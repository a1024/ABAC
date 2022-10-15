//aac.c - Arithmetic Coder implementation, CPU-side
//Copyright (C) 2022  Ayman Wagih Mohsen
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include		"awm_ac.h"
#include		<stdio.h>
#include		<stdlib.h>
#include		<string.h>
#define _USE_MATH_DEFINES
#include		<math.h>
#include		<immintrin.h>
//#include		<tmmintrin.h>
#ifdef __GNUC__
#include		<x86intrin.h>
#else
#include		<intrin.h>
#endif
#include		<stdarg.h>
static const char file[]=__FILE__;


//	#define		ENABLE_GRAYCODE
//	#define		ENABLE_SYMBOL_ESTIMATION
	#define		ENABLE_CONFIDENCE//FIXME commenting this breaks abac0a

//	#define		PRINT_SSE2


#define			AC_MEASURE_PREDICTION
#define			LOG_WINDOW_SIZE		16	//[2, 16]	do not change
#define			LOG_CONFBOOST		14
#define			ABAC2_CONF_MSB_RELATION
#define			WINDOW_SIZE			(1<<LOG_WINDOW_SIZE)
#define			PROB_MASK			(WINDOW_SIZE-1)
#define			PROB_MAX			(WINDOW_SIZE-2)
#define			PROB_HALF			(1<<(LOG_WINDOW_SIZE-1))
#define			PROB_INIT			(PROB_HALF-1)
#define			BOOST_POWER			4.
#define			MIN_CONF			0.55
#define			FAIL(MSG, ...)		return LOG_ERROR(MSG, ##__VA_ARGS__), 0
#define			PROB_CLAMP(PROB)	PROB=PROB<1?1:(PROB_MAX<PROB?PROB_MAX:PROB)


typedef unsigned long long u64;
const int		magic_ac04='A'|'C'<<8|'0'<<16|'4'<<24,//ABAC
				magic_ac05='A'|'C'<<8|'0'<<16|'5'<<24,//ANS
				magic_ac06='A'|'C'<<8|'0'<<16|'6'<<24,//64bit ANS
				magic_ac07='A'|'C'<<8|'0'<<16|'7'<<24,//SSE2 ANS
				magic_ac08='A'|'C'<<8|'0'<<16|'8'<<24,//AVX2 ANS

				magic_ac09='A'|'C'<<8|'0'<<16|'9'<<24,//OpenCL ABAC
				magic_an09='A'|'N'<<8|'0'<<16|'9'<<24,//OpenCL ANS

				magic_ac0A='A'|'C'<<8|'0'<<16|'A'<<24,//8-bit SSE2 ABAC
				magic_ac0T='A'|'C'<<8|'0'<<16|'T'<<24;//testbench


#ifdef PRINT_SSE2
static void		print_reg(__m128i const *r, int n, const char *format, ...)
{
	va_list args;

	if(n==16)
	{
		for(int k=0;k<8;++k)
			printf("%04X ", r->m128i_i16[7-k]);
	}
	else if(n==32)
	{
		for(int k2=1;k2>=0;--k2)
			for(int k=3;k>=0;--k)
				printf("%08X ", r[k2].m128i_i32[k]);
	}
	if(format)
	{
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
	}
	printf("\n");
}
static void		print_prob(int prob, int w, const char *format, ...)
{
	va_list args;

	for(int k=0;k<w;++k)
		printf("%c", '1'-(k*0x7FFF<prob*w));
	if(format)
	{
		printf(" ");
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
	}
	printf("\n");
}
#else
#define			print_reg(...)
#define			print_prob(...)
#endif


//aligned malloc	//https://stackoverflow.com/questions/38088732/explanation-to-aligned-malloc-implementation
void* aligned_malloc(size_t size, size_t align)
{
    void *p0, **pA;
	size_t offset;

    offset=align-1+sizeof(void*);
	p0=(void*)malloc(size+offset);
    if(!p0)
       return 0;
    pA=(void**)(((size_t)p0+offset)&~(align-1));
    pA[-1]=p0;
    return pA;
}
void aligned_free(void *p)
{
	if(p)
	    free(((void**)p)[-1]);
}


//JPEG

//https://en.wikipedia.org/wiki/YCbCr#JPEG_conversion
void			YCbCr_fwd(const int *image, float *Y, float *Cb, float *Cr, int count)
{
	int k, R, G, B;
	const int *v;

	for(k=0;k<count;++k)
	{
		v=image+k;
		R=*v&0xFF;
		G=*v>>8&0xFF;
		B=*v>>16&0xFF;
		Y[k]=0.299f*R+0.587f*G+0.114f*B;
		Cb[k]=128-0.168736f*R-0.331264f*G+0.5f*B;
		Cr[k]=128+0.5f*R+0.418688f*G+0.081312f*B;
	}
}
void			YCbCr_inv(int *image, const float *Y, const float *Cb, const float *Cr, int count)
{
	int k, R, G, B;
	float Cb2, Cr2;

	for(k=0;k<count;++k)
	{
		Cb2=Cb[k]-128;
		Cr2=Cr[k]-128;
		R=(int)(Y[k]+1.402f*Cr2);
		G=(int)(Y[k]-0.344136f*Cb2-0.714136f*Cr2);
		B=(int)(Y[k]+1.772f*Cb2);
		image[k]=B<<16|G<<8|R;
	}
}

//DCT 8x8 float
void			DCT8_ref(double *data, int inv)
{
#define C0_5	0.5*0.9807852804032304491261822361342
#define C1_0	0.5*0.9238795325112867561281831893968
#define C1_5	0.5*0.8314696123025452370787883776179
#define C2_0	0.5*0.7071067811865475244008443621048
#define C2_5	0.5*0.5555702330196022247428308139485
#define C3_0	0.5*0.3826834323650897717284599840304
#define C3_5	0.5*0.195090322016128267848284868477
	//https://en.wikipedia.org/wiki/Discrete_cosine_transform#DCT-II
	static const double f_coeff[]=//X[k] = sqrt(2/N) * sum x[n] * cos(pi/N*(n+1/2)*k)
	{//k \ n
		0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
		C0_5,  C1_5,  C2_5,  C3_5, -C3_5, -C2_5, -C1_5, -C0_5,
		C1_0,  C3_0, -C3_0, -C1_0, -C1_0, -C3_0,  C3_0,  C1_0,
		C1_5, -C3_5, -C0_5, -C2_5,  C2_5,  C0_5,  C3_5, -C1_5,
		C2_0, -C2_0, -C2_0,  C2_0,  C2_0, -C2_0, -C2_0,  C2_0,
		C2_5, -C0_5,  C3_5,  C1_5, -C1_5, -C3_5,  C0_5, -C2_5,
		C3_0, -C1_0,  C1_0, -C3_0, -C3_0,  C1_0, -C1_0,  C3_0,
		C3_5, -C2_5,  C1_5, -C0_5,  C0_5, -C1_5,  C2_5, -C3_5,
	};
	//https://en.wikipedia.org/wiki/Discrete_cosine_transform#DCT-III
	static const double i_coeff[]=//X[k] = sqrt(2/N) * (x0/2 + sum x[n] * cos(pi/N*(k+1/2)*n))
	{//k \ n
		0.25,  C0_5,  C1_0,  C1_5,  C2_0,  C2_5,  C3_0,  C3_5,
		0.25,  C1_5,  C3_0, -C3_5, -C2_0, -C0_5, -C1_0, -C2_5,
		0.25,  C2_5, -C3_0, -C0_5, -C2_0,  C3_5,  C1_0,  C1_5,
		0.25,  C3_5, -C1_0, -C2_5,  C2_0,  C1_5, -C3_0, -C0_5,
		0.25, -C3_5, -C1_0,  C2_5,  C2_0, -C1_5, -C3_0,  C0_5,
		0.25, -C2_5, -C3_0,  C0_5, -C2_0, -C3_5,  C1_0, -C1_5,
		0.25, -C1_5,  C3_0,  C3_5, -C2_0,  C0_5, -C1_0,  C2_5,
		0.25, -C0_5,  C1_0, -C1_5,  C2_0, -C2_5,  C3_0, -C3_5,
	};
#undef	C0_5
#undef	C1_0
#undef	C1_5
#undef	C2_0
#undef	C2_5
#undef	C3_0
#undef	C3_5
	const double *coeff=inv?i_coeff:f_coeff;
	double d2[8];
	int k, n;

	for(k=0;k<8;++k)
	{
		d2[k]=0;
		for(n=0;n<8;++n)
			d2[k]+=coeff[k<<3|n]*data[n];
	}
	memcpy(data, d2, 8*sizeof(double));
}
#if 0
static const double S[]=//scaling coefficients
{
	0.353553390593273762200422/2,
	0.254897789552079584470970,
	0.270598050073098492199862,
	0.300672443467522640271861,
	0.353553390593273762200422,
	0.449988111568207852319255,
	0.653281482438188263928322,
	1.281457723870753089398043,
};

static const double A[]=
{
	NAN,//0
	0.707106781186547524400844,//1
	0.541196100146196984399723,//2
	0.707106781186547524400844,//3
	1.306562964876376527856643,//4
	0.382683432365089771728460,//5
};
#endif
void			DCT8f_fwd(double *data)
{
#if 1
	static const double coeff[]=
	{
		0.382683432365089771728460,//1
		0.541196100146196984399723,//2
		0.707106781186547524400844,//0
		1.306562964876376527856643,
	};

	double
		a0=data[0]+data[7],
		a1=data[1]+data[6],
		a2=data[2]+data[5],
		a3=data[3]+data[4],
		a4=data[3]-data[4],
		a5=data[2]-data[5],
		a6=data[1]-data[6],
		a7=data[0]-data[7];
	
	double
		b0=a0+a3,
		b1=a1+a2,
		b2=a1-a2,
		b3=a0-a3,
		b4=-a4-a5,
		b5=(a5+a6)*coeff[2],
		b6=a6+a7;
	
	double c0=(b4+b6)*coeff[0];

	double
		d0=(b2+b3)*coeff[2],
		d1=-b4*coeff[1]-c0,
		d2=b6*coeff[3]-c0,
		d3=b5+a7,
		d4=a7-b5;
	
	data[0]=(b0+b1)*0.125;
	data[1]=(d3+d2)*0.125;
	data[2]=(d0+b3)*0.125;
	data[3]=(d4-d1)*0.125;
	data[4]=(b0-b1)*0.125;
	data[5]=(d1+d4)*0.125;
	data[6]=(b3-d0)*0.125;
	data[7]=(d3-d2)*0.125;

	//data[0]=(b0+b1)*S[0];
	//data[1]=(d3+d2)*S[1];
	//data[2]=(d0+b3)*S[2];
	//data[3]=(d4-d1)*S[3];
	//data[4]=(b0-b1)*S[4];
	//data[5]=(d1+d4)*S[5];
	//data[6]=(b3-d0)*S[6];
	//data[7]=(d3-d2)*S[7];
#endif

	//https://www.nayuki.io/page/fast-discrete-cosine-transform-algorithms
#if 0
	const double v0=data[0]+data[7];
	const double v1=data[1]+data[6];
	const double v2=data[2]+data[5];
	const double v3=data[3]+data[4];
	const double v4=data[3]-data[4];
	const double v5=data[2]-data[5];
	const double v6=data[1]-data[6];
	const double v7=data[0]-data[7];
	
	const double v8=v0+v3;
	const double v9=v1+v2;
	const double v10=v1-v2;
	const double v11=v0-v3;
	const double v12=-v4-v5;
	const double v13=(v5+v6)*A[3];
	const double v14=v6+v7;
	
	const double v15=v8+v9;
	const double v16=v8-v9;
	const double v17=(v10+v11)*A[1];
	const double v18=(v12+v14)*A[5];
	
	const double v19=-v12*A[2]-v18;
	const double v20=v14*A[4]-v18;
	
	const double v21=v17+v11;
	const double v22=v11-v17;
	const double v23=v13+v7;
	const double v24=v7-v13;
	
	const double v25=v19+v24;
	const double v26=v23+v20;
	const double v27=v23-v20;
	const double v28=v24-v19;
	
	data[0]=S[0]*v15;
	data[1]=S[1]*v26;
	data[2]=S[2]*v21;
	data[3]=S[3]*v28;
	data[4]=S[4]*v16;
	data[5]=S[5]*v25;
	data[6]=S[6]*v22;
	data[7]=S[7]*v27;
#endif

	//Plonka 2002 DCT-II	cosk = cos(k*pi/2N), here N=8		//https://fgiesen.wordpress.com/2013/11/04/bink-2-2-integer-dct-design-part-1/
#if 0
	const float
		cos7=0.195090322016128267848284868477f,
		sin7=0.9807852804032304491261822361342f,
		cos5=0.5555702330196022247428308139485f,
		sin5=0.8314696123025452370787883776179f,
		cos6=0.3826834323650897717284599840304f,
		sin6=0.9238795325112867561281831893968f;

	float
		a0=data[0]+data[7],
		a1=data[1]+data[6],
		a2=data[2]+data[5],
		a3=data[3]+data[4],
		a4=data[3]-data[4],
		a5=data[2]-data[5],
		a6=data[1]-data[6],
		a7=data[0]-data[7];

	float
		b0=a0+a3,
		b1=a1+a2,
		b2=a1-a2,
		b3=a0-a3,
		b4=a4*cos7+a7*sin7,
		b5=a5*cos5+a6*sin5,
		b6=a6*cos5-a5*sin5,
		b7=a7*cos7-a4*sin7;

	a0=b0+b1;
	a1=b0-b1;
	a2=b2*cos6+b3*sin6;
	a3=b3*cos6-b2*sin6;
	a4=b4+b5;
	a5=b4-b5;
	a6=b6+b7;
	a7=b6-b7;

	data[0]=a0*0.125;
	data[1]=a4*(0.125*M_SQRT2);
	data[2]=a2*(0.125*M_SQRT2);
	data[3]=a5*0.125;
	data[4]=a1*0.125;
	data[5]=a6*0.125;
	data[6]=a3*(0.125*M_SQRT2);
	data[7]=a7*(0.125*-M_SQRT2);
#endif
}
void			DCT8f_inv(double *data)
{
#if 1
	static const double coeff[]=
	{
		 0.76536686473017954345692,
		-2.613125929752753055713286,
		 1.082392200292393968799446,
		 1.41421356237309504880168872421,
	};
	//data[0]/=S[0];//0 BRP is its own inverse
	//data[1]/=S[1];//4
	//data[2]/=S[2];//2
	//data[3]/=S[3];//6
	//data[4]/=S[4];//1
	//data[5]/=S[5];//5
	//data[6]/=S[6];//3
	//data[7]/=S[7];//7
	
	double
		a0=data[0]+data[4],
		a1=data[2]+data[6],
		a2=data[2]-data[6],
		a3=data[0]-data[4],
		a4=data[1]+data[7],
		a5=data[5]+data[3],
		a7=data[5]-data[3],
		a6=data[1]-data[7];

	const double v0=a0+a1;
	const double v3=a0-a1;
	const double b0=a4+a5;
	const double b1=a4-a5;
	const double v18=(a7-a6)*coeff[0];
	const double v10=a2*coeff[3]-a1;

	const double v1=a3+v10;
	const double v2=a3-v10;
	const double v12=a7*coeff[1]+v18;
	const double v14=a6*coeff[2]-v18;
	
	const double v6=v14-b0;

	const double v5=b1*coeff[3]-v6;
	
	const double v4=-v5-v12;
	
	data[0]=v0+b0;
	data[1]=v1+v6;
	data[2]=v2+v5;
	data[3]=v3+v4;
	data[4]=v3-v4;
	data[5]=v2-v5;
	data[6]=v1-v6;
	data[7]=v0-b0;
#endif

	//https://www.nayuki.io/page/fast-discrete-cosine-transform-algorithms
#if 0
	const double v15=data[0]/S[0];
	const double v26=data[1]/S[1];
	const double v21=data[2]/S[2];
	const double v28=data[3]/S[3];
	const double v16=data[4]/S[4];
	const double v25=data[5]/S[5];
	const double v22=data[6]/S[6];
	const double v27=data[7]/S[7];
	
	const double v19=(v25-v28)/2;
	const double v20=(v26-v27)/2;
	const double v23=(v26+v27)/2;
	const double v24=(v25+v28)/2;
	
	const double v7 =(v23+v24)/2;
	const double v11=(v21+v22)/2;
	const double v13=(v23-v24)/2;
	const double v17=(v21-v22)/2;
	
	const double v8=(v15+v16)/2;
	const double v9=(v15-v16)/2;
	
	const double v18=(v19-v20)*A[5];//Different from original
	const double v12=(v19*A[4]-v18)/(A[2]*A[5]-A[2]*A[4]-A[4]*A[5]);//den=-1
	const double v14=(v18-v20*A[2])/(A[2]*A[5]-A[2]*A[4]-A[4]*A[5]);
	
	const double v6=v14-v7;
	const double v5=v13/A[3]-v6;
	const double v4=-v5-v12;
	const double v10=v17/A[1]-v11;
	
	const double v0=(v8+v11)/2;
	const double v1=(v9+v10)/2;
	const double v2=(v9-v10)/2;
	const double v3=(v8-v11)/2;
	
	data[0]=(v0+v7)*0.5;
	data[1]=(v1+v6)*0.5;
	data[2]=(v2+v5)*0.5;
	data[3]=(v3+v4)*0.5;
	data[4]=(v3-v4)*0.5;
	data[5]=(v2-v5)*0.5;
	data[6]=(v1-v6)*0.5;
	data[7]=(v0-v7)*0.5;
#endif
}
void			DCT8x8f_fwd(float *data, int bw, int bh)
{
}


//https://en.wikipedia.org/wiki/YCoCg#The_lifting-based_YCoCg-R_variation
void			YCoCg_fwd(int *image, int count)
{
	int *v, R, G, B, Y, Co, Cg;
	for(int k=0;k<count;++k)
	{
		v=image+k;
		R=*v&0xFF;
		G=*v>>8&0xFF;
		B=*v>>16&0xFF;

		Co=R-B;
		B+=Co>>1;
		Cg=G-B;
		Y=B+(Cg>>1);

		*v=Cg<<17|Co<<8|Y;
	}
}
void			YCoCg_inv(int *image, int count)
{
	int *v, R, G, B, Y, Co, Cg;
	for(int k=0;k<count;++k)
	{
		v=image+k;
		Y=*v&0xFF;
		Co=*v>>8&0x1FF;
		Cg=*v>>17&0x1FF;

		R=Y-(Cg>>1);
		G=Cg+R;
		B=R-(Co>>1);
		R=B+Co;

		*v=B<<16|G<<8|R;
	}
}


//test: encodes YCoCg26 pixels stored in 32bit integers, no alpha.
//src pixel: MSB {6 zeros}{9 Cg}{9 Co}{8 Y} LSB
int				test_encode(const int *src, int iw, int ih, int bpp, Predictor predictor, ArrayHandle *output, int loud)
{
	DList list[32];
	unsigned start[32]={0}, range[32];
	unsigned short prob[32]={0}, hit_hist[32];
	u64 predicted=0;
	for(int kb=0;kb<bpp;++kb)
	{
		range[kb]=0xFFFFFFFF;
		prob[kb]=0x8000;
		hit_hist[kb]=0x8000;
		dlist_init(list+kb, 1, 128, 0);
	}
	for(int ky=0;ky<ih;++ky)
	{
		printf("\r%lf%%\t\t", 100.*ky/(ih-1));
		for(int kx=0;kx<iw;++kx)
		{
			predictor(src, iw, ih, bpp, kx, ky, hit_hist, prob);
			for(int kb=0;kb<bpp;++kb)
			{
				unsigned r2=(unsigned)((long long)range[kb]*prob[kb]>>16);
				r2+=(r2==0)-(r2==range[kb]);

				int bit=src[iw*ky+kx]>>kb&1;

				if(bit)
				{
					++r2;
					start[kb]+=r2, range[kb]-=r2;
				}
				else
					range[kb]=r2-1;
				
				int correct=bit^(prob[kb]>=PROB_HALF);
				predicted+=correct;
				hit_hist[kb]=correct<<15|hit_hist[kb]>>1;

				//renormalize
				while((start[kb]^(start[kb]+(unsigned)range[kb]))<0x1000000)//most significant byte has stabilized			zpaq 1.10
				{
					dlist_push_back1(list+kb, (char*)(start+kb)+3);
					start[kb]<<=8;
					range[kb]=range[kb]<<8|0xFF;
				}
				if(range[kb]<3)
				{
					dlist_push_back1(list+kb, (char*)(start+kb)+3);//big endian
					dlist_push_back1(list+kb, (char*)(start+kb)+2);
					dlist_push_back1(list+kb, (char*)(start+kb)+1);
					dlist_push_back1(list+kb, start+kb);

					start[kb]=0, range[kb]=0xFFFFFFFF;//because 1=0.9999...
				}
			}
		}
	}
	printf("\n");
	size_t offset0;
	if(*output)
	{
		if(output[0]->esize!=1)
			FAIL("Output array element size must be 1 byte");
		offset0=output[0]->count;
		ARRAY_APPEND(*output, 0, (1+8)*sizeof(int), 1, 0);
	}
	else
	{
		offset0=0;
		ARRAY_ALLOC(char, *output, 0, (1+32)*sizeof(int), 0, 0);//magic + sizes[bitdepth]
	}
	memcpy(output[0]->data+offset0, &magic_ac0T, sizeof(int));
	for(int kb=0;kb<bpp;++kb)
		memcpy(output[0]->data+offset0+(1+kb)*sizeof(int), &list[kb].nobj, sizeof(int));
	if(loud)
	{
		int symcount=iw*ih;
		u64 original_bitsize=(u64)symcount*bpp, compressed_bitsize=0;
		for(int kb=0;kb<bpp;++kb)
			compressed_bitsize+=list[kb].nobj;
		compressed_bitsize<<=3;
		printf("TEST ENC: (%dx%dx%d/8) %lld -> %lld,  ratio: %lf,  %lf bpp\n", iw, ih, bpp, original_bitsize>>3, compressed_bitsize>>3, (double)original_bitsize/compressed_bitsize, (double)compressed_bitsize/symcount);
		printf("Predicted: %6lld / %6lld = %lf%%\n", predicted, original_bitsize, 100.*predicted/original_bitsize);
		printf("Bit\tbytes\tratio,\tbytes/bitplane = %d\n", (symcount+7)>>3);
		for(int kb=0;kb<bpp;++kb)
			printf("%2d\t%5d\t%lf\n", kb, list[kb].nobj, (double)symcount/(list[kb].nobj<<3));
		
		printf("Preview:\n");
		int kprint=(int)(output[0]->count-offset0);
		if(kprint>200)
			kprint=200;
		for(int k=0;k<kprint;++k)
			printf("%02X-", output[0]->data[offset0+k]&0xFF);
		printf("\n");
	}
	for(int kb=0;kb<bpp;++kb)
		dlist_clear(list+kb);
#if 0
	//3 2
	//1[0]	causal neighbors
	int *ptr, wnd[4], val, Y[4], Co[4], Cg[4];
	ptr=src;
	for(int ky=0;ky<ih;++ky)
	{
		wnd[3]=0;
		wnd[1]=0;
		for(int kx=0;kx<iw;++kx, ++ptr)
		{
			wnd[2]=ky?ptr[-iw]:0;
			wnd[0]=*ptr;

			Y[0]=wnd[0]&0xFF, Co[0]=wnd[0]>>8&0x1FF, Cg[0]=wnd[0]>>17&0x1FF;
			Y[1]=wnd[1]&0xFF, Co[1]=wnd[1]>>8&0x1FF, Cg[1]=wnd[1]>>17&0x1FF;
			Y[2]=wnd[2]&0xFF, Co[2]=wnd[2]>>8&0x1FF, Cg[2]=wnd[2]>>17&0x1FF;
			Y[3]=wnd[3]&0xFF, Co[3]=wnd[3]>>8&0x1FF, Cg[3]=wnd[3]>>17&0x1FF;

			//prediction
			Y[1]+=Y[3]-Y[2];
			Co[1]+=Co[3]-Co[2];
			Cg[1]+=Cg[3]-Cg[2];

			//clamp
			Y[1]=Y[1]>0x1FF?0x1FF:(Y[1]<0?0:Y[1]);
			Co[1]=Co[1]>0x1FF?0x1FF:(Co[1]<0?0:Co[1]);
			Cg[1]=Cg[1]>0x1FF?0x1FF:(Cg[1]<0?0:Cg[1]);

			//Y[0]-=Y[1]+Y[3]-Y[2];//this can be done beforehand with non-causal neighbors
			//Co[0]-=Co[1]+Co[3]-Co[2];
			//Cg[0]-=Cg[1]+Cg[3]-Cg[2];
			//
			//Y[0]=abs(Y[0])<<1|(Y[0]<0);//how many bits is this now?
			//Co[0]=abs(Co[0])<<1|(Co[0]<0);
			//Cg[0]=abs(Cg[0])<<1|(Cg[0]<0);

			wnd[3]=wnd[2];
			wnd[1]=wnd[0];
		}
	}
#endif
	return 1;
}
//const void*		test_encode(const void *src_start, const void *src_end, int *dst, int iw, int ih, int loud)
//{
//}

//abac0a_encode(): Encodes 8-bit symbols. Uses SSSE3 with 15-bit probability.
//If output was initialized, output->esize must be 1.
int				abac0a_encode(const unsigned char *src, int count, int bytestride, ArrayHandle *output, int loud)
{
	DList list[8];
	const unsigned char *buffer;
	unsigned char *header;
	size_t offset0;
	u64 t1, t2;
	__m128i
		ones=_mm_set1_epi32(-1), ones32=_mm_set1_epi32(1), three=_mm_set1_epi32(3),
		prob_min=_mm_set1_epi16(1), prob_half=_mm_set1_epi16(0x4000), prob_max=_mm_set1_epi16(0x7FFE),
		prob=prob_half,
#ifdef ENABLE_CONFIDENCE
		prob_correct=prob_half,
#endif
		limit=_mm_set1_epi32(0x1000000);
	__m128i start[2], range[2], r2[2];
	__m128i shuf_pack=_mm_set_epi8(15, 14, 11, 10, 7, 6, 3, 2, 13, 12, 9, 8, 5, 4, 1, 0);

	if(!src||!count||!bytestride)
		FAIL("abac4_encode(src=%p, count=%d, stride=%d)", src, count, bytestride);

	t1=__rdtsc();
	if(*output)
	{
		if(output[0]->esize!=1)
			FAIL("Output array element size must be 1 byte");
		offset0=output[0]->count;
		ARRAY_APPEND(*output, 0, (1+8)*sizeof(int), 1, 0);
	}
	else
	{
		offset0=0;
		ARRAY_ALLOC(char, *output, 0, (1+8)*sizeof(int), 0, 0);//magic + sizes[bitdepth]
	}

	header=(unsigned char*)array_back(output)+1-(1+8)*sizeof(int);
	memcpy(header, &magic_ac0A, sizeof(int));
	header+=sizeof(int);

	for(int k=0;k<8;++k)
		dlist_init(list+k, 1, 1024, 0);

	buffer=(const unsigned char*)src;
	start[0]=_mm_setzero_si128();
	start[1]=_mm_setzero_si128();
	range[0]=_mm_set1_epi32(0x7FFFFFFF);//the MSB is the 2's comp sign, all SSE2 cmp instructions are signed
	range[1]=_mm_set1_epi32(0x7FFFFFFF);
	//unsigned start[8]={0}, r2[8], range[8];
	//for(int k=0;k<8;++k)
	//	range[k]=0xFFFFFFFF;
	for(int ks=0, ks2=0;ks<count;++ks, ks2+=bytestride)
	{
		print_reg(start, 32, "%d start", ks);
		print_reg(range, 32, "%d range", ks);

		//calculate probability of zero bit
#ifdef ENABLE_CONFIDENCE
		__m128i p0=_mm_sub_epi16(prob, prob_half);

		p0=_mm_slli_epi16(p0, 1);
		prob_correct=_mm_slli_epi16(prob_correct, 1);
		p0=_mm_mulhi_epu16(p0, prob_correct);
		p0=_mm_mulhi_epu16(p0, prob_correct);
		p0=_mm_srli_epi16(p0, 1);
		prob_correct=_mm_srli_epi16(prob_correct, 1);

		p0=_mm_add_epi16(p0, prob_half);
#else
		__m128i p0=prob;
#endif

		//clamp probability			eg: {0, 1, half, FFFF}
		__m128i cmp=_mm_cmplt_epi16(prob_min, p0);//{0, FFFF, FFFF, FFFF}
		p0=_mm_and_si128(p0, cmp);
		cmp=_mm_xor_si128(cmp, ones);
		cmp=_mm_and_si128(cmp, prob_min);
		p0=_mm_or_si128(p0, cmp);

		cmp=_mm_cmplt_epi16(p0, prob_max);//{FFFF, FFFF, FFFF, 0}
		p0=_mm_and_si128(p0, cmp);
		cmp=_mm_xor_si128(cmp, ones);
		cmp=_mm_and_si128(cmp, prob_max);
		p0=_mm_or_si128(p0, cmp);
		
		//calculate middle (r2)
		__m128i t0=_mm_shuffle_epi8(range[0], shuf_pack);
		__m128i t1=_mm_shuffle_epi8(range[1], shuf_pack);
		__m128i r_lo=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(t0), _mm_castsi128_ps(t1), _MM_SHUFFLE(1, 0, 1, 0)));
		__m128i r_hi=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(t0), _mm_castsi128_ps(t1), _MM_SHUFFLE(3, 2, 3, 2)));
		p0=_mm_slli_epi16(p0, 1);
		__m128i a=_mm_mulhi_epu16(r_lo, p0);
		__m128i b0=_mm_mullo_epi16(r_hi, p0);
		__m128i b1=_mm_mulhi_epu16(r_hi, p0);
		p0=_mm_srli_epi16(p0, 1);

		r2[0]=_mm_unpacklo_epi16(b0, b1);
		r2[1]=_mm_unpackhi_epi16(b0, b1);

		__m128i d0=_mm_unpacklo_epi16(a, _mm_setzero_si128());
		__m128i d1=_mm_unpackhi_epi16(a, _mm_setzero_si128());

		r2[0]=_mm_add_epi32(r2[0], d0);
		r2[1]=_mm_add_epi32(r2[1], d1);

		//avoid degenerate range		r2[k]+=(r2[k]==0)-(r2[k]==range[k]);
		cmp=_mm_cmpeq_epi32(r2[0], _mm_setzero_si128());
		cmp=_mm_and_si128(cmp, ones32);
		__m128i cmp2=_mm_cmpeq_epi32(r2[0], range[0]);
		cmp2=_mm_and_si128(cmp2, ones32);
		cmp=_mm_sub_epi32(cmp, cmp2);
		r2[0]=_mm_add_epi32(r2[0], cmp);
		
		cmp=_mm_cmpeq_epi32(r2[0], _mm_setzero_si128());
		cmp=_mm_and_si128(cmp, ones32);
		cmp2=_mm_cmpeq_epi32(r2[0], range[1]);
		cmp2=_mm_and_si128(cmp2, ones32);
		cmp=_mm_sub_epi32(cmp, cmp2);
		r2[1]=_mm_add_epi32(r2[1], cmp);

		print_reg(r2, 32, "%d r2", ks);
		
		print_reg(&prob, 16, "%d prob ", ks);
		print_prob(prob.m128i_u16[7], 64, "bit 7");
		print_reg(&p0, 16, "%d p0 ", ks);
		print_prob(p0.m128i_u16[7], 64, "bit 7");


		unsigned char sym=src[ks2];
#ifdef ENABLE_GRAYCODE
		sym^=sym>>1;
#endif
		__m128i bit=_mm_set_epi16(sym>>7&1, sym>>6&1, sym>>5&1, sym>>4&1, sym>>3&1, sym>>2&1, sym>>1&1, sym>>0&1);
		
#ifdef ENABLE_CONFIDENCE
		//update prob_correct			correct=bit^(p0>=PROB_HALF), prob_correct=correct<<15|prob_correct>>1;
		cmp=_mm_cmplt_epi16(p0, prob_half);
		cmp=_mm_xor_si128(cmp, ones);
		cmp=_mm_and_si128(cmp, prob_min);
		cmp=_mm_xor_si128(cmp, bit);

		cmp=_mm_slli_epi16(cmp, 14);//set pre-MSB, because setting MSB would negate the 2's comp value
		prob_correct=_mm_srli_epi16(prob_correct, 1);
		prob_correct=_mm_or_si128(prob_correct, cmp);
#endif

		//update prob zero				prob=!bit<<15|prob>>1;
		__m128i bit2=_mm_xor_si128(bit, prob_min);
		prob=_mm_srli_epi16(prob, 1);
		prob=_mm_or_si128(prob, _mm_slli_epi16(bit2, 14));
		
		print_reg(&prob, 16, "%d prob2 ", ks);
		print_reg(&bit, 16, "%d bit ", ks);


		//update range
		__m128i bitlo=_mm_unpacklo_epi16(bit, _mm_setzero_si128());
		__m128i bithi=_mm_unpackhi_epi16(bit, _mm_setzero_si128());

		//r2[k]+=bit.m128i_i16[k]-!bit.m128i_i16[k];
		r2[0]=_mm_add_epi32(r2[0], bitlo);
		r2[1]=_mm_add_epi32(r2[1], bithi);
		r2[0]=_mm_sub_epi32(r2[0], _mm_sub_epi32(ones32, bitlo));
		r2[1]=_mm_sub_epi32(r2[1], _mm_sub_epi32(ones32, bithi));

		//start[k]+=r2[k]&-bit.m128i_i16[k];
		bitlo=_mm_cmpeq_epi32(bitlo, ones32);
		bithi=_mm_cmpeq_epi32(bithi, ones32);
		t0=_mm_and_si128(bitlo, r2[0]);
		t1=_mm_and_si128(bithi, r2[1]);
		start[0]=_mm_add_epi32(start[0], t0);
		start[1]=_mm_add_epi32(start[1], t1);
		//__m128i start_lo=_mm_loadu_si128((const __m128i*)start);
		//__m128i start_hi=_mm_loadu_si128((const __m128i*)(start+4));

		//range[k]=((range[k]-r2[k])&-bit.m128i_i16[k])|(r2[k]&-!bit.m128i_i16[k]);
		range[0]=_mm_sub_epi32(range[0], r2[0]);
		range[1]=_mm_sub_epi32(range[1], r2[1]);
		range[0]=_mm_and_si128(range[0], bitlo);
		range[1]=_mm_and_si128(range[1], bithi);
		bitlo=_mm_xor_si128(bitlo, ones);
		bithi=_mm_xor_si128(bithi, ones);
		r2[0]=_mm_and_si128(r2[0], bitlo);
		r2[1]=_mm_and_si128(r2[1], bithi);
		range[0]=_mm_or_si128(range[0], r2[0]);
		range[1]=_mm_or_si128(range[1], r2[1]);

		//renormalize
		unsigned condition;
		for(;;)
		{
			t0=_mm_add_epi32(start[0], range[0]);
			t1=_mm_add_epi32(start[1], range[1]);
			t0=_mm_xor_si128(t0, start[0]);
			t1=_mm_xor_si128(t1, start[1]);
			t0=_mm_cmpeq_epi32(t0, limit);
			t1=_mm_cmpeq_epi32(t1, limit);
			condition=_mm_movemask_epi8(t1)<<16|_mm_movemask_epi8(t0);
			if(!condition)
				break;
			for(int k=0;k<8;++k)
			{
				if(condition>>(k<<2)&15)
				{
					int reg_idx=k>>2, idx=k&3;
					dlist_push_back1(list+k, (char*)&start[reg_idx].m128i_i32[idx]+3);
#ifdef PRINT_SSE2
					printf("RENORM bit %d <-[%02X]%02X%04X\n", k, ((char*)&start[reg_idx].m128i_i32[idx])[3], ((char*)&start[reg_idx].m128i_i32[idx])[2], ((short*)&start[reg_idx].m128i_i32[idx])[0]);
#endif
					start[reg_idx].m128i_i32[idx]<<=8;
					range[reg_idx].m128i_i32[idx]=range[reg_idx].m128i_i32[idx]<<8|0xFF;
				}
			}
		}
		t0=_mm_cmplt_epi32(range[0], three);
		t1=_mm_cmplt_epi32(range[1], three);
		condition=_mm_movemask_epi8(t1)<<16|_mm_movemask_epi8(t0);
		if(condition)
		{
			for(int k=0;k<8;++k)
			{
				if(condition>>(k<<2)&15)
				{
					int reg_idx=k>>2, idx=k&3;
					dlist_push_back1(list+k, (char*)&start[reg_idx].m128i_i32[idx]+3);//big endian
					dlist_push_back1(list+k, (char*)&start[reg_idx].m128i_i32[idx]+2);
					dlist_push_back1(list+k, (char*)&start[reg_idx].m128i_i32[idx]+1);
					dlist_push_back1(list+k, &start[reg_idx].m128i_i32[idx]);
#ifdef PRINT_SSE2
					printf("RENORM bit %d <-[%08X]\n", k, start[reg_idx].m128i_i32[idx]);
#endif

					start[reg_idx].m128i_i32[idx]=0;
					range[reg_idx].m128i_i32[idx]=0x7FFFFFFF;//because 1=0.9999...
				}
			}
		}
#ifdef PRINT_SSE2
		printf("\n");
#endif
	}
	for(int k=0;k<8;++k)
	{
		int reg_idx=k>>2, idx=k&3;
		dlist_push_back1(list+k, (char*)&start[reg_idx].m128i_i32[idx]+3);//big endian
		dlist_push_back1(list+k, (char*)&start[reg_idx].m128i_i32[idx]+2);
		dlist_push_back1(list+k, (char*)&start[reg_idx].m128i_i32[idx]+1);
		dlist_push_back1(list+k, &start[reg_idx].m128i_i32[idx]);
	}
	for(int k=0;k<8;++k)
		memcpy(header+k*sizeof(int), &list[k].nobj, sizeof(int));
	for(int k=0;k<8;++k)//header pointer is now invalid
		dlist_appendtoarray(list+k, output);
	t2=__rdtsc();

	if(loud)
	{
		size_t total_csize=output[0]->count-offset0;
		printf("AC encode:  %lld cycles, %lf c/byte\n", t2-t1, (double)(t2-t1)/count);
		printf("Size: %d -> %d, ratio: %lf, %lf bpp\n", count, total_csize, (double)count/total_csize, (double)(total_csize<<3)/count);
#ifdef AC_MEASURE_PREDICTION
		//printf("Predicted: %6lld / %6lld = %lf%%\n", hitnum, hitden, 100.*hitnum/hitden);
#endif
		printf("Bit\tbytes\tratio,\tbytes/bitplane = %d\n", (count+7)>>3);
		for(int k=0;k<8;++k)
			printf("%2d\t%5d\t%lf\n", k, list[k].nobj, (double)count/(list[k].nobj<<3));
		
		printf("Preview:\n");
		int kprint=total_csize<200?(int)total_csize:200;
		for(int k=0;k<kprint;++k)
			printf("%02X-", output[0]->data[offset0+k]&0xFF);
		printf("\n");
	}

	for(int k=0;k<8;++k)
		dlist_clear(list+k);

	return 1;
}
const void*		abac0a_decode(const void *src_start, const void *src_end, unsigned char *dst, int count, int bytestride, int loud)
{
	const unsigned char *sizes, *data;
	unsigned char *buffer;
	const unsigned char *ptr[8], *end[8];
	u64 t1, t2;
	__m128i
		ones=_mm_set1_epi32(-1), ones32=_mm_set1_epi32(1), three=_mm_set1_epi32(3),
		prob_min=_mm_set1_epi16(1), prob_half=_mm_set1_epi16(0x4000), prob_max=_mm_set1_epi16(0x7FFE),
		prob=prob_half,
#ifdef ENABLE_CONFIDENCE
		prob_correct=prob_half,
#endif
		limit=_mm_set1_epi32(0x1000000);
	__m128i code[2], start[2], range[2], r2[2];
	__m128i shuf_pack=_mm_set_epi8(15, 14, 11, 10, 7, 6, 3, 2, 13, 12, 9, 8, 5, 4, 1, 0);

	if(!src_start||!src_end||!dst||!count||!bytestride)
		FAIL("abac4_decode(src_start=%p, src_end=%p, dst=%p, count=%d, bytestride=%d)", src_start, src_end, dst, count, bytestride);

	t1=__rdtsc();
	data=(const unsigned char*)src_start;
	buffer=(unsigned char*)dst;
	
	int headercount=1+8;
	if((int)(headercount*sizeof(int))>=(unsigned char*)src_end-data)
		FAIL("File is %d bytes < %d", (int)((unsigned char*)src_end-data), headercount*sizeof(int));
	int magic;
	memcpy(&magic, data, sizeof(int));
	if(magic!=magic_ac0A)
		FAIL("Invalid magic number 0x%08X, expected 0x%08X", magic, magic_ac0A);
	sizes=data+sizeof(int);
	data=sizes+8*sizeof(int);

	size_t offset=0;
	for(int k=0;k<8;++k)
	{
		int size;
		memcpy(&size, sizes+k*sizeof(int), sizeof(int));
		ptr[k]=data+offset;
		offset+=size;
		end[k]=data+offset;
	}
	if((unsigned char*)src_end<end[7])
		FAIL("File is %d bytes, header requires %d bytes", (int)((size_t)src_end-(size_t)src_start), (int)((size_t)end[7]-(size_t)src_start));
	
	start[0]=_mm_setzero_si128();
	start[1]=_mm_setzero_si128();
	range[0]=_mm_set1_epi32(0x7FFFFFFF);//the MSB is the 2's comp sign, all SSE2 cmp instructions are signed
	range[1]=_mm_set1_epi32(0x7FFFFFFF);
	for(int k=0;k<8;++k)
	{
		int reg_idx=k>>2, idx=k&3;
		ptr[k]+=sizeof(int);
		if(end[k]<ptr[k])
			FAIL("Decode OOB: bit %d ptr %d size %d", k, (int)((size_t)ptr[k]-(size_t)src_start), (int)((size_t)src_end-(size_t)src_start));
		code[reg_idx].m128i_i32[idx]=ptr[k][-4]<<24|ptr[k][-3]<<16|ptr[k][-2]<<8|ptr[k][-1];//big endian
	}
	for(int ks=0;ks<count;++ks, dst+=bytestride)
	{
		print_reg(code, 32, "%d code", ks);
		print_reg(start, 32, "%d start", ks);
		print_reg(range, 32, "%d range", ks);
		
		//calculate probability of zero bit
#ifdef ENABLE_CONFIDENCE
		__m128i p0=_mm_sub_epi16(prob, prob_half);

		p0=_mm_slli_epi16(p0, 1);
		prob_correct=_mm_slli_epi16(prob_correct, 1);
		p0=_mm_mulhi_epu16(p0, prob_correct);
		p0=_mm_mulhi_epu16(p0, prob_correct);
		p0=_mm_srli_epi16(p0, 1);
		prob_correct=_mm_srli_epi16(prob_correct, 1);

		p0=_mm_add_epi16(p0, prob_half);
#else
		__m128i p0=prob;
#endif

		//clamp probability			eg: {0, 1, half, FFFF}
		__m128i cmp=_mm_cmplt_epi16(prob_min, p0);//{0, FFFF, FFFF, FFFF}
		p0=_mm_and_si128(p0, cmp);
		cmp=_mm_xor_si128(cmp, ones);
		cmp=_mm_and_si128(cmp, prob_min);
		p0=_mm_or_si128(p0, cmp);

		cmp=_mm_cmplt_epi16(p0, prob_max);//{FFFF, FFFF, FFFF, 0}
		p0=_mm_and_si128(p0, cmp);
		cmp=_mm_xor_si128(cmp, ones);
		cmp=_mm_and_si128(cmp, prob_max);
		p0=_mm_or_si128(p0, cmp);

		//calculate middle (r2)
		__m128i t0=_mm_shuffle_epi8(range[0], shuf_pack);
		__m128i t1=_mm_shuffle_epi8(range[1], shuf_pack);
		__m128i r_lo=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(t0), _mm_castsi128_ps(t1), _MM_SHUFFLE(1, 0, 1, 0)));
		__m128i r_hi=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(t0), _mm_castsi128_ps(t1), _MM_SHUFFLE(3, 2, 3, 2)));
		p0=_mm_slli_epi16(p0, 1);
		__m128i a=_mm_mulhi_epu16(r_lo, p0);
		__m128i b0=_mm_mullo_epi16(r_hi, p0);
		__m128i b1=_mm_mulhi_epu16(r_hi, p0);
		p0=_mm_srli_epi16(p0, 1);

		r2[0]=_mm_unpacklo_epi16(b0, b1);
		r2[1]=_mm_unpackhi_epi16(b0, b1);

		__m128i d0=_mm_unpacklo_epi16(a, _mm_setzero_si128());
		__m128i d1=_mm_unpackhi_epi16(a, _mm_setzero_si128());

		r2[0]=_mm_add_epi32(r2[0], d0);
		r2[1]=_mm_add_epi32(r2[1], d1);

		//avoid degenerate range		r2[k]+=(r2[k]==0)-(r2[k]==range[k]);
		cmp=_mm_cmpeq_epi32(r2[0], _mm_setzero_si128());
		cmp=_mm_and_si128(cmp, ones32);
		__m128i cmp2=_mm_cmpeq_epi32(r2[0], range[0]);
		cmp2=_mm_and_si128(cmp2, ones32);
		cmp=_mm_sub_epi32(cmp, cmp2);
		r2[0]=_mm_add_epi32(r2[0], cmp);
		
		cmp=_mm_cmpeq_epi32(r2[0], _mm_setzero_si128());
		cmp=_mm_and_si128(cmp, ones32);
		cmp2=_mm_cmpeq_epi32(r2[0], range[1]);
		cmp2=_mm_and_si128(cmp2, ones32);
		cmp=_mm_sub_epi32(cmp, cmp2);
		r2[1]=_mm_add_epi32(r2[1], cmp);

		__m128i bitlo=_mm_add_epi32(start[0], r2[0]);
		__m128i bithi=_mm_add_epi32(start[1], r2[1]);
		bitlo=_mm_cmplt_epi32(bitlo, code[0]);
		bithi=_mm_cmplt_epi32(bithi, code[1]);

		//pack low 16-bit halves of each 32-bit int
		bitlo=_mm_shuffle_epi8(bitlo, shuf_pack);
		bithi=_mm_shuffle_epi8(bithi, shuf_pack);
		__m128i bit=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(bitlo), _mm_castsi128_ps(bithi), _MM_SHUFFLE(1, 0, 1, 0)));

		//store symbol
		int mask=_mm_movemask_epi8(bit);
		unsigned char sym=mask>>8&128|mask>>7&64|mask>>6&32|mask>>5&16|mask>>4&8|mask>>3&4|mask>>1&2|mask&1;
#ifdef ENABLE_GRAYCODE
		sym^=sym>>4;
		sym^=sym>>2;
		sym^=sym>>1;
#endif
		*dst=sym;
		bit=_mm_and_si128(bit, prob_min);//leave only the LSB

		print_reg(r2, 32, "%d r2", ks);
		
		print_reg(&prob, 16, "%d prob ", ks);
		print_reg(&p0, 16, "%d p0 ", ks);
		
#ifdef ENABLE_CONFIDENCE
		//update prob_correct			prob_correct=correct<<15|prob_correct>>1;
		cmp=_mm_cmplt_epi16(p0, prob_half);
		cmp=_mm_xor_si128(cmp, ones);
		cmp=_mm_and_si128(cmp, prob_min);
		cmp=_mm_xor_si128(cmp, bit);

		cmp=_mm_slli_epi16(cmp, 14);//set pre-MSB, because setting MSB would negate the 2's comp value
		prob_correct=_mm_srli_epi16(prob_correct, 1);
		prob_correct=_mm_or_si128(prob_correct, cmp);
#endif

		//update prob zero				prob=!bit<<15|prob>>1;
		__m128i bit2=_mm_xor_si128(bit, prob_min);
		prob=_mm_srli_epi16(prob, 1);
		prob=_mm_or_si128(prob, _mm_slli_epi16(bit2, 14));

		print_reg(&prob, 16, "%d prob2 ", ks);
		print_reg(&bit, 16, "%d bit ", ks);

#if 1
		//update range
		bitlo=_mm_unpacklo_epi16(bit, _mm_setzero_si128());
		bithi=_mm_unpackhi_epi16(bit, _mm_setzero_si128());

		//r2[k]+=bit.m128i_i16[k]-!bit.m128i_i16[k];
		r2[0]=_mm_add_epi32(r2[0], bitlo);
		r2[1]=_mm_add_epi32(r2[1], bithi);
		r2[0]=_mm_sub_epi32(r2[0], _mm_sub_epi32(ones32, bitlo));
		r2[1]=_mm_sub_epi32(r2[1], _mm_sub_epi32(ones32, bithi));

		//start[k]+=r2[k]&-bit.m128i_i16[k];
		bitlo=_mm_cmpeq_epi32(bitlo, ones32);
		bithi=_mm_cmpeq_epi32(bithi, ones32);
		t0=_mm_and_si128(bitlo, r2[0]);
		t1=_mm_and_si128(bithi, r2[1]);
		start[0]=_mm_add_epi32(start[0], t0);
		start[1]=_mm_add_epi32(start[1], t1);
		
		//range[k]=((range[k]-r2[k])&-bit.m128i_i16[k])|(r2[k]&-!bit.m128i_i16[k]);
		range[0]=_mm_sub_epi32(range[0], r2[0]);
		range[1]=_mm_sub_epi32(range[1], r2[1]);
		range[0]=_mm_and_si128(range[0], bitlo);
		range[1]=_mm_and_si128(range[1], bithi);
		bitlo=_mm_xor_si128(bitlo, ones);
		bithi=_mm_xor_si128(bithi, ones);
		r2[0]=_mm_and_si128(r2[0], bitlo);
		r2[1]=_mm_and_si128(r2[1], bithi);
		range[0]=_mm_or_si128(range[0], r2[0]);
		range[1]=_mm_or_si128(range[1], r2[1]);

		//renormalize
		unsigned condition;
		for(;;)
		{
			t0=_mm_add_epi32(start[0], range[0]);
			t1=_mm_add_epi32(start[1], range[1]);
			t0=_mm_xor_si128(t0, start[0]);
			t1=_mm_xor_si128(t1, start[1]);
			t0=_mm_cmpeq_epi32(t0, limit);
			t1=_mm_cmpeq_epi32(t1, limit);
			condition=_mm_movemask_epi8(t1)<<16|_mm_movemask_epi8(t0);
			if(!condition)
				break;
			for(int k=0;k<8;++k)
			{
				if(condition>>(k<<2)&15)
				{
					int reg_idx=k>>2, idx=k&3;

					++ptr[k];
					if(end[k]<ptr[k])
						FAIL("Decode OOB: bit %d ptr %d size %d", k, (int)((size_t)ptr[k]-(size_t)src_start), (int)((size_t)src_end-(size_t)src_start));
#ifdef PRINT_SSE2
					printf("RENORM bit %d %08X<-[%02X]\n", k, code[reg_idx].m128i_i32[idx], (int)ptr[k][-1]);
#endif
					code[reg_idx].m128i_i32[idx]<<=8;
					code[reg_idx].m128i_i32[idx]|=ptr[k][-1];

					start[reg_idx].m128i_i32[idx]<<=8;
					range[reg_idx].m128i_i32[idx]=range[reg_idx].m128i_i32[idx]<<8|0xFF;
				}
			}
		}
		t0=_mm_cmplt_epi32(range[0], three);
		t1=_mm_cmplt_epi32(range[1], three);
		condition=_mm_movemask_epi8(t1)<<16|_mm_movemask_epi8(t0);
		if(condition)
		{
			for(int k=0;k<8;++k)
			{
				if(condition>>(k<<2)&15)
				{
					int reg_idx=k>>2, idx=k&3;

					ptr[k]+=sizeof(int);
					if(end[k]<ptr[k])
						FAIL("Decode OOB: bit %d ptr %d size %d", k, (int)((size_t)ptr[k]-(size_t)src_start), (int)((size_t)src_end-(size_t)src_start));
#ifdef PRINT_SSE2
					printf("RENORM bit %d %08X<-[%08X]\n", k, code[reg_idx].m128i_i32[idx], ptr[k][-4]<<24|ptr[k][-3]<<16|ptr[k][-2]<<8|ptr[k][-1]);
#endif
					code[reg_idx].m128i_i32[idx]=ptr[k][-4]<<24|ptr[k][-3]<<16|ptr[k][-2]<<8|ptr[k][-1];//big endian

					start[reg_idx].m128i_i32[idx]=0;
					range[reg_idx].m128i_i32[idx]=0x7FFFFFFF;//because 1=0.9999...
				}
			}
		}
#ifdef PRINT_SSE2
		printf("\n");
#endif
#endif
#if 0
		for(int k=0;k<8;++k)
		{
			if(bit.m128i_i16[k])
			{
				++r2[k];
				start[k]+=r2[k], range[k]-=r2[k];
			}
			else
				range[k]=r2[k]-1;
			
			while((start[k]^(start[k]+(unsigned)range[k]))<0x1000000)//shift-out identical bytes			zpaq 1.10
			{
#ifdef DEBUG_ABAC2
				if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
					printf("range 0x%08X byte-out 0x%02X\n", (int)range, code>>24);
#endif
				++ptr[k];
				if(end[k]<ptr[k])
					FAIL("Decode OOB: bit %d ptr %d size %d", k, (int)((size_t)ptr[k]-(size_t)src_start), (int)((size_t)src_end-(size_t)src_start));
				code[k]=code[k]<<8|ptr[k][-1];
				start[k]<<=8;
				range[k]=range[k]<<8|0xFF;
			}

			if(range[k]<3)
			{
				ptr[k]+=sizeof(int);
				if(end[k]<ptr[k])
					FAIL("Decode OOB: bit %d ptr %d size %d", k, (int)((size_t)ptr[k]-(size_t)src_start), (int)((size_t)src_end-(size_t)src_start));
				code[k]=ptr[k][-4]<<24|ptr[k][-3]<<16|ptr[k][-2]<<8|ptr[k][-1];//big endian

				start[k]=0, range[k]=0xFFFFFFFF;//because 1=0.9999...
			}
		}
#endif
	}
	t2=__rdtsc();

	if(loud)
		printf("AC decode:  %lld cycles, %lf c/byte\n", t2-t1, (double)(t2-t1)/count);
	return ptr[7];
}

int				abac4_encode(const void *src, int symcount, int bitoffset, int bitdepth, int bytestride, ArrayHandle *output, int loud)
{
	DList list;
	const unsigned char *buffer;
	unsigned *header, *sizes, *conf;
	u64 t1, t2;

	if(!src||!symcount||!bitdepth||!bytestride)
		FAIL("abac4_encode(src=%p, symcount=%d, bitdepth=%d, stride=%d)", src, symcount, bitdepth, bytestride);
	t1=__rdtsc();

	dlist_init(&list, 1, 1024, 0);//bitdepth < 127 bit			for max bitdepth of 32: objpernode >= 260
	header=dlist_push_back(&list, 0, (1+bitdepth*2)*sizeof(int));//magic + sizes[bitdepth] + conf[bitdepth]		all must fit in list node
	header[0]=magic_ac04;
	sizes=header+1;//don't worry about alignment here
	conf=sizes+bitdepth;

	buffer=(const unsigned char*)src;

#ifdef AC_MEASURE_PREDICTION
	u64 hitnum=0, hitden=0;//prediction efficiency
#endif

	for(int kp=bitdepth-1;kp>=0;--kp)//bit-plane loop		encode MSB first
	{
		size_t out_planestart=list.nobj;
		int bit_offset=(bitoffset+kp)>>3, bit_shift=(bitoffset+kp)&7;
		int bit_offset2=(kp+1)>>3, bit_shift2=(kp+1)&7;
		unsigned prob=PROB_HALF, prob_correct=PROB_HALF;//cheap weighted average predictor

		u64 hitcount=1;
		
#ifdef PRINT_ABAC
		printf("\n");//
#endif
		for(int kb=0, kb2=0;kb<symcount;++kb, kb2+=bytestride)//analyze bitplane
		{
			int bit=buffer[kb2+bit_offset]>>bit_shift&1;
			int p0=prob-PROB_HALF;
			p0=((long long)p0*prob_correct>>16);
			p0+=PROB_HALF;
			//unsigned p0=PROB_HALF+(long long)(prob-PROB_HALF)*hitcount/(kb+1);
			PROB_CLAMP(p0);
			int correct=bit^(p0>=PROB_HALF);
#ifdef PRINT_ABAC
			printf("%d", bit);//
#endif
			//if(kp==1)
			//	printf("%d", bit);//actual bits
			//	printf("%d", p0<PROB_HALF);//predicted bits
			//	printf("%d", !correct);//prediction error
			hitcount+=correct;
			prob=!bit<<15|prob>>1;
			prob_correct=correct<<15|prob_correct>>1;
		}
#ifdef PRINT_ABAC
		printf("\n");//
#endif
		conf[kp]=(int)hitcount;

		if(hitcount<symcount*MIN_CONF)//incompressible, bypass
		{
#ifdef PRINT_ABAC
			printf("BYPASS\n");//
#endif
			for(int kb=0, kb2=bit_offset;kb+7<symcount;kb+=8)
			{
				unsigned char frag=0;
				frag|= buffer[kb2]>>bit_shift&1, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<1, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<2, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<3, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<4, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<5, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<6, kb2+=bytestride;
				frag|=(buffer[kb2]>>bit_shift&1)<<7, kb2+=bytestride;
				dlist_push_back1(&list, &frag);
			}
		}
		else
		{
			int hitratio_sure=(int)(0x10000*pow((double)hitcount/symcount, 1/BOOST_POWER)), hitratio_notsure=(int)(0x10000*pow((double)hitcount/symcount, BOOST_POWER));
			int hitratio_delta=hitratio_sure-hitratio_notsure;
		//	int hitratio_sure=int(0x10000*cbrt((double)hitcount/symcount)), hitratio_notsure=int(0x10000*(double)hitcount*hitcount*hitcount/((double)symcount*symcount*symcount));
		//	int hitratio_sure=int(0x10000*sqrt((double)hitcount/symcount)), hitratio_notsure=int(0x10000*(double)hitcount*hitcount/((double)symcount*symcount));
			hitcount=(hitcount<<16)/symcount;

			//hitcount=unsigned(((u64)hitcount<<16)/symcount);
			//hitcount=abac2_normalize16(hitcount, logimsize);
			//hitcount*=invimsize;

			prob_correct=prob=PROB_HALF;

#ifdef ABAC2_CONF_MSB_RELATION
			int prevbit0=0;
#endif
			
			unsigned start=0;
			u64 range=0xFFFFFFFF;
			for(int kb=0, kb2=0;kb<symcount;kb2+=bytestride)//bit-pixel loop		http://mattmahoney.net/dc/dce.html#Section_32
			{
				int bit=buffer[kb2+bit_offset]>>bit_shift&1;
#ifdef ABAC2_CONF_MSB_RELATION
				int prevbit=buffer[kb2+bit_offset2]>>bit_shift2&1;
			//	int prevbit=buffer[kb]>>(kp+1)&1;
#endif
				
				if(range<3)
				{
					dlist_push_back1(&list, (char*)&start+3);//big endian
					dlist_push_back1(&list, (char*)&start+2);
					dlist_push_back1(&list, (char*)&start+1);
					dlist_push_back1(&list, &start);

					start=0, range=0xFFFFFFFF;//because 1=0.9999...
				}
				
				//if(kp==7&&kb==5)
				//	int LOL_1=0;
				int p0=prob-PROB_HALF;
				p0=(int)((long long)p0*prob_correct>>16);
				p0=(int)((long long)p0*prob_correct>>16);
				int sure=-(prevbit==prevbit0)&-(kp!=bitdepth-1);
				p0=(int)((long long)p0*(hitratio_notsure+(hitratio_delta&sure))>>16);
				//p0=p0*(prevbit==prevbit0?hitratio_sure:hitratio_notsure)>>16;
				//p0=(long long)p0*hitcount>>16;
				p0+=PROB_HALF;
				//if(prevbit!=prevbit0)
				//	p0=PROB_HALF;
				//	p0=0xFFFF-p0;

				//unsigned p0=PROB_HALF+((long long)(prob-PROB_HALF)*(prevbit==prevbit0?hitratio_sure:hitratio_notsure)>>16);

				//unsigned p0=(long long)(prob-PROB_HALF)*sqrthitcount>>16;
				//if(prevbit==prevbit0)
				//	p0=(long long)p0*hitcount>>16;
				//p0+=PROB_HALF;

				//int confboost=prevbit==prevbit0;
				//confboost-=!confboost;
				//confboost<<=LOG_CONFBOOST;
				//unsigned p0=PROB_HALF+((long long)(prob-PROB_HALF)*(hitcount+confboost)>>16);

			//	unsigned p0=PROB_HALF+(int)((prob-PROB_HALF)*(prevbit==prevbit0?sqrt((double)test_conf[kp]/symcount):(double)test_conf[kp]*test_conf[kp]/((double)symcount*symcount)));
			//	unsigned p0=prevbit==prevbit0?prob:PROB_HALF;
			//	unsigned p0=PROB_HALF+(long long)(prob-PROB_HALF)*test_conf[kp]/symcount;
			//	unsigned p0=PROB_HALF+(long long)(prob-PROB_HALF)*hitcount/(kb+1);
				PROB_CLAMP(p0);
				unsigned r2=(unsigned)(range*p0>>16);
				r2+=(r2==0)-(r2==range);
#ifdef DEBUG_ABAC2
				if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
					printf("%6d %6d %d %08X+%08X %08X %08X\n", kp, kb, bit, start, (int)range, r2, start+r2);
#endif

				int correct=bit^(p0>=PROB_HALF);
			//	hitcount+=correct;
				prob=!bit<<15|prob>>1;
				prob_correct=correct<<15|prob_correct>>1;
#ifdef ABAC2_CONF_MSB_RELATION
				prevbit0=prevbit;
#endif
#ifdef AC_MEASURE_PREDICTION
				hitnum+=correct, ++hitden;
#endif
				unsigned start0=start;
#ifdef PRINT_ABAC
				u64 r0=range;
#endif
				if(bit)
				{
					++r2;
					start+=r2, range-=r2;
				}
				//	start=middle+1;
				else
					range=r2-1;
				//	end=middle-1;
#ifdef PRINT_ABAC
				if(kp==0)
					printf("%d %d bit %d p0 %04X %08X+%08X -> %08X+%08X\n", kp, kb, bit, p0, start0, r0, start, range);
#endif
				if(start<start0)//
				{
					FAIL("AC OVERFLOW: start = %08X -> %08X, r2 = %08X", start0, start, r2);
					//printf("OVERFLOW\nstart = %08X -> %08X, r2 = %08X", start0, start, r2);
					//int k=0;
					//scanf_s("%d", &k);
				}
				++kb;
				
				while((start^(start+(unsigned)range))<0x1000000)//most significant byte has stabilized			zpaq 1.10
				{
#ifdef DEBUG_ABAC2
					if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
						printf("range 0x%08X byte-out 0x%02X\n", (int)range, start>>24);
#endif
					dlist_push_back1(&list, (char*)&start+3);
					start<<=8;
					range=range<<8|0xFF;
				}
			}
			//start+=range>>1;
			dlist_push_back1(&list, (char*)&start+3);//big endian
			dlist_push_back1(&list, (char*)&start+2);
			dlist_push_back1(&list, (char*)&start+1);
			dlist_push_back1(&list, &start);
		}
		if(loud)
			printf("bit %d: conf = %6d / %6d = %lf%%\n", kp, conf[kp]-1, symcount, 100.*(conf[kp]-1)/symcount);
		sizes[kp]=(int)(list.nobj-out_planestart);
	}
	size_t offset0=*output?output[0]->count*output[0]->esize:0;
	dlist_appendtoarray(&list, output);
	t2=__rdtsc();

	if(loud)
	{
		int original_bitsize=symcount*bitdepth, compressed_bitsize=(int)(list.nobj)<<3;
		printf("AC encode:  %lld cycles, %lf c/byte\n", t2-t1, (double)(t2-t1)/(original_bitsize>>3));
		printf("Size: %d -> %d, ratio: %lf, %lf bpp\n", original_bitsize>>3, compressed_bitsize>>3, (double)original_bitsize/compressed_bitsize, (double)compressed_bitsize/symcount);
#ifdef AC_MEASURE_PREDICTION
		printf("Predicted: %6lld / %6lld = %lf%%\n", hitnum, hitden, 100.*hitnum/hitden);
#endif
		printf("Bit\tbytes\tratio,\tbytes/bitplane = %d\n", (symcount+7)>>3);
		for(int k=0;k<bitdepth;++k)
			printf("%2d\t%5d\t%lf\n", bitdepth-1-k, sizes[k], (double)symcount/(sizes[k]<<3));
		
		printf("Preview:\n");
		int kprint=list.nobj<200?(int)list.nobj:200;
		for(int k=0;k<kprint;++k)
			printf("%02X-", output[0]->data[offset0+k]&0xFF);
		printf("\n");
	}
	dlist_clear(&list);
	return 1;
}
const void*		abac4_decode(const void *in_start, const void *in_end, void *dst, int imsize, int bitoffset, int bitdepth, int bytestride, int loud)//set the dst buffer to zero
{
	const unsigned char *sizes, *conf, *data;
	unsigned char *buffer;
	u64 t1, t2;

	if(!in_start||!in_end||!dst||!imsize||!bitdepth||!bytestride)
		FAIL("abac4_decode(in_start=%p, in_end=%p, dst=%p, imsize=%d, bitdepth=%d, stride=%d)", in_start, in_end, dst, imsize, bitdepth, bytestride);
	t1=__rdtsc();
	data=(const unsigned char*)in_start;
	buffer=(unsigned char*)dst;
	//memset(buffer, 0, imsize*bytestride);

	int headercount=1+(bitdepth<<1);
	if((int)(headercount*sizeof(int))>=(unsigned char*)in_end-data)
		FAIL("File is %d bytes < %d", (int)((unsigned char*)in_end-data), headercount*sizeof(int));
	int magic;
	memcpy(&magic, data, sizeof(int));
	if(magic!=magic_ac04)
		FAIL("Invalid magic number 0x%08X, expected 0x%08X", magic, magic_ac04);
	sizes=data+sizeof(int);
	conf=sizes+bitdepth*sizeof(int);
	data=conf+bitdepth*sizeof(int);
	
	for(int kp=bitdepth-1;kp>=0;--kp)//bit-plane loop
	{
		int kp2=kp+bitoffset;
		int bit_offset=kp2>>3, bit_shift=kp2&7;
		int bit_offset2=(kp2+1)>>3, bit_shift2=(kp2+1)&7;
		int ncodes;
		memcpy(&ncodes, sizes+kp*sizeof(int), sizeof(int));
		
		unsigned prob=PROB_HALF, prob_correct=PROB_HALF;
#if 1
		u64 hitcount=0;
		memcpy(&hitcount, conf+kp*sizeof(int), sizeof(int));
		if(hitcount<imsize*MIN_CONF)
		{
			for(int kb=0, kb2=0;kb<imsize;++kb, kb2+=bytestride)
			{
				int byte_idx=kb>>3, bit_idx=kb&7;
				const unsigned char *b=data+byte_idx;
				if((const unsigned char*)in_end<b)
				{
					FAIL("Decode error (bypass): ptr %p > in_end %p", b, in_end);
					return 0;
				}
				int bit=*b>>bit_idx&1;
				buffer[kb2+bit_offset]|=bit<<bit_shift;
			//	buffer[kb]|=bit<<kp;
#ifdef ENABLE_GUIDE
				int original_bit=((unsigned char*)guide)[kb2+bit_offset]>>bit_shift&1;
				if(bit!=original_bit)
					FAIL("Decode error (bypass): data %d expected %d got %d", kp, original_bit, bit);
#endif
			}
			data+=ncodes;
			//cusize+=ncodes;
			continue;
		}
#ifdef ABAC2_CONF_MSB_RELATION
		int prevbit0=0;
#endif
		int hitratio_sure=(int)(0x10000*pow((double)hitcount/imsize, 1/BOOST_POWER)), hitratio_notsure=(int)(0x10000*pow((double)hitcount/imsize, BOOST_POWER));
		int hitratio_delta=hitratio_sure-hitratio_notsure;
		//int hitratio_sure=int(0x10000*cbrt((double)hitcount/imsize)), hitratio_notsure=int(0x10000*(double)hitcount*hitcount*hitcount/((double)imsize*imsize*imsize));
		//int hitratio_sure=int(0x10000*sqrt((double)hitcount/imsize)), hitratio_notsure=int(0x10000*(double)hitcount*hitcount/((double)imsize*imsize));
		hitcount=(hitcount<<16)/imsize;
		//hitcount=unsigned(((u64)hitcount<<16)/imsize);
		//hitcount=abac2_normalize16(hitcount, logimsize);
		//hitcount*=invimsize;
#endif

		unsigned code, start=0;
		u64 range=0xFFFFFFFF;

		code=data[0]<<24|data[1]<<16|data[2]<<8|data[3];
		data+=sizeof(int);
		for(int kb=0, kb2=0;kb<imsize;kb2+=bytestride)//bit-pixel loop
		{
			if(range<3)
			{
				code=data[0]<<24|data[1]<<16|data[2]<<8|data[3];
				data+=sizeof(int);
				start=0, range=0xFFFFFFFF;//because 1=0.9999...
			}
#ifdef ABAC2_CONF_MSB_RELATION
			int prevbit=0;
			if(kp+1<bitdepth)
				prevbit=buffer[kb2+bit_offset2]>>bit_shift2&1;
		//	int prevbit=buffer[kb]>>(kp+1)&1;
#endif
			//if(kp==7&&kb==5)
			//	int LOL_1=0;
			int p0=prob-PROB_HALF;
			p0=(long long)p0*prob_correct>>16;
			p0=(long long)p0*prob_correct>>16;
			int sure=-(prevbit==prevbit0)&-(kp!=bitdepth-1);
			p0=(long long)p0*(hitratio_notsure+(hitratio_delta&sure))>>16;
			//p0=p0*(prevbit==prevbit0?hitratio_sure:hitratio_notsure)>>16;
			//p0=(long long)p0*hitcount>>16;
			p0+=PROB_HALF;
			//if(prevbit!=prevbit0)
			//	p0=PROB_HALF;
			//	p0=0xFFFF-p0;

			//unsigned p0=PROB_HALF+((long long)(prob-PROB_HALF)*(prevbit==prevbit0?hitratio_sure:hitratio_notsure)>>16);

			//unsigned p0=(long long)(prob-PROB_HALF)*sqrthitcount>>16;
			//if(prevbit==prevbit0)
			//	p0=(long long)p0*hitcount>>16;
			//p0+=PROB_HALF;

			//int confboost=prevbit==prevbit0;
			//confboost-=!confboost;
			//confboost<<=LOG_CONFBOOST;
			//unsigned p0=PROB_HALF+((long long)(prob-PROB_HALF)*(hitcount+confboost)>>16);

		//	unsigned p0=PROB_HALF+(int)((prob-PROB_HALF)*(prevbit==prevbit0?sqrt((double)test_conf[kp]/imsize):(double)test_conf[kp]*test_conf[kp]/((double)imsize*imsize)));
		//	unsigned p0=prevbit==prevbit0?prob:PROB_HALF;
		//	unsigned p0=PROB_HALF+(long long)(prob-PROB_HALF)*test_conf[kp]/imsize;
		//	unsigned p0=PROB_HALF+(long long)(prob-PROB_HALF)*hitcount/(kb+1);
			PROB_CLAMP(p0);
			unsigned r2=(unsigned)(range*p0>>16);
			r2+=(r2==0)-(r2==range);
			unsigned middle=start+r2;
			int bit=code>middle;
#ifdef DEBUG_ABAC2
			if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
				printf("%6d %6d %d %08X+%08X %08X %08X %08X\n", kp, kb, bit, start, (int)range, r2, middle, code);
#endif
			
			int correct=bit^(p0>=PROB_HALF);
		//	hitcount+=correct;
			prob=!bit<<15|prob>>1;
			prob_correct=correct<<15|prob_correct>>1;
#ifdef ABAC2_CONF_MSB_RELATION
			prevbit0=prevbit;
#endif
			
#ifdef PRINT_ABAC
			auto start0=start;
			auto r0=range;
#endif
			if(bit)
			{
				++r2;
				start+=r2, range-=r2;
			}
			//	start=middle+1;
			else
				range=r2-1;
			//	end=middle-1;
#ifdef PRINT_ABAC
			if(kp==0)
				printf("%d %d bit %d p0 %04X %08X+%08X -> %08X+%08X c %08X\n", kp, kb, bit, p0, start0, r0, start, range, code);
#endif
			
			buffer[kb2+bit_offset]|=bit<<bit_shift;
		//	buffer[kb]|=bit<<kp;
#ifdef ENABLE_GUIDE
			int original_bit=((unsigned char*)guide)[kb2+bit_offset]>>bit_shift&1;
			if(bit!=original_bit)
				FAIL("Decode error (AC): data %d expected %d got %d", kp, original_bit, bit);
#endif
			++kb;
			
			while((start^(start+(unsigned)range))<0x1000000)//shift-out identical bytes			zpaq 1.10
			{
#ifdef DEBUG_ABAC2
				if(kp==examined_plane&&kb>=examined_start&&kb<examined_end)
					printf("range 0x%08X byte-out 0x%02X\n", (int)range, code>>24);
#endif
				code=code<<8|*data;
				++data;
				start<<=8;
				range=range<<8|0xFF;
			}
		}
		//data+=ncodes;
		//cusize+=ncodes;
	}
	t2=__rdtsc();

	if(loud)
		printf("AC decode:  %lld cycles, %lf c/byte\n", t2-t1, (double)(t2-t1)/(imsize*bitdepth>>3));
	return data;
}