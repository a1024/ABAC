#define _USE_MATH_DEFINES
#include<math.h>
#include<immintrin.h>
#include"awm_util.h"

//color transforms

//YCoCg > YCbCr	https://en.wikipedia.org/wiki/YCoCg
void			YCoCg_f32_fwd(const int *image, size_t size, float *bufY, float *bufCo, float *bufCg)
{
	ptrdiff_t k;

	for(k=0;k<(ptrdiff_t)size;++k)
	{
		unsigned char *p;
		int R, G, B, Y, Co, Cg;

		p=(unsigned char*)(image+k);
		R=p[0], G=p[1], B=p[2];

		Co=R-B;
		Y=B+(Co>>1);
		Cg=G-Y;
		Y+=(Cg>>1);

//#define MINMAX(MIN, MAX, VAR)	if(MIN>VAR)MIN=VAR; if(MAX<VAR)MAX=VAR;
//		MINMAX(minY, maxY, Y)
//		MINMAX(minCo, maxCo, Co)
//		MINMAX(minCg, maxCg, Cg)
//#undef	MINMAX

		//if(Y<0||Co<0||Cg<0)
		//	Y=Y;

		bufY[k]=(float)Y, bufCo[k]=(float)Co, bufCg[k]=(float)Cg;//u8+s9+s9 bit
	}
}
void			YCoCg_f32_inv(int *image, size_t size, const float *bufY, const float *bufCo, const float *bufCg)
{
	ptrdiff_t k;

	for(k=0;k<(ptrdiff_t)size;++k)
	{
		unsigned char *p;
		float R, G, B, tmp;

		p=(unsigned char*)(image+k);
		tmp=bufY[k]-bufCg[k]/2;
		G=bufCg[k]+tmp;
		B=tmp-bufCo[k]/2;
		R=B+bufCo[k];
#define CLAMP(X)	X=X>255?255:(X<0?0:X)
		CLAMP(R);
		CLAMP(G);
		CLAMP(B);
#undef	CLAMP
		p[0]=(unsigned char)R, p[1]=(unsigned char)G, p[2]=(unsigned char)B, p[3]=0xFF;
	}
}

//https://en.wikipedia.org/wiki/YCoCg#The_lifting-based_YCoCg-R_variation
void			YCoCg_i32_fwd(int *image, int count)
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
void			YCoCg_i32_inv(int *image, int count)
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


//spatial transforms


//DCT-II/III naive, arbitrary size
void init_fwdDCT(float *matrix, int lgdim)
{
	int dim=1<<lgdim;
	float scale=sqrtf(2.f/dim), freq=(float)M_PI/dim;
	for(int ky=0;ky<dim;++ky)
	{
		float *row=matrix+dim*ky;
		for(int kx=0;kx<dim;++kx)
			row[kx]=scale*cosf(freq*(kx+0.5f)*ky);
	}
}
void init_invDCT(float *matrix, int lgdim)
{
	int dim=1<<lgdim;
	float scale=sqrtf(2.f/dim), freq=(float)M_PI/dim;
	for(int ky=0;ky<dim;++ky)
	{
		float *row=matrix+dim*ky;
		for(int kx=0;kx<dim;++kx)
		{
			if(!kx)
				row[kx]=scale*0.5f;
			else
				row[kx]=scale*cosf(freq*kx*(ky+0.5f));
		}
	}
}
void apply_DCT_1D(const float *matrix, float *data, int count, int stride, float *temp)
{
	const float *row=matrix;
	for(int ks=0, kd=0;kd<count;ks+=stride, ++kd)
		temp[kd]=data[ks];
	for(int ks=0, kd=0;ks<count;++ks, kd+=stride, row+=count)
	{
		float sum=0;
		for(int k2=0;k2<count;++k2)
			sum+=row[k2]*temp[k2];
		data[kd]=sum;
	}
}
void apply_DCT_2D(const float *hmatrix, const float *vmatrix, float *data, int bw, int bh, float *temp)
{
	int size=bw*bh;
	for(int ky=0;ky<size;ky+=bw)
		apply_DCT_1D(hmatrix, data+ky, bw, 1, temp);
	for(int kx=0;kx<bw;++kx)
		apply_DCT_1D(vmatrix, data+kx, bh, bw, temp);
}


//DCT-II/III fast, P.O.T. size
#define BIT_REVERSE(N, NBITS)\
	N=N<<1&0xAAAA|N>>1&0x5555;\
	N=N<<2&0xCCCC|N>>2&0x3333;\
	N=N<<4&0xF0F0|N>>4&0x0F0F;\
	N=N<<8&0xFF00|N>>8&0x00FF;\
	N>>=16-NBITS;
void fft_dit_ps(float *re, float *im, int lgsize, int forward)//DIT vs DIF: simply reverse the order of stage matrices along with the permutation, and transpose the stages
{
	int size=1<<lgsize, sign=!forward-(forward!=0);
	
	for(int k=0;k<size;k+=2)//separate mul-free sum & difference (m==2)
	{
		float
			re0=re[k],
			im0=im[k],
			re1=re[k+1],
			im1=im[k+1];

		re[k]=re0+re1;
		im[k]=im0+im1;
		re[k+1]=re0-re1;
		im[k+1]=im0-im1;
	}
	for(int m=4;m<=size;m<<=1)
	{
		int m2=m>>1;
		float angle=(float)((2*M_PI)*sign/m);
		float rewm=(float)cos(angle), imwm=(float)sin(angle);//exp(-+ 2*pi*i*(N/m)/N)
		for(int k=0;k<size;k+=m)
		{
			float rew=1, imw=0;

			for(int j=0;;)
			{
				int idx=k+j;
				
				float
					re0=re[idx],
					im0=im[idx],
					re1=re[idx+m2],
					im1=im[idx+m2],
					re2=rew*re1-imw*im1,
					im2=rew*im1+imw*re1;

				re[idx]=re0+re2;
				im[idx]=im0+im2;
				re[idx+m2]=re0-re2;
				im[idx+m2]=im0-im2;

				++j;
				if(j>=m2)
					break;

				re0=rewm*rew-imwm*imw;//w*=wm
				im0=rewm*imw+imwm*rew;
				rew=re0;
				imw=im0;
			}
		}
	}
}
void fft_dif_ps(float *re, float *im, int lgsize, int forward)
{
	int size=1<<lgsize, sign=!forward-(forward!=0);
	
	for(int m=size;m>2;m>>=1)
	{
		int m2=m>>1;
		float angle=(float)((2*M_PI)*sign/m);
		float rewm=(float)cos(angle), imwm=(float)sin(angle);//exp(-+ 2*pi*i*(N/m)/N)
		for(int k=0;k<size;k+=m)
		{
			float rew=1, imw=0;
			for(int j=0;;)
			{
				int idx=k+j;
				float
					re0=re[idx],
					im0=im[idx],
					re1=re[idx+m2],
					im1=im[idx+m2];

				re[idx]=re0+re1;
				im[idx]=im0+im1;

				re0-=re1;
				im0-=im1;
				
				re[idx+m2]=rew*re0-imw*im0;
				im[idx+m2]=rew*im0+imw*re0;
				
				++j;
				if(j>=m2)
					break;
				
				re0=rewm*rew-imwm*imw;//w*=wm
				im0=rewm*imw+imwm*rew;
				rew=re0;
				imw=im0;
			}
		}
	}
	for(int k=0;k<size;k+=2)//separate mul-free sum & difference (m==2)
	{
		float
			re0=re[k],
			im0=im[k],
			re1=re[k+1],
			im1=im[k+1];

		re[k]=re0+re1;
		im[k]=im0+im1;
		re[k+1]=re0-re1;
		im[k+1]=im0-im1;
	}
}

//A fast cosine transform (FCT) in one and two dimensions - page 103
typedef struct FCT_ParamsStruct
{
	int lgsize;
	int *fctp;
	float *re, *im, *rew4N, *imw4N;
} FCT_Params;
void gen_FCT(int lgsize, FCT_Params *p)//up to size 65536
{
	int size=1<<lgsize, sh=16-lgsize, mask=(1<<lgsize)-1;
	p->lgsize=lgsize;
	p->fctp=(int*)malloc(size*sizeof(int));
	p->re=(float*)malloc(size*sizeof(float));
	p->im=(float*)malloc(size*sizeof(float));
	p->rew4N=(float*)malloc(size*sizeof(float));
	p->imw4N=(float*)malloc(size*sizeof(float));
	for(int k=0;k<size;++k)//FCT permutation: f(n) = (bit_reverse(n>>1)^-(n&1)) & ((1<<lgsize)-1)
	{
		unsigned val=k>>1;
		BIT_REVERSE(val, lgsize)
		val^=-(k&1);
		val&=mask;
		p->fctp[k]=val;
	}
	float invsqrtN=(float)(1/sqrt(size));
	for(int k=0;k<size;++k)
	{
		float temp=-2*M_PI*k/(size<<2);
		p->rew4N[k]=invsqrtN*cos(temp);
		p->imw4N[k]=invsqrtN*sin(temp);
	}
}
void free_FCT(FCT_Params *p)
{
	free(p->fctp);
	free(p->re);
	free(p->im);
	free(p->rew4N);
	free(p->imw4N);
	memset(p, 0, sizeof(*p));
}
void apply_FCT_1D_ps(FCT_Params *p, float *data, int stride)//data, re, im, rew, imw are all of size (1<<lgsize), result is in data buffer
{
	int size=1<<p->lgsize;

	for(int ks=0, kd=0;kd<size;ks+=stride, ++kd)
		p->re[p->fctp[kd]]=data[ks];
	memset(p->im, 0, size*sizeof(float));

	fft_dit_ps(p->re, p->im, p->lgsize, 1);//V[k] = sum n=0 to N-1: v[n] WN^nk

	for(int ks=0, kd=0;ks<size;++ks, kd+=stride)//C[k] = 2 Re{W4N^k V[k]}
	{
		float r=p->re[ks], i=p->im[ks];
		//data[kd]=p->rew4N[ks]*r-p->imw4N[ks]*i;
		data[kd]=2*(p->rew4N[ks]*r-p->imw4N[ks]*i);
		//imaginary part would be 2(imw4N[ks]*r+rew4N[ks]*i) = data[(N*stride-kd)%N*stride] or data[(N-1)*stride-kd] ?
	}
}
void apply_IFCT_1D_ps(FCT_Params *p, float *data, int stride)//data, re, im, rew, imw are all of size (1<<lgsize), result is in data buffer
{
	int size=1<<p->lgsize, offset=size*stride;

	for(int ks=0, kd=0;kd<size;ks+=stride, ++kd)
	{
		float r=data[ks], i=ks?-data[offset-ks]:0;//V[k] = (1/2) W4N^-k (C[k] - i*C[N-k])
		//p->re[kd]= p->rew4N[kd]*r-p->imw4N[kd]*i;
		//p->im[kd]=-p->rew4N[kd]*i-p->imw4N[kd]*r;
		p->re[kd]=0.5f*(p->rew4N[kd]*r+p->imw4N[kd]*i);	//(1/2) (cos(2pik/4N) *  C[k]   - -sin(2pik/4N) *(-C[N-k]))
		p->im[kd]=0.5f*(p->rew4N[kd]*i-p->imw4N[kd]*r);	//(1/2) (cos(2pik/4N) *(-C[N-k])+ -sin(2pik/4N) *  C[k])
	}

	fft_dif_ps(p->re, p->im, p->lgsize, 0);

	for(int ks=0, kd=0;kd<size;ks+=stride, ++kd)
		data[p->fctp[kd]]=p->re[ks];
}


//DCT-II/III 8x8, optimized, fixed precision 16.16
void DCT8_fix_fwd(int *data)
{
	static const long long coeff[]=//fixed 32.32
	{
		0x061F78A9BLL,//0.382683432365089771728460,
		0x08A8BD3DFLL,//0.541196100146196984399723,
		0x0B504F334LL,//0.707106781186547524400844,
		0x14E7AE914LL,//1.306562964876376527856643,
	};
	const int
		a0=data[0]+data[7],
		a1=data[1]+data[6],
		a2=data[2]+data[5],
		a3=data[3]+data[4],
		a4=data[3]-data[4],
		a5=data[2]-data[5],
		a6=data[1]-data[6],
		a7=data[0]-data[7];

	const int
		e0=a4+a5,
		e1=(a5+a6)*coeff[2]>>32,
		e2=a6+a7,
		e3=a0+a3,
		e4=a1+a2,
		e6=a0-a3,
		e5=a1-a2;
	
	const int
		c0=(e2-e0)*coeff[0]>>32,
		d3=a7+e1,
		d4=a7-e1;

	const int
		d0=(e5+e6)*coeff[2]>>32,
		d1=(e0*coeff[1]>>32)-c0,
		d2=(e2*coeff[3]>>32)-c0;
	
	data[0]=(e3+e4)>>3;
	data[1]=(d3+d2)>>3;
	data[2]=(d0+e6)>>3;
	data[3]=(d4-d1)>>3;
	data[4]=(e3-e4)>>3;
	data[5]=(d1+d4)>>3;
	data[6]=(e6-d0)>>3;
	data[7]=(d3-d2)>>3;
}
void DCT8_fix_inv(int *data)
{
	static const long long coeff[]=//fixed 32.32
	{
		 0x0C3EF1535,// 0.76536686473017954345692,
		-0x29CF5D229,//-2.613125929752753055713286,
		 0x11517A7BE,// 1.082392200292393968799446,
		 0x16A09E668,// 1.41421356237309504880168872421,
	};
	const int
		a0=data[0]+data[4],
		a1=data[2]+data[6],
		a2=data[2]-data[6],
		a3=data[0]-data[4],
		a4=data[1]+data[7],
		a5=data[5]+data[3],
		a7=data[5]-data[3],
		a6=data[1]-data[7];

	const int
		m0=a0+a1,
		m1=a0-a1,
		m2=a4+a5,
		m3=a4-a5,
		m4=(a7-a6)*coeff[0]>>32,
		m5=(a2*coeff[3]>>32)-a1;

	const int
		n0=a3+m5,
		n1=a3-m5,
		n2=(a7*coeff[1]>>32)+m4,
		n3=(a6*coeff[2]>>32)-m4;
	
	const int o0=n3-m2;

	const int p0=(m3*coeff[3]>>32)-o0;
	
	const int q0=-p0-n2;
	
	data[0]=m0+m2;
	data[1]=n0+o0;
	data[2]=n1+p0;
	data[3]=m1+q0;
	data[4]=m1-q0;
	data[5]=n1-p0;
	data[6]=n0-o0;
	data[7]=m0-m2;
}


//DWT2 (lifting scheme)
void dwt2_1d_scale(float *even, float *odd, int halfsize, float ce, float co)
{
	int xrem=halfsize&3, xround=halfsize-xrem;
	__m128 factor=_mm_set1_ps(ce);
	for(int k=0;k<xround;k+=4)//x[2k] *= ce
	{
		__m128 val=_mm_loadu_ps(even+k);
		val=_mm_mul_ps(val, factor);
		_mm_storeu_ps(even+k, val);
	}
	for(int k=xround;k<halfsize;++k)
		even[k]*=ce;
	factor=_mm_set1_ps(co);
	for(int k=0;k<xround;k+=4)//x[2k+1] *= co
	{
		__m128 val=_mm_loadu_ps(odd+k);
		val=_mm_mul_ps(val, factor);
		_mm_storeu_ps(odd+k, val);
	}
	for(int k=xround;k<halfsize;++k)
		odd[k]*=co;
}
void dwt2_1d_predict_next(float *even, float *odd, int halfsize, float z10)//symmetric padding at boundary
{
	--halfsize;
	int xround=halfsize-(halfsize&3);
	__m128 f0=_mm_set1_ps(z10);
	for(int k=0;k<xround;k+=4)//x[2k+1] += z10*(x[2k]+x[k2+2])
	{
		__m128 vo=_mm_loadu_ps(odd+k);
		__m128 ve0=_mm_loadu_ps(even+k);
		__m128 ve1=_mm_loadu_ps(even+k+1);
		ve0=_mm_add_ps(ve0, ve1);
		ve0=_mm_mul_ps(ve0, f0);
		ve0=_mm_add_ps(ve0, vo);
		_mm_storeu_ps(odd+k, ve0);
	}
	for(int k=xround;k<halfsize;++k)
		odd[k]+=z10*(even[k]+even[k+1]);
	odd[halfsize]+=z10*(even[halfsize]+even[halfsize-1]);
}
void dwt2_1d_update_prev(float *even, float *odd, int halfsize, float z0m1)//symmetric padding at boundary
{
	__m128 f0=_mm_set1_ps(z0m1);
	int k=halfsize-5;
	for(;k>0;k-=4)//x[2k] += z0m1*(x[2k-1]+x[k2+1])
	{
		__m128 ve=_mm_loadu_ps(even+k);
		__m128 vo0=_mm_loadu_ps(odd+k);
		__m128 vo1=_mm_loadu_ps(odd+k-1);
		vo0=_mm_add_ps(vo0, vo1);
		vo0=_mm_mul_ps(vo0, f0);
		vo0=_mm_add_ps(vo0, ve);
		_mm_storeu_ps(even+k, vo0);
	}
	k+=4;
	for(;k>0;--k)
		even[k]+=z0m1*(odd[k]+odd[k-1]);
	even[0]+=z0m1*(odd[0]+odd[1]);
}
const float dwt97_coeff[]=
{
	-1.58613434342059f,	//alpha
	-0.0529801185729f,	//beta
	0.8829110755309f,	//gamma
	0.4435068520439f,	//delta
	1.1496043988602f,	//zeta
};
void dwt2_1d_9_7(float *even, float *odd, int halfsize)
{
	//'factring wavelet transforms into lifting steps' - page 19
	//'a wavelet tour of signal processing - the sparse way' - page 376
	dwt2_1d_scale(even, odd, halfsize, dwt97_coeff[4], 1/dwt97_coeff[4]);
	dwt2_1d_predict_next(even, odd, halfsize, dwt97_coeff[3]);
	dwt2_1d_update_prev(even, odd, halfsize, dwt97_coeff[2]);
	dwt2_1d_predict_next(even, odd, halfsize, dwt97_coeff[1]);
	dwt2_1d_update_prev(even, odd, halfsize, dwt97_coeff[0]);
}
void dwt2_1d_9_7_inv(float *even, float *odd, int halfsize)
{
	dwt2_1d_update_prev(even, odd, halfsize, -dwt97_coeff[0]);
	dwt2_1d_predict_next(even, odd, halfsize, -dwt97_coeff[1]);
	dwt2_1d_update_prev(even, odd, halfsize, -dwt97_coeff[2]);
	dwt2_1d_predict_next(even, odd, halfsize, -dwt97_coeff[3]);
	dwt2_1d_scale(even, odd, halfsize, 1/dwt97_coeff[4], dwt97_coeff[4]);
}
void dwt2_1d(float *buffer, int size, int stride, float *b2)//size is even
{
	int halfsize=size>>1;
	float *even=b2+halfsize, *odd=b2;
	
	for(int k=0, ks=0;k<halfsize;++k, ks+=stride<<1)//lazy wavelet: split into even (high frequency) & odd (low frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}

	dwt2_1d_9_7(even, odd, halfsize);

	for(int k=0, ks=0;k<size;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt2_1d_inv(float *buffer, int size, int stride, float *b2)//size is even
{
	int halfsize=size>>1;
	float *even=b2+halfsize, *odd=b2;

	for(int k=0, ks=0;k<size;++k, ks+=stride)
		b2[k]=buffer[ks];

	dwt2_1d_9_7_inv(even, odd, halfsize);

	for(int k=0, ks=0;k<halfsize;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
}
void dwt2_2d_fwd(float *buffer, int bw, int bh, int nstages)
{
	int tsize=maximum(bw, bh);
	float *temp=(float*)malloc(tsize*sizeof(float));
	for(int w2=bw, h2=bh, it=0;w2>=3&&h2>=3&&(!nstages||it<nstages);++it)
	{
		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt2_1d(buffer+bw*ky, w2, 1, temp);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt2_1d(buffer+kx, h2, bw, temp);

		w2-=w2>>1;//w=ceil(w/2)
		h2-=h2>>1;//h=ceil(h/2)
	}
	free(temp);
}
void dwt2_2d_inv(float *buffer, int bw, int bh, int nstages)
{
	int lw=floor_log2(bw), lh=floor_log2(bh);
	short *sizes=(short*)malloc((minimum(lw, lh)<<1)*sizeof(short));
	int nsizes=0;
	for(int w2=bw, h2=bh;w2>=3&&h2>=3&&(!nstages||nsizes<nstages);++nsizes)//calculate dimensions of each stage
	{
		sizes[nsizes<<1]=w2;
		sizes[nsizes<<1|1]=h2;
		w2-=w2>>1;//w=ceil(w/2)
		h2-=h2>>1;//h=ceil(h/2)
	}

	int tsize=maximum(bw, bh);
	float *temp=(float*)malloc(tsize*sizeof(float));
	for(int it=nsizes-1;it>=0;--it)
	{
		int w2=sizes[it<<1], h2=sizes[(it<<1)+1];

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt2_1d_inv(buffer+kx, h2, bw, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt2_1d_inv(buffer+bw*ky, w2, 1, temp);
	}
	free(sizes);
	free(temp);
}