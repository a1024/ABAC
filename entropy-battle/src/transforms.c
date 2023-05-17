#include"battle.h"
#include<stdio.h>//for debugging
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#ifdef __AVX__
#include<immintrin.h>
#else
#include<tmmintrin.h>
#endif
static const char file[]=__FILE__;

static void print_sbuf(short *buf, int bw, int bh, int kc, int bytestride)
{
	static int call=0;
	printf("call %d ch %d\n", call, kc);
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
			printf("%d\t", buf[bytestride*(bw*ky+kx)+kc]);
		printf("\n");
	}
	printf("\n");
	++call;
}
static void print_fbuf(float *buf, int bw, int bh, int nplaces)
{
	static int call=0;
	printf("call %d\n", call);
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
			printf("%.*f\t", nplaces, buf[bw*ky+kx]);
		printf("\n");
	}
	printf("\n");
	++call;
}
void compare_bufs_uint8(unsigned char *b1, unsigned char *b0, int iw, int ih, int symbytes, int bytestride, const char *name, int backward)
{
	ptrdiff_t len=(ptrdiff_t)bytestride*iw*ih;
	int inc=bytestride*(1-(backward<<1));
	for(ptrdiff_t k=backward?len-bytestride:0;k>=0&&k<len;k+=inc)
	{
		//{//
		//	ptrdiff_t idx=k>>2, kx=idx%1920, ky=idx/1920;
		//	if(kx==5&&ky==1)
		//		kx=5;
		//}//
		if(memcmp(b1+k, b0+k, symbytes))
		{
			ptrdiff_t idx=k/bytestride, kx=idx%iw, ky=idx/iw;
			printf("%s error XY (%5lld, %5lld) / %5d x %5d  b1 != b0\n", name, kx, ky, iw, ih);
			for(int kc=0;kc<symbytes;++kc)
				printf("C%d  0x%02X != 0x%02X    %d != %d\n", kc, (unsigned)b1[k+kc], (unsigned)b0[k+kc], (unsigned)b1[k+kc], (unsigned)b0[k+kc]);
			//if(backward)
			//	printf("%s error at %d - %d: 0x%02X != 0x%02X\n", name, (int)len-1, (int)(len-1-k), b1[k], b0[k]);
			//else
			//	printf("%s error at %d: 0x%02X != 0x%02X\n", name, (int)k, b1[k], b0[k]);
			return;
		}
	}
	printf("%s:\tSUCCESS\n", name);
}
void compare_bufs_ps(float *b1, float *b0, int iw, int ih, const char *name, int backward)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	int inc=1-(backward<<1);
	for(ptrdiff_t k=backward?res-1:0;k>=0&&k<res;k+=inc)
	{
		if(fabsf(b1[k]-b0[k])>1e-6)
		{
			ptrdiff_t kx=k%iw, ky=k/iw;
			printf("%s error XY (%5lld, %5lld) / %5d x %5d  b1 != b0  %f != %f  error=%f\n", name, kx, ky, iw, ih, b1[k], b0[k], b1[k]-b0[k]);
			return;
		}
	}
	printf("%s:\tSUCCESS\n", name);
}


#if 0
void print_matrix_pd(double *matrix, int bw, int bh)
{
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
			printf("\t%lf", matrix[bw*ky+kx]);
		printf("\n");
	}
	printf("\n");
}
#if 0
void impl_ref_pd(double *m, short dx, short dy)
{
#ifdef _DEBUG
	double pivot;
#endif
	double coeff;
	int mindim=dx<dy?dx:dy, it, ky, kx, npivots, kpivot;
	for(it=0, npivots=0;it<mindim;++it)//iteration
	{
		for(ky=npivots;ky<dy;++ky)//find pivot
		{
			if(m[dx*ky+it])
			{
#ifdef _DEBUG
				pivot=m[dx*ky+it];
#endif
				kpivot=ky;
				++npivots;
				break;
			}
		}
		if(ky<dy)
		{
			if(ky>it)
				for(kx=0;kx<dx;++kx)//swap rows
					coeff=m[dx*it+kx], m[dx*it+kx]=m[dx*ky+kx], m[dx*ky+kx]=coeff;
			for(++ky;ky<dy;++ky)//subtract pivot row
			{
				coeff=m[dx*ky+it]/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*ky+kx]-=coeff*m[dx*kpivot+kx];
			}
		}
	}
}
void impl_rref_pd(double *m, short dx, short dy)
{
#ifdef _DEBUG
	double pivot;
#endif
	double coeff;
	int mindim=dx<dy?dx:dy, it, ky, kx, npivots, kpivot;
	for(it=0, npivots=0;it<mindim;++it)//iteration
	{
		kpivot=-1;
		for(ky=npivots;ky<dy;++ky)//find pivot
		{
			if(m[dx*ky+it])
			{
#ifdef _DEBUG
				pivot=m[dx*ky+it];
#endif
				kpivot=ky;
				++npivots;
				break;
			}
		}
		if(kpivot==-1)
			continue;
		if(kpivot>npivots-1)
		{
			for(kx=0;kx<dx;++kx)//swap rows
				coeff=m[dx*kpivot+kx], m[dx*kpivot+kx]=m[dx*(npivots-1)+kx], m[dx*(npivots-1)+kx]=coeff;
			kpivot=npivots-1;
		}
		for(ky=0;ky<dy;++ky)
		{
			if(ky==kpivot)//normalize pivot row
			{
				coeff=0x100000000/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
				{
					int idx=dx*kpivot+kx;
					m[idx]*=coeff;
					//m[idx]>>=16;
				}
			}
			else//subtract pivot row from all other rows
			{
				coeff=(m[dx*ky+it])/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*ky+kx]-=coeff*m[dx*kpivot+kx];
			}
			//print_matrix_fixed(m, dx, dy, 16);
		}
	}
}
double impl_det_pd(double *m, int dx)//m is destroyed
{
	int k, dxplus1=dx+1;
	double result;

	//print_matrix_debug(m, dx, dx);//
	impl_ref(m, dx, dx);
	//print_matrix_debug(m, dx, dx);//

	result=m[0];//accumulate diagonal
	for(k=1;k<dx;++k)
		result=result*m[dxplus1*k];
	return result;
}
void impl_matinv_pd(double *m, short dx)//resize m to (dy * 2dx) temporarily,		dx==dy always
{
	int k, dy=dx, size=dx*dy;
			//print_matrix_fixed(m, dx<<1, dy, 16);
	for(k=size-dx;k>=0;k-=dx)//expand M into [M, 0]
	{
		memcpy(m+((size_t)k<<1), m+k, dx*sizeof(double));
		memset(m+((size_t)k<<1)+dx, 0, dx*sizeof(double));
	}
			//print_matrix_fixed(m, dx<<1, dy, 16);

	for(k=0;k<dx;++k)//add identity: [M, I]
		m[(dx<<1)*k+dx+k]=1;
			//print_matrix_fixed(m, dx<<1, dy, 16);

	impl_rref(m, dx<<1, dy);//[I, M^-1]
			//print_matrix_fixed(m, dx<<1, dy, 16);

	for(k=0;k<size;k+=dx)//pack M^-1
		memcpy(m+k, m+((size_t)k<<1)+dx, dx*sizeof(double));
			//print_matrix_fixed(m, dx<<1, dy, 16);
}
#endif
#endif
//float jxlparams_ps[33]=
//{
//		0.78f,    0.71f,    0.63f,   0.7f ,		-0.08f,   -0.01f,    0.59f,   0.12f,    -0.11f,   0.28f,    0.67f,
//		0.63f,    0.51f,    1.33f,   0.79f,		 0.28f,    0.02f,   -0.07f,   0.f  ,     0.01f,   0.39f,    0.15f,
//		0.7f ,    0.76f,    0.86f,   1.1f ,		-0.08f,   -0.06f,    0.38f,   0.04f,    -0.03f,   0.1f ,    0.91f,
//};
void apply_transforms_fwd(unsigned char *buf, int bw, int bh)
{
	ArrayHandle sizes=dwt2d_gensizes(bw, bh, 7, 7, 0);
	unsigned char *temp=(unsigned char*)malloc(MAXVAR(bw, bh));

	addbuf(buf, bw, bh, 3, 4, 128);//unsigned char -> signed char
	
	colortransform_ycocgt_fwd((char*)buf, bw, bh);
	//colortransform_ycocg_fwd((char*)buf, bw, bh);
	//colortransform_xgz_fwd((char*)buf, bw, bh);
	//colortransform_xyz_fwd((char*)buf, bw, bh);
	
	pred_jxl_apply((char*)buf, bw, bh, jxlparams_i16, 1);

	//pred_jxl((char*)buf, bw, bh, 0, 4, jxlpred_params   , 1);
	//pred_jxl((char*)buf, bw, bh, 1, 4, jxlpred_params+11, 1);
	//pred_jxl((char*)buf, bw, bh, 2, 4, jxlpred_params+22, 1);

	//pred_grad_fwd((char*)buf, bw, bh, 3, 4);
	//image_differentiate((char*)buf, bw, bh, 3, 4);
	//for(int kc=0;kc<3;++kc)
	//	//dwt2d_haar_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	//	//dwt2d_squeeze_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	//	dwt2d_cdf53_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	//	//dwt2d_cdf97_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);

	addbuf(buf, bw, bh, 3, 4, 128);
	
	//save_DWT_int8_all("kodim21-cdf53", buf, (DWTSize*)sizes->data, (int)sizes->count);
	//save_DWT_int8("kodim21-squeeze-stage", buf, (DWTSize*)sizes->data, (int)sizes->count, 4);//
	//save_DWT_int8("kodim21-cubic-stage", buf, (DWTSize*)sizes->data, (int)sizes->count, 4);//

	free(temp);
	array_free(&sizes);
}
void apply_transforms_inv(unsigned char *buf, int bw, int bh)
{
	ArrayHandle sizes=dwt2d_gensizes(bw, bh, 3, 3, 0);
	unsigned char *temp=(unsigned char*)malloc(MAXVAR(bw, bh));
	
	addbuf(buf, bw, bh, 3, 4, 128);
	
	pred_jxl_apply((char*)buf, bw, bh, jxlparams_i16, 0);

	//pred_jxl((char*)buf, bw, bh, 0, 4, jxlpred_params   , 0);
	//pred_jxl((char*)buf, bw, bh, 1, 4, jxlpred_params+11, 0);
	//pred_jxl((char*)buf, bw, bh, 2, 4, jxlpred_params+22, 0);

	//pred_grad_inv((char*)buf, bw, bh, 3, 4);
	//image_integrate((char*)buf, bw, bh, 3, 4);
	//for(int kc=0;kc<3;++kc)
	//	//dwt2d_haar_inv((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	//	dwt2d_squeeze_inv((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	//	//dwt2d_cdf53_inv((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	//	//dwt2d_cdf97_inv((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, (char*)temp);
	
	colortransform_ycocgt_inv((char*)buf, bw, bh);
	//colortransform_ycocg_inv((char*)buf, bw, bh);
	//colortransform_xgz_inv((char*)buf, bw, bh);
	//colortransform_xyz_inv((char*)buf, bw, bh);

	addbuf(buf, bw, bh, 3, 4, 128);//unsigned char -> signed char

	free(temp);
	array_free(&sizes);
}
void addbuf(unsigned char *buf, int iw, int ih, int nch, int bytestride, int ammount)
{
	for(int kp=0, len=iw*ih*bytestride;kp<len;kp+=bytestride)
	{
		for(int kc=0;kc<nch;++kc)
			buf[kp+kc]+=ammount;
	}
}


void image_differentiate(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					sub=left+top-topleft;
				buf[idx]-=sub;
			}
		}
	}
}
void image_integrate(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					sub=left+top-topleft;
				buf[idx]+=sub;
			}
		}
	}
}
void pred_grad_fwd(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					sub;

				if(topleft>MAXVAR(left, top))//planar prediction (unplane)
					sub=MINVAR(left, top);
				else if(topleft<MINVAR(left, top))
					sub=MAXVAR(left, top);
				else
					sub=left+top-topleft;

				buf[idx]-=sub;
			}
		}
	}
}
void pred_grad_inv(char *buf, int iw, int ih, int nch, int bytestride)
{
	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				char
					left=kx?buf[idx-bytestride]:0,
					top=ky?buf[idx-rowlen]:0,
					topleft=kx&&ky?buf[idx-rowlen-bytestride]:0,
					sub;

				if(topleft>MAXVAR(left, top))
					sub=MINVAR(left, top);
				else if(topleft<MINVAR(left, top))
					sub=MAXVAR(left, top);
				else
					sub=left+top-topleft;

				buf[idx]+=sub;
			}
		}
	}
}


double jxlpred_params[33]=
{
	 0.78,    0.71,    0.63,   0.7 ,
	-0.08,   -0.01,    0.59,   0.12,    -0.11,   0.28,    0.67,

	 0.63,    0.51,    1.33,   0.79,
	 0.28,    0.02,   -0.07,   0   ,     0.01,   0.39,    0.15,

	 0.7 ,    0.76,    0.86,   1.1 ,
	-0.08,   -0.06,    0.38,   0.04,    -0.03,   0.1 ,    0.91,

	// 0.93,    0.82,   0.71,   0.81,
	// 0.72,    0.51,   1.23,   0.89,
	// 0.84,    0.8 ,   0.89,   1.19,
	//-0.08,   -0.01,   0.58,   0.07,    -0.07,    0.32,    0.63,
	// 0.35,    0.01,  -0.03,   0.02,     0   ,    0.54,    0.17,
	//-0.09,   -0.06,   0.25,   0   ,     0.01,    0.27,    0.82,

	//0.860000, 0.860000, 0.830000, 0.830000,
	//0.790000, 0.600000, 0.970000, 0.970000,
	//1.200000, 0.800000, 0.920000, 1.170000,
	//-0.030000, -0.010000, 0.360000, 0.000000, -0.010000, 0.450000, 0.540000,
	//0.320000, -0.020000, 0.070000, -0.010000, 0.000000, 0.870000, 0.220000,
	//-0.010000, -0.030000, 0.150000, -0.010000, -0.000000, 0.600000, 0.560000,

	//1, 1, 1, 1,
	//1, 1, 1, 1,
	//1, 1, 1, 1,
};
short jxlparams_i16[33]=//signed fixed 7.8 bit
{
	0x0B37,  0x110B,  0x121B,  0x0BFC, -0x0001,  0x000E, -0x0188, -0x00E7, -0x00BB, -0x004A,  0x00BA,
	0x0DB8,  0x0E22,  0x181F,  0x0BF3, -0x005C, -0x005B,  0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,
	0x064C,  0x0F31,  0x1040,  0x0BF8, -0x0007, -0x000D, -0x0085, -0x0063, -0x00A2, -0x0017,  0x00F2,

	//0x0BD6, 0x10E7, 0x11F9, 0x0BEC,    0x0003,  0x0016, -0x01AC, -0x0127, -0x00C3, -0x0048,  0x00C3,
	//0x0DB8, 0x0E22, 0x181F, 0x0BF3,   -0x005C, -0x005B,  0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,
	//0x0684, 0x0F1F, 0x100D, 0x0C40,   -0x0007, -0x000D, -0x0087, -0x0063, -0x009F, -0x0018,  0x00F2,

	//0x0EEE, 0x103B, 0x0FAC, 0x0FF6, -0x0005,  0x0006, -0x0125, -0x007A, -0x004C,  0x0013,  0x0032,
	//0x0F92, 0x0FF3, 0x0FDF, 0x1006, -0x007A, -0x00C2, -0x000B, -0x006E,  0x007F, -0x004E, -0x00A0,
	//0x0E7E, 0x0FD7, 0x1024, 0x0FF4, -0x0011, -0x0048, -0x0054, -0x0053, -0x0038,  0x0018,  0x007D,

	//0x1000, 0x1000, 0x1000, 0x1000, 0, 0, 0, 0, 0, 0, 0,
	//0x1000, 0x1000, 0x1000, 0x1000, 0, 0, 0, 0, 0, 0, 0,
	//0x1000, 0x1000, 0x1000, 0x1000, 0, 0, 0, 0, 0, 0, 0,

	//0x0C7B, 0x0B5C, 0x0A14, 0x0B33,  -0x0148,-0x0029, 0, 0, 0, 0, 0,
	//0x0A14, 0x0829, 0x1548, 0x0CA4,   0x047B, 0x0052, 0, 0, 0, 0, 0,
	//0x0B33, 0x0C29, 0x0DC3, 0x119A,  -0x0148,-0x00F6, 0, 0, 0, 0, 0,

	//0x0C7B, 0x0B5C, 0x0A14, 0x0B33,-0x0148,-0x0029, 0x0971, 0x01EC,-0x01C3, 0x047B, 0x0AB8,
	//0x0A14, 0x0829, 0x1548, 0x0CA4, 0x047B, 0x0052,-0x011F, 0x0000, 0x0029, 0x063D, 0x0266,
	//0x0B33, 0x0C29, 0x0DC3, 0x119A,-0x0148,-0x00F6, 0x0614, 0x00A4,-0x007B, 0x019A, 0x0E8F,

	//0x63, 0x5A, 0x50, 0x59,    -0x0A, -0x01,  0x4B, 0x0F, -0x0E, -0x24, 0x55,
	//0x50, 0x41, 0xA9, 0x64,     0x24,  0x03, -0x09,    0,  0x01,  0x32, 0x13,
	//0x59, 0x61, 0x6D, 0x8C,    -0x0A, -0x08,  0x30, 0x05, -0x04,  0x0D, 0x74,
};
//float jxlparams_ps[33]=
//{
//	0.78f,    0.71f,    0.63f,   0.7f ,		-0.08f,   -0.01f,    0.59f,   0.12f,    -0.11f,   0.28f,    0.67f,
//	0.63f,    0.51f,    1.33f,   0.79f,		 0.28f,    0.02f,   -0.07f,   0.f  ,     0.01f,   0.39f,    0.15f,
//	0.7f ,    0.76f,    0.86f,   1.1f ,		-0.08f,   -0.06f,    0.38f,   0.04f,    -0.03f,   0.1f ,    0.91f,
//};
void pred_jxl_prealloc(const char *src, int iw, int ih, int kc, const short *params, int fwd, char *dst, int *temp_w10)
{
	int res=iw*ih, errorbuflen=iw*2, rowlen=iw<<2;
	int *error=temp_w10, *pred_errors[]=
	{
		temp_w10+errorbuflen,
		temp_w10+errorbuflen*2,
		temp_w10+errorbuflen*3,
		temp_w10+errorbuflen*4,
	};
	int idx=kc;
	for(int ky=0;ky<ih;++ky)
	{
		int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
		for(int kx=0;kx<iw;++kx, idx+=4)
		{
			int pred, curr;
			
			const char *src2=fwd?src:dst;
			char
				ctt      =         ky-2>=0?src2[idx-rowlen*2]:0,
				ctopleft =kx-1>=0&&ky-1>=0?src2[idx-rowlen-4]:0,
				ctop     =kx  <iw&&ky-1>=0?src2[idx-rowlen  ]:0,
				ctopright=kx+1<iw&&ky-1>=0?src2[idx-rowlen+4]:0,
				cleft    =kx-1>=0         ?src2[idx       -4]:0;

			//if(kx==(iw>>1)&&ky==(ih>>1))
			//	kx=iw>>1;

			//w0   w1   w2   w3
			//p3Ca p3Cb p3Cc p3Cd p3Ce
			//p1C  p2c
			
			int weights[4];//fixed 23.8 bit
			for(int k=0;k<4;++k)
			{
				int w=(ky-1>=0?pred_errors[k][prevrow+kx]:0)+(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0)+(ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0);
				weights[k]=(params[k]<<8)/(w+1);
			}

			int
				etop=ky-1>=0?error[prevrow+kx]:0,
				eleft=kx-1>=0?error[currrow+kx-1]:0,
				etopleft=ky-1>=0&&kx-1>=0?error[prevrow+kx-1]:0,
				etopright=ky-1>=0&&kx+1<iw?error[prevrow+kx+1]:0,
				sumtopleft=etop+eleft;
			int predictions[]=//fixed 23.8 bit
			{
				(cleft+ctopright-ctop)<<8,
				(ctop<<8)-((sumtopleft+etopright)*params[4]>>8),
				(cleft<<8)-((sumtopleft+etopleft)*params[5]>>8),

				//(ctop<<8)-(((etopleft*params[6])>>8)+ctop*params[7]+ctopright*params[8]+(ctt-ctop)*params[9]+(ctopleft-cleft)*params[10]),
				(ctop<<8)-(((etopleft*params[6]+etop*params[7]+etopright*params[8])>>8)+(ctt-ctop)*params[9]+(ctopleft-cleft)*params[10]),
			};

			int sum=weights[0]+weights[1]+weights[2]+weights[3];
			if(sum)
				pred=(predictions[0]*weights[0]+predictions[1]*weights[1]+predictions[2]*weights[2]+predictions[3]*weights[3]+(sum>>1)-1)/sum;
			else
				pred=predictions[0];

			int vmin=cleft, vmax=cleft;
			if(vmin>ctopright)
				vmin=ctopright;
			if(vmin>ctop)
				vmin=ctop;

			if(vmax<ctopright)
				vmax=ctopright;
			if(vmax<ctop)
				vmax=ctop;

			vmin<<=8;
			vmax<<=8;

			pred=CLAMP(vmin, pred, vmax);
			if(fwd)
			{
				curr=src[idx]<<8;
				dst[idx]=src[idx]-((pred+127)>>8);
			}
			else
			{
				dst[idx]=src[idx]+((pred+127)>>8);
				curr=dst[idx]<<8;
			}

			error[currrow+kx]=curr-pred;
			for(int k=0;k<4;++k)
			{
				int e=abs(curr-predictions[k]);
				pred_errors[k][currrow+kx]=e;
				if(kx+1<iw)
					pred_errors[k][prevrow+kx+1]+=e;
			}
		}
	}
}
double pred_jxl_calccsize(const char *src, int iw, int ih, int kc, const short *params, int *temp, char *dst, int *hist, int it, int loud)
{
	int res=iw*ih;
	pred_jxl_prealloc(src, iw, ih, kc, params, 1, dst, temp);
	//addbuf(dst+kc, iw, ih, 1, 4, 128);//just slows down the optimization
	calc_histogram(dst+kc, (ptrdiff_t)res<<2, 4, hist);

	double entropy=0;
	int freq;
	double p;
	for(int k=0;k<256;++k)
	{
		freq=hist[k];
		if(freq)
		{
			p=(double)freq/res;
			p*=0x10000-255;
			++p;
			p/=0x10000;
			entropy-=p*log2(p);
		}
	}
	double invCR=entropy/8, csize=res*invCR;
	//if(loud)
	//	printf("%4d %14lf\r", it, csize);
	return csize;
}
double pred_jxl_optimize(const char *src, int iw, int ih, int kc, short *params, int step, int pidx, char *dst, int loud)
{
	int *temp=(int*)malloc((size_t)iw*10*sizeof(int));
	int *hist=(int*)malloc(256*sizeof(int));
	if(!temp||!hist)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	
	int res=iw*ih;
	double csize0, csize, csize00;
	short p0=params[pidx];
	csize=pred_jxl_calccsize(src, iw, ih, kc, params, temp, dst, hist, pidx, loud);
	csize00=csize;

	int subit;
	for(subit=0;subit<20;++subit)
	{
		csize0=csize;
		params[pidx]+=step;
		csize=pred_jxl_calccsize(src, iw, ih, kc, params, temp, dst, hist, pidx, loud);
		if(csize>=csize0)//cancel last change and break
		{
			params[pidx]-=step;
			csize=csize0;
			break;
		}
	}
		
	for(subit=0;subit<20;++subit)
	{
		csize0=csize;
		params[pidx]-=step;
		csize=pred_jxl_calccsize(src, iw, ih, kc, params, temp, dst, hist, pidx, loud);
		if(csize>=csize0)
		{
			params[pidx]+=step;
			csize=csize0;
			break;
		}
	}
	free(hist);
	free(temp);

	if(csize>=csize00)//prevent CR from worsening
	{
		params[pidx]=p0;
		csize=csize00;
	}
	//if(loud)
	//	printf("%4d %14lf\r", pidx, csize);
	return csize;
}
void pred_jxl_apply(char *buf, int iw, int ih, short *allparams, int fwd)
{
	int res=iw*ih;
	int *temp=(int*)malloc((size_t)iw*10*sizeof(int));
	char *buf2=(char*)malloc((size_t)res<<2);
	if(!temp||!buf2)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	
	pred_jxl_prealloc(buf, iw, ih, 0, allparams   , fwd, buf2, temp);
	pred_jxl_prealloc(buf, iw, ih, 1, allparams+11, fwd, buf2, temp);
	pred_jxl_prealloc(buf, iw, ih, 2, allparams+22, fwd, buf2, temp);

	for(int k=0;k<res;++k)
	{
		buf[k<<2  ]=buf2[k<<2  ];
		buf[k<<2|1]=buf2[k<<2|1];
		buf[k<<2|2]=buf2[k<<2|2];
	}

	free(temp);
	free(buf2);
}
#if 0
void pred_jxl(char *buf, int iw, int ih, int kc, int bytestride, double *params, int fwd)
{
	int res=iw*ih;
	char *b2=(char*)malloc((size_t)res*bytestride);
	int errorbuflen=iw*2;
	char *pred_errors[]=
	{
		(char*)malloc((size_t)errorbuflen),
		(char*)malloc((size_t)errorbuflen),
		(char*)malloc((size_t)errorbuflen),
		(char*)malloc((size_t)errorbuflen),
	};
	char *error=(char*)malloc((size_t)errorbuflen);
	if(!b2||!pred_errors[0]||!pred_errors[1]||!pred_errors[2]||!pred_errors[3]||!error)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(pred_errors[0], 0, errorbuflen);
	memset(pred_errors[1], 0, errorbuflen);
	memset(pred_errors[2], 0, errorbuflen);
	memset(pred_errors[3], 0, errorbuflen);
	memset(error, 0, errorbuflen);
	memset(b2, 0, (size_t)res*bytestride);
	int rowlen=iw*bytestride;

	int idx=kc;
	for(int ky=0;ky<ih;++ky)
	{
		int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
		for(int kx=0;kx<iw;++kx, idx+=bytestride)
		{
			int pred, curr;
			
			char *src=fwd?buf:b2;
			char
				ctt      =         ky-2>=0?src[idx-rowlen*2]:0,

				ctopleft =kx-1>=0&&ky-1>=0?src[idx-rowlen-bytestride  ]:0,
				ctop     =kx  <iw&&ky-1>=0?src[idx-rowlen             ]:0,
				ctopright=kx+1<iw&&ky-1>=0?src[idx-rowlen+bytestride  ]:0,
				ctrr     =kx+2<iw&&ky-1>=0?src[idx-rowlen+bytestride*2]:0,
				
				cl3      =kx-3>=0?src[idx-bytestride*3]:0,
				cll      =kx-2>=0?src[idx-bytestride*2]:0,
				cleft    =kx-1>=0?src[idx-bytestride  ]:0;


			//w0   w1   w2   w3
			//p3Ca p3Cb p3Cc p3Cd p3Ce
			//p1C  p2c

			double weights[4];
			for(int k=0;k<4;++k)
			{
				int w=(ky-1>=0?pred_errors[k][prevrow+kx]:0)+(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0)+(ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0);
				weights[k]=params[k]/(w+1);
			}

			char
				etop=ky-1>=0?error[prevrow+kx]:0,
				eleft=kx-1>=0?error[currrow+kx-1]:0,
				etopleft=ky-1>=0&&kx-1>=0?error[prevrow+kx-1]:0,
				etopright=ky-1>=0&&kx+1<iw?error[prevrow+kx+1]:0;
			double predictions[]=
			{
				cleft+ctopright-ctop,
				ctop-(int)((etop+eleft+etopright)*params[4]),
				cleft-(int)((etop+eleft+etopleft)*params[5]),
				ctop-(int)(etopleft*params[6]+ctop*params[7]+ctopright*params[8]+(ctt-ctop)*params[9]+(ctopleft-cleft)*params[10]),
			};

			double sum=weights[0]+weights[1]+weights[2]+weights[3];
			if(sum)
				pred=(int)round((predictions[0]*weights[0]+predictions[1]*weights[1]+predictions[2]*weights[2]+predictions[3]*weights[3])/sum);
			else
				pred=(int)round(predictions[0]);

			int vmin=cleft, vmax=cleft;
			if(vmin>ctopright)
				vmin=ctopright;
			if(vmin>ctop)
				vmin=ctop;

			if(vmax<ctopright)
				vmax=ctopright;
			if(vmax<ctop)
				vmax=ctop;

			pred=CLAMP(vmin, pred, vmax);
			//pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
			if(fwd)
			{
				curr=buf[idx];
				b2[idx]=buf[idx]-pred;
			}
			else
			{
				b2[idx]=buf[idx]+pred;
				curr=b2[idx];
			}

			error[currrow+kx]=pred-curr;
			for(int k=0;k<4;++k)
			{
				int e=(int)round(fabs(curr-predictions[k]));
				pred_errors[k][currrow+kx]=e;
				if(kx+1<iw)
					pred_errors[k][prevrow+kx+1]+=e;
			}
		}
	}
	for(int k=0;k<res;++k)
	{
		int idx=k*bytestride+kc;
		buf[idx]=b2[idx];
	}
	free(b2);
	free(pred_errors[0]);
	free(pred_errors[1]);
	free(pred_errors[2]);
	free(pred_errors[3]);
	free(error);
}
#endif


//YCoCg-R (lifting-based YCoCg)		8 bit <-> ps,  3 channels, pixel stride 4 bytes,  Y in [0, 1], Co/Cg in [-1, 1],  used by AVC/HEVC/VVC
#ifndef CLAMP
#define CLAMP(LO, X, HI)	(X<LO?LO:(X>HI?HI:X))
#endif
void YCoCg_8bit_ps_fwd(const unsigned char *src, ptrdiff_t res, float *bufY, float *bufCo, float *bufCg)
{
	const float gain=1.f/255;
	for(ptrdiff_t k=0;k<res;++k)
	{
		ptrdiff_t idx=k<<2;
		unsigned R=src[idx], G=src[idx|1], B=src[idx|2];
		int Co=R-B,
			temp=B+(Co>>1),
			Cg=G-temp,
			Y=temp+(Cg>>1);
		bufY[k]=(float)Y*gain;
		bufCo[k]=(float)Co*gain;
		bufCg[k]=(float)Cg*gain;
	}
}
void YCoCg_8bit_ps_inv(const float *bufY, const float *bufCo, const float *bufCg, ptrdiff_t res, unsigned char *dst)
{
	const float gain=255.f;
	for(ptrdiff_t k=0;k<res;++k)
	{
		ptrdiff_t idx=k<<2;
		int Y=(int)(gain*bufY[k]), Co=(int)(gain*bufCo[k]), Cg=(int)(gain*bufCg[k]);
		int temp=Y-(Cg>>1), G=Cg+temp, B=temp-(Co>>1), R=B+Co;
		R=CLAMP(0, R, 255);
		G=CLAMP(0, G, 255);
		B=CLAMP(0, B, 255);
		dst[idx]=R, dst[idx|1]=G, dst[idx|3]=B;
	}
}

void colortransform_xgz_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;
		b-=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_xgz_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		b+=g;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_xyz_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;//what is this transform?
		b-=g;
		g+=(r+b)>>1;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;

		//buf[k  ]=r-g+128;//XGZ
		//buf[k|1]=g;
		//buf[k|2]=b-g+128;

		//buf[k  ]=r;//RYZ
		//buf[k|1]=g-r+128;
		//buf[k|2]=b-r+128;

		//buf[k  ]=r-b+128;//XYBdash
		//buf[k|1]=g-b+128;
		//buf[k|2]=b;
	}
}
void colortransform_xyz_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		g-=(r+b)>>1;//what is this transform???
		b+=g;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocg_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=b;
		b+=r>>1;
		g-=b;
		b+=g>>1;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocg_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];
		
		b-=g>>1;
		g+=b;
		b-=r>>1;
		r+=b;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocgt_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;
		g+=r>>1;
		b-=g;
		g+=b>>1;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocgt_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];
		
		g-=b>>1;
		b+=g;
		g-=r>>1;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}

float lrct[]=
{
	-0.235968455672264100f, -0.005472748540341854f,  0.546232640743255600f,
	 0.042025841772556305f, -0.281145840883255000f, -0.041713438928127290f,
	-0.270333796739578250f,  0.382541775703430200f,  0.122593350708484650f,
	-0.654174447059631300f, -0.111721098423004150f,  0.049731694161891940f,
	-0.180651888251304630f, -0.665180742740631100f, -0.449275195598602300f,
	-0.634068489074707000f, -0.438444226980209350f, -0.079102315008640290f,
	-0.492376208305358900f, -0.011669249273836613f,  0.612039744853973400f,
	-0.140227079391479500f,  0.255208671092987060f, -0.405863404273986800f,
	-0.133800372481346130f, -0.049969632178545000f, -0.006972887087613344f,
	-0.498686790466308600f, -0.417643666267395000f,  0.129942402243614200f,
	-0.584400236606597900f, -0.273861020803451540f,  0.104215130209922790f,
	 0.297003507614135740f, -0.688771069049835200f,  0.330031186342239400f,
};
char lift(char v1, char v2, float *coeff)
{
	return (char)(127*(v1*coeff[0]/127+v2*coeff[1]/127+coeff[2]));
}
void colortransform_learned_fwd(char *buf, int iw, int ih)
{
	for(int k=0, res=iw*ih;k<res;++k)
	{
		char r=buf[k<<2], g=buf[k<<2|1], b=buf[k<<2|2];
		r-=lift(g, b, lrct+3* 0);
		g-=lift(r, b, lrct+3* 1);
		b-=lift(r, g, lrct+3* 2);
		r+=lift(g, b, lrct+3* 3);
		g+=lift(r, b, lrct+3* 4);
		b+=lift(r, g, lrct+3* 5);
		r-=lift(g, b, lrct+3* 6);
		g-=lift(r, b, lrct+3* 7);
		b-=lift(r, g, lrct+3* 8);
		r+=lift(g, b, lrct+3* 9);
		g+=lift(r, b, lrct+3*10);
		b+=lift(r, g, lrct+3*11);
		buf[k<<2]=r, buf[k<<2|1]=g, buf[k<<2|2]=b;
	}
}
void colortransform_learned_inv(char *buf, int iw, int ih)
{
	for(int k=0, res=iw*ih;k<res;++k)
	{
		char r=buf[k<<2], g=buf[k<<2|1], b=buf[k<<2|2];
		b-=lift(r, g, lrct+3*11);
		g-=lift(r, b, lrct+3*10);
		r-=lift(g, b, lrct+3* 9);
		b+=lift(r, g, lrct+3* 8);
		g+=lift(r, b, lrct+3* 7);
		r+=lift(g, b, lrct+3* 6);
		b-=lift(r, g, lrct+3* 5);
		g-=lift(r, b, lrct+3* 4);
		r-=lift(g, b, lrct+3* 3);
		b+=lift(r, g, lrct+3* 2);
		g+=lift(r, b, lrct+3* 1);
		r+=lift(g, b, lrct+3* 0);
		buf[k<<2]=r, buf[k<<2|1]=g, buf[k<<2|2]=b;
	}
}
#if 0
#if 1
float biases[]=
{
	-0.007593977730721235f,
	-0.20463228225708008f,
	0.048144787549972534f,
	-0.052841395139694214f,
	-0.229848250746727f,
	-0.10854046046733856f,
	-0.27017149329185486f,
	-0.28025728464126587f,
	0.14076033234596252f,
	0.11415114998817444f,
	-0.24117304384708405f,
	-0.22909477353096008f,
};
float lrt_c01[]=
{
	 0.21996481716632843f,
	 0.3447851538658142f ,
	 0.4309721887111664f ,
	-0.3914514482021332f ,
	 0.1072743684053421f ,
	 0.09292115271091461f,
	 0.09527700394392014f,
	-0.7720031142234802f ,
	-0.03891817852854729f,
	 0.18422403931617737f,
	 0.6045548319816589f ,
};
float lrt_c02[]=
{
	0.16338559985160828f,
	-0.21537575125694275f,
	-0.01003316044807434f,
	0.0806010365486145f,
	0.1516939401626587f,
	0.22873657941818237f,
	-0.13795271515846252f,
	0.13615399599075317f,
	-0.07778285443782806f,
	0.291059672832489f,
	0.22079020738601685f,
};
float lrt_c03[]=
{
	0.010352730751037598f,
	0.1270458996295929f,
	-0.18215587735176086f,
	-0.164699524641037f,
	0.23722952604293823f,
	0.27438366413116455f,
	-0.05265103280544281f,
	0.18550443649291992f,
	-0.21034197509288788f,
	-0.13812671601772308f,
	-0.015580296516418457f,
};
float lrt_c04[]=
{
	-0.06830897927284241f,
	-0.0477299690246582f,
	-0.2479163110256195f,
	-0.04657846689224243f,
	-0.0802721232175827f,
	-0.004736065864562988f,
	-0.07743765413761139f,
	-0.26771998405456543f,
	0.013397663831710815f,
	-0.23979492485523224f,
	0.12400856614112854f,
};
float lrt_c05[]=
{
	-0.04351922869682312f,
	-0.11087171733379364f,
	-0.11081892251968384f,
	0.301426351070404050f,
	0.090425968170166020f,
	0.171547293663024900f,
	0.166505038738250730f,
	0.086036622524261470f,
	-0.25268214941024780f,
	0.168854206800460820f,
	0.221817135810852050f,
};
float lrt_c06[]=
{
	0.216941952705383300f,
	0.145165383815765380f,
	-0.30066022276878357f,
	0.062808305025100710f,
	-0.11871317028999329f,
	-0.03821510076522827f,
	0.057682424783706665f,
	-0.15345700085163116f,
	-0.10404434800148010f,
	0.159864872694015500f,
	0.030692487955093384f,
};
float lrt_c07[]=
{
	-0.11159592866897583f,
	0.267521142959594700f,
	-0.04237860441207886f,
	-0.12383817136287689f,
	0.221899330615997310f,
	0.031740218400955200f,
	-0.11875738203525543f,
	0.109959483146667480f,
	-0.17683275043964386f,
	-0.12659740447998047f,
	-0.21873277425765990f,
};
float lrt_c08[]=
{
	-0.18569713830947876f,
	0.224896728992462160f,
	-0.28754058480262756f,
	0.111558407545089720f,
	-0.18325410783290863f,
	-0.27818745374679565f,
	0.241845905780792240f,
	0.179279416799545300f,
	-0.05782620608806610f,
	-0.10470718145370483f,
	0.128715872764587400f,
};
float lrt_c09[]=
{
	0.0736332833766937300f,
	-0.278686881065368650f,
	-0.062273472547531130f,
	0.2923494577407837000f,
	0.1254649162292480500f,
	0.0480349659919738800f,
	0.2714382410049438500f,
	-0.016680717468261720f,
	0.0569746196269989000f,
	-0.021338820457458496f,
	0.0394396781921386700f,
};
float lrt_c10[]=
{
	0.068611204624176030f,
	-0.12615610659122467f,
	-0.07715198397636414f,
	0.252162516117095950f,
	-0.02322089672088623f,
	0.062368571758270264f,
	0.138979107141494750f,
	0.039735138416290280f,
	0.171300053596496580f,
	-0.13656683266162872f,
	-0.07565732300281525f,
};
float lrt_c11[]=
{
	0.1174532473087310800f,
	0.0668951570987701400f,
	-0.071375116705894470f,
	-0.240586921572685240f,
	0.1582727432250976600f,
	0.1294519901275634800f,
	0.2438962459564209000f,
	-0.014519184827804565f,
	0.1907848119735717800f,
	-0.215079009532928470f,
	0.0485148429870605500f,
};
float lrt_c12[]=
{
	0.086784452199935910f,
	0.276734590530395500f,
	-0.24387010931968690f,
	0.102736681699752810f,
	-0.28351089358329773f,
	0.072889894247055050f,
	-0.18153974413871765f,
	0.022794693708419800f,
	-0.27390283346176150f,
	-0.27018001675605774f,
	0.102353841066360470f,
};
#endif
void learnedtransform_fwd(char *buf, int iw, int ih)
{
	for(int ky=0;ky<ih-1;ky+=2)
	{
		for(int kx=0;kx<iw-1;kx+=2)
		{
			int idx=iw*ky+kx;
			char
				x01=buf[idx<<2  ], x02=buf[(idx+1)<<2  ], x03=buf[(idx+iw)<<2  ], x04=buf[(idx+iw+1)<<2  ],
				x05=buf[idx<<2|1], x06=buf[(idx+1)<<2|1], x07=buf[(idx+iw)<<2|1], x08=buf[(idx+iw+1)<<2|1],
				x09=buf[idx<<2|2], x10=buf[(idx+1)<<2|2], x11=buf[(idx+iw)<<2|2], x12=buf[(idx+iw+1)<<2|2];
			x01-=(char)(lrt_c01[0]*x02/255+lrt_c01[1]*x03/255+lrt_c01[2]*x04/255+lrt_c01[3]*x05/255+lrt_c01[4]*x06/255+lrt_c01[5]*x07/255+lrt_c01[6]*x08/255+lrt_c01[7]*x09/255+lrt_c01[8]*x10/255+lrt_c01[9]*x11/255+lrt_c01[10]*x12/255);
			x02-=(char)(lrt_c02[0]*x02/255+lrt_c02[1]*x03/255+lrt_c02[2]*x04/255+lrt_c02[3]*x05/255+lrt_c02[4]*x06/255+lrt_c02[5]*x07/255+lrt_c02[6]*x08/255+lrt_c02[7]*x09/255+lrt_c02[8]*x10/255+lrt_c02[9]*x11/255+lrt_c02[10]*x12/255);
			x03-=(char)(lrt_c03[0]*x02/255+lrt_c03[1]*x03/255+lrt_c03[2]*x04/255+lrt_c03[3]*x05/255+lrt_c03[4]*x06/255+lrt_c03[5]*x07/255+lrt_c03[6]*x08/255+lrt_c03[7]*x09/255+lrt_c03[8]*x10/255+lrt_c03[9]*x11/255+lrt_c03[10]*x12/255);
			x04-=(char)(lrt_c04[0]*x02/255+lrt_c04[1]*x03/255+lrt_c04[2]*x04/255+lrt_c04[3]*x05/255+lrt_c04[4]*x06/255+lrt_c04[5]*x07/255+lrt_c04[6]*x08/255+lrt_c04[7]*x09/255+lrt_c04[8]*x10/255+lrt_c04[9]*x11/255+lrt_c04[10]*x12/255);
			x05-=(char)(lrt_c05[0]*x02/255+lrt_c05[1]*x03/255+lrt_c05[2]*x04/255+lrt_c05[3]*x05/255+lrt_c05[4]*x06/255+lrt_c05[5]*x07/255+lrt_c05[6]*x08/255+lrt_c05[7]*x09/255+lrt_c05[8]*x10/255+lrt_c05[9]*x11/255+lrt_c05[10]*x12/255);
			x06-=(char)(lrt_c06[0]*x02/255+lrt_c06[1]*x03/255+lrt_c06[2]*x04/255+lrt_c06[3]*x05/255+lrt_c06[4]*x06/255+lrt_c06[5]*x07/255+lrt_c06[6]*x08/255+lrt_c06[7]*x09/255+lrt_c06[8]*x10/255+lrt_c06[9]*x11/255+lrt_c06[10]*x12/255);
			x07-=(char)(lrt_c07[0]*x02/255+lrt_c07[1]*x03/255+lrt_c07[2]*x04/255+lrt_c07[3]*x05/255+lrt_c07[4]*x06/255+lrt_c07[5]*x07/255+lrt_c07[6]*x08/255+lrt_c07[7]*x09/255+lrt_c07[8]*x10/255+lrt_c07[9]*x11/255+lrt_c07[10]*x12/255);
			x08-=(char)(lrt_c08[0]*x02/255+lrt_c08[1]*x03/255+lrt_c08[2]*x04/255+lrt_c08[3]*x05/255+lrt_c08[4]*x06/255+lrt_c08[5]*x07/255+lrt_c08[6]*x08/255+lrt_c08[7]*x09/255+lrt_c08[8]*x10/255+lrt_c08[9]*x11/255+lrt_c08[10]*x12/255);
			x09-=(char)(lrt_c09[0]*x02/255+lrt_c09[1]*x03/255+lrt_c09[2]*x04/255+lrt_c09[3]*x05/255+lrt_c09[4]*x06/255+lrt_c09[5]*x07/255+lrt_c09[6]*x08/255+lrt_c09[7]*x09/255+lrt_c09[8]*x10/255+lrt_c09[9]*x11/255+lrt_c09[10]*x12/255);
			x10-=(char)(lrt_c10[0]*x02/255+lrt_c10[1]*x03/255+lrt_c10[2]*x04/255+lrt_c10[3]*x05/255+lrt_c10[4]*x06/255+lrt_c10[5]*x07/255+lrt_c10[6]*x08/255+lrt_c10[7]*x09/255+lrt_c10[8]*x10/255+lrt_c10[9]*x11/255+lrt_c10[10]*x12/255);
			x11-=(char)(lrt_c11[0]*x02/255+lrt_c11[1]*x03/255+lrt_c11[2]*x04/255+lrt_c11[3]*x05/255+lrt_c11[4]*x06/255+lrt_c11[5]*x07/255+lrt_c11[6]*x08/255+lrt_c11[7]*x09/255+lrt_c11[8]*x10/255+lrt_c11[9]*x11/255+lrt_c11[10]*x12/255);
			x12-=(char)(lrt_c12[0]*x02/255+lrt_c12[1]*x03/255+lrt_c12[2]*x04/255+lrt_c12[3]*x05/255+lrt_c12[4]*x06/255+lrt_c12[5]*x07/255+lrt_c12[6]*x08/255+lrt_c12[7]*x09/255+lrt_c12[8]*x10/255+lrt_c12[9]*x11/255+lrt_c12[10]*x12/255);
			
			buf[idx<<2  ]=x01, buf[(idx+1)<<2  ]=x02, buf[(idx+iw)<<2  ]=x03, buf[(idx+iw+1)<<2  ]=x04,
			buf[idx<<2|1]=x05, buf[(idx+1)<<2|1]=x06, buf[(idx+iw)<<2|1]=x07, buf[(idx+iw+1)<<2|1]=x08,
			buf[idx<<2|2]=x09, buf[(idx+1)<<2|2]=x10, buf[(idx+iw)<<2|2]=x11, buf[(idx+iw+1)<<2|2]=x12;
		}
	}
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 1);
	for(int kc=0;kc<3;++kc)
		dwt2d_lazy_fwd(buf, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
	array_free(&sizes);
	free(temp);
}
void learnedtransform_inv(char *buf, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 1);
	for(int kc=0;kc<3;++kc)
		dwt2d_lazy_inv(buf, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
	array_free(&sizes);
	free(temp);

	for(int ky=0;ky<ih-1;ky+=2)
	{
		for(int kx=0;kx<iw-1;kx+=2)
		{
			int idx=iw*ky+kx;
			char
				x01=buf[idx<<2  ], x02=buf[(idx+1)<<2  ], x03=buf[(idx+iw)<<2  ], x04=buf[(idx+iw+1)<<2  ],
				x05=buf[idx<<2|1], x06=buf[(idx+1)<<2|1], x07=buf[(idx+iw)<<2|1], x08=buf[(idx+iw+1)<<2|1],
				x09=buf[idx<<2|2], x10=buf[(idx+1)<<2|2], x11=buf[(idx+iw)<<2|2], x12=buf[(idx+iw+1)<<2|2];
			x12+=(char)(lrt_c12[0]*x02/255+lrt_c12[1]*x03/255+lrt_c12[2]*x04/255+lrt_c12[3]*x05/255+lrt_c12[4]*x06/255+lrt_c12[5]*x07/255+lrt_c12[6]*x08/255+lrt_c12[7]*x09/255+lrt_c12[8]*x10/255+lrt_c12[9]*x11/255+lrt_c12[10]*x12/255);
			x11+=(char)(lrt_c11[0]*x02/255+lrt_c11[1]*x03/255+lrt_c11[2]*x04/255+lrt_c11[3]*x05/255+lrt_c11[4]*x06/255+lrt_c11[5]*x07/255+lrt_c11[6]*x08/255+lrt_c11[7]*x09/255+lrt_c11[8]*x10/255+lrt_c11[9]*x11/255+lrt_c11[10]*x12/255);
			x10+=(char)(lrt_c10[0]*x02/255+lrt_c10[1]*x03/255+lrt_c10[2]*x04/255+lrt_c10[3]*x05/255+lrt_c10[4]*x06/255+lrt_c10[5]*x07/255+lrt_c10[6]*x08/255+lrt_c10[7]*x09/255+lrt_c10[8]*x10/255+lrt_c10[9]*x11/255+lrt_c10[10]*x12/255);
			x09+=(char)(lrt_c09[0]*x02/255+lrt_c09[1]*x03/255+lrt_c09[2]*x04/255+lrt_c09[3]*x05/255+lrt_c09[4]*x06/255+lrt_c09[5]*x07/255+lrt_c09[6]*x08/255+lrt_c09[7]*x09/255+lrt_c09[8]*x10/255+lrt_c09[9]*x11/255+lrt_c09[10]*x12/255);
			x08+=(char)(lrt_c08[0]*x02/255+lrt_c08[1]*x03/255+lrt_c08[2]*x04/255+lrt_c08[3]*x05/255+lrt_c08[4]*x06/255+lrt_c08[5]*x07/255+lrt_c08[6]*x08/255+lrt_c08[7]*x09/255+lrt_c08[8]*x10/255+lrt_c08[9]*x11/255+lrt_c08[10]*x12/255);
			x07+=(char)(lrt_c07[0]*x02/255+lrt_c07[1]*x03/255+lrt_c07[2]*x04/255+lrt_c07[3]*x05/255+lrt_c07[4]*x06/255+lrt_c07[5]*x07/255+lrt_c07[6]*x08/255+lrt_c07[7]*x09/255+lrt_c07[8]*x10/255+lrt_c07[9]*x11/255+lrt_c07[10]*x12/255);
			x06+=(char)(lrt_c06[0]*x02/255+lrt_c06[1]*x03/255+lrt_c06[2]*x04/255+lrt_c06[3]*x05/255+lrt_c06[4]*x06/255+lrt_c06[5]*x07/255+lrt_c06[6]*x08/255+lrt_c06[7]*x09/255+lrt_c06[8]*x10/255+lrt_c06[9]*x11/255+lrt_c06[10]*x12/255);
			x05+=(char)(lrt_c05[0]*x02/255+lrt_c05[1]*x03/255+lrt_c05[2]*x04/255+lrt_c05[3]*x05/255+lrt_c05[4]*x06/255+lrt_c05[5]*x07/255+lrt_c05[6]*x08/255+lrt_c05[7]*x09/255+lrt_c05[8]*x10/255+lrt_c05[9]*x11/255+lrt_c05[10]*x12/255);
			x04+=(char)(lrt_c04[0]*x02/255+lrt_c04[1]*x03/255+lrt_c04[2]*x04/255+lrt_c04[3]*x05/255+lrt_c04[4]*x06/255+lrt_c04[5]*x07/255+lrt_c04[6]*x08/255+lrt_c04[7]*x09/255+lrt_c04[8]*x10/255+lrt_c04[9]*x11/255+lrt_c04[10]*x12/255);
			x03+=(char)(lrt_c03[0]*x02/255+lrt_c03[1]*x03/255+lrt_c03[2]*x04/255+lrt_c03[3]*x05/255+lrt_c03[4]*x06/255+lrt_c03[5]*x07/255+lrt_c03[6]*x08/255+lrt_c03[7]*x09/255+lrt_c03[8]*x10/255+lrt_c03[9]*x11/255+lrt_c03[10]*x12/255);
			x02+=(char)(lrt_c02[0]*x02/255+lrt_c02[1]*x03/255+lrt_c02[2]*x04/255+lrt_c02[3]*x05/255+lrt_c02[4]*x06/255+lrt_c02[5]*x07/255+lrt_c02[6]*x08/255+lrt_c02[7]*x09/255+lrt_c02[8]*x10/255+lrt_c02[9]*x11/255+lrt_c02[10]*x12/255);
			x01+=(char)(lrt_c01[0]*x02/255+lrt_c01[1]*x03/255+lrt_c01[2]*x04/255+lrt_c01[3]*x05/255+lrt_c01[4]*x06/255+lrt_c01[5]*x07/255+lrt_c01[6]*x08/255+lrt_c01[7]*x09/255+lrt_c01[8]*x10/255+lrt_c01[9]*x11/255+lrt_c01[10]*x12/255);
			
			buf[idx<<2  ]=x01, buf[(idx+1)<<2  ]=x02, buf[(idx+iw)<<2  ]=x03, buf[(idx+iw+1)<<2  ]=x04,
			buf[idx<<2|1]=x05, buf[(idx+1)<<2|1]=x06, buf[(idx+iw)<<2|1]=x07, buf[(idx+iw+1)<<2|1]=x08,
			buf[idx<<2|2]=x09, buf[(idx+1)<<2|2]=x10, buf[(idx+iw)<<2|2]=x11, buf[(idx+iw+1)<<2|2]=x12;
		}
	}
}
#endif




//lifting-based 8bit lazy DWT
void dwt1d_lazy_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_lazy_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_lazy_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_lazy_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_lazy_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_lazy_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_lazy_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_lazy_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit Haar DWT
void dwt1d_haar_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(int k=0;k<nodd;++k)//const predictor (difference)
		even[k]-=odd[k];
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(int k=0;k<nodd;++k)//update (O+(E-O)/2 = average)
		odd[k]+=even[k]>>1;


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_haar_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	
	for(int k=0;k<nodd;++k)//unupdate
		odd[k]-=even[k]>>1;
	
	for(int k=0;k<nodd;++k)//unpredict
		even[k]+=odd[k];
	if(extraeven)
		even[nodd]+=odd[nodd-1];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_haar_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_haar_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_haar_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_haar_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_haar_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_haar_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit squeeze DWT
char smoothtendency(char B, char a, char n)
{
	char diff=0;
	if(B>=a&&a>=n)
	{
		diff=(4*B-3*n-a+6)/12;
		//      2C = a<<1 + diff - diff&1 <= 2B  so diff - diff&1 <= 2B - 2a
		//      2D = a<<1 - diff - diff&1 >= 2n  so diff + diff&1 <= 2a - 2n
		if(diff-(diff&1)>2*(B-a))
			diff=2*(B-a)+1;
		if(diff+(diff&1)>2*(a-n))
			diff=2*(a-n);
	}
	else if(B<=a&&a<=n)
	{
		diff=(4*B-3*n-a-6)/12;
		//      2C = a<<1 + diff + diff&1 >= 2B  so diff + diff&1 >= 2B - 2a
		//      2D = a<<1 - diff + diff&1 <= 2n  so diff - diff&1 >= 2a - 2n
		if(diff+(diff&1)<2*(B-a))
			diff=2*(B-a)-1;
		if(diff-(diff&1)<2*(a-n))
			diff=2*(a-n);
	}
	return diff;
}
void dwt1d_squeeze_fwd(char *buffer, int count, int stride, char *b2)//nOdd <= nEven			nOdd = nEven - (count&1)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int ks=0, kd=0;kd<nodd+extraeven;ks+=stride<<1, ++kd)
	{
		char
			o1=kd>0?buffer[ks-stride]:0,
			e =buffer[ks],
			o =(kd<<1)+1<count?buffer[ks+stride]:0,
			e2=(kd<<1)+2<count?buffer[ks+(stride<<1)]:0,//n-1-(idx-(n-1)) = ((n-1)<<1)-idx
			o2=(kd<<1)+3<count?buffer[ks+ stride*3  ]:0;
		//char
		//	o1=kd>0?buffer[ks-stride]:buffer[ks+stride],
		//	e =buffer[ks],
		//	o =(kd<<1)+1<count?buffer[ks+stride]:buffer[((count-1)*stride<<1)-(ks+stride)],
		//	e2=(kd<<1)+2<count?buffer[ks+(stride<<1)]:buffer[((count-1)*stride<<1)-(ks+(stride<<1))],//n-1-(idx-(n-1)) = ((n-1)<<1)-idx
		//	o2=(kd<<1)+3<count?buffer[ks+ stride*3  ]:buffer[((count-1)*stride<<1)-(ks+ stride*3  )];

		//if(kd==512)//
		//	printf("fwd [%d] before B %3d  e %3d o %3d  e2 %3d o2 %3d\n", kd<<1, o1, e, o, e2, o2);//

		e-=o;		//diff
		o+=e>>1;	//av
		e2-=o2;
		o2+=e2>>1;
		char st=smoothtendency(o1, o, o2);
		e-=st;

		//if(kd==512)//
		//	printf("fwd [%d] after  B %3d  diff %3d av %3d  diff2 %3d av2 %3d,  st %3d\n", kd<<1, o1, e, o, e2, o2, st);//

		//char diff=0;
		//if(o1>=o&&o>=o2)
		//{
		//	diff=(4*o1+3*o2-o+6)/12;
		//	if(diff-(diff&1)>2*(o1-o))
		//		diff=2*(o1-o)+1;
		//	if(diff+(diff&1)>2*(o-o2))
		//		diff=2*(o-o2);
		//}
		//else if(o1<=o&&o<=o2)
		//{
		//	diff=(4*o1+3*o2-o-6)/12;
		//	if(diff+(diff&1)>2*(o1-o))
		//		diff=2*(o1-o)-1;
		//	if(diff-(diff&1)>2*(o-o2))
		//		diff=2*(o-o2);
		//}
		//e-=diff;
		
		if(kd<nodd)
			odd[kd]=o;
		even[kd]=e;
	}
	
	//for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	//{
	//	even[k]=buffer[ks];
	//	odd[k]=buffer[ks+stride];
	//}
	//if(extraeven)
	//	even[nodd]=buffer[stride*(count-1)];

	//for(int k=0;k<nodd;++k)//const predictor (difference)
	//	even[k]-=odd[k];
	//if(extraeven)
	//	even[nodd]-=odd[nodd-1];
	//
	//for(int k=0;k<nodd;++k)//update (O+(E-O)/2 = average)
	//	odd[k]+=even[k]>>1;
	//
	//for(int k=0;k<nodd;++k)//predict
	//	even[k]-=


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_squeeze_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	for(int ks=0, kd=0;ks<nodd+extraeven;++ks, kd+=stride<<1)
	{
		char
			o1=ks>0?buffer[kd-stride]:0,
			e =even[ks],
			o=ks<nodd?odd[ks]:0,
			o2=ks+1<nodd?odd[ks+1]:0;
		//char
		//	o1=ks>0?buffer[kd-stride]:buffer[kd+stride],
		//	e =even[ks],
		//	o=ks<nodd?odd[ks]:odd[((nodd-1)<<1)-ks],
		//	o2=ks+1<nodd?odd[ks+1]:odd[((nodd-1)<<1)-(ks+1)];

		//if(ks==512)//
		//	printf("inv [%d] before B %3d  diff %3d av %3d  av2 %3d\n", ks<<1, o1, e, o, o2);//

		char st=smoothtendency(o1, o, o2);
		e+=st;
		o-=e>>1;
		e+=o;

		//if(ks==512)//
		//	printf("inv [%d] after  B %3d  e %3d o %3d  av2 %3d,  st %3d\n", ks<<1, o1, e, o, o2, st);//

		buffer[kd]=e;
		if(ks<nodd)
			buffer[kd+stride]=o;
	}
	
	
	//for(int k=0;k<nodd;++k)//unupdate
	//	odd[k]-=even[k]>>1;
	//
	//for(int k=0;k<nodd;++k)//unpredict
	//	even[k]+=odd[k];
	//if(extraeven)
	//	even[nodd]+=odd[nodd-1];


	//for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	//{
	//	buffer[ks]=even[k];
	//	buffer[ks+stride]=odd[k];
	//}
	//if(extraeven)
	//	buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_squeeze_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
		{
			dwt1d_squeeze_fwd(buffer+rowlen*ky, w2, stride, temp);
			dwt1d_squeeze_inv(buffer+rowlen*ky, w2, stride, temp);//
			dwt1d_squeeze_fwd(buffer+rowlen*ky, w2, stride, temp);//
		}

		//save_channel(buffer, iw, ih, 4, 128, "squeeze-stage%02dA.PNG", it);

		//for(int kx=0;kx<w2;++kx)//vertical DWT
		//	dwt1d_squeeze_fwd(buffer+stride*kx, h2, rowlen, temp);

		//save_channel(buffer, iw, ih, 4, 128, "squeeze-stage%02dB.PNG", it);
	}
}
void dwt2d_squeeze_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		//for(int kx=0;kx<w2;++kx)//vertical IDWT
		//	dwt1d_squeeze_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_squeeze_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit CDF 5/3 DWT
void dwt1d_cdf53_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

	//cubic prediction
#if 0
	for(int k=0;k<nodd+extraeven;++k)
	{
		char
			v0=k>=2?odd[k-2]:odd[2-k],
			v1=k>=1?odd[k-1]:odd[1-k],
			v2=k<nodd?odd[k]:odd[k-1],
			v3=k+1<nodd?odd[k+1]:odd[k-2];
		even[k]-=((v1+v2)*9-(v0+v3))>>4;
	}
	//for(int k=0;k<nodd;++k)
	//{
	//	char
	//		v0=k>=2?even[k-1]:even[1-k],
	//		v1=even[k],
	//		v2=k+1<nodd+extraeven?even[k+1]:even[k-1],
	//		v3=k+2<nodd+extraeven?even[k+2]:even[k-2];
	//	odd[k]+=((v1+v2)*9-(v0+v3))>>5;
	//}
#endif

	//linear prediction
#if 1
	even[0]-=odd[0];
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]>>1;
#endif


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_cdf53_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];

	
	for(int k=0;k<nodd-!extraeven;++k)//un-update
		odd[k]-=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]>>1;
	
	even[0]+=odd[0];
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]+=odd[nodd-1];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf53_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf53_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf53_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_cdf53_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_cdf53_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_cdf53_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit CDF 9/7 DWT
static const int cdf97_coeff[]=
{
	//'factring wavelet transforms into lifting steps' - page 19
	//'a wavelet tour of signal processing - the sparse way' - page 376
	//Lifting-based Discrete Wavelet Transform for Real-Time Signal Detection 2015-10
	 0x1960C,	//-1.58613434342059f,	//alpha
	 0x00D90,	//-0.0529801185729f,	//beta
	-0x0E206,	// 0.8829110755309f,	//gamma
	-0x07189,	// 0.4435068520439f,	//delta
	 0x1264C,	// 1.1496043988602f,	//zeta		output gain is 1.89
};
static void dwt1d_u8_predict(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	even[0]-=odd[0]*coeff>>15;
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(odd[k-1]+odd[k])*coeff>>16;
	if(extraeven)
		even[nodd]-=odd[nodd-1]*coeff>>15;
}
static void dwt1d_u8_update(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(even[k]+even[k+1])*coeff>>16;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]*coeff>>15;
}
static void dwt1d_u8_unpredict(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	even[0]+=odd[0]*coeff>>15;
	for(int k=1;k<nodd;++k)//unpredict
		even[k]+=(odd[k-1]+odd[k])*coeff>>16;
	if(extraeven)
		even[nodd]+=odd[nodd-1]*coeff>>15;
}
static void dwt1d_u8_unupdate(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//unupdate
		odd[k]-=(even[k]+even[k+1])*coeff>>16;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]*coeff>>15;
}
static void dwt1d_u8_scale(char *buf, int count, int coeff)
{
	for(int k=0;k<count;++k)
		buf[k]=buf[k]*coeff>>16;
}
static void dwt1d_u8_unscale(char *buf, int count, int coeff)
{
	for(int k=0;k<count;++k)
		buf[k]=(buf[k]<<16)/coeff;
}
void dwt1d_cdf97_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

	dwt1d_u8_predict(odd, even, nodd, extraeven, cdf97_coeff[0]);
	dwt1d_u8_update (odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_u8_predict(odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_u8_update (odd, even, nodd, extraeven, cdf97_coeff[3]);
	//dwt1d_u8_scale(b2, count, cdf97_coeff[4]);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_cdf97_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	//dwt1d_u8_unscale(b2, count, cdf97_coeff[4]);
	dwt1d_u8_unupdate (odd, even, nodd, extraeven, cdf97_coeff[3]);
	dwt1d_u8_unpredict(odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_u8_unupdate (odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_u8_unpredict(odd, even, nodd, extraeven, cdf97_coeff[0]);

	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf97_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf97_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dA.PNG", it);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf97_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_cdf97_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_cdf97_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_cdf97_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}


//lifting-based 9/7 DWT
ArrayHandle dwt2d_gensizes(int iw, int ih, int wstop, int hstop, int nstages_override)//calculate dimensions of each DWT stage in descending order
{
	ArrayHandle sizes;
	int lw=floor_log2(iw), lh=floor_log2(ih), lmin=MINVAR(lw, lh);
	ARRAY_ALLOC(DWTSize, sizes, 0, 0, lmin, 0);
	if(wstop<3)
		wstop=3;
	if(hstop<3)
		hstop=3;
	int nstages=0;
	DWTSize *p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);
	p->w=iw;
	p->h=ih;
	for(int w2=iw, h2=ih;w2>=wstop&&h2>=hstop&&(!nstages_override||nstages<nstages_override);++nstages)
	{
		p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);

		w2>>=1;//w=floor(w/2)
		h2>>=1;//h=floor(h/2)
		//w2-=w2>>1;//w=ceil(w/2)
		//h2-=h2>>1;//h=ceil(h/2)

		p->w=w2;
		p->h=h2;
	}
	return sizes;
}
static const float dwt97_coeff[]=
{
	//'factring wavelet transforms into lifting steps' - page 19
	//'a wavelet tour of signal processing - the sparse way' - page 376
	-1.58613434342059f,	//alpha
	-0.0529801185729f,	//beta
	 0.8829110755309f,	//gamma
	 0.4435068520439f,	//delta
	 1.1496043988602f,	//zeta		output gain is 1.89
};
void        dwt1d_ps_scale(float *odd, float *even, int nodd, int neven, float ce, float co)
{
	int k, c2;
	__m128 factor;

	factor=_mm_set1_ps(ce);
	for(k=0, c2=nodd-3;k<c2;k+=4)//x[2k] *= ce
	{
		__m128 val=_mm_load_ps(odd+k);
		val=_mm_mul_ps(val, factor);
		_mm_store_ps(odd+k, val);
	}
	for(;k<nodd;++k)
		odd[k]*=ce;

	factor=_mm_set1_ps(co);
	for(k=0, c2=neven-3;k<c2;k+=4)//x[2k+1] *= co
	{
		__m128 val=_mm_load_ps(even+k);
		val=_mm_mul_ps(val, factor);
		_mm_store_ps(even+k, val);
	}
	for(;k<neven;++k)
		even[k]*=co;
}
void        dwt1d_ps_predict_next(float *odd, float *even, int nodd, int neven, float z10)
{
	int k, c2;
	__m128 f0;
	__m128i sh1, sh2;
	int extraeven=nodd<neven;

	//  _mm_set_epi8(15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0);
	sh1=_mm_set_epi8(-1, -1, -1, -1, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4);
	sh2=_mm_set_epi8( 3,  2,  1,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	f0=_mm_set1_ps(z10);

	//c2=neven-4;//same
	//c2=nodd-3-!extraeven;

	for(k=0, c2=neven-4;k<c2;k+=4)//x[2k+1] += z10*(x[2k]+x[2k+2])		odd[k]+=z10*(even[k]+even[k+1])
	{
		__m128 vo=_mm_load_ps(odd+k);
		__m128 ve0=_mm_load_ps(even+k);
		
		__m128 ve1=_mm_load_ss(even+k+4);
		__m128 temp=_mm_castsi128_ps(_mm_shuffle_epi8(_mm_castps_si128(ve0), sh1));
		ve1=_mm_castsi128_ps(_mm_shuffle_epi8(_mm_castps_si128(ve1), sh2));
		ve1=_mm_or_ps(ve1, temp);
		
		//__m128 ve1=_mm_loadu_ps(even+k+1);//unaligned

		ve0=_mm_add_ps(ve0, ve1);
		ve0=_mm_mul_ps(ve0, f0);
		ve0=_mm_add_ps(ve0, vo);
		_mm_store_ps(odd+k, ve0);
	}
	for(;k<nodd;++k)
		odd[k]+=z10*(even[k]+even[k+(k+1<neven)]);

	//for(c2=nodd-1-!extraeven;k<c2;++k)
	//	odd[k]+=z10*(even[k]+even[k+1]);
	//--nodd;
	//odd[nodd]+=z10*(even[nodd]+even[nodd+extraeven]);//symmetric pad if no extra even

	//if(nodd==neven)//count was even: use symmetric padding
	//	odd[nodd-1]+=z10*(even[nodd-1]+even[nodd-1]);
	//else//count was odd
	//	odd[nodd-1]+=z10*(even[nodd-1]+even[nodd]);
}
void        dwt1d_ps_update_prev(float *odd, float *even, int nodd, int neven, float z0m1)
{
	int k, c2;
	__m128 f0;
	__m128i shift=_mm_set_epi8(11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0, -1, -1, -1, -1);

	//even[0]+=z0m1*(odd[0]+odd[0]);//always pad at start
	//for(k=1;k<3;++k)
	//	even[k]+=z0m1*(odd[k-1]+odd[k]);
	
	f0=_mm_set1_ps(z0m1);
	for(k=0, c2=nodd-3;k<c2;k+=4)//x[2k] += z0m1*(x[2k-1]+x[2k+1])		even[k]+=z0m1*(odd[k-1]+odd[k])
	{
		__m128 ve=_mm_load_ps(even+k);
		__m128 vo0=_mm_load_ps(odd+k);

		__m128 vo1=_mm_load_ss(odd+k-(k>0));//always pad symmetric at start
		__m128 temp=_mm_castsi128_ps(_mm_shuffle_epi8(_mm_castps_si128(vo0), shift));
		vo1=_mm_or_ps(vo1, temp);

		//__m128 vo1=_mm_loadu_ps(odd+k-1);//unaligned

		vo0=_mm_add_ps(vo0, vo1);
		vo0=_mm_mul_ps(vo0, f0);
		vo0=_mm_add_ps(vo0, ve);
		_mm_storeu_ps(even+k, vo0);
	}
	for(;k<neven;++k)
		even[k]+=z0m1*(odd[k-1]+odd[k-(k>=nodd)]);

#if 0
	__m128 f0=_mm_set1_ps(z0m1);
	int k=halfsize-5;
	for(;k>0;k-=4)//x[2k] += z0m1*(x[2k-1]+x[2k+1])
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
#endif
}
void        dwt1d_ps_fwd(float *buffer, int count, int stride, float *b2)
{
	int nodd=count>>1, neven=count-nodd;//counting from 0:  nOdd <= nEven			nOdd = nEven - (count&1)
	float *odd=b2, *even=odd+nodd;
	even=(float*)(((size_t)even+15)&~15);//16 byte alignment
	
	int k, ks;
	for(k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(count&1)//extra even if count is odd
		even[k]=buffer[ks];
	
	dwt1d_ps_scale       (odd, even, nodd, neven, dwt97_coeff[4], 1/dwt97_coeff[4]);
	dwt1d_ps_predict_next(odd, even, nodd, neven, dwt97_coeff[3]);
	dwt1d_ps_update_prev (odd, even, nodd, neven, dwt97_coeff[2]);
	dwt1d_ps_predict_next(odd, even, nodd, neven, dwt97_coeff[1]);
	dwt1d_ps_update_prev (odd, even, nodd, neven, dwt97_coeff[0]);

	for(k=0, ks=0;k<nodd;++k, ks+=stride)
		buffer[ks]=odd[k];
	for(k=0;k<neven;++k, ks+=stride)//skip align pad
		buffer[ks]=even[k];
}
void        dwt1d_ps_inv(float *buffer, int count, int stride, float *b2)
{
	int nodd=count>>1, neven=count-nodd;
	float *odd=b2, *even=b2+nodd;
	even=(float*)(((size_t)even+15)&~15);//16 byte alignment

	int k, ks;
	for(k=0, ks=0;k<nodd;++k, ks+=stride)
		odd[k]=buffer[ks];
	for(k=0;k<neven;++k, ks+=stride)//skip align pad
		even[k]=buffer[ks];
	
	dwt1d_ps_update_prev (odd, even, nodd, neven, -dwt97_coeff[0]);
	dwt1d_ps_predict_next(odd, even, nodd, neven, -dwt97_coeff[1]);
	dwt1d_ps_update_prev (odd, even, nodd, neven, -dwt97_coeff[2]);
	dwt1d_ps_predict_next(odd, even, nodd, neven, -dwt97_coeff[3]);
	dwt1d_ps_scale       (odd, even, nodd, neven, 1/dwt97_coeff[4], dwt97_coeff[4]);

	for(k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(count&1)//extra even if count is odd
		buffer[ks]=even[k];
}
void        dwt2d_ps_fwd(float *buffer, DWTSize *sizes, int nsizes)//sizes in descending order
{
	if(!nsizes)
		return;
	int iw=sizes->w, ih=sizes->h,
		tsize=MAXVAR(iw, ih)+16;//16 is for align pad
	float *temp=(float*)_mm_malloc(tsize*sizeof(float), 16);
	for(int it=0;it<nsizes-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_ps_fwd(buffer+iw*ky, w2, 1, temp);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_ps_fwd(buffer+kx, h2, iw, temp);
	}
	_mm_free(temp);
}
void        dwt2d_ps_inv(float *buffer, DWTSize *sizes, int nsizes)//sizes in descending order
{
	if(!nsizes)
		return;
	int iw=sizes->w, ih=sizes->h,
		tsize=MAXVAR(iw, ih)+16;
	float *temp=(float*)_mm_malloc(tsize*sizeof(float), 16);
	for(int it=nsizes-2;it>=0;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_ps_inv(buffer+kx, h2, iw, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_ps_inv(buffer+iw*ky, w2, 1, temp);
	}
	_mm_free(temp);
}


//Haar (squeeze) DWT: n bit -> {n, n+1, n+1, n+2} bit
void squeeze_1d_fwd(short *buffer, int count, int stride, int vmax, short *b2)//nOdd <= nEven			nOdd = nEven - (count&1)
{
	int floorhalf=count>>1;
	short *odd=b2, *even=b2+floorhalf;
	
	for(int k=0, ks=0;k<floorhalf;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(count&1)
		even[floorhalf]=buffer[stride*(count-1)];

	for(int k=0;k<floorhalf;++k)
	{
		short
			av=(even[k]+odd[k]+(even[k]>odd[k]))>>1,//average is rounded towards the even (first) value		from JPEG XL
			diff=even[k]-odd[k]+vmax;

		//int neg=diff<0;
		//diff^=-neg;
		//diff+=neg;
		//diff<<=1;
		//diff|=neg;

		//if(diff<0)
		//	LOG_ERROR("Wrong vmax %d", vmax);

		odd[k]=av;
		even[k]=diff;
	}
	//dwt2_1d_9_7_fwd(even, odd, halfsize);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void squeeze_2d_fwd(short *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int vmax, short *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			squeeze_1d_fwd(buffer+rowlen*ky, w2, stride, vmax, temp);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			squeeze_1d_fwd(buffer+stride*kx, h2, rowlen, vmax<<(kx>=(w2>>1)), temp);
	}
}
void squeeze_1d_inv(short *buffer, int count, int stride, int vmax, short *b2)
{
	int floorhalf=count>>1;//ceil(count/2)
	short *odd=b2, *even=b2+floorhalf;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];

	for(int k=0;k<floorhalf;++k)
	{
		short av=odd[k], diff=even[k]-vmax;
		
		//int neg=diff&1;
		//diff>>=1;
		//diff^=-neg;
		//diff+=neg;
		//diff|=neg<<7&-!diff;

		even[k]=((av<<1)+diff+(diff>0?-(diff&1):(diff&1)))>>1;//JPEG XL
		odd[k]=even[k]-diff;
	}
	//dwt2_1d_9_7_inv(even, odd, ceilhalf);

	for(int k=0, ks=0;k<floorhalf;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(count&1)
		buffer[stride*(count-1)]=even[floorhalf];
}
void squeeze_2d_inv(short *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int vmax, short *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			squeeze_1d_inv(buffer+stride*kx, h2, rowlen, vmax, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			squeeze_1d_inv(buffer+rowlen*ky, w2, stride, vmax, temp);
	}
}


//Haar DWT: 8 bit -> 9 bit
void haar_1d_fwd(short *buffer, int count, int stride, short *b2)//nOdd <= nEven			nOdd = nEven - (count&1)
{
	int floorhalf=count>>1,
		ceilhalf=count-floorhalf;//ceil(count/2)
	short *even=b2, *odd=b2+ceilhalf;
	
	for(int k=0, ks=0;k<floorhalf;++k, ks+=stride<<1)//lazy wavelet: split into even (low frequency) & odd (high frequency), contrary to DWT notation			TODO _mm_shuffle_epi8?
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(count&1)
		even[floorhalf]=buffer[stride*(count-1)];

	for(int k=0;k<floorhalf;++k)
	{
		short
			av=(even[k]+odd[k]+(even[k]>odd[k]))>>1,//average is rounded towards the even (first) value		from JPEG XL
			diff=even[k]-odd[k];

		//int neg=diff<0;
		//diff^=-neg;
		//diff+=neg;
		//diff<<=1;
		//diff|=neg;

		even[k]=av;
		odd[k]=diff;
	}
	//dwt2_1d_9_7_fwd(even, odd, halfsize);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void haar_2d_fwd(const unsigned char *buf, int bw, int bh, int nch, int bytestride, int nstages, short **ret)//lossless DWT, don't forget to free ret if was zero
{
	int tsize=bw>bh?bw:bh;
	short *temp=(short*)malloc(tsize*sizeof(short));

	size_t dstlen=(size_t)bw*bh*bytestride;
	if(!*ret)
	{
		*ret=(short*)malloc(dstlen*sizeof(short));
		if(!*ret)
		{
			free(temp);
			return;
		}
	}
	for(int k=0;k<dstlen;++k)
		ret[0][k]=buf[k];
	//if(*ret!=buf)
	//	memcpy(*ret, buf, dstlen);

	int rowlen=bytestride*bw;
	for(int kc=0;kc<nch;++kc)
	{
		//print_sbuf(*ret, bw, bh, kc, bytestride);//

		for(int w2=bw, h2=bh, it=0;w2>=3&&h2>=3&&(!nstages||it<nstages);++it)
		{
			for(int ky=0;ky<h2;++ky)//horizontal DWT
				haar_1d_fwd(*ret+rowlen*ky+kc, w2, bytestride, temp);
			
			//print_sbuf(*ret, bw, bh, kc, bytestride);//

			for(int kx=0;kx<w2;++kx)//vertical DWT
				haar_1d_fwd(*ret+bytestride*kx+kc, h2, rowlen, temp);
			
			//print_sbuf(*ret, bw, bh, kc, bytestride);//

			w2-=w2>>1;//w=ceil(w/2)
			h2-=h2>>1;//h=ceil(h/2)
		}
	}
	free(temp);
}

void haar_1d_inv(short *buffer, int count, int stride, short *b2)
{
	int floorhalf=count>>1, ceilhalf=count-floorhalf;//ceil(count/2)
	short *even=b2, *odd=b2+ceilhalf;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];

	for(int k=0;k<floorhalf;++k)
	{
		short av=even[k], diff=odd[k];
		
		//int neg=diff&1;
		//diff>>=1;
		//diff^=-neg;
		//diff+=neg;
		//diff|=neg<<7&-!diff;

		even[k]=((av<<1)+diff+(diff>0?-(diff&1):(diff&1)))>>1;
		odd[k]=even[k]-diff;
	}
	//dwt2_1d_9_7_inv(even, odd, ceilhalf);

	for(int k=0, ks=0;k<floorhalf;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(count&1)
		buffer[stride*(count-1)]=even[floorhalf];
}
void haar_2d_inv(short *buf, int bw, int bh, int nch, int bytestride, int nstages, unsigned char **ret)//buf is destroyed
{
	int lw=floor_log2(bw), lh=floor_log2(bh);
	int *sizes=(int*)malloc(((size_t)(lw<lh?lw:lh)<<1)*sizeof(int));
	if(!sizes)
		return;
	int nsizes=0;
	for(int w2=bw, h2=bh;w2>=3&&h2>=3&&(!nstages||nsizes<nstages);++nsizes)//calculate dimensions of each stage
	{
		sizes[nsizes<<1]=w2;
		sizes[nsizes<<1|1]=h2;
		w2-=w2>>1;//w=ceil(w/2)
		h2-=h2>>1;//h=ceil(h/2)
	}

	int tsize=maximum(bw, bh);
	short *temp=(short*)malloc(tsize*sizeof(short));
	if(!temp)
	{
		free(sizes);
		return;
	}

	int rowlen=bytestride*bw;
	for(int kc=0;kc<nch;++kc)
	{
		//print_sbuf(buf, bw, bh, kc, bytestride);//

		for(int it=nsizes-1;it>=0;--it)
		{
			int w2=sizes[it<<1], h2=sizes[it<<1|1];

			for(int kx=0;kx<w2;++kx)//vertical IDWT
				haar_1d_inv(buf+bytestride*kx+kc, h2, rowlen, temp);

			//print_sbuf(buf, bw, bh, kc, bytestride);//

			for(int ky=0;ky<h2;++ky)//horizontal IDWT
				haar_1d_inv(buf+rowlen*ky+kc, w2, bytestride, temp);

			//print_sbuf(buf, bw, bh, kc, bytestride);//
		}
	}
	free(sizes);
	free(temp);

	size_t dstlen=(size_t)bw*bh*bytestride;
	if(!*ret)
	{
		*ret=(unsigned char*)malloc(dstlen);
		if(!*ret)
			return;
	}
	for(int k=0;k<dstlen;++k)
		ret[0][k]=(unsigned char)buf[k];
}


//integer DCT-II		TODO benchmark vs FCT, where is the inverse?
void transpose_block8x8_ps(float *ptr)//32 byte aligned
{
	//AVX implementation is broken
#if 0
//#ifdef __AVX__
	//https://stackoverflow.com/questions/25622745/transpose-an-8x8-float-using-avx-avx2
	__m256
		r0=_mm256_load_ps(ptr   ),
		r4=_mm256_load_ps(ptr+ 8),
		r6=_mm256_load_ps(ptr+16),
		r2=_mm256_load_ps(ptr+24),
		r7=_mm256_load_ps(ptr+32),
		r3=_mm256_load_ps(ptr+40),
		r5=_mm256_load_ps(ptr+48),
		r1=_mm256_load_ps(ptr+56);

	__m256 t0, t1, t2, t3, t4, t5, t6, t7;
	__m256 tt0, tt1, tt2, tt3, tt4, tt5, tt6, tt7;
	t0=_mm256_unpacklo_ps(r0, r1);
	t1=_mm256_unpackhi_ps(r0, r1);
	t2=_mm256_unpacklo_ps(r2, r3);
	t3=_mm256_unpackhi_ps(r2, r3);
	t4=_mm256_unpacklo_ps(r4, r5);
	t5=_mm256_unpackhi_ps(r4, r5);
	t6=_mm256_unpacklo_ps(r6, r7);
	t7=_mm256_unpackhi_ps(r6, r7);
	tt0=_mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(1,0,1,0));
	tt1=_mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(3,2,3,2));
	tt2=_mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1,0,1,0));
	tt3=_mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(3,2,3,2));
	tt4=_mm256_shuffle_ps(t4, t6, _MM_SHUFFLE(1,0,1,0));
	tt5=_mm256_shuffle_ps(t4, t6, _MM_SHUFFLE(3,2,3,2));
	tt6=_mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(1,0,1,0));
	tt7=_mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(3,2,3,2));
	r0=_mm256_permute2f128_ps(tt0, tt4, 0x20);
	r1=_mm256_permute2f128_ps(tt1, tt5, 0x20);
	r2=_mm256_permute2f128_ps(tt2, tt6, 0x20);
	r3=_mm256_permute2f128_ps(tt3, tt7, 0x20);
	r4=_mm256_permute2f128_ps(tt0, tt4, 0x31);
	r5=_mm256_permute2f128_ps(tt1, tt5, 0x31);
	r6=_mm256_permute2f128_ps(tt2, tt6, 0x31);
	r7=_mm256_permute2f128_ps(tt3, tt7, 0x31);

	_mm256_store_ps((float*)ptr   , r0);
	_mm256_store_ps((float*)ptr+ 8, r1);
	_mm256_store_ps((float*)ptr+16, r2);
	_mm256_store_ps((float*)ptr+24, r3);
	_mm256_store_ps((float*)ptr+32, r4);
	_mm256_store_ps((float*)ptr+40, r5);
	_mm256_store_ps((float*)ptr+48, r6);
	_mm256_store_ps((float*)ptr+56, r7);
#else
	__m128
		a0lo=_mm_load_ps(ptr   ), a0hi=_mm_load_ps(ptr+ 4),
		a1lo=_mm_load_ps(ptr+ 8), a1hi=_mm_load_ps(ptr+12),
		a2lo=_mm_load_ps(ptr+16), a2hi=_mm_load_ps(ptr+20),
		a3lo=_mm_load_ps(ptr+24), a3hi=_mm_load_ps(ptr+28),
		a4lo=_mm_load_ps(ptr+32), a4hi=_mm_load_ps(ptr+36),
		a5lo=_mm_load_ps(ptr+40), a5hi=_mm_load_ps(ptr+44),
		a6lo=_mm_load_ps(ptr+48), a6hi=_mm_load_ps(ptr+52),
		a7lo=_mm_load_ps(ptr+56), a7hi=_mm_load_ps(ptr+60);

	_MM_TRANSPOSE4_PS(a0lo, a1lo, a2lo, a3lo);
	_MM_TRANSPOSE4_PS(a0hi, a1hi, a2hi, a3hi);
	_MM_TRANSPOSE4_PS(a4lo, a5lo, a6lo, a7lo);
	_MM_TRANSPOSE4_PS(a4hi, a5hi, a6hi, a7hi);
	
	_mm_store_ps(ptr   , a0lo); _mm_store_ps(ptr+ 4, a4lo);
	_mm_store_ps(ptr+ 8, a1lo); _mm_store_ps(ptr+12, a5lo);
	_mm_store_ps(ptr+16, a2lo); _mm_store_ps(ptr+20, a6lo);
	_mm_store_ps(ptr+24, a3lo); _mm_store_ps(ptr+28, a7lo);
	_mm_store_ps(ptr+32, a0hi); _mm_store_ps(ptr+36, a4hi);
	_mm_store_ps(ptr+40, a1hi); _mm_store_ps(ptr+44, a5hi);
	_mm_store_ps(ptr+48, a2hi); _mm_store_ps(ptr+52, a6hi);
	_mm_store_ps(ptr+56, a3hi); _mm_store_ps(ptr+60, a7hi);
#endif
}
void intDCT2_8x8_p32_rows(int *buf)//buf is aligned packed 32 bit as fixed precision
{
	//https://fgiesen.wordpress.com/2013/11/04/bink-2-2-integer-dct-design-part-1/
	__m128i half=_mm_set1_epi32(1<<5);

	__m128i
		a0lo=_mm_load_si128((__m128i*)buf+ 0), a0hi=_mm_load_si128((__m128i*)buf+ 1),
		a1lo=_mm_load_si128((__m128i*)buf+ 2), a1hi=_mm_load_si128((__m128i*)buf+ 3),
		a2lo=_mm_load_si128((__m128i*)buf+ 4), a2hi=_mm_load_si128((__m128i*)buf+ 5),
		a3lo=_mm_load_si128((__m128i*)buf+ 6), a3hi=_mm_load_si128((__m128i*)buf+ 7),
		a4lo=_mm_load_si128((__m128i*)buf+ 8), a4hi=_mm_load_si128((__m128i*)buf+ 9),
		a5lo=_mm_load_si128((__m128i*)buf+10), a5hi=_mm_load_si128((__m128i*)buf+11),
		a6lo=_mm_load_si128((__m128i*)buf+12), a6hi=_mm_load_si128((__m128i*)buf+13),
		a7lo=_mm_load_si128((__m128i*)buf+14), a7hi=_mm_load_si128((__m128i*)buf+15);

	__m128i
		b0lo=_mm_add_epi32(a0lo, a7lo), b0hi=_mm_add_epi16(a0hi, a7hi),
		b1lo=_mm_add_epi32(a1lo, a6lo), b1hi=_mm_add_epi32(a1hi, a6hi),
		b2lo=_mm_add_epi32(a2lo, a5lo), b2hi=_mm_add_epi32(a2hi, a5hi),
		b3lo=_mm_add_epi32(a3lo, a4lo), b3hi=_mm_add_epi32(a3hi, a4hi),
		b4lo=_mm_sub_epi32(a3lo, a4lo), b4hi=_mm_sub_epi32(a3hi, a4hi),
		b5lo=_mm_sub_epi32(a2lo, a5lo), b5hi=_mm_sub_epi32(a2hi, a5hi),
		b6lo=_mm_sub_epi32(a1lo, a6lo), b6hi=_mm_sub_epi32(a1hi, a6hi),
		b7lo=_mm_sub_epi32(a0lo, a7lo), b7hi=_mm_sub_epi32(a0hi, a7hi);

#define MUL_C7(REG)		_mm_add_epi32(REG, _mm_add_epi32(_mm_slli_epi32(REG, 2), _mm_slli_epi32(REG, 3)))
#define MUL_S7(REG)		_mm_add_epi32(REG, _mm_slli_epi32(REG, 6))
#define MUL_C5(REG)		_mm_add_epi32(REG, _mm_add_epi32(_mm_slli_epi32(REG, 2), _mm_slli_epi32(REG, 5)))
#define MUL_S5(REG)		_mm_sub_epi32(_mm_slli_epi32(REG, 6), _mm_add_epi32(_mm_slli_epi32(REG, 3), REG))
#define MUL_C6(REG)		_mm_sub_epi32(_mm_slli_epi32(REG, 5), _mm_add_epi32(_mm_slli_epi32(REG, 1), REG))
#define MUL_S6(REG)		_mm_add_epi32(_mm_slli_epi32(REG, 6), _mm_add_epi32(_mm_slli_epi32(REG, 2), _mm_slli_epi32(REG, 1)))
#define ROUND(REG)		_mm_srai_epi32(_mm_add_epi32(REG, half), 6)
	a0lo=_mm_add_epi32(b0lo, b3lo); a0hi=_mm_add_epi32(b0hi, b3hi);
	a1lo=_mm_add_epi32(b1lo, b2lo); a1hi=_mm_add_epi32(b1hi, b2hi);
	a2lo=_mm_sub_epi32(b1lo, b2lo); a2hi=_mm_sub_epi32(b1hi, b2hi);
	a3lo=_mm_sub_epi32(b0lo, b3lo); a3hi=_mm_sub_epi32(b0hi, b3hi);
	a4lo=_mm_add_epi32(MUL_C7(b4lo), MUL_S7(b7lo)); a4hi=_mm_add_epi32(MUL_C7(b4hi), MUL_S7(b7hi));//c7*b4+s7*b7
	a5lo=_mm_add_epi32(MUL_C5(b5lo), MUL_S5(b6lo)); a5hi=_mm_add_epi32(MUL_C5(b5hi), MUL_S5(b6hi));//c5*b5+s5*b6
	a6lo=_mm_sub_epi32(MUL_C5(b6lo), MUL_S7(b5lo)); a6hi=_mm_sub_epi32(MUL_C5(b6hi), MUL_S5(b5hi));//c5*b6-s5*b5
	a7lo=_mm_sub_epi32(MUL_C7(b7lo), MUL_S5(b4lo)); a7hi=_mm_sub_epi32(MUL_C7(b7hi), MUL_S7(b4hi));//c7*b7-s7*b4

	b0lo=_mm_add_epi32(a0lo, a1lo); b0hi=_mm_add_epi32(a0hi, a1hi);
	b1lo=_mm_sub_epi32(a0lo, a1lo); b1hi=_mm_sub_epi32(a0hi, a1hi);
	b2lo=_mm_add_epi32(MUL_C6(a2lo), MUL_S6(a3lo)); b2hi=_mm_add_epi32(MUL_C6(a2hi), MUL_S6(a3hi));//c6*a2+s6*a3
	b3lo=_mm_sub_epi32(MUL_C6(a3lo), MUL_S6(a2lo)); b3hi=_mm_sub_epi32(MUL_C6(a3hi), MUL_S6(a2hi));//c6*a3-s6*a2
	b4lo=_mm_add_epi32(a4lo, a5lo); b4hi=_mm_add_epi32(a4hi, a5hi);
	b5lo=_mm_sub_epi32(a4lo, a5lo); b5hi=_mm_sub_epi32(a4hi, a5hi);
	b6lo=_mm_add_epi32(a6lo, a7lo); b6hi=_mm_add_epi32(a6hi, a7hi);
	b7lo=_mm_sub_epi32(a6lo, a7lo); b7hi=_mm_sub_epi32(a6hi, a7hi);

	b2lo=ROUND(b2lo); b2hi=ROUND(b2hi);//divide by 64
	b3lo=ROUND(b3lo); b3hi=ROUND(b3hi);
	b4lo=ROUND(b4lo); b4hi=ROUND(b4hi);
	b5lo=ROUND(b5lo); b5hi=ROUND(b5hi);
	b6lo=ROUND(b6lo); b6hi=ROUND(b6hi);
	b7lo=ROUND(b7lo); b7hi=ROUND(b7hi);

	a5lo=_mm_add_epi32(b5lo, b6lo); a5hi=_mm_add_epi32(b5hi, b6hi);
	a6lo=_mm_sub_epi32(b5lo, b6lo); a6hi=_mm_sub_epi32(b5hi, b6hi);

	_mm_store_si128((__m128i*)buf+ 0, b0lo); _mm_store_si128((__m128i*)buf+ 1, b0hi);
	_mm_store_si128((__m128i*)buf+ 2, b4lo); _mm_store_si128((__m128i*)buf+ 3, b4hi);
	_mm_store_si128((__m128i*)buf+ 4, b2lo); _mm_store_si128((__m128i*)buf+ 5, b2hi);
	_mm_store_si128((__m128i*)buf+ 6, a6lo); _mm_store_si128((__m128i*)buf+ 7, a6hi);
	_mm_store_si128((__m128i*)buf+ 8, b1lo); _mm_store_si128((__m128i*)buf+ 9, b1hi);
	_mm_store_si128((__m128i*)buf+10, b3lo); _mm_store_si128((__m128i*)buf+11, b3hi);
	_mm_store_si128((__m128i*)buf+12, a5lo); _mm_store_si128((__m128i*)buf+13, a5hi);
	_mm_store_si128((__m128i*)buf+14, b7lo); _mm_store_si128((__m128i*)buf+15, b7hi);
#undef	MUL_C7
#undef	MUL_S7
#undef	MUL_C5
#undef	MUL_S5
#undef	MUL_C6
#undef	MUL_S6
}
void intDCT2_8x8_p32(const unsigned char *buf, int bw, int bh, int nch, int bytestride, unsigned short **ret)
{
	if(bw&7||bh&7)
		return;
	size_t len=(size_t)bw*bh*bytestride;
	if(!*ret)
	{
		*ret=(unsigned short*)malloc(len*sizeof(short));
		if(!*ret)
			return;
	}
	ALIGN(16) int temp[64];
	for(int ky=0;ky<bh;ky+=8)
	{
		for(int kx=0;kx<bw;kx+=8)
		{
			for(int kc=0;kc<nch;++kc)
			{
				for(int ky2=0;ky2<8;++ky2)
				{
					for(int kx2=0;kx2<8;++kx2)
						temp[ky2<<3|kx2]=buf[bytestride*(bw*(ky|ky2)+(kx|kx2))+kc];
					intDCT2_8x8_p32_rows(temp);
					transpose_block8x8_ps((float*)temp);
					intDCT2_8x8_p32_rows(temp);
					transpose_block8x8_ps((float*)temp);
					for(int kx2=0;kx2<8;++kx2)
						ret[0][bytestride*(bw*(ky|ky2)+(kx|kx2))+kc]=temp[ky2<<3|kx2];
				}
			}
		}
	}
}
#if 0
void DCTtest()
{
	ALIGN(16) int block[64];
	int b2[64];
	long long result[64];
	for(int ky=0;ky<8;++ky)
	{
		for(int kx=0;kx<8;++kx)
			block[ky<<3|kx]=ky==kx;
	}
	intDCT2_8x8_p32_rows(block);
	memcpy(b2, block, 64*sizeof(int));
	transpose_block8x8_ps((float*)block);

	for(int ky=0;ky<8;++ky)//matmul
	{
		for(int kx=0;kx<8;++kx)
		{
			long long sum=0;
			for(int ki=0;ki<8;++ki)
				sum+=(long long)b2[ky<<3|ki]*block[ki<<3|kx];
			result[ky<<3|kx]=sum;
		}
	}
	for(int ky=0;ky<8;++ky)
	{
		for(int kx=0;kx<8;++kx)
			printf("%lld\t", result[ky<<3|kx]);
		printf("\n");
	}
}
#endif


//	#define DCT_8x8_PS_SCALED

//fast 8 point DCT-II/III, call transpose then call again, ptr to 64 packed floats is aligned by 32 bytes
void DCT2_8x8_ps(float *ptr)
{
	//Loeffler's factorization using lifting		http://www.sfu.ca/~jiel/papers/c003-bindct-ieee.pdf
#ifdef __AVX__
	__m256
		half=_mm256_set1_ps(0.5f),
		p1=_mm256_set1_ps(13.f/32),//binDCT-L5
		u1=_mm256_set1_ps(11.f/32),
		p2=_mm256_set1_ps(19.f/64),
		u2=_mm256_set1_ps( 9.f/16),
		p3=_mm256_set1_ps(19.f/64),
		p4=_mm256_set1_ps( 3.f/32),
		u3=_mm256_set1_ps( 3.f/16),
		p5=_mm256_set1_ps( 3.f/32);
	
	__m256
		x0=_mm256_load_ps(ptr   ),
		x1=_mm256_load_ps(ptr+ 8),
		x2=_mm256_load_ps(ptr+16),
		x3=_mm256_load_ps(ptr+24),
		x4=_mm256_load_ps(ptr+32),
		x5=_mm256_load_ps(ptr+40),
		x6=_mm256_load_ps(ptr+48),
		x7=_mm256_load_ps(ptr+56);
	
	//new stage after each extra newline
	x7=_mm256_sub_ps(x0, x7);
	x6=_mm256_sub_ps(x1, x6);
	x5=_mm256_sub_ps(x2, x5);
	x4=_mm256_sub_ps(x3, x4);

	x0=_mm256_fnmadd_ps(x7, half, x0);
	x1=_mm256_fnmadd_ps(x6, half, x1);
	x2=_mm256_fnmadd_ps(x5, half, x2);
	x3=_mm256_fnmadd_ps(x4, half, x3);

	x3=_mm256_sub_ps(x0, x3);
	x2=_mm256_sub_ps(x1, x2);
	x7=_mm256_fnmadd_ps(x4, p2, x7);
	x6=_mm256_fnmadd_ps(x5, p4, x6);
	
	x0=_mm256_fnmadd_ps(x3, half, x0);
	x1=_mm256_fnmadd_ps(x2, half, x1);
	x4=_mm256_fmadd_ps(x7, u2, x4);
	x5=_mm256_fmadd_ps(x6, u3, x5);

	x0=_mm256_add_ps(x0, x1);
	x2=_mm256_fmsub_ps(x3, p1, x2);
	x7=_mm256_fnmadd_ps(x4, p3, x7);
	x6=_mm256_fnmadd_ps(x5, p5, x6);

	x1=_mm256_fmsub_ps(x0, half, x1);
	x3=_mm256_fnmadd_ps(x2, u1, x3);
	x6=_mm256_sub_ps(x4, x6);
	x5=_mm256_sub_ps(x7, x5);

	x4=_mm256_fnmadd_ps(x6, half, x4);
	x7=_mm256_fnmadd_ps(x5, half, x7);

	x7=_mm256_add_ps(x7, x4);
	
	x4=_mm256_fmsub_ps(x7, half, x4);
	
	_mm256_store_ps(ptr   , x0);
	_mm256_store_ps(ptr+ 8, x4);
	_mm256_store_ps(ptr+16, x6);
	_mm256_store_ps(ptr+24, x2);
	_mm256_store_ps(ptr+32, x7);
	_mm256_store_ps(ptr+40, x3);
	_mm256_store_ps(ptr+48, x5);
	_mm256_store_ps(ptr+56, x1);
#else
#error TODO
#endif

#if 0
	//Chen's factorization using lifting		http://www.sfu.ca/~jiel/papers/c003-bindct-ieee.pdf
#ifdef __AVX__
	__m256
		half=_mm256_set1_ps(0.5f),
		p1=_mm256_set1_ps(13.f/32),//binDCT-C7
		u1=_mm256_set1_ps(11.f/32),
		p2=_mm256_set1_ps(11.f/16),
		u2=_mm256_set1_ps(15.f/32),
		p3=_mm256_set1_ps( 3.f/16),
		u3=_mm256_set1_ps( 3.f/16),
		p4=_mm256_set1_ps(13.f/32),
		u4=_mm256_set1_ps(11.f/16),
		p5=_mm256_set1_ps(13.f/32);
	
	__m256
		x0=_mm256_load_ps(ptr   ),
		x1=_mm256_load_ps(ptr+ 8),
		x2=_mm256_load_ps(ptr+16),
		x3=_mm256_load_ps(ptr+24),
		x4=_mm256_load_ps(ptr+32),
		x5=_mm256_load_ps(ptr+40),
		x6=_mm256_load_ps(ptr+48),
		x7=_mm256_load_ps(ptr+56);

	//new stage after each extra newline
	x7=_mm256_sub_ps(x0, x7);
	x6=_mm256_sub_ps(x1, x6);
	x5=_mm256_sub_ps(x2, x5);
	x4=_mm256_sub_ps(x3, x4);

	x0=_mm256_sub_ps(x0, _mm256_mul_ps(x7, half));
	x1=_mm256_sub_ps(x1, _mm256_mul_ps(x6, half));
	x2=_mm256_sub_ps(x2, _mm256_mul_ps(x5, half));
	x3=_mm256_sub_ps(x3, _mm256_mul_ps(x4, half));

	x3=_mm256_sub_ps(x0, x3);
	x2=_mm256_sub_ps(x1, x2);
	x5=_mm256_sub_ps(x5, _mm256_mul_ps(x6, p4));

	x0=_mm256_sub_ps(x0, _mm256_mul_ps(x3, half));
	x1=_mm256_sub_ps(x1, _mm256_mul_ps(x2, half));
	x6=_mm256_add_ps(x6, _mm256_mul_ps(x5, u4));

	x5=_mm256_sub_ps(_mm256_mul_ps(x6, p5), x5);

	x5=_mm256_sub_ps(x4, x5);
	x6=_mm256_sub_ps(x7, x6);
	x1=_mm256_sub_ps(x0, x1);
	x2=_mm256_sub_ps(_mm256_mul_ps(x3, p1), x2);

	x4=_mm256_sub_ps(x4, _mm256_mul_ps(x5, half));
	x7=_mm256_sub_ps(x7, _mm256_mul_ps(x6, half));

	x4=_mm256_sub_ps(_mm256_mul_ps(x7, p3), x4);
	x5=_mm256_add_ps(_mm256_mul_ps(x6, p2), x5);
	
	x0=_mm256_sub_ps(x0, _mm256_mul_ps(x1, half));
	x2=_mm256_sub_ps(x2, _mm256_mul_ps(x1, half));

#else
	__m128
		half=_mm_set1_ps(0.5f),
		p1=_mm_set1_ps(13.f/32),//binDCT-C7
		u1=_mm_set1_ps(11.f/32),
		p2=_mm_set1_ps(11.f/16),
		u2=_mm_set1_ps(15.f/32),
		p3=_mm_set1_ps( 3.f/16),
		u3=_mm_set1_ps( 3.f/16),
		p4=_mm_set1_ps(13.f/32),
		u4=_mm_set1_ps(11.f/16),
		p5=_mm_set1_ps(13.f/32);
	
	__m128
		x0lo=_mm_load_ps(ptr   ), x0hi=_mm_load_ps(ptr+ 4),
		x1lo=_mm_load_ps(ptr+ 8), x1hi=_mm_load_ps(ptr+12),
		x2lo=_mm_load_ps(ptr+16), x2hi=_mm_load_ps(ptr+20),
		x3lo=_mm_load_ps(ptr+24), x3hi=_mm_load_ps(ptr+28),
		x4lo=_mm_load_ps(ptr+32), x4hi=_mm_load_ps(ptr+36),
		x5lo=_mm_load_ps(ptr+40), x5hi=_mm_load_ps(ptr+44),
		x6lo=_mm_load_ps(ptr+48), x6hi=_mm_load_ps(ptr+52),
		x7lo=_mm_load_ps(ptr+56), x7hi=_mm_load_ps(ptr+60);

	x7lo=_mm_sub_ps(x0lo, x7lo), x7hi=_mm_sub_ps(x0hi, x7hi);
#endif
#endif

#if 0
	//Loeffler's factorization using lifting		http://www.sfu.ca/~jiel/papers/c003-bindct-ieee.pdf
	__m128
		half=_mm_set1_ps(0.5f),
		p1=_mm_set1_ps(13.f/32),//binDCT-L5
		u1=_mm_set1_ps(11.f/32),
		p2=_mm_set1_ps(19.f/64),
		u2=_mm_set1_ps( 9.f/16),
		p3=_mm_set1_ps(19.f/64),
		p4=_mm_set1_ps( 3.f/32),
		u3=_mm_set1_ps( 3.f/16),
		p5=_mm_set1_ps( 3.f/32);
	
#ifdef DCT_8x8_PS_SCALED
	__m128
		g0=_mm_set1_ps(0.3535533905932737622004221810524f),//sin(pi/4)/2
		g1=_mm_set1_ps(0.7071067811865475244008443621048f),//sin(pi/4)
		g2=_mm_set1_ps(0.4619397662556433780640915946984f),//sin(3pi/8)/2
		g3=_mm_set1_ps(0.5411961001461969843997232053664f),//1/(2sin(3pi/8))
		g4=_mm_set1_ps(0.7071067811865475244008443621048f),//1/sqrt2
		g5=half,//1/2
		g6=half,//1/2
		g7=_mm_set1_ps(0.3535533905932737622004221810524f);//1/sqrt8
#endif

	__m128
		a0lo=_mm_load_ps(ptr   ), a0hi=_mm_load_ps(ptr+ 4),
		a1lo=_mm_load_ps(ptr+ 8), a1hi=_mm_load_ps(ptr+12),
		a2lo=_mm_load_ps(ptr+16), a2hi=_mm_load_ps(ptr+20),
		a3lo=_mm_load_ps(ptr+24), a3hi=_mm_load_ps(ptr+28),
		a4lo=_mm_load_ps(ptr+32), a4hi=_mm_load_ps(ptr+36),
		a5lo=_mm_load_ps(ptr+40), a5hi=_mm_load_ps(ptr+44),
		a6lo=_mm_load_ps(ptr+48), a6hi=_mm_load_ps(ptr+52),
		a7lo=_mm_load_ps(ptr+56), a7hi=_mm_load_ps(ptr+60);

	__m128
		b0lo=_mm_add_ps(a0lo, a7lo), b0hi=_mm_add_ps(a0hi, a7hi),
		b1lo=_mm_add_ps(a1lo, a6lo), b1hi=_mm_add_ps(a1hi, a6hi),
		b2lo=_mm_add_ps(a2lo, a5lo), b2hi=_mm_add_ps(a2hi, a5hi),
		b3lo=_mm_add_ps(a3lo, a4lo), b3hi=_mm_add_ps(a3hi, a4hi),
		b4lo=_mm_sub_ps(a0lo, a7lo), b4hi=_mm_sub_ps(a0hi, a7hi),
		b5lo=_mm_sub_ps(a1lo, a6lo), b5hi=_mm_sub_ps(a1hi, a6hi),
		b6lo=_mm_sub_ps(a2lo, a5lo), b6hi=_mm_sub_ps(a2hi, a5hi),
		b7lo=_mm_sub_ps(a3lo, a4lo), b7hi=_mm_sub_ps(a3hi, a4hi);

	__m128 t0lo, t0hi;

	t0lo=_mm_mul_ps(b4lo, p2  ), t0hi=_mm_mul_ps(b4hi, p2  );//pred2:   x[7]-=p2*x[4]
	b7lo=_mm_sub_ps(b7lo, t0lo), b7hi=_mm_sub_ps(b7hi, t0hi);

	b7lo=_mm_mul_ps(b7lo, u2  ), b7hi=_mm_mul_ps(b7hi, u2  );//update2: x[4]+=u2*x[7]
	b4lo=_mm_add_ps(b4lo, b7lo), b4hi=_mm_add_ps(b4hi, b7hi);

	t0lo=_mm_mul_ps(b4lo, p3  ), t0hi=_mm_mul_ps(b4hi, p3  );//pred3:   x[7]-=p3*x[4]
	b7lo=_mm_sub_ps(b7lo, t0lo), b7hi=_mm_sub_ps(b7hi, t0hi);


	t0lo=_mm_mul_ps(b5lo, p4  ), t0hi=_mm_mul_ps(b5hi, p4  );//pred4:   x[6]-=p4*x[5]
	b6lo=_mm_sub_ps(b6lo, t0lo), b6hi=_mm_sub_ps(b6hi, t0hi);

	t0lo=_mm_mul_ps(b6lo, u3  ), t0hi=_mm_mul_ps(b6hi, u3  );//update3: x[5]+=u3*x[6]
	b5lo=_mm_add_ps(b5lo, t0lo), b5hi=_mm_add_ps(b5hi, t0hi);

	t0lo=_mm_mul_ps(b5lo, p5  ), t0hi=_mm_mul_ps(b5hi, p5  );//pred5:   x[6]-=p5*x[5]
	b6lo=_mm_sub_ps(b6lo, t0lo), b6hi=_mm_sub_ps(b6hi, t0hi);

	a0lo=_mm_add_ps(b0lo, b3lo), a0hi=_mm_add_ps(b0hi, b3hi);
	a1lo=_mm_add_ps(b1lo, b2lo), a1hi=_mm_add_ps(b1hi, b2hi);
	a2lo=_mm_sub_ps(b0lo, b3lo), a2hi=_mm_sub_ps(b0hi, b3hi);
	a3lo=_mm_sub_ps(b1lo, b2lo), a3hi=_mm_sub_ps(b1hi, b2hi);
	a4lo=_mm_add_ps(b4lo, b6lo), a4hi=_mm_add_ps(b4hi, b6hi);
	a5lo=_mm_add_ps(b5lo, b7lo), a5hi=_mm_add_ps(b5hi, b7hi);
	a6lo=_mm_sub_ps(b4lo, b6lo), a6hi=_mm_sub_ps(b4hi, b6hi);
	a7lo=_mm_sub_ps(b5lo, b7lo), a7hi=_mm_sub_ps(b5hi, b7hi);

	a0lo=_mm_add_ps(a0lo, a1lo), a0hi=_mm_add_ps(a0hi, a1hi);//x[0]+=x[1]
	a7lo=_mm_add_ps(a7lo, a4lo), a7hi=_mm_add_ps(a7hi, a4hi);//x[7]+=x[4]

	t0lo=_mm_mul_ps(a0lo, half), t0hi=_mm_mul_ps(a0hi, half);//x[1] = x[0]/2 - x[1]
	a1lo=_mm_sub_ps(t0lo, a1lo), a1hi=_mm_sub_ps(t0hi, a1hi);

	t0lo=_mm_mul_ps(a3lo, p1  ), t0hi=_mm_mul_ps(a3hi, p1  );//pred1: x[2] = p1*x[3] - x[2]
	a2lo=_mm_sub_ps(t0lo, a2lo), a2hi=_mm_sub_ps(t0hi, a2hi);

	t0lo=_mm_mul_ps(a2lo, u1  ), t0hi=_mm_mul_ps(a2hi, u1  );//update1: x[3]-=u1*x[2]
	a3lo=_mm_sub_ps(a3lo, t0lo), a3hi=_mm_sub_ps(a3hi, t0hi);

	t0lo=_mm_mul_ps(a7lo, half), t0hi=_mm_mul_ps(a7hi, half);//x[4] = x[7]/2 - x[4]
	a4lo=_mm_sub_ps(t0lo, a4lo), a4hi=_mm_sub_ps(t0hi, a4hi);

#ifdef DCT_8x8_PS_SCALED
	a0lo=_mm_mul_ps(a0lo, g0), a0hi=_mm_mul_ps(a0hi, g0);
	a1lo=_mm_mul_ps(a1lo, g1), a1hi=_mm_mul_ps(a1hi, g1);
	a2lo=_mm_mul_ps(a2lo, g2), a2hi=_mm_mul_ps(a2hi, g2);
	a3lo=_mm_mul_ps(a3lo, g3), a3hi=_mm_mul_ps(a3hi, g3);
	a4lo=_mm_mul_ps(a4lo, g4), a4hi=_mm_mul_ps(a4hi, g4);
	a5lo=_mm_mul_ps(a5lo, g5), a5hi=_mm_mul_ps(a5hi, g5);
	a6lo=_mm_mul_ps(a6lo, g6), a6hi=_mm_mul_ps(a6hi, g6);
	a7lo=_mm_mul_ps(a7lo, g7), a7hi=_mm_mul_ps(a7hi, g7);
#endif
	
	_mm_store_ps(ptr   , a0lo), _mm_load_ps(ptr+ 4, a0hi);
	_mm_store_ps(ptr+ 8, a4lo), _mm_load_ps(ptr+12, a4hi);
	_mm_store_ps(ptr+16, a6lo), _mm_load_ps(ptr+20, a6hi);
	_mm_store_ps(ptr+24, a2lo), _mm_load_ps(ptr+28, a2hi);
	_mm_store_ps(ptr+32, a7lo), _mm_load_ps(ptr+36, a7hi);
	_mm_store_ps(ptr+40, a3lo), _mm_load_ps(ptr+44, a3hi);
	_mm_store_ps(ptr+48, a5lo), _mm_load_ps(ptr+52, a5hi);
	_mm_store_ps(ptr+56, a1lo), _mm_load_ps(ptr+60, a1hi);
#endif
}
void DCT3_8x8_ps(float *ptr)
{
	//INVERSE OF Loeffler's factorization using lifting		http://www.sfu.ca/~jiel/papers/c003-bindct-ieee.pdf
#ifdef __AVX__
	__m256
		half=_mm256_set1_ps(0.5f),
		p1=_mm256_set1_ps(13.f/32),//binDCT-L5
		u1=_mm256_set1_ps(11.f/32),
		p2=_mm256_set1_ps(19.f/64),
		u2=_mm256_set1_ps( 9.f/16),
		p3=_mm256_set1_ps(19.f/64),
		p4=_mm256_set1_ps( 3.f/32),
		u3=_mm256_set1_ps( 3.f/16),
		p5=_mm256_set1_ps( 3.f/32);
	
	__m256
		x0=_mm256_load_ps(ptr   ),
		x4=_mm256_load_ps(ptr+ 8),
		x6=_mm256_load_ps(ptr+16),
		x2=_mm256_load_ps(ptr+24),
		x7=_mm256_load_ps(ptr+32),
		x3=_mm256_load_ps(ptr+40),
		x5=_mm256_load_ps(ptr+48),
		x1=_mm256_load_ps(ptr+56);
	
	//new stage after each extra newline, to undo lifting:
	//- permuted load instead of store
	//- reverse operation order
	//- replace each operation with its inverse
	x4=_mm256_fmsub_ps(x7, half, x4);

	x7=_mm256_sub_ps(x7, x4);//add

	x4=_mm256_fmadd_ps(x6, half, x4);//_mm256_fnmadd_ps
	x7=_mm256_fmadd_ps(x5, half, x7);//_mm256_fnmadd_ps

	x1=_mm256_fmsub_ps(x0, half, x1);
	x3=_mm256_fmadd_ps(x2, u1, x3);//_mm256_fnmadd_ps
	x6=_mm256_sub_ps(x4, x6);
	x5=_mm256_sub_ps(x7, x5);

	x0=_mm256_sub_ps(x0, x1);//add
	x2=_mm256_fmsub_ps(x3, p1, x2);
	x7=_mm256_fmadd_ps(x4, p3, x7);//_mm256_fnmadd_ps
	x6=_mm256_fmadd_ps(x5, p5, x6);//_mm256_fnmadd_ps
	
	x0=_mm256_fmadd_ps(x3, half, x0);//_mm256_fnmadd_ps
	x1=_mm256_fmadd_ps(x2, half, x1);//_mm256_fnmadd_ps
	x4=_mm256_fnmadd_ps(x7, u2, x4);//_mm256_fmadd_ps
	x5=_mm256_fnmadd_ps(x6, u3, x5);//_mm256_fmadd_ps

	x3=_mm256_sub_ps(x0, x3);
	x2=_mm256_sub_ps(x1, x2);
	x7=_mm256_fmadd_ps(x4, p2, x7);//_mm256_fnmadd_ps
	x6=_mm256_fmadd_ps(x5, p4, x6);//_mm256_fnmadd_ps

	x0=_mm256_fmadd_ps(x7, half, x0);//_mm256_fnmadd_ps
	x1=_mm256_fmadd_ps(x6, half, x1);//_mm256_fnmadd_ps
	x2=_mm256_fmadd_ps(x5, half, x2);//_mm256_fnmadd_ps
	x3=_mm256_fmadd_ps(x4, half, x3);//_mm256_fnmadd_ps

	x7=_mm256_sub_ps(x0, x7);
	x6=_mm256_sub_ps(x1, x6);
	x5=_mm256_sub_ps(x2, x5);
	x4=_mm256_sub_ps(x3, x4);
	
	_mm256_store_ps(ptr   , x0);
	_mm256_store_ps(ptr+ 8, x1);
	_mm256_store_ps(ptr+16, x2);
	_mm256_store_ps(ptr+24, x3);
	_mm256_store_ps(ptr+32, x4);
	_mm256_store_ps(ptr+40, x5);
	_mm256_store_ps(ptr+48, x6);
	_mm256_store_ps(ptr+56, x7);
#else
#error TODO
#endif
}
#if 1
void DCTtest2()
{
	ALIGN(32) float buf[64], b2[64];
	
	float
		srcmin=0, srcmax=0,//src [0.000000 1.000000]		Y/Co/Cg	src [-1.000000 1.000061]
		dstmin=0, dstmax=0;//DCT [-4.664061 3.825694]				DCT [-7.948018 10.167978]
	for(int it=0;;++it)
	{
		int k;
		srand((unsigned)__rdtsc());

		for(k=0;k<64;++k)
			//buf[k]=roundf((float)rand()*(10000.f/RAND_MAX))*0.0001f;//[0, 1] 4 digits
			//buf[k]=(float)rand()/RAND_MAX;//[0, 1]
			buf[k]=(float)(rand()-(RAND_MAX>>1))/(RAND_MAX>>1);//[-1, 1]
		memcpy(b2, buf, 64*sizeof(float));
		for(k=0;k<64;++k)
		{
			if(!it&&!k||srcmin>b2[k])
				srcmin=b2[k];
			if(!it&&!k||srcmax<b2[k])
				srcmax=b2[k];
		}

		DCT2_8x8_ps(b2);
		transpose_block8x8_ps(b2);
		DCT2_8x8_ps(b2);
		
		for(k=0;k<64;++k)
		{
			if(!it&&!k||dstmin>b2[k])
				dstmin=b2[k];
			if(!it&&!k||dstmax<b2[k])
				dstmax=b2[k];
		}

		DCT3_8x8_ps(b2);
		transpose_block8x8_ps(b2);
		DCT3_8x8_ps(b2);

		printf("src [%f %f] DCT [%f %f]\t\t", srcmin, srcmax, dstmin, dstmax);
		//int nplaces=4;
		//print_fbuf(buf, 8, 8, nplaces);
		//print_fbuf(b2, 8, 8, nplaces);
		compare_bufs_ps(b2, buf, 8, 8, "DCT", 0);
		pause1();
	}

	printf("Done.\n");
	pause();
	exit(0);
}
#endif


//FCT-II/III, O(N*lg(N)), POT size (uses FFT size N)
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
		float rewm=(float)cosf(angle), imwm=(float)sinf(angle);//exp(-+ 2*pi*i*(N/m)/N)
		for(int k=0;k<size;k+=m)
		{

			float
				re0=re[k],
				im0=im[k],
				re1=re[k+m2],
				im1=im[k+m2];

			re[k]=re0+re1;
			im[k]=im0+im1;
			re[k+m2]=re0-re1;
			im[k+m2]=im0-im1;

			float rew=rewm, imw=imwm;
			for(int j=1;;)
			{
				int idx=k+j;
				re0=re[idx];
				im0=im[idx];
				re1=re[idx+m2];
				im1=im[idx+m2];
				
				float
					re2=rew*re1-imw*im1,
					im2=rew*im1+imw*re1;

				re[idx]=re0+re2;
				im[idx]=im0+im2;
				re[idx+m2]=re0-re2;
				im[idx+m2]=im0-im2;

				//nmuls+=4;//

				++j;
				if(j>=m2)
					break;

				re0=rewm*rew-imwm*imw;//w*=wm
				im0=rewm*imw+imwm*rew;
				rew=re0;
				imw=im0;

				//nmuls+=4;//
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
		float rewm=(float)cosf(angle), imwm=(float)sinf(angle);//exp(-+ 2*pi*i*(N/m)/N)
		for(int k=0;k<size;k+=m)
		{
			float
				re0=re[k],
				im0=im[k],
				re1=re[k+m2],
				im1=im[k+m2];

			re[k]=re0+re1;
			im[k]=im0+im1;
			re[k+m2]=re0-re1;
			im[k+m2]=im0-im1;
			
			float rew=rewm, imw=imwm;
			for(int j=1;;)
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

				//nmuls+=4;//
				
				++j;
				if(j>=m2)
					break;
				
				re0=rewm*rew-imwm*imw;//w*=wm
				im0=rewm*imw+imwm*rew;
				rew=re0;
				imw=im0;

				//nmuls+=4;//
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
#if 0
typedef struct FCT1D_PS_Params
{
	int lgsize;
	int *fctp;
	float *re, *im, *rew4N, *imw4N;
} FCT1D_PS_Params;
#endif
#define BIT_REVERSE(N, NBITS)\
	N=N<<1&0xAAAA|N>>1&0x5555;\
	N=N<<2&0xCCCC|N>>2&0x3333;\
	N=N<<4&0xF0F0|N>>4&0x0F0F;\
	N=N<<8&0xFF00|N>>8&0x00FF;\
	N>>=16-NBITS;
void FCT1D_ps_gen(int lgsize, FCT1D_PS_Params *p)//up to size 65536
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
		float temp=(float)(-2*M_PI*k/((size_t)size<<2));
		p->rew4N[k]=invsqrtN*cosf(temp);
		p->imw4N[k]=invsqrtN*sinf(temp);
	}
}
void FCT1D_ps_free(FCT1D_PS_Params *p)
{
	free(p->fctp);
	free(p->re);
	free(p->im);
	free(p->rew4N);
	free(p->imw4N);
	memset(p, 0, sizeof(*p));
}
void FCT1D_ps_fwd(FCT1D_PS_Params *p, float *data, int stride)//data, re, im, rew, imw are all of size (1<<lgsize), result is in data buffer
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
void FCT1D_ps_inv(FCT1D_PS_Params *p, float *data, int stride)//data, re, im, rew, imw are all of size (1<<lgsize), result is in data buffer
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


void printmat_pd(double *matrix, int bw, int bh, int prec)
{
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
			printf("%*.*lf", prec+5, prec, matrix[bw*ky+kx]);
		printf("\n");
	}
	printf("\n");
}
void matmul_pd(double *res, const double *A, const double *B, int h1, int w1h2, int w2)
{
	for(int ky=0;ky<h1;++ky)
	{
		for(int kx=0;kx<w2;++kx)
		{
			double sum=0;
			for(int ki=0;ki<w1h2;++ki)
				sum+=A[w1h2*ky+ki]*B[w2*ki+kx];
			res[w2*ky+kx]=sum;
		}
	}
}
void calc_DCT_matrix(int lgn, double *DCT2, double *DCT3)//arrays size (1<<(lgn<<1))
{
	int n=1<<lgn;
	double gain=sqrt(2./n), pi_n=M_PI/n;
	for(int ko=0;ko<n;++ko)
	{
		for(int ki=0;ki<n;++ki)
		{
			DCT2[ko<<2|ki]=gain*cos(pi_n*(ki+0.5)*ko);
			if(ki)
				DCT3[ko<<2|ki]=gain*cos(pi_n*(ko+0.5)*ki);
			else
				DCT3[ko<<2|ki]=gain*0.5;
		}
	}
}
#define lgN 1
#define N (1<<lgN)
double DCT2[N*N], DCT3[N*N], prod[N*N];
void test4()
{
	calc_DCT_matrix(lgN, DCT2, DCT3);

	printf("DCT2 (fwd):\n");
	printmat_pd(DCT2, N, N, 10);

	printf("DCT3 (inv):\n");
	printmat_pd(DCT3, N, N, 10);

	matmul_pd(DCT3, DCT2, prod, N, N, N);

	printf("Product:\n");
	printmat_pd(prod, N, N, 10);

	printf("Done.\n");
	pause();
	exit(0);
}