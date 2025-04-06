#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;

#define BLOCKX 16
#define BLOCKY 16

#define OLS6_REACH 1
#define OLS6_NPARAMS0 (2*(OLS6_REACH+1)*OLS6_REACH*3+0+1)
#define OLS6_NPARAMS1 (2*(OLS6_REACH+1)*OLS6_REACH*3+1+1)
#define OLS6_NPARAMS2 (2*(OLS6_REACH+1)*OLS6_REACH*3+2+1)
#define OLS6_NPARAMST (2*(OLS6_REACH+1)*OLS6_REACH*3*3+3+3)

#define PRECBITS 12


//	#define ENABLE_MA


typedef enum _RCTInfoIdx
{
	II_OCH_Y,
	II_OCH_U,
	II_OCH_V,

	II_PERM_Y,
	II_PERM_U,
	II_PERM_V,

	II_COEFF_U_SUB_Y,
	II_COEFF_V_SUB_Y,
	II_COEFF_V_SUB_U,

	II_COUNT,
} RCTInfoIdx;
typedef enum _OCHIndex
{
	OCH_R,
	OCH_G,
	OCH_B,
	OCH_RG,
	OCH_GB,
	OCH_BR,
	OCH_COUNT,
	OCH_GR=OCH_RG,
	OCH_BG=OCH_GB,
	OCH_RB=OCH_BR,
} OCHIndex;
#define RCTLIST\
	RCT(_400_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_400_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4)\
	RCT(_400_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0)\
	RCT(_040_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0)\
	RCT(_040_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4)\
	RCT(_004_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4)\
	RCT(_004_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0)\
	RCT(_040_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0)\
	RCT(_040_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4)\
	RCT(_040_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4)\
	RCT(_004_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0)\
	RCT(_004_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4)\
	RCT(_004_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4)\
	RCT(_400_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0)\
	RCT(_400_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4)\
	RCT(_400_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)
static const char *rct_names[]=
{
#define RCT(LABEL, ...) #LABEL,
	RCTLIST
#undef  RCT
};
typedef enum _RCTIndex
{
#define RCT(LABEL, ...) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][II_COUNT]=
{
#define RCT(LABEL, ...) {__VA_ARGS__},
	RCTLIST
#undef  RCT
};
static int ch_decorrelate(unsigned char *image, int iw, int ih)
{
	long long counters[6]={0};
	int W[6]={0};
	for(ptrdiff_t k=0, len=(ptrdiff_t)3*iw*ih;k<len;k+=3)
	{
		int
			r=image[k+0],
			g=image[k+1],
			b=image[k+2],
			rg=r-g,
			gb=g-b,
			br=b-r;
		counters[0]+=abs(r -W[0]);
		counters[1]+=abs(g -W[1]);
		counters[2]+=abs(b -W[2]);
		counters[3]+=abs(rg-W[3]);
		counters[4]+=abs(gb-W[4]);
		counters[5]+=abs(br-W[5]);
		W[0]=r;
		W[1]=g;
		W[2]=b;
		W[3]=rg;
		W[4]=gb;
		W[5]=br;
	}
	int bestrct=0;
	long long minerr=0;
	for(int kt=0;kt<RCT_COUNT;++kt)
	{
		const unsigned char *rct=rct_combinations[kt];
		long long currerr=
			+counters[rct[0]]
			+counters[rct[1]]
			+counters[rct[2]]
		;
		if(!kt||minerr>currerr)
		{
			minerr=currerr;
			bestrct=kt;
		}
	}
	const unsigned char *combination=rct_combinations[bestrct];
	int
		yidx=combination[II_PERM_Y],
		uidx=combination[II_PERM_U],
		vidx=combination[II_PERM_V];
	for(ptrdiff_t k=0, len=(ptrdiff_t)3*iw*ih;k<len;k+=3)
	{
		char yuv[]=
		{
			image[k+yidx]-128,
			image[k+uidx]-128,
			image[k+vidx]-128,
		};
#ifdef ENABLE_MA
		yuv[2]-=(combination[II_COEFF_V_SUB_Y]*yuv[0]+combination[II_COEFF_V_SUB_U]*yuv[1])>>2;
		yuv[1]-=yuv[0]&-(combination[II_COEFF_U_SUB_Y]!=0);
#endif
		image[k+0]=yuv[0]+128;
		image[k+1]=yuv[1]+128;
		image[k+2]=yuv[2]+128;
	}
	return bestrct;
}
static void ch_reconstruct(unsigned char *image, int iw, int ih, int rct)
{
	const unsigned char *combination=rct_combinations[rct];
	int
		yidx=combination[II_PERM_Y],
		uidx=combination[II_PERM_U],
		vidx=combination[II_PERM_V];
	for(ptrdiff_t k=0, len=(ptrdiff_t)3*iw*ih;k<len;k+=3)
	{
		char yuv[]=
		{
			image[k+0]-128,
			image[k+1]-128,
			image[k+2]-128,
		};
#ifdef ENABLE_MA
		yuv[1]+=yuv[0]&-(combination[II_COEFF_U_SUB_Y]!=0);
		yuv[2]+=(combination[II_COEFF_V_SUB_Y]*yuv[0]+combination[II_COEFF_V_SUB_U]*yuv[1])>>2;
#endif
		image[k+yidx]=yuv[0]+128;
		image[k+uidx]=yuv[1]+128;
		image[k+vidx]=yuv[2]+128;
	}
}

static double calc_entropy(const int *hist)
{
	double e=0;
	int sum=0;
	for(int k=0;k<256;++k)
		sum+=hist[k];
	for(int k=0;k<256;++k)
	{
		int freq=hist[k];
		if(freq)
			e-=freq*log2((double)freq/sum);
	}
	//printf(" %d", sum);//
	return e/8;
}
static int solve_Mx_v_cholesky(double *matrix, const double *vec, int n, double *solution)
{
	int success=1;
	double sum;
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<i;++j)
		{
			sum=matrix[i*n+j];
			for(int k=0;k<j;++k)
				sum-=matrix[i*n+k]*matrix[j*n+k];
			matrix[i*n+j]=sum/matrix[j*n+j];
		}
		sum=matrix[i*n+i];
		for(int k=0;k<i;++k)
			sum-=matrix[i*n+k]*matrix[i*n+k];
		if(sum<=1e-8)
		{
			success=0;
			break;
		}
		matrix[i*n+i]=sqrt(sum);
	}
	if(success)
	{
		for(int ky=0;ky<n;++ky)
		{
			//solution[0] = vec[0]/matrix[0][0]
			//solution[1] = (vec[1] - solution[0]*matrix[1][0])/matrix[1][1]
			//solution[2] = (vec[2] - solution[0]*matrix[2][0] - solution[1]*matrix[2][1])/matrix[2][2]
			//...
			sum=vec[ky];
			for(int kx=0;kx<ky;++kx)//lower triangle downwards	_*\ |
				sum-=matrix[ky*n+kx]*solution[kx];
			solution[ky]=sum/matrix[(n+1)*ky];
		}
		for(int kx=n-1;kx>=0;--kx)
		{
			sum=solution[kx];
			for(int ky=kx+1;ky<n;++ky)//upper triangle upwards	_ \*|
				sum-=matrix[ky*n+kx]*solution[ky];
			solution[kx]=sum/matrix[(n+1)*kx];
		}
	}
	return success;
}
static int solveblock(const unsigned char *image, int iw, int ih, int x1, int x2, int y1, int y2, short *coeffs)//13+14+15 coeffs
{
	//const int msize0=OLS6_NPARAMS0*OLS6_NPARAMS0;
	//const int msize1=OLS6_NPARAMS1*OLS6_NPARAMS1;
	//const int msize2=OLS6_NPARAMS2*OLS6_NPARAMS2;
	double matrix0[OLS6_NPARAMS0*OLS6_NPARAMS0]={0}, vec0[OLS6_NPARAMS0]={0};
	double matrix1[OLS6_NPARAMS1*OLS6_NPARAMS1]={0}, vec1[OLS6_NPARAMS1]={0};
	double matrix2[OLS6_NPARAMS2*OLS6_NPARAMS2]={0}, vec2[OLS6_NPARAMS2]={0};
	double solution[OLS6_NPARAMST]={0};
	int nb[16]={0};
	const unsigned char *ptr=image+3*(iw*y1+x1);
	int rowstride=3*iw;
	for(int ky2=y1;ky2<y2;++ky2)
	{
		const unsigned char *ptr2=ptr;
		memset(nb, 0, sizeof(nb));
		for(int kx2=x1;kx2<x2;++kx2, ptr2+=3)
		{
			nb[0]=128;
			if(ky2)
			{
				if(kx2)
				{
					nb[1]=ptr2[0-1*3-rowstride]-128;
					nb[2]=ptr2[1-1*3-rowstride]-128;
					nb[3]=ptr2[2-1*3-rowstride]-128;
				}
				nb[4]=ptr2[0+0*3-rowstride]-128;
				nb[5]=ptr2[1+0*3-rowstride]-128;
				nb[6]=ptr2[2+0*3-rowstride]-128;
				if(kx2<iw-1)
				{
					nb[7]=ptr2[0+1*3-rowstride]-128;
					nb[8]=ptr2[1+1*3-rowstride]-128;
					nb[9]=ptr2[2+1*3-rowstride]-128;
				}
				else
				{
					nb[7]=0;
					nb[8]=0;
					nb[9]=0;
				}
			}
			if(kx2)
			{
				nb[10]=ptr2[0-1*3]-128;
				nb[11]=ptr2[1-1*3]-128;
				nb[12]=ptr2[2-1*3]-128;
			}
			nb[13]=ptr2[0+0*3]-128;
			nb[14]=ptr2[1+0*3]-128;
			nb[15]=ptr2[2+0*3]-128;
			for(int ky3=0, idx3=0;ky3<OLS6_NPARAMS0;++ky3)
			{
				double val=nb[ky3];
				for(int kx3=0;kx3<OLS6_NPARAMS0;++kx3, ++idx3)
					matrix0[idx3]+=val*nb[kx3];
				vec0[ky3]+=val*nb[OLS6_NPARAMS0];
			}
			for(int ky3=0, idx3=0;ky3<OLS6_NPARAMS1;++ky3)
			{
				double val=nb[ky3];
				for(int kx3=0;kx3<OLS6_NPARAMS1;++kx3, ++idx3)
					matrix1[idx3]+=val*nb[kx3];
				vec1[ky3]+=val*nb[OLS6_NPARAMS1];
			}
			for(int ky3=0, idx3=0;ky3<OLS6_NPARAMS2;++ky3)
			{
				double val=nb[ky3];
				for(int kx3=0;kx3<OLS6_NPARAMS2;++kx3, ++idx3)
					matrix2[idx3]+=val*nb[kx3];
				vec2[ky3]+=val*nb[OLS6_NPARAMS2];
			}
		}
		ptr+=rowstride;
	}
	int success=1;
	success&=solve_Mx_v_cholesky(matrix0, vec0, OLS6_NPARAMS0, solution);
	success&=solve_Mx_v_cholesky(matrix1, vec1, OLS6_NPARAMS1, solution+OLS6_NPARAMS0);
	success&=solve_Mx_v_cholesky(matrix2, vec2, OLS6_NPARAMS2, solution+OLS6_NPARAMS0+OLS6_NPARAMS1);
	for(int k=0;k<OLS6_NPARAMS0+OLS6_NPARAMS1+OLS6_NPARAMS2;++k)
	{
		double val=solution[k];
		val*=1<<PRECBITS;
		if(val<-0x8000||val>0x7FFF)
			success=0;
		coeffs[k]=(short)CVTFP64_I64(val);
	}
	return success;
}
static void ols_predict(unsigned char *image, int iw, int ih, const short *allcoeffs, int fwd)
{
	int psize=(iw+16LL)*(int)sizeof(short[4*4]);//4 padded rows * 4 channels max * pixels
	short *pixels=(short*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	int xblocks=iw/BLOCKX;//floor
	unsigned char *ptr=image;
	for(int ky=0;ky<ih;++ky)
	{
		ALIGN(32) short *rows[]=
		{
			pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4,
		};
		int by=ky/BLOCKY;
		if((by+1)*BLOCKY>ih)
			--by;
		const short *coeffrow=allcoeffs+OLS6_NPARAMST*xblocks*by;
		int kc=0, pred=0, curr=0;
		int vmin=0, vmax=0;
		for(int kx=0;kx<iw;++kx)
		{
			int bx=kx/BLOCKX;
			if((bx+1)*BLOCKX>iw)
				--bx;
			const short *coeffs=coeffrow+OLS6_NPARAMST*bx;
			short nb[]=
			{
				128,
				rows[1][0-1*4],//NW	1~3
				rows[1][1-1*4],
				rows[1][2-1*4],
				rows[1][0+0*4],//N	4~6
				rows[1][1+0*4],
				rows[1][2+0*4],
				rows[1][0+1*4],//NE	7~9
				rows[1][1+1*4],
				rows[1][2+1*4],
				rows[0][0-1*4],//W	10~12
				rows[0][1-1*4],
				rows[0][2-1*4],
				0,//curr	13~15
				0,
				0,
			};

			//predict Y
			kc=0;
			pred=0;
			for(int k=0;k<OLS6_NPARAMS0;++k)
				pred+=*coeffs++*nb[k];
			pred=(pred+(1<<PRECBITS>>1))>>PRECBITS;
			vmax=nb[4+kc]; vmin=nb[10+kc];
			if(nb[4+kc]<nb[10+kc])vmin=nb[4+kc], vmax=nb[10+kc];
			if(vmin>nb[7+kc])vmin=nb[7+kc];
			if(vmax<nb[7+kc])vmax=nb[7+kc];
			CLAMP2(pred, vmin, vmax);
			if(fwd)
			{
				curr=*ptr-128;
				*ptr=(char)(curr-pred)+128;
			}
			else
			{
				curr=(char)(*ptr+pred-128);
				*ptr=curr+128;
			}
			rows[0][kc]=nb[13+kc]=curr;
			++ptr;

			//predict U
			kc=1;
			pred=0;
			for(int k=0;k<OLS6_NPARAMS1;++k)
				pred+=*coeffs++*nb[k];
			pred=(pred+(1<<PRECBITS>>1))>>PRECBITS;
			vmax=nb[4+kc]; vmin=nb[10+kc];
			if(nb[4+kc]<nb[10+kc])vmin=nb[4+kc], vmax=nb[10+kc];
			if(vmin>nb[7+kc])vmin=nb[7+kc];
			if(vmax<nb[7+kc])vmax=nb[7+kc];
			CLAMP2(pred, vmin, vmax);
			if(fwd)
			{
				curr=*ptr-128;
				*ptr=(char)(curr-pred)+128;
			}
			else
			{
				curr=(char)(*ptr+pred-128);
				*ptr=curr+128;
			}
			rows[0][kc]=nb[13+kc]=curr;
			++ptr;

			//predict V
			kc=2;
			pred=0;
			for(int k=0;k<OLS6_NPARAMS2;++k)
				pred+=*coeffs++*nb[k];
			pred=(pred+(1<<PRECBITS>>1))>>PRECBITS;
			vmax=nb[4+kc]; vmin=nb[10+kc];
			if(nb[4+kc]<nb[10+kc])vmin=nb[4+kc], vmax=nb[10+kc];
			if(vmin>nb[7+kc])vmin=nb[7+kc];
			if(vmax<nb[7+kc])vmax=nb[7+kc];
			CLAMP2(pred, vmin, vmax);
			if(fwd)
			{
				curr=*ptr-128;
				*ptr=(char)(curr-pred)+128;
			}
			else
			{
				curr=(char)(*ptr+pred-128);
				*ptr=curr+128;
			}
			rows[0][kc]=nb[13+kc]=curr;
			++ptr;

			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	_mm_free(pixels);
}

#define NCTX 17
#define NPREDS 8
#define PREDLIST\
	PRED(240, N)\
	PRED(240, W)\
	PRED(180, 3*(N-NN)+NNN)\
	PRED(180, 3*(W-WW)+WWW)\
	PRED(140, W+NE-N)\
	PRED(160, (WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)\
	PRED(120, N+W-NW)\
	PRED(120, N+NE-NNE)

AWM_INLINE int wp_mix_naive(const int *preds, const int *errors, const int *iweights)//faster
{
	long long pred=0, wsum=0;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
	for(int k=0;k<NPREDS;++k)
	{
		int coeff=(unsigned)(0x800000*iweights[k])/(unsigned)(errors[k]+1);
		pred+=(long long)preds[k]*coeff;
		wsum+=coeff;
	}
//	pred+=wsum>>1;
	pred/=wsum;
	return (int)pred;
}
AWM_INLINE int wp_mix_eBench(const int *preds, const int *errors, const float *weights)
{
	__m256i me=_mm256_add_epi32(_mm256_load_si256((__m256i*)errors), _mm256_set1_epi32(1));
	__m256 fp=_mm256_cvtepi32_ps(_mm256_load_si256((__m256i*)preds));
	__m256 fe=_mm256_cvtepi32_ps(me);
	fe=_mm256_div_ps(_mm256_load_ps(weights), fe);

	fp=_mm256_mul_ps(fp, fe);
	fe=_mm256_hadd_ps(fe, fe);
	fp=_mm256_hadd_ps(fp, fp);
	fe=_mm256_hadd_ps(fe, fe);
	fp=_mm256_hadd_ps(fp, fp);
	__m128 fe4=_mm_add_ps(_mm256_castps256_ps128(fe), _mm256_extractf128_ps(fe, 1));
	__m128 fp4=_mm_add_ps(_mm256_castps256_ps128(fp), _mm256_extractf128_ps(fp, 1));
	fp4=_mm_div_ss(fp4, fe4);
	int pred=_mm_cvt_ss2si(fp4);
	return pred;
}
AWM_INLINE int wp_mix(const int *preds, const int *errors, const float *weights)
{
	__m256i mp=_mm256_load_si256((__m256i*)preds);
	__m256i me=_mm256_load_si256((__m256i*)errors);
	me=_mm256_srli_epi32(me, 1);
	me=_mm256_add_epi32(me, _mm256_set1_epi32(1));
	__m256 fe=_mm256_cvtepi32_ps(me);
	__m256 fp=_mm256_cvtepi32_ps(mp);
	fe=_mm256_rcp_ps(fe);
	fe=_mm256_mul_ps(fe, _mm256_load_ps(weights));
	fp=_mm256_mul_ps(fp, fe);
	fe=_mm256_hadd_ps(fe, fe);
	fp=_mm256_hadd_ps(fp, fp);
	fe=_mm256_hadd_ps(fe, fe);
	fp=_mm256_hadd_ps(fp, fp);
	__m128 fe4=_mm_add_ps(_mm256_castps256_ps128(fe), _mm256_extractf128_ps(fe, 1));
	__m128 fp4=_mm_add_ps(_mm256_castps256_ps128(fp), _mm256_extractf128_ps(fp, 1));
	fp4=_mm_div_ss(fp4, fe4);
	int pred=_mm_cvt_ss2si(fp4);
	return pred;
}
static void wp_estimate(const unsigned char *image, int iw, int ih, int rct, int usewp, int blockx, int blocky, double *esizes, double *elapsed)
{
	double t=time_sec();
#ifndef ENABLE_MA
	const unsigned char *combination=rct_combinations[rct];
	//int
	//	yidx=combination[II_PERM_Y],
	//	uidx=combination[II_PERM_U],
	//	vidx=combination[II_PERM_V];
	int vfromy=-(combination[II_COEFF_U_SUB_Y]!=0);
#endif
	ALIGN(32) static const float weights[]=
	{
#define PRED(WEIGHT, EXPR) WEIGHT,
		PREDLIST
#undef  PRED
	};
	ALIGN(32) static const int iweights[]=
	{
#define PRED(WEIGHT, EXPR) WEIGHT,
		PREDLIST
#undef  PRED
	};
	int psize=(int)sizeof(short[4*3*(2+NPREDS)])*(2*blockx+16);//4 padded rows * 3 channels * {yuv, ctx, prederrors...}
	short *pixels=(short*)malloc(psize);
	int hsize=(int)sizeof(int[3][NCTX][256]);
	int *hist=(int*)malloc(hsize);
	//int hist0[3][256]={0};
	if(!pixels||!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(hist, 0, hsize);
	for(int y1=0;y1<ih;y1+=blocky)
	{
		int y2=y1+blocky;
		if(y2+blocky>ih)
			y2=ih;
		for(int x1=0;x1<iw;x1+=blockx)
		{
			int x2=x1+blockx;
			if(x2+blockx>iw)
				x2=iw;
			ALIGN(32) int errors[NPREDS]={0};
			memset(pixels, 0, psize);
			for(int ky=y1;ky<y2;++ky)
			{
				ALIGN(32) short *rows[]=
				{
					pixels+((x2-x1+16LL)*((ky-0LL)&3)+8LL)*3*(2+NPREDS),
					pixels+((x2-x1+16LL)*((ky-1LL)&3)+8LL)*3*(2+NPREDS),
					pixels+((x2-x1+16LL)*((ky-2LL)&3)+8LL)*3*(2+NPREDS),
					pixels+((x2-x1+16LL)*((ky-3LL)&3)+8LL)*3*(2+NPREDS),
				};
				const unsigned char *ptr=image+3*(iw*ky+x1);
				ALIGN(32) int preds[NPREDS]={0}, e2[NPREDS]={0};
				int pred=0, vmin=0, vmax=0;
				int p2=0;
				(void)p2;
				for(int kx=x1;kx<x2;++kx)
				{
#ifndef ENABLE_MA
					int y=0, offset=0;
#endif
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 3
#endif
					for(int kc=0;kc<3;++kc)
					{
						short
							NNN	=rows[3][0+0*3*(2+NPREDS)],
							NN	=rows[2][0+0*3*(2+NPREDS)],
							NNE	=rows[2][0+1*3*(2+NPREDS)],
							NW	=rows[1][0-1*3*(2+NPREDS)],
							N	=rows[1][0+0*3*(2+NPREDS)],
							NE	=rows[1][0+1*3*(2+NPREDS)],
							NEE	=rows[1][0+2*3*(2+NPREDS)],
							NEEE	=rows[1][0+3*3*(2+NPREDS)],
							NEEEE	=rows[1][0+4*3*(2+NPREDS)],
							WWWW	=rows[0][0-4*3*(2+NPREDS)],
							WWW	=rows[0][0-3*3*(2+NPREDS)],
							WW	=rows[0][0-2*3*(2+NPREDS)],
							W	=rows[0][0-1*3*(2+NPREDS)],
							eW	=rows[0][1-1*3*(2+NPREDS)];
						int ctx=FLOOR_LOG2(eW*eW+1);
						if(ctx>NCTX-1)
							ctx=NCTX-1;

						vmax=N, vmin=W;
						if(N<W)vmin=N, vmax=W;
						if(!usewp)
							pred=N+W-NW;
						else
						{
							int j=0;
#define PRED(WEIGHT, EXPR) preds[j++]=EXPR;
							PREDLIST
#undef  PRED
							switch(usewp)
							{
							case 1:
								pred=wp_mix(preds, errors, weights);
								//{
								//	p2=wp_mix_naive(preds, errors, iweights);
								//	if(abs(p2-pred)>5)
								//		printf("");
								//}
								break;
							case 2:
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
								for(int k=0;k<NPREDS;++k)
									e2[k]=errors[k]
										+rows[1][k+2-1*3*(2+NPREDS)]
										+2*(rows[1][k+2+0*3*(2+NPREDS)]+rows[0][k+2-1*3*(2+NPREDS)])
										+rows[1][k+2+1*3*(2+NPREDS)]
										+rows[2][k+2+1*3*(2+NPREDS)];
								pred=wp_mix(preds, e2, weights);
								break;
							case 3:
								pred=wp_mix_naive(preds, errors, iweights);
								break;
							case 4:
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
								for(int k=0;k<NPREDS;++k)
									e2[k]=errors[k]
										+rows[1][k+2-1*3*(2+NPREDS)]
										+2*(rows[1][k+2+0*3*(2+NPREDS)]+rows[0][k+2-1*3*(2+NPREDS)])
										+rows[1][k+2+1*3*(2+NPREDS)]
										+rows[2][k+2+1*3*(2+NPREDS)];
								pred=wp_mix_naive(preds, e2, iweights);
								break;
							case 5:
								pred=wp_mix_eBench(preds, errors, weights);
								break;
							case 6:
#if 1
								//pI+eWW+eNW+2*(eN+eW)+eNE+eNNE
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
								for(int k=0;k<NPREDS;++k)
									e2[k]=(errors[k]>>2)
										+rows[0][k+2-2*3*(2+NPREDS)]
										+rows[1][k+2-1*3*(2+NPREDS)]
										+3*(rows[1][k+2+0*3*(2+NPREDS)]+rows[0][k+2-1*3*(2+NPREDS)])
										+rows[1][k+2+1*3*(2+NPREDS)]
										+rows[2][k+2+0*3*(2+NPREDS)]
										+rows[2][k+2+1*3*(2+NPREDS)]
										+rows[1][k+2+2*3*(2+NPREDS)]
										+rows[1][k+2+3*3*(2+NPREDS)];
								pred=wp_mix_eBench(preds, e2, weights);
#endif
							//	for(int k=0;k<NPREDS;++k)
							//	{
							//		e2[k]=errors[k]+rows[1][k+2-1*3*(2+8)]+2*rows[1][k+2+0*3*(2+8)]+rows[1][k+2+1*3*(2+8)]+rows[2][k+2+1*3*(2+8)];
							//		//if(e2[k]<0)//
							//		//	LOG_ERROR("");
							//	}
							//	pred=(e2[1]*preds[0]+e2[0]*preds[1])/(e2[0]+e2[1]+1);
							//	pred=(preds[0]+preds[1]+preds[2]+preds[3]+preds[4]+preds[5]+preds[6]+preds[7]+4)>>3;//
							//	pred=preds[3];//
								break;
							case 7:
								{
									int p1=wp_mix(preds, errors, weights);		//3176982.10
									int p2=wp_mix_naive(preds, errors, iweights);	//3216902.64	worst
									int p3=wp_mix_eBench(preds, errors, weights);	//3176849.18
									pred=p3;
									(void)p1;
									(void)p2;
									(void)p3;
								}
								break;
							case 8:
								{
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
									for(int k=0;k<NPREDS;++k)
										e2[k]=errors[k]+rows[1][k+2-1*3*(2+NPREDS)]+2*rows[1][k+2+0*3*(2+NPREDS)]+rows[1][k+2+1*3*(2+NPREDS)]+rows[2][k+2+1*3*(2+NPREDS)];
									int p1=wp_mix(preds, e2, weights);		//3148157.01
									int p2=wp_mix_naive(preds, e2, iweights);	//3188495.31	worst
									int p3=wp_mix_eBench(preds, e2, weights);	//3148094.15
									pred=p3;
									(void)p1;
									(void)p2;
									(void)p3;
								}
								break;
							}
							//if(usewp==1)
							//	pred=wp_mix(preds, errors, weights);
							//else
							//{
							//	//if(ky==10&&kx==10)//
							//	//	printf("");
							//
							//	for(int k=0;k<NPREDS;++k)
							//		e2[k]=errors[k]+rows[1][k+2-1*3*(2+NPREDS)]+2*rows[1][k+2+0*3*(2+NPREDS)]+rows[1][k+2+1*3*(2+NPREDS)]+rows[2][k+2+1*3*(2+NPREDS)];
							//	pred=wp_mix(preds, e2, weights);
							//}
							if(vmin>NE)vmin=NE;
							if(vmax<NE)vmax=NE;
						}
						CLAMP2(pred, vmin, vmax);
#ifndef ENABLE_MA
						pred+=offset;
						CLAMP2(pred, -128, 127);
#endif

						int curr=*ptr++-128;
						int error=(char)(curr-pred);

						++hist[(NCTX*kc+ctx)<<8|(error+128)];
						
#ifdef ENABLE_MA
						rows[0][0]=curr;
#else
						rows[0][0]=curr-offset;
#endif
						rows[0][1]=(2*eW+((error<<1^error>>31)<<3)+rows[1][1+2*3*(2+NPREDS)])>>2;
						if(usewp)
						{
							//static const int factors[]={97, 99, 99};
							//int factor=factors[kc];
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
							for(int k=0;k<NPREDS;++k)
							{
#ifdef ENABLE_MA
								int e=abs(curr-preds[k])<<1;
#else
								int e=abs(curr-offset-preds[k])<<1;
#endif
							//	errors[k]=(errors[k]+e)*factor>>7;
								errors[k]+=((e<<2)-errors[k]+(1<<3>>1))>>3;

								//(2*eW+e+min(eN, eNEE))>>2
								rows[0][k+2+0*3*(2+NPREDS)]=(2*rows[0][k+2-1*3*(2+NPREDS)]+e+MINVAR(rows[1][k+2+0*3*(2+NPREDS)], rows[1][k+2+2*3*(2+NPREDS)]))>>2;
								
								//(eW+min(eN+eW)+e+min(eN, eNEE))>>2
								//rows[0][k+2+0*3*(2+NPREDS)]=(
								//	+rows[0][k+2-1*3*(2+NPREDS)]
								//	+MINVAR(rows[1][k+2+0*3*(2+NPREDS)], rows[0][k+2-1*3*(2+NPREDS)])
								//	+e
								//	+MINVAR(rows[1][k+2+0*3*(2+NPREDS)], rows[1][k+2+2*3*(2+NPREDS)])
								//)>>2;

								rows[1][k+2+1*3*(2+NPREDS)]=MAXVAR(rows[1][k+2+1*3*(2+NPREDS)], e);
							//	rows[1][k+2+1*3*(2+NPREDS)]=(rows[1][k+2+1*3*(2+NPREDS)]+3*e)>>2;
							//	if(rows[1][k+2+1*3*(2+8)]>256*16)//
							//		LOG_ERROR("");
							}
						}
#ifndef ENABLE_MA
						switch(kc)
						{
						case 0:
							offset=curr&vfromy;
							y=curr;
							break;
						case 1:
							offset=(combination[II_COEFF_V_SUB_Y]*y+combination[II_COEFF_V_SUB_U]*curr)>>2;
							break;
						}
#endif
						rows[0]+=2+NPREDS;
						rows[1]+=2+NPREDS;
						rows[2]+=2+NPREDS;
						rows[3]+=2+NPREDS;
					}
				}
			}
			if(x2==iw)
				break;
		}
		if(y2==ih)
			break;
	}
	free(pixels);
	memset(esizes, 0, sizeof(double[3]));
	for(int k=0;k<NCTX;++k)
	{
		esizes[0]+=calc_entropy(hist+((ptrdiff_t)0*NCTX+k)*256);
		esizes[1]+=calc_entropy(hist+((ptrdiff_t)1*NCTX+k)*256);
		esizes[2]+=calc_entropy(hist+((ptrdiff_t)2*NCTX+k)*256);
		//printf("\n");//
	}
	free(hist);
	t=time_sec()-t;
	if(elapsed)*elapsed=t;
}

static void calc_esize(const unsigned char *image, int iw, int ih, double *esizes0, double *esizes1)//esizes0[3], esizes1[3][NCTX]
{
	int rowstride=3*iw;
	const unsigned char *ptr=image;
	int hist0[3][256]={0};
	int hist1size=(int)sizeof(int[3][NCTX][256]);
	int *hist1=(int*)malloc(hist1size);
	if(!hist1)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(hist1, 0, hist1size);
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, ptr+=3)
		{
			int y=ptr[0];
			int u=ptr[1];
			int v=ptr[2];
			int yctx=0, uctx=0, vctx=0;
			if(ky)
			{
				if(kx)
				{
					yctx+=abs(ptr[0-1*3-rowstride]-128);
					uctx+=abs(ptr[1-1*3-rowstride]-128);
					vctx+=abs(ptr[2-1*3-rowstride]-128);
				}
				yctx+=abs(ptr[0+0*3-rowstride]-128);
				uctx+=abs(ptr[1+0*3-rowstride]-128);
				vctx+=abs(ptr[2+0*3-rowstride]-128);
				if(kx<iw-1)
				{
					yctx+=abs(ptr[0+1*3-rowstride]-128);
					uctx+=abs(ptr[1+1*3-rowstride]-128);
					vctx+=abs(ptr[2+1*3-rowstride]-128);
				}
			}
			if(kx)
			{
				yctx+=abs(ptr[0-1*3]-128);
				uctx+=abs(ptr[1-1*3]-128);
				vctx+=abs(ptr[2-1*3]-128);
			}
			yctx=FLOOR_LOG2(yctx*yctx+1);
			uctx=FLOOR_LOG2(uctx*uctx+1);
			vctx=FLOOR_LOG2(vctx*vctx+1);
			if(yctx>NCTX-1)yctx=NCTX-1;
			if(uctx>NCTX-1)uctx=NCTX-1;
			if(vctx>NCTX-1)vctx=NCTX-1;
			++hist0[0][y];
			++hist0[1][u];
			++hist0[2][v];
			++hist1[(0*NCTX+yctx)*256+y];
			++hist1[(1*NCTX+uctx)*256+u];
			++hist1[(2*NCTX+vctx)*256+v];
		}
	}
	esizes0[0]=calc_entropy(hist0[0]);
	esizes0[1]=calc_entropy(hist0[1]);
	esizes0[2]=calc_entropy(hist0[2]);
	for(int k=0;k<NCTX;++k)
	{
		esizes1[0*NCTX+k]=calc_entropy(hist1+((ptrdiff_t)0*NCTX+k)*256);
		esizes1[1*NCTX+k]=calc_entropy(hist1+((ptrdiff_t)1*NCTX+k)*256);
		esizes1[2*NCTX+k]=calc_entropy(hist1+((ptrdiff_t)2*NCTX+k)*256);
	}
	free(hist1);
}
static void print_yuvsizes(double *sizes, ptrdiff_t res, double elapsed)
{
	double ctotal=sizes[0]+sizes[1]+sizes[2];
	printf("U %9td  ", res*3);
	printf("TYUV %12.2lf %12.2lf %12.2lf %12.2lf"
		, ctotal
		, sizes[0]
		, sizes[1]
		, sizes[2]
	);
	printf(" TYUV %8.4lf%% %8.4lf%% %8.4lf%% %8.4lf%%  %12.6lf sec %12.6lf MB/s\n"
		, ctotal/(3*res)*100
		, sizes[0]/res*100
		, sizes[1]/res*100
		, sizes[2]/res*100
		, elapsed
		, 3*res/(elapsed*1024*1024)
	);
	//printf(
	//	"T %8.4lf%% %12.2lf / %12td\n"
	//	"Y %8.4lf%% %12.2lf / %12td\n"
	//	"U %8.4lf%% %12.2lf\n"
	//	"V %8.4lf%% %12.2lf\n"
	//	, ctotal/(3*res)*100, ctotal, res*3
	//	, sizes[0]/res*100, sizes[0], res
	//	, sizes[1]/res*100, sizes[1]
	//	, sizes[2]/res*100, sizes[2]
	//);
}

static unsigned char* load_ppm(const char *fn, int *ret_iw, int *ret_ih)
{
	FILE *fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		LOG_ERROR("Cannot open \"%s\"", fn);
		return 0;
	}
	int tag=0;
	fread(&tag, 1, 2, fsrc);
	if(tag!=('P'|'6'<<8))
	{
		LOG_ERROR("Unsupported file \"%s\"", fn);
		return 0;
	}
#ifdef LOUD
	print_timestamp("%Y-%m-%d_%H%M%S\n");
#endif
	int temp=fgetc(fsrc);
	if(temp!='\n')
	{
		LOG_ERROR("Invalid PPM file");
		return 0;
	}
	int iw=0, ih=0;
	int nread=fscanf(fsrc, "%d %d", &iw, &ih);
	if(nread!=2)
	{
		LOG_ERROR("Unsupported PPM file");
		return 0;
	}
	int vmax=0;
	nread=fscanf(fsrc, "%d", &vmax);
	if(nread!=1||vmax!=255)
	{
		LOG_ERROR("Unsupported PPM file");
		return 0;
	}
	temp=fgetc(fsrc);
	if(temp!='\n')
	{
		LOG_ERROR("Invalid PPM file");
		return 0;
	}
	if(iw<1||ih<1)
	{
		LOG_ERROR("Unsupported source file");
		return 0;
	}
	ptrdiff_t size=(ptrdiff_t)3*iw*ih;
	unsigned char *image=(unsigned char*)malloc(size+sizeof(__m256i));
	fread(image, 1, size, fsrc);//read image
	fclose(fsrc);
	if(ret_iw)*ret_iw=iw;
	if(ret_ih)*ret_ih=ih;
	return image;
}
static void save_pgm(const char *fn, const unsigned char *image, int iw, int ih)
{
	FILE *fdst=fopen(fn, "wb");
	if(!fdst)
	{
		LOG_ERROR("Cannot open \"%s\" for writing", fn);
		return;
	}
	fprintf(fdst, "P5\n%d %d\n255\n", iw, ih);
	fwrite(image, 1, (ptrdiff_t)iw*ih, fdst);
	fclose(fdst);
}
static void save_ppm(const char *fn, const unsigned char *image, int iw, int ih)
{
	FILE *fdst=fopen(fn, "wb");
	if(!fdst)
	{
		LOG_ERROR("Cannot open \"%s\" for writing", fn);
		return;
	}
	fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	fwrite(image, 1, (ptrdiff_t)3*iw*ih, fdst);
	fclose(fdst);
}
int c30_codec(const char *srcfn, const char *dstfn, int nthreads0)
{
	const char *ext[]=
	{
		"ppm",
	};
	ArrayHandle fns=get_filenames(srcfn, ext, _countof(ext), 1);
	if(!fns)
	{
		LOG_ERROR("No images in \"%s\"", srcfn);
		return 1;
	}
	double total=0;
	long long usize=0;
	double csize=0;
	for(int k=0;k<(int)fns->count;++k)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&fns, k);
		int iw=0, ih=0;
		unsigned char *image=load_ppm((char*)fn[0]->data, &iw, &ih);
		if(!image)
			continue;
		int rct=ch_decorrelate(image, iw, ih);

		ptrdiff_t res=(ptrdiff_t)iw*ih;
		double esizes[3]={0}, t=0;
		wp_estimate(image, iw, ih, rct, 6, iw/4, ih/8, esizes, &t);
		printf("%4d ", k+1);
		print_yuvsizes(esizes, res, t);

		usize+=3LL*iw*ih;
		csize+=esizes[0]+esizes[1]+esizes[2];
		total+=t;

		free(image);
	}
	printf("%16.2lf/%16lld  %12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB\n"
		, csize
		, usize
		, total
		, usize/(total*1024*1024)
		, (total*1000*1024*1024)/usize
	);
#if 0
	ptrdiff_t usize=0;
	int fwd=0;
	int iw=0, ih=0;
	unsigned char *image=0;
	{
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		int tag=0;
		fread(&tag, 1, 2, fsrc);
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('3'|'0'<<8))
		{
			LOG_ERROR("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		if(!fwd)
		{
			LOG_ERROR("This is not a codec");
			return 1;
		}
#ifdef LOUD
		print_timestamp("%Y-%m-%d_%H%M%S\n");
#endif
		int temp=fgetc(fsrc);
		if(temp!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
		int nread=fscanf(fsrc, "%d %d", &iw, &ih);
		if(nread!=2)
		{
			LOG_ERROR("Unsupported PPM file");
			return 1;
		}
		int vmax=0;
		nread=fscanf(fsrc, "%d", &vmax);
		if(nread!=1||vmax!=255)
		{
			LOG_ERROR("Unsupported PPM file");
			return 1;
		}
		temp=fgetc(fsrc);
		if(temp!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
		if(iw<1||ih<1)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		usize=(ptrdiff_t)3*iw*ih;
		image=(unsigned char*)malloc(usize+sizeof(__m256i));
		fread(image, 1, usize, fsrc);//read image
		fclose(fsrc);
	}
	
	int rct=ch_decorrelate(image, iw, ih);

	printf("%s\n", srcfn);//
	printf("%s\n", rct_names[rct]);//

#if 1
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	double esizes[3]={0};
	for(int wp=0;wp<5;++wp)
	{
		double t=0;

		wp=6;

		printf("B01 WP%d ", wp);
		wp_estimate(image, iw, ih, rct, wp, iw, ih, esizes, &t);
		print_yuvsizes(esizes, res, t);
		
		//printf("B16 WP%d ", wp);
		//wp_estimate(image, iw, ih, rct, wp, iw/4, ih/4, esizes, &t);
		//print_yuvsizes(esizes, res, t);
		//
		//printf("B32 WP%d ", wp);
		//wp_estimate(image, iw, ih, rct, wp, iw/8, ih/4, esizes, &t);
		//print_yuvsizes(esizes, res, t);
		//
		//printf("B64 WP%d ", wp);
		//wp_estimate(image, iw, ih, rct, wp, iw/8, ih/8, esizes, &t);
		//print_yuvsizes(esizes, res, t);

		printf("\n");
		break;
	}
#endif
#if 0
	int xblocks=iw/BLOCKX;//floor
	int yblocks=ih/BLOCKY;
	int nblocks=xblocks*yblocks;
	int coeffssize=(int)sizeof(short[OLS6_NPARAMST])*nblocks;
	short *coeffs=(short*)malloc(coeffssize);
	if(!coeffs)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	for(int ky=0, blockidx=0;ky<ih;ky+=BLOCKY)
	{
		int y2=ky+BLOCKY;
		if(y2+BLOCKX>ih)
			y2=ih;
		for(int kx=0;kx<iw;kx+=BLOCKX)
		{
			int x2=kx+BLOCKX;
			if(x2+BLOCKX>iw)
				x2=iw;
			int success=solveblock(image, iw, ih, kx, x2, ky, y2, coeffs+OLS6_NPARAMST*blockidx++);
			printf("%d", success);//
			if(x2==iw)
				break;
		}
		printf("\n");//
		if(y2==ih)
			break;
	}

	ols_predict(image, iw, ih, coeffs, 1);

	double esizes0[3]={0}, esizes1[3][NCTX]={0}, tsizes1[3]={0};
	calc_esize(image, iw, ih, esizes0, (double*)esizes1);
	for(int k=0;k<NCTX;++k)
	{
		tsizes1[0]+=esizes1[0][k];
		tsizes1[1]+=esizes1[1][k];
		tsizes1[2]+=esizes1[2][k];
	}
	printf("Order-0:  ");//
	print_yuvsizes(esizes0, (ptrdiff_t)iw*ih);//
	printf("Order-C:  ");//
	print_yuvsizes(tsizes1, (ptrdiff_t)iw*ih);//
	
	ols_predict(image, iw, ih, coeffs, 0);
	ch_reconstruct(image, iw, ih, rct);
	save_ppm("20250402_0800_recon.ppm", image, iw, ih);//

	//preview coeffs
#if 1
	//AAABBBCCC	6*9 blocks
	//Ab0B00C00
	//DDDEEEFFF
	//DD0Eb0F00
	//GGGHHHIII
	//GG0HH0I00
	#define PREVX (3*3)
	#define PREVY (2*3)
//	#define PREVX 7
//	#define PREVY (OLS6_NPARAMST/PREVX)
	int ctilessize=(int)sizeof(char[PREVY*PREVX])*nblocks;
	unsigned char *ctiles=(unsigned char*)malloc(ctilessize);
	if(!ctiles)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(ctiles, 0, ctilessize);
	for(int k=0;k<OLS6_NPARAMST;++k)
	{
		int bx=0, by=0;
		switch(k)
		{
		case           0:by=  1; bx=  1;break;
		case     0+1+0+0:by=0+0; bx=0+0;break;
		case     0+1+0+1:by=0+0; bx=1+0;break;
		case     0+1+0+2:by=0+0; bx=2+0;break;
		case     0+1+0+3:by=1+0; bx=0+0;break;
		case     0+1+4+0:by=0+0; bx=0+3;break;
		case     0+1+4+1:by=0+0; bx=1+3;break;
		case     0+1+4+2:by=0+0; bx=2+3;break;
		case     0+1+4+3:by=1+0; bx=0+3;break;
		case     0+1+8+0:by=0+0; bx=0+6;break;
		case     0+1+8+1:by=0+0; bx=1+6;break;
		case     0+1+8+2:by=0+0; bx=2+6;break;
		case     0+1+8+3:by=1+0; bx=0+6;break;
		case          13:by=  3; bx=  5;break;
		case    13+1+0+0:by=0+2; bx=0+0;break;
		case    13+1+0+1:by=0+2; bx=1+0;break;
		case    13+1+0+2:by=0+2; bx=2+0;break;
		case    13+1+0+3:by=1+2; bx=0+0;break;
		case    13+1+0+4:by=1+2; bx=1+0;break;
		case    13+1+5+0:by=0+2; bx=0+3;break;
		case    13+1+5+1:by=0+2; bx=1+3;break;
		case    13+1+5+2:by=0+2; bx=2+3;break;
		case    13+1+5+3:by=1+2; bx=0+3;break;
		case    13+1+9+0:by=0+2; bx=0+6;break;
		case    13+1+9+1:by=0+2; bx=1+6;break;
		case    13+1+9+2:by=0+2; bx=2+6;break;
		case    13+1+9+3:by=1+2; bx=0+6;break;
		case       13+14:by=  5; bx=  7;break;
		case 13+14+1+0+0:by=0+4; bx=0+0;break;
		case 13+14+1+0+1:by=0+4; bx=1+0;break;
		case 13+14+1+0+2:by=0+4; bx=2+0;break;
		case 13+14+1+0+3:by=1+4; bx=0+0;break;
		case 13+14+1+0+4:by=1+4; bx=1+0;break;
		case 13+14+1+5+0:by=0+4; bx=0+3;break;
		case 13+14+1+5+1:by=0+4; bx=1+3;break;
		case 13+14+1+5+2:by=0+4; bx=2+3;break;
		case 13+14+1+5+3:by=1+4; bx=0+3;break;
		case 13+14+1+9+0:by=0+4; bx=0+6;break;
		case 13+14+1+9+1:by=0+4; bx=1+6;break;
		case 13+14+1+9+2:by=0+4; bx=2+6;break;
		case 13+14+1+9+3:by=1+4; bx=0+6;break;
		}
		//int by=k/PREVX;
		//int bx=k%PREVX;
		for(int ky=0;ky<yblocks;++ky)
		{
			for(int kx=0;kx<xblocks;++kx)
			{
				int val=coeffs[OLS6_NPARAMST*(xblocks*ky+kx)+k];
				val>>=5;
				val+=128;
				CLAMP2(val, 0, 255);
				ctiles[(PREVX*xblocks)*(yblocks*by+ky)+xblocks*bx+kx]=val;
			}
		}
	}
	save_pgm("20250402_0644_coeffs.pgm", ctiles, PREVX*xblocks, PREVY*yblocks);//
	free(ctiles);
#endif
	free(coeffs);
#endif
	free(image);
#endif
#ifndef __GNUC__
	exit(0);
#endif
//	LOG_ERROR("This is not a codec");
	(void)ols_predict;
	(void)ch_reconstruct;
	(void)solveblock;
	(void)calc_esize;
	(void)save_pgm;
	(void)save_ppm;
	(void)rct_names;
	return 0;
}
