#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
//#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>
static const char file[]=__FILE__;


	#define ALLOW_AVX2
	#define OLS4_OPTIMAL_CLAMP
//	#define OLS4_DEBUG
	#define BASE_OFFSET_ADDRESSING	//faster


#define ALLOCASSERT(C)\
	do\
		if(C)\
		{\
			LOG_ERROR("Alloc error");\
			return;\
		}\
	while(0)

//#define OLS4_RMAX 4
//#define OLS4_CTXSIZE (2*(OLS4_RMAX+1)*OLS4_RMAX)
#define PADX 8
#define PADY 8	//must be a power of two
int ols4_period=512;
double ols4_lr[4]={0.0018, 0.003, 0.003, 0.003};
//double ols4_lr[4]={0.0003, 0.0003, 0.0003, 0.0003};
unsigned char ols4_mask[4][OLS4_CTXSIZE+1]=//MSB {E3 E2 E1 E0  P3 P2 P1 P0} LSB,  last element can't exceed 1<<(kc<<1) for causality
{
#if 1
	{
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x01,0x01,0x01,0x01,0x01,0x01,0x00,
		0x00,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x00,
		0x00,0x01,0x01,0x01,0x00,
	},
	{
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x00,0x02,0x02,0x02,0x02,0x00,0x00,
		0x00,0x00,0x02,0x03,0x03,0x03,0x02,0x00,0x00,
		0x00,0x00,0x02,0x03,0x01,
	},
	{
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x00,0x04,0x04,0x04,0x04,0x00,0x00,
		0x00,0x00,0x04,0x05,0x05,0x05,0x04,0x00,0x00,
		0x00,0x00,0x04,0x05,0x01,
	},
	{
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x00,0x00,0x00,
	},
#endif
#if 0
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00,
		0x00, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00,
		0x00, 0x01, 0x01, 0x01, 0x00,
	},
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x02, 0x02, 0x02, 0x02, 0x00, 0x00,
		0x00, 0x00, 0x02, 0x03, 0x03, 0x03, 0x02, 0x02, 0x00,
		0x00, 0x00, 0x02, 0x03, 0x01,
	},
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x04, 0x04, 0x04, 0x04, 0x00, 0x00,
		0x00, 0x00, 0x04, 0x05, 0x05, 0x05, 0x04, 0x00, 0x00,
		0x00, 0x00, 0x04, 0x05, 0x01,
	},
#endif
#if 0
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x01, 0x01, 0x07, 0x07, 0x01, 0x00, 0x00,
		0x00, 0x00, 0x01, 0x07, 0x00,
	},
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x02, 0x02, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x02, 0x02, 0x07, 0x07, 0x02, 0x00, 0x00,
		0x00, 0x00, 0x02, 0x07, 0x01,
	},
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x04, 0x04, 0x07, 0x07, 0x04, 0x00, 0x00,
		0x00, 0x00, 0x04, 0x07, 0x03,
	},
#endif
#if 0
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00,
		0x00, 0x00, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
		0x00, 0x01, 0x01, 0x07, 0x07, 0x07, 0x01, 0x01, 0x01,
		0x01, 0x01, 0x01, 0x07, 0x00,
	},
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x02, 0x02, 0x02, 0x02, 0x02, 0x00,
		0x00, 0x00, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02,
		0x00, 0x02, 0x02, 0x03, 0x07, 0x07, 0x02, 0x02, 0x02,
		0x02, 0x02, 0x02, 0x07, 0x01,
	},
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x04, 0x04, 0x04, 0x04, 0x04, 0x00,
		0x00, 0x00, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04,
		0x00, 0x04, 0x04, 0x05, 0x07, 0x05, 0x04, 0x04, 0x04,
		0x04, 0x04, 0x04, 0x07, 0x03,
	},
#endif
#if 0
	{
		0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
		0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
		0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
		0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
		0x07, 0x07, 0x07, 0x07, 0x00,
	},
	{
		0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
		0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
		0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
		0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
		0x07, 0x07, 0x07, 0x07, 0x01,
	},
	{
		0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
		0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
		0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
		0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
		0x07, 0x07, 0x07, 0x07, 0x03,
	},
#endif
#if 0
	{
		0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
		0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
		0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
		0x01, 0x01, 0x01, 0x01, 0x07, 0x07, 0x01, 0x01, 0x01,
		0x01, 0x01, 0x01, 0x07, 0x00,
	},
	{
		0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02,
		0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02,
		0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02,
		0x02, 0x02, 0x02, 0x02, 0x07, 0x07, 0x02, 0x02, 0x02,
		0x02, 0x02, 0x02, 0x07, 0x01,
	},
	{
		0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04,
		0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04,
		0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04,
		0x04, 0x04, 0x04, 0x04, 0x07, 0x07, 0x04, 0x04, 0x04,
		0x04, 0x04, 0x04, 0x07, 0x03,
	},
#endif
#if 0
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x07, 0x07, 0x07, 0x07, 0x07, 0x00, 0x00,
		0x00, 0x00, 0x07, 0x07, 0x07, 0x07, 0x07, 0x00, 0x00,
		0x00, 0x00, 0x07, 0x07, 0x00,
	},
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x07, 0x07, 0x07, 0x07, 0x07, 0x00, 0x00,
		0x00, 0x00, 0x07, 0x07, 0x07, 0x07, 0x07, 0x00, 0x00,
		0x00, 0x00, 0x07, 0x07, 0x01,
	},
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x07, 0x07, 0x07, 0x07, 0x07, 0x00, 0x00,
		0x00, 0x00, 0x07, 0x07, 0x07, 0x07, 0x07, 0x00, 0x00,
		0x00, 0x00, 0x07, 0x07, 0x03,
	},
#endif
#if 0
	{
		0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x01, 0x07, 0x07, 0x07, 0x01, 0x00, 0x00,
		0x00, 0x00, 0x01, 0x07, 0x00,
	},
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x02, 0x02, 0x02, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x02, 0x07, 0x07, 0x07, 0x02, 0x00, 0x00,
		0x02, 0x02, 0x02, 0x07, 0x01,
	},
	{
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x04, 0x04, 0x04, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x04, 0x07, 0x07, 0x07, 0x04, 0x04, 0x04,
		0x00, 0x00, 0x04, 0x07, 0x03,
	},
#endif
};
void pred_ols4(Image *src, int period, double *lrs, unsigned char *mask0, unsigned char *mask1, unsigned char *mask2, unsigned char *mask3, int fwd)
{
	double t_start=time_sec();
	if(loud_transforms)
		DisableProcessWindowsGhosting();
	unsigned char *masks[]=
	{
		mask0,
		mask1,
		mask2,
		mask3,
	};
	int ctxsize[4]={0}, maxcount=0;
	for(int kc=0;kc<4;++kc)
	{
		for(int km=0;km<OLS4_CTXSIZE+1;++km)
		{
			int val=masks[kc][km];
			//if(km==OLS4_CTXSIZE)//causality mask
			//{
			//	int cmask=((1<<(kc<<1))-1);
			//	val&=cmask<<4|cmask;
			//}
			val-=val>>1&0x55;
			val=(val>>2&0x33)+(val&0x33);
			val=(val>>4)+(val&15);
			ctxsize[kc]+=val;
			//ctxsize[kc]+=hammingweight16(masks[kc][km]);
		}
		UPDATE_MAX(maxcount, ctxsize[kc]);
	}
	double **ctx[4]={0}, *vec[4]={0};
	int matsize[4]={0};
	double *cov[4]={0}, *cholesky[4]={0};
	double *params[4]={0};
	for(int kc=0;kc<4;++kc)
	{
		if(ctxsize[kc])
		{
			matsize[kc]=ctxsize[kc]*ctxsize[kc];
			ctx[kc]=(double**)_mm_malloc(sizeof(double*)*ctxsize[kc], sizeof(__m256i));
			vec[kc]=(double*)_mm_malloc(sizeof(double)*(ctxsize[kc]+4LL), sizeof(__m256d));
			cov[kc]=(double*)_mm_malloc(sizeof(double)*(matsize[kc]+4LL), sizeof(__m256d));
			cholesky[kc]=(double*)malloc(sizeof(double)*matsize[kc]);
			params[kc]=(double*)_mm_malloc(sizeof(double)*ctxsize[kc], sizeof(__m256d));
			ALLOCASSERT(!ctx[kc]||!vec[kc]||!cov[kc]||!cholesky[kc]||!params[kc]);
			memset(ctx[kc], 0, sizeof(double*)*ctxsize[kc]);
			memset(vec[kc], 0, sizeof(double)*(ctxsize[kc]+4LL));
			memset(cov[kc], 0, sizeof(double)*(matsize[kc]+4LL));
			memset(cholesky[kc], 0, sizeof(double)*matsize[kc]);
			memset(params[kc], 0, sizeof(double)*ctxsize[kc]);
		}
	}
	double *ctxtemp=(double*)_mm_malloc((maxcount+4LL)*sizeof(double), sizeof(__m256d));

	size_t bufsize=(src->iw+PADX*2LL)*sizeof(double[4*PADY*2]);//PADY padded rows * 4 channels * {pixels, errors}
	double *pixels=(double*)malloc(bufsize);
	ALLOCASSERT(!pixels||!ctxtemp);
	memset(pixels, 0, bufsize);
	memset(ctxtemp, 0, (maxcount+4LL)*sizeof(double));

	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int half[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int successcount=0;
#ifdef ALLOW_AVX2
	__m256d mlr[]=
	{
		_mm256_set1_pd(lrs[0]),
		_mm256_set1_pd(lrs[1]),
		_mm256_set1_pd(lrs[2]),
		_mm256_set1_pd(lrs[3]),
	};
	__m256d mlr_comp[]=
	{
		_mm256_set1_pd(1-lrs[0]),
		_mm256_set1_pd(1-lrs[1]),
		_mm256_set1_pd(1-lrs[2]),
		_mm256_set1_pd(1-lrs[3]),
	};
#ifndef BASE_OFFSET_ADDRESSING
	__m256i stride=_mm256_set1_epi64x(sizeof(double[8]));
#endif
#endif
	ALIGN(16) int ipred[4]={0};
	int olsnextpoint=1;
	for(int ky=0, idx=0, olsidx=1;ky<src->ih;++ky)
	{
		double *rows[PADY];
		for(int k=0;k<PADY;++k)
			rows[k]=pixels+(((src->iw+PADX*2LL)*(((size_t)ky-k)&(PADY-1))+PADX)<<3);
		for(int kc=0;kc<4;++kc)
		{
			int nparams=ctxsize[kc];
			if(!nparams)
				continue;
			const unsigned char *ctxcell=masks[kc];
			double **curr_ctx=ctx[kc];
			for(int ky2=-OLS4_RMAX, loadidx=0;ky2<=0;++ky2)
			{
				for(int kx2=-OLS4_RMAX;kx2<=OLS4_RMAX;++kx2, ++ctxcell)
				{
					double *comp=rows[-ky2]+((size_t)kx2<<3);
					for(int kc2=0;kc2<8;++kc2)
					{
#ifdef OLS4_DEBUG
						if(ctxcell-masks[kc]>=OLS4_CTXSIZE+1)
							LOG_ERROR("");
#endif
						if(*ctxcell>>kc2&1)
						{
#ifdef OLS4_DEBUG
							if(loadidx>=nparams)//
								LOG_ERROR("");
#endif
							curr_ctx[loadidx]=comp+kc2;
							++loadidx;
							if(loadidx>=nparams)
								goto finish_loading;
						}
					}
#ifdef OLS4_DEBUG
					if(!ky2&&!kx2)
					{
						if(loadidx!=nparams)
							LOG_ERROR("");
						goto finish_loading;
					}
#endif
				}
			}
	finish_loading:
			;
		}
		for(int kx=0;kx<src->iw;++kx, ++olsidx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				int nparams=ctxsize[kc];
				if(!nparams)
					continue;
				double **curr_ctx=ctx[kc];
				double *curr_vec=vec[kc];
				double *curr_cov=cov[kc];
				double *curr_cholesky=cholesky[kc];
				double *curr_params=params[kc];
				double fpred;
#ifdef OLS4_OPTIMAL_CLAMP
				double temp=rows[0][kc-8];//W
				double cmin=temp, cmax=temp;
				temp=rows[1][kc];//N
				UPDATE_MIN(cmin, temp);
				UPDATE_MAX(cmax, temp);
				temp=rows[1][kc+8];//NE
				UPDATE_MIN(cmin, temp);
				UPDATE_MAX(cmax, temp);
#endif
				fpred=0;
#if 0
				//THIS IS SLOWER
				for(int k=0;k<nparams;++k)
					ctxtemp[k]=*curr_ctx[k];
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");
				{
					int k=0;
					__m256d msum=_mm256_setzero_pd();
					for(;k<nparams-3;k+=4)
					{
						__m256d mp=_mm256_load_pd(curr_params+k);
						__m256d mc=_mm256_load_pd(ctxtemp+k);
						mp=_mm256_mul_pd(mp, mc);
						msum=_mm256_add_pd(msum, mp);
					}
					ALIGN(32) double asum[4];
					msum=_mm256_hadd_pd(msum, _mm256_setzero_pd());
					_mm256_store_pd(asum, msum);
					for(;k<nparams;++k)
						fpred+=curr_params[k]*ctxtemp[k];
					fpred+=asum[0]+asum[2];
				}
#elif 1
				//FASTEST
				{
					int k=nparams-1;
					//int k=0, n=nparams-1;
					double **p1=curr_ctx;
					double *p2=curr_params;
					double *p3=ctxtemp;
					double fpred2=0;
					for(;k>=1;k-=2)
					//for(;k<n;k+=2)
					{
#ifdef BASE_OFFSET_ADDRESSING
						double nb0=p1[0][kx<<3], nb1=p1[1][kx<<3];
#else
						double nb0=*p1[0], nb1=*p1[1];
#endif
						fpred+=nb0*p2[0];
						fpred2+=nb1*p2[1];
						p3[0]=nb0;
						p3[1]=nb1;
						p1+=2;
						p2+=2;
						p3+=2;
					}
					if(k>=0)
					//if(k<nparams)
					{
#ifdef BASE_OFFSET_ADDRESSING
						double nb0=p1[0][kx<<3];
#else
						double nb0=*p1[0];
#endif
						fpred+=nb0*p2[0];
						p3[0]=nb0;
					}
					fpred+=fpred2;
				}
#else
				for(int k=0;k<nparams;++k)
				{
					double nb=*curr_ctx[k];
					fpred+=curr_params[k]*nb;
					ctxtemp[k]=nb;
				}
#endif
#ifdef OLS4_DEBUG
				double fpred0=fpred;//
#endif
#ifdef OLS4_OPTIMAL_CLAMP
				fpred=CLAMP(cmin, fpred, cmax);
#else
				fpred=CLAMP(-half[kc], fpred, half[kc]-1);
#endif
				__m128d fp=_mm_set_sd(fpred);
				__m128i ip=_mm_cvtpd_epi32(fp);
				_mm_store_si128((__m128i*)ipred, ip);
				int pred=ipred[0];
				//int pred=(int)round(fpred);

				int val=src->data[idx];
#ifdef OLS4_DEBUG
				if(!kc&&kx==4&&ky==4)//
					printf("");
#endif

				if(fwd)
				{
					rows[0][kc+0]=val;//pixel
					val-=pred;
					val+=half[kc];
					val&=nlevels[kc]-1;
					val-=half[kc];
				}
				else
				{
					val+=pred;
					val+=half[kc];
					val&=nlevels[kc]-1;
					val-=half[kc];
					rows[0][kc+0]=val;//pixel
				}
				src->data[idx]=val;
				rows[0][kc+4]=rows[0][kc+0]-fpred;//high-res error
				if(olsidx==olsnextpoint)
				//if(!(olsidx&63))
				{
				double lr=lrs[kc];
#ifdef ALLOW_AVX2
				//if((nparams>>1<<1)!=nparams)
				//	printf("");
				for(int ky2=0, midx=0;ky2<nparams;++ky2)
				{
					__m256d yctx=_mm256_set1_pd(ctxtemp[ky2]);
					int kx2=0;
					for(;kx2<nparams-3;kx2+=4, midx+=4)
					{
						__m256d mcov=_mm256_load_pd(curr_cov+midx);
						__m256d xctx=_mm256_load_pd(ctxtemp+kx2);
						xctx=_mm256_mul_pd(xctx, yctx);
						xctx=_mm256_sub_pd(xctx, mcov);
						xctx=_mm256_mul_pd(xctx, mlr[kc]);
						xctx=_mm256_add_pd(xctx, mcov);
						_mm256_store_pd(curr_cov+midx, xctx);
						//curr_cov[midx+0]+=(ctxtemp[kx2+0]*ctxtemp[ky2]-curr_cov[midx+0])*lr;
						//curr_cov[midx+1]+=(ctxtemp[kx2+1]*ctxtemp[ky2]-curr_cov[midx+1])*lr;
						//curr_cov[midx+2]+=(ctxtemp[kx2+2]*ctxtemp[ky2]-curr_cov[midx+2])*lr;
						//curr_cov[midx+3]+=(ctxtemp[kx2+3]*ctxtemp[ky2]-curr_cov[midx+3])*lr;
						//for(int k=0;k<4;++k)
						//{
						//	if(abs(curr_cov[midx+k]-xctx.m256d_f64[k])>1e-6)
						//		printf("");
						//}
					}
					for(;kx2<nparams;++kx2, ++midx)
						curr_cov[midx]+=(ctxtemp[kx2]*ctxtemp[ky2]-curr_cov[midx])*lr;
				}
#else
				for(int ky2=0, midx=0;ky2<nparams;++ky2)
				{
					for(int kx2=0;kx2<nparams;++kx2, ++midx)
						curr_cov[midx]+=(ctxtemp[kx2]*ctxtemp[ky2]-curr_cov[midx])*lr;
				}
#endif
#if 1
#ifdef ALLOW_AVX2
				__m256d mval=_mm256_set1_pd(rows[0][kc+0]*lr);
				double lval=rows[0][kc+0]*lr, lr_comp=1-lr;//
				{
					int k=0;
					for(;k<nparams-3;k+=4)
					{
						__m256d mvec=_mm256_load_pd(curr_vec+k);
						__m256d mtmp=_mm256_load_pd(ctxtemp+k);
						mtmp=_mm256_mul_pd(mval, mtmp);
						mvec=_mm256_mul_pd(mvec, mlr_comp[kc]);
						mvec=_mm256_add_pd(mvec, mtmp);
						_mm256_store_pd(curr_vec+k, mvec);
						//curr_vec[k+0]=lval*ctxtemp[k+0]+lr_comp*curr_vec[k+0];
						//curr_vec[k+1]=lval*ctxtemp[k+1]+lr_comp*curr_vec[k+1];
						//curr_vec[k+2]=lval*ctxtemp[k+2]+lr_comp*curr_vec[k+2];
						//curr_vec[k+3]=lval*ctxtemp[k+3]+lr_comp*curr_vec[k+3];
					}
					for(;k<nparams;++k)
						curr_vec[k]=lval*ctxtemp[k]+lr_comp*curr_vec[k];
				}
#else
				double lval=rows[0][kc+0]*lr, lr_comp=1-lr;
				for(int k=0;k<nparams;++k)
					curr_vec[k]=lval*ctxtemp[k]+lr_comp*curr_vec[k];
#endif
#else
				double lval=rows[0][kc+0];
				for(int k=0;k<nparams;++k)
					curr_vec[k]+=(lval*ctxtemp[k]-curr_vec[k])*lr;
#endif
				//if(olsidx==period)//OLS solver by Cholesky decomposition from paq8px
				{
					int success=1;
					double sum;
					int n=nparams;
					memcpy(curr_cholesky, curr_cov, sizeof(double)*matsize[kc]);
					for(int k=0;k<matsize[kc];k+=n+1)
						curr_cholesky[k]+=0.0075;
					for(int i=0;i<n;++i)
					{
						for(int j=0;j<i;++j)
						{
							sum=curr_cholesky[i*n+j];
							for(int k=0;k<j;++k)
								sum-=curr_cholesky[i*n+k]*curr_cholesky[j*n+k];
							curr_cholesky[i*n+j]=sum/curr_cholesky[j*n+j];
						}
						sum=curr_cholesky[i*n+i];
						for(int k=0;k<i;++k)
							sum-=curr_cholesky[i*n+k]*curr_cholesky[i*n+k];
						if(sum<=1e-8)
						{
							success=0;
							break;
						}
						curr_cholesky[i*n+i]=sqrt(sum);
					}
					if(success)
					{
						for(int i=0;i<n;++i)
						{
							sum=curr_vec[i];
							for(int j=0;j<i;++j)
								sum-=curr_cholesky[i*n+j]*curr_params[j];
							curr_params[i]=sum/curr_cholesky[i*n+i];
						}
						for(int i=n-1;i>=0;--i)
						{
							sum=curr_params[i];
							for(int j=i+1;j<n;++j)
								sum-=curr_cholesky[j*n+i]*curr_params[j];
							curr_params[i]=sum/curr_cholesky[i*n+i];
						}
						//for(int k=0;k<n;++k)//X
						//	curr_params[k]+=(ctxtemp[k]-curr_params[k])*0.9;
						++successcount;
					}
				}
				}
#ifndef BASE_OFFSET_ADDRESSING
#ifdef ALLOW_AVX2
				{
					int k=nparams;
					while(k>=32)
					{
						__m256i mp0=_mm256_load_si256((__m256i*)curr_ctx+0);
						__m256i mp1=_mm256_load_si256((__m256i*)curr_ctx+1);
						__m256i mp2=_mm256_load_si256((__m256i*)curr_ctx+2);
						__m256i mp3=_mm256_load_si256((__m256i*)curr_ctx+3);
						__m256i mp4=_mm256_load_si256((__m256i*)curr_ctx+4);
						__m256i mp5=_mm256_load_si256((__m256i*)curr_ctx+5);
						__m256i mp6=_mm256_load_si256((__m256i*)curr_ctx+6);
						__m256i mp7=_mm256_load_si256((__m256i*)curr_ctx+7);
						mp0=_mm256_add_epi64(mp0, stride);
						mp1=_mm256_add_epi64(mp1, stride);
						mp2=_mm256_add_epi64(mp2, stride);
						mp3=_mm256_add_epi64(mp3, stride);
						mp4=_mm256_add_epi64(mp4, stride);
						mp5=_mm256_add_epi64(mp5, stride);
						mp6=_mm256_add_epi64(mp6, stride);
						mp7=_mm256_add_epi64(mp7, stride);
						_mm256_store_si256((__m256i*)curr_ctx+0, mp0);
						_mm256_store_si256((__m256i*)curr_ctx+1, mp1);
						_mm256_store_si256((__m256i*)curr_ctx+2, mp2);
						_mm256_store_si256((__m256i*)curr_ctx+3, mp3);
						_mm256_store_si256((__m256i*)curr_ctx+4, mp4);
						_mm256_store_si256((__m256i*)curr_ctx+5, mp5);
						_mm256_store_si256((__m256i*)curr_ctx+6, mp6);
						_mm256_store_si256((__m256i*)curr_ctx+7, mp7);
						k-=32;
						curr_ctx+=32;
					}
					if(k>=16)
					{
						__m256i mp0=_mm256_load_si256((__m256i*)curr_ctx+0);
						__m256i mp1=_mm256_load_si256((__m256i*)curr_ctx+1);
						__m256i mp2=_mm256_load_si256((__m256i*)curr_ctx+2);
						__m256i mp3=_mm256_load_si256((__m256i*)curr_ctx+3);
						mp0=_mm256_add_epi64(mp0, stride);
						mp1=_mm256_add_epi64(mp1, stride);
						mp2=_mm256_add_epi64(mp2, stride);
						mp3=_mm256_add_epi64(mp3, stride);
						_mm256_store_si256((__m256i*)curr_ctx+0, mp0);
						_mm256_store_si256((__m256i*)curr_ctx+1, mp1);
						_mm256_store_si256((__m256i*)curr_ctx+2, mp2);
						_mm256_store_si256((__m256i*)curr_ctx+3, mp3);
						k-=16;
						curr_ctx+=16;
					}
					if(k>=8)
					{
						__m256i mp0=_mm256_load_si256((__m256i*)curr_ctx+0);
						__m256i mp1=_mm256_load_si256((__m256i*)curr_ctx+1);
						mp0=_mm256_add_epi64(mp0, stride);
						mp1=_mm256_add_epi64(mp1, stride);
						_mm256_store_si256((__m256i*)curr_ctx+0, mp0);
						_mm256_store_si256((__m256i*)curr_ctx+1, mp1);
						k-=8;
						curr_ctx+=8;
					}
					if(k>=4)
					{
						__m256i mp0=_mm256_load_si256((__m256i*)curr_ctx+0);
						mp0=_mm256_add_epi64(mp0, stride);
						_mm256_store_si256((__m256i*)curr_ctx+0, mp0);
						k-=4;
						curr_ctx+=4;
					}
					if(k>=2)
					{
						curr_ctx[0]+=8;
						curr_ctx[1]+=8;
						k-=2;
						curr_ctx+=2;
					}
					if(k>=1)
						*curr_ctx+=8;
				}
#else
				for(int k=0;k<nparams;++k)
					curr_ctx[k]+=8;
#endif
#endif
			}
			if(olsidx==olsnextpoint)
			{
				olsnextpoint<<=olsnextpoint<period;
				olsidx=0;
			}
			//olsidx&=-(olsidx<period);
			//if(olsidx==period)
			//	olsidx=0;
			for(int k=0;k<PADY;++k)
				rows[k]+=8;
		}
		if(loud_transforms&&(ky&63))
		{
			double elapsed=time_sec()-t_start;
			if(elapsed>1)
				set_window_title("%d/%d = %7.3lf%%  OLS4 rate %lf%%  %lf sec", ky+1, src->ih, 100.*(ky+1)/src->ih, 100.*successcount*period/(src->nch*src->iw*(ky+1)), elapsed);
		}
	}
	if(loud_transforms)
		set_window_title("OLS4 %lf%%  %lf sec {%d %d %d %d} params",
			100.*successcount*period/(src->nch*src->iw*src->ih),
			time_sec()-t_start,
			ctxsize[0],
			ctxsize[1],
			ctxsize[2],
			ctxsize[3]
		);
	for(int kc=0;kc<4;++kc)
	{
		if(ctxsize[kc])
		{
			_mm_free(ctx[kc]);
			_mm_free(vec[kc]);
			_mm_free(cov[kc]);
			free(cholesky[kc]);
			_mm_free(params[kc]);
		}
	}
	_mm_free(ctxtemp);
	free(pixels);
}