#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>


	#define ALLOW_AVX2
	#define OLS4_OPTIMAL_CLAMP
	#define CLAMPGRAD
#if 0
	#define OLS4_DEBUG
#endif
	#define BASE_OFFSET_ADDRESSING	/*faster*/


#define ALLOCASSERT(C)\
	do\
		if(C)\
		{\
			CRASH("Alloc error");\
			return;\
		}\
	while(0)

#if 0
#	define OLS4_RMAX 4
#	define OLS4_CTXSIZE (2*(OLS4_RMAX+1)*OLS4_RMAX)
#endif
#define PADX 8
#define PADY 8	/*must be a power of two*/
int ols4_period=512;
double ols4_lr[4]=
{
	0.0018, 0.003, 0.003, 0.003
#if 0
	0.0003, 0.0003, 0.0003, 0.0003
#endif
};
uint8_t ols4_cache=0;
uint8_t ols4_mask[4][OLS4_CTXSIZE+1]=/*MSB {E3 E2 E1 E0  P3 P2 P1 P0} LSB,  last element can't exceed 1<<(kc<<1) for causality*/
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
#if 0
static int solve_Mx_v(double *matrix, double *vec, int size)
{
	int it, success, k, kx, ky;

	success=1;
	for(it=0;it<size;++it)
	{
		int kp;
		double pivot;

		kp=it;
		for(;kp<size;++kp)
		{
			if(fabs(matrix[size*kp+it])>1e-6)
				break;
		}
		if(kp==size)
		{
			success=0;
			break;
		}
		if(kp!=it)
		{
			double temp;
			for(k=it;k<size;++k)
				SWAPVAR(matrix[size*it+k], matrix[size*kp+k], temp);
			SWAPVAR(vec[it], vec[kp], temp);
		}
		pivot=matrix[size*it+it];
		for(kx=it;kx<size;++kx)
			matrix[size*it+kx]/=pivot;
		vec[it]/=pivot;
		for(ky=0;ky<size;++ky)
		{
			double factor;

			if(ky==it)
				continue;
			factor=matrix[size*ky+it];
			if(fabs(factor)>1e-6)
			{
				for(kx=it;kx<size;++kx)
					matrix[size*ky+kx]-=matrix[size*it+kx]*factor;
				vec[ky]-=vec[it]*factor;
			}
		}
	}
	return success;
}
#endif
void pred_ols4(Image *src, int period, double *lrs, uint8_t *mask0, uint8_t *mask1, uint8_t *mask2, uint8_t *mask3, int fwd)
{
	double t_start=time_sec();
	uint8_t *masks[4];
	int ctxsize[4]={0}, maxcount=0;
	double **ctx[4]={0}, *vec[4]={0};
#if 0
	double *vec2[4]={0};
#endif
	int matsize[4]={0};
	double *cov[4]={0}, *cholesky[4]={0};
	double *params[4]={0};
	double *ctxtemp;
	size_t bufsize;
	double *pixels;
	int nlevels[4], half[4], successcount;
	__m256d mlr[4], mlr_comp[4];
#ifndef BASE_OFFSET_ADDRESSING
	__m256i stride;
#endif
	ALIGN(16) int ipred[4]={0};
	int olsnextpoint;
	double *rows[PADY]={0};
	int kc, km, kx, ky, idx, olsidx, i, j, k, ky2, kx2, loadidx, kc2, midx;


	masks[0]=mask0;
	masks[1]=mask1;
	masks[2]=mask2;
	masks[3]=mask3;
	nlevels[0]=1<<src->depth[0];
	nlevels[1]=1<<src->depth[1];
	nlevels[2]=1<<src->depth[2];
	nlevels[3]=1<<src->depth[3];
	half[0]=nlevels[0]>>1;
	half[1]=nlevels[1]>>1;
	half[2]=nlevels[2]>>1;
	half[3]=nlevels[3]>>1;
	successcount=0;
#ifdef ALLOW_AVX2
	mlr[0]=_mm256_set1_pd(lrs[0]);
	mlr[1]=_mm256_set1_pd(lrs[1]);
	mlr[2]=_mm256_set1_pd(lrs[2]);
	mlr[3]=_mm256_set1_pd(lrs[3]);
	mlr_comp[0]=_mm256_set1_pd(1-lrs[0]);
	mlr_comp[1]=_mm256_set1_pd(1-lrs[1]);
	mlr_comp[2]=_mm256_set1_pd(1-lrs[2]);
	mlr_comp[3]=_mm256_set1_pd(1-lrs[3]);
#ifndef BASE_OFFSET_ADDRESSING
	stride=_mm256_set1_epi64x(sizeof(double[8]));
#endif
#endif
	olsnextpoint=1;
	if(loud_transforms)
		DisableProcessWindowsGhosting();
	for(kc=0;kc<4;++kc)
	{
		for(km=0;km<OLS4_CTXSIZE+1;++km)
		{
			int val=masks[kc][km];
			val-=val>>1&0x55;
			val=(val>>2&0x33)+(val&0x33);
			val=(val>>4)+(val&15);
			ctxsize[kc]+=val;
		}
#ifdef CLAMPGRAD
		ctxsize[kc]+=src->depth[kc]!=0;
#endif
		UPDATE_MAX(maxcount, ctxsize[kc]);
	}
	for(kc=0;kc<4;++kc)
	{
		if(ctxsize[kc])
		{
			matsize[kc]=ctxsize[kc]*ctxsize[kc];
			ctx[kc]=(double**)_mm_malloc(sizeof(double*)*ctxsize[kc], sizeof(__m256i));
			vec[kc]=(double*)_mm_malloc(sizeof(double)*(ctxsize[kc]+4), sizeof(__m256d));
#if 0
			vec2[kc]=(double*)_mm_malloc(sizeof(double)*(ctxsize[kc]+4), sizeof(__m256d));
#endif
			cov[kc]=(double*)_mm_malloc(sizeof(double)*(matsize[kc]+4), sizeof(__m256d));
			cholesky[kc]=(double*)malloc(sizeof(double)*matsize[kc]);
			params[kc]=(double*)_mm_malloc(sizeof(double)*ctxsize[kc], sizeof(__m256d));
			ALLOCASSERT(!ctx[kc]||!vec[kc]||!cov[kc]||!cholesky[kc]||!params[kc]);
			memset(ctx[kc], 0, sizeof(double*)*ctxsize[kc]);
			memset(vec[kc], 0, sizeof(double)*(ctxsize[kc]+4));
#if 0
			memset(vec2[kc], 0, sizeof(double)*(ctxsize[kc]+4));
#endif
			memset(cov[kc], 0, sizeof(double)*(matsize[kc]+4));
			memset(cholesky[kc], 0, sizeof(double)*matsize[kc]);
			memset(params[kc], 0, sizeof(double)*ctxsize[kc]);
		}
	}
	ctxtemp=(double*)_mm_malloc((maxcount+4)*sizeof(double), sizeof(__m256d));

	bufsize=(src->iw+PADX*2)*sizeof(double[4*PADY*2]);/*PADY padded rows * 4 channels * {pixels, errors}*/
	pixels=(double*)malloc(bufsize);
	ALLOCASSERT(!pixels||!ctxtemp);
	memset(pixels, 0, bufsize);
	memset(ctxtemp, 0, (maxcount+4)*sizeof(double));

	for(ky=0, idx=0, olsidx=1;ky<src->ih;++ky)
	{
		for(k=0;k<PADY;++k)
			rows[k]=pixels+(((src->iw+PADX*2)*(((size_t)ky-k)&(PADY-1))+PADX)<<3);
		for(kc=0;kc<4;++kc)
		{
			const uint8_t *ctxcell;
			double **curr_ctx;
			int nparams;

			nparams=ctxsize[kc];
			if(!nparams)
				continue;
			ctxcell=masks[kc];
			curr_ctx=ctx[kc];
			for(ky2=-OLS4_RMAX, loadidx=0;ky2<=0;++ky2)
			{
				for(kx2=-OLS4_RMAX;kx2<=OLS4_RMAX;++kx2, ++ctxcell)
				{
					double *comp=rows[-ky2]+((size_t)kx2<<3);
					for(kc2=0;kc2<8;++kc2)
					{
#ifdef OLS4_DEBUG
						if(ctxcell-masks[kc]>=OLS4_CTXSIZE+1)
							CRASH("");
#endif
						if(*ctxcell>>kc2&1)
						{
#ifdef OLS4_DEBUG
							if(loadidx>=nparams)/**/
								CRASH("");
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
							CRASH("");
						goto finish_loading;
					}
#endif
				}
			}
	finish_loading:
			;
		}
		for(kx=0;kx<src->iw;++kx, ++olsidx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				int nparams=ctxsize[kc];
				double
					**curr_ctx,
					*curr_vec,
					*curr_cov,
					*curr_cholesky,
					*curr_params,
					fpred;
#ifdef CLAMPGRAD
				ALIGN(16) double cgrad[2], cmin[2], cmax[2];
#endif

				if(!nparams)
					continue;
				curr_ctx=ctx[kc];
				curr_vec=vec[kc];
				curr_cov=cov[kc];
				curr_cholesky=cholesky[kc];
				curr_params=params[kc];
#ifdef CLAMPGRAD
				{
					__m128d a=_mm_set_pd(0, rows[1][kc]);/*N*/
					__m128d b=_mm_set_pd(0, rows[0][kc-8]);/*W*/
					__m128d c=_mm_set_pd(0, rows[1][kc+8]);/*NE*/
					__m128d x=_mm_set_pd(0, rows[1][kc]+rows[0][kc-8]-rows[1][kc-8]);/*N+W-NW*/
					__m128d vmin=_mm_min_pd(a, b);
					__m128d vmax=_mm_max_pd(a, b);
					x=_mm_max_pd(x, vmin);
					x=_mm_min_pd(x, vmax);
					vmin=_mm_min_pd(vmin, c);
					vmax=_mm_max_pd(vmax, c);
					_mm_store_pd(cgrad, x);
					_mm_store_pd(cmin, vmin);
					_mm_store_pd(cmax, vmax);
				}
#if 0
				cgrad=rows[1][kc]+rows[0][kc-8]-rows[1][kc-8];
				if(rows[1][kc]<rows[0][kc-8])
					cmin=rows[1][kc], cmax=rows[0][kc-8];/*N, W*/
				else
					cmin=rows[0][kc-8], cmax=rows[1][kc];
				cgrad=CLAMP(cmin, cgrad, cmax);
				{
					double temp=rows[1][kc+8];/*NE*/
					UPDATE_MIN(cmin, temp);
					UPDATE_MAX(cmax, temp);
				}
				{
					double temp=rows[1][kc-8];/*NW	X*/
					UPDATE_MIN(cmin, temp);
					UPDATE_MAX(cmax, temp);
				}
#endif
				--nparams;
#else
#ifdef OLS4_OPTIMAL_CLAMP
				double temp=rows[0][kc-8];/*W*/
				double cmin=temp, cmax=temp;
				temp=rows[1][kc];/*N*/
				UPDATE_MIN(cmin, temp);
				UPDATE_MAX(cmax, temp);
				temp=rows[1][kc+8];/*NE*/
				UPDATE_MIN(cmin, temp);
				UPDATE_MAX(cmax, temp);
#endif
#endif
				fpred=0;
#if 0
				/*THIS IS SLOWER*/
				for(k=0;k<nparams;++k)
					ctxtemp[k]=*curr_ctx[k];
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
				/*FASTEST*/
				{
					double **p1, *p2, *p3, fpred2;

					k=nparams-1;
					p1=curr_ctx;
					p2=curr_params;
					p3=ctxtemp;
					fpred2=0;

					for(;k>=1;k-=2)
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
					{
						double nb0;
#ifdef BASE_OFFSET_ADDRESSING
						nb0=p1[0][kx<<3];
#else
						nb0=*p1[0];
#endif
						fpred+=nb0*p2[0];
						p3[0]=nb0;
					}
					fpred+=fpred2;
				}
#else
				for(k=0;k<nparams;++k)
				{
					double nb=*curr_ctx[k];
					fpred+=curr_params[k]*nb;
					ctxtemp[k]=nb;
				}
#endif
#ifdef CLAMPGRAD
				ctxtemp[nparams]=cgrad[0];
				fpred+=cgrad[0]*curr_params[nparams];
				++nparams;
#endif
#ifdef OLS4_DEBUG
				double fpred0=fpred;
#endif
#ifdef OLS4_OPTIMAL_CLAMP
				CLAMP2_fp64(fpred, fpred, cmin[0], cmax[0]);
#else
				fpred=CLAMP(-half[kc], fpred, half[kc]-1);
#endif
				{
					__m128d fp;
					__m128i ip;
					int pred, val;

					fp=_mm_set_sd(fpred);
					ip=_mm_cvtpd_epi32(fp);
					_mm_store_si128((__m128i*)ipred, ip);
					pred=ipred[0];
					val=((int32_t*)HDR_DATA(src))[idx];

					if(fwd)
					{
						rows[0][kc+0]=val;/*pixel*/
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
						rows[0][kc+0]=val;/*pixel*/
					}
					((int32_t*)HDR_DATA(src))[idx]=val;
				}
				rows[0][kc+4]=rows[0][kc+0]-fpred;/*high-res error*/
#if 0
				rows[0][kc+4]=fpred;/*high-res pred*/
#endif
				if(olsidx==olsnextpoint)
				{
				double lr=lrs[kc];
#ifdef ALLOW_AVX2
				for(ky2=0, midx=0;ky2<nparams;++ky2)
				{
					__m256d yctx=_mm256_set1_pd(ctxtemp[ky2]);

					for(kx2=0;kx2<nparams-3;kx2+=4, midx+=4)
					{
						__m256d mcov=_mm256_load_pd(curr_cov+midx);
						__m256d xctx=_mm256_load_pd(ctxtemp+kx2);
						xctx=_mm256_mul_pd(xctx, yctx);
						xctx=_mm256_sub_pd(xctx, mcov);
						xctx=_mm256_mul_pd(xctx, mlr[kc]);
						xctx=_mm256_add_pd(xctx, mcov);
						_mm256_store_pd(curr_cov+midx, xctx);
#if 0
						curr_cov[midx+0]+=(ctxtemp[kx2+0]*ctxtemp[ky2]-curr_cov[midx+0])*lr;
						curr_cov[midx+1]+=(ctxtemp[kx2+1]*ctxtemp[ky2]-curr_cov[midx+1])*lr;
						curr_cov[midx+2]+=(ctxtemp[kx2+2]*ctxtemp[ky2]-curr_cov[midx+2])*lr;
						curr_cov[midx+3]+=(ctxtemp[kx2+3]*ctxtemp[ky2]-curr_cov[midx+3])*lr;
						for(k=0;k<4;++k)
						{
							if(abs(curr_cov[midx+k]-xctx.m256d_f64[k])>1e-6)
								printf("");
						}
#endif
					}
					for(;kx2<nparams;++kx2, ++midx)
						curr_cov[midx]+=(ctxtemp[kx2]*ctxtemp[ky2]-curr_cov[midx])*lr;
				}
#else
				for(ky2=0, midx=0;ky2<nparams;++ky2)
				{
					for(kx2=0;kx2<nparams;++kx2, ++midx)
						curr_cov[midx]+=(ctxtemp[kx2]*ctxtemp[ky2]-curr_cov[midx])*lr;
				}
#endif
#if 1
#ifdef ALLOW_AVX2
				{
					__m256d mval;
					double lval, lr_comp;

					lval=rows[0][kc+0]*lr;
					lr_comp=1-lr;
					mval=_mm256_set1_pd(lval);
					for(k=0;k<nparams-3;k+=4)
					{
						__m256d mvec=_mm256_load_pd(curr_vec+k);
						__m256d mtmp=_mm256_load_pd(ctxtemp+k);
						mtmp=_mm256_mul_pd(mval, mtmp);
						mvec=_mm256_mul_pd(mvec, mlr_comp[kc]);
						mvec=_mm256_add_pd(mvec, mtmp);
						_mm256_store_pd(curr_vec+k, mvec);
#if 0
						curr_vec[k+0]=lval*ctxtemp[k+0]+lr_comp*curr_vec[k+0];
						curr_vec[k+1]=lval*ctxtemp[k+1]+lr_comp*curr_vec[k+1];
						curr_vec[k+2]=lval*ctxtemp[k+2]+lr_comp*curr_vec[k+2];
						curr_vec[k+3]=lval*ctxtemp[k+3]+lr_comp*curr_vec[k+3];
#endif
					}
					for(;k<nparams;++k)
						curr_vec[k]=lval*ctxtemp[k]+lr_comp*curr_vec[k];
				}
#else
				double lval=rows[0][kc+0]*lr, lr_comp=1-lr;
				for(k=0;k<nparams;++k)
					curr_vec[k]=lval*ctxtemp[k]+lr_comp*curr_vec[k];
#endif
#else
				double lval=rows[0][kc+0];
				for(k=0;k<nparams;++k)
					curr_vec[k]+=(lval*ctxtemp[k]-curr_vec[k])*lr;
#endif
#if 0
				if(olsidx==period)/*OLS solver by Cholesky decomposition from paq8px*/
#endif
				{
					double sum;
					int success, n;

					success=1;
					n=nparams;
					memcpy(curr_cholesky, curr_cov, sizeof(double)*matsize[kc]);
					for(k=0;k<matsize[kc];k+=n+1)
						curr_cholesky[k]+=0.0075;
#if 1
					for(i=0;i<n;++i)
					{
						for(j=0;j<i;++j)
						{
							sum=curr_cholesky[i*n+j];
							for(k=0;k<j;++k)
								sum-=curr_cholesky[i*n+k]*curr_cholesky[j*n+k];
							curr_cholesky[i*n+j]=sum/curr_cholesky[j*n+j];
						}
						sum=curr_cholesky[i*n+i];
						for(k=0;k<i;++k)
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
						for(i=0;i<n;++i)
						{
							sum=curr_vec[i];
							for(j=0;j<i;++j)
								sum-=curr_cholesky[i*n+j]*curr_params[j];
							curr_params[i]=sum/curr_cholesky[i*n+i];
						}
						for(i=n-1;i>=0;--i)
						{
							sum=curr_params[i];
							for(j=i+1;j<n;++j)
								sum-=curr_cholesky[j*n+i]*curr_params[j];
							curr_params[i]=sum/curr_cholesky[i*n+i];
						}
#if 0
						for(k=0;k<n;++k)/*X*/
							curr_params[k]+=(ctxtemp[k]-curr_params[k])*0.9;
#endif
						++successcount;
					}
#endif
#if 0
					double *curr_vec2=vec2[kc];
					memcpy(curr_vec2, curr_vec, nparams*sizeof(double));
					success=solve_Mx_v(curr_cholesky, curr_vec2, nparams);
					if(success)
					{
						memcpy(curr_params, curr_vec2, nparams*sizeof(double));
#if 0
						for(k=0;k<n;++k)
							curr_params[k]+=(curr_vec2[k]-curr_params[k])*0.2;
#endif
					}
#endif
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
				for(k=0;k<nparams;++k)
					curr_ctx[k]+=8;
#endif
#endif
			}
			if(olsidx==olsnextpoint)
			{
				olsnextpoint<<=olsnextpoint<period;
				olsidx=0;
			}
#if 0
			olsidx&=-(olsidx<period);
			if(olsidx==period)
				olsidx=0;
#endif
			for(k=0;k<PADY;++k)
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
	for(kc=0;kc<4;++kc)
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