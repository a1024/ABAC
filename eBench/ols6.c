#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif

static int solve_Mx_v_cholesky(double *matrix, const double *vec, int n, double *solution)
{
	double sum;
	int success, i, j, k, kx, ky;

	success=1;
	for(i=0;i<n;++i)
	{
		for(j=0;j<i;++j)
		{
			sum=matrix[i*n+j];
			for(k=0;k<j;++k)
				sum-=matrix[i*n+k]*matrix[j*n+k];
			matrix[i*n+j]=sum/matrix[j*n+j];
		}
		sum=matrix[i*n+i];
		for(k=0;k<i;++k)
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
		for(ky=0;ky<n;++ky)
		{
			/*
			solution[0] = vec[0]/matrix[0][0]
			solution[1] = (vec[1] - solution[0]*matrix[1][0])/matrix[1][1]
			solution[2] = (vec[2] - solution[0]*matrix[2][0] - solution[1]*matrix[2][1])/matrix[2][2]
			...
			*/
			sum=vec[ky];
			for(kx=0;kx<ky;++kx)/*lower triangle downwards	_*\ |	*/
				sum-=matrix[ky*n+kx]*solution[kx];
			solution[ky]=sum/matrix[(n+1)*ky];
		}
		for(kx=n-1;kx>=0;--kx)
		{
			sum=solution[kx];
			for(ky=kx+1;ky<n;++ky)/*upper triangle upwards	_ \*|	*/
				sum-=matrix[ky*n+kx]*solution[ky];
			solution[kx]=sum/matrix[(n+1)*kx];
		}
	}
	return success;
}
#if 0
static int solve_Mx_v(double *matrix, double *vec, int size)
{
	int success, it, k, kx, ky;

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

typedef struct _V4Context
{
	int64_t *vbuf, *isum, *ptr;
	double *sum, *solution, *params, gain;
	int iw, nparams, cellsize, init, successcount, totalcalls, txid;
} V4Context;
static int v4_init(V4Context *ctx, int iw, int nparams, int depth)
{
	int cellsize=nparams*(nparams+1);
	int vbufsize=iw*cellsize;
	memset(ctx, 0, sizeof(*ctx));
	ctx->iw=iw;
	ctx->nparams=nparams;
	ctx->cellsize=cellsize;
	ctx->init=0;
	ctx->successcount=0;
	ctx->totalcalls=0;
	ctx->gain=1<<depth>>1;
	ctx->vbuf=(int64_t*)malloc(vbufsize*sizeof(int64_t));
	ctx->isum=(int64_t*)_mm_malloc(cellsize*sizeof(int64_t), sizeof(__m256i));
	ctx->sum=(double*)_mm_malloc(cellsize*sizeof(double), sizeof(__m256d));
	ctx->solution=(double*)malloc(nparams*sizeof(double));
	ctx->params=(double*)malloc((nparams+(intptr_t)1)*sizeof(double[256]));
	if(!ctx->vbuf||!ctx->isum||!ctx->sum||!ctx->solution||!ctx->params)
	{
		CRASH("Alloc error");
		return 1;
	}
	memset(ctx->vbuf, 0, vbufsize*sizeof(int64_t));
	memset(ctx->isum, 0, cellsize*sizeof(int64_t));
	memset(ctx->sum, 0, cellsize*sizeof(double));
	memset(ctx->solution, 0, nparams*sizeof(double));
	memset(ctx->params, 0, (nparams+(intptr_t)1)*sizeof(double[256]));

	return 0;
}
static void v4_finish(V4Context *ctx)
{
	free(ctx->vbuf);
	_mm_free(ctx->isum);
	_mm_free(ctx->sum);
	free(ctx->solution);
	free(ctx->params);
}
static void v4_addsample(V4Context *ctx, int64_t *dst, const int *nb)
{
	int kx, ky, msize, idx;

	msize=ctx->nparams*ctx->nparams;
	for(ky=0, idx=0;ky<ctx->nparams;++ky)
	{
		int val=nb[ky];
		for(kx=0;kx<ctx->nparams;++kx)
			dst[idx++]+=(int64_t)val*nb[kx];
		dst[msize+ky]+=(int64_t)val*nb[ctx->nparams];
	}
}
static void v4_subsample(V4Context *ctx, int64_t *dst, const int *nb)
{
	int kx, ky, msize, idx;

	msize=ctx->nparams*ctx->nparams;
	for(ky=0, idx=0;ky<ctx->nparams;++ky)
	{
		int val=nb[ky];
		for(kx=0;kx<ctx->nparams;++kx)
			dst[idx++]-=(int64_t)val*nb[kx];
		dst[msize+ky]-=(int64_t)val*nb[ctx->nparams];
	}
}
static void v4_addmv(V4Context *ctx, int64_t *dst, const int64_t *src)
{
	int k;

	for(k=0;k<ctx->cellsize;++k)
		dst[k]+=src[k];
}
static void v4_submv(V4Context *ctx, int64_t *dst, const int64_t *src)
{
	int k;

	for(k=0;k<ctx->cellsize;++k)
		dst[k]-=src[k];
}
static int v4_predict(V4Context *ctx, const int *nb, int W, int N, int NW, int NE, int NEEE, int txid)
{
	double *pWW, *pW, *params, fpred, reg;
	int k, msize, success;
	ALIGN(16) int pred[4];
	__m128i mpred, mN, mW, mNW, mNE, mNEEE, vmin, vmax;

	ctx->txid=txid;
	pWW	=ctx->params+(ctx->nparams+(intptr_t)1)*((ctx->successcount-(intptr_t)2)&3);
	pW	=ctx->params+(ctx->nparams+(intptr_t)1)*((ctx->successcount-(intptr_t)1)&3);
	params	=ctx->params+(ctx->nparams+(intptr_t)1)*((ctx->successcount+(intptr_t)0)&3);
	{
		msize=ctx->nparams*ctx->nparams;
		for(k=0;k<ctx->cellsize;++k)
			ctx->sum[k]=(double)ctx->isum[k];
		reg=ctx->gain*0.00005;
		for(k=0;k<ctx->nparams;++k)
			ctx->sum[(ctx->nparams+1)*k]+=reg;
		success=solve_Mx_v_cholesky(ctx->sum, ctx->sum+msize, ctx->nparams, ctx->solution);
		if(success)
		{
			if(ctx->init)
			{
				for(k=0;k<ctx->nparams;++k)
					params[k]=(2*pWW[k]+pW[k]+ctx->solution[k])*0.25;
			}
			else
			{
				memcpy(params, ctx->solution, ctx->nparams*sizeof(double));
				ctx->init=1;
			}
			++ctx->successcount;
		}
		++ctx->totalcalls;
	}
	fpred=0;
	for(k=0;k<ctx->nparams;++k)
		fpred+=params[k]*nb[k];

	mpred=_mm_cvtpd_epi32(_mm_set_sd(fpred));
	mN	=_mm_set_epi32(0, 0, 0, N);
	mW	=_mm_set_epi32(0, 0, 0, W);
	mNW	=_mm_set_epi32(0, 0, 0, NW);
	mNE	=_mm_set_epi32(0, 0, 0, NE);
	mNEEE	=_mm_set_epi32(0, 0, 0, NEEE);
	vmin=_mm_min_epi32(mN, mW);
	vmax=_mm_max_epi32(mN, mW);
	vmin=_mm_min_epi32(vmin, mNW);
	vmax=_mm_max_epi32(vmax, mNW);
	vmin=_mm_min_epi32(vmin, mNE);
	vmax=_mm_max_epi32(vmax, mNE);
	vmin=_mm_min_epi32(vmin, mNEEE);
	vmax=_mm_max_epi32(vmax, mNEEE);
	mpred=_mm_max_epi32(mpred, vmin);
	mpred=_mm_min_epi32(mpred, vmax);
	_mm_store_si128((__m128i*)pred, mpred);
	return *pred;
}
static void v4_newrow(V4Context *ctx)
{
	ctx->ptr=ctx->vbuf;
	memset(ctx->isum, 0, ctx->cellsize*sizeof(int64_t));
	ctx->init=0;
}


#define OLS6_REACH 1
#define OLS6_STEP 1
#define OLS6_PADX 16	/*64*/
#define OLS6_PADY 16	/*multiple of 4*/
#define OLS6_SAMPLEREACHX 8	/*60*/
#define OLS6_SAMPLEREACHY 8	/*12*/
#define OLS6_NPARAMS0 (2*(OLS6_REACH+1)*OLS6_REACH*3+0+1)
#define OLS6_NPARAMS1 (2*(OLS6_REACH+1)*OLS6_REACH*3+1+1)
#define OLS6_NPARAMS2 (2*(OLS6_REACH+1)*OLS6_REACH*3+2+1)
static void load_nb(int16_t **rows, int half, int *nb, int x, int y)
{
	int j=0;
	nb[j++]=half;
#if 0
	nb[j++]=rows[2-y][0+(x-2)*4];/*NNWW*/
	nb[j++]=rows[2-y][1+(x-2)*4];
	nb[j++]=rows[2-y][2+(x-2)*4];
	nb[j++]=rows[2-y][0+(x-1)*4];/*NNW*/
	nb[j++]=rows[2-y][1+(x-1)*4];
	nb[j++]=rows[2-y][2+(x-1)*4];
	nb[j++]=rows[2-y][0+(x+0)*4];/*NN*/
	nb[j++]=rows[2-y][1+(x+0)*4];
	nb[j++]=rows[2-y][2+(x+0)*4];
	nb[j++]=rows[2-y][0+(x+1)*4];/*NNE*/
	nb[j++]=rows[2-y][1+(x+1)*4];
	nb[j++]=rows[2-y][2+(x+1)*4];
	nb[j++]=rows[2-y][0+(x+2)*4];/*NNEE*/
	nb[j++]=rows[2-y][1+(x+2)*4];
	nb[j++]=rows[2-y][2+(x+2)*4];
	nb[j++]=rows[1-y][0+(x-2)*4];/*NWW*/
	nb[j++]=rows[1-y][1+(x-2)*4];
	nb[j++]=rows[1-y][2+(x-2)*4];
#endif
	nb[j++]=rows[1-y][0+(x-1)*4];/*NW*/
	nb[j++]=rows[1-y][1+(x-1)*4];
	nb[j++]=rows[1-y][2+(x-1)*4];
	nb[j++]=rows[1-y][0+(x+0)*4];/*N*/
	nb[j++]=rows[1-y][1+(x+0)*4];
	nb[j++]=rows[1-y][2+(x+0)*4];
	nb[j++]=rows[1-y][0+(x+1)*4];/*NE*/
	nb[j++]=rows[1-y][1+(x+1)*4];
	nb[j++]=rows[1-y][2+(x+1)*4];
#if 0
	nb[j++]=rows[1-y][0+(x+2)*4];/*NEE*/
	nb[j++]=rows[1-y][1+(x+2)*4];
	nb[j++]=rows[1-y][2+(x+2)*4];
	nb[j++]=rows[0-y][0+(x-2)*4];/*WW*/
	nb[j++]=rows[0-y][1+(x-2)*4];
	nb[j++]=rows[0-y][2+(x-2)*4];
#endif
	nb[j++]=rows[0-y][0+(x-1)*4];/*W*/
	nb[j++]=rows[0-y][1+(x-1)*4];
	nb[j++]=rows[0-y][2+(x-1)*4];
	nb[j++]=rows[0-y][0+(x+0)*4];/*curr*/
	nb[j++]=rows[0-y][1+(x+0)*4];
	nb[j++]=rows[0-y][2+(x+0)*4];
}
void pred_ols6(Image *src, int fwd)
{
	V4Context v4[3]={0};
	int half;
	int bufsize;
	int16_t *pixels;
	int ky, kx, k, kc;
	double t;

	t=time_sec();
	half=1<<src->depth[0]>>1;
	v4_init(v4+0, src->iw, OLS6_NPARAMS0, src->depth[0]);
	v4_init(v4+1, src->iw, OLS6_NPARAMS1, src->depth[1]);
	v4_init(v4+2, src->iw, OLS6_NPARAMS2, src->depth[2]);
	bufsize=(src->iw+OLS6_PADX*2)*(int)sizeof(int16_t[OLS6_PADY*4]);/*16 padded rows * 4 channels max*/
	pixels=(int16_t*)malloc(bufsize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	if(loud_transforms)
		DisableProcessWindowsGhosting();
	memset(pixels, 0, bufsize);
	for(ky=0;ky<src->ih;++ky)
	{
		__m256i rowstride=_mm256_set1_epi64x(sizeof(int16_t[4]));
		ALIGN(32) int16_t *rows[OLS6_PADY];

		for(k=0;k<(int)_countof(rows);++k)
			rows[k]=pixels+((src->iw+OLS6_PADX*2)*((ky-(size_t)k)%OLS6_PADY)+OLS6_PADX)*4;
		v4_newrow(v4+0);
		v4_newrow(v4+1);
		v4_newrow(v4+2);
		for(kx=0;kx<src->iw;++kx)
		{
			int future_bottom	[OLS6_NPARAMS2+1];
			int future_top		[OLS6_NPARAMS2+1];
			int present		[OLS6_NPARAMS2+1];
			int past		[OLS6_NPARAMS2+1];
			int nb			[OLS6_NPARAMS2+1];
			load_nb(rows, half, future_bottom	, OLS6_SAMPLEREACHX	, -1			);
			load_nb(rows, half, future_top		, OLS6_SAMPLEREACHX	, -OLS6_SAMPLEREACHY-1	);
			load_nb(rows, half, present		, -1			, 0			);
			load_nb(rows, half, past		, -OLS6_SAMPLEREACHX-1	, 0			);
			load_nb(rows, half, nb			, 0			, 0			);
			for(kc=0;kc<3;++kc)
			{
				V4Context *ctx;
				int pred, idx, curr;

				ctx=v4+kc;
				if(kx<src->iw-OLS6_SAMPLEREACHX)
				{
					if(ky>1)
						v4_addsample(ctx, ctx->ptr+OLS6_SAMPLEREACHX*ctx->cellsize, future_bottom);
					if(ky>OLS6_SAMPLEREACHY+1)
						v4_subsample(ctx, ctx->ptr+OLS6_SAMPLEREACHX*ctx->cellsize, future_top);
					v4_addmv(ctx, ctx->isum, ctx->ptr+OLS6_SAMPLEREACHX*ctx->cellsize);
				}
				if(kx>=OLS6_SAMPLEREACHX)
					v4_submv(ctx, ctx->isum, ctx->ptr-OLS6_SAMPLEREACHX*ctx->cellsize);
				if(kx>1)
					v4_addsample(ctx, ctx->isum, present);
				if(kx>OLS6_SAMPLEREACHX+1)
					v4_subsample(ctx, ctx->isum, past);
				if(ky>=OLS6_SAMPLEREACHY&&kx>=OLS6_SAMPLEREACHX&&kx<src->iw-OLS6_SAMPLEREACHX)
				{
					int txid=(((rows[1][kc+(0+0)*4]>>(src->depth[kc]-4))+8)&15)<<4|(((rows[0][kc+(0-1)*4]>>(src->depth[kc]-4))+8)&15);
					pred=v4_predict(ctx, nb, rows[0][kc-1*4], rows[1][kc+0*4], rows[1][kc-1*4], rows[1][kc+1*4], rows[1][kc+3*4], txid);
				}
				else
					MEDIAN3_32(pred, rows[1][kc+0*4], rows[0][kc-1*4], rows[1][kc+0*4]+rows[0][kc-1*4]-rows[1][kc-1*4]);

				idx=(src->iw*ky+kx)<<2;
				curr=0;
				if(fwd)
				{
					curr=((int32_t*)HDR_DATA(src))[idx|kc];
					pred=curr-pred;
					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];
				}
				else
				{
					pred+=((int32_t*)HDR_DATA(src))[idx|kc];
					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];
					curr=pred;
				}
				((int32_t*)HDR_DATA(src))[idx|kc]=pred;
				rows[0][kc]=curr;
				nb[_countof(nb)-3+kc]=curr;

				ctx->ptr+=ctx->cellsize;
			}
			for(k=0;k<(int)(sizeof(rows)/(sizeof(__m256i)));++k)
			{
				__m256i mr=_mm256_load_si256((__m256i*)rows+k);
				mr=_mm256_add_epi64(mr, rowstride);
				_mm256_store_si256((__m256i*)rows+k, mr);
			}
		}
		if(loud_transforms&&(!((ky+1)&63)||ky>=src->ih-1))
			set_window_title(
				"%d/%d %lf%%  %lf%% %lf%% %lf%%  %lf sec",
				ky+1, src->ih,
				100.*(ky+1)/src->ih,
				100.*v4[0].successcount/v4[0].totalcalls,
				100.*v4[1].successcount/v4[1].totalcalls,
				100.*v4[2].successcount/v4[2].totalcalls,
				time_sec()-t
			);
	}
	v4_finish(v4+0);
	v4_finish(v4+1);
	v4_finish(v4+2);
	free(pixels);
}