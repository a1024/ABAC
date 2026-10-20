#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>


	#define OLS5_POPUP


#define ALLOCASSERT(C)\
	do\
		if(C)\
		{\
			CRASH("Alloc error");\
			return;\
		}\
	while(0)

typedef struct OLS5Struct
{
	int n;
	double *nb, *vec, *cov, *cholesky, *params;
} OLS5;
static void ols5_init(OLS5 *p, int n)
{
	int msize=n*n;
	memset(p, 0, sizeof(*p));
	p->n=n;
	p->nb=(double*)_mm_malloc(sizeof(double)*n, sizeof(__m256d));
	p->vec=(double*)_mm_malloc(sizeof(double)*n, sizeof(__m256d));
	p->cov=(double*)_mm_malloc(sizeof(double)*msize, sizeof(__m256d));
	p->cholesky=(double*)_mm_malloc(sizeof(double)*msize, sizeof(__m256d));
	p->params=(double*)_mm_malloc(sizeof(double)*n, sizeof(__m256d));
	ALLOCASSERT(!p->nb||!p->vec||!p->cov||!p->cholesky||!p->params);
	memset(p->nb, 0, sizeof(double)*n);
	memset(p->vec, 0, sizeof(double)*n);
	memset(p->cov, 0, sizeof(double)*msize);
	memset(p->cholesky, 0, sizeof(double)*msize);
	memset(p->params, 0, sizeof(double)*n);
}
static void ols5_clear(OLS5 *p)
{
	_mm_free(p->nb);
	_mm_free(p->vec);
	_mm_free(p->cov);
	_mm_free(p->cholesky);
	_mm_free(p->params);
}
static double ols5_predict(OLS5 *p, double *nb)
{
	double fpred;
	int k;

	fpred=0;
	for(k=0;k<p->n;++k)
		fpred+=p->params[k]*nb[k];
	return fpred;
}
static int ols5_update(OLS5 *p, double *nb, double curr, double lr)/*OLS solver by Cholesky decomposition from paq8px*/
{
	double lval, lr_comp, sum;
	int ky2, kx2, midx, i, j, k, success, msize;

	for(ky2=0, midx=0;ky2<p->n;++ky2)
	{
		for(kx2=0;kx2<p->n;++kx2, ++midx)
			p->cov[midx]+=(nb[kx2]*nb[ky2]-p->cov[midx])*lr;
	}
	lval=curr*lr;
	lr_comp=1-lr;
	for(k=0;k<p->n;++k)
		p->vec[k]=lval*nb[k]+lr_comp*p->vec[k];
	success=1;
	msize=p->n*p->n;
	memcpy(p->cholesky, p->cov, sizeof(double)*msize);
	for(k=0;k<msize;k+=p->n+1)
		p->cholesky[k]+=0.0075;
	for(i=0;i<p->n;++i)
	{
		for(j=0;j<i;++j)
		{
			sum=p->cholesky[i*p->n+j];
			for(k=0;k<j;++k)
				sum-=p->cholesky[i*p->n+k]*p->cholesky[j*p->n+k];
			p->cholesky[i*p->n+j]=sum/p->cholesky[j*p->n+j];
		}
		sum=p->cholesky[i*p->n+i];
		for(k=0;k<i;++k)
			sum-=p->cholesky[i*p->n+k]*p->cholesky[i*p->n+k];
		if(sum<=1e-8)
		{
			success=0;
			break;
		}
		p->cholesky[i*p->n+i]=sqrt(sum);
	}
	if(success)
	{
		for(i=0;i<p->n;++i)
		{
			sum=p->vec[i];
			for(j=0;j<i;++j)
				sum-=p->cholesky[i*p->n+j]*p->params[j];
			p->params[i]=sum/p->cholesky[i*p->n+i];
		}
		for(i=p->n-1;i>=0;--i)
		{
			sum=p->params[i];
			for(j=i+1;j<p->n;++j)
				sum-=p->cholesky[j*p->n+i]*p->params[j];
			p->params[i]=sum/p->cholesky[i*p->n+i];
		}
	}
	return success;
}
static int ols5_clamp4(double fpred, double a, double b, double c, double d)
{
	ALIGN(16) int pred[4];
	__m128i ipred;
	__m128d mpred=_mm_set_sd(fpred);
	__m128d ma=_mm_set_sd(a);
	__m128d mb=_mm_set_sd(b);
	__m128d mc=_mm_set_sd(c);
	__m128d md=_mm_set_sd(d);
	__m128d vmin=_mm_min_pd(ma, mb);
	__m128d vmax=_mm_max_pd(ma, mb);

	vmin=_mm_min_pd(vmin, mc);
	vmax=_mm_max_pd(vmax, mc);
	vmin=_mm_min_pd(vmin, md);
	vmax=_mm_max_pd(vmax, md);
	mpred=_mm_max_pd(mpred, vmin);
	mpred=_mm_min_pd(mpred, vmax);
	ipred=_mm_cvtpd_epi32(mpred);
	_mm_store_si128((__m128i*)pred, ipred);
	return pred[0];
#if 0
	{
		double vmin=MINVAR(a, b), vmax=MAXVAR(a, b);
		int pred;

		UPDATE_MIN(vmin, c);
		UPDATE_MAX(vmax, c);
		UPDATE_MIN(vmin, d);
		UPDATE_MAX(vmax, d);
		fpred=CLAMP(vmin, fpred, vmax);
		pred=(int)round(fpred);
		return pred;
	}
#endif
}
#define LOAD(BUF, C, X, Y) ((uint32_t)(ky+(Y))<(uint32_t)src->ih&&(uint32_t)(kx+(X))<(uint32_t)src->iw?(double)BUF[(src->iw*(ky+(Y))+kx+(X))<<2|C]:0)
void pred_ols5(Image *src, int fwd)
{
	const double lr=0.005;
	double t_start;
	Image *dst;
	OLS5 *ols, *curr_ols;
	int nlevels[4], half[4], fwdmask, successcount, kc, kx, ky, ks, j, k;
	const int *pixels;

	nlevels[0]=1<<src->depth[0];
	nlevels[1]=1<<src->depth[1];
	nlevels[2]=1<<src->depth[2];
	nlevels[3]=1<<src->depth[3];
	half[0]=nlevels[0]>>1;
	half[1]=nlevels[1]>>1;
	half[2]=nlevels[2]>>1;
	half[3]=nlevels[3]>>1;
	fwdmask=-fwd;
	successcount=0;

	t_start=time_sec();
	if(loud_transforms)
		DisableProcessWindowsGhosting();
	dst=0;
	image_copy(&dst, src);
	ols=(OLS5*)malloc(sizeof(OLS5[16]));
	ALLOCASSERT(!dst||!ols);
	for(kc=0;kc<4;++kc)
	{
		ols5_init(ols+kc+ 0, 12);
		ols5_init(ols+kc+ 4, 16);

		ols5_init(ols+kc+ 8, 6+6);
#if 0
		ols5_init(ols+kc+ 8, 12+12);
#endif

		ols5_init(ols+kc+12, 8+8);
#if 0
		ols5_init(ols+kc+12, 16+12+12);
#endif
	}

	pixels=fwd?((int32_t*)HDR_DATA(dst)):((int32_t*)HDR_DATA(src));
	curr_ols=ols+0;
	for(ky=0;ky<src->ih-1;ky+=2)
	{
		for(kx=0;kx<src->iw-1;kx+=2)
		{
			for(kc=0;kc<4;++kc)
			{
				double nb[12], fpred, vmin, vmax;
				int pred, idx;

				if(!src->depth[kc])
					continue;
				j=0;
				nb[j++]=LOAD(pixels, kc, -4, -4);
				nb[j++]=LOAD(pixels, kc, -2, -4);
				nb[j++]=LOAD(pixels, kc,  0, -4);
				nb[j++]=LOAD(pixels, kc,  2, -4);
				nb[j++]=LOAD(pixels, kc,  4, -4);
				nb[j++]=LOAD(pixels, kc, -4, -2);
				nb[j++]=LOAD(pixels, kc, -2, -2);
				nb[j++]=LOAD(pixels, kc,  0, -2);/*7:NN*/
				nb[j++]=LOAD(pixels, kc,  2, -2);
				nb[j++]=LOAD(pixels, kc,  4, -2);
				nb[j++]=LOAD(pixels, kc, -4,  0);
				nb[j++]=LOAD(pixels, kc, -2,  0);/*11:WW*/
				fpred=ols5_predict(curr_ols+kc, nb);
				vmin=MINVAR(nb[7], nb[11]);
				vmax=MAXVAR(nb[7], nb[11]);
				fpred=CLAMP(vmin, fpred, vmax);
				pred=(int)round(fpred);

				idx=(src->iw*ky+kx)<<2|kc;
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=((int32_t*)HDR_DATA(src))[idx];
				pred+=half[kc];
				pred&=nlevels[kc]-1;
				pred-=half[kc];
				((int32_t*)HDR_DATA(src))[idx]=pred;
				successcount+=ols5_update(curr_ols+kc, nb, pixels[idx], lr);
			}
		}
	}
	if(loud_transforms)
		set_window_title("OLS5 1/4  %f sec", time_sec()-t_start);
	curr_ols=ols+4;
	for(ky=1;ky<src->ih-1;ky+=2)
	{
		for(kx=1;kx<src->iw-1;kx+=2)
		{
			for(kc=0;kc<4;++kc)
			{
				double nb[16], fpred;
				int pred, idx;

				if(!src->depth[kc])
					continue;
				j=0;
				nb[j++]=LOAD(pixels, kc, -3, -3);
				nb[j++]=LOAD(pixels, kc, -1, -3);
				nb[j++]=LOAD(pixels, kc,  1, -3);
				nb[j++]=LOAD(pixels, kc,  3, -3);
				nb[j++]=LOAD(pixels, kc, -3, -1);
				nb[j++]=LOAD(pixels, kc, -1, -1);/* 5*/
				nb[j++]=LOAD(pixels, kc,  1, -1);/* 6*/
				nb[j++]=LOAD(pixels, kc,  3, -1);
				nb[j++]=LOAD(pixels, kc, -3,  1);
				nb[j++]=LOAD(pixels, kc, -1,  1);/* 9*/
				nb[j++]=LOAD(pixels, kc,  1,  1);/*10*/
				nb[j++]=LOAD(pixels, kc,  3,  1);
				nb[j++]=LOAD(pixels, kc, -3,  3);
				nb[j++]=LOAD(pixels, kc, -1,  3);
				nb[j++]=LOAD(pixels, kc,  1,  3);
				nb[j++]=LOAD(pixels, kc,  3,  3);
				fpred=ols5_predict(curr_ols+kc, nb);
				pred=ols5_clamp4(fpred, nb[5], nb[6], nb[9], nb[10]);
				idx=(src->iw*ky+kx)<<2|kc;
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=((int32_t*)HDR_DATA(src))[idx];
				pred+=half[kc];
				pred&=nlevels[kc]-1;
				pred-=half[kc];
				((int32_t*)HDR_DATA(src))[idx]=pred;
				successcount+=ols5_update(curr_ols+kc, nb, pixels[idx], lr);
			}
		}
	}
	if(loud_transforms)
		set_window_title("OLS5 2/4  %f sec", time_sec()-t_start);
	curr_ols=ols+8;
	for(ky=0;ky<src->ih-1;ky+=2)
	{
		for(kx=1;kx<src->iw-1;kx+=2)
		{
			for(kc=0;kc<4;++kc)
			{
				double nb[12], fpred;
				int pred, idx;

				if(!src->depth[kc])
					continue;
				j=0;
				nb[j++]=LOAD(pixels, kc, -1, -2);
				nb[j++]=LOAD(pixels, kc,  1, -2);
				nb[j++]=LOAD(pixels, kc, -1,  0);/* 2*/
				nb[j++]=LOAD(pixels, kc,  1,  0);/* 3*/
				nb[j++]=LOAD(pixels, kc, -1,  2);
				nb[j++]=LOAD(pixels, kc,  1,  2);

				nb[j++]=LOAD(pixels, kc, -2, -1);
				nb[j++]=LOAD(pixels, kc,  0, -1);/* 7*/
				nb[j++]=LOAD(pixels, kc,  2, -1);
				nb[j++]=LOAD(pixels, kc, -2,  1);
				nb[j++]=LOAD(pixels, kc,  0,  1);/*10*/
				nb[j++]=LOAD(pixels, kc,  2,  1);
#if 0
				nb[j++]=LOAD(pixels, kc, -3, -2);
				nb[j++]=LOAD(pixels, kc, -1, -2);
				nb[j++]=LOAD(pixels, kc,  1, -2);
				nb[j++]=LOAD(pixels, kc,  3, -2);
				nb[j++]=LOAD(pixels, kc, -3,  0);
				nb[j++]=LOAD(pixels, kc, -1,  0);/* 5*/
				nb[j++]=LOAD(pixels, kc,  1,  0);/* 6*/
				nb[j++]=LOAD(pixels, kc,  3,  0);
				nb[j++]=LOAD(pixels, kc, -3,  2);
				nb[j++]=LOAD(pixels, kc, -1,  2);
				nb[j++]=LOAD(pixels, kc,  1,  2);
				nb[j++]=LOAD(pixels, kc,  3,  2);
				
				nb[j++]=LOAD(pixels, kc, -2, -3);
				nb[j++]=LOAD(pixels, kc,  0, -3);
				nb[j++]=LOAD(pixels, kc,  2, -3);
				nb[j++]=LOAD(pixels, kc, -2, -1);
				nb[j++]=LOAD(pixels, kc,  0, -1);/*16*/
				nb[j++]=LOAD(pixels, kc,  2, -1);
				nb[j++]=LOAD(pixels, kc, -2,  1);
				nb[j++]=LOAD(pixels, kc,  0,  1);/*19*/
				nb[j++]=LOAD(pixels, kc,  2,  1);
				nb[j++]=LOAD(pixels, kc, -2,  3);
				nb[j++]=LOAD(pixels, kc,  0,  3);
				nb[j++]=LOAD(pixels, kc,  2,  3);
#endif
				fpred=ols5_predict(curr_ols+kc, nb);
				pred=ols5_clamp4(fpred, nb[2], nb[3], nb[7], nb[10]);
				idx=(src->iw*ky+kx)<<2|kc;
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=((int32_t*)HDR_DATA(src))[idx];
				pred+=half[kc];
				pred&=nlevels[kc]-1;
				pred-=half[kc];
				((int32_t*)HDR_DATA(src))[idx]=pred;
				successcount+=ols5_update(curr_ols+kc, nb, pixels[idx], lr);
			}
		}
	}
	if(loud_transforms)
		set_window_title("OLS5 3/4  %f sec", time_sec()-t_start);
	curr_ols=ols+12;
	for(ky=1;ky<src->ih-1;ky+=2)
	{
		for(kx=0;kx<src->iw-1;kx+=2)
		{
			for(kc=0;kc<4;++kc)
			{
				double nb[16], fpred;
				int pred, idx;

				if(!src->depth[kc])
					continue;
				j=0;
#if 1
				nb[j++]=LOAD(pixels, kc, -1, -1);
				nb[j++]=LOAD(pixels, kc,  0, -1);/* 1*/
				nb[j++]=LOAD(pixels, kc,  1, -1);
				nb[j++]=LOAD(pixels, kc, -1,  0);/* 3*/
				nb[j++]=LOAD(pixels, kc,  1,  0);/* 4*/
				nb[j++]=LOAD(pixels, kc, -1,  1);
				nb[j++]=LOAD(pixels, kc,  0,  1);/* 6*/
				nb[j++]=LOAD(pixels, kc,  1,  1);
				
				nb[j++]=LOAD(pixels, kc, -1, -2);
				nb[j++]=LOAD(pixels, kc,  1, -2);
				nb[j++]=LOAD(pixels, kc, -2, -1);
				nb[j++]=LOAD(pixels, kc,  2, -1);
				nb[j++]=LOAD(pixels, kc, -2,  1);
				nb[j++]=LOAD(pixels, kc,  2,  1);
				nb[j++]=LOAD(pixels, kc, -1,  2);
				nb[j++]=LOAD(pixels, kc,  1,  2);
#endif
#if 0
				nb[j++]=LOAD(pixels, kc, -2, -3);
				nb[j++]=LOAD(pixels, kc,  0, -3);
				nb[j++]=LOAD(pixels, kc,  2, -3);
				nb[j++]=LOAD(pixels, kc, -2, -1);
				nb[j++]=LOAD(pixels, kc,  0, -1);/* 4*/
				nb[j++]=LOAD(pixels, kc,  2, -1);
				nb[j++]=LOAD(pixels, kc, -2,  1);
				nb[j++]=LOAD(pixels, kc,  0,  1);/* 7*/
				nb[j++]=LOAD(pixels, kc,  2,  1);
				nb[j++]=LOAD(pixels, kc, -2,  3);
				nb[j++]=LOAD(pixels, kc,  0,  3);
				nb[j++]=LOAD(pixels, kc,  2,  3);
				
				nb[j++]=LOAD(pixels, kc, -3, -2);
				nb[j++]=LOAD(pixels, kc, -1, -2);
				nb[j++]=LOAD(pixels, kc,  1, -2);
				nb[j++]=LOAD(pixels, kc,  3, -2);
				nb[j++]=LOAD(pixels, kc, -3,  0);
				nb[j++]=LOAD(pixels, kc, -1,  0);/*17*/
				nb[j++]=LOAD(pixels, kc,  1,  0);/*18*/
				nb[j++]=LOAD(pixels, kc,  3,  0);
				nb[j++]=LOAD(pixels, kc, -3,  2);
				nb[j++]=LOAD(pixels, kc, -1,  2);
				nb[j++]=LOAD(pixels, kc,  1,  2);
				nb[j++]=LOAD(pixels, kc,  3,  2);
				
				nb[j++]=LOAD(pixels, kc, -3, -3);
				nb[j++]=LOAD(pixels, kc, -1, -3);
				nb[j++]=LOAD(pixels, kc,  1, -3);
				nb[j++]=LOAD(pixels, kc,  3, -3);
				nb[j++]=LOAD(pixels, kc, -3, -1);
				nb[j++]=LOAD(pixels, kc, -1, -1);
				nb[j++]=LOAD(pixels, kc,  1, -1);
				nb[j++]=LOAD(pixels, kc,  3, -1);
				nb[j++]=LOAD(pixels, kc, -3,  1);
				nb[j++]=LOAD(pixels, kc, -1,  1);
				nb[j++]=LOAD(pixels, kc,  1,  1);
				nb[j++]=LOAD(pixels, kc,  3,  1);
				nb[j++]=LOAD(pixels, kc, -3,  3);
				nb[j++]=LOAD(pixels, kc, -1,  3);
				nb[j++]=LOAD(pixels, kc,  1,  3);
				nb[j++]=LOAD(pixels, kc,  3,  3);
#endif
				fpred=ols5_predict(curr_ols+kc, nb);
				pred=ols5_clamp4(fpred, nb[1], nb[3], nb[4], nb[6]);
				idx=(src->iw*ky+kx)<<2|kc;
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=((int32_t*)HDR_DATA(src))[idx];
				pred+=half[kc];
				pred&=nlevels[kc]-1;
				pred-=half[kc];
				((int32_t*)HDR_DATA(src))[idx]=pred;
				successcount+=ols5_update(curr_ols+kc, nb, pixels[idx], lr);
			}
		}
	}
	if(loud_transforms)
		set_window_title("OLS5 4/4  %f sec", time_sec()-t_start);
#ifdef OLS5_POPUP
	if(loud_transforms)
	{
		enum
		{
			BUFSIZE=0x10000
		};
		int maxdepth, *hist;
		int qres;
		char *str;
		double csizes[16]={0};
		int printed;
		double total;

		maxdepth=MAXVAR(src->depth[0], src->depth[1]);
		UPDATE_MAX(maxdepth, src->depth[2]);
		UPDATE_MAX(maxdepth, src->depth[3]);
		hist=(int*)malloc(sizeof(int[16])<<maxdepth);
		ALLOCASSERT(!hist);
		memset(hist, 0, sizeof(int[16])<<maxdepth);
		for(ky=0;ky<src->ih-1;ky+=2)
		{
			for(kx=0;kx<src->iw-1;kx+=2)
			{
				for(kc=0;kc<4;++kc)
				{
					int v0=(((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+0)<<2|kc]+half[kc])&(nlevels[kc]-1);
					int v1=(((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+1)<<2|kc]+half[kc])&(nlevels[kc]-1);
					int v2=(((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|kc]+half[kc])&(nlevels[kc]-1);
					int v3=(((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+0)<<2|kc]+half[kc])&(nlevels[kc]-1);
					++hist[(kc+ 0)<<maxdepth|v0];
					++hist[(kc+ 4)<<maxdepth|v1];
					++hist[(kc+ 8)<<maxdepth|v2];
					++hist[(kc+12)<<maxdepth|v3];
				}
			}
		}
		qres=(src->iw>>1)*(src->ih>>1);
		for(kc=0;kc<16;++kc)
		{
			double psum=0;
			for(ks=0;ks<nlevels[kc&3];++ks)
			{
				int freq=hist[kc<<maxdepth|ks];
				if(freq)
				{
					double p=(double)freq/qres;
					csizes[kc]-=freq*log(p);
					psum+=p;
				}
			}
			if(fabs(psum-1)>1e-6)
				crash(__FILE__, __LINE__, 1, "Wrong probability denominator:  C%d  P=%f", kc, psum);
			csizes[kc]*=1./(8*M_LN2);
		}
		str=(char*)malloc(BUFSIZE);
		ALLOCASSERT(!str);
		printed=0;
		total=0;
		printed+=snprintf(str+printed, (size_t)BUFSIZE-printed-1, "U     %13d\n"
			, qres
		);
		for(kc=0;kc<16;++kc)
		{
			printed+=snprintf(str+printed, (size_t)BUFSIZE-printed-1, "C%dQ%d  %16.2f  %10.6f%%  %10.6f\n"
				, kc&3, kc>>2, csizes[kc], 100.*csizes[kc]/qres, qres/csizes[kc]
			);
			total+=csizes[kc];
		}
		qres*=src->nch;
		for(kc=0;kc<16-3;kc+=4)
		{
			double csize=csizes[kc+0]+csizes[kc+1]+csizes[kc+2]+csizes[kc+3];
			printed+=snprintf(str+printed, (size_t)BUFSIZE-printed-1, "AvQ%d  %16.2f  %10.6f%%  %10.6f\n"
				, kc>>2, csize, 100.*csize/qres, qres/csize
			);
		}
		{
			int usize=src->iw*src->ih*(src->src_depth[0]+src->src_depth[1]+src->src_depth[2]+src->src_depth[3])>>3;
			printed+=snprintf(str+printed, (size_t)BUFSIZE-printed-1, "T     %16.2f/%d  %10.6f%%  %10.6f\n"
				, total, usize, 100.*total/usize, usize/total
			);
		}
		copy_to_clipboard(str, printed);
		messagebox(MBOX_OK, "Copied to clipboard", "%s", str);
		free(str);
		free(hist);
	}
#endif
	free(dst);
	for(k=0;k<16;++k)
		ols5_clear(ols+k);
	(void)successcount;
}