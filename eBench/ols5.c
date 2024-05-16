#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
//#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>
static const char file[]=__FILE__;


	#define OLS5_POPUP


#define ALLOCASSERT(C)\
	do\
		if(C)\
		{\
			LOG_ERROR("Alloc error");\
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
	double fpred=0;
	for(int k=0;k<p->n;++k)
		fpred+=p->params[k]*nb[k];
	return fpred;
}
static int ols5_update(OLS5 *p, double *nb, double curr, double lr)//OLS solver by Cholesky decomposition from paq8px
{
	for(int ky2=0, midx=0;ky2<p->n;++ky2)
	{
		for(int kx2=0;kx2<p->n;++kx2, ++midx)
			p->cov[midx]+=(nb[kx2]*nb[ky2]-p->cov[midx])*lr;
	}
	double lval=curr*lr, lr_comp=1-lr;
	for(int k=0;k<p->n;++k)
		p->vec[k]=lval*nb[k]+lr_comp*p->vec[k];
	int success=1;
	double sum;
	int msize=p->n*p->n;
	memcpy(p->cholesky, p->cov, sizeof(double)*msize);
	for(int k=0;k<msize;k+=p->n+1)
		p->cholesky[k]+=0.0075;
	for(int i=0;i<p->n;++i)
	{
		for(int j=0;j<i;++j)
		{
			sum=p->cholesky[i*p->n+j];
			for(int k=0;k<j;++k)
				sum-=p->cholesky[i*p->n+k]*p->cholesky[j*p->n+k];
			p->cholesky[i*p->n+j]=sum/p->cholesky[j*p->n+j];
		}
		sum=p->cholesky[i*p->n+i];
		for(int k=0;k<i;++k)
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
		for(int i=0;i<p->n;++i)
		{
			sum=p->vec[i];
			for(int j=0;j<i;++j)
				sum-=p->cholesky[i*p->n+j]*p->params[j];
			p->params[i]=sum/p->cholesky[i*p->n+i];
		}
		for(int i=p->n-1;i>=0;--i)
		{
			sum=p->params[i];
			for(int j=i+1;j<p->n;++j)
				sum-=p->cholesky[j*p->n+i]*p->params[j];
			p->params[i]=sum/p->cholesky[i*p->n+i];
		}
	}
	return success;
}
static int ols5_clamp4(double fpred, double a, double b, double c, double d)
{
	double vmin=MINVAR(a, b), vmax=MAXVAR(a, b);
	UPDATE_MIN(vmin, c);
	UPDATE_MAX(vmax, c);
	UPDATE_MIN(vmin, d);
	UPDATE_MAX(vmax, d);
	fpred=CLAMP(vmin, fpred, vmax);
	int pred=(int)round(fpred);
	return pred;
}
#define LOAD(BUF, C, X, Y) ((unsigned)(ky+(Y))<(unsigned)src->ih&&(unsigned)(kx+(X))<(unsigned)src->iw?(double)BUF[(src->iw*(ky+(Y))+kx+(X))<<2|C]:0)
void pred_ols5(Image *src, int fwd)
{
	double t_start=time_sec();
	if(loud_transforms)
		DisableProcessWindowsGhosting();
	Image *dst=0;
	image_copy(&dst, src);
	OLS5 *ols=(OLS5*)malloc(sizeof(OLS5[16]));
	ALLOCASSERT(!dst||!ols);
	for(int kc=0;kc<4;++kc)
	{
		ols5_init(ols+kc+ 0, 12);
		ols5_init(ols+kc+ 4, 16);

		ols5_init(ols+kc+ 8, 6+6);
		//ols5_init(ols+kc+ 8, 12+12);

		ols5_init(ols+kc+12, 8+8);
		//ols5_init(ols+kc+12, 16+12+12);
	}

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
	double lr=0.005;
	int fwdmask=-fwd;
	int successcount=0;
	const int *pixels=fwd?dst->data:src->data;
	OLS5 *curr_ols=ols+0;
	for(int ky=0;ky<src->ih-1;ky+=2)
	{
		for(int kx=0;kx<src->iw-1;kx+=2)
		{
			for(int kc=0;kc<4;++kc)
			{
				if(!src->depth[kc])
					continue;
				double nb[]=
				{
					LOAD(pixels, kc, -4, -4),
					LOAD(pixels, kc, -2, -4),
					LOAD(pixels, kc,  0, -4),
					LOAD(pixels, kc,  2, -4),
					LOAD(pixels, kc,  4, -4),
					LOAD(pixels, kc, -4, -2),
					LOAD(pixels, kc, -2, -2),
					LOAD(pixels, kc,  0, -2),//7:NN
					LOAD(pixels, kc,  2, -2),
					LOAD(pixels, kc,  4, -2),
					LOAD(pixels, kc, -4,  0),
					LOAD(pixels, kc, -2,  0),//11:WW
				};
				double fpred=ols5_predict(curr_ols+kc, nb);
				double vmin=MINVAR(nb[7], nb[11]), vmax=MAXVAR(nb[7], nb[11]);
				fpred=CLAMP(vmin, fpred, vmax);
				int pred=(int)round(fpred);

				int idx=(src->iw*ky+kx)<<2|kc;
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=src->data[idx];
				pred+=half[kc];
				pred&=nlevels[kc]-1;
				pred-=half[kc];
				src->data[idx]=pred;
				
				successcount+=ols5_update(curr_ols+kc, nb, pixels[idx], lr);
			}
		}
	}
	if(loud_transforms)
		set_window_title("OLS5 1/4  %lf sec", time_sec()-t_start);
	curr_ols=ols+4;
	for(int ky=1;ky<src->ih-1;ky+=2)
	{
		for(int kx=1;kx<src->iw-1;kx+=2)
		{
			for(int kc=0;kc<4;++kc)
			{
				if(!src->depth[kc])
					continue;
				double nb[]=
				{
					LOAD(pixels, kc, -3, -3),
					LOAD(pixels, kc, -1, -3),
					LOAD(pixels, kc,  1, -3),
					LOAD(pixels, kc,  3, -3),
					LOAD(pixels, kc, -3, -1),
					LOAD(pixels, kc, -1, -1),//5
					LOAD(pixels, kc,  1, -1),//6
					LOAD(pixels, kc,  3, -1),
					LOAD(pixels, kc, -3,  1),
					LOAD(pixels, kc, -1,  1),//9
					LOAD(pixels, kc,  1,  1),//10
					LOAD(pixels, kc,  3,  1),
					LOAD(pixels, kc, -3,  3),
					LOAD(pixels, kc, -1,  3),
					LOAD(pixels, kc,  1,  3),
					LOAD(pixels, kc,  3,  3),
				};
				double fpred=ols5_predict(curr_ols+kc, nb);
				int pred=ols5_clamp4(fpred, nb[5], nb[6], nb[9], nb[10]);

				int idx=(src->iw*ky+kx)<<2|kc;
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=src->data[idx];
				pred+=half[kc];
				pred&=nlevels[kc]-1;
				pred-=half[kc];
				src->data[idx]=pred;
				
				successcount+=ols5_update(curr_ols+kc, nb, pixels[idx], lr);
			}
		}
	}
	if(loud_transforms)
		set_window_title("OLS5 2/4  %lf sec", time_sec()-t_start);
	curr_ols=ols+8;
	for(int ky=0;ky<src->ih-1;ky+=2)
	{
		for(int kx=1;kx<src->iw-1;kx+=2)
		{
			for(int kc=0;kc<4;++kc)
			{
				if(!src->depth[kc])
					continue;
				double nb[]=
				{
					LOAD(pixels, kc, -1, -2),
					LOAD(pixels, kc,  1, -2),
					LOAD(pixels, kc, -1,  0),//2
					LOAD(pixels, kc,  1,  0),//3
					LOAD(pixels, kc, -1,  2),
					LOAD(pixels, kc,  1,  2),

					LOAD(pixels, kc, -2, -1),
					LOAD(pixels, kc,  0, -1),//7
					LOAD(pixels, kc,  2, -1),
					LOAD(pixels, kc, -2,  1),
					LOAD(pixels, kc,  0,  1),//10
					LOAD(pixels, kc,  2,  1),
#if 0
					LOAD(pixels, kc, -3, -2),
					LOAD(pixels, kc, -1, -2),
					LOAD(pixels, kc,  1, -2),
					LOAD(pixels, kc,  3, -2),
					LOAD(pixels, kc, -3,  0),
					LOAD(pixels, kc, -1,  0),//5
					LOAD(pixels, kc,  1,  0),//6
					LOAD(pixels, kc,  3,  0),
					LOAD(pixels, kc, -3,  2),
					LOAD(pixels, kc, -1,  2),
					LOAD(pixels, kc,  1,  2),
					LOAD(pixels, kc,  3,  2),
					
					LOAD(pixels, kc, -2, -3),
					LOAD(pixels, kc,  0, -3),
					LOAD(pixels, kc,  2, -3),
					LOAD(pixels, kc, -2, -1),
					LOAD(pixels, kc,  0, -1),//16
					LOAD(pixels, kc,  2, -1),
					LOAD(pixels, kc, -2,  1),
					LOAD(pixels, kc,  0,  1),//19
					LOAD(pixels, kc,  2,  1),
					LOAD(pixels, kc, -2,  3),
					LOAD(pixels, kc,  0,  3),
					LOAD(pixels, kc,  2,  3),
#endif
				};
				double fpred=ols5_predict(curr_ols+kc, nb);
				int pred=ols5_clamp4(fpred, nb[2], nb[3], nb[7], nb[10]);
				//int pred=ols5_clamp4(fpred, nb[5], nb[6], nb[16], nb[19]);

				int idx=(src->iw*ky+kx)<<2|kc;
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=src->data[idx];
				pred+=half[kc];
				pred&=nlevels[kc]-1;
				pred-=half[kc];
				src->data[idx]=pred;
				
				successcount+=ols5_update(curr_ols+kc, nb, pixels[idx], lr);
			}
		}
	}
	if(loud_transforms)
		set_window_title("OLS5 3/4  %lf sec", time_sec()-t_start);
	curr_ols=ols+12;
	for(int ky=1;ky<src->ih-1;ky+=2)
	{
		for(int kx=0;kx<src->iw-1;kx+=2)
		{
			for(int kc=0;kc<4;++kc)
			{
				if(!src->depth[kc])
					continue;
				double nb[]=
				{
#if 1
					LOAD(pixels, kc, -1, -1),
					LOAD(pixels, kc,  0, -1),//1
					LOAD(pixels, kc,  1, -1),
					LOAD(pixels, kc, -1,  0),//3
					LOAD(pixels, kc,  1,  0),//4
					LOAD(pixels, kc, -1,  1),
					LOAD(pixels, kc,  0,  1),//6
					LOAD(pixels, kc,  1,  1),
					
					LOAD(pixels, kc, -1, -2),
					LOAD(pixels, kc,  1, -2),
					LOAD(pixels, kc, -2, -1),
					LOAD(pixels, kc,  2, -1),
					LOAD(pixels, kc, -2,  1),
					LOAD(pixels, kc,  2,  1),
					LOAD(pixels, kc, -1,  2),
					LOAD(pixels, kc,  1,  2),
#endif
#if 0
					LOAD(pixels, kc, -2, -3),
					LOAD(pixels, kc,  0, -3),
					LOAD(pixels, kc,  2, -3),
					LOAD(pixels, kc, -2, -1),
					LOAD(pixels, kc,  0, -1),//4
					LOAD(pixels, kc,  2, -1),
					LOAD(pixels, kc, -2,  1),
					LOAD(pixels, kc,  0,  1),//7
					LOAD(pixels, kc,  2,  1),
					LOAD(pixels, kc, -2,  3),
					LOAD(pixels, kc,  0,  3),
					LOAD(pixels, kc,  2,  3),
					
					LOAD(pixels, kc, -3, -2),
					LOAD(pixels, kc, -1, -2),
					LOAD(pixels, kc,  1, -2),
					LOAD(pixels, kc,  3, -2),
					LOAD(pixels, kc, -3,  0),
					LOAD(pixels, kc, -1,  0),//17
					LOAD(pixels, kc,  1,  0),//18
					LOAD(pixels, kc,  3,  0),
					LOAD(pixels, kc, -3,  2),
					LOAD(pixels, kc, -1,  2),
					LOAD(pixels, kc,  1,  2),
					LOAD(pixels, kc,  3,  2),
					
					LOAD(pixels, kc, -3, -3),
					LOAD(pixels, kc, -1, -3),
					LOAD(pixels, kc,  1, -3),
					LOAD(pixels, kc,  3, -3),
					LOAD(pixels, kc, -3, -1),
					LOAD(pixels, kc, -1, -1),
					LOAD(pixels, kc,  1, -1),
					LOAD(pixels, kc,  3, -1),
					LOAD(pixels, kc, -3,  1),
					LOAD(pixels, kc, -1,  1),
					LOAD(pixels, kc,  1,  1),
					LOAD(pixels, kc,  3,  1),
					LOAD(pixels, kc, -3,  3),
					LOAD(pixels, kc, -1,  3),
					LOAD(pixels, kc,  1,  3),
					LOAD(pixels, kc,  3,  3),
#endif
				};
				double fpred=ols5_predict(curr_ols+kc, nb);
				int pred=ols5_clamp4(fpred, nb[1], nb[3], nb[4], nb[6]);
				//int pred=ols5_clamp4(fpred, nb[4], nb[7], nb[17], nb[18]);

				int idx=(src->iw*ky+kx)<<2|kc;
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=src->data[idx];
				pred+=half[kc];
				pred&=nlevels[kc]-1;
				pred-=half[kc];
				src->data[idx]=pred;
				
				successcount+=ols5_update(curr_ols+kc, nb, pixels[idx], lr);
			}
		}
	}
	if(loud_transforms)
		set_window_title("OLS5 4/4  %lf sec", time_sec()-t_start);
#ifdef OLS5_POPUP
	if(loud_transforms)
	{
		int maxdepth=MAXVAR(src->depth[0], src->depth[1]);
		UPDATE_MAX(maxdepth, src->depth[2]);
		UPDATE_MAX(maxdepth, src->depth[3]);
		int *hist=(int*)malloc(sizeof(int[16])<<maxdepth);
		ALLOCASSERT(!hist);
		memset(hist, 0, sizeof(int[16])<<maxdepth);
		for(int ky=0;ky<src->ih-1;ky+=2)
		{
			for(int kx=0;kx<src->iw-1;kx+=2)
			{
				for(int kc=0;kc<4;++kc)
				{
					int v0=(src->data[(src->iw*(ky+0)+kx+0)<<2|kc]+half[kc])&(nlevels[kc]-1);
					int v1=(src->data[(src->iw*(ky+1)+kx+1)<<2|kc]+half[kc])&(nlevels[kc]-1);
					int v2=(src->data[(src->iw*(ky+0)+kx+1)<<2|kc]+half[kc])&(nlevels[kc]-1);
					int v3=(src->data[(src->iw*(ky+1)+kx+0)<<2|kc]+half[kc])&(nlevels[kc]-1);
					++hist[(kc+ 0)<<maxdepth|v0];
					++hist[(kc+ 4)<<maxdepth|v1];
					++hist[(kc+ 8)<<maxdepth|v2];
					++hist[(kc+12)<<maxdepth|v3];
				}
			}
		}
		double csizes[16]={0};
		int qres=(src->iw>>1)*(src->ih>>1);
		for(int kc=0;kc<16;++kc)
		{
			double psum=0;
			for(int ks=0;ks<nlevels[kc&3];++ks)
			{
				int freq=hist[kc<<maxdepth|ks];
				if(freq)
				{
					double p=(double)freq/qres;
					csizes[kc]-=freq*log2(p);
					psum+=p;
				}
			}
			if(fabs(psum-1)>1e-6)
				LOG_ERROR("Wrong probability denominator:  C%d  P=%lf", kc, psum);
			csizes[kc]/=8;
		}
		char *str=(char*)malloc(0x10000);
		ALLOCASSERT(!str);
		int printed=0;
		double total=0;
		printed+=snprintf(str+printed, 0x10000LL-printed-1, "U     %13d\n", qres);
		for(int kc=0;kc<16;++kc)
		{
			printed+=snprintf(str+printed, 0x10000LL-printed-1, "C%dQ%d  %16.2lf  %10.6lf%%  %10.6lf\n",
				kc&3, kc>>2, csizes[kc], 100.*csizes[kc]/qres, qres/csizes[kc]
			);
			total+=csizes[kc];
		}
		qres*=src->nch;
		for(int kc=0;kc<16-3;kc+=4)
		{
			double csize=csizes[kc+0]+csizes[kc+1]+csizes[kc+2]+csizes[kc+3];
			printed+=snprintf(str+printed, 0x10000LL-printed-1, "AvQ%d  %16.2lf  %10.6lf%%  %10.6lf\n",
				kc>>2, csize, 100.*csize/qres, qres/csize
			);
		}
		int usize=src->iw*src->ih*(src->src_depth[0]+src->src_depth[1]+src->src_depth[2]+src->src_depth[3])>>3;
		printed+=snprintf(str+printed, 0x10000LL-printed-1, "T     %16.2lf/%d  %10.6lf%%  %10.6lf\n", total, usize, 100.*total/usize, usize/total);
		copy_to_clipboard(str, printed);
		messagebox(MBOX_OK, "Copied to clipboard", "%s", str);
		free(str);
		free(hist);
	}
#endif
	free(dst);
	for(int k=0;k<16;++k)
		ols5_clear(ols+k);
}