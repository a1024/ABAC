#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
//#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;

//WG:
#define WG_DECAY_NUM	493
#define WG_DECAY_SH	9

#define WG_NPREDS	16
#define WG_PREDLIST\
	WG_PRED(25, spred)\
	WG_PRED(50, wgrad)\
	WG_PRED(50, 3*N-W-NE)\
	WG_PRED(30, 2*NE-W)\
	WG_PRED(132, N+W-NW)\
	WG_PRED(100, N+W-NW+((eN+eW-eNW+16)>>5))\
	WG_PRED(150, N+eN)\
	WG_PRED(100, N)\
	WG_PRED(150, W+eW)\
	WG_PRED(100, W)\
	WG_PRED(165, W+NE-N)\
	WG_PRED(220, N+NE-NNE)\
	WG_PRED(25, N+NE-NNE+((eN+eNE-eNNE)>>2))\
	WG_PRED(165, 3*(N-NN)+NNN)\
	WG_PRED(165, 3*(W-WW)+WWW)\
	WG_PRED(150, (W+NEEE)/2)
static void wg_init(int *weights)
{
	int j=0;
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
	WG_PREDLIST
#undef  WG_PRED
}
static int wg_predict(
	const int *weights,
	int **rows, const int stride, int kc2,
	int spred, const int *perrors, int *preds
)
{
	double pred2=0, wsum=0;
	int j=0;
	ALIGN(16) int pred[4];
	int
		NNNWWWW	=rows[3][kc2-4*stride+0],
		NNWWWW	=rows[2][kc2-4*stride+0],
		NNN	=rows[3][kc2+0*stride+0],
		NN	=rows[2][kc2+0*stride+0],
		NNE	=rows[2][kc2+1*stride+0],
		NNEEE	=rows[2][kc2+3*stride+0],
		NW	=rows[1][kc2-1*stride+0],
		N	=rows[1][kc2+0*stride+0],
		NE	=rows[1][kc2+1*stride+0],
		NEE	=rows[1][kc2+2*stride+0],
		NEEE	=rows[1][kc2+3*stride+0],
		WWW	=rows[0][kc2-3*stride+0],
		WW	=rows[0][kc2-2*stride+0],
		W	=rows[0][kc2-1*stride+0],
		eNNE	=rows[2][kc2+1*stride+1],
		eNW	=rows[1][kc2-1*stride+1],
		eN	=rows[1][kc2+0*stride+1],
		eNE	=rows[1][kc2+1*stride+1],
		eW	=rows[0][kc2-1*stride+1];

	int gy=abs(eN)+1, gx=abs(eW)+1;
	int wgrad=(N*gy+W*gx)/(gy+gx);

#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST
#undef  WG_PRED
	
	for(int k=0;k<WG_NPREDS;++k)
	{
		double weight=(double)weights[k]/(perrors[k]+1);
		pred2+=weight*preds[k];
		wsum+=weight;
	}
	pred2/=wsum;
#if 1
	__m128i mp=_mm_cvtpd_epi32(_mm_set_pd(0, pred2));
	__m128i mN	=_mm_set_epi32(0, 0, 0, N);
	__m128i mW	=_mm_set_epi32(0, 0, 0, W);
	__m128i mNE	=_mm_set_epi32(0, 0, 0, NE);
	__m128i vmin=_mm_min_epi32(mN, mW);
	__m128i vmax=_mm_max_epi32(mN, mW);
	vmin=_mm_min_epi32(vmin, mNE);
	vmax=_mm_max_epi32(vmax, mNE);
	mp=_mm_max_epi32(mp, vmin);
	mp=_mm_min_epi32(mp, vmax);
	_mm_store_si128((__m128i*)pred, mp);
#else
#ifdef _MSC_VER
	pred=_cvt_dtoi_fast(pred2);
#else
	pred=(int)pred2;
#endif
	CLAMP3_32(pred, pred, N, W, NE);
#endif
	return pred[0];
}
static void wg_update(int curr, const int *preds, int *perrors, int *weights)
{
	for(int k=0;k<WG_NPREDS;++k)
	{
		int e2=abs(curr-preds[k]);
		perrors[k]=(perrors[k]+e2)*WG_DECAY_NUM>>WG_DECAY_SH;
	}
}
void pred_wgrad4(Image *src, int fwd)
{
	ALIGN(32) int wg_weights[4][WG_NPREDS]={0}, wg_perrors[4][WG_NPREDS]={0}, wg_preds[WG_NPREDS]={0};
	int nch;
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int fwdmask=-fwd;
	int bufsize;
	int *pixels;

	bufsize=(src->iw+16LL)*sizeof(int[4*4*2]);//4 padded rows * 4 channels max
	pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(int kc=0;kc<nch;++kc)
		wg_init(wg_weights[kc]);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(int kc=0;kc<src->nch;++kc)
			{
				int pred;
				int
					kc2=kc<<1,
					NW	=rows[1][kc2-1*4*2],
					N	=rows[1][kc2+0*4*2],
					W	=rows[0][kc2-1*4*2];
				//if(ky==10&&kx==10)//
				//	printf("");
				MEDIAN3_32(pred, N, W, N+W-NW);
				pred=wg_predict(wg_weights[kc], rows, 4*2, kc2, pred, wg_perrors[kc], wg_preds);
				{
					int curr=src->data[idx+kc], pred0=pred;
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=curr;

					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];

					src->data[idx+kc]=pred;
					rows[0][kc2+0]=fwd?curr:pred;
					rows[0][kc2+1]=rows[0][kc2]-pred0;
				}
				wg_update(rows[0][kc2], wg_preds, wg_perrors[kc], wg_weights[kc]);
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	free(pixels);
}