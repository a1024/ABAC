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
#define WG_NBITS 0

#define WG_DECAY_NUM	493
#define WG_DECAY_SH	9

#define WG_NPREDS	12	//multiple of 4
#if 1
#define WG_PREDLIST\
	WG_PRED(280,	N+eN/(kc2?2:4))\
	WG_PRED(280,	W+eW/(kc2?2:4))\
	WG_PRED(160,	3*(N-NN)+NNN)\
	WG_PRED(160,	3*(W-WW)+WWW)\
	WG_PRED(30,	4*(N+NNN)-6*NN+NNNW+NNNE-(NNWW+NNEE)/2)\
	WG_PRED(30,	4*(W+WWW)-6*WW+WWWW)\
	WG_PRED(130,	N+W-NW+((eN+eW-eNW+8)>>4))\
	WG_PRED(130,	N+NE-NNE+((eN+eNE-eNNE+8)>>4))\
	WG_PRED(180,	W+NE-N)\
	WG_PRED(150,	(W+NEEE-eW)/2)\
	WG_PRED(4,	NE+eNE/4)\
	WG_PRED(4,	NW+eNW/4)
//	WG_PRED(0, (W+2*NE-NNE+(eW+2*eNE-eNNE)/4)/2)
//	WG_PRED(0, N+W-NW)
//	WG_PRED(0, (-WW+3*W+N)/3)
//	WG_PRED(0, (4*(N+W)-WW-NN)/6)
//	WG_PRED(0, (10*W+5*(N-WW)-NN+WWW)/10)
//	WG_PRED(0, (15*(N+W)-6*(NN+WW)+NNN+WWW)/20)
//	WG_PRED(0, N+NE-NNE+((eN+eNE-eNNE)>>2))
#endif
#if 0
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
#endif
static void wg_init(double *weights, int kc)
{
	int j=0;
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
	WG_PREDLIST
#undef  WG_PRED
}
static int wg_predict(
	const double *weights,
	int **rows, const int stride, int kc2,
	int spred, const int *perrors, const int *Werrors, const int *Nerrors, const int *NEerrors, int *preds
)
{
	//double fpred=0, wsum=0;
	int j=0;
	ALIGN(16) int pred[4];
	int
		NNNWWWW	=rows[3][kc2-4*stride+0],
		NNNWWW	=rows[3][kc2-3*stride+0],
		NNNW	=rows[3][kc2-1*stride+0],
		NNN	=rows[3][kc2+0*stride+0],
		NNNE	=rows[3][kc2+1*stride+0],
		NNNEEE	=rows[3][kc2+3*stride+0],
		NNWWWW	=rows[2][kc2-4*stride+0],
		NNWW	=rows[2][kc2-2*stride+0],
		NNW	=rows[2][kc2-1*stride+0],
		NN	=rows[2][kc2+0*stride+0],
		NNE	=rows[2][kc2+1*stride+0],
		NNEE	=rows[2][kc2+2*stride+0],
		NNEEE	=rows[2][kc2+3*stride+0],
		NWW	=rows[1][kc2-2*stride+0],
		NW	=rows[1][kc2-1*stride+0],
		N	=rows[1][kc2+0*stride+0],
		NE	=rows[1][kc2+1*stride+0],
		NEE	=rows[1][kc2+2*stride+0],
		NEEE	=rows[1][kc2+3*stride+0],
		WWWW	=rows[0][kc2-4*stride+0],
		WWW	=rows[0][kc2-3*stride+0],
		WW	=rows[0][kc2-2*stride+0],
		W	=rows[0][kc2-1*stride+0],
		eNN	=rows[2][kc2+0*stride+1],
		eNNE	=rows[2][kc2+1*stride+1],
		eNW	=rows[1][kc2-1*stride+1],
		eN	=rows[1][kc2+0*stride+1],
		eNE	=rows[1][kc2+1*stride+1],
		eWW	=rows[0][kc2-2*stride+1],
		eW	=rows[0][kc2-1*stride+1];
	int sh=kc2?1:2;
	//int gy=abs(eN)+1, gx=abs(eW)+1;
	//int wgrad=(N*gy+W*gx)/(gy+gx);

	//int cgrad, cgrad2;
	//MEDIAN3_32(cgrad, N, W, N+W-NW);
	//MEDIAN3_32(cgrad2, N, NE, N+NE-NNE);
	
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST
#undef  WG_PRED

#if 1
	__m128i one=_mm_set1_epi32(1);
	__m128i me0=_mm_load_si128((__m128i*)perrors+0);
	__m128i me1=_mm_load_si128((__m128i*)perrors+1);
	__m128i me2=_mm_load_si128((__m128i*)perrors+2);
//	__m128i me3=_mm_load_si128((__m128i*)perrors+3);
	__m256d w0=_mm256_load_pd(weights+0*4);
	__m256d w1=_mm256_load_pd(weights+1*4);
	__m256d w2=_mm256_load_pd(weights+2*4);
//	__m256d w3=_mm256_load_pd(weights+3*4);
	//me0=_mm_srai_epi32(me0, 1);
	//me1=_mm_srai_epi32(me1, 1);
	//me2=_mm_srai_epi32(me2, 1);
	//me3=_mm_srai_epi32(me3, 1);
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)Werrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)Werrors+1));
	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)Werrors+2));
//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)Werrors+3));
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)Nerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)Nerrors+1));
	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)Nerrors+2));
//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)Nerrors+3));
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)NEerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)NEerrors+1));
	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)NEerrors+2));
//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)NEerrors+3));
	me0=_mm_add_epi32(me0, one);
	me1=_mm_add_epi32(me1, one);
	me2=_mm_add_epi32(me2, one);
//	me3=_mm_add_epi32(me3, one);
	__m256d de0=_mm256_cvtepi32_pd(me0);
	__m256d de1=_mm256_cvtepi32_pd(me1);
	__m256d de2=_mm256_cvtepi32_pd(me2);
//	__m256d de3=_mm256_cvtepi32_pd(me3);
	w0=_mm256_div_pd(w0, de0);
	w1=_mm256_div_pd(w1, de1);
	w2=_mm256_div_pd(w2, de2);
//	w3=_mm256_div_pd(w3, de3);
	de0=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+0));
	de1=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+1));
	de2=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+2));
//	de3=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+3));
	de0=_mm256_mul_pd(de0, w0);
	de1=_mm256_mul_pd(de1, w1);
	de2=_mm256_mul_pd(de2, w2);
//	de3=_mm256_mul_pd(de3, w3);
	w0=_mm256_add_pd(w0, w1);
	w0=_mm256_add_pd(w0, w2);
//	w0=_mm256_add_pd(w0, w3);
	de0=_mm256_add_pd(de0, de1);
	de0=_mm256_add_pd(de0, de2);
//	de0=_mm256_add_pd(de0, de3);
	
	//[num3 num2 num1 num0]
	//[den3 den2 den1 den0]
	//r = hadd(num, den) = [den3+den2 num3+num2 den1+den0 num1+num0]
	//lo=_mm256_extractf128_pd(r, 0)
	//hi=_mm256_extractf128_pd(r, 1)
	//hi+lo = [den3+den2+den1+den0 num3+num2+num1+num0]
	w0=_mm256_hadd_pd(de0, w0);
	__m128d dp=_mm_add_pd(_mm256_extractf128_pd(w0, 1), _mm256_extractf128_pd(w0, 0));
	dp=_mm_div_pd(dp, _mm_permute_pd(dp, 3));
	__m128i mp=_mm_cvtpd_epi32(dp);
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
	for(int k=0;k<WG_NPREDS;++k)
	{
		double weight=(double)weights[k]/(perrors[k]+1);
		fpred+=weight*preds[k];
		wsum+=weight;
	}
	fpred/=wsum;
#ifdef _MSC_VER
	pred=_cvt_dtoi_fast(fpred);
#else
	pred=(int)fpred;
#endif
	CLAMP3_32(pred, pred, N, W, NE);
#endif
	return pred[0];
}
static void wg_update(int curr, const int *preds, int *perrors, int *Werrors, int *currerrors, int *NEerrors)
{
	for(int k=0;k<WG_NPREDS;++k)
	{
		int e2=abs(curr-preds[k]);
		perrors[k]=(perrors[k]+e2)*WG_DECAY_NUM>>WG_DECAY_SH;
		e2<<=3;
		currerrors[k]=(2*Werrors[k]+e2+NEerrors[k])>>2;
		NEerrors[k]+=e2;
	}
}
void pred_wgrad4(Image *src, int fwd)
{
	ALIGN(32) double wg_weights[4][WG_NPREDS]={0};
	ALIGN(32) int wg_perrors[4][WG_NPREDS]={0}, wg_preds[WG_NPREDS]={0};
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

	int pesize=(src->iw+16)*(int)sizeof(int[4][4][WG_NPREDS]);//4 padded rows * 3 channels * WG_NPREDS
	int *ebuf=(int*)_mm_malloc(pesize, sizeof(__m256i));
	bufsize=(src->iw+16LL)*sizeof(int[4*4*2]);//4 padded rows * 4 channels max
	pixels=(int*)malloc(bufsize);
	if(!pixels||!ebuf)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(ebuf, 0, pesize);
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(int kc=0;kc<nch;++kc)
		wg_init(wg_weights[kc], kc);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		int *erows[]=
		{
			ebuf+((src->iw+16LL)*((ky-0LL)&3)+8LL)*4*WG_NPREDS,
			ebuf+((src->iw+16LL)*((ky-1LL)&3)+8LL)*4*WG_NPREDS,
			ebuf+((src->iw+16LL)*((ky-2LL)&3)+8LL)*4*WG_NPREDS,
			ebuf+((src->iw+16LL)*((ky-3LL)&3)+8LL)*4*WG_NPREDS,
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(int kc=0;kc<src->nch;++kc)
			{
				int pred;
				int
					kc2=kc<<1,
					kc3=kc*WG_NPREDS,
					NW	=rows[1][kc2-1*4*2],
					N	=rows[1][kc2+0*4*2],
					W	=rows[0][kc2-1*4*2],
					*eNW	=erows[1]+kc3-1*4*WG_NPREDS,
					*eN	=erows[1]+kc3+0*4*WG_NPREDS,
					*eNE	=erows[1]+kc3+1*4*WG_NPREDS,
					*eW	=erows[0]+kc3-1*4*WG_NPREDS,
					*ecurr	=erows[0]+kc3+0*4*WG_NPREDS;
				//if(ky==10&&kx==10)//
				//	printf("");
				//MEDIAN3_32(pred, N, W, N+W-NW);
				pred=wg_predict(wg_weights[kc], rows, 4*2, kc2, 0, wg_perrors[kc], eNW, eN, eNE, wg_preds);
				{
					int curr=src->data[idx+kc], pred0=pred;
					pred+=1<<WG_NBITS>>1;
					pred>>=WG_NBITS;
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=curr;

					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];

					src->data[idx+kc]=pred;
					rows[0][kc2+0]=(fwd?curr:pred)<<WG_NBITS;
					rows[0][kc2+1]=rows[0][kc2]-pred0;
				}
				wg_update(rows[0][kc2], wg_preds, wg_perrors[kc], eW, ecurr, eNE);
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
			erows[0]+=4*WG_NPREDS;
			erows[1]+=4*WG_NPREDS;
			erows[2]+=4*WG_NPREDS;
			erows[3]+=4*WG_NPREDS;
		}
	}
	free(pixels);
	_mm_free(ebuf);
}