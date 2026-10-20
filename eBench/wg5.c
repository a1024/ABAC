#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<process.h>
#include<immintrin.h>

/*WG:*/
#if 1
#define WG_PREDLIST\
	WG_PRED(240,	N)\
	WG_PRED(240,	W)\
	WG_PRED(180,	3*(N-NN)+NNN)\
	WG_PRED(180,	3*(W-WW)+WWW)\
	WG_PRED(140,	W+NE-N)\
	WG_PRED(160,	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)\
	WG_PRED(120,	N+W-NW)\
	WG_PRED(120,	N+NE-NNE)\

#endif
#if 0
#define WG_PREDLIST\
	WG_PRED(210,	N)\
	WG_PRED(210,	W)\
	WG_PRED(215,	3*(N-NN)+NNN)\
	WG_PRED(215,	3*(W-WW)+WWW)\
	WG_PRED(140,	W+NE-N)\
	WG_PRED(230,	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)\
	WG_PRED(120,	N+W-NW)\
	WG_PRED(140,	N+NE-NNE)\

#endif
enum
{
	WG_NBITS=0,		/*why 0 is best?*/
#define WG_PRED(W, E) +1
	WG_NPREDS=WG_PREDLIST	/*multiple of 4*/
#undef  WG_PRED
};
static void wg_init(double *weights, int kc)
{
	int j;

#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
	j=0;
	WG_PREDLIST
#undef  WG_PRED

	(void)kc;
}
static int wg_predict(
	const double *weights, int **rows, const int stride, int kc2,
	const int *perrors, const int *NWerrors, const int *Nerrors, const int *NEerrors, const int *NNEerrors, int *preds
)
{
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
		NNEEEE	=rows[2][kc2+4*stride+0],
		NWWW	=rows[1][kc2-3*stride+0],
		NWW	=rows[1][kc2-2*stride+0],
		NW	=rows[1][kc2-1*stride+0],
		N	=rows[1][kc2+0*stride+0],
		NE	=rows[1][kc2+1*stride+0],
		NEE	=rows[1][kc2+2*stride+0],
		NEEE	=rows[1][kc2+3*stride+0],
		NEEEE	=rows[1][kc2+4*stride+0],
		WWWWW	=rows[0][kc2-5*stride+0],
		WWWW	=rows[0][kc2-4*stride+0],
		WWW	=rows[0][kc2-3*stride+0],
		WW	=rows[0][kc2-2*stride+0],
		W	=rows[0][kc2-1*stride+0],
		eNNN	=rows[3][kc2+0*stride+1],
		eNN	=rows[2][kc2+0*stride+1],
		eNNE	=rows[2][kc2+1*stride+1],
		eNW	=rows[1][kc2-1*stride+1],
		eN	=rows[1][kc2+0*stride+1],
		eNE	=rows[1][kc2+1*stride+1],
		eNEE	=rows[1][kc2+2*stride+1],
		eNEEE	=rows[1][kc2+3*stride+1],
		eWWWW	=rows[0][kc2-4*stride+1],
		eWWW	=rows[0][kc2-3*stride+1],
		eWW	=rows[0][kc2-2*stride+1],
		eW	=rows[0][kc2-1*stride+1];
	__m256d w0, w1, de0, de1;
	__m128d dp;
	__m128i one, me0, me1, mp, mN, mW, mNE, vmin, vmax;
	ALIGN(16) int pred[4];
	int j;

#if 1
	(void)NNNWWWW;
	(void)NNNWWW;
	(void)NNNW;
	(void)NNN;
	(void)NNNE;
	(void)NNNEEE;
	(void)NNWWWW;
	(void)NNWW;
	(void)NNW;
	(void)NN;
	(void)NNE;
	(void)NNEE;
	(void)NNEEE;
	(void)NNEEEE;
	(void)NWWW;
	(void)NWW;
	(void)NW;
	(void)N;
	(void)NE;
	(void)NEE;
	(void)NEEE;
	(void)NEEEE;
	(void)WWWWW;
	(void)WWWW;
	(void)WWW;
	(void)WW;
	(void)W;
	(void)eNNN;
	(void)eNN;
	(void)eNNE;
	(void)eNW;
	(void)eN;
	(void)eNE;
	(void)eNEE;
	(void)eNEEE;
	(void)eWWWW;
	(void)eWWW;
	(void)eWW;
	(void)eW;
#endif
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	j=0;
	WG_PREDLIST
#undef  WG_PRED
#if 1
	one=_mm_set1_epi32(1);
	me0=_mm_load_si128((const __m128i*)perrors+0);
	me1=_mm_load_si128((const __m128i*)perrors+1);
	w0=_mm256_load_pd(weights+0*4);
	w1=_mm256_load_pd(weights+1*4);
	me0=_mm_add_epi32(me0, _mm_load_si128((const __m128i*)NWerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((const __m128i*)NWerrors+1));
	me0=_mm_add_epi32(me0, _mm_slli_epi32(_mm_load_si128((const __m128i*)Nerrors+0), 1));
	me1=_mm_add_epi32(me1, _mm_slli_epi32(_mm_load_si128((const __m128i*)Nerrors+1), 1));
	me0=_mm_add_epi32(me0, _mm_load_si128((const __m128i*)NEerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((const __m128i*)NEerrors+1));
	me0=_mm_add_epi32(me0, _mm_load_si128((const __m128i*)NNEerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((const __m128i*)NNEerrors+1));
	me0=_mm_add_epi32(me0, one);
	me1=_mm_add_epi32(me1, one);
	de0=_mm256_cvtepi32_pd(me0);
	de1=_mm256_cvtepi32_pd(me1);
	w0=_mm256_div_pd(w0, de0);
	w1=_mm256_div_pd(w1, de1);
	de0=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+0));
	de1=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+1));
	de0=_mm256_mul_pd(de0, w0);
	de1=_mm256_mul_pd(de1, w1);
	w0=_mm256_add_pd(w0, w1);
	de0=_mm256_add_pd(de0, de1);
	
	/*
	[num3 num2 num1 num0]
	[den3 den2 den1 den0]
	r = hadd(num, den) = [den3+den2 num3+num2 den1+den0 num1+num0]
	lo=_mm256_extractf128_pd(r, 0)
	hi=_mm256_extractf128_pd(r, 1)
	hi+lo = [den3+den2+den1+den0 num3+num2+num1+num0]
	*/
	w0=_mm256_hadd_pd(de0, w0);
	dp=_mm_add_pd(_mm256_extractf128_pd(w0, 1), _mm256_extractf128_pd(w0, 0));
	dp=_mm_div_pd(dp, _mm_permute_pd(dp, 3));
	mp=_mm_cvtpd_epi32(dp);
	mN	=_mm_set_epi32(0, 0, 0, N);
	mW	=_mm_set_epi32(0, 0, 0, W);
	mNE	=_mm_set_epi32(0, 0, 0, NE);
	vmin=_mm_min_epi32(mN, mW);
	vmax=_mm_max_epi32(mN, mW);
	vmin=_mm_min_epi32(vmin, mNE);
	vmax=_mm_max_epi32(vmax, mNE);
	mp=_mm_max_epi32(mp, vmin);
	mp=_mm_min_epi32(mp, vmax);
	_mm_store_si128((__m128i*)pred, mp);
#else
	{
		int k;

		for(k=0;k<WG_NPREDS;++k)
		{
			double weight=(double)weights[k]/(perrors[k]+1);
			fpred+=weight*preds[k];
			wsum+=weight;
		}
	}
	fpred/=wsum;
	pred=(int)TRUNC64(fpred);
	CLAMP3_32(pred, pred, N, W, NE);
#endif
	return pred[0];
}
static void wg_update(int curr, int kc, const int *preds, int *perrors, int *Werrors, int *currerrors, int *NEerrors)
{
	static const int factors[]={97, 99, 99};
	int factor, k, e2;

	factor=factors[kc];
	for(k=0;k<WG_NPREDS;++k)
	{
		e2=abs(curr-preds[k])<<1;
		perrors[k]=(perrors[k]+e2)*factor>>7;
		currerrors[k]=(2*Werrors[k]+e2+NEerrors[k])>>2;
		NEerrors[k]+=e2;
	}
}
void pred_wgrad5(Image *src, int fwd)
{
	ALIGN(32) double wg_weights[4][WG_NPREDS]={0};
	ALIGN(32) int wg_perrors[4][WG_NPREDS]={0}, wg_preds[WG_NPREDS]={0};
	int amin[4], amax[4], nch, bufsize, *pixels, *rows[4], invdist, pesize, *ebuf, *erows[4], ky, kx, kc, idx;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	invdist=((1<<16)+g_dist-1)/g_dist;
	pesize=(src->iw+16)*(int)sizeof(int[4][4][WG_NPREDS]);/*4 padded rows * 3 channels * WG_NPREDS*/
	ebuf=(int*)_mm_malloc(pesize, sizeof(__m256i));
	bufsize=(src->iw+16)*(int)sizeof(int[4*4*2]);/*4 padded rows * 4 channels max*/
	pixels=(int*)malloc(bufsize);
	if(!pixels||!ebuf)
	{
		CRASH("Alloc error");
		return;
	}
	memset(ebuf, 0, pesize);
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(kc=0;kc<nch;++kc)
		wg_init(wg_weights[kc], kc);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4*2;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4*2;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4*2;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4*2;
		erows[0]=ebuf+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4*WG_NPREDS;
		erows[1]=ebuf+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4*WG_NPREDS;
		erows[2]=ebuf+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4*WG_NPREDS;
		erows[3]=ebuf+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4*WG_NPREDS;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<src->nch;++kc)
			{
				int
					kc2=kc<<1,
					kc3=kc*WG_NPREDS,
					*eNNE	=erows[2]+kc3+1*4*WG_NPREDS,
					*eNW	=erows[1]+kc3-1*4*WG_NPREDS,
					*eN	=erows[1]+kc3+0*4*WG_NPREDS,
					*eNE	=erows[1]+kc3+1*4*WG_NPREDS,
					*eW	=erows[0]+kc3-1*4*WG_NPREDS,
					*ecurr	=erows[0]+kc3+0*4*WG_NPREDS;
				int pred, curr, e;

				pred=wg_predict(wg_weights[0], rows, 4*2, kc2, wg_perrors[kc], eNW, eN, eNE, eNNE, wg_preds);
				curr=((int32_t*)HDR_DATA(src))[idx+kc];
				if(g_dist>1)
				{
					if(fwd)
					{
						curr-=(int)pred;
						curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));/*curr/=g_dist*/
						((int32_t*)HDR_DATA(src))[idx+kc]=curr;

						curr=g_dist*curr+(int)pred;
						CLAMP2(curr, amin[kc], amax[kc]);
					}
					else
					{
						curr=g_dist*curr+(int)pred;
						CLAMP2(curr, amin[kc], amax[kc]);

						((int32_t*)HDR_DATA(src))[idx+kc]=curr;
					}
				}
				else
				{
					if(fwd)
					{
						e=curr-pred;
						e<<=32-src->depth[kc];
						e>>=32-src->depth[kc];
						((int32_t*)HDR_DATA(src))[idx+kc]=e;
					}
					else
					{
						curr+=pred;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						((int32_t*)HDR_DATA(src))[idx+kc]=curr;
					}
				}
				rows[0][kc2+0]=curr;
				rows[0][kc2+1]=curr-pred;
				wg_update(rows[0][kc2], kc, wg_preds, wg_perrors[kc], eW, ecurr, eNE);
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

#define WM_PREDLIST\
	WM_PRED(210,	N)\
	WM_PRED(210,	W)\
	WM_PRED(215,	3*(N-NN)+NNN)\
	WM_PRED(215,	3*(W-WW)+WWW)\
	WM_PRED(140,	W+NE-N)\
	WM_PRED(230,	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)\
	WM_PRED(120,	N+W-NW)\
	WM_PRED(140,	N+NE-NNE)\

enum
{
#define WM_PRED(W, E) +1
	WM_NPREDS=WM_PREDLIST	/*multiple of 4*/
#undef  WM_PRED
};
void pred_wmix(Image *src, int fwd)
{
	ALIGN(32) int alphas[4][WM_NPREDS]={0}, preds[WM_NPREDS]={0};
	int amin[4], amax[4], nch, bufsize, *pixels, *rows[4], ky, kx, kc, idx;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	FILLMEM((int*)alphas, 0x1000/WM_NPREDS, sizeof(alphas), sizeof(int));
	bufsize=(src->iw+16)*(int)sizeof(int[4*4*2]);/*4 padded rows * 4 channels max*/
	pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4*2;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4*2;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4*2;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4*2;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<src->nch;++kc)
			{
				int
					kc2	=kc<<1,
					NNN	=rows[3][kc2+0*4*2+0],
					NN	=rows[2][kc2+0*4*2+0],
					NNE	=rows[2][kc2+1*4*2+0],
					NW	=rows[1][kc2-1*4*2+0],
					N	=rows[1][kc2+0*4*2+0],
					NE	=rows[1][kc2+1*4*2+0],
					NEE	=rows[1][kc2+2*4*2+0],
					NEEE	=rows[1][kc2+3*4*2+0],
					NEEEE	=rows[1][kc2+4*4*2+0],
					WWWW	=rows[0][kc2-4*4*2+0],
					WWW	=rows[0][kc2-3*4*2+0],
					WW	=rows[0][kc2-2*4*2+0],
					W	=rows[0][kc2-1*4*2+0];
				int pred, j, k, vmin, vmax, curr, val, esum;

#define WM_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
				j=0;
				WM_PREDLIST;
#undef  WM_PRED
				pred=1<<12>>1;
				for(k=1;k<WM_NPREDS;++k)
					pred+=preds[k]*alphas[kc][k];
				pred>>=12;
				vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				CLAMP2(pred, vmin, vmax);

				curr=((int32_t*)HDR_DATA(src))[idx+kc];
				if(fwd)
				{
					val=(curr-(int)pred+g_dist/2)/g_dist;
					curr=g_dist*val+(int)pred;
				}
				else
				{
					val=g_dist*curr+(int)pred;
					curr=val;
					CLAMP2(val, amin[kc], amax[kc]);
				}
				((int32_t*)HDR_DATA(src))[idx+kc]=keyboard[KEY_ALT]?pred:val;
				rows[0][kc2+0]=curr;
				rows[0][kc2+1]=curr-pred;
				esum=0;
				for(k=0;k<WM_NPREDS;++k)
					esum+=preds[k]=abs(curr-preds[k]);
				for(k=0;k<WM_NPREDS;++k)
					alphas[kc][k]+=(((esum-WM_NPREDS*preds[k])<<2)-alphas[kc][k]+0x1000/WM_NPREDS+(1<<8>>1))>>8;
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	free(pixels);
}
