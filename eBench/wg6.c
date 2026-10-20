#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<process.h>
#include<immintrin.h>


#if 1
#	define USE_DIVTABLE
#endif


#if 1
#define PREDLIST\
	PRED(210, N+(23*eN-2*(eNN+eNW)+9*eW)/110)\
	PRED(210, W+(23*eW-2*(eWW+eNW)+9*eN)/110)\
	PRED(215, 3*(N-NN)+NNN+(eN/2+eNN/2+eNNN)/3)\
	PRED(215, 3*(W-WW)+WWW+(eW/2+eWW/2+eWWW)/3)\
	PRED(140, W+NE-N+((-13*eN)>>4)+eW*13/32-(eW>>7))\
	PRED(230, (WWW+NN-2*NW+NEE+NEEE+NEEEE+(N-W)/6+(NNN-NE-(5*(eN+eW)+eWW))/2)/3)\
	PRED(120, N+W-NW+(eN+eW-eNW)/3)\
	PRED(140, N+NE-NNE+((eN+eNE+eNNE+4)>>3))\

#endif
#if 0
#define PREDLIST\
	PRED(210, N-eN/16)\
	PRED(210, W-eW/16)\
	PRED(215, 3*(N-NN)+NNN)\
	PRED(215, 3*(W-WW)+WWW)\
	PRED(140, W+NE-N)\
	PRED(230, (WWWW+WWW+NNN-2*NW+NEE+NEEE+NEEEE)/4)\
	PRED(120, N+W-NW)\
	PRED(140, N+NE-NNE)\

#endif
#if 0
#define PREDLIST\
	PRED(210, N)\
	PRED(210, W)\
	PRED(215, 3*(N-NN)+NNN)\
	PRED(215, 3*(W-WW)+WWW)\
	PRED(140, W+NE-N)\
	PRED(230, (WWWW+WWW+NNN-2*NW+NEE+NEEE+NEEEE)/4)\
	PRED(120, N+W-NW)\
	PRED(140, N+NE-NNE)\

#endif

enum
{
#define PRED(WEIGHT, EXPR) +1
	NPREDS=PREDLIST
#undef  PRED
};
static const int wg_weights[]=
{
#define PRED(WEIGHT, EXPR) WEIGHT,
	PREDLIST
#undef  PRED
};
static int wg_predict(
	const int *errors,
	const int *coeffs, const uint16_t *cN, const uint16_t *ccurr,
	const int16_t *NNNptr, const int16_t *NNptr, const int16_t *Nptr, const int16_t *currptr,
	int *subpreds, const int *mixer, const int *sse,
	int *preds, int depth, int kc
)
{
#ifdef USE_DIVTABLE
	enum
	{
		NUMBITS=15,
		DENBITS=7,
		DIVLUTSIZE=(1<<DENBITS)
	};
	static int divlookup[DIVLUTSIZE]={0};
#endif
	int
		NNNN	=currptr[0+0*4*2],
		NNNWWWW	=NNNptr	[0-4*4*2],
		NNNWWW	=NNNptr	[0-3*4*2],
		NNNW	=NNNptr	[0-1*4*2],
		NNN	=NNNptr	[0+0*4*2],
		NNNE	=NNNptr	[0+1*4*2],
		NNNEEE	=NNNptr	[0+3*4*2],
		NNWWWW	=NNptr	[0-4*4*2],
		NNWW	=NNptr	[0-2*4*2],
		NNW	=NNptr	[0-1*4*2],
		NN	=NNptr	[0+0*4*2],
		NNE	=NNptr	[0+1*4*2],
		NNEE	=NNptr	[0+2*4*2],
		NNEEE	=NNptr	[0+3*4*2],
		NNEEEE	=NNptr	[0+4*4*2],
		NWWW	=Nptr	[0-3*4*2],
		NWW	=Nptr	[0-2*4*2],
		NW	=Nptr	[0-1*4*2],
		N	=Nptr	[0+0*4*2],
		NE	=Nptr	[0+1*4*2],
		NEE	=Nptr	[0+2*4*2],
		NEEE	=Nptr	[0+3*4*2],
		NEEEE	=Nptr	[0+4*4*2],
		WWWWW	=currptr[0-5*4*2],
		WWWW	=currptr[0-4*4*2],
		WWW	=currptr[0-3*4*2],
		WW	=currptr[0-2*4*2],
		W	=currptr[0-1*4*2],
		eNNN	=NNNptr	[4+0*4*2],
		eNN	=NNptr	[4+0*4*2],
		eNNE	=NNptr	[4+1*4*2],
		eNW	=Nptr	[4-1*4*2],
		eN	=Nptr	[4+0*4*2],
		eNE	=Nptr	[4+1*4*2],
		eNEE	=Nptr	[4+2*4*2],
		eNEEE	=Nptr	[4+3*4*2],
		eWWWW	=currptr[4-4*4*2],
		eWWW	=currptr[4-3*4*2],
		eWW	=currptr[4-2*4*2],
		eW	=currptr[4-1*4*2];
	uint32_t coeff[NPREDS], wsum=0;
	int j, k, kp, ipred, pred, sh[NPREDS], vmin, vmax;

#if 1
	(void)NNNN;
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
#define PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	j=0;
	PREDLIST;
#undef  PRED
#ifdef USE_DIVTABLE
	if(!*divlookup)
	{
		for(k=0;k<DIVLUTSIZE;++k)
			divlookup[k]=(1<<NUMBITS)/(k+1);
	}
	ipred=0;
	pred=0;
	for(kp=0;kp<NPREDS;++kp)
	{
		/*
		|			1	1
		|		2	2	1
		|	1	2	?
		*/
		coeff[kp]=/*maxerror = 18<<depth*/
			+(errors[kp]>>3)		/*+I	*8*/
			+cN	[kp-1*4*NPREDS]*2	/*+eNW	*2*/
			+cN	[kp+0*4*NPREDS]*2	/*+eN	*2*/
			+ccurr	[kp-1*4*NPREDS]*2	/*+eW	*2*/
			+ccurr	[kp+0*4*NPREDS]		/*+eNN	*1*/
			+ccurr	[kp+1*4*NPREDS]		/*+eNNE	*1*/
			+cN	[kp+1*4*NPREDS]		/*+eNE	*1*/
			+ccurr	[kp-2*4*NPREDS]		/*+eWW	*1*/
			+1
		;
	}
	for(k=0;k<NPREDS;++k)
		sh[k]=31-(DENBITS-1)-_lzcnt_u32(coeff[k]+1);/*invert errros*/
	for(k=0;k<NPREDS;++k)
		if(sh[k]<0)sh[k]=0;
	for(k=0;k<NPREDS;++k)
		coeff[k]=(wg_weights[k]*divlookup[coeff[k]>>sh[k]]>>sh[k])+(1<<DENBITS>>2)/NPREDS;
	wsum=0;
	for(k=0;k<NPREDS;++k)
		wsum+=coeff[k];
	sh[0]=31-(DENBITS-2)-LZCNT32(wsum);/*invert coeff sum*/
	wsum=0;
	for(k=0;k<NPREDS;++k)
		coeff[k]>>=sh[0];
	wsum=0;
	for(k=0;k<NPREDS;++k)
		wsum+=coeff[k];
	ipred=0;
	for(k=0;k<NPREDS;++k)
		ipred+=coeff[k]*preds[k];
	pred=ipred*divlookup[wsum-1]>>NUMBITS;
	ipred=(ipred+16)>>5;
#endif

	vmax=N; vmin=W;
	if(N<W)vmin=N, vmax=W;
	if(vmin>NE)vmin=NE;
	if(vmax<NE)vmax=NE;
	CLAMP2(pred, vmin, vmax);

	(void)coeffs;
	(void)subpreds;
	(void)mixer;
	(void)sse;
	(void)depth;
	(void)kc;
	return (int)pred;
}
static void wg_update(
	int curr, int pred,
	int *subpreds, int *mixer, int *sse,
	const int *preds,
	int *coeffs, uint16_t *cN, uint16_t *ccurr,
	int *errors, int depth
)
{
	int kp, e2[NPREDS], best;

	best=0x7FFFFFFF;
	for(kp=0;kp<NPREDS;++kp)
	{
		int e=abs(curr-preds[kp]);
		if(best>e)
			best=e;
		e2[kp]=e;
	}
	for(kp=0;kp<NPREDS;++kp)
	{
		int e=e2[kp]-best;
		errors[kp]+=((e<<6)-errors[kp]+(1<<3>>1))>>3;
		ccurr[kp]=e;
	}

	(void)pred;
	(void)coeffs;
	(void)cN;
	(void)subpreds;
	(void)mixer;
	(void)sse;
	(void)depth;
}
void pred_wgrad6(Image *src, int fwd)
{
	ALIGN(32) int wg_coeffs[4][NPREDS]={0}, wg_errors[4][NPREDS]={0}, wg_preds[NPREDS*2]={0};
	int amin[4], amax[4], subpreds[2]={0}, mixer[4][NPREDS]={0}, sse[4][NPREDS]={0}, nch, pesize, bufsize, ky, kx, kc, idx;
	uint16_t *ebuf, *erows[2];
	int16_t *pixels, *rows[4];

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;

	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	bufsize=(src->iw+16)*(int)sizeof(int[4*4*2]);/*4 padded rows * 4 channels max*/
	pixels=(int16_t*)malloc(bufsize);
	pesize=(src->iw+16)*(int)sizeof(int16_t[2][4][NPREDS]);/*2 padded rows * 4 channels * NPREDS*/
	ebuf=(uint16_t*)_mm_malloc(pesize, sizeof(__m256i));
	if(!pixels||!ebuf)
	{
		CRASH("Alloc error");
		return;
	}
	memset(ebuf, 0, pesize);
	memset(pixels, 0, bufsize);
	FILLMEM((int*)wg_coeffs, 0x10000/NPREDS, sizeof(wg_coeffs), sizeof(int));
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int pred, curr, val;

		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4*2;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4*2;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4*2;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4*2;
		erows[0]=ebuf+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&1)+8)*4*NPREDS;
		erows[1]=ebuf+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&1)+8)*4*NPREDS;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<src->nch;++kc)
			{
				pred=wg_predict(wg_errors[kc]
					, wg_coeffs[kc]
					, erows[1]+(ptrdiff_t)NPREDS*kc
					, erows[0]+(ptrdiff_t)NPREDS*kc
					, rows[3]+kc
					, rows[2]+kc
					, rows[1]+kc
					, rows[0]+kc
					, subpreds
					, mixer[kc], sse[kc]
					, wg_preds
					, src->depth[kc]
					, kc
				);
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
				((int32_t*)HDR_DATA(src))[idx+kc]=val;
				rows[0][kc+0]=curr;
				rows[0][kc+4]=curr-pred;

				wg_update(rows[0][kc+0]
					, pred
					, subpreds
					, mixer[kc]
					, sse[kc]
					, wg_preds
					, wg_coeffs[kc]
					, erows[1]+(ptrdiff_t)NPREDS*kc
					, erows[0]+(ptrdiff_t)NPREDS*kc
					, wg_errors[kc]
					, src->depth[kc]
				);
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
			erows[0]+=4*NPREDS;
			erows[1]+=4*NPREDS;
		}
	}
	free(pixels);
	_mm_free(ebuf);
}
