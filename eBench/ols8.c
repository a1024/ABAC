#include"ebench.h"
#include<stdint.h>
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


#if 1
#define PREDLIST\
	PRED(100000, N+W-NW)\
	PRED(200000, N+W-NW)\
	PRED(100000, N)\
	PRED(100000, W)\
	PRED(100000, 3*(N-NN)+NNN)\
	PRED(100000, 3*(W-WW)+WWW)\
	PRED(100000, W+NE-N)\
	PRED(100000, N+NE-NNE)\
	PRED(100000, W+((NEEE+NEEEEE-N-W)>>3))\
	PRED( 50000, W+NW-NWW)\
	PRED( 50000, N+NW-NNW)\
	PRED( 50000, NE+NEE-NNEEE)\
	PRED( 50000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)\
	PRED( 40000, NW)\
	PRED( 40000, NE)\
	PRED( 40000, NEEE)\
	PRED( 40000, NN)\
	PRED( 40000, WW)\
	PRED( 40000, NNN)\
	PRED( 40000, WWW)\

#endif

#define PREDLIST2 PREDLIST
#if 0
#define PREDLIST2\
	PRED(100000, N+W-NW)\
	PRED(200000, N+W-NW)\
	PRED(100000, N)\
	PRED(100000, W)\
	PRED(100000, 3*(N-NN)+NNN)\
	PRED(100000, 3*(W-WW)+WWW)\
	PRED(100000, W+NE-N)\
	PRED(100000, N+NE-NNE)\
	PRED(100000, W+((NEEE+NEEEEE-N-W)>>3))\
	PRED( 50000, W+NW-NWW)\
	PRED( 50000, N+NW-NNW)\
	PRED( 50000, NE+NEE-NNEEE)\
	PRED( 50000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)\
	PRED( 40000, NW)\
	PRED( 40000, NE)\
	PRED( 40000, NN)\
	PRED( 40000, WW)\
	PRED( 40000, NNN)\
	PRED( 40000, WWW)\

#endif
#if 0
#define PREDLIST2\
	PRED(100000, eN)\
	PRED(100000, eW)\
	PRED(100000, eNW)\
	PRED(100000, eNE)\
	PRED(100000, eNN)\
	PRED(100000, eWW)\

#endif
#if 0
#define PREDLIST2\
	PRED(100000, eN+eW-eNW)\
	PRED(100000, eN)\
	PRED(100000, eW)\
	PRED(100000, 3*(eN-eNN)+eNNN)\
	PRED(100000, 3*(eW-eWW)+eWWW)\
	PRED(100000, eW+eNE-eN)\
	PRED(100000, eN+eNE-eNNE)\
	PRED(100000, eW+((eNEEE+eNEEEEE-eN-eW)>>3))\
	PRED( 50000, eW+eNW-eNWW)\
	PRED( 50000, eN+eNW-eNNW)\
	PRED( 50000, eNE+eNEE-eNNEEE)\
	PRED( 50000, (eWWWWW+eWW-eW+eNNN+eN+eNEEEEE)>>2)
#endif

enum
{
	CTXBITS1=0,
	CTXBITS2=4,

	L1SH=19,
	L1SH2=17,
#define PRED(W, E) +1
	NPREDS=PREDLIST,
	NPREDS2=PREDLIST2
#undef  PRED
};
void pred_ols8(Image *src, int fwd)
{
	const int wsize=sizeof(int64_t[4][1<<CTXBITS1][NPREDS+1]);
	const int w2size=sizeof(int64_t[4][1<<CTXBITS2][NPREDS2+1]);
	int64_t *weights, *weights2;
	int amin[4], amax[4], invdist, kx, ky, kc, idx;
	int bufsize;
	int16_t *pixels, *rows[6];

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[0]=(1<<src->depth[1]>>1)-1;
	amax[0]=(1<<src->depth[2]>>1)-1;
	amax[0]=(1<<src->depth[3]>>1)-1;
	invdist=((1<<16)+g_dist-1)/g_dist;
	weights=(int64_t*)malloc(wsize);
	weights2=(int64_t*)malloc(w2size);
	bufsize=(src->iw+8*2)*(int)sizeof(int16_t[6*4*3]);/*6 padded rows * 4 channels max * {pixels, residuals1, residuals2}*/
	pixels=(int16_t*)malloc(bufsize);
	if(!pixels||!weights||!weights2)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	memset(weights, 0, wsize);
	memset(weights2, 0, w2size);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)0+6)%6)+8)*4-1)*3;/*FIXME FIXME FIXME FIXME*/
		rows[1]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)1+6)%6)+8)*4-1)*3;
		rows[2]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)2+6)%6)+8)*4-1)*3;
		rows[3]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)3+6)%6)+8)*4-1)*3;
		rows[4]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)4+6)%6)+8)*4-1)*3;
		rows[5]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)5+6)%6)+8)*4-1)*3;
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				int64_t *currw, pred1, *currw2, pred2, pred2s, predc;
				int preds[NPREDS], preds2[NPREDS2], e, j, k, vmax, vmin, curr;

				rows[0]+=3;
				rows[1]+=3;
				rows[2]+=3;
				rows[3]+=3;
				rows[4]+=3;
				rows[5]+=3;
				if(!src->depth[kc])
					continue;
				{
					int
						NNNNN		=rows[5][+0*4*3],
						NNNNWW		=rows[4][-2*4*3],
						NNNNW		=rows[4][-1*4*3],
						NNNN		=rows[4][+0*4*3],
						NNNNE		=rows[4][+1*4*3],
						NNNNEE		=rows[4][+2*4*3],
						NNNNEEEE	=rows[4][+4*4*3],
						NNNWWW		=rows[3][-3*4*3],
						NNNW		=rows[3][-1*4*3],
						NNN		=rows[3][+0*4*3],
						NNNE		=rows[3][+1*4*3],
						NNNEE		=rows[3][+2*4*3],
						NNNEEE		=rows[3][+3*4*3],
						NNNEEEE		=rows[3][+4*4*3],
						NNWWWW		=rows[2][-4*4*3],
						NNWWW		=rows[2][-3*4*3],
						NNWW		=rows[2][-2*4*3],
						NNW		=rows[2][-1*4*3],
						NN		=rows[2][+0*4*3],
						NNE		=rows[2][+1*4*3],
						NNEE		=rows[2][+2*4*3],
						NNEEE		=rows[2][+3*4*3],
						NNEEEE		=rows[2][+4*4*3],
						NWWWW		=rows[1][-4*4*3],
						NWWW		=rows[1][-3*4*3],
						NWW		=rows[1][-2*4*3],
						NW		=rows[1][-1*4*3],
						N		=rows[1][+0*4*3],
						NE		=rows[1][+1*4*3],
						NEE		=rows[1][+2*4*3],
						NEEE		=rows[1][+3*4*3],
						NEEEE		=rows[1][+4*4*3],
						NEEEEE		=rows[1][+5*4*3],
						NEEEEEE		=rows[1][+6*4*3],
						NEEEEEEE	=rows[1][+7*4*3],
						NEEEEEEEE	=rows[1][+8*4*3],
						WWWWWWWWW	=rows[0][-9*4*3],
						WWWWWWWW	=rows[0][-8*4*3],
						WWWWWWW		=rows[0][-7*4*3],
						WWWWWW		=rows[0][-6*4*3],
						WWWWW		=rows[0][-5*4*3],
						WWWW		=rows[0][-4*4*3],
						WWW		=rows[0][-3*4*3],
						WW		=rows[0][-2*4*3],
						W		=rows[0][-1*4*3];

#if 1
					(void)NNNNN;
					(void)NNNNWW;
					(void)NNNNW;
					(void)NNNN;
					(void)NNNNE;
					(void)NNNNEE;
					(void)NNNNEEEE;
					(void)NNNWWW;
					(void)NNNW;
					(void)NNN;
					(void)NNNE;
					(void)NNNEE;
					(void)NNNEEE;
					(void)NNNEEEE;
					(void)NNWWWW;
					(void)NNWWW;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
					(void)NNEEE;
					(void)NNEEEE;
					(void)NWWWW;
					(void)NWWW;
					(void)NWW;
					(void)NW;
					(void)N;
					(void)NE;
					(void)NEE;
					(void)NEEE;
					(void)NEEEE;
					(void)NEEEEE;
					(void)NEEEEEE;
					(void)NEEEEEEE;
					(void)NEEEEEEEE;
					(void)WWWWWWWWW;
					(void)WWWWWWWW;
					(void)WWWWWWW;
					(void)WWWWWW;
					(void)WWWWW;
					(void)WWWW;
					(void)WWW;
					(void)WW;
					(void)W;
#endif
#define PRED(W0, E) preds[j++]=E;
					j=0;
					PREDLIST;
#undef  PRED
					vmax=N; vmin=W;
					if(N<W)vmin=N, vmax=W;
					CLAMP2(preds[0], vmin, vmax);
					currw=weights+(NPREDS+1)*kc;
					pred1=currw[NPREDS];
					for(k=0;k<NPREDS;++k)
						pred1+=currw[k]*preds[k];
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					currw2=weights2+(NPREDS2+1)*((1<<CTXBITS2)*kc+((N+W)/2<<CTXBITS2>>src->depth[kc]&((1<<CTXBITS2)-1)));
				}
				{
					int
						NNNNN		=rows[5][1+0*4*3],
						NNNNWW		=rows[4][1-2*4*3],
						NNNNW		=rows[4][1-1*4*3],
						NNNN		=rows[4][1+0*4*3],
						NNNNE		=rows[4][1+1*4*3],
						NNNNEE		=rows[4][1+2*4*3],
						NNNNEEEE	=rows[4][1+4*4*3],
						NNNWWW		=rows[3][1-3*4*3],
						NNNW		=rows[3][1-1*4*3],
						NNN		=rows[3][1+0*4*3],
						NNNE		=rows[3][1+1*4*3],
						NNNEE		=rows[3][1+2*4*3],
						NNNEEE		=rows[3][1+3*4*3],
						NNNEEEE		=rows[3][1+4*4*3],
						NNWWWW		=rows[2][1-4*4*3],
						NNWWW		=rows[2][1-3*4*3],
						NNWW		=rows[2][1-2*4*3],
						NNW		=rows[2][1-1*4*3],
						NN		=rows[2][1+0*4*3],
						NNE		=rows[2][1+1*4*3],
						NNEE		=rows[2][1+2*4*3],
						NNEEE		=rows[2][1+3*4*3],
						NNEEEE		=rows[2][1+4*4*3],
						NWWWW		=rows[1][1-4*4*3],
						NWWW		=rows[1][1-3*4*3],
						NWW		=rows[1][1-2*4*3],
						NW		=rows[1][1-1*4*3],
						N		=rows[1][1+0*4*3],
						NE		=rows[1][1+1*4*3],
						NEE		=rows[1][1+2*4*3],
						NEEE		=rows[1][1+3*4*3],
						NEEEE		=rows[1][1+4*4*3],
						NEEEEE		=rows[1][1+5*4*3],
						NEEEEEE		=rows[1][1+6*4*3],
						NEEEEEEE	=rows[1][1+7*4*3],
						NEEEEEEEE	=rows[1][1+8*4*3],
						WWWWWWWWW	=rows[0][1-9*4*3],
						WWWWWWWW	=rows[0][1-8*4*3],
						WWWWWWW		=rows[0][1-7*4*3],
						WWWWWW		=rows[0][1-6*4*3],
						WWWWW		=rows[0][1-5*4*3],
						WWWW		=rows[0][1-4*4*3],
						WWW		=rows[0][1-3*4*3],
						WW		=rows[0][1-2*4*3],
						W		=rows[0][1-1*4*3];

#if 1
					(void)NNNNN;
					(void)NNNNWW;
					(void)NNNNW;
					(void)NNNN;
					(void)NNNNE;
					(void)NNNNEE;
					(void)NNNNEEEE;
					(void)NNNWWW;
					(void)NNNW;
					(void)NNN;
					(void)NNNE;
					(void)NNNEE;
					(void)NNNEEE;
					(void)NNNEEEE;
					(void)NNWWWW;
					(void)NNWWW;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
					(void)NNEEE;
					(void)NNEEEE;
					(void)NWWWW;
					(void)NWWW;
					(void)NWW;
					(void)NW;
					(void)N;
					(void)NE;
					(void)NEE;
					(void)NEEE;
					(void)NEEEE;
					(void)NEEEEE;
					(void)NEEEEEE;
					(void)NEEEEEEE;
					(void)NEEEEEEEE;
					(void)WWWWWWWWW;
					(void)WWWWWWWW;
					(void)WWWWWWW;
					(void)WWWWWW;
					(void)WWWWW;
					(void)WWWW;
					(void)WWW;
					(void)WW;
					(void)W;
#endif
#define PRED(W0, E) preds2[j++]=E;
					j=0;
					PREDLIST2;
#undef  PRED
				}
				pred2=currw2[NPREDS2];
				for(k=0;k<NPREDS2;++k)
					pred2+=currw2[k]*preds2[k];
				pred1+=pred2<<3;/*[sic]*/
				pred1+=1<<L1SH>>1;
				pred1>>=L1SH;
				pred2+=1<<L1SH2>>1;
				pred2>>=L1SH2;

				pred2s=pred1+pred2;
				predc=pred2s;
				CLAMP2(predc, vmin, vmax);
				
				curr=((int32_t*)HDR_DATA(src))[idx];
				if(g_dist>1)
				{
					if(fwd)
					{
						curr-=(int)predc;
						curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));/*curr/=g_dist*/
						((int32_t*)HDR_DATA(src))[idx]=curr;

						curr=g_dist*curr+(int)predc;
						CLAMP2(curr, amin[kc], amax[kc]);
					}
					else
					{
						curr=g_dist*curr+(int)predc;
						CLAMP2(curr, amin[kc], amax[kc]);

						((int32_t*)HDR_DATA(src))[idx]=curr;
					}
				}
				else
				{
					if(fwd)
					{
						e=curr-(int)predc;
						e<<=32-src->depth[kc];
						e>>=32-src->depth[kc];
						((int32_t*)HDR_DATA(src))[idx]=e;
					}
					else
					{
						curr+=(int)predc;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						((int32_t*)HDR_DATA(src))[idx]=curr;
					}
				}
				rows[0][0]=curr;
				rows[0][1]=curr-(int)pred1;
				rows[0][2]=curr-(int)pred2s;

				/*update*/
				e=(curr>pred1)-(curr<pred1);
				currw[NPREDS]+=e;/*bias*/
				for(k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];/*coeffs*/

				e=(curr>pred2s)-(curr<pred2s);
				currw2[NPREDS2]+=e;
				for(k=0;k<NPREDS2;++k)
					currw2[k]+=e*preds2[k];
			}
		}
	}
	free(pixels);
	free(weights);
	free(weights2);
}
void pred_ols8_crct(Image *src, int fwd)
{
	enum
	{
		XPAD=16,
		NROWS=8,
		NCH=4,
		NVAL=3
	};
	const int wsize=(int)sizeof(int32_t[4][1<<CTXBITS1][NPREDS+1]);
	const int w2size=(int)sizeof(int32_t[4][1<<CTXBITS2][NPREDS2+1]);
	int32_t *weights, *weights2;
	int amin[4], amax[4], invdist, ky, kx, kc, idx, yuv[3];
	int psize;
	int32_t *pixels, *rows[8];
	int uc0, vc0, vc1;
	const uint8_t *perm=0;
	int a_dist[4]={0}, dist2symx[4]={0};
#if 0
	uint64_t getrand=0xC2B2AE3D27D4EB4F;
#endif

	(void)invdist;
	(void)dist2symx;
	a_dist[0]=g_dist;
	a_dist[1]=g_dist;
	a_dist[2]=g_dist;
	a_dist[3]=g_dist;
	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	invdist=((1<<16)+g_dist-1)/g_dist;
	weights=(int32_t*)malloc(wsize);
	weights2=(int32_t*)malloc(w2size);
	psize=(src->iw+2*XPAD)*(int)sizeof(int32_t[NROWS*NCH*NVAL]);
	pixels=(int32_t*)malloc(psize);
	uc0=vc0=vc1=0;
	perm=0;
	if(fwd)
		src->rct=crct2_analysis(src);
	perm=crct2_unpack(src->rct, &uc0, &vc0, &vc1);
	if(!pixels||!weights||!weights2)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	memset(weights, 0, wsize);
	memset(weights2, 0, w2size);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;/*sub 1 channel for pre-increment*/
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		rows[4]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)4+NROWS)%NROWS)*NVAL;
		rows[5]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)5+NROWS)%NROWS)*NVAL;
		rows[6]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)6+NROWS)%NROWS)*NVAL;
		rows[7]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)7+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			int offset=0;

			yuv[0]=((int32_t*)HDR_DATA(src))[idx+perm[0]];
			yuv[1]=((int32_t*)HDR_DATA(src))[idx+perm[1]];
			yuv[2]=((int32_t*)HDR_DATA(src))[idx+perm[2]];
			for(kc=0;kc<4;++kc)
			{
				int preds[NPREDS], preds2[NPREDS2], j, k, vmax, vmin, *currw, pred1, *currw2, pred2, pred2s, predc, curr, e;
				
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				rows[4]+=NROWS*NVAL;
				rows[5]+=NROWS*NVAL;
				rows[6]+=NROWS*NVAL;
				rows[7]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				{
					int
						NNNNN		=rows[5][0+0*NCH*NROWS*NVAL],
						NNNNWW		=rows[4][0-2*NCH*NROWS*NVAL],
						NNNNW		=rows[4][0-1*NCH*NROWS*NVAL],
						NNNN		=rows[4][0+0*NCH*NROWS*NVAL],
						NNNNE		=rows[4][0+1*NCH*NROWS*NVAL],
						NNNNEE		=rows[4][0+2*NCH*NROWS*NVAL],
						NNNNEEEE	=rows[4][0+4*NCH*NROWS*NVAL],
						NNNWWW		=rows[3][0-3*NCH*NROWS*NVAL],
						NNNW		=rows[3][0-1*NCH*NROWS*NVAL],
						NNN		=rows[3][0+0*NCH*NROWS*NVAL],
						NNNE		=rows[3][0+1*NCH*NROWS*NVAL],
						NNNEE		=rows[3][0+2*NCH*NROWS*NVAL],
						NNNEEE		=rows[3][0+3*NCH*NROWS*NVAL],
						NNNEEEE		=rows[3][0+4*NCH*NROWS*NVAL],
						NNWWWW		=rows[2][0-4*NCH*NROWS*NVAL],
						NNWWW		=rows[2][0-3*NCH*NROWS*NVAL],
						NNWW		=rows[2][0-2*NCH*NROWS*NVAL],
						NNW		=rows[2][0-1*NCH*NROWS*NVAL],
						NN		=rows[2][0+0*NCH*NROWS*NVAL],
						NNE		=rows[2][0+1*NCH*NROWS*NVAL],
						NNEE		=rows[2][0+2*NCH*NROWS*NVAL],
						NNEEE		=rows[2][0+3*NCH*NROWS*NVAL],
						NNEEEE		=rows[2][0+4*NCH*NROWS*NVAL],
						NWWWW		=rows[1][0-4*NCH*NROWS*NVAL],
						NWWW		=rows[1][0-3*NCH*NROWS*NVAL],
						NWW		=rows[1][0-2*NCH*NROWS*NVAL],
						NW		=rows[1][0-1*NCH*NROWS*NVAL],
						N		=rows[1][0+0*NCH*NROWS*NVAL],
						NE		=rows[1][0+1*NCH*NROWS*NVAL],
						NEE		=rows[1][0+2*NCH*NROWS*NVAL],
						NEEE		=rows[1][0+3*NCH*NROWS*NVAL],
						NEEEE		=rows[1][0+4*NCH*NROWS*NVAL],
						NEEEEE		=rows[1][0+5*NCH*NROWS*NVAL],
						NEEEEEE		=rows[1][0+6*NCH*NROWS*NVAL],
						NEEEEEEE	=rows[1][0+7*NCH*NROWS*NVAL],
						NEEEEEEEE	=rows[1][0+8*NCH*NROWS*NVAL],
						WWWWWWWWW	=rows[0][0-9*NCH*NROWS*NVAL],
						WWWWWWWW	=rows[0][0-8*NCH*NROWS*NVAL],
						WWWWWWW		=rows[0][0-7*NCH*NROWS*NVAL],
						WWWWWW		=rows[0][0-6*NCH*NROWS*NVAL],
						WWWWW		=rows[0][0-5*NCH*NROWS*NVAL],
						WWWW		=rows[0][0-4*NCH*NROWS*NVAL],
						WWW		=rows[0][0-3*NCH*NROWS*NVAL],
						WW		=rows[0][0-2*NCH*NROWS*NVAL],
						W		=rows[0][0-1*NCH*NROWS*NVAL];

#if 1
					(void)NNNNN;
					(void)NNNNWW;
					(void)NNNNW;
					(void)NNNN;
					(void)NNNNE;
					(void)NNNNEE;
					(void)NNNNEEEE;
					(void)NNNWWW;
					(void)NNNW;
					(void)NNN;
					(void)NNNE;
					(void)NNNEE;
					(void)NNNEEE;
					(void)NNNEEEE;
					(void)NNWWWW;
					(void)NNWWW;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
					(void)NNEEE;
					(void)NNEEEE;
					(void)NWWWW;
					(void)NWWW;
					(void)NWW;
					(void)NW;
					(void)N;
					(void)NE;
					(void)NEE;
					(void)NEEE;
					(void)NEEEE;
					(void)NEEEEE;
					(void)NEEEEEE;
					(void)NEEEEEEE;
					(void)NEEEEEEEE;
					(void)WWWWWWWWW;
					(void)WWWWWWWW;
					(void)WWWWWWW;
					(void)WWWWWW;
					(void)WWWWW;
					(void)WWWW;
					(void)WWW;
					(void)WW;
					(void)W;
#endif
#define PRED(W0, E) preds[j++]=E;
					j=0;
					PREDLIST;
#undef  PRED
					vmax=N; vmin=W;
					if(N<W)vmin=N, vmax=W;
					CLAMP2(preds[0], vmin, vmax);
					currw=weights+(NPREDS+1)*kc;
					pred1=currw[NPREDS];
					for(k=0;k<NPREDS;++k)
						pred1+=currw[k]*preds[k];
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					currw2=weights2+(NPREDS2+1)*((1<<CTXBITS2)*kc+((N+W)/2<<CTXBITS2>>src->depth[kc]&((1<<CTXBITS2)-1)));
				}
				{
					int
						NNNNN		=rows[5][1+0*NCH*NROWS*NVAL],
						NNNNWW		=rows[4][1-2*NCH*NROWS*NVAL],
						NNNNW		=rows[4][1-1*NCH*NROWS*NVAL],
						NNNN		=rows[4][1+0*NCH*NROWS*NVAL],
						NNNNE		=rows[4][1+1*NCH*NROWS*NVAL],
						NNNNEE		=rows[4][1+2*NCH*NROWS*NVAL],
						NNNNEEEE	=rows[4][1+4*NCH*NROWS*NVAL],
						NNNWWW		=rows[3][1-3*NCH*NROWS*NVAL],
						NNNW		=rows[3][1-1*NCH*NROWS*NVAL],
						NNN		=rows[3][1+0*NCH*NROWS*NVAL],
						NNNE		=rows[3][1+1*NCH*NROWS*NVAL],
						NNNEE		=rows[3][1+2*NCH*NROWS*NVAL],
						NNNEEE		=rows[3][1+3*NCH*NROWS*NVAL],
						NNNEEEE		=rows[3][1+4*NCH*NROWS*NVAL],
						NNWWWW		=rows[2][1-4*NCH*NROWS*NVAL],
						NNWWW		=rows[2][1-3*NCH*NROWS*NVAL],
						NNWW		=rows[2][1-2*NCH*NROWS*NVAL],
						NNW		=rows[2][1-1*NCH*NROWS*NVAL],
						NN		=rows[2][1+0*NCH*NROWS*NVAL],
						NNE		=rows[2][1+1*NCH*NROWS*NVAL],
						NNEE		=rows[2][1+2*NCH*NROWS*NVAL],
						NNEEE		=rows[2][1+3*NCH*NROWS*NVAL],
						NNEEEE		=rows[2][1+4*NCH*NROWS*NVAL],
						NWWWW		=rows[1][1-4*NCH*NROWS*NVAL],
						NWWW		=rows[1][1-3*NCH*NROWS*NVAL],
						NWW		=rows[1][1-2*NCH*NROWS*NVAL],
						NW		=rows[1][1-1*NCH*NROWS*NVAL],
						N		=rows[1][1+0*NCH*NROWS*NVAL],
						NE		=rows[1][1+1*NCH*NROWS*NVAL],
						NEE		=rows[1][1+2*NCH*NROWS*NVAL],
						NEEE		=rows[1][1+3*NCH*NROWS*NVAL],
						NEEEE		=rows[1][1+4*NCH*NROWS*NVAL],
						NEEEEE		=rows[1][1+5*NCH*NROWS*NVAL],
						NEEEEEE		=rows[1][1+6*NCH*NROWS*NVAL],
						NEEEEEEE	=rows[1][1+7*NCH*NROWS*NVAL],
						NEEEEEEEE	=rows[1][1+8*NCH*NROWS*NVAL],
						WWWWWWWWW	=rows[0][1-9*NCH*NROWS*NVAL],
						WWWWWWWW	=rows[0][1-8*NCH*NROWS*NVAL],
						WWWWWWW		=rows[0][1-7*NCH*NROWS*NVAL],
						WWWWWW		=rows[0][1-6*NCH*NROWS*NVAL],
						WWWWW		=rows[0][1-5*NCH*NROWS*NVAL],
						WWWW		=rows[0][1-4*NCH*NROWS*NVAL],
						WWW		=rows[0][1-3*NCH*NROWS*NVAL],
						WW		=rows[0][1-2*NCH*NROWS*NVAL],
						W		=rows[0][1-1*NCH*NROWS*NVAL];

#if 1
					(void)NNNNN;
					(void)NNNNWW;
					(void)NNNNW;
					(void)NNNN;
					(void)NNNNE;
					(void)NNNNEE;
					(void)NNNNEEEE;
					(void)NNNWWW;
					(void)NNNW;
					(void)NNN;
					(void)NNNE;
					(void)NNNEE;
					(void)NNNEEE;
					(void)NNNEEEE;
					(void)NNWWWW;
					(void)NNWWW;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
					(void)NNEEE;
					(void)NNEEEE;
					(void)NWWWW;
					(void)NWWW;
					(void)NWW;
					(void)NW;
					(void)N;
					(void)NE;
					(void)NEE;
					(void)NEEE;
					(void)NEEEE;
					(void)NEEEEE;
					(void)NEEEEEE;
					(void)NEEEEEEE;
					(void)NEEEEEEEE;
					(void)WWWWWWWWW;
					(void)WWWWWWWW;
					(void)WWWWWWW;
					(void)WWWWWW;
					(void)WWWWW;
					(void)WWWW;
					(void)WWW;
					(void)WW;
					(void)W;
#endif
#define PRED(W0, E) preds2[j++]=E;
					j=0;
					PREDLIST2;
#undef  PRED
				}
				pred2=currw2[NPREDS2];
				for(k=0;k<NPREDS2;++k)
					pred2+=currw2[k]*preds2[k];

				pred1+=pred2<<3;/*[sic]*/
				pred1+=1<<L1SH>>1;
				pred1>>=L1SH;
				pred2+=1<<L1SH2>>1;
				pred2>>=L1SH2;

				pred2s=pred1+pred2;
				predc=pred2s;
				CLAMP2(predc, vmin, vmax);

				offset=0;
				if(kc==1)offset=uc0*yuv[0]>>RCTBITS;
				if(kc==2)offset=(vc0*yuv[0]+vc1*yuv[1])>>RCTBITS;
				predc+=offset;
				CLAMP2(predc, amin[kc], amax[kc]);
				
				curr=yuv[kc];
				if(g_dist>1)
				{
					if(fwd)
					{
						curr-=predc;
						curr/=a_dist[kc];
					/*	curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));*/ /*curr/=g_dist*/
						((int32_t*)HDR_DATA(src))[idx+kc]=curr;
					}
					else
						curr=((int32_t*)HDR_DATA(src))[idx+kc];
					{
						int r0=curr;
#if 0
						int inc=0;

						rows[0][2]=curr;
						if(curr)
							inc=1;
						else if(
							!rows[0][2+0*NCH*NROWS*NVAL]&&
							!rows[0][2+1*NCH*NROWS*NVAL]&&
							!rows[0][2+2*NCH*NROWS*NVAL]&&
							!rows[0][2+3*NCH*NROWS*NVAL]&&
							!rows[0][2+4*NCH*NROWS*NVAL]&&
							!rows[0][2+5*NCH*NROWS*NVAL]&&
							!rows[0][2+6*NCH*NROWS*NVAL]&&

							!rows[5][2-7*NCH*NROWS*NVAL]&&
							!rows[5][2-6*NCH*NROWS*NVAL]&&
							!rows[5][2-5*NCH*NROWS*NVAL]&&
							!rows[5][2-4*NCH*NROWS*NVAL]&&
							!rows[5][2-3*NCH*NROWS*NVAL]&&
							!rows[5][2-2*NCH*NROWS*NVAL]&&
							!rows[5][2-1*NCH*NROWS*NVAL]&&
							!rows[5][2+0*NCH*NROWS*NVAL]&&
							!rows[5][2+1*NCH*NROWS*NVAL]&&
							!rows[5][2+2*NCH*NROWS*NVAL]&&
							!rows[5][2+3*NCH*NROWS*NVAL]&&
							!rows[5][2+4*NCH*NROWS*NVAL]&&
							!rows[5][2+5*NCH*NROWS*NVAL]&&
							!rows[5][2+6*NCH*NROWS*NVAL]&&
							!rows[5][2+7*NCH*NROWS*NVAL]&&

							!rows[4][2-7*NCH*NROWS*NVAL]&&
							!rows[4][2-6*NCH*NROWS*NVAL]&&
							!rows[4][2-5*NCH*NROWS*NVAL]&&
							!rows[4][2-4*NCH*NROWS*NVAL]&&
							!rows[4][2-3*NCH*NROWS*NVAL]&&
							!rows[4][2-2*NCH*NROWS*NVAL]&&
							!rows[4][2-1*NCH*NROWS*NVAL]&&
							!rows[4][2+0*NCH*NROWS*NVAL]&&
							!rows[4][2+1*NCH*NROWS*NVAL]&&
							!rows[4][2+2*NCH*NROWS*NVAL]&&
							!rows[4][2+3*NCH*NROWS*NVAL]&&
							!rows[4][2+4*NCH*NROWS*NVAL]&&
							!rows[4][2+5*NCH*NROWS*NVAL]&&
							!rows[4][2+6*NCH*NROWS*NVAL]&&
							!rows[4][2+7*NCH*NROWS*NVAL]&&

							!rows[3][2-8*NCH*NROWS*NVAL]&&
							!rows[3][2-7*NCH*NROWS*NVAL]&&
							!rows[3][2-6*NCH*NROWS*NVAL]&&
							!rows[3][2-5*NCH*NROWS*NVAL]&&
							!rows[3][2-4*NCH*NROWS*NVAL]&&
							!rows[3][2-3*NCH*NROWS*NVAL]&&
							!rows[3][2-2*NCH*NROWS*NVAL]&&
							!rows[3][2-1*NCH*NROWS*NVAL]&&
							!rows[3][2+0*NCH*NROWS*NVAL]&&
							!rows[3][2+1*NCH*NROWS*NVAL]&&
							!rows[3][2+2*NCH*NROWS*NVAL]&&
							!rows[3][2+3*NCH*NROWS*NVAL]&&
							!rows[3][2+4*NCH*NROWS*NVAL]&&
							!rows[3][2+5*NCH*NROWS*NVAL]&&
							!rows[3][2+6*NCH*NROWS*NVAL]&&
							!rows[3][2+7*NCH*NROWS*NVAL]&&
							!rows[3][2+8*NCH*NROWS*NVAL]&&

							!rows[2][2-8*NCH*NROWS*NVAL]&&
							!rows[2][2-7*NCH*NROWS*NVAL]&&
							!rows[2][2-6*NCH*NROWS*NVAL]&&
							!rows[2][2-5*NCH*NROWS*NVAL]&&
							!rows[2][2-4*NCH*NROWS*NVAL]&&
							!rows[2][2-3*NCH*NROWS*NVAL]&&
							!rows[2][2-2*NCH*NROWS*NVAL]&&
							!rows[2][2-1*NCH*NROWS*NVAL]&&
							!rows[2][2+0*NCH*NROWS*NVAL]&&
							!rows[2][2+1*NCH*NROWS*NVAL]&&
							!rows[2][2+2*NCH*NROWS*NVAL]&&
							!rows[2][2+3*NCH*NROWS*NVAL]&&
							!rows[2][2+4*NCH*NROWS*NVAL]&&
							!rows[2][2+5*NCH*NROWS*NVAL]&&
							!rows[2][2+6*NCH*NROWS*NVAL]&&
							!rows[2][2+7*NCH*NROWS*NVAL]&&
							!rows[2][2+8*NCH*NROWS*NVAL]&&

							!rows[1][2-8*NCH*NROWS*NVAL]&&
							!rows[1][2-7*NCH*NROWS*NVAL]&&
							!rows[1][2-6*NCH*NROWS*NVAL]&&
							!rows[1][2-5*NCH*NROWS*NVAL]&&
							!rows[1][2-4*NCH*NROWS*NVAL]&&
							!rows[1][2-3*NCH*NROWS*NVAL]&&
							!rows[1][2-2*NCH*NROWS*NVAL]&&
							!rows[1][2-1*NCH*NROWS*NVAL]&&
							!rows[1][2+0*NCH*NROWS*NVAL]&&
							!rows[1][2+1*NCH*NROWS*NVAL]&&
							!rows[1][2+2*NCH*NROWS*NVAL]&&
							!rows[1][2+3*NCH*NROWS*NVAL]&&
							!rows[1][2+4*NCH*NROWS*NVAL]&&
							!rows[1][2+5*NCH*NROWS*NVAL]&&
							!rows[1][2+6*NCH*NROWS*NVAL]&&
							!rows[1][2+7*NCH*NROWS*NVAL]&&
							!rows[1][2+8*NCH*NROWS*NVAL]&&

							!rows[0][2-9*NCH*NROWS*NVAL]&&
							!rows[0][2-8*NCH*NROWS*NVAL]&&
							!rows[0][2-7*NCH*NROWS*NVAL]&&
							!rows[0][2-6*NCH*NROWS*NVAL]&&
							!rows[0][2-5*NCH*NROWS*NVAL]&&
							!rows[0][2-4*NCH*NROWS*NVAL]&&
							!rows[0][2-3*NCH*NROWS*NVAL]&&
							!rows[0][2-2*NCH*NROWS*NVAL]&&
							!rows[0][2-1*NCH*NROWS*NVAL]
						)
							inc=-1;
#endif
#if 0
						enum
						{
							MAXRUN=8
						};
						int reset=0;
						int dist2symy=rows[1][2+0*NCH*NROWS*NVAL];
						if(!curr)
						{
							++dist2symx[kc];
							++dist2symy;
						}
						else
							reset=1;
#endif
						curr*=a_dist[kc];

#if 0
						if(curr>0)
							curr+=a_dist[kc]>>1;
						if(curr<0)
							curr-=a_dist[kc]>>1;
#endif
#if 0
						if(curr)
						{
							int inc=getrand%a_dist[kc];
							getrand=getrand<<31|getrand>>(64-31);
							getrand*=0x9E3779B185EBCA87;

							curr+=inc;
						}
#endif

						curr+=predc;
						CLAMP2(curr, amin[kc], amax[kc]);
						yuv[kc]=curr;

#if 1
					/*	rows[0][2+0*NCH*NROWS*NVAL]=rows[1][2+0*NCH*NROWS*NVAL]<<!(ky&1)|(r0!=0);*/ /*x2 register length?*/
						rows[0][2+0*NCH*NROWS*NVAL]=rows[1][2+0*NCH*NROWS*NVAL]<<1|(r0!=0);
						if(!((((uint64_t)1<<16)-1)&(
							rows[1][2+ 1*NCH*NROWS*NVAL]|	//transposed yM reach N
							rows[1][2+ 2*NCH*NROWS*NVAL]|
							rows[1][2+ 3*NCH*NROWS*NVAL]|
							rows[1][2+ 4*NCH*NROWS*NVAL]|
							rows[1][2+ 5*NCH*NROWS*NVAL]|
							rows[1][2+ 6*NCH*NROWS*NVAL]|
							rows[1][2+ 7*NCH*NROWS*NVAL]|
							rows[1][2+ 8*NCH*NROWS*NVAL]|
							rows[0][2- 8*NCH*NROWS*NVAL]|
							rows[0][2- 7*NCH*NROWS*NVAL]|
							rows[0][2- 6*NCH*NROWS*NVAL]|
							rows[0][2- 5*NCH*NROWS*NVAL]|
							rows[0][2- 4*NCH*NROWS*NVAL]|
							rows[0][2- 3*NCH*NROWS*NVAL]|
							rows[0][2- 2*NCH*NROWS*NVAL]|
							rows[0][2- 1*NCH*NROWS*NVAL]|
							rows[0][2+ 0*NCH*NROWS*NVAL]
						)))
						{
							a_dist[kc]-=1;
							if(a_dist[kc]<1)
								a_dist[kc]=1;
						/*	a_dist[kc]-=a_dist[kc]>1;*/
						}
						else if(r0)
						{
							a_dist[kc]+=3;
							if(a_dist[kc]>g_dist)
								a_dist[kc]=g_dist;
						}
						//	a_dist[kc]=g_dist;
#endif
#if 0
					/*	rows[0][2+0*NCH*NROWS*NVAL]=rows[1][2+0*NCH*NROWS*NVAL]<<!(ky&1)|(r0!=0);*/ /*x2 register length?*/
						rows[0][2+0*NCH*NROWS*NVAL]=rows[1][2+0*NCH*NROWS*NVAL]<<1|(r0!=0);
						if(!((((uint64_t)1<<24)-1)&(
							rows[1][2+ 1*NCH*NROWS*NVAL]|	//transposed yM reach N		GDCC/photo-69 d10
							rows[1][2+ 2*NCH*NROWS*NVAL]|
							rows[1][2+ 3*NCH*NROWS*NVAL]|
							rows[1][2+ 4*NCH*NROWS*NVAL]|
							rows[1][2+ 5*NCH*NROWS*NVAL]|
							rows[1][2+ 6*NCH*NROWS*NVAL]|
							rows[1][2+ 7*NCH*NROWS*NVAL]|	//				351736.30 37.810059 ty24r7NE	369832.26 37.990536 y16tr7NE
							rows[1][2+ 8*NCH*NROWS*NVAL]|
						//	rows[1][2+ 9*NCH*NROWS*NVAL]|
						//	rows[1][2+10*NCH*NROWS*NVAL]|
						//	rows[0][2-10*NCH*NROWS*NVAL]|	//319118.90 37.232280 ty32r10
						//	rows[0][2- 9*NCH*NROWS*NVAL]|	//322163.83 37.384869 ty32r09
							rows[0][2- 8*NCH*NROWS*NVAL]|	//325551.29 37.440283 ty32r08
							rows[0][2- 7*NCH*NROWS*NVAL]|	//332067.35 37.542593 ty32r07	342936.67 37.672720 ty24r7	362164.74 37.897771 y16tr7
							rows[0][2- 6*NCH*NROWS*NVAL]|	//				351132.21 37.784647 ty24r6	374346.84 38.012554 y16tr6
							rows[0][2- 5*NCH*NROWS*NVAL]|	//				361990.51 37.873900 ty24r5	388855.65 38.131213 y16tr5
							rows[0][2- 4*NCH*NROWS*NVAL]|	//				377841.88 38.053066 ty24r4	409910.08 38.301898 y16tr4
							rows[0][2- 3*NCH*NROWS*NVAL]|	//								439559.23 38.496528 y16tr3	472932.49 38.700900 y12tr3	530157.13 38.981538 y08tr3
							rows[0][2- 2*NCH*NROWS*NVAL]|	//								485172.92 38.789841 y16tr2	527451.33 39.005420 y12tr2	594467.73 39.313925 y08tr2
							rows[0][2- 1*NCH*NROWS*NVAL]|	//												625410.77 39.489361 y12tr1	708195.64 39.810866 y08tr1		558840.90 39.190841 tr1y12 X
							rows[0][2+ 0*NCH*NROWS*NVAL]	//												893590.28 40.594033 y12tr0	984029.78 40.914294 y08tr0
						)))
							a_dist[kc]-=a_dist[kc]>1;
						else if(r0)
							a_dist[kc]=g_dist;
#endif
#if 0
						rows[0][2+0*NCH*NROWS*NVAL]=rows[0][2-1*NCH*NROWS*NVAL]<<1|(r0!=0);
						if(!(
							(rows[0][2+5*NCH*NROWS*NVAL]&0x3FF)|	//  456922.53 38.594845		GDCC/photo-69 d10
							(rows[7][2+5*NCH*NROWS*NVAL]&0x3FF)|	//  470579.56 38.677359
							(rows[6][2+5*NCH*NROWS*NVAL]&0x3FF)|	//  490922.82 38.766315
							(rows[5][2+5*NCH*NROWS*NVAL]&0x3FF)|	//  513889.88 38.886249
							(rows[4][2+5*NCH*NROWS*NVAL]&0x3FF)|	//  543352.58 39.019337
							(rows[3][2+5*NCH*NROWS*NVAL]&0x3FF)|	//  586799.06 39.185988
							(rows[2][2+5*NCH*NROWS*NVAL]&0x3FF)|	//  650573.70 39.416613
							(rows[1][2+5*NCH*NROWS*NVAL]&0x3FF)|	//  764659.76 39.781853
							(rows[0][2+0*NCH*NROWS*NVAL]&0x1F)	// 1044751.03 40.614498
						))
							a_dist[kc]-=a_dist[kc]>1;
						else if(r0)
							a_dist[kc]=g_dist;
#endif
#if 0
						rows[0][2+0*NCH*NROWS*NVAL]=rows[0][2-1*NCH*NROWS*NVAL]<<1|(r0!=0);
						if(!(
						//	(rows[0][2+7*NCH*NROWS*NVAL]&0x7FFF)|	//405665.82 38.221620 [best]	GDCC/photo-69 d10
						//	(rows[7][2+7*NCH*NROWS*NVAL]&0x7FFF)|	//411168.52 38.252767
						//	(rows[6][2+7*NCH*NROWS*NVAL]&0x7FFF)|	//429789.93 38.389618	circle 449074.84 38.514027
						//	(rows[5][2+7*NCH*NROWS*NVAL]&0x7FFF)|	//445989.62 38.452232	circle 459858.86 38.558775
						//	(rows[4][2+7*NCH*NROWS*NVAL]&0x7FFF)|	//313903.93 35.958377
							(rows[3][2+7*NCH*NROWS*NVAL]&0x7FFF)|	//503858.39 38.720748
							(rows[2][2+7*NCH*NROWS*NVAL]&0x7FFF)|	//558880.59 38.928703
							(rows[1][2+7*NCH*NROWS*NVAL]&0x7FFF)|	//658957.86 39.255624
							(rows[0][2+0*NCH*NROWS*NVAL]&0xFF)	//920677.08 40.045000
						))						//298458.15 36.196825	disabled
							a_dist[kc]-=a_dist[kc]>1;
						else if(r0)
							a_dist[kc]=g_dist;
#endif
#if 0
						rows[0][2+0*NCH*NROWS*NVAL]=rows[0][2-1*NCH*NROWS*NVAL]<<1|(r0!=0);
						if(!(/*XPAD=16*/
							(rows[3][2+9*NCH*NROWS*NVAL]&0x7FFFF)|	//464992.41 38.453202	GDCC/photo-69 d10
							(rows[2][2+9*NCH*NROWS*NVAL]&0x7FFFF)|	//512390.47 38.650679
							(rows[1][2+9*NCH*NROWS*NVAL]&0x7FFFF)|	//602956.41 38.985118
							(rows[0][2+0*NCH*NROWS*NVAL]&0x3FF)	//858994.71 39.787507
						))
							a_dist[kc]-=a_dist[kc]>1;
						else if(r0)
							a_dist[kc]=g_dist;
#endif
#if 0
						rows[0][2+0*NCH*NROWS*NVAL]=rows[0][2-1*NCH*NROWS*NVAL]<<1|(r0!=0);
						if(!(
							(rows[0][2+11*NCH*NROWS*NVAL]&0x7FFFF)|	//380184.38 37.962585		GDCC/photo-69 d10
							(rows[7][2+11*NCH*NROWS*NVAL]&0x7FFFF)|	//389472.16 38.032595
							(rows[6][2+11*NCH*NROWS*NVAL]&0x7FFFF)|	//401211.33 38.106955
							(rows[5][2+11*NCH*NROWS*NVAL]&0x7FFFF)|	//417009.16 38.218038
							(rows[4][2+11*NCH*NROWS*NVAL]&0x7FFFF)|	//438699.09 38.333200
							(rows[3][2+11*NCH*NROWS*NVAL]&0x7FFFF)|	//470224.43 38.501410
							(rows[2][2+11*NCH*NROWS*NVAL]&0x7FFFF)|	//521492.35 38.722718
							(rows[1][2+11*NCH*NROWS*NVAL]&0x7FFFF)|	//623100.37 39.103427
							(rows[0][2+0*NCH*NROWS*NVAL]&0xFF)	//920677.08 40.045000
						))						//
							a_dist[kc]-=a_dist[kc]>1;
						else if(r0)
							a_dist[kc]=g_dist;
#endif
#if 0
						++dist2symx[kc];
						if(dist2symx[kc]>=8)
						{
							a_dist[kc]-=a_dist[kc]>1;
							dist2symx[kc]-=8;
						}
						else if(r0)
						{
							a_dist[kc]=g_dist;
							dist2symx[kc]=0;
						}
#endif
#if 0
						{
							int dN=rows[1][2+0*NCH*NROWS*NVAL]+1;
							int dW=rows[0][2-1*NCH*NROWS*NVAL]+1;
							int dcurr=dN;

						//	++dist2symx[kc];
						//	if(r0)dist2symx[kc]=0;
							if(dcurr>dW)dcurr=dW;
							if(r0)dcurr=0;
							rows[0][2+0*NCH*NROWS*NVAL]=dcurr;
							if(dcurr>=g_dist<<1)
								a_dist[kc]-=a_dist[kc]>1;
							else if(r0)
								a_dist[kc]=g_dist;
						}
#endif
#if 0
						if(inc<0)
							a_dist[kc]-=a_dist[kc]>1;
						else if(inc>0)
							a_dist[kc]=g_dist;
#endif
#if 0
						if(reset)
						{
							a_dist[kc]=g_dist;
							dist2symx[kc]=0;
							dist2symy=0;
						}
						else if(dist2symx[kc]>=MAXRUN||dist2symy>=MAXRUN)
						{
							a_dist[kc]-=a_dist[kc]>1;
							dist2symx[kc]=0;
							dist2symy=0;
						}
						rows[0][2+0*NCH*NROWS*NVAL]=dist2symy;
#endif
					}
				}
				else
				{
					if(fwd)
					{
						e=curr-predc;
						e<<=32-src->depth[kc];
						e>>=32-src->depth[kc];
						((int32_t*)HDR_DATA(src))[idx+kc]=e;
					}
					else
					{
						curr=((int32_t*)HDR_DATA(src))[idx+kc]+predc;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						yuv[kc]=curr;
					}
				}
				curr-=offset;
				rows[0][0]=curr;
				rows[0][1]=curr-(int)pred1;
			/*	rows[0][2]=curr-(int)pred2s;*/

				/*update*/
				e=(curr>pred1)-(curr<pred1);
				currw[NPREDS]+=e;/*bias*/
				for(k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];/*coeffs*/

				e=(curr>pred2s)-(curr<pred2s);
				currw2[NPREDS2]+=e;
				for(k=0;k<NPREDS2;++k)
					currw2[k]+=e*preds2[k];
			}
			if(!fwd)
			{
				((int32_t*)HDR_DATA(src))[idx+perm[0]]=yuv[0];
				((int32_t*)HDR_DATA(src))[idx+perm[1]]=yuv[1];
				((int32_t*)HDR_DATA(src))[idx+perm[2]]=yuv[2];
			}
		}
	}
	free(pixels);
	free(weights);
	free(weights2);
}
