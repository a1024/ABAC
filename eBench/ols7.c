#include"ebench.h"
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<immintrin.h>




#if 0
#define L1SH 19
#define PREDLIST\
	PRED(100000, N)\
	PRED(100000, W)\
	PRED( 80000, 3*(N-NN)+NNN)\
	PRED( 80000, 3*(W-WW)+WWW)\
	PRED( 40000, (WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)>>2)\
	PRED( 50000, W+NE-N)\
	PRED(150000, N+W-NW)\
	PRED( 50000, N+NE-NNE)\

#endif
#if 1
#define PREDLIST\
	PRED(100000, N)\
	PRED(100000, W)\
	PRED( 40000, NNN)\
	PRED( 40000, WWW)\
	PRED( 40000, NEEE)\
	PRED( 80000, 3*(N-NN)+NNN)\
	PRED( 80000, 3*(W-WW)+WWW)\
	PRED( 50000, W+NE-N)\
	PRED(150000, N+W-NW)\
	PRED( 50000, N+NE-NNE)\
	PRED( 40000, (WWWW+WWW+NNN+NNEE+NEEE+NEEEE-(N+W))>>2)\

#if 0
	PRED( 40000, (WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)>>2)
#endif
#endif
#if 0
#define L1SH 19
#define PREDLIST\
	PRED(100000, N)\
	PRED(100000, W)\
	PRED( 80000, 3*(N-NN)+NNN)\
	PRED( 80000, 3*(W-WW)+WWW)\
	PRED( 50000, W+NE-N)\
	PRED( 50000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)\
	PRED(150000, N+W-NW)\
	PRED( 50000, N+NE-NNE)\
	PRED( 40000, N+NW-NNW)\
	PRED( 40000, W+NW-NWW)\
	PRED( 40000, NEEE)\
	PRED( 40000, NW)\
	PRED( 40000, NE)\
	PRED( 40000, NN)\
	PRED( 40000, WW)\

#endif
#if 0
#define PREDLIST\
	PRED( 40000, N)\
	PRED( 40000, W)\
	PRED( 40000, 3*(N-NN)+NNN)\
	PRED( 40000, 3*(W-WW)+WWW)\
	PRED( 40000, W+NE-N)\
	PRED(160000, N+W-NW)\
	PRED( 40000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)\
	PRED( 40000, N+NE-NNE)\
	PRED( 40000, W+NW-NWW)\
	PRED( 40000, NEEE)\

#endif
#if 0
#define PREDLIST\
	PRED( 38000, N)\
	PRED( 69000, W)\
	PRED( 41000, 3*(N-NN)+NNN)\
	PRED( 72000, 3*(W-WW)+WWW)\
	PRED( 70000, W+NE-N)\
	PRED( 83000, N+W-NW)\
	PRED(-10000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)\
	PRED( 61000, N+NE-NNE)\
	PRED( 81000, W+NW-NWW)\
	PRED( 18000, NEEE)\

#endif
enum
{
	L1SH=21,
#define PRED(WEIGHT, EXPR) +1
	NPREDS=PREDLIST
#undef  PRED
};
void pred_ols7(Image *src, int fwd)
{
	int amin[4], amax[4];
	int weights[4][NPREDS]={0}, bias[4]={1<<L1SH>>1, 1<<L1SH>>1, 1<<L1SH>>1, 1<<L1SH>>1};
	int invdist;
	int psize;
	int16_t *pixels;
	int kx, ky, kc, idx;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	invdist=((1<<16)+g_dist-1)/g_dist;
	psize=(src->iw+16*2)*(int)sizeof(int16_t[4*4*1]);/*4 padded rows * 4 channels max * {pixels}		FIXME FIXME FIXME FIXME*/
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	FILLMEM((int*)weights, (1<<L1SH)/NPREDS, sizeof(weights), sizeof(int));
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int16_t *rows[4];

		rows[0]=pixels+(((src->iw+(intptr_t)16*2)*((ky-(intptr_t)0+4)%4)+16)*4-1)*1;
		rows[1]=pixels+(((src->iw+(intptr_t)16*2)*((ky-(intptr_t)1+4)%4)+16)*4-1)*1;
		rows[2]=pixels+(((src->iw+(intptr_t)16*2)*((ky-(intptr_t)2+4)%4)+16)*4-1)*1;
		rows[3]=pixels+(((src->iw+(intptr_t)16*2)*((ky-(intptr_t)3+4)%4)+16)*4-1)*1;
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				++rows[0];
				++rows[1];
				++rows[2];
				++rows[3];
				if(!src->depth[kc])
					continue;
				{
					int
						NNNWWW		=rows[3][-3*4*1],
						NNNW		=rows[3][-1*4*1],
						NNN		=rows[3][+0*4*1],
						NNNE		=rows[3][+1*4*1],
						NNNEE		=rows[3][+2*4*1],
						NNNEEE		=rows[3][+3*4*1],
						NNNEEEE		=rows[3][+4*4*1],
						NNWWWW		=rows[2][-4*4*1],
						NNWWW		=rows[2][-3*4*1],
						NNWW		=rows[2][-2*4*1],
						NNW		=rows[2][-1*4*1],
						NN		=rows[2][+0*4*1],
						NNE		=rows[2][+1*4*1],
						NNEE		=rows[2][+2*4*1],
						NNEEE		=rows[2][+3*4*1],
						NNEEEE		=rows[2][+4*4*1],
						NWWWW		=rows[1][-4*4*1],
						NWWW		=rows[1][-3*4*1],
						NWW		=rows[1][-2*4*1],
						NW		=rows[1][-1*4*1],
						N		=rows[1][+0*4*1],
						NE		=rows[1][+1*4*1],
						NEE		=rows[1][+2*4*1],
						NEEE		=rows[1][+3*4*1],
						NEEEE		=rows[1][+4*4*1],
						NEEEEE		=rows[1][+5*4*1],
						NEEEEEE		=rows[1][+6*4*1],
						NEEEEEEE	=rows[1][+7*4*1],
						NEEEEEEEE	=rows[1][+8*4*1],
						WWWWWWWWW	=rows[0][-9*4*1],
						WWWWWWWW	=rows[0][-8*4*1],
						WWWWWWW		=rows[0][-7*4*1],
						WWWWWW		=rows[0][-6*4*1],
						WWWWW		=rows[0][-5*4*1],
						WWWW		=rows[0][-4*4*1],
						WWW		=rows[0][-3*4*1],
						WW		=rows[0][-2*4*1],
						W		=rows[0][-1*4*1];
					int preds[NPREDS], j, *currw, predc, vmax, vmin, p0, k, e, curr;

#if 1
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
					currw=weights[kc];
					predc=bias[kc];
					for(k=0;k<NPREDS;++k)
						predc+=currw[k]*preds[k];
					predc>>=L1SH;
#if 0
					p0-=p0>>31;/*X  deadzone bad with advanced pred*/
#endif
					p0=predc;

					vmax=N; vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(predc, vmin, vmax);

					curr=((int32_t*)HDR_DATA(src))[idx];
					if(g_dist>1)
					{
						if(fwd)
						{
							curr-=predc;
#if 0
							curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));
#endif
							curr=(curr*invdist>>16)-(curr>>31);/*curr/=g_dist*/
							((int32_t*)HDR_DATA(src))[idx]=curr;
						}
						curr=g_dist*curr+predc;
						CLAMP2(curr, amin[kc], amax[kc]);
						if(!fwd)
							((int32_t*)HDR_DATA(src))[idx]=curr;
					}
					else
					{
						if(fwd)
						{
							e=curr-predc;
							e<<=32-src->depth[kc];
							e>>=32-src->depth[kc];
							((int32_t*)HDR_DATA(src))[idx]=e;
						}
						else
						{
							curr+=predc;
							curr<<=32-src->depth[kc];
							curr>>=32-src->depth[kc];
							((int32_t*)HDR_DATA(src))[idx]=curr;
						}
					}
					rows[0][0]=curr;

					/*update*/
					e=(curr>p0)-(curr<p0);/*L1*/
					bias[kc]+=e<<4;
					for(k=0;k<NPREDS;++k)
						currw[k]+=e*preds[k];
				}
			}
		}
	}
	free(pixels);
}

void pred_mixN(Image *src, int fwd)
{
	enum
	{
		MIXPREDS=6,

		SHIFT=20,

		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=2
	};
	ALIGN(16) int32_t coeffs[4][MIXPREDS*3]={0}, bias[4]={1<<SHIFT>>1}, estims[MIXPREDS*3]={0};
	int amin[4], amax[4], rmin[4], rmax[4], invdist, kx, ky, kc, idx;
	int psize;
	int32_t *pixels;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	rmin[0]=amin[0]/g_dist;
	rmin[1]=amin[1]/g_dist;
	rmin[2]=amin[2]/g_dist;
	rmin[3]=amin[3]/g_dist;
	rmax[0]=amax[0]/g_dist,
	rmax[1]=amax[1]/g_dist,
	rmax[2]=amax[2]/g_dist,
	rmax[3]=amax[3]/g_dist,
	invdist=((1<<16)+g_dist-1)/g_dist;
	/*
	cache-friendly layout:
	...
	A(NNN NN N C)WW

	Y(NNN NN N C)W
	U(NNN NN N C)W
	V(NNN NN N C)W
	A(NNN NN N C)W

	Y(NNN NN N C)curr
	...
	*/
	psize=(src->iw+2*XPAD)*(int)sizeof(int32_t[NROWS*NCH*NVAL]);
	pixels=(int32_t*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	FILLMEM((int32_t*)coeffs, (1<<SHIFT)/MIXPREDS, sizeof(coeffs), sizeof(int32_t));
	bias[3]=bias[2]=bias[1]=bias[0];
	memset(pixels, 0, psize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int32_t *rows[4];

		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;/*sub 1 channel for pre-increment*/
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				{
					int32_t
						NNN	=rows[3][0+0*NCH*NROWS*NVAL],
						NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
						NNW	=rows[2][0-1*NCH*NROWS*NVAL],
						NN	=rows[2][0+0*NCH*NROWS*NVAL],
						NNE	=rows[2][0+1*NCH*NROWS*NVAL],
						NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
						NWW	=rows[1][0-2*NCH*NROWS*NVAL],
						NW	=rows[1][0-1*NCH*NROWS*NVAL],
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						NE	=rows[1][0+1*NCH*NROWS*NVAL],
						NEE	=rows[1][0+2*NCH*NROWS*NVAL],
						NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
						NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
						WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
						WWW	=rows[0][0-3*NCH*NROWS*NVAL],
						WW	=rows[0][0-2*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL],
						cN	=rows[1][1+0*NCH*NROWS*NVAL],
						cW	=rows[0][1-1*NCH*NROWS*NVAL];
					int curr, j, o1, o2, pred, e;
					int64_t p1;

#if 1
					(void)NNN;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
					(void)NWW;
					(void)NW;
					(void)N;
					(void)NE;
					(void)NEE;
					(void)NEEE;
					(void)NEEEE;
					(void)WWWW;
					(void)WWW;
					(void)WW;
					(void)W;
					(void)cN;
					(void)cW;
#endif
					curr=((int32_t*)HDR_DATA(src))[idx];
					j=0;
					estims[j++]=W;
					estims[j++]=NE;
					estims[j++]=2*N-NN;
					estims[j++]=N+W-NW;

					o1=kc==2?-2:1;
					o2=kc==0?2:-1;
					estims[j++]=rows[0][0+(o1-1*NCH)*NROWS*NVAL];/*176 MB/s*/
					estims[j++]=rows[0][0+(o2-1*NCH)*NROWS*NVAL];

					p1=bias[kc];
					for(j=0;j<MIXPREDS;++j)
						p1+=(int64_t)coeffs[kc][j]*estims[j];
					p1>>=SHIFT;
					pred=(int)p1;
#if 0
					{
						int vmax=N, vmin=W;

						if(N<W)vmin=N, vmax=W;
						if(vmin>NE)vmin=NE;
						if(vmax<NE)vmax=NE;
						if(vmin>NEEE)vmin=NEEE;
						if(vmax<NEEE)vmax=NEEE;
						CLAMP2(pred, vmin, vmax);
					}
#endif
					if(g_dist>1)
					{
						if(fwd)
						{
							curr-=(int)pred;
#if 0
							curr=curr<0?-(-curr>>1):curr>>1;/*IMG0008 d3  17.9% 33.3 dB  19.19% 34.0 dB*/

							curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));
#endif
							curr=(curr*invdist>>16)-(curr>>31);/*curr/=g_dist*/
							CLAMP2(curr, rmin[kc], rmax[kc]);
							((int32_t*)HDR_DATA(src))[idx]=curr;
						}
						curr=g_dist*curr+(int)pred;
						CLAMP2(curr, amin[kc], amax[kc]);
						if(!fwd)
							((int32_t*)HDR_DATA(src))[idx]=curr;
					}
					else
					{
						if(fwd)
						{
							e=curr-pred;
							e<<=32-src->depth[kc];
							e>>=32-src->depth[kc];
							((int32_t*)HDR_DATA(src))[idx]=e;
						}
						else
						{
							curr+=pred;
							curr<<=32-src->depth[kc];
							curr>>=32-src->depth[kc];
							((int32_t*)HDR_DATA(src))[idx]=curr;
						}
					}
					rows[0][0]=curr;
#if 0
					rows[0][1]=curr-p1;
#endif

					/*150 MB/s  6.65 ms/MB  i5-1145G7*/
#if 0
					{
						__m128i p=_mm_load_si128((__m128i*)estims);
						__m128i c=_mm_load_si128((__m128i*)coeffs[kc]);
						p=_mm_sign_epi32(p, _mm_set1_epi32(curr-p1));
						c=_mm_add_epi32(c, p);
						_mm_store_si128((__m128i*)coeffs[kc], c);
					}
#endif

					/*147 MB/s  6.77 ms/MB  i5-1145G7*/
#if 1
					e=(curr>(int)p1)-(curr<(int)p1);
#if 0
					e=((curr-p1)>>31)-((p1-curr)>>31);
					e=curr-p1; CLAMP2(e, -1, 1);/*jump?*/
#endif
					bias[kc]+=e;
					for(j=0;j<MIXPREDS;++j)
						coeffs[kc][j]+=(int16_t)((int16_t)e*(int16_t)estims[j]);/*casts prevent pmulld*/
#endif
				}
			}
		}
	}
	_mm_free(pixels);
}
void pred_mixR(Image *src, int fwd)
{
#define ESTIMLIST\
	ESTIM(W)\
	ESTIM(2*N-NN)\
	ESTIM(N+W-NW)\
	ESTIM(NE)\

#if 0
#define ESTIMLIST\
	ESTIM(N)\
	ESTIM(W)\
	ESTIM(2*N-NN)\
	ESTIM(2*W-WW)\
	ESTIM(3*(N-NN)+NNN)\
	ESTIM(3*(W-WW)+WWW)\
	ESTIM(N+W-NW)\
	ESTIM(W+NE-N)\

#endif
	enum
	{
#if 0
		RCTBITS2=23,
#endif
		SHIFT=21,
#define ESTIM(E) +1
		MIXPREDS=ESTIMLIST,
#undef  ESTIM

		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=2
	};
	ALIGN(16) int32_t coeffs[4][MIXPREDS*3]={0}, bias[4]={0}, estims[MIXPREDS*3]={0};
	int psize;
	int32_t *pixels, *rows[4];
	int amin[4], amax[4], rmin[4], rmax[4], invdist, ky, kx, kc, idx;
	const uint8_t *perm;
	int uc0, vc0, vc1;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	rmin[0]=amin[0]/g_dist;
	rmin[1]=amin[1]/g_dist;
	rmin[2]=amin[2]/g_dist;
	rmin[3]=amin[3]/g_dist;
	rmax[0]=amax[0]/g_dist;
	rmax[1]=amax[1]/g_dist;
	rmax[2]=amax[2]/g_dist;
	rmax[3]=amax[3]/g_dist;
	invdist=((1<<16)+g_dist-1)/g_dist;
	uc0=vc0=vc1=0;
	psize=(src->iw+2*XPAD)*(int)sizeof(int32_t[NROWS*NCH*NVAL]);
	pixels=(int32_t*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	if(fwd)
		src->rct=crct2_analysis(src);
	perm=crct2_unpack(src->rct, &uc0, &vc0, &vc1);
#if 0
	if(fwd)
		src->rct=crct_analysis2(src);
	combination=rct_combinations[src->rct];
	yidx=combination[II_PERM_Y];
	uidx=combination[II_PERM_U];
	vidx=combination[II_PERM_V];
#endif
	FILLMEM((int32_t*)coeffs, (1<<SHIFT)/MIXPREDS, sizeof(coeffs), sizeof(int32_t));
	bias[3]=bias[2]=bias[1]=bias[0]=1<<SHIFT>>1;
	memset(pixels, 0, psize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;/*sub 1 channel for pre-increment*/
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			int yuv[3];

			if(fwd)
			{
				yuv[0]=((int32_t*)(src+1))[idx+perm[0]];
				yuv[1]=((int32_t*)(src+1))[idx+perm[1]];
				yuv[2]=((int32_t*)(src+1))[idx+perm[2]];
			}
			for(kc=0;kc<4;++kc)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(kc==3)
					continue;
				{
					int32_t
						NNN	=rows[3][0+0*NCH*NROWS*NVAL],
						NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
						NNW	=rows[2][0-1*NCH*NROWS*NVAL],
						NN	=rows[2][0+0*NCH*NROWS*NVAL],
						NNE	=rows[2][0+1*NCH*NROWS*NVAL],
						NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
						NWW	=rows[1][0-2*NCH*NROWS*NVAL],
						NW	=rows[1][0-1*NCH*NROWS*NVAL],
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						NE	=rows[1][0+1*NCH*NROWS*NVAL],
						NEE	=rows[1][0+2*NCH*NROWS*NVAL],
						NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
						NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
						WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
						WWW	=rows[0][0-3*NCH*NROWS*NVAL],
						WW	=rows[0][0-2*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL],
						cN	=rows[1][1+0*NCH*NROWS*NVAL],
						cW	=rows[0][1-1*NCH*NROWS*NVAL];
					int offset, curr, j, pred, e;
					int64_t p1;

#if 1
					(void)NNN;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
					(void)NWW;
					(void)NW;
					(void)N;
					(void)NE;
					(void)NEE;
					(void)NEEE;
					(void)NEEEE;
					(void)WWWW;
					(void)WWW;
					(void)WW;
					(void)W;
					(void)cN;
					(void)cW;
#endif
					if(!fwd)
						yuv[kc]=((int32_t*)(src+1))[idx+kc];
					curr=yuv[kc];
#define ESTIM(E) estims[j++]=E;
					j=0;
					ESTIMLIST;
#undef  ESTIM
#if 0
					estims[j++]=W;
					estims[j++]=NE;
					estims[j++]=2*N-NN;
					estims[j++]=N+W-NW;
#endif

					p1=bias[kc];
					for(j=0;j<MIXPREDS;++j)
						p1+=(int64_t)coeffs[kc][j]*estims[j];
					p1>>=SHIFT;
					pred=(int)p1;
				/*	p2=(int)p1;*/

#if 0
					int vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(pred, vmin, vmax);
#endif
					offset=0;
					if(kc==1)offset=uc0*yuv[0];
					if(kc==2)offset=vc0*yuv[0]+vc1*yuv[1];
					offset>>=RCTBITS;
					pred+=offset;
#if 0
					if(kc==1)pred+=(int)((uc0*errors[0])>>RCTBITS);
#if 0
					if(kc==2)
					{
						eblend=(int)(errors[1]+((((int64_t)errors[0]-errors[1])*vc0+((int64_t)1<<RCTBITS2>>1))>>RCTBITS2));
						pred+=(int)((eblend*vc1+((int64_t)1<<RCTBITS2>>1))>>RCTBITS2);
					}
#endif
					if(kc==2)pred+=(int)((vc0*errors[0]+vc1*errors[1])>>RCTBITS);
#endif
					CLAMP2(pred, amin[kc], amax[kc]);

					if(g_dist>1)
					{
						if(fwd)
						{
							curr-=(int)pred;
#if 0
							curr=curr<0?-(-curr>>1):curr>>1;/*IMG0008 d3  17.9% 33.3 dB  19.19% 34.0 dB*/

							curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));
#endif
							curr=(curr*invdist>>16)-(curr>>31);/*curr/=g_dist*/
							CLAMP2(curr, rmin[kc], rmax[kc]);
							((int32_t*)(src+1))[idx+kc]=curr;
						}
						curr=g_dist*curr+(int)pred;
						CLAMP2(curr, amin[kc], amax[kc]);
						yuv[kc]=curr;
					}
					else
					{
						if(fwd)
						{
							e=curr-pred;
							e<<=32-src->depth[kc];
							e>>=32-src->depth[kc];
							((int32_t*)(src+1))[idx+kc]=e;
						}
						else
						{
							curr+=pred;
							curr<<=32-src->depth[kc];
							curr>>=32-src->depth[kc];
							yuv[kc]=curr;
						}
					}
					curr-=offset;
					rows[0][0]=curr;
#if 0
					errors[kc]=curr-p2;
#endif

					/*147 MB/s  6.77 ms/MB  i5-1145G7*/
#if 1
					e=(curr>(int)p1)-(curr<(int)p1);/*cmp setgt setlt sub*/
#if 0
					e=((curr-p1)>>31)-((p1-curr)>>31);/*sub sar sub sar sub*/
					e=curr-p1; CLAMP2(e, -1, 1);/*sub cmp cmovlt cmovgt*/
#endif
					bias[kc]+=e;
					for(j=0;j<MIXPREDS;++j)
						coeffs[kc][j]+=(int16_t)((int16_t)e*(int16_t)estims[j]);
#endif
				}
			}
#if 0
			uc0+=((errors[1]>0)-(errors[1]<0))*errors[0];

#if 0
			vc0+=((int64_t)((errors[2]>0)-(errors[2]<0))*(errors[0]-errors[1])*vc1+((int64_t)1<<RCTBITS2>>1))>>RCTBITS2;
			vc1+=((errors[2]>0)-(errors[2]<0))*eblend;
#endif
			vc0+=((errors[2]>0)-(errors[2]<0))*errors[0];
			vc1+=((errors[2]>0)-(errors[2]<0))*errors[1];
#endif
			if(!fwd)
			{
				((int32_t*)(src+1))[idx+perm[0]]=yuv[0];
				((int32_t*)(src+1))[idx+perm[1]]=yuv[1];
				((int32_t*)(src+1))[idx+perm[2]]=yuv[2];
			}
		}
	}
	_mm_free(pixels);
}

enum
{
	LSTM_NESTIMS=6,
	LSTM_NFEATURES=14,
	LSTM_SIZE=7
};
typedef struct _LSTMState
{
	float W[4*LSTM_SIZE][LSTM_NFEATURES+LSTM_SIZE];
	float bias[4*LSTM_SIZE];
	float cell[LSTM_SIZE];
	float hidden[LSTM_SIZE];
	float Wmix[LSTM_NESTIMS][LSTM_SIZE];
	float bmix[LSTM_NESTIMS];
} LSTMState;
void pred_mixNC(Image *src, int fwd)
{
	enum
	{
		MIXPREDS=4,

		SHIFT=20,
		SHIFT2=16,

		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=2
	};
	ALIGN(16) LSTMState lstms[3]={0};
	float estims[LSTM_NESTIMS]={0}, xh[LSTM_NFEATURES+LSTM_SIZE]={0}, epreds[3][LSTM_NESTIMS]={0}, eprev[3]={0};
	int amin[4], amax[4], rmin[4], rmax[4], invdist, ky, kx, kc, idx;
	int psize;
	int32_t *pixels;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	rmin[0]=amin[0]/g_dist;
	rmin[1]=amin[1]/g_dist;
	rmin[2]=amin[2]/g_dist;
	rmin[3]=amin[3]/g_dist;
	rmax[0]=amax[0]/g_dist;
	rmax[1]=amax[1]/g_dist;
	rmax[2]=amax[2]/g_dist;
	rmax[3]=amax[3]/g_dist;
	invdist=((1<<16)+g_dist-1)/g_dist;
	psize=(src->iw+2*XPAD)*(int)sizeof(int32_t[NROWS*NCH*NVAL]);
	pixels=(int32_t*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int32_t *rows[4];

		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc]||kc==3)
					continue;
				{
					int32_t
						NNN	=rows[3][0+0*NCH*NROWS*NVAL],
						NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
						NNW	=rows[2][0-1*NCH*NROWS*NVAL],
						NN	=rows[2][0+0*NCH*NROWS*NVAL],
						NNE	=rows[2][0+1*NCH*NROWS*NVAL],
						NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
						NWW	=rows[1][0-2*NCH*NROWS*NVAL],
						NW	=rows[1][0-1*NCH*NROWS*NVAL],
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						NE	=rows[1][0+1*NCH*NROWS*NVAL],
						NEE	=rows[1][0+2*NCH*NROWS*NVAL],
						NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
						NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
						WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
						WWW	=rows[0][0-3*NCH*NROWS*NVAL],
						WW	=rows[0][0-2*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL];
					int curr, j, k, k2, vmax, vmin, pred, e;
					float yhat, fmax, sum, z[4*LSTM_SIZE]={0}, c0[LSTM_SIZE], tanhc[LSTM_SIZE], logits[LSTM_NESTIMS]={0}, w[LSTM_NESTIMS]={0};
					LSTMState *lstm;

					(void)NNN;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
					(void)NWW;
					(void)NW;
					(void)N;
					(void)NE;
					(void)NEE;
					(void)NEEE;
					(void)NEEEE;
					(void)WWWW;
					(void)WWW;
					(void)WW;
					(void)W;
					curr=((int32_t*)HDR_DATA(src))[idx];
					lstm=lstms+kc;
					j=0;
					estims[j++]=(float)(W);
					estims[j++]=(float)(NE);
					estims[j++]=(float)(2*N-NN);
					estims[j++]=(float)(N+W-NW);
					estims[j++]=(float)(W+NE-N);
					estims[j++]=(float)(N+NE-NNE);
					j=0;
					xh[j++]=eprev[kc];
					xh[j++]=(float)abs(N-NW);
					xh[j++]=(float)abs(W-NW);
					xh[j++]=epreds[kc][0];
					xh[j++]=epreds[kc][1];
					xh[j++]=epreds[kc][2];
					xh[j++]=epreds[kc][3];
					xh[j++]=epreds[kc][4];
					xh[j++]=epreds[kc][5];
					xh[j++]=estims[0]-estims[1];
					xh[j++]=estims[0]-estims[2];
					xh[j++]=estims[0]-estims[3];
					xh[j++]=estims[0]-estims[4];
					xh[j++]=estims[0]-estims[5];
					for(k=0;k<LSTM_SIZE;++k)
						xh[j++]=lstm->hidden[k];
					for(k=0;k<LSTM_NFEATURES+LSTM_SIZE;++k)
					{
						float feature=xh[k];
						for(k2=0;k2<4*LSTM_SIZE;++k2)
							z[k2]+=lstm->W[k2][k]*feature;
					}
					for(k2=0;k2<4*LSTM_SIZE;++k2)
						z[k2]+=lstm->bias[k2];
					for(k2=0;k2<3*LSTM_SIZE;++k2)
						z[k2]=1/(1+(float)exp(-z[k2]));
					for(k2=3*LSTM_SIZE;k2<4*LSTM_SIZE;++k2)
						z[k2]=(float)tanh(z[k2]);
					memcpy(c0, lstm->cell, sizeof(c0));
					for(k2=0;k2<LSTM_SIZE;++k2)
						lstm->cell[k2]=z[k2+1*LSTM_SIZE]*lstm->cell[k2]+z[k2+0*LSTM_SIZE]*z[k2+3*LSTM_SIZE];
					for(k2=0;k2<LSTM_SIZE;++k2)
					{
						tanhc[k2]=(float)tanh(lstm->cell[k2]);
						lstm->hidden[k2]=z[k2+2*LSTM_SIZE]*tanhc[k2];
					}
					for(k=0;k<LSTM_SIZE;++k)
					{
						for(k2=0;k2<LSTM_NESTIMS;++k2)
							logits[k2]+=lstm->Wmix[k2][k]*lstm->hidden[k];
					}
					for(k2=0;k2<LSTM_NESTIMS;++k2)
						logits[k2]+=lstm->bmix[k2];
					fmax=logits[0];
					sum=0;
					for(k2=1;k2<LSTM_NESTIMS;++k2)
					{
						if(fmax<logits[k2])
							fmax=logits[k2];
					}
					for(k2=0;k2<LSTM_NESTIMS;++k2)
						sum+=w[k2]=(float)exp(logits[k2]-fmax);
					sum=1.f/sum;
					for(k2=0;k2<LSTM_NESTIMS;++k2)
						w[k2]*=sum;
					yhat=0;
					for(k2=0;k2<LSTM_NESTIMS;++k2)
						yhat+=w[k2]*estims[k2];
					pred=(int)round(yhat);

					vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(pred, vmin, vmax);
					if(g_dist>1)
					{
						if(fwd)
						{
							curr-=(int)pred;
							curr=(curr*invdist>>16)-(curr>>31);/*curr/=g_dist*/
							CLAMP2(curr, rmin[kc], rmax[kc]);
							((int32_t*)HDR_DATA(src))[idx]=curr;
						}
						curr=g_dist*curr+(int)pred;
						CLAMP2(curr, amin[kc], amax[kc]);
						if(!fwd)
							((int32_t*)HDR_DATA(src))[idx]=curr;
					}
					else
					{
						if(fwd)
						{
							e=curr-pred;
							e<<=32-src->depth[kc];
							e>>=32-src->depth[kc];
							((int32_t*)HDR_DATA(src))[idx]=e;
						}
						else
						{
							curr+=pred;
							curr<<=32-src->depth[kc];
							curr>>=32-src->depth[kc];
							((int32_t*)HDR_DATA(src))[idx]=curr;
						}
					}
					rows[0][0]=curr;
					eprev[kc]=(float)curr-yhat;
					/*
					fwd:  estims p[4]; features x[10]
					i = sigmoid(Wi*x + Ui*hprev + bi)
					f = sigmoid(Wf*x + Uf*hprev + bf)
					o = sigmoid(Wo*x + Uo*hprev + bo)
					g =    tanh(Wg*x + Ug*hprev + bg)
					c = f.*cprev + i.*g
					h = o.*tanh(c)
					logits = Wmix*h + bmix
					w = softmax(logits)
					yhat = w.p
					L = abs(y - yhat)
					
					rules:
					sigmoid' = sigmoid*(1-sigmoid)
					tanh' = 1-tanh^2

					backprop:
					dL/dyhat = -sgn(y-yhat) = e
					dL/dw[k] = e*p[k]
					s = e*p.w = e*yhat
					dL/dlogits[k] = w[k]*(dL/dw[k] - s) = w[k]*e*(p[k] - yhat)
					dL/Wmix = outer(dL/dlogits, h)
					dL/bmix = dL/dlogits
					dL/dh = WmixT * dL/dlogits
					dL/dc = dL/dh .* o .* (1-tanh(c)^2) + dL_future/dc    (BPTT)
					dL/di = dL/dc .* g
					dL/df = dL/dc .* cprev
					dL/do = dL/dh .* tanh(c)
					dL/dg = i .* dL/dc
					dL/dcprev = dL/dc .* f    (BPTT)
					
					dai = dL/di .* i .* (1-i)	= dL/dc .* g      .* i .* (1-i)
					daf = dL/df .* f .* (1-f)	= dL/dc .* cprev  .* f .* (1-f)
					dao = dL/do .* o .* (1-o)	= dL/dh .* tanh(c) .* o .* (1-o)
					dag = dL/dg .* (1-g.*g)		= i .* dL/dc .* (1-g.*g)
					dL/bi += dai
					dL/bf += daf
					dL/bo += dao
					dL/bg += dag
					dL/dhprev = UiT*dL/bi + UfT*dL/bf + UoT*dL/bo + UgT*dL/bg
					dL/dWi = outer(dai, x)
					dL/dWf = outer(daf, x)
					dL/dWo = outer(dao, x)
					dL/dWg = outer(dag, x)
					dL/dUi = outer(dai, hprev)
					dL/dUf = outer(daf, hprev)
					dL/dUo = outer(dao, hprev)
					dL/dUg = outer(dag, hprev)

					dL/dx = WiT*dai + WfT*daf + WoT*dao + WgT*dag  (unused)
					*/
					if(pred!=curr)
					{
						const float LR=0.02f;
					/*	float s=e*yhat;*/
						float dlogits[LSTM_NESTIMS]={0};
						float dh[LSTM_SIZE]={0};
						float dc[LSTM_SIZE]={0};
						float *i=z+0*LSTM_SIZE;
						float *f=z+1*LSTM_SIZE;
						float *o=z+2*LSTM_SIZE;
						float *g=z+3*LSTM_SIZE;
						float grad[4*LSTM_SIZE];
						float *dai=grad+0*LSTM_SIZE;
						float *daf=grad+1*LSTM_SIZE;
						float *dao=grad+2*LSTM_SIZE;
						float *dag=grad+3*LSTM_SIZE;
						
						e=-((curr>pred)-(curr<pred));/*-sgn(y-yhat)*/
						for(k=0;k<LSTM_NESTIMS;++k)
							dlogits[k]=w[k]*e*(estims[k]-yhat);
						for(k=0;k<LSTM_SIZE;++k)
						{
							for(k2=0;k2<LSTM_NESTIMS;++k2)
								dh[k]+=lstm->Wmix[k2][k]*dlogits[k2];
						}
						for(k=0;k<LSTM_NESTIMS;++k)
						{
							lstm->bmix[k]-=LR*dlogits[k];/*dL/bmix = dL/dlogits*/
							for(k2=0;k2<LSTM_SIZE;++k2)/*dL/Wmix = outer(dL/dlogits, h)*/
								lstm->Wmix[k][k2]-=LR*dlogits[k]*lstm->hidden[k2];
						}
						for(k=0;k<LSTM_SIZE;++k)
						{
							dc[k]=dh[k]*o[k]*(1-tanhc[k]*tanhc[k]);
							dai[k]=dc[k]*g[k]*i[k]*(1-i[k]);
							daf[k]=dc[k]*c0[k]*f[k]*(1-f[k]);
							dao[k]=dh[k]*tanhc[k]*o[k]*(1-o[k]);
							dag[k]=i[k]*dc[k]*(1-g[k]*g[k]);
						}
						for(k=0;k<4*LSTM_SIZE;++k)
							lstm->bias[k]-=LR*grad[k];
						for(k=0;k<4*LSTM_SIZE;++k)
						{
							for(k2=0;k2<LSTM_NFEATURES+LSTM_SIZE;++k2)
								lstm->W[k][k2]-=LR*grad[k]*xh[k2];
						}
					}
					for(k=0;k<LSTM_NESTIMS;++k)
						epreds[kc][k]+=(fabsf(curr-estims[k])-epreds[kc][k])*(1.f/32);
				}
			}
		}
	}
	_mm_free(pixels);
}
void pred_bestN(Image *src, int fwd)
{
	enum
	{
		MIXPREDS=4,

		SHIFT=18+6,

		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=1
	};
	ALIGN(16) int32_t coeffs[4][MIXPREDS]={0}, bias[4]={1<<SHIFT>>1}, estims[MIXPREDS]={0};
	int amin[4], amax[4], rmin[4], rmax[4], invdist, kx, ky, kc, idx;
	int psize;
	int32_t *pixels;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	rmin[0]=amin[0]/g_dist;
	rmin[1]=amin[1]/g_dist;
	rmin[2]=amin[2]/g_dist;
	rmin[3]=amin[3]/g_dist;
	rmax[0]=amax[0]/g_dist;
	rmax[1]=amax[1]/g_dist;
	rmax[2]=amax[2]/g_dist;
	rmax[3]=amax[3]/g_dist;
	invdist=((1<<16)+g_dist-1)/g_dist;

	psize=(src->iw+2*XPAD)*(int)sizeof(int32_t[NROWS*NCH*NVAL]);
	pixels=(int32_t*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	FILLMEM((int32_t*)coeffs, (1<<SHIFT)/MIXPREDS, sizeof(coeffs), sizeof(int32_t));
	bias[3]=bias[2]=bias[1]=bias[0];
	memset(pixels, 0, psize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int32_t *rows[4];

		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;/*sub 1 channel for pre-increment*/
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				{
					int32_t
						NNN	=rows[3][0+0*NCH*NROWS*NVAL],
						NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
						NNW	=rows[2][0-1*NCH*NROWS*NVAL],
						NN	=rows[2][0+0*NCH*NROWS*NVAL],
						NNE	=rows[2][0+1*NCH*NROWS*NVAL],
						NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
						NWW	=rows[1][0-2*NCH*NROWS*NVAL],
						NW	=rows[1][0-1*NCH*NROWS*NVAL],
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						NE	=rows[1][0+1*NCH*NROWS*NVAL],
						NEE	=rows[1][0+2*NCH*NROWS*NVAL],
						NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
						NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
						WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
						WWW	=rows[0][0-3*NCH*NROWS*NVAL],
						WW	=rows[0][0-2*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL];
					int curr, j, p1, pred, vmax, vmin, e;

					(void)NNN;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
					(void)NWW;
					(void)NW;
					(void)N;
					(void)NE;
					(void)NEE;
					(void)NEEE;
					(void)NEEEE;
					(void)WWWW;
					(void)WWW;
					(void)WW;
					(void)W;
					curr=((int32_t*)HDR_DATA(src))[idx];
					j=0;
					estims[j++]=W;
					estims[j++]=NE;
					estims[j++]=2*N-NN;
					estims[j++]=N+W-NW;

					p1=(int)((bias[kc]
						+(int64_t)coeffs[kc][0]*estims[0]
						+(int64_t)coeffs[kc][1]*estims[1]
						+(int64_t)coeffs[kc][2]*estims[2]
						+(int64_t)coeffs[kc][3]*estims[3]
					)>>SHIFT);
					pred=p1;
					
					vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(pred, vmin, vmax);
					if(g_dist>1)
					{
						if(fwd)
						{
							curr-=(int)pred;
							curr=(curr*invdist>>16)-(curr>>31);
							CLAMP2(curr, rmin[kc], rmax[kc]);
							((int32_t*)HDR_DATA(src))[idx]=curr;
						}
						curr=g_dist*curr+(int)pred;
						CLAMP2(curr, amin[kc], amax[kc]);
						if(!fwd)
							((int32_t*)HDR_DATA(src))[idx]=curr;
					}
					else
					{
						if(fwd)
						{
							e=curr-pred;
							e<<=32-src->depth[kc];
							e>>=32-src->depth[kc];
							((int32_t*)HDR_DATA(src))[idx]=e;
						}
						else
						{
							curr+=pred;
							curr<<=32-src->depth[kc];
							curr>>=32-src->depth[kc];
							((int32_t*)HDR_DATA(src))[idx]=curr;
						}
					}
					rows[0][0]=curr;
					e=(curr>p1)-(curr<p1);
					bias[kc]+=e;
					coeffs[kc][0]+=(int16_t)((int16_t)e*(int16_t)estims[0]);
					coeffs[kc][1]+=(int16_t)((int16_t)e*(int16_t)estims[1]);
					coeffs[kc][2]+=(int16_t)((int16_t)e*(int16_t)estims[2]);
					coeffs[kc][3]+=(int16_t)((int16_t)e*(int16_t)estims[3]);
				}
			}
		}
	}
	_mm_free(pixels);
}
void pred_rls(Image *src, int fwd)
{
	enum
	{
		MIXPREDS=4,
		SHIFT=18,

		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=1
	};
	ALIGN(32) double params[4][MIXPREDS]={0}, invcov[4][MIXPREDS*MIXPREDS]={0};
	int estims[MIXPREDS]={0};
	int amin[4], amax[4], invdist, kx, ky, kc, idx;
	int psize;
	int16_t *pixels;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	invdist=((1<<16)+g_dist-1)/g_dist;
	psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	invcov[0][(MIXPREDS+1)*0]=1;
	invcov[0][(MIXPREDS+1)*1]=1;
	invcov[0][(MIXPREDS+1)*2]=1;
	invcov[0][(MIXPREDS+1)*3]=1;
	memcpy(invcov[1], invcov[0], sizeof(double[MIXPREDS*MIXPREDS]));
	memcpy(invcov[2], invcov[0], sizeof(double[MIXPREDS*MIXPREDS]));
	memcpy(invcov[3], invcov[0], sizeof(double[MIXPREDS*MIXPREDS]));
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int16_t *rows[4];

		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				{
					int16_t
						NNN	=rows[3][0+0*NCH*NROWS*NVAL],
						NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
						NNW	=rows[2][0-1*NCH*NROWS*NVAL],
						NN	=rows[2][0+0*NCH*NROWS*NVAL],
						NNE	=rows[2][0+1*NCH*NROWS*NVAL],
						NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
						NWW	=rows[1][0-2*NCH*NROWS*NVAL],
						NW	=rows[1][0-1*NCH*NROWS*NVAL],
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						NE	=rows[1][0+1*NCH*NROWS*NVAL],
						NEE	=rows[1][0+2*NCH*NROWS*NVAL],
						NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
						NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
						WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
						WWW	=rows[0][0-3*NCH*NROWS*NVAL],
						WW	=rows[0][0-2*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL];
					int curr, vmax, vmin, j, pred;
					double fpred, vec[MIXPREDS], norm, gain[MIXPREDS], e2;

					(void)NNN;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
					(void)NWW;
					(void)NW;
					(void)N;
					(void)NE;
					(void)NEE;
					(void)NEEE;
					(void)NEEEE;
					(void)WWWW;
					(void)WWW;
					(void)WW;
					(void)W;
					curr=((int32_t*)HDR_DATA(src))[idx];
					/*
					|mix 4
					|		NN
					|	NW	N	NE
					|	W	?		198 MB/s
					*/
					j=0;
					estims[j++]=W;
					estims[j++]=N+W-NW;
					estims[j++]=2*N-NN;
					estims[j++]=NE;
					fpred=
						+params[kc][0]*estims[0]
						+params[kc][1]*estims[1]
						+params[kc][2]*estims[2]
						+params[kc][3]*estims[3]
					;
					pred=(int)ROUND64(fpred);
					vmax=N; vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(pred, vmin, vmax);
					if(g_dist>1)
					{
						if(fwd)
						{
							curr-=(int)pred;
#if 0
							curr=curr<0?-(-curr>>1):curr>>1;/*IMG0008 d3  17.9% 33.3 dB  19.19% 34.0 dB*/

							curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));
#endif
							curr=(curr*invdist>>16)-(curr>>31);/*curr/=g_dist*/
							((int32_t*)HDR_DATA(src))[idx]=curr;
						}
						curr=g_dist*curr+(int)pred;
						CLAMP2(curr, amin[kc], amax[kc]);
						if(!fwd)
							((int32_t*)HDR_DATA(src))[idx]=curr;
					}
					else
					{
						if(fwd)
						{
							int e=curr-pred;
							e<<=32-src->depth[kc];
							e>>=32-src->depth[kc];
							((int32_t*)HDR_DATA(src))[idx]=e;
						}
						else
						{
							curr+=pred;
							curr<<=32-src->depth[kc];
							curr>>=32-src->depth[kc];
							((int32_t*)HDR_DATA(src))[idx]=curr;
						}
					}
					rows[0][0]=curr;

					vec[0]=
						+invcov[kc][MIXPREDS*0+0]*estims[0]
						+invcov[kc][MIXPREDS*0+1]*estims[1]
						+invcov[kc][MIXPREDS*0+2]*estims[2]
						+invcov[kc][MIXPREDS*0+3]*estims[3]
					;
					vec[1]=
						+invcov[kc][MIXPREDS*1+0]*estims[0]
						+invcov[kc][MIXPREDS*1+1]*estims[1]
						+invcov[kc][MIXPREDS*1+2]*estims[2]
						+invcov[kc][MIXPREDS*1+3]*estims[3]
					;
					vec[2]=
						+invcov[kc][MIXPREDS*2+0]*estims[0]
						+invcov[kc][MIXPREDS*2+1]*estims[1]
						+invcov[kc][MIXPREDS*2+2]*estims[2]
						+invcov[kc][MIXPREDS*2+3]*estims[3]
					;
					vec[3]=
						+invcov[kc][MIXPREDS*3+0]*estims[0]
						+invcov[kc][MIXPREDS*3+1]*estims[1]
						+invcov[kc][MIXPREDS*3+2]*estims[2]
						+invcov[kc][MIXPREDS*3+3]*estims[3]
					;
					norm=1/(
						+estims[0]*vec[0]
						+estims[1]*vec[1]
						+estims[2]*vec[2]
						+estims[3]*vec[3]
						+1
					);
					gain[0]=vec[0]*norm;
					gain[1]=vec[1]*norm;
					gain[2]=vec[2]*norm;
					gain[3]=vec[3]*norm;
					e2=curr-fpred;
					params[kc][0]+=gain[0]*e2;
					params[kc][1]+=gain[1]*e2;
					params[kc][2]+=gain[2]*e2;
					params[kc][3]+=gain[3]*e2;
					invcov[kc][MIXPREDS*0+0]-=gain[0]*vec[0];
					invcov[kc][MIXPREDS*0+1]-=gain[0]*vec[1];
					invcov[kc][MIXPREDS*0+2]-=gain[0]*vec[2];
					invcov[kc][MIXPREDS*0+3]-=gain[0]*vec[3];

					invcov[kc][MIXPREDS*1+0]-=gain[1]*vec[0];
					invcov[kc][MIXPREDS*1+1]-=gain[1]*vec[1];
					invcov[kc][MIXPREDS*1+2]-=gain[1]*vec[2];
					invcov[kc][MIXPREDS*1+3]-=gain[1]*vec[3];

					invcov[kc][MIXPREDS*2+0]-=gain[2]*vec[0];
					invcov[kc][MIXPREDS*2+1]-=gain[2]*vec[1];
					invcov[kc][MIXPREDS*2+2]-=gain[2]*vec[2];
					invcov[kc][MIXPREDS*2+3]-=gain[2]*vec[3];

					invcov[kc][MIXPREDS*3+0]-=gain[3]*vec[0];
					invcov[kc][MIXPREDS*3+1]-=gain[3]*vec[1];
					invcov[kc][MIXPREDS*3+2]-=gain[3]*vec[2];
					invcov[kc][MIXPREDS*3+3]-=gain[3]*vec[3];
				}
			}
		}
	}
	_mm_free(pixels);
}

enum
{
	HISTMATCHNORM=8192
};
static void normalize(int32_t *hist, int nlevels, int res)
{
	int ks, c;

	for(ks=0, c=0;ks<nlevels;++ks)
	{
		int f=hist[ks];
		hist[ks]=(int)((int64_t)c*HISTMATCHNORM/res);
		c+=f;
	}
}
void ct_histmatch(Image *image, int fwd)
{
	intptr_t k, len, res;
	int hsize;
	int32_t *hist, *CDF_G, *CDF1, *CDF2, *invCDF_R, *invCDF_G, *invCDF_B, h2[256];
	int ks;

	if(!fwd)
		return;
	hsize=(int)sizeof(int32_t[256+512+512+HISTMATCHNORM*3]);
	hist=(int32_t*)malloc(hsize);
	if(!hist)
	{
		CRASH("Alloc error");
		return;
	}
	CDF_G=hist;
	CDF1=CDF_G+256;
	CDF2=CDF1+512;
	invCDF_R=CDF2+512;
	invCDF_G=invCDF_R+HISTMATCHNORM;
	invCDF_B=invCDF_G+HISTMATCHNORM;
	memset(hist, 0, hsize);
	res=(intptr_t)image->iw*image->ih;
	for(k=0, len=4*res;k<len;k+=4)
	{
		int r=((int32_t*)HDR_DATA(image))[k+0];
		int g=((int32_t*)HDR_DATA(image))[k+1];
		int b=((int32_t*)HDR_DATA(image))[k+2];
		++CDF_G   [(uint8_t)(g+128)];
		++invCDF_R[(uint8_t)(r+128)];
		++invCDF_B[(uint8_t)(b+128)];
	}
	normalize(CDF_G, 256, (int)res);
	memcpy(h2, invCDF_R, sizeof(h2));
	normalize(h2, 256, (int)res);
	for(ks=0;;)
	{
		int x1=h2[ks], x2=ks<255?h2[ks+1]:HISTMATCHNORM;
		int mid=(x1+x2+1)>>1;
		for(k=x1;k<mid;++k)
			invCDF_R[k]=ks;
		++ks;
		if(ks>=256)
		{
			for(k=mid;k<x2;++k)
				invCDF_R[k]=255;
			break;
		}
		for(k=mid;k<x2;++k)
			invCDF_R[k]=ks;
	}
	memcpy(h2, invCDF_B, sizeof(h2));
	normalize(h2, 256, (int)res);
	for(ks=0;;)
	{
		int x1=h2[ks], x2=ks<255?h2[ks+1]:HISTMATCHNORM;
		int mid=(x1+x2+1)>>1;
		for(k=x1;k<mid;++k)
			invCDF_B[k]=ks;
		++ks;
		if(ks>=256)
		{
			for(k=mid;k<x2;++k)
				invCDF_B[k]=255;
			break;
		}
		for(k=mid;k<x2;++k)
			invCDF_B[k]=ks;
	}
	for(k=0, len=(intptr_t)4*image->iw*image->ih;k<len;k+=4)
	{
		int r=((int32_t*)HDR_DATA(image))[k+0];
		int g=((int32_t*)HDR_DATA(image))[k+1];
		int b=((int32_t*)HDR_DATA(image))[k+2];
		int ug=CDF_G[(uint8_t)(g+128)];
		r-=invCDF_R[ug]-128;
		b-=invCDF_B[ug]-128;
		((int32_t*)HDR_DATA(image))[k+0]=g;
		((int32_t*)HDR_DATA(image))[k+1]=b;
		((int32_t*)HDR_DATA(image))[k+2]=r;
	}
	free(hist);
	image->depth[0]=8;
	image->depth[1]=9;
	image->depth[2]=9;
}
int crct_analysis(Image *src)
{
	int64_t counters[OCH_COUNT]={0}, minerr;
	int prev[OCH_COUNT]={0}, rct, kt;
	intptr_t k, len;

	for(k=0, len=(intptr_t)src->iw*src->ih*4;k<len;k+=4)
	{
		int
			r=((int32_t*)HDR_DATA(src))[k+0]<<2,
			g=((int32_t*)HDR_DATA(src))[k+1]<<2,
			b=((int32_t*)HDR_DATA(src))[k+2]<<2,
			rg=r-g,
			gb=g-b,
			br=b-r;

		counters[0]+=abs(r -prev[0]);
		counters[1]+=abs(g -prev[1]);
		counters[2]+=abs(b -prev[2]);
		counters[3]+=abs(rg-prev[3]);
		counters[4]+=abs(gb-prev[4]);
		counters[5]+=abs(br-prev[5]);
		prev[0]=r;
		prev[1]=g;
		prev[2]=b;
		prev[3]=rg;
		prev[4]=gb;
		prev[5]=br;
#ifdef ENABLE_EXTENDED_RCT
#define UPDATE(IDXA, A0, IDXB, B0, IDXC, C0)\
	do\
	{\
		int a0=A0, b0=B0, c0=C0;\
		counters[IDXA]+=abs(a0-prev[IDXA]);\
		counters[IDXB]+=abs(b0-prev[IDXB]);\
		counters[IDXC]+=abs(c0-prev[IDXC]);\
		prev[IDXA]=a0;\
		prev[IDXB]=b0;\
		prev[IDXC]=c0;\
	}while(0)

		UPDATE(
			OCH_CX31, rg+(gb>>2),/*r-(3*g+b)/4 = r-g-(b-g)/4*/
			OCH_C3X1, rg+(br>>2),/*g-(3*r+b)/4 = g-r-(b-r)/4*/
			OCH_C31X, br+(rg>>2) /*b-(3*r+g)/4 = b-r-(g-r)/4*/
		);
		UPDATE(
			OCH_CX13, br+(gb>>2),/*r-(g+3*b)/4 = r-b-(g-b)/4*/
			OCH_C1X3, gb+(br>>2),/*g-(r+3*b)/4 = g-b-(r-b)/4*/
			OCH_C13X, gb+(rg>>2) /*b-(r+3*g)/4 = b-g-(r-g)/4*/
		);
		UPDATE(
			OCH_CX22, (rg-br)>>1,/*r-(g+b)/2 = (r-g + r-b)/2*/
			OCH_C2X2, (gb-rg)>>1,/*g-(r+b)/2 = (g-r + g-b)/2*/
			OCH_C22X, (br-gb)>>1 /*b-(r+g)/2 = (b-r + b-g)/2*/
		);
#undef  UPDATE
#endif
	}
	rct=0;
	minerr=0;
	for(kt=0;kt<RCT_COUNT;++kt)
	{
		const uint8_t *group=rct_combinations[kt];
		int64_t currerr=
			+counters[group[0]]
			+counters[group[1]]
			+counters[group[2]]
		;
		if(!kt||minerr>currerr)
		{
			minerr=currerr;
			rct=kt;
		}
	}
	return rct;
}
static int32_t crct2_hist[OCH_COUNT][256];
int crct_analysis2(Image *src)
{
	enum
	{
		XPAD=8,
		NROWS=2,
		NCH=OCH_COUNT,
		NVAL=1
	};
	int psize;
	int16_t *pixels;
	double csizes[OCH_COUNT]={0}, minerr;
	int kt, ky, kx, kc, idx, rct;

	psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		CRASH("Alloc error");
		return 0;
	}
	memset(pixels, 0, psize);
	memset(crct2_hist, 0, sizeof(crct2_hist));
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int16_t *rows[4];

		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			int r=((int32_t*)HDR_DATA(src))[idx+0];
			int g=((int32_t*)HDR_DATA(src))[idx+1];
			int b=((int32_t*)HDR_DATA(src))[idx+2];
			int rg=r-g;
			int gb=g-b;
			int br=b-r;
			int yuv[OCH_COUNT], j;

			j=0;
			yuv[j++]=r;
			yuv[j++]=g;
			yuv[j++]=b;
			yuv[j++]=rg;
			yuv[j++]=gb;
			yuv[j++]=br;
			yuv[j++]=rg+(gb>>2);/*r-(3*g+b)/4 = r-g-(b-g)/4*/
			yuv[j++]=rg+(br>>2);/*g-(3*r+b)/4 = g-r-(b-r)/4*/
			yuv[j++]=br+(rg>>2);/*b-(3*r+g)/4 = b-r-(g-r)/4*/
			yuv[j++]=br+(gb>>2);/*r-(g+3*b)/4 = r-b-(g-b)/4*/
			yuv[j++]=gb+(br>>2);/*g-(r+3*b)/4 = g-b-(r-b)/4*/
			yuv[j++]=gb+(gb>>2);/*b-(r+3*g)/4 = b-g-(r-g)/4*/
			yuv[j++]=(rg-br)>>1;/*r-(g+b)/2 = (r-g + r-b)/2*/
			yuv[j++]=(rg-gb)>>1;/*g-(r+b)/2 = (g-r + g-b)/2*/
			yuv[j++]=(br-gb)>>1;/*b-(r+g)/2 = (b-r + b-g)/2*/
			for(kc=0;kc<OCH_COUNT;++kc)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				{
					int
						NW	=rows[1][0-1*NCH*NROWS*NVAL],
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL];
					int pred;

					MEDIAN3_CMOV(pred, N, W, N+W-NW);
#if 0
					pred=(N+W)>>1;
#endif
					++crct2_hist[kc][(uint8_t)(yuv[kc]-pred+128)];
					rows[0][0]=yuv[kc];
				}
			}
		}
	}
	_mm_free(pixels);
	for(kc=0;kc<OCH_COUNT;++kc)
	{
		int32_t *hist=crct2_hist[kc], sum=0;
		double e=0, norm;
		int ks;

		for(ks=0;ks<256;++ks)
			sum+=hist[ks];
		norm=1./sum;
		for(ks=0;ks<256;++ks)
		{
			int f=hist[ks];
			if(f)
				e-=f*log(f*norm);
		}
		csizes[kc]=e*(1./(8*M_LN2));
	}
	rct=0;
	minerr=0;
	for(kt=0;kt<RCT_COUNT;++kt)
	{
		const uint8_t *group=rct_combinations[kt];
		double currerr=
			+csizes[group[0]]
			+csizes[group[1]]
			+csizes[group[2]]
		;
		if(!kt||minerr>currerr)
		{
			minerr=currerr;
			rct=kt;
		}
	}
	return rct;
}
int crct2_analysis(Image *src)
{
	enum
	{
		STRIDE=7
	};
	int64_t counters[NPERMS*(1+RCTLEVELSU+RCTLEVELSV)]={0}, bestsum=0;
	RCTInfo rct={0};
	int rstr=0, rctdata=0;
	int ky, kx, kp, uc, vc;

	rstr=3*src->iw;
	for(ky=1;ky<=src->ih-STRIDE;ky+=STRIDE)
	{
		int *imptr=((int32_t*)HDR_DATA(src))+rstr*ky+4;
		for(kx=1;kx<=src->iw-STRIDE;kx+=STRIDE, imptr+=4*STRIDE)
		{
			int64_t *ctrptr=counters;
			int rgb[3], yuv[3];
			
			rgb[0]=(imptr[0]-imptr[-4+0]-imptr[-rstr+0]+imptr[-rstr-4+0])<<RCTBITS;
			rgb[1]=(imptr[1]-imptr[-4+1]-imptr[-rstr+1]+imptr[-rstr-4+1])<<RCTBITS;
			rgb[2]=(imptr[2]-imptr[-4+2]-imptr[-rstr+2]+imptr[-rstr-4+2])<<RCTBITS;
			for(kp=0;kp<NPERMS;++kp)
			{
				yuv[0]=rgb[perms[3*kp+0]];
				yuv[1]=rgb[perms[3*kp+1]];
				yuv[2]=rgb[perms[3*kp+2]];

				*ctrptr++ += abs(yuv[0]);
#if 1
				for(uc=0;uc<RCTLEVELSU;++uc)
					*ctrptr++ += abs(yuv[1]-((crct2_uc[uc]*yuv[0]+RCTROUND)>>RCTBITS));
				for(vc=0;vc<RCTLEVELSV;++vc)
					*ctrptr++ += abs(yuv[2]-((crct2_vc[2*vc+0]*yuv[0]+crct2_vc[2*vc+1]*yuv[1]+RCTROUND)>>RCTBITS));
#else
				for(uc0=0;uc0<RCTLEVELS;++uc0)
					*ctrptr++ += abs(yuv[1]-(uc0*yuv[0]>>RCTBITS));
				for(vc0=0;vc0<RCTLEVELS;++vc0)
				{
					for(vc1=0;vc0+vc1<RCTLEVELS;++vc1)
						*ctrptr++ += abs(yuv[2]-((vc0*yuv[0]+vc1*yuv[1])>>RCTBITS));
				}
#endif
			}
		}
	}
	bestsum=0;
#if 1
	for(kp=0;kp<NPERMS;++kp)
	{
		int64_t *currctrs=counters+(1+RCTLEVELSU+RCTLEVELSV)*kp;
		for(uc=0;uc<RCTLEVELSU;++uc)
		{
			for(vc=0;vc<RCTLEVELSV;++vc)
			{
				int64_t sum=currctrs[0]+currctrs[1+uc]+currctrs[1+RCTLEVELSU+vc];
				if(!bestsum||bestsum>sum)
				{
					bestsum=sum;
					rct.pidx=kp;
					rct.uc0=uc;
					rct.vc0=vc;
					rct.vc1=0;
				}
			}
		}
	}
	rctdata=rct.pidx;
	rctdata=RCTLEVELSU*rctdata+rct.uc0;
	rctdata=RCTLEVELSV*rctdata+rct.vc0;
	return rctdata|0x8000;
#else
	for(kp=0;kp<NPERMS;++kp)
	{
		int64_t *currctrs=counters+(1+RCTLEVELS+RCTLEVELS*(RCTLEVELS+1)/2)*kp;
		for(uc0=0;uc0<RCTLEVELS;++uc0)
		{
			for(vc0=0, idx=0;vc0<RCTLEVELS;++vc0)
			{
				for(vc1=0;vc0+vc1<RCTLEVELS;++vc1, ++idx)
				{
					int64_t sum=currctrs[0]+currctrs[1+uc0]+currctrs[1+RCTLEVELS+idx];
					if(uc0==1||vc0+vc1==1)
						continue;
					if(!bestsum||bestsum>sum)
					{
						bestsum=sum;
						rct.pidx=kp;
						rct.uc0=uc0;
						rct.vc0=vc0;
						rct.vc1=vc1;
					}
				}
			}
		}
	}
	rctdata=rct.pidx;
	rctdata=RCTLEVELS*rctdata+rct.uc0;
	rctdata=RCTLEVELS*rctdata+rct.vc0;
	rctdata=RCTLEVELS*rctdata+rct.vc1;
	return rctdata|0x8000;
#endif
}
const uint8_t* crct2_unpack(int rctdata, int *ret_uc0, int *ret_vc0, int *ret_vc1)
{
#if 1
	rctdata&=0x7FFF;
	*ret_vc0=crct2_vc[rctdata%RCTLEVELSV*2+0];
	*ret_vc1=crct2_vc[rctdata%RCTLEVELSV*2+1];
	rctdata/=RCTLEVELSV;
	*ret_uc0=crct2_uc[rctdata%RCTLEVELSU];
	rctdata/=RCTLEVELSU;
	if((uint32_t)rctdata>=(uint32_t)NPERMS)
		crash(__FILE__, __LINE__, 1, "Invalid RCT perm %d", rctdata);
	return perms+3*rctdata;
#else
	rctdata&=0x7FFF;
	*vc1=rctdata%RCTLEVELS;
	rctdata/=RCTLEVELS;
	*vc0=rctdata%RCTLEVELS;
	rctdata/=RCTLEVELS;
	*uc0=rctdata%RCTLEVELS;
	rctdata/=RCTLEVELS;
	if((uint32_t)rctdata>=(uint32_t)NPERMS)
		CRASH("Invalid RCT perm %d", rctdata);
	return perms+3*rctdata;
#endif
}

void pred_l1crct(Image *src, int fwd)
{
	enum
	{
		SHIFT=18,

		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=2
	};
	int amin[4], amax[4];
	int32_t weights[4][NPREDS]={0}, bias[4]={0};
	int psize;
	int16_t *pixels;
	int invdist;
	int uc0, vc0, vc1;
	const uint8_t *perm;
	int ky, kx, kc, idx;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i));
	invdist=((1<<16)+g_dist-1)/g_dist;
	uc0=vc0=vc1=0;
	if(fwd)
		src->rct=crct2_analysis(src);
	perm=crct2_unpack(src->rct, &uc0, &vc0, &vc1);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	FILLMEM((int*)weights, (1<<L1SH)/NPREDS, sizeof(weights), sizeof(int));
	bias[0]=1<<L1SH>>1;
	bias[1]=1<<L1SH>>1;
	bias[2]=1<<L1SH>>1;
	bias[3]=1<<L1SH>>1;
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int16_t *rows[4];

		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;/*sub 1 channel for pre-increment*/
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			int offset=0, yuv[3];

			yuv[0]=((int32_t*)HDR_DATA(src))[idx+perm[0]];
			yuv[1]=((int32_t*)HDR_DATA(src))[idx+perm[1]];
			yuv[2]=((int32_t*)HDR_DATA(src))[idx+perm[2]];
			for(kc=0;kc<4;++kc)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				{
					int
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
						W		=rows[0][0-1*NCH*NROWS*NVAL],
						aNW		=rows[1][1-1*NCH*NROWS*NVAL],
						aN		=rows[1][1+0*NCH*NROWS*NVAL],
						aNE		=rows[1][1+1*NCH*NROWS*NVAL],
						aNEEE		=rows[1][1+3*NCH*NROWS*NVAL],
						aW		=rows[0][1-1*NCH*NROWS*NVAL];
					int preds[NPREDS], j, k, *currw, p0, predc, vmax, vmin, curr, e;

#if 1
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
					(void)aNW;
					(void)aN;
					(void)aNE;
					(void)aNEEE;
					(void)aW;
#endif
#define PRED(W0, E) preds[j++]=E;
					j=0;
					PREDLIST;
#undef  PRED
					currw=weights[kc];
					p0=bias[kc];

					for(k=0;k<NPREDS;++k)
						p0+=currw[k]*preds[k];
					p0>>=L1SH;
					predc=p0;
					vmax=N; vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(predc, vmin, vmax);
					
					offset=0;
					if(kc==1)offset=uc0*yuv[0];
					if(kc==2)offset=vc0*yuv[0]+vc1*yuv[1];
					offset=(offset+RCTROUND)>>RCTBITS;
					predc+=offset;
					CLAMP2(predc, amin[kc], amax[kc]);

					curr=yuv[kc];
					if(g_dist>1)
					{
						if(fwd)
						{
							curr-=predc;
							curr=(curr*invdist>>16)-(curr>>31);
							((int32_t*)HDR_DATA(src))[idx+kc]=curr;

							curr=g_dist*curr+predc;
						}
						else
							curr=g_dist*((int32_t*)HDR_DATA(src))[idx+kc]+predc;
						CLAMP2(curr, amin[kc], amax[kc]);
						yuv[kc]=curr;
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
					rows[0][1]=curr;
					curr-=offset;
					rows[0][0]=curr;

					/*update*/
					e=(curr>p0)-(curr<p0);/*L1*/
#if 0
					e=curr-p0;/*L2 (faster rise, worse steady state)*/
#endif
					bias[kc]+=e;
					for(k=0;k<NPREDS;++k)
						currw[k]+=e*preds[k];
				}
			}
			if(!fwd)
			{
				((int32_t*)HDR_DATA(src))[idx+perm[0]]=yuv[0];
				((int32_t*)HDR_DATA(src))[idx+perm[1]]=yuv[1];
				((int32_t*)HDR_DATA(src))[idx+perm[2]]=yuv[2];
			}
		}
	}
	_mm_free(pixels);
}


void pred_grfilt(Image *src, int fwd)
{
	int psize;
	int16_t *pixels;
	int kx, ky, kc, idx;

	(void)fwd;
	psize=(src->iw+8*2)*(int)sizeof(int16_t[4*4*3]);/*4 padded rows * 4 channels max * {pixels, error}*/
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int16_t *rows[4];

		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0+4)%4)+8)*4*3;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1+4)%4)+8)*4*3;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2+4)%4)+8)*4*3;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3+4)%4)+8)*4*3;
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=3;
				rows[1]+=3;
				rows[2]+=3;
				rows[3]+=3;
				if(!src->depth[kc])
					continue;
				{
					int
						NNNWWW		=rows[3][0-3*4*3],
						NNNW		=rows[3][0-1*4*3],
						NNN		=rows[3][0+0*4*3],
						NNNE		=rows[3][0+1*4*3],
						NNNEE		=rows[3][0+2*4*3],
						NNNEEE		=rows[3][0+3*4*3],
						NNNEEEE		=rows[3][0+4*4*3],
						NNWWWW		=rows[2][0-4*4*3],
						NNWWW		=rows[2][0-3*4*3],
						NNWW		=rows[2][0-2*4*3],
						NNW		=rows[2][0-1*4*3],
						NN		=rows[2][0+0*4*3],
						NNE		=rows[2][0+1*4*3],
						NNEE		=rows[2][0+2*4*3],
						NNEEE		=rows[2][0+3*4*3],
						NNEEEE		=rows[2][0+4*4*3],
						NWWWW		=rows[1][0-4*4*3],
						NWWW		=rows[1][0-3*4*3],
						NWW		=rows[1][0-2*4*3],
						NW		=rows[1][0-1*4*3],
						N		=rows[1][0+0*4*3],
						NE		=rows[1][0+1*4*3],
						NEE		=rows[1][0+2*4*3],
						NEEE		=rows[1][0+3*4*3],
						NEEEE		=rows[1][0+4*4*3],
						NEEEEE		=rows[1][0+5*4*3],
						NEEEEEE		=rows[1][0+6*4*3],
						NEEEEEEE	=rows[1][0+7*4*3],
						NEEEEEEEE	=rows[1][0+8*4*3],
						WWWWWWWWW	=rows[0][0-9*4*3],
						WWWWWWWW	=rows[0][0-8*4*3],
						WWWWWWW		=rows[0][0-7*4*3],
						WWWWWW		=rows[0][0-6*4*3],
						WWWWW		=rows[0][0-5*4*3],
						WWWW		=rows[0][0-4*4*3],
						WWW		=rows[0][0-3*4*3],
						WW		=rows[0][0-2*4*3],
						W		=rows[0][0-1*4*3],
						eNE		=rows[1][1+1*4*3],
						eNEE		=rows[1][1+2*4*3],
						eNEEE		=rows[1][1+3*4*3],
						eW		=rows[0][1-1*4*3],
						e2NE		=rows[1][2+1*4*3],
						e2NEE		=rows[1][2+2*4*3],
						e2NEEE		=rows[1][2+3*4*3],
						e2W		=rows[0][2-1*4*3];
					int curr, e;

#if 1
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
					(void)eNE;
					(void)eNEE;
					(void)eNEEE;
					(void)eW;
					(void)e2NE;
					(void)e2NEE;
					(void)e2NEEE;
					(void)e2W;
#endif
					curr=((int32_t*)HDR_DATA(src))[idx];
					((int32_t*)HDR_DATA(src))[idx]=eW>>6;
					rows[0][0]=curr;

					e=abs(curr);
#if 0
					e=curr<<1^curr>>31;
#endif
					rows[0][1]=(eW+(eW<eNE?eW:eNE)+(e<<6)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
#if 0
					rows[0][1]=(2*eW+e+eNEEE)>>2;
#endif

					rows[0][2]=(e2W+(e2W<e2NE?e2W:e2NE)+eW+(e2NEE>e2NEEE?e2NEE:e2NEEE))>>2;
				}
			}
		}
	}
	free(pixels);
}

void pred_adaquant(Image *src, int fwd)
{
	enum
	{
		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=2,

		NPREDS2=4,
		SHIFT=18
	};
	int amin[4], amax[4];
	int32_t weights[NCH][NPREDS2]={0}, bias[NCH]={0}, estims[NPREDS2]={0};
	int psize;
	int16_t *pixels;
	int kx, ky, kc, idx;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	FILLMEM((int*)weights, (1<<L1SH)/NPREDS2, sizeof(weights), sizeof(int));
	bias[3]=bias[2]=bias[1]=bias[0]=1<<SHIFT>>1;
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int16_t *rows[4];

		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				{
					int16_t
						NNN	=rows[3][0+0*NCH*NROWS*NVAL],
						NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
						NNW	=rows[2][0-1*NCH*NROWS*NVAL],
						NN	=rows[2][0+0*NCH*NROWS*NVAL],
						NNE	=rows[2][0+1*NCH*NROWS*NVAL],
						NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
						NWW	=rows[1][0-2*NCH*NROWS*NVAL],
						NW	=rows[1][0-1*NCH*NROWS*NVAL],
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						NE	=rows[1][0+1*NCH*NROWS*NVAL],
						NEE	=rows[1][0+2*NCH*NROWS*NVAL],
						NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
						NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
						WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
						WWW	=rows[0][0-3*NCH*NROWS*NVAL],
						WW	=rows[0][0-2*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL],
						eN	=rows[1][1+0*NCH*NROWS*NVAL],
						eNE	=rows[1][1+1*NCH*NROWS*NVAL],
						eNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
						eW	=rows[0][1-1*NCH*NROWS*NVAL];
					int pred, p1, j, curr, error2, vmax, vmin, qden, e;

					j=0;

					(void)NNN;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
					(void)NWW;
					(void)NW;
					(void)N;
					(void)NE;
					(void)NEE;
					(void)NEEE;
					(void)NEEEE;
					(void)WWWW;
					(void)WWW;
					(void)WW;
					(void)W;
					(void)eN;
					(void)eNE;
					(void)eNEEE;
					(void)eW;
					j=0;
					estims[j++]=W;
					estims[j++]=N+W-NW;
					estims[j++]=2*N-NN;
					estims[j++]=NE;
					pred=(
						+bias[kc]
						+weights[kc][0]*estims[0]
						+weights[kc][1]*estims[1]
						+weights[kc][2]*estims[2]
						+weights[kc][3]*estims[3]
					)>>SHIFT;
					p1=pred;
					vmax=N; vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(pred, vmin, vmax);

					qden=eW+1;
					if(qden>(1<<g_dist))/*just to decrease penalty for easy regions*/
						qden=(1<<g_dist);
					curr=((int32_t*)HDR_DATA(src))[idx];
					if(fwd)
					{
						curr-=pred;
						curr/=qden;
						error2=curr;
						curr*=qden;
						curr+=pred;
						CLAMP2(curr, amin[kc], amax[kc]);
						((int32_t*)HDR_DATA(src))[idx]=error2;
					}
					else
					{
						error2=curr;
						curr*=qden;
						curr+=pred;
						CLAMP2(curr, amin[kc], amax[kc]);
						((int32_t*)HDR_DATA(src))[idx]=curr;
					}
					rows[0][0]=curr;
					rows[0][1]=(2*eW+((error2<<1^error2>>31)<<g_dist)+eNEEE)>>2;

					e=(curr>p1)-(curr<p1);/*L1*/
					weights[kc][0]+=e*estims[0];
					weights[kc][1]+=e*estims[1];
					weights[kc][2]+=e*estims[2];
					weights[kc][3]+=e*estims[3];
				}
			}
		}
	}
	_mm_free(pixels);
}

void pred_gray(Image *src, int fwd)
{
	int half[4], ky, kx, kc, idx, val;

	half[0]=1<<src->depth[0]>>1;
	half[1]=1<<src->depth[1]>>1;
	half[2]=1<<src->depth[2]>>1;
	half[3]=1<<src->depth[3]>>1;
	if(fwd)
	{
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx)
			{
				for(kc=0;kc<4;++kc, ++idx)
				{
					if(!src->depth[kc])
						continue;
					val=((int32_t*)HDR_DATA(src))[idx]+half[kc];
					val^=val>>1;
					((int32_t*)HDR_DATA(src))[idx]=val-half[kc];
				}
			}
		}
	}
	else
	{
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx)
			{
				for(kc=0;kc<4;++kc, ++idx)
				{
					if(!src->depth[kc])
						continue;
					val=((int32_t*)HDR_DATA(src))[idx]+half[kc];
					val^=val>>16;
					val^=val>> 8;
					val^=val>> 4;
					val^=val>> 2;
					val^=val>> 1;
					((int32_t*)HDR_DATA(src))[idx]=val-half[kc];
				}
			}
		}
	}
}

void pred_awav(Image *src, int fwd)
{
	int niter, rstr, kc, ky, start;

	niter=31^LZCNT32(((src->iw|src->ih)>>5)+1);
	rstr=4*src->iw;
	for(kc=0;kc<3;++kc)
	{
		int amin=-(1<<src->depth[kc]>>1), amax=(1<<src->depth[kc]>>1)-1;
		int it=fwd?0:niter-1;

		for(;;it+=fwd?1:-1)
		{
			int step=1<<it;
			int it2=fwd?0:3;
			if((uint32_t)it>(uint32_t)(niter-1))
				break;
			for(;;it2+=fwd?1:-1)
			{
				if((uint32_t)it2>3)
					break;
				for(ky=it2==2?step:0, start=(it2&1)^1;ky<=src->ih-step;ky+=step<<(it2>1), start^=it2<=1)
				{
					int kx=start?step:0;
					int *ptr=((int32_t*)HDR_DATA(src))+rstr*ky+4*kx+kc;
					for(;kx<=src->iw-step;kx+=2*step)
					{
						int pred=0, curr=0;
						int W=-step, E=+step, N=-step, S=+step;

						if(kx+W<0)W=+step;
						if(kx+E>=src->iw)E=-step;
						if(ky+N<0)N=+step;
						if(ky+S>=src->ih)S=-step;
						if(it2<=1)
						{
							pred+=ptr[rstr*N];
							pred+=ptr[rstr*S];
							pred+=ptr[4*W];
							pred+=ptr[4*E];
						}
						else
						{
							pred+=ptr[rstr*N+4*E];
							pred+=ptr[rstr*N+4*W];
							pred+=ptr[rstr*S+4*E];
							pred+=ptr[rstr*S+4*W];
						}
						pred>>=2;
						CLAMP2(pred, amin, amax);
						if(!fwd)
							pred=-pred;
						curr=*ptr;
						if(it2&1)
							curr+=pred;
						else
							curr-=pred;
						*ptr=curr;
						ptr+=4*2*step;
					}
				}
			}
		}
	}
}

static int squash(int x)
{
	enum
	{
		PROBBITS_SQUASHFUNC=12
	};
	static const int t[33]=/*2^5 table elements, table amplitude 2^12*/
	{
		   1,    2,    3,    6,   10,   16,   27,   45,   73,  120,  194,
		 310,  488,  747, 1101, 1546, 2047, 2549, 2994, 3348, 3607, 3785,
		3901, 3975, 4022, 4050, 4068, 4079, 4085, 4089, 4092, 4093, 4094,
	};
	int w;

	w=x&((1<<(PROBBITS_SQUASHFUNC-5))-1);
	x=(x>>(PROBBITS_SQUASHFUNC-5))+16;
	if(x>31)
		return (1<<PROBBITS_SQUASHFUNC)-1;
	if(x<0)
		return 1;
	x=(t[x]*((1<<(PROBBITS_SQUASHFUNC-5))-w)+t[x+1]*w+64)>>(12-5);
	return x;
}
void pred_extrap(Image *src)
{
#define PREDLIST2\
	PRED(N)\
	PRED(W)\
	PRED(2*N-NN)\
	PRED(2*W-WW)\
	PRED(N+W-NW)\
	PRED(W+NE-N)\
	PRED(N+NE-NNE)\
	PRED(W+NW-NWW)\

	enum
	{
		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=2,
#define PRED(E) +1
		NESTIM=PREDLIST2,
#undef  PRED
		USEBITS=12
	};
	int statssize=0;
	int32_t *stats=0;
	int psize=0;
	int32_t *pixels=0;
	uint32_t state=0x01234567;
	int ky, kx, kc, idx;

	#define ESTIM_CSIZE
#if defined ESTIM_CSIZE
	double csize=0;
#endif
	
	psize=(src->iw+2*XPAD)*(int)sizeof(int32_t[NROWS*NCH*NVAL]);
	pixels=(int32_t*)_mm_malloc(psize, sizeof(__m128i));
	statssize=sizeof(int32_t[NCH*NESTIM*256*256]);
	stats=(int32_t*)malloc(statssize);
	if(!pixels||!stats)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	memset(stats, 0, statssize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int32_t *rows[4];

		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(kc==3)
					continue;
				{
					int32_t
						NNN	=rows[3][0+0*NCH*NROWS*NVAL],
						NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
						NNW	=rows[2][0-1*NCH*NROWS*NVAL],
						NN	=rows[2][0+0*NCH*NROWS*NVAL],
						NNE	=rows[2][0+1*NCH*NROWS*NVAL],
						NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
						NWW	=rows[1][0-2*NCH*NROWS*NVAL],
						NW	=rows[1][0-1*NCH*NROWS*NVAL],
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						NE	=rows[1][0+1*NCH*NROWS*NVAL],
						NEE	=rows[1][0+2*NCH*NROWS*NVAL],
						NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
						NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
						WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
						WWW	=rows[0][0-3*NCH*NROWS*NVAL],
						WW	=rows[0][0-2*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL],
						cN	=rows[1][1+0*NCH*NROWS*NVAL],
						cW	=rows[0][1-1*NCH*NROWS*NVAL];
					int32_t *currstats[NESTIM];
					int j, k, tidx, kb, curr;

					(void)NNN;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
					(void)NWW;
					(void)NW;
					(void)N;
					(void)NE;
					(void)NEE;
					(void)NEEE;
					(void)NEEEE;
					(void)WWWW;
					(void)WWW;
					(void)WW;
					(void)W;
					(void)cN;
					(void)cW;
#define PRED(E) currstats[j]=stats+0x10000*j+256*(uint8_t)(E); ++j;
					j=0;
					PREDLIST2;
#undef  PRED
#if 1
					state^=state<<13;
					state^=state>>17;
					state^=state<<5;
					if(state>0x1FFFFFFF)/*synthesize*/
#else
					if((ky^kx)>>1&1)
#endif
					{
						for(kb=7, tidx=1;kb>=0;--kb)
						{
							int bit=0;
							int p0=0;
							for(k=0;k<NESTIM;++k)
								p0+=currstats[k][tidx];
							bit=p0<0;
							tidx=2*tidx+bit;
						}
						curr=(int8_t)tidx;
					}
					else/*learn*/
					{
						curr=((int32_t*)HDR_DATA(src))[idx];
						for(kb=7, tidx=1;kb>=0;--kb)
						{
							int bit, p0, proberror;

							bit=curr>>kb&1;
							p0=0;
							for(k=0;k<NESTIM;++k)
								p0+=currstats[k][tidx];
							p0>>=11;
							p0=squash(p0);
#if 0
							p0+=1<<USEBITS>>1;
							CLAMP2(p0, 1, (1<<USEBITS)-1);
#endif
#if defined ESTIM_CSIZE
							csize+=USEBITS-log((double)(bit?(1<<USEBITS)-p0:p0))*(1./M_LN2);
#endif
							proberror=((bit^1)<<USEBITS)-p0;
							for(k=0;k<NESTIM;++k)
								currstats[k][tidx]+=proberror;
#if 0
							{
								int truth=((bit^1)<<STOREBITS)-(1<<STOREBITS>>1)+(1<<7>>1);
								for(k=0;k<NESTIM;++k)
								{
									int32_t p=currstats[k][tidx];
									p+=(truth-p)>>7;
									currstats[k][tidx]=p;
								}
							}
#endif
							tidx=2*tidx+bit;
						}
					}
					rows[0][0]=((int32_t*)HDR_DATA(src))[idx]=curr;
				}
			}
		}
	}
	_mm_free(pixels);
	free(stats);
#if defined ESTIM_CSIZE
	messagebox(MBOX_OK, "Info", "%12.2lf bytes", csize/8);
#endif
}
void pred_l1dither(Image *src)
{
	enum
	{
		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=2,

		SHIFT=17
	};
	int weights[4][NPREDS]={0}, bias[4]={1<<SHIFT>>1, 1<<SHIFT>>1, 1<<SHIFT>>1, 1<<SHIFT>>1};
	int psize=0;
	int16_t *pixels=0;
	uint32_t state=0x01234567, synth=0;
	int ky, kx, kc, idx;

	psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	FILLMEM((int*)weights, (1<<SHIFT)/NPREDS, sizeof(weights), sizeof(int));
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int16_t *rows[4];

		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;/*sub 1 channel for pre-increment*/
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				++rows[0];
				++rows[1];
				++rows[2];
				++rows[3];
				if(!src->depth[kc])
					continue;
				{
					int32_t
						NNN	=rows[3][0+0*NCH*NROWS*NVAL],
						NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
						NNW	=rows[2][0-1*NCH*NROWS*NVAL],
						NN	=rows[2][0+0*NCH*NROWS*NVAL],
						NNE	=rows[2][0+1*NCH*NROWS*NVAL],
						NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
						NWW	=rows[1][0-2*NCH*NROWS*NVAL],
						NW	=rows[1][0-1*NCH*NROWS*NVAL],
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						NE	=rows[1][0+1*NCH*NROWS*NVAL],
						NEE	=rows[1][0+2*NCH*NROWS*NVAL],
						NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
						NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
						WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
						WWW	=rows[0][0-3*NCH*NROWS*NVAL],
						WW	=rows[0][0-2*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL],
						cN	=rows[1][1+0*NCH*NROWS*NVAL],
						cW	=rows[0][1-1*NCH*NROWS*NVAL];
					int preds[NPREDS], e, j, k, *currw, predc, p0, vmax, vmin, curr;

#if 1
					(void)NNN;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
					(void)NWW;
					(void)NW;
					(void)N;
					(void)NE;
					(void)NEE;
					(void)NEEE;
					(void)NEEEE;
					(void)WWWW;
					(void)WWW;
					(void)WW;
					(void)W;
					(void)cN;
					(void)cW;
#endif
#define PRED(W0, E) preds[j++]=E;
					j=0;
					PREDLIST;
#undef  PRED
					currw=weights[kc];
					predc=bias[kc];
					for(k=0;k<NPREDS;++k)
						predc+=currw[k]*preds[k];
					predc>>=SHIFT;
					p0=predc;
					vmax=N; vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(predc, vmin, vmax);

					curr=((int32_t*)HDR_DATA(src))[idx];
#if 1
					state^=state<<13;
					state^=state>>17;
					state^=state<<5;
					synth=state>0x7FFFFFFF;
					if(synth)
						((int32_t*)HDR_DATA(src))[idx]=curr=predc;
#else
					if(g_dist>1)
					{
						if(fwd)
						{
							curr-=predc;
						/*	curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));*/
							curr=(curr*invdist>>16)-(curr>>31);/*curr/=g_dist*/
							src->data[idx]=curr;
						}
						curr=g_dist*curr+predc;
						CLAMP2(curr, amin[kc], amax[kc]);
						if(!fwd)
							src->data[idx]=curr;
					}
					else
					{
						if(fwd)
						{
							int e=curr-predc;
							e<<=32-src->depth[kc];
							e>>=32-src->depth[kc];
							src->data[idx]=e;
						}
						else
						{
							curr+=predc;
							curr<<=32-src->depth[kc];
							curr>>=32-src->depth[kc];
							src->data[idx]=curr;
						}
					}
#endif
					rows[0][0]=curr;

					if(!synth)/*update*/
					{
						e=(curr>p0)-(curr<p0);/*L1*/
						bias[kc]+=e<<4;
						for(k=0;k<NPREDS;++k)
							currw[k]+=e*preds[k];
					}
				}
			}
		}
	}
	free(pixels);
}
