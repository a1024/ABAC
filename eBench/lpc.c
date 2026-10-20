#include"ebench.h"
#include<stdint.h>
#include<stdio.h>
#include<immintrin.h>

#if 1
#define PREDLIST0\
	PRED(WWWWW)\
	PRED(WWWW)\
	PRED(WWW)\
	PRED(WW)\
	PRED(W)\
	PRED(0)\

#define PREDLIST1\
	PRED(NWWW)\
	PRED(NWW)\
	PRED(NW)\
	PRED(N)\
	PRED(NE)\
	PRED(NEE)\
	PRED(NEEE)\

#define PREDLIST2\
	PRED(NNWW)\
	PRED(NNW)\
	PRED(NN)\
	PRED(NNE)\
	PRED(NNEE)\

#define PREDLIST3\
	PRED(NNNW)\
	PRED(NNN)\
	PRED(NNNE)\

#endif
#if 0
#define PREDLIST\
	PRED(N+W-NW)\
	PRED(N)\
	PRED(W)\
	PRED(NE)\
	PRED(2*N-NN)\
	PRED(2*W-WW)\
	PRED(2*(NE-NNE)+NN)\
	PRED(W+NEE-NE)\

#endif
enum
{
#define PRED(E) +1
	NPREDS3	=PREDLIST3,
	NPREDS2	=PREDLIST2,
	NPREDS1	=PREDLIST1,
	NPREDS0	=PREDLIST0,
#undef  PRED
	XPAD	=8,
	NROWS	=4,
	NCH	=4,
	NVAL	=1
};
#define NMAX 16
typedef struct _LPCInfo
{
	int estim[NMAX];
	double autocorr[NMAX+1], coeffs[NMAX];
} LPCInfo;
static void lpc_update(LPCInfo *lpc, int count, double lr)
{
	/*HALAC 1.9.0*/
#if 1
	double coeffs[NMAX+1]={0}, temp[NMAX+1];
	double e=lpc->autocorr[0];
	int i, j, k;

	coeffs[0]=1;
	if(fabs(e)<1e-9)
		e=1;
	for(i=0;i<count;++i)
	{
		double r=-coeffs[0]*lpc->autocorr[i+1];
		for(j=1;j<=i;++j)
			r-=coeffs[j]*lpc->autocorr[i+1-j];
		r/=e;
		for(j=0;j<=i+1;++j)
			temp[j]=coeffs[j]+r*coeffs[i+1-j];
		for(j=0;j<=i+1;++j)
			coeffs[j]=temp[j];
#if 0
		memcpy(coeffs, temp, sizeof(double)*(i+1));
#endif
		e-=e*r*r;
	}
	for(k=0;k<count;++k)
		lpc->coeffs[k]+=(coeffs[k+1]-lpc->coeffs[k])*lr;
#endif

	/*vorbis*/
#if 0
	double e=lpc->autocorr[0]*(1+1e-10);
	double epsilon=lpc->autocorr[0]*1e-9+1e-10;
	int i, j, k;

	for(i=0;i<count;++i)
	{
		double r=-lpc->autocorr[i+1];/*reflection*/
		if(e<epsilon)
			break;
		for(j=0;j<i;++j)
			r-=lpc[j]*autocorr[i-j];
		r/=e;
		lpc[i]=r;
		for(j=0;j<i/2;j++)
		{
			double tmp=lpc[j];

			lpc[j]+=r*lpc[i-1-j];
			lpc[i-1-j]+=r*tmp;
		}
		if(i&1)
			lpc[j]+=lpc[j]*r;

		e-=e*r*r;
	}
	for(k=0;k<NPREDS;++k)
		coeffs[k]+=(lpc[k]-coeffs[k])*1;
#endif
}
void pred_lpc(Image *src, int fwd)
{
	double coeffs[4][4]={0};
	int psize;
	int16_t *pixels;
	int lpcsize;
	LPCInfo *lpc;
	int kc, kx, ky, idx;

	psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i*));
	lpcsize=(int)sizeof(LPCInfo[4][4]);
	lpc=(LPCInfo*)malloc(lpcsize);
	if(!pixels||!lpc)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	memset(lpc, 0, lpcsize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int j=0, pred=0, curr=0;
		int16_t *rows[4];
		
		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky+(intptr_t)NROWS-0)%NROWS)*NVAL;
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky+(intptr_t)NROWS-1)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky+(intptr_t)NROWS-2)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky+(intptr_t)NROWS-3)%NROWS)*NVAL;
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
						NNNW	=rows[3][0-1*NROWS*NCH*NVAL],
						NNN	=rows[3][0+0*NROWS*NCH*NVAL],
						NNNE	=rows[3][0+1*NROWS*NCH*NVAL],
						NNWW	=rows[2][0-2*NROWS*NCH*NVAL],
						NNW	=rows[2][0-1*NROWS*NCH*NVAL],
						NN	=rows[2][0+0*NROWS*NCH*NVAL],
						NNE	=rows[2][0+1*NROWS*NCH*NVAL],
						NNEE	=rows[2][0+2*NROWS*NCH*NVAL],
						NWWW	=rows[1][0-3*NROWS*NCH*NVAL],
						NWW	=rows[1][0-2*NROWS*NCH*NVAL],
						NW	=rows[1][0-1*NROWS*NCH*NVAL],
						N	=rows[1][0+0*NROWS*NCH*NVAL],
						NE	=rows[1][0+1*NROWS*NCH*NVAL],
						NEE	=rows[1][0+2*NROWS*NCH*NVAL],
						NEEE	=rows[1][0+3*NROWS*NCH*NVAL],
						NEEEE	=rows[1][0+4*NROWS*NCH*NVAL],
						WWWWW	=rows[0][0-5*NROWS*NCH*NVAL],
						WWWW	=rows[0][0-4*NROWS*NCH*NVAL],
						WWW	=rows[0][0-3*NROWS*NCH*NVAL],
						WW	=rows[0][0-2*NROWS*NCH*NVAL],
						W	=rows[0][0-1*NROWS*NCH*NVAL];
					double estim[4]={0}, fpred, p0;
					LPCInfo *currlpc;
					int vmax, vmin;

#if 1
					(void)NNNW;
					(void)NNN;
					(void)NNNE;
					(void)NNWW;
					(void)NNW;
					(void)NN;
					(void)NNE;
					(void)NNEE;
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
#endif
					vmax=N; vmin=W;
					if(N<W)vmin=N, vmax=W;
#if 0
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
#endif

					/*predict*/
					currlpc=lpc+4*kc;
#define PRED(E) estim[0]+=currlpc[0].coeffs[j]*(currlpc[0].estim[j]=E); ++j;
					j=0;
					PREDLIST0
#undef  PRED
#if 0
#define PRED(E) estim[1]+=currlpc[1].coeffs[j]*(currlpc[1].estim[j]=E); ++j;
					j=0;
					PREDLIST1
#undef  PRED
#define PRED(E) estim[2]+=currlpc[2].coeffs[j]*(currlpc[2].estim[j]=E); ++j;
					j=0;
					PREDLIST2
#undef  PRED
#define PRED(E) estim[3]+=currlpc[3].coeffs[j]*(currlpc[3].estim[j]=E); ++j;
					j=0;
					PREDLIST3
#undef  PRED
#endif
					fpred=
						+coeffs[kc][0]*estim[0]
						+coeffs[kc][1]*estim[1]
						+coeffs[kc][2]*estim[2]
						+coeffs[kc][3]*estim[3]
					;
					p0=fpred;
#if 0
					CLAMP2(fpred, vmin, vmax);	//TUNE
					pred=(int)ROUND64(estim[0]);
#endif
					pred=N+W-NW;
					CLAMP2(pred, vmin, vmax);

					curr=((int32_t*)HDR_DATA(src))[idx];
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
					rows[0][0]=curr;

					/*update*/
					{
						double e=(double)((curr>p0)-(curr<p0))/0x40000;
						coeffs[kc][0]+=e*estim[0];
						coeffs[kc][1]+=e*estim[1];
						coeffs[kc][2]+=e*estim[2];
						coeffs[kc][3]+=e*estim[3];
					}

					currlpc[0].estim[NPREDS0-1]=curr;/*special*/
#define PRED(E) currlpc[0].autocorr[j]+=currlpc[0].estim[j]*curr-currlpc[0].autocorr[j]*0.00; ++j;
					j=0;
					PREDLIST0
#undef  PRED
#if 0
#define PRED(E) currlpc[1].autocorr[j]+=currlpc[1].estim[j]*curr-currlpc[1].autocorr[j]*0.01; ++j;
					j=0;
					PREDLIST1
#undef  PRED
#define PRED(E) currlpc[2].autocorr[j]+=currlpc[2].estim[j]*curr-currlpc[2].autocorr[j]*0.01; ++j;
					j=0;
					PREDLIST2
#undef  PRED
#define PRED(E) currlpc[3].autocorr[j]+=currlpc[3].estim[j]*curr-currlpc[3].autocorr[j]*0.01; ++j;
					j=0;
					PREDLIST3
#undef  PRED
#endif
					lpc_update(currlpc+0, NPREDS0, 1);
#if 0
					lpc_update(currlpc+1, NPREDS1, 0.001);
					lpc_update(currlpc+2, NPREDS2, 0.001);
					lpc_update(currlpc+3, NPREDS3, 0.001);
#endif
				}
			}
		}
	}
	_mm_free(pixels);
}
