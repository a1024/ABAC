#include"ebench.h"
#include<stdint.h>
#include<stdio.h>
#include<immintrin.h>
static const char file[]=__FILE__;

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
#define PRED(...) +1
	NPREDS3	=PREDLIST3,
	NPREDS2	=PREDLIST2,
	NPREDS1	=PREDLIST1,
	NPREDS0	=PREDLIST0,
#undef  PRED
	XPAD	=8,
	NROWS	=4,
	NCH	=4,
	NVAL	=1,
};
#define NMAX 16
typedef struct _LPCInfo
{
	int estim[NMAX];
	double autocorr[NMAX+1], coeffs[NMAX];
} LPCInfo;
static void lpc_update(LPCInfo *lpc, int count, double lr)
{
	//HALAC 1.9.0
#if 1
	double coeffs[NMAX+1]={1}, temp[NMAX+1];
	double error=lpc->autocorr[0];
	if(!error)
		error=1;
	for(int i=0;i<count;++i)
	{
		double r=-coeffs[0]*lpc->autocorr[i+1];
		for(int j=1;j<=i;++j)
			r-=coeffs[j]*lpc->autocorr[i+1-j];
		r/=error;
		for(int j=0;j<=i+1;++j)
			temp[j]=coeffs[j]+r*coeffs[i+1-j];
		for(int j=0;j<=i+1;++j)
			coeffs[j]=temp[j];
		//memcpy(coeffs, temp, sizeof(double)*(i+1));

		error-=error*r*r;
	}
	for(int k=0;k<count;++k)
		lpc->coeffs[k]+=(coeffs[k+1]-lpc->coeffs[k])*lr;
#endif

	//vorbis
#if 0
	double error=lpc->autocorr[0]*(1+1e-10);
	double epsilon=lpc->autocorr[0]*1e-9+1e-10;
	for(int i=0;i<count;++i)
	{
		double r=-lpc->autocorr[i+1];//reflection
		if(error<epsilon)
			break;
		for(j=0;j<i;++j)
			r-=lpc[j]*autocorr[i-j];
		r/=error;
		lpc[i]=r;
		for(j=0;j<i/2;j++)
		{
			double tmp=lpc[j];

			lpc[j]+=r*lpc[i-1-j];
			lpc[i-1-j]+=r*tmp;
		}
		if(i&1)
			lpc[j]+=lpc[j]*r;

		error-=error*r*r;
	}
	for(int k=0;k<NPREDS;++k)
		coeffs[k]+=(lpc[k]-coeffs[k])*1;
#endif
}
void pred_lpc(Image *src, int fwd)
{
	double coeffs[4][4]={0};
	int psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	int16_t *pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i*));
	const int lpcsize=(int)sizeof(LPCInfo[4][4]);
	LPCInfo *lpc=(LPCInfo*)malloc(lpcsize);
	if(!pixels||!lpc)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	memset(lpc, 0, lpcsize);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int j=0, pred=0, curr=0;
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS-NROWS+(ky+NROWS-0LL)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky+NROWS-1LL)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky+NROWS-2LL)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky+NROWS-3LL)%NROWS)*NVAL,
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
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
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				//if(vmin>NE)vmin=NE;
				//if(vmax<NE)vmax=NE;
				//if(vmin>NEEE)vmin=NEEE;
				//if(vmax<NEEE)vmax=NEEE;

				//predict
				double estim[4]={0};
				LPCInfo *currlpc=lpc+4*kc;
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
				double fpred=
					+coeffs[kc][0]*estim[0]
					+coeffs[kc][1]*estim[1]
					+coeffs[kc][2]*estim[2]
					+coeffs[kc][3]*estim[3]
				;
				double p0=fpred;
				//CLAMP2(fpred, vmin, vmax);	//TUNE
				//pred=(int)CVTFP64_I64(estim[0]);
				pred=N+W-NW;
				CLAMP2(pred, vmin, vmax);

				curr=src->data[idx];
				if(fwd)
				{
					int error=curr-pred;
					error<<=32-src->depth[kc];
					error>>=32-src->depth[kc];
					src->data[idx]=error;
				}
				else
				{
					curr+=pred;
					curr<<=32-src->depth[kc];
					curr>>=32-src->depth[kc];
					src->data[idx]=curr;
				}
				rows[0][0]=curr;
				
				//if(ky>=src->ih/2&&kx>=src->iw/2&&!kc)//
				//	printf("");//

				//update
				double e=(double)((curr>p0)-(curr<p0))/0x40000;
				coeffs[kc][0]+=e*estim[0];
				coeffs[kc][1]+=e*estim[1];
				coeffs[kc][2]+=e*estim[2];
				coeffs[kc][3]+=e*estim[3];

				currlpc[0].estim[NPREDS0-1]=curr;//special
#define PRED(...) currlpc[0].autocorr[j]+=currlpc[0].estim[j]*curr-currlpc[0].autocorr[j]*0.00; ++j;
				j=0;
				PREDLIST0
#undef  PRED
#if 0
#define PRED(...) currlpc[1].autocorr[j]+=currlpc[1].estim[j]*curr-currlpc[1].autocorr[j]*0.01; ++j;
				j=0;
				PREDLIST1
#undef  PRED
#define PRED(...) currlpc[2].autocorr[j]+=currlpc[2].estim[j]*curr-currlpc[2].autocorr[j]*0.01; ++j;
				j=0;
				PREDLIST2
#undef  PRED
#define PRED(...) currlpc[3].autocorr[j]+=currlpc[3].estim[j]*curr-currlpc[3].autocorr[j]*0.01; ++j;
				j=0;
				PREDLIST3
#undef  PRED
#endif
				lpc_update(currlpc+0, NPREDS0, 1);
				//lpc_update(currlpc+1, NPREDS1, 0.001);
				//lpc_update(currlpc+2, NPREDS2, 0.001);
				//lpc_update(currlpc+3, NPREDS3, 0.001);
			}
		}
	}
	_mm_free(pixels);
}
