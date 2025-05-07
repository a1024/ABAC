#include"ebench.h"
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
static const char file[]=__FILE__;


#define L1SH 19
#define NPREDS 10
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
	PRED( 40000, NEEE)

void pred_ols7(Image *src, int fwd)
{
	const int wsize=sizeof(long long[4][NPREDS]);
	long long *weights=(long long*)malloc(wsize);
	int bufsize=(src->iw+8*2)*(int)sizeof(short[4*4*1]);//4 padded rows * 4 channels max * {pixels}
	short *pixels=(short*)malloc(bufsize);
	if(!pixels||!weights)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	static const int w0[]=
	{
#define PRED(W0, EXPR) W0,
		PREDLIST
#undef  PRED
	};
	for(int kc=0;kc<4;++kc)
	{
		for(int kp=0;kp<NPREDS;++kp)
			weights[kc*NPREDS+kp]=w0[kp];
	}
	//memset(weights, 0, wsize);
	//FILLMEM(weights, (1<<L1SH)/NPREDS, wsize, sizeof(*weights));
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((src->iw+16LL)*((ky-0LL+4)%4)+8)*4-1)*1,
			pixels+(((src->iw+16LL)*((ky-1LL+4)%4)+8)*4-1)*1,
			pixels+(((src->iw+16LL)*((ky-2LL+4)%4)+8)*4-1)*1,
			pixels+(((src->iw+16LL)*((ky-3LL+4)%4)+8)*4-1)*1,
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				++rows[0];
				++rows[1];
				++rows[2];
				++rows[3];
				if(!src->depth[kc])
					continue;
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
				int preds[]=
				{
#define PRED(W0, EXPR) EXPR,
					PREDLIST
#undef  PRED
				};
				long long *currw=weights+NPREDS*kc;
				long long pred=0;
				for(int k=0;k<NPREDS;++k)
					pred+=currw[k]*preds[k];
				long long p0=pred;
				pred+=1LL<<L1SH>>1;
				pred>>=L1SH;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);

				int curr=src->data[idx];
				int val=(int)pred;
				if(!keyboard[KEY_SHIFT])
				{
					val^=-fwd;
					val+=fwd;
					val+=curr;
					val<<=32-src->depth[kc];
					val>>=32-src->depth[kc];
				}
				src->data[idx]=val;
				if(!fwd)
					curr=val;
				rows[0][0]=curr;

				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");

				//update
				int e=(curr>pred)-(curr<pred);
				for(int k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];
			}
		}
	}
	free(pixels);
	free(weights);
}
