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


#if 0
#define L1SH 18
#define PREDLIST\
	PRED(N)\
	PRED(W)\
	PRED(NNN)\
	PRED(WWW)\
	PRED(NEEE)\
	PRED(3*(N-NN)+NNN)\
	PRED(3*(W-WW)+WWW)\
	PRED(W+NE-N)\
	PRED(N+W-NW)\
	PRED(N+NE-NNE)\
	PRED((WWWW+WWW+NNN+NNEE+NEEE+NEEEE-(N+W))>>2)\

#endif
#if 1
#define L1SH 26
#define NBITS 4
#define PREDLIST\
	PRED(4*(WWWW+WWW+NNN)+3*NEEEE+NW)\
	PRED(24*W-8*WW-4*NW+3*N+NE)\
	PRED(6*W+8*NE+2*NEE)\
	PRED(24*N-8*NN+4*(W-NW+NE-NNE))\
	PRED(8*(N+W)+2*(NE-NW))\
	PRED(20*N-4*NN)\
	PRED(20*W-4*WW)\

#endif
#if 0
#define L1SH 16
#define NBITS 0
#define PREDLIST\
	PRED(N)\
	PRED(W)\
	PRED(N+W-NW)\
	PRED(NE)\
	PRED(NEE)\
	PRED(NN)\
	PRED(WW)\

#endif
enum
{
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED
};


void pred_ols9(Image *src, int fwd)
{
	int weights[4][NPREDS]={0};
	//int weights[4][2][2][NPREDS]={0};
	int bufsize=(src->iw+8*2)*(int)sizeof(short[4*4*(NPREDS+1)]);//4 padded rows * 4 channels max * {pixel, prederrors...}
	short *pixels=(short*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	FILLMEM((int*)weights, (1<<L1SH>>NBITS)/NPREDS, sizeof(weights), sizeof(int));
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((src->iw+16LL)*((ky-0LL+4)%4)+8)*4-1)*(NPREDS+1),
			pixels+(((src->iw+16LL)*((ky-1LL+4)%4)+8)*4-1)*(NPREDS+1),
			pixels+(((src->iw+16LL)*((ky-2LL+4)%4)+8)*4-1)*(NPREDS+1),
			pixels+(((src->iw+16LL)*((ky-3LL+4)%4)+8)*4-1)*(NPREDS+1),
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NPREDS+1;
				rows[1]+=NPREDS+1;
				rows[2]+=NPREDS+1;
				rows[3]+=NPREDS+1;
				if(!src->depth[kc])
					continue;
				int
					NNNWWW		=rows[3][-3*4*(NPREDS+1)],
					NNNW		=rows[3][-1*4*(NPREDS+1)],
					NNN		=rows[3][+0*4*(NPREDS+1)],
					NNNE		=rows[3][+1*4*(NPREDS+1)],
					NNNEE		=rows[3][+2*4*(NPREDS+1)],
					NNNEEE		=rows[3][+3*4*(NPREDS+1)],
					NNNEEEE		=rows[3][+4*4*(NPREDS+1)],
					NNWWWW		=rows[2][-4*4*(NPREDS+1)],
					NNWWW		=rows[2][-3*4*(NPREDS+1)],
					NNWW		=rows[2][-2*4*(NPREDS+1)],
					NNW		=rows[2][-1*4*(NPREDS+1)],
					NN		=rows[2][+0*4*(NPREDS+1)],
					NNE		=rows[2][+1*4*(NPREDS+1)],
					NNEE		=rows[2][+2*4*(NPREDS+1)],
					NNEEE		=rows[2][+3*4*(NPREDS+1)],
					NNEEEE		=rows[2][+4*4*(NPREDS+1)],
					NWWWW		=rows[1][-4*4*(NPREDS+1)],
					NWWW		=rows[1][-3*4*(NPREDS+1)],
					NWW		=rows[1][-2*4*(NPREDS+1)],
					NW		=rows[1][-1*4*(NPREDS+1)],
					N		=rows[1][+0*4*(NPREDS+1)],
					NE		=rows[1][+1*4*(NPREDS+1)],
					NEE		=rows[1][+2*4*(NPREDS+1)],
					NEEE		=rows[1][+3*4*(NPREDS+1)],
					NEEEE		=rows[1][+4*4*(NPREDS+1)],
					NEEEEE		=rows[1][+5*4*(NPREDS+1)],
					NEEEEEE		=rows[1][+6*4*(NPREDS+1)],
					NEEEEEEE	=rows[1][+7*4*(NPREDS+1)],
					NEEEEEEEE	=rows[1][+8*4*(NPREDS+1)],
					WWWWWWWWW	=rows[0][-9*4*(NPREDS+1)],
					WWWWWWWW	=rows[0][-8*4*(NPREDS+1)],
					WWWWWWW		=rows[0][-7*4*(NPREDS+1)],
					WWWWWW		=rows[0][-6*4*(NPREDS+1)],
					WWWWW		=rows[0][-5*4*(NPREDS+1)],
					WWWW		=rows[0][-4*4*(NPREDS+1)],
					WWW		=rows[0][-3*4*(NPREDS+1)],
					WW		=rows[0][-2*4*(NPREDS+1)],
					W		=rows[0][-1*4*(NPREDS+1)];
				int preds[]=
				{
#define PRED(EXPR) EXPR,
					PREDLIST
#undef  PRED
				};

				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");
#if 0
				int wp=0, wsum=0;
				for(int k=0;k<NPREDS;++k)
				{
					int e=(
						+rows[2][k+1+0*(NPREDS+1)]	//NN
						+rows[1][k+1-1*(NPREDS+1)]	//NW
						+rows[1][k+1+0*(NPREDS+1)]	//N
						+rows[1][k+1+1*(NPREDS+1)]	//NE
						+rows[1][k+1+2*(NPREDS+1)]	//NEE
						+rows[0][k+1-2*(NPREDS+1)]	//WW
						+rows[0][k+1-1*(NPREDS+1)]	//W
					);
					e=0x10000/(e+1);
					wp+=e*preds[k];
					wsum+=e;
				}
				wsum<<=NBITS;
				wp+=wsum>>1;
				wp/=wsum+1;
#endif
#if 1
				int *currw=weights[kc];
				//int *currw=weights[kc][ky&1][kx&1];
				long long p0=1<<L1SH>>1;
				for(int k=0;k<NPREDS;++k)
					p0+=(long long)currw[k]*preds[k];
				p0>>=L1SH;
#endif
			//	int pred=((int)p0+wp)>>1;
				int pred=(int)p0;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);

				//if(ky==src->ih/2&&kx>10&&!kc)//
				//	printf("");

				int curr=src->data[idx];
				if(fwd)
				{
					int res=curr-pred;
					res<<=32-src->depth[kc];
					res>>=32-src->depth[kc];
					src->data[idx]=res;
				}
				else
				{
					curr+=pred;
					curr<<=32-src->depth[kc];
					curr>>=32-src->depth[kc];
					src->data[idx]=curr;
				}
				rows[0][0]=curr;
#if 0
				for(int k=0;k<NPREDS;++k)
					rows[0][k+1+0*(NPREDS+1)]=abs(curr-((preds[k]+(1<<NBITS>>1))>>NBITS));
#endif
#if 1
				int e=(curr>p0)-(curr<p0);
				for(int k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];
#endif
			}
		}
	}
	free(pixels);
}
