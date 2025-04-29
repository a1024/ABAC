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


//	#define SIMULATE_INTERLEAVE
//	#define PRINT_L1_BOUNDS


#define NPREDS 9
void pred_ols7(Image *src, int fwd)
{
#ifdef PRINT_L1_BOUNDS
	int cmin=0, cmax=0;
	int bmin=0, bmax=0;
#endif
#ifdef SIMULATE_INTERLEAVE
	int weights[4][4][4][NPREDS+1]={0};//coeffs: 18 bit	bias: 15 bit
#else
	int weights[4][NPREDS+1]={0};//coeffs: 18 bit	bias: 15 bit
#endif
	int bufsize=(src->iw+8*2)*(int)sizeof(int[4*4]);//4 padded rows * 4 channels max
	int *pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	//if(loud_transforms)
	//	DisableProcessWindowsGhosting();
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8)*4-1,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8)*4-1,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8)*4-1,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8)*4-1,
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
					NNNN	=rows[0][+0*4],
					NNNW	=rows[3][-1*4],
					NNN	=rows[3][+0*4],
					NNNE	=rows[3][+1*4],
					NN	=rows[2][+0*4],
					NNE	=rows[2][+1*4],
					NW	=rows[1][-1*4],
					N	=rows[1][+0*4],
					NE	=rows[1][+1*4],
					NEE	=rows[1][+2*4],
					NEEE	=rows[1][+3*4],
					NEEEE	=rows[1][+4*4],
					NEEEEE	=rows[1][+5*4],
					WWWWWW	=rows[0][-6*4],
					WWWWW	=rows[0][-5*4],
					WWWW	=rows[0][-4*4],
					WWW	=rows[0][-3*4],
					WW	=rows[0][-2*4],
					W	=rows[0][-1*4];
				int preds[]=
				{
					N,
					W,
					3*(N-NN)+NNN,
					3*(W-WW)+WWW,
					W+NE-N,
					N+W-NW,
					(WWWWW+WW-W+NNN+N+NEEEEE)>>2,
					N+NE-NNE,
					NEEE,

				//	N,
				//	W,
				//	3*(N-NN)+NNN,
				//	3*(W-WW)+WWW,
				//	W+NE-N,
				//	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4,
				//	N+W-NW,
				//	N+NE-NNE,
				};
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");
#ifdef SIMULATE_INTERLEAVE
				int *currw=weights[kc][ky*4/src->ih][kx*4/src->iw];
#else
				int *currw=weights[kc];
#endif
				int pred=currw[NPREDS];
#ifdef PRINT_L1_BOUNDS
				if(bmin>currw[NPREDS])bmin=currw[NPREDS];
				if(bmax<currw[NPREDS])bmax=currw[NPREDS];
#endif
				for(int k=0;k<NPREDS;++k)
				{
					pred+=currw[k]*preds[k];
#ifdef PRINT_L1_BOUNDS
					if(cmin>currw[k])cmin=currw[k];
					if(cmax<currw[k])cmax=currw[k];
#endif
				}
				pred+=1<<(src->depth[kc]+11)>>1;
				pred>>=src->depth[kc]+11;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NW)vmin=NW;
				if(vmax<NW)vmax=NW;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);

				int curr=src->data[idx];
				int pred2=(int)pred;
				if(!keyboard[KEY_SHIFT])
				{
					pred2^=-fwd;
					pred2+=fwd;
					pred2+=curr;
					pred2<<=32-src->depth[kc];
					pred2>>=32-src->depth[kc];
				}
				src->data[idx]=(int)pred2;
				if(!fwd)
					curr=pred2;
				rows[0][0]=curr;

				//update
				int e=curr-(int)pred;
				//	... -8 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7
				//	... -1 -1 -1 -1 -2 -2 -2 -2  0 +2 +2 +2 +2 +1 +1 +1 ...
				e=(((e>0)-(e<0))<<1)-((e>4)-(e<-4));
				currw[NPREDS]+=e;
				for(int k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];
			}
		}
	}
#ifdef PRINT_L1_BOUNDS
	if(loud_transforms)
		messagebox(MBOX_OK, "Info",
			"Coeffs %8d ~ %8d\n"
			"Bias   %8d ~ %8d\n"
			, cmin, cmax
			, bmin, bmax
		);
#endif
}
