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

#define NPREDS 8
void pred_ols7(Image *src, int fwd)
{
	int weights[4][NPREDS+2]={0};
	weights[0][NPREDS+1]=1<<src->depth[0];
	weights[1][NPREDS+1]=1<<src->depth[1];
	weights[2][NPREDS+1]=1<<src->depth[2];
	weights[3][NPREDS+1]=1<<src->depth[3];
	int bufsize=(src->iw+8*2)*(int)sizeof(int[4*4*(1+NPREDS+1)]);//4 padded rows * 4 channels max * {pixel, 8 weights, bias}
	int *pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	if(loud_transforms)
		DisableProcessWindowsGhosting();
	memset(pixels, 0, bufsize);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8)*4*(1+NPREDS+1)-1*(1+NPREDS+1),
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8)*4*(1+NPREDS+1)-1*(1+NPREDS+1),
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8)*4*(1+NPREDS+1)-1*(1+NPREDS+1),
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8)*4*(1+NPREDS+1)-1*(1+NPREDS+1),
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=1+NPREDS+1;
				rows[1]+=1+NPREDS+1;
				rows[2]+=1+NPREDS+1;
				rows[3]+=1+NPREDS+1;
				if(!src->depth[kc])
					continue;
				int
					NNNN	=rows[0][+0*4*(1+NPREDS+1)],
					NNN	=rows[3][+0*4*(1+NPREDS+1)],
					NNNE	=rows[3][+1*4*(1+NPREDS+1)],
					NN	=rows[2][+0*4*(1+NPREDS+1)],
					NNE	=rows[2][+1*4*(1+NPREDS+1)],
					NW	=rows[1][-1*4*(1+NPREDS+1)],
					N	=rows[1][+0*4*(1+NPREDS+1)],
					NE	=rows[1][+1*4*(1+NPREDS+1)],
					NEE	=rows[1][+2*4*(1+NPREDS+1)],
					NEEE	=rows[1][+3*4*(1+NPREDS+1)],
					NEEEE	=rows[1][+4*4*(1+NPREDS+1)],
					WWWWW	=rows[0][-5*4*(1+NPREDS+1)],
					WWWW	=rows[0][-4*4*(1+NPREDS+1)],
					WWW	=rows[0][-3*4*(1+NPREDS+1)],
					WW	=rows[0][-2*4*(1+NPREDS+1)],
					W	=rows[0][-1*4*(1+NPREDS+1)];
				int preds[]=
				{
#if 1
					N,
					W,
					3*(N-NN)+NNN,
					3*(W-WW)+WWW,
				//	3*NN-2*NNN,//X
				//	3*WW-2*WWW,//X
					W+NE-N,
					N+W-NW,
					(WWWW+NNN+NEE+NEEEE)>>2,
				//	(WWWW+WWW+NNN+NEEE+NEEEE-NW)>>2,//X
				//	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4,//X
					N+NE-NNE,
#endif
#if 0
					N,
					W,
					NW,
					NE,
					2*N-NN,
					2*W-WW,
					3*(N-NN)+NNN,
					3*(W-WW)+WWW,
					NEE,
					NEEE,
					(WWWW+NNN+NEE+NEEEE)>>2,
					(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4,
#endif
				};
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");
				int *currw=weights[kc];
			//	int *currwW	=rows[0]-1*4*(1+NPREDS+1)+1;
			//	int *currNW	=rows[1]-1*4*(1+NPREDS+1)+1;
				int *currN	=rows[1]+0*4*(1+NPREDS+1)+1;
			//	int *currNE	=rows[1]+1*4*(1+NPREDS+1)+1;
				long long pred=currw[NPREDS];
				for(int k=0;k<NPREDS;++k)
				{
					int w=currw[k]+currN[k];
					pred+=(long long)w*preds[k];
				}
			//	int pred0=pred;
			//	pred*=currw[NPREDS+1];
				pred+=1LL<<(src->depth[kc]*2+4)>>1;
				pred>>=src->depth[kc]*2+4;
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
				int *currC=rows[0]+0*4*(1+NPREDS+1)+1;
				int e=curr-(int)pred;
				//int ex=e>>1, ey=e+(e>>1);
				//if(abs(curr-N)<abs(curr-W))
				//	ex=e+(e>>1), ey=e>>1;
			//	e*=currw[NPREDS+1];
				currw[NPREDS]+=e;
				currC[NPREDS]=currN[NPREDS]+e;
				for(int k=0;k<NPREDS;++k)
				{
					int e2=e*preds[k];
					currw[k]+=e2;
					currC[k]=currN[k]+e2;
				}
				//{
				//	currw[k]+=ex*preds[k];
				//	currC[k]=currN[k]+ey*preds[k];
				//}
			//	currw[NPREDS+1]+=(curr-(int)pred)*pred0;
			}
		}
	}
}
