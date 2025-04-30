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


//	#define PRINTINFO
//	#define PRINT_L1_BOUNDS


#define NPREDS 13
void pred_ols7(Image *src, int fwd)
{
#ifdef PRINTINFO
	static int silence=0;
#endif
#ifdef PRINT_L1_BOUNDS
	int cmin=0, cmax=0;
	int bmin=0, bmax=0;
#endif
	int weights[4][NPREDS+1]={0};//coeffs: 18 bit	bias: 15 bit
	int bufsize=(src->iw+8*2)*(int)sizeof(int[4*4]);//4 padded rows * 4 channels max
	int *pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
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
					NNW	=rows[2][-1*4],
					NN	=rows[2][+0*4],
					NNE	=rows[2][+1*4],
					NNEEE	=rows[2][+3*4],
					NWW	=rows[1][-2*4],
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
					N+W-NW,
					N+W-NW,
					N,
					W,
					3*(N-NN)+NNN,
					3*(W-WW)+WWW,	//west is more important that north in DIV2K & GDCC
					W+NE-N,
					(WWWWW+WW-W+NNN+N+NEEEEE)>>2,
					N+NE-NNE,
					W+NW-NWW,
					N+NW-NNW,
					NE+NEE-NNEEE,
					W+((NEEE+NEEEEE-N-W)>>3),

				//	N,
				//	W,
				//	3*(N-NN)+NNN,
				//	3*(W-WW)+WWW,
				//	W+NE-N,
				//	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4,
				//	N+W-NW,
				//	N+NE-NNE,
				};
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				CLAMP2(preds[0], vmin, vmax);
				int *currw=weights[kc];
				int pred=currw[NPREDS];
				//int wsum=0;
#ifdef PRINT_L1_BOUNDS
				if(bmin>currw[NPREDS])bmin=currw[NPREDS];
				if(bmax<currw[NPREDS])bmax=currw[NPREDS];
#endif
				for(int k=0;k<NPREDS;++k)
				{
					//wsum+=currw[k];
					pred+=currw[k]*preds[k];
#ifdef PRINT_L1_BOUNDS
					if(cmin>currw[k])cmin=currw[k];
					if(cmax<currw[k])cmax=currw[k];
#endif
				}
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");
	#define L1SH 20	//DIV2K
//	#define L1SH 20	//GDCC
//	#define L1SH 15	//synth
				//if(wsum)
				//	pred/=wsum;
				//const int sh[]={20, 20, 17, 16};
				//pred+=1<<sh[kc]>>1;
				//pred>>=sh[kc];
				pred+=1<<L1SH>>1;
				pred>>=L1SH;
				//int p0=pred;
				if(vmin>NW)vmin=NW;
				if(vmax<NW)vmax=NW;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);
#ifdef PRINTINFO
				if(!(ky&0x7F)&&kx==src->iw/2&&!silence&&!kc)
				{
					static const char *prednames[]=
					{
						"N+W-NW",
						"N+W-NW",
						"N",
						"W",
						"3*(N-NN)+NNN",
						"3*(W-WW)+WWW",
						"W+NE-N",
						"(WWWWW+WW-W+NNN+N+NEEEEE)>>2",
						"N+NE-NNE",
						"W+NW-NWW",
						"N+NW-NNW",
						"NE+NEE-NNEEE",
						"W+((NEEE+NEEEEE-N-W)>>3)",
						"bias",
					};
					static char msg[2048]={0};
					int printed=0;
					printed+=snprintf(msg+printed, sizeof(msg)-1-printed,
						"YXC %d %d %d%c  (%6.2lf%%)\n"
						, ky, kx, kc, "YUVA"[kc], 100.*ky/src->ih
					);
					for(int k=0;k<NPREDS+1;++k)
					{
						int nstars=abs(currw[k]/6144);
						printed+=snprintf(msg+printed, sizeof(msg)-1-printed,
							"%02d  %+08d  %s\n"
							, k
							, currw[k]
							, prednames[k]
						);
						for(int k2=0;k2<nstars;++k2)
							printed+=snprintf(msg+printed, sizeof(msg)-1-printed, "%c", currw[k]>0?'*':'v');
						printed+=snprintf(msg+printed, sizeof(msg)-1-printed, "\n\n");
					}
					//printed+=snprintf(msg+printed, sizeof(msg)-1-printed, "\n|");
					//for(int k=0;k<NPREDS;++k)
					//{
					//	int nstars=abs(currw[k]/32768);
					//	printed+=snprintf(msg+printed, sizeof(msg)-1-printed, "%02d", k);
					//	for(int k2=0;k2<nstars;++k2)
					//		printed+=snprintf(msg+printed, sizeof(msg)-1-printed, "%c", currw[k]>0?'*':'v');
					//	printed+=snprintf(msg+printed, sizeof(msg)-1-printed, "|");
					//}
					//printed+=snprintf(msg+printed, sizeof(msg)-1-printed, "\n");
					int choice=messagebox(MBOX_OKCANCEL, "Continue?", "%s", msg);
					if(choice)
						silence=1;
				}
#endif

				int curr=src->data[idx];
				int val=pred;
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

				//update
				int e=curr-pred;
				CLAMP2(e, 1-8, 15-8);
				static const int etable[]=
				{
					//   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
					-1, -1, -1, -1, -2, -2, -2, -2,  0, +2, +2, +2, +2, +1, +1, +1,
				//	-4, -4, -4, -5, -6, -8, -9, -6,  0, +6, +9, +8, +6, +5, +4, +4,
				};
				e=etable[e+8];
				//	... -8 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7
				//	... -1 -1 -1 -1 -2 -2 -2 -2  0 +2 +2 +2 +2 +1 +1 +1 ...
				//e=(((e>0)-(e<0))<<1)-((e>4)-(e<-4));
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
#ifdef PRINTINFO
	silence=0;
#endif
}
