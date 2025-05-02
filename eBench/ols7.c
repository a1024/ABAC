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

#if 1
#define BIAS0 0
#define NPREDS 13
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
	PRED( 50000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)
#endif
#if 0
#define BIAS0 0
#define NPREDS 8
#define PREDLIST\
	PRED(100000, N)\
	PRED(100000, W)\
	PRED( 90000, 3*(N-NN)+NNN)\
	PRED( 90000, 3*(W-WW)+WWW)\
	PRED( 80000, W+NE-N)\
	PRED( 80000, N+W-NW)\
	PRED( 70000, (WWWW+NNN+NEE+NEEEE)>>2)\
	PRED( 70000, N+NE-NNE)
#endif
#if 0
#define BIAS0 0
#define NPREDS 12
#define PREDLIST\
	PRED( 90000, N)\
	PRED( 90000, W)\
	PRED( 90000, 3*(N-NN)+NNN)\
	PRED( 60000, 3*(W-WW)+WWW)\
	PRED(160000, N+W-NW)\
	PRED( 90000, W+NE-N)\
	PRED(100000, N+NE-NNE)\
	PRED( 20000, W+NW-NWW)\
	PRED( 50000, N+NW-NNW)\
	PRED( 70000, NE+NEE-NNEEE)\
	PRED( 80000, W+((NEEE+NEEEEE-N-W)>>3))\
	PRED( 90000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)
#endif

void pred_ols7(Image *src, int fwd)
{
#ifdef PRINTINFO
	static int silence=0;
#endif
#ifdef PRINT_L1_BOUNDS
	long long cmin=0, cmax=0;
	long long bmin=0, bmax=0;
#endif
	int wsize=sizeof(long long[4][NPREDS+1]);
	long long *weights=(long long*)malloc(wsize);
	int bufsize=(src->iw+8*2)*(int)sizeof(int[6*4]);//6 padded rows * 4 channels max
	int *pixels=(int*)malloc(bufsize);
	if(!pixels||!weights)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	{
		int j=0;
#define PRED(W0, EXPR) weights[0*(NPREDS+1)+j]=weights[1*(NPREDS+1)+j]=weights[2*(NPREDS+1)+j]=weights[3*(NPREDS+1)+j]=W0; ++j;
		PREDLIST
#undef  PRED
		weights[0*(NPREDS+1)+NPREDS]=weights[1*(NPREDS+1)+NPREDS]=weights[2*(NPREDS+1)+NPREDS]=weights[3*(NPREDS+1)+NPREDS]=BIAS0;
	}
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL+6)%6)+8)*4-1,
			pixels+((src->iw+16LL)*((ky-1LL+6)%6)+8)*4-1,
			pixels+((src->iw+16LL)*((ky-2LL+6)%6)+8)*4-1,
			pixels+((src->iw+16LL)*((ky-3LL+6)%6)+8)*4-1,
			pixels+((src->iw+16LL)*((ky-4LL+6)%6)+8)*4-1,
			pixels+((src->iw+16LL)*((ky-5LL+6)%6)+8)*4-1,
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				++rows[0];
				++rows[1];
				++rows[2];
				++rows[3];
				++rows[4];
				++rows[5];
				if(!src->depth[kc])
					continue;
				int
					NNNNN		=rows[5][+0*4],
					NNNNWW		=rows[4][-2*4],
					NNNNW		=rows[4][-1*4],
					NNNN		=rows[4][+0*4],
					NNNNE		=rows[4][+1*4],
					NNNNEE		=rows[4][+2*4],
					NNNNEEEE	=rows[4][+4*4],
					NNNWWW		=rows[3][-3*4],
					NNNW		=rows[3][-1*4],
					NNN		=rows[3][+0*4],
					NNNE		=rows[3][+1*4],
					NNNEE		=rows[3][+2*4],
					NNNEEE		=rows[3][+3*4],
					NNNEEEE		=rows[3][+4*4],
					NNWWWW		=rows[2][-4*4],
					NNWWW		=rows[2][-3*4],
					NNWW		=rows[2][-2*4],
					NNW		=rows[2][-1*4],
					NN		=rows[2][+0*4],
					NNE		=rows[2][+1*4],
					NNEE		=rows[2][+2*4],
					NNEEE		=rows[2][+3*4],
					NNEEEE		=rows[2][+4*4],
					NWWWW		=rows[1][-4*4],
					NWWW		=rows[1][-3*4],
					NWW		=rows[1][-2*4],
					NW		=rows[1][-1*4],
					N		=rows[1][+0*4],
					NE		=rows[1][+1*4],
					NEE		=rows[1][+2*4],
					NEEE		=rows[1][+3*4],
					NEEEE		=rows[1][+4*4],
					NEEEEE		=rows[1][+5*4],
					NEEEEEE		=rows[1][+6*4],
					NEEEEEEE	=rows[1][+7*4],
					NEEEEEEEE	=rows[1][+8*4],
					WWWWWWWWW	=rows[0][-9*4],
					WWWWWWWW	=rows[0][-8*4],
					WWWWWWW		=rows[0][-7*4],
					WWWWWW		=rows[0][-6*4],
					WWWWW		=rows[0][-5*4],
					WWWW		=rows[0][-4*4],
					WWW		=rows[0][-3*4],
					WW		=rows[0][-2*4],
					W		=rows[0][-1*4];
				int preds[]=
				{
#define PRED(W0, EXPR) EXPR,
					PREDLIST
#undef  PRED

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
				long long *currw=weights+(NPREDS+1)*kc;
				long long pred=currw[NPREDS];
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
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");

#define L1SH 19
//13 preds:
//	#define L1SH 20	//DIV2K
//	#define L1SH 20	//GDCC
//	#define L1SH 15	//synth

//89 preds:
//	#define L1SH 26
				pred+=1<<L1SH>>1;
				pred>>=L1SH;
				long long p0=pred;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);
#ifdef PRINTINFO
			//	if(!(ky&0x7F)&&kx==src->iw/2&&!silence&&!kc)
				if(ky==src->ih-1&&kx==src->iw-1&&!silence&&!kc)
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

					copy_to_clipboard(msg, printed);//
					int choice=messagebox(MBOX_OKCANCEL, "Continue?", "%s", msg);
					if(choice)
						silence=1;
				}
#endif

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

				//update
				long long e=(curr>p0)-(curr<p0);
				currw[NPREDS]+=e;//bias
				for(int k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];//coeffs
			}
		}
	}
#ifdef PRINT_L1_BOUNDS
	if(loud_transforms)
		messagebox(MBOX_OK, "Info",
			"Coeffs %8lld ~ %8lld\n"
			"Bias   %8lld ~ %8lld\n"
			, cmin, cmax
			, bmin, bmax
		);
#endif
#ifdef PRINTINFO
	silence=0;
#endif
	free(pixels);
	free(weights);
}
