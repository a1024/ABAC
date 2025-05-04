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

#define NSTAGES 3
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

#if 1
#define NPREDS2 12
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
#if 0
#define NPREDS2 8
#define PREDLIST2\
	PRED(100000, eN)\
	PRED(100000, eW)\
	PRED(100000, 3*(eN-eNN)+eNNN)\
	PRED(100000, 3*(eW-eWW)+eWWW)\
	PRED(200000, eN+eW-eNW)\
	PRED(100000, eW+eNE-eN)\
	PRED(100000, eN+eNE-eNNE)\
	PRED( 50000, (eWWWWW+eWW-eW+eNNN+eN+eNEEEEE)>>2)
#endif

#if 0
#define NPREDS3 8
#define PREDLIST3\
	PRED(100000, e3N)\
	PRED(100000, e3W)\
	PRED(100000, 3*(e3N-e3NN)+e3NNN)\
	PRED(100000, 3*(e3W-e3WW)+e3WWW)\
	PRED(200000, e3N+e3W-e3NW)\
	PRED(100000, e3W+e3NE-e3N)\
	PRED(100000, e3N+e3NE-e3NNE)\
	PRED( 50000, (e3WWWWW+e3WW-e3W+e3NNN+e3N+e3NEEEEE)>>2)
#endif
#if 0
#define NPREDS3 12
#define PREDLIST3\
	PRED(100000, e3N+e3W-e3NW)\
	PRED(100000, e3N)\
	PRED(100000, e3W)\
	PRED(100000, 3*(e3N-e3NN)+e3NNN)\
	PRED(100000, 3*(e3W-e3WW)+e3WWW)\
	PRED(100000, e3W+e3NE-e3N)\
	PRED(100000, e3N+e3NE-e3NNE)\
	PRED(100000, e3W+((e3NEEE+e3NEEEEE-e3N-e3W)>>3))\
	PRED( 50000, e3W+e3NW-e3NWW)\
	PRED( 50000, e3N+e3NW-e3NNW)\
	PRED( 50000, e3NE+e3NEE-e3NNEEE)\
	PRED( 50000, (e3WWWWW+e3WW-e3W+e3NNN+e3N+e3NEEEEE)>>2)
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
	long long c2min=0, c2max=0;
	long long b2min=0, b2max=0;
	long long c3min=0, c3max=0;
	long long b3min=0, b3max=0;
#endif
	const int wsize=sizeof(long long[4][NPREDS+1]);
	long long *weights=(long long*)malloc(wsize);
	const int w2size=sizeof(long long[4][NPREDS2+1]);
	long long *weights2=(long long*)malloc(w2size);
	const int w3size=sizeof(long long[4][NPREDS2+1]);
	long long *weights3=(long long*)malloc(w2size);
	long long mixweights[4][NSTAGES]={0};
	int bufsize=(src->iw+8*2)*(int)sizeof(int[6*4*3]);//6 padded rows * 4 channels max * {pixels, residuals1, residuals2}
	int *pixels=(int*)malloc(bufsize);
	if(!pixels||!weights||!weights2||!weights3)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	//FILLMEM((long long*)mixweights, (1LL<<16)/NSTAGES, sizeof(mixweights), sizeof(long long));
	memset(pixels, 0, bufsize);
	{
		int j=0;
#define PRED(W0, EXPR) weights[0*(NPREDS+1)+j]=weights[1*(NPREDS+1)+j]=weights[2*(NPREDS+1)+j]=weights[3*(NPREDS+1)+j]=W0; ++j;
		PREDLIST
#undef  PRED
		weights[0*(NPREDS+1)+NPREDS]=weights[1*(NPREDS+1)+NPREDS]=weights[2*(NPREDS+1)+NPREDS]=weights[3*(NPREDS+1)+NPREDS]=BIAS0;
	}
	memset(weights2, 0, w2size);
	memset(weights3, 0, w3size);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+(((src->iw+16LL)*((ky-0LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-1LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-2LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-3LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-4LL+6)%6)+8)*4-1)*3,
			pixels+(((src->iw+16LL)*((ky-5LL+6)%6)+8)*4-1)*3,
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=3;
				rows[1]+=3;
				rows[2]+=3;
				rows[3]+=3;
				rows[4]+=3;
				rows[5]+=3;
				if(!src->depth[kc])
					continue;
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
				long long pred1=currw[NPREDS];
#ifdef PRINT_L1_BOUNDS
				if(bmin>currw[NPREDS])bmin=currw[NPREDS];
				if(bmax<currw[NPREDS])bmax=currw[NPREDS];
#endif
				for(int k=0;k<NPREDS;++k)
				{
					pred1+=currw[k]*preds[k];
#ifdef PRINT_L1_BOUNDS
					if(cmin>currw[k])cmin=currw[k];
					if(cmax<currw[k])cmax=currw[k];
#endif
				}
				int
					eNNNNN		=rows[5][1+0*4*3],
					eNNNNWW		=rows[4][1-2*4*3],
					eNNNNW		=rows[4][1-1*4*3],
					eNNNN		=rows[4][1+0*4*3],
					eNNNNE		=rows[4][1+1*4*3],
					eNNNNEE		=rows[4][1+2*4*3],
					eNNNNEEEE	=rows[4][1+4*4*3],
					eNNNWWW		=rows[3][1-3*4*3],
					eNNNW		=rows[3][1-1*4*3],
					eNNN		=rows[3][1+0*4*3],
					eNNNE		=rows[3][1+1*4*3],
					eNNNEE		=rows[3][1+2*4*3],
					eNNNEEE		=rows[3][1+3*4*3],
					eNNNEEEE	=rows[3][1+4*4*3],
					eNNWWWW		=rows[2][1-4*4*3],
					eNNWWW		=rows[2][1-3*4*3],
					eNNWW		=rows[2][1-2*4*3],
					eNNW		=rows[2][1-1*4*3],
					eNN		=rows[2][1+0*4*3],
					eNNE		=rows[2][1+1*4*3],
					eNNEE		=rows[2][1+2*4*3],
					eNNEEE		=rows[2][1+3*4*3],
					eNNEEEE		=rows[2][1+4*4*3],
					eNWWWW		=rows[1][1-4*4*3],
					eNWWW		=rows[1][1-3*4*3],
					eNWW		=rows[1][1-2*4*3],
					eNW		=rows[1][1-1*4*3],
					eN		=rows[1][1+0*4*3],
					eNE		=rows[1][1+1*4*3],
					eNEE		=rows[1][1+2*4*3],
					eNEEE		=rows[1][1+3*4*3],
					eNEEEE		=rows[1][1+4*4*3],
					eNEEEEE		=rows[1][1+5*4*3],
					eNEEEEEE	=rows[1][1+6*4*3],
					eNEEEEEEE	=rows[1][1+7*4*3],
					eNEEEEEEEE	=rows[1][1+8*4*3],
					eWWWWWWWWW	=rows[0][1-9*4*3],
					eWWWWWWWW	=rows[0][1-8*4*3],
					eWWWWWWW	=rows[0][1-7*4*3],
					eWWWWWW		=rows[0][1-6*4*3],
					eWWWWW		=rows[0][1-5*4*3],
					eWWWW		=rows[0][1-4*4*3],
					eWWW		=rows[0][1-3*4*3],
					eWW		=rows[0][1-2*4*3],
					eW		=rows[0][1-1*4*3];
				int preds2[]=
				{
#define PRED(W0, EXPR) EXPR,
					PREDLIST2
#undef  PRED
				};
				long long *currw2=weights2+(NPREDS2+1)*kc;
				long long pred2=currw2[NPREDS2];
#ifdef PRINT_L1_BOUNDS
				if(b2min>currw2[NPREDS2])b2min=currw2[NPREDS2];
				if(b2max<currw2[NPREDS2])b2max=currw2[NPREDS2];
#endif
				for(int k=0;k<NPREDS2;++k)
				{
					pred2+=currw2[k]*preds2[k];
#ifdef PRINT_L1_BOUNDS
					if(c2min>currw2[k])c2min=currw2[k];
					if(c2max<currw2[k])c2max=currw2[k];
#endif
				}
#if 0
				int
					e3NNNNN		=rows[5][2+0*4*3],
					e3NNNNWW	=rows[4][2-2*4*3],
					e3NNNNW		=rows[4][2-1*4*3],
					e3NNNN		=rows[4][2+0*4*3],
					e3NNNNE		=rows[4][2+1*4*3],
					e3NNNNEE	=rows[4][2+2*4*3],
					e3NNNNEEEE	=rows[4][2+4*4*3],
					e3NNNWWW	=rows[3][2-3*4*3],
					e3NNNW		=rows[3][2-1*4*3],
					e3NNN		=rows[3][2+0*4*3],
					e3NNNE		=rows[3][2+1*4*3],
					e3NNNEE		=rows[3][2+2*4*3],
					e3NNNEEE	=rows[3][2+3*4*3],
					e3NNNEEEE	=rows[3][2+4*4*3],
					e3NNWWWW	=rows[2][2-4*4*3],
					e3NNWWW		=rows[2][2-3*4*3],
					e3NNWW		=rows[2][2-2*4*3],
					e3NNW		=rows[2][2-1*4*3],
					e3NN		=rows[2][2+0*4*3],
					e3NNE		=rows[2][2+1*4*3],
					e3NNEE		=rows[2][2+2*4*3],
					e3NNEEE		=rows[2][2+3*4*3],
					e3NNEEEE	=rows[2][2+4*4*3],
					e3NWWWW		=rows[1][2-4*4*3],
					e3NWWW		=rows[1][2-3*4*3],
					e3NWW		=rows[1][2-2*4*3],
					e3NW		=rows[1][2-1*4*3],
					e3N		=rows[1][2+0*4*3],
					e3NE		=rows[1][2+1*4*3],
					e3NEE		=rows[1][2+2*4*3],
					e3NEEE		=rows[1][2+3*4*3],
					e3NEEEE		=rows[1][2+4*4*3],
					e3NEEEEE	=rows[1][2+5*4*3],
					e3NEEEEEE	=rows[1][2+6*4*3],
					e3NEEEEEEE	=rows[1][2+7*4*3],
					e3NEEEEEEEE	=rows[1][2+8*4*3],
					e3WWWWWWWWW	=rows[0][2-9*4*3],
					e3WWWWWWWW	=rows[0][2-8*4*3],
					e3WWWWWWW	=rows[0][2-7*4*3],
					e3WWWWWW	=rows[0][2-6*4*3],
					e3WWWWW		=rows[0][2-5*4*3],
					e3WWWW		=rows[0][2-4*4*3],
					e3WWW		=rows[0][2-3*4*3],
					e3WW		=rows[0][2-2*4*3],
					e3W		=rows[0][2-1*4*3];
				int preds3[]=
				{
#define PRED(W0, EXPR) EXPR,
					PREDLIST3
#undef  PRED
				};
				long long *currw3=weights3+(NPREDS2+1)*kc;
				long long pred3=currw3[NPREDS2];
#ifdef PRINT_L1_BOUNDS
				if(b3min>currw2[NPREDS2])b3min=currw2[NPREDS2];
				if(b3max<currw2[NPREDS2])b3max=currw2[NPREDS2];
#endif
				for(int k=0;k<NPREDS3;++k)
				{
					pred3+=currw3[k]*preds3[k];
#ifdef PRINT_L1_BOUNDS
					if(c3min>currw2[k])c3min=currw2[k];
					if(c3max<currw2[k])c3max=currw2[k];
#endif
				}
#endif
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");

#define L1SH 19
//13 preds:
//	#define L1SH 20	//DIV2K
//	#define L1SH 20	//GDCC
//	#define L1SH 15	//synth

//89 preds:
//	#define L1SH 26
				pred1+=pred2;
				pred1+=1<<L1SH>>1;
				pred1>>=L1SH;
				pred2+=1<<L1SH>>1;
				pred2>>=L1SH;
#if 0
				pred3+=1<<L1SH>>1;
				pred3>>=L1SH;

				//long long pred2s=pred1+pred2, pred3s=pred2s+pred3,
				//	predfs=(mixweights[kc][0]*pred1+mixweights[kc][1]*(pred1+pred2)+mixweights[kc][2]*(pred1+pred2+pred3))>>16,
				//	predc=predfs;

				//long long pred2s=pred1+((pred2-pred1)*mixweights[kc][0]>>20);
				//long long pred3s=pred2s+((pred3-pred2s)*mixweights[kc][1]>>20);
				//long long predc=pred3s;

				long long pred2s=pred1+pred2, pred3s=pred2s+pred3, predc=pred3s;
#endif
				//long long pred2s=(mixweights[kc][0]*pred1+mixweights[kc][1]*pred2)>>16, predc=pred2s;

				//long long pred1s=pred1*mixweights[kc][0]>>12;
				//long long pred2s=pred1s+((pred2-pred1s)*mixweights[kc][1]>>12);
				//long long predc=pred2s;

				long long pred2s=pred1+pred2, predc=pred2s;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				//if(vmin>NEEEE)vmin=NEEEE;
				//if(vmax<NEEEE)vmax=NEEEE;
				//if(vmin>NW)vmin=NW;
				//if(vmax<NW)vmax=NW;
				//if(vmin>WWW)vmin=WWW;
				//if(vmax<WWW)vmax=WWW;
				CLAMP2(predc, vmin, vmax);
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
				int val=(int)predc;
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
				rows[0][1]=curr-(int)pred1;
				rows[0][2]=curr-(int)pred2s;

				//update
				int e=(curr>pred1)-(curr<pred1);
				currw[NPREDS]+=e;//bias
				for(int k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];//coeffs

				//e=(curr>pred1s)-(curr<pred1s);
				e=(curr>pred2s)-(curr<pred2s);
				currw2[NPREDS2]+=e;//bias
				for(int k=0;k<NPREDS2;++k)
					currw2[k]+=e*preds2[k];//coeffs

				//mixweights[kc][0]+=e*(pred2-pred1);
#if 0
				e=(curr>pred3s)-(curr<pred3s);
				currw3[NPREDS3]+=e;//bias
				for(int k=0;k<NPREDS3;++k)
					currw3[k]+=e*preds3[k];//coeffs
#endif
				//mixweights[kc][1]+=e*(pred3-pred2s);

				//mixweights[kc][0]+=((curr>pred1s)-(curr<pred1s))*pred1;
				//mixweights[kc][1]+=((curr>pred2s)-(curr<pred2s))*pred2;

				//e=(curr>predfs)-(curr<predfs);
				//mixweights[kc][0]+=e*pred1;
				//mixweights[kc][1]+=e*pred2;
				//mixweights[kc][2]+=e*pred3;
			}
		}
	}
#ifdef PRINT_L1_BOUNDS
	if(loud_transforms)
		messagebox(MBOX_OK, "Info",
			"Coeffs  %8lld ~ %8lld\n"
			"Bias    %8lld ~ %8lld\n"
			"Coeffs2 %8lld ~ %8lld\n"
			"Bias2   %8lld ~ %8lld\n"
			, cmin, cmax
			, bmin, bmax
			, c2min, c2max
			, b2min, b2max
		);
#endif
#ifdef PRINTINFO
	silence=0;
#endif
	free(pixels);
	free(weights);
	free(weights2);
	free(weights3);
}
