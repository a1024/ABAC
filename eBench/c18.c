#include"ebench.h"
#include"c18.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<process.h>
#include<immintrin.h>
static const char file[]=__FILE__;

#define CODE_N(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_N	=N	##NBIDX[KC];\
		int pred=nb_N;\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_N)<<maxdepth|pred];\
	}while(0)
#define CODE_W(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_W	=W	##NBIDX[KC];\
		int pred=nb_W;\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_W)<<maxdepth|pred];\
	}while(0)
#define CODE_AV2(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_N	=N	##NBIDX[KC],\
			nb_W	=W	##NBIDX[KC];\
		int pred=(nb_N+nb_W)/2;\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_AV2)<<maxdepth|pred];\
	}while(0)
#define CODE_WG(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_NN	=NN	##NBIDX[KC],\
			nb_NNE	=NNE	##NBIDX[KC],\
			nb_NW	=NW	##NBIDX[KC],\
			nb_N	=N	##NBIDX[KC],\
			nb_NE	=NE	##NBIDX[KC],\
			nb_WW	=WW	##NBIDX[KC],\
			nb_W	=W	##NBIDX[KC];\
		int gx=abs(nb_W-nb_WW)+abs(nb_N-nb_NW)+abs(nb_NE-nb_N)+1;\
		int gy=abs(nb_W-nb_NW)+abs(nb_N-nb_NN)+abs(nb_NE-nb_NNE)+1;\
		int pred=(gx*nb_N+gy*nb_W)/(gx+gy);\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_WG)<<maxdepth|pred];\
	}while(0)
#define CODE_CG(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_NW	=NW	##NBIDX[KC],\
			nb_N	=N	##NBIDX[KC],\
			nb_W	=W	##NBIDX[KC];\
		int pred=nb_N+nb_W-nb_NW;\
		pred=MAXVAR(pred, vmin[KC]);\
		pred=MINVAR(pred, vmax[KC]);\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_CG)<<maxdepth|pred];\
	}while(0)
#define CODE_AV3(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_NW	=NW	##NBIDX[KC],\
			nb_N	=N	##NBIDX[KC],\
			nb_W	=W	##NBIDX[KC];\
		int pred=(3*(nb_N+nb_W)-2*nb_NW+2)>>2;\
		pred=MAXVAR(pred, vmin[KC]);\
		pred=MINVAR(pred, vmax[KC]);\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_AV3)<<maxdepth|pred];\
	}while(0)
#define CODE_AV4(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_NW	=NW	##NBIDX[KC],\
			nb_N	=N	##NBIDX[KC],\
			nb_NE	=NE	##NBIDX[KC],\
			nb_W	=W	##NBIDX[KC];\
		int pred=(4*(nb_N+nb_W)+nb_NE-nb_NW+4)>>3;\
		pred=MAXVAR(pred, vmin[KC]);\
		pred=MINVAR(pred, vmax[KC]);\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_AV4)<<maxdepth|pred];\
	}while(0)
#define CODE_AV5(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_NW	=NW	##NBIDX[KC],\
			nb_N	=N	##NBIDX[KC],\
			nb_NE	=NE	##NBIDX[KC],\
			nb_WW	=WW	##NBIDX[KC],\
			nb_W	=W	##NBIDX[KC];\
		int pred=nb_W+((5*(nb_N-nb_NW)+nb_NE-nb_WW+4)>>3);\
		pred=MAXVAR(pred, vmin[KC]);\
		pred=MINVAR(pred, vmax[KC]);\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_AV5)<<maxdepth|pred];\
	}while(0)
#define CODE_AV6(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_NN	=NN	##NBIDX[KC],\
			nb_NW	=NW	##NBIDX[KC],\
			nb_N	=N	##NBIDX[KC],\
			nb_NE	=NE	##NBIDX[KC],\
			nb_WW	=WW	##NBIDX[KC],\
			nb_W	=W	##NBIDX[KC];\
		int pred=nb_W+((6*nb_N-5*nb_NW+nb_NE-nb_NN-nb_WW+4)>>3);\
		pred=MAXVAR(pred, vmin[KC]);\
		pred=MINVAR(pred, vmax[KC]);\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_AV6)<<maxdepth|pred];\
	}while(0)
#define CODE_AV9(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_NNW	=NNW	##NBIDX[KC],\
			nb_NN	=NN	##NBIDX[KC],\
			nb_NNE	=NNE	##NBIDX[KC],\
			nb_NWW	=NWW	##NBIDX[KC],\
			nb_NW	=NW	##NBIDX[KC],\
			nb_N	=N	##NBIDX[KC],\
			nb_NE	=NE	##NBIDX[KC],\
			nb_WW	=WW	##NBIDX[KC],\
			nb_W	=W	##NBIDX[KC];\
		int pred=nb_W+((10*nb_N-9*nb_NW+4*nb_NE-2*(nb_NN+nb_WW)-nb_NNE+nb_NNW-nb_NWW+8)>>4);\
		pred=MAXVAR(pred, vmin[KC]);\
		pred=MINVAR(pred, vmax[KC]);\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_AV9)<<maxdepth|pred];\
	}while(0)
#define CODE_AVB(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_NNWW	=NNWW	##NBIDX[KC],\
			nb_NNW	=NNW	##NBIDX[KC],\
			nb_NN	=NN	##NBIDX[KC],\
			nb_NNE	=NNE	##NBIDX[KC],\
			nb_NWW	=NWW	##NBIDX[KC],\
			nb_NW	=NW	##NBIDX[KC],\
			nb_N	=N	##NBIDX[KC],\
			nb_NE	=NE	##NBIDX[KC],\
			nb_NEE	=NEE	##NBIDX[KC],\
			nb_WW	=WW	##NBIDX[KC],\
			nb_W	=W	##NBIDX[KC];\
		int pred=((\
			+0x04*nb_NNWW	+0x03*nb_NNW	-0x1F*nb_NN	-0x26*nb_NNE\
			+0x07*nb_NWW	-0x9E*nb_NW	+0xDB*nb_N	+0x1E*nb_NE	+0x13*nb_NEE\
			-0x2A*nb_WW	+0xF3*nb_W\
		+128)>>8);\
		pred=MAXVAR(pred, vmin[KC]);\
		pred=MINVAR(pred, vmax[KC]);\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_AVB)<<maxdepth|pred];\
	}while(0)

void c18_analyze(Image const *src, int x1, int x2, int y1, int y2, C18Info *info)
{
	double t_start=time_sec();
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int half[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int mask[]=
	{
		nlevels[0]-1,
		nlevels[1]-1,
		nlevels[2]-1,
		nlevels[3]-1,
	};
	int maxdepth=calc_maxdepth(src, 0), maxlevels=1<<maxdepth;
	int histsize=sizeof(int[OCH_COUNT*PRED_COUNT])<<maxdepth;
	int *hist=(int*)malloc(histsize);
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(hist, 0, histsize);
	int count=0;
	for(int ky=y1+2;ky<y2;++ky)//analysis loop
	{
		for(int kx=x1+2;kx<x2-2;++kx, ++count)
		{
			const int
				*NNWW0	=src->data+4*(src->iw*(ky-2LL)+kx-2LL),
				*NNW0	=src->data+4*(src->iw*(ky-2LL)+kx-1LL),
				*NN0	=src->data+4*(src->iw*(ky-2LL)+kx+0LL),
				*NNE0	=src->data+4*(src->iw*(ky-2LL)+kx+1LL),
				*NNEE0	=src->data+4*(src->iw*(ky-2LL)+kx+2LL),
				*NWW0	=src->data+4*(src->iw*(ky-1LL)+kx-2LL),
				*NW0	=src->data+4*(src->iw*(ky-1LL)+kx-1LL),
				*N0	=src->data+4*(src->iw*(ky-1LL)+kx+0LL),
				*NE0	=src->data+4*(src->iw*(ky-1LL)+kx+1LL),
				*NEE0	=src->data+4*(src->iw*(ky-1LL)+kx+2LL),
				*WW0	=src->data+4*(src->iw*(ky+0LL)+kx-2LL),
				*W0	=src->data+4*(src->iw*(ky+0LL)+kx-1LL),
				*target	=src->data+4*(src->iw*(ky+0LL)+kx+0LL),
				offset0	[]={0, 0, 0};
#define DECL_NB(NB) {NB##0[0]-NB##0[1], NB##0[1]-NB##0[2], NB##0[2]-NB##0[0]}
			int
				NNWW1	[]=DECL_NB(NNWW),
				NNW1	[]=DECL_NB(NNW),
				NN1	[]=DECL_NB(NN),
				NNE1	[]=DECL_NB(NNE),
				NNEE1	[]=DECL_NB(NNEE),
				NWW1	[]=DECL_NB(NWW),
				NW1	[]=DECL_NB(NW),
				N1	[]=DECL_NB(N),
				NE1	[]=DECL_NB(NE),
				NEE1	[]=DECL_NB(NEE),
				WW1	[]=DECL_NB(WW),
				W1	[]=DECL_NB(W),
				offset1	[]={target[1], target[2], target[0]};
#undef  DECL_NB
#define DECL_NB(NB) {NB##0[0]-(3*NB##0[1]+NB##0[2])/4, NB##0[1]-(3*NB##0[2]+NB##0[0])/4, NB##0[2]-(3*NB##0[0]+NB##0[1])/4}
			int
				NNWW2	[]=DECL_NB(NNWW),
				NNW2	[]=DECL_NB(NNW),
				NN2	[]=DECL_NB(NN),
				NNE2	[]=DECL_NB(NNE),
				NNEE2	[]=DECL_NB(NNEE),
				NWW2	[]=DECL_NB(NWW),
				NW2	[]=DECL_NB(NW),
				N2	[]=DECL_NB(N),
				NE2	[]=DECL_NB(NE),
				NEE2	[]=DECL_NB(NEE),
				WW2	[]=DECL_NB(WW),
				W2	[]=DECL_NB(W),
				offset2	[]={(3*target[1]+target[2])/4, (3*target[2]+target[0])/4, (3*target[0]+target[1])/4};
#undef  DECL_NB
#define DECL_NB(NB) {NB##0[0]-(NB##0[1]+NB##0[2])/2, NB##0[1]-(NB##0[2]+NB##0[0])/2, NB##0[2]-(NB##0[0]+NB##0[1])/2}
			int
				NNWW3	[]=DECL_NB(NNWW),
				NNW3	[]=DECL_NB(NNW),
				NN3	[]=DECL_NB(NN),
				NNE3	[]=DECL_NB(NNE),
				NNEE3	[]=DECL_NB(NNEE),
				NWW3	[]=DECL_NB(NWW),
				NW3	[]=DECL_NB(NW),
				N3	[]=DECL_NB(N),
				NE3	[]=DECL_NB(NE),
				NEE3	[]=DECL_NB(NEE),
				WW3	[]=DECL_NB(WW),
				W3	[]=DECL_NB(W),
				offset3	[]={(target[1]+target[2])/2, (target[2]+target[0])/2, (target[0]+target[1])/2};
#undef  DECL_NB
#define DECL_NB(NB) {NB##0[0]-(NB##0[1]+3*NB##0[2])/4, NB##0[1]-(NB##0[2]+3*NB##0[0])/4, NB##0[2]-(NB##0[0]+3*NB##0[1])/4}
			int
				NNWW4	[]=DECL_NB(NNWW),
				NNW4	[]=DECL_NB(NNW),
				NN4	[]=DECL_NB(NN),
				NNE4	[]=DECL_NB(NNE),
				NNEE4	[]=DECL_NB(NNEE),
				NWW4	[]=DECL_NB(NWW),
				NW4	[]=DECL_NB(NW),
				N4	[]=DECL_NB(N),
				NE4	[]=DECL_NB(NE),
				NEE4	[]=DECL_NB(NEE),
				WW4	[]=DECL_NB(WW),
				W4	[]=DECL_NB(W),
				offset4	[]={(target[1]+3*target[2])/4, (target[2]+3*target[0])/4, (target[0]+3*target[1])/4};
#undef  DECL_NB
#define DECL_NB(NB) {NB##0[0]-NB##0[2], NB##0[1]-NB##0[0], NB##0[2]-NB##0[1]}
			int
				NNWW5	[]=DECL_NB(NNWW),
				NNW5	[]=DECL_NB(NNW),
				NN5	[]=DECL_NB(NN),
				NNE5	[]=DECL_NB(NNE),
				NNEE5	[]=DECL_NB(NNEE),
				NWW5	[]=DECL_NB(NWW),
				NW5	[]=DECL_NB(NW),
				N5	[]=DECL_NB(N),
				NE5	[]=DECL_NB(NE),
				NEE5	[]=DECL_NB(NEE),
				WW5	[]=DECL_NB(WW),
				W5	[]=DECL_NB(W),
				offset5	[]={target[2], target[0], target[1]};
#undef  DECL_NB
			int vmin[3], vmax[3];
			
			CODE_N(0, 0, OCH_R);
			CODE_N(0, 1, OCH_G);
			CODE_N(0, 2, OCH_B);
			CODE_W(0, 0, OCH_R);
			CODE_W(0, 1, OCH_G);
			CODE_W(0, 2, OCH_B);
			CODE_AV2(0, 0, OCH_R);
			CODE_AV2(0, 1, OCH_G);
			CODE_AV2(0, 2, OCH_B);
			CODE_WG(0, 0, OCH_R);
			CODE_WG(0, 1, OCH_G);
			CODE_WG(0, 2, OCH_B);
			vmin[0]=MINVAR(N0[0], W0[0]);
			vmin[1]=MINVAR(N0[1], W0[1]);
			vmin[2]=MINVAR(N0[2], W0[2]);
			vmax[0]=MAXVAR(N0[0], W0[0]);
			vmax[1]=MAXVAR(N0[1], W0[1]);
			vmax[2]=MAXVAR(N0[2], W0[2]);
			CODE_CG(0, 0, OCH_R);
			CODE_CG(0, 1, OCH_G);
			CODE_CG(0, 2, OCH_B);
			vmin[0]=MINVAR(vmin[0], NE0[0]);
			vmin[1]=MINVAR(vmin[1], NE0[1]);
			vmin[2]=MINVAR(vmin[2], NE0[2]);
			vmax[0]=MAXVAR(vmax[0], NE0[0]);
			vmax[1]=MAXVAR(vmax[1], NE0[1]);
			vmax[2]=MAXVAR(vmax[2], NE0[2]);
			CODE_AV3(0, 0, OCH_R);
			CODE_AV3(0, 1, OCH_G);
			CODE_AV3(0, 2, OCH_B);
			CODE_AV4(0, 0, OCH_R);
			CODE_AV4(0, 1, OCH_G);
			CODE_AV4(0, 2, OCH_B);
			CODE_AV5(0, 0, OCH_R);
			CODE_AV5(0, 1, OCH_G);
			CODE_AV5(0, 2, OCH_B);
			CODE_AV6(0, 0, OCH_R);
			CODE_AV6(0, 1, OCH_G);
			CODE_AV6(0, 2, OCH_B);
			CODE_AV9(0, 0, OCH_R);
			CODE_AV9(0, 1, OCH_G);
			CODE_AV9(0, 2, OCH_B);
			CODE_AVB(0, 0, OCH_R);
			CODE_AVB(0, 1, OCH_G);
			CODE_AVB(0, 2, OCH_B);
			
			CODE_N(1, 0, OCH_RG);
			CODE_N(1, 1, OCH_GB);
			CODE_N(1, 2, OCH_BR);
			CODE_W(1, 0, OCH_RG);
			CODE_W(1, 1, OCH_GB);
			CODE_W(1, 2, OCH_BR);
			CODE_AV2(1, 0, OCH_RG);
			CODE_AV2(1, 1, OCH_GB);
			CODE_AV2(1, 2, OCH_BR);
			CODE_WG(1, 0, OCH_RG);
			CODE_WG(1, 1, OCH_GB);
			CODE_WG(1, 2, OCH_BR);
			vmin[0]=MINVAR(N1[0], W1[0]);
			vmin[1]=MINVAR(N1[1], W1[1]);
			vmin[2]=MINVAR(N1[2], W1[2]);
			vmax[0]=MAXVAR(N1[0], W1[0]);
			vmax[1]=MAXVAR(N1[1], W1[1]);
			vmax[2]=MAXVAR(N1[2], W1[2]);
			CODE_CG(1, 0, OCH_RG);
			CODE_CG(1, 1, OCH_GB);
			CODE_CG(1, 2, OCH_BR);
			vmin[0]=MINVAR(vmin[0], NE1[0]);
			vmin[1]=MINVAR(vmin[1], NE1[1]);
			vmin[2]=MINVAR(vmin[2], NE1[2]);
			vmax[0]=MAXVAR(vmax[0], NE1[0]);
			vmax[1]=MAXVAR(vmax[1], NE1[1]);
			vmax[2]=MAXVAR(vmax[2], NE1[2]);
			CODE_AV3(1, 0, OCH_RG);
			CODE_AV3(1, 1, OCH_GB);
			CODE_AV3(1, 2, OCH_BR);
			CODE_AV4(1, 0, OCH_RG);
			CODE_AV4(1, 1, OCH_GB);
			CODE_AV4(1, 2, OCH_BR);
			CODE_AV5(1, 0, OCH_RG);
			CODE_AV5(1, 1, OCH_GB);
			CODE_AV5(1, 2, OCH_BR);
			CODE_AV6(1, 0, OCH_RG);
			CODE_AV6(1, 1, OCH_GB);
			CODE_AV6(1, 2, OCH_BR);
			CODE_AV9(1, 0, OCH_RG);
			CODE_AV9(1, 1, OCH_GB);
			CODE_AV9(1, 2, OCH_BR);
			CODE_AVB(1, 0, OCH_RG);
			CODE_AVB(1, 1, OCH_GB);
			CODE_AVB(1, 2, OCH_BR);
			
			CODE_N(2, 0, OCH_R1);
			CODE_N(2, 1, OCH_G1);
			CODE_N(2, 2, OCH_B1);
			CODE_W(2, 0, OCH_R1);
			CODE_W(2, 1, OCH_G1);
			CODE_W(2, 2, OCH_B1);
			CODE_AV2(2, 0, OCH_R1);
			CODE_AV2(2, 1, OCH_G1);
			CODE_AV2(2, 2, OCH_B1);
			CODE_WG(2, 0, OCH_R1);
			CODE_WG(2, 1, OCH_G1);
			CODE_WG(2, 2, OCH_B1);
			vmin[0]=MINVAR(N2[0], W2[0]);
			vmin[1]=MINVAR(N2[1], W2[1]);
			vmin[2]=MINVAR(N2[2], W2[2]);
			vmax[0]=MAXVAR(N2[0], W2[0]);
			vmax[1]=MAXVAR(N2[1], W2[1]);
			vmax[2]=MAXVAR(N2[2], W2[2]);
			CODE_CG(2, 0, OCH_R1);
			CODE_CG(2, 1, OCH_G1);
			CODE_CG(2, 2, OCH_B1);
			vmin[0]=MINVAR(vmin[0], NE2[0]);
			vmin[1]=MINVAR(vmin[1], NE2[1]);
			vmin[2]=MINVAR(vmin[2], NE2[2]);
			vmax[0]=MAXVAR(vmax[0], NE2[0]);
			vmax[1]=MAXVAR(vmax[1], NE2[1]);
			vmax[2]=MAXVAR(vmax[2], NE2[2]);
			CODE_AV3(2, 0, OCH_R1);
			CODE_AV3(2, 1, OCH_G1);
			CODE_AV3(2, 2, OCH_B1);
			CODE_AV4(2, 0, OCH_R1);
			CODE_AV4(2, 1, OCH_G1);
			CODE_AV4(2, 2, OCH_B1);
			CODE_AV5(2, 0, OCH_R1);
			CODE_AV5(2, 1, OCH_G1);
			CODE_AV5(2, 2, OCH_B1);
			CODE_AV6(2, 0, OCH_R1);
			CODE_AV6(2, 1, OCH_G1);
			CODE_AV6(2, 2, OCH_B1);
			CODE_AV9(2, 0, OCH_R1);
			CODE_AV9(2, 1, OCH_G1);
			CODE_AV9(2, 2, OCH_B1);
			CODE_AVB(2, 0, OCH_R1);
			CODE_AVB(2, 1, OCH_G1);
			CODE_AVB(2, 2, OCH_B1);
			
			CODE_N(3, 0, OCH_R2);
			CODE_N(3, 1, OCH_G2);
			CODE_N(3, 2, OCH_B2);
			CODE_W(3, 0, OCH_R2);
			CODE_W(3, 1, OCH_G2);
			CODE_W(3, 2, OCH_B2);
			CODE_AV2(3, 0, OCH_R2);
			CODE_AV2(3, 1, OCH_G2);
			CODE_AV2(3, 2, OCH_B2);
			CODE_WG(3, 0, OCH_R2);
			CODE_WG(3, 1, OCH_G2);
			CODE_WG(3, 2, OCH_B2);
			vmin[0]=MINVAR(N3[0], W3[0]);
			vmin[1]=MINVAR(N3[1], W3[1]);
			vmin[2]=MINVAR(N3[2], W3[2]);
			vmax[0]=MAXVAR(N3[0], W3[0]);
			vmax[1]=MAXVAR(N3[1], W3[1]);
			vmax[2]=MAXVAR(N3[2], W3[2]);
			CODE_CG(3, 0, OCH_R2);
			CODE_CG(3, 1, OCH_G2);
			CODE_CG(3, 2, OCH_B2);
			vmin[0]=MINVAR(vmin[0], NE3[0]);
			vmin[1]=MINVAR(vmin[1], NE3[1]);
			vmin[2]=MINVAR(vmin[2], NE3[2]);
			vmax[0]=MAXVAR(vmax[0], NE3[0]);
			vmax[1]=MAXVAR(vmax[1], NE3[1]);
			vmax[2]=MAXVAR(vmax[2], NE3[2]);
			CODE_AV3(3, 0, OCH_R2);
			CODE_AV3(3, 1, OCH_G2);
			CODE_AV3(3, 2, OCH_B2);
			CODE_AV4(3, 0, OCH_R2);
			CODE_AV4(3, 1, OCH_G2);
			CODE_AV4(3, 2, OCH_B2);
			CODE_AV5(3, 0, OCH_R2);
			CODE_AV5(3, 1, OCH_G2);
			CODE_AV5(3, 2, OCH_B2);
			CODE_AV6(3, 0, OCH_R2);
			CODE_AV6(3, 1, OCH_G2);
			CODE_AV6(3, 2, OCH_B2);
			CODE_AV9(3, 0, OCH_R2);
			CODE_AV9(3, 1, OCH_G2);
			CODE_AV9(3, 2, OCH_B2);
			CODE_AVB(3, 0, OCH_R2);
			CODE_AVB(3, 1, OCH_G2);
			CODE_AVB(3, 2, OCH_B2);
			
			CODE_N(4, 0, OCH_R3);
			CODE_N(4, 1, OCH_G3);
			CODE_N(4, 2, OCH_B3);
			CODE_W(4, 0, OCH_R3);
			CODE_W(4, 1, OCH_G3);
			CODE_W(4, 2, OCH_B3);
			CODE_AV2(4, 0, OCH_R3);
			CODE_AV2(4, 1, OCH_G3);
			CODE_AV2(4, 2, OCH_B3);
			CODE_WG(4, 0, OCH_R3);
			CODE_WG(4, 1, OCH_G3);
			CODE_WG(4, 2, OCH_B3);
			vmin[0]=MINVAR(N4[0], W4[0]);
			vmin[1]=MINVAR(N4[1], W4[1]);
			vmin[2]=MINVAR(N4[2], W4[2]);
			vmax[0]=MAXVAR(N4[0], W4[0]);
			vmax[1]=MAXVAR(N4[1], W4[1]);
			vmax[2]=MAXVAR(N4[2], W4[2]);
			CODE_CG(4, 0, OCH_R3);
			CODE_CG(4, 1, OCH_G3);
			CODE_CG(4, 2, OCH_B3);
			vmin[0]=MINVAR(vmin[0], NE4[0]);
			vmin[1]=MINVAR(vmin[1], NE4[1]);
			vmin[2]=MINVAR(vmin[2], NE4[2]);
			vmax[0]=MAXVAR(vmax[0], NE4[0]);
			vmax[1]=MAXVAR(vmax[1], NE4[1]);
			vmax[2]=MAXVAR(vmax[2], NE4[2]);
			CODE_AV3(4, 0, OCH_R3);
			CODE_AV3(4, 1, OCH_G3);
			CODE_AV3(4, 2, OCH_B3);
			CODE_AV4(4, 0, OCH_R3);
			CODE_AV4(4, 1, OCH_G3);
			CODE_AV4(4, 2, OCH_B3);
			CODE_AV5(4, 0, OCH_R3);
			CODE_AV5(4, 1, OCH_G3);
			CODE_AV5(4, 2, OCH_B3);
			CODE_AV6(4, 0, OCH_R3);
			CODE_AV6(4, 1, OCH_G3);
			CODE_AV6(4, 2, OCH_B3);
			CODE_AV9(4, 0, OCH_R3);
			CODE_AV9(4, 1, OCH_G3);
			CODE_AV9(4, 2, OCH_B3);
			CODE_AVB(4, 0, OCH_R3);
			CODE_AVB(4, 1, OCH_G3);
			CODE_AVB(4, 2, OCH_B3);
			
			CODE_N(5, 0, OCH_RB);
			CODE_N(5, 1, OCH_GR);
			CODE_N(5, 2, OCH_BG);
			CODE_W(5, 0, OCH_RB);
			CODE_W(5, 1, OCH_GR);
			CODE_W(5, 2, OCH_BG);
			CODE_AV2(5, 0, OCH_RB);
			CODE_AV2(5, 1, OCH_GR);
			CODE_AV2(5, 2, OCH_BG);
			CODE_WG(5, 0, OCH_RB);
			CODE_WG(5, 1, OCH_GR);
			CODE_WG(5, 2, OCH_BG);
			vmin[0]=MINVAR(N5[0], W5[0]);
			vmin[1]=MINVAR(N5[1], W5[1]);
			vmin[2]=MINVAR(N5[2], W5[2]);
			vmax[0]=MAXVAR(N5[0], W5[0]);
			vmax[1]=MAXVAR(N5[1], W5[1]);
			vmax[2]=MAXVAR(N5[2], W5[2]);
			CODE_CG(5, 0, OCH_RB);
			CODE_CG(5, 1, OCH_GR);
			CODE_CG(5, 2, OCH_BG);
			vmin[0]=MINVAR(vmin[0], NE5[0]);
			vmin[1]=MINVAR(vmin[1], NE5[1]);
			vmin[2]=MINVAR(vmin[2], NE5[2]);
			vmax[0]=MAXVAR(vmax[0], NE5[0]);
			vmax[1]=MAXVAR(vmax[1], NE5[1]);
			vmax[2]=MAXVAR(vmax[2], NE5[2]);
			CODE_AV3(5, 0, OCH_RB);
			CODE_AV3(5, 1, OCH_GR);
			CODE_AV3(5, 2, OCH_BG);
			CODE_AV4(5, 0, OCH_RB);
			CODE_AV4(5, 1, OCH_GR);
			CODE_AV4(5, 2, OCH_BG);
			CODE_AV5(5, 0, OCH_RB);
			CODE_AV5(5, 1, OCH_GR);
			CODE_AV5(5, 2, OCH_BG);
			CODE_AV6(5, 0, OCH_RB);
			CODE_AV6(5, 1, OCH_GR);
			CODE_AV6(5, 2, OCH_BG);
			CODE_AV9(5, 0, OCH_RB);
			CODE_AV9(5, 1, OCH_GR);
			CODE_AV9(5, 2, OCH_BG);
			CODE_AVB(5, 0, OCH_RB);
			CODE_AVB(5, 1, OCH_GR);
			CODE_AVB(5, 2, OCH_BG);
		}
	}
	if(!count)
		return;
	double norm=1./count;
	//int res=(x2-x1-4)*(y2-y1-2);
	for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)
	{
		int *curr_hist=hist+((size_t)kc<<maxdepth);
		double e=0;
		for(int ks=0;ks<maxlevels;++ks)
		{
			int freq=curr_hist[ks];
			if(freq)
			{
				double p=freq*norm;
				e-=p*log2(p);
			}
		}
		info->esizes[kc]=e/8;
	}
	int predsel[OCH_COUNT]={0};
	for(int kc=0;kc<OCH_COUNT;++kc)//select best predictors
	{
		int bestpred=0;
		for(int kp=1;kp<PRED_COUNT;++kp)
		{
			if(info->esizes[kc*PRED_COUNT+bestpred]>info->esizes[kc*PRED_COUNT+kp])
				bestpred=kp;
		}
		predsel[kc]=bestpred;
	}
	double bestsize=0;
	for(int kt=0;kt<RCT_COUNT;++kt)//select best RCT
	{
		const unsigned char *group=rct_combinations[kt];
		double csize=
			+info->esizes[group[0]*PRED_COUNT+predsel[group[0]]]
			+info->esizes[group[1]*PRED_COUNT+predsel[group[1]]]
			+info->esizes[group[2]*PRED_COUNT+predsel[group[2]]];
		info->rctsizes[kt]=csize/3;
		if(!kt||bestsize>csize)
		{
			bestsize=csize;
			info->bestrct=kt;
		}
	}
	const unsigned char *group=rct_combinations[info->bestrct];
	info->predidx[0]=predsel[group[0]];
	info->predidx[1]=predsel[group[1]];
	info->predidx[2]=predsel[group[2]];
	free(hist);
	info->t_analysis=time_sec()-t_start;
}