#include"ebench.h"
#include"c18.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>

#define CODE_ZERO(NBIDX, KC, OCHIDX)\
	do\
	{\
		int pred=0;\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_ZERO)<<maxdepth|pred];\
	}while(0)
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
		CLAMP2(pred, vmin[KC], vmax[KC]);\
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
		CLAMP2(pred, vmin[KC], vmax[KC]);\
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
		CLAMP2(pred, vmin[KC], vmax[KC]);\
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
		CLAMP2(pred, vmin[KC], vmax[KC]);\
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
		CLAMP2(pred, vmin[KC], vmax[KC]);\
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
		CLAMP2(pred, vmin[KC], vmax[KC]);\
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
		CLAMP2(pred, vmin[KC], vmax[KC]);\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -half[KC], half[KC]-1);\
		pred=(target[KC]-pred+half[KC])&mask[KC];\
		++hist[(OCHIDX*PRED_COUNT+PRED_AVB)<<maxdepth|pred];\
	}while(0)

void c18_analyze(Image const *src, int x1, int x2, int y1, int y2, C18Info *info)
{
	double t_start=time_sec();
	int nlevels[4], half[4], mask[4], maxdepth, maxlevels, histsize, *hist, count, kx, ky;
	const int32_t *dataptr;
	double norm, csize, bestsize;
	int predsel[OCH2_COUNT]={0}, kc, kt, ks, kp;
	const unsigned char *group;

	dataptr=(const int32_t*)HDR_DATC(src);
	nlevels[0]=1<<src->depth[0];
	nlevels[1]=1<<src->depth[1];
	nlevels[2]=1<<src->depth[2];
	nlevels[3]=1<<src->depth[3];
	half[0]=nlevels[0]>>1;
	half[1]=nlevels[1]>>1;
	half[2]=nlevels[2]>>1;
	half[3]=nlevels[3]>>1;
	mask[0]=nlevels[0]-1;
	mask[1]=nlevels[1]-1;
	mask[2]=nlevels[2]-1;
	mask[3]=nlevels[3]-1;
	maxdepth=calc_maxdepth(src, 0);
	maxlevels=1<<maxdepth;
	histsize=sizeof(int[OCH2_COUNT*PRED_COUNT])<<maxdepth;
	hist=(int*)malloc(histsize);
	if(!hist)
	{
		CRASH("Alloc error");
		return;
	}
	memset(hist, 0, histsize);
	count=0;
	for(ky=y1+2;ky<y2;++ky)/*analysis loop*/
	{
		for(kx=x1+2;kx<x2-2;++kx, ++count)
		{
			const int
				*NNWW0	=dataptr+4*(src->iw*(ky-(intptr_t)2)+kx-(intptr_t)2),
				*NNW0	=dataptr+4*(src->iw*(ky-(intptr_t)2)+kx-(intptr_t)1),
				*NN0	=dataptr+4*(src->iw*(ky-(intptr_t)2)+kx+(intptr_t)0),
				*NNE0	=dataptr+4*(src->iw*(ky-(intptr_t)2)+kx+(intptr_t)1),
			/*	*NNEE0	=dataptr+4*(src->iw*(ky-(intptr_t)2)+kx+(intptr_t)2),*/
				*NWW0	=dataptr+4*(src->iw*(ky-(intptr_t)1)+kx-(intptr_t)2),
				*NW0	=dataptr+4*(src->iw*(ky-(intptr_t)1)+kx-(intptr_t)1),
				*N0	=dataptr+4*(src->iw*(ky-(intptr_t)1)+kx+(intptr_t)0),
				*NE0	=dataptr+4*(src->iw*(ky-(intptr_t)1)+kx+(intptr_t)1),
				*NEE0	=dataptr+4*(src->iw*(ky-(intptr_t)1)+kx+(intptr_t)2),
				*WW0	=dataptr+4*(src->iw*(ky+(intptr_t)0)+kx-(intptr_t)2),
				*W0	=dataptr+4*(src->iw*(ky+(intptr_t)0)+kx-(intptr_t)1),
				*target	=dataptr+4*(src->iw*(ky+(intptr_t)0)+kx+(intptr_t)0);
			int offset0	[3]={0};
			int
				NNWW1	[3]={0},
				NNW1	[3]={0},
				NN1	[3]={0},
				NNE1	[3]={0},
				NWW1	[3]={0},
				NW1	[3]={0},
				N1	[3]={0},
				NE1	[3]={0},
				NEE1	[3]={0},
				WW1	[3]={0},
				W1	[3]={0},
				offset1	[3]={0};
			int
				NNWW2	[3]={0},
				NNW2	[3]={0},
				NN2	[3]={0},
				NNE2	[3]={0},
				NWW2	[3]={0},
				NW2	[3]={0},
				N2	[3]={0},
				NE2	[3]={0},
				NEE2	[3]={0},
				WW2	[3]={0},
				W2	[3]={0},
				offset2	[3]={0};
			int
				NNWW3	[3]={0},
				NNW3	[3]={0},
				NN3	[3]={0},
				NNE3	[3]={0},
				NWW3	[3]={0},
				NW3	[3]={0},
				N3	[3]={0},
				NE3	[3]={0},
				NEE3	[3]={0},
				WW3	[3]={0},
				W3	[3]={0},
				offset3	[3]={0};
			int
				NNWW4	[3]={0},
				NNW4	[3]={0},
				NN4	[3]={0},
				NNE4	[3]={0},
				NWW4	[3]={0},
				NW4	[3]={0},
				N4	[3]={0},
				NE4	[3]={0},
				NEE4	[3]={0},
				WW4	[3]={0},
				W4	[3]={0},
				offset4	[3]={0};
			int
				NNWW5	[3]={0},
				NNW5	[3]={0},
				NN5	[3]={0},
				NNE5	[3]={0},
				NWW5	[3]={0},
				NW5	[3]={0},
				N5	[3]={0},
				NE5	[3]={0},
				NEE5	[3]={0},
				WW5	[3]={0},
				W5	[3]={0},
				offset5	[3]={0};
			int vmin[3], vmax[3];

#define DECL_NB(NB) NB##1[0]=NB##0[0]-NB##0[1], NB##1[1]=NB##0[1]-NB##0[2], NB##1[2]=NB##0[2]-NB##0[0]
			DECL_NB(NNWW);
			DECL_NB(NNW);
			DECL_NB(NN);
			DECL_NB(NNE);
		/*	DECL_NB(NNEE);*/
			DECL_NB(NWW);
			DECL_NB(NW);
			DECL_NB(N);
			DECL_NB(NE);
			DECL_NB(NEE);
			DECL_NB(WW);
			DECL_NB(W);
			offset1[0]=target[1], offset1[1]=target[2], offset1[2]=target[0];
#undef  DECL_NB
#define DECL_NB(NB) NB##2[0]=NB##0[0]-(3*NB##0[1]+NB##0[2])/4, NB##2[1]=NB##0[1]-(3*NB##0[2]+NB##0[0])/4, NB##2[2]=NB##0[2]-(3*NB##0[0]+NB##0[1])/4
			DECL_NB(NNWW);
			DECL_NB(NNW);
			DECL_NB(NN);
			DECL_NB(NNE);
		/*	DECL_NB(NNEE);*/
			DECL_NB(NWW);
			DECL_NB(NW);
			DECL_NB(N);
			DECL_NB(NE);
			DECL_NB(NEE);
			DECL_NB(WW);
			DECL_NB(W);
			offset2[0]=(3*target[1]+target[2])/4;
			offset2[1]=(3*target[2]+target[0])/4;
			offset2[2]=(3*target[0]+target[1])/4;
#undef  DECL_NB
#define DECL_NB(NB) NB##3[0]=NB##0[0]-(NB##0[1]+NB##0[2])/2, NB##3[1]=NB##0[1]-(NB##0[2]+NB##0[0])/2, NB##3[2]=NB##0[2]-(NB##0[0]+NB##0[1])/2
			DECL_NB(NNWW);
			DECL_NB(NNW);
			DECL_NB(NN);
			DECL_NB(NNE);
		/*	DECL_NB(NNEE);*/
			DECL_NB(NWW);
			DECL_NB(NW);
			DECL_NB(N);
			DECL_NB(NE);
			DECL_NB(NEE);
			DECL_NB(WW);
			DECL_NB(W);
			offset3[0]=(target[1]+target[2])/2;
			offset3[1]=(target[2]+target[0])/2;
			offset3[2]=(target[0]+target[1])/2;
#undef  DECL_NB
#define DECL_NB(NB) NB##4[0]=NB##0[0]-(NB##0[1]+3*NB##0[2])/4, NB##4[1]=NB##0[1]-(NB##0[2]+3*NB##0[0])/4, NB##4[2]=NB##0[2]-(NB##0[0]+3*NB##0[1])/4
			DECL_NB(NNWW);
			DECL_NB(NNW);
			DECL_NB(NN);
			DECL_NB(NNE);
		/*	DECL_NB(NNEE);*/
			DECL_NB(NWW);
			DECL_NB(NW);
			DECL_NB(N);
			DECL_NB(NE);
			DECL_NB(NEE);
			DECL_NB(WW);
			DECL_NB(W);
			offset4[0]=(target[1]+3*target[2])/4;
			offset4[1]=(target[2]+3*target[0])/4;
			offset4[2]=(target[0]+3*target[1])/4;
#undef  DECL_NB
#define DECL_NB(NB) NB##5[0]=NB##0[0]-NB##0[2], NB##5[1]=NB##0[1]-NB##0[0], NB##5[2]=NB##0[2]-NB##0[1]
			DECL_NB(NNWW);
			DECL_NB(NNW);
			DECL_NB(NN);
			DECL_NB(NNE);
		/*	DECL_NB(NNEE);*/
			DECL_NB(NWW);
			DECL_NB(NW);
			DECL_NB(N);
			DECL_NB(NE);
			DECL_NB(NEE);
			DECL_NB(WW);
			DECL_NB(W);
			offset5[0]=target[2];
			offset5[1]=target[0];
			offset5[2]=target[1];
#undef  DECL_NB

			CODE_ZERO	(0, 0, OCH2_R);
			CODE_ZERO	(0, 1, OCH2_G);
			CODE_ZERO	(0, 2, OCH2_B);
			CODE_N		(0, 0, OCH2_R);
			CODE_N		(0, 1, OCH2_G);
			CODE_N		(0, 2, OCH2_B);
			CODE_W		(0, 0, OCH2_R);
			CODE_W		(0, 1, OCH2_G);
			CODE_W		(0, 2, OCH2_B);
			CODE_AV2	(0, 0, OCH2_R);
			CODE_AV2	(0, 1, OCH2_G);
			CODE_AV2	(0, 2, OCH2_B);
			CODE_WG		(0, 0, OCH2_R);
			CODE_WG		(0, 1, OCH2_G);
			CODE_WG		(0, 2, OCH2_B);
			vmin[0]=MINVAR(N0[0], W0[0]);
			vmin[1]=MINVAR(N0[1], W0[1]);
			vmin[2]=MINVAR(N0[2], W0[2]);
			vmax[0]=MAXVAR(N0[0], W0[0]);
			vmax[1]=MAXVAR(N0[1], W0[1]);
			vmax[2]=MAXVAR(N0[2], W0[2]);
			CODE_CG		(0, 0, OCH2_R);
			CODE_CG		(0, 1, OCH2_G);
			CODE_CG		(0, 2, OCH2_B);
			vmin[0]=MINVAR(vmin[0], NE0[0]);
			vmin[1]=MINVAR(vmin[1], NE0[1]);
			vmin[2]=MINVAR(vmin[2], NE0[2]);
			vmax[0]=MAXVAR(vmax[0], NE0[0]);
			vmax[1]=MAXVAR(vmax[1], NE0[1]);
			vmax[2]=MAXVAR(vmax[2], NE0[2]);
			CODE_AV3	(0, 0, OCH2_R);
			CODE_AV3	(0, 1, OCH2_G);
			CODE_AV3	(0, 2, OCH2_B);
			CODE_AV4	(0, 0, OCH2_R);
			CODE_AV4	(0, 1, OCH2_G);
			CODE_AV4	(0, 2, OCH2_B);
			CODE_AV5	(0, 0, OCH2_R);
			CODE_AV5	(0, 1, OCH2_G);
			CODE_AV5	(0, 2, OCH2_B);
			CODE_AV6	(0, 0, OCH2_R);
			CODE_AV6	(0, 1, OCH2_G);
			CODE_AV6	(0, 2, OCH2_B);
			CODE_AV9	(0, 0, OCH2_R);
			CODE_AV9	(0, 1, OCH2_G);
			CODE_AV9	(0, 2, OCH2_B);
			CODE_AVB	(0, 0, OCH2_R);
			CODE_AVB	(0, 1, OCH2_G);
			CODE_AVB	(0, 2, OCH2_B);
			
			CODE_ZERO	(1, 0, OCH2_RG);
			CODE_ZERO	(1, 1, OCH2_GB);
			CODE_ZERO	(1, 2, OCH2_BR);
			CODE_N		(1, 0, OCH2_RG);
			CODE_N		(1, 1, OCH2_GB);
			CODE_N		(1, 2, OCH2_BR);
			CODE_W		(1, 0, OCH2_RG);
			CODE_W		(1, 1, OCH2_GB);
			CODE_W		(1, 2, OCH2_BR);
			CODE_AV2	(1, 0, OCH2_RG);
			CODE_AV2	(1, 1, OCH2_GB);
			CODE_AV2	(1, 2, OCH2_BR);
			CODE_WG		(1, 0, OCH2_RG);
			CODE_WG		(1, 1, OCH2_GB);
			CODE_WG		(1, 2, OCH2_BR);
			vmin[0]=MINVAR(N1[0], W1[0]);
			vmin[1]=MINVAR(N1[1], W1[1]);
			vmin[2]=MINVAR(N1[2], W1[2]);
			vmax[0]=MAXVAR(N1[0], W1[0]);
			vmax[1]=MAXVAR(N1[1], W1[1]);
			vmax[2]=MAXVAR(N1[2], W1[2]);
			CODE_CG		(1, 0, OCH2_RG);
			CODE_CG		(1, 1, OCH2_GB);
			CODE_CG		(1, 2, OCH2_BR);
			vmin[0]=MINVAR(vmin[0], NE1[0]);
			vmin[1]=MINVAR(vmin[1], NE1[1]);
			vmin[2]=MINVAR(vmin[2], NE1[2]);
			vmax[0]=MAXVAR(vmax[0], NE1[0]);
			vmax[1]=MAXVAR(vmax[1], NE1[1]);
			vmax[2]=MAXVAR(vmax[2], NE1[2]);
			CODE_AV3	(1, 0, OCH2_RG);
			CODE_AV3	(1, 1, OCH2_GB);
			CODE_AV3	(1, 2, OCH2_BR);
			CODE_AV4	(1, 0, OCH2_RG);
			CODE_AV4	(1, 1, OCH2_GB);
			CODE_AV4	(1, 2, OCH2_BR);
			CODE_AV5	(1, 0, OCH2_RG);
			CODE_AV5	(1, 1, OCH2_GB);
			CODE_AV5	(1, 2, OCH2_BR);
			CODE_AV6	(1, 0, OCH2_RG);
			CODE_AV6	(1, 1, OCH2_GB);
			CODE_AV6	(1, 2, OCH2_BR);
			CODE_AV9	(1, 0, OCH2_RG);
			CODE_AV9	(1, 1, OCH2_GB);
			CODE_AV9	(1, 2, OCH2_BR);
			CODE_AVB	(1, 0, OCH2_RG);
			CODE_AVB	(1, 1, OCH2_GB);
			CODE_AVB	(1, 2, OCH2_BR);
			
			CODE_ZERO	(2, 0, OCH2_R1);
			CODE_ZERO	(2, 1, OCH2_G1);
			CODE_ZERO	(2, 2, OCH2_B1);
			CODE_N		(2, 0, OCH2_R1);
			CODE_N		(2, 1, OCH2_G1);
			CODE_N		(2, 2, OCH2_B1);
			CODE_W		(2, 0, OCH2_R1);
			CODE_W		(2, 1, OCH2_G1);
			CODE_W		(2, 2, OCH2_B1);
			CODE_AV2	(2, 0, OCH2_R1);
			CODE_AV2	(2, 1, OCH2_G1);
			CODE_AV2	(2, 2, OCH2_B1);
			CODE_WG		(2, 0, OCH2_R1);
			CODE_WG		(2, 1, OCH2_G1);
			CODE_WG		(2, 2, OCH2_B1);
			vmin[0]=MINVAR(N2[0], W2[0]);
			vmin[1]=MINVAR(N2[1], W2[1]);
			vmin[2]=MINVAR(N2[2], W2[2]);
			vmax[0]=MAXVAR(N2[0], W2[0]);
			vmax[1]=MAXVAR(N2[1], W2[1]);
			vmax[2]=MAXVAR(N2[2], W2[2]);
			CODE_CG		(2, 0, OCH2_R1);
			CODE_CG		(2, 1, OCH2_G1);
			CODE_CG		(2, 2, OCH2_B1);
			vmin[0]=MINVAR(vmin[0], NE2[0]);
			vmin[1]=MINVAR(vmin[1], NE2[1]);
			vmin[2]=MINVAR(vmin[2], NE2[2]);
			vmax[0]=MAXVAR(vmax[0], NE2[0]);
			vmax[1]=MAXVAR(vmax[1], NE2[1]);
			vmax[2]=MAXVAR(vmax[2], NE2[2]);
			CODE_AV3	(2, 0, OCH2_R1);
			CODE_AV3	(2, 1, OCH2_G1);
			CODE_AV3	(2, 2, OCH2_B1);
			CODE_AV4	(2, 0, OCH2_R1);
			CODE_AV4	(2, 1, OCH2_G1);
			CODE_AV4	(2, 2, OCH2_B1);
			CODE_AV5	(2, 0, OCH2_R1);
			CODE_AV5	(2, 1, OCH2_G1);
			CODE_AV5	(2, 2, OCH2_B1);
			CODE_AV6	(2, 0, OCH2_R1);
			CODE_AV6	(2, 1, OCH2_G1);
			CODE_AV6	(2, 2, OCH2_B1);
			CODE_AV9	(2, 0, OCH2_R1);
			CODE_AV9	(2, 1, OCH2_G1);
			CODE_AV9	(2, 2, OCH2_B1);
			CODE_AVB	(2, 0, OCH2_R1);
			CODE_AVB	(2, 1, OCH2_G1);
			CODE_AVB	(2, 2, OCH2_B1);
			
			CODE_ZERO	(3, 0, OCH2_R2);
			CODE_ZERO	(3, 1, OCH2_G2);
			CODE_ZERO	(3, 2, OCH2_B2);
			CODE_N		(3, 0, OCH2_R2);
			CODE_N		(3, 1, OCH2_G2);
			CODE_N		(3, 2, OCH2_B2);
			CODE_W		(3, 0, OCH2_R2);
			CODE_W		(3, 1, OCH2_G2);
			CODE_W		(3, 2, OCH2_B2);
			CODE_AV2	(3, 0, OCH2_R2);
			CODE_AV2	(3, 1, OCH2_G2);
			CODE_AV2	(3, 2, OCH2_B2);
			CODE_WG		(3, 0, OCH2_R2);
			CODE_WG		(3, 1, OCH2_G2);
			CODE_WG		(3, 2, OCH2_B2);
			vmin[0]=MINVAR(N3[0], W3[0]);
			vmin[1]=MINVAR(N3[1], W3[1]);
			vmin[2]=MINVAR(N3[2], W3[2]);
			vmax[0]=MAXVAR(N3[0], W3[0]);
			vmax[1]=MAXVAR(N3[1], W3[1]);
			vmax[2]=MAXVAR(N3[2], W3[2]);
			CODE_CG		(3, 0, OCH2_R2);
			CODE_CG		(3, 1, OCH2_G2);
			CODE_CG		(3, 2, OCH2_B2);
			vmin[0]=MINVAR(vmin[0], NE3[0]);
			vmin[1]=MINVAR(vmin[1], NE3[1]);
			vmin[2]=MINVAR(vmin[2], NE3[2]);
			vmax[0]=MAXVAR(vmax[0], NE3[0]);
			vmax[1]=MAXVAR(vmax[1], NE3[1]);
			vmax[2]=MAXVAR(vmax[2], NE3[2]);
			CODE_AV3	(3, 0, OCH2_R2);
			CODE_AV3	(3, 1, OCH2_G2);
			CODE_AV3	(3, 2, OCH2_B2);
			CODE_AV4	(3, 0, OCH2_R2);
			CODE_AV4	(3, 1, OCH2_G2);
			CODE_AV4	(3, 2, OCH2_B2);
			CODE_AV5	(3, 0, OCH2_R2);
			CODE_AV5	(3, 1, OCH2_G2);
			CODE_AV5	(3, 2, OCH2_B2);
			CODE_AV6	(3, 0, OCH2_R2);
			CODE_AV6	(3, 1, OCH2_G2);
			CODE_AV6	(3, 2, OCH2_B2);
			CODE_AV9	(3, 0, OCH2_R2);
			CODE_AV9	(3, 1, OCH2_G2);
			CODE_AV9	(3, 2, OCH2_B2);
			CODE_AVB	(3, 0, OCH2_R2);
			CODE_AVB	(3, 1, OCH2_G2);
			CODE_AVB	(3, 2, OCH2_B2);
			
			CODE_ZERO	(4, 0, OCH2_R3);
			CODE_ZERO	(4, 1, OCH2_G3);
			CODE_ZERO	(4, 2, OCH2_B3);
			CODE_N		(4, 0, OCH2_R3);
			CODE_N		(4, 1, OCH2_G3);
			CODE_N		(4, 2, OCH2_B3);
			CODE_W		(4, 0, OCH2_R3);
			CODE_W		(4, 1, OCH2_G3);
			CODE_W		(4, 2, OCH2_B3);
			CODE_AV2	(4, 0, OCH2_R3);
			CODE_AV2	(4, 1, OCH2_G3);
			CODE_AV2	(4, 2, OCH2_B3);
			CODE_WG		(4, 0, OCH2_R3);
			CODE_WG		(4, 1, OCH2_G3);
			CODE_WG		(4, 2, OCH2_B3);
			vmin[0]=MINVAR(N4[0], W4[0]);
			vmin[1]=MINVAR(N4[1], W4[1]);
			vmin[2]=MINVAR(N4[2], W4[2]);
			vmax[0]=MAXVAR(N4[0], W4[0]);
			vmax[1]=MAXVAR(N4[1], W4[1]);
			vmax[2]=MAXVAR(N4[2], W4[2]);
			CODE_CG		(4, 0, OCH2_R3);
			CODE_CG		(4, 1, OCH2_G3);
			CODE_CG		(4, 2, OCH2_B3);
			vmin[0]=MINVAR(vmin[0], NE4[0]);
			vmin[1]=MINVAR(vmin[1], NE4[1]);
			vmin[2]=MINVAR(vmin[2], NE4[2]);
			vmax[0]=MAXVAR(vmax[0], NE4[0]);
			vmax[1]=MAXVAR(vmax[1], NE4[1]);
			vmax[2]=MAXVAR(vmax[2], NE4[2]);
			CODE_AV3	(4, 0, OCH2_R3);
			CODE_AV3	(4, 1, OCH2_G3);
			CODE_AV3	(4, 2, OCH2_B3);
			CODE_AV4	(4, 0, OCH2_R3);
			CODE_AV4	(4, 1, OCH2_G3);
			CODE_AV4	(4, 2, OCH2_B3);
			CODE_AV5	(4, 0, OCH2_R3);
			CODE_AV5	(4, 1, OCH2_G3);
			CODE_AV5	(4, 2, OCH2_B3);
			CODE_AV6	(4, 0, OCH2_R3);
			CODE_AV6	(4, 1, OCH2_G3);
			CODE_AV6	(4, 2, OCH2_B3);
			CODE_AV9	(4, 0, OCH2_R3);
			CODE_AV9	(4, 1, OCH2_G3);
			CODE_AV9	(4, 2, OCH2_B3);
			CODE_AVB	(4, 0, OCH2_R3);
			CODE_AVB	(4, 1, OCH2_G3);
			CODE_AVB	(4, 2, OCH2_B3);
			
			CODE_ZERO	(5, 0, OCH2_RB);
			CODE_ZERO	(5, 1, OCH2_GR);
			CODE_ZERO	(5, 2, OCH2_BG);
			CODE_N		(5, 0, OCH2_RB);
			CODE_N		(5, 1, OCH2_GR);
			CODE_N		(5, 2, OCH2_BG);
			CODE_W		(5, 0, OCH2_RB);
			CODE_W		(5, 1, OCH2_GR);
			CODE_W		(5, 2, OCH2_BG);
			CODE_AV2	(5, 0, OCH2_RB);
			CODE_AV2	(5, 1, OCH2_GR);
			CODE_AV2	(5, 2, OCH2_BG);
			CODE_WG		(5, 0, OCH2_RB);
			CODE_WG		(5, 1, OCH2_GR);
			CODE_WG		(5, 2, OCH2_BG);
			vmin[0]=MINVAR(N5[0], W5[0]);
			vmin[1]=MINVAR(N5[1], W5[1]);
			vmin[2]=MINVAR(N5[2], W5[2]);
			vmax[0]=MAXVAR(N5[0], W5[0]);
			vmax[1]=MAXVAR(N5[1], W5[1]);
			vmax[2]=MAXVAR(N5[2], W5[2]);
			CODE_CG		(5, 0, OCH2_RB);
			CODE_CG		(5, 1, OCH2_GR);
			CODE_CG		(5, 2, OCH2_BG);
			vmin[0]=MINVAR(vmin[0], NE5[0]);
			vmin[1]=MINVAR(vmin[1], NE5[1]);
			vmin[2]=MINVAR(vmin[2], NE5[2]);
			vmax[0]=MAXVAR(vmax[0], NE5[0]);
			vmax[1]=MAXVAR(vmax[1], NE5[1]);
			vmax[2]=MAXVAR(vmax[2], NE5[2]);
			CODE_AV3	(5, 0, OCH2_RB);
			CODE_AV3	(5, 1, OCH2_GR);
			CODE_AV3	(5, 2, OCH2_BG);
			CODE_AV4	(5, 0, OCH2_RB);
			CODE_AV4	(5, 1, OCH2_GR);
			CODE_AV4	(5, 2, OCH2_BG);
			CODE_AV5	(5, 0, OCH2_RB);
			CODE_AV5	(5, 1, OCH2_GR);
			CODE_AV5	(5, 2, OCH2_BG);
			CODE_AV6	(5, 0, OCH2_RB);
			CODE_AV6	(5, 1, OCH2_GR);
			CODE_AV6	(5, 2, OCH2_BG);
			CODE_AV9	(5, 0, OCH2_RB);
			CODE_AV9	(5, 1, OCH2_GR);
			CODE_AV9	(5, 2, OCH2_BG);
			CODE_AVB	(5, 0, OCH2_RB);
			CODE_AVB	(5, 1, OCH2_GR);
			CODE_AVB	(5, 2, OCH2_BG);
		}
	}
	if(!count)
		return;
	norm=1./count;
	for(kc=0;kc<OCH2_COUNT*PRED_COUNT;++kc)
	{
		int *curr_hist=hist+((size_t)kc<<maxdepth);
		double e=0;
		for(ks=0;ks<maxlevels;++ks)
		{
			int freq=curr_hist[ks];
			if(freq)
			{
				double p=freq*norm;
				e-=p*log(p);
			}
		}
		info->esizes[kc]=e*(1./(8*M_LN2));
	}
	for(kc=0;kc<OCH2_COUNT;++kc)/*select best predictors*/
	{
		int bestpred=0;
		const double *curr_esizes=info->esizes+kc*PRED_COUNT;
		for(kp=1;kp<PRED_COUNT;++kp)
		{
			if(curr_esizes[bestpred]>curr_esizes[kp])
				bestpred=kp;
		}
		predsel[kc]=bestpred;
	}
	bestsize=0;
	for(kt=0;kt<RCT2_COUNT;++kt)/*select best RCT*/
	{
		group=rct_combinations[kt];
		csize=
			+info->esizes[group[0]*PRED_COUNT+predsel[group[0]]]
			+info->esizes[group[1]*PRED_COUNT+predsel[group[1]]]
			+info->esizes[group[2]*PRED_COUNT+predsel[group[2]]]
		;
		info->rctsizes[kt]=csize/3;
		if(!kt||bestsize>csize)
		{
			bestsize=csize;
			info->bestrct=kt;
		}
	}
	group=rct_combinations[info->bestrct];
	info->predidx[0]=predsel[group[0]];
	info->predidx[1]=predsel[group[1]];
	info->predidx[2]=predsel[group[2]];
	free(hist);
	info->t_analysis=time_sec()-t_start;
}