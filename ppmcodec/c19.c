#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define ESTIMATE_SIZE


#define AC3_PREC
#include"entropy.h"


#define OCHLIST\
	OCH(R)\
	OCH(G)\
	OCH(B)\
	OCH(RG)\
	OCH(GB)\
	OCH(BR)\
	OCH(R1)\
	OCH(G1)\
	OCH(B1)\
	OCH(R2)\
	OCH(G2)\
	OCH(B2)\
	OCH(R3)\
	OCH(G3)\
	OCH(B3)\
	OCH(RB)\
	OCH(GR)\
	OCH(BG)
typedef enum _OCHIndex
{
#define OCH(LABEL) OCH_##LABEL,
	OCHLIST
#undef  OCH
	OCH_COUNT,
} OCHIndex;
static const char *och_names[OCH_COUNT]=
{
#define OCH(LABEL) #LABEL,
	OCHLIST
#undef  OCH
};

//rct_name, output_channels, permutation, 2*3 coeffs
typedef enum _RCTInfo
{
	IDX_OCH_Y,
	IDX_OCH_U,
	IDX_OCH_V,
	IDX_PERM_Y,
	IDX_PERM_U,
	IDX_PERM_V,
	IDX_YC0,
	IDX_UC0,
	IDX_VC0,
	IDX_YC1,
	IDX_UC1,
	IDX_VC1,

	IDX_COUNT,
} RCTInfo;
#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0, 0,	0, 0,	0, 0)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0, 0,	0, 0,	0, 4)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0, 0,	0, 0,	4, 0)\
	RCT(R_G_B1,	OCH_R,		OCH_G,		OCH_B1,		0, 1, 2,	0, 0,	0, 0,	3, 1)\
	RCT(R_G_B2,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0, 0,	0, 0,	2, 2)\
	RCT(R_G_B3,	OCH_R,		OCH_G,		OCH_B3,		0, 1, 2,	0, 0,	0, 0,	1, 3)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	0, 0,	4, 0,	4, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	0, 0,	4, 0,	0, 4)\
	RCT(R_GR_B1,	OCH_R,		OCH_GR,		OCH_B1,		0, 1, 2,	0, 0,	4, 0,	3, 1)\
	RCT(R_GR_B2,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	0, 0,	4, 0,	2, 2)\
	RCT(R_GR_B3,	OCH_R,		OCH_GR,		OCH_B3,		0, 1, 2,	0, 0,	4, 0,	1, 3)\
	RCT(R_B_G1,	OCH_R,		OCH_B,		OCH_G1,		0, 2, 1,	0, 0,	0, 0,	1, 3)\
	RCT(R_B_G2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0, 0,	0, 0,	2, 2)\
	RCT(R_B_G3,	OCH_R,		OCH_B,		OCH_G3,		0, 2, 1,	0, 0,	0, 0,	3, 1)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	0, 0,	4, 0,	0, 4)\
	RCT(R_BR_G1,	OCH_R,		OCH_BR,		OCH_G1,		0, 2, 1,	0, 0,	4, 0,	1, 3)\
	RCT(R_BR_G2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	0, 0,	4, 0,	2, 2)\
	RCT(R_BR_G3,	OCH_R,		OCH_BR,		OCH_G3,		0, 2, 1,	0, 0,	4, 0,	3, 1)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0, 0,	0, 0,	4, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0, 0,	0, 0,	0, 4)\
	RCT(G_B_R1,	OCH_G,		OCH_B,		OCH_R1,		1, 2, 0,	0, 0,	0, 0,	3, 1)\
	RCT(G_B_R2,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0, 0,	0, 0,	2, 2)\
	RCT(G_B_R3,	OCH_G,		OCH_B,		OCH_R3,		1, 2, 0,	0, 0,	0, 0,	1, 3)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	0, 0,	4, 0,	4, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	0, 0,	4, 0,	0, 4)\
	RCT(G_BG_R1,	OCH_G,		OCH_BG,		OCH_R1,		1, 2, 0,	0, 0,	4, 0,	3, 1)\
	RCT(G_BG_R2,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	0, 0,	4, 0,	2, 2)\
	RCT(G_BG_R3,	OCH_G,		OCH_BG,		OCH_R3,		1, 2, 0,	0, 0,	4, 0,	1, 3)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	0, 0,	4, 0,	0, 4)\
	RCT(G_RG_B1,	OCH_G,		OCH_RG,		OCH_B1,		1, 0, 2,	0, 0,	4, 0,	1, 3)\
	RCT(G_RG_B2,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	0, 0,	4, 0,	2, 2)\
	RCT(G_RG_B3,	OCH_G,		OCH_RG,		OCH_B3,		1, 0, 2,	0, 0,	4, 0,	3, 1)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0, 0,	0, 0,	0, 4)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0, 0,	0, 0,	4, 0)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	0, 0,	4, 0,	4, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	0, 0,	4, 0,	0, 4)\
	RCT(B_RB_G1,	OCH_B,		OCH_RB,		OCH_G1,		2, 0, 1,	0, 0,	4, 0,	1, 3)\
	RCT(B_RB_G2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	0, 0,	4, 0,	2, 2)\
	RCT(B_RB_G3,	OCH_B,		OCH_RB,		OCH_G3,		2, 0, 1,	0, 0,	4, 0,	3, 1)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	0, 0,	4, 0,	0, 4)\
	RCT(B_GB_R1,	OCH_B,		OCH_GB,		OCH_R1,		2, 1, 0,	0, 0,	4, 0,	1, 3)\
	RCT(B_GB_R2,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	0, 0,	4, 0,	2, 2)\
	RCT(B_GB_R3,	OCH_B,		OCH_GB,		OCH_R3,		2, 1, 0,	0, 0,	4, 0,	3, 1)
typedef enum _RCTIndex
{
#define RCT(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  YC0, YC1,  UC0, UC1,  VC0, VC1) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][IDX_COUNT]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  YC0, YC1,  UC0, UC1,  VC0, VC1)\
		  {YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  YC0, UC0, VC0,  YC1, UC1, VC1},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  YC0, YC1,  UC0, UC1,  VC0, VC1) #LABEL,
	RCTLIST
#undef  RCT
};

#define PREDLIST\
	PRED(N, "N")\
	PRED(W, "W")\
	PRED(AV2, "(N+W)/2")\
	PRED(WG, "(gx*N+gy*W)/(gx+gy)")\
	PRED(CG, "median(N, W, N+W-NW)")\
	PRED(AV3, "[-2 [3];3 ?]/4")\
	PRED(AV4, "[-1 [4] 1;4 ?]/8")\
	PRED(AV5, "[-5 [5] 1;-1 8 ?]/8")\
	PRED(AV6, "[[-1];-5 [6] 1;-1 8 ?]/8")\
	PRED(AV9, "[1 [-2] -1;-1 -9 [10] 4;-2 16 ?]/16")\
	PRED(AVB, "[4 3 [-31] -38;7 -158 [219] 30 19;-42 243 ?]/256")
static const char *pred_desc[]=
{
#define PRED(LABEL, DESC) DESC,
	PREDLIST
#undef  PRED
};
typedef enum _PredIndex
{
#define PRED(LABEL, DESC) PRED_##LABEL,
	PREDLIST
#undef  PRED
	PRED_COUNT,
} PredIndex;
static const char *pred_names[PRED_COUNT]=
{
#define PRED(LABEL, DESC) #LABEL,
	PREDLIST
#undef  PRED
};
typedef enum _NBIndex
{
	NB_NNWW,	NB_NNW,		NB_NN,		NB_NNE,		NB_NNEE,
	NB_NWW,		NB_NW,		NB_N,		NB_NE,		NB_NEE,
	NB_WW,		NB_W,		NB_curr,

	NB_COUNT,
} NBIndex;
//static const short av12_icoeffs[12]=
//{
//	 0x04,	 0x03,	-0x1F,	-0x26,	 0x00,
//	 0x07,	-0x9E,	 0xDB,	 0x1E,	 0x13,
//	-0x2A,	 0xF3,
//};
#define CODE_N(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_N	=N	##NBIDX[KC];\
		int pred=nb_N;\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -128, 127);\
		pred=(target[KC]-pred+128)&255;\
		++hist[(OCHIDX*PRED_COUNT+PRED_N)<<8|pred];\
	}while(0)
#define CODE_W(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_W	=W	##NBIDX[KC];\
		int pred=nb_W;\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -128, 127);\
		pred=(target[KC]-pred+128)&255;\
		++hist[(OCHIDX*PRED_COUNT+PRED_W)<<8|pred];\
	}while(0)
#define CODE_AV2(NBIDX, KC, OCHIDX)\
	do\
	{\
		int\
			nb_N	=N	##NBIDX[KC],\
			nb_W	=W	##NBIDX[KC];\
		int pred=(nb_N+nb_W)/2;\
		pred+=offset##NBIDX[KC];\
		CLAMP2(pred, -128, 127);\
		pred=(target[KC]-pred+128)&255;\
		++hist[(OCHIDX*PRED_COUNT+PRED_AV2)<<8|pred];\
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
		CLAMP2(pred, -128, 127);\
		pred=(target[KC]-pred+128)&255;\
		++hist[(OCHIDX*PRED_COUNT+PRED_WG)<<8|pred];\
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
		CLAMP2(pred, -128, 127);\
		pred=(target[KC]-pred+128)&255;\
		++hist[(OCHIDX*PRED_COUNT+PRED_CG)<<8|pred];\
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
		CLAMP2(pred, -128, 127);\
		pred=(target[KC]-pred+128)&255;\
		++hist[(OCHIDX*PRED_COUNT+PRED_AV3)<<8|pred];\
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
		CLAMP2(pred, -128, 127);\
		pred=(target[KC]-pred+128)&255;\
		++hist[(OCHIDX*PRED_COUNT+PRED_AV4)<<8|pred];\
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
		CLAMP2(pred, -128, 127);\
		pred=(target[KC]-pred+128)&255;\
		++hist[(OCHIDX*PRED_COUNT+PRED_AV5)<<8|pred];\
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
		CLAMP2(pred, -128, 127);\
		pred=(target[KC]-pred+128)&255;\
		++hist[(OCHIDX*PRED_COUNT+PRED_AV6)<<8|pred];\
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
		CLAMP2(pred, -128, 127);\
		pred=(target[KC]-pred+128)&255;\
		++hist[(OCHIDX*PRED_COUNT+PRED_AV9)<<8|pred];\
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
		CLAMP2(pred, -128, 127);\
		pred=(target[KC]-pred+128)&255;\
		++hist[(OCHIDX*PRED_COUNT+PRED_AVB)<<8|pred];\
	}while(0)
typedef struct _AnalysisInfo
{
	int bestrct, predidx[3];
	int hist[3][256];
} AnalysisInfo;
static void c19_analyze(const unsigned char *src, int iw, int ih, AnalysisInfo *info)
{
#ifdef _MSC_VER
	double t_start=time_sec();
#endif
	int histsize=sizeof(int[OCH_COUNT*PRED_COUNT<<8]);
	int *hist=(int*)malloc(histsize);
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(hist, 0, histsize);
	int count=0;
	int vmin[3]={0}, vmax[3]={0};
	for(int ky=0+2;ky<ih;++ky)//analysis loop
	{
		const unsigned char
			*NNptr	=src+3LL*(iw*(ky-2LL)+2),
			*Nptr	=src+3LL*(iw*(ky-1LL)+2),
			*currptr=src+3LL*(iw*(ky+0LL)+2);
		for(int kx=0+2;kx<iw-2;++kx, ++count)
		{
#define DECL_NN(X)   {NNptr	[3*X+0]-128, NNptr	[3*X+1]-128, NNptr	[3*X+2]-128}
#define DECL_N(X)    {Nptr	[3*X+0]-128, Nptr	[3*X+1]-128, Nptr	[3*X+2]-128}
#define DECL_curr(X) {currptr	[3*X+0]-128, currptr	[3*X+1]-128, currptr	[3*X+2]-128}
			int
				NNWW0	[]=DECL_NN(-2),
				NNW0	[]=DECL_NN(-1),
				NN0	[]=DECL_NN( 0),
				NNE0	[]=DECL_NN( 1),
			//	NNEE0	[]=DECL_NN( 2),
				NWW0	[]=DECL_N(-2),
				NW0	[]=DECL_N(-1),
				N0	[]=DECL_N( 0),
				NE0	[]=DECL_N( 1),
				NEE0	[]=DECL_N( 2),
				WW0	[]=DECL_curr(-2),
				W0	[]=DECL_curr(-1),
				target	[]=DECL_curr( 0),
				offset0	[]={0, 0, 0};
#undef  DECL_NN
#undef  DECL_N
#undef  DECL_curr
			//const unsigned char
			//	*NNWW0	=src+3*(iw*(ky-2LL)+kx-2LL),
			//	*NNW0	=src+3*(iw*(ky-2LL)+kx-1LL),
			//	*NN0	=src+3*(iw*(ky-2LL)+kx+0LL),
			//	*NNE0	=src+3*(iw*(ky-2LL)+kx+1LL),
			//	*NNEE0	=src+3*(iw*(ky-2LL)+kx+2LL),
			//	*NWW0	=src+3*(iw*(ky-1LL)+kx-2LL),
			//	*NW0	=src+3*(iw*(ky-1LL)+kx-1LL),
			//	*N0	=src+3*(iw*(ky-1LL)+kx+0LL),
			//	*NE0	=src+3*(iw*(ky-1LL)+kx+1LL),
			//	*NEE0	=src+3*(iw*(ky-1LL)+kx+2LL),
			//	*WW0	=src+3*(iw*(ky+0LL)+kx-2LL),
			//	*W0	=src+3*(iw*(ky+0LL)+kx-1LL),
			//	*target	=src+3*(iw*(ky+0LL)+kx+0LL),
			//	offset0	[]={0, 0, 0};
#define DECL_NB(NB) {NB##0[0]-NB##0[1], NB##0[1]-NB##0[2], NB##0[2]-NB##0[0]}
			int
				NNWW1	[]=DECL_NB(NNWW),
				NNW1	[]=DECL_NB(NNW),
				NN1	[]=DECL_NB(NN),
				NNE1	[]=DECL_NB(NNE),
			//	NNEE1	[]=DECL_NB(NNEE),
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
			//	NNEE2	[]=DECL_NB(NNEE),
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
			//	NNEE3	[]=DECL_NB(NNEE),
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
			//	NNEE4	[]=DECL_NB(NNEE),
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
			//	NNEE5	[]=DECL_NB(NNEE),
				NWW5	[]=DECL_NB(NWW),
				NW5	[]=DECL_NB(NW),
				N5	[]=DECL_NB(N),
				NE5	[]=DECL_NB(NE),
				NEE5	[]=DECL_NB(NEE),
				WW5	[]=DECL_NB(WW),
				W5	[]=DECL_NB(W),
				offset5	[]={target[2], target[0], target[1]};
#undef  DECL_NB
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

			NNptr+=3;
			Nptr+=3;
			currptr+=3;
		}
	}
	if(!count)
		return;
	double norm=1./count;
	double esizes[OCH_COUNT*PRED_COUNT]={0};
	for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)
	{
		int *curr_hist=hist+((size_t)kc<<8);
		double e=0;
		for(int ks=0;ks<256;++ks)
		{
			int freq=curr_hist[ks];
			if(freq)
			{
				double p=freq*norm;
				e-=p*log2(p);
			}
		}
		esizes[kc]=e/8;
#ifdef _MSC_VER
		printf("%-4s %-4s %10.2lf %8.4lf%%\n",
			och_names[kc/PRED_COUNT],
			pred_names[kc%PRED_COUNT],
			esizes[kc]*count,
			esizes[kc]*100
		);
#endif
	}
	int predsel[OCH_COUNT]={0};
	for(int kc=0;kc<OCH_COUNT;++kc)//select best predictors
	{
		int bestpred=0;
		for(int kp=1;kp<PRED_COUNT;++kp)
		{
			if(esizes[kc*PRED_COUNT+bestpred]>esizes[kc*PRED_COUNT+kp])
				bestpred=kp;
		}
		predsel[kc]=bestpred;
	}
	double bestsize=0;
	int bestrct=0;
	for(int kt=0;kt<RCT_COUNT;++kt)//select best RCT
	{
		const unsigned char *group=rct_combinations[kt];
		double csize=
			+esizes[group[0]*PRED_COUNT+predsel[group[0]]]
			+esizes[group[1]*PRED_COUNT+predsel[group[1]]]
			+esizes[group[2]*PRED_COUNT+predsel[group[2]]];
		if(!kt||bestsize>csize)
		{
			bestsize=csize;
			bestrct=kt;
		}
#ifdef _MSC_VER
		printf("%s\t%s %s %s\t%10.2lf %8.4lf%%\n",
			rct_names[kt],
			pred_names[predsel[group[0]]],
			pred_names[predsel[group[1]]],
			pred_names[predsel[group[2]]],
			csize*count,
			csize*100/3
		);
#endif
	}
	const unsigned char *group=rct_combinations[bestrct];
	info->bestrct=bestrct;
	info->predidx[0]=predsel[group[0]];
	info->predidx[1]=predsel[group[1]];
	info->predidx[2]=predsel[group[2]];
	memcpy(info->hist[0], hist+(((size_t)group[0]*PRED_COUNT+predsel[group[0]])<<8), sizeof(int[256]));
	memcpy(info->hist[1], hist+(((size_t)group[1]*PRED_COUNT+predsel[group[1]])<<8), sizeof(int[256]));
	memcpy(info->hist[2], hist+(((size_t)group[2]*PRED_COUNT+predsel[group[2]])<<8), sizeof(int[256]));
	free(hist);
#ifdef _MSC_VER
	printf("\n");
	printf("%s %s %s %s  %lf sec\n",
		rct_names[bestrct],
		pred_names[info->predidx[0]],
		pred_names[info->predidx[1]],
		pred_names[info->predidx[2]],
		time_sec()-t_start
	);
#endif
}


static void update_CDF(const int *hist, unsigned *CDF)
{
	int hsum=0, nused=0;
	for(int ks=0;ks<256;++ks)
	{
		int freq=hist[ks];
		hsum+=freq;
		nused+=freq!=0;
	}
	for(int ks=0, c=0, guard=0;ks<256;++ks)
	{
		int freq=hist[ks];
	//	CDF[ks]=(int)(c*((1LL<<16)-nused)/hsum)+guard;
		CDF[ks]=(int)(c*((1LL<<16)-256)/hsum)+ks;
		c+=freq;
		guard+=freq!=0;
	}
	CDF[256]=1<<16;
	//const int hsum=hist[nlevels];
	//for(int ks=0, c=0;ks<nlevels;++ks)
	//{
	//	int freq=hist[ks];
	//	CDF[ks]=(int)(c*((1LL<<16)-nlevels)/hsum)+ks;
	//	c+=freq;
	//}
	//CDF[nlevels]=1<<16;
}
int c19_codec(const char *srcfn, const char *dstfn)
{
#ifdef ESTIMATE_SIZE
	double csizes[3]={0};
#endif
#ifdef _MSC_VER
	double t_start=time_sec();
	unsigned long long c_start=__rdtsc();
#endif
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	ptrdiff_t srcsize;
	unsigned char *srcbuf=0, *srcptr=0, *srcend=0;
	{
		srcsize=get_filesize(srcfn);
		if(srcsize<1)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		srcbuf=(unsigned char*)malloc(srcsize+1);
		if(!srcbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		fread(srcbuf, 1, srcsize, fsrc);
		fclose(fsrc);
		srcbuf[srcsize]=0;
	}
	srcptr=srcbuf;
	srcend=srcbuf+srcsize;

	if(srcsize<2)
	{
		LOG_ERROR("File is empty");
		return 1;
	}
	int tag=*(unsigned short*)srcptr;
	srcptr+=2;
	int fwd=tag==('P'|'6'<<8);
	if(!fwd&&tag!=('C'|'H'<<8))
	{
		LOG_ERROR("Unsupported source file");
		return 1;
	}
	int iw=0, ih=0;
	ptrdiff_t dstbufsize;
	unsigned char *dstbuf=0, *image=0;
	if(fwd)//parse PPM header
	{
		if(*srcptr++ != '\n')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		while((unsigned)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++ - '0';
		if(iw<1||*srcptr++ != ' ')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		while((unsigned)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++ - '0';
		if(ih<1||*srcptr++ != '\n')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		int temp=0;
		while((unsigned)(*srcptr-'0')<10)
			temp=10*temp+*srcptr++ - '0';
		if(temp!=255||*srcptr++ != '\n')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		image=srcptr;
	}
	else
	{
		if(srcptr+8>=srcend)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		memcpy(&iw, srcptr, 4); srcptr+=4;
		memcpy(&ih, srcptr, 4); srcptr+=4;

		dstbufsize=4LL*iw*ih;
		dstbuf=(unsigned char*)malloc(dstbufsize);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		image=dstbuf;
	}

	int psize=(iw+16LL)*sizeof(short[4*4]);//4 padded rows * 4 channels max
	short *pixels=(short*)_mm_malloc(psize, sizeof(__m128i));
	int hsize=(int)sizeof(unsigned[3][257]);
	unsigned *CDF=(unsigned*)malloc(hsize);
	unsigned char *CDF2sym=0;
	if(!pixels||!CDF)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	BList list={0};
	AC3 ec={0};
	memset(pixels, 0, psize);
	int bestrct=0;
	const unsigned char *rctinfo=0;
	int predidx[3]={0};
	if(fwd)//encode
	{
		AnalysisInfo info={0};
		c19_analyze(image, iw, ih, &info);
		update_CDF(info.hist[0], CDF+257*0);
		update_CDF(info.hist[1], CDF+257*1);
		update_CDF(info.hist[2], CDF+257*2);
		bestrct=info.bestrct;
		rctinfo=rct_combinations[info.bestrct];
		predidx[0]=info.predidx[0];
		predidx[1]=info.predidx[1];
		predidx[2]=info.predidx[2];
		
		blist_init(&list);
		ac3_enc_init(&ec, &list);
		ac3_enc_bypass_NPOT(&ec, bestrct, RCT_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[0], PRED_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[1], PRED_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[2], PRED_COUNT);
		for(int ks=1;ks<256;++ks)
		{
			ac3_enc_bypass_NPOT(&ec, CDF[257*0+ks]-CDF[257*0+ks-1], 0x10000-CDF[257*0+ks-1]);
			ac3_enc_bypass_NPOT(&ec, CDF[257*1+ks]-CDF[257*1+ks-1], 0x10000-CDF[257*1+ks-1]);
			ac3_enc_bypass_NPOT(&ec, CDF[257*2+ks]-CDF[257*2+ks-1], 0x10000-CDF[257*2+ks-1]);
		}
	}
	else
	{
		CDF2sym=(unsigned char*)malloc(sizeof(char[0x30000]));
		if(!CDF2sym)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		ac3_dec_init(&ec, srcptr, srcend);
		bestrct=ac3_dec_bypass_NPOT(&ec, RCT_COUNT);
		rctinfo=rct_combinations[bestrct];
		predidx[0]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[1]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[2]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		CDF[257*0+0]=0;
		CDF[257*1+0]=0;
		CDF[257*2+0]=0;
		for(int ks=1;ks<256;++ks)
		{
			CDF[257*0+ks]=CDF[257*0+ks-1]+ac3_dec_bypass_NPOT(&ec, 0x10000-CDF[257*0+ks-1]);
			CDF[257*1+ks]=CDF[257*1+ks-1]+ac3_dec_bypass_NPOT(&ec, 0x10000-CDF[257*1+ks-1]);
			CDF[257*2+ks]=CDF[257*2+ks-1]+ac3_dec_bypass_NPOT(&ec, 0x10000-CDF[257*2+ks-1]);
		}
		CDF[257*0+256]=0x10000;
		CDF[257*1+256]=0x10000;
		CDF[257*2+256]=0x10000;
		for(int ks=0;ks<256;++ks)
		{
			int cdf, freq;

			cdf=CDF[257*0+ks];
			freq=CDF[257*0+ks+1]-cdf;
			if(freq)
				memset(CDF2sym+0x10000*0+cdf, ks, freq);

			cdf=CDF[257*1+ks];
			freq=CDF[257*1+ks+1]-cdf;
			if(freq)
				memset(CDF2sym+0x10000*1+cdf, ks, freq);

			cdf=CDF[257*2+ks];
			freq=CDF[257*2+ks+1]-cdf;
			if(freq)
				memset(CDF2sym+0x10000*2+cdf, ks, freq);
		}
	}
	for(int ky=0, idx=0;ky<ih;++ky)
	{
		ALIGN(32) short *rows[]=
		{
			pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4,
		};
		int pred=0;
		char yuv[3]={0};
		unsigned char sym=0;
		for(int kx=0;kx<iw;++kx, idx+=3)
		{
			short
				*NNWW	=rows[2]-2*4,
				*NNW	=rows[2]-1*4,
				*NN	=rows[2]+0*4,
				*NNE	=rows[2]+1*4,
			//	*NNEE	=rows[2]+2*4,
				*NWW	=rows[1]-2*4,
				*NW	=rows[1]-1*4,
				*N	=rows[1]+0*4,
				*NE	=rows[1]+1*4,
				*NEE	=rows[1]+2*4,
				*WW	=rows[0]-2*4,
				*W	=rows[0]-1*4,
				*curr	=rows[0]+0*4;
			if(fwd)
			{
				yuv[0]=image[idx+rctinfo[IDX_PERM_Y]]-128;
				yuv[1]=image[idx+rctinfo[IDX_PERM_U]]-128;
				yuv[2]=image[idx+rctinfo[IDX_PERM_V]]-128;
			}
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
			for(int kc=0;kc<3;++kc)
			{
				int offset=(rctinfo[kc+IDX_YC0]*yuv[0]+rctinfo[kc+IDX_YC1]*yuv[1])>>2;
				switch(predidx[kc])
				{
				case PRED_N:
					pred=N[kc];
					break;
				case PRED_W:
					pred=W[kc];
					break;
				case PRED_AV2:
					pred=(N[kc]+W[kc]+1)>>1;
					break;
				case PRED_WG:
					{
						int gx=abs(W[kc]-WW[kc])+abs(N[kc]-NW[kc])+abs(NE[kc]-N[kc])+1;
						int gy=abs(N[kc]-NN[kc])+abs(W[kc]-NW[kc])+abs(NE[kc]-NNE[kc])+1;
						pred=(gx*N[kc]+gy*W[kc])/(gx+gy);
					}
					break;
				case PRED_CG:
					MEDIAN3_32(pred, N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
					break;
				case PRED_AV3:
					CLAMP3_32(pred, (3*(N[kc]+W[kc])-2*NW[kc]+2)>>2, N[kc], W[kc], NE[kc]);
					break;
				case PRED_AV4:
					CLAMP3_32(pred, (4*(N[kc]+W[kc])+NE[kc]-NW[kc]+4)>>3, N[kc], W[kc], NE[kc]);
					break;
				case PRED_AV5:
					CLAMP3_32(pred, W[kc]+((5*(N[kc]-NW[kc])+NE[kc]-WW[kc]+4)>>3), N[kc], W[kc], NE[kc]);
					break;
				case PRED_AV6:
					CLAMP3_32(pred, W[kc]+((6*N[kc]-5*NW[kc]+NE[kc]-NN[kc]-WW[kc]+4)>>3), N[kc], W[kc], NE[kc]);
					break;
				case PRED_AV9:
					pred=W[kc]+((10*N[kc]-9*NW[kc]+4*NE[kc]-2*(NN[kc]+WW[kc])-NNE[kc]+NNW[kc]-NWW[kc]+8)>>4);
					CLAMP3_32(pred, pred, N[kc], W[kc], NE[kc]);
					break;
				case PRED_AVB:
					pred=(
						+0x04*NNWW	[kc]+0x03*NNW	[kc]-0x1F*NN	[kc]-0x26*NNE	[kc]
						+0x07*NWW	[kc]-0x9E*NW	[kc]+0xDB*N	[kc]+0x1E*NE	[kc]+0x13*NEE	[kc]
						-0x2A*WW	[kc]+0xF3*W	[kc]
					+128)>>8;
					CLAMP3_32(pred, pred, N[kc], W[kc], NE[kc]);
					break;
				}
				pred+=offset;
				CLAMP2(pred, -128, 127);

				int cdf, freq;
				if(fwd)
				{
					curr[kc]=yuv[kc];
					sym=(yuv[kc]-pred+128)&255;
					cdf=CDF[257*kc+sym];
					freq=CDF[257*kc+sym+1]-cdf;
					ac3_enc_update(&ec, cdf, freq);
#ifdef ESTIMATE_SIZE
					csizes[kc]-=log2((double)freq/0x10000);
#endif
				}
				else
				{
					int code=ac3_dec_getcdf(&ec);
					sym=CDF2sym[kc<<16|code];
					cdf=CDF[257*kc+sym];
					freq=CDF[257*kc+sym+1]-cdf;
					ac3_dec_update(&ec, cdf, freq);
					yuv[kc]=((sym+pred)&255)-128;
					curr[kc]=yuv[kc];
				}
				curr[kc]-=offset;
			}
			if(!fwd)
			{
				image[idx+rctinfo[IDX_PERM_Y]]=yuv[0]+128;
				image[idx+rctinfo[IDX_PERM_U]]=yuv[1]+128;
				image[idx+rctinfo[IDX_PERM_V]]=yuv[2]+128;
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	ptrdiff_t csize=0;
	FILE *fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
	ptrdiff_t imsize=0;
	if(fwd)
	{
		ArrayHandle stream;
		ac3_enc_flush(&ec);
		csize=list.nbytes+10;
		ARRAY_ALLOC(char, stream, 0, 0, csize, 0);

		fwrite("CH", 1, 2, fdst);
		fwrite(&iw, 1, 4, fdst);
		fwrite(&ih, 1, 4, fdst);
		blist_appendtofile(&list, fdst);
		blist_clear(&list);

		fwrite(stream->data, 1, stream->count, fdst);
		array_free(&stream);
		imsize=srcsize;
	}
	else
	{
		csize=srcsize;
		imsize+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		fwrite(image, 1, 3LL*iw*ih, fdst);
		free(CDF2sym);
		imsize+=3LL*iw*ih;
	}
	fclose(fdst);
	_mm_free(pixels);
	free(CDF);
#ifdef _MSC_VER
	{
		unsigned long long c_elapsed=__rdtsc()-c_start;
		double elapsed=time_sec()-t_start;
#ifdef ESTIMATE_SIZE
		if(fwd)
		{
			printf("Y %10.2lf\n", csizes[0]/8);
			printf("U %10.2lf\n", csizes[1]/8);
			printf("V %10.2lf\n", csizes[2]/8);
		}
#endif
		printf("%td -> %td bytes  %lf:1  %12lf sec  %12lf MB/s  %12lf C/B\n",
			fwd?imsize:csize,
			fwd?csize:imsize,
			(double)imsize/csize,
			elapsed,
			imsize/elapsed/(1024*1024),
			(double)c_elapsed/imsize
		);
	}
#endif
	return 0;
}