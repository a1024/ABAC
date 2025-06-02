#pragma once
#ifndef INC_C18_H
#define INC_C18_H


#define OCH2LIST\
	OCH2(R)\
	OCH2(G)\
	OCH2(B)\
	OCH2(RG)\
	OCH2(GB)\
	OCH2(BR)\
	OCH2(R1)\
	OCH2(G1)\
	OCH2(B1)\
	OCH2(R2)\
	OCH2(G2)\
	OCH2(B2)\
	OCH2(R3)\
	OCH2(G3)\
	OCH2(B3)\
	OCH2(RB)\
	OCH2(GR)\
	OCH2(BG)
typedef enum _OCH2Index
{
#define OCH2(LABEL) OCH2_##LABEL,
	OCH2LIST
#undef  OCH2
	OCH2_COUNT,
} OCH2Index;
static const char *och_names[OCH_COUNT]=
{
#define OCH(LABEL) #LABEL,
	OCHLIST
#undef  OCH
};

#define RCT2LIST\
	RCT2(R_G_B,	OCH2_R,		OCH2_G,		OCH2_B,		0, 1, 2,	3,	3, 3, 0)\
	RCT2(R_G_BG,	OCH2_R,		OCH2_G,		OCH2_BG,	0, 1, 2,	3,	1, 3, 0)\
	RCT2(R_G_BR,	OCH2_R,		OCH2_G,		OCH2_BR,	0, 1, 2,	3,	0, 3, 0)\
	RCT2(R_G_B1,	OCH2_R,		OCH2_G,		OCH2_B1,	0, 1, 2,	3,	0, 1, 1)\
	RCT2(R_G_B2,	OCH2_R,		OCH2_G,		OCH2_B2,	0, 1, 2,	3,	0, 1, 1)\
	RCT2(R_G_B3,	OCH2_R,		OCH2_G,		OCH2_B3,	0, 1, 2,	3,	0, 1, 1)\
	RCT2(R_GR_BR,	OCH2_R,		OCH2_GR,	OCH2_BR,	0, 1, 2,	0,	0, 3, 0)\
	RCT2(R_GR_BG,	OCH2_R,		OCH2_GR,	OCH2_BG,	0, 1, 2,	0,	1, 3, 0)\
	RCT2(R_GR_B1,	OCH2_R,		OCH2_GR,	OCH2_B1,	0, 1, 2,	0,	0, 1, 1)\
	RCT2(R_GR_B2,	OCH2_R,		OCH2_GR,	OCH2_B2,	0, 1, 2,	0,	0, 1, 1)\
	RCT2(R_GR_B3,	OCH2_R,		OCH2_GR,	OCH2_B3,	0, 1, 2,	0,	0, 1, 1)\
	RCT2(R_B_G1,	OCH2_R,		OCH2_B,		OCH2_G1,	0, 2, 1,	3,	0, 1, 1)\
	RCT2(R_B_G2,	OCH2_R,		OCH2_B,		OCH2_G2,	0, 2, 1,	3,	0, 1, 1)\
	RCT2(R_B_G3,	OCH2_R,		OCH2_B,		OCH2_G3,	0, 2, 1,	3,	0, 1, 1)\
	RCT2(R_BR_GB,	OCH2_R,		OCH2_BR,	OCH2_GB,	0, 2, 1,	0,	1, 3, 0)\
	RCT2(R_BR_G1,	OCH2_R,		OCH2_BR,	OCH2_G1,	0, 2, 1,	0,	0, 1, 1)\
	RCT2(R_BR_G2,	OCH2_R,		OCH2_BR,	OCH2_G2,	0, 2, 1,	0,	0, 1, 1)\
	RCT2(R_BR_G3,	OCH2_R,		OCH2_BR,	OCH2_G3,	0, 2, 1,	0,	0, 1, 1)\
	RCT2(G_B_RG,	OCH2_G,		OCH2_B,		OCH2_RG,	1, 2, 0,	3,	0, 3, 0)\
	RCT2(G_B_RB,	OCH2_G,		OCH2_B,		OCH2_RB,	1, 2, 0,	3,	1, 3, 0)\
	RCT2(G_B_R1,	OCH2_G,		OCH2_B,		OCH2_R1,	1, 2, 0,	3,	0, 1, 1)\
	RCT2(G_B_R2,	OCH2_G,		OCH2_B,		OCH2_R2,	1, 2, 0,	3,	0, 1, 1)\
	RCT2(G_B_R3,	OCH2_G,		OCH2_B,		OCH2_R3,	1, 2, 0,	3,	0, 1, 1)\
	RCT2(G_BG_RG,	OCH2_G,		OCH2_BG,	OCH2_RG,	1, 2, 0,	0,	0, 3, 0)\
	RCT2(G_BG_RB,	OCH2_G,		OCH2_BG,	OCH2_RB,	1, 2, 0,	0,	1, 3, 0)\
	RCT2(G_BG_R1,	OCH2_G,		OCH2_BG,	OCH2_R1,	1, 2, 0,	0,	0, 1, 1)\
	RCT2(G_BG_R2,	OCH2_G,		OCH2_BG,	OCH2_R2,	1, 2, 0,	0,	0, 1, 1)\
	RCT2(G_BG_R3,	OCH2_G,		OCH2_BG,	OCH2_R3,	1, 2, 0,	0,	0, 1, 1)\
	RCT2(G_RG_BR,	OCH2_G,		OCH2_RG,	OCH2_BR,	1, 0, 2,	0,	1, 3, 0)\
	RCT2(G_RG_B1,	OCH2_G,		OCH2_RG,	OCH2_B1,	1, 0, 2,	0,	0, 1, 1)\
	RCT2(G_RG_B2,	OCH2_G,		OCH2_RG,	OCH2_B2,	1, 0, 2,	0,	0, 1, 1)\
	RCT2(G_RG_B3,	OCH2_G,		OCH2_RG,	OCH2_B3,	1, 0, 2,	0,	0, 1, 1)\
	RCT2(B_R_GR,	OCH2_B,		OCH2_R,		OCH2_GR,	2, 0, 1,	3,	1, 3, 0)\
	RCT2(B_R_GB,	OCH2_B,		OCH2_R,		OCH2_GB,	2, 0, 1,	3,	0, 3, 0)\
	RCT2(B_RB_GB,	OCH2_B,		OCH2_RB,	OCH2_GB,	2, 0, 1,	0,	0, 3, 0)\
	RCT2(B_RB_GR,	OCH2_B,		OCH2_RB,	OCH2_GR,	2, 0, 1,	0,	1, 3, 0)\
	RCT2(B_RB_G1,	OCH2_B,		OCH2_RB,	OCH2_G1,	2, 0, 1,	0,	0, 1, 1)\
	RCT2(B_RB_G2,	OCH2_B,		OCH2_RB,	OCH2_G2,	2, 0, 1,	0,	0, 1, 1)\
	RCT2(B_RB_G3,	OCH2_B,		OCH2_RB,	OCH2_G3,	2, 0, 1,	0,	0, 1, 1)\
	RCT2(B_GB_RG,	OCH2_B,		OCH2_GB,	OCH2_RG,	2, 1, 0,	0,	1, 3, 0)\
	RCT2(B_GB_R1,	OCH2_B,		OCH2_GB,	OCH2_R1,	2, 1, 0,	0,	0, 1, 1)\
	RCT2(B_GB_R2,	OCH2_B,		OCH2_GB,	OCH2_R2,	2, 1, 0,	0,	0, 1, 1)\
	RCT2(B_GB_R3,	OCH2_B,		OCH2_GB,	OCH2_R3,	2, 1, 0,	0,	0, 1, 1)
typedef enum _RCT2Index
{
#define RCT2(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2) RCT2_##LABEL,
	RCT2LIST
#undef  RCT2
	RCT2_COUNT,
} RCT2Index;
static const unsigned char rct2_combinations[RCT2_COUNT][10]=
{
#define RCT2(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2)\
	{YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2},
	RCT2LIST
#undef  RCT2
};
static const char *rct2_names[RCT2_COUNT]=
{
#define RCT2(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2) #LABEL,
	RCT2LIST
#undef  RCT2
};

#define PREDLIST\
	PRED(ZERO, "ZERO")\
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
static const short av12_icoeffs[12]=
{
	 0x04,	 0x03,	-0x1F,	-0x26,	 0x00,
	 0x07,	-0x9E,	 0xDB,	 0x1E,	 0x13,
	-0x2A,	 0xF3,
};
typedef struct _C18Info
{
	double esizes[OCH_COUNT*PRED_COUNT], rctsizes[RCT_COUNT];
	int bestrct, predidx[3];
//	double ACsize[3], GRsize[3];
	double t_analysis;
} C18Info;

void c18_analyze(Image const *src, int x1, int x2, int y1, int y2, C18Info *info);

#endif//INC_C18_H
