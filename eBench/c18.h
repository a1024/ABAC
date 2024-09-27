#pragma once
#ifndef INC_C18_H
#define INC_C18_H


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

#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	3,	3, 3, 0)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	3,	1, 3, 0)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	3,	0, 3, 0)\
	RCT(R_G_B1,	OCH_R,		OCH_G,		OCH_B1,		0, 1, 2,	3,	0, 1, 1)\
	RCT(R_G_B2,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	3,	0, 1, 1)\
	RCT(R_G_B3,	OCH_R,		OCH_G,		OCH_B3,		0, 1, 2,	3,	0, 1, 1)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	0,	0, 3, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	0,	1, 3, 0)\
	RCT(R_GR_B1,	OCH_R,		OCH_GR,		OCH_B1,		0, 1, 2,	0,	0, 1, 1)\
	RCT(R_GR_B2,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	0,	0, 1, 1)\
	RCT(R_GR_B3,	OCH_R,		OCH_GR,		OCH_B3,		0, 1, 2,	0,	0, 1, 1)\
	RCT(R_B_G1,	OCH_R,		OCH_B,		OCH_G1,		0, 2, 1,	3,	0, 1, 1)\
	RCT(R_B_G2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	3,	0, 1, 1)\
	RCT(R_B_G3,	OCH_R,		OCH_B,		OCH_G3,		0, 2, 1,	3,	0, 1, 1)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	0,	1, 3, 0)\
	RCT(R_BR_G1,	OCH_R,		OCH_BR,		OCH_G1,		0, 2, 1,	0,	0, 1, 1)\
	RCT(R_BR_G2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	0,	0, 1, 1)\
	RCT(R_BR_G3,	OCH_R,		OCH_BR,		OCH_G3,		0, 2, 1,	0,	0, 1, 1)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	3,	0, 3, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	3,	1, 3, 0)\
	RCT(G_B_R1,	OCH_G,		OCH_B,		OCH_R1,		1, 2, 0,	3,	0, 1, 1)\
	RCT(G_B_R2,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	3,	0, 1, 1)\
	RCT(G_B_R3,	OCH_G,		OCH_B,		OCH_R3,		1, 2, 0,	3,	0, 1, 1)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	0,	0, 3, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	0,	1, 3, 0)\
	RCT(G_BG_R1,	OCH_G,		OCH_BG,		OCH_R1,		1, 2, 0,	0,	0, 1, 1)\
	RCT(G_BG_R2,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	0,	0, 1, 1)\
	RCT(G_BG_R3,	OCH_G,		OCH_BG,		OCH_R3,		1, 2, 0,	0,	0, 1, 1)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	0,	1, 3, 0)\
	RCT(G_RG_B1,	OCH_G,		OCH_RG,		OCH_B1,		1, 0, 2,	0,	0, 1, 1)\
	RCT(G_RG_B2,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	0,	0, 1, 1)\
	RCT(G_RG_B3,	OCH_G,		OCH_RG,		OCH_B3,		1, 0, 2,	0,	0, 1, 1)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	3,	1, 3, 0)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	3,	0, 3, 0)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	0,	0, 3, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	0,	1, 3, 0)\
	RCT(B_RB_G1,	OCH_B,		OCH_RB,		OCH_G1,		2, 0, 1,	0,	0, 1, 1)\
	RCT(B_RB_G2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	0,	0, 1, 1)\
	RCT(B_RB_G3,	OCH_B,		OCH_RB,		OCH_G3,		2, 0, 1,	0,	0, 1, 1)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	0,	1, 3, 0)\
	RCT(B_GB_R1,	OCH_B,		OCH_GB,		OCH_R1,		2, 1, 0,	0,	0, 1, 1)\
	RCT(B_GB_R2,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	0,	0, 1, 1)\
	RCT(B_GB_R3,	OCH_B,		OCH_GB,		OCH_R3,		2, 1, 0,	0,	0, 1, 1)
typedef enum _RCTIndex
{
#define RCT(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][10]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2)\
	{YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2) #LABEL,
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
