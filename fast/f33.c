#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#ifdef _MSC_VER
#include<intrin.h>
#endif
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT

	#define ALLOW_AVX2
	#define ENABLE_OLS
	#define ENABLE_SSE


#define AC3_PREC
#include"ac.h"
#define BLOCKSIZE 384
#define MAXPRINTEDBLOCKS 0	//20

#define OLS_STRIDE 15
#ifdef ENABLE_SSE
#define SSEBITS 6
#endif
//#define USE_PROB_BITS PROB_BITS
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif

//RCT1~7:
//Cr-=y		Cr = [ 1	-1	0].RGB
//y+=Cr>>1	y  = [ 1/2	 1/2	0].RGB
//Cb-=y		Cb = [-1/2	-1/2	1].RGB
//y+=LUMA_UPDATE(Cb, Cr)
	#define LUMA_UPDATE_1(Cb, Cr) ((Cb)>>1)			//RCT1	y = [1/4	1/4	 1/2].RGB
	#define LUMA_UPDATE_2(Cb, Cr) ((2*(Cb)-(Cr)+4)>>3)	//RCT2	y = [1/4	1/2	 1/4].RGB	coincidentally, A710 from GFWX
	#define LUMA_UPDATE_3(Cb, Cr) ((2*(Cb)+(Cr)+4)>>3)	//RCT3a	y = [1/2	1/4	 1/4].RGB	never used
//	#define LUMA_UPDATE_3(Cb, Cr) (((Cb)+(Cr)+2)>>2)	//RCT3b	y = [5/8	1/8	 1/4].RGB
	#define LUMA_UPDATE_4(Cb, Cr) ((Cb)/3)			//RCT4	y = [1/3	1/3	 1/3].RGB
	#define LUMA_UPDATE_5(Cb, Cr) ((3*(Cb)+4)>>3)		//RCT5	y = [5/16	5/16	 6/16].RGB
	#define LUMA_UPDATE_6(Cb, Cr) ((7*(Cb)+8)>>4)		//RCT6	y = [9/32	9/32	14/32].RGB
	#define LUMA_UPDATE_7(Cb, Cr) ((10*(Cb)-(Cr)+16)>>5)	//RCT7	y = [5/16	6/16	 5/16].RGB
static const int rgb2yuv_permutations[6][3]=
{
	//YUV
	1, 2, 0,
	2, 0, 1,
	0, 1, 2,
	2, 1, 0,
	1, 0, 2,
	0, 2, 1,
};
typedef enum _IChannelType
{
	ICH_R,
	ICH_G,
	ICH_B,

	//YCbCr-R "RCT1~7"		FIXME deduplicate channels
	ICH_R1_120Y,
	ICH_R2_120Y,
	ICH_R3_120Y,
	ICH_R4_120Y,
	ICH_R5_120Y,
	ICH_R6_120Y,
	ICH_R7_120Y,
	ICH_R1_120U,
	ICH_R1_120V,

	ICH_R1_201Y,
	ICH_R2_201Y,
	ICH_R3_201Y,
	ICH_R4_201Y,
	ICH_R5_201Y,
	ICH_R6_201Y,
	ICH_R7_201Y,
	ICH_R1_201U,
	ICH_R1_201V,

	ICH_R1_012Y,
	ICH_R2_012Y,
	ICH_R3_012Y,
	ICH_R4_012Y,
	ICH_R5_012Y,
	ICH_R6_012Y,
	ICH_R7_012Y,
	ICH_R1_012U,
	ICH_R1_012V,

	ICH_R1_210Y,
	ICH_R2_210Y,
	ICH_R3_210Y,
	ICH_R4_210Y,
	ICH_R5_210Y,
	ICH_R6_210Y,
	ICH_R7_210Y,
	ICH_R1_210U,
	ICH_R1_210V,

	ICH_R1_102Y,
	ICH_R2_102Y,
	ICH_R3_102Y,
	ICH_R4_102Y,
	ICH_R5_102Y,
	ICH_R6_102Y,
	ICH_R7_102Y,
	ICH_R1_102U,
	ICH_R1_102V,

	ICH_R1_021Y,
	ICH_R2_021Y,
	ICH_R3_021Y,
	ICH_R4_021Y,
	ICH_R5_021Y,
	ICH_R6_021Y,
	ICH_R7_021Y,
	ICH_R1_021U,
	ICH_R1_021V,

	//RCT_Pei09
	ICH_P9_120Y,
	ICH_P9_120U,
	ICH_P9_120V,
	ICH_P9_201Y,
	ICH_P9_201U,
	ICH_P9_201V,
	ICH_P9_012Y,
	ICH_P9_012U,
	ICH_P9_012V,
	ICH_P9_210Y,
	ICH_P9_210U,
	ICH_P9_210V,
	ICH_P9_102Y,
	ICH_P9_102U,
	ICH_P9_102V,
	ICH_P9_021Y,
	ICH_P9_021U,
	ICH_P9_021V,

	//RCT_JPEG2000
	ICH_J2_120Y,
	ICH_J2_120U,
	ICH_J2_120V,
	ICH_J2_201Y,
	ICH_J2_201U,
	ICH_J2_201V,
	ICH_J2_012Y,
	ICH_J2_012U,
	ICH_J2_012V,

	ICH_ZERO,
	ICH_COUNT,
} IChannelType;
typedef enum _InfoIndex
{
	II_TARGET,
	II_HELPER1,
	II_HELPER2,
	II_HSHIFT,
	II_INFLATION,

	II_COUNT,
} InfoIndex;
#define OCHLIST\
	OCH(OCH_R,		0,	ICH_R,		ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_G,		0,	ICH_G,		ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_B,		0,	ICH_B,		ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_RG,		0,	ICH_R,		ICH_G,		ICH_ZERO,	0)\
	OCH(OCH_RB,		0,	ICH_R,		ICH_B,		ICH_ZERO,	0)\
	OCH(OCH_GR,		0,	ICH_G,		ICH_R,		ICH_ZERO,	0)\
	OCH(OCH_GB,		0,	ICH_G,		ICH_B,		ICH_ZERO,	0)\
	OCH(OCH_BG,		0,	ICH_B,		ICH_G,		ICH_ZERO,	0)\
	OCH(OCH_BR,		0,	ICH_B,		ICH_R,		ICH_ZERO,	0)\
	OCH(OCH_R2,		0,	ICH_R,		ICH_G,		ICH_B,		1)\
	OCH(OCH_G2,		0,	ICH_G,		ICH_B,		ICH_R,		1)\
	OCH(OCH_B2,		0,	ICH_B,		ICH_R,		ICH_G,		1)\
	OCH(OCH_R1_120Y,	0,	ICH_R1_120Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R2_120Y,	0,	ICH_R2_120Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R3_120Y,	0,	ICH_R3_120Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R4_120Y,	0,	ICH_R4_120Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R5_120Y,	0,	ICH_R5_120Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R6_120Y,	0,	ICH_R6_120Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R7_120Y,	0,	ICH_R7_120Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_120U,	1,	ICH_R1_120U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_120V,	1,	ICH_R1_120V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_201Y,	0,	ICH_R1_201Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R2_201Y,	0,	ICH_R2_201Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R3_201Y,	0,	ICH_R3_201Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R4_201Y,	0,	ICH_R4_201Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R5_201Y,	0,	ICH_R5_201Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R6_201Y,	0,	ICH_R6_201Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R7_201Y,	0,	ICH_R7_201Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_201U,	1,	ICH_R1_201U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_201V,	1,	ICH_R1_201V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_012Y,	0,	ICH_R1_012Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R2_012Y,	0,	ICH_R2_012Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R3_012Y,	0,	ICH_R3_012Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R4_012Y,	0,	ICH_R4_012Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R5_012Y,	0,	ICH_R5_012Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R6_012Y,	0,	ICH_R6_012Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R7_012Y,	0,	ICH_R7_012Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_012U,	1,	ICH_R1_012U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_012V,	1,	ICH_R1_012V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_210Y,	0,	ICH_R1_210Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R2_210Y,	0,	ICH_R2_210Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R3_210Y,	0,	ICH_R3_210Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R4_210Y,	0,	ICH_R4_210Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R5_210Y,	0,	ICH_R5_210Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R6_210Y,	0,	ICH_R6_210Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R7_210Y,	0,	ICH_R7_210Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_210U,	1,	ICH_R1_210U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_210V,	1,	ICH_R1_210V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_102Y,	0,	ICH_R1_102Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R2_102Y,	0,	ICH_R2_102Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R3_102Y,	0,	ICH_R3_102Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R4_102Y,	0,	ICH_R4_102Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R5_102Y,	0,	ICH_R5_102Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R6_102Y,	0,	ICH_R6_102Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R7_102Y,	0,	ICH_R7_102Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_102U,	1,	ICH_R1_102U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_102V,	1,	ICH_R1_102V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_021Y,	0,	ICH_R1_021Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R2_021Y,	0,	ICH_R2_021Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R3_021Y,	0,	ICH_R3_021Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R4_021Y,	0,	ICH_R4_021Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R5_021Y,	0,	ICH_R5_021Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R6_021Y,	0,	ICH_R6_021Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R7_021Y,	0,	ICH_R7_021Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_021U,	1,	ICH_R1_021U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_R1_021V,	1,	ICH_R1_021V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_120Y,	0,	ICH_P9_120Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_120U,	1,	ICH_P9_120U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_120V,	1,	ICH_P9_120V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_201Y,	0,	ICH_P9_201Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_201U,	1,	ICH_P9_201U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_201V,	1,	ICH_P9_201V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_012Y,	0,	ICH_P9_012Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_012U,	1,	ICH_P9_012U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_012V,	1,	ICH_P9_012V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_210Y,	0,	ICH_P9_210Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_210U,	1,	ICH_P9_210U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_210V,	1,	ICH_P9_210V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_102Y,	0,	ICH_P9_102Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_102U,	1,	ICH_P9_102U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_102V,	1,	ICH_P9_102V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_021Y,	0,	ICH_P9_021Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_021U,	1,	ICH_P9_021U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_P9_021V,	1,	ICH_P9_021V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_J2_120Y,	0,	ICH_J2_120Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_J2_120U,	1,	ICH_J2_120U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_J2_120V,	1,	ICH_J2_120V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_J2_201Y,	0,	ICH_J2_201Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_J2_201U,	1,	ICH_J2_201U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_J2_201V,	1,	ICH_J2_201V,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_J2_012Y,	0,	ICH_J2_012Y,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_J2_012U,	1,	ICH_J2_012U,	ICH_ZERO,	ICH_ZERO,	0)\
	OCH(OCH_J2_012V,	1,	ICH_J2_012V,	ICH_ZERO,	ICH_ZERO,	0)
typedef enum _OChannelType
{
#define OCH(ONAME, OINF, TARGET, HELPER1, HELPER2, HSHIFT) ONAME,
	OCHLIST
#undef  OCH
	OCH_COUNT,
} OChannelType;
static const int och_info[OCH_COUNT][II_COUNT]=
{
#define OCH(ONAME, OINF, TARGET, HELPER1, HELPER2, HSHIFT) {TARGET, HELPER1, HELPER2, HSHIFT, OINF},
	OCHLIST
#undef  OCH
};
static const char *och_names[OCH_COUNT]=
{
#define OCH(ONAME, OINF, TARGET, HELPER1, HELPER2, HSHIFT) #ONAME,
	OCHLIST
#undef  OCH
};
#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		4, 4, 4,	4, 4, 4)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		4, 4, 1,	4, 4, 4)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		4, 4, 0,	4, 4, 4)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		4, 4, 0,	4, 4, 4)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		4, 4, 1,	4, 4, 4)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		4, 4, 1,	4, 4, 4)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		4, 4, 0,	4, 4, 4)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		4, 0, 0,	4, 4, 4)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		4, 0, 1,	4, 4, 4)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		4, 0, 1,	4, 4, 4)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		4, 0, 0,	4, 4, 4)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		4, 0, 1,	4, 4, 4)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		4, 0, 1,	4, 4, 4)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		4, 0, 0,	4, 4, 4)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		4, 0, 1,	4, 4, 4)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		4, 0, 1,	4, 4, 4)\
	RCT(R_G_B2,	OCH_R,		OCH_G,		OCH_B2,		4, 4, 0,	4, 4, 1)\
	RCT(R_B_G2,	OCH_R,		OCH_B,		OCH_G2,		4, 4, 0,	4, 4, 1)\
	RCT(G_B_R2,	OCH_G,		OCH_B,		OCH_R2,		4, 4, 0,	4, 4, 1)\
	RCT(R_GR_B2,	OCH_R,		OCH_GR,		OCH_B2,		4, 0, 0,	4, 4, 1)\
	RCT(R_BR_G2,	OCH_R,		OCH_BR,		OCH_G2,		4, 0, 0,	4, 4, 1)\
	RCT(G_BG_R2,	OCH_G,		OCH_BG,		OCH_R2,		4, 0, 0,	4, 4, 1)\
	RCT(G_RG_B2,	OCH_G,		OCH_RG,		OCH_B2,		4, 0, 0,	4, 4, 1)\
	RCT(B_RB_G2,	OCH_B,		OCH_RB,		OCH_G2,		4, 0, 0,	4, 4, 1)\
	RCT(B_GB_R2,	OCH_B,		OCH_GB,		OCH_R2,		4, 0, 0,	4, 4, 1)\
	RCT(RCT1_120,	OCH_R1_120Y,	OCH_R1_120U,	OCH_R1_120V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT1_201,	OCH_R1_201Y,	OCH_R1_201U,	OCH_R1_201V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT1_012,	OCH_R1_012Y,	OCH_R1_012U,	OCH_R1_012V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT1_210,	OCH_R1_210Y,	OCH_R1_210U,	OCH_R1_210V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT1_102,	OCH_R1_102Y,	OCH_R1_102U,	OCH_R1_102V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT1_021,	OCH_R1_021Y,	OCH_R1_021U,	OCH_R1_021V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT2_120,	OCH_R2_120Y,	OCH_R1_120U,	OCH_R1_120V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT2_201,	OCH_R2_201Y,	OCH_R1_201U,	OCH_R1_201V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT2_012,	OCH_R2_012Y,	OCH_R1_012U,	OCH_R1_012V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT2_210,	OCH_R2_210Y,	OCH_R1_210U,	OCH_R1_210V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT2_102,	OCH_R2_102Y,	OCH_R1_102U,	OCH_R1_102V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT2_021,	OCH_R2_021Y,	OCH_R1_021U,	OCH_R1_021V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT3_120,	OCH_R3_120Y,	OCH_R1_120U,	OCH_R1_120V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT3_201,	OCH_R3_201Y,	OCH_R1_201U,	OCH_R1_201V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT3_012,	OCH_R3_012Y,	OCH_R1_012U,	OCH_R1_012V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT3_210,	OCH_R3_210Y,	OCH_R1_210U,	OCH_R1_210V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT3_102,	OCH_R3_102Y,	OCH_R1_102U,	OCH_R1_102V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT3_021,	OCH_R3_021Y,	OCH_R1_021U,	OCH_R1_021V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT4_120,	OCH_R4_120Y,	OCH_R1_120U,	OCH_R1_120V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT4_201,	OCH_R4_201Y,	OCH_R1_201U,	OCH_R1_201V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT4_012,	OCH_R4_012Y,	OCH_R1_012U,	OCH_R1_012V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT4_210,	OCH_R4_210Y,	OCH_R1_210U,	OCH_R1_210V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT4_102,	OCH_R4_102Y,	OCH_R1_102U,	OCH_R1_102V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT4_021,	OCH_R4_021Y,	OCH_R1_021U,	OCH_R1_021V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT5_120,	OCH_R5_120Y,	OCH_R1_120U,	OCH_R1_120V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT5_201,	OCH_R5_201Y,	OCH_R1_201U,	OCH_R1_201V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT5_012,	OCH_R5_012Y,	OCH_R1_012U,	OCH_R1_012V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT5_210,	OCH_R5_210Y,	OCH_R1_210U,	OCH_R1_210V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT5_102,	OCH_R5_102Y,	OCH_R1_102U,	OCH_R1_102V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT5_021,	OCH_R5_021Y,	OCH_R1_021U,	OCH_R1_021V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT6_120,	OCH_R6_120Y,	OCH_R1_120U,	OCH_R1_120V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT6_201,	OCH_R6_201Y,	OCH_R1_201U,	OCH_R1_201V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT6_012,	OCH_R6_012Y,	OCH_R1_012U,	OCH_R1_012V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT6_210,	OCH_R6_210Y,	OCH_R1_210U,	OCH_R1_210V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT6_102,	OCH_R6_102Y,	OCH_R1_102U,	OCH_R1_102V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT6_021,	OCH_R6_021Y,	OCH_R1_021U,	OCH_R1_021V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT7_120,	OCH_R7_120Y,	OCH_R1_120U,	OCH_R1_120V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT7_201,	OCH_R7_201Y,	OCH_R1_201U,	OCH_R1_201V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT7_012,	OCH_R7_012Y,	OCH_R1_012U,	OCH_R1_012V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT7_210,	OCH_R7_210Y,	OCH_R1_210U,	OCH_R1_210V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT7_102,	OCH_R7_102Y,	OCH_R1_102U,	OCH_R1_102V,	4, 4, 4,	4, 4, 4)\
	RCT(RCT7_021,	OCH_R7_021Y,	OCH_R1_021U,	OCH_R1_021V,	4, 4, 4,	4, 4, 4)\
	RCT(Pei09_120,	OCH_P9_120Y,	OCH_P9_120U,	OCH_P9_120V,	4, 4, 4,	4, 4, 4)\
	RCT(Pei09_201,	OCH_P9_201Y,	OCH_P9_201U,	OCH_P9_201V,	4, 4, 4,	4, 4, 4)\
	RCT(Pei09_012,	OCH_P9_012Y,	OCH_P9_012U,	OCH_P9_012V,	4, 4, 4,	4, 4, 4)\
	RCT(Pei09_210,	OCH_P9_210Y,	OCH_P9_210U,	OCH_P9_210V,	4, 4, 4,	4, 4, 4)\
	RCT(Pei09_102,	OCH_P9_102Y,	OCH_P9_102U,	OCH_P9_102V,	4, 4, 4,	4, 4, 4)\
	RCT(Pei09_021,	OCH_P9_021Y,	OCH_P9_021U,	OCH_P9_021V,	4, 4, 4,	4, 4, 4)\
	RCT(J2K_120,	OCH_J2_120Y,	OCH_J2_120U,	OCH_J2_120V,	4, 4, 4,	4, 4, 4)\
	RCT(J2K_201,	OCH_J2_201Y,	OCH_J2_201U,	OCH_J2_201V,	4, 4, 4,	4, 4, 4)\
	RCT(J2K_012,	OCH_J2_012Y,	OCH_J2_012U,	OCH_J2_012V,	4, 4, 4,	4, 4, 4)
typedef enum _RCTType
{
#define RCT(NAME, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2) RCT_##NAME,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTType;
static const int rct_combinations[RCT_COUNT][12]=
{
#define RCT(NAME, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2) {YIDX, UIDX, VIDX, 0, YOFF1, UOFF1, VOFF1, 4, YOFF2, UOFF2, VOFF2, 4},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(NAME, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2) #NAME,
	RCTLIST
#undef  RCT
};


typedef struct _OLS4Context
{
	int nparams;
	double *nb, *vec, *cov, *cholesky, *params;
} OLS4Context;
static void ols4_init(OLS4Context *ctx, int nparams)
{
	int msize=nparams*nparams;
	memset(ctx, 0, sizeof(*ctx));
	ctx->nparams=nparams;
	ctx->nb=(double*)_mm_malloc(nparams*sizeof(double), sizeof(__m256d));
	ctx->vec=(double*)_mm_malloc(nparams*sizeof(double), sizeof(__m256d));
	ctx->cov=(double*)_mm_malloc(msize*sizeof(double), sizeof(__m256d));
	ctx->cholesky=(double*)malloc(msize*sizeof(double));
	ctx->params=(double*)_mm_malloc(nparams*sizeof(double), sizeof(__m256d));
	if(!ctx->nb||!ctx->vec||!ctx->cov||!ctx->cholesky||!ctx->params)
	{
		LOG_ERROR("Alloc error");
		return;
	}
}
static void ols4_reset(OLS4Context *ctx)
{
	int msize=ctx->nparams*ctx->nparams;
	memset(ctx->nb, 0, ctx->nparams*sizeof(double));
	memset(ctx->vec, 0, ctx->nparams*sizeof(double));
	memset(ctx->cov, 0, msize*sizeof(double));
	memset(ctx->cholesky, 0, msize*sizeof(double));
	memset(ctx->params, 0, ctx->nparams*sizeof(double));
}
static void ols4_free(OLS4Context *ctx)
{
	_mm_free(ctx->nb);
	_mm_free(ctx->vec);
	_mm_free(ctx->cov);
	free(ctx->cholesky);
	_mm_free(ctx->params);
	memset(ctx, 0, sizeof(*ctx));
}
static int ols4_predict(OLS4Context *ctx, const int *nb, int cmin, int cmax)
{
	double fpred=0;
#ifdef ALLOW_AVX2
	int k=0;
	double temp;
	ALIGN(32) double asum[4];
	__m256d msum=_mm256_setzero_pd();
	for(;k<ctx->nparams-3;k+=4)
	{
		__m128i mnb=_mm_load_si128((__m128i*)(nb+k));
		__m256d mp=_mm256_load_pd(ctx->params+k);
		__m256d dnb=_mm256_cvtepi32_pd(mnb);
		_mm256_store_pd(ctx->nb+k, dnb);
		dnb=_mm256_mul_pd(dnb, mp);
		msum=_mm256_add_pd(msum, dnb);
	}
	_mm256_store_pd(asum, msum);
	fpred=asum[0]+asum[1]+asum[2]+asum[3];
	switch(ctx->nparams-k)
	{
	case 3:temp=nb[k]; ctx->nb[k]=temp; fpred+=ctx->params[k]*temp; ++k;
	case 2:temp=nb[k]; ctx->nb[k]=temp; fpred+=ctx->params[k]*temp; ++k;
	case 1:temp=nb[k]; ctx->nb[k]=temp; fpred+=ctx->params[k]*temp;
		break;
	}
#else
	for(int k=0;k<ctx->nparams;++k)
	{
		double val=nb[k];
		ctx->nb[k]=val;
		fpred+=ctx->params[k]*val;
	}
#endif
	{
		ALIGN(16) int pred[4];
		__m128i x=_mm_cvtpd_epi32(_mm_set_sd(fpred));
		__m128i mmin=_mm_set_epi32(0, 0, 0, cmin);
		__m128i mmax=_mm_set_epi32(0, 0, 0, cmax);
		x=_mm_max_epi16(x, mmin);
		x=_mm_min_epi16(x, mmax);
		_mm_store_si128((__m128i*)pred, x);
		return pred[0];
	}
}
static void ols4_update(OLS4Context *ctx, int target, double lr)
{
	double lval=target*lr, lr_comp=1-lr;//
#ifdef ALLOW_AVX2
	__m256d mlr=_mm256_set1_pd(lr);
	for(int ky2=0, midx=0;ky2<ctx->nparams;++ky2)
	{
		__m256d yctx=_mm256_set1_pd(ctx->nb[ky2]);
		int kx2=0;
		for(;kx2<ctx->nparams-3;kx2+=4, midx+=4)
		{
			__m256d mcov=_mm256_load_pd(ctx->cov+midx);
			__m256d xctx=_mm256_load_pd(ctx->nb+kx2);
			xctx=_mm256_mul_pd(xctx, yctx);
			xctx=_mm256_sub_pd(xctx, mcov);
			xctx=_mm256_mul_pd(xctx, mlr);
			xctx=_mm256_add_pd(xctx, mcov);
			_mm256_store_pd(ctx->cov+midx, xctx);
			//curr_cov[midx+0]+=(ctx->nb[kx2+0]*ctx->nb[ky2]-curr_cov[midx+0])*lr;
			//curr_cov[midx+1]+=(ctx->nb[kx2+1]*ctx->nb[ky2]-curr_cov[midx+1])*lr;
			//curr_cov[midx+2]+=(ctx->nb[kx2+2]*ctx->nb[ky2]-curr_cov[midx+2])*lr;
			//curr_cov[midx+3]+=(ctx->nb[kx2+3]*ctx->nb[ky2]-curr_cov[midx+3])*lr;
			//for(int k=0;k<4;++k)
			//{
			//	if(abs(curr_cov[midx+k]-xctx.m256d_f64[k])>1e-6)
			//		printf("");
			//}
		}
		for(;kx2<ctx->nparams;++kx2, ++midx)
			ctx->cov[midx]+=(ctx->nb[kx2]*ctx->nb[ky2]-ctx->cov[midx])*lr;
	}
	{
		__m256d mlr_comp=_mm256_set1_pd(lr_comp);
		__m256d mval=_mm256_set1_pd(lval);
		int k=0;
		for(;k<ctx->nparams-3;k+=4)
		{
			__m256d mvec=_mm256_load_pd(ctx->vec+k);
			__m256d mtmp=_mm256_load_pd(ctx->nb+k);
			mtmp=_mm256_mul_pd(mval, mtmp);
			mvec=_mm256_mul_pd(mvec, mlr_comp);
			mvec=_mm256_add_pd(mvec, mtmp);
			_mm256_store_pd(ctx->vec+k, mvec);
			//ctx->vec[k+0]=lval*ctx->nb[k+0]+lr_comp*ctx->vec[k+0];
			//ctx->vec[k+1]=lval*ctx->nb[k+1]+lr_comp*ctx->vec[k+1];
			//ctx->vec[k+2]=lval*ctx->nb[k+2]+lr_comp*ctx->vec[k+2];
			//ctx->vec[k+3]=lval*ctx->nb[k+3]+lr_comp*ctx->vec[k+3];
		}
		for(;k<ctx->nparams;++k)
			ctx->vec[k]=lval*ctx->nb[k]+lr_comp*ctx->vec[k];
	}
#else
	for(int ky2=0, midx=0;ky2<ctx->nparams;++ky2)
	{
		for(int kx2=0;kx2<ctx->nparams;++kx2, ++midx)
			ctx->cov[midx]+=(ctx->nb[kx2]*ctx->nb[ky2]-ctx->cov[midx])*lr;
	}
	for(int k=0;k<ctx->nparams;++k)
		ctx->vec[k]=lval*ctx->nb[k]+lr_comp*ctx->vec[k];
#endif
	int n=ctx->nparams, msize=n*n, success=1;
	memcpy(ctx->cholesky, ctx->cov, sizeof(double)*msize);
	for(int k=0;k<msize;k+=n+1)
		ctx->cholesky[k]+=0.0075;
	double sum;
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<i;++j)
		{
			sum=ctx->cholesky[i*n+j];
			for(int k=0;k<j;++k)
				sum-=ctx->cholesky[i*n+k]*ctx->cholesky[j*n+k];
			ctx->cholesky[i*n+j]=sum/ctx->cholesky[j*n+j];
		}
		sum=ctx->cholesky[i*n+i];
		for(int k=0;k<i;++k)
			sum-=ctx->cholesky[i*n+k]*ctx->cholesky[i*n+k];
		if(sum<=1e-8)
		{
			success=0;
			break;
		}
		ctx->cholesky[i*n+i]=sqrt(sum);
	}
	if(success)
	{
		for(int i=0;i<n;++i)
		{
			sum=ctx->vec[i];
			for(int j=0;j<i;++j)
				sum-=ctx->cholesky[i*n+j]*ctx->params[j];
			ctx->params[i]=sum/ctx->cholesky[i*n+i];
		}
		for(int i=n-1;i>=0;--i)
		{
			sum=ctx->params[i];
			for(int j=i+1;j<n;++j)
				sum-=ctx->cholesky[j*n+i]*ctx->params[j];
			ctx->params[i]=sum/ctx->cholesky[i*n+i];
		}
	}
}

#define OLS4_NPARAMS 16		//multiple of 4,  NPARAMS is unified because of interleaving


//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
#define CONFIG_EXP 4
#define CONFIG_MSB 1
#define CONFIG_LSB 0
static void quantize_pixel(int val, int *token, int *bypass, int *nbits)
{
	if(val<(1<<CONFIG_EXP))
	{
		*token=val;//token
		*nbits=0;
		*bypass=0;
	}
	else
	{
		int lgv=FLOOR_LOG2((unsigned)val);
		int mantissa=val-(1<<lgv);
		*token = (1<<CONFIG_EXP) + (
				(lgv-CONFIG_EXP)<<(CONFIG_MSB+CONFIG_LSB)|
				(mantissa>>(lgv-CONFIG_MSB))<<CONFIG_LSB|
				(mantissa&((1<<CONFIG_LSB)-1))
			);
		*nbits=lgv-(CONFIG_MSB+CONFIG_LSB);
		*bypass=val>>CONFIG_LSB&((1LL<<*nbits)-1);
	}
}

typedef struct _ThreadArgs
{
	const Image *src;
	Image *dst;
	int fwd, loud, x1, x2, y1, y2;
	
	short *pixels;
	int bufsize;

	int *hist, histsize, histindices[OCH_COUNT+1];

	//AC
	int tlevels, clevels, statssize;
	//unsigned *stats;
	DList list;
	const unsigned char *decstart, *decend;
	
#ifdef ENABLE_OLS
	OLS4Context ols4[3];
#endif
#ifdef ENABLE_SSE
	long long sse[3][1<<SSEBITS][1<<SSEBITS];
#endif

	//aux
	int blockidx, bestrct;
	double usizes[3], csizes[3], bestsize;
} ThreadArgs;
#if 0
#define CDFSTRIDE 16
static void update_CDF(const int *hist, unsigned *CDF, int tlevels)
{
	int sum=hist[tlevels], c=0;
	for(int ks=0;ks<tlevels;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*((1LL<<USE_PROB_BITS)-tlevels)/sum)+ks;
		c+=freq;
	}
	CDF[tlevels]=1<<USE_PROB_BITS;

	//for(int ks=0;ks<tlevels;++ks)//
	//{
	//	if(CDF[ks]>CDF[ks+1])
	//		LOG_ERROR("");
	//}
}
#endif
static void block_thread(void *param)
{
	AC3 ec;
	ThreadArgs *args=(ThreadArgs*)param;
	Image const *image=args->fwd?args->src:args->dst;
	int depths[OCH_COUNT]={0}, halfs[OCH_COUNT]={0};
	double bestsize=0;
	int bestrct=0;
	const int *combination=0;
	int nctx=args->clevels*args->clevels, cdfstride=args->tlevels+1, chsize=nctx*cdfstride;
	//ALIGN(16) int ols4_ctx[OLS4_NPARAMS]={0};//
	ALIGN(16) int ols4_ctx0[OLS4_NPARAMS]={0}, ols4_ctx1[OLS4_NPARAMS]={0}, ols4_ctx2[OLS4_NPARAMS]={0};
	static const double lrs[]={0.0018, 0.003, 0.003};

	for(int kc=0;kc<OCH_COUNT;++kc)
	{
		int depth=image->depth+och_info[kc][II_INFLATION];
		UPDATE_MIN(depth, 16);
		depths[kc]=depth;
		halfs[kc]=1<<depth>>1;
	}
	if(args->fwd)//encode
	{
		if(image->nch>=3)
		{
			int res=((args->x2-args->x1)>>2)*((args->y2-args->y1)>>2);
			int nlevels[OCH_COUNT]={0};
			double csizes[OCH_COUNT]={0};

			for(int kc=0;kc<OCH_COUNT;++kc)
				nlevels[kc]=1<<depths[kc];
			
			//for(int k=0;k<image->nch&&k<3;++k)
			//	ols4_reset(args->ols4+k);
			memset(args->hist, 0, args->histsize);
			memset(args->pixels, 0, args->bufsize);
			for(int ky=args->y1, ay=0;ky<args->y2;ky+=4, ++ay)
			{
				ALIGN(16) short *rows[]=
				{
					args->pixels+((image->iw+16LL)*((ay-0LL)&3)+8LL)*OCH_COUNT,
					args->pixels+((image->iw+16LL)*((ay-1LL)&3)+8LL)*OCH_COUNT,
					args->pixels+((image->iw+16LL)*((ay-2LL)&3)+8LL)*OCH_COUNT,
					args->pixels+((image->iw+16LL)*((ay-3LL)&3)+8LL)*OCH_COUNT,
				};
				int idx=image->nch*(image->iw*ky+args->x1);
				short input[ICH_COUNT]={0};
				ALIGN(32) short preds[OCH_COUNT]={0};
				for(int kx=args->x1;kx<args->x2;kx+=4, idx+=image->nch*4)
				{
					int cidx=3;
					short
						*NNN	=rows[3]+0*OCH_COUNT,
						*NNWW	=rows[2]-2*OCH_COUNT,
						*NNW	=rows[2]-1*OCH_COUNT,
						*NN	=rows[2]+0*OCH_COUNT,
						*NNE	=rows[2]+1*OCH_COUNT,
						*NNEE	=rows[2]+2*OCH_COUNT,
						*NNEEE	=rows[2]+3*OCH_COUNT,
						*NWWW	=rows[1]-3*OCH_COUNT,
						*NWW	=rows[1]-2*OCH_COUNT,
						*NW	=rows[1]-1*OCH_COUNT,
						*N	=rows[1]+0*OCH_COUNT,
						*NE	=rows[1]+1*OCH_COUNT,
						*NEE	=rows[1]+2*OCH_COUNT,
						*NEEE	=rows[1]+3*OCH_COUNT,
						*WWWW	=rows[0]-4*OCH_COUNT,
						*WWW	=rows[0]-3*OCH_COUNT,
						*WW	=rows[0]-2*OCH_COUNT,
						*W	=rows[0]-1*OCH_COUNT,
						*curr	=rows[0]+0*OCH_COUNT;
#if 0
					if(ky<=args->y1+2)
					{
						if(ky<=args->y1+1)
						{
							if(ky==args->y1)
								NEEE=NEE=NE=NWW=NW=N=W;
							NNWW=NWW;
							NNW=NW;
							NN=N;
							NNE=NE;
							NNEE=NEE;
							NNEEE=NEEE;
						}
						NNN=NN;
					}
					if(kx<=args->x1+3)
					{
						if(kx<=args->x1+2)
						{
							if(kx<=args->x1+1)
							{
								if(kx<=args->x1)
									NW=W=N;
								WW=W;
								NWW=NW;
							}
							WWW=WW;
						}
						WWWW=WWW;
					}
					if(kx>=args->x2-3)
					{
						if(kx>=args->x2-2)
						{
							if(kx>=args->x2-1)
							{
								NNE=NN;
								NE=N;
							}
							NNEE=NNE;
							NEE=NE;
						}
						NEEE=NEE;
					}
#endif
					input[ICH_R]=image->data[idx+0];//r
					input[ICH_G]=image->data[idx+1];//g
					input[ICH_B]=image->data[idx+2];//b
					for(int kp=0;kp<6;++kp)//RCT1~7		r-=g; g+=r>>1; b-=g; g+=(A*r+B*b+C)>>SH;
					{
						//int j=0;
						input[cidx+0]=input[rgb2yuv_permutations[kp][0]];//Y
						input[cidx+7]=input[rgb2yuv_permutations[kp][1]];//U
						input[cidx+8]=input[rgb2yuv_permutations[kp][2]];//V
						input[cidx+8]-=input[cidx+0];				//	Cr = [ 1	-1	0].RGB
						input[cidx+0]+=input[cidx+8]>>1;			//	y = [1/2	1/2	0].RGB
						input[cidx+7]-=input[cidx+0];				//	Cb = [-1/2	-1/2	1].RGB
						input[cidx+6]=input[cidx+5]=input[cidx+4]=input[cidx+3]=input[cidx+2]=input[cidx+1]=input[cidx+0];
						input[cidx+0]+=LUMA_UPDATE_1(input[cidx+7], input[cidx+8]);
						input[cidx+1]+=LUMA_UPDATE_2(input[cidx+7], input[cidx+8]);
						input[cidx+2]+=LUMA_UPDATE_3(input[cidx+7], input[cidx+8]);
						input[cidx+3]+=LUMA_UPDATE_4(input[cidx+7], input[cidx+8]);
						input[cidx+4]+=LUMA_UPDATE_5(input[cidx+7], input[cidx+8]);
						input[cidx+5]+=LUMA_UPDATE_6(input[cidx+7], input[cidx+8]);
						input[cidx+6]+=LUMA_UPDATE_7(input[cidx+7], input[cidx+8]);
						cidx+=9;
					}
					for(int kp=0;kp<6;++kp)//Pei09		b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
					{
						input[cidx+0]=input[rgb2yuv_permutations[kp][0]];//Y
						input[cidx+1]=input[rgb2yuv_permutations[kp][1]];//U
						input[cidx+2]=input[rgb2yuv_permutations[kp][2]];//V
						input[cidx+1]-=(87*input[cidx+2]+169*input[cidx+0]+128)>>8;
						input[cidx+2]-=input[cidx+0];
						input[cidx+0]+=(86*input[cidx+2]+29*input[cidx+1]+128)>>8;
						cidx+=3;
					}
					for(int kp=0;kp<3;++kp)//JPEG2000	r-=g; b-=g; g+=(r+b)>>2;
					{
						input[cidx+0]=input[rgb2yuv_permutations[kp][0]];//Y
						input[cidx+1]=input[rgb2yuv_permutations[kp][1]];//U
						input[cidx+2]=input[rgb2yuv_permutations[kp][2]];//V
						input[cidx+2]-=input[cidx+0];
						input[cidx+1]-=input[cidx+0];
						input[cidx+0]+=(input[cidx+2]+input[cidx+1])>>2;
						cidx+=3;
					}
					{
						int kc=0;
						for(;kc<OCH_COUNT-15;kc+=16)
						{
							__m256i mNW	=_mm256_loadu_si256((__m256i*)(NW	+kc));
							__m256i mN	=_mm256_loadu_si256((__m256i*)(N	+kc));
							__m256i mW	=_mm256_loadu_si256((__m256i*)(W	+kc));
							__m256i mp=_mm256_sub_epi16(_mm256_add_epi16(mN, mW), mNW);
							__m256i vmin=_mm256_min_epi16(mN, mW);
							__m256i vmax=_mm256_max_epi16(mN, mW);
							mp=_mm256_max_epi16(mp, vmin);
							mp=_mm256_min_epi16(mp, vmax);
							_mm256_store_si256((__m256i*)(preds+kc), mp);
						}
						if(kc<OCH_COUNT-7)
						{
							__m128i mNW	=_mm_loadu_si128((__m128i*)(NW	+kc));
							__m128i mN	=_mm_loadu_si128((__m128i*)(N	+kc));
							__m128i mW	=_mm_loadu_si128((__m128i*)(W	+kc));
							__m128i mp=_mm_sub_epi16(_mm_add_epi16(mN, mW), mNW);
							__m128i vmin=_mm_min_epi16(mN, mW);
							__m128i vmax=_mm_max_epi16(mN, mW);
							mp=_mm_max_epi16(mp, vmin);
							mp=_mm_min_epi16(mp, vmax);
							_mm_store_si128((__m128i*)(preds+kc), mp);
							kc+=8;
						}
						for(;kc<OCH_COUNT;++kc)
							MEDIAN3_32(preds[kc], N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
					}
					for(int kc=0;kc<OCH_COUNT;++kc)
					{
						int
							pred,
							offset=(input[och_info[kc][II_HELPER1]]+input[och_info[kc][II_HELPER2]])>>och_info[kc][II_HSHIFT],
							target=input[och_info[kc][II_TARGET]];

						pred=preds[kc];
						//MEDIAN3_32(pred, N[kc], W[kc], N[kc]+W[kc]-NW[kc]);

						//int j=0, cmin, cmax;
						//{
						//	ALIGN(16) int arr[4];
						//	__m128i va=_mm_set_epi32(0, 0, 0, N[kc]);
						//	__m128i vb=_mm_set_epi32(0, 0, 0, W[kc]);
						//	__m128i vc=_mm_set_epi32(0, 0, 0, N[kc]+W[kc]-NW[kc]);
						//	__m128i _vmin=_mm_min_epi32(va, vb);
						//	__m128i _vmax=_mm_max_epi32(va, vb);
						//	vc=_mm_max_epi32(vc, _vmin);
						//	vc=_mm_min_epi32(vc, _vmax);
						//	_mm_store_si128((__m128i*)arr, vc);
						//	pred=arr[0];
						//
						//	vc=_mm_set_epi32(0, 0, 0, NE[kc]);
						//	_vmin=_mm_min_epi32(_vmin, vc);
						//	_vmax=_mm_max_epi32(_vmax, vc);
						//	_mm_store_si128((__m128i*)arr, _vmin);
						//	cmin=arr[0];
						//	_mm_store_si128((__m128i*)arr, _vmax);
						//	cmax=arr[0];
						//}
						//ols4_ctx[j++]=pred;
						//ols4_ctx[j++]=NNWW	[kc];
						//ols4_ctx[j++]=NNW	[kc];
						//ols4_ctx[j++]=NN	[kc];
						//ols4_ctx[j++]=NNE	[kc];
						//ols4_ctx[j++]=NNEE	[kc];
						//ols4_ctx[j++]=NNEEE	[kc];
						//ols4_ctx[j++]=NWWW	[kc];
						//ols4_ctx[j++]=NWW	[kc];
						//ols4_ctx[j++]=NW	[kc];
						//ols4_ctx[j++]=N		[kc];
						//ols4_ctx[j++]=NE	[kc];
						//ols4_ctx[j++]=NEE	[kc];
						//ols4_ctx[j++]=NEEE	[kc];
						//ols4_ctx[j++]=WWW	[kc];
						//ols4_ctx[j++]=WW	[kc];
						//ols4_ctx[j++]=W		[kc];
						//pred=ols4_predict(args->ols4+0, ols4_ctx, cmin, cmax);

						pred+=offset;
						CLAMP2_32(pred, pred, -halfs[kc], halfs[kc]-1);

						int *curr_hist=args->hist+args->histindices[kc];
						int val=target-pred;
						val+=halfs[kc];
						val&=nlevels[kc]-1;
						++curr_hist[val];

						target-=offset;
						curr[kc]=target;
						//if(kc<3&&!(kx&63))
						//	ols4_update(args->ols4+kc, target, lrs[kc]);
					}
					rows[0]+=OCH_COUNT;
					rows[1]+=OCH_COUNT;
					rows[2]+=OCH_COUNT;
					rows[3]+=OCH_COUNT;
				}
			}
			for(int kc=0;kc<OCH_COUNT;++kc)//calculate channel sizes
			{
				int *curr_hist=args->hist+args->histindices[kc];
				int nlevels2=nlevels[kc];
				for(int ks=0;ks<nlevels2;++ks)
				{
					int freq=curr_hist[ks];
					if(freq)
						csizes[kc]-=freq*log2((double)freq/res);
				}
				csizes[kc]/=8;
			}
			for(int kt=0;kt<RCT_COUNT;++kt)//select best RCT
			{
				const int *group=rct_combinations[kt];
				double csize=
					csizes[group[0]]+
					csizes[group[1]]+
					csizes[group[2]];
				if(!kt||bestsize>csize)
					bestsize=csize, bestrct=kt;
			}
			combination=rct_combinations[bestrct];

			args->bestrct=bestrct;
			args->usizes[0]=(res*image->depth+7)>>3;
			args->usizes[1]=(res*image->depth+7)>>3;
			args->usizes[2]=(res*image->depth+7)>>3;
			args->csizes[0]=csizes[combination[0]];
			args->csizes[1]=csizes[combination[1]];
			args->csizes[2]=csizes[combination[2]];
			args->bestsize=bestsize;
			if(args->loud)
			{
				printf("X %5d~%5d  Y %5d~%5d  best %12.2lf bytes  %s\n",
					args->x1, args->x2, args->y1, args->y2,
					bestsize,
					rct_names[bestrct]
				);

				for(int kc=0;kc<OCH_COUNT;++kc)
					printf(" %12.2lf  %s\n", csizes[kc], och_names[kc]);

				for(int kt=0;kt<RCT_COUNT;++kt)
				{
					const int *group=rct_combinations[kt];
					double csize=
						csizes[group[0]]+
						csizes[group[1]]+
						csizes[group[2]];
					printf("%12.2lf %c  %-10s\n",
						csize,
						kt==bestrct?'*':' ',
						rct_names[kt]
					);
				}
			}
		}
		dlist_init(&args->list, 1, 1024, 0);
		ac3_enc_init(&ec, &args->list);
		if(image->nch>=3)
			ac3_enc_bypass_NPOT(&ec, bestrct, RCT_COUNT);

	}
	else//decode
	{
		const unsigned char *srcstart=args->decstart, *srcend=args->decend;

		ac3_dec_init(&ec, srcstart, srcend);
		if(image->nch>=3)
		{
			bestrct=ac3_dec_bypass_NPOT(&ec, RCT_COUNT);
			combination=rct_combinations[bestrct];
		}
	}
	if(image->nch<3)
	{
		bestrct=0;
		combination=rct_combinations[bestrct];
	}
	memset(args->hist, 0, args->statssize);
	//for(int kc=0;kc<image->nch;++kc)
	//{
	//	int *curr_hist=args->hist+chsize*kc;
	//	unsigned *curr_CDF=args->stats+chsize*kc;
	//	
	//	*curr_hist=1;//init bypass
	//	memfill(curr_hist+1, curr_hist, sizeof(int)*(args->tlevels-1LL), sizeof(int));
	//	curr_hist[args->tlevels]=args->tlevels;
	//	update_CDF(curr_hist, curr_CDF, args->tlevels);
	//	memfill(curr_hist+cdfstride, curr_hist, ((size_t)chsize-cdfstride)*sizeof(int), cdfstride*sizeof(int));
	//	memfill(curr_CDF+cdfstride, curr_CDF, ((size_t)chsize-cdfstride)*sizeof(int), cdfstride*sizeof(int));
	//}
#ifdef ENABLE_OLS
	int kols=0, kols2=OLS_STRIDE;
	for(int k=0;k<image->nch&&k<3;++k)
		ols4_reset(args->ols4+k);
#endif
#ifdef ENABLE_SSE
	memset(args->sse, 0, sizeof(args->sse));
#endif
	memset(args->pixels, 0, args->bufsize);
	int depths2[]=
	{
		depths[combination[0]],
		depths[combination[1]],
		depths[combination[2]],
	};
	int depths3[]=
	{
		depths2[0]-8,
		depths2[1]-8,
		depths2[2]-8,
	};
	int halfs2[]=
	{
		halfs[combination[0]],
		halfs[combination[1]],
		halfs[combination[2]],
	};
	__m128i swap16=_mm_set_epi8(
		13, 12, 15, 14,  9,  8, 11, 10,  5,  4,  7,  6,  1,  0,  3,  2
	//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
	);
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) short *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		int idx=image->nch*(image->iw*ky+args->x1);
		short yuv[5]={0};
		int token=0, bypass=0, nbits=0;

		for(int kx=args->x1;kx<args->x2;++kx, idx+=image->nch)
		{
			short
				*NNN	=rows[3]+0*4*2,
				*NNWW	=rows[2]-2*4*2,
				*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NNEE	=rows[2]+2*4*2,
				*NNEEE	=rows[2]+3*4*2,
				*NWWW	=rows[1]-3*4*2,
				*NWW	=rows[1]-2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*NEE	=rows[1]+2*4*2,
				*NEEE	=rows[1]+3*4*2,
				*WWWW	=rows[0]-4*4*2,
				*WWW	=rows[0]-3*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
			if(ky<=args->y1+2)
			{
				if(ky<=args->y1+1)
				{
					if(ky==args->y1)
						NEEE=NEE=NE=NWWW=NWW=NW=N=W;
					NNWW=NWW;
					NNW=NW;
					NN=N;
					NNE=NE;
					NNEE=NEE;
					NNEEE=NEEE;
				}
				NNN=NN;
			}
			if(kx<=args->x1+3)
			{
				if(kx<=args->x1+2)
				{
					if(kx<=args->x1+1)
					{
						if(kx<=args->x1)
							NW=W=N;
						WW=W;
						NWW=NW;
					}
					WWW=WW;
				}
				WWWW=WWW;
			}
			if(kx>=args->x2-3)
			{
				if(kx>=args->x2-2)
				{
					if(kx>=args->x2-1)
					{
						NNE=NN;
						NE=N;
					}
					NNEE=NNE;
					NEE=NE;
				}
				NEEE=NEE;
			}
			if(args->fwd)
			{
				int kc=0;
				if(image->nch>=3)
				{
					switch(bestrct)
					{
					case RCT_R_G_B:
					case RCT_R_G_BG:
					case RCT_R_G_BR:
					case RCT_R_GR_BR:
					case RCT_R_GR_BG:
					case RCT_R_G_B2:
					case RCT_R_GR_B2:
						yuv[0]=image->data[idx+0];
						yuv[1]=image->data[idx+1];
						yuv[2]=image->data[idx+2];
						break;
					case RCT_G_B_RG:
					case RCT_G_B_RB:
					case RCT_G_BG_RG:
					case RCT_G_BG_RB:
					case RCT_G_B_R2:
					case RCT_G_BG_R2:
						yuv[0]=image->data[idx+1];
						yuv[1]=image->data[idx+2];
						yuv[2]=image->data[idx+0];
						break;
					case RCT_B_R_GR:
					case RCT_B_R_GB:
					case RCT_B_RB_GB:
					case RCT_B_RB_GR:
					case RCT_B_RB_G2:
						yuv[0]=image->data[idx+2];
						yuv[1]=image->data[idx+0];
						yuv[2]=image->data[idx+1];
						break;
					case RCT_G_RG_BR:
					case RCT_G_RG_B2:
						yuv[0]=image->data[idx+1];
						yuv[1]=image->data[idx+0];
						yuv[2]=image->data[idx+2];
						break;
					case RCT_B_GB_RG:
					case RCT_B_GB_R2:
						yuv[0]=image->data[idx+2];
						yuv[1]=image->data[idx+1];
						yuv[2]=image->data[idx+0];
						break;
					case RCT_R_BR_GB:
					case RCT_R_B_G2:
					case RCT_R_BR_G2:
						yuv[0]=image->data[idx+0];
						yuv[1]=image->data[idx+2];
						yuv[2]=image->data[idx+1];
						break;
					case RCT_RCT1_120:
					case RCT_RCT1_201:
					case RCT_RCT1_012:
					case RCT_RCT1_210:
					case RCT_RCT1_102:
					case RCT_RCT1_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT1_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT1_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT1_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=LUMA_UPDATE_1(yuv[1], yuv[2]);
						break;
					case RCT_RCT2_120:
					case RCT_RCT2_201:
					case RCT_RCT2_012:
					case RCT_RCT2_210:
					case RCT_RCT2_102:
					case RCT_RCT2_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT2_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT2_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT2_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=LUMA_UPDATE_2(yuv[1], yuv[2]);
						break;
					case RCT_RCT3_120:
					case RCT_RCT3_201:
					case RCT_RCT3_012:
					case RCT_RCT3_210:
					case RCT_RCT3_102:
					case RCT_RCT3_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT3_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT3_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT3_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=LUMA_UPDATE_3(yuv[1], yuv[2]);
						break;
					case RCT_RCT4_120:
					case RCT_RCT4_201:
					case RCT_RCT4_012:
					case RCT_RCT4_210:
					case RCT_RCT4_102:
					case RCT_RCT4_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT4_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT4_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT4_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=LUMA_UPDATE_4(yuv[1], yuv[2]);
						break;
					case RCT_RCT5_120:
					case RCT_RCT5_201:
					case RCT_RCT5_012:
					case RCT_RCT5_210:
					case RCT_RCT5_102:
					case RCT_RCT5_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT5_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT5_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT5_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=LUMA_UPDATE_5(yuv[1], yuv[2]);
						break;
					case RCT_RCT6_120:
					case RCT_RCT6_201:
					case RCT_RCT6_012:
					case RCT_RCT6_210:
					case RCT_RCT6_102:
					case RCT_RCT6_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT6_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT6_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT6_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=LUMA_UPDATE_6(yuv[1], yuv[2]);
						break;
					case RCT_RCT7_120:
					case RCT_RCT7_201:
					case RCT_RCT7_012:
					case RCT_RCT7_210:
					case RCT_RCT7_102:
					case RCT_RCT7_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT7_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT7_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT7_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=LUMA_UPDATE_7(yuv[1], yuv[2]);
						break;
					case RCT_Pei09_120:
					case RCT_Pei09_201:
					case RCT_Pei09_012:
					case RCT_Pei09_210:
					case RCT_Pei09_102:
					case RCT_Pei09_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_Pei09_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_Pei09_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_Pei09_120][2]];
						yuv[1]-=(87*yuv[2]+169*yuv[0]+128)>>8;
						yuv[2]-=yuv[0];
						yuv[0]+=(86*yuv[2]+29*yuv[1]+128)>>8;
						break;
					case RCT_J2K_120:
					case RCT_J2K_201:
					case RCT_J2K_012:
						yuv[0]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_J2K_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_J2K_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[bestrct-RCT_J2K_120][2]];
						yuv[1]-=yuv[0];
						yuv[2]-=yuv[0];
						yuv[0]+=(yuv[1]+yuv[2])>>2;
						break;
					}
					kc=3;
				}
				if(kc<image->nch)
					yuv[kc]=image->data[idx+kc];
			}
#ifdef ENABLE_OLS
			++kols;
#endif
			//if(ky==10&&kx==10)//
			//	printf("");
			__m128i mNNN	=_mm_load_si128((__m128i*)NNN);
			__m128i mNN	=_mm_load_si128((__m128i*)NN);
			__m128i mNNE	=_mm_load_si128((__m128i*)NNE);
			__m128i mNW	=_mm_load_si128((__m128i*)NW);
			__m128i mN	=_mm_load_si128((__m128i*)N);
			__m128i mNE	=_mm_load_si128((__m128i*)NE);
			__m128i mWWW	=_mm_load_si128((__m128i*)WWW);
			__m128i mWW	=_mm_load_si128((__m128i*)WW);
			__m128i mW	=_mm_load_si128((__m128i*)W);
			__m128i t0=_mm_abs_epi16(_mm_sub_epi16(mW, mWW));
			__m128i t1=_mm_abs_epi16(_mm_sub_epi16(mN, mNN));
			t0=_mm_add_epi16(t0, _mm_abs_epi16(_mm_sub_epi16(mN, mNW)));
			t1=_mm_add_epi16(t1, _mm_abs_epi16(_mm_sub_epi16(mW, mNW)));
			t0=_mm_add_epi16(t0, _mm_abs_epi16(_mm_sub_epi16(mNE, mN)));
			t1=_mm_add_epi16(t1, _mm_abs_epi16(_mm_sub_epi16(mNE, mNNE)));
			t0=_mm_add_epi16(t0, _mm_shuffle_epi8(_mm_abs_epi16(mWWW), swap16));
			t1=_mm_add_epi16(t1, _mm_shuffle_epi8(_mm_abs_epi16(mNNN), swap16));
			t0=_mm_add_epi16(t0, _mm_shuffle_epi8(_mm_abs_epi16(mWW), swap16));
			t1=_mm_add_epi16(t1, _mm_shuffle_epi8(_mm_abs_epi16(mNN), swap16));
			t0=_mm_add_epi16(t0, _mm_slli_epi16(_mm_shuffle_epi8(_mm_abs_epi16(mW), swap16), 1));
			t1=_mm_add_epi16(t1, _mm_slli_epi16(_mm_shuffle_epi8(_mm_abs_epi16(mN), swap16), 1));
			ALIGN(16) short grads[16];
			_mm_store_si128((__m128i*)grads+0, t0);
			_mm_store_si128((__m128i*)grads+1, t1);
			grads[ 0]>>=depths3[0], grads[ 0]=FLOOR_LOG2(grads[ 0]+1);
			grads[ 2]>>=depths3[1], grads[ 2]=FLOOR_LOG2(grads[ 2]+1);
			grads[ 4]>>=depths3[2], grads[ 4]=FLOOR_LOG2(grads[ 4]+1);
			grads[ 8]>>=depths3[0], grads[ 8]=FLOOR_LOG2(grads[ 8]+1);
			grads[10]>>=depths3[1], grads[10]=FLOOR_LOG2(grads[10]+1);
			grads[12]>>=depths3[2], grads[12]=FLOOR_LOG2(grads[12]+1);
			int ctxidx[]=
			{
				cdfstride*(nctx*0+args->clevels*MINVAR(grads[ 8], 8)+MINVAR(grads[0], 8)),
				cdfstride*(nctx*1+args->clevels*MINVAR(grads[10], 8)+MINVAR(grads[2], 8)),
				cdfstride*(nctx*2+args->clevels*MINVAR(grads[12], 8)+MINVAR(grads[4], 8)),
			};
			__m128i vmin=_mm_min_epi16(mN, mW);
			__m128i vmax=_mm_max_epi16(mN, mW);
			__m128i mg=_mm_sub_epi16(_mm_add_epi16(mN, mW), mNW);
			mg=_mm_max_epi16(mg, vmin);
			mg=_mm_min_epi16(mg, vmax);
			ALIGN(16) short cgrads[8];
			_mm_store_si128((__m128i*)cgrads, mg);
#if defined ENABLE_OLS && defined ALLOW_AVX2
			ALIGN(16) int preds[4];
			vmin=_mm_min_epi16(vmin, mNE);
			vmax=_mm_max_epi16(vmax, mNE);
			vmin=_mm_slli_epi32(vmin, 16);//need low signed 16 bits -> 32 bit
			vmax=_mm_slli_epi32(vmax, 16);
			vmin=_mm_srai_epi32(vmin, 16);
			vmax=_mm_srai_epi32(vmax, 16);
			int j0=0, j1=0, j2=0;
			ols4_ctx0[j0++]=cgrads	[0];
		//	ols4_ctx0[j0++]=NNWW	[0];
			ols4_ctx0[j0++]=NNW	[0];
			ols4_ctx0[j0++]=NN	[0];
			ols4_ctx0[j0++]=NNE	[0];
			ols4_ctx0[j0++]=NNEE	[0];
			ols4_ctx0[j0++]=NNEEE	[0];
			ols4_ctx0[j0++]=NWWW	[0];
			ols4_ctx0[j0++]=NWW	[0];
			ols4_ctx0[j0++]=NW	[0];
			ols4_ctx0[j0++]=N	[0];
			ols4_ctx0[j0++]=NE	[0];
			ols4_ctx0[j0++]=NEE	[0];
			ols4_ctx0[j0++]=NEEE	[0];
			ols4_ctx0[j0++]=WWW	[0];
			ols4_ctx0[j0++]=WW	[0];
			ols4_ctx0[j0++]=W	[0];

			ols4_ctx1[j1++]=cgrads	[2];
		//	ols4_ctx1[j1++]=NNWW	[2];
			ols4_ctx1[j1++]=NNW	[2];
			ols4_ctx1[j1++]=NN	[2];
			ols4_ctx1[j1++]=NNE	[2];
			ols4_ctx1[j1++]=NNEE	[2];
			ols4_ctx1[j1++]=NWW	[2];
			ols4_ctx1[j1++]=NW	[2];
			ols4_ctx1[j1++]=NW	[0];
			ols4_ctx1[j1++]=N	[2];
			ols4_ctx1[j1++]=N	[0];
			ols4_ctx1[j1++]=NE	[2];
			ols4_ctx1[j1++]=NE	[0];
			ols4_ctx1[j1++]=NEE	[2];
			ols4_ctx1[j1++]=WW	[2];
			ols4_ctx1[j1++]=W	[2];
			ols4_ctx1[j1++]=W	[0];
		//	ols4_ctx1[j1++]=curr	[0];//X  unavailable

			ols4_ctx2[j2++]=cgrads	[4];
		//	ols4_ctx2[j2++]=NNWW	[4];
			ols4_ctx2[j2++]=NNW	[4];
			ols4_ctx2[j2++]=NN	[4];
			ols4_ctx2[j2++]=NNE	[4];
			ols4_ctx2[j2++]=NNEE	[4];
			ols4_ctx2[j2++]=NWW	[4];
			ols4_ctx2[j2++]=NW	[4];
			ols4_ctx2[j2++]=NW	[0];
			ols4_ctx2[j2++]=N	[4];
			ols4_ctx2[j2++]=N	[0];
			ols4_ctx2[j2++]=NE	[4];
			ols4_ctx2[j2++]=NE	[0];
			ols4_ctx2[j2++]=NEE	[4];
			ols4_ctx2[j2++]=WW	[4];
			ols4_ctx2[j2++]=W	[4];
			ols4_ctx2[j2++]=W	[0];
		//	ols4_ctx2[j2++]=curr	[0];//X  unavailable
			{
				double *pp0=args->ols4[0].params;
				double *pp1=args->ols4[1].params;
				double *pp2=args->ols4[2].params;
				double *pnb0=args->ols4[0].nb;
				double *pnb1=args->ols4[1].nb;
				double *pnb2=args->ols4[2].nb;

				__m128i mnb0=_mm_load_si128((__m128i*)ols4_ctx0+0);
				__m128i mnb1=_mm_load_si128((__m128i*)ols4_ctx1+0);
				__m128i mnb2=_mm_load_si128((__m128i*)ols4_ctx2+0);
				__m256d dnb0=_mm256_cvtepi32_pd(mnb0);
				__m256d dnb1=_mm256_cvtepi32_pd(mnb1);
				__m256d dnb2=_mm256_cvtepi32_pd(mnb2);
				_mm256_store_pd(pnb0+0*4, dnb0);
				_mm256_store_pd(pnb1+0*4, dnb1);
				_mm256_store_pd(pnb2+0*4, dnb2);
				__m256d sum0=_mm256_mul_pd(dnb0, _mm256_load_pd(pp0+0*4));
				__m256d sum1=_mm256_mul_pd(dnb1, _mm256_load_pd(pp1+0*4));
				__m256d sum2=_mm256_mul_pd(dnb2, _mm256_load_pd(pp2+0*4));

				mnb0=_mm_load_si128((__m128i*)ols4_ctx0+1);
				mnb1=_mm_load_si128((__m128i*)ols4_ctx1+1);
				mnb2=_mm_load_si128((__m128i*)ols4_ctx2+1);
				dnb0=_mm256_cvtepi32_pd(mnb0);
				dnb1=_mm256_cvtepi32_pd(mnb1);
				dnb2=_mm256_cvtepi32_pd(mnb2);
				_mm256_store_pd(pnb0+1*4, dnb0);
				_mm256_store_pd(pnb1+1*4, dnb1);
				_mm256_store_pd(pnb2+1*4, dnb2);
				sum0=_mm256_add_pd(sum0, _mm256_mul_pd(dnb0, _mm256_load_pd(pp0+1*4)));
				sum1=_mm256_add_pd(sum1, _mm256_mul_pd(dnb1, _mm256_load_pd(pp1+1*4)));
				sum2=_mm256_add_pd(sum2, _mm256_mul_pd(dnb2, _mm256_load_pd(pp2+1*4)));

				mnb0=_mm_load_si128((__m128i*)ols4_ctx0+2);
				mnb1=_mm_load_si128((__m128i*)ols4_ctx1+2);
				mnb2=_mm_load_si128((__m128i*)ols4_ctx2+2);
				dnb0=_mm256_cvtepi32_pd(mnb0);
				dnb1=_mm256_cvtepi32_pd(mnb1);
				dnb2=_mm256_cvtepi32_pd(mnb2);
				_mm256_store_pd(pnb0+2*4, dnb0);
				_mm256_store_pd(pnb1+2*4, dnb1);
				_mm256_store_pd(pnb2+2*4, dnb2);
				sum0=_mm256_add_pd(sum0, _mm256_mul_pd(dnb0, _mm256_load_pd(pp0+2*4)));
				sum1=_mm256_add_pd(sum1, _mm256_mul_pd(dnb1, _mm256_load_pd(pp1+2*4)));
				sum2=_mm256_add_pd(sum2, _mm256_mul_pd(dnb2, _mm256_load_pd(pp2+2*4)));

				mnb0=_mm_load_si128((__m128i*)ols4_ctx0+3);
				mnb1=_mm_load_si128((__m128i*)ols4_ctx1+3);
				mnb2=_mm_load_si128((__m128i*)ols4_ctx2+3);
				dnb0=_mm256_cvtepi32_pd(mnb0);
				dnb1=_mm256_cvtepi32_pd(mnb1);
				dnb2=_mm256_cvtepi32_pd(mnb2);
				_mm256_store_pd(pnb0+3*4, dnb0);
				_mm256_store_pd(pnb1+3*4, dnb1);
				_mm256_store_pd(pnb2+3*4, dnb2);
				sum0=_mm256_add_pd(sum0, _mm256_mul_pd(dnb0, _mm256_load_pd(pp0+3*4)));
				sum1=_mm256_add_pd(sum1, _mm256_mul_pd(dnb1, _mm256_load_pd(pp1+3*4)));
				sum2=_mm256_add_pd(sum2, _mm256_mul_pd(dnb2, _mm256_load_pd(pp2+3*4)));

				ALIGN(32) double sums[12];
				_mm256_store_pd(sums+0*4, sum0);
				_mm256_store_pd(sums+1*4, sum1);
				_mm256_store_pd(sums+2*4, sum2);

				sums[ 0]+=sums[ 1]+sums[ 2]+sums[ 3];
				sums[ 4]+=sums[ 5]+sums[ 6]+sums[ 7];
				sums[ 8]+=sums[ 9]+sums[10]+sums[11];
#if OLS_NPARAMS&15
				for(int k=16;k<OLS_NPARAMS;++k)
				{
					pnb0[k]=(double)ols4_ctx0[k];
					pnb1[k]=(double)ols4_ctx1[k];
					pnb2[k]=(double)ols4_ctx2[k];
					sums[ 0]+=pp0[k]*pnb0[k];
					sums[ 4]+=pp1[k]*pnb1[k];
					sums[ 8]+=pp2[k]*pnb2[k];
				}
#endif
				sums[1]=sums[4];
				sums[2]=sums[8];

				__m128i mp=_mm256_cvtpd_epi32(_mm256_load_pd(sums));//rounded
				mp=_mm_min_epi32(mp, vmax);
				mp=_mm_max_epi32(mp, vmin);
				_mm_store_si128((__m128i*)preds, mp);
			}
#endif
			for(int kc=0;kc<image->nch;++kc)
			{
				int ch=combination[kc], depth=depths[ch], half=halfs2[kc];
				int kc2=kc<<1;
				int offset=(yuv[combination[kc+4]]+yuv[combination[kc+8]])>>och_info[ch][II_HSHIFT];
				int pred=0, error, sym;
				//int
				//	vx=(abs(W[kc2]-WW[kc2])+abs(N[kc2]-NW[kc2])+abs(NE[kc2]-N  [kc2])+abs(WWW[kc2+1])+abs(WW[kc2+1])+abs(W[kc2+1])*2)<<8>>depth,
				//	vy=(abs(N[kc2]-NN[kc2])+abs(W[kc2]-NW[kc2])+abs(NE[kc2]-NNE[kc2])+abs(NNN[kc2+1])+abs(NN[kc2+1])+abs(N[kc2+1])*2)<<8>>depth;
				//int qeN=FLOOR_LOG2(vy+1);
				//int qeW=FLOOR_LOG2(vx+1);
				//int cidx=cdfstride*(nctx*kc+args->clevels*MINVAR(qeN, 8)+MINVAR(qeW, 8));
				//int *curr_hist=args->hist+cidx;
				int *curr_hist=args->hist+ctxidx[kc];
				//unsigned *curr_CDF=args->stats+cidx;

				int den=curr_hist[args->tlevels]+args->tlevels, cdf=0, freq;
#ifdef ENABLE_OLS
#ifdef ALLOW_AVX2
				pred=preds[kc];
#else
				int j=0, cmin, cmax;
				{
					ALIGN(16) int arr[4];
					__m128i va=_mm_set_epi32(0, 0, 0, N[kc2]);
					__m128i vb=_mm_set_epi32(0, 0, 0, W[kc2]);
					__m128i vc=_mm_set_epi32(0, 0, 0, N[kc2]+W[kc2]-NW[kc2]);
					__m128i _vmin=_mm_min_epi32(va, vb);
					__m128i _vmax=_mm_max_epi32(va, vb);
					vc=_mm_max_epi32(vc, _vmin);
					vc=_mm_min_epi32(vc, _vmax);
					_mm_store_si128((__m128i*)arr, vc);
					pred=arr[0];

					vc=_mm_set_epi32(0, 0, 0, NE[kc2]);
					_vmin=_mm_min_epi32(_vmin, vc);
					_vmax=_mm_max_epi32(_vmax, vc);
					_mm_store_si128((__m128i*)arr, _vmin);
					cmin=arr[0];
					_mm_store_si128((__m128i*)arr, _vmax);
					cmax=arr[0];
				}
				switch(kc)
				{
				case 0:
					ols4_ctx[j++]=pred;
					ols4_ctx[j++]=NNWW	[kc2];
					ols4_ctx[j++]=NNW	[kc2];
					ols4_ctx[j++]=NN	[kc2];
					ols4_ctx[j++]=NNE	[kc2];
					ols4_ctx[j++]=NNEE	[kc2];
					ols4_ctx[j++]=NNEEE	[kc2];
					ols4_ctx[j++]=NWWW	[kc2];
					ols4_ctx[j++]=NWW	[kc2];
					ols4_ctx[j++]=NW	[kc2];
					ols4_ctx[j++]=N		[kc2];
					ols4_ctx[j++]=NE	[kc2];
					ols4_ctx[j++]=NEE	[kc2];
					ols4_ctx[j++]=NEEE	[kc2];
					ols4_ctx[j++]=WWW	[kc2];
					ols4_ctx[j++]=WW	[kc2];
					ols4_ctx[j++]=W		[kc2];
					pred=ols4_predict(args->ols4+0, ols4_ctx, cmin, cmax);
					break;
				case 1:
					ols4_ctx[j++]=pred;
					ols4_ctx[j++]=NNWW	[kc2];
					ols4_ctx[j++]=NNW	[kc2];
					ols4_ctx[j++]=NN	[kc2];
					ols4_ctx[j++]=NNE	[kc2];
					ols4_ctx[j++]=NNEE	[kc2];
					ols4_ctx[j++]=NWW	[kc2];
					ols4_ctx[j++]=NW	[kc2];
					ols4_ctx[j++]=NW	[kc2-2];
					ols4_ctx[j++]=N		[kc2];
					ols4_ctx[j++]=N		[kc2-2];
					ols4_ctx[j++]=NE	[kc2];
					ols4_ctx[j++]=NE	[kc2-2];
					ols4_ctx[j++]=NEE	[kc2];
					ols4_ctx[j++]=WW	[kc2];
					ols4_ctx[j++]=W		[kc2];
					ols4_ctx[j++]=W		[kc2-2];
				//	ols4_ctx[j++]=curr	[kc2-2];
					pred=ols4_predict(args->ols4+1, ols4_ctx, cmin, cmax);
					break;
				case 2:
					ols4_ctx[j++]=pred;
					ols4_ctx[j++]=NNWW	[kc2];
					ols4_ctx[j++]=NNW	[kc2];
					ols4_ctx[j++]=NN	[kc2];
					ols4_ctx[j++]=NNE	[kc2];
					ols4_ctx[j++]=NNEE	[kc2];
					ols4_ctx[j++]=NWW	[kc2];
					ols4_ctx[j++]=NW	[kc2];
					ols4_ctx[j++]=NW	[kc2-4];
					ols4_ctx[j++]=N		[kc2];
					ols4_ctx[j++]=N		[kc2-4];
					ols4_ctx[j++]=NE	[kc2];
					ols4_ctx[j++]=NE	[kc2-4];
					ols4_ctx[j++]=NEE	[kc2];
					ols4_ctx[j++]=WW	[kc2];
					ols4_ctx[j++]=W		[kc2];
					ols4_ctx[j++]=W		[kc2-4];
				//	ols4_ctx[j++]=curr	[kc2-4];
					pred=ols4_predict(args->ols4+2, ols4_ctx, cmin, cmax);
					break;
				}
#endif
#else
				pred=cgrads[kc2];
				//MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
#endif
#ifdef ENABLE_SSE
				int pred0=pred;
				long long *ssecell=0, ssesum=0;
				int ssecount=0;
				int ssecorr=0;
				if(kc<3)
				{
					int sseidx1=(N[kc2]-NW[kc2]+half)<<SSEBITS>>depth;
					int sseidx2=(W[kc2]-NW[kc2]+half)<<SSEBITS>>depth;
					CLAMP2(sseidx1, 0, (1<<SSEBITS)-1);
					CLAMP2(sseidx2, 0, (1<<SSEBITS)-1);
					ssecell=&args->sse[kc][sseidx1][sseidx2];
					ssesum=*ssecell>>12;
					ssecount=(int)(*ssecell&0xFFF);
					ssecorr=(int)(ssesum/(ssecount+1));
					pred+=ssecorr;
				}
#endif
				pred+=offset;
				CLAMP2_32(pred, pred, -half, half-1);

				//if(ky==13&&kx==623)//
				//if(ky==80&&kx==510&&kc==2)//
				//if(ky==16&&kx==640&&kc==1)//
				//if(ky==0&&kx==258&&kc==1)//
				//if(ky==17&&kx==458&&kc==2)//
				//if(ky==255&&kx==767&&kc==0)//
				//if(ky==0&&kx==125&&kc==0)//
				//if(ky==0&&kx==102&&kc==1)//
				//if(ky==1&&kx==111&&kc==1)//
				//if(ky==22&&kx==242&&kc==3)//
				//if(ky==0&&kx==4&&kc==0)//
				//	printf("");

				if(args->fwd)
				{
					curr[kc2+0]=yuv[kc];
					curr[kc2+1]=error=yuv[kc]-pred;
					{
						int upred=half-abs(pred), aval=abs(error);
						if(aval<=upred)
						{
							sym=error;
#ifdef ENABLE_SSE
							{
								int negmask=-((ssecorr<0)&(sym!=-half));//sign is flipped if SSE correction was negative, to skew the histogram
								sym^=negmask;
								sym-=negmask;
							}
#endif
							sym=sym<<1^sym>>31;//pack sign
						}
						else
							sym=upred+aval;//error sign is known
					}
					quantize_pixel(sym, &token, &bypass, &nbits);

					freq=curr_hist[token]+1;
					for(int k=0;k<token;++k)
						cdf+=curr_hist[k]+1;
					ac3_enc_update_NPOT(&ec, cdf, freq, den);
					//ac3_enc(&ec, token, curr_CDF);
					if(nbits)
						ac3_enc_bypass(&ec, bypass, nbits);//up to 16 bits
				}
				else
				{
					int code=ac3_dec_getcdf_NPOT(&ec, den);
					token=0;
					for(;;)
					{
						int cdf2;

						freq=curr_hist[token]+1;
						cdf2=cdf+freq;
						if(cdf2>code)
							break;
						cdf=cdf2;
						++token;
					}
					//if(token>=args->tlevels)//
					//	LOG_ERROR("");
					ac3_dec_update_NPOT(&ec, cdf, freq, den);
					//token=ac3_dec(&ec, curr_CDF, args->tlevels);//try ac_dec_packedsign()
					sym=token;
					if(sym>=(1<<CONFIG_EXP))
					{
						int lsb, msb;

						sym-=1<<CONFIG_EXP;
						lsb=sym&((1<<CONFIG_LSB)-1);
						sym>>=CONFIG_LSB;
						msb=sym&((1<<CONFIG_MSB)-1);
						sym>>=CONFIG_MSB;
						nbits=sym+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB);
						bypass=ac3_dec_bypass(&ec, nbits);
						sym=1;
						sym<<=CONFIG_MSB;
						sym|=msb;
						sym<<=nbits;
						sym|=bypass;
						sym<<=CONFIG_LSB;
						sym|=lsb;
					}
					{
						int upred=half-abs(pred), negmask=0;
						if(sym<=(upred<<1))
						{
							error=sym>>1^-(sym&1);
#ifdef ENABLE_SSE
							negmask=-((ssecorr<0)&(error!=-half));
#endif
						}
						else
						{
							error=sym-upred;
							negmask=-(pred>0);
						}
						error^=negmask;
						error-=negmask;
					}
					curr[kc2+0]=yuv[kc]=error+pred;
					curr[kc2+1]=error;
				}
				++curr_hist[token];
				++curr_hist[args->tlevels];
				
				if(curr_hist[args->tlevels]>=6144)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist[ks]>>=1;
					curr_hist[args->tlevels]=sum;
				}
#if 0
				//if((kx&(CDFSTRIDE-1))==CDFSTRIDE-1)
				if(curr_hist[args->tlevels]>=args->tlevels&&!(curr_hist[args->tlevels]&(CDFSTRIDE-1)))
					update_CDF(curr_hist, curr_CDF, args->tlevels);
				if(curr_hist[args->tlevels]>=6144)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist[ks]>>=1;
					curr_hist[args->tlevels]=sum;
				}
#endif
#ifdef ENABLE_SSE
				if(ssecell)
				{
					++ssecount;
					ssesum+=(long long)curr[kc2]-pred;
					if(ssecount>=256)
					{
						ssecount>>=1;
						ssesum>>=1;
					}
					*ssecell=ssesum<<12|ssecount;
				}
#endif
				curr[kc2]-=offset;
#if defined ENABLE_OLS && !defined ALLOW_AVX2
				if(kc<3&&(kx<5||kols>=kols2))
					ols4_update(args->ols4+kc, curr[kc2], lrs[kc]);
#endif
			}
#ifdef ENABLE_OLS
			if(kx<5||kols>=kols2)
			{
				double *pnb0=args->ols4[0].nb;
				double *pnb1=args->ols4[1].nb;
				double *pnb2=args->ols4[2].nb;
				double *cov0=args->ols4[0].cov;
				double *cov1=args->ols4[1].cov;
				double *cov2=args->ols4[2].cov;
				double lval0=curr[0]*lrs[0], lr_comp0=1-lrs[0];
				double lval1=curr[2]*lrs[1], lr_comp1=1-lrs[1];
				double lval2=curr[4]*lrs[2], lr_comp2=1-lrs[2];
				__m256d mlr0=_mm256_set1_pd(lrs[0]);
				__m256d mlr1=_mm256_set1_pd(lrs[1]);
				__m256d mlr2=_mm256_set1_pd(lrs[2]);
				for(int ky2=0, midx=0;ky2<OLS4_NPARAMS;++ky2)
				{
					__m256d vy0=_mm256_set1_pd(pnb0[ky2]);
					__m256d vy1=_mm256_set1_pd(pnb1[ky2]);
					__m256d vy2=_mm256_set1_pd(pnb2[ky2]);
					int kx2=0;
					for(;kx2<OLS4_NPARAMS-3;kx2+=4, midx+=4)
					{
						__m256d mcov0=_mm256_load_pd(cov0+midx);
						__m256d mcov1=_mm256_load_pd(cov1+midx);
						__m256d mcov2=_mm256_load_pd(cov2+midx);
						__m256d vx0=_mm256_load_pd(pnb0+kx2);
						__m256d vx1=_mm256_load_pd(pnb1+kx2);
						__m256d vx2=_mm256_load_pd(pnb2+kx2);
						vx0=_mm256_mul_pd(vx0, vy0);
						vx1=_mm256_mul_pd(vx1, vy1);
						vx2=_mm256_mul_pd(vx2, vy2);
						vx0=_mm256_sub_pd(vx0, mcov0);
						vx1=_mm256_sub_pd(vx1, mcov1);
						vx2=_mm256_sub_pd(vx2, mcov2);
						vx0=_mm256_mul_pd(vx0, mlr0);
						vx1=_mm256_mul_pd(vx1, mlr1);
						vx2=_mm256_mul_pd(vx2, mlr2);
						vx0=_mm256_add_pd(vx0, mcov0);
						vx1=_mm256_add_pd(vx1, mcov1);
						vx2=_mm256_add_pd(vx2, mcov2);
						_mm256_store_pd(cov0+midx, vx0);
						_mm256_store_pd(cov1+midx, vx1);
						_mm256_store_pd(cov2+midx, vx2);
					}
#if OLS4_NPARAMS&15
					for(;kx2<OLS4_NPARAMS;++kx2, ++midx)
					{
						cov0[midx]+=(pnb0[kx2]*pnb0[ky2]-cov0[midx])*lrs[0];
						cov1[midx]+=(pnb1[kx2]*pnb1[ky2]-cov1[midx])*lrs[1];
						cov2[midx]+=(pnb2[kx2]*pnb2[ky2]-cov2[midx])*lrs[2];
					}
#endif
				}
				double *vec0=args->ols4[0].vec;
				double *vec1=args->ols4[1].vec;
				double *vec2=args->ols4[2].vec;
				{
					__m256d mlr_comp0=_mm256_set1_pd(lr_comp0);
					__m256d mlr_comp1=_mm256_set1_pd(lr_comp1);
					__m256d mlr_comp2=_mm256_set1_pd(lr_comp2);
					__m256d mval0=_mm256_set1_pd(lval0);
					__m256d mval1=_mm256_set1_pd(lval1);
					__m256d mval2=_mm256_set1_pd(lval2);
					int k=0;
					for(;k<OLS4_NPARAMS-3;k+=4)
					{
						__m256d mvec0=_mm256_load_pd(vec0+k);
						__m256d mvec1=_mm256_load_pd(vec1+k);
						__m256d mvec2=_mm256_load_pd(vec2+k);
						__m256d mtmp0=_mm256_load_pd(pnb0+k);
						__m256d mtmp1=_mm256_load_pd(pnb1+k);
						__m256d mtmp2=_mm256_load_pd(pnb2+k);
						mvec0=_mm256_mul_pd(mvec0, mlr_comp0);
						mvec1=_mm256_mul_pd(mvec1, mlr_comp1);
						mvec2=_mm256_mul_pd(mvec2, mlr_comp2);
						mtmp0=_mm256_mul_pd(mtmp0, mval0);
						mtmp1=_mm256_mul_pd(mtmp1, mval1);
						mtmp2=_mm256_mul_pd(mtmp2, mval2);
						mvec0=_mm256_add_pd(mvec0, mtmp0);
						mvec1=_mm256_add_pd(mvec1, mtmp1);
						mvec2=_mm256_add_pd(mvec2, mtmp2);
						_mm256_store_pd(vec0+k, mvec0);
						_mm256_store_pd(vec1+k, mvec1);
						_mm256_store_pd(vec2+k, mvec2);
					}
#if OLS4_NPARAMS&15
					for(;k<OLS4_NPARAMS;++k)
					{
						vec0[k]=lval0*pnb0[k]+lr_comp0*vec0[k];
						vec1[k]=lval1*pnb1[k]+lr_comp1*vec1[k];
						vec2[k]=lval2*pnb2[k]+lr_comp2*vec2[k];
					}
#endif
				}
				double *chol0=args->ols4[0].cholesky;
				double *chol1=args->ols4[1].cholesky;
				double *chol2=args->ols4[2].cholesky;
				double *params0=args->ols4[0].params;
				double *params1=args->ols4[1].params;
				double *params2=args->ols4[2].params;
				int success0=1;
				int success1=1;
				int success2=1;
				memcpy(chol0, cov0, sizeof(double[OLS4_NPARAMS*OLS4_NPARAMS]));
				memcpy(chol1, cov1, sizeof(double[OLS4_NPARAMS*OLS4_NPARAMS]));
				memcpy(chol2, cov2, sizeof(double[OLS4_NPARAMS*OLS4_NPARAMS]));
				for(int k=0;k<OLS4_NPARAMS*OLS4_NPARAMS;k+=OLS4_NPARAMS+1)
				{
					chol0[k]+=0.0075;
					chol1[k]+=0.0075;
					chol2[k]+=0.0075;
				}
				double sum;
				for(int i=0;i<OLS4_NPARAMS;++i)
				{
					for(int j=0;j<i;++j)
					{
						sum=chol0[i*OLS4_NPARAMS+j];
						for(int k=0;k<j;++k)
							sum-=chol0[i*OLS4_NPARAMS+k]*chol0[j*OLS4_NPARAMS+k];
						chol0[i*OLS4_NPARAMS+j]=sum/chol0[j*OLS4_NPARAMS+j];
					}
					sum=chol0[i*OLS4_NPARAMS+i];
					for(int k=0;k<i;++k)
						sum-=chol0[i*OLS4_NPARAMS+k]*chol0[i*OLS4_NPARAMS+k];
					if(sum<=1e-8)
					{
						success0=0;
						break;
					}
					chol0[i*OLS4_NPARAMS+i]=sqrt(sum);
				}
				for(int i=0;i<OLS4_NPARAMS;++i)
				{
					for(int j=0;j<i;++j)
					{
						sum=chol1[i*OLS4_NPARAMS+j];
						for(int k=0;k<j;++k)
							sum-=chol1[i*OLS4_NPARAMS+k]*chol1[j*OLS4_NPARAMS+k];
						chol1[i*OLS4_NPARAMS+j]=sum/chol1[j*OLS4_NPARAMS+j];
					}
					sum=chol1[i*OLS4_NPARAMS+i];
					for(int k=0;k<i;++k)
						sum-=chol1[i*OLS4_NPARAMS+k]*chol1[i*OLS4_NPARAMS+k];
					if(sum<=1e-8)
					{
						success1=0;
						break;
					}
					chol1[i*OLS4_NPARAMS+i]=sqrt(sum);
				}
				for(int i=0;i<OLS4_NPARAMS;++i)
				{
					for(int j=0;j<i;++j)
					{
						sum=chol2[i*OLS4_NPARAMS+j];
						for(int k=0;k<j;++k)
							sum-=chol2[i*OLS4_NPARAMS+k]*chol2[j*OLS4_NPARAMS+k];
						chol2[i*OLS4_NPARAMS+j]=sum/chol2[j*OLS4_NPARAMS+j];
					}
					sum=chol2[i*OLS4_NPARAMS+i];
					for(int k=0;k<i;++k)
						sum-=chol2[i*OLS4_NPARAMS+k]*chol2[i*OLS4_NPARAMS+k];
					if(sum<=1e-8)
					{
						success2=0;
						break;
					}
					chol2[i*OLS4_NPARAMS+i]=sqrt(sum);
				}
				if(success0)
				{
					for(int i=0;i<OLS4_NPARAMS;++i)
					{
						sum=vec0[i];
						for(int j=0;j<i;++j)
							sum-=chol0[i*OLS4_NPARAMS+j]*params0[j];
						params0[i]=sum/chol0[i*OLS4_NPARAMS+i];
					}
					for(int i=OLS4_NPARAMS-1;i>=0;--i)
					{
						sum=params0[i];
						for(int j=i+1;j<OLS4_NPARAMS;++j)
							sum-=chol0[j*OLS4_NPARAMS+i]*params0[j];
						params0[i]=sum/chol0[i*OLS4_NPARAMS+i];
					}
				}
				if(success1)
				{
					for(int i=0;i<OLS4_NPARAMS;++i)
					{
						sum=vec1[i];
						for(int j=0;j<i;++j)
							sum-=chol1[i*OLS4_NPARAMS+j]*params1[j];
						params1[i]=sum/chol1[i*OLS4_NPARAMS+i];
					}
					for(int i=OLS4_NPARAMS-1;i>=0;--i)
					{
						sum=params1[i];
						for(int j=i+1;j<OLS4_NPARAMS;++j)
							sum-=chol1[j*OLS4_NPARAMS+i]*params1[j];
						params1[i]=sum/chol1[i*OLS4_NPARAMS+i];
					}
				}
				if(success2)
				{
					for(int i=0;i<OLS4_NPARAMS;++i)
					{
						sum=vec2[i];
						for(int j=0;j<i;++j)
							sum-=chol2[i*OLS4_NPARAMS+j]*params2[j];
						params2[i]=sum/chol2[i*OLS4_NPARAMS+i];
					}
					for(int i=OLS4_NPARAMS-1;i>=0;--i)
					{
						sum=params2[i];
						for(int j=i+1;j<OLS4_NPARAMS;++j)
							sum-=chol2[j*OLS4_NPARAMS+i]*params2[j];
						params2[i]=sum/chol2[i*OLS4_NPARAMS+i];
					}
				}
			}
			if(kols>=kols2)
				kols2+=OLS_STRIDE;
#endif
			if(!args->fwd)
			{
				int kc=0;
				if(image->nch>=3)
				{
					switch(bestrct)
					{
					case RCT_R_G_B:
					case RCT_R_G_BG:
					case RCT_R_G_BR:
					case RCT_R_GR_BR:
					case RCT_R_GR_BG:
					case RCT_R_G_B2:
					case RCT_R_GR_B2:
						image->data[idx+0]=yuv[0];
						image->data[idx+1]=yuv[1];
						image->data[idx+2]=yuv[2];
						break;
					case RCT_G_B_RG:
					case RCT_G_B_RB:
					case RCT_G_BG_RG:
					case RCT_G_BG_RB:
					case RCT_G_B_R2:
					case RCT_G_BG_R2:
						image->data[idx+1]=yuv[0];
						image->data[idx+2]=yuv[1];
						image->data[idx+0]=yuv[2];
						break;
					case RCT_B_R_GR:
					case RCT_B_R_GB:
					case RCT_B_RB_GB:
					case RCT_B_RB_GR:
					case RCT_B_RB_G2:
						image->data[idx+2]=yuv[0];
						image->data[idx+0]=yuv[1];
						image->data[idx+1]=yuv[2];
						break;
					case RCT_G_RG_BR:
					case RCT_G_RG_B2:
						image->data[idx+1]=yuv[0];
						image->data[idx+0]=yuv[1];
						image->data[idx+2]=yuv[2];
						break;
					case RCT_B_GB_RG:
					case RCT_B_GB_R2:
						image->data[idx+2]=yuv[0];
						image->data[idx+1]=yuv[1];
						image->data[idx+0]=yuv[2];
						break;
					case RCT_R_BR_GB:
					case RCT_R_B_G2:
					case RCT_R_BR_G2:
						image->data[idx+0]=yuv[0];
						image->data[idx+2]=yuv[1];
						image->data[idx+1]=yuv[2];
						break;
					case RCT_RCT1_120:
					case RCT_RCT1_201:
					case RCT_RCT1_012:
					case RCT_RCT1_210:
					case RCT_RCT1_102:
					case RCT_RCT1_021:
						yuv[0]-=LUMA_UPDATE_1(yuv[1], yuv[2]);
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT1_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT1_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT1_120][2]]=yuv[2];
						break;
					case RCT_RCT2_120:
					case RCT_RCT2_201:
					case RCT_RCT2_012:
					case RCT_RCT2_210:
					case RCT_RCT2_102:
					case RCT_RCT2_021:
						yuv[0]-=LUMA_UPDATE_2(yuv[1], yuv[2]);
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT2_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT2_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT2_120][2]]=yuv[2];
						break;
					case RCT_RCT3_120:
					case RCT_RCT3_201:
					case RCT_RCT3_012:
					case RCT_RCT3_210:
					case RCT_RCT3_102:
					case RCT_RCT3_021:
						yuv[0]-=LUMA_UPDATE_3(yuv[1], yuv[2]);
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT3_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT3_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT3_120][2]]=yuv[2];
						break;
					case RCT_RCT4_120:
					case RCT_RCT4_201:
					case RCT_RCT4_012:
					case RCT_RCT4_210:
					case RCT_RCT4_102:
					case RCT_RCT4_021:
						yuv[0]-=LUMA_UPDATE_4(yuv[1], yuv[2]);
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT4_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT4_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT4_120][2]]=yuv[2];
						break;
					case RCT_RCT5_120:
					case RCT_RCT5_201:
					case RCT_RCT5_012:
					case RCT_RCT5_210:
					case RCT_RCT5_102:
					case RCT_RCT5_021:
						yuv[0]-=LUMA_UPDATE_5(yuv[1], yuv[2]);
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT5_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT5_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT5_120][2]]=yuv[2];
						break;
					case RCT_RCT6_120:
					case RCT_RCT6_201:
					case RCT_RCT6_012:
					case RCT_RCT6_210:
					case RCT_RCT6_102:
					case RCT_RCT6_021:
						yuv[0]-=LUMA_UPDATE_6(yuv[1], yuv[2]);
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT6_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT6_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT6_120][2]]=yuv[2];
						break;
					case RCT_RCT7_120:
					case RCT_RCT7_201:
					case RCT_RCT7_012:
					case RCT_RCT7_210:
					case RCT_RCT7_102:
					case RCT_RCT7_021:
						yuv[0]-=LUMA_UPDATE_7(yuv[1], yuv[2]);
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT7_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT7_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_RCT7_120][2]]=yuv[2];
						break;
					case RCT_Pei09_120:
					case RCT_Pei09_201:
					case RCT_Pei09_012:
					case RCT_Pei09_210:
					case RCT_Pei09_102:
					case RCT_Pei09_021:
						yuv[0]-=(86*yuv[2]+29*yuv[1]+128)>>8;
						yuv[2]+=yuv[0];
						yuv[1]+=(87*yuv[2]+169*yuv[0]+128)>>8;
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_Pei09_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_Pei09_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_Pei09_120][2]]=yuv[2];
						break;
					case RCT_J2K_120:
					case RCT_J2K_201:
					case RCT_J2K_012:
						yuv[0]-=(yuv[1]+yuv[2])>>2;
						yuv[2]+=yuv[0];
						yuv[1]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_J2K_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_J2K_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[bestrct-RCT_J2K_120][2]]=yuv[2];
						break;
					}
					kc=3;
				}
				if(kc<image->nch)
					image->data[idx+kc]=yuv[kc];
#ifdef ENABLE_GUIDE
				if(memcmp(image->data+idx, guide->data+idx, sizeof(short)*image->nch))
				{
					short orig[4]={0};
					memcpy(orig, guide->data+idx, image->nch*sizeof(short));
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	if(args->fwd)
		ac3_enc_flush(&ec);
}
int f33_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	int xblocks, yblocks, nblocks, ncores, nthreads, coffset;
	ptrdiff_t start=0;
	ptrdiff_t memusage=0;
	ptrdiff_t argssize;
	ThreadArgs *args;
	int histindices[OCH_COUNT+1]={0}, histsize=0;
	int tlevels=0, clevels=0, statssize=0;
	double bestsize=0;
	
	ncores=query_cpu_cores();
	xblocks=(image->iw+BLOCKSIZE-1)/BLOCKSIZE;
	yblocks=(image->ih+BLOCKSIZE-1)/BLOCKSIZE;
	nblocks=xblocks*yblocks;
	nthreads=MINVAR(nblocks, ncores);
	coffset=sizeof(int)*nblocks;
	argssize=nthreads*sizeof(ThreadArgs);
	args=(ThreadArgs*)malloc(argssize);
	if(fwd)
	{
#ifdef ENABLE_GUIDE
		guide=image;
#endif
		histsize=0;
		for(int kc=0;kc<OCH_COUNT;++kc)
		{
			int depth=image->depth+och_info[kc][II_INFLATION];
			UPDATE_MIN(depth, 16);
			int nlevels=1<<depth;
			histindices[kc]=histsize;
			histsize+=nlevels;
		}
		histindices[OCH_COUNT]=histsize;
		start=array_append(data, 0, 1, coffset, 1, 0, 0);
	}
	else//integrity check
	{
		start=coffset;
		for(int kt=0;kt<nblocks;++kt)
		{
			int size=0;
			memcpy(&size, cbuf+sizeof(int)*kt, sizeof(int));
			start+=size;
		}
		if(start!=(ptrdiff_t)clen)
			LOG_ERROR("Corrupt file");
		start=coffset;
	}
	{
		int nlevels=2<<image->depth;//chroma-inflated
		int token=0, bypass=0, nbits=0;

		quantize_pixel(nlevels, &token, &bypass, &nbits);
		tlevels=token+1;
		clevels=9;
		statssize=clevels*clevels*(tlevels+1)*image->nch*(int)sizeof(int);
		if(fwd)
		{
			UPDATE_MAX(histsize, statssize);
		}
		else
			histsize=statssize;
	}

	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memusage+=argssize;
	memset(args, 0, argssize);
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=src;
		arg->dst=dst;
		arg->bufsize=sizeof(short[4*OCH_COUNT])*(image->iw+16LL);//4 padded rows * OCH_COUNT * {pixels, wg_errors}
		arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));

		arg->histsize=histsize*sizeof(int);
		arg->hist=(int*)malloc(arg->histsize);

		arg->statssize=statssize;
		//arg->stats=(unsigned*)malloc(statssize);
#ifdef ENABLE_OLS
		for(int k=0;k<image->nch&&k<3;++k)
			ols4_init(arg->ols4+k, OLS4_NPARAMS);
#endif
		if(!arg->pixels||!arg->hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=arg->bufsize;
		memusage+=arg->histsize;
		//memusage+=arg->statssize;
		
		arg->tlevels=tlevels;
		arg->clevels=clevels;
		memcpy(arg->histindices, histindices, sizeof(arg->histindices));
		arg->fwd=fwd;
#ifdef DISABLE_MT
		arg->loud=loud&&nblocks<MAXPRINTEDBLOCKS;
#else
		arg->loud=0;
#endif
	}
	for(int kt=0;kt<nblocks;kt+=nthreads)
	{
		int nthreads2=MINVAR(kt+nthreads, nblocks)-kt;
		for(int kt2=0;kt2<nthreads2;++kt2)
		{
			ThreadArgs *arg=args+kt2;
			int kx, ky;

			arg->blockidx=kx=kt+kt2;
			ky=kx/xblocks;
			kx%=xblocks;
			arg->x1=BLOCKSIZE*kx;
			arg->y1=BLOCKSIZE*ky;
			arg->x2=MINVAR(arg->x1+BLOCKSIZE, image->iw);
			arg->y2=MINVAR(arg->y1+BLOCKSIZE, image->ih);
			if(!fwd)
			{
				int size=0;
				memcpy(&size, cbuf+sizeof(int)*((ptrdiff_t)kt+kt2), sizeof(int));
				arg->decstart=cbuf+start;
				start+=size;
				arg->decend=cbuf+start;
			}
		}
#ifdef DISABLE_MT
		for(int k=0;k<nthreads2;++k)
			block_thread(args+k);
#else
		void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
		mt_finish(ctx);
#endif
		if(fwd)
		{
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				int blocksize=((arg->y2-arg->y1)*(arg->x2-arg->x1)*image->nch*image->depth+7)>>3;
				if(loud)
				{
					if(nblocks<MAXPRINTEDBLOCKS)
					{
						//if(!(kt+kt2))
						//	printf("block,  nrows,  usize,     best  ->  actual,  (actual-best)\n");
						printf(
							"[%3d]  %4d*%4d  %8d  %11.2lf -> %8zd bytes  (%+11.2lf)  %s\n",
							arg->blockidx,
							arg->x2-arg->x1,
							arg->y2-arg->y1,
							blocksize,
							arg->bestsize,
							arg->list.nobj,
							(double)arg->list.nobj-arg->bestsize,
							rct_names[arg->bestrct]
						);
					}
					bestsize+=arg->bestsize;
				}
				memcpy(data[0]->data+start+sizeof(int)*((ptrdiff_t)kt+kt2), &arg->list.nobj, sizeof(int));
				dlist_appendtoarray(&arg->list, data);
				dlist_clear(&arg->list);
			}
		}
	}
	if(loud)
	{
		ptrdiff_t usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
		t0=time_sec()-t0;
		if(fwd)
		{
			ptrdiff_t csize=data[0]->count-start;
			printf("Best %15.2lf (%+13.2lf) bytes\n", bestsize, csize-bestsize);
			printf("Size %12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
			printf("Mem usage: ");
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
		}
		printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		_mm_free(arg->pixels);
		free(arg->hist);
		//free(arg->stats);
#ifdef ENABLE_OLS
		for(int k=0;k<image->nch;++k)
			ols4_free(arg->ols4+k);
#endif
	}
	free(args);
	return 0;
}