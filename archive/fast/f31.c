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


	#define DISABLE_OLS

#include"ac.h"
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
#define BLOCKDX 768
#define BLOCKDY 768

//RCT1~7:
//Cr-=y		Cr = [ 1	-1	0].RGB
//y+=Cr>>1	y  = [ 1/2	 1/2	0].RGB
//Cb-=y		Cb = [-1/2	-1/2	1].RGB
	#define LUMA_UPDATE_1(Cb, Cr) ((Cb)>>1)			//RCT1	y = [1/4	1/4	 1/2].RGB
	#define LUMA_UPDATE_2(Cb, Cr) ((2*(Cb)-(Cr)+4)>>3)	//RCT2	y = [1/4	1/2	 1/4].RGB
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

//#define PREDLIST\
//	PRED(W)\
//	PRED(cgrad)\
//	PRED(wgrad)
//typedef enum _PredType
//{
//#define PRED(NAME) PRED_##NAME,
//	PREDLIST
//#undef  PRED
//	PRED_COUNT,
//} PredType;
//static const char *pred_names[PRED_COUNT]=
//{
//#define PRED(NAME) #NAME,
//	PREDLIST
//#undef  PRED
//};


//WG
#if 0
#define WG_DECAY_NUM	493
#define WG_DECAY_SH	9

#define WG_NPREDS	11
#define WG_PREDLIST\
	WG_PRED(132, N+W-NW)\
	WG_PRED(176, N+W-NW+((eN+eW-eNW+16)>>5))\
	WG_PRED(176, N+eN)\
	WG_PRED(100, N)\
	WG_PRED(176, W+eW)\
	WG_PRED(100, W)\
	WG_PRED(165, W+NE-N)\
	WG_PRED(220, N+NE-NNE)\
	WG_PRED(165, 3*(N-NN)+NNN)\
	WG_PRED(165, 3*(W-WW)+WWW)\
	WG_PRED(176, (W+NEEE)/2)
//	WG_PRED(  0, N+NE-NNE+((eN+eNE-eNNE)>>2))
static void wg_init(int *weights)
{
	int j=0;
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
	WG_PREDLIST
#undef  WG_PRED
}
static int wg_predict(
	const int *weights,
	short **rows, const int stride, int kc2, int depth,
	const int *perrors, int *preds
)
{
	double pred2=0, wsum=0;
	int j=0, pred;
	short
		NNNWWWW	=rows[3][kc2-4*stride+0],
		NNWWWW	=rows[2][kc2-4*stride+0],
		NNN	=rows[3][kc2+0*stride+0],
		NN	=rows[2][kc2+0*stride+0],
		NNE	=rows[2][kc2+1*stride+0],
		NNEEE	=rows[2][kc2+3*stride+0],
		NW	=rows[1][kc2-1*stride+0],
		N	=rows[1][kc2+0*stride+0],
		NE	=rows[1][kc2+1*stride+0],
		NEE	=rows[1][kc2+2*stride+0],
		NEEE	=rows[1][kc2+3*stride+0],
		WWW	=rows[0][kc2-3*stride+0],
		WW	=rows[0][kc2-2*stride+0],
		W	=rows[0][kc2-1*stride+0],
		eNNE	=rows[2][kc2+1*stride+1],
		eNW	=rows[1][kc2-1*stride+1],
		eN	=rows[1][kc2+0*stride+1],
		eNE	=rows[1][kc2+1*stride+1],
		eW	=rows[0][kc2-1*stride+1];

#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST
#undef  WG_PRED
	
	for(int k=0;k<WG_NPREDS;++k)
	{
		double weight=(double)weights[k]/(perrors[k]+1);
		pred2+=weight*preds[k];
		wsum+=weight;
	}
	pred2/=wsum;
#ifdef _MSC_VER
	pred=_cvt_dtoi_fast(pred2);
#else
	pred=(int)pred2;
#endif
	CLAMP3_32(pred, pred, N, W, NE);
	return pred;
}
#if 0
static int wg_predict(
	const double *weights,
	int NNN, int NN, int NNE,
	int NW, int N, int NE, int NEE, int NEEE,
	int WWW, int WW, int W,
	const int *perrors, int *preds
)
{
	int pred;
	double pred2=0, wsum=0;
	int j=0;
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST
#undef  WG_PRED
	
	for(int k=0;k<WG_NPREDS;++k)
	{
		double weight=weights[k]/(perrors[k]+1);
		pred2+=weight*preds[k];
		wsum+=weight;
	}
	pred2/=wsum;
#ifdef _MSC_VER
	pred=_cvt_dtoi_fast(pred2);
#else
	pred=(int)pred2;
#endif
	CLAMP3_32(pred, pred, N, W, NE);
	return pred;
}
#endif
static void wg_update(int curr, const int *preds, int *perrors, int *weights)
{
	for(int k=0;k<WG_NPREDS;++k)
		perrors[k]=(perrors[k]+abs(curr-preds[k]))*WG_DECAY_NUM>>WG_DECAY_SH;
}
#endif


#define OLS_NPARAMS 10
#define OLS_MEM (1+OLS_NPARAMS+OLS_NPARAMS*OLS_NPARAMS)
#define OLS_ALPHA 5
#define OLS_BETA 3
#define FB1 12
#define FB2 2
#define FB3 (FB1-FB2)
#define BIAS_INIT (2<<FB2)
//#define N_MAPPER 20
#define K_STEP 3
#define K_MAX (15/K_STEP)
static int ols_predict(long long *bufe, long long *buff, const long long *vec, long long bias, long long *pred, int half)
{
	long long buf[OLS_MEM], *vecb=buf+1, *matA=vecb+OLS_NPARAMS;
	for(int k=1;k<OLS_MEM;++k)
		buf[k]=bufe[k]+buff[k];
	for(int k=0;k<OLS_NPARAMS;++k)
	{
		vecb[k]+=bias<<FB3;
		matA[(OLS_NPARAMS+1)*k]+=bias*OLS_NPARAMS;
	}
	for(int it=0;it<OLS_NPARAMS-1;++it)//Gaussian elimination	A*x = b -> U*x = b
	{
		int kp=it;
		for(int ky=it+1;ky<OLS_NPARAMS;++ky)
		{
			if(llabs(matA[OLS_NPARAMS*ky+it])>llabs(matA[OLS_NPARAMS*kp+it]))
				kp=ky;
		}
		if(kp!=it)
		{
			long long temp;
			for(int kx=it;kx<OLS_NPARAMS;++kx)
				SWAPVAR(matA[OLS_NPARAMS*kp+kx], matA[OLS_NPARAMS*it+kx], temp);
			SWAPVAR(vecb[kp], vecb[it], temp);
		}

		long long pivot=matA[(OLS_NPARAMS+1)*it];
		if(!pivot)
			return 0;
		for(int ky=it+1;ky<OLS_NPARAMS;++ky)
		{
			long long other=matA[OLS_NPARAMS*ky+it];
			matA[OLS_NPARAMS*ky+it]=0;
			for(int kx=it+1;kx<OLS_NPARAMS;++kx)
				matA[OLS_NPARAMS*ky+kx]-=matA[OLS_NPARAMS*it+kx]*other/pivot;
			vecb[ky]-=vecb[it]*other/pivot;
		}
	}
	for(int k=OLS_NPARAMS-1;k>=0;--k)//solve U*x = b
	{
		long long diag=matA[(OLS_NPARAMS+1)*k];
		if(!diag)
			return 0;
		for(int k2=0;k2<k;++k2)
		{
			long long other=matA[OLS_NPARAMS*k2+k];
			matA[OLS_NPARAMS*k2+k]=0;
			if(other)
				vecb[k2]-=vecb[k]*other/diag;
		}
	}
	long long x=0;
	for(int k=0;k<OLS_NPARAMS;++k)
	{
		long long diag=matA[(OLS_NPARAMS+1)*k];
		x+=((vecb[k]*vec[k]<<FB2)+(diag>>1))/diag;
	}
	CLAMP2(x, -half, half-1);
	*pred=x;
	return 1;
}
static void entropy_update(ArithmeticCoder *ec, unsigned *ca, unsigned *cb, int qw, int *bit, int fwd)
{
	int inc;
	int p0a=((ca[0]+1)<<16)/(ca[0]+ca[1]+2);
	int p0b=((cb[0]+1)<<16)/(cb[0]+cb[1]+2);
	int p0=p0a+(((p0b-p0a)*qw+(1<<5>>1))>>5);
	CLAMP2_32(p0, p0, 1, 0xFFFF);
	
	if(fwd)
		ac_enc_bin(ec, p0, *bit);
	else
		*bit=ac_dec_bin(ec, p0);

	inc=(1<<5)-qw;
	if(*bit^(ca[0]<ca[1])||ca[0]==ca[1])
		ca[*bit]+=inc;
	else
	{
		//ca[!*bit]=(ca[!*bit]+1)>>1;
		ca[*bit]+=inc;
		//if(ca[!*bit]>inc>>1)
		//	ca[!*bit]-=inc>>1;
		//else
		//	ca[!*bit]=0;
		//ca[*bit]+=(1<<5)-qw;
	}
	inc=qw;
	if(*bit^(cb[0]<cb[1])||cb[0]==cb[1])
		cb[*bit]+=inc;
	else
	{
		//cb[!*bit]=(cb[!*bit]+1)>>1;
		cb[*bit]+=inc;
	}
	//ca[*bit]+=(1<<5)-qw;
	//cb[*bit]+=qw;
	if(ca[0]+ca[1]>256<<5)
	{
		ca[0]=(ca[0]+1)>>1;
		ca[1]=(ca[1]+1)>>1;
	}
	if(cb[0]+cb[1]>256<<5)
	{
		cb[0]=(cb[0]+1)>>1;
		cb[1]=(cb[1]+1)>>1;
	}
}


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

	long long ols_bufe[4][OLS_MEM], *ols_bufb, *ols_buff;
	int sse[4][2048];

	//AC
	int tlevels, clevels, statssize;
	unsigned *stats;
	DList list;
	const unsigned char *decstart, *decend;

	//aux
	int blockidx;
	int bestrct;
	//int predidx[3];
	//double usizes[3];
	double csizes[3], bestsize;

	//WG
	//int wg_weights[OCH_COUNT*WG_NPREDS];
	//int wg_errors[OCH_COUNT*WG_NPREDS];
} ThreadArgs;
static void block_thread(void *param)
{
	ArithmeticCoder ec;
	ThreadArgs *args=(ThreadArgs*)param;
	Image const *image=args->fwd?args->src:args->dst;
	int depths[OCH_COUNT]={0}, halfs[OCH_COUNT]={0};
	double bestsize=0;
	int bestrct=0;
	const int *combination;
	//int combination[12]={0}, predidx[4]={0}, flag=0;
	int res=image->iw*(args->y2-args->y1);
	int nctx=args->clevels*args->clevels, cdfstride=args->tlevels+1, chsize=nctx*cdfstride;

	for(int kc=0;kc<OCH_COUNT;++kc)
	{
		int depth=image->depth+och_info[kc][II_INFLATION];
		UPDATE_MIN(depth, 16);
		depths[kc]=depth;
		halfs[kc]=1<<depth>>1;
	}
	if(args->fwd)//encode
	{
		int nlevels[OCH_COUNT]={0};
		double csizes[OCH_COUNT]={0};
		int predsel[OCH_COUNT]={0};

		for(int kc=0;kc<OCH_COUNT;++kc)
			nlevels[kc]=1<<depths[kc];

		//for(int kc=0;kc<OCH_COUNT;++kc)
		//	wg_init(args->wg_weights+WG_NPREDS*kc);
		//memset(args->wg_errors, 0, sizeof(args->wg_errors));

		memset(args->hist, 0, args->histsize);
		memset(args->pixels, 0, args->bufsize);
		for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)//analysis loop
		{
			ALIGN(16) short *rows[]=
			{
				args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*OCH_COUNT*2,
				args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*OCH_COUNT*2,
				args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*OCH_COUNT*2,
				args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*OCH_COUNT*2,
			};
			short input[ICH_COUNT]={0};
			//int preds[PRED_COUNT]={0};
			//int wg_preds[WG_NPREDS]={0};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				int cidx=3;
				short
					*NNN	=rows[3]+0*OCH_COUNT*2,
					*NNW	=rows[2]-1*OCH_COUNT*2,
					*NN	=rows[2]+0*OCH_COUNT*2,
					*NNE	=rows[2]+1*OCH_COUNT*2,
					*NNEE	=rows[2]+2*OCH_COUNT*2,
					*NW	=rows[1]-1*OCH_COUNT*2,
					*N	=rows[1]+0*OCH_COUNT*2,
					*NE	=rows[1]+1*OCH_COUNT*2,
					*NEE	=rows[1]+2*OCH_COUNT*2,
					*NEEE	=rows[1]+3*OCH_COUNT*2,
					*WWW	=rows[0]-3*OCH_COUNT*2,
					*WW	=rows[0]-2*OCH_COUNT*2,
					*W	=rows[0]-1*OCH_COUNT*2,
					*curr	=rows[0]+0*OCH_COUNT*2;
				if(ky<=args->y1+2)
				{
					if(ky<=args->y1+1)
					{
						if(ky==args->y1)
							NEEE=NEE=NE=NW=N=W;
						NN=N;
						NNE=NE;
					}
					NNN=NN;
				}
				if(kx<=2)
				{
					if(kx<=1)
					{
						if(!kx)
							NW=W=N;
						WW=W;
					}
					WWW=WW;
				}
				if(kx>=image->iw-3)
				{
					if(kx>=image->iw-1)
					{
						NNE=NN;
						NE=N;
					}
					NEEE=NE;
				}
				input[ICH_R]=image->data[idx+0];//r
				input[ICH_G]=image->data[idx+1];//g
				input[ICH_B]=image->data[idx+2];//b
				for(int kp=0;kp<6;++kp)//RCT1~7		r-=g; g+=r>>1; b-=g; g+=(A*r+B*b+C)>>SH;
				{
					int j=0;
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
#if 0
					input[cidx+0]+=input[cidx+7]>>1;			//RCT1	y = [1/4	1/4	 1/2].RGB
					input[cidx+1]+=(2*input[cidx+7]-input[cidx+8]+4)>>3;	//RCT2	y = [1/4	1/2	 1/4].RGB
				//	input[cidx+2]+=(2*input[cidx+7]+input[cidx+8]+4)>>3;	//RCT3a	y = [1/2	1/4	 1/4].RGB	never used
					input[cidx+2]+=(input[cidx+7]+input[cidx+8]+2)>>2;	//RCT3b	y = [5/8	1/8	 1/4].RGB
					input[cidx+3]+=input[cidx+7]/3;				//RCT4	y = [1/3	1/3	 1/3].RGB
					input[cidx+4]+=(3*input[cidx+7]+4)>>3;			//RCT5	y = [5/16	5/16	 6/16].RGB
					input[cidx+5]+=(7*input[cidx+7]+8)>>4;			//RCT6	y = [9/32	9/32	14/32].RGB
					input[cidx+6]+=(10*input[cidx+7]-input[cidx+8]+16)>>5;	//RCT7	y = [5/16	6/16	 5/16].RGB
#endif
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
				for(int kc=0;kc<OCH_COUNT;++kc)
				{
					int
						kc2=kc<<1,
						target=input[och_info[kc][II_TARGET]],
						offset=(input[och_info[kc][II_HELPER1]]+input[och_info[kc][II_HELPER2]])>>och_info[kc][II_HSHIFT];
					
					int pred;
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);

					//preds[PRED_W]=W[kc2];
					//MEDIAN3_32(preds[PRED_cgrad], N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					//preds[PRED_wgrad]=wg_predict(
					//	args->wg_weights+WG_NPREDS*kc,
					//	rows, OCH_COUNT*2, kc2, depths[kc],
					//	args->wg_errors+WG_NPREDS*kc, wg_preds
					//);

					if(offset)
					{
						pred+=offset;
						CLAMP2_32(pred, pred, -halfs[kc], halfs[kc]-1);

						//for(int kp=0;kp<PRED_COUNT;++kp)
						//{
						//	preds[kp]+=offset;
						//	CLAMP2_32(preds[kp], preds[kp], -halfs[kc], halfs[kc]-1);
						//}
					}
					int *curr_hist=args->hist+args->histindices[kc];
					int val=target-pred;
					val+=halfs[kc];
					val&=nlevels[kc]-1;
					++curr_hist[val];

					target-=offset;
					curr[kc2+0]=target;
					//for(int kp=0;kp<PRED_COUNT;++kp)
					//{
					//	int *curr_hist=args->hist+args->histindices[kc*PRED_COUNT+kp];
					//	int val=target-preds[kp];
					//	val+=halfs[kc];
					//	val&=nlevels[kc]-1;
					//	++curr_hist[val];
					//}
					//curr[kc2+1]=target-preds[PRED_wgrad];
					//target-=offset;
					//wg_update(target, wg_preds, args->wg_errors+WG_NPREDS*kc, args->wg_weights+WG_NPREDS*kc);
					//curr[kc2+0]=target;
				}
				rows[0]+=OCH_COUNT*2;
				rows[1]+=OCH_COUNT*2;
				rows[2]+=OCH_COUNT*2;
				rows[3]+=OCH_COUNT*2;
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
		//for(int kc=0;kc<OCH_COUNT;++kc)//select best predictors
		//{
		//	int bestpred=0;
		//	for(int kp=1;kp<PRED_COUNT;++kp)
		//	{
		//		if(csizes[kc*PRED_COUNT+bestpred]>csizes[kc*PRED_COUNT+kp])
		//			bestpred=kp;
		//	}
		//	predsel[kc]=bestpred;
		//}
		for(int kt=0;kt<RCT_COUNT;++kt)//select best R.C.T
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
		//memcpy(combination, rct_combinations[bestrct], sizeof(combination));
		//predidx[0]=predsel[combination[0]];
		//predidx[1]=predsel[combination[1]];
		//predidx[2]=predsel[combination[2]];
		//flag=bestrct;
		//flag=flag*PRED_COUNT+predidx[2];
		//flag=flag*PRED_COUNT+predidx[1];
		//flag=flag*PRED_COUNT+predidx[0];//19*3*3*3 = 513

		args->bestrct=bestrct;
		//args->predidx[0]=predidx[0];
		//args->predidx[1]=predidx[1];
		//args->predidx[2]=predidx[2];
		//args->usizes[0]=(res*image->depth+7)>>3;
		//args->usizes[1]=(res*image->depth+7)>>3;
		//args->usizes[2]=(res*image->depth+7)>>3;
		args->csizes[0]=csizes[combination[0]];
		args->csizes[1]=csizes[combination[1]];
		args->csizes[2]=csizes[combination[2]];
		args->bestsize=bestsize;
		if(args->loud)
		{
			printf("Y %5d~%5d  best %12.2lf bytes  %s\n",
				args->y1, args->y2,
				bestsize,
				rct_names[bestrct]
			);

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
		dlist_init(&args->list, 1, 1024, 0);
		//dlist_push_back(&args->list, &flag, sizeof(char[2]));
		ac_enc_init(&ec, &args->list);
		ac_enc_bypass_NPOT(&ec, bestrct, RCT_COUNT);
	}
	else//decode
	{
		ac_dec_init(&ec, args->decstart, args->decend);
		bestrct=ac_dec_bypass_NPOT(&ec, RCT_COUNT);
		combination=rct_combinations[bestrct];
	}
	long long bias=2LL<<FB2;
	static const int c_thresholds[]={1, 3, 9, 20, 50, 110, 300, 800};
	static const int q_bins[]=
	{
		0, 2, 4, 7, 10, 14, 20, 26, 34, 42, 52, 64, 78, 95, 135, 200,
	};
	memset(args->stats, 0, args->statssize);
#ifndef DISABLE_OLS
	memset(args->ols_bufb, 0, image->iw*sizeof(long long[OLS_MEM])*image->nch);
	memset(args->ols_buff, 0, image->iw*sizeof(long long[OLS_MEM])*image->nch);
#endif
	//for(int kc=0;kc<3;++kc)
	//{
	//	if(predidx[kc]==PRED_wgrad)
	//		wg_init(args->wg_weights+WG_NPREDS*kc);
	//}
	//memset(args->wg_errors, 0, sizeof(args->wg_errors));
	memset(args->pixels, 0, args->bufsize);
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) short *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		short yuv[5]={0};
		//int wg_preds[WG_NPREDS]={0};
		int token=0, bypass=0, nbits=0;

		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int idx=image->nch*(image->iw*ky+kx);
			short
				*NNN	=rows[3]+0*4*2,
				*NNWW	=rows[2]-2*4*2,
				*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NNEE	=rows[2]+2*4*2,
				*NWW	=rows[1]-2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*NEE	=rows[1]+2*4*2,
				*NEEE	=rows[1]+3*4*2,
				*WWW	=rows[0]-3*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
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
				//	NNEEE=NEEE;
				}
				NNN=NN;
			}
			//if(kx<=args->x1+3)
			//{
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
			//	WWWW=WWW;
			//}
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
						//yuv[0]+=yuv[1]>>1;
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
						//yuv[0]+=(2*yuv[1]-yuv[2]+4)>>3;//A710
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
						//yuv[0]+=(2*yuv[1]+yuv[2]+4)>>3;
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
						//yuv[0]+=yuv[1]/3;
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
						//yuv[0]+=(3*yuv[1]+4)>>3;
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
						//yuv[0]+=(7*yuv[1]+8)>>4;
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
						//yuv[0]+=(10*yuv[1]-yuv[2]+16)>>5;
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
			//if(idx==12)//
			//if(idx==760758)
			//	printf("");
			for(int kc=0;kc<image->nch;++kc)
			{
				int ch=combination[kc];
				int kc2=kc<<1;
				int offset=(yuv[combination[kc+4]]+yuv[combination[kc+8]])>>och_info[ch][II_HSHIFT];
				
				int qu, qv, qw;
				int ctx;

				int depth=depths[ch], half=halfs[ch];
				long long nb[OLS_NPARAMS]=
				{
				//	NNWW	[kc2],
					NNW	[kc2],
					NN	[kc2],
					NNE	[kc2],
				//	NNEE	[kc2],
					NWW	[kc2],
					NW	[kc2],
					N	[kc2],
					NE	[kc2],
					NEE	[kc2],
					WW	[kc2],
					W	[kc2],
				};
				int pred, pred2, sym, error;
				long long preds[2];
#ifndef DISABLE_OLS
				int success[2];
				long long
					*bufe=args->ols_bufe[kc],
					*bufb=args->ols_bufb+OLS_MEM*(image->nch*kx+kc),
					*buff=args->ols_buff+OLS_MEM*(image->nch*kx+kc);
				long long bias1, bias2;
				bias1=bias*21/22;
				bias2=bias*22/21;
				CLAMP2(bias1, 0, bias-1);
				CLAMP2(bias2, bias+1, 1024<<FB2);
				UPDATE_MIN(bias1, 1024<<FB2);
				UPDATE_MAX(bias2, 0);
				success[0]=ols_predict(bufe, buff, nb, bias1, preds+0, half);
				success[1]=ols_predict(bufe, buff, nb, bias2, preds+1, half);
				if(success[0])
					pred=(int)((preds[0]+(1LL<<FB1>>1))>>FB1);
				else//GAP from NBLIC
#endif
				{
					int csum=0, cmin=0xFFFFFF, cost, ehalf=half<<4, preda=0, wt;
					pred=9*(N[kc2]+W[kc2])+2*(NE[kc2]-NW[kc2])-(NN[kc2]+WW[kc2]);
					CLAMP2_32(pred, pred, -ehalf, ehalf-1);

					cost=2*(abs(W[kc2]-WW[kc2])+abs(NW[kc2]-NWW[kc2])+abs(N[kc2]-NW[kc2])+abs(NE[kc2]-N[kc2]));//180 degrees
					csum+=cost;
					if(cmin>cost)
					{
						cmin=cost;
						preda=W[kc2]*2;
					}

					cost=2*(abs(W[kc2]-NW[kc2])+abs(NW[kc2]-NNW[kc2])+abs(N[kc2]-NN[kc2])+abs(NE[kc2]-NNE[kc2]));//90 degrees
					csum+=cost;
					if(cmin>cost)
					{
						cmin=cost;
						preda=N[kc2]*2;
					}

					cost=2*(abs(W[kc2]-NWW[kc2])+abs(NW[kc2]-NNWW[kc2])+abs(N[kc2]-NNW[kc2])+abs(NE[kc2]-NN[kc2]));//135 degrees
					csum+=cost;
					if(cmin>cost)
					{
						cmin=cost;
						preda=NW[kc2]*2;
					}

					cost=2*(abs(W[kc2]-N[kc2])+abs(NW[kc2]-NN[kc2])+abs(N[kc2]-NNE[kc2])+abs(NE[kc2]-NNEE[kc2]));//45 degrees
					csum+=cost;
					if(cmin>cost)
					{
						cmin=cost;
						preda=NE[kc2]*2;
					}

					cost=abs(2*W[kc2]-WW[kc2]-NWW[kc2])+abs(2*NW[kc2]-NWW[kc2]-NNWW[kc2])+abs(2*N[kc2]-NW[kc2]-NNW[kc2])+abs(2*NE[kc2]-N[kc2]-NN[kc2]);//157.5 degrees
					csum+=cost;
					if(cmin>cost)
					{
						cmin=cost;
						preda=W[kc2]+NW[kc2];
					}

					cost=abs(2*W[kc2]-NWW[kc2]-NW[kc2])+abs(2*NW[kc2]-NNWW[kc2]-NNW[kc2])+abs(2*N[kc2]-NNW[kc2]-NN[kc2])+abs(2*NE[kc2]-NN[kc2]-NNE[kc2]);//112.5 degrees
					csum+=cost;
					if(cmin>cost)
					{
						cmin=cost;
						preda=NW[kc2]+N[kc2];
					}

					cost=abs(2*W[kc2]-NW[kc2]-N[kc2])+abs(2*NW[kc2]-NNW[kc2]-NN[kc2])+abs(2*N[kc2]-NN[kc2]-NNE[kc2])+abs(2*NE[kc2]-NNE[kc2]-NNEE[kc2]);//67.5 degrees
					csum+=cost;
					if(cmin>cost)
					{
						cmin=cost;
						preda=N[kc2]+NE[kc2];
					}

					csum-=7*cmin;
					for(wt=0;wt<(int)_countof(c_thresholds)&&c_thresholds[wt]*half<csum*8;++wt);
					pred=(8*wt*preda+(8-wt)*pred+64)>>7;
					preds[0]=(long long)pred<<FB1;
				}
				{
					int delta=abs(W[kc2]-WW[kc2])+abs(N[kc2]-NW[kc2])+abs(N[kc2]-NE[kc2])+abs(W[kc2]-NW[kc2])+abs(N[kc2]-NN[kc2])+abs(NE[kc2]-NNE[kc2])+abs(W[kc2+1])+abs(N[kc2+1]), qd;
					for(qd=0;qd<(int)_countof(q_bins)-1&&q_bins[qd]<<depth<delta<<8;++qd);
					qu=qv=qd;
					qw=0;
					if(delta<q_bins[qd])
					{
						qw=((delta-q_bins[qd-1])<<5)/(q_bins[qd]-q_bins[qd-1]);
						if(qw<1<<5>>1)
							qu=qd-1;
						else
						{
							qv=qd-1;
							qw=(1<<5)-qw;
						}
					}
				}
				ctx=qu>>1;
				ctx=ctx<<1|(pred>2*N[kc2]-NN[kc2]);
				ctx=ctx<<1|(pred>2*W[kc2]-WW[kc2]);
				ctx=ctx<<1|(pred>NN[kc2]);
				ctx=ctx<<1|(pred>WW[kc2]);
				ctx=ctx<<1|(pred>NE[kc2]);
				ctx=ctx<<1|(pred>NW[kc2]);
				ctx=ctx<<1|(pred>N[kc2]);
				ctx=ctx<<1|(pred>W[kc2]);
				int *pitem=args->sse[kc]+ctx, item=*pitem;
				int sign=item>>7&1;
				pred2=pred+(item>>8)+sign;
				CLAMP2_32(pred2, pred2, -half, half-1);
				int pred3=pred+offset;
				CLAMP2_32(pred3, pred3, -half, half-1);
				//if(idx==760758)//
				//	printf("");
				if(args->fwd)
				{
					curr[kc2]=yuv[kc];
					curr[kc2+1]=error=yuv[kc]-pred3;
#if 1
					{
						int upred=half-abs(pred3), aval=abs(error);
						if(aval<=upred)
						{
							sym=error;
							{
								int negmask=-((pred2<pred)&(sym!=-half));//sign is flipped if SSE correction was negative, to skew the histogram
								sym^=negmask;
								sym-=negmask;
							}
							sym=sym<<1^sym>>31;//pack sign
						}
						else
							sym=upred+aval;//error sign is known
					}
				}
				else
					sym=0, error=0;
#endif
				int i, bit;
				unsigned *curr_stats=args->stats+(16LL*2*kc<<(image->depth+1));
				if(qv/K_STEP!=qu/K_STEP)
					qv=qu;
				i=0;
				for(int kb=depth;kb>=0;--kb)
				{
					if(args->fwd)
						bit=sym>>kb&1;
					entropy_update(&ec, curr_stats+(qu<<depth|i)*2, curr_stats+(qv<<depth|i)*2, qw, &bit, args->fwd);
					if(!args->fwd)
						sym|=bit<<kb;
					if(!i&&bit)
						break;
					i+=i+!bit;
				}
#if 0
				int k;
				for(i=0;;)
				{
					k=qu/K_STEP;
					if(args->fwd)
						bit=i>>K_MAX<sym>>k;
					entropy_update(&ec, curr_stats+(qu<<depth|i)*2, curr_stats+(qv<<depth|i)*2, qw, &bit, args->fwd);
					if(!bit)
						break;
					i+=1<<K_MAX;
					if(i>=half<<1)
					{
						i>>=1;
						qv=qu=(k+1)*K_STEP;
					}
				}
				if(!args->fwd)
					sym=i>>K_MAX<<k;
				for(++i, --k;k>=0;--k)
				{
					if(args->fwd)
						bit=sym>>k&1;
					entropy_update(&ec, curr_stats+(qu<<depth|i)*2, curr_stats+(qv<<depth|i)*2, qw, &bit, args->fwd);
					if(!args->fwd)
						sym|=bit<<k;
					i+=bit?1<<k:1;
				}
#endif
				if(!args->fwd)
				{
#if 1
					{
						int upred=half-abs(pred3), negmask=0;
						if(sym<=(upred<<1))
						{
							error=sym>>1^-(sym&1);
							negmask=-((pred2<pred)&(error!=-half));
						}
						else
						{
							error=sym-upred;
							negmask=-(pred3>0);
						}
						error^=negmask;
						error-=negmask;
					}
#endif
					curr[kc2+0]=yuv[kc]=error+pred3;
					curr[kc2+1]=error;
				}
				curr[kc2]-=offset;
				*pitem=(item*((1<<7)-1)+(error<<8)+64)>>7;
#ifndef DISABLE_OLS
				{
					int x=curr[kc2];
					long long s_curr=llabs(((long long)x<<FB1)-preds[0]);
					long long sum=(bufe[0]+buff[0])+s_curr*OLS_BETA/(OLS_BETA-1), half_sum;
					long long buf[OLS_MEM], *vecb=buf+1, *matA=vecb+OLS_NPARAMS;
					buf[0]=s_curr;

					sum+=1LL<<FB1;
					CLAMP2(sum, 1<<FB1, 16<<FB1);
					half_sum=sum>>1;
					for(int k=0;k<OLS_NPARAMS;++k)//b = x*n
						vecb[k]=((x*nb[k]<<(4+FB1+FB1))+half_sum)/sum;
					for(int ky=0;ky<OLS_NPARAMS;++ky)//A = n*nT
					{
						for(int kx=0;kx<OLS_NPARAMS;++kx)
							matA[OLS_NPARAMS*ky+kx]=((nb[ky]*nb[kx]<<(4+FB1+FB1))+half_sum)/sum;
					}
					int ab=OLS_BETA;
					for(int k=0;k<OLS_MEM;++k)
					{
						bufb[k]=(bufb[k]*(ab-1LL)+(ab>>1))/ab+buf[k];
						bufe[k]=(bufe[k]*(ab-1LL)+(ab>>1))/ab+buf[k];
						ab=OLS_ALPHA;
					}
					if(success[0]&&success[1])
					{
						long long
							e1=llabs(((long long)x<<FB1)-preds[0]),
							e2=llabs(((long long)x<<FB1)-preds[1]);
						bias=e1<=e2?bias1:bias2;
					}
				}
#endif
				//if(predidx[kc]==PRED_wgrad)
				//	wg_update(curr[kc2], wg_preds, args->wg_errors+WG_NPREDS*kc, args->wg_weights+WG_NPREDS*kc);
			}
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
						//yuv[0]-=yuv[1]>>1;
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
						//yuv[0]-=(2*yuv[1]-yuv[2]+4)>>3;//A710
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
						//yuv[0]-=(2*yuv[1]+yuv[2]+4)>>3;
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
						//yuv[0]-=yuv[1]/3;
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
						//yuv[0]-=(3*yuv[1]+4)>>3;
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
						//yuv[0]-=(7*yuv[1]+8)>>4;
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
						//yuv[0]-=(10*yuv[1]-yuv[2]+16)>>5;
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
	{
		ac_enc_flush(&ec);
		//if(args->loud)
		//	printf("Actual %8zd bytes (%+11.2lf)\n\n", args->list.nobj, args->list.nobj-bestsize);
	}
}
int f31_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
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
	size_t usize;
	
	ncores=query_cpu_cores();
	usize=((size_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
	xblocks=(image->iw+BLOCKDX-1)/BLOCKDX;
	yblocks=(image->ih+BLOCKDY-1)/BLOCKDY;
	nblocks=xblocks*yblocks, nthreads=MINVAR(nblocks, ncores);
	coffset=sizeof(int)*nblocks;
	argssize=sizeof(ThreadArgs)*nthreads;
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
		statssize=image->nch*sizeof(int[16*2])<<(image->depth+1);//chroma inflated
		//statssize=clevels*clevels*(tlevels+1)*image->nch*(int)sizeof(int);
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
		arg->bufsize=sizeof(short[4*OCH_COUNT*2])*(image->iw+16LL);//4 padded rows * OCH_COUNT * {pixels, wg_errors}
		arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));

		arg->histsize=histsize*sizeof(int);
		arg->hist=(int*)malloc(arg->histsize);

		arg->statssize=statssize;
		arg->stats=(unsigned*)malloc(statssize);
		
#ifndef DISABLE_OLS
		int ols_size=image->iw*sizeof(long long[OLS_MEM])*image->nch;
		arg->ols_bufb=(long long*)malloc(ols_size);
		arg->ols_buff=(long long*)malloc(ols_size);
#endif
		if(!arg->pixels||!arg->hist||!arg->stats
#ifndef DISABLE_OLS
			||!arg->ols_bufb||!arg->ols_buff
#endif
		)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=arg->bufsize;
		memusage+=arg->histsize;
		memusage+=arg->statssize;
#ifndef DISABLE_OLS
		memusage+=ols_size*2LL;
#endif
		
		arg->tlevels=tlevels;
		arg->clevels=clevels;
		memcpy(arg->histindices, histindices, sizeof(arg->histindices));
		arg->fwd=fwd;
#ifdef DISABLE_MT
		arg->loud=loud;
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
			arg->x1=BLOCKDX*kx;
			arg->y1=BLOCKDY*ky;
			arg->x2=MINVAR(arg->x1+BLOCKDX, image->iw);
			arg->y2=MINVAR(arg->y1+BLOCKDY, image->ih);
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
				int blocksize=(arg->y2-arg->y1)*image->iw*image->nch*image->depth/8;
				if(loud)
				{
					if(!(kt+kt2))
						printf("block,  nrows,  usize,     best  ->  actual,  (actual-best)\n");
					printf(
						"[%3d]  %4d  %8d  %11.2lf -> %8zd bytes  (%+11.2lf)  %s\n",
						kt+kt2,
						arg->y2-arg->y1,
						blocksize,
						arg->bestsize,
						arg->list.nobj,
						(double)arg->list.nobj-arg->bestsize,
						rct_names[arg->bestrct]
					);
					//printf("[%3d]  %13.2lf -> %10zd bytes (%+13.2lf)\n", kt+kt2, arg->bestsize, arg->list.nobj, arg->list.nobj-arg->bestsize);
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
			printf("Mem usage: ");
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
			printf("Best %15.2lf (%+13.2lf) bytes\n", bestsize, csize-bestsize);
			printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
		}
		printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		_mm_free(arg->pixels);
		free(arg->hist);
		free(arg->stats);
#ifndef DISABLE_OLS
		free(arg->ols_bufb);
		free(arg->ols_buff);
#endif
	}
	free(args);
	return 0;
}