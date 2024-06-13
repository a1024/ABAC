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


#include"ac.h"
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
#define BLOCKSIZE 64

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
static const int rct_combinations[RCT_COUNT][9]=
{
#define RCT(NAME, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2) {YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(NAME, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2) #NAME,
	RCTLIST
#undef  RCT
};

#define PREDLIST\
	PRED(W)\
	PRED(cgrad)\
	PRED(wgrad)
typedef enum _PredType
{
#define PRED(NAME) PRED_##NAME,
	PREDLIST
#undef  PRED
	PRED_COUNT,
} PredType;
static const char *pred_names[PRED_COUNT]=
{
#define PRED(NAME) #NAME,
	PREDLIST
#undef  PRED
};


//WG:

//	#define WG_UPDATE
#define WG_DECAY_NUM	493
#define WG_DECAY_SH	9

#define WG_NPREDS	8
#if 1
#define WG_PREDLIST\
	WG_PRED(1.2, N+W-NW)\
	WG_PRED(1.5, N)\
	WG_PRED(1.5, W)\
	WG_PRED(1.0, W+NE-N)\
	WG_PRED(1.0, N+NE-NNE)\
	WG_PRED(1.0, 3*(N-NN)+NNN)\
	WG_PRED(1.0, 3*(W-WW)+WWW)\
	WG_PRED(1.0, (W+NEEE)/2)
#endif
#if 0
#define WG_PREDLIST\
	WG_PRED(10000, N+W-NW)\
	WG_PRED(10000, N)\
	WG_PRED(10000, W)\
	WG_PRED(10000, W+NE-N)\
	WG_PRED(10000, N+NE-NNE)\
	WG_PRED(10000, 3*(N-NN)+NNN)\
	WG_PRED(10000, 3*(W-WW)+WWW)\
	WG_PRED(10000, (W+NEEE)/2)
#endif
#if 0
#define WG_PREDLIST\
	WG_PRED(211.2, N+W-NW)\
	WG_PRED(264.0, N)\
	WG_PRED(264.0, W)\
	WG_PRED(176.0, W+NE-N)\
	WG_PRED(176.0, N+NE-NNE)\
	WG_PRED(176.0, 3*(N-NN)+NNN)\
	WG_PRED(176.0, 3*(W-WW)+WWW)\
	WG_PRED(176.0, (W+NEEE)/2)
#endif
//	WG_PRED(0.5, NW)
//	WG_PRED(0.5, NE)
static void wg_init(double *weights)
{
	int j=0;
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
	WG_PREDLIST
#undef  WG_PRED
}
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
static void wg_update(int curr, const int *preds, int *perrors, double *weights)
{
#ifdef WG_UPDATE
	int kbest=0, ebest=0;
	//int kworst=0, eworst=0;
#endif
	for(int k=0;k<WG_NPREDS;++k)
	{
		int e2=abs(curr-preds[k]);
		perrors[k]=(perrors[k]+e2)*WG_DECAY_NUM>>WG_DECAY_SH;
#ifdef WG_UPDATE
		if(!k||ebest>e2)
			kbest=k, ebest=e2;
		//if(!k||eworst<e2)
		//	kworst=k, eworst=e2;
		weights[k]+=1/sqrt(e2+1./256);
		//weights[k]+=1./(e2+1);
#endif
	}
#ifdef WG_UPDATE
	if(weights[kbest]>10)
	{
		for(int k=0;k<WG_NPREDS;++k)
			weights[k]*=0.125;
	}
	//if(weights[kbest]>10)//100	352
	//{
	//	for(int k=0;k<WG_NPREDS;++k)
	//		weights[k]*=0.5;
	//}
#endif
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
		*nbits=lgv-CONFIG_MSB+CONFIG_LSB;
		*bypass=val>>CONFIG_LSB&((1LL<<*nbits)-1);
	}
}

static double g_rct_usizes[RCT_COUNT]={0};
static double g_usizes[OCH_COUNT*PRED_COUNT]={0};
static double g_csizes[OCH_COUNT*PRED_COUNT]={0};
//static double g_utotal[OCH_COUNT*PRED_COUNT]={0};
//static double g_ctotal[OCH_COUNT*PRED_COUNT]={0};
typedef struct _ThreadArgs
{
	const Image *image;
//	int x1, x2;//FIXME later
	int y1, y2;
	int fwd, loud;//fwd always true
	
	short *pixels;
	int bufsize;

	double wg_weights[OCH_COUNT*WG_NPREDS];
	int wg_errors[OCH_COUNT*WG_NPREDS];

	int *hist, histsize, histindices[OCH_COUNT*PRED_COUNT+1];
	
	double csizes[OCH_COUNT*PRED_COUNT], bestsize;
	int predsel[OCH_COUNT], bestrct;
} ThreadArgs;
#define CDFSTRIDE 32	//power-of-two
static void update_CDF(const int *hist, unsigned *CDF, int tlevels)
{
	int sum=hist[tlevels], c=0;
	for(int ks=0;ks<tlevels;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*((1LL<<PROB_BITS)-tlevels)/sum)+ks;
		c+=freq;
	}
	CDF[tlevels]=1<<PROB_BITS;
}
static void block_thread(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	Image const *image=args->image;
	int depths[OCH_COUNT]={0}, halfs[OCH_COUNT]={0};
	double bestsize=0;
	int res=image->iw*(args->y2-args->y1);
	int nlevels[OCH_COUNT]={0};

	for(int kc=0;kc<OCH_COUNT;++kc)
	{
		int depth=image->depth+och_info[kc][II_INFLATION];
		UPDATE_MIN(depth, 16);
		depths[kc]=depth;
		halfs[kc]=1<<depth>>1;
	}
	memset(args->pixels, 0, args->bufsize);

	for(int kc=0;kc<OCH_COUNT;++kc)
		nlevels[kc]=1<<depths[kc];

	for(int kc=0;kc<OCH_COUNT;++kc)
		wg_init(args->wg_weights+WG_NPREDS*kc);
	memset(args->wg_errors, 0, sizeof(args->wg_errors));

	memset(args->hist, 0, args->histsize);
	for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
	{
		ALIGN(16) short *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*OCH_COUNT,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*OCH_COUNT,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*OCH_COUNT,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*OCH_COUNT,
		};
		short input[ICH_COUNT]={0};
		int preds[PRED_COUNT]={0};
		int wg_preds[WG_NPREDS]={0};
		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			int cidx=3;
			short
				*NNN	=rows[3]+0*OCH_COUNT,
				*NNW	=rows[2]-1*OCH_COUNT,
				*NN	=rows[2]+0*OCH_COUNT,
				*NNE	=rows[2]+1*OCH_COUNT,
				*NNEE	=rows[2]+2*OCH_COUNT,
				*NW	=rows[1]-1*OCH_COUNT,
				*N	=rows[1]+0*OCH_COUNT,
				*NE	=rows[1]+1*OCH_COUNT,
				*NEE	=rows[1]+2*OCH_COUNT,
				*NEEE	=rows[1]+3*OCH_COUNT,
				*WWW	=rows[0]-3*OCH_COUNT,
				*WW	=rows[0]-2*OCH_COUNT,
				*W	=rows[0]-1*OCH_COUNT,
				*curr	=rows[0]+0*OCH_COUNT;
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
				input[cidx+0]=input[rgb2yuv_permutations[kp][0]];//Y
				input[cidx+7]=input[rgb2yuv_permutations[kp][1]];//U
				input[cidx+8]=input[rgb2yuv_permutations[kp][2]];//V
				input[cidx+8]-=input[cidx+0];
				input[cidx+0]+=input[cidx+8]>>1;
				input[cidx+7]-=input[cidx+0];
				input[cidx+6]=input[cidx+5]=input[cidx+4]=input[cidx+3]=input[cidx+2]=input[cidx+1]=input[cidx+0];
				input[cidx+0]+=input[cidx+7]>>1;
				input[cidx+1]+=(2*input[cidx+7]-input[cidx+8]+4)>>3;
				input[cidx+2]+=(2*input[cidx+7]+input[cidx+8]+4)>>3;
				input[cidx+3]+=input[cidx+7]/3;
				input[cidx+4]+=(3*input[cidx+7]+4)>>3;
				input[cidx+5]+=(7*input[cidx+7]+8)>>4;
				input[cidx+6]+=(10*input[cidx+7]-input[cidx+8]+16)>>5;
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
					target=input[och_info[kc][II_TARGET]],
					offset=(input[och_info[kc][II_HELPER1]]+input[och_info[kc][II_HELPER2]])>>och_info[kc][II_HSHIFT];

				preds[PRED_W]=W[kc];
				MEDIAN3_32(preds[PRED_cgrad], N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
				preds[PRED_wgrad]=wg_predict(
					args->wg_weights+WG_NPREDS*kc,
					NNN[kc], NN[kc], NNE[kc],
					NW[kc], N[kc], NE[kc], NEE[kc], NEEE[kc],
					WWW[kc], WW[kc], W[kc],
					args->wg_errors+WG_NPREDS*kc, wg_preds
				);

				if(offset)
				{
					for(int kp=0;kp<PRED_COUNT;++kp)
					{
						preds[kp]+=offset;
						CLAMP2_32(preds[kp], preds[kp], -halfs[kc], halfs[kc]-1);
					}
				}
				for(int kp=0;kp<PRED_COUNT;++kp)
				{
					int *curr_hist=args->hist+args->histindices[kc*PRED_COUNT+kp];
					int val=target-preds[kp];
					val+=halfs[kc];
					val&=nlevels[kc]-1;
					++curr_hist[val];
				}
				target-=offset;
				wg_update(target, wg_preds, args->wg_errors+WG_NPREDS*kc, args->wg_weights+WG_NPREDS*kc);
				curr[kc]=target;
			}
			rows[0]+=OCH_COUNT;
			rows[1]+=OCH_COUNT;
			rows[2]+=OCH_COUNT;
			rows[3]+=OCH_COUNT;
		}
	}
	memset(args->csizes, 0, sizeof(args->csizes));
	for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)//calculate channel sizes
	{
		int *curr_hist=args->hist+args->histindices[kc];
		int nlevels2=nlevels[kc/PRED_COUNT];
		for(int ks=0;ks<nlevels2;++ks)
		{
			int freq=curr_hist[ks];
			if(freq)
				args->csizes[kc]-=freq*log2((double)freq/res);
		}
		args->csizes[kc]/=8;
	}
	for(int kc=0;kc<OCH_COUNT;++kc)//select best predictors
	{
		int bestpred=0;
		for(int kp=1;kp<PRED_COUNT;++kp)
		{
			if(args->csizes[kc*PRED_COUNT+bestpred]>args->csizes[kc*PRED_COUNT+kp])
				bestpred=kp;
		}
		args->predsel[kc]=bestpred;
	}
	for(int kt=0;kt<RCT_COUNT;++kt)//select best R.C.T
	{
		const int *group=rct_combinations[kt];
		double csize=
			args->csizes[group[0]*PRED_COUNT+args->predsel[group[0]]]+
			args->csizes[group[1]*PRED_COUNT+args->predsel[group[1]]]+
			args->csizes[group[2]*PRED_COUNT+args->predsel[group[2]]];
		if(!kt||bestsize>csize)
			bestsize=csize, args->bestrct=kt;
	}
}
typedef struct _DecorrInfo
{
	unsigned char bestrct, predidx[3];
} DecorrInfo;
int f27_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	ArithmeticCoder ec;
	DList list={0};
	int nblocks=(image->ih+BLOCKSIZE-1)/BLOCKSIZE;
	int bufsize=sizeof(short[4*OCH_COUNT*1])*(image->iw+16LL);//4 padded rows * OCH_COUNT * {pixels}
	ptrdiff_t memusage=0;
	int tlevels=0, clevels=0, statssize=0;
	double bestsize=0;
	DecorrInfo *flags=(DecorrInfo*)malloc(nblocks*sizeof(DecorrInfo));
	int nctx, cdfstride, chsize;
	int *hist=0;
	unsigned *stats=0;
	double wg_weights[4*WG_NPREDS]={0};
	int wg_errors[4*WG_NPREDS]={0};
	short *pixels=0;
	int depths[OCH_COUNT]={0}, halfs[OCH_COUNT]={0};

	if(!flags)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(flags, 0, nblocks*sizeof(DecorrInfo));
	if(fwd)
	{
		int histindices[OCH_COUNT*PRED_COUNT+1]={0}, histsize=0;
		int ncores=query_cpu_cores(), nthreads=MINVAR(nblocks, ncores);
		ptrdiff_t argssize=nthreads*sizeof(ThreadArgs);
		ThreadArgs *args=(ThreadArgs*)malloc(argssize);

		if(!args)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=argssize;
		dlist_init(&list, 1, 1024, 0);
		ac_enc_init(&ec, &list);
#ifdef ENABLE_GUIDE
		guide=image;
#endif
		histsize=0;
		for(int kc=0;kc<OCH_COUNT;++kc)
		{
			int depth=image->depth+och_info[kc][II_INFLATION], nlevels;
			UPDATE_MIN(depth, 16);
			depths[kc]=depth;
			nlevels=1<<depth;
			halfs[kc]=nlevels>>1;
			for(int kp=0;kp<PRED_COUNT;++kp)
			{
				histindices[kc*PRED_COUNT+kp]=histsize;
				histsize+=nlevels;
			}
		}
		histindices[OCH_COUNT*PRED_COUNT]=histsize;

		memset(args, 0, argssize);
		for(int k=0;k<nthreads;++k)
		{
			ThreadArgs *arg=args+k;
			arg->image=image;
			arg->bufsize=bufsize;
			arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));

			arg->histsize=histsize*sizeof(int);
			arg->hist=(int*)malloc(arg->histsize);

			if(!arg->pixels||!arg->hist)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			memusage+=arg->bufsize;
			memusage+=arg->histsize;

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
				arg->y1=BLOCKSIZE*(kt+kt2);
				arg->y2=MINVAR(arg->y1+BLOCKSIZE, image->ih);
			}
#ifdef DISABLE_MT
			for(int k=0;k<nthreads2;++k)
				block_thread(args+k);
#else
			void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
			mt_finish(ctx);
#endif
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				DecorrInfo *flag=flags+kt+kt2;
				const int *comb=rct_combinations[arg->bestrct];
				flag->bestrct=arg->bestrct;//all this just to return 4 bytes
				flag->predidx[0]=arg->predsel[comb[0]];
				flag->predidx[1]=arg->predsel[comb[1]];
				flag->predidx[2]=arg->predsel[comb[2]];
				ac_enc_bypass_NPOT(&ec, flag->bestrct, RCT_COUNT);
				ac_enc_bypass_NPOT(&ec, flag->predidx[0], PRED_COUNT);
				ac_enc_bypass_NPOT(&ec, flag->predidx[1], PRED_COUNT);
				ac_enc_bypass_NPOT(&ec, flag->predidx[2], PRED_COUNT);

				//aux info
				{
					int usize=(image->iw*(arg->y2-arg->y1)*image->depth+7)>>3;
					int csize=(int)(arg->csizes[0]+arg->csizes[1]+arg->csizes[2]);
					int cidx[]=
					{
						comb[0]*PRED_COUNT+flag->predidx[0],
						comb[1]*PRED_COUNT+flag->predidx[1],
						comb[2]*PRED_COUNT+flag->predidx[2],
					};
					g_rct_usizes[arg->bestrct]+=csize;
					g_usizes[cidx[0]]+=usize;
					g_usizes[cidx[1]]+=usize;
					g_usizes[cidx[2]]+=usize;
					g_csizes[cidx[0]]+=arg->csizes[cidx[0]];
					g_csizes[cidx[1]]+=arg->csizes[cidx[1]];
					g_csizes[cidx[2]]+=arg->csizes[cidx[2]];
					if(loud)
					{
						printf("Y %5d~%5d  best %12.2lf bytes  %s [YUV: %s %s %s]\n",
							arg->y1, arg->y2,
							arg->bestsize,
							rct_names[flag->bestrct],
							pred_names[flag->predidx[0]],
							pred_names[flag->predidx[1]],
							pred_names[flag->predidx[2]]
						);

						for(int kp=0;kp<PRED_COUNT;++kp)
							printf("%14s", pred_names[kp]);
						printf("\n");
						for(int kc=0;kc<OCH_COUNT;++kc)
						{
							for(int kp=0;kp<PRED_COUNT;++kp)
								printf(" %12.2lf%c", arg->csizes[kc*PRED_COUNT+kp], kp==arg->predsel[kc]?'*':' ');
							printf("  %s\n", och_names[kc]);
						}

						for(int kt=0;kt<RCT_COUNT;++kt)
						{
							const int *group=rct_combinations[kt];
							double csize=
								arg->csizes[group[0]*PRED_COUNT+arg->predsel[group[0]]]+
								arg->csizes[group[1]*PRED_COUNT+arg->predsel[group[1]]]+
								arg->csizes[group[2]*PRED_COUNT+arg->predsel[group[2]]];
							printf("%12.2lf %c  %-10s %-10s %-10s %-10s\n",
								csize,
								kt==arg->bestrct?'*':' ',
								rct_names[kt],
								pred_names[arg->predsel[group[0]]],
								pred_names[arg->predsel[group[1]]],
								pred_names[arg->predsel[group[2]]]
							);
						}
#if 0
						if(!(kt+kt2))
							printf("block,  nrows,  usize,     best  ->  actual,  (actual-best)\n");
						printf(
							"[%3d]  %4d  %8d  %11.2lf -> %8zd bytes  (%+11.2lf)  %s %s %s %s\n",
							kt+kt2,
							arg->y2-arg->y1,
							csize,
							arg->bestsize,
							arg->list.nobj,
							(double)arg->list.nobj-arg->bestsize,
							rct_names[arg->bestrct],
							pred_names[arg->predidx[0]],
							pred_names[arg->predidx[1]],
							pred_names[arg->predidx[2]]
						);
#endif
						bestsize+=arg->bestsize;
					}
				}
			}
		}
		for(int k=0;k<nthreads;++k)
		{
			ThreadArgs *arg=args+k;
			_mm_free(arg->pixels);
			free(arg->hist);
		}
		free(args);
	}
	else//decode
	{
		for(int kc=0;kc<OCH_COUNT;++kc)
		{
			int depth=image->depth+och_info[kc][II_INFLATION], nlevels;
			UPDATE_MIN(depth, 16);
			depths[kc]=depth;
			nlevels=1<<depth;
			halfs[kc]=nlevels>>1;
		}
		ac_dec_init(&ec, cbuf, cbuf+clen);
		for(int kb=0;kb<nblocks;++kb)
		{
			DecorrInfo *flag=flags+kb;
			flag->bestrct=ac_dec_bypass_NPOT(&ec, RCT_COUNT);
			flag->predidx[0]=ac_dec_bypass_NPOT(&ec, PRED_COUNT);
			flag->predidx[1]=ac_dec_bypass_NPOT(&ec, PRED_COUNT);
			flag->predidx[2]=ac_dec_bypass_NPOT(&ec, PRED_COUNT);
		}
	}

	{
		int nlevels=2<<image->depth;//chroma-inflated for token
		int token=0, bypass=0, nbits=0;

		quantize_pixel(nlevels, &token, &bypass, &nbits);
		tlevels=token+1;
		clevels=9;
		statssize=clevels*clevels*(tlevels+1)*image->nch*(int)sizeof(int);
	}
	nctx=clevels*clevels;
	cdfstride=tlevels+1;
	chsize=nctx*cdfstride;
	
//#define PIXELSTRIDE 4*(1+WG_NPREDS)
#define PIXELSTRIDE 4*1
	bufsize=sizeof(short[4*PIXELSTRIDE])*(image->iw+16LL);//4 padded rows * max 4 channels * {pixels, pred errors...}
	pixels=(short*)_mm_malloc(bufsize, sizeof(__m128i));
	stats=(unsigned*)malloc(statssize);
	hist=(int*)malloc(statssize);
	if(!pixels||!hist||!stats)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	for(int kc=0;kc<image->nch;++kc)
	{
		int *curr_hist=hist+chsize*kc;
		unsigned *curr_CDF=stats+chsize*kc;
		
		*curr_hist=1;
		memfill(curr_hist+1, curr_hist, sizeof(int)*(tlevels-1LL), sizeof(int));
		curr_hist[tlevels]=tlevels;
		update_CDF(curr_hist, curr_CDF, tlevels);
		memfill(curr_hist+cdfstride, curr_hist, ((size_t)chsize-cdfstride)*sizeof(int), cdfstride*sizeof(int));
		memfill(curr_CDF+cdfstride, curr_CDF, ((size_t)chsize-cdfstride)*sizeof(int), cdfstride*sizeof(int));
	}
	for(int kc=0;kc<image->nch;++kc)
		wg_init(wg_weights+WG_NPREDS*kc);
	memset(pixels, 0, bufsize);
	for(int ky=0, idx=0;ky<image->ih;++ky)//codec loop
	{
		ALIGN(16) short *rows[]=
		{
			pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*PIXELSTRIDE,
			pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*PIXELSTRIDE,
			pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*PIXELSTRIDE,
			pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*PIXELSTRIDE,
		};
		short yuv[5]={0};
		int wg_preds[WG_NPREDS]={0};
		int token=0, bypass=0, nbits=0;
		int kb=ky/BLOCKSIZE;
		DecorrInfo *flag=flags+kb;
		const int *combination=rct_combinations[flag->bestrct];

		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			short
				*NNN	=rows[3]+0*PIXELSTRIDE,
				*NNW	=rows[2]-1*PIXELSTRIDE,
				*NN	=rows[2]+0*PIXELSTRIDE,
				*NNE	=rows[2]+1*PIXELSTRIDE,
				*NNEE	=rows[2]+2*PIXELSTRIDE,
				*NW	=rows[1]-1*PIXELSTRIDE,
				*N	=rows[1]+0*PIXELSTRIDE,
				*NE	=rows[1]+1*PIXELSTRIDE,
				*NEE	=rows[1]+2*PIXELSTRIDE,
				*NEEE	=rows[1]+3*PIXELSTRIDE,
				*WWW	=rows[0]-3*PIXELSTRIDE,
				*WW	=rows[0]-2*PIXELSTRIDE,
				*W	=rows[0]-1*PIXELSTRIDE,
				*curr	=rows[0]+0*PIXELSTRIDE;
			//           NNN
			//       NNW NN  NNE NNEE
			//       NW  N   NE  NEE  NEEE
			//WWW WW W  [?]
			if(ky<=2)
			{
				if(ky<=1)
				{
					if(ky==0)
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
				if(kx>=image->iw-2)
				{
					if(kx>=image->iw-1)
					{
						NNE=NN;
						NE=N;
					}
					NEE=NE;
				}
				NEEE=NEE;
			}
			if(fwd)
			{
				int kc=0;
				if(image->nch>=3)
				{
					switch(flag->bestrct)
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
						yuv[0]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT1_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT1_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT1_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=yuv[1]>>1;
						break;
					case RCT_RCT2_120:
					case RCT_RCT2_201:
					case RCT_RCT2_012:
					case RCT_RCT2_210:
					case RCT_RCT2_102:
					case RCT_RCT2_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT2_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT2_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT2_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=(2*yuv[1]-yuv[2]+4)>>3;//A710 from GWFX
						break;
					case RCT_RCT3_120:
					case RCT_RCT3_201:
					case RCT_RCT3_012:
					case RCT_RCT3_210:
					case RCT_RCT3_102:
					case RCT_RCT3_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT3_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT3_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT3_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=(2*yuv[1]+yuv[2]+4)>>3;
						break;
					case RCT_RCT4_120:
					case RCT_RCT4_201:
					case RCT_RCT4_012:
					case RCT_RCT4_210:
					case RCT_RCT4_102:
					case RCT_RCT4_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT4_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT4_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT4_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=yuv[1]/3;
						break;
					case RCT_RCT5_120:
					case RCT_RCT5_201:
					case RCT_RCT5_012:
					case RCT_RCT5_210:
					case RCT_RCT5_102:
					case RCT_RCT5_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT5_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT5_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT5_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=(3*yuv[1]+4)>>3;
						break;
					case RCT_RCT6_120:
					case RCT_RCT6_201:
					case RCT_RCT6_012:
					case RCT_RCT6_210:
					case RCT_RCT6_102:
					case RCT_RCT6_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT6_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT6_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT6_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=(7*yuv[1]+8)>>4;
						break;
					case RCT_RCT7_120:
					case RCT_RCT7_201:
					case RCT_RCT7_012:
					case RCT_RCT7_210:
					case RCT_RCT7_102:
					case RCT_RCT7_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT7_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT7_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT7_120][2]];
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=(10*yuv[1]-yuv[2]+16)>>5;
						break;
					case RCT_Pei09_120:
					case RCT_Pei09_201:
					case RCT_Pei09_012:
					case RCT_Pei09_210:
					case RCT_Pei09_102:
					case RCT_Pei09_021:
						yuv[0]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_Pei09_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_Pei09_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_Pei09_120][2]];
						yuv[1]-=(87*yuv[2]+169*yuv[0]+128)>>8;
						yuv[2]-=yuv[0];
						yuv[0]+=(86*yuv[2]+29*yuv[1]+128)>>8;
						break;
					case RCT_J2K_120:
					case RCT_J2K_201:
					case RCT_J2K_012:
						yuv[0]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_J2K_120][0]];
						yuv[1]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_J2K_120][1]];
						yuv[2]=image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_J2K_120][2]];
						yuv[1]-=yuv[0];
						yuv[2]-=yuv[0];
						yuv[0]+=(yuv[1]+yuv[2])>>2;
						break;
					}
					kc=3;
				}
				if(kc<image->nch)
					yuv[kc]=image->data[kc];
			}
			for(int kc=0;kc<image->nch;++kc)
			{
				int ch=combination[kc];
				int kc2=kc;
				//int kc2=kc*(1+WG_NPREDS);
				int offset=(yuv[combination[kc+3]]+yuv[combination[kc+6]])>>och_info[ch][II_HSHIFT];
				int pred=0, error, sym;
				int
					vx=(abs(W[kc2]-WW[kc2])+abs(N[kc2]-NW[kc2])+abs(NE[kc2]-N[kc2]))<<8>>depths[ch],
					vy=(abs(W[kc2]-NW[kc2])+abs(N[kc2]-NN[kc2])+abs(NE[kc2]-NNE[kc2]))<<8>>depths[ch];
				int qeN=FLOOR_LOG2_P1(vy);
				int qeW=FLOOR_LOG2_P1(vx);
				int cidx=cdfstride*(nctx*kc+clevels*MINVAR(qeN, 8)+MINVAR(qeW, 8));
				int *curr_hist=hist+cidx;
				unsigned *curr_CDF=stats+cidx;

				switch(flag->predidx[kc])
				{
				case PRED_W:
					pred=W[kc2];
					break;
				case PRED_cgrad:
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					break;
				case PRED_wgrad:
					//for(int kp=0;kp<WG_NPREDS;++kp)
					//{
					//	//       [NNW]   [NN]   [NNE] [NNEE]
					//	//NWW  2*[NW]  4*[N]  2*[NE]  [NEE]
					//	//WW   3* W       ?
					//	wg_errors[kp]=
					//		NNW	[kc2+kp+1]+//NNW+NWW
					//		NN	[kc2+kp+1]+//NN+NW
					//		NNE	[kc2+kp+1]+//NNE+N
					//		NNEE	[kc2+kp+1]+//NNEE+NE
					//		NW	[kc2+kp+1]+//NW+WW
					//		N	[kc2+kp+1]*3+//N+W
					//		NE	[kc2+kp+1]+//NE
					//		NEE	[kc2+kp+1];//NEE
					//}
					pred=wg_predict(
						wg_weights+WG_NPREDS*kc,
						NNN[kc2], NN[kc2], NNE[kc2],
						NW[kc2], N[kc2], NE[kc2], NEE[kc2], NEEE[kc2],
						WWW[kc2], WW[kc2], W[kc2],
						wg_errors+WG_NPREDS*kc, wg_preds
					);
					break;
				}
				pred+=offset;
				CLAMP2_32(pred, pred, -halfs[ch], halfs[ch]-1);
				if(fwd)
				{
					curr[kc2]=yuv[kc];
					error=yuv[kc]-pred;
					{
						int upred=halfs[ch]-abs(pred), aval=abs(error);
						if(aval<=upred)
						{
							sym=error;
#ifdef ENABLE_BIASCORR
							{
								int negmask=-((ibias_corr<0)&(sym!=-halfs[ch]));//sign is flipped if SSE correction was negative, to skew the histogram
								sym^=negmask;
								sym-=negmask;
							}
#endif
							sym=sym<<1^(sym>>31);//pack sign
						}
						else
							sym=upred+aval;//error sign is known
					}
					quantize_pixel(sym, &token, &bypass, &nbits);
					ac_enc(&ec, token, curr_CDF);
					if(nbits)
						ac_enc_bypass(&ec, bypass, nbits);//up to 16 bits
				}
				else
				{
					token=ac_dec(&ec, curr_CDF, tlevels);//try ac_dec_packedsign()
					sym=token;
					if(sym>=(1<<CONFIG_EXP))
					{
						sym-=1<<CONFIG_EXP;
						int lsb=sym&((1<<CONFIG_LSB)-1);
						sym>>=CONFIG_LSB;
						int msb=sym&((1<<CONFIG_MSB)-1);
						sym>>=CONFIG_MSB;
						nbits=sym+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB);
						bypass=ac_dec_bypass(&ec, nbits);
						sym=1;
						sym<<=CONFIG_MSB;
						sym|=msb;
						sym<<=nbits;
						sym|=bypass;
						sym<<=CONFIG_LSB;
						sym|=lsb;
					}
					{
						int upred=halfs[ch]-abs(pred), negmask=0;
						if(sym<=(upred<<1))
						{
							error=sym>>1^-(sym&1);
#ifdef ENABLE_BIASCORR
							negmask=-((ibias_corr<0)&(error!=-halfs[ch]));
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
					curr[kc2]=yuv[kc]=error+pred;
				}
				++curr_hist[token];
				++curr_hist[tlevels];
				if(curr_hist[tlevels]>=tlevels&&!(curr_hist[tlevels]&(CDFSTRIDE-1)))
					update_CDF(curr_hist, curr_CDF, tlevels);
				if(curr_hist[tlevels]>=6144)
				{
					int sum=0;
					for(int ks=0;ks<tlevels;++ks)
						sum+=curr_hist[ks]>>=1;
					curr_hist[tlevels]=sum;
				}
				curr[kc2]-=offset;
				if(flag->predidx[kc]==PRED_wgrad)
				{
					wg_update(curr[kc2], wg_preds, wg_errors+WG_NPREDS*kc, wg_weights+WG_NPREDS*kc);
					//for(int kp=0;kp<WG_NPREDS;++kp)
					//{
					//	//       [NNW]   [NN]   [NNE] [NNEE]
					//	//NWW  2*[NW]  4*[N]  2*[NE]  [NEE]
					//	//WW   3* W       ?
					//	curr[kc2+kp+1]=(curr[kc2+kp+1]+3*wg_errors[kp])>>2;
					//	NE[kc2+kp+1]+=wg_errors[kp];
					//}
				}
			}
			if(!fwd)
			{
				Image *image=dst;
				int kc=0;
				if(image->nch>=3)
				{
					switch(flag->bestrct)
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
						yuv[0]-=yuv[1]>>1;
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT1_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT1_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT1_120][2]]=yuv[2];
						break;
					case RCT_RCT2_120:
					case RCT_RCT2_201:
					case RCT_RCT2_012:
					case RCT_RCT2_210:
					case RCT_RCT2_102:
					case RCT_RCT2_021:
						yuv[0]-=(2*yuv[1]-yuv[2]+4)>>3;//A710
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT2_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT2_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT2_120][2]]=yuv[2];
						break;
					case RCT_RCT3_120:
					case RCT_RCT3_201:
					case RCT_RCT3_012:
					case RCT_RCT3_210:
					case RCT_RCT3_102:
					case RCT_RCT3_021:
						yuv[0]-=(2*yuv[1]+yuv[2]+4)>>3;
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT3_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT3_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT3_120][2]]=yuv[2];
						break;
					case RCT_RCT4_120:
					case RCT_RCT4_201:
					case RCT_RCT4_012:
					case RCT_RCT4_210:
					case RCT_RCT4_102:
					case RCT_RCT4_021:
						yuv[0]-=yuv[1]/3;
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT4_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT4_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT4_120][2]]=yuv[2];
						break;
					case RCT_RCT5_120:
					case RCT_RCT5_201:
					case RCT_RCT5_012:
					case RCT_RCT5_210:
					case RCT_RCT5_102:
					case RCT_RCT5_021:
						yuv[0]-=(3*yuv[1]+4)>>3;
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT5_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT5_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT5_120][2]]=yuv[2];
						break;
					case RCT_RCT6_120:
					case RCT_RCT6_201:
					case RCT_RCT6_012:
					case RCT_RCT6_210:
					case RCT_RCT6_102:
					case RCT_RCT6_021:
						yuv[0]-=(7*yuv[1]+8)>>4;
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT6_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT6_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT6_120][2]]=yuv[2];
						break;
					case RCT_RCT7_120:
					case RCT_RCT7_201:
					case RCT_RCT7_012:
					case RCT_RCT7_210:
					case RCT_RCT7_102:
					case RCT_RCT7_021:
						yuv[0]-=(10*yuv[1]-yuv[2]+16)>>5;
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT7_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT7_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_RCT7_120][2]]=yuv[2];
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
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_Pei09_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_Pei09_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_Pei09_120][2]]=yuv[2];
						break;
					case RCT_J2K_120:
					case RCT_J2K_201:
					case RCT_J2K_012:
						yuv[0]-=(yuv[1]+yuv[2])>>2;
						yuv[2]+=yuv[0];
						yuv[1]+=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_J2K_120][0]]=yuv[0];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_J2K_120][1]]=yuv[1];
						image->data[idx+rgb2yuv_permutations[flag->bestrct-RCT_J2K_120][2]]=yuv[2];
						break;
					}
					kc=3;
				}
				if(kc<image->nch)
					image->data[kc]=yuv[kc];
			}
#ifdef ENABLE_GUIDE
			if(memcmp(image->data+idx, guide->data+idx, sizeof(short)*image->nch))
			{
				short orig[4]={0};
				memcpy(orig, guide->data+idx, image->nch*sizeof(short));
				LOG_ERROR("Guide error XY %d %d", kx, ky);
				printf("");//
			}
#endif
			rows[0]+=PIXELSTRIDE;
			rows[1]+=PIXELSTRIDE;
			rows[2]+=PIXELSTRIDE;
			rows[3]+=PIXELSTRIDE;
		}
	}
	if(fwd)
	{
		ac_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
	}
	if(loud)
	{
		ptrdiff_t usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
		t0=time_sec()-t0;
		if(fwd)
		{
			ptrdiff_t csize=list.nobj;
			printf("Best %15.2lf (%+13.2lf) bytes\n", bestsize, csize-bestsize);
			printf("Size %12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
			printf("Mem usage: ");
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
		}
		printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
	}
	if(fwd)
		dlist_clear(&list);
	free(hist);
	free(stats);
	_mm_free(pixels);
	free(flags);
	return 0;
}
void f27_curiosity()
{
	double rowusum[OCH_COUNT]={0}, colusum[PRED_COUNT]={0};
	double rowcsum[OCH_COUNT]={0}, colcsum[PRED_COUNT]={0};
	double utotal=0, ctotal=0;
	for(int kc=0;kc<OCH_COUNT;++kc)
	{
		for(int kp=0;kp<PRED_COUNT;++kp)
		{
			rowusum[kc]+=g_usizes[kc*PRED_COUNT+kp];
			rowcsum[kc]+=g_csizes[kc*PRED_COUNT+kp];
			colusum[kp]+=g_usizes[kc*PRED_COUNT+kp];
			colcsum[kp]+=g_csizes[kc*PRED_COUNT+kp];
		}
	}
	for(int kp=0;kp<PRED_COUNT;++kp)
	{
		utotal+=colusum[kp];
		ctotal+=colcsum[kp];
	}

	printf("Output channel compressed sizes:\n");
	for(int kp=0;kp<PRED_COUNT;++kp)
		printf("%*s%*s", 15, pred_names[kp], 14, "");
	printf("\n");

	for(int kp=0;kp<PRED_COUNT;++kp)
	{
		double p=100.*colusum[kp]/utotal;
		printf("%*s", 9, "");
		print_nan(p, 10, 6);
		printf("%c%*s", p!=p?' ':'%', 10, "");
	}
	printf("  coverage\n");

	for(int kp=0;kp<PRED_COUNT;++kp)
	{
		double p=100.*colcsum[kp]/colusum[kp];
		printf("%*s", 9, "");
		print_nan(p, 10, 6);
		printf("%c%*s", p!=p?' ':'%', 10, "");
	}
	printf("  ratio\n");
	for(int kc=0;kc<OCH_COUNT;++kc)
	{
		double p1=100.*rowusum[kc]/utotal, p2=100.*rowcsum[kc]/rowusum[kc];
		for(int kp=0;kp<PRED_COUNT;++kp)
			printf(" %12.2lf->%12.2lf", g_usizes[kc*PRED_COUNT+kp], g_csizes[kc*PRED_COUNT+kp]);
		printf("  ");
		print_nan(p1, 10, 6);
		printf("%c ", p1!=p1?' ':'%');
		print_nan(p2, 10, 6);
		printf("%c  %s\n", p2!=p2?' ':'%', och_names[kc]);
	}

	printf("RCT contributions:\n");
	for(int kt=0;kt<RCT_COUNT;++kt)
		printf("%12.2lf  %8.4lf%%  %s\n", g_rct_usizes[kt], 100.*g_rct_usizes[kt]/utotal, rct_names[kt]);
}