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

//	#define ENABLE_CURIOSITY

//	#define USE_AC2//bad
//	#define MIXCDF//bad
//	#define WG_T47
//	#define WG_UPDATE//horrible


//#define PROB_BITS 24//X
//#define EMIT_BITS 1
#include"ac.h"
#ifdef USE_AC2
#define USE_PROB_BITS AC2_PROB_BITS
#else
#define USE_PROB_BITS PROB_BITS
#endif
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
#define BLOCKSIZE 256

#define ANALYSIS_STRIDE 3

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
#define WG_DECAY_NUM	493
#define WG_DECAY_SH	9

#if 1
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
#endif
#if 0
#define WG_NPREDS	13
#define WG_PREDLIST\
	WG_PRED(132, N+W-NW)\
	WG_PRED(176, N+W-NW+((eN+eW-eNW+16)>>5))\
	WG_PRED(120, N+eN)\
	WG_PRED( 90, N+(eN>>2))\
	WG_PRED( 70, N)\
	WG_PRED(120, W+eW)\
	WG_PRED( 90, W+(eW>>2))\
	WG_PRED( 70, W)\
	WG_PRED(165, W+NE-N)\
	WG_PRED(220, N+NE-NNE)\
	WG_PRED(165, 3*(N-NN)+NNN)\
	WG_PRED(165, 3*(W-WW)+WWW)\
	WG_PRED(176, (W+NEEE)/2)
//	WG_PRED(  0, N+NE-NNE+((eN+eNE-eNNE)>>2))
#endif
#if 0
#define WG_NPREDS	13
#define WG_PREDLIST\
	WG_PRED(132, N+W-NW)\
	WG_PRED(176, N+W-NW+((eN+eW-eNW+16)>>5))\
	WG_PRED(176, N+eN)\
	WG_PRED( 88, N)\
	WG_PRED(176, W+eW)\
	WG_PRED( 88, W)\
	WG_PRED(165, W+NE-N)\
	WG_PRED(200, N+NE-NNE)\
	WG_PRED( 10, NE+NEE-NNEEE)\
	WG_PRED( 10, 3*NE-(NEE+NEEE))\
	WG_PRED(165, 3*(N-NN)+NNN)\
	WG_PRED(165, 3*(W-WW)+WWW)\
	WG_PRED(176, (W+NEEE)/2)
//	WG_PRED(  0, N+NE-NNE+((eN+eNE-eNNE)>>2))
#endif
#if 0
#define WG_NPREDS	11
#define WG_PREDLIST\
	WG_PRED(0.75, N+W-NW)\
	WG_PRED(1.0, N+W-NW+((eN+eW-eNW+16)>>5))\
	WG_PRED(1.0, N+eN)\
	WG_PRED(1.0, W+eW)\
	WG_PRED(0.5, N)\
	WG_PRED(0.5, W)\
	WG_PRED(0.9375, W+NE-N)\
	WG_PRED(1.25, N+NE-NNE)\
	WG_PRED(0.9375, 3*(N-NN)+NNN)\
	WG_PRED(0.9375, 3*(W-WW)+WWW)\
	WG_PRED(1.0, (W+NEEE)/2)
//	WG_PRED(0.5, NW+(eNW>>2))//X
//	WG_PRED(0.5, NE+(eNE>>2))
#endif
#if 0
#define WG_NPREDS	8
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
#if 0
#define WG_NPREDS	8
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
#ifdef WG_T47
#define WG_NPREDS	32
#define WG_PREDLIST\
	WG_PRED(176.0, W+NE-N-((2*(eN+eW)+eNE-eNW+4)>>3))\
	WG_PRED(176.0, N-(int)(((long long)eN+eW+eNE)*-0x05C>>8))\
	WG_PRED(176.0, W-(int)(((long long)eN+eW+eNW)*-0x05B>>8))\
	WG_PRED(176.0, N+(int)((-eNN*0x0DFLL-eN*0x051LL-eNE*0x0BDLL+((long long)N-NN)*0x05C+((long long)NW-W)*0x102)>>8))\
	WG_PRED(176.0, 3*(N-NN)+NNN)\
	WG_PRED(176.0, (N+W)>>1)\
	WG_PRED(176.0, N+W-NW)\
	WG_PRED(176.0, (W+NEE)>>1)\
	WG_PRED(176.0, (3*W+NEEE)>>2)\
	WG_PRED(176.0, (3*(3*W+NE+NEE)-10*N+2)/5)\
	WG_PRED(176.0, (3*(3*W+NE+NEE)-10*N)/5)\
	WG_PRED(176.0, (4*N-2*NN+NW+NE)>>2)\
	WG_PRED(176.0, N+NE-NNE-eNNE)\
	WG_PRED(176.0, (4*(N+W+NW+NE)-(NN+WW+NNWW+NNEE)+6)/12)\
	WG_PRED(176.0, W+((eW-eWW)>>1))\
	WG_PRED(176.0, paper_GAP)\
	WG_PRED(176.0, calic_GAP)\
	WG_PRED(176.0, N+W-((NW+NN+WW+NE)>>2))\
	WG_PRED(176.0, ((2*(N+W)-(NW+NN+WW+NE))*9+(WWW+NWW+NNW+NNN+NNE+NEE)*2)/12)\
	WG_PRED(176.0, 3*(N+W-NW-(NN+WW-NNWW))+NNN+WWW-NNNWWW)\
	WG_PRED(176.0, 2*(W+NE-N)-(WW+NNEE-NN))\
	WG_PRED(176.0, (2*W+NEE-N)>>1)\
	WG_PRED(176.0, NW+NWW-NNWWW)\
	WG_PRED(176.0, (14*NE-(NNEE+NNNEE+NNEEE))/11)\
	WG_PRED(176.0, (NEEE+NEEEE)>>1)\
	WG_PRED(176.0, (NNNEEEE+NNEEE)>>1)\
	WG_PRED(176.0, NNEEEE)\
	WG_PRED(176.0, (NNWWWW+NNNWWWW)>>1)\
	WG_PRED(176.0, (WWW+WWWW)>>1)\
	WG_PRED(176.0, (N+NN)>>1)\
	WG_PRED(176.0, (NE+NNEE)>>1)\
	WG_PRED(176.0, (NE+NNE+NEE+NNEE)>>2)
#endif
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
#ifdef WG_T47
	short
		NNNWWWW	=rows[3][kc2-4*stride+0],
		NNNWWW	=rows[3][kc2-3*stride+0],
		NNN	=rows[3][kc2+0*stride+0],
		NNNEE	=rows[3][kc2+2*stride+0],
		NNNEEEE	=rows[3][kc2+4*stride+0],
		NNWWWW	=rows[2][kc2-4*stride+0],
		NNWWW	=rows[2][kc2-3*stride+0],
		NNWW	=rows[2][kc2-2*stride+0],
		NNW	=rows[2][kc2-1*stride+0],
		NN	=rows[2][kc2+0*stride+0],
		NNE	=rows[2][kc2+1*stride+0],
		NNEE	=rows[2][kc2+2*stride+0],
		NNEEE	=rows[2][kc2+3*stride+0],
		NNEEEE	=rows[2][kc2+4*stride+0],
		NWW	=rows[1][kc2-2*stride+0],
		NW	=rows[1][kc2-1*stride+0],
		N	=rows[1][kc2+0*stride+0],
		NE	=rows[1][kc2+1*stride+0],
		NEE	=rows[1][kc2+2*stride+0],
		NEEE	=rows[1][kc2+3*stride+0],
		NEEEE	=rows[1][kc2+4*stride+0],
		WWWW	=rows[0][kc2-4*stride+0],
		WWW	=rows[0][kc2-3*stride+0],
		WW	=rows[0][kc2-2*stride+0],
		W	=rows[0][kc2-1*stride+0],
		eNN	=rows[2][kc2+0*stride+1],
		eNNE	=rows[2][kc2+1*stride+1],
		eNW	=rows[1][kc2-1*stride+1],
		eN	=rows[1][kc2+0*stride+1],
		eNE	=rows[1][kc2+1*stride+1],
		eWW	=rows[0][kc2-2*stride+1],
		eW	=rows[0][kc2-1*stride+1];
	int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
	int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
	int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
	int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
	int diff=(dy-dx)<<8>>depth, diff2=(d45-d135)<<8>>depth, diff3=NE-NW;
	int paper_GAP, calic_GAP;

	if(dy+dx>32)//[sic]
		paper_GAP=(int)(((long long)dx*N+(long long)dy*W)/((long long)dy+dx));
	else if(diff>12)
		paper_GAP=(N+2*W)/3;
	else if(diff<-12)
		paper_GAP=(2*N+W)/3;
	else
		paper_GAP=(N+W)>>1;

	if(diff2>32)
		paper_GAP+=diff3>>2;
	else if(diff2>16)
		paper_GAP+=diff3*3>>4;
	else if(diff2>=-16)
		paper_GAP+=diff3>>3;
	else if(diff2>=-32)
		paper_GAP+=diff3>>4;

	if(diff>80)
		calic_GAP=W;
	else if(diff<-80)
		calic_GAP=N;
	else if(diff>32)
		calic_GAP=(2*N+6*W+NE-NW)>>3;		//c1	[1/4  3/4  1/8  -1/8].[N W NE NW]
	else if(diff>8)
		calic_GAP=(6*N+10*W+3*(NE-NW))>>4;	//c2	[3/8  5/8  3/16  -3/16]
	else if(diff<-32)
		calic_GAP=(6*N+2*W+NE-NW)>>3;		//c3	[3/4  1/4  1/8  -1/8]
	else if(diff<-8)
		calic_GAP=(10*N+6*W+3*(NE-NW))>>4;	//c4	[5/8  3/8  3/16  -3/16]
	else
		calic_GAP=(((N+W)<<1)+NE-NW)>>2;	//c5	[1/2  1/2  1/4  -1/4]
#else
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
#endif

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
#ifdef WG_UPDATE
	double wsum=0;
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
		
		//{
		//	//https://stackoverflow.com/questions/11644441/fast-inverse-square-root-on-x64
		//	double t1=e2+1, t2=t1*0.5;
		//	long long t3=*(long long*)&t1;
		//	t3=0x5FE6EB50C7B537A9-(t3>>1);
		//	t1=*(double*)&t3;
		//	t1=t1*(1.5-(t2*t1*t1));
		//	wsum+=weights[k]+=t1;
		//}
		//wsum+=weights[k]+=1./sqrt(e2+1);
		//weights[k]+=1./(e2+1);
#endif
	}
#ifdef WG_UPDATE
	//wsum=1/wsum;
	//for(int k=0;k<WG_NPREDS;++k)
	//	weights[k]*=wsum;

	//weights[kbest]+=1./256;
	++weights[kbest];
	//if(weights[kbest]>22)//352/16.
	if(weights[kbest]>352)
	{
		for(int k=0;k<WG_NPREDS;++k)
			weights[k]>>=1;
			//weights[k]*=0.5;
	}
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
		*nbits=lgv-(CONFIG_MSB+CONFIG_LSB);
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
	const Image *src;
	Image *dst;
	int fwd, loud, y1, y2;
	
	short *pixels;
	int bufsize;

	int *hist, histsize, histindices[OCH_COUNT*PRED_COUNT+1];

	//AC
	int tlevels, clevels, statssize;
	unsigned *stats;
	DList list;
	const unsigned char *decstart, *decend;

	//aux
	int bestrct, predidx[3];
	double usizes[3], csizes[3], bestsize;

	//WG
	int wg_weights[OCH_COUNT*WG_NPREDS];
	int wg_errors[OCH_COUNT*WG_NPREDS];
} ThreadArgs;
#ifndef MIXCDF
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
#ifdef USE_AC2
	AC2 ec;
#else
	ArithmeticCoder ec;
#endif
	ThreadArgs *args=(ThreadArgs*)param;
	Image const *image=args->fwd?args->src:args->dst;
	int depths[OCH_COUNT]={0}, halfs[OCH_COUNT]={0};
	double bestsize=0;
	int bestrct=0, combination[12]={0}, predidx[4]={0}, flag=0;
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
		int res=(image->iw/ANALYSIS_STRIDE)*((args->y2-args->y1)/ANALYSIS_STRIDE);
		int nlevels[OCH_COUNT]={0};
		double csizes[OCH_COUNT*PRED_COUNT]={0};
		int predsel[OCH_COUNT]={0};

		for(int kc=0;kc<OCH_COUNT;++kc)
			nlevels[kc]=1<<depths[kc];

		for(int kc=0;kc<OCH_COUNT;++kc)
			wg_init(args->wg_weights+WG_NPREDS*kc);
		memset(args->wg_errors, 0, sizeof(args->wg_errors));

		memset(args->hist, 0, args->histsize);
		memset(args->pixels, 0, args->bufsize);
		for(int ky=args->y1, ay=0;ky<args->y2-3;ky+=ANALYSIS_STRIDE, ++ay)
		{
			ALIGN(16) short *rows[]=
			{
				args->pixels+((image->iw+16LL)*((ay-0LL)&3)+8LL)*OCH_COUNT*2,
				args->pixels+((image->iw+16LL)*((ay-1LL)&3)+8LL)*OCH_COUNT*2,
				args->pixels+((image->iw+16LL)*((ay-2LL)&3)+8LL)*OCH_COUNT*2,
				args->pixels+((image->iw+16LL)*((ay-3LL)&3)+8LL)*OCH_COUNT*2,
			};
			int idx=image->nch*image->iw*ky;
			short input[ICH_COUNT]={0};
			int preds[PRED_COUNT]={0};
			int wg_preds[WG_NPREDS]={0};
			for(int kx=0;kx<image->iw;kx+=ANALYSIS_STRIDE, idx+=image->nch*ANALYSIS_STRIDE)
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

					preds[PRED_W]=W[kc2];
					MEDIAN3_32(preds[PRED_cgrad], N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					preds[PRED_wgrad]=wg_predict(
						args->wg_weights+WG_NPREDS*kc,
						rows, OCH_COUNT*2, kc2, depths[kc],
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
					curr[kc2+1]=target-preds[PRED_wgrad];
					target-=offset;
					wg_update(target, wg_preds, args->wg_errors+WG_NPREDS*kc, args->wg_weights+WG_NPREDS*kc);
					curr[kc2+0]=target;
				}
				rows[0]+=OCH_COUNT*2;
				rows[1]+=OCH_COUNT*2;
				rows[2]+=OCH_COUNT*2;
				rows[3]+=OCH_COUNT*2;
			}
		}
		for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)//calculate channel sizes
		{
			int *curr_hist=args->hist+args->histindices[kc];
			int nlevels2=nlevels[kc/PRED_COUNT];
			for(int ks=0;ks<nlevels2;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
					csizes[kc]-=freq*log2((double)freq/res);
			}
			csizes[kc]/=8;
		}
		for(int kc=0;kc<OCH_COUNT;++kc)//select best predictors
		{
			int bestpred=0;
			for(int kp=1;kp<PRED_COUNT;++kp)
			{
				if(csizes[kc*PRED_COUNT+bestpred]>csizes[kc*PRED_COUNT+kp])
					bestpred=kp;
			}
			predsel[kc]=bestpred;
		}
		for(int kt=0;kt<RCT_COUNT;++kt)//select best R.C.T
		{
			const int *group=rct_combinations[kt];
			double csize=
				csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
				csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
				csizes[group[2]*PRED_COUNT+predsel[group[2]]];
			if(!kt||bestsize>csize)
				bestsize=csize, bestrct=kt;
		}
		memcpy(combination, rct_combinations[bestrct], sizeof(combination));
		predidx[0]=predsel[combination[0]];
		predidx[1]=predsel[combination[1]];
		predidx[2]=predsel[combination[2]];
		flag=bestrct;
		flag=flag*PRED_COUNT+predidx[2];
		flag=flag*PRED_COUNT+predidx[1];
		flag=flag*PRED_COUNT+predidx[0];//19*3*3*3 = 513

		args->bestrct=bestrct;
		args->predidx[0]=predidx[0];
		args->predidx[1]=predidx[1];
		args->predidx[2]=predidx[2];
		args->usizes[0]=(res*image->depth+7)>>3;
		args->usizes[1]=(res*image->depth+7)>>3;
		args->usizes[2]=(res*image->depth+7)>>3;
		{
			const int *group=rct_combinations[bestrct];
			args->csizes[0]=csizes[group[0]*PRED_COUNT+predsel[group[0]]];
			args->csizes[1]=csizes[group[1]*PRED_COUNT+predsel[group[1]]];
			args->csizes[2]=csizes[group[2]*PRED_COUNT+predsel[group[2]]];
		}
		args->bestsize=bestsize;
		if(args->loud)
		{
			printf("Y %5d~%5d  best %12.2lf bytes  %s [YUV: %s %s %s]\n",
				args->y1, args->y2,
				bestsize,
				rct_names[bestrct],
				pred_names[predidx[0]],
				pred_names[predidx[1]],
				pred_names[predidx[2]]
			);

			for(int kp=0;kp<PRED_COUNT;++kp)
				printf("%14s", pred_names[kp]);
			printf("\n");
			for(int kc=0;kc<OCH_COUNT;++kc)
			{
				for(int kp=0;kp<PRED_COUNT;++kp)
					printf(" %12.2lf%c", csizes[kc*PRED_COUNT+kp], kp==predsel[kc]?'*':' ');
				printf("  %s\n", och_names[kc]);
			}

			for(int kt=0;kt<RCT_COUNT;++kt)
			{
				const int *group=rct_combinations[kt];
				double csize=
					csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
					csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
					csizes[group[2]*PRED_COUNT+predsel[group[2]]];
				printf("%12.2lf %c  %-10s %-10s %-10s %-10s\n",
					csize,
					kt==bestrct?'*':' ',
					rct_names[kt],
					pred_names[predsel[group[0]]],
					pred_names[predsel[group[1]]],
					pred_names[predsel[group[2]]]
				);
			}
		}
		dlist_init(&args->list, 1, 1024, 0);
		dlist_push_back(&args->list, &flag, sizeof(char[2]));
#ifdef USE_AC2
		ac2_enc_init(&ec, &args->list);
#else
		ac_enc_init(&ec, &args->list);
#endif
	}
	else//decode
	{
		const unsigned char *srcstart=args->decstart, *srcend=args->decend;
			
		memcpy(&flag, srcstart, sizeof(char[2]));
		srcstart+=sizeof(char[2]);
		bestrct=flag;
		predidx[0]=bestrct%PRED_COUNT;	bestrct/=PRED_COUNT;
		predidx[1]=bestrct%PRED_COUNT;	bestrct/=PRED_COUNT;
		predidx[2]=bestrct%PRED_COUNT;	bestrct/=PRED_COUNT;
		if((unsigned)bestrct>=(unsigned)RCT_COUNT)
		{
			LOG_ERROR("Corrupt file");
			return;
		}
		memcpy(combination, rct_combinations[bestrct], sizeof(combination));
#ifdef USE_AC2
		ac2_dec_init(&ec, srcstart, srcend);
#else
		ac_dec_init(&ec, srcstart, srcend);
#endif
	}
	predidx[3]=PRED_cgrad;
#ifdef MIXCDF
	for(int ks=0;ks<=args->tlevels;++ks)//init bypass
		args->stats[ks]=(unsigned)(((1LL<<USE_PROB_BITS)-args->tlevels)*ks/args->tlevels);
	memfill(args->stats+cdfstride, args->stats, sizeof(int)*((size_t)chsize*image->nch-cdfstride), sizeof(int)*cdfstride);
#else
	for(int kc=0;kc<image->nch;++kc)
	{
		int *curr_hist=args->hist+chsize*kc;
		unsigned *curr_CDF=args->stats+chsize*kc;
		
		*curr_hist=1;//init bypass
		memfill(curr_hist+1, curr_hist, sizeof(int)*(args->tlevels-1LL), sizeof(int));
		curr_hist[args->tlevels]=args->tlevels;
		update_CDF(curr_hist, curr_CDF, args->tlevels);
		memfill(curr_hist+cdfstride, curr_hist, ((size_t)chsize-cdfstride)*sizeof(int), cdfstride*sizeof(int));
		memfill(curr_CDF+cdfstride, curr_CDF, ((size_t)chsize-cdfstride)*sizeof(int), cdfstride*sizeof(int));
	}
#endif
	for(int kc=0;kc<3;++kc)
	{
		if(predidx[kc]==PRED_wgrad)
			wg_init(args->wg_weights+WG_NPREDS*kc);
	}
	memset(args->wg_errors, 0, sizeof(args->wg_errors));
	memset(args->pixels, 0, args->bufsize);
	for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) short *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		short yuv[5]={0};
		int wg_preds[WG_NPREDS]={0};
		int token=0, bypass=0, nbits=0;

		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			short
				*NNN	=rows[3]+0*4*2,
				*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NNEE	=rows[2]+2*4*2,
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
			for(int kc=0;kc<image->nch;++kc)
			{
				int ch=combination[kc];
				int kc2=kc<<1;
				int offset=(yuv[combination[kc+4]]+yuv[combination[kc+8]])>>och_info[ch][II_HSHIFT];
				int pred=0, error, sym;
				int
					vx=(abs(W[kc2]-WW[kc2])+abs(N[kc2]-NW[kc2])+abs(NE[kc2]-N  [kc2])+abs(WWW[kc2+1])+abs(WW[kc2+1])+abs(W[kc2+1])*2)<<8>>depths[ch],
					vy=(abs(W[kc2]-NW[kc2])+abs(N[kc2]-NN[kc2])+abs(NE[kc2]-NNE[kc2])+abs(NNN[kc2+1])+abs(NN[kc2+1])+abs(N[kc2+1])*2)<<8>>depths[ch];
				int qeN=FLOOR_LOG2(vy+1);
				int qeW=FLOOR_LOG2(vx+1);
				int cidx=cdfstride*(nctx*kc+args->clevels*MINVAR(qeN, 8)+MINVAR(qeW, 8));
#ifndef MIXCDF
				int *curr_hist=args->hist+cidx;
#endif
				unsigned *curr_CDF=args->stats+cidx;

				switch(predidx[kc])
				{
				case PRED_W:
					pred=W[kc2];
					break;
				case PRED_cgrad:
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					break;
				case PRED_wgrad:
					pred=wg_predict(
						args->wg_weights+WG_NPREDS*kc,
						rows, 4*2, kc2, depths[kc],
						args->wg_errors+WG_NPREDS*kc, wg_preds
					);
					break;
				}
				pred+=offset;
				CLAMP2_32(pred, pred, -halfs[ch], halfs[ch]-1);

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
				//	printf("");

				if(args->fwd)
				{
					curr[kc2+0]=yuv[kc];
					curr[kc2+1]=error=yuv[kc]-pred;
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
							sym=sym<<1^sym>>31;//pack sign
						}
						else
							sym=upred+aval;//error sign is known
					}
					quantize_pixel(sym, &token, &bypass, &nbits);
#ifdef MIXCDF
					ac2_enc_update(&ec, curr_CDF[token]+token, curr_CDF[token+1]+token+1);
#elif defined USE_AC2
					ac2_enc_update(&ec, curr_CDF[token], curr_CDF[token+1]);
#else
					ac_enc(&ec, token, curr_CDF);
#endif
					if(nbits)
#ifdef USE_AC2
						ac2_enc_bypass(&ec, bypass, nbits);
#else
						ac_enc_bypass(&ec, bypass, nbits);//up to 16 bits
#endif
				}
				else
				{
#ifdef MIXCDF
					unsigned cdf=ac2_dec_getcdf(&ec);
					int range=args->tlevels;

					token=0;
					while(range)//binary search
					{
						int floorhalf=range>>1, ks=token+floorhalf+1;
						if(cdf>=curr_CDF[ks]+ks)
							token+=range-floorhalf;
						range=floorhalf;
					}
					//if(token>=args->tlevels)
					//	LOG_ERROR("");
					ac2_dec_update(&ec, curr_CDF[token]+token, curr_CDF[token+1]+token+1);
#elif defined USE_AC2
					unsigned cdf=ac2_dec_getcdf(&ec);
					int range=args->tlevels;

					token=0;
					while(range)//binary search
					{
						int floorhalf=range>>1, ks=token+floorhalf+1;
						if(cdf>=curr_CDF[ks])
							token+=range-floorhalf;
						range=floorhalf;
					}
					//if(token>=args->tlevels)
					//	LOG_ERROR("");
					ac2_dec_update(&ec, curr_CDF[token], curr_CDF[token+1]);
#else
					token=ac_dec(&ec, curr_CDF, args->tlevels);//try ac_dec_packedsign()
#endif
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
#ifdef USE_AC2
						bypass=ac2_dec_bypass(&ec, nbits);
#else
						bypass=ac_dec_bypass(&ec, nbits);
#endif
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
					curr[kc2+0]=yuv[kc]=error+pred;
					curr[kc2+1]=error;
				}
#ifdef MIXCDF
				for(int ks=0, mixin=0, mag=(1<<USE_PROB_BITS)-args->tlevels;ks<args->tlevels;++ks)
				{
					mixin|=mag&(token-ks)>>31;
					curr_CDF[ks]+=(int)(mixin-curr_CDF[ks])>>8;
				}
				//for(int ks=0;ks<args->tlevels;++ks)
				//	curr_CDF[ks]+=((((1<<USE_PROB_BITS)-args->tlevels)&-(ks>token))+ks-curr_CDF[ks])>>8;

				//for(int ks=0;ks<args->tlevels-1;++ks)//
				//{
				//	if(curr_CDF[ks]>curr_CDF[ks+1])
				//		LOG_ERROR("");
				//}
#else
				++curr_hist[token];
				++curr_hist[args->tlevels];
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
				curr[kc2]-=offset;
				if(predidx[kc]==PRED_wgrad)
					wg_update(curr[kc2], wg_preds, args->wg_errors+WG_NPREDS*kc, args->wg_weights+WG_NPREDS*kc);
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
#ifdef USE_AC2
		ac2_enc_flush(&ec);
#else
		ac_enc_flush(&ec);
#endif
		if(args->loud)
			printf("Actual %8zd bytes (%+11.2lf)\n\n", args->list.nobj, args->list.nobj-bestsize);
	}
}
int f26_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	int ncores=query_cpu_cores();
	int nblocks=(image->ih+BLOCKSIZE-1)/BLOCKSIZE, nthreads=MINVAR(nblocks, ncores);
	int coffset=sizeof(int)*nblocks;
	ptrdiff_t start=0;
	ptrdiff_t memusage=0;
	ptrdiff_t argssize=nthreads*sizeof(ThreadArgs);
	ThreadArgs *args=(ThreadArgs*)malloc(argssize);
	int histindices[OCH_COUNT*PRED_COUNT+1]={0}, histsize=0;
	int tlevels=0, clevels=0, statssize=0;
	double bestsize=0;
	
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
			for(int kp=0;kp<PRED_COUNT;++kp)
			{
				int nlevels=1<<depth;
				histindices[kc*PRED_COUNT+kp]=histsize;
				histsize+=nlevels;
			}
		}
		histindices[OCH_COUNT*PRED_COUNT]=histsize;
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
		arg->bufsize=sizeof(short[4*OCH_COUNT*2])*(image->iw+16LL);//4 padded rows * OCH_COUNT * {pixels, wg_errors}
		arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));

		arg->histsize=histsize*sizeof(int);
		arg->hist=(int*)malloc(arg->histsize);

		arg->statssize=statssize;
		arg->stats=(unsigned*)malloc(statssize);
		if(!arg->pixels||!arg->hist||!arg->stats)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=arg->bufsize;
		memusage+=arg->histsize;
		memusage+=arg->statssize;
		
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
			arg->y1=BLOCKSIZE*(kt+kt2);
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
				int blocksize=(arg->y2-arg->y1)*image->iw*image->nch*image->depth/8;
				{
					const int *group=rct_combinations[arg->bestrct];
					g_rct_usizes[arg->bestrct]+=arg->bestsize;
					g_usizes[group[0]*PRED_COUNT+arg->predidx[0]]+=arg->usizes[0];
					g_usizes[group[1]*PRED_COUNT+arg->predidx[1]]+=arg->usizes[1];
					g_usizes[group[2]*PRED_COUNT+arg->predidx[2]]+=arg->usizes[2];
					g_csizes[group[0]*PRED_COUNT+arg->predidx[0]]+=arg->csizes[0];
					g_csizes[group[1]*PRED_COUNT+arg->predidx[1]]+=arg->csizes[1];
					g_csizes[group[2]*PRED_COUNT+arg->predidx[2]]+=arg->csizes[2];
				}
				if(loud)
				{
					if(!(kt+kt2))
						printf("block,  nrows,  usize,     best  ->  actual,  (actual-best)\n");
					printf(
						"[%3d]  %4d  %8d  %11.2lf -> %8zd bytes  (%+11.2lf)  %s %s %s %s\n",
						kt+kt2,
						arg->y2-arg->y1,
						blocksize,
						arg->bestsize,
						arg->list.nobj,
						(double)arg->list.nobj-arg->bestsize,
						rct_names[arg->bestrct],
						pred_names[arg->predidx[0]],
						pred_names[arg->predidx[1]],
						pred_names[arg->predidx[2]]
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
		free(arg->stats);
	}
	free(args);
	return 0;
}
void f26_curiosity()
{
#ifdef ENABLE_CURIOSITY
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
#endif
}