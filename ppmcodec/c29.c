#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define PROFILER

#ifdef _MSC_VER
	#define LOUD		//size & time
	#define PROFILE_SIZE
	#define ESTIMATE_SIZE

//	#define STREAM_STATSONLY
	#define ENABLE_GUIDE
	#define ANS_VAL
#define VAL_NCODERS 1
#endif

	#define WG_DISABLE_eW


//3*17+3=54 contexts
#define GRBITS 3
#define NCTX 17

#define NCODERS 32

#define PROBBITS 12	//12 bit max	James Bonfield's CDF2sym: {freq<<20 | bias<<8 | sym}

#define RANS_STATE_BITS 31
#define RANS_RENORM_BITS 16

#include"entropy.h"
#if 0
static ArrayHandle debugsimdstate=0;
static void anssimdval_push_states(__m256i *mstate)
{
	if(!debugsimdstate)
		ARRAY_ALLOC(int[NCODERS], debugsimdstate, 0, 0, 1024, 0);
	ARRAY_APPEND(debugsimdstate, mstate, 1, 1, 0);
}
static void anssimdval_check_states(__m256i *mstate)
{
	__m256i *mstate0=(__m256i*)array_back(&debugsimdstate);
	if(memcmp(mstate, mstate0, sizeof(int[NCODERS])))
	{
		unsigned original[NCODERS];
		unsigned wrong[NCODERS];
		memcpy(original, mstate0, sizeof(original));
		memcpy(wrong, mstate, sizeof(wrong));
		for(int k=0;k<32;++k)
			printf(" 0x%08X 0x%08X\n", original[k], wrong[k]);
		LOG_ERROR("");
	}
	--debugsimdstate->count;
}
static ArrayHandle debugstate=0;
static void ansval_push_state(unsigned state)
{
	if(!debugstate)
		ARRAY_ALLOC(int, debugstate, 0, 0, 1024, 0);
	ARRAY_APPEND(debugstate, &state, 1, 1, 0);
}
static void ansval_check_state(unsigned state)
{
	unsigned state0=*(unsigned*)array_back(&debugstate);
	if(state!=state0)
	{
		printf(" 0x%08X 0x%08X\n", state0, state);
		LOG_ERROR("");
	}
	--debugstate->count;
}
#endif
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static void guide_save(const unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(const unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		LOG_ERROR("");
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif
#ifdef PROFILE_SIZE
#include<stdarg.h>
static void profile_size(const unsigned char *dstbwdptr, const char *msg, ...)
{
	static ptrdiff_t size=0;
	static const unsigned char *prev=0;
	if(prev)
	{
		ptrdiff_t diff=prev-dstbwdptr;
		size+=diff;
		printf("%9td (%+9td) bytes", size, diff);
		if(msg)
		{
			va_list args;
			va_start(args, msg);
			printf("  ");
			vprintf(msg, args);
			va_end(args);
		}
		printf("\n");
	}
	prev=dstbwdptr;
}
#else
#define profile_size(...)
#endif
#ifdef PROFILER
static void prof_checkpoint(const char *msg)
{
	static double t=0;
	double t2=time_sec();
	if(t)
	{
		printf("%16.12lf sec", t2-t);
		if(msg)
			printf(" %s", msg);
		printf("\n");
	}
	t=t2;
}
#else
#define prof_checkpoint(...)
#endif

#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
	OCH(CX40) OCH(C0X4) OCH(C40X)\
	OCH(CX31) OCH(C3X1) OCH(C31X)\
	OCH(CX13) OCH(C1X3) OCH(C13X)\
	OCH(CX22) OCH(C2X2) OCH(C22X)
#if 0
#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
	OCH(Y310) OCH(Y031) OCH(Y103)\
	OCH(Y301) OCH(Y130) OCH(Y013)\
	OCH(Y211) OCH(Y121) OCH(Y112)\
	OCH(CX40) OCH(C0X4) OCH(C40X)\
	OCH(CX31) OCH(C3X1) OCH(C31X)\
	OCH(CX13) OCH(C1X3) OCH(C13X)\
	OCH(CX22) OCH(C2X2) OCH(C22X)
#endif
typedef enum _OCHIndex
{
#define OCH(X) OCH_##X,
	OCHLIST
#undef  OCH
	OCH_COUNT,
	OCH_C4X0=OCH_CX40,
	OCH_C04X=OCH_C0X4,
	OCH_CX04=OCH_C40X,
	OCH_R=OCH_Y400,
	OCH_G=OCH_Y040,
	OCH_B=OCH_Y004,
	OCH_BG=OCH_C04X,
	OCH_BR=OCH_C40X,
	OCH_RG=OCH_CX40,
	OCH_RB=OCH_CX04,
	OCH_GB=OCH_C0X4,
	OCH_GR=OCH_C4X0,
	OCH_R1=OCH_CX13,
	OCH_G1=OCH_C3X1,
	OCH_B1=OCH_C13X,
	OCH_R2=OCH_CX22,
	OCH_G2=OCH_C2X2,
	OCH_B2=OCH_C22X,
	OCH_R3=OCH_CX31,
	OCH_G3=OCH_C1X3,
	OCH_B3=OCH_C31X,
} OCHIndex;
static const char *och_names[]=
{
#define OCH(X) #X,
	OCHLIST
#undef  OCH
};
typedef enum _RCTInfoIdx
{
	II_OCH_Y,
	II_OCH_U,
	II_OCH_V,

	II_PERM_Y,
	II_PERM_U,
	II_PERM_V,

	II_COEFF_U_SUB_Y,
	II_COEFF_V_SUB_Y,
	II_COEFF_V_SUB_U,

//	II_COEFF_Y_SUB_U,
//	II_COEFF_Y_SUB_V,
//	II_COEFF_U_SUB_V_NBLI,
//	II_COEFF_V_SUB_U_NBLI,

	II_COUNT,
} RCTInfoIdx;
//YUV = RCT * RGB	watch out for permutation in last row
//luma: averaging	chroma: subtraction
//example: _400_40X_3X1 == [1 0 0; -1 0 1; -3/4 1 -1/4]
#define RCTLIST\
	RCT(_400_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_400_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4)\
	RCT(_400_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0)\
	RCT(_040_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0)\
	RCT(_040_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4)\
	RCT(_004_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4)\
	RCT(_004_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0)\
	RCT(_040_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0)\
	RCT(_040_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4)\
	RCT(_040_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4)\
	RCT(_004_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0)\
	RCT(_004_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4)\
	RCT(_004_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4)\
	RCT(_400_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0)\
	RCT(_400_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4)\
	RCT(_400_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)\
	RCT(_400_0X0_13X,	OCH_R,		OCH_G,		OCH_B1,		0, 1, 2,	0,  1, 3)\
	RCT(_400_4X0_13X,	OCH_R,		OCH_GR,		OCH_B1,		0, 1, 2,	4,  1, 3)\
	RCT(_400_00X_3X1,	OCH_R,		OCH_B,		OCH_G1,		0, 2, 1,	0,  3, 1)\
	RCT(_400_40X_3X1,	OCH_R,		OCH_BR,		OCH_G1,		0, 2, 1,	4,  3, 1)\
	RCT(_040_00X_X13,	OCH_G,		OCH_B,		OCH_R1,		1, 2, 0,	0,  1, 3)\
	RCT(_040_04X_X13,	OCH_G,		OCH_BG,		OCH_R1,		1, 2, 0,	4,  1, 3)\
	RCT(_040_X40_13X,	OCH_G,		OCH_RG,		OCH_B1,		1, 0, 2,	4,  3, 1)\
	RCT(_004_X04_3X1,	OCH_B,		OCH_RB,		OCH_G1,		2, 0, 1,	4,  1, 3)\
	RCT(_004_04X_X13,	OCH_B,		OCH_GB,		OCH_R1,		2, 1, 0,	4,  3, 1)\
	RCT(_400_0X0_22X,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,  2, 2)\
	RCT(_400_4X0_22X,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	4,  2, 2)\
	RCT(_400_00X_2X2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,  2, 2)\
	RCT(_400_40X_2X2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	4,  2, 2)\
	RCT(_040_00X_X22,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,  2, 2)\
	RCT(_040_04X_X22,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	4,  2, 2)\
	RCT(_040_X40_22X,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	4,  2, 2)\
	RCT(_004_X04_2X2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	4,  2, 2)\
	RCT(_004_0X4_X22,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	4,  2, 2)\
	RCT(_400_0X0_31X,	OCH_R,		OCH_G,		OCH_B3,		0, 1, 2,	0,  3, 1)\
	RCT(_400_4X0_31X,	OCH_R,		OCH_GR,		OCH_B3,		0, 1, 2,	4,  3, 1)\
	RCT(_400_00X_1X3,	OCH_R,		OCH_B,		OCH_G3,		0, 2, 1,	0,  1, 3)\
	RCT(_400_40X_1X3,	OCH_R,		OCH_BR,		OCH_G3,		0, 2, 1,	4,  1, 3)\
	RCT(_040_00X_X31,	OCH_G,		OCH_B,		OCH_R3,		1, 2, 0,	0,  3, 1)\
	RCT(_040_04X_X31,	OCH_G,		OCH_BG,		OCH_R3,		1, 2, 0,	4,  3, 1)\
	RCT(_040_X40_31X,	OCH_G,		OCH_RG,		OCH_B3,		1, 0, 2,	4,  1, 3)\
	RCT(_004_X04_1X3,	OCH_B,		OCH_RB,		OCH_G3,		2, 0, 1,	4,  3, 1)\
	RCT(_004_0X4_X31,	OCH_B,		OCH_GB,		OCH_R3,		2, 1, 0,	4,  1, 3)
#if 0
	RCT(_211_4X0_40X,	OCH_Y211,	OCH_C4X0,	OCH_C40X,	0, 1, 2,	4,  4, 0,	1, 1, 0, 0)\
	RCT(_211_4X0_31X,	OCH_Y211,	OCH_C4X0,	OCH_C31X,	0, 1, 2,	4,  4, 0,	1, 1, 0, 1)\
	RCT(_211_3X1_40X,	OCH_Y211,	OCH_C3X1,	OCH_C40X,	0, 1, 2,	4,  4, 0,	1, 1, 1, 0)\
	RCT(_310_4X0_40X,	OCH_Y310,	OCH_C4X0,	OCH_C40X,	0, 1, 2,	4,  4, 0,	1, 0, 0, 0)\
	RCT(_310_4X0_31X,	OCH_Y310,	OCH_C4X0,	OCH_C31X,	0, 1, 2,	4,  4, 0,	1, 0, 0, 1)\
	RCT(_310_3X1_40X,	OCH_Y310,	OCH_C3X1,	OCH_C40X,	0, 1, 2,	4,  4, 0,	1, 0, 1, 0)\
	RCT(_301_4X0_40X,	OCH_Y301,	OCH_C4X0,	OCH_C40X,	0, 1, 2,	4,  4, 0,	0, 1, 0, 0)\
	RCT(_301_4X0_31X,	OCH_Y301,	OCH_C4X0,	OCH_C31X,	0, 1, 2,	4,  4, 0,	0, 1, 0, 1)\
	RCT(_301_3X1_40X,	OCH_Y301,	OCH_C3X1,	OCH_C40X,	0, 1, 2,	4,  4, 0,	0, 1, 1, 0)\
	RCT(_121_04X_X40,	OCH_Y121,	OCH_C04X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	1, 1, 0, 0)\
	RCT(_121_04X_X31,	OCH_Y121,	OCH_C04X,	OCH_CX31,	1, 2, 0,	4,  4, 0,	1, 1, 0, 1)\
	RCT(_121_13X_X40,	OCH_Y121,	OCH_C13X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	1, 1, 1, 0)\
	RCT(_031_04X_X40,	OCH_Y031,	OCH_C04X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	0, 1, 0, 0)\
	RCT(_031_04X_X31,	OCH_Y031,	OCH_C04X,	OCH_CX31,	1, 2, 0,	4,  4, 0,	0, 1, 0, 1)\
	RCT(_031_13X_X40,	OCH_Y031,	OCH_C13X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	0, 1, 1, 0)\
	RCT(_130_40X_X40,	OCH_Y130,	OCH_C04X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	1, 0, 0, 0)\
	RCT(_130_40X_X31,	OCH_Y130,	OCH_C04X,	OCH_CX31,	1, 2, 0,	4,  4, 0,	1, 0, 0, 1)\
	RCT(_130_31X_X40,	OCH_Y130,	OCH_C13X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	1, 0, 1, 0)\
	RCT(_112_X04_0X4,	OCH_Y112,	OCH_CX04,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	1, 1, 0, 0)\
	RCT(_112_X04_1X3,	OCH_Y112,	OCH_CX04,	OCH_C1X3,	2, 0, 1,	4,  4, 0,	1, 1, 0, 1)\
	RCT(_112_X13_0X4,	OCH_Y112,	OCH_CX13,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	1, 1, 1, 0)\
	RCT(_013_X04_0X4,	OCH_Y013,	OCH_CX04,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	0, 1, 0, 0)\
	RCT(_013_X04_1X3,	OCH_Y013,	OCH_CX04,	OCH_C1X3,	2, 0, 1,	4,  4, 0,	0, 1, 0, 1)\
	RCT(_013_X13_0X4,	OCH_Y013,	OCH_CX13,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	0, 1, 1, 0)\
	RCT(_103_X40_0X4,	OCH_Y103,	OCH_CX04,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	1, 0, 0, 0)\
	RCT(_103_X40_1X3,	OCH_Y103,	OCH_CX04,	OCH_C1X3,	2, 0, 1,	4,  4, 0,	1, 0, 0, 1)\
	RCT(_103_X31_0X4,	OCH_Y103,	OCH_CX13,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	1, 0, 1, 0)
#endif
typedef enum _RCTIndex
{
#define RCT(LABEL, ...) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][II_COUNT]=
{
#define RCT(LABEL, ...) {__VA_ARGS__},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, ...) #LABEL,
	RCTLIST
#undef  RCT
};


#define WG_NPREDS	8	//multiple of 4
AWM_INLINE void wg_mix_pt2(const __m256i *wgpreds, const __m256i *prederrors, __m256 *result)//process halflane
{
	//mix(const short wgpreds[8, stride 6][16], short prederrors[8][16]) -> float[8][16]
	__m256 pr[16], coeff[16];

	//convert 16 shorts * 8 regs -> 8 floats * 16 regs		each 2 float regs make half a wave
	{
		__m128i hi[8];
		__m256i i32[8];
		hi[0]=_mm256_extracti128_si256(prederrors[0], 1);
		hi[1]=_mm256_extracti128_si256(prederrors[1], 1);
		hi[2]=_mm256_extracti128_si256(prederrors[2], 1);
		hi[3]=_mm256_extracti128_si256(prederrors[3], 1);
		hi[4]=_mm256_extracti128_si256(prederrors[4], 1);
		hi[5]=_mm256_extracti128_si256(prederrors[5], 1);
		hi[6]=_mm256_extracti128_si256(prederrors[6], 1);
		hi[7]=_mm256_extracti128_si256(prederrors[7], 1);
		i32[0]=_mm256_cvtepi16_epi32(hi[0]);
		i32[1]=_mm256_cvtepi16_epi32(hi[1]);
		i32[2]=_mm256_cvtepi16_epi32(hi[2]);
		i32[3]=_mm256_cvtepi16_epi32(hi[3]);
		i32[4]=_mm256_cvtepi16_epi32(hi[4]);
		i32[5]=_mm256_cvtepi16_epi32(hi[5]);
		i32[6]=_mm256_cvtepi16_epi32(hi[6]);
		i32[7]=_mm256_cvtepi16_epi32(hi[7]);
		coeff[0x1]=_mm256_cvtepi32_ps(i32[0]);
		coeff[0x3]=_mm256_cvtepi32_ps(i32[1]);
		coeff[0x5]=_mm256_cvtepi32_ps(i32[2]);
		coeff[0x7]=_mm256_cvtepi32_ps(i32[3]);
		coeff[0x9]=_mm256_cvtepi32_ps(i32[4]);
		coeff[0xB]=_mm256_cvtepi32_ps(i32[5]);
		coeff[0xD]=_mm256_cvtepi32_ps(i32[6]);
		coeff[0xF]=_mm256_cvtepi32_ps(i32[7]);
		i32[0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[0]));
		i32[1]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[1]));
		i32[2]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[2]));
		i32[3]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[3]));
		i32[4]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[4]));
		i32[5]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[5]));
		i32[6]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[6]));
		i32[7]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[7]));
		coeff[0x0]=_mm256_cvtepi32_ps(i32[0]);
		coeff[0x2]=_mm256_cvtepi32_ps(i32[1]);
		coeff[0x4]=_mm256_cvtepi32_ps(i32[2]);
		coeff[0x6]=_mm256_cvtepi32_ps(i32[3]);
		coeff[0x8]=_mm256_cvtepi32_ps(i32[4]);
		coeff[0xA]=_mm256_cvtepi32_ps(i32[5]);
		coeff[0xC]=_mm256_cvtepi32_ps(i32[6]);
		coeff[0xE]=_mm256_cvtepi32_ps(i32[7]);
	}
	
	//convert 16 shorts * 8 regs -> 8 floats * 16 regs
	{
		__m256i t0[16];
		t0[0*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[0*6+0], 1));
		t0[1*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[1*6+0], 1));
		t0[2*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[2*6+0], 1));
		t0[3*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[3*6+0], 1));
		t0[4*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[4*6+0], 1));
		t0[5*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[5*6+0], 1));
		t0[6*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[6*6+0], 1));
		t0[7*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[7*6+0], 1));
		t0[0*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[0*6+0]));
		t0[1*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[1*6+0]));
		t0[2*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[2*6+0]));
		t0[3*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[3*6+0]));
		t0[4*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[4*6+0]));
		t0[5*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[5*6+0]));
		t0[6*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[6*6+0]));
		t0[7*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[7*6+0]));

		pr[0x0]=_mm256_cvtepi32_ps(t0[0x0]);
		pr[0x1]=_mm256_cvtepi32_ps(t0[0x1]);
		pr[0x2]=_mm256_cvtepi32_ps(t0[0x2]);
		pr[0x3]=_mm256_cvtepi32_ps(t0[0x3]);
		pr[0x4]=_mm256_cvtepi32_ps(t0[0x4]);
		pr[0x5]=_mm256_cvtepi32_ps(t0[0x5]);
		pr[0x6]=_mm256_cvtepi32_ps(t0[0x6]);
		pr[0x7]=_mm256_cvtepi32_ps(t0[0x7]);
		pr[0x8]=_mm256_cvtepi32_ps(t0[0x8]);
		pr[0x9]=_mm256_cvtepi32_ps(t0[0x9]);
		pr[0xA]=_mm256_cvtepi32_ps(t0[0xA]);
		pr[0xB]=_mm256_cvtepi32_ps(t0[0xB]);
		pr[0xC]=_mm256_cvtepi32_ps(t0[0xC]);
		pr[0xD]=_mm256_cvtepi32_ps(t0[0xD]);
		pr[0xE]=_mm256_cvtepi32_ps(t0[0xE]);
		pr[0xF]=_mm256_cvtepi32_ps(t0[0xF]);
	}

	//coeff = 1/error
	coeff[0x0]=_mm256_rcp_ps(coeff[0x0]);//FIXME max error < (1.5 / 0x1000)
	coeff[0x1]=_mm256_rcp_ps(coeff[0x1]);
	coeff[0x2]=_mm256_rcp_ps(coeff[0x2]);
	coeff[0x3]=_mm256_rcp_ps(coeff[0x3]);
	coeff[0x4]=_mm256_rcp_ps(coeff[0x4]);
	coeff[0x5]=_mm256_rcp_ps(coeff[0x5]);
	coeff[0x6]=_mm256_rcp_ps(coeff[0x6]);
	coeff[0x7]=_mm256_rcp_ps(coeff[0x7]);
	coeff[0x8]=_mm256_rcp_ps(coeff[0x8]);
	coeff[0x9]=_mm256_rcp_ps(coeff[0x9]);
	coeff[0xA]=_mm256_rcp_ps(coeff[0xA]);
	coeff[0xB]=_mm256_rcp_ps(coeff[0xB]);
	coeff[0xC]=_mm256_rcp_ps(coeff[0xC]);
	coeff[0xD]=_mm256_rcp_ps(coeff[0xD]);
	coeff[0xE]=_mm256_rcp_ps(coeff[0xE]);
	coeff[0xF]=_mm256_rcp_ps(coeff[0xF]);

	//pr = (sum: coeff[i]*wg[i])/(sum: coeff[i])	NPREDS=8  (*16 interleaved floats)
	pr[0x0]=_mm256_mul_ps(pr[0x0], coeff[0x0]);
	pr[0x1]=_mm256_mul_ps(pr[0x1], coeff[0x1]);
	pr[0x2]=_mm256_mul_ps(pr[0x2], coeff[0x2]);
	pr[0x3]=_mm256_mul_ps(pr[0x3], coeff[0x3]);
	pr[0x4]=_mm256_mul_ps(pr[0x4], coeff[0x4]);
	pr[0x5]=_mm256_mul_ps(pr[0x5], coeff[0x5]);
	pr[0x6]=_mm256_mul_ps(pr[0x6], coeff[0x6]);
	pr[0x7]=_mm256_mul_ps(pr[0x7], coeff[0x7]);
	pr[0x8]=_mm256_mul_ps(pr[0x8], coeff[0x8]);
	pr[0x9]=_mm256_mul_ps(pr[0x9], coeff[0x9]);
	pr[0xA]=_mm256_mul_ps(pr[0xA], coeff[0xA]);
	pr[0xB]=_mm256_mul_ps(pr[0xB], coeff[0xB]);
	pr[0xC]=_mm256_mul_ps(pr[0xC], coeff[0xC]);
	pr[0xD]=_mm256_mul_ps(pr[0xD], coeff[0xD]);
	pr[0xE]=_mm256_mul_ps(pr[0xE], coeff[0xE]);
	pr[0xF]=_mm256_mul_ps(pr[0xF], coeff[0xF]);

	coeff[0x0]=_mm256_add_ps(coeff[0x0], coeff[0x2]);
	coeff[0x1]=_mm256_add_ps(coeff[0x1], coeff[0x3]);
	coeff[0x4]=_mm256_add_ps(coeff[0x4], coeff[0x6]);
	coeff[0x5]=_mm256_add_ps(coeff[0x5], coeff[0x7]);
	coeff[0x8]=_mm256_add_ps(coeff[0x8], coeff[0xA]);
	coeff[0x9]=_mm256_add_ps(coeff[0x9], coeff[0xB]);
	coeff[0xC]=_mm256_add_ps(coeff[0xC], coeff[0xE]);
	coeff[0xD]=_mm256_add_ps(coeff[0xD], coeff[0xF]);

	coeff[0x0]=_mm256_add_ps(coeff[0x0], coeff[0x4]);
	coeff[0x1]=_mm256_add_ps(coeff[0x1], coeff[0x5]);
	coeff[0x8]=_mm256_add_ps(coeff[0x8], coeff[0xC]);
	coeff[0x9]=_mm256_add_ps(coeff[0x9], coeff[0xD]);
	
	coeff[0x0]=_mm256_add_ps(coeff[0x0], coeff[0x8]);
	coeff[0x1]=_mm256_add_ps(coeff[0x1], coeff[0x9]);

	pr[0x0]=_mm256_add_ps(pr[0x0], pr[0x2]);
	pr[0x1]=_mm256_add_ps(pr[0x1], pr[0x3]);
	pr[0x4]=_mm256_add_ps(pr[0x4], pr[0x6]);
	pr[0x5]=_mm256_add_ps(pr[0x5], pr[0x7]);
	pr[0x8]=_mm256_add_ps(pr[0x8], pr[0xA]);
	pr[0x9]=_mm256_add_ps(pr[0x9], pr[0xB]);
	pr[0xC]=_mm256_add_ps(pr[0xC], pr[0xE]);
	pr[0xD]=_mm256_add_ps(pr[0xD], pr[0xF]);

	pr[0x0]=_mm256_add_ps(pr[0x0], pr[0x4]);
	pr[0x1]=_mm256_add_ps(pr[0x1], pr[0x5]);
	pr[0x8]=_mm256_add_ps(pr[0x8], pr[0xC]);
	pr[0x9]=_mm256_add_ps(pr[0x9], pr[0xD]);
	
	pr[0x0]=_mm256_add_ps(pr[0x0], pr[0x8]);
	pr[0x1]=_mm256_add_ps(pr[0x1], pr[0x9]);

	result[0]=_mm256_div_ps(pr[0], coeff[0]);
	result[1]=_mm256_div_ps(pr[4], coeff[4]);
}
AWM_INLINE void wg_mix(int kc, const short *wgWerrors, const short *rowscurrNNptr, const short *rowsNptr, const __m256i *wgpreds, __m256i *preds)//process channel wave (NCODERS lanes)
{
	//pred/errors layout:  {Y[NCODERS], U[NCODERS], V[NCODERS]} * NPREDS
	const short *chWerrors=wgWerrors+(ptrdiff_t)NCODERS*kc;
	const short *Nptr=rowsNptr+(ptrdiff_t)NCODERS*kc;
	const short *currNNptr=rowscurrNNptr+(ptrdiff_t)NCODERS*kc;
	__m256i prederrors[WG_NPREDS*NCODERS*sizeof(short)/sizeof(__m256i)];//2 regs per lane  *  NPREDS

	//errors = pI+pNW+2*pN+pNE+pNNE		(__m256i*)(erows[-Y]+((P+X*NPREDS)*3+C)*NCODERS)+R
	prederrors[0x0]=_mm256_load_si256((const __m256i*)(Nptr+((0+0*WG_NPREDS)*3+0)*NCODERS)+0);//pN
	prederrors[0x1]=_mm256_load_si256((const __m256i*)(Nptr+((0+0*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0x2]=_mm256_load_si256((const __m256i*)(Nptr+((1+0*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0x3]=_mm256_load_si256((const __m256i*)(Nptr+((1+0*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0x4]=_mm256_load_si256((const __m256i*)(Nptr+((2+0*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0x5]=_mm256_load_si256((const __m256i*)(Nptr+((2+0*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0x6]=_mm256_load_si256((const __m256i*)(Nptr+((3+0*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0x7]=_mm256_load_si256((const __m256i*)(Nptr+((3+0*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0x8]=_mm256_load_si256((const __m256i*)(Nptr+((4+0*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0x9]=_mm256_load_si256((const __m256i*)(Nptr+((4+0*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0xA]=_mm256_load_si256((const __m256i*)(Nptr+((5+0*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0xB]=_mm256_load_si256((const __m256i*)(Nptr+((5+0*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0xC]=_mm256_load_si256((const __m256i*)(Nptr+((6+0*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0xD]=_mm256_load_si256((const __m256i*)(Nptr+((6+0*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0xE]=_mm256_load_si256((const __m256i*)(Nptr+((7+0*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0xF]=_mm256_load_si256((const __m256i*)(Nptr+((7+0*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0x0]=_mm256_slli_epi16(prederrors[0x0], 1);
	prederrors[0x1]=_mm256_slli_epi16(prederrors[0x1], 1);
	prederrors[0x2]=_mm256_slli_epi16(prederrors[0x2], 1);
	prederrors[0x3]=_mm256_slli_epi16(prederrors[0x3], 1);
	prederrors[0x4]=_mm256_slli_epi16(prederrors[0x4], 1);
	prederrors[0x5]=_mm256_slli_epi16(prederrors[0x5], 1);
	prederrors[0x6]=_mm256_slli_epi16(prederrors[0x6], 1);
	prederrors[0x7]=_mm256_slli_epi16(prederrors[0x7], 1);
	prederrors[0x8]=_mm256_slli_epi16(prederrors[0x8], 1);
	prederrors[0x9]=_mm256_slli_epi16(prederrors[0x9], 1);
	prederrors[0xA]=_mm256_slli_epi16(prederrors[0xA], 1);
	prederrors[0xB]=_mm256_slli_epi16(prederrors[0xB], 1);
	prederrors[0xC]=_mm256_slli_epi16(prederrors[0xC], 1);
	prederrors[0xD]=_mm256_slli_epi16(prederrors[0xD], 1);
	prederrors[0xE]=_mm256_slli_epi16(prederrors[0xE], 1);
	prederrors[0xF]=_mm256_slli_epi16(prederrors[0xF], 1);
#ifdef WG_DISABLE_eW
	(void)chWerrors;
#else
	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], _mm256_load_si256((const __m256i*)(chWerrors+(0*3+0)*NCODERS)+0));//pI	(__m256i*)(wgWerrors+(P*3+C)*NCODERS)+R
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], _mm256_load_si256((const __m256i*)(chWerrors+(0*3+0)*NCODERS)+1));
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], _mm256_load_si256((const __m256i*)(chWerrors+(1*3+0)*NCODERS)+0));
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], _mm256_load_si256((const __m256i*)(chWerrors+(1*3+0)*NCODERS)+1));
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], _mm256_load_si256((const __m256i*)(chWerrors+(2*3+0)*NCODERS)+0));
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], _mm256_load_si256((const __m256i*)(chWerrors+(2*3+0)*NCODERS)+1));
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], _mm256_load_si256((const __m256i*)(chWerrors+(3*3+0)*NCODERS)+0));
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], _mm256_load_si256((const __m256i*)(chWerrors+(3*3+0)*NCODERS)+1));
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], _mm256_load_si256((const __m256i*)(chWerrors+(4*3+0)*NCODERS)+0));
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], _mm256_load_si256((const __m256i*)(chWerrors+(4*3+0)*NCODERS)+1));
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], _mm256_load_si256((const __m256i*)(chWerrors+(5*3+0)*NCODERS)+0));
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], _mm256_load_si256((const __m256i*)(chWerrors+(5*3+0)*NCODERS)+1));
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], _mm256_load_si256((const __m256i*)(chWerrors+(6*3+0)*NCODERS)+0));
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], _mm256_load_si256((const __m256i*)(chWerrors+(6*3+0)*NCODERS)+1));
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], _mm256_load_si256((const __m256i*)(chWerrors+(7*3+0)*NCODERS)+0));
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], _mm256_load_si256((const __m256i*)(chWerrors+(7*3+0)*NCODERS)+1));
#endif

	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], _mm256_load_si256((__m256i*)(Nptr+((0-1*WG_NPREDS)*3+0)*NCODERS)+0));//pNW
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], _mm256_load_si256((__m256i*)(Nptr+((0-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], _mm256_load_si256((__m256i*)(Nptr+((1-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], _mm256_load_si256((__m256i*)(Nptr+((1-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], _mm256_load_si256((__m256i*)(Nptr+((2-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], _mm256_load_si256((__m256i*)(Nptr+((2-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], _mm256_load_si256((__m256i*)(Nptr+((3-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], _mm256_load_si256((__m256i*)(Nptr+((3-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], _mm256_load_si256((__m256i*)(Nptr+((4-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], _mm256_load_si256((__m256i*)(Nptr+((4-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], _mm256_load_si256((__m256i*)(Nptr+((5-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], _mm256_load_si256((__m256i*)(Nptr+((5-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], _mm256_load_si256((__m256i*)(Nptr+((6-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], _mm256_load_si256((__m256i*)(Nptr+((6-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], _mm256_load_si256((__m256i*)(Nptr+((7-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], _mm256_load_si256((__m256i*)(Nptr+((7-1*WG_NPREDS)*3+0)*NCODERS)+1));

	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], _mm256_load_si256((__m256i*)(Nptr+((0+1*WG_NPREDS)*3+0)*NCODERS)+0));//pNE
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], _mm256_load_si256((__m256i*)(Nptr+((0+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], _mm256_load_si256((__m256i*)(Nptr+((1+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], _mm256_load_si256((__m256i*)(Nptr+((1+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], _mm256_load_si256((__m256i*)(Nptr+((2+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], _mm256_load_si256((__m256i*)(Nptr+((2+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], _mm256_load_si256((__m256i*)(Nptr+((3+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], _mm256_load_si256((__m256i*)(Nptr+((3+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], _mm256_load_si256((__m256i*)(Nptr+((4+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], _mm256_load_si256((__m256i*)(Nptr+((4+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], _mm256_load_si256((__m256i*)(Nptr+((5+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], _mm256_load_si256((__m256i*)(Nptr+((5+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], _mm256_load_si256((__m256i*)(Nptr+((6+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], _mm256_load_si256((__m256i*)(Nptr+((6+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], _mm256_load_si256((__m256i*)(Nptr+((7+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], _mm256_load_si256((__m256i*)(Nptr+((7+1*WG_NPREDS)*3+0)*NCODERS)+1));

	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], _mm256_load_si256((__m256i*)(currNNptr+((0+1*WG_NPREDS)*3+0)*NCODERS)+0));//pNNE
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], _mm256_load_si256((__m256i*)(currNNptr+((0+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], _mm256_load_si256((__m256i*)(currNNptr+((1+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], _mm256_load_si256((__m256i*)(currNNptr+((1+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], _mm256_load_si256((__m256i*)(currNNptr+((2+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], _mm256_load_si256((__m256i*)(currNNptr+((2+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], _mm256_load_si256((__m256i*)(currNNptr+((3+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], _mm256_load_si256((__m256i*)(currNNptr+((3+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], _mm256_load_si256((__m256i*)(currNNptr+((4+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], _mm256_load_si256((__m256i*)(currNNptr+((4+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], _mm256_load_si256((__m256i*)(currNNptr+((5+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], _mm256_load_si256((__m256i*)(currNNptr+((5+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], _mm256_load_si256((__m256i*)(currNNptr+((6+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], _mm256_load_si256((__m256i*)(currNNptr+((6+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], _mm256_load_si256((__m256i*)(currNNptr+((7+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], _mm256_load_si256((__m256i*)(currNNptr+((7+1*WG_NPREDS)*3+0)*NCODERS)+1));

	__m256i one=_mm256_set1_epi16(1);
	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], one);
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], one);
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], one);
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], one);
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], one);
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], one);
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], one);
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], one);
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], one);
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], one);
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], one);
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], one);
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], one);
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], one);
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], one);
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], one);

	__m256 result[4];
	wg_mix_pt2(wgpreds+(ptrdiff_t)2*kc+0, prederrors+0, result+0);//P*6+2*C+R
	wg_mix_pt2(wgpreds+(ptrdiff_t)2*kc+1, prederrors+8, result+2);

	__m256i ires[4];
	ires[0]=_mm256_cvtps_epi32(result[0]);
	ires[1]=_mm256_cvtps_epi32(result[1]);
	ires[2]=_mm256_cvtps_epi32(result[2]);
	ires[3]=_mm256_cvtps_epi32(result[3]);
	ires[0]=_mm256_packs_epi32(ires[0], ires[1]);
	ires[2]=_mm256_packs_epi32(ires[2], ires[3]);
	preds[0]=_mm256_permute4x64_epi64(ires[0], _MM_SHUFFLE(3, 1, 2, 0));
	preds[1]=_mm256_permute4x64_epi64(ires[2], _MM_SHUFFLE(3, 1, 2, 0));
}


//https://github.com/rygorous/ryg_rans
//https://github.com/samtools/htscodecs
typedef struct _rANS_SIMD_SymInfo	//16 bytes/level	4KB/ctx = 1<<12 bytes
{
	unsigned smax, invf, cdf;
	unsigned short negf, sh;
} rANS_SIMD_SymInfo;
static void enc_hist2stats(int *hist, rANS_SIMD_SymInfo *syminfo, long long *ctxmask, int ctxidx)
{
	int sum=0, count=0;
	for(int ks=0;ks<256;++ks)
	{
		int freq=hist[ks];
		sum+=freq;
		count+=freq!=0;
	}
	int rare=sum<12*256/8;
	*ctxmask|=(unsigned long long)rare<<ctxidx;
	if(rare)
	{
		for(int ks=0;ks<256;++ks)//bypass for now	FIXME try emitting mask with present symbols, and handle degenerate distribution
			hist[ks]=1;
		sum=256;
		count=256;
	}
	int sum2=0;
	for(int ks=0, ks2=0;ks<256;++ks)//absent symbols get zero freqs
	{
		int freq=hist[ks];
		hist[ks]=(int)(sum2*((1ULL<<PROBBITS)-count)/sum)+ks2;
		ks2+=freq!=0;
		sum2+=freq;
	}
	//for(int ks=0;ks<256;++ks)//never allows zero freqs	INEFFICIENT
	//{
	//	int freq=hist[ks];
	//	hist[ks]=(int)(sum2*((1ULL<<PROBBITS)-256)/sum)+ks;
	//	sum2+=freq;
	//}
#ifdef LOUD
	double e=0, norm=1./0x1000;
	for(int ks=0;ks<256;++ks)//estimate
	{
		int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
		if(freq)
		{
			double p=freq*norm;
			e-=p*log2(p);
		}
		if(e!=e)
			LOG_ERROR("");
	}
	e*=sum/8.;
	printf("%12.2lf bytes\n", e);
#endif
	int next=1<<PROBBITS;
	for(int ks=255;ks>=0;--ks)
	{
		rANS_SIMD_SymInfo *info=syminfo+ks;
		int curr=hist[ks];
		int freq=next-curr;
		next=curr;
		hist[ks]=freq;
		info->smax=(freq<<(RANS_STATE_BITS-PROBBITS))-1;//rescale freq to match the rANS state, and decrement to use _mm256_cmpgt_epi32 instead of '>='
		info->cdf=curr;
		info->negf=(1<<PROBBITS)-freq;
		//encoding:  state  =  q<<16|(cdf+r)
		//div-free:  state  =  q*M+cdf+state-q*freq  =  state+q*(M-freq)+cdf  =  state+(state*invf>>sh)*(M-freq)+cdf
		//sh = FLOOR_LOG2(freq)+32
		//invf = ceil(2^sh/freq)		state is 31 bits
		if(freq<2)
		{
			//ideally  q = state*inv(1)>>sh(1) = state*2^32>>32
			//here  q' = state*(2^32-1)>>32 = floor(state-state/2^32) = state-1  if  1 <= x < 2^32
			//enc  state = (state/1)*M+cdf+state%1  =  state+q*(M-1)+cdf
			//but  q' = state-1
			//so  state = state+(state-1+1)*(M-1)+cdf  =  state+q'*(M-1)+(cdf+M-1)
			info->sh=0;
			info->invf=0xFFFFFFFF;
			info->cdf+=(1<<PROBBITS)-1;
		}
		else
		{
			info->sh=FLOOR_LOG2(freq);
			info->invf=(int)(((0x100000000<<info->sh)+freq-1)/freq);
		}
	}
}
static void enc_packhist(BitPackerLIFO *ec, const int *hist, unsigned long long ctxmask, int ctxidx)//histogram must be normalized to PROBBITS, with spike at 128
{
	if(ctxmask>>ctxidx&1)
		return;
	int sum=0;
	unsigned short CDF[257];
	for(int ks=0;ks<256;++ks)//integrage to zigzag CDF to be packed backwards
	{
		int sym=((ks>>1^-(ks&1))+128)&255;
		int freq=hist[sym];
		CDF[ks]=sum;//separate buffer for faster access in 2nd loop
		sum+=freq;
	}
	CDF[256]=1<<PROBBITS;
	
//#ifdef _DEBUG
//	size_t checkpoint=(size_t)ec->dstbwdptr;
//	int count=0;
//	for(int ks=0;ks<256;++ks)//
//		count+=CDF[ks]<CDF[ks+1];
//#endif
	int cdfW=CDF[0];
	int CDFlevels=1<<PROBBITS;
	int startsym=0;
	for(int ks=1;ks<=256;++ks)//push nbypass
	{
		int next=CDF[ks], freq=next-cdfW, nlevels=(1<<PROBBITS)-cdfW;
		int nbypass=FLOOR_LOG2(CDFlevels);
		if(ks>1)
			nbypass-=7;
		if(nbypass<0)
			nbypass=0;
		CDF[ks]=nbypass<<PROBBITS|freq;
		cdfW=next;
		CDFlevels-=freq;
		startsym=ks;
		if(!CDFlevels)
			break;
	}
	for(int ks=startsym;ks>0;--ks)//encode GR
	{
		int freq=CDF[ks], nbypass=freq>>PROBBITS;
		freq&=(1<<PROBBITS)-1;
		int nzeros=freq>>nbypass, bypass=freq&((1<<nbypass)-1);
		if(nbypass)
			bitpacker_enc(ec, nbypass, bypass);
		bitpacker_enc(ec, 1, 1);
		while(nzeros)
		{
			bitpacker_enc(ec, 1, 0);
			--nzeros;
		}
#ifdef ANS_VAL
		ansval_push(&ks, sizeof(ks));
		//if(ctxidx==5||ctxidx==4)
		//	printf("PUSH [%3d] freq %8d  nzeros %8d  bypass %08X/%4d  size %8.2lf  state %016llX\n", ks, freq, nzeros, bypass, nbypass, (nzeros+1+nbypass)/8., ec->state);//
#endif
	}
//#ifdef _DEBUG
//	printf("ctx %3d  count %4d  size %8td  state %016llX\n", ctxidx, count, checkpoint-(size_t)ec->dstbwdptr, ec->state);//
//#endif
}
static void dec_unpackhist(BitPackerLIFO *ec, unsigned *CDF2sym, unsigned long long ctxmask, int ctxidx)
{
	unsigned short hist[257];
	if(ctxmask>>ctxidx&1)//rare context
	{
		for(int ks=0;ks<256;++ks)//bypass for now
			hist[ks]=(1<<PROBBITS)/256;
	}
	else
	{
		unsigned short CDF[257]={0};
		int CDFlevels=1<<PROBBITS;
		int cdfW=0;
		CDF[0]=0;
		for(int ks=0;ks<256;++ks)//decode GR
		{
			int freq=-1;//stop bit doesn't count
			int nbypass=FLOOR_LOG2(CDFlevels);
			int ks2=ks+1;
			if(ks2>1)
				nbypass-=7;
			if(nbypass<0)
				nbypass=0;
#ifdef ANS_VAL
			ansval_check(&ks2, sizeof(ks2));
#endif
			int bit=0;
			do
			{
				bit=bitpacker_dec(ec, 1);
				++freq;
			}while(!bit);
			if(nbypass)
				freq=freq<<nbypass|bitpacker_dec(ec, nbypass);

			cdfW+=freq;
			CDF[ks]=freq;
			CDFlevels-=freq;
			if(CDFlevels<=0)
			{
#ifdef _DEBUG
				if(CDFlevels<0)
					LOG_ERROR("CDF unpack error");
#endif
				break;
			}
		}
		if(CDFlevels)
			LOG_ERROR("CDF unpack error");
		for(int ks=0;ks<256;++ks)//undo zigzag
		{
			int sym=((ks>>1^-(ks&1))+128)&255;
			hist[sym]=CDF[ks];
		}
	}
	int sum=0;
	for(int ks=0;ks<256;++ks)//integrate
	{
		int freq=hist[ks];
		hist[ks]=sum;
		sum+=freq;
	}
	hist[256]=1<<PROBBITS;
	for(int ks=0;ks<256;++ks)//CDF2sym contains {freq, (state&0xFFF)-cdf, sym}
	{
		int cdf=hist[ks], next=hist[ks+1], freq=next-cdf;
		int val=(freq<<PROBBITS|0)<<8|ks;
		for(int ks2=cdf, bias=0;ks2<next;++ks2, ++bias)
			CDF2sym[ks2]=val|bias<<8;
	}
}
AWM_INLINE void dec_yuv(__m256i *mstate, __m256i *ctx0, __m256i *ctx1, const int *CDF2syms, const int *ans_permute, unsigned char **pstreamptr, const unsigned char *streamend, __m256i *syms)
{
	const unsigned char *streamptr=*pstreamptr;
	__m256i decctx[4];
	{
		__m128i t0=_mm256_extracti128_si256(*ctx0, 1);
		__m128i t1=_mm256_extracti128_si256(*ctx1, 1);
		decctx[0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(*ctx0));
		decctx[2]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(*ctx1));
		decctx[1]=_mm256_cvtepi16_epi32(t0);
		decctx[3]=_mm256_cvtepi16_epi32(t1);
	}
	decctx[0]=_mm256_slli_epi32(decctx[0], PROBBITS);
	decctx[1]=_mm256_slli_epi32(decctx[1], PROBBITS);
	decctx[2]=_mm256_slli_epi32(decctx[2], PROBBITS);
	decctx[3]=_mm256_slli_epi32(decctx[3], PROBBITS);
#ifdef ANS_VAL
	ansval_check(mstate, sizeof(int[VAL_NCODERS]));
#endif
	{
		__m256i mprobmask=_mm256_set1_epi32((1<<PROBBITS)-1);
		__m256i rem0=_mm256_and_si256(mstate[0], mprobmask);
		__m256i rem1=_mm256_and_si256(mstate[1], mprobmask);
		__m256i rem2=_mm256_and_si256(mstate[2], mprobmask);
		__m256i rem3=_mm256_and_si256(mstate[3], mprobmask);
		decctx[0]=_mm256_add_epi32(decctx[0], rem0);
		decctx[1]=_mm256_add_epi32(decctx[1], rem1);
		decctx[2]=_mm256_add_epi32(decctx[2], rem2);
		decctx[3]=_mm256_add_epi32(decctx[3], rem3);
	}
	decctx[0]=_mm256_i32gather_epi32(CDF2syms, decctx[0], sizeof(int));//FIXME try not using gather
	decctx[1]=_mm256_i32gather_epi32(CDF2syms, decctx[1], sizeof(int));
	decctx[2]=_mm256_i32gather_epi32(CDF2syms, decctx[2], sizeof(int));
	decctx[3]=_mm256_i32gather_epi32(CDF2syms, decctx[3], sizeof(int));

	//update		state = (state>>12)*freq+(rem-cdf)	rem-cdf is prebaked
	mstate[0]=_mm256_srli_epi32(mstate[0], PROBBITS);
	mstate[1]=_mm256_srli_epi32(mstate[1], PROBBITS);
	mstate[2]=_mm256_srli_epi32(mstate[2], PROBBITS);
	mstate[3]=_mm256_srli_epi32(mstate[3], PROBBITS);
	{
		__m256i mfreq0=_mm256_srli_epi32(decctx[0], PROBBITS+8);//1 <= freq <= 0xF01
		__m256i mfreq1=_mm256_srli_epi32(decctx[1], PROBBITS+8);
		__m256i mfreq2=_mm256_srli_epi32(decctx[2], PROBBITS+8);
		__m256i mfreq3=_mm256_srli_epi32(decctx[3], PROBBITS+8);
		mstate[0]=_mm256_mullo_epi32(mstate[0], mfreq0);//10 cycles
		mstate[1]=_mm256_mullo_epi32(mstate[1], mfreq1);
		mstate[2]=_mm256_mullo_epi32(mstate[2], mfreq2);
		mstate[3]=_mm256_mullo_epi32(mstate[3], mfreq3);
	}
	__m256i lobyte=_mm256_set1_epi32(255);
	__m256i msym0=_mm256_and_si256(decctx[0], lobyte);
	__m256i msym1=_mm256_and_si256(decctx[1], lobyte);
	__m256i msym2=_mm256_and_si256(decctx[2], lobyte);
	__m256i msym3=_mm256_and_si256(decctx[3], lobyte);
	msym0=_mm256_packus_epi32(msym0, msym1);
	msym1=_mm256_packus_epi32(msym2, msym3);
	msym0=_mm256_permute4x64_epi64(msym0, _MM_SHUFFLE(3, 1, 2, 0));
	msym1=_mm256_permute4x64_epi64(msym1, _MM_SHUFFLE(3, 1, 2, 0));
	{
		__m256i mbias0=_mm256_slli_epi32(decctx[0], PROBBITS);
		__m256i mbias1=_mm256_slli_epi32(decctx[1], PROBBITS);
		__m256i mbias2=_mm256_slli_epi32(decctx[2], PROBBITS);
		__m256i mbias3=_mm256_slli_epi32(decctx[3], PROBBITS);
		mbias0=_mm256_srli_epi32(mbias0, 32-PROBBITS);
		mbias1=_mm256_srli_epi32(mbias1, 32-PROBBITS);
		mbias2=_mm256_srli_epi32(mbias2, 32-PROBBITS);
		mbias3=_mm256_srli_epi32(mbias3, 32-PROBBITS);
		mstate[0]=_mm256_add_epi32(mstate[0], mbias0);
		mstate[1]=_mm256_add_epi32(mstate[1], mbias1);
		mstate[2]=_mm256_add_epi32(mstate[2], mbias2);
		mstate[3]=_mm256_add_epi32(mstate[3], mbias3);
	}
	//renorm
	{
#ifdef _DEBUG
		if(streamptr>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__m256i lo0=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr));
		__m256i smin=_mm256_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
		__m256i cond0=_mm256_cmpgt_epi32(smin, mstate[0]);//FIXME this is signed comparison
		__m256i cond1=_mm256_cmpgt_epi32(smin, mstate[1]);
		int mask0=_mm256_movemask_ps(_mm256_castsi256_ps(cond0));
		int mask1=_mm256_movemask_ps(_mm256_castsi256_ps(cond1));
		__m256i idx0=_mm256_load_si256((const __m256i*)ans_permute+mask0);
		__m256i idx1=_mm256_load_si256((const __m256i*)ans_permute+mask1);
		mask0=_mm_popcnt_u32(mask0);
		mask1=_mm_popcnt_u32(mask1);
		__m256i renorm0=_mm256_slli_epi32(mstate[0], 16);
		__m256i renorm1=_mm256_slli_epi32(mstate[1], 16);
		streamptr+=mask0*sizeof(short);
		lo0=_mm256_permutevar8x32_epi32(lo0, idx0);
#ifdef _DEBUG
		if(streamptr>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__m256i lo1=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr));
		streamptr+=mask1*sizeof(short);
		lo1=_mm256_permutevar8x32_epi32(lo1, idx1);
		renorm0=_mm256_or_si256(renorm0, lo0);
		renorm1=_mm256_or_si256(renorm1, lo1);
		
#ifdef _DEBUG
		if(streamptr>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__m256i lo2=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr));
		__m256i cond2=_mm256_cmpgt_epi32(smin, mstate[2]);
		__m256i cond3=_mm256_cmpgt_epi32(smin, mstate[3]);
		int mask2=_mm256_movemask_ps(_mm256_castsi256_ps(cond2));
		int mask3=_mm256_movemask_ps(_mm256_castsi256_ps(cond3));
		__m256i idx2=_mm256_load_si256((const __m256i*)ans_permute+mask2);
		__m256i idx3=_mm256_load_si256((const __m256i*)ans_permute+mask3);
		mask2=_mm_popcnt_u32(mask2);
		mask3=_mm_popcnt_u32(mask3);
		__m256i renorm2=_mm256_slli_epi32(mstate[2], 16);
		__m256i renorm3=_mm256_slli_epi32(mstate[3], 16);
		streamptr+=mask2*sizeof(short);
		lo2=_mm256_permutevar8x32_epi32(lo2, idx2);
#ifdef _DEBUG
		if(streamptr>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__m256i lo3=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr));
		streamptr+=mask3*sizeof(short);
		lo3=_mm256_permutevar8x32_epi32(lo3, idx3);
		renorm2=_mm256_or_si256(renorm2, lo2);
		renorm3=_mm256_or_si256(renorm3, lo3);

		mstate[0]=_mm256_blendv_epi8(mstate[0], renorm0, cond0);
		mstate[1]=_mm256_blendv_epi8(mstate[1], renorm1, cond1);
		mstate[2]=_mm256_blendv_epi8(mstate[2], renorm2, cond2);
		mstate[3]=_mm256_blendv_epi8(mstate[3], renorm3, cond3);
	}
	*pstreamptr=(unsigned char*)(size_t)streamptr;
}
static void interleave_fwd(unsigned char *image, int iw, int ih, unsigned char *temprow)
{
	int blockw=iw/NCODERS, blockstride=3*blockw;
	int rowstride=3*iw;
	unsigned char *tptrs[NCODERS]={0}, *tptrs0[NCODERS]={0};
	{
		unsigned char *temp=temprow;
		for(int k=0;k<NCODERS;++k)
		{
			tptrs[k]=temp;
			temp+=blockstride;
		}
	}
	memcpy(tptrs0, tptrs, sizeof(tptrs0));
	unsigned char *imptr=image;
	for(int ky=0;ky<ih;++ky)//interleave coder lanes
	{
		unsigned char *dstptr=imptr;
		memcpy(temprow, imptr, (ptrdiff_t)blockstride*NCODERS);
		memcpy(tptrs, tptrs0, sizeof(tptrs));
		for(int kx=0;kx<blockw;++kx)
		{
			//toy example with 2 coders:
			// [0] [2] [4] ...                                 [1] [3] [5] ...
			//|r00 g00 b00|r01 g01 b01|r02 g02 b02|r03 g03 b03|r04 g04 b04|r05 g05 b05|r06 g06 b06|r07 g07 b07|	src
			//|r00 r04|g00 g04|b00 b04|r01 r05|g01 g05|b01 b05|r02 r06|g02 g06|b02 b06|r03 r07|g03 g07|b03 b07|	dst
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
			for(int k=0;k<NCODERS;++k)
				*dstptr++=*tptrs[k]++;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
			for(int k=0;k<NCODERS;++k)
				*dstptr++=*tptrs[k]++;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
			for(int k=0;k<NCODERS;++k)
				*dstptr++=*tptrs[k]++;
		}
		imptr+=rowstride;
	}
}
static void interleave_inv(unsigned char *image, int iw, int ih, unsigned char *temprow)
{
	int blockw=iw/NCODERS, blockstride=3*blockw;
	int rowstride=3*iw;
	unsigned char *tptrs[NCODERS]={0}, *tptrs0[NCODERS]={0};
	{
		unsigned char *temp=temprow;
		for(int k=0;k<NCODERS;++k)
		{
			tptrs[k]=temp;
			temp+=blockstride;
		}
	}
	memcpy(tptrs0, tptrs, sizeof(tptrs0));
	unsigned char *imptr=image;
	for(int ky=0;ky<ih;++ky)//deinterleave coder lanes
	{
		unsigned char *dstptr=imptr;
		memcpy(tptrs, tptrs0, sizeof(tptrs));
		for(int kx=0;kx<blockw;++kx)
		{
			//toy example with 2 coders:
			// [0] [2] [4] ...                                 [1] [3] [5] ...
			//|r00 g00 b00|r01 g01 b01|r02 g02 b02|r03 g03 b03|r04 g04 b04|r05 g05 b05|r06 g06 b06|r07 g07 b07|	dst
			//|r00 r04|g00 g04|b00 b04|r01 r05|g01 g05|b01 b05|r02 r06|g02 g06|b02 b06|r03 r07|g03 g07|b03 b07|	src
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
			for(int k=0;k<NCODERS;++k)
				*tptrs[k]++=*dstptr++;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
			for(int k=0;k<NCODERS;++k)
				*tptrs[k]++=*dstptr++;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
			for(int k=0;k<NCODERS;++k)
				*tptrs[k]++=*dstptr++;
		}
		memcpy(imptr, temprow, (ptrdiff_t)blockstride*NCODERS);
		imptr+=rowstride;
	}
}
static void save_ppm(const char *fn, const unsigned char *image, int iw, int ih)
{
	FILE *fdst=fopen(fn, "wb");
	if(!fdst)
	{
		LOG_ERROR("Cannot open \"%s\" for writing", fdst);
		return;
	}
	fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	fwrite(image, 1, (ptrdiff_t)3*iw*ih, fdst);
	fclose(fdst);
}
int c29_codec(const char *srcfn, const char *dstfn, int nthreads0)
{
#ifdef ESTIMATE_SIZE
	double esize[3*NCODERS]={0};
#endif
#ifdef LOUD
	ptrdiff_t usize2=get_filesize(srcfn);
	double t=time_sec();
#endif
	prof_checkpoint(0);
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	int fwd=0, iw=0, ih=0, rowstride=0, blockw=0, remwidth=0, blockstride=0, mainwidth=0, mainbytes=0;
	ptrdiff_t usize=0, cap=0;
	unsigned char *image=0, *imptr=0, *streamptr=0, *streamstart=0, *streamend=0;
	unsigned char *context=0;
	int psize=0;
	short *pixels=0;
	int wgsize=0;
	short *wgerrors=0;
	//ptrdiff_t uheadersize=0;
	ptrdiff_t cheadersize=0, csize=0;
	{
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		int tag=0;
		fread(&tag, 1, 2, fsrc);
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('2'|'9'<<8))
		{
			LOG_ERROR("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
		{
			int temp=fgetc(fsrc);
			if(temp!='\n')
			{
				LOG_ERROR("Invalid PPM file");
				return 1;
			}
			int nread=fscanf(fsrc, "%d %d", &iw, &ih);
			if(nread!=2)
			{
				LOG_ERROR("Unsupported PPM file");
				return 1;
			}
			int vmax=0;
			nread=fscanf(fsrc, "%d", &vmax);
			if(nread!=1||vmax!=255)
			{
				LOG_ERROR("Unsupported PPM file");
				return 1;
			}
			temp=fgetc(fsrc);
			if(temp!='\n')
			{
				LOG_ERROR("Invalid PPM file");
				return 1;
			}
			//uheadersize=ftell(fsrc);
		}
		else
		{
			fread(&iw, 1, 4, fsrc);
			fread(&ih, 1, 4, fsrc);
			cheadersize=ftell(fsrc);
		}
		if(iw<1||ih<1)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		blockw=iw/NCODERS;
		remwidth=iw%NCODERS;
		blockstride=3*blockw;
		mainwidth=NCODERS*blockw;
		mainbytes=3*mainwidth;
		rowstride=3*iw;
		usize=(ptrdiff_t)3*iw*ih;
		cap=(ptrdiff_t)4*iw*ih;
		image=(unsigned char*)malloc(cap+sizeof(__m256i));
		if(!image)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		if(fwd)
		{
			fread(image, 1, usize, fsrc);//read image
			streamptr=streamstart=image+cap;//bwd-bwd ANS encoding
			profile_size(streamptr, "start");
		}
		else
		{
			csize=get_filesize(srcfn);
			streamptr=streamstart=image+cap-(csize-cheadersize)-sizeof(__m256i);
			streamend=image+cap-sizeof(__m256i);
			fread(streamstart, 1, csize-cheadersize, fsrc);//read stream
		}
		fclose(fsrc);
	}
	prof_checkpoint("fread");
	int bestrct=0, use_wg4=0;
	unsigned long long ctxmask=0;//3*NCTX+3 = 54 flags	0: rare context (bypass)  1: emit stats
	const int hsize=(int)sizeof(int[3*NCTX<<8]);//3 channels
	int *hists=fwd?(int*)malloc(hsize):0;//fwd-only
	const int rhsize=(int)sizeof(int[3*256]);
	int *rhist=fwd?(int*)malloc(rhsize):0;

	int CDF2syms_size=(int)sizeof(int[3*NCTX<<PROBBITS]);
	if(fwd)//DIV-free rANS encoder reuses these as SIMD symbol info
		CDF2syms_size=(int)sizeof(rANS_SIMD_SymInfo[3*NCTX<<8]);
	unsigned *CDF2syms=(unsigned*)_mm_malloc(CDF2syms_size, sizeof(__m128i*));

	int rCDF2syms_size=(int)sizeof(int[3<<PROBBITS]);
	if(fwd)
		rCDF2syms_size=(int)sizeof(rANS_SIMD_SymInfo[3<<8]);
	unsigned *rCDF2syms=(unsigned*)_mm_malloc(rCDF2syms_size, sizeof(__m128i*));

	int ans_permute_size=sizeof(__m256i[256]);
	int *ans_permute=(int*)_mm_malloc(ans_permute_size, sizeof(__m256i));

	psize=(int)sizeof(short[4*6*NCODERS])*(blockw+16);//4 padded rows  *  {Y*NCODERS, U*NCODERS, V*NCODERS,  eY*NCODERS, eU*NCODERS, eV*NCODERS} = 2*3*32 = 192 channels  ~48*iw bytes
	pixels=(short*)_mm_malloc(psize, sizeof(__m256i));//~188 KB for 4K/12MP
	wgsize=(int)sizeof(short[2*NCODERS*3*WG_NPREDS])*(blockw+16);//2 padded rows  *  {WGY*NCODERS, WGU*NCODERS, WGV*NCODERS} * NPREDS = 3*32*8 = 768 channels  ~96*iw bytes
	wgerrors=(short*)_mm_malloc(wgsize, sizeof(__m256i));//~375 KB for 4K/12MP		NNEerrors = currerrors
	const int wgstatesize=(int)sizeof(short[2*NCODERS*3*WG_NPREDS]);//{preds, Wprederrors} * NPREDS * 3 channels * NCODERS
	short *wgstate=(short*)_mm_malloc(wgstatesize, sizeof(__m256i));
	if((fwd&&(!hists||!rhist))||!CDF2syms||!rCDF2syms||!ans_permute||!pixels||!wgerrors||!wgstate)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(ans_permute, 0, ans_permute_size);
	if(fwd)
	{
		memset(hists, 0, hsize);
		memset(rhist, 0, rhsize);
#ifdef TEST_INTERLEAVE
		guide_save(image, iw, ih);
		save_ppm("20250227_1225AM_original.PPM", image, iw, ih);
		interleave_fwd(image, iw, ih, (unsigned char*)pixels);
		save_ppm("20250226_1153PM_interleaved.PPM", image, iw, ih);
		interleave_inv(image, iw, ih, (unsigned char*)pixels);
		save_ppm("20250227_1244AM_deinterleaved1.PPM", image, iw, ih);
		//interleave_fwd(image, iw, ih, (unsigned char*)pixels);
		//interleave_inv(image, iw, ih, (unsigned char*)pixels);
		{
			if(memcmp(image, g_image, usize))
			{
				const char *fn="20250226_0230AM_corrupt.PPM";
				save_ppm(fn, image, iw, ih);
				printf("corrupt  %s\n", fn);
				LOG_ERROR("");
			}
			printf("SUCCESS\n");
			exit(0);
		}
#endif
		interleave_fwd(image, iw, ih, (unsigned char*)pixels);
		prof_checkpoint("interleave");
		guide_save(image, iw, ih);
		{
			ALIGN(32) long long counters[OCH_COUNT]={0};
			__m256i mcounters[OCH_COUNT];//64-bit
			__m128i half8=_mm_set1_epi8(-128);
			__m256i wordmask=_mm256_set1_epi64x(0xFFFF);
			memset(mcounters, 0, sizeof(mcounters));
			imptr=image;
			for(int ky=0;ky<ih;++ky)//analysis
			{
				__m256i prev[OCH_COUNT*2];//16-bit
				memset(prev, 0, sizeof(prev));
				for(int kx=0;kx<rowstride-3*NCODERS;kx+=3*NCODERS)
				{
					__m256i r0=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx)+0), half8));
					__m256i r1=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx)+1), half8));
					__m256i g0=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx)+2), half8));
					__m256i g1=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx)+3), half8));
					__m256i b0=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx)+4), half8));
					__m256i b1=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx)+5), half8));
					r0=_mm256_slli_epi16(r0, 2);
					r1=_mm256_slli_epi16(r1, 2);
					g0=_mm256_slli_epi16(g0, 2);
					g1=_mm256_slli_epi16(g1, 2);
					b0=_mm256_slli_epi16(b0, 2);
					b1=_mm256_slli_epi16(b1, 2);
					__m256i rg0=_mm256_sub_epi16(r0, g0);
					__m256i rg1=_mm256_sub_epi16(r1, g1);
					__m256i gb0=_mm256_sub_epi16(g0, b0);
					__m256i gb1=_mm256_sub_epi16(g1, b1);
					__m256i br0=_mm256_sub_epi16(b0, r0);
					__m256i br1=_mm256_sub_epi16(b1, r1);
					__m256i t0, t1, t2, t3, t4, t5;
#define UPDATE(IDXA, IDXB, IDXC, A0, A1, B0, B1, C0, C1)\
	do\
	{\
		__m256i a0=_mm256_sub_epi16(A0, prev[IDXA*2+0]);\
		__m256i a1=_mm256_sub_epi16(A1, prev[IDXA*2+1]);\
		__m256i b0=_mm256_sub_epi16(B0, prev[IDXB*2+0]);\
		__m256i b1=_mm256_sub_epi16(B1, prev[IDXB*2+1]);\
		__m256i c0=_mm256_sub_epi16(C0, prev[IDXC*2+0]);\
		__m256i c1=_mm256_sub_epi16(C1, prev[IDXC*2+1]);\
		prev[IDXA*2+0]=A0;\
		prev[IDXA*2+1]=A1;\
		prev[IDXB*2+0]=B0;\
		prev[IDXB*2+1]=B1;\
		prev[IDXC*2+0]=C0;\
		prev[IDXC*2+1]=C1;\
		a0=_mm256_abs_epi16(a0);\
		a1=_mm256_abs_epi16(a1);\
		b0=_mm256_abs_epi16(b0);\
		b1=_mm256_abs_epi16(b1);\
		c0=_mm256_abs_epi16(c0);\
		c1=_mm256_abs_epi16(c1);\
		a0=_mm256_add_epi16(a0, a1);\
		b0=_mm256_add_epi16(b0, b1);\
		c0=_mm256_add_epi16(c0, c1);\
		a0=_mm256_add_epi16(a0, _mm256_srli_epi64(a0, 32));\
		b0=_mm256_add_epi16(b0, _mm256_srli_epi64(b0, 32));\
		c0=_mm256_add_epi16(c0, _mm256_srli_epi64(c0, 32));\
		a0=_mm256_add_epi16(a0, _mm256_srli_epi64(a0, 16));\
		b0=_mm256_add_epi16(b0, _mm256_srli_epi64(b0, 16));\
		c0=_mm256_add_epi16(c0, _mm256_srli_epi64(c0, 16));\
		mcounters[IDXA]=_mm256_add_epi64(mcounters[IDXA], _mm256_and_si256(a0, wordmask));\
		mcounters[IDXB]=_mm256_add_epi64(mcounters[IDXB], _mm256_and_si256(b0, wordmask));\
		mcounters[IDXC]=_mm256_add_epi64(mcounters[IDXC], _mm256_and_si256(c0, wordmask));\
	}while(0)
					UPDATE(OCH_Y400, OCH_Y040, OCH_Y004, r0, r1, g0, g1, b0, b1);
					UPDATE(OCH_CX40, OCH_C0X4, OCH_C40X, rg0, rg1, gb0, gb1, br0, br1);
					t0=_mm256_add_epi16(rg0, _mm256_srai_epi16(gb0, 2));//r-(3*g+b)/4 = r-g-(b-g)/4
					t1=_mm256_add_epi16(rg1, _mm256_srai_epi16(gb1, 2));
					t2=_mm256_add_epi16(rg0, _mm256_srai_epi16(br0, 2));//g-(3*r+b)/4 = g-r-(b-r)/4
					t3=_mm256_add_epi16(rg1, _mm256_srai_epi16(br1, 2));
					t4=_mm256_add_epi16(br0, _mm256_srai_epi16(rg0, 2));//b-(3*r+g)/4 = b-r-(g-r)/4
					t5=_mm256_add_epi16(br1, _mm256_srai_epi16(rg1, 2));
					UPDATE(OCH_CX31, OCH_C3X1, OCH_C31X, t0, t1, t2, t3, t4, t5);
					t0=_mm256_add_epi16(br0, _mm256_srai_epi16(gb0, 2));//r-(g+3*b)/4 = r-b-(g-b)/4
					t1=_mm256_add_epi16(br1, _mm256_srai_epi16(gb1, 2));
					t2=_mm256_add_epi16(gb0, _mm256_srai_epi16(br0, 2));//g-(r+3*b)/4 = g-b-(r-b)/4
					t3=_mm256_add_epi16(gb1, _mm256_srai_epi16(br1, 2));
					t4=_mm256_add_epi16(gb0, _mm256_srai_epi16(rg0, 2));//b-(r+3*g)/4 = b-g-(r-g)/4
					t5=_mm256_add_epi16(gb1, _mm256_srai_epi16(rg1, 2));
					UPDATE(OCH_CX13, OCH_C1X3, OCH_C13X, t0, t1, t2, t3, t4, t5);
					t0=_mm256_srai_epi16(_mm256_sub_epi16(rg0, br0), 1);//r-(g+b)/2 = (r-g + r-b)/2
					t1=_mm256_srai_epi16(_mm256_sub_epi16(rg1, br1), 1);
					t2=_mm256_srai_epi16(_mm256_sub_epi16(gb0, rg0), 1);//g-(r+b)/2 = (g-r + g-b)/2
					t3=_mm256_srai_epi16(_mm256_sub_epi16(gb1, rg1), 1);
					t4=_mm256_srai_epi16(_mm256_sub_epi16(br0, gb0), 1);//b-(r+g)/2 = (b-r + b-g)/2
					t5=_mm256_srai_epi16(_mm256_sub_epi16(br1, gb1), 1);
					UPDATE(OCH_CX22, OCH_C2X2, OCH_C22X, t0, t1, t2, t3, t4, t5);
				}
				imptr+=rowstride;
			}
			for(int k=0;k<OCH_COUNT;++k)
			{
				ALIGN(32) long long temp[8]={0};
				_mm256_store_si256((__m256i*)temp, mcounters[k]);
				counters[k]=temp[0]+temp[1]+temp[2]+temp[3]+temp[4]+temp[5]+temp[6]+temp[7];
			}
			long long minerr=0;
			for(int kt=0;kt<RCT_COUNT;++kt)
			{
				const unsigned char *rct=rct_combinations[kt];
				long long currerr=
					+counters[rct[0]]
					+counters[rct[1]]
					+counters[rct[2]]
				;
#ifdef LOUD
				//printf("%-14s %12lld + %12lld + %12lld = %12lld\n", rct_names[kt], counters[rct[0]], counters[rct[1]], counters[rct[2]], currerr);
#endif
				if(!kt||minerr>currerr)
				{
					minerr=currerr;
					bestrct=kt;
				}
			}
			const unsigned char *rct=rct_combinations[bestrct];
			long long count=(long long)mainbytes*ih * 2;//FIXME tune
			use_wg4=counters[rct[0]]+counters[rct[1]]+counters[rct[2]]>count;
		}
		prof_checkpoint("analysis");
#ifdef LOUD
		printf("%s  WG4=%d  %td bytes\n", rct_names[bestrct], use_wg4, usize);
#endif
		context=(unsigned char*)malloc(usize+sizeof(__m256i));
		if(!context)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		//generate encode permutations		eg: mask = MSB 0b11000101 LSB  ->  LO {9, 9, 9, 9, 0, 2, 6, 7} HI
		for(int km=0;km<256;++km)
		{
			int *curr=ans_permute+((ptrdiff_t)km<<3);
			int kb2=7;
			for(int kb=7;kb>=0;--kb)
			{
				int bit=km>>kb&1;
				if(bit)
					curr[kb2--]=kb;
			}
		}
		prof_checkpoint("gen permutation");
	}
	else
	{
		//decode flags, stats
		int flags=*streamptr++;
		use_wg4=flags&1;
		bestrct=flags>>1;

		ctxmask=*(unsigned long long*)streamptr;
		streamptr+=8;
		BitPackerLIFO ec;
		bitpacker_dec_init(&ec, streamptr, streamend);
		for(int kc=0;kc<3*NCTX;++kc)
			dec_unpackhist(&ec, CDF2syms+((ptrdiff_t)kc<<PROBBITS), ctxmask, kc);
		dec_unpackhist(&ec, rCDF2syms+((ptrdiff_t)0<<PROBBITS), ctxmask, 3*NCTX+0);
		dec_unpackhist(&ec, rCDF2syms+((ptrdiff_t)1<<PROBBITS), ctxmask, 3*NCTX+1);
		dec_unpackhist(&ec, rCDF2syms+((ptrdiff_t)2<<PROBBITS), ctxmask, 3*NCTX+2);
		streamptr=(unsigned char*)(size_t)ec.srcfwdptr;
		prof_checkpoint("unpack histograms");

		//generate decode permutations		eg: mask = MSB 0b11000101 LSB  ->  LO {0, 9, 1, 9, 9, 9, 2, 3} HI
		for(int km=0;km<256;++km)
		{
			int *curr=ans_permute+((ptrdiff_t)km<<3);
			int idx=0;
			for(int kb=0;kb<8;++kb)
			{
				int bit=km>>kb&1;
				if(bit)
					curr[kb]=idx++;
			}
		}
		prof_checkpoint("gen permutation");
	}
	//use_wg4=0;
	const unsigned char *combination=rct_combinations[bestrct];
	int
		yidx=combination[II_PERM_Y]*NCODERS,
		uidx=combination[II_PERM_U]*NCODERS,
		vidx=combination[II_PERM_V]*NCODERS;
	__m256i uhelpmask=_mm256_set1_epi16(-(combination[II_COEFF_U_SUB_Y]!=0));
	__m256i vc0=_mm256_set1_epi16(combination[II_COEFF_V_SUB_Y]);
	__m256i vc1=_mm256_set1_epi16(combination[II_COEFF_V_SUB_U]);
	int paddedwidth=blockw+16;
	memset(pixels, 0, psize);
	memset(wgerrors, 0, wgsize);
	__m256i mctxmax=_mm256_set1_epi16(NCTX-1);
	__m256i mctxuoffset=_mm256_set1_epi16(NCTX);
	__m256i mctxvoffset=_mm256_set1_epi16(NCTX*2);
	__m256i amin=_mm256_set1_epi16(-128);
	__m256i amax=_mm256_set1_epi16(127);
	__m128i half8=_mm_set1_epi8(-128);
	__m256i bytemask=_mm256_set1_epi16(255);
	__m256i wordmask=_mm256_set1_epi32(0xFFFF);
	__m256i myuv[6];
	memset(myuv, 0, sizeof(myuv));
	imptr=image;
	unsigned char *ctxptr=fwd?context:0;
	__m256i mstate[4];
	short *wgWerrors=wgstate+0*(ptrdiff_t)NCODERS*3*WG_NPREDS;
	__m256i *wgpreds=(__m256i*)(wgstate+1*(ptrdiff_t)NCODERS*3*WG_NPREDS);
	memset(wgstate, 0, wgstatesize);
	if(!fwd)
	{
#ifdef _DEBUG
		if(streamptr>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		memcpy(mstate, streamptr, sizeof(mstate));
		streamptr+=sizeof(mstate);
	}
	for(int ky=0;ky<ih;++ky)//main coding loop
	{
		ALIGN(32) short *rows[]=
		{
			pixels+(paddedwidth*((ky-0LL)&3)+8LL)*6*NCODERS,
			pixels+(paddedwidth*((ky-1LL)&3)+8LL)*6*NCODERS,
			pixels+(paddedwidth*((ky-2LL)&3)+8LL)*6*NCODERS,
			pixels+(paddedwidth*((ky-3LL)&3)+8LL)*6*NCODERS,
		};
		ALIGN(32) short *erows[]=
		{
			wgerrors+(paddedwidth*((ky-0LL)&1)+8LL)*NCODERS*3*WG_NPREDS,
			wgerrors+(paddedwidth*((ky-1LL)&1)+8LL)*NCODERS*3*WG_NPREDS,
		};
		ALIGN(32) unsigned short syms[32]={0};
		__m256i NW[6], N[6], W[6];
		__m256i eW[6], ecurr[6], eNEE[6], eNEEE[6];
		memset(NW, 0, sizeof(NW));
		memset(N, 0, sizeof(N));
		memset(W, 0, sizeof(W));
		memset(eW, 0, sizeof(eW));
		memset(ecurr, 0, sizeof(ecurr));
		memset(eNEE, 0, sizeof(eNEE));
		memset(eNEEE, 0, sizeof(eNEEE));
		memset(wgpreds, 0, sizeof(short[NCODERS*3*WG_NPREDS]));
		//                       (__m256i*)(rows[-Y]+(E+C+X*6)*NCODERS)+R
		eNEE[0]=_mm256_load_si256((__m256i*)(rows[1]+(3+0+2*6)*NCODERS)+0);
		eNEE[1]=_mm256_load_si256((__m256i*)(rows[1]+(3+0+2*6)*NCODERS)+1);
		eNEE[2]=_mm256_load_si256((__m256i*)(rows[1]+(3+1+2*6)*NCODERS)+0);
		eNEE[3]=_mm256_load_si256((__m256i*)(rows[1]+(3+1+2*6)*NCODERS)+1);
		eNEE[4]=_mm256_load_si256((__m256i*)(rows[1]+(3+2+2*6)*NCODERS)+0);
		eNEE[5]=_mm256_load_si256((__m256i*)(rows[1]+(3+2+2*6)*NCODERS)+1);
		for(int kx=0;kx<rowstride-3*NCODERS;kx+=3*NCODERS)
		{
			//						(E+C+X*6)*NCODERS
			N	[0]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+0*6)*NCODERS)+0);//y0 neighbors
			N	[1]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+0*6)*NCODERS)+1);//y1
			N	[2]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+0*6)*NCODERS)+0);//u0
			N	[3]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+0*6)*NCODERS)+1);//u1
			N	[4]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+0*6)*NCODERS)+0);//v0
			N	[5]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+0*6)*NCODERS)+1);//v1
			__m256i
				predY0, predY1, ctxY0, ctxY1,
				predU0, predU1, ctxU0, ctxU1,
				predV0, predV1, ctxV0, ctxV1;
			{
				__m256i one=_mm256_set1_epi32(1);
				__m256i cy0=_mm256_and_si256(eW[0], wordmask), cy1=_mm256_srli_epi32(eW[0], 16);
				__m256i cy2=_mm256_and_si256(eW[1], wordmask), cy3=_mm256_srli_epi32(eW[1], 16);
				__m256i cu0=_mm256_and_si256(eW[2], wordmask), cu1=_mm256_srli_epi32(eW[2], 16);
				__m256i cu2=_mm256_and_si256(eW[3], wordmask), cu3=_mm256_srli_epi32(eW[3], 16);
				__m256i cv0=_mm256_and_si256(eW[4], wordmask), cv1=_mm256_srli_epi32(eW[4], 16);
				__m256i cv2=_mm256_and_si256(eW[5], wordmask), cv3=_mm256_srli_epi32(eW[5], 16);
				cy0=_mm256_mullo_epi32(cy0, cy0);
				cy1=_mm256_mullo_epi32(cy1, cy1);
				cy2=_mm256_mullo_epi32(cy2, cy2);
				cy3=_mm256_mullo_epi32(cy3, cy3);
				cu0=_mm256_mullo_epi32(cu0, cu0);
				cu1=_mm256_mullo_epi32(cu1, cu1);
				cu2=_mm256_mullo_epi32(cu2, cu2);
				cu3=_mm256_mullo_epi32(cu3, cu3);
				cv0=_mm256_mullo_epi32(cv0, cv0);
				cv1=_mm256_mullo_epi32(cv1, cv1);
				cv2=_mm256_mullo_epi32(cv2, cv2);
				cv3=_mm256_mullo_epi32(cv3, cv3);
				cy0=_mm256_add_epi32(cy0, one);
				cy1=_mm256_add_epi32(cy1, one);
				cy2=_mm256_add_epi32(cy2, one);
				cy3=_mm256_add_epi32(cy3, one);
				cu0=_mm256_add_epi32(cu0, one);
				cu1=_mm256_add_epi32(cu1, one);
				cu2=_mm256_add_epi32(cu2, one);
				cu3=_mm256_add_epi32(cu3, one);
				cv0=_mm256_add_epi32(cv0, one);
				cv1=_mm256_add_epi32(cv1, one);
				cv2=_mm256_add_epi32(cv2, one);
				cv3=_mm256_add_epi32(cv3, one);
				//FLOOR_LOG2_32x8(X) = _mm256_sub_epi32(_mm256_srli_epi32(_mm256_castps_si256(_mm256_cvtepi32_ps(X)), 23), _mm256_set1_epi32(127))
				cy0=_mm256_castps_si256(_mm256_cvtepi32_ps(cy0));
				cy1=_mm256_castps_si256(_mm256_cvtepi32_ps(cy1));
				cy2=_mm256_castps_si256(_mm256_cvtepi32_ps(cy2));
				cy3=_mm256_castps_si256(_mm256_cvtepi32_ps(cy3));
				cu0=_mm256_castps_si256(_mm256_cvtepi32_ps(cu0));
				cu1=_mm256_castps_si256(_mm256_cvtepi32_ps(cu1));
				cu2=_mm256_castps_si256(_mm256_cvtepi32_ps(cu2));
				cu3=_mm256_castps_si256(_mm256_cvtepi32_ps(cu3));
				cv0=_mm256_castps_si256(_mm256_cvtepi32_ps(cv0));
				cv1=_mm256_castps_si256(_mm256_cvtepi32_ps(cv1));
				cv2=_mm256_castps_si256(_mm256_cvtepi32_ps(cv2));
				cv3=_mm256_castps_si256(_mm256_cvtepi32_ps(cv3));
				cy0=_mm256_srli_epi32(cy0, 23);
				cy1=_mm256_srli_epi32(cy1, 23);
				cy2=_mm256_srli_epi32(cy2, 23);
				cy3=_mm256_srli_epi32(cy3, 23);
				cu0=_mm256_srli_epi32(cu0, 23);
				cu1=_mm256_srli_epi32(cu1, 23);
				cu2=_mm256_srli_epi32(cu2, 23);
				cu3=_mm256_srli_epi32(cu3, 23);
				cv0=_mm256_srli_epi32(cv0, 23);
				cv1=_mm256_srli_epi32(cv1, 23);
				cv2=_mm256_srli_epi32(cv2, 23);
				cv3=_mm256_srli_epi32(cv3, 23);
				__m256i expbias=_mm256_set1_epi32(127);
				cy0=_mm256_sub_epi32(cy0, expbias);
				cy1=_mm256_sub_epi32(cy1, expbias);
				cy2=_mm256_sub_epi32(cy2, expbias);
				cy3=_mm256_sub_epi32(cy3, expbias);
				cu0=_mm256_sub_epi32(cu0, expbias);
				cu1=_mm256_sub_epi32(cu1, expbias);
				cu2=_mm256_sub_epi32(cu2, expbias);
				cu3=_mm256_sub_epi32(cu3, expbias);
				cv0=_mm256_sub_epi32(cv0, expbias);
				cv1=_mm256_sub_epi32(cv1, expbias);
				cv2=_mm256_sub_epi32(cv2, expbias);
				cv3=_mm256_sub_epi32(cv3, expbias);
				cy1=_mm256_slli_epi32(cy1, 16);
				cy3=_mm256_slli_epi32(cy3, 16);
				cu1=_mm256_slli_epi32(cu1, 16);
				cu3=_mm256_slli_epi32(cu3, 16);
				cv1=_mm256_slli_epi32(cv1, 16);
				cv3=_mm256_slli_epi32(cv3, 16);
				ctxY0=_mm256_or_si256(cy0, cy1);
				ctxY1=_mm256_or_si256(cy2, cy3);
				ctxU0=_mm256_or_si256(cu0, cu1);
				ctxU1=_mm256_or_si256(cu2, cu3);
				ctxV0=_mm256_or_si256(cv0, cv1);
				ctxV1=_mm256_or_si256(cv2, cv3);
				ctxY0=_mm256_min_epi16(ctxY0, mctxmax);
				ctxY1=_mm256_min_epi16(ctxY1, mctxmax);
				ctxU0=_mm256_min_epi16(ctxU0, mctxmax);
				ctxU1=_mm256_min_epi16(ctxU1, mctxmax);
				ctxV0=_mm256_min_epi16(ctxV0, mctxmax);
				ctxV1=_mm256_min_epi16(ctxV1, mctxmax);
			}
			{
				__m256i ymin0=_mm256_min_epi16(N[0], W[0]);
				__m256i ymax0=_mm256_max_epi16(N[0], W[0]);
				__m256i ymin1=_mm256_min_epi16(N[1], W[1]);
				__m256i ymax1=_mm256_max_epi16(N[1], W[1]);
				__m256i umin0=_mm256_min_epi16(N[2], W[2]);
				__m256i umax0=_mm256_max_epi16(N[2], W[2]);
				__m256i umin1=_mm256_min_epi16(N[3], W[3]);
				__m256i umax1=_mm256_max_epi16(N[3], W[3]);
				__m256i vmin0=_mm256_min_epi16(N[4], W[4]);
				__m256i vmax0=_mm256_max_epi16(N[4], W[4]);
				__m256i vmin1=_mm256_min_epi16(N[5], W[5]);
				__m256i vmax1=_mm256_max_epi16(N[5], W[5]);
				predY0=_mm256_add_epi16(N[0], W[0]);
				predY1=_mm256_add_epi16(N[1], W[1]);
				predU0=_mm256_add_epi16(N[2], W[2]);
				predU1=_mm256_add_epi16(N[3], W[3]);
				predV0=_mm256_add_epi16(N[4], W[4]);
				predV1=_mm256_add_epi16(N[5], W[5]);
				predY0=_mm256_sub_epi16(predY0, NW[0]);
				predY1=_mm256_sub_epi16(predY1, NW[1]);
				predU0=_mm256_sub_epi16(predU0, NW[2]);
				predU1=_mm256_sub_epi16(predU1, NW[3]);
				predV0=_mm256_sub_epi16(predV0, NW[4]);
				predV1=_mm256_sub_epi16(predV1, NW[5]);
				if(use_wg4)//predict
				{
					/*
					0	N
					1	W
					2	3*(N-NN)+NNN
					3	3*(W-WW)+WWW
					4	W+NE-N
					5	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4
					6	N+W-NW
					7	N+NE-NNE
					*/
					__m256i cache[6];

					//N
					wgpreds[0*6+0]=N[0];
					wgpreds[0*6+1]=N[1];
					wgpreds[0*6+2]=N[2];
					wgpreds[0*6+3]=N[3];
					wgpreds[0*6+4]=N[4];
					wgpreds[0*6+5]=N[5];

					//W
					wgpreds[1*6+0]=W[0];
					wgpreds[1*6+1]=W[1];
					wgpreds[1*6+2]=W[2];
					wgpreds[1*6+3]=W[3];
					wgpreds[1*6+4]=W[4];
					wgpreds[1*6+5]=W[5];

					//3*(N-NN)+NNN
					cache[0]=_mm256_load_si256((__m256i*)(rows[2]+(0+0+0*6)*NCODERS)+0);//NN
					cache[1]=_mm256_load_si256((__m256i*)(rows[2]+(0+0+0*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[2]+(0+1+0*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[2]+(0+1+0*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[2]+(0+2+0*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[2]+(0+2+0*6)*NCODERS)+1);
					wgpreds[2*6+0]=_mm256_sub_epi16(N[0], cache[0]);
					wgpreds[2*6+1]=_mm256_sub_epi16(N[1], cache[1]);
					wgpreds[2*6+2]=_mm256_sub_epi16(N[2], cache[2]);
					wgpreds[2*6+3]=_mm256_sub_epi16(N[3], cache[3]);
					wgpreds[2*6+4]=_mm256_sub_epi16(N[4], cache[4]);
					wgpreds[2*6+5]=_mm256_sub_epi16(N[5], cache[5]);
					cache[0]=_mm256_slli_epi16(wgpreds[2*6+0], 1);
					cache[1]=_mm256_slli_epi16(wgpreds[2*6+1], 1);
					cache[2]=_mm256_slli_epi16(wgpreds[2*6+2], 1);
					cache[3]=_mm256_slli_epi16(wgpreds[2*6+3], 1);
					cache[4]=_mm256_slli_epi16(wgpreds[2*6+4], 1);
					cache[5]=_mm256_slli_epi16(wgpreds[2*6+5], 1);
					wgpreds[2*6+0]=_mm256_add_epi16(wgpreds[2*6+0], cache[0]);
					wgpreds[2*6+1]=_mm256_add_epi16(wgpreds[2*6+1], cache[1]);
					wgpreds[2*6+2]=_mm256_add_epi16(wgpreds[2*6+2], cache[2]);
					wgpreds[2*6+3]=_mm256_add_epi16(wgpreds[2*6+3], cache[3]);
					wgpreds[2*6+4]=_mm256_add_epi16(wgpreds[2*6+4], cache[4]);
					wgpreds[2*6+5]=_mm256_add_epi16(wgpreds[2*6+5], cache[5]);
					cache[0]=_mm256_load_si256((__m256i*)(rows[3]+(0+0+0*6)*NCODERS)+0);//NNN
					cache[1]=_mm256_load_si256((__m256i*)(rows[3]+(0+0+0*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[3]+(0+1+0*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[3]+(0+1+0*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[3]+(0+2+0*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[3]+(0+2+0*6)*NCODERS)+1);
					wgpreds[2*6+0]=_mm256_add_epi16(wgpreds[2*6+0], cache[0]);
					wgpreds[2*6+1]=_mm256_add_epi16(wgpreds[2*6+1], cache[1]);
					wgpreds[2*6+2]=_mm256_add_epi16(wgpreds[2*6+2], cache[2]);
					wgpreds[2*6+3]=_mm256_add_epi16(wgpreds[2*6+3], cache[3]);
					wgpreds[2*6+4]=_mm256_add_epi16(wgpreds[2*6+4], cache[4]);
					wgpreds[2*6+5]=_mm256_add_epi16(wgpreds[2*6+5], cache[5]);

					//3*(W-WW)+WWW
					cache[0]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-2*6)*NCODERS)+0);//WW
					cache[1]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-2*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-2*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-2*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-2*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-2*6)*NCODERS)+1);
					wgpreds[3*6+0]=_mm256_sub_epi16(W[0], cache[0]);
					wgpreds[3*6+1]=_mm256_sub_epi16(W[1], cache[1]);
					wgpreds[3*6+2]=_mm256_sub_epi16(W[2], cache[2]);
					wgpreds[3*6+3]=_mm256_sub_epi16(W[3], cache[3]);
					wgpreds[3*6+4]=_mm256_sub_epi16(W[4], cache[4]);
					wgpreds[3*6+5]=_mm256_sub_epi16(W[5], cache[5]);
					cache[0]=_mm256_slli_epi16(wgpreds[3*6+0], 1);
					cache[1]=_mm256_slli_epi16(wgpreds[3*6+1], 1);
					cache[2]=_mm256_slli_epi16(wgpreds[3*6+2], 1);
					cache[3]=_mm256_slli_epi16(wgpreds[3*6+3], 1);
					cache[4]=_mm256_slli_epi16(wgpreds[3*6+4], 1);
					cache[5]=_mm256_slli_epi16(wgpreds[3*6+5], 1);
					wgpreds[3*6+0]=_mm256_add_epi16(wgpreds[3*6+0], cache[0]);
					wgpreds[3*6+1]=_mm256_add_epi16(wgpreds[3*6+1], cache[1]);
					wgpreds[3*6+2]=_mm256_add_epi16(wgpreds[3*6+2], cache[2]);
					wgpreds[3*6+3]=_mm256_add_epi16(wgpreds[3*6+3], cache[3]);
					wgpreds[3*6+4]=_mm256_add_epi16(wgpreds[3*6+4], cache[4]);
					wgpreds[3*6+5]=_mm256_add_epi16(wgpreds[3*6+5], cache[5]);
					cache[0]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-3*6)*NCODERS)+0);//WWW
					cache[1]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-3*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-3*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-3*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-3*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-3*6)*NCODERS)+1);
					wgpreds[3*6+0]=_mm256_add_epi16(wgpreds[3*6+0], cache[0]);
					wgpreds[3*6+1]=_mm256_add_epi16(wgpreds[3*6+1], cache[1]);
					wgpreds[3*6+2]=_mm256_add_epi16(wgpreds[3*6+2], cache[2]);
					wgpreds[3*6+3]=_mm256_add_epi16(wgpreds[3*6+3], cache[3]);
					wgpreds[3*6+4]=_mm256_add_epi16(wgpreds[3*6+4], cache[4]);
					wgpreds[3*6+5]=_mm256_add_epi16(wgpreds[3*6+5], cache[5]);

					//W+NE-N
					cache[0]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+0);//NE
					cache[1]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+1);
					wgpreds[4*6+0]=_mm256_sub_epi16(W[0], N[0]);
					wgpreds[4*6+1]=_mm256_sub_epi16(W[1], N[1]);
					wgpreds[4*6+2]=_mm256_sub_epi16(W[2], N[2]);
					wgpreds[4*6+3]=_mm256_sub_epi16(W[3], N[3]);
					wgpreds[4*6+4]=_mm256_sub_epi16(W[4], N[4]);
					wgpreds[4*6+5]=_mm256_sub_epi16(W[5], N[5]);
					wgpreds[4*6+0]=_mm256_add_epi16(wgpreds[4*6+0], cache[0]);
					wgpreds[4*6+1]=_mm256_add_epi16(wgpreds[4*6+1], cache[1]);
					wgpreds[4*6+2]=_mm256_add_epi16(wgpreds[4*6+2], cache[2]);
					wgpreds[4*6+3]=_mm256_add_epi16(wgpreds[4*6+3], cache[3]);
					wgpreds[4*6+4]=_mm256_add_epi16(wgpreds[4*6+4], cache[4]);
					wgpreds[4*6+5]=_mm256_add_epi16(wgpreds[4*6+5], cache[5]);

					//(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4
					wgpreds[5*6+0]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-4*6)*NCODERS)+0);//WWWW
					wgpreds[5*6+1]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-4*6)*NCODERS)+1);
					wgpreds[5*6+2]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-4*6)*NCODERS)+0);
					wgpreds[5*6+3]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-4*6)*NCODERS)+1);
					wgpreds[5*6+4]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-4*6)*NCODERS)+0);
					wgpreds[5*6+5]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-4*6)*NCODERS)+1);
					cache[0]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-3*6)*NCODERS)+0);//WWW
					cache[1]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-3*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-3*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-3*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-3*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-3*6)*NCODERS)+1);
					wgpreds[5*6+0]=_mm256_add_epi16(wgpreds[5*6+0], cache[0]);
					wgpreds[5*6+1]=_mm256_add_epi16(wgpreds[5*6+1], cache[1]);
					wgpreds[5*6+2]=_mm256_add_epi16(wgpreds[5*6+2], cache[2]);
					wgpreds[5*6+3]=_mm256_add_epi16(wgpreds[5*6+3], cache[3]);
					wgpreds[5*6+4]=_mm256_add_epi16(wgpreds[5*6+4], cache[4]);
					wgpreds[5*6+5]=_mm256_add_epi16(wgpreds[5*6+5], cache[5]);
					cache[0]=_mm256_load_si256((__m256i*)(rows[3]+(0+0+0*6)*NCODERS)+0);//NNN
					cache[1]=_mm256_load_si256((__m256i*)(rows[3]+(0+0+0*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[3]+(0+1+0*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[3]+(0+1+0*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[3]+(0+2+0*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[3]+(0+2+0*6)*NCODERS)+1);
					wgpreds[5*6+0]=_mm256_add_epi16(wgpreds[5*6+0], cache[0]);
					wgpreds[5*6+1]=_mm256_add_epi16(wgpreds[5*6+1], cache[1]);
					wgpreds[5*6+2]=_mm256_add_epi16(wgpreds[5*6+2], cache[2]);
					wgpreds[5*6+3]=_mm256_add_epi16(wgpreds[5*6+3], cache[3]);
					wgpreds[5*6+4]=_mm256_add_epi16(wgpreds[5*6+4], cache[4]);
					wgpreds[5*6+5]=_mm256_add_epi16(wgpreds[5*6+5], cache[5]);
					cache[0]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+2*6)*NCODERS)+0);//NEE
					cache[1]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+2*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+2*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+2*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+2*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+2*6)*NCODERS)+1);
					wgpreds[5*6+0]=_mm256_add_epi16(wgpreds[5*6+0], cache[0]);
					wgpreds[5*6+1]=_mm256_add_epi16(wgpreds[5*6+1], cache[1]);
					wgpreds[5*6+2]=_mm256_add_epi16(wgpreds[5*6+2], cache[2]);
					wgpreds[5*6+3]=_mm256_add_epi16(wgpreds[5*6+3], cache[3]);
					wgpreds[5*6+4]=_mm256_add_epi16(wgpreds[5*6+4], cache[4]);
					wgpreds[5*6+5]=_mm256_add_epi16(wgpreds[5*6+5], cache[5]);
					cache[0]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+3*6)*NCODERS)+0);//NEEE
					cache[1]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+3*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+3*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+3*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+3*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+3*6)*NCODERS)+1);
					wgpreds[5*6+0]=_mm256_add_epi16(wgpreds[5*6+0], cache[0]);
					wgpreds[5*6+1]=_mm256_add_epi16(wgpreds[5*6+1], cache[1]);
					wgpreds[5*6+2]=_mm256_add_epi16(wgpreds[5*6+2], cache[2]);
					wgpreds[5*6+3]=_mm256_add_epi16(wgpreds[5*6+3], cache[3]);
					wgpreds[5*6+4]=_mm256_add_epi16(wgpreds[5*6+4], cache[4]);
					wgpreds[5*6+5]=_mm256_add_epi16(wgpreds[5*6+5], cache[5]);
					cache[0]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+4*6)*NCODERS)+0);//NEEEE
					cache[1]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+4*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+4*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+4*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+4*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+4*6)*NCODERS)+1);
					wgpreds[5*6+0]=_mm256_add_epi16(wgpreds[5*6+0], cache[0]);
					wgpreds[5*6+1]=_mm256_add_epi16(wgpreds[5*6+1], cache[1]);
					wgpreds[5*6+2]=_mm256_add_epi16(wgpreds[5*6+2], cache[2]);
					wgpreds[5*6+3]=_mm256_add_epi16(wgpreds[5*6+3], cache[3]);
					wgpreds[5*6+4]=_mm256_add_epi16(wgpreds[5*6+4], cache[4]);
					wgpreds[5*6+5]=_mm256_add_epi16(wgpreds[5*6+5], cache[5]);
					cache[0]=_mm256_slli_epi16(NW[0], 1);
					cache[1]=_mm256_slli_epi16(NW[1], 1);
					cache[2]=_mm256_slli_epi16(NW[2], 1);
					cache[3]=_mm256_slli_epi16(NW[3], 1);
					cache[4]=_mm256_slli_epi16(NW[4], 1);
					cache[5]=_mm256_slli_epi16(NW[5], 1);
					wgpreds[5*6+0]=_mm256_sub_epi16(wgpreds[5*6+0], cache[0]);
					wgpreds[5*6+1]=_mm256_sub_epi16(wgpreds[5*6+1], cache[1]);
					wgpreds[5*6+2]=_mm256_sub_epi16(wgpreds[5*6+2], cache[2]);
					wgpreds[5*6+3]=_mm256_sub_epi16(wgpreds[5*6+3], cache[3]);
					wgpreds[5*6+4]=_mm256_sub_epi16(wgpreds[5*6+4], cache[4]);
					wgpreds[5*6+5]=_mm256_sub_epi16(wgpreds[5*6+5], cache[5]);
					wgpreds[5*6+0]=_mm256_srai_epi16(wgpreds[5*6+0], 2);
					wgpreds[5*6+1]=_mm256_srai_epi16(wgpreds[5*6+1], 2);
					wgpreds[5*6+2]=_mm256_srai_epi16(wgpreds[5*6+2], 2);
					wgpreds[5*6+3]=_mm256_srai_epi16(wgpreds[5*6+3], 2);
					wgpreds[5*6+4]=_mm256_srai_epi16(wgpreds[5*6+4], 2);
					wgpreds[5*6+5]=_mm256_srai_epi16(wgpreds[5*6+5], 2);

					//N+W-NW
					wgpreds[6*6+0]=predY0;
					wgpreds[6*6+1]=predY1;
					wgpreds[6*6+2]=predU0;
					wgpreds[6*6+3]=predU1;
					wgpreds[6*6+4]=predV0;
					wgpreds[6*6+5]=predV1;

					//N+NE-NNE
					cache[0]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+0);//NE
					cache[1]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+1);
					wgpreds[7*6+0]=_mm256_add_epi16(N[0], cache[0]);
					wgpreds[7*6+1]=_mm256_add_epi16(N[1], cache[1]);
					wgpreds[7*6+2]=_mm256_add_epi16(N[2], cache[2]);
					wgpreds[7*6+3]=_mm256_add_epi16(N[3], cache[3]);
					wgpreds[7*6+4]=_mm256_add_epi16(N[4], cache[4]);
					wgpreds[7*6+5]=_mm256_add_epi16(N[5], cache[5]);
					cache[0]=_mm256_load_si256((__m256i*)(rows[2]+(0+0+1*6)*NCODERS)+0);//NNE
					cache[1]=_mm256_load_si256((__m256i*)(rows[2]+(0+0+1*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[2]+(0+1+1*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[2]+(0+1+1*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[2]+(0+2+1*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[2]+(0+2+1*6)*NCODERS)+1);
					wgpreds[7*6+0]=_mm256_sub_epi16(wgpreds[7*6+0], cache[0]);
					wgpreds[7*6+1]=_mm256_sub_epi16(wgpreds[7*6+1], cache[1]);
					wgpreds[7*6+2]=_mm256_sub_epi16(wgpreds[7*6+2], cache[2]);
					wgpreds[7*6+3]=_mm256_sub_epi16(wgpreds[7*6+3], cache[3]);
					wgpreds[7*6+4]=_mm256_sub_epi16(wgpreds[7*6+4], cache[4]);
					wgpreds[7*6+5]=_mm256_sub_epi16(wgpreds[7*6+5], cache[5]);

					//mix
					__m256i result[6];
					wg_mix(0, wgWerrors, erows[0], erows[1], wgpreds, result+0);
					wg_mix(1, wgWerrors, erows[0], erows[1], wgpreds, result+2);
					wg_mix(2, wgWerrors, erows[0], erows[1], wgpreds, result+4);
					predY0=result[0];
					predY1=result[1];
					predU0=result[2];
					predU1=result[3];
					predV0=result[4];
					predV1=result[5];

#if 0
					__m256d de[16];//de[P*2+highhalf]		prederrors[P*2+halflaneidx]
					de[0*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[0*2+0], 1));
					de[1*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[1*2+0], 1));
					de[2*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[2*2+0], 1));
					de[3*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[3*2+0], 1));
					de[4*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[4*2+0], 1));
					de[5*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[5*2+0], 1));
					de[6*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[6*2+0], 1));
					de[7*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[7*2+0], 1));
					de[0*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[0*2+0]));
					de[1*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[1*2+0]));
					de[2*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[2*2+0]));
					de[3*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[3*2+0]));
					de[4*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[4*2+0]));
					de[5*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[5*2+0]));
					de[6*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[6*2+0]));
					de[7*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[7*2+0]));
#endif

					//loosen pred range
					cache[0]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+0);//NE
					cache[1]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+1);
					ymin0=_mm256_min_epi16(ymin0, cache[0]);
					ymax0=_mm256_max_epi16(ymax0, cache[0]);
					ymin1=_mm256_min_epi16(ymin1, cache[1]);
					ymax1=_mm256_max_epi16(ymax1, cache[1]);
					umin0=_mm256_min_epi16(umin0, cache[2]);
					umax0=_mm256_max_epi16(umax0, cache[2]);
					umin1=_mm256_min_epi16(umin1, cache[3]);
					umax1=_mm256_max_epi16(umax1, cache[3]);
					vmin0=_mm256_min_epi16(vmin0, cache[4]);
					vmax0=_mm256_max_epi16(vmax0, cache[4]);
					vmin1=_mm256_min_epi16(vmin1, cache[5]);
					vmax1=_mm256_max_epi16(vmax1, cache[5]);
				}
				predY0=_mm256_max_epi16(predY0, ymin0);
				predY1=_mm256_max_epi16(predY1, ymin1);
				predU0=_mm256_max_epi16(predU0, umin0);
				predU1=_mm256_max_epi16(predU1, umin1);
				predV0=_mm256_max_epi16(predV0, vmin0);
				predV1=_mm256_max_epi16(predV1, vmin1);
				predY0=_mm256_min_epi16(predY0, ymax0);
				predY1=_mm256_min_epi16(predY1, ymax1);
				predU0=_mm256_min_epi16(predU0, umax0);
				predU1=_mm256_min_epi16(predU1, umax1);
				predV0=_mm256_min_epi16(predV0, vmax0);
				predV1=_mm256_min_epi16(predV1, vmax1);
			}
			__m256i msyms0, msyms1, msyms8, mctx, moffset0, moffset1;
			if(fwd)
			{
				myuv[0]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx+yidx)+0), half8));//load yuv
				myuv[1]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx+yidx)+1), half8));
				myuv[2]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx+uidx)+0), half8));
				myuv[3]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx+uidx)+1), half8));
				myuv[4]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx+vidx)+0), half8));
				myuv[5]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx+vidx)+1), half8));
				
				msyms0=_mm256_sub_epi16(myuv[0], predY0);//sub pred
				msyms1=_mm256_sub_epi16(myuv[1], predY1);
				ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));//ecurr = pack_sign(yuv-pred)		FIXME try abs
				ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms0=_mm256_sub_epi16(msyms0, amin);
				msyms1=_mm256_sub_epi16(msyms1, amin);
				msyms0=_mm256_and_si256(msyms0, bytemask);//sym = (yuv-pred+128)&255
				msyms1=_mm256_and_si256(msyms1, bytemask);
				msyms8=_mm256_packus_epi16(msyms0, msyms1);//FIXME optimize: shift-or
				mctx=_mm256_packus_epi16(ctxY0, ctxY1);
				ctxY0=_mm256_slli_epi16(ctxY0, 8);
				ctxY1=_mm256_slli_epi16(ctxY1, 8);
				ctxY0=_mm256_or_si256(ctxY0, msyms0);
				ctxY1=_mm256_or_si256(ctxY1, msyms1);
				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
				mctx=_mm256_permute4x64_epi64(mctx, _MM_SHUFFLE(3, 1, 2, 0));
				_mm256_storeu_si256((__m256i*)(imptr+kx+yidx), msyms8);//store residuals
				_mm256_storeu_si256((__m256i*)(ctxptr+kx+yidx), mctx);//store contexts
				_mm256_store_si256((__m256i*)syms+0, ctxY0);
				_mm256_store_si256((__m256i*)syms+1, ctxY1);
				W[0]=myuv[0];
				W[1]=myuv[1];
				_mm256_store_si256((__m256i*)(rows[0]+(0+0+0*6)*NCODERS)+0, myuv[0]);//store Y neighbors
				_mm256_store_si256((__m256i*)(rows[0]+(0+0+0*6)*NCODERS)+1, myuv[1]);
				++hists[syms[0x00]];
				++hists[syms[0x01]];
				++hists[syms[0x02]];
				++hists[syms[0x03]];
				++hists[syms[0x04]];
				++hists[syms[0x05]];
				++hists[syms[0x06]];
				++hists[syms[0x07]];
				++hists[syms[0x08]];
				++hists[syms[0x09]];
				++hists[syms[0x0A]];
				++hists[syms[0x0B]];
				++hists[syms[0x0C]];
				++hists[syms[0x0D]];
				++hists[syms[0x0E]];
				++hists[syms[0x0F]];
				++hists[syms[0x10]];
				++hists[syms[0x11]];
				++hists[syms[0x12]];
				++hists[syms[0x13]];
				++hists[syms[0x14]];
				++hists[syms[0x15]];
				++hists[syms[0x16]];
				++hists[syms[0x17]];
				++hists[syms[0x18]];
				++hists[syms[0x19]];
				++hists[syms[0x1A]];
				++hists[syms[0x1B]];
				++hists[syms[0x1C]];
				++hists[syms[0x1D]];
				++hists[syms[0x1E]];
				++hists[syms[0x1F]];

				moffset0=_mm256_and_si256(myuv[0], uhelpmask);
				moffset1=_mm256_and_si256(myuv[1], uhelpmask);
				predU0=_mm256_add_epi16(predU0, moffset0);
				predU1=_mm256_add_epi16(predU1, moffset1);
				predU0=_mm256_max_epi16(predU0, amin);
				predU1=_mm256_max_epi16(predU1, amin);
				predU0=_mm256_min_epi16(predU0, amax);
				predU1=_mm256_min_epi16(predU1, amax);
				
				msyms0=_mm256_sub_epi16(myuv[2], predU0);
				msyms1=_mm256_sub_epi16(myuv[3], predU1);
				ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));
				ecurr[3]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms0=_mm256_sub_epi16(msyms0, amin);
				msyms1=_mm256_sub_epi16(msyms1, amin);
				msyms0=_mm256_and_si256(msyms0, bytemask);
				msyms1=_mm256_and_si256(msyms1, bytemask);
				msyms8=_mm256_packus_epi16(msyms0, msyms1);
				mctx=_mm256_packus_epi16(ctxU0, ctxU1);
				ctxU0=_mm256_add_epi16(ctxU0, mctxuoffset);
				ctxU1=_mm256_add_epi16(ctxU1, mctxuoffset);
				ctxU0=_mm256_slli_epi16(ctxU0, 8);
				ctxU1=_mm256_slli_epi16(ctxU1, 8);
				ctxU0=_mm256_or_si256(ctxU0, msyms0);
				ctxU1=_mm256_or_si256(ctxU1, msyms1);
				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
				mctx=_mm256_permute4x64_epi64(mctx, _MM_SHUFFLE(3, 1, 2, 0));
				_mm256_storeu_si256((__m256i*)(imptr+kx+uidx), msyms8);
				_mm256_storeu_si256((__m256i*)(ctxptr+kx+uidx), mctx);
				msyms0=_mm256_sub_epi16(myuv[2], moffset0);
				msyms1=_mm256_sub_epi16(myuv[3], moffset1);
				_mm256_store_si256((__m256i*)syms+0, ctxU0);
				_mm256_store_si256((__m256i*)syms+1, ctxU1);
				W[2]=msyms0;
				W[3]=msyms1;
				_mm256_store_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+0, msyms0);//store U neighbors
				_mm256_store_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+1, msyms1);
				++hists[syms[0x00]];
				++hists[syms[0x01]];
				++hists[syms[0x02]];
				++hists[syms[0x03]];
				++hists[syms[0x04]];
				++hists[syms[0x05]];
				++hists[syms[0x06]];
				++hists[syms[0x07]];
				++hists[syms[0x08]];
				++hists[syms[0x09]];
				++hists[syms[0x0A]];
				++hists[syms[0x0B]];
				++hists[syms[0x0C]];
				++hists[syms[0x0D]];
				++hists[syms[0x0E]];
				++hists[syms[0x0F]];
				++hists[syms[0x10]];
				++hists[syms[0x11]];
				++hists[syms[0x12]];
				++hists[syms[0x13]];
				++hists[syms[0x14]];
				++hists[syms[0x15]];
				++hists[syms[0x16]];
				++hists[syms[0x17]];
				++hists[syms[0x18]];
				++hists[syms[0x19]];
				++hists[syms[0x1A]];
				++hists[syms[0x1B]];
				++hists[syms[0x1C]];
				++hists[syms[0x1D]];
				++hists[syms[0x1E]];
				++hists[syms[0x1F]];

				moffset0=_mm256_mullo_epi16(vc0, myuv[0]);
				moffset1=_mm256_mullo_epi16(vc0, myuv[1]);
				moffset0=_mm256_add_epi16(moffset0, _mm256_mullo_epi16(vc1, myuv[2]));
				moffset1=_mm256_add_epi16(moffset1, _mm256_mullo_epi16(vc1, myuv[3]));
				predV0=_mm256_add_epi16(predV0, moffset0);
				predV1=_mm256_add_epi16(predV1, moffset1);
				predV0=_mm256_srai_epi16(predV0, 2);
				predV1=_mm256_srai_epi16(predV1, 2);
				predV0=_mm256_max_epi16(predV0, amin);
				predV1=_mm256_max_epi16(predV1, amin);
				predV0=_mm256_min_epi16(predV0, amax);
				predV1=_mm256_min_epi16(predV1, amax);
				
				msyms0=_mm256_sub_epi16(myuv[4], predV0);
				msyms1=_mm256_sub_epi16(myuv[5], predV1);
				ecurr[4]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));
				ecurr[5]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms0=_mm256_sub_epi16(msyms0, amin);
				msyms1=_mm256_sub_epi16(msyms1, amin);
				msyms0=_mm256_and_si256(msyms0, bytemask);
				msyms1=_mm256_and_si256(msyms1, bytemask);
				msyms8=_mm256_packus_epi16(msyms0, msyms1);
				mctx=_mm256_packus_epi16(ctxV0, ctxV1);
				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
				mctx=_mm256_permute4x64_epi64(mctx, _MM_SHUFFLE(3, 1, 2, 0));
				ctxV0=_mm256_add_epi16(ctxV0, mctxvoffset);
				ctxV1=_mm256_add_epi16(ctxV1, mctxvoffset);
				ctxV0=_mm256_slli_epi16(ctxV0, 8);
				ctxV1=_mm256_slli_epi16(ctxV1, 8);
				ctxV0=_mm256_or_si256(ctxV0, msyms0);
				ctxV1=_mm256_or_si256(ctxV1, msyms1);
				_mm256_storeu_si256((__m256i*)(imptr+kx+vidx), msyms8);
				_mm256_storeu_si256((__m256i*)(ctxptr+kx+vidx), mctx);
				_mm256_store_si256((__m256i*)syms+0, ctxV0);
				_mm256_store_si256((__m256i*)syms+1, ctxV1);
				msyms0=_mm256_slli_epi16(myuv[4], 2);
				msyms1=_mm256_slli_epi16(myuv[5], 2);
				msyms0=_mm256_sub_epi16(msyms0, moffset0);
				msyms1=_mm256_sub_epi16(msyms1, moffset1);
				W[4]=msyms0;
				W[5]=msyms1;
				_mm256_store_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+0, msyms0);//store V neighbors
				_mm256_store_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+1, msyms1);
				++hists[syms[0x00]];
				++hists[syms[0x01]];
				++hists[syms[0x02]];
				++hists[syms[0x03]];
				++hists[syms[0x04]];
				++hists[syms[0x05]];
				++hists[syms[0x06]];
				++hists[syms[0x07]];
				++hists[syms[0x08]];
				++hists[syms[0x09]];
				++hists[syms[0x0A]];
				++hists[syms[0x0B]];
				++hists[syms[0x0C]];
				++hists[syms[0x0D]];
				++hists[syms[0x0E]];
				++hists[syms[0x0F]];
				++hists[syms[0x10]];
				++hists[syms[0x11]];
				++hists[syms[0x12]];
				++hists[syms[0x13]];
				++hists[syms[0x14]];
				++hists[syms[0x15]];
				++hists[syms[0x16]];
				++hists[syms[0x17]];
				++hists[syms[0x18]];
				++hists[syms[0x19]];
				++hists[syms[0x1A]];
				++hists[syms[0x1B]];
				++hists[syms[0x1C]];
				++hists[syms[0x1D]];
				++hists[syms[0x1E]];
				++hists[syms[0x1F]];
			}
			else
			{
				//decode main

				//decode Y
				dec_yuv(mstate, &ctxY0, &ctxY1, (int*)CDF2syms, ans_permute, &streamptr, streamend, myuv+0*2);

				//yuv = (sym+pred-128)&255
				myuv[0]=_mm256_add_epi16(myuv[0], predY0);
				myuv[1]=_mm256_add_epi16(myuv[1], predY1);
				myuv[0]=_mm256_add_epi16(myuv[0], amin);
				myuv[1]=_mm256_add_epi16(myuv[1], amin);
				myuv[0]=_mm256_and_si256(myuv[0], bytemask);
				myuv[1]=_mm256_and_si256(myuv[1], bytemask);
				msyms0=_mm256_sub_epi16(myuv[0], predY0);//sub pred
				msyms1=_mm256_sub_epi16(myuv[1], predY1);
				ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));
				ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms8=_mm256_packus_epi16(myuv[0], myuv[1]);
				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
				msyms8=_mm256_xor_si256(msyms8, _mm256_set1_epi8(-128));
				_mm256_storeu_si256((__m256i*)(imptr+kx+yidx), msyms8);//store Y bytes
				W[0]=myuv[0];
				W[1]=myuv[1];
				_mm256_store_si256((__m256i*)(rows[0]+(0+0+0*6)*NCODERS)+0, myuv[0]);//store Y neighbors
				_mm256_store_si256((__m256i*)(rows[0]+(0+0+0*6)*NCODERS)+1, myuv[1]);

				//decode U
				moffset0=_mm256_and_si256(myuv[0], uhelpmask);
				moffset1=_mm256_and_si256(myuv[1], uhelpmask);
				predU0=_mm256_add_epi16(predU0, moffset0);
				predU1=_mm256_add_epi16(predU1, moffset1);
				predU0=_mm256_max_epi16(predU0, amin);
				predU1=_mm256_max_epi16(predU1, amin);
				predU0=_mm256_min_epi16(predU0, amax);
				predU1=_mm256_min_epi16(predU1, amax);
				
				dec_yuv(mstate, &ctxU0, &ctxU1, (int*)CDF2syms, ans_permute, &streamptr, streamend, myuv+1*2);
				
				myuv[2]=_mm256_add_epi16(myuv[2], predU0);
				myuv[3]=_mm256_add_epi16(myuv[3], predU1);
				myuv[2]=_mm256_add_epi16(myuv[2], amin);
				myuv[3]=_mm256_add_epi16(myuv[3], amin);
				myuv[2]=_mm256_and_si256(myuv[2], bytemask);
				myuv[3]=_mm256_and_si256(myuv[3], bytemask);
				msyms0=_mm256_sub_epi16(myuv[2], predU0);//sub pred
				msyms1=_mm256_sub_epi16(myuv[3], predU1);
				ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));
				ecurr[3]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms8=_mm256_packus_epi16(myuv[2], myuv[3]);
				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
				msyms8=_mm256_xor_si256(msyms8, _mm256_set1_epi8(-128));
				_mm256_storeu_si256((__m256i*)(imptr+kx+uidx), msyms8);//store U bytes
				myuv[2]=_mm256_sub_epi16(myuv[2], moffset0);//subtract Uoffset from U
				myuv[3]=_mm256_sub_epi16(myuv[3], moffset1);
				W[2]=myuv[2];
				W[3]=myuv[3];
				_mm256_store_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+0, myuv[2]);//store U neighbors
				_mm256_store_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+1, myuv[3]);
				
				//decode V
				moffset0=_mm256_mullo_epi16(vc0, myuv[0]);
				moffset1=_mm256_mullo_epi16(vc0, myuv[1]);
				moffset0=_mm256_add_epi16(moffset0, _mm256_mullo_epi16(vc1, myuv[2]));
				moffset1=_mm256_add_epi16(moffset1, _mm256_mullo_epi16(vc1, myuv[3]));
				predV0=_mm256_add_epi16(predV0, moffset0);
				predV1=_mm256_add_epi16(predV1, moffset1);
				predV0=_mm256_srai_epi16(predV0, 2);
				predV1=_mm256_srai_epi16(predV1, 2);
				predV0=_mm256_max_epi16(predV0, amin);
				predV1=_mm256_max_epi16(predV1, amin);
				predV0=_mm256_min_epi16(predV0, amax);
				predV1=_mm256_min_epi16(predV1, amax);
				
				dec_yuv(mstate, &ctxV0, &ctxV1, (int*)CDF2syms, ans_permute, &streamptr, streamend, myuv+2*2);
				
				myuv[4]=_mm256_add_epi16(myuv[4], predV0);
				myuv[5]=_mm256_add_epi16(myuv[5], predV1);
				myuv[4]=_mm256_add_epi16(myuv[4], amin);
				myuv[5]=_mm256_add_epi16(myuv[5], amin);
				myuv[4]=_mm256_and_si256(myuv[4], bytemask);
				myuv[5]=_mm256_and_si256(myuv[5], bytemask);
				msyms0=_mm256_sub_epi16(myuv[4], predV0);
				msyms1=_mm256_sub_epi16(myuv[5], predV1);
				ecurr[4]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));
				ecurr[5]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms8=_mm256_packus_epi16(myuv[4], myuv[5]);
				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
				msyms8=_mm256_xor_si256(msyms8, _mm256_set1_epi8(-128));
				_mm256_storeu_si256((__m256i*)(imptr+kx+vidx), msyms8);//store V bytes
				myuv[4]=_mm256_slli_epi16(myuv[4], 2);
				myuv[5]=_mm256_slli_epi16(myuv[5], 2);
				myuv[4]=_mm256_sub_epi16(myuv[4], moffset0);//subtract Voffset from V
				myuv[5]=_mm256_sub_epi16(myuv[5], moffset1);
				W[4]=myuv[4];
				W[5]=myuv[5];
				_mm256_store_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+0, myuv[4]);//store V neighbors
				_mm256_store_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+1, myuv[5]);
				
#ifdef ENABLE_GUIDE
				for(int kx2=0;kx2<3*NCODERS;kx2+=3)
					guide_check(image, kx+kx2, ky);
#endif
			}
			if(use_wg4)//update
			{
				for(int kp=0;kp<8;++kp)
				{
					__m256i t[8];

					//e = abs(curr-wgpred[i])<<1
					__m256i me0=_mm256_sub_epi16(myuv[0], wgpreds[kp*6+0]);
					__m256i me1=_mm256_sub_epi16(myuv[1], wgpreds[kp*6+1]);
					__m256i me2=_mm256_sub_epi16(myuv[2], wgpreds[kp*6+2]);
					__m256i me3=_mm256_sub_epi16(myuv[3], wgpreds[kp*6+3]);
					__m256i me4=_mm256_sub_epi16(myuv[4], wgpreds[kp*6+4]);
					__m256i me5=_mm256_sub_epi16(myuv[5], wgpreds[kp*6+5]);
					me0=_mm256_abs_epi16(me0);
					me1=_mm256_abs_epi16(me1);
					me2=_mm256_abs_epi16(me2);
					me3=_mm256_abs_epi16(me3);
					me4=_mm256_abs_epi16(me4);
					me5=_mm256_abs_epi16(me5);
					me0=_mm256_slli_epi16(me0, 1);
					me1=_mm256_slli_epi16(me1, 1);
					me2=_mm256_slli_epi16(me2, 1);
					me3=_mm256_slli_epi16(me3, 1);
					me4=_mm256_slli_epi16(me4, 1);
					me5=_mm256_slli_epi16(me5, 1);

					//eW += (e-eW+(1<<7>>1))>>7
#ifdef WG_DISABLE_eW
					(void)wgWerrors;
#else
					__m256i bias=_mm256_set1_epi16(1<<7>>1);
					t[0]=_mm256_sub_epi16(me0, wgWerrors[kp*6+0]);
					t[1]=_mm256_sub_epi16(me1, wgWerrors[kp*6+1]);
					t[2]=_mm256_sub_epi16(me2, wgWerrors[kp*6+2]);
					t[3]=_mm256_sub_epi16(me3, wgWerrors[kp*6+3]);
					t[4]=_mm256_sub_epi16(me4, wgWerrors[kp*6+4]);
					t[5]=_mm256_sub_epi16(me5, wgWerrors[kp*6+5]);
					t[0]=_mm256_add_epi16(t[0], bias);
					t[1]=_mm256_add_epi16(t[1], bias);
					t[2]=_mm256_add_epi16(t[2], bias);
					t[3]=_mm256_add_epi16(t[3], bias);
					t[4]=_mm256_add_epi16(t[4], bias);
					t[5]=_mm256_add_epi16(t[5], bias);
					t[0]=_mm256_srai_epi16(t[0], 1);
					t[1]=_mm256_srai_epi16(t[1], 1);
					t[2]=_mm256_srai_epi16(t[2], 1);
					t[3]=_mm256_srai_epi16(t[3], 1);
					t[4]=_mm256_srai_epi16(t[4], 1);
					t[5]=_mm256_srai_epi16(t[5], 1);
					wgWerrors[kp*6+0]=_mm256_add_epi16(wgWerrors[kp*6+0], t[0]);
					wgWerrors[kp*6+1]=_mm256_add_epi16(wgWerrors[kp*6+1], t[1]);
					wgWerrors[kp*6+2]=_mm256_add_epi16(wgWerrors[kp*6+2], t[2]);
					wgWerrors[kp*6+3]=_mm256_add_epi16(wgWerrors[kp*6+3], t[3]);
					wgWerrors[kp*6+4]=_mm256_add_epi16(wgWerrors[kp*6+4], t[4]);
					wgWerrors[kp*6+5]=_mm256_add_epi16(wgWerrors[kp*6+5], t[5]);
#endif
					//curr = (2*eW+e+eNE)>>1		(__m256i*)(erows[-Y]+((P+X*NPREDS)*3+C)*NCODERS)+R
					t[0]=_mm256_slli_epi16(_mm256_load_si256((__m256i*)(erows[0]+((kp-1*WG_NPREDS)*3+0)*NCODERS)+0), 1);//eW
					t[1]=_mm256_slli_epi16(_mm256_load_si256((__m256i*)(erows[0]+((kp-1*WG_NPREDS)*3+0)*NCODERS)+1), 1);
					t[2]=_mm256_slli_epi16(_mm256_load_si256((__m256i*)(erows[0]+((kp-1*WG_NPREDS)*3+1)*NCODERS)+0), 1);
					t[3]=_mm256_slli_epi16(_mm256_load_si256((__m256i*)(erows[0]+((kp-1*WG_NPREDS)*3+1)*NCODERS)+1), 1);
					t[4]=_mm256_slli_epi16(_mm256_load_si256((__m256i*)(erows[0]+((kp-1*WG_NPREDS)*3+2)*NCODERS)+0), 1);
					t[5]=_mm256_slli_epi16(_mm256_load_si256((__m256i*)(erows[0]+((kp-1*WG_NPREDS)*3+2)*NCODERS)+1), 1);
					t[0]=_mm256_add_epi16(t[0], me0);
					t[1]=_mm256_add_epi16(t[1], me1);
					t[2]=_mm256_add_epi16(t[2], me2);
					t[3]=_mm256_add_epi16(t[3], me3);
					t[4]=_mm256_add_epi16(t[4], me4);
					t[5]=_mm256_add_epi16(t[5], me5);
					__m256i eNE0=_mm256_load_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+0)*NCODERS)+0);//eNE
					__m256i eNE1=_mm256_load_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+0)*NCODERS)+1);
					__m256i eNE2=_mm256_load_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+1)*NCODERS)+0);
					__m256i eNE3=_mm256_load_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+1)*NCODERS)+1);
					__m256i eNE4=_mm256_load_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+2)*NCODERS)+0);
					__m256i eNE5=_mm256_load_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+2)*NCODERS)+1);
					t[0]=_mm256_add_epi16(t[0], eNE0);
					t[1]=_mm256_add_epi16(t[1], eNE1);
					t[2]=_mm256_add_epi16(t[2], eNE2);
					t[3]=_mm256_add_epi16(t[3], eNE3);
					t[4]=_mm256_add_epi16(t[4], eNE4);
					t[5]=_mm256_add_epi16(t[5], eNE5);
					t[0]=_mm256_srli_epi16(t[0], 2);//eW
					t[1]=_mm256_srli_epi16(t[1], 2);
					t[2]=_mm256_srli_epi16(t[2], 2);
					t[3]=_mm256_srli_epi16(t[3], 2);
					t[4]=_mm256_srli_epi16(t[4], 2);
					t[5]=_mm256_srli_epi16(t[5], 2);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+0*WG_NPREDS)*3+0)*NCODERS)+0, t[0]);//ecurr
					_mm256_store_si256((__m256i*)(erows[1]+((kp+0*WG_NPREDS)*3+0)*NCODERS)+1, t[1]);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+0*WG_NPREDS)*3+1)*NCODERS)+0, t[2]);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+0*WG_NPREDS)*3+1)*NCODERS)+1, t[3]);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+0*WG_NPREDS)*3+2)*NCODERS)+0, t[4]);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+0*WG_NPREDS)*3+2)*NCODERS)+1, t[5]);

					//eNE += e
					eNE0=_mm256_add_epi16(eNE0, me0);
					eNE1=_mm256_add_epi16(eNE1, me1);
					eNE2=_mm256_add_epi16(eNE2, me2);
					eNE3=_mm256_add_epi16(eNE3, me3);
					eNE4=_mm256_add_epi16(eNE4, me4);
					eNE5=_mm256_add_epi16(eNE5, me5);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+0)*NCODERS)+0, eNE0);//eNE
					_mm256_store_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+0)*NCODERS)+1, eNE1);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+1)*NCODERS)+0, eNE2);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+1)*NCODERS)+1, eNE3);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+2)*NCODERS)+0, eNE4);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+2)*NCODERS)+1, eNE5);
				}
				//for(int kr=0;kr<6;++kr)//reuse register but fragmented access
				//{
				//	__m256i me0=_mm256_sub_epi16(myuv[kr], wgpreds[0*6+kr]);
				//	__m256i me1=_mm256_sub_epi16(myuv[kr], wgpreds[1*6+kr]);
				//	__m256i me2=_mm256_sub_epi16(myuv[kr], wgpreds[2*6+kr]);
				//	__m256i me3=_mm256_sub_epi16(myuv[kr], wgpreds[3*6+kr]);
				//	__m256i me4=_mm256_sub_epi16(myuv[kr], wgpreds[4*6+kr]);
				//	__m256i me5=_mm256_sub_epi16(myuv[kr], wgpreds[5*6+kr]);
				//	__m256i me6=_mm256_sub_epi16(myuv[kr], wgpreds[6*6+kr]);
				//	__m256i me7=_mm256_sub_epi16(myuv[kr], wgpreds[7*6+kr]);
				//}
				erows[0]+=NCODERS*3*WG_NPREDS;
				erows[1]+=NCODERS*3*WG_NPREDS;
			}
			eNEEE[0]=_mm256_load_si256((__m256i*)(rows[1]+(3+0+3*6)*NCODERS)+0);
			eNEEE[1]=_mm256_load_si256((__m256i*)(rows[1]+(3+0+3*6)*NCODERS)+1);
			eNEEE[2]=_mm256_load_si256((__m256i*)(rows[1]+(3+1+3*6)*NCODERS)+0);
			eNEEE[3]=_mm256_load_si256((__m256i*)(rows[1]+(3+1+3*6)*NCODERS)+1);
			eNEEE[4]=_mm256_load_si256((__m256i*)(rows[1]+(3+2+3*6)*NCODERS)+0);
			eNEEE[5]=_mm256_load_si256((__m256i*)(rows[1]+(3+2+3*6)*NCODERS)+1);
			eW[0]=_mm256_slli_epi16(eW[0], 1);
			eW[1]=_mm256_slli_epi16(eW[1], 1);
			eW[2]=_mm256_slli_epi16(eW[2], 1);
			eW[3]=_mm256_slli_epi16(eW[3], 1);
			eW[4]=_mm256_slli_epi16(eW[4], 1);
			eW[5]=_mm256_slli_epi16(eW[5], 1);
			eW[0]=_mm256_add_epi16(eW[0], _mm256_max_epi16(eNEE[0], eNEEE[0]));
			eW[1]=_mm256_add_epi16(eW[1], _mm256_max_epi16(eNEE[1], eNEEE[1]));
			eW[2]=_mm256_add_epi16(eW[2], _mm256_max_epi16(eNEE[2], eNEEE[2]));
			eW[3]=_mm256_add_epi16(eW[3], _mm256_max_epi16(eNEE[3], eNEEE[3]));
			eW[4]=_mm256_add_epi16(eW[4], _mm256_max_epi16(eNEE[4], eNEEE[4]));
			eW[5]=_mm256_add_epi16(eW[5], _mm256_max_epi16(eNEE[5], eNEEE[5]));
			ecurr[0]=_mm256_slli_epi16(ecurr[0], GRBITS);
			ecurr[1]=_mm256_slli_epi16(ecurr[1], GRBITS);
			ecurr[2]=_mm256_slli_epi16(ecurr[2], GRBITS);
			ecurr[3]=_mm256_slli_epi16(ecurr[3], GRBITS);
			ecurr[4]=_mm256_slli_epi16(ecurr[4], GRBITS);
			ecurr[5]=_mm256_slli_epi16(ecurr[5], GRBITS);
			eW[0]=_mm256_add_epi16(eW[0], ecurr[0]);
			eW[1]=_mm256_add_epi16(eW[1], ecurr[1]);
			eW[2]=_mm256_add_epi16(eW[2], ecurr[2]);
			eW[3]=_mm256_add_epi16(eW[3], ecurr[3]);
			eW[4]=_mm256_add_epi16(eW[4], ecurr[4]);
			eW[5]=_mm256_add_epi16(eW[5], ecurr[5]);
			eW[0]=_mm256_srli_epi16(eW[0], 2);
			eW[1]=_mm256_srli_epi16(eW[1], 2);
			eW[2]=_mm256_srli_epi16(eW[2], 2);
			eW[3]=_mm256_srli_epi16(eW[3], 2);
			eW[4]=_mm256_srli_epi16(eW[4], 2);
			eW[5]=_mm256_srli_epi16(eW[5], 2);
			_mm256_store_si256((__m256i*)(rows[0]+(3+0+0*6)*NCODERS)+0, eW[0]);//store current contexts
			_mm256_store_si256((__m256i*)(rows[0]+(3+0+0*6)*NCODERS)+1, eW[1]);
			_mm256_store_si256((__m256i*)(rows[0]+(3+1+0*6)*NCODERS)+0, eW[2]);
			_mm256_store_si256((__m256i*)(rows[0]+(3+1+0*6)*NCODERS)+1, eW[3]);
			_mm256_store_si256((__m256i*)(rows[0]+(3+2+0*6)*NCODERS)+0, eW[4]);
			_mm256_store_si256((__m256i*)(rows[0]+(3+2+0*6)*NCODERS)+1, eW[5]);
			eNEE[0]=eNEEE[0];
			eNEE[1]=eNEEE[1];
			eNEE[2]=eNEEE[2];
			eNEE[3]=eNEEE[3];
			eNEE[4]=eNEEE[4];
			eNEE[5]=eNEEE[5];
			NW[0]=N[0];
			NW[1]=N[1];
			NW[2]=N[2];
			NW[3]=N[3];
			NW[4]=N[4];
			NW[5]=N[5];
			rows[0]+=6*NCODERS;
			rows[1]+=6*NCODERS;
			rows[2]+=6*NCODERS;
			rows[3]+=6*NCODERS;
		}
		imptr+=rowstride;
		ctxptr+=rowstride;
	}
	prof_checkpoint("main");
	_mm_free(wgerrors);
	_mm_free(wgstate);
	if(mainwidth<iw)
	{
		int sym=0;
		int vfromy=-(combination[II_COEFF_U_SUB_Y]!=0);
		int vc0=combination[II_COEFF_V_SUB_Y];
		int vc1=combination[II_COEFF_V_SUB_U];
		int offset=0;
		short N[32*3]={0}, yuv[32*3]={0};
		int rembytes=3*remwidth;
		unsigned state=0;
		if(!fwd)
		{
#ifdef _DEBUG
			if(streamptr>streamend)
				LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=*(unsigned*)streamptr;
			streamptr+=4;
		}
		imptr=image+mainbytes;
		yidx=combination[II_PERM_Y],
		uidx=combination[II_PERM_U],
		vidx=combination[II_PERM_V];
		for(int ky=0;ky<ih;++ky)//remainder coding loop		simple differentiation, static-o0 coding
		{
			for(int kx=0;kx<rembytes;kx+=3)
			{
				if(fwd)
				{
					//sym = (yuv-N+128)&255
					yuv[kx+0]=imptr[kx+yidx];
					yuv[kx+1]=imptr[kx+uidx];
					yuv[kx+2]=imptr[kx+vidx];
					imptr[kx+0]=sym=(unsigned char)(yuv[kx+0]-N[kx+0]+128);
					N[kx+0]=yuv[kx+0];
					++rhist[256*0+sym];

					offset=yuv[kx+0]&vfromy;
					N[kx+1]+=offset;
					CLAMP2(N[kx+1], -128, 127);
					imptr[kx+1]=sym=(unsigned char)(yuv[kx+1]-N[kx+1]+128);
					N[kx+1]=yuv[kx+1]-offset;
					++rhist[256*1+sym];

					offset=vc0*yuv[kx+0]+vc1*yuv[kx+1];
					int vpred=(N[kx+2]+offset)>>2;
					CLAMP2(vpred, -128, 127);
					imptr[kx+2]=sym=(unsigned char)(yuv[kx+2]-vpred+128);
					N[kx+2]=4*yuv[kx+2]-offset;
					++rhist[256*2+sym];
				}
				else
				{
					int info;
					info=rCDF2syms[state&((1<<PROBBITS)-1)];
					yuv[kx+0]=(unsigned char)(info+N[kx+0]-128);
					N[kx+0]=yuv[kx+0];
#ifdef ANS_VAL
					ansval_check(&state, sizeof(state));
#endif
					state=(state>>PROBBITS)*(info>>(PROBBITS+8))+(info<<PROBBITS>>(32-PROBBITS));
					if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
					{
#ifdef _DEBUG
						if(streamptr>streamend)
							LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
						state=state<<16|*(unsigned*)streamptr;
						streamptr+=4;
					}

					offset=yuv[kx+0]&vfromy;
					N[kx+1]+=offset;
					CLAMP2(N[kx+1], -128, 127);
					info=rCDF2syms[state&((1<<PROBBITS)-1)];
					yuv[kx+1]=(unsigned char)(info+N[kx+1]-128);
					N[kx+1]=yuv[kx+1]-offset;
#ifdef ANS_VAL
					ansval_check(&state, sizeof(state));
#endif
					state=(state>>PROBBITS)*(info>>(PROBBITS+8))+(info<<PROBBITS>>(32-PROBBITS));
					if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
					{
#ifdef _DEBUG
						if(streamptr>streamend)
							LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
						state=state<<16|*(unsigned*)streamptr;
						streamptr+=4;
					}

					offset=vc0*yuv[kx+0]+vc1*yuv[kx+1];
					N[kx+2]+=offset;
					CLAMP2(N[kx+2], -128, 127);
					info=rCDF2syms[state&((1<<PROBBITS)-1)];
					yuv[kx+2]=(unsigned char)(info+N[kx+2]-128);
					N[kx+2]=yuv[kx+2]-offset;
#ifdef ANS_VAL
					ansval_check(&state, sizeof(state));
#endif
					state=(state>>PROBBITS)*(info>>(PROBBITS+8))+(info<<PROBBITS>>(32-PROBBITS));
					if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
					{
#ifdef _DEBUG
						if(streamptr>streamend)
							LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
						state=state<<16|*(unsigned*)streamptr;
						streamptr+=4;
					}
					imptr[kx+yidx]=(unsigned char)yuv[kx+0];
					imptr[kx+uidx]=(unsigned char)yuv[kx+1];
					imptr[kx+vidx]=(unsigned char)yuv[kx+2];
#ifdef ENABLE_GUIDE
					guide_check(image, kx+mainbytes, ky);
#endif
				}
			}
			imptr+=rowstride;
		}
	}
	prof_checkpoint("remainder");
	if(fwd)//all rANS encoding is bwd-bwd
	{
		//normalize/integrate hists
		rANS_SIMD_SymInfo *syminfo=(rANS_SIMD_SymInfo*)CDF2syms;
		rANS_SIMD_SymInfo *rsyminfo=(rANS_SIMD_SymInfo*)rCDF2syms;
		for(int kc=0;kc<3*NCTX;++kc)
			enc_hist2stats(hists+(ptrdiff_t)256*kc, syminfo+(ptrdiff_t)256*kc, &ctxmask, kc);
		if(blockw*NCODERS<iw)
		{
			enc_hist2stats(rhist+(ptrdiff_t)256*0, rsyminfo+(ptrdiff_t)256*0, &ctxmask, 3*NCTX+0);
			enc_hist2stats(rhist+(ptrdiff_t)256*1, rsyminfo+(ptrdiff_t)256*1, &ctxmask, 3*NCTX+1);
			enc_hist2stats(rhist+(ptrdiff_t)256*2, rsyminfo+(ptrdiff_t)256*2, &ctxmask, 3*NCTX+2);
		}

		//encode remainder
		if(blockw*NCODERS<iw)
		{
			unsigned state=1<<(RANS_STATE_BITS-RANS_RENORM_BITS);
			rANS_SIMD_SymInfo *info=0;
			int rembytes=3*remwidth;
			imptr=image+usize-rembytes;
			for(int ky=ih-1;ky>=0;--ky)
			{
				for(int kx=rembytes-3;kx>=0;kx-=3)
				{
					info=rsyminfo+imptr[kx+2]+256*2;
					if(state>info->smax)
					{
						streamptr-=2;
#ifdef _DEBUG
						if(streamptr<=image)
							LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
						*(unsigned short*)streamptr=(unsigned short)state;
						state>>=RANS_RENORM_BITS;
					}
					state+=((unsigned long long)state*info->invf>>info->sh)*info->negf+info->cdf;
#ifndef STREAM_STATSONLY
#ifdef ANS_VAL
					ansval_push(&state, sizeof(state));
#endif
#endif

					info=rsyminfo+imptr[kx+1]+256*1;
					if(state>info->smax)
					{
						streamptr-=2;
#ifdef _DEBUG
						if(streamptr<=image)
							LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
						*(unsigned short*)streamptr=(unsigned short)state;
						state>>=RANS_RENORM_BITS;
					}
					state+=((unsigned long long)state*info->invf>>info->sh)*info->negf+info->cdf;
#ifndef STREAM_STATSONLY
#ifdef ANS_VAL
					ansval_push(&state, sizeof(state));
#endif
#endif

					info=rsyminfo+imptr[kx+0]+256*0;
					if(state>info->smax)
					{
						streamptr-=2;
#ifdef _DEBUG
						if(streamptr<=image)
							LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
						*(unsigned short*)streamptr=(unsigned short)state;
						state>>=RANS_RENORM_BITS;
					}
					state+=((unsigned long long)state*info->invf>>info->sh)*info->negf+info->cdf;
#ifndef STREAM_STATSONLY
#ifdef ANS_VAL
					ansval_push(&state, sizeof(state));
#endif
#endif
				}
				imptr-=rowstride;
			}
			//flush
			streamptr-=4;
#ifdef _DEBUG
			if(streamptr<=image)
				LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
			*(unsigned*)streamptr=state;
		}
		prof_checkpoint("encode remainder");
		profile_size(streamptr, "/ %9d bytes remainders", 3*remwidth*ih);

		//encode main
		mstate[0]=_mm256_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
		mstate[1]=mstate[0];
		mstate[2]=mstate[0];
		mstate[3]=mstate[0];
		imptr=image+usize-rowstride;
		ctxptr=context+usize-rowstride;
		for(int ky=ih-1;ky>=0;--ky)
		{
			__m256i mkctxinit=_mm256_set1_epi32((NCTX<<8)*2);
			__m256i mkctxstep=_mm256_set1_epi32(NCTX<<8);
			__m256i mkctx=mkctxinit;
#ifdef ESTIMATE_SIZE
			int kc=2;
#endif
			ALIGN(32) int indices[32]={0};
			//indices are permuted:
			//00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 0F	loaded bytes
			//00 04 08 0C 10 14 18 1C 01 05 09 0D 11 15 19 1D 02 06 0A 0E 12 16 1A 1E 03 07 0B 0F 13 17 1B 1F	saved int32s
			//lane 0bEDCBA -> 0xBAEDC	indices[0bBA * 8 + 0bEDC]
			
			//00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 0F	loaded bytes	X
			//00 02 04 06 08 0A 0C 0E 10 12 14 16 18 1A 1C 1E 01 03 05 07 09 0B 0D 0F 11 13 15 16 19 1B 1D 1F	saved int16s	X
			//lane 0bEDCBA -> 0bDCBAE
			for(int kx=mainbytes-NCODERS;kx>=0;kx-=NCODERS)
			{
				//idx = ctx<<8|subpixel		FIXME encoder should use 16-bit buffer
				//32*3 -> 32+2 accesses
#if 1
					__m256i mctxmask=_mm256_set1_epi32(0x0000FF00);
					__m256i lanesyms=_mm256_loadu_si256((__m256i*)(imptr+kx));
					__m256i lanectxs=_mm256_loadu_si256((__m256i*)(ctxptr+kx));
					__m256i midx0=_mm256_slli_epi32(lanesyms, 3*8);
					__m256i midx1=_mm256_slli_epi32(lanesyms, 2*8);
					__m256i midx2=_mm256_slli_epi32(lanesyms, 1*8);
					__m256i midx3=lanesyms;
					__m256i cidx0=_mm256_slli_epi32(lanectxs, 1*8);
					__m256i cidx1=lanectxs;
					__m256i cidx2=_mm256_srli_epi32(lanectxs, 1*8);
					__m256i cidx3=_mm256_srli_epi32(lanectxs, 2*8);
				{
					midx0=_mm256_srli_epi32(midx0, 3*8);
					midx1=_mm256_srli_epi32(midx1, 3*8);
					midx2=_mm256_srli_epi32(midx2, 3*8);
					midx3=_mm256_srli_epi32(midx3, 3*8);
					cidx0=_mm256_and_si256(cidx0, mctxmask);
					cidx1=_mm256_and_si256(cidx1, mctxmask);
					cidx2=_mm256_and_si256(cidx2, mctxmask);
					cidx3=_mm256_and_si256(cidx3, mctxmask);
					midx0=_mm256_or_si256(midx0, cidx0);
					midx1=_mm256_or_si256(midx1, cidx1);
					midx2=_mm256_or_si256(midx2, cidx2);
					midx3=_mm256_or_si256(midx3, cidx3);
					midx0=_mm256_add_epi32(midx0, mkctx);
					midx1=_mm256_add_epi32(midx1, mkctx);
					midx2=_mm256_add_epi32(midx2, mkctx);
					midx3=_mm256_add_epi32(midx3, mkctx);
					midx0=_mm256_slli_epi32(midx0, 4);//multiply by sizeof(__m128i)
					midx1=_mm256_slli_epi32(midx1, 4);
					midx2=_mm256_slli_epi32(midx2, 4);
					midx3=_mm256_slli_epi32(midx3, 4);
					_mm256_store_si256((__m256i*)indices+0, midx0);
					_mm256_store_si256((__m256i*)indices+1, midx1);
					_mm256_store_si256((__m256i*)indices+2, midx2);
					_mm256_store_si256((__m256i*)indices+3, midx3);
				}
#endif
				//2 uloads	8 shift/ORs	1 store		32 adds		32 loads
#if 0
				__m256i lanesyms=_mm256_loadu_si256((__m256i*)(imptr+kx));
				__m256i lanectxs=_mm256_loadu_si256((__m256i*)(ctxptr+kx));
				__m256i midx0=_mm256_slli_epi16(lanesyms, 8);
				__m256i midx1=_mm256_srli_epi16(lanesyms, 8);
				midx0=_mm256_srli_epi16(midx0, 8);
				__m256i cidx1=_mm256_srli_epi16(lanectxs, 8);
				__m256i cidx0=_mm256_slli_epi16(lanectxs, 8);
				cidx1=_mm256_slli_epi16(cidx1, 8);
				midx0=_mm256_or_si256(midx0, cidx0);
				midx1=_mm256_or_si256(midx1, cidx1);
				_mm256_store_si256((__m256i*)indices+0, midx0);
				_mm256_store_si256((__m256i*)indices+1, midx1);
#endif
				//unsigned char *imptrx=imptr+kx, *ctxptrx=ctxptr+kx;

				__m256i mmax[4], minvf[4], mcdf[4], mnegf_sh[4];
				{
					__m256i s0, s1, s2, s3;
					__m256i t0, t1, t2, t3;
#define SHUFFLE_PS(LO, HI, IMM8_HHLL) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(LO), _mm256_castsi256_ps(HI), IMM8_HHLL))
					
				//	//original: indices[wave[2]|reg[2]|hi[1]]	permuted: indices[(reg&1)[1] | hi[1] | wave[2] | reg/2[1]]
					s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|0<<3|0<<1|0])));//_mm_extract_epi16 x32	4/1
					s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|0<<3|0<<1|0])));
					s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|0<<3|0<<1|1])));
					s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|0<<3|0<<1|1])));
					s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|1<<3|0<<1|0])), 1);
					s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|1<<3|0<<1|0])), 1);
					s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|1<<3|0<<1|1])), 1);
					s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|1<<3|0<<1|1])), 1);
					//original: indices[wave[2]|reg[2]|hi[1]]	permuted: indices[(wave&1)[1] | reg[2] | hi[1] | wave/2[1]]	X
				//	s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[0*8+0*2+0]<<8|imptrx[0*8+0*2+0])));//32 shift-ORs
				//	s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[0*8+1*2+0]<<8|imptrx[0*8+1*2+0])));
				//	s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[0*8+2*2+0]<<8|imptrx[0*8+2*2+0])));
				//	s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[0*8+3*2+0]<<8|imptrx[0*8+3*2+0])));
				//	s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[0*8+0*2+1]<<8|imptrx[0*8+0*2+1])), 1);
				//	s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[0*8+1*2+1]<<8|imptrx[0*8+1*2+1])), 1);
				//	s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[0*8+2*2+1]<<8|imptrx[0*8+2*2+1])), 1);
				//	s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[0*8+3*2+1]<<8|imptrx[0*8+3*2+1])), 1);
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));//_MM_TRANSPOSE4_PS
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[0]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[0]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[0]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[0]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));

					s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|0<<3|1<<1|0])));
					s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|0<<3|1<<1|0])));
					s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|0<<3|1<<1|1])));
					s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|0<<3|1<<1|1])));
					s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|1<<3|1<<1|0])), 1);
					s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|1<<3|1<<1|0])), 1);
					s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|1<<3|1<<1|1])), 1);
					s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|1<<3|1<<1|1])), 1);
				//	s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[1*8+0*2+0]<<8|imptrx[1*8+0*2+0])));
				//	s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[1*8+1*2+0]<<8|imptrx[1*8+1*2+0])));
				//	s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[1*8+2*2+0]<<8|imptrx[1*8+2*2+0])));
				//	s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[1*8+3*2+0]<<8|imptrx[1*8+3*2+0])));
				//	s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[1*8+0*2+1]<<8|imptrx[1*8+0*2+1])), 1);
				//	s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[1*8+1*2+1]<<8|imptrx[1*8+1*2+1])), 1);
				//	s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[1*8+2*2+1]<<8|imptrx[1*8+2*2+1])), 1);
				//	s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[1*8+3*2+1]<<8|imptrx[1*8+3*2+1])), 1);
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[1]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[1]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[1]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[1]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));

					s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|0<<3|2<<1|0])));
					s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|0<<3|2<<1|0])));
					s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|0<<3|2<<1|1])));
					s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|0<<3|2<<1|1])));
					s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|1<<3|2<<1|0])), 1);
					s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|1<<3|2<<1|0])), 1);
					s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|1<<3|2<<1|1])), 1);
					s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|1<<3|2<<1|1])), 1);
				//	s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[2*8+0*2+0]<<8|imptrx[2*8+0*2+0])));
				//	s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[2*8+1*2+0]<<8|imptrx[2*8+1*2+0])));
				//	s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[2*8+2*2+0]<<8|imptrx[2*8+2*2+0])));
				//	s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[2*8+3*2+0]<<8|imptrx[2*8+3*2+0])));
				//	s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[2*8+0*2+1]<<8|imptrx[2*8+0*2+1])), 1);
				//	s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[2*8+1*2+1]<<8|imptrx[2*8+1*2+1])), 1);
				//	s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[2*8+2*2+1]<<8|imptrx[2*8+2*2+1])), 1);
				//	s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[2*8+3*2+1]<<8|imptrx[2*8+3*2+1])), 1);
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[2]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[2]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[2]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[2]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));

					s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|0<<3|3<<1|0])));
					s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|0<<3|3<<1|0])));
					s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|0<<3|3<<1|1])));
					s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|0<<3|3<<1|1])));
					s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|1<<3|3<<1|0])), 1);
					s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|1<<3|3<<1|0])), 1);
					s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)((char*)syminfo+indices[0<<4|1<<3|3<<1|1])), 1);
					s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)((char*)syminfo+indices[1<<4|1<<3|3<<1|1])), 1);
				//	s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[3*8+0*2+0]<<8|imptrx[3*8+0*2+0])));
				//	s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[3*8+1*2+0]<<8|imptrx[3*8+1*2+0])));
				//	s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[3*8+2*2+0]<<8|imptrx[3*8+2*2+0])));
				//	s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[3*8+3*2+0]<<8|imptrx[3*8+3*2+0])));
				//	s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[3*8+0*2+1]<<8|imptrx[3*8+0*2+1])), 1);
				//	s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[3*8+1*2+1]<<8|imptrx[3*8+1*2+1])), 1);
				//	s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[3*8+2*2+1]<<8|imptrx[3*8+2*2+1])), 1);
				//	s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)syminfo+((ptrdiff_t)ctxptrx[3*8+3*2+1]<<8|imptrx[3*8+3*2+1])), 1);
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[3]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[3]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[3]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[3]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));
				}
#ifdef ESTIMATE_SIZE
				{
					ALIGN(32) int anegf[NCODERS]={0};
					memcpy(anegf, mnegf_sh, sizeof(anegf));
					const double norm=1./(1<<PROBBITS);
					for(int k=0;k<NCODERS;++k)
					{
						int freq=(1<<PROBBITS)-(anegf[k]&0xFFFF);
						if((unsigned)(freq-1)>=(unsigned)((1<<PROBBITS)-1))
							LOG_ERROR("freq = %d", freq);
						esize[kc*NCODERS+k]-=log2(freq*norm)*0.125;
					}
				}
				--kc;
				if(kc<0)
					kc=2;
#endif

				//enc renorm
				__m256i cond0=_mm256_cmpgt_epi32(mstate[0], mmax[0]);
				__m256i cond1=_mm256_cmpgt_epi32(mstate[1], mmax[1]);
				__m256i cond2=_mm256_cmpgt_epi32(mstate[2], mmax[2]);
				__m256i cond3=_mm256_cmpgt_epi32(mstate[3], mmax[3]);
				int mask0=_mm256_movemask_ps(_mm256_castsi256_ps(cond0));
				int mask1=_mm256_movemask_ps(_mm256_castsi256_ps(cond1));
				int mask2=_mm256_movemask_ps(_mm256_castsi256_ps(cond2));
				int mask3=_mm256_movemask_ps(_mm256_castsi256_ps(cond3));
				__m256i idx0=_mm256_load_si256((__m256i*)ans_permute+mask0);
				__m256i idx1=_mm256_load_si256((__m256i*)ans_permute+mask1);
				__m256i idx2=_mm256_load_si256((__m256i*)ans_permute+mask2);
				__m256i idx3=_mm256_load_si256((__m256i*)ans_permute+mask3);
				__m256i emit0=_mm256_permutevar8x32_epi32(_mm256_and_si256(mstate[0], cond0), idx0);
				__m256i emit1=_mm256_permutevar8x32_epi32(_mm256_and_si256(mstate[1], cond1), idx1);
				__m256i emit2=_mm256_permutevar8x32_epi32(_mm256_and_si256(mstate[2], cond2), idx2);
				__m256i emit3=_mm256_permutevar8x32_epi32(_mm256_and_si256(mstate[3], cond3), idx3);
				emit0=_mm256_and_si256(emit0, _mm256_set1_epi32(0xFFFF));
				emit1=_mm256_and_si256(emit1, _mm256_set1_epi32(0xFFFF));
				emit2=_mm256_and_si256(emit2, _mm256_set1_epi32(0xFFFF));
				emit3=_mm256_and_si256(emit3, _mm256_set1_epi32(0xFFFF));
				emit0=_mm256_packus_epi32(emit0, emit1);
				emit2=_mm256_packus_epi32(emit2, emit3);
				emit0=_mm256_permute4x64_epi64(emit0, _MM_SHUFFLE(3, 1, 2, 0));
				emit2=_mm256_permute4x64_epi64(emit2, _MM_SHUFFLE(3, 1, 2, 0));
				__m128i e3=_mm256_extractf128_si256(emit2, 1);
				__m128i e2=_mm256_castsi256_si128(emit2);
				__m128i e1=_mm256_extractf128_si256(emit0, 1);
				__m128i e0=_mm256_castsi256_si128(emit0);
				mask3=_mm_popcnt_u32(mask3);
				mask2=_mm_popcnt_u32(mask2);
				mask1=_mm_popcnt_u32(mask1);
				mask0=_mm_popcnt_u32(mask0);
#ifdef _DEBUG
				if(streamptr-(2*((ptrdiff_t)mask0+mask1+mask2+mask3)+sizeof(__m128i))<=image)
					LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
				_mm_storeu_si128((__m128i*)streamptr-1, e3); streamptr-=mask3*sizeof(short);
				_mm_storeu_si128((__m128i*)streamptr-1, e2); streamptr-=mask2*sizeof(short);
				_mm_storeu_si128((__m128i*)streamptr-1, e1); streamptr-=mask1*sizeof(short);
				_mm_storeu_si128((__m128i*)streamptr-1, e0); streamptr-=mask0*sizeof(short);
				{
					__m256i state0=_mm256_srli_epi32(mstate[0], 16);
					__m256i state1=_mm256_srli_epi32(mstate[1], 16);
					__m256i state2=_mm256_srli_epi32(mstate[2], 16);
					__m256i state3=_mm256_srli_epi32(mstate[3], 16);
					mstate[0]=_mm256_blendv_epi8(mstate[0], state0, cond0);
					mstate[1]=_mm256_blendv_epi8(mstate[1], state1, cond1);
					mstate[2]=_mm256_blendv_epi8(mstate[2], state2, cond2);
					mstate[3]=_mm256_blendv_epi8(mstate[3], state3, cond3);
				}

				//enc update		state += (state*invf>>sh)*negf+cdf
				{
					__m256i lo0=_mm256_mul_epu32(mstate[0], minvf[0]);//q = mulhi32(state, invf)
					__m256i lo1=_mm256_mul_epu32(mstate[1], minvf[1]);
					__m256i lo2=_mm256_mul_epu32(mstate[2], minvf[2]);
					__m256i lo3=_mm256_mul_epu32(mstate[3], minvf[3]);
					__m256i hi0=_mm256_mul_epu32(_mm256_srli_epi64(mstate[0], 32), _mm256_srli_epi64(minvf[0], 32));
					__m256i hi1=_mm256_mul_epu32(_mm256_srli_epi64(mstate[1], 32), _mm256_srli_epi64(minvf[1], 32));
					__m256i hi2=_mm256_mul_epu32(_mm256_srli_epi64(mstate[2], 32), _mm256_srli_epi64(minvf[2], 32));
					__m256i hi3=_mm256_mul_epu32(_mm256_srli_epi64(mstate[3], 32), _mm256_srli_epi64(minvf[3], 32));
					minvf[0]=_mm256_blend_epi32(_mm256_srli_epi64(lo0, 32), hi0, 0xAA);
					minvf[1]=_mm256_blend_epi32(_mm256_srli_epi64(lo1, 32), hi1, 0xAA);
					minvf[2]=_mm256_blend_epi32(_mm256_srli_epi64(lo2, 32), hi2, 0xAA);
					minvf[3]=_mm256_blend_epi32(_mm256_srli_epi64(lo3, 32), hi3, 0xAA);
				}
				{
					__m256i sh0=_mm256_srli_epi32(mnegf_sh[0], 16);
					__m256i sh1=_mm256_srli_epi32(mnegf_sh[1], 16);
					__m256i sh2=_mm256_srli_epi32(mnegf_sh[2], 16);
					__m256i sh3=_mm256_srli_epi32(mnegf_sh[3], 16);
					_mm256_srlv_epi32(minvf[0], sh0);
					_mm256_srlv_epi32(minvf[1], sh1);
					_mm256_srlv_epi32(minvf[2], sh2);
					_mm256_srlv_epi32(minvf[3], sh3);
				}
				mstate[0]=_mm256_add_epi32(mstate[0], mcdf[0]);
				mstate[1]=_mm256_add_epi32(mstate[1], mcdf[1]);
				mstate[2]=_mm256_add_epi32(mstate[2], mcdf[2]);
				mstate[3]=_mm256_add_epi32(mstate[3], mcdf[3]);
				{
					__m256i lomask=_mm256_set1_epi32(0xFFFF);
					__m256i negf0=_mm256_and_si256(mnegf_sh[0], lomask);
					__m256i negf1=_mm256_and_si256(mnegf_sh[1], lomask);
					__m256i negf2=_mm256_and_si256(mnegf_sh[2], lomask);
					__m256i negf3=_mm256_and_si256(mnegf_sh[3], lomask);
					minvf[0]=_mm256_mullo_epi32(minvf[0], negf0);
					minvf[1]=_mm256_mullo_epi32(minvf[1], negf1);
					minvf[2]=_mm256_mullo_epi32(minvf[2], negf2);
					minvf[3]=_mm256_mullo_epi32(minvf[3], negf3);
				}
				mstate[0]=_mm256_add_epi32(mstate[0], minvf[0]);
				mstate[1]=_mm256_add_epi32(mstate[1], minvf[1]);
				mstate[2]=_mm256_add_epi32(mstate[2], minvf[2]);
				mstate[3]=_mm256_add_epi32(mstate[3], minvf[3]);
#ifndef STREAM_STATSONLY
#ifdef ANS_VAL
				ansval_push(mstate, sizeof(int[VAL_NCODERS]));
#endif
#endif

				mkctx=_mm256_sub_epi32(mkctx, mkctxstep);
				mkctx=_mm256_blendv_epi8(mkctx, mkctxinit, _mm256_srai_epi32(mkctx, 31));
			}
			imptr-=rowstride;
			ctxptr-=rowstride;
		}
		//flush
		streamptr-=sizeof(mstate);
#ifdef _DEBUG
		if(streamptr<=image)
			LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
		memcpy(streamptr, mstate, sizeof(mstate));
		//_mm256_storeu_si256((__m256i*)streamptr+0, mstate[0]);
		//_mm256_storeu_si256((__m256i*)streamptr+1, mstate[1]);
		//_mm256_storeu_si256((__m256i*)streamptr+2, mstate[2]);
		//_mm256_storeu_si256((__m256i*)streamptr+3, mstate[3]);
		prof_checkpoint("encode main");
		profile_size(streamptr, "/ %9td bytes main", (ptrdiff_t)mainbytes*ih);

		//pack hists
		{
#ifdef STREAM_STATSONLY
			streamptr=streamstart;
#endif
			BitPackerLIFO ec;
			bitpacker_enc_init(&ec, image, streamptr);
			enc_packhist(&ec, rhist+256*2, ctxmask, 3*NCTX+2);
			enc_packhist(&ec, rhist+256*1, ctxmask, 3*NCTX+1);
			enc_packhist(&ec, rhist+256*0, ctxmask, 3*NCTX+0);
			for(int kc=3*NCTX-1;kc>=0;--kc)
				enc_packhist(&ec, hists+(ptrdiff_t)256*kc, ctxmask, kc);
			bitpacker_enc_flush(&ec);
			streamptr=ec.dstbwdptr;

			streamptr-=8;
			*(unsigned long long*)streamptr=ctxmask;
		}
		prof_checkpoint("pack histograms");
		profile_size(streamptr, "/ %9d bytes overhead", (3*NCTX+3)*12<<8>>3);

		//save compressed file
		{
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", fdst);
				return 1;
			}
			ptrdiff_t csize2=0;
			csize2+=fwrite("29", 1, 2, fdst);
			csize2+=fwrite(&iw, 1, 4, fdst);
			csize2+=fwrite(&ih, 1, 4, fdst);
			int flags=bestrct<<1|use_wg4;
			csize2+=fwrite(&flags, 1, 1, fdst);
#ifdef _DEBUG
			if(streamptr>streamstart)
				LOG_ERROR("OOB ptr %016zX > %016zX", streamptr, streamstart);
			if(streamptr<image)
				LOG_ERROR("OOB ptr %016zX < %016zX", streamptr, image);
#endif
			csize2+=fwrite(streamptr, 1, streamstart-streamptr, fdst);
			fclose(fdst);
			
#ifdef ESTIMATE_SIZE
			double etotal=0;
			for(int k=0;k<3*NCODERS;++k)
				etotal+=esize[k];
			printf("Total estimate  %12.2lf bytes\n", etotal);
#endif
#ifdef LOUD
			printf("%s\n", srcfn);
			printf("%8td/%8td bytes\n", csize2, usize2);
#endif
			(void)csize2;
		}
		prof_checkpoint("fwrite");
		free(hists);
		free(rhist);
	}
	else
	{
		//deinterleave
		{
			unsigned char *temprow=(unsigned char*)pixels;
			unsigned char *tptrs[NCODERS]={0}, *tptrs0[NCODERS]={0};
			{
				unsigned char *temp=temprow;
				for(int k=0;k<NCODERS;++k)
				{
					tptrs[k]=temp;
					temp+=blockw;
				}
			}
			memcpy(tptrs0, tptrs, sizeof(tptrs0));
			imptr=image;
			for(int ky=0;ky<ih;++ky)//deinterleave coder lanes
			{
				unsigned char *dstptr=imptr;
				memcpy(tptrs, tptrs0, sizeof(tptrs));
				int kx=0;
				for(;kx<=rowstride-3*NCODERS;kx+=3*NCODERS)
				{
					//toy example with 2 coders:
					// p0  p2  p4                                      p1  p3  p5
					//|r00 g00 b00|r01 g01 b01|r02 g02 b02|r03 g03 b03|r04 g04 b04|r05 g05 b05|r06 g06 b06|r07 g07 b07|	dst
					//|r00 r04|g00 g04|b00 b04|r01 r05|g01 g05|b01 b05|r02 r06|g02 g06|b02 b06|r03 r07|g03 g07|b03 b07|	src
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
					for(int k=0;k<NCODERS;++k)
						*tptrs[k]++=*dstptr++;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
					for(int k=0;k<NCODERS;++k)
						*tptrs[k]++=*dstptr++;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
					for(int k=0;k<NCODERS;++k)
						*tptrs[k]++=*dstptr++;
				}
				memcpy(imptr, temprow, rowstride);
				imptr+=rowstride;
			}
		}
		prof_checkpoint("deinterleave");

		//save PPM file
		save_ppm(dstfn, image, iw, ih);
		prof_checkpoint("fwrite");
	}
	_mm_free(pixels);
	_mm_free(CDF2syms);
	_mm_free(rCDF2syms);
	free(image);

#ifdef LOUD
	t=time_sec()-t;
	printf("%c  %12.6lf sec  %12.6lf MB/s\n", 'D'+fwd, t, usize2/(t*1024*1024));
#endif
	(void)och_names;
	(void)rct_names;
	return 0;
}