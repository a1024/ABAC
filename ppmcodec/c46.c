#ifdef _MSC_VER
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif
#elif defined __linux__ && !defined _GNU_SOURCE
#	define _GNU_SOURCE
#	include<stddef.h>//ptrdiff_t
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<math.h>
#include<sys/stat.h>
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#else
#include<time.h>
#endif
#include<immintrin.h>


#ifdef _MSC_VER
	#define LOUD
	#define ESTIMATE_SIZES
	#define ENABLE_GUIDE
#endif


//	#define TRACK_TRAVEL
//	#define ENABLE_ANALYSIS2
//	#define MIXTEST	//X
//	#define MIXLEAK	//X
#ifndef ENABLE_ANALYSIS2
	#define USE_WP
#endif

//	#define MATCH_HISTOGRAMS
//	#define MATCH_HISTOGRAMS_TEST
//	#define USE_W
//	#define USE_CG
	#define USE_L1

	#define USE_AC
//	#define GR_L1
//	#define ANALYSIS_SIMD


#ifdef USE_L1
#define PREDLIST\
	PRED(N+W-NW)\
	PRED(W)\
	PRED(2*N-NN)\
	PRED(NE)\
	PRED(2*W-WW)\
	PRED(3*(N-NN)+NNN)\
	PRED(3*(W-WW)+WWW)\
	PRED(W+NE-N)\
	PRED(W+NEE-NE)\
	PRED(N+NE-NNE)\
	PRED(W+NW-NWW)\
	PRED(N+NW-NNW)\

#endif

#ifdef GR_L1
#define GRESTIMLIST\
	GRESTIM(eNW)\
	GRESTIM(eN)\
	GRESTIM(eNE)\
	GRESTIM(eW)\

#endif

enum
{
#ifdef USE_L1
	SHIFT=18,
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED
#endif
	
#ifdef GR_L1
	GRSHIFT=18,
#define GRESTIM(...) +1
	NGRESTIMS=GRESTIMLIST,
#undef  GRESTIM
#else
	GRBITS=4,
#endif
#ifdef USE_AC
	NCTX=18,
	NLEVELS=256,
#endif

	XPAD=8,
	NCH=3,
	NROWS=4,
	NVAL=3,
#ifdef ENABLE_ANALYSIS2
	A2XSTRIDE=1,
	A2YSTRIDE=1,
#endif
};

//runtime
#if 1
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
#ifdef _MSC_VER
#	define	ALIGN(N) __declspec(align(N))
#	define AWM_INLINE __forceinline static
#else
#	define	ALIGN(N) __attribute__((aligned(N)))
#	define AWM_INLINE __attribute__((always_inline)) inline static
#	ifndef _countof
#		define _countof(A) (sizeof(A)/sizeof(*(A)))
#	endif
#endif
#if defined _M_X64 || defined __x86_64__
#define FLOOR_LOG2(X)\
	(sizeof(X)==8?63-(int32_t)_lzcnt_u64(X):31-_lzcnt_u32((uint32_t)(X)))
#else
AWM_INLINE int floor_log2_64(uint64_t n)
{
	int	logn=-!n;
	int	sh=(n>=1ULL<<32)<<5;	logn+=sh, n>>=sh;
		sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
}
AWM_INLINE int floor_log2_32(uint32_t n)
{
	int	logn=-!n;
	int	sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
}
#define FLOOR_LOG2(X)\
	(sizeof(X)==8?floor_log2_64(X):floor_log2_32((uint32_t)(X)))
#endif
#define CVTFP32_I32(X) _mm_cvt_ss2si(_mm_set_ss(X))
#define CVTTFP32_I32(X) _mm_cvtt_ss2si(_mm_set_ss(X))
#define CVTFP64_I64(X) _mm_cvtsd_si64(_mm_set_sd(X))
#define CVTTFP64_I64(X) _mm_cvttsd_si64(_mm_set_sd(X))
static void crash(const char *file, int line, const char *format, ...)
{
	printf("%s(%d):\n", file, line);
	if(format)
	{
		va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
	}
	printf("\n");
	exit(1);
}
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT,##__VA_ARGS__)
static double time_sec(void)
{
#ifdef _WIN32
	static long long t0=0;
	LARGE_INTEGER li;
	double t;
	QueryPerformanceCounter(&li);
	if(!t0)
		t0=li.QuadPart;
	t=(double)(li.QuadPart-t0);
	QueryPerformanceFrequency(&li);
	t/=(double)li.QuadPart;
	return t;
#else
	struct timespec t;
	clock_gettime(CLOCK_REALTIME, &t);//<time.h>
	return t.tv_sec+t.tv_nsec*1e-9;
#endif
}
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static uint8_t *g_image=0;
static double g_sqe[3]={0};
static void guide_save(uint8_t *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(uint8_t*)malloc(size);
	if(!g_image)
	{
		CRASH("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(uint8_t *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		CRASH("Guide error  XY %d %d", kx, ky);
		printf("");
	}
}
//static void guide_update(uint8_t *image, int kx, int ky)
//{
//	int idx=3*(g_iw*ky+kx), diff;
//	diff=g_image[idx+0]-image[idx+0]; g_sqe[0]+=diff*diff; if(abs(diff)>96)CRASH("");
//	diff=g_image[idx+1]-image[idx+1]; g_sqe[1]+=diff*diff; if(abs(diff)>96)CRASH("");
//	diff=g_image[idx+2]-image[idx+2]; g_sqe[2]+=diff*diff; if(abs(diff)>96)CRASH("");
//}
#else
#define guide_save(...)
#define guide_check(...)
#define guide_update(...)
#endif//ENABLE_GUIDE
#endif


#ifdef MIXTEST
double best_mad[3];
#endif

//cRCT
#if 1
	#define ENABLE_EXTENDED_RCT
#ifndef ENABLE_EXTENDED_RCT
#define OCHLIST\
	OCH(YX00) OCH(Y0X0) OCH(Y00X)\
	OCH(CX40) OCH(C0X4) OCH(C40X)
#endif
#ifdef ENABLE_EXTENDED_RCT
#define OCHLIST\
	OCH(YX00) OCH(Y0X0) OCH(Y00X)\
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
#endif
typedef enum _OCHIndex
{
#define OCH(X) OCH_##X,
	OCHLIST
#undef  OCH
	OCH_COUNT,

	OCH_R=OCH_YX00,
	OCH_G=OCH_Y0X0,
	OCH_B=OCH_Y00X,
	OCH_C4X0=OCH_CX40,
	OCH_C04X=OCH_C0X4,
	OCH_CX04=OCH_C40X,
	OCH_BG=OCH_C04X,
	OCH_BR=OCH_C40X,
	OCH_RG=OCH_CX40,
	OCH_RB=OCH_CX04,
	OCH_GB=OCH_C0X4,
	OCH_GR=OCH_C4X0,
#ifdef ENABLE_EXTENDED_RCT
	OCH_R1=OCH_CX13,
	OCH_G1=OCH_C3X1,
	OCH_B1=OCH_C13X,
	OCH_R2=OCH_CX22,
	OCH_G2=OCH_C2X2,
	OCH_B2=OCH_C22X,
	OCH_R3=OCH_CX31,
	OCH_G3=OCH_C1X3,
	OCH_B3=OCH_C31X,
#endif
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
//	II_COEFF_U_SUB_V2,
//	II_COEFF_V_SUB_U2,

	II_COUNT,
} RCTInfoIdx;
//YUV = RCT * RGB	watch out for permutation in last row
//luma: averaging	chroma: subtraction
//example: _X00_40X_3X1 == [1 0 0; -1 0 1; -3/4 1 -1/4]
#ifndef ENABLE_EXTENDED_RCT
#define RCTLIST\
	RCT(_X00_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_X00_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4)\
	RCT(_X00_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0)\
	RCT(_0X0_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0)\
	RCT(_0X0_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4)\
	RCT(_00X_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4)\
	RCT(_00X_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0)\
	RCT(_0X0_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0)\
	RCT(_0X0_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4)\
	RCT(_0X0_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4)\
	RCT(_00X_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0)\
	RCT(_00X_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4)\
	RCT(_00X_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4)\
	RCT(_X00_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0)\
	RCT(_X00_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4)\
	RCT(_X00_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)
#endif
#ifdef ENABLE_EXTENDED_RCT
#define RCTLIST\
	RCT(_X00_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_X00_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4)\
	RCT(_X00_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0)\
	RCT(_0X0_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0)\
	RCT(_0X0_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4)\
	RCT(_00X_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4)\
	RCT(_00X_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0)\
	RCT(_0X0_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0)\
	RCT(_0X0_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4)\
	RCT(_0X0_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4)\
	RCT(_00X_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0)\
	RCT(_00X_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4)\
	RCT(_00X_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4)\
	RCT(_X00_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0)\
	RCT(_X00_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4)\
	RCT(_X00_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)\
	RCT(_X00_0X0_13X,	OCH_R,		OCH_G,		OCH_B1,		0, 1, 2,	0,  1, 3)\
	RCT(_X00_4X0_13X,	OCH_R,		OCH_GR,		OCH_B1,		0, 1, 2,	4,  1, 3)\
	RCT(_X00_00X_3X1,	OCH_R,		OCH_B,		OCH_G1,		0, 2, 1,	0,  3, 1)\
	RCT(_X00_40X_3X1,	OCH_R,		OCH_BR,		OCH_G1,		0, 2, 1,	4,  3, 1)\
	RCT(_0X0_00X_X13,	OCH_G,		OCH_B,		OCH_R1,		1, 2, 0,	0,  1, 3)\
	RCT(_0X0_04X_X13,	OCH_G,		OCH_BG,		OCH_R1,		1, 2, 0,	4,  1, 3)\
	RCT(_0X0_X40_13X,	OCH_G,		OCH_RG,		OCH_B1,		1, 0, 2,	4,  3, 1)\
	RCT(_00X_X04_3X1,	OCH_B,		OCH_RB,		OCH_G1,		2, 0, 1,	4,  1, 3)\
	RCT(_00X_04X_X13,	OCH_B,		OCH_GB,		OCH_R1,		2, 1, 0,	4,  3, 1)\
	RCT(_X00_0X0_22X,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,  2, 2)\
	RCT(_X00_4X0_22X,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	4,  2, 2)\
	RCT(_X00_00X_2X2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,  2, 2)\
	RCT(_X00_40X_2X2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	4,  2, 2)\
	RCT(_0X0_00X_X22,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,  2, 2)\
	RCT(_0X0_04X_X22,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	4,  2, 2)\
	RCT(_0X0_X40_22X,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	4,  2, 2)\
	RCT(_00X_X04_2X2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	4,  2, 2)\
	RCT(_00X_0X4_X22,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	4,  2, 2)\
	RCT(_X00_0X0_31X,	OCH_R,		OCH_G,		OCH_B3,		0, 1, 2,	0,  3, 1)\
	RCT(_X00_4X0_31X,	OCH_R,		OCH_GR,		OCH_B3,		0, 1, 2,	4,  3, 1)\
	RCT(_X00_00X_1X3,	OCH_R,		OCH_B,		OCH_G3,		0, 2, 1,	0,  1, 3)\
	RCT(_X00_40X_1X3,	OCH_R,		OCH_BR,		OCH_G3,		0, 2, 1,	4,  1, 3)\
	RCT(_0X0_00X_X31,	OCH_G,		OCH_B,		OCH_R3,		1, 2, 0,	0,  3, 1)\
	RCT(_0X0_04X_X31,	OCH_G,		OCH_BG,		OCH_R3,		1, 2, 0,	4,  3, 1)\
	RCT(_0X0_X40_31X,	OCH_G,		OCH_RG,		OCH_B3,		1, 0, 2,	4,  1, 3)\
	RCT(_00X_X04_1X3,	OCH_B,		OCH_RB,		OCH_G3,		2, 0, 1,	4,  3, 1)\
	RCT(_00X_0X4_X31,	OCH_B,		OCH_GB,		OCH_R3,		2, 1, 0,	4,  1, 3)
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
#endif
typedef enum _RCTIndex
{
#define RCT(LABEL, ...) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const uint8_t rct_combinations[RCT_COUNT][II_COUNT]=
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
#endif
static int crct_analysis(uint8_t *image, int iw, int ih)
{
	long long counters[OCH_COUNT]={0};
#ifdef ANALYSIS_SIMD
	__m256i mprev=_mm256_setzero_si256();
#else
	int prev[OCH_COUNT]={0};
#endif
	for(uint8_t *ptr=image, *end=image+(ptrdiff_t)3*iw*ih;ptr<end;ptr+=3)
	{
		int
			r=ptr[0]<<2,
			g=ptr[1]<<2,
			b=ptr[2]<<2,
			rg=r-g,
			gb=g-b,
			br=b-r;
#ifdef ANALYSIS_SIMD
		__m256i mval=_mm256_setr_epi16(
			r, g, b,
			rg, gb, br,
			rg+(gb>>2), rg+(br>>2), br+(rg>>2),
			br+(gb>>2), gb+(br>>2), gb+(rg>>2),
			(rg-br)>>1, (gb-rg)>>1, (br-gb)>>1,
			0
		);
		__m256i mdelta=_mm256_sub_epi16(mval, mprev);
		mprev=mval;
		mdelta=_mm256_abs_epi16(mdelta);
		ALIGN(32) int16_t deltas[16];
		_mm256_store_si256((__m256i*)deltas, mdelta);
		counters[0x0]+=deltas[0x0];
		counters[0x1]+=deltas[0x1];
		counters[0x2]+=deltas[0x2];
		counters[0x3]+=deltas[0x3];
		counters[0x4]+=deltas[0x4];
		counters[0x5]+=deltas[0x5];
		counters[0x6]+=deltas[0x6];
		counters[0x7]+=deltas[0x7];
		counters[0x8]+=deltas[0x8];
		counters[0x9]+=deltas[0x9];
		counters[0xA]+=deltas[0xA];
		counters[0xB]+=deltas[0xB];
		counters[0xC]+=deltas[0xC];
		counters[0xD]+=deltas[0xD];
		counters[0xE]+=deltas[0xE];
#else
		counters[0]+=abs(r -prev[0]);
		counters[1]+=abs(g -prev[1]);
		counters[2]+=abs(b -prev[2]);
		counters[3]+=abs(rg-prev[3]);
		counters[4]+=abs(gb-prev[4]);
		counters[5]+=abs(br-prev[5]);
		prev[0]=r;
		prev[1]=g;
		prev[2]=b;
		prev[3]=rg;
		prev[4]=gb;
		prev[5]=br;
#ifdef ENABLE_EXTENDED_RCT
#define UPDATE(IDXA, IDXB, IDXC, A0, B0, C0)\
	do\
	{\
		int a0=A0, b0=B0, c0=C0;\
		counters[IDXA]+=abs(a0-prev[IDXA]);\
		counters[IDXB]+=abs(b0-prev[IDXB]);\
		counters[IDXC]+=abs(c0-prev[IDXC]);\
		prev[IDXA]=a0;\
		prev[IDXB]=b0;\
		prev[IDXC]=c0;\
	}while(0)
		//r-(3*g+b)/4 = r-g-(b-g)/4
		//g-(3*r+b)/4 = g-r-(b-r)/4
		//b-(3*r+g)/4 = b-r-(g-r)/4
		UPDATE(OCH_CX31, OCH_C3X1, OCH_C31X, rg+(gb>>2), rg+(br>>2), br+(rg>>2));

		//r-(g+3*b)/4 = r-b-(g-b)/4
		//g-(r+3*b)/4 = g-b-(r-b)/4
		//b-(r+3*g)/4 = b-g-(r-g)/4
		UPDATE(OCH_CX13, OCH_C1X3, OCH_C13X, br+(gb>>2), gb+(br>>2), gb+(rg>>2));

		//r-(g+b)/2 = (r-g + r-b)/2
		//g-(r+b)/2 = (g-r + g-b)/2
		//b-(r+g)/2 = (b-r + b-g)/2
		UPDATE(OCH_CX22, OCH_C2X2, OCH_C22X, (rg-br)>>1, (gb-rg)>>1, (br-gb)>>1);
#undef  UPDATE
#endif
#endif
	}
	int bestrct=0;
	long long minerr=0;
	for(int kt=0;kt<RCT_COUNT;++kt)
	{
		const uint8_t *rct=rct_combinations[kt];
		long long currerr=
			+counters[rct[0]]
			+counters[rct[1]]
			+counters[rct[2]]
		;
		if(!kt||minerr>currerr)
		{
			minerr=currerr;
			bestrct=kt;
		}
	}
#ifdef MIXTEST
	best_mad[0]=(double)counters[rct_combinations[bestrct][0]]/((double)iw*ih);
	best_mad[1]=(double)counters[rct_combinations[bestrct][1]]/((double)iw*ih);
	best_mad[2]=(double)counters[rct_combinations[bestrct][2]]/((double)iw*ih);
#endif
	return bestrct;
}

#ifdef ESTIMATE_SIZES
static const char *bsize_labels[]=
{
	"unary",
	"bypass",
};
static int g_kc;
static int64_t bsizes[3][2];
#endif
typedef struct _RiceCoder
{
	uint64_t cache, nbits;
	uint8_t *ptr, *end;
} RiceCoder;
AWM_INLINE void rice_init(RiceCoder *ec, uint8_t *start, uint8_t *end)
{
	ec->cache=0;
	ec->nbits=64;
	ec->ptr=start;
	ec->end=end;
}
AWM_INLINE void rice_flush(RiceCoder *ec)
{
	*(uint64_t*)ec->ptr=ec->cache;
	ec->ptr+=8;
}
AWM_INLINE void rice_enc(RiceCoder *ec, int nbypass, int sym)
{
	//buffer: {c,c,c,b,b,a,a,a, f,f,f,e,e,e,d,c}, cache: MSB gg[hhh]000 LSB	nbits is number of ASSIGNED bits
	//written 64-bit words are byte-reversed because the CPU is little-endian

	int nzeros=sym>>nbypass;
	int bypass=sym&0x7FFFFFFF>>(31-nbypass);
//	int bypass=sym&((1<<nbypass)-1);
#ifdef ESTIMATE_SIZES
	bsizes[g_kc][0]+=nzeros+1;
	bsizes[g_kc][1]+=nbypass;
#endif
	if(nzeros>=ec->nbits)//fill the rest of cache with zeros, and flush
	{
		nzeros-=(int)ec->nbits;
		*(uint64_t*)ec->ptr=ec->cache;
		ec->ptr+=8;
		ec->cache=0;
		while(nzeros>=64)//just flush zeros
		{
			nzeros-=64;
			*(uint64_t*)ec->ptr=0;
			ec->ptr+=8;
		}
		ec->nbits=64;
	}
	//now there is room for zeros:  0 <= nzeros < nbits <= 64
	ec->nbits-=nzeros;//emit remaining zeros to cache

	bypass|=1<<nbypass;//append 1 stop bit
	++nbypass;
	if(nbypass>=ec->nbits)//cache would overflow:  fill, flush, and repeat
	{
		nbypass-=(int)ec->nbits;
		ec->cache|=(uint64_t)bypass>>nbypass;
		bypass&=0x7FFFFFFF>>(31-nbypass);
	//	bypass&=(1<<nbypass)-1;
		*(uint64_t*)ec->ptr=ec->cache;
		ec->ptr+=8;
		ec->cache=0;
		ec->nbits=64;
	}
	//now there is room for bypass:  0 <= nbypass < nbits <= 64
	if(nbypass)
	{
		ec->nbits-=nbypass;//emit remaining bypass to cache
		ec->cache|=(uint64_t)bypass<<ec->nbits;
	}
}
AWM_INLINE int rice_dec(RiceCoder *ec, int nbypass)
{
	//cache: MSB 00[hhh]ijj LSB	nbits is number of CLEARED bits (past codes must be cleared from cache)
	
	int sym;

	sym=-(int)ec->nbits;
	while(!ec->cache)
	{
		sym+=64;
		ec->cache=*(uint64_t*)ec->ptr;
		ec->ptr+=8;
	}
	ec->nbits=_lzcnt_u64(ec->cache);
	sym+=(int)ec->nbits;

	sym<<=nbypass;
	ec->cache&=0x7FFFFFFFFFFFFFFF>>ec->nbits;//remove stop bit
	ec->nbits+=(uint64_t)nbypass+1;
	if(ec->nbits>=64)//nbits = nbits0+nbypass > N
	{
		//example: 000000[11 1]1010010	nbits=6, nbypass=3	6+3-8 = 1
		ec->nbits-=64;
		sym|=(int)(ec->cache<<ec->nbits);
		ec->cache=*(uint64_t*)ec->ptr;
		ec->ptr+=8;
		nbypass=(int)ec->nbits;
	}
	if(nbypass)
	{
		sym|=(int)(ec->cache>>(64-ec->nbits));
		ec->cache&=0xFFFFFFFFFFFFFFFF>>ec->nbits;//nbits=61 -> cache&=7;

		//ec->cache&=(1ULL<<sh)-1;//sh=3 -> cache&=7;
	}
	return sym;
}

enum
{
	SHIFTSTART=12,
	NSHIFTS=15,//last one is WP

	NVAL2=NSHIFTS+3,

	YBLOCKS=1,
};
#ifdef ENABLE_ANALYSIS2
static int32_t testhist[3][YBLOCKS][NSHIFTS][NCTX][NLEVELS];
static int64_t weights[3][YBLOCKS][NSHIFTS][NPREDS];
static void analysis2(const char *srcfn, uint8_t *image, int iw, int ih, int bestrct, uint8_t *ret_sh)
{
	const int fwd=1;
//	int bestrct=crct_analysis(image, iw, ih);
	int yidx=rct_combinations[bestrct][II_PERM_Y];
	int uidx=rct_combinations[bestrct][II_PERM_U];
	int vidx=rct_combinations[bestrct][II_PERM_V];
	int uc0=rct_combinations[bestrct][II_COEFF_U_SUB_Y];
	int vc0=rct_combinations[bestrct][II_COEFF_V_SUB_Y];
	int vc1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
	int psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL2]);
	int16_t *pixels=(int16_t*)malloc(psize);
	int rowstride=3*iw;
#if defined MIXTEST || defined LOUD
	double t=time_sec();
#endif
	if(!pixels)
	{
		CRASH("Alloc error");
		exit(1);
		return;
	}
	memset(testhist, 0, sizeof(testhist));
	memset(pixels, 0, psize);
	memset(weights, 0, sizeof(weights));
	for(int kc=0;kc<3;++kc)
	{
		for(int qy=0;qy<YBLOCKS;++qy)
		{
			for(int ks=0;ks<NSHIFTS-1;++ks)
			{
				for(int j=0;j<NPREDS;++j)
					weights[kc][qy][ks][j]=(1LL<<(ks+SHIFTSTART)>>1)/NPREDS;
			}
		}
	}
	for(int ky=0;ky<ih;ky+=A2YSTRIDE)
	{
		uint8_t *imptr=image+rowstride*ky;
#ifdef USE_L1
		int estim[NPREDS]={0}, j;
#endif
		int yuv[3]={0};
		int error=0, sym=0, curr=0;
		int qy=ky*YBLOCKS/ih;
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL2,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL2,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL2,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL2,
		};
		for(int kx=0;kx<iw;kx+=A2XSTRIDE, imptr+=3*A2XSTRIDE)
		{
			int offset=0;
			if(fwd)
			{
				yuv[0]=imptr[yidx];
				yuv[1]=imptr[uidx];
				yuv[2]=imptr[vidx];
			}
			for(int kc=0;kc<3;++kc)
			{
				int16_t
					NNN	=rows[3][0+0*NCH*NROWS*NVAL2],
					NNW	=rows[2][0-1*NCH*NROWS*NVAL2],
					NN	=rows[2][0+0*NCH*NROWS*NVAL2],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL2],
					NWW	=rows[1][0-2*NCH*NROWS*NVAL2],
					NW	=rows[1][0-1*NCH*NROWS*NVAL2],
					N	=rows[1][0+0*NCH*NROWS*NVAL2],
					NE	=rows[1][0+1*NCH*NROWS*NVAL2],
					NEE	=rows[1][0+2*NCH*NROWS*NVAL2],
					NEEE	=rows[1][0+3*NCH*NROWS*NVAL2],
					WWWW	=rows[0][0-4*NCH*NROWS*NVAL2],
					WWW	=rows[0][0-3*NCH*NROWS*NVAL2],
					WW	=rows[0][0-2*NCH*NROWS*NVAL2],
					W	=rows[0][0-1*NCH*NROWS*NVAL2],
					cN	=rows[1][1+0*NCH*NROWS*NVAL2],
					cW	=rows[0][1-1*NCH*NROWS*NVAL2];
				int64_t p1[NSHIFTS];
				int vmax, vmin;
				
				vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
#define PRED(E) estim[j++]=E;
				j=0;
				PREDLIST
#undef  PRED
				curr=yuv[kc]-offset;
				for(int ks=0;ks<NSHIFTS;++ks)
				{
					uint16_t
					//	eNW	=rows[1][ks+1-1*NCH*NROWS*NVAL2],
					//	eN	=rows[1][ks+1+0*NCH*NROWS*NVAL2],
					//	eNE	=rows[1][ks+1+1*NCH*NROWS*NVAL2],
						eNEE	=rows[1][ks+2+2*NCH*NROWS*NVAL2],
						eNEEE	=rows[1][ks+2+3*NCH*NROWS*NVAL2],
						eW	=rows[0][ks+2-1*NCH*NROWS*NVAL2];
					int ctx=FLOOR_LOG2(eW*eW+1);
					if(ctx>NCTX-1)
						ctx=NCTX-1;
					if(ks==NSHIFTS-1)
					{
						int32_t wsum=0;
						int32_t wp[NPREDS]={0};
#define PRED(...) wsum+=wp[j]=0x100000/((int32_t)weights[kc][qy][ks][j]+1); p1[ks]+=wp[j]*estim[j]; ++j;
						j=0;
						PREDLIST
#undef  PRED
						//if(ky==ih/2&&kx==iw/2)//
						//	printf("");

						p1[ks]/=wsum;
					}
					else
					{
						int sh=SHIFTSTART+ks;
						p1[ks]=1LL<<sh>>1;
#define PRED(...) p1[ks]+=weights[kc][qy][ks][j]*estim[j]; ++j;
						j=0;
						PREDLIST
#undef  PRED
						p1[ks]>>=sh;
					}
					int prob=(int)p1[ks];
					CLAMP2(prob, vmin, vmax);
					prob+=offset;
					CLAMP2(prob, 0, 255);

					error=(int8_t)(yuv[kc]-prob);
					sym=error<<1^error>>31;
					++testhist[kc][qy][ks][ctx][sym];
					
					if(ks==NSHIFTS-1)
					{
#define PRED(...) weights[kc][qy][ks][j]+=(((int64_t)abs(curr-estim[j])<<8)-weights[kc][qy][ks][j])>>3; ++j;
						j=0;
						PREDLIST
#undef  PRED
					}
					else
					{
						int e=(curr>p1[ks])-(curr<p1[ks]);
#define PRED(...) weights[kc][qy][ks][j]+=e*estim[j]; ++j;
						j=0;
						PREDLIST
#undef  PRED
#ifdef MIXLEAK
#define PRED(...) weights[kc][qy][ks][j]-=weights[kc][qy][ks][j]>>12; ++j;
						j=0;
						PREDLIST
#undef  PRED
#endif
					}
					rows[0][ks+2]=(2*eW+(sym<<4)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
				}
				//if(ky==ih/2&&kx==iw/2)//
				//	printf("");
				rows[0][0]=curr;
				rows[0][1]=curr-(int)p1[18-SHIFTSTART];
				offset=(kc?vc0*yuv[0]+vc1*yuv[1]:uc0*yuv[0])>>2;
				rows[0]+=NROWS*NVAL2;
				rows[1]+=NROWS*NVAL2;
				rows[2]+=NROWS*NVAL2;
				rows[3]+=NROWS*NVAL2;
			}
		}
	}
	uint8_t kbest[3*YBLOCKS]={0};
	double csizes[3][YBLOCKS][NSHIFTS]={0};
	for(int kc=0;kc<3;++kc)
	{
		for(int qy=0;qy<YBLOCKS;++qy)
		{
			for(int sh=0;sh<NSHIFTS;++sh)
			{
				for(int ctx=0;ctx<NCTX;++ctx)
				{
					int32_t *currhist=testhist[kc][qy][sh][ctx];
					int sum=0;
					for(int ks=0;ks<NLEVELS;++ks)
						sum+=currhist[ks];
					double invsum=1./sum, e=0;
					for(int ks=0;ks<NLEVELS;++ks)
					{
						int freq=currhist[ks];
						if(freq)
							e-=freq*log2((double)freq*invsum);
					}
					csizes[kc][qy][sh]+=e/8;
				}
				if(!sh||csizes[kc][qy][kbest[YBLOCKS*kc+qy]]>csizes[kc][qy][sh])
					kbest[YBLOCKS*kc+qy]=sh;
			}
		}
	}
	if(ret_sh)
	{
		memcpy(ret_sh, kbest, sizeof(kbest));
		//ret_sh[0]=SHIFTSTART+kbest[0];
		//ret_sh[1]=SHIFTSTART+kbest[1];
		//ret_sh[2]=SHIFTSTART+kbest[2];
	}
#if defined MIXTEST || defined LOUD
	t=time_sec()-t;
	printf("WH %d*%d  %lld bytes  \"%s\"\n"
		, iw, ih
		, (int64_t)3*iw*ih
		, srcfn
	);
	printf("channel usize %lld\n", (int64_t)iw*ih);
	printf("%12.6lf sec  %12.6lf MB/s\n"
		, t, (double)3*iw*ih/(t*1024*1024)
	);
#ifdef MIXTEST
	printf("YUV  mad=%12.4lf %12.4lf %12.4lf\n", best_mad[0], best_mad[1], best_mad[2]);
#endif
	for(int qy=0;qy<YBLOCKS;++qy)
	{
		printf("QY = %2d\n", qy);
		for(int sh=0;sh<NSHIFTS;++sh)
		{
			if(sh==NSHIFTS-1)
				printf("  WP:");
			else
				printf("  %2d:"
					, sh+SHIFTSTART
				);
			printf(" %12.2lf%s %12.2lf%s %12.2lf%s\n"
				, csizes[0][qy][sh], sh==kbest[YBLOCKS*0+qy]?" <-":"   "
				, csizes[1][qy][sh], sh==kbest[YBLOCKS*1+qy]?" <-":"   "
				, csizes[2][qy][sh], sh==kbest[YBLOCKS*2+qy]?" <-":"   "
			);
		}
		printf("\n");
	}
#endif
	free(pixels);
#ifdef MIXTEST
	exit(0);
#endif
}
#endif

#ifdef MATCH_HISTOGRAMS
int32_t matchhist[4][256];
uint16_t matchtable[2][256];
#ifdef MATCH_HISTOGRAMS_TEST
static double calc_csize(int32_t *hist)
{
	int histsum=0;
	for(int ks=0;ks<256;++ks)
		histsum+=hist[ks];
	double invsum=1./histsum, e=0;
	for(int ks=0;ks<256;++ks)
	{
		int freq=hist[ks];
		if(freq)
			e-=freq*log2((double)freq*invsum);
	}
	return e/8;
}
#endif
#endif
#ifdef USE_AC
uint16_t hists[3][NCTX][NLEVELS];
uint16_t hcounts[3][NCTX];
#endif
int c46_codec(int argc, char **argv)
{
	const uint16_t tag='4'|'6'<<8;

	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0;
	uint64_t c=0;
	int fwd=0, iw=0, ih=0, bestrct=0;
	int64_t usize=0, ccap=0, csize=0;
	int psize=0;
	int16_t *pixels=0;
	int wpsize=0;
	uint16_t *wperrors=0;
	uint8_t *image=0, *stream=0, *imptr=0;
	int yidx=0, uidx=0, vidx=0, uc0=0, vc0=0, vc1=0;
#ifdef USE_L1
	int64_t weights[NCH][NPREDS]={0};
	//int32_t weights[NCH][NPREDS]=
	//{
	//	{110945, 45985, 18752, 97238, 24553, 2839, -3204, -45897, 11426},
	//	{135091, 35463, 8137, 110673, 10225, 10071, 8608, -61836, 4727},
	//	{150376, 50763, -37948, 96009, -18017, 26950, 21077, -56707, 32104},
	//};
#endif
#ifdef GR_L1
	int64_t grweights[NCH][NGRESTIMS]={0};
#endif
	uint8_t l1sh[3*YBLOCKS]={0}, prevsh[3]={0};
#ifdef USE_AC
	uint64_t low=0, range=0xFFFFFFFFFFFF, code=0;
	uint8_t *streamptr=0;
#else
	RiceCoder ec;
#endif
#ifdef LOUD
	double t=0;
#endif

	if(argc!=3)
	{
		printf(
			"Usage:  \"%s\"  src  dst\n"
			"Only for 24-bit PPM images\n"
			"Built on %s %s\n"
			, argv[0]
			, __DATE__, __TIME__
		);
		return 1;
	}
#ifdef LOUD
	t=time_sec();
#endif
	srcfn=argv[1];
	dstfn=argv[2];
	
	fsrc=fopen(srcfn, "rb");
	if(!fsrc)
	{
		CRASH("Cannot open \"%s\"", srcfn);
		return 1;
	}
	fread(&c, 1, 2, fsrc);
	fwd=c==('P'|'6'<<8);
	if(!fwd&&c!=tag)
	{
		CRASH("Unsupported file \"%s\"", srcfn);
		return 1;
	}
	if(fwd)//parse header
	{
		c=fgetc(fsrc);
		if(c!='\n')
		{
			CRASH("Unsupported PPM file");
			return 1;
		}
		c=fgetc(fsrc);
		while(c=='#')
		{
			c=fgetc(fsrc);
			while(c!='\n')
				c=fgetc(fsrc);
			c=fgetc(fsrc);
		}
		iw=0;
		while((uint32_t)(c-'0')<10)
		{
			iw=10*iw+(int32_t)c-'0';
			c=fgetc(fsrc);
		}
		while(c<=' ')
			c=fgetc(fsrc);
		ih=0;
		while((uint32_t)(c-'0')<10)
		{
			ih=10*ih+(int32_t)c-'0';
			c=fgetc(fsrc);
		}
		while(c=='#')
		{
			c=fgetc(fsrc);
			while(c!='\n')
				c=fgetc(fsrc);
			c=fgetc(fsrc);
		}
		c|=(int64_t)fgetc(fsrc)<<8*1;
		c|=(int64_t)fgetc(fsrc)<<8*2;
		c|=(int64_t)fgetc(fsrc)<<8*3;
		c|=(int64_t)fgetc(fsrc)<<8*4;
		if(c!=(
			(uint64_t)'\n'<<8*0|
			(uint64_t) '2'<<8*1|
			(uint64_t) '5'<<8*2|
			(uint64_t) '5'<<8*3|
			(uint64_t)'\n'<<8*4
		))
		{
			CRASH("Unsupported PPM file");
			return 1;
		}
		ccap=(int64_t)4*iw*ih;
	}
	else
	{
		iw=0;
		ih=0;
		fread(&iw, 1, 3, fsrc);
		fread(&ih, 1, 3, fsrc);
		fread(&bestrct, 1, 1, fsrc);
#ifdef ENABLE_ANALYSIS2
		fread(l1sh, 1, (size_t)3*YBLOCKS, fsrc);
#endif
		{
			struct stat info={0};

			stat(srcfn, &info);
			ccap=(int64_t)info.st_size-ftell(fsrc);
		}
	}
	if(iw<1||ih<1)
	{
		CRASH("Unsupported source file");
		return 1;
	}
	usize=(int64_t)3*iw*ih;
	image=(uint8_t*)malloc(usize);
	stream=(uint8_t*)malloc(ccap);
	psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!image||!stream||!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	wpsize=(iw+2*XPAD)*(int)sizeof(uint16_t[NCH*NROWS*NPREDS]);
	wperrors=(uint16_t*)malloc(wpsize);
	if(!wperrors)
	{
		CRASH("Alloc error");
		return 1;
	}
	memset(wperrors, 0, wpsize);
	if(fwd)
	{
		fread(image, 1, usize, fsrc);
		guide_save(image, iw, ih);
		bestrct=crct_analysis(image, iw, ih);
#ifdef ENABLE_ANALYSIS2
		analysis2(srcfn, image, iw, ih, bestrct, l1sh);
#endif
#ifdef MATCH_HISTOGRAMS
		yidx=rct_combinations[bestrct][II_PERM_Y];
		uidx=rct_combinations[bestrct][II_PERM_U];
		vidx=rct_combinations[bestrct][II_PERM_V];
		uc0=rct_combinations[bestrct][II_COEFF_U_SUB_Y];
		vc0=rct_combinations[bestrct][II_COEFF_V_SUB_Y];
		vc1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
		memset(matchhist, 0, sizeof(matchhist));
		for(uint8_t *ptr=image, *end=image+usize;ptr<end;ptr+=3)
		{
			int y=ptr[yidx], u=ptr[uidx], v=ptr[vidx];
			int offset1=uc0*y>>2;
			int offset2=(vc0*y+vc1*u)>>2;
			++matchhist[0][offset1];
			++matchhist[1][u];
			++matchhist[2][offset2];
			++matchhist[3][v];
		}
		int sum0=0, sum1=0, sum2=0, sum3=0;
		for(int ks=0;ks<256;++ks)//integrate CDFs
		{
			int freq0=matchhist[0][ks];
			int freq1=matchhist[1][ks];
			int freq2=matchhist[2][ks];
			int freq3=matchhist[3][ks];
			matchhist[0][ks]=sum0;
			matchhist[1][ks]=sum1;
			matchhist[2][ks]=sum2;
			matchhist[3][ks]=sum3;
			sum0+=freq0;
			sum1+=freq1;
			sum2+=freq2;
			sum3+=freq3;
		}
		for(int kc=0;kc<2;++kc)
		{
			int32_t *srccdf=matchhist[2*kc+0], *dstcdf=matchhist[2*kc+1];
			uint16_t *table=matchtable[kc];
			int start=0;
			for(int ks=0;ks<256;++ks)
			{
				int s0=start;
				int v1=srccdf[ks];
				while(dstcdf[start]<v1)
					++start;
			//	table[ks]=ks;
				table[ks]=s0;
			//	table[ks]=(2*ks+s0+start+1)>>2;
			//	table[ks]=start;
			}
		}
#ifdef MATCH_HISTOGRAMS_TEST
		int64_t ctr[4]={0};//{ubefore vbefore uafter vafter}
		int prev[4]={0};
		memset(matchhist, 0, sizeof(matchhist));
		for(uint8_t *ptr=image, *end=image+usize;ptr<end;ptr+=3)
		{
			int y=ptr[yidx], u=ptr[uidx], v=ptr[vidx];

			int offset1=uc0*y>>2;
			int pu=prev[0]+offset1;
			CLAMP2(pu, 0, 255);
			ctr[0]+=abs(u-pu);
			++matchhist[0][(u-pu+128)&255];

			int offset2=(vc0*y+vc1*u)>>2;
			int pv=prev[1]+offset2;
			CLAMP2(pv, 0, 255);
			ctr[1]+=abs(v-pv);
			++matchhist[1][(v-pv+128)&255];

			int offset1b=matchtable[0][offset1];
			int pu2=prev[2]+offset1b;
			CLAMP2(pu2, 0, 255);
			ctr[2]+=abs(u-pu2);
			++matchhist[2][(u-pu2+128)&255];

			int offset2b=matchtable[1][offset2];
			int pv2=prev[3]+offset2b;
			CLAMP2(pv2, 0, 255);
			ctr[3]+=abs(v-pv2);
			++matchhist[3][(v-pv2+128)&255];
			
			prev[0]=u-offset1;
			prev[1]=v-offset2;
			prev[2]=u-offset1b;
			prev[3]=v-offset2b;
		}
		double csizes[4]={0};
		csizes[0]=calc_csize(matchhist[0]);
		csizes[1]=calc_csize(matchhist[1]);
		csizes[2]=calc_csize(matchhist[2]);
		csizes[3]=calc_csize(matchhist[3]);
		printf("U before %12.2lf %12lld  ->  U after  %12.2lf %12lld\n", csizes[0], ctr[0], csizes[2], ctr[2]);
		printf("V before %12.2lf %12lld  ->  V after  %12.2lf %12lld\n", csizes[1], ctr[1], csizes[3], ctr[3]);
		exit(0);
#endif
#if 0
		int64_t res=(int64_t)iw*ih;
		int64_t norm=(0xFFFFLL<<32)/res;
		int64_t c0=0, c1=0, c2=0, c3=0;
		for(int ks=0;ks<256;++ks)//normalize
		{
			int freq0=matchhist[0][ks];
			int freq1=matchhist[1][ks];
			int freq2=matchhist[2][ks];
			int freq3=matchhist[3][ks];
			matchhist[0][ks]=(int32_t)(c0*norm>>32);
			matchhist[1][ks]=(int32_t)(c1*norm>>32);
			matchhist[2][ks]=(int32_t)(c2*norm>>32);
			matchhist[3][ks]=(int32_t)(c3*norm>>32);
			c0+=freq0;
			c1+=freq1;
			c2+=freq2;
			c3+=freq3;
		}
#endif
#endif
	}
	else
	{
		fread(stream, 1, ccap, fsrc);
	}
	fclose(fsrc);
	
#ifndef ENABLE_ANALYSIS2
#ifdef USE_WP
	memset(l1sh, NSHIFTS-1, sizeof(l1sh));//WP
#else
	memset(l1sh, 18-SHIFTSTART, sizeof(l1sh));//L1
#endif
#endif
#ifdef USE_AC
	streamptr=stream;
	if(!fwd)
	{
		code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);//load
		code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);
	}
	memset(hists, 0, sizeof(hists));
	memset(hcounts, 0, sizeof(hcounts));
#else
	rice_init(&ec, stream, stream+ccap);
#endif
	yidx=rct_combinations[bestrct][II_PERM_Y];
	uidx=rct_combinations[bestrct][II_PERM_U];
	vidx=rct_combinations[bestrct][II_PERM_V];
	uc0=rct_combinations[bestrct][II_COEFF_U_SUB_Y];
	vc0=rct_combinations[bestrct][II_COEFF_V_SUB_Y];
	vc1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
#if defined USE_L1 && !defined USE_WP
	for(int kc=0;kc<NCH;++kc)
	{
		if(l1sh[kc]!=SHIFTSTART+NSHIFTS-1)
		{
			for(int kp=0;kp<NPREDS;++kp)
				weights[kc][kp]=(1<<SHIFT)/NPREDS;
		}
	}
#endif
#ifdef RICE_L1
	for(int k=0;k<NCH*2;++k)
		((int32_t*)riceweights)[k]=(1<<SHIFT)/2;
#endif
	memset(pixels, 0, psize);
	imptr=image;
	prevsh[0]=l1sh[YBLOCKS*0+0]+SHIFTSTART;
	prevsh[1]=l1sh[YBLOCKS*1+0]+SHIFTSTART;
	prevsh[2]=l1sh[YBLOCKS*2+0]+SHIFTSTART;
	for(int ky=0;ky<ih;++ky)
	{
#ifdef USE_L1
		int estim[NPREDS]={0}, j;
#endif
#ifdef GR_L1
		int grestims[NGRESTIMS]={0};
#endif
		int qy=ky*YBLOCKS/ih;
		uint8_t sh[]=
		{
			l1sh[YBLOCKS*0+qy]+SHIFTSTART,
			l1sh[YBLOCKS*1+qy]+SHIFTSTART,
			l1sh[YBLOCKS*2+qy]+SHIFTSTART,
		};
		int yuv[3]={0};
		int error=0, sym=0, curr=0;
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		uint16_t *wprows[]=
		{
			wperrors+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NPREDS,
			wperrors+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NPREDS,
			wperrors+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NPREDS,
			wperrors+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NPREDS,
		};
#ifdef ENABLE_ANALYSIS2
		if(qy!=(ky-1)*YBLOCKS/ih)
		{
			for(int kc=0;kc<3;++kc)
			{
				if(prevsh[kc]!=sh[kc])
				{
					if(prevsh[kc]==SHIFTSTART+NSHIFTS-1||sh[kc]==SHIFTSTART+NSHIFTS-1)//switch to L1 or WP
						memset(weights[kc], 0, sizeof(int64_t[NPREDS]));
					else if(prevsh[kc]>sh[kc])//shift has decreased
					{
						int s=prevsh[kc]-sh[kc];
						for(int j=0;j<NPREDS;++j)
							weights[kc][j]>>=s;
					}
					else//shift has increased
					{
						int s=sh[kc]-prevsh[kc];
						for(int j=0;j<NPREDS;++j)
							weights[kc][j]<<=s;
					}
				}
			}
		}
#endif
		for(int kx=0;kx<iw;++kx, imptr+=3)
		{
			int offset=0;
			if(fwd)
			{
				yuv[0]=imptr[yidx];
				yuv[1]=imptr[uidx];
				yuv[2]=imptr[vidx];
			}
			for(int kc=0;kc<3;++kc)
			{
				int
					NNN	=rows[3][0+0*NCH*NROWS*NVAL],
					NNW	=rows[2][0-1*NCH*NROWS*NVAL],
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL],
					NWW	=rows[1][0-2*NCH*NROWS*NVAL],
					NW	=rows[1][0-1*NCH*NROWS*NVAL],
					N	=rows[1][0+0*NCH*NROWS*NVAL],
					NE	=rows[1][0+1*NCH*NROWS*NVAL],
					NEE	=rows[1][0+2*NCH*NROWS*NVAL],
					NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
					WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
					WWW	=rows[0][0-3*NCH*NROWS*NVAL],
					WW	=rows[0][0-2*NCH*NROWS*NVAL],
					W	=rows[0][0-1*NCH*NROWS*NVAL],
					
					cNW	=rows[1][1-1*NCH*NROWS*NVAL],
					cN	=rows[1][1+0*NCH*NROWS*NVAL],
					cNE	=rows[1][1+1*NCH*NROWS*NVAL],
					cW	=rows[0][1-1*NCH*NROWS*NVAL],

					eNEE	=rows[1][2+2*NCH*NROWS*NVAL],
					eNEEE	=rows[1][2+3*NCH*NROWS*NVAL],
					eW	=rows[0][2-1*NCH*NROWS*NVAL];
#ifdef GR_L1
				int64_t grestim=1<<GRSHIFT>>1;
#define GRESTIM(E) grestims[j]=E; grestim+=grweights[kc][j]*grestims[j]; ++j;
				j=0;
				GRESTIMLIST
#undef  GRESTIM
				grestim>>=GRSHIFT;
				//int riceestim=(int)((
				//	+(int64_t)grweights[kc][0]*eW
				//	+(int64_t)grweights[kc][1]*eNE
				//	+(int64_t)grweights[kc][2]*eN
				//)>>SHIFT);
				int nbypass=FLOOR_LOG2((int)(grestim<0?0:grestim)+1);
#else
#ifdef USE_AC
				int ctx=FLOOR_LOG2(eW*eW+1);
				if(ctx>NCTX-1)
					ctx=NCTX-1;
#else
				int nbypass=FLOOR_LOG2(((eW+eNE)>>(GRBITS+1))+1);
			//	int nbypass=FLOOR_LOG2((eW>>GRBITS)+1);
#endif
#endif
#ifdef USE_W
				int pred=W;
#endif
#ifdef USE_CG
				int pred=N+W-NW;
				int vmax=N, vmin=W;
				
				if(N<W)vmin=N, vmax=W;
				CLAMP2(pred, vmin, vmax);
#endif
#ifdef USE_L1
				int64_t p1;
				int pred, e;
				int vmax=N, vmin=W;

				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				if(sh[kc]==SHIFTSTART+NSHIFTS-1)
				{
					int32_t wp[NPREDS]={0};
					//e = 2*(N+W+NW)+NN+NNE+NE+WW+I/4
#define PRED(E)\
	wp[j]=\
		+wprows[0][j-1*NCH*NROWS*NPREDS]*2\
		+wprows[1][j+0*NCH*NROWS*NPREDS]*2\
		+wprows[1][j-1*NCH*NROWS*NPREDS]*2\
		+wprows[2][j+0*NCH*NROWS*NPREDS]\
		+wprows[2][j+1*NCH*NROWS*NPREDS]\
		+wprows[1][j+1*NCH*NROWS*NPREDS]\
		+wprows[0][j-1*NCH*NROWS*NPREDS]\
		+(int32_t)(weights[kc][j]>>2)\
	;++j;
					j=0;
					PREDLIST
#undef  PRED
					//if(ky==ih/2&&kx==iw/2)//
					//	printf("");
#if 0
					float coeff=0, wsum=0, psum=0;
#define PRED(E) estim[j]=E; coeff=65536.f/((int32_t)weights[kc][j]+1.0f); wsum+=coeff; psum+=coeff*estim[j]; ++j;
					j=0;
					PREDLIST
#undef  PRED
					//if(wsum)
					//	p1=CVTFP32_I32(psum/wsum);
					//else
					//	p1=estim[0];
					p1=CVTFP32_I32(psum/(wsum+0.5f));
#endif
#if 1
					int32_t wsum=0;
					int32_t coeff=0;
					p1=0;
#define PRED(E) estim[j]=E; coeff=0x100000/(wp[j]+1); wsum+=coeff; p1+=coeff*estim[j]; ++j;
					j=0;
					PREDLIST
#undef  PRED
					int64_t sign=p1>>63;
					p1^=sign;
					p1-=sign;
					p1=(p1+(wsum>>1))/wsum;//towards nearest
					p1^=sign;
					p1-=sign;

				//	if(p1<0)
				//		p1=-((-p1+(wsum>>1))/wsum);//towards nearest
				//	else
				//		p1=(p1+(wsum>>1))/wsum;

				//	p1/=wsum;//towards zero		X
#endif
				}
				else
				{
					p1=1LL<<sh[kc]>>1;
#define PRED(E) estim[j]=E; p1+=weights[kc][j]*estim[j]; ++j;
					j=0;
					PREDLIST
#undef  PRED
					p1>>=sh[kc];
				}
				pred=(int)p1;
				CLAMP2(pred, vmin, vmax);
#endif
				pred+=offset;
				CLAMP2(pred, 0, 255);
#ifdef ESTIMATE_SIZES
				g_kc=kc;
#endif
#ifdef USE_AC
				int den=hcounts[kc][ctx]+NLEVELS, cdf=0, freq;
				uint16_t *currhist=hists[kc][ctx];
				if(fwd)
				{
					error=(int8_t)(yuv[kc]-pred);
					sym=error<<1^error>>31;
					if(range<=0xFFFF)
					{
						*(uint32_t*)streamptr=(uint32_t)(low>>32);
						streamptr+=sizeof(uint32_t);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					for(int t=0;;++t)
					{
						freq=currhist[t]+1;
						if(t>=sym)
							break;
						cdf+=freq;
					}
					low+=range*cdf/den;
					range=range*freq/den-1;
				}
				else
				{
					if(range<=0xFFFF)
					{
						code=code<<32|*(uint32_t*)streamptr;
						streamptr+=sizeof(uint32_t);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					int c=(int)(((code-low+1)*den-1)/range);
					for(sym=0;;++sym)
					{
						freq=currhist[sym]+1;
						if(cdf+freq>c)
							break;
						cdf+=freq;
					}
					low+=range*cdf/den;
					range=range*freq/den-1;
					error=sym>>1^-(sym&1);
					yuv[kc]=(uint8_t)(error+pred);
				}
				++currhist[sym];
				++hcounts[kc][ctx];
				if(hcounts[kc][ctx]>=0xFFFF-2*NLEVELS)
				{
					den=0;
					for(int ks=0;ks<NLEVELS;++ks)
						den+=currhist[ks]>>=1;
					hcounts[kc][ctx]=den;
				}
#else
				if(fwd)
				{
					error=(int8_t)(yuv[kc]-pred);
					sym=error<<1^error>>31;
					rice_enc(&ec, nbypass, sym);
				}
				else
				{
					sym=rice_dec(&ec, nbypass);
					error=sym>>1^-(sym&1);
					yuv[kc]=(uint8_t)(error+pred);
				}
#endif
				curr=yuv[kc]-offset;
				//if(ky==ih/2&&kx==iw/2)//
				//	printf("");
#ifdef USE_L1
				if(sh[kc]==SHIFTSTART+NSHIFTS-1)
				{
					int best=0x7FFFFFFF;
					int32_t errors[NPREDS];
#define PRED(...) errors[j]=abs(curr-estim[j]); if(best>errors[j])best=errors[j]; ++j;
					j=0;
					PREDLIST;
#undef  PRED
	#define PRED(...) wprows[0][j]=abs(curr-estim[j])-best; weights[kc][j]+=(((int64_t)wprows[0][j]<<5)-weights[kc][j]+(1<<3>>1))>>3; ++j;
//	#define PRED(...) weights[kc][j]+=(((int64_t)abs(curr-estim[j])<<8)-weights[kc][j])>>3; ++j;
					j=0;
					PREDLIST
#undef  PRED
				}
				else
				{
					e=(curr>p1)-(curr<p1);
				//	e+=(curr-p1)>>0;
				//	e=curr-(int)p1;
				//	CLAMP2(e, -2, 2);
#define PRED(...) weights[kc][j]+=e*estim[j]; ++j;
					j=0;
					PREDLIST
#undef  PRED
#ifdef MIXLEAK
#define PRED(...) weights[kc][j]-=weights[kc][j]>>12; ++j;
					j=0;
					PREDLIST
#undef  PRED
#endif
				}
#endif
				rows[0][0]=curr;
				rows[0][1]=curr-pred;
#ifdef GR_L1
				int gre=(sym>grestim)-(sym<grestim);
#define GRESTIM(...) grweights[kc][j]+=gre*grestims[j]; ++j;
				j=0;
				GRESTIMLIST
#undef  GRESTIM
				rows[0][1]=sym;
#else
#ifdef USE_AC
				rows[0][2]=(2*eW+(sym<<GRBITS)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
#else
				rows[0][1]=(2*eW+(sym<<GRBITS)+eNEEE)>>2;
#endif
#endif
				offset=(kc?vc0*yuv[0]+vc1*yuv[1]:uc0*yuv[0])>>2;
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				(void)NNN	;
				(void)NN	;
				(void)NNE	;
				(void)NW	;
				(void)N		;
				(void)NE	;
				(void)NEEE	;
				(void)WWW	;
				(void)WW	;
				(void)W		;
				(void)eNEEE	;
				(void)eW	;
			}
			if(!fwd)
			{
				imptr[yidx]=yuv[0];
				imptr[uidx]=yuv[1];
				imptr[vidx]=yuv[2];
				guide_check(image, kx, ky);
			}
		}
	}
	free(pixels);
	free(wperrors);
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing", dstfn);
			return 1;
		}
		if(fwd)
		{
#ifdef USE_AC
			*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;//flush
			*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;
#else
			rice_flush(&ec);
			uint8_t *streamptr=ec.ptr;
#endif

			csize=0;
			csize+=fwrite(&tag, 1, 2, fdst);
			csize+=fwrite(&iw, 1, 3, fdst);
			csize+=fwrite(&ih, 1, 3, fdst);
			csize+=fwrite(&bestrct, 1, 1, fdst);
#ifdef ENABLE_ANALYSIS2
			csize+=fwrite(l1sh, 1, (size_t)3*YBLOCKS, fdst);
#endif
			csize+=fwrite(stream, 1, streamptr-stream, fdst);
		}
		else
		{
			int headersize=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(image, 1, usize, fdst);
			usize+=headersize;
		}
		fclose(fdst);
	}
	free(image);
	free(stream);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
#if defined ESTIMATE_SIZES && !defined USE_AC
		int64_t btotal=0;
		for(int k=0;k<2;++k)
			btotal+=bsizes[0][k]+bsizes[1][k]+bsizes[2][k];
		for(int k=0;k<2;++k)
			printf("%12.2lf %12.2lf %12.2lf    %9.5lf%% %9.5lf%% %9.5lf%%  %s\n"
				, (double)bsizes[0][k]/8.
				, (double)bsizes[1][k]/8.
				, (double)bsizes[2][k]/8.
				, 100.*bsizes[0][k]/btotal
				, 100.*bsizes[1][k]/btotal
				, 100.*bsizes[2][k]/btotal
				, bsize_labels[k]
			);
#endif
		printf("CWH=3*%d*%d  RCT %2d %s  \"%s\"\n", iw, ih, bestrct, rct_names[bestrct], srcfn);
		printf("%10td->%10td  %8.4lf%%  %12.6lf:1  BPD %12.6lf\n"
			, usize
			, csize
			, 100.*csize/usize
			, (double)usize/csize
			, 8.*csize/usize
		);
	}
	printf("%c  %12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB\n"
		, 'D'+fwd
		, t
		, usize/(t*1024*1024)
		, t*1024*1024*1000/usize
	);
#endif
	(void)csize;
	(void)&time_sec;
	(void)och_names;
	(void)rct_names;
	return 0;
}
