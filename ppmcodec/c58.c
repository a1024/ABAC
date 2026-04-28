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
#define _USE_MATH_DEFINES
#include<math.h>
#include<sys/stat.h>
#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>//QueryPerformanceCounter
#elif defined __GNUC__
#include<time.h>
#endif
#include<immintrin.h>//_lzcnt_u32
#ifdef PROFILER
#include"util.h"
#endif


#if defined _MSC_VER && !defined RELEASE
	#define LOUD
//	#define PRINT_RCT
//	#define ESTIMATE_BITSIZE
//	#define PRINT_STATS

	#define ENABLE_GUIDE
	#define FIFOVAL
#endif


	#define USE_CTXCTR
//	#define USE_HIST


#if 1
#define PREDLIST\
	PRED(N)\
	PRED(NNN)\
	PRED(NNNN)\
	PRED(N+N-NN)\
	PRED(N+NE-NNE)\
	PRED(N+NW-NNW)\
	PRED(N+W-NW)\
	PRED(NEEEE)\
	PRED(W)\
	PRED(WWW)\
	PRED(WWWW)\
	PRED(W+W-WW)\
	PRED(W+NE-N)\
	PRED(W+NEE-NE)\
	PRED(W+NW-NWW)\
	PRED(3*(W-WW)+WWW)\

#endif
#if 0
#define PREDLIST\
	PRED(N)\
	PRED(W)\
	PRED(NNN)\
	PRED(WWW)\
	PRED(NEEE)\
	PRED(NEEEE)\
	PRED(2*N-NN)\
	PRED(2*W-WW)\
	PRED(N+W-NW)\
	PRED(W+NE-N)\
	PRED(W+NEE-NE)\
	PRED(N+NE-NNE)\
	PRED(N+NW-NNW)\
	PRED(W+NW-NWW)\

#endif
#if 0
	PRED(N+2*yN-yNN)\
	PRED(NEEE)\
	PRED(4*(N+NNN)-6*NN-NNNN)\
	PRED(4*(W+WWW)-6*WW-WWWW)\
	PRED(W+NEEE-NEE)\

#endif
enum
{
	L1SH=20,
#define PRED(EXPR) +1
	L1NPREDS=PREDLIST,
#undef  PRED

	GRBITS=6,
	NCTX=24,
	GRLIMIT=16,

	PROBBITS_STORE=20,
	PROBBITS_USE=12,
	PROBSHIFT=PROBBITS_STORE-PROBBITS_USE,

	XPAD=8,
	NCH=3,
	NROWS=4,
	NVAL=2,
};




//runtime
#if 1
#ifdef _MSC_VER
#define INLINE __forceinline static
#else
#define INLINE __attribute__((always_inline)) inline static
#endif
#ifndef ALIGN
#ifdef _MSC_VER
#define	ALIGN(N) __declspec(align(N))
#else
#define	ALIGN(N) __attribute__((aligned(N)))
#endif
#endif
#define CLAMP2(X, LO, HI) X=X>LO?X:LO, X=X<HI?X:HI
#define CVTFP32_I32(X) _mm_cvt_ss2si(_mm_set_ss(X))
#define CVTTFP32_I32(X) _mm_cvtt_ss2si(_mm_set_ss(X))
#define CVTFP64_I64(X) _mm_cvtsd_si64(_mm_set_sd(X))
#define CVTTFP64_I64(X) _mm_cvttsd_si64(_mm_set_sd(X))
static void memfill_s(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
{
	size_t copied;
	char *d=(char*)dst;
	const char *s=(const char*)src;
#ifdef _DEBUG
	if(!dstbytes||!srcbytes)
		return;
#endif
	if(dstbytes<srcbytes)
	{
		memcpy(dst, src, dstbytes);
		return;
	}
	copied=srcbytes;
	memcpy(d, s, srcbytes);
	while((copied<<1)<=dstbytes)
	{
		memcpy(d+copied, d, copied);
		copied<<=1;
	}
	if(copied<dstbytes)
		memcpy(d+copied, d, dstbytes-copied);
}
#define FILLMEM_S(PTR, DATA, ASIZE, ESIZE)\
	do\
	{\
		*(PTR)=(DATA);\
		memfill_s((PTR)+1, PTR, (ASIZE)-(ESIZE), ESIZE);\
	}while(0)
static double time_sec2(void)
{
#ifdef _MSC_VER
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
static void crash(const char *file, int line, const char *format, ...)
{
	printf("%s(%d):\n", file, line);
	if(format)
	{
		va_list args;
		printf("\n");
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
		printf("\n");
	}
	printf("CRASH\n");
	exit(1);
}
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT, ##__VA_ARGS__)
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
		CRASH("guide  XY %d %d", kx, ky);
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
#endif
#ifdef FIFOVAL
static ptrdiff_t fifoidx=0, fifocap=0, fifoidx2=0;
static uint32_t *fifoval=0;
static void fifoval_enqueue(uint32_t val)
{
	if(fifoidx+1>=fifocap)
	{
		void *p=0;

		if(!fifocap)
			fifocap=1;
		fifocap<<=1;
		p=realloc(fifoval, fifocap*sizeof(uint32_t));
		if(!p)
		{
			CRASH("Alloc error");
			return;
		}
		fifoval=(uint32_t*)p;
	}
	fifoval[fifoidx++]=val;
}
static void fifoval_check(uint32_t val)
{
	uint32_t val0=fifoval[fifoidx2++];
	if(val!=val0)
	{
		--fifoidx2;
		printf(
			"\n"
			"FIFO Error  at %10lld,  remaining %10lld\n"
			"    0x%08X  !=  original 0x%08X\n"
			"\n"
			, fifoidx2
			, fifoidx-fifoidx2//current element was not decoded successfully
			, val, val0
		);
		for(int k=-32;k<32;++k)
		{
			ptrdiff_t idx=fifoidx2+k;
			if((size_t)idx<(size_t)fifoidx)
			{
				printf(
					"%10td  0x%08X"
					, idx
					, fifoval[idx]
				);
				if(idx<fifoidx2)
					printf("  OK");
				if(idx==fifoidx2)
					printf("  !=  corrupt 0x%08X", val);
				printf("\n");
			}
		}
		CRASH("");
	}
}
#endif
#endif


//cRCT
#if 1
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
//	II_COEFF_U_SUB_V2,
//	II_COEFF_V_SUB_U2,

	II_COUNT,
} RCTInfoIdx;
//YUV = RCT * RGB	example: _400_40X_3X1 == [1 0 0; -1 0 1; -3/4 1 -1/4]
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

#ifdef USE_HIST
#if 1
static uint32_t squash(int32_t d)
{
	static const int t[33]=
	{
		1,2,3,6,10,16,27,45,73,120,194,310,488,747,1101,
		1546,2047,2549,2994,3348,3607,3785,3901,3975,4022,
		4050,4068,4079,4085,4089,4092,4093,4094
	};
	if(d>2047)return 4095;
	if(d<-2047)return 1;
	int w=d&127;
	d=(d>>7)+16;
	return (t[d]*(128-w)+t[(d+1)]*w+64)>>7;
}
#endif
static int32_t stats[3][256][NCTX][257];
#else
static int32_t stats1[3][1024][NCTX][GRLIMIT];//unary
static int32_t stats2[3][1024][256];//remainder
static int32_t stats3[3][8][256];//bypass on GRLIMIT
static const size_t memusage=sizeof(stats1)+sizeof(stats2)+sizeof(stats3);
#endif
typedef struct _ACState
{
	uint64_t low, range, code;
	uint8_t *ptr, *end;

	uint64_t bitidx, totalbits;
	uint64_t n[2];
} ACState;
#ifdef ESTIMATE_BITSIZE
static uint32_t unary_count=0, binary_count=0;
static double shannontable[1<<PROBBITS_USE];
static double bitsizes[3][GRLIMIT+8];
static uint32_t bitctr[3][GRLIMIT+8][2];
static uint32_t winctr[3][GRLIMIT+8];
static int ekc, eidx;
#endif
#ifdef USE_HIST
INLINE void codebit(ACState *ac, int *bit, int p1, const int fwd)
{
	uint64_t r2, mid;
	int rbit;
	//int p1;

#ifdef _MSC_VER
	//if(z1>=z2||t1>=t2)
	//	CRASH("%d %d  %d %d", z1, z2, t1, t2);
	++ac->bitidx;
#endif
	//p1=((CDF[z2]-CDF[z1]+z2-z1)<<PROBBITS_USE)/(CDF[t2]-CDF[t1]+t2-t1);
	if(ac->range<=0xFFFF)
	{
		if(ac->ptr>=ac->end)
		{
#ifdef _MSC_VER
			CRASH("inflation  %d/%d  %8.4lf%%\n"
				, (int32_t)ac->bitidx
				, (int32_t)ac->totalbits
				, 100.*ac->totalbits/ac->bitidx
			);
#endif
			exit(1);
		}
		if(fwd)
			*(uint32_t*)ac->ptr=(uint32_t)(ac->low>>32);
		else
			ac->code=ac->code<<32|*(uint32_t*)ac->ptr;
		ac->ptr+=4;
		ac->low<<=32;
		ac->range=ac->range<<32|0xFFFFFFFF;
		if(ac->range>~ac->low)
			ac->range=~ac->low;
	}
	r2=ac->range*p1>>PROBBITS_USE;
//	r2=(ac->range>>PROBBITS_USE)*p1;
	mid=ac->low+r2;
	ac->range-=r2;
	--r2;
	rbit=*bit;
	rbit=fwd?rbit:ac->code<mid;
	*bit=rbit;
#ifdef FIFOVAL
	if(p1<1||p1>(1<<PROBBITS_USE)-1)
		CRASH("");
	if(fwd)
		fifoval_enqueue(rbit<<PROBBITS_USE^p1);
	else
		fifoval_check(rbit<<PROBBITS_USE^p1);
#endif
	ac->range=rbit?r2:ac->range;
	ac->low=rbit?ac->low:mid;
#ifdef ESTIMATE_BITSIZE
	bitsizes[ekc][eidx]+=shannontable[rbit?(1<<PROBBITS_USE)-p1:p1];
	++bitctr[ekc][eidx][rbit];
	winctr[ekc][eidx]+=rbit==(p1>=1<<PROBBITS_USE);
#endif
}
#else
INLINE void codebit(ACState *ac, int32_t *pcell, int *bit, const int fwd)
{
	uint64_t r2, mid;
	int32_t cell=*pcell;
//	uint32_t p1=cell>>PROBSHIFT;
	uint32_t p1=cell>>(PROBSHIFT+11);
	int rbit;
	
	p1+=p1<1<<PROBBITS_USE>>1;
	//p1+=1<<PROBBITS_USE>>1;
	//CLAMP2(p1, 1, (1<<PROBBITS_USE)-1);

#ifdef _MSC_VER
	++ac->bitidx;
#endif
	if(ac->range<=0xFFFF)
	{
		if(ac->ptr>=ac->end)
		{
#ifdef _MSC_VER
			CRASH("inflation  %d/%d  %8.4lf%%\n"
				, (int32_t)ac->bitidx
				, (int32_t)ac->totalbits
				, 100.*ac->totalbits/ac->bitidx
			);
#endif
			exit(1);
		}
		if(fwd)
			*(uint32_t*)ac->ptr=(uint32_t)(ac->low>>32);
		else
			ac->code=ac->code<<32|*(uint32_t*)ac->ptr;
		ac->ptr+=4;
		ac->low<<=32;
		ac->range=ac->range<<32|0xFFFFFFFF;
		if(ac->range>~ac->low)
			ac->range=~ac->low;
	}
	r2=ac->range*p1>>PROBBITS_USE;
//	r2=(ac->range>>PROBBITS_USE)*p1;
	mid=ac->low+r2;
	ac->range-=r2;
	--r2;
	rbit=*bit;
	rbit=fwd?rbit:ac->code<mid;
	*bit=rbit;
#ifdef FIFOVAL
	//if(p1<1||p1>(1<<PROBBITS_USE)-1)
	//	CRASH("");
	if(fwd)
		fifoval_enqueue(rbit<<PROBBITS_STORE^p1);
	else
		fifoval_check(rbit<<PROBBITS_STORE^p1);
#endif
	//*pcell=(((rbit<<PROBBITS_STORE)-cell)>>7)+cell;
	//*pcell=((rbit<<PROBBITS_STORE)-(1<<PROBBITS_STORE>>1)-(cell>>11))+cell;
	int count=cell&0x7FF;
	*pcell=(((rbit<<PROBBITS_STORE)-(1<<PROBBITS_STORE>>1)-(cell>>11))*(40000/(256+3*count))&0xFFFFF800)+cell+(count<0x7FF);
	ac->range=rbit?r2:ac->range;
	ac->low=rbit?ac->low:mid;
#ifdef ESTIMATE_BITSIZE
	bitsizes[ekc][eidx]+=shannontable[rbit?p1:(1<<PROBBITS_USE)-p1];
	++bitctr[ekc][eidx][rbit];
	winctr[ekc][eidx]+=rbit==(p1>=1<<PROBBITS_USE);
#endif
}
#endif
INLINE void mainloop(int iw, int ih, int bestrct, uint8_t *image, uint8_t *stream, ACState *ac, const int fwd)
{
	int
		yidx=rct_combinations[bestrct][II_PERM_Y],
		uidx=rct_combinations[bestrct][II_PERM_U],
		vidx=rct_combinations[bestrct][II_PERM_V],
		cu0=rct_combinations[bestrct][II_COEFF_U_SUB_Y],
		cv0=rct_combinations[bestrct][II_COEFF_V_SUB_Y],
		cv1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
	int32_t ky, kx;
	int32_t psize=0;
	int16_t *pixels=0;
	ALIGN(32) int32_t coeffs[3][L1NPREDS]={0}, bias[3]={0};
	uint8_t *imptr=image;
#ifndef USE_HIST
#ifdef PRINT_STATS
	int64_t unarysum=0, unarycount=0;
	uint32_t statctr[GRLIMIT+1][2]={0};
#endif
#endif
	
#ifndef USE_HIST
	(void)memusage;
#endif
	psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error\n");
		free(image);
		free(stream);
		return;
	}
	memset(pixels, 0, psize);
#ifdef USE_HIST
	memset(stats, 0, sizeof(stats));
#else
	FILLMEM_S((uint32_t*)stats1, 1<<PROBBITS_STORE>>1, sizeof(stats1), sizeof(int32_t));
	FILLMEM_S((uint32_t*)stats2, 1<<PROBBITS_STORE>>1, sizeof(stats2), sizeof(int32_t));
	FILLMEM_S((uint32_t*)stats3, 1<<PROBBITS_STORE>>1, sizeof(stats3), sizeof(int32_t));
	//memset(stats1, 0, sizeof(stats1));
	//memset(stats2, 0, sizeof(stats2));
	//memset(stats3, 0, sizeof(stats3));
#endif
	FILLMEM_S((int32_t*)coeffs, (1<<L1SH)/L1NPREDS, sizeof(coeffs), sizeof(int32_t));
	bias[0]=1<<L1SH>>1;
	bias[1]=1<<L1SH>>1;
	bias[2]=1<<L1SH>>1;
#ifdef ESTIMATE_BITSIZE
	memset(bitsizes, 0, sizeof(bitsizes));
	memset(bitctr, 0, sizeof(bitctr));
	memset(winctr, 0, sizeof(winctr));
#endif

	for(ky=0;ky<ih;++ky)
	{
		uint32_t yuv[4]={0};
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(kx=0;kx<iw;++kx, imptr+=3)
		{
			int kc;
			int offset, offset0;

			if(fwd)
			{
				yuv[0]=imptr[yidx];
				yuv[1]=imptr[uidx];
				yuv[2]=imptr[vidx];
			}
			offset0=offset=0;
			for(kc=0;kc<3;++kc)
			{
				int32_t
					NNNN	=rows[0][0+0*NCH*NROWS*NVAL],
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
					NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
					WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
					WWW	=rows[0][0-3*NCH*NROWS*NVAL],
					WW	=rows[0][0-2*NCH*NROWS*NVAL],
					W	=rows[0][0-1*NCH*NROWS*NVAL],
					eN	=rows[1][1+0*NCH*NROWS*NVAL],
					eNE	=rows[1][1+1*NCH*NROWS*NVAL],
					eNEE	=rows[1][1+2*NCH*NROWS*NVAL],
					eNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
					eW	=rows[0][1-1*NCH*NROWS*NVAL];
				int32_t pred;
#ifndef USE_HIST
				int32_t upred2;
				int epred;
#endif
				int32_t vmax, vmin, pred0;
				int32_t error;
#ifndef USE_HIST
				int32_t nbypass, nbypass0, grflag;
#endif
				int32_t nzeros=-1;
				int32_t tidx=0;
				int32_t *statsptr;
#ifndef USE_HIST
				int32_t bit=0;
				int32_t ctx;
#endif
				ALIGN(32) int32_t preds[L1NPREDS];
				
#if 1
				(void)NNNN	;
				(void)NNN	;
				(void)NN	;
				(void)NNE	;
				(void)NW	;
				(void)N		;
				(void)NE	;
				(void)NEE	;
				(void)NEEE	;
				(void)NEEEE	;
				(void)WWWW	;
				(void)WWW	;
				(void)WW	;
				(void)W		;
				(void)eN	;
				(void)eNE	;
				(void)eNEE	;
				(void)eNEEE	;
				(void)eW	;
#endif
				{
					int j;

					pred=bias[kc];
#define PRED(EXPR) preds[j]=EXPR; pred+=coeffs[kc][j]*preds[j]; ++j;
					j=0;
					PREDLIST;
#undef  PRED
#ifndef USE_HIST
					upred2=pred>>(L1SH-2);
#endif
					pred>>=L1SH;
				}
#ifndef USE_HIST
				upred2+=offset0&~3;
				CLAMP2(upred2, 0, 1023);
#endif
				pred0=(int32_t)pred;
				vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);
				pred+=offset;
				CLAMP2(pred, 0, 255);

				ctx=31-_lzcnt_u32(eW*eW+2);
				nbypass=(ctx>>1)-GRBITS;
				if(ctx>NCTX-1)
					ctx=NCTX-1;
				CLAMP2(nbypass, 0, 7);
#ifdef USE_HIST
				statsptr=stats[kc][pred][0];
			//	statsptr=stats[kc][pred][ctx];
				{
#if 1
					int x1=0;
					for(int kb=7, tidx=1;kb>=0;--kb)
					{
						int bit=yuv[kc]>>kb&1;
						int32_t cell=statsptr[tidx];
						int p1=cell>>(PROBSHIFT+11);
						//p1=squash(p1>>PROBSHIFT);
						//p1>>=PROBSHIFT;
						p1+=p1<0;
						p1+=1<<PROBBITS_USE>>1;
						//if(p1<1||p1>(1<<PROBBITS_USE)-1)
						//	CRASH("");
						codebit(ac, &bit, p1, fwd);
						x1+=bit<<kb;
						int count=cell&0x7FF;
						statsptr[tidx]=(((bit<<PROBBITS_STORE)-(1<<PROBBITS_STORE>>1)-(cell>>11))*(30000/(512+3*count))&0xFFFFF800)+cell+(count<0x7FF);
					//	statsptr[tidx]+=((bit<<PROBBITS_STORE)-(int16_t)statsptr[tidx]-(1<<PROBBITS_STORE>>1)+(1<<3>>1))>>3;
					//	statsptr[tidx]+=((bit<<PROBBITS_STORE)-(int16_t)statsptr[tidx]-(1<<PROBBITS_STORE>>1)+(1<<6>>1))>>6;
						tidx=2*tidx+bit;
					}
#else
					int x1=0, x2=256;
					while(x2-x1>1)
					{
						int bit, p1, c0, c1, c2, mid;

						c0=currCDF[x1]+x1;
						c2=currCDF[x2]+x2;

						mid=(x1+x2)>>1;
						//c1=(c0+c2)>>1;
						//for(mid=x1+1;mid<x2-1&&currCDF[mid]+mid<c1;++mid);

						c1=currCDF[mid]+mid;
						p1=((c2-c1)<<PROBBITS_USE)/(c2-c0);
						CLAMP2(p1, 1, (1<<PROBBITS_USE)-1);
						//p1+=p1<1<<PROBBITS_USE>>1;
						bit=yuv[kc]>=(uint32_t)mid;
						codebit(ac, &bit, p1, fwd);
						x1=bit?mid:x1;
						x2=bit?x2:mid;
					}
					for(int k=x1+1;k<257;++k)
						++currCDF[k];
					if(currCDF[256]>=0xFFFF-256)
					{
						for(int k=0;k<257;++k)
							currCDF[k]>>=1;
					}
#endif
#ifdef _MSC_VER
					if(fwd&&x1!=yuv[kc])
						CRASH("");
#endif
					yuv[kc]=x1;
#ifdef ENABLE_GUIDE
					if(!fwd)
					{
						uint8_t *pval=&g_image[3*(iw*ky+kx)+rct_combinations[bestrct][II_PERM_Y+kc]];
						uint8_t val=*pval;
						uint8_t pixel=yuv[kc];
						if(pixel!=val)
							CRASH("guide  YXC %d %d %d", ky, kx, kc);
					}
#endif
				}
#else
				epred=128-abs(pred-128);
				if(fwd)
				{
					error=yuv[kc]-(int32_t)pred;
					int negmask=error>>31;
					int abserr=(error^negmask)-negmask;
					error=error<<1^negmask;
					if(epred<abserr)
						error=epred+abserr;
					if(error==256)
					{
						error=(int8_t)(yuv[kc]-pred);
						error=error<<1^error>>31;
					}
					nzeros=error>>nbypass;
#ifdef PRINT_STATS
					unarysum+=nzeros;
					++unarycount;
#endif
				}
				else
					error=0;
				statsptr=stats1[kc][upred2][ctx];
				tidx=0;
#ifdef ESTIMATE_BITSIZE
				ekc=kc;
#endif
				do
				{
					bit=tidx>=nzeros;
#ifdef ESTIMATE_BITSIZE
					eidx=tidx;
					++unary_count;
#endif
					codebit(ac, statsptr+tidx, &bit, fwd);
#ifdef PRINT_STATS
					++statctr[tidx][bit];
#endif
					if(bit)
						break;
					++tidx;
				}while(tidx<GRLIMIT);
				nbypass0=nbypass;
				grflag=tidx==GRLIMIT;
				if(grflag)
				{
					error-=(GRLIMIT-1)<<nbypass;
					statsptr=stats3[kc][nbypass];
					tidx=1;
					nbypass=8;
				}
				else
				{
					statsptr=stats2[kc][upred2];
					tidx=(256>>nbypass)+tidx;//bit coding:  tidx=2*tidx+bit  tidx=0b1XX
				}
				{
					int32_t kb=nbypass-1;

					for(;kb>=0;--kb)
					{
						if(fwd)
							bit=error>>kb&1;
#ifdef ESTIMATE_BITSIZE
						eidx=GRLIMIT+8-nbypass+kb;
						++binary_count;
#endif
						codebit(ac, statsptr+tidx, &bit, fwd);
						tidx=2*tidx+bit;
					}
				}
				if(grflag)
					tidx+=(GRLIMIT-1)<<nbypass0;
#ifdef _MSC_VER
				if(fwd&&grflag)
					error+=(GRLIMIT-1)<<nbypass0;
				if(fwd&&tidx!=error+256)
					CRASH("");
#endif
				if(!fwd)
				{
					error=(uint8_t)tidx;
					if(2*pred+error==512)
					{
						error=error>>1^-(error&1);
						yuv[kc]=(uint8_t)(error+pred);
					}
					else
					{
						int negmask=((int32_t)pred-128)>>31;
						int sym=error;
						int e2=epred-sym;
						error=sym>>1^-(sym&1);
						e2=(e2^negmask)-negmask;
						if((epred<<1)<sym)
							error=e2;
						yuv[kc]=error+(int32_t)pred;
					}
#ifdef ENABLE_GUIDE
					{
						uint8_t *pval=&g_image[3*(iw*ky+kx)+rct_combinations[bestrct][II_PERM_Y+kc]];
						uint8_t val=*pval;
						uint8_t pixel=yuv[kc];
						if(pixel!=val)
							CRASH("GUIDE YXC %d %d %d", ky, kx, kc);
					}
#endif
				}
#endif
				{
					int32_t curr=yuv[kc]-offset;
					int32_t k, e;

					e=(curr>pred0)-(curr<pred0);

					bias[kc]+=e<<8;
#define PRED(EXPR) coeffs[kc][k]+=e*preds[k]; ++k;
					k=0;
					PREDLIST;
#undef  PRED

					error=yuv[kc]-(int32_t)pred;
					error=error<<1^error>>31;
					rows[0][0]=curr;

					rows[0][1]=(2*eW+(error<<GRBITS)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
					offset0=kc ? cv0*yuv[0]+cv1*yuv[1] : cu0*yuv[0];
					offset=offset0>>2;
				}
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
			}
			if(!fwd)
			{
				imptr[yidx]=yuv[0];
				imptr[uidx]=yuv[1];
				imptr[vidx]=yuv[2];
#ifdef ENABLE_GUIDE
				guide_check(image, kx, ky);
#endif
			}
		}
	}
	free(pixels);
#ifndef USE_HIST
#ifdef PRINT_STATS
	if(fwd)
	{
		printf("unary: zeros %lld vs ones %lld  %8.4lf%%\n", unarysum, unarycount, 100.*unarysum/(unarysum+unarycount));
		for(int k=0;k<GRLIMIT+1;++k)
		{
			int sum=statctr[k][0]+statctr[k][1];
			double p0=(double)statctr[k][0]/sum;
			double csize0=-(statctr[k][0]*log2(p0)+statctr[k][1]*log2(1-p0));
			printf("%3d  %9d vs %9d  %8.4lf%%  %12.2lf -> %12.2lf  %8.4lf%%\n"
				, k
				, statctr[k][0]
				, statctr[k][1]
				, 100.*p0
				, sum/8.
				, csize0/8
				, 100.*csize0/sum
			);
		}
	}
#endif
#endif
}
int c58_codec(int argc, char **argv)
{
	static const uint16_t tag='5'|'8'<<8;

	const char *srcfn=0, *dstfn=0;
	ptrdiff_t srcsize=0, dstsize=0;
	int fwd=0;
	int32_t iw=0, ih=0;
	ptrdiff_t res=0, usize=0, csize=0;
	uint8_t *image=0, *imptr=0, *imend=0, *stream=0, *streamptr=0, *streamend=0;
	ptrdiff_t streamsize=0;
	int bestrct=0;
	ACState ac=
	{
		0, 0xFFFFFFFFFFFF, 0,
		0, 0,
	};
#ifdef LOUD
	double t=time_sec2();
#endif
#ifdef PROFILER
	void *prof_ctx=prof_start();
#endif

	if(argc!=3&&argc!=4)
	{
		printf(
			"Usage:  \"%s\"  input  output\n"
			, argv[0]
		);
		return 1;
	}
	srcfn=argv[1];
	dstfn=argv[2];
	
#ifdef ESTIMATE_BITSIZE
	for(int k=0;k<1<<PROBBITS_USE;++k)
	{
		double val=(double)(k+(k<1<<PROBBITS_USE>>1))*(1./(1<<PROBBITS_USE));
		val=-log(val)*(1./(8*M_LN2));
		shannontable[k]=val;
	}
#endif
	//read source
	{
		struct stat info={0};
		int error=stat(srcfn, &info);
		if(error)
		{
			CRASH("Cannot stat \"%s\"", srcfn);
			return 1;
		}
		srcsize=info.st_size;
	}
	{
		FILE *fsrc;
		ptrdiff_t nread;
		int c;
		
		fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			CRASH("Cannot open \"%s\"", srcfn);
			return 1;
		}
		c=0;
		fread(&c, 1, 2, fsrc);
		fwd=c==('P'|'6'<<8);
		if(!fwd&&c!=tag)
		{
			CRASH("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
		{
			c=fgetc(fsrc);
			if(c!='\n')
			{
				CRASH("Invalid PPM file");
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
				iw=10*iw+c-'0';
				c=fgetc(fsrc);
			}
			while(c<=' ')
				c=fgetc(fsrc);
			ih=0;
			while((uint32_t)(c-'0')<10)
			{
				ih=10*ih+c-'0';
				c=fgetc(fsrc);
			}
			while(c<=' ')
				c=fgetc(fsrc);
			while(c=='#')
			{
				c=fgetc(fsrc);
				while(c!='\n')
					c=fgetc(fsrc);
				c=fgetc(fsrc);
			}
			c=c<<8|fgetc(fsrc);
			c=c<<8|fgetc(fsrc);
			c=c<<8|fgetc(fsrc);
			if(c!=('2'<<24|'5'<<16|'5'<<8|'\n'))
			{
				CRASH("Unsupported PPM file");
				return 1;
			}
		}
		else
		{
			iw=0;
			ih=0;
			bestrct=0;
			fread(&iw, 1, 3, fsrc);
			fread(&ih, 1, 3, fsrc);
			fread(&bestrct, 1, 1, fsrc);
			nread=ftell(fsrc);
			streamsize=srcsize-nread;

		}
		if(iw<1||ih<1)
		{
			CRASH("Unsupported image dimensions  WH %d*%d", iw, ih);
			return 1;
		}
		res=(ptrdiff_t)iw*ih;
		usize=3*res;
		if(fwd)
			streamsize=usize;
		image=(uint8_t*)malloc(usize);
		stream=(uint8_t*)malloc(streamsize+sizeof(char[32]));
		if(!image||!stream)
		{
			CRASH("Alloc error");
			return 1;
		}
		imend=image+usize;
		{
			ptrdiff_t expected=0;
			if(fwd)
			{
				expected=usize;
				nread=fread(image, 1, usize, fsrc);
				guide_save(image, iw, ih);
			}
			else
			{
				expected=streamsize;
				nread=fread(stream, 1, streamsize, fsrc);
			}
			if(nread!=expected)
				printf("Truncated  expected %td  read %td", expected, nread);
		}
		fclose(fsrc);
	}
	if(fwd)
	{
		//analysis
		int64_t counters[OCH_COUNT]={0}, minerr=0;
		int prev[OCH_COUNT]={0};
		int rowstride=3*iw;

		imptr=image+rowstride;
		while(imptr<imend)
		{
			int r, g, b, rg, gb, br;

			r=(imptr[0]-imptr[0-rowstride])<<2;
			g=(imptr[1]-imptr[1-rowstride])<<2;
			b=(imptr[2]-imptr[2-rowstride])<<2;
			imptr+=3;
			rg=r-g;
			gb=g-b;
			br=b-r;
#define UPDATE(I0, E0, I1, E1, I2, E2)\
	do\
	{\
		int t0=E0;\
		int t1=E1;\
		int t2=E2;\
		counters[I0]+=abs(t0-prev[I0]);\
		counters[I1]+=abs(t1-prev[I1]);\
		counters[I2]+=abs(t2-prev[I2]);\
		prev[I0]=t0;\
		prev[I1]=t1;\
		prev[I2]=t2;\
	}while(0)

			UPDATE(
				OCH_Y400, r,
				OCH_Y040, g,
				OCH_Y004, b
			);
			UPDATE(
				OCH_CX40, rg,
				OCH_C0X4, gb,
				OCH_C40X, br
			);
			UPDATE(
				OCH_CX31, rg+(gb>>2),//r-(3*g+b)/4 = r-g-(b-g)/4
				OCH_C3X1, rg+(br>>2),//g-(3*r+b)/4 = g-r-(b-r)/4
				OCH_C31X, br+(rg>>2) //b-(3*r+g)/4 = b-r-(g-r)/4
			);
			UPDATE(
				OCH_CX13, br+(gb>>2),//r-(g+3*b)/4 = r-b-(g-b)/4
				OCH_C1X3, gb+(br>>2),//g-(r+3*b)/4 = g-b-(r-b)/4
				OCH_C13X, gb+(rg>>2) //b-(r+3*g)/4 = b-g-(r-g)/4
			);
			UPDATE(
				OCH_CX22, (rg-br)>>1,//r-(g+b)/2 = (r-g + r-b)/2
				OCH_C2X2, (rg-gb)>>1,//g-(r+b)/2 = (g-r + g-b)/2
				OCH_C22X, (br-gb)>>1 //b-(r+g)/2 = (b-r + b-g)/2
			);
#undef  UPDATE
		}
		{
			int kt;

#ifdef PRINT_RCT
			for(kt=0;kt<OCH_COUNT;++kt)
				printf("%d %16lld\n", kt, counters[kt]);
			printf("\n");
#endif
			for(kt=0;kt<RCT_COUNT;++kt)
			{
				const uint8_t *combination=rct_combinations[kt];
				long long currerr=
					+counters[combination[0]]
					+counters[combination[1]]
					+counters[combination[2]]
				;
				if(!kt||minerr>currerr)
				{
					minerr=currerr;
					bestrct=kt;
				}
#ifdef PRINT_RCT
				printf("RCT%02d %16lld%s\n", kt, currerr, kt==bestrct?" <-":"");
#endif
			}
#ifdef LOUD
			printf("WH %d*%d  %lld B  RCT %d %s  \"%s\"\n"
				, iw, ih, usize
				, bestrct, rct_names[bestrct]
				, srcfn
			);
#else
			(void)och_names;
			(void)rct_names;
#endif
		}
		streamptr=stream;
		streamend=stream+usize;

		ac.ptr=streamptr;
		ac.end=streamend;
	}
	else
	{
		streamptr=stream;
		streamend=stream+srcsize;

		ac.code=*(uint64_t*)streamptr;//load
		streamptr+=sizeof(uint64_t);
		ac.code=ac.code<<32|ac.code>>32;
		ac.ptr=streamptr;
		ac.end=streamend;

		csize=srcsize;
	}
#ifdef _MSC_VER
	ac.totalbits=(int64_t)24*iw*ih;
#endif
	if(fwd)
		mainloop(iw, ih, bestrct, image, stream, &ac, 1);
	else
		mainloop(iw, ih, bestrct, image, stream, &ac, 0);
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing\n", dstfn);
			free(image);
			free(stream);
			return 1;
		}
		if(fwd)
		{
			*(uint64_t*)ac.ptr=ac.low<<32|ac.low>>32;//flush
			ac.ptr+=sizeof(uint64_t);

			csize=ac.ptr-stream;
			dstsize+=fwrite(&tag, 1, 2, fdst);
			dstsize+=fwrite(&iw, 1, 3, fdst);
			dstsize+=fwrite(&ih, 1, 3, fdst);
			dstsize+=fwrite(&bestrct, 1, 1, fdst);
			dstsize+=fwrite(stream, 1, csize, fdst);
			csize=dstsize;
		}
		else
		{
			dstsize+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			dstsize+=fwrite(image, 1, usize, fdst);
			usize=dstsize;
		}
		fclose(fdst);
	}
	free(image);
	free(stream);
#ifdef LOUD
	t=time_sec2()-t;
	if(fwd)
	{
		usize=srcsize;
#ifdef ESTIMATE_BITSIZE
		double total[3]={0}, total_u=0, total_b=0;
		printf("plane  csize / usize = invCR  nzeros%%\n");
		for(int kv=0;kv<GRLIMIT+8;++kv)
		{
			printf("%3d %12.2lf/%12.2lf=%8.4lf%% %8.4lf%%  %12.2lf/%12.2lf=%8.4lf%% %8.4lf%%  %12.2lf/%12.2lf=%8.4lf%% %8.4lf%%  %8.4lf%% %8.4lf%% %8.4lf%%\n"
				, kv
				, bitsizes[0][kv], (bitctr[0][kv][0]+bitctr[0][kv][1])/8., 800.*bitsizes[0][kv]/(bitctr[0][kv][0]+bitctr[0][kv][1]), 100.*bitctr[0][kv][1]/(bitctr[0][kv][0]+bitctr[0][kv][1])
				, bitsizes[1][kv], (bitctr[1][kv][0]+bitctr[1][kv][1])/8., 800.*bitsizes[1][kv]/(bitctr[1][kv][0]+bitctr[1][kv][1]), 100.*bitctr[1][kv][1]/(bitctr[1][kv][0]+bitctr[1][kv][1])
				, bitsizes[2][kv], (bitctr[2][kv][0]+bitctr[2][kv][1])/8., 800.*bitsizes[2][kv]/(bitctr[2][kv][0]+bitctr[2][kv][1]), 100.*bitctr[2][kv][1]/(bitctr[2][kv][0]+bitctr[2][kv][1])
				, 100.*winctr[0][kv]/(bitctr[0][kv][0]+bitctr[0][kv][1])
				, 100.*winctr[1][kv]/(bitctr[1][kv][0]+bitctr[1][kv][1])
				, 100.*winctr[2][kv]/(bitctr[2][kv][0]+bitctr[2][kv][1])
			);
			total[0]+=bitsizes[0][kv];
			total[1]+=bitsizes[1][kv];
			total[2]+=bitsizes[2][kv];
			if(kv==GRLIMIT-1||kv==GRLIMIT+8-1)
			{
				printf("TYUV %12.2lf %12.2lf %12.2lf %12.2lf\n\n"
					, total[0]+total[1]+total[2]
					, total[0]
					, total[1]
					, total[2]
				);
				if(kv==GRLIMIT-1)
					total_u=total[0]+total[1]+total[2];
				else
					total_b=total[0]+total[1]+total[2];
				total[2]=total[1]=total[0]=0;
			}
		}
		printf("%12.2lf  %12.6lf  ->  %12.2lf  %12.6lf  unary\n", unary_count/8., (double)unary_count/usize, total_u, total_u*8/usize);
		printf("%12.2lf  %12.6lf  ->  %12.2lf  %12.6lf  binary\n", binary_count/8., (double)binary_count/usize, total_b, total_b*8/usize);
		printf("%12.2lf  %12.6lf  ->  %12.2lf  %12.6lf  total\n", (unary_count+binary_count)/8., (double)(unary_count+binary_count)/usize, total_u+total_b, (total_u+total_b)*8/usize);
		printf("\n");
#endif
		printf("%9td->%9td  %8.4lf%%  %12.6lf:1  BPD %12.6lf\n"
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
	(void)dstsize;
	(void)csize;
	(void)&time_sec2;
	//(void)&squash;
#ifdef PROFILER
	prof_end(prof_ctx);
#endif
	return 0;
}
