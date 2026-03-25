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
	#define ESTIMATE_BITSIZE
//	#define PRINTBITS

	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif

//	#define NEW_SYSTEM//X
	#define SIMD_L1
	#define USE_TABLES
	#define USE_COUNTERS
//	#define ZIPF_VIEW

	#define UNSIGNED_PIXEL


#if 1
#define PREDLIST_LOSSY\
	PRED(N)\
	PRED(NNN)\
	PRED(NNNN)\
	PRED(N+yN)\
	PRED(N+yNE)\
	PRED(N+yNW)\
	PRED(N+yW)\
	PRED(W)\
	PRED(WWW)\
	PRED(WWWW)\
	PRED(W+xW)\
	PRED(W+xNE)\
	PRED(W+xNEE)\
	PRED(W+xNW)\
	PRED(N+2*yN-yNN)\
	PRED(W+2*xW-xWW)\
	PRED(NE)\
	PRED(NEEE)\
	PRED(NEEEE)\

#define PREDLIST\
	PRED(N)\
	PRED(NNN)\
	PRED(NNNN)\
	PRED(N+yN)\
	PRED(N+yNE)\
	PRED(N+yNW)\
	PRED(N+yW)\
	PRED(W)\
	PRED(WWW)\
	PRED(WWWW)\
	PRED(W+xW)\
	PRED(W+xNE)\
	PRED(W+xNEE)\
	PRED(W+xNW)\
	PRED(W+2*xW-xWW)\
	PRED(NEEEE)\

#if 0
	PRED(N+2*yN-yNN)\
	PRED(NEEE)\

#endif
#endif
#if 0
#define PREDLIST\
	PRED(N)\
	PRED(W)\
	PRED(NNN)\
	PRED(WWW)\
	PRED(NNNN)\
	PRED(WWWW)\
	PRED(NEEE)\
	PRED(NEEEE)\
	PRED(N+yN)\
	PRED(W+xW)\
	PRED(N+yW)\
	PRED(W+xNE)\
	PRED(W+xNEE)\
	PRED(N+yNE)\
	PRED(N+yNW)\
	PRED(W+xNW)\

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
enum
{
	L1SH_LOSSY=18,
	L1SH=20,
#define PRED(EXPR) +1
	L1NPREDS=PREDLIST,
	L1NPREDS_LOSSY=PREDLIST_LOSSY,
#undef  PRED

	GRBITS=6,
	NCTX=24,
	GRLIMIT=18,
	PROBBITS_STORE=15,
	PROBBITS_USE=14,
//	PROBBITS_STORE=17,
//	PROBBITS_USE=12,
	PROBSHIFT=PROBBITS_STORE-PROBBITS_USE,

	XPAD=8,
	NCH=3,
	NROWS=4,
	NVAL=4,
	NVAL0=2,
#ifdef USE_COUNTERS
	CTRBITS=9,
	CTRMASK=(1<<CTRBITS)-1,

	HISTBITS=22,// <= 64-CTRBITS*2
	HISTMASK=(1<<HISTBITS)-1,

//	CTRFBITS=(32-CTRBITS)>>1,
//	CTRFBITS=4,
//	CTRFMASK=(1<<CTRFBITS)-1,
#endif
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
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
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

#ifdef USE_COUNTERS
static uint64_t stats1[3][1024][NCTX][GRLIMIT];//unary
static uint64_t stats2[3][1024][256];//remainder
static uint64_t stats3[3][8][256];//bypass on GRLIMIT
static const size_t memusage=sizeof(stats1)+sizeof(stats2)+sizeof(stats3);
#else
#ifdef NEW_SYSTEM
static int32_t stats1[3][1024][NCTX][GRLIMIT+2];//unary
#else
static int32_t stats1[3][1024][NCTX][GRLIMIT];//unary
#endif
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
static double shannontable[1<<PROBBITS_USE];
static double bitsizes[3][GRLIMIT+8];
static uint32_t bitctr[3][GRLIMIT+8][2];
static uint32_t winctr[3][GRLIMIT+8];
static int ekc, eidx;
#endif
#ifdef ZIPF_VIEW
static double zsize=0;
static uint8_t *zptr=0;
#endif
#ifdef USE_TABLES
static uint8_t epredtable[256];
static uint8_t clamptable[512];
static uint16_t errortable[512][2];
static uint8_t ctxtable[(2<<NCTX/2)/3][2];
#endif
#ifdef USE_COUNTERS
//static int32_t table;
#endif
//#ifdef _MSC_VER
//static uint64_t unary_count=0, binary_count=0;
//#endif
static int squash(int32_t d)
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
INLINE void codebit(ACState *ac, uint64_t *pcell, int32_t *bit, const int fwd)
{
	int rbit;
	uint64_t r2, mid;
	
#ifdef USE_COUNTERS
#if 1
	uint64_t cell=*pcell;
	int32_t n[]={cell&CTRMASK, cell>>CTRBITS&CTRMASK};
	uint64_t hist=cell>>CTRBITS*2&HISTMASK;
	int32_t alpha=(int32_t)(cell>>(HISTBITS+CTRBITS*2)&0xFFFF);
	int n0e=n[0]+2;
	int n1e=n[1]+2;
	int sum=n0e+n1e;
	int hwt=(((int)_mm_popcnt_u64(hist)+1)<<PROBBITS_USE)/(HISTBITS+2);
	int32_t p1=((n1e<<PROBBITS_USE)+(sum>>1))/sum;
	int32_t x=hwt-p1;
	p1+=x*alpha>>16;
#endif
#else
	int32_t p1=*pp1a>>PROBSHIFT;
//	int32_t p10a=*pp1a;
//	int32_t p1=p10a>>sh;

	p1+=p1<(1<<PROBBITS_USE>>1);
//	CLAMP2(p1, 1, (1<<PROBBITS_USE)-1);
//	p1=squash(p1-(1<<PROBBITS_USE>>1));
#endif
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
	if(fwd)
		rbit=*bit;
	else
		*bit=rbit=ac->code<mid;
#ifdef FIFOVAL
	//if(p1<1||p1>(1<<PROBBITS_USE)-1)
	//	CRASH("");
	if(fwd)
		fifoval_enqueue(rbit<<PROBBITS_STORE^p1);
	else
		fifoval_check(rbit<<PROBBITS_STORE^p1);
#endif
	if(rbit)
		ac->range=r2;
	else
		ac->low=mid;

#ifdef ZIPF_VIEW
	if(fwd)
		zsize-=log2((double)(rbit?p1:(1<<PROBBITS_USE)-p1)/(1<<PROBBITS_USE));
#endif
#ifdef USE_COUNTERS
#if defined _MSC_VER && 1
	static int ctrctr=0;
	++ctrctr;
	if((uint32_t)(ctrctr-4000000)<5000)
	{
		printf(" %d %d %04X %04X %3d %3d ", ctrctr, rbit, p1, alpha, n[0], n[1]);
		for(int kb=HISTBITS-1;kb>=0;--kb)
			printf("%d", (int)(hist>>kb&1));
		if((ctrctr&3)==3)
			printf("\n");
	}
#endif
#if 1
	++n[rbit];
	if(n[rbit]>CTRMASK)
	{
		n[0]>>=1;
		n[1]>>=1;
	}
	hist=(uint64_t)rbit<<(HISTBITS-1)|hist>>1;
	alpha+=(x<<9)/(p1-(!rbit<<PROBBITS_USE));
	CLAMP2(alpha, 1, 0xFFFF);
	*pcell=(uint64_t)alpha<<(HISTBITS+CTRBITS*2)|(uint64_t)hist<<CTRBITS*2|(uint64_t)n[1]<<CTRBITS|n[0];
#endif
#else
//	p10a+=(int32_t)((rbit<<PROBBITS_STORE)-p10a)>>sh;
//	p10a+=(int32_t)((rbit<<PROBBITS_STORE)-p10a)>>7;
//	*pp1a=p10a;
#endif
#ifdef ESTIMATE_BITSIZE
	bitsizes[ekc][eidx]+=shannontable[rbit?p1:(1<<PROBBITS_USE)-p1];
	//if(isinf(bitsizes[ekc][eidx]))
	//	printf("");
	++bitctr[ekc][eidx][rbit];
	winctr[ekc][eidx]+=rbit==(p1>=1<<PROBBITS_USE);
	//if(bitsizes[ekc][eidx]>ac->bitidx/8.*1.2)
	//	CRASH("");
#endif
#ifdef _MSC_VER
	++ac->n[rbit];
#endif
	//return p1;
}
INLINE void mainloop(int iw, int ih, int bestrct, int dist, uint8_t *image, uint8_t *stream, ACState *ac, int32_t *gaintable2, const int lossy, const int fwd)
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
	ALIGN(32) int32_t lcoeffs[3][L1NPREDS_LOSSY+1]={0}, coeffs[3][L1NPREDS]={0}, bias[3]={0};
	int32_t invdist=((1<<16)+dist-1)/dist;
	uint8_t *imptr=image;
#ifdef PRINTBITS
	ptrdiff_t idx=0, usize=(ptrdiff_t)3*iw*ih;
#endif
	
	//srand(10);
	(void)memusage;
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
#ifdef USE_COUNTERS
	memset(stats1, 0, sizeof(stats1));
	memset(stats2, 0, sizeof(stats2));
	memset(stats3, 0, sizeof(stats3));
	//memset(hist1, 0xAA, sizeof(hist1));
	//memset(hist2, 0xAA, sizeof(hist2));
	//memset(hist3, 0xAA, sizeof(hist3));
#else
	FILLMEM_S((uint32_t*)stats1, 1<<PROBBITS_STORE>>1, sizeof(stats1), sizeof(int32_t));
	FILLMEM_S((uint32_t*)stats2, 1<<PROBBITS_STORE>>1, sizeof(stats2), sizeof(int32_t));
	FILLMEM_S((uint32_t*)stats3, 1<<PROBBITS_STORE>>1, sizeof(stats3), sizeof(int32_t));
#endif
	if(lossy)
	{
		FILLMEM_S((int32_t*)lcoeffs, (1<<L1SH_LOSSY)/L1NPREDS_LOSSY, sizeof(lcoeffs), sizeof(int32_t));
	}
	else
	{
		FILLMEM_S((int32_t*)coeffs, (1<<L1SH)/L1NPREDS, sizeof(coeffs), sizeof(int32_t));
		bias[0]=1<<L1SH>>1;
		bias[1]=1<<L1SH>>1;
		bias[2]=1<<L1SH>>1;
	}
#ifdef USE_TABLES
	uint8_t *const epredptr=epredtable;
	for(int k=0;k<256;++k)
		epredtable[k]=lossy ? k>=128?255-k:k : 128-abs(k-128);
	uint8_t *const clampptr=clamptable+128;
	for(int k=0;k<_countof(clamptable);++k)
	{
		int val=k-128;
		CLAMP2(val, 0, 255);
		clamptable[k]=val;
	}
	uint16_t (*const errorptr)[2]=errortable+256;
	for(int k=0;k<512;++k)
	{
		int error=k-256;
		errortable[k][0]=error<<1^error>>31;
		errortable[k][1]=abs(error);
	}
	for(int k=0;k<_countof(ctxtable);++k)
	{
		int ctx=31-_lzcnt_u32(k*k+2);
		int nbypass=(ctx>>1)-GRBITS;
		CLAMP2(nbypass, 0, 7);
		if(ctx>NCTX-1)
			ctx=NCTX-1;
		ctxtable[k][0]=ctx;
		ctxtable[k][1]=nbypass;
	}
#endif
#ifdef ESTIMATE_BITSIZE
	memset(bitsizes, 0, sizeof(bitsizes));
	memset(bitctr, 0, sizeof(bitctr));
	memset(winctr, 0, sizeof(winctr));
#endif
#ifdef ZIPF_VIEW
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	uint8_t *zimage=0;
	if(fwd)
	{
		zimage=(uint8_t*)malloc(res);
		if(!zimage)
		{
			CRASH("Alloc error");
			return;
		}
		memset(zimage, 0, res);
		zptr=zimage;
	}
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
#ifdef UNSIGNED_PIXEL
				yuv[0]=imptr[yidx];
				yuv[1]=imptr[uidx];
				yuv[2]=imptr[vidx];
#else
				yuv[0]=imptr[yidx]-128;
				yuv[1]=imptr[uidx]-128;
				yuv[2]=imptr[vidx]-128;
#endif
			}
			offset0=offset=0;
			for(kc=0;kc<3;++kc)
			{
				int32_t
					NNNN	=rows[0][0+0*NCH*NROWS*NVAL],
					NNN	=rows[3][0+0*NCH*NROWS*NVAL],
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
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
					xNW	=rows[1][1-1*NCH*NROWS*NVAL],
					xNE	=rows[1][1+1*NCH*NROWS*NVAL],
					xNEE	=rows[1][1+2*NCH*NROWS*NVAL],
					xNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
					xWW	=rows[0][1-2*NCH*NROWS*NVAL],
					xW	=rows[0][1-1*NCH*NROWS*NVAL],
					yNN	=rows[2][2+0*NCH*NROWS*NVAL],
					yNW	=rows[1][2-1*NCH*NROWS*NVAL],
					yN	=rows[1][2+0*NCH*NROWS*NVAL],
					yNE	=rows[1][2+1*NCH*NROWS*NVAL],
					yW	=rows[0][2-1*NCH*NROWS*NVAL],
					eNE	=rows[1][3+1*NCH*NROWS*NVAL],
					eNEE	=rows[1][3+2*NCH*NROWS*NVAL],
					eNEEE	=rows[1][3+3*NCH*NROWS*NVAL],
					eW	=rows[0][3-1*NCH*NROWS*NVAL];
#if 0
				int32_t
					NNNN	=rows[0][0+0*NCH*NROWS*NVAL],
				//	NNN	=rows[3][0+0*NCH*NROWS*NVAL],
					NNW	=rows[2][0-1*NCH*NROWS*NVAL],
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL],
					NWW	=rows[1][0-2*NCH*NROWS*NVAL],
					NW	=rows[1][0-1*NCH*NROWS*NVAL],
				//	N	=rows[1][0+0*NCH*NROWS*NVAL],
				//	NE	=rows[1][0+1*NCH*NROWS*NVAL],
					NEE	=rows[1][0+2*NCH*NROWS*NVAL],
				//	NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
				//	NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
					WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
				//	WWW	=rows[0][0-3*NCH*NROWS*NVAL],
					WW	=rows[0][0-2*NCH*NROWS*NVAL];
				//	W	=rows[0][0-1*NCH*NROWS*NVAL],
				//	eN	=rows[1][1+0*NCH*NROWS*NVAL],
				//	eNE	=rows[1][1+1*NCH*NROWS*NVAL],
				//	eNEE	=rows[1][1+2*NCH*NROWS*NVAL],
				//	eNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
				//	eW	=rows[0][1-1*NCH*NROWS*NVAL];
#endif
				int32_t pred;
				//int32_t upred;
				int32_t upred2;
				int32_t vmax, vmin, pred0;
				int32_t error;
				int32_t nbypass, nbypass0;
				int32_t nzeros=-1, grflag;
				int32_t tidx=0;
				uint64_t *statsptr;
				//uint64_t *histptr;
				int32_t bit=0;
				int32_t ctx;
				ALIGN(32) int32_t preds[L1NPREDS>L1NPREDS_LOSSY?L1NPREDS:L1NPREDS_LOSSY];

				{
					int j;

					if(lossy)
					{
						pred=1<<L1SH_LOSSY>>1;
#define PRED(EXPR) preds[j]=EXPR; pred+=lcoeffs[kc][j]*preds[j]; ++j;
						j=0;
						PREDLIST_LOSSY;
#undef  PRED
						upred2=pred>>(L1SH_LOSSY-2);
						pred>>=L1SH_LOSSY;
					}
					else
					{
						pred=bias[kc];
#ifdef SIMD_L1
#define PRED(EXPR) preds[j++]=EXPR;
						j=0;
						PREDLIST;
#undef  PRED
						__m256i mc0=_mm256_load_si256((__m256i*)coeffs[kc]+0);
						__m256i mc1=_mm256_load_si256((__m256i*)coeffs[kc]+1);
						__m256i mp0=_mm256_load_si256((__m256i*)preds+0);
						__m256i mp1=_mm256_load_si256((__m256i*)preds+1);

						__m256i lo0=_mm256_mul_epi32(mc0, mp0);
						__m256i lo1=_mm256_mul_epi32(mc1, mp1);
						__m256i hi0=_mm256_mul_epi32(_mm256_srli_epi64(mc0, 32), _mm256_srli_epi64(mp0, 32));
						__m256i hi1=_mm256_mul_epi32(_mm256_srli_epi64(mc1, 32), _mm256_srli_epi64(mp1, 32));
						mc0=_mm256_add_epi64(lo0, lo1);
						mc0=_mm256_add_epi64(mc0, hi0);
						mc0=_mm256_add_epi64(mc0, hi1);
						__m128i xc0=_mm_add_epi64(_mm256_extracti128_si256(mc0, 1), _mm256_castsi256_si128(mc0));
						xc0=_mm_add_epi64(xc0, _mm_shuffle_epi32(xc0, _MM_SHUFFLE(1, 0, 3, 2)));
						pred+=_mm_extract_epi32(xc0, 0);

#if 0
						int LOL_1=mc0.m256i_i32[0], LOL_2=mp0.m256i_i32[0];
						if(ky==ih/2&&kx==iw/2&&kc==1)//
							printf("");
						mp0=_mm256_blend_epi16(mp0, _mm256_slli_epi32(mp0, 16), 0xAA);
						mp1=_mm256_blend_epi16(mp1, _mm256_slli_epi32(mp1, 16), 0xAA);
						__m256i lo0=_mm256_mullo_epi16(mc0, mp0);
						__m256i hi0=_mm256_mulhi_epi16(mc0, mp0);
						__m256i lo1=_mm256_mullo_epi16(mc1, mp1);
						__m256i hi1=_mm256_mulhi_epi16(mc1, mp1);
						mc0=_mm256_add_epi16(lo0, _mm256_slli_epi32(hi0, 16));
						mc1=_mm256_add_epi16(lo1, _mm256_slli_epi32(hi1, 16));
						if(mc0.m256i_i32[0]!=LOL_1*LOL_2)
							CRASH("");
						mc0=_mm256_add_epi32(mc0, mc1);
						__m128i xc0=_mm_add_epi32(_mm256_extracti128_si256(mc0, 1), _mm256_castsi256_si128(mc0));
						xc0=_mm_add_epi32(xc0, _mm_shuffle_epi32(xc0, _MM_SHUFFLE(1, 0, 3, 2)));
						xc0=_mm_add_epi32(xc0, _mm_shuffle_epi32(xc0, _MM_SHUFFLE(2, 3, 0, 1)));
						pred+=_mm_extract_epi32(xc0, 0);
#endif
#if 0
						mc0=_mm256_mullo_epi32(mc0, mp0);
						mc1=_mm256_mullo_epi32(mc1, mp1);
						mc0=_mm256_add_epi32(mc0, mc1);
						__m128i xc0=_mm_add_epi32(_mm256_extracti128_si256(mc0, 1), _mm256_castsi256_si128(mc0));
						xc0=_mm_add_epi32(xc0, _mm_shuffle_epi32(xc0, _MM_SHUFFLE(1, 0, 3, 2)));
						xc0=_mm_add_epi32(xc0, _mm_shuffle_epi32(xc0, _MM_SHUFFLE(2, 3, 0, 1)));
						pred+=_mm_extract_epi32(xc0, 0);
#endif
#else
#define PRED(EXPR) preds[j]=EXPR; pred+=coeffs[kc][j]*preds[j]; ++j;
						j=0;
						PREDLIST;
#undef  PRED
#endif
						upred2=pred>>(L1SH-2);
						pred>>=L1SH;
					}
				}
				upred2+=offset0&~3;
				CLAMP2(upred2, 0, 1023);
				pred0=(int32_t)pred;
				vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);

				pred+=offset;
#ifdef USE_TABLES
				pred=clampptr[pred];
#else
#ifdef UNSIGNED_PIXEL
				CLAMP2(pred, 0, 255);
#else
				CLAMP2(pred, -128, 127);
#endif
#endif
#ifdef USE_TABLES
				uint8_t *ctxptr=ctxtable[eW<_countof(ctxtable)-1?eW:_countof(ctxtable)-1];
				ctx=ctxptr[0];
				nbypass=ctxptr[1];
#else
				nbypass=31-_lzcnt_u32(eW*eW+2);
				ctx=nbypass;
				nbypass>>=1;
				nbypass-=GRBITS;
				CLAMP2(nbypass, 0, 7);
				if(ctx>NCTX-1)
					ctx=NCTX-1;
#endif
#if 1
				(void)NNNN	;
				(void)NNN	;
				(void)NN	;
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
				(void)eNE	;
				(void)eNEE	;
				(void)eNEEE	;
				(void)eW	;
				(void)xNW	;
				(void)xNE	;
				(void)xNEE	;
				(void)xNEEE	;
				(void)xW	;
				(void)yNW	;
				(void)yN	;
				(void)yNE	;
				(void)yW	;
#endif
#ifdef NEW_SYSTEM
				int32_t p1, p1u[GRLIMIT+1];
				int nonzeroflag=0, signflag=0, sym=0;
				error=0;
				(void)p1;
				(void)p1u;
				if(fwd)
				{
					error=(int8_t)(yuv[kc]-pred);
					nonzeroflag=error!=0;
				}
				statsptr=stats[kc][upred2][ctx]+2;
				p1=codebit(ac, statsptr-2, &nonzeroflag, PROBSHIFT, fwd);
				statsptr[-2]+=(int32_t)((nonzeroflag<<PROBBITS_STORE)-statsptr[-2])>>7;
			//	statsptr[-2]+=(int32_t)((nonzeroflag<<PROBBITS_USE)-p1);
				if(nonzeroflag)
				{
					int32_t *statsptr2=0;

					nzeros=0;
					if(fwd)
					{
						signflag=error>>31&1;
						sym=abs(error);
						nzeros=sym>>nbypass;
					}
					p1=codebit(ac, statsptr-1, &signflag, PROBSHIFT, fwd);
					statsptr[-1]+=(int32_t)((signflag<<PROBBITS_STORE)-statsptr[-1])>>12;
				//	statsptr[-1]+=(int32_t)((signflag<<PROBBITS_USE)-p1);
					for(tidx=0;;)
					{
						bit=tidx>=nzeros;
						p1u[tidx]=codebit(ac, statsptr+tidx, &bit, PROBSHIFT, fwd);
						if(bit)
							break;
						++tidx;
						if(tidx>=GRLIMIT)
							break;
					}
					{
						//static const uint8_t shifttable[GRLIMIT]=
						//{//	0	1	2	3	4	5	6	7	8	9
						//	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,
						//	5,	5,	5,	5,	5,	5,	5,	5,	4,
						//};
						//static const int16_t gaintable[GRLIMIT]=//<256!
						//{
						//	8,	15,	20,	36,	52,	64,	72,	80,	96,	128,
						//	192,	192,	192,	192,	192,	192,	192,	192,
						//};
						//int sh=shifttable[tidx];
						for(int k=0;k<tidx+(tidx<GRLIMIT);++k)
							statsptr[k]+=(int32_t)(((k>=tidx)<<PROBBITS_STORE)-statsptr[k])>>7;
						//	statsptr[k]+=(int32_t)(((k>=tidx)<<PROBBITS_USE)-p1u[k])*gaintable[k]>>sh;
						//	statsptr[k]+=(int32_t)(((k>=tidx)<<PROBBITS_STORE)-statsptr[k])*gaintable[k]>>sh;
					}
					nbypass0=nbypass;
					grflag=tidx>=GRLIMIT;
					if(grflag)
					{
						if(fwd)
							sym-=(GRLIMIT-1)<<nbypass;
						statsptr2=stats3[kc][nbypass];
						tidx=1;
						nbypass=8;
					}
					else
					{
						statsptr2=stats2[kc][upred2];
						tidx=(256>>nbypass)+tidx;//bit coding:  tidx=2*tidx+bit  tidx=0b1XX
					}
					for(int32_t kb=nbypass-1;kb>=0;--kb)
					{
						bit=sym>>kb&1;
						p1=codebit(ac, statsptr2+tidx, &bit, PROBSHIFT, fwd);
						statsptr2[tidx]+=(int32_t)((bit<<PROBBITS_STORE)-statsptr2[tidx])>>7;
					//	statsptr2[tidx]+=(int32_t)((bit<<PROBBITS_USE)-p1+(1<<1>>1))>>1;
						tidx=2*tidx+bit;
					}
					tidx-=256;
#ifdef _MSC_VER
					if(fwd&&tidx!=sym)
						CRASH("");
#endif
					if(!fwd)
					{
						if(grflag)
							tidx+=(GRLIMIT-1)<<nbypass0;
						error=tidx;
						if(signflag)
							error=-tidx;
					}
				}
				if(!fwd)
				{
					yuv[kc]=(uint8_t)(error+pred);
#ifdef ENABLE_GUIDE
					{
						uint8_t *pval=&g_image[3*(iw*ky+kx)+rct_combinations[bestrct][II_PERM_Y+kc]];
						uint8_t val=*pval;
						uint8_t pixel=yuv[kc];
#ifndef UNSIGNED_PIXEL
						pixel+=128;
#endif
						if(pixel!=val)
							CRASH("guide  YXC %d %d %d", ky, kx, kc);
					}
#endif
				}
#else
				int epred;
#ifdef USE_TABLES
				epred=epredptr[pred];
#else
#ifdef UNSIGNED_PIXEL
				if(lossy)
					epred=pred>=128?255-pred:pred;
				else
					epred=128-abs(pred-128);
#else
				if(lossy)
					epred=pred>=0?127-pred:pred+128;
				else
					epred=128-abs(pred);
#endif
#endif
				if(fwd)
				{
					if(lossy)
					{
						int pixel;
						error=yuv[kc]-pred;
						epred=epred*invdist>>16;
						error=(error*invdist>>16)-(error>>31);
						pixel=error*dist+pred;
#ifdef UNSIGNED_PIXEL
#ifdef USE_TABLES
						pixel=clampptr[pixel];
#else
						CLAMP2(pixel, 0, 255);
#endif
						yuv[kc]=pixel;
#else
						CLAMP2(pixel, -128, 127);
						yuv[kc]=pixel;
						pixel+=128;
#endif
#ifdef ENABLE_GUIDE
						{
							uint8_t *pval=&g_image[3*(iw*ky+kx)+rct_combinations[bestrct][II_PERM_Y+kc]];
							uint8_t val=*pval;
							int diff=(int)val-pixel;
							g_sqe[kc]+=diff*diff;
							*pval=pixel;
						}
#endif
						{
							int negmask=error>>31;
							int abserr=(error^negmask)-negmask;
							error=error<<1^negmask;
							if(epred<abserr)
								error=epred+abserr;
						}
					}
					else
					{
#ifdef USE_TABLES
						error=yuv[kc]-(int32_t)pred;
						uint16_t *p=errorptr[error];
						int e2=(int8_t)error;
						error=p[0];
						if(epred<p[1])
							error=epred+p[1];
						if(error==256)
							error=errorptr[e2][0];
#else
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
#endif
					}
					nzeros=error>>nbypass;
				}
				else
					error=0;
//#ifdef UNSIGNED_PIXEL
//				upred=pred;
//#else
//				upred=(uint8_t)(pred+128);
//#endif
				statsptr=stats1[kc][upred2][ctx];
				//histptr=hist1[kc][upred2][ctx];
				tidx=0;
#ifdef ESTIMATE_BITSIZE
				ekc=kc;
#endif
			//	static const shift_unary[GRLIMIT]=
			//	{//	0	1	2	3	4	5	6	7	8	9
			//		9,	8,	8,	7,	7,	7,	7,	7,	7,	7,
			//		7,	7,	7,	7,	7,	7,	7,	7,
			//	};
				//int bits[GRLIMIT]={0};
				
				//if(ky==0&&kx==1044&&kc==2)//
				//	printf("");

				do
				{
					bit=tidx>=nzeros;
#ifdef ESTIMATE_BITSIZE
					eidx=tidx;
#endif
					codebit(ac, statsptr+tidx, &bit, fwd);
					//bits[tidx]=bit;
#ifdef PRINTBITS
					if(fwd&&(uint32_t)(idx-10000)<1000)printf("%c", '0'+bit);//
#endif
//#ifdef _MSC_VER
//					++unary_count;
//#endif
					if(bit)
						break;
					++tidx;
				}while(tidx<GRLIMIT);
#ifndef USE_COUNTERS
				{
					static const uint8_t shifttable[GRLIMIT+1]=
					{//	0	1	2	3	4	5	6	7	8	9
						7+4,	7+4,	7+4,	7+4,	7+4,	7+4,	7+4,	7+4,	7+4,	7+4,
						7+4,	7+4,	7+4,	7+4,	7+4,	7+4,	7+4,	7+4,	4+4,
					};
					static const int16_t gaintable[GRLIMIT]=//<256!
					{
						8,	15,	20,	36,	52,	64,	72,	80,	96,	128,
						192,	192,	192,	192,	192,	192,	192,	192,
					};
					int sh=shifttable[tidx];
					for(int k=0;k<tidx+(tidx<GRLIMIT);++k)
					//	statsptr[k]+=(int32_t)(((k>=tidx)<<PROBBITS_STORE)-statsptr[k])*gaintable[kc*(GRLIMIT+1)+k]>>sh;
						statsptr[k]+=(int32_t)(((k>=tidx)<<PROBBITS_STORE)-statsptr[k])*gaintable[k]>>sh;
				}
#endif
				//if(tidx<GRLIMIT)
				//{
				//	for(int k=0;k<tidx+(tidx<GRLIMIT);++k)
				//	{
				//		//int sh=shift_unary[k]+5;
				//		//statsptr[k]+=(int32_t)(((k>=tidx)<<PROBBITS_STORE)-statsptr[k])*(k-tidx+50)>>sh;
				//
				//		//statsptr[k]+=(int32_t)(((k>=tidx)<<PROBBITS_STORE)-statsptr[k])>>shift_unary[k];
				//		//statsptr[k]+=(int32_t)(((k>=tidx)<<PROBBITS_STORE)-statsptr[k])>>7;
				//		
				//		//statsptr[k]+=(int32_t)((bits[k]<<PROBBITS_STORE)-statsptr[k])*(80-tidx)>>(7+6);
				//		statsptr[k]+=(int32_t)(((k>=tidx)<<PROBBITS_STORE)-statsptr[k])>>7;
				//
				//		//if(statsptr[k]<1||statsptr[k]>(1<<PROBBITS_STORE)-1)
				//		//	CRASH("");
				//	}
				//}
				//else
				//{
				//	for(int k=0;k<tidx+(tidx<GRLIMIT);++k)
				//		statsptr[k]+=(int32_t)(((k>=tidx)<<PROBBITS_STORE)-statsptr[k])>>1;
				//}
				nbypass0=nbypass;
				grflag=tidx==GRLIMIT;
				if(grflag)
				{
					error-=(GRLIMIT-1)<<nbypass;
					statsptr=stats3[kc][nbypass];
					//histptr=hist3[kc][nbypass];
					tidx=1;
					nbypass=8;
				}
				else
				{
					statsptr=stats2[kc][upred2];
					//histptr=hist2[kc][upred2];
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
#endif
						codebit(ac, statsptr+tidx, &bit, fwd);
#ifndef USE_COUNTERS
						statsptr[tidx]+=(int32_t)((bit<<PROBBITS_STORE)-statsptr[tidx])>>7;
#endif
#ifdef PRINTBITS
						if(fwd&&(uint32_t)(idx-10000)<1000)printf("%c", '0'+bit);//
#endif
						tidx=2*tidx+bit;
//#ifdef _MSC_VER
//						++binary_count;
//#endif
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
					if(lossy)
					{
						int pixel, negmask, sym, e2;

						epred=epred*invdist>>16;
#ifdef UNSIGNED_PIXEL
						negmask=((int32_t)pred-128)>>31;
#else
						negmask=(int32_t)pred>>31;
#endif
						sym=error;
						e2=epred-sym;
						error=sym>>1^-(sym&1);
						e2=(e2^negmask)-negmask;
						if((epred<<1)<sym)
							error=e2;

						pixel=error*dist+(int32_t)pred;
#ifdef UNSIGNED_PIXEL
						CLAMP2(pixel, 0, 255);
						yuv[kc]=pixel;
#else
						CLAMP2(pixel, -128, 127);
						yuv[kc]=pixel;
						pixel+=128;
#endif
					}
					else
					{
#ifdef UNSIGNED_PIXEL
						//if(2*(pred-128)+error==256)
						if(2*pred+error==512)
#else
						if(2*pred+error==256)
#endif
						{
							error=error>>1^-(error&1);
#ifdef UNSIGNED_PIXEL
							yuv[kc]=(uint8_t)(error+pred);
#else
							yuv[kc]=(int8_t)(error+pred);
#endif
						}
						else
						{
#ifdef UNSIGNED_PIXEL
							int negmask=((int32_t)pred-128)>>31;
#else
							int negmask=(int32_t)pred>>31;
#endif
							int sym=error;
							int e2=epred-sym;
							error=sym>>1^-(sym&1);
							e2=(e2^negmask)-negmask;
							if((epred<<1)<sym)
								error=e2;
							yuv[kc]=error+(int32_t)pred;
						}
					}
#ifdef ENABLE_GUIDE
					{
						uint8_t *pval=&g_image[3*(iw*ky+kx)+rct_combinations[bestrct][II_PERM_Y+kc]];
						uint8_t val=*pval;
						uint8_t pixel=yuv[kc];
#ifndef UNSIGNED_PIXEL
						pixel+=128;
#endif
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

					k=0;
					if(lossy)
					{
#define PRED(EXPR) lcoeffs[kc][k]+=e*preds[k]; ++k;
						//currw[L1NPREDS_LOSSY]+=e;
						PREDLIST_LOSSY;
#undef  PRED
					}
					else
					{
						bias[kc]+=e<<8;
#define PRED(EXPR) coeffs[kc][k]+=e*preds[k]; ++k;
						PREDLIST;
#undef  PRED
					}

					error=yuv[kc]-(int32_t)pred;
					if(lossy)
						error=abs(error);
					else
						error=error<<1^error>>31;
					rows[0][0]=curr;
					rows[0][1]=curr-W;//gx
					rows[0][2]=curr-N;//gy

					rows[0][3]=(eW+(eW<eNE?eW:eNE)+(error<<GRBITS)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
					offset0=kc ? cv0*yuv[0]+cv1*yuv[1] : cu0*yuv[0];
					offset=offset0>>2;
				}
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
#ifdef PRINTBITS
				++idx;
				(void)idx;
#endif
			}
			if(!fwd)
			{
#if defined UNSIGNED_PIXEL
				imptr[yidx]=yuv[0];
				imptr[uidx]=yuv[1];
				imptr[vidx]=yuv[2];
#else
				imptr[yidx]=yuv[0]+128;
				imptr[uidx]=yuv[1]+128;
				imptr[vidx]=yuv[2]+128;
#endif
#ifdef ENABLE_GUIDE
				guide_check(image, kx, ky);
#endif
			}
#ifdef ZIPF_VIEW
			if(fwd)
			{
				zsize*=255./24;
				if(zsize>255)
					zsize=255;
				*zptr++=(uint8_t)zsize;
				zsize=0;
			}
#endif
			//printf("XY %5d %5d\r", kx, ky);
		}
	}
	free(pixels);
#ifdef ZIPF_VIEW
	if(fwd)
	{
		const char fn[]="c12zipf-20251203_0302am.pgm";
		FILE *fdst=fopen(fn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing\n", fn);
			return;
		}
		fprintf(fdst, "P5\n%d %d\n255\n", iw, ih);
		fwrite(zimage, 1, res, fdst);
		fclose(fdst);
		printf("Saved Zipf file \"%s\"\n", fn);
		free(zimage);
	}
#endif
}
int c12_codec(int argc, char **argv)
{
	static const uint16_t tag='1'|'2'<<8;

	const char *srcfn=0, *dstfn=0;
	int dist=0;
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
	int32_t gaintable[3*(GRLIMIT+1)]={0};

#if 0
#define NBITS 4
#define HALF (1<<NBITS>>1)
	{
		printf("pred\\target  naive modular arithmetic sign packing\n");
		printf("\t");
		for(int k=0;k<2*HALF;++k)
			printf(" %4d", k-HALF);
		printf("\n\n");
		for(int kp=-HALF;kp<HALF;++kp)
		{
			printf(" %4d\t", kp);
			for(int kt=-HALF;kt<HALF;++kt)
			{
				int e=kt-kp, e0;

				e0=e<<(32-NBITS)>>(32-NBITS);
				e0=e0<<1^e0>>31;
				printf(" %4d", e0);
			}
			printf("\n");
		}
		printf("\n");

		printf("pred\\target  CALIC sign deduction\n");
		printf("\t");
		for(int k=0;k<2*HALF;++k)
			printf(" %4d", k-HALF);
		printf("\n\n");
		for(int kp=-HALF;kp<HALF;++kp)
		{
			printf(" %4d\t", kp);
			for(int kt=-HALF;kt<HALF;++kt)
			{
				int e=kt-kp, e1;

				if(kt==-HALF&&kp>0)
				{
					e1=e<<(32-NBITS)>>(32-NBITS);
					e1=e1<<1^e1>>31;
				}
				else
				{
					int upred=HALF-abs(kp);
					int negmask=e>>31;
					int abse=(e^negmask)-negmask;
					e1=e<<1^negmask;
					if(upred<abse)
						e1=upred+abse;
				}
				printf(" %4d", e1);

				//deduce kt from e1 and kp
				{
					int kt2=e1>>1^-(e1&1);
					kt2+=kp;
					kt2=kt2<<(32-NBITS)>>(32-NBITS);
					if(!(kp>0&&kt2==-HALF))
					{
						int upred=HALF-abs(kp);
						int negmask=kp>>31;
						int e2=upred-e1;
						kt2=e1>>1^-(e1&1);
						e2=(e2^negmask)-negmask;
						if(2*upred<e1)
							kt2=e2;
						kt2+=kp;
					}
					if(kt2!=kt)
						CRASH("ERROR");
				}
			}
			printf("\n");
		}
		printf("\n");
		exit(0);
	}
#endif
	if(argc!=3&&argc!=4)
	{
		printf(
			"Usage:  \"%s\"  input  output  [dist]\n"
			"  dist=1 for lossless (default).  Or 3 <= dist <= 17 for lossy.\n"
			, argv[0]
		);
		return 1;
	}
	srcfn=argv[1];
	dstfn=argv[2];
	dist=argc<4?1:atoi(argv[3]);
	if(dist!=1)
		CLAMP2(dist, 3, 17);
	
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
			fread(&iw, 1, 4, fsrc);
			fread(&ih, 1, 4, fsrc);
			bestrct=0;
			fread(&bestrct, 1, 1, fsrc);
			dist=1;
			fread(&dist, 1, 1, fsrc);
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
#if 0
		{
			int yidx=rct_combinations[bestrct][II_PERM_Y];
			int uidx=rct_combinations[bestrct][II_PERM_U];
			int vidx=rct_combinations[bestrct][II_PERM_V];
			int cu0=rct_combinations[bestrct][II_COEFF_U_SUB_Y];
			int cv0=rct_combinations[bestrct][II_COEFF_V_SUB_Y];
			int cv1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
			int psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL0]);
			int16_t *pixels=(int16_t*)malloc(psize);

			if(!pixels)
			{
				CRASH("Alloc error\n");
				return 1;
			}
			memset(pixels, 0, psize);
			imptr=image;
			for(int ky=0;ky<ih;++ky)
			{
				int16_t *rows[]=
				{
					pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL0,
					pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL0,
					pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL0,
					pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL0,
				};
				for(int kx=0;kx<iw;++kx, imptr+=3)
				{
					int y, u, v, py, pu, pv, uoffset, voffset, nby, nbu, nbv;

					y=imptr[yidx];
					u=imptr[uidx];
					v=imptr[vidx];
					py=(rows[1][(0+0*NCH)*NROWS*NVAL0]+rows[0][(0-1*NCH)*NROWS*NVAL0])>>1;
					pu=(rows[1][(1+0*NCH)*NROWS*NVAL0]+rows[0][(1-1*NCH)*NROWS*NVAL0])>>1;
					pv=(rows[1][(2+0*NCH)*NROWS*NVAL0]+rows[0][(2-1*NCH)*NROWS*NVAL0])>>1;
					uoffset=cu0*y>>2;
					voffset=(cv0*y+cv1*v)>>2;
					pu+=uoffset; CLAMP2(pu, 0, 255);
					pv+=voffset; CLAMP2(pv, 0, 255);
					rows[0][(0+0*NCH)*NROWS*NVAL0]=y;
					rows[0][(1+0*NCH)*NROWS*NVAL0]=u-uoffset;
					rows[0][(2+0*NCH)*NROWS*NVAL0]=v-voffset;
					y-=py;
					u-=pu;
					v-=pv;
					y=(int8_t)y;
					u=(int8_t)u;
					v=(int8_t)v;
					y=y<<1^y>>31;
					u=u<<1^u>>31;
					v=v<<1^v>>31;
					nby=31-_lzcnt_u32((rows[0][1+(0-1*NCH)*NROWS*NVAL0]>>GRBITS)+1);
					nbu=31-_lzcnt_u32((rows[0][1+(1-1*NCH)*NROWS*NVAL0]>>GRBITS)+1);
					nbv=31-_lzcnt_u32((rows[0][1+(2-1*NCH)*NROWS*NVAL0]>>GRBITS)+1);
					if(nby>7)nby=7;
					if(nbu>7)nbu=7;
					if(nbv>7)nbv=7;
					rows[0][1+(0+0*NCH)*NROWS*NVAL0]=(2*rows[0][1+(0-1*NCH)*NROWS*NVAL0]+(y<<GRBITS)+rows[1][1+(0+3*NCH)*NROWS*NVAL0])>>2;
					rows[0][1+(1+0*NCH)*NROWS*NVAL0]=(2*rows[0][1+(1-1*NCH)*NROWS*NVAL0]+(u<<GRBITS)+rows[1][1+(1+3*NCH)*NROWS*NVAL0])>>2;
					rows[0][1+(2+0*NCH)*NROWS*NVAL0]=(2*rows[0][1+(2-1*NCH)*NROWS*NVAL0]+(v<<GRBITS)+rows[1][1+(2+3*NCH)*NROWS*NVAL0])>>2;
					y>>=nby;
					u>>=nbu;
					v>>=nbv;
					if(y>GRLIMIT)y=GRLIMIT;
					if(u>GRLIMIT)u=GRLIMIT;
					if(v>GRLIMIT)v=GRLIMIT;
					++gaintable[0*(GRLIMIT+1)+y];
					++gaintable[1*(GRLIMIT+1)+u];
					++gaintable[2*(GRLIMIT+1)+v];
					rows[0]+=NCH*NROWS*NVAL0;
					rows[1]+=NCH*NROWS*NVAL0;
					rows[2]+=NCH*NROWS*NVAL0;
					rows[3]+=NCH*NROWS*NVAL0;
				}
			}
			for(int kc=0;kc<3;++kc)
			{
				int32_t vmax=0;
				for(int ks=0;ks<GRLIMIT+1;++ks)
				{
					int f=gaintable[kc*(GRLIMIT+1)+ks]+1;
					f=(f*127>>7)+(iw*ih*1>>7);
					f=(0x7FFFFFFF / (GRLIMIT+1))/f;
					if(vmax<f)
						vmax=f;
					gaintable[kc*(GRLIMIT+1)+ks]=f;
				}
				for(int ks=0;ks<GRLIMIT+1;++ks)
					gaintable[kc*(GRLIMIT+1)+ks]=(uint32_t)(((uint64_t)gaintable[kc*(GRLIMIT+1)+ks]*255)/vmax);
			}
			free(pixels);
		}
#endif
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
			printf("WH %d*%d  %lld B  RCT %d %s  dist %d  \"%s\"\n"
				, iw, ih, usize
				, bestrct, rct_names[bestrct]
				, dist
				, srcfn
			);
#else
			(void)och_names;
			(void)rct_names;
#endif
		}
		streamptr=stream;
		streamend=stream+usize;
#if 0
		for(int kc=0;kc<3;++kc)
		{
			for(int ks=0;ks<GRLIMIT+1;++ks)
				*streamptr++=(uint8_t)gaintable[kc*(GRLIMIT+1)+ks];
		}
#endif
		ac.ptr=streamptr;
		ac.end=streamend;
	}
	else
	{
		streamptr=stream;
		streamend=stream+srcsize;
#if 0
		for(int kc=0;kc<3;++kc)
		{
			for(int ks=0;ks<GRLIMIT+1;++ks)
				gaintable[kc*(GRLIMIT+1)+ks]=*streamptr++;
		}
#endif
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
	if(dist>1)
	{
		if(fwd)
			mainloop(iw, ih, bestrct, dist, image, stream, &ac, gaintable, 1, 1);
		else
			mainloop(iw, ih, bestrct, dist, image, stream, &ac, gaintable, 1, 0);
	}
	else
	{
		if(fwd)
			mainloop(iw, ih, bestrct, dist, image, stream, &ac, gaintable, 0, 1);
		else
			mainloop(iw, ih, bestrct, dist, image, stream, &ac, gaintable, 0, 0);
	}
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
			dstsize+=fwrite(&iw, 1, 4, fdst);
			dstsize+=fwrite(&ih, 1, 4, fdst);
			dstsize+=fwrite(&bestrct, 1, 1, fdst);
			dstsize+=fwrite(&dist, 1, 1, fdst);
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
#ifdef PRINTBITS
		printf("\n");
#endif
//#ifdef _MSC_VER
//		printf("%12.2lf bytes unary  %8.4lf%%\n", unary_count/8., 100.*unary_count/(unary_count+binary_count));
//		printf("%12.2lf bytes binary\n", binary_count/8.);
//		printf("%8.4lf%% GR CR\n", 100.*csize/((unary_count+binary_count)/8.));
//#endif
		//printf("%12.2lf B zeros\n%12.2lf B ones\n", ac.n[0]/8., ac.n[1]/8.);
		//printf("%12.2lf /%12.2lf bytes GR  %12.6lf bit/sym\n", ac.bitidx/8., ac.totalbits/8., ac.bitidx/(3.*iw*ih));
#ifdef ESTIMATE_BITSIZE
		double total[3]={0};
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
				total[2]=total[1]=total[0]=0;
			}
		}
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
#ifdef ENABLE_GUIDE
	if(fwd&&dist>1)
	{
		double rmse[]=
		{
			sqrt((g_sqe[0]+g_sqe[1]+g_sqe[2])/((double)3*iw*ih)),
			sqrt(g_sqe[0]/((double)iw*ih)),
			sqrt(g_sqe[1]/((double)iw*ih)),
			sqrt(g_sqe[2]/((double)iw*ih)),
		};
		double psnr[]=
		{
			20*log10(255/rmse[0]),
			20*log10(255/rmse[1]),
			20*log10(255/rmse[2]),
			20*log10(255/rmse[3]),
		};
		printf("RMSE  PSNR\n");
		printf("T %12.6lf  %12.6lf\n", rmse[0], psnr[0]);
		printf("Y %12.6lf  %12.6lf\n", rmse[1], psnr[1]);
		printf("U %12.6lf  %12.6lf\n", rmse[2], psnr[2]);
		printf("V %12.6lf  %12.6lf\n", rmse[3], psnr[3]);
	}
#endif
#endif
#ifdef ZIPF_VIEW
	exit(0);
#endif
	(void)dstsize;
	(void)csize;
	(void)&time_sec2;
	(void)&squash;
#ifdef PROFILER
	prof_end(prof_ctx);
#endif
	return 0;
}
