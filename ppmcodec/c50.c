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
#if defined _MSC_VER || defined _WIN32
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
	#define FIFOVAL
#endif


	#define USE_FELICS
//	#define USE_W
//	#define USE_CG
//	#define USE_L1


#ifdef USE_L1
#define PREDLIST\
	PRED(W)\
	PRED(N+W-NW)\
	PRED(2*N-NN)\
	PRED(NE)\
	PRED(2*W-WW)\
	PRED(3*(N-NN)+NNN)\
	PRED(3*(W-WW)+WWW)\
	PRED(W+NE-N)\
	PRED(N+NE-NNE)\

#endif

enum
{
#ifdef USE_L1
	SHIFT=18,
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED
#endif

	SBUFSIZE=512*1024,
	PBUFSIZE=SBUFSIZE/3*3,

	GRBITS=3,
	
	NCTX=18,

	XPAD=8,
	NCH=3,
	NROWS=4,
	NVAL=2,

	GRPREC=3,
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

#ifdef FIFOVAL
static ptrdiff_t fifoidx=0, fifocap=0, fifoidx2=0;
static uint32_t *fifoval=0;
static void valfifo_enqueue(uint32_t val)
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
static void valfifo_check(uint32_t val)
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
	int prev[OCH_COUNT]={0};
	for(ptrdiff_t k=0, len=(ptrdiff_t)3*iw*ih;k<len;k+=3)
	{
		int
			r=image[k+0]<<2,
			g=image[k+1]<<2,
			b=image[k+2]<<2,
			rg=r-g,
			gb=g-b,
			br=b-r;
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
	return bestrct;
}


#ifdef ESTIMATE_SIZES
#define BSIZELIST\
	BSIZE(UNARY)\
	BSIZE(BYPASS)\
	BSIZE(TRUNC)\
	BSIZE(BITS)\

enum
{
#define BSIZE(X) SIZE_##X,
	BSIZELIST
#undef  BSIZE
	SIZE_COUNT,
};
static const char *bsize_labels[]=
{
#define BSIZE(X) #X,
	BSIZELIST
#undef  BSIZE
};
static int g_kc=0;
static int64_t bsizes[3][SIZE_COUNT]={0};
#endif
typedef struct _RiceCoder
{
	uint64_t cache;
	int64_t nbits;
	uint8_t *ptr, *end;
} RiceCoder;
AWM_INLINE void rice_init(RiceCoder *ec, uint8_t *start, uint8_t *end, int fwd)
{
	ec->cache=0;
	ec->nbits=64;
	ec->ptr=start;
	ec->end=end;
	if(!fwd)
	{
		ec->cache=*(uint64_t*)ec->ptr;
		ec->ptr+=8;
		ec->nbits=0;
	}
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
	bsizes[g_kc][SIZE_UNARY]+=nzeros+1LL;
	bsizes[g_kc][SIZE_BYPASS]+=nbypass;
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
AWM_INLINE void bit_pack(RiceCoder *ec, int bit)
{
#ifdef ESTIMATE_SIZES
	++bsizes[g_kc][SIZE_BITS];
#endif
	--ec->nbits;
	ec->cache|=(uint64_t)bit<<ec->nbits;
	if(ec->nbits<=0)
	{
		*(uint64_t*)ec->ptr=ec->cache;
		ec->ptr+=8;
		ec->cache=0;
		ec->nbits=64;
	}
}
AWM_INLINE int bit_unpack(RiceCoder *ec)
{
	int bit;

	++ec->nbits;
	bit=(int)(ec->cache>>(64-ec->nbits));
	ec->cache&=0xFFFFFFFFFFFFFFFF>>ec->nbits;
	if(ec->nbits>=64)
	{
		ec->nbits-=64;
		ec->cache=*(uint64_t*)ec->ptr;
		ec->ptr+=8;
	}
	return bit;
}
AWM_INLINE void truncbin_enc(RiceCoder *ec, int nlevels, int val)
{
	/*
	truncaed binary code:	nlevels>=2
	nlevels	6	5	4	3	2

	0	00	00	00	0	0
	1	01	01	01	10	1
	2	100	10	10	11
	3	101	110	11
	4	110	111
	5	111
	*/
	int k=FLOOR_LOG2(nlevels);
	int nunused=(1<<(k+1))-nlevels;
	int bypass, nbypass;
	if(val<nunused)
		bypass=val, nbypass=k;
	else
		bypass=val+nunused, nbypass=k+1;
#ifdef ESTIMATE_SIZES
	bsizes[g_kc][SIZE_TRUNC]+=bypass;
#endif
	
	if(nbypass>=ec->nbits)
	{
		nbypass-=(int)ec->nbits;
		ec->cache|=(uint64_t)bypass>>nbypass;
		bypass&=0x7FFFFFFF>>(31-nbypass);
		*(uint64_t*)ec->ptr=ec->cache;
		ec->ptr+=8;
		ec->cache=0;
		ec->nbits=64;
	}
	if(nbypass)
	{
		ec->nbits-=nbypass;
		ec->cache|=(uint64_t)bypass<<ec->nbits;
	}
}
AWM_INLINE int truncbin_dec(RiceCoder *ec, int nlevels)
{
	int k=FLOOR_LOG2(nlevels);
	int nunused=(1<<(k+1))-nlevels;
	int val=0;
	
	ec->nbits+=k;
	if(ec->nbits>=64)
	{
		ec->nbits-=64;
		val|=(int)(ec->cache<<ec->nbits);
		ec->cache=*(uint64_t*)ec->ptr;
		ec->ptr+=8;
		k=(int)ec->nbits;
	}
	if(k)
	{
		val|=(int)(ec->cache>>(64-ec->nbits));
		ec->cache&=0xFFFFFFFFFFFFFFFF>>ec->nbits;
	}
	if(val>=nunused)
	{
		++ec->nbits;
		val<<=1;
		val|=(int)(ec->cache>>(64-ec->nbits));
		ec->cache&=0xFFFFFFFFFFFFFFFF>>ec->nbits;
		val-=nunused;
		if(ec->nbits>=64)
		{
			ec->nbits-=64;
			ec->cache=*(uint64_t*)ec->ptr;
			ec->ptr+=8;
		}
	}
	return val;
}

int c50_codec(int argc, char **argv)
{
	const uint16_t tag='5'|'0'<<8;

	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0;
	uint64_t c=0;
	int fwd=0, iw=0, ih=0, bestrct=0;
	int64_t usize=0, ccap=0, csize=0;
	int psize=0;
	int16_t *pixels=0;
	uint8_t *image=0, *stream=0, *imptr=0;
	int yidx=0, uidx=0, vidx=0, uc0=0, vc0=0, vc1=0;
#ifdef USE_L1
	int32_t weights[NCH][NPREDS]={0};
#endif
#ifdef USE_FELICS
	int amin[]=
	{
		0,
		-255,
		-255,
	};
	int amax[]=
	{
		255,
		255,
		255,
	};
	int estim=1<<GRPREC;
#endif
	RiceCoder ec;
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
		ccap=(int64_t)6*iw*ih;
	}
	else
	{
		iw=0;
		ih=0;
		fread(&iw, 1, 3, fsrc);
		fread(&ih, 1, 3, fsrc);
		fread(&bestrct, 1, 1, fsrc);
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
	if(fwd)
	{
		fread(image, 1, usize, fsrc);
		guide_save(image, iw, ih);
		bestrct=crct_analysis(image, iw, ih);
	}
	else
	{
		fread(stream, 1, ccap, fsrc);
	}
	fclose(fsrc);

	rice_init(&ec, stream, stream+ccap, fwd);
	yidx=rct_combinations[bestrct][II_PERM_Y];
	uidx=rct_combinations[bestrct][II_PERM_U];
	vidx=rct_combinations[bestrct][II_PERM_V];
	uc0=rct_combinations[bestrct][II_COEFF_U_SUB_Y];
	vc0=rct_combinations[bestrct][II_COEFF_V_SUB_Y];
	vc1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
#ifdef USE_L1
	for(int k=0;k<NCH*NPREDS;++k)
		((int32_t*)weights)[k]=(1<<SHIFT)/NPREDS;
#endif
	memset(pixels, 0, psize);
	imptr=image;
	for(int ky=0;ky<ih;++ky)
	{
#ifdef USE_FELICS
		int inside=0, below=0, nlevels=0, val=0;
#else
		int error=0, sym=0;
#endif
#ifdef USE_L1
		int estim[NPREDS]={0};
#endif
		int yuv[3]={0};
		int curr=0;
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=0;kx<iw;++kx, imptr+=3)
		{
#ifndef USE_FELICS
			int offset=0;
#endif
			if(fwd)
			{
				yuv[0]=imptr[yidx];
				yuv[1]=imptr[uidx];
				yuv[2]=imptr[vidx];
#ifdef USE_FELICS
				yuv[2]-=(vc0*yuv[0]+vc1*yuv[1])>>2;
				yuv[1]-=uc0*yuv[0]>>2;
				if(vc0+vc1&&uc0)
					yuv[0]+=(yuv[1]+yuv[2])>>2;
#endif
			}
			for(int kc=0;kc<3;++kc)
			{
				int
					NNN	=rows[3][0+0*NCH*NROWS*NVAL],
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL],
					NW	=rows[1][0-1*NCH*NROWS*NVAL],
					N	=rows[1][0+0*NCH*NROWS*NVAL],
					NE	=rows[1][0+1*NCH*NROWS*NVAL],
					NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
					WWW	=rows[0][0-3*NCH*NROWS*NVAL],
					WW	=rows[0][0-2*NCH*NROWS*NVAL],
					W	=rows[0][0-1*NCH*NROWS*NVAL],
					eNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
					eW	=rows[0][1-1*NCH*NROWS*NVAL];
#ifdef ESTIMATE_SIZES
				g_kc=kc;
#endif
#ifdef USE_FELICS
				if(!kx)
					W=NE;
				if(!ky)
					N=WW;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;

				//if(ky==0&&kx==16)//
				//if(ky==0&&kx==5)//
				//if(ky==0&&kx==7&&kc==0)//
				//if(ky==0&&kx==11&&kc==0)//
				//if(ky==0&&kx==184&&kc==2)//
				//if(ky==0&&kx==0&&kc==1)//
				//if(ky==0&&kx==0&&kc==0)//
				//if(ky==0&&kx==0&&kc==1)//
				//if(ky==0&&kx==2&&kc==1)//
				//if(ky==1&&kx==440&&kc==2)//
				//	printf("");

				if(fwd)
				{
					curr=yuv[kc];
					inside=(uint32_t)(curr-vmin)<=(uint32_t)(vmax-vmin);
					bit_pack(&ec, inside);
					if(inside)
					{
						nlevels=vmax-vmin+1;
						if(nlevels>1)
						{
							val=curr-((vmax+vmin+1)>>1);
							truncbin_enc(&ec, vmax-vmin+1, val<<1^val>>31);
						}
					}
					else
					{
						below=curr<vmin;
						bit_pack(&ec, below);
						if(below)
						{
							nlevels=vmin-amin[kc];
							val=vmin-1-curr;
						}
						else
						{
							nlevels=amax[kc]-vmax;
							val=curr-1-vmax;
						}
						if(nlevels>1)
						{
							int e2=estim>>GRPREC;
							rice_enc(&ec, FLOOR_LOG2(e2+1), val);
							estim+=((val<<GRPREC)-estim)>>GRPREC;
						//	estim+=val-e2;
						}
					}
#ifdef FIFOVAL
					valfifo_enqueue(below<<25^inside<<24^estim<<16^nlevels<<8^curr);
#endif
				}
				else
				{
					inside=bit_unpack(&ec);
					if(inside)
					{
						nlevels=vmax-vmin+1;
						val=0;
						if(nlevels>1)
							val=truncbin_dec(&ec, nlevels);
						curr=(val>>1^-(val&1))+((vmax+vmin+1)>>1);
					}
					else
					{
						below=bit_unpack(&ec);
						if(below)
							nlevels=vmin-amin[kc];
						else
							nlevels=amax[kc]-vmax;
						val=0;
						if(nlevels>1)
						{
							int e2=estim>>GRPREC;
							val=rice_dec(&ec, FLOOR_LOG2(e2+1));
							estim+=((val<<GRPREC)-estim)>>GRPREC;
						//	estim+=val-e2;
						}
						if(below)
							curr=vmin-1-val;
						else
							curr=vmax+1+val;
					}
#ifdef FIFOVAL
					valfifo_check(below<<25^inside<<24^estim<<16^nlevels<<8^curr);
#endif
					yuv[kc]=curr;
				}
#else
				int nbypass=FLOOR_LOG2((eW>>GRBITS)+1);
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
				int pred=1<<SHIFT>>1, j=0, p1, e;
				int vmax=N, vmin=W;

				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
#define PRED(E) estim[j]=E; pred+=weights[kc][j]*estim[j]; ++j;
				j=0;
				PREDLIST
#undef  PRED
				pred>>=SHIFT;
				p1=pred;
				CLAMP2(pred, vmin, vmax);
#endif
				pred+=offset;
				CLAMP2(pred, 0, 255);
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
				curr=yuv[kc]-offset;
#ifdef USE_L1
				e=(curr>p1)-(curr<p1);
#define PRED(...) weights[kc][j]+=e*estim[j]; ++j;
				j=0;
				PREDLIST
#undef  PRED
#endif
				rows[0][1]=(2*eW+(sym<<GRBITS)+eNEEE)>>2;
				offset=(kc?vc0*yuv[0]+vc1*yuv[1]:uc0*yuv[0])>>2;
#endif
				rows[0][0]=curr;
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
#ifdef USE_FELICS
				if(vc0+vc1&&uc0)
					yuv[0]-=(yuv[1]+yuv[2])>>2;
				yuv[1]+=uc0*yuv[0]>>2;
				yuv[2]+=(vc0*yuv[0]+vc1*yuv[1])>>2;
#endif
				imptr[yidx]=yuv[0];
				imptr[uidx]=yuv[1];
				imptr[vidx]=yuv[2];
				guide_check(image, kx, ky);
			}
		}
	}
	free(pixels);
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing", dstfn);
			return 1;
		}
		if(fwd)
		{
			rice_flush(&ec);

			csize=0;
			csize+=fwrite(&tag, 1, 2, fdst);
			csize+=fwrite(&iw, 1, 3, fdst);
			csize+=fwrite(&ih, 1, 3, fdst);
			csize+=fwrite(&bestrct, 1, 1, fdst);
			csize+=fwrite(stream, 1, ec.ptr-stream, fdst);
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
#ifdef ESTIMATE_SIZES
		int64_t btotal=0;
		for(int k=0;k<SIZE_COUNT;++k)
			btotal+=bsizes[0][k]+bsizes[1][k]+bsizes[2][k];
		for(int k=0;k<SIZE_COUNT;++k)
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
