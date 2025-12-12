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


#ifdef _MSC_VER
	#define LOUD
	#define ENABLE_GUIDE
#endif


#define GRBITS 3

#define RUNMIN 4
#define RUNBITS 9
#define RUNMAX (1<<RUNBITS)


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
uint8_t *g_im1=0, *g_im2=0;
double g_sqe[3]={0};
static uint8_t* guide_save(const uint8_t *image, int iw, int ih)
{
	uint8_t *im2=0;
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	im2=(uint8_t*)malloc(size);
	if(!im2)
	{
		CRASH("Alloc error");
		return 0;
	}
	memcpy(im2, image, size);
	return im2;
}
static void guide_check(const uint8_t *image, const uint8_t *im0, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, im0+idx, 3))
	{
		printf("\n\nGuide error  X %5d  Y %5d  0x%02X%02X%02X != 0x%02X%02X%02X\n\n"
			, kx
			, ky
			, image[idx+0]
			, image[idx+1]
			, image[idx+2]
			, im0[idx+0]
			, im0[idx+1]
			, im0[idx+2]
		);
		CRASH("Guide error");
		printf("\n");//trick for old debuggers
	}
}
#endif
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


typedef struct _RiceCoder
{
	uint64_t cache, nbits;
	uint8_t *ptr, *end;
} RiceCoder;
AWM_INLINE void rice_enc_init(RiceCoder *ec, uint8_t *start, uint8_t *end)
{
	ec->cache=0;
	ec->nbits=64;
	ec->ptr=start;
	ec->end=end;
}
AWM_INLINE void rice_dec_init(RiceCoder *ec, uint8_t *start, uint8_t *end)
{
	ec->cache=0;
	ec->nbits=64;
	ec->ptr=start;
	ec->end=end;
}
AWM_INLINE void rice_enc_flush(RiceCoder *ec)
{
	*(uint64_t*)ec->ptr=ec->cache;
	ec->ptr+=8;
}
AWM_INLINE void rice_enc(RiceCoder *ec, int nbypass, int sym)
{
	//buffer: {c,c,c,b,b,a,a,a, f,f,f,e,e,e,d,c}, cache: MSB gg[hhh]000 LSB	nbits is number of ASSIGNED bits
	//written 64-bit words are byte-reversed because the CPU is little-endian

	int nzeros=sym>>nbypass, bypass=sym&((1<<nbypass)-1);
	if(nzeros>=ec->nbits)//fill the rest of cache with zeros, and flush
	{
		nzeros-=(int)ec->nbits;
		*(uint64_t*)ec->ptr=ec->cache;
		ec->ptr+=8;
		ec->cache=0;
		while(nzeros>=64)//just flush zeros
		{
			nzeros-=64;
			*(uint64_t*)ec->ptr=ec->cache;
			ec->ptr+=8;
		}
		ec->nbits=64;
	}
	//now there is room for zeros:  0 <= nzeros < nbits <= 64
	ec->nbits-=nzeros;//emit remaining zeros to cache

	bypass|=1<<nbypass;//append 1 stop bit
	++nbypass;
	if(nbypass>=ec->nbits)//not enough free bits in cache:  fill cache, write to list, and repeat
	{
		nbypass-=(int)ec->nbits;
		ec->cache|=(uint64_t)bypass>>nbypass;
		bypass&=(1<<nbypass)-1;
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

	if(ec->cache)
	{
		uint64_t lz=_lzcnt_u64(ec->cache);
		sym=(int)(lz-ec->nbits);
		ec->nbits=lz;
	}
	else
	{
		sym=64-(int)ec->nbits;
		for(;;)
		{
			ec->cache=*(uint64_t*)ec->ptr;
			ec->ptr+=8;
			if(ec->cache)
			{
				ec->nbits=_lzcnt_u64(ec->cache);
				sym+=(int)ec->nbits;
				break;
			}
			sym+=64;
		}
	}
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
		if(ec->nbits)
		{
			sym|=(int)(ec->cache>>(64-ec->nbits));
			ec->cache&=0xFFFFFFFFFFFFFFFF>>ec->nbits;//nbits=61 -> cache&=7;

			//ec->cache&=(1ULL<<sh)-1;//sh=3 -> cache&=7;
		}
		return sym;
	}
	if(nbypass)
	{
		sym|=(int)(ec->cache>>(64-ec->nbits));
		ec->cache&=0xFFFFFFFFFFFFFFFF>>ec->nbits;
	}
	return sym;
}

typedef struct _Symbol
{
	uint8_t nbypass, sym;
} Symbol;
typedef struct _RLECoder
{
	RiceCoder *ec;
	int32_t symbolmode;

	int32_t qstart, qend, qcount, runstart, runcount;
	Symbol *queue;
} RLECoder;
AWM_INLINE void rle_enc_init(RLECoder *rc, RiceCoder *ec)
{
	int queuesize=sizeof(Symbol[RUNMAX]);
	memset(rc, 0, sizeof(*rc));
	rc->ec=ec;
	rc->symbolmode=1;
	rc->queue=(Symbol*)malloc(queuesize);
	if(!rc->queue)
	{
		CRASH("Alloc error");
		return;
	}
	memset(rc->queue, 0, queuesize);
}
AWM_INLINE void rle_enc_finish(RLECoder *rc)
{
	if(rc->qcount)
	{
		rice_enc(rc->ec, 0, rc->symbolmode);
		rice_enc(rc->ec, RUNBITS, rc->qcount);
		if(rc->symbolmode)
		{
			while(rc->qstart!=rc->qend)
			{
				Symbol *sym=rc->queue+rc->qstart;
				rice_enc(rc->ec, sym->nbypass, sym->sym);
				rc->qstart=(rc->qstart+1)%RUNMAX;
			}
		}
	}
	rice_enc_flush(rc->ec);
	free(rc->queue);
}
AWM_INLINE void rle_enc(RLECoder *rc, int nbypass, int sym)
{
	if(rc->qcount>=RUNMAX)
	{
		rice_enc(rc->ec, 0, rc->symbolmode);
		rice_enc(rc->ec, RUNBITS, rc->qcount-1);
		if(rc->symbolmode)
		{
			while(rc->qcount)
			{
				Symbol *sym=rc->queue+rc->qstart;
				rice_enc(rc->ec, sym->nbypass, sym->sym);
				rc->qstart=(rc->qstart+1)%RUNMAX;
				--rc->qcount;
			}
		}
		rc->qstart=rc->qend;
		rc->qcount=0;
		rc->runstart=rc->qend;
		rc->runcount=0;
		rc->symbolmode=1;
	}
	Symbol *ptr=rc->queue+rc->qend;

	//enqueue symbol
	ptr->nbypass=nbypass;
	ptr->sym=sym;
	rc->qend=(rc->qend+1)%RUNMAX;
	++rc->qcount;

	if(sym)
	{
		if(!rc->symbolmode)//end of run
		{
			int end0=(rc->qend-1+RUNMAX)%RUNMAX;
			int count=(rc->qend-1-rc->qstart+RUNMAX)%RUNMAX;
			if(count>=RUNMIN)
			{
				rice_enc(rc->ec, 0, 0);
				rice_enc(rc->ec, RUNBITS, count-1);
				rc->qstart=end0;
				rc->qcount=1;
			}
			rc->symbolmode=1;
		}
		rc->runstart=rc->qend;
		rc->runcount=0;
	}
	else
	{
		++rc->runcount;
		if(rc->symbolmode)
		{
			if(rc->runcount>=RUNMIN)//run confirmed: backtrack
			{
				int count=(rc->runstart-rc->qstart+RUNMAX)%RUNMAX;
				if(count)
				{
					rice_enc(rc->ec, 0, 1);
					rice_enc(rc->ec, RUNBITS, count-1);
					while(rc->qstart!=rc->runstart)
					{
						Symbol *sym=rc->queue+rc->qstart;
						rice_enc(rc->ec, sym->nbypass, sym->sym);
						rc->qstart=(rc->qstart+1)%RUNMAX;
					}
					rc->qcount-=count;
				}
				rc->symbolmode=0;
			}
		}
	}
}
AWM_INLINE void rle_dec_init(RLECoder *rc, RiceCoder *ec)
{
	memset(rc, 0, sizeof(*rc));
	rc->ec=ec;
}
AWM_INLINE int rle_dec(RLECoder *rc, int nbypass)
{
	int sym=0;
	if(rc->qcount<=0)
	{
		rc->symbolmode=rice_dec(rc->ec, 0);
		rc->qcount=rice_dec(rc->ec, RUNBITS)+1;
	}
	if(rc->symbolmode)
		sym=rice_dec(rc->ec, nbypass);
	--rc->qcount;
	return sym;
}

static void predict(uint8_t *image, int iw, int ih, int16_t *pixels, int psize, int rct, int fwd)
{
	uint8_t *imptr=image;
	int yidx=rct_combinations[rct][II_PERM_Y];
	int uidx=rct_combinations[rct][II_PERM_U];
	int vidx=rct_combinations[rct][II_PERM_V];
	int umask=-(rct_combinations[rct][II_COEFF_U_SUB_Y]!=0);
	int vc0=rct_combinations[rct][II_COEFF_V_SUB_Y];
	int vc1=rct_combinations[rct][II_COEFF_V_SUB_U];
#ifdef ENABLE_GUIDE
	int perm[]={yidx, uidx, vidx};
#endif

	memset(pixels, 0, psize);
	for(int ky=0;ky<ih;++ky)
	{
		ALIGN(32) int16_t *rows[]=
		{
			pixels+(8*3*2+(ky-0LL+2)%2)*1,//base + (XPAD*NCH*NROWS + (CY-NY+NROWS)%NROWS)*NVAL
			pixels+(8*3*2+(ky-1LL+2)%2)*1,
		};
		uint8_t yuv[3]={0};
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
					NW	=rows[1][-1*3*2*1],//NCH*NROWS*NVAL
					N	=rows[1][+0*3*2*1],
					W	=rows[0][-1*3*2*1];
				int error, epred, negmask, sym;

				int pred=abs(N-NW)>abs(W-NW)?N:W;
				//int vmax=N, vmin=W, pred=N+W-NW;
				//if(N<W)vmin=N, vmax=W;
				//CLAMP2(pred, vmin, vmax);

				pred+=offset;
				CLAMP2(pred, 0, 255);
				
				epred=128-abs(pred-128);
				if(fwd)
				{
					int e0, abserr;

					error=yuv[kc]-pred;
					e0=(int8_t)error;
					negmask=error>>31;
					abserr=(error^negmask)-negmask;
					sym=error<<1^negmask;
					if(epred<abserr)
						sym=epred+abserr;
					if(sym==256)
						sym=e0<<1^e0>>31;
					imptr[kc]=sym;
				}
				else
				{
					sym=imptr[kc];
					
					error=sym>>1^-(sym&1);
					yuv[kc]=(int8_t)(error+pred);
					if(2*pred+sym!=512)
					{
						negmask=(pred-128)>>31;
						int e2=epred-sym;
						error=sym>>1^-(sym&1);
						e2=(e2^negmask)-negmask;
						if((epred<<1)<sym)
							error=e2;
						yuv[kc]=error+pred;
					}
#ifdef ENABLE_GUIDE
					if(g_im1[imptr-image+perm[kc]]!=yuv[kc])
					{
						printf("\n\nGuide error  X %5d  Y %5d C%d  0x%02X != 0x%02X\n\n"
							, kx
							, ky
							, kc
							, yuv[kc]
							, g_im1[3*(iw*ky+kx)+perm[kc]]
						);
						CRASH("Guide error");
					}
#endif
				}
				rows[0][0]=yuv[kc]-offset;
				offset=kc?(vc0*yuv[0]+vc1*yuv[1])>>2:yuv[0]&umask;
				rows[0]+=2*1;//NROWS*NVAL
				rows[1]+=2*1;
			}
			if(!fwd)
			{
				imptr[yidx]=yuv[0];
				imptr[uidx]=yuv[1];
				imptr[vidx]=yuv[2];
#ifdef ENABLE_GUIDE
				guide_check(image, g_im1, kx, ky);
#endif
			}
		}
	}
}
int c41_codec(int argc, char **argv)
{
	if(argc!=3)
	{
		printf(
			"Usage:  \"%s\"  input  output    To encode/decode.\n"
			"Built on %s %s\n"
			, argv[0]
			, __DATE__, __TIME__
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];
#ifdef LOUD
	double t=time_sec();
#endif
	int fwd=0, iw=0, ih=0;
	int bestrct=0;
	ptrdiff_t usize=0, csize=0, headersize=0, cap=0;
	uint8_t *buf=0, *image=0, *streamstart=0, *streamend=0;
	int padw=0, psize=0;
	int16_t *pixels=0;
	RiceCoder ec;
	RLECoder rc;
#ifdef ENABLE_GUIDE
	static uint8_t *im0=0;
#endif

#ifdef LOUD
	t=time_sec();
#endif
	{
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			CRASH("Cannot open \"%s\"", srcfn);
			return 1;
		}
		int64_t c=0;
		fread(&c, 1, 2, fsrc);
		fwd=c==('P'|'6'<<8);
		if(!fwd&&c!=('4'|'1'<<8))
		{
			CRASH("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
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

			//nread=fscanf(fsrc, "%d %d", &iw, &ih);
			//if(nread!=2)
			//{
			//	CRASH("Unsupported PPM file");
			//	return 1;
			//}
			//nread=fscanf(fsrc, "%d", &vmax);
			//if(nread!=1||vmax!=255)
			//{
			//	CRASH("Unsupported PPM file");
			//	return 1;
			//}
			//c=fgetc(fsrc);
			//if(c!='\n')
			//{
			//	CRASH("Invalid PPM file");
			//	return 1;
			//}
		}
		else
		{
			iw=0;
			ih=0;
			fread(&iw, 1, 3, fsrc);
			fread(&ih, 1, 3, fsrc);
			fread(&bestrct, 1, 1, fsrc);
		}
		if(iw<1||ih<1)
		{
			CRASH("Unsupported source file");
			return 1;
		}
		headersize=ftell(fsrc);
		usize=(ptrdiff_t)3*iw*ih;
		cap=(ptrdiff_t)7*iw*ih;
		buf=(uint8_t*)malloc(cap);
		
		padw=iw+8*2;
		psize=padw*(int)sizeof(int16_t[3*2*1]);//int16_t[iw+2*XPAD][NCH][NROWS][NVAL]

		//psize=padw*(int)sizeof(int16_t[4*3*2]);//4 padded rows * 3 channels * {pixel, error}
		pixels=(int16_t*)malloc(psize);
		if(!buf||!pixels)
		{
			CRASH("Alloc error");
			return 1;
		}
		if(fwd)
		{
			image=buf+cap-usize-sizeof(uint64_t);
			streamstart=buf;
			streamend=buf+cap;

			fread(image, 1, usize, fsrc);//read image
#ifdef ENABLE_GUIDE
			g_im1=guide_save(image, iw, ih);
#endif
			bestrct=crct_analysis(image, iw, ih);

			rice_enc_init(&ec, streamstart, streamend);
			rle_enc_init(&rc, &ec);
		}
		else
		{
			struct stat info={0};
			stat(srcfn, &info);
			csize=info.st_size;

			image=buf;
			streamstart=buf+cap-csize-sizeof(uint64_t);
			streamend=buf+cap;

			fread(streamstart, 1, csize-headersize, fsrc);//read stream
			
			rice_dec_init(&ec, streamstart, streamend);
			rle_dec_init(&rc, &ec);
		}
		fclose(fsrc);
	}
	if(fwd)
	{
		predict(image, iw, ih, pixels, psize, bestrct, 1);
#ifdef ENABLE_GUIDE
		g_im2=guide_save(image, iw, ih);
#endif
	}
	for(int k=0;k<psize/sizeof(int16_t);k+=2)
	{
		pixels[k+0]=0;
		pixels[k+1]=128;
	}
	for(int kc=0;kc<3;++kc)
	{
		uint8_t *imptr=image+kc;
		for(int ky=0;ky<ih;++ky)
		{
			ALIGN(32) int16_t *rows[]=
			{
				pixels+(8*1*2+(ky-0LL+2)%2)*1,//base + (XPAD*NCH*NROWS + (CY-NY+NROWS)%NROWS)*NVAL
				pixels+(8*1*2+(ky-1LL+2)%2)*1,
			};
			for(int kx=0;kx<iw;++kx, imptr+=3)
			{
				int
					eNEE	=rows[1][+2*1*2*1],//NCH*NROWS*NVAL
					eNEEE	=rows[1][+3*1*2*1],
					eW	=rows[0][-1*1*2*1];
				int nbypass=eW>>GRBITS;
				int error;

				nbypass=FLOOR_LOG2(nbypass+1);
				if(fwd)
				{
					error=*imptr;

					rle_enc(&rc, nbypass, error);
				}
				else
				{
					error=rle_dec(&rc, nbypass);
					*imptr=error;
#ifdef ENABLE_GUIDE
					if(g_im2[imptr-image]!=*imptr)
					{
						printf("\n\nGuide error  X %5d  Y %5d C%d  0x%02X != 0x%02X\n\n"
							, kx
							, ky
							, kc
							, *imptr
							, g_im2[3*(iw*ky+kx)+kc]
						);
						CRASH("Guide error");
					}
#endif
				}
				rows[0][0]=(2*eW+(error<<GRBITS)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
				rows[0]+=2*1;//NROWS*NVAL
				rows[1]+=2*1;
			}
		}
	}
	if(!fwd)
	{
#ifdef ENABLE_GUIDE
		free(g_im2);
#endif
		predict(image, iw, ih, pixels, psize, bestrct, 0);
	}
	free(pixels);
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing\n", dstfn);
			free(buf);
			return 1;
		}
		if(fwd)
		{
			rle_enc_finish(&rc);

			csize=0;
			csize+=fwrite("41", 1, 2, fdst);
			csize+=fwrite(&iw, 1, 3, fdst);
			csize+=fwrite(&ih, 1, 3, fdst);
			csize+=fwrite(&bestrct, 1, 1, fdst);
			csize+=fwrite(streamstart, 1, ec.ptr-streamstart, fdst);
		}
		else
		{
			headersize=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(image, 1, usize, fdst);
		}
		fclose(fdst);
	}
	free(buf);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		usize+=headersize;
		printf("\"%s\"  CWH=3*%d*%d  %s\n", srcfn, iw, ih, rct_names[bestrct]);
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
	(void)time_sec;
	(void)och_names;
	(void)rct_names;
	return 0;
}
