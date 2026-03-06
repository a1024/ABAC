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
#include<math.h>
#include<sys/stat.h>
#if defined _MSC_VER || defined _WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#else
#include<time.h>
#endif
#ifdef _MSC_VER
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif


#ifdef _MSC_VER
	#define LOUD
//	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif

	#define USE_MIX4
	#define ENABLE_CRCT
	#define ENABLE_GRCTX


enum
{
	NCTX=18,
	PROBBITS=14,
	BUFSIZE=512*1024,

	XPAD=8,
	NROWS=4,
	NVAL=8,
#ifdef USE_MIX4
	NPREDS=4,
	SHIFT=20,
#endif

	MIXNSHIFT=7,
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
#	define INLINE __forceinline static
#else
#	define	ALIGN(N) __attribute__((aligned(N)))
#	define INLINE __attribute__((always_inline)) inline static
#	ifndef _countof
#		define _countof(A) (sizeof(A)/sizeof(*(A)))
#	endif
#endif
#if defined _M_X64 || defined __x86_64__
#define FLOOR_LOG2(X)\
	(sizeof(X)==8?63-(int32_t)_lzcnt_u64(X):31-_lzcnt_u32((uint32_t)(X)))
#else
INLINE int floor_log2_64(uint64_t n)
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
INLINE int floor_log2_32(uint32_t n)
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
#define FLOOR_LOG2_32x4(X) _mm_sub_epi32(_mm_srli_epi32(_mm_castps_si128(_mm_cvtepi32_ps(X)), 23), _mm_set1_epi32(127))
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
	static int64_t t0=0;
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
static void guide_save(FILE *f, int iw, int ih)
{
	ptrdiff_t idx=0, size=0;
	
	size=(ptrdiff_t)3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(uint8_t*)malloc(size);
	if(!g_image)
	{
		CRASH("Alloc error");
		return;
	}
	idx=ftell(f);
	fread(g_image, 1, size, f);
	fseek(f, (long)idx, SEEK_SET);
}
static void guide_check(uint8_t *rgb, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(rgb, g_image+idx, 3))
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


static uint8_t rdbuf[BUFSIZE+sizeof(uint64_t)], wtbuf[BUFSIZE+sizeof(uint64_t)];
INLINE uint64_t acme_read(uint8_t **pptr, ptrdiff_t size, FILE *f)
{
	uint8_t *ptr=*pptr;
	uint64_t data=*(uint64_t*)ptr;
	ptrdiff_t left;

	/*
	overflow:
	|                    ______left______   ______right_____
	|                   /                \ /                \
	|buf1start ... ... [datastart  buf1end|buf2start  dataend] ...
	|                   \________________    _______________/
	|                                    size
	*/
	
	left=rdbuf+BUFSIZE-ptr;
	ptr+=size;
	if(left<size)
	{
		fread(rdbuf, 1, BUFSIZE, f);
		ptr=(rdbuf+size)-left;
		left<<=3;
		data&=0xFFFFFFFFFFFFFFFF>>(64-left);
		data|=*(uint64_t*)rdbuf<<left;
	}
	*pptr=ptr;
	return data;
}
INLINE void acme_write(uint8_t **pptr, ptrdiff_t size, FILE *f, uint64_t data)
{
	uint8_t *ptr=*pptr;
	ptrdiff_t left;
	
	/*
	overflow:
	|                    ______left______   ______right_____
	|                   /                \ /                \
	|buf1start ... ... [datastart  buf1end|buf2start  dataend] ...
	|                   \________________    _______________/
	|                                    size
	*/
	
	*(uint64_t*)ptr=data;
	left=wtbuf+BUFSIZE-ptr;
	ptr+=size;
	if(left<size)
	{
		fwrite(wtbuf, 1, BUFSIZE, f);
		ptr=(wtbuf+size)-left;
		left<<=3;
		data>>=left;
		*(uint64_t*)wtbuf=data;
	}
	*pptr=ptr;
}


//cRCT
#ifdef ENABLE_CRCT
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
static int crct_analysis(FILE *f, int iw, int ih)
{
	int64_t counters[OCH_COUNT]={0};
	int prev[OCH_COUNT]={0};
	uint8_t *ptr=rdbuf+BUFSIZE;
	long idx=ftell(f);

	for(ptrdiff_t k=0, size=(ptrdiff_t)3*iw*ih;k<size;k+=3)
	{
		int r, g, b, rg, gb, br;

		uint64_t data=acme_read(&ptr, 3, f);
		r=data>> 0&255;
		g=data>> 8&255;
		b=data>>16&255;
		rg=r-g;
		gb=g-b;
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
#define UPDATE(IDXA, A0, IDXB, B0, IDXC, C0)\
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
			OCH_C2X2, (gb-rg)>>1,//g-(r+b)/2 = (g-r + g-b)/2
			OCH_C22X, (br-gb)>>1 //b-(r+g)/2 = (b-r + b-g)/2
		);
#undef  UPDATE
#endif
	}
	fseek(f, idx, SEEK_SET);
	{
		int bestrct=0;
		int64_t minerr=0;
		for(int kt=0;kt<RCT_COUNT;++kt)
		{
			const uint8_t *rct=rct_combinations[kt];
			int64_t currerr=
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
}
#endif


ALIGN(32) static uint16_t stats0[3][NCTX][32], stats1[3][NCTX][16][32], mixinCDFs[16][16];
int c11_codec(int argc, char **argv)
{
	const uint16_t tag='1'|'1'<<8;

	if(argc!=3)
	{
		printf(
			"Usage: \"%s\"  input  output    Encode/decode.\n"
			, argv[0]
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];
#ifdef LOUD
	double t=time_sec();
#endif
	if(!srcfn||!dstfn)
	{
		CRASH("Codec requires both source and destination filenames");
		return 1;
	}
	FILE *fsrc=fopen(srcfn, "rb");
	int c=0;
	size_t nread=fread(&c, 1, 2, fsrc);
	if(nread!=2)
	{
		CRASH("File is empty");
		return 1;
	}
	int fwd=c==('P'|'6'<<8);
	//int bypass=c==('B'|'P'<<8);
	if(!fwd&&c!=tag)
	{
		CRASH("Unsupported source file");
		return 1;
	}
	int iw=0, ih=0;
	int64_t usize=0;
	int bestrct=0;
	if(fwd)
	{
		int temp=fgetc(fsrc);
		if(temp!='\n')
		{
			CRASH("Unsupported PPM file");
			return 1;
		}
		nread=fscanf(fsrc, "%d %d", &iw, &ih);
		if(nread!=2)
		{
			CRASH("Unsupported PPM file");
			return 1;
		}
		int vmax=0;
		nread=fscanf(fsrc, "%d", &vmax);
		if(nread!=1||vmax!=255)
		{
			CRASH("Unsupported PPM file");
			return 1;
		}
		temp=fgetc(fsrc);
		if(temp!='\n')
		{
			CRASH("Unsupported PPM file");
			return 1;
		}
	}
	else
	{
		iw=0;
		ih=0;
		nread=fread(&iw, 1, 3, fsrc);
		nread+=fread(&ih, 1, 3, fsrc);
		nread+=fread(&bestrct, 1, 1, fsrc);
		if(nread!=(size_t)3+3+1)
		{
			CRASH("Unsupported archive");
			return 1;
		}
	}
	usize=(int64_t)3*iw*ih;
	FILE *fdst=fopen(dstfn, "wb");
	uint64_t low=0, range=0xFFFFFFFFFFFF, code=0;
	uint8_t *rdptr=rdbuf+BUFSIZE, *wtptr=wtbuf;
	if(fwd)
	{
		bestrct=crct_analysis(fsrc, iw, ih);
		guide_save(fsrc, iw, ih);

		fwrite(&tag, 1, 2, fdst);
		fwrite(&iw, 1, 3, fdst);
		fwrite(&ih, 1, 3, fdst);
		fwrite(&bestrct, 1, 1, fdst);
	}
	else
	{
		code=acme_read(&rdptr, sizeof(code), fsrc);
		code=code<<32|code>>32;

		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	}
	int psize=(iw+2*XPAD)*(int)sizeof(int32_t[NROWS*NVAL]);
	int32_t *pixels=(int32_t*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	
	for(int k=0;k<16;++k)
		stats0[0][0][k]=k<<(PROBBITS-4);
	for(int k=0;k<16;++k)
		stats0[0][0][k+16]=0x4000;
	for(int k=1;k<NCTX;++k)
		memcpy(stats0[0][k], stats0[0][0], sizeof(uint16_t[32]));
	for(int kc=1;kc<3;++kc)
		memcpy(stats0[kc], stats0[0], sizeof(uint16_t[NCTX][32]));
	for(int kc=0;kc<3;++kc)
	{
		for(int kctx=0;kctx<NCTX;++kctx)
		{
			for(int ks=0;ks<16;++ks)
			{
				for(int k=0;k<16;++k)
					stats1[kc][kctx][ks][k]=k<<(PROBBITS-4);
				for(int k=0;k<16;++k)
					stats1[kc][kctx][ks][k+16]=0x4000;
			}
		}
	}
	for(int ks=0;ks<16;++ks)
	{
		for(int kv=0;kv<16;++kv)
			mixinCDFs[ks][kv]=((0x4000-16)&-(ks<kv))+kv;
	}
	memset(pixels, 0, psize);
	int
		yidx=rct_combinations[bestrct][II_PERM_Y],
		uidx=rct_combinations[bestrct][II_PERM_U],
		vidx=rct_combinations[bestrct][II_PERM_V],
		cu0=rct_combinations[bestrct][II_COEFF_U_SUB_Y],
		cv0=rct_combinations[bestrct][II_COEFF_V_SUB_Y],
		cv1=rct_combinations[bestrct][II_COEFF_V_SUB_U],
		yshift=yidx*8,
		ushift=uidx*8,
		vshift=vidx*8;
	ALIGN(16) int32_t coeffs[NPREDS][4]={0};
	for(int k=0;k<NPREDS*4;++k)
		((uint32_t*)coeffs)[k]=(1<<SHIFT)/NPREDS;
	for(int ky=0, idx=0;ky<ih;++ky)
	{
		ALIGN(32) int32_t *rows[]=
		{
			pixels+(XPAD*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		ALIGN(16) int32_t preds[4]={0}, ctx[4]={0};
		ALIGN(16) int32_t errors[4]={0}, syms[6]={0}, cdfs[6]={0}, freqs[6]={0};
		uint8_t yuv[3]={0};
		int uoffset, voffset;
		for(int kx=0;kx<iw;++kx, idx+=3)
		{
			int32_t
				*NNN	=rows[3]+0*NROWS*NVAL+0,
				*NN	=rows[2]+0*NROWS*NVAL+0,
				*NNE	=rows[2]+1*NROWS*NVAL+0,
				*NW	=rows[1]-1*NROWS*NVAL+0,
				*N	=rows[1]+0*NROWS*NVAL+0,
				*NE	=rows[1]+1*NROWS*NVAL+0,
				*NEEE	=rows[1]+3*NROWS*NVAL+0,
				*WWW	=rows[0]-3*NROWS*NVAL+0,
				*WW	=rows[0]-2*NROWS*NVAL+0,
				*W	=rows[0]-1*NROWS*NVAL+0,
				*curr	=rows[0]+0*NROWS*NVAL+0;
			int32_t
				*eNEE	=rows[1]+2*NROWS*NVAL+4,
				*eNEEE	=rows[1]+3*NROWS*NVAL+4,
				*eW	=rows[0]-1*NROWS*NVAL+4,
				*ecurr	=rows[0]+0*NROWS*NVAL+4;
#if 1
			(void)NNN;
			(void)NN;
			(void)NNE;
			(void)NW;
			(void)N;
			(void)NE;
			(void)NEEE;
			(void)WWW;
			(void)WW;
			(void)W;
#endif
#ifdef USE_MIX4
			/*
			0	N
			1	N+W-NW
			2	2*N-NN
			3	NE
			*/
			__m128i e0=_mm_load_si128((__m128i*)W);
			__m128i e1=_mm_load_si128((__m128i*)N);
			__m128i e2=_mm_load_si128((__m128i*)NN);
			__m128i e3=_mm_load_si128((__m128i*)NE);
			e2=_mm_add_epi32(e2, e2);
			e2=_mm_sub_epi32(e2, e1);
			e1=_mm_add_epi32(e1, e0);
			e1=_mm_sub_epi32(e1, _mm_load_si128((__m128i*)NW));
			__m128i co0=_mm_load_si128((__m128i*)coeffs[0]);
			__m128i co1=_mm_load_si128((__m128i*)coeffs[1]);
			__m128i co2=_mm_load_si128((__m128i*)coeffs[2]);
			__m128i co3=_mm_load_si128((__m128i*)coeffs[3]);
			__m128i t0=_mm_mullo_epi32(co0, e0);
			__m128i t1=_mm_mullo_epi32(co1, e1);
			__m128i t2=_mm_mullo_epi32(co2, e2);
			__m128i t3=_mm_mullo_epi32(co3, e3);
			
			__m128i s0=_mm_load_si128((__m128i*)N);
			__m128i s1=_mm_load_si128((__m128i*)W);
			__m128i vmin=_mm_min_epi32(s0, s1);
			__m128i vmax=_mm_max_epi32(s0, s1);
			s0=_mm_load_si128((__m128i*)NE);
			vmin=_mm_min_epi32(vmin, s0);
			vmax=_mm_max_epi32(vmax, s0);

			__m128i mp=_mm_set1_epi32(1<<SHIFT>>1);
			mp=_mm_add_epi32(mp, t0);
			mp=_mm_add_epi32(mp, t1);
			mp=_mm_add_epi32(mp, t2);
			mp=_mm_add_epi32(mp, t3);
			mp=_mm_srai_epi32(mp, SHIFT);
			mp=_mm_max_epi32(mp, vmin);
			mp=_mm_min_epi32(mp, vmax);
			_mm_store_si128((__m128i*)preds, mp);
#else
			__m128i mNW	=_mm_load_si128((__m128i*)NW);
			__m128i mN	=_mm_load_si128((__m128i*)N);
			__m128i mW	=_mm_load_si128((__m128i*)W);
			__m128i vmin=_mm_min_epi32(mN, mW);
			__m128i vmax=_mm_max_epi32(mN, mW);
			__m128i mp=_mm_sub_epi32(_mm_add_epi32(mN, mW), mNW);
			mp=_mm_max_epi32(mp, vmin);
			mp=_mm_min_epi32(mp, vmax);
			_mm_store_si128((__m128i*)preds, mp);
#endif
#ifdef ENABLE_GRCTX
			__m128i mctx=_mm_load_si128((__m128i*)eW);
			mctx=_mm_mullo_epi32(mctx, mctx);
			mctx=_mm_add_epi32(mctx, _mm_set1_epi32(1));
			mctx=FLOOR_LOG2_32x4(mctx);
			mctx=_mm_min_epi32(mctx, _mm_set1_epi32(NCTX-1));
			_mm_store_si128((__m128i*)ctx, mctx);
#endif
			uint16_t *currstats[6]=
			{
				stats0[0][ctx[0]],
				stats0[1][ctx[1]],
				stats0[2][ctx[2]],
			};
			//if(ky==1&&kx==39)//
			//if(ky==1&&kx==5)//
			//if(ky==1&&kx==743)//
			//if(ky==10&&kx==10)//
			//if(ky==ih/2&&kx==iw/2)//
			//	printf("");
			if(fwd)
			{
				uint64_t data=acme_read(&rdptr, 3, fsrc);
				yuv[0]=data>>yshift&255;
				yuv[1]=data>>ushift&255;
				yuv[2]=data>>vshift&255;
				
				uoffset=cu0*yuv[0]>>2;
				voffset=(cv0*yuv[0]+cv1*yuv[1])>>2;
				preds[1]+=uoffset; CLAMP2(preds[1], 0, 255);
				preds[2]+=voffset; CLAMP2(preds[2], 0, 255);
				errors[0]=yuv[0]-preds[0];
				errors[1]=yuv[1]-preds[1];
				errors[2]=yuv[2]-preds[2];
				curr[0]=yuv[0];
				curr[1]=yuv[1]-uoffset;
				curr[2]=yuv[2]-voffset;

				syms[0]=errors[0]>>4&15;
				syms[1]=errors[1]>>4&15;
				syms[2]=errors[2]>>4&15;
				syms[3]=errors[0]&15;
				syms[4]=errors[1]&15;
				syms[5]=errors[2]&15;
				currstats[3]=stats1[0][ctx[0]][syms[0]];
				currstats[4]=stats1[1][ctx[1]][syms[1]];
				currstats[5]=stats1[2][ctx[2]][syms[2]];
				cdfs[0]=currstats[0][syms[0]];
				cdfs[1]=currstats[1][syms[1]];
				cdfs[2]=currstats[2][syms[2]];
				cdfs[3]=currstats[3][syms[3]];
				cdfs[4]=currstats[4][syms[4]];
				cdfs[5]=currstats[5][syms[5]];
				freqs[0]=currstats[0][syms[0]+1]-cdfs[0];
				freqs[1]=currstats[1][syms[1]+1]-cdfs[1];
				freqs[2]=currstats[2][syms[2]+1]-cdfs[2];
				freqs[3]=currstats[3][syms[3]+1]-cdfs[3];
				freqs[4]=currstats[4][syms[4]+1]-cdfs[4];
				freqs[5]=currstats[5][syms[5]+1]-cdfs[5];
#ifdef FIFOVAL
				fifoval_enqueue(cdfs[0]<<16|freqs[0]);
				fifoval_enqueue(cdfs[1]<<16|freqs[1]);
				fifoval_enqueue(cdfs[2]<<16|freqs[2]);
				fifoval_enqueue(cdfs[3]<<16|freqs[3]);
				fifoval_enqueue(cdfs[4]<<16|freqs[4]);
				fifoval_enqueue(cdfs[5]<<16|freqs[5]);
#endif
				for(int k=0;k<6;++k)
				{
					while(range<=0xFFFF)
					{
						acme_write(&wtptr, sizeof(uint32_t), fdst, (uint32_t)(low>>32));
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					low+=range*cdfs[k]>>PROBBITS;
					range=(range*freqs[k]>>PROBBITS)-1;
				}
			}
			else
			{
				for(int kc=0;kc<3;++kc)
				{
					unsigned cdf, freq;

					while(range<=0xFFFF)
					{
						code=code<<32|(acme_read(&rdptr, sizeof(uint32_t), fsrc)&0xFFFFFFFF);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					__m256i ms0=_mm256_load_si256((__m256i*)currstats[kc]);
					__m256i mc2=_mm256_set1_epi16((uint16_t)(((code-low)<<PROBBITS|((1LL<<PROBBITS)-1))/range));
					int mask=_mm256_movemask_epi8(_mm256_cmpgt_epi16(ms0, mc2));
					syms[kc]=FLOOR_LOG2(~mask)>>1;
					cdf=currstats[kc][syms[kc]];
					freq=currstats[kc][syms[kc]+1]-cdf;
#ifdef FIFOVAL
					fifoval_check(cdf<<16|freq);
#endif
					low+=range*cdf>>PROBBITS;
					range=(range*freq>>PROBBITS)-1;
				}
				currstats[3]=stats1[0][ctx[0]][syms[0]];
				currstats[4]=stats1[1][ctx[1]][syms[1]];
				currstats[5]=stats1[2][ctx[2]][syms[2]];
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
				for(int kc=0;kc<3;++kc)
				{
					unsigned cdf, freq;
					
					while(range<=0xFFFF)
					{
						code=code<<32|(acme_read(&rdptr, sizeof(uint32_t), fsrc)&0xFFFFFFFF);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					__m256i ms0=_mm256_load_si256((__m256i*)currstats[kc+3]);
					__m256i mc2=_mm256_set1_epi16((uint16_t)(((code-low)<<PROBBITS|((1LL<<PROBBITS)-1))/range));
					int mask=_mm256_movemask_epi8(_mm256_cmpgt_epi16(ms0, mc2));
					syms[kc+3]=FLOOR_LOG2(~mask)>>1;
					cdf=currstats[kc+3][syms[kc+3]];
					freq=currstats[kc+3][syms[kc+3]+1]-cdf;
#ifdef FIFOVAL
					fifoval_check(cdf<<16|freq);
#endif
					low+=range*cdf>>PROBBITS;
					range=(range*freq>>PROBBITS)-1;
				}
				yuv[0]=errors[0]=syms[0]<<4|syms[3];
				yuv[1]=errors[1]=syms[1]<<4|syms[4];
				yuv[2]=errors[2]=syms[2]<<4|syms[5];
				yuv[0]+=preds[0];
				uoffset=cu0*yuv[0]>>2;
				preds[1]+=uoffset; CLAMP2(preds[1], 0, 255);
				yuv[1]+=preds[1];
				voffset=(cv0*yuv[0]+cv1*yuv[1])>>2;
				preds[2]+=voffset; CLAMP2(preds[2], 0, 255);
				yuv[2]+=preds[2];
				curr[0]=yuv[0];
				curr[1]=yuv[1]-uoffset;
				curr[2]=yuv[2]-voffset;
				acme_write(&wtptr, 3, fdst, (uint64_t)yuv[0]<<yshift|(uint64_t)yuv[1]<<ushift|(uint64_t)yuv[2]<<vshift);
				guide_check(buffer+idx, kx, ky);
			}
#ifdef ENABLE_GRCTX
			{
				__m128i me=_mm_load_si128((__m128i*)errors);
				me=_mm_slli_epi32(me, 24);
				me=_mm_srai_epi32(me, 24);
				me=_mm_xor_si128(_mm_add_epi32(me, me), _mm_srai_epi32(me, 31));
				me=_mm_slli_epi32(me, 3);
				__m128i meNEE	=_mm_load_si128((__m128i*)eNEE);
				__m128i meNEEE	=_mm_load_si128((__m128i*)eNEEE);
				__m128i meW	=_mm_load_si128((__m128i*)eW);
				meNEE=_mm_max_epi32(meNEE, meNEEE);
				meW=_mm_slli_epi32(meW, 1);
				me=_mm_add_epi32(me, meW);
				me=_mm_add_epi32(me, meNEE);
				me=_mm_srli_epi32(me, 2);
				_mm_store_si128((__m128i*)ecurr, me);
			}
#endif
#ifdef USE_MIX4
			{
				__m128i me=_mm_load_si128((__m128i*)curr);
				me=_mm_sub_epi32(me, mp);
				co0=_mm_add_epi32(co0, _mm_sign_epi32(e0, me));
				co1=_mm_add_epi32(co1, _mm_sign_epi32(e1, me));
				co2=_mm_add_epi32(co2, _mm_sign_epi32(e2, me));
				co3=_mm_add_epi32(co3, _mm_sign_epi32(e3, me));
				_mm_store_si128((__m128i*)coeffs[0], co0);
				_mm_store_si128((__m128i*)coeffs[1], co1);
				_mm_store_si128((__m128i*)coeffs[2], co2);
				_mm_store_si128((__m128i*)coeffs[3], co3);
			}
#endif
			__m256i mixin0=_mm256_load_si256((__m256i*)mixinCDFs[syms[0]]);
			__m256i mixin1=_mm256_load_si256((__m256i*)mixinCDFs[syms[1]]);
			__m256i mixin2=_mm256_load_si256((__m256i*)mixinCDFs[syms[2]]);
			__m256i mixin3=_mm256_load_si256((__m256i*)mixinCDFs[syms[3]]);
			__m256i mixin4=_mm256_load_si256((__m256i*)mixinCDFs[syms[4]]);
			__m256i mixin5=_mm256_load_si256((__m256i*)mixinCDFs[syms[5]]);
			__m256i ms0=_mm256_load_si256((__m256i*)currstats[0]);
			__m256i ms1=_mm256_load_si256((__m256i*)currstats[1]);
			__m256i ms2=_mm256_load_si256((__m256i*)currstats[2]);
			__m256i ms3=_mm256_load_si256((__m256i*)currstats[3]);
			__m256i ms4=_mm256_load_si256((__m256i*)currstats[4]);
			__m256i ms5=_mm256_load_si256((__m256i*)currstats[5]);
			mixin0=_mm256_sub_epi16(mixin0, ms0);
			mixin1=_mm256_sub_epi16(mixin1, ms1);
			mixin2=_mm256_sub_epi16(mixin2, ms2);
			mixin3=_mm256_sub_epi16(mixin3, ms3);
			mixin4=_mm256_sub_epi16(mixin4, ms4);
			mixin5=_mm256_sub_epi16(mixin5, ms5);
			__m256i tmp=_mm256_set1_epi16(1<<MIXNSHIFT>>1);
			mixin0=_mm256_add_epi16(mixin0, tmp);
			mixin1=_mm256_add_epi16(mixin1, tmp);
			mixin2=_mm256_add_epi16(mixin2, tmp);
			mixin3=_mm256_add_epi16(mixin3, tmp);
			mixin4=_mm256_add_epi16(mixin4, tmp);
			mixin5=_mm256_add_epi16(mixin5, tmp);
			mixin0=_mm256_srai_epi16(mixin0, MIXNSHIFT);
			mixin1=_mm256_srai_epi16(mixin1, MIXNSHIFT);
			mixin2=_mm256_srai_epi16(mixin2, MIXNSHIFT);
			mixin3=_mm256_srai_epi16(mixin3, MIXNSHIFT);
			mixin4=_mm256_srai_epi16(mixin4, MIXNSHIFT);
			mixin5=_mm256_srai_epi16(mixin5, MIXNSHIFT);
			ms0=_mm256_add_epi16(ms0, mixin0);
			ms1=_mm256_add_epi16(ms1, mixin1);
			ms2=_mm256_add_epi16(ms2, mixin2);
			ms3=_mm256_add_epi16(ms3, mixin3);
			ms4=_mm256_add_epi16(ms4, mixin4);
			ms5=_mm256_add_epi16(ms5, mixin5);
			_mm256_store_si256((__m256i*)currstats[0], ms0);
			_mm256_store_si256((__m256i*)currstats[1], ms1);
			_mm256_store_si256((__m256i*)currstats[2], ms2);
			_mm256_store_si256((__m256i*)currstats[3], ms3);
			_mm256_store_si256((__m256i*)currstats[4], ms4);
			_mm256_store_si256((__m256i*)currstats[5], ms5);
			rows[0]+=NROWS*NVAL;
			rows[1]+=NROWS*NVAL;
			rows[2]+=NROWS*NVAL;
			rows[3]+=NROWS*NVAL;
		}
	}
	if(fwd)
		acme_write(&wtptr, sizeof(low), fdst, low<<32|low>>32);
	if(wtptr>wtbuf)
		fwrite(wtbuf, 1, wtptr-wtbuf, fdst);
	fclose(fsrc);
	fclose(fdst);
	_mm_free(pixels);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		int64_t csize=0;
		struct stat info={0};

		stat(dstfn, &info);
		csize=info.st_size;
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
	(void)usize;
	(void)rct_names;
	(void)och_names;
	(void)&time_sec;
	return 0;
}