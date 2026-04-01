#ifdef _MSC_VER
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif
#elif defined __linux__ && !defined _GNU_SOURCE
#	define _GNU_SOURCE
#endif
#include<stddef.h>
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
#ifdef PROFILER
#include"util.h"
#endif


#ifdef _MSC_VER
	#define LOUD
//	#define PROFILE_SIZE
//	#define RICE_OPT
	#define ANALYSIS2

	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif


enum
{
	BLOCKX=300,
	BLOCKY=300,

	GRBITS=1,

	XPAD=8,
	NCH=3,
	NROWS=1,
	NVAL=2,

	BUFSIZE=512*1024,

	R3DEPTH=8,
	R3LIMIT=10,
	R3CODEMAX=R3LIMIT+R3DEPTH,//CODEMAX*3 <= 64-8
};

//runtime
#if 1

//clobbers A B C
#define MEDIAN3V_CLOB(M, A, B, C)\
	M=A, A=A<B?A:B, B=B<M?B:M,\
	M=B, B=C<B?C:B, C=T>C?T:C,\
	M=A, A=B<A?B:A, M=T>B?T:B

#define CLAMP2(X, LO, HI) X=X>LO?X:LO, X=X<HI?X:HI
#ifdef _MSC_VER
#	define ALIGN(N) __declspec(align(N))
#	define INLINE __forceinline static
#	define LIKELY(C) C
#else
#	define ALIGN(N) __attribute__((aligned(N)))
#	define INLINE __attribute__((always_inline)) inline static
#	define LIKELY(C) __builtin_expect(C, 1)
#	ifndef _countof
#		define _countof(A) (sizeof(A)/sizeof(*(A)))
#	endif
#endif
#ifndef FLOOR_LOG2
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
#endif
#if 1
INLINE int median3i(int a, int b, int c)
{
	int t;

	t=a;
	a=a>b?b:a;
	b=t>b?t:b;

	t=b;
	b=b>c?c:b;
	c=t>c?t:c;

	t=a;
	a=a>b?b:a;
	b=t>b?t:b;

	return b;
}
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
static double time_sec2(void)
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


/*
overflow:
|                    ______left______   ______right_____
|                   /                \ /                \
|buf1start ... ... [datastart  buf1end|buf2start  dataend] ...
|                   \________________    _______________/
|                                    size
*/
static uint8_t rdbuf[BUFSIZE+sizeof(uint64_t)]={0}, wtbuf[BUFSIZE+sizeof(uint64_t)]={0};
INLINE uint64_t acme_read(uint8_t **pptr, ptrdiff_t size, FILE *f)
{
	uint8_t *ptr=*pptr;
	uint64_t data=*(uint64_t*)ptr;
	ptrdiff_t left;

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
	int64_t counters[OCH_COUNT]={0};
	int prev[OCH_COUNT]={0};
	int rowstride=3*iw;
	uint8_t *ptr=image+rowstride, *end=image+(ptrdiff_t)rowstride*ih;

	while(ptr<end)
	{
		int
			r=(ptr[0]-ptr[0-rowstride])<<2,
			g=(ptr[1]-ptr[1-rowstride])<<2,
			b=(ptr[2]-ptr[2-rowstride])<<2,
			rg=r-g,
			gb=g-b,
			br=b-r;
		ptr+=3;
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

#ifdef ANALYSIS2
static void analysis2(uint8_t *image, int iw, int ih)
{
	enum
	{
		A2_NCH=OCH_COUNT,
		A2_NROWS=2,
		A2_NVAL=2,
	};
	int xblocks=iw/BLOCKX, yblocks=ih/BLOCKY;
	if(xblocks<1)xblocks=1;
	if(yblocks<1)yblocks=1;
	int nblocks=xblocks*yblocks;
	int cssize=nblocks*sizeof(uint32_t[OCH_COUNT]);
	uint32_t *csizes=(uint32_t*)malloc(cssize);
	int psize=(2*BLOCKX+2*XPAD)*(int)sizeof(int16_t[A2_NCH*A2_NROWS*A2_NVAL]);
	int16_t *pixels=(int16_t*)malloc(psize);
	if(!csizes||!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(csizes, 0, cssize);
	memset(pixels, 0, psize);
	for(int by=0, y1=0;by<yblocks;++by)
	{
		int y2=y1+BLOCKY;
		if(by==yblocks-1)
			y2=ih;
		int dy=y2-y1;
		for(int bx=0, x1=0;bx<xblocks;++bx)
		{
			int x2=x1+BLOCKX;
			if(bx==xblocks-1)
				x2=iw;
			int dx=x2-x1;
			int bidx=xblocks*by+bx;
			for(int ky=0;ky<dy;++ky)
			{
				int16_t *rptr0=pixels+(XPAD*A2_NCH*A2_NROWS+(ky+0LL+A2_NROWS)%A2_NROWS)*A2_NVAL;
				int16_t *rptr1=pixels+(XPAD*A2_NCH*A2_NROWS+(ky-1LL+A2_NROWS)%A2_NROWS)*A2_NVAL;
				uint8_t *imptr=image+3*(iw*(y1+ky)+x1);
				for(int kx=0;kx<dx;++kx, imptr+=3)
				{
					int r=imptr[0];
					int g=imptr[1];
					int b=imptr[2];
					int rg=r-g;
					int gb=g-b;
					int br=b-r;
					int16_t curr[OCH_COUNT]=
					{
						r,
						g,
						b,
						rg,
						gb,
						br,
						rg+(gb>>2),//r-(3*g+b)/4 = r-g-(b-g)/4
						rg+(br>>2),//g-(3*r+b)/4 = g-r-(b-r)/4
						br+(rg>>2),//b-(3*r+g)/4 = b-r-(g-r)/4
						br+(gb>>2),//r-(g+3*b)/4 = r-b-(g-b)/4
						gb+(br>>2),//g-(r+3*b)/4 = g-b-(r-b)/4
						gb+(rg>>2),//b-(r+3*g)/4 = b-g-(r-g)/4
						(rg-br)>>1,//r-(g+b)/2 = (r-g + r-b)/2
						(gb-rg)>>1,//g-(r+b)/2 = (g-r + g-b)/2
						(br-gb)>>1,//b-(r+g)/2 = (b-r + b-g)/2
					};
					for(int kc=0;kc<OCH_COUNT;++kc)
					{
						int16_t
							NW	=rptr1[0-1*A2_NCH*A2_NROWS*A2_NVAL],
							N	=rptr1[0+0*A2_NCH*A2_NROWS*A2_NVAL],
							W	=rptr0[0-1*A2_NCH*A2_NROWS*A2_NVAL],
							eNEE	=rptr1[1+2*A2_NCH*A2_NROWS*A2_NVAL],
							eW	=rptr0[1-1*A2_NCH*A2_NROWS*A2_NVAL];
						int nbypass=FLOOR_LOG2((eW>>GRBITS)+1);
						int pred=median3i(N, W, N+W-NW);
						int sym=(int8_t)(curr[kc]-pred);
						sym=sym<<1^sym>>31;
						int nzeros=sym>>nbypass;
						int stopbit=nzeros<R3LIMIT;
						int nbypass2=stopbit?nbypass:R3DEPTH;
						csizes[OCH_COUNT*bidx+kc]+=nzeros+stopbit+nbypass2;

						rptr0[0+0*A2_NCH*A2_NROWS*A2_NVAL]=curr[kc];
						rptr0[1+0*A2_NCH*A2_NROWS*A2_NVAL]=(2*eW+(sym<<GRBITS)+eNEE)>>2;
						rptr0+=A2_NROWS*A2_NVAL;
						rptr1+=A2_NROWS*A2_NVAL;
					}
				}
			}
			x1=x2;
		}
		y1=y2;
	}
	int rsize=nblocks*(int)sizeof(uint32_t[2]);
	uint32_t *results=(int*)malloc(rsize);
	if(!results)
	{
		CRASH("Alloc error");
		return;
	}
	memset(results, 0, rsize);
	uint32_t total=0;
	for(int kb=0;kb<nblocks;++kb)
	{
		int x1=kb%xblocks*BLOCKX, y1=kb/xblocks*BLOCKY;
		uint32_t *csizes2=csizes+OCH_COUNT*kb;
		int bestrct=0;
		uint32_t bestsize=0;
		for(int kt=0;kt<RCT_COUNT;++kt)
		{
			const uint8_t *rct=rct_combinations[kt];
			uint32_t currsize=
				+csizes2[rct[0]]
				+csizes2[rct[1]]
				+csizes2[rct[2]]
			;
			if(!kt||bestsize>currsize)
			{
				bestsize=currsize;
				bestrct=kt;
			}
		}
		printf("%3d  (%5d, %5d)  %12.2lf  %2d %s\n", kb, x1, y1, bestsize/8., bestrct, rct_names[bestrct]);
		total+=bestsize+1;
		results[2*kb+0]=bestsize;
		results[2*kb+1]=bestrct;
	}
	printf("%12.2lf  Total\n", total/8.);
	free(csizes);
	free(pixels);
	free(results);//
}
#endif


#ifdef RICE_OPT
enum
{
	OPT_QSIZE=256,
};
static int8_t ropt_queue[3][OPT_QSIZE];
#endif
static uint8_t clamptable[1024];
int c57_codec(int argc, char **argv)
{
	const uint16_t tag='5'|'7'<<8;

	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	uint64_t c=0;
	int fwd=0, iw=0, ih=0, bestrct=0;
	int64_t usize=0, csize=0;
	int psize=0;
	int16_t *pixels=0;
	uint8_t *image=0, *imptr=0;
	int yidx=0, uidx=0, vidx=0, uc0=0, vc0=0, vc1=0;
	int ysh=0, ush=0, vsh=0;
	uint8_t *rdptr=rdbuf+BUFSIZE, *wtptr=wtbuf;
	uint64_t cache=0;
	int nbits=0;
	uint8_t *const clampptr=clamptable+256;
#ifdef PROFILER
	void *prof_ctx=prof_start();
#endif
#ifdef LOUD
	double t=0;
#endif

	(void)csize;
	(void)&time_sec2;
	(void)och_names;
	(void)rct_names;
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
	t=time_sec2();
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
	usize=(int64_t)3*iw*ih;
	psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	if(fwd)
	{
		image=(uint8_t*)malloc(usize);
		if(!image)
		{
			CRASH("Alloc error");
			return 1;
		}
		fread(image, 1, usize, fsrc);
		guide_save(image, iw, ih);
		bestrct=crct_analysis(image, iw, ih);
#ifdef ANALYSIS2
		analysis2(image, iw, ih);
#endif
	}
	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
	if(fwd)
	{
		csize=0;
		csize+=fwrite(&tag, 1, 2, fdst);
		csize+=fwrite(&iw, 1, 3, fdst);
		csize+=fwrite(&ih, 1, 3, fdst);
		csize+=fwrite(&bestrct, 1, 1, fdst);
	}
	else
	{
		usize+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	}
	yidx=rct_combinations[bestrct][II_PERM_Y];
	uidx=rct_combinations[bestrct][II_PERM_U];
	vidx=rct_combinations[bestrct][II_PERM_V];
	uc0=rct_combinations[bestrct][II_COEFF_U_SUB_Y];
	vc0=rct_combinations[bestrct][II_COEFF_V_SUB_Y];
	vc1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
	for(int k=0;k<_countof(clamptable);++k)
	{
		int val=k-256;
		CLAMP2(val, 0, 255);
		clamptable[k]=val;
	}
	memset(pixels, 0, psize);
	if(fwd)
	{
#ifdef PROFILE_SIZE
		int symcount=0;
		int unarysize=0, binarysize=0;
#endif
#ifdef RICE_OPT
		int qlen=0;
		int roptusize=0, roptbsize=0;
#endif
		imptr=image;
		for(int ky=0;ky<ih;++ky)
		{
			int16_t *rptr=pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL;
			uint8_t *end=imptr+3*iw-2;
			int
				NW0=0, W0=0, eW0=4<<GRBITS,
				NW1=0, W1=0, eW1=4<<GRBITS,
				NW2=0, W2=0, eW2=4<<GRBITS;
#ifdef _MSC_VER
			int kx=0;
#endif
			while(imptr<end)
			{
				int uoffset, voffset,
					N0, pred0, vmin0, vmax0, x0,
					N1, pred1, vmin1, vmax1, x1,
					N2, pred2, vmin2, vmax2, x2;
				
				int nbypass0=FLOOR_LOG2((eW0>>GRBITS)+1);
				int nbypass1=FLOOR_LOG2((eW1>>GRBITS)+1);
				int nbypass2=FLOOR_LOG2((eW2>>GRBITS)+1);
				N0=rptr[0+(0+0*NCH)*NROWS*NVAL];
				N1=rptr[0+(1+0*NCH)*NROWS*NVAL];
				N2=rptr[0+(2+0*NCH)*NROWS*NVAL];
				pred0=N0+W0-NW0;
				pred1=N1+W1-NW1;
				pred2=N2+W2-NW2;
				vmax0=N0, vmin0=W0;
				vmax1=N1, vmin1=W1;
				vmax2=N2, vmin2=W2;
				if(N0<W0)vmin0=N0, vmax0=W0;
				if(N1<W1)vmin1=N1, vmax1=W1;
				if(N2<W2)vmin2=N2, vmax2=W2;
				CLAMP2(pred0, vmin0, vmax0);
				CLAMP2(pred1, vmin1, vmax1);
				CLAMP2(pred2, vmin2, vmax2);

				x0=imptr[yidx];
				x1=imptr[uidx];
				x2=imptr[vidx];
				uoffset=uc0*x0>>2;
				voffset=(vc0*x0+vc1*x1)>>2;
				rptr[0+(0+0*NCH)*NROWS*NVAL]=W0=x0;
				rptr[0+(1+0*NCH)*NROWS*NVAL]=W1=x1-uoffset;
				rptr[0+(2+0*NCH)*NROWS*NVAL]=W2=x2-voffset;
				pred1+=uoffset;
				pred2+=voffset;
				CLAMP2(pred1, 0, 255);
				CLAMP2(pred2, 0, 255);
				//pred1=clampptr[pred1];
				//pred2=clampptr[pred2];
				x0=(int8_t)(x0-pred0);
				x1=(int8_t)(x1-pred1);
				x2=(int8_t)(x2-pred2);
				x0=x0<<1^x0>>31;
				x1=x1<<1^x1>>31;
				x2=x2<<1^x2>>31;

#ifdef FIFOVAL
				fifoval_enqueue(x2<<16^x1<<8^x0);
#endif
				{
					uint64_t code;
					int len;

					int s0=x0, s1=x1, s2=x2;
					int nzeros0=x0>>nbypass0;
					int nzeros1=x1>>nbypass1;
					int nzeros2=x2>>nbypass2;
					int stopbit0=nzeros0<R3LIMIT;
					nzeros0=stopbit0?nzeros0:R3LIMIT;
					nbypass0=stopbit0?nbypass0:R3DEPTH;

					int stopbit1=nzeros1<R3LIMIT;
					nzeros1=stopbit1?nzeros1:R3LIMIT;
					nbypass1=stopbit1?nbypass1:R3DEPTH;

					int stopbit2=nzeros2<R3LIMIT;
					nzeros2=stopbit2?nzeros2:R3LIMIT;
					nbypass2=stopbit2?nbypass2:R3DEPTH;
					int mask0=255<<nbypass0>>8;
					int mask1=255<<nbypass1>>8;
					int mask2=255<<nbypass2>>8;
					s0&=mask0;
					s1&=mask1;
					s2&=mask2;

					code=s2;			len=nbypass2;
					code=code<<stopbit2|stopbit2;	len+=stopbit2;
					code<<=nzeros2;			len+=nzeros2;
					code=code<<nbypass1|s1;		len+=nbypass1;
					code=code<<stopbit1|stopbit1;	len+=stopbit1;
					code<<=nzeros1;			len+=nzeros1;
					code=code<<nbypass0|s0;		len+=nbypass0;
					code=code<<stopbit0|stopbit0;	len+=stopbit0;
					code<<=nzeros0;			len+=nzeros0;

#ifdef PROFILE_SIZE
					symcount+=3;
					unarysize+=nzeros0+stopbit0+nzeros1+stopbit1+nzeros2+stopbit2;
					binarysize+=nbypass0+nbypass1+nbypass2;
#endif
#ifdef RICE_OPT
					ropt_queue[0][qlen]=x0;
					ropt_queue[1][qlen]=x1;
					ropt_queue[2][qlen]=x2;
					++qlen;
					if(qlen>=OPT_QSIZE)
					{
						qlen=0;
						for(int kc=0;kc<3;++kc)
						{
							int ricek=0;
							for(int k=0;k<OPT_QSIZE;++k)
								ricek+=ropt_queue[kc][k];
							ricek/=OPT_QSIZE;
							ricek=FLOOR_LOG2(ricek+1);
							for(int k=0;k<OPT_QSIZE;++k)
							{
								int val=ropt_queue[kc][k];
								int nzeros=val>>ricek;
								roptusize+=nzeros+1;
								roptbsize+=ricek;
							}
						}
						roptbsize+=3;
					}
#endif

					cache|=code<<nbits;
					nbits+=len;
					int spill=nbits>=64;
					if(spill)
						acme_write(&wtptr, 8, fdst, cache);
					int rem=nbits-64;
					uint64_t next=code>>(len-rem);
					cache=spill?next:cache;
					nbits=spill?rem:nbits;
				}

				rptr[1+(0+0*NCH)*NROWS*NVAL]=eW0=(2*eW0+(x0<<GRBITS)+rptr[1+(0+2*NCH)*NROWS*NVAL])>>2;
				rptr[1+(1+0*NCH)*NROWS*NVAL]=eW1=(2*eW1+(x1<<GRBITS)+rptr[1+(1+2*NCH)*NROWS*NVAL])>>2;
				rptr[1+(2+0*NCH)*NROWS*NVAL]=eW2=(2*eW2+(x2<<GRBITS)+rptr[1+(2+2*NCH)*NROWS*NVAL])>>2;

				imptr+=3;
				rptr+=NCH*NROWS*NVAL;
				NW0=N0;
				NW1=N1;
				NW2=N2;
#ifdef _MSC_VER
				++kx;
#endif
			}
		}
		acme_write(&wtptr, 8, fdst, cache);
#ifdef PROFILE_SIZE
		printf("%12.2lf bytes  %12.4lf BPD  unary\n", (double)unarysize/8, (double)unarysize/symcount);
		printf("%12.2lf bytes  %12.4lf BPD  binary\n", (double)binarysize/8, (double)binarysize/symcount);
		printf("%12.2lf bytes  %12.4lf BPD  Rice\n", (double)(unarysize+binarysize)/8, (double)(unarysize+binarysize)/symcount);

		printf("%12.2lf bytes  %12.4lf BPD  ropt_unary\n", (double)roptusize/8, (double)roptusize/symcount);
		printf("%12.2lf bytes  %12.4lf BPD  ropt_binary\n", (double)roptbsize/8, (double)roptbsize/symcount);
		printf("%12.2lf bytes  %12.4lf BPD  ropt_total\n", (double)(roptusize+roptbsize)/8, (double)(roptusize+roptbsize)/symcount);
#endif
	}
	else//dec
	{
		ysh=yidx*8;
		ush=uidx*8;
		vsh=vidx*8;
		acme_read(&rdptr, 8, fsrc);
		cache=((volatile uint64_t*)rdptr)[-1];
		nbits=0;//64&7
		for(int ky=0;ky<ih;++ky)
		{
			int16_t *rptr=pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL;
			int
				NW0=0, W0=0, eW0=4<<GRBITS,
				NW1=0, W1=0, eW1=4<<GRBITS,
				NW2=0, W2=0, eW2=4<<GRBITS;
			int kx=0;
			while(kx<iw)
			{
				int uoffset, voffset,
					N0, pred0, vmin0, vmax0, x0,
					N1, pred1, vmin1, vmax1, x1,
					N2, pred2, vmin2, vmax2, x2;
				
				int nbypass0=FLOOR_LOG2((eW0>>GRBITS)+1);
				int nbypass1=FLOOR_LOG2((eW1>>GRBITS)+1);
				int nbypass2=FLOOR_LOG2((eW2>>GRBITS)+1);
				N0=rptr[0+(0+0*NCH)*NROWS*NVAL];
				N1=rptr[0+(1+0*NCH)*NROWS*NVAL];
				N2=rptr[0+(2+0*NCH)*NROWS*NVAL];
				pred0=N0+W0-NW0;
				pred1=N1+W1-NW1;
				pred2=N2+W2-NW2;
				vmax0=N0, vmin0=W0;
				vmax1=N1, vmin1=W1;
				vmax2=N2, vmin2=W2;
				if(N0<W0)vmin0=N0, vmax0=W0;
				if(N1<W1)vmin1=N1, vmax1=W1;
				if(N2<W2)vmin2=N2, vmax2=W2;
				CLAMP2(pred0, vmin0, vmax0);
				CLAMP2(pred1, vmin1, vmax1);
				CLAMP2(pred2, vmin2, vmax2);
#if 1
				int s0, s1, s2;
				{
					int nzeros, stopbit;
					uint64_t code=cache>>nbits;

					nzeros=(int)_tzcnt_u64(code);
					stopbit=nzeros<R3LIMIT;
					nbypass0=stopbit?nbypass0:R3DEPTH;
					nzeros=stopbit?nzeros:R3LIMIT;
					code>>=nzeros+stopbit;
					s0=(uint8_t)(nzeros<<nbypass0|((uint8_t)code&255<<nbypass0>>8));
					code>>=nbypass0;
					nbits+=(uint64_t)nzeros+stopbit+nbypass0;

					nzeros=(int)_tzcnt_u64(code);
					stopbit=nzeros<R3LIMIT;
					nbypass1=stopbit?nbypass1:R3DEPTH;
					nzeros=stopbit?nzeros:R3LIMIT;
					code>>=nzeros+stopbit;
					s1=(uint8_t)(nzeros<<nbypass1|((uint8_t)code&255<<nbypass1>>8));
					code>>=nbypass1;
					nbits+=(uint64_t)nzeros+stopbit+nbypass1;

					nzeros=(int)_tzcnt_u64(code);
					stopbit=nzeros<R3LIMIT;
					nbypass2=stopbit?nbypass2:R3DEPTH;
					nzeros=stopbit?nzeros:R3LIMIT;
					code>>=nzeros+stopbit;
					s2=(uint8_t)(nzeros<<nbypass2|((uint8_t)code&255<<nbypass2>>8));
					code>>=nbypass2;
					nbits+=(uint64_t)nzeros+stopbit+nbypass2;
#ifdef FIFOVAL
					fifoval_check(s2<<16^s1<<8^s0);
#endif
					int nbytes=nbits>>3;
					if(nbytes)
					{
						int sh=nbytes<<3;
						cache=(acme_read(&rdptr, nbytes, fsrc)&((1LL<<sh)-1))<<(64-sh)|cache>>sh;
						nbits&=7;
					}
				}
#else
				int s0=r2_dec(&ec, nbypass0);
				int s1=r2_dec(&ec, nbypass1);
				int s2=r2_dec(&ec, nbypass2);
#endif
				x0=s0>>1^-(s0&1);
				x1=s1>>1^-(s1&1);
				x2=s2>>1^-(s2&1);
				x0=(uint8_t)(x0+pred0);
				uoffset=uc0*x0>>2;
				pred1+=uoffset;
				//CLAMP2(pred1, 0, 255);
				pred1=clampptr[pred1];
				x1=(uint8_t)(x1+pred1);
				voffset=(vc0*x0+vc1*x1)>>2;
				pred2+=voffset;
				//CLAMP2(pred2, 0, 255);
				pred2=clampptr[pred2];
				x2=(uint8_t)(x2+pred2);
				acme_write(&wtptr, 3, fdst, (uint64_t)x2<<vsh|(uint64_t)x1<<ush|(uint64_t)x0<<ysh);
				rptr[0+(0+0*NCH)*NROWS*NVAL]=W0=x0;
				rptr[0+(1+0*NCH)*NROWS*NVAL]=W1=x1-uoffset;
				rptr[0+(2+0*NCH)*NROWS*NVAL]=W2=x2-voffset;
				
				rptr[1+(0+0*NCH)*NROWS*NVAL]=eW0=(2*eW0+(s0<<GRBITS)+rptr[1+(0+2*NCH)*NROWS*NVAL])>>2;
				rptr[1+(1+0*NCH)*NROWS*NVAL]=eW1=(2*eW1+(s1<<GRBITS)+rptr[1+(1+2*NCH)*NROWS*NVAL])>>2;
				rptr[1+(2+0*NCH)*NROWS*NVAL]=eW2=(2*eW2+(s2<<GRBITS)+rptr[1+(2+2*NCH)*NROWS*NVAL])>>2;
				rptr+=NCH*NROWS*NVAL;
				NW0=N0;
				NW1=N1;
				NW2=N2;
				++kx;
			}
		}
	}
	if(wtptr>wtbuf)
		fwrite(wtbuf, 1, wtptr-wtbuf, fdst);
	fclose(fsrc);
	fclose(fdst);
	free(pixels);
	if(fwd)
		free(image);
#ifdef LOUD
	t=time_sec2()-t;
	if(fwd)
	{
		struct stat info={0};
		stat(dstfn, &info);
		csize=info.st_size;
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
#ifdef PROFILER
	prof_end(prof_ctx);
#endif
	return 0;
}
