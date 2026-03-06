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
#ifdef PROFILER
void* prof_start();
void prof_end(void *prof_ctx);
#endif


#ifdef _MSC_VER
	#define LOUD
	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif


//	#define USE_GAMMA	//X
//	#define USE_LUT
	#define USE_L1
	#define ENABLE_CRCT


#ifdef USE_L1
#define PREDLIST\
	PRED(W)\
	PRED(N+W-NW)\
	PRED(2*N-NN)\
	PRED(NE)\

#endif
enum
{
#ifdef USE_L1
	SHIFT=24,
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED
#endif
	GRBITS=5,

	XPAD=8,
	NCH=3,
	NROWS=4,
	NVAL=2,

	BUFSIZE=512*1024,
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
static uint16_t *g_image=0;
static double g_sqe[3]={0};
static void guide_save(FILE *f, int iw, int ih)
{
	ptrdiff_t idx=0, size=0;
	
	size=(ptrdiff_t)6*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(uint16_t*)malloc(size);
	if(!g_image)
	{
		CRASH("Alloc error");
		return;
	}
	idx=ftell(f);
	fread(g_image, 1, size, f);
	fseek(f, (long)idx, SEEK_SET);
	for(ptrdiff_t k=0, res=(ptrdiff_t)3*iw*ih-2;k<res;k+=3)
	{
		g_image[k+0]=(uint16_t)(g_image[k+0]<<8|g_image[k+0]>>8);
		g_image[k+1]=(uint16_t)(g_image[k+1]<<8|g_image[k+1]>>8);
		g_image[k+2]=(uint16_t)(g_image[k+2]<<8|g_image[k+2]>>8);
	}
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
static int crct_analysis(FILE *f, int iw, int ih)
{
	int64_t counters[OCH_COUNT]={0};
	int prev[OCH_COUNT]={0};
	uint8_t *ptr=rdbuf+BUFSIZE;
	long idx=ftell(f);

	for(ptrdiff_t k=0, size=(ptrdiff_t)6*iw*ih;k<size;k+=6)
	{
		int r, g, b, rg, gb, br;

		uint64_t data=acme_read(&ptr, 6, f);
		r=data>> 0&65535;
		g=data>>16&65535;
		b=data>>32&65535;
		r=(uint16_t)(r<<8|r>>8);
		g=(uint16_t)(g<<8|g>>8);
		b=(uint16_t)(b<<8|b>>8);
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


#ifdef USE_LUT
typedef struct _CSymInfo
{
	uint16_t sym;
	uint8_t nzeros, bypass;
} CSymInfo;
static CSymInfo csymtable[8][256];
static int8_t dsymtable[256];
#endif
int c55_codec(int argc, char **argv)
{
	const uint16_t tag='5'|'5'<<8;

	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	uint64_t c=0;
	int fwd=0, iw=0, ih=0;
	int64_t usize=0, csize=0;
	int psize=0;
	int32_t *pixels=0;
	uint64_t cache=0;
	int nbits=0;
	uint8_t *rdptr=0, *wtptr=0;
#ifdef ENABLE_CRCT
	int bestrct=0, yidx=0, uidx=0, vidx=0, uc0=0, vc0=0, vc1=0;
#endif
#ifdef USE_L1
	int32_t coeffs[3][NPREDS]={0}, estim[NPREDS]={0};
	int64_t p1=0;
	int j=0;
#endif
#ifdef LOUD
	double t=0;
#endif
#ifdef PROFILER
	void *prof_ctx=prof_start();
#endif

	if(argc!=3)
	{
		printf(
			"Usage:  \"%s\"  src  dst\n"
			"Only for 48-bit PPM images\n"
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
		c|=(int64_t)fgetc(fsrc)<<8*5;
		c|=(int64_t)fgetc(fsrc)<<8*6;
		if(c!=(
			(uint64_t)'\n'<<8*0|
			(uint64_t) '6'<<8*1|
			(uint64_t) '5'<<8*2|
			(uint64_t) '5'<<8*3|
			(uint64_t) '3'<<8*4|
			(uint64_t) '5'<<8*5|
			(uint64_t)'\n'<<8*6
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
#ifdef ENABLE_CRCT
		fread(&bestrct, 1, 1, fsrc);
#endif
	}
	if(iw<1||ih<1)
	{
		CRASH("Unsupported source file");
		return 1;
	}
	usize=(int64_t)6*iw*ih;
	psize=(iw+2*XPAD)*(int)sizeof(int32_t[NCH*NROWS*NVAL]);
	pixels=(int32_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	cache=0;
	nbits=64;
	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
	if(fwd)
	{
#ifdef ENABLE_GUIDE
		guide_save(fsrc, iw, ih);
#endif
#ifdef ENABLE_CRCT
		bestrct=crct_analysis(fsrc, iw, ih);
#endif
		fwrite(&tag, 1, 2, fdst);
		fwrite(&iw, 1, 3, fdst);
		fwrite(&ih, 1, 3, fdst);
#ifdef ENABLE_CRCT
		fwrite(&bestrct, 1, 1, fdst);
#endif
#ifdef USE_LUT
		for(int nbypass=0;nbypass<8;++nbypass)
		{
			for(int ks=0;ks<256;++ks)
			{
				CSymInfo *p=csymtable[nbypass]+ks;
				int sym=ks-128;
				sym=sym<<1^sym>>31;
				p->sym=sym<<GRBITS;
				p->nzeros=sym>>nbypass;
				p->bypass=sym&((1<<nbypass)-1);
			}
		}
#endif
	}
	else
	{
		fread(&cache, 1, sizeof(uint64_t), fsrc);
		nbits=0;

		fprintf(fdst, "P6\n%d %d\n65535\n", iw, ih);
#ifdef USE_LUT
		for(int k=0;k<256;++k)
			dsymtable[k]=k>>1^-(k&1);
#endif
	}
#ifdef ENABLE_CRCT
	yidx=rct_combinations[bestrct][II_PERM_Y]*16;
	uidx=rct_combinations[bestrct][II_PERM_U]*16;
	vidx=rct_combinations[bestrct][II_PERM_V]*16;
	uc0=rct_combinations[bestrct][II_COEFF_U_SUB_Y];
	vc0=rct_combinations[bestrct][II_COEFF_V_SUB_Y];
	vc1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
#endif
#ifdef USE_L1
	for(int kc=0;kc<3;++kc)
	{
		for(int kp=0;kp<NPREDS;++kp)
			coeffs[kc][kp]=(1<<SHIFT)/NPREDS;
	}
#endif
	memset(pixels, 0, psize);
	for(int k=0;k<psize/sizeof(int32_t);++k)
		pixels[k]=256;
	rdptr=rdbuf+BUFSIZE;
	wtptr=wtbuf;
	for(int ky=0;ky<ih;++ky)
	{
		int yuv[3]={0};
		int pred=0;
		int32_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=0;kx<iw;++kx)
		{
#ifdef ENABLE_CRCT
			int offset=0;
#endif
			if(fwd)
			{
				uint64_t data=acme_read(&rdptr, 6, fsrc);
#ifdef ENABLE_CRCT
				yuv[0]=data>>yidx&65535;
				yuv[1]=data>>uidx&65535;
				yuv[2]=data>>vidx&65535;
				yuv[0]=(uint16_t)(yuv[0]<<8|yuv[0]>>8);
				yuv[1]=(uint16_t)(yuv[1]<<8|yuv[1]>>8);
				yuv[2]=(uint16_t)(yuv[2]<<8|yuv[2]>>8);
#else
				yuv[0]=data>> 0&65535;
				yuv[1]=data>>16&65535;
				yuv[2]=data>>32&65535;
				yuv[0]-=yuv[1];
				yuv[2]-=yuv[1];
#endif
			}
			for(int kc=0;kc<3;++kc)
			{
				int32_t
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
					eN	=rows[1][1+0*NCH*NROWS*NVAL],
					eNE	=rows[1][1+1*NCH*NROWS*NVAL],
					eNEE	=rows[1][1+2*NCH*NROWS*NVAL],
					eNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
					eWW	=rows[0][1-2*NCH*NROWS*NVAL],
					eW	=rows[0][1-1*NCH*NROWS*NVAL];
				int nbypass;
				int vmin, vmax;
				int error, sym;
				int nzeros, bypass;
#ifdef ENABLE_CRCT
				if(kc==1)
					offset=uc0*yuv[0]>>2;
				if(kc==2)
					offset=(vc0*yuv[0]+vc1*yuv[1])>>2;
#endif
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
				(void)eN;
				(void)eNE;
				(void)eNEE;
				(void)eNEEE;
				(void)eWW;
				(void)eW;
#ifndef USE_GAMMA
				nbypass=FLOOR_LOG2((eW+(1<<GRBITS>>1))>>GRBITS);
				if(nbypass<8)
					nbypass=8;
#endif
#ifdef USE_L1
				p1=1LL<<SHIFT>>1;
#define PRED(E) estim[j]=E; p1+=(int64_t)coeffs[kc][j]*estim[j]; ++j;
				j=0;
				PREDLIST;
#undef  PRED
				p1>>=SHIFT;
				pred=(int)p1;
				vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);
#else
				pred=N+W-NW;
				vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				CLAMP2(pred, vmin, vmax);
				//pred=(N+W)>>1;
#endif
#ifdef ENABLE_CRCT
				pred+=offset;
				CLAMP2(pred, 0, 65535);
#endif
				//if(ky==ih/2&&kx==iw/2)//
				//	printf("");
				if(fwd)
				{
					error=(int16_t)(yuv[kc]-pred);
					sym=error<<1^error>>31;
#ifdef USE_GAMMA
					++sym;
					bypass=sym;
					nzeros=nbypass=31-_lzcnt_u32(sym);
#else
					nzeros=sym>>nbypass;
					bypass=sym&0x7FFFFFFF>>(31-nbypass);
#endif
					if(nzeros>=nbits)//fill the rest of cache with zeros, and flush
					{
						nzeros-=nbits;
						acme_write(&wtptr, sizeof(cache), fdst, cache);
						//fwrite(&cache, 1, sizeof(cache), fdst);
						cache=0;
						while(nzeros>=64)//just flush zeros
						{
							nzeros-=64;
							acme_write(&wtptr, sizeof(cache), fdst, cache);
							//fwrite(&cache, 1, sizeof(cache), fdst);
						}
						nbits=64;
					}
					//now there is room for zeros:  0 <= nzeros < nbits <= 64
					nbits-=nzeros;//emit remaining zeros to cache

					bypass|=1<<nbypass;//append 1 stop bit
					++nbypass;
					if(nbypass>=nbits)//cache would overflow:  fill, flush, and repeat
					{
						nbypass-=nbits;
						cache|=(uint64_t)bypass>>nbypass;
						bypass&=0x7FFFFFFF>>(31-nbypass);
						acme_write(&wtptr, sizeof(cache), fdst, cache);
						//fwrite(&cache, 1, sizeof(cache), fdst);
						cache=0;
						nbits=64;
					}
					//now there is room for bypass:  0 <= nbypass < nbits <= 64
					if(nbypass)
					{
						nbits-=nbypass;//emit remaining bypass to cache
						cache|=(uint64_t)bypass<<nbits;
					}
				}
				else
				{
					sym=-nbits;
					while(!cache)
					{
						sym+=64;
						cache=acme_read(&rdptr, sizeof(cache), fsrc);
						//fread(&cache, 1, sizeof(cache), fsrc);
					}
					nbits=(int)_lzcnt_u64(cache);
					sym+=nbits;
					
#ifdef USE_GAMMA
					nbypass=sym;
					sym=1<<nbypass;
#endif
					sym<<=nbypass;
					cache&=0x7FFFFFFFFFFFFFFF>>nbits;//remove stop bit
					nbits+=nbypass+1;
					if(nbits>=64)//nbits = nbits0+nbypass > N
					{
						//example: 000000[11 1]1010010	nbits=6, nbypass=3	6+3-8 = 1
						nbits-=64;
						sym|=(int)(cache<<nbits);
						cache=acme_read(&rdptr, sizeof(cache), fsrc);
						//fread(&cache, 1, sizeof(cache), fsrc);
						nbypass=nbits;
					}
					if(nbypass)
					{
						sym|=(int)(cache>>(64-nbits));
						cache&=0xFFFFFFFFFFFFFFFF>>nbits;//nbits=61 -> cache&=7;
					}
#ifdef USE_GAMMA
					--sym;
#endif
					error=sym>>1^-(sym&1);
					yuv[kc]=(uint16_t)(error+pred);
#ifdef ENABLE_GUIDE
					int perm[]=
					{
						rct_combinations[bestrct][II_PERM_Y],
						rct_combinations[bestrct][II_PERM_U],
						rct_combinations[bestrct][II_PERM_V],
					};
					if(yuv[kc]!=g_image[3*(iw*ky+kx)+perm[kc]])
					{
						CRASH("Guide  X%d Y%d C%d  %d != %d", kx, ky, kc, yuv[kc], g_image[3*(iw*ky+kx)+perm[kc]]);
						return 1;
					}
#endif
				}
#ifdef ENABLE_CRCT
				rows[0][0]=yuv[kc]-offset;
#else
				rows[0][0]=yuv[kc];
#endif
#ifdef USE_L1
				{
					int e=(rows[0][0]>p1)-(rows[0][0]<p1);

#define PRED(...) coeffs[kc][j]+=e*estim[j]; ++j;
					j=0;
					PREDLIST;
#undef  PRED
				}
#endif
#ifdef USE_LUT
				rows[0][1]=(2*eW+sym+eNEEE)>>2;
#else
				rows[0][1]=(2*eW+(sym<<GRBITS)+eNEEE)>>2;
#endif
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
			}
			if(!fwd)
			{
#ifdef ENABLE_CRCT
				yuv[0]=(uint16_t)(yuv[0]<<8|yuv[0]>>8);
				yuv[1]=(uint16_t)(yuv[1]<<8|yuv[1]>>8);
				yuv[2]=(uint16_t)(yuv[2]<<8|yuv[2]>>8);
				acme_write(&wtptr, 6, fdst, (uint64_t)yuv[2]<<vidx|(uint64_t)yuv[1]<<uidx|(uint64_t)yuv[0]<<yidx);
#else
				yuv[2]+=yuv[1];
				yuv[0]+=yuv[1];
				acme_write(&wtptr, 3, fdst, (uint64_t)yuv[2]<<16|(uint64_t)yuv[1]<<8|yuv[0]);
#endif
			}
		}
	}
	if(fwd)
		acme_write(&wtptr, sizeof(cache), fdst, cache);
		//fwrite(&cache, 1, sizeof(cache), fdst);

	if(wtptr>wtbuf)
		fwrite(wtbuf, 1, wtptr-wtbuf, fdst);
	free(pixels);
	fclose(fsrc);
	fclose(fdst);
#ifdef PROFILER
	prof_end(prof_ctx);
#endif
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		struct stat info={0};
		stat(dstfn, &info);
		csize=info.st_size;
		printf("CWH=3*%d*%d  \"%s\"\n", iw, ih, srcfn);
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
#ifdef ENABLE_CRCT
	(void)rct_names;
	(void)och_names;
#endif
	(void)usize;
	(void)csize;
	(void)&time_sec;
	return 0;
}
