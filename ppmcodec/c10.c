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
#include<intrin.h>
#else
#include<time.h>
#endif
#include<immintrin.h>


#ifdef _MSC_VER
	#define LOUD
	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif


	#define SIGNED_PIXEL


#define PREDLIST\
	PRED(N+W-NW)\
	PRED(N)\
	PRED(W)\
	PRED(W+NE-N)\
	PRED(3*(N-NN)+NNN)\
	PRED(3*(W-WW)+WWW)\
	PRED(N+NE-NNE)\
	PRED(NEE)\
	PRED(NN)\
	PRED(WW)\
	PRED(2*N-NN)\
	PRED(2*W-WW)\
	PRED(NEEE)\
	PRED(NEEEE)\
	PRED(NNWW)\
	PRED(NNEE)\
	PRED(N+NW-NNW)\
	PRED(W+NW-NWW)\
	PRED((WWWW+NEEEE)>>1)\
	PRED((WWW+NNN+NEEE-NW)>>1)\


enum
{
	L1SH=20,
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED

	NCTX=18,

	XPAD=8,
	NROWS=4,
	NCH=3,
	NVAL=2,

	NLEVELS=256,
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
#define CVTFP32_I32(X)  _mm_cvt_ss2si(_mm_set_ss(X))
#define CVTTFP32_I32(X) _mm_cvtt_ss2si(_mm_set_ss(X))
#define CVTFP64_I64(X)  _mm_cvtsd_si64(_mm_set_sd(X))
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
#endif
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

	II_COUNT,
} RCTInfoIdx;
#ifdef ENABLE_EXTENDED_RCT
#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
	OCH(CX40) OCH(C0X4) OCH(C40X)\
	OCH(CX31) OCH(C3X1) OCH(C31X)\
	OCH(CX13) OCH(C1X3) OCH(C13X)\
	OCH(CX22) OCH(C2X2) OCH(C22X)
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
#endif
#ifndef ENABLE_EXTENDED_RCT
typedef enum _OCHIndex
{
	OCH_R,
	OCH_G,
	OCH_B,
	OCH_RG,
	OCH_GB,
	OCH_BR,
	OCH_COUNT,
	OCH_GR=OCH_RG,
	OCH_BG=OCH_GB,
	OCH_RB=OCH_BR,
} OCHIndex;
#endif
#ifdef ENABLE_EXTENDED_RCT
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
#endif
#ifndef ENABLE_EXTENDED_RCT
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
	RCT(_400_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)
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
static const char *rct_names[]=
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
#if 0
	uint8_t *ptr=image;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, ptr+=3)
		{
			int
				r=ptr[0]<<2,
				g=ptr[1]<<2,
				b=ptr[2]<<2,
				rg=r-g,
				gb=g-b,
				br=b-r,
				c6=rg+(gb>>2),//r-(3*g+b)/4 = r-g-(b-g)/4
				c7=rg+(br>>2),//g-(3*r+b)/4 = g-r-(b-r)/4
				c8=br+(rg>>2),//b-(3*r+g)/4 = b-r-(g-r)/4
				c9=br+(gb>>2),//r-(g+3*b)/4 = r-b-(g-b)/4
				cA=gb+(br>>2),//g-(r+3*b)/4 = g-b-(r-b)/4
				cB=gb+(rg>>2),//b-(r+3*g)/4 = b-g-(r-g)/4
				cC=(rg-br)>>1,//r-(g+b)/2 = (r-g + r-b)/2
				cD=(gb-rg)>>1,//g-(r+b)/2 = (g-r + g-b)/2
				cE=(br-gb)>>1;//b-(r+g)/2 = (b-r + b-g)/2
			int pred[OCH_COUNT];
			memcpy(pred, prev, sizeof(pred));
			if(ky)
			{
				int
					Nr=ptr[0]<<2,
					Ng=ptr[1]<<2,
					Nb=ptr[2]<<2,
					Nrg=Nr-Ng,
					Ngb=Ng-Nb,
					Nbr=Nb-Nr,
					Nc6=Nrg+(Ngb>>2),//r-(3*g+b)/4 = r-g-(b-g)/4
					Nc7=Nrg+(Nbr>>2),//g-(3*r+b)/4 = g-r-(b-r)/4
					Nc8=Nbr+(Nrg>>2),//b-(3*r+g)/4 = b-r-(g-r)/4
					Nc9=Nbr+(Ngb>>2),//r-(g+3*b)/4 = r-b-(g-b)/4
					NcA=Ngb+(Nbr>>2),//g-(r+3*b)/4 = g-b-(r-b)/4
					NcB=Ngb+(Nrg>>2),//b-(r+3*g)/4 = b-g-(r-g)/4
					NcC=(Nrg-Nbr)>>1,//r-(g+b)/2 = (r-g + r-b)/2
					NcD=(Ngb-Nrg)>>1,//g-(r+b)/2 = (g-r + g-b)/2
					NcE=(Nbr-Ngb)>>1;//b-(r+g)/2 = (b-r + b-g)/2
				pred[0x0]=(pred[0x0]+Nr)>>1;
				pred[0x1]=(pred[0x1]+Ng)>>1;
				pred[0x2]=(pred[0x2]+Nb)>>1;
				pred[0x3]=(pred[0x3]+Nrg)>>1;
				pred[0x4]=(pred[0x4]+Ngb)>>1;
				pred[0x5]=(pred[0x5]+Nbr)>>1;
				pred[0x6]=(pred[0x6]+Nc6)>>1;
				pred[0x7]=(pred[0x7]+Nc7)>>1;
				pred[0x8]=(pred[0x8]+Nc8)>>1;
				pred[0x9]=(pred[0x9]+Nc9)>>1;
				pred[0xA]=(pred[0xA]+NcA)>>1;
				pred[0xB]=(pred[0xB]+NcB)>>1;
				pred[0xC]=(pred[0xC]+NcC)>>1;
				pred[0xD]=(pred[0xD]+NcD)>>1;
				pred[0xE]=(pred[0xE]+NcE)>>1;
			}
			counters[0x0]+=abs(r	-pred[0x0]);
			counters[0x1]+=abs(g	-pred[0x1]);
			counters[0x2]+=abs(b	-pred[0x2]);
			counters[0x3]+=abs(rg	-pred[0x3]);
			counters[0x4]+=abs(gb	-pred[0x4]);
			counters[0x5]+=abs(br	-pred[0x5]);
			counters[0x6]+=abs(c6	-pred[0x6]);
			counters[0x7]+=abs(c7	-pred[0x7]);
			counters[0x8]+=abs(c8	-pred[0x8]);
			counters[0x9]+=abs(c9	-pred[0x9]);
			counters[0xA]+=abs(cA	-pred[0xA]);
			counters[0xB]+=abs(cB	-pred[0xB]);
			counters[0xC]+=abs(cC	-pred[0xC]);
			counters[0xD]+=abs(cD	-pred[0xD]);
			counters[0xE]+=abs(cE	-pred[0xE]);
			prev[0x0]=r;
			prev[0x1]=g;
			prev[0x2]=b;
			prev[0x3]=rg;
			prev[0x4]=gb;
			prev[0x5]=br;
			prev[0x6]=c6;
			prev[0x7]=c7;
			prev[0x8]=c8;
			prev[0x9]=c9;
			prev[0xA]=cA;
			prev[0xB]=cB;
			prev[0xC]=cC;
			prev[0xD]=cD;
			prev[0xE]=cE;
		}
	}
#endif
#if 1
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
#endif
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

static int32_t hist[3][NCTX][NLEVELS+1];
int c10_codec(int argc, char **argv)
{
	const uint16_t tag='1'|'0'<<8;
	
	const char *srcfn=0, *dstfn=0;
	ptrdiff_t usize=0, overhead=0, csize=0, streamsize=0;
	uint8_t *buf=0, *image=0, *stream=0;
	int iw=0, ih=0, fwd=0;
	uint8_t *streamptr=0;
	uint64_t low=0, range=0, code=0;
	int bestrct=0, yidx=0, uidx=0, vidx=0, uc0=0, vc0=0, vc1=0;
	uint8_t *NNptr=0, *Nptr=0, *ptr=0;
	long long weights[3*NPREDS]={0};
	int psize=0;
	int16_t *pixels=0;
#ifdef LOUD
	double t=time_sec();
#endif

	if(argc!=3)
	{
		printf(
			"Usage: \"%s\"  input  output    Encode/decode.\n"
			, argv[0]
		);
		return 1;
	}
	srcfn=argv[1];
	dstfn=argv[2];
	{
		int64_t c=0, nread=0;
		FILE *fsrc=fopen(srcfn, "rb");
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
			usize=(ptrdiff_t)3*iw*ih;
			overhead=(ptrdiff_t)iw*ih;
			buf=(uint8_t*)malloc(usize+overhead);
			if(!buf)
			{
				CRASH("Alloc error");
				fclose(fsrc);
				return 1;
			}
			memset(buf, 128, overhead);
			stream=buf;
			image=buf+overhead;
			nread=fread(image, 1, usize, fsrc);
			if(nread!=usize)
				printf("Warning: PPM truncated at %td/%td\n", nread, usize);
		}
		else
		{
			struct stat info={0};

			fread(&iw, 1, 3, fsrc);
			fread(&ih, 1, 3, fsrc);
			fread(&bestrct, 1, 1, fsrc);
			stat(srcfn, &info);
			csize=info.st_size;
			streamsize=csize-ftell(fsrc);
			usize=(ptrdiff_t)3*iw*ih;
			overhead=(ptrdiff_t)iw*ih;
			buf=(uint8_t*)malloc(usize+overhead);
			if(!buf)
			{
				CRASH("Alloc error");
				fclose(fsrc);
				return 1;
			}
			memset(buf, 128, overhead);
			image=buf+(ptrdiff_t)9*iw;
			stream=buf+usize+overhead-streamsize;
			nread=fread(stream, 1, streamsize, fsrc);
			if(nread!=streamsize)
				printf("Warning: stream truncated at %td/%td\n", nread, streamsize);
		}
		fclose(fsrc);
	}
	streamptr=stream;
	low=0, range=0xFFFFFFFFFFFF, code=0;
	if(fwd)//analysis
	{
#ifdef ENABLE_GUIDE
		guide_save(image, iw, ih);
#endif
		bestrct=crct_analysis(image, iw, ih);
	}
	else
	{
		code=0;
		code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);
		code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);
	}
	//short *ebuf=(short*)malloc((iw+16LL)*sizeof(short[3]));
	//if(!ebuf)
	//{
	//	CRASH("Alloc error");
	//	return 1;
	//}
	//memset(ebuf, 0, (iw+16LL)*sizeof(short[3]));
	yidx=rct_combinations[bestrct][II_PERM_Y];
	uidx=rct_combinations[bestrct][II_PERM_U];
	vidx=rct_combinations[bestrct][II_PERM_V];
	uc0=rct_combinations[bestrct][II_COEFF_U_SUB_Y];
	vc0=rct_combinations[bestrct][II_COEFF_V_SUB_Y];
	vc1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
	NNptr=image-(ptrdiff_t)6*iw;
	Nptr=image-(ptrdiff_t)3*iw;
	ptr=image;
	for(int k=0;k<3*NPREDS;++k)
		weights[k]=(1<<L1SH)/NPREDS;
	psize=(iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	memset(pixels, 0, psize);
	memset(hist, 0, sizeof(hist));
	for(int ky=0;ky<ih;++ky)
	{
		short *rows[]=
		{
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,//sub 1 channel for pre-increment
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS-NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		int estim[NPREDS]={0};
		int yuv[4]={0};
		int sym, cdf, freq, den;
		for(int kx=0;kx<iw;++kx)
		{
			int offset=0, error;
			if(fwd)
			{
				yuv[0]=ptr[yidx];
				yuv[1]=ptr[uidx];
				yuv[2]=ptr[vidx];
#ifdef SIGNED_PIXEL
				yuv[0]-=128;
				yuv[1]-=128;
				yuv[2]-=128;
#endif
			}
			for(int kc=0;kc<3;++kc)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				int16_t
					NNN	=rows[3][0+0*NCH*NROWS*NVAL],
					NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
					NNW	=rows[2][0-1*NCH*NROWS*NVAL],
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL],
					NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
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
					eNEE	=rows[1][1+2*NCH*NROWS*NVAL],
					eNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
					eW	=rows[0][1-1*NCH*NROWS*NVAL];
				int ctx=FLOOR_LOG2(eW*eW+1);
				int32_t *currhist;
				int64_t *currweights=weights+NPREDS*kc, p1=1LL<<L1SH>>1;
				int pred, curr, vmin, vmax, j, e;

				if(ctx>NCTX-1)
					ctx=NCTX-1;
				currhist=hist[kc][ctx];
#define PRED(EXPR) estim[j]=EXPR; p1+=currweights[j]*estim[j]; ++j;
				j=0;
				PREDLIST
#undef  PRED
				p1>>=L1SH;
				pred=(int)p1;
				vmax=N; vmin=W;
				if(N<W)vmin=N, vmax=W;
				if((uint32_t)(kx-2)>(uint32_t)(iw-2)||ky<2)//CG fallback
					pred=estim[0];
				else
				{
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
				}
				CLAMP2(pred, vmin, vmax);
				pred+=offset;
#ifdef SIGNED_PIXEL
				CLAMP2(pred, -128, 127);
#else
				CLAMP2(pred, 0, 255);
#endif
				den=currhist[NLEVELS]+NLEVELS;
				if(fwd)
				{
					int t;

					error=(int8_t)(yuv[kc]-pred);
					sym=error<<1^error>>31;

					for(t=0, cdf=0;;++t)
					{
						freq=currhist[t]+1;
						if(t>=sym)
							break;
						cdf+=freq;
					}
					if(range<=0xFFFF)
					{
						*(uint32_t*)streamptr=(uint32_t)(low>>32);
						streamptr+=sizeof(uint32_t);
#ifdef _DEBUG
						if(streamptr>NNptr)
							CRASH("");
#endif
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
#ifdef FIFOVAL
					valfifo_enqueue(freq<<16|cdf);
#endif
					low+=range*cdf/den;//DIV takes 31.55% E, 31.46% D
					range=range*freq/den-1;
				}
				else
				{
					int c, cdf2;
					
					if(range<=0xFFFF)
					{
						code=code<<32|*(unsigned*)streamptr;
						streamptr+=sizeof(uint32_t);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					c=(int)(((code-low)*den+den-1)/range);
					for(sym=0, cdf=0;;++sym)
					{
						freq=currhist[sym]+1;
						cdf2=cdf+freq;
						if(cdf2>c)
							break;
						cdf=cdf2;
					}
#ifdef FIFOVAL
					valfifo_check(freq<<16|cdf);
#endif
					low+=range*cdf/den;
					range=range*freq/den-1;

					error=sym>>1^-(sym&1);
#ifdef SIGNED_PIXEL
					yuv[kc]=(int8_t)(error+pred);
#else
					yuv[kc]=(uint8_t)(error+pred);
#endif
#ifdef ENABLE_GUIDE
					int perm[]={yidx, uidx, vidx};
#ifdef SIGNED_PIXEL
					int val=(uint8_t)(yuv[kc]+128);
#else
					int val=yuv[kc];
#endif
					if(g_image[ptr-image+perm[kc]]!=val)
						CRASH("");
#endif
				}
				++currhist[sym];
				++currhist[NLEVELS];
				if(currhist[NLEVELS]>=0xFFFF-2*NLEVELS)
				{
					den=0;
					for(int ks=0;ks<NLEVELS;++ks)
						den+=currhist[ks]>>=1;
					currhist[NLEVELS]=den;
				}
				curr=yuv[kc]-offset, e=(curr>(int32_t)p1)-(curr<(int32_t)p1);
#define PRED(EXPR) currweights[j]+=(int64_t)e*estim[j]; ++j;
				j=0;
				PREDLIST
#undef  PRED
				rows[0][0]=curr;
				rows[0][1]=(2*eW+(sym<<3)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
				offset=(kc?vc0*yuv[0]+vc1*yuv[1]:uc0*yuv[0])>>2;
			}
			if(!fwd)
			{
#ifdef SIGNED_PIXEL
				yuv[0]+=128;
				yuv[1]+=128;
				yuv[2]+=128;
#endif
				ptr[yidx]=yuv[0];
				ptr[uidx]=yuv[1];
				ptr[vidx]=yuv[2];
			}
			NNptr+=3;
			Nptr+=3;
			ptr+=3;
		}
	}
	_mm_free(pixels);
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing", dstfn);
			return 1;
		}
		if(fwd)
		{
			*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;//flush
			*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;

			fwrite("10", 1, 2, fdst);
			fwrite(&iw, 1, 3, fdst);
			fwrite(&ih, 1, 3, fdst);
			fwrite(&bestrct, 1, 1, fdst);
			fwrite(stream, 1, streamptr-stream, fdst);
#ifdef LOUD
			csize=streamptr-stream;
#endif
		}
		else
		{
			fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(image, 1, usize, fdst);
		}
		fclose(fdst);
	}
	free(buf);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		printf("WH %d*%d   RCT %2d %s  \"%s\"\n"
			, iw, ih
			, bestrct, rct_names[bestrct]
			, srcfn
		);
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
	(void)&time_sec;
	(void)&rct_names;
	(void)NNptr;
	(void)Nptr;
	return 0;
}
