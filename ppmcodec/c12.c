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
	#define ESTIMATE_BITSIZE
//	#define PRINT_STATS
	#define PRINT_RCT
//	#define PRINTBITS
//	#define PRINTBITS2
//	#define EXPERIMENT
//	#define PRINT_ERRORCOUNTS

	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif

//	#define ANALYSIS_NB
//	#define ANALYSIS_L1
//	#define ANALYSIS_BY_ENTROPY

	#define RCT_BASE8
//	#define RCT_BASE64

//	#define NEW_SYSTEM//X
//	#define SIMD_L1
//	#define USE_TABLES
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

#endif


#if 0
#define PREDLIST\
	PRED(N)\
	PRED(NNN)\
	PRED(NNNN)\
	PRED(N+yN)\
	PRED(N+yNE)\
	PRED(N+yNW)\
	PRED(N+yW)\
	PRED(NEEEE)\
	PRED(W)\
	PRED(WWW)\
	PRED(WWWW)\
	PRED(W+xW)\
	PRED(W+xNE)\
	PRED(W+xNEE)\
	PRED(W+xNW)\
	PRED(W+2*xW-xWW)\

#endif
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
	PRED((N+W)>>1)\
	PRED(N+W-NW)\
	PRED((3*(N+W)-2*NW)>>2)\
	PRED((4*(N+W)+NE-NW)>>3)\
	PRED(W+((5*(N-NW)+NE-WW)>>3))\
	PRED(W+((6*N-5*NW+NE-WW-NN)>>3))\
	PRED(W+((NNW-2*NN-NNE-NWW-9*NW+10*N+4*NE-2*WW)>>3))\

#endif
#if 0
#define PREDLIST\
	PRED(NNNN	)\
	PRED(NNN	)\
	PRED(NN		)\
	PRED(NNE	)\
	PRED(NW		)\
	PRED(N		)\
	PRED(NE		)\
	PRED(NEE	)\
	PRED(NEEE	)\
	PRED(NEEEE	)\
	PRED(WWWW	)\
	PRED(WWW	)\
	PRED(WW		)\
	PRED(W		)\
	PRED(xNW	)\
	PRED(xNE	)\
	PRED(xNEE	)\
	PRED(xNEEE	)\
	PRED(xWW	)\
	PRED(xW		)\
	PRED(yNN	)\
	PRED(yNW	)\
	PRED(yN		)\
	PRED(yNE	)\
	PRED(yW		)\

#endif
#if 0
	PRED(N+2*yN-yNN)\
	PRED(NEEE)\
	PRED(4*(N+NNN)-6*NN-NNNN)\
	PRED(4*(W+WWW)-6*WW-WWWW)\
	PRED(W+NEEE-NEE)\

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
	GRLIMIT=16,
#ifdef USE_COUNTERS
	PROBBITS_USE=12,

	CTRBITS=9,
	CTRMASK=(1<<CTRBITS)-1,
	
	PROBBITS_STORE=24,
	PROBSHIFT=PROBBITS_STORE-PROBBITS_USE,

//	HISTBITS=22,// <= 64-CTRBITS*2
//	HISTMASK=(1<<HISTBITS)-1,
//	CTRFBITS=(32-CTRBITS)>>1,
//	CTRFBITS=4,
//	CTRFMASK=(1<<CTRFBITS)-1,
#else
	PROBBITS_STORE=15,
	PROBBITS_USE=14,
//	PROBBITS_STORE=17,
//	PROBBITS_USE=12,
	PROBSHIFT=PROBBITS_STORE-PROBBITS_USE,
#endif
	XPAD=8,
	NCH=3,
	NROWS=4,
	NVAL=4,
//	NVAL0=2,
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
#if 0
ALIGN(32) static const int16_t cRCT_out_coeffs1[]=
{
	0,
	1, 2, 3, 4, 5, 6, 7, 8,
	1,
	1, 2,
	1, 2, 3,
	1, 2, 3, 4,
	1, 2, 3, 4, 5,
	1, 2, 3, 4, 5, 6,
	1, 2, 3, 4, 5, 6, 7,

	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};
ALIGN(32) static const int16_t cRCT_out_coeffs2[]=
{
	0,
	0, 0, 0, 0, 0, 0, 0, 0,
	1,
	2, 1,
	3, 2, 1,
	4, 3, 2, 1,
	5, 4, 3, 2, 1,
	6, 5, 4, 3, 2, 1,
	7, 6, 5, 4, 3, 2, 1,
	
	1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0,
};
#endif
#define OCHLIST\
	OCH(YX00) OCH(Y0X0) OCH(Y00X)\
	\
	OCH(CX10) OCH(C0X1) OCH(C10X)\
	OCH(CX20) OCH(C0X2) OCH(C20X)\
	OCH(CX30) OCH(C0X3) OCH(C30X)\
	OCH(CX40) OCH(C0X4) OCH(C40X)\
	OCH(CX50) OCH(C0X5) OCH(C50X)\
	OCH(CX60) OCH(C0X6) OCH(C60X)\
	OCH(CX70) OCH(C0X7) OCH(C70X)\
	OCH(CX80) OCH(C0X8) OCH(C80X)\
	\
	OCH(CX01) OCH(C1X0) OCH(C01X)\
	OCH(CX02) OCH(C2X0) OCH(C02X)\
	OCH(CX03) OCH(C3X0) OCH(C03X)\
	OCH(CX04) OCH(C4X0) OCH(C04X)\
	OCH(CX05) OCH(C5X0) OCH(C05X)\
	OCH(CX06) OCH(C6X0) OCH(C06X)\
	OCH(CX07) OCH(C7X0) OCH(C07X)\
	OCH(CX08) OCH(C8X0) OCH(C08X)\
	\
	OCH(CX11) OCH(C1X1) OCH(C11X)\
	\
	OCH(CX12) OCH(C2X1) OCH(C12X)\
	OCH(CX21) OCH(C1X2) OCH(C21X)\
	\
	OCH(CX13) OCH(C3X1) OCH(C13X)\
	OCH(CX22) OCH(C2X2) OCH(C22X)\
	OCH(CX31) OCH(C1X3) OCH(C31X)\
	\
	OCH(CX14) OCH(C4X1) OCH(C14X)\
	OCH(CX23) OCH(C3X2) OCH(C23X)\
	OCH(CX32) OCH(C2X3) OCH(C32X)\
	OCH(CX41) OCH(C1X4) OCH(C41X)\
	\
	OCH(CX15) OCH(C5X1) OCH(C15X)\
	OCH(CX24) OCH(C4X2) OCH(C24X)\
	OCH(CX33) OCH(C3X3) OCH(C33X)\
	OCH(CX42) OCH(C2X4) OCH(C42X)\
	OCH(CX51) OCH(C1X5) OCH(C51X)\
	\
	OCH(CX16) OCH(C6X1) OCH(C16X)\
	OCH(CX25) OCH(C5X2) OCH(C25X)\
	OCH(CX34) OCH(C4X3) OCH(C34X)\
	OCH(CX43) OCH(C3X4) OCH(C43X)\
	OCH(CX52) OCH(C2X5) OCH(C52X)\
	OCH(CX61) OCH(C1X6) OCH(C61X)\
	\
	OCH(CX17) OCH(C7X1) OCH(C17X)\
	OCH(CX26) OCH(C6X2) OCH(C26X)\
	OCH(CX35) OCH(C5X3) OCH(C35X)\
	OCH(CX44) OCH(C4X4) OCH(C44X)\
	OCH(CX53) OCH(C3X5) OCH(C53X)\
	OCH(CX62) OCH(C2X6) OCH(C62X)\
	OCH(CX71) OCH(C1X7) OCH(C71X)\

typedef enum _OCHIndex
{
#define OCH(X) OCH_##X,
	OCHLIST
#undef  OCH
	OCH_COUNT,
#if 0
	OCH_C1X0=OCH_CX10,
	OCH_C01X=OCH_C0X1,
	OCH_CX01=OCH_C10X,
	OCH_C2X0=OCH_CX20,
	OCH_C02X=OCH_C0X2,
	OCH_CX02=OCH_C20X,
	OCH_C3X0=OCH_CX30,
	OCH_C03X=OCH_C0X3,
	OCH_CX03=OCH_C30X,
	OCH_C4X0=OCH_CX40,
	OCH_C04X=OCH_C0X4,
	OCH_CX04=OCH_C40X,
	OCH_C5X0=OCH_CX50,
	OCH_C05X=OCH_C0X5,
	OCH_CX05=OCH_C50X,
	OCH_C6X0=OCH_CX60,
	OCH_C06X=OCH_C0X6,
	OCH_CX06=OCH_C60X,
	OCH_C7X0=OCH_CX70,
	OCH_C07X=OCH_C0X7,
	OCH_CX07=OCH_C70X,
	OCH_C8X0=OCH_CX80,
	OCH_C08X=OCH_C0X8,
	OCH_CX08=OCH_C80X,
	OCH_R=OCH_YX00,
	OCH_G=OCH_Y0X0,
	OCH_B=OCH_Y00X,
	OCH_BG=OCH_C08X,
	OCH_BR=OCH_C80X,
	OCH_RG=OCH_CX80,
	OCH_RB=OCH_CX08,
	OCH_GB=OCH_C0X8,
	OCH_GR=OCH_C8X0,
	OCH_R1=OCH_CX26,
	OCH_G1=OCH_C6X2,
	OCH_B1=OCH_C26X,
	OCH_R2=OCH_CX44,
	OCH_G2=OCH_C4X4,
	OCH_B2=OCH_C44X,
	OCH_R3=OCH_CX62,
	OCH_G3=OCH_C2X6,
	OCH_B3=OCH_C62X,
#endif
} OCHIndex;
static const char *och_names[]=
{
#define OCH(X) #X,
	OCHLIST
#undef  OCH
};
typedef struct _RCTInfo
{
	uint8_t och[3], perm[3], cu0, cv0, cv1;
//	uint8_t cy0, cy1, cu2, cv2;
} RCTInfo;
static void rct2str(char *str, int size, const RCTInfo *rct)//17 bytes
{
	const char rgb[]="RGB";

	snprintf(str, size, "%c_%c%03d_%c%03d_%03d"
		, rgb[rct->perm[0]]
		, rgb[rct->perm[1]]
		, rct->cu0
		, rgb[rct->perm[2]]
		, rct->cv0
		, rct->cv1
	);
}
static void print_rct2(const RCTInfo *rct)//17 bytes
{
	char str[20]={0};
	rct2str(str, sizeof(str)-1, rct);
	printf("%s", str);
}
static void rct1str(char *str, const RCTInfo *rct)//>=13 chars
{
	memset(str, '0', 12);
	str[0]='_';
	str[rct->perm[0]+1]='X';
	str[4]='_';
	str[rct->perm[0]+5]=rct->cu0+(rct->cu0<10?'0':'A'-10);
	str[rct->perm[1]+5]='X';
	str[8]='_';
	str[rct->perm[0]+9]=rct->cv0+(rct->cv0<10?'0':'A'-10);
	str[rct->perm[1]+9]=rct->cv1+(rct->cv1<10?'0':'A'-10);
	str[rct->perm[2]+9]='X';
	str[12]=0;
}
static void print_rct(const RCTInfo *rct)//up to 4-bit coeffs
{
	char info[20]={0};

	rct1str(info, rct);
	printf("RCT__%.3s_%.3s_%.3s"
		, info+0
		, info+4
		, info+8
	);
}
static int och_getidx(const char *label)
{
	if(label[1]=='X')
		return 0;
	if(label[2]=='X')
		return 1;
	if(label[3]=='X')
		return 2;
	CRASH("");
	return 0;
}
static void crct_get(RCTInfo *rct, int c0, int c1, int c2)
{
	const char *n0=och_names[c0];
	const char *n1=och_names[c1];
	const char *n2=och_names[c2];

	rct->perm[0]=och_getidx(n0);
	rct->perm[1]=och_getidx(n1);
	rct->perm[2]=och_getidx(n2);
	rct->cu0=n1[rct->perm[0]+1]-'0';
	rct->cv0=n2[rct->perm[0]+1]-'0';
	rct->cv1=n2[rct->perm[1]+1]-'0';
	if((uint32_t)rct->cu0>8||(uint32_t)rct->cv0>8||(uint32_t)rct->cv1>8)
		CRASH("");
}
#if 0
enum
{
	OCH0_COUNT=3,
	OCH1_COUNT=3*8,
	OCH2_COUNT=3*(7*(7+1)/2),
};
static void och_get(const char *och, int yidx, int uidx, int *ret_ch, int *ret_c0, int *ret_c1)
{
	if(och[1]=='X')
	{
		*ret_ch=0;
		*ret_c0=yidx==1?och[3]-'0':och[2]-'0';
		*ret_c1=uidx==1?och[3]-'0':och[2]-'0';
	}
	else if(och[2]=='X')
	{
		*ret_ch=1;
		*ret_c0=yidx==2?och[3]-'0':och[1]-'0';
		*ret_c1=uidx==2?och[3]-'0':och[1]-'0';
	}
	else if(och[3]=='X')
	{
		*ret_ch=2;
		*ret_c0=yidx==1?och[2]-'0':och[1]-'0';
		*ret_c1=uidx==1?och[2]-'0':och[1]-'0';
	}
	else
		CRASH("");
	if(yidx==*ret_ch||uidx==*ret_ch)
		CRASH("");
}
static int crct_isvalid(RCTInfo *rct)
{
	int valid=1;

	valid&=(uint32_t)rct->perm[0]<=3;
	valid&=(uint32_t)rct->perm[1]<=3;
	valid&=(uint32_t)rct->perm[2]<=3;
	valid&=rct->perm[0]==rct->perm[1];
	valid&=rct->perm[1]==rct->perm[2];
	valid&=rct->perm[2]==rct->perm[0];
	return valid;
}
static void crct_getparams(RCTInfo *rct, int k0, int k1, int k2)
{
	static const uint8_t coeffs[]=
	{
		1, 1,

		1, 2,
		2, 1,

		1, 3,
		2, 2,
		3, 1,

		1, 4,
		2, 3,
		3, 2,
		4, 1,

		1, 5,
		2, 4,
		3, 3,
		4, 2,
		5, 1,

		1, 6,
		2, 5,
		3, 4,
		4, 3,
		5, 2,
		6, 1,

		1, 7,
		2, 6,
		3, 5,
		4, 4,
		5, 3,
		6, 2,
		7, 1,
	};
	int yidx=k0, uidx=k1%3, vidx=k2%3;
	int cu0=k1/3;
	int cv0=0, cv1=0;

	if(k2>=OCH0_COUNT)
	{
		if(k2>=OCH0_COUNT+OCH1_COUNT)
		{
			int idx=(k2-(OCH0_COUNT+OCH1_COUNT))/3+1;
			cv1=coeffs[2*idx+0];
			cv0=coeffs[2*idx+1];
		}
		else
			cv1=(k2-OCH0_COUNT)/3+1;
		if((yidx+2)%3==vidx)
		{
			int tmp=cv1;
			cv1=cv0;
			cv0=tmp;
		}
	}
	rct->perm[0]=yidx;
	rct->perm[1]=uidx;
	rct->perm[2]=vidx;
	rct->cu0=cu0;
	rct->cv0=cv0;
	rct->cv1=cv1;
#if 0
	if(crct_isvalid(rct))
	{
		char str[20]={0}, truth[20]={0};
		int printed=0;

		rct1str(str, rct);
		printed+=snprintf(truth+printed, sizeof(truth)-1-printed, "_%.3s_%.3s_%.3s"
			, och_names[k0]+1
			, och_names[k1]+1
			, och_names[k2]+1
		);
		if(memcmp(str, truth, 9))
			CRASH("");
	}
#endif
}
//YUV = RCT * RGB	example: _X00_80X_6X2 == [1 0 0; -1 0 1; -3/4 1 -1/4]
#define RCTLIST\
	RCT(_X00_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_X00_0X0_08X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 8)\
	RCT(_X00_0X0_80X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  8, 0)\
	RCT(_0X0_00X_X80,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  8, 0)\
	RCT(_0X0_00X_X08,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 8)\
	RCT(_00X_X00_8X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 8)\
	RCT(_00X_X00_0X8,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  8, 0)\
	RCT(_0X0_08X_X80,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	8,  8, 0)\
	RCT(_0X0_08X_X08,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	8,  0, 8)\
	RCT(_0X0_X80_80X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	8,  0, 8)\
	RCT(_00X_X08_0X8,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	8,  8, 0)\
	RCT(_00X_X08_8X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	8,  0, 8)\
	RCT(_00X_0X8_X80,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	8,  0, 8)\
	RCT(_X00_8X0_80X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	8,  8, 0)\
	RCT(_X00_8X0_08X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	8,  0, 8)\
	RCT(_X00_80X_0X8,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	8,  0, 8)\
	RCT(_X00_0X0_26X,	OCH_R,		OCH_G,		OCH_B1,		0, 1, 2,	0,  2, 6)\
	RCT(_X00_8X0_26X,	OCH_R,		OCH_GR,		OCH_B1,		0, 1, 2,	8,  2, 6)\
	RCT(_X00_00X_6X2,	OCH_R,		OCH_B,		OCH_G1,		0, 2, 1,	0,  6, 2)\
	RCT(_X00_80X_6X2,	OCH_R,		OCH_BR,		OCH_G1,		0, 2, 1,	8,  6, 2)\
	RCT(_0X0_00X_X26,	OCH_G,		OCH_B,		OCH_R1,		1, 2, 0,	0,  2, 6)\
	RCT(_0X0_08X_X26,	OCH_G,		OCH_BG,		OCH_R1,		1, 2, 0,	8,  2, 6)\
	RCT(_0X0_X80_26X,	OCH_G,		OCH_RG,		OCH_B1,		1, 0, 2,	8,  6, 2)\
	RCT(_00X_X08_6X2,	OCH_B,		OCH_RB,		OCH_G1,		2, 0, 1,	8,  2, 6)\
	RCT(_00X_08X_X26,	OCH_B,		OCH_GB,		OCH_R1,		2, 1, 0,	8,  6, 2)\
	RCT(_X00_0X0_44X,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,  4, 4)\
	RCT(_X00_8X0_44X,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	8,  4, 4)\
	RCT(_X00_00X_4X4,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,  4, 4)\
	RCT(_X00_80X_4X4,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	8,  4, 4)\
	RCT(_0X0_00X_X44,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,  4, 4)\
	RCT(_0X0_08X_X44,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	8,  4, 4)\
	RCT(_0X0_X80_44X,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	8,  4, 4)\
	RCT(_00X_X08_4X4,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	8,  4, 4)\
	RCT(_00X_0X8_X44,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	8,  4, 4)\
	RCT(_X00_0X0_62X,	OCH_R,		OCH_G,		OCH_B3,		0, 1, 2,	0,  6, 2)\
	RCT(_X00_8X0_62X,	OCH_R,		OCH_GR,		OCH_B3,		0, 1, 2,	8,  6, 2)\
	RCT(_X00_00X_2X6,	OCH_R,		OCH_B,		OCH_G3,		0, 2, 1,	0,  2, 6)\
	RCT(_X00_80X_2X6,	OCH_R,		OCH_BR,		OCH_G3,		0, 2, 1,	8,  2, 6)\
	RCT(_0X0_00X_X62,	OCH_G,		OCH_B,		OCH_R3,		1, 2, 0,	0,  6, 2)\
	RCT(_0X0_08X_X62,	OCH_G,		OCH_BG,		OCH_R3,		1, 2, 0,	8,  6, 2)\
	RCT(_0X0_X80_62X,	OCH_G,		OCH_RG,		OCH_B3,		1, 0, 2,	8,  2, 6)\
	RCT(_00X_X08_2X6,	OCH_B,		OCH_RB,		OCH_G3,		2, 0, 1,	8,  6, 2)\
	RCT(_00X_0X8_X62,	OCH_B,		OCH_GB,		OCH_R3,		2, 1, 0,	8,  2, 6)\
	RCT(_X00_0X0_17X,	OCH_R,		OCH_G,		OCH_C17X,	0, 1, 2,	0,  1, 7)\
	RCT(_X00_8X0_17X,	OCH_R,		OCH_GR,		OCH_C17X,	0, 1, 2,	8,  1, 7)\
	RCT(_X00_00X_7X1,	OCH_R,		OCH_B,		OCH_C7X1,	0, 2, 1,	0,  7, 1)\
	RCT(_X00_80X_7X1,	OCH_R,		OCH_BR,		OCH_C7X1,	0, 2, 1,	8,  7, 1)\
	RCT(_0X0_00X_X17,	OCH_G,		OCH_B,		OCH_CX17,	1, 2, 0,	0,  1, 7)\
	RCT(_0X0_08X_X17,	OCH_G,		OCH_BG,		OCH_CX17,	1, 2, 0,	8,  1, 7)\
	RCT(_0X0_X80_17X,	OCH_G,		OCH_RG,		OCH_C17X,	1, 0, 2,	8,  7, 1)\
	RCT(_00X_X08_7X1,	OCH_B,		OCH_RB,		OCH_C7X1,	2, 0, 1,	8,  1, 7)\
	RCT(_00X_08X_X17,	OCH_B,		OCH_GB,		OCH_CX17,	2, 1, 0,	8,  7, 1)\
	RCT(_X00_0X0_71X,	OCH_R,		OCH_G,		OCH_C71X,	0, 1, 2,	0,  7, 1)\
	RCT(_X00_8X0_71X,	OCH_R,		OCH_GR,		OCH_C71X,	0, 1, 2,	8,  7, 1)\
	RCT(_X00_00X_1X7,	OCH_R,		OCH_B,		OCH_C1X7,	0, 2, 1,	0,  1, 7)\
	RCT(_X00_80X_1X7,	OCH_R,		OCH_BR,		OCH_C1X7,	0, 2, 1,	8,  1, 7)\
	RCT(_0X0_00X_X71,	OCH_G,		OCH_B,		OCH_CX71,	1, 2, 0,	0,  7, 1)\
	RCT(_0X0_08X_X71,	OCH_G,		OCH_BG,		OCH_CX71,	1, 2, 0,	8,  7, 1)\
	RCT(_0X0_X80_71X,	OCH_G,		OCH_RG,		OCH_C71X,	1, 0, 2,	8,  1, 7)\
	RCT(_00X_X08_1X7,	OCH_B,		OCH_RB,		OCH_C1X7,	2, 0, 1,	8,  7, 1)\
	RCT(_00X_0X8_X71,	OCH_B,		OCH_GB,		OCH_CX71,	2, 1, 0,	8,  1, 7)\
	RCT(_X00_0X0_33X,	OCH_R,		OCH_G,		OCH_C33X,	0, 1, 2,	0,  3, 3)\
	RCT(_X00_8X0_33X,	OCH_R,		OCH_GR,		OCH_C33X,	0, 1, 2,	8,  3, 3)\
	RCT(_X00_00X_3X3,	OCH_R,		OCH_B,		OCH_C3X3,	0, 2, 1,	0,  3, 3)\
	RCT(_X00_80X_3X3,	OCH_R,		OCH_BR,		OCH_C3X3,	0, 2, 1,	8,  3, 3)\
	RCT(_0X0_00X_X33,	OCH_G,		OCH_B,		OCH_CX33,	1, 2, 0,	0,  3, 3)\
	RCT(_0X0_08X_X33,	OCH_G,		OCH_BG,		OCH_CX33,	1, 2, 0,	8,  3, 3)\
	RCT(_0X0_X80_33X,	OCH_G,		OCH_RG,		OCH_C33X,	1, 0, 2,	8,  3, 3)\
	RCT(_00X_X08_3X3,	OCH_B,		OCH_RB,		OCH_C3X3,	2, 0, 1,	8,  3, 3)\
	RCT(_00X_08X_X33,	OCH_B,		OCH_GB,		OCH_CX33,	2, 1, 0,	8,  3, 3)\
	RCT(_X00_0X0_35X,	OCH_R,		OCH_G,		OCH_C35X,	0, 1, 2,	0,  3, 5)\
	RCT(_X00_8X0_35X,	OCH_R,		OCH_GR,		OCH_C35X,	0, 1, 2,	8,  3, 5)\
	RCT(_X00_00X_5X3,	OCH_R,		OCH_B,		OCH_C5X3,	0, 2, 1,	0,  5, 3)\
	RCT(_X00_80X_5X3,	OCH_R,		OCH_BR,		OCH_C5X3,	0, 2, 1,	8,  5, 3)\
	RCT(_0X0_00X_X35,	OCH_G,		OCH_B,		OCH_CX35,	1, 2, 0,	0,  3, 5)\
	RCT(_0X0_08X_X35,	OCH_G,		OCH_BG,		OCH_CX35,	1, 2, 0,	8,  3, 5)\
	RCT(_0X0_X80_35X,	OCH_G,		OCH_RG,		OCH_C35X,	1, 0, 2,	8,  5, 3)\
	RCT(_00X_X08_5X3,	OCH_B,		OCH_RB,		OCH_C5X3,	2, 0, 1,	8,  3, 5)\
	RCT(_00X_08X_X35,	OCH_B,		OCH_GB,		OCH_CX35,	2, 1, 0,	8,  5, 3)\
	RCT(_X00_0X0_53X,	OCH_R,		OCH_G,		OCH_C53X,	0, 1, 2,	0,  5, 3)\
	RCT(_X00_8X0_53X,	OCH_R,		OCH_GR,		OCH_C53X,	0, 1, 2,	8,  5, 3)\
	RCT(_X00_00X_3X5,	OCH_R,		OCH_B,		OCH_C3X5,	0, 2, 1,	0,  3, 5)\
	RCT(_X00_80X_3X5,	OCH_R,		OCH_BR,		OCH_C3X5,	0, 2, 1,	8,  3, 5)\
	RCT(_0X0_00X_X53,	OCH_G,		OCH_B,		OCH_CX53,	1, 2, 0,	0,  5, 3)\
	RCT(_0X0_08X_X53,	OCH_G,		OCH_BG,		OCH_CX53,	1, 2, 0,	8,  5, 3)\
	RCT(_0X0_X80_53X,	OCH_G,		OCH_RG,		OCH_C53X,	1, 0, 2,	8,  3, 5)\
	RCT(_00X_X08_3X5,	OCH_B,		OCH_RB,		OCH_C3X5,	2, 0, 1,	8,  5, 3)\
	RCT(_00X_0X8_X53,	OCH_B,		OCH_GB,		OCH_CX53,	2, 1, 0,	8,  3, 5)\
	RCT(_X00_0X0_06X,	OCH_R,		OCH_G,		OCH_C06X,	0, 1, 2,	0,  0, 6)\
	RCT(_X00_8X0_06X,	OCH_R,		OCH_GR,		OCH_C06X,	0, 1, 2,	8,  0, 6)\
	RCT(_X00_00X_6X0,	OCH_R,		OCH_B,		OCH_C6X0,	0, 2, 1,	0,  6, 0)\
	RCT(_X00_80X_6X0,	OCH_R,		OCH_BR,		OCH_C6X0,	0, 2, 1,	8,  6, 0)\
	RCT(_0X0_00X_X06,	OCH_G,		OCH_B,		OCH_CX06,	1, 2, 0,	0,  0, 6)\
	RCT(_0X0_08X_X06,	OCH_G,		OCH_BG,		OCH_CX06,	1, 2, 0,	8,  0, 6)\
	RCT(_0X0_X80_06X,	OCH_G,		OCH_RG,		OCH_C06X,	1, 0, 2,	8,  6, 0)\
	RCT(_00X_X08_6X0,	OCH_B,		OCH_RB,		OCH_C6X0,	2, 0, 1,	8,  0, 6)\
	RCT(_00X_08X_X06,	OCH_B,		OCH_GB,		OCH_CX06,	2, 1, 0,	8,  6, 0)\
	RCT(_X00_0X0_60X,	OCH_R,		OCH_G,		OCH_C60X,	0, 1, 2,	0,  6, 0)\
	RCT(_X00_8X0_60X,	OCH_R,		OCH_GR,		OCH_C60X,	0, 1, 2,	8,  6, 0)\
	RCT(_X00_00X_0X6,	OCH_R,		OCH_B,		OCH_C0X6,	0, 2, 1,	0,  0, 6)\
	RCT(_X00_80X_0X6,	OCH_R,		OCH_BR,		OCH_C0X6,	0, 2, 1,	8,  0, 6)\
	RCT(_0X0_00X_X60,	OCH_G,		OCH_B,		OCH_CX60,	1, 2, 0,	0,  6, 0)\
	RCT(_0X0_08X_X60,	OCH_G,		OCH_BG,		OCH_CX60,	1, 2, 0,	8,  6, 0)\
	RCT(_0X0_X80_60X,	OCH_G,		OCH_RG,		OCH_C60X,	1, 0, 2,	8,  0, 6)\
	RCT(_00X_X08_0X6,	OCH_B,		OCH_RB,		OCH_C0X6,	2, 0, 1,	8,  6, 0)\
	RCT(_00X_0X8_X60,	OCH_B,		OCH_GB,		OCH_CX60,	2, 1, 0,	8,  0, 6)\
	RCT(_X00_0X0_15X,	OCH_R,		OCH_G,		OCH_C15X,	0, 1, 2,	0,  1, 5)\
	RCT(_X00_8X0_15X,	OCH_R,		OCH_GR,		OCH_C15X,	0, 1, 2,	8,  1, 5)\
	RCT(_X00_00X_5X1,	OCH_R,		OCH_B,		OCH_C5X1,	0, 2, 1,	0,  5, 1)\
	RCT(_X00_80X_5X1,	OCH_R,		OCH_BR,		OCH_C5X1,	0, 2, 1,	8,  5, 1)\
	RCT(_0X0_00X_X15,	OCH_G,		OCH_B,		OCH_CX15,	1, 2, 0,	0,  1, 5)\
	RCT(_0X0_08X_X15,	OCH_G,		OCH_BG,		OCH_CX15,	1, 2, 0,	8,  1, 5)\
	RCT(_0X0_X80_15X,	OCH_G,		OCH_RG,		OCH_C15X,	1, 0, 2,	8,  5, 1)\
	RCT(_00X_X08_5X1,	OCH_B,		OCH_RB,		OCH_C5X1,	2, 0, 1,	8,  1, 5)\
	RCT(_00X_08X_X15,	OCH_B,		OCH_GB,		OCH_CX15,	2, 1, 0,	8,  5, 1)\
	RCT(_X00_0X0_51X,	OCH_R,		OCH_G,		OCH_C51X,	0, 1, 2,	0,  5, 1)\
	RCT(_X00_8X0_51X,	OCH_R,		OCH_GR,		OCH_C51X,	0, 1, 2,	8,  5, 1)\
	RCT(_X00_00X_1X5,	OCH_R,		OCH_B,		OCH_C1X5,	0, 2, 1,	0,  1, 5)\
	RCT(_X00_80X_1X5,	OCH_R,		OCH_BR,		OCH_C1X5,	0, 2, 1,	8,  1, 5)\
	RCT(_0X0_00X_X51,	OCH_G,		OCH_B,		OCH_CX51,	1, 2, 0,	0,  5, 1)\
	RCT(_0X0_08X_X51,	OCH_G,		OCH_BG,		OCH_CX51,	1, 2, 0,	8,  5, 1)\
	RCT(_0X0_X80_51X,	OCH_G,		OCH_RG,		OCH_C51X,	1, 0, 2,	8,  1, 5)\
	RCT(_00X_X08_1X5,	OCH_B,		OCH_RB,		OCH_C1X5,	2, 0, 1,	8,  5, 1)\
	RCT(_00X_0X8_X51,	OCH_B,		OCH_GB,		OCH_CX51,	2, 1, 0,	8,  1, 5)\
	RCT(_X00_0X0_24X,	OCH_R,		OCH_G,		OCH_C24X,	0, 1, 2,	0,  2, 4)\
	RCT(_X00_8X0_24X,	OCH_R,		OCH_GR,		OCH_C24X,	0, 1, 2,	8,  2, 4)\
	RCT(_X00_00X_4X2,	OCH_R,		OCH_B,		OCH_C4X2,	0, 2, 1,	0,  4, 2)\
	RCT(_X00_80X_4X2,	OCH_R,		OCH_BR,		OCH_C4X2,	0, 2, 1,	8,  4, 2)\
	RCT(_0X0_00X_X24,	OCH_G,		OCH_B,		OCH_CX24,	1, 2, 0,	0,  2, 4)\
	RCT(_0X0_08X_X24,	OCH_G,		OCH_BG,		OCH_CX24,	1, 2, 0,	8,  2, 4)\
	RCT(_0X0_X80_24X,	OCH_G,		OCH_RG,		OCH_C24X,	1, 0, 2,	8,  4, 2)\
	RCT(_00X_X08_4X2,	OCH_B,		OCH_RB,		OCH_C4X2,	2, 0, 1,	8,  2, 4)\
	RCT(_00X_08X_X24,	OCH_B,		OCH_GB,		OCH_CX24,	2, 1, 0,	8,  4, 2)\
	RCT(_X00_0X0_42X,	OCH_R,		OCH_G,		OCH_C42X,	0, 1, 2,	0,  4, 2)\
	RCT(_X00_8X0_42X,	OCH_R,		OCH_GR,		OCH_C42X,	0, 1, 2,	8,  4, 2)\
	RCT(_X00_00X_2X4,	OCH_R,		OCH_B,		OCH_C2X4,	0, 2, 1,	0,  2, 4)\
	RCT(_X00_80X_2X4,	OCH_R,		OCH_BR,		OCH_C2X4,	0, 2, 1,	8,  2, 4)\
	RCT(_0X0_00X_X42,	OCH_G,		OCH_B,		OCH_CX42,	1, 2, 0,	0,  4, 2)\
	RCT(_0X0_08X_X42,	OCH_G,		OCH_BG,		OCH_CX42,	1, 2, 0,	8,  4, 2)\
	RCT(_0X0_X80_42X,	OCH_G,		OCH_RG,		OCH_C42X,	1, 0, 2,	8,  2, 4)\
	RCT(_00X_X08_2X4,	OCH_B,		OCH_RB,		OCH_C2X4,	2, 0, 1,	8,  4, 2)\
	RCT(_00X_0X8_X42,	OCH_B,		OCH_GB,		OCH_CX42,	2, 1, 0,	8,  2, 4)\
	RCT(_X00_0X0_07X,	OCH_R,		OCH_G,		OCH_C07X,	0, 1, 2,	0,  0, 7)\
	RCT(_X00_8X0_07X,	OCH_R,		OCH_GR,		OCH_C07X,	0, 1, 2,	8,  0, 7)\
	RCT(_X00_00X_7X0,	OCH_R,		OCH_B,		OCH_C7X0,	0, 2, 1,	0,  7, 0)\
	RCT(_X00_80X_7X0,	OCH_R,		OCH_BR,		OCH_C7X0,	0, 2, 1,	8,  7, 0)\
	RCT(_0X0_00X_X07,	OCH_G,		OCH_B,		OCH_CX07,	1, 2, 0,	0,  0, 7)\
	RCT(_0X0_08X_X07,	OCH_G,		OCH_BG,		OCH_CX07,	1, 2, 0,	8,  0, 7)\
	RCT(_0X0_X80_07X,	OCH_G,		OCH_RG,		OCH_C07X,	1, 0, 2,	8,  7, 0)\
	RCT(_00X_X08_7X0,	OCH_B,		OCH_RB,		OCH_C7X0,	2, 0, 1,	8,  0, 7)\
	RCT(_00X_08X_X07,	OCH_B,		OCH_GB,		OCH_CX07,	2, 1, 0,	8,  7, 0)\
	RCT(_X00_0X0_70X,	OCH_R,		OCH_G,		OCH_C70X,	0, 1, 2,	0,  7, 0)\
	RCT(_X00_8X0_70X,	OCH_R,		OCH_GR,		OCH_C70X,	0, 1, 2,	8,  7, 0)\
	RCT(_X00_00X_0X7,	OCH_R,		OCH_B,		OCH_C0X7,	0, 2, 1,	0,  0, 7)\
	RCT(_X00_80X_0X7,	OCH_R,		OCH_BR,		OCH_C0X7,	0, 2, 1,	8,  0, 7)\
	RCT(_0X0_00X_X70,	OCH_G,		OCH_B,		OCH_CX70,	1, 2, 0,	0,  7, 0)\
	RCT(_0X0_08X_X70,	OCH_G,		OCH_BG,		OCH_CX70,	1, 2, 0,	8,  7, 0)\
	RCT(_0X0_X80_70X,	OCH_G,		OCH_RG,		OCH_C70X,	1, 0, 2,	8,  0, 7)\
	RCT(_00X_X08_0X7,	OCH_B,		OCH_RB,		OCH_C0X7,	2, 0, 1,	8,  7, 0)\
	RCT(_00X_0X8_X70,	OCH_B,		OCH_GB,		OCH_CX70,	2, 1, 0,	8,  0, 7)\
	RCT(_X00_0X0_05X,	OCH_R,		OCH_G,		OCH_C05X,	0, 1, 2,	0,  0, 5)\
	RCT(_X00_8X0_05X,	OCH_R,		OCH_GR,		OCH_C05X,	0, 1, 2,	8,  0, 5)\
	RCT(_X00_00X_5X0,	OCH_R,		OCH_B,		OCH_C5X0,	0, 2, 1,	0,  5, 0)\
	RCT(_X00_80X_5X0,	OCH_R,		OCH_BR,		OCH_C5X0,	0, 2, 1,	8,  5, 0)\
	RCT(_0X0_00X_X05,	OCH_G,		OCH_B,		OCH_CX05,	1, 2, 0,	0,  0, 5)\
	RCT(_0X0_08X_X05,	OCH_G,		OCH_BG,		OCH_CX05,	1, 2, 0,	8,  0, 5)\
	RCT(_0X0_X80_05X,	OCH_G,		OCH_RG,		OCH_C05X,	1, 0, 2,	8,  5, 0)\
	RCT(_00X_X08_5X0,	OCH_B,		OCH_RB,		OCH_C5X0,	2, 0, 1,	8,  0, 5)\
	RCT(_00X_08X_X05,	OCH_B,		OCH_GB,		OCH_CX05,	2, 1, 0,	8,  5, 0)\
	RCT(_X00_0X0_50X,	OCH_R,		OCH_G,		OCH_C50X,	0, 1, 2,	0,  5, 0)\
	RCT(_X00_8X0_50X,	OCH_R,		OCH_GR,		OCH_C50X,	0, 1, 2,	8,  5, 0)\
	RCT(_X00_00X_0X5,	OCH_R,		OCH_B,		OCH_C0X5,	0, 2, 1,	0,  0, 5)\
	RCT(_X00_80X_0X5,	OCH_R,		OCH_BR,		OCH_C0X5,	0, 2, 1,	8,  0, 5)\
	RCT(_0X0_00X_X50,	OCH_G,		OCH_B,		OCH_CX50,	1, 2, 0,	0,  5, 0)\
	RCT(_0X0_08X_X50,	OCH_G,		OCH_BG,		OCH_CX50,	1, 2, 0,	8,  5, 0)\
	RCT(_0X0_X80_50X,	OCH_G,		OCH_RG,		OCH_C50X,	1, 0, 2,	8,  0, 5)\
	RCT(_00X_X08_0X5,	OCH_B,		OCH_RB,		OCH_C0X5,	2, 0, 1,	8,  5, 0)\
	RCT(_00X_0X8_X50,	OCH_B,		OCH_GB,		OCH_CX50,	2, 1, 0,	8,  0, 5)\

#if 0
	RCT(_211_4X0_80X,	OCH_Y211,	OCH_C8X0,	OCH_C80X,	0, 1, 2,	4,  4, 0,	1, 1, 0, 0)\
	RCT(_211_4X0_31X,	OCH_Y211,	OCH_C8X0,	OCH_C62X,	0, 1, 2,	4,  4, 0,	1, 1, 0, 1)\
	RCT(_211_3X1_80X,	OCH_Y211,	OCH_C6X2,	OCH_C80X,	0, 1, 2,	4,  4, 0,	1, 1, 1, 0)\
	RCT(_310_4X0_80X,	OCH_Y310,	OCH_C8X0,	OCH_C80X,	0, 1, 2,	4,  4, 0,	1, 0, 0, 0)\
	RCT(_310_4X0_31X,	OCH_Y310,	OCH_C8X0,	OCH_C62X,	0, 1, 2,	4,  4, 0,	1, 0, 0, 1)\
	RCT(_310_3X1_80X,	OCH_Y310,	OCH_C6X2,	OCH_C80X,	0, 1, 2,	4,  4, 0,	1, 0, 1, 0)\
	RCT(_301_4X0_80X,	OCH_Y301,	OCH_C8X0,	OCH_C80X,	0, 1, 2,	4,  4, 0,	0, 1, 0, 0)\
	RCT(_301_4X0_31X,	OCH_Y301,	OCH_C8X0,	OCH_C62X,	0, 1, 2,	4,  4, 0,	0, 1, 0, 1)\
	RCT(_301_3X1_80X,	OCH_Y301,	OCH_C6X2,	OCH_C80X,	0, 1, 2,	4,  4, 0,	0, 1, 1, 0)\
	RCT(_121_04X_X80,	OCH_Y121,	OCH_C08X,	OCH_CX80,	1, 2, 0,	4,  4, 0,	1, 1, 0, 0)\
	RCT(_121_04X_X31,	OCH_Y121,	OCH_C08X,	OCH_CX62,	1, 2, 0,	4,  4, 0,	1, 1, 0, 1)\
	RCT(_121_13X_X80,	OCH_Y121,	OCH_C26X,	OCH_CX80,	1, 2, 0,	4,  4, 0,	1, 1, 1, 0)\
	RCT(_031_04X_X80,	OCH_Y031,	OCH_C08X,	OCH_CX80,	1, 2, 0,	4,  4, 0,	0, 1, 0, 0)\
	RCT(_031_04X_X31,	OCH_Y031,	OCH_C08X,	OCH_CX62,	1, 2, 0,	4,  4, 0,	0, 1, 0, 1)\
	RCT(_031_13X_X80,	OCH_Y031,	OCH_C26X,	OCH_CX80,	1, 2, 0,	4,  4, 0,	0, 1, 1, 0)\
	RCT(_130_80X_X80,	OCH_Y130,	OCH_C08X,	OCH_CX80,	1, 2, 0,	4,  4, 0,	1, 0, 0, 0)\
	RCT(_130_80X_X31,	OCH_Y130,	OCH_C08X,	OCH_CX62,	1, 2, 0,	4,  4, 0,	1, 0, 0, 1)\
	RCT(_130_31X_X80,	OCH_Y130,	OCH_C26X,	OCH_CX80,	1, 2, 0,	4,  4, 0,	1, 0, 1, 0)\
	RCT(_112_X04_0X8,	OCH_Y112,	OCH_CX08,	OCH_C0X8,	2, 0, 1,	4,  4, 0,	1, 1, 0, 0)\
	RCT(_112_X04_1X3,	OCH_Y112,	OCH_CX08,	OCH_C2X6,	2, 0, 1,	4,  4, 0,	1, 1, 0, 1)\
	RCT(_112_X13_0X8,	OCH_Y112,	OCH_CX26,	OCH_C0X8,	2, 0, 1,	4,  4, 0,	1, 1, 1, 0)\
	RCT(_013_X04_0X8,	OCH_Y013,	OCH_CX08,	OCH_C0X8,	2, 0, 1,	4,  4, 0,	0, 1, 0, 0)\
	RCT(_013_X04_1X3,	OCH_Y013,	OCH_CX08,	OCH_C2X6,	2, 0, 1,	4,  4, 0,	0, 1, 0, 1)\
	RCT(_013_X13_0X8,	OCH_Y013,	OCH_CX26,	OCH_C0X8,	2, 0, 1,	4,  4, 0,	0, 1, 1, 0)\
	RCT(_103_X80_0X8,	OCH_Y103,	OCH_CX08,	OCH_C0X8,	2, 0, 1,	4,  4, 0,	1, 0, 0, 0)\
	RCT(_103_X80_1X3,	OCH_Y103,	OCH_CX08,	OCH_C2X6,	2, 0, 1,	4,  4, 0,	1, 0, 0, 1)\
	RCT(_103_X31_0X8,	OCH_Y103,	OCH_CX26,	OCH_C0X8,	2, 0, 1,	4,  4, 0,	1, 0, 1, 0)
#endif
typedef enum _RCTIndex
{
#define RCT(LABEL, ...) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const RCTInfo rct_combinations[RCT_COUNT]=
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
#endif

#ifdef USE_COUNTERS
static uint64_t ctrtable[1<<CTRBITS*2];//{p1, next_0(n0, n1}, next_1(n0, n1)}
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
static uint32_t unary_count=0, binary_count=0;
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
//#ifdef _MSC_VER
//static uint64_t unary_count=0, binary_count=0;
//#endif
#if 0
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
#endif
//void inflation(void)
//{
//#ifdef _MSC_VER
//	CRASH("inflation");
//#endif
//	exit(1);
//}
#ifdef EXPERIMENT
static uint16_t ehists[3][NCTX][257];
static uint32_t estats[3][NCTX][256];
static uint32_t erice[3][NCTX][GRLIMIT];
static uint32_t ectr[3][NCTX][GRLIMIT];
static uint16_t er2[3][NCTX][GRLIMIT+2];
static uint64_t ehistory[3][NCTX][GRLIMIT*4];
static void experiment(uint8_t *image, int iw, int ih, int rct)
{
	enum
	{
		EPROBBITS_STORE=24,
		EPROBBITS_USE=12,
		EPROBSHIFT=EPROBBITS_STORE-EPROBBITS_USE,
	};
	int psize;
	int16_t *pixels;
	uint8_t *imptr=image;
	int yidx=rct_combinations[rct][II_PERM_Y];
	int uidx=rct_combinations[rct][II_PERM_U];
	int vidx=rct_combinations[rct][II_PERM_V];
	int cu0=rct_combinations[rct][II_COEFF_U_SUB_Y];
	int cv0=rct_combinations[rct][II_COEFF_V_SUB_Y];
	int cv1=rct_combinations[rct][II_COEFF_V_SUB_U];
	int32_t coeffs[3][L1NPREDS], bias[3];
	double hsize=0, ssize=0, rsize=0, csize=0, r2size=0;
	
	memset(ehists, 0, sizeof(ehists));
	memset(estats, 0, sizeof(estats));
	memset(erice, 0, sizeof(erice));
	memset(ectr, 0, sizeof(ectr));
//	FILLMEM_S((uint32_t*)ectr, 0x00800080, sizeof(ectr), sizeof(uint32_t));
	memset(er2, 0, sizeof(er2));
	memset(ehistory, 0, sizeof(ehistory));

	FILLMEM_S((int32_t*)coeffs, (1<<L1SH)/L1NPREDS, sizeof(coeffs), sizeof(int32_t));
	bias[0]=1<<L1SH>>1;
	bias[1]=1<<L1SH>>1;
	bias[2]=1<<L1SH>>1;
	psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error\n");
		return;
	}
	memset(pixels, 0, psize);
	for(int ky=0;ky<ih;++ky)
	{
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=0;kx<iw;++kx, imptr+=3)
		{
			int yuv[]=
			{
				imptr[yidx],
				imptr[uidx],
				imptr[vidx],
			};
			int offset=0;
			for(int kc=0;kc<3;++kc)
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
				ALIGN(32) int32_t estim[L1NPREDS];

				int32_t pred=0;
				int j, vmax, vmin, e, pred0, curr, error;
				int ctx=31-_lzcnt_u32(eW*eW+1);

				if(ctx>NCTX-1)
					ctx=NCTX-1;
				{
					pred=bias[kc];
#define PRED(EXPR) estim[j++]=EXPR;
					j=0;
					PREDLIST;
#undef  PRED
					__m256i mc0=_mm256_load_si256((__m256i*)coeffs[kc]+0);
					__m256i mc1=_mm256_load_si256((__m256i*)coeffs[kc]+1);
					__m256i mp0=_mm256_load_si256((__m256i*)estim+0);
					__m256i mp1=_mm256_load_si256((__m256i*)estim+1);

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
				}
				pred>>=L1SH;
				pred0=pred;
				vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);
				pred+=offset;
				CLAMP2(pred, 0, 255);

				curr=yuv[kc];

				error=(int8_t)(curr-pred);

				uint16_t *currhist=ehists[kc][ctx];
				int32_t *currstats=estats[kc][ctx];
				int32_t *currrice=erice[kc][ctx];
				int32_t *currctr=ectr[kc][ctx];
				uint16_t *currr2=er2[kc][ctx];
				uint64_t *currhistory=ehistory[kc][ctx];
				
				int den=currhist[256]+256;
				int sym=error<<1^error>>31;
				hsize-=log2((double)(currhist[sym]+1)/den);//estimate doesn't need cdf
				++currhist[sym];
				++currhist[256];
				if(currhist[256]>=0xFFFF-256)
				{
					den=0;
					for(int k=0;k<256;++k)
						den+=currhist[k]>>=1;
					currhist[256]=den;
				}
				int nbypass=(ctx>>2)-GRBITS;
				CLAMP2(nbypass, 0, 7);
				int nzeros=sym>>nbypass;
				if(nzeros>GRLIMIT)
					nzeros=GRLIMIT;
				int den2=currr2[GRLIMIT+1]+GRLIMIT+1;
				//if(ky==ih/2&&kx==iw/2)//
				//	printf("");
				for(int k=0;k<=nzeros&&k<GRLIMIT;++k)
				{
					int recent=(int)(currhistory[k]>>60);
					//int change=currhistory[k]>>60;
					//change=change==0x5||change==0xA;
					//int change=((currhistory[k]>>63)^(currhistory[k]>>62&1))&&((currhistory[k]>>62&1)^(currhistory[k]>>61&1));
					int bit=k>=nzeros&&k<GRLIMIT;
					int p0=currrice[k]>>EPROBSHIFT;
					p0+=p0<0;
					p0+=1<<EPROBBITS_USE>>1;
					//if(recent==15)
					//	p0+=((0<<EPROBBITS_USE)-p0+(1<<4>>1))>>4;
					//else if(recent==0)
					//	p0+=((1<<EPROBBITS_USE)-p0+(1<<4>>1))>>4;
					//if(change)
					//	p0=(1<<EPROBBITS_USE)-p0;
					//p0=squash(p0-2048);
					//double bitsize=-log2((double)(bit?(1<<EPROBBITS_USE)-p0:p0)*(1./(1<<EPROBBITS_USE)));
					//if(isinf(bitsize))
					//	CRASH("");
					rsize-=log2((double)(bit?(1<<EPROBBITS_USE)-p0:p0)*(1./(1<<EPROBBITS_USE)));
					{
						uint16_t *pcell=(uint16_t*)(currctr+k);
						int sum=pcell[0]+pcell[1]+2;
						int p0c=(((pcell[0]+1)<<EPROBBITS_USE)+(sum>>1))/sum;
						//if(change)
						//	p0=(1<<EPROBBITS_USE)-p0;
						csize-=log2((double)(bit?(1<<EPROBBITS_USE)-p0c:p0c)*(1./(1<<EPROBBITS_USE)));
						++pcell[bit];
						if(pcell[bit]>=511)
						{
							pcell[0]>>=1;
							pcell[1]>>=1;
						}
					}
					{
						int p1d=((currr2[k]+1)<<EPROBBITS_USE)/den2;
						//p1d=squash(p1d-2058);
						CLAMP2(p1d, 1, (1<<EPROBBITS_USE)-1);
						r2size-=log2((double)(bit?p1d:(1<<EPROBBITS_USE)-p1d)*(1./(1<<EPROBBITS_USE)));
						//r2size-=log2((double)(1<<EPROBBITS_USE>>1)*(1./(1<<EPROBBITS_USE)));
						//currr2[k]+=!bit;
						//currr2[GRLIMIT]+=!bit;
						//if(currr2[GRLIMIT]>=0xFFFF-GRLIMIT)
						//{
						//	den2=0;
						//	for(int k=0;k<GRLIMIT;++k)
						//		den2+=currr2[k]>>=1;
						//	currr2[GRLIMIT]=den2;
						//}
					}
					{
						currhistory[4*k+0]=(currhistory[4*k+1]&1)<<63|currhistory[4*k+0]>>1;
						currhistory[4*k+1]=(currhistory[4*k+2]&1)<<63|currhistory[4*k+1]>>1;
						currhistory[4*k+2]=(currhistory[4*k+3]&1)<<63|currhistory[4*k+2]>>1;
						currhistory[4*k+3]=(uint64_t)bit<<63|currhistory[4*k+3]>>1;
					}
					currrice[k]+=((!bit<<EPROBBITS_STORE)-(1<<EPROBBITS_STORE>>1)-currrice[k])>>8;
				}
#if 1
				if(ky==ih/2&&kx==iw/2)//
				{
					for(int k=0;k<GRLIMIT;++k)
					{
						for(int kb=64*4-1;kb>=0;--kb)
							printf("%d", (int)(currhistory[4*k+(kb>>6)]>>(kb&63)&1));
						printf("\n");
					}
					printf("\n");
				}
#endif
#if 0
				if(ky==ih/2&&kx==iw/2)//
				{
					//double hist[GRLIMIT+1];
					//for(int k=0;k<GRLIMIT+1;++k)
					//	hist[k]=1;
					//memset(hist, 0, sizeof(hist));
					double pc=1;
					for(int k=0;k<GRLIMIT+1;++k)
					{
						uint16_t *pcell=(uint16_t*)(currctr+k);
						double p0=(double)(pcell[0]+1)/(pcell[0]+pcell[1]+2);
						//hist[k]=p0*pc;
						printf("%3d  %12.6lf  %4d %4d  %12.6lf\n", k, p0*pc, pcell[0], pcell[1], p0);
						pc*=1-p0;
					}
				}
#endif
				//for(int ks=nzeros;ks<=GRLIMIT;++ks)
				//{
				//	++currr2[ks];
				//	++currr2[GRLIMIT+1];
				//	if(currr2[GRLIMIT+1]>=0xFFFF-GRLIMIT)
				//	{
				//		den2=0;
				//		for(int k=0;k<=GRLIMIT;++k)
				//			den2+=currr2[k]>>=1;
				//		currr2[GRLIMIT+1]=den2;
				//	}
				//}
				++currr2[nzeros];
				++currr2[GRLIMIT+1];
				if(currr2[GRLIMIT+1]>=0xFFFF-GRLIMIT)
				{
					den2=0;
					for(int k=0;k<=GRLIMIT;++k)
						den2+=currr2[k]>>=1;
					currr2[GRLIMIT+1]=den2;
				}
				//if(ky==ih/2&&kx==iw/2)//
				//	printf("");
				rsize+=nzeros>=GRLIMIT?8:nbypass;
				for(int kb=7, tidx=1;kb>=0;--kb)
				{
					int bit=error>>kb&1;
					int p0=currstats[tidx]>>EPROBSHIFT;
					p0+=p0<0;
					p0+=1<<EPROBBITS_USE>>1;
					ssize-=log2((double)(bit?(1<<EPROBBITS_USE)-p0:p0)*(1./(1<<EPROBBITS_USE)));
					currstats[tidx]+=((!bit<<EPROBBITS_STORE)-(1<<EPROBBITS_STORE>>1)-currstats[tidx])>>8;
					tidx=2*tidx+bit;
				}
				
				e=curr-pred;
				rows[0][3+0*NCH*NROWS*NVAL]=(2*eW+((e<<1^e>>31)<<GRBITS)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
				curr-=offset;
				e=(curr>pred0)-(curr<pred0);
				bias[kc]+=e<<8;
#define PRED(EXPR) coeffs[kc][j]+=e*estim[j]; ++j;
				j=0;
				PREDLIST;
#undef  PRED
				rows[0][0+0*NCH*NROWS*NVAL]=curr;
				rows[0][1+0*NCH*NROWS*NVAL]=curr-W;
				rows[0][2+0*NCH*NROWS*NVAL]=curr-N;
				offset=(kc ? cv0*yuv[0]+cv1*yuv[1] : cu0*yuv[0])>>2;
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
			}
		}
	}
	free(pixels);
	int64_t usize=(int64_t)3*iw*ih;
	printf("%12.2lf bytes %12.6lf BPD hists\n"	, hsize/8, hsize/usize);
	printf("%12.2lf bytes %12.6lf BPD stats\n"	, ssize/8, ssize/usize);
	printf("%12.2lf bytes %12.6lf BPD Rice\n"	, rsize/8, rsize/usize);
	printf("%12.2lf bytes %12.6lf BPD Rice ctr\n"	, csize/8, csize/usize);
	printf("%12.2lf bytes %12.6lf BPD Rice ctr2\n"	, r2size/8, r2size/usize);
}
#endif
#ifdef PRINT_ERRORCOUNTS
int64_t g_n[2];
enum
{
	ERRORHISTBITS=7,
};
int32_t errorhist[1<<ERRORHISTBITS]={0};
#endif
INLINE void codebit(ACState *ac, uint64_t *pcell, int32_t *bit, const int fwd, int unary)
{
	uint64_t r2, mid;
	int rbit, sh;
	uint64_t entry=ctrtable[*pcell];
	int32_t p1=(int32_t)(entry&((1ULL<<PROBBITS_USE)-1));
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
	r2=ac->range*(uint32_t)p1>>PROBBITS_USE;
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
	ac->range=rbit?r2:ac->range;
	ac->low=rbit?ac->low:mid;
	//uint64_t mask=-(uint64_t)rbit;
	//ac->range^=(r2^ac->range)&mask;
	//sh=(uint8_t)(PROBBITS_USE+(CTRBITS*2LL&mask));
	//ac->low^=(mid^ac->low)&~mask;
	sh=rbit?PROBBITS_USE+CTRBITS*2:PROBBITS_USE;
	*pcell=entry>>sh&((1ULL<<CTRBITS*2)-1);
#ifdef ESTIMATE_BITSIZE
	bitsizes[ekc][eidx]+=shannontable[rbit?p1:(1<<PROBBITS_USE)-p1];
	++bitctr[ekc][eidx][rbit];
	winctr[ekc][eidx]+=rbit==(p1>=1<<PROBBITS_USE);
#endif
#ifdef PRINT_ERRORCOUNTS
	if(unary)
	{
		static int ctr=0;
		if(fwd&&(uint32_t)(ctr-10000000)<100000)
		{
			int errorgrade=abs((rbit<<PROBBITS_USE)-p1)>>(PROBBITS_USE-ERRORHISTBITS);
			++errorhist[errorgrade];
			printf("%c", "0123456789ABCDEF"[errorgrade>>(ERRORHISTBITS-4)]);
			//char c=rbit?'+':'-';
			//if(abs(p1-(rbit<<PROBBITS_USE))>1<<PROBBITS_USE>>1)
			//	c=rbit?'1':'0';
			//printf("%c", c);
		}
		if(fwd&&abs(p1-(rbit<<PROBBITS_USE))>1<<PROBBITS_USE>>1)
			++g_n[rbit];
		++ctr;
	}
#endif
}
INLINE void mainloop(int iw, int ih, RCTInfo *rct, int dist, uint8_t *image, uint8_t *stream, ACState *ac, const int lossy, const int fwd)
{
	int
		yidx=rct->perm[0],
		uidx=rct->perm[1],
		vidx=rct->perm[2],
		cu0=rct->cu0,
		cv0=rct->cv0,
		cv1=rct->cv1;
	int32_t ky, kx;
	int32_t psize=0;
	int16_t *pixels=0;
	ALIGN(32) int32_t lcoeffs[3][L1NPREDS_LOSSY+1]={0}, coeffs[3][L1NPREDS]={0}, bias[3]={0};
	int32_t invdist=((1<<16)+dist-1)/dist;
	uint8_t *imptr=image;
#ifdef PRINTBITS
	ptrdiff_t idx=0, usize=(ptrdiff_t)3*iw*ih;
#endif
#ifdef PRINT_STATS
	int64_t unarysum=0, unarycount=0;
	uint32_t statctr[GRLIMIT+1][2]={0};
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
#ifdef USE_COUNTERS
	for(int k=0;k<1<<CTRBITS*2;++k)//MSB {next_1(n0, n1), next_0(n0, n1}, p1} LSB
	{
		int32_t n[]={k&CTRMASK, k>>CTRBITS&CTRMASK};
		int n0e=n[0]+2;
		int n1e=n[1]+2;
		int sum=n0e+n1e;
		int32_t p1=((n1e<<PROBBITS_USE)+(sum>>1))/sum;
		uint64_t cell=0;
		int32_t n0[]={k&CTRMASK, k>>CTRBITS&CTRMASK};
		int32_t n1[]={k&CTRMASK, k>>CTRBITS&CTRMASK};
		CLAMP2(p1, 1, ((1<<PROBBITS_USE)-1));
#if 1
		if(n0[0]>=CTRMASK)
		{
			n0[0]>>=1;
			n0[1]>>=1;
		}
		++n0[0];
		if(n1[1]>=CTRMASK)
		{
			n1[0]>>=1;
			n1[1]>>=1;
		}
		++n1[1];
#endif
#if 0
		++n0[0];
		if(n0[0]>CTRMASK)
		{
			n0[0]>>=1;
			n0[1]>>=1;
		}
		++n1[1];
		if(n1[1]>CTRMASK)
		{
			n1[0]>>=1;
			n1[1]>>=1;
		}
#endif
		cell=cell<<CTRBITS|n1[1];
		cell=cell<<CTRBITS|n1[0];
		cell=cell<<CTRBITS|n0[1];
		cell=cell<<CTRBITS|n0[0];
		cell=cell<<PROBBITS_USE|p1;
		ctrtable[k]=cell;
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
					eN	=rows[1][2+0*NCH*NROWS*NVAL],
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
				int32_t upred2, pred0;
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
					//	pred>>=L1SH;
					}
				}
				pred0=upred2>>2;
#ifdef RCT_BASE8
				upred2+=offset0>>1;
#elif defined RCT_BASE64
				upred2+=offset0>>4;
#else
#error
#endif
				CLAMP2(upred2, 0, 1023);
#if 0
				if(!kc)
				{
					int vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(pred, vmin, vmax);
				}
				else
				{
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
				}
#endif
				pred=upred2>>2;
#ifdef USE_TABLES
				uint8_t *ctxptr=ctxtable[eW<_countof(ctxtable)-1?eW:_countof(ctxtable)-1];
				ctx=ctxptr[0];
				nbypass=ctxptr[1];
				//ctx=abs(N-W)>>3;
				//if(ctx>NCTX-1)
				//	ctx=NCTX-1;
#else
				nbypass=ctx=31-_lzcnt_u32(eW*eW+2);
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
			//	(void)eNE	;
			//	(void)eNEE	;
			//	(void)eNEEE	;
			//	(void)eW	;
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
							uint8_t *pval=&g_image[3*(iw*ky+kx)+rct->perm[kc]];
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
						int e2=(int8_t)error;
						int negmask=error>>31;
						int abserr=(error^negmask)-negmask;
						error=error<<1^negmask;
						e2=e2<<1^e2>>31;
						if(epred<abserr)
							error=epred+abserr;
						if(error==256)
							error=e2;
#endif
					}
					nzeros=error>>nbypass;
#ifdef PRINT_STATS
					unarysum+=nzeros;
					++unarycount;
#endif
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
					++unary_count;
#endif
					codebit(ac, statsptr+tidx, &bit, fwd, 1);
					//bits[tidx]=bit;
#ifdef PRINT_STATS
					++statctr[tidx][bit];
#endif
#ifdef PRINTBITS2
					//if(fwd&&!tidx)
					//{
					//	static int fctr=0, bitctr=0, reg=0;
					//	static FILE *f=0;
					//	++fctr;
					//	++bitctr;
					//	reg=reg<<1|bit;
					//	if(fctr==1)
					//	{
					//		f=fopen("20260404_0531pm.dat", "wb");
					//		if(!f)
					//		{
					//			CRASH("Cannot write");
					//			return;
					//		}
					//	}
					//	if(bitctr==32)
					//	{
					//		fwrite(&reg, 1, 4, f);
					//		reg=0;
					//		bitctr=0;
					//	}
					//}
					//if(fwd&&!tidx)
					//{
					//	static int nbits=0, val=0;
					//	val|=bit<<nbits++;
					//	if(nbits>=32)
					//	{
					//		printf("%08X", val);
					//		val=0;
					//		nbits=0;
					//	}
					//	//printf("%c", '0'+bit);
					//}
#endif
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
						++binary_count;
#endif
						codebit(ac, statsptr+tidx, &bit, fwd, 0);
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
						uint8_t *pval=&g_image[3*(iw*ky+kx)+rct->perm[kc]];
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
				{
					int32_t curr=yuv[kc]-offset;//FIXME prec?
					int32_t k, e;

					e=(curr>pred0)-(curr<pred0);

					if(lossy)
					{
						//currw[L1NPREDS_LOSSY]+=e;
#define PRED(EXPR) lcoeffs[kc][k]+=e*preds[k]; ++k;
						k=0;
						PREDLIST_LOSSY;
#undef  PRED
					}
					else
					{
						bias[kc]+=e<<8;
#define PRED(EXPR) coeffs[kc][k]+=e*preds[k]; ++k;
						k=0;
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
					
				//	rows[0][3]=(16*eW+7*(error<<GRBITS)+9*(eNEE>eNEEE?eNEE:eNEEE))>>5;
				//	rows[0][3]=(2*eW+(error<<GRBITS)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
					rows[0][3]=(eW+(eW<eNE?eW:eNE)+(error<<GRBITS)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
					offset0=kc ? cv0*yuv[0]+cv1*yuv[1] : cu0*yuv[0];
#ifdef RCT_BASE8
					offset=offset0>>3;
#elif defined RCT_BASE64
					offset=offset0>>6;
#else
#error
#endif
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
#ifdef PRINT_ERRORCOUNTS
	if(fwd)
	{
		printf("\n%lld vs %lld\n", g_n[0], g_n[1]);
		int vmax=0;
		for(int k=0;k<1<<ERRORHISTBITS;++k)
		{
			if(vmax<errorhist[k])
				vmax=errorhist[k];
		}
		for(int k=0;k<1<<ERRORHISTBITS;++k)
		{
			printf("%4d  %9d  ", k, errorhist[k]);
			for(int k2=0, nstars=errorhist[k]*64/vmax;k2<nstars;++k2)
				printf("*");
			printf("\n");
		}
	}
#endif
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
	RCTInfo rctinfo={0};
//	int bestrct=0;
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
	//int32_t gaintable[3*(GRLIMIT+1)]={0};

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
			iw=0;
			ih=0;
		//	bestrct=0;
			dist=1;
			fread(&iw, 1, 3, fsrc);
			fread(&ih, 1, 3, fsrc);
		//	fread(&bestrct, 1, 1, fsrc);
			fread(&rctinfo.perm[0], 1, 1, fsrc);
			fread(&rctinfo.perm[1], 1, 1, fsrc);
			fread(&rctinfo.perm[2], 1, 1, fsrc);
			fread(&rctinfo.cu0, 1, 1, fsrc);
			fread(&rctinfo.cv0, 1, 1, fsrc);
			fread(&rctinfo.cv1, 1, 1, fsrc);
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
#ifdef ANALYSIS_BY_ENTROPY
		int hsize=(int)sizeof(int32_t[OCH_C1X0*256]);
		int32_t *hists=(int32_t*)malloc(hsize);
		if(!hists)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(hists, 0, hsize);
#else
		int64_t counters[OCH_COUNT]={0};
#endif
#ifdef ANALYSIS_NB
		enum
		{
			AXPAD=8,
			ACH=OCH_COUNT,
			AROWS=4,
			AVAL=1,
		};
		int psize=(iw+2*XPAD)*(int)sizeof(int16_t[ACH*AROWS*AVAL]);
		int16_t *pixels=(int16_t*)malloc(psize);
#ifdef ANALYSIS_L1
		int wsize=(int)sizeof(int32_t[ACH*L1NPREDS]);
		int32_t *weights=(int32_t*)malloc(wsize);
		if(!weights)
		{
			CRASH("Alloc error");
			return 1;
		}
		for(int k=0;k<ACH*L1NPREDS;++k)
			weights[k]=(1<<L1SH)/L1NPREDS;
#endif

		if(!pixels)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(pixels, 0, psize);
		imptr=image;
		for(int ky=0;ky<ih;++ky)
		{
			int16_t *rows[]=
			{
				pixels+(AXPAD*ACH*AROWS+(ky-0LL+AROWS)%AROWS)*AVAL,
				pixels+(AXPAD*ACH*AROWS+(ky-1LL+AROWS)%AROWS)*AVAL,
				pixels+(AXPAD*ACH*AROWS+(ky-2LL+AROWS)%AROWS)*AVAL,
				pixels+(AXPAD*ACH*AROWS+(ky-3LL+AROWS)%AROWS)*AVAL,
			};
			ALIGN(16) int16_t rramp[8]={0}, gramp[8]={0}, bramp[8]={0};
			for(int kx=0;kx<iw;++kx)
			{
				int r, g, b;

				r=imptr[0]<<3;
				g=imptr[1]<<3;
				b=imptr[2]<<3;
				imptr+=3;

				rramp[1-1]=r*1>>3;
				gramp[1-1]=g*1>>3;
				bramp[1-1]=b*1>>3;
				rramp[2-1]=r*2>>3;
				gramp[2-1]=g*2>>3;
				bramp[2-1]=b*2>>3;
				rramp[3-1]=r*3>>3;
				gramp[3-1]=g*3>>3;
				bramp[3-1]=b*3>>3;
				rramp[4-1]=r*4>>3;
				gramp[4-1]=g*4>>3;
				bramp[4-1]=b*4>>3;
				rramp[5-1]=r*5>>3;
				gramp[5-1]=g*5>>3;
				bramp[5-1]=b*5>>3;
				rramp[6-1]=r*6>>3;
				gramp[6-1]=g*6>>3;
				bramp[6-1]=b*6>>3;
				rramp[7-1]=r*7>>3;
				gramp[7-1]=g*7>>3;
				bramp[7-1]=b*7>>3;
				rramp[8-1]=r*8>>3;
				gramp[8-1]=g*8>>3;
				bramp[8-1]=b*8>>3;
				int och[OCH_COUNT]=
				{
					r,
					g,
					b,

					r-gramp[1-1],
					g-bramp[1-1],
					b-rramp[1-1],
					r-gramp[2-1],
					g-bramp[2-1],
					b-rramp[2-1],
					r-gramp[3-1],
					g-bramp[3-1],
					b-rramp[3-1],
					r-gramp[4-1],
					g-bramp[4-1],
					b-rramp[4-1],
					r-gramp[5-1],
					g-bramp[5-1],
					b-rramp[5-1],
					r-gramp[6-1],
					g-bramp[6-1],
					b-rramp[6-1],
					r-gramp[7-1],
					g-bramp[7-1],
					b-rramp[7-1],
					r-gramp[8-1],
					g-bramp[8-1],
					b-rramp[8-1],

					r-bramp[1-1],
					g-rramp[1-1],
					b-gramp[1-1],
					r-bramp[2-1],
					g-rramp[2-1],
					b-gramp[2-1],
					r-bramp[3-1],
					g-rramp[3-1],
					b-gramp[3-1],
					r-bramp[4-1],
					g-rramp[4-1],
					b-gramp[4-1],
					r-bramp[5-1],
					g-rramp[5-1],
					b-gramp[5-1],
					r-bramp[6-1],
					g-rramp[6-1],
					b-gramp[6-1],
					r-bramp[7-1],
					g-rramp[7-1],
					b-gramp[7-1],
					r-bramp[8-1],
					g-rramp[8-1],
					b-gramp[8-1],

					r-(gramp[1-1]+bramp[1-1]),
					g-(bramp[1-1]+rramp[1-1]),
					b-(rramp[1-1]+gramp[1-1]),

					r-(gramp[1-1]+bramp[2-1]),
					g-(bramp[1-1]+rramp[2-1]),
					b-(rramp[1-1]+gramp[2-1]),
					r-(gramp[2-1]+bramp[1-1]),
					g-(bramp[2-1]+rramp[1-1]),
					b-(rramp[2-1]+gramp[1-1]),

					r-(gramp[1-1]+bramp[3-1]),
					g-(bramp[1-1]+rramp[3-1]),
					b-(rramp[1-1]+gramp[3-1]),
					r-(gramp[2-1]+bramp[2-1]),
					g-(bramp[2-1]+rramp[2-1]),
					b-(rramp[2-1]+gramp[2-1]),
					r-(gramp[3-1]+bramp[1-1]),
					g-(bramp[3-1]+rramp[1-1]),
					b-(rramp[3-1]+gramp[1-1]),

					r-(gramp[1-1]+bramp[4-1]),
					g-(bramp[1-1]+rramp[4-1]),
					b-(rramp[1-1]+gramp[4-1]),
					r-(gramp[2-1]+bramp[3-1]),
					g-(bramp[2-1]+rramp[3-1]),
					b-(rramp[2-1]+gramp[3-1]),
					r-(gramp[3-1]+bramp[2-1]),
					g-(bramp[3-1]+rramp[2-1]),
					b-(rramp[3-1]+gramp[2-1]),
					r-(gramp[4-1]+bramp[1-1]),
					g-(bramp[4-1]+rramp[1-1]),
					b-(rramp[4-1]+gramp[1-1]),

					r-(gramp[1-1]+bramp[5-1]),
					g-(bramp[1-1]+rramp[5-1]),
					b-(rramp[1-1]+gramp[5-1]),
					r-(gramp[2-1]+bramp[4-1]),
					g-(bramp[2-1]+rramp[4-1]),
					b-(rramp[2-1]+gramp[4-1]),
					r-(gramp[3-1]+bramp[3-1]),
					g-(bramp[3-1]+rramp[3-1]),
					b-(rramp[3-1]+gramp[3-1]),
					r-(gramp[4-1]+bramp[2-1]),
					g-(bramp[4-1]+rramp[2-1]),
					b-(rramp[4-1]+gramp[2-1]),
					r-(gramp[5-1]+bramp[1-1]),
					g-(bramp[5-1]+rramp[1-1]),
					b-(rramp[5-1]+gramp[1-1]),

					r-(gramp[1-1]+bramp[6-1]),
					g-(bramp[1-1]+rramp[6-1]),
					b-(rramp[1-1]+gramp[6-1]),
					r-(gramp[2-1]+bramp[5-1]),
					g-(bramp[2-1]+rramp[5-1]),
					b-(rramp[2-1]+gramp[5-1]),
					r-(gramp[3-1]+bramp[4-1]),
					g-(bramp[3-1]+rramp[4-1]),
					b-(rramp[3-1]+gramp[4-1]),
					r-(gramp[4-1]+bramp[3-1]),
					g-(bramp[4-1]+rramp[3-1]),
					b-(rramp[4-1]+gramp[3-1]),
					r-(gramp[5-1]+bramp[2-1]),
					g-(bramp[5-1]+rramp[2-1]),
					b-(rramp[5-1]+gramp[2-1]),
					r-(gramp[6-1]+bramp[1-1]),
					g-(bramp[6-1]+rramp[1-1]),
					b-(rramp[6-1]+gramp[1-1]),

					r-(gramp[1-1]+bramp[7-1]),
					g-(bramp[1-1]+rramp[7-1]),
					b-(rramp[1-1]+gramp[7-1]),
					r-(gramp[2-1]+bramp[6-1]),
					g-(bramp[2-1]+rramp[6-1]),
					b-(rramp[2-1]+gramp[6-1]),
					r-(gramp[3-1]+bramp[5-1]),
					g-(bramp[3-1]+rramp[5-1]),
					b-(rramp[3-1]+gramp[5-1]),
					r-(gramp[4-1]+bramp[4-1]),
					g-(bramp[4-1]+rramp[4-1]),
					b-(rramp[4-1]+gramp[4-1]),
					r-(gramp[5-1]+bramp[3-1]),
					g-(bramp[5-1]+rramp[3-1]),
					b-(rramp[5-1]+gramp[3-1]),
					r-(gramp[6-1]+bramp[2-1]),
					g-(bramp[6-1]+rramp[2-1]),
					b-(rramp[6-1]+gramp[2-1]),
					r-(gramp[7-1]+bramp[1-1]),
					g-(bramp[7-1]+rramp[1-1]),
					b-(rramp[7-1]+gramp[1-1]),
				};
				for(int kc=0;kc<OCH_COUNT;++kc)
				{
					int
						NNNN	=rows[0][(kc+0*ACH)*AROWS*AVAL],
						NNN	=rows[3][(kc+0*ACH)*AROWS*AVAL],
						NNW	=rows[2][(kc-1*ACH)*AROWS*AVAL],
						NN	=rows[2][(kc+0*ACH)*AROWS*AVAL],
						NNE	=rows[2][(kc+1*ACH)*AROWS*AVAL],
						NWW	=rows[1][(kc-2*ACH)*AROWS*AVAL],
						NW	=rows[1][(kc-1*ACH)*AROWS*AVAL],
						N	=rows[1][(kc+0*ACH)*AROWS*AVAL],
						NE	=rows[1][(kc+1*ACH)*AROWS*AVAL],
						NEE	=rows[1][(kc+2*ACH)*AROWS*AVAL],
						NEEEE	=rows[1][(kc+4*ACH)*AROWS*AVAL],
						WWWW	=rows[0][(kc-4*ACH)*AROWS*AVAL],
						WWW	=rows[0][(kc-3*ACH)*AROWS*AVAL],
						WW	=rows[0][(kc-2*ACH)*AROWS*AVAL],
						W	=rows[0][(kc-1*ACH)*AROWS*AVAL];
#if 1
					(void)NNNN	;
					(void)NNN	;
					(void)NNW	;
					(void)NN	;
					(void)NNE	;
					(void)NWW	;
					(void)NW	;
					(void)N		;
					(void)NE	;
					(void)NEE	;
					(void)NEEEE	;
					(void)WWWW	;
					(void)WWW	;
					(void)WW	;
					(void)W		;
#endif
#ifdef ANALYSIS_L1
					int estim[]=
					{
#define PRED(E) E,
						PREDLIST
#undef  PRED
					};
					int32_t *coeffs=weights+L1NPREDS*kc;
					int64_t pred=1LL<<L1SH>>1;
					for(int k=0;k<L1NPREDS;++k)
						pred+=coeffs[k]*estim[k];
					pred>>=L1SH;
					int e=(och[kc]>pred)-(och[kc]<pred);
					for(int k=0;k<L1NPREDS;++k)
						coeffs[k]+=e*estim[k];
#else
				//	int pred=N+W-NW;
				//	int vmax=N, vmin=W;
				//	if(N<W)vmin=N, vmax=W;
				//	CLAMP2(pred, vmin, vmax);

					int pred=(N+W)>>1;

				//	int pred=N+W-NW;

				//	int pred=W;
#endif
#ifdef ANALYSIS_BY_ENTROPY
					++hists[kc*256+((((och[kc]-pred)>>3)+128)&255)];
#else
					counters[kc]+=abs(och[kc]-(int32_t)pred);
#endif
					rows[0][(kc+0*ACH)*AROWS*AVAL]=och[kc];
				}
				rows[0]+=ACH*AROWS*AVAL;
				rows[1]+=ACH*AROWS*AVAL;
			}
		}
#ifdef ANALYSIS_L1
		free(weights);
#endif
		free(pixels);
#else
		ALIGN(32) int prev[OCH_COUNT]={0};
		int rowstride=3*iw;

		imptr=image+rowstride;
		{
			//int res=iw*ih, idx=0;//

			ALIGN(16) int16_t rramp[8]={0}, gramp[8]={0}, bramp[8]={0};
			while(imptr<imend)
			{
				int r, g, b;

				r=(imptr[0]-imptr[0-rowstride])<<3;
				g=(imptr[1]-imptr[1-rowstride])<<3;
				b=(imptr[2]-imptr[2-rowstride])<<3;
				imptr+=3;

				rramp[1-1]=r*1>>3;
				gramp[1-1]=g*1>>3;
				bramp[1-1]=b*1>>3;
				rramp[2-1]=r*2>>3;
				gramp[2-1]=g*2>>3;
				bramp[2-1]=b*2>>3;
				rramp[3-1]=r*3>>3;
				gramp[3-1]=g*3>>3;
				bramp[3-1]=b*3>>3;
				rramp[4-1]=r*4>>3;
				gramp[4-1]=g*4>>3;
				bramp[4-1]=b*4>>3;
				rramp[5-1]=r*5>>3;
				gramp[5-1]=g*5>>3;
				bramp[5-1]=b*5>>3;
				rramp[6-1]=r*6>>3;
				gramp[6-1]=g*6>>3;
				bramp[6-1]=b*6>>3;
				rramp[7-1]=r*7>>3;
				gramp[7-1]=g*7>>3;
				bramp[7-1]=b*7>>3;
				rramp[8-1]=r*8>>3;
				gramp[8-1]=g*8>>3;
				bramp[8-1]=b*8>>3;

#if 0
				if(idx>=res/4)//
					printf("");
				++idx;
#endif
				int och[OCH_COUNT]=
				{
					r,
					g,
					b,

					r-gramp[1-1],
					g-bramp[1-1],
					b-rramp[1-1],
					r-gramp[2-1],
					g-bramp[2-1],
					b-rramp[2-1],
					r-gramp[3-1],
					g-bramp[3-1],
					b-rramp[3-1],
					r-gramp[4-1],
					g-bramp[4-1],
					b-rramp[4-1],
					r-gramp[5-1],
					g-bramp[5-1],
					b-rramp[5-1],
					r-gramp[6-1],
					g-bramp[6-1],
					b-rramp[6-1],
					r-gramp[7-1],
					g-bramp[7-1],
					b-rramp[7-1],
					r-gramp[8-1],
					g-bramp[8-1],
					b-rramp[8-1],

					r-bramp[1-1],
					g-rramp[1-1],
					b-gramp[1-1],
					r-bramp[2-1],
					g-rramp[2-1],
					b-gramp[2-1],
					r-bramp[3-1],
					g-rramp[3-1],
					b-gramp[3-1],
					r-bramp[4-1],
					g-rramp[4-1],
					b-gramp[4-1],
					r-bramp[5-1],
					g-rramp[5-1],
					b-gramp[5-1],
					r-bramp[6-1],
					g-rramp[6-1],
					b-gramp[6-1],
					r-bramp[7-1],
					g-rramp[7-1],
					b-gramp[7-1],
					r-bramp[8-1],
					g-rramp[8-1],
					b-gramp[8-1],

					r-(gramp[1-1]+bramp[1-1]),
					g-(bramp[1-1]+rramp[1-1]),
					b-(rramp[1-1]+gramp[1-1]),

					r-(gramp[1-1]+bramp[2-1]),
					g-(bramp[1-1]+rramp[2-1]),
					b-(rramp[1-1]+gramp[2-1]),
					r-(gramp[2-1]+bramp[1-1]),
					g-(bramp[2-1]+rramp[1-1]),
					b-(rramp[2-1]+gramp[1-1]),

					r-(gramp[1-1]+bramp[3-1]),
					g-(bramp[1-1]+rramp[3-1]),
					b-(rramp[1-1]+gramp[3-1]),
					r-(gramp[2-1]+bramp[2-1]),
					g-(bramp[2-1]+rramp[2-1]),
					b-(rramp[2-1]+gramp[2-1]),
					r-(gramp[3-1]+bramp[1-1]),
					g-(bramp[3-1]+rramp[1-1]),
					b-(rramp[3-1]+gramp[1-1]),

					r-(gramp[1-1]+bramp[4-1]),
					g-(bramp[1-1]+rramp[4-1]),
					b-(rramp[1-1]+gramp[4-1]),
					r-(gramp[2-1]+bramp[3-1]),
					g-(bramp[2-1]+rramp[3-1]),
					b-(rramp[2-1]+gramp[3-1]),
					r-(gramp[3-1]+bramp[2-1]),
					g-(bramp[3-1]+rramp[2-1]),
					b-(rramp[3-1]+gramp[2-1]),
					r-(gramp[4-1]+bramp[1-1]),
					g-(bramp[4-1]+rramp[1-1]),
					b-(rramp[4-1]+gramp[1-1]),

					r-(gramp[1-1]+bramp[5-1]),
					g-(bramp[1-1]+rramp[5-1]),
					b-(rramp[1-1]+gramp[5-1]),
					r-(gramp[2-1]+bramp[4-1]),
					g-(bramp[2-1]+rramp[4-1]),
					b-(rramp[2-1]+gramp[4-1]),
					r-(gramp[3-1]+bramp[3-1]),
					g-(bramp[3-1]+rramp[3-1]),
					b-(rramp[3-1]+gramp[3-1]),
					r-(gramp[4-1]+bramp[2-1]),
					g-(bramp[4-1]+rramp[2-1]),
					b-(rramp[4-1]+gramp[2-1]),
					r-(gramp[5-1]+bramp[1-1]),
					g-(bramp[5-1]+rramp[1-1]),
					b-(rramp[5-1]+gramp[1-1]),

					r-(gramp[1-1]+bramp[6-1]),
					g-(bramp[1-1]+rramp[6-1]),
					b-(rramp[1-1]+gramp[6-1]),
					r-(gramp[2-1]+bramp[5-1]),
					g-(bramp[2-1]+rramp[5-1]),
					b-(rramp[2-1]+gramp[5-1]),
					r-(gramp[3-1]+bramp[4-1]),
					g-(bramp[3-1]+rramp[4-1]),
					b-(rramp[3-1]+gramp[4-1]),
					r-(gramp[4-1]+bramp[3-1]),
					g-(bramp[4-1]+rramp[3-1]),
					b-(rramp[4-1]+gramp[3-1]),
					r-(gramp[5-1]+bramp[2-1]),
					g-(bramp[5-1]+rramp[2-1]),
					b-(rramp[5-1]+gramp[2-1]),
					r-(gramp[6-1]+bramp[1-1]),
					g-(bramp[6-1]+rramp[1-1]),
					b-(rramp[6-1]+gramp[1-1]),

					r-(gramp[1-1]+bramp[7-1]),
					g-(bramp[1-1]+rramp[7-1]),
					b-(rramp[1-1]+gramp[7-1]),
					r-(gramp[2-1]+bramp[6-1]),
					g-(bramp[2-1]+rramp[6-1]),
					b-(rramp[2-1]+gramp[6-1]),
					r-(gramp[3-1]+bramp[5-1]),
					g-(bramp[3-1]+rramp[5-1]),
					b-(rramp[3-1]+gramp[5-1]),
					r-(gramp[4-1]+bramp[4-1]),
					g-(bramp[4-1]+rramp[4-1]),
					b-(rramp[4-1]+gramp[4-1]),
					r-(gramp[5-1]+bramp[3-1]),
					g-(bramp[5-1]+rramp[3-1]),
					b-(rramp[5-1]+gramp[3-1]),
					r-(gramp[6-1]+bramp[2-1]),
					g-(bramp[6-1]+rramp[2-1]),
					b-(rramp[6-1]+gramp[2-1]),
					r-(gramp[7-1]+bramp[1-1]),
					g-(bramp[7-1]+rramp[1-1]),
					b-(rramp[7-1]+gramp[1-1]),
				};
#ifdef ANALYSIS_BY_ENTROPY
				for(int k=0;k<OCH_COUNT;++k)
					++hists[k*256+((((och[k]-prev[k])>>3)+128)&255)];
#else
				for(int k=0;k<OCH_COUNT;++k)
					counters[k]+=abs(och[k]-prev[k]);
#endif
				for(int k=0;k<OCH_COUNT;++k)
					prev[k]=och[k];
			}
		}
#endif
#ifdef ANALYSIS_BY_ENTROPY
		double csizes[OCH_COUNT]={0};
		for(int kc=0;kc<OCH_COUNT;++kc)
		{
			int32_t *hist=hists+kc*256;
			int32_t sum=0;
			for(int ks=0;ks<256;++ks)
				sum+=hist[ks];
			double norm=1./sum, e=0;
			for(int ks=0;ks<256;++ks)
			{
				int freq=hist[ks];
				if(freq)
					e-=freq*log2(freq*norm);
			}
			csizes[kc]=e/8;
		}
		//memcpy(csizes+OCH_C1X0, csizes+OCH_CX10, sizeof(double[OCH_COUNT-OCH_C1X0]));
#ifdef PRINT_RCT
#if 0
		//memcpy(counters+OCH_C1X0, counters+OCH_CX10, sizeof(int64_t[OCH_COUNT-OCH_C1X0]));
		for(int kc=0;kc<OCH_COUNT-2;kc+=3)
			printf("%4d  %s  %s  %s  %12.2lf  %12.2lf  %12.2lf\n"
				, kc
				, och_names[kc+0]
				, och_names[kc+1]
				, och_names[kc+2]
				, counters[kc+0]/csizes[kc+0]
				, counters[kc+1]/csizes[kc+1]
				, counters[kc+2]/csizes[kc+2]
			);
		printf("\n");
#endif

		for(int kc=0;kc<OCH_COUNT-2;kc+=3)
			printf("%4d  %s  %s  %s  %12.2lf  %12.2lf  %12.2lf\n"
				, kc
				, och_names[kc+0]
				, och_names[kc+1]
				, och_names[kc+2]
				, csizes[kc+0]
				, csizes[kc+1]
				, csizes[kc+2]
			);
#endif
#else
		//memcpy(counters+OCH_C1X0, counters+OCH_CX10, sizeof(int64_t[OCH_COUNT-OCH_C1X0]));
#ifdef PRINT_RCT
		for(int kc=0;kc<OCH_COUNT-2;kc+=3)
			printf("%4d  %s  %s  %s  %12lld  %12lld  %12lld\n"
				, kc
				, och_names[kc+0]
				, och_names[kc+1]
				, och_names[kc+2]
				, counters[kc+0]
				, counters[kc+1]
				, counters[kc+2]
			);
#endif
#endif
		{
			static const int perms[]=
			{
				0, 1, 2,	0,
				1, 2, 0,	0,
				2, 0, 1,	0,
				1, 0, 2,	1,
				0, 2, 1,	1,
				2, 1, 0,	1,
			};
#ifdef ANALYSIS_BY_ENTROPY
			double bestval=0;
#else
			int64_t bestval=0;
#endif
			int it_count=0, it_select=0;
			int best0=0, best1=0, best2=0;
#ifdef PRINT_RCT
			RCTInfo rctinfo2={0};
#endif

			for(int kp=0;kp<_countof(perms)/4;++kp)
			{
				int yidx=perms[4*kp+0];
				int uidx=perms[4*kp+1];
				int vidx=perms[4*kp+2];
				int c0=yidx;
				for(int k1=0;k1<OCH_CX11/3;++k1)
				{
					int c1=k1*3+uidx;
					const char *label=och_names[c1];
					if(label[vidx+1]!='0')
						continue;
					for(int k2=0;k2<OCH_COUNT/3;++k2)
					{
						int c2=k2*3+vidx;
#ifdef ANALYSIS_BY_ENTROPY
						double val=
							+csizes[c0]
							+csizes[c1]
							+csizes[c2]
						;
#else
						int64_t val=
							+counters[c0]
							+counters[c1]
							+counters[c2]
						;
#endif

						//if(it_count==867)//
						//if(it_count==37)//
						//if(it_count==562)//
						//if(it_count==45)//
						//if(it_count==360)//
						//if(it_count==2428)//
						//	printf("");

						if(!it_count||bestval>val)
						{
							//int mirror=perms[4*kp+3];

							//if(mirror!=((uidx+1)%3==yidx)||mirror!=((vidx+1)%3!=yidx))
							//	CRASH("");

							bestval=val;
							best0=c0;
							best1=c1;
							best2=c2;
							//rctinfo.perm[0]=yidx;
							//rctinfo.perm[1]=uidx;
							//rctinfo.perm[2]=vidx;
							//rctinfo.cu0=(uint8_t)cRCT_out_coeffs1[k1];
							//rctinfo.cv0=(uint8_t)cRCT_out_coeffs1[k2];
							//rctinfo.cv1=(uint8_t)cRCT_out_coeffs2[k2];
							//rctinfo.cu0=(uint8_t)((uidx+1)%3!=yidx?cRCT_out_coeffs1:cRCT_out_coeffs2)[k1];
							//if((vidx+1)%3!=yidx)
							//{
							//	rctinfo.cv0=(uint8_t)cRCT_out_coeffs2[k2];
							//	rctinfo.cv1=(uint8_t)cRCT_out_coeffs1[k2];
							//}
							//else
							//{
							//	rctinfo.cv0=(uint8_t)cRCT_out_coeffs1[k2];
							//	rctinfo.cv1=(uint8_t)cRCT_out_coeffs2[k2];
							//}
							it_select=it_count;
						}
#ifdef PRINT_RCT
#if 1
						crct_get(&rctinfo2, c0, c1, c2);
						//rctinfo2.perm[0]=yidx;
						//rctinfo2.perm[1]=uidx;
						//rctinfo2.perm[2]=vidx;
						//rctinfo2.cu0=(uint8_t)cRCT_out_coeffs1[k1];
						//rctinfo2.cv0=(uint8_t)cRCT_out_coeffs1[k2];
						//rctinfo2.cv1=(uint8_t)cRCT_out_coeffs2[k2];
						//rctinfo2.cu0=(uint8_t)((uidx+1)%3!=yidx?cRCT_out_coeffs1:cRCT_out_coeffs2)[k1];
						//if((vidx+1)%3==yidx)
						//{
						//	rctinfo2.cv0=(uint8_t)cRCT_out_coeffs1[k2];
						//	rctinfo2.cv1=(uint8_t)cRCT_out_coeffs2[k2];
						//}
						//else
						//{
						//	rctinfo2.cv0=(uint8_t)cRCT_out_coeffs2[k2];
						//	rctinfo2.cv1=(uint8_t)cRCT_out_coeffs1[k2];
						//}
						print_rct2(&rctinfo2);
						printf(" %d %3d %3d _%s_%s_%s  "
							, kp, k1, k2
							, och_names[c0]
							, och_names[c1]
							, och_names[c2]
						);
#endif
#ifdef ANALYSIS_BY_ENTROPY
						printf("%4d  %12.2lf +%12.2lf +%12.2lf =%12.2lf%s\n"
							, it_count
							, csizes[c0]
							, csizes[c1]
							, csizes[c2]
							, val
							, it_count==it_select?" <-":""
						);
#else
						printf("%4d  %12lld +%12lld +%12lld =%12lld%s\n"
							, it_count
							, counters[c0]
							, counters[c1]
							, counters[c2]
							, val
							, it_count==it_select?" <-":""
						);
#endif
#endif
						++it_count;
					}
				}
			}
			crct_get(&rctinfo, best0, best1, best2);
#if 0
			//PIA13912	B_G008_R002_006
			rctinfo.perm[0]=2;
			rctinfo.perm[1]=1;
			rctinfo.perm[2]=0;
			rctinfo.cu0=8;
			rctinfo.cv0=2;
			rctinfo.cv1=6;
#endif
#ifdef ANALYSIS_BY_ENTROPY
			free(hists);
#endif
			(void)it_count;
			(void)it_select;
			(void)&print_rct2;
			(void)&print_rct;
			//(void)&crct_getparams;
			//(void)&rct_combinations;
			//(void)&rct_names;
#if 0
			int best0=0, best1=0, best2=0;
			for(int k0=0;k0<OCH0_COUNT;++k0)
			{
				for(int k1=0;k1<OCH0_COUNT+OCH1_COUNT;++k1)
				{
					int uch=k1%3;
					if(uch==k0)//disallow _X00_0X1_00X
						continue;
					for(int k2=0;k2<OCH0_COUNT+OCH1_COUNT+OCH2_COUNT;++k2)
					{
						int vch=k2%3;
						if(vch==k0||vch==uch)
						{
#if 0
							crct_getparams(&rctinfo, k0, k1, k2);
							print_rct2(&rctinfo);
							printf("  forbidden\n");
#endif
							continue;
						}
						int64_t val=
							+counters[k0]
							+counters[k1]
							+counters[k2]
						;
						if(!it_count||bestval>val)
							bestval=val, best0=k0, best1=k1, best2=k2, it_select=it_count;
#ifdef PRINT_RCT
#if 1
						crct_getparams(&rctinfo, k0, k1, k2);
						print_rct2(&rctinfo);
						printf("  ");
#endif
						printf("%4d  %12lld +%12lld +%12lld =%12lld%s\n"
							, it_count
							, counters[k0]
							, counters[k1]
							, counters[k2]
							, val
							, it_count==it_select?" <-":""
						);
#endif
						++it_count;
					}
				}
			}
			crct_getparams(&rctinfo, best0, best1, best2);
#endif
#ifdef LOUD
			printf("WH %d*%d  %lld bytes"
				, iw, ih, usize
			);

			printf("  RCT %4d/%4d ", it_select, it_count);
			print_rct2(&rctinfo);

			printf("  dist %d  \"%s\"\n"
				, dist
				, srcfn
			);
#endif
		}
#if 0
		int64_t minerr=0;
		while(imptr<imend)
		{
			int r, g, b, rg, gb, br;

			r=(imptr[0]-imptr[0-rowstride])<<3;
			g=(imptr[1]-imptr[1-rowstride])<<3;
			b=(imptr[2]-imptr[2-rowstride])<<3;
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
				OCH_YX00, r,
				OCH_Y0X0, g,
				OCH_Y00X, b
			);
			UPDATE(
				OCH_CX80, rg,
				OCH_C0X8, gb,
				OCH_C80X, br
			);
			UPDATE(
				OCH_CX44, (rg-br)>>1,//r-(g+b)/2 = (r-g + r-b)/2
				OCH_C4X4, (rg-gb)>>1,//g-(r+b)/2 = (g-r + g-b)/2
				OCH_C44X, (br-gb)>>1 //b-(r+g)/2 = (b-r + b-g)/2
			);
			UPDATE(
				OCH_CX62, rg+(gb>>2),//r-(3*g+b)/4 = r-g-(b-g)/4
				OCH_C6X2, rg+(br>>2),//g-(3*r+b)/4 = g-r-(b-r)/4
				OCH_C62X, br+(rg>>2) //b-(3*r+g)/4 = b-r-(g-r)/4
			);
			UPDATE(
				OCH_CX26, br+(gb>>2),//r-(g+3*b)/4 = r-b-(g-b)/4
				OCH_C2X6, gb+(br>>2),//g-(r+3*b)/4 = g-b-(r-b)/4
				OCH_C26X, gb+(rg>>2) //b-(r+3*g)/4 = b-g-(r-g)/4
			);
			UPDATE(
				OCH_CX71, rg+(gb>>3),//r-(7*g+b)/8 = r-g-(b-g)/8
				OCH_C7X1, rg+(br>>3),//g-(7*r+b)/8 = g-r-(b-r)/8
				OCH_C71X, br+(rg>>3) //b-(7*r+g)/8 = b-r-(g-r)/8
			);
			UPDATE(
				OCH_CX17, br+(gb>>3),//r-(g+7*b)/8 = r-b-(g-b)/8
				OCH_C1X7, gb+(br>>3),//g-(r+7*b)/8 = g-b-(r-b)/8
				OCH_C17X, gb+(rg>>3) //b-(r+7*g)/8 = b-g-(r-g)/8
			);
			UPDATE(
				OCH_CX53, (rg-br-(gb>>2))>>1,//r-(5*g+3*b)/8 =  (5*rg-3*br)/8 = (4*(rg-br)+rg+br)/8 = (((rg+br)>>2)+rg-br)>>1
				OCH_C5X3, (rg-gb-(br>>2))>>1,//g-(5*r+3*b)/8 = -(5*rg-3*gb)/8
				OCH_C53X, (br-gb-(rg>>2))>>1 //b-(5*r+3*g)/8 =  (5*br-3*gb)/8
			);
			UPDATE(
				OCH_CX35, (rg-br+(gb>>2))>>1,//r-(3*g+5*b)/8 =  (3*rg-5*br)/8 = (rg-br-((rg+br)>>2))>>1
				OCH_C3X5, (rg-gb+(br>>2))>>1,//g-(3*r+5*b)/8 = -(3*rg-5*gb)/8
				OCH_C35X, (br-gb+(rg>>2))>>1 //b-(3*r+5*g)/8 =  (3*br-5*gb)/8
			);
			UPDATE(
				OCH_CX33, r-((g+b)*3>>3),//r-3*(g+b)>>3
				OCH_C3X3, g-((r+b)*3>>3),//g-3*(r+b)>>3
				OCH_C33X, b-((r+g)*3>>3) //b-3*(r+g)>>3
			);
			UPDATE(
				OCH_CX60, r-(3*g>>2),
				OCH_C6X0, g-(3*r>>2),
				OCH_C60X, b-(3*r>>2)
			);
			UPDATE(
				OCH_CX06, r-(3*b>>2),
				OCH_C0X6, g-(3*b>>2),
				OCH_C06X, b-(3*g>>2)
			);
			UPDATE(
				OCH_CX70, r-(7*g>>3),
				OCH_C7X0, g-(7*r>>3),
				OCH_C70X, b-(7*r>>3)
			);
			UPDATE(
				OCH_CX07, r-(7*b>>3),
				OCH_C0X7, g-(7*b>>3),
				OCH_C07X, b-(7*g>>3)
			);
			UPDATE(
				OCH_CX50, r-(5*g>>3),
				OCH_C5X0, g-(5*r>>3),
				OCH_C50X, b-(5*r>>3)
			);
			UPDATE(
				OCH_CX05, r-(5*b>>3),
				OCH_C0X5, g-(5*b>>3),
				OCH_C05X, b-(5*g>>3)
			);
			UPDATE(
				OCH_CX51, r-((5*g+b)>>3),
				OCH_C5X1, g-((5*r+b)>>3),
				OCH_C51X, b-((5*r+g)>>3)
			);
			UPDATE(
				OCH_CX15, r-((g+5*b)>>3),
				OCH_C1X5, g-((r+5*b)>>3),
				OCH_C15X, b-((r+5*g)>>3)
			);
			UPDATE(
				OCH_CX42, r-((2*g+b)>>2),
				OCH_C4X2, g-((2*r+b)>>2),
				OCH_C42X, b-((2*r+g)>>2)
			);
			UPDATE(
				OCH_CX24, r-((g+2*b)>>2),
				OCH_C2X4, g-((r+2*b)>>2),
				OCH_C24X, b-((r+2*g)>>2)
			);
#undef  UPDATE
		}
#endif
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
#if 0
		{
			int kt;

#ifdef PRINT_RCT
			for(kt=0;kt<OCH_COUNT;++kt)
				printf("%2d %s %16lld\n", kt, och_names[kt], counters[kt]);
			printf("\n");
#endif
		//	for(kt=0;kt<43;++kt)
		//	for(kt=0;kt<61;++kt)
		//	for(kt=0;kt<79;++kt)
			for(kt=0;kt<RCT_COUNT;++kt)
			{
				const RCTInfo *combination=rct_combinations+kt;
				long long currerr=
					+counters[combination->och[0]]
					+counters[combination->och[1]]
					+counters[combination->och[2]]
				;
				if(!kt||minerr>currerr)
				{
					minerr=currerr;
					bestrct=kt;
				}
#ifdef PRINT_RCT
				printf("RCT%03d %s %16lld +%16lld +%16lld =%16lld%s\n"
					, kt
					, rct_names[kt]
					, counters[combination->och[0]]
					, counters[combination->och[1]]
					, counters[combination->och[2]]
					, currerr
					, kt==bestrct?" <-":""
				);
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
#endif
#ifdef EXPERIMENT
		experiment(image, iw, ih, bestrct);
#endif
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
	//RCTInfo rctinfo=rct_combinations[bestrct];
	//rctinfo.cu0<<=3;
	//rctinfo.cv0<<=3;
	//rctinfo.cv1<<=3;
	if(dist>1)
	{
		if(fwd)
			mainloop(iw, ih, &rctinfo, dist, image, stream, &ac, 1, 1);
		else
			mainloop(iw, ih, &rctinfo, dist, image, stream, &ac, 1, 0);
	}
	else
	{
		if(fwd)
#if 0
		{
			int64_t bestsize=0;
			ACState ac2=ac;
			for(int kt=0;kt<RCT_COUNT;++kt)
			{
				ac=ac2;
				mainloop(iw, ih, kt, dist, image, stream, &ac, 0, 1);
				int64_t csize=ac.ptr-stream;
				if(!kt||bestsize>csize)
					bestsize=csize, bestrct=kt;
#ifdef LOUD
				printf("%03d %s %12lld%s\n", kt, rct_names[kt], csize, bestrct==kt?" <-":"");
#endif
			}
			mainloop(iw, ih, bestrct, dist, image, stream, &ac, 0, 1);
#ifdef LOUD
			printf("\nRCT %03d %s is best at %12lld\n", bestrct, rct_names[bestrct], bestsize);
#endif
		}
#elif 0
		{
			const int perms[]=
			{
				0, 1, 2,
				1, 2, 0,
				2, 0, 1,
				2, 1, 0,
				0, 2, 1,
				1, 0, 2,
			};
			const int cus[]={0, 24, 28, 32, 36, 48, 64};
			const int cvs[]=
			{
				 0,  0,	//sum=0

				 0, 32,	//sum=1/2
				 8, 24,
				16, 16,
				24,  8,
				32,  0,
				
				 0, 40,	//sum=5/8
				10, 30,
				20, 20,
				30, 10,
				40,  0,

				 0, 48,	//sum=3/4
				12, 36,
				24, 24,
				30, 18,
				36, 12,
				42,  6,
				48,  0,

				 0, 56,	//sum=7/8
				14, 42,
				28, 28,
				35, 21,
				42, 14,
				49,  7,
				56,  0,

				 0, 64,	//sum=1
				16, 48,
				32, 32,
				40, 24,
				48, 16,
				56,  8,
				64,  0,
			};
			const int rct_count=(_countof(perms)/3)*_countof(cus)*(_countof(cvs)/2);
			int64_t bestsize=0;
			int bestrctidx=0;
			ACState ac2=ac;
			double t2=time_sec2();
#if 1
			for(int kp=0, it=0;kp<_countof(perms)/3;++kp)
			{
				const int *perm=perms+kp*3;
				for(int ku=0;ku<_countof(cus);++ku)
				{
					for(int kv=0;kv<_countof(cvs)/2;++kv, ++it)
					{
						RCTInfo rct3=
						{
							{0, 0, 0},
							{perm[0], perm[1], perm[2]},
							cus[ku], cvs[2*kv+0], cvs[2*kv+1],
						};
						int64_t csize;

						ac=ac2;
						mainloop(iw, ih, &rct3, dist, image, stream, &ac, 0, 1);
						csize=ac.ptr-stream;
						if(!it||bestsize>csize)
							bestsize=csize, rctinfo=rct3, bestrctidx=it;

						printf("%4d/%4d ", it, rct_count);
						print_rct2(&rct3);
						printf(" %12lld %s  ETA %12.6lf mins\n"
							, csize
							, bestrctidx==it?"<-":"  "
							, (time_sec2()-t2)/(60*(it+1))*(rct_count-(it+1))
						);
					}
				}
			}
#else
			rctinfo.perm[0]=2;
			rctinfo.perm[1]=1;
			rctinfo.perm[2]=0;
			rctinfo.cu0=64;
			rctinfo.cv0=16;
			rctinfo.cv1=48;
#endif
			ac=ac2;
			mainloop(iw, ih, &rctinfo, dist, image, stream, &ac, 0, 1);
			printf("\nBest: ");
			print_rct2(&rctinfo);
			printf(" %12lld\n", bestsize);
		}
#else
			mainloop(iw, ih, &rctinfo, dist, image, stream, &ac, 0, 1);
#endif
		else
			mainloop(iw, ih, &rctinfo, dist, image, stream, &ac, 0, 0);
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
			dstsize+=fwrite(&iw, 1, 3, fdst);
			dstsize+=fwrite(&ih, 1, 3, fdst);
		//	dstsize+=fwrite(&bestrct, 1, 1, fdst);
			dstsize+=fwrite(&rctinfo.perm[0], 1, 1, fdst);
			dstsize+=fwrite(&rctinfo.perm[1], 1, 1, fdst);
			dstsize+=fwrite(&rctinfo.perm[2], 1, 1, fdst);
			dstsize+=fwrite(&rctinfo.cu0, 1, 1, fdst);
			dstsize+=fwrite(&rctinfo.cv0, 1, 1, fdst);
			dstsize+=fwrite(&rctinfo.cv1, 1, 1, fdst);
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
		double total[3]={0}, total2[3]={0}, total_u=0, total_b=0;
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
			total2[0]+=bitsizes[0][kv];
			total2[1]+=bitsizes[1][kv];
			total2[2]+=bitsizes[2][kv];
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
		printf("TYUV %12.2lf %12.2lf %12.2lf %12.2lf  total\n\n"
			, total2[0]+total2[1]+total2[2]
			, total2[0]
			, total2[1]
			, total2[2]
		);
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
	//(void)&squash;
#ifdef PROFILER
	prof_end(prof_ctx);
#endif
	return 0;
}
