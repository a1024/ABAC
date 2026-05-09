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
//	#define PRINT_RCT
//	#define PRINT_ANALYSISSPEED
	#define PRINTBITS4 2

	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif


	#define SUB_LUMAMEAN
//	#define USE_TABLES

//	#define USE_RCPTABLE	//X  slow
//	#define USE_DIV		//X  slow
	#define USE_COUNTERS

//	#define SIGNED_PIXEL


//lossless estims
#if 1
#define PREDLIST\
	PRED(N - 4*dW)\
	PRED(NNWW)\
	PRED(NNW)\
	PRED(NN)\
	PRED(NNE)\
	PRED(NNEE)\
	PRED(NNN)\
	PRED(NNNN)\
	PRED(NNNNN)\
	PRED(NNNNNN)\
	PRED(NNNNNNN)\
	PRED(3*(N-NN)+NNN)\
	PRED(2*N-NN - dNW)\
	PRED(N+NE-NNE)\
	PRED(N+NW-NNW)\
	PRED(N+W-NW + 2*dNW)\
	PRED(N+W-NW + 5*dN)\
	PRED(N+W-NW + 2*dNE)\
	PRED(N+W-NW + 5*dW)\
	PRED(W - 4*dN)\
	PRED(WW)\
	PRED(WWW)\
	PRED(WWWW)\
	PRED(WWWWW)\
	PRED(WWWWWW)\
	PRED(WWWWWWW)\
	PRED(W+NE-N)\
	PRED(W+NE - dNE)\
	PRED(W+NW-NWW)\
	PRED(2*W-WW - dWW)\
	PRED(3*WW-2*WWW)\
	PRED(4*WWW-3*WWWW)\
	PRED(5*W-4*WW)\
	PRED(NWWWW)\
	PRED(NWWW)\
	PRED(NWW)\
	PRED(NW + dNWW)\
	PRED(NE + dN)\
	PRED(NEE + dNEE)\
	PRED(NEEE)\
	PRED(NEEEE)\
	PRED(NEEEEE)\
	PRED(NEEEEEE)\
	PRED(NN+WW-NW)\
	PRED(8*(dN+dW+dNW+dNE))\

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
	PRED(2*W-WW+dW)\

#endif


//lossy estims
#define PREDLIST_LOSSY PREDLIST
#if 0
#define PREDLIST_LOSSY\
	PRED(N)\
	PRED(NNN)\
	PRED(NNNN)\
	PRED(N+N-NN)\
	PRED(N+NE-NNE)\
	PRED(N+NW-NNW)\
	PRED(N+W-NW)\
	PRED(W)\
	PRED(WWW)\
	PRED(WWWW)\
	PRED(W+W-WW)\
	PRED(W+NE-N)\
	PRED(W+NEE-NE)\
	PRED(W+NW-NWW)\
	PRED(NE)\
	PRED(NEEE)\
	PRED(NEEEE)\

#endif
enum
{
	RCT_BITS=3,
	ADDBITS=2,
	
	L1SH=20,
	L1SH_LOSSY=18,
#define PRED(EXPR) +1
	L1NPREDS=PREDLIST,
	L1NPREDS_LOSSY=PREDLIST_LOSSY,
#undef  PRED

	NCTX=24,
	GRBITS=6,
	GRLIMIT=16,

	PROBBITS_USE=14,
#ifdef USE_COUNTERS
	CTRBITS=9,
	CTRMASK=(1<<CTRBITS)-1,
#else
	PROBBITS_STORE=24,
	PROBSHIFT=PROBBITS_STORE-PROBBITS_USE,
#endif

	XPAD=8,
	NCH=3,
	NROWS=8,
	NVAL=4,
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
#define CLAMP2(X, LO, HI) X=X>(LO)?X:LO, X=X<(HI)?X:HI
#define ROUND32(X) _mm_cvt_ss2si(_mm_set_ss(X))
#define ROUND64(X) _mm_cvtsd_si64(_mm_set_sd(X))
#define TRUNC32(X) _mm_cvtt_ss2si(_mm_set_ss(X))
#define TRUNC64(X) _mm_cvttsd_si64(_mm_set_sd(X))
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
} OCHIndex;
static const char *och_names[]=
{
#define OCH(X) #X,
	OCHLIST
#undef  OCH
};
typedef struct _RCTInfo
{
	uint8_t och[3], perm[3];
	int16_t cu0, cv0, cv1;
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
#endif

#if defined USE_COUNTERS && !defined USE_DIV && !defined USE_RCPTABLE
static uint64_t ctrtable[1<<CTRBITS*2];//{p1, next_0(n0, n1}, next_1(n0, n1)}
#endif
#ifdef USE_RCPTABLE
static uint32_t rcptable[2<<CTRBITS];
#endif
typedef uint32_t Cell_t;
static Cell_t stats1[3][256<<ADDBITS][NCTX][GRLIMIT];//unary
static Cell_t stats2[3][256<<ADDBITS][256];//remainder
static Cell_t stats3[3][8][256];//bypass on GRLIMIT
static const size_t memusage=sizeof(stats1)+sizeof(stats2)+sizeof(stats3);
typedef struct _ACState
{
	uint64_t lo, hi, code;
	uint8_t *ptr, *end;

	int64_t renormctr, byteidx;
} ACState;
#ifdef ESTIMATE_BITSIZE
static uint32_t unary_count=0, binary_count=0;
static double shannontable[1<<PROBBITS_USE];
static double bitsizes[3][GRLIMIT+8];
static uint32_t bitctr[3][GRLIMIT+8][2];
static uint32_t winctr[3][GRLIMIT+8];
static int ekc, eidx;
#endif
#ifdef USE_TABLES
static uint8_t epredtable[256];
static uint8_t clamptable[512];
static uint16_t errortable[512][2];
static uint8_t ctxtable[(2<<NCTX/2)/3][2];
#endif
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
INLINE void codebit(ACState *ac, Cell_t *pcell, int32_t *bit, const int fwd)
{
	int rbit;
	Cell_t cell=*pcell;
#ifdef USE_RCPTABLE
	int n[]={cell&CTRMASK, cell>>CTRBITS};
	int sum=n[0]+n[1];
	int rcp=rcptable[sum];
	int num=(n[1]+1)<<PROBBITS_USE;
	int32_t p1=(uint32_t)((uint64_t)num*rcp>>30);
	p1+=p1<1<<PROBBITS_USE>>1;
	//CLAMP2(p1, 1, ((1<<PROBBITS_USE)-1));
#elif defined USE_DIV
	int n[]={cell&CTRMASK, cell>>CTRBITS};
	int m[]={n[0]*8+19, n[1]*8+19};
	int den=m[0]+m[1];
	int num=(m[1]<<PROBBITS_USE)+(den>>1);
	int32_t p1=(uint32_t)num/(uint32_t)den;
	CLAMP2(p1, 1, ((1<<PROBBITS_USE)-1));
#elif defined USE_COUNTERS
	int sh;
	uint64_t entry=ctrtable[cell];
#else
	int32_t p1=cell>>PROBSHIFT;
	p1+=p1<1<<PROBBITS_USE>>1;
#endif
	uint64_t range=ac->hi-ac->lo;
	if(range<=0xFFFF)
	{
#ifdef _MSC_VER
		++ac->renormctr;
#endif
		if(ac->ptr>=ac->end)
		{
#ifdef _MSC_VER
			printf("\n%d/%d  %8.4lf%%\n"
				, (int32_t)(4*ac->renormctr)
				, (int32_t)ac->byteidx
				, 100.*4*ac->renormctr/ac->byteidx
			);
#endif
			CRASH("inflation");
		}
		if(fwd)
			*(uint32_t*)ac->ptr=(uint32_t)(ac->lo>>32);
		else
			ac->code=ac->code<<32|*(uint32_t*)ac->ptr;
		ac->ptr+=4;
		ac->lo<<=32;
		range=range<<32|0xFFFFFFFF;
		
		if(range>~ac->lo)
			range=~ac->lo;
		ac->hi=ac->lo+range;
	}
#ifdef USE_COUNTERS
	int32_t p1=(int32_t)(entry&((1ULL<<PROBBITS_USE)-1));
#endif
#ifdef _MSC_VER
	if((uint32_t)(p1-1)>(uint32_t)((1<<PROBBITS_USE)-2))
		CRASH("Invalid p1 0x%08X / %d bit", p1, PROBBITS_USE);
#endif
	{
		uint64_t mid=ac->lo+(range*(uint32_t)p1>>PROBBITS_USE);
		rbit=*bit;
		rbit=fwd?rbit:ac->code<mid;
		*bit=rbit;
#ifdef FIFOVAL
		if(fwd)
			fifoval_enqueue(rbit<<24^p1);
		else
			fifoval_check(rbit<<24^p1);
#endif
		*(rbit?&ac->hi:&ac->lo)=mid-rbit;
	}
	
#if defined USE_DIV || defined USE_RCPTABLE
	if(n[rbit]>=CTRMASK)
	{
		n[0]>>=1;
		n[1]>>=1;
	}
	++n[rbit];
	cell=n[1]<<CTRBITS|n[0];
#elif defined USE_COUNTERS
	sh=rbit?PROBBITS_USE+CTRBITS*2:PROBBITS_USE;
	cell=(Cell_t)(entry>>sh&((1ULL<<CTRBITS*2)-1));
#else
	cell+=(int32_t)((rbit<<PROBBITS_STORE)-cell+(1<<7>>1))>>7;
#endif
	*pcell=cell;
#ifdef ESTIMATE_BITSIZE
	bitsizes[ekc][eidx]+=shannontable[rbit?p1:(1<<PROBBITS_USE)-p1];
	++bitctr[ekc][eidx][rbit];
	winctr[ekc][eidx]+=rbit==(p1>=1<<PROBBITS_USE>>1);
#endif
}
INLINE void mainloop(int iw, int ih, RCTInfo *rct, int dist, uint8_t *image, uint8_t *stream, ACState *ac, const int lossy, const int fwd, uint16_t *lumamean)
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
	ALIGN(32) int32_t lcoeffs[3][L1NPREDS_LOSSY+1]={0};
	ALIGN(32) int32_t coeffs[3][L1NPREDS]={0}, bias[3]={0};
	int32_t invdist=((1<<16)+dist-1)/dist;
	uint8_t *imptr=image;
#ifdef PRINTBITS4
	int32_t printidx=0;
#endif
	uint16_t offset_bias[]=
	{
					lumamean[rct->perm[0]],
		rct->cu0?0:		lumamean[rct->perm[1]],
		rct->cv0+rct->cv1?0:	lumamean[rct->perm[2]],
	//				lumamean[rct->perm[0]]&~((1<<RCT_BITS>>1)-1),
	//	rct->cu0?0:		lumamean[rct->perm[1]]&~((1<<RCT_BITS>>1)-1),
	//	rct->cv0+rct->cv1?0:	lumamean[rct->perm[2]]&~((1<<RCT_BITS>>1)-1),
	};
	
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
#ifdef USE_RCPTABLE
	for(int k=0;k<2<<CTRBITS;++k)
		rcptable[k]=(1<<30)/(k+2);
#endif
#if defined USE_COUNTERS && !defined USE_DIV && !defined USE_RCPTABLE
	for(int k=0;k<1<<CTRBITS*2;++k)//MSB {next_1(n0, n1), next_0(n0, n1}, p1} LSB
	{
		int32_t n[]={k&CTRMASK, k>>CTRBITS&CTRMASK};
		int n0e=n[0]*8+19;
		int n1e=n[1]*8+19;
		int sum=n0e+n1e;
		int32_t p1=((n1e<<PROBBITS_USE)+(sum>>1))/sum;
		uint64_t cell=0;
		int32_t n0[]={k&CTRMASK, k>>CTRBITS&CTRMASK};
		int32_t n1[]={k&CTRMASK, k>>CTRBITS&CTRMASK};
		CLAMP2(p1, 1, ((1<<PROBBITS_USE)-1));
		
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
	for(ky=0;ky<ih;++ky)
	{
		int yuv[3]={0};
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-4LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-5LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-6LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-7LL+NROWS)%NROWS)*NVAL,
		};
		for(kx=0;kx<iw;++kx, imptr+=3)
		{
			int kc;
			int offset, offset0;

			if(fwd)
			{
#ifdef SIGNED_PIXEL
				yuv[0]=imptr[yidx]-128;
				yuv[1]=imptr[uidx]-128;
				yuv[2]=imptr[vidx]-128;
#else
				yuv[0]=imptr[yidx];
				yuv[1]=imptr[uidx];
				yuv[2]=imptr[vidx];
#endif
			}
			offset0=offset=0;
			for(kc=0;kc<3;++kc)
			{
				int32_t
					NNNNNNN	=rows[7][0+0*NCH*NROWS*NVAL],
					NNNNNN	=rows[6][0+0*NCH*NROWS*NVAL],
					NNNNN	=rows[5][0+0*NCH*NROWS*NVAL],
					NNNN	=rows[4][0+0*NCH*NROWS*NVAL],
					NNNNEE	=rows[4][0+2*NCH*NROWS*NVAL],
					NNNW	=rows[3][0-1*NCH*NROWS*NVAL],
					NNN	=rows[3][0+0*NCH*NROWS*NVAL],
					NNNE	=rows[3][0+1*NCH*NROWS*NVAL],
					NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
					NNW	=rows[2][0-1*NCH*NROWS*NVAL],
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL],
					NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
					NWWWWW	=rows[1][0-5*NCH*NROWS*NVAL],
					NWWWW	=rows[1][0-4*NCH*NROWS*NVAL],
					NWWW	=rows[1][0-3*NCH*NROWS*NVAL],
					NWW	=rows[1][0-2*NCH*NROWS*NVAL],
					NW	=rows[1][0-1*NCH*NROWS*NVAL],
					N	=rows[1][0+0*NCH*NROWS*NVAL],
					NE	=rows[1][0+1*NCH*NROWS*NVAL],
					NEE	=rows[1][0+2*NCH*NROWS*NVAL],
					NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
					NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
					NEEEEE	=rows[1][0+5*NCH*NROWS*NVAL],
					NEEEEEE	=rows[1][0+6*NCH*NROWS*NVAL],
					WWWWWWW	=rows[0][0-7*NCH*NROWS*NVAL],
					WWWWWW	=rows[0][0-6*NCH*NROWS*NVAL],
					WWWWW	=rows[0][0-5*NCH*NROWS*NVAL],
					WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
					WWW	=rows[0][0-3*NCH*NROWS*NVAL],
					WW	=rows[0][0-2*NCH*NROWS*NVAL],
					W	=rows[0][0-1*NCH*NROWS*NVAL],
					
					dNNNN	=rows[4][1+0*NCH*NROWS*NVAL],
					dNNN	=rows[3][1+0*NCH*NROWS*NVAL],
					dNNWW	=rows[2][1-2*NCH*NROWS*NVAL],
					dNNW	=rows[2][1-1*NCH*NROWS*NVAL],
					dNN	=rows[2][1+0*NCH*NROWS*NVAL],
					dNNE	=rows[2][1+1*NCH*NROWS*NVAL],
					dNNEE	=rows[2][1+2*NCH*NROWS*NVAL],
					dNWWW	=rows[1][1-3*NCH*NROWS*NVAL],
					dNWW	=rows[1][1-2*NCH*NROWS*NVAL],
					dNW	=rows[1][1-1*NCH*NROWS*NVAL],
					dN	=rows[1][1+0*NCH*NROWS*NVAL],
					dNE	=rows[1][1+1*NCH*NROWS*NVAL],
					dNEE	=rows[1][1+2*NCH*NROWS*NVAL],
					dNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
					dNEEEE	=rows[1][1+4*NCH*NROWS*NVAL],
					dWWWW	=rows[0][1-4*NCH*NROWS*NVAL],
					dWWW	=rows[0][1-3*NCH*NROWS*NVAL],
					dWW	=rows[0][1-2*NCH*NROWS*NVAL],
					dW	=rows[0][1-1*NCH*NROWS*NVAL],

					xN	=rows[1][2+0*NCH*NROWS*NVAL],
					xNE	=rows[1][2+1*NCH*NROWS*NVAL],
					xNEE	=rows[1][2+2*NCH*NROWS*NVAL],
					xNEEE	=rows[1][2+3*NCH*NROWS*NVAL],
					xW	=rows[0][2-1*NCH*NROWS*NVAL],

					yN	=rows[1][3+0*NCH*NROWS*NVAL],
					yNE	=rows[1][3+1*NCH*NROWS*NVAL],
					yNEE	=rows[1][3+2*NCH*NROWS*NVAL],
					yNEEE	=rows[1][3+3*NCH*NROWS*NVAL],
					yW	=rows[0][3-1*NCH*NROWS*NVAL];
				int64_t pred;
				int32_t tidx=0;
				int32_t bit=0;
				ALIGN(32) int32_t preds[L1NPREDS>L1NPREDS_LOSSY?L1NPREDS:L1NPREDS_LOSSY];
				int32_t pred0, upred2;
				int32_t error;
				int32_t nbypass, nbypass0;
				int32_t nzeros=0, grflag;
				Cell_t *statsptr;
				int32_t ctx;
				int epred;
				int j;

				if(kc==0)
					offset0=offset_bias[0];
				else if(kc==1)
					offset0=cu0*yuv[0]+offset_bias[1];
				else 
					offset0=cv0*yuv[0]+cv1*yuv[1]+offset_bias[2];
				offset=offset0>>RCT_BITS;
				{
					if(lossy)
					{
						pred=1<<L1SH_LOSSY>>1;
#define PRED(EXPR) preds[j]=EXPR; pred+=lcoeffs[kc][j]*preds[j]; ++j;
						j=0;
						PREDLIST_LOSSY;
#undef  PRED
						upred2=(int32_t)(pred>>(L1SH_LOSSY-ADDBITS));
					}
					else
					{
						pred=(int64_t)bias[kc];
#define PRED(EXPR) preds[j]=EXPR; pred+=coeffs[kc][j]*preds[j]; ++j;
						j=0;
						PREDLIST;
#undef  PRED
						upred2=(int32_t)(pred>>(L1SH-ADDBITS));
					}
				}
				pred0=upred2>>ADDBITS;
				upred2+=offset0<<ADDBITS>>RCT_BITS;
				CLAMP2(upred2, 0, (256<<ADDBITS)-1);
				pred=upred2>>ADDBITS;
#ifdef USE_TABLES
				uint8_t *ctxptr=ctxtable[eW<_countof(ctxtable)-1?eW:_countof(ctxtable)-1];
				ctx=ctxptr[0];
				nbypass=ctxptr[1];
				//ctx=abs(N-W)>>3;
				//if(ctx>NCTX-1)
				//	ctx=NCTX-1;
#else
				ctx=31-_lzcnt_u32(xW*xW+1);
			//	ctx=31-_lzcnt_u32(xW*yW+1);
			//	ctx=31-_lzcnt_u32(xN*yW+1);
				nbypass=(ctx>>1)-GRBITS;
				CLAMP2(nbypass, 0, 7);
				if(ctx>NCTX-1)
					ctx=NCTX-1;
#endif
#if 1
				(void)NNNNNNN	;
				(void)NNNNNN	;
				(void)NNNNN	;
				(void)NNNN	;
				(void)NNNNEE	;
				(void)NNNW	;
				(void)NNN	;
				(void)NNNE	;
				(void)NNWW	;
				(void)NNW	;
				(void)NN	;
				(void)NNE	;
				(void)NNEE	;
				(void)NWWWWW	;
				(void)NWWWW	;
				(void)NWWW	;
				(void)NWW	;
				(void)NW	;
				(void)N		;
				(void)NE	;
				(void)NEE	;
				(void)NEEE	;
				(void)NEEEE	;
				(void)NEEEEE	;
				(void)NEEEEEE	;
				(void)WWWWWWW	;
				(void)WWWWWW	;
				(void)WWWWW	;
				(void)WWWW	;
				(void)WWW	;
				(void)WW	;
				(void)W		;

				(void)dNNNN	;
				(void)dNNN	;
				(void)dNNWW	;
				(void)dNNW	;
				(void)dNN	;
				(void)dNNE	;
				(void)dNNEE	;
				(void)dNWWW	;
				(void)dNWW	;
				(void)dNW	;
				(void)dN	;
				(void)dNE	;
				(void)dNEE	;
				(void)dNEEE	;
				(void)dNEEEE	;
				(void)dWWWW	;
				(void)dWWW	;
				(void)dWW	;
				(void)dW	;

				(void)xN	;
				(void)xNE	;
				(void)xNEE	;
				(void)xNEEE	;
				(void)xW	;

				(void)yN	;
				(void)yNE	;
				(void)yNEE	;
				(void)yNEEE	;
				(void)yW	;
#endif
#ifdef USE_TABLES
				epred=epredptr[pred];
#else
#ifdef SIGNED_PIXEL
				if(lossy)
					epred=pred>=0?127-pred:pred+128;
				else
					epred=128-abs(pred);
#else
				if(lossy)
					epred=(int32_t)pred>=128?255-(int32_t)pred:(int32_t)pred;
				else
					epred=128-abs((int32_t)pred-128);
#endif
#endif
				if(fwd)
				{
					//if(ky==ih/2&&kx==iw/2&&kc==0)//
					//	printf("");
					if(lossy)
					{
						int pixel;
						error=yuv[kc]-(int32_t)pred;
						epred=epred*invdist>>16;
						error=(error*invdist>>16)-(error>>31);
						pixel=error*dist+(int32_t)pred;
#ifdef SIGNED_PIXEL
						CLAMP2(pixel, -128, 127);
						yuv[kc]=pixel;
						pixel+=128;
#else
#ifdef USE_TABLES
						pixel=clampptr[pixel];
#else
						CLAMP2(pixel, 0, 255);
#endif
						yuv[kc]=pixel;
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
#ifdef PRINTBITS4
					if(fwd&&!kc&&(uint32_t)(printidx-100000)<100000&&nbypass==PRINTBITS4)//
						printf("%c", '0'+bit);
#endif
					if(bit)
						break;
					++tidx;
				}while(tidx<GRLIMIT);
				nbypass0=nbypass;
				grflag=tidx==GRLIMIT;
				if(grflag)
				{
					error-=GRLIMIT<<nbypass;
					statsptr=stats3[kc][nbypass];
					switch(nbypass)
					{
					case 0:nbypass=8;break;
					case 1:nbypass=8;break;
					case 2:nbypass=8;break;
					case 3:nbypass=7;break;
					case 4:nbypass=1;break;
					case 5:nbypass=1;break;
					case 6:nbypass=1;break;
					case 7:nbypass=1;break;
					}
				//	nbypass=8;
					tidx=0;
				}
				else
					statsptr=stats2[kc][upred2];
				tidx+=256>>nbypass;//bit coding:  tidx=2*tidx+bit  tidx=0b1XX
#ifdef PRINTBITS4
				if(fwd&&!kc&&(uint32_t)(printidx-100000)<100000&&nbypass&&nbypass==PRINTBITS4)//
					printf(" ");
#endif
				{
					int32_t kb=nbypass-1;

					for(;kb>=0;--kb)
					{
						//if(fwd)
							bit=error>>kb&1;
#ifdef ESTIMATE_BITSIZE
						eidx=GRLIMIT+8-nbypass+kb;
						++binary_count;
#endif
						codebit(ac, statsptr+tidx, &bit, fwd);
#ifdef PRINTBITS4
						if(fwd&&!kc&&(uint32_t)(printidx-100000)<100000&&nbypass==PRINTBITS4)//
							printf("%c", bit?'T':'F');
#endif
						tidx=2*tidx+bit;
					}
				}
#ifdef PRINTBITS4
				if(fwd&&!kc&&(uint32_t)(printidx-100000)<100000&&nbypass==PRINTBITS4)//
					printf("%c", printidx&63?' ':'\n');
				printidx+=!kc;
#endif
				if(grflag)
					tidx+=GRLIMIT<<nbypass0;
#ifdef _MSC_VER
				if(fwd&&grflag)
					error+=GRLIMIT<<nbypass0;
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
#ifdef SIGNED_PIXEL
						negmask=(int32_t)pred>>31;
#else
						negmask=((int32_t)pred-128)>>31;
#endif
						sym=error;
						e2=epred-sym;
						error=sym>>1^-(sym&1);
						e2=(e2^negmask)-negmask;
						if((epred<<1)<sym)
							error=e2;

						pixel=error*dist+(int32_t)pred;
#ifdef SIGNED_PIXEL
						CLAMP2(pixel, -128, 127);
						yuv[kc]=pixel;
						pixel+=128;
#else
						CLAMP2(pixel, 0, 255);
						yuv[kc]=pixel;
#endif
					}
					else
					{
#ifdef SIGNED_PIXEL
						if(2*pred+error==256)
#else
						//if(2*(pred-128)+error==256)
						if(2*pred+error==512)
#endif
						{
							error=error>>1^-(error&1);
#ifdef SIGNED_PIXEL
							yuv[kc]=(int8_t)(error+pred);
#else
							yuv[kc]=(uint8_t)(error+pred);
#endif
						}
						else
						{
#ifdef SIGNED_PIXEL
							int negmask=(int32_t)pred>>31;
#else
							int negmask=((int32_t)pred-128)>>31;
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
#ifdef SIGNED_PIXEL
						pixel+=128;
#endif
						if(pixel!=val)
							CRASH("GUIDE YXC %d %d %d", ky, kx, kc);
					}
#endif
				}
				{
					int32_t curr=yuv[kc]-offset;
					
					rows[0][0]=curr;

					error=(yuv[kc]<<ADDBITS)-upred2;
					rows[0][1]=error;

					error=yuv[kc]-(int32_t)pred;
					if(lossy)
						error=abs(error);
					else
						error=error<<1^error>>31;
					error<<=GRBITS;
					rows[0][2]=(xW+(xW<xNE?xW:xNE)+error+(xNEE>xNEEE?xNEE:xNEEE))>>2;
				//	rows[0][2]=xW+((error-xW)>>2);
				//	rows[0][3]=yNE+((error-yNE)>>2);
				//	rows[0][2]=(16*eW+7*(error<<GRBITS)+9*(eNEE>eNEEE?eNEE:eNEEE))>>5;
				//	rows[0][2]=(2*eW+(error<<GRBITS)+(eNEE>eNEEE?eNEE:eNEEE))>>2;

					//if(curr-pred0!=yuv[kc]-(int32_t)pred)//
					//	printf("");
					
					error=(curr>pred0)-(curr<pred0);
					if(lossy)
					{
						//currw[L1NPREDS_LOSSY]+=e;
#define PRED(EXPR) lcoeffs[kc][j]+=error*preds[j]; ++j;
						j=0;
						PREDLIST_LOSSY;
#undef  PRED
					}
					else
					{
						bias[kc]+=error<<9;
					//	bias[kc]+=(curr-pred0)<<12;
#define PRED(EXPR) coeffs[kc][j]+=error*preds[j]; ++j;
						j=0;
						PREDLIST;
#undef  PRED
					}
				}
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				rows[4]+=NROWS*NVAL;
				rows[5]+=NROWS*NVAL;
				rows[6]+=NROWS*NVAL;
				rows[7]+=NROWS*NVAL;
#ifdef _MSC_VER
				++ac->byteidx;
#endif
			}
			if(!fwd)
			{
#if defined SIGNED_PIXEL
				imptr[yidx]=yuv[0]+128;
				imptr[uidx]=yuv[1]+128;
				imptr[vidx]=yuv[2]+128;
#else
				imptr[yidx]=yuv[0];
				imptr[uidx]=yuv[1];
				imptr[vidx]=yuv[2];
#endif
#ifdef ENABLE_GUIDE
				guide_check(image, kx, ky);
#endif
			}
		}
	}
	free(pixels);
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
	ACState ac={0};
#ifdef SUB_LUMAMEAN
	uint16_t lumamean[3]={0};
#endif
#ifdef LOUD
	double t=time_sec2();
#endif
#ifdef PROFILER
	void *prof_ctx=prof_start();
#endif

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
	
	ac.hi=0xFFFFFFFFFFFF;
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
			dist=1;
			fread(&iw, 1, 3, fsrc);
			fread(&ih, 1, 3, fsrc);
			fread(&rctinfo.perm[0], 1, 1, fsrc);
			fread(&rctinfo.perm[1], 1, 1, fsrc);
			fread(&rctinfo.perm[2], 1, 1, fsrc);
			fread(&rctinfo.cu0, 1, 2, fsrc);
			fread(&rctinfo.cv0, 1, 2, fsrc);
			fread(&rctinfo.cv1, 1, 2, fsrc);
#ifdef SUB_LUMAMEAN
			fread(&lumamean, 1, sizeof(lumamean), fsrc);
#endif
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
#ifdef PRINT_ANALYSISSPEED
		double t3=time_sec2();
#endif
		int64_t counters[OCH_COUNT]={0};
		ALIGN(32) int prev[OCH_COUNT]={0};
		int rowstride=3*iw;
#ifdef SUB_LUMAMEAN
		int64_t mean[3]={0};
		for(imptr=image;imptr<imend;imptr+=3)
		{
			mean[0]+=imptr[0];
			mean[1]+=imptr[1];
			mean[2]+=imptr[2];
		}
		uint64_t res=(uint64_t)iw*ih;
		lumamean[0]=(uint16_t)(((mean[0]<<RCT_BITS)+(res>>1))/res);
		lumamean[1]=(uint16_t)(((mean[1]<<RCT_BITS)+(res>>1))/res);
		lumamean[2]=(uint16_t)(((mean[2]<<RCT_BITS)+(res>>1))/res);
		//lumamean[0]&=~((1<<RCT_BITS>>1)-1);
		//lumamean[1]&=~((1<<RCT_BITS>>1)-1);
		//lumamean[2]&=~((1<<RCT_BITS>>1)-1);
#endif

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
				for(int k=0;k<OCH_COUNT;++k)
					counters[k]+=abs(och[k]-prev[k]);
				for(int k=0;k<OCH_COUNT;++k)
					prev[k]=och[k];
			}
		}
#ifdef PRINT_RCT
		for(int kc=0;kc<OCH_COUNT;kc+=3)
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
		{
			static const int perms[]=
			{
				0, 1, 2,
				1, 2, 0,
				2, 0, 1,
				1, 0, 2,
				0, 2, 1,
				2, 1, 0,
			};
			int64_t bestval=0;
			int it_count=0, it_select=0;
			int best0=0, best1=0, best2=0;
#ifdef PRINT_RCT
			RCTInfo rctinfo2={0};
#endif

			for(int kp=0;kp<_countof(perms)/3;++kp)
			{
				int yidx=perms[3*kp+0];
				int uidx=perms[3*kp+1];
				int vidx=perms[3*kp+2];
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
						int64_t val=
							+counters[c0]
							+counters[c1]
							+counters[c2]
						;
						if(!it_count||bestval>val)
						{
							bestval=val;
							best0=c0;
							best1=c1;
							best2=c2;
							it_select=it_count;
						}
#ifdef PRINT_RCT
#if 1
						crct_get(&rctinfo2, c0, c1, c2);
						print_rct2(&rctinfo2);
						printf(" %d %3d %3d _%s_%s_%s  "
							, kp, k1, k2
							, och_names[c0]
							, och_names[c1]
							, och_names[c2]
						);
#endif
						printf("%4d  %12lld +%12lld +%12lld =%12lld%s\n"
							, it_count
							, counters[c0]
							, counters[c1]
							, counters[c2]
							, val
							, it_count==it_select?" <-":""
						);
#endif
						++it_count;
					}
				}
			}
			crct_get(&rctinfo, best0, best1, best2);
#ifdef PRINT_ANALYSISSPEED
			t3=time_sec2()-t3;
			printf("Analysis 1  %12.6lf sec  %12.6lf MB/s\n"
				, t3
				, usize/(t3*1024*1024)
			);
#endif
#if 0
			//PIA13912	B_G008_R002_006
			rctinfo.perm[0]=2;
			rctinfo.perm[1]=1;
			rctinfo.perm[2]=0;
			rctinfo.cu0=8;
			rctinfo.cv0=2;
			rctinfo.cv1=6;
#endif
			(void)it_count;
			(void)it_select;
			(void)&print_rct2;
			(void)&print_rct;
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
		
		streamptr=stream;
		streamend=stream+usize;
	}
	else
	{
		streamptr=stream;
		streamend=stream+srcsize;

		ac.code=*(uint64_t*)streamptr;//load
		streamptr+=sizeof(uint64_t);
		ac.code=ac.code<<32|ac.code>>32;

		csize=srcsize;
	}
	ac.ptr=streamptr;
	ac.end=streamend;
	if(dist>1)
	{
		if(fwd)
			mainloop(iw, ih, &rctinfo, dist, image, stream, &ac, 1, 1, lumamean);
		else
			mainloop(iw, ih, &rctinfo, dist, image, stream, &ac, 1, 0, lumamean);
	}
	else
	{
		if(fwd)
#if 0
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
			mainloop(iw, ih, &rctinfo, dist, image, stream, &ac, 0, 1, lumamean);
#endif
		else
			mainloop(iw, ih, &rctinfo, dist, image, stream, &ac, 0, 0, lumamean);
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
			*(uint64_t*)ac.ptr=ac.lo<<32|ac.lo>>32;//flush
			ac.ptr+=sizeof(uint64_t);

			csize=ac.ptr-stream;

			dstsize+=fwrite(&tag, 1, 2, fdst);
			dstsize+=fwrite(&iw, 1, 3, fdst);
			dstsize+=fwrite(&ih, 1, 3, fdst);
			dstsize+=fwrite(&rctinfo.perm[0], 1, 1, fdst);
			dstsize+=fwrite(&rctinfo.perm[1], 1, 1, fdst);
			dstsize+=fwrite(&rctinfo.perm[2], 1, 1, fdst);
			dstsize+=fwrite(&rctinfo.cu0, 1, 2, fdst);
			dstsize+=fwrite(&rctinfo.cv0, 1, 2, fdst);
			dstsize+=fwrite(&rctinfo.cv1, 1, 2, fdst);
#ifdef SUB_LUMAMEAN
			dstsize+=fwrite(&lumamean, 1, sizeof(lumamean), fdst);
#endif
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
#ifdef ESTIMATE_BITSIZE
		double total[3]={0}, total2[3]={0}, total_u=0, total_b=0;
		printf("plane  csize / usize = invCR  nzeros%%\n");
		for(int kv=0;kv<GRLIMIT+8;++kv)
		{
			uint32_t sum0=bitctr[0][kv][0]+bitctr[0][kv][1];
			uint32_t sum1=bitctr[1][kv][0]+bitctr[1][kv][1];
			uint32_t sum2=bitctr[2][kv][0]+bitctr[2][kv][1];
			printf("%3d %12.2lf/%12.2lf=%8.4lf%% %8.4lf%%  %12.2lf/%12.2lf=%8.4lf%% %8.4lf%%  %12.2lf/%12.2lf=%8.4lf%% %8.4lf%%  %8.4lf%% %8.4lf%% %8.4lf%%\n"
				, kv
				, bitsizes[0][kv], sum0/8., 800.*bitsizes[0][kv]/(sum0?sum0:1), 100.*bitctr[0][kv][1]/(sum0?sum0:1)
				, bitsizes[1][kv], sum1/8., 800.*bitsizes[1][kv]/(sum1?sum1:1), 100.*bitctr[1][kv][1]/(sum1?sum1:1)
				, bitsizes[2][kv], sum2/8., 800.*bitsizes[2][kv]/(sum2?sum2:1), 100.*bitctr[2][kv][1]/(sum2?sum2:1)
				, 100.*winctr[0][kv]/(sum0?sum0:1)
				, 100.*winctr[1][kv]/(sum1?sum1:1)
				, 100.*winctr[2][kv]/(sum2?sum2:1)
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
	(void)dstsize;
	(void)csize;
	(void)&time_sec2;
	//(void)prednames;
	//(void)&print_weights;
	//(void)&squash;
#ifdef PROFILER
	prof_end(prof_ctx);
#endif
	return 0;
}
