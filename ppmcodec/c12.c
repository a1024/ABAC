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

	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif


#define PREDLIST\
	PRED(3*dN-14*dW)\
	PRED(3*dW-14*dN)\
	PRED(-8*dNE)\
	PRED(-8*dNW)\
	PRED(N+W-NW + 2*dNW)\
	PRED(N+W-NW + 5*dN)\
	PRED(N+W-NW + 2*dNE)\
	PRED(N+W-NW + 5*dW)\
	PRED(2*N-NN - dNW)\
	PRED(W+NE-N - dNE)\
	PRED(2*W-WW + 2*dWW)\
	PRED(NW + dNWW)\
	PRED(NE + 2*dN)\
	PRED(NEE + dNEE)\
	PRED(NNNNN)\
	PRED(NNNN)\
	PRED(NNN)\
	PRED(NNWWWW)\
	PRED(NNWWW)\
	PRED(NNWW)\
	PRED(NNW)\
	PRED(NNE)\
	PRED(NNEE)\
	PRED(NNEEE)\
	PRED(NWWWW)\
	PRED(NWWW)\
	PRED(NEEE)\
	PRED(NEEEE)\
	PRED(NEEEEE)\
	PRED(WWWWW)\
	PRED(WWWW)\
	PRED(WWW)\
	PRED(3*(N-NN)+NNN)\
	PRED(W+NW-NWW)\
	PRED(3*WW-2*WWW)\
	PRED(4*WWW-3*WWWW)\
	PRED(5*W-4*WW)\
	PRED(N+NE-NNE)\
	PRED(N+NW-NNW)\
	PRED(NN+WW-NW)\


enum
{
	RCT_BITS=3,
	ADDBITS=2,
	
#define PRED(EXPR) +1
	L1NPREDS=PREDLIST,
#undef  PRED

	CTXMIN=9,
	CTXMAX=18,
	NCTX=CTXMAX+2-CTXMIN,
	GRBITS=5,
	GRLIMIT=18,

	XPAD=8,
	NCH=3,
	NROWS=8,
	NVAL=4,
};


//runtime
#if 1
#ifdef _MSC_VER
#	define INLINE __forceinline static
#else
#	define INLINE __attribute__((always_inline)) inline static
#endif
#ifndef ALIGN
#	ifdef _MSC_VER
#		define	ALIGN(N) __declspec(align(N))
#	else
#		define	ALIGN(N) __attribute__((aligned(N)))
#	endif
#endif
#ifdef _MSC_VER
#	define LZCNT32 _lzcnt_u32
#	define LZCNT64 _lzcnt_u64
#	define TZCNT32 _tzcnt_u32
#	define TZCNT64 _tzcnt_u64
#else
#	define LZCNT32 __builtin_clz
#	define LZCNT64 __builtin_clzll
#	define TZCNT32 __builtin_ctz
#	define TZCNT64 __builtin_ctzll
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


typedef int32_t Cell_t;
ALIGN(64) static Cell_t stats1[3][NCTX][256<<ADDBITS][GRLIMIT];//unary
ALIGN(64) static Cell_t stats2[3][256<<ADDBITS][256];//remainder
ALIGN(64) static Cell_t stats3[3][8][256];//bypass on GRLIMIT
static const size_t memusage=sizeof(stats1)+sizeof(stats2)+sizeof(stats3);
typedef struct _ACState
{
	uint64_t lo, hi, code;
	uint8_t *ptr, *end, *start;

	int64_t byteidx;
} ACState;
static int squash(int x)
{
	enum
	{
		USEBITS=12,
	};
	static const int t[33]=//2^5 table elements, table amplitude 2^12
	{
		   1,    2,    3,    6,   10,   16,   27,   45,   73,  120,  194,
		 310,  488,  747, 1101, 1546, 2047, 2549, 2994, 3348, 3607, 3785,
		3901, 3975, 4022, 4050, 4068, 4079, 4085, 4089, 4092, 4093, 4094,
	};
	int w=x&((1<<(USEBITS-5))-1);
	x=(x>>(USEBITS-5))+16;
	if(x>31)
		return (1<<USEBITS)-1;
	if(x<0)
		return 1;
	x=(t[x]*((1<<(USEBITS-5))-w)+t[x+1]*w+64)>>(12-5);
	return x;
}
static int stretch(int x)//ln(x/(1-x))		probs -> logits
{
	enum
	{
		USEBITS=12,
	};
	static short t[4096];
	static int initialized=0;
	if(!initialized)
	{
		initialized=1;
		
		int pi=0;
		for(int k=-2047;k<=2047;++k)//invert squash()
		{
			int i=squash(k<<(USEBITS-12))>>(USEBITS-12);
			for(int j=pi;j<=i;++j)
				t[j]=k<<(USEBITS-12);
			pi=i+1;
		}
		t[4095]=(1<<USEBITS>>1)-1;
	}
	x=t[x>>(USEBITS-12)];
	return x;
}
INLINE void codebit(ACState *ac, Cell_t *pcell, int32_t *pbit, const int fwd)
{
	enum
	{
		USEBITS=10,
		STOREBITS=20,
		AGEBITS=9,
		AGEMAX=(1<<AGEBITS)-1,
	};
	uint64_t x;
	int32_t prob, count, p1, sh;
	Cell_t cell;
	int bit;

	cell=*pcell;
	x=ac->hi-ac->lo;
	sh=p1=(uint32_t)(cell>>(STOREBITS-USEBITS+AGEBITS));
	count=(uint32_t)(cell&AGEMAX);
	CLAMP2(sh, -48, 48);
	prob=(uint32_t)(cell>>AGEBITS);
	p1-=sh>>4;
	sh=31^LZCNT32(count+10);
	count+=count<AGEMAX;
	p1+=1<<USEBITS>>1;
	if(x<=0xFFFF)
	{
		if(ac->ptr>=ac->end)
		{
#ifdef _MSC_VER
			int64_t usize=ac->ptr-ac->start;

			printf("\n%d/%d  %8.4lf%%\n"
				, (int32_t)usize
				, (int32_t)ac->byteidx
				, 100.*usize/ac->byteidx
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
		ac->hi=ac->hi<<32|~0U;
		if(ac->hi<ac->lo)ac->hi=~0ULL;
		x=ac->hi-ac->lo;
	}
	x=ac->lo+(x*(uint32_t)p1>>USEBITS);
	bit=*pbit;
	bit=fwd?bit:ac->code<x;
	*pbit=bit;
#ifdef _MSC_VER
	if((uint32_t)(p1-1)>(uint32_t)((1<<USEBITS)-2))
		CRASH("Invalid p1 0x%08X / %d bit", p1, USEBITS);
#ifdef FIFOVAL
	if(fwd)
		fifoval_enqueue(rbit<<24^p1);
	else
		fifoval_check(rbit<<24^p1);
#endif
#endif
	*(bit?&ac->hi:&ac->lo)=x-bit;
	prob+=(((2*bit-1)<<(STOREBITS-1))-prob)>>sh;
	*pcell=prob<<AGEBITS|count;
	(void)&squash;
	(void)&stretch;
}
INLINE void mainloop(int iw, int ih, RCTInfo *rct, uint8_t *image, uint8_t *stream, ACState *ac, const int fwd, uint16_t *lumamean, uint8_t *lumamin,  uint8_t *lumamax)
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
	ALIGN(32) int32_t coeffs[3][L1NPREDS]={0}, bias[3]={0};
	uint8_t *imptr=image;
	uint16_t offset_bias[]=
	{
					lumamean[rct->perm[0]],
		rct->cu0?0:		lumamean[rct->perm[1]],
		rct->cv0+rct->cv1?0:	lumamean[rct->perm[2]],
	};
	int64_t res=0;
	int sh=0;

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
	memset(stats1, 0, sizeof(stats1));
	memset(stats2, 0, sizeof(stats2));
	memset(stats3, 0, sizeof(stats3));
	res=(int64_t)iw*ih;
//	if(res<        10000)	sh=17;
//	else if(res<  300000)	sh=18;
//	else if(res< 3000000)	sh=19;
	if(res<     20000000)	sh=20;
	else if(res<40000000)	sh=21;
	else			sh=22;
	FILLMEM_S((int32_t*)coeffs, (1<<sh)/L1NPREDS, sizeof(coeffs), sizeof(int32_t));
	bias[0]=1<<sh>>1;
	bias[1]=1<<sh>>1;
	bias[2]=1<<sh>>1;
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
				yuv[0]=imptr[yidx];
				yuv[1]=imptr[uidx];
				yuv[2]=imptr[vidx];
			}
			offset0=offset=0;
			for(kc=0;kc<3;++kc)
			{
				int32_t
					NNNNNNN	=rows[7][0+0*NCH*NROWS*NVAL],
					NNNNNNNE=rows[7][0+1*NCH*NROWS*NVAL],
					NNNNNN	=rows[6][0+0*NCH*NROWS*NVAL],
					NNNNNNE	=rows[6][0+1*NCH*NROWS*NVAL],
					NNNNN	=rows[5][0+0*NCH*NROWS*NVAL],
					NNNNNE	=rows[5][0+1*NCH*NROWS*NVAL],
					NNNNWW	=rows[4][0-2*NCH*NROWS*NVAL],
					NNNN	=rows[4][0+0*NCH*NROWS*NVAL],
					NNNNE	=rows[4][0+1*NCH*NROWS*NVAL],
					NNNNEE	=rows[4][0+2*NCH*NROWS*NVAL],
					NNNWWW	=rows[3][0-3*NCH*NROWS*NVAL],
					NNNWW	=rows[3][0-2*NCH*NROWS*NVAL],
					NNNW	=rows[3][0-1*NCH*NROWS*NVAL],
					NNN	=rows[3][0+0*NCH*NROWS*NVAL],
					NNNE	=rows[3][0+1*NCH*NROWS*NVAL],
					NNNEE	=rows[3][0+2*NCH*NROWS*NVAL],
					NNNEEE	=rows[3][0+3*NCH*NROWS*NVAL],
					NNWWWW	=rows[2][0-4*NCH*NROWS*NVAL],
					NNWWW	=rows[2][0-3*NCH*NROWS*NVAL],
					NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
					NNW	=rows[2][0-1*NCH*NROWS*NVAL],
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL],
					NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
					NNEEE	=rows[2][0+3*NCH*NROWS*NVAL],
					NNEEEE	=rows[2][0+4*NCH*NROWS*NVAL],
					NNEEEEE	=rows[2][0+5*NCH*NROWS*NVAL],
					NWWWWWW	=rows[1][0-6*NCH*NROWS*NVAL],
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

					dNNNNNNN=rows[7][2+0*NCH*NROWS*NVAL],
					dNNNNNN	=rows[6][2+0*NCH*NROWS*NVAL],
					dNNNNN	=rows[5][2+0*NCH*NROWS*NVAL],
					dNNNNWW	=rows[4][2-2*NCH*NROWS*NVAL],
					dNNNN	=rows[4][2+0*NCH*NROWS*NVAL],
					dNNNNEE	=rows[4][2+2*NCH*NROWS*NVAL],
					dNNNWWW	=rows[3][2-3*NCH*NROWS*NVAL],
					dNNNWW	=rows[3][2-2*NCH*NROWS*NVAL],
					dNNNW	=rows[3][2-1*NCH*NROWS*NVAL],
					dNNN	=rows[3][2+0*NCH*NROWS*NVAL],
					dNNNE	=rows[3][2+1*NCH*NROWS*NVAL],
					dNNNEE	=rows[3][2+2*NCH*NROWS*NVAL],
					dNNNEEE	=rows[3][2+3*NCH*NROWS*NVAL],
					dNNWWW	=rows[2][2-3*NCH*NROWS*NVAL],
					dNNWW	=rows[2][2-2*NCH*NROWS*NVAL],
					dNNW	=rows[2][2-1*NCH*NROWS*NVAL],
					dNN	=rows[2][2+0*NCH*NROWS*NVAL],
					dNNE	=rows[2][2+1*NCH*NROWS*NVAL],
					dNNEE	=rows[2][2+2*NCH*NROWS*NVAL],
					dNNEEE	=rows[2][2+3*NCH*NROWS*NVAL],
					dNNEEEE	=rows[2][2+4*NCH*NROWS*NVAL],
					dNNEEEEE=rows[2][2+5*NCH*NROWS*NVAL],
					dNWWWWWW=rows[1][2-6*NCH*NROWS*NVAL],
					dNWWWWW	=rows[1][2-5*NCH*NROWS*NVAL],
					dNWWWW	=rows[1][2-4*NCH*NROWS*NVAL],
					dNWWW	=rows[1][2-3*NCH*NROWS*NVAL],
					dNWW	=rows[1][2-2*NCH*NROWS*NVAL],
					dNW	=rows[1][2-1*NCH*NROWS*NVAL],
					dN	=rows[1][2+0*NCH*NROWS*NVAL],
					dNE	=rows[1][2+1*NCH*NROWS*NVAL],
					dNEE	=rows[1][2+2*NCH*NROWS*NVAL],
					dNEEE	=rows[1][2+3*NCH*NROWS*NVAL],
					dNEEEE	=rows[1][2+4*NCH*NROWS*NVAL],
					dNEEEEE	=rows[1][2+5*NCH*NROWS*NVAL],
					dNEEEEEE=rows[1][2+6*NCH*NROWS*NVAL],
					dWWWWWWW=rows[0][2-7*NCH*NROWS*NVAL],
					dWWWWWW	=rows[0][2-6*NCH*NROWS*NVAL],
					dWWWWW	=rows[0][2-5*NCH*NROWS*NVAL],
					dWWWW	=rows[0][2-4*NCH*NROWS*NVAL],
					dWWW	=rows[0][2-3*NCH*NROWS*NVAL],
					dWW	=rows[0][2-2*NCH*NROWS*NVAL],
					dW	=rows[0][2-1*NCH*NROWS*NVAL],

					xN	=rows[1][3+0*NCH*NROWS*NVAL],
					xNE	=rows[1][3+1*NCH*NROWS*NVAL],
					xNEE	=rows[1][3+2*NCH*NROWS*NVAL],
					xNEEE	=rows[1][3+3*NCH*NROWS*NVAL],
					xW	=rows[0][3-1*NCH*NROWS*NVAL];
				ALIGN(32) int32_t estim[L1NPREDS];
				int64_t pred;
				int32_t hpred1, hpred2, hpredf;
				int32_t tidx=0, bit=0;
				int32_t error;
				int32_t nbypass, nbypass0;
				int32_t nzeros=0, grflag;
				Cell_t *statsptr;
				int32_t ctx;
				int epred;
				int j;
				
#if 1
				(void)NNNNNNN	;
				(void)NNNNNNNE	;
				(void)NNNNNN	;
				(void)NNNNNNE	;
				(void)NNNNN	;
				(void)NNNNNE	;
				(void)NNNNWW	;
				(void)NNNN	;
				(void)NNNNE	;
				(void)NNNNEE	;
				(void)NNNWWW	;
				(void)NNNWW	;
				(void)NNNW	;
				(void)NNN	;
				(void)NNNE	;
				(void)NNNEE	;
				(void)NNNEEE	;
				(void)NNWWWW	;
				(void)NNWWW	;
				(void)NNWW	;
				(void)NNW	;
				(void)NN	;
				(void)NNE	;
				(void)NNEE	;
				(void)NNEEE	;
				(void)NNEEEE	;
				(void)NNEEEEE	;
				(void)NWWWWWW	;
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

				(void)dNNNNNNN	;
				(void)dNNNNNN	;
				(void)dNNNNN	;
				(void)dNNNNWW	;
				(void)dNNNN	;
				(void)dNNNNEE	;
				(void)dNNNWWW	;
				(void)dNNNWW	;
				(void)dNNNW	;
				(void)dNNN	;
				(void)dNNNE	;
				(void)dNNNEE	;
				(void)dNNNEEE	;
				(void)dNNWWW	;
				(void)dNNWW	;
				(void)dNNW	;
				(void)dNN	;
				(void)dNNE	;
				(void)dNNEE	;
				(void)dNNEEE	;
				(void)dNNEEEE	;
				(void)dNNEEEEE	;
				(void)dNWWWWWW	;
				(void)dNWWWWW	;
				(void)dNWWWW	;
				(void)dNWWW	;
				(void)dNWW	;
				(void)dNW	;
				(void)dN	;
				(void)dNE	;
				(void)dNEE	;
				(void)dNEEE	;
				(void)dNEEEE	;
				(void)dNEEEEE	;
				(void)dNEEEEEE	;
				(void)dWWWWWWW	;
				(void)dWWWWWW	;
				(void)dWWWWW	;
				(void)dWWWW	;
				(void)dWWW	;
				(void)dWW	;
				(void)dW	;

				(void)xN	;
				(void)xNE	;
				(void)xNEE	;
				(void)xNEEE	;
				(void)xW	;
#endif
				if(kc==0)offset0=offset_bias[0];
				if(kc==1)offset0=cu0*yuv[0]+offset_bias[1];
				if(kc==2)offset0=cv0*yuv[0]+cv1*yuv[1]+offset_bias[2];
				offset=offset0>>RCT_BITS;
				pred=(int64_t)bias[kc];
#define PRED(EXPR) estim[j]=EXPR; pred+=coeffs[kc][j]*estim[j]; ++j;
				j=0;
				PREDLIST;
#undef  PRED
				hpred2=hpred1=(int32_t)(pred>>(sh-ADDBITS));
				ctx=31^LZCNT32(xW*xW+1);
				nbypass=(ctx>>1)-GRBITS;
				CLAMP2(nbypass, 0, 7);
				{
					int nzctx=ctx;
					if(nzctx>1)
						nzctx=1;
					CLAMP2(ctx, CTXMIN, CTXMAX);
					ctx-=CTXMIN-nzctx;
				}
			//	if(ctx>NCTX-1)ctx=NCTX-1;
				hpredf=hpred2+(offset0<<ADDBITS>>RCT_BITS);
				CLAMP2(hpredf, 0, 255<<ADDBITS);
			//	CLAMP2(hpredf, lumamin[kc]<<ADDBITS, lumamax[kc]<<ADDBITS);
				pred=hpredf>>ADDBITS;
				epred=128-abs((int32_t)pred-128);
				error=0;
				if(fwd)
				{
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
					nzeros=error>>nbypass;
				}
				statsptr=stats1[kc][ctx][hpredf];
				tidx=0;
				_mm_prefetch((char*)statsptr, _MM_HINT_T0);
				do
				{
					bit=tidx>=nzeros;
					//if(fwd&&!kc&&ctx==NCTX/2&&hpredf==60)
					//	printf("%d", bit);
					codebit(ac, statsptr+tidx, &bit, fwd);
					if(bit)
						break;
					++tidx;
				}while(tidx<GRLIMIT);
				nbypass0=nbypass;
				grflag=tidx==GRLIMIT;
				statsptr=stats2[kc][hpredf];
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
				tidx+=256>>nbypass;//bit coding:  tidx=2*tidx+bit  tidx=0b1XX
				{
					int32_t kb=nbypass-1;
					
					for(;kb>=0;--kb)
					{
						bit=error>>kb&1;
						codebit(ac, statsptr+tidx, &bit, fwd);
						tidx=2*tidx+bit;
					}
				}
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
					int negmask=((int32_t)pred-128)>>31;
					int sym=error;
					int e2=epred-sym;
					error=sym>>1^-(sym&1);
					e2=(e2^negmask)-negmask;
					if(2*epred<sym&&2*pred+sym!=512)
						error=e2;
					yuv[kc]=(uint8_t)(error+(int32_t)pred);
#ifdef ENABLE_GUIDE
					uint8_t *pval=&g_image[3*(iw*ky+kx)+rct->perm[kc]];
					uint8_t val=*pval;
					uint8_t pixel=yuv[kc];
					if(pixel!=val)
						CRASH("GUIDE YXC %d %d %d", ky, kx, kc);
#endif
				}
				{
					int32_t curr=yuv[kc]-offset, e2;
					
					rows[0][0]=curr;
					rows[0][1]=curr-(hpred1>>ADDBITS);

					error=(yuv[kc]<<ADDBITS)-hpredf;
					rows[0][2]=error;

					e2=error=yuv[kc]-(int32_t)pred;
					error=error<<1^error>>31;
					error<<=GRBITS;
					rows[0][3]=(xW+(xW<xNE?xW:xNE)+error+(xNEE>xNEEE?xNEE:xNEEE))>>2;
					
					hpred1>>=ADDBITS;
					error=(curr>hpred1)-(curr<hpred1);
					bias[kc]+=(error<<6)+e2;
#define PRED(EXPR) coeffs[kc][j]+=error*estim[j]; ++j;
					j=0;
					PREDLIST;
#undef  PRED
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
}
int c12_codec(int argc, char **argv)
{
	static const uint16_t tag='1'|'2'<<8;

	const char *srcfn=0, *dstfn=0;
	ptrdiff_t srcsize=0, dstsize=0;
	int fwd=0;
	int32_t iw=0, ih=0;
	ptrdiff_t res=0, usize=0, csize=0;
	uint8_t *image=0, *imptr=0, *imend=0, *stream=0, *streamptr=0, *streamend=0;
	ptrdiff_t streamsize=0;
	RCTInfo rctinfo={0};
	ACState ac={0};
	uint16_t lumamean[3]={0};
	uint8_t lumamin[3]={0}, lumamax[3]={0};
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
		printf("Usage:  \"%s\"  input  output\n", argv[0]);
		return 1;
	}
	srcfn=argv[1];
	dstfn=argv[2];
	ac.hi=0xFFFFFFFFFFFF;
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
			fread(&iw, 1, 3, fsrc);
			fread(&ih, 1, 3, fsrc);
			fread(&rctinfo.perm[0], 1, 1, fsrc);
			fread(&rctinfo.perm[1], 1, 1, fsrc);
			fread(&rctinfo.perm[2], 1, 1, fsrc);
			fread(&rctinfo.cu0, 1, 2, fsrc);
			fread(&rctinfo.cv0, 1, 2, fsrc);
			fread(&rctinfo.cv1, 1, 2, fsrc);
			fread(&lumamean, 1, sizeof(lumamean), fsrc);
			fread(&lumamin, 1, sizeof(lumamin), fsrc);
			fread(&lumamax, 1, sizeof(lumamax), fsrc);
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
		int64_t counters[OCH_COUNT]={0};
		int rowstride=3*iw;
		int64_t mean[3]={0};

		lumamin[0]=255;
		lumamin[1]=255;
		lumamin[2]=255;
		for(imptr=image;imptr<imend;imptr+=3)
		{
			int r=imptr[0];
			int g=imptr[1];
			int b=imptr[2];
			mean[0]+=r;
			mean[1]+=g;
			mean[2]+=b;
			if(lumamin[0]>r)lumamin[0]=r;
			if(lumamin[1]>g)lumamin[1]=g;
			if(lumamin[2]>b)lumamin[2]=b;
			if(lumamax[0]<r)lumamax[0]=r;
			if(lumamax[1]<g)lumamax[1]=g;
			if(lumamax[2]<b)lumamax[2]=b;
		}
		lumamean[0]=(uint16_t)(((mean[0]<<RCT_BITS)+(res>>1))/res);
		lumamean[1]=(uint16_t)(((mean[1]<<RCT_BITS)+(res>>1))/res);
		lumamean[2]=(uint16_t)(((mean[2]<<RCT_BITS)+(res>>1))/res);
		imptr=image+rowstride;
		{
			ALIGN(16) int16_t rramp[8]={0}, gramp[8]={0}, bramp[8]={0};
			int prev[3]={0};
			while(imptr<imend)
			{
				int r, g, b;
				int r2, g2, b2;

				r=imptr[0]-imptr[0-rowstride];
				g=imptr[1]-imptr[1-rowstride];
				b=imptr[2]-imptr[2-rowstride];
				imptr+=3;
				r2=r;
				g2=g;
				b2=b;
				r-=prev[0];
				g-=prev[1];
				b-=prev[2];
				prev[0]=r2;
				prev[1]=g2;
				prev[2]=b2;

				r<<=RCT_BITS;
				g<<=RCT_BITS;
				b<<=RCT_BITS;

				rramp[1-1]=r*1>>RCT_BITS;
				gramp[1-1]=g*1>>RCT_BITS;
				bramp[1-1]=b*1>>RCT_BITS;
				rramp[2-1]=r*2>>RCT_BITS;
				gramp[2-1]=g*2>>RCT_BITS;
				bramp[2-1]=b*2>>RCT_BITS;
				rramp[3-1]=r*3>>RCT_BITS;
				gramp[3-1]=g*3>>RCT_BITS;
				bramp[3-1]=b*3>>RCT_BITS;
				rramp[4-1]=r*4>>RCT_BITS;
				gramp[4-1]=g*4>>RCT_BITS;
				bramp[4-1]=b*4>>RCT_BITS;
				rramp[5-1]=r*5>>RCT_BITS;
				gramp[5-1]=g*5>>RCT_BITS;
				bramp[5-1]=b*5>>RCT_BITS;
				rramp[6-1]=r*6>>RCT_BITS;
				gramp[6-1]=g*6>>RCT_BITS;
				bramp[6-1]=b*6>>RCT_BITS;
				rramp[7-1]=r*7>>RCT_BITS;
				gramp[7-1]=g*7>>RCT_BITS;
				bramp[7-1]=b*7>>RCT_BITS;
				rramp[8-1]=r*8>>RCT_BITS;
				gramp[8-1]=g*8>>RCT_BITS;
				bramp[8-1]=b*8>>RCT_BITS;

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
					counters[k]+=abs(och[k]);
			}
		}
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
						crct_get(&rctinfo2, c0, c1, c2);
						print_rct2(&rctinfo2);
						printf(" %d %3d %3d _%s_%s_%s  "
							, kp, k1, k2
							, och_names[c0]
							, och_names[c1]
							, och_names[c2]
						);
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
			(void)it_count;
			(void)it_select;
			(void)&print_rct2;
			(void)&print_rct;
#ifdef LOUD
			printf("WH %d*%d  %lld bytes  RCT %4d/%4d ", iw, ih, usize, it_select, it_count);
			print_rct2(&rctinfo);
			printf("  \"%s\"\n", srcfn);
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
	ac.ptr=ac.start=streamptr;
	ac.end=streamend;
	if(fwd)	mainloop(iw, ih, &rctinfo, image, stream, &ac, 1, lumamean, lumamin, lumamax);
	else	mainloop(iw, ih, &rctinfo, image, stream, &ac, 0, lumamean, lumamin, lumamax);
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
			dstsize+=fwrite(&lumamean, 1, sizeof(lumamean), fdst);
			dstsize+=fwrite(&lumamin, 1, sizeof(lumamin), fdst);
			dstsize+=fwrite(&lumamax, 1, sizeof(lumamax), fdst);
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
#ifdef _MSC_VER
	if(fwd)
	{
		extern int64_t g_csize;
		g_csize=csize;
	}
#endif
	(void)dstsize;
	(void)csize;
	(void)&time_sec2;
	return 0;
}
