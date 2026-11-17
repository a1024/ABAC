#ifdef _MSC_VER
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif
#elif defined __linux__ && !defined _GNU_SOURCE
#	define _GNU_SOURCE
#endif
#include<stddef.h>
#include<stdint.h>
#include<inttypes.h>
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<sys/stat.h>
#if defined _WIN32 || defined WIN32
#	define WIN32_LEAN_AND_MEAN
#	include<Windows.h>
#else
#	include<time.h>
#endif
#if defined _M_X64 || defined __x86_64__
#	ifdef _MSC_VER
#		include<intrin.h>
#	elif defined __GNUC__
#		include<x86intrin.h>
#	endif
#endif
#if defined _MSC_VER && defined _WIN32
#	include<process.h>
#	define THREAD_CALL __stdcall
typedef unsigned THREAD_RET;
#else
#	include<pthread.h>
#	define THREAD_CALL
typedef void *THREAD_RET;
#endif


#ifdef _MSC_VER
	#define LOUD
	#define ENABLE_GUIDE
#ifdef _DEBUG
	#define MT_DISABLE
//	#define PROFILE_SIZE
//	#define FIFOVAL
#endif
//	#define PRINT_SQUASH
#endif


#define SHIFT_EXPR 7+(ky>64)
#if 1
//LPCB 1013463642 B       43.189465      46.577634 sec     76.457671    70.895956 MB/s	ETAP7
//LPCB 1015939746 B       37.780346      43.329843 sec     87.404331    76.209967 MB/s	ETAP5

//weakest to strongest contributor
#define PREDLIST\
	PRED(N, N, N)\
	PRED(W, W, W)\
	PRED(N+NE-NNE, N+NE-NNE+((m1-(N1+NE1-NNE1))*rct.uc0>>RCTBITS), N+NE-NNE+(((m1-(N1+NE1-NNE1))*rct.vc0+(m2-(N2+NE2-NNE2))*rct.vc1)>>RCTBITS))\
	PRED(W+NE-N, W+NE-N+((m1-(W1+NE1-N1))*rct.uc0>>RCTBITS), W+NE-N+(((m1-(W1+NE1-N1))*rct.vc0+(m2-(W2+NE2-N2))*rct.vc1)>>RCTBITS))\
	PRED((WWWW+WWW+NEE+NEEE)>>2, m1, m1)\
	PRED((N+NN)>>1, N+m1-N1, N+(m1-N1+m2-N2)/2)\
	PRED((W+WW)>>1, W+m1-W1, W+(m1-W1+m2-W2)/2)\
	PRED(N+W-NW, N+W-NW+((m1-(N1+W1-NW1))*rct.uc0>>RCTBITS), N+W-NW+(((m1-(N1+W1-NW1))*rct.vc0+(m2-(N2+W2-NW2))*rct.vc1)>>RCTBITS))\
	PRED(2*N-NN, 2*N-NN+((m1-(2*N1-NN1))*rct.uc0>>RCTBITS), 2*N-NN+(((m1-(2*N1-NN1))*rct.vc0+(m2-(2*N2-NN2))*rct.vc1)>>RCTBITS))\
	PRED(2*W-WW, 2*W-WW+((m1-(2*W1-WW1))*rct.uc0>>RCTBITS), 2*W-WW+(((m1-(2*W1-WW1))*rct.vc0+(m2-(2*W2-WW2))*rct.vc1)>>RCTBITS))\
	PRED(3*(N-NN)+NNN, 3*(N-NN)+NNN+((m1-(3*(N1-NN1)+NNN1))*rct.uc0>>RCTBITS), 3*(N-NN)+NNN+(((m1-(3*(N1-NN1)+NNN1))*rct.vc0+(m2-(3*(N2-NN2)+NNN2))*rct.vc1)>>RCTBITS))\
	PRED(3*(W-WW)+WWW, 3*(W-WW)+WWW+((m1-(3*(W1-WW1)+WWW1))*rct.uc0>>RCTBITS), 3*(W-WW)+WWW+(((m1-(3*(W1-WW1)+WWW1))*rct.vc0+(m2-(3*(W2-WW2)+WWW2))*rct.vc1)>>RCTBITS))\

#endif
#if 0
//LPCB 1017638958 B       38.603066      43.045737 sec     85.541544    76.712960 MB/s	ETAP7
//LPCB 1019837346 B       34.305270      40.837117 sec     96.258267    80.861875 MB/s	ETAP5
//LPCB 1024692998 B       32.338126      36.990283 sec    102.113706    89.271170 MB/s	ETAP3
#define PREDLIST\
	PRED(N, N, N)\
	PRED(W, W, W)\
	PRED(NN, N+m1-N1, N+(m1-N1+m2-N2)/2)\
	PRED(WW, W+m1-W1, W+(m1-W1+m2-W2)/2)\
	PRED(W+NE-N, W+NE-N+m1-(W1+NE1-N1), W+NE-N+(m1-(W1+NE1-N1)+m2-(W2+NE2-N2))/2)\
	PRED((WWWW+WWW+NEE+NEEE)>>2, m1, m1)\
	PRED(N+W-NW, N+W-NW+((m1-(N1+W1-NW1))*rct.uc0>>RCTBITS), N+W-NW+(((m1-(N1+W1-NW1))*rct.vc0+(m2-(N2+W2-NW2))*rct.vc1)>>RCTBITS))\
	PRED((W+2*NE-NNE)/2, (W+2*NE-NNE)/2+((m1-(W1+2*NE1-NNE1)/2)*rct.uc0>>RCTBITS), (W+2*NE-NNE)/2+(((m1-(W1+2*NE1-NNE1)/2)*rct.vc0+(m2-(W2+2*NE2-NNE2)/2)*rct.vc1)>>RCTBITS))\
	PRED(2*N-NN, 2*N-NN+((m1-(2*N1-NN1))*rct.uc0>>RCTBITS), 2*N-NN+(((m1-(2*N1-NN1))*rct.vc0+(m2-(2*N2-NN2))*rct.vc1)>>RCTBITS))\
	PRED(2*W-WW, 2*W-WW+((m1-(2*W1-WW1))*rct.uc0>>RCTBITS), 2*W-WW+(((m1-(2*W1-WW1))*rct.vc0+(m2-(2*W2-WW2))*rct.vc1)>>RCTBITS))\

#endif
#if 0//LPCB 1075810202 B       21.874057      26.616999 sec    150.962663   124.062291 MB/s	NPREDS 5
#define PREDLIST\
	PRED(W, W, W)\
	PRED(NE, NE+((m1-NE1)*rct.uc0>>RCTBITS), NE+(((m1-NE1)*rct.vc0+(m2-NE2)*rct.vc1)>>RCTBITS))\
	PRED(N+W-NW, N+W-NW+((m1-(N1+W1-NW1))*rct.uc0>>RCTBITS), N+W-NW+(((m1-(N1+W1-NW1))*rct.vc0+(m2-(N2+W2-NW2))*rct.vc1)>>RCTBITS))\
	PRED(2*N-NN, 2*N-NN+((m1-(2*N1-NN1))*rct.uc0>>RCTBITS), 2*N-NN+(((m1-(2*N1-NN1))*rct.vc0+(m2-(2*N2-NN2))*rct.vc1)>>RCTBITS))\
	PRED((N+W)>>1, m1, m1)\

#endif
#if 0//LPCB 1086245706 B       22.882644      28.024728 sec    144.308756   117.830433 MB/s
#define PREDLIST\
	PRED(W, W, W)\
	PRED(NE, NE+((m1-NE1)*rct.uc0>>RCTBITS), NE+(((m1-NE1)*rct.vc0+(m2-NE2)*rct.vc1)>>RCTBITS))\
	PRED(N+W-NW, N+W-NW+((m1-(N1+W1-NW1))*rct.uc0>>RCTBITS), N+W-NW+(((m1-(N1+W1-NW1))*rct.vc0+(m2-(N2+W2-NW2))*rct.vc1)>>RCTBITS))\
	PRED(2*N-NN, 2*N-NN+((m1-(2*N1-NN1))*rct.uc0>>RCTBITS), 2*N-NN+(((m1-(2*N1-NN1))*rct.vc0+(m2-(2*N2-NN2))*rct.vc1)>>RCTBITS))\

#endif


#ifdef _MSC_VER
#	define ALIGN(N) __declspec(align(N))
#	define INLINE __forceinline static
#	if _MSC_VER<1900
#		define snprintf sprintf_s
#	endif
#else
#	define ALIGN(N) __attribute__((aligned(N)))
#	define INLINE __attribute__((always_inline)) inline static
#	ifndef _countof
#		define _countof(A) (sizeof(A)/sizeof(*(A)))
#	endif
#endif
#define CLAMP(X, L, H) X=X<(L)?L:X, X=X>(H)?H:X
#if _MSC_VER
#define LZCNT32 _lzcnt_u32
#define LZCNT64 _lzcnt_u64
#define TZCNT32 _tzcnt_u32
#define TZCNT64 _tzcnt_u64
#else
#define LZCNT32(X) (X?__builtin_clz(X):32)
#define LZCNT64(X) (X?__builtin_clzll(X):64)
#define TZCNT32(X) (X?__builtin_ctz(X):32)
#define TZCNT64(X) (X?__builtin_ctzll(X):64)
#endif
#define ROUND32(X) _mm_cvt_ss2si(_mm_set_ss(X))
#define ROUND64(X) _mm_cvtsd_si64(_mm_set_sd(X))
#define TRUNC32(X) _mm_cvtt_ss2si(_mm_set_ss(X))
#define TRUNC64(X) _mm_cvttsd_si64(_mm_set_sd(X))
//#define ROUND32(X) roundf(X)
//#define ROUND64(X) round(X)
//#define TRUNC32(X) truncf(X)
//#define TRUNC64(X) trunc(X)
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
static double time_sec(void)//time since first call
{
#if defined _WIN32 || defined WIN32
	static long long t0=0;
	LARGE_INTEGER li;
	double t;
	QueryPerformanceCounter(&li);//computer up time
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
	void *p=0;
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	p=realloc(g_image, size);
	if(!p)
	{
		CRASH("Alloc error");
		return;
	}
	g_image=(uint8_t*)p;
	memcpy(g_image, image, size);
}
static void guide_check(uint8_t *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		CRASH("guide  XY %d %d", kx, ky);
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
static int query_cpu_cores(void)
{
#ifdef _WIN32
	SYSTEM_INFO info;
	GetNativeSystemInfo(&info);
	return info.dwNumberOfProcessors;
#else
	return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}
#if 0
static int clamp(int x, int lo, int hi)
{
	CLAMP(x, lo, hi);
	return x;
}
static int clip(int x, int a, int b, int c)
{
	int vmin, vmax;

	vmax=a; vmin=b;
	if(a<b)vmin=a, vmax=b;
	if(vmax<c)vmax=c;
	if(vmin>c)vmin=c;
	CLAMP(x, vmin, vmax);
	return x;
}
#endif

//RCT
static const uint8_t perms[]=
{
	1, 0, 2,	//good
	1, 2, 0,

	//0, 1, 2,	2, 1, 0,	//bad
	//2, 0, 1,	1, 0, 2,
	//1, 2, 0,	0, 2, 1,
};
enum
{
	ETAP=5,			// <- ETAP		3 good, 1,5 bad
	PLEVELS=256,		//POT			256 best
	ELEVELS=PLEVELS-1+ETAP,
	PROBST1=7,	MIXST1=0,
	PROBST3=29,	MIXST3=17,
	PROBST5=29,	MIXST5=17,
	PROBST7=29,	MIXST7=17,
	PROBSTB=29,	MIXSTB=17,
	TABLEBITS=10,
	PROBBITS=13,
	PROBBITS0=12,

	BLOCKDIM=1024,//0xFFFF
	NCTX=12,
	RLIMIT=11,
	PBITS=5,
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED

	RCTBITS=3,
	RCTMAX=1<<RCTBITS,
	RCTLEVELS=RCTMAX+1,//+unity
	NRCTS=RCTLEVELS*RCTLEVELS*(RCTLEVELS+1)/2,
	NPERMS=_countof(perms)/3,

	XPAD=8,
	NCH=3,
	NROWS=4,
	NVAL=3+NPREDS,
};
typedef struct _RCTInfo
{
	uint8_t pidx, uc0, vc0, vc1;
} RCTInfo;
typedef struct _AnalysisCtrs
{
	uint32_t yctr, uctrs[RCTLEVELS], vctrs[RCTLEVELS*(RCTLEVELS+1)/2];
} AnalysisCtrs;
typedef struct _CodecCtx
{
	uint8_t *image, *streams;
	int *csizes;
	size_t *cbounds;//nblocks+1
	int iw, ih, blockdim, fwd, nthreads, ncores;//base
	int xblocks, yblocks, nblocks;//derived
	size_t usize, ccap, csize;
	int blockcap;
} CodecCtx;
typedef int16_t Pixel_t;
typedef int16_t Cell_t;
typedef struct _ThreadCtx
{
	CodecCtx *cc;
	int tid, b1, b2;
	uint8_t *block;
	Cell_t stats[NPREDS*256*ELEVELS];
} ThreadCtx;
static void print_rct(RCTInfo *rct)
{
	printf("RCT%c%c%c_%c_%c%c/%c"
		, '0'+perms[rct->pidx*3+0]
		, '0'+perms[rct->pidx*3+1]
		, '0'+perms[rct->pidx*3+2]
		, rct->uc0+(rct->uc0<9?'0':'A'-10)
		, rct->vc0+(rct->vc0<9?'0':'A'-10)
		, rct->vc1+(rct->vc1<9?'0':'A'-10)
		, RCTMAX+(RCTMAX<9?'0':'A'-10)
	);
}
typedef struct _ACState
{
	uint64_t lo, hi, code;
	uint8_t *ptr, *end, *start;

	int64_t byteidx;
} ACState;
#ifdef PROFILE_SIZE
static int g_kc;
static double esizes[24];
static double shannontable[1<<PROBBITS];
#endif
static uint16_t squashtable[1<<TABLEBITS]={0};
#if defined _MSC_VER &&0
static int g_debug=0;
#endif
static int codesym(ACState *ac, Cell_t **stats, int64_t *mixer, int shift, int sym, int fwd)
{
	const int learngain[]=
	{//kb	0	1	2	3	4	5	6	7
		6,	8,	12,	16,	21,	25,	27,	32,
	};
	uint64_t x;
	int64_t p0;
	int32_t probs[64];
	int16_t pupdate[64];
	int kb, tidx, bit, proberror;
	
	if(!fwd)
		sym=0;
	tidx=ELEVELS;
	for(kb=7;kb>=0;--kb)
	{
		switch((int)ETAP)
		{
		case 1:
			p0=1LL<<PROBST1>>1;
			for(int k=0;k<NPREDS;++k)
			{
				p0+=stats[k][tidx+0];
				pupdate[k]=stats[k][tidx+ETAP/2];
			}
			p0>>=PROBST1;// *1
			break;
		case 3:
			probs[0]=0;
			probs[1]=0;
			probs[2]=0;
			for(int k=0;k<NPREDS;++k)
			{
				probs[0]+=stats[k][tidx+0];
				probs[1]+=stats[k][tidx+1];
				probs[2]+=stats[k][tidx+2];
				pupdate[k]=stats[k][tidx+ETAP/2];
			}
			p0=1LL<<PROBST3>>1;
			for(int k=0;k<ETAP;++k)
				p0+=mixer[k]*probs[k];
			p0>>=PROBST3;// *4
			break;
		case 5:
			probs[0]=0;
			probs[1]=0;
			probs[2]=0;
			probs[3]=0;
			probs[4]=0;
			for(int k=0;k<NPREDS;++k)
			{
				probs[0]+=stats[k][tidx+0];
				probs[1]+=stats[k][tidx+1];
				probs[2]+=stats[k][tidx+2];
				probs[3]+=stats[k][tidx+3];
				probs[4]+=stats[k][tidx+4];
				pupdate[k]=stats[k][tidx+ETAP/2];
			}
			p0=1LL<<PROBST5>>1;
			for(int k=0;k<ETAP;++k)
				p0+=mixer[k]*probs[k];
			p0>>=PROBST5;// *16
			break;
		case 7:
			probs[0]=0;
			probs[1]=0;
			probs[2]=0;
			probs[3]=0;
			probs[4]=0;
			probs[5]=0;
			probs[6]=0;
			for(int k=0;k<NPREDS;++k)
			{
				probs[0]+=stats[k][tidx+0];
				probs[1]+=stats[k][tidx+1];
				probs[2]+=stats[k][tidx+2];
				probs[3]+=stats[k][tidx+3];
				probs[4]+=stats[k][tidx+4];
				probs[5]+=stats[k][tidx+5];
				probs[6]+=stats[k][tidx+6];
				pupdate[k]=stats[k][tidx+ETAP/2];
			}
			p0=1LL<<PROBST7>>1;
			for(int k=0;k<ETAP;++k)
				p0+=mixer[k]*probs[k];
			p0>>=PROBST7;// *16
			break;
		case 11:
			probs[0x0]=0;
			probs[0x1]=0;
			probs[0x2]=0;
			probs[0x3]=0;
			probs[0x4]=0;
			probs[0x5]=0;
			probs[0x6]=0;
			probs[0x7]=0;
			probs[0x8]=0;
			probs[0x9]=0;
			probs[0xA]=0;
			for(int k=0;k<NPREDS;++k)
			{
				probs[0x0]+=stats[k][tidx+0x0];
				probs[0x1]+=stats[k][tidx+0x1];
				probs[0x2]+=stats[k][tidx+0x2];
				probs[0x3]+=stats[k][tidx+0x3];
				probs[0x4]+=stats[k][tidx+0x4];
				probs[0x5]+=stats[k][tidx+0x5];
				probs[0x6]+=stats[k][tidx+0x6];
				probs[0x7]+=stats[k][tidx+0x7];
				probs[0x8]+=stats[k][tidx+0x8];
				probs[0x9]+=stats[k][tidx+0x9];
				probs[0xA]+=stats[k][tidx+0xA];
				pupdate[k]=stats[k][tidx+ETAP/2];
			}
			p0=1LL<<PROBSTB>>1;
			for(int k=0;k<ETAP;++k)
				p0+=mixer[k]*probs[k];
			p0>>=PROBSTB;// *16
			break;
		}
		CLAMP(p0, -((int)_countof(squashtable)/2-1), (int)_countof(squashtable)/2-1);
		p0=(squashtable+_countof(squashtable)/2)[p0];
		x=ac->hi-ac->lo;
		if(x<=0xFFFF)
		{
#if defined _MSC_VER &&0
			uint64_t c0=ac->code;
#endif
			if(ac->ptr>=ac->end)
				CRASH("inflation");
			if(fwd)
				*(uint32_t*)ac->ptr=(uint32_t)(ac->lo>>32);
			else
				ac->code=ac->code<<32|*(uint32_t*)ac->ptr;
#if defined _MSC_VER &&0
			if(!*(uint32_t*)ac->ptr)//
				printf(" %c0 %016llX %016llX %016llX %016llX %d\n"
					, 'D'+fwd
					, ac->lo
					, ac->hi
					, c0
					, x
					, (int)(ac->end-ac->ptr)
				);
#endif
			ac->ptr+=sizeof(uint32_t);
			ac->lo<<=32;
			ac->hi=ac->hi<<32|(uint32_t)~0;
			if(ac->hi<ac->lo)ac->hi=~(uint64_t)0;
			x=ac->hi-ac->lo;
		}
		x=ac->lo+(x*p0>>PROBBITS);
		bit=fwd?sym>>kb&1:ac->code>=x;
#if defined FIFOVAL &&1
		if(fwd)
			fifoval_enqueue(bit<<31^(int)p0);
		else
			fifoval_check(bit<<31^(int)p0);
#endif
		*(bit?&ac->lo:&ac->hi)=x-(bit^1);
#ifdef PROFILE_SIZE
		if(fwd)
			esizes[8*g_kc+kb]+=shannontable[bit?(1LL<<PROBBITS)-p0:p0];
#endif
		proberror=((bit^1)<<PROBBITS)-(int)p0;
#if defined _MSC_VER &&0
		if(g_debug)//
			printf("[%d] %5d %d  %016llX  %016llX  %016llX\n"
				, kb
				, (int)p0
				, bit
				, ac->lo
				, ac->hi
				, ac->code
			);//
#endif
		if(ETAP>1)
		{
			int sh=0;
			switch((int)ETAP)
			{
			case  1:sh=MIXST1;break;
			case  3:sh=MIXST3;break;
			case  5:sh=MIXST5;break;
			case  7:sh=MIXST7;break;
		//	case  9:sh=MIXST9;break;
			case 11:sh=MIXSTB;break;
		//	case 13:sh=MIXSTD;break;
		//	case 15:sh=MIXSTF;break;
			}
			for(int k=0;k<ETAP;++k)
				mixer[k]+=((int64_t)proberror*probs[k]+(1LL<<sh>>1))>>sh;
		}
		proberror=(proberror*learngain[kb]+(1<<shift>>1))>>shift;
		if(proberror)
		{
			__m128i u=_mm_set1_epi16(proberror);
			for(int k=0;k<(NPREDS+7)>>3;++k)
			{
				__m128i v=_mm_load_si128((__m128i*)pupdate+k);
				v=_mm_adds_epi16(v, u);
				_mm_store_si128((__m128i*)pupdate+k, v);
			}
			for(int k=0;k<NPREDS;++k)
				stats[k][tidx+ETAP/2]=pupdate[k];
			//for(int k=0;k<NPREDS;++k)
			//{
			//	int32_t cell=stats[k][tidx+ETAP/2];
			//	cell+=proberror;
			//	CLAMP(cell, -0x8000, 0x7FFF);
			//	stats[k][tidx+ETAP/2]=(Cell_t)cell;
			//}
		}
		tidx=2*tidx+ELEVELS*bit;
		sym|=bit<<kb;
	}
	return sym;
}
static THREAD_RET THREAD_CALL codec_thread(void *p)
{
	ThreadCtx *a=(ThreadCtx*)p;
	CodecCtx *cc=a->cc;
	uint8_t *block=a->block+3*8*(cc->iw<cc->blockdim?cc->iw:cc->blockdim);
	Cell_t *stats[NPREDS]={0};

	if((int)(block-a->block)+3*(cc->iw<cc->blockdim?cc->iw:cc->blockdim)*(cc->ih<cc->blockdim?cc->ih:cc->blockdim)>cc->blockcap)
		CRASH("small block");
	for(int kb=a->b1;kb<a->b2;++kb)
	{
		int bx, by, x1, x2, y1, y2, dx, dy, rstr;
		RCTInfo rct={0};
		const uint8_t *perm=0;
		ACState ac={0};
		int64_t mixer[ETAP]={0};
#ifdef _DEBUG
		int debug[NPREDS];
#endif
		
		by=kb/cc->xblocks;
		bx=kb%cc->xblocks;

		y1=by*cc->blockdim;
		y2=y1+cc->blockdim;
	//	if(y2+cc->blockdim>cc->ih)y2=cc->ih;
		if(y2>cc->ih)y2=cc->ih;
		dy=y2-y1;

		x1=bx*cc->blockdim;
		x2=x1+cc->blockdim;
	//	if(x2+cc->blockdim>cc->iw)x2=cc->iw;
		if(x2>cc->iw)x2=cc->iw;
		dx=x2-x1;

		rstr=3*dx;
		ac.hi=0xFFFFFFFFFFFF;
		ac.ptr=ac.start=cc->streams+cc->cbounds[kb];
		ac.end=cc->streams+cc->cbounds[kb+1];
		if(cc->fwd)//analysis
		{
			enum
			{
				STRIDE=3,
			};
			AnalysisCtrs counters[NPERMS]={0};
			uint32_t bestscore=0;
			for(int ky=y1<1?1:y1;ky<y2;ky+=STRIDE)
			{
				int kx=x1<1?1:x1;
				uint8_t *imptr=cc->image+3*(cc->iw*ky+kx);
				for(;kx<x2;kx+=STRIDE)
				{
					int rgb[]=
					{
						(imptr[0]-imptr[-rstr+0]-imptr[-3+0]+imptr[-rstr-3+0])<<RCTBITS,
						(imptr[1]-imptr[-rstr+1]-imptr[-3+1]+imptr[-rstr-3+1])<<RCTBITS,
						(imptr[2]-imptr[-rstr+2]-imptr[-3+2]+imptr[-rstr-3+2])<<RCTBITS,
					};
					imptr+=3*STRIDE;
					for(int kp=0;kp<NPERMS;++kp)
					{
						AnalysisCtrs *currctrs=counters+kp;
						int yuv[]=
						{
							rgb[perms[3*kp+0]],
							rgb[perms[3*kp+1]],
							rgb[perms[3*kp+2]],
						};
						currctrs->yctr+=abs(yuv[0]);
						for(int kp1=0, idx=0;kp1<=RCTMAX;++kp1)
						{
							int u=yuv[1]-(kp1*yuv[0]>>RCTBITS);
							currctrs->uctrs[kp1]+=abs(u);
							for(int kp2=0;kp1+kp2<=RCTMAX;++kp2, ++idx)
							{
								int v=yuv[2]-((kp1*yuv[0]+kp2*yuv[1])>>RCTBITS);
								currctrs->vctrs[idx]+=abs(v);
							}
						}
					}
				}
			}
			for(int kp=0;kp<NPERMS;++kp)
			{
				AnalysisCtrs *currctrs=counters+kp;
				for(int uc0=0;uc0<=RCTMAX;++uc0)
				{
					for(int vc0=0, idx=0;vc0<=RCTMAX;++vc0)
					{
						for(int vc1=0;vc0+vc1<=RCTMAX;++vc1, ++idx)
						{
							uint32_t score=currctrs->yctr+currctrs->uctrs[uc0]+currctrs->vctrs[idx];
							if(!bestscore||bestscore>score)
							{
								bestscore=score;
								rct.pidx=kp;
								rct.uc0=uc0;
								rct.vc0=vc0;
								rct.vc1=vc1;
							}
						}
					}
				}
			}
#if defined _MSC_VER && defined MT_DISABLE
			print_rct(&rct);
			printf("\n");
#endif
			memcpy(ac.ptr, &rct, sizeof(rct));
			ac.ptr+=sizeof(rct);
		}
		else
		{
			memcpy(&rct, ac.ptr, sizeof(rct));
			ac.ptr+=sizeof(rct);
			ac.code=*(uint64_t*)ac.ptr;
			ac.code=ac.code<<32|ac.code>>32;
			ac.ptr+=sizeof(uint64_t);
		}
		perm=perms+3*rct.pidx;
		memset(a->stats, 0, sizeof(a->stats));
		memset(a->block, 0, cc->blockcap);
#ifdef PROFILE_SIZE
		g_kc=0;
#endif
		for(int ky=0;ky<dy;++ky)
		{
			int shift=SHIFT_EXPR;
			uint8_t *imptr, *bptr;

			imptr=cc->image+3*(cc->iw*(y1+ky)+x1)+perm[0];
			bptr=block+rstr*ky;
			bptr+=perm[0];
			for(int kx=0;kx<dx;++kx, imptr+=3, bptr+=3)
			{
				uint8_t *Nptr		=bptr	-rstr;
				uint8_t *NNptr		=Nptr	-rstr;
				uint8_t *NNNptr		=NNptr	-rstr;
				uint8_t *NNNNptr	=NNNptr	-rstr;
				int
					NNNNWW	=NNNNptr	[-2*3],//y
					NNNN	=NNNNptr	[+0*3],
					NNNWWW	=NNNptr		[-3*3],
					NNN	=NNNptr		[+0*3],
					NNNE	=NNNptr		[+1*3],
					NNWW	=NNptr		[-2*3],
					NNW	=NNptr		[-1*3],
					NN	=NNptr		[+0*3],
					NNE	=NNptr		[+1*3],
					NWWW	=Nptr		[-3*3],
					NWW	=Nptr		[-2*3],
					NW	=Nptr		[-1*3],
					N	=Nptr		[+0*3],
					NE	=Nptr		[+1*3],
					NEE	=Nptr		[+2*3],
					NEEE	=Nptr		[+3*3],
					NEEEE	=Nptr		[+4*3],
					WWWWW	=bptr		[-5*3],
					WWWW	=bptr		[-4*3],
					WWW	=bptr		[-3*3],
					WW	=bptr		[-2*3],
					W	=bptr		[-1*3];
				int j;
#if 1
				(void)NNNNWW;
				(void)NNNN;
				(void)NNNWWW;
				(void)NNN;
				(void)NNNE;
				(void)NNWW;
				(void)NNW;
				(void)NN;
				(void)NNE;
				(void)NWWW;
				(void)NWW;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)NEEEE;
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
#endif
#ifdef _DEBUG
#define PRED(E0, E1, E2) stats[j]=a->stats+j*256*ELEVELS+(debug[j]=(E0)&(PLEVELS-1)); ++j;
#else
#define PRED(E0, E1, E2) stats[j]=a->stats+j*256*ELEVELS+((E0)&(PLEVELS-1)); ++j;
#endif
				j=0;
				PREDLIST;
#undef  PRED
#if defined _MSC_VER &&0
				if(ky==2962&&kx==2751)//
				{
					printf(" %d %d %d %d %d %d %d %d %d %d %d %d\n"
						, debug[0x0]
						, debug[0x1]
						, debug[0x2]
						, debug[0x3]
						, debug[0x4]
						, debug[0x5]
						, debug[0x6]
						, debug[0x7]
						, debug[0x8]
						, debug[0x9]
						, debug[0xA]
						, debug[0xB]
					);
					g_debug=1;//
				}
#endif
				*bptr=*imptr=codesym(&ac, stats, mixer, shift, *imptr, cc->fwd);
#if defined ENABLE_GUIDE
#if defined _MSC_VER &&0
					g_debug=0;//
#endif
				if(!cc->fwd)
				{
					uint8_t val=g_image[3*(cc->iw*(y1+ky)+x1+kx)+perm[0]];
					if(*bptr!=val)
						CRASH("guide  Y %d  X %d  C %d  0x%02X 0x%02X [original]", ky, kx, 0, *bptr, val);
				}
#endif
			}
		}
#ifdef PROFILE_SIZE
		g_kc=1;
#endif
		for(int ky=0;ky<dy;++ky)
		{
			int shift=SHIFT_EXPR;
			uint8_t *imptr, *bptr, *bptr1;
			
			imptr=cc->image+3*(cc->iw*(y1+ky)+x1)+perm[1];
			bptr=block+rstr*ky;
			bptr1=bptr+perm[0];
			bptr+=perm[1];
			for(int kx=0;kx<dx;++kx, imptr+=3, bptr+=3, bptr1+=3)
			{
				uint8_t *Nptr		=bptr	-rstr;
				uint8_t *NNptr		=Nptr	-rstr;
				uint8_t *NNNptr		=NNptr	-rstr;
				uint8_t *NNNNptr	=NNNptr	-rstr;
				uint8_t *Nptr1		=bptr1	-rstr;
				uint8_t *NNptr1		=Nptr1	-rstr;
				uint8_t *NNNptr1	=NNptr1	-rstr;
				uint8_t *NNNNptr1	=NNNptr1-rstr;
				int
					NNNNWW1	=NNNNptr1	[-2*3],//y
					NNNN1	=NNNNptr1	[+0*3],
					NNNWWW1	=NNNptr1	[-3*3],
					NNN1	=NNNptr1	[+0*3],
					NNNE1	=NNNptr1	[+1*3],
					NNWW1	=NNptr1		[-2*3],
					NNW1	=NNptr1		[-1*3],
					NN1	=NNptr1		[+0*3],
					NNE1	=NNptr1		[+1*3],
					NWWW1	=Nptr1		[-3*3],
					NWW1	=Nptr1		[-2*3],
					NW1	=Nptr1		[-1*3],
					N1	=Nptr1		[+0*3],
					NE1	=Nptr1		[+1*3],
					NEE1	=Nptr1		[+2*3],
					NEEE1	=Nptr1		[+3*3],
					NEEEE1	=Nptr1		[+4*3],
					WWWWW1	=bptr1		[-5*3],
					WWWW1	=bptr1		[-4*3],
					WWW1	=bptr1		[-3*3],
					WW1	=bptr1		[-2*3],
					W1	=bptr1		[-1*3],
					m1	=bptr1		[+0*3];
				int
					NNNNWW	=NNNNptr	[-2*3],//u
					NNNN	=NNNNptr	[+0*3],
					NNNWWW	=NNNptr		[-3*3],
					NNN	=NNNptr		[+0*3],
					NNNE	=NNNptr		[+1*3],
					NNWW	=NNptr		[-2*3],
					NNW	=NNptr		[-1*3],
					NN	=NNptr		[+0*3],
					NNE	=NNptr		[+1*3],
					NWWW	=Nptr		[-3*3],
					NWW	=Nptr		[-2*3],
					NW	=Nptr		[-1*3],
					N	=Nptr		[+0*3],
					NE	=Nptr		[+1*3],
					NEE	=Nptr		[+2*3],
					NEEE	=Nptr		[+3*3],
					NEEEE	=Nptr		[+4*3],
					WWWWW	=bptr		[-5*3],
					WWWW	=bptr		[-4*3],
					WWW	=bptr		[-3*3],
					WW	=bptr		[-2*3],
					W	=bptr		[-1*3];
				int j;
				
#if 1
				(void)NNNNWW1;
				(void)NNNN1;
				(void)NNNWWW1;
				(void)NNN1;
				(void)NNNE1;
				(void)NNWW1;
				(void)NNW1;
				(void)NN1;
				(void)NNE1;
				(void)NWWW1;
				(void)NWW1;
				(void)NW1;
				(void)N1;
				(void)NE1;
				(void)NEE1;
				(void)NEEE1;
				(void)NEEEE1;
				(void)WWWWW1;
				(void)WWWW1;
				(void)WWW1;
				(void)WW1;
				(void)W1;
				(void)NNNNWW;
				(void)NNNN;
				(void)NNNWWW;
				(void)NNN;
				(void)NNNE;
				(void)NNWW;
				(void)NNW;
				(void)NN;
				(void)NNE;
				(void)NWWW;
				(void)NWW;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)NEEEE;
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
#endif
#ifdef _DEBUG
#define PRED(E0, E1, E2) stats[j]=a->stats+j*256*ELEVELS+(debug[j]=(E1)&(PLEVELS-1)); ++j;
#else
#define PRED(E0, E1, E2) stats[j]=a->stats+j*256*ELEVELS+((E1)&(PLEVELS-1)); ++j;
#endif
				j=0;
				PREDLIST;
#undef  PRED
				*bptr=*imptr=codesym(&ac, stats, mixer, shift, *imptr, cc->fwd);
#if defined ENABLE_GUIDE
				if(!cc->fwd)
				{
					uint8_t val=g_image[3*(cc->iw*(y1+ky)+x1+kx)+perm[1]];
					if(*bptr!=val)
						CRASH("guide  Y %d  X %d  C %d  0x%02X 0x%02X [original]", ky, kx, 1, *bptr, val);
				}
#endif
			}
		}
#ifdef PROFILE_SIZE
		g_kc=2;
#endif
		for(int ky=0;ky<dy;++ky)
		{
			int shift=SHIFT_EXPR;
			uint8_t *imptr, *bptr, *bptr1, *bptr2;
			
			imptr=cc->image+3*(cc->iw*(y1+ky)+x1)+perm[2];
			bptr=block+rstr*ky;
			bptr1=bptr+perm[0];
			bptr2=bptr+perm[1];
			bptr+=perm[2];
			for(int kx=0;kx<dx;++kx, imptr+=3, bptr+=3, bptr1+=3, bptr2+=3)
			{
				uint8_t *Nptr		=bptr	-rstr;
				uint8_t *NNptr		=Nptr	-rstr;
				uint8_t *NNNptr		=NNptr	-rstr;
				uint8_t *NNNNptr	=NNNptr	-rstr;
				uint8_t *Nptr1		=bptr1	-rstr;
				uint8_t *NNptr1		=Nptr1	-rstr;
				uint8_t *NNNptr1	=NNptr1	-rstr;
				uint8_t *NNNNptr1	=NNNptr1-rstr;
				uint8_t *Nptr2		=bptr2	-rstr;
				uint8_t *NNptr2		=Nptr2	-rstr;
				uint8_t *NNNptr2	=NNptr2	-rstr;
				uint8_t *NNNNptr2	=NNNptr2-rstr;
				int
					NNNNWW1	=NNNNptr1	[-2*3],//y
					NNNN1	=NNNNptr1	[+0*3],
					NNNWWW1	=NNNptr1	[-3*3],
					NNN1	=NNNptr1	[+0*3],
					NNNE1	=NNNptr1	[+1*3],
					NNWW1	=NNptr1		[-2*3],
					NNW1	=NNptr1		[-1*3],
					NN1	=NNptr1		[+0*3],
					NNE1	=NNptr1		[+1*3],
					NWWW1	=Nptr1		[-3*3],
					NWW1	=Nptr1		[-2*3],
					NW1	=Nptr1		[-1*3],
					N1	=Nptr1		[+0*3],
					NE1	=Nptr1		[+1*3],
					NEE1	=Nptr1		[+2*3],
					NEEE1	=Nptr1		[+3*3],
					NEEEE1	=Nptr1		[+4*3],
					WWWWW1	=bptr1		[-5*3],
					WWWW1	=bptr1		[-4*3],
					WWW1	=bptr1		[-3*3],
					WW1	=bptr1		[-2*3],
					W1	=bptr1		[-1*3],
					m1	=bptr1		[+0*3];
				int
					NNNNWW2	=NNNNptr2	[-2*3],//u
					NNNN2	=NNNNptr2	[+0*3],
					NNNWWW2	=NNNptr2	[-3*3],
					NNN2	=NNNptr2	[+0*3],
					NNNE2	=NNNptr2	[+1*3],
					NNWW2	=NNptr2		[-2*3],
					NNW2	=NNptr2		[-1*3],
					NN2	=NNptr2		[+0*3],
					NNE2	=NNptr2		[+1*3],
					NWWW2	=Nptr2		[-3*3],
					NWW2	=Nptr2		[-2*3],
					NW2	=Nptr2		[-1*3],
					N2	=Nptr2		[+0*3],
					NE2	=Nptr2		[+1*3],
					NEE2	=Nptr2		[+2*3],
					NEEE2	=Nptr2		[+3*3],
					NEEEE2	=Nptr2		[+4*3],
					WWWWW2	=bptr2		[-5*3],
					WWWW2	=bptr2		[-4*3],
					WWW2	=bptr2		[-3*3],
					WW2	=bptr2		[-2*3],
					W2	=bptr2		[-1*3],
					m2	=bptr2		[+0*3];
				int
					NNNNWW	=NNNNptr	[-2*3],//v
					NNNN	=NNNNptr	[+0*3],
					NNNWWW	=NNNptr		[-3*3],
					NNN	=NNNptr		[+0*3],
					NNNE	=NNNptr		[+1*3],
					NNWW	=NNptr		[-2*3],
					NNW	=NNptr		[-1*3],
					NN	=NNptr		[+0*3],
					NNE	=NNptr		[+1*3],
					NWWW	=Nptr		[-3*3],
					NWW	=Nptr		[-2*3],
					NW	=Nptr		[-1*3],
					N	=Nptr		[+0*3],
					NE	=Nptr		[+1*3],
					NEE	=Nptr		[+2*3],
					NEEE	=Nptr		[+3*3],
					NEEEE	=Nptr		[+4*3],
					WWWWW	=bptr		[-5*3],
					WWWW	=bptr		[-4*3],
					WWW	=bptr		[-3*3],
					WW	=bptr		[-2*3],
					W	=bptr		[-1*3];
				int j;
				
#if 1
				(void)NNNNWW1;
				(void)NNNN1;
				(void)NNNWWW1;
				(void)NNN1;
				(void)NNNE1;
				(void)NNWW1;
				(void)NNW1;
				(void)NN1;
				(void)NNE1;
				(void)NWWW1;
				(void)NWW1;
				(void)NW1;
				(void)N1;
				(void)NE1;
				(void)NEE1;
				(void)NEEE1;
				(void)NEEEE1;
				(void)WWWWW1;
				(void)WWWW1;
				(void)WWW1;
				(void)WW1;
				(void)W1;
				(void)NNNNWW2;
				(void)NNNN2;
				(void)NNNWWW2;
				(void)NNN2;
				(void)NNNE2;
				(void)NNWW2;
				(void)NNW2;
				(void)NN2;
				(void)NNE2;
				(void)NWWW2;
				(void)NWW2;
				(void)NW2;
				(void)N2;
				(void)NE2;
				(void)NEE2;
				(void)NEEE2;
				(void)NEEEE2;
				(void)WWWWW2;
				(void)WWWW2;
				(void)WWW2;
				(void)WW2;
				(void)W2;
				(void)m2;
				(void)NNNNWW;
				(void)NNNN;
				(void)NNNWWW;
				(void)NNN;
				(void)NNNE;
				(void)NNWW;
				(void)NNW;
				(void)NN;
				(void)NNE;
				(void)NWWW;
				(void)NWW;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)NEEEE;
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
#endif
#ifdef _DEBUG
#define PRED(E0, E1, E2) stats[j]=a->stats+j*256*ELEVELS+(debug[j]=(E2)&(PLEVELS-1)); ++j;
#else
#define PRED(E0, E1, E2) stats[j]=a->stats+j*256*ELEVELS+((E2)&(PLEVELS-1)); ++j;
#endif
				j=0;
				PREDLIST;
#undef  PRED
				*bptr=*imptr=codesym(&ac, stats, mixer, shift, *imptr, cc->fwd);
#if defined ENABLE_GUIDE
				if(!cc->fwd)
				{
					uint8_t val=g_image[3*(cc->iw*(y1+ky)+x1+kx)+perm[2]];
					if(*bptr!=val)
						CRASH("guide  Y %d  X %d  C %d  0x%02X 0x%02X [original]", ky, kx, 2, *bptr, val);
				}
#endif
			}
		}
		if(cc->fwd)
		{
			*(uint64_t*)ac.ptr=ac.lo<<32|ac.lo>>32;//flush
			ac.ptr+=sizeof(uint64_t);
			cc->csizes[kb]=(int)(ac.ptr-ac.start);
#if defined _MSC_VER && defined MT_DISABLE
			if(ac.ptr>ac.end)
				CRASH("inflation  kb %d", kb);
#endif
		}
	}
	return 0;
}
int c64_codec(int argc, char **argv)
{
	const uint16_t tag='6'|'4'<<8;
	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	uint32_t c=0;
	CodecCtx _cc={0}, *cc=&_cc;
	ThreadCtx *args=0;
	void **thandles=0;
	size_t bsum=0;
	int blockdim=0;
#ifdef LOUD
	double t;
#endif

	if((uint32_t)(argc-3)>(uint32_t)(5-3))
	{
		printf(
			"C64: High-Efficiency Multithreaded Lossless Image Compressor.  AWM 2026.\n"
			"Usage:  \"%s\"  src  dst  [nThreads]  [blockSize]\n"
			"  nThreads:   Default 0 - number of logical cores,  pass 1 for single thread.\n"
			"  blockSize:  Default %d, minimum 16,  pass 0 to disable blocking.\n"
			"Examples:\n"
			"  \"%s\"  src  dst            Multithreaded, default block size.\n"
			"  \"%s\"  src  dst  1         Single thread, default block size.\n"
			"  \"%s\"  src  dst  1  0      No blocking.\n"
			, argv[0]
			, BLOCKDIM
			, argv[0]
			, argv[0]
			, argv[0]
		);
		return 1;
	}
#ifdef LOUD
	t=time_sec();
#endif
	srcfn=argv[1];
	dstfn=argv[2];
	cc->nthreads=argc>3?atoi(argv[3]):0;
	blockdim=argc>4?atoi(argv[4]):BLOCKDIM;
	if(!blockdim)
		blockdim=0x7FFFFFFF;
	if(blockdim<16)
	{
		printf("Warning: Invalid block dimension %d\n", blockdim);
		blockdim=BLOCKDIM;
	}
	cc->ncores=query_cpu_cores();
	if(!cc->nthreads)
		cc->nthreads=cc->ncores;
	if(cc->nthreads>cc->ncores)
		cc->nthreads=cc->ncores;
#ifdef MT_DISABLE
	cc->nthreads=1;
#endif
	
	fsrc=fopen(srcfn, "rb");
	if(!fsrc)
	{
		CRASH("Cannot open \"%s\"", srcfn);
		return 1;
	}
	fread(&c, 1, 2, fsrc);
	cc->fwd=c==('P'|'6'<<8);
	if(!cc->fwd&&c!=tag)
	{
		CRASH("Unsupported file \"%s\"", srcfn);
		return 1;
	}
	if(cc->fwd)//parse header
	{
		c=fgetc(fsrc);
		if(c!='\n')
		{
			CRASH("Unsupported file \"%s\"", srcfn);
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
		cc->iw=0;
		while((uint32_t)(c-'0')<10)
		{
			cc->iw=10*cc->iw+c-'0';
			c=fgetc(fsrc);
		}
		while(c<=' ')
			c=fgetc(fsrc);
		cc->ih=0;
		while((uint32_t)(c-'0')<10)
		{
			cc->ih=10*cc->ih+c-'0';
			c=fgetc(fsrc);
		}
		while(c=='#')
		{
			c=fgetc(fsrc);
			while(c!='\n')
				c=fgetc(fsrc);
			c=fgetc(fsrc);
		}
		c|=fgetc(fsrc)<<8*1;
		c|=fgetc(fsrc)<<8*2;
		c|=fgetc(fsrc)<<8*3;
		if(c!=(
			'\n'<<8*0|
			'2' <<8*1|
			'5' <<8*2|
			'5' <<8*3
		)||fgetc(fsrc)!='\n')
		{
			CRASH("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		cc->blockdim=cc->iw;
		if(cc->blockdim<cc->ih)
			cc->blockdim=cc->ih;
		if(cc->blockdim>blockdim)
			cc->blockdim=blockdim;
	}
	else
	{
		cc->iw=0;
		cc->ih=0;
		cc->blockdim=0;
		fread(&cc->iw, 1, 3, fsrc);
		fread(&cc->ih, 1, 3, fsrc);
		fread(&cc->blockdim, 1, 2, fsrc);
	}
	if(cc->iw<1||cc->ih<1)
	{
		CRASH("Invalid image %d*%d", cc->iw, cc->ih);
		return 1;
	}
	cc->usize=(size_t)3*cc->iw*cc->ih;
	cc->ccap=(size_t)4*cc->iw*cc->ih;
	cc->xblocks=(cc->iw+cc->blockdim-1)/cc->blockdim;
	cc->yblocks=(cc->ih+cc->blockdim-1)/cc->blockdim;
//	if(!cc->xblocks)cc->xblocks=1;
//	if(!cc->yblocks)cc->yblocks=1;
	cc->nblocks=cc->xblocks*cc->yblocks;
	if(cc->nthreads>cc->nblocks)
		cc->nthreads=cc->nblocks;
	cc->csizes=(int*)malloc(cc->nblocks*sizeof(int));
	cc->cbounds=(size_t*)malloc(((size_t)cc->nblocks+1)*sizeof(size_t));
	if(!cc->csizes||!cc->cbounds)
	{
		CRASH("Alloc error");
		return 1;
	}
	if(!cc->fwd)
	{
		bsum=0;
		for(int k=0;k<cc->nblocks;++k)
		{
			cc->csizes[k]=0;
			fread(cc->csizes+k, 1, 4, fsrc);
			cc->cbounds[k]=bsum;
			bsum+=cc->csizes[k];
		}
		cc->csize=cc->cbounds[cc->nblocks]=bsum;
	}
	cc->blockcap=(int)sizeof(uint8_t[3])*(cc->blockdim<cc->iw?cc->blockdim:cc->iw)*((cc->blockdim<cc->ih?cc->blockdim:cc->ih)+16);//zero-pad 16 rows
#if defined _MSC_VER &&1
	if(cc->fwd)
		printf("WH %d*%d  %.2f MB  %.2f MP  blockdim %d  blockcap %d  csize %d\n"
			, cc->iw
			, cc->ih
			, 3.*cc->iw*cc->ih/(1<<20)
			, cc->iw*cc->ih*1e-6
			, cc->blockdim
			, cc->blockcap
			, (int)cc->csize
		);
#endif
	cc->image=(uint8_t*)malloc(cc->usize);
	cc->streams=(uint8_t*)malloc(cc->ccap);
	args=(ThreadCtx*)malloc(cc->nthreads*sizeof(ThreadCtx));
	thandles=(void**)malloc(cc->nthreads*sizeof(void*));
	if(!cc->image||!cc->streams||!cc->cbounds||!args)
	{
		free(cc->image);
		free(cc->streams);
		free(cc->csizes);
		free(cc->cbounds);
		free(args);
		CRASH("Alloc error");
		return 1;
	}
#ifdef PROFILE_SIZE
	for(int kp=0;kp<1<<PROBBITS;++kp)
		shannontable[kp]=-log2((kp+!kp)*(1./(1<<PROBBITS)))/8;
	memset(esizes, 0, sizeof(esizes));
#endif
	for(int k=0;k<_countof(squashtable);++k)
	{
		int p0=k-_countof(squashtable)/2;

		int signmask=p0>>31;
		p0<<=1;
		p0^=signmask;
		p0-=signmask;
		p0=_countof(squashtable)-p0;
		p0=(uint64_t)p0*p0>>TABLEBITS;
		p0=(uint64_t)p0*p0>>(2*TABLEBITS-PROBBITS);
		p0=(1<<PROBBITS)-p0;
		p0^=signmask;
		p0-=signmask;
		p0/=2;
		p0+=1<<PROBBITS>>1;
		CLAMP(p0, 1, (1<<PROBBITS)-1);
		squashtable[k]=p0;
#if defined _MSC_VER && defined PRINT_SQUASH
		if(cc->fwd)//
			printf(" %5d:%6d,", k, p0);
#endif
	}
	if(cc->fwd)
	{
#if defined _MSC_VER && defined PRINT_SQUASH
		printf("\n");//
#endif
		fread(cc->image, 1, cc->usize, fsrc);
		guide_save(cc->image, cc->iw, cc->ih);
		memset(cc->csizes, 0, cc->nblocks*sizeof(int));
		bsum=0;
		for(int ky=0, y1=0, idx=0;ky<cc->yblocks;++ky)
		{
			int y2=y1+cc->blockdim, dy;

		//	if(y2+cc->blockdim>cc->ih)y2=cc->ih;
			if(y2>cc->ih)y2=cc->ih;
			dy=y2-y1;
			for(int kx=0, x1=0;kx<cc->xblocks;++kx, ++idx)
			{
				int x2=x1+cc->blockdim, dx, bsize;

			//	if(x2+cc->blockdim>cc->iw)x2=cc->iw;
				if(x2>cc->iw)x2=cc->iw;
				dx=x2-x1;

				bsize=4*dx*dy;
				cc->cbounds[idx]=bsum;
				bsum+=bsize;

				x1=x2;
			}
			y1=y2;
		}
		cc->cbounds[cc->nblocks]=bsum;
	}
	else
		fread(cc->streams, 1, cc->csize, fsrc);
	{
		int bq=cc->nblocks/cc->nthreads, br=cc->nblocks%cc->nthreads, bstop=0;

		memset(args, 0, cc->nthreads*sizeof(ThreadCtx));
		for(int kt=0;kt<cc->nthreads;++kt)
		{
			ThreadCtx *arg=args+kt;
			int tblocks=bq+(kt<br);

			arg->cc=cc;
			arg->tid=kt;
			arg->b1=bstop;
			bstop+=tblocks;
			arg->b2=bstop;
			arg->block=(uint8_t*)malloc(cc->blockcap);
			if(!arg->block)
			{
				CRASH("Alloc error");
				return 1;
			}
			memset(arg->block, 0, cc->blockcap);
		}
	}
	if(cc->nthreads==1)
		codec_thread(args);
	else
	{
		for(int kt=0;kt<cc->nthreads;++kt)
		{
			ThreadCtx *arg=args+kt;
			int error;

#ifdef _MSC_VER
			thandles[kt]=(void*)_beginthreadex(0, 0, codec_thread, arg, 0, 0);
			error=!thandles[kt];
#else
			error=pthread_create((pthread_t*)(thandles+kt), 0, codec_thread, arg);
#endif
			if(error)
			{
				CRASH("Alloc error");
				return 1;
			}
		}
#ifdef _MSC_VER
		WaitForMultipleObjects(cc->nthreads, (HANDLE*)thandles, TRUE, INFINITE);
		for(int kt=0;kt<cc->nthreads;++kt)
			CloseHandle(thandles[kt]);
#else
		for(int kt=0;kt<cc->nthreads;++kt)
			pthread_join((pthread_t)thandles[kt], 0);
#endif
	}
	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing");
		return 1;
	}
	if(cc->fwd)
	{
		cc->csize=0;
		cc->csize+=fwrite(&tag, 1, 2, fdst);
		cc->csize+=fwrite(&cc->iw, 1, 3, fdst);
		cc->csize+=fwrite(&cc->ih, 1, 3, fdst);
		cc->csize+=fwrite(&cc->blockdim, 1, 2, fdst);
		for(int kb=0;kb<cc->nblocks;++kb)
			cc->csize+=fwrite(cc->csizes+kb, 1, 4, fdst);
		for(int kb=0;kb<cc->nblocks;++kb)
			cc->csize+=fwrite(cc->streams+cc->cbounds[kb], 1, cc->csizes[kb], fdst);
	}
	else
	{
		fprintf(fdst, "P6\n%d %d\n255\n", cc->iw, cc->ih);
		fwrite(cc->image, 1, cc->usize, fdst);
	}
	fclose(fdst);
	for(int kt=0;kt<cc->nthreads;++kt)
	{
		ThreadCtx *arg=args+kt;

		free(arg->block);
	}
	free(cc->image);
	free(cc->streams);
	free(cc->csizes);
	free(cc->cbounds);
	free(thandles);
	free(args);
#ifdef ENABLE_GUIDE
	if(!cc->fwd)
		free(g_image);
#endif
#ifdef LOUD
	t=time_sec()-t;
#ifdef PROFILE_SIZE
	if(cc->fwd)
	{
		double ctotal[3]={0};

		printf("bitplane  TYUV\n");
		for(int kb=7;kb>=0;--kb)
		{
			ctotal[0]+=esizes[8*0+kb];
			ctotal[1]+=esizes[8*1+kb];
			ctotal[2]+=esizes[8*2+kb];
		}
		printf("T   %12.2lf  %12.2lf  %12.2lf  %12.2lf\n\n"
			, ctotal[0]+ctotal[1]+ctotal[2]
			, ctotal[0]
			, ctotal[1]
			, ctotal[2]
		);
		for(int kb=7;kb>=0;--kb)
		{
			printf("B%d  %12.2lf  %12.2lf  %12.2lf  %12.2lf\n"
				, kb
				, esizes[8*0+kb]+esizes[8*1+kb]+esizes[8*2+kb]
				, esizes[8*0+kb]
				, esizes[8*1+kb]
				, esizes[8*2+kb]
			);
		}
	}
#endif
	{
		struct stat info={0};
		int error=0;

		error=stat(srcfn, &info);
		cc->usize=info.st_size;
		error=stat(dstfn, &info);
		cc->csize=info.st_size;
		if(!cc->fwd)
		{
			size_t tmp=cc->usize;
			cc->usize=cc->csize;
			cc->csize=tmp;
		}
	}
	if(cc->fwd)
	{
		printf("%9td->%9td  %8.4lf%%  %12.6lf:1  BPD %12.6lf\n"
			, cc->usize
			, cc->csize
			, 100.*cc->csize/cc->usize
			, (double)cc->usize/cc->csize
			, 8.*cc->csize/cc->usize
		);
	}
	printf("%c  %12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB\n"
		, 'D'+cc->fwd
		, t
		, cc->usize/(t*1024*1024)
		, t*1024*1024*1000/cc->usize
	);
#endif
	(void)&time_sec;
	(void)&print_rct;
	return 0;
}
