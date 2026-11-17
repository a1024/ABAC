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

//config
#ifdef _MSC_VER
	#define LOUD
	#define ENABLE_GUIDE

//	#define ANS_VAL
#endif


//	#define PROFILE_SIZE


enum
{
	//LBITS:  log2 lookup table size
	//EBITS:  log2 number of collision cells
#if 1
	//lz+fse
	LBITS=14, EBITS=8,	//16 MB
//	LBITS=10, EBITS=8,	//256 KB	prev
//	LBITS=12, EBITS=10,
#elif 1
	//lzac fast
	LBITS=14,
	EBITS=3,
#else
	//lzac high
	LBITS=22,
	EBITS=6,
#endif

	LZ_MINRUN=3*4,

	CONFIG_EXP=7,
	CONFIG_MSB=3,
	CONFIG_LSB=0,

	LITCTX=2,
	NHISTS=3*LITCTX+3,//[3 ch]*LITCTX + [litlen  offset  matchlen]

	PROBBITS=12,
};
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
#ifdef ANS_VAL
#define ANS_VAL_HISTSIZE 128
typedef struct _ANSVALHeader
{
	unsigned short esize, count;
	unsigned idx;
	struct _ANSVALHeader *above, *below;
	uint8_t data[];
} ANSVALNode;
static ANSVALNode *debugstack=0;
static int ansvalidx=0, ansvalmax=0;
static void ansval_push(const void *data, int esize, int count)
{
	int size=count*esize;
	ANSVALNode *node=(ANSVALNode*)malloc(sizeof(ANSVALNode)+size);
	if(!node)
	{
		printf("Alloc error\n");
		exit(1);
	}
	memset(node, 0, sizeof(ANSVALNode)+size);
	node->esize=esize;
	node->count=count;
	node->idx=ansvalidx++;
	node->above=0;
	node->below=debugstack;
	if(debugstack)
		debugstack->above=node;
	memcpy(node->data, data, size);
	debugstack=node;
	++ansvalmax;
}
static void ansval_printr(const void *data, int esize, int count, const void *xdata)//print elements in reverse because little-endian
{
	const uint8_t *p=(const uint8_t*)data, *p2=(const uint8_t*)xdata;
	int size=count*esize, k;
	for(k=0;k<size;k+=esize)
	{
		int k2=esize-1;
		printf(" ");
		for(;k2>=0;--k2)
		{
			int val=p[k+k2];
			if(p2)
				val^=p2[k+k2];
			if(p2&&!val)
				printf("--");
			else
				printf("%02X", val);
		}
	}
	printf("\n");
}
static void ansval_check(const void *data, int esize, int count)
{
	--ansvalidx;
	if(!debugstack)
	{
		printf("Debug stack is empty\n");
		ansval_printr(data, esize, count, 0);
		CRASH("");
	}
	else if(debugstack->esize!=esize||debugstack->count!=count||memcmp(data, debugstack->data, esize*count))
	{
		printf("\n\nValidation Error  [enc ^ | v dec]\n");
		if(debugstack->above)
		{
			ANSVALNode *node=debugstack->above;
			if(node->above)
			{
				ANSVALNode *node2=node->above;
				printf("[%10d] Verified:   ", node2->idx);
				ansval_printr(node2->data, node2->esize, node2->count, 0);
			}
			printf("[%10d] Verified:   ", node->idx);
			ansval_printr(node->data, node->esize, node->count, 0);
			printf("\n");
		}

		printf("[%10d] Original:   ", debugstack->idx);
		ansval_printr(debugstack->data, esize, count, 0);

		printf("[%10d] Corrupt:    ", debugstack->idx);
		ansval_printr(data, esize, count, 0);
		
		if(debugstack->esize==esize&&debugstack->count==count)
		{
			printf("[%10d] XOR:        ", debugstack->idx);
			ansval_printr(debugstack->data, esize, count, data);
		}
		if(debugstack->below)
		{
			ANSVALNode *node=debugstack->below;
			printf("\n");
			printf("[%10d] Below:      ", node->idx);
			ansval_printr(node->data, node->esize, node->count, 0);
			if(node->below)
			{
				node=node->below;
				printf("[%10d] Below:      ", node->idx);
				ansval_printr(node->data, node->esize, node->count, 0);
			}
		}
		printf("\n\n");
		CRASH("");
	}
	if(debugstack->below)
		debugstack=debugstack->below;
}
#else
#define ansval_push(...)
#define ansval_check(...)
#endif

//RCT
#if 1
typedef struct _RCTInfo
{
	uint8_t pidx, uc0, vc0, vc1;
} RCTInfo;
static const uint8_t perms[]=
{
	0, 1, 2,	0, 2, 1,
	1, 2, 0,	1, 0, 2,
	2, 0, 1,	2, 1, 0,
};
enum
{
	RCTBITS=2,
	RCTMAX=1<<RCTBITS,
	RCTLEVELS=RCTMAX+1,//+unity
	NPERMS=_countof(perms)/3,
	OCH_COUNT=NPERMS*(1+RCTLEVELS+RCTLEVELS*(RCTLEVELS+1)/2),
};
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
static void crct_analysis(uint8_t *image, int iw, int ih, RCTInfo *ret_rct)
{
	enum
	{
		STRIDE=7,
	};
	int64_t counters[OCH_COUNT]={0}, bestsum=0;
	RCTInfo rct={0};
	int rstr=0;
#ifdef LOUD
	double t=time_sec();
#endif

	rstr=3*iw;
	for(int ky=1;ky<=ih-STRIDE;ky+=STRIDE)
	{
		uint8_t *imptr=image+rstr*ky+3;
		for(int kx=1;kx<=iw-STRIDE;kx+=STRIDE, imptr+=3*STRIDE)
		{
			int64_t *ctrptr=counters;
			int rgb[]=
			{
				(imptr[0]-imptr[-3+0]-imptr[-rstr+0]+imptr[-rstr-3+0])<<RCTBITS,
				(imptr[1]-imptr[-3+1]-imptr[-rstr+1]+imptr[-rstr-3+1])<<RCTBITS,
				(imptr[2]-imptr[-3+2]-imptr[-rstr+2]+imptr[-rstr-3+2])<<RCTBITS,
			};
			for(int kp=0;kp<NPERMS;++kp)
			{
				int yuv[]=
				{
					rgb[perms[3*kp+0]],
					rgb[perms[3*kp+1]],
					rgb[perms[3*kp+2]],
				};
				*ctrptr++ += abs(yuv[0]);
				for(int uc0=0;uc0<=RCTMAX;++uc0)
					*ctrptr++ += abs(yuv[1]-(yuv[0]*uc0>>RCTBITS));
				for(int vc0=0;vc0<=RCTMAX;++vc0)
				{
					for(int vc1=0;vc0+vc1<=RCTMAX;++vc1)
						*ctrptr++ += abs(yuv[2]-((yuv[0]*vc0+yuv[1]*vc1)>>RCTBITS));
				}
			}
		}
	}

	bestsum=0;
	for(int kp=0;kp<NPERMS;++kp)
	{
		int64_t *currctrs=counters+(1+RCTLEVELS+RCTLEVELS*(RCTLEVELS+1)/2)*kp;
		for(int uc0=0;uc0<=RCTMAX;++uc0)
		{
			for(int vc0=0, idx=0;vc0<=RCTMAX;++vc0)
			{
				for(int vc1=0;vc0+vc1<=RCTMAX;++vc1, ++idx)
				{
					int64_t sum=currctrs[0]+currctrs[1+uc0]+currctrs[1+RCTLEVELS+idx];
					if(!bestsum||bestsum>sum)
					{
						bestsum=sum;
						rct.pidx=kp;
						rct.uc0=uc0;
						rct.vc0=vc0;
						rct.vc1=vc1;
#ifdef PRINT_ANALYSIS
						print_rct(&rct);
						printf("  %16lld\n", sum);
#endif
					}
				}
			}
		}
	}
	*ret_rct=rct;
#ifdef LOUD
	t=time_sec()-t;
	{
		int64_t usize=(int64_t)3*iw*ih;
		printf("A  %12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB\n"
			, t
			, usize/(t*1024*1024)
			, t*1024*1024*1000/usize
		);
	}
#endif
}
#endif

static const uint64_t g_bitmask[]=
{
	0x0,			0x1,			0x3,			0x7,
	0xF,			0x1F,			0x3F,			0x7F,
	0xFF,			0x1FF,			0x3FF,			0x7FF,
	0xFFF,			0x1FFF,			0x3FFF,			0x7FFF,
	0xFFFF,			0x1FFFF,		0x3FFFF,		0x7FFFF,
	0xFFFFF,		0x1FFFFF,		0x3FFFFF,		0x7FFFFF,
	0xFFFFFF,		0x1FFFFFF,		0x3FFFFFF,		0x7FFFFFF,
	0xFFFFFFF,		0x1FFFFFFF,		0x3FFFFFFF,		0x7FFFFFFF,
	0xFFFFFFFF,		0x1FFFFFFFF,		0x3FFFFFFFF,		0x7FFFFFFFF,
	0xFFFFFFFFF,		0x1FFFFFFFFF,		0x3FFFFFFFFF,		0x7FFFFFFFFF,
	0xFFFFFFFFFF,		0x1FFFFFFFFFF,		0x3FFFFFFFFFF,		0x7FFFFFFFFFF,
	0xFFFFFFFFFFF,		0x1FFFFFFFFFFF,		0x3FFFFFFFFFFF,		0x7FFFFFFFFFFF,
	0xFFFFFFFFFFFF,		0x1FFFFFFFFFFFF,	0x3FFFFFFFFFFFF,	0x7FFFFFFFFFFFF,
	0xFFFFFFFFFFFFF,	0x1FFFFFFFFFFFFF,	0x3FFFFFFFFFFFFF,	0x7FFFFFFFFFFFFF,
	0xFFFFFFFFFFFFFF,	0x1FFFFFFFFFFFFFF,	0x3FFFFFFFFFFFFFF,	0x7FFFFFFFFFFFFFF,
	0xFFFFFFFFFFFFFFF,	0x1FFFFFFFFFFFFFFF,	0x3FFFFFFFFFFFFFFF,	0x7FFFFFFFFFFFFFFF,
};
//LIFO Bit Packer inspired by FPC	https://github.com/kagiannis/FPC
typedef struct _BitPackerLIFO
{
	int64_t bitidx;
	uint64_t cache;
	uint8_t *start, *ptr, *end;//fwd enc / bwd dec
} BitPackerLIFO;
INLINE void bitpacker_enc(BitPackerLIFO *ec, int64_t bitcount, uint64_t sym)//1 <= bitcount <= 56
{
	ansval_push(&sym, sizeof(sym), 1);
	ec->cache|=sym<<ec->bitidx;
	ec->bitidx+=bitcount;
	*(uint64_t*)ec->ptr=ec->cache;
	ec->cache>>=ec->bitidx&56;
	ec->ptr+=ec->bitidx>>3;
	ec->bitidx&=7;
}
INLINE uint64_t bitpacker_dec(BitPackerLIFO *ec, int64_t bitcount)
{
	uint64_t sym;

	ec->bitidx-=bitcount;
	sym=((uint64_t*)ec->ptr)[-1]>>ec->bitidx&((1ULL<<bitcount)-1);
	ansval_check(&sym, sizeof(sym), 1);
	ec->ptr-=ec->bitidx>>3^7;
	ec->bitidx|=56;
	return sym;
}
static void bitpacker_enc_init(BitPackerLIFO *ec, uint8_t *start, uint8_t *end)
{
	ec->bitidx=0;
	ec->cache=0;
	ec->start=start;
	ec->ptr=start+sizeof(uint64_t);
	ec->end=end;
	*(uint64_t*)start=0;
}
static uint8_t* bitpacker_flush(BitPackerLIFO *ec)
{
	bitpacker_enc(ec, 1, 1);
	++ec->ptr;
	return ec->ptr;
}
static void bitpacker_dec_init(BitPackerLIFO *ec, uint8_t *start, uint8_t *end)
{
	ec->start=start;
	ec->ptr=end;
	ec->end=end;
	ec->bitidx=64-LZCNT64(((uint64_t*)ec->ptr)[-1]);
	bitpacker_dec(ec, 1);
}

//FSE		https://github.com/Cyan4973/FiniteStateEntropy
#if 1
typedef struct _FSE_CState_t
{
	ptrdiff_t value;
	const void *stateTable, *symbolTT;
	uint32_t stateLog;
} FSE_CState_t;
typedef struct _FSE_DState_t
{
	size_t state;
	const void *table;
} FSE_DState_t;
typedef struct _FSE_symbolCompressionTransform
{
	int32_t deltaFindState;
	uint32_t deltaNbBits;
} FSE_symbolCompressionTransform;
typedef uint32_t FSE_CTable;
typedef uint32_t FSE_DTable;
typedef struct _FSE_DTableHeader
{
	uint16_t tableLog, fastMode;
} FSE_DTableHeader;
typedef struct _FSE_decode_t
{
	uint16_t newState;
	uint8_t symbol, nbBits;
} FSE_decode_t;
#define FSE_MAX_SYMBOL_VALUE 255
#define FSE_FUNCTION_TYPE uint8_t
//#define FSE_TABLELOG_ABSOLUTE_MAX 15
#define FSE_MAX_MEMORY_USAGE 14
//#define FSE_DEFAULT_MEMORY_USAGE 13
//#define FSE_MIN_TABLELOG 5
#define FSE_MAX_TABLELOG (FSE_MAX_MEMORY_USAGE-2)
//#define FSE_DEFAULT_TABLELOG (FSE_DEFAULT_MEMORY_USAGE-2)
#define FSE_TABLESTEP(tableSize) (((tableSize)>>1) + ((tableSize)>>3) + 3)
//#define FSE_CTABLE_SIZE_U32(maxTableLog, maxSymbolValue)	(1+(1<<((maxTableLog)-1))+(((maxSymbolValue)+1)*2))
//#define FSE_DTABLE_SIZE_U32(maxTableLog)			(1+(1<<(maxTableLog)))
#define FSE_DECODE_TYPE FSE_decode_t
#endif
//end of copied FSE declarations

static void normalizehist(const uint32_t *hist, uint16_t *nhist, int probbits)
{
	int hsum=0, nusedlevels=0;
	for(int ks=0;ks<256;++ks)
	{
		int freq=hist[ks];
		hsum+=freq;
		nusedlevels+=freq!=0;
	}
	int64_t rsum=(((1LL<<probbits)-nusedlevels)<<24)/hsum;
	uint16_t CDF[257]={0};
	for(int ks=0, c=0, c2=0;ks<256;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*rsum>>24)+c2;
		c+=freq;
		c2+=freq!=0;
	}
	CDF[256]=1<<probbits;
	for(int ks=0;ks<256;++ks)
		nhist[ks]=CDF[ks+1]-CDF[ks];
}
static void enc_packhist(BitPackerLIFO *ec, const int *hist, uint64_t bypassmask, int ctxidx, int userans, int probbits)//histogram must be normalized to PROBBITS, with spike at 128
{
	int sum, cdfW, CDFlevels, startsym, ks;
	uint16_t hist2[256], CDF[257];

	if(bypassmask>>ctxidx&1)
		return;
	sum=0;
	if(userans)//32-bit hist
	{
		for(ks=0;ks<256;++ks)
			hist2[ks]=hist[ks];
	}
	else//16-bit hist
		memcpy(hist2, hist, sizeof(hist2));
	for(ks=0;ks<256;++ks)//integrage to zigzag CDF to be packed backwards
	{
		int sym=(uint8_t)((ks>>1^ks<<31>>31)+128);
		int freq=hist2[sym];
		CDF[ks]=sum;//separate buffer for faster access in 2nd loop
		sum+=freq;
	}
	CDF[256]=1<<probbits;
	
	cdfW=CDF[0];
	CDFlevels=1<<probbits;
	startsym=0;
	for(ks=1;ks<=256;++ks)//push GR.k
	{
		int next=CDF[ks], freq=next-cdfW;
		int nbypass=31-LZCNT32(CDFlevels);
		if(ks>1)
			nbypass-=7;
		if(nbypass<0)
			nbypass=0;
		CDF[ks]=nbypass<<probbits|freq;
		cdfW=next;
		CDFlevels-=freq;
		startsym=ks;
		if(!CDFlevels)
			break;
	}
	for(ks=startsym;ks>0;--ks)//encode GR
	{
		int freq=CDF[ks], nbypass=freq>>probbits;
		freq&=(1<<probbits)-1;
		int nzeros=freq>>nbypass, bypass=freq&((1<<nbypass)-1);
#ifdef ANS_VAL
		//ansval_push(&freq, sizeof(freq), 1);
#endif
		if(nbypass)
			bitpacker_enc(ec, nbypass, bypass);
		bitpacker_enc(ec, 1, 1);
		while(nzeros)
		{
			bitpacker_enc(ec, 1, 0);
			--nzeros;
		}
#ifdef ANS_VAL
		//ansval_push(&ks, sizeof(ks), 1);
#endif
	}
}
static void dec_unpackhist(BitPackerLIFO *ec, uint64_t bypassmask, int ctxidx, uint16_t *hist, int probbits)
{
	if(bypassmask>>ctxidx&1)//rare context
	{
		for(int ks=0;ks<256;++ks)//bypass
			hist[ks]=(1<<probbits)/256;
	}
	else
	{
		uint16_t CDF[257]={0};
		int CDFlevels=1<<probbits;
		CDF[0]=0;
		for(int ks=0;ks<256;++ks)//decode GR
		{
			int freq=-1;//stop bit doesn't count
			int nbypass=31-LZCNT32(CDFlevels);
			int ks2=ks+1;
			if(ks2>1)
				nbypass-=7;
			if(nbypass<0)
				nbypass=0;
#ifdef ANS_VAL
			//ansval_check(&ks2, sizeof(ks2), 1);
#endif
			int bit=0;
			do
			{
				bit=(int)bitpacker_dec(ec, 1);
				++freq;
			}while(!bit);
			if(nbypass)
				freq=freq<<nbypass|(int)bitpacker_dec(ec, nbypass);
#ifdef ANS_VAL
			//ansval_check(&freq, sizeof(freq), 1);
#endif

			CDF[ks]=freq;
			CDFlevels-=freq;
			if(CDFlevels<=0)
			{
#ifdef _DEBUG
				if(CDFlevels<0)
					CRASH("CDF unpack error");
#endif
				break;
			}
		}
		if(CDFlevels)
			CRASH("CDF unpack error");
		for(int ks=0;ks<256;++ks)//undo zigzag
		{
			int sym=(uint8_t)((ks>>1^ks<<31>>31)+128);
			hist[sym]=CDF[ks];
		}
	}
	//int sum=0;
	//for(int ks=0;ks<256;++ks)//integrate
	//{
	//	int freq=hist[ks];
	//	hist[ks]=sum;
	//	sum+=freq;
	//}
	//hist[256]=1<<PROBBITS;
	//for(int ks=0;ks<256;++ks)//CDF2sym contains {freq, (state&0xFFF)-cdf, sym}
	//{
	//	int cdf=hist[ks], next=hist[ks+1], freq=next-cdf;
	//	int val=(freq<<PROBBITS|0)<<8|ks;
	//	for(int ks2=cdf;ks2<next;++ks2, val+=1<<8)
	//		CDF2sym[ks2]=val;
	//}
}

typedef struct _ETable
{
	int32_t etable[1<<EBITS], estart, eend, ecount;
} ETable;
typedef struct _LZPacket
{
	int64_t nliterals, matchlen, offset;
} LZPacket;
typedef struct _LZFSETables
{
	FSE_symbolCompressionTransform symbolTT[NHISTS][256];
	uint16_t stateTable[NHISTS][1<<PROBBITS];
	int32_t hists[NHISTS][256];
	int16_t nhists[NHISTS][256];
	uint8_t wksp[1<<14>>2];
} LZFSETables;
typedef struct _LZFSDTables
{
	FSE_decode_t decodeTable[NHISTS][1<<PROBBITS];
} LZFSDTables;
static void lz_findmatch(ETable *tables, uint8_t *image, uint8_t *imptr, uint8_t *searchend, uint32_t *ret_matchidx, uint32_t *ret_matchlen)
{
	ETable *table;
	int32_t matchidx=0, matchlen=0;
	int ctr;

	{
		uint64_t a;

		a=0;
		a^=((uint64_t*)imptr)[0]*0x9E3779B97F4A7C15>>(64-LBITS);//24 bytes
		a^=((uint64_t*)imptr)[1]*0xC2B2AE3D27D4EB4F>>(64-LBITS);
		a^=((uint64_t*)imptr)[2]*0x165667B19E3779F9>>(64-LBITS);
	//	a^=((uint64_t*)imptr)[3]*0x85EBCA77C2B2AE63>>(64-LBITS);
		table=tables+a;
	}
//	table=tables+(*(uint64_t*)imptr*0x9E3779B97F4A7C15>>(64-LBITS));
	ctr=table->ecount;
	while(ctr>0)
	{
		int tidx=(table->estart+ctr-1)&((1<<EBITS)-1);
		int idx=table->etable[tidx], len=0;
		uint8_t *search1=image+idx, *search2=imptr;
		uint64_t cmp1=~0ULL, c1;

		if(*(uint32_t*)search1==*(uint32_t*)search2)
		{
			enum
			{
				WINDOWSIZE=0x400000,
			};
			for(;;)
			{
				cmp1=*(uint64_t*)search1^*(uint64_t*)search2;
				if(cmp1||search2>=searchend)
					break;
				search1+=sizeof(uint64_t);
				search2+=sizeof(uint64_t);
			}
			c1=TZCNT64(cmp1)>>3;
			len+=(int)(search2-imptr+c1);
			if(len>searchend-imptr+sizeof(uint64_t))
				len=(int)(searchend-imptr+sizeof(uint64_t));

			if(imptr-image-idx<WINDOWSIZE&&(len>matchlen||(len==matchlen&&idx>matchidx)))
				matchidx=idx, matchlen=len;
		}
		--ctr;
	}
	table->etable[table->eend]=(int32_t)(imptr-image);
	table->eend=(table->eend+1)&((1<<EBITS)-1);
	if(table->ecount>=1<<EBITS)
		table->estart=(table->estart+1)&((1<<EBITS)-1);
	else
		++table->ecount;
	*ret_matchidx=matchidx;
	*ret_matchlen=matchlen;
}
INLINE int split_uint(int64_t val, int *ret_nbits, int64_t *ret_bypass)
{
	int token, nbits, bypass;

	if(val<(int64_t)1<<CONFIG_EXP)
	{
		token=(int)val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		int64_t mantissa;
		int msb;

		msb=63^(int)LZCNT64(val);
		mantissa=val-((int64_t)1<<msb);
		token = (1<<CONFIG_EXP) + (
			(msb-CONFIG_EXP)<<(CONFIG_MSB+CONFIG_LSB)|
			(int)(mantissa>>(msb-CONFIG_MSB))<<CONFIG_LSB|
			(int)(mantissa&((1<<CONFIG_LSB)-1))
		);
		nbits=msb-(CONFIG_MSB+CONFIG_LSB);
		bypass=val>>CONFIG_LSB&(((int64_t)1<<nbits)-1);
	}
	if(ret_nbits)*ret_nbits=nbits;
	if(ret_bypass)*ret_bypass=bypass;
	return token;
}
int c63_codec(int argc, char **argv)
{
	enum
	{
		DEPTH=8,

		XPAD=8,
		NCH=3,
		NROWS=4,
		NVAL=4,
	};
	typedef int16_t Pixel_t;
	const uint16_t tag='6'|'3'<<8;
	
	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	uint32_t c=0;
	int fwd=0, iw=0, ih=0;
	int64_t usize=0, csize=0, streamsize=0;
	uint8_t *image=0, *imptr=0, *imptr0, *imend, *stream=0, *streamend=0;
	int psize=0;
	Pixel_t *pixels=0;
	RCTInfo rct={0};
	const uint8_t *perm=0;
	int kx, ky, kc;
	int yuv[3]={0};
	int degenmask=0, bypassmask=0;
	BitPackerLIFO ec={0};
	uint32_t states[3]={0};
	LZPacket *packets, *packet;
	int npackets;
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
		csize=(int64_t)3*iw*ih;
	}
	else
	{
		struct stat info={0};

		stat(srcfn, &info);
		csize=info.st_size;
		iw=0;
		ih=0;
		fread(&iw, 1, 3, fsrc);
		fread(&ih, 1, 3, fsrc);
		fread(&rct, 1, sizeof(rct), fsrc);
	}
	if(iw<1||ih<1)
	{
		CRASH("Unsupported file \"%s\"", srcfn);
		return 1;
	}
	usize=(int64_t)3*iw*ih;
	image=(uint8_t*)malloc(usize);
	stream=(uint8_t*)malloc(csize);
	psize=(iw+2*XPAD)*(int)sizeof(int32_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	streamend=stream+csize;
	if(fwd)
	{
		fread(image, 1, usize, fsrc);
		guide_save(image, iw, ih);
	}
	else
	{
		streamsize=csize-ftell(fsrc);
		fread(stream, 1, streamsize, fsrc);
		streamend=stream+streamsize;
	}
	fclose(fsrc);
	//analysis
	if(fwd)
		crct_analysis(image, iw, ih, &rct);
	perm=perms+3*rct.pidx;
	memset(pixels, 0, psize);
	if(fwd)
	{
		uint64_t usize;
		uint8_t *searchend;
		ETable *tables;
		uint32_t matchidx, matchlen;
		LZFSETables *fse;
		int packetscap;
		int rstr;
		int encbytes;
		FSE_symbolCompressionTransform symbolTT0;
		FSE_symbolCompressionTransform symbolTT1;
		FSE_symbolCompressionTransform symbolTT2;
		uint32_t nbBitsOut0;
		uint32_t nbBitsOut1;
		uint32_t nbBitsOut2;
		uint64_t sym0;
		uint64_t sym1;
		uint64_t sym2;
#ifdef PROFILE_SIZE
		uint64_t nbypass_total=0;
#endif

		rstr=3*iw;
		usize=(uint64_t)3*iw*ih;
		tables=(ETable*)malloc(sizeof(ETable[1<<LBITS]));
		packetscap=(int)(sizeof(LZPacket)*((usize+31)>>5));
		packets=(LZPacket*)malloc(packetscap);
		fse=(LZFSETables*)malloc(sizeof(*fse));
		if(!tables||!packets||!fse)
		{
			CRASH("Alloc error");
			return 0;
		}
		memset(tables, 0, sizeof(ETable[1<<LBITS]));
		npackets=0;
		memset(fse, 0, sizeof(*fse));
		imptr0=imptr=image;
		imend=image+usize;
		searchend=imend-sizeof(uint64_t);
		matchidx=0;
		matchlen=0;
		for(;;)
		{
			int tokena, tokenb, tokenc;

			if(imptr<imend)
			{
				matchidx=0;
				matchlen=0;
				lz_findmatch(tables, image, imptr, searchend, &matchidx, &matchlen);
				matchlen/=3;
				if(matchlen<(LZ_MINRUN+2)/3)
				{
					imptr+=3;
					continue;
				}

				//lazy parser
#if 0
				{
					uint32_t matchidx2, matchlen2;
					uint8_t *imptr2, *searchend;

					imptr2=imptr+3;
					searchend=imptr+matchlen*3;
				//	while(imptr2<searchend)
					{
						matchidx2=0;
						matchlen2=0;
						lz_findmatch(tables, image, imptr2, searchend, &matchidx2, &matchlen2);
						matchlen2/=3;
						if(matchlen2+1>matchlen)
						{
							imptr=imptr2;
							matchlen=matchlen2;
							matchidx=matchidx2;
					//		imptr2+=3;
					//		continue;
						}
					//	imptr2+=3;
					//	break;
					}
				}
#endif
			}
#if defined _MSC_VER &&1
			if(npackets>=packetscap/sizeof(LZPacket))
				CRASH("packets OoM");
#endif
			packet=packets+npackets++;
			packet->nliterals=(imptr-imptr0)/3;
			imptr0=imptr;
			if(imptr<imend)
			{
				packet->matchlen=matchlen;
				packet->offset=(imptr-image-matchidx)/3;
				imptr+=matchlen*3;
				imptr0=imptr;
			}
			else
			{
				packet->matchlen=0;
				packet->offset=0;
			}
#ifdef PROFILE_SIZE
			{
				int na=0, nb=0, nc=0;

				tokena=split_uint(packet->nliterals, &na, 0);
				tokenb=split_uint(packet->matchlen, &nb, 0);
				tokenc=split_uint(packet->offset, &nc, 0);
				nbypass_total+=na+nb+nc;
			}
#else
			tokena=split_uint(packet->nliterals, 0, 0);
			tokenb=split_uint(packet->matchlen, 0, 0);
			tokenc=split_uint(packet->offset, 0, 0);
#endif
#if defined _MSC_VER &&0
			if(tokena>>8||tokenb>>8||tokenc>>8)//
				CRASH("LZ token %d %d %d", tokena, tokenb, tokenc);
#endif
			++fse->hists[0][tokena];
			++fse->hists[1][tokenb];
			++fse->hists[2][tokenc];
			if(imptr>=imend)
				break;
		}
#if defined _MSC_VER &&0
		{
			int64_t sumlit=0, summatch=0;

			packet=packets;
			while(packet<packets+npackets)
			{
				sumlit+=packet->nliterals;
				summatch+=packet->matchlen;
				++packet;
			}
			if(sumlit+summatch!=(int64_t)iw*ih)
				CRASH("LZ error");
		}
#endif
		free(tables);

		imptr=image;
		packet=packets;
		imptr0=image+packet->nliterals*3;
		encbytes=1;
		for(ky=0;ky<ih;++ky)
		{
			int16_t *rows[]=
			{
				pixels+(XPAD*NCH*NROWS+(ky+0LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
			};
			int offset[2]={0}, N[3]={0}, W[3]={0}, pred[3]={0}, error[3]={0}, curr[3]={0}, ctx[3]={0};
			for(kx=0;kx<iw;++kx, imptr+=3)
			{
				yuv[0]=imptr[perm[0]];
				yuv[1]=imptr[perm[1]];
				yuv[2]=imptr[perm[2]];
				while(imptr>=imptr0)
				{
					if(encbytes)
						imptr0+=packet->matchlen*3;
					else
					{
						if(packet-packets>=npackets)
							CRASH("LZFSE predictor  out of LZ packets");
						++packet;
						imptr0+=packet->nliterals*3;
					}
					encbytes^=1;
				}
				{
					int epred, e0, negmask, abserr;

					offset[0]=rct.uc0*yuv[0]>>RCTBITS;
					offset[1]=(rct.vc0*yuv[0]+rct.vc1*yuv[1])>>RCTBITS;
					N[0]=rows[1][0+(0+0*NCH)*NROWS*NVAL];
					W[0]=rows[0][0+(0-1*NCH)*NROWS*NVAL];
					N[1]=rows[1][0+(1+0*NCH)*NROWS*NVAL];
					W[1]=rows[0][0+(1-1*NCH)*NROWS*NVAL];
					N[2]=rows[1][0+(2+0*NCH)*NROWS*NVAL];
					W[2]=rows[0][0+(2-1*NCH)*NROWS*NVAL];
					ctx[0]=0;
					ctx[1]=0;
					ctx[2]=0;
					imptr[0]=imptr[1]=imptr[2]=0;
					if(encbytes)
					{
						pred[0]=
							+rows[1][1+(0-1*NCH)*NROWS*NVAL]
							+rows[1][1+(0+0*NCH)*NROWS*NVAL]
							+rows[0][1+(0-1*NCH)*NROWS*NVAL]
							>
							+rows[1][2+(0-1*NCH)*NROWS*NVAL]
							+rows[1][2+(0+0*NCH)*NROWS*NVAL]
							+rows[0][2+(0-1*NCH)*NROWS*NVAL]
							? N[0] : W[0];
						pred[1]=
							+rows[1][1+(1-1*NCH)*NROWS*NVAL]
							+rows[1][1+(1+0*NCH)*NROWS*NVAL]
							+rows[0][1+(1-1*NCH)*NROWS*NVAL]
							>
							+rows[1][2+(1-1*NCH)*NROWS*NVAL]
							+rows[1][2+(1+0*NCH)*NROWS*NVAL]
							+rows[0][2+(1-1*NCH)*NROWS*NVAL]
							? N[1] : W[1];
						pred[2]=
							+rows[1][1+(2-1*NCH)*NROWS*NVAL]
							+rows[1][1+(2+0*NCH)*NROWS*NVAL]
							+rows[0][1+(2-1*NCH)*NROWS*NVAL]
							>
							+rows[1][2+(2-1*NCH)*NROWS*NVAL]
							+rows[1][2+(2+0*NCH)*NROWS*NVAL]
							+rows[0][2+(2-1*NCH)*NROWS*NVAL]
							? N[2] : W[2];
						pred[1]+=offset[0];
						pred[2]+=offset[1];
						CLAMP(pred[1], 0, 255);
						CLAMP(pred[2], 0, 255);
						for(kc=0;kc<3;++kc)
						{
							error[kc]=yuv[kc]-pred[kc];
							epred=pred[kc]>128?256-pred[kc]:pred[kc];
							e0=(int8_t)error[kc];
							negmask=error[kc]>>31;
							abserr=(error[kc]^negmask)-negmask;
							error[kc]=(uint32_t)error[kc]<<1^negmask;
							if(epred<abserr)
								error[kc]=epred+abserr;
							if(error[kc]==256)
								error[kc]=e0<<1^e0>>31;
							imptr[kc]=error[kc];
						}
						++fse->hists[2*0+(rows[1][3+(0+0*NCH)*NROWS*NVAL]|rows[0][3+(0-1*NCH)*NROWS*NVAL])+3][error[0]];
						++fse->hists[2*1+(rows[1][3+(1+0*NCH)*NROWS*NVAL]|rows[0][3+(1-1*NCH)*NROWS*NVAL])+3][error[1]];
						++fse->hists[2*2+(rows[1][3+(2+0*NCH)*NROWS*NVAL]|rows[0][3+(2-1*NCH)*NROWS*NVAL])+3][error[2]];
						ctx[0]=error[0]!=0;
						ctx[1]=error[1]!=0;
						ctx[2]=error[2]!=0;
					}
					curr[0]=yuv[0];
					curr[1]=yuv[1]-offset[0];
					curr[2]=yuv[2]-offset[1];
					rows[0][0+(0+0*NCH)*NROWS*NVAL]=curr[0];
					rows[0][1+(0+0*NCH)*NROWS*NVAL]=abs(curr[0]-W[0]);
					rows[0][2+(0+0*NCH)*NROWS*NVAL]=abs(curr[0]-N[0]);
					rows[0][3+(0+0*NCH)*NROWS*NVAL]=ctx[0];
					rows[0][0+(1+0*NCH)*NROWS*NVAL]=curr[1];
					rows[0][1+(1+0*NCH)*NROWS*NVAL]=abs(curr[1]-W[1]);
					rows[0][2+(1+0*NCH)*NROWS*NVAL]=abs(curr[1]-N[1]);
					rows[0][3+(1+0*NCH)*NROWS*NVAL]=ctx[1];
					rows[0][0+(2+0*NCH)*NROWS*NVAL]=curr[2];
					rows[0][1+(2+0*NCH)*NROWS*NVAL]=abs(curr[2]-W[2]);
					rows[0][2+(2+0*NCH)*NROWS*NVAL]=abs(curr[2]-N[2]);
					rows[0][3+(2+0*NCH)*NROWS*NVAL]=ctx[2];
					rows[0]+=NROWS*NVAL*NCH;
					rows[1]+=NROWS*NVAL*NCH;
					rows[2]+=NROWS*NVAL*NCH;
					rows[3]+=NROWS*NVAL*NCH;
				}
#if 0
				int offset=0;
				for(kc=0;kc<3;++kc)
				{
					int
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL],
					
						xNW	=rows[1][1-1*NCH*NROWS*NVAL],
						xN	=rows[1][1+0*NCH*NROWS*NVAL],
						xW	=rows[0][1-1*NCH*NROWS*NVAL],
					
						yNW	=rows[1][2-1*NCH*NROWS*NVAL],
						yN	=rows[1][2+0*NCH*NROWS*NVAL],
						yW	=rows[0][2-1*NCH*NROWS*NVAL],
						
						rN	=rows[1][3+0*NCH*NROWS*NVAL],
						rW	=rows[0][3-1*NCH*NROWS*NVAL];
					int ctx, pred, epred, error, curr, negmask, abserr, e0;
					
					offset=0;
					if(kc==1)offset=rct.uc0*yuv[0];
					if(kc==2)offset=rct.vc0*yuv[0]+rct.vc1*yuv[1];
					offset>>=RCTBITS;
					ctx=rN|rW;
					pred = xNW+xN+xW>yNW+yN+yW ? N : W;
					pred+=offset;
					CLAMP(pred, 0, 255);
					error=yuv[kc]-pred;
					epred=pred>128?256-pred:pred;
					e0=(int8_t)error;
					negmask=error>>31;
					abserr=(error^negmask)-negmask;
					error=error<<1^negmask;
					if(epred<abserr)
						error=epred+abserr;
					if(error==256)
						error=e0<<1^e0>>31;
					imptr[kc]=error;

					if(encbytes)
						++fse->hists[2*kc+ctx+3][error];

					curr=yuv[kc]-offset;
					rows[0][0]=curr;
					rows[0][1]=abs(curr-W);
					rows[0][2]=abs(curr-N);
					rows[0][3]=error!=0;
					rows[0]+=NROWS*NVAL;
					rows[1]+=NROWS*NVAL;
					rows[2]+=NROWS*NVAL;
					rows[3]+=NROWS*NVAL;
				}
#endif
			}
		}
#ifdef PROFILE_SIZE
		{
			double esizes[NHISTS]={0}, esizes2[2]={0};
			for(kc=0;kc<NHISTS;++kc)
			{
				int16_t *currnhist=fse->nhists[kc];
				int32_t *hcurr=fse->hists[kc];
				int ks;

				normalizehist((uint32_t*)hcurr, (uint16_t*)currnhist, PROBBITS);
				for(ks=0;ks<256;++ks)
				{
					int c=hcurr[ks];
					int f=currnhist[ks];
					if(c)
						esizes[kc]-=c*log2(f*(1./(1<<PROBBITS)));
				}
				printf("%5d  %12.2lf\n", kc, esizes[kc]/8.);
			}
			printf("\n");
			esizes2[0]=(esizes[0]+esizes[1]+esizes[2]+nbypass_total)/8.;
			esizes2[1]=(esizes[3]+esizes[4]+esizes[5]+esizes[6]+esizes[7]+esizes[8])/8.;
			printf("%12.2lf\n", esizes2[0]);
			printf("%12.2lf\n", esizes2[1]);
		}
#endif
		for(kc=0;kc<NHISTS;++kc)
		{
			enum
			{
				TABLESIZE=1<<PROBBITS,
				TABLEMASK=TABLESIZE-1,
				STEP=FSE_TABLESTEP(TABLESIZE),
				VMAX=255,
			};
			int16_t *currnhist=fse->nhists[kc];
			int32_t *hcurr=fse->hists[kc];
			int sum=0, nlevels=0, vmax=0, fmax=0;
			FSE_FUNCTION_TYPE *const tableSymbol=(FSE_FUNCTION_TYPE*)fse->wksp;
			uint16_t *tableU16=fse->stateTable[kc];
			FSE_symbolCompressionTransform *const symbolTT=fse->symbolTT[kc];
			uint32_t CDF[FSE_MAX_SYMBOL_VALUE+2];
			int highthreshold=TABLESIZE-1;
			
			if(PROBBITS>15||((size_t)1<<PROBBITS)*sizeof(FSE_FUNCTION_TYPE)>sizeof(fse->wksp))//PROBBITS<=15 required for threshold strategy to work
			{
				CRASH("tableLog too large");
				return 1;
			}
			for(int ks=0;ks<256;++ks)
			{
				int freq=hcurr[ks];
				sum+=freq;
				nlevels+=freq!=0;
				if(fmax<freq)
				{
					fmax=freq;
					vmax=ks;
				}
			}
			if(nlevels==1)//degenerate distribution
			{
				memset(currnhist, 0, sizeof(int16_t[256]));
				currnhist[vmax]=(1<<PROBBITS)-1;
				currnhist[(uint8_t)(vmax+1)]=1;
				degenmask|=1LL<<kc;
			}
			else if(sum<=256)//bypass distribution
			{
				for(int ks=0;ks<256;++ks)
					currnhist[ks]=1<<PROBBITS>>8;
				bypassmask|=1LL<<kc;
			}
			else//ordinary distribution
				normalizehist((uint32_t*)hcurr, (uint16_t*)currnhist, PROBBITS);

			{//Spread symbols
				uint32_t position=0, symbol;
				for(symbol=0;symbol<=VMAX;++symbol)
				{
					int freq=currnhist[symbol], nbOccurrences;
					for(nbOccurrences=0;nbOccurrences<freq;++nbOccurrences)
					{
						tableSymbol[position]=(FSE_FUNCTION_TYPE)symbol;
						position=(position+STEP)&TABLEMASK;
						while((int)position>highthreshold)
							position=(position+STEP)&TABLEMASK;//Low probability area
					}
				}
				if(position)//Must have initialized all positions
				{
					CRASH("Must have initialized all positions");
					return 1;
				}
			}

#ifdef __clang_analyzer__
			memset(fse->wksp, 0, sizeof(fse->wksp));//useless initialization, just to keep scan-build happy
#endif
			//For explanations on how to distribute symbol values over the table
			//http://fastcompression.blogspot.fr/2014/02/fse-distributing-symbol-values.html
			{//symbol start positions
				uint32_t u;
				CDF[0]=0;
				for(u=1;u<=VMAX+1;++u)
				{
					if(currnhist[u-1]==-1)//Low probability symbol
					{
						CDF[u]=CDF[u-1]+1;
						tableSymbol[highthreshold--]=(FSE_FUNCTION_TYPE)(u-1);
					}
					else
						CDF[u]=CDF[u-1]+currnhist[u-1];
				}
				CDF[VMAX+1]=TABLESIZE+1;
			}
			{//Build table
				uint32_t u;
				for(u=0;u<TABLESIZE;++u)
				{
					FSE_FUNCTION_TYPE s=tableSymbol[u];//note : static analyzer may not understand tableSymbol is properly initialized
					tableU16[CDF[s]++]=(uint16_t)(TABLESIZE+u);//TableU16 : sorted by symbol order; gives next state value
				}
			}
			{//Build Symbol Transformation Table
				uint32_t total=0, s;
				for(s=0;s<=VMAX;++s)
				{
					switch(currnhist[s])
					{
					case  0://filling nonetheless, for compatibility with FSE_getMaxNbBits()
						symbolTT[s].deltaNbBits=((PROBBITS+1)<<16)-(1<<PROBBITS);
						break;
					case -1:
					case  1:
						symbolTT[s].deltaNbBits=(PROBBITS<<16)-(1<<PROBBITS);
						symbolTT[s].deltaFindState=total-1;
						++total;
						break;
					default:
						{
							uint32_t const maxBitsOut=PROBBITS-(31-LZCNT32(currnhist[s]-1));
							uint32_t const minStatePlus=currnhist[s]<<maxBitsOut;
							symbolTT[s].deltaNbBits=(maxBitsOut<<16)-minStatePlus;
							symbolTT[s].deltaFindState=total-currnhist[s];
							total+=currnhist[s];
						}
						break;
					}
				}
			}
#if defined _DEBUG && 0		//debug : symbol costs
			printf("table statistics\n");
			{
				uint32_t symbol;
				for(symbol=0;symbol<=maxSymbolValue;++symbol)
				{
					printf("%3u: w=%3i,   maxBits=%u, fracBits=%.2f"
						, symbol
						, normalizedCounter[symbol]
						, FSE_getMaxNbBits(symbolTT, symbol)
						, (double)FSE_bitCost(symbolTT, tableLog, symbol, 8)/256
					);
				}
			}
#endif
			if(nlevels==1)//RLE
				currnhist[0]=vmax;
		}
		bitpacker_enc_init(&ec, stream, streamend);
		for(kc=0;kc<3;++kc)
			states[kc]=1<<PROBBITS;
		imptr=imend;
		packet=packets+npackets-1;
		imptr0=imptr-packet->matchlen*3;
		encbytes=0;
		kc=0;
		{
			int c0, c1, c2;

			kx=iw;
			while(imptr>image)
			{
				imptr-=3;
				--kx;
				if(kx<0)
					kx=iw-1;
				while(imptr<imptr0)
				{
					if(encbytes)
					{
						if(packet<=packets)
							CRASH("FSE encoder  out of LZ packets");
						--packet;
						imptr0-=packet->matchlen*3;
					}
					else
						imptr0-=packet->nliterals*3;
					encbytes^=1;
					if(!encbytes)
					{
					//	memset(imptr0, 0, imptr-imptr0);
						imptr=imptr0;
						kx=(imptr-image)/3%iw;
					}
				}
				if(encbytes)
				{
					c0=0;
					c1=0;
					c2=0;
					if(kx>=1)
					{
						c0|=imptr[0-3]!=0;
						c1|=imptr[1-3]!=0;
						c2|=imptr[2-3]!=0;
					}
					if(imptr>=image+rstr)
					{
						c0|=imptr[0-rstr]!=0;
						c1|=imptr[1-rstr]!=0;
						c2|=imptr[2-rstr]!=0;
					}
					symbolTT2=fse->symbolTT[c2+2*2+3][imptr[2]];
					symbolTT1=fse->symbolTT[c1+1*2+3][imptr[1]];
					symbolTT0=fse->symbolTT[c0+0*2+3][imptr[0]];
					nbBitsOut2=(states[2]+symbolTT2.deltaNbBits)>>16;
					nbBitsOut1=(states[1]+symbolTT1.deltaNbBits)>>16;
					nbBitsOut0=(states[0]+symbolTT0.deltaNbBits)>>16;
					sym2=states[2]&g_bitmask[nbBitsOut2];
					sym1=states[1]&g_bitmask[nbBitsOut1];
					sym0=states[0]&g_bitmask[nbBitsOut0];
					states[2]=(states[2]>>nbBitsOut2)+symbolTT2.deltaFindState;
					states[1]=(states[1]>>nbBitsOut1)+symbolTT1.deltaFindState;
					states[0]=(states[0]>>nbBitsOut0)+symbolTT0.deltaFindState;
					sym2|=sym1<<nbBitsOut2;
					sym2|=sym0<<(nbBitsOut2+nbBitsOut1);
					bitpacker_enc(&ec, (int64_t)nbBitsOut2+nbBitsOut1+nbBitsOut0, sym2);
#ifdef ANS_VAL
					{
						int c=c2<<2|c1<<1|c0;
						ansval_push(&c, sizeof(c), 1);
					}
#endif
					states[2]=fse->stateTable[c2+2*2+3][states[2]];
					states[1]=fse->stateTable[c1+1*2+3][states[1]];
					states[0]=fse->stateTable[c0+0*2+3][states[0]];
				}
			}
		}
		{
			int tokens[3]={0}, nbypass[3]={0};
			int64_t bypass[3]={0};

			packet=packets+npackets;
			while(packet>packets)
			{
				--packet;
				tokens[2]=split_uint(packet->offset	, nbypass+2, bypass+2);
				tokens[1]=split_uint(packet->matchlen	, nbypass+1, bypass+1);
				tokens[0]=split_uint(packet->nliterals	, nbypass+0, bypass+0);

				symbolTT2=fse->symbolTT[2][tokens[2]];
				symbolTT1=fse->symbolTT[1][tokens[1]];
				symbolTT0=fse->symbolTT[0][tokens[0]];
				nbBitsOut2=(states[2]+symbolTT2.deltaNbBits)>>16;
				nbBitsOut1=(states[1]+symbolTT1.deltaNbBits)>>16;
				nbBitsOut0=(states[0]+symbolTT0.deltaNbBits)>>16;
#if defined _MSC_VER &&0
				if(nbBitsOut2>PROBBITS||nbBitsOut1>PROBBITS||nbBitsOut0>PROBBITS)
					CRASH("FSE nbBitsOut %d %d %d", nbBitsOut2, nbBitsOut1, nbBitsOut0);
#endif
				sym2=states[2]&g_bitmask[nbBitsOut2];
				sym1=states[1]&g_bitmask[nbBitsOut1];
				sym0=states[0]&g_bitmask[nbBitsOut0];
				states[2]=(states[2]>>nbBitsOut2)+symbolTT2.deltaFindState;
				states[1]=(states[1]>>nbBitsOut1)+symbolTT1.deltaFindState;
				states[0]=(states[0]>>nbBitsOut0)+symbolTT0.deltaFindState;
				if(nbypass[2]>0)
					bitpacker_enc(&ec, nbypass[2], bypass[2]);
				if(nbypass[1]>0)
					bitpacker_enc(&ec, nbypass[1], bypass[1]);
				if(nbypass[0]>0)
					bitpacker_enc(&ec, nbypass[0], bypass[0]);
				sym2|=sym1<<nbBitsOut2;
				sym2|=sym0<<(nbBitsOut2+nbBitsOut1);
				bitpacker_enc(&ec, (int64_t)nbBitsOut2+nbBitsOut1+nbBitsOut0, sym2);
				states[2]=fse->stateTable[2][states[2]];
				states[1]=fse->stateTable[1][states[1]];
				states[0]=fse->stateTable[0][states[0]];
			}
		}
		for(kc=3-1;kc>=0;--kc)
			bitpacker_enc(&ec, PROBBITS, states[kc]&((1<<PROBBITS)-1));
		bitpacker_enc(&ec, 24, npackets);
		for(kc=NHISTS-1;kc>=0;--kc)
		{
			int16_t *currnhist=fse->nhists[kc];
			if(degenmask>>kc&1)//degenerate distribution
				bitpacker_enc(&ec, 8, (uint8_t)currnhist[0]);
			else if(!(bypassmask>>kc&1))//ordinary distribution
				enc_packhist(&ec, (int*)currnhist, bypassmask, kc, 0, PROBBITS);
		}
		bitpacker_enc(&ec, (int64_t)NHISTS, bypassmask);
		bitpacker_enc(&ec, (int64_t)NHISTS, degenmask);
		streamend=bitpacker_flush(&ec);
		csize=streamend-stream;
		free(packets);
		free(fse);
	}
	else//dec
	{
		LZFSDTables *fsd=0;
		int decbytes=0;
		uint32_t DInfo[3]={0};
		int syms[3]={0}, nbits[3]={0}, c1=0, c2=0;
		int preds[3]={0}, curr[3]={0};

		fsd=(LZFSDTables*)malloc(sizeof(*fsd));
		if(!fsd)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(fsd, 0, sizeof(*fsd));
		bitpacker_dec_init(&ec, stream, streamend);
		degenmask=(int)bitpacker_dec(&ec, (int64_t)NHISTS);
		bypassmask=(int)bitpacker_dec(&ec, (int64_t)NHISTS);
		for(kc=0;kc<NHISTS;++kc)
		{
			int16_t nhist[256]={0};
			if(degenmask>>kc&1)//degenerate distribution
			{
				int sym=(int)bitpacker_dec(&ec, 8);
				memset(nhist, 0, sizeof(nhist));
				nhist[sym]=(1<<PROBBITS)-1;
				nhist[(uint8_t)(sym+1)]=1;
			}
			else if(bypassmask>>kc&1)//bypass distribution
			{
				for(int ks=0;ks<256;++ks)
					nhist[ks]=1<<PROBBITS>>8;
			}
			else//ordinary distribution
				dec_unpackhist(&ec, bypassmask, kc, (uint16_t*)nhist, PROBBITS);

			{
				enum
				{
					TABLESIZE=1<<PROBBITS,
					TABLEMASK=TABLESIZE-1,
					STEP=FSE_TABLESTEP(TABLESIZE),
					VMAX=255,
					NLEVELS=VMAX+1,
				//	LARGELIMIT=1<<(PROBBITS-1),
				};
				FSE_DECODE_TYPE *const tableDecode=(FSE_DECODE_TYPE*)fsd->decodeTable[kc];
				uint16_t symbolNext[FSE_MAX_SYMBOL_VALUE+1];
				uint32_t highThreshold=TABLESIZE-1;

				if(VMAX>FSE_MAX_SYMBOL_VALUE)//Sanity Checks
				{
					CRASH("maxSymbolValue too large");
					return 1;
				}
				if(PROBBITS>FSE_MAX_TABLELOG)
				{
					CRASH("tableLog too large");
					return 1;
				}
				{//Init, lay down lowprob symbols
					uint32_t s;
					for(s=0;s<NLEVELS;++s)
					{
						if(nhist[s]==-1)
						{
							tableDecode[highThreshold--].symbol=(FSE_FUNCTION_TYPE)s;
							symbolNext[s]=1;
						}
						else
							symbolNext[s]=nhist[s];
					}
				}
				{//Spread symbols
					uint32_t s, position=0;
					for(s=0;s<NLEVELS;++s)
					{
						int i;
						for(i=0;i<nhist[s];++i)
						{
							tableDecode[position].symbol=(FSE_FUNCTION_TYPE)s;
							position=(position+STEP)&TABLEMASK;
							while(position>highThreshold)
								position=(position+STEP)&TABLEMASK;//lowprob area
						}
					}
					if(position!=0)//position must reach all cells once, otherwise normalizedCounter is incorrect
					{
						CRASH("position must reach all cells once, otherwise normalizedCounter is incorrect");
						return 1;
					}
				}
				{//Build Decoding table
					uint32_t u;
					for(u=0;u<TABLESIZE;++u)
					{
						FSE_FUNCTION_TYPE const symbol=(FSE_FUNCTION_TYPE)tableDecode[u].symbol;
						uint32_t const nextState=symbolNext[symbol]++;
						tableDecode[u].nbBits=(uint8_t)(PROBBITS-31+LZCNT32(nextState));
						tableDecode[u].newState=(uint16_t)((nextState<<tableDecode[u].nbBits)-TABLESIZE);
					}
				}
			}
		}
		npackets=(int)bitpacker_dec(&ec, 24);
		packets=(LZPacket*)malloc(npackets*sizeof(LZPacket));
		if(!packets)
		{
			CRASH("Alloc error");
			return 1;
		}
		for(kc=0;kc<3;++kc)
			states[kc]=(uint32_t)bitpacker_dec(&ec, PROBBITS);
		packet=packets;
		{
			uint64_t bypass[3]={0}, val[3]={0};
			uint64_t sym=0;

			while(packet<packets+npackets)
			{
				DInfo[0]=*(uint32_t*)&fsd->decodeTable[0][states[0]];
				DInfo[1]=*(uint32_t*)&fsd->decodeTable[1][states[1]];
				DInfo[2]=*(uint32_t*)&fsd->decodeTable[2][states[2]];
				syms[0]=(uint8_t)(DInfo[0]>>16);
				syms[1]=(uint8_t)(DInfo[1]>>16);
				syms[2]=(uint8_t)(DInfo[2]>>16);
				nbits[0]=(uint8_t)(DInfo[0]>>24);
				nbits[1]=(uint8_t)(DInfo[1]>>24);
				nbits[2]=(uint8_t)(DInfo[2]>>24);
				states[0]=(uint16_t)DInfo[0];
				states[1]=(uint16_t)DInfo[1];
				states[2]=(uint16_t)DInfo[2];
				c1=nbits[2];
				c2=nbits[2]+nbits[1];
				sym=bitpacker_dec(&ec, (int64_t)nbits[0]+c2);
				states[0]+=(uint32_t)(sym>>c2&g_bitmask[nbits[0]]);
				states[1]+=(uint32_t)(sym>>c1&g_bitmask[nbits[1]]);
				states[2]+=(uint32_t)(sym>> 0&g_bitmask[nbits[2]]);

				for(kc=0;kc<3;++kc)
				{
					val[kc]=syms[kc];
					if(syms[kc]>=1<<CONFIG_EXP)
					{
						nbits[kc]=(syms[kc]>>(CONFIG_MSB+CONFIG_LSB)) - ((1<<(CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB)))-(CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB)));
						bypass[kc]=bitpacker_dec(&ec, nbits[kc]);
						val[kc]=
							+((((uint64_t)1<<(CONFIG_MSB+CONFIG_LSB))+(syms[kc]&(((uint64_t)1<<CONFIG_MSB)-1)<<CONFIG_LSB))<<nbits[kc])//7 instructions
							+(bypass[kc]<<CONFIG_LSB)
							+(syms[kc]&((1<<CONFIG_LSB)-1))
						;
					}
				}
				packet->nliterals=val[0];
				packet->matchlen=val[1];
				packet->offset=val[2];
				++packet;
			}
		}
		imptr=image;
		packet=packets;
		decbytes=1;
		imptr0=imptr+packet->nliterals*3;
		for(ky=0;ky<ih;++ky)
		{
			int16_t *rows[]=
			{
				pixels+(XPAD*NCH*NROWS+(ky+0LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
			};
			for(kx=0;kx<iw;++kx, imptr+=3)
			{
				int offsets[2]={0};
				int
					N0	=rows[1][0+(0+0*NCH)*NROWS*NVAL],
					W0	=rows[0][0+(0-1*NCH)*NROWS*NVAL],
					N1	=rows[1][0+(1+0*NCH)*NROWS*NVAL],
					W1	=rows[0][0+(1-1*NCH)*NROWS*NVAL],
					N2	=rows[1][0+(2+0*NCH)*NROWS*NVAL],
					W2	=rows[0][0+(2-1*NCH)*NROWS*NVAL];
				uint64_t sym;
				int negmask, e2, epred, error;

				while(imptr>=imptr0)
				{
					if(decbytes)
						imptr0+=packet->matchlen*3;
					else
					{
						if(packet-packets>=npackets)
							CRASH("LZFSE decoder  out of LZ packets");
						++packet;
						imptr0+=packet->nliterals*3;
					}
					decbytes^=1;
				}
				if(!decbytes)
				{
					yuv[0]=imptr[perm[0]-packet->offset*3];
					yuv[1]=imptr[perm[1]-packet->offset*3];
					yuv[2]=imptr[perm[2]-packet->offset*3];
					offsets[0]=rct.uc0*yuv[0]>>RCTBITS;
					offsets[1]=(rct.vc0*yuv[0]+rct.vc1*yuv[1])>>RCTBITS;
					syms[0]=0;
					syms[1]=0;
					syms[2]=0;
				}
				else
				{
					int
						xNW0	=rows[1][1+(0-1*NCH)*NROWS*NVAL],
						xN0	=rows[1][1+(0+0*NCH)*NROWS*NVAL],
						xW0	=rows[0][1+(0-1*NCH)*NROWS*NVAL],
						yNW0	=rows[1][2+(0-1*NCH)*NROWS*NVAL],
						yN0	=rows[1][2+(0+0*NCH)*NROWS*NVAL],
						yW0	=rows[0][2+(0-1*NCH)*NROWS*NVAL],
						rN0	=rows[1][3+(0+0*NCH)*NROWS*NVAL],
						rW0	=rows[0][3+(0-1*NCH)*NROWS*NVAL],
						xNW1	=rows[1][1+(1-1*NCH)*NROWS*NVAL],
						xN1	=rows[1][1+(1+0*NCH)*NROWS*NVAL],
						xW1	=rows[0][1+(1-1*NCH)*NROWS*NVAL],
						yNW1	=rows[1][2+(1-1*NCH)*NROWS*NVAL],
						yN1	=rows[1][2+(1+0*NCH)*NROWS*NVAL],
						yW1	=rows[0][2+(1-1*NCH)*NROWS*NVAL],
						rN1	=rows[1][3+(1+0*NCH)*NROWS*NVAL],
						rW1	=rows[0][3+(1-1*NCH)*NROWS*NVAL],
						xNW2	=rows[1][1+(2-1*NCH)*NROWS*NVAL],
						xN2	=rows[1][1+(2+0*NCH)*NROWS*NVAL],
						xW2	=rows[0][1+(2-1*NCH)*NROWS*NVAL],
						yNW2	=rows[1][2+(2-1*NCH)*NROWS*NVAL],
						yN2	=rows[1][2+(2+0*NCH)*NROWS*NVAL],
						yW2	=rows[0][2+(2-1*NCH)*NROWS*NVAL],
						rN2	=rows[1][3+(2+0*NCH)*NROWS*NVAL],
						rW2	=rows[0][3+(2-1*NCH)*NROWS*NVAL];
					preds[0] = xNW0+xN0+xW0>yNW0+yN0+yW0 ? N0 : W0;
					preds[1] = xNW1+xN1+xW1>yNW1+yN1+yW1 ? N1 : W1;
					preds[2] = xNW2+xN2+xW2>yNW2+yN2+yW2 ? N2 : W2;
					DInfo[0]=*(uint32_t*)&fsd->decodeTable[(rN0|rW0)+0*2+3][states[0]];
					DInfo[1]=*(uint32_t*)&fsd->decodeTable[(rN1|rW1)+1*2+3][states[1]];
					DInfo[2]=*(uint32_t*)&fsd->decodeTable[(rN2|rW2)+2*2+3][states[2]];
					syms[0]=(uint8_t)(DInfo[0]>>16);
					syms[1]=(uint8_t)(DInfo[1]>>16);
					syms[2]=(uint8_t)(DInfo[2]>>16);
					nbits[0]=(uint8_t)(DInfo[0]>>24);
					nbits[1]=(uint8_t)(DInfo[1]>>24);
					nbits[2]=(uint8_t)(DInfo[2]>>24);
					states[0]=(uint16_t)DInfo[0];
					states[1]=(uint16_t)DInfo[1];
					states[2]=(uint16_t)DInfo[2];
					c1=nbits[2];
					c2=nbits[2]+nbits[1];
#ifdef ANS_VAL
					{
						int c=(rN2|rW2)<<2|(rN1|rW1)<<1|(rN0|rW0);
						ansval_check(&c, sizeof(c), 1);
					}
#endif
					sym=bitpacker_dec(&ec, (int64_t)nbits[0]+c2);
					states[0]+=(uint32_t)(sym>>c2&g_bitmask[nbits[0]]);
					states[1]+=(uint32_t)(sym>>c1&g_bitmask[nbits[1]]);
					states[2]+=(uint32_t)(sym>> 0&g_bitmask[nbits[2]]);
					
					epred=preds[0]-128;
					negmask=epred>>31;
					epred=128-((epred^negmask)-negmask);
					e2=epred-syms[0];
					error=syms[0]>>1^syms[0]<<31>>31;
					yuv[0]=(uint8_t)(error+preds[0]);
					e2=(e2^negmask)-negmask;
					if(epred<<1<syms[0])
						error=e2;
					if(2*preds[0]+syms[0]!=2<<8)
						yuv[0]=error+(int32_t)preds[0];

					offsets[0]=rct.uc0*yuv[0]>>RCTBITS;
					preds[1]+=offsets[0];
					CLAMP(preds[1], 0, 255);
					epred=preds[1]-128;
					negmask=epred>>31;
					epred=128-((epred^negmask)-negmask);
					e2=epred-syms[1];
					error=syms[1]>>1^syms[1]<<31>>31;
					yuv[1]=(uint8_t)(error+preds[1]);
					e2=(e2^negmask)-negmask;
					if(epred<<1<syms[1])
						error=e2;
					if(2*preds[1]+syms[1]!=2<<8)
						yuv[1]=error+(int32_t)preds[1];

					offsets[1]=(rct.vc0*yuv[0]+rct.vc1*yuv[1])>>RCTBITS;
					preds[2]+=offsets[1];
					CLAMP(preds[2], 0, 255);
					epred=preds[2]-128;
					negmask=epred>>31;
					epred=128-((epred^negmask)-negmask);
					e2=epred-syms[2];
					error=syms[2]>>1^syms[2]<<31>>31;
					yuv[2]=(uint8_t)(error+preds[2]);
					e2=(e2^negmask)-negmask;
					if(epred<<1<syms[2])
						error=e2;
					if(2*preds[2]+syms[2]!=2<<8)
						yuv[2]=error+(int32_t)preds[2];

					syms[0]=syms[0]!=0;
					syms[1]=syms[1]!=0;
					syms[2]=syms[2]!=0;
				}
				imptr[perm[0]]=yuv[0];
				imptr[perm[1]]=yuv[1];
				imptr[perm[2]]=yuv[2];
#if defined ENABLE_GUIDE &&1
				if(memcmp(imptr, g_image+(imptr-image), 3))
				{
					printf("\n\nguide  X %5d  Y %5d  %02X%02X%02X  != OG %06X\n\n"
						, kx
						, ky
						, imptr[2]
						, imptr[1]
						, imptr[0]
						, *(uint32_t*)(g_image+(imptr-image))<<8>>8
					);
					CRASH("");
				}
#endif
				curr[0]=yuv[0];
				curr[1]=yuv[1]-offsets[0];
				curr[2]=yuv[2]-offsets[1];
				rows[0][0+(0+0*NCH)*NROWS*NVAL]=curr[0];
				rows[0][1+(0+0*NCH)*NROWS*NVAL]=abs(curr[0]-W0);
				rows[0][2+(0+0*NCH)*NROWS*NVAL]=abs(curr[0]-N0);
				rows[0][3+(0+0*NCH)*NROWS*NVAL]=syms[0];
				
				rows[0][0+(1+0*NCH)*NROWS*NVAL]=curr[1];
				rows[0][1+(1+0*NCH)*NROWS*NVAL]=abs(curr[1]-W1);
				rows[0][2+(1+0*NCH)*NROWS*NVAL]=abs(curr[1]-N1);
				rows[0][3+(1+0*NCH)*NROWS*NVAL]=syms[1];
				
				rows[0][0+(2+0*NCH)*NROWS*NVAL]=curr[2];
				rows[0][1+(2+0*NCH)*NROWS*NVAL]=abs(curr[2]-W2);
				rows[0][2+(2+0*NCH)*NROWS*NVAL]=abs(curr[2]-N2);
				rows[0][3+(2+0*NCH)*NROWS*NVAL]=syms[2];
				rows[0]+=NCH*NROWS*NVAL;
				rows[1]+=NCH*NROWS*NVAL;
				rows[2]+=NCH*NROWS*NVAL;
				rows[3]+=NCH*NROWS*NVAL;
			}
		}
		free(packets);
		free(fsd);
	}
	free(pixels);

	//save
	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing\n", dstfn);
		free(image);
		free(stream);
		return 1;
	}
	if(fwd)
	{
		csize=0;
		csize+=fwrite(&tag, 1, 2, fdst);
		csize+=fwrite(&iw, 1, 3, fdst);
		csize+=fwrite(&ih, 1, 3, fdst);
		csize+=fwrite(&rct, 1, sizeof(rct), fdst);
		csize+=fwrite(stream, 1, streamend-stream, fdst);
	}
	else
	{
		int64_t usize2=0;

		usize2+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		usize2+=fwrite(image, 1, usize, fdst);
		usize=usize2;
	}
	fclose(fdst);
	free(image);
	free(stream);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
#ifdef TRAIN_HIST
		for(int k=0;k<THSIZE;++k)
			printf("%12d %12d %12d\n", thist[0][k], thist[1][k], thist[2][k]);
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
#endif
	(void)&time_sec;
	(void)&print_rct;
	return 0;
}
