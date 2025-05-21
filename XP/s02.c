#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
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
#include<intrin.h>
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>//QueryPerformanceCounter
#elif defined __GNUC__
#include<x86intrin.h>
#include<time.h>
#endif


#ifdef _MSC_VER
	#define LOUD
#ifdef _DEBUG
	#define ESTIMATE_SIZE
	#define ENABLE_GUIDE
//	#define ANS_VAL
//	#define PRINTBITS
//	#define PRINTGR
#endif
#endif

	#define USE_L1
//	#define USE_CASCADE


//3*17+3=54 contexts
#define GRBITS 3
#define NCTX 18		//18*3+3 = 57 total

#define PROBBITS 12	//12 bit max	James Bonfield's CDF2sym: {freq<<20 | bias<<8 | sym}
#define RANS_STATE_BITS 31
#define RANS_RENORM_BITS 16

#if 1
#define L1NPREDS 8
#define PREDLIST\
	PRED(100000, N)\
	PRED(100000, W)\
	PRED(100000, 3*(W-WW)+WWW)\
	PRED(100000, 3*(N-NN)+NNN)\
	PRED(100000, N+W-NW)\
	PRED(100000, W+NE-N)\
	PRED(100000, N+NE-NNE)\
	PRED(100000, (WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)>>2)
#endif
#if 0
#define L1NPREDS 13
#define PREDLIST\
	PRED(100000, N+W-NW)\
	PRED(200000, N+W-NW)\
	PRED(100000, N)\
	PRED(100000, W)\
	PRED(100000, 3*(N-NN)+NNN)\
	PRED(100000, 3*(W-WW)+WWW)\
	PRED(100000, W+NE-N)\
	PRED(100000, N+NE-NNE)\
	PRED(100000, W+((NEEE+NEEEEE-N-W)>>3))\
	PRED( 50000, W+NW-NWW)\
	PRED( 50000, N+NW-NNW)\
	PRED( 50000, NE+NEE-NNEEE)\
	PRED( 50000, (WWWWW+WW-W+NNN+N+NEEEEE)>>2)
#endif
#ifdef USE_CASCADE
#define L1NPREDS2 12
#define PREDLIST2\
	PRED(eN+eW-eNW)\
	PRED(eN)\
	PRED(eW)\
	PRED(3*(eN-eNN)+eNNN)\
	PRED(3*(eW-eWW)+eWWW)\
	PRED(eW+eNE-eN)\
	PRED(eN+eNE-eNNE)\
	PRED(eW+((eNEEE+eNEEEEE-eN-eW)>>3))\
	PRED(eW+eNW-eNWW)\
	PRED(eN+eNW-eNNW)\
	PRED(eNE+eNEE-eNNEEE)\
	PRED((eWWWWW+eWW-eW+eNNN+eN+eNEEEEE)>>2)
#endif

#ifdef _MSC_VER
#define AWM_INLINE __forceinline static
#else
#define AWM_INLINE __attribute__((always_inline)) inline static
#endif
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
static void memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
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
#define FILLMEM(PTR, DATA, ASIZE, ESIZE)\
	do\
	{\
		*(PTR)=(DATA);\
		memfill((PTR)+1, PTR, (ASIZE)-(ESIZE), ESIZE);\
	}while(0)
static double time_sec(void)
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
AWM_INLINE uint32_t floor_log2(uint32_t n)//n >= 1
{
#ifdef M_LZCNT
	return 31-_lzcnt_u32(n);
#else
	if(n>0x7FFFFFFF)
		return 31;
	{
		__m128i t0=_mm_castpd_si128(_mm_cvtsi32_sd(_mm_setzero_pd(), n));
		t0=_mm_srli_epi64(t0, 52);
		return _mm_cvtsi128_si32(t0)-1023;
	}
#if 0
	float x=(float)n;
	size_t addr=(size_t)&x;
	uint32_t bits=*(uint32_t*)addr;
	return (bits>>23)-127;
#endif
#endif
}
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
	exit(1);
}
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT,##__VA_ARGS__)
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static void guide_save(unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		CRASH("");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		CRASH("");
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif
static unsigned char* load_file(const char *fn, ptrdiff_t *ret_size)
{
	struct stat info={0};
	int error;
	ptrdiff_t size, nread;
	FILE *fsrc;
	unsigned char *buf;

	error=stat(fn, &info);
	if(error)
	{
		printf("Cannot stat \"%s\"\n", fn);
		return 0;
	}
	size=info.st_size;
	fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		printf("Cannot open \"%s\"\n", fn);
		return 0;
	}
	buf=(unsigned char*)malloc(size+16);
	if(!buf)
	{
		printf("Alloc error\n");
		exit(1);
		return 0;
	}
	nread=fread(buf, 1, size, fsrc);
	if(nread!=size)
		printf("Error reading \"%s\"\n", fn);
	fclose(fsrc);
	*ret_size=size;
	return buf;
}
typedef enum _OCHType
{
	OCH_R,
	OCH_G,
	OCH_B,
	OCH_RG,
	OCH_GB,
	OCH_BR,
//	OCH_R2,
//	OCH_G2,
//	OCH_B2,

	OCH_COUNT,

	OCH_RB=OCH_BR,
	OCH_GR=OCH_RG,
	OCH_BG=OCH_GB,

} OCHType;
static const unsigned char rct_indices[][8]=
{//	output channels			permutation	helper index
	{OCH_R,	OCH_G,	OCH_B,		0, 1, 2,	3, 3},// 0
	{OCH_R,	OCH_G,	OCH_BG,		0, 1, 2,	3, 1},// 1
	{OCH_R,	OCH_G,	OCH_BR,		0, 1, 2,	3, 0},// 2
	{OCH_R,	OCH_GR,	OCH_BR,		0, 1, 2,	0, 0},// 3
	{OCH_R,	OCH_GR,	OCH_BG,		0, 1, 2,	0, 1},// 4
	{OCH_R,	OCH_BR,	OCH_GB,		0, 2, 1,	0, 1},// 5
	{OCH_G,	OCH_B,	OCH_RG,		1, 2, 0,	3, 0},// 6
	{OCH_G,	OCH_B,	OCH_RB,		1, 2, 0,	3, 1},// 7
	{OCH_G,	OCH_BG,	OCH_RG,		1, 2, 0,	0, 0},// 8
	{OCH_G,	OCH_BG,	OCH_RB,		1, 2, 0,	0, 1},// 9
	{OCH_G,	OCH_RG,	OCH_BR,		1, 0, 2,	0, 1},//10
	{OCH_B,	OCH_R,	OCH_GR,		2, 0, 1,	3, 1},//11
	{OCH_B,	OCH_R,	OCH_GB,		2, 0, 1,	3, 0},//12
	{OCH_B,	OCH_RB,	OCH_GB,		2, 0, 1,	0, 0},//13
	{OCH_B,	OCH_RB,	OCH_GR,		2, 0, 1,	0, 1},//14
	{OCH_B,	OCH_GB,	OCH_RG,		2, 1, 0,	0, 1},//15
};

#ifdef ANS_VAL
#define ANS_VAL_HISTSIZE 128
typedef struct _ANSVALHeader
{
	unsigned short esize, count;
	unsigned idx;
	struct _ANSVALHeader *above, *below;
	unsigned char data[];
} ANSVALNode;
static ANSVALNode *debugstack=0;
static int ansvalidx=0;
static void ansval_push(const void *data, int esize, int count)
{
	int size=sizeof(ANSVALNode)+esize*count;
	ANSVALNode *node=(ANSVALNode*)malloc(size);
	if(!node)
	{
		printf("Alloc error\n");
		exit(1);
	}
	memset(node, 0, size);
	node->esize=esize;
	node->count=count;
	node->idx=ansvalidx++;
	node->above=0;
	node->below=debugstack;
	if(debugstack)
		debugstack->above=node;
	memcpy(node->data, data, size-sizeof(ANSVALNode));
	debugstack=node;
}
static void ansval_printr(const void *data, int esize, int count, const void *xdata)//print elements in reverse because little-endian
{
	const unsigned char *p=(const unsigned char*)data, *p2=(const unsigned char*)xdata;
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
	//for(int k=size-1;k>=0;--k)
	//{
	//	if((k&3)==3)
	//		printf(" ");
	//	printf("%02X", ((unsigned char*)data)[k]);
	//}
	//printf("\n");
}
static void* ansval_ptrguard(const void *start, const void *end, const void *ptr, ptrdiff_t nbytes)
{
	size_t istart=(size_t)start, iend=(size_t)end;
	ptrdiff_t size=iend-istart;
	size_t ip1=(size_t)ptr, ip2=ip1+nbytes;
	int problems[]=
	{
		size<0,
		(size_t)(ip1-istart)>=(size_t)size,
		(size_t)(ip2-istart)>=(size_t)size,
	};
	if(problems[0]||problems[1]||problems[2])
	{
		printf("\nOOB\n");
		printf("  inc     %+16td bytes\n", nbytes);
		printf("  start   %016zd  %16d\n", istart, 0);
		if(nbytes<0)
		{
			printf("  after   %016zd  %16td%s\n", ip2, ip2-istart, problems[2]?"  <-":"");
			printf("  before  %016zd  %16td%s\n", ip1, ip1-istart, problems[1]?"  <-":"");
		}
		else
		{
			printf("  before  %016zd  %16td%s\n", ip1, ip1-istart, problems[1]?"  <-":"");
			printf("  after   %016zd  %16td%s\n", ip2, ip2-istart, problems[2]?"  <-":"");
		}
		printf("  end     %016zd  %16td%s\n", iend, size, problems[0]?"  <-":"");
		CRASH("\n");
		return 0;
	}
	return (void*)(nbytes<0?ip2:ip1);
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
#if 0
	static const unsigned char *endptr=0;
	static int totalcount=0, popcount=0;
	int firstpop=!endptr;
	if(firstpop)
		endptr=debugstack->data+debugstack->count;
	const ANSVALHeader *hiheader=(const ANSVALHeader*)(debugstack->data+debugstack->count)-1;
	debugstack->count-=sizeof(ANSVALHeader);
	const unsigned char *data0=debugstack->data+debugstack->count-hiheader->count*hiheader->esize;
	debugstack->count-=hiheader->count*hiheader->esize;
	const ANSVALHeader *loheader=(const ANSVALHeader*)(debugstack->data+debugstack->count)-1;
	debugstack->count-=sizeof(ANSVALHeader);
	if(firstpop)
		totalcount=hiheader->idx+1;
	if(memcmp(hiheader, loheader, sizeof(*hiheader)))
	{
		printf("\n");
		printf("Validation Header Mismatch  idx,esize,count: loheader %d,%d,%d != hiheader %d,%d,%d\n",
			loheader->idx, loheader->esize, loheader->count,
			hiheader->idx, hiheader->esize, hiheader->count
		);
		CRASH("\n");
	}
	if(esize!=loheader->esize||count!=loheader->count||memcmp(data, data0, count*esize))
	{
		printf("\n");
		printf("Validation Error    pop #%d / total %d,  remaining %d,  using %8.2lf/%8.2lf MB\n",
			popcount,
			totalcount,
			totalcount-popcount-1,
			(double)debugstack->count/(1024*1024),
			(double)debugstack->cap/(1024*1024)
		);
		printf("\n");

		const unsigned char *verptr=debugstack->data+debugstack->count+loheader->count*loheader->esize+sizeof(ANSVALHeader[2]);
		const unsigned char *ptrstack[ANS_VAL_HISTSIZE]={0};
		const ANSVALHeader *loheader2=0, *hiheader2=0;
		const unsigned char *verdata=0, *unverdata=0;
		int nptrs=0;
		printf("Verified pops:\n");
		for(int k=0;k<ANS_VAL_HISTSIZE;++k)
		{
			if(verptr>=endptr)
				break;

			loheader2=(const ANSVALHeader*)ansval_ptrguard(debugstack->data, endptr, verptr, +sizeof(ANSVALHeader));
			verptr+=sizeof(ANSVALHeader);
			verdata=(const unsigned char*)ansval_ptrguard(debugstack->data, endptr, verptr, +loheader2->count*loheader2->esize);
			verptr+=loheader2->count*loheader2->esize;
			hiheader2=(const ANSVALHeader*)verptr;
			verptr+=sizeof(ANSVALHeader);

			ptrstack[nptrs++]=(const unsigned char*)loheader2;
			(void)verdata;
			(void)hiheader2;
		}
		for(int k=nptrs-1;k>=0;--k)
		{
			const unsigned char *ptr=ptrstack[k];
			loheader2=(const ANSVALHeader*)ptr;
			ptr+=sizeof(ANSVALHeader);
			verdata=ptr;
			ptr+=loheader2->count*loheader2->esize;
			hiheader2=(const ANSVALHeader*)ptr;
			ptr+=sizeof(ANSVALHeader);

			printf("  [%7d] %7d B    ", loheader2->idx, loheader2->count*loheader2->esize);
			ansval_printr(verdata, loheader2->esize, loheader2->count, 0);
			(void)hiheader2;
		}
		if(!nptrs)
			printf("  No data\n");
		printf("\n");

		printf("The error:\n");
		printf("  [%7d] Original %7d B    ", loheader->idx, loheader->count*loheader->esize);
		ansval_printr(data0, loheader->esize, loheader->count, 0);
		printf("  [%7d] Corrupt  %7d B    ", loheader->idx, count*esize);
		ansval_printr(data, esize, count, 0);
		printf("  [%7d] XOR      %7s      ", loheader->idx, "");
		ansval_printr(data0, esize, count, data);
		printf("\n");
		
		const unsigned char *unverptr=debugstack->data+debugstack->count;
		printf("Remaining pops:\n");
		for(int k=0;k<ANS_VAL_HISTSIZE;++k)
		{
			if(unverptr<=debugstack->data)
			{
				printf("  No data\n");
				break;
			}

			unverptr=(const unsigned char*)ansval_ptrguard(debugstack->data, endptr, unverptr, -(ptrdiff_t)sizeof(ANSVALHeader));
			hiheader2=(const ANSVALHeader*)unverptr;
			unverptr=(const unsigned char*)ansval_ptrguard(debugstack->data, endptr, unverptr, -(ptrdiff_t)hiheader2->count*hiheader2->esize);
			unverdata=unverptr;
			unverptr=(const unsigned char*)ansval_ptrguard(debugstack->data, endptr, unverptr, -(ptrdiff_t)sizeof(ANSVALHeader));
			loheader2=(const ANSVALHeader*)unverptr;
			
			if(memcmp(hiheader2, loheader2, sizeof(*hiheader2)))
			{
				printf("  Header Mismatch:  (idx,esize,count)\n");
				printf("    hiheader  0x%08X, 0x%08X, 0x%08X    %d, %d, %d\n",
					hiheader2->idx, hiheader2->esize, hiheader2->count,
					hiheader2->idx, hiheader2->esize, hiheader2->count
				);
				printf("    loheader  0x%08X, 0x%08X, 0x%08X    %d, %d, %d\n",
					loheader2->idx, loheader2->esize, loheader2->count,
					loheader2->idx, loheader2->esize, loheader2->count
				);
				break;
			}
			printf("  [%7d] %7d B    ", hiheader2->idx, hiheader2->count*hiheader2->esize);
			ansval_printr(unverdata, hiheader2->esize, hiheader2->count, 0);
		}
		printf("\n");
		CRASH("");
	}
	++popcount;
#endif
}
#endif

//LIFO Bypass Coder
typedef struct _BitPackerLIFO//bwd enc / fwd dec
{
	uint32_t state;
	int enc_nwritten, dec_navailable;//bitcounts, only for tracking renorms
	unsigned char *dstbwdptr;
	const unsigned char *srcfwdptr, *streamend;
} BitPackerLIFO;
AWM_INLINE void bitpacker_enc_init(BitPackerLIFO *ec, const unsigned char *bufstart, unsigned char *bufptr0_OOB)
{
	memset(ec, 0, sizeof(*ec));
	ec->state=1ULL<<16;
	ec->enc_nwritten=17;
	ec->streamend=bufstart;
	ec->dstbwdptr=bufptr0_OOB;
}
AWM_INLINE void bitpacker_dec_init(BitPackerLIFO *ec, const unsigned char *bufptr0_start, const unsigned char *bufend)
{
	memset(ec, 0, sizeof(*ec));
	ec->srcfwdptr=bufptr0_start+4;
	ec->streamend=bufend;
	ec->state=*(const uint32_t*)bufptr0_start;
	ec->dec_navailable=floor_log2(ec->state)+1;
}
AWM_INLINE void bitpacker_enc_flush(BitPackerLIFO *ec)
{
	ec->dstbwdptr-=4;
#ifdef _DEBUG
	if(ec->dstbwdptr<ec->streamend)
		CRASH("IntPacker Encoder OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
	*(uint32_t*)ec->dstbwdptr=ec->state;
}
AWM_INLINE void bitpacker_enc(BitPackerLIFO *ec, int inbits, int sym)
{
	//uint32_t s0=ec->state, s1;//
#ifdef _DEBUG
	if(!inbits)
		CRASH("BitPacker inbits=0");
#endif
	//renorm then push inbits
	ec->enc_nwritten+=inbits;
	if(ec->enc_nwritten>32)//renorm on overflow
	{
		ec->enc_nwritten-=16;
		ec->dstbwdptr-=2;
#ifdef _DEBUG
		if(ec->dstbwdptr<ec->streamend)
			CRASH("IntPacker OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
		*(uint16_t*)ec->dstbwdptr=(uint16_t)ec->state;
		ec->state>>=16;
#ifdef ANS_VAL
		ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
	}
	//s1=ec->state;//
	ec->state=ec->state<<inbits|sym;

	//if(ec->state==0x8124E830)//
	//	printf("");

#ifdef ANS_VAL
	ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
}
AWM_INLINE int bitpacker_dec(BitPackerLIFO *ec, int outbits)
{
	//uint32_t s0=ec->state, s1;//
	int sym=ec->state&((1ULL<<outbits)-1);
	
#ifdef _DEBUG
	if(!outbits)
		CRASH("BitPacker outbits=0");
#endif
	//pop outbits then renorm
#ifdef ANS_VAL
	ansval_check(&ec->state, sizeof(ec->state), 1);
#endif
	ec->dec_navailable-=outbits;
	ec->state>>=outbits;
	//s1=ec->state;//
	if(ec->dec_navailable<=16)
	{
#ifdef ANS_VAL
		ansval_check(&ec->state, sizeof(ec->state), 1);
#endif
		ec->dec_navailable+=16;
#ifdef _DEBUG
		if(ec->srcfwdptr>=ec->streamend)
			CRASH("IntPacker OOB:  srcfwdptr = 0x%016llX >= 0x%016llX", (uint64_t)ec->srcfwdptr, (uint64_t)ec->streamend);
#endif
		ec->state=ec->state<<16|*(const uint16_t*)ec->srcfwdptr;
		ec->srcfwdptr+=2;
	}
	//if(ec->state==0x81242049)//
	//	printf("");
	return sym;
}

typedef struct _rANS_SIMD_SymInfo	//16 bytes/level	4KB/ctx = 1<<12 bytes
{
	uint32_t smax, invf, cdf;
	uint16_t negf, sh;
} rANS_SIMD_SymInfo;
static void enc_hist2stats(int *hist, rANS_SIMD_SymInfo *syminfo, uint32_t *ctxmask, int ctxidx)
{
	int sum=0, count=0, ks, rare;
#ifdef ESTIMATE_SIZE
	int count0, sum0;
#endif
	for(ks=0;ks<256;++ks)
	{
		int freq=hist[ks];
		sum+=freq;
		count+=freq!=0;
	}
	rare=sum<12*256/8;
	ctxmask[ctxidx>>5]|=(uint32_t)rare<<(ctxidx&31);
#ifdef ESTIMATE_SIZE
	count0=count;
	sum0=sum;
#endif
	if(rare)
	{
		for(ks=0;ks<256;++ks)//bypass
			hist[ks]=1;
		sum=256;
		count=256;
	}
	else if(count==1)//disallow degenerate distribution
	{
		for(ks=0;ks<256;++ks)
		{
			int freq=hist[ks];
			if(freq==(1<<PROBBITS))
			{
				--freq;
				if(!ks)
					++hist[ks+1];
				else
					++hist[ks-1];
				break;
			}
		}
		count=2;
	}
	{
		int sum2=0, ks2;
		for(ks=0, ks2=0;ks<256;++ks)//absent symbols get zero freqs
		{
			int freq=hist[ks];
			hist[ks]=(int)(sum2*((1ULL<<PROBBITS)-count)/sum)+ks2;
			ks2+=freq!=0;
			sum2+=freq;
		}
		//for(ks=0;ks<256;++ks)//never allows zero freqs	INEFFICIENT
		//{
		//	int freq=hist[ks];
		//	hist[ks]=(int)(sum2*((1ULL<<PROBBITS)-256)/sum)+ks;
		//	sum2+=freq;
		//}
	}
#ifdef ESTIMATE_SIZE
	{
		double e=sum0;
		if(count==count0)
		{
			double norm=1./0x1000;
			e=0;
			for(ks=0;ks<256;++ks)//estimate
			{
				int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
				if(freq)
				{
					double p=freq*norm;
					e-=p*log(p);
				}
				if(e!=e)
					CRASH("");
			}
			e*=sum/8./M_LN2;//log2(x) = ln(x)/ln2
		}
		if(ctxidx&&!(ctxidx%NCTX))
			printf("\n");
		printf("%c  ctx %3d  %12.2lf / %9d bytes%8.2lf%%  %3d %s",
			ctxidx<3*NCTX?"YUV"[ctxidx/NCTX]:"yuv"[ctxidx-3*NCTX],
			ctxidx, e, sum0, 100.*e/sum0, count0, count==count0?"levels":"bypass"
		);
		if(count==count0&&count<256)
		{
			int fmax=0, ks;
			printf(" %3d", count);
			for(ks=0;ks<256;++ks)
			{
				int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
				if(fmax<freq)
					fmax=freq;
			}
			for(ks=0;ks<256;++ks)
			{
				int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
				if(!(ks&15))
					printf(" ");

			//	int shade=48+freq*(255-48)/fmax;
			//	colorprintf(shade<<16|shade<<8|shade, freq?0x808080:COLORPRINTF_BK_DEFAULT, "%c", "0123456789ABCDEF"[ks&15]);
				//int shade=freq*255/fmax;
				//colorprintf(freq?0xFFFFFF:0x808080, shade<<16|0<<8|shade, "%c", "0123456789ABCDEF"[ks&15]);

				printf("%c", freq?"0123456789ABCDEF"[ks&15]:'-');
			}
		}
		printf("\n");
	}
#endif
	{
		int next=1<<PROBBITS, ks;
		for(ks=255;ks>=0;--ks)
		{
			rANS_SIMD_SymInfo *info=syminfo+ks;
			int curr=hist[ks];
			int freq=next-curr;
			next=curr;
			hist[ks]=freq;
			info->smax=(freq<<(RANS_STATE_BITS-PROBBITS))-1;//rescale freq to match the rANS state, and decrement to use _mm256_cmpgt_epi32 instead of '>='
			info->cdf=curr;
			info->negf=(1<<PROBBITS)-freq;
			//encoding:  state  =  q<<16|(cdf+r)
			//div-free:  state  =  q*M+cdf+state-q*freq  =  state+q*(M-freq)+cdf  =  state+(state*invf>>sh)*(M-freq)+cdf
			//sh = floor_log2(freq)+32
			//invf = ceil(2^sh/freq)		state is 31 bits
			if(freq<2)
			{
				//freq=1
				//ideally  q = state*inv(1)>>sh(1) = state*2^32>>32
				//here  q' = state*(2^32-1)>>32 = floor(state-state/2^32) = state-1  if  1 <= x < 2^32
				//enc  state = (state/1)*M+cdf+state%1  =  state+q*(M-1)+cdf
				//but  q' = state-1
				//so  state = state+(state-1+1)*(M-1)+cdf  =  state+q'*(M-1)+(cdf+M-1)
				info->sh=0;
				info->invf=0xFFFFFFFF;
				info->cdf+=(1<<PROBBITS)-1;
			}
			else
			{
				uint64_t inv;

				info->sh=floor_log2(freq);//eg: x/2 = x*0x80000000>>32>>0
				inv=((0x100000000ULL<<info->sh)+freq-1)/freq;
				info->invf=(uint32_t)inv;
				if(inv>>32)
				{
					--info->sh;
					info->invf=(uint32_t)(inv>>1);
				}
			}
		}
	}
}
static void enc_packhist(BitPackerLIFO *ec, const int *hist, uint32_t *ctxmask, int ctxidx)//histogram must be normalized to PROBBITS, with spike at 128
{
	int sum=0;
	unsigned short CDF[257];
	int ks;

	if(ctxmask[ctxidx>>5]>>(ctxidx&31)&1)
		return;

	for(ks=0;ks<256;++ks)//integrage to zigzag CDF to be packed backwards
	{
		int sym=((ks>>1^-(ks&1))+128)&255;
		int freq=hist[sym];
		CDF[ks]=sum;//separate buffer for faster access in 2nd loop
		sum+=freq;
	}
	CDF[256]=1<<PROBBITS;
	
	{
		int cdfW=CDF[0];
		int CDFlevels=1<<PROBBITS;
		int startsym=0;
		int ks;
		for(ks=1;ks<=256;++ks)//push GR.k
		{
			int next=CDF[ks], freq=next-cdfW;
			int nbypass=floor_log2(CDFlevels);
			if(ks>1)
				nbypass-=7;
			if(nbypass<0)
				nbypass=0;
			CDF[ks]=nbypass<<PROBBITS|freq;
			cdfW=next;
			CDFlevels-=freq;
			startsym=ks;
			if(!CDFlevels)
				break;
		}
		for(ks=startsym;ks>0;--ks)//encode GR
		{
			int freq=CDF[ks], nbypass=freq>>PROBBITS;
			freq&=(1<<PROBBITS)-1;
			{
				int nzeros=freq>>nbypass, bypass=freq&((1<<nbypass)-1);
#ifdef ANS_VAL
				ansval_push(&freq, sizeof(freq), 1);
#endif
				if(nbypass)
					bitpacker_enc(ec, nbypass, bypass);
				bitpacker_enc(ec, 1, 1);
				while(nzeros--)
					bitpacker_enc(ec, 1, 0);
#ifdef ANS_VAL
				//ansval_push(&ks, sizeof(ks), 1);
#endif
			}
		}
	}
}
static void dec_unpackhist(BitPackerLIFO *ec, unsigned *CDF2sym, uint32_t *ctxmask, int ctxidx)
{
	unsigned short hist[257];
	int ks;
	int sum=0;

	if(ctxmask[ctxidx>>5]>>(ctxidx&31)&1)//rare context
	{
		for(ks=0;ks<256;++ks)//bypass
			hist[ks]=(1<<PROBBITS)/256;
	}
	else
	{
		unsigned short CDF[257]={0};
		int CDFlevels=1<<PROBBITS;
		int cdfW=0;
		CDF[0]=0;
		for(ks=0;ks<256;++ks)//decode GR
		{
			int freq=-1;//stop bit doesn't count
			int nbypass=floor_log2(CDFlevels);
			int ks2=ks+1;
			int bit=0;

			if(ks2>1)
				nbypass-=7;
			if(nbypass<0)
				nbypass=0;
#ifdef ANS_VAL
			//ansval_check(&ks2, sizeof(ks2), 1);
#endif
			do
			{
				bit=bitpacker_dec(ec, 1);
				++freq;
			}while(!bit);
			if(nbypass)
				freq=freq<<nbypass|bitpacker_dec(ec, nbypass);
#ifdef ANS_VAL
			ansval_check(&freq, sizeof(freq), 1);
#endif

			cdfW+=freq;
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
		for(ks=0;ks<256;++ks)//undo zigzag
		{
			int sym=((ks>>1^-(ks&1))+128)&255;
			hist[sym]=CDF[ks];
		}
	}
	for(ks=0;ks<256;++ks)//integrate
	{
		int freq=hist[ks];
		hist[ks]=sum;
		sum+=freq;
	}
	hist[256]=1<<PROBBITS;
	for(ks=0;ks<256;++ks)//CDF2sym contains {freq, (state&0xFFF)-cdf, sym}
	{
		int cdf=hist[ks], next=hist[ks+1], freq=next-cdf;
		int val=(freq<<PROBBITS|0)<<8|ks;
		int ks2;
		for(ks2=cdf;ks2<next;++ks2, val+=1<<8)
			CDF2sym[ks2]=val;
	}
}

#ifdef ESTIMATE_SIZE
static int32_t ehist[3][256]={0};
#endif
static int32_t hist[3][NCTX][256]={0};
int s02_codec(const char *command, const char *srcfn, const char *dstfn)
{
	unsigned char *srcbuf=0, *srcptr=0, *srcend=0;
	ptrdiff_t srcsize=0, dstsize=0;
	int fwd=0;
	int32_t iw=0, ih=0;
	ptrdiff_t res=0, usize=0, csize=0;
	unsigned char *image=0;
	unsigned char *dstbuf=0;
	unsigned char *streamptr=0, *streamend=0;
	unsigned char *imptr=0;
	int bestrct=0;
#ifdef LOUD
	double t=time_sec();
#endif
#ifdef PRINTGR
	long long gr_symlen=0, gr_bypsum=0;
#endif
	uint32_t ctxmask[2]={0};
	const int32_t statssize=(int)sizeof(rANS_SIMD_SymInfo[3*NCTX<<8]);
	rANS_SIMD_SymInfo *stats=0;
	uint16_t *ctxbuf=0, *ctxptr=0;
	const int32_t CDF2syms_size=(int)sizeof(int[3*NCTX<<PROBBITS]);
	uint32_t *CDF2syms=0;

	srcbuf=load_file(srcfn, &srcsize);
	if(!srcbuf)
	{
		CRASH("");
		return 1;
	}
	srcptr=srcbuf;
	srcend=srcbuf+srcsize;
	{
		int tag=*(uint16_t*)srcptr;
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('0'|'2'<<8))
		{
			printf("Unsupported file  \"%s\"\n", srcfn);
			CRASH("");
			free(srcbuf);
			return 1;
		}
		srcptr+=2;
	}
	if(fwd)
	{
		if(*srcptr++!='\n')
		{
			printf("Unsupported file\n");
			CRASH("");
			free(srcbuf);
			return 1;
		}
		while(*srcptr=='#')//skip comment(s)
		{
			++srcptr;
			while(*srcptr++!='\n');
		}
		iw=0;
		while((uint32_t)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++-'0';
		while(*srcptr<=' ')++srcptr;//skip space
		ih=0;
		while((uint32_t)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++-'0';
		if(*srcptr++!='\n')
		{
			printf("Unsupported file\n");
			CRASH("");
			free(srcbuf);
			return 1;
		}
		while(*srcptr=='#')//skip comment(s)
		{
			++srcptr;
			while(*srcptr++!='\n');
		}
		if(memcmp(srcptr, "255\n", 4))
		{
			printf("Unsupported file\n");
			CRASH("");
			free(srcbuf);
			return 1;
		}
		srcptr+=4;
		image=srcptr;
		guide_save(image, iw, ih);
	}
	else
	{
		iw=((uint32_t*)srcptr)[0];
		ih=((uint32_t*)srcptr)[1];
		srcptr+=sizeof(uint32_t[2]);
	}
	if(iw<1||ih<1)
	{
		printf("Invalid file\n");
		CRASH("");
		free(srcbuf);
		return 1;
	}
	res=(ptrdiff_t)iw*ih;
	usize=3*res;
	dstbuf=(unsigned char*)malloc(usize);
	if(!dstbuf)
	{
		printf("Alloc error");
		CRASH("");
		free(srcbuf);
		return 1;
	}
	if(fwd)
	{
		int64_t counters[OCH_COUNT]={0}, minerr=0;
		int prev[OCH_COUNT]={0};

		//analysis
		imptr=image;
		while(imptr<srcend)
		{
			int r, g, b, rg, gb, br;

			r=imptr[0];
			g=imptr[1];
			b=imptr[2];
			imptr+=3*3;//stride 3 pixels
			rg=r-g;
			gb=g-b;
			br=b-r;
			counters[0]+=abs(r	-prev[0]);
			counters[1]+=abs(g	-prev[1]);
			counters[2]+=abs(b	-prev[2]);
			counters[3]+=abs(rg	-prev[3]);
			counters[4]+=abs(gb	-prev[4]);
			counters[5]+=abs(br	-prev[5]);
			prev[0]=r;
			prev[1]=g;
			prev[2]=b;
			prev[3]=rg;
			prev[4]=gb;
			prev[5]=br;
		}
		imptr=image;
		{
			int kt;

#ifdef ESTIMATE_SIZE
			for(kt=0;kt<OCH_COUNT;++kt)
				printf("%d %16lld\n", kt, counters[kt]);
			printf("\n");
#endif
			for(kt=0;kt<_countof(rct_indices);++kt)
			{
				const unsigned char *rct=rct_indices[kt];
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
#ifdef ESTIMATE_SIZE
				printf("RCT%02d %16lld%s\n", kt, currerr, kt==bestrct?" <-":"");
#endif
			}
		}
		streamptr=dstbuf+usize;
		streamend=streamptr;

		stats=(rANS_SIMD_SymInfo*)malloc(statssize);
		ctxbuf=(uint16_t*)malloc(usize<<1);
		if(!stats||!ctxbuf)
		{
			CRASH("Alloc error");
			return 1;
		}
		ctxptr=ctxbuf;
	}
	else
	{
		imptr=dstbuf;
		streamptr=srcptr;
		streamend=srcbuf+srcsize;

		bestrct=*streamptr++;

		memcpy(ctxmask, streamptr, sizeof(ctxmask));
		streamptr+=sizeof(ctxmask);
		

		CDF2syms=(uint32_t*)malloc(CDF2syms_size);
		if(!CDF2syms)
		{
			CRASH("Alloc error");
			return 1;
		}
		{
			BitPackerLIFO ec;
			int kc;

			bitpacker_dec_init(&ec, streamptr, streamend);
			for(kc=0;kc<3*NCTX;++kc)
				dec_unpackhist(&ec, CDF2syms+((ptrdiff_t)kc<<PROBBITS), ctxmask, kc);
			streamptr=(unsigned char*)(size_t)ec.srcfwdptr;
		}

		csize=srcsize;
	}
	{
		int
			yidx=rct_indices[bestrct][3+0],
			uidx=rct_indices[bestrct][3+1],
			vidx=rct_indices[bestrct][3+2],
			uhelpidx=rct_indices[bestrct][6+0],
			vhelpidx=rct_indices[bestrct][6+1];
		int32_t ky, kx, idx;
		int32_t psize=0;
		int16_t *pixels=0;
		int32_t paddedwidth=iw+16;
		uint32_t state=0;
#ifdef USE_L1
		int32_t coeffs[3][L1NPREDS+1];
		//{
		//	{
		//		0,
		//		//100000,	//0	N
		//		//100000,	//1	W
		//		// 80000,	//2	3*(N-NN)+NNN
		//		// 80000,	//3	3*(W-WW)+WWW
		//		//150000,	//4	N+W-NW
		//		// 50000,	//5	W+NE-N
		//		// 50000,	//6	N+NE-NNE
		//		// 50000,	//7	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4
		//		//0,		//bias
		//	}
		//};
#ifdef USE_CASCADE
#define ECTXBITS 4
		const int32_t ecoeffsize=(int32_t)sizeof(int32_t[3<<ECTXBITS][L1NPREDS2+1]);
		int32_t *ecoeffs=(int32_t*)malloc(ecoeffsize);
		//int32_t coeffs2[3][L1NPREDS+1]={0};

		if(!ecoeffs)
		{
			printf("Alloc error");
			CRASH("");
		}
		memset(ecoeffs, 0, ecoeffsize);
#endif
		FILLMEM((int32_t*)coeffs, 100000, sizeof(coeffs), sizeof(int));
		//memfill(coeffs[1], coeffs[0], sizeof(int[2][L1NPREDS+1]), sizeof(int[L1NPREDS+1]));
#endif

		if(!fwd)
		{
#ifdef _DEBUG
			if(streamptr>streamend)
				CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=*(uint32_t*)streamptr;
			streamptr+=sizeof(state);
		}
		else
			memset(hist, 0, sizeof(hist));

		psize=(int32_t)sizeof(int16_t[4*3*3])*(iw+16);//4 padded rows * 3 channels * {pixels, residuals, ctx}
		pixels=(int16_t*)malloc(psize);
		if(!pixels)
		{
			printf("Alloc error\n");
			CRASH("");
			free(srcbuf);
			free(dstbuf);
			return 1;
		}
		memset(pixels, 0, psize);
		//FILLMEM((uint16_t*)stats, 0x8000, sizeof(stats), sizeof(int16_t));
		for(ky=0, idx=0;ky<ih;++ky)
		{
			char yuv[4]={0};
			int16_t *rows[]=
			{
				pixels+(paddedwidth*((ky-0)&3)+8)*3*3,
				pixels+(paddedwidth*((ky-1)&3)+8)*3*3,
				pixels+(paddedwidth*((ky-2)&3)+8)*3*3,
				pixels+(paddedwidth*((ky-3)&3)+8)*3*3,
			};
			for(kx=0;kx<iw;++kx, ++idx)
			{
				int offset=0, kc;

				if(fwd)
				{
					yuv[0]=imptr[yidx]-128;
					yuv[1]=imptr[uidx]-128;
					yuv[2]=imptr[vidx]-128;
				}
				for(kc=0;kc<3;++kc)
				{
					int32_t
						NNN	=rows[3][0+0*3*3],
						NNNE	=rows[3][0+1*3*3],
						NNW	=rows[2][0-1*3*3],
						NN	=rows[2][0+0*3*3],
						NNE	=rows[2][0+1*3*3],
						NNEEE	=rows[2][0+3*3*3],
						NWW	=rows[1][0-2*3*3],
						NW	=rows[1][0-1*3*3],
						N	=rows[1][0+0*3*3],
						NE	=rows[1][0+1*3*3],
						NEE	=rows[1][0+2*3*3],
						NEEE	=rows[1][0+3*3*3],
						NEEEE	=rows[1][0+4*3*3],
						NEEEEE	=rows[1][0+5*3*3],
						WWWWW	=rows[0][0-5*3*3],
						WWWW	=rows[0][0-4*3*3],
						WWW	=rows[0][0-3*3*3],
						WW	=rows[0][0-2*3*3],
						W	=rows[0][0-1*3*3],
#ifdef USE_CASCADE
						eNNN	=rows[3][2+0*3*3],
						eNNW	=rows[2][2-1*3*3],
						eNN	=rows[2][2+0*3*3],
						eNNE	=rows[2][2+1*3*3],
						eNNEEE	=rows[2][2+3*3*3],
						eNWW	=rows[1][2-2*3*3],
						eNW	=rows[1][2-1*3*3],
						eN	=rows[1][2+0*3*3],
						eNE	=rows[1][2+1*3*3],
						eNEE	=rows[1][2+2*3*3],
						eNEEE	=rows[1][2+3*3*3],
						eNEEEE	=rows[1][2+4*3*3],
						eNEEEEE	=rows[1][2+5*3*3],
						eWWWWW	=rows[0][2-5*3*3],
						eWWWW	=rows[0][2-4*3*3],
						eWWW	=rows[0][2-3*3*3],
						eWW	=rows[0][2-2*3*3],
						eW	=rows[0][2-1*3*3],
#endif
						cNEEE	=rows[1][1+3*3*3],
						cW	=rows[0][1-1*3*3];
					int32_t predc, vmax, vmin;
					int32_t error;
					int32_t ctx;
#ifdef USE_L1
#ifdef USE_CASCADE
					int32_t ectx=((N+W)/2+128)>>(8-ECTXBITS+(kc!=0))&((1<<ECTXBITS)-1);//try removing (kc!=0)
					int32_t *curr_ecoeffs=ecoeffs+(L1NPREDS2+1)*(kc<<ECTXBITS|ectx);
#endif
					int32_t pred1, pred2;
					int32_t preds[]=
					{
#define PRED(W0, EXPR) EXPR,
						PREDLIST
#undef  PRED
					};
#ifdef USE_CASCADE
					int32_t preds2[]=
					{
#define PRED(EXPR) EXPR,
						PREDLIST2
#undef  PRED
					};
					(void)NEEEE;
					(void)WWWW;
					(void)eNEEEE;
					(void)eWWWW;
#endif
					(void)NNNE;
					(void)NNW;
					(void)NNEEE;
					(void)NWW;
					(void)NEEEEE;
					(void)WWWWW;
					{
						int j;
#ifdef USE_CASCADE
						int32_t tmp;
#endif
	#define L1SH 19
	#define L1SH2 19
//	#define L1SH 15
//	#define L1SH2 15
						pred1=0;
						//pred1=coeffs[kc][L1NPREDS];
						for(j=0;j<L1NPREDS;++j)
						{
							pred1+=coeffs[kc][j]*preds[j];
							//if(abs(coeffs[kc][j]>>4)>0x7FFF)//never hit on LPCB
							//	CRASH("");
							//pred1+=(coeffs[kc][j]>>4)*preds[j];
							//pred1+=(coeffs[kc][j]>>10)*preds[j]>>2;
							//pred1+=(coeffs[kc][j]>>12)*preds[j];//X
							//pred1+=((coeffs[kc][j]>>3)*preds[j]+(1<<15>>1))>>15;//X
							//pred1+=((coeffs[kc][j]>>4)*preds[j]+(1<<15>>1))>>15;//X
						}
#ifdef USE_CASCADE
						tmp=curr_ecoeffs[L1NPREDS2];
						for(j=0;j<L1NPREDS2;++j)
							tmp+=curr_ecoeffs[j]*preds2[j];
#endif
						pred1+=1<<L1SH>>1;
						pred1>>=L1SH;
#ifdef USE_CASCADE
						tmp+=1<<L1SH2>>1;
						tmp>>=L1SH2;
						pred2=pred1+tmp;
#else
						pred2=pred1;
#endif
					}
					predc=pred2;
					vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(predc, vmin, vmax);
#else
					predc=N+W-NW;
					vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					CLAMP2(predc, vmin, vmax);
#endif
					predc+=offset;
					CLAMP2(predc, -128, 127);

					{
						float fval=(float)(cW*cW+1);
						size_t addr=(size_t)&fval;
						int32_t bits=*(int32_t*)addr;
						ctx=(bits>>23)-127;
						if(ctx>NCTX-1)
							ctx=NCTX-1;
#ifdef _DEBUG
						if(ctx<0)
						{
							printf("OOB  ctx %d / %d\n", ctx, NCTX);
							CRASH("");
						}
						if((uint32_t)kc>=3)
						{
							printf("OOB  kc %d / 3\n", kc);
							CRASH("");
						}
#endif
					}

					if(fwd)
					{
						error=(unsigned char)(yuv[kc]-predc+128);
						++hist[kc][ctx][error];
#ifdef _DEBUG
						if(ctxptr>ctxbuf+usize)
						{
							printf("OOB ctxidx %016llX / %016llX"
								, ctxptr-ctxbuf
								, usize
							);
							CRASH("");
						}
#endif
						*ctxptr++=ctx<<8|error;
					}
					else
					{
						//uint32_t s0=state;//
						uint32_t info=CDF2syms[(NCTX*kc+ctx)<<PROBBITS|(state&((1<<PROBBITS)-1))];
						int32_t bias=info>>8&((1<<PROBBITS)-1), freq=info>>(PROBBITS+8);
#ifdef ANS_VAL
						ansval_check(&state, sizeof(state), 1);
#endif
						state=(state>>PROBBITS)*freq+bias;//state = (state>>12)*freq+(rem-cdf)
#ifdef ANS_VAL
						//if(ansvalidx==16597437)//
						//	printf("");

						ansval_check(&state, sizeof(state), 1);
#endif
						if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
						{
							state=state<<RANS_RENORM_BITS|*(uint16_t*)streamptr;
							streamptr+=2;
						}
						yuv[kc]=(char)(info+predc-128);
					}
					error=abs(yuv[kc]-predc);
					rows[0][0]=yuv[kc]-offset;
					rows[0][1]=(2*cW+(error<<GRBITS)+cNEEE)>>2;
#ifdef USE_L1
					rows[0][2]=rows[0][0]-pred1;
					{
						int32_t k, curr=rows[0][0], e;

						e=(curr>pred1)-(curr<pred1);
						//coeffs[kc][L1NPREDS]+=e;
						for(k=0;k<L1NPREDS;++k)
							coeffs[kc][k]+=e*preds[k];
#ifdef USE_CASCADE
						e=(curr>pred2)-(curr<pred2);
						curr_ecoeffs[L1NPREDS2]+=e;
						for(k=0;k<L1NPREDS2;++k)
							curr_ecoeffs[k]+=e*preds2[k];
#endif
					}
#endif
					offset=yuv[kc?vhelpidx:uhelpidx];
					//if(kc)
					//	offset=yuv[vhelpidx];
					//else
					//	offset=yuv[uhelpidx];
					rows[0]+=3;
					rows[1]+=3;
					rows[2]+=3;
					rows[3]+=3;
				}
				if(!fwd)
				{
					imptr[yidx]=yuv[0]+128;
					imptr[uidx]=yuv[1]+128;
					imptr[vidx]=yuv[2]+128;
					guide_check(dstbuf, kx, ky);
				}
				imptr+=3;
			}
		}
#ifdef USE_CASCADE
		free(ecoeffs);
#endif
		free(pixels);
	}
	if(fwd)
	{
		int kc;
		uint32_t state;

		//normalize/integrate hists
		for(kc=0;kc<3*NCTX;++kc)
			enc_hist2stats((int32_t*)hist+(ptrdiff_t)256*kc, stats+(ptrdiff_t)256*kc, ctxmask, kc);
		
		kc=2*NCTX*256;
		state=1<<(RANS_STATE_BITS-RANS_RENORM_BITS);
		while(ctxptr-->ctxbuf)
		{
			//uint32_t s0=state;//
			rANS_SIMD_SymInfo *info=stats+kc+*ctxptr;
#if _DEBUG
			if((size_t)(info-stats)>=(size_t)statssize)
			{
				printf("OOB  info idx %016llX / %016llX\n"
					, info-stats
					, statssize
				);
				CRASH("");
				return 1;
			}
#endif
			if(state>info->smax)
			{
				streamptr-=2;
#ifdef _DEBUG
				if(streamptr<dstbuf)
				{
					printf("OOB  ptr %016llX < start %016llX  inflation %8.4lf%%\n"
						, streamptr
						, dstbuf
						, (double)(dstbuf+usize-streamptr)/(ctxbuf+usize-ctxptr)
					);
					CRASH("");
					return 1;
				}
#endif
				*(uint16_t*)streamptr=(uint16_t)state;
				state>>=RANS_RENORM_BITS;
			}

			//if(state==0x0000055E)//
			//if(ansvalidx==16597436)//
			//	printf("");
#ifdef ANS_VAL
			ansval_push(&state, sizeof(state), 1);
#endif
			state+=(((uint64_t)state*info->invf>>32)>>info->sh)*info->negf+info->cdf;//state = state/freq<<12|(cdf+state%freq)
#ifdef ANS_VAL
			ansval_push(&state, sizeof(state), 1);
#endif
			kc-=NCTX*256;
			if(kc<0)
				kc=2*NCTX*256;
		}
		//flush
		streamptr-=4;
		*(uint32_t*)streamptr=state;

		//pack hists
		{
			BitPackerLIFO ec;

			bitpacker_enc_init(&ec, image, streamptr);
			for(kc=3*NCTX-1;kc>=0;--kc)
				enc_packhist(&ec, (int32_t*)hist+(ptrdiff_t)256*kc, ctxmask, kc);
			bitpacker_enc_flush(&ec);
			streamptr=ec.dstbwdptr;
		}
	}
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			printf("Cannot open \"%s\" for writing\n", dstfn);
			CRASH("");
			free(srcbuf);
			free(dstbuf);
			return 1;
		}
		if(fwd)
		{
			dstsize+=fwrite("02", 1, 2, fdst);
			dstsize+=fwrite(&iw, 1, 4, fdst);
			dstsize+=fwrite(&ih, 1, 4, fdst);
			dstsize+=fwrite(&bestrct, 1, 1, fdst);
			dstsize+=fwrite(ctxmask, 1, sizeof(ctxmask), fdst);
#ifdef _DEBUG
			if(streamptr>streamend)
				CRASH("OOB ptr %016zX > %016zX", streamptr, streamend);
			if(streamptr<dstbuf)
				CRASH("OOB ptr %016zX < %016zX", streamptr, image);
#endif
			dstsize+=fwrite(streamptr, 1, streamend-streamptr, fdst);
		}
		else
		{
			dstsize+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			dstsize+=fwrite(dstbuf, 1, usize, fdst);
		}
		fclose(fdst);
	}
	free(srcbuf);
	free(dstbuf);
	if(fwd)
	{
		free(stats);
		free(ctxbuf);
	}
	else
		free(CDF2syms);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
#ifdef PRINTBITS
		printf("\n");
#endif
#ifdef ESTIMATE_SIZE
		{
			double csizes[3]={0};
			int kc, ks;
			for(kc=0;kc<3;++kc)
			{
				double norm;
				int32_t sum=0;
				for(ks=0;ks<256;++ks)
					sum+=ehist[kc][ks];
				norm=1./sum;
				for(ks=0;ks<256;++ks)
				{
					int32_t freq=ehist[kc][ks];
					if(freq)
						csizes[kc]-=freq*log(freq*norm);
				}
				csizes[kc]*=1./(M_LN2*8);//convert ln->log2
			}
			printf("TYUV %12.2lf %12.2lf %12.2lf %12.2lf\n"
				, csizes[0]+csizes[1]+csizes[2]
				, csizes[0]
				, csizes[1]
				, csizes[2]
			);
		}
#endif
#ifdef PRINTGR
		printf("%12lld bypass %12.6lf bits\n", gr_bypsum, (double)gr_bypsum/usize);
		printf("%12lld symlen %12.6lf bits\n", gr_symlen, (double)gr_symlen/usize);
#endif
		printf("%9d->%9d  %8.4lf%%  %12.6lf\n"
			, srcsize
			, dstsize
			, 100.*dstsize/srcsize
			, (double)srcsize/dstsize
		);
	}
	printf("%c  %12.6lf  %12.6lf MB/s\n"
		, 'D'+fwd
		, t
		, usize/(t*1024*1024)
	);
#endif
	(void)csize;
	(void)time_sec;
	(void)memfill;
	return 0;
}
