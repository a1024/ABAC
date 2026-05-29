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
#if defined _M_X64 || defined __x86_64__
#include<immintrin.h>
#endif
#ifdef PROFILER
#include"util.h"
#endif


#ifdef _MSC_VER
	#define LOUD

	#define RICE_EXPERIMENT
//	#define MEASURE_PSNR

//	#define FIFOVAL
	#define ANS_VAL
#endif

	#define USE_RANS
//	#define NO_RCT
//	#define NEARLOSSLESS
	#define USE_ROWS
//	#define USE_CG


enum
{
	DEPTH=9,
//	RLIMIT=23,
	RLIMIT=12,
#ifdef USE_ROWS
	RSHIFT=3,
#else
	RSHIFT=1,
#endif

	BUFSIZE=512<<10,
	
#ifdef USE_ROWS
	XPAD=8,
	NCH=3,
	NROWS=2,
	NVAL=2,
#endif
#ifdef USE_RANS
	NCTX=8,
	PROBBITS=12,
	RANS_STATE_BITS=31,
	RANS_RENORM_BITS=16,
#endif
};


//runtime
#if 1
#ifndef MEDIAN3V_CLOB
//clobbers A B C
#define MEDIAN3V_CLOB(M, A, B, C)\
	M=A, A=A>B?B:A, B=M>B?M:B,\
	M=B, B=B>C?C:B, C=M>C?M:C,\
	M=A, A=A>B?B:A, M=M>B?M:B
#endif
#ifndef CLAMP2
#define CLAMP2(X, L, H) X=(X)>(L)?X:L, X=(X)<(H)?X:H
#endif
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
#define ROUND32(X) _mm_cvt_ss2si(_mm_set_ss(X))
#define ROUND64(X) _mm_cvtsd_si64(_mm_set_sd(X))
#define TRUNC32(X) _mm_cvtt_ss2si(_mm_set_ss(X))
#define TRUNC64(X) _mm_cvttsd_si64(_mm_set_sd(X))
#if defined _M_X64 || defined __x86_64__
#define LZCNT32 _lzcnt_u32
#define LZCNT64 _lzcnt_u64
#else
#define LZCNT32 __builtin_clz
#define LZCNT64 __builtin_clzll
#endif
#ifndef FLOOR_LOG2
#if defined _M_X64 || defined __x86_64__
#define FLOOR_LOG2(X) (sizeof(X)==8?63-(int32_t)_lzcnt_u64(X):31-_lzcnt_u32((uint32_t)(X)))
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
static int fifoval_check(uint32_t val)
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
		return 1;
		//CRASH("");
	}
	return 0;
}
#endif
#endif


static uint8_t rdbuf[BUFSIZE+sizeof(uint64_t[4])]={0}, wtbuf[BUFSIZE+sizeof(uint64_t[4])]={0};
#if 0
/*
overflow:
|                    ______left______   ______right_____
|                   /                \ /                \
|buf1start ... ... [datastart  buf1end|buf2start  dataend] ...
|                   \________________    _______________/
|                                    size
*/
INLINE uint64_t acme_read(uint8_t **pptr, const ptrdiff_t size, FILE *f)
{
	uint8_t *ptr=*pptr;
	uint64_t data;
	
	memcpy(&data, ptr, sizeof(uint64_t));
	if(ptr>=rdbuf+sizeof(uint64_t)+BUFSIZE-size)
	{
		uint64_t d2;

		fread(rdbuf+sizeof(uint64_t), 1, BUFSIZE, f);
		ptr-=BUFSIZE;
		memcpy(&d2, ptr, sizeof(uint64_t));
		data|=d2;
	}
	*pptr=ptr+size;
	return data;
}
INLINE void acme_write(uint8_t **pptr, const ptrdiff_t size, FILE *f, uint64_t data)
{
	uint8_t *ptr=*pptr;
	
	memcpy(ptr, &data, sizeof(uint64_t));
	if(ptr>=wtbuf+sizeof(uint64_t)+BUFSIZE-size)
	{
		fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, f);
		ptr-=BUFSIZE;
		memcpy(ptr, &data, sizeof(uint64_t));
	}
	*pptr=ptr+size;
}
#endif

#ifdef RICE_EXPERIMENT
static void rice_experiment(void)
{
	enum
	{
		USIZE=48<<20,
	};
	int64_t csize=0;
	int run=0, estim=0;
	int hist[256]={0};
	double esize=0;
	int runestim=0;
	//int runk=(int)round(log2(-1/log2(1-p1)));
	int csize2=0, estim2=0;
	for(int it=0;it<USIZE;++it)
	{
		int x=0, bit, nbypass, limit;
		//double p1=cos(0.001*it);
		//p1*=p1;
		//p1*=p1;
		//p1*=0.25;
		//p1+=0.001;
		double p1=0.0001;
		limit=(int)ROUND64(RAND_MAX*p1);
		for(int it2=0;it2<256;++it2)
		{
			bit=rand()<limit;
			x+=bit;
			if(!bit)
				break;
		}
		esize-=log2((double)(hist[x]+1)/(it+256));
		++hist[x];

		nbypass=31-_lzcnt_u32((estim2>>1)+1);
		csize2+=(x>>nbypass)+1+nbypass;
		estim2+=(2*x-estim2)>>2;

		if(!x)
			++run;
		else
		{
			++csize;//run flag
			if(run)
			{
				nbypass=31-_lzcnt_u32((runestim>>1)+1);
				csize+=(int64_t)(run>>nbypass)+1+nbypass;
				runestim+=(2*run-runestim)>>2;
				run=0;
			}
			nbypass=31-_lzcnt_u32((estim>>1)+1);
			csize+=(int64_t)(x>>nbypass)+1+nbypass;
			estim+=(2*x-estim)>>2;
		}
	}
	csize+=64;
	printf("%9d->%12.2lf  %12.2lf  %12.2lf\n", USIZE, csize/8., csize2/8., esize/8.);
#if 0
	enum
	{
		USIZE=48<<20,
	};
	static const double p1a[]={0.001, 0.002, 0.005, 0.01, 0.015, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5};
	for(int it0=0;it0<_countof(p1a);++it0)
	{
		double p1=p1a[it0];
		int64_t csize=0;
		int run=0, estim=0;
		int hist[256]={0};
		double esize=0;
		int runk=(int)round(log2(-1/log2(1-p1)));
		int csize2=0, estim2=0;
		for(int it=0;it<USIZE;++it)
		{
			int x=0, bit, nbypass;
			for(int it2=0;it2<256;++it2)
			{
				bit=rand()<RAND_MAX*p1;
				x+=bit;
				if(!bit)
					break;
			}
			esize-=log2((double)(hist[x]+1)/(it+256));
			++hist[x];

			nbypass=31-_lzcnt_u32((estim2>>1)+1);
			csize2+=(x>>nbypass)+1+nbypass;
			estim2+=(2*x-estim2)>>2;

			if(!x)
			{
				++run;
				continue;
			}
			++csize;
			if(run)
			{
				csize+=(int64_t)(run>>runk)+1+runk;
				run=0;
			}
			nbypass=31-_lzcnt_u32((estim>>1)+1);
			csize+=(int64_t)(x>>nbypass)+1+nbypass;
			estim+=(2*x-estim)>>2;
		}
		csize+=64;
		printf("%9d->%12.2lf  %12.2lf  %12.2lf  [%8.4lf%%]  %12.3lf\n", USIZE, csize/8., csize2/8., esize/8., 100.*esize/csize, p1);
	}
#endif
	exit(0);
}
#endif
#ifdef USE_RANS
//ANS validation
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
		printf("  inc     %+16lld bytes\n", (uint64_t)nbytes);
		printf("  start   %016zd  %16d\n", istart, 0);
		if(nbytes<0)
		{
			printf("  after   %016lld  %16lld%s\n", (uint64_t)ip2, (uint64_t)(ip2-istart), problems[2]?"  <-":"");
			printf("  before  %016lld  %16lld%s\n", (uint64_t)ip1, (uint64_t)(ip1-istart), problems[1]?"  <-":"");
		}
		else
		{
			printf("  before  %016lld  %16lld%s\n", (uint64_t)ip1, (uint64_t)(ip1-istart), problems[1]?"  <-":"");
			printf("  after   %016lld  %16lld%s\n", (uint64_t)ip2, (uint64_t)(ip2-istart), problems[2]?"  <-":"");
		}
		printf("  end     %016lld  %16lld%s\n", (uint64_t)iend, (uint64_t)size, problems[0]?"  <-":"");
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
}
#else
#define ansval_push(...)
#define ansval_check(...)
#endif

//LIFO Bypass Coder
#if 1
#define BITPACKERMAX 32
typedef struct _BitPackerLIFO//bwd enc / fwd dec
{
	uint64_t state;
	int32_t enc_nwritten, dec_navailable;//bitcounts, only for tracking renorms
	uint8_t *dstbwdptr;
	const uint8_t *srcfwdptr, *streamend;
} BitPackerLIFO;
INLINE void bitpacker_enc_init(BitPackerLIFO *ec, const uint8_t *bufstart, uint8_t *bufptr0_OOB)
{
	memset(ec, 0, sizeof(*ec));
	ec->state=1ULL<<32;
	ec->enc_nwritten=33;
	ec->streamend=bufstart;
	ec->dstbwdptr=bufptr0_OOB;
}
INLINE void bitpacker_dec_init(BitPackerLIFO *ec, const uint8_t *bufptr0_start, const uint8_t *bufend)
{
	memset(ec, 0, sizeof(*ec));
	ec->srcfwdptr=bufptr0_start+8;
	ec->streamend=bufend;
	ec->state=*(const uint64_t*)bufptr0_start;
	ec->dec_navailable=FLOOR_LOG2(ec->state)+1;
}
INLINE void bitpacker_enc_flush(BitPackerLIFO *ec)
{
	ec->dstbwdptr-=8;
#ifdef _DEBUG
	if(ec->dstbwdptr<ec->streamend)
		CRASH("IntPacker Encoder OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
	*(uint64_t*)ec->dstbwdptr=ec->state;
}
INLINE void bitpacker_enc(BitPackerLIFO *ec, int inbits, int sym)
{
#ifdef _DEBUG
	if(inbits>BITPACKERMAX)
		CRASH("BitPacker inbits %d", inbits);
#endif
	//renorm then push inbits
	ec->enc_nwritten+=inbits;
	if(ec->enc_nwritten>64)//renorm on overflow
	{
		ec->enc_nwritten-=32;
		ec->dstbwdptr-=4;
#ifdef _DEBUG
		if(ec->dstbwdptr<ec->streamend)
			CRASH("IntPacker OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
		*(unsigned*)ec->dstbwdptr=(unsigned)ec->state;
		ec->state>>=32;
#ifdef ANS_VAL
		ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
	}
	ec->state=ec->state<<inbits|(uint64_t)sym;
#ifdef ANS_VAL
	ansval_push(&inbits, sizeof(inbits), 1);
	ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
}
INLINE int bitpacker_dec(BitPackerLIFO *ec, int outbits)
{
	int sym;

#ifdef _DEBUG
	if(outbits>BITPACKERMAX)
		CRASH("BitPacker outbits %d", outbits);
#endif
	sym=(int)(ec->state&((1ULL<<outbits)-1));

	//pop outbits then renorm
#ifdef ANS_VAL
	ansval_check(&ec->state, sizeof(ec->state), 1);
	ansval_check(&outbits, sizeof(outbits), 1);
#endif
	ec->dec_navailable-=outbits;
	ec->state>>=outbits;
	if(ec->dec_navailable<=32)
	{
#ifdef ANS_VAL
		ansval_check(&ec->state, sizeof(ec->state), 1);
#endif
		ec->dec_navailable+=32;
#ifdef _DEBUG
		if(ec->srcfwdptr>=ec->streamend)
			CRASH("IntPacker OOB:  srcfwdptr = 0x%016zX >= 0x%016zX", ec->srcfwdptr, ec->streamend);
#endif
		ec->state=ec->state<<32|*(const unsigned*)ec->srcfwdptr;
		ec->srcfwdptr+=4;
	}
	return sym;
}
#endif

//SIMD static-o1 rANS	https://github.com/rygorous/ryg_rans	https://github.com/samtools/htscodecs
typedef struct _rANS_SIMD_SymInfo	//16 bytes/level	4KB/ctx = 1<<12 bytes
{
	uint32_t smax, invf, cdf;
	uint16_t negf, sh;
} rANS_SIMD_SymInfo;
static void enc_hist2stats(int *hist, rANS_SIMD_SymInfo *syminfo, uint64_t *bypassmask, int ctxidx, int sse41, int oneshift)
{
#ifdef ESTIMATE_SIZE
	int count0=0, sum0=0;
#endif
	int sum=0, count=0, ks, rare;
	for(ks=0;ks<256;++ks)
	{
		int freq=hist[ks];
		sum+=freq;
		count+=freq!=0;
	}
	rare=sum<12*256/8;
	*bypassmask|=(uint64_t)rare<<ctxidx;
#ifdef ESTIMATE_SIZE
	count0=count; sum0=sum;
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
		int sum2, ks2;

		for(ks=0, ks2=0, sum2=0;ks<256;++ks)//absent symbols get zero freqs
		{
			int freq=hist[ks];
			hist[ks]=(int)((uint64_t)sum2*((1ULL<<PROBBITS)-(uint64_t)count)/(uint64_t)sum)+ks2;
			ks2+=freq!=0;
			sum2+=freq;
		}
		//for(ks=0, sum2=0;ks<256;++ks)//never allows zero freqs	INEFFICIENT
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
			int ks;
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
			e*=sum/(8.*M_LN2);
		}
		if(ctxidx&&!(ctxidx%NCTX))
			printf("\n");
		printf("%c  ctx %3d  %12.2lf / %9d bytes%10.2lf%%  %3d %s",
			ctxidx<3*NCTX?"YUV"[ctxidx/NCTX]:"yuv"[ctxidx-3*NCTX],
			ctxidx, e, sum0, 100.*e/sum0, count0, count==count0?"levels":"bypass"
		);
		if(count==count0&&count<256)
		{
			int fmax, ks;

			printf(" %3d", count);
			fmax=0;
			for(ks=0;ks<256;++ks)
			{
				int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
				if(fmax<freq)
					fmax=freq;
			}
			for(ks=0;ks<256;++ks)
			{
				int freq, shade;

				freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
				if(!(ks&15))
					printf(" ");

				shade=48+freq*(255-48)/fmax;
				colorprintf(shade<<16|shade<<8|shade, freq?0x808080:COLORPRINTF_BK_DEFAULT, "%c", "0123456789ABCDEF"[ks&15]);
				//int shade=freq*255/fmax;
				//colorprintf(freq?0xFFFFFF:0x808080, shade<<16|0<<8|shade, "%c", "0123456789ABCDEF"[ks&15]);

				//printf("%c", freq?"0123456789ABCDEF"[ks&15]:'-');
			}
		}
		printf("\n");
	}
#ifdef DEBUG_HIST
	if(ctxidx==0)
	{
		const int amplitude=512;
		int ks;

		printf("Context %d: (1 star = %d steps)\n", ctxidx, (1<<PROBBITS)/amplitude);
		for(ks=0;ks<256;++ks)
		{
			int freq, nstars, k;

			freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
			nstars=freq*amplitude>>PROBBITS;
			printf("%3d %4d ", ks, freq);
			for(k=0;k<nstars;++k)
				printf("*");
			printf("\n");
		}
	}
#endif
#endif
	{
		int next=1<<PROBBITS;
		for(ks=255;ks>=0;--ks)
		{
			rANS_SIMD_SymInfo *info=syminfo+ks;
			int curr=hist[ks];
			int freq=next-curr;
			next=curr;
			hist[ks]=freq;
			info->smax=(uint32_t)((freq<<(RANS_STATE_BITS-PROBBITS))-1);//rescale freq to match the rANS state, and decrement to use _mm_cmpgt_epi32 instead of '>='
			info->cdf=(uint32_t)curr;
			info->negf=(uint16_t)((1<<PROBBITS)-freq);
			//encoding:  state  =  q<<16|(cdf+r)
			//div-free:  state  =  q*M+cdf+state-q*freq  =  state+q*(M-freq)+cdf  =  state+(state*invf>>sh)*(M-freq)+cdf
			//sh = FLOOR_LOG2(freq)+32
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

				info->sh=(uint16_t)FLOOR_LOG2(freq);//eg: x/2 = x*0x80000000>>32>>0
				inv=((0x100000000ULL<<info->sh)+(uint64_t)freq-1)/(uint64_t)freq;
				info->invf=(uint32_t)inv;
				if(inv>0xFFFFFFFF)
				{
					--info->sh;
					info->invf=(uint32_t)(inv>>1);
				}
			}
#ifdef PRINT_SHIFTBOUNDS
			if(minsh>info->sh)
				minsh=info->sh;
			if(maxsh<info->sh)
				maxsh=info->sh;
#endif
			if(sse41)
				info->sh=(uint16_t)(1<<(PROBBITS-1-info->sh));
			else if(oneshift)
				info->sh+=32;
		}
	}
}
static void enc_packhist(BitPackerLIFO *ec, const int *hist, uint64_t bypassmask, int ctxidx)//histogram must be normalized to PROBBITS, with spike at 128
{
	uint16_t CDF[257];
	int ks;

	if(bypassmask>>ctxidx&1)
		return;
	{
		int sum=0;
		for(ks=0;ks<256;++ks)//integrage to zigzag CDF to be packed backwards
		{
			int sym=((ks>>1^-(ks&1))+128)&255;
			int freq=hist[sym];
			CDF[ks]=(uint16_t)sum;//separate buffer for faster access in 2nd loop
			sum+=freq;
		}
		CDF[256]=1<<PROBBITS;
	}
	{
		int cdfW=CDF[0];
		int CDFlevels=1<<PROBBITS;
		int startsym=0;
		for(ks=1;ks<=256;++ks)//push GR.k
		{
			int next=CDF[ks], freq=next-cdfW;
			int nbypass=FLOOR_LOG2(CDFlevels);
			if(ks>1)
				nbypass-=7;
			if(nbypass<0)
				nbypass=0;
			CDF[ks]=(uint16_t)(nbypass<<PROBBITS|freq);
			cdfW=next;
			CDFlevels-=freq;
			startsym=ks;
			if(!CDFlevels)
				break;
		}
		for(ks=startsym;ks>0;--ks)//encode GR
		{
			int freq, nbypass, nzeros, bypass;

			freq=CDF[ks];
			nbypass=freq>>PROBBITS;
			freq&=(1<<PROBBITS)-1;
			nzeros=freq>>nbypass, bypass=freq&((1<<nbypass)-1);
#ifdef ANS_VAL
			ansval_push(&freq, sizeof(freq), 1);
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
			ansval_push(&ks, sizeof(ks), 1);
#endif
		}
	}
}
static void dec_unpackhist(BitPackerLIFO *ec, uint32_t *CDF2sym, uint64_t bypassmask, int ctxidx)
{
	uint16_t hist[257];
	int ks;

	if(bypassmask>>ctxidx&1)//rare context
	{
		for(ks=0;ks<256;++ks)//bypass
			hist[ks]=(1<<PROBBITS)/256;
	}
	else
	{
		uint16_t CDF[257]={0};
		int CDFlevels=1<<PROBBITS;
		CDF[0]=0;
		for(ks=0;ks<256;++ks)//decode GR
		{
			int freq, nbypass, ks2, bit;

			freq=-1;//stop bit doesn't count
			nbypass=FLOOR_LOG2(CDFlevels);
			ks2=ks+1;
			if(ks2>1)
				nbypass-=7;
			if(nbypass<0)
				nbypass=0;
#ifdef ANS_VAL
			ansval_check(&ks2, sizeof(ks2), 1);
#endif
			bit=0;
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

			CDF[ks]=(uint16_t)freq;
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
	{
		int sum=0;
		for(ks=0;ks<256;++ks)//integrate
		{
			int freq=hist[ks];
			hist[ks]=(uint16_t)sum;
			sum+=freq;
		}
	}
	hist[256]=1<<PROBBITS;
	for(ks=0;ks<256;++ks)//CDF2sym contains {freq, (state&0xFFF)-cdf, sym}
	{
		int cdf, next, freq, val, ks2;

		cdf=hist[ks];
		next=hist[ks+1];
		freq=next-cdf;
		val=(freq<<PROBBITS|0)<<8|ks;
		for(ks2=cdf;ks2<next;++ks2, val+=1<<8)
			CDF2sym[ks2]=(uint16_t)val;
	}
#ifdef DEBUG_HIST
	if(ctxidx==DEBUG_HIST)
	{
		const int amplitude=512;
		int ks;

		printf("Context %d: (1 star = %d steps)\n", ctxidx, (1<<PROBBITS)/amplitude);
		for(ks=0;ks<256;++ks)
		{
			int freq, nstars, k;

			freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks], nstars=freq*amplitude>>PROBBITS;
			printf("%3d %4d ", ks, freq);
			for(k=0;k<nstars;++k)
				printf("*");
			printf("\n");
		}
	}
#endif
}

static int32_t hists[3][NCTX][256];
#endif
#ifdef MEASURE_PSNR
static void print_psnr(const int64_t *esum, int64_t res, const double *mixweights)
{
	double rmse[4], psnr[4], invres=1./res;

	rmse[0]=(double)esum[0]*invres;
	rmse[1]=(double)esum[1]*invres;
	rmse[2]=(double)esum[2]*invres;
	rmse[3]=(mixweights[0]*rmse[0]+mixweights[1]*rmse[1]+mixweights[2]*rmse[2])/(mixweights[0]+mixweights[1]+mixweights[2]);
	rmse[0]=sqrt(rmse[0]);
	rmse[1]=sqrt(rmse[1]);
	rmse[2]=sqrt(rmse[2]);
	rmse[3]=sqrt(rmse[3]);
	psnr[0]=-20*log10(rmse[0]*(1./255));
	psnr[1]=-20*log10(rmse[1]*(1./255));
	psnr[2]=-20*log10(rmse[2]*(1./255));
	psnr[3]=-20*log10(rmse[3]*(1./255));
	printf(
		"PSNR  %12.6lf %12.6lf %12.6lf %12.6lf\n"
		"RMSE  %12.6lf %12.6lf %12.6lf %12.6lf\n"
		, psnr[3], psnr[0], psnr[1], psnr[2]
		, rmse[3], rmse[0], rmse[1], rmse[2]
	);
}
#endif
static uint8_t log2table[1<<(DEPTH+RSHIFT)];
#ifndef USE_RANS
static uint32_t enctable[1<<(DEPTH+3)];
#endif
static int16_t packsign[1024], *const packsignptr=packsign+512;
int c59_codec(int argc, char **argv)
{
	const uint16_t tag='5'|'9'<<8;
	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	int64_t c=0;
	int fwd=0, iw=0, ih=0;
	uint8_t *rdptr=0, *wtptr=0;
	int pred[3]={0}, yuv[3]={0};
	int sym[3]={0};
	int estim[3]={0};
	int64_t usize=0;
#ifdef USE_ROWS
	int psize=0;
	int16_t *pixels=0;
#endif
#ifdef LOUD
	double t=time_sec2();
#endif
#ifdef PROFILER
	void *prof_ctx=prof_start();
#endif
#ifdef MEASURE_PSNR
	int64_t esum_yuv[3]={0}, esum_rgb[3]={0};
#endif

#ifdef RICE_EXPERIMENT
	rice_experiment();
#endif
	if(argc!=3)
	{
		printf("Usage:  \"%s\"  src  dst\n", argv[0]);
		return 1;
	}
	srcfn=argv[1];
	dstfn=argv[2];
	fsrc=fopen(srcfn, "rb");
	if(!fsrc)
	{
		CRASH("Cannot open \"%s\"", srcfn);
		return 1;
	}
	c=fgetc(fsrc);
	c|=(int64_t)fgetc(fsrc)<<8;
	if(c==('P'|'6'<<8))
		fwd=1;
	else if(c!=tag)
	{
		CRASH("Unsupported file \"%s\"", srcfn);
		return 1;
	}
	if(fwd)
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
			iw=10*iw+(int)c-'0';
			c=fgetc(fsrc);
		}
		while(c<=' ')
			c=fgetc(fsrc);
		ih=0;
		while((uint32_t)(c-'0')<10)
		{
			ih=10*ih+(int)c-'0';
			c=fgetc(fsrc);
		}
		if(c!='\n')
		{
			CRASH("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		while(c<=' ')
			c=fgetc(fsrc);
		c|=(int64_t)fgetc(fsrc)<< 8;
		c|=(int64_t)fgetc(fsrc)<<16;
		c|=(int64_t)fgetc(fsrc)<<24;
		if(c!=('2'|'5'<<8|'5'<<16|'\n'<<24))
		{
			CRASH("Unsupported file \"%s\"", srcfn);
			return 1;
		}
	}
	else
	{
		iw=0;
		ih=0;
		fread(&iw, 1, 3, fsrc);
		fread(&ih, 1, 3, fsrc);
	}
	if(iw<1||ih<1)
	{
		CRASH("Unsupported file \"%s\"", srcfn);
		return 1;
	}
	usize=(int64_t)3*iw*ih;
	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
#ifdef USE_ROWS
	psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error\n");
		fclose(fsrc);
		fclose(fdst);
		return 1;
	}
	memset(pixels, 0, psize);
#endif
	memset(rdbuf, 0, sizeof(rdbuf));
	memset(wtbuf, 0, sizeof(wtbuf));
	rdptr=rdbuf+sizeof(uint64_t)+BUFSIZE;
	wtptr=wtbuf+sizeof(uint64_t);
	for(int ks=0;ks<1<<(DEPTH+RSHIFT);++ks)
	{
		int val=(ks>>RSHIFT)+1;
		val=LZCNT32(val);
		val^=31;
		if(val>7)
			val=7;
		log2table[ks]=val;
	}
	if(fwd)
	{
		static uint8_t *const rdend=rdbuf+sizeof(uint64_t)+BUFSIZE-3;
		static uint8_t *const wtend=wtbuf+sizeof(uint64_t)+BUFSIZE-sizeof(uint64_t);
		uint64_t cache=0;
		int nbits=0;

		for(int k=0;k<1024;++k)
		{
			int val=k-512;
			val=val<<(32-DEPTH)>>(32-DEPTH-1);
			val=val^val>>31;
			packsign[k]=val;
		}
#ifndef USE_RANS
		for(int ks=0;ks<1<<DEPTH;++ks)
		{
			for(int kb=0;kb<8;++kb)
			{
				uint32_t code;
				int nbypass, codelen, nzeros;
		
				nzeros=ks>>kb;
				codelen=nzeros+1+kb;
				code=ks;
				nbypass=kb^31;
				code<<=nbypass;
				code|=1<<31;
				code>>=nbypass;
				if(nzeros>RLIMIT-1)
					code=ks, codelen=RLIMIT+DEPTH;
				enctable[ks<<3|kb]=code<<8|codelen;
			}
		}
#endif
		fwrite(&tag, 1, 2, fdst);
		fwrite(&iw, 1, 3, fdst);
		fwrite(&ih, 1, 3, fdst);
		for(int ky=0;ky<ih;++ky)
		{
#ifdef USE_ROWS
			int16_t *rows[]=
			{
				pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			};
#endif
			for(int kx=0;kx<iw;++kx)
			{
				uint64_t reg;
#ifndef USE_RANS
				uint64_t code[3];
				int codelen[3];
#endif
				int nbypass[3];
				memcpy(&reg, rdptr, sizeof(uint64_t));
				if(rdptr>=rdend)
				{
					uint64_t d2;

					fread(rdbuf+sizeof(uint64_t), 1, BUFSIZE, fsrc);
					rdptr-=BUFSIZE;
					memcpy(&d2, rdptr, sizeof(uint64_t));
					reg|=d2;
				}
				rdptr+=3;
				//reg=acme_read(&rdptr, 3, fsrc);
#ifdef USE_ROWS
#ifdef USE_CG
				{
					int N[3], W[3], NW[3];
					
					NW[0]	=rows[1][0+(0-1*NCH)*NROWS*NVAL];
					NW[1]	=rows[1][0+(1-1*NCH)*NROWS*NVAL];
					NW[2]	=rows[1][0+(2-1*NCH)*NROWS*NVAL];
					N[0]	=rows[1][0+(0+0*NCH)*NROWS*NVAL];
					N[1]	=rows[1][0+(1+0*NCH)*NROWS*NVAL];
					N[2]	=rows[1][0+(2+0*NCH)*NROWS*NVAL];
					W[0]	=rows[0][0+(0-1*NCH)*NROWS*NVAL];
					W[1]	=rows[0][0+(1-1*NCH)*NROWS*NVAL];
					W[2]	=rows[0][0+(2-1*NCH)*NROWS*NVAL];
					NW[0]=N[0]+W[0]-NW[0];
					NW[1]=N[1]+W[1]-NW[1];
					NW[2]=N[2]+W[2]-NW[2];
					MEDIAN3V_CLOB(pred[0], N[0], W[0], NW[0]);
					MEDIAN3V_CLOB(pred[1], N[1], W[1], NW[1]);
					MEDIAN3V_CLOB(pred[2], N[2], W[2], NW[2]);
				}
#else
				pred[0]=(rows[1][0+(0+0*NCH)*NROWS*NVAL]+rows[0][0+(0-1*NCH)*NROWS*NVAL])>>1;
				pred[1]=(rows[1][0+(1+0*NCH)*NROWS*NVAL]+rows[0][0+(1-1*NCH)*NROWS*NVAL])>>1;
				pred[2]=(rows[1][0+(2+0*NCH)*NROWS*NVAL]+rows[0][0+(2-1*NCH)*NROWS*NVAL])>>1;
#endif
				estim[0]=rows[0][1+(0-1*NCH)*NROWS*NVAL];
				estim[1]=rows[0][1+(1-1*NCH)*NROWS*NVAL];
				estim[2]=rows[0][1+(2-1*NCH)*NROWS*NVAL];
#endif
				nbypass[0]=log2table[estim[0]];
				nbypass[1]=log2table[estim[1]];
				nbypass[2]=log2table[estim[2]];
#ifdef USE_RANS
				yuv[0]=(uint8_t)(reg>> 8);
				yuv[1]=(uint8_t)(reg>>16);
				yuv[2]=(uint8_t)(reg>> 0);
#else
				yuv[0]=(uint8_t)(reg>> 0);
				yuv[1]=(uint8_t)(reg>> 8);
				yuv[2]=(uint8_t)(reg>>16);
#ifdef MEASURE_PSNR
				int rgb0[3];
				rgb0[0]=yuv[0];
				rgb0[1]=yuv[1];
				rgb0[2]=yuv[2];
#endif
#ifndef NO_RCT
				yuv[0]-=yuv[1];
				yuv[2]-=yuv[1];
				yuv[1]+=(yuv[0]+yuv[2])>>2;
				yuv[2]-=yuv[0]>>2;
#endif
#ifdef NEARLOSSLESS
				sym[0]=yuv[0]-pred[0];
				sym[1]=yuv[1]-pred[1];
				sym[2]=yuv[2]-pred[2];
				sym[0]>>=2;
				sym[1]>>=2;
				sym[2]>>=2;
#ifdef MEASURE_PSNR
				int yuv0[3];
				yuv0[0]=yuv[0];
				yuv0[1]=yuv[1];
				yuv0[2]=yuv[2];
#endif
				yuv[0]=sym[0]<<2;
				yuv[1]=sym[1]<<2;
				yuv[2]=sym[2]<<2;
				yuv[0]+=pred[0];
				yuv[1]+=pred[1];
				yuv[2]+=pred[2];
				//CLAMP2(yuv[0], -255, 255);
				//CLAMP2(yuv[1],    0, 255);
				//CLAMP2(yuv[2], -255, 255);
#ifdef MEASURE_PSNR
				yuv0[0]-=yuv[0];
				yuv0[1]-=yuv[1];
				yuv0[2]-=yuv[2];
				esum_yuv[0]+=(int64_t)yuv0[0]*yuv0[0];
				esum_yuv[1]+=(int64_t)yuv0[1]*yuv0[1];
				esum_yuv[2]+=(int64_t)yuv0[2]*yuv0[2];
				yuv0[0]=yuv[0];
				yuv0[1]=yuv[1];
				yuv0[2]=yuv[2];
#ifndef NO_RCT
				yuv0[2]+=yuv0[0]>>2;
				yuv0[1]-=(yuv0[0]+yuv0[2])>>2;
				yuv0[2]+=yuv0[1];
				yuv0[0]+=yuv0[1];
#endif
				rgb0[0]-=yuv0[0];
				rgb0[1]-=yuv0[1];
				rgb0[2]-=yuv0[2];
				esum_rgb[0]+=(int64_t)rgb0[0]*rgb0[0];
				esum_rgb[1]+=(int64_t)rgb0[1]*rgb0[1];
				esum_rgb[2]+=(int64_t)rgb0[2]*rgb0[2];
#endif
				sym[0]=packsignptr[sym[0]];
				sym[1]=packsignptr[sym[1]];
				sym[2]=packsignptr[sym[2]];
#else
				sym[0]=packsignptr[yuv[0]-pred[0]];
				sym[1]=packsignptr[yuv[1]-pred[1]];
				sym[2]=packsignptr[yuv[2]-pred[2]];
#endif
				//sym[0]=sym[0]<<23>>23;
				//sym[1]=(int8_t)sym[1];
				//sym[2]=sym[2]<<23>>23;
				//sym[0]=sym[0]<<1^sym[0]>>31;
				//sym[1]=sym[1]<<1^sym[1]>>31;
				//sym[2]=sym[2]<<1^sym[2]>>31;
#endif
#ifdef USE_ROWS
				rows[0][0+(0+0*NCH)*NROWS*NVAL]=yuv[0];
				rows[0][0+(1+0*NCH)*NROWS*NVAL]=yuv[1];
				rows[0][0+(2+0*NCH)*NROWS*NVAL]=yuv[2];
				rows[0][1+(0+0*NCH)*NROWS*NVAL]=(2*rows[0][1+(0-1*NCH)*NROWS*NVAL]+(sym[0]<<RSHIFT)+rows[1][1+(0+3*NCH)*NROWS*NVAL])>>2;
				rows[0][1+(1+0*NCH)*NROWS*NVAL]=(2*rows[0][1+(1-1*NCH)*NROWS*NVAL]+(sym[1]<<RSHIFT)+rows[1][1+(1+3*NCH)*NROWS*NVAL])>>2;
				rows[0][1+(2+0*NCH)*NROWS*NVAL]=(2*rows[0][1+(2-1*NCH)*NROWS*NVAL]+(sym[2]<<RSHIFT)+rows[1][1+(2+3*NCH)*NROWS*NVAL])>>2;
				rows[0]+=NCH*NROWS*NVAL;
				rows[1]+=NCH*NROWS*NVAL;
#else
				pred[0]=yuv[0];
				pred[1]=yuv[1];
				pred[2]=yuv[2];
				estim[0]+=(int)((sym[0]<<RSHIFT)-estim[0])>>(RSHIFT+1);
				estim[1]+=(int)((sym[1]<<RSHIFT)-estim[1])>>(RSHIFT+1);
				estim[2]+=(int)((sym[2]<<RSHIFT)-estim[2])>>(RSHIFT+1);
#endif
#ifdef FIFOVAL
				fifoval_enqueue(sym[2]<<DEPTH*2^sym[1]<<DEPTH^sym[0]);
				fifoval_enqueue(yuv[2]<<DEPTH*2^yuv[1]<<DEPTH^yuv[0]);
#endif
#ifdef USE_RANS

#else
				code[0]=enctable[8*sym[0]+nbypass[0]];
				code[1]=enctable[8*sym[1]+nbypass[1]];
				code[2]=enctable[8*sym[2]+nbypass[2]];
				codelen[0]=(uint8_t)code[0];
				codelen[1]=(uint8_t)code[1];
				codelen[2]=(uint8_t)code[2];
				code[0]>>=8;
				code[1]>>=8;
				code[2]>>=8;
				{
					uint8_t s0=64-codelen[0];
					uint8_t s1=64-codelen[0]-codelen[1];
					uint8_t s2=64-codelen[0]-codelen[1]-codelen[2];

					code[0]=
						code[0]<<s0
					|	code[1]<<s1
					|	code[2]<<s2
					;
					codelen[2]+=codelen[0]+codelen[1];
				}
				codelen[2]+=nbits;
				cache|=code[0]>>nbits;
				if(codelen[2]>=64)
				{
					memcpy(wtptr, &cache, sizeof(uint64_t));
					if(wtptr>=wtend)
					{
						fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
						wtptr-=BUFSIZE;
						memcpy(wtptr, &cache, sizeof(uint64_t));
					}
					wtptr+=sizeof(uint64_t);
					//acme_write(&wtptr, 8, fdst, cache);
					cache=code[0]<<(64-nbits);
					codelen[2]-=64;
					if(!nbits)
						cache=0;
				}
				nbits=codelen[2];
#endif
			}
		}
		memcpy(wtptr, &cache, sizeof(uint64_t));
		if(wtptr>=wtend)
		{
			fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
			wtptr-=BUFSIZE;
			memcpy(wtptr, &cache, sizeof(uint64_t));
		}
		wtptr+=sizeof(uint64_t);
		//acme_write(&wtptr, 8, fdst, cache);
#ifdef MEASURE_PSNR
		{
			static const double mix_rgb[]={1, 1, 1};
			static const double mix_yuv[]={1, 6, 1};
			printf("TRGB:\n");
			print_psnr(esum_rgb, (int64_t)iw*ih, mix_rgb);
			printf("TVYU/CrYCb:\n");
			print_psnr(esum_yuv, (int64_t)iw*ih, mix_yuv);
		}
#endif
	}
	else//dec
	{
		static uint8_t *const rdend=rdbuf+sizeof(uint64_t)+BUFSIZE-sizeof(uint64_t);
		static uint8_t *const wtend=wtbuf+sizeof(uint64_t)+BUFSIZE-3;
		uint64_t cache=0, cache2=0;
		int nbits=0;

		for(int k=0;k<512;++k)
		{
			int val=k;
			val=val>>1^-(val&1);
			packsign[k]=val;
		}
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		fread(&cache, 1, 8, fsrc);
		fread(&cache2, 1, 8, fsrc);
	//	cache	=acme_read(&rdptr, 8, fsrc);
	//	cache2	=acme_read(&rdptr, 8, fsrc);
		for(int ky=0;ky<ih;++ky)
		{
#ifdef USE_ROWS
			int16_t *rows[]=
			{
				pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			};
#endif
			for(int kx=0;kx<iw;++kx)
			{
				uint64_t code, reg;
				int nzeros, prefix, nbypass[3];

#ifdef USE_ROWS
#ifdef USE_CG
				{
					int N[3], W[3], NW[3];
					
					NW[0]	=rows[1][0+(0-1*NCH)*NROWS*NVAL];
					NW[1]	=rows[1][0+(1-1*NCH)*NROWS*NVAL];
					NW[2]	=rows[1][0+(2-1*NCH)*NROWS*NVAL];
					N[0]	=rows[1][0+(0+0*NCH)*NROWS*NVAL];
					N[1]	=rows[1][0+(1+0*NCH)*NROWS*NVAL];
					N[2]	=rows[1][0+(2+0*NCH)*NROWS*NVAL];
					W[0]	=rows[0][0+(0-1*NCH)*NROWS*NVAL];
					W[1]	=rows[0][0+(1-1*NCH)*NROWS*NVAL];
					W[2]	=rows[0][0+(2-1*NCH)*NROWS*NVAL];
					NW[0]=N[0]+W[0]-NW[0];
					NW[1]=N[1]+W[1]-NW[1];
					NW[2]=N[2]+W[2]-NW[2];
					MEDIAN3V_CLOB(pred[0], N[0], W[0], NW[0]);
					MEDIAN3V_CLOB(pred[1], N[1], W[1], NW[1]);
					MEDIAN3V_CLOB(pred[2], N[2], W[2], NW[2]);
				}
#else
				pred[0]=(rows[1][0+(0+0*NCH)*NROWS*NVAL]+rows[0][0+(0-1*NCH)*NROWS*NVAL])>>1;
				pred[1]=(rows[1][0+(1+0*NCH)*NROWS*NVAL]+rows[0][0+(1-1*NCH)*NROWS*NVAL])>>1;
				pred[2]=(rows[1][0+(2+0*NCH)*NROWS*NVAL]+rows[0][0+(2-1*NCH)*NROWS*NVAL])>>1;
#endif
				estim[0]=rows[0][1+(0-1*NCH)*NROWS*NVAL];
				estim[1]=rows[0][1+(1-1*NCH)*NROWS*NVAL];
				estim[2]=rows[0][1+(2-1*NCH)*NROWS*NVAL];
#endif
				nbypass[0]=log2table[estim[0]];
				nbypass[1]=log2table[estim[1]];
				nbypass[2]=log2table[estim[2]];
				if(nbits>=64)
				{
					cache=cache2;
					memcpy(&cache2, rdptr, sizeof(uint64_t));
					if(rdptr>=rdend)
					{
						uint64_t d2;

						fread(rdbuf+sizeof(uint64_t), 1, BUFSIZE, fsrc);
						rdptr-=BUFSIZE;
						memcpy(&d2, rdptr, sizeof(uint64_t));
						cache2|=d2;
					}
					rdptr+=8;
					//cache2=acme_read(&rdptr, 8, fsrc);
					nbits-=64;
				}
				code=cache<<nbits;
				if(nbits)
					code|=cache2>>(64-nbits);

				nzeros=(int)_lzcnt_u64(code);
				prefix=nzeros<<nbypass[0];
				if(nzeros>RLIMIT-1)
					nzeros=RLIMIT-1, nbypass[0]=DEPTH, prefix=0;
				code<<=nzeros+1;
				sym[0]=(int)(code>>(64-nbypass[0]));
				if(!nbypass[0])
					sym[0]=0;
				code<<=nbypass[0];
				sym[0]|=prefix;
				nbits+=nzeros+1+nbypass[0];

				nzeros=(int)_lzcnt_u64(code);
				prefix=nzeros<<nbypass[1];
				if(nzeros>RLIMIT-1)
					nzeros=RLIMIT-1, nbypass[1]=DEPTH, prefix=0;
				code<<=nzeros+1;
				sym[1]=(int)(code>>(64-nbypass[1]));
				if(!nbypass[1])
					sym[1]=0;
				code<<=nbypass[1];
				sym[1]|=prefix;
				nbits+=nzeros+1+nbypass[1];

				nzeros=(int)_lzcnt_u64(code);
				prefix=nzeros<<nbypass[2];
				if(nzeros>RLIMIT-1)
					nzeros=RLIMIT-1, nbypass[2]=DEPTH, prefix=0;
				code<<=nzeros+1;
				sym[2]=(int)(code>>(64-nbypass[2]));
				if(!nbypass[2])
					sym[2]=0;
			//	code<<=nbypass[2];
				sym[2]|=prefix;
#ifdef FIFOVAL
				if(fifoval_check(sym[2]<<DEPTH*2^sym[1]<<DEPTH^sym[0]))
				{
					printf("%016llX\n", cache);
					printf("%016llX\n", cache2);
					CRASH("");
				}
#endif
				nbits+=nzeros+1+nbypass[2];
#ifdef USE_ROWS
				//rows[0][1+(0+0*NCH)*NROWS*NVAL]=(estim[0]+(sym[0]<<RSHIFT))>>2;
				//rows[0][1+(1+0*NCH)*NROWS*NVAL]=(estim[1]+(sym[1]<<RSHIFT))>>2;
				//rows[0][1+(2+0*NCH)*NROWS*NVAL]=(estim[2]+(sym[2]<<RSHIFT))>>2;
				rows[0][1+(0+0*NCH)*NROWS*NVAL]=(2*rows[0][1+(0-1*NCH)*NROWS*NVAL]+(sym[0]<<RSHIFT)+rows[1][1+(0+3*NCH)*NROWS*NVAL])>>2;
				rows[0][1+(1+0*NCH)*NROWS*NVAL]=(2*rows[0][1+(1-1*NCH)*NROWS*NVAL]+(sym[1]<<RSHIFT)+rows[1][1+(1+3*NCH)*NROWS*NVAL])>>2;
				rows[0][1+(2+0*NCH)*NROWS*NVAL]=(2*rows[0][1+(2-1*NCH)*NROWS*NVAL]+(sym[2]<<RSHIFT)+rows[1][1+(2+3*NCH)*NROWS*NVAL])>>2;
#else
				estim[0]+=(int)((sym[0]<<RSHIFT)-estim[0])>>(RSHIFT+1);
				estim[1]+=(int)((sym[1]<<RSHIFT)-estim[1])>>(RSHIFT+1);
				estim[2]+=(int)((sym[2]<<RSHIFT)-estim[2])>>(RSHIFT+1);
#endif
				sym[0]=packsign[sym[0]];
				sym[1]=packsign[sym[1]];
				sym[2]=packsign[sym[2]];
#ifdef NEARLOSSLESS
				sym[0]<<=2;
				sym[1]<<=2;
				sym[2]<<=2;
#endif
				sym[0]+=pred[0];
				sym[1]+=pred[1];
				sym[2]+=pred[2];
				sym[0]<<=32-DEPTH;
				sym[1]<<=32-DEPTH;
				sym[2]<<=32-DEPTH;
				sym[0]>>=32-DEPTH;
				sym[1]>>=32-DEPTH;
				sym[2]>>=32-DEPTH;
#ifdef FIFOVAL
				if(fifoval_check(sym[2]<<DEPTH*2^sym[1]<<DEPTH^sym[0]))
				{
					printf("%016llX\n", cache);
					printf("%016llX\n", cache2);
					CRASH("");
				}
#endif
#ifdef USE_ROWS
				rows[0][0+(0+0*NCH)*NROWS*NVAL]=sym[0];
				rows[0][0+(1+0*NCH)*NROWS*NVAL]=sym[1];
				rows[0][0+(2+0*NCH)*NROWS*NVAL]=sym[2];
				rows[0]+=NCH*NROWS*NVAL;
				rows[1]+=NCH*NROWS*NVAL;
#else
				pred[0]=sym[0];
				pred[1]=sym[1];
				pred[2]=sym[2];
#endif
#ifndef NO_RCT
				sym[2]+=sym[0]>>2;
				sym[1]-=(sym[0]+sym[2])>>2;
				sym[2]+=sym[1];
				sym[0]+=sym[1];
#endif
#ifdef NEARLOSSLESS
				CLAMP2(sym[0], 0, 255);
				CLAMP2(sym[1], 0, 255);
				CLAMP2(sym[2], 0, 255);
#endif
				reg=(uint64_t)sym[2]<<16|(uint64_t)sym[1]<<8|sym[0];
				//sym[1]<<=8;
				//sym[2]<<=16;
				//sym[0]|=sym[1];
				//sym[0]|=sym[2];
			//	reg=sym[0];
				memcpy(wtptr, &reg, sizeof(uint64_t));
				if(wtptr>=wtend)
				{
					fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
					wtptr-=BUFSIZE;
					memcpy(wtptr, &reg, sizeof(uint64_t));
				}
				wtptr+=3;
				//acme_write(&wtptr, 3, fdst, (uint64_t)sym[2]<<16|(uint64_t)sym[1]<<8|sym[0]);
			}
		}
	}
	if(wtptr>wtbuf+sizeof(uint64_t))
		fwrite(wtbuf+sizeof(uint64_t), 1, wtptr-(wtbuf+sizeof(uint64_t)), fdst);
	fclose(fsrc);
	fclose(fdst);
#ifdef USE_ROWS
	free(pixels);
#endif
#ifdef LOUD
	{
		t=time_sec2()-t;
		if(fwd)
		{
			int64_t csize=0;
			struct stat info={0};

			stat(dstfn, &info);
			csize=info.st_size;
			printf("%12lld->%12lld  %12.6lf:1\n", usize, csize, (double)usize/csize);
		}
		printf("%12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB\n"
			, t
			, usize/(t*1024*1024)
			, t*1024*1024*1000/usize
		);
	}
#endif
#ifdef PROFILER
	prof_end(prof_ctx);
#endif
	(void)usize;
	(void)&time_sec2;
	(void)packsignptr;
	return 0;
}
