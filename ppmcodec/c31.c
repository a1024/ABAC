//C31: Cross-platform C29 reference
//standard headers-only
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>
#include<math.h>
#include<sys/stat.h>
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#else
#include<time.h>
#endif
static const char file[]=__FILE__;


#ifdef _MSC_VER
#include"util.h"//
#include<immintrin.h>//
	#define LOUD			//size & time
	#define ESTIMATE_SIZE
	#define ENABLE_GUIDE		//DEBUG		checks interleaved pixels
	#define ANS_VAL			//DEBUG

//	#define C29_VAL			//DESIGN
//	#define PRINTRESIDUALS		//DESIGN

//	#define DISABLE_WG
#else
#define colorprintf(TXTCOLOR, BKCOLOR, FORMAT, ...) printf(FORMAT, ##__VA_ARGS__)
#endif


//3*17+3=54 contexts
#define GRBITS 3
#define NCTX 17

#define XCODERS 4	//xrem 1~3 cols
#define YCODERS 8	//yrem 1~7 rows

#define NCODERS 32

#define PROBBITS 12	//12 bit max	James Bonfield's CDF2sym: {freq<<20 | bias<<8 | sym}

#define RANS_STATE_BITS 31
#define RANS_RENORM_BITS 16

#define WG_NPREDS	8	//multiple of 4

unsigned char* c31_encode(const unsigned char *image, int iw, int ih, size_t *csize, unsigned char **rawptr, int flags);
unsigned char* c31_decode(const unsigned char *data, size_t csize, int *ret_iw, int *ret_ih);

#ifndef ALIGN
#ifdef _MSC_VER
#define	ALIGN(N) __declspec(align(N))
#elif defined __GNUC__
#define	ALIGN(N) __attribute__((aligned(N)))
#endif
#endif
#ifndef AWM_INLINE
#ifdef _MSC_VER
#define AWM_INLINE __forceinline static
#else
#define AWM_INLINE __attribute__((always_inline)) inline static
#endif
#endif
#ifndef MINVAR
#define MINVAR(A, B) ((A)<(B)?(A):(B))
#endif
#ifndef MAXVAR
#define MAXVAR(A, B) ((A)>(B)?(A):(B))
#endif
#ifndef CLAMP2
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
#endif
#ifndef FILLMEM
static void memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
{
	size_t copied;
	char *d=(char*)dst;
	const char *s=(const char*)src;
	//printf("MEMFILL  %016zX %016zX %016zX %016zX\n", (size_t)dst, (size_t)src, dstbytes, srcbytes);//
#ifdef _DEBUG
	if(!dstbytes||!srcbytes)
	{
		//LOG_ERROR("memfill:  dstbytes %zd  srcbytes %zd", dstbytes, srcbytes);
		return;
	}
#endif
	if(dstbytes<srcbytes)
	{
		//printf("  MEMCPY1  %016zX %016zX %016zX\n", (size_t)dst, (size_t)src, dstbytes);//
		memcpy(dst, src, dstbytes);
		return;
	}
	copied=srcbytes;
	//printf("  MEMCPY2  %016zX %016zX %016zX\n", (size_t)d, (size_t)s, srcbytes);//
	memcpy(d, s, srcbytes);
	while((copied<<1)<=dstbytes)
	{
		//printf("  MEMCPY3  %016zX %016zX %016zX\n", (size_t)(d+copied), (size_t)d, copied);//
		memcpy(d+copied, d, copied);
		copied<<=1;
	}
	if(copied<dstbytes)
	{
		//printf("  MEMCPY4  %016zX %016zX %016zX\n", (size_t)(d+copied), (size_t)d, dstbytes-copied);//
		memcpy(d+copied, d, dstbytes-copied);
	}
}
#define FILLMEM(PTR, DATA, ASIZE, ESIZE)\
	do\
	{\
		*(PTR)=(DATA);\
		memfill((PTR)+1, PTR, (ASIZE)-(ESIZE), ESIZE);\
	}while(0)
#endif
#ifndef FLOOR_LOG2
static int floor_log2(unsigned n)
{
	double f=(double)n;
	unsigned long long x=*(unsigned long long*)&f;
	return (x>>52)-1023;
#if 0
	int logn, sh;

	logn=-!n;
	sh=(n>=1ULL<<32)<<5; logn+=sh, n>>=sh;
	sh=(n>=1ULL<<16)<<4; logn+=sh, n>>=sh;
	sh=(n>=1ULL<< 8)<<3; logn+=sh, n>>=sh;
	sh=(n>=1ULL<< 4)<<2; logn+=sh, n>>=sh;
	sh=(n>=1ULL<< 2)<<1; logn+=sh, n>>=sh;
	sh=(n>=1ULL<< 1)<<0; logn+=sh;
	return logn;
#endif
}
#define FLOOR_LOG2 floor_log2
#define FLOOR_LOG2_P1(X) (floor_log2(X)+1)
#endif
#ifndef LOG_ERROR
static int log_error(const char *fn, int line, int quit, const char *format, ...)
{
	ptrdiff_t size=strlen(fn), start=size-1;
	for(;start>=0&&fn[start]!='/'&&fn[start]!='\\';--start);
	start+=start==-1||fn[start]=='/'||fn[start]=='\\';

	printf("\n%s(%d): ", fn+start, line);
	va_list args;
	va_start(args, format);
	vprintf(format, args);
	va_end(args);

	if(quit)
		exit(0);
	return 0;
}
#define LOG_ERROR(MSG, ...) log_error(file, __LINE__, 1, MSG, ##__VA_ARGS__)
#endif
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


#ifdef ANS_VAL
typedef struct _Buffer
{
	size_t count, cap;//in bytes
	unsigned char data[];
} Buffer, *BufferHandle;
static void buffer_realloc(BufferHandle *buf, size_t count, size_t pad)//CANNOT be nullptr, array must be initialized with array_alloc()
{
	size_t dstsize, newcap;

	dstsize=count+pad, newcap=1;
	for(;newcap<dstsize;newcap<<=1);
	if(newcap>buf[0]->cap)
	{
		void *p2=realloc(*buf, sizeof(BufferHandle)+newcap);
		if(!p2)
		{
			LOG_ERROR("Alloc error %016llX", sizeof(BufferHandle)+newcap);
			return;
		}
		*buf=(BufferHandle)p2;
		if(buf[0]->cap<newcap)
			memset(buf[0]->data+buf[0]->cap, 0, newcap-buf[0]->cap);
		buf[0]->cap=newcap;
	}
	buf[0]->count=count;
}
static BufferHandle buffer_construct(const void *src, size_t count, size_t pad)
{
	BufferHandle buf;
	size_t cap;
	
	cap=count+pad;
	buf=(BufferHandle)malloc(sizeof(BufferHandle)+cap);
	if(!buf)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	buf->count=count;
	buf->cap=cap;
	if(src)
		memcpy(buf->data, src, count);
	else
		memset(buf->data, 0, count);
		
	if(cap-count>0)//zero pad
		memset(buf->data+count, 0, cap-count);
	return buf;
}
static void* buffer_insert(BufferHandle *buf, size_t idx, const void *data, size_t count, size_t pad)//cannot be nullptr
{
	size_t start, srcsize, dstsize, movesize;
	
	start=idx;
	srcsize=count;
	dstsize=srcsize;
	movesize=buf[0]->count-start;
	buffer_realloc(buf, buf[0]->count+count, pad);
	memmove(buf[0]->data+start+dstsize, buf[0]->data+start, movesize);
	if(data)
		memfill(buf[0]->data+start, data, dstsize, srcsize);
	else
		memset(buf[0]->data+start, 0, dstsize);
	return buf[0]->data+start;
}
#define BUFFER_APPEND(BUF, DATA, COUNT, PAD) buffer_insert(&(BUF), (BUF)->count, DATA, COUNT, PAD)

#define ANS_VAL_HISTSIZE 128
typedef struct _ANSVALHeader
{
	unsigned short esize, count;
	unsigned idx;
} ANSVALHeader;
static BufferHandle debugstack=0;
static void ansval_push(const void *data, int esize, int count)
{
	static int idx=0;
	ANSVALHeader header={esize, count, idx};
	++idx;
	if(!debugstack)
		debugstack=buffer_construct(0, 0, 1024);
	BUFFER_APPEND(debugstack, &header, sizeof(header), 0);//lo header
	BUFFER_APPEND(debugstack, data, (ptrdiff_t)count*esize, 0);
	BUFFER_APPEND(debugstack, &header, sizeof(header), 0);//hi header
}
static void ansval_printr(const void *data, int esize, int count, const void *xdata)//print elements in reverse because little-endian
{
	const unsigned char *p=(const unsigned char*)data;
	const unsigned char *p2=(const unsigned char*)xdata;
	int size=count*esize;
	for(int k=0;k<size;k+=esize)
	{
		printf(" ");
		for(int k2=esize-1;k2>=0;--k2)
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
		LOG_ERROR("\n");
		return 0;
	}
	return (void*)(nbytes<0?ip2:ip1);
}
static void ansval_check(const void *data, int esize, int count)
{
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
		LOG_ERROR("\n");
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
		LOG_ERROR("");
	}
	++popcount;
}
#endif


//LIFO Bypass Coder
typedef struct _BitPackerLIFO//bwd enc / fwd dec
{
	unsigned long long state;
	int enc_nwritten, dec_navailable;//bitcounts, only for tracking renorms
	unsigned char *dstbwdptr;
	const unsigned char *srcfwdptr, *streamend;
} BitPackerLIFO;
AWM_INLINE void bitpacker_enc_init(BitPackerLIFO *ec, const unsigned char *bufstart, unsigned char *bufptr0_OOB)
{
	memset(ec, 0, sizeof(*ec));
	ec->state=1ULL<<32;
	ec->enc_nwritten=33;
	ec->streamend=bufstart;
	ec->dstbwdptr=bufptr0_OOB;
}
AWM_INLINE void bitpacker_dec_init(BitPackerLIFO *ec, const unsigned char *bufptr0_start, const unsigned char *bufend)
{
	memset(ec, 0, sizeof(*ec));
	ec->srcfwdptr=bufptr0_start+8;
	ec->streamend=bufend;
	ec->state=*(const unsigned long long*)bufptr0_start;
	ec->dec_navailable=FLOOR_LOG2_P1(ec->state);
}
AWM_INLINE void bitpacker_enc_flush(BitPackerLIFO *ec)
{
	ec->dstbwdptr-=8;
#ifdef _DEBUG
	if(ec->dstbwdptr<ec->streamend)
		LOG_ERROR("IntPacker Encoder OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
	*(unsigned long long*)ec->dstbwdptr=ec->state;
}
AWM_INLINE void bitpacker_enc(BitPackerLIFO *ec, int inbits, int sym)
{
#ifdef _DEBUG
	if(!inbits)
		LOG_ERROR("BitPacker inbits=0");
#endif
	//renorm then push inbits
	ec->enc_nwritten+=inbits;
	if(ec->enc_nwritten>64)//renorm on overflow
	{
		ec->enc_nwritten-=32;
		ec->dstbwdptr-=4;
#ifdef _DEBUG
		if(ec->dstbwdptr<ec->streamend)
			LOG_ERROR("IntPacker OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
		*(unsigned*)ec->dstbwdptr=(unsigned)ec->state;
		ec->state>>=32;
#ifdef ANS_VAL
		ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
	}
	ec->state=ec->state<<inbits|sym;
#ifdef ANS_VAL
	ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
}
AWM_INLINE int bitpacker_dec(BitPackerLIFO *ec, int outbits)
{
#ifdef _DEBUG
	if(!outbits)
		LOG_ERROR("BitPacker outbits=0");
#endif
	int sym=ec->state&((1ULL<<outbits)-1);

	//pop outbits then renorm
#ifdef ANS_VAL
	ansval_check(&ec->state, sizeof(ec->state), 1);
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
			LOG_ERROR("IntPacker OOB:  srcfwdptr = 0x%016zX >= 0x%016zX", ec->srcfwdptr, ec->streamend);
#endif
		ec->state=ec->state<<32|*(const unsigned*)ec->srcfwdptr;
		ec->srcfwdptr+=4;
	}
	return sym;
}


#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static void guide_save(const unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(const unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		LOG_ERROR("");
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif
#ifdef PROFILE_SIZE
#include<stdarg.h>
static void profile_size(const unsigned char *dstbwdptr, const char *msg, ...)
{
	static ptrdiff_t size=0;
	static const unsigned char *prev=0;
	if(prev)
	{
		ptrdiff_t diff=prev-dstbwdptr;
		size+=diff;
		printf("%10td (%+10td) bytes", size, diff);
		if(msg)
		{
			va_list args;
			va_start(args, msg);
			printf("  ");
			vprintf(msg, args);
			va_end(args);
		}
		printf("\n");
	}
	prev=dstbwdptr;
}
#else
#define profile_size(...)
#endif
#ifdef PROFILE_TIME
typedef struct _SpeedProfilerInfo
{
	double t;
	ptrdiff_t size;
	const char *msg;
} SpeedProfilerInfo;
#define PROF_CAP 128
static double prof_timestamp=0;
static SpeedProfilerInfo prof_data[PROF_CAP]={0};
static int prof_count=0;
static void prof_checkpoint(ptrdiff_t size, const char *msg)
{
	double t2=time_sec();
	if(prof_timestamp)
	{
		SpeedProfilerInfo *info=prof_data+prof_count++;
		if(prof_count>PROF_CAP)
		{
			LOG_ERROR("Profiler OOB");
			return;
		}
		info->t=t2-prof_timestamp;
		info->size=size;
		info->msg=msg;
	}
	prof_timestamp=t2;
}
static void prof_print(ptrdiff_t usize)
{
	static char buf[2048]={0};
	double timesum=0, tmax=0;
	for(int k=0;k<prof_count;++k)
	{
		double t=prof_data[k].t;
		timesum+=t;
		if(tmax<t)
			tmax=t;
	}
	int printed=0;
	int prev=0;
	double csum=0;
	//printf("| ");
	int colors[128]={0};
	srand((unsigned)__rdtsc());
	colorgen(colors, prof_count, 64, 300, 100);
	//colorgen0(colors, prof_count, 0xC0C0C0);
	printf("1 char = 4 ms\n");
	printf("|");
	for(int k=0;k<prof_count;++k)
	{
		SpeedProfilerInfo *info=prof_data+k;
		csum+=info->t;
		int curr=(int)(csum*250);//fixed scale
		int space=curr-prev;
		int len=0;
		if(info->msg)
			len=(int)strlen(info->msg);
		if(space>2047)//printf("%*s", HUGE, ""); CRASHES
			space=2047;
		if(info->msg&&space>=len)
		{
			int labelstart=(space-len)>>1;
			int labelend=labelstart+len;

			memset(buf, '-', labelstart);
			buf[labelstart]=0;
			colorprintf(colors[k], colors[k], buf);
			colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%s", info->msg);
			memset(buf, '-', (ptrdiff_t)space-labelend);
			buf[space-labelend]=0;
			colorprintf(colors[k], colors[k], buf);
		}
		else
		{
			memset(buf, '-', space);
			buf[space]=0;
			colorprintf(colors[k], colors[k], buf);
		}
		printf("|");
		prev=curr;
		printed+=space+2+(k<prof_count-1);
	}
	printf("\n");
	for(int k=0;k<prof_count;++k)
	{
		SpeedProfilerInfo *info=prof_data+k;
		printf("%16.7lf ms %8.4lf%% ", info->t*1000, 100.*info->t/timesum);
		if(info->size)
			printf(" %12.6lf MB/s %10td bytes ", info->size/(info->t*1024*1024), info->size);
		if(info->msg)
			colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%-20s", info->msg);
		else// if(nstars)
			colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%-20s", "");
		printf("\n");
	}
	printf("\n");
	printf("%16.7lf ms %12.6lf MB/s Total\n", timesum*1000, usize/(timesum*1024*1024));
	printf("\n");
	prof_count=0;
	prof_timestamp=0;
}
#else
#define prof_checkpoint(...)
#define prof_print(...)
#endif
#ifdef C29_VAL
void c29val_init(int iw, int ih);
void c29val_assign(unsigned char pred);
void c29val_check(unsigned char pred);

static int g2_iw=0, g2_ih=0;
static unsigned char *predimage=0;
static ptrdiff_t g2_idx=0, g2_idx2=0;
void c29val_init(int iw, int ih)
{
	ptrdiff_t size=(ptrdiff_t)3*iw*ih;
	g2_iw=iw;
	g2_ih=ih;
	predimage=(unsigned char*)malloc(size);
	if(!predimage)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(predimage, 0, size);
}
void c29val_assign(unsigned char pred)
{
	predimage[g2_idx++]=pred;
}
void c29val_check(unsigned char pred)
{
	if(pred!=predimage[g2_idx2])
	{
		LOG_ERROR("C29val Error  pred %02X vs OG %02X", pred, predimage[g2_idx2]);
		printf("");
	}
	++g2_idx2;
}
#endif


#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
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
//	II_COEFF_U_SUB_V_NBLI,
//	II_COEFF_V_SUB_U_NBLI,

	II_COUNT,
} RCTInfoIdx;
//YUV = RCT * RGB	watch out for permutation in last row
//luma: averaging	chroma: subtraction
//example: _400_40X_3X1 == [1 0 0; -1 0 1; -3/4 1 -1/4]
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
typedef enum _RCTIndex
{
#define RCT(LABEL, ...) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][II_COUNT]=
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

//https://github.com/rygorous/ryg_rans
//https://github.com/samtools/htscodecs
typedef struct _rANS_SIMD_SymInfo	//16 bytes/level	4KB/ctx = 1<<12 bytes
{
	unsigned smax, invf, cdf;
	unsigned short negf, sh;
} rANS_SIMD_SymInfo;
#ifdef ESTIMATE_SIZE
static double g_esize=0;
#endif
static void enc_hist2stats(int *hist, rANS_SIMD_SymInfo *syminfo, unsigned long long *ctxmask, int ctxidx)
{
	int sum=0, count=0;
	for(int ks=0;ks<256;++ks)
	{
		int freq=hist[ks];
		sum+=freq;
		count+=freq!=0;
	}
	int rare=sum<12*256/8;
	*ctxmask|=(unsigned long long)rare<<ctxidx;
#ifdef ESTIMATE_SIZE
	static double esize=0;
	int count0=count, sum0=sum;
#endif
	if(rare)
	{
		for(int ks=0;ks<256;++ks)//bypass
			hist[ks]=1;
		sum=256;
		count=256;
	}
	else if(count==1)//disallow degenerate distribution
	{
		for(int ks=0;ks<256;++ks)
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
	int sum2=0;
	for(int ks=0, ks2=0;ks<256;++ks)//absent symbols get zero freqs
	{
		int freq=hist[ks];
		hist[ks]=(int)(sum2*((1ULL<<PROBBITS)-count)/sum)+ks2;
		ks2+=freq!=0;
		sum2+=freq;
	}
	//for(int ks=0;ks<256;++ks)//never allows zero freqs	INEFFICIENT
	//{
	//	int freq=hist[ks];
	//	hist[ks]=(int)(sum2*((1ULL<<PROBBITS)-256)/sum)+ks;
	//	sum2+=freq;
	//}
#ifdef ESTIMATE_SIZE
	double e=sum0;
	if(count==count0)
	{
		double norm=1./0x1000;
		e=0;
		for(int ks=0;ks<256;++ks)//estimate
		{
			int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
			if(freq)
			{
				double p=freq*norm;
				e-=p*log2(p);
			}
			//	e-=freq*log2(freq*norm);
			//if(e!=e)
			//	LOG_ERROR("");
		}
		e*=sum/8.;
	//	e/=8;
	}
	g_esize+=e;
	if(ctxidx&&!(ctxidx%NCTX))
		printf("\n");
	printf("%c  ctx %3d  %12.2lf / %9d bytes%8.2lf%%  %3d %s",
		ctxidx<3*NCTX?"YUV"[ctxidx/NCTX]:"yuv"[ctxidx-3*NCTX],
		ctxidx, e, sum0, 100.*e/sum0, count0, count==count0?"levels":"bypass"
	);
	if(count==count0&&count<256)
	{
		printf(" %3d", count);
		int fmax=0;
		for(int ks=0;ks<256;++ks)
		{
			int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
			if(fmax<freq)
				fmax=freq;
		}
		for(int ks=0;ks<256;++ks)
		{
			int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
			if(!(ks&15))
				printf(" ");

			int shade=48+freq*(255-48)/fmax;
			colorprintf(shade<<16|shade<<8|shade, freq?0x808080:COLORPRINTF_BK_DEFAULT, "%c", "0123456789ABCDEF"[ks&15]);

			//int shade=freq*255/fmax;
			//colorprintf(freq?0xFFFFFF:0x808080, shade<<16|0<<8|shade, "%c", "0123456789ABCDEF"[ks&15]);

			//printf("%c", freq?"0123456789ABCDEF"[ks&15]:'-');
		}
	}
	printf("\n");
#endif
	int next=1<<PROBBITS;
	for(int ks=255;ks>=0;--ks)
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
			info->sh=FLOOR_LOG2(freq);//eg: x/2 = x*0x80000000>>32>>0
			unsigned long long inv=((0x100000000ULL<<info->sh)+freq-1)/freq;
			info->invf=(unsigned)inv;
			if(inv>0xFFFFFFFF)
			{
				--info->sh;
				info->invf=(unsigned)(inv>>1);
			}
		}
	}
}
static void enc_packhist(BitPackerLIFO *ec, const int *hist, unsigned long long ctxmask, int ctxidx)//histogram must be normalized to PROBBITS, with spike at 128
{
	if(ctxmask>>ctxidx&1)
		return;
	int sum=0;
	unsigned short CDF[257];
	for(int ks=0;ks<256;++ks)//integrage to zigzag CDF to be packed backwards
	{
		int sym=((ks>>1^-(ks&1))+128)&255;
		int freq=hist[sym];
		CDF[ks]=sum;//separate buffer for faster access in 2nd loop
		sum+=freq;
	}
	CDF[256]=1<<PROBBITS;
	
	int cdfW=CDF[0];
	int CDFlevels=1<<PROBBITS;
	int startsym=0;
	for(int ks=1;ks<=256;++ks)//push GR.k
	{
		int next=CDF[ks], freq=next-cdfW;
		int nbypass=FLOOR_LOG2(CDFlevels);
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
	for(int ks=startsym;ks>0;--ks)//encode GR
	{
		int freq=CDF[ks], nbypass=freq>>PROBBITS;
		freq&=(1<<PROBBITS)-1;
		int nzeros=freq>>nbypass, bypass=freq&((1<<nbypass)-1);
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
static void dec_unpackhist(BitPackerLIFO *ec, unsigned *CDF2sym, unsigned long long ctxmask, int ctxidx)
{
	unsigned short hist[257];
	if(ctxmask>>ctxidx&1)//rare context
	{
		for(int ks=0;ks<256;++ks)//bypass
			hist[ks]=(1<<PROBBITS)/256;
	}
	else
	{
		unsigned short CDF[257]={0};
		int CDFlevels=1<<PROBBITS;
		int cdfW=0;
		CDF[0]=0;
		for(int ks=0;ks<256;++ks)//decode GR
		{
			int freq=-1;//stop bit doesn't count
			int nbypass=FLOOR_LOG2(CDFlevels);
			int ks2=ks+1;
			if(ks2>1)
				nbypass-=7;
			if(nbypass<0)
				nbypass=0;
#ifdef ANS_VAL
			ansval_check(&ks2, sizeof(ks2), 1);
#endif
			int bit=0;
			do
			{
				bit=bitpacker_dec(ec, 1);
				++freq;
			}while(!bit);
			if(nbypass)
				freq=freq<<nbypass|bitpacker_dec(ec, nbypass);

			cdfW+=freq;
			CDF[ks]=freq;
			CDFlevels-=freq;
			if(CDFlevels<=0)
			{
#ifdef _DEBUG
				if(CDFlevels<0)
					LOG_ERROR("CDF unpack error");
#endif
				break;
			}
		}
		if(CDFlevels)
			LOG_ERROR("CDF unpack error");
		for(int ks=0;ks<256;++ks)//undo zigzag
		{
			int sym=((ks>>1^-(ks&1))+128)&255;
			hist[sym]=CDF[ks];
		}
	}
	int sum=0;
	for(int ks=0;ks<256;++ks)//integrate
	{
		int freq=hist[ks];
		hist[ks]=sum;
		sum+=freq;
	}
	hist[256]=1<<PROBBITS;
	for(int ks=0;ks<256;++ks)//CDF2sym contains {freq, (state&0xFFF)-cdf, sym}
	{
		int cdf=hist[ks], next=hist[ks+1], freq=next-cdf;
		int val=(freq<<PROBBITS|0)<<8|ks;
		for(int ks2=cdf;ks2<next;++ks2, val+=1<<8)
			CDF2sym[ks2]=val;
	}
}
static void decorr1d(const unsigned char *data, int count, int bytestride, int bestrct, int *rhist, unsigned char **pdstptr)
{
	const unsigned char *combination=rct_combinations[bestrct];
	int yidx=combination[II_PERM_Y];
	int uidx=combination[II_PERM_U];
	int vidx=combination[II_PERM_V];
	int ufromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	int vc0=combination[II_COEFF_V_SUB_Y];
	int vc1=combination[II_COEFF_V_SUB_U];

	const unsigned char *ptr=data;
	unsigned char *dstptr=*pdstptr;
	int prevy=0, prevu=0, prevv=0, offset=0;
	for(int k=0;k<count;++k)
	{
		int y=ptr[yidx]-128;
		int u=ptr[uidx]-128;
		int v=ptr[vidx]-128;
		int sym;
		*dstptr++=sym=(unsigned char)(y-prevy+128);
		++rhist[256*0+sym];
		prevy=y;

		offset=y&ufromy;
		prevu+=offset;
		CLAMP2(prevu, -128, 127);
		*dstptr++=sym=(unsigned char)(u-prevu+128);
		++rhist[256*1+sym];
		prevu=u-offset;

		offset=vc0*y+vc1*u;
		int vpred=(prevv+offset)>>2;
		CLAMP2(vpred, -128, 127);
		*dstptr++=sym=(unsigned char)(v-vpred+128);
		++rhist[256*2+sym];
		prevv=4*v-offset;
		ptr+=bytestride;
	}
	*pdstptr=dstptr;
}
static void encode1d(int count, unsigned *pstate, unsigned char **pstreamptr, const unsigned char *streamend, const rANS_SIMD_SymInfo *rsyminfo, unsigned char **pdptr)
{
	unsigned char *dptr=*pdptr;
	unsigned char *streamptr=*pstreamptr;
	unsigned state=*pstate;
	//unsigned char *ptr=data+(count-(ptrdiff_t)1)*bytestride;
	const rANS_SIMD_SymInfo *info=0;
	for(int k=0;k<count;++k)
	{
		--dptr;
		info=rsyminfo+*dptr+256*2;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)//"streamend" is buffer start
				LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(unsigned short*)streamptr=(unsigned short)state;
			state>>=RANS_RENORM_BITS;
		}
		state+=((unsigned long long)state*info->invf>>32>>info->sh)*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif
		
		--dptr;
		info=rsyminfo+*dptr+256*1;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)
				LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(unsigned short*)streamptr=(unsigned short)state;
			state>>=RANS_RENORM_BITS;
		}
		state+=((unsigned long long)state*info->invf>>32>>info->sh)*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif
		
		--dptr;
		info=rsyminfo+*dptr+256*0;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)
				LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(unsigned short*)streamptr=(unsigned short)state;
			state>>=RANS_RENORM_BITS;
		}
		state+=((unsigned long long)state*info->invf>>32>>info->sh)*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif
		//ptr-=bytestride;
	}
	*pstreamptr=streamptr;
	*pstate=state;
	*pdptr=dptr;
}
static void decode1d(unsigned char *data, int count, int bytestride, int bestrct, unsigned *pstate, const unsigned char **pstreamptr, const unsigned char *streamend, unsigned *rCDF2syms)
{
	const unsigned char *combination=rct_combinations[bestrct];
	int yidx=combination[II_PERM_Y];
	int uidx=combination[II_PERM_U];
	int vidx=combination[II_PERM_V];
	int ufromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	int vc0=combination[II_COEFF_V_SUB_Y];
	int vc1=combination[II_COEFF_V_SUB_U];

	const unsigned char *streamptr=*pstreamptr;
	unsigned state=*pstate;
	unsigned char *ptr=data;
	int prevy=0, prevu=0, prevv=0, offset=0;
	int y=0, u=0, v=0;
	for(int k=0;k<count;++k)
	{
		unsigned info;

		//yuv = (char)(error+N-128)
		info=rCDF2syms[0<<PROBBITS|(state&((1<<PROBBITS)-1))];
		y=(char)(info+prevy-128);
		prevy=y;
#ifdef ANS_VAL
		ansval_check(&state, sizeof(state), 1);
#endif
		state=(state>>PROBBITS)*(info>>(PROBBITS+8))+(info<<PROBBITS>>(32-PROBBITS));
		if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
		{
#ifdef _DEBUG
			if(streamptr>streamend)
				LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(unsigned short*)streamptr;
			streamptr+=2;
		}

		offset=y&ufromy;
		prevu+=offset;
		CLAMP2(prevu, -128, 127);
		info=rCDF2syms[1<<PROBBITS|(state&((1<<PROBBITS)-1))];
		u=(char)(info+prevu-128);
		prevu=u-offset;
#ifdef ANS_VAL
		ansval_check(&state, sizeof(state), 1);
#endif
		state=(state>>PROBBITS)*(info>>(PROBBITS+8))+(info<<PROBBITS>>(32-PROBBITS));
		if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
		{
#ifdef _DEBUG
			if(streamptr>streamend)
				LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(unsigned short*)streamptr;
			streamptr+=2;
		}

		offset=vc0*y+vc1*u;
		int vpred=(prevv+offset)>>2;
		CLAMP2(vpred, -128, 127);
		info=rCDF2syms[2<<PROBBITS|(state&((1<<PROBBITS)-1))];
		v=(char)(info+vpred-128);
		prevv=4*v-offset;
#ifdef ANS_VAL
		ansval_check(&state, sizeof(state), 1);
#endif
		state=(state>>PROBBITS)*(info>>(PROBBITS+8))+(info<<PROBBITS>>(32-PROBBITS));
		if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
		{
#ifdef _DEBUG
			if(streamptr>streamend)
				LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(unsigned short*)streamptr;
			streamptr+=2;
		}
		ptr[yidx]=y+128;
		ptr[uidx]=u+128;
		ptr[vidx]=v+128;
		ptr+=bytestride;
	}
	*pstreamptr=streamptr;
	*pstate=state;
}
static void interleave_blocks_fwd(const unsigned char *original, int iw, int ih, unsigned char *interleaved)
{
	//original[ih][iw][3]
	//interleaved[ih/YCODERS][iw/XCODERS][3][NCODERS]	contiguous & aligned
	//xrem[ih%YCODERS][iw][3]
	//yrem[ih-iw%YCODERS][iw%XCODERS][3]

	//only difference between fwd and inv:		swap assignments (const slow->fast)
	int rowstride=3*iw;
	int ixyblockw=iw/XCODERS;
	int ixyblockh=ih/YCODERS;
	int blockxbytes=3*NCODERS*ixyblockw;
#ifdef INTERLEAVESIMD
	int SIMDxcount=blockxbytes&~((int)sizeof(__m256i[NCODERS])-1);
	__m256i slowinc=_mm256_set1_epi64x(sizeof(__m256i));
#endif
	unsigned char *fastptr=interleaved;
	ALIGN(32) const unsigned char *slowptrs[NCODERS]={0}, *slowptrs0[NCODERS]={0};
	for(int ky=0;ky<YCODERS;++ky)//spread slow pointers
	{
		for(int kx=0;kx<XCODERS;++kx)
			slowptrs0[XCODERS*ky+kx]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);
	}
	for(int ky=0;ky<ixyblockh;++ky)//interleave
	{
		int kx=0;
		memcpy((void*)slowptrs, slowptrs0, sizeof(slowptrs));
#if 1
		for(;kx<blockxbytes;kx+=NCODERS)
		{
			/*
			toy example
			CWH 3*4*5 original
			4 coders
			2*2 blocks
			4 Y's then 4 U's then 4 V's

			original	channel[block][pixel]
			r00 g00 b00 r01 g01 b01  r10 g10 b10 r11 g11 b11
			r02 g02 b02 r03 g03 b03  r12 g12 b12 r13 g13 b13
			r20 g20 b20 r21 g21 b21  r30 g30 b30 r31 g31 b31
			r22 g22 b22 r23 g23 b23  r32 g32 b32 r33 g33 b33

			interleaved
			r00 r10 r20 r30 g00 g10 g20 g30 b00 b10 b20 b30  r01 r11 r21 r31 g01 g11 g21 g31 b01 b11 b21 b31
			r02 r12 r22 r32 g02 g12 g22 g32 b02 b12 b22 b32  r03 r13 r23 r33 g03 g13 g23 g33 b03 b13 b23 b33
			*/
			fastptr[0x00]=*slowptrs[0x00]++;
			fastptr[0x01]=*slowptrs[0x01]++;
			fastptr[0x02]=*slowptrs[0x02]++;
			fastptr[0x03]=*slowptrs[0x03]++;
			fastptr[0x04]=*slowptrs[0x04]++;
			fastptr[0x05]=*slowptrs[0x05]++;
			fastptr[0x06]=*slowptrs[0x06]++;
			fastptr[0x07]=*slowptrs[0x07]++;
			fastptr[0x08]=*slowptrs[0x08]++;
			fastptr[0x09]=*slowptrs[0x09]++;
			fastptr[0x0A]=*slowptrs[0x0A]++;
			fastptr[0x0B]=*slowptrs[0x0B]++;
			fastptr[0x0C]=*slowptrs[0x0C]++;
			fastptr[0x0D]=*slowptrs[0x0D]++;
			fastptr[0x0E]=*slowptrs[0x0E]++;
			fastptr[0x0F]=*slowptrs[0x0F]++;
			fastptr[0x10]=*slowptrs[0x10]++;
			fastptr[0x11]=*slowptrs[0x11]++;
			fastptr[0x12]=*slowptrs[0x12]++;
			fastptr[0x13]=*slowptrs[0x13]++;
			fastptr[0x14]=*slowptrs[0x14]++;
			fastptr[0x15]=*slowptrs[0x15]++;
			fastptr[0x16]=*slowptrs[0x16]++;
			fastptr[0x17]=*slowptrs[0x17]++;
			fastptr[0x18]=*slowptrs[0x18]++;
			fastptr[0x19]=*slowptrs[0x19]++;
			fastptr[0x1A]=*slowptrs[0x1A]++;
			fastptr[0x1B]=*slowptrs[0x1B]++;
			fastptr[0x1C]=*slowptrs[0x1C]++;
			fastptr[0x1D]=*slowptrs[0x1D]++;
			fastptr[0x1E]=*slowptrs[0x1E]++;
			fastptr[0x1F]=*slowptrs[0x1F]++;
			fastptr+=NCODERS;
//#if defined __GNUC__ && !defined PROFILER
//#pragma GCC unroll 32
//#endif
//			for(int k=0;k<NCODERS;++k)
//				*fastptr++=*slowptrs[k]++;
		}
#endif
		for(int k=0;k<NCODERS;++k)
			slowptrs0[k]+=rowstride;
	}
}
static void interleave_blocks_inv(const unsigned char *interleaved, int iw, int ih, unsigned char *original)
{
	//original[ih][iw][3]
	//interleaved[ih/YCODERS][iw/XCODERS][3][NCODERS]	contiguous & aligned
	//xrem[ih%YCODERS][iw][3]
	//yrem[ih-iw%YCODERS][iw%XCODERS][3]

	//only difference between inv and fwd:		swap assignments (slow<-const fast)
	int rowstride=3*iw;
	int ixyblockw=iw/XCODERS;
	int ixyblockh=ih/YCODERS;
	int blockxbytes=3*NCODERS*ixyblockw;
#ifdef INTERLEAVESIMD
	int SIMDxcount=blockxbytes&~((int)sizeof(__m256i[NCODERS])-1);
	__m256i slowinc=_mm256_set1_epi64x(sizeof(__m256i));
#endif
	const unsigned char *fastptr=interleaved;
	unsigned char *slowptrs[NCODERS]={0}, *slowptrs0[NCODERS]={0};
	for(int ky=0;ky<YCODERS;++ky)//spread slow pointers
	{
		for(int kx=0;kx<XCODERS;++kx)
			slowptrs0[XCODERS*ky+kx]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);
	}
	for(int ky=0;ky<ixyblockh;++ky)//deinterleave
	{
		int kx=0;
		memcpy(slowptrs, slowptrs0, sizeof(slowptrs));
#if 1
		for(;kx<blockxbytes;kx+=NCODERS)
		{
			*slowptrs[0x00]++=fastptr[0x00];
			*slowptrs[0x01]++=fastptr[0x01];
			*slowptrs[0x02]++=fastptr[0x02];
			*slowptrs[0x03]++=fastptr[0x03];
			*slowptrs[0x04]++=fastptr[0x04];
			*slowptrs[0x05]++=fastptr[0x05];
			*slowptrs[0x06]++=fastptr[0x06];
			*slowptrs[0x07]++=fastptr[0x07];
			*slowptrs[0x08]++=fastptr[0x08];
			*slowptrs[0x09]++=fastptr[0x09];
			*slowptrs[0x0A]++=fastptr[0x0A];
			*slowptrs[0x0B]++=fastptr[0x0B];
			*slowptrs[0x0C]++=fastptr[0x0C];
			*slowptrs[0x0D]++=fastptr[0x0D];
			*slowptrs[0x0E]++=fastptr[0x0E];
			*slowptrs[0x0F]++=fastptr[0x0F];
			*slowptrs[0x10]++=fastptr[0x10];
			*slowptrs[0x11]++=fastptr[0x11];
			*slowptrs[0x12]++=fastptr[0x12];
			*slowptrs[0x13]++=fastptr[0x13];
			*slowptrs[0x14]++=fastptr[0x14];
			*slowptrs[0x15]++=fastptr[0x15];
			*slowptrs[0x16]++=fastptr[0x16];
			*slowptrs[0x17]++=fastptr[0x17];
			*slowptrs[0x18]++=fastptr[0x18];
			*slowptrs[0x19]++=fastptr[0x19];
			*slowptrs[0x1A]++=fastptr[0x1A];
			*slowptrs[0x1B]++=fastptr[0x1B];
			*slowptrs[0x1C]++=fastptr[0x1C];
			*slowptrs[0x1D]++=fastptr[0x1D];
			*slowptrs[0x1E]++=fastptr[0x1E];
			*slowptrs[0x1F]++=fastptr[0x1F];
			fastptr+=NCODERS;
//#if defined __GNUC__ && !defined PROFILER
//#pragma GCC unroll 32
//#endif
//			for(int k=0;k<NCODERS;++k)
//				*slowptrs[k]++=*fastptr++;
		}
#endif
		for(int k=0;k<NCODERS;++k)
			slowptrs0[k]+=rowstride;
	}
}
static void save_ppm(const char *fn, const unsigned char *image, int iw, int ih)
{
	FILE *fdst=fopen(fn, "wb");
	if(!fdst)
	{
		LOG_ERROR("Cannot open \"%s\" for writing", fn);
		return;
	}
	fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	fwrite(image, 1, (ptrdiff_t)3*iw*ih, fdst);
	fclose(fdst);
}
typedef struct _Context
{
	int iw, ih, blockw, blockh, qxbytes, ixcount, ixbytes, xremw, yremh, xrembytes;
	ptrdiff_t usize, isize, interleavedsize;

	int hsize, *hists;
	int rhsize, *rhist;
	//int ans_permute_size, *ans_permute;
	int CDF2syms_size, *CDF2syms;
	int rCDF2syms_size, *rCDF2syms;

	int psize;
	short *pixels;

	int wgsize;
	short *wgerrors;

	//int wgstatesize;
	//short *wgstate;
} Context;
static int ctx_alloc(Context *ctx, int iw, int ih, int fwd)
{
	memset(ctx, 0, sizeof(*ctx));
	ctx->iw=iw;
	ctx->ih=ih;
	ctx->usize=(ptrdiff_t)3*iw*ih;
	ctx->blockw=iw/XCODERS;
	ctx->blockh=ih/YCODERS;
	ctx->qxbytes=ctx->blockw*XCODERS*3;//iw/XCODERS*XCODERS*3
	ctx->ixcount=ctx->blockw*NCODERS;
	ctx->ixbytes=3*ctx->ixcount;//ix = interleaved circular buffer width		iw/XCODERS*NCODERS
	ctx->xremw=iw-ctx->blockw*XCODERS;
	ctx->yremh=ih-ctx->blockh*YCODERS;
	ctx->xrembytes=3*ctx->xremw;
	ctx->isize=(ptrdiff_t)ctx->ixbytes*ctx->blockh;
	ctx->interleavedsize=ctx->isize<<fwd;//fwd ? interleave residuals & context : pack residuals
	if(fwd)
	{
		ctx->hsize=(int)sizeof(int[3*NCTX<<8]);//3 channels
		ctx->hists=(int*)malloc(ctx->hsize);
		ctx->rhsize=(int)sizeof(int[3*256]);
		ctx->rhist=(int*)malloc(ctx->rhsize);
		if(!ctx->hists||!ctx->rhist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(ctx->hists, 0, ctx->hsize);
		memset(ctx->rhist, 0, ctx->rhsize);
	}
	ctx->CDF2syms_size=(int)sizeof(int[3*NCTX<<PROBBITS]);
	if(fwd)//DIV-free rANS encoder reuses these as SIMD symbol info
		ctx->CDF2syms_size=(int)sizeof(rANS_SIMD_SymInfo[3*NCTX<<8]);
	ctx->CDF2syms=(unsigned*)malloc(ctx->CDF2syms_size);

	ctx->rCDF2syms_size=(int)sizeof(int[3<<PROBBITS]);
	if(fwd)
		ctx->rCDF2syms_size=(int)sizeof(rANS_SIMD_SymInfo[3<<8]);
	ctx->rCDF2syms=(unsigned*)malloc(ctx->rCDF2syms_size);

	//ctx->ans_permute_size=sizeof(char[32][256]);
	//ctx->ans_permute=(int*)malloc(ctx->ans_permute_size);

	ctx->psize=(int)sizeof(short[4*6*NCODERS])*(ctx->blockw+16);//4 padded rows  *  {Y*NCODERS, U*NCODERS, V*NCODERS,  eY*NCODERS, eU*NCODERS, eV*NCODERS} = 2*3*32 = 192 channels  ~48*iw bytes
	ctx->pixels=(short*)malloc(ctx->psize);//~188 KB for 4K/12MP
	ctx->wgsize=(int)sizeof(short[2*NCODERS*3*WG_NPREDS])*(ctx->blockw+16);//2 padded rows  *  {WGY*NCODERS, WGU*NCODERS, WGV*NCODERS} * NPREDS = 3*32*8 = 768 channels  ~96*iw bytes
	ctx->wgerrors=(short*)malloc(ctx->wgsize);//~375 KB for 4K/12MP		NNEerrors = currerrors
	//ctx->wgstatesize=(int)sizeof(short[2*NCODERS*3*WG_NPREDS]);//{preds, Wprederrors} * NPREDS * 3 channels * NCODERS
	//ctx->wgstate=(short*)malloc(ctx->wgstatesize);
	if(!ctx->CDF2syms||!ctx->rCDF2syms||!ctx->pixels||!ctx->wgerrors)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	return 0;
}
static void ctx_free(Context *ctx, int fwd)
{
	if(fwd)
	{
		free(ctx->hists);
		free(ctx->rhist);
	}
	free(ctx->CDF2syms);
	free(ctx->rCDF2syms);
	//free(ctx->ans_permute);
	free(ctx->pixels);
	free(ctx->wgerrors);
	//free(ctx->wgstate);
}
static void c31_main(Context *ctx, unsigned char *interleaved, int rct, int use_wg4, int fwd, const unsigned char **pstreamptr, const unsigned char *streamend)
{
	const unsigned char *streamptr=0;
	const unsigned char *combination=rct_combinations[rct];
	//int
	//	yidx=combination[II_PERM_Y]*NCODERS,
	//	uidx=combination[II_PERM_U]*NCODERS,
	//	vidx=combination[II_PERM_V]*NCODERS;
	int uhelpmask=-(combination[II_COEFF_U_SUB_Y]!=0);
	int vc0=combination[II_COEFF_V_SUB_Y];
	int vc1=combination[II_COEFF_V_SUB_U];
	int paddedwidth=ctx->blockw+16;
	unsigned char *imptr=interleaved+(fwd?ctx->isize:0);
	unsigned short *ctxptr=(unsigned short*)interleaved;
	unsigned state[NCODERS]={0};
	if(!fwd)
	{
		streamptr=*pstreamptr;
#ifdef _DEBUG
		if(streamptr>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		memcpy(state, streamptr, sizeof(state));
		streamptr+=sizeof(state);
	}
	memset(ctx->pixels, 0, ctx->psize);
	if(use_wg4)
		memset(ctx->wgerrors, 0, ctx->wgsize);
	for(int ky=0;ky<ctx->blockh;++ky)//main coding loop
	{
		ALIGN(32) short *rows[]=
		{
			ctx->pixels+(paddedwidth*((ky-0LL)&3)+8LL)*6*NCODERS,
			ctx->pixels+(paddedwidth*((ky-1LL)&3)+8LL)*6*NCODERS,
			ctx->pixels+(paddedwidth*((ky-2LL)&3)+8LL)*6*NCODERS,
			ctx->pixels+(paddedwidth*((ky-3LL)&3)+8LL)*6*NCODERS,
		};
		ALIGN(32) short *erows[]=
		{
			ctx->wgerrors+(paddedwidth*((ky-0LL)&1)+8LL)*NCODERS*3*WG_NPREDS,
			ctx->wgerrors+(paddedwidth*((ky-1LL)&1)+8LL)*NCODERS*3*WG_NPREDS,
		};
		int pred=0;
		int wgpreds[WG_NPREDS]={0};
		int offsets[NCODERS]={0};
		char yuv[3][NCODERS]={0};
		int val=0, sym=0, error=0;
		for(int kx=0;kx<ctx->ixbytes;kx+=3*NCODERS)
		{
			for(int kc=0;kc<3;++kc)
			{
				for(int kl=0;kl<NCODERS;++kl)
				{
					//context = FLOOR_LOG2(eW*eW+1)
					int eW=rows[0][(kc+3-1*6)*NCODERS+kl];
					int ctx2=FLOOR_LOG2(eW*eW+1);
					if(ctx2>NCTX-1)
						ctx2=NCTX-1;
					//	rows[-Y][(kc+E+X*6)*NCODERS+kl]
					int NW	=rows[1][(kc+0-1*6)*NCODERS+kl];
					int N	=rows[1][(kc+0+0*6)*NCODERS+kl];
					int W	=rows[0][(kc+0-1*6)*NCODERS+kl];
					int vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					pred=N+W-NW;
					if(use_wg4)
					{
						int NNN		=rows[3][(kc+0+0*6)*NCODERS+kl];
						int NN		=rows[2][(kc+0+0*6)*NCODERS+kl];
						int NNE		=rows[2][(kc+0+1*6)*NCODERS+kl];
						int NE		=rows[1][(kc+0+1*6)*NCODERS+kl];
						int NEE		=rows[1][(kc+0+2*6)*NCODERS+kl];
						int NEEE	=rows[1][(kc+0+3*6)*NCODERS+kl];
						int NEEEE	=rows[1][(kc+0+4*6)*NCODERS+kl];
						int WWWW	=rows[0][(kc+0-4*6)*NCODERS+kl];
						int WWW		=rows[0][(kc+0-3*6)*NCODERS+kl];
						int WW		=rows[0][(kc+0-2*6)*NCODERS+kl];

						wgpreds[0]=N;
						wgpreds[1]=W;
						wgpreds[2]=3*(N-NN)+NNN;
						wgpreds[3]=3*(W-WW)+WWW;
						wgpreds[4]=W+NE-N;
						wgpreds[5]=(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)>>2;
						wgpreds[6]=pred;
						wgpreds[7]=N+NE-NNE;

						//if(fwd&&ky==0&&kx<4*3*NCODERS&&kl==0)//
						if(fwd&&ky==0&&kx==12000&&kc==0&&kl==24)//
							printf("");

						float fpred=0, csum=0;
						//errors = [pI]+pNW+2*pN+pNE+pNNE	(__m256i*)(erows[-Y]+((P+X*NPREDS)*3+C)*NCODERS)+R
						for(int kp=0;kp<WG_NPREDS;++kp)
						{
							static const float weights[]=
							{
								240,
								240,
								180,
								180,
								140,
								160,
								120,
								120,
							};
							int e=1//(__m256i*)(erows[-Y]+((P+X*NPREDS)*3+C)*NCODERS)+R
								+erows[1][((kp+0*WG_NPREDS)*3+kc)*NCODERS+kl]*2	//eN*2
								+erows[1][((kp-1*WG_NPREDS)*3+kc)*NCODERS+kl]	//eNW
								+erows[1][((kp+1*WG_NPREDS)*3+kc)*NCODERS+kl]	//eNE
								+erows[0][((kp+1*WG_NPREDS)*3+kc)*NCODERS+kl]	//eNNE
							;
							float coeff=1.f/(float)e;
							*(unsigned*)&coeff&=0xFFFF8000;//use 9 bit like in C29		FIXME emulate with integers and FLOOR_LOG2
							float fwp=(float)wgpreds[kp];
							coeff*=weights[kp];
							fpred+=coeff*fwp;
							csum+=coeff;

							if(fwd&&ky==0&&kx==12000&&kc==0&&kl==24)//
								printf("0x%08X 0x%08X  %20.15lf %20.15lf %20.15lf\n",
									*(unsigned*)&coeff,
									*(unsigned*)&fwp,
									coeff,
									coeff*fwp,
									fwp
								);
						}
						fpred/=csum;
						pred=(int)round(fpred);
						
						if(vmin>NW)vmin=NW;
						if(vmax<NW)vmax=NW;
						if(vmin>NE)vmin=NE;
						if(vmax<NE)vmax=NE;
						if(vmin>NEEE)vmin=NEEE;
						if(vmax<NEEE)vmax=NEEE;
					}
					CLAMP2(pred, vmin, vmax);
					if(kc)
					{
						pred+=offsets[kl];
						if(kc==2)
							pred>>=2;
						CLAMP2(pred, -128, 127);
					}

					//if(fwd&&ky==0&&kx<4*3*NCODERS&&kl==0)//
					//	printf("");

					if(fwd)
					{
#ifdef C29_VAL
						c29val_check(pred);
#endif
						val=imptr[combination[II_PERM_Y+kc]*NCODERS+kl];
						sym=(unsigned char)(val-pred);
						int cell=(NCTX*kc+ctx2)<<8|sym;
						*ctxptr++=cell;
						//if((unsigned)cell>=(unsigned)ctx->hsize)//
						//	LOG_ERROR("");
						++ctx->hists[cell];
#ifdef PRINTRESIDUALS
						if(!ky&&!kl)
							printf("%02X %02X %04X%s", val, pred+128, cell, kc==2?"\n":"  ");
#endif
					}
					else
					{
#ifdef ANS_VAL
						ansval_check(&state[kl], sizeof(int), 1);
#endif
						//update
						int rem=state[kl]&((1<<PROBBITS)-1);
						unsigned info=ctx->CDF2syms[(NCTX*kc+ctx2)<<PROBBITS|rem];//James Bonfield's CDF2sym: {freq<<20 | bias<<8 | sym}
						int freq=info>>20, bias=info>>8&((1<<PROBBITS)-1);//1 <= freq <= 0xF01
						sym=info&255;
						state[kl]=(state[kl]>>PROBBITS)*freq+bias;//state = (state>>12)*freq+(rem-cdf)	rem-cdf is prebaked

#ifdef ANS_VAL
						ansval_check(&state[kl], sizeof(int), 1);
#endif
						//renorm
						if(state[kl]<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
						{
#ifdef _DEBUG
							if(streamptr>=streamend)
								LOG_ERROR("OOB ptr %016zX >= start %016zX", streamptr, streamend);
#endif
							state[kl]=state[kl]<<16|*(unsigned short*)streamptr;
							streamptr+=2;
						}

						val=(unsigned char)(sym+pred);
						imptr[combination[II_PERM_Y+kc]*NCODERS+kl]=val;
					}
					val-=128;
					yuv[kc][kl]=val;
					error=val-pred;
					error=error<<1^error>>31;
					if(kc==2)
						val<<=2;
					val-=offsets[kl];
					//	(kc+E+X*6)*NCODERS+kl
					rows[0][(kc+0+0*6)*NCODERS+kl]=val;
					//context update = (2*eW+(e<<3)+max(eNEE, eNEEE))>>2
					rows[0][(kc+3+0*6)*NCODERS+kl]=(
						+2*eW
						+(error<<3)
						+MAXVAR(rows[1][(kc+3+2*6)*NCODERS+kl], rows[1][(kc+3+3*6)*NCODERS+kl])
					)>>2;
					if(use_wg4)
					{
						for(int kp=0;kp<WG_NPREDS;++kp)
						{
							int p=abs(val-wgpreds[kp])<<1;

							//curr = (eW + min(eN, eW) + e + min(eN, eNEE))>>2
							int pN	=erows[1][((kp+0*WG_NPREDS)*3+kc)*NCODERS+kl];
							int pNEE=erows[1][((kp+2*WG_NPREDS)*3+kc)*NCODERS+kl];
							int pW	=erows[0][((kp-1*WG_NPREDS)*3+kc)*NCODERS+kl];
							//if(eN<0||eNEE<0||eW<0||e<0)//
							//	LOG_ERROR("");
							erows[0][((kp+0*WG_NPREDS)*3+kc)*NCODERS+kl]=(
								+pW	+MINVAR(pN, pW)
								+p	+MINVAR(pN, pNEE)
							)>>2;
							
							//eNE += e
							erows[1][((kp+1*WG_NPREDS)*3+kc)*NCODERS+kl]+=p;
							//if(erows[1][((kp+1*WG_NPREDS)*3+kc)*NCODERS+kl]<0)//
							//	LOG_ERROR("");
						}
					}
					if(!kc)
						offsets[kl]=yuv[0][kl]&uhelpmask;
					else if(kc==1)
						offsets[kl]=vc0*yuv[0][kl]+vc1*yuv[1][kl];
					else
						offsets[kl]=0;
				}
			}
			imptr+=3*NCODERS;
			rows[0]+=6*NCODERS;
			rows[1]+=6*NCODERS;
			rows[2]+=6*NCODERS;
			rows[3]+=6*NCODERS;
			erows[0]+=3*NCODERS*WG_NPREDS;
			erows[1]+=3*NCODERS*WG_NPREDS;
		}
	}
	if(!fwd)
		*pstreamptr=streamptr;
}
unsigned char* c31_encode(const unsigned char *image, int iw, int ih, size_t *csize, unsigned char **rawptr, int flags)
{
	Context ctx={0};
	ctx_alloc(&ctx, iw, ih, 1);
	unsigned char *interleaved=(unsigned char*)malloc(ctx.interleavedsize);
	if(!interleaved)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	interleave_blocks_fwd(image, iw, ih, interleaved+ctx.isize);
	guide_save(interleaved+ctx.isize, ctx.ixcount, ctx.blockh);
	prof_checkpoint(ctx.usize, "interleave");
	int bestrct=0, use_wg4=0;
	{
		long long counters[OCH_COUNT]={0};
		unsigned char *imptr=interleaved+ctx.isize;
		for(int ky=0;ky<ctx.blockh;++ky)//analysis
		{
			__m256i mprev[OCH_COUNT*2];//16-bit
			memset(mprev, 0, sizeof(mprev));//

			int prev[OCH_COUNT][NCODERS]={0};
			for(int kx=0;kx<ctx.ixbytes-3*NCODERS;kx+=3*NCODERS)
			{
				for(int kl=0;kl<NCODERS;++kl)
				{
					//if(ky>=ctx.blockh/2&&kx>=ctx.ixbytes/2)//
					//	printf("");

					int r=imptr[kl+0*NCODERS]-128;
					int g=imptr[kl+1*NCODERS]-128;
					int b=imptr[kl+2*NCODERS]-128;
					r<<=2;
					g<<=2;
					b<<=2;
					int rg=r-g;
					int gb=g-b;
					int br=b-r;
#define UPDATE(IDX0, IDX1, IDX2, EXPR0, EXPR1, EXPR2)\
	do\
	{\
		int curr0=EXPR0;\
		int curr1=EXPR1;\
		int curr2=EXPR2;\
		counters[IDX0]+=abs(curr0-prev[IDX0][kl]);\
		counters[IDX1]+=abs(curr1-prev[IDX1][kl]);\
		counters[IDX2]+=abs(curr2-prev[IDX2][kl]);\
		prev[IDX0][kl]=curr0;\
		prev[IDX1][kl]=curr1;\
		prev[IDX2][kl]=curr2;\
	}while(0)
					UPDATE(OCH_Y400, OCH_Y040, OCH_Y004, r, g, b);
					UPDATE(OCH_CX40, OCH_C0X4, OCH_C40X, rg, gb, br);
					UPDATE(OCH_CX31, OCH_C3X1, OCH_C31X, rg+(gb>>2), rg+(br>>2), br+(rg>>2));//r-(3*g+b)/4,  g-(3*r+b)/4,  b-(3*r+g)/4
					UPDATE(OCH_CX13, OCH_C1X3, OCH_C13X, br+(gb>>2), gb+(br>>2), gb+(rg>>2));//r-(g+3*b)/4,  g-(r+3*b)/4,  b-(r+3*g)/4
					UPDATE(OCH_CX22, OCH_C2X2, OCH_C22X, (rg-br)>>1, (gb-rg)>>1, (br-gb)>>1);//r-(g+b)/2,  g-(r+b)/2,  b-(r+g)/2
#undef  UPDATE
				}
				imptr+=3*NCODERS;
			}
		}
		long long minerr=0;
		for(int kt=0;kt<RCT_COUNT;++kt)
		{
			const unsigned char *rct=rct_combinations[kt];
			long long currerr=
				+counters[rct[0]]
				+counters[rct[1]]
				+counters[rct[2]]
			;
//#ifdef LOUD
//			printf("%-14s %12lld + %12lld + %12lld = %12lld%s\n",
//				rct_names[kt],
//				counters[rct[0]],
//				counters[rct[1]],
//				counters[rct[2]],
//				currerr,
//				!kt||minerr>currerr?" <-":""
//			);
//#endif
			if(!kt||minerr>currerr)
			{
				minerr=currerr;
				bestrct=kt;
			}
		}
		const unsigned char *rct=rct_combinations[bestrct];
		use_wg4=counters[rct[0]]+counters[rct[1]]+counters[rct[2]] > ctx.isize*2;//FIXME tune
		if(flags&1)//force CG
			use_wg4=0;
		else if(flags&2)//force WG4
			use_wg4=1;
#ifdef DISABLE_WG
		use_wg4=0;//DISABLE WG4
#endif
		prof_checkpoint(ctx.usize, "analysis");
#ifdef LOUD
		printf("%s  WG4=%d  %td bytes\n", rct_names[bestrct], use_wg4, (ptrdiff_t)3*iw*ih);
#endif
		//bestrct=0;//MARKER
	}

	//generate encode permutations		eg: mask = MSB 0b11000101 LSB  ->  LO {x, x, x, x, 0, 2, 6, 7} HI
#if 0
	for(int km=0;km<256;++km)
	{
		int *curr=ctx.ans_permute+((ptrdiff_t)km<<3);
		int kb2=7;
		for(int kb=7;kb>=0;--kb)
		{
			int bit=km>>kb&1;
			if(bit)
				curr[kb2--]=kb;
		}
	}
	prof_checkpoint(ctx.ans_permute_size, "gen permutation");
#endif

	c31_main(&ctx, interleaved, bestrct, use_wg4, 1, 0, 0);
	prof_checkpoint(ctx.isize, "main");
	
	rANS_SIMD_SymInfo *syminfo=(rANS_SIMD_SymInfo*)ctx.CDF2syms;
	rANS_SIMD_SymInfo *rsyminfo=(rANS_SIMD_SymInfo*)ctx.rCDF2syms;
	
	//normalize/integrate hists
	unsigned long long ctxmask=0;
	for(int kc=0;kc<3*NCTX;++kc)
		enc_hist2stats(ctx.hists+(ptrdiff_t)256*kc, syminfo+(ptrdiff_t)256*kc, &ctxmask, kc);
	ptrdiff_t cap=(ptrdiff_t)4*iw*ih;
	unsigned char *data=(unsigned char*)malloc(cap);
	if(!data)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	unsigned char *streamptr=data+cap;
	if(ctx.xremw||ctx.yremh)//encode remainder
	{
		int rembufsize=3*(iw*ctx.yremh+ih*ctx.xremw);
		unsigned char *rembuf=(unsigned char*)malloc(rembufsize);
		if(!rembuf)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		unsigned char *dptr=rembuf;
		int rowstride=3*ctx.iw;
		memset(ctx.rhist, 0, ctx.rhsize);
		for(int ky=0;ky<ctx.yremh;++ky)
			decorr1d(image+rowstride*(ctx.blockh*YCODERS+ky), iw, 3, bestrct, ctx.rhist, &dptr);
		for(int kx=0;kx<ctx.xremw;++kx)
			decorr1d(image+ctx.qxbytes+3*kx, ctx.blockh*YCODERS, rowstride, bestrct, ctx.rhist, &dptr);
		enc_hist2stats(ctx.rhist+(ptrdiff_t)256*0, rsyminfo+(ptrdiff_t)256*0, &ctxmask, 3*NCTX+0);
		enc_hist2stats(ctx.rhist+(ptrdiff_t)256*1, rsyminfo+(ptrdiff_t)256*1, &ctxmask, 3*NCTX+1);
		enc_hist2stats(ctx.rhist+(ptrdiff_t)256*2, rsyminfo+(ptrdiff_t)256*2, &ctxmask, 3*NCTX+2);
#ifdef ESTIMATE_SIZE
		printf("%12.2lf/%8td bytes\n", g_esize, (ptrdiff_t)3*iw*ih);
#endif
		unsigned state=1<<(RANS_STATE_BITS-RANS_RENORM_BITS);
		for(int kx=ctx.xremw-1;kx>=0;--kx)
			encode1d(ctx.blockh*YCODERS, &state, &streamptr, image, rsyminfo, &dptr);
		//	encode1d(image+ctx.qxbytes+3*kx, ctx.blockh*YCODERS, rowstride, &state, &streamptr, image, rsyminfo);
		for(int ky=ctx.yremh-1;ky>=0;--ky)
			encode1d(iw, &state, &streamptr, image, rsyminfo, &dptr);
		//	encode1d(image+rowstride*(ctx.blockh*YCODERS+ky), iw, 3, &state, &streamptr, image, rsyminfo);
		free(rembuf);
		//flush
		streamptr-=4;
#ifdef _DEBUG
		if(streamptr<=data)
			LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
		*(unsigned*)streamptr=state;
		prof_checkpoint(ctx.usize-ctx.isize, "encode remainder");
		profile_size(streamptr, "/ %9td bytes remainders", usize-isize);
	}
	unsigned state[NCODERS];
	for(int k=0;k<NCODERS;++k)
		state[k]=1<<(RANS_STATE_BITS-RANS_RENORM_BITS);
	unsigned short *ctxptr2=(unsigned short*)(interleaved+(ctx.isize<<1)-sizeof(short[NCODERS]));
	for(int ky=ctx.blockh-1;ky>=0;--ky)
	{
		for(int kx=ctx.ixbytes-NCODERS;kx>=0;kx-=NCODERS)//ixbytes = iw/XCODERS*NCODERS*3
		{
			if(ky==ctx.blockh/2)//
				printf("");
			for(int kl=NCODERS-1;kl>=0;--kl)
			{
				rANS_SIMD_SymInfo *info=syminfo+ctxptr2[kl];

				//enc renorm
				if(state[kl]>info->smax)
				{
					streamptr-=2;
#ifdef _DEBUG
					if(streamptr<data)
						LOG_ERROR("OOB ptr %016zX < start %016zX", streamptr, data);
#endif
					*(unsigned short*)streamptr=state[kl];
				}
#ifdef ANS_VAL
				ansval_push(&state[kl], sizeof(int), 1);
#endif
				//enc update		state += (state*invf>>32>>sh)*negf+cdf		state = state/freq<<12|(cdf+state%freq)
				state[kl]+=(unsigned)((unsigned long long)state[kl]*info->invf>>32>>info->sh)*info->negf+info->cdf;
#ifdef ANS_VAL
				ansval_push(&state[kl], sizeof(int), 1);
#endif
			}
			ctxptr2-=NCODERS;
		}
	}
	streamptr-=sizeof(state);
#ifdef _DEBUG
	if(streamptr<=data)
		LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, data);
#endif
	memcpy(streamptr, state, sizeof(state));
	prof_checkpoint(ctx.isize, "encode main");
	profile_size(streamptr, "/ %9td bytes main", isize);
	
	//pack hists
	{
		BitPackerLIFO ec;
		bitpacker_enc_init(&ec, image, streamptr);
		if(ctx.xremw||ctx.yremh)
		{
			enc_packhist(&ec, ctx.rhist+256*2, ctxmask, 3*NCTX+2);
			enc_packhist(&ec, ctx.rhist+256*1, ctxmask, 3*NCTX+1);
			enc_packhist(&ec, ctx.rhist+256*0, ctxmask, 3*NCTX+0);
		}
		for(int kc=3*NCTX-1;kc>=0;--kc)
			enc_packhist(&ec, ctx.hists+(ptrdiff_t)256*kc, ctxmask, kc);
		bitpacker_enc_flush(&ec);
		streamptr=ec.dstbwdptr;

		streamptr-=8;
		*(unsigned long long*)streamptr=ctxmask;
	}
	prof_checkpoint(((ptrdiff_t)3*NCTX+((ctx.xremw||ctx.yremh)?3:0))*256, "pack histograms");
	profile_size(streamptr, "/ %9d bytes overhead", (3*NCTX+3)*12<<8>>3);

	streamptr-=4;
	*(int*)streamptr=ih;
	streamptr-=4;
	*(int*)streamptr=iw;
	streamptr-=2;
	*(short*)streamptr='2'|'9'<<8;

	free(interleaved);
	ctx_free(&ctx, 1);
	*csize=data+cap-streamptr;
	*rawptr=data;
	return streamptr;
}
unsigned char* c31_decode(const unsigned char *data, size_t csize, int *ret_iw, int *ret_ih)
{
	const unsigned char *streamptr=data, *streamend=data+csize;
	int iw, ih, rct, use_wg4;
	unsigned long long ctxmask;
	{
		int flags;

		if(*(unsigned short*)streamptr!=('2'|'9'<<8))
		{
			LOG_ERROR("Decode error");
			return 0;
		}
		streamptr+=2;

		iw=*(int*)streamptr;
		streamptr+=4;

		ih=*(int*)streamptr;
		streamptr+=4;

		flags=*streamptr++;
		use_wg4=flags&1;
		rct=flags>>1;

		ctxmask=*(unsigned long long*)streamptr;
		streamptr+=8;
	}
	Context ctx={0};
	ctx_alloc(&ctx, iw, ih, 0);
	unsigned char *image=malloc(ctx.usize+16);
	if(!image)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	{
		BitPackerLIFO ec;
		bitpacker_dec_init(&ec, streamptr, streamend);
		for(int kc=0;kc<3*NCTX;++kc)
			dec_unpackhist(&ec, ctx.CDF2syms+((ptrdiff_t)kc<<PROBBITS), ctxmask, kc);
		if(ctx.xremw||ctx.yremh)
		{
			dec_unpackhist(&ec, ctx.rCDF2syms+((ptrdiff_t)0<<PROBBITS), ctxmask, 3*NCTX+0);
			dec_unpackhist(&ec, ctx.rCDF2syms+((ptrdiff_t)1<<PROBBITS), ctxmask, 3*NCTX+1);
			dec_unpackhist(&ec, ctx.rCDF2syms+((ptrdiff_t)2<<PROBBITS), ctxmask, 3*NCTX+2);
		}
		streamptr=(unsigned char*)(size_t)ec.srcfwdptr;
		prof_checkpoint((ptrdiff_t)ctx.CDF2syms_size+ctx.rCDF2syms_size, "unpack histograms");
	}

	//generate decode permutations		eg: mask = MSB 0b11000101 LSB  ->  LO {0, x, 1, x, x, x, 2, 3} HI
#if 0
	for(int km=0;km<256;++km)
	{
		int *curr=ctx.ans_permute+((ptrdiff_t)km<<3);
		int idx=0;
		for(int kb=0;kb<8;++kb)
		{
			int bit=km>>kb&1;
			if(bit)
				curr[kb]=idx++;
		}
	}
	prof_checkpoint(ctx.ans_permute_size, "gen permutation");
#endif
	unsigned char *interleaved=(unsigned char*)malloc(ctx.interleavedsize);
	if(!interleaved)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	c31_main(&ctx, interleaved, rct, use_wg4, 0, &streamptr, streamend);
	prof_checkpoint(ctx.isize, "main");
	
	interleave_blocks_inv(interleaved, iw, ih, image);
	prof_checkpoint(usize, "deinterleave");

	if(ctx.xremw||ctx.yremh)
	{
		int rowstride=3*ctx.iw;
#ifdef _DEBUG
		if(streamptr>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		unsigned state=*(unsigned*)streamptr;
		streamptr+=4;
		for(int ky=0;ky<ctx.yremh;++ky)
			decode1d(image+rowstride*(ctx.blockh*YCODERS+ky), iw, 3, rct, &state, (const unsigned char**)&streamptr, streamend, ctx.rCDF2syms);
		for(int kx=0;kx<ctx.xremw;++kx)
			decode1d(image+ctx.qxbytes+3*kx, ctx.blockh*YCODERS, rowstride, rct, &state, (const unsigned char**)&streamptr, streamend, ctx.rCDF2syms);
		prof_checkpoint(usize-isize, "remainder");
	}

	free(interleaved);
	ctx_free(&ctx, 0);
	*ret_iw=iw;
	*ret_ih=ih;
	return image;
}
int c31_codec(const char *srcfn, const char *dstfn, int nthreads0)
{
#ifdef LOUD
	double t=time_sec();
#endif
	ptrdiff_t srcsize=0;
	unsigned char *srcbuf=0;
	{
		struct stat info={0};
		int e2=stat(srcfn, &info);
		if(e2)
		{
			printf("File not found  \"%s\"", srcfn);
			return 1;
		}
#if defined _MSC_VER || defined _WIN32
		if((info.st_mode&S_IFMT)==S_IFREG)
#else
		if(S_ISREG(info.st_mode))
#endif
			srcsize=info.st_size;
		else
		{
			printf("Not a regular file  \"%s\"", srcfn);
			return 1;
		}
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		srcbuf=(unsigned char*)malloc(srcsize+16);
		if(!srcbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		ptrdiff_t nread=fread(srcbuf, 1, srcsize, fsrc);
		if(nread!=srcsize)
			printf("File size %td, but %td bytes were read.\n", srcsize, nread);
		fclose(fsrc);
	}
	int iw=0, ih=0, fwd=0;
	unsigned char *srcptr=srcbuf;
	unsigned char *image=0, *data=0;
	ptrdiff_t csize=0;
	{
		int tag=*(unsigned short*)srcptr;
		srcptr+=2;
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('2'|'9'<<8))
		{
			printf("Unsupported file  \"%s\"", srcfn);
			free(srcbuf);
			return 1;
		}
	}
	if(fwd)
	{
#ifdef LOUD
		print_timestamp("%Y-%m-%d_%H%M%S\n");
#endif
		if(*srcptr++!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
		while(*srcptr=='#')
		{
			++srcptr;//skip '#'
			while(*srcptr!='\n')++srcptr;
			++srcptr;//skip '\n'
		}
		iw=0;
		while((unsigned)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++-'0';
		while(*srcptr<'!')++srcptr;//skip space
		ih=0;
		while((unsigned)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++-'0';
		if(memcmp(srcptr, "\n255\n", 5))
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
		srcptr+=5;
		image=srcptr;
		unsigned char *rawptr=0;
		data=c31_encode(image, iw, ih, &csize, &rawptr, nthreads0);

		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			LOG_ERROR("Cannot open \"%s\" for writing", fdst);
			return 1;
		}
		fwrite(data, 1, csize, fdst);
		fclose(fdst);
		free(rawptr);
#ifdef LOUD
		printf("%s  WH %d*%d\n", srcfn, iw, ih);
		printf("%8td/%8td bytes\n", csize, (ptrdiff_t)3*iw*ih);
#endif
	}
	else
	{
		image=c31_decode(srcbuf, srcsize, &iw, &ih);
		save_ppm(dstfn, image, iw, ih);
		free(image);
	}
#ifdef LOUD
	t=time_sec()-t;
	printf("%c  %12.6lf sec  %12.6lf MB/s\n", 'D'+fwd, t, (ptrdiff_t)3*iw*ih/(t*1024*1024));
#endif
	return 0;
}