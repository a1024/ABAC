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


#ifdef _MSC_VER
	#define LOUD
	#define ENABLE_GUIDE
#endif

#define GRBITS 6


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
static ptrdiff_t g_size=0;
static uint8_t *g_buf=0;
static uint8_t* guide_save(const uint8_t *buf, ptrdiff_t size)
{
	uint8_t *buf2=0;

	g_size=size;
	buf2=(uint8_t*)malloc(size);
	if(!buf2)
	{
		CRASH("Alloc error");
		return 0;
	}
	memcpy(buf2, buf, size);
	return buf2;
}
#endif
#endif


typedef struct _RiceCoder
{
	uint64_t cache, nbits;
	uint8_t *ptr, *end;
} RiceCoder;
AWM_INLINE void rice_init(RiceCoder *ec, uint8_t *start, uint8_t *end)
{
	ec->cache=0;
	ec->nbits=64;
	ec->ptr=start;
	ec->end=end;
}
AWM_INLINE void rice_enc_flush(RiceCoder *ec)
{
	*(uint64_t*)ec->ptr=ec->cache;
	ec->ptr+=8;
}
AWM_INLINE void rice_enc(RiceCoder *ec, int nbypass, int sym)
{
	//buffer: {c,c,c,b,b,a,a,a, f,f,f,e,e,e,d,c}, cache: MSB gg[hhh]000 LSB	nbits is number of ASSIGNED bits
	//written 64-bit words are byte-reversed because the CPU is little-endian

	int nzeros=sym>>nbypass, bypass=sym&((1<<nbypass)-1);
	if(nzeros>=ec->nbits)//fill the rest of cache with zeros, and flush
	{
		nzeros-=(int)ec->nbits;
		*(uint64_t*)ec->ptr=ec->cache;
		ec->ptr+=8;
		ec->cache=0;
		while(nzeros>=64)//just flush zeros
		{
			nzeros-=64;
			*(uint64_t*)ec->ptr=ec->cache;
			ec->ptr+=8;
		}
		ec->nbits=64;
	}
	//now there is room for zeros:  0 <= nzeros < nbits <= 64
	ec->nbits-=nzeros;//emit remaining zeros to cache

	bypass|=1<<nbypass;//append 1 stop bit
	++nbypass;
	if(nbypass>=ec->nbits)//not enough free bits in cache:  fill cache, write to list, and repeat
	{
		nbypass-=(int)ec->nbits;
		ec->cache|=(uint64_t)bypass>>nbypass;
		bypass&=(1<<nbypass)-1;
		*(uint64_t*)ec->ptr=ec->cache;
		ec->ptr+=8;
		ec->cache=0;
		ec->nbits=64;
	}
	//now there is room for bypass:  0 <= nbypass < nbits <= 64
	if(nbypass)
	{
		ec->nbits-=nbypass;//emit remaining bypass to cache
		ec->cache|=(uint64_t)bypass<<ec->nbits;
	}
}
AWM_INLINE int rice_dec(RiceCoder *ec, int nbypass)
{
	//cache: MSB 00[hhh]ijj LSB	nbits is number of CLEARED bits (past codes must be cleared from cache)
	
	int sym;

	if(ec->cache)
	{
		uint64_t lz=_lzcnt_u64(ec->cache);
		sym=(int)(lz-ec->nbits);
		ec->nbits=lz;
	}
	else
	{
		sym=64-(int)ec->nbits;
		for(;;)
		{
			ec->cache=*(uint64_t*)ec->ptr;
			ec->ptr+=8;
			if(ec->cache)
			{
				ec->nbits=_lzcnt_u64(ec->cache);
				sym+=(int)ec->nbits;
				break;
			}
			sym+=64;
		}
	}
	sym<<=nbypass;
	ec->cache&=0x7FFFFFFFFFFFFFFF>>ec->nbits;//remove stop bit
	ec->nbits+=(uint64_t)nbypass+1;
	if(ec->nbits>=64)//nbits = nbits0+nbypass > N
	{
		//example: 000000[11 1]1010010	nbits=6, nbypass=3	6+3-8 = 1
		ec->nbits-=64;
		sym|=(int)(ec->cache<<ec->nbits);
		ec->cache=*(uint64_t*)ec->ptr;
		ec->ptr+=8;
		if(ec->nbits)
		{
			sym|=(int)(ec->cache>>(64-ec->nbits));
			ec->cache&=0xFFFFFFFFFFFFFFFF>>ec->nbits;//nbits=61 -> cache&=7;

			//ec->cache&=(1ULL<<sh)-1;//sh=3 -> cache&=7;
		}
		return sym;
	}
	if(nbypass)
	{
		sym|=(int)(ec->cache>>(64-ec->nbits));
		ec->cache&=0xFFFFFFFFFFFFFFFF>>ec->nbits;
	}
	return sym;
}

#define LZMIN 5
#define LZLENBITS 9
#define LZMAX (1<<LZLENBITS)
#define LZWBITS 11
#define LZWSIZE (1<<LZWBITS)

#define EBITS 18
#define ESIZE (1<<EBITS)
typedef struct _ETable
{
	int32_t etable[ESIZE], estart, eend, ecount;
} ETable;
static ETable tables[256];
int c43_codec(int argc, char **argv)
{
	const uint16_t tag='4'|'3'<<8;
	if(argc!=4||(uint32_t)((argv[1][0]&0xDF)-'D')>=2||argv[1][1])
	{
		printf(
			"Usage:  \"%s\"  e|d  input  output    To encode/decode.\n"
			"Built on %s %s\n"
			, argv[0]
			, __DATE__, __TIME__
		);
		return 1;
	}
	int fwd=(argv[1][0]&0xDF)=='E';
	const char *srcfn=argv[2], *dstfn=argv[3];
#ifdef LOUD
	double t=time_sec(), t2=0;
#endif
	ptrdiff_t usize=0, csize=0, cap=0;
	uint8_t *buf=0, *bufptr=0, *bufend=0;
	uint8_t *streamstart=0, *streamend=0;
	RiceCoder ec;
	const int nbypass=7;
#ifdef ENABLE_GUIDE
	static uint8_t *im0=0;
#endif

#ifdef LOUD
	t=time_sec();
#endif
	{
		struct stat info={0};
		int error=stat(srcfn, &info);
		if(error||!info.st_size)
		{
			CRASH("Cannot stat \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
			usize=info.st_size;
		else
			csize=info.st_size;
	}
	{
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			CRASH("Cannot open \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
		{
			cap=usize*4/3;
		}
		else
		{
			int h=0;

			cap=csize;
			fread(&h, 1, 2, fsrc);
			if(h!=tag)
			{
				CRASH("Unupported file \"%s\"", srcfn);
				return 1;
			}
			fread(&usize, 1, 4, fsrc);
		}
		buf=(uint8_t*)malloc(usize);
		streamstart=(uint8_t*)malloc(cap);
		if(!buf||!streamstart)
		{
			CRASH("Alloc error");
			return 1;
		}
		streamend=buf+cap;
		if(fwd)
		{

			fread(buf, 1, usize, fsrc);//read image
#ifdef ENABLE_GUIDE
			g_buf=guide_save(buf, usize);
#endif
		}
		else
		{
			fread(streamstart, 1, csize-ftell(fsrc), fsrc);//read stream
		}
		fclose(fsrc);
	}
	rice_init(&ec, streamstart, streamend);
	bufptr=buf;
	bufend=buf+usize;
	if(fwd)
	{
		const int queuesize=(int)sizeof(uint8_t[LZMAX]);
		uint8_t *queue=(uint8_t*)malloc(queuesize);
		int qcount=0, qstart=0, qend=0;
		if(!queue)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(tables, 0, sizeof(tables));
		memset(queue, 0, queuesize);
		while(bufptr<bufend)
		{
			uint8_t sym=*bufptr;
			int matchidx=-1, matchlen=0;
			
			ETable *table=tables+sym;
			int ctr=table->ecount;
			while(ctr)
			{
				int tidx=(table->estart+ctr)%ESIZE;
				int idx=table->etable[tidx], len=0;
				while(bufptr+len<bufend&&bufptr[len]==buf[idx+len])
					++len;
				if(len>matchlen)
					matchidx=idx, matchlen=len;

				//check for other previous encounters
				--ctr;
			}
			//add current position to table
			if(table->ecount>=ESIZE)//table is full
			{
				table->etable[table->eend]=(int32_t)(bufptr-buf);
				table->eend=(table->eend+1)%ESIZE;
				table->estart=(table->estart+1)%ESIZE;
			}
			else
			{
				table->etable[table->eend]=(int32_t)(bufptr-buf);
				table->eend=(table->eend+1)%ESIZE;
				++table->ecount;
			}
#if 0
			for(int k=1;k<LZWSIZE;++k)//search the window
			{
				uint8_t *src=bufptr-k;
				int len;

				if(src<buf)
					break;
				if(*src!=sym)
					continue;
				for(len=1;len<LZMAX;++len)
				{
					if(src[len]!=bufptr[len])
						break;
				}
				if(len>matchlen)
				{
					matchlen=len;
					matchidx=(int32_t)(src-buf);
				}
			}
#endif
			if(qcount>=LZMAX)//queue is full
			{
				rice_enc(&ec, 0, 1);
				rice_enc(&ec, LZLENBITS, qcount-1);
				while(qcount)
				{
					uint8_t sym=queue[qstart];
					rice_enc(&ec, nbypass, sym);
					qstart=(qstart+1)%LZMAX;
					--qcount;
				}
				qstart=qend;
			}
			if(matchlen>=LZMIN)//emit match
			{
				if(qcount)
				{
					rice_enc(&ec, 0, 1);
					rice_enc(&ec, LZLENBITS, qcount-1);
					while(qcount)
					{
						int sym=queue[qstart];
						rice_enc(&ec, nbypass, sym);
						qstart=(qstart+1)%LZMAX;
						--qcount;
					}
				}
				rice_enc(&ec, 0, 0);
				rice_enc(&ec, LZLENBITS, matchlen-1);
				rice_enc(&ec, LZWBITS, (int)(bufptr-buf-matchidx));

				//jump
				bufptr+=matchlen;
			}
			else//append byte to stream
			{
				//enqueue symbol
				queue[qend]=sym;
				qend=(qend+1)%LZMAX;
				++qcount;

				++bufptr;
			}
#ifdef LOUD
			static int printctr=0;
			if(!((printctr+1)<<(32-16)))
				printf("\rIDX %12lld  CR %8.4lf%%"
					, bufptr-buf
					, 100.*(ec.ptr-streamstart)/(bufptr-buf)
				);
			++printctr;
#endif
		}
#ifdef LOUD
		printf("\n");
#endif
		if(qcount)
		{
			rice_enc(&ec, 0, 1);
			rice_enc(&ec, LZLENBITS, qcount-1);
			while(qcount)
			{
				uint8_t sym=queue[qstart];
				rice_enc(&ec, nbypass, sym);
				qstart=(qstart+1)%LZMAX;
				--qcount;
			}
			qstart=qend;
		}
		free(queue);
		rice_enc_flush(&ec);
	}
	else
	{
		for(;bufptr<bufend;)
		{
			int mode=rice_dec(&ec, 0);
			int count=rice_dec(&ec, LZLENBITS)+1;
			uint8_t *fillend=bufptr+count;
			if(fillend>bufend)//guard
				fillend=bufend;
			if(mode)//symbols
			{
				while(bufptr<fillend)
				{
					int sym=rice_dec(&ec, 7);
					*bufptr++=sym;
#ifdef ENABLE_GUIDE
					if(bufptr[-1]!=g_buf[bufptr-buf-1])
					{
						printf("Guide error  IDX %10d  0x%02X != 0x%02X\n"
							, (int)(bufptr-buf-1)
							, bufptr[-1]
							, g_buf[bufptr-buf-1]
						);
						CRASH("Guide error");
					}
#endif
				}
			}
			else//copy
			{
				int backtrack=rice_dec(&ec, LZWBITS);
				uint8_t *src=bufptr-backtrack;
				while(bufptr<fillend)
				{
					*bufptr++=*src++;
#ifdef ENABLE_GUIDE
					if(bufptr[-1]!=g_buf[bufptr-buf-1])
					{
						printf("Guide error  IDX %10d  0x%02X != 0x%02X\n"
							, (int)(bufptr-buf-1)
							, bufptr[-1]
							, g_buf[bufptr-buf-1]
						);
						CRASH("Guide error");
					}
#endif
				}
			}
		}
	}
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing\n", dstfn);
			free(buf);
			return 1;
		}
		if(fwd)
		{
			csize=0;
			csize+=fwrite("42", 1, 2, fdst);
			csize+=fwrite(&usize, 1, 4, fdst);
			csize+=fwrite(streamstart, 1, ec.ptr-streamstart, fdst);
		}
		else
			fwrite(buf, 1, usize, fdst);
		fclose(fdst);
	}
	free(buf);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		printf("%10lld bytes  \"%s\"\n", (int64_t)usize, srcfn);
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
	(void)time_sec;
	return 0;
}
