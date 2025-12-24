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


#ifdef _MSC_VER
	#define LOUD
	#define ENABLE_GUIDE
#endif


#define GRBITS 6
#define LRSHIFT 4


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
static int g_iw=0, g_ih=0;
static uint8_t *g_im1=0, *g_im2=0;
static double g_sqe[3]={0};
static uint8_t* guide_save(const uint8_t *image, int iw, int ih)
{
	uint8_t *im2=0;
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	im2=(uint8_t*)malloc(size);
	if(!im2)
	{
		CRASH("Alloc error");
		return 0;
	}
	memcpy(im2, image, size);
	return im2;
}
static void guide_check(const uint8_t *image, const uint8_t *im0, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, im0+idx, 3))
	{
		printf("\n\nGuide error  X %5d  Y %5d  0x%02X%02X%02X != 0x%02X%02X%02X\n\n"
			, kx
			, ky
			, image[idx+0]
			, image[idx+1]
			, image[idx+2]
			, im0[idx+0]
			, im0[idx+1]
			, im0[idx+2]
		);
		CRASH("Guide error");
		printf("\n");//trick for old debuggers
	}
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

#define LZMIN 6
#define LZLENBITS 9
#define LZMAX (1<<LZLENBITS)
#define LZBACKBITS 24

#define EBITS 14
#define ESIZE (1<<EBITS)
typedef struct _ETable
{
	int32_t etable[ESIZE], estart, eend, ecount;
} ETable;
static ETable tables[0x100];
int c42_codec(int argc, char **argv)
{
	const char *srcfn=argv[1], *dstfn=argv[2];
	int fwd=0, iw=0, ih=0;
	ptrdiff_t usize=0, csize=0, headersize=0, cap=0;
	uint8_t *buf=0, *image=0, *streamstart=0, *streamend=0;
	RiceCoder ec;
	uint8_t *imptr=0, *imend=0;
	int prevdata=0, prevsym=128;
#ifdef ENABLE_GUIDE
	static uint8_t *im0=0;
#endif
#ifdef LOUD
	double t=0, t2=0;
#endif

	if(argc!=3)
	{
		printf(
			"Usage:  \"%s\"  input  output    To encode/decode.\n"
			"Built on %s %s\n"
			, argv[0]
			, __DATE__, __TIME__
		);
		return 1;
	}
#ifdef LOUD
	t=time_sec();
#endif
	{
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			CRASH("Cannot open \"%s\"", srcfn);
			return 1;
		}
		int64_t c=0;
		fread(&c, 1, 2, fsrc);
		fwd=c==('P'|'6'<<8);
		if(!fwd&&c!=('4'|'2'<<8))
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
			CRASH("Unsupported source file");
			return 1;
		}
		headersize=ftell(fsrc);
		//rowstride=3*iw;
		usize=(ptrdiff_t)3*iw*ih;
		cap=(ptrdiff_t)7*iw*ih;
		buf=(uint8_t*)malloc(cap);
		if(!buf)
		{
			CRASH("Alloc error");
			return 1;
		}
		if(fwd)
		{
			image=buf+cap-usize-sizeof(uint64_t);
			streamstart=buf;
			streamend=buf+cap;

			fread(image, 1, usize, fsrc);//read image
#ifdef ENABLE_GUIDE
			g_im1=guide_save(image, iw, ih);
#endif
		}
		else
		{
			struct stat info={0};
			stat(srcfn, &info);
			csize=info.st_size;

			image=buf;
			streamstart=buf+cap-csize-sizeof(uint64_t);
			streamend=buf+cap;

			fread(streamstart, 1, csize-headersize, fsrc);//read stream
		}
		fclose(fsrc);
	}
	rice_init(&ec, streamstart, streamend);
	imptr=image;
	imend=image+usize;
	if(fwd)
	{
		const int queuesize=sizeof(uint8_t[LZMAX]);
		uint8_t *queue=(uint8_t*)malloc(queuesize);
		int qcount=0, qstart=0, qend=0;
//#ifdef LOUD
//		int64_t avback=0, avlen=0, codecount=0;
//#endif

		if(!queue)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(tables, 0, sizeof(tables));
		memset(queue, 0, queuesize);
		*imend=0;
		while(imptr<imend)
		{
			int matchidx=-1, matchlen=0;

			uint32_t sym0=*(uint8_t*)imptr;
		//	uint32_t sym0=*(uint16_t*)imptr;
			ETable *table=tables+sym0;
			int ctr=table->ecount;
			while(ctr)
			{
				int tidx=(table->estart+ctr)%ESIZE;
				int idx=table->etable[tidx], len=0;
				uint8_t *search1=image+idx;
				uint8_t *search2=imptr;
				ptrdiff_t searchend=image+usize-imptr-sizeof(uint64_t);
				while(len<searchend&&*(uint64_t*)(search1+len)==*(uint64_t*)(search2+len))
					len+=sizeof(uint64_t);
				searchend=image+usize-imptr-1;
				while(len<searchend&&search1[len]==search2[len])
					++len;
				if(len>matchlen)
					matchidx=idx, matchlen=len;

				//check for other previous encounters
				--ctr;
			}
			//add current position to table
			if(table->ecount>=ESIZE)//table is full
			{
				table->etable[table->eend]=(int32_t)(imptr-image);
				table->eend=(table->eend+1)%ESIZE;
				table->estart=(table->estart+1)%ESIZE;
			}
			else
			{
				table->etable[table->eend]=(int32_t)(imptr-image);
				table->eend=(table->eend+1)%ESIZE;
				++table->ecount;
			}
			if(qcount>=LZMAX)//queue is full
			{
				rice_enc(&ec, 0, 1);
				rice_enc(&ec, LZLENBITS, qcount-1);
				while(qcount)
				{
					int data=queue[qstart];
					int sym=data-prevdata;
					prevdata=data;
					sym=sym<<1^sym>>31;
					rice_enc(&ec, FLOOR_LOG2((prevsym>>GRBITS)+1), sym);
					qstart=(qstart+1)%LZMAX;
					--qcount;
					prevsym+=((sym<<GRBITS)-prevsym)>>LRSHIFT;
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
						int data=queue[qstart];
						int sym=data-prevdata;
						prevdata=data;
						sym=sym<<1^sym>>31;
						rice_enc(&ec, FLOOR_LOG2((prevsym>>GRBITS)+1), sym);
						qstart=(qstart+1)%LZMAX;
						--qcount;
						prevsym+=((sym<<GRBITS)-prevsym)>>LRSHIFT;
					}
				}
				rice_enc(&ec, 0, 0);
				rice_enc(&ec, LZLENBITS, matchlen-1);
				rice_enc(&ec, LZBACKBITS, (int)(imptr-image-matchidx));
				
#if defined LOUD && 1
#if 0
				static int ctr=0;
				++ctr;
				if(ctr<<(32-9))
				{
					printf("nentries %10d back %10d  len %10d\n"
						, nentries
						, (int)(imptr-image)-matchidx
						, matchlen
					);
				}
#endif
				//avback+=(int)(imptr-image)-matchidx;
				//avlen+=matchlen;
				//++codecount;
#endif

				//jump
				imptr+=matchlen;
			}
			else//append byte to stream
			{
				//enqueue symbol
				queue[qend]=sym0;
			//	queue[qend]=sym0&255;
				qend=(qend+1)%LZMAX;
				++qcount;

				++imptr;
			}
		}
		if(qcount)
		{
			rice_enc(&ec, 0, 1);
			rice_enc(&ec, LZLENBITS, qcount-1);
			while(qcount)
			{
				int data=queue[qstart];
				int sym=data-prevdata;
				prevdata=data;
				sym=sym<<1^sym>>31;
				rice_enc(&ec, FLOOR_LOG2((prevsym>>GRBITS)+1), sym);
				qstart=(qstart+1)%LZMAX;
				--qcount;
				prevsym+=((sym<<GRBITS)-prevsym)>>LRSHIFT;
			}
			qstart=qend;
		}
//#ifdef LOUD
//		printf("\n");
//		printf("avback %12.2lf  avlen %12.2lf\n"
//			, (double)avback/codecount
//			, (double)avlen/codecount
//		);
//#endif
		//free(table);
		free(queue);
		rice_enc_flush(&ec);
	}
	else
	{
		for(;imptr<imend;)
		{
			int mode=rice_dec(&ec, 0);
			int count=rice_dec(&ec, LZLENBITS)+1;
			uint8_t *fillend=imptr+count;
			if(fillend>imend)//guard
				fillend=imend;
			if(mode)//symbols
			{
				while(imptr<fillend)
				{
					int sym=rice_dec(&ec, FLOOR_LOG2((prevsym>>GRBITS)+1));
					int data=sym>>1^-(sym&1);
					prevsym+=((sym<<GRBITS)-prevsym)>>LRSHIFT;
					data+=prevdata;
					prevdata=data;
					*imptr++=data;
#ifdef ENABLE_GUIDE
					if(imptr[-1]!=g_im1[imptr-image-1])
					{
						printf("Guide error  IDX %10d  0x%02X != 0x%02X\n"
							, (int)(imptr-image-1)
							, imptr[-1]
							, g_im1[imptr-image-1]
						);
						CRASH("Guide error");
					}
#endif
				}
			}
			else//copy
			{
				int backtrack=rice_dec(&ec, LZBACKBITS);
				uint8_t *src=imptr-backtrack;
				while(imptr<fillend)
				{
					*imptr++=*src++;
#ifdef ENABLE_GUIDE
					if(imptr[-1]!=g_im1[imptr-image-1])
					{
						printf("Guide error  IDX %10d  0x%02X != 0x%02X\n"
							, (int)(imptr-image-1)
							, imptr[-1]
							, g_im1[imptr-image-1]
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
			csize+=fwrite(&iw, 1, 3, fdst);
			csize+=fwrite(&ih, 1, 3, fdst);
			csize+=fwrite(streamstart, 1, ec.ptr-streamstart, fdst);
		}
		else
		{
			headersize=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(image, 1, usize, fdst);
		}
		fclose(fdst);
	}
	free(buf);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		printf("Mem usage %lld bytes\n", sizeof(tables));
		usize+=headersize;
		printf("WH %d*%d  \"%s\"\n", iw, ih, srcfn);
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
