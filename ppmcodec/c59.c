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
#include<immintrin.h>
#ifdef PROFILER
#include"util.h"
#endif


#ifdef _MSC_VER
	#define LOUD
//	#define PROFILE_SIZE
//	#define RICE_OPT
	#define ANALYSIS2

	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif

	#define USE_ROWS


enum
{
	DEPTH_LUMA=8,
	DEPTH_CHROMA=9,
	RLIMIT=23,
#ifdef USE_ROWS
	RSHIFT=3,
#else
	RSHIFT=1,
#endif

	BUFSIZE=512<<10,

	NCTX=32,
#ifdef USE_ROWS
	XPAD=8,
	NCH=3,
	NROWS=2,
	NVAL=2,
#endif
};


//runtime
#if 1

//clobbers A B C
#define MEDIAN3V_CLOB(M, A, B, C)\
	M=A, A=A<B?A:B, B=B<M?B:M,\
	M=B, B=C<B?C:B, C=T>C?T:C,\
	M=A, A=B<A?B:A, M=T>B?T:B

#ifndef CLAMP2
#define CLAMP2(X, LO, HI) X=X>LO?X:LO, X=X<HI?X:HI
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
#ifndef FLOOR_LOG2
#if defined _M_X64 || defined __x86_64__
#define FLOOR_LOG2(X)\
	(sizeof(X)==8?63-(int32_t)_lzcnt_u64(X):31-_lzcnt_u32((uint32_t)(X)))
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
#if 0
INLINE int median3i(int a, int b, int c)
{
	int t;

	t=a;
	a=a>b?b:a;
	b=t>b?t:b;

	t=b;
	b=b>c?c:b;
	c=t>c?t:c;

	t=a;
	a=a>b?b:a;
	b=t>b?t:b;

	return b;
}
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
#endif//ENABLE_GUIDE
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
static uint8_t rdbuf[BUFSIZE+sizeof(uint64_t)], wtbuf[BUFSIZE+sizeof(uint64_t)];
INLINE uint64_t acme_read(uint8_t **pptr, ptrdiff_t size, FILE *f)
{
	uint8_t *ptr=*pptr;
	uint64_t data=*(uint64_t*)ptr;
	ptrdiff_t left;

	/*
	overflow:
	|                    ______left______   ______right_____
	|                   /                \ /                \
	|buf1start ... ... [datastart  buf1end|buf2start  dataend] ...
	|                   \________________    _______________/
	|                                    size
	*/
	
	left=rdbuf+BUFSIZE-ptr;
	ptr+=size;
	if(left<size)
	{
		fread(rdbuf, 1, BUFSIZE, f);
		ptr=(rdbuf+size)-left;
		left<<=3;
		data&=0xFFFFFFFFFFFFFFFF>>(64-left);
		data|=*(uint64_t*)rdbuf<<left;
	}
	*pptr=ptr;
	return data;
}
INLINE void acme_write(uint8_t **pptr, ptrdiff_t size, FILE *f, uint64_t data)
{
	uint8_t *ptr=*pptr;
	ptrdiff_t left;
	
	/*
	overflow:
	|                    ______left______   ______right_____
	|                   /                \ /                \
	|buf1start ... ... [datastart  buf1end|buf2start  dataend] ...
	|                   \________________    _______________/
	|                                    size
	*/
	
	*(uint64_t*)ptr=data;
	left=wtbuf+BUFSIZE-ptr;
	ptr+=size;
	if(left<size)
	{
		fwrite(wtbuf, 1, BUFSIZE, f);
		ptr=(wtbuf+size)-left;
		left<<=3;
		data>>=left;
		*(uint64_t*)wtbuf=data;
	}
	*pptr=ptr;
}
typedef struct _RiceCoder
{
	uint64_t cache, nbits;
	uint8_t *ptr;
	FILE *f;
} RiceCoder;
INLINE void rice_init(RiceCoder *ec, uint8_t *ptr, FILE *f)
{
	ec->cache=0;
	ec->nbits=64;
	ec->ptr=ptr;
	ec->f=f;
}
INLINE void rice_flush(RiceCoder *ec)
{
	acme_write(&ec->ptr, 8, ec->f, ec->cache);
}
INLINE void rice_enc(RiceCoder *ec, int nbypass, int sym)
{
	//buffer: {c,c,c,b,b,a,a,a, f,f,f,e,e,e,d,c}, cache: MSB gg[hhh]000 LSB	nbits is number of ASSIGNED bits
	//written 64-bit words are byte-reversed because the CPU is little-endian

	int nzeros=sym>>nbypass;
	int bypass=sym&0x7FFFFFFF>>(31-nbypass);
//	int bypass=sym&((1<<nbypass)-1);
#ifdef ESTIMATE_SIZES
	bsizes[g_kc][0]+=nzeros+1;
	bsizes[g_kc][1]+=nbypass;
#endif
	if(nzeros>=ec->nbits)//fill the rest of cache with zeros, and flush
	{
		nzeros-=(int)ec->nbits;
		acme_write(&ec->ptr, 8, ec->f, ec->cache);
		ec->cache=0;
		while(nzeros>=64)//just flush zeros
		{
			nzeros-=64;
			acme_write(&ec->ptr, 8, ec->f, 0);
		}
		ec->nbits=64;
	}
	//now there is room for zeros:  0 <= nzeros < nbits <= 64
	ec->nbits-=nzeros;//emit remaining zeros to cache

	bypass|=1<<nbypass;//append 1 stop bit
	++nbypass;
	if(nbypass>=ec->nbits)//cache would overflow:  fill, flush, and repeat
	{
		nbypass-=(int)ec->nbits;
		ec->cache|=(uint64_t)bypass>>nbypass;
		bypass&=0x7FFFFFFF>>(31-nbypass);
		acme_write(&ec->ptr, 8, ec->f, ec->cache);
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
INLINE int rice_dec(RiceCoder *ec, int nbypass)
{
	//cache: MSB 00[hhh]ijj LSB	nbits is number of CLEARED bits (past codes must be cleared from cache)
	
	int sym;

	sym=-(int)ec->nbits;
	while(!ec->cache)
	{
		sym+=64;
		ec->cache=acme_read(&ec->ptr, 8, ec->f);
	}
	ec->nbits=_lzcnt_u64(ec->cache);
	sym+=(int)ec->nbits;

	sym<<=nbypass;
	ec->cache&=0x7FFFFFFFFFFFFFFF>>ec->nbits;//remove stop bit
	ec->nbits+=(uint64_t)nbypass+1;
	if(ec->nbits>=64)//nbits = nbits0+nbypass > N
	{
		//example: 000000[11 1]1010010	nbits=6, nbypass=3	6+3-8 = 1
		ec->nbits-=64;
		sym|=(int)(ec->cache<<ec->nbits);
		ec->cache=acme_read(&ec->ptr, 8, ec->f);
		nbypass=(int)ec->nbits;
	}
	if(nbypass)
	{
		sym|=(int)(ec->cache>>(64-ec->nbits));
		ec->cache&=0xFFFFFFFFFFFFFFFF>>ec->nbits;//nbits=61 -> cache&=7;
	}
	return sym;
}

static int16_t packsign[1024], *const packsignptr=packsign+512;
int c59_codec(int argc, char **argv)
{
	const uint16_t tag='5'|'9'<<8;
	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	int64_t c=0;
	int fwd=0, iw=0, ih=0;
	uint8_t *rdptr=0, *wtptr=0;
#ifndef USE_ROWS
	int prev[3]={0};
#endif
	int yuv[3]={0}, sym[3]={0}, estim[3]={0};
	RiceCoder ec={0};
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
	rdptr=rdbuf+BUFSIZE;
	wtptr=wtbuf;
	if(fwd)
	{
		for(int k=0;k<1024;++k)
		{
			int val=k-512;
			val=val<<23>>22;
			val=val^val>>31;
			packsign[k]=val;
		}
		fwrite(&tag, 1, 2, fdst);
		fwrite(&iw, 1, 3, fdst);
		fwrite(&ih, 1, 3, fdst);
		rice_init(&ec, wtptr, fdst);
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
#ifdef USE_ROWS
				int pred[3];

				estim[0]=rows[0][1+(0-1*NCH)*NROWS*NVAL];
				estim[1]=rows[0][1+(1-1*NCH)*NROWS*NVAL];
				estim[2]=rows[0][1+(2-1*NCH)*NROWS*NVAL];
				pred[0]=(rows[1][0+(0+0*NCH)*NROWS*NVAL]+rows[0][0+(0-1*NCH)*NROWS*NVAL])>>1;
				pred[1]=(rows[1][0+(1+0*NCH)*NROWS*NVAL]+rows[0][0+(1-1*NCH)*NROWS*NVAL])>>1;
				pred[2]=(rows[1][0+(2+0*NCH)*NROWS*NVAL]+rows[0][0+(2-1*NCH)*NROWS*NVAL])>>1;
#endif
				uint64_t x=acme_read(&rdptr, 3, fsrc);

				int nbypass0=(estim[0]>>RSHIFT)+1;
				int nbypass1=(estim[1]>>RSHIFT)+1;
				int nbypass2=(estim[2]>>RSHIFT)+1;

				
#ifdef _MSC_VER
				_BitScanReverse(&nbypass0, nbypass0);
				_BitScanReverse(&nbypass1, nbypass1);
				_BitScanReverse(&nbypass2, nbypass2);
#else
				nbypass0=__builtin_clz(nbypass0);
				nbypass1=__builtin_clz(nbypass1);
				nbypass2=__builtin_clz(nbypass2);
				nbypass0^=31;
				nbypass1^=31;
				nbypass2^=31;
#endif
				yuv[0]=x>> 0&255;
				yuv[1]=x>> 8&255;
				yuv[2]=x>>16&255;
				yuv[0]-=yuv[1];
				yuv[2]-=yuv[1];
				yuv[1]+=(yuv[0]+yuv[2])>>2;
				yuv[2]-=yuv[0]>>2;
#ifdef USE_ROWS
				rows[0][0+(0+0*NCH)*NROWS*NVAL]=yuv[0];
				rows[0][0+(1+0*NCH)*NROWS*NVAL]=yuv[1];
				rows[0][0+(2+0*NCH)*NROWS*NVAL]=yuv[2];
				sym[0]=packsignptr[yuv[0]-pred[0]];
				sym[1]=packsignptr[yuv[1]-pred[1]];
				sym[2]=packsignptr[yuv[2]-pred[2]];
#else
				sym[0]=packsignptr[yuv[0]-prev[0]];
				sym[1]=packsignptr[yuv[1]-prev[1]];
				sym[2]=packsignptr[yuv[2]-prev[2]];
				//sym[0]=sym[0]<<23>>23;
				//sym[1]=(int8_t)sym[1];
				//sym[2]=sym[2]<<23>>23;
				//sym[0]=sym[0]<<1^sym[0]>>31;
				//sym[1]=sym[1]<<1^sym[1]>>31;
				//sym[2]=sym[2]<<1^sym[2]>>31;
				prev[0]=yuv[0];
				prev[1]=yuv[1];
				prev[2]=yuv[2];
#endif
				
#ifdef FIFOVAL
				fifoval_enqueue(sym[2]<<18|sym[1]<<9|sym[0]);
#endif
				{
					uint64_t code[3];
					int nzeros[3], codelen[3];

					//new rice coder 2026
					nzeros[0]=sym[0]>>nbypass0;
					nzeros[1]=sym[1]>>nbypass1;
					nzeros[2]=sym[2]>>nbypass2;
					codelen[0]=nzeros[0]+1+nbypass0;
					codelen[1]=nzeros[1]+1+nbypass1;
					codelen[2]=nzeros[2]+1+nbypass2;
					code[0]=sym[0];
					code[1]=sym[1];
					code[2]=sym[2];
					code[0]<<=31-nbypass0;
					code[1]<<=31-nbypass1;
					code[2]<<=31-nbypass2;
					code[0]=(uint32_t)code[0];
					code[1]=(uint32_t)code[1];
					code[2]=(uint32_t)code[2];
					code[0]|=0x80000000;
					code[1]|=0x80000000;
					code[2]|=0x80000000;
					code[0]>>=31-nbypass0;
					code[1]>>=31-nbypass1;
					code[2]>>=31-nbypass2;
					if(nzeros[0]>RLIMIT-1)code[0]=sym[0], codelen[0]=32;
					if(nzeros[1]>RLIMIT-1)code[1]=sym[1], codelen[1]=32;
					if(nzeros[2]>RLIMIT-1)code[2]=sym[2], codelen[2]=32;
					for(int k=0;k<3;++k)
					{
						if(codelen[k]>=ec.nbits)//cache would overflow:  fill, flush, and repeat
						{
							codelen[k]-=(int)ec.nbits;
							ec.cache|=(uint64_t)code[k]>>codelen[k];
							code[k]&=0x7FFFFFFFULL>>(31-codelen[k]);
							acme_write(&ec.ptr, 8, ec.f, ec.cache);
							ec.cache=0;
							ec.nbits=64;
						}
						//now there is room for bypass:  0 <= nbypass < nbits <= 64
						if(codelen[k])
						{
							ec.nbits-=codelen[k];//emit remaining bypass to cache
							ec.cache|=code[k]<<ec.nbits;
						}
					}
				}
				//rice_enc(&ec, nbypass0, sym[0]);
				//rice_enc(&ec, nbypass1, sym[1]);
				//rice_enc(&ec, nbypass2, sym[2]);
#ifdef USE_ROWS
				rows[0][1+(0+0*NCH)*NROWS*NVAL]=(2*rows[0][1+(0-1*NCH)*NROWS*NVAL]+(sym[0]<<RSHIFT)+rows[1][1+(0+3*NCH)*NROWS*NVAL])>>2;
				rows[0][1+(1+0*NCH)*NROWS*NVAL]=(2*rows[0][1+(1-1*NCH)*NROWS*NVAL]+(sym[1]<<RSHIFT)+rows[1][1+(1+3*NCH)*NROWS*NVAL])>>2;
				rows[0][1+(2+0*NCH)*NROWS*NVAL]=(2*rows[0][1+(2-1*NCH)*NROWS*NVAL]+(sym[2]<<RSHIFT)+rows[1][1+(2+3*NCH)*NROWS*NVAL])>>2;
				rows[0]+=NCH*NROWS*NVAL;
				rows[1]+=NCH*NROWS*NVAL;
#else
				estim[0]+=((sym[0]<<RSHIFT)-estim[0])>>(RSHIFT+1);
				estim[1]+=((sym[1]<<RSHIFT)-estim[1])>>(RSHIFT+1);
				estim[2]+=((sym[2]<<RSHIFT)-estim[2])>>(RSHIFT+1);
#endif
			}
		}
		rice_flush(&ec);
		wtptr=ec.ptr;
	}
	else//dec
	{
		uint64_t cache2=0;
		int nbits=0;

		for(int k=0;k<512;++k)
		{
			int val=k;
			val=val>>1^-(val&1);
			packsign[k]=val;
		}
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		rice_init(&ec, rdptr, fsrc);
		fread(&ec.cache, 1, 8, fsrc);
		fread(&cache2, 1, 8, fsrc);
		
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
				int nbypass[3];
#ifdef USE_ROWS
				int pred[3];

				estim[0]=rows[0][1+(0-1*NCH)*NROWS*NVAL];
				estim[1]=rows[0][1+(1-1*NCH)*NROWS*NVAL];
				estim[2]=rows[0][1+(2-1*NCH)*NROWS*NVAL];
				pred[0]=(rows[1][0+(0+0*NCH)*NROWS*NVAL]+rows[0][0+(0-1*NCH)*NROWS*NVAL])>>1;
				pred[1]=(rows[1][0+(1+0*NCH)*NROWS*NVAL]+rows[0][0+(1-1*NCH)*NROWS*NVAL])>>1;
				pred[2]=(rows[1][0+(2+0*NCH)*NROWS*NVAL]+rows[0][0+(2-1*NCH)*NROWS*NVAL])>>1;
#endif

				nbypass[0]=(estim[0]>>RSHIFT)+1;
				nbypass[1]=(estim[1]>>RSHIFT)+1;
				nbypass[2]=(estim[2]>>RSHIFT)+1;
				
#ifdef _MSC_VER
				_BitScanReverse(&nbypass[0], nbypass[0]);
				_BitScanReverse(&nbypass[1], nbypass[1]);
				_BitScanReverse(&nbypass[2], nbypass[2]);
#else
				nbypass[0]=__builtin_clz(nbypass[0]);
				nbypass[1]=__builtin_clz(nbypass[1]);
				nbypass[2]=__builtin_clz(nbypass[2]);
				nbypass[0]^=31;
				nbypass[1]^=31;
				nbypass[2]^=31;
#endif
				for(int k=0;k<3;++k)
				{
					uint64_t code;
					int nzeros;
					int prefix;

					if(nbits>=64)
					{
						ec.cache=cache2;
						cache2=acme_read(&rdptr, 8, fsrc);
						nbits-=64;
					}
					code=ec.cache<<nbits;
					if(nbits)
						code|=cache2>>(64-nbits);
					nzeros=(int)_lzcnt_u64(code);
					prefix=nzeros<<nbypass[k];
					if(nzeros>RLIMIT-1)
						nzeros=RLIMIT-1, nbypass[k]=9, prefix=0;
					code<<=nzeros+1;
					code>>=64-nbypass[k];
					if(!nbypass[k])
						code=0;
					sym[k]=(int)code|prefix;
					nbits+=nzeros+1+nbypass[k];
				}
#ifdef FIFOVAL
				fifoval_check(sym[2]<<18|sym[1]<<9|sym[0]);
#endif
#ifdef USE_ROWS
				rows[0][1+(0+0*NCH)*NROWS*NVAL]=(2*rows[0][1+(0-1*NCH)*NROWS*NVAL]+(sym[0]<<RSHIFT)+rows[1][1+(0+3*NCH)*NROWS*NVAL])>>2;
				rows[0][1+(1+0*NCH)*NROWS*NVAL]=(2*rows[0][1+(1-1*NCH)*NROWS*NVAL]+(sym[1]<<RSHIFT)+rows[1][1+(1+3*NCH)*NROWS*NVAL])>>2;
				rows[0][1+(2+0*NCH)*NROWS*NVAL]=(2*rows[0][1+(2-1*NCH)*NROWS*NVAL]+(sym[2]<<RSHIFT)+rows[1][1+(2+3*NCH)*NROWS*NVAL])>>2;
#else
				estim[0]+=((sym[0]<<RSHIFT)-estim[0])>>(RSHIFT+1);
				estim[1]+=((sym[1]<<RSHIFT)-estim[1])>>(RSHIFT+1);
				estim[2]+=((sym[2]<<RSHIFT)-estim[2])>>(RSHIFT+1);
#endif
				sym[0]=packsign[sym[0]];
				sym[1]=packsign[sym[1]];
				sym[2]=packsign[sym[2]];
#ifdef USE_ROWS
				sym[0]+=pred[0];
				sym[1]+=pred[1];
				sym[2]+=pred[2];
#else
				sym[0]+=prev[0];
				sym[1]+=prev[1];
				sym[2]+=prev[2];
#endif
				sym[0]=sym[0]<<23>>23;
				sym[1]=sym[1]<<23>>23;
				sym[2]=sym[2]<<23>>23;
#ifdef USE_ROWS
				rows[0][0+(0+0*NCH)*NROWS*NVAL]=sym[0];
				rows[0][0+(1+0*NCH)*NROWS*NVAL]=sym[1];
				rows[0][0+(2+0*NCH)*NROWS*NVAL]=sym[2];
				rows[0]+=NCH*NROWS*NVAL;
				rows[1]+=NCH*NROWS*NVAL;
#else
				prev[0]=sym[0];
				prev[1]=sym[1];
				prev[2]=sym[2];
#endif
				sym[2]+=sym[0]>>2;
				sym[1]-=(sym[0]+sym[2])>>2;
				sym[2]+=sym[1];
				sym[0]+=sym[1];
				acme_write(&wtptr, 3, fdst, (uint64_t)sym[2]<<16|(uint64_t)sym[1]<<8|sym[0]);
			}
		}
	}
	if(wtptr>wtbuf)
		fwrite(wtbuf, 1, wtptr-wtbuf, fdst);
	fclose(fsrc);
	fclose(fdst);
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
	(void)&rice_enc;
	(void)&rice_dec;
	return 0;
}
