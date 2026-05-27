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

//	#define FIFOVAL
#endif

//	#define USE_RGB
	#define USE_LOG2TABLE	// 1 KB, packsign: 2 KB
	#define USE_ENCTABLE	//16 KB
//	#define USE_ENCTABLE2	//32 KB, slow
	#define USE_ROWS
//	#define USE_CG


enum
{
#ifdef USE_RGB
	DEPTH=8,
#else
	DEPTH=9,
#endif
//	RLIMIT=23,
	RLIMIT=12,
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
#ifndef MEDIAN3V_CLOB
//clobbers A B C
#define MEDIAN3V_CLOB(M, A, B, C)\
	M=A, A=A>B?B:A, B=M>B?M:B,\
	M=B, B=B>C?C:B, C=M>C?M:C,\
	M=A, A=A>B?B:A, M=M>B?M:B
#endif
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
/*
overflow:
|                    ______left______   ______right_____
|                   /                \ /                \
|buf1start ... ... [datastart  buf1end|buf2start  dataend] ...
|                   \________________    _______________/
|                                    size
*/
static uint8_t rdbuf[BUFSIZE+sizeof(uint64_t[4])]={0}, wtbuf[BUFSIZE+sizeof(uint64_t[4])]={0};
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

#ifdef USE_LOG2TABLE
static uint8_t log2table[1<<(DEPTH+RSHIFT)];
#else
static const int nbypasstable[]={0, 1, 2, 3, 4, 5, 6, 7, 7, 7};
#endif
#ifdef USE_ENCTABLE
#ifdef USE_ENCTABLE2
static uint32_t enctable2[1<<(1+DEPTH+3)], *const encptr2=enctable2+(1<<(1+DEPTH+3)>>1);
#endif
#if !defined USE_ENCTABLE2
static uint32_t enctable[1<<(DEPTH+3)];
#endif
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
#ifdef USE_ENCTABLE2
	uint32_t sym[3]={0};
#else
	int sym[3]={0};
#endif
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
#ifdef USE_LOG2TABLE
	for(int ks=0;ks<1<<(DEPTH+RSHIFT);++ks)
	{
		int val=31-_lzcnt_u32((ks>>RSHIFT)+1);
		if(val>7)
			val=7;
		log2table[ks]=val;
	}
#endif
	if(fwd)
	{
		uint64_t cache=0;
		int nbits=0;

		for(int k=0;k<1024;++k)
		{
			int val=k-512;
			val=val<<(32-DEPTH)>>(32-DEPTH-1);
			val=val^val>>31;
			packsign[k]=val;
		}
#ifdef USE_ENCTABLE
#ifdef USE_ENCTABLE2
		for(int k=0;k<_countof(enctable2);++k)
		{
			uint32_t code;
			int sym, nbypass, nzeros, codelen;

			sym=(k-_countof(enctable2)/2)>>3;
			nbypass=k&7;
			sym=sym<<(32-DEPTH)>>(32-DEPTH-1);
			sym=sym^sym>>31;
			nzeros=sym>>nbypass;
			codelen=nzeros+1+nbypass;
			code=sym;
			nbypass^=31;
			code<<=nbypass;
			code|=1<<31;
			code>>=nbypass;
			if(nzeros>RLIMIT-1)
				code=sym, codelen=RLIMIT+DEPTH;
			//if(code<<(32-DEPTH-5)>>(32-DEPTH-5)!=code)
			//	CRASH("");
			enctable2[k]=((32*code+codelen)<<DEPTH)+sym;
		}
#else
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
				uint64_t reg, code[3];
				int nbypass[3], codelen[3];
#ifndef USE_ENCTABLE
				int nzeros[3];
#endif

				reg=acme_read(&rdptr, 3, fsrc);
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
#ifdef USE_LOG2TABLE
				nbypass[0]=log2table[estim[0]];
				nbypass[1]=log2table[estim[1]];
				nbypass[2]=log2table[estim[2]];
#else
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
				nbypass[0]=nbypasstable[nbypass[0]];
				nbypass[1]=nbypasstable[nbypass[1]];
				nbypass[2]=nbypasstable[nbypass[2]];
#endif
				yuv[0]=(uint8_t)(reg>> 0);
				yuv[1]=(uint8_t)(reg>> 8);
				yuv[2]=(uint8_t)(reg>>16);
#ifndef USE_RGB
				yuv[0]-=yuv[1];
				yuv[2]-=yuv[1];
				yuv[1]+=(yuv[0]+yuv[2])>>2;
				yuv[2]-=yuv[0]>>2;
#endif
#ifdef USE_ENCTABLE2
				code[0]=encptr2[8*(yuv[0]-pred[0])+nbypass[0]];
				code[1]=encptr2[8*(yuv[1]-pred[1])+nbypass[1]];
				code[2]=encptr2[8*(yuv[2]-pred[2])+nbypass[2]];
				sym[0]=(uint32_t)code[0]<<(32-DEPTH);
				sym[1]=(uint32_t)code[1]<<(32-DEPTH);
				sym[2]=(uint32_t)code[2]<<(32-DEPTH);
				codelen[0]=(uint32_t)code[0]>>DEPTH;
				codelen[1]=(uint32_t)code[1]>>DEPTH;
				codelen[2]=(uint32_t)code[2]>>DEPTH;
				code[0]>>=DEPTH+5;
				code[1]>>=DEPTH+5;
				code[2]>>=DEPTH+5;
				sym[0]>>=32-DEPTH;
				sym[1]>>=32-DEPTH;
				sym[2]>>=32-DEPTH;
#ifdef FIFOVAL
				fifoval_enqueue(sym[2]<<DEPTH*2^sym[1]<<DEPTH^sym[0]);
				fifoval_enqueue(yuv[2]<<DEPTH*2^yuv[1]<<DEPTH^yuv[0]);
#endif
				codelen[0]&=31;
				codelen[1]&=31;
				codelen[2]&=31;
				estim[0]+=(int)((sym[0]<<RSHIFT)-estim[0])>>(RSHIFT+1);
				estim[1]+=(int)((sym[1]<<RSHIFT)-estim[1])>>(RSHIFT+1);
				estim[2]+=(int)((sym[2]<<RSHIFT)-estim[2])>>(RSHIFT+1);
#if 0
				if(	sym[0]!=packsignptr[yuv[0]-pred[0]]
				||	sym[1]!=packsignptr[yuv[1]-pred[1]]
				||	sym[2]!=packsignptr[yuv[2]-pred[2]]
				)
					CRASH("");
				if(	256*code[0]+codelen[0]!=enctable[8*sym[0]+nbypass[0]]
				||	256*code[1]+codelen[1]!=enctable[8*sym[1]+nbypass[1]]
				||	256*code[2]+codelen[2]!=enctable[8*sym[2]+nbypass[2]]
				)
					CRASH("");
#endif
				pred[0]=yuv[0];
				pred[1]=yuv[1];
				pred[2]=yuv[2];
#else
				sym[0]=packsignptr[yuv[0]-pred[0]];
				sym[1]=packsignptr[yuv[1]-pred[1]];
				sym[2]=packsignptr[yuv[2]-pred[2]];
				//sym[0]=sym[0]<<23>>23;
				//sym[1]=(int8_t)sym[1];
				//sym[2]=sym[2]<<23>>23;
				//sym[0]=sym[0]<<1^sym[0]>>31;
				//sym[1]=sym[1]<<1^sym[1]>>31;
				//sym[2]=sym[2]<<1^sym[2]>>31;
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
#ifdef USE_ENCTABLE
				code[0]=enctable[8*sym[0]+nbypass[0]];
				code[1]=enctable[8*sym[1]+nbypass[1]];
				code[2]=enctable[8*sym[2]+nbypass[2]];
				codelen[0]=(uint8_t)code[0];
				codelen[1]=(uint8_t)code[1];
				codelen[2]=(uint8_t)code[2];
				code[0]>>=8;
				code[1]>>=8;
				code[2]>>=8;
#else
				nzeros[0]=sym[0]>>nbypass[0];
				nzeros[1]=sym[1]>>nbypass[1];
				nzeros[2]=sym[2]>>nbypass[2];
				codelen[0]=nzeros[0]+1+nbypass[0];
				codelen[1]=nzeros[1]+1+nbypass[1];
				codelen[2]=nzeros[2]+1+nbypass[2];
				code[0]=sym[0];
				code[1]=sym[1];
				code[2]=sym[2];
				nbypass[0]^=63;
				nbypass[1]^=63;
				nbypass[2]^=63;
				code[0]<<=nbypass[0];
				code[1]<<=nbypass[1];
				code[2]<<=nbypass[2];
				code[0]|=1ULL<<63;
				code[1]|=1ULL<<63;
				code[2]|=1ULL<<63;
				code[0]>>=nbypass[0];
				code[1]>>=nbypass[1];
				code[2]>>=nbypass[2];
				if(nzeros[0]>RLIMIT-1)code[0]=sym[0], codelen[0]=RLIMIT+DEPTH;
				if(nzeros[1]>RLIMIT-1)code[1]=sym[1], codelen[1]=RLIMIT+DEPTH;
				if(nzeros[2]>RLIMIT-1)code[2]=sym[2], codelen[2]=RLIMIT+DEPTH;
#endif
#endif
#if 1
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
#else
				codelen[1]+=codelen[0];
				codelen[2]+=codelen[1];
				code[0]<<=64-codelen[0];
				code[1]<<=64-codelen[1];
				code[2]<<=64-codelen[2];
				//code[0]=_rotr64(code[0], codelen[0]);//no inline
				//code[1]=_rotr64(code[1], codelen[1]);
				//code[2]=_rotr64(code[2], codelen[2]);
				code[0]|=code[1];
				code[0]|=code[2];
#endif
				codelen[2]+=nbits;
				cache|=code[0]>>nbits;
				if(codelen[2]>=64)
				{
					acme_write(&wtptr, 8, fdst, cache);
					cache=code[0]<<(64-nbits);
					codelen[2]-=64;
					if(!nbits)
						cache=0;
				}
				nbits=codelen[2];
			}
		}
		acme_write(&wtptr, 8, fdst, cache);
	}
	else//dec
	{
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
				uint64_t code;
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
				
#ifdef USE_LOG2TABLE
				nbypass[0]=log2table[estim[0]];
				nbypass[1]=log2table[estim[1]];
				nbypass[2]=log2table[estim[2]];
#else
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
				nbypass[0]=nbypasstable[nbypass[0]];
				nbypass[1]=nbypasstable[nbypass[1]];
				nbypass[2]=nbypasstable[nbypass[2]];
#endif
				if(nbits>=64)
				{
					cache=cache2;
					cache2=acme_read(&rdptr, 8, fsrc);
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
				nbits+=nzeros+1+nbypass[2];
#ifdef FIFOVAL
				if(fifoval_check(sym[2]<<DEPTH*2^sym[1]<<DEPTH^sym[0]))
				{
					printf("%016llX\n", cache);
					printf("%016llX\n", cache2);
					CRASH("");
				}
#endif
#ifdef USE_ROWS
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
				sym[0]+=pred[0];
				sym[1]+=pred[1];
				sym[2]+=pred[2];
#ifdef USE_RGB
				sym[0]=(uint8_t)sym[0];
				sym[1]=(uint8_t)sym[1];
				sym[2]=(uint8_t)sym[2];
#else
				sym[0]<<=32-DEPTH;
				sym[1]<<=32-DEPTH;
				sym[2]<<=32-DEPTH;
				sym[0]>>=32-DEPTH;
				sym[1]>>=32-DEPTH;
				sym[2]>>=32-DEPTH;
#endif
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
#ifndef USE_RGB
				sym[2]+=sym[0]>>2;
				sym[1]-=(sym[0]+sym[2])>>2;
				sym[2]+=sym[1];
				sym[0]+=sym[1];
#endif
				acme_write(&wtptr, 3, fdst, (uint64_t)sym[2]<<16|(uint64_t)sym[1]<<8|sym[0]);
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
#ifdef USE_ENCTABLE2
	(void)encptr2;
#endif
	return 0;
}
