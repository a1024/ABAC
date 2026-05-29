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

//	#define RICE_EXPERIMENT
	#define MEASURE_PSNR

//	#define FIFOVAL 2
#endif

	#define USE_RLE
	#define NO_RCT
	#define NEARLOSSLESS
	#define USE_ROWS
//	#define USE_SELECT	//incompatible with near
//	#define USE_CG		//


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
#ifdef NEARLOSSLESS
	NEARSHIFT=4,
#endif
#ifdef USE_RLE
	MAXRUNBITS=30,
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
		CRASH("");
		//return 1;
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
static uint32_t enctable[1<<(DEPTH+3)];
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
#ifdef USE_RLE
	int64_t res=0, csize=0;
	uint8_t *cbuf=0, *cptr[3]={0};
	int run[3]={0}, runestim[3]={0};
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
#ifdef USE_RLE
	cbuf=(uint8_t*)malloc(usize);
	if(!cbuf)
	{
		CRASH("Alloc error");
		return 1;
	}
	if(fwd)
	{
		res=(int64_t)iw*ih;
		cptr[0]=cbuf+0*res;
		cptr[1]=cbuf+1*res;
		cptr[2]=cbuf+2*res;
	}
	else
	{
		struct stat info={0};
		int64_t csizes[3]={0}, nread;
		int idx;

		fread(csizes+0, 1, 4, fsrc);
		fread(csizes+1, 1, 4, fsrc);
		fread(csizes+2, 1, 4, fsrc);
		idx=(int)ftell(fsrc);
		stat(srcfn, &info);
		nread=fread(cbuf, 1, (size_t)info.st_size-idx, fsrc);
		cptr[0]=cbuf;
		cptr[1]=cptr[0]+csizes[0];
		cptr[2]=cptr[1]+csizes[1];
		if(cptr[2]+csizes[2]!=cbuf+nread)
			printf("size expected %lld  got %lld+%lld+%lld = %lld\n"
				, nread
				, csizes[0]
				, csizes[1]
				, csizes[2]
				, csizes[0]+csizes[1]+csizes[2]
			);
	}
#endif
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
	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
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
	//	static uint8_t *const wtend=wtbuf+sizeof(uint64_t)+BUFSIZE-sizeof(uint64_t);
#ifdef USE_RLE
		uint64_t cache[3]={0};
		int nbits[3]={0};
#else
		uint64_t cache=0;
		int nbits=0;
#endif
		for(int k=0;k<1024;++k)
		{
			int val=k-512;
			val=val<<(32-DEPTH)>>(32-DEPTH-1);
			val=val^val>>31;
			packsign[k]=val;
		}
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
		fwrite(&tag, 1, 2, fdst);
		fwrite(&iw, 1, 3, fdst);
		fwrite(&ih, 1, 3, fdst);
		for(int ky=0, idx=0;ky<ih;++ky)
		{
#ifdef USE_ROWS
			int16_t *rows[]=
			{
				pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
				pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			};
#endif
			for(int kx=0;kx<iw;++kx, ++idx)
			{
				uint64_t reg;
#ifndef USE_RLE
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
#elif defined USE_SELECT
				{
					int NW[]=
					{
						rows[1][0+(0-1*NCH)*NROWS*NVAL],
						rows[1][0+(1-1*NCH)*NROWS*NVAL],
						rows[1][0+(2-1*NCH)*NROWS*NVAL],
					};
					int N[]=
					{
						rows[1][0+(0+0*NCH)*NROWS*NVAL],
						rows[1][0+(1+0*NCH)*NROWS*NVAL],
						rows[1][0+(2+0*NCH)*NROWS*NVAL],
					};
					int W[]=
					{
						rows[0][0+(0-1*NCH)*NROWS*NVAL],
						rows[0][0+(1-1*NCH)*NROWS*NVAL],
						rows[0][0+(2-1*NCH)*NROWS*NVAL],
					};
					int t0[]=
					{
						N[0]-NW[0],
						N[1]-NW[1],
						N[2]-NW[2],
					};
					int t1[]=
					{
						W[0]-NW[0],
						W[1]-NW[1],
						W[2]-NW[2],
					};
					pred[0]=(t0[0]-t1[0])*(t0[0]+t1[0])>0?N[0]:W[0];
					pred[1]=(t0[1]-t1[1])*(t0[1]+t1[1])>0?N[1]:W[1];
					pred[2]=(t0[2]-t1[2])*(t0[2]+t1[2])>0?N[2]:W[2];
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
				sym[0]>>=NEARSHIFT;
				sym[1]>>=NEARSHIFT;
				sym[2]>>=NEARSHIFT;
#ifdef MEASURE_PSNR
				int yuv0[3];
				yuv0[0]=yuv[0];
				yuv0[1]=yuv[1];
				yuv0[2]=yuv[2];
#endif
				yuv[0]=sym[0]<<NEARSHIFT;
				yuv[1]=sym[1]<<NEARSHIFT;
				yuv[2]=sym[2]<<NEARSHIFT;
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
#if defined FIFOVAL && !defined USE_RLE
				fifoval_enqueue(sym[2]<<DEPTH*2^sym[1]<<DEPTH^sym[0]);
				fifoval_enqueue(yuv[2]<<DEPTH*2^yuv[1]<<DEPTH^yuv[0]);
#endif
#ifdef USE_RLE
				for(int kc=0;kc<3;++kc)
				{
//#ifdef FIFOVAL
//					if((uint32_t)(idx-1600)<200&&kc==FIFOVAL)//
//						printf("enc 0x%016zX  %2d\n", (size_t)cache[kc], nbits[kc]);
//#endif
					//if(ky==2&&kx==162&&kc==2)//
					//if(kc==2)
					//{
					//	if(ky==2&&(kx==160||kx==161||kx==162))//
					//		printf("");
					//}
					if(!sym[kc]&&idx<res-1)
						++run[kc];
					else
					{
						//uint64_t code=(uint64_t)(run[kc]!=0)<<63;
						//int nb2=nbits[kc]+1;
						//cache[kc]|=code>>nbits[kc];
						//if(nb2>=64)
						//{
						//	*(uint64_t*)cptr[kc]=cache[kc];
						//	cptr[kc]+=sizeof(uint64_t);
						//	cache[kc]=code<<(64-nbits[kc]);
						//	nb2-=64;
						//	if(!nbits[kc])
						//		cache[kc]=0;
						//}
						//nbits[kc]=nb2;
						
						//if(kc==2)
						//{
						//	if(ky==2&&kx>=160)//
						//		printf("");
						//}
						if(run[kc])
						{
							int nbypass2=31-_lzcnt_u32((runestim[kc]>>1)+1);
							int nzeros=run[kc]>>nbypass2, stopbit=1, bypassmask=(1<<nbypass2)-1, codelen=nzeros+stopbit+nbypass2;
							uint64_t code=run[kc];
							if(nzeros>MAXRUNBITS-1)
								nzeros=MAXRUNBITS, stopbit=0, nbypass2=MAXRUNBITS, bypassmask=~0, codelen=2*MAXRUNBITS;
							code&=(uint32_t)bypassmask;
							code|=(uint64_t)stopbit<<nbypass2;
							code<<=63-codelen;
							++codelen;//prepend true bit (run)
							code|=1ULL<<63;
#ifdef FIFOVAL
							if(kc==FIFOVAL)
								fifoval_enqueue(1<<MAXRUNBITS|run[kc]);
#endif
							
							codelen+=nbits[kc];
							cache[kc]|=code>>nbits[kc];
							if(codelen>=64)
							{
								*(uint64_t*)cptr[kc]=cache[kc];
								cptr[kc]+=sizeof(uint64_t);
								cache[kc]=code<<(64-nbits[kc]);
								codelen-=64;
								if(!nbits[kc])
									cache[kc]=0;
							}
							nbits[kc]=codelen;

							runestim[kc]+=(2*run[kc]-runestim[kc])>>2;
						}
						{
							uint64_t code=enctable[8*sym[kc]+nbypass[kc]];
							int codelen=(uint8_t)code;

							code>>=8;
							if(!run[kc])
								++codelen;//prepend zero bit (not run)
							run[kc]=0;
							code<<=64-codelen;
#ifdef FIFOVAL
							if(kc==FIFOVAL)
								fifoval_enqueue(sym[kc]);
#endif

							codelen+=nbits[kc];
							cache[kc]|=code>>nbits[kc];
							if(codelen>=64)
							{
								*(uint64_t*)cptr[kc]=cache[kc];
								cptr[kc]+=sizeof(uint64_t);
								cache[kc]=code<<(64-nbits[kc]);
								codelen-=64;
								if(!nbits[kc])
									cache[kc]=0;
							}
							nbits[kc]=codelen;
						}
					}
				}
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
#ifdef USE_RLE
		*(uint64_t*)cptr[0]=cache[0]; cptr[0]+=sizeof(uint64_t);
		*(uint64_t*)cptr[1]=cache[1]; cptr[1]+=sizeof(uint64_t);
		*(uint64_t*)cptr[2]=cache[2]; cptr[2]+=sizeof(uint64_t);
		cptr[0]=(uint8_t*)(cptr[0]-(cbuf+0*res));
		cptr[1]=(uint8_t*)(cptr[1]-(cbuf+1*res));
		cptr[2]=(uint8_t*)(cptr[2]-(cbuf+2*res));
		csize=0;
		csize+=fwrite(cptr+0, 1, 4, fdst);
		csize+=fwrite(cptr+1, 1, 4, fdst);
		csize+=fwrite(cptr+2, 1, 4, fdst);
		csize+=fwrite(cbuf+0*res, 1, (size_t)cptr[0], fdst);
		csize+=fwrite(cbuf+1*res, 1, (size_t)cptr[1], fdst);
		csize+=fwrite(cbuf+2*res, 1, (size_t)cptr[2], fdst);
#else
		memcpy(wtptr, &cache, sizeof(uint64_t));
		if(wtptr>=wtend)
		{
			fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
			wtptr-=BUFSIZE;
			memcpy(wtptr, &cache, sizeof(uint64_t));
		}
		wtptr+=sizeof(uint64_t);
		//acme_write(&wtptr, 8, fdst, cache);
#endif
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
	//	static uint8_t *const rdend=rdbuf+sizeof(uint64_t)+BUFSIZE-sizeof(uint64_t);
		static uint8_t *const wtend=wtbuf+sizeof(uint64_t)+BUFSIZE-3;
#ifdef USE_RLE
		uint64_t cache[3]={0}, cache2[3]={0};
		int nbits[3]={0};
		int run0[3]={0};
#else
		uint64_t cache=0, cache2=0;
		int nbits=0;
#endif
		for(int k=0;k<512;++k)
		{
			int val=k;
			val=val>>1^-(val&1);
			packsign[k]=val;
		}
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
#ifdef USE_RLE
		cache[0]=*(uint64_t*)cptr[0], cptr[0]+=sizeof(uint64_t), cache2[0]=*(uint64_t*)cptr[0], cptr[0]+=sizeof(uint64_t);
		cache[1]=*(uint64_t*)cptr[1], cptr[1]+=sizeof(uint64_t), cache2[1]=*(uint64_t*)cptr[1], cptr[1]+=sizeof(uint64_t);
		cache[2]=*(uint64_t*)cptr[2], cptr[2]+=sizeof(uint64_t), cache2[2]=*(uint64_t*)cptr[2], cptr[2]+=sizeof(uint64_t);
#else
		fread(&cache, 1, 8, fsrc);
		fread(&cache2, 1, 8, fsrc);
	//	cache	=acme_read(&rdptr, 8, fsrc);
	//	cache2	=acme_read(&rdptr, 8, fsrc);
#endif
#ifdef _MSC_VER
		int idx=0;
#endif
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
#elif defined USE_SELECT
				{
					int NW[]=
					{
						rows[1][0+(0-1*NCH)*NROWS*NVAL],
						rows[1][0+(1-1*NCH)*NROWS*NVAL],
						rows[1][0+(2-1*NCH)*NROWS*NVAL],
					};
					int N[]=
					{
						rows[1][0+(0+0*NCH)*NROWS*NVAL],
						rows[1][0+(1+0*NCH)*NROWS*NVAL],
						rows[1][0+(2+0*NCH)*NROWS*NVAL],
					};
					int W[]=
					{
						rows[0][0+(0-1*NCH)*NROWS*NVAL],
						rows[0][0+(1-1*NCH)*NROWS*NVAL],
						rows[0][0+(2-1*NCH)*NROWS*NVAL],
					};
					int t0[]=
					{
						N[0]-NW[0],
						N[1]-NW[1],
						N[2]-NW[2],
					};
					int t1[]=
					{
						W[0]-NW[0],
						W[1]-NW[1],
						W[2]-NW[2],
					};
					pred[0]=(t0[0]-t1[0])*(t0[0]+t1[0])>0?N[0]:W[0];
					pred[1]=(t0[1]-t1[1])*(t0[1]+t1[1])>0?N[1]:W[1];
					pred[2]=(t0[2]-t1[2])*(t0[2]+t1[2])>0?N[2]:W[2];
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
#ifdef USE_RLE
				for(int kc=0;kc<3;++kc)
				{
//#ifdef FIFOVAL
//					if((uint32_t)(idx-1600)<200&&kc==FIFOVAL)//
//						printf("dec 0x%016zX  0x%016zX  %2d\n", (size_t)cache[kc], (uint64_t)cache2[kc], nbits[kc]);
//#endif
//					if(ky==2&&(kx==160||kx==161||kx==162)&&kc==2)//
//						printf("");
					if(run[kc])
						--run[kc], sym[kc]=0;
					else
					{
						if(nbits[kc]>=64)
						{
							cache[kc]=cache2[kc];
							cache2[kc]=*(uint64_t*)cptr[kc];
							cptr[kc]+=sizeof(uint64_t);
							nbits[kc]-=64;
						}
						code=cache[kc]<<nbits[kc];
						if(nbits[kc])
							code|=cache2[kc]>>(64-nbits[kc]);

						if(!run0[kc])//previous symbol was not from a run
						{
							run[kc]=code>>63;
							code<<=1;
							++nbits[kc];
						}
						if(run[kc])
						{
							int nbypass2=31-_lzcnt_u32((runestim[kc]>>1)+1);

							nzeros=(int)_lzcnt_u64(code);
							prefix=nzeros<<nbypass2;
							if(nzeros>MAXRUNBITS-1)
								nzeros=MAXRUNBITS-1, nbypass2=MAXRUNBITS, prefix=0;
							code<<=nzeros+1;
							run[kc]=(int)(code>>(64-nbypass2));
							if(!nbypass2)
								run[kc]=0;
							//code<<=nbypass2;
							run[kc]|=prefix;
							run0[kc]=run[kc];
#ifdef FIFOVAL
							if(kc==FIFOVAL)
								fifoval_check(1<<MAXRUNBITS|run[kc]);
#endif
							nbits[kc]+=nzeros+1+nbypass2;
							
							runestim[kc]+=(2*run[kc]-runestim[kc])>>2;
							--run[kc];
							sym[kc]=0;
						}
						else
						{
							nzeros=(int)_lzcnt_u64(code);
							prefix=nzeros<<nbypass[kc];
							if(nzeros>RLIMIT-1)
								nzeros=RLIMIT-1, nbypass[kc]=DEPTH, prefix=0;
							code<<=nzeros+1;
							sym[kc]=(int)(code>>(64-nbypass[kc]));
							if(!nbypass[kc])
								sym[kc]=0;
							//code<<=nbypass[kc];
							sym[kc]|=prefix;
#ifdef FIFOVAL
							if(kc==FIFOVAL)
								fifoval_check(sym[kc]);
#endif
							nbits[kc]+=nzeros+1+nbypass[kc];
							run0[kc]=0;
						}
					}
				}
#else
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
#endif
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
				sym[0]<<=NEARSHIFT;
				sym[1]<<=NEARSHIFT;
				sym[2]<<=NEARSHIFT;
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
#if defined FIFOVAL && !defined USE_RLE
				if(fifoval_check(sym[2]<<DEPTH*2^sym[1]<<DEPTH^sym[0]))
				{
					printf("%016zX\n", (size_t)cache);
					printf("%016zX\n", (size_t)cache2);
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
#ifdef _MSC_VER
				++idx;
#endif
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
	(void)csize;
	(void)&time_sec2;
	(void)packsignptr;
	return 0;
}
