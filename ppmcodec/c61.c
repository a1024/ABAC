//gcc -Wall -O3 -mavx2 -mlzcnt -static -s -Dc61_codec=main c61.c -o c61
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
#ifdef _MSC_VER
#	include<intrin.h>
#elif defined __GNUC__
#	include<x86intrin.h>
#endif
#if defined _MSC_VER || defined _WIN32
#	define WIN32_LEAN_AND_MEAN
#	include<Windows.h>
#else
#	include<time.h>
#endif


#ifdef _MSC_VER
	#define LOUD
//	#define PRINT_RCT
//	#define PRINT_CODES
//	#define ESTIM_ENTROPY
//	#define GOLOMB_TEST

	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif


	#define SYNTH_PATH


enum
{
	BUFSIZE=512*1024,

	NCTX=7,
	DEPTH=8,
	RLIMIT=16,//DEPTH + RLIMIT <= 24
};

#if 1
#ifdef _MSC_VER
#	define ALIGN(N) __declspec(align(N))
#	define INLINE __forceinline static
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
static void guide_save(FILE *f, int iw, int ih)
{
	ptrdiff_t idx=0, size=0;
	
	size=(ptrdiff_t)3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(uint8_t*)malloc(size);
	if(!g_image)
	{
		CRASH("Alloc error");
		return;
	}
	idx=ftell(f);
	fread(g_image, 1, size, f);
	fseek(f, (long)idx, SEEK_SET);
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


//cRCT
#if 1
static const uint8_t perms[]=
{
	//1, 0, 2,
	//1, 2, 0,

	0, 1, 2,	2, 1, 0,
	2, 0, 1,	1, 0, 2,
	1, 2, 0,	0, 2, 1,
};
enum
{
	RCTBITS=2,
	RCTMAX=1<<RCTBITS,
	RCTLEVELS=RCTMAX+1,//+unity
	NRCTS=RCTLEVELS*RCTLEVELS*(RCTLEVELS+1)/2,
	NPERMS=_countof(perms)/3,
};
typedef struct _RCTInfo
{
	uint8_t pidx, uc0, vc0, vc1;
} RCTInfo;
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
typedef struct _AnalysisCtrs
{
	uint32_t yctr, uctrs[RCTLEVELS], vctrs[RCTLEVELS*(RCTLEVELS+1)/2];
} AnalysisCtrs;
static void crct_analysis(FILE *fsrc, int iw, int ih, RCTInfo *ret_rct)
{
	enum
	{
		STRIDE=5,
#if defined SYNTH_PATH
		STRIDE_S=5,
#endif
	};
	AnalysisCtrs counters[NPERMS]={0};
	uint8_t *image=0, *imptr=0;
	int32_t size=0;
	int rstr=0, kx=0, ky=0;
	long fidx=0;
	uint32_t bestscore=0;
	RCTInfo rct={0};
#if defined SYNTH_PATH
	uint32_t synthcounter=0;
	int issynth=0;
#endif
#ifdef PRINT_RCT
	int it=0;
#endif
	
	fidx=ftell(fsrc);
	size=(int32_t)3*STRIDE*iw;
	image=(uint8_t*)malloc(size);
	if(!image)
	{
		CRASH("Alloc error");
		return;
	}
	rstr=3*iw;
#if defined SYNTH_PATH
	for(ky=1;ky<=ih-STRIDE_S;ky+=STRIDE_S)
	{
		kx=1;
		imptr=image+rstr*(ky-ky/STRIDE*STRIDE)+3*kx;
		fread(image, 1, size, fsrc);
		for(;kx<=iw-STRIDE_S;kx+=STRIDE_S)
		{
			int rgbX[]=
			{
				imptr[0]-imptr[-3+0],
				imptr[1]-imptr[-3+1],
				imptr[2]-imptr[-3+2],
			//	abs((int8_t)(imptr[0]-(abs(imptr[-rstr+0]-imptr[-rstr-3+0])>abs(imptr[-3+0]-imptr[-rstr-3+0])?imptr[-rstr+0]:imptr[-3+0]))),
			//	abs((int8_t)(imptr[1]-(abs(imptr[-rstr+1]-imptr[-rstr-3+1])>abs(imptr[-3+1]-imptr[-rstr-3+1])?imptr[-rstr+1]:imptr[-3+1]))),
			//	abs((int8_t)(imptr[2]-(abs(imptr[-rstr+2]-imptr[-rstr-3+2])>abs(imptr[-3+2]-imptr[-rstr-3+2])?imptr[-rstr+2]:imptr[-3+2]))),
			};
			int rgbY[]=
			{
				imptr[0]-imptr[-rstr+0],
				imptr[1]-imptr[-rstr+1],
				imptr[2]-imptr[-rstr+2],
			//	abs((int8_t)(imptr[0]-imptr[-rstr+0]-imptr[-3+0]+imptr[-rstr-3+0])),
			//	abs((int8_t)(imptr[1]-imptr[-rstr+1]-imptr[-3+1]+imptr[-rstr-3+1])),
			//	abs((int8_t)(imptr[2]-imptr[-rstr+2]-imptr[-3+2]+imptr[-rstr-3+2])),
			};
			imptr+=3*STRIDE_S;
			synthcounter+=!rgbX[0]+!rgbY[0]+!rgbX[1]+!rgbY[1]+!rgbX[2]+!rgbY[2];
		//	synthcounters[0]+=rgb0[0]+rgb0[1]+rgb0[2];
		//	synthcounters[1]+=rgb1[0]+rgb1[1]+rgb1[2];
		}
	}
	issynth=(double)synthcounter>=(double)iw*ih*(3.*2*0.50/(STRIDE_S*STRIDE_S));// >= 50% flat

//	issynth=synthcounters[0]<=synthcounters[1];
	fseek(fsrc, fidx, SEEK_SET);
#endif
	for(ky=1;ky<=ih-STRIDE;ky+=STRIDE)
	{
		kx=1;
		imptr=image+rstr+3*kx;
		fread(image, 1, size, fsrc);
		for(;kx<=iw-STRIDE;kx+=STRIDE)
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
	free(image);
	fseek(fsrc, fidx, SEEK_SET);
	
#ifdef PRINT_RCT
	printf("Y:\n");
	for(int kp=0;kp<NPERMS;++kp)
		printf("  %12u", counters[kp].yctr);
	printf("\n");
	printf("U:\n");
	for(int k=0;k<RCTLEVELS;++k)
	{
		for(int kp=0;kp<NPERMS;++kp)
			printf("  %12u", counters[kp].uctrs[k]);
		printf("\n");
	}
	printf("V:\n");
	for(int k=0;k<RCTLEVELS*(RCTLEVELS+1)/2;++k)
	{
		for(int kp=0;kp<NPERMS;++kp)
			printf("  %12u", counters[kp].vctrs[k]);
		printf("\n");
	}
#endif
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
#ifdef PRINT_RCT
					++it;
#endif
					if(!bestscore||bestscore>score)
					{
						bestscore=score;
						rct.pidx=kp;
						rct.uc0=uc0;
						rct.vc0=vc0;
						rct.vc1=vc1;
#ifdef PRINT_RCT
						printf("%5d  ", it);
						print_rct(&rct);
						printf("  %10d  =  %10d + %10d + %10d\n"
							, score
							, currctrs->yctr
							, currctrs->uctrs[uc0]
							, currctrs->vctrs[idx]
						);
#endif
					}
				}
			}
		}
	}
#if defined SYNTH_PATH
//	if((double)bestscore*STRIDE*STRIDE<(double)iw*ih*(3<<RCTBITS))//synth test
	if(issynth)
		rct.pidx+=NPERMS;
#endif
#ifdef PRINT_RCT
	print_rct(&rct);
	printf("\n");
#endif
	*ret_rct=rct;
}
#endif


static uint8_t logtable[1<<DEPTH];
static uint32_t enctable[DEPTH<<DEPTH];
static uint8_t signpack[1<<DEPTH];
static uint8_t rdbuf[BUFSIZE+sizeof(uint32_t[2])], wtbuf[BUFSIZE+sizeof(uint32_t[2])];
#if defined SYNTH_PATH
static void largerice_enc(FILE *fdst, uint8_t *wtend, uint8_t **pwtptr, uint32_t *pcache, int *pnbits, uint32_t sym, int rk)
{
	uint8_t *wtptr;
	uint32_t cache;
	int nbits, inc;

	wtptr=*pwtptr;
	cache=*pcache;
	nbits=*pnbits;
	nbits+=sym>>rk;
	if(nbits>7)//nbits = number of written bits
	{
		do
		{
			*wtptr++=(uint8_t)cache;
			if(wtptr>=wtend)
			{
				fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
				wtptr-=BUFSIZE;
			}
			cache>>=8;
			nbits-=8;
		}while(nbits>7);
#if 0
		inc=nbits>>3&3;
		nbits+=32-nb0;
		do
		{
			*(uint32_t*)wtptr=cache;
			wtptr+=inc;
			if(wtptr>=wtend)
			{
				fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
				wtptr-=BUFSIZE;
				*(uint32_t*)(wtptr-inc)=cache;
			}
			cache>>=nbits&24;
			inc=sizeof(uint32_t);
			nbits-=32;
		}while(nbits>31);
#endif
	}
	cache|=((sym&((1<<rk)-1))<<1|1)<<nbits;		//1. cache |= code<<nbits
	nbits+=rk+1;					//2. nbits += len
	*(uint32_t*)wtptr=cache;			//3. store(cache)
	inc=nbits>>3;					//4. advance(nbits>>3)
	wtptr+=inc;
	if(wtptr>=wtend)
	{
		fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
		wtptr-=BUFSIZE;
		*(uint32_t*)(wtptr-inc)=cache;
	}
	cache>>=nbits&24;				//5. cache >>= nbits&24
	nbits&=7;					//6. nbits &= 7
	*pnbits=nbits;
	*pcache=cache;
	*pwtptr=wtptr;
}
static uint32_t largerice_dec(FILE *fsrc, uint8_t *rdend, uint8_t **prdptr, uint32_t *pcache, int *pnbits, int rk)
{
	uint8_t *rdptr;
	uint32_t cache, q, sym;
	int nbits, inc, nzeros;

	rdptr=*prdptr;
	cache=*pcache;
	nbits=*pnbits;
	q=0;
	do
	{
		cache|=*(uint32_t*)rdptr<<nbits;//nbits = number of read bits
		inc=nbits>>3^3;
		rdptr+=inc;
		if(rdptr>=rdend)
		{
			fread(rdbuf+sizeof(uint32_t), 1, BUFSIZE, fsrc);
			rdptr-=BUFSIZE;
			cache|=*(uint32_t*)(rdptr-inc)<<nbits;
		}
		nbits|=24;
		nzeros=(int)TZCNT32(cache);
		if(nzeros>24)
			nzeros=24;
		cache>>=nzeros;
		nbits-=nzeros;
#if defined _MSC_VER
		//if(nbits<0)//
		//	CRASH("");
#endif
		q+=nzeros;
	}while(nzeros>=24);
	sym=q<<rk;
	cache|=*(uint32_t*)rdptr<<nbits;		//1. cache |= read32()<<nbits
	inc=nbits>>3^3;					//2. advance(nbits>>3^3)
	rdptr+=inc;
	if(rdptr>=rdend)
	{
		fread(rdbuf+sizeof(uint32_t), 1, BUFSIZE, fsrc);
		rdptr-=BUFSIZE;
		cache|=*(uint32_t*)(rdptr-inc)<<nbits;
	}
	nbits|=24;					//3. nbits |= 24
	cache>>=1;//stopbit
	sym|=(int)(cache&((1<<rk)-1));			//4. sym = cache&mask
	cache>>=rk;					//5. cache >>= len
	nbits-=rk+1;					//6. nbits -= len
	*prdptr=rdptr;
	*pcache=cache;
	*pnbits=nbits;
	return sym;
}
#endif
int c61_codec(int argc, char **argv)
{
	enum
	{
		XPAD=8,
		NCH=4,
		NROWS=2,
		NVAL=4,

		PREDADD=4,
		TOTALADD=PREDADD+RCTBITS,
		RBITS=3,
	};
	typedef int32_t Pixel_t;
	const uint16_t tag='6'|'1'<<8;

	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	uint32_t c=0;
	int fwd=0, iw=0, ih=0;
	int64_t usize=0, csize=0;
	int psize=0;
	Pixel_t *pixels=0;
	uint32_t cache=0;
	int nbits=0, inc=0;
	uint8_t *rdptr=0, *wtptr=0, *rdend=0, *wtend=0;
	RCTInfo rct={0};
	int kx=0, ky=0, kc=0;
	int c0[3]={0}, c1[3]={0};
	int ysh=0, ush=0, vsh=0;
#if defined SYNTH_PATH
	int issynth=0;
	int rs_estim[4]={0};
	uint32_t *symbols=0;
#endif
	ALIGN(32) int yuv[4]={0}, nbypass[4]={0}, pred[4]={0}, offset[4]={0}, sym[4]={0};
#ifdef ESTIM_ENTROPY
	double esize=0;
#endif
#ifdef GOLOMB_TEST
	int gsizes[3]={0}, rsizes[3]={0};
#endif
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
	}
	else
	{
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
	psize=(iw+2*XPAD)*(int)sizeof(Pixel_t[NCH*NROWS*NVAL]);
	pixels=(Pixel_t*)_mm_malloc(psize, 32);
	if(!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	cache=0;
	nbits=0;
	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
	if(fwd)
	{
#ifdef ENABLE_GUIDE
		guide_save(fsrc, iw, ih);
#endif
		crct_analysis(fsrc, iw, ih, &rct);
#if defined LOUD && 0
		{
			double t2=time_sec()-t;
			print_rct(&rct);
			printf("   analysis  %12.6lf ms/MB\n", t2*1024*1024*1000/usize);
		}
#endif
		fwrite(&tag, 1, 2, fdst);
		fwrite(&iw, 1, 3, fdst);
		fwrite(&ih, 1, 3, fdst);
		fwrite(&rct, 1, sizeof(rct), fdst);
	}
	else
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	{
		int pidx=rct.pidx%NPERMS*3;

		ysh=perms[pidx+0]*8;
		ush=perms[pidx+1]*8;
		vsh=perms[pidx+2]*8;
#if defined SYNTH_PATH
		issynth=rct.pidx>=NPERMS;
#endif
	}
	for(int ks=0;ks<1<<DEPTH;++ks)
	{
		int val=31^LZCNT32(ks+1);
		if(val>NCTX-1)val=NCTX-1;
		logtable[ks]=val;
	}
	if(fwd)
	{
		//for(int ks=0;ks<1<<DEPTH;++ks)
		//	logtable[ks]<<=DEPTH;
		for(int ks=0;ks<1<<DEPTH;++ks)
		{
			int sym=(int8_t)(ks-(1<<DEPTH>>1));
			sym=sym<<1^sym>>31;
			for(int kb=0;kb<DEPTH;++kb)
			{
				uint32_t code;
				int nzeros, stopbit, nbypass, codelen;
				
				nbypass=kb;
				nzeros=sym>>kb;
				stopbit=nzeros<RLIMIT;
				codelen=nzeros+1+nbypass;
				if(nzeros>=RLIMIT)
					nzeros=RLIMIT, stopbit=0, nbypass=DEPTH, codelen=RLIMIT+DEPTH;
				code=(sym&((1<<nbypass)-1))<<(nzeros+stopbit)|stopbit<<nzeros;
				enctable[kb<<DEPTH|ks]=code<<8|codelen;
			}
		}
		for(int ks=0;ks<1<<DEPTH;++ks)//s->u
		{
			int sym=ks-(1<<DEPTH>>1);
			signpack[ks]=sym<<1^sym>>31;
		}
	}
	else
	{
		for(int ks=0;ks<1<<DEPTH;++ks)//u->s
			signpack[ks]=(uint8_t)(ks>>1^ks<<31>>31);
	}
	memset(pixels, 0, psize);
	rdptr=rdbuf+sizeof(uint32_t)+BUFSIZE;
	wtptr=wtbuf+sizeof(uint32_t);
	rdend=rdbuf+sizeof(uint32_t)+BUFSIZE;
	wtend=wtbuf+sizeof(uint32_t)+BUFSIZE;
	if(fwd)
		rdend-=3;
	else
		wtend-=3;
	c0[0]=0;
	c0[1]=rct.uc0<<PREDADD;
	c0[2]=rct.vc0<<PREDADD;
	c1[0]=0;
	c1[1]=0;
	c1[2]=rct.vc1<<PREDADD;
#if defined SYNTH_PATH
	if(issynth&&fwd)
		symbols=(uint32_t*)malloc(iw*sizeof(uint32_t[4]));
#endif
	for(ky=0;ky<ih;++ky)
	{
		const uint32_t masks[]={0, 1, 3, 7, 15, 31, 63, 127, 255};
		Pixel_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky+0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
		};
		uint32_t data, code[3]={0};

#if defined SYNTH_PATH &&1
		if(issynth)
		{
#define PREDICT_S(C)\
	do\
	{\
		N[C]=rows[1][0+(C+0*NCH)*NROWS*NVAL];\
		W[C]=rows[0][0+(C-1*NCH)*NROWS*NVAL];\
		pred[C]=\
		+rows[1][1+(C-1*NCH)*NROWS*NVAL]\
		+rows[1][1+(C+0*NCH)*NROWS*NVAL]\
		+rows[0][1+(C-1*NCH)*NROWS*NVAL]\
		>\
		+rows[1][2+(C-1*NCH)*NROWS*NVAL]\
		+rows[1][2+(C+0*NCH)*NROWS*NVAL]\
		+rows[0][2+(C-1*NCH)*NROWS*NVAL]\
		? N[C] : W[C];\
	}while(0)

			if(fwd)
			{
				uint32_t *symptr=symbols;
				int N[3]={0}, W[3]={0}, error[3]={0}, curr[3]={0}, run=0, rsum=0, ssum[3]={0};

				for(kx=0;kx<iw;++kx)
				{
					PREDICT_S(0);
					PREDICT_S(1);
					PREDICT_S(2);

					data=*(uint32_t*)rdptr;
					if(rdptr>=rdend)
					{
						fread(rdbuf+sizeof(uint32_t), 1, BUFSIZE, fsrc);
						rdptr-=BUFSIZE;
						data|=*(uint32_t*)rdptr;
					}
					rdptr+=3;

					yuv[0]=(uint8_t)(data>>ysh);
					yuv[1]=(uint8_t)(data>>ush);
					yuv[2]=(uint8_t)(data>>vsh);
					offset[1]=(rct.uc0*yuv[0])>>RCTBITS;
					offset[2]=(rct.vc0*yuv[0]+rct.vc1*yuv[1])>>RCTBITS;
					pred[1]+=offset[1];
					pred[2]+=offset[2];
					CLAMP(pred[1], 0, 255);
					CLAMP(pred[2], 0, 255);
					error[0]=(int8_t)(yuv[0]-pred[0]);
					error[1]=(int8_t)(yuv[1]-pred[1]);
					error[2]=(int8_t)(yuv[2]-pred[2]);

				//	if(ky==1&&kx==1240)//
				//		printf("");

					if(!(error[0]|error[1]|error[2]))
						++run;
					else
					{
						*symptr++=run;
						rsum+=run;
						run=0;
						sym[0]=(signpack+(1<<DEPTH>>1))[error[0]];
						sym[1]=(signpack+(1<<DEPTH>>1))[error[1]];
						sym[2]=(signpack+(1<<DEPTH>>1))[error[2]];
						*symptr++=sym[2]<<16|sym[1]<<8|sym[0];
						ssum[0]+=sym[0];
						ssum[1]+=sym[1];
						ssum[2]+=sym[2];
					}
					curr[0]=yuv[0];
					curr[1]=yuv[1]-offset[1];
					curr[2]=yuv[2]-offset[2];
					rows[0][0+(0+0*NCH)*NROWS*NVAL]=curr[0];
					rows[0][1+(0+0*NCH)*NROWS*NVAL]=abs(curr[0]-W[0]);
					rows[0][2+(0+0*NCH)*NROWS*NVAL]=abs(curr[0]-N[0]);
					rows[0][0+(1+0*NCH)*NROWS*NVAL]=curr[1];
					rows[0][1+(1+0*NCH)*NROWS*NVAL]=abs(curr[1]-W[1]);
					rows[0][2+(1+0*NCH)*NROWS*NVAL]=abs(curr[1]-N[1]);
					rows[0][0+(2+0*NCH)*NROWS*NVAL]=curr[2];
					rows[0][1+(2+0*NCH)*NROWS*NVAL]=abs(curr[2]-W[2]);
					rows[0][2+(2+0*NCH)*NROWS*NVAL]=abs(curr[2]-N[2]);
					rows[0]+=NVAL*NROWS*NCH;
					rows[1]+=NVAL*NROWS*NCH;
				}
				if(run)
				{
					*symptr++=run;
					rsum+=run;
					run=0;
				}
				{
					int ntotal, nsyms, nruns, xk[4]={0};

					ntotal=(int)(symptr-symbols);
					nsyms=ntotal>>1;
					nruns=ntotal-nsyms;
					xk[0]=nruns?31^LZCNT32(rsum/nruns+1):0;
					xk[1]=nsyms?31^LZCNT32(ssum[0]/nsyms+1):0;
					xk[2]=nsyms?31^LZCNT32(ssum[1]/nsyms+1):0;
					xk[3]=nsyms?31^LZCNT32(ssum[2]/nsyms+1):0;
					for(kc=0;kc<4;++kc)
					{
						uint32_t code;
						int s, nbypass, nzeros, stopbit, codelen;

						nbypass=31^LZCNT32((rs_estim[kc]>>RBITS)+1);

						s=xk[kc];
						nzeros=s>>nbypass;
						stopbit=nzeros<RLIMIT;
						codelen=nzeros+1+nbypass;
						if(nzeros>=RLIMIT)
							nzeros=RLIMIT, stopbit=0, nbypass=DEPTH, codelen=RLIMIT+DEPTH;
						code=(s&((1<<nbypass)-1))<<(nzeros+stopbit)|stopbit<<nzeros;
						cache|=code<<nbits;
						nbits+=codelen;
						*(uint32_t*)wtptr=cache;
						inc=nbits>>3;
						wtptr+=inc;
						if(wtptr>=wtend)
						{
							fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
							wtptr-=BUFSIZE;
							*(uint32_t*)(wtptr-inc)=cache;
						}
						cache>>=nbits&24;
						nbits&=7;

						rs_estim[kc]+=((s<<RBITS)-rs_estim[kc])>>(RBITS+1);
					}
					symptr=symbols;
					for(kx=0;kx<ntotal;++kx)
					{
						int s[4], s2, nzeros[4], stopbit, codelen, nbypass;

						s[0]=*symptr++;
						s2=*symptr++;
						s[1]=(uint8_t)(s2>>0);
						s[2]=(uint8_t)(s2>>8);
						s[3]=(uint8_t)(s2>>16);

					//	nzeros[0]=s[0]>>xk[0];
						nzeros[1]=s[1]>>xk[1];
						nzeros[2]=s[2]>>xk[2];
						nzeros[3]=s[3]>>xk[3];

					//	if(ky==1&&s[0]==1003)//
					//		printf("");//

						largerice_enc(fdst, wtend, &wtptr, &cache, &nbits, s[0], xk[0]);
						++kx;
						if(kx>=ntotal)
							break;
						for(kc=1;kc<4;++kc)
						{
							nbypass=xk[kc];
							s2=s[kc];
							stopbit=nzeros[kc]<RLIMIT;
							codelen=nzeros[kc]+1+nbypass;
							if(nzeros[kc]>=RLIMIT)
								nzeros[kc]=RLIMIT, stopbit=0, nbypass=DEPTH, codelen=RLIMIT+DEPTH;
							cache|=((s2&((1<<nbypass)-1))<<(nzeros[kc]+stopbit)|stopbit<<nzeros[kc])<<nbits;
							nbits+=codelen;
							*(uint32_t*)wtptr=cache;
							inc=nbits>>3;
							wtptr+=inc;
							if(wtptr>=wtend)
							{
								fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
								wtptr-=BUFSIZE;
								*(uint32_t*)(wtptr-inc)=cache;
							}
							cache>>=nbits&24;
							nbits&=7;
						}
#if 0
						for(kc=0;kc<4;++kc)
						{
							s2=s[kc];
							nzeros=s2>>xk[kc];
							stopbit=nzeros<RLIMIT;
							codelen=nzeros+1+nbypass;
							if(nzeros>=RLIMIT)
								nzeros=RLIMIT, stopbit=0, nbypass=DEPTH, codelen=RLIMIT+DEPTH;
							code=(s2&((1<<nbypass)-1))<<(nzeros+stopbit)|stopbit<<nzeros;
							cache|=code<<nbits;
							nbits+=codelen;
							*(uint32_t*)wtptr=cache;
							inc=nbits>>3;
							wtptr+=inc;
							if(wtptr>=wtend)
							{
								fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
								wtptr-=BUFSIZE;
								*(uint32_t*)(wtptr-inc)=cache;
							}
							cache>>=nbits&24;
							nbits&=7;
						}
#endif
					}
				}
			}
			else//dec
			{
				uint32_t run;
#ifdef _MSC_VER
				uint32_t r0;
#endif
				int N[3]={0}, W[3]={0}, curr[3]={0}, xk[4]={0}, nbypass, nzeros;

				for(kc=0;kc<4;++kc)
				{
					nbypass=31^LZCNT32((rs_estim[kc]>>RBITS)+1);
					
					cache|=*(uint32_t*)rdptr<<nbits;
					inc=nbits>>3^3;
					rdptr+=inc;
					if(rdptr>=rdend)
					{
						fread(rdbuf+sizeof(uint32_t), 1, BUFSIZE, fsrc);
						rdptr-=BUFSIZE;
						cache|=*(uint32_t*)(rdptr-inc)<<nbits;
					}
					nbits|=24;
					nzeros=(int)TZCNT32(cache);
					xk[kc]=nzeros<<nbypass;
					if(nzeros>RLIMIT-1)
						nzeros=RLIMIT-1, nbypass=DEPTH, xk[kc]=0;
					cache>>=nzeros+1;
					xk[kc]|=(int)(cache&masks[nbypass]);
					cache>>=nbypass;
					nbits-=nzeros+1+nbypass;
					
					rs_estim[kc]+=((xk[kc]<<RBITS)-rs_estim[kc])>>(RBITS+1);
				}
				for(kx=0;kx<iw;++kx)
				{

				//	if(ky==1&&kx==1240)//
				//	if(ky==1&&kx==248)//
				//		printf("");
#ifdef _MSC_VER
					r0=
#endif
					run=largerice_dec(fsrc, rdend, &rdptr, &cache, &nbits, xk[0]);
					while(run>0)
					{
						PREDICT_S(0);
						PREDICT_S(1);
						PREDICT_S(2);
						yuv[0]=pred[0];
						offset[1]=(rct.uc0*yuv[0])>>RCTBITS;
						pred[1]+=offset[1];
						CLAMP(pred[1], 0, 255);
						yuv[1]=pred[1];
						offset[2]=(rct.vc0*yuv[0]+rct.vc1*yuv[1])>>RCTBITS;
						pred[2]+=offset[2];
						CLAMP(pred[2], 0, 255);
						yuv[2]=pred[2];

						data=yuv[2]<<vsh|yuv[1]<<ush|yuv[0]<<ysh;
#ifdef ENABLE_GUIDE
						if(memcmp(&data, g_image+3*(iw*ky+kx), 3))
						{
							printf("Y %d  X %d\n", ky, kx);
							CRASH("guide");
							return 1;
						}
#endif
						*(uint32_t*)wtptr=data;
						if(wtptr>=wtend)
						{
							fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
							wtptr-=BUFSIZE;
							*(uint32_t*)wtptr=data;
						}
						wtptr+=3;

						++kx;
						--run;
						curr[0]=yuv[0];
						curr[1]=yuv[1]-offset[1];
						curr[2]=yuv[2]-offset[2];
						rows[0][0+(0+0*NCH)*NROWS*NVAL]=curr[0];
						rows[0][1+(0+0*NCH)*NROWS*NVAL]=abs(curr[0]-W[0]);
						rows[0][2+(0+0*NCH)*NROWS*NVAL]=abs(curr[0]-N[0]);
						rows[0][0+(1+0*NCH)*NROWS*NVAL]=curr[1];
						rows[0][1+(1+0*NCH)*NROWS*NVAL]=abs(curr[1]-W[1]);
						rows[0][2+(1+0*NCH)*NROWS*NVAL]=abs(curr[1]-N[1]);
						rows[0][0+(2+0*NCH)*NROWS*NVAL]=curr[2];
						rows[0][1+(2+0*NCH)*NROWS*NVAL]=abs(curr[2]-W[2]);
						rows[0][2+(2+0*NCH)*NROWS*NVAL]=abs(curr[2]-N[2]);
						rows[0]+=NVAL*NROWS*NCH;
						rows[1]+=NVAL*NROWS*NCH;
					}
					if(kx>=iw)
						break;

					for(kc=0;kc<3;++kc)
					{
						nbypass=xk[kc+1];
						cache|=*(uint32_t*)rdptr<<nbits;
						inc=nbits>>3^3;
						rdptr+=inc;
						if(rdptr>=rdend)
						{
							fread(rdbuf+sizeof(uint32_t), 1, BUFSIZE, fsrc);
							rdptr-=BUFSIZE;
							cache|=*(uint32_t*)(rdptr-inc)<<nbits;
						}
						nbits|=24;
						nzeros=(int)TZCNT32(cache);
						sym[kc]=nzeros<<nbypass;
						if(nzeros>RLIMIT-1)
							nzeros=RLIMIT-1, nbypass=DEPTH, sym[kc]=0;
						cache>>=nzeros+1;
						sym[kc]|=(int)(cache&masks[nbypass]);
						cache>>=nbypass;
						nbits-=nzeros+1+nbypass;
					}
					sym[0]=(int8_t)signpack[sym[0]];
					sym[1]=(int8_t)signpack[sym[1]];
					sym[2]=(int8_t)signpack[sym[2]];

					PREDICT_S(0);
					PREDICT_S(1);
					PREDICT_S(2);
					yuv[0]=(uint8_t)(sym[0]+pred[0]);
					offset[1]=(rct.uc0*yuv[0])>>RCTBITS;
					pred[1]+=offset[1];
					CLAMP(pred[1], 0, 255);
					yuv[1]=(uint8_t)(sym[1]+pred[1]);
					offset[2]=(rct.vc0*yuv[0]+rct.vc1*yuv[1])>>RCTBITS;
					pred[2]+=offset[2];
					CLAMP(pred[2], 0, 255);
					yuv[2]=(uint8_t)(sym[2]+pred[2]);
					
					data=yuv[2]<<vsh|yuv[1]<<ush|yuv[0]<<ysh;
#ifdef ENABLE_GUIDE
					if(memcmp(&data, g_image+3*(iw*ky+kx), 3))
					{
						printf("Y %d  X %d\n", ky, kx);
						CRASH("guide");
						return 1;
					}
#endif
					*(uint32_t*)wtptr=data;
					if(wtptr>=wtend)
					{
						fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
						wtptr-=BUFSIZE;
						*(uint32_t*)wtptr=data;
					}
					wtptr+=3;

					curr[0]=yuv[0];
					curr[1]=yuv[1]-offset[1];
					curr[2]=yuv[2]-offset[2];
					rows[0][0+(0+0*NCH)*NROWS*NVAL]=curr[0];
					rows[0][1+(0+0*NCH)*NROWS*NVAL]=abs(curr[0]-W[0]);
					rows[0][2+(0+0*NCH)*NROWS*NVAL]=abs(curr[0]-N[0]);
					rows[0][0+(1+0*NCH)*NROWS*NVAL]=curr[1];
					rows[0][1+(1+0*NCH)*NROWS*NVAL]=abs(curr[1]-W[1]);
					rows[0][2+(1+0*NCH)*NROWS*NVAL]=abs(curr[1]-N[1]);
					rows[0][0+(2+0*NCH)*NROWS*NVAL]=curr[2];
					rows[0][1+(2+0*NCH)*NROWS*NVAL]=abs(curr[2]-W[2]);
					rows[0][2+(2+0*NCH)*NROWS*NVAL]=abs(curr[2]-N[2]);
					rows[0]+=NVAL*NROWS*NCH;
					rows[1]+=NVAL*NROWS*NCH;
				}
			}
			continue;
		}
#endif
		//				-2
		//			-5	10	2
		//	1	-3	13	[?]			284.457069   267.943686 MB/s
		//
		//	+8*(+N+W)
		//	+5*(+aW)
		//	+2*(+aN+NE-WW)
		//	+1*(-WW+WWW)
#if 1
#define PREDICT(CH) pred[CH]=\
	+8*(+rows[1][+CH+(+0*NVAL+0)*NROWS*NCH]+rows[0][+CH+(-1*NVAL+0)*NROWS*NCH])\
	+5*(+rows[0][+CH+(-1*NVAL+2)*NROWS*NCH])\
	+2*(+rows[1][+CH+(+0*NVAL+2)*NROWS*NCH]+rows[1][+CH+(+1*NVAL+0)*NROWS*NCH]-rows[0][+CH+(-2*NVAL+0)*NROWS*NCH])\
	+1*(-rows[0][+CH+(-2*NVAL+0)*NROWS*NCH]+rows[0][+CH+(-3*NVAL+0)*NROWS*NCH])\

#define UPDATE_AUX(CH) rows[0][+CH+(+0*NVAL+2)*NROWS*NCH]=\
	rows[0][+CH+(+0*NVAL+0)*NROWS*NCH]-rows[1][+CH+(+0*NVAL+0)*NROWS*NCH]\

#endif

		//				-2
		//			-5	10	2
		//	1	-3	13	[?]
#if 0
#define PREDICT(CH) pred[CH]=\
	+10*rows[1][+CH+(+0*NVAL+0)*NROWS*NCH]\
	+13*rows[0][+CH+(-1*NVAL+0)*NROWS*NCH]\
	- 5*rows[1][+CH+(-1*NVAL+0)*NROWS*NCH]\
	- 3*rows[0][+CH+(-2*NVAL+0)*NROWS*NCH]\
	+ 2*rows[1][+CH+(+1*NVAL+0)*NROWS*NCH]\
	- 2*rows[0][+CH+(+0*NVAL+0)*NROWS*NCH]\
	+ 1*rows[0][+CH+(-3*NVAL+0)*NROWS*NCH]\

#define UPDATE_AUX(...)
#endif


		//			8	12	8
		//	20	[17]	?
		//
		//	+16*(eW+sym)
		//	+ 8*(eNE+eNEE+eNEEE)
		//	+ 4*(eW+eNEE)
		//	+ 1*(sym)
		//	+ 9
#if 1
#define UPDATE(CH)\
	do\
	{\
		int eW=rows[0][+CH+(-1*NVAL+1)*NROWS*NCH], eNEE=rows[1][+CH+(+2*NVAL+1)*NROWS*NCH];\
		\
		nbypass[CH]=rows[0][+CH+(+0*NCH+1)*NROWS*NVAL]=(\
			+16*(eW+sym[CH])\
			+ 8*(rows[1][+CH+(+1*NCH+1)*NROWS*NVAL]+eNEE+rows[1][+CH+(+3*NCH+1)*NROWS*NVAL])\
			+ 4*(eW+eNEE)\
			+ 1*sym[CH]\
		+9)>>6;\
	}while(0)

#endif

		//			8	12	8
		//	20	[17]	?
#if 0
#define UPDATE(CH) nbypass[CH]=\
	rows[0][+CH+(+0*NVAL+1)*NROWS*NCH]=(\
		+20*rows[0][+CH+(-1*NVAL+1)*NROWS*NCH]\
		+17*sym[CH]\
		+ 8*rows[1][+CH+(+1*NVAL+1)*NROWS*NCH]\
		+12*rows[1][+CH+(+2*NVAL+1)*NROWS*NCH]\
		+ 8*rows[1][+CH+(+3*NVAL+1)*NROWS*NCH]\
	+9)>>6

#endif
		nbypass[0]=4;
		nbypass[1]=4;
		nbypass[2]=4;
		PREDICT(0);
		PREDICT(1);
		PREDICT(2);
		if(fwd)
		{
			for(kx=0;kx<iw;++kx)
			{
				data=*(uint32_t*)rdptr;
				if(rdptr>=rdend)
				{
					fread(rdbuf+sizeof(uint32_t), 1, BUFSIZE, fsrc);
					rdptr-=BUFSIZE;
					data|=*(uint32_t*)rdptr;
				}
				rdptr+=3;
				yuv[0]=(uint8_t)(data>>ysh);
				yuv[1]=(uint8_t)(data>>ush);
				yuv[2]=(uint8_t)(data>>vsh);
				offset[0]=0;
				offset[1]=c0[1]*yuv[0];
				offset[2]=c0[2]*yuv[0]+c1[2]*yuv[1];
				if(c0[1]<1<<TOTALADD)offset[1]-=offset[1]>>6;
				if(c0[2]+c1[2]<1<<TOTALADD)offset[2]-=offset[2]>>6;
				pred[0]=(pred[0]+offset[0]+((1<<TOTALADD>>1)+2))>>TOTALADD;
				pred[1]=(pred[1]+offset[1]+((1<<TOTALADD>>1)+2))>>TOTALADD;
				pred[2]=(pred[2]+offset[2]+((1<<TOTALADD>>1)+2))>>TOTALADD;
				CLAMP(pred[0], 0, 255);
				CLAMP(pred[1], 0, 255);
				CLAMP(pred[2], 0, 255);
				sym[0]=(int8_t)(yuv[0]-pred[0]);
				sym[1]=(int8_t)(yuv[1]-pred[1]);
				sym[2]=(int8_t)(yuv[2]-pred[2]);

				code[0]=(enctable+(1<<DEPTH>>1))[(nbypass[0]<<DEPTH)+sym[0]];
				code[1]=(enctable+(1<<DEPTH>>1))[(nbypass[1]<<DEPTH)+sym[1]];
				code[2]=(enctable+(1<<DEPTH>>1))[(nbypass[2]<<DEPTH)+sym[2]];
				sym[0]=(signpack+(1<<DEPTH>>1))[sym[0]];
				sym[1]=(signpack+(1<<DEPTH>>1))[sym[1]];
				sym[2]=(signpack+(1<<DEPTH>>1))[sym[2]];
#ifdef GOLOMB_TEST
				for(kc=0;kc<3;++kc)
				{
					int x, m, k, u, b;

					x=sym[kc];

					m=rows[0][+kc+(-1*NVAL+1)*NROWS*NCH]+1;
					k=31^LZCNT32(m);
					u=x/m;
					b=x%m;
					gsizes[0]+=u>=RLIMIT?0:u+1;
					gsizes[1]+=u>=RLIMIT?0:k+(b>=(1<<(k+1)));
					gsizes[2]+=u>=RLIMIT?RLIMIT+DEPTH:0;
				//	gsize+=u>=RLIMIT?RLIMIT+DEPTH:u+1+k+(x>=(1<<(k+1)));
					
					k=31^LZCNT32(rows[0][+kc+(-1*NVAL+1)*NROWS*NCH]+1);
					u=x>>k;
					rsizes[0]+=u>=RLIMIT?0:u+1;
					rsizes[1]+=u>=RLIMIT?0:k;
					rsizes[1]+=u>=RLIMIT?RLIMIT+DEPTH:0;
				//	rsize+=u>=RLIMIT?RLIMIT+DEPTH:u+1+k;
				}
#endif
				for(kc=0;kc<3;++kc)
				{
#ifdef ESTIM_ENTROPY
					{
						static int totalbits=0, n0=0;
						uint32_t c=code[kc]>>8;
						int len=(uint8_t)code[kc];
						for(int kb=0;kb<len;++kb)
						{
							int bit=c>>kb&1;
							if(!kb)
							{
								esize+=log2(totalbits+2)-log2((bit?totalbits-n0:n0)+1);
								++totalbits;
								n0+=bit^1;
								if(totalbits>511)
								{
									totalbits>>=1;
									n0>>=1;
								}
							}
							else
							{
								++esize;
							}
						}
					}
#endif
#if defined PRINT_CODES &&0
					{
						static int prev=0;
						uint32_t c=code[kc]>>8;
						int len=(uint8_t)code[kc];
						printf(" ");
						for(int kb=0;kb<len;++kb)
						{
							int bit=c>>kb&1;
							//if(prev^bit)
							//	printf(" ");
							printf("%c", '0'+bit);
							prev=bit;
						}
					}
#endif
					cache|=code[kc]>>8<<nbits;
					nbits+=(uint8_t)code[kc];
					*(uint32_t*)wtptr=cache;
					inc=nbits>>3;
					wtptr+=inc;
					if(wtptr>=wtend)
					{
						fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
						wtptr-=BUFSIZE;
						*(uint32_t*)(wtptr-inc)=cache;
					}
					cache>>=nbits&24;
					nbits&=7;
				}
				rows[0][+0+(+0*NVAL+0)*NROWS*NCH]=(yuv[0]<<RCTBITS)-((offset[0]-7)>>PREDADD);
				rows[0][+1+(+0*NVAL+0)*NROWS*NCH]=(yuv[1]<<RCTBITS)-((offset[1]-7)>>PREDADD);
				rows[0][+2+(+0*NVAL+0)*NROWS*NCH]=(yuv[2]<<RCTBITS)-((offset[2]-7)>>PREDADD);
				UPDATE(0);
				UPDATE(1);
				UPDATE(2);
				UPDATE_AUX(0);
				UPDATE_AUX(1);
				UPDATE_AUX(2);
				rows[0]+=NVAL*NROWS*NCH;
				rows[1]+=NVAL*NROWS*NCH;
				PREDICT(0);
				PREDICT(1);
				PREDICT(2);
				nbypass[0]=logtable[nbypass[0]];//next nW
				nbypass[1]=logtable[nbypass[1]];
				nbypass[2]=logtable[nbypass[2]];
			}
		}
		else//dec
		{
			int nzeros;

			for(kx=0;kx<iw;++kx)
			{
				for(kc=0;kc<3;++kc)
				{
					cache|=*(uint32_t*)rdptr<<nbits;
					inc=nbits>>3^3;
					rdptr+=inc;
					if(rdptr>=rdend)
					{
						fread(rdbuf+sizeof(uint32_t), 1, BUFSIZE, fsrc);
						rdptr-=BUFSIZE;
						cache|=*(uint32_t*)(rdptr-inc)<<nbits;
					}
					nbits|=24;
					nzeros=(int)TZCNT32(cache);
					sym[kc]=nzeros<<nbypass[kc];
					if(nzeros>RLIMIT-1)
						nzeros=RLIMIT-1, nbypass[kc]=DEPTH, sym[kc]=0;
					cache>>=nzeros+1;
					sym[kc]|=(int)(cache&masks[nbypass[kc]]);
					cache>>=nbypass[kc];
					nbits-=nzeros+1+nbypass[kc];
				}
				UPDATE(0);
				UPDATE(1);
				UPDATE(2);
				sym[0]=(int8_t)signpack[sym[0]];
				sym[1]=(int8_t)signpack[sym[1]];
				sym[2]=(int8_t)signpack[sym[2]];

				offset[0]=0;
				pred[0]=(pred[0]+offset[0]+((1<<TOTALADD>>1)+2))>>TOTALADD;
				CLAMP(pred[0], 0, 255);
				yuv[0]=(uint8_t)(sym[0]+pred[0]);
#ifdef ENABLE_GUIDE
				if(yuv[0]!=g_image[3*(iw*ky+kx)+perms[rct.pidx*3+0]])
				{
					printf("Y %d  X %d  C0\n", ky, kx);
					CRASH("guide");
					return 1;
				}
#endif
				
				offset[1]=c0[1]*yuv[0];
				if(c0[1]<1<<TOTALADD)offset[1]-=offset[1]>>6;
				pred[1]=(pred[1]+offset[1]+((1<<TOTALADD>>1)+2))>>TOTALADD;
				CLAMP(pred[1], 0, 255);
				yuv[1]=(uint8_t)(sym[1]+pred[1]);
#ifdef ENABLE_GUIDE
				if(yuv[1]!=g_image[3*(iw*ky+kx)+perms[rct.pidx*3+1]])
				{
					printf("Y %d  X %d  C1\n", ky, kx);
					CRASH("guide");
					return 1;
				}
#endif
				
				offset[2]=c0[2]*yuv[0]+c1[2]*yuv[1];
				if(c0[2]+c1[2]<1<<TOTALADD)offset[2]-=offset[2]>>6;
				pred[2]=(pred[2]+offset[2]+((1<<TOTALADD>>1)+2))>>TOTALADD;
				CLAMP(pred[2], 0, 255);
				yuv[2]=(uint8_t)(sym[2]+pred[2]);
#ifdef ENABLE_GUIDE
				if(yuv[2]!=g_image[3*(iw*ky+kx)+perms[rct.pidx*3+2]])
				{
					printf("Y %d  X %d  C2\n", ky, kx);
					CRASH("guide");
					return 1;
				}
#endif
				data=yuv[2]<<vsh|yuv[1]<<ush|yuv[0]<<ysh;
				*(uint32_t*)wtptr=data;
				if(wtptr>=wtend)
				{
					fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
					wtptr-=BUFSIZE;
					*(uint32_t*)wtptr=data;
				}
				wtptr+=3;
				rows[0][+0+(+0*NVAL+0)*NROWS*NCH]=(yuv[0]<<RCTBITS)-((offset[0]-7)>>PREDADD);
				rows[0][+1+(+0*NVAL+0)*NROWS*NCH]=(yuv[1]<<RCTBITS)-((offset[1]-7)>>PREDADD);
				rows[0][+2+(+0*NVAL+0)*NROWS*NCH]=(yuv[2]<<RCTBITS)-((offset[2]-7)>>PREDADD);
				UPDATE_AUX(0);
				UPDATE_AUX(1);
				UPDATE_AUX(2);
				rows[0]+=NCH*NROWS*NVAL;
				rows[1]+=NCH*NROWS*NVAL;
				PREDICT(0);
				PREDICT(1);
				PREDICT(2);
				nbypass[0]=logtable[nbypass[0]];//next nW
				nbypass[1]=logtable[nbypass[1]];
				nbypass[2]=logtable[nbypass[2]];
			}
		}
	}
	if(fwd)
	{
		*(uint32_t*)wtptr=cache;
		if(wtptr>=wtend)
		{
			fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
			wtptr-=BUFSIZE;
			*(uint32_t*)wtptr=cache;
		}
		++wtptr;
#if defined SYNTH_PATH
		if(issynth)
			free(symbols);
#endif
	}
	if(wtptr>wtbuf+sizeof(uint32_t))
		fwrite(wtbuf+sizeof(uint32_t), 1, wtptr-(wtbuf+sizeof(uint32_t)), fdst);
	_mm_free(pixels);
	fclose(fsrc);
	fclose(fdst);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		struct stat info={0};
		stat(dstfn, &info);
		csize=info.st_size;
#ifdef GOLOMB_TEST
		printf("%12.2f gsize.t\n", (gsizes[0]+gsizes[1]+gsizes[2])/8.);
		printf("%12.2f gsize.u\n", gsizes[0]/8.);
		printf("%12.2f gsize.b\n", gsizes[1]/8.);
		printf("%12.2f gsize.o\n", gsizes[2]/8.);
		printf("%12.2f rsize.t\n", (rsizes[0]+rsizes[1]+rsizes[2])/8.);
		printf("%12.2f rsize.u\n", rsizes[0]/8.);
		printf("%12.2f rsize.b\n", rsizes[1]/8.);
		printf("%12.2f rsize.o\n", rsizes[2]/8.);
	//	printf("%12.2f gsize\n", gsize/8.);
	//	printf("%12.2f rsize\n", rsize/8.);
#endif
#ifdef ESTIM_ENTROPY
		printf("%12.2lf\n", esize/8.);
#endif
		printf("CWH=3*%d*%d  \"%s\"\n", iw, ih, srcfn);
		printf("%10td->%10td  %8.4f%%  %12.6f:1  BPD %12.6f\n"
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
	(void)usize;
	(void)csize;
	(void)&time_sec;
	(void)&print_rct;
	return 0;
}
