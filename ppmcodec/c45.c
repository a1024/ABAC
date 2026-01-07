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
//	#define FIFOVAL
#endif


#if 0
#define PREDLIST\
	PRED(N+W-NW)\
	PRED(N)\
	PRED(W)\
	PRED(NE)\

#endif
enum
{
#if 0
	SHIFT=18,
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED
#endif
	
	SBUFSIZE=512*1024,
	PBUFSIZE=SBUFSIZE/3*3,

	GRBITS=3,
	
	NCTX=18,

	XPAD=8,
	NCH=3,
	NROWS=4,
	NVAL=2,

	PROBBITS=12,
};


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
static uint8_t *g_buf=0;
static uint8_t* guide_save(FILE *f, ptrdiff_t size)
{
	ptrdiff_t idx=0;
	uint8_t *buf2=0;

	buf2=(uint8_t*)malloc(size);
	if(!buf2)
	{
		CRASH("Alloc error");
		return 0;
	}
	idx=ftell(f);
	fread(buf2, 1, size, f);
	fseek(f, (long)idx, SEEK_SET);
	return buf2;
}
#endif
#endif


#ifdef FIFOVAL
ptrdiff_t fifoidx=0, fifocap=0, fifoidx2=0;
uint32_t *fifoval=0;
static void valfifo_enqueue(uint32_t val)
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
static void valfifo_check(uint32_t val)
{
	uint32_t val0=fifoval[fifoidx2++];
	if(val!=val0)
	{
		--fifoidx2;
		printf(
			"\n"
			"FIFO Error\n"
			"    0x%08X  !=  original 0x%08X\n"
			"\n"
			, val, val0
		);
		for(int k=-32;k<32;++k)
		{
			ptrdiff_t idx=fifoidx2+k;
			if((size_t)idx<fifoidx)
			{
				printf(
					"%8td 0x%08X"
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


#ifdef LOUD
static volatile double tread=0, twrite=0;
#endif

uint8_t imbuf[PBUFSIZE+3], streambuf[SBUFSIZE];
uint32_t hweight[3*NCTX];
uint32_t hists[3*NCTX*256];

uint16_t enccdf[3*NCTX*256];//enc

uint32_t CDF2sym[3*NCTX<<PROBBITS];//dec
static void normalize_enc(int kc, int rescale)
{
	uint64_t invsum=0;
	uint32_t cdf=0, sum2=0, bypass=0;
	uint32_t *currh=hists+256*kc;
	uint16_t *currs=enccdf+256*kc;
	int ks=0;

	invsum=hweight[kc];
	if(invsum<256)
	{
		invsum+=256;
		bypass=1;
	}
	invsum=(0x100000000ULL+invsum-1)/invsum*((1ULL<<PROBBITS)-256);

	for(ks=0, cdf=0, sum2=0;ks<=255;++ks)
	{
		uint32_t freq=currh[ks];
		currs[ks]=(uint32_t)(cdf*invsum>>32)+ks;
		cdf+=freq+bypass;
		if(rescale)
		{
			freq>>=1;
			currh[ks]=freq;
			sum2+=freq;
		}
	}
	if(rescale)
		hweight[kc]=sum2;
}
static void normalize_dec(int kc, int rescale)
{
	uint64_t invsum=0;
	uint32_t cdf=0, sum2=0, bypass=0;
	uint32_t *currh=hists+256*kc, *currCDF2sym=CDF2sym+((ptrdiff_t)kc<<PROBBITS);
	uint16_t CDF[257];
	int ks=0;
	
	invsum=hweight[kc];
	if(invsum<256)
	{
		invsum+=256;
		bypass=1;
	}
	invsum=(0x100000000ULL+invsum-1)/invsum*((1ULL<<PROBBITS)-256);

	for(ks=0, cdf=0, sum2=0;ks<=255;++ks)
	{
		uint32_t freq=currh[ks];
		CDF[ks]=(uint32_t)(cdf*invsum>>32)+ks;
#ifdef _MSC_VER
		if(CDF[ks]>(1<<PROBBITS))
			CRASH("");
#endif
		cdf+=freq+bypass;
		if(rescale)
		{
			freq>>=1;
			currh[ks]=freq;
			sum2+=freq;
		}
	}
	CDF[256]=1<<PROBBITS;
	if(rescale)
		hweight[kc]=sum2;
	for(ks=0, cdf=0;ks<256;++ks)
	{
		uint32_t val=0;
		int freq=CDF[ks], end=CDF[ks+1];
		for(val=(end-freq)<<(PROBBITS+8)|freq<<8|ks;freq<end;)//for AC
			currCDF2sym[freq++]=val;
	}
}
int c45_codec(int argc, char **argv)
{
	const uint16_t tag='4'|'5'<<8;

	const char *srcfn=0, *dstfn=0;
	int64_t c=0;
	FILE *fsrc=0, *fdst=0;
	int fwd=0, iw=0, ih=0;
	ptrdiff_t usize=0, csize=0;
	uint8_t *imptr=0, *imend=0, *streamptr=0;
	int psize=0;
	int16_t *pixels=0;
	int kx=0, ky=0, kc=0, idx=0;
#ifdef LOUD
	double t=0;
#endif
	uint64_t nread=0;
	uint64_t low=0, range=0xFFFFFFFFFFFF, code=0;

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
		csize=2;
		iw=0;
		ih=0;
		csize+=fread(&iw, 1, 3, fsrc);
		csize+=fread(&ih, 1, 3, fsrc);
	}
	if(iw<1||ih<1)
	{
		CRASH("Unsupported source file");
		return 1;
	}
	usize=(ptrdiff_t)3*iw*ih;
	psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
	memset(pixels, 0, psize);
	memset(hists, 0, sizeof(hists));
	memset(hweight, 0, sizeof(hweight));
	if(fwd)//write header
	{
#ifdef ENABLE_GUIDE
		g_buf=guide_save(fsrc, usize);
#endif
		for(kc=0;kc<3*NCTX;++kc)
			normalize_enc(kc, 0);

		imend=imbuf+PBUFSIZE;
		if(imend>imbuf+usize)
			imend=imbuf+usize;
		imptr=imend;

		csize+=fwrite(&tag, 1, 2, fdst);
		csize+=fwrite(&iw, 1, 3, fdst);
		csize+=fwrite(&ih, 1, 3, fdst);

		streamptr=streambuf;
	}
	else
	{
		for(kc=0;kc<3*NCTX;++kc)
			normalize_dec(kc, 0);
		imend=imbuf+PBUFSIZE;
		if(imend>imbuf+usize)
			imend=imbuf+usize;
		imptr=imbuf;
		{
			ptrdiff_t nact=0, chunksize=0;
#ifdef LOUD
			double t2=0;
#endif
			
#ifdef LOUD
			t2=time_sec();
#endif
			nact=fread(&chunksize, 1, 3, fsrc);
			csize+=nact;
			if(nact!=3||(size_t)chunksize>SBUFSIZE)
			{
				CRASH("Corrupt archive \"%s\"", srcfn);
				return 1;
			}
			nact=fread(streambuf, 1, chunksize, fsrc);
			csize+=nact;
			if(nact!=chunksize)
			{
				CRASH("Truncated archive \"%s\"", srcfn);
				return 1;
			}
			streamptr=streambuf;
			
#ifdef LOUD
			tread+=time_sec()-t2;
#endif
			code=0;
			code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);
			code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);
		}

		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	}
	for(ky=0, idx=0;ky<ih;++ky)
	{
#if 0
		int estim[NPREDS]={0};
#endif
		uint8_t rgb[3]={0}, yuv[3]={0};
		int sym=0, error=0;
		int freq=0, cdf=0;
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(kx=0;kx<iw;++kx, idx+=3)
		{
			int offset=0;
			if(fwd)
			{
				uint8_t *imptr0=imptr;

				rgb[0]=imptr[0];
				rgb[1]=imptr[1];
				rgb[2]=imptr[2];
				imptr+=3;
				if(imptr>imend)//enc reload
				{
					ptrdiff_t nreq=0, nact=0;
					int k=0;
#ifdef LOUD
					double t2=0;
#endif

					imptr=imptr0;

					if(idx)
					{
						ptrdiff_t chunksize=0;
						
#ifdef LOUD
						t2=time_sec();
#endif
						*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;//flush
						*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;
						chunksize=streamptr-streambuf;
						csize+=fwrite(&chunksize, 1, 3, fdst);
						csize+=fwrite(streambuf, 1, chunksize, fdst);
						
						streamptr=streambuf;
						low=0;
						range=0xFFFFFFFFFFFF;
#ifdef LOUD
						twrite+=time_sec()-t2;
#endif
					}
					
#ifdef LOUD
					t2=time_sec();
#endif
					for(k=0;imptr<imend;)
						rgb[k++]=*imptr++;
					nreq=PBUFSIZE;
					if((size_t)nreq>(size_t)(usize-nread))
						nreq=usize-nread;
					nact=fread(imbuf, 1, nreq, fsrc);
					nread+=nact;
					if(nact!=nreq)
					{
						CRASH("Truncated PPM \"%s\"", srcfn);
						return 1;
					}
					imend=imbuf+nreq;
					imptr=imbuf;
					while(k<3)
						rgb[k++]=*imptr++;
#ifdef LOUD
					tread+=time_sec()-t2;
#endif
				}
				yuv[0]=rgb[1];
				yuv[1]=rgb[0];
				yuv[2]=rgb[2];
			}
			for(kc=0;kc<NCH;++kc)
			{
				int
					NW	=rows[1][0-1*NROWS*NCH*NVAL],
					N	=rows[1][0+0*NROWS*NCH*NVAL],
					W	=rows[0][0-1*NROWS*NCH*NVAL],
					eNEEE	=rows[1][1+3*NROWS*NCH*NVAL],
					eW	=rows[0][1-1*NROWS*NCH*NVAL];

				int pred=N+W-NW, ctx=FLOOR_LOG2(eW*eW+1);
				int vmax=N, vmin=W;
				if(ctx>NCTX-1)
					ctx=NCTX-1;
				ctx+=NCTX*kc;
				if(N<W)vmin=N, vmax=W;
#if 0
				int pred=1<<SHIFT>>1, j=0;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
#define PRED(E) estim[j++]=E;
				j=0;
				PREDLIST
#undef  PRED
				pred>>=SHIFT;
#endif
				CLAMP2(pred, vmin, vmax);
				pred+=offset;
				CLAMP2(pred, 0, 255);

				if(fwd)
				{
					error=yuv[kc]-pred;
					sym=(uint8_t)error;
					cdf=enccdf[ctx<<8|sym];
					freq=1<<PROBBITS;
					if(sym<255)
						freq=enccdf[ctx<<8|(sym+1)];
					freq-=cdf;
#ifdef FIFOVAL
					valfifo_enqueue(symptr[kc]);
#endif
					if(range<=0xFFFF)
					{
						*(uint32_t*)streamptr=(uint32_t)(low>>32);
						streamptr+=sizeof(uint32_t);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					low+=range*cdf>>PROBBITS;
					range=(range*freq>>PROBBITS)-1;
				}
				else
				{
					if(range<=0xFFFF)//stall: unpredictable branch
					{
#ifdef _MSC_VER
						if(streamptr>streambuf+SBUFSIZE-sizeof(uint16_t))
							CRASH("");
#endif
						code=code<<32|*(uint32_t*)streamptr;
						streamptr+=sizeof(uint32_t);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					int c2=(int)(((code-low)<<PROBBITS|((1ULL<<PROBBITS)-1))/range);//stall: DIV64
#ifdef _MSC_VER
					if(c2>>PROBBITS)
						CRASH("Dec error X %d  Y %d  C%d", kx, ky, kc);
#endif
					uint32_t info=CDF2sym[ctx<<PROBBITS|c2];//stall: cache miss
					sym=(uint8_t)info;
					cdf=info<<PROBBITS>>(32-PROBBITS);
					freq=info>>(PROBBITS+8);
#ifdef FIFOVAL
					valfifo_check(freq<<16|cdf);
#endif
					low+=range*cdf>>PROBBITS;
					range=(range*freq>>PROBBITS)-1;

					yuv[kc]=sym+pred;
					error=yuv[kc]-pred;
				}
				++hists[ctx<<8|sym];
				++hweight[ctx];
				if(!(hweight[ctx]&(hweight[ctx]-1)))
				{
					if(fwd)
						normalize_enc(ctx, hweight[ctx]>=0x1000);
					else
						normalize_dec(ctx, hweight[ctx]>=0x1000);
				}
				rows[0][0]=yuv[kc]-offset;
				rows[0][1]=(2*eW+((error<<1^error>>31)<<GRBITS)+eNEEE)>>2;
				offset=yuv[0];
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
			}
			if(!fwd)
			{
				uint8_t *imptr0=imptr;

				rgb[1]=yuv[0];
				rgb[0]=yuv[1];
				rgb[2]=yuv[2];
				imptr[0]=rgb[0];
				imptr[1]=rgb[1];
				imptr[2]=rgb[2];
#ifdef ENABLE_GUIDE
				if(memcmp(imptr, g_buf+idx, 3))
				{
					printf(
						"\n\n"
						"Guide error  X %5d  Y %5d\n"
						"  R  decoded 0x%02X %4d != original 0x%02X %4d\n"
						"  G  decoded 0x%02X %4d != original 0x%02X %4d\n"
						"  B  decoded 0x%02X %4d != original 0x%02X %4d\n"
						"\n"
						, kx
						, ky
						, imptr[0], imptr[0], g_buf[idx+0], g_buf[idx+0]
						, imptr[1], imptr[1], g_buf[idx+1], g_buf[idx+1]
						, imptr[2], imptr[2], g_buf[idx+2], g_buf[idx+2]
					);
					CRASH("");
				}
#endif
				imptr+=3;
				if(imptr>=imend-2)//dec flush
				{
					ptrdiff_t nact=0;
					int k=0;
#ifdef LOUD
					double t2=0;
#endif
					
#ifdef LOUD
					t2=time_sec();
#endif
					imptr=imptr0;
					for(k=0;imptr<imend;)
						*imptr++=rgb[k++];
					nact=fwrite(imbuf, 1, PBUFSIZE, fdst);
					if(nact!=PBUFSIZE)
					{
						CRASH("Cannot write to \"%s\"", dstfn);
						return 1;
					}
					imptr=imbuf;
					while(k<3)
						*imptr++=rgb[k++];
#ifdef LOUD
					twrite+=time_sec()-t2;
#endif

					//dec reload
					{
						ptrdiff_t chunksize=0, nact=0;
						
#ifdef LOUD
						t2=time_sec();
#endif
						nact=fread(&chunksize, 1, 3, fsrc);
						csize+=nact;
						if(nact!=3||(size_t)chunksize>SBUFSIZE)
						{
							CRASH("Corrupt archive \"%s\"", srcfn);
							return 1;
						}
						nact=fread(streambuf, 1, chunksize, fsrc);
						csize+=nact;
						if(nact!=chunksize)
						{
							CRASH("Truncated archive \"%s\"", srcfn);
							return 1;
						}
						streamptr=streambuf;
#ifdef LOUD
						tread+=time_sec()-t2;
#endif
						
						streamptr=streambuf;
						low=0;
						range=0xFFFFFFFFFFFF;
						code=0;
						code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);
						code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);
					}
				}
			}
		}
	}
	free(pixels);
	if(fwd)
	{
		if(imptr>imbuf)
		{
			ptrdiff_t chunksize=0;
#ifdef LOUD
			volatile double t2=0;
#endif
			
#ifdef LOUD
			t2=time_sec();
#endif
			*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;//flush
			*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;
			chunksize=streamptr-streambuf;
			csize+=fwrite(&chunksize, 1, 3, fdst);
			csize+=fwrite(streambuf, 1, chunksize, fdst);

			low=0;
			range=0xFFFFFFFFFFFF;
#ifdef LOUD
			twrite+=time_sec()-t2;
#endif
		}
	}
	else if(imptr>imbuf)
	{
		ptrdiff_t nreq=imptr-imbuf, nact=0;
		nact=fwrite(imbuf, 1, nreq, fdst);
		if(nact!=nreq)
		{
			CRASH("Cannot write to \"%s\"", dstfn);
			return 1;
		}
	}
	fclose(fsrc);
	fclose(fdst);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		double tother=t-twrite-tread;

		printf("WH %5d*%5d  \"%s\"\n", iw, ih, srcfn);
		printf("%9td->%9td  %8.4lf%%  %12.6lf:1  BPD %12.6lf\n"
			, usize
			, csize
			, 100.*csize/usize
			, (double)usize/csize
			, 8.*csize/usize
		);
		printf("%12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB  read\n"
			, tread
			, usize/(1024*1024*tread)
			, 1024*1024*1000*tread/usize
		);
		printf("%12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB  other\n"
			, tother
			, usize/(1024*1024*tother)
			, 1024*1024*1000*tother/usize
		);
		printf("%12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB  write\n"
			, twrite
			, usize/(1024*1024*twrite)
			, 1024*1024*1000*twrite/usize
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
	(void)csize;
	return 0;
}
