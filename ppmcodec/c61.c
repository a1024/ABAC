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
#	define WIN32_LEAN_AND_MEAN
#	include<Windows.h>
#else
#	include<time.h>
#endif
#ifdef _MSC_VER
#	include<intrin.h>
#elif defined __GNUC__
#	include<x86intrin.h>
#endif


#ifdef _MSC_VER
	#define LOUD
	#define PRINT_RCT

	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif


//	#define ENABLE_LUMA_UPDATE	//X  bad
//	#define TEST_LUMA_UPDATE	//X  bad


enum
{
	BUFSIZE=512*1024,

	NCTX=7,
	DEPTH=8,
	RLIMIT=13,//<=16
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
	uint32_t yctrs[
#if defined ENABLE_LUMA_UPDATE
		RCTLEVELS*(RCTLEVELS+1)/2
#else
		1
#endif
	][
#if defined ENABLE_LUMA_UPDATE
		RCTLEVELS*(RCTLEVELS+1)/2
#else
		1
#endif
	], uctrs[RCTLEVELS], vctrs[RCTLEVELS*(RCTLEVELS+1)/2];
} AnalysisCtrs;
static void crct_analysis(FILE *fsrc, int iw, int ih, RCTInfo *ret_rct)
{
	enum
	{
		STRIDE=11,
	};
	AnalysisCtrs counters[NPERMS]={0};
	uint8_t *image=0, *imptr=0;
	int32_t size=0;
	int rstr=0;
	long fidx=0;
	uint32_t bestscore=0;
	RCTInfo rct={0};
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
	for(int ky=1;ky<=ih-STRIDE;ky+=STRIDE)
	{
		int kx=1;
		imptr=image+rstr+3*kx;
		fread(image, 1, size, fsrc);
		for(;kx<=iw-STRIDE;kx+=STRIDE)
		{
			int rgb[]=
			{
				(imptr[0]-imptr[-3+0]-imptr[-rstr+0]+imptr[-rstr-3+0])<<RCTBITS,
				(imptr[1]-imptr[-3+1]-imptr[-rstr+1]+imptr[-rstr-3+1])<<RCTBITS,
				(imptr[2]-imptr[-3+2]-imptr[-rstr+2]+imptr[-rstr-3+2])<<RCTBITS,
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
#if defined ENABLE_LUMA_UPDATE
				for(int kp1=0, idx=0;kp1<=RCTMAX;++kp1)
				{
					int u=yuv[1]-(kp1*yuv[0]>>RCTBITS);
					currctrs->uctrs[kp1]+=abs(u);
					for(int kp2=0;kp1+kp2<=RCTMAX;++kp2, ++idx)
					{
						int v=yuv[2]-((kp1*yuv[0]+kp2*yuv[1])>>RCTBITS);
						currctrs->vctrs[idx]+=abs(v);
						for(int kp3=0, idx2=0;kp3<=RCTMAX;++kp3)
						{
							for(int kp4=0;kp3+kp4<=RCTMAX;++kp4, ++idx2)
							{
								int y=yuv[0]+((kp3*u+kp4*v+(1<<RCTBITS))>>(RCTBITS+1));
							//	y<<=32-(8+RCTBITS);
							//	y>>=32-(8+RCTBITS);
								currctrs->yctrs[idx][idx2]+=abs(y);
							}
						}
					}
				}
#else
				currctrs->yctrs[0][0]+=abs(yuv[0]);
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
#endif
			}
		}
	}
	free(image);
	fseek(fsrc, fidx, SEEK_SET);
	
#ifdef PRINT_RCT
	printf("Y:\n");
	for(int k=0;k<RCTLEVELS*(RCTLEVELS+1)/2;++k)
	{
		for(int k2=0;k2<RCTLEVELS*(RCTLEVELS+1)/2;++k2)
		{
			for(int kp=0;kp<NPERMS;++kp)
				printf("  %12u", counters[kp].yctrs[k][k2]);
			printf("\n");
		}
		printf("\n");
	}
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
					for(int yc0=0, idx2=0;yc0<=RCTMAX;++yc0)
					{
						for(int yc1=0;yc0+yc1<=RCTMAX;++yc1, ++idx2)
						{
#if defined ENABLE_LUMA_UPDATE
							uint32_t score=currctrs->yctrs[idx][idx2]+currctrs->uctrs[uc0]+currctrs->vctrs[idx];
#else
							uint32_t score=currctrs->yctrs[0][0]+currctrs->uctrs[uc0]+currctrs->vctrs[idx];
#endif
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
							//	rct.yc0=yc0;
							//	rct.yc1=yc1;
#ifdef PRINT_RCT
								printf("%5d  ", it);
								print_rct(&rct);
								printf("  %10d  =  %10d + %10d + %10d\n"
									, score
									, currctrs->yctrs[idx][idx2]
									, currctrs->uctrs[uc0]
									, currctrs->vctrs[idx]
								);
#endif
							}
						}
					}
				}
			}
		}
	}
#ifdef PRINT_RCT
	print_rct(&rct);
	printf("\n");
#endif
	memcpy(ret_rct, &rct, sizeof(*ret_rct));
}
#endif


static uint8_t logtable[1<<DEPTH];
static uint32_t enctable[DEPTH<<DEPTH];
static uint8_t signpack[1<<DEPTH];
//static uint8_t clamptable[1024];
static uint8_t rdbuf[BUFSIZE+sizeof(uint32_t[2])], wtbuf[BUFSIZE+sizeof(uint32_t[2])];
int c61_codec(int argc, char **argv)
{
	enum
	{
		XPAD=8,
		NCH=3,
		NROWS=4,
		NVAL=3,

		PREDADD=4,
		TOTALADD=PREDADD+RCTBITS,
	};
	const uint16_t tag='6'|'1'<<8;

	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	uint32_t c=0;
	int fwd=0, iw=0, ih=0;
	int64_t usize=0, csize=0;
	int psize=0;
	int16_t *pixels=0;
	uint32_t cache=0;
	uint32_t inc=0;
	int nbits=0;
	uint8_t *rdptr=0, *wtptr=0, *rdend=0, *wtend=0;
	RCTInfo rct={0};
	int c0[3]={0}, c1[3]={0};
	int ysh=0, ush=0, vsh=0;
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
	psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)malloc(psize);
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
	ysh=perms[rct.pidx*3+0]*8;
	ush=perms[rct.pidx*3+1]*8;
	vsh=perms[rct.pidx*3+2]*8;
	for(int ks=0;ks<1<<DEPTH;++ks)
	{
		int val=31^LZCNT32(ks+1);
		if(val>NCTX-1)
			val=NCTX-1;
		logtable[ks]=val;
	}
	//for(int k=0;k<_countof(clamptable);++k)
	//{
	//	int val=k-(_countof(clamptable)*3>>3);
	//	CLAMP(val, 0, 255);
	//	clamptable[k]=val;
	//}
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
				codelen=nzeros+1+kb;
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
	for(int ky=0;ky<ih;++ky)
	{
		static const uint8_t masks[]={0, 1, 3, 7, 15, 31, 63, 127, 255};
		int yuv[3]={0};
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};

		//interleaved
#if 1
		uint32_t data;
		int nbypass[3], pred[3], offset[3];
		//int vmin[3], vmax[3], t[6];
		int sym[3];
		uint32_t code[3];
		
		//					1/8
		//					-2
		//				-5	10	2	1/8
		//	-2/8	1	-3	13	[?]
#if 0
#define PREDICT(CH) pred[CH]=\
	+10*rows[1][0+(+0*NCH+CH)*NROWS*NVAL]\
	+13*rows[0][0+(-1*NCH+CH)*NROWS*NVAL]\
	- 5*rows[1][0+(-1*NCH+CH)*NROWS*NVAL]\
	+ 2*rows[1][0+(+1*NCH+CH)*NROWS*NVAL]\
	- 2*rows[2][0+(+0*NCH+CH)*NROWS*NVAL]\
	- 3*rows[0][0+(-2*NCH+CH)*NROWS*NVAL]\
	+   rows[0][0+(-3*NCH+CH)*NROWS*NVAL]\
	+((\
		-2*rows[0][0+(-4*NCH+CH)*NROWS*NVAL]\
		+  rows[3][0+(+0*NCH+CH)*NROWS*NVAL]\
		+  rows[1][0+(+2*NCH+CH)*NROWS*NVAL]\
	)>>3)

#endif
		
		//					1/8
		//					-2
		//				-5	10	2
		//	-1/8	1	-3	13	[?]
#if 0
#define PREDICT(CH) pred[CH]=\
	+10*rows[1][0+(+0*NCH+CH)*NROWS*NVAL]\
	+13*rows[0][0+(-1*NCH+CH)*NROWS*NVAL]\
	- 5*rows[1][0+(-1*NCH+CH)*NROWS*NVAL]\
	+ 2*rows[1][0+(+1*NCH+CH)*NROWS*NVAL]\
	- 2*rows[2][0+(+0*NCH+CH)*NROWS*NVAL]\
	- 3*rows[0][0+(-2*NCH+CH)*NROWS*NVAL]\
	+   rows[0][0+(-3*NCH+CH)*NROWS*NVAL]\
	+((\
		-rows[0][0+(-4*NCH+CH)*NROWS*NVAL]\
		+rows[3][0+(+0*NCH+CH)*NROWS*NVAL]\
	)>>3)

#endif
		
		//				-2
		//			-5	10	2
		//	1	-3	13	[?]			opt
#if 1
#define PREDICT(CH) pred[CH]=\
	+ 5*rows[1][2+(+0*NCH+CH)*NROWS*NVAL]\
	+11*rows[0][0+(-1*NCH+CH)*NROWS*NVAL]\
	+ 1*rows[0][2+(-1*NCH+CH)*NROWS*NVAL]\
	+ 2*rows[1][0+(+1*NCH+CH)*NROWS*NVAL]\
	- 2*rows[2][0+(+0*NCH+CH)*NROWS*NVAL]\
	- 1*rows[0][2+(-2*NCH+CH)*NROWS*NVAL]\

#endif
		
		//				-2
		//			-5	10	2
		//	1	-3	13	[?]
#if 0
#define PREDICT(CH) pred[CH]=\
	+10*rows[1][0+(+0*NCH+CH)*NROWS*NVAL]\
	+13*rows[0][0+(-1*NCH+CH)*NROWS*NVAL]\
	- 5*rows[1][0+(-1*NCH+CH)*NROWS*NVAL]\
	+ 2*rows[1][0+(+1*NCH+CH)*NROWS*NVAL]\
	- 2*rows[2][0+(+0*NCH+CH)*NROWS*NVAL]\
	- 3*rows[0][0+(-2*NCH+CH)*NROWS*NVAL]\
	+   rows[0][0+(-3*NCH+CH)*NROWS*NVAL]\

#endif

		//			-2
		//		-4	10	2
		//	-2	12	[?]
#if 0
#define PREDICT(CH) pred[CH]=\
	+10*rows[1][0+(+0*NCH+CH)*NROWS*NVAL]\
	+12*rows[0][0+(-1*NCH+CH)*NROWS*NVAL]\
	- 4*rows[1][0+(-1*NCH+CH)*NROWS*NVAL]\
	+ 2*rows[1][0+(+1*NCH+CH)*NROWS*NVAL]\
	- 2*rows[2][0+(+0*NCH+CH)*NROWS*NVAL]\
	- 2*rows[0][0+(-2*NCH+CH)*NROWS*NVAL]\

#endif
		
		//			-1
		//		-2	+5	1
		//	-1	+6	[?]
#if 0
#define PREDICT(CH) pred[CH]=(\
		+5*rows[1][0+(+0*NCH+CH)*NROWS*NVAL]\
		+6*rows[0][0+(-1*NCH+CH)*NROWS*NVAL]\
		-2*rows[1][0+(-1*NCH+CH)*NROWS*NVAL]\
		+1*rows[1][0+(+1*NCH+CH)*NROWS*NVAL]\
		-1*rows[2][0+(+0*NCH+CH)*NROWS*NVAL]\
		-1*rows[0][0+(-2*NCH+CH)*NROWS*NVAL]\
	)<<1

#endif
		
		//	-2	+3
		//	+3	[?]
#if 0
#define PREDICT(CH) pred[CH]=(\
		+3*rows[1][0+(+0*NCH+CH)*NROWS*NVAL]\
		+3*rows[0][0+(-1*NCH+CH)*NROWS*NVAL]\
		-2*rows[1][0+(-1*NCH+CH)*NROWS*NVAL]\
	)<<2

#endif
		
		//			8	12	8
		//	20	[17]	?
#if 1
#define UPDATE(CH) nbypass[CH]=\
	rows[0][1+(+0*NCH+CH)*NROWS*NVAL]=(\
		+20*rows[0][1+(-1*NCH+CH)*NROWS*NVAL]\
		+17*sym[CH]\
		+ 8*rows[1][1+(+1*NCH+CH)*NROWS*NVAL]\
		+12*rows[1][1+(+2*NCH+CH)*NROWS*NVAL]\
		+ 8*rows[1][1+(+3*NCH+CH)*NROWS*NVAL]\
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
			for(int kx=0;kx<iw;++kx)
			{
#if 0
				t[0]=vmax[0]=rows[1][0+(+0*NCH+0)*NROWS*NVAL]; t[0+3]=vmin[0]=rows[0][0+(-1*NCH+0)*NROWS*NVAL];//N W
				t[1]=vmax[1]=rows[1][0+(+0*NCH+1)*NROWS*NVAL]; t[1+3]=vmin[1]=rows[0][0+(-1*NCH+1)*NROWS*NVAL];
				t[2]=vmax[2]=rows[1][0+(+0*NCH+2)*NROWS*NVAL]; t[2+3]=vmin[2]=rows[0][0+(-1*NCH+2)*NROWS*NVAL];
				vmin[0]=t[0]<t[0+3]?t[0]:vmin[0]; vmax[0]=t[0]<t[0+3]?t[0+3]:vmax[0];
				vmin[1]=t[1]<t[1+3]?t[1]:vmin[1]; vmax[1]=t[1]<t[1+3]?t[1+3]:vmax[1];
				vmin[2]=t[2]<t[2+3]?t[2]:vmin[2]; vmax[2]=t[2]<t[2+3]?t[2+3]:vmax[2];
				t[0]=rows[1][0+(+1*NCH+0)*NROWS*NVAL];//NE
				t[1]=rows[1][0+(+1*NCH+1)*NROWS*NVAL];
				t[2]=rows[1][0+(+1*NCH+2)*NROWS*NVAL];
				t[3]=rows[1][0+(+3*NCH+0)*NROWS*NVAL];//NEEE
				t[4]=rows[1][0+(+3*NCH+1)*NROWS*NVAL];
				t[5]=rows[1][0+(+3*NCH+2)*NROWS*NVAL];
				vmin[0]=vmin[0]>t[0]?t[0]:vmin[0];
				vmin[1]=vmin[1]>t[1]?t[1]:vmin[1];
				vmin[2]=vmin[2]>t[2]?t[2]:vmin[2];
				vmax[0]=vmax[0]<t[0]?t[0]:vmax[0];
				vmax[1]=vmax[1]<t[1]?t[1]:vmax[1];
				vmax[2]=vmax[2]<t[2]?t[2]:vmax[2];
				vmin[0]=vmin[0]>t[0+3]?t[0+3]:vmin[0];
				vmin[1]=vmin[1]>t[1+3]?t[1+3]:vmin[1];
				vmin[2]=vmin[2]>t[2+3]?t[2+3]:vmin[2];
				vmax[0]=vmax[0]<t[0+3]?t[0+3]:vmax[0];
				vmax[1]=vmax[1]<t[1+3]?t[1+3]:vmax[1];
				vmax[2]=vmax[2]<t[2+3]?t[2+3]:vmax[2];
				vmin[0]<<=PREDADD;
				vmax[0]<<=PREDADD;
				vmin[1]<<=PREDADD;
				vmax[1]<<=PREDADD;
				vmin[2]<<=PREDADD;
				vmax[2]<<=PREDADD;
				CLAMP(pred[0], vmin[0]-(1<<PREDADD>>1), vmax[0]+(1<<PREDADD>>1));
				CLAMP(pred[1], vmin[1]-(1<<PREDADD>>1), vmax[1]+(1<<PREDADD>>1));
				CLAMP(pred[2], vmin[2]-(1<<PREDADD>>1), vmax[2]+(1<<PREDADD>>1));
#endif				
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
				//pred[0]=(clamptable+(_countof(clamptable)*3>>3))[pred[0]];
				//pred[1]=(clamptable+(_countof(clamptable)*3>>3))[pred[1]];
				//pred[2]=(clamptable+(_countof(clamptable)*3>>3))[pred[2]];
				CLAMP(pred[0], 0, 255);
				CLAMP(pred[1], 0, 255);
				CLAMP(pred[2], 0, 255);

				sym[0]=(int8_t)(yuv[0]-pred[0]);
				sym[1]=(int8_t)(yuv[1]-pred[1]);
				sym[2]=(int8_t)(yuv[2]-pred[2]);
				
#if defined ENABLE_LUMA_UPDATE || defined TEST_LUMA_UPDATE
				sym[1]-=(sym[0]+(1<<5>>1))>>5;		sym[1]=(int8_t)sym[1];
				sym[2]-=(sym[0]+(1<<5>>1))>>5;		sym[2]=(int8_t)sym[2];
				sym[0]+=(sym[1]+sym[2]+(1<<6>>1))>>6;	sym[0]=(int8_t)sym[0];

			//	sym[0]-=(rct.uc0*sym[1]+(rct.vc0+rct.vc1)*sym[2]+16)>>5;//
			//	sym[0]+=(3*sym[1]+sym[2]+16)>>5;
			//	sym[0]+=(rct.yc0*sym[1]+rct.yc1*sym[2]+(1<<RCTBITS))>>(RCTBITS+1);//
			//	sym[0]=(int8_t)sym[0];
#endif

				code[0]=(enctable+(1<<DEPTH>>1))[(nbypass[0]<<DEPTH)+sym[0]];
				code[1]=(enctable+(1<<DEPTH>>1))[(nbypass[1]<<DEPTH)+sym[1]];
				code[2]=(enctable+(1<<DEPTH>>1))[(nbypass[2]<<DEPTH)+sym[2]];
				sym[0]=(signpack+(1<<DEPTH>>1))[sym[0]];
				sym[1]=(signpack+(1<<DEPTH>>1))[sym[1]];
				sym[2]=(signpack+(1<<DEPTH>>1))[sym[2]];

				for(int kc=0;kc<3;++kc)
				{
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
				rows[0][0+(+0*NCH+0)*NROWS*NVAL]=(yuv[0]<<RCTBITS)-((offset[0]-7)>>PREDADD);
				rows[0][0+(+0*NCH+1)*NROWS*NVAL]=(yuv[1]<<RCTBITS)-((offset[1]-7)>>PREDADD);
				rows[0][0+(+0*NCH+2)*NROWS*NVAL]=(yuv[2]<<RCTBITS)-((offset[2]-7)>>PREDADD);
				UPDATE(0);
				UPDATE(1);
				UPDATE(2);
				rows[0][2+(+0*NCH+0)*NROWS*NVAL]=2*rows[0][0+(+0*NCH+0)*NROWS*NVAL]-rows[0][0+(-1*NCH+0)*NROWS*NVAL];
				rows[0][2+(+0*NCH+1)*NROWS*NVAL]=2*rows[0][0+(+0*NCH+1)*NROWS*NVAL]-rows[0][0+(-1*NCH+1)*NROWS*NVAL];
				rows[0][2+(+0*NCH+2)*NROWS*NVAL]=2*rows[0][0+(+0*NCH+2)*NROWS*NVAL]-rows[0][0+(-1*NCH+2)*NROWS*NVAL];
				rows[0]+=NCH*NROWS*NVAL;
				rows[1]+=NCH*NROWS*NVAL;
				rows[2]+=NCH*NROWS*NVAL;
				rows[3]+=NCH*NROWS*NVAL;
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

			for(int kx=0;kx<iw;++kx)
			{
#if 0
				t[0]=vmax[0]=rows[1][0+(+0*NCH+0)*NROWS*NVAL]; t[0+3]=vmin[0]=rows[0][0+(-1*NCH+0)*NROWS*NVAL];//N W
				t[1]=vmax[1]=rows[1][0+(+0*NCH+1)*NROWS*NVAL]; t[1+3]=vmin[1]=rows[0][0+(-1*NCH+1)*NROWS*NVAL];
				t[2]=vmax[2]=rows[1][0+(+0*NCH+2)*NROWS*NVAL]; t[2+3]=vmin[2]=rows[0][0+(-1*NCH+2)*NROWS*NVAL];
				vmin[0]=t[0]<t[0+3]?t[0]:vmin[0]; vmax[0]=t[0]<t[0+3]?t[0+3]:vmax[0];
				vmin[1]=t[1]<t[1+3]?t[1]:vmin[1]; vmax[1]=t[1]<t[1+3]?t[1+3]:vmax[1];
				vmin[2]=t[2]<t[2+3]?t[2]:vmin[2]; vmax[2]=t[2]<t[2+3]?t[2+3]:vmax[2];
				t[0]=rows[1][0+(+1*NCH+0)*NROWS*NVAL];//NE
				t[1]=rows[1][0+(+1*NCH+1)*NROWS*NVAL];
				t[2]=rows[1][0+(+1*NCH+2)*NROWS*NVAL];
				t[3]=rows[1][0+(+3*NCH+0)*NROWS*NVAL];//NEEE
				t[4]=rows[1][0+(+3*NCH+1)*NROWS*NVAL];
				t[5]=rows[1][0+(+3*NCH+2)*NROWS*NVAL];
				vmin[0]=vmin[0]>t[0]?t[0]:vmin[0];
				vmin[1]=vmin[1]>t[1]?t[1]:vmin[1];
				vmin[2]=vmin[2]>t[2]?t[2]:vmin[2];
				vmax[0]=vmax[0]<t[0]?t[0]:vmax[0];
				vmax[1]=vmax[1]<t[1]?t[1]:vmax[1];
				vmax[2]=vmax[2]<t[2]?t[2]:vmax[2];
				vmin[0]=vmin[0]>t[0+3]?t[0+3]:vmin[0];
				vmin[1]=vmin[1]>t[1+3]?t[1+3]:vmin[1];
				vmin[2]=vmin[2]>t[2+3]?t[2+3]:vmin[2];
				vmax[0]=vmax[0]<t[0+3]?t[0+3]:vmax[0];
				vmax[1]=vmax[1]<t[1+3]?t[1+3]:vmax[1];
				vmax[2]=vmax[2]<t[2+3]?t[2+3]:vmax[2];
				vmin[0]<<=PREDADD;
				vmax[0]<<=PREDADD;
				vmin[1]<<=PREDADD;
				vmax[1]<<=PREDADD;
				vmin[2]<<=PREDADD;
				vmax[2]<<=PREDADD;
				CLAMP(pred[0], vmin[0]-(1<<PREDADD>>1), vmax[0]+(1<<PREDADD>>1));
				CLAMP(pred[1], vmin[1]-(1<<PREDADD>>1), vmax[1]+(1<<PREDADD>>1));
				CLAMP(pred[2], vmin[2]-(1<<PREDADD>>1), vmax[2]+(1<<PREDADD>>1));
#endif
				for(int kc=0;kc<3;++kc)
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
				
#if defined ENABLE_LUMA_UPDATE || defined TEST_LUMA_UPDATE
				sym[0]-=(sym[1]+sym[2]+(1<<6>>1))>>6;		sym[0]=(int8_t)sym[0];
				sym[2]+=(sym[0]+(1<<5>>1))>>5;			sym[2]=(int8_t)sym[2];
				sym[1]+=(sym[0]+(1<<5>>1))>>5;			sym[1]=(int8_t)sym[1];

			//	sym[0]+=(rct.uc0*sym[1]+(rct.vc0+rct.vc1)*sym[2]+16)>>5;//
			//	sym[0]-=(3*sym[1]+sym[2]+16)>>5;
			//	sym[0]-=(rct.yc0*sym[1]+rct.yc1*sym[2]+(1<<RCTBITS))>>(RCTBITS+1);
			//	sym[0]=(int8_t)sym[0];
#endif

				offset[0]=0;
				pred[0]=(pred[0]+offset[0]+((1<<TOTALADD>>1)+2))>>TOTALADD;
				//pred[0]=(clamptable+(_countof(clamptable)*3>>3))[pred[0]];
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
				//pred[1]=(clamptable+(_countof(clamptable)*3>>3))[pred[1]];
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
				//pred[2]=(clamptable+(_countof(clamptable)*3>>3))[pred[2]];
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
				data=(uint32_t)yuv[2]<<vsh|(uint32_t)yuv[1]<<ush|(uint32_t)yuv[0]<<ysh;
				*(uint32_t*)wtptr=data;
				if(wtptr>=wtend)
				{
					fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
					wtptr-=BUFSIZE;
					*(uint32_t*)wtptr=data;
				}
				wtptr+=3;
				rows[0][0+(+0*NCH+0)*NROWS*NVAL]=(yuv[0]<<RCTBITS)-((offset[0]-7)>>PREDADD);
				rows[0][0+(+0*NCH+1)*NROWS*NVAL]=(yuv[1]<<RCTBITS)-((offset[1]-7)>>PREDADD);
				rows[0][0+(+0*NCH+2)*NROWS*NVAL]=(yuv[2]<<RCTBITS)-((offset[2]-7)>>PREDADD);
				rows[0][2+(+0*NCH+0)*NROWS*NVAL]=2*rows[0][0+(+0*NCH+0)*NROWS*NVAL]-rows[0][0+(-1*NCH+0)*NROWS*NVAL];
				rows[0][2+(+0*NCH+1)*NROWS*NVAL]=2*rows[0][0+(+0*NCH+1)*NROWS*NVAL]-rows[0][0+(-1*NCH+1)*NROWS*NVAL];
				rows[0][2+(+0*NCH+2)*NROWS*NVAL]=2*rows[0][0+(+0*NCH+2)*NROWS*NVAL]-rows[0][0+(-1*NCH+2)*NROWS*NVAL];
				rows[0]+=NCH*NROWS*NVAL;
				rows[1]+=NCH*NROWS*NVAL;
				rows[2]+=NCH*NROWS*NVAL;
				rows[3]+=NCH*NROWS*NVAL;
				PREDICT(0);
				PREDICT(1);
				PREDICT(2);
				nbypass[0]=logtable[nbypass[0]];//next nW
				nbypass[1]=logtable[nbypass[1]];
				nbypass[2]=logtable[nbypass[2]];
			}
		}
#endif
		//reference
#if 0
		int pred=0;
		for(int kx=0;kx<iw;++kx)
		{
			int offset=0;
			if(fwd)
			{
				uint32_t data=*(uint32_t*)rdptr;
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
			}
			for(int kc=0;kc<3;++kc)
			{
				int
					NNN	=rows[3][0+0*NCH*NROWS*NVAL],
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL],
					NWW	=rows[1][0-2*NCH*NROWS*NVAL],
					NW	=rows[1][0-1*NCH*NROWS*NVAL],
					N	=rows[1][0+0*NCH*NROWS*NVAL],
					NE	=rows[1][0+1*NCH*NROWS*NVAL],
					NEE	=rows[1][0+2*NCH*NROWS*NVAL],
					NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
					NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
					WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
					WWW	=rows[0][0-3*NCH*NROWS*NVAL],
					WW	=rows[0][0-2*NCH*NROWS*NVAL],
					W	=rows[0][0-1*NCH*NROWS*NVAL],
					
					nNN	=rows[2][1+0*NCH*NROWS*NVAL],
					nNNE	=rows[2][1+1*NCH*NROWS*NVAL],
					nNNEE	=rows[2][1+2*NCH*NROWS*NVAL],
					nNNEEE	=rows[2][1+3*NCH*NROWS*NVAL],
					nNW	=rows[1][1-1*NCH*NROWS*NVAL],
					nN	=rows[1][1+0*NCH*NROWS*NVAL],
					nNE	=rows[1][1+1*NCH*NROWS*NVAL],
					nNEE	=rows[1][1+2*NCH*NROWS*NVAL],
					nNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
					nNEEEE	=rows[1][1+4*NCH*NROWS*NVAL],
					nNEEEEE	=rows[1][1+5*NCH*NROWS*NVAL],
					nWW	=rows[0][1-2*NCH*NROWS*NVAL],
					nW	=rows[0][1-1*NCH*NROWS*NVAL];
				int nbypass, sym, nzeros;

#if 1
				(void)NNN;
				(void)NN;
				(void)NNE;
				(void)NWW;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)NEEEE;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;

				(void)nNN;
				(void)nNNE;
				(void)nNNEE;
				(void)nNNEEE;
				(void)nNW;
				(void)nN;
				(void)nNE;
				(void)nNEE;
				(void)nNEEE;
				(void)nNEEEE;
				(void)nNEEEEE;
				(void)nWW;
				(void)nW;
#endif
				offset=c0[kc]*yuv[0]+c1[kc]*yuv[1];
				if(c0[kc]+c1[kc]<1<<TOTALADD)
					offset-=offset>>6;
				nbypass=logtable[nW];
				//					1/8
				//					-2
				//				-5	10	2	1/8
				//	-2/8	1	-3	13	[?]
				pred=13*W+10*N-5*NW+2*NE-2*NN-3*WW+WWW+((-2*WWWW+NNN+NEE)>>3);
				{
					int vmin, vmax;
					vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					vmin<<=PREDADD;
					vmax<<=PREDADD;
					CLAMP(pred, vmin-(1<<PREDADD>>1), vmax+(1<<PREDADD>>1));
				}
				pred=(pred+offset+((1<<TOTALADD>>1)+2))>>TOTALADD;
				CLAMP(pred, 0, 255);
				offset=(offset-7)>>PREDADD;
				if(fwd)
				{
					sym=(int8_t)(yuv[kc]-pred);
					code=(enctable+(1<<DEPTH>>1))[(nbypass<<DEPTH)+sym];
					sym=(signpack+(1<<DEPTH>>1))[sym];
					//sym=sym<<1^sym>>31;
					cache|=code>>8<<nbits;
					nbits+=(uint8_t)code;
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
				else
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
					sym=nzeros<<nbypass;
					if(nzeros>RLIMIT-1)
						nzeros=RLIMIT-1, nbypass=DEPTH, sym=0;
					cache>>=nzeros+1;
					sym|=(int)(cache&masks[nbypass]);
				//	sym|=(int)(cache&((1ULL<<nbypass)-1));
					cache>>=nbypass;
					nbits-=nzeros+1+nbypass;
					yuv[kc]=(uint8_t)(signpack[sym]+pred);
				//	yuv[kc]=(uint8_t)((sym>>1^sym<<31>>31)+pred);
#ifdef ENABLE_GUIDE
					if(yuv[kc]!=g_image[3*(iw*ky+kx)+perms[rct.pidx*3+kc]])
					{
						printf("Y %d  X %d  C%d\n", ky, kx, kc);
						CRASH("guide");
						return 1;
					}
#endif
				}
				rows[0][0]=(yuv[kc]<<RCTBITS)-offset;
				//			8	12	8
				//	20	[17]	?
				rows[0][1]=(20*nW+17*sym+8*nNE+12*nNEE+8*nNEEE+9)>>6;
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
			}
			if(!fwd)
			{
				uint32_t data=(uint32_t)yuv[2]<<vsh|(uint32_t)yuv[1]<<ush|(uint32_t)yuv[0]<<ysh;
				*(uint32_t*)wtptr=data;
				if(wtptr>=wtend)
				{
					fwrite(wtbuf+sizeof(uint32_t), 1, BUFSIZE, fdst);
					wtptr-=BUFSIZE;
					*(uint32_t*)wtptr=data;
				}
				wtptr+=3;
			}
		}
#endif
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
	}

	if(wtptr>wtbuf+sizeof(uint32_t))
		fwrite(wtbuf+sizeof(uint32_t), 1, wtptr-(wtbuf+sizeof(uint32_t)), fdst);
	free(pixels);
	fclose(fsrc);
	fclose(fdst);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		struct stat info={0};
		stat(dstfn, &info);
		csize=info.st_size;
		printf("CWH=3*%d*%d  \"%s\"\n", iw, ih, srcfn);
		printf("%10td->%10td  %8.4lf%%  %12.6lf:1  BPD %12.6lf\n"
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
