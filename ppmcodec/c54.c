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
//	#define PRINT_RCT
	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif


enum
{
	BUFSIZE=512*1024,

	NCTX=6,
	DEPTH=8,
	RLIMIT=11,
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
#define CLAMP2(X, L, H) X=X<(L)?L:X, X=X>(H)?H:X
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


static uint8_t rdbuf[BUFSIZE+sizeof(uint64_t[2])], wtbuf[BUFSIZE+sizeof(uint64_t[2])];


//cRCT
#if 1
static const int perms[]=
{
	1, 0, 2,
	1, 2, 0,

	//0, 1, 2,	2, 1, 0,
	//2, 0, 1,	1, 0, 2,
	//1, 2, 0,	0, 2, 1,
};
enum
{
	RCTBITS=2,
	RCTMAX=1<<RCTBITS,
	RCTLEVELS=RCTMAX+1,//+zero
	NRCTS=RCTLEVELS*RCTLEVELS*(RCTLEVELS+1)/2,
	NPERMS=_countof(perms)/3,
};
typedef struct _RCTInfo
{
	uint8_t pidx, uc0, vc0, vc1;
} RCTInfo;
static void print_rct(RCTInfo *rct, int tidx, int64_t score)
{
	printf("[%7d]  RCT%c%c%c_%c_%c%c/%c  %16lld"
		, tidx
		, '0'+perms[rct->pidx*3+0]
		, '0'+perms[rct->pidx*3+1]
		, '0'+perms[rct->pidx*3+2]
		, rct->uc0+(rct->uc0<9?'0':'A'-10)
		, rct->vc0+(rct->vc0<9?'0':'A'-10)
		, rct->vc1+(rct->vc1<9?'0':'A'-10)
		, RCTMAX+(RCTMAX<9?'0':'A'-10)
		, score
	);
}
typedef struct _AnalysisCtrs
{
	int64_t yctr, uctrs[RCTLEVELS], vctrs[RCTLEVELS*(RCTLEVELS+1)/2+1];
} AnalysisCtrs;
static void crct_analysis(FILE *fsrc, int iw, int ih, RCTInfo *ret_rct)
{
	enum
	{
		XPAD=8,
		NCH=3,
		NROWS=1,
		NVAL=1,
	};
	AnalysisCtrs counters[NPERMS]={0};
	uint8_t *const rdend=rdbuf+sizeof(uint64_t)+BUFSIZE, *rdptr=0;
	long fidx=0;
	int prev[3]={0}, rgb[3]={0}, yuv[3]={0};
	int64_t bestscore=0;
	RCTInfo rct={0};
#ifdef PRINT_RCT
	int it=0;
#endif
	int psize=0;
	int16_t *pixels=0;
	
	rdptr=rdend;
	fidx=ftell(fsrc);
	psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	for(int ky=0;ky<ih;++ky)
	{
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=0;kx<iw;++kx)
		{
			uint64_t data=*(uint64_t*)rdptr;
			if(rdptr>=rdend)
			{
				fread(rdbuf+sizeof(uint64_t), 1, BUFSIZE, fsrc);
				rdptr-=BUFSIZE;
				data|=*(uint64_t*)rdptr;
			}
			rdptr+=3;
			rgb[0]=(uint8_t)(data>>0*8);
			rgb[1]=(uint8_t)(data>>1*8);
			rgb[2]=(uint8_t)(data>>2*8);
			yuv[0]=rgb[0]-prev[0];
			yuv[1]=rgb[1]-prev[1];
			yuv[2]=rgb[2]-prev[2];
			prev[0]=rgb[0];
			prev[1]=rgb[1];
			prev[2]=rgb[2];
			rgb[0]=yuv[0]-rows[0][0];
			rgb[1]=yuv[1]-rows[0][1];
			rgb[2]=yuv[2]-rows[0][2];
			rows[0][0]=yuv[0];
			rows[0][1]=yuv[1];
			rows[0][2]=yuv[2];
			rows[0]+=NROWS*NVAL*NCH;
			rgb[0]<<=RCTBITS;
			rgb[1]<<=RCTBITS;
			rgb[2]<<=RCTBITS;
			for(int kp=0;kp<NPERMS;++kp)
			{
				AnalysisCtrs *currctrs=counters+kp;
				yuv[0]=rgb[perms[3*kp+0]];
				yuv[1]=rgb[perms[3*kp+1]];
				yuv[2]=rgb[perms[3*kp+2]];
				currctrs->yctr+=abs(yuv[0]);
				for(int kp1=0, idx=0;kp1<=RCTMAX;++kp1)
				{
					currctrs->uctrs[kp1]+=abs(yuv[1]-(kp1*yuv[0]>>RCTBITS));
					for(int kp2=0;kp1+kp2<=RCTMAX;++kp2, ++idx)
						currctrs->vctrs[idx]+=abs(yuv[2]-((kp1*yuv[0]+kp2*yuv[1])>>RCTBITS));
				}
			}
		}
	}
	fseek(fsrc, fidx, SEEK_SET);
	free(pixels);
	for(int kp=0;kp<NPERMS;++kp)
	{
		AnalysisCtrs *currctrs=counters+kp;
		for(int uc0=0;uc0<=RCTMAX;++uc0)
		{
			for(int vc0=0, idx=0;vc0<=RCTMAX;++vc0)
			{
				for(int vc1=0;vc0+vc1<=RCTMAX;++vc1, ++idx)
				{
					int64_t score=currctrs->yctr+currctrs->uctrs[uc0]+currctrs->vctrs[idx];
#ifdef PRINT_RCT
					{
						RCTInfo rct2={0};
						rct2.pidx=kp;
						rct2.uc0=uc0;
						rct2.vc0=vc0;
						rct2.vc1=vc1;
						print_rct(&rct2, it++, score);
						if(!bestscore||bestscore>score)
							printf(" <-");
						printf("\n");
					}
#endif
					if(!bestscore||bestscore>score)
					{
						bestscore=score;
						rct.pidx=kp;
						rct.uc0=uc0;
						rct.vc0=vc0;
						rct.vc1=vc1;
					}
				}
			}
		}
	}
	memcpy(ret_rct, &rct, sizeof(*ret_rct));
}
#endif


static uint8_t logtable[1<<DEPTH];
static uint32_t enctable[DEPTH<<DEPTH];
int c54_codec(int argc, char **argv)
{
	enum
	{
		XPAD=8,
		NCH=3,
		NROWS=4,
		NVAL=3,

		TOTALADD=4+RCTBITS,
	};
	const uint16_t tag='5'|'4'<<8;

	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	uint64_t c=0;
	int fwd=0, iw=0, ih=0;
	int64_t usize=0, csize=0;
	int psize=0;
	int16_t *pixels=0;
	uint64_t cache=0;
	uint32_t code=0, inc=0;
	int nbits=0;
	uint8_t *rdptr=0, *wtptr=0, *rdend=0, *wtend=0;
	RCTInfo rct={0};
	int ysh=0, ush=0, vsh=0;
	int csums[3]={0};
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
	if(fwd)
	{
		//for(int ks=0;ks<1<<DEPTH;++ks)
		//	logtable[ks]<<=8;
		for(int ks=0;ks<1<<DEPTH;++ks)
		{
			for(int kb=0;kb<DEPTH;++kb)
			{
				uint32_t code;
				int nzeros, stopbit, nbypass, codelen;
				
				nbypass=kb;
				nzeros=ks>>kb;
				stopbit=nzeros<RLIMIT;
				codelen=nzeros+1+kb;
				if(nzeros>=RLIMIT)
					nzeros=RLIMIT, stopbit=0, nbypass=DEPTH, codelen=RLIMIT+DEPTH;
				code=(ks&((1<<nbypass)-1))<<(nzeros+stopbit)|stopbit<<nzeros;
				enctable[kb<<DEPTH|ks]=code<<8|codelen;
			}
		}
	}
	memset(pixels, 0, psize);
	rdptr=rdbuf+sizeof(uint64_t)+BUFSIZE;
	wtptr=wtbuf+sizeof(uint64_t);
	rdend=rdbuf+sizeof(uint64_t)+BUFSIZE;
	wtend=wtbuf+sizeof(uint64_t)+BUFSIZE;
	csums[0]=0;
	csums[1]=rct.uc0;
	csums[2]=rct.vc0+rct.vc1;
	for(int ky=0;ky<ih;++ky)
	{
		int yuv[3]={0};
		int pred=0;
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=0;kx<iw;++kx)
		{
			int offset=0, csum=0;
			if(fwd)
			{
				uint64_t data=*(uint64_t*)rdptr;
				if(rdptr+3>=rdend)
				{
					fread(rdbuf+sizeof(uint64_t), 1, BUFSIZE, fsrc);
					rdptr-=BUFSIZE;
					data|=*(uint64_t*)rdptr;
				}
				rdptr+=3;
				yuv[0]=data>>ysh&255;
				yuv[1]=data>>ush&255;
				yuv[2]=data>>vsh&255;
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

				offset=0;
				if(kc==1)offset=rct.uc0*yuv[0];
				if(kc==2)offset=rct.vc0*yuv[0]+rct.vc1*yuv[1];
				csum=csums[kc];
				offset<<=4;
				if(csum<4)
					offset+=offset>>4;
#if 1
				(void)NNN;
				(void)NN;
				(void)NNE;
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
				nbypass=logtable[nW];
				//					1/8
				//					-2
				//				-5	10	2	1/8
				//	-2/8	1	-3	13	[?]
				int p1=13*W+10*N-5*NW+2*NE-2*NN-3*WW+WWW+((-2*WWWW+NNN+NEE)>>3);
				{
					int vmin, vmax;
					vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					vmin<<=4;
					vmax<<=4;
					CLAMP2(p1, vmin, vmax);
				}
				pred=(p1+offset+(1<<TOTALADD>>1))>>TOTALADD;
				CLAMP2(pred, 0, 255);
				offset=(offset-7)>>4;
				if(fwd)
				{
					nbypass<<=8;
					sym=(int8_t)(yuv[kc]-pred);
					sym=sym<<1^sym>>31;
					code=enctable[nbypass|sym];
					cache|=(uint64_t)code>>8<<nbits;
					nbits+=(uint8_t)code;
					*(uint64_t*)wtptr=cache;
					inc=nbits>>3;
					wtptr+=inc;
					if(wtptr>=wtend)
					{
						fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
						wtptr-=BUFSIZE;
						*(uint64_t*)(wtptr-inc)=cache;
					}
					cache>>=nbits&56;
					nbits&=7;
				}
				else
				{
					cache|=*(uint64_t*)rdptr<<nbits;
					inc=nbits>>3^7;
					rdptr+=inc;
					if(rdptr>=rdend)
					{
						fread(rdbuf+sizeof(uint64_t), 1, BUFSIZE, fsrc);
						rdptr-=BUFSIZE;
						cache|=*(uint64_t*)(rdptr-inc)<<nbits;
					}
					nbits|=56;
					nzeros=(int)TZCNT64(cache);
					sym=nzeros<<nbypass;
					if(nzeros>RLIMIT-1)
						nzeros=RLIMIT-1, nbypass=DEPTH, sym=0;
					cache>>=nzeros+1;
					sym|=(int)(cache&((1ULL<<nbypass)-1));
					cache>>=nbypass;
					nbits-=nzeros+1+nbypass;
					yuv[kc]=(uint8_t)((sym>>1^sym<<31>>31)+pred);
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
				uint64_t data=(uint64_t)yuv[2]<<vsh|(uint64_t)yuv[1]<<ush|(uint64_t)yuv[0]<<ysh;
				*(uint64_t*)wtptr=data;
				if(wtptr+3>=wtend)
				{
					fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
					wtptr-=BUFSIZE;
					*(uint64_t*)wtptr=data;
				}
				wtptr+=3;
			}
		}
	}
	if(fwd)
	{
		*(uint64_t*)wtptr=cache;
		if(wtptr>=wtend)
		{
			fwrite(wtbuf+sizeof(uint64_t), 1, BUFSIZE, fdst);
			wtptr-=BUFSIZE;
			*(uint64_t*)wtptr=cache;
		}
		++wtptr;
	}

	if(wtptr>wtbuf+sizeof(uint64_t))
		fwrite(wtbuf+sizeof(uint64_t), 1, wtptr-(wtbuf+sizeof(uint64_t)), fdst);
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
