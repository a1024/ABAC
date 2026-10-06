//gcc -Wall -O3 -mavx2 -mlzcnt -static -s -Dc62_codec=main c62.c -o c62
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
	#define PRINT_ANALYSIS
	#define TRAIN_HIST

	#define ENABLE_GUIDE
//	#define FIFOVAL
#endif


	#define USE_AGETABLE



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
#elif defined __GNUC__//zero -> UB
#define LZCNT32 __builtin_clz
#define LZCNT64 __builtin_clzll
#define TZCNT32 __builtin_ctz
#define TZCNT64 __builtin_ctzll
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
static void guide_save(const uint8_t *image, int iw, int ih)
{
	int64_t size=(int64_t)3*iw*ih;

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
static void guide_save_file(FILE *f, int iw, int ih)
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
#endif//FIFOVAL
#endif//1


#define PREDLIST\
	PRED(W)\
	PRED(NE)\
	PRED(2*N-NN)\
	PRED(N+W-NW)\

static const uint8_t perms[]=
{
	0, 1, 2,	2, 1, 0,
	2, 0, 1,	1, 0, 2,
	1, 2, 0,	0, 2, 1,
};
enum
{
	RCTBITS=2,
	RCTMAX=1<<RCTBITS,
	RCTLEVELS=RCTMAX+1,//+unity
	//NRCTS=RCTLEVELS*RCTLEVELS*(RCTLEVELS+1)/2,
	NPERMS=_countof(perms)/3,
	OCH_COUNT=NPERMS*(1+RCTLEVELS+RCTLEVELS*(RCTLEVELS+1)/2),
	
	SHIFT=18,
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED

	XPAD=8,
	NROWS=3,
	NCH=3,
	NVAL=2,

	GRBITS=5,
	CTXMIN=9,
	CTXMAX=18,
	NCTX=CTXMAX+2-CTXMIN,

	UDEPTH=16,
	BDEPTH=8,
	
	USEBITS=10,
	STOREBITS=20,
	AGEBITS=9,
	AGEMAX=(1<<AGEBITS)-1,
};
typedef struct _RCTInfo
{
	uint8_t pidx, uc0, vc0, vc1;
} RCTInfo;

typedef struct _AC
{
	uint64_t lo, hi, code;
	uint8_t *ptr, *end, *start;

	int64_t byteidx;
} AC;
typedef int16_t Pixel_t;
typedef int32_t Cell_t;
static Cell_t stats[NCH][NCTX][256][UDEPTH+BDEPTH];
#ifdef USE_AGETABLE
static uint16_t agetable[1<<AGEBITS];
#endif

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
static void crct_analysis(uint8_t *image, int iw, int ih, RCTInfo *ret_rct)
{
	enum
	{
		STRIDE=7,
	};
	int64_t counters[OCH_COUNT]={0}, bestsum=0;
	RCTInfo rct={0};
	int rstr=0;
#ifdef LOUD
	double t=time_sec();
#endif

	rstr=3*iw;
	for(int ky=1;ky<=ih-STRIDE;ky+=STRIDE)
	{
		uint8_t *imptr=image+rstr*ky+3;
		for(int kx=1;kx<=iw-STRIDE;kx+=STRIDE, imptr+=3*STRIDE)
		{
			int64_t *ctrptr=counters;
			int rgb[]=
			{
				(imptr[0]-imptr[-3+0]-imptr[-rstr+0]+imptr[-rstr-3+0])<<RCTBITS,
				(imptr[1]-imptr[-3+1]-imptr[-rstr+1]+imptr[-rstr-3+1])<<RCTBITS,
				(imptr[2]-imptr[-3+2]-imptr[-rstr+2]+imptr[-rstr-3+2])<<RCTBITS,
			};
			for(int kp=0;kp<NPERMS;++kp)
			{
				int yuv[]=
				{
					rgb[perms[3*kp+0]],
					rgb[perms[3*kp+1]],
					rgb[perms[3*kp+2]],
				};
				*ctrptr++ += abs(yuv[0]);
				for(int uc0=0;uc0<=RCTMAX;++uc0)
					*ctrptr++ += abs(yuv[1]-(yuv[0]*uc0>>RCTBITS));
				for(int vc0=0;vc0<=RCTMAX;++vc0)
				{
					for(int vc1=0;vc0+vc1<=RCTMAX;++vc1)
						*ctrptr++ += abs(yuv[2]-((yuv[0]*vc0+yuv[1]*vc1)>>RCTBITS));
				}
			}
		}
	}

	bestsum=0;
	for(int kp=0;kp<NPERMS;++kp)
	{
		int64_t *currctrs=counters+(1+RCTLEVELS+RCTLEVELS*(RCTLEVELS+1)/2)*kp;
		for(int uc0=0;uc0<=RCTMAX;++uc0)
		{
			for(int vc0=0, idx=0;vc0<=RCTMAX;++vc0)
			{
				for(int vc1=0;vc0+vc1<=RCTMAX;++vc1, ++idx)
				{
					int64_t sum=currctrs[0]+currctrs[1+uc0]+currctrs[1+RCTLEVELS+idx];
					if(!bestsum||bestsum>sum)
					{
						bestsum=sum;
						rct.pidx=kp;
						rct.uc0=uc0;
						rct.vc0=vc0;
						rct.vc1=vc1;
#ifdef PRINT_ANALYSIS
						print_rct(&rct);
						printf("  %16lld\n", sum);
#endif
					}
				}
			}
		}
	}
	*ret_rct=rct;
#ifdef LOUD
	t=time_sec()-t;
	{
		int64_t usize=(int64_t)3*iw*ih;
		printf("A  %12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB\n"
			, t
			, usize/(t*1024*1024)
			, t*1024*1024*1000/usize
		);
	}
#endif
}
INLINE int codebit(AC *ac, Cell_t *pcell, int bit, const int fwd)
{
	uint64_t x;
	int32_t prob, count, p1, sh;
	Cell_t cell;
	
#ifdef USE_AGETABLE
	cell=*pcell;
	x=ac->hi-ac->lo;
	sh=p1=(uint32_t)(cell>>(STOREBITS-USEBITS+AGEBITS));
	count=agetable[cell&AGEMAX];
	CLAMP(sh, -48, 48);
	prob=(uint32_t)(cell>>AGEBITS);
	p1-=sh>>4;
	sh=count&15;
	count>>=4;
	p1+=1<<USEBITS>>1;
#else
	cell=*pcell;
	x=ac->hi-ac->lo;
	sh=p1=(uint32_t)(cell>>(STOREBITS-USEBITS+AGEBITS));
	count=(uint32_t)(cell&AGEMAX);
	CLAMP(sh, -48, 48);
	prob=(uint32_t)(cell>>AGEBITS);
	p1-=sh>>4;
	sh=31^LZCNT32(count+10);
	count+=count<AGEMAX;
	p1+=1<<USEBITS>>1;
#endif
	if(x<=0xFFFF)
	{
		if(ac->ptr>=ac->end)
		{
#ifdef _MSC_VER
			int64_t usize=ac->ptr-ac->start;

			printf("\n%d/%d  %8.4lf%%\n"
				, (int32_t)usize
				, (int32_t)ac->byteidx
				, 100.*usize/ac->byteidx
			);
#endif
			CRASH("inflation");
		}
		if(fwd)
			*(uint32_t*)ac->ptr=(uint32_t)(ac->lo>>32);
		else
			ac->code=ac->code<<32|*(uint32_t*)ac->ptr;
		ac->ptr+=4;
		ac->lo<<=32;
		ac->hi=ac->hi<<32|~0U;
		if(ac->hi<ac->lo)ac->hi=~0ULL;
		x=ac->hi-ac->lo;
	}
	x=ac->lo+(x*(uint32_t)p1>>USEBITS);
	bit=fwd?bit:ac->code<x;
#ifdef _MSC_VER
	if((uint32_t)(p1-1)>(uint32_t)((1<<USEBITS)-2))
		CRASH("Invalid p1 0x%08X / %d bit", p1, USEBITS);
#ifdef FIFOVAL
	if(fwd)
		fifoval_enqueue(bit<<31^p1);
	else
		fifoval_check(bit<<31^p1);
#endif
#endif
	*(bit?&ac->hi:&ac->lo)=x-bit;
	prob+=(((2*bit-1)<<(STOREBITS-1))-prob)>>sh;
	*pcell=prob<<AGEBITS|count;
	return bit;
}
int c62_codec(int argc, char **argv)
{
	const uint16_t tag='6'|'2'<<8;

	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	uint32_t c=0;
	int fwd=0, iw=0, ih=0;
	int64_t usize=0, csize=0, streamsize=0;
	uint8_t *image=0, *imptr=0, *stream=0;
	int psize=0;
	Pixel_t *pixels=0;
	RCTInfo rct={0};
	const uint8_t *perm=0;
	ALIGN(16) int coeffs[3][NPREDS]={0}, bias[3]={0};
	int yuv[3]={0};
	AC ac={0};
#ifdef TRAIN_HIST
	enum
	{
	//	THSIZE=11,
		THSIZE=3,
	};
	int thist[3][THSIZE]={0};
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
		csize=(int64_t)3*iw*ih;
	}
	else
	{
		struct stat info={0};

		stat(srcfn, &info);
		csize=info.st_size;
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
	image=(uint8_t*)malloc(usize);
	stream=(uint8_t*)malloc(csize);
	psize=(iw+2*XPAD)*(int)sizeof(Pixel_t[NCH*NROWS*NVAL]);
	pixels=(Pixel_t*)malloc(psize);
	if(!image||!stream||!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	streamsize=usize;
	if(fwd)
	{
		fread(image, 1, usize, fsrc);
		guide_save(image, iw, ih);
	}
	else
	{
		streamsize=csize-ftell(fsrc);
		fread(stream, 1, streamsize, fsrc);
	}
	fclose(fsrc);

	//analysis
	if(fwd)
		crct_analysis(image, iw, ih, &rct);
	
	memset(stats, 0, sizeof(stats));
	ac.hi=0xFFFFFFFFFFFF;
	ac.ptr=ac.start=stream;
	ac.end=stream+streamsize;
	if(!fwd)
	{
		ac.code=*(uint64_t*)ac.ptr;//load
		ac.ptr+=sizeof(uint64_t);
		ac.code=ac.code<<32|ac.code>>32;
	}
	{
		int k;
	
		for(k=0;k<3*NPREDS;++k)
			((int32_t*)coeffs)[k]=(1<<SHIFT)/NPREDS;
		bias[0]=1<<NPREDS>>1;
		bias[1]=1<<NPREDS>>1;
		bias[2]=1<<NPREDS>>1;
#ifdef USE_AGETABLE
		for(k=0;k<=AGEMAX;++k)
		{
			int next=k+(k<AGEMAX);
			int sh=31^LZCNT32(k+10);
			agetable[k]=next<<4|sh;
		}
#endif
	}
	perm=perms+3*rct.pidx;
	memset(pixels, 0, psize);
	imptr=image;
	for(int ky=0;ky<ih;++ky)
	{
		Pixel_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky+0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=0;kx<iw;++kx, imptr+=3)
		{
			int offset=0;
			if(fwd)
			{
				yuv[0]=imptr[perm[0]];
				yuv[1]=imptr[perm[1]];
				yuv[2]=imptr[perm[2]];
			}
			for(int kc=0;kc<NCH;++kc)
			{
				int32_t
					NNWWWW	=rows[2][0-4*NCH*NROWS*NVAL],
					NNWWW	=rows[2][0-3*NCH*NROWS*NVAL],
					NNWW	=rows[2][0-2*NCH*NROWS*NVAL],
					NNW	=rows[2][0-1*NCH*NROWS*NVAL],
					NN	=rows[2][0+0*NCH*NROWS*NVAL],
					NNE	=rows[2][0+1*NCH*NROWS*NVAL],
					NNEE	=rows[2][0+2*NCH*NROWS*NVAL],
					NNEEE	=rows[2][0+3*NCH*NROWS*NVAL],
					NNEEEE	=rows[2][0+4*NCH*NROWS*NVAL],
					NNEEEEE	=rows[2][0+5*NCH*NROWS*NVAL],
					NWWWWWW	=rows[1][0-6*NCH*NROWS*NVAL],
					NWWWWW	=rows[1][0-5*NCH*NROWS*NVAL],
					NWWWW	=rows[1][0-4*NCH*NROWS*NVAL],
					NWWW	=rows[1][0-3*NCH*NROWS*NVAL],
					NWW	=rows[1][0-2*NCH*NROWS*NVAL],
					NW	=rows[1][0-1*NCH*NROWS*NVAL],
					N	=rows[1][0+0*NCH*NROWS*NVAL],
					NE	=rows[1][0+1*NCH*NROWS*NVAL],
					NEE	=rows[1][0+2*NCH*NROWS*NVAL],
					NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
					NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
					NEEEEE	=rows[1][0+5*NCH*NROWS*NVAL],
					NEEEEEE	=rows[1][0+6*NCH*NROWS*NVAL],
					WWWWWWW	=rows[0][0-7*NCH*NROWS*NVAL],
					WWWWWW	=rows[0][0-6*NCH*NROWS*NVAL],
					WWWWW	=rows[0][0-5*NCH*NROWS*NVAL],
					WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
					WWW	=rows[0][0-3*NCH*NROWS*NVAL],
					WW	=rows[0][0-2*NCH*NROWS*NVAL],
					W	=rows[0][0-1*NCH*NROWS*NVAL],

					eNEE	=rows[0][1-3*NCH*NROWS*NVAL],
					eNEEE	=rows[0][1-2*NCH*NROWS*NVAL],
					eW	=rows[0][1-1*NCH*NROWS*NVAL];
				int pred, ctx, nbypass, sym, nzeros, tidx, bit, kb;
				Cell_t *statsptr;
				int estims[]=
				{
#define PRED(E) E,
					PREDLIST
#undef  PRED
				};
				int pred0, k, curr;

#if 1
				(void)NNWWWW;
				(void)NNWWW;
				(void)NNWW;
				(void)NNW;
				(void)NN;
				(void)NNE;
				(void)NNEE;
				(void)NNEEE;
				(void)NNEEEE;
				(void)NNEEEEE;
				(void)NWWWWWW;
				(void)NWWWWW;
				(void)NWWWW;
				(void)NWWW;
				(void)NWW;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)NEEEE;
				(void)NEEEEE;
				(void)NEEEEEE;
				(void)WWWWWWW;
				(void)WWWWWW;
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
				(void)eNEE;
				(void)eNEEE;
				(void)eW;
#endif
				pred=bias[kc];
				for(k=0;k<NPREDS;++k)
					pred+=coeffs[kc][k]*(int16_t)estims[k];
				pred>>=SHIFT;
				pred0=pred;

				if(kc==1)offset=yuv[0]*rct.uc0;
				if(kc==2)offset=yuv[0]*rct.vc0+yuv[1]*rct.vc1;
				offset>>=RCTBITS;
				pred+=offset;
				CLAMP(pred, 0, 255);

				ctx=31^LZCNT32(eW*eW+1);
				nbypass=(ctx>>1)-GRBITS;
				CLAMP(nbypass, 0, 6);
				{
					int nzctx=ctx;
					if(nzctx>1)nzctx=1;
					CLAMP(ctx, CTXMIN, CTXMAX);
					ctx-=CTXMIN-nzctx;
				}
				statsptr=stats[kc][ctx][pred];
				
				sym=0;
				nzeros=0;
				if(fwd)
				{
					sym=(int8_t)(yuv[kc]-pred);
					sym=sym<<1^sym>>31;
					nzeros=sym>>nbypass;
				}
				tidx=0;
				_mm_prefetch((char*)statsptr, _MM_HINT_T0);
				for(;;)
				{
					bit=codebit(&ac, statsptr+(tidx<UDEPTH-1?tidx:UDEPTH-1), tidx>=nzeros, fwd);
					if(bit)
						break;
					++tidx;
				}
				for(kb=nbypass-1;kb>=0;--kb)
				{
					bit=codebit(&ac, statsptr+UDEPTH+kb, sym>>kb&1, fwd);
					sym|=bit<<kb;
				}
				if(!fwd)
				{
					sym|=tidx<<nbypass;
					yuv[kc]=(uint8_t)((sym>>1^sym<<31>>31)+pred);
#ifdef ENABLE_GUIDE
					if(yuv[kc]!=g_image[imptr-image+perm[kc]])
					{
						printf("guide  Y %d  X %d  C%d  decoded 0x%02X != original 0x%02X\n"
							, ky
							, kx
							, kc
							, yuv[kc]
							, g_image[imptr-image+perm[kc]]
						);
						CRASH("");
					}
#endif
				}
				curr=yuv[kc]-offset;
#ifdef TRAIN_HIST
				if(fwd)
				{
					int e=curr-pred0;
					CLAMP(e, -THSIZE/2, THSIZE/2);
					++thist[kc][e+THSIZE/2];
				}
#endif
				pred0=(curr>pred0)-(curr<pred0);
				bias[kc]+=pred0<<2;
				for(k=0;k<NPREDS;++k)
					coeffs[kc][k]+=pred0*estims[k];
				rows[0][0+0*NCH*NROWS*NVAL]=curr;
				rows[0][1+0*NCH*NROWS*NVAL]=(2*eW+(sym<<GRBITS)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
			}
			if(!fwd)
			{
				imptr[perm[0]]=yuv[0];
				imptr[perm[1]]=yuv[1];
				imptr[perm[2]]=yuv[2];
			}
		}
	}
	free(pixels);

	//save
	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing\n", dstfn);
		free(image);
		free(stream);
		return 1;
	}
	if(fwd)
	{
		*(uint64_t*)ac.ptr=ac.lo<<32|ac.lo>>32;//flush
		ac.ptr+=sizeof(uint64_t);

		csize=0;
		csize+=fwrite(&tag, 1, 2, fdst);
		csize+=fwrite(&iw, 1, 3, fdst);
		csize+=fwrite(&ih, 1, 3, fdst);
		csize+=fwrite(&rct, 1, sizeof(rct), fdst);
		csize+=fwrite(stream, 1, ac.ptr-stream, fdst);
	}
	else
	{
		int64_t usize2=0;

		usize2+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		usize2+=fwrite(image, 1, usize, fdst);
		usize=usize2;
	}
	fclose(fdst);
	free(image);
	free(stream);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
#ifdef TRAIN_HIST
		for(int k=0;k<THSIZE;++k)
			printf("%12d %12d %12d\n", thist[0][k], thist[1][k], thist[2][k]);
#endif
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
	(void)&time_sec;
	(void)&print_rct;
	return 0;
}
