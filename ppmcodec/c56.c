#ifdef _MSC_VER
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif
#elif defined __linux__
#	ifndef !defined _GNU_SOURCE
#		define _GNU_SOURCE
#	endif
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
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif


#ifdef _MSC_VER
	#define LOUD
//	#define ENABLE_GUIDE
//	#define DEBUG_LZ

//	#define FIFOVAL
#endif


//	#define LZ_SINGLECHANNEL


//runtime
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
#define CLAMP2(X, L, H) X=(X)<(L)?L:X, X=(X)>(H)?H:X
#if defined _M_X64 || defined __x86_64__
#define LZCNT32 _lzcnt_u32
#define LZCNT64 _lzcnt_u64
#define TZCNT32 _tzcnt_u32
#define TZCNT64 _tzcnt_u64
#else
#define LZCNT32 __builtin_clz
#define LZCNT64 __builtin_clzll
#define TZCNT32 __builtin_ctz
#define TZCNT64 __builtin_ctzll
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
static uint8_t *g_image=0;
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

enum
{
//	LBITS=9,
//	EBITS=14,
	LBITS=12,//log2 lookup table size
	EBITS=10,//log2 number of collision cells
	ESIZE=(1<<EBITS),
//	EMASK=ESIZE-1,

	CTRBITS=7,

	LZAC_MIN=5,		//tune
	LZAC_QBITS=6,		//tune
	LZAC_QMAX=(1<<LZAC_QBITS),

	LZAC_STOREBITS=20,
	LZAC_USEBITS=10,	//tune
//	LZAC_SHIFT=LZAC_STOREBITS-LZAC_USEBITS,

//	LZAC_ESTIMBITS=3,
//	LZAC_ESTIMLR=5,		//tune
	
	LZAC_CTXBITS=12,
//	LZAC_CTXBITS=14,
//	LZAC_MODEHIST=8,
//	LZAC_GAMMA_L=28,
//	LZAC_NCTX=8,

//	LZAC_LR_MODE=6,		//tune
//	LZAC_LR_FIXED=5,
//	LZAC_LR_GAMMA_U=6,
//	LZAC_LR_GAMMA_B=7,

	RICEBITS=3,
	RICELIMIT=28,
};
typedef struct _LZAC2Stats
{
	int32_t
		litlen_ustats[RICELIMIT*RICELIMIT], litlen_bstats[(RICELIMIT+1)*512],
		matchlen_ustats[RICELIMIT*RICELIMIT], matchlen_bstats[(RICELIMIT+1)*512],
		backtrack_ustats[RICELIMIT*RICELIMIT], backtrack_bstats[(RICELIMIT+1)*512],
		sym_ustats[RICELIMIT*RICELIMIT], sym_bstats[(RICELIMIT+1)*512];
	//	sym_stats[3][1<<LZAC_CTXBITS][256];
} LZAC2Stats;
static LZAC2Stats stats;
//typedef struct _SymCtx
//{
//	uint8_t sym, ctx;
//} SymCtx;
typedef struct _ETable
{
	int32_t etable[ESIZE], estart, eend, ecount;
} ETable;
static ETable tables[1<<LBITS];
typedef struct _ACState
{
	uint64_t lo, hi, code;
	uint8_t *ptr, *end;
} ACState;
static uint32_t ctrtable[1<<CTRBITS];
INLINE void codebit(ACState *ac, int32_t *pcell, int32_t *pbit, const int fwd)
{
	uint64_t x;
	int32_t cell, ctr, prob, p0;
	int bit;
	
	cell=*pcell;
	x=ac->hi-ac->lo;
	ctr=ctrtable[cell&((1<<CTRBITS)-1)];
	prob=cell>>CTRBITS;
	p0=cell>>(CTRBITS+LZAC_STOREBITS-LZAC_USEBITS);
	p0+=p0<0;
	p0+=1<<LZAC_USEBITS>>1;
	if(x<=0xFFFF)
	{
		if(ac->ptr>=ac->end)
		{
#ifdef _MSC_VER
			CRASH("inflation\n");
#endif
			return;
		}
		if(fwd)
			*(uint32_t*)ac->ptr=(uint32_t)(ac->lo>>32);
		else
			ac->code=ac->code<<32|*(uint32_t*)ac->ptr;
		ac->ptr+=sizeof(uint32_t);
		ac->lo<<=32;
		ac->hi=ac->hi<<32|0xFFFFFFFF;
		if(ac->hi<ac->lo)
			ac->hi=~0ULL;
		x=ac->hi-ac->lo;
	}
	x=ac->lo+(x*p0>>LZAC_USEBITS);
	bit=*pbit;
	bit=fwd?bit:ac->code>=x;
	*pbit=bit;
#ifdef FIFOVAL
	//const int breakidx=45;
	//if(fifoidx==breakidx)//
	//	printf("");
	//if(fifoidx2==breakidx)//
	//	printf("");
	if((uint32_t)(p0-1)>=(1<<LZAC_USEBITS))
		CRASH("");
	if(fwd)
		fifoval_enqueue(cell<<1^bit);
	else
		fifoval_check(cell<<1^bit);
#endif
	*(bit?&ac->lo:&ac->hi)=x+bit-1;
	prob+=((((1-2*bit)<<(LZAC_STOREBITS-1))-prob)>>(uint8_t)ctr);
	*pcell=prob<<CTRBITS|ctr>>16;
}
#if 0
static void codefixed(ACState *ac, int32_t *stats, int *psym, const int nbits, const int fwd)//stats[1<<NBITS]
{
	int sym=0;
	if(fwd)
		sym=*psym;
	for(int kb=nbits-1, tidx=1;kb>=0;--kb)
	{
		int bit=sym>>kb&1;
		codebit(ac, stats+tidx, &bit, fwd, LZAC_LR_FIXED+(kb<1));
		sym|=bit<<kb;
		tidx=2*tidx+bit;
	}
#ifdef _MSC_VER
	if(fwd&&sym!=*psym)
		CRASH("");
#endif
	if(!fwd)
		*psym=sym;
}
static void codegamma(ACState *ac, int32_t *stats, int *psym, const int fwd)//sym>0	stats[2*GAMMA_L+1]
{
	int sym=0, nbits=0, tidx=0, bit;
	if(fwd)
	{
		sym=*psym;
#ifdef _MSC_VER
		if(sym<1)
			CRASH("");
#endif
		nbits=31^LZCNT32(sym);
	}
	do
	{
		bit=tidx>=nbits;
		codebit(ac, stats+tidx, &bit, fwd, LZAC_LR_GAMMA_U);
		++tidx;
	}while(!bit);
	--tidx;
#ifdef _MSC_VER
	if(fwd&&tidx!=nbits)
		CRASH("");
#endif
	nbits=tidx;
	sym|=1<<nbits;
	for(tidx=nbits-1;tidx>=0;--tidx)
	{
		bit=sym>>tidx&1;
		codebit(ac, stats+LZAC_GAMMA_L+1+tidx, &bit, fwd, LZAC_LR_GAMMA_B);
		sym|=bit<<tidx;
	}
	if(!fwd)
		*psym=sym;
}
//static int32_t
//	stats_mode[LZAC_MODEHIST],
//	stats_matchlen[2*LZAC_GAMMA_L+1],
//	stats_backtrack[2*LZAC_GAMMA_L+1],
//	stats_qlen[LZAC_QMAX],
//	stats_sym[3*LZAC_NCTX][256];
#endif
static int codefixed(ACState *ac, int32_t *stats, int sym, const int nbits, const int fwd)//stats[1<<NBITS]
{
	for(int kb=nbits-1, tidx=1;kb>=0;--kb)
	{
		int bit=sym>>kb&1;
		codebit(ac, stats+tidx, &bit, fwd);
		sym|=bit<<kb;
		tidx=2*tidx+bit;
	}
#ifdef FIFOVAL
	if(fwd)
		fifoval_enqueue((uint32_t)(size_t)stats^sym);
	else
		fifoval_check((uint32_t)(size_t)stats^sym);
#endif
	return sym;
}
static int coderice(ACState *ac, int32_t *ustats, int32_t *bstats, int sym, int *pestim, int ricelimit, int maxbits, const int fwd)
{
	int ricek, q, tidx, bit, limitflag, ricek2, tidx2, nbits;

	//ustats[RICELIMIT_SYM*RICELIMIT_SYM], bstats[(RICELIMIT_SYM+1)*MIN(NLEVELS, 512)]
	ricek=31-LZCNT32((*pestim>>RICEBITS)+1);
	q=sym>>ricek;
	tidx=0;
	ustats+=ricek*ricelimit;
	do
	{
		bit=tidx>=q;
		codebit(ac, ustats+tidx, &bit, fwd);
		if(bit)
			break;
		++tidx;
	}
	while(tidx<ricelimit);
	q=tidx<<ricek;
	limitflag=0;
	nbits=maxbits;
	if(maxbits>9)
		nbits=9;
	bstats+=(ptrdiff_t)ricek<<nbits;//ctx
	ricek2=ricek;
	if(tidx>=ricelimit)
	{
		limitflag=1;
		sym-=q;
		ricek2=(1<<maxbits)-q;
		if(ricek2<=0)
			CRASH("ricek2<=0");
		ricek2=32-LZCNT32(ricek2);
		tidx=0;
	}
	tidx+=1<<maxbits>>ricek2;//tidx=0b1xx...
	tidx2=1;
	if(tidx<=511)
		tidx2=tidx;
	for(int kb=ricek2-1;kb>=0;--kb)
	{
		bit=sym>>kb&1;
		codebit(ac, bstats+tidx2, &bit, fwd);
		tidx=2*tidx+bit;
		if(tidx<=511)
			tidx2=tidx;
	}
#ifdef _MSC_VER
	if(fwd&&tidx-(1<<maxbits)!=sym)
		CRASH("enc rice");
#endif
	if(limitflag)
		tidx+=q;
	sym=tidx-(1<<maxbits);
	*pestim+=((sym<<RICEBITS)-*pestim)>>(RICEBITS+1);
#ifdef FIFOVAL
	if(fwd)
		fifoval_enqueue(sym);
	else
		fifoval_check(sym);
#endif
	return sym;
}
static void findmatch(uint8_t *imptr, uint8_t *imstart, uint8_t *searchend, uint8_t *imend, int *ret_matchlen, int *ret_matchidx)
{
	uint64_t x;
	ETable *table;
	int matchlen=0, matchidx=0, ctr;

	x=*(uint64_t*)imptr+1;
	x^=x>>33;
	x*=0xFF51AFD7ED558CCDULL;
	x^=x>>33;
//	x*=0xC4CEB9FE1A85EC53ULL;
//	x^=x>>33;
	table=tables+(x>>(64-LBITS));
	ctr=table->ecount;
	while(ctr)
	{
		int tidx=(table->estart+ctr-1)&((1<<EBITS)-1);
		int idx=table->etable[tidx], len=0;
		uint8_t *search1=imstart+idx;
		uint8_t *search2=imptr;
		uint64_t cmp1=~0ULL, cmp2=~0ULL;
		uint64_t c1, c2;

		if(*(uint32_t*)search1==*(uint32_t*)search2)
		{
			for(;;)
			{
				cmp1=((uint64_t*)search1)[0]^((uint64_t*)search2)[0];//22%
				cmp2=((uint64_t*)search1)[1]^((uint64_t*)search2)[1];
				if((cmp1|cmp2)||search2>=searchend)
					break;
				search1+=sizeof(uint64_t[2]);
				search2+=sizeof(uint64_t[2]);
			}
			c1=TZCNT64(cmp1)>>3;
			c2=TZCNT64(cmp2)>>3;
			if(!cmp1)
				len+=(int)c2;
			len+=(int)(search2-imptr+c1);
			if(len>imend-imptr)
				len=(int)(imend-imptr);
			if(len>matchlen||(len==matchlen&&idx>matchidx))
				matchidx=idx, matchlen=len;
		}
		--ctr;
	}
	table->etable[table->eend]=(int32_t)(imptr-imstart);
	table->eend=(table->eend+1)&((1<<EBITS)-1);
	if(table->ecount>=1<<EBITS)
		table->estart=(table->estart+1)&((1<<EBITS)-1);
	else
		++table->ecount;
	*ret_matchlen=matchlen;
	*ret_matchidx=matchidx;
}
typedef struct _LZRectangle
{
	int x, y, dx, dy;
} LZRectangle;
typedef struct _Instruction
{
	LZRectangle;
	int x0, y0;
} Instruction;
static void findredundancy2D(uint8_t *image, int iw, int ih, int kx, int ky, LZRectangle *ret)
{
	uint64_t x;
	ETable *table;
//	int matchlen=0, matchidx=0;
	LZRectangle rect={0};
	int bestsum=1, bestcount=0;
	int ctr=0;
	int rowstride=3*iw;
	uint8_t *imptr=image+rowstride*ky+3*kx;
	uint8_t *rowend=image+rowstride*(ky+1);
	uint8_t *imend=image+3*iw*ih;

	x=*(uint64_t*)imptr+1;
	x^=x>>33;
	x*=0xFF51AFD7ED558CCDULL;
	x^=x>>33;
//	x*=0xC4CEB9FE1A85EC53ULL;
//	x^=x>>33;
	table=tables+(x>>(64-LBITS));
	ctr=table->ecount;
	while(ctr)
	{
		int tidx=(table->estart+ctr-1)&((1<<EBITS)-1);
		int idx=table->etable[tidx], len=0;

#if 1

		if((*(uint32_t*)(image+3*idx)&0xFFFFFF)==(*(uint32_t*)(image+3*(iw*ky+kx))&0xFFFFFF))
		{
			int sy=idx/iw, sx=idx%iw, dx=1, dy=1;
			uint8_t *p1=image+3*(iw*sy+sx);
			uint8_t *p2=image+3*(iw*ky+kx);
			int sum=1, count=0, stop=0;

		again:
			while(stop!=3)
			{
				if(!(stop&1))
				{
					if(kx+dx>=iw||sx+dx>=iw)
						stop|=1;
					else
					{
						for(int y2=0, yidx=3*dx;y2<dy;++y2, yidx+=rowstride)//expand right
						{
							int sum2=sum
								+abs(p2[yidx+0]-p1[yidx+0])
								+abs(p2[yidx+1]-p1[yidx+1])
								+abs(p2[yidx+2]-p1[yidx+2])
							;
							int count2=count+3;
							if(sum2>count2)
							{
								stop|=1;
								break;
							}
							sum=sum2;
							count=count2;
						}
						if(!(stop&1))
							++dx;
					}
				}
				if(!(stop&2))
				{
					if(ky+dy>=ih)
						stop|=2;
					else
					{
						for(int x2=0, xidx=dy*rowstride;x2<dx;++x2, xidx+=3)//expand down
						{
							int sum2=sum
								+abs(p2[xidx+0]-p1[xidx+0])
								+abs(p2[xidx+1]-p1[xidx+1])
								+abs(p2[xidx+2]-p1[xidx+2])
							;
							int count2=count+3;
							if(sum2>count2)
							{
								stop|=2;
								break;
							}
							sum=sum2;
							count=count2;
						}
						if(!(stop&2))
							++dy;
					}
				}
			}
#if 1
			if(sx+dx>iw||sy+dy>ih||kx+dx>iw||ky+dy>ih)//
			{
				stop=0;
				dx=1;
				dy=1;
				goto again;
			}
#endif
			if((int64_t)bestsum*count>(int64_t)sum*bestcount)
			{
				bestcount=count;
				bestsum=sum;
				rect.x=sx;
				rect.y=sy;
				rect.dx=dx;
				rect.dy=dy;
			}
		}
#endif
#if 0
		int sy=idx/iw, sx=idx%iw;
		int sum=0, count=0, xmin=0;

		xmin=kx;
		if(xmin<sx)
			xmin=sx;
		xmin=iw-xmin;
		for(;;)
		{
			uint8_t *p1=image+3*(iw*sy+sx);
			uint8_t *p2=image+3*(iw*ky+kx);
			uint8_t *end2=image+3*(iw*ky+xmin);
			int width=0;
			int sum0=sum, count0=count;

			for(;;)
			{
				int sum2=sum+abs(p2[0]-p1[0])+abs(p2[1]-p1[1])+abs(p2[2]-p1[2]);
				int count2=count+3;
				if(p2>=end2||sum2>=count2||width>xmin)
				{
					if(width<xmin)
					{
					}
					break;
				}
				sum=sum2;
				count=count2;
				p1+=3;
				p2+=3;
				++width;
			}
		}
#endif
#if 0
		uint8_t *search1=image+idx;
		uint8_t *search2=imptr;
		uint64_t cmp1=~0ULL, cmp2=~0ULL;
		uint64_t c1, c2;

		if(*(uint32_t*)search1==*(uint32_t*)search2)
		{
#if 0
			int sum=0, count=0;

			for(;;)
			{
				uint8_t *start1=search1, *start2=search2;

				if(start2+rowstride>imend||sum>=count)
					break;
				for(;;)
				{
					int sum2=sum+abs(search2[0]-search1[0])+abs(search2[1]-search1[1])+abs(search2[2]-search1[2]);
					int count2=count+3;
					if(search2>=rowend||sum2>=count2)
						break;
					sum=sum2;
					count=count2;
					search1+=3;
					search2+=3;
				}
				search1=start1+rowstride;
				search2=start2+rowstride;
			}
#endif
#if 0
			for(;;)
			{
				cmp1=((uint64_t*)search1)[0]^((uint64_t*)search2)[0];//22%
				cmp2=((uint64_t*)search1)[1]^((uint64_t*)search2)[1];
				if((cmp1|cmp2)||search2>=searchend)
					break;
				search1+=sizeof(uint64_t[2]);
				search2+=sizeof(uint64_t[2]);
			}
			c1=TZCNT64(cmp1)>>3;
			c2=TZCNT64(cmp2)>>3;
			if(!cmp1)
				len+=(int)c2;
			len+=(int)(search2-imptr+c1);
			if(len>imend-imptr)
				len=(int)(imend-imptr);
			if(len>matchlen||(len==matchlen&&idx>matchidx))
				matchidx=idx, matchlen=len;
#endif
		}
#endif
		--ctr;
	}
	table->etable[table->eend]=(int32_t)((imptr-image)/3);
	table->eend=(table->eend+1)&((1<<EBITS)-1);
	if(table->ecount>=1<<EBITS)
		table->estart=(table->estart+1)&((1<<EBITS)-1);
	else
		++table->ecount;
	memcpy(ret, &rect, sizeof(*ret));
	//*ret_matchlen=matchlen;
	//*ret_matchidx=matchidx;
}
#if 0
static void memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
{
	size_t copied;
	char *d=(char*)dst;
	const char *s=(const char*)src;

	if(dstbytes<srcbytes)
	{
		memcpy(dst, src, dstbytes);
		return;
	}
	copied=srcbytes;
	memcpy(d, s, srcbytes);
	while((copied<<1)<=dstbytes)
	{
		memcpy(d+copied, d, copied);
		copied<<=1;
	}
	if(copied<dstbytes)
		memcpy(d+copied, d, dstbytes-copied);
}
#endif
static void ppm_save(const char *fn, const uint8_t *image, int iw, int ih)
{
	FILE *fdst=fopen(fn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", fn);
		return;
	}
	fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	fwrite(image, 1, (ptrdiff_t)3*iw*ih, fdst);
	fclose(fdst);
}
int c56_codec(int argc, char **argv)
{
	const int16_t tag='5'|'6'<<8;

	const char *srcfn=0, *dstfn=0;
	int fwd=0, iw=0, ih=0;
	int rowstride=0;
//	int estim[]={32, 32, 32};
//	int prevmode=1, prevqlen=0, prevmatchlen=0, prevbacktrack=0;
//#ifdef LZ_SINGLECHANNEL
//	const
//#endif
	int kc=0;
	ptrdiff_t usize=0, csize=0, headersize=0, cap=0;
	uint8_t *buf=0, *image=0, *streamstart=0, *streamend=0;
	ACState ac={0};
	uint8_t *imptr=0, *imend=0;
	int32_t litlen_estim=1<<RICEBITS, matchlen_estim=1<<RICEBITS, backtrack_estim=1<<RICEBITS;
	int32_t sym_estim[3]={1<<RICEBITS, 1<<RICEBITS, 1<<RICEBITS};
#ifdef ENABLE_GUIDE
	static uint8_t *im0=0;
#endif
#ifdef LOUD
	double t=0, t2=0;
#endif
#ifdef DEBUG_LZ
	int debug_target=
		0
	//	71868
	//	71878
	//	138173
	//	6382425
	;
	int debug_block=0;
#endif

	//(void)prevqlen;
	//(void)prevmatchlen;
	//(void)prevbacktrack;
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
	srcfn=argv[1];
	dstfn=argv[2];
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
		if(!fwd&&c!=tag)
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
		rowstride=3*iw;
		usize=(ptrdiff_t)3*iw*ih;
		cap=(ptrdiff_t)7*iw*ih;
		buf=(uint8_t*)malloc(cap);
		if(!buf)
			CRASH("Alloc error");
#ifdef _MSC_VER
		memset(buf, 0, cap);
#endif
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

			memset(image-(ptrdiff_t)3*rowstride, 0, (ptrdiff_t)3*rowstride);
			fread(image, 1, usize, fsrc);//read image
#ifdef ENABLE_GUIDE
			g_image=guide_save(image, iw, ih);
#endif
		}
		else
		{
			struct stat info={0};
			stat(srcfn, &info);
			csize=info.st_size;

			image=buf+(ptrdiff_t)3*rowstride;
			streamstart=buf+cap-csize-sizeof(uint64_t);
			streamend=buf+cap;
			memset(image-(ptrdiff_t)3*rowstride, 0, (ptrdiff_t)3*rowstride);

			fread(streamstart, 1, csize-headersize, fsrc);//read stream
		}
		fclose(fsrc);
	}
	//memset(stats_mode, 0, sizeof(stats_mode));
	//memset(stats_matchlen, 0, sizeof(stats_matchlen));
	//memset(stats_backtrack, 0, sizeof(stats_backtrack));
	//memset(stats_qlen, 0, sizeof(stats_qlen));
	//memset(stats_sym, 0, sizeof(stats_sym));
	memset(&stats, 0, sizeof(stats));
	ac.hi=0xFFFFFFFFFFFF;
	ac.ptr=streamstart;
	ac.end=streamend;
	imptr=image;
	imend=image+usize;
	for(int k=0;k<1<<CTRBITS;++k)
	{
		int next=k+1, sh=31^LZCNT32(k+2);
		if(next>(1<<CTRBITS)-1)
			next=(1<<CTRBITS)-1;
		ctrtable[k]=next<<16|sh;
	}
	if(fwd)
	{
#if 1
		uint8_t *imptr0=imptr, *searchend=image+usize-sizeof(uint64_t[2]);
		
#if 0
		{
			int saved=0;
			int64_t res=(int64_t)iw*ih;
			int cap=(int)res, rcount=0;
			Instruction *rects=(Instruction*)malloc(cap*sizeof(Instruction));
			uint8_t *mask=(uint8_t*)malloc(res);
			double t=time_sec();
			if(!rects||!mask)
				CRASH("Alloc error");
			memset(mask, 0, res);
			for(int ky=0;ky<ih;++ky)
			{
				for(int kx=0;kx<iw;++kx)
				{
					if(!mask[iw*ky+kx])
					{
						LZRectangle rect={0};

						findredundancy2D(image, iw, ih, kx, ky, &rect);

						if(rect.dx>=2&&rect.dy>=1&&rcount<cap)
						{
							for(int ky2=0;ky2<rect.dy;++ky2)
							{
								for(int kx2=0;kx2<rect.dx;++kx2)
								{
									LZRectangle rect2={0};
									findredundancy2D(image, iw, ih, kx, ky, &rect2);
								}
							}
							Instruction *instr=rects+rcount++;
							memcpy(instr, &rect, sizeof(rect));
							instr->x0=kx;
							instr->y0=ky;
							for(int ky2=0;ky2<rect.dy;++ky2)
							{
								for(int kx2=0;kx2<rect.dx;++kx2)
									mask[iw*(ky+ky2)+kx+kx2]=1;
							}
#if 1
							if(instr->x0+instr->dx>iw||instr->x+instr->dx>iw||instr->y0+instr->dy>ih||instr->y+instr->dy>ih)//
								CRASH("");
#endif
						}
					}
				}
				printf("%5d  %12.6lf mins\r", ky+1, (time_sec()-t)*(ih-(ky+1))/(60*(ky+1)));
			}
			printf("\n");
			for(int k=0;k<res;++k)
				saved+=mask[k];
			printf("Saved %d/%lld bytes  %8.4lf%%\n", 3*saved, 3*res, 100.*(double)(res-saved)/res);
			for(int k=rcount-1;k>=0;--k)
			{
				Instruction *instr=rects+k;

				for(int ky=instr->dy-1;ky>=0;--ky)
				{
					for(int kx=instr->dx-1;kx>=0;--kx)
					{
						int idx1=3*(iw*(instr->y0+ky)+instr->x0+kx);
						int idx2=3*(iw*(instr->y+ky)+instr->x+kx);
						//image[idx2+0]=0xFF;
						//image[idx2+1]=0x00;
						//image[idx2+2]=0xFF;
						image[idx2+0]-=image[idx1+0];
						image[idx2+1]-=image[idx1+1];
						image[idx2+2]-=image[idx1+2];
					}
				}
			}
			ppm_save("20260630Tu_0705pm.ppm", image, iw, ih);
			//for(int ky=0;ky<ih;++ky)
			//{
			//	for(int kx=0;kx<iw;++kx)
			//	{
			//	}
			//}
			exit(0);
		}
#endif
		//[lit_len/Rice][literals/fixed8...][match_len/Rice][offset/Rice]
		memset(tables, 0, sizeof(tables));
		for(;;)
		{
			int matchlen=0, matchidx=0;

			if(imptr<imend)
			{
				findmatch(imptr, image, searchend, imend, &matchlen, &matchidx);
				if(matchlen<=LZAC_MIN)
				{
					++imptr;
					continue;
				}
				{
					uint8_t *imptr1=imptr;
					int matchlen1=0, matchidx1=0;
					int matchlen2=0, matchidx2=0;
					int matchlen3=0, matchidx3=0;
					int matchlen4=0, matchidx4=0;

					findmatch(imptr+1, image, searchend, imend, &matchlen1, &matchidx1);
					findmatch(imptr+2, image, searchend, imend, &matchlen2, &matchidx2);
					findmatch(imptr+3, image, searchend, imend, &matchlen3, &matchidx3);
					findmatch(imptr+4, image, searchend, imend, &matchlen4, &matchidx4);
					if(matchlen1>matchlen+LZAC_MIN)
					{
						imptr=imptr1+1;
						matchlen=matchlen1;
						matchidx=matchidx1;
					}
					if(matchlen2>matchlen+LZAC_MIN)
					{
						imptr=imptr1+2;
						matchlen=matchlen2;
						matchidx=matchidx2;
					}
					if(matchlen3>matchlen+LZAC_MIN)
					{
						imptr=imptr1+3;
						matchlen=matchlen3;
						matchidx=matchidx3;
					}
					if(matchlen4>matchlen+LZAC_MIN)
					{
						imptr=imptr1+4;
						matchlen=matchlen4;
						matchidx=matchidx4;
					}
				}
			}
			coderice(&ac, stats.litlen_ustats, stats.litlen_bstats, (int)(imptr-imptr0), &litlen_estim, RICELIMIT, RICELIMIT, 1);
			kc=(uint32_t)(imptr0-image)%3;
			while(imptr0<imptr)
			{
				int NW, N, W, pred, sym;

				NW=imptr0[-rowstride-3];
				N=imptr0[-rowstride];
				W=imptr0[-3];
				if(kc)
				{
					NW-=imptr0[-rowstride-3-1];
					N-=imptr0[-rowstride-1];
					W-=imptr0[-3-1];
				}
				pred=abs(N-NW)>abs(W-NW)?N:W;
				if(kc)
				{
					pred+=imptr0[-1];
					CLAMP2(pred, 0, 255);
				}
				sym=(int8_t)(*imptr0-pred);
				sym=sym<<1^sym>>31;
				coderice(&ac, stats.sym_ustats, stats.sym_bstats, sym, sym_estim+kc, RICELIMIT, 8, 1);
				++imptr0;
				++kc;
				if(kc>2)
					kc=0;
			}
			if(imptr>=imend)
				break;
			coderice(&ac, stats.matchlen_ustats, stats.matchlen_bstats, matchlen, &matchlen_estim, RICELIMIT, RICELIMIT, 1);
			coderice(&ac, stats.backtrack_ustats, stats.backtrack_bstats, (int)(imptr-image-matchidx), &backtrack_estim, RICELIMIT, RICELIMIT, 1);
			imptr+=matchlen;
			imptr0=imptr;
			if(imptr>=imend)
				break;
		}
		*(uint64_t*)ac.ptr=ac.lo<<32|ac.lo>>32;//flush
		ac.ptr+=sizeof(uint64_t);
		streamend=ac.ptr;
#endif
#if 0
		const int queuesize=sizeof(SymCtx[LZAC_QMAX]);
		SymCtx *queue=(SymCtx*)malloc(queuesize);
		int qcount=0, qstart=0, qend=0;
		uint8_t *matchend=image+usize-sizeof(uint64_t[4]), *matchend2=matchend-(sizeof(uint64_t[4])+1);
		if(!queue)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(tables, 0, sizeof(tables));
		memset(queue, 0, queuesize);
		while(imptr<imend)
		{
			int matchidx=-1, matchlen=0;
			uint8_t pixel=*imptr;
			//int lookup=pixel;
			//if(imptr>=image+rowstride)
			//	lookup=(lookup+imptr[-rowstride])>>1;
			ETable *table=tables+(((*(uint16_t*)imptr+1)*0x9E3779B9)>>(32-LBITS));
		//	ETable *table=tables+pixel;
			int ctr=table->ecount;
			while(ctr)
			{
				int tidx=(table->estart+ctr)&EMASK;
				int idx=table->etable[tidx];
				uint8_t *match2=imptr;

#if 0
				int64_t len=0;
				ptrdiff_t offset1=image+idx-match2;
				__m256i ones=_mm256_set1_epi8(-1);
				__m256i cand=ones;
				if(match2<matchend)
				{
					for(;;)
					{
						__m256i truth=_mm256_loadu_si256((__m256i*)match2);
						cand=_mm256_loadu_si256((__m256i*)(match2+offset1));
						cand=_mm256_cmpeq_epi8(cand, truth);
						cand=_mm256_xor_si256(cand, ones);
						int mask=_mm256_movemask_epi8(cand);
						int64_t cond=(matchend2-match2)>>63;
						cond=mask|cond>>63;
						if(cond)
							break;
						match2+=sizeof(__m256i);
					}
				}
				int64_t cmp0=_mm256_extract_epi64(cand, 0);
				int64_t cmp1=_mm256_extract_epi64(cand, 1);
				int64_t cmp2=_mm256_extract_epi64(cand, 2);
				int64_t cmp3=_mm256_extract_epi64(cand, 3);
				cmp0=_tzcnt_u64(cmp0);
				cmp1=_tzcnt_u64(cmp1);
				cmp2=_tzcnt_u64(cmp2);
				cmp3=_tzcnt_u64(cmp3);
				cmp0>>=3;
				cmp1>>=3;
				cmp2>>=3;
				cmp3>>=3;
				cmp1=cmp0<8?0:cmp1;
				len+=(int)cmp0;
				cmp2=cmp1<8?0:cmp2;
				len+=(int)cmp1;
				cmp3=cmp2<8?0:cmp3;
				len+=(int)cmp2;
				len+=(int)cmp3;
				len+=(int)(match2-imptr);
#ifdef _MSC_VER
				if(imptr+len>imend)//
					CRASH("");
#endif
#endif
#if 1
				int64_t len=0;
				ptrdiff_t offset1=image+idx-match2;
				uint64_t cmp0=~0ULL;
				uint64_t cmp1=~0ULL;
				uint64_t cmp2=~0ULL;
				uint64_t cmp3=~0ULL;
				//_mm_prefetch((char*)(match2+offset1), _MM_HINT_T0);
				if(match2<matchend)
				{
					for(;;)
					{
						int64_t cond=(matchend2-match2)>>63;
						cmp0=*(uint64_t*)(match2+offset1+0*sizeof(uint64_t))^*(uint64_t*)(match2+0*sizeof(uint64_t));
						cmp1=*(uint64_t*)(match2+offset1+1*sizeof(uint64_t))^*(uint64_t*)(match2+1*sizeof(uint64_t));
						cmp2=*(uint64_t*)(match2+offset1+2*sizeof(uint64_t))^*(uint64_t*)(match2+2*sizeof(uint64_t));
						cmp3=*(uint64_t*)(match2+offset1+3*sizeof(uint64_t))^*(uint64_t*)(match2+3*sizeof(uint64_t));
						cond=cmp0|cmp1|cmp2|cmp3|cond>>63;
						if(cond)
							break;
						match2+=sizeof(uint64_t[4]);
					}
				}
				cmp0=_tzcnt_u64(cmp0);
				cmp1=_tzcnt_u64(cmp1);
				cmp2=_tzcnt_u64(cmp2);
				cmp3=_tzcnt_u64(cmp3);
				cmp0>>=3;
				cmp1>>=3;
				cmp2>>=3;
				cmp3>>=3;
				cmp1=cmp0<8?0:cmp1;
				len+=(int)cmp0;
				cmp2=cmp1<8?0:cmp2;
				len+=(int)cmp1;
				cmp3=cmp2<8?0:cmp3;
				len+=(int)cmp2;
				len+=(int)cmp3;
				len+=(int)(match2-imptr);
#endif

				//int len=0;
				//ptrdiff_t offset1=image+idx-match2;
				//uint64_t cmp=~0ULL;
				//if(match2<matchend)
				//{
				//	for(;;)
				//	{
				//		int64_t cond=(matchend2-match2)>>63;
				//		cmp=*(uint64_t*)(match2+offset1)^*(uint64_t*)match2;
				//		cond=cmp|cond>>63;
				//		if(cond)
				//			break;
				//		match2+=sizeof(uint64_t);
				//	}
				//}
				//len=(int)(match2-imptr+(_tzcnt_u64(cmp)>>3));
				
				//int len=0;
				//ptrdiff_t searchend=image+usize-sizeof(uint64_t);
				//uint8_t *search1=image+idx;
				//while(len<searchend&&*(uint64_t*)(search1+len)==*(uint64_t*)(match2+len))
				//	len+=sizeof(uint64_t);
				//searchend=image+usize-imptr-1;
				//while(len<searchend&&search1[len]==match2[len])
				//	++len;


				if(len>matchlen)
					matchidx=idx, matchlen=(int)len;

				//check for other previous encounters
				--ctr;
			}
			//add current position to table
			if(table->ecount>=ESIZE)//table is full
			{
				table->etable[table->eend]=(int32_t)(imptr-image);
				table->eend=(table->eend+1)&EMASK;
				table->estart=(table->estart+1)&EMASK;
			}
			else
			{
				table->etable[table->eend]=(int32_t)(imptr-image);
				table->eend=(table->eend+1)&EMASK;
				++table->ecount;
			}

			if(qcount>=LZAC_QMAX)//queue is full
			{
#ifdef DEBUG_LZ
				++debug_block;
				if(debug_block==debug_target)//
					__debugbreak();
#endif
				{
					int mode=1;
					codebit(&ac, stats_mode+(prevmode&(LZAC_MODEHIST-1)), &mode, fwd, LZAC_LR_MODE);
					prevmode=prevmode<<1|mode;
				}
				{
					int c1=qcount-1;
					codefixed(&ac, stats_qlen, &c1, LZAC_QBITS, fwd);
					prevqlen=c1;
				}
				while(qcount)
				{
					SymCtx *p=queue+qstart;
					int s=p->sym;
					codefixed(&ac, stats_sym[p->ctx], &s, 8, fwd);
					qstart=(qstart+1)%LZAC_QMAX;
					--qcount;
				}
				qstart=qend;
				qcount=0;
			}
			if(matchlen>=LZAC_MIN)//emit match
			{
				if(qcount)
				{
#ifdef DEBUG_LZ
					++debug_block;
					if(debug_block==debug_target)//
						__debugbreak();
#endif
					{
						int mode=1;
						codebit(&ac, stats_mode+(prevmode&(LZAC_MODEHIST-1)), &mode, fwd, LZAC_LR_MODE);
						prevmode=prevmode<<1|mode;
					}
					{
						int c1=qcount-1;//[0 ~ LZAC_QMAX-1]
						codefixed(&ac, stats_qlen, &c1, LZAC_QBITS, fwd);
						prevqlen=c1;
					}
					while(qcount)
					{
						SymCtx *p=queue+qstart;
						int s=p->sym;
						codefixed(&ac, stats_sym[p->ctx], &s, 8, fwd);
						qstart=(qstart+1)%LZAC_QMAX;
						--qcount;
					}
				}
#ifdef DEBUG_LZ
				++debug_block;
				if(debug_block==debug_target)//
					__debugbreak();
#endif
				{
					int mode=0;
					codebit(&ac, stats_mode+(prevmode&(LZAC_MODEHIST-1)), &mode, fwd, LZAC_LR_MODE);
					prevmode=prevmode<<1|mode;
				}
				codegamma(&ac, stats_matchlen, &matchlen, fwd);
				prevmatchlen=matchlen;
				{
					int backtrack=(int)(imptr-image-matchidx);
					codegamma(&ac, stats_backtrack, &backtrack, fwd);
					prevbacktrack=backtrack;
				}

				//jump
				imptr+=matchlen;
#ifndef LZ_SINGLECHANNEL
				kc+=matchlen;
				kc%=3;
#endif
			}
			else//encode pixel normally
			{
				int ctx=31^LZCNT32(estim[kc]*estim[kc]+1);
				if(ctx>LZAC_NCTX-1)
					ctx=LZAC_NCTX-1;
				int
					NW	=imptr[-rowstride-3],
					N	=imptr[          -3],
					W	=imptr[-rowstride  ];
				if(kc)//cRCT
				{
					NW	-=imptr[-rowstride-3-1];
					N	-=imptr[          -3-1];
					W	-=imptr[-rowstride  -1];
				}

				int pred=abs(N-NW)>abs(W-NW)?N:W;
				//int pred=N+W-NW, vmax=N, vmin=W;
				//if(N<W)vmin=N, vmax=W;
				//CLAMP2(pred, vmin, vmax);

				if(kc)
				{
					pred+=imptr[-1];
					CLAMP2(pred, 0, 255);
				}
				int sym=(int8_t)(pixel-pred);
				sym=sym<<1^sym>>31;
				SymCtx *p=queue+qend;
				p->sym=sym;//enqueue symbol
				p->ctx=3*ctx+kc;
				estim[kc]+=((sym<<LZAC_ESTIMBITS)-estim[kc]+(1<<LZAC_ESTIMLR>>1))>>LZAC_ESTIMLR;
				qend=(qend+1)%LZAC_QMAX;
				++qcount;

				++imptr;
#ifndef LZ_SINGLECHANNEL
				++kc;
				if(kc>=3)
					kc=0;
#endif
			}
		}
		if(qcount)
		{
#ifdef DEBUG_LZ
			++debug_block;
			if(debug_block==debug_target)//
				__debugbreak();
#endif
			{
				int mode=1;
				codebit(&ac, stats_mode+(prevmode&(LZAC_MODEHIST-1)), &mode, fwd, LZAC_LR_MODE);
				prevmode=prevmode<<1|mode;
			}
			{
				int c1=qcount-1;
				codefixed(&ac, stats_qlen, &c1, LZAC_QBITS, fwd);
				prevqlen=c1;
			}
			while(qcount)
			{
				SymCtx *p=queue+qstart;
				int s=p->sym;
				codefixed(&ac, stats_sym[p->ctx], &s, 8, fwd);
				qstart=(qstart+1)%LZAC_QMAX;
				--qcount;
#ifndef LZ_SINGLECHANNEL
				++kc;
				if(kc>=3)
					kc=0;
#endif
			}
		}
		free(queue);
		ac.low=ac.low<<32|ac.low>>32;
		*(uint64_t*)ac.ptr=ac.low;//flush
		ac.ptr+=sizeof(uint64_t);
		streamend=ac.ptr;
#endif
	}
	else//dec
	{
		int litlen, matchlen, bt;
		uint8_t *srcptr, *matchend;

		ac.code=*(uint64_t*)ac.ptr;//load
		ac.ptr+=sizeof(uint64_t);
		ac.code=ac.code<<32|ac.code>>32;
		while(imptr<imend)
		{
			litlen=coderice(&ac, stats.litlen_ustats, stats.litlen_bstats, 0, &litlen_estim, RICELIMIT, RICELIMIT, 0);
			kc=(uint32_t)(imptr-image)%3;
			while(litlen--)
			{
				int NW, N, W, pred, sym;

				NW=imptr[-rowstride-3];
				N=imptr[-rowstride];
				W=imptr[-3];
				if(kc)
				{
					NW-=imptr[-rowstride-3-1];
					N-=imptr[-rowstride-1];
					W-=imptr[-3-1];
				}
				pred=abs(N-NW)>abs(W-NW)?N:W;
				if(kc)
				{
					pred+=imptr[-1];
					CLAMP2(pred, 0, 255);
				}
				sym=coderice(&ac, stats.sym_ustats, stats.sym_bstats, 0, sym_estim+kc, RICELIMIT, 8, 0);
				sym=sym>>1^sym<<31>>31;
				sym=(uint8_t)(sym+pred);
				*imptr=sym;
#ifdef ENABLE_GUIDE
				if(*imptr!=g_image[imptr-image])
					CRASH("guide  %d", imptr-image);
#endif
				++imptr;
				++kc;
				if(kc>2)
					kc=0;
			}
			if(imptr>=imend)
				break;
			matchlen=coderice(&ac, stats.matchlen_ustats, stats.matchlen_bstats, 0, &matchlen_estim, RICELIMIT, RICELIMIT, 0);
			bt=coderice(&ac, stats.backtrack_ustats, stats.backtrack_bstats, 0, &backtrack_estim, RICELIMIT, RICELIMIT, 0);
			srcptr=imptr-bt;
			matchend=imptr+matchlen;

#if 1
			if(bt<matchlen)//overlap
			{
				size_t copied=bt;

				while((size_t)(copied-1)<((size_t)matchlen+bt-1))
				{
					memcpy(imptr, srcptr, copied);
					imptr+=copied;
					copied+=copied;
					if(copied>(size_t)(matchend-imptr))
						copied=matchend-imptr;
				}
				//if(imptr!=matchend)
				//	CRASH("");
			}
			else//disjoint
				memcpy(imptr, srcptr, matchlen);
			imptr=matchend;
#endif
			//while(imptr<matchend)
			//{
			//	*imptr=imptr[-bt];
			//	++imptr;
			//}
#if 0
			if(bt==1)
			{
				memset(imptr, imptr[-1], matchlen);
				imptr+=matchlen;
			}
			else
			{
				while(imptr<matchend)
					*imptr++=*srcptr++;
			}
#endif

			//X
#if 0
			//if(bt==1)
			//	printf("");
			if(matchlen>bt)
			{
				uint8_t *fillend=imptr+bt, *s0=srcptr;

				while(imptr<fillend)
					*imptr++=*srcptr++;
				srcptr=s0;
			}
			while(imptr<matchend)
			{
				*(uint64_t*)imptr=*(uint64_t*)srcptr;
#ifdef ENABLE_GUIDE
				if(*imptr!=g_image[imptr-image])
					CRASH("guide  %d", imptr-image);
#endif
				srcptr+=sizeof(uint64_t);
				imptr+=sizeof(uint64_t);
			}
			imptr=matchend;
#endif
			//while(imptr<matchend)
			//	*imptr++=*srcptr++;
		}
#if 0
		ac.code=*(uint64_t*)ac.ptr;//load
		ac.ptr+=sizeof(uint64_t);
		ac.code=ac.code<<32|ac.code>>32;
		while(imptr<imend)
		{
			uint8_t *fillend;
			int mode, count;

#ifdef DEBUG_LZ
			++debug_block;
			if(debug_block==debug_target)//
				__debugbreak();
#endif
			{
				mode=0;
				codebit(&ac, stats_mode+(prevmode&(LZAC_MODEHIST-1)), &mode, fwd, LZAC_LR_MODE);
				prevmode=prevmode<<1|mode;
			}
			if(mode)
			{
				count=0;
				codefixed(&ac, stats_qlen, &count, LZAC_QBITS, fwd);
				prevqlen=count;
				++count;
			}
			else
			{
				codegamma(&ac, stats_matchlen, &count, fwd);
				prevmatchlen=count;
			}
			fillend=imptr+count;
			if(fillend>imend)//end guard
			{
#ifdef _MSC_VER
				CRASH("");
#endif
				fillend=imend;
			}
			if(mode)//symbols
			{
				while(imptr<fillend)
				{
					int ctx=31^LZCNT32(estim[kc]*estim[kc]+1);
					if(ctx>LZAC_NCTX-1)
						ctx=LZAC_NCTX-1;
					int
						NW	=imptr[-rowstride-3],
						N	=imptr[          -3],
						W	=imptr[-rowstride  ];
					if(kc)//cRCT
					{
						NW	-=imptr[-rowstride-3-1];
						N	-=imptr[          -3-1];
						W	-=imptr[-rowstride  -1];
					}

					int pred=abs(N-NW)>abs(W-NW)?N:W;
					//int pred=N+W-NW, vmax=N, vmin=W;
					//if(N<W)vmin=N, vmax=W;
					//CLAMP2(pred, vmin, vmax);

					if(kc)
					{
						pred+=imptr[-1];
						CLAMP2(pred, 0, 255);
					}
					int sym=0;
					codefixed(&ac, stats_sym[3*ctx+kc], &sym, 8, fwd);
					int pix2=sym>>1^-(sym&1);
					pix2=(uint8_t)(pix2+pred);
					*imptr++=pix2;
#ifdef ENABLE_GUIDE
					if(imptr[-1]!=g_image[imptr-image-1])
					{
						int idx=(int)(imptr-image-1), kc2, ky, kx;
						{
							int idx2=idx/3;
							kc2=idx%3;
							kx=idx2%iw;
							ky=idx2/iw;
						}
						printf("guide  IDX %d  X %d  Y %d  C %d  0x%02X != 0x%02X\n"
							, idx
							, kx
							, ky
							, kc2
							, imptr[-1]
							, g_image[imptr-image-1]
						);
						CRASH("");
					}
#endif
					estim[kc]+=((sym<<LZAC_ESTIMBITS)-estim[kc]+(1<<LZAC_ESTIMLR>>1))>>LZAC_ESTIMLR;
#ifndef LZ_SINGLECHANNEL
					++kc;
					if(kc>=3)
						kc=0;
#endif
				}
			}
			else//copy
			{
				uint8_t *src;
				int backtrack=0;
				codegamma(&ac, stats_backtrack, &backtrack, fwd);
				prevbacktrack=backtrack;
				src=imptr-backtrack;
				if(src<image)//start guard
				{
#ifdef _MSC_VER
					CRASH("");
#endif
					src=image;
				}
				while(imptr<fillend)
				{
					*imptr++=*src++;
#ifdef ENABLE_GUIDE
					if(imptr[-1]!=g_image[imptr-image-1])
					{
						int idx=(int)(imptr-image-1), kc2, ky, kx;
						int srcidx=(int)(src-image-1), srcc, srcx, srcy;
						{
							int idx2=idx/3;
							kc2=idx%3;
							kx=idx2%iw;
							ky=idx2/iw;
						}
						{
							int idx2=srcidx/3;
							srcc=srcidx%3;
							srcx=idx2%iw;
							srcy=idx2/iw;
						}
						printf(
							"guide  IDX %d  Y %d  X %d  C %d  0x%02X\n"
							"src    IDX %d  Y %d  X %d  C %d  0x%02X\n"
							, idx
							, ky
							, kx
							, kc2
							, imptr[-1]
							, srcidx
							, srcy
							, srcx
							, srcc
							, g_image[imptr-image-1]
						);
						CRASH("");
					}
#endif
				}
#ifndef LZ_SINGLECHANNEL
				kc+=count;
				kc%=3;
#endif
			}
		}
#endif
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
			csize+=fwrite(&tag, 1, 2, fdst);
			csize+=fwrite(&iw, 1, 3, fdst);
			csize+=fwrite(&ih, 1, 3, fdst);
			csize+=fwrite(streamstart, 1, ac.ptr-streamstart, fdst);
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
	(void)codefixed;
	return 0;
}
