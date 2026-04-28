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
	#define LOUD
	#define ENABLE_GUIDE
#endif


//	#define LZ_SINGLECHANNEL


//runtime
#if 1
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
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
#endif

enum
{
	EBITS=14,
	ESIZE=(1<<EBITS),

	LZAC_MIN=5,		//tune
	LZAC_QBITS=6,		//tune
	LZAC_QMAX=(1<<LZAC_QBITS),

	LZAC_STOREBITS=24,
	LZAC_USEBITS=10,	//tune
	LZAC_SHIFT=LZAC_STOREBITS-LZAC_USEBITS,

	LZAC_ESTIMBITS=3,
	LZAC_ESTIMLR=5,		//tune

	LZAC_MODEHIST=8,
	LZAC_GAMMA_L=28,
	LZAC_NCTX=8,

	LZAC_LR_MODE=6,		//tune
	LZAC_LR_FIXED=5,
	LZAC_LR_GAMMA_U=6,
	LZAC_LR_GAMMA_B=7,

//	BUFSIZE=512*1024,	//buffering/streaming is incompatible with LZ
};
typedef struct _SymCtx
{
	uint8_t sym, ctx;
} SymCtx;
typedef struct _ETable
{
	int32_t etable[ESIZE], estart, eend, ecount;
} ETable;
static ETable tables[0x100];
typedef struct _ACState
{
	uint64_t low, range, code;
	uint8_t *ptr, *end;
} ACState;
#if 0
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
#endif
INLINE void codebit(ACState *ac, int32_t *pp0, int32_t *bit, const int fwd, const int lr)
{
	int rbit;
	uint64_t r2, mid;
	
	int32_t p0r=*pp0;
	int32_t p0=p0r>>LZAC_SHIFT;

	p0+=p0<0;
	p0+=1<<LZAC_USEBITS>>1;
	if(ac->range<=0xFFFF)
	{
		if(ac->ptr>=ac->end)
		{
#ifdef _MSC_VER
			CRASH("inflation\n");
#endif
			return;
		}
		if(fwd)
			*(uint32_t*)ac->ptr=(uint32_t)(ac->low>>32);
		else
			ac->code=ac->code<<32|*(uint32_t*)ac->ptr;
		ac->ptr+=sizeof(uint32_t);
		ac->low<<=32;
		ac->range=ac->range<<32|0xFFFFFFFF;
		if(ac->range>~ac->low)
			ac->range=~ac->low;
	}
	r2=ac->range*p0>>LZAC_USEBITS;
	mid=ac->low+r2;
	ac->range-=r2;
	--r2;
	if(fwd)
		rbit=*bit;
	else
		*bit=rbit=ac->code>=mid;
	if(rbit)
		ac->low=mid;
	else
		ac->range=r2;

	p0r+=((!rbit<<LZAC_STOREBITS)-(1<<LZAC_STOREBITS>>1)-p0r)>>lr;
	*pp0=p0r;
}
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
		nbits=FLOOR_LOG2(sym);
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
static int32_t
	stats_mode[LZAC_MODEHIST],
	stats_matchlen[2*LZAC_GAMMA_L+1],
	stats_backtrack[2*LZAC_GAMMA_L+1],
	stats_qlen[LZAC_QMAX],
	stats_sym[3*LZAC_NCTX][256];
int c56_codec(int argc, char **argv)
{
	const int16_t tag='5'|'6'<<8;

	const char *srcfn=0, *dstfn=0;
	int fwd=0, iw=0, ih=0;
	int rowstride=0;
	int estim[]={32, 32, 32};
	int prevmode=1, prevqlen=0, prevmatchlen=0, prevbacktrack=0;
#ifdef LZ_SINGLECHANNEL
	const
#endif
	int kc=0;
	ptrdiff_t usize=0, csize=0, headersize=0, cap=0;
	uint8_t *buf=0, *image=0, *streamstart=0, *streamend=0;
	ACState ac={0};
	uint8_t *imptr=0, *imend=0;
#ifdef ENABLE_GUIDE
	static uint8_t *im0=0;
#endif
#ifdef LOUD
	double t=0, t2=0;
#endif

	(void)prevqlen;
	(void)prevmatchlen;
	(void)prevbacktrack;
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
	memset(stats_mode, 0, sizeof(stats_mode));
	memset(stats_matchlen, 0, sizeof(stats_matchlen));
	memset(stats_backtrack, 0, sizeof(stats_backtrack));
	memset(stats_qlen, 0, sizeof(stats_qlen));
	memset(stats_sym, 0, sizeof(stats_sym));
	ac.range=0xFFFFFFFFFFFF;
	ac.ptr=streamstart;
	ac.end=streamend;
	imptr=image;
	imend=image+usize;
	if(fwd)
	{
		const int queuesize=sizeof(SymCtx[LZAC_QMAX]);
		SymCtx *queue=(SymCtx*)malloc(queuesize);
		int qcount=0, qstart=0, qend=0;
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
			ETable *table=tables+pixel;
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
				int ctx=FLOOR_LOG2(estim[kc]*estim[kc]+1);
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
	}
	else
	{
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
					int ctx=FLOOR_LOG2(estim[kc]*estim[kc]+1);
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
	return 0;
}
