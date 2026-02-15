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
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#else
#include<time.h>
#endif
#include<immintrin.h>


#ifdef _MSC_VER
	#define LOUD
//	#define ESTIMATE_SIZE
//	#define ENABLE_GUIDE
#endif


//	#define USE_CASCADE
	#define USE_L1
//	#define ANALYSIS_SIMD


#ifdef USE_L1
#define PREDLIST\
	PRED(240, N)\
	PRED(240, W)\
	PRED(120, N+W-NW)\
	PRED(100, 2*N-NN)\
	PRED(100, NE)\
	PRED(100, NEEE)\
	PRED(100, 2*W-WW)\
	PRED(180, NNN)\
	PRED(180, WWW)\
	PRED(140, W+NE-N)\
	PRED(100, W+NEE-NE)\
	PRED(120, N+NE-NNE)\
	PRED(100, W+NW-NWW)\
	PRED(100, N+NW-NNW)\

#endif

#ifdef USE_CASCADE
#define CPREDLIST\
	PRED(240, N)\
	PRED(240, W)\
	PRED(120, N+W-NW)\
	PRED(100, 2*N-NN)\
	PRED(100, NE)\
	PRED(100, NEEE)\
	PRED(100, 2*W-WW)\
	PRED(180, NNN)\
	PRED(180, WWW)\
	PRED(140, (W+NE+N)/3)\
	PRED(100, NEE)\
	PRED(120, NNE)\
	PRED(100, NWW)\
	PRED(100, NNW)\

#endif

enum
{
#ifdef USE_L1
	SHIFT=22,
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED
#endif
#ifdef USE_CASCADE
	CSHIFT=18,
#define PRED(...) +1
	NCASCADE=CPREDLIST,
#undef  PRED
#endif
	
	GRBITS=5,
	NCTX=18,
	PREDBITS=8,
	NLEVELS=256,

	XPAD=8,
	NROWS=4,
	NVAL=3,
};

//runtime
#if 1
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
#define SWAPVAR(A, B, T) T=A, A=B, B=T
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
#define CVTFP32_I32(X) _mm_cvt_ss2si(_mm_set_ss(X))
#define CVTTFP32_I32(X) _mm_cvtt_ss2si(_mm_set_ss(X))
#define CVTFP64_I64(X) _mm_cvtsd_si64(_mm_set_sd(X))
#define CVTTFP64_I64(X) _mm_cvttsd_si64(_mm_set_sd(X))
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
static int strcmp_i(const char *s1, const char *s2)
{
	while(*s1&&*s2&&tolower(*s1)==tolower(*s2))++s1, ++s2;
	return *s1-*s2;
}
#ifdef ENABLE_GUIDE
static int g_nch=0, g_iw=0, g_ih=0;
static uint8_t *g_image=0;
static double *g_sqe=0;
static void guide_save(uint8_t *image, int nch, int iw, int ih)
{
	ptrdiff_t size=(ptrdiff_t)nch*iw*ih;
	g_nch=nch;
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
	ptrdiff_t idx=g_nch*((ptrdiff_t)g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, g_nch))
	{
		CRASH("Guide error  XY %d %d", kx, ky);
		printf("");
	}
}
//static void guide_update(uint8_t *image, int kx, int ky)
//{
//	int idx=g_nch*(g_iw*ky+kx), diff;
//	diff=g_image[idx+0]-image[idx+0]; g_sqe[0]+=diff*diff; if(abs(diff)>96)CRASH("");
//	diff=g_image[idx+1]-image[idx+1]; g_sqe[1]+=diff*diff; if(abs(diff)>96)CRASH("");
//	diff=g_image[idx+2]-image[idx+2]; g_sqe[2]+=diff*diff; if(abs(diff)>96)CRASH("");
//	diff=g_image[idx+3]-image[idx+3]; g_sqe[3]+=diff*diff; if(abs(diff)>96)CRASH("");
//}
#else
#define guide_save(...)
#define guide_check(...)
#define guide_update(...)
#endif//ENABLE_GUIDE
#endif

static int parse_files(int argc, char **argv, int nch, uint8_t *buf, int *ret_nch, int *ret_iw, int *ret_ih)
{
	int iw=-1, ih=-1, kc=0, k=0;
	const char *prevfn=0;

	for(k=3, kc=0;k<argc;++k)
	{
		FILE *f=0;
		int64_t c=0;
		int iw2=0, ih2=0;
			
		f=fopen(argv[k], "rb");
		if(!f)
		{
			CRASH("Cannot open \"%s\"", argv[k]);
			return 1;
		}
		c=fgetc(f);
		c|=(int64_t)fgetc(f)<<8;
		c|=(int64_t)fgetc(f)<<16;
		if(c==('P'|'6'<<8|'\n'<<16)||c==('P'|'5'<<8|'\n'<<16))
		{
			int color=c==('P'|'6'<<8|'\n'<<16);
			c=fgetc(f);
			while(c=='#')
			{
				c=fgetc(f);
				while(c!='\n')
					c=fgetc(f);
				c=fgetc(f);
			}
			iw2=0;
			while((uint32_t)(c-'0')<10)
			{
				iw2=10*iw2+(int32_t)c-'0';
				c=fgetc(f);
			}
			while(c<=' ')
				c=fgetc(f);
			ih2=0;
			while((uint32_t)(c-'0')<10)
			{
				ih2=10*ih2+(int32_t)c-'0';
				c=fgetc(f);
			}
			while(c=='#')
			{
				c=fgetc(f);
				while(c!='\n')
					c=fgetc(f);
				c=fgetc(f);
			}
			c|=(int64_t)fgetc(f)<<8*1;
			c|=(int64_t)fgetc(f)<<8*2;
			c|=(int64_t)fgetc(f)<<8*3;
			c|=(int64_t)fgetc(f)<<8*4;
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
			if((iw!=-1&&iw!=iw2)||(ih!=-1&&ih!=ih2))
			{
				CRASH(
					"Resolution mismatch %d*%d != %d*%d\n"
					"  \"%s\"\n"
					"  \"%s\"\n"
					, iw, ih, iw2, ih2
					, prevfn
					, argv[k]
				);
			}
			if(buf)
			{
				uint8_t *ptr=buf+kc, *end=buf+(ptrdiff_t)nch*iw2*ih2;

				if(color)
				{
					while(ptr<end)
					{
						ptr[0]=fgetc(f);
						ptr[1]=fgetc(f);
						ptr[2]=fgetc(f);
						ptr+=nch;
					}
					kc+=3;
				}
				else
				{
					while(ptr<end)
					{
						*ptr=fgetc(f);
						ptr+=nch;
					}
					++kc;
				}
			}
			else
				nch+=color?3:1;
			iw=iw2;
			ih=ih2;
			prevfn=argv[k];
		}
		else
		{
			CRASH("Unsupported file \"%s\"", argv[k]);
			return 1;
		}
		fclose(f);
	}
	if(ret_nch)*ret_nch=nch;
	if(ret_iw)*ret_iw=iw;
	if(ret_ih)*ret_ih=ih;
	return 0;
}
int c52_codec(int argc, char **argv)
{
	const uint16_t tag='5'|'2'<<8;

	const char *arname=0;
	int fwd=0, nch=0, iw=0, ih=0;
	int rowstride=0;
	ptrdiff_t usize=0, csize=0;
	uint8_t *image=0, *imptr=0, *imend=0;
	uint8_t *stream=0, *streamptr=0, *streamend=0;
	int rctsize=0;
	int8_t *rctidx=0;
	int psize=0;
	int16_t *pixels=0;
	int64_t *weights=0;
	int hsize=0;
	int16_t *hists=0;
	uint64_t low=0, range=0xFFFFFFFFFFFF, code=0;
#ifdef LOUD
	double t=0;
#endif
#ifdef ESTIMATE_SIZE
	double *csizes=0;
#endif

	if(argc<3||argv[1][1]||(uint32_t)(argv[1][0]-'c')>=2)
	{
		printf(
			"Usage:  \"%s\"  c|d  archive  files.ppm/pgm...\n"
			"Built on %s %s\n"
			, argv[0]
			, __DATE__, __TIME__
		);
		return 1;
	}
#ifdef LOUD
	t=time_sec();
#endif
	arname=argv[2];
	fwd=(argv[1][0]&0xDF)=='C';
	if(fwd)
	{
		parse_files(argc, argv, 0, 0, &nch, &iw, &ih);
		
		usize=(ptrdiff_t)nch*iw*ih;
		stream=(uint8_t*)malloc(usize);
		if(!stream)
		{
			CRASH("Alloc error");
			return 1;
		}
		streamend=stream+usize;
		streamptr=stream;
	}
	else
	{
		int error=0;
		struct stat info={0};
		FILE *f=0;
		int64_t c=0;

		error=stat(arname, &info);
		if(error)
		{
			CRASH("Canot stat \"%s\"", arname);
			return 1;
		}
		f=fopen(arname, "rb");
		if(!f)
		{
			CRASH("Cannot open \"%s\"", arname);
			return 1;
		}
		c=0;
		fread(&c, 1, 2, f);
		if(c!=tag)
		{
			CRASH("Unsupported file \"%s\"", arname);
			return 1;
		}
		nch=0;
		iw=0;
		ih=0;
		fread(&nch, 1, 2, f);
		fread(&iw, 1, 3, f);
		fread(&ih, 1, 3, f);
		csize=info.st_size;
		stream=(uint8_t*)malloc(csize);
		rctsize=nch*sizeof(int8_t);
		rctidx=(int8_t*)malloc(rctsize);
		if(!stream||!rctidx)
		{
			CRASH("Alloc error");
			return 1;
		}
		rctidx[0]=-1;
		if(nch>1)
			fread(rctidx+1, 1, rctsize-1LL, f);
		streamend=stream+csize-ftell(f);
		streamptr=stream;
		fread(stream, 1, streamend-stream, f);
		fclose(f);

		code=0;
		code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);//load
		code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);
	}
	rowstride=nch*iw;
	usize=(ptrdiff_t)nch*iw*ih;
	image=(uint8_t*)malloc(usize);
#ifdef ESTIMATE_SIZE
	csizes=(double*)malloc(nch*sizeof(double));
	if(!csizes)
	{
		CRASH("Alloc error");
		return 1;
	}
	memset(csizes, 0, nch*sizeof(double));
#endif
	if(!image)
	{
		CRASH("Alloc error");
		return 1;
	}
	imend=image+usize;
	if(fwd)
	{
		rctsize=nch*sizeof(int8_t);
		rctidx=(int8_t*)malloc(rctsize);
		if(!rctidx)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(rctidx, 0, rctsize);
		parse_files(argc, argv, nch, image, 0, 0, 0);
#ifdef ENABLE_GUIDE
		guide_save(image, nch, iw, ih);
#endif

		//simple analysis
#if 1
		int ncandidates=nch*(nch+1LL)>>1;
		int64_t *ctrs=(int64_t*)malloc(sizeof(int64_t)*ncandidates);
		int16_t *prev=(int16_t*)malloc(sizeof(int16_t)*ncandidates);
		if(!ctrs)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(ctrs, 0, sizeof(int64_t)*ncandidates);
		memset(prev, 0, sizeof(int16_t)*ncandidates);
		imptr=image+rowstride;
		while(imptr<imend)
		{
			int j=0;
			for(int kc=0;kc<nch;++kc)
			{
				int curr=imptr[kc]-imptr[kc-rowstride];
				ctrs[j]+=abs(curr-prev[j]);
				prev[j]=curr;
				++j;
			}
			for(int k1=0;k1<nch-1;++k1)
			{
				for(int k2=k1+1;k2<nch;++k2)
				{
					int curr=prev[k2]-prev[k1];
					ctrs[j]+=abs(curr-prev[j]);
					prev[j]=curr;
					++j;
				}
			}
			imptr+=nch;
		}
		rctidx[0]=-1;
		for(int kc=1;kc<nch;++kc)
		{
			int64_t ebest=0;
			int kbest=0;
			for(int kc2=-1, idx=kc, stride=nch;kc2<kc;++kc2, --stride, idx+=stride)
			{
				int64_t curr=ctrs[idx];
#ifdef LOUD
				printf("%2d  %16lld%s\n", kc, curr, kc2==-1||ebest>curr?" <-":"");
#endif
				if(kc2==-1||ebest>curr)
					ebest=curr, kbest=kc2;
			}
#ifdef LOUD
			printf("\t%2d\n", kbest);
#endif
			rctidx[kc]=kbest;
		}
#ifdef LOUD
		printf("\n");
		for(int kc=0;kc<nch;++kc)
		{
			printf("%2d", kc);
			if(rctidx[kc]>=0)
				printf("  %2d", rctidx[kc]);
			else
				printf("  %2s", "");
			printf("\n");
		}
#endif
#endif
		free(ctrs);
		free(prev);
	}
	psize=(iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NVAL])*nch;
	pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i));
	weights=(int64_t*)malloc(nch*sizeof(int64_t[NPREDS]));
	hsize=nch*sizeof(int16_t[NCTX*(NLEVELS+1)<<PREDBITS]);
	hists=(int16_t*)malloc(hsize);
	if(!pixels||!weights||!hists)
	{
		CRASH("Alloc error");
		return 1;
	}
	memset(pixels, 0, psize);
	for(int k=0;k<nch*NPREDS;++k)
		weights[k]=(1LL<<SHIFT)/NPREDS;
	memset(hists, 0, hsize);
	imptr=image;
	for(int ky=0;ky<ih;++ky)
	{
		int estim[NPREDS]={0};
		int16_t *rows[]=
		{
			pixels+(XPAD*nch*NROWS+(ky-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*nch*NROWS+(ky-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*nch*NROWS+(ky-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*nch*NROWS+(ky-3LL+NROWS)%NROWS)*NVAL,
		};
		for(int kx=0;kx<iw;++kx, imptr+=nch)
		{
			int offset;
			for(int kc=0;kc<nch;++kc)
			{
				int16_t
					NNN	=rows[3][0+0*nch*NROWS*NVAL],
					NNWW	=rows[2][0-2*nch*NROWS*NVAL],
					NNW	=rows[2][0-1*nch*NROWS*NVAL],
					NN	=rows[2][0+0*nch*NROWS*NVAL],
					NNE	=rows[2][0+1*nch*NROWS*NVAL],
					NNEE	=rows[2][0+2*nch*NROWS*NVAL],
					NWW	=rows[1][0-2*nch*NROWS*NVAL],
					NW	=rows[1][0-1*nch*NROWS*NVAL],
					N	=rows[1][0+0*nch*NROWS*NVAL],
					NE	=rows[1][0+1*nch*NROWS*NVAL],
					NEE	=rows[1][0+2*nch*NROWS*NVAL],
					NEEE	=rows[1][0+3*nch*NROWS*NVAL],
					NEEEE	=rows[1][0+4*nch*NROWS*NVAL],
					WWWW	=rows[0][0-4*nch*NROWS*NVAL],
					WWW	=rows[0][0-3*nch*NROWS*NVAL],
					WW	=rows[0][0-2*nch*NROWS*NVAL],
					W	=rows[0][0-1*nch*NROWS*NVAL],

					eNEE	=rows[1][2+2*nch*NROWS*NVAL],
					eNEEE	=rows[1][2+3*nch*NROWS*NVAL],
					eW	=rows[0][2-1*nch*NROWS*NVAL];
				int64_t p1, *currweights;
				int j, ctx, pred, vmax, vmin, error, sym, curr;
				
				ctx=FLOOR_LOG2(eW*eW+1);
				if(ctx>NCTX-1)
					ctx=NCTX-1;

				offset=0;
				if(rctidx[kc]>=0)
					offset=imptr[rctidx[kc]];
				currweights=weights+NPREDS*kc;
				p1=1LL<<SHIFT>>1;
#define PRED(W0, E) estim[j]=E; p1+=currweights[j]*estim[j]; ++j;
				j=0;
				PREDLIST;
#undef  PRED
				p1>>=SHIFT;
				pred=(int)p1;
				vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);
				pred+=offset;
				CLAMP2(pred, 0, 255);

				uint16_t *currhist=hists+(NLEVELS+1LL)*((NCTX*(ptrdiff_t)kc+ctx)<<PREDBITS|pred>>(8-PREDBITS));
				int den=currhist[NLEVELS]+NLEVELS;
				int cdf=0, freq, tmp;
				if(fwd)
				{
					error=(int8_t)(imptr[kc]-pred);
					sym=error<<1^error>>31;
					if(range<=0xFFFF)
					{
#ifdef _MSC_VER
						if(streamptr+sizeof(uint32_t)>streamend)
						{
							CRASH("Inflation %8.4lf%%", 100.*ih/(ky+1));
							return 1;
						}
#endif
						*(uint32_t*)streamptr=(uint32_t)(low>>32);
						streamptr+=sizeof(uint32_t);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					for(tmp=0;;++tmp)
					{
						freq=currhist[tmp]+1;
						if(tmp>=sym)
							break;
						cdf+=freq;
					}
					low+=range*cdf/den;
					range=range*freq/den-1;
#ifdef ESTIMATE_SIZE
					csizes[kc]-=log2((double)freq/den);
#endif
				}
				else
				{
					if(range<=0xFFFF)
					{
						code=code<<32|*(uint32_t*)streamptr;
						streamptr+=sizeof(uint32_t);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					tmp=(int)(((code-low+1)*den-1)/range);
					for(sym=0;;++sym)
					{
						int cdf2;
						
						freq=currhist[sym]+1;
						cdf2=cdf+freq;
						if(cdf2>tmp)
							break;
						cdf=cdf2;
					}
					low+=range*cdf/den;
					range=range*freq/den-1;
					error=sym>>1^-(sym&1);
					imptr[kc]=(uint8_t)(error+pred);
#ifdef ENABLE_GUIDE
					if(g_image[imptr+kc-image]!=imptr[kc])
					{
						CRASH("Guide error");
						return 1;
					}
#endif
				}
				++currhist[sym];
				++currhist[NLEVELS];
				if(currhist[NLEVELS]>=0xFFFF-2*NLEVELS)
				{
					den=0;
					for(int ks=0;ks<NLEVELS;++ks)
						den+=currhist[ks]>>=1;
					currhist[NLEVELS]=den;
				}

				rows[0][0]=curr=imptr[kc]-offset;
				rows[0][2]=(2*eW+(sym<<GRBITS)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
				{
					error=(curr>p1)-(curr<p1);
#define PRED(...) currweights[j]+=error*estim[j]; ++j;
					j=0;
					PREDLIST;
#undef  PRED
				}
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
			}
		}
	}
	_mm_free(pixels);
	free(weights);
	free(hists);
	if(fwd)
	{
		FILE *fdst=fopen(arname, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing", arname);
			return 1;
		}
		*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;//flush
		*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;

		csize=0;
		csize+=fwrite(&tag, 1, 2, fdst);
		csize+=fwrite(&nch, 1, 2, fdst);
		csize+=fwrite(&iw, 1, 3, fdst);
		csize+=fwrite(&ih, 1, 3, fdst);
		if(nch>1)
			csize+=fwrite(rctidx+1, 1, rctsize-1LL, fdst);
		csize+=fwrite(stream, 1, streamptr-stream, fdst);
		fclose(fdst);
	}
	else
	{
		int channelidx=0;
		int noutputchannels=0;
		for(int k=3;k<argc;++k)
		{
			int nch2=0;
			const char *fn=argv[k];
			char *ext=strrchr(fn, '.');
			if(!ext)
			{
				CRASH("Output files must have one of these extensions: \'.PGM\' or \'.PPM\'");
				return 1;
			}
			if(!strcmp_i(ext, ".PGM"))nch2=1;
			else if(!strcmp_i(ext, ".PPM"))nch2=3;
			else
			{
				CRASH("Output files must have one of these extensions: \'.PGM\' or \'.PPM\'");
				return 1;
			}
			noutputchannels+=nch2;
		}
		if(noutputchannels!=nch)
		{
			CRASH("Output channel count (%d) != original decoded channel count (%d)", noutputchannels, nch);
			return 1;
		}
		for(int k=3;k<argc;++k)
		{
			int nch2=0;
			const char *fn=argv[k];
			char *ext=strrchr(fn, '.');
			if(!strcmp_i(ext, ".PGM"))nch2=1;
			else if(!strcmp_i(ext, ".PPM"))nch2=3;

			FILE *fdst=fopen(fn, "wb");
			if(!fdst)
			{
				CRASH("Cannot open \"%s\" for writing", fn);
				return 1;
			}
			imptr=image+channelidx;
			if(nch2==1)
			{
				fprintf(fdst, "P5\n%d %d\n255\n", iw, ih);
				while(imptr<imend)
				{
					fputc(*imptr, fdst);
					imptr+=nch;
				}
			}
			else
			{
				fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
				while(imptr<imend)
				{
					fputc(imptr[0], fdst);
					fputc(imptr[1], fdst);
					fputc(imptr[2], fdst);
					imptr+=nch;
				}
			}
			fclose(fdst);
			channelidx+=nch2;
		}
	}
	free(rctidx);
	free(image);
	free(stream);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
#ifdef ESTIMATE_SIZE
		for(int kc=0;kc<nch;++kc)
			printf("%2d  %12.2lf\n", kc, csizes[kc]/8);
#endif
		printf("CWH=3*%d*%d\n", iw, ih);
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
	(void)csize;
	(void)&time_sec;
	return 0;
}
