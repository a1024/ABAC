#include"e2.h"
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctype.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;


//debug
//	#define T51_ENABLE_GUIDE

//efficiency
//	#define DYNAMIC_DEN
//	#define ENABLE_RCT_MA
//	#define RECT_BUFFER
//	#define ENABLE_SSE

//speed
//	#define PROFILER


#define NROWS 4
#define PAD_SIZE 4
//#define MAXQBITS 6

#if 1
#define T51_NPREDS 26
#define T51_PREDLIST_SELF(C)\
	T51_PRED(N(C))\
	T51_PRED(W(C))\
	T51_PRED(NW(C))\
	T51_PRED(NE(C))\
	T51_PRED(NEE(C))\
	T51_PRED((N(C)+W(C))>>1)\
	T51_PRED(N(C)+W(C)-NW(C))\
	T51_PRED(W(C)+NE(C)-N(C))\
	T51_PRED((W(C)+NEE(C))>>1)\
	T51_PRED((3*W(C)+NEEE(C))>>2)\
	T51_PRED((4*N(C)-2*NN(C)+NW(C)+NE(C))>>2)\
	T51_PRED(N(C)+NE(C)-NNE(C))\
	T51_PRED((4*(N(C)+W(C)+NW(C)+NE(C))-(NN(C)+WW(C)+NNWW(C)+NNEE(C))+6)/12)\
	T51_PRED(N(C)+W(C)+NNWW(C)-NNW(C)-NWW(C))\
	T51_PRED(geomean)\
	T51_PRED(paper_GAP)\
	T51_PRED(calic_GAP)

#define T51_PREDLIST_C0_FROM_C1
//	T51_PRED(N(0)+W(1)-NW(1))
//	T51_PRED(N(1))
//	T51_PRED(W(1))
//	T51_PRED(NW(1))
//	T51_PRED(NE(1))
//	T51_PRED((N(0)+W(0)+NW(0)+NE(0)+N(1)+W(1)+NW(1)+NE(1))>>3)
#define T51_PREDLIST_C0_FROM_C2
//	T51_PRED(N(0)+W(2)-NW(2))
//	T51_PRED(N(2))
//	T51_PRED(W(2))
//	T51_PRED(NW(2))
//	T51_PRED(NE(2))
//	T51_PRED((N(0)+W(0)+NW(0)+NE(0)+N(2)+W(2)+NW(2)+NE(2))>>3)
#define T51_PREDLIST_C0_FROM_C3
//	T51_PRED(N(0)+W(3)-NW(3))
//	T51_PRED(N(3))
//	T51_PRED(W(3))
//	T51_PRED(NW(3))
//	T51_PRED(NE(3))
//	T51_PRED((N(0)+W(0)+NW(0)+NE(0)+N(3)+W(3)+NW(3)+NE(3))>>3)

#define T51_PREDLIST_C1_FROM_C0\
	T51_PRED(CURR(0)+W(1)-W(0))\
	T51_PRED(CURR(0)+N(1)-N(0))\
	T51_PRED(CURR(0))
//	T51_PRED(CURR(0)+N(1)+W(1)+NW(0)-(N(0)+W(0)+NW(1)))
//	T51_PRED(CURR(0)+((NW(1)+NE(1)-NW(0)-NE(0))>>1))
//	T51_PRED(N(0))
//	T51_PRED(W(0))
//	T51_PRED(NW(0))
//	T51_PRED(NE(0))
//	T51_PRED((N(1)+W(1)+NW(1)+NE(1)+N(0)+W(0)+NW(0)+NE(0))>>3)
#define T51_PREDLIST_C1_FROM_C2
//	T51_PRED(N(1)+W(2)-NW(2))
//	T51_PRED(N(2))
//	T51_PRED(W(2))
//	T51_PRED(NW(2))
//	T51_PRED(NE(2))
//	T51_PRED((N(1)+W(1)+NW(1)+NE(1)+N(2)+W(2)+NW(2)+NE(2))>>3)
#define T51_PREDLIST_C1_FROM_C3
//	T51_PRED(N(1)+W(3)-NW(3))
//	T51_PRED(N(3))
//	T51_PRED(W(3))
//	T51_PRED(NW(3))
//	T51_PRED(NE(3))
//	T51_PRED((N(1)+W(1)+NW(1)+NE(1)+N(3)+W(3)+NW(3)+NE(3))>>3)

#define T51_PREDLIST_C2_FROM_C0_C1\
	T51_PRED(CURR(0)+W(2)-W(0))\
	T51_PRED(CURR(1)+W(2)-W(1))\
	T51_PRED(CURR(0)+N(2)-N(0))\
	T51_PRED(CURR(1)+N(2)-N(1))\
	T51_PRED(CURR(0))\
	T51_PRED(CURR(1))\
	T51_PRED((CURR(0)+CURR(1))>>1)
//	T51_PRED(CURR(0)+N(2)+W(2)+NW(0)-(N(0)+W(0)+NW(2)))
//	T51_PRED(CURR(1)+N(2)+W(2)+NW(1)-(N(1)+W(1)+NW(2)))
//	T51_PRED(N(0))
//	T51_PRED(W(0))
//	T51_PRED(NW(0))
//	T51_PRED(NE(0))
//	T51_PRED((N(2)+W(2)+NW(2)+NE(2)+N(0)+W(0)+NW(0)+NE(0))>>3)
//	T51_PRED(N(1))
//	T51_PRED(W(1))
//	T51_PRED(NW(1))
//	T51_PRED(NE(1))
//	T51_PRED((N(2)+W(2)+NW(2)+NE(2)+N(1)+W(1)+NW(1)+NE(1))>>3)
#define T51_PREDLIST_C2_FROM_C3
//	T51_PRED(N(2)+W(3)-NW(3))
//	T51_PRED((N(2)+W(2)+NW(2)+NE(2)+N(3)+W(3)+NW(3)+NE(3))>>3)
//	T51_PRED(N(3))
//	T51_PRED(W(3))
//	T51_PRED(NW(3))
//	T51_PRED(NE(3))

#define T51_PREDLIST_C3_FROM_C0_C1_C2\
	T51_PRED(W(3)+CURR(0)-W(0))\
	T51_PRED(W(3)+CURR(1)-W(1))\
	T51_PRED(W(3)+CURR(2)-W(2))\
	T51_PRED(N(3)+CURR(0)-N(0))\
	T51_PRED(N(3)+CURR(1)-N(1))\
	T51_PRED(N(3)+CURR(2)-N(2))\
	T51_PRED(CURR(0))\
	T51_PRED(CURR(1))\
	T51_PRED(CURR(2))
//	T51_PRED((N(3)+W(3)+NW(3)+NE(3)+N(0)+W(0)+NW(0)+NE(0))>>3)
//	T51_PRED(N(0))
//	T51_PRED(W(0))
//	T51_PRED(NW(0))
//	T51_PRED(NE(0))
//	T51_PRED((N(3)+W(3)+NW(3)+NE(3)+N(1)+W(1)+NW(1)+NE(1))>>3)
//	T51_PRED(N(1))
//	T51_PRED(W(1))
//	T51_PRED(NW(1))
//	T51_PRED(NE(1))
//	T51_PRED((N(3)+W(3)+NW(3)+NE(3)+N(2)+W(2)+NW(2)+NE(2))>>3)
//	T51_PRED(N(2))
//	T51_PRED(W(2))
//	T51_PRED(NW(2))
//	T51_PRED(NE(2))
#endif

#ifdef PROFILER
#define CHECKPOINTLIST\
	CHECKPOINT(OUTSIDE)\
	CHECKPOINT(GETCTX)\
	CHECKPOINT(PRED_MEM)\
	CHECKPOINT(PRED_SQRT)\
	CHECKPOINT(UPDATE)\
	CHECKPOINT(CSIZE)
typedef enum ProfilerLabelEnum
{
#define CHECKPOINT(X) PROF_##X,
	CHECKPOINTLIST
#undef  CHECKPOINT
	PROF_COUNT,
} ProfilerLabel;
static const char *prof_labels[]=
{
#define CHECKPOINT(X) #X,
	CHECKPOINTLIST
#undef  CHECKPOINT
};
static long long prof_timestamp=0, prof_cycles[PROF_COUNT]={0};
#define PROF_START() memset(prof_cycles, 0, _countof(prof_cycles)), prof_timestamp=__rdtsc()
#define PROF(X) prof_cycles[PROF_##X]+=__rdtsc()-prof_timestamp, prof_timestamp=__rdtsc()
static void prof_print()
{
	long long sum=0;
	int maxlen=0;
	for(int k=0;k<_countof(prof_labels);++k)
	{
		int len=(int)strlen(prof_labels[k]);
		UPDATE_MAX(maxlen, len);
		sum+=prof_cycles[k];
	}
	printf("Profiler:\n");
	for(int k=0;k<_countof(prof_labels);++k)
	{
		double percent=100.*prof_cycles[k]/sum;
		printf("%-*s %16lld %6.2lf%%  ", maxlen, prof_labels[k], prof_cycles[k], percent);
		for(int k2=0, npoints=(int)percent;k2<npoints;++k2)
			printf("*");
		printf("\n");
	}
}
#else
#define PROF_START()
#define PROF(...)
#define prof_print()
#endif
typedef short StatNode;
typedef struct T51CtxStruct
{
	int iw, ih, nch;
	char depths[4];
	int
		maxdepth, maxlevels,
		nlevels[4], half[4],
#ifdef RECT_BUFFER
		partialmask,
#else
		nqlevels[4], lgdepths[4], maxqlevels, (*quantize[4])(int),
		//shiftq[4],
#endif
		shift8[4];
	int *pixels;

	int preds[T51_NPREDS], qpreds[T51_NPREDS], npreds;
#ifdef ENABLE_SSE
	long long *sse;//nch * npreds * nqlevels
	long long *sse_cells[T51_NPREDS];
#endif
	StatNode *stats;//nch * npreds << maxdepth
	StatNode *stat_cells[T51_NPREDS];
	int *mixer;
#ifdef DYNAMIC_DEN
	long long csize_approx[128];
	//double csize_actual[128];
#endif
	
#ifdef ENABLE_SSE
	long long *p0_bias;
#endif
	long long p0;
	int lglrden;

	int kc, kx, ky, bitidx;
#ifndef RECT_BUFFER
	int treeidx;
#endif
	int partial, kym[4];
	ArithmeticCoder *ec;
} T51Ctx;
#ifndef RECT_BUFFER
static int quantize0(int x)
{
	int negmask=-(x<0);
	x=abs(x);
	int qx=floor_log2_32(x)+1;
	qx^=negmask;
	qx-=negmask;
	return qx;
}
static int quantize1(int x)
{
	int negmask=-(x<0);
	x=abs(x);
	int qx=floor_log2_32(x)+1;
	x>>=qx-2;
	x&=1;
	qx<<=1;
	qx|=x;
	qx^=negmask;
	qx-=negmask;
	return qx;
}
static int quantize2(int x)
{
	int negmask=-(x<0);
	x=abs(x);
	int qx=floor_log2_32(x)+1;
	qx=floor_log2_32(x)+1;
	qx^=negmask;
	qx-=negmask;
	return qx;
}
static int quantize3(int x)
{
	int negmask=-(x<0);
	x=abs(x);
	int qx=floor_log2_32(x)+1;
	qx=floor_log2_32(x)+1;
	qx=floor_log2_32(x)+1;
	qx^=negmask;
	qx-=negmask;
	return qx;
}
static int (*const quantize[])(int)=
{
	quantize0,
	quantize0,
	quantize1,
	quantize2,
	quantize3,
};
#endif
static ptrdiff_t t51_init(T51Ctx *pr, int iw, int ih, int nch, const char *depths)//returns memory usage or 0 on failure
{
	PROF_START();
	if(iw<1||ih<1||nch<1)
	{
		LOG_ERROR("Invalid image");
		return 0;
	}
	memset(pr, 0, sizeof(*pr));
	pr->iw=iw;
	pr->ih=ih;
	pr->nch=nch;
	memcpy(pr->depths, depths, nch);
	pr->maxdepth=0;
#ifndef RECT_BUFFER
	pr->maxqlevels=0;
#endif
	for(int kc=0;kc<nch;++kc)
	{
		pr->nlevels[kc]=1<<pr->depths[kc];
		pr->half[kc]=pr->nlevels[kc]>>1;
		pr->shift8[kc]=pr->depths[kc]-8;
		UPDATE_MAX(pr->maxdepth, pr->depths[kc]);
#ifndef RECT_BUFFER
		pr->lgdepths[kc]=floor_log2_32(pr->depths[kc]);
		UPDATE_MIN(pr->lgdepths[kc], 4);
		pr->quantize[kc]=quantize[pr->lgdepths[kc]];
		pr->nqlevels[kc]=pr->quantize[kc](pr->half[kc])<<1|1;
		//pr->nqlevels[kc]=quantize(MINVAR(pr->half[kc], (1<<MAXQBITS)>>1))<<1|1;
		//pr->shiftq[kc]=MAXVAR(pr->depths[kc], MAXQBITS)-MAXQBITS;
		UPDATE_MAX(pr->maxqlevels, pr->nqlevels[kc]);
#endif
	}
	pr->maxlevels=1<<pr->maxdepth;

	size_t pixels_size=(iw+PAD_SIZE*2LL)*sizeof(int[4*NROWS]);
#ifdef RECT_BUFFER
	size_t stats_size=(size_t)pr->nch*pr->maxdepth*sizeof(StatNode[T51_NPREDS])<<pr->maxdepth;
#else
	size_t stats_size=(size_t)pr->nch*pr->maxqlevels*sizeof(StatNode[T51_NPREDS])<<pr->maxdepth;
#endif
	size_t mixer_size=(size_t)pr->nch*pr->maxdepth*sizeof(int[T51_NPREDS]);

	size_t memusage=0, size;
#define ALLOC(PTR, PTR_TYPE, SIZE) size=SIZE, memusage+=size, PTR=(PTR_TYPE)malloc(size)
	ALLOC(pr->pixels, int*, pixels_size);
	ALLOC(pr->stats, StatNode*, stats_size);
	ALLOC(pr->mixer, int*, mixer_size);
#ifdef ENABLE_SSE
	ALLOC(pr->sse, long long*, (size_t)pr->nch*pr->maxqlevels*sizeof(long long[T51_NPREDS]));
	ALLOC(pr->p0_bias, long long*, (size_t)pr->nch*pr->maxdepth*sizeof(long long));
#endif
#undef  ALLOC
	if(!pr->pixels||!pr->stats||!pr->mixer
#ifdef ENABLE_SSE
		||!pr->sse||!pr->p0_bias
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(pr->pixels, 0, pixels_size);
	memset(pr->stats, 0, stats_size);
	memset(pr->mixer, 0, mixer_size);
#ifdef ENABLE_SSE
	memset(pr->sse, 0, (size_t)pr->nch*pr->maxqlevels*sizeof(long long[T51_NPREDS]));
	memset(pr->p0_bias, 0, (size_t)pr->nch*pr->maxdepth*sizeof(long long));
#endif
	return memusage;
}
static void t51_free(T51Ctx *pr)
{
	free(pr->pixels);
	free(pr->stats);
	free(pr->mixer);
#ifdef ENABLE_SSE
	free(pr->sse);
	free(pr->p0_bias);
#endif
}
static void t51_nextrow(T51Ctx *pr, int ky)
{
	pr->ky=ky;
	pr->kym[0]=(pr->iw+PAD_SIZE*2)*(ky&3);
	pr->kym[1]=(pr->iw+PAD_SIZE*2)*((ky-1)&3);
	pr->kym[2]=(pr->iw+PAD_SIZE*2)*((ky-2)&3);
	pr->kym[3]=(pr->iw+PAD_SIZE*2)*((ky-3)&3);
}
#define LOAD(BUF, X, Y) BUF[(pr->kym[Y]+pr->kx+PAD_SIZE-(X))<<2|pr->kc]
#if 1
#define NNNWWWW(C)	NNNrow[C-16]
#define NNNWWW(C)	NNNrow[C-12]
#define NNNWW(C)	NNNrow[C-8]
#define NNNW(C)		NNNrow[C-4]
#define NNN(C)		NNNrow[C]
#define NNNE(C)		NNNrow[C+4]
#define NNNEE(C)	NNNrow[C+8]
#define NNNEEE(C)	NNNrow[C+13]
#define NNNEEEE(C)	NNNrow[C+17]
#define NNWWWW(C)	NNrow[C-16]
#define NNWWW(C)	NNrow[C-12]
#define NNWW(C)		NNrow[C-8]
#define NNW(C)		NNrow[C-4]
#define NN(C)		NNrow[C]
#define NNE(C)		NNrow[C+4]
#define NNEE(C)		NNrow[C+8]
#define NNEEE(C)	NNrow[C+13]
#define NNEEEE(C)	NNrow[C+17]
#define NWWWW(C)	Nrow[C-16]
#define NWWW(C)		Nrow[C-12]
#define NWW(C)		Nrow[C-8]
#define NW(C)		Nrow[C-4]
#define N(C)		Nrow[C]
#define NE(C)		Nrow[C+4]
#define NEE(C)		Nrow[C+8]
#define NEEE(C)		Nrow[C+13]
#define NEEEE(C)	Nrow[C+17]
#define WWWW(C)		currrow[C-16]
#define WWW(C)		currrow[C-12]
#define WW(C)		currrow[C-8]
#define W(C)		currrow[C-4]
#define CURR(C)		currrow[C]
//#define WWWW(X) X[-16]
//#define WWW(X) X[-12]
//#define WW(X) X[-8]
//#define W(X) X[-4]
//#define CURR(X) X[0]
//#define E(X) X[4]
//#define EE(X) X[8]
//#define EEE(X) X[12]
//#define EEEE(X) X[16]
#endif
static void t51_getctx(T51Ctx *pr, int kc, int kx)
{
	PROF(OUTSIDE);
	//int idx=(pr->iw*pr->ky+kx)<<2|kc;
	pr->kc=kc;
	pr->kx=kx;
	//int NNNrow[9*4], NNrow[9*4], Nrow[9*4], currrow[4*4];
	//memcpy(NNNrow,	pr->pixels+((pr->kym[3]+pr->kx)<<2), sizeof(int[9*4]));
	//memcpy(NNrow,		pr->pixels+((pr->kym[2]+pr->kx)<<2), sizeof(int[9*4]));
	//memcpy(Nrow,		pr->pixels+((pr->kym[1]+pr->kx)<<2), sizeof(int[9*4]));
	//memcpy(currrow,	pr->pixels+((pr->kym[0]+pr->kx)<<2), sizeof(int[4*4]));
	int
	//	*NNNrow	=pr->pixels+(((size_t)pr->kym[3]+pr->kx+PAD_SIZE)<<2),
		*NNrow	=pr->pixels+(((size_t)pr->kym[2]+pr->kx+PAD_SIZE)<<2),
		*Nrow	=pr->pixels+(((size_t)pr->kym[1]+pr->kx+PAD_SIZE)<<2),
		*currrow=pr->pixels+(((size_t)pr->kym[0]+pr->kx+PAD_SIZE)<<2);
	//int
	//	*cNNN	=pr->pixels+(((size_t)pr->kym[3]+pr->kx+PAD_SIZE)<<2),
	//	*cNN	=pr->pixels+(((size_t)pr->kym[2]+pr->kx+PAD_SIZE)<<2),
	//	*cN	=pr->pixels+(((size_t)pr->kym[1]+pr->kx+PAD_SIZE)<<2),
	//	*ccurr	=pr->pixels+(((size_t)pr->kym[0]+pr->kx+PAD_SIZE)<<2),
	//	*NNN	=cNNN+kc,
	//	*NN	=cNN+kc,
	//	*N	=cN+kc,
	//	*curr	=ccurr+kc;
#if 0
	int
#if PAD_SIZE>=4
		NNNWWWW	=LOAD(pr->pixels,  4, 3),
		NNNWWW	=LOAD(pr->pixels,  3, 3),
		NNNWW	=LOAD(pr->pixels,  2, 3),
		NNNW	=LOAD(pr->pixels,  1, 3),
		NNN	=LOAD(pr->pixels,  0, 3),
		NNNE	=LOAD(pr->pixels, -1, 3),
		NNNEE	=LOAD(pr->pixels, -2, 3),
		NNNEEE	=LOAD(pr->pixels, -3, 3),
		NNNEEEE	=LOAD(pr->pixels, -4, 3),
		
		NNWWWW	=LOAD(pr->pixels,  4, 2),
		NNWWW	=LOAD(pr->pixels,  3, 2),
#endif
		NNWW	=LOAD(pr->pixels,  2, 2),
		NNW	=LOAD(pr->pixels,  1, 2),
		NN	=LOAD(pr->pixels,  0, 2),
		NNE	=LOAD(pr->pixels, -1, 2),
		NNEE	=LOAD(pr->pixels, -2, 2),
#if PAD_SIZE>=4
		NNEEE	=LOAD(pr->pixels, -3, 2),
		NNEEEE	=LOAD(pr->pixels, -4, 2),
		NWWWW	=LOAD(pr->pixels,  4, 1),
		NWWW	=LOAD(pr->pixels,  3, 1),
#endif
		NWW	=LOAD(pr->pixels,  2, 1),
		NW	=LOAD(pr->pixels,  1, 1),
		N	=LOAD(pr->pixels,  0, 1),
		NE	=LOAD(pr->pixels, -1, 1),
		NEE	=LOAD(pr->pixels, -2, 1),
#if PAD_SIZE>=4
		NEEE	=LOAD(pr->pixels, -3, 1),
		NEEEE	=LOAD(pr->pixels, -4, 1),
		WWWW	=LOAD(pr->pixels,  4, 0),
		WWW	=LOAD(pr->pixels,  3, 0),
#endif
		WW	=LOAD(pr->pixels,  2, 0),
		W	=LOAD(pr->pixels,  1, 0);
#endif
	int sh=pr->shift8[kc];

	int dx=abs(W(kc)-WW(kc))+abs(N(kc)-NW(kc))+abs(NE(kc)-N(kc));
	int dy=abs(W(kc)-NW(kc))+abs(N(kc)-NN(kc))+abs(NE(kc)-NNE(kc));
	int d45=abs(W(kc)-NWW(kc))+abs(NW(kc)-NNWW(kc))+abs(N(kc)-NNW(kc));
	int d135=abs(NE(kc)-NNEE(kc))+abs(N(kc)-NNE(kc))+abs(W(kc)-N(kc));
	int diff=SHIFT_RIGHT_SIGNED(dy-dx, sh), diff2=SHIFT_RIGHT_SIGNED(d45-d135, sh), diff3=NE(kc)-NW(kc);
	int paper_GAP, calic_GAP;
	if(dy+dx>32)
		paper_GAP=(int)(((long long)dx*N(kc)+(long long)dy*W(kc))/((long long)dy+dx));
	else if(diff>12)
		paper_GAP=(N(kc)+2*W(kc))/3;
	else if(diff<-12)
		paper_GAP=(2*N(kc)+W(kc))/3;
	else
		paper_GAP=(N(kc)+W(kc))>>1;

	if(diff2>32)
		paper_GAP+=diff3>>2;
	else if(diff2>16)
		paper_GAP+=diff3*3>>4;
	else if(diff2>=-16)
		paper_GAP+=diff3>>3;
	else if(diff2>=-32)
		paper_GAP+=diff3>>4;

	if(diff>80)
		calic_GAP=W(kc);
	else if(diff<-80)
		calic_GAP=N(kc);
	else if(diff>32)
		calic_GAP=(2*N(kc)+6*W(kc)+NE(kc)-NW(kc))>>3;		//c1	[1/4  3/4  1/8  -1/8].[N W NE NW]
	else if(diff>8)
		calic_GAP=(6*N(kc)+10*W(kc)+3*(NE(kc)-NW(kc)))>>4;	//c2	[3/8  5/8  3/16  -3/16]
	else if(diff<-32)
		calic_GAP=(6*N(kc)+2*W(kc)+NE(kc)-NW(kc))>>3;		//c3	[3/4  1/4  1/8  -1/8]
	else if(diff<-8)
		calic_GAP=(10*N(kc)+6*W(kc)+3*(NE(kc)-NW(kc)))>>4;	//c4	[5/8  3/8  3/16  -3/16]
	else
		calic_GAP=(((N(kc)+W(kc))<<1)+NE(kc)-NW(kc))>>2;	//c5	[1/2  1/2  1/4  -1/4]
	
	int geomean;
	long long aN=N(kc)+0x800000LL, aW=W(kc)+0x800000LL;
	aN=CLAMP(0, aN, 0xFFFFFF);
	aW=CLAMP(0, aW, 0xFFFFFF);
	geomean=(int)floor_sqrt(aN*aW)-0x800000;

	int j=-1;
#define T51_PRED(X) pr->preds[++j]=X;
	T51_PREDLIST_SELF(kc)
	switch(kc)
	{
	case 0:
		if(pr->nch>=2)
		{
			T51_PREDLIST_C0_FROM_C1
			if(pr->nch>=3)
			{
				T51_PREDLIST_C0_FROM_C2
				if(pr->nch>=4)
				{
					T51_PREDLIST_C0_FROM_C3
				}
			}
		}
		break;
	case 1:
		T51_PREDLIST_C1_FROM_C0
		if(pr->nch>=3)
		{
			T51_PREDLIST_C1_FROM_C2
			if(pr->nch>=4)
			{
				T51_PREDLIST_C1_FROM_C3
			}
		}
		break;
	case 2:
		T51_PREDLIST_C2_FROM_C0_C1
		if(pr->nch>=4)
		{
			T51_PREDLIST_C2_FROM_C3
		}
		break;
	case 3:
		T51_PREDLIST_C3_FROM_C0_C1_C2
		break;
	}
#undef  T51_PRED
	pr->npreds=j+1;
	
#ifdef RECT_BUFFER
	for(int k=0;k<pr->npreds;++k)
		pr->qpreds[k]=CLAMP(-pr->half[pr->kc], pr->preds[k], pr->half[pr->kc]-1)+pr->half[pr->kc];
#else
	for(int k=0;k<pr->npreds;++k)
	{
		//pr->preds[k]=CLAMP(-pr->half[pr->kc], pr->preds[k], pr->half[pr->kc]-1)>>pr->shiftq[kc];
		pr->preds[k]=CLAMP(-pr->half[pr->kc], pr->preds[k], pr->half[pr->kc]-1);
#ifdef ENABLE_SSE
		pr->qpreds[k]=quantize(pr->preds[k])+(pr->maxqlevels>>1);
		pr->sse_cells[k]=pr->sse+pr->maxqlevels*((size_t)T51_NPREDS*pr->kc+k)+pr->qpreds[k];
		long long sum=*pr->sse_cells[k]>>12;
		int count=(int)(*pr->sse_cells[k]&0xFFF);
		pr->preds[k]+=(int)(sum/(count+1LL));
#endif
		pr->qpreds[k]=pr->quantize[kc](pr->preds[k])+(pr->maxqlevels>>1);
		pr->stat_cells[k]=pr->stats+((pr->maxqlevels*((size_t)T51_NPREDS*pr->kc+k)+pr->qpreds[k])<<pr->maxdepth);
	}
	pr->treeidx=1;
#endif
	pr->bitidx=pr->depths[pr->kc]-1;
	pr->partial=0;
	pr->lglrden=9;
	PROF(GETCTX);
}
static int t51_predict(T51Ctx *pr)
{
	pr->p0=0;
	for(int k=0;k<pr->npreds;++k)
	{
#ifdef RECT_BUFFER
		pr->stat_cells[k]=pr->stats+((pr->maxdepth*((size_t)T51_NPREDS*pr->kc+k)+pr->bitidx)<<pr->maxdepth|pr->qpreds[k]);
		pr->p0+=(long long)pr->mixer[pr->maxdepth*(T51_NPREDS*pr->kc+k)+pr->bitidx]**pr->stat_cells[k];
#else
		pr->p0+=(long long)pr->mixer[pr->maxdepth*(T51_NPREDS*pr->kc+k)+pr->bitidx]*pr->stat_cells[k][pr->treeidx];
#endif
	}
	PROF(PRED_MEM);
	int negmask=-(pr->p0<0);
	pr->p0=floor_sqrt(llabs(pr->p0));
	pr->p0^=negmask;
	pr->p0-=negmask;
	pr->p0+=0x8000;
#ifdef ENABLE_SSE
	long long *cell=pr->p0_bias+pr->maxdepth*pr->kc+pr->bitidx;
	long long sum=*cell>>12;
	int count=*cell&0xFFF;
	pr->p0+=sum/(count+1LL);
#endif
	pr->p0=CLAMP(1, pr->p0, 0xFFFF);
	PROF(PRED_SQRT);
	return (int)pr->p0;
}
static void t51_update(T51Ctx *pr, int bit)
{
	int proberror=((!bit<<16)-(int)pr->p0);

#ifdef DYNAMIC_DEN
	unsigned long long prob=bit?0x10000-pr->p0:pr->p0;

	//pr->csize_actual[pr->kc<<5|pr->bitidx]-=log2((double)prob/0x10000);
	int est_csize=-log2_fix24(prob<<8);
	//prob*=prob;
	//prob*=prob;
	//unsigned long long prob_hi=0, prob_lo;
	//prob_lo=_mul128(prob, prob, &prob_hi);
	//int est_csize=prob_hi?64-floor_log2(prob_hi):128-floor_log2(prob_lo);//[1~128]/8, fixed 5.3 bit

	//int est_csize=64-floor_log2(prob);//upper bound for csize, [1~64]/4, fixed 5.2 bit
	long long *csize=pr->csize_approx+(pr->kc<<5|pr->bitidx);
	*csize+=est_csize;
	pr->lglrden=*csize/((long long)pr->iw*pr->ky+pr->kx+1)>>(24-4);
	pr->lglrden=CLAMP(7, pr->lglrden, 9);
	//UPDATE_MIN(pr->lglrden, 8);
	int update=proberror>>pr->lglrden;
#else
	int update=proberror>>8;
#endif
	for(int k=0;k<pr->npreds;++k)
	{
		int *w=pr->mixer+pr->maxdepth*(T51_NPREDS*pr->kc+k)+pr->bitidx;
#ifdef RECT_BUFFER
		short *c=pr->stat_cells[k];
		int idx=pr->depths[pr->kc]-1-pr->bitidx;
		pr->qpreds[k]&=~(1<<idx);
		pr->qpreds[k]|=bit<<idx;
#else
		short *c=pr->stat_cells[k]+pr->treeidx;
#endif

		long long delta=(long long)update**c>>16;
		delta+=*w;
		delta=CLAMP(-0x8000, delta, 0x7FFF);
		*w=(int)delta;

		delta=(long long)*c+update;
		delta=CLAMP(-0x8000, delta, 0x7FFF);
		*c=(short)delta;
	}
#ifdef ENABLE_SSE
	long long *cell=pr->p0_bias+pr->maxdepth*pr->kc+pr->bitidx;
	long long sum=*cell>>12;
	int count=*cell&0xFFF;
	sum+=proberror;
	++count;
	if(count>=0x400)
	{
		count>>=1;
		sum>>=1;
	}
	*cell=sum<<12|count;
#endif

	pr->partial|=bit<<pr->bitidx;
	if(!pr->bitidx)
	{
		LOAD(pr->pixels, 0, 0)=pr->partial;
#ifdef ENABLE_SSE
		for(int k=0;k<pr->npreds;++k)
		{
			int error=pr->partial-pr->preds[k];
			sum=*pr->sse_cells[k]>>12;
			count=(int)(*pr->sse_cells[k]&0xFFF);
			sum+=error;
			++count;
			if(count>0x400)
			{
				count>>=1;
				sum>>=1;
			}
			*pr->sse_cells[k]=sum<<12|count;
		}
#endif
	}
	//pr->lglrden+=(pr->bitidx&2)!=0;
#ifndef RECT_BUFFER
	pr->treeidx=pr->treeidx<<1|bit;
#endif
	--pr->bitidx;
	PROF(UPDATE);
}
#undef  LOAD
#undef  NNNWWWW
#undef  NNNWWW
#undef  NNNWW
#undef  NNNW
#undef  NNN
#undef  NNNE
#undef  NNNEE
#undef  NNNEEE
#undef  NNNEEEE
#undef  NNWWWW
#undef  NNWWW
#undef  NNWW
#undef  NNW
#undef  NN
#undef  NNE
#undef  NNEE
#undef  NNEEE
#undef  NNEEEE
#undef  NWWWW
#undef  NWWW
#undef  NWW
#undef  NW
#undef  N
#undef  NE
#undef  NEE
#undef  NEEE
#undef  NEEEE
#undef  WWWW
#undef  WWW
#undef  WW
#undef  W
#ifdef T51_ENABLE_GUIDE
static Image *guide=0;
#endif
int t51_encode(Image const *src, ArrayHandle *data, int loud)
{
#ifdef T51_ENABLE_GUIDE
	//guide=src;
#endif
	double t_start=time_sec();
	int nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);//TODO
	UPDATE_MIN(nch, src->nch);
	if(loud)
	{
		int maxdepth=calc_maxdepth(src, 0);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T51  Enc %s  CWHD %d*%d*%d*%d/8\n", g_buf, nch, src->iw, src->ih, maxdepth);
	}
	char depths[4]=
	{
#ifdef ENABLE_RCT_MA
		src->depth[1],
		src->depth[2],
		src->depth[0],
		src->depth[3],
#else
		src->depth[1],
		src->depth[2]+1,
		src->depth[0]+1,
		src->depth[3],
#endif
	};
	double csizes[17*4]={0};
	DList list;
	ArithmeticCoder ec;
	T51Ctx pr;
	dlist_init(&list, 1, 1024, 0);
	ac_enc_init(&ec, &list);
	ptrdiff_t memusage=t51_init(&pr, src->iw, src->ih, nch, depths);
	Image *im2=0;
	image_copy(&im2, src);
	if(!memusage||!im2)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	rct_JPEG2000_32(im2, 1);
	pred_clampgrad(im2, 1, depths);
#ifdef T51_ENABLE_GUIDE
	free(guide);
	image_copy(&guide, im2);
#endif
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		t51_nextrow(&pr, ky);
		for(int kx=0;kx<src->iw;++kx, ++idx)
		{
			for(int kc=0;kc<nch;++kc)
			{
				//if(kc==0&&kx==1&&ky==0)//
				//if(kc==0&&kx==10&&ky==1)//
				//if(kc==1&&ky==0&&kx==1)//
				//if(kc==0&&kx==256&&ky==256)//
				//	printf("");

				t51_getctx(&pr, kc, kx);
				for(int kb=depths[kc]-1;kb>=0;--kb)
				{
					int p0=t51_predict(&pr);

					int bit=im2->data[idx<<2|kc]>>kb&1;
					ac_enc_bin(&ec, p0, bit);

					if(loud)
					{
						double prob=(double)(bit?0x10000-p0:p0)/0x10000;
						double bitsize=-log2(prob);
						csizes[17*kc+kb]+=bitsize;
						PROF(CSIZE);
					}

					t51_update(&pr, bit);
				}
			}
		}
	}
	ac_enc_flush(&ec);
	dlist_appendtoarray(&list, data);
	if(loud)
	{
		double usize=image_getBMPsize(src);
		printf("\n");
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		
		printf("csize,invCR:");
		for(int kc=0;kc<nch;++kc)
			printf("\t\tC%d%c", kc, kc<nch-1?'\t':'\n');
		int res=src->iw*src->ih;
		for(int kb=pr.maxdepth-1;kb>=0;--kb)
		{
			printf("B%02d ", kb);
			for(int kc=0;kc<nch;++kc)
			{
				double csize=csizes[17*kc+kb];
				printf("%15.0lf%8.3lf%%", csize/8, 100.*csize/res);
			}
			printf("\n");
		}
		printf("\n");

		printf("csize %8d  invCR %10.6lf%%  used %lf MB\n",
			(int)list.nobj, 100.*list.nobj/usize, (double)memusage/(1024*1024)
		);

		prof_print();
	}
	t51_free(&pr);
	dlist_clear(&list);
	free(im2);
	return 1;
}
int t51_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	double t_start=time_sec();
	int dst_nch=(dst->depth[0]!=0)+(dst->depth[1]!=0)+(dst->depth[2]!=0)+(dst->depth[3]!=0);
	int nch=MINVAR(dst_nch, dst->nch);
	char depths[4]=
	{
#ifdef ENABLE_RCT_MA
		dst->depth[1],
		dst->depth[2],
		dst->depth[0],
		dst->depth[3],
#else
		dst->depth[1],
		dst->depth[2]+1,
		dst->depth[0]+1,
		dst->depth[3],
#endif
	};
	ArithmeticCoder ec;
	T51Ctx pr;
	ac_dec_init(&ec, data, data+srclen);
	ptrdiff_t memusage=t51_init(&pr, dst->iw, dst->ih, nch, depths);
	if(!memusage)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int ky=0, idx=0;ky<dst->ih;++ky)
	{
		t51_nextrow(&pr, ky);
		for(int kx=0;kx<dst->iw;++kx, ++idx)
		{
			for(int kc=0;kc<nch;++kc)
			{
				int depth=depths[kc];
				//if(kc==1&&ky==0&&kx==1)//
				//	printf("");

				int sp=0;
				t51_getctx(&pr, kc, kx);
				for(int kb=depth-1;kb>=0;--kb)
				{
					int p0=t51_predict(&pr);

					int bit=ac_dec_bin(&ec, p0);
					sp|=bit<<kb;

					t51_update(&pr, bit);
				}
				sp<<=32-depth;//set sign
				sp>>=32-depth;
#ifdef T51_ENABLE_GUIDE
				if(guide&&sp!=guide->data[idx<<2|kc])
					LOG_ERROR("Guide error CXY %d %d %d  0x%08X!=0x%08X", kc, kx, ky, sp, guide->data[idx<<2|kc]);
#endif
				dst->data[idx<<2|kc]=sp;
			}
		}
	}
	pred_clampgrad(dst, 0, depths);
	rct_JPEG2000_32(dst, 0);
	if(loud)
	{
		printf("\n");
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	t51_free(&pr);
	return 1;
}

static void enc_bin(ArithmeticCoder *ec, int x, int depth, unsigned short *stats)
{
	int treeidx=1;
	for(int kb=depth-1;kb>=0;--kb)
	{
		int p0=stats[treeidx];

		int bit=x>>kb&1;
		ac_enc_bin(ec, p0, bit);

		p0+=((!bit<<16)-p0)>>5;
		p0=CLAMP(1, p0, 0xFFFF);
		stats[treeidx]=p0;
		treeidx=treeidx<<1|bit;
	}
}
#define SLIC5_CONFIG_EXP 5
#define SLIC5_CONFIG_MSB 2
#define SLIC5_CONFIG_LSB 0
typedef struct HybridUintStruct
{
	unsigned short token, nbits;
	unsigned bypass;
} HybridUint;
static void test_quantize(unsigned val, HybridUint *hu)
{
	int token, bypass, nbits;
	if(val<(1<<SLIC5_CONFIG_EXP))
	{
		token=val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=floor_log2_32((unsigned)val);
		int mantissa=val-(1<<lgv);
		token = (1<<SLIC5_CONFIG_EXP) + (
				(lgv-SLIC5_CONFIG_EXP)<<(SLIC5_CONFIG_MSB+SLIC5_CONFIG_LSB)|
				(mantissa>>(lgv-SLIC5_CONFIG_MSB))<<SLIC5_CONFIG_LSB|
				(mantissa&((1<<SLIC5_CONFIG_LSB)-1))
			);
		nbits=lgv-(SLIC5_CONFIG_MSB+SLIC5_CONFIG_LSB);
		bypass=val>>SLIC5_CONFIG_LSB&((1LL<<nbits)-1);
	}
	hu->token=token;
	hu->bypass=bypass;
	hu->nbits=nbits;
}
#if 0
static void enc_alpha(ArithmeticCoder *ec, int x, unsigned short *hist, int *hcount, unsigned *CDF, int cdfsize)
{
	x=x<<1^-(x<0);//pack sign
	HybridUint hu={0};
	test_quantize(x, &hu);
	ac_enc(ec, hu.token, CDF, cdfsize, 0);
	if(hu.nbits)
	{
		int bypass=hu.bypass, nbits=hu.nbits;
		while(nbits>8)
		{
			ac_enc(ec, bypass>>(nbits-8)&0xFF, 0, 1<<8, 16-8);
			nbits-=8;
		}
		ac_enc(ec, bypass&((1<<nbits)-1), 0, 1<<nbits, 16-nbits);
	}
	++*hcount;
	++hist[hu.token];
	if(*hcount>cdfsize&&!(*hcount&63))//update CDF
	{
		int sum=*hcount;
		long long c=0;
		for(int ks=0;ks<cdfsize;++ks)
		{
			long long freq=hist[ks];
			CDF[ks]=(int)(c*((1LL<<16)-cdfsize)/sum)+ks;
			c+=freq;
		}
		CDF[cdfsize]=1<<16;
	}
	if(*hcount>=0x1800)//rescale hist
	{
		int sum=0;
		for(int k=0;k<cdfsize;++k)
		{
			hist[k]=(hist[k]+1)>>1;
			sum+=hist[k];
		}
		*hcount=sum;
	}
}
#endif
void test_alphaVSbin(Image const *src)
{
	printf("\nAlphabetic vs Binary Coding\n");
	int inflation[]={0, 1, 1, 0};
	int maxdepth=calc_maxdepth(src, inflation);
	char depths[]=
	{
		src->depth[1],
		src->depth[2]+1,
		src->depth[0]+1,
		src->depth[3],
	};
	Image *im2=0;
	double elapsed=time_sec();
	image_copy(&im2, src);
	if(!im2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	rct_JPEG2000_32(im2, 1);
	pred_clampgrad(im2, 1, depths);
	elapsed=time_sec()-elapsed;
	printf("Decorrelation %lf sec\n", elapsed);

	//alphabetic coding
	{
		elapsed=time_sec();
		HybridUint hu;
		test_quantize(1<<maxdepth, &hu);
		int cdfsize=hu.token+1;
		unsigned short *hist=(unsigned short*)malloc(cdfsize*sizeof(short[3]));
		unsigned *CDF=(unsigned*)malloc((cdfsize+1LL)*sizeof(int[3]));
		if(!hist||!CDF)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		memset(hist, 0, cdfsize*sizeof(short[3]));
		int sum=0;
		for(int ks=0;ks<cdfsize;++ks)
		{
			int freq=0x10000/(ks+1);
			CDF[ks]=freq;
			sum+=freq;
			//freq=(int)((long long)freq*freq>>16);
		}
		for(int ks=0, c=0;ks<cdfsize;++ks)
		{
			int freq=CDF[ks];
			CDF[ks+2*(cdfsize+1)]=CDF[ks+cdfsize+1]=CDF[ks]=(int)((long long)c*((1LL<<16)-cdfsize)/sum)+ks;
			c+=freq;
		}
		CDF[cdfsize]=1<<16;
		
		ArithmeticCoder ec;
		DList list;
		dlist_init(&list, 1, 0x10000, 0);
		ac_enc_init(&ec, &list);
		int hcount[3]={0};
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				for(int kc=0;kc<3;++kc)
				{
					unsigned short *curr_hist=hist+kc*cdfsize;
					unsigned *curr_CDF=CDF+kc*(cdfsize+1);
					int x=im2->data[idx<<2|kc];
					x=x<<1^-(x<0);//pack sign
					HybridUint hu={0};
					test_quantize(x, &hu);
					ac_enc(&ec, hu.token, curr_CDF, cdfsize, 0);
					if(hu.nbits)
					{
						int bypass=hu.bypass, nbits=hu.nbits;
						while(nbits>8)
						{
							ac_enc(&ec, bypass>>(nbits-8)&0xFF, 0, 1<<8, 16-8);
							nbits-=8;
						}
						ac_enc(&ec, bypass&((1<<nbits)-1), 0, 1<<nbits, 16-nbits);
					}
					++hcount[kc];
					++curr_hist[hu.token];
					if(hcount[kc]>cdfsize&&!(hcount[kc]&63))//update CDF
					{
						int sum=hcount[kc];
						long long c=0;
						for(int ks=0;ks<cdfsize;++ks)
						{
							long long freq=curr_hist[ks];
							curr_CDF[ks]=(int)(c*((1LL<<16)-cdfsize)/sum)+ks;
							c+=freq;
						}
						curr_CDF[cdfsize]=1<<16;
					}
					if(hcount[kc]>=0x1800)//rescale hist
					{
						int sum=0;
						for(int k=0;k<cdfsize;++k)
						{
							curr_hist[k]=(curr_hist[k]+1)>>1;
							sum+=curr_hist[k];
						}
						hcount[kc]=sum;
					}
				}
				//int r=im2->data[idx<<2|0], g=im2->data[idx<<2|1], b=im2->data[idx<<2|2];
				//enc_alpha(&ec, r, hist+0*cdfsize, hcount+0, CDF+0*(cdfsize+1), cdfsize);
				//enc_alpha(&ec, g, hist+1*cdfsize, hcount+1, CDF+1*(cdfsize+1), cdfsize);
				//enc_alpha(&ec, b, hist+2*cdfsize, hcount+2, CDF+2*(cdfsize+1), cdfsize);
			}
		}
		ac_enc_flush(&ec);
		elapsed=time_sec()-elapsed;

		double usize=image_getBMPsize(src);
		printf("ALPHA %10lf sec  %.0lf -> %zd  invCR %6.2lf%%\n", elapsed, usize, list.nobj, 100.*list.nobj/usize);
		dlist_clear(&list);
		free(hist);
		free(CDF);
	}

	//binary coding
	{
		elapsed=time_sec();
		unsigned short *stats=(unsigned short*)malloc(sizeof(short[3])<<maxdepth);
		if(!stats)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		unsigned short fillval=0x8000;
		memfill(stats, &fillval, sizeof(short[3])<<maxdepth, sizeof(short));
		//memset(stats, 0, sizeof(short[3])<<maxdepth);

		ArithmeticCoder ec;
		DList list;
		dlist_init(&list, 1, 0x10000, 0);
		ac_enc_init(&ec, &list);
		//int prevW[3]={0};
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				int r=im2->data[idx<<2|0], g=im2->data[idx<<2|1], b=im2->data[idx<<2|2];
				enc_bin(&ec, r, depths[0], stats+(0LL<<maxdepth));
				enc_bin(&ec, g, depths[1], stats+(1LL<<maxdepth));
				enc_bin(&ec, b, depths[2], stats+(2LL<<maxdepth));
			}
		}
		ac_enc_flush(&ec);
		elapsed=time_sec()-elapsed;

		double usize=image_getBMPsize(src);
		printf("BIN   %10lf sec  %.0lf -> %zd  invCR %6.2lf%%\n", elapsed, usize, list.nobj, 100.*list.nobj/usize);
		dlist_clear(&list);
		free(stats);
	}
	free(im2);
}