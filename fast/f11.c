#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;

#include"ac.h"
#include"profiler.h"


static int CG(int N, int W, int NW)
{
	int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
	int pred=N+W-NW;
	pred=CLAMP(vmin, pred, vmax);
	return pred;
}
static int clamp(int x, int half)
{
	return CLAMP(-half, x, half-1);
}
#define NPREDS 3
#define NRCH (NPREDS*3)
#define NGCH NPREDS
#define NBCH (NPREDS*2)
#define CHLIST\
	CHANNEL(g)\
	CHANNEL(g-(gN+gW)/2)\
	CHANNEL(g-CG(gN, gW, gNW))\
	CHANNEL(b)\
	CHANNEL(b-(bN+bW)/2)\
	CHANNEL(b-CG(bN, bW, bNW))\
	CHANNEL(b-g)\
	CHANNEL(b-clamp((bN-gN+bW-gW)/2+g, half))\
	CHANNEL(b-clamp(CG(bN-gN, bW-gW, bNW-gNW)+g, half))\
	CHANNEL(r)\
	CHANNEL(r-(rN+rW)/2)\
	CHANNEL(r-CG(rN, rW, rNW))\
	CHANNEL(r-g)\
	CHANNEL(r-clamp((rN-gN+rW-gW)/2+g, half))\
	CHANNEL(r-clamp(CG(rN-gN, rW-gW, rNW-gNW)+g, half))\
	CHANNEL(r-b)\
	CHANNEL(r-clamp((rN-bN+rW-bW)/2+b, half))\
	CHANNEL(r-clamp(CG(rN-bN, rW-bW, rNW-bNW)+b, half))

//	CHANNEL(g-gN)
//	CHANNEL(g-gW)
//	CHANNEL(b-bN)
//	CHANNEL(b-bW)
//	CHANNEL(b-clamp(bN-gN+g, half))
//	CHANNEL(b-clamp(bW-gW+g, half))
//	CHANNEL(r-rN)
//	CHANNEL(r-rW)
//	CHANNEL(r-clamp(rN-gN+g, half))
//	CHANNEL(r-clamp(rW-gW+g, half))
//	CHANNEL(r-clamp(rN-bN+b, half))
//	CHANNEL(r-clamp(rW-bW+b, half))
static const char *chnames[]=
{
#define CHANNEL(X) #X,
	CHLIST
#undef  CHANNEL
};
int f11_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;

	if(image->nch!=3)
	{
		LOG_ERROR("Expected 3 channels, got %d", image->nch);
		return 0;
	}
	int depth=image->depth, nlevels=1<<depth, half=nlevels>>1;
	size_t histsize=sizeof(int[30])<<image->depth;
	int *hist=(int*)malloc(histsize);
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(hist, 0, histsize);
	int rowstride=image->iw*3;
	double t=time_sec();
	for(int ky=1;ky<image->ih;++ky)
	{
		for(int kx=1;kx<image->iw;++kx)
		{
#define LOAD(C, X, Y) image->data[idx-rowstride*Y-3*X+C]
			int
				idx=(image->iw*ky+kx)*3,
				r	=LOAD(0, 0, 0),
				rN	=LOAD(0, 0, 1),
				rW	=LOAD(0, 1, 0),
				rNW	=LOAD(0, 1, 1),
				g	=LOAD(1, 0, 0),
				gN	=LOAD(1, 0, 1),
				gW	=LOAD(1, 1, 0),
				gNW	=LOAD(1, 1, 1),
				b	=LOAD(2, 0, 0),
				bN	=LOAD(2, 0, 1),
				bW	=LOAD(2, 1, 0),
				bNW	=LOAD(2, 1, 1);
			int preds[]=
			{
#define CHANNEL(X) X,
				CHLIST
#undef  CHANNEL
			};
			for(int kt=0;kt<NRCH+NGCH+NBCH;++kt)
			{
				preds[kt]+=half;
				preds[kt]&=nlevels-1;
				++hist[kt<<depth|preds[kt]];
			}
		}
	}
	double t1=time_sec();
	ptrdiff_t field=(image->iw-1LL)*(image->ih-1LL);
	double csizes[NRCH+NGCH+NBCH]={0};
	for(int kt=0;kt<NRCH+NGCH+NBCH;++kt)
	{
		int *curr_hist=hist+((size_t)kt<<depth);
		for(int ks=0;ks<nlevels;++ks)
		{
			int freq=curr_hist[ks];
			if(freq)
				csizes[kt]-=freq*log2((double)freq/field);
		}
		csizes[kt]/=8;
	}
	double t2=time_sec();
	int ridx=0, gidx=0, bidx=0;
	for(int kt=1;kt<NGCH;++kt)
	{
		if(csizes[gidx]>csizes[kt])
			gidx=kt;
	}
	bidx=NGCH;
	for(int kt=NGCH+1;kt<NGCH+NBCH;++kt)
	{
		if(csizes[bidx]>csizes[kt])
			bidx=kt;
	}
	ridx=NGCH+NBCH;
	for(int kt=NGCH+NBCH+1;kt<NRCH+NGCH+NBCH;++kt)
	{
		if(csizes[ridx]>csizes[kt])
			ridx=kt;
	}
	ptrdiff_t usize=field*depth>>3;
	printf("U   %11lld\n", usize);
	for(int kt=0;kt<NRCH+NGCH+NBCH;++kt)
	{
		char c=' ';
		if(kt==ridx)
			c='r';
		if(kt==gidx)
			c='g';
		if(kt==bidx)
			c='b';
		printf("%2d  %14.2lf  %6.2lf%%  %14.6lf  %c  %s\n", kt, csizes[kt], 100.*csizes[kt]/usize, usize/csizes[kt], c, chnames[kt]);
	}
	printf("Sampling %14.6lf sec  %lf B/s\n", t1-t, usize/(t1-t));
	printf("Analysis %14.6lf sec\n", t2-t);

	free(hist);
	LOG_ERROR("This isn't a codec");
	return 1;
}