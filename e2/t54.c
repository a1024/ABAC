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

#if 0
typedef struct T54CtxStruct
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
} T54Ctx;
int t54_encode(Image const *src, ArrayHandle *data, int loud)
{
#ifdef ENABLE_GUIDE
	//guide=src;
#endif
	double t_start=time_sec();
	int nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);//TODO
	UPDATE_MIN(nch, src->nch);
	if(loud)
	{
		int maxdepth=calc_maxdepth(src, 0);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T54 Bitplanes  Enc %s  CWHD %d*%d*%d*%d/8\n", g_buf, nch, src->iw, src->ih, maxdepth);
	}
	char depths[4]=
	{
		src->depth[1],
		src->depth[2]+1,
		src->depth[0]+1,
		src->depth[3],
	};
	double csizes[17*4]={0};
	DList list;
	ArithmeticCoder ec;
	T54Ctx pr;
	dlist_init(&list, 1, 1024, 0);
	ac_enc_init(&ec, &list);
	ptrdiff_t memusage=t54_init(&pr, src->iw, src->ih, nch, depths, &ec);
	Image *im2=0;
	image_copy(&im2, src);
	if(!memusage||!im2)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	rct_JPEG2000_32(im2, 1);
	pred_clampgrad(im2, 1, depths);
	for(int kc=0;kc<nch;++kc)
	{
		int depth=depths[kc];
		for(int kb=depth-1;kb>=0;--kb)
		{
			for(int ky=0, idx=0;ky<src->ih;++ky)
			{
				for(int kx=0;kx<src->iw;++kx, ++idx)
				{
					int bit=im2->data[idx<<2|kc]>>kb&1;
				}
			}
		}
	}
}
int t54_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
}
#endif