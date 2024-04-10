#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"ac.h"
static const char file[]=__FILE__;

	#define ESTIMATE_CSIZES
//	#define ENABLE_GUIDE

//#ifdef ENABLE_GUIDE
//static const Image *guide=0;
//#endif
static int quantize(int val, int clevels)
{
	int negmask=-(val<0);
	val=abs(val);
	val=floor_log2_32(val)+1;
	val>>=1;
	val^=negmask;
	val-=negmask;
	val+=clevels>>1;
	return val;
}
int f04_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
//#ifdef ENABLE_GUIDE
//	if(fwd)
//		guide=image;
//#endif
	char depths[4]={0};
	memfill(depths, &image->depth, sizeof(char[4]), sizeof(char));
	if(image->nch>=3)
	{
		++depths[1];
		++depths[2];
	}
	int depth=MINVAR(image->depth+(image->nch>=3), 16), nlevels=1<<depth, half=nlevels>>1;
	int clevels=quantize(1<<depth, 0)<<1|1;
	int nctx=clevels;
	for(int k=1;k<image->nch;++k)
		nctx*=clevels;
	int cdfsize=1<<image->nch, CDFoffset=cdfsize+1, ctxsize=CDFoffset<<1;
	Image im2={0};
	image_copy_nodata(&im2, image);
	unsigned *stats=(unsigned*)malloc(sizeof(int)*ctxsize*nctx);//hist & CDF
	if(!im2.data||!stats)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	if(fwd)
	{
		memcpy(im2.data, image->data, sizeof(short)*im2.iw*im2.ih*im2.nch);
		rct_JPEG2000_32(&im2, 1);
		pred_clampgrad(&im2, 1, depths);
	}
	else
		memset(im2.data, 0, sizeof(short)*image->iw*image->ih*image->nch);
#ifdef ESTIMATE_CSIZES
	double csizes[16]={0};
#endif
	DList list;
	dlist_init(&list, 1, 256, 0);
	ArithmeticCoder ec;
	if(fwd)
		ac_enc_init(&ec, &list);
	else
		ac_dec_init(&ec, cbuf, cbuf+clen);
	int dx=image->nch;
	int dy=image->nch*image->iw;
	for(int k=0;k<CDFoffset;++k)
	{
		stats[k]=1;
		stats[CDFoffset+k]=(k<<16)/cdfsize;
	}
	stats[cdfsize]=cdfsize;//histsum
	memfill(stats+ctxsize, stats, sizeof(int)*ctxsize*(nctx-1LL), sizeof(int)*ctxsize);
	for(int kb=depth-1;kb>=0;--kb)
	{
		int MSBoffset=half&-(kb==depth-1);
		unsigned
			pmaska=~0U<< kb   , pmaskb=~pmaska>>1,
			fmaska=~0U<<(kb+1), fmaskb=~fmaska>>1;
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				int ctx=0, sym=0;
				for(int kc=0;kc<image->nch;++kc)
				{
#define LOAD(MA, MB, X, Y) ((unsigned)(ky+(Y))<(unsigned)image->ih&&(unsigned)(kx+(X))<(unsigned)image->iw?(im2.data[idx+(Y)*dy+(X)*dx+kc]+half)&MA:0)|MB
					int
						NW	=LOAD(pmaska, pmaskb, -1, -1),
						N	=LOAD(pmaska, pmaskb,  0, -1),
						NE	=LOAD(pmaska, pmaskb,  1, -1),
						W	=LOAD(pmaska, pmaskb, -1,  0),
						curr	=LOAD(fmaska, fmaskb,  0,  0),
						E	=LOAD(fmaska, fmaskb,  1,  0),
						SW	=LOAD(fmaska, fmaskb, -1,  1),
						S	=LOAD(fmaska, fmaskb,  0,  1),
						SE	=LOAD(fmaska, fmaskb,  1,  1);
#undef  LOAD
					int predNW=N+W-NW;	predNW=MEDIAN3(N, W, predNW);
					int predSW=S+W-SW;	predSW=MEDIAN3(S, W, predSW);
					int predSE=S+E-SE;	predSE=MEDIAN3(S, E, predSE);
					int predNE=N+E-NE;	predNE=MEDIAN3(N, E, predNE);
					int pred=(8*curr-(predNW+predSW+predSE+predNE))>>2;
					pred=CLAMP(0, pred, nlevels-1);
					ctx=clevels*ctx+quantize(pred, clevels);

					if(fwd)
						sym|=((im2.data[idx+kc]+half)>>kb&1)<<kc;
				}

				unsigned *curr_hist=stats+ctxsize*ctx, *curr_CDF=curr_hist+CDFoffset;
				if(fwd)
				{
					ac_enc(&ec, sym, curr_CDF, cdfsize);
#ifdef ESTIMATE_CSIZES
					csizes[kb]-=log2((double)(curr_CDF[sym+1]-curr_CDF[sym])/0x10000);
#endif
				}
				else
					sym=ac_dec(&ec, curr_CDF, cdfsize);

				++curr_hist[sym];
				++curr_hist[cdfsize];
				if(curr_hist[cdfsize]>=640)
				{
					int sum=0;
					for(int ks=0;ks<cdfsize;++ks)
						sum+=curr_hist[ks]=(curr_hist[ks]+1)>>1;
					curr_hist[cdfsize]=sum;
				}
				if(!(idx&63))
				{
					int sum=curr_hist[cdfsize];
					for(int ks=0, c=0;ks<cdfsize;++ks)
					{
						long long freq=curr_hist[ks];
						curr_CDF[ks]=(int)(c*(0x10000LL-cdfsize)/sum)+ks;
						c+=(int)freq;
					}
				}
				if(!fwd)
				{
					for(int kc=0;kc<image->nch;++kc)
					{
						int bit=sym>>kc&1;
						int val=im2.data[idx+kc];
						val|=bit<<kb;
						val-=MSBoffset;
						im2.data[idx+kc]=val;
//#ifdef ENABLE_GUIDE
//						if(guide&&(image->data[idx+kc]>>kb&1)!=(guide->data[idx+kc]>>kb&1))
//							LOG_ERROR("Guide error CXYB %d %d %d %d", kc, kx, ky, kb);
//#endif
					}
				}
			}
		}
	}
	if(fwd)
	{
		ac_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
	}
	else
	{
		memcpy(dst->data, im2.data, sizeof(short)*im2.iw*im2.ih*im2.nch);
		pred_clampgrad(dst, 0, depths);
		rct_JPEG2000_32(dst, 0);
	}
	if(loud)
	{
		t0=time_sec()-t0;
		if(fwd)
		{
			ptrdiff_t usize=image_getBMPsize(image);
			ptrdiff_t csize=list.nobj;
#ifdef ESTIMATE_CSIZES
			double uchsize=(double)image->iw*image->ih*3/8;
			double ctotal=0;
			printf("Bx        %12.3lf\n", uchsize);
			for(int kb=0;kb<depth;++kb)
			{
				double csize=csizes[kb]/8.;
				printf("B%d        %12.3lf  %15.6lf%%  %9.6lf\n", kb, csize, 100.*csize/uchsize, uchsize/csize);
				ctotal+=csizes[kb];
			}
			uchsize*=image->depth;
			ctotal/=8;
			printf("Estimated %12.3lf  %15.6lf%%  %9.6lf\n\n", ctotal, 100.*ctotal/uchsize, uchsize/ctotal);
#endif
			printf("csize %12lld      %15.6lf%%  %9.6lf  stats %lld bytes\n",
				csize, 100.*csize/usize, (double)usize/csize, sizeof(int)*ctxsize*nctx
			);
		}
		printf("%c %15.6lf sec\n", 'D'+fwd, t0);
	}
	dlist_clear(&list);
	image_clear(&im2);
	free(stats);
	return 1;
}