#include"best.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
//#define EC_USE_ARRAY
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE

//	#define USE_ANS
//	#define ANS_DEC_FWD
	#define USE_GOLOMB
//	#define USE_ABAC
	#define USE_STATIC_PROB
	#define USE_CDF2SYM


//THIS CODEC WAS MODIFIED FOR  Image struct with int[4]
#define AC_IMPLEMENTATION
#include"ac.h"

#ifdef USE_ANS
static const char ecname[]="ANS";
#elif defined USE_GOLOMB
static const char ecname[]="Golomb-Rice";
#else
static const char ecname[]="AC";
#endif
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
int f02_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	const int pixelstride=4;
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif
#ifdef USE_STATIC_PROB
	unsigned *CDF=(unsigned*)malloc(sizeof(int[257]));
	unsigned *pCDF=(unsigned*)malloc(sizeof(int[257]));
#ifdef USE_CDF2SYM
	unsigned char *CDF2sym=(unsigned char*)malloc(sizeof(char[0x10000]));
#endif
	if(!CDF||!pCDF
#ifdef USE_CDF2SYM
		||!CDF2sym
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int k=0;k<257;++k)
		CDF[k]=k<<8;
	if(!fwd)
	{
		int nlevels=256;
		for(int step=nlevels, ks2=1;step>1;step>>=1)
		{
			for(int ks=0;ks<nlevels;ks+=step, ++ks2)
				pCDF[ks2]=CDF[ks|step>>1];
		}
		pCDF[0]=CDF[0];
		pCDF[nlevels]=CDF[nlevels];
#ifdef USE_CDF2SYM
		int ks=0;
		for(unsigned c=0;c<0x10000;++c)
		{
			ks+=c>=CDF[ks+1];
			CDF2sym[c]=(unsigned char)ks;
		}
#endif
	}
#endif
	DList list;
	dlist_init(&list, 1, 256, 0);
#ifdef USE_ANS
	ANSCoder ec;
	if(fwd)
		ans_enc_init(&ec, &list);
	else
		ans_dec_init(&ec, cbuf, cbuf+clen);
	int step=!fwd-fwd;
#elif defined USE_GOLOMB
	const int step=1;
	GolombRiceCoder ec;
#ifdef EC_USE_ARRAY
	size_t dststart=0;
#endif
	if(fwd)
#ifdef EC_USE_ARRAY
	{
		dststart=array_append(data, 0, 1, (ptrdiff_t)image->iw*image->ih*image->nch*sizeof(short[2]), 1, 0, 0);
		gr_enc_init(&ec, data[0]->data, data[0]->data+data[0]->count);
	}
#else
		gr_enc_init(&ec, &list);
#endif
	else
		gr_dec_init(&ec, cbuf, cbuf+clen);
#else
	ArithmeticCoder ec;
	if(fwd)
		ac_enc_init(&ec, &list);
	else
		ac_dec_init(&ec, cbuf, cbuf+clen);
	const int step=1;
#endif
	int bypass;
	for(ptrdiff_t ky=step<0?(ptrdiff_t)image->ih-1:0, cstep=(ptrdiff_t)pixelstride*step, idx=cstep<0?((size_t)image->iw*image->ih-1LL)*image->nch:0;(size_t)ky<(size_t)image->ih;ky+=step)
	{
		for(ptrdiff_t kx=step<0?(ptrdiff_t)image->iw-1:0;(size_t)kx<(size_t)image->iw;kx+=step, idx+=cstep)
		{
			const int *comp=image->data+idx;
			for(int kc=step<0?image->nch-1:0;(unsigned)kc<(unsigned)image->nch;kc+=step)
			{
				//if(ky==50&&kx==395&&kc==2)
				//if(ky==2&&kx==506&&kc==0)
				//if(idx==2316&&kc==2)//
				//	printf("");
				if(fwd)
				{
					bypass=comp[kc]<<1^-(comp[kc]<0);
#ifdef USE_GOLOMB
					gr_enc_POT(&ec, bypass, image->depth[kc]-2);
#elif defined USE_ANS
					ans_enc(&ec, bypass, CDF, 1<<image->depth);//up to 16 bit
#elif defined USE_ABAC
					for(int kb=image->depth-1;kb>=0;--kb)
						ac_enc_bin(&ec, 0x8000, bypass>>kb&1);
#else
					ac_enc(&ec, bypass, CDF);
					//ac_enc_bypass(&ec, bypass, 1<<image->depth);
					//while(nbits>8)
					//{
					//	ac_enc(&ec, bypass>>(nbits-8)&0xFF, 0, 1<<8, 16-8);
					//	nbits-=8;
					//}
					//ac_enc(&ec, bypass&((1<<nbits)-1), 0, 1<<nbits, 16-nbits);
#endif
				}
				else
				{
#ifdef USE_GOLOMB
					bypass=gr_dec_POT(&ec, image->depth[kc]-2);
#elif defined USE_ANS
#ifdef USE_CDF2SYM
					bypass=ans_dec_CDF2sym(&ec, CDF, CDF2sym);
#else
					bypass=ans_dec_POT(&ec, CDF, image->depth);
					//bypass=ans_dec(&ec, 0, 1<<image->depth);
#endif
#elif defined USE_ABAC
					bypass=0;
					for(int kb=image->depth-1;kb>=0;--kb)
						bypass|=ac_dec_bin(&ec, 0x8000)<<kb;
#else
#ifdef USE_CDF2SYM
					bypass=ac_dec_CDF2sym(&ec, CDF, CDF2sym);
#else
					bypass=ac_dec_POT_permuted(&ec, pCDF, CDF, 8);
					//bypass=ac_dec_POT(&ec, CDF, 8);
					//bypass=ac_dec(&ec, CDF, 256);
#endif
					//bypass=ac_dec_bypass(&ec, 1<<image->depth);

					//bypass=0;
					//while(nbits>8)
					//{
					//	nbits-=8;
					//	bypass|=ac_dec(&ec, 0, 1<<8, 16-8)<<nbits;
					//}
					//bypass|=ac_dec(&ec, 0, 1<<nbits, 16-nbits);
#endif
					dst->data[idx+kc]=(short)(bypass>>1^-(bypass&1));
				}
			}
#ifdef ENABLE_GUIDE
			if(!fwd)
			{
				if(guide&&memcmp(comp, guide->data+idx, image->nch*sizeof(short)))
				{
					short comp0[4]={0};
					memcpy(comp0, guide->data+idx, image->nch*sizeof(short));
					//for(int k=1;k<image->nch;++k)
					//{
					//	pixels[k]-=pixels[0];
					//	comp0[k]-=comp0[0];
					//}
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
			}
#endif
		}
	}
#ifdef USE_ANS
	volatile double tr=0;
#endif
	if(fwd)
	{
#ifdef USE_ANS
		ans_enc_flush(&ec);
#elif defined USE_GOLOMB
		gr_enc_flush(&ec);
#else
		ac_enc_flush(&ec);
#endif
#ifdef EC_USE_ARRAY
		data[0]->count=dststart+ec.srcptr-ec.srcstart;
#else
#if defined USE_ANS && defined ANS_DEC_FWD
		size_t dststart=dlist_appendtoarray(&list, data);
#if 0
		{
			unsigned short xdata[]=
			{
				 0,  1,  2,  3,  4,  5,  6,  7,
				 8,  9, 10, 11, 12, 13, 14, 15,
				16, 17, 18, 19, 20, 21, 22, 23,
				24, 25, 26, 27, 28, 29, 30, 31,
				32, 33, 34, 35, 36, 37, 38, 39,
				40,
			};
			reverse16(xdata, xdata+_countof(xdata));
			printf("idx  reversed\n");
			for(int k=0;k<_countof(xdata);++k)
				printf("%2d %2d\n", k, xdata[k]);
			LOG_ERROR("Test");
		}
#endif
		tr=time_sec();
		reverse16(data[0]->data+dststart, data[0]->data+data[0]->count);
		tr=time_sec()-tr;
#else
		dlist_appendtoarray(&list, data);
#endif
#endif
	}
	if(loud)
	{
		t0=time_sec()-t0;
		ptrdiff_t usize=image_getBMPsize(image);
		if(fwd)
		{
#ifdef EC_USE_ARRAY
			ptrdiff_t csize=ec.srcptr-ec.srcstart;
#else
			ptrdiff_t csize=list.nobj;
#endif
			printf("csize %12td  %10.6lf%%  CR %8.6lf\n", csize, 100.*csize/usize, (double)usize/csize);
#if defined USE_ANS && defined ANS_DEC_FWD
			printf("reverse %15.6lf sec  %15.6lf MB/s\n", tr, usize/(1024.*1024.*tr));
#endif
		}
		printf("F02-%s  %c %15.6lf sec  %15.6lf MB/s\n", ecname, 'D'+fwd, t0, usize/(1024.*1024.*t0));
	}
	dlist_clear(&list);
#ifdef USE_STATIC_PROB
	free(CDF);
#endif
	return 1;
}