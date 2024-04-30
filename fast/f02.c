#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
//#define EC_USE_ARRAY
#include"ac.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE

//	#define USE_ANS
//	#define USE_GOLOMB
//	#define USE_ABAC
	#define USE_STATIC_PROB
	#define USE_CDF2SYM


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
		for(int c=0;c<0x10000;++c)
		{
			ks+=c>=CDF[ks+1];
			CDF2sym[c]=ks;
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
	for(ptrdiff_t ky=step<0?(ptrdiff_t)image->ih-1:0, cstep=(ptrdiff_t)image->nch*step, idx=cstep<0?((size_t)image->iw*image->ih-1LL)*image->nch:0;(size_t)ky<(size_t)image->ih;ky+=step)
	{
		for(ptrdiff_t kx=step<0?(ptrdiff_t)image->iw-1:0;(size_t)kx<(size_t)image->iw;kx+=step, idx+=cstep)
		{
			const short *comp=image->data+idx;
			for(int kc=step<0?image->nch-1:0;(unsigned)kc<(unsigned)image->nch;kc+=step)
			{
				//if(ky==50&&kx==395&&kc==2)
				//if(ky==2&&kx==506&&kc==0)
				//	printf("");
				if(fwd)
				{
					int bypass=comp[kc]<<1^-(comp[kc]<0);
#ifdef USE_GOLOMB
					gr_enc(&ec, bypass, 1<<image->depth>>1);
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
					int bypass=gr_dec(&ec, 1<<image->depth>>1);
#elif defined USE_ANS
					int bypass=ans_dec_POT(&ec, CDF, image->depth);
					//int bypass=ans_dec(&ec, 0, 1<<image->depth);
#elif defined USE_ABAC
					int bypass=0;
					for(int kb=image->depth-1;kb>=0;--kb)
						bypass|=ac_dec_bin(&ec, 0x8000)<<kb;
#else
#ifdef USE_CDF2SYM
					int bypass=ac_dec_CDF2sym(&ec, CDF, CDF2sym, 8);
#else
					int bypass=ac_dec_POT_permuted(&ec, pCDF, CDF, 8);
					//int bypass=ac_dec_POT(&ec, CDF, 8);
					//int bypass=ac_dec(&ec, CDF, 256);
#endif
					//int bypass=ac_dec_bypass(&ec, 1<<image->depth);

					//int bypass=0;
					//while(nbits>8)
					//{
					//	nbits-=8;
					//	bypass|=ac_dec(&ec, 0, 1<<8, 16-8)<<nbits;
					//}
					//bypass|=ac_dec(&ec, 0, 1<<nbits, 16-nbits);
#endif
					dst->data[idx+kc]=bypass>>1^-(bypass&1);
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
		dlist_appendtoarray(&list, data);
#endif
	}
	if(loud)
	{
		t0=time_sec()-t0;
		if(fwd)
		{
			ptrdiff_t usize=image_getBMPsize(image);
#ifdef EC_USE_ARRAY
			ptrdiff_t csize=ec.srcptr-ec.srcstart;
#else
			ptrdiff_t csize=list.nobj;
#endif
			printf("csize %12lld  %10.6lf%%  CR %8.6lf\n", csize, 100.*csize/usize, (double)usize/csize);
		}
		printf("F02-%s  %c %15.6lf sec\n", ecname, 'D'+fwd, t0);
	}
	dlist_clear(&list);
#ifdef USE_STATIC_PROB
	free(CDF);
#endif
	return 1;
}