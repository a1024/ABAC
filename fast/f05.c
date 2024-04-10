#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
//#define EC_USE_ARRAY
#include"ac.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define PROFILER 0
//	#define USE_ANS
//	#define USE_GOLOMB


#ifdef PROFILER
#define CHECKPOINTLIST\
	CHECKPOINT(INIT)\
	CHECKPOINT(RCT)\
	CHECKPOINT(PRED)\
	CHECKPOINT(EC)\
	CHECKPOINT(HIST)\
	CHECKPOINT(CDF)\
	CHECKPOINT(FINISH)
#endif
#include"profiler.h"
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
static void update_CDF(unsigned *CDF, const unsigned short *hist, int nlevels)
{
	int sum=hist[nlevels];
	for(int ks=0, c=0;ks<nlevels;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*(0x10000LL-nlevels)/sum)+ks;
		//if(CDF[ks]>0x10000)
		//	LOG_ERROR("");
		c+=freq;
	}
	CDF[nlevels]=0x10000;
}
static void update_hist(unsigned short *hist, int nlevels, int sym)
{
	++hist[sym];
	++hist[nlevels];
	if(hist[nlevels]>8192)
	{
		int sum=0;
		for(int ks=0;ks<nlevels;++ks)
			sum+=hist[ks]=(hist[ks]+1)>>1;
		hist[nlevels]=sum;
	}
}
int f05_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	PROF_START();
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif
	DList list;
	dlist_init(&list, 1, 256, 0);
	ArithmeticCoder ec;
	if(image->depth>8)
		LOG_ERROR("Unsupported depth %d", image->depth);
	int nlevels=1<<image->depth, clevels=nlevels<<1, half=nlevels>>1, chalf=nlevels;
	short fillval=1;
	unsigned short *hist=(unsigned short*)malloc(sizeof(short[513*3]));
	unsigned *CDF=(unsigned*)malloc(sizeof(int[513*3]));
	if(!hist||!CDF)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memfill(hist, &fillval, sizeof(short[513*3]), sizeof(short));
	hist[513*0+nlevels]=nlevels;
	hist[513*1+clevels]=clevels;
	hist[513*2+clevels]=clevels;
	update_CDF(CDF+513*0, hist+513*0, 256);
	update_CDF(CDF+513*1, hist+513*1, 512);
	update_CDF(CDF+513*2, hist+513*2, 512);
	PROF(INIT);
	if(fwd)
	{
		ac_enc_init(&ec, &list);
		if(image->nch==3)
		{
			for(int ky=0, idx=0;ky<image->ih;++ky)
			{
				short W[3]={0};
				for(int kx=0;kx<image->iw;++kx, idx+=3)
				{
					//if(idx==2319)//
					//if(idx==2304)//
					//if(idx==0x06fbb434-3)//
					//	printf("");

					short curr[3], val[3];
					memcpy(curr, image->data+idx, sizeof(short[3]));
					curr[0]-=curr[1];
					curr[2]-=curr[1];
					curr[1]+=(curr[0]+curr[2])>>2;
					PROF(RCT);

					val[1]=(curr[1]-W[1]+half)&(nlevels-1);
					val[2]=(curr[2]-W[2]+chalf)&(clevels-1);
					val[0]=(curr[0]-W[0]+chalf)&(clevels-1);
					PROF(PRED);

					ac_enc(&ec, val[1], CDF+513*0, nlevels);
					ac_enc(&ec, val[2], CDF+513*1, clevels);
					ac_enc(&ec, val[0], CDF+513*2, clevels);
					PROF(EC);

					update_hist(hist+513*0, nlevels, val[1]);
					update_hist(hist+513*1, clevels, val[2]);
					update_hist(hist+513*2, clevels, val[0]);
					PROF(HIST);

					memcpy(W, curr, sizeof(short[3]));
				}
				update_CDF(CDF+513*0, hist+513*0, nlevels);
				update_CDF(CDF+513*1, hist+513*1, clevels);
				update_CDF(CDF+513*2, hist+513*2, clevels);
				PROF(CDF);
			}
		}
		else if(image->nch==1)
		{
			for(int ky=0, idx=0;ky<image->ih;++ky)
			{
				int W=0;
				for(int kx=0;kx<image->iw;++kx, ++idx)
				{
					int val=(image->data[idx]-W+half)&(nlevels-1);
					ac_enc(&ec, val, CDF+513*0, nlevels);
					update_hist(hist+513*0, nlevels, val);
				}
				update_CDF(CDF+513*0, hist+513*0, nlevels);
			}
		}
		else
			LOG_ERROR("Unsupported number of channels %d", image->nch);
		ac_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
		PROF(FINISH);
	}
	else
	{
		ac_dec_init(&ec, cbuf, cbuf+clen);
		if(image->nch==3)
		{
			for(int ky=0, idx=0;ky<image->ih;++ky)
			{
				short W[3]={0};
				for(int kx=0;kx<image->iw;++kx, idx+=3)
				{
					//if(idx==2319)//
					//if(idx==2304)//
					//if(idx==0x06fbb434-3)//
					//	printf("");

					short curr[3], val[3];

					val[1]=ac_dec(&ec, CDF+513*0, nlevels);
					val[2]=ac_dec(&ec, CDF+513*1, clevels);
					val[0]=ac_dec(&ec, CDF+513*2, clevels);
					PROF(EC);

					curr[1]=((val[1]+W[1])&(nlevels-1))-half;
					curr[2]=((val[2]+W[2])&(clevels-1))-chalf;
					curr[0]=((val[0]+W[0])&(clevels-1))-chalf;
					PROF(PRED);

					update_hist(hist+513*0, nlevels, val[1]);
					update_hist(hist+513*1, clevels, val[2]);
					update_hist(hist+513*2, clevels, val[0]);
					PROF(HIST);

					memcpy(W, curr, sizeof(short[3]));
					curr[1]-=(curr[0]+curr[2])>>2;
					curr[2]+=curr[1];
					curr[0]+=curr[1];
					memcpy(dst->data+idx, curr, sizeof(short[3]));
					PROF(RCT);
#ifdef ENABLE_GUIDE
					if(guide&&memcmp(curr, guide->data+idx, image->nch*sizeof(short)))
					{
						short comp0[4]={0};
						memcpy(comp0, guide->data+idx, image->nch*sizeof(short));
						curr[0]-=curr[1];
						curr[2]-=curr[1];
						curr[1]+=(curr[0]+curr[2])>>2;
						comp0[0]-=comp0[1];
						comp0[2]-=comp0[1];
						comp0[1]+=(comp0[0]+comp0[2])>>2;
						LOG_ERROR("Guide error XY %d %d", kx, ky);
						printf("");//
					}
#endif
				}
				update_CDF(CDF+513*0, hist+513*0, nlevels);
				update_CDF(CDF+513*1, hist+513*1, clevels);
				update_CDF(CDF+513*2, hist+513*2, clevels);
				PROF(CDF);
			}
		}
		else if(image->nch==1)
		{
			for(int ky=0, idx=0;ky<image->ih;++ky)
			{
				int W=0;
				for(int kx=0;kx<image->iw;++kx, ++idx)
				{
					int val=ac_dec(&ec, CDF+513*0, nlevels);
					dst->data[idx]=((val+W)&(nlevels-1))-half;
					update_hist(hist+513*0, nlevels, val);
				}
				update_CDF(CDF+513*0, hist+513*0, nlevels);
			}
		}
		else
			LOG_ERROR("Unsupported number of channels %d", image->nch);
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
		printf("F05  %c %15.6lf sec\n", 'D'+fwd, t0);
		prof_print();
	}
	dlist_clear(&list);
	free(hist);
	free(CDF);
	return 1;
}