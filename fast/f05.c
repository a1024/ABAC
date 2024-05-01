#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
//#define EC_USE_ARRAY
#include"ac.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define PROFILER 1
	#define USE_CLAMPGRAD
	#define ALLOW_SIMD


#ifdef PROFILER
#define CHECKPOINTLIST\
	CHECKPOINT(INIT)\
	CHECKPOINT(RCT)\
	CHECKPOINT(PRED)\
	CHECKPOINT(DUMMY)\
	CHECKPOINT(EC)\
	CHECKPOINT(HIST)\
	CHECKPOINT(CDF)\
	CHECKPOINT(FINISH)
#endif
#include"profiler.h"
#ifdef ALLOW_SIMD
#include<tmmintrin.h>
#endif
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
		hist[nlevels]=(short)sum;
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
	DList list[3];
	dlist_init(list+0, 1, 256, 0);
	dlist_init(list+1, 1, 256, 0);
	dlist_init(list+2, 1, 256, 0);
	ArithmeticCoder ec[3];
	if(image->depth>8)
		LOG_ERROR("Unsupported depth %d", image->depth);
	int nlevels=1<<image->depth, clevels=nlevels<<1, half=nlevels>>1, chalf=nlevels;
	short fillval=1;
	unsigned short *hist=(unsigned short*)malloc(sizeof(short[513*3]));
	unsigned *CDF=(unsigned*)malloc(sizeof(int[513*3]));
#ifdef USE_CLAMPGRAD
	short *pixels=(short*)malloc((image->iw+4LL)*sizeof(short[2*4]));//2 padded rows * 4 channels max
#endif
	if(!hist||!CDF
#ifdef USE_CLAMPGRAD
		||!pixels
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memfill(hist, &fillval, sizeof(short[513*3]), sizeof(short));
	hist[513*0+nlevels]=(short)nlevels;
	hist[513*1+clevels]=(short)clevels;
	hist[513*2+clevels]=(short)clevels;
	update_CDF(CDF+513*0, hist+513*0, 256);
	update_CDF(CDF+513*1, hist+513*1, 512);
	update_CDF(CDF+513*2, hist+513*2, 512);
#ifdef USE_CLAMPGRAD
	memset(pixels, 0, (image->iw+4LL)*sizeof(short[2*4]));
#endif
	PROF(INIT);
	if(fwd)
	{
		if(image->nch==3)
		{
			ac_enc_init(ec+0, list+0);
			ac_enc_init(ec+1, list+1);
			ac_enc_init(ec+2, list+2);
			for(int ky=0, idx=0;ky<image->ih;++ky)
			{
#ifdef USE_CLAMPGRAD
				short *rows[2]=
				{
					pixels+(((image->iw+4LL)*(ky&1)+1)<<2),
					pixels+(((image->iw+4LL)*((ky-1)&1)+1)<<2),
				};
#else
				short W[3]={0};
#endif
				for(int kx=0;kx<image->iw;++kx, idx+=3)
				{
					//if(idx==2319)//
					//if(idx==2304)//
					//if(idx==0x06fbb434-3)//
					//	printf("");
#ifdef USE_CLAMPGRAD
					__m128i mNW	=_mm_loadu_si128((const __m128i*)(rows[1]-4));
					__m128i mN	=_mm_loadu_si128((const __m128i*)(rows[1]+0));
					__m128i mW	=_mm_loadu_si128((const __m128i*)(rows[0]-4));
					__m128i vmin=_mm_min_epi16(mN, mW);
					__m128i vmax=_mm_max_epi16(mN, mW);
					__m128i mpred=_mm_add_epi16(mN, mW);
					mpred=_mm_sub_epi16(mpred, mNW);
					mpred=_mm_min_epi16(mpred, vmax);
					mpred=_mm_max_epi16(mpred, vmin);
					ALIGN(16) short pred[8];
					_mm_store_si128((__m128i*)pred, mpred);

					//short
					//	NW[]={rows[1][-4], rows[1][-3], rows[1][-2]},
					//	N []={rows[1][ 0], rows[1][ 1], rows[1][ 2]},
					//	W []={rows[0][-4], rows[0][-3], rows[0][-2]};

					//short pred[]=
					//{
					//	N[0]+W[0]-NW[0],//Cr
					//	N[1]+W[1]-NW[1],//Y
					//	N[2]+W[2]-NW[2],//Cb
					//};
					//pred[0]=MEDIAN3(N[0], W[0], pred[0]);
					//pred[1]=MEDIAN3(N[1], W[1], pred[1]);
					//pred[2]=MEDIAN3(N[2], W[2], pred[2]);

					//short pred[3], vmin[3], vmax[3];	//X  slower
					//for(int k=0;k<3;++k)
					//{
					//	if(N[k]<W[k])
					//		vmin[k]=N[k], vmax[k]=W[k];
					//	else
					//		vmin[k]=W[k], vmax[k]=N[k];
					//	if((unsigned)(NW[k]-vmin[k])<(unsigned)(vmax[k]-vmin[k]))
					//		pred[k]=N[k]+W[k]-NW[k];
					//	else if(NW[k]<vmin[k])
					//		pred[k]=vmax[k];
					//	else
					//		pred[k]=vmin[k];
					//}
					
					short *curr=rows[0];
					rows[0]+=4;
					rows[1]+=4;
					PROF(PRED);
					PROF(DUMMY);

					memcpy(curr, image->data+idx, sizeof(short[3]));
					curr[0]-=curr[1];
					curr[2]-=curr[1];
					curr[1]+=(curr[0]+curr[2])>>2;

					short val[]=
					{
						(short)((curr[0]-pred[0]+chalf)&(clevels-1)),
						(short)((curr[1]-pred[1]+ half)&(nlevels-1)),
						(short)((curr[2]-pred[2]+chalf)&(clevels-1)),
					};
					PROF(RCT);
#else
					short curr[3], val[3];
					memcpy(curr, image->data+idx, sizeof(short[3]));
					curr[0]-=curr[1];
					curr[2]-=curr[1];
					curr[1]+=(curr[0]+curr[2])>>2;
					PROF(RCT);

					val[1]=(curr[1]-W[1]+half)&(nlevels-1);
					val[2]=(curr[2]-W[2]+chalf)&(clevels-1);
					val[0]=(curr[0]-W[0]+chalf)&(clevels-1);

					memcpy(W, curr, sizeof(short[3]));
					PROF(PRED);
#endif

					ac_enc(ec+0, val[1], CDF+513*0);
					ac_enc(ec+1, val[2], CDF+513*1);
					ac_enc(ec+2, val[0], CDF+513*2);
					PROF(EC);

					update_hist(hist+513*0, nlevels, val[1]);
					update_hist(hist+513*1, clevels, val[2]);
					update_hist(hist+513*2, clevels, val[0]);
					PROF(HIST);
				}
				update_CDF(CDF+513*0, hist+513*0, nlevels);
				update_CDF(CDF+513*1, hist+513*1, clevels);
				update_CDF(CDF+513*2, hist+513*2, clevels);
				PROF(CDF);
			}
			ac_enc_flush(ec+0);
			ac_enc_flush(ec+1);
			ac_enc_flush(ec+2);
			unsigned bm[]=
			{
				(unsigned)list[0].nobj,
				(unsigned)list[1].nobj,
				(unsigned)list[2].nobj,
			};
			array_append(data, bm, 1, sizeof(bm), 1, 0, 0);
			dlist_appendtoarray(list+0, data);
			dlist_appendtoarray(list+1, data);
			dlist_appendtoarray(list+2, data);
		}
		else if(image->nch==1)
		{
			ac_enc_init(ec+0, list+0);
			for(int ky=0, idx=0;ky<image->ih;++ky)
			{
				int W=0;
				for(int kx=0;kx<image->iw;++kx, ++idx)
				{
					int val=(image->data[idx]-W+half)&(nlevels-1);
					ac_enc(ec, val, CDF+513*0);
					update_hist(hist+513*0, nlevels, val);
				}
				update_CDF(CDF+513*0, hist+513*0, nlevels);
			}
			dlist_appendtoarray(list+0, data);
		}
		else
			LOG_ERROR("Unsupported number of channels %d", image->nch);
		PROF(FINISH);
	}
	else
	{
		if(image->nch==3)
		{
			unsigned bm[3];
			memcpy(bm, cbuf, sizeof(bm));
			cbuf+=sizeof(bm);
			clen-=sizeof(bm);
			const unsigned char *ptr=cbuf;
			ac_dec_init(ec+0, ptr, ptr+bm[0]);	ptr+=bm[0];
			ac_dec_init(ec+1, ptr, ptr+bm[1]);	ptr+=bm[1];
			ac_dec_init(ec+2, ptr, ptr+bm[2]);	ptr+=bm[2];
			for(int ky=0, idx=0;ky<image->ih;++ky)
			{
#ifdef USE_CLAMPGRAD
				short *rows[2]=
				{
					pixels+(((image->iw+4LL)*(ky&1)+1)<<2),
					pixels+(((image->iw+4LL)*((ky-1)&1)+1)<<2),
				};
#else
				short W[3]={0};
#endif
				for(int kx=0;kx<image->iw;++kx, idx+=3)
				{
					//if(idx==2319)//
					//if(idx==2304)//
					//if(idx==0x06fbb434-3)//
					//	printf("");
					
#ifdef USE_CLAMPGRAD
					__m128i mNW	=_mm_loadu_si128((const __m128i*)(rows[1]-4));
					__m128i mN	=_mm_loadu_si128((const __m128i*)(rows[1]+0));
					__m128i mW	=_mm_loadu_si128((const __m128i*)(rows[0]-4));
					__m128i vmin=_mm_min_epi16(mN, mW);
					__m128i vmax=_mm_max_epi16(mN, mW);
					__m128i mpred=_mm_add_epi16(mN, mW);
					mpred=_mm_sub_epi16(mpred, mNW);
					mpred=_mm_min_epi16(mpred, vmax);
					mpred=_mm_max_epi16(mpred, vmin);
					ALIGN(16) short pred[8];
					_mm_store_si128((__m128i*)pred, mpred);
					
					short *curr=rows[0];
					//short
					//	NW[]={rows[1][-4], rows[1][-3], rows[1][-2]},
					//	N []={rows[1][ 0], rows[1][ 1], rows[1][ 2]},
					//	W []={rows[0][-4], rows[0][-3], rows[0][-2]};

					//short pred[]=
					//{
					//	N[0]+W[0]-NW[0],
					//	N[1]+W[1]-NW[1],
					//	N[2]+W[2]-NW[2],
					//};
					//pred[0]=MEDIAN3(N[0], W[0], pred[0]);
					//pred[1]=MEDIAN3(N[1], W[1], pred[1]);
					//pred[2]=MEDIAN3(N[2], W[2], pred[2]);
					
					//short pred[3], vmin[3], vmax[3];		//X  slower
					//for(int k=0;k<3;++k)
					//{
					//	if(N[k]<W[k])
					//		vmin[k]=N[k], vmax[k]=W[k];
					//	else
					//		vmin[k]=W[k], vmax[k]=N[k];
					//	if((unsigned)(NW[k]-vmin[k])<(unsigned)(vmax[k]-vmin[k]))
					//		pred[k]=N[k]+W[k]-NW[k];
					//	else if(NW[k]<vmin[k])
					//		pred[k]=vmax[k];
					//	else
					//		pred[k]=vmin[k];
					//}

					rows[0]+=4;
					rows[1]+=4;
					PROF(PRED);
					PROF(DUMMY);

					short val[3];
#else
					short curr[3], val[3];
#endif
					
					val[1]=(short)ac_dec(ec+0, CDF+513*0, nlevels);
					val[2]=(short)ac_dec(ec+1, CDF+513*1, clevels);
					val[0]=(short)ac_dec(ec+2, CDF+513*2, clevels);
					PROF(EC);

#ifdef USE_CLAMPGRAD
					curr[0]=(short)(((val[0]+pred[0])&(clevels-1))-chalf);
					curr[1]=(short)(((val[1]+pred[1])&(nlevels-1))- half);
					curr[2]=(short)(((val[2]+pred[2])&(clevels-1))-chalf);
					
					short *rgb=dst->data+idx;
					memcpy(rgb, curr, sizeof(short[3]));
					rgb[1]-=(rgb[0]+rgb[2])>>2;
					rgb[2]+=rgb[1];
					rgb[0]+=rgb[1];
					PROF(RCT);
#else
					curr[1]=((val[1]+W[1])&(nlevels-1))-half;
					curr[2]=((val[2]+W[2])&(clevels-1))-chalf;
					curr[0]=((val[0]+W[0])&(clevels-1))-chalf;
					memcpy(W, curr, sizeof(short[3]));
					PROF(PRED);

					curr[1]-=(curr[0]+curr[2])>>2;
					curr[2]+=curr[1];
					curr[0]+=curr[1];
					memcpy(dst->data+idx, curr, sizeof(short[3]));
					PROF(RCT);
#endif

					update_hist(hist+513*0, nlevels, val[1]);
					update_hist(hist+513*1, clevels, val[2]);
					update_hist(hist+513*2, clevels, val[0]);
					PROF(HIST);
#ifdef ENABLE_GUIDE
#ifdef USE_CLAMPGRAD
					curr=rgb;
#endif
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
			ac_dec_init(ec+0, cbuf, cbuf+clen);
			for(int ky=0, idx=0;ky<image->ih;++ky)
			{
				int W=0;
				for(int kx=0;kx<image->iw;++kx, ++idx)
				{
					int val=ac_dec(ec, CDF+513*0, nlevels);
					dst->data[idx]=(short)(((val+W)&(nlevels-1))-half);
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
			ptrdiff_t csize=list[0].nobj+list[1].nobj+list[2].nobj;
			printf("YUV %12zd %12zd %12zd\n",
				list[0].nobj,
				list[1].nobj,
				list[2].nobj
			);
#endif
			printf("csize %12td  %10.6lf%%  CR %8.6lf\n", csize, 100.*csize/usize, (double)usize/csize);
		}
		printf("F05  %c %15.6lf sec\n", 'D'+fwd, t0);
		prof_print();
	}
	dlist_clear(list+0);
	dlist_clear(list+1);
	dlist_clear(list+2);
	free(hist);
	free(CDF);
#ifdef USE_CLAMPGRAD
	free(pixels);
#endif
	return 1;
}