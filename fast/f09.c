#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
	#define SIMPLE_PREDS
//	#define PROFILER 1


#ifdef PROFILER
#define CHECKPOINTLIST\
	CHECKPOINT(INIT)\
	CHECKPOINT(WP)\
	CHECKPOINT(CTX)\
	CHECKPOINT(EC)\
	CHECKPOINT(DUMMY)\
	CHECKPOINT(HIST)\
	CHECKPOINT(CDF)\
	CHECKPOINT(FINISH)
#endif
#include"ac.h"
#include"profiler.h"
#define CDF_UPDATE_MASK 1023
#define QLEVELS 7
#ifdef SIMPLE_PREDS
#define DEPTH_C0 8
#define DEPTH_C1 9
#define DEPTH_C2 9
#define NLEVELS_C0 256
#define NLEVELS_C1 512
#define NLEVELS_C2 512
#else
#define DEPTH_C0 9
#define DEPTH_C1 8
#define DEPTH_C2 9
#define NLEVELS_C0 512
#define NLEVELS_C1 256
#define NLEVELS_C2 512
#endif
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
static void update_CDF(unsigned short *CDF, const unsigned short *hist, int nlevels)
{
	int sum=hist[nlevels], c=0;
	if(CDF[nlevels]!=sum)
	{
		for(int ks=0;ks<nlevels;++ks)
		{
			int freq=hist[ks];
			//if(c>sum)
			//	LOG_ERROR("");
			CDF[ks]=(int)(c*(0x10000LL-nlevels)/sum)+ks;
			//int LOL_1=(int)(c*(0x10000LL-nlevels)/sum)+ks;
			//if(LOL_1>0x10000)
			//	LOG_ERROR("");
			//CDF[ks]=LOL_1;
			c+=freq;
			//if(ks&&CDF[ks]<CDF[ks-1])
			//	LOG_ERROR("");
		}
		CDF[nlevels]=sum;//0x10000 doesn't fit in uint16 anyway, so we store last histsum here
	}
}
int f09_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	PROF_START();
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	if(image->depth!=8)
		LOG_ERROR("Unsupported bit depth %d", image->depth);
	if(image->nch!=3)
		LOG_ERROR("Unsupported number of channels %d", image->nch);
#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif
	char depths[]=
	{
#ifdef SIMPLE_PREDS
		image->depth,
		image->depth+1,
		image->depth+1,
#else
		image->depth+1,
		image->depth,
		image->depth+1,
#endif
	};
	int nlevels[]=
	{
		1<<depths[0],
		1<<depths[1],
		1<<depths[2],
	};
	Image im2={0}, *im=0;
	if(fwd)
	{
		image_copy(&im2, src);
		im=&im2;
		if(!im->data)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		PROF(INIT);
		
#ifndef ENABLE_GUIDE
#ifdef SIMPLE_PREDS
		rct_JPEG2000_32(im, 1);
		pred_clampgrad(im, 1, depths);
#else
		pred_wp_deferred(im, 1);
#endif
		PROF(WP);
#endif
	}
	else
		im=dst;
	DList list[3]={0};
	dlist_init(list+0, 1, 1024, 0);
	dlist_init(list+1, 1, 1024, 0);
	dlist_init(list+2, 1, 1024, 0);
	ArithmeticCoder ec[3]={0};
	//ANSCoder ec[3]={0};
	unsigned short *stats[]=
	{
		(unsigned short*)_mm_malloc(sizeof(short[(NLEVELS_C0+1LL)*2*QLEVELS*QLEVELS]), sizeof(__m256i)),//nlevels * {hist, CDF} * QLEVELS bins ^ {eN, eW}
		(unsigned short*)_mm_malloc(sizeof(short[(NLEVELS_C1+1LL)*2*QLEVELS*QLEVELS]), sizeof(__m256i)),
		(unsigned short*)_mm_malloc(sizeof(short[(NLEVELS_C2+1LL)*2*QLEVELS*QLEVELS]), sizeof(__m256i)),
	};
	if(!stats[0]||!stats[1]||!stats[2])
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int kc=0;kc<3;++kc)//initialize stats
	{
		int n=nlevels[kc], depth=depths[kc];
		unsigned short *curr_stats=stats[kc];
		curr_stats[0]=1;
		memfill(curr_stats+1, curr_stats, (n+1LL)*sizeof(short)-sizeof(*curr_stats)*2, sizeof(*curr_stats));
		curr_stats[n]=n;//histsum
		for(int ks=0;ks<=n;++ks)
			curr_stats[n+1+ks]=ks<<(16-depth);
		memfill(curr_stats+2*(n+1LL), curr_stats, (n+1LL)*(sizeof(short[2*QLEVELS*QLEVELS])-sizeof(short[2])), (n+1LL)*sizeof(short[2]));
	}
	ALIGN(16) int ctx[4]={0};
	ALIGN(16) int sym[4]={0};
	int w3=im->iw*3;
	__m128i qlevels[]=
	{
		//_mm_set1_epi32(-128),
		_mm_set1_epi32( -32),
		_mm_set1_epi32(  -8),
		_mm_set1_epi32(  -2),
		_mm_set1_epi32(   2),
		_mm_set1_epi32(   8),
		_mm_set1_epi32(  32),
		//_mm_set1_epi32( 128),
	};
	PROF(INIT);
	if(fwd)
	{
#if 1
		ac_enc_init(ec+0, list+0);
		ac_enc_init(ec+1, list+1);
		ac_enc_init(ec+2, list+2);
#else
		ans_enc_init(ec+0, list+0);
		ans_enc_init(ec+1, list+1);
		ans_enc_init(ec+2, list+2);
#endif
		//for(int ky=im->ih-1, idx=im->iw*im->ih*3-3;ky>=0;--ky)
		//	for(int kx=im->iw-1;kx>=0;--kx, idx-=3)
		for(int ky=0, idx=0;ky<im->ih;++ky)
		{
			for(int kx=0;kx<im->iw;++kx, idx+=3)
			{
				//if(idx==2319)//
				//if(idx==148554)//
				//	printf("");

				__m128i eN=_mm_set_epi32(
					0,
					ky?im->data[idx-w3+2]:0,
					ky?im->data[idx-w3+1]:0,
					ky?im->data[idx-w3+0]:0
				);
				__m128i eW=_mm_set_epi32(
					0,
					kx?im->data[idx-3+2]:0,
					kx?im->data[idx-3+1]:0,
					kx?im->data[idx-3+0]:0
				);
				__m128i qN=_mm_cmpgt_epi32(eN, qlevels[0]);
				__m128i qW=_mm_cmpgt_epi32(eW, qlevels[0]);
				for(int k=1;k<_countof(qlevels);++k)
				{
					qN=_mm_add_epi32(qN, _mm_cmpgt_epi32(eN, qlevels[k]));
					qW=_mm_add_epi32(qW, _mm_cmpgt_epi32(eW, qlevels[k]));
				}
				qN=_mm_sub_epi32(_mm_setzero_si128(), qN);
				qW=_mm_sub_epi32(_mm_setzero_si128(), qW);
				__m128i mctx=_mm_add_epi32(_mm_mullo_epi32(qN, _mm_set1_epi32(QLEVELS)), qW);
				_mm_store_si128((__m128i*)ctx, mctx);
				PROF(CTX);
				unsigned short *curr_stats[]=
				{
					stats[0]+(NLEVELS_C0+1LL)*2*ctx[0],
					stats[1]+(NLEVELS_C1+1LL)*2*ctx[1],
					stats[2]+(NLEVELS_C2+1LL)*2*ctx[2],
				};
				const unsigned short *CDFs[]=
				{
					curr_stats[0]+NLEVELS_C0+1LL,
					curr_stats[1]+NLEVELS_C1+1LL,
					curr_stats[2]+NLEVELS_C2+1LL,
				};

				__m128i mcurr=_mm_set_epi32(
					0,
					im->data[idx+2],
					im->data[idx+1],
					im->data[idx+0]
				);
				
				//if(!idx)//
				//if(idx==195)//
				//	printf("");

				mcurr=_mm_xor_si128(_mm_slli_epi32(mcurr, 1), _mm_srai_epi32(mcurr, 31));//pack sign		x = x<<1^-(x<0)
				_mm_store_si128((__m128i*)sym, mcurr);
#if 1
				ac_enc_packedCDF(ec+0, sym[0], CDFs[0], NLEVELS_C0);
				ac_enc_packedCDF(ec+1, sym[1], CDFs[1], NLEVELS_C1);
				ac_enc_packedCDF(ec+2, sym[2], CDFs[2], NLEVELS_C2);
#else
				unsigned cdf[]=
				{
					CDFs[0][sym[0]],
					CDFs[1][sym[1]],
					CDFs[2][sym[2]],
				};
				int freq[]=
				{
					//CDFs[0][sym[0]+1]-cdf[0],//X  0x10000 doesn't fit in unsigned short!
					//CDFs[1][sym[1]+1]-cdf[1],
					//CDFs[2][sym[2]+1]-cdf[2],
					(sym[0]==511?0x10000:CDFs[0][sym[0]+1])-cdf[0],
					(sym[1]==255?0x10000:CDFs[1][sym[1]+1])-cdf[1],
					(sym[2]==511?0x10000:CDFs[2][sym[2]+1])-cdf[2],
				};
				if((ec[2].state>>16)>=(unsigned)freq[2])
				{
					dlist_push_back(ec[2].list, &ec[2].state, 2);
					ec[2].state>>=16;
				}
				if((ec[1].state>>16)>=(unsigned)freq[1])
				{
					dlist_push_back(ec[1].list, &ec[1].state, 2);
					ec[1].state>>=16;
				}
				if((ec[0].state>>16)>=(unsigned)freq[0])//renorm
				{
					dlist_push_back(ec[0].list, &ec[0].state, 2);
					ec[0].state>>=16;
				}
				debug_enc_update(ec[2].state, cdf[2], freq[2], 0, 0, 0, 0, sym[2]);
				debug_enc_update(ec[1].state, cdf[1], freq[1], 0, 0, 0, 0, sym[1]);
				debug_enc_update(ec[0].state, cdf[0], freq[0], 0, 0, 0, 0, sym[0]);
				div_t qr[]=
				{
					div(ec[0].state, (unsigned short)freq[0]),
					div(ec[1].state, (unsigned short)freq[1]),
					div(ec[2].state, (unsigned short)freq[2]),
				};
				ec[2].state=qr[2].quot<<16|(cdf[2]+qr[2].rem);//update
				ec[1].state=qr[1].quot<<16|(cdf[1]+qr[1].rem);
				ec[0].state=qr[0].quot<<16|(cdf[0]+qr[0].rem);
#endif
				PROF(EC);
				PROF(DUMMY);

				for(int kc=0;kc<3;++kc)//update hist
				{
					int n=nlevels[kc];
					unsigned short *hist=curr_stats[kc];
					++hist[sym[kc]];
					++hist[n];
					if(hist[n]>=1024)
					{
						int sum=0;
						for(int ks=0;ks<n;++ks)
							sum+=hist[ks]=(hist[ks]+1)>>1;
						hist[n]=sum;
					}
				}
				PROF(HIST);
				if(idx&&!(idx&CDF_UPDATE_MASK))
				{
					unsigned short *stats2[]=
					{
						stats[0],
						stats[1],
						stats[2],
					};
					for(int k=0;k<QLEVELS*QLEVELS;++k)
					{
						update_CDF(stats2[0]+NLEVELS_C0+1LL, stats2[0], NLEVELS_C0);
						update_CDF(stats2[1]+NLEVELS_C1+1LL, stats2[1], NLEVELS_C1);
						update_CDF(stats2[2]+NLEVELS_C2+1LL, stats2[2], NLEVELS_C2);
						stats2[0]+=(NLEVELS_C0+1LL)*2;
						stats2[1]+=(NLEVELS_C1+1LL)*2;
						stats2[2]+=(NLEVELS_C2+1LL)*2;
					}
					PROF(CDF);
				}
			}
		}
#if 1
		ac_enc_flush(ec+0);
		ac_enc_flush(ec+1);
		ac_enc_flush(ec+2);
#else
		ans_enc_flush(ec+0);
		ans_enc_flush(ec+1);
		ans_enc_flush(ec+2);
#endif
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
		PROF(FINISH);
	}
	else
	{
		unsigned bm[3];
		memcpy(bm, cbuf, sizeof(bm));
		cbuf+=sizeof(bm);
		clen-=sizeof(bm);
		const unsigned char *ptr=cbuf;
#if 1
		ac_dec_init(ec+0, ptr, ptr+bm[0]);	ptr+=bm[0];
		ac_dec_init(ec+1, ptr, ptr+bm[1]);	ptr+=bm[1];
		ac_dec_init(ec+2, ptr, ptr+bm[2]);	ptr+=bm[2];
#else
		ans_dec_init(ec+0, ptr, ptr+bm[0]);	ptr+=bm[0];
		ans_dec_init(ec+1, ptr, ptr+bm[1]);	ptr+=bm[1];
		ans_dec_init(ec+2, ptr, ptr+bm[2]);	ptr+=bm[2];
#endif
		for(int ky=0, idx=0;ky<im->ih;++ky)
		{
			for(int kx=0;kx<im->iw;++kx, idx+=3)
			{
				//if(idx==2319)//
				//if(idx==148554)//
				//	printf("");

				__m128i eN=_mm_set_epi32(
					0,
					ky?im->data[idx-w3+2]:0,
					ky?im->data[idx-w3+1]:0,
					ky?im->data[idx-w3+0]:0
				);
				__m128i eW=_mm_set_epi32(
					0,
					kx?im->data[idx-3+2]:0,
					kx?im->data[idx-3+1]:0,
					kx?im->data[idx-3+0]:0
				);
				__m128i qN=_mm_cmpgt_epi32(eN, qlevels[0]);
				__m128i qW=_mm_cmpgt_epi32(eW, qlevels[0]);
				for(int k=1;k<_countof(qlevels);++k)
				{
					qN=_mm_add_epi32(qN, _mm_cmpgt_epi32(eN, qlevels[k]));
					qW=_mm_add_epi32(qW, _mm_cmpgt_epi32(eW, qlevels[k]));
				}
				qN=_mm_sub_epi32(_mm_setzero_si128(), qN);
				qW=_mm_sub_epi32(_mm_setzero_si128(), qW);
				__m128i mctx=_mm_add_epi32(_mm_mullo_epi32(qN, _mm_set1_epi32(QLEVELS)), qW);
				_mm_store_si128((__m128i*)ctx, mctx);
				unsigned short *curr_stats[]=
				{
					stats[0]+(NLEVELS_C0+1LL)*2*ctx[0],
					stats[1]+(NLEVELS_C1+1LL)*2*ctx[1],
					stats[2]+(NLEVELS_C2+1LL)*2*ctx[2],
				};
				const unsigned short *CDFs[]=
				{
					curr_stats[0]+NLEVELS_C0+1LL,
					curr_stats[1]+NLEVELS_C1+1LL,
					curr_stats[2]+NLEVELS_C2+1LL,
				};
				
#if 1
				sym[0]=ac_dec_packedCDF_POT(ec+0, CDFs[0], DEPTH_C0);
				sym[1]=ac_dec_packedCDF_POT(ec+1, CDFs[1], DEPTH_C1);
				sym[2]=ac_dec_packedCDF_POT(ec+2, CDFs[2], DEPTH_C2);
#else
				unsigned short c[]=
				{
					(unsigned short)ec[0].state,
					(unsigned short)ec[1].state,
					(unsigned short)ec[2].state,
				};
				sym[0] =(CDFs[0][       256]<=c[0])<<8;
				sym[2] =(CDFs[2][       256]<=c[2])<<8;

				sym[0]|=(CDFs[0][sym[0]|128]<=c[0])<<7;
				sym[1] =(CDFs[1][       128]<=c[1])<<7;
				sym[2]|=(CDFs[2][sym[2]|128]<=c[2])<<7;

				sym[0]|=(CDFs[0][sym[0]| 64]<=c[0])<<6;
				sym[1]|=(CDFs[1][sym[1]| 64]<=c[1])<<6;
				sym[2]|=(CDFs[2][sym[2]| 64]<=c[2])<<6;
				sym[0]|=(CDFs[0][sym[0]| 32]<=c[0])<<5;
				sym[1]|=(CDFs[1][sym[1]| 32]<=c[1])<<5;
				sym[2]|=(CDFs[2][sym[2]| 32]<=c[2])<<5;
				sym[0]|=(CDFs[0][sym[0]| 16]<=c[0])<<4;
				sym[1]|=(CDFs[1][sym[1]| 16]<=c[1])<<4;
				sym[2]|=(CDFs[2][sym[2]| 16]<=c[2])<<4;
				sym[0]|=(CDFs[0][sym[0]|  8]<=c[0])<<3;
				sym[1]|=(CDFs[1][sym[1]|  8]<=c[1])<<3;
				sym[2]|=(CDFs[2][sym[2]|  8]<=c[2])<<3;
				sym[0]|=(CDFs[0][sym[0]|  4]<=c[0])<<2;
				sym[1]|=(CDFs[1][sym[1]|  4]<=c[1])<<2;
				sym[2]|=(CDFs[2][sym[2]|  4]<=c[2])<<2;
				sym[0]|=(CDFs[0][sym[0]|  2]<=c[0])<<1;
				sym[1]|=(CDFs[1][sym[1]|  2]<=c[1])<<1;
				sym[2]|=(CDFs[2][sym[2]|  2]<=c[2])<<1;

				sym[0]|= CDFs[0][sym[0]|  1]<=c[0];
				sym[1]|= CDFs[1][sym[1]|  1]<=c[1];
				sym[2]|= CDFs[2][sym[2]|  1]<=c[2];

				unsigned cdf[]=
				{
					CDFs[0][sym[0]],
					CDFs[1][sym[1]],
					CDFs[2][sym[2]],
				};
				int freq[]=
				{
					(sym[0]==511?0x10000:CDFs[0][sym[0]+1])-cdf[0],
					(sym[1]==255?0x10000:CDFs[1][sym[1]+1])-cdf[1],
					(sym[2]==511?0x10000:CDFs[2][sym[2]+1])-cdf[2],
				};

				//if(freq[0]==0xFFFF||freq[1]==0xFFFF||freq[2]==0xFFFF)//
				//if(idx==195)//
				//	printf("");

				debug_dec_update(ec[0].state, cdf[0], freq[0], 0, 0, 0, 0, sym[0]);
				debug_dec_update(ec[1].state, cdf[1], freq[1], 0, 0, 0, 0, sym[1]);
				debug_dec_update(ec[2].state, cdf[2], freq[2], 0, 0, 0, 0, sym[2]);
				ec[0].state=(unsigned short)freq[0]*(unsigned short)(ec[0].state>>16)+c[0]-cdf[0];//update
				ec[1].state=(unsigned short)freq[1]*(unsigned short)(ec[1].state>>16)+c[1]-cdf[1];
				ec[2].state=(unsigned short)freq[2]*(unsigned short)(ec[2].state>>16)+c[2]-cdf[2];
				if(ec[0].state<0x10000)//renorm
				{
					ec[0].state<<=16;
					if(ec[0].srcptr-2>=ec[0].srcstart)
					{
						ec[0].srcptr-=2;
						memcpy(&ec[0].state, ec[0].srcptr, 2);
					}
				}
				if(ec[1].state<0x10000)
				{
					ec[1].state<<=16;
					if(ec[1].srcptr-2>=ec[1].srcstart)
					{
						ec[1].srcptr-=2;
						memcpy(&ec[1].state, ec[1].srcptr, 2);
					}
				}
				if(ec[2].state<0x10000)
				{
					ec[2].state<<=16;
					if(ec[2].srcptr-2>=ec[2].srcstart)
					{
						ec[2].srcptr-=2;
						memcpy(&ec[2].state, ec[2].srcptr, 2);
					}
				}
#endif
				__m128i mcurr=_mm_load_si128((__m128i*)sym);
				//unpack sign		x = x>>1^-(x&1)
				__m128i one=_mm_set1_epi32(1);
				mcurr=_mm_xor_si128(_mm_srli_epi32(mcurr, 1), _mm_cmpeq_epi32(_mm_and_si128(mcurr, one), one));
				_mm_store_si128((__m128i*)ctx, mcurr);
				dst->data[idx+0]=ctx[0];
				dst->data[idx+1]=ctx[1];
				dst->data[idx+2]=ctx[2];
				for(int kc=0;kc<3;++kc)//update hist
				{
					int n=nlevels[kc];
					unsigned short *hist=curr_stats[kc];
					++hist[sym[kc]];
					++hist[n];
					if(hist[n]>=1024)
					{
						int sum=0;
						for(int ks=0;ks<n;++ks)
							sum+=hist[ks]=(hist[ks]+1)>>1;
						hist[n]=sum;
					}
				}
				PROF(HIST);
				if(idx&&!(idx&CDF_UPDATE_MASK))
				{
					unsigned short *stats2[]=
					{
						stats[0],
						stats[1],
						stats[2],
					};
					for(int k=0;k<QLEVELS*QLEVELS;++k)
					{
						update_CDF(stats2[0]+NLEVELS_C0+1LL, stats2[0], NLEVELS_C0);
						update_CDF(stats2[1]+NLEVELS_C1+1LL, stats2[1], NLEVELS_C1);
						update_CDF(stats2[2]+NLEVELS_C2+1LL, stats2[2], NLEVELS_C2);
						stats2[0]+=(NLEVELS_C0+1LL)*2;
						stats2[1]+=(NLEVELS_C1+1LL)*2;
						stats2[2]+=(NLEVELS_C2+1LL)*2;
					}
					PROF(CDF);
				}
#ifdef ENABLE_GUIDE
				if(guide&&memcmp(dst->data+idx, guide->data+idx, im->nch*sizeof(short)))
				{
					//short orig[4]={0};
					//memcpy(orig, guide->data+idx, im->nch*sizeof(short));
					//curr[0]-=curr[1];
					//curr[2]-=curr[1];
					//curr[1]+=(curr[0]+curr[2])>>2;
					//orig[0]-=orig[1];
					//orig[2]-=orig[1];
					//orig[1]+=(orig[0]+orig[2])>>2;
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
			}
		}
#ifndef ENABLE_GUIDE
#ifdef SIMPLE_PREDS
		pred_clampgrad(im, 0, depths);
		rct_JPEG2000_32(im, 0);
#else
		pred_wp_deferred(dst, 0);
#endif
		PROF(WP);
#endif
	}
	if(loud)
	{
		t0=time_sec()-t0;
		if(fwd)
		{
			ptrdiff_t usize=image_getBMPsize(image);
			ptrdiff_t csize=list[0].nobj+list[1].nobj+list[2].nobj;
			printf("YUV %12lld %12lld %12lld\n",
				list[1].nobj,
				list[2].nobj,
				list[0].nobj
			);
			printf("csize %12lld  %10.6lf%%  CR %8.6lf\n", csize, 100.*csize/usize, (double)usize/csize);
		}
		printf("F08  %c %15.6lf sec\n", 'D'+fwd, t0);
		prof_print();
	}
	if(fwd)
		image_clear(&im2);
	dlist_clear(list+0);
	dlist_clear(list+1);
	dlist_clear(list+2);
	_mm_free(stats[0]);
	_mm_free(stats[1]);
	_mm_free(stats[2]);
	return 1;
}