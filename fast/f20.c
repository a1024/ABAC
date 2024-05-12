#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;


//	#define PROFILE_CSIZE
//	#define ENABLE_GUIDE
//	#define ALLOW_IDIV


#include"ac.h"
#ifdef ENABLE_GUIDE
static Image *guide=0;
#endif
static void calc_ctx(const __m128i *N, const __m128i *W, const __m128i *NW, const __m128i *eN, const __m128i *eW, const __m128i *mminmag, __m128i *pred, __m128i *mag)
{
	__m128i mmin=_mm_min_epi32(*N, *W);
	__m128i mmax=_mm_max_epi32(*N, *W);
	*pred=_mm_add_epi32(*N, *W);
	*mag=_mm_add_epi32(*eN, *eW);
	*pred=_mm_sub_epi32(*pred, *NW);
	*mag=_mm_srai_epi32(*mag, 1);//mag = (eN+eW)>>1		//mag = (4*(eN+eW)+eNE-eNW)>>2
	*pred=_mm_max_epi32(*pred, mmin);
	*mag=_mm_max_epi32(*mag, *mminmag);
	*pred=_mm_min_epi32(*pred, mmax);
}
static void pack_sign(__m128i *data)
{
	//x = curr-pred
	//x = x<<1^-(x<0)
	__m128i negmask=_mm_cmpgt_epi32(_mm_setzero_si128(), *data);
	*data=_mm_slli_epi32(*data, 1);
	*data=_mm_xor_si128(*data, negmask);
}
static void unpack_sign(__m128i *curr)
{
	//curr = (curr>>1^-(curr&1)) + pred
	__m128i negmask=_mm_and_si128(*curr, _mm_set1_epi32(1));
	negmask=_mm_cmpeq_epi32(negmask, _mm_set1_epi32(1));
	*curr=_mm_srai_epi32(*curr, 1);
	*curr=_mm_xor_si128(*curr, negmask);
}
static void calc_update(const __m128i *eN, const __m128i *eW, const __m128i *eNEEE, __m128i *ecurr)//writes to hi64 bits
{
	//x = (eN+eW+eNEEE+ecurr)>>2
	*ecurr=_mm_add_epi32(*ecurr, *eN);
	*ecurr=_mm_add_epi32(*ecurr, *eW);
	*ecurr=_mm_add_epi32(*ecurr, *eNEEE);
	*ecurr=_mm_srli_epi32(*ecurr, 2);
}
int f20_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;

	size_t ebufsize=sizeof(int[4*8])*(image->iw+8LL);//4 padded rows * 4 channels max * {pixels, errors}
	int *pixels=(int*)_mm_malloc(ebufsize, sizeof(__m128i));
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	int minmag=1<<image->depth>>7;
	pixels[0]=0;
	pixels[1]=0;
	pixels[2]=0;
	pixels[3]=0;
	pixels[4]=minmag;
	pixels[5]=minmag;
	pixels[6]=minmag;
	pixels[7]=minmag;
	memfill(pixels+8, pixels, ebufsize-sizeof(int[8]), sizeof(int[8]));
	__m128i mminmag=_mm_set1_epi32(minmag);
	if(image->depth==8)
	{
		mminmag=_mm_srli_epi32(mminmag, 2);
		minmag>>=2;
	}
	ptrdiff_t usize=(ptrdiff_t)image->iw*image->ih*image->nch*image->depth>>3;
	if(fwd)
	{
#ifdef ENABLE_GUIDE
		guide=image;
#endif
#ifdef PROFILE_CSIZE
		ptrdiff_t csizes[4*3]={0};//{unary, bypass}
#endif
		GolombRiceCoder ec;
		DList list;
		dlist_init(&list, 1, 1024, 0);
		gr_enc_init(&ec, &list);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			int *rows[]=
			{
				pixels+(((image->iw+8LL)*((ky-0LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-1LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-2LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-3LL)&3)+4)<<3),
			};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				__m128i
					NW	=_mm_load_si128((__m128i*)rows[1]-1*2+0),
					N	=_mm_load_si128((__m128i*)rows[1]+0*2+0),
					W	=_mm_load_si128((__m128i*)rows[0]-1*2+0),
					eN	=_mm_load_si128((__m128i*)rows[1]+0*2+1),
					eNEEE	=_mm_load_si128((__m128i*)rows[1]+3*2+1),
					eW	=_mm_load_si128((__m128i*)rows[0]-1*2+1);
				int	*curr	=rows[0]+0*8;
				__m128i mcurr, mpred, mmag;
				calc_ctx(&N, &W, &NW, &eN, &eW, &mminmag, &mpred, &mmag);
				//if(ky==1&&kx==5)//
				//if(idx==1170)//
				//if(idx==6936)//
				//	printf("");
				ALIGN(16) int mag[4];
				ALIGN(16) short c2[8];
				_mm_store_si128((__m128i*)mag, mmag);
				memcpy(c2, image->data+idx, sizeof(short)*image->nch);
				if(image->nch>=3)
				{
					//JPEG2000 RCT
					c2[0]-=c2[1];
					c2[2]-=c2[1];
					c2[1]+=(c2[0]+c2[2])>>2;
				}
				mcurr=_mm_load_si128((__m128i*)c2);
				mcurr=_mm_cvtepi16_epi32(mcurr);
				__m128i delta=_mm_sub_epi32(mcurr, mpred);
				pack_sign(&delta);
				_mm_store_si128((__m128i*)curr+0, mcurr);
				_mm_store_si128((__m128i*)curr+1, delta);
				switch(image->nch)
				{
#ifdef ALLOW_IDIV
				case 4:gr_enc(&ec, curr[4+3], mag[3]+1);
				case 3:gr_enc(&ec, curr[4+2], mag[2]+1);
				case 2:gr_enc(&ec, curr[4+1], mag[1]+1);
				case 1:gr_enc(&ec, curr[4+0], mag[0]+1);
#else
				case 4:gr_enc_POT(&ec, curr[4+3], FLOOR_LOG2(mag[3]+1));
				case 3:gr_enc_POT(&ec, curr[4+2], FLOOR_LOG2(mag[2]+1));
				case 2:gr_enc_POT(&ec, curr[4+1], FLOOR_LOG2(mag[1]+1));
				case 1:gr_enc_POT(&ec, curr[4+0], FLOOR_LOG2(mag[0]+1));
#endif
					break;
				}
#ifdef PROFILE_CSIZE
				for(int kc=image->nch-1;kc>=0;--kc)
				{
					int sym=curr[kc+4];
					int magnitude=mag[kc]+1;
					csizes[0<<2|kc]+=sym/magnitude;
					++csizes[1<<2|kc];
					int bypass=sym%magnitude;
					int nbypass=floor_log2_32(magnitude)+1;
					csizes[2<<2|kc]+=nbypass-(bypass<(1<<nbypass)-magnitude);
				}
#endif
				calc_update(&eN, &eW, &eNEEE, &delta);
				_mm_store_si128((__m128i*)curr+1, delta);
				//if(curr[4+1]>=(4<<image->depth))
				//	LOG_ERROR("");
				rows[0]+=8;
				rows[1]+=8;
				rows[2]+=8;
				rows[3]+=8;
			}
		}
		gr_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
		if(loud)
		{
			ptrdiff_t csize=list.nobj;
			t0=time_sec()-t0;
			printf("%14td/%14td = %10.6lf%%  CR %lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
#ifdef PROFILE_CSIZE
			printf("C         unary         stop_bit          bypass\n");
			double total_unary=0, total_stop=0, total_bypass=0;
			for(int kc=0;kc<image->nch;++kc)
			{
				double size_unary=csizes[0<<2|kc]/8., size_stop=csizes[1<<2|kc]/8., size_bypass=csizes[2<<2|kc]/8.;
				printf("%d %16.2lf %16.2lf %16.2lf\n", kc, size_unary, size_stop, size_bypass);
				total_unary+=size_unary;
				total_stop+=size_stop;
				total_bypass+=size_bypass;
			}
			printf("T %16.2lf %16.2lf %16.2lf\n", total_unary, total_stop, total_bypass);
#endif
			printf("E  %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
		dlist_clear(&list);
	}
	else
	{
		GolombRiceCoder ec;
		gr_dec_init(&ec, cbuf, cbuf+clen);
		__m128i packpixels=_mm_set_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1, 13, 12,  9,  8,  5,  4,  1,  0
			//15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
		);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			int *rows[]=
			{
				pixels+(((image->iw+8LL)*((ky-0LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-1LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-2LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-3LL)&3)+4)<<3),
			};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				__m128i
					NW	=_mm_load_si128((__m128i*)rows[1]-1*2+0),
					N	=_mm_load_si128((__m128i*)rows[1]+0*2+0),
					W	=_mm_load_si128((__m128i*)rows[0]-1*2+0),
					eN	=_mm_load_si128((__m128i*)rows[1]+0*2+1),
					eNEEE	=_mm_load_si128((__m128i*)rows[1]+3*2+1),
					eW	=_mm_load_si128((__m128i*)rows[0]-1*2+1);
				int	*curr	=rows[0]+0*8;
				__m128i mpred, mmag;
				calc_ctx(&N, &W, &NW, &eN, &eW, &mminmag, &mpred, &mmag);
				//if(ky==1&&kx==5)//
				//if(idx==1161)//
				//if(idx==6936)//
				//	printf("");
				ALIGN(16) int mag[4];
				_mm_store_si128((__m128i*)mag, mmag);
				switch(image->nch)
				{
#ifdef ALLOW_IDIV
				case 4:curr[3+4]=gr_dec(&ec, mag[3]+1);
				case 3:curr[2+4]=gr_dec(&ec, mag[2]+1);
				case 2:curr[1+4]=gr_dec(&ec, mag[1]+1);
				case 1:curr[0+4]=gr_dec(&ec, mag[0]+1);
#else
				case 4:curr[3+4]=gr_dec_POT(&ec, FLOOR_LOG2(mag[3]+1));
				case 3:curr[2+4]=gr_dec_POT(&ec, FLOOR_LOG2(mag[2]+1));
				case 2:curr[1+4]=gr_dec_POT(&ec, FLOOR_LOG2(mag[1]+1));
				case 1:curr[0+4]=gr_dec_POT(&ec, FLOOR_LOG2(mag[0]+1));
#endif
					break;
				}
				__m128i delta=_mm_load_si128((__m128i*)curr+1);
				__m128i mcurr=delta;
				unpack_sign(&mcurr);
				calc_update(&eN, &eW, &eNEEE, &delta);
				mcurr=_mm_add_epi32(mcurr, mpred);
				_mm_store_si128((__m128i*)curr+0, mcurr);
				_mm_store_si128((__m128i*)curr+1, delta);

				mcurr=_mm_shuffle_epi8(mcurr, packpixels);
				short *c2=image->data+idx;
				ALIGN(16) short c3[8];
				_mm_store_si128((__m128i*)c3, mcurr);
				memcpy(c2, c3, sizeof(short)*image->nch);
				if(image->nch>=3)
				{
					//JPEG2000 RCT
					c2[1]-=(c2[0]+c2[2])>>2;
					c2[2]+=c2[1];
					c2[0]+=c2[1];
				}
#ifdef ENABLE_GUIDE
				if(guide&&memcmp(c2, guide->data+idx, sizeof(short)*image->nch))
				{
					short c0[4]={0};
					memcpy(c0, guide->data+idx, sizeof(short)*image->nch);
					c0[0]-=c0[1];
					c0[2]-=c0[1];
					c0[1]+=(c0[0]+c0[2])>>2;
					curr[0]-=curr[1];
					curr[2]-=curr[1];
					curr[1]+=(curr[0]+curr[2])>>2;
					LOG_ERROR("Guide error IDX %d/%d", idx, image->nch*image->iw*image->ih);
					printf("");//
				}
#endif
				rows[0]+=8;
				rows[1]+=8;
				rows[2]+=8;
				rows[3]+=8;
			}
		}
		if(loud)
		{
			t0=time_sec()-t0;
			printf("D  %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
	}
	_mm_free(pixels);
	return 0;
}