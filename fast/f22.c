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


//	#define ENABLE_GUIDE
	#define ALLOW_AVX512


#include"ac.h"
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
#ifdef ALLOW_AVX512
static void calc_ctx_512(const __m512i *N, const __m512i *W, const __m512i *NW, const __m512i *eW, const __m512i *mminmag, __m512i *pred, __m512i *mag)
{
	//pred = median3(N, W, N+W-NW)
	//mag = eW
	//mag = FLOOR_LOG2(mag+1)
	__m512i mmin=_mm512_min_epi32(*N, *W);
	__m512i mmax=_mm512_max_epi32(*N, *W);
	*pred=_mm512_add_epi32(*N, *W);
	*pred=_mm512_sub_epi32(*pred, *NW);
	*pred=_mm512_max_epi32(*pred, mmin);
	*pred=_mm512_min_epi32(*pred, mmax);

	*mag=_mm512_max_epi32(*eW, *mminmag);
	*mag=_mm512_add_epi32(*mag, _mm512_set1_epi32(1));
	*mag=_mm512_lzcnt_epi32(*mag);
	*mag=_mm512_sub_epi32(_mm512_set1_epi32(31), *mag);
}
static void pack_sign_512(__m512i *data)
{
	//x = curr-pred
	//x = x<<1 ^ x>>31
	__m512i negmask=_mm512_srai_epi32(*data, 31);
	*data=_mm512_slli_epi32(*data, 1);
	*data=_mm512_xor_si512(*data, negmask);
}
static void unpack_sign_512(__m512i *curr)
{
	//curr = (curr>>1^-(curr&1)) + pred
	__mmask16 negmask=_mm512_test_epi32_mask(*curr, _mm512_set1_epi32(1));
	*curr=_mm512_srai_epi32(*curr, 1);
	*curr=_mm512_mask_xor_epi32(*curr, negmask, *curr, _mm512_set1_epi32(-1));
}
static void calc_update_512(const __m512i *eW, const __m512i *eNEEE, __m512i *ecurr)//writes to hi64 bits
{
	//x = (2*eW+eNEEE+ecurr)>>2
	*ecurr=_mm512_add_epi32(*ecurr, _mm512_slli_epi32(*eW, 1));
	*ecurr=_mm512_add_epi32(*ecurr, *eNEEE);
	*ecurr=_mm512_srai_epi32(*ecurr, 2);//doesn't matter logic or arithmetic shift
}
#endif
static void calc_ctx_256x2(const __m256i *N, const __m256i *W, const __m256i *NW, const __m256i *eW, const __m256i *mminmag, __m256i *pred, __m256i *mag)
{
	//pred = median(N, W, N+W-NW)
	//mag = eW
	//mag = FLOOR_LOG2(PACK_SIGN(mag)+1)
	__m256i mmin0=_mm256_min_epi32(N[0], W[0]);
	__m256i mmin1=_mm256_min_epi32(N[1], W[1]);
	__m256i mmax0=_mm256_max_epi32(N[0], W[0]);
	__m256i mmax1=_mm256_max_epi32(N[1], W[1]);
	pred[0]=_mm256_add_epi32(N[0], W[0]);
	pred[1]=_mm256_add_epi32(N[1], W[1]);
	pred[0]=_mm256_sub_epi32(pred[0], NW[0]);
	pred[1]=_mm256_sub_epi32(pred[1], NW[1]);
	pred[0]=_mm256_max_epi32(pred[0], mmin0);
	pred[1]=_mm256_max_epi32(pred[1], mmin1);
	pred[0]=_mm256_min_epi32(pred[0], mmax0);
	pred[1]=_mm256_min_epi32(pred[1], mmax1);

	mag[0]=_mm256_max_epi32(eW[0], *mminmag);
	mag[1]=_mm256_max_epi32(eW[1], *mminmag);
	mag[0]=_mm256_add_epi32(mag[0], _mm256_set1_epi32(1));
	mag[1]=_mm256_add_epi32(mag[1], _mm256_set1_epi32(1));
}
static void pack_sign_256x2(__m256i *data)
{
	//x = curr-pred
	//x = x<<1 ^ x>>31
	__m256i negmask[]=
	{
		_mm256_srai_epi32(data[0], 31),
		_mm256_srai_epi32(data[1], 31),
	};
	data[0]=_mm256_slli_epi32(data[0], 1);
	data[1]=_mm256_slli_epi32(data[1], 1);
	data[0]=_mm256_xor_si256(data[0], negmask[0]);
	data[1]=_mm256_xor_si256(data[1], negmask[1]);
}
static void unpack_sign_256x2(__m256i *curr)
{
	//curr = (curr>>1^-(curr&1)) + pred
	__m256i negmask[]=
	{
		_mm256_and_si256(curr[0], _mm256_set1_epi32(1)),
		_mm256_and_si256(curr[1], _mm256_set1_epi32(1)),
	};
	negmask[0]=_mm256_cmpeq_epi32(negmask[0], _mm256_set1_epi32(1));
	negmask[1]=_mm256_cmpeq_epi32(negmask[1], _mm256_set1_epi32(1));
	curr[0]=_mm256_srai_epi32(curr[0], 1);
	curr[1]=_mm256_srai_epi32(curr[1], 1);
	curr[0]=_mm256_xor_si256(curr[0], negmask[0]);
	curr[1]=_mm256_xor_si256(curr[1], negmask[1]);
}
static void calc_update_256x2(const __m256i *eW, const __m256i *eNEEE, __m256i *ecurr)//writes to hi64 bits
{
	//x = (2*eW+eNEEE+ecurr)>>2
	ecurr[0]=_mm256_add_epi32(ecurr[0], _mm256_slli_epi32(eW[0], 1));
	ecurr[1]=_mm256_add_epi32(ecurr[1], _mm256_slli_epi32(eW[1], 1));
	ecurr[0]=_mm256_add_epi32(ecurr[0], eNEEE[0]);
	ecurr[1]=_mm256_add_epi32(ecurr[1], eNEEE[1]);
	ecurr[0]=_mm256_srai_epi32(ecurr[0], 2);
	ecurr[1]=_mm256_srai_epi32(ecurr[1], 2);
}
int f22_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif
	size_t ebufsize=sizeof(int[2*4*2*4])*(image->iw+16LL);//2 padded rows * 4 channels max * {pixels, errors} * 4 threads (interleaved)
	int *pixels=(int*)_mm_malloc(ebufsize, sizeof(__m512i));
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	int minmag=1<<image->depth>>7;
	memset(pixels, 0, sizeof(int[16]));
	pixels[16]=minmag;
	memfill(pixels+17, pixels+16, sizeof(int[15]), sizeof(int));
	memfill(pixels+32, pixels, ebufsize-sizeof(int[32]), sizeof(int[32]));
	ptrdiff_t nvals=(ptrdiff_t)image->iw*image->ih*image->nch;
	ptrdiff_t usize=nvals*image->depth>>3;
	int ceilqheight=(image->ih+3)>>2;
	int blocksize=ceilqheight*image->iw*image->nch;
	const short *endptr=image->data+nvals;
	GolombRiceCoder ec[4];
	__m256i stride=_mm256_set1_epi64x(sizeof(short)*image->nch);
#ifdef ALLOW_AVX512
	int cpu_features=get_cpu_features();
	if(cpu_features>>1&1)
	{
		__m512i mminmag=_mm512_set1_epi32(minmag);
		if(image->depth==8)
		{
			mminmag=_mm512_srli_epi32(mminmag, 2);
			minmag>>=2;
		}
		if(fwd)	//AVX-512	enc
		{
			DList list[4];
			dlist_init(list+0, 1, 1024, 0);
			dlist_init(list+1, 1, 1024, 0);
			dlist_init(list+2, 1, 1024, 0);
			dlist_init(list+3, 1, 1024, 0);
			gr_enc_init(ec+0, list+0);
			gr_enc_init(ec+1, list+1);
			gr_enc_init(ec+2, list+2);
			gr_enc_init(ec+3, list+3);
			ALIGN(32) const short *srcptr[]=
			{
				image->data+0*blocksize,
				image->data+1*blocksize,
				image->data+2*blocksize,
				image->data+3*blocksize,
			};
			for(int ky=0;ky<ceilqheight;++ky)
			{
				ALIGN(16) int *rows[]=
				{
					pixels+(((image->iw+16LL)*((ky-0LL)&1)+8)<<5),
					pixels+(((image->iw+16LL)*((ky-1LL)&1)+8)<<5),
				};
				for(int kx=0;kx<image->iw;++kx)
				{
					__m512i
						NW	=_mm512_load_si512((__m512i*)rows[1+0]-1*2+0),
						N	=_mm512_load_si512((__m512i*)rows[1+0]+0*2+0),
						W	=_mm512_load_si512((__m512i*)rows[0+0]-1*2+0),
						eNEEE	=_mm512_load_si512((__m512i*)rows[1+0]+3*2+1),
						eW	=_mm512_load_si512((__m512i*)rows[0+0]-1*2+1);
					int	*curr	=rows[0+0]+0*32;
					__m512i mcurr, mpred, mmag;
					calc_ctx_512(&N, &W, &NW, &eW, &mminmag, &mpred, &mmag);
					ALIGN(64) int mag[16];
					_mm512_store_si512((__m512i*)mag, mmag);
					int workmask=srcptr[3]<endptr;
					//if(!workmask)
					//	printf("");
					switch(image->nch)	//try separate switch instead of CMOVs
					{
					case 4:
						curr[3*4+0]=		srcptr[0][3];
						curr[3*4+1]=		srcptr[1][3];
						curr[3*4+2]=		srcptr[2][3];
						curr[3*4+3]=workmask?	srcptr[3][3]:0;
					case 3:
						curr[2*4+0]=		srcptr[0][2];
						curr[2*4+1]=		srcptr[1][2];
						curr[2*4+2]=		srcptr[2][2];
						curr[2*4+3]=workmask?	srcptr[3][2]:0;
					case 2:
						curr[1*4+0]=		srcptr[0][1];
						curr[1*4+1]=		srcptr[1][1];
						curr[1*4+2]=		srcptr[2][1];
						curr[1*4+3]=workmask?	srcptr[3][1]:0;
					case 1:
						curr[0*4+0]=		srcptr[0][0];
						curr[0*4+1]=		srcptr[1][0];
						curr[0*4+2]=		srcptr[2][0];
						curr[0*4+3]=workmask?	srcptr[3][0]:0;
						break;
					}
					if(image->nch>=3)
					{
						//JPEG2000 RCT
						__m128i mr=_mm_load_si128((__m128i*)curr+0);
						__m128i mg=_mm_load_si128((__m128i*)curr+1);
						__m128i mb=_mm_load_si128((__m128i*)curr+2);
						mr=_mm_sub_epi32(mr, mg);
						mb=_mm_sub_epi32(mb, mg);
						mg=_mm_add_epi32(mg, _mm_srai_epi32(_mm_add_epi32(mr, mb), 2));
						_mm_store_si128((__m128i*)curr+0, mr);
						_mm_store_si128((__m128i*)curr+1, mg);
						_mm_store_si128((__m128i*)curr+2, mb);
					}
					mcurr=_mm512_load_si512((__m128i*)curr+0);
					__m512i delta=_mm512_sub_epi32(mcurr, mpred);
					pack_sign_512(&delta);
					//_mm512_store_si512((__m512i*)curr+0, mcurr);
					_mm512_store_si512((__m512i*)curr+1, delta);
					if(workmask)
					{
						switch(image->nch)
						{
						case 4:
							gr_enc_POT(ec+3, curr[16+3*4+3], mag[3*4+3]);
							gr_enc_POT(ec+2, curr[16+3*4+2], mag[3*4+2]);
							gr_enc_POT(ec+1, curr[16+3*4+1], mag[3*4+1]);
							gr_enc_POT(ec+0, curr[16+3*4+0], mag[3*4+0]);
						case 3:
							gr_enc_POT(ec+3, curr[16+2*4+3], mag[2*4+3]);
							gr_enc_POT(ec+2, curr[16+2*4+2], mag[2*4+2]);
							gr_enc_POT(ec+1, curr[16+2*4+1], mag[2*4+1]);
							gr_enc_POT(ec+0, curr[16+2*4+0], mag[2*4+0]);
						case 2:
							gr_enc_POT(ec+3, curr[16+1*4+3], mag[1*4+3]);
							gr_enc_POT(ec+2, curr[16+1*4+2], mag[1*4+2]);
							gr_enc_POT(ec+1, curr[16+1*4+1], mag[1*4+1]);
							gr_enc_POT(ec+0, curr[16+1*4+0], mag[1*4+0]);
						case 1:
							gr_enc_POT(ec+3, curr[16+0*4+3], mag[0*4+3]);
							gr_enc_POT(ec+2, curr[16+0*4+2], mag[0*4+2]);
							gr_enc_POT(ec+1, curr[16+0*4+1], mag[0*4+1]);
							gr_enc_POT(ec+0, curr[16+0*4+0], mag[0*4+0]);
							break;
						}
					}
					else
					{
						switch(image->nch)
						{
						case 4:
							gr_enc_POT(ec+2, curr[16+3*4+2], mag[3*4+2]);
							gr_enc_POT(ec+1, curr[16+3*4+1], mag[3*4+1]);
							gr_enc_POT(ec+0, curr[16+3*4+0], mag[3*4+0]);
						case 3:
							gr_enc_POT(ec+2, curr[16+2*4+2], mag[2*4+2]);
							gr_enc_POT(ec+1, curr[16+2*4+1], mag[2*4+1]);
							gr_enc_POT(ec+0, curr[16+2*4+0], mag[2*4+0]);
						case 2:
							gr_enc_POT(ec+2, curr[16+1*4+2], mag[1*4+2]);
							gr_enc_POT(ec+1, curr[16+1*4+1], mag[1*4+1]);
							gr_enc_POT(ec+0, curr[16+1*4+0], mag[1*4+0]);
						case 1:
							gr_enc_POT(ec+2, curr[16+0*4+2], mag[0*4+2]);
							gr_enc_POT(ec+1, curr[16+0*4+1], mag[0*4+1]);
							gr_enc_POT(ec+0, curr[16+0*4+0], mag[0*4+0]);
							break;
						}
					}
					calc_update_512(&eW, &eNEEE, &delta);
					_mm512_store_si512((__m512i*)curr+1, delta);
				
					__m128i mrows=_mm_load_si128((__m128i*)rows);
					mrows=_mm_add_epi64(mrows, _mm_set1_epi64x(sizeof(__m512i[2])));
					_mm_store_si128((__m128i*)rows, mrows);

					__m256i mptr=_mm256_load_si256((__m256i*)srcptr);
					mptr=_mm256_add_epi64(mptr, stride);
					_mm256_store_si256((__m256i*)srcptr, mptr);
				}
			}
			gr_enc_flush(ec+0);
			gr_enc_flush(ec+1);
			gr_enc_flush(ec+2);
			gr_enc_flush(ec+3);
			size_t dststart=
			array_append(data, &list[0].nobj, 1, 4, 1, 0, 0);
			array_append(data, &list[1].nobj, 1, 4, 1, 0, 0);
			array_append(data, &list[2].nobj, 1, 4, 1, 0, 0);
			dlist_appendtoarray(list+0, data);
			dlist_appendtoarray(list+1, data);
			dlist_appendtoarray(list+2, data);
			dlist_appendtoarray(list+3, data);
			if(loud)
			{
				ptrdiff_t csize=data[0]->count-dststart;
				t0=time_sec()-t0;
				printf("%14td/%14td = %10.6lf%%  CR %lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
				printf("E (AVX-512)  %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
			}
			dlist_clear(list+0);
			dlist_clear(list+1);
			dlist_clear(list+2);
			dlist_clear(list+3);
		}
		else	//AVX-512	dec
		{
			unsigned listsizes[3];
			memcpy(listsizes, cbuf, sizeof(listsizes));
			cbuf+=sizeof(listsizes);
			clen-=sizeof(listsizes);
			gr_dec_init(ec+0, cbuf, cbuf+listsizes[0]);	cbuf+=listsizes[0], clen-=listsizes[0];
			gr_dec_init(ec+1, cbuf, cbuf+listsizes[1]);	cbuf+=listsizes[1], clen-=listsizes[1];
			gr_dec_init(ec+2, cbuf, cbuf+listsizes[2]);	cbuf+=listsizes[2], clen-=listsizes[2];
			gr_dec_init(ec+3, cbuf, cbuf+clen);
			ALIGN(32) short *dstptr[]=
			{
				dst->data+0*blocksize,
				dst->data+1*blocksize,
				dst->data+2*blocksize,
				dst->data+3*blocksize,
			};
			for(int ky=0;ky<ceilqheight;++ky)
			{
				ALIGN(16) int *rows[]=
				{
					pixels+(((image->iw+16LL)*((ky-0LL)&1)+8)<<5),
					pixels+(((image->iw+16LL)*((ky-1LL)&1)+8)<<5),
				};
				for(int kx=0;kx<image->iw;++kx)
				{
					__m512i
						NW	=_mm512_load_si512((__m512i*)rows[1+0]-1*2+0),
						N	=_mm512_load_si512((__m512i*)rows[1+0]+0*2+0),
						W	=_mm512_load_si512((__m512i*)rows[0+0]-1*2+0),
						eNEEE	=_mm512_load_si512((__m512i*)rows[1+0]+3*2+1),
						eW	=_mm512_load_si512((__m512i*)rows[0+0]-1*2+1);
					int	*curr	=rows[0+0]+0*32;
					__m512i mcurr, mpred, mmag;
					calc_ctx_512(&N, &W, &NW, &eW, &mminmag, &mpred, &mmag);
					ALIGN(64) int mag[16];
					_mm512_store_si512((__m512i*)mag, mmag);
					int workmask=dstptr[3]<endptr;
					if(workmask)
					{
						switch(image->nch)
						{
						case 4:
							curr[16+3*4+3]=gr_dec_POT(ec+3, mag[3*4+3]);
							curr[16+3*4+2]=gr_dec_POT(ec+2, mag[3*4+2]);
							curr[16+3*4+1]=gr_dec_POT(ec+1, mag[3*4+1]);
							curr[16+3*4+0]=gr_dec_POT(ec+0, mag[3*4+0]);
						case 3:
							curr[16+2*4+3]=gr_dec_POT(ec+3, mag[2*4+3]);
							curr[16+2*4+2]=gr_dec_POT(ec+2, mag[2*4+2]);
							curr[16+2*4+1]=gr_dec_POT(ec+1, mag[2*4+1]);
							curr[16+2*4+0]=gr_dec_POT(ec+0, mag[2*4+0]);
						case 2:
							curr[16+1*4+3]=gr_dec_POT(ec+3, mag[1*4+3]);
							curr[16+1*4+2]=gr_dec_POT(ec+2, mag[1*4+2]);
							curr[16+1*4+1]=gr_dec_POT(ec+1, mag[1*4+1]);
							curr[16+1*4+0]=gr_dec_POT(ec+0, mag[1*4+0]);
						case 1:
							curr[16+0*4+3]=gr_dec_POT(ec+3, mag[0*4+3]);
							curr[16+0*4+2]=gr_dec_POT(ec+2, mag[0*4+2]);
							curr[16+0*4+1]=gr_dec_POT(ec+1, mag[0*4+1]);
							curr[16+0*4+0]=gr_dec_POT(ec+0, mag[0*4+0]);
							break;
						}
					}
					else
					{
						switch(image->nch)
						{
						case 4:
							curr[16+3*4+2]=gr_dec_POT(ec+2, mag[3*4+2]);
							curr[16+3*4+1]=gr_dec_POT(ec+1, mag[3*4+1]);
							curr[16+3*4+0]=gr_dec_POT(ec+0, mag[3*4+0]);
						case 3:
							curr[16+2*4+2]=gr_dec_POT(ec+2, mag[2*4+2]);
							curr[16+2*4+1]=gr_dec_POT(ec+1, mag[2*4+1]);
							curr[16+2*4+0]=gr_dec_POT(ec+0, mag[2*4+0]);
						case 2:
							curr[16+1*4+2]=gr_dec_POT(ec+2, mag[1*4+2]);
							curr[16+1*4+1]=gr_dec_POT(ec+1, mag[1*4+1]);
							curr[16+1*4+0]=gr_dec_POT(ec+0, mag[1*4+0]);
						case 1:
							curr[16+0*4+2]=gr_dec_POT(ec+2, mag[0*4+2]);
							curr[16+0*4+1]=gr_dec_POT(ec+1, mag[0*4+1]);
							curr[16+0*4+0]=gr_dec_POT(ec+0, mag[0*4+0]);
							break;
						}
					}
					__m512i delta=_mm512_load_si512((__m512i*)curr+1);
					mcurr=delta;
					unpack_sign_512(&mcurr);
					calc_update_512(&eW, &eNEEE, &delta);
					mcurr=_mm512_add_epi32(mcurr, mpred);
					_mm512_store_si512((__m512i*)curr+0, mcurr);
					_mm512_store_si512((__m512i*)curr+1, delta);

					ALIGN(64) int c2[16];
					_mm512_store_si512((__m512i*)c2, mcurr);
					if(image->nch>=3)
					{
						//JPEG2000 RCT
						__m128i mr=_mm_load_si128((__m128i*)c2+0);
						__m128i mg=_mm_load_si128((__m128i*)c2+1);
						__m128i mb=_mm_load_si128((__m128i*)c2+2);
						mg=_mm_sub_epi32(mg, _mm_srai_epi32(_mm_add_epi32(mr, mb), 2));
						mr=_mm_add_epi32(mr, mg);
						mb=_mm_add_epi32(mb, mg);
						_mm_store_si128((__m128i*)c2+0, mr);
						_mm_store_si128((__m128i*)c2+1, mg);
						_mm_store_si128((__m128i*)c2+2, mb);
					}
					switch(image->nch)
					{
					case 4:
				if(workmask)	dstptr[3][3]=(short)c2[3*4+3];
						dstptr[2][3]=(short)c2[3*4+2];
						dstptr[1][3]=(short)c2[3*4+1];
						dstptr[0][3]=(short)c2[3*4+0];
					case 3:
				if(workmask)	dstptr[3][2]=(short)c2[2*4+3];
						dstptr[2][2]=(short)c2[2*4+2];
						dstptr[1][2]=(short)c2[2*4+1];
						dstptr[0][2]=(short)c2[2*4+0];
					case 2:
				if(workmask)	dstptr[3][1]=(short)c2[1*4+3];
						dstptr[2][1]=(short)c2[1*4+2];
						dstptr[1][1]=(short)c2[1*4+1];
						dstptr[0][1]=(short)c2[1*4+0];
					case 1:
				if(workmask)	dstptr[3][0]=(short)c2[0*4+3];
						dstptr[2][0]=(short)c2[0*4+2];
						dstptr[1][0]=(short)c2[0*4+1];
						dstptr[0][0]=(short)c2[0*4+0];
						break;
					}
#ifdef ENABLE_GUIDE
					ptrdiff_t idx0=dstptr[0]-image->data;
					if(guide&&memcmp(dstptr[0], guide->data+idx0, sizeof(short)*image->nch))
					{
						short c0[4]={0}, c1[4]={0};
						memcpy(c0, guide->data+idx0, sizeof(short)*image->nch);
						memcpy(c1, dstptr[0], sizeof(short)*image->nch);
						c0[0]-=c0[1];
						c0[2]-=c0[1];
						c0[1]+=(c0[0]+c0[2])>>2;
						c1[0]-=c1[1];
						c1[2]-=c1[1];
						c1[1]+=(c1[0]+c1[2])>>2;
						LOG_ERROR("Guide error IDX %d/%d", idx0, image->nch*image->iw*image->ih);
						printf("");//
					}
#endif
					__m128i mrows=_mm_load_si128((__m128i*)rows);
					mrows=_mm_add_epi64(mrows, _mm_set1_epi64x(sizeof(__m512i[2])));
					_mm_store_si128((__m128i*)rows, mrows);

					__m256i mptr=_mm256_load_si256((__m256i*)dstptr);
					mptr=_mm256_add_epi64(mptr, stride);
					_mm256_store_si256((__m256i*)dstptr, mptr);
				}
			}
			if(loud)
			{
				t0=time_sec()-t0;
				printf("D (AVX-512)  %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
			}
		}
	}
	else
#endif
	{
		__m256i mminmag=_mm256_set1_epi32(minmag);
		if(image->depth==8)
		{
			mminmag=_mm256_srli_epi32(mminmag, 2);
			minmag>>=2;
		}
		if(fwd)	//AVX2		enc
		{
			DList list[4];
			dlist_init(list+0, 1, 1024, 0);
			dlist_init(list+1, 1, 1024, 0);
			dlist_init(list+2, 1, 1024, 0);
			dlist_init(list+3, 1, 1024, 0);
			gr_enc_init(ec+0, list+0);
			gr_enc_init(ec+1, list+1);
			gr_enc_init(ec+2, list+2);
			gr_enc_init(ec+3, list+3);
			ALIGN(32) const short *srcptr[]=
			{
				image->data+0*blocksize,
				image->data+1*blocksize,
				image->data+2*blocksize,
				image->data+3*blocksize,
			};
			for(int ky=0;ky<ceilqheight;++ky)
			{
				ALIGN(16) int *rows[]=
				{
					pixels+(((image->iw+16LL)*((ky-0LL)&1)+8)<<5),
					pixels+(((image->iw+16LL)*((ky-1LL)&1)+8)<<5),
				};
				for(int kx=0;kx<image->iw;++kx)
				{
					__m256i
						NW[]	={_mm256_load_si256((__m256i*)rows[1+0]-1*4+0+0), _mm256_load_si256((__m256i*)rows[1+0]-1*4+0+1)},
						N[]	={_mm256_load_si256((__m256i*)rows[1+0]+0*4+0+0), _mm256_load_si256((__m256i*)rows[1+0]+0*4+0+1)},
						W[]	={_mm256_load_si256((__m256i*)rows[0+0]-1*4+0+0), _mm256_load_si256((__m256i*)rows[0+0]-1*4+0+1)},
						eNEEE[]	={_mm256_load_si256((__m256i*)rows[1+0]+3*4+2+0), _mm256_load_si256((__m256i*)rows[1+0]+3*4+2+1)},
						eW[]	={_mm256_load_si256((__m256i*)rows[0+0]-1*4+2+0), _mm256_load_si256((__m256i*)rows[0+0]-1*4+2+1)};
					int	*curr	=rows[0+0]+0*32;
					__m256i mcurr[2], mpred[2], mmag[2];
					calc_ctx_256x2(N, W, NW, eW, &mminmag, mpred, mmag);
					ALIGN(32) int mag[16];
					_mm256_store_si256((__m256i*)mag+0, mmag[0]);
					_mm256_store_si256((__m256i*)mag+1, mmag[1]);
					int workmask=srcptr[3]<endptr;
					//if(ky==3&&kx==130)//
					//	printf("");
					switch(image->nch)
					{
					case 4:
						curr[3*4+0]=		srcptr[0][3];
						curr[3*4+1]=		srcptr[1][3];
						curr[3*4+2]=		srcptr[2][3];
						curr[3*4+3]=workmask?	srcptr[3][3]:0;
					case 3:
						curr[2*4+0]=		srcptr[0][2];
						curr[2*4+1]=		srcptr[1][2];
						curr[2*4+2]=		srcptr[2][2];
						curr[2*4+3]=workmask?	srcptr[3][2]:0;
					case 2:
						curr[1*4+0]=		srcptr[0][1];
						curr[1*4+1]=		srcptr[1][1];
						curr[1*4+2]=		srcptr[2][1];
						curr[1*4+3]=workmask?	srcptr[3][1]:0;
					case 1:
						curr[0*4+0]=		srcptr[0][0];
						curr[0*4+1]=		srcptr[1][0];
						curr[0*4+2]=		srcptr[2][0];
						curr[0*4+3]=workmask?	srcptr[3][0]:0;
						break;
					}
					if(image->nch>=3)
					{
						//JPEG2000 RCT
						__m128i mr=_mm_load_si128((__m128i*)curr+0);
						__m128i mg=_mm_load_si128((__m128i*)curr+1);
						__m128i mb=_mm_load_si128((__m128i*)curr+2);
						mr=_mm_sub_epi32(mr, mg);
						mb=_mm_sub_epi32(mb, mg);
						mg=_mm_add_epi32(mg, _mm_srai_epi32(_mm_add_epi32(mr, mb), 2));
						_mm_store_si128((__m128i*)curr+0, mr);
						_mm_store_si128((__m128i*)curr+1, mg);
						_mm_store_si128((__m128i*)curr+2, mb);
					}
					mcurr[0]=_mm256_load_si256((__m256i*)curr+0+0);
					mcurr[1]=_mm256_load_si256((__m256i*)curr+0+1);
					__m256i delta[]=
					{
						_mm256_sub_epi32(mcurr[0], mpred[0]),
						_mm256_sub_epi32(mcurr[1], mpred[1]),
					};
					pack_sign_256x2(delta);
					_mm256_store_si256((__m256i*)curr+2+0, delta[0]);
					_mm256_store_si256((__m256i*)curr+2+1, delta[1]);
					if(workmask)
					{
						switch(image->nch)
						{
						case 4:	gr_enc_POT(ec+3, curr[16+3*4+3], FLOOR_LOG2(mag[3*4+3]));
							gr_enc_POT(ec+2, curr[16+3*4+2], FLOOR_LOG2(mag[3*4+2]));
							gr_enc_POT(ec+1, curr[16+3*4+1], FLOOR_LOG2(mag[3*4+1]));
							gr_enc_POT(ec+0, curr[16+3*4+0], FLOOR_LOG2(mag[3*4+0]));
						case 3:	gr_enc_POT(ec+3, curr[16+2*4+3], FLOOR_LOG2(mag[2*4+3]));
							gr_enc_POT(ec+2, curr[16+2*4+2], FLOOR_LOG2(mag[2*4+2]));
							gr_enc_POT(ec+1, curr[16+2*4+1], FLOOR_LOG2(mag[2*4+1]));
							gr_enc_POT(ec+0, curr[16+2*4+0], FLOOR_LOG2(mag[2*4+0]));
						case 2:	gr_enc_POT(ec+3, curr[16+1*4+3], FLOOR_LOG2(mag[1*4+3]));
							gr_enc_POT(ec+2, curr[16+1*4+2], FLOOR_LOG2(mag[1*4+2]));
							gr_enc_POT(ec+1, curr[16+1*4+1], FLOOR_LOG2(mag[1*4+1]));
							gr_enc_POT(ec+0, curr[16+1*4+0], FLOOR_LOG2(mag[1*4+0]));
						case 1:	gr_enc_POT(ec+3, curr[16+0*4+3], FLOOR_LOG2(mag[0*4+3]));
							gr_enc_POT(ec+2, curr[16+0*4+2], FLOOR_LOG2(mag[0*4+2]));
							gr_enc_POT(ec+1, curr[16+0*4+1], FLOOR_LOG2(mag[0*4+1]));
							gr_enc_POT(ec+0, curr[16+0*4+0], FLOOR_LOG2(mag[0*4+0]));
							break;
						}
					}
					else
					{
						switch(image->nch)
						{
						case 4:	gr_enc_POT(ec+2, curr[16+3*4+2], FLOOR_LOG2(mag[3*4+2]));
							gr_enc_POT(ec+1, curr[16+3*4+1], FLOOR_LOG2(mag[3*4+1]));
							gr_enc_POT(ec+0, curr[16+3*4+0], FLOOR_LOG2(mag[3*4+0]));
						case 3:	gr_enc_POT(ec+2, curr[16+2*4+2], FLOOR_LOG2(mag[2*4+2]));
							gr_enc_POT(ec+1, curr[16+2*4+1], FLOOR_LOG2(mag[2*4+1]));
							gr_enc_POT(ec+0, curr[16+2*4+0], FLOOR_LOG2(mag[2*4+0]));
						case 2:	gr_enc_POT(ec+2, curr[16+1*4+2], FLOOR_LOG2(mag[1*4+2]));
							gr_enc_POT(ec+1, curr[16+1*4+1], FLOOR_LOG2(mag[1*4+1]));
							gr_enc_POT(ec+0, curr[16+1*4+0], FLOOR_LOG2(mag[1*4+0]));
						case 1:	gr_enc_POT(ec+2, curr[16+0*4+2], FLOOR_LOG2(mag[0*4+2]));
							gr_enc_POT(ec+1, curr[16+0*4+1], FLOOR_LOG2(mag[0*4+1]));
							gr_enc_POT(ec+0, curr[16+0*4+0], FLOOR_LOG2(mag[0*4+0]));
							break;
						}
					}
					calc_update_256x2(eW, eNEEE, delta);
					_mm256_store_si256((__m256i*)curr+0+2, delta[0]);
					_mm256_store_si256((__m256i*)curr+1+2, delta[1]);

					__m128i mrows=_mm_load_si128((__m128i*)rows);
					mrows=_mm_add_epi64(mrows, _mm_set1_epi64x(sizeof(__m256i[4])));
					_mm_store_si128((__m128i*)rows, mrows);

					__m256i mptr=_mm256_load_si256((__m256i*)srcptr);
					mptr=_mm256_add_epi64(mptr, stride);
					_mm256_store_si256((__m256i*)srcptr, mptr);
				}
			}
			gr_enc_flush(ec+0);
			gr_enc_flush(ec+1);
			gr_enc_flush(ec+2);
			gr_enc_flush(ec+3);
			size_t dststart=
			array_append(data, &list[0].nobj, 1, 4, 1, 0, 0);
			array_append(data, &list[1].nobj, 1, 4, 1, 0, 0);
			array_append(data, &list[2].nobj, 1, 4, 1, 0, 0);
			dlist_appendtoarray(list+0, data);
			dlist_appendtoarray(list+1, data);
			dlist_appendtoarray(list+2, data);
			dlist_appendtoarray(list+3, data);
			if(loud)
			{
				ptrdiff_t csize=data[0]->count-dststart;
				t0=time_sec()-t0;
				printf("%14td/%14td = %10.6lf%%  CR %lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
				printf("E  %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
			}
			dlist_clear(list+0);
			dlist_clear(list+1);
			dlist_clear(list+2);
			dlist_clear(list+3);
		}
		else	//AVX2		dec
		{
			unsigned listsizes[3];
			memcpy(listsizes, cbuf, sizeof(listsizes));
			cbuf+=sizeof(listsizes);
			clen-=sizeof(listsizes);
			gr_dec_init(ec+0, cbuf, cbuf+listsizes[0]);	cbuf+=listsizes[0], clen-=listsizes[0];
			gr_dec_init(ec+1, cbuf, cbuf+listsizes[1]);	cbuf+=listsizes[1], clen-=listsizes[1];
			gr_dec_init(ec+2, cbuf, cbuf+listsizes[2]);	cbuf+=listsizes[2], clen-=listsizes[2];
			gr_dec_init(ec+3, cbuf, cbuf+clen);
			ALIGN(32) short *dstptr[]=
			{
				dst->data+0*blocksize,
				dst->data+1*blocksize,
				dst->data+2*blocksize,
				dst->data+3*blocksize,
			};
			for(int ky=0;ky<ceilqheight;++ky)
			{
				ALIGN(32) int *rows[]=
				{
					pixels+(((image->iw+16LL)*((ky-0LL)&1)+8)<<5),
					pixels+(((image->iw+16LL)*((ky-1LL)&1)+8)<<5),
				};
				for(int kx=0;kx<image->iw;++kx)
				{
					__m256i
						NW[]	={_mm256_load_si256((__m256i*)rows[1+0]-1*4+0+0), _mm256_load_si256((__m256i*)rows[1+0]-1*4+0+1)},
						N[]	={_mm256_load_si256((__m256i*)rows[1+0]+0*4+0+0), _mm256_load_si256((__m256i*)rows[1+0]+0*4+0+1)},
						W[]	={_mm256_load_si256((__m256i*)rows[0+0]-1*4+0+0), _mm256_load_si256((__m256i*)rows[0+0]-1*4+0+1)},
						eNEEE[]	={_mm256_load_si256((__m256i*)rows[1+0]+3*4+2+0), _mm256_load_si256((__m256i*)rows[1+0]+3*4+2+1)},
						eW[]	={_mm256_load_si256((__m256i*)rows[0+0]-1*4+2+0), _mm256_load_si256((__m256i*)rows[0+0]-1*4+2+1)};
					int	*curr	=rows[0+0]+0*32;
					__m256i mcurr[2], mpred[2], mmag[2];
					//if(ky==3&&kx==130)//
					//	printf("");
					calc_ctx_256x2(N, W, NW, eW, &mminmag, mpred, mmag);
					ALIGN(32) int mag[16];
					_mm256_store_si256((__m256i*)mag+0, mmag[0]);
					_mm256_store_si256((__m256i*)mag+1, mmag[1]);
					int workmask=dstptr[3]<endptr;
					if(workmask)
					{
						switch(image->nch)
						{
						case 4:	curr[16+3*4+3]=gr_dec_POT(ec+3, FLOOR_LOG2(mag[3*4+3]));
							curr[16+3*4+2]=gr_dec_POT(ec+2, FLOOR_LOG2(mag[3*4+2]));
							curr[16+3*4+1]=gr_dec_POT(ec+1, FLOOR_LOG2(mag[3*4+1]));
							curr[16+3*4+0]=gr_dec_POT(ec+0, FLOOR_LOG2(mag[3*4+0]));
						case 3:	curr[16+2*4+3]=gr_dec_POT(ec+3, FLOOR_LOG2(mag[2*4+3]));
							curr[16+2*4+2]=gr_dec_POT(ec+2, FLOOR_LOG2(mag[2*4+2]));
							curr[16+2*4+1]=gr_dec_POT(ec+1, FLOOR_LOG2(mag[2*4+1]));
							curr[16+2*4+0]=gr_dec_POT(ec+0, FLOOR_LOG2(mag[2*4+0]));
						case 2:	curr[16+1*4+3]=gr_dec_POT(ec+3, FLOOR_LOG2(mag[1*4+3]));
							curr[16+1*4+2]=gr_dec_POT(ec+2, FLOOR_LOG2(mag[1*4+2]));
							curr[16+1*4+1]=gr_dec_POT(ec+1, FLOOR_LOG2(mag[1*4+1]));
							curr[16+1*4+0]=gr_dec_POT(ec+0, FLOOR_LOG2(mag[1*4+0]));
						case 1:	curr[16+0*4+3]=gr_dec_POT(ec+3, FLOOR_LOG2(mag[0*4+3]));
							curr[16+0*4+2]=gr_dec_POT(ec+2, FLOOR_LOG2(mag[0*4+2]));
							curr[16+0*4+1]=gr_dec_POT(ec+1, FLOOR_LOG2(mag[0*4+1]));
							curr[16+0*4+0]=gr_dec_POT(ec+0, FLOOR_LOG2(mag[0*4+0]));
							break;
						}
					}
					else
					{
						switch(image->nch)
						{
						case 4:	curr[16+3*4+2]=gr_dec_POT(ec+2, FLOOR_LOG2(mag[3*4+2]));
							curr[16+3*4+1]=gr_dec_POT(ec+1, FLOOR_LOG2(mag[3*4+1]));
							curr[16+3*4+0]=gr_dec_POT(ec+0, FLOOR_LOG2(mag[3*4+0]));
						case 3:	curr[16+2*4+2]=gr_dec_POT(ec+2, FLOOR_LOG2(mag[2*4+2]));
							curr[16+2*4+1]=gr_dec_POT(ec+1, FLOOR_LOG2(mag[2*4+1]));
							curr[16+2*4+0]=gr_dec_POT(ec+0, FLOOR_LOG2(mag[2*4+0]));
						case 2:	curr[16+1*4+2]=gr_dec_POT(ec+2, FLOOR_LOG2(mag[1*4+2]));
							curr[16+1*4+1]=gr_dec_POT(ec+1, FLOOR_LOG2(mag[1*4+1]));
							curr[16+1*4+0]=gr_dec_POT(ec+0, FLOOR_LOG2(mag[1*4+0]));
						case 1:	curr[16+0*4+2]=gr_dec_POT(ec+2, FLOOR_LOG2(mag[0*4+2]));
							curr[16+0*4+1]=gr_dec_POT(ec+1, FLOOR_LOG2(mag[0*4+1]));
							curr[16+0*4+0]=gr_dec_POT(ec+0, FLOOR_LOG2(mag[0*4+0]));
							break;
						}
					}
					__m256i delta[]=
					{
						_mm256_load_si256((__m256i*)curr+0+2),
						_mm256_load_si256((__m256i*)curr+1+2),
					};
					mcurr[0]=delta[0];
					mcurr[1]=delta[1];
					unpack_sign_256x2(mcurr);
					calc_update_256x2(eW, eNEEE, delta);
					mcurr[0]=_mm256_add_epi32(mcurr[0], mpred[0]);
					mcurr[1]=_mm256_add_epi32(mcurr[1], mpred[1]);
					_mm256_store_si256((__m256i*)curr+0+0, mcurr[0]);
					_mm256_store_si256((__m256i*)curr+1+0, mcurr[1]);
					_mm256_store_si256((__m256i*)curr+0+2, delta[0]);
					_mm256_store_si256((__m256i*)curr+1+2, delta[1]);

					ALIGN(32) int c2[16];
					_mm256_store_si256((__m256i*)c2+0, mcurr[0]);
					_mm256_store_si256((__m256i*)c2+1, mcurr[1]);
					if(image->nch>=3)
					{
						//JPEG2000 RCT
						__m128i mr=_mm_load_si128((__m128i*)c2+0);
						__m128i mg=_mm_load_si128((__m128i*)c2+1);
						__m128i mb=_mm_load_si128((__m128i*)c2+2);
						mg=_mm_sub_epi32(mg, _mm_srai_epi32(_mm_add_epi32(mr, mb), 2));
						mr=_mm_add_epi32(mr, mg);
						mb=_mm_add_epi32(mb, mg);
						_mm_store_si128((__m128i*)c2+0, mr);
						_mm_store_si128((__m128i*)c2+1, mg);
						_mm_store_si128((__m128i*)c2+2, mb);
					}
					switch(image->nch)
					{
					case 4:
				if(workmask)	dstptr[3][3]=(short)c2[3*4+3];
						dstptr[2][3]=(short)c2[3*4+2];
						dstptr[1][3]=(short)c2[3*4+1];
						dstptr[0][3]=(short)c2[3*4+0];
					case 3:
				if(workmask)	dstptr[3][2]=(short)c2[2*4+3];
						dstptr[2][2]=(short)c2[2*4+2];
						dstptr[1][2]=(short)c2[2*4+1];
						dstptr[0][2]=(short)c2[2*4+0];
					case 2:
				if(workmask)	dstptr[3][1]=(short)c2[1*4+3];
						dstptr[2][1]=(short)c2[1*4+2];
						dstptr[1][1]=(short)c2[1*4+1];
						dstptr[0][1]=(short)c2[1*4+0];
					case 1:
				if(workmask)	dstptr[3][0]=(short)c2[0*4+3];
						dstptr[2][0]=(short)c2[0*4+2];
						dstptr[1][0]=(short)c2[0*4+1];
						dstptr[0][0]=(short)c2[0*4+0];
						break;
					}
#ifdef ENABLE_GUIDE
					ptrdiff_t idx0=dstptr[0]-image->data;
					if(guide&&memcmp(dstptr[0], guide->data+idx0, sizeof(short)*image->nch))
					{
						short c0[4]={0}, c1[4]={0};
						memcpy(c0, guide->data+idx0, sizeof(short)*image->nch);
						memcpy(c1, dstptr[0], sizeof(short)*image->nch);
						c0[0]-=c0[1];
						c0[2]-=c0[1];
						c0[1]+=(c0[0]+c0[2])>>2;
						c1[0]-=c1[1];
						c1[2]-=c1[1];
						c1[1]+=(c1[0]+c1[2])>>2;
						LOG_ERROR("Guide error IDX %d/%d", idx0, image->nch*image->iw*image->ih);
						printf("");//
					}
#endif
					__m128i mrows=_mm_load_si128((__m128i*)rows);
					mrows=_mm_add_epi64(mrows, _mm_set1_epi64x(sizeof(__m256i[4])));
					_mm_store_si128((__m128i*)rows, mrows);
					
					__m256i mptr=_mm256_load_si256((__m256i*)dstptr);
					mptr=_mm256_add_epi64(mptr, stride);
					_mm256_store_si256((__m256i*)dstptr, mptr);
				}
			}
			if(loud)
			{
				t0=time_sec()-t0;
				printf("D  %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
			}
		}
	}

	_mm_free(pixels);
	return 0;
}