#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define PROFILER 1
	#define N_CODER
//	#define DISABLE_RCT


#ifdef PROFILER
#define CHECKPOINTLIST\
	CHECKPOINT(INIT)\
	CHECKPOINT(RCT)\
	CHECKPOINT(PRED)\
	CHECKPOINT(DUMMY)\
	CHECKPOINT(EC)\
	CHECKPOINT(CDF)\
	CHECKPOINT(FINISH)
#endif
#include"ac.h"
#include"profiler.h"
#define QLEVELS 15
#define RATE_SR 8
#ifdef N_CODER
#define EC_IDX_MASK ~0
#else
#define EC_IDX_MASK 0
#endif
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
static void get_ctx(__m128i const *nb, short *pred, int *ctx)//nb={N, W, NW}, 16-bit
{
	__m128i hmasks[]=
	{
		_mm_set1_epi16(0x5555),
		_mm_set1_epi16(0x3333),
		_mm_set1_epi16(0x0F0F),
	};
	__m128i qhalf=_mm_set1_epi16(QLEVELS>>1);
	__m128i qmax=_mm_set1_epi16(QLEVELS-1);
	__m128i factor=_mm_set_epi16(QLEVELS, QLEVELS, QLEVELS, QLEVELS, 1, 1, 1, 1);
	__m128i swap64=_mm_set_epi8(
		 7,  6,  5,  4,  3,  2,  1,  0, 15, 14, 13, 12, 11, 10,  9,  8
	);
	__m128i cvt32=_mm_set_epi8(
		-1, -1,  7,  6, -1, -1,  5,  4, -1, -1,  3,  2, -1, -1,  1,  0
	);
	__m128i offset=_mm_set_epi32(
		256*QLEVELS*QLEVELS*3,
		256*QLEVELS*QLEVELS*2,
		256*QLEVELS*QLEVELS*1,
		256*QLEVELS*QLEVELS*0
	);

	__m128i vmin=_mm_min_epi16(nb[0], nb[1]);
	__m128i vmax=_mm_max_epi16(nb[0], nb[1]);
	__m128i mpred=_mm_add_epi16(nb[0], nb[1]);
	mpred=_mm_sub_epi16(mpred, nb[2]);
	mpred=_mm_min_epi16(mpred, vmax);
	mpred=_mm_max_epi16(mpred, vmin);
	_mm_store_si128((__m128i*)pred, mpred);

	//get ctx
	__m128i errors=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(nb[0]), _mm_castsi128_ps(nb[1]), _MM_SHUFFLE(3, 2, 3, 2)));//hi halves contain errors
	__m128i negmask=_mm_cmplt_epi16(errors, _mm_setzero_si128());
	errors=_mm_abs_epi16(errors);//remove sign bit (8 -> 7-bit)

	errors=_mm_or_si128(errors, _mm_srli_epi16(errors, 1));//set LSBs
	errors=_mm_or_si128(errors, _mm_srli_epi16(errors, 2));
	errors=_mm_or_si128(errors, _mm_srli_epi16(errors, 4));

	errors=_mm_sub_epi16(errors, _mm_and_si128(_mm_srli_epi16(errors, 1), hmasks[0]));//7-bit hamming weight
	errors=_mm_add_epi16(_mm_and_si128(errors, hmasks[1]), _mm_and_si128(_mm_srli_epi16(errors, 2), hmasks[1]));
	errors=_mm_add_epi16(_mm_and_si128(errors, hmasks[2]), _mm_and_si128(_mm_srli_epi16(errors, 4), hmasks[2]));

	errors=_mm_xor_si128(errors, negmask);//negate if was negative
	errors=_mm_sub_epi16(errors, negmask);
	errors=_mm_add_epi16(errors, qhalf);//add half
	errors=_mm_min_epi16(errors, qmax);//clamp
	errors=_mm_max_epi16(errors, _mm_setzero_si128());

	errors=_mm_mullo_epi16(errors, factor);
	errors=_mm_add_epi16(errors, _mm_shuffle_epi8(errors, swap64));

	errors=_mm_mullo_epi16(errors, _mm_set1_epi16(256));//5 cycles
	errors=_mm_shuffle_epi8(errors, cvt32);
	//errors=_mm_cvtepi16_epi32(errors);
	//errors=_mm_mullo_epi32(errors, _mm_set1_epi32(256));//10 cycles

	errors=_mm_add_epi32(errors, offset);
	_mm_store_si128((__m128i*)ctx, errors);
}
static void update_CDFs(unsigned char *val, unsigned short *stats, int *ctx)
{
	for(int kc=0;kc<3;++kc)
	{
		unsigned short *curr_CDF=stats+ctx[kc];
#if 0
		for(int ks=0;ks<256;++ks)
			curr_CDF[ks]+=((0xFF00&-(ks>val[kc]))+ks-curr_CDF[ks])>>RATE_SR;
#else
		__m256i sym=_mm256_set1_epi16(val[kc]);
		__m256i r0=_mm256_set_epi16(15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0);
		__m256i r1=_mm256_set_epi16(31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16);
		__m256i r2=_mm256_set_epi16(47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32);
		__m256i r3=_mm256_set_epi16(63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48);
		__m256i stride=_mm256_set1_epi16(64);
		__m256i amplitude=_mm256_set1_epi16(0xFF00);
		__m256i expandlo=_mm256_set_epi8(
			-1, -1,  7,  6, -1, -1,  5,  4, -1, -1,  3,  2, -1, -1,  1,  0,
			-1, -1,  7,  6, -1, -1,  5,  4, -1, -1,  3,  2, -1, -1,  1,  0
		);
		__m256i expandhi=_mm256_set_epi8(
			-1, -1, 15, 14, -1, -1, 13, 12, -1, -1, 11, 10, -1, -1,  9,  8,
			-1, -1, 15, 14, -1, -1, 13, 12, -1, -1, 11, 10, -1, -1,  9,  8
		);
		__m256i compactlo=_mm256_set_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1, 13, 12,  9,  8,  5,  4,  1,  0,
			-1, -1, -1, -1, -1, -1, -1, -1, 13, 12,  9,  8,  5,  4,  1,  0
		);
		__m256i compacthi=_mm256_set_epi8(
			13, 12,  9,  8,  5,  4,  1,  0, -1, -1, -1, -1, -1, -1, -1, -1,
			13, 12,  9,  8,  5,  4,  1,  0, -1, -1, -1, -1, -1, -1, -1, -1
		);
		for(int k=0;k<256;k+=64)
		{
			__m256i c0=_mm256_load_si256((__m256i*)curr_CDF+0);
			__m256i c1=_mm256_load_si256((__m256i*)curr_CDF+1);
			__m256i c2=_mm256_load_si256((__m256i*)curr_CDF+2);
			__m256i c3=_mm256_load_si256((__m256i*)curr_CDF+3);
			__m256i t0=_mm256_cmpgt_epi16(r0, sym);
			__m256i t1=_mm256_cmpgt_epi16(r1, sym);
			__m256i t2=_mm256_cmpgt_epi16(r2, sym);
			__m256i t3=_mm256_cmpgt_epi16(r3, sym);
			t0=_mm256_and_si256(amplitude, t0);
			t1=_mm256_and_si256(amplitude, t1);
			t2=_mm256_and_si256(amplitude, t2);
			t3=_mm256_and_si256(amplitude, t3);
			t0=_mm256_add_epi16(t0, r0);
			t1=_mm256_add_epi16(t1, r1);
			t2=_mm256_add_epi16(t2, r2);
			t3=_mm256_add_epi16(t3, r3);

			//need 17 bits for this (including sign bit)
			__m256i tlo, thi, clo, chi;
			tlo=_mm256_shuffle_epi8(t0, expandlo);
			thi=_mm256_shuffle_epi8(t0, expandhi);
			clo=_mm256_shuffle_epi8(c0, expandlo);
			chi=_mm256_shuffle_epi8(c0, expandhi);
			clo=_mm256_sub_epi32(tlo, clo);
			chi=_mm256_sub_epi32(thi, chi);
			clo=_mm256_srai_epi32(clo, RATE_SR);
			chi=_mm256_srai_epi32(chi, RATE_SR);
			clo=_mm256_shuffle_epi8(clo, compactlo);
			chi=_mm256_shuffle_epi8(chi, compacthi);
			t0=_mm256_or_si256(clo, chi);
			
			tlo=_mm256_shuffle_epi8(t1, expandlo);
			thi=_mm256_shuffle_epi8(t1, expandhi);
			clo=_mm256_shuffle_epi8(c1, expandlo);
			chi=_mm256_shuffle_epi8(c1, expandhi);
			clo=_mm256_sub_epi32(tlo, clo);
			chi=_mm256_sub_epi32(thi, chi);
			clo=_mm256_srai_epi32(clo, RATE_SR);
			chi=_mm256_srai_epi32(chi, RATE_SR);
			clo=_mm256_shuffle_epi8(clo, compactlo);
			chi=_mm256_shuffle_epi8(chi, compacthi);
			t1=_mm256_or_si256(clo, chi);
			
			tlo=_mm256_shuffle_epi8(t2, expandlo);
			thi=_mm256_shuffle_epi8(t2, expandhi);
			clo=_mm256_shuffle_epi8(c2, expandlo);
			chi=_mm256_shuffle_epi8(c2, expandhi);
			clo=_mm256_sub_epi32(tlo, clo);
			chi=_mm256_sub_epi32(thi, chi);
			clo=_mm256_srai_epi32(clo, RATE_SR);
			chi=_mm256_srai_epi32(chi, RATE_SR);
			clo=_mm256_shuffle_epi8(clo, compactlo);
			chi=_mm256_shuffle_epi8(chi, compacthi);
			t2=_mm256_or_si256(clo, chi);
			
			tlo=_mm256_shuffle_epi8(t3, expandlo);
			thi=_mm256_shuffle_epi8(t3, expandhi);
			clo=_mm256_shuffle_epi8(c3, expandlo);
			chi=_mm256_shuffle_epi8(c3, expandhi);
			clo=_mm256_sub_epi32(tlo, clo);
			chi=_mm256_sub_epi32(thi, chi);
			clo=_mm256_srai_epi32(clo, RATE_SR);
			chi=_mm256_srai_epi32(chi, RATE_SR);
			clo=_mm256_shuffle_epi8(clo, compactlo);
			chi=_mm256_shuffle_epi8(chi, compacthi);
			t3=_mm256_or_si256(clo, chi);

			c0=_mm256_add_epi16(c0, t0);
			c1=_mm256_add_epi16(c1, t1);
			c2=_mm256_add_epi16(c2, t2);
			c3=_mm256_add_epi16(c3, t3);
			_mm256_store_si256((__m256i*)curr_CDF+0, c0);
			_mm256_store_si256((__m256i*)curr_CDF+1, c1);
			_mm256_store_si256((__m256i*)curr_CDF+2, c2);
			_mm256_store_si256((__m256i*)curr_CDF+3, c3);
			curr_CDF+=64;
			r0=_mm256_add_epi16(r0, stride);
			r1=_mm256_add_epi16(r1, stride);
			r2=_mm256_add_epi16(r2, stride);
			r3=_mm256_add_epi16(r3, stride);
		}
#endif
		//curr_CDF=stats+ctx[kc];//
		//for(int ks=0;ks<255;++ks)
		//{
		//	if(curr_CDF[ks]>curr_CDF[ks+1])
		//		LOG_ERROR("");
		//}
	}
}
int f07_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
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
	DList list[3];
	dlist_init(list+0, 1, 1024, 0);
	dlist_init(list+1, 1, 1024, 0);
	dlist_init(list+2, 1, 1024, 0);
	ArithmeticCoder ec[3];
	unsigned short *stats=(unsigned short*)_mm_malloc(sizeof(short[256*QLEVELS*QLEVELS*4]), sizeof(__m256i));//(CDFSIZE+1) * nodes_in_tree * 4 channels max
	char *pixels=(char*)malloc((image->iw+4LL)*sizeof(char[2*4*2]));//2 padded rows * 4 channels max * {pixels, errors}
	if(!stats||!pixels)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int ks=0;ks<256;++ks)
		stats[ks]=(short)(ks<<8);
	memfill(stats+256, stats, sizeof(short[256*QLEVELS*QLEVELS*4])-sizeof(short[256]), sizeof(short[256]));
	memset(pixels, 0, (image->iw+4LL)*sizeof(char[2*4*2]));
	PROF(INIT);
	if(fwd)
	{
		ac_enc_init(ec+0, list+0);
		ac_enc_init(ec+1, list+1);
		ac_enc_init(ec+2, list+2);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			char *rows[2]=
			{
				pixels+(((image->iw+4LL)*(ky&1)+1)<<3),
				pixels+(((image->iw+4LL)*((ky-1)&1)+1)<<3),
			};
			for(int kx=0;kx<image->iw;++kx, idx+=3)
			{
				__m128i nb[]=
				{
					_mm_cvtepi8_epi16(_mm_loadu_si128((const __m128i*)(rows[1]+0))),//N
					_mm_cvtepi8_epi16(_mm_loadu_si128((const __m128i*)(rows[0]-8))),//W
					_mm_cvtepi8_epi16(_mm_loadu_si128((const __m128i*)(rows[1]-8))),//NW
				};
				ALIGN(16) short pred[8];
				ALIGN(16) int ctx[4];
				get_ctx(nb, pred, ctx);

				char *curr=rows[0];
				rows[0]+=8;
				rows[1]+=8;
				PROF(PRED);
				PROF(DUMMY);

				curr[0]=(char)image->data[idx+0];
				curr[1]=(char)image->data[idx+1];
				curr[2]=(char)image->data[idx+2];
#ifndef DISABLE_RCT
				curr[0]-=curr[1];
				curr[2]-=curr[1];
				curr[1]+=(curr[0]+curr[2])>>2;
#endif
				curr[4]=(char)(curr[0]-pred[0]);
				curr[5]=(char)(curr[1]-pred[1]);
				curr[6]=(char)(curr[2]-pred[2]);

				unsigned char val[]=
				{
					(unsigned char)(curr[4]+128),
					(unsigned char)(curr[5]+128),
					(unsigned char)(curr[6]+128),
				};
				PROF(RCT);
				//if(!ky&&kx==42)//
				//if(idx==84945-3)//
				//if(idx==2319)//
				//if(idx==9)//
				//if(idx==9621)//
				//if(idx==7299)//
				//	printf("");
#ifdef N_CODER
				ac_enc_packedCDF_8x3(ec, val, stats+ctx[0], stats+ctx[1], stats+ctx[2]);
#else
				ac_enc_packedCDF(ec+(0&EC_IDX_MASK), val[0], stats+ctx[0], 256);
				ac_enc_packedCDF(ec+(1&EC_IDX_MASK), val[1], stats+ctx[1], 256);
				ac_enc_packedCDF(ec+(2&EC_IDX_MASK), val[2], stats+ctx[2], 256);
#endif
				PROF(EC);

				update_CDFs(val, stats, ctx);
				PROF(CDF);
			}
		}
		ac_enc_flush(ec+0);
#ifdef N_CODER
		ac_enc_flush(ec+1);
		ac_enc_flush(ec+2);
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
		ac_dec_init(ec+0, ptr, ptr+bm[0]);	ptr+=bm[0];
#ifdef N_CODER
		ac_dec_init(ec+1, ptr, ptr+bm[1]);	ptr+=bm[1];
		ac_dec_init(ec+2, ptr, ptr+bm[2]);	ptr+=bm[2];
#endif
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			char *rows[2]=
			{
				pixels+(((image->iw+4LL)*(ky&1)+1)<<3),
				pixels+(((image->iw+4LL)*((ky-1)&1)+1)<<3),
			};
			for(int kx=0;kx<image->iw;++kx, idx+=3)
			{
				__m128i nb[]=
				{
					_mm_cvtepi8_epi16(_mm_loadu_si128((const __m128i*)(rows[1]+0))),//N
					_mm_cvtepi8_epi16(_mm_loadu_si128((const __m128i*)(rows[0]-8))),//W
					_mm_cvtepi8_epi16(_mm_loadu_si128((const __m128i*)(rows[1]-8))),//NW
				};
				ALIGN(16) short pred[8];
				ALIGN(16) int ctx[4];
				get_ctx(nb, pred, ctx);
					
				char *curr=rows[0];

				rows[0]+=8;
				rows[1]+=8;
				PROF(PRED);
				PROF(DUMMY);

				//if(idx==1095)//
				//if(idx==2319)//
				//if(idx==9)//
				//if(idx==9621)//
				//if(idx==7299)//
				//	printf("");
					
				unsigned char val[3];
#ifdef N_CODER
				ac_dec_packedCDF_8x3(ec, stats+ctx[0], stats+ctx[1], stats+ctx[2], val);
#else
				val[0]=ac_dec_packedCDF_POT(ec+(0&EC_IDX_MASK), stats+ctx[0], 8);
				val[1]=ac_dec_packedCDF_POT(ec+(1&EC_IDX_MASK), stats+ctx[1], 8);
				val[2]=ac_dec_packedCDF_POT(ec+(2&EC_IDX_MASK), stats+ctx[2], 8);
#endif
				PROF(EC);

				curr[4]=val[0]-128;
				curr[5]=val[1]-128;
				curr[6]=val[2]-128;

				char c2[]=
				{
					curr[0]=(char)(curr[4]+pred[0]),
					curr[1]=(char)(curr[5]+pred[1]),
					curr[2]=(char)(curr[6]+pred[2]),
				};
#ifndef DISABLE_RCT
				c2[1]-=(c2[0]+c2[2])>>2;
				c2[2]+=c2[1];
				c2[0]+=c2[1];
#endif
				short *rgb=dst->data+idx;
				rgb[0]=c2[0];
				rgb[1]=c2[1];
				rgb[2]=c2[2];
				PROF(RCT);
					
				update_CDFs(val, stats, ctx);
				PROF(CDF);
#ifdef ENABLE_GUIDE
				if(guide&&memcmp(rgb, guide->data+idx, image->nch*sizeof(short)))
				{
					short comp0[4]={0};
					memcpy(comp0, guide->data+idx, image->nch*sizeof(short));
					rgb[0]-=rgb[1];
					rgb[2]-=rgb[1];
					rgb[1]+=(rgb[0]+rgb[2])>>2;
					comp0[0]-=comp0[1];
					comp0[2]-=comp0[1];
					comp0[1]+=(comp0[0]+comp0[2])>>2;
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
			}
		}
	}
	if(loud)
	{
		t0=time_sec()-t0;
		if(fwd)
		{
			ptrdiff_t usize=image_getBMPsize(image);
			ptrdiff_t csize=list[0].nobj+list[1].nobj+list[2].nobj;
			printf("YUV %12zd %12zd %12zd\n",
				list[0].nobj,
				list[1].nobj,
				list[2].nobj
			);
			printf("csize %12td  %10.6lf%%  CR %8.6lf\n", csize, 100.*csize/usize, (double)usize/csize);
		}
		printf("F07  %c %15.6lf sec\n", 'D'+fwd, t0);
		prof_print();
	}
	dlist_clear(list+0);
	dlist_clear(list+1);
	dlist_clear(list+2);
	_mm_free(stats);
	free(pixels);
	return 1;
}