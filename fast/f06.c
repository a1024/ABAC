#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define PROFILER 1

	#define N_CODER
//	#define EC_USE_ARRAY
//	#define LINEAR_SEARCH
	#define DEDICATED_DECODER


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
#define QLEVELS 17
#ifdef N_CODER
#define EC_IDX_MASK ~0
#else
#define EC_IDX_MASK 0
#endif
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
static void update_CDFs(short *val, unsigned *stats, int *ctx, const unsigned *mixin_CDFs)
{
	unsigned *curr_CDF, sym;
	const unsigned *mcdf;
	for(int kc=0;kc<3;++kc)
	{
		curr_CDF=stats+ctx[kc];
		sym=val[kc]>>4;
		mcdf=mixin_CDFs+32*sym;
		__m256i c0=_mm256_loadu_si256((__m256i*)curr_CDF+0);
		__m256i c1=_mm256_loadu_si256((__m256i*)curr_CDF+1);
		__m256i c2=_mm256_loadu_si256((__m256i*)curr_CDF+2);
		__m256i c3=_mm256_loadu_si256((__m256i*)curr_CDF+3);
		__m256i m0=_mm256_load_si256((__m256i const*)mcdf+0);
		__m256i m1=_mm256_load_si256((__m256i const*)mcdf+1);
		__m256i m2=_mm256_load_si256((__m256i const*)mcdf+2);
		__m256i m3=_mm256_load_si256((__m256i const*)mcdf+3);
		__m256i d0=_mm256_sub_epi32(m0, c0);
		__m256i d1=_mm256_sub_epi32(m1, c1);
		__m256i d2=_mm256_sub_epi32(m2, c2);
		__m256i d3=_mm256_sub_epi32(m3, c3);
		d0=_mm256_srai_epi32(d0, 7);
		d1=_mm256_srai_epi32(d1, 7);
		d2=_mm256_srai_epi32(d2, 7);
		d3=_mm256_srai_epi32(d3, 7);
		c0=_mm256_add_epi32(c0, d0);
		c1=_mm256_add_epi32(c1, d1);
		c2=_mm256_add_epi32(c2, d2);
		c3=_mm256_add_epi32(c3, d3);
		_mm256_storeu_si256((__m256i*)curr_CDF+0, c0);
		_mm256_storeu_si256((__m256i*)curr_CDF+1, c1);
		_mm256_storeu_si256((__m256i*)curr_CDF+2, c2);
		_mm256_storeu_si256((__m256i*)curr_CDF+3, c3);
#ifdef ENABLE_GUIDE
		for(int k=1;k<33;++k)
		{
			if(curr_CDF[k]>0x10000||curr_CDF[k-1]>curr_CDF[k])
				LOG_ERROR("");
		}
#endif

		curr_CDF+=33+17*sym;
		sym=val[kc]&15;
		mcdf=mixin_CDFs+1024+16*sym;
		//mcdf=mixin_CDFs+32*sym;
		c0=_mm256_loadu_si256((__m256i*)curr_CDF+0);
		c1=_mm256_loadu_si256((__m256i*)curr_CDF+1);
		m0=_mm256_load_si256((__m256i*)mcdf+0);
		m1=_mm256_load_si256((__m256i*)mcdf+1);
		d0=_mm256_sub_epi32(m0, c0);
		d1=_mm256_sub_epi32(m1, c1);
		d0=_mm256_srai_epi32(d0, 8);
		d1=_mm256_srai_epi32(d1, 8);
		c0=_mm256_add_epi32(c0, d0);
		c1=_mm256_add_epi32(c1, d1);
		_mm256_storeu_si256((__m256i*)curr_CDF+0, c0);
		_mm256_storeu_si256((__m256i*)curr_CDF+1, c1);
	}
}
static void get_ctx(__m128i const *nb, short *pred, int *ctx)//nb={N, W, NW}
{
	__m128i hmasks[]=
	{
		_mm_set1_epi32(0x55555555),
		_mm_set1_epi32(0x33333333),
		_mm_set1_epi32(0x0F0F0F0F),
		_mm_set1_epi32(0x00FF00FF),
	};
	__m128i qhalf=_mm_set1_epi16(QLEVELS>>1);
	__m128i qmax=_mm_set1_epi16(QLEVELS-1);
	__m128i factor=_mm_set_epi16(QLEVELS, QLEVELS, QLEVELS, QLEVELS, 1, 1, 1, 1);
	__m128i shuf=_mm_set_epi8(
		 7,  6,  5,  4,  3,  2,  1,  0, 15, 14, 13, 12, 11, 10,  9,  8
	);
	//__m128i shuf2=_mm_set_epi8(
	//	-1, -1,  7,  6, -1, -1,  5,  4, -1, -1,  3,  2, -1, -1,  1,  0
	//);
	__m128i stride=_mm_set1_epi32(33+17*32);
	__m128i offset=_mm_set_epi32(
		(33+17*32)*QLEVELS*QLEVELS*3,
		(33+17*32)*QLEVELS*QLEVELS*2,
		(33+17*32)*QLEVELS*QLEVELS*1,
		(33+17*32)*QLEVELS*QLEVELS*0
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
	errors=_mm_abs_epi16(errors);//remove sign bit (9 -> 8-bit)
	errors=_mm_or_si128(errors, _mm_srli_epi16(errors, 1));//set LSBs
	errors=_mm_or_si128(errors, _mm_srli_epi16(errors, 2));
	errors=_mm_or_si128(errors, _mm_srli_epi16(errors, 4));
	errors=_mm_or_si128(errors, _mm_srli_epi16(errors, 8));
	errors=_mm_sub_epi16(errors, _mm_and_si128(_mm_srli_epi16(errors, 1), hmasks[0]));//hamming weight
	errors=_mm_add_epi16(_mm_and_si128(errors, hmasks[1]), _mm_and_si128(_mm_srli_epi16(errors, 2), hmasks[1]));
	errors=_mm_add_epi16(_mm_and_si128(errors, hmasks[2]), _mm_and_si128(_mm_srli_epi16(errors, 4), hmasks[2]));
	errors=_mm_add_epi16(_mm_and_si128(errors, hmasks[3]), _mm_and_si128(_mm_srli_epi16(errors, 8), hmasks[3]));//optimize: *0x0101

	//errors=_mm_srli_epi16(errors, 1);
	//errors=_mm_setzero_si128();//

	errors=_mm_xor_si128(errors, negmask);
	errors=_mm_sub_epi16(errors, negmask);
	errors=_mm_add_epi16(errors, qhalf);
	errors=_mm_min_epi16(errors, qmax);
	errors=_mm_max_epi16(errors, _mm_setzero_si128());
	errors=_mm_mullo_epi16(errors, factor);
	errors=_mm_add_epi16(errors, _mm_shuffle_epi8(errors, shuf));
	errors=_mm_cvtepi16_epi32(errors);
	errors=_mm_mullo_epi32(errors, stride);
	errors=_mm_add_epi32(errors, offset);
	_mm_store_si128((__m128i*)ctx, errors);
}
int f06_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
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
	dlist_init(list+0, 1, 1024, 0);
	dlist_init(list+1, 1, 1024, 0);
	dlist_init(list+2, 1, 1024, 0);
	ArithmeticCoder ec[3];
	if(image->depth!=8)
		LOG_ERROR("Unsupported depth %d", image->depth);
	int nlevels=1<<image->depth, clevels=nlevels<<1, half=nlevels>>1, chalf=nlevels;
	unsigned *mixin_CDFs=(unsigned*)_mm_malloc(sizeof(int[32*32+16*16]), sizeof(__m256i));
	unsigned *stats=(unsigned*)malloc(sizeof(int[(33+17*32)*QLEVELS*QLEVELS*4]));//(CDFSIZE+1) * nodes_in_tree * 4 channels max
	short *pixels=(short*)malloc((image->iw+4LL)*sizeof(short[2*4*2]));//2 padded rows * 4 channels max * {pixels, errors}
	if(!mixin_CDFs||!stats||!pixels)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	//initialize mixin_CDFs
	for(int kt=0;kt<32;++kt)
	{
		unsigned *curr_CDF=mixin_CDFs+32*kt;
		for(int ks=0;ks<32;++ks)
			curr_CDF[ks]=(ks>kt)*(0x10000-32)+ks;
		//{
		//	int amplitude;//X
		//	if(ks<kt)
		//		amplitude=0;
		//	else if(ks==kt)
		//		amplitude=(0x10000-32)/3;
		//	else if(ks==kt+1)
		//		amplitude=(0x10000-32)*2/3;
		//	else
		//		amplitude=(0x10000-32);
		//	curr_CDF[ks]=amplitude+ks;
		//}
	}
	for(int kt=0;kt<16;++kt)
	{
		unsigned *curr_CDF=mixin_CDFs+1024+16*kt;
		for(int ks=0;ks<16;++ks)
			curr_CDF[ks]=(ks>kt)*(0x10000-16)+ks;
	}
	//initialize to bypass
	for(int ks=0;ks<33;++ks)
		stats[ks]=(ks<<16)/32;
	for(int k=0;k<32;++k)
	{
		unsigned *curr_CDF=stats+33+17*k;
		for(int ks=0;ks<17;++ks)
			curr_CDF[ks]=(ks<<16)/16;
	}
	memfill(stats+33+17*32, stats, sizeof(int[(33+17*32)*QLEVELS*QLEVELS*4])-sizeof(int[33+17*32]), sizeof(int[33+17*32]));
	memset(pixels, 0, (image->iw+4LL)*sizeof(short[2*4*2]));
	//unsigned *CDF[4]=
	//{
	//	stats+(33+17*32)*9*9*0,
	//	stats+(33+17*32)*9*9*1,
	//	stats+(33+17*32)*9*9*2,
	//	stats+(33+17*32)*9*9*3,
	//};
	PROF(INIT);
	//ALIGN(32) const int ramp[]=
	//{
	//	 0,  1,  2,  3,  4,  5,  6,  7,
	//	 8,  9, 10, 11, 12, 13, 14, 15,
	//	16, 17, 18, 19, 20, 21, 22, 23,
	//	24, 25, 26, 27, 28, 29, 30, 31,
	//};
	if(fwd)
	{
		if(image->nch==3)
		{
			ac_enc_init(ec+0, list+0);
			ac_enc_init(ec+1, list+1);
			ac_enc_init(ec+2, list+2);
			for(int ky=0, idx=0;ky<image->ih;++ky)
			{
				short *rows[2]=
				{
					pixels+(((image->iw+4LL)*(ky&1)+1)<<3),
					pixels+(((image->iw+4LL)*((ky-1)&1)+1)<<3),
				};
				for(int kx=0;kx<image->iw;++kx, idx+=3)
				{
					//if(ky==256&&kx==256)
					//	printf("");
					__m128i nb[]=
					{
						_mm_loadu_si128((const __m128i*)(rows[1]+0)),//N
						_mm_loadu_si128((const __m128i*)(rows[0]-8)),//W
						_mm_loadu_si128((const __m128i*)(rows[1]-8)),//NW
					};
					ALIGN(16) short pred[8];
					ALIGN(16) int ctx[4];
					get_ctx(nb, pred, ctx);
#if 0
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

					//get ctx
					__m128i errors=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(mN), _mm_castsi128_ps(mW), _MM_SHUFFLE(3, 2, 3, 2)));//hi halves contain errors
					__m128i negmask=_mm_cmplt_epi16(errors, _mm_setzero_si128());
					errors=_mm_abs_epi16(errors);//remove sign bit (9 -> 8-bit)
					errors=_mm_or_si128(errors, _mm_srli_epi16(errors, 1));//set LSBs
					errors=_mm_or_si128(errors, _mm_srli_epi16(errors, 2));
					errors=_mm_or_si128(errors, _mm_srli_epi16(errors, 4));
					errors=_mm_or_si128(errors, _mm_srli_epi16(errors, 8));
					errors=_mm_sub_epi16(errors, _mm_and_si128(_mm_srli_epi16(errors, 1), hmasks[0]));//hamming weight
					errors=_mm_sub_epi16(_mm_and_si128(errors, hmasks[1]), _mm_and_si128(_mm_srli_epi16(errors, 2), hmasks[1]));
					errors=_mm_sub_epi16(_mm_and_si128(errors, hmasks[2]), _mm_and_si128(_mm_srli_epi16(errors, 4), hmasks[2]));
					errors=_mm_sub_epi16(_mm_and_si128(errors, hmasks[3]), _mm_and_si128(_mm_srli_epi16(errors, 8), hmasks[3]));//optimize: *0x0101
					errors=_mm_add_epi16(errors, qhalf);
					errors=_mm_min_epi16(errors, qmax);
					errors=_mm_max_epi16(errors, _mm_setzero_si128());
					errors=_mm_mullo_epi16(errors, factor);
					errors=_mm_add_epi16(errors, _mm_shuffle_epi8(errors, shuf));
					errors=_mm_cvtepi16_epi32(errors);
					errors=_mm_mullo_epi32(errors, stride);
					errors=_mm_add_epi32(errors, offset);
					ALIGN(16) int ctx[4];
					_mm_store_si128((__m128i*)ctx, errors);
#endif
					short *curr=rows[0];
					rows[0]+=8;
					rows[1]+=8;
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
					curr[4]=(short)(val[0]-chalf);
					curr[5]=(short)(val[1]- half);
					curr[6]=(short)(val[2]-chalf);
					PROF(RCT);
					//if(!ky&&kx==42)//
					//if(idx==84945-3)//
					//if(idx==2319)//
					//if(idx==9)//
					//if(idx==9621)//
					//if(idx==7299)//
					//	printf("");
					ac_enc(ec+(0&EC_IDX_MASK), val[1]>>4, stats+ctx[1]);
					ac_enc(ec+(1&EC_IDX_MASK), val[2]>>4, stats+ctx[2]);
					ac_enc(ec+(2&EC_IDX_MASK), val[0]>>4, stats+ctx[0]);
					ac_enc(ec+(0&EC_IDX_MASK), val[1]&15, stats+ctx[1]+33+17*(val[1]>>4));
					ac_enc(ec+(1&EC_IDX_MASK), val[2]&15, stats+ctx[2]+33+17*(val[2]>>4));
					ac_enc(ec+(2&EC_IDX_MASK), val[0]&15, stats+ctx[0]+33+17*(val[0]>>4));
					PROF(EC);

					update_CDFs(val, stats, ctx, mixin_CDFs);
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
#ifdef N_CODER
			ac_dec_init(ec+1, ptr, ptr+bm[1]);	ptr+=bm[1];
			ac_dec_init(ec+2, ptr, ptr+bm[2]);	ptr+=bm[2];
#endif
			for(int ky=0, idx=0;ky<image->ih;++ky)
			{
				short *rows[2]=
				{
					pixels+(((image->iw+4LL)*(ky&1)+1)<<3),
					pixels+(((image->iw+4LL)*((ky-1)&1)+1)<<3),
				};
				for(int kx=0;kx<image->iw;++kx, idx+=3)
				{
					__m128i nb[]=
					{
						_mm_loadu_si128((const __m128i*)(rows[1]+0)),//N
						_mm_loadu_si128((const __m128i*)(rows[0]-8)),//W
						_mm_loadu_si128((const __m128i*)(rows[1]-8)),//NW
					};
					ALIGN(16) short pred[8];
					ALIGN(16) int ctx[4];
					get_ctx(nb, pred, ctx);
					//__m128i mNW	=_mm_loadu_si128((const __m128i*)(rows[1]-4));
					//__m128i mN	=_mm_loadu_si128((const __m128i*)(rows[1]+0));
					//__m128i mW	=_mm_loadu_si128((const __m128i*)(rows[0]-4));
					//__m128i vmin=_mm_min_epi16(mN, mW);
					//__m128i vmax=_mm_max_epi16(mN, mW);
					//__m128i mpred=_mm_add_epi16(mN, mW);
					//mpred=_mm_sub_epi16(mpred, mNW);
					//mpred=_mm_min_epi16(mpred, vmax);
					//mpred=_mm_max_epi16(mpred, vmin);
					//ALIGN(16) short pred[8];
					//_mm_store_si128((__m128i*)pred, mpred);
					
					short *curr=rows[0];

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
					
					short val[3];
#ifdef DEDICATED_DECODER
					val[1]=(short)ac_dec_5bit(ec+(0&EC_IDX_MASK), stats+ctx[1]);
					val[2]=(short)ac_dec_5bit(ec+(1&EC_IDX_MASK), stats+ctx[2]);
					val[0]=(short)ac_dec_5bit(ec+(2&EC_IDX_MASK), stats+ctx[0]);
					val[1]=(short)(val[1]<<4|ac_dec_4bit(ec+(0&EC_IDX_MASK), stats+ctx[1]+33+17*val[1]));
					val[2]=(short)(val[2]<<4|ac_dec_4bit(ec+(1&EC_IDX_MASK), stats+ctx[2]+33+17*val[2]));
					val[0]=(short)(val[0]<<4|ac_dec_4bit(ec+(2&EC_IDX_MASK), stats+ctx[0]+33+17*val[0]));
#else
					val[1]=ac_dec(ec+(0&EC_IDX_MASK), stats+ctx[1], 32);
					val[2]=ac_dec(ec+(1&EC_IDX_MASK), stats+ctx[2], 32);
					val[0]=ac_dec(ec+(2&EC_IDX_MASK), stats+ctx[0], 32);
					val[1]=val[1]<<4|ac_dec(ec+(0&EC_IDX_MASK), stats+ctx[1]+33+17*val[1], 16);
					val[2]=val[2]<<4|ac_dec(ec+(1&EC_IDX_MASK), stats+ctx[2]+33+17*val[2], 16);
					val[0]=val[0]<<4|ac_dec(ec+(2&EC_IDX_MASK), stats+ctx[0]+33+17*val[0], 16);
#endif
					PROF(EC);

					curr[4]=(short)(val[0]-chalf);
					curr[5]=(short)(val[1]- half);
					curr[6]=(short)(val[2]-chalf);

					curr[0]=(short)(((val[0]+pred[0])&(clevels-1))-chalf);
					curr[1]=(short)(((val[1]+pred[1])&(nlevels-1))- half);
					curr[2]=(short)(((val[2]+pred[2])&(clevels-1))-chalf);
					
					short *rgb=dst->data+idx;
					memcpy(rgb, curr, sizeof(short[3]));
					rgb[1]-=(rgb[0]+rgb[2])>>2;
					rgb[2]+=rgb[1];
					rgb[0]+=rgb[1];
					PROF(RCT);
					
					update_CDFs(val, stats, ctx, mixin_CDFs);
					PROF(CDF);
#ifdef ENABLE_GUIDE
					curr=rgb;
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
		printf("F06  %c %15.6lf sec\n", 'D'+fwd, t0);
		prof_print();
	}
	dlist_clear(list+0);
	dlist_clear(list+1);
	dlist_clear(list+2);
	_mm_free(mixin_CDFs);
	free(stats);
	free(pixels);
	return 1;
}