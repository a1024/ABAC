#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define PROFILER 1
	#define N_CODER
	#define DEDICATED_AC
//	#define SHOW_PRED_ERRORS


#ifdef PROFILER
#define CHECKPOINTLIST\
	CHECKPOINT(INIT)\
	CHECKPOINT(RCT)\
	CHECKPOINT(PRED)\
	CHECKPOINT(DUMMY)\
	CHECKPOINT(EC)\
	CHECKPOINT(CDF)\
	CHECKPOINT(WP_ADD)\
	CHECKPOINT(WP_UPDATE)\
	CHECKPOINT(FINISH)
#endif
#include"ac.h"
#include"profiler.h"
#define QLEVELS 17
#define WP_UPDATE_MASK 0	//7
#define NWP 8
#define LGBLOCKSIZE 4
#define BLOCKSIZE (1<<LGBLOCKSIZE)
#define WP_BITS 16
#ifdef N_CODER
#define EC_IDX_MASK ~0
#else
#define EC_IDX_MASK 0
#endif
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
typedef enum NeighborIndexEnum
{
	NB_N,
	NB_W,
	NB_NW,
	NB_NE,
	NB_NEE,
	NB_NEEE,
	NB_NN,
	NB_WW,
	NB_NNWW,
	NB_NNE,
} NeighborIndex;
static void get_ctx(__m128i const *nb, const int *weights, int *pred, int *ctx)//nb={N, W, NW, NE, NEE, NEEE, NN, WW, NNWW, NNE}
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
	__m128i stride=_mm_set1_epi32(33+17*32);
	__m128i offset=_mm_set_epi32(
		(33+17*32)*QLEVELS*QLEVELS*3,
		(33+17*32)*QLEVELS*QLEVELS*2,
		(33+17*32)*QLEVELS*QLEVELS*1,
		(33+17*32)*QLEVELS*QLEVELS*0
	);
	//__m128i three=_mm_set1_epi16(3);
	
	__m128i eN=_mm_shuffle_epi32(nb[NB_N], _MM_SHUFFLE(1, 0, 3, 2));
	__m128i eW=_mm_shuffle_epi32(nb[NB_W], _MM_SHUFFLE(1, 0, 3, 2));
	__m128i wt2=_mm_slli_epi16(nb[NB_W], 2);
	//__m128i wt2=_mm_slli_epi16(nb[NB_W], 1);
	__m128i weighted=_mm_add_epi16(nb[NB_N], nb[NB_NE]);
	__m128i weird=_mm_add_epi16(nb[NB_W], nb[NB_NE]);
	__m128i cgrad=_mm_add_epi16(nb[NB_N], nb[NB_W]);
	__m128i vmin=_mm_min_epi16(nb[NB_N], nb[NB_W]);
	__m128i vmax=_mm_max_epi16(nb[NB_N], nb[NB_W]);
	__m128i mp2=_mm_srai_epi16(cgrad, 1);
	__m128i funky=_mm_add_epi16(nb[NB_W], nb[NB_NEE]);
	__m128i extrapN=_mm_slli_epi16(nb[NB_N], 1);
	__m128i top=_mm_add_epi16(nb[NB_N], nb[NB_NN]);
	wt2=_mm_add_epi16(wt2, nb[NB_N]);
	//wt2=_mm_add_epi16(wt2, nb[NB_NEEE]);
	weighted=_mm_srai_epi16(weighted, 1);
	//weighted=_mm_add_epi16(weighted, nb[NB_NNE]);
	weird=_mm_sub_epi16(weird, nb[NB_N]);
	cgrad=_mm_sub_epi16(cgrad, nb[NB_NW]);
	cgrad=_mm_sub_epi16(cgrad, _mm_srai_epi16(_mm_add_epi16(eN, eW), 5));
	wt2=_mm_add_epi16(wt2, nb[NB_NE]);
	//wt2=_mm_sub_epi16(wt2, nb[NB_N]);
	cgrad=_mm_min_epi16(cgrad, vmax);
	cgrad=_mm_max_epi16(cgrad, vmin);
	wt2=_mm_add_epi16(wt2, nb[NB_NEE]);
	funky=_mm_srai_epi16(funky, 1);
	extrapN=_mm_sub_epi16(extrapN, nb[NB_NN]);
	top=_mm_srai_epi16(top, 1);
	wt2=_mm_add_epi16(wt2, nb[NB_NEEE]);
	wt2=_mm_srai_epi16(wt2, 3);
	//wt2=_mm_srai_epi16(wt2, 1);

	vmin=_mm_min_epi16(vmin, nb[NB_NE]);
	vmax=_mm_max_epi16(vmax, nb[NB_NE]);
	cgrad	=_mm_cvtepi16_epi32(cgrad);
	mp2	=_mm_cvtepi16_epi32(mp2);
	funky	=_mm_cvtepi16_epi32(funky);
	weird	=_mm_cvtepi16_epi32(weird);
	extrapN	=_mm_cvtepi16_epi32(extrapN);
	wt2	=_mm_cvtepi16_epi32(wt2);
	top	=_mm_cvtepi16_epi32(top);
	weighted=_mm_cvtepi16_epi32(weighted);
	vmin=_mm_cvtepi16_epi32(vmin);
	vmax=_mm_cvtepi16_epi32(vmax);
	__m128i wp1=_mm_mullo_epi32(cgrad	, _mm_load_si128((__m128i const*)weights+0));//median3(N, W, N+W-NW)
	__m128i wp2=_mm_mullo_epi32(mp2		, _mm_load_si128((__m128i const*)weights+1));//(N+W)>>1
	__m128i wp3=_mm_mullo_epi32(funky	, _mm_load_si128((__m128i const*)weights+2));//(W+NEE)>>1
	__m128i wp4=_mm_mullo_epi32(weird	, _mm_load_si128((__m128i const*)weights+3));//W+NE-N
	__m128i wp5=_mm_mullo_epi32(extrapN	, _mm_load_si128((__m128i const*)weights+4));//2*N-NN
	__m128i wp6=_mm_mullo_epi32(wt2		, _mm_load_si128((__m128i const*)weights+5));//(4*W+N+NE+NEE+NEEE)>>3		//(2*W+NEEE-N)>>1
	__m128i wp7=_mm_mullo_epi32(top		, _mm_load_si128((__m128i const*)weights+6));//(N+NN)>>1
	__m128i wp8=_mm_mullo_epi32(weighted	, _mm_load_si128((__m128i const*)weights+7));//(N+NE)>>1			//N+NE-NNE
	wp1=_mm_add_epi32(wp1, wp2);
	wp1=_mm_add_epi32(wp1, wp3);
	wp1=_mm_add_epi32(wp1, wp4);
	wp1=_mm_add_epi32(wp1, wp5);
	wp1=_mm_add_epi32(wp1, wp6);
	wp1=_mm_add_epi32(wp1, wp7);
	wp1=_mm_add_epi32(wp1, wp8);
	wp1=_mm_srai_epi32(wp1, WP_BITS);//weights must add up to 1<<WP_BITS
	wp1=_mm_min_epi32(wp1, vmax);
	wp1=_mm_max_epi32(wp1, vmin);
	_mm_store_si128((__m128i*)pred+0, wp1);
	_mm_store_si128((__m128i*)pred+1, cgrad);
	_mm_store_si128((__m128i*)pred+2, mp2);
	_mm_store_si128((__m128i*)pred+3, funky);
	_mm_store_si128((__m128i*)pred+4, weird);
	_mm_store_si128((__m128i*)pred+5, extrapN);
	_mm_store_si128((__m128i*)pred+6, wt2);
	_mm_store_si128((__m128i*)pred+7, top);
	_mm_store_si128((__m128i*)pred+8, weighted);
#if 0
	__m128i wp	=_mm_mullo_epi16(cgrad		, _mm_load_si128((__m128i const*)weights+0));//median3(N, W, N+W-NW)
	__m128i wp2	=_mm_mullo_epi16(mp2		, _mm_load_si128((__m128i const*)weights+1));//(N+W)>>1
	__m128i wp3	=_mm_mullo_epi16(funky		, _mm_load_si128((__m128i const*)weights+2));//(W+NEE)>>1
	__m128i wp4	=_mm_mullo_epi16(weird		, _mm_load_si128((__m128i const*)weights+3));//W+NE-N
	__m128i wp5	=_mm_mullo_epi16(extrapN	, _mm_load_si128((__m128i const*)weights+4));//2*N-NN
	__m128i wp6	=_mm_mullo_epi16(extrapW	, _mm_load_si128((__m128i const*)weights+5));//2*W-WW
	__m128i wp7	=_mm_mullo_epi16(top		, _mm_load_si128((__m128i const*)weights+6));//(N+NN)>>1
	__m128i wp8	=_mm_mullo_epi16(weighted	, _mm_load_si128((__m128i const*)weights+7));//(3*W+NEEE)>>1
	wp=_mm_add_epi16(wp, wp2);
	wp=_mm_add_epi16(wp, wp3);
	wp=_mm_add_epi16(wp, wp4);
	wp=_mm_add_epi16(wp, wp5);
	wp=_mm_add_epi16(wp, wp6);
	wp=_mm_add_epi16(wp, wp7);
	wp=_mm_add_epi16(wp, wp8);
	wp=_mm_srai_epi16(wp, WP_BITS);//weights must add up to 1<<WP_BITS
	wp=_mm_min_epi16(wp, vmax);
	wp=_mm_max_epi16(wp, vmin);
	_mm_store_si128((__m128i*)pred+0, wp);
	_mm_store_si128((__m128i*)pred+1, cgrad);
	_mm_store_si128((__m128i*)pred+2, mp2);
	_mm_store_si128((__m128i*)pred+3, funky);
	_mm_store_si128((__m128i*)pred+4, weird);
	_mm_store_si128((__m128i*)pred+5, extrapN);
	_mm_store_si128((__m128i*)pred+6, extrapW);
	_mm_store_si128((__m128i*)pred+7, top);
	_mm_store_si128((__m128i*)pred+8, weighted);
#endif

	//get ctx
	__m128i errors=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(nb[NB_N]), _mm_castsi128_ps(nb[NB_W]), _MM_SHUFFLE(3, 2, 3, 2)));//hi halves contain errors
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
static void update_wp(int *weights, int *errors)
{
	for(int kc=0;kc<3;++kc)
	{
		long long
			w0=0x400000LL/(errors[kc+4*0]+1LL),
			w1=0x400000LL/(errors[kc+4*1]+1LL),
			w2=0x400000LL/(errors[kc+4*2]+1LL),
			w3=0x400000LL/(errors[kc+4*3]+1LL),
			w4=0x400000LL/(errors[kc+4*4]+1LL),
			w5=0x400000LL/(errors[kc+4*5]+1LL),
			w6=0x400000LL/(errors[kc+4*6]+1LL),
			w7=0x400000LL/(errors[kc+4*7]+1LL),
			sum=w0+w1+w2+w3+w4+w5+w6+w7;
		sum+=!sum;
		w0=((w0<<WP_BITS)+(sum>>1))/sum;
		w1=((w1<<WP_BITS)+(sum>>1))/sum;
		w2=((w2<<WP_BITS)+(sum>>1))/sum;
		w3=((w3<<WP_BITS)+(sum>>1))/sum;
		w4=((w4<<WP_BITS)+(sum>>1))/sum;
		w5=((w5<<WP_BITS)+(sum>>1))/sum;
		w6=((w6<<WP_BITS)+(sum>>1))/sum;
		w7=((w7<<WP_BITS)+(sum>>1))/sum;
#if 0
			e0=errors[kc+ 0],
			e1=errors[kc+ 4],
			e2=errors[kc+ 8],
			e3=errors[kc+12],
			sum=e0+e1+e2+e3+1,
			w3=(((e0+e1+e2)<<WP_BITS)+(sum>>1))/sum,
			w2=(((e0+e1+e3)<<WP_BITS)+(sum>>1))/sum,
			w1=(((e0+e2+e3)<<WP_BITS)+(sum>>1))/sum,
			w0=(1<<WP_BITS)-(w1+w2+w3);
#endif
		weights[kc+4*0]=(int)w0;
		weights[kc+4*1]=(int)w1;
		weights[kc+4*2]=(int)w2;
		weights[kc+4*3]=(int)w3;
		weights[kc+4*4]=(int)w4;
		weights[kc+4*5]=(int)w5;
		weights[kc+4*6]=(int)w6;
		weights[kc+4*7]=(int)w7;
		//errors[kc+4*0]>>=3;
		//errors[kc+4*1]>>=3;
		//errors[kc+4*2]>>=3;
		//errors[kc+4*3]>>=3;
		//errors[kc+4*4]>>=3;
		//errors[kc+4*5]>>=3;
		//errors[kc+4*6]>>=3;
		//errors[kc+4*7]>>=3;
	}
	memset(errors, 0, sizeof(int[4*NWP]));
}
int f08_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
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
	int nlevels=1<<image->depth, clevels=nlevels<<1, half=nlevels>>1, chalf=nlevels;
	__m128i mhalf=_mm_set_epi16(0, 0, 0, 0, 0, chalf, half, chalf);
	__m128i mhalf32=_mm_set_epi32(0, chalf, half, chalf);
	__m128i pxmask=_mm_set_epi16(0, 0, 0, 0, 0, clevels-1, nlevels-1, clevels-1);
	__m128i pxmask32=_mm_set_epi32(0, clevels-1, nlevels-1, clevels-1);
	__m128i pack16=_mm_set_epi8(
		-1, -1, -1, -1, -1, -1, -1, -1, 13, 12,  9,  8,  5,  4,  1,  0
	);
	unsigned *mixin_CDFs=(unsigned*)_mm_malloc(sizeof(int[32*32+16*16]), sizeof(__m256i));
	unsigned *stats=(unsigned*)malloc(sizeof(int[(33+17*32)*QLEVELS*QLEVELS*4]));//(CDFSIZE+1) * nodes_in_tree * 4 channels max
	short *pixels=(short*)malloc((image->iw+16LL)*sizeof(short[4*4*2]));//4 padded rows * 4 channels max * {pixels, errors}
	int nblocks=(image->iw+BLOCKSIZE-1)>>LGBLOCKSIZE;
	int *weights=(int*)_mm_malloc(nblocks*sizeof(__m128i[NWP]), sizeof(__m128i));
	if(!mixin_CDFs||!stats||!pixels||!weights)
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
	memset(pixels, 0, (image->iw+16LL)*sizeof(short[4*4*2]));
	weights[0]=(1<<WP_BITS)/NWP;
	memfill(weights+1, weights, nblocks*sizeof(__m128i[NWP])-sizeof(*weights), sizeof(*weights));
	PROF(INIT);
	ALIGN(32) int errors[4*NWP]={0};
	ALIGN(32) int pred[4*(1+NWP)];
	ALIGN(32) int ctx[4];
	ALIGN(16) short val[8]={0};
#ifdef SHOW_PRED_ERRORS
	long long total_errors[_countof(errors)]={0};
#endif
	//__m128i weights[]=
	//{
	//	_mm_set1_epi16((1<<WP_BITS)/NWP),
	//	_mm_set1_epi16((1<<WP_BITS)/NWP),
	//	_mm_set1_epi16((1<<WP_BITS)/NWP),
	//	_mm_set1_epi16((1<<WP_BITS)/NWP),
	//	_mm_set1_epi16((1<<WP_BITS)/NWP),
	//	_mm_set1_epi16((1<<WP_BITS)/NWP),
	//};
	if(fwd)
	{
		ac_enc_init(ec+0, list+0);
		ac_enc_init(ec+1, list+1);
		ac_enc_init(ec+2, list+2);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			int updatewp=!(ky&WP_UPDATE_MASK);
			int *wp=weights;
			short *rows[]=
			{
				pixels+(((image->iw+16LL)*((ky-0)&1)+2)<<3),
				pixels+(((image->iw+16LL)*((ky-1)&1)+2)<<3),
				pixels+(((image->iw+16LL)*((ky-2)&1)+2)<<3),
				pixels+(((image->iw+16LL)*((ky-3)&1)+2)<<3),
			};
			for(int kx=0;kx<image->iw;++kx, idx+=3)
			{
				//if(ky==256&&kx==256)
				//if(idx==3462)//
				//if(ky==1&&kx==386)//
				//if(idx==29889)//
				//if(idx==2115)//
				//if(idx==48321)//
				//	printf("");
				__m128i nb[]=
				{
					_mm_loadu_si128((const __m128i*)(rows[1]+0)),//N
					_mm_loadu_si128((const __m128i*)(rows[0]-8)),//W
					_mm_loadu_si128((const __m128i*)(rows[1]-8)),//NW
					_mm_loadu_si128((const __m128i*)(rows[1]+8)),//NE
					_mm_loadu_si128((const __m128i*)(rows[1]+16)),//NEE
					_mm_loadu_si128((const __m128i*)(rows[1]+24)),//NEEE
					_mm_loadu_si128((const __m128i*)(rows[2]+0)),//NN
					_mm_loadu_si128((const __m128i*)(rows[0]-16)),//WW
					_mm_loadu_si128((const __m128i*)(rows[2]-16)),//NNWW
					_mm_loadu_si128((const __m128i*)(rows[2]+8)),//NNE
				};
				get_ctx(nb, wp, pred, ctx);

				short *curr=rows[0];
				PROF(PRED);
				PROF(DUMMY);

				memcpy(curr, image->data+idx, sizeof(short[3]));
				curr[0]-=curr[1];
				curr[2]-=curr[1];
				curr[1]+=(curr[0]+curr[2])>>2;

				__m128i mc=_mm_loadu_si128((__m128i*)curr);
				mc=_mm_cvtepi16_epi32(mc);
				__m128i wp0=_mm_load_si128((__m128i*)pred+0);
				wp0=_mm_sub_epi32(mc, wp0);
				wp0=_mm_shuffle_epi8(wp0, pack16);
				wp0=_mm_add_epi16(wp0, mhalf);
				wp0=_mm_and_si128(wp0, pxmask);
				_mm_store_si128((__m128i*)val, wp0);
				wp0=_mm_sub_epi16(wp0, mhalf);
				_mm_storeu_si128((__m128i*)(curr+4), wp0);

				//short val[]=
				//{
				//	(curr[0]-pred[0]+chalf)&(clevels-1),
				//	(curr[1]-pred[1]+ half)&(nlevels-1),
				//	(curr[2]-pred[2]+chalf)&(clevels-1),
				//};
				//curr[4]=val[0]-chalf;
				//curr[5]=val[1]- half;
				//curr[6]=val[2]-chalf;
				PROF(RCT);
				//if(!ky&&kx==42)//
				//if(idx==84945-3)//
				//if(idx==2319)//
				//if(idx==9)//
				//if(idx==9621)//
				//if(idx==7299)//
				//if(ky==1&&kx==1)//
				//if(ky==1&&kx==386)//
				//if(idx==3462)//
				//	printf("");

#ifdef DEDICATED_AC
				int sym[]=
				{
					val[0]>>4,
					val[1]>>4,
					val[2]>>4,
					val[0]&15,
					val[1]&15,
					val[2]&15,
				};
				unsigned *CDFs[]=
				{
					stats+ctx[0],
					stats+ctx[1],
					stats+ctx[2],
					stats+ctx[0]+33+17*sym[0],
					stats+ctx[1]+33+17*sym[1],
					stats+ctx[2]+33+17*sym[2],
				};
				unsigned cdf[]=
				{
					CDFs[0][sym[0]],
					CDFs[1][sym[1]],
					CDFs[2][sym[2]],
					CDFs[3][sym[3]],
					CDFs[4][sym[4]],
					CDFs[5][sym[5]],
				};
				int freq[]=
				{
					CDFs[0][sym[0]+1]-cdf[0],
					CDFs[1][sym[1]+1]-cdf[1],
					CDFs[2][sym[2]+1]-cdf[2],
					CDFs[3][sym[3]+1]-cdf[3],
					CDFs[4][sym[4]+1]-cdf[4],
					CDFs[5][sym[5]+1]-cdf[5],
				};
				ec[0].low+=ec[0].range*cdf[0]>>16;
				ec[1].low+=ec[1].range*cdf[1]>>16;
				ec[2].low+=ec[2].range*cdf[2]>>16;
				ec[0].range=(ec[0].range*freq[0]>>16)-1;//must decrement hi because decoder fails when code == hi2
				ec[1].range=(ec[1].range*freq[1]>>16)-1;
				ec[2].range=(ec[2].range*freq[2]>>16)-1;
				while(ec[0].range<(1LL<<PROB_BITS))//only when freq=1 -> range=0, this loop runs twice
					ac_enc_renorm(ec+0);
				while(ec[1].range<(1LL<<PROB_BITS))
					ac_enc_renorm(ec+1);
				while(ec[2].range<(1LL<<PROB_BITS))
					ac_enc_renorm(ec+2);

				ec[0].low+=ec[0].range*cdf[3]>>16;
				ec[1].low+=ec[1].range*cdf[4]>>16;
				ec[2].low+=ec[2].range*cdf[5]>>16;
				ec[0].range=(ec[0].range*freq[3]>>16)-1;
				ec[1].range=(ec[1].range*freq[4]>>16)-1;
				ec[2].range=(ec[2].range*freq[5]>>16)-1;
				while(ec[0].range<(1LL<<PROB_BITS))
					ac_enc_renorm(ec+0);
				while(ec[1].range<(1LL<<PROB_BITS))
					ac_enc_renorm(ec+1);
				while(ec[2].range<(1LL<<PROB_BITS))
					ac_enc_renorm(ec+2);
#else
				ac_enc(ec+(0&EC_IDX_MASK), val[0]>>4, stats+ctx[0]);
				ac_enc(ec+(1&EC_IDX_MASK), val[1]>>4, stats+ctx[1]);
				ac_enc(ec+(2&EC_IDX_MASK), val[2]>>4, stats+ctx[2]);
				ac_enc(ec+(0&EC_IDX_MASK), val[0]&15, stats+ctx[0]+33+17*(val[0]>>4));
				ac_enc(ec+(1&EC_IDX_MASK), val[1]&15, stats+ctx[1]+33+17*(val[1]>>4));
				ac_enc(ec+(2&EC_IDX_MASK), val[2]&15, stats+ctx[2]+33+17*(val[2]>>4));
#endif
				PROF(EC);

				update_CDFs(val, stats, ctx, mixin_CDFs);
				PROF(CDF);
				if(updatewp)
				{
				__m128i p0=_mm_load_si128((__m128i*)pred+1);//calculate pred errors
				__m128i p1=_mm_load_si128((__m128i*)pred+2);
				__m128i p2=_mm_load_si128((__m128i*)pred+3);
				__m128i p3=_mm_load_si128((__m128i*)pred+4);
				__m128i p4=_mm_load_si128((__m128i*)pred+5);
				__m128i p5=_mm_load_si128((__m128i*)pred+6);
				__m128i p6=_mm_load_si128((__m128i*)pred+7);
				__m128i p7=_mm_load_si128((__m128i*)pred+8);
				__m128i e0=_mm_load_si128((__m128i*)errors+0);
				__m128i e1=_mm_load_si128((__m128i*)errors+1);
				__m128i e2=_mm_load_si128((__m128i*)errors+2);
				__m128i e3=_mm_load_si128((__m128i*)errors+3);
				__m128i e4=_mm_load_si128((__m128i*)errors+4);
				__m128i e5=_mm_load_si128((__m128i*)errors+5);
				__m128i e6=_mm_load_si128((__m128i*)errors+6);
				__m128i e7=_mm_load_si128((__m128i*)errors+7);
				p0=_mm_sub_epi32(mc, p0);
				p1=_mm_sub_epi32(mc, p1);
				p2=_mm_sub_epi32(mc, p2);
				p3=_mm_sub_epi32(mc, p3);
				p4=_mm_sub_epi32(mc, p4);
				p5=_mm_sub_epi32(mc, p5);
				p6=_mm_sub_epi32(mc, p6);
				p7=_mm_sub_epi32(mc, p7);
				p0=_mm_abs_epi32(p0);
				p1=_mm_abs_epi32(p1);
				p2=_mm_abs_epi32(p2);
				p3=_mm_abs_epi32(p3);
				p4=_mm_abs_epi32(p4);
				p5=_mm_abs_epi32(p5);
				p6=_mm_abs_epi32(p6);
				p7=_mm_abs_epi32(p7);
				e0=_mm_add_epi32(e0, p0);
				e1=_mm_add_epi32(e1, p1);
				e2=_mm_add_epi32(e2, p2);
				e3=_mm_add_epi32(e3, p3);
				e4=_mm_add_epi32(e4, p4);
				e5=_mm_add_epi32(e5, p5);
				e6=_mm_add_epi32(e6, p6);
				e7=_mm_add_epi32(e7, p7);
				_mm_store_si128((__m128i*)errors+0, e0);
				_mm_store_si128((__m128i*)errors+1, e1);
				_mm_store_si128((__m128i*)errors+2, e2);
				_mm_store_si128((__m128i*)errors+3, e3);
				_mm_store_si128((__m128i*)errors+4, e4);
				_mm_store_si128((__m128i*)errors+5, e5);
				_mm_store_si128((__m128i*)errors+6, e6);
				_mm_store_si128((__m128i*)errors+7, e7);
				PROF(WP_ADD);
				}
				rows[0]+=8;
				rows[1]+=8;
				rows[2]+=8;
				rows[3]+=8;
				if(updatewp&&kx&&!(kx&(BLOCKSIZE-1)))
				{
#ifdef SHOW_PRED_ERRORS
					for(int k=0;k<_countof(errors);++k)
						total_errors[k]+=errors[k];
#endif
					update_wp(wp, errors);
					wp+=sizeof(__m128i[NWP])/(sizeof(*wp));
					PROF(WP_UPDATE);
				}
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
			int updatewp=!(ky&WP_UPDATE_MASK);
			int *wp=weights;
			short *rows[]=
			{
				pixels+(((image->iw+16LL)*((ky-0)&1)+2)<<3),
				pixels+(((image->iw+16LL)*((ky-1)&1)+2)<<3),
				pixels+(((image->iw+16LL)*((ky-2)&1)+2)<<3),
				pixels+(((image->iw+16LL)*((ky-3)&1)+2)<<3),
			};
			for(int kx=0;kx<image->iw;++kx, idx+=3)
			{
				//if(idx==3462)//
				//if(ky==1&&kx==386)//
				//if(idx==29889)//
				//if(idx==2115)//
				//if(idx==48321)//
				//	printf("");
				__m128i nb[]=
				{
					_mm_loadu_si128((const __m128i*)(rows[1]+0)),//N
					_mm_loadu_si128((const __m128i*)(rows[0]-8)),//W
					_mm_loadu_si128((const __m128i*)(rows[1]-8)),//NW
					_mm_loadu_si128((const __m128i*)(rows[1]+8)),//NE
					_mm_loadu_si128((const __m128i*)(rows[1]+16)),//NEE
					_mm_loadu_si128((const __m128i*)(rows[1]+24)),//NEEE
					_mm_loadu_si128((const __m128i*)(rows[2]+0)),//NN
					_mm_loadu_si128((const __m128i*)(rows[0]-16)),//WW
					_mm_loadu_si128((const __m128i*)(rows[2]-16)),//NNWW
					_mm_loadu_si128((const __m128i*)(rows[2]+8)),//NNE
				};
				get_ctx(nb, wp, pred, ctx);
					
				short *curr=rows[0];
				PROF(PRED);
				PROF(DUMMY);

				//if(idx==1095)//
				//if(idx==2319)//
				//if(idx==9)//
				//if(idx==9621)//
				//if(idx==7299)//
				//if(ky==1&&kx==1)//
				//if(ky==1&&kx==386)//
				//if(idx==3462)//
				//	printf("");
#ifdef DEDICATED_AC
				unsigned *CDFs[]=
				{
					stats+ctx[0],
					stats+ctx[1],
					stats+ctx[2],
					//stats+ctx[0]+33+17*(val[0]>>4),
					//stats+ctx[1]+33+17*(val[1]>>4),
					//stats+ctx[2]+33+17*(val[2]>>4),
				};
				int sym[3];
				{
					unsigned long long range[]=
					{
						ec[0].code-ec[0].low,
						ec[1].code-ec[1].low,
						ec[2].code-ec[2].low,
					};
					sym[0] =(ec[0].range*CDFs[0][16]>>16<=range[0])<<4;
					sym[1] =(ec[1].range*CDFs[1][16]>>16<=range[1])<<4;
					sym[2] =(ec[2].range*CDFs[2][16]>>16<=range[2])<<4;
					sym[0]|=(ec[0].range*CDFs[0][sym[0]|8]>>16<=range[0])<<3;
					sym[1]|=(ec[1].range*CDFs[1][sym[1]|8]>>16<=range[1])<<3;
					sym[2]|=(ec[2].range*CDFs[2][sym[2]|8]>>16<=range[2])<<3;
					sym[0]|=(ec[0].range*CDFs[0][sym[0]|4]>>16<=range[0])<<2;
					sym[1]|=(ec[1].range*CDFs[1][sym[1]|4]>>16<=range[1])<<2;
					sym[2]|=(ec[2].range*CDFs[2][sym[2]|4]>>16<=range[2])<<2;
					sym[0]|=(ec[0].range*CDFs[0][sym[0]|2]>>16<=range[0])<<1;
					sym[1]|=(ec[1].range*CDFs[1][sym[1]|2]>>16<=range[1])<<1;
					sym[2]|=(ec[2].range*CDFs[2][sym[2]|2]>>16<=range[2])<<1;
					sym[0]|= ec[0].range*CDFs[0][sym[0]|1]>>16<=range[0];
					sym[1]|= ec[1].range*CDFs[1][sym[1]|1]>>16<=range[1];
					sym[2]|= ec[2].range*CDFs[2][sym[2]|1]>>16<=range[2];
					unsigned cdf[]=
					{
						CDFs[0][sym[0]],
						CDFs[1][sym[1]],
						CDFs[2][sym[2]],
					};
					int freq[]=
					{
						CDFs[0][sym[0]+1]-cdf[0],
						CDFs[1][sym[1]+1]-cdf[1],
						CDFs[2][sym[2]+1]-cdf[2],
					};
					ec[0].low+=ec[0].range*cdf[0]>>16;
					ec[1].low+=ec[1].range*cdf[1]>>16;
					ec[2].low+=ec[2].range*cdf[2]>>16;
					ec[0].range=(ec[0].range*freq[0]>>16)-1;
					ec[1].range=(ec[1].range*freq[1]>>16)-1;
					ec[2].range=(ec[2].range*freq[2]>>16)-1;
					while(ec[0].range<(1LL<<PROB_BITS))
						ac_dec_renorm(ec+0);
					while(ec[1].range<(1LL<<PROB_BITS))
						ac_dec_renorm(ec+1);
					while(ec[2].range<(1LL<<PROB_BITS))
						ac_dec_renorm(ec+2);
				}
				{
					unsigned long long range[]=
					{
						ec[0].code-ec[0].low,
						ec[1].code-ec[1].low,
						ec[2].code-ec[2].low,
					};
					CDFs[0]+=33+17*sym[0];
					CDFs[1]+=33+17*sym[1];
					CDFs[2]+=33+17*sym[2];
					val[0]=sym[0]<<4;
					val[1]=sym[1]<<4;
					val[2]=sym[2]<<4;
					sym[0] =(ec[0].range*CDFs[0][8]>>16<=range[0])<<3;
					sym[1] =(ec[1].range*CDFs[1][8]>>16<=range[1])<<3;
					sym[2] =(ec[2].range*CDFs[2][8]>>16<=range[2])<<3;
					sym[0]|=(ec[0].range*CDFs[0][sym[0]|4]>>16<=range[0])<<2;
					sym[1]|=(ec[1].range*CDFs[1][sym[1]|4]>>16<=range[1])<<2;
					sym[2]|=(ec[2].range*CDFs[2][sym[2]|4]>>16<=range[2])<<2;
					sym[0]|=(ec[0].range*CDFs[0][sym[0]|2]>>16<=range[0])<<1;
					sym[1]|=(ec[1].range*CDFs[1][sym[1]|2]>>16<=range[1])<<1;
					sym[2]|=(ec[2].range*CDFs[2][sym[2]|2]>>16<=range[2])<<1;
					sym[0]|= ec[0].range*CDFs[0][sym[0]|1]>>16<=range[0];
					sym[1]|= ec[1].range*CDFs[1][sym[1]|1]>>16<=range[1];
					sym[2]|= ec[2].range*CDFs[2][sym[2]|1]>>16<=range[2];
					unsigned cdf[]=
					{
						CDFs[0][sym[0]],
						CDFs[1][sym[1]],
						CDFs[2][sym[2]],
					};
					int freq[]=
					{
						CDFs[0][sym[0]+1]-cdf[0],
						CDFs[1][sym[1]+1]-cdf[1],
						CDFs[2][sym[2]+1]-cdf[2],
					};
					ec[0].low+=ec[0].range*cdf[0]>>16;
					ec[1].low+=ec[1].range*cdf[1]>>16;
					ec[2].low+=ec[2].range*cdf[2]>>16;
					ec[0].range=(ec[0].range*freq[0]>>16)-1;
					ec[1].range=(ec[1].range*freq[1]>>16)-1;
					ec[2].range=(ec[2].range*freq[2]>>16)-1;
					while(ec[0].range<(1LL<<PROB_BITS))
						ac_dec_renorm(ec+0);
					while(ec[1].range<(1LL<<PROB_BITS))
						ac_dec_renorm(ec+1);
					while(ec[2].range<(1LL<<PROB_BITS))
						ac_dec_renorm(ec+2);
				}
				val[0]|=sym[0];
				val[1]|=sym[1];
				val[2]|=sym[2];
#else
				val[0]=ac_dec_5bit(ec+(0&EC_IDX_MASK), stats+ctx[0]);
				val[1]=ac_dec_5bit(ec+(1&EC_IDX_MASK), stats+ctx[1]);
				val[2]=ac_dec_5bit(ec+(2&EC_IDX_MASK), stats+ctx[2]);
				val[0]=val[0]<<4|ac_dec_4bit(ec+(0&EC_IDX_MASK), stats+ctx[0]+33+17*val[0]);
				val[1]=val[1]<<4|ac_dec_4bit(ec+(1&EC_IDX_MASK), stats+ctx[1]+33+17*val[1]);
				val[2]=val[2]<<4|ac_dec_4bit(ec+(2&EC_IDX_MASK), stats+ctx[2]+33+17*val[2]);
#endif
				PROF(EC);

				__m128i mc=_mm_load_si128((__m128i*)val);
				__m128i wp0=_mm_load_si128((__m128i*)pred+0);
				__m128i me=_mm_sub_epi16(mc, mhalf);
				mc=_mm_cvtepi16_epi32(mc);
				mc=_mm_add_epi32(mc, wp0);
				mc=_mm_and_si128(mc, pxmask32);
				mc=_mm_sub_epi32(mc, mhalf32);
				wp0=_mm_shuffle_epi8(mc, pack16);
				_mm_storeu_si128((__m128i*)curr, wp0);
				_mm_storeu_si128((__m128i*)(curr+4), me);
				//curr[4]=val[0]-chalf;
				//curr[5]=val[1]- half;
				//curr[6]=val[2]-chalf;
				//
				//curr[0]=((val[0]+pred[0])&(clevels-1))-chalf;
				//curr[1]=((val[1]+pred[1])&(nlevels-1))- half;
				//curr[2]=((val[2]+pred[2])&(clevels-1))-chalf;
					
				short *rgb=dst->data+idx;
				memcpy(rgb, curr, sizeof(short[3]));
				rgb[1]-=(rgb[0]+rgb[2])>>2;
				rgb[2]+=rgb[1];
				rgb[0]+=rgb[1];
				PROF(RCT);
					
				update_CDFs(val, stats, ctx, mixin_CDFs);
				PROF(CDF);
				if(updatewp)
				{
				__m128i p0=_mm_load_si128((__m128i*)pred+1);//calculate pred errors
				__m128i p1=_mm_load_si128((__m128i*)pred+2);
				__m128i p2=_mm_load_si128((__m128i*)pred+3);
				__m128i p3=_mm_load_si128((__m128i*)pred+4);
				__m128i p4=_mm_load_si128((__m128i*)pred+5);
				__m128i p5=_mm_load_si128((__m128i*)pred+6);
				__m128i p6=_mm_load_si128((__m128i*)pred+7);
				__m128i p7=_mm_load_si128((__m128i*)pred+8);
				__m128i e0=_mm_load_si128((__m128i*)errors+0);
				__m128i e1=_mm_load_si128((__m128i*)errors+1);
				__m128i e2=_mm_load_si128((__m128i*)errors+2);
				__m128i e3=_mm_load_si128((__m128i*)errors+3);
				__m128i e4=_mm_load_si128((__m128i*)errors+4);
				__m128i e5=_mm_load_si128((__m128i*)errors+5);
				__m128i e6=_mm_load_si128((__m128i*)errors+6);
				__m128i e7=_mm_load_si128((__m128i*)errors+7);
				p0=_mm_sub_epi32(mc, p0);
				p1=_mm_sub_epi32(mc, p1);
				p2=_mm_sub_epi32(mc, p2);
				p3=_mm_sub_epi32(mc, p3);
				p4=_mm_sub_epi32(mc, p4);
				p5=_mm_sub_epi32(mc, p5);
				p6=_mm_sub_epi32(mc, p6);
				p7=_mm_sub_epi32(mc, p7);
				p0=_mm_abs_epi32(p0);
				p1=_mm_abs_epi32(p1);
				p2=_mm_abs_epi32(p2);
				p3=_mm_abs_epi32(p3);
				p4=_mm_abs_epi32(p4);
				p5=_mm_abs_epi32(p5);
				p6=_mm_abs_epi32(p6);
				p7=_mm_abs_epi32(p7);
				e0=_mm_add_epi32(e0, p0);
				e1=_mm_add_epi32(e1, p1);
				e2=_mm_add_epi32(e2, p2);
				e3=_mm_add_epi32(e3, p3);
				e4=_mm_add_epi32(e4, p4);
				e5=_mm_add_epi32(e5, p5);
				e6=_mm_add_epi32(e6, p6);
				e7=_mm_add_epi32(e7, p7);
				_mm_store_si128((__m128i*)errors+0, e0);
				_mm_store_si128((__m128i*)errors+1, e1);
				_mm_store_si128((__m128i*)errors+2, e2);
				_mm_store_si128((__m128i*)errors+3, e3);
				_mm_store_si128((__m128i*)errors+4, e4);
				_mm_store_si128((__m128i*)errors+5, e5);
				_mm_store_si128((__m128i*)errors+6, e6);
				_mm_store_si128((__m128i*)errors+7, e7);
				PROF(WP_ADD);
				}
				rows[0]+=8;
				rows[1]+=8;
				rows[2]+=8;
				rows[3]+=8;
				
				if(updatewp&&kx&&!(kx&(BLOCKSIZE-1)))
				{
					update_wp(wp, errors);
					wp+=sizeof(__m128i[NWP])/(sizeof(*wp));
					PROF(WP_UPDATE);
				}
#ifdef ENABLE_GUIDE
				curr=rgb;
				if(guide&&memcmp(curr, guide->data+idx, image->nch*sizeof(short)))
				{
					short orig[4]={0};
					memcpy(orig, guide->data+idx, image->nch*sizeof(short));
					curr[0]-=curr[1];
					curr[2]-=curr[1];
					curr[1]+=(curr[0]+curr[2])>>2;
					orig[0]-=orig[1];
					orig[2]-=orig[1];
					orig[1]+=(orig[0]+orig[2])>>2;
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
#ifdef SHOW_PRED_ERRORS
			for(int k=0;k<_countof(errors);++k)
				printf("%2d  %12lld\n", k, total_errors[k]);
#endif
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
	dlist_clear(list+0);
	dlist_clear(list+1);
	dlist_clear(list+2);
	_mm_free(mixin_CDFs);
	free(stats);
	free(pixels);
	_mm_free(weights);
	return 1;
}