#include"codec.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
#ifndef DISABLE_MT
	#define ENABLE_MT
#endif
//	#define AC_DISABLE_WRITE
//	#define PRINT_PREDHIST
//	#define ENABLE_EDGECASES

//OG preds: W/NE/CG/AV4/AV9/WG

	#define ENABLE_W	//good with GDCC
//	#define ENABLE_N	//useless
//	#define ENABLE_NW	//good
	#define ENABLE_NE	//?
//	#define ENABLE_NWE	//bad
//	#define ENABLE_CG45	//bad
//	#define ENABLE_AV2	//good, but redundant with WG
//	#define ENABLE_GRAD	//weaker
	#define ENABLE_CG
//	#define ENABLE_CGv2	//useless
//	#define ENABLE_SELECT	//bad for synthetic due to quantization
//	#define ENABLE_IZ	//weak with CG/AV4
//	#define ENABLE_AV3	//weak with CG/AV4
	#define ENABLE_AV4	//good with GDCC
//	#define ENABLE_AV4v2	//?
//	#define ENABLE_AV5	//redundant with AV4
	#define ENABLE_AV9	//good				GDCC 0.31% slower, 0.37% smaller
//	#define ENABLE_AV9v2	//weak				GDCC 3~5% slower, 0.03% smaller
	#define ENABLE_WG	//for noisy areas		GDCC 5.69% slower, 0.07% smaller
//	#define ENABLE_WG2	//bad
//	#define ENABLE_WG3	//bad

	#define SIMD_CTX	//0.8% faster
//	#define SIMD_PREDS	//3.8% slower

//	#define ENABLE_MIX1	//inefficient			GDCC 9.7~6% faster, 0.6% larger
//	#define ENABLE_MIX2	//bad
//	#define ENABLE_MIX3	//slow				GDCC 16% slower, 0.2% smaller
//	#define ENABLE_MIX4	//slow				GDCC 28% slower, 0.3% smaller

//	#define AC3_ENC_BRANCHLESSRENORM	//slow
//	#define AC3_DEC_BRANCHLESSRENORM	//insignificantly (0.2%) faster
//	#define BYPASS_ENC_BRANCHLESSRENORM	//insignificantly (0.4%) faster
//	#define BYPASS_DEC_BRANCHLESSRENORM	//slow

//PROB_BITS <= 16
//PROB_BITS + PROB_EBITS + 1 <= 31
#define PROB_BITS 16
#define PROB_EBITS 14

#define ANALYSIS_XSTRIDE 2
#define ANALYSIS_YSTRIDE 3

//	#define AC_VALIDATE
//	#define AC_IMPLEMENTATION
#include"entropy.h"

#define BLOCKX 512	//384
#define BLOCKY 512

#define MAXPRINTEDBLOCKS 500
//11<<6 = 704 CDFs size 32
#define ELEVELS 11
#define CBITS 6
#define CLEVELS (1<<CBITS)

#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static void guide_save(const unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(const unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		LOG_ERROR("");
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif
#define OCHLIST\
	OCH(R)\
	OCH(G)\
	OCH(B)\
	OCH(RG)\
	OCH(RB)\
	OCH(GB)\
	OCH(GR)\
	OCH(BR)\
	OCH(BG)
typedef enum _OCHIndex
{
#define OCH(LABEL) OCH_##LABEL,
	OCHLIST
#undef  OCH
	OCH_COUNT,
} OCHIndex;
//static const char *och_names[OCH_COUNT]=
//{
//#define OCH(LABEL) #LABEL,
//	OCHLIST
//#undef  OCH
//};

#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	3, 3, 3)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	3, 3, 1)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	3, 3, 0)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	3, 3, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	3, 3, 1)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	3, 3, 1)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	3, 3, 0)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	3, 0, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	3, 0, 1)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	3, 0, 1)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	3, 0, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	3, 0, 1)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	3, 0, 1)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	3, 0, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	3, 0, 1)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	3, 0, 1)
typedef enum _RCTIndex
{
#define RCT(LABEL, ...) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][9]=
{
#define RCT(LABEL, ...) {__VA_ARGS__},
	RCTLIST
#undef  RCT
};
//static const char *rct_names[RCT_COUNT]=
//{
//#define RCT(LABEL, ...) #LABEL,
//	RCTLIST
//#undef  RCT
//};

#define PREDLIST\
	PRED(W)\
	PRED(NE)\
	PRED(CG)\
	PRED(AV4)\
	PRED(AV9)\
	PRED(WG)
typedef enum _PredIndex
{
#define PRED(LABEL) PRED_##LABEL,
	PREDLIST
#undef  PRED
	PRED_COUNT,
} PredIndex;
//static const char *pred_names[PRED_COUNT]=
//{
//#define PRED(LABEL) #LABEL,
//	PREDLIST
//#undef  PRED
//};
#ifdef ENABLE_WG
#define MIX2_16x16(DST, A0, A1, C0, C1)\
	do\
	{\
		__m256 _a0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(A0, 16), 16));\
		__m256 _a0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(A0, 16));\
		__m256 _a1lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(A1, 16), 16));\
		__m256 _a1hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(A1, 16));\
		__m256 _c0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(C0, 16), 16));\
		__m256 _c0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(C0, 16));\
		__m256 _c1lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(C1, 16), 16));\
		__m256 _c1hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(C1, 16));\
		_a0lo=_mm256_mul_ps(_a0lo, _c0lo);\
		_a0hi=_mm256_mul_ps(_a0hi, _c0hi);\
		_a1lo=_mm256_mul_ps(_a1lo, _c1lo);\
		_a1hi=_mm256_mul_ps(_a1hi, _c1hi);\
		_c0lo=_mm256_add_ps(_c0lo, _c1lo);\
		_c0hi=_mm256_add_ps(_c0hi, _c1hi);\
		_a0lo=_mm256_add_ps(_a0lo, _a1lo);\
		_a0hi=_mm256_add_ps(_a0hi, _a1hi);\
		_a0lo=_mm256_div_ps(_a0lo, _c0lo);\
		_a0hi=_mm256_div_ps(_a0hi, _c0hi);\
		DST=_mm256_and_si256(_mm256_cvttps_epi32(_a0lo), _mm256_set1_epi32(0xFFFF));\
		DST=_mm256_or_si256(DST, _mm256_slli_epi32(_mm256_cvttps_epi32(_a0hi), 16));\
	}while(0)
#define MIX2_16x8(DST, A0, A1, C0, C1)\
	do\
	{\
		__m128 _a0lo=_mm_cvtepi32_ps(_mm_srai_epi32(_mm_slli_epi32(A0, 16), 16));\
		__m128 _a0hi=_mm_cvtepi32_ps(_mm_srai_epi32(A0, 16));\
		__m128 _a1lo=_mm_cvtepi32_ps(_mm_srai_epi32(_mm_slli_epi32(A1, 16), 16));\
		__m128 _a1hi=_mm_cvtepi32_ps(_mm_srai_epi32(A1, 16));\
		__m128 _c0lo=_mm_cvtepi32_ps(_mm_srai_epi32(_mm_slli_epi32(C0, 16), 16));\
		__m128 _c0hi=_mm_cvtepi32_ps(_mm_srai_epi32(C0, 16));\
		__m128 _c1lo=_mm_cvtepi32_ps(_mm_srai_epi32(_mm_slli_epi32(C1, 16), 16));\
		__m128 _c1hi=_mm_cvtepi32_ps(_mm_srai_epi32(C1, 16));\
		_a0lo=_mm_mul_ps(_a0lo, _c0lo);\
		_a0hi=_mm_mul_ps(_a0hi, _c0hi);\
		_a1lo=_mm_mul_ps(_a1lo, _c1lo);\
		_a1hi=_mm_mul_ps(_a1hi, _c1hi);\
		_c0lo=_mm_add_ps(_c0lo, _c1lo);\
		_c0hi=_mm_add_ps(_c0hi, _c1hi);\
		_a0lo=_mm_add_ps(_a0lo, _a1lo);\
		_a0hi=_mm_add_ps(_a0hi, _a1hi);\
		_a0lo=_mm_div_ps(_a0lo, _c0lo);\
		_a0hi=_mm_div_ps(_a0hi, _c0hi);\
		DST=_mm_and_si128(_mm_cvttps_epi32(_a0lo), _mm_set1_epi32(0xFFFF));\
		DST=_mm_or_si128(DST, _mm_slli_epi32(_mm_cvttps_epi32(_a0hi), 16));\
	}while(0)
#endif
#ifdef ENABLE_WG2
#define MIX2_16x16_LOG2(DST, A0, A1, C0, C1)\
	do\
	{\
		__m256i _a0lo=_mm256_srai_epi32(_mm256_slli_epi32(A0, 16), 16);\
		__m256i _a0hi=_mm256_srai_epi32(A0, 16);\
		__m256i _a1lo=_mm256_srai_epi32(_mm256_slli_epi32(A1, 16), 16);\
		__m256i _a1hi=_mm256_srai_epi32(A1, 16);\
		__m256i _c0lo=FLOOR_LOG2_32x8(_mm256_srai_epi32(_mm256_slli_epi32(C0, 16), 16));\
		__m256i _c0hi=FLOOR_LOG2_32x8(_mm256_srai_epi32(C0, 16));\
		__m256i _c1lo=FLOOR_LOG2_32x8(_mm256_srai_epi32(_mm256_slli_epi32(C1, 16), 16));\
		__m256i _c1hi=FLOOR_LOG2_32x8(_mm256_srai_epi32(C1, 16));\
		__m256i t0lo=_mm256_sub_epi32(_a0lo, _mm256_srav_epi32(_mm256_sub_epi32(_a1lo, _a0lo), _mm256_add_epi32(_mm256_sub_epi32(_c0lo, _c1lo), _mm256_set1_epi32(1))));\
		__m256i t1lo=_mm256_sub_epi32(_a1lo, _mm256_srav_epi32(_mm256_sub_epi32(_a0lo, _a1lo), _mm256_add_epi32(_mm256_sub_epi32(_c1lo, _c0lo), _mm256_set1_epi32(1))));\
		__m256i t0hi=_mm256_sub_epi32(_a0hi, _mm256_srav_epi32(_mm256_sub_epi32(_a1hi, _a0hi), _mm256_add_epi32(_mm256_sub_epi32(_c0hi, _c1hi), _mm256_set1_epi32(1))));\
		__m256i t1hi=_mm256_sub_epi32(_a1hi, _mm256_srav_epi32(_mm256_sub_epi32(_a0hi, _a1hi), _mm256_add_epi32(_mm256_sub_epi32(_c1hi, _c0hi), _mm256_set1_epi32(1))));\
		t0lo=_mm256_blendv_epi8(t0lo, t1lo, _mm256_cmpgt_epi32(_c0lo, _c1lo));\
		t0hi=_mm256_blendv_epi8(t0hi, t1hi, _mm256_cmpgt_epi32(_c0hi, _c1hi));\
		DST=_mm256_or_si256(_mm256_and_si256(t0lo, _mm256_set1_epi32(0xFFFF)), _mm256_slli_epi32(t0hi, 16));\
	}while(0)
#endif
#ifdef ENABLE_WG3
AWM_INLINE void pred_wg3(short *pred, const short *gx, const short *gy, const short *N, const short *W, const short *NW)
{
	for(int k=0;k<16;++k)
	{
		pred[k]=(gx[k]*N[k]+gy[k]*W[k])/(gx[k]+gy[k]);
		pred[k]=pred[k]+(pred[k]-NW[k])*(gx[k]+gy[k])/(768*2);
	}
}
#endif

//from libjxl		sym = packsign(delta) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of FLOOR_LOG2(sym)-E,  bypass = 0bBB...BB  FLOOR_LOG2(sym)-(M+L) bits
#define CONFIG_EXP 4	//revise AVX2 CDF search to change config
#define CONFIG_MSB 1	//410->1+24+1	411->1+32+1	421->1+48+1
#define CONFIG_LSB 1	//		511->1+44+1
AWM_INLINE void quantize_pixel(int val, int *token, int *bypass, int *nbits)
{
	if(val<(1<<CONFIG_EXP))
	{
		*token=val;//token
		*nbits=0;
		*bypass=0;
	}
	else
	{
		int msb=FLOOR_LOG2((unsigned)val);
		int mantissa=val-(1<<msb);
		*token = (1<<CONFIG_EXP) + (
				(msb-CONFIG_EXP)<<(CONFIG_MSB+CONFIG_LSB)|
				(mantissa>>(msb-CONFIG_MSB))<<CONFIG_LSB|
				(mantissa&((1<<CONFIG_LSB)-1))
			);
		*nbits=msb-(CONFIG_MSB+CONFIG_LSB);
		*bypass=val>>CONFIG_LSB&((1LL<<*nbits)-1);
	}
}
typedef struct _ThreadArgs
{
	const unsigned char *src;
	unsigned char *dst;
	int iw, ih;

	int fwd, test, loud, b1, b2, xblocks, nblocks, x1, x2, y1, y2;
	int bufsize, histsize;
	short *pixels;
	int *hist;

	//tokens and bypass are encoded separately:
	// 0      enctokenoffsets[0]      [1]             [2]      [N-2]             [N-1]=tokentotal
	// | tokenstream 0 | tokenstream 1 | tokenstream 2 |   ...   | tokenstream N-1 |
	//
	// 0     encbypassoffsets[0]        [1]              [2]      [N-2]              [N-1]=bypasstotal
	// | bypassstream 0 | bypassstream 1 | bypassstream 2 |   ...   | bypassstream N-1 |
	int *enctokenoffsets, *encbypassoffsets, encbufsize;
	unsigned char *enctokenbuf, *encbypassbuf;
	const int *decoffsets;
	const unsigned char *decstream;

	int tlevels;

	const unsigned *mixincdfs;
	int ctxctrs[3][ELEVELS][CLEVELS];//good

	//aux
	int blockidx;
	double bestsize;
	int bestrct, predidx[3];
} ThreadArgs;
static void block_func(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	AC3 ec;
	const int half=128;
	BypassCoder bc;
	int cdfstride=args->tlevels+1;
	const unsigned char *image=args->fwd?args->src:args->dst;
	unsigned char bestrct=0, predidx[3]={0};
	const unsigned char *combination=0;
	int entropylevel=0, entropyidx=0;

	if(args->fwd)
	{
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		unsigned char predsel[OCH_COUNT]={0};
		int dx=(args->x2-args->x1-3)/ANALYSIS_XSTRIDE/5*5, dy=(args->y2-args->y1-2)/ANALYSIS_YSTRIDE;
		int count=0;
		int ystride=args->iw*3;
		if(dx<=0||dy<=0)
		{
			bestrct=RCT_G_BG_RG;
			combination=rct_combinations[bestrct];
#ifdef ENABLE_CG
			predidx[0]=PRED_CG;
			predidx[1]=PRED_CG;
			predidx[2]=PRED_CG;
#else
			predidx[0]=0;
			predidx[1]=0;
			predidx[2]=0;
#endif
			entropyidx=3;
			goto skip_analysis;
		}
		
		memset(args->hist, 0, args->histsize);
		for(int ky=args->y1+2;ky<args->y2;ky+=ANALYSIS_YSTRIDE)//analysis loop
		{
			int kx=args->x1+2;
			const unsigned char *ptr=image+3*(args->iw*ky+kx);

			__m256i amin=_mm256_set1_epi16(-128);
			__m256i amax=_mm256_set1_epi16(127);
			__m256i amask=_mm256_set1_epi16(255);
			__m128i half8=_mm_set1_epi8(128);
			__m128i shuf=_mm_set_epi8(
				-1,
				12, 14, 13,
				 9, 11, 10,
				 6,  8,  7,
				 3,  5,  4,
				 0,  2,  1
				//15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
			);
		//	__m256i four=_mm256_set1_epi16(4);
		//	__m256i eight=_mm256_set1_epi16(8);
			ALIGN(32) short result[16]={0};
			for(;kx<=args->x2-(5*ANALYSIS_XSTRIDE+1+1);kx+=5*ANALYSIS_XSTRIDE, ptr+=15*ANALYSIS_XSTRIDE, count+=5)
			{
				//		NNW	NN	NNE
				//	NWW	NW	N	NE
				//	WW	W	?
				__m256i
					NNW,	NNW2,	NNW3,	NNW4,
					NN,	NN2,	NN3,	NN4,
					NNE,	NNE2,	NNE3,	NNE4,
					NWW,	NWW2,	NWW3,	NWW4,
					NW,	NW2,	NW3,	NW4,
					N,	N2,	N3,	N4,
					NE,	NE2,	NE3,	NE4,
					WW,	WW2,	WW3,	WW4,
					W,	W2,	W3,	W4,
					curr,	curr2;
				__m256i vmin[3], vmax[3], pred;
				{
					__m128i NNW8	=_mm_loadu_si128((__m128i*)(ptr-2*ystride-1*3+0));
					__m128i NN8	=_mm_loadu_si128((__m128i*)(ptr-2*ystride+0*3+0));
					__m128i NNE8	=_mm_loadu_si128((__m128i*)(ptr-2*ystride+1*3+0));
					__m128i NWW8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride-2*3+0));
					__m128i NW8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride-1*3+0));
					__m128i N8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride+0*3+0));
					__m128i NE8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride+1*3+0));
					__m128i WW8	=_mm_loadu_si128((__m128i*)(ptr+0*ystride-2*3+0));
					__m128i W8	=_mm_loadu_si128((__m128i*)(ptr+0*ystride-1*3+0));
					__m128i curr8	=_mm_loadu_si128((__m128i*)(ptr+0*ystride+0*3+0));
					NNW8	=_mm_xor_si128(NNW8	, half8);
					NN8	=_mm_xor_si128(NN8	, half8);
					NNE8	=_mm_xor_si128(NNE8	, half8);
					NWW8	=_mm_xor_si128(NWW8	, half8);
					NW8	=_mm_xor_si128(NW8	, half8);
					N8	=_mm_xor_si128(N8	, half8);
					NE8	=_mm_xor_si128(NE8	, half8);
					WW8	=_mm_xor_si128(WW8	, half8);
					W8	=_mm_xor_si128(W8	, half8);
					curr8	=_mm_xor_si128(curr8	, half8);
					__m128i NNW82	=_mm_shuffle_epi8(NNW8	, shuf);
					__m128i NN82	=_mm_shuffle_epi8(NN8	, shuf);
					__m128i NNE82	=_mm_shuffle_epi8(NNE8	, shuf);
					__m128i NWW82	=_mm_shuffle_epi8(NWW8	, shuf);
					__m128i NW82	=_mm_shuffle_epi8(NW8	, shuf);
					__m128i N82	=_mm_shuffle_epi8(N8	, shuf);
					__m128i NE82	=_mm_shuffle_epi8(NE8	, shuf);
					__m128i WW82	=_mm_shuffle_epi8(WW8	, shuf);
					__m128i W82	=_mm_shuffle_epi8(W8	, shuf);
					__m128i curr82	=_mm_shuffle_epi8(curr8	, shuf);
					NNW	=_mm256_cvtepi8_epi16(NNW8);
					NN	=_mm256_cvtepi8_epi16(NN8);
					NNE	=_mm256_cvtepi8_epi16(NNE8);
					NWW	=_mm256_cvtepi8_epi16(NWW8);
					NW	=_mm256_cvtepi8_epi16(NW8);
					N	=_mm256_cvtepi8_epi16(N8);
					NE	=_mm256_cvtepi8_epi16(NE8);
					WW	=_mm256_cvtepi8_epi16(WW8);
					W	=_mm256_cvtepi8_epi16(W8);
					curr	=_mm256_cvtepi8_epi16(curr8);
					NNW2	=_mm256_cvtepi8_epi16(NNW82);
					NN2	=_mm256_cvtepi8_epi16(NN82);
					NNE2	=_mm256_cvtepi8_epi16(NNE82);
					NWW2	=_mm256_cvtepi8_epi16(NWW82);
					NW2	=_mm256_cvtepi8_epi16(NW82);
					N2	=_mm256_cvtepi8_epi16(N82);
					NE2	=_mm256_cvtepi8_epi16(NE82);
					WW2	=_mm256_cvtepi8_epi16(WW82);
					W2	=_mm256_cvtepi8_epi16(W82);
					curr2	=_mm256_cvtepi8_epi16(curr82);
					NNW3	=_mm256_sub_epi16(NNW	, NNW2	);
					NN3	=_mm256_sub_epi16(NN	, NN2	);
					NNE3	=_mm256_sub_epi16(NNE	, NNE2	);
					NWW3	=_mm256_sub_epi16(NWW	, NWW2	);
					NW3	=_mm256_sub_epi16(NW	, NW2	);
					N3	=_mm256_sub_epi16(N	, N2	);
					NE3	=_mm256_sub_epi16(NE	, NE2	);
					WW3	=_mm256_sub_epi16(WW	, WW2	);
					W3	=_mm256_sub_epi16(W	, W2	);
				//	curr3	=_mm256_sub_epi16(curr	, curr2	);
					NNW4	=_mm256_sub_epi16(NNW2	, NNW	);
					NN4	=_mm256_sub_epi16(NN2	, NN	);
					NNE4	=_mm256_sub_epi16(NNE2	, NNE	);
					NWW4	=_mm256_sub_epi16(NWW2	, NWW	);
					NW4	=_mm256_sub_epi16(NW2	, NW	);
					N4	=_mm256_sub_epi16(N2	, N	);
					NE4	=_mm256_sub_epi16(NE2	, NE	);
					WW4	=_mm256_sub_epi16(WW2	, WW	);
					W4	=_mm256_sub_epi16(W2	, W	);
				//	curr4	=_mm256_sub_epi16(curr2	, curr	);
				}
#define UPDATE(PREDIDX, IDX0, IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7, IDX8, IDX9, IDXA, IDXB, IDXC, IDXD, IDXE)\
	do\
	{\
		pred=_mm256_sub_epi16(pred, amin);\
		pred=_mm256_and_si256(pred, amask);\
		_mm256_store_si256((__m256i*)result, pred);\
		++args->hist[(IDX0*PRED_COUNT+PREDIDX)<<8|result[0x0]];\
		++args->hist[(IDX1*PRED_COUNT+PREDIDX)<<8|result[0x1]];\
		++args->hist[(IDX2*PRED_COUNT+PREDIDX)<<8|result[0x2]];\
		++args->hist[(IDX3*PRED_COUNT+PREDIDX)<<8|result[0x3]];\
		++args->hist[(IDX4*PRED_COUNT+PREDIDX)<<8|result[0x4]];\
		++args->hist[(IDX5*PRED_COUNT+PREDIDX)<<8|result[0x5]];\
		++args->hist[(IDX6*PRED_COUNT+PREDIDX)<<8|result[0x6]];\
		++args->hist[(IDX7*PRED_COUNT+PREDIDX)<<8|result[0x7]];\
		++args->hist[(IDX8*PRED_COUNT+PREDIDX)<<8|result[0x8]];\
		++args->hist[(IDX9*PRED_COUNT+PREDIDX)<<8|result[0x9]];\
		++args->hist[(IDXA*PRED_COUNT+PREDIDX)<<8|result[0xA]];\
		++args->hist[(IDXB*PRED_COUNT+PREDIDX)<<8|result[0xB]];\
		++args->hist[(IDXC*PRED_COUNT+PREDIDX)<<8|result[0xC]];\
		++args->hist[(IDXD*PRED_COUNT+PREDIDX)<<8|result[0xD]];\
		++args->hist[(IDXE*PRED_COUNT+PREDIDX)<<8|result[0xE]];\
	}while(0)
				vmin[0]=_mm256_min_epi16(N, W);
				vmax[0]=_mm256_max_epi16(N, W);
				vmin[1]=_mm256_min_epi16(N3, W3);
				vmax[1]=_mm256_max_epi16(N3, W3);
				vmin[2]=_mm256_min_epi16(N4, W4);
				vmax[2]=_mm256_max_epi16(N4, W4);
#if defined ENABLE_AV9v2 || defined ENABLE_AV4v2
				__m256i vmin0[3]=
				{
					vmin[0],
					vmin[1],
					vmin[2],
				};
				__m256i vmax0[3]=
				{
					vmax[0],
					vmax[1],
					vmax[2],
				};
#endif
#ifdef ENABLE_AV4v2
				__m256i av4pred[3];
#endif
#ifdef ENABLE_AV9v2
				__m256i av9pred[3];
#endif
#ifdef ENABLE_CGv2
				__m256i cgpred[3];
#endif

				//N
#ifdef ENABLE_N
				pred=_mm256_sub_epi16(curr, N);
				UPDATE(
					PRED_N,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(N3, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_N,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(N4, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_N,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif
				
				//W
#ifdef ENABLE_W
				pred=_mm256_sub_epi16(curr, W);
				UPDATE(
					PRED_W,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(W3, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_W,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(W4, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_W,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif

				//NW
#ifdef ENABLE_NW
				pred=NW;

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_NW,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=NW3;

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_NW,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=NW4;

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_NW,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif

				//NE
#ifdef ENABLE_NE
				pred=NE;

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_NE,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=NE3;

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_NE,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=NE4;

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_NE,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif

				//(NW+NE)/2
#ifdef ENABLE_NWE
				pred=_mm256_add_epi16(NW, NE);
				pred=_mm256_srai_epi16(pred, 1);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_NWE,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(NW3, NE3);
				pred=_mm256_srai_epi16(pred, 1);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_NWE,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(NW4, NE4);
				pred=_mm256_srai_epi16(pred, 1);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_NWE,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif

				//AV2 = (N+W)>>1
#ifdef ENABLE_AV2
				pred=_mm256_srai_epi16(_mm256_add_epi16(N, W), 1);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV2,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_srai_epi16(_mm256_add_epi16(N3, W3), 1);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV2,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_srai_epi16(_mm256_add_epi16(N4, W4), 1);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_AV2,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif

				//grad = N+W-NW
#ifdef ENABLE_GRAD
				pred=_mm256_sub_epi16(_mm256_add_epi16(N, W), NW);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_grad,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_sub_epi16(_mm256_add_epi16(N3, W3), NW3);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_grad,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_sub_epi16(_mm256_add_epi16(N4, W4), NW4);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_grad,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif

				//Select = abs(N-NW)>abs(W-NW) ? N : W
#ifdef ENABLE_SELECT
				pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N, NW)), _mm256_abs_epi16(_mm256_sub_epi16(W, NW)));
				pred=_mm256_blendv_epi8(W, N, pred);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_Paeth2,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N3, NW3)), _mm256_abs_epi16(_mm256_sub_epi16(W3, NW3)));
				pred=_mm256_blendv_epi8(W3, N3, pred);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_Paeth2,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N4, NW4)), _mm256_abs_epi16(_mm256_sub_epi16(W4, NW4)));
				pred=_mm256_blendv_epi8(W4, N4, pred);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_Paeth2,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif

				//CG45 = median(NW, NE, NW+NE-NN)
#ifdef ENABLE_CG45
				vmin[0]=_mm256_min_epi16(NW, NE);
				vmax[0]=_mm256_max_epi16(NW, NE);
				vmin[1]=_mm256_min_epi16(NW3, NE3);
				vmax[1]=_mm256_max_epi16(NW3, NE3);
				vmin[2]=_mm256_min_epi16(NW4, NE4);
				vmax[2]=_mm256_max_epi16(NW4, NE4);

				pred=_mm256_sub_epi16(_mm256_add_epi16(NW, NE), NN);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CG45,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_sub_epi16(_mm256_add_epi16(NW3, NE3), NN3);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CG45,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_sub_epi16(_mm256_add_epi16(NW4, NE4), NN4);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_CG45,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif

				//CG = median(N, W, N+W-NW)
#ifdef ENABLE_CG
				pred=_mm256_sub_epi16(_mm256_add_epi16(N, W), NW);
#ifdef ENABLE_CGv2
				cgpred[0]=pred;
#endif
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CG,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_sub_epi16(_mm256_add_epi16(N3, W3), NW3);
#ifdef ENABLE_CGv2
				cgpred[1]=pred;
#endif
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CG,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_sub_epi16(_mm256_add_epi16(N4, W4), NW4);
#ifdef ENABLE_CGv2
				cgpred[2]=pred;
#endif
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_CG,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif
				vmin[0]=_mm256_min_epi16(vmin[0], NE);
				vmax[0]=_mm256_max_epi16(vmax[0], NE);
				vmin[1]=_mm256_min_epi16(vmin[1], NE3);
				vmax[1]=_mm256_max_epi16(vmax[1], NE3);
				vmin[2]=_mm256_min_epi16(vmin[2], NE4);
				vmax[2]=_mm256_max_epi16(vmax[2], NE4);

				//CGv2 = CLAMP(N+W-NW, N,W,NE)
#ifdef ENABLE_CGv2
				pred=cgpred[0];
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CGv2,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=cgpred[1];
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CGv2,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=cgpred[2];
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_CGv2,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif
				
				//AV3 = clamp((5*(N+W)-2*NW)>>3, N,W,NE)
#ifdef ENABLE_AV3
				pred=_mm256_add_epi16(N, W);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_sub_epi16(pred, _mm256_slli_epi16(NW, 1));
				pred=_mm256_srai_epi16(pred, 3);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV3,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(N3, W3);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_sub_epi16(pred, _mm256_slli_epi16(NW3, 1));
				pred=_mm256_srai_epi16(pred, 3);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV3,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(N4, W4);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_sub_epi16(pred, _mm256_slli_epi16(NW4, 1));
				pred=_mm256_srai_epi16(pred, 3);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_AV3,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif

				//IZ = clamp((3*(N+W)-2*NW)>>2, N,W,NE)
#ifdef ENABLE_IZ
				pred=_mm256_add_epi16(N, W);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(_mm256_sub_epi16(pred, NW), 1));
				pred=_mm256_srai_epi16(pred, 2);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_IZ,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(N3, W3);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(_mm256_sub_epi16(pred, NW3), 1));
				pred=_mm256_srai_epi16(pred, 2);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_IZ,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(N4, W4);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(_mm256_sub_epi16(pred, NW4), 1));
				pred=_mm256_srai_epi16(pred, 2);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_IZ,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif

				//AV4 = clamp((4*(N+W)+NE-NW)>>3, N,W,NE)
#ifdef ENABLE_AV4
				pred=_mm256_add_epi16(N, W);
				pred=_mm256_slli_epi16(pred, 2);
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE, NW));
				pred=_mm256_srai_epi16(pred, 3);
#ifdef ENABLE_AV4v2
				av4pred[0]=pred;
#endif
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV4,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(N3, W3);
				pred=_mm256_slli_epi16(pred, 2);
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE3, NW3));
				pred=_mm256_srai_epi16(pred, 3);
#ifdef ENABLE_AV4v2
				av4pred[1]=pred;
#endif
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV4,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(N4, W4);
				pred=_mm256_slli_epi16(pred, 2);
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE4, NW4));
				pred=_mm256_srai_epi16(pred, 3);
#ifdef ENABLE_AV4v2
				av4pred[2]=pred;
#endif
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_AV4,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif

				//AV4 = clamp((4*(N+W)+NE-NW)>>3, N,W,NE)
#ifdef ENABLE_AV4v2
				pred=av4pred[0];
				pred=_mm256_max_epi16(pred, vmin0[0]);
				pred=_mm256_min_epi16(pred, vmax0[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV4v2,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=av4pred[1];
				pred=_mm256_max_epi16(pred, vmin0[1]);
				pred=_mm256_min_epi16(pred, vmax0[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV4v2,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=av4pred[2];
				pred=_mm256_max_epi16(pred, vmin0[2]);
				pred=_mm256_min_epi16(pred, vmax0[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_AV4v2,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif
				
				//AV5 = clamp(W+((5*(N-NW)+NE-WW)>>3), N,W,NE)
#ifdef ENABLE_AV5
				pred=_mm256_sub_epi16(N, NW);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE, WW));
			//	pred=_mm256_add_epi16(pred, four);
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), W);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV5,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_sub_epi16(N3, NW3);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE3, WW3));
			//	pred=_mm256_add_epi16(pred, four);
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), W3);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV5,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_sub_epi16(N4, NW4);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE4, WW4));
			//	pred=_mm256_add_epi16(pred, four);
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), W4);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_AV5,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif
				
				//AV9
				//		1	-2	-1
				//	-1	-9	10	4
				//	-2	16	[?]>>4		clamp(N,W,NE)
#ifdef ENABLE_AV9
				pred=_mm256_add_epi16(N, _mm256_slli_epi16(N, 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(NN, WW));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(NE, 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(NW, 3), NW));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NNW, _mm256_add_epi16(NNE, NWW)));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
			//	pred=_mm256_add_epi16(pred, eight);
				pred=_mm256_add_epi16(W, _mm256_srai_epi16(pred, 4));
#ifdef ENABLE_AV9v2
				av9pred[0]=pred;
#endif
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV9,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(N3, _mm256_slli_epi16(N3, 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(NN3, WW3));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(NE3, 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(NW3, 3), NW3));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NNW3, _mm256_add_epi16(NNE3, NWW3)));
			//	pred=_mm256_add_epi16(pred, eight);
				pred=_mm256_add_epi16(W3, _mm256_srai_epi16(pred, 4));
#ifdef ENABLE_AV9v2
				av9pred[1]=pred;
#endif
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV9,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(N4, _mm256_slli_epi16(N4, 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(NN4, WW4));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(NE4, 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(NW4, 3), NW4));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NNW4, _mm256_add_epi16(NNE4, NWW4)));
			//	pred=_mm256_add_epi16(pred, eight);
				pred=_mm256_add_epi16(W4, _mm256_srai_epi16(pred, 4));
#ifdef ENABLE_AV9v2
				av9pred[2]=pred;
#endif
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_AV9,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif
				
				//AV9v2
				//		1	-2	-1
				//	-1	-9	10	4
				//	-2	16	[?]>>4		clamp(N,W)
#ifdef ENABLE_AV9v2
				pred=av9pred[0];
				//pred=_mm256_add_epi16(N, _mm256_slli_epi16(N, 2));//5*N
				//pred=_mm256_sub_epi16(pred, _mm256_add_epi16(NN, WW));//5*N - (NN+WW)
				//pred=_mm256_add_epi16(pred, _mm256_slli_epi16(NE, 1));//5*N-NN-WW + 2*NE
				//pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(NW, 3), NW));//2*(5*N-NN-WW+2*NE) - 9*NW
				//pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NNW, _mm256_add_epi16(NNE, NWW)));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				//pred=_mm256_add_epi16(W, _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin0[0]);
				pred=_mm256_min_epi16(pred, vmax0[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV9v2,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=av9pred[1];
				//pred=_mm256_add_epi16(N3, _mm256_slli_epi16(N3, 2));
				//pred=_mm256_sub_epi16(pred, _mm256_add_epi16(NN3, WW3));
				//pred=_mm256_add_epi16(pred, _mm256_slli_epi16(NE3, 1));
				//pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(NW3, 3), NW3));
				//pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NNW3, _mm256_add_epi16(NNE3, NWW3)));
				//pred=_mm256_add_epi16(W3, _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin0[1]);
				pred=_mm256_min_epi16(pred, vmax0[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV9v2,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=av9pred[2];
				//pred=_mm256_add_epi16(N4, _mm256_slli_epi16(N4, 2));
				//pred=_mm256_sub_epi16(pred, _mm256_add_epi16(NN4, WW4));
				//pred=_mm256_add_epi16(pred, _mm256_slli_epi16(NE4, 1));
				//pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(NW4, 3), NW4));
				//pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NNW4, _mm256_add_epi16(NNE4, NWW4)));
				//pred=_mm256_add_epi16(W4, _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin0[2]);
				pred=_mm256_min_epi16(pred, vmax0[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_AV9v2,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif
				
				//WG
				//gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1
				//gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1
				//pred=(gx*N+gy*W)/(gx+gy)
#ifdef ENABLE_WG
				__m256i gx, gy;
			//	gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W, WW)), 1);
				gx=_mm256_abs_epi16(_mm256_sub_epi16(W, WW));
				gy=_mm256_abs_epi16(_mm256_sub_epi16(W, NW));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N, NW)));
			//	gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N, NN)), 1));
				gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(N, NN)));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE, N)));
				gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE, NNE)));
				gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
				gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
				MIX2_16x16(pred, N, W, gx, gy);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_WG,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
			//	gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W3, WW3)), 1);
				gx=_mm256_abs_epi16(_mm256_sub_epi16(W3, WW3));
				gy=_mm256_abs_epi16(_mm256_sub_epi16(W3, NW3));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N3, NW3)));
			//	gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N3, NN3)), 1));
				gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(N3, NN3)));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE3, N3)));
				gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE3, NNE3)));
				gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
				gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
				MIX2_16x16(pred, N3, W3, gx, gy);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_WG,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
			//	gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W4, WW4)), 1);
				gx=_mm256_abs_epi16(_mm256_sub_epi16(W4, WW4));
				gy=_mm256_abs_epi16(_mm256_sub_epi16(W4, NW4));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N4, NW4)));
			//	gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N4, NN4)), 1));
				gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(N4, NN4)));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE4, N4)));
				gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE4, NNE4)));
				gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
				gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
				MIX2_16x16(pred, N4, W4, gx, gy);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_WG,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif
				
				//WG2
				//gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1
				//gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1
				//pred = gx<gy ? W+((N-W)>>(gy-gx+1)) : N-((W-N)>>(gx-gy+1))
#ifdef ENABLE_WG2
				__m256i gx, gy;
				gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W, WW)), 1);
				gy=_mm256_abs_epi16(_mm256_sub_epi16(W, NW));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N, NW)));
				gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N, NN)), 1));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE, N)));
				gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE, NNE)));
				gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
				gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
				MIX2_16x16_LOG2(pred, N, W, gx, gy);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_WG2,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W3, WW3)), 1);
				gy=_mm256_abs_epi16(_mm256_sub_epi16(W3, NW3));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N3, NW3)));
				gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N3, NN3)), 1));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE3, N3)));
				gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE3, NNE3)));
				gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
				gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
				MIX2_16x16_LOG2(pred, N3, W3, gx, gy);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_WG2,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W4, WW4)), 1);
				gy=_mm256_abs_epi16(_mm256_sub_epi16(W4, NW4));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N4, NW4)));
				gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N4, NN4)), 1));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE4, N4)));
				gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE4, NNE4)));
				gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
				gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
				MIX2_16x16_LOG2(pred, N4, W4, gx, gy);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_WG2,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif
				
				//WG3
				//gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1
				//gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1
				//pred=(gx*N+gy*W)/(gx+gy);
				//pred=pred+(pred-NW)*(gx+gy)/(2*gmax);
#ifdef ENABLE_WG3
				__m256i gx, gy;
				gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W, WW)), 1);
				gy=_mm256_abs_epi16(_mm256_sub_epi16(W, NW));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N, NW)));
				gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N, NN)), 1));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE, N)));
				gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE, NNE)));
				gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
				gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
				pred_wg3((short*)&pred, (short*)&gx, (short*)&gy, (short*)&N, (short*)&W, (short*)&NW);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_WG3,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W3, WW3)), 1);
				gy=_mm256_abs_epi16(_mm256_sub_epi16(W3, NW3));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N3, NW3)));
				gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N3, NN3)), 1));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE3, N3)));
				gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE3, NNE3)));
				gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
				gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
				pred_wg3((short*)&pred, (short*)&gx, (short*)&gy, (short*)&N3, (short*)&W3, (short*)&NW3);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_WG3,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W4, WW4)), 1);
				gy=_mm256_abs_epi16(_mm256_sub_epi16(W4, NW4));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N4, NW4)));
				gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N4, NN4)), 1));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE4, N4)));
				gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE4, NNE4)));
				gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
				gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
				pred_wg3((short*)&pred, (short*)&gx, (short*)&gy, (short*)&N4, (short*)&W4, (short*)&NW4);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_WG3,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#endif
			}
		}
		double gain=1./count;
		for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)
		{
			int *curr_hist=args->hist+((size_t)kc<<8);
			for(int ks=0;ks<256;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
					csizes[kc]-=freq*log2((double)freq*gain);
			}
			csizes[kc]/=8;
		}
		for(int kc=0;kc<OCH_COUNT;++kc)//select best predictors
		{
			int bestpred=0;
			for(int kp=1;kp<PRED_COUNT;++kp)
			{
				if(csizes[kc*PRED_COUNT+bestpred]>csizes[kc*PRED_COUNT+kp])
					bestpred=kp;
			}
			predsel[kc]=bestpred;
		}
		for(int kt=0;kt<RCT_COUNT;++kt)//select best RCT
		{
			const unsigned char *group=rct_combinations[kt];
			double csize=
				csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
				csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
				csizes[group[2]*PRED_COUNT+predsel[group[2]]];
			if(!kt||bestsize>csize)
				bestsize=csize, bestrct=kt;
		}
		combination=rct_combinations[bestrct];
		predidx[0]=predsel[combination[0]];
		predidx[1]=predsel[combination[1]];
		predidx[2]=predsel[combination[2]];
#if 0
		if(csizes[combination[0]*PRED_COUNT+predidx[0]]<0.1*count)predidx[0]=PRED_CG;	//X
		if(csizes[combination[1]*PRED_COUNT+predidx[1]]<0.1*count)predidx[1]=PRED_CG;
		if(csizes[combination[2]*PRED_COUNT+predidx[2]]<0.1*count)predidx[2]=PRED_CG;
#endif
		args->bestsize=bestsize;
		entropylevel=(int)((
			csizes[combination[0]*PRED_COUNT+predidx[0]]+
			csizes[combination[1]*PRED_COUNT+predidx[1]]+
			csizes[combination[2]*PRED_COUNT+predidx[2]]
		)*gain*(100/3.));//percent
		CLAMP2(entropylevel, 0, 99);
	skip_analysis:
		{
			int prevblockidx=args->blockidx-args->b1-1;
			ac3_encbuf_init(&ec, args->enctokenbuf+(prevblockidx>=0?args->enctokenoffsets[prevblockidx]:0), args->enctokenbuf+args->encbufsize);
			bypass_encbuf_init(&bc, args->encbypassbuf+(prevblockidx>=0?args->encbypassoffsets[prevblockidx]:0), args->encbypassbuf+args->encbufsize);
		}
		ac3_encbuf_bypass_NPOT(&ec, entropylevel, 100);
		ac3_encbuf_bypass_NPOT(&ec, bestrct, RCT_COUNT);
		ac3_encbuf_bypass_NPOT(&ec, predidx[0], PRED_COUNT);
		ac3_encbuf_bypass_NPOT(&ec, predidx[1], PRED_COUNT);
		ac3_encbuf_bypass_NPOT(&ec, predidx[2], PRED_COUNT);

		args->bestrct=bestrct;
		args->predidx[0]=predidx[0];
		args->predidx[1]=predidx[1];
		args->predidx[2]=predidx[2];
#if 0
		if(args->loud)
		{
			printf("Y %5d~%5d  best %12.2lf bytes  %s [YUV: %s %s %s]\n",
				args->y1, args->y2,
				bestsize,
				rct_names[bestrct],
				pred_names[predidx[0]],
				pred_names[predidx[1]],
				pred_names[predidx[2]]
			);

			for(int kp=0;kp<PRED_COUNT;++kp)
				printf("%14s", pred_names[kp]);
			printf("\n");
			for(int kc=0;kc<OCH_COUNT;++kc)
			{
				for(int kp=0;kp<PRED_COUNT;++kp)
					printf(" %12.2lf%c", csizes[kc*PRED_COUNT+kp], kp==predsel[kc]?'*':' ');
				printf("  %s\n", och_names[kc]);
			}

			for(int kt=0;kt<RCT_COUNT;++kt)
			{
				const unsigned char *group=rct_combinations[kt];
				double csize=
					csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
					csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
					csizes[group[2]*PRED_COUNT+predsel[group[2]]];
				printf("%12.2lf %c  %-10s %-10s %-10s %-10s\n",
					csize,
					kt==bestrct?'*':' ',
					rct_names[kt],
					pred_names[predsel[group[0]]],
					pred_names[predsel[group[1]]],
					pred_names[predsel[group[2]]]
				);
			}
		}
#endif
	}
	else//decode header
	{
		ac3_dec_init(&ec, args->decstream+(args->blockidx>0?args->decoffsets[args->blockidx-1]:0), args->decstream+args->decoffsets[args->blockidx]);
		bypass_decbuf_init(&bc, args->decstream+args->decoffsets[args->nblocks+args->blockidx-1], args->decstream+args->decoffsets[args->nblocks+args->blockidx]);
		entropylevel=ac3_dec_bypass_NPOT(&ec, 100);
		bestrct=ac3_dec_bypass_NPOT(&ec, RCT_COUNT);
		combination=rct_combinations[bestrct];
		predidx[0]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[1]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[2]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
	}
	if(entropylevel>25)		entropyidx=3;
	else if(entropylevel>17)	entropyidx=2;
	else if(entropylevel>10)	entropyidx=1;
	else				entropyidx=0;
#if 0
	predidx[0]=PRED_CG;//synth2 0.8% smaller
	predidx[1]=PRED_CG;
	predidx[2]=PRED_CG;
	//predidx[0]=PRED_WG;//GDCC 0.9% larger
	//predidx[1]=PRED_WG;
	//predidx[2]=PRED_WG;
#endif
	int shiftgain=(entropylevel*entropylevel>>4)+15;
	unsigned *histptr=(unsigned*)args->hist;
	for(int ks=0;ks<cdfstride;++ks)
	{
		histptr[ks]=(unsigned)((ks*((1LL<<PROB_BITS)-args->tlevels)<<PROB_EBITS)/args->tlevels);
#ifdef _DEBUG
		if(ks&&histptr[ks-1]>histptr[ks])
			LOG_ERROR("Init error");
#endif
	}
	//histptr[cdfstride-1]=((1LL<<PROB_BITS)-args->tlevels)<<PROB_EBITS;
	memfill(histptr+cdfstride, histptr, args->histsize-cdfstride*sizeof(int), cdfstride*sizeof(int));
	*(int*)args->ctxctrs=1;
	memfill((int*)args->ctxctrs+1, args->ctxctrs, sizeof(args->ctxctrs)-sizeof(int), sizeof(int));
#ifdef SIMD_PREDS
#ifdef ENABLE_NE
	__m128i predmaskNE	=_mm_set_epi16(0, 0, 0, 0, 0,
		-(predidx[2]==PRED_NE),
		-(predidx[1]==PRED_NE),
		-(predidx[0]==PRED_NE)
	);
#endif
#ifdef ENABLE_CG
	int enableCG=predidx[2]==PRED_CG||predidx[1]==PRED_CG||predidx[0]==PRED_CG;
	__m128i predmaskCG	=_mm_set_epi16(0, 0, 0, 0, 0,
		-(predidx[2]==PRED_CG),
		-(predidx[1]==PRED_CG),
		-(predidx[0]==PRED_CG)
	);
#endif
#ifdef ENABLE_AV4
	int enableAV4=predidx[2]==PRED_AV4||predidx[1]==PRED_AV4||predidx[0]==PRED_AV4;
	__m128i predmaskAV4	=_mm_set_epi16(0, 0, 0, 0, 0,
		-(predidx[2]==PRED_AV4),
		-(predidx[1]==PRED_AV4),
		-(predidx[0]==PRED_AV4)
	);
#endif
#ifdef ENABLE_AV9
	int enableAV9=predidx[2]==PRED_AV9||predidx[1]==PRED_AV9||predidx[0]==PRED_AV9;
	__m128i predmaskAV9	=_mm_set_epi16(0, 0, 0, 0, 0,
		-(predidx[2]==PRED_AV9),
		-(predidx[1]==PRED_AV9),
		-(predidx[0]==PRED_AV9)
	);
#endif
#ifdef ENABLE_WG
	int enableWG=predidx[2]==PRED_WG||predidx[1]==PRED_WG||predidx[0]==PRED_WG;
	__m128i predmaskWG	=_mm_set_epi16(0, 0, 0, 0, 0,
		-(predidx[2]==PRED_WG),
		-(predidx[1]==PRED_WG),
		-(predidx[0]==PRED_WG)
	);
#endif
#endif
	int paddedblockwidth=args->x2-args->x1+16;
	memset(args->pixels, 0, args->bufsize);
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(32) short *rows[]=
		{
			args->pixels+(paddedblockwidth*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+(paddedblockwidth*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+(paddedblockwidth*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+(paddedblockwidth*((ky-3LL)&3)+8LL)*4*2,
		};
		int yuv[4]={0};
		int *helpers[3]=
		{
			yuv+combination[6+0],
			yuv+combination[6+1],
			yuv+combination[6+2],
		};
		int token=0, bypass=0, nbits=0;
		int pred=0, error=0, sym=0;
		int yidx=3*(args->iw*ky+args->x1)+combination[3+0];
		int uidx=3*(args->iw*ky+args->x1)+combination[3+1];
		int vidx=3*(args->iw*ky+args->x1)+combination[3+2];
		ALIGN(16) short regW[12]={0};
		ALIGN(16) short preds[8]={0};
		ALIGN(16) int grads[8]={0};
		//ALIGN(16) char apred8[16];
		ALIGN(32) int CDF[8*6]={0};//{8 zeros, CDF size 8*4, 8 CDFmax}
		CDF[8*5+0]=1<<PROB_BITS;
		CDF[8*5+1]=1<<PROB_BITS;
		CDF[8*5+2]=1<<PROB_BITS;
		CDF[8*5+3]=1<<PROB_BITS;
		CDF[8*5+4]=1<<PROB_BITS;
		CDF[8*5+5]=1<<PROB_BITS;
		CDF[8*5+6]=1<<PROB_BITS;
		CDF[8*5+7]=1<<PROB_BITS;
		int *ctxctrs[12]={0};
		unsigned *cdfptrs[12]={0};
		const unsigned *mixinptrs[3]={0};
#if defined SIMD_CTX || defined SIMD_PREDS
		__m128i
			mNNW	=_mm_setzero_si128(),
			mNN	=_mm_loadu_si128((__m128i*)(rows[2]+0*4*2)),
			mNNE	=_mm_loadu_si128((__m128i*)(rows[2]+1*4*2)),
			mNWW	=_mm_setzero_si128(),
			mNW	=_mm_setzero_si128(),
			mN	=_mm_loadu_si128((__m128i*)(rows[1]+0*4*2)),
			mNE	=_mm_loadu_si128((__m128i*)(rows[1]+1*4*2)),
			mNEE	=_mm_loadu_si128((__m128i*)(rows[1]+2*4*2)),
			mWW	=_mm_setzero_si128(),
			mW	=_mm_setzero_si128();
		(void)mNNW;
		(void)mNWW;
#endif
		for(int kx=args->x1;kx<args->x2;++kx, yidx+=3, uidx+=3, vidx+=3)
		{
#ifdef SIMD_CTX
			__m128i msum=_mm_add_epi16(mN, mW);
			//__m128i mqp=_mm_add_epi16(msum, _mm_slli_epi16(_mm_sub_epi16(msum, mNW), 1));
			//mqp=_mm_srai_epi16(mqp, 11-CBITS);
			//mqp=_mm_and_si128(mqp, _mm_set1_epi16((1<<CBITS)-1));

			__m128i mqe=_mm_add_epi16(mNW, mNE);
			mqe=_mm_add_epi16(mqe, _mm_add_epi16(mWW, mNN));
			mqe=_mm_srli_epi16(_mm_add_epi16(mqe, mNEE), 1);
			mqe=_mm_add_epi16(mqe, msum);
			mqe=_mm_shuffle_epi32(mqe, _MM_SHUFFLE(1, 0, 3, 2));
			mqe=_mm_add_epi16(mqe, _mm_abs_epi16(_mm_sub_epi16(mN, mW)));

			//mqp=_mm_cvtepi16_epi32(mqp);
			//_mm_store_si128((__m128i*)grads+1, mqp);

			mqe=_mm_add_epi16(mqe, _mm_set1_epi16(1));
			mqe=FLOOR_LOG2_32x4(_mm_cvtepi16_epi32(mqe));//7 cycles vs (_lzcnt_u32 3 cycles)*3 (actually 5 cycles)		2.45%
			mqe=_mm_min_epi32(mqe, _mm_set1_epi32(ELEVELS-1));
			_mm_store_si128((__m128i*)grads+0, mqe);
#endif
#ifndef SIMD_PREDS
			short
				*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NWW	=rows[1]-2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
			//	*NEE	=rows[1]+2*4*2,
				*WW	=rows[0]-2*4*2,
			//	*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 3
#endif
			for(int kc=0;kc<3;++kc)
			{
				switch(predidx[kc])
				{
#ifdef ENABLE_N
				case PRED_N:
					pred=N[kc];
					break;
#endif
#ifdef ENABLE_W
				case PRED_W:
					pred=regW[kc];
					break;
#endif
#ifdef ENABLE_AV2
				case PRED_AV2:
					pred=(N[kc]+regW[kc])>>1;
					break;
#endif
#ifdef ENABLE_GRAD
				case PRED_grad:
					pred=N[kc]+regW[kc]-NW[kc];
					break;
#endif
#ifdef ENABLE_SELECT
				case PRED_Paeth2:
					pred=abs(regW[kc]-NW[kc])<=abs(N[kc]-NW[kc])?N[kc]:regW[kc];
					break;
#endif
#ifdef ENABLE_NW
				case PRED_NW:
					pred=NW[kc];
					break;
#endif
#ifdef ENABLE_NE
				case PRED_NE:
					pred=NE[kc];
					break;
#endif
#ifdef ENABLE_NWE
				case PRED_NWE:
					pred=(NW[kc]+NE[kc])>>1;
					break;
#endif
#ifdef ENABLE_CG45
				case PRED_CG45:
					MEDIAN3_32(pred, NW[kc], NE[kc], NW[kc]+NE[kc]-NN[kc]);
					break;
#endif
#ifdef ENABLE_CG
				case PRED_CG:
					MEDIAN3_32(pred, N[kc], regW[kc], N[kc]+regW[kc]-NW[kc]);
					//{
					//	int vmax=N[kc], vmin=regW[kc];
					//	if(N[kc]<regW[kc])
					//		vmin=N[kc], vmax=regW[kc];
					//	pred=N[kc]+regW[kc]-NW[kc];
					//	CLAMP2(pred, vmin, vmax);
					//}
					break;
#endif
#ifdef ENABLE_CGv2
				case PRED_CGv2:
					CLAMP3_32(pred, N[kc]+regW[kc]-NW[kc], N[kc], regW[kc], NE[kc]);
					break;
#endif
#ifdef ENABLE_AV3
				case PRED_AV3:
					CLAMP3_32(pred,
						(5*(N[kc]+regW[kc])-2*NW[kc])>>3,
						N[kc], regW[kc], NE[kc]
					);
					break;
#endif
#ifdef ENABLE_IZ
				case PRED_IZ:
					{
						int sum=N[kc]+regW[kc];
						CLAMP3_32(pred,
							(2*(sum-NW[kc])+sum)>>2,
							N[kc], regW[kc], NE[kc]
						);
					}
					break;
#endif
#ifdef ENABLE_AV4
				case PRED_AV4:
					CLAMP3_32(pred,
						(4*(N[kc]+regW[kc])+NE[kc]-NW[kc])>>3,
						N[kc], regW[kc], NE[kc]
					);
					//{
					//	int vmax=N[kc], vmin=regW[kc];
					//	if(N[kc]<regW[kc])
					//		vmin=N[kc], vmax=regW[kc];
					//	if(vmin>NE[kc])vmin=NE[kc];
					//	if(vmax<NE[kc])vmax=NE[kc];
					//	pred=(4*(N[kc]+regW[kc])+NE[kc]-NW[kc])>>3;
					//	CLAMP2(pred, vmin, vmax);
					//}
					break;
#endif
#ifdef ENABLE_AV4v2
				case PRED_AV4v2:
					MEDIAN3_32(pred,
						(4*(N[kc]+regW[kc])+NE[kc]-NW[kc])>>3,
						N[kc], regW[kc]
					);
					break;
#endif
#ifdef ENABLE_AV5
				case PRED_AV5:
					CLAMP3_32(pred,
						regW[kc]+((5*(N[kc]-NW[kc])+NE[kc]-WW[kc])>>3),
						N[kc], regW[kc], NE[kc]
					);
					break;
#endif
#ifdef ENABLE_AV9
				case PRED_AV9:
					CLAMP3_32(pred,
						regW[kc]+((10*N[kc]-9*NW[kc]+4*NE[kc]-2*(NN[kc]+WW[kc])+NNW[kc]-(NNE[kc]+NWW[kc]))>>4),
						N[kc], regW[kc], NE[kc]
					);
					//{
					//	int vmax=N[kc], vmin=regW[kc];
					//	if(N[kc]<regW[kc])
					//		vmin=N[kc], vmax=regW[kc];
					//	if(vmin>NE[kc])vmin=NE[kc];
					//	if(vmax<NE[kc])vmax=NE[kc];
					//	pred=regW[kc]+((10*N[kc]-9*NW[kc]+4*NE[kc]-2*(NN[kc]+WW[kc])+NNW[kc]-(NNE[kc]+NWW[kc]))>>4);
					//	CLAMP2(pred, vmin, vmax);
					//}
					break;
#endif
#ifdef ENABLE_AV9v2
				case PRED_AV9v2:
					MEDIAN3_32(pred,
						regW[kc]+((10*N[kc]-9*NW[kc]+4*NE[kc]-2*(NN[kc]+WW[kc])+NNW[kc]-(NNE[kc]+NWW[kc]))>>4),
						N[kc], regW[kc]
					);
					break;
#endif
#ifdef ENABLE_WG
				case PRED_WG:
					{
						int gx=abs(regW[kc]-WW[kc])+abs(N[kc]-NW[kc])+abs(NE[kc]-N[kc])+1;//don't add errors
						int gy=abs(regW[kc]-NW[kc])+abs(N[kc]-NN[kc])+abs(NE[kc]-NNE[kc])+1;
						pred=(gx*N[kc]+gy*regW[kc])/(gx+gy);
					}
					break;
#endif
#ifdef ENABLE_WG2
				case PRED_WG2:
					{
						int gx=abs(regW[kc]-WW[kc])+abs(N[kc]-NW[kc])+abs(NE[kc]-N[kc])+1;//1~769
						int gy=abs(regW[kc]-NW[kc])+abs(N[kc]-NN[kc])+abs(NE[kc]-NNE[kc])+1;
						gx=FLOOR_LOG2(gx);
						gy=FLOOR_LOG2(gy);
						pred=gx<gy ? regW[kc]+((N[kc]-regW[kc])>>(gy-gx+1)) : N[kc]-((regW[kc]-N[kc])>>(gx-gy+1));
					}
					break;
#endif
#ifdef ENABLE_WG3
				case PRED_WG3:
					{
						int gx=abs(regW[kc]-WW[kc])+abs(N[kc]-NW[kc])+abs(NE[kc]-N[kc])+1;//1~769
						int gy=abs(regW[kc]-NW[kc])+abs(N[kc]-NN[kc])+abs(NE[kc]-NNE[kc])+1;
						pred=(gx*N[kc]+gy*regW[kc])/(gx+gy);
						pred=pred+(pred-NW[kc])*(gx+gy)/(2*768);
					}
					break;
#endif
				}
				preds[kc]=pred;
			}
#else
			short *curr=rows[0];
			//if(ky==0&&kx==1)//
			//	printf("");

			//predictors
			__m128i vmin=_mm_min_epi16(mN, mW);
			__m128i vmax=_mm_max_epi16(mN, mW);
#ifdef ENABLE_CG
			__m128i predCG=_mm_setzero_si128();
			if(enableCG)
			{
				predCG=_mm_sub_epi16(msum, mNW);
				predCG=_mm_max_epi16(predCG, vmin);
				predCG=_mm_min_epi16(predCG, vmax);
			}
#endif
			vmin=_mm_min_epi16(vmin, mNE);
			vmax=_mm_max_epi16(vmax, mNE);
#ifdef ENABLE_AV4
			__m128i predAV4=_mm_setzero_si128();
			if(enableAV4)
			{
				predAV4=_mm_slli_epi16(msum, 2);
				predAV4=_mm_add_epi16(predAV4, _mm_sub_epi16(mNE, mNW));
				predAV4=_mm_srai_epi16(predAV4, 3);
				predAV4=_mm_max_epi16(predAV4, vmin);
				predAV4=_mm_min_epi16(predAV4, vmax);
			}
#endif
#ifdef ENABLE_AV9
			__m128i predAV9=_mm_setzero_si128();
			if(enableAV9)
			{
				predAV9=_mm_add_epi16(mN, _mm_slli_epi16(mN, 2));//5*N
				predAV9=_mm_sub_epi16(predAV9, _mm_add_epi16(mNN, mWW));//5*N - (NN+WW)
				predAV9=_mm_add_epi16(predAV9, _mm_slli_epi16(mNE, 1));//5*N-NN-WW + 2*NE
				predAV9=_mm_sub_epi16(_mm_slli_epi16(predAV9, 1), _mm_add_epi16(_mm_slli_epi16(mNW, 3), mNW));//2*(5*N-NN-WW+2*NE) - 9*NW
				predAV9=_mm_add_epi16(predAV9, _mm_sub_epi16(mNNW, _mm_add_epi16(mNNE, mNWW)));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				//predAV9=_mm_add_epi16(predAV9, eight);
				predAV9=_mm_add_epi16(mW, _mm_srai_epi16(predAV9, 4));
				predAV9=_mm_max_epi16(predAV9, vmin);
				predAV9=_mm_min_epi16(predAV9, vmax);
			}
#endif
#ifdef ENABLE_WG
			__m128i predWG=_mm_setzero_si128();
			if(enableWG)
			{
				__m128i gx, gy;
			//	gx=_mm_slli_epi16(_mm_abs_epi16(_mm_sub_epi16(mW, mWW)), 1);
				gx=_mm_abs_epi16(_mm_sub_epi16(mW, mWW));
				gy=_mm_abs_epi16(_mm_sub_epi16(mW, mNW));
				gx=_mm_add_epi16(gx, _mm_abs_epi16(_mm_sub_epi16(mN, mNW)));
			//	gy=_mm_add_epi16(gy, _mm_slli_epi16(_mm_abs_epi16(_mm_sub_epi16(mN, mNN)), 1));
				gy=_mm_add_epi16(gy, _mm_abs_epi16(_mm_sub_epi16(mN, mNN)));
				gx=_mm_add_epi16(gx, _mm_abs_epi16(_mm_sub_epi16(mNE, mN)));
				gy=_mm_add_epi16(gy, _mm_abs_epi16(_mm_sub_epi16(mNE, mNNE)));
				gx=_mm_add_epi16(gx, _mm_set1_epi16(1));
				gy=_mm_add_epi16(gy, _mm_set1_epi16(1));
				MIX2_16x8(predWG, mN, mW, gx, gy);
			}
#endif
			__m128i mpreds=mW;
#ifdef ENABLE_NE
			mpreds=_mm_blendv_epi8(mpreds, mNE, predmaskNE);
#endif
#ifdef ENABLE_CG
			mpreds=_mm_blendv_epi8(mpreds, predCG, predmaskCG);
#endif
#ifdef ENABLE_AV4
			mpreds=_mm_blendv_epi8(mpreds, predAV4, predmaskAV4);
#endif
#ifdef ENABLE_AV9
			mpreds=_mm_blendv_epi8(mpreds, predAV9, predmaskAV9);
#endif
#ifdef ENABLE_WG
			mpreds=_mm_blendv_epi8(mpreds, predWG, predmaskWG);
#endif
			_mm_store_si128((__m128i*)preds, mpreds);
#endif
			if(args->fwd)
			{
				yuv[0]=args->src[yidx]-128;
				yuv[1]=args->src[uidx]-128;
				yuv[2]=args->src[vidx]-128;
			}
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 3
#endif
			for(int kc=0;kc<3;++kc)
			{
				int offset=0;
				pred=preds[kc];
				//if(ky==480&&kx==109&&kc==2)//
				//if(ky==147&&kx==317&&kc==1)//
				//if(ky==1&&kx==455&&kc==0)//
				//if(ky==0&&kx==5&&kc==1)//
				//if(ky==0&&kx==0&&kc==0)//
				//if(ky==747&&kx==2443&&kc==1)//
				//if(ky==0&&kx==25&&kc==2)//
				//if(ky==0&&kx==0&&kc==2)//
				//if(ky==747&&kx==2444&&kc==0)//
				//if(ky==0&&kx==0&&kc==0)//
				//if(ky==7&&kx==342&&kc==2)//
				//if(ky==0&&kx==2&&kc==1)//
				//if(ky==0&&kx==12&&kc==2)//
				//if(ky==0&&kx==13&&kc==0)//
				//if(ky==0&&kx==14&&kc==0)//
				//if(ky==13&&kx==19&&kc==0)//
				//if(ky==12&&kx==17&&kc==1)//
				//if(ky==12&&kx==16&&kc==1)//
				//	printf("");
#if 1
				//int qp;
				int upred;
				__m128i mp=_mm_set_epi32(0, 0, 0, pred);
				__m128i mhalf=_mm_set_epi32(0, 0, 0, 128);
				if(!kc)
				{
					__m128i mhelp=_mm_sub_epi32(mhalf, _mm_abs_epi32(mp));
					upred=_mm_extract_epi32(mhelp, 0);
				}
				else
				{
					__m128i mhelp=_mm_set_epi32(0, 0, 0, offset=*helpers[kc]);
					mp=_mm_add_epi32(mp, mhelp);
					mp=_mm_max_epi32(mp, _mm_set_epi32(0, 0, 0, -128));
					mp=_mm_min_epi32(mp, _mm_set_epi32(0, 0, 0, 127));
					mhelp=_mm_sub_epi32(mhalf, _mm_abs_epi32(mp));
					pred=_mm_extract_epi32(mp, 0);
					upred=_mm_extract_epi32(mhelp, 0);
				}
				//mp=_mm_add_epi32(mp, mhalf);
				//mp=_mm_srli_epi32(mp, 8-CBITS);
				//qp=_mm_extract_epi32(mp, 0);
				(void)half;
#else
				if(kc)
				{
					pred+=offset=*helpers[kc];
					CLAMP2_32(pred, pred, -128, 127);	//MOV-MOV-MOV-MIN-MAX-MOV (6 instr)
				//	CLAMP2(pred, -128, 127);		//MOV-CMP-CMOV (6 instr) slower?
				}
				int upred=half-abs(pred);
#endif
#ifdef SIMD_CTX
				int qp=(pred+128)>>(8-CBITS);
				//int qe=FLOOR_LOG2(grads[kc+0]+1);
				//if(qe>ELEVELS-1)
				//	qe=ELEVELS-1;
				int qe=grads[kc+0];
#else
				int qp=(3*(N[kc]+regW[kc])-2*NW[kc])>>(11-CBITS)&((1<<CBITS)-1);
				int qe=abs(N[kc]-regW[kc])+N[kc+4]+regW[kc+4]+((NW[kc+4]+NE[kc+4]+WW[kc+4]+NN[kc+4]+NEE[kc+4])>>1);
				qe=FLOOR_LOG2(qe+1);
				if(qe>ELEVELS-1)
					qe=ELEVELS-1;
#endif
				unsigned *curr_hist0=cdfptrs[kc+0]=histptr+cdfstride*(ELEVELS*(CLEVELS*kc+qp)+qe);
				ctxctrs[kc+0]=&args->ctxctrs[kc][qe][qp];
				int qe2=-(qe<ELEVELS-1);
				unsigned *curr_hist1=cdfptrs[kc+3]=curr_hist0+(cdfstride&qe2);
				ctxctrs[kc+3]=&args->ctxctrs[kc][qe-qe2][qp];

			//	unsigned *curr_hist0=cdfptrs[kc+0]=histptr+cdfstride*(CLEVELS*(ELEVELS*kc+qe)+qp);	//synth2 1.7% larger, GDCC 0.04% smaller
			//	ctxctrs[kc+0]=&args->ctxctrs[kc][qe][qp];
			//	int qp2=qp+(qp<CLEVELS-1);
			//	unsigned *curr_hist1=cdfptrs[kc+3]=histptr+cdfstride*(CLEVELS*(ELEVELS*kc+qe)+qp2);
			//	ctxctrs[kc+3]=&args->ctxctrs[kc][qe][qp2];
#ifdef ENABLE_MIX3
				int qp2=qp+(qp<CLEVELS-1);
				unsigned *curr_hist2=cdfptrs[kc+6]=histptr+cdfstride*(ELEVELS*(CLEVELS*kc+qp2)+qe);
				ctxctrs[kc+6]=&args->ctxctrs[kc][qe][qp2];
	#define GETCDF(X) (((curr_hist0[X]+curr_hist0[X]+curr_hist1[X]+curr_hist2[X])>>(PROB_EBITS+2))+(X))
#elif defined ENABLE_MIX4
				int qp2=qp+(entropyidx>1&&qp<CLEVELS-1);
				unsigned *curr_hist2=cdfptrs[kc+6]=histptr+cdfstride*(ELEVELS*(CLEVELS*kc+qp2)+qe);
				ctxctrs[kc+6]=&args->ctxctrs[kc][qe][qp2];
				unsigned *curr_hist3=cdfptrs[kc+9]=histptr+cdfstride*(ELEVELS*(CLEVELS*kc+qp2)+qe2);
				ctxctrs[kc+9]=&args->ctxctrs[kc][qe2][qp2];
	#define GETCDF(X) (((curr_hist0[X]+curr_hist1[X]+curr_hist2[X]+curr_hist3[X])>>(PROB_EBITS+2))+(X))
#elif defined ENABLE_MIX2
	#define GETCDF(X) (((curr_hist0[X]+curr_hist0[X]+curr_hist0[X]+curr_hist1[X])>>(PROB_EBITS+2))+(X))
#elif defined ENABLE_MIX1
	#define GETCDF(X) ((curr_hist0[X]>>PROB_EBITS)+(X))
#else
	#define GETCDF(X) (((curr_hist0[X]+curr_hist1[X])>>(PROB_EBITS+1))+(X))
#endif
//	#define GETCDF(X) ((curr_hist0[X]>>PROB_EBITS)+(X))

				unsigned cdf, freq;
				if(args->fwd)
				{
					error=yuv[kc]-pred;
					{
						int abserr=abs(error);
						sym=upred+abserr;
						if(abserr<=upred)//CMOV
							sym=error<<1^error>>31;

						//if(abserr<=upred)
						//	sym=error<<1^error>>31;//pack sign
						//else
						//	sym=upred+abserr;//error sign is known
					}
					token=sym;
					if(sym>=(1<<CONFIG_EXP))
					{
						//int nbits=FLOOR_LOG2((unsigned)sym)-(CONFIG_MSB+CONFIG_LSB);
						int nbits=31-(CONFIG_MSB+CONFIG_LSB)-_lzcnt_u32(sym);
						token = (1<<CONFIG_EXP)-((CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB))<<(CONFIG_MSB+CONFIG_LSB)) + (
								nbits<<(CONFIG_MSB+CONFIG_LSB)|
								(sym>>nbits&((1<<CONFIG_MSB)-1)<<CONFIG_LSB)|
								(sym&((1<<CONFIG_LSB)-1))
							);
						bypass=sym>>CONFIG_LSB&((1LL<<nbits)-1);
						//bypass=_bextr_u32(sym>>CONFIG_LSB, 0, nbits);	//slow
						bypass_encbuf(&bc, bypass, nbits);
					}
#ifdef _DEBUG
					if(token>=args->tlevels)
						LOG_ERROR("YXC %d %d %d  token %d/%d", ky, kx, kc, token, args->tlevels);
#endif
					cdf=GETCDF(token);
					freq=GETCDF(token+1)-cdf;
					ac3_encbuf_update_N(&ec, cdf, freq, PROB_BITS);
#ifdef _DEBUG
					if(ec.dstptr>=ec.dstend)
						LOG_ERROR("Encoder out of memory at YXC %d %d %d", ky, kx, kc);
#endif
				}
				else
				{
					if(ec.range<(1ULL<<PROB_BITS))
					{
#ifdef _DEBUG
						if(ec.srcptr>=ec.srcend)
							LOG_ERROR("Decoder out of memory at YXC %d %d %d", ky, kx, kc);
#endif
						ec.code=ec.code<<32|*(unsigned*)ec.srcptr;
						ec.srcptr+=4;
						ec.range=ec.range<<32|0xFFFFFFFF;
						ec.low<<=32;
						if(ec.range>~ec.low)
							ec.range=~ec.low;
					}
#if 1
					if(entropyidx)
					{
						__m256i mr2=_mm256_set1_epi32((unsigned short)(((ec.code-ec.low)<<PROB_BITS|((1LL<<PROB_BITS)-1))/ec.range));
						__m256i mc0=_mm256_loadu_si256((__m256i*)(curr_hist0+1+0*8));
						__m256i mc1=_mm256_loadu_si256((__m256i*)(curr_hist0+1+1*8));
						__m256i mc2=_mm256_loadu_si256((__m256i*)(curr_hist0+1+2*8));
						__m256i mc3=_mm256_loadu_si256((__m256i*)(curr_hist0+1+3*8));
					//	__m256i mc4=_mm256_loadu_si256((__m256i*)(curr_hist0+1+4*8));
					//	__m256i mc5=_mm256_loadu_si256((__m256i*)(curr_hist0+1+5*8));
						__m256i md0=_mm256_loadu_si256((__m256i*)(curr_hist1+1+0*8));
						__m256i md1=_mm256_loadu_si256((__m256i*)(curr_hist1+1+1*8));
						__m256i md2=_mm256_loadu_si256((__m256i*)(curr_hist1+1+2*8));
						__m256i md3=_mm256_loadu_si256((__m256i*)(curr_hist1+1+3*8));
					//	__m256i md4=_mm256_loadu_si256((__m256i*)(curr_hist1+1+4*8));
					//	__m256i md5=_mm256_loadu_si256((__m256i*)(curr_hist1+1+5*8));
#ifdef ENABLE_MIX3
						__m256i me0=_mm256_loadu_si256((__m256i*)(curr_hist2+1+0*8));
						__m256i me1=_mm256_loadu_si256((__m256i*)(curr_hist2+1+1*8));
						__m256i me2=_mm256_loadu_si256((__m256i*)(curr_hist2+1+2*8));
						__m256i me3=_mm256_loadu_si256((__m256i*)(curr_hist2+1+3*8));
					//	__m256i me4=_mm256_loadu_si256((__m256i*)(curr_hist2+1+4*8));
					//	__m256i me5=_mm256_loadu_si256((__m256i*)(curr_hist2+1+5*8));
						__m256i mx0=_mm256_add_epi32(mc0, mc0);
						__m256i mx1=_mm256_add_epi32(mc1, mc1);
						__m256i mx2=_mm256_add_epi32(mc2, mc2);
						__m256i mx3=_mm256_add_epi32(mc3, mc3);
					//	__m256i mx4=_mm256_add_epi32(mc4, mc4);
					//	__m256i mx5=_mm256_add_epi32(mc5, mc5);
						mx0=_mm256_add_epi32(mx0, md0);
						mx1=_mm256_add_epi32(mx1, md1);
						mx2=_mm256_add_epi32(mx2, md2);
						mx3=_mm256_add_epi32(mx3, md3);
					//	mx4=_mm256_add_epi32(mx4, md4);
					//	mx5=_mm256_add_epi32(mx5, md5);
						mx0=_mm256_add_epi32(mx0, me0);
						mx1=_mm256_add_epi32(mx1, me1);
						mx2=_mm256_add_epi32(mx2, me2);
						mx3=_mm256_add_epi32(mx3, me3);
					//	mx4=_mm256_add_epi32(mx4, me4);
					//	mx5=_mm256_add_epi32(mx5, me5);
						mx0=_mm256_srli_epi32(mx0, PROB_EBITS+2);
						mx1=_mm256_srli_epi32(mx1, PROB_EBITS+2);
						mx2=_mm256_srli_epi32(mx2, PROB_EBITS+2);
						mx3=_mm256_srli_epi32(mx3, PROB_EBITS+2);
					//	mx4=_mm256_srli_epi32(mx4, PROB_EBITS+2);
					//	mx5=_mm256_srli_epi32(mx5, PROB_EBITS+2);
#elif defined ENABLE_MIX4
						__m256i me0=_mm256_loadu_si256((__m256i*)(curr_hist2+1+0*8));
						__m256i me1=_mm256_loadu_si256((__m256i*)(curr_hist2+1+1*8));
						__m256i me2=_mm256_loadu_si256((__m256i*)(curr_hist2+1+2*8));
						__m256i me3=_mm256_loadu_si256((__m256i*)(curr_hist2+1+3*8));
					//	__m256i me4=_mm256_loadu_si256((__m256i*)(curr_hist2+1+4*8));
					//	__m256i me5=_mm256_loadu_si256((__m256i*)(curr_hist2+1+5*8));
						__m256i mf0=_mm256_loadu_si256((__m256i*)(curr_hist3+1+0*8));
						__m256i mf1=_mm256_loadu_si256((__m256i*)(curr_hist3+1+1*8));
						__m256i mf2=_mm256_loadu_si256((__m256i*)(curr_hist3+1+2*8));
						__m256i mf3=_mm256_loadu_si256((__m256i*)(curr_hist3+1+3*8));
					//	__m256i mf4=_mm256_loadu_si256((__m256i*)(curr_hist3+1+4*8));
					//	__m256i mf5=_mm256_loadu_si256((__m256i*)(curr_hist3+1+5*8));
						__m256i mx0=_mm256_add_epi32(mc0, md0);
						__m256i mx1=_mm256_add_epi32(mc1, md1);
						__m256i mx2=_mm256_add_epi32(mc2, md2);
						__m256i mx3=_mm256_add_epi32(mc3, md3);
					//	__m256i mx4=_mm256_add_epi32(mc4, md4);
					//	__m256i mx5=_mm256_add_epi32(mc5, md5);
						mx0=_mm256_add_epi32(mx0, me0);
						mx1=_mm256_add_epi32(mx1, me1);
						mx2=_mm256_add_epi32(mx2, me2);
						mx3=_mm256_add_epi32(mx3, me3);
					//	mx4=_mm256_add_epi32(mx4, me4);
					//	mx5=_mm256_add_epi32(mx5, me5);
						mx0=_mm256_add_epi32(mx0, mf0);
						mx1=_mm256_add_epi32(mx1, mf1);
						mx2=_mm256_add_epi32(mx2, mf2);
						mx3=_mm256_add_epi32(mx3, mf3);
					//	mx4=_mm256_add_epi32(mx4, mf4);
					//	mx5=_mm256_add_epi32(mx5, mf5);
						mx0=_mm256_srli_epi32(mx0, PROB_EBITS+2);
						mx1=_mm256_srli_epi32(mx1, PROB_EBITS+2);
						mx2=_mm256_srli_epi32(mx2, PROB_EBITS+2);
						mx3=_mm256_srli_epi32(mx3, PROB_EBITS+2);
					//	mx4=_mm256_srli_epi32(mx4, PROB_EBITS+2);
					//	mx5=_mm256_srli_epi32(mx5, PROB_EBITS+2);
#elif defined ENABLE_MIX2
						__m256i mx0=_mm256_add_epi32(mc0, md0);
						__m256i mx1=_mm256_add_epi32(mc1, md1);
						__m256i mx2=_mm256_add_epi32(mc2, md2);
						__m256i mx3=_mm256_add_epi32(mc3, md3);
					//	__m256i mx4=_mm256_add_epi32(mc4, md4);
					//	__m256i mx5=_mm256_add_epi32(mc5, md5);
						mx0=_mm256_add_epi32(mx0, mc0);
						mx1=_mm256_add_epi32(mx1, mc1);
						mx2=_mm256_add_epi32(mx2, mc2);
						mx3=_mm256_add_epi32(mx3, mc3);
					//	mx4=_mm256_add_epi32(mx4, mc4);
					//	mx5=_mm256_add_epi32(mx5, mc5);
						mx0=_mm256_add_epi32(mx0, mc0);
						mx1=_mm256_add_epi32(mx1, mc1);
						mx2=_mm256_add_epi32(mx2, mc2);
						mx3=_mm256_add_epi32(mx3, mc3);
					//	mx4=_mm256_add_epi32(mx4, mc4);
					//	mx5=_mm256_add_epi32(mx5, mc5);
						mx0=_mm256_srli_epi32(mx0, PROB_EBITS+2);
						mx1=_mm256_srli_epi32(mx1, PROB_EBITS+2);
						mx2=_mm256_srli_epi32(mx2, PROB_EBITS+2);
						mx3=_mm256_srli_epi32(mx3, PROB_EBITS+2);
					//	mx4=_mm256_srli_epi32(mx4, PROB_EBITS+2);
					//	mx5=_mm256_srli_epi32(mx5, PROB_EBITS+2);
#elif defined ENABLE_MIX1
						__m256i mx0=_mm256_srli_epi32(mc0, PROB_EBITS);
						__m256i mx1=_mm256_srli_epi32(mc1, PROB_EBITS);
						__m256i mx2=_mm256_srli_epi32(mc2, PROB_EBITS);
						__m256i mx3=_mm256_srli_epi32(mc3, PROB_EBITS);
					//	__m256i mx4=_mm256_srli_epi32(mc4, PROB_EBITS);
					//	__m256i mx5=_mm256_srli_epi32(mc5, PROB_EBITS);
#else
						__m256i mx0=_mm256_add_epi32(mc0, md0);
						__m256i mx1=_mm256_add_epi32(mc1, md1);
						__m256i mx2=_mm256_add_epi32(mc2, md2);
						__m256i mx3=_mm256_add_epi32(mc3, md3);
					//	__m256i mx4=_mm256_add_epi32(mc4, md4);
					//	__m256i mx5=_mm256_add_epi32(mc5, md5);
						mx0=_mm256_srli_epi32(mx0, PROB_EBITS+1);
						mx1=_mm256_srli_epi32(mx1, PROB_EBITS+1);
						mx2=_mm256_srli_epi32(mx2, PROB_EBITS+1);
						mx3=_mm256_srli_epi32(mx3, PROB_EBITS+1);
					//	mx4=_mm256_srli_epi32(mx4, PROB_EBITS+1);
					//	mx5=_mm256_srli_epi32(mx5, PROB_EBITS+1);
#endif
						mx0=_mm256_add_epi32(mx0, _mm256_set_epi32( 8,  7,  6,  5,  4,  3,  2,  1));
						mx1=_mm256_add_epi32(mx1, _mm256_set_epi32(16, 15, 14, 13, 12, 11, 10,  9));
						mx2=_mm256_add_epi32(mx2, _mm256_set_epi32(24, 23, 22, 21, 20, 19, 18, 17));
						mx3=_mm256_add_epi32(mx3, _mm256_set_epi32(32, 31, 30, 29, 28, 27, 26, 25));
					//	mx4=_mm256_add_epi32(mx4, _mm256_set_epi32(40, 39, 38, 37, 36, 35, 34, 33));
					//	mx5=_mm256_add_epi32(mx5, _mm256_set_epi32(48, 47, 46, 45, 44, 43, 42, 41));
						_mm256_store_si256((__m256i*)CDF+1, mx0);
						_mm256_store_si256((__m256i*)CDF+2, mx1);
						_mm256_store_si256((__m256i*)CDF+3, mx2);
						_mm256_store_si256((__m256i*)CDF+4, mx3);
					//	_mm256_store_si256((__m256i*)CDF+5, mx4);
					//	_mm256_store_si256((__m256i*)CDF+6, mx5);
						__m256i cond0=_mm256_cmpgt_epi32(mx0, mr2);
						__m256i cond1=_mm256_cmpgt_epi32(mx1, mr2);
						__m256i cond2=_mm256_cmpgt_epi32(mx2, mr2);
						__m256i cond3=_mm256_cmpgt_epi32(mx3, mr2);
					//	__m256i cond4=_mm256_cmpgt_epi32(mx4, mr2);
					//	__m256i cond5=_mm256_cmpgt_epi32(mx5, mr2);
						unsigned long long mask=_mm256_movemask_ps(_mm256_castsi256_ps(cond0));
						mask|=(unsigned long long)_mm256_movemask_ps(_mm256_castsi256_ps(cond1))<<8;
						mask|=(unsigned long long)_mm256_movemask_ps(_mm256_castsi256_ps(cond2))<<16;
						mask|=(unsigned long long)_mm256_movemask_ps(_mm256_castsi256_ps(cond3))<<24;
					//	mask|=(unsigned long long)_mm256_movemask_ps(_mm256_castsi256_ps(cond4))<<32;
					//	mask|=(unsigned long long)_mm256_movemask_ps(_mm256_castsi256_ps(cond5))<<40;
						mask|=1ULL<<32;
						token=(int)_tzcnt_u64(mask);//BMI1 (2013, with LZCNT & AVX2)
						
						cdf=CDF[token+7];
						freq=CDF[token+8]-cdf;
					}
					else//low entropy
#endif
					{
						unsigned long long r2=(ec.code-ec.low)<<PROB_BITS|((1LL<<PROB_BITS)-1);
						cdf=0;
						token=0;
						for(;;)
						{
							freq=GETCDF(token+1);
							if(freq*ec.range>r2)
								break;
#ifdef _DEBUG
							if(token>=args->tlevels)
								LOG_ERROR("YXC %d %d %d  token %d/%d", ky, kx, kc, token, args->tlevels);
#endif
							++token;
							cdf=freq;
						}
						freq-=cdf;
					}
					ac3_dec_update_N(&ec, cdf, freq, PROB_BITS);
					sym=token;
					if(token>=(1<<CONFIG_EXP))
					{
						nbits=(sym>>(CONFIG_MSB+CONFIG_LSB))-((1<<(CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB)))-(CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB)));//2 instructions
						bypass=bypass_decbuf(&bc, nbits);
						sym=
							+(((1<<(CONFIG_MSB+CONFIG_LSB))+(sym&((1<<CONFIG_MSB)-1)<<CONFIG_LSB))<<nbits)//7 instructions
							+(bypass<<CONFIG_LSB)
							+(sym&((1<<CONFIG_LSB)-1))
						;
					}
					error=upred-sym;
					if(pred<0)
						error=sym-upred;
					if(sym<=(upred<<1))//CMOV
						error=sym>>1^-(sym&1);

					//if(sym<=(upred<<1))
					//	error=sym>>1^-(sym&1);//unpack sign
					//else
					//{
					//	int neg=pred>>31;
					//	error=upred-sym;
					//	error^=neg;
					//	error-=neg;
					//
					//	//error=sym-upred;
					//	//if(pred>0)
					//	//	error=-error;
					//}
					yuv[kc]=error+pred;
				}
				curr[kc+0]=regW[kc]=yuv[kc]-offset;
				curr[kc+4]=regW[kc+4]=abs(error);
			//	regW[kc+8]=error;
				mixinptrs[kc]=args->mixincdfs+args->tlevels*token;
			}
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 3
#endif
			for(int kc=0;kc<3;++kc)
			{
				int sh=*ctxctrs[kc+0]+=shiftgain;
				sh=FLOOR_LOG2(sh)>>1;
#if !defined ENABLE_MT && 0
				{
					static int shist[SMAX+1-SMIN]={0};
					++shist[sh-SMIN];
					if(args->fwd&&ky==args->y2-1&&kx==args->x2-1&&kc==2)
					{
						printf("block %d\n", args->blockidx);
						for(int k=0;k<SMAX+1-SMIN;++k)
							printf("%2d %8d\n", k+SMIN, shist[k]);
					}
				}
#endif
				__m128i ms=_mm_set_epi32(0, 0, 0, sh);
				__m256i my0=_mm256_loadu_si256((__m256i*)(mixinptrs[kc]+1+0*8));
				__m256i my1=_mm256_loadu_si256((__m256i*)(mixinptrs[kc]+1+1*8));
				__m256i my2=_mm256_loadu_si256((__m256i*)(mixinptrs[kc]+1+2*8));
				__m256i my3=_mm256_loadu_si256((__m256i*)(mixinptrs[kc]+1+3*8));
			//	__m256i my4=_mm256_loadu_si256((__m256i*)(mixinptrs[kc]+1+4*8));
			//	__m256i my5=_mm256_loadu_si256((__m256i*)(mixinptrs[kc]+1+5*8));
				__m256i mc0=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+0]+1+0*8));
				__m256i mc1=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+0]+1+1*8));
				__m256i mc2=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+0]+1+2*8));
				__m256i mc3=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+0]+1+3*8));
			//	__m256i mc4=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+0]+1+4*8));
			//	__m256i mc5=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+0]+1+5*8));
				__m256i mx0=_mm256_sub_epi32(my0, mc0);
				__m256i mx1=_mm256_sub_epi32(my1, mc1);
				__m256i mx2=_mm256_sub_epi32(my2, mc2);
				__m256i mx3=_mm256_sub_epi32(my3, mc3);
			//	__m256i mx4=_mm256_sub_epi32(my4, mc4);
			//	__m256i mx5=_mm256_sub_epi32(my5, mc5);
				mx0=_mm256_sra_epi32(mx0, ms);
				mx1=_mm256_sra_epi32(mx1, ms);
				mx2=_mm256_sra_epi32(mx2, ms);
				mx3=_mm256_sra_epi32(mx3, ms);
			//	mx4=_mm256_sra_epi32(mx4, ms);
			//	mx5=_mm256_sra_epi32(mx5, ms);
				mx0=_mm256_add_epi32(mx0, mc0);
				mx1=_mm256_add_epi32(mx1, mc1);
				mx2=_mm256_add_epi32(mx2, mc2);
				mx3=_mm256_add_epi32(mx3, mc3);
			//	mx4=_mm256_add_epi32(mx4, mc4);
			//	mx5=_mm256_add_epi32(mx5, mc5);
				_mm256_storeu_si256((__m256i*)(cdfptrs[kc+0]+1+0*8), mx0);
				_mm256_storeu_si256((__m256i*)(cdfptrs[kc+0]+1+1*8), mx1);
				_mm256_storeu_si256((__m256i*)(cdfptrs[kc+0]+1+2*8), mx2);
				_mm256_storeu_si256((__m256i*)(cdfptrs[kc+0]+1+3*8), mx3);
			//	_mm256_storeu_si256((__m256i*)(cdfptrs[kc+0]+1+4*8), mx4);
			//	_mm256_storeu_si256((__m256i*)(cdfptrs[kc+0]+1+5*8), mx5);
				
#ifndef ENABLE_MIX1
				//if(cdfptrs[kc+3]!=cdfptrs[kc+0])
				{
					sh=*ctxctrs[kc+3]+=shiftgain;
					sh=FLOOR_LOG2(sh)>>1;
					ms=_mm_set_epi32(0, 0, 0, sh);
					mc0=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+3]+1+0*8));
					mc1=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+3]+1+1*8));
					mc2=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+3]+1+2*8));
					mc3=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+3]+1+3*8));
				//	mc4=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+3]+1+4*8));
				//	mc5=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+3]+1+5*8));
					mx0=_mm256_sub_epi32(my0, mc0);
					mx1=_mm256_sub_epi32(my1, mc1);
					mx2=_mm256_sub_epi32(my2, mc2);
					mx3=_mm256_sub_epi32(my3, mc3);
				//	mx4=_mm256_sub_epi32(my4, mc4);
				//	mx5=_mm256_sub_epi32(my5, mc5);
					mx0=_mm256_sra_epi32(mx0, ms);
					mx1=_mm256_sra_epi32(mx1, ms);
					mx2=_mm256_sra_epi32(mx2, ms);
					mx3=_mm256_sra_epi32(mx3, ms);
				//	mx4=_mm256_sra_epi32(mx4, ms);
				//	mx5=_mm256_sra_epi32(mx5, ms);
					mx0=_mm256_add_epi32(mx0, mc0);
					mx1=_mm256_add_epi32(mx1, mc1);
					mx2=_mm256_add_epi32(mx2, mc2);
					mx3=_mm256_add_epi32(mx3, mc3);
				//	mx4=_mm256_add_epi32(mx4, mc4);
				//	mx5=_mm256_add_epi32(mx5, mc5);
					_mm256_storeu_si256((__m256i*)(cdfptrs[kc+3]+1+0*8), mx0);
					_mm256_storeu_si256((__m256i*)(cdfptrs[kc+3]+1+1*8), mx1);
					_mm256_storeu_si256((__m256i*)(cdfptrs[kc+3]+1+2*8), mx2);
					_mm256_storeu_si256((__m256i*)(cdfptrs[kc+3]+1+3*8), mx3);
				//	_mm256_storeu_si256((__m256i*)(cdfptrs[kc+3]+1+4*8), mx4);
				//	_mm256_storeu_si256((__m256i*)(cdfptrs[kc+3]+1+5*8), mx5);
				}
#endif
#if defined ENABLE_MIX3 || defined ENABLE_MIX4
				//if(cdfptrs[kc+6]!=cdfptrs[kc+0])
				{
					sh=*ctxctrs[kc+6]+=shiftgain;
					sh=FLOOR_LOG2(sh)>>1;
					ms=_mm_set_epi32(0, 0, 0, sh);
					mc0=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+6]+1+0*8));
					mc1=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+6]+1+1*8));
					mc2=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+6]+1+2*8));
					mc3=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+6]+1+3*8));
				//	mc4=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+6]+1+4*8));
				//	mc5=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+6]+1+5*8));
					mx0=_mm256_sub_epi32(my0, mc0);
					mx1=_mm256_sub_epi32(my1, mc1);
					mx2=_mm256_sub_epi32(my2, mc2);
					mx3=_mm256_sub_epi32(my3, mc3);
				//	mx4=_mm256_sub_epi32(my4, mc4);
				//	mx5=_mm256_sub_epi32(my5, mc5);
					mx0=_mm256_sra_epi32(mx0, ms);
					mx1=_mm256_sra_epi32(mx1, ms);
					mx2=_mm256_sra_epi32(mx2, ms);
					mx3=_mm256_sra_epi32(mx3, ms);
				//	mx4=_mm256_sra_epi32(mx4, ms);
				//	mx5=_mm256_sra_epi32(mx5, ms);
					mx0=_mm256_add_epi32(mx0, mc0);
					mx1=_mm256_add_epi32(mx1, mc1);
					mx2=_mm256_add_epi32(mx2, mc2);
					mx3=_mm256_add_epi32(mx3, mc3);
				//	mx4=_mm256_add_epi32(mx4, mc4);
				//	mx5=_mm256_add_epi32(mx5, mc5);
					_mm256_storeu_si256((__m256i*)(cdfptrs[kc+6]+1+0*8), mx0);
					_mm256_storeu_si256((__m256i*)(cdfptrs[kc+6]+1+1*8), mx1);
					_mm256_storeu_si256((__m256i*)(cdfptrs[kc+6]+1+2*8), mx2);
					_mm256_storeu_si256((__m256i*)(cdfptrs[kc+6]+1+3*8), mx3);
				//	_mm256_storeu_si256((__m256i*)(cdfptrs[kc+6]+1+4*8), mx4);
				//	_mm256_storeu_si256((__m256i*)(cdfptrs[kc+6]+1+5*8), mx5);
				}
#endif
#ifdef ENABLE_MIX4
				//if(cdfptrs[kc+9]!=cdfptrs[kc+0])
				{
					sh=*ctxctrs[kc+9]+=shiftgain;
					sh=FLOOR_LOG2(sh)>>1;
					ms=_mm_set_epi32(0, 0, 0, sh);
					mc0=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+9]+1+0*8));
					mc1=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+9]+1+1*8));
					mc2=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+9]+1+2*8));
					mc3=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+9]+1+3*8));
				//	mc4=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+9]+1+4*8));
				//	mc5=_mm256_loadu_si256((__m256i*)(cdfptrs[kc+9]+1+5*8));
					mx0=_mm256_sub_epi32(my0, mc0);
					mx1=_mm256_sub_epi32(my1, mc1);
					mx2=_mm256_sub_epi32(my2, mc2);
					mx3=_mm256_sub_epi32(my3, mc3);
				//	mx4=_mm256_sub_epi32(my4, mc4);
				//	mx5=_mm256_sub_epi32(my5, mc5);
					mx0=_mm256_sra_epi32(mx0, ms);
					mx1=_mm256_sra_epi32(mx1, ms);
					mx2=_mm256_sra_epi32(mx2, ms);
					mx3=_mm256_sra_epi32(mx3, ms);
				//	mx4=_mm256_sra_epi32(mx4, ms);
				//	mx5=_mm256_sra_epi32(mx5, ms);
					mx0=_mm256_add_epi32(mx0, mc0);
					mx1=_mm256_add_epi32(mx1, mc1);
					mx2=_mm256_add_epi32(mx2, mc2);
					mx3=_mm256_add_epi32(mx3, mc3);
				//	mx4=_mm256_add_epi32(mx4, mc4);
				//	mx5=_mm256_add_epi32(mx5, mc5);
					_mm256_storeu_si256((__m256i*)(cdfptrs[kc+9]+1+0*8), mx0);
					_mm256_storeu_si256((__m256i*)(cdfptrs[kc+9]+1+1*8), mx1);
					_mm256_storeu_si256((__m256i*)(cdfptrs[kc+9]+1+2*8), mx2);
					_mm256_storeu_si256((__m256i*)(cdfptrs[kc+9]+1+3*8), mx3);
				//	_mm256_storeu_si256((__m256i*)(cdfptrs[kc+9]+1+4*8), mx4);
				//	_mm256_storeu_si256((__m256i*)(cdfptrs[kc+9]+1+5*8), mx5);
				}
#endif
			}
			if(!args->fwd)
			{
				args->dst[yidx]=yuv[0]+128;
				args->dst[uidx]=yuv[1]+128;
				args->dst[vidx]=yuv[2]+128;
				guide_check(args->dst, kx, ky);
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
#if defined SIMD_CTX || defined SIMD_PREDS
			mNNW	=mNN;
			mNN	=mNNE;
			mNNE	=_mm_loadu_si128((__m128i*)(rows[2]+1*4*2));
			mNWW	=mNW;
			mNW	=mN;
			mN	=mNE;
			mNE	=mNEE;
			mNEE	=_mm_loadu_si128((__m128i*)(rows[1]+2*4*2));
			mWW	=mW;
			mW	=_mm_load_si128((__m128i*)regW);
		//	mW	=_mm_loadu_si128((__m128i*)(rows[0]-1*4*2));
#endif
		}
	}
	if(args->fwd)
	{
		int currblockidx=args->blockidx-args->b1;
		unsigned char *ptr;

		ptr=ac3_encbuf_flush(&ec);
		args->enctokenoffsets[currblockidx]=(int)(ptr-args->enctokenbuf);
#ifdef _DEBUG
		if((ptr-args->enctokenbuf)>>32)
			LOG_ERROR("Integer overflow");
#endif
		
		ptr=bypass_encbuf_flush(&bc);
		args->encbypassoffsets[currblockidx]=(int)(ptr-args->encbypassbuf);
#ifdef _DEBUG
		if((ptr-args->encbypassbuf)>>32)
			LOG_ERROR("Integer overflow");
#endif
	}
}
static void block_thread(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	for(int kb=args->b1;kb<args->b2;++kb)
	{
		int kx, ky;

		args->blockidx=kb;
		ky=args->blockidx/args->xblocks;
		kx=args->blockidx%args->xblocks;
		args->x1=BLOCKX*kx;
		args->y1=BLOCKY*ky;
		args->x2=MINVAR(args->x1+BLOCKX, args->iw);
		args->y2=MINVAR(args->y1+BLOCKY, args->ih);
		//if(!args->fwd)
		//{
		//	args->decstart=args->decbuf+args->offsets[kb*2+0];
		//	args->decend=args->decbuf+args->offsets[kb+2+1];
		//	args->bypassstart=args->decbuf+args->offsets[kb*2+1];
		//	args->bypassend=args->decbuf+args->offsets[kb*2+2];
		//}
		block_func(param);
	}
}
int c24_codec(const char *srcfn, const char *dstfn)
{
	double t0;
	ArrayHandle src, dst;
	int headersize, printed;
	int iw, ih;
	const unsigned char *image, *imageend;
	unsigned char *image2;
	CodecID codec;
	int xblocks, yblocks, nblocks, ncores, nthreads, coffset;
	ptrdiff_t start, memusage, argssize;
	ThreadArgs *args;
	int test, fwd;
	int histsize;
	double esize;
	int usize;
	int maxwidth;
	int tlevels, statssize;
	unsigned *mixincdfs=0;
#ifdef PRINT_PREDHIST
	int predhist[PRED_COUNT]={0};
#endif
	int *streamoffsets=0;
	int encbufsize=0;
	
	t0=time_sec();
	src=load_file(srcfn, 1, 3, 1);
	headersize=header_read(src->data, (int)src->count, &iw, &ih, &codec);
	image=src->data+headersize;
	imageend=src->data+src->count;
	if(codec==CODEC_INVALID||codec==CODEC_PGM)
	{
		LOG_ERROR("Unsupported codec %d.\n", codec);
		array_free(&src);
		return 1;
	}
	else if(codec==CODEC_C01&&!dstfn)
	{
		LOG_ERROR(
			"Test mode expects PPM source.\n"
			"Decode mode expects destination filename."
		);
		return 1;
	}
	test=!dstfn;
	fwd=codec==CODEC_PPM;
	
	usize=iw*ih*3;
	xblocks=(iw+BLOCKX-1)/BLOCKX;
	yblocks=(ih+BLOCKY-1)/BLOCKY;
	nblocks=xblocks*yblocks;
	ncores=query_cpu_cores();
	nthreads=MINVAR(nblocks, ncores);
	coffset=(int)sizeof(int[2])*nblocks;
	start=0;
	memusage=0;
	argssize=nthreads*sizeof(ThreadArgs);
	args=(ThreadArgs*)malloc(argssize);
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	esize=0;
	memusage+=argssize;
	memset(args, 0, argssize);
	streamoffsets=(int*)malloc(sizeof(int[2])*nblocks);
	if(!streamoffsets)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(streamoffsets, 0, sizeof(int[2])*nblocks);
	if(fwd)
	{
		guide_save(image, iw, ih);

		dst=0;
		printed=snprintf(g_buf, G_BUF_SIZE-1, "C01\n%d %d\n", iw, ih);
		array_append(&dst, g_buf, 1, printed, 1, 0, 0);
		start=dst->count;
		//start=array_append(&dst, 0, 1, coffset, 1, 0, 0);
		
		image2=0;
		encbufsize=(int)((long long)4*iw*ih/nthreads);
	}
	else//integrity check
	{
		dst=0;
		printed=snprintf(g_buf, G_BUF_SIZE-1, "P6\n%d %d\n255\n", iw, ih);
		array_append(&dst, g_buf, 1, printed, 1, 0, 0);
		array_append(&dst, 0, 1, usize, 1, 0, 0);

		memcpy(streamoffsets, image, sizeof(int[2])*nblocks);
		int sum=0;
		for(int k=0;k<nblocks*2;++k)
		{
			int size=streamoffsets[k];
			sum+=size;
			streamoffsets[k]=sum;
		}
		start=coffset;
		if(image+coffset+sum!=imageend)
			LOG_ERROR("Corrupt file");
		//start=coffset;
		//for(int kt=0;kt<nblocks;++kt)
		//{
		//	int size[2]={0};
		//	memcpy(size, image+sizeof(int[2])*kt, sizeof(int[2]));
		//	start+=size[0];
		//	start+=size[1];
		//}
		//if(image+start!=imageend)
		//	LOG_ERROR("Corrupt file");
		//start=coffset;

		image2=(unsigned char*)malloc(usize);
		if(!image2)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		memset(image2, 0, usize);
	}
	{
		int token=0, bypass=0, nbits=0;

		quantize_pixel(256, &token, &bypass, &nbits);//256 is a valid symbol due to CALIC-like sign packing
		tlevels=token+1;
		statssize=(tlevels+1)*(int)sizeof(int[3*ELEVELS*CLEVELS]);//CDF padding, contains (1<<PROB_BITS)
		histsize=(int)sizeof(int[OCH_COUNT*PRED_COUNT<<8]);
		if(histsize<statssize)
			histsize=statssize;
		maxwidth=iw;
		if(maxwidth>BLOCKX)
			maxwidth=BLOCKX;
		mixincdfs=(unsigned*)malloc(sizeof(int)*tlevels*tlevels);
		if(!mixincdfs)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		for(int ks=0;ks<tlevels;++ks)
		{
			for(int ks2=0;ks2<tlevels;++ks2)
				mixincdfs[tlevels*ks+ks2]=(((1<<PROB_BITS)-tlevels)<<PROB_EBITS)&-(ks<ks2);
		}
	}
	for(int k=0, kb=0;k<nthreads;++k)//initialization
	{
		ThreadArgs *arg=args+k;
		arg->src=image;
		arg->dst=fwd?0:dst->data+printed;
		arg->iw=iw;
		arg->ih=ih;
		
		arg->b1=kb;
		kb=(k+1)*nblocks/nthreads;
		arg->b2=kb;
		arg->xblocks=xblocks;
		arg->nblocks=nblocks;

		arg->bufsize=sizeof(short[4*OCH_COUNT*2])*(maxwidth+16LL);//4 padded rows * OCH_COUNT * {pixels, wg_errors}
		arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));
		arg->histsize=histsize;
		arg->hist=(int*)malloc(histsize);
		if(!arg->pixels||!arg->hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=arg->bufsize;
		memusage+=arg->histsize;
		if(fwd)
		{
			arg->encbufsize=encbufsize;
			arg->enctokenbuf=(unsigned char*)malloc(encbufsize);
			arg->encbypassbuf=(unsigned char*)malloc(encbufsize);
			if(!arg->enctokenbuf||!arg->encbypassbuf)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			memusage+=2LL*encbufsize;

			arg->enctokenoffsets=streamoffsets+arg->b1;
			arg->encbypassoffsets=arg->enctokenoffsets+nblocks;
		}
		else
		{
			arg->decstream=image+start;
			arg->decoffsets=streamoffsets;
		}
		//if(fwd)
		//{
		//	arg->tokenbufsize=(ptrdiff_t)4*BLOCKX*BLOCKY;//actually 2.5*BLOCKX*BLOCKY
		//	arg->tokenbuf=(unsigned char*)malloc(arg->tokenbufsize);
		//	arg->bypassbufsize=(ptrdiff_t)4*BLOCKX*BLOCKY;//actually 2.7*BLOCKX*BLOCKY
		//	arg->bypassbuf=(unsigned char*)malloc(arg->bypassbufsize+4);//4 byte pad for branchless renorm
		//	if(!arg->tokenbuf||!arg->bypassbuf)
		//	{
		//		LOG_ERROR("Alloc error");
		//		return 1;
		//	}
		//	memusage+=arg->tokenbufsize;
		//	memusage+=arg->bypassbufsize;
		//}
		
		arg->tlevels=tlevels;
		arg->fwd=fwd;
		arg->test=test;
		arg->mixincdfs=mixincdfs;
#ifdef ENABLE_MT
		arg->loud=0;
#else
		arg->loud=test&&nblocks<MAXPRINTEDBLOCKS;
#endif
	}
	for(int k2=0;k2<=test;++k2)
	{
#ifdef ENABLE_MT
		void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads);
		mt_finish(ctx);
#else
		for(int k=0;k<nthreads;++k)
			block_thread(args+k);
#endif
		if(fwd)
		{
			{
				int *streamsizes=(int*)ARRAY_APPEND(dst, 0, sizeof(int[2])*nblocks, 1, 0);
				if(!streamsizes)
				{
					LOG_ERROR("Alloc error");
					return 1;
				}
				for(int k=0;k<nthreads;++k)//append token stream sizes
				{
					ThreadArgs *arg=args+k;
					for(int kb=arg->b1;kb<arg->b2;++kb)
						streamsizes[kb]=streamoffsets[kb]-(kb>arg->b1?streamoffsets[kb-1]:0);
				}
				for(int k=0;k<nthreads;++k)//append bypass stream sizes
				{
					ThreadArgs *arg=args+k;
					for(int kb=arg->b1;kb<arg->b2;++kb)
					{
						int kb2=nblocks+kb;
						streamsizes[kb2]=streamoffsets[kb2]-(kb>arg->b1?streamoffsets[kb2-1]:0);
					}
				}
				//for(int k=0;k<nblocks*2;++k)
				//	printf("%3d  %8d\n", k, streamsizes[k]);
			}

			for(int k=0;k<nthreads;++k)//append token streams
			{
				ThreadArgs *arg=args+k;
				ARRAY_APPEND(dst, arg->enctokenbuf, streamoffsets[arg->b2-1], 1, 0);
			}
			for(int k=0;k<nthreads;++k)//append bypass streams
			{
				ThreadArgs *arg=args+k;
				ARRAY_APPEND(dst, arg->encbypassbuf, streamoffsets[nblocks+arg->b2-1], 1, 0);
			}
		}
#if 0
		for(int kt=0;kt<nblocks;kt+=nthreads)
		{
			int nthreads2=MINVAR(kt+nthreads, nblocks)-kt;
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				int kx, ky;
			
				arg->blockidx=kx=kt+kt2;
				ky=kx/xblocks;
				kx%=xblocks;
				arg->x1=BLOCKX*kx;
				arg->y1=BLOCKY*ky;
				arg->x2=MINVAR(arg->x1+BLOCKX, iw);
				arg->y2=MINVAR(arg->y1+BLOCKY, ih);
				//if(!fwd)
				//{
				//	int size[2]={0};
				//	memcpy(size, image+sizeof(int[2])*((ptrdiff_t)kt+kt2), sizeof(int[2]));
				//
				//	arg->decstart=image+start;
				//	start+=size[0];
				//	arg->decend=image+start;
				//
				//	arg->bypassstart=image+start;
				//	start+=size[1];
				//	arg->bypassend=image+start;
				//}
			}
#ifdef ENABLE_MT
			void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
			mt_finish(ctx);
#else
			for(int k=0;k<nthreads2;++k)
				block_thread(args+k);
#endif
			if(fwd)
			{
				for(int kt2=0;kt2<nthreads2;++kt2)
				{
					ThreadArgs *arg=args+kt2;
					if(test)
					{
#if 0
						int blocksize=(3*8*(arg->x2-arg->x1)*(arg->y2-arg->y1)+7)>>3;
						int kx, ky;

						kx=kt+kt2;
						ky=kx/xblocks;
						kx%=xblocks;
						if(nblocks<MAXPRINTEDBLOCKS)
						{
							ptrdiff_t nbytes=arg->tokenstream.nbytes+args->bypassstream.nbytes;
							//if(!(kt+kt2))
							//	printf("block,  nrows,  usize,     best  ->  actual,  (actual-best)\n");
							printf(
								"block %4d/%4d  XY %3d %3d  %4d*%4d:  %8d->%16lf->%8zd bytes (%+10.2lf)  %10.6lf%%  CR %10lf  %s %s %s %s\n",
								kt+kt2+1, nblocks,
								kx, ky,
								arg->y2-arg->y1,
								arg->x2-arg->x1,
								blocksize,
								arg->bestsize,
								nbytes,
								nbytes-arg->bestsize,
								100.*nbytes/blocksize,
								(double)blocksize/nbytes,
								rct_names[arg->bestrct],
								pred_names[arg->predidx[0]],
								pred_names[arg->predidx[1]],
								pred_names[arg->predidx[2]]
							);
						}
#endif
						esize+=arg->bestsize;
					}
					//ARRAY_APPEND(dst, arg->tokenbuf, arg->tokenptr-arg->tokenbuf, 1, 0);
					//ARRAY_APPEND(dst, arg->bypassbuf, arg->bypassptr-arg->bypassbuf, 1, 0);
					//*(unsigned*)(dst->data+start+8LL*((ptrdiff_t)kt+kt2)+4*0)=(unsigned)(arg->tokenptr-arg->tokenbuf);
					//*(unsigned*)(dst->data+start+8LL*((ptrdiff_t)kt+kt2)+4*1)=(unsigned)(arg->bypassptr-arg->bypassbuf);
#ifdef PRINT_PREDHIST
					++predhist[arg->predidx[0]];
					++predhist[arg->predidx[1]];
					++predhist[arg->predidx[2]];
#endif
				}
			}
		}
#endif
		if(test)
		{
			ptrdiff_t usize=((ptrdiff_t)3*8*iw*ih+7)>>3;
			ptrdiff_t csize=dst->count;
			t0=time_sec()-t0;
			if(fwd)
			{
				printf("Best %15.2lf (%+13.2lf) bytes\n", esize, csize-esize);
				printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
				printf("Mem usage: ");
				print_size((double)memusage, 8, 4, 0, 0);
				printf("\n");
			}
			printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
			if(!fwd)
				compare_bufs_8(image2, src->data+headersize, iw, ih, 3, 3, "C01", 0, 1);
		}
		if(!k2&&test)
		{
			int usize=iw*ih*3;
			fwd=0;
			image2=(unsigned char*)malloc(usize);
			if(!image2)
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			memset(image2, 0, usize);
			for(int kt=0;kt<nthreads;++kt)
			{
				ThreadArgs *arg=args+kt;
				arg->dst=image2;
				arg->fwd=0;
			}
			image=dst->data+printed;
			start=coffset;
		}
		t0=time_sec();
	}
	if(!test)
		save_file(dstfn, dst->data, dst->count, 1);
	if(image2)
		free(image2);
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		_mm_free(arg->pixels);
		free(arg->hist);
		if(fwd||test)
		{
			free(arg->enctokenbuf);
			free(arg->encbypassbuf);
		}
	}
	free(args);
	array_free(&src);
	array_free(&dst);
	free(mixincdfs);
	free(streamoffsets);
#ifdef PRINT_PREDHIST
	if(fwd)
	{
		for(int k=0;k<PRED_COUNT;++k)
			printf("\t%8d %s", predhist[k], pred_names[k]);
		printf("\n");
	}
#endif
	return 0;
}