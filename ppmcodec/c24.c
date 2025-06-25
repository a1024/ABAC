static const char file[]=__FILE__;
#include"util.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<sys/stat.h>


#ifdef _MSC_VER
	#define LOUD
	#define ENABLE_GUIDE
//	#define ENABLE_VALIDATION
#endif
#if !defined DISABLE_MT && !defined ENABLE_VALIDATION
	#define ENABLE_MT
#endif

//	#define AC_DISABLE_WRITE
//	#define PRINT_PREDHIST
//	#define PRINT_CSIZES

//OG preds: W/NE/CG/AV4/AV9/WG

	#define ENABLE_ZERO	//good with astro, requires static-o0 to be fast
//	#define ENABLE_CLEARTYPE//slightly compensates ZERO-stealing with text
//	#define ENABLE_N	//bad with GDCC
	#define ENABLE_W	//good with GDCC
	#define ENABLE_NW	//good
//	#define ENABLE_NE	//?
//	#define ENABLE_NWE	//bad
//	#define ENABLE_CG45	//bad
	#define ENABLE_AV2	//good, but redundant with WG
//	#define ENABLE_GRAD	//weaker
//	#define ENABLE_CGX	//bad
//	#define ENABLE_CGY	//bad
	#define ENABLE_CG
//	#define ENABLE_CGv2	//useless
	#define ENABLE_SELECT	//good with synth
	#define ENABLE_IZ	//weak with CG/AV4
	#define ENABLE_AV3	//weak with CG/AV4
	#define ENABLE_AV4	//good with GDCC
//	#define ENABLE_AV4v2	//?
//	#define ENABLE_AV5	//bad, redundant with AV4
	#define ENABLE_AV6	//good
	#define ENABLE_AV8	//?		LPF
	#define ENABLE_AV9	//good				GDCC 0.31% slower, 0.37% smaller
//	#define ENABLE_AV9v2	//weak				GDCC 3~5% slower, 0.03% smaller
//	#define ENABLE_AV11	//bad
	#define ENABLE_WG	//for noisy areas		GDCC 5.69% slower, 0.07% smaller
//	#define ENABLE_WGB	//bad, steals from CG?
//	#define ENABLE_WG2	//bad
//	#define ENABLE_WG3	//bad
//	#define ENABLE_WG4	//?		WIP

	#define SIMD_CTX	//0.8% faster
//	#define SIMD_PREDS	//3.8% slower

//	#define STATIC_ESTIMATE
//	#define STATIC_ESTIMATE2
#define PBITS 1
#define PLEVELS (1<<PBITS)

	#define STATIC_ZERO
//	#define ENABLE_MIX1	//inefficient			GDCC 9.7~6% faster, 0.6% larger
//	#define ENABLE_MIX2	//bad
//	#define ENABLE_MIX3	//slow				GDCC 16% slower, 0.2% smaller
	#define ENABLE_MIX4	//worse with synth		GDCC 28% slower, 0.3% smaller

//	#define AC3_ENC_BRANCHLESSRENORM	//slow
//	#define AC3_DEC_BRANCHLESSRENORM	//slow	//insignificantly (0.2%) faster?
//	#define BYPASS_ENC_BRANCHLESSRENORM	//slow	//insignificantly (0.4%) faster?
//	#define BYPASS_DEC_BRANCHLESSRENORM	//slow

//PROB_BITS <= 16
//PROB_BITS + PROB_EBITS + 1 <= 31
#define PROB_BITS 16
#define PROB_EBITS 14

#define ANALYSIS_XSTRIDE 2
#define ANALYSIS_YSTRIDE 3

#ifdef ENABLE_VALIDATION
	#define AC_VALIDATE
	#define AC_IMPLEMENTATION
#endif
#include"entropy.h"

#define BLOCKX 448	//384
#define BLOCKY 448

#define MAXPRINTEDBLOCKS 500
//11<<6 = 704 CDFs size 32
#define ELEVELS 11
#define CBITS_EFFECTIVE 4
#define CBITS 6
#define CLEVELS (1<<CBITS)
#define WG4_PREC 5
#define WG4_LR 5

#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static double g_sqe=0;
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
static void guide_update_sqe(const unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	double diff;
	diff=g_image[idx+0]-image[idx+0]; g_sqe+=diff*diff;
	diff=g_image[idx+1]-image[idx+1]; g_sqe+=diff*diff;
	diff=g_image[idx+2]-image[idx+2]; g_sqe+=diff*diff;
}
#else
#define guide_save(...)
#define guide_check(...)
#define guide_update_sqe(...)
#endif
#define OCHLIST\
	OCH(R)\
	OCH(G)\
	OCH(B)\
	OCH(RG)\
	OCH(GB)\
	OCH(BR)\
	OCH(R2)\
	OCH(G2)\
	OCH(B2)\
	OCH(RB)\
	OCH(GR)\
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

typedef enum _RCTInfoIdx
{
	II_OCH_Y,
	II_OCH_U,
	II_OCH_V,

	II_PERM_Y,
	II_PERM_U,
	II_PERM_V,

	II_HELP_U,
	II_HELP_V0,
	II_HELP_V1,

	II_COUNT,
} RCTInfoIdx;
#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 2)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  2, 0)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  2, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 2)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 2)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  2, 0)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	2,  2, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	2,  0, 2)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	2,  0, 2)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	2,  2, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	2,  0, 2)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	2,  0, 2)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	2,  2, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	2,  0, 2)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	2,  0, 2)\
	RCT(R_G_B2,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,  1, 1)\
	RCT(R_GR_B2,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	2,  1, 1)\
	RCT(R_B_G2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,  1, 1)\
	RCT(R_BR_G2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	2,  1, 1)\
	RCT(G_B_R2,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,  1, 1)\
	RCT(G_BG_R2,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	2,  1, 1)\
	RCT(G_RG_B2,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	2,  1, 1)\
	RCT(B_RB_G2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	2,  1, 1)\
	RCT(B_GB_R2,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	2,  1, 1)
typedef enum _RCTIndex
{
#define RCT(LABEL, ...) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][II_COUNT]=
{
#define RCT(LABEL, ...) {__VA_ARGS__},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, ...) #LABEL,
	RCTLIST
#undef  RCT
};

#define PREDLIST\
	PRED(ZERO)\
	PRED(W)\
	PRED(NW)\
	PRED(SELECT)\
	PRED(CG)\
	PRED(AV2)\
	PRED(IZ)\
	PRED(AV3)\
	PRED(AV4)\
	PRED(AV6)\
	PRED(AV8)\
	PRED(AV9)\
	PRED(WG)
typedef enum _PredIndex
{
#define PRED(LABEL) PRED_##LABEL,
	PREDLIST
#undef  PRED
	PRED_COUNT,
} PredIndex;
static const char *pred_names[PRED_COUNT]=
{
#define PRED(LABEL) #LABEL,
	PREDLIST
#undef  PRED
};
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
	unsigned char *imagebuf, *streambuf;
	int iw, ih;
	int *offsets, *chunksizes;//[nblocks*2+1]

	int fwd, loud, b1, b2, xblocks, nblocks, x1, x2, y1, y2, dist;
	int bufsize, histsize;
	short *pixels;
	int *hist;

	int tlevels;

	const unsigned *mixincdfs;
	int ctxctrs[3][ELEVELS][CLEVELS];//good

	//aux
	int blockidx;
	double bestsize;
	int bestrct, predidx[3];
#ifdef STATIC_ESTIMATE
	int ehistsize, *ehist;
	double esizes[3*(ELEVELS<<PBITS)*2];//3 channels  *  ELEVELS contexts  *  {streamsize, statsize}
#endif
#ifdef STATIC_ESTIMATE2
	int ehist2[3][ELEVELS][PLEVELS][33], ezerohist2[3][256];
	double esize2, hsize2;
#endif
#ifdef PRINT_CSIZES
	double csizes[3];
#endif
} ThreadArgs;
static void block_func(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	BypassCoder bc;
	AC3 ec;
	int cdfstride=args->tlevels+1;
	const unsigned char *image=args->imagebuf;
	unsigned char bestrct=0, predidx[3]={0};
	const unsigned char *combination=0;
	int entropylevel=0, entropyidx=0;
#ifdef PRINT_CSIZES
	double esizes[3]={0};
#endif

	if(args->fwd)
	{
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		unsigned char predsel[OCH_COUNT]={0};
		int dx=(args->x2-args->x1-3)/ANALYSIS_XSTRIDE/5*5, dy=(args->y2-args->y1-2)/ANALYSIS_YSTRIDE;
		int count=0;
		int ystride=args->iw*3;
#ifdef ENABLE_WG4
		int wg4_prederrors[4][OCH_COUNT]={0};
#endif
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
			entropylevel=99;
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
			__m128i half8=_mm_set1_epi8(-128);
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
				//1		rgb
				//2		gbr			unused exept for curr2
				//3		rgb - gbr
				//4		gbr - rgb
				//5 curr	(gbr+brg)/2
				//5 ...		rgb - (gbr+brg)/2
				__m256i
				//	NNWW,	NNWW3,	NNWW4,	NNWW5,
					NNW,	NNW3,	NNW4,	NNW5,
					NN,	NN3,	NN4,	NN5,
					NNE,	NNE3,	NNE4,	NNE5,
					NWW,	NWW3,	NWW4,	NWW5,
					NW,	NW3,	NW4,	NW5,
					N,	N3,	N4,	N5,
					NE,	NE3,	NE4,	NE5,
					NEE,	NEE3,	NEE4,	NEE5,
					WW,	WW3,	WW4,	WW5,
					W,	W3,	W4,	W5,
					curr,	curr2,		curr5;
				__m256i vmin[4], vmax[4], pred;
				{
				//	__m128i NNWW8	=_mm_loadu_si128((__m128i*)(ptr-2*ystride-2*3+0));
					__m128i NNW8	=_mm_loadu_si128((__m128i*)(ptr-2*ystride-1*3+0));
					__m128i NN8	=_mm_loadu_si128((__m128i*)(ptr-2*ystride+0*3+0));
					__m128i NNE8	=_mm_loadu_si128((__m128i*)(ptr-2*ystride+1*3+0));
					__m128i NWW8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride-2*3+0));
					__m128i NW8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride-1*3+0));
					__m128i N8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride+0*3+0));
					__m128i NE8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride+1*3+0));
					__m128i NEE8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride+2*3+0));
					__m128i WW8	=_mm_loadu_si128((__m128i*)(ptr+0*ystride-2*3+0));
					__m128i W8	=_mm_loadu_si128((__m128i*)(ptr+0*ystride-1*3+0));
					__m128i curr8	=_mm_loadu_si128((__m128i*)(ptr+0*ystride+0*3+0));
				//	NNWW8	=_mm_xor_si128(NNWW8	, half8);
					NNW8	=_mm_xor_si128(NNW8	, half8);
					NN8	=_mm_xor_si128(NN8	, half8);
					NNE8	=_mm_xor_si128(NNE8	, half8);
					NWW8	=_mm_xor_si128(NWW8	, half8);
					NW8	=_mm_xor_si128(NW8	, half8);
					N8	=_mm_xor_si128(N8	, half8);
					NE8	=_mm_xor_si128(NE8	, half8);
					NEE8	=_mm_xor_si128(NEE8	, half8);
					WW8	=_mm_xor_si128(WW8	, half8);
					W8	=_mm_xor_si128(W8	, half8);
					curr8	=_mm_xor_si128(curr8	, half8);
				//	__m128i NNWW82	=_mm_shuffle_epi8(NNWW8	, shuf);
					__m128i NNW82	=_mm_shuffle_epi8(NNW8	, shuf);
					__m128i NN82	=_mm_shuffle_epi8(NN8	, shuf);
					__m128i NNE82	=_mm_shuffle_epi8(NNE8	, shuf);
					__m128i NWW82	=_mm_shuffle_epi8(NWW8	, shuf);
					__m128i NW82	=_mm_shuffle_epi8(NW8	, shuf);
					__m128i N82	=_mm_shuffle_epi8(N8	, shuf);
					__m128i NE82	=_mm_shuffle_epi8(NE8	, shuf);
					__m128i NEE82	=_mm_shuffle_epi8(NEE8	, shuf);
					__m128i WW82	=_mm_shuffle_epi8(WW8	, shuf);
					__m128i W82	=_mm_shuffle_epi8(W8	, shuf);
					__m128i curr82	=_mm_shuffle_epi8(curr8	, shuf);
				//	NNWW	=_mm256_cvtepi8_epi16(NNWW8);
					NNW	=_mm256_cvtepi8_epi16(NNW8);
					NN	=_mm256_cvtepi8_epi16(NN8);
					NNE	=_mm256_cvtepi8_epi16(NNE8);
					NWW	=_mm256_cvtepi8_epi16(NWW8);
					NW	=_mm256_cvtepi8_epi16(NW8);
					N	=_mm256_cvtepi8_epi16(N8);
					NE	=_mm256_cvtepi8_epi16(NE8);
					NEE	=_mm256_cvtepi8_epi16(NEE8);
					WW	=_mm256_cvtepi8_epi16(WW8);
					W	=_mm256_cvtepi8_epi16(W8);
					curr	=_mm256_cvtepi8_epi16(curr8);
				//	__m256i NNWW2	=_mm256_cvtepi8_epi16(NNWW82);
					__m256i NNW2	=_mm256_cvtepi8_epi16(NNW82);
					__m256i NN2	=_mm256_cvtepi8_epi16(NN82);
					__m256i NNE2	=_mm256_cvtepi8_epi16(NNE82);
					__m256i NWW2	=_mm256_cvtepi8_epi16(NWW82);
					__m256i NW2	=_mm256_cvtepi8_epi16(NW82);
					__m256i N2	=_mm256_cvtepi8_epi16(N82);
					__m256i NE2	=_mm256_cvtepi8_epi16(NE82);
					__m256i NEE2	=_mm256_cvtepi8_epi16(NEE82);
					__m256i WW2	=_mm256_cvtepi8_epi16(WW82);
					__m256i W2	=_mm256_cvtepi8_epi16(W82);
					curr2	=_mm256_cvtepi8_epi16(curr82);
				//	NNWW3	=_mm256_sub_epi16(NNWW	, NNWW2	);
					NNW3	=_mm256_sub_epi16(NNW	, NNW2	);
					NN3	=_mm256_sub_epi16(NN	, NN2	);
					NNE3	=_mm256_sub_epi16(NNE	, NNE2	);
					NWW3	=_mm256_sub_epi16(NWW	, NWW2	);
					NW3	=_mm256_sub_epi16(NW	, NW2	);
					N3	=_mm256_sub_epi16(N	, N2	);
					NE3	=_mm256_sub_epi16(NE	, NE2	);
					NEE3	=_mm256_sub_epi16(NEE	, NEE2	);
					WW3	=_mm256_sub_epi16(WW	, WW2	);
					W3	=_mm256_sub_epi16(W	, W2	);
				//	curr3	=_mm256_sub_epi16(curr	, curr2	);
				//	NNWW4	=_mm256_sub_epi16(NNWW2	, NNWW	);
					NNW4	=_mm256_sub_epi16(NNW2	, NNW	);
					NN4	=_mm256_sub_epi16(NN2	, NN	);
					NNE4	=_mm256_sub_epi16(NNE2	, NNE	);
					NWW4	=_mm256_sub_epi16(NWW2	, NWW	);
					NW4	=_mm256_sub_epi16(NW2	, NW	);
					N4	=_mm256_sub_epi16(N2	, N	);
					NE4	=_mm256_sub_epi16(NE2	, NE	);
					NEE4	=_mm256_sub_epi16(NEE2	, NEE	);
					WW4	=_mm256_sub_epi16(WW2	, WW	);
					W4	=_mm256_sub_epi16(W2	, W	);
				//	curr4	=_mm256_sub_epi16(curr2	, curr	);
					
				//	NNWW5	=_mm256_add_epi16(NNWW2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NNWW82	, shuf)));
					NNW5	=_mm256_add_epi16(NNW2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NNW82	, shuf)));
					NN5	=_mm256_add_epi16(NN2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NN82	, shuf)));
					NNE5	=_mm256_add_epi16(NNE2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NNE82	, shuf)));
					NWW5	=_mm256_add_epi16(NWW2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NWW82	, shuf)));
					NW5	=_mm256_add_epi16(NW2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NW82	, shuf)));
					N5	=_mm256_add_epi16(N2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(N82	, shuf)));
					NE5	=_mm256_add_epi16(NE2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NE82	, shuf)));
					NEE5	=_mm256_add_epi16(NEE2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NEE82	, shuf)));
					WW5	=_mm256_add_epi16(WW2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(WW82	, shuf)));
					W5	=_mm256_add_epi16(W2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(W82	, shuf)));
					curr5	=_mm256_add_epi16(curr2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(curr82	, shuf)));
				//	NNWW5	=_mm256_srai_epi16(NNWW5, 1);
					NNW5	=_mm256_srai_epi16(NNW5	, 1);
					NN5	=_mm256_srai_epi16(NN5	, 1);
					NNE5	=_mm256_srai_epi16(NNE5	, 1);
					NWW5	=_mm256_srai_epi16(NWW5	, 1);
					NW5	=_mm256_srai_epi16(NW5	, 1);
					N5	=_mm256_srai_epi16(N5	, 1);
					NE5	=_mm256_srai_epi16(NE5	, 1);
					NEE5	=_mm256_srai_epi16(NEE5	, 1);
					WW5	=_mm256_srai_epi16(WW5	, 1);
					W5	=_mm256_srai_epi16(W5	, 1);
					curr5	=_mm256_srai_epi16(curr5, 1);
				//	NNWW5	=_mm256_sub_epi16(NNWW	, NNWW5);
					NNW5	=_mm256_sub_epi16(NNW	, NNW5);
					NN5	=_mm256_sub_epi16(NN	, NN5);
					NNE5	=_mm256_sub_epi16(NNE	, NNE5);
					NWW5	=_mm256_sub_epi16(NWW	, NWW5);
					NW5	=_mm256_sub_epi16(NW	, NW5);
					N5	=_mm256_sub_epi16(N	, N5);
					NE5	=_mm256_sub_epi16(NE	, NE5);
					NEE5	=_mm256_sub_epi16(NEE	, NEE5);
					WW5	=_mm256_sub_epi16(WW	, WW5);
					W5	=_mm256_sub_epi16(W	, W5);
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
				vmin[3]=_mm256_min_epi16(N5, W5);
				vmax[3]=_mm256_max_epi16(N5, W5);
#if defined ENABLE_AV9v2 || defined ENABLE_AV4v2 || defined ENABLE_WGB
				__m256i vmin0[3]=
				{
					vmin[0],
					vmin[1],
					vmin[2],
					vmin[3],
				};
				__m256i vmax0[3]=
				{
					vmax[0],
					vmax[1],
					vmax[2],
					vmax[3],
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
				
				//ZERO
#ifdef ENABLE_ZERO
				pred=curr;
				UPDATE(
					PRED_ZERO,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_sub_epi16(curr, curr2);
				UPDATE(
					PRED_ZERO,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_sub_epi16(curr2, curr);
				UPDATE(
					PRED_ZERO,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_sub_epi16(curr, curr5);
				UPDATE(
					PRED_ZERO,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);
#endif

				//CLEARTYPE
				//rNW gNW bNW rN gN bN
				//rW  gW  bW  r  g  b
				//pred(r) = SELECT(rN, bW, bNW)
				//pred(g) = SELECT(gN, r, rN)
				//pred(b) = SELECT(bN, g, gN)
#ifdef ENABLE_CLEARTYPE
				{
					__m256i blendmask=_mm256_set_epi32(
						0x0000FFFF,
						0xFFFF0000,
						0xFFFFFFFF,
						0x0000FFFF,
						0xFFFF0000,
						0xFFFFFFFF,
						0x0000FFFF,
						0xFFFF0000
					);
					__m256i sNW, sW;

					sNW=_mm256_blendv_epi8(W, curr, blendmask);
					sW=_mm256_blendv_epi8(W, NW, N);
					pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N, sNW)), _mm256_abs_epi16(_mm256_sub_epi16(sW, sNW)));
					pred=_mm256_blendv_epi8(sW, N, pred);

					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_CLEARTYPE,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					sNW=_mm256_blendv_epi8(W3, curr, blendmask);
					sW=_mm256_blendv_epi8(W3, NW3, N3);
					pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N3, sNW)), _mm256_abs_epi16(_mm256_sub_epi16(sW, sNW)));
					pred=_mm256_blendv_epi8(sW, N3, pred);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_CLEARTYPE,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					sNW=_mm256_blendv_epi8(W4, curr, blendmask);
					sW=_mm256_blendv_epi8(W4, NW4, N4);
					pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N4, sNW)), _mm256_abs_epi16(_mm256_sub_epi16(sW, sNW)));
					pred=_mm256_blendv_epi8(sW, N4, pred);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_CLEARTYPE,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					sNW=_mm256_blendv_epi8(W5, curr, blendmask);
					sW=_mm256_blendv_epi8(W5, NW5, N5);
					pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N5, sNW)), _mm256_abs_epi16(_mm256_sub_epi16(sW, sNW)));
					pred=_mm256_blendv_epi8(sW, N5, pred);

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_CLEARTYPE,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
				}
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
				pred=_mm256_add_epi16(W5, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_W,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
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
				pred=NW5;

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_NW,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
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
				pred=NE5;

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_NE,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
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
				pred=_mm256_srai_epi16(_mm256_add_epi16(N5, W5), 1);

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
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

				//SELECT = abs(N-NW)>abs(W-NW) ? N : W
#ifdef ENABLE_SELECT
				pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N, NW)), _mm256_abs_epi16(_mm256_sub_epi16(W, NW)));
				pred=_mm256_blendv_epi8(W, N, pred);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_SELECT,
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
					PRED_SELECT,
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
					PRED_SELECT,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N5, NW5)), _mm256_abs_epi16(_mm256_sub_epi16(W5, NW5)));
				pred=_mm256_blendv_epi8(W5, N5, pred);

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_SELECT,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
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

				//CGX = median(N, W, W+(N-NW)/2)
#ifdef ENABLE_CGX
				pred=_mm256_add_epi16(W, _mm256_srai_epi16(_mm256_sub_epi16(N, NW), 1));
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CGX,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(W3, _mm256_srai_epi16(_mm256_sub_epi16(N3, NW3), 1));
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CGX,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(W4, _mm256_srai_epi16(_mm256_sub_epi16(N4, NW4), 1));
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_CGX,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_add_epi16(W5, _mm256_srai_epi16(_mm256_sub_epi16(N5, NW5), 1));
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CGX,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);
#endif

				//CGY = median(N, W, N+(W-NW)/2)
#ifdef ENABLE_CGY
				pred=_mm256_add_epi16(N, _mm256_srai_epi16(_mm256_sub_epi16(W, NW), 1));
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CGY,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(N3, _mm256_srai_epi16(_mm256_sub_epi16(W3, NW3), 1));
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CGY,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(N4, _mm256_srai_epi16(_mm256_sub_epi16(W4, NW4), 1));
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_CGY,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_add_epi16(N5, _mm256_srai_epi16(_mm256_sub_epi16(W5, NW5), 1));
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CGY,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
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
				pred=_mm256_sub_epi16(_mm256_add_epi16(N5, W5), NW5);
#ifdef ENABLE_CGv2
				cgpred[1]=pred;
#endif
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CG,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);
#endif
				vmin[0]=_mm256_min_epi16(vmin[0], NE);
				vmax[0]=_mm256_max_epi16(vmax[0], NE);
				vmin[1]=_mm256_min_epi16(vmin[1], NE3);
				vmax[1]=_mm256_max_epi16(vmax[1], NE3);
				vmin[2]=_mm256_min_epi16(vmin[2], NE4);
				vmax[2]=_mm256_max_epi16(vmax[2], NE4);
				vmin[3]=_mm256_min_epi16(vmin[3], NE5);
				vmax[3]=_mm256_max_epi16(vmax[3], NE5);

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
				pred=_mm256_add_epi16(N5, W5);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_sub_epi16(pred, _mm256_slli_epi16(NW5, 1));
				pred=_mm256_srai_epi16(pred, 3);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV3,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
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
				pred=_mm256_add_epi16(N5, W5);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(_mm256_sub_epi16(pred, NW5), 1));
				pred=_mm256_srai_epi16(pred, 2);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_IZ,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
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
				pred=_mm256_add_epi16(N5, W5);
				pred=_mm256_slli_epi16(pred, 2);
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE5, NW5));
				pred=_mm256_srai_epi16(pred, 3);
#ifdef ENABLE_AV4v2
				av4pred[1]=pred;
#endif
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV4,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
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
				pred=_mm256_sub_epi16(N5, NW5);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE5, WW5));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), W5);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, curr5);
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
#endif
				
				//AV6
				//			-1
				//		-5	6	1
				//	-1	8	[?]>>3		clamp(N,W,NE)
#ifdef ENABLE_AV6
				pred=_mm256_sub_epi16(N, NW);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));//5*(N-NW)
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE, WW));//5*(N-NW)+NE-WW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(N, NN));//6*N-5*NW-NN-WW+NE
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), W);//W+(6*N-5*NW-NN-WW+NE)/8
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV6,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_sub_epi16(N3, NW3);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));//5*(N-NW)
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE3, WW3));//5*(N-NW)+NE-WW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(N3, NN3));//6*N-5*NW-NN-WW+NE
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), W3);//W+(6*N-5*NW-NN-WW+NE)/8
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV6,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_sub_epi16(N4, NW4);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));//5*(N-NW)
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE4, WW4));//5*(N-NW)+NE-WW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(N4, NN4));//6*N-5*NW-NN-WW+NE
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), W4);//W+(6*N-5*NW-NN-WW+NE)/8
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_AV6,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_sub_epi16(N5, NW5);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));//5*(N-NW)
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE5, WW5));//5*(N-NW)+NE-WW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(N5, NN5));//6*N-5*NW-NN-WW+NE
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), W5);//W+(6*N-5*NW-NN-WW+NE)/8
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV6,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);
#endif
				
				//AV8
				//			1	1
				//		1	1	1	1
				//	1	1	[?]>>3
#ifdef ENABLE_AV8
				pred=_mm256_add_epi16(N, W);
				pred=_mm256_add_epi16(pred, _mm256_add_epi16(NN, WW));
				pred=_mm256_add_epi16(pred, _mm256_add_epi16(NW, NE));
				pred=_mm256_add_epi16(pred, _mm256_add_epi16(NNE, NEE));
				pred=_mm256_srai_epi16(pred, 3);

				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV8,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(N3, W3);
				pred=_mm256_add_epi16(pred, _mm256_add_epi16(NN3, WW3));
				pred=_mm256_add_epi16(pred, _mm256_add_epi16(NW3, NE3));
				pred=_mm256_add_epi16(pred, _mm256_add_epi16(NNE3, NEE3));
				pred=_mm256_srai_epi16(pred, 3);

				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV8,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(N4, W4);
				pred=_mm256_add_epi16(pred, _mm256_add_epi16(NN4, WW4));
				pred=_mm256_add_epi16(pred, _mm256_add_epi16(NW4, NE4));
				pred=_mm256_add_epi16(pred, _mm256_add_epi16(NNE4, NEE4));
				pred=_mm256_srai_epi16(pred, 3);

				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_AV8,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_add_epi16(N5, W5);
				pred=_mm256_add_epi16(pred, _mm256_add_epi16(NN5, WW5));
				pred=_mm256_add_epi16(pred, _mm256_add_epi16(NW5, NE5));
				pred=_mm256_add_epi16(pred, _mm256_add_epi16(NNE5, NEE5));
				pred=_mm256_srai_epi16(pred, 3);

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV8,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
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
				pred=_mm256_add_epi16(N5, _mm256_slli_epi16(N5, 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(NN5, WW5));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(NE5, 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(NW5, 3), NW5));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NNW5, _mm256_add_epi16(NNE5, NWW5)));
				pred=_mm256_add_epi16(W5, _mm256_srai_epi16(pred, 4));
#ifdef ENABLE_AV9v2
				av9pred[1]=pred;
#endif
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, curr5);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV9,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
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

				//AV11 = clamp((
				//	+4*NNWW	+3*NNW	-31*NN	-38*NNE
				//	+7*NWW	-158*NW	+219*N	+30*NE	+19*NEE
				//	-42*WW	+243*W
				//+128)>>8, N,W,NE)
				//
				//at analysis:
				//AV11 = clamp((
				//			-2*NN	-2*NNE
				//		-10*NW	+14*N	+2*NE	+1*NEE
				//	-3*WW	+16*W
				//)>>4, N,W,NE)
#ifdef ENABLE_AV11
				{
					__m256i threeW, t0, t1;
					
					//			-2	-2
					//		-8-2	16-2	+2	+1
					//	-2-1	+16	?
					threeW=_mm256_add_epi16(W, _mm256_srai_epi16(W, 1));
					pred=_mm256_add_epi16(W, NEE);
					pred=_mm256_add_epi16(pred, NN);
					pred=_mm256_sub_epi16(pred, NNW);
					pred=_mm256_sub_epi16(pred, NWW);
					pred=_mm256_add_epi16(pred, threeW);
					t0=_mm256_add_epi16(W, NW);
					t0=_mm256_sub_epi16(t0, NE);
					t0=_mm256_sub_epi16(t0, NNE);
					t0=_mm256_add_epi16(t0, NEE);
					t0=_mm256_sub_epi16(t0, NNW);
					t0=_mm256_add_epi16(t0, WW);
					t1=_mm256_add_epi16(NNWW, NNE);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t0);	//	pred/128
					t0=_mm256_add_epi16(WW, NWW);
					t0=_mm256_add_epi16(t0, threeW);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t1);	//	pred/64
					t1=_mm256_sub_epi16(NEE, W);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t0);	//	pred/32
					t0=_mm256_add_epi16(NN, NNE);
					t0=_mm256_add_epi16(t0, NW);
					t0=_mm256_add_epi16(t0, WW);
					t0=_mm256_sub_epi16(NE, t0);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t1);		//	pred/16
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t0);		//	pred/8
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), threeW);	//	pred/4
					pred=_mm256_sub_epi16(_mm256_srai_epi16(pred, 1), NW);		//	pred/2
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), W);		//	pred/1

					//pred=_mm256_add_epi16(N, W);
					//pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), NW);
					//t0=_mm256_add_epi16(NW, N);
					//t1=_mm256_add_epi16(NN, NNE);
					//t0=_mm256_sub_epi16(NE, t0);
					//t1=_mm256_add_epi16(t1, WW);
					//pred=_mm256_add_epi16(_mm256_slli_epi16(pred, 2), _mm256_sub_epi16(t0, t1));
					//pred=_mm256_add_epi16(_mm256_slli_epi16(pred, 1), _mm256_sub_epi16(NEE, WW));
					pred=_mm256_max_epi16(pred, vmin[0]);
					pred=_mm256_min_epi16(pred, vmax[0]);

					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV11,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					threeW=_mm256_add_epi16(W3, _mm256_srai_epi16(W3, 1));
					pred=_mm256_add_epi16(W3, NEE3);
					pred=_mm256_add_epi16(pred, NN3);
					pred=_mm256_sub_epi16(pred, NNW3);
					pred=_mm256_sub_epi16(pred, NWW3);
					pred=_mm256_add_epi16(pred, threeW);
					t0=_mm256_add_epi16(W3, NW3);
					t0=_mm256_sub_epi16(t0, NE3);
					t0=_mm256_sub_epi16(t0, NNE3);
					t0=_mm256_add_epi16(t0, NEE3);
					t0=_mm256_sub_epi16(t0, NNW3);
					t0=_mm256_add_epi16(t0, WW3);
					t1=_mm256_add_epi16(NNWW3, NNE3);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t0);
					t0=_mm256_add_epi16(WW3, NWW3);
					t0=_mm256_add_epi16(t0, threeW);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t1);
					t1=_mm256_sub_epi16(NEE3, W3);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t0);
					t0=_mm256_add_epi16(NN3, NNE3);
					t0=_mm256_add_epi16(t0, NW3);
					t0=_mm256_add_epi16(t0, WW3);
					t0=_mm256_sub_epi16(NE3, t0);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t1);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t0);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), threeW);
					pred=_mm256_sub_epi16(_mm256_srai_epi16(pred, 1), NW3);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), W3);

					//pred=_mm256_add_epi16(N3, W3);
					//pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), NW3);
					//t0=_mm256_add_epi16(NW3, N3);
					//t1=_mm256_add_epi16(NN3, NNE3);
					//t0=_mm256_sub_epi16(NE3, t0);
					//t1=_mm256_add_epi16(t1, WW3);
					//pred=_mm256_add_epi16(_mm256_slli_epi16(pred, 2), _mm256_sub_epi16(t0, t1));
					//pred=_mm256_add_epi16(_mm256_slli_epi16(pred, 1), _mm256_sub_epi16(NEE3, WW3));
					pred=_mm256_max_epi16(pred, vmin[1]);
					pred=_mm256_min_epi16(pred, vmax[1]);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV11,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					threeW=_mm256_add_epi16(W4, _mm256_srai_epi16(W4, 1));
					pred=_mm256_add_epi16(W4, NEE4);
					pred=_mm256_add_epi16(pred, NN4);
					pred=_mm256_sub_epi16(pred, NNW4);
					pred=_mm256_sub_epi16(pred, NWW4);
					pred=_mm256_add_epi16(pred, threeW);
					t0=_mm256_add_epi16(W4, NW4);
					t0=_mm256_sub_epi16(t0, NE4);
					t0=_mm256_sub_epi16(t0, NNE4);
					t0=_mm256_add_epi16(t0, NEE4);
					t0=_mm256_sub_epi16(t0, NNW4);
					t0=_mm256_add_epi16(t0, WW4);
					t1=_mm256_add_epi16(NNWW4, NNE4);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t0);
					t0=_mm256_add_epi16(WW4, NWW4);
					t0=_mm256_add_epi16(t0, threeW);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t1);
					t1=_mm256_sub_epi16(NEE4, W4);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t0);
					t0=_mm256_add_epi16(NN4, NNE4);
					t0=_mm256_add_epi16(t0, NW4);
					t0=_mm256_add_epi16(t0, WW4);
					t0=_mm256_sub_epi16(NE4, t0);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t1);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t0);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), threeW);
					pred=_mm256_sub_epi16(_mm256_srai_epi16(pred, 1), NW4);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), W4);
					//pred=_mm256_add_epi16(N4, W4);
					//pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), NW4);
					//t0=_mm256_add_epi16(NW4, N4);
					//t1=_mm256_add_epi16(NN4, NNE4);
					//t0=_mm256_sub_epi16(NE4, t0);
					//t1=_mm256_add_epi16(t1, WW4);
					//pred=_mm256_add_epi16(_mm256_slli_epi16(pred, 2), _mm256_sub_epi16(t0, t1));
					//pred=_mm256_add_epi16(_mm256_slli_epi16(pred, 1), _mm256_sub_epi16(NEE4, WW4));
					pred=_mm256_max_epi16(pred, vmin[2]);
					pred=_mm256_min_epi16(pred, vmax[2]);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_AV11,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					threeW=_mm256_add_epi16(W5, _mm256_srai_epi16(W5, 1));
					pred=_mm256_add_epi16(W5, NEE5);
					pred=_mm256_add_epi16(pred, NN5);
					pred=_mm256_sub_epi16(pred, NNW5);
					pred=_mm256_sub_epi16(pred, NWW5);
					pred=_mm256_add_epi16(pred, threeW);
					t0=_mm256_add_epi16(W5, NW5);
					t0=_mm256_sub_epi16(t0, NE5);
					t0=_mm256_sub_epi16(t0, NNE5);
					t0=_mm256_add_epi16(t0, NEE5);
					t0=_mm256_sub_epi16(t0, NNW5);
					t0=_mm256_add_epi16(t0, WW5);
					t1=_mm256_add_epi16(NNWW5, NNE5);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t0);
					t0=_mm256_add_epi16(WW5, NWW5);
					t0=_mm256_add_epi16(t0, threeW);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t1);
					t1=_mm256_sub_epi16(NEE, W5);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t0);
					t0=_mm256_add_epi16(NN5, NNE5);
					t0=_mm256_add_epi16(t0, NW5);
					t0=_mm256_add_epi16(t0, WW5);
					t0=_mm256_sub_epi16(NE5, t0);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t1);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), t0);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), threeW);
					pred=_mm256_sub_epi16(_mm256_srai_epi16(pred, 1), NW5);
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 1), W5);
					//pred=_mm256_add_epi16(N5, W5);
					//pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), NW5);
					//t0=_mm256_add_epi16(NW5, N5);
					//t1=_mm256_add_epi16(NN5, NNE5);
					//t0=_mm256_sub_epi16(NE5, t0);
					//t1=_mm256_add_epi16(t1, WW5);
					//pred=_mm256_add_epi16(_mm256_slli_epi16(pred, 2), _mm256_sub_epi16(t0, t1));
					//pred=_mm256_add_epi16(_mm256_slli_epi16(pred, 1), _mm256_sub_epi16(NEE5, WW5));
					pred=_mm256_max_epi16(pred, vmin[3]);
					pred=_mm256_min_epi16(pred, vmax[3]);

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV11,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
				}
#endif
				
				//WG
				//gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1
				//gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1
				//pred=(gx*N+gy*W)/(gx+gy)
#ifdef ENABLE_WG
				{
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
					gx=_mm256_abs_epi16(_mm256_sub_epi16(W5, WW5));
					gy=_mm256_abs_epi16(_mm256_sub_epi16(W5, NW5));
					gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N5, NW5)));
					gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(N5, NN5)));
					gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE5, N5)));
					gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE5, NNE5)));
					gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
					gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
					MIX2_16x16(pred, N5, W5, gx, gy);

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_WG,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
				}
#endif
				
				//WGB
				//gx=abs(N-NW)+abs(W-WW)+abs(NE-N)+abs(NW-NWW)+1
				//gy=abs(N-NN)+abs(W-NW)+abs(NE-NNE)+abs(NW-NNW)+1
				//sum=gx+gy
				//sh=(depth+3-FLOOR_LOG2(sum))>>1
				//h1=gx>>sh
				//h2=gy>>sh
				//pred=CLAMP(((gx+h1)*N+(gy+h1)*W-(h1+h2)*NW)/sum, N,W)
#ifdef ENABLE_WGB
				{
					__m256i gx, gy, hx, hy;
					gx=_mm256_abs_epi16(_mm256_sub_epi16(W, WW));
					gy=_mm256_abs_epi16(_mm256_sub_epi16(W, NW));
					gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N, NW)));
					gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(N, NN)));
					gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE, N)));
					gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE, NNE)));
					gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NW, NWW)));
					gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NW, NNW)));
					gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
					gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
					hx=_mm256_srli_epi16(gx, 3);
					hy=_mm256_srli_epi16(gy, 3);
					gx=_mm256_add_epi16(gx, hx);
					gy=_mm256_add_epi16(gy, hy);
					hx=_mm256_add_epi16(hx, hy);
					{
						__m256 a0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(N, 16), 16));
						__m256 a0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(N, 16));
						__m256 a1lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(W, 16), 16));
						__m256 a1hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(W, 16));
						__m256 a2lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(NW, 16), 16));
						__m256 a2hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(NW, 16));
						__m256 c0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(gx, 16), 16));
						__m256 c0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(gx, 16));
						__m256 c1lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(gy, 16), 16));
						__m256 c1hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(gy, 16));
						__m256 c2lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(hx, 16), 16));
						__m256 c2hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(hx, 16));
						a0lo=_mm256_mul_ps(a0lo, c0lo);
						a0hi=_mm256_mul_ps(a0hi, c0hi);
						a1lo=_mm256_mul_ps(a1lo, c1lo);
						a1hi=_mm256_mul_ps(a1hi, c1hi);
						a2lo=_mm256_mul_ps(a2lo, c2lo);
						a2hi=_mm256_mul_ps(a2hi, c2hi);
						c0lo=_mm256_add_ps(c0lo, c1lo);
						c0hi=_mm256_add_ps(c0hi, c1hi);
						a0lo=_mm256_add_ps(a0lo, a1lo);
						a0hi=_mm256_add_ps(a0hi, a1hi);
						a0lo=_mm256_sub_ps(a0lo, a2lo);
						a0hi=_mm256_sub_ps(a0hi, a2hi);
						a0lo=_mm256_div_ps(a0lo, c0lo);
						a0hi=_mm256_div_ps(a0hi, c0hi);
						pred=_mm256_and_si256(_mm256_cvttps_epi32(a0lo), _mm256_set1_epi32(0xFFFF));
						pred=_mm256_or_si256(pred, _mm256_slli_epi32(_mm256_cvttps_epi32(a0hi), 16));
					}
					pred=_mm256_max_epi16(pred, vmin0[0]);
					pred=_mm256_min_epi16(pred, vmax0[0]);

					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_WGB,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					gx=_mm256_abs_epi16(_mm256_sub_epi16(W3, WW3));
					gy=_mm256_abs_epi16(_mm256_sub_epi16(W3, NW3));
					gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N3, NW3)));
					gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(N3, NN3)));
					gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE3, N3)));
					gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE3, NNE3)));
					gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NW3, NWW3)));
					gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NW3, NNW3)));
					gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
					gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
					hx=_mm256_srli_epi16(gx, 3);
					hy=_mm256_srli_epi16(gy, 3);
					gx=_mm256_add_epi16(gx, hx);
					gy=_mm256_add_epi16(gy, hy);
					hx=_mm256_add_epi16(hx, hy);
					{
						__m256 a0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(N3, 16), 16));
						__m256 a0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(N3, 16));
						__m256 a1lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(W3, 16), 16));
						__m256 a1hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(W3, 16));
						__m256 a2lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(NW3, 16), 16));
						__m256 a2hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(NW3, 16));
						__m256 c0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(gx, 16), 16));
						__m256 c0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(gx, 16));
						__m256 c1lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(gy, 16), 16));
						__m256 c1hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(gy, 16));
						__m256 c2lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(hx, 16), 16));
						__m256 c2hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(hx, 16));
						a0lo=_mm256_mul_ps(a0lo, c0lo);
						a0hi=_mm256_mul_ps(a0hi, c0hi);
						a1lo=_mm256_mul_ps(a1lo, c1lo);
						a1hi=_mm256_mul_ps(a1hi, c1hi);
						a2lo=_mm256_mul_ps(a2lo, c2lo);
						a2hi=_mm256_mul_ps(a2hi, c2hi);
						c0lo=_mm256_add_ps(c0lo, c1lo);
						c0hi=_mm256_add_ps(c0hi, c1hi);
						a0lo=_mm256_add_ps(a0lo, a1lo);
						a0hi=_mm256_add_ps(a0hi, a1hi);
						a0lo=_mm256_sub_ps(a0lo, a2lo);
						a0hi=_mm256_sub_ps(a0hi, a2hi);
						a0lo=_mm256_div_ps(a0lo, c0lo);
						a0hi=_mm256_div_ps(a0hi, c0hi);
						pred=_mm256_and_si256(_mm256_cvttps_epi32(a0lo), _mm256_set1_epi32(0xFFFF));
						pred=_mm256_or_si256(pred, _mm256_slli_epi32(_mm256_cvttps_epi32(a0hi), 16));
					}
					pred=_mm256_max_epi16(pred, vmin0[1]);
					pred=_mm256_min_epi16(pred, vmax0[1]);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_WGB,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					gx=_mm256_abs_epi16(_mm256_sub_epi16(W4, WW4));
					gy=_mm256_abs_epi16(_mm256_sub_epi16(W4, NW4));
					gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N4, NW4)));
					gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(N4, NN4)));
					gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE4, N4)));
					gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE4, NNE4)));
					gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NW4, NWW4)));
					gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NW4, NNW4)));
					gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
					gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
					hx=_mm256_srli_epi16(gx, 3);
					hy=_mm256_srli_epi16(gy, 3);
					gx=_mm256_add_epi16(gx, hx);
					gy=_mm256_add_epi16(gy, hy);
					hx=_mm256_add_epi16(hx, hy);
					{
						__m256 a0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(N4, 16), 16));
						__m256 a0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(N4, 16));
						__m256 a1lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(W4, 16), 16));
						__m256 a1hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(W4, 16));
						__m256 a2lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(NW4, 16), 16));
						__m256 a2hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(NW4, 16));
						__m256 c0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(gx, 16), 16));
						__m256 c0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(gx, 16));
						__m256 c1lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(gy, 16), 16));
						__m256 c1hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(gy, 16));
						__m256 c2lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(hx, 16), 16));
						__m256 c2hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(hx, 16));
						a0lo=_mm256_mul_ps(a0lo, c0lo);
						a0hi=_mm256_mul_ps(a0hi, c0hi);
						a1lo=_mm256_mul_ps(a1lo, c1lo);
						a1hi=_mm256_mul_ps(a1hi, c1hi);
						a2lo=_mm256_mul_ps(a2lo, c2lo);
						a2hi=_mm256_mul_ps(a2hi, c2hi);
						c0lo=_mm256_add_ps(c0lo, c1lo);
						c0hi=_mm256_add_ps(c0hi, c1hi);
						a0lo=_mm256_add_ps(a0lo, a1lo);
						a0hi=_mm256_add_ps(a0hi, a1hi);
						a0lo=_mm256_sub_ps(a0lo, a2lo);
						a0hi=_mm256_sub_ps(a0hi, a2hi);
						a0lo=_mm256_div_ps(a0lo, c0lo);
						a0hi=_mm256_div_ps(a0hi, c0hi);
						pred=_mm256_and_si256(_mm256_cvttps_epi32(a0lo), _mm256_set1_epi32(0xFFFF));
						pred=_mm256_or_si256(pred, _mm256_slli_epi32(_mm256_cvttps_epi32(a0hi), 16));
					}
					pred=_mm256_max_epi16(pred, vmin0[2]);
					pred=_mm256_min_epi16(pred, vmax0[2]);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_WGB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
				}
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

				//WG4
				//analysis = mix(
				//	N,
				//	W,
				//	W+NE-N,
				//	(NN+NNE+WW+NEE)>>2
				//)
				//pred = mix(
				//	N,
				//	W,
				//	3*(N-NN)+NNN,
				//	3*(W-WW)+WWW,
				//	W+NE-N,
				//	N+W-NW,
				//	N+NE-NNE,
				//	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4
				//)
#ifdef ENABLE_WG4
#define WG4_MIX4_16x16(DST, C0, C1, C2, C3, V0, V1, V2, V3)\
	do\
	{\
		__m256 c0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(C0, 16), 16));\
		__m256 c0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(C0, 16));\
		__m256 c0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(C1, 16), 16));\
		__m256 c0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(C1, 16));\
		__m256 c0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(C2, 16), 16));\
		__m256 c0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(C2, 16));\
		__m256 c0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(C3, 16), 16));\
		__m256 c0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(C3, 16));\
		__m256 v0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(V0, 16), 16));\
		__m256 v0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(V0, 16));\
		__m256 v0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(V1, 16), 16));\
		__m256 v0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(V1, 16));\
		__m256 v0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(V2, 16), 16));\
		__m256 v0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(V2, 16));\
		__m256 v0lo=_mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_slli_epi32(V3, 16), 16));\
		__m256 v0hi=_mm256_cvtepi32_ps(_mm256_srai_epi32(V3, 16));\
		v0lo=_mm256_mul_ps(v0lo, c0lo);\
		v0hi=_mm256_mul_ps(v0hi, c0hi);\
		v1lo=_mm256_mul_ps(v1lo, c1lo);\
		v1hi=_mm256_mul_ps(v1hi, c1hi);\
		v2lo=_mm256_mul_ps(v2lo, c2lo);\
		v2hi=_mm256_mul_ps(v2hi, c2hi);\
		v3lo=_mm256_mul_ps(v3lo, c3lo);\
		v3hi=_mm256_mul_ps(v3hi, c3hi);\
		c0lo=_mm256_add_ps(c0lo, c1lo);\
		c0hi=_mm256_add_ps(c0hi, c1hi);\
		c0lo=_mm256_add_ps(c0lo, c2lo);\
		c0hi=_mm256_add_ps(c0hi, c2hi);\
		c0lo=_mm256_add_ps(c0lo, c3lo);\
		c0hi=_mm256_add_ps(c0hi, c3hi);\
		v0lo=_mm256_add_ps(v0lo, v1lo);\
		v0hi=_mm256_add_ps(v0hi, v1hi);\
		v0lo=_mm256_add_ps(v0lo, v2lo);\
		v0hi=_mm256_add_ps(v0hi, v2hi);\
		v0lo=_mm256_add_ps(v0lo, v3lo);\
		v0hi=_mm256_add_ps(v0hi, v3hi);\
		v0lo=_mm256_div_ps(v0lo, c0lo);\
		v0hi=_mm256_div_ps(v0hi, c0hi);\
		DST=_mm256_blend_epi32(_mm256_cvttps_epi32(_a0lo), _mm256_slli_epi32(_mm256_cvttps_epi32(_a0hi), 16), 0xAA);\
	}while(0)
				{
					__m256i e0, e1, e2, e3;

					e0=N;
					e1=W;
					e2=_mm256_sub_epi16(_mm256_add_epi16(W, NE), N);
					e3=_mm256_srai_epi16(_mm256_add_epi16(_mm256_add_epi16(NN, NNE), _mm256_add_epi16(WW, NEE)), 2);
					//__m256i c0, c1, c2, c3;
					//
					//c0=_mm256_set_epi32(
					//	0, 0,
					//	wg4_prederrors[0][2],
					//	wg4_prederrors[0][1],
					//	wg4_prederrors[0][0],
					//	wg4_prederrors[0][2],
					//	wg4_prederrors[0][1],
					//	wg4_prederrors[0][0]
					//);
				}
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
			if(args->dist>1)
				bestpred=1;
			for(int kp=bestpred+1;kp<PRED_COUNT;++kp)
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
#ifdef ENABLE_CLEARTYPE
		if(predidx[0]==PRED_CLEARTYPE)
			predidx[0]=PRED_SELECT;
		if(predidx[1]==PRED_CLEARTYPE)
			predidx[1]=PRED_SELECT;
		if(predidx[2]==PRED_CLEARTYPE)
			predidx[2]=PRED_SELECT;
#endif
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
		ac3_encbuf_init(&ec, args->streambuf+args->offsets[args->blockidx*2+0], args->streambuf+args->offsets[args->blockidx*2+1]);
		bypass_encbuf_init(&bc, args->streambuf+args->offsets[args->blockidx*2+1], args->streambuf+args->offsets[args->blockidx*2+2]);
		//{
		//	int prevblockidx=args->blockidx-args->b1-1;
		//	ac3_encbuf_init(&ec, args->enctokenbuf+(prevblockidx>=0?args->enctokenoffsets[prevblockidx]:0), args->enctokenbuf+args->encbufsize);
		//	bypass_encbuf_init(&bc, args->encbypassbuf+(prevblockidx>=0?args->encbypassoffsets[prevblockidx]:0), args->encbypassbuf+args->encbufsize);
		//}
	//	ac3_encbuf_bypass_NPOT(&ec, args->dist, 2);
		ac3_encbuf_bypass_NPOT(&ec, entropylevel, 100);
		ac3_encbuf_bypass_NPOT(&ec, bestrct, RCT_COUNT);
		ac3_encbuf_bypass_NPOT(&ec, predidx[0], PRED_COUNT);
		ac3_encbuf_bypass_NPOT(&ec, predidx[1], PRED_COUNT);
		ac3_encbuf_bypass_NPOT(&ec, predidx[2], PRED_COUNT);

		args->bestrct=bestrct;
		args->predidx[0]=predidx[0];
		args->predidx[1]=predidx[1];
		args->predidx[2]=predidx[2];
#if defined ENABLE_ZERO && defined STATIC_ZERO
		if(predidx[0]==PRED_ZERO||predidx[1]==PRED_ZERO||predidx[2]==PRED_ZERO)
		{
			int hist0[3][256]={0};
			int yidx=combination[II_PERM_Y];
			int uidx=combination[II_PERM_U];
			int vidx=combination[II_PERM_V];
			int helpcoeffs[]=
			{
				0,
				combination[II_HELP_U],
				combination[II_HELP_V0],

				0,
				0,
				combination[II_HELP_V1],
			};
			for(int ky=args->y1;ky<args->y2;++ky)
			{
				int kx=args->x1;
				const unsigned char *ptr=image+3*(args->iw*ky+kx);
				for(;kx<args->x2;++kx, ptr+=3)
				{
					char yuv[]=
					{
						ptr[yidx]-128,
						ptr[uidx]-128,
						ptr[vidx]-128,
					};
					yuv[2]-=(helpcoeffs[2]*yuv[0]+helpcoeffs[5]*yuv[1])>>1;
					yuv[1]-=helpcoeffs[1]*yuv[0]>>1;
					++hist0[0][(yuv[0]+128)&255];
					++hist0[1][(yuv[1]+128)&255];
					++hist0[2][(yuv[2]+128)&255];
				}
			}
			for(int kc=0;kc<3;++kc)
			{
				if(predidx[kc]==PRED_ZERO)
				{
					int *curr_hist0=hist0[kc];
					int sum=0;
					for(int ks=0;ks<256;++ks)
						sum+=curr_hist0[ks];
					int sum2=0;
					int prev=1;
					for(int ks=0;ks<256;++ks)
					{
						int freq=curr_hist0[ks];
						curr_hist0[ks]=(int)(((long long)sum2*((1ULL<<PROB_BITS)-256))/sum)+ks;//adaptive formula because blocksize>0x10000	FIXME use reciprocal because sum is const
						sum2+=freq;

						if(ks)//(NLEVELS-1) increasing values
						{
							int nlevels=(1<<PROB_BITS)-prev;
							int sym=curr_hist0[ks]-prev;
							//if(nlevels>=0x10000)
							//{
							//	ac3_encbuf_bypass_NPOT(&ec, sym>>8, nlevels>>8);
							//	sym&=255;
							//	nlevels=256;
							//}
							ac3_encbuf_bypass_NPOT(&ec, sym, nlevels);
						}
						prev=curr_hist0[ks]+1;
					}
					unsigned *dsthist=(unsigned*)args->hist+ELEVELS*CLEVELS*cdfstride*kc;
					memcpy(dsthist, curr_hist0, sizeof(int[256]));
					dsthist[256]=1<<PROB_BITS;
				}
			}
		}
#endif
#ifdef PRINT_CSIZES
		printf("BLOCK %3d XY %5d %5d  %-10s %-10s %-10s %-10s  %11.2lf %11.2lf %11.2lf  ",
			args->blockidx,
			args->x1, args->y1,
			rct_names[bestrct],
			pred_names[predidx[0]],
			pred_names[predidx[1]],
			pred_names[predidx[2]],
			csizes[combination[0]*PRED_COUNT+predidx[0]],
			csizes[combination[1]*PRED_COUNT+predidx[1]],
			csizes[combination[2]*PRED_COUNT+predidx[2]]
		);
		//if(!args->blockidx)
		//{
		//	printf("Block0\n%8d count\n%11.2lf Y est.\n%11.2lf U est.\n%11.2lf V est.\n",
		//		count,
		//		csizes[combination[0]*PRED_COUNT+predidx[0]],
		//		csizes[combination[1]*PRED_COUNT+predidx[1]],
		//		csizes[combination[2]*PRED_COUNT+predidx[2]]
		//	);
		//}
#endif
		(void)rct_names;
		(void)pred_names;
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
		ac3_dec_init(&ec, args->streambuf+args->offsets[args->blockidx*2+0], args->streambuf+args->offsets[args->blockidx*2+1]);
		bypass_decbuf_init(&bc, args->streambuf+args->offsets[args->blockidx*2+1], args->streambuf+args->offsets[args->blockidx*2+2]);
	//	ac3_dec_init(&ec, args->decstream+(args->blockidx>0?args->decoffsets[args->blockidx-1]:0), args->decstream+args->decoffsets[args->blockidx]);
	//	bypass_decbuf_init(&bc, args->decstream+args->decoffsets[args->nblocks+args->blockidx-1], args->decstream+args->decoffsets[args->nblocks+args->blockidx]);
	//	args->dist=ac3_dec_bypass_NPOT(&ec, 2);
		entropylevel=ac3_dec_bypass_NPOT(&ec, 100);
		bestrct=ac3_dec_bypass_NPOT(&ec, RCT_COUNT);
		combination=rct_combinations[bestrct];
		predidx[0]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[1]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[2]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
#if defined ENABLE_ZERO && defined STATIC_ZERO
		for(int kc=0;kc<3;++kc)
		{
			if(predidx[kc]==PRED_ZERO)
			{
				unsigned *dsthist=(unsigned*)args->hist+ELEVELS*CLEVELS*cdfstride*kc;
				int prev=1;
				dsthist[0]=0;
				for(int ks=1;ks<256;++ks)
				{
					int nlevels=(1<<PROB_BITS)-prev;
					int sym=0;
					//if(nlevels>=0x10000)
					//{
					//	sym=ac3_dec_bypass_NPOT(&ec, nlevels>>8)<<8;
					//	nlevels=256;
					//}
					sym|=ac3_dec_bypass_NPOT(&ec, nlevels);
					dsthist[ks]=sym+prev;
					prev=dsthist[ks]+1;
				}
				dsthist[256]=1<<PROB_BITS;
			}
		}
#endif
	}
	if(entropylevel>25)		entropyidx=3;
	else if(entropylevel>17)	entropyidx=2;
	else if(entropylevel>10)	entropyidx=1;
	else				entropyidx=0;
	int shiftgain=(entropylevel*entropylevel>>4)+15;//up to 640, with blocksize 512, counters should not exceed 0xA000000
#if 0
	predidx[0]=PRED_CG;//synth2 0.8% smaller
	predidx[1]=PRED_CG;
	predidx[2]=PRED_CG;
	//predidx[0]=PRED_WG;//GDCC 0.9% larger
	//predidx[1]=PRED_WG;
	//predidx[2]=PRED_WG;
#endif
	unsigned *histptr=(unsigned*)args->hist;
#if defined ENABLE_ZERO
	//int disable_quant=0;
	//if(predidx[0]==PRED_ZERO||predidx[1]==PRED_ZERO||predidx[2]==PRED_ZERO)
	//{
	//	disable_quant=1;
	//	//predidx[0]=PRED_ZERO;
	//	//predidx[1]=PRED_ZERO;
	//	//predidx[2]=PRED_ZERO;
	//}
//#ifdef ENABLE_SELECT
//	if(predidx[0]==PRED_SELECT||predidx[1]==PRED_SELECT||predidx[2]==PRED_SELECT)
//		disable_quant=1;
//#endif
	unsigned char *CDF2sym[]=
	{
		(unsigned char*)(histptr+ELEVELS*CLEVELS*cdfstride*0+257),
		(unsigned char*)(histptr+ELEVELS*CLEVELS*cdfstride*1+257),
		(unsigned char*)(histptr+ELEVELS*CLEVELS*cdfstride*2+257),
	};
	for(int kc=0;kc<3;++kc)
	{
		unsigned *curr_hist=histptr+ELEVELS*CLEVELS*cdfstride*kc;
		if(predidx[kc]==PRED_ZERO)
		{
			if(!args->fwd)
			{
				unsigned char *curr_CDF2sym=CDF2sym[kc];
				for(int ks=0;ks<256;++ks)
				{
					int start=curr_hist[ks], end=curr_hist[ks+1];
					memset(curr_CDF2sym+start, ks, (size_t)end-start);
				}
			}
			continue;
		}
		//if(disable_quant)
		//{
		//	for(int ks=0;ks<258;++ks)
		//	{
		//		curr_hist[ks]=(unsigned)(ks*((((1LL<<PROB_BITS)-257)<<PROB_EBITS) / 257));
		//		//if(ks==257)//
		//		//	printf("");
		//	}
		//	memfill(curr_hist+258, curr_hist, sizeof(int[ELEVELS*258-258]), sizeof(int[258]));
		//}
		//else
		{
			for(int ks=0;ks<cdfstride;++ks)
				curr_hist[ks]=(unsigned)((ks*((1LL<<PROB_BITS)-args->tlevels)<<PROB_EBITS)/args->tlevels);
			memfill(curr_hist+cdfstride, curr_hist, (sizeof(int[ELEVELS*CLEVELS])-sizeof(int))*cdfstride, sizeof(int)*cdfstride);
		}
	}
	int invdist=((1<<16)-1+args->dist)/args->dist;
#else
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
#endif
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
	//int ysh=combination[3+0];
	//int ush=combination[3+1];
	//int vsh=combination[3+2];
	//__m128i shufload=_mm_set_epi8(
	//	-1, -1, -1, -1,
	//	-1, -1, -1, combination[3+2],
	//	-1, -1, -1, combination[3+1],
	//	-1, -1, -1, combination[3+0]
	//);
	//int invpermutation[3]={0};
	//invpermutation[combination[3+0]]=0*4;
	//invpermutation[combination[3+1]]=1*4;
	//invpermutation[combination[3+2]]=2*4;
	//__m128i shufstore=_mm_set_epi8(
	//	-1, -1, -1, -1,
	//	-1, -1, -1, -1,
	//	-1, -1, -1, -1,
	//	-1, invpermutation[2], invpermutation[1], invpermutation[0]
	//);
	//__m128i mhalf=_mm_set_epi32(0, 128, 128, 128);
#ifdef STATIC_ESTIMATE
	memset(args->ehist, 0, args->ehistsize);
#endif
#ifdef STATIC_ESTIMATE2
	memset(args->ehist2, 0, sizeof(args->ehist2));
	memset(args->ezerohist2, 0, sizeof(args->ezerohist2));
#endif
	int paddedblockwidth=args->x2-args->x1+16;
	memset(args->pixels, 0, args->bufsize);
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(32) short *rows[]=
		{
			args->pixels+(paddedblockwidth*((ky-0LL)&3)+8LL)*4*4,
			args->pixels+(paddedblockwidth*((ky-1LL)&3)+8LL)*4*4,
			args->pixels+(paddedblockwidth*((ky-2LL)&3)+8LL)*4*4,
			args->pixels+(paddedblockwidth*((ky-3LL)&3)+8LL)*4*4,
		};
		ALIGN(16) int yuv[4]={0};
		int helpcoeffs[]=
		{
			0,
			combination[II_HELP_U],
			combination[II_HELP_V0],

			0,
			0,
			combination[II_HELP_V1],
		};
		//int *helpers[3]=
		//{
		//	yuv+combination[6+0],
		//	yuv+combination[6+1],
		//	yuv+combination[6+2],
		//};
		int token=0, bypass=0, nbits=0;
		int pred=0, error=0, sym=0;
		unsigned char *imptr=args->imagebuf+3*(args->iw*ky+args->x1);
		//const unsigned char *srcptr=args->fwd?args->src+3*(args->iw*ky+args->x1):0;
		//unsigned char *dstptr=args->fwd?0:args->dst+3*(args->iw*ky+args->x1);
		int yidx=combination[II_PERM_Y];
		int uidx=combination[II_PERM_U];
		int vidx=combination[II_PERM_V];
		//int idx=3*(args->iw*ky+args->x1);
		//int yidx=3*(args->iw*ky+args->x1)+combination[3+0];
		//int uidx=3*(args->iw*ky+args->x1)+combination[3+1];
		//int vidx=3*(args->iw*ky+args->x1)+combination[3+2];
		ALIGN(16) short regW[16]={0};
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
		//unsigned char *effectiveCDF=(unsigned char*)(CDF+7);
		int *ctxctrs[12]={0};
		unsigned *cdfptrs[12]={0};
		const unsigned *mixinptrs[3]={0};
//#ifdef ENABLE_ZERO
//		int syms[3]={0};
//#endif
#if defined SIMD_CTX || defined SIMD_PREDS
		__m128i
			mNNW	=_mm_setzero_si128(),
			mNN	=_mm_loadu_si128((__m128i*)(rows[2]+0*4*4)),
			mNNE	=_mm_loadu_si128((__m128i*)(rows[2]+1*4*4)),
			mNWW	=_mm_setzero_si128(),
			mNW	=_mm_setzero_si128(),
			mN	=_mm_loadu_si128((__m128i*)(rows[1]+0*4*4)),
			mNE	=_mm_loadu_si128((__m128i*)(rows[1]+1*4*4)),
			mNEE	=_mm_loadu_si128((__m128i*)(rows[1]+2*4*4)),
			mWW	=_mm_setzero_si128(),
			mW	=_mm_setzero_si128();
		(void)mNNW;
		(void)mNWW;
#endif
		//for(int kx0=args->x2-args->x1;kx0--;yidx+=3, uidx+=3, vidx+=3)
		for(int kx0=args->x2-args->x1;kx0--;)
		{
			int kx=args->x2-1-kx0;
			(void)kx;
#if defined SIMD_CTX || defined SIMD_PREDS
			__m128i msum=_mm_add_epi16(mN, mW);
#endif
#ifdef SIMD_CTX
			//qe = FLOOR_LOG2(abs(N-W)+eN+eW+((eNW+eNE+eNEE+eNN+eWW)>>1)+1);  qe=CLAMP(qe, NLEVELS-1)

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
			//	*NNWW	=rows[2]-2*4*4,
				*NNW	=rows[2]-1*4*4,
				*NN	=rows[2]+0*4*4,
				*NNE	=rows[2]+1*4*4,
				*NWW	=rows[1]-2*4*4,
				*NW	=rows[1]-1*4*4,
				*N	=rows[1]+0*4*4,
				*NE	=rows[1]+1*4*4,
				*NEE	=rows[1]+2*4*4,
				*WW	=rows[0]-2*4*4,
			//	*W	=rows[0]-1*4*4,
				*curr	=rows[0]+0*4*4;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 3
#endif
			for(int kc=0;kc<3;++kc)
			{
				switch(predidx[kc])
				{
#ifdef ENABLE_ZERO
				case PRED_ZERO:
					pred=0;
					break;
#endif
#if 0
				case PRED_CLEARTYPE:
					//NW0 NW1 NW2 N0 N1 N2
					//W0  W1  W2  c0 c1 c2
					//pred(r) = SELECT(N0, W2, NW2)
					//pred(g) = SELECT(N1, c0, N0)
					//pred(b) = SELECT(N2, c1, N1)
					switch(kc)
					{
					case 0:
						pred=(abs(N[0]-NW[2])<abs(regW[2]-NW[2])?regW[2]:N[0]);
						CLAMP2(pred, -128, 127);
						break;
					case 1:pred=(abs(N[1]-N [0])<abs(regW[0]-N [0])?regW[0]:N[1]);break;
					case 2:pred=(abs(N[2]-N [1])<abs(regW[1]-N [1])?regW[1]:N[2]);break;
					}
					//if(pred<-256||pred>255||!kc&&(pred<-128||pred>127))
					//	LOG_ERROR("");
					break;
#endif
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
				case PRED_SELECT:
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
#ifdef ENABLE_CGX
				case PRED_CGX:
					MEDIAN3_32(pred, N[kc], regW[kc], regW[kc]+((N[kc]-NW[kc])>>1));
					break;
#endif
#ifdef ENABLE_CGY
				case PRED_CGY:
					MEDIAN3_32(pred, N[kc], regW[kc], N[kc]+((regW[kc]-NW[kc])>>1));
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
#ifdef ENABLE_AV6
				case PRED_AV6:
					CLAMP3_32(pred,
						regW[kc]+((6*N[kc]-5*NW[kc]-NN[kc]-WW[kc]+NE[kc])>>3),
						N[kc], regW[kc], NE[kc]
					);
					break;
#endif
#ifdef ENABLE_AV8
				case PRED_AV8:
					pred=(
						+NN[kc]+NNE[kc]
						+NW[kc]+N[kc]+NE[kc]+NEE[kc]
						+WW[kc]+regW[kc]
					)>>3;
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
#ifdef ENABLE_AV11
				case PRED_AV11:
					CLAMP3_32(pred, (
						+4*NNWW[kc]	-3*NNW[kc]	-31*NN[kc]	-38*NNE[kc]
						+7*NWW[kc]	-158*NW[kc]	+219*N[kc]	+30*NE[kc]	+19*NEE[kc]
						-42*WW[kc]	+243*regW[kc]
						+128)>>8, N[kc], regW[kc], NE[kc]
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
#ifdef ENABLE_WGB
				case PRED_WGB:
					{
						int cN=N[kc], cW=regW[kc];
						int gy=abs(cN-NN[kc])+abs(cW-NW[kc])+abs(NE[kc]-NNE[kc])+abs(NW[kc]-NNW[kc])+1;
						int gx=abs(cN-NW[kc])+abs(cW-WW[kc])+abs(NE[kc]-N[kc])+abs(NW[kc]-NWW[kc])+1;
						int sum=gx+gy;
						//int sh=(8+3-FLOOR_LOG2(sum))>>1;
						//int h1=gx>>sh;
						//int h2=gy>>sh;
						int h1=gx>>3;
						int h2=gy>>3;
						gx+=h1;
						gy+=h2;
						int pred=(gx*cN+gy*cW-(h1+h2)*NW[kc])/sum;
						int vmax=cN, vmin=cW;
						if(cN<cW)vmin=cN, vmax=cW;
						CLAMP2(pred, vmin, vmax);
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
				yuv[0]=imptr[yidx]-128;
				yuv[1]=imptr[uidx]-128;
				yuv[2]=imptr[vidx]-128;
				imptr+=3;

				//unsigned rgb=*(unsigned*)(args->src+idx);
				//rgb^=0x00808080;
				//yuv[0]=(char)(rgb>>ysh);
				//yuv[1]=(char)(rgb>>ush);
				//yuv[2]=(char)(rgb>>vsh);

				//__m128i rgb=_mm_set_epi32(0, 0, 0, *(unsigned*)(args->src+idx));
				//rgb=_mm_shuffle_epi8(rgb, shufload);
				//rgb=_mm_sub_epi32(rgb, mhalf);
				//_mm_store_si128((__m128i*)yuv, rgb);

				//yuv[0]=args->src[yidx]-128;
				//yuv[1]=args->src[uidx]-128;
				//yuv[2]=args->src[vidx]-128;
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
					offset=(helpcoeffs[kc]*yuv[0]+helpcoeffs[kc+3]*yuv[1])>>1;
					//offset=*helpers[kc];
					__m128i mhelp=_mm_set_epi32(0, 0, 0, offset);
					mp=_mm_add_epi32(mp, mhelp);
					mp=_mm_max_epi32(mp, _mm_set_epi32(0, 0, 0, -128));
					mp=_mm_min_epi32(mp, _mm_set_epi32(0, 0, 0, 127));
					mhelp=_mm_sub_epi32(mhalf, _mm_abs_epi32(mp));
					pred=_mm_extract_epi32(mp, 0);
					upred=_mm_extract_epi32(mhelp, 0);
				}
#else
				if(kc)
				{
					pred+=offset=*helpers[kc];
					CLAMP2_32(pred, pred, -128, 127);	//MOV-MOV-MOV-MIN-MAX-MOV (6 instr)
				//	CLAMP2(pred, -128, 127);		//MOV-CMP-CMOV (6 instr) slower?
				}
				int upred=128-abs(pred);
#endif
				unsigned cdf, freq;
#ifdef SIMD_CTX
				int qp=(pred+128)>>(8-CBITS_EFFECTIVE)&((1<<CBITS_EFFECTIVE)-1);
				//int qp=(pred+128)>>(8-CBITS)&((1<<CBITS)-1);
				//int qp=(pred+128)>>(8-CBITS);
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
#ifdef STATIC_ESTIMATE
				if(args->fwd)
				{
					int ectx=qe<<1|qp<<PBITS>>CBITS_EFFECTIVE;
					error=(unsigned char)(yuv[kc]-pred+128);
					++args->ehist[ectx<<8|error];
				}
#endif
#ifdef ENABLE_ZERO
				if(predidx[kc]==PRED_ZERO)
				//if(disable_quant)
				{
					unsigned *curr_hist0=cdfptrs[kc+0]=histptr+(cdfstride*ELEVELS*CLEVELS)*kc;
					if(args->fwd)
					{
						sym=(yuv[kc]-pred+128)&255;
						cdf=curr_hist0[sym];
						freq=curr_hist0[sym+1]-cdf;
						ac3_encbuf_update_N(&ec, cdf, freq, PROB_BITS);
#ifdef PRINT_CSIZES
						esizes[kc]-=log2(freq/(double)(1<<PROB_BITS));
#endif
#ifdef STATIC_ESTIMATE2
						++args->ezerohist2[kc][sym];
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
						
						unsigned c=(unsigned)(((ec.code-ec.low)<<PROB_BITS|((1LL<<PROB_BITS)-1))/ec.range);
						sym=CDF2sym[kc][c];
						cdf=curr_hist0[sym];
						freq=curr_hist0[sym+1]-cdf;
						ac3_dec_update_N(&ec, cdf, freq, PROB_BITS);
						yuv[kc]=(char)(sym-128+pred);
					}
					curr[kc+0]=regW[kc]=yuv[kc]-offset;
#if 0
					//qe=0;//o0	X  bad
					unsigned *curr_hist=cdfptrs[kc+0]=histptr+(cdfstride*ELEVELS*CLEVELS)*kc+258*qe;
					//if(curr_hist[0])//
					//	LOG_ERROR("");
	#define GETCDF0(X) ((curr_hist[X]>>PROB_EBITS)+(X))
					ctxctrs[kc+0]=&args->ctxctrs[kc][qe][0];
					if(args->fwd)
					{
						error=yuv[kc]-pred;
						{
							int negmask=error>>31;	//8 cycles
							int abserr=(error^negmask)-negmask;
							sym=error<<1^negmask;
							if(upred<abserr)//CMOV
								sym=upred+abserr;
						
							//int abserr=abs(error);
							//if(abserr<=upred)
							//	sym=error<<1^error>>31;//pack sign
							//else
							//	sym=upred+abserr;//error sign is known
						}
						
						cdf=GETCDF0(sym);
						freq=GETCDF0(sym+1)-cdf;
						ac3_encbuf_update_N(&ec, cdf, freq, PROB_BITS);
#ifdef PRINT_CSIZES
						esizes[kc]-=log2(freq/(double)(1<<PROB_BITS));
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
						
						unsigned long long r2=(ec.code-ec.low)<<PROB_BITS|((1LL<<PROB_BITS)-1);
						sym =(r2>=GETCDF0(    256)*ec.range)<<8;//256 is a valid symbol
						if(!sym)
						{
							sym =(r2>=GETCDF0(sym|128)*ec.range)<<7;
							sym|=(r2>=GETCDF0(sym| 64)*ec.range)<<6;
							sym|=(r2>=GETCDF0(sym| 32)*ec.range)<<5;
							sym|=(r2>=GETCDF0(sym| 16)*ec.range)<<4;
							sym|=(r2>=GETCDF0(sym|  8)*ec.range)<<3;
							sym|=(r2>=GETCDF0(sym|  4)*ec.range)<<2;
							sym|=(r2>=GETCDF0(sym|  2)*ec.range)<<1;
							sym|= r2>=GETCDF0(sym|  1)*ec.range;
						}
						//if((unsigned)sym>256)
						//	LOG_ERROR("");
						cdf=GETCDF0(sym);
						freq=GETCDF0(sym+1)-cdf;
						ac3_dec_update_N(&ec, cdf, freq, PROB_BITS);
						{
							int negmask=pred>>31;	//11 cycles
							int e2=upred-sym;
							error=sym>>1^-(sym&1);
							e2=(e2^negmask)-negmask;
							if((upred<<1)<sym)//CMOV
								error=e2;
						}

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
					syms[kc]=sym;
#endif
					continue;
				}
#endif
				unsigned *curr_hist0=cdfptrs[kc+0]=histptr+cdfstride*(ELEVELS*(CLEVELS*kc+qp)+qe);
				ctxctrs[kc+0]=&args->ctxctrs[kc][qe][qp];
				int qe2=-(qe<ELEVELS-1);
				unsigned *curr_hist1=cdfptrs[kc+3]=curr_hist0+(cdfstride&qe2);
				qe2=qe-qe2;
				ctxctrs[kc+3]=&args->ctxctrs[kc][qe2][qp];

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

				if(args->fwd)
				{
				//again:
					if(args->dist>1)
					{
						error=(yuv[kc]-pred);
						//naive deadzone
						error=(error*invdist>>16)+((unsigned)error>>31);
					//	error=(yuv[kc]-pred)/args->dist;
						sym=error<<1^error>>31;
						yuv[kc]=error*args->dist+pred;
						CLAMP2(yuv[kc], -128, 127);
						error=yuv[kc]-pred;
					}
					else
					{
						error=yuv[kc]-pred;
						int negmask=error>>31;	//8 cycles
						int abserr=(error^negmask)-negmask;
						sym=error<<1^negmask;
						if(upred<abserr)//CMOV
							sym=upred+abserr;

						//int abserr=abs(error);
						//sym=upred+abserr;
						//if(abserr<=upred)//CMOV
						//	sym=error<<1^error>>31;
						
						//int abserr=abs(error);
						//if(abserr<=upred)
						//	sym=error<<1^error>>31;//pack sign
						//else
						//	sym=upred+abserr;//error sign is known
					}
					{
						token=sym;
						if(sym>=(1<<CONFIG_EXP))
						{
							//nbits=FLOOR_LOG2((unsigned)sym)-(CONFIG_MSB+CONFIG_LSB);
							nbits=(31-(CONFIG_MSB+CONFIG_LSB))-_lzcnt_u32(sym);
							token = (1<<CONFIG_EXP)-((CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB))<<(CONFIG_MSB+CONFIG_LSB)) + (
									nbits<<(CONFIG_MSB+CONFIG_LSB)|
									(sym>>nbits&((1<<CONFIG_MSB)-1)<<CONFIG_LSB)|
									(sym&((1<<CONFIG_LSB)-1))
								);
							bypass=sym>>CONFIG_LSB&((1LL<<nbits)-1);
							//bypass=_bextr_u32(sym>>CONFIG_LSB, 0, nbits);	//slow
							bypass_encbuf(&bc, bypass, nbits);
#ifdef PRINT_CSIZES
							esizes[kc]+=nbits;
#endif
#ifdef STATIC_ESTIMATE2
							args->esize2+=nbits;
#endif
						}
#ifdef _DEBUG
						if(token>=args->tlevels)
							LOG_ERROR("YXC %d %d %d  token %d/%d", ky, kx, kc, token, args->tlevels);
#endif
						cdf=GETCDF(token);
						freq=GETCDF(token+1)-cdf;
						ac3_encbuf_update_N(&ec, cdf, freq, PROB_BITS);
#ifdef PRINT_CSIZES
						esizes[kc]-=log2(freq/(double)(1<<PROB_BITS));
#endif
#ifdef _DEBUG
						if(ec.dstptr>=ec.dstend)
							LOG_ERROR("Encoder out of memory at YXC %d %d %d  inflation %8.4lf%%", ky, kx, kc,
								100.*3*((args->x2-args->x1)*ky+kx)/((args->x2-args->x1)*(args->y2-args->y1))
							);
#endif
#ifdef STATIC_ESTIMATE2
						{
							++args->ehist2[kc][qe][qp<<PBITS>>CBITS_EFFECTIVE][token];
						}
#endif
					}
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
					if(entropyidx&&args->dist<4)
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

						//__m256i mtoken=_mm256_castsi128_si256(_mm_broadcastd_epi32(_mm_set_epi32(0, 0, 0, token)));//X  need freq & CDF
						//mc0=_mm256_permutevar8x32_epi32(mx0, mtoken);
						//mc1=_mm256_permutevar8x32_epi32(mx1, mtoken);
						//mc2=_mm256_permutevar8x32_epi32(mx2, mtoken);
						//mc3=_mm256_permutevar8x32_epi32(mx3, mtoken);
						//mtoken=_mm256_srli_epi32(mtoken, 3);
						//mc0=_mm256_blend_epi32(mc0, mc1, 0xFE);
						//mc2=_mm256_blend_epi32(mc2, mc3, 0xF8);
						//mc0=_mm256_blend_epi32(mc0, mc2, 0xFC);
						//mc0=_mm256_permutevar8x32_epi32(mc0, mtoken);
						//_mm256_extract_epi32(mc0, 0);
						
						cdf=CDF[token+7];//5%
						freq=CDF[token+8]-cdf;
						
						//unsigned long long reg=*(unsigned long long*)(effectiveCDF+4LL*token);//12%
						//cdf=(unsigned)reg;
						//freq=(reg>>32)-cdf;
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
					if(args->dist>=4)
					{
						//naive deadzone
						error=(sym>>1^-(sym&1))*args->dist;
						yuv[kc]=error+pred;
						CLAMP2(yuv[kc], -128, 127);
						error=yuv[kc]-pred;
					}
					else
					{
						int negmask=pred>>31;	//11 cycles
						int e2=upred-sym;
						error=sym>>1^-(sym&1);
						e2=(e2^negmask)-negmask;
						if((upred<<1)<sym)//CMOV
							error=e2;

						//error=upred-sym;
						//if(pred<0)
						//	error=sym-upred;
						//if(sym<=(upred<<1))//CMOV
						//	error=sym>>1^-(sym&1);

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
				}
				curr[kc+0]=regW[kc]=yuv[kc]-offset;
			//	curr[kc+4]=regW[kc+4]=(2*rows[0][kc-1*4*4+4]+((error<<1^error>>31)<<6)+MAXVAR(rows[1][kc+2*4*4+4], rows[1][kc+3*4*4+4]))>>2;
			//	curr[kc+4]=regW[kc+4]=(2*rows[0][kc-1*4*4+4]+(abs(error)<<6)+MAXVAR(rows[1][kc+2*4*4+4], rows[1][kc+3*4*4+4]))>>2;
			//	curr[kc+4]=regW[kc+4]=(2*rows[0][kc-1*4*4+4]+(token<<6)+MAXVAR(rows[1][kc+2*4*4+4], rows[1][kc+3*4*4+4]))>>2;
				curr[kc+4]=regW[kc+4]=abs(error);
			//	regW[kc+8]=error;
				mixinptrs[kc]=args->mixincdfs+args->tlevels*token;
			}
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 3
#endif
			for(int kc=0;kc<3;++kc)
			{
#ifdef ENABLE_ZERO
				//if(disable_quant)
				if(predidx[kc]==PRED_ZERO)
				{
#if 0
					int sh=*ctxctrs[kc+0]+=shiftgain;
					sh=FLOOR_LOG2(sh)>>1;
					unsigned *curr_cdf=cdfptrs[kc];
					token=syms[kc];
					//if(curr_cdf[0])//
					//	LOG_ERROR("");
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
					for(int ks=1;ks<257;++ks)
						curr_cdf[ks]+=(((((1LL<<PROB_BITS)-257)<<PROB_EBITS)&-(ks>token))-curr_cdf[ks])>>sh;
					//	curr_cdf[ks]+=(((((1LL<<PROB_BITS)-257)<<PROB_EBITS)&-(ks>token))-curr_cdf[ks])>>10;
					//	curr_cdf[ks]+=(((((1LL<<PROB_BITS)-257)<<PROB_EBITS)&-(ks>token))-curr_cdf[ks])>>7;
#endif
					continue;
				}
#endif
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
				imptr[yidx]=yuv[0]+128;
				imptr[uidx]=yuv[1]+128;
				imptr[vidx]=yuv[2]+128;
				imptr+=3;

				//unsigned *ptr=(unsigned*)(args->dst+idx);		//X  UB
				//int rgb=(yuv[0]&255)<<ysh|(yuv[1]&255)<<ush|(yuv[2]&255)<<vsh;
				//unsigned v0=*ptr&0xFF00000;
				//rgb^=0x808080;
				//*ptr=v0|rgb;

				//unsigned *ptr=(unsigned*)(args->dst+idx);		//X  UB
				//__m128i rgb=_mm_load_si128((__m128i*)yuv);
				//unsigned val=*ptr;
				//rgb=_mm_add_epi32(rgb, mhalf);
				//val&=0xFF000000;
				//rgb=_mm_shuffle_epi8(rgb, shufstore);
				//*ptr=val|_mm_extract_epi32(rgb, 0);

				//args->dst[yidx]=yuv[0]+128;
				//args->dst[uidx]=yuv[1]+128;
				//args->dst[vidx]=yuv[2]+128;
				if(args->dist>=4)
					guide_update_sqe(args->imagebuf, kx, ky);
				else
					guide_check(args->imagebuf, kx, ky);
			}
			rows[0]+=4*4;
			rows[1]+=4*4;
			rows[2]+=4*4;
			rows[3]+=4*4;
#if defined SIMD_CTX || defined SIMD_PREDS
			mNNW	=mNN;
			mNN	=mNNE;
			mNNE	=_mm_loadu_si128((__m128i*)(rows[2]+1*4*4));
			mNWW	=mNW;
			mNW	=mN;
			mN	=mNE;
			mNE	=mNEE;
			mNEE	=_mm_loadu_si128((__m128i*)(rows[1]+2*4*4));
			mWW	=mW;
			mW	=_mm_load_si128((__m128i*)regW);
		//	mW	=_mm_loadu_si128((__m128i*)(rows[0]-1*4*4));
#endif
		}
	}
	if(args->fwd)
	{
		ac3_encbuf_flush(&ec);
		args->chunksizes[args->blockidx*2+0]=(int)(ec.dstptr-ec.dststart);
		//args->chunksizes[args->blockidx*2+0]=(int)(ptr-(args->streambuf+args->offsets[args->blockidx*2+0]));
		bypass_encbuf_flush(&bc);
		args->chunksizes[args->blockidx*2+1]=(int)(bc.dstptr-bc.dststart);
		//args->chunksizes[args->blockidx*2+1]=(int)(ptr-(args->streambuf+args->offsets[args->blockidx*2+1]));
#ifdef PRINT_CSIZES
		printf(" Act %11.2lf %11.2lf %11.2lf\n",
			esizes[0]/8,
			esizes[1]/8,
			esizes[2]/8
		);
		args->csizes[0]+=esizes[0];
		args->csizes[1]+=esizes[1];
		args->csizes[2]+=esizes[2];
		//if(!args->blockidx)
		//{
		//	printf("%11.2lf Y act.\n", args->csizes[0]);
		//	printf("%11.2lf U act.\n", args->csizes[1]);
		//	printf("%11.2lf V act.\n", args->csizes[2]);
		//	printf("%8d cblock\n%8d Tokens\n%8d Bypass\n",
		//		args->enctokenoffsets[currblockidx]+args->encbypassoffsets[currblockidx]+8,
		//		args->enctokenoffsets[currblockidx],
		//		args->encbypassoffsets[currblockidx]
		//	);
		//}
#endif
#ifdef STATIC_ESTIMATE
		for(int kctx=0;kctx<3*(ELEVELS<<PBITS);++kctx)
		{
			int *curr_hist=args->ehist+((ptrdiff_t)kctx<<8);
			int sum=0;
			for(int ks=0;ks<256;++ks)
				sum+=curr_hist[ks];
			if(!sum)
				continue;
			double e=0, norm=1./sum;
			for(int ks=0;ks<256;++ks)//calc entropy
			{
				int freq=curr_hist[ks];
				if(freq)
					e-=freq*log2(freq*norm);
			}
			args->esizes[kctx<<1|0]+=e/8;
			
			double hsize=0;
			int cdfW=0;
			int sum2=0;
			for(int ks=0;ks<256;++ks)//calc model stat overhead
			{
				int sym=(unsigned char)((ks>>1^-(ks&1))+128);
				int freq=curr_hist[sym];
				int cdf=sum2*((1ULL<<16)-256)/sum+ks;
				int csym=cdf-cdfW;
				hsize+=log2(csym+1);
				cdfW=cdf;
				sum2+=freq;
			}
			args->esizes[kctx<<1|1]+=(hsize+256-1)/8;
		}
#endif
#ifdef STATIC_ESTIMATE2
		for(int kc=0;kc<3;++kc)
		{
			int sum=0;
			for(int ks=0;ks<256;++ks)
				sum+=args->ezerohist2[kc][ks];
			if(sum)
			{
				double e=0, norm=1./sum;
				for(int ks=0;ks<256;++ks)
				{
					int freq=args->ezerohist2[kc][ks];
					if(freq)
						e-=freq*log2(freq*norm);
				}
				//if(e!=e)
				//	LOG_ERROR("sum %d", sum);
				args->esize2+=e/8;
			}
			for(int ke=0;ke<ELEVELS;++ke)
			{
				for(int kp=0;kp<PLEVELS;++kp)
				{
					int *curr_hist=args->ehist2[kc][ke][kp];
					sum=0;
					for(int ks=0;ks<33;++ks)
						sum+=curr_hist[ks];
					if(!sum)
						continue;
					double e=0, norm=1./sum;
					for(int ks=0;ks<33;++ks)//calc entropy
					{
						int freq=curr_hist[ks];
						if(freq)
							e-=freq*log2(freq*norm);
					}
					//if(e!=e)
					//	LOG_ERROR("sum %d", sum);
					args->esize2+=e/8;
					double hsize=0;
					int cdfW=0;
					int sum2=0;
					for(int ks=0;ks<33;++ks)//calc model stat overhead
					{
						int freq=curr_hist[ks];
						int cdf=sum2*((1ULL<<16)-33)/sum+ks;
						int csym=cdf-cdfW;
						//double LOL_1=log2(csym+1);
						//if(LOL_1!=LOL_1)//
						//	LOG_ERROR("edible thermonuclear");
						hsize+=log2(csym+1);
						cdfW=cdf;
						sum2+=freq;
					}
					args->hsize2+=(hsize+33-1)/8;
				}
			}
		}
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
int c24_codec(int argc, char **argv)
{
#ifdef LOUD
	double t=time_sec();
#endif
	if(argc!=2&&argc!=3&&argc!=4&&argc!=5)
	{
		printf(
			"Usage: \"%s\"  input  output  [maxthreads]  [dist]    Encode/decode.\n"
			"       \"%s\"  input                                  Test without saving.\n"
			"[maxthreads]:\n"
			"  0: nthreads = number of cores (default)\n"
			"  1: Single thread\n"
			, argv[0]
			, argv[0]
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argc<3?0:argv[2];
	int maxthreads=argc<4?0:atoi(argv[3]), dist=argc<5?1:atoi(argv[4]);
	if(dist>1)
		CLAMP2(dist, 4, 64);
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	int fwd=0, iw=0, ih=0;
	ptrdiff_t srcsize=0, streamsize=0, usize=0;
	unsigned char *imagebuf=0, *streambuf=0;
	unsigned char *streamptr=0, *streamstart=0, *streamend=0;
	(void)streamstart;
	{
		struct stat info={0};
		int error=stat(srcfn, &info);
		if(error)
		{
			LOG_ERROR("Cannot stat \"%s\"", srcfn);
			return 1;
		}
		srcsize=info.st_size;
	}
	{
		ptrdiff_t nread;
		int tag=0;
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		fread(&tag, 1, 2, fsrc);
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('2'|'4'<<8))
		{
			LOG_ERROR("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
		{
			int c;
#ifdef LOUD
			print_timestamp("%Y-%m-%d_%H%M%S\n");
#endif
			c=fgetc(fsrc);
			if(c!='\n')
			{
				LOG_ERROR("Invalid PPM file");
				return 1;
			}
			c=fgetc(fsrc);
			while(c=='#')
			{
				c=fgetc(fsrc);
				while(c!='\n')
					c=fgetc(fsrc);
				c=fgetc(fsrc);
			}
			iw=0;
			while((unsigned)(c-'0')<10)
			{
				iw=10*iw+c-'0';
				c=fgetc(fsrc);
			}
			while(c<=' ')
				c=fgetc(fsrc);
			ih=0;
			while((unsigned)(c-'0')<10)
			{
				ih=10*ih+c-'0';
				c=fgetc(fsrc);
			}
			while(c<=' ')
				c=fgetc(fsrc);
			while(c=='#')
			{
				c=fgetc(fsrc);
				while(c!='\n')
					c=fgetc(fsrc);
				c=fgetc(fsrc);
			}
			c=c<<8|fgetc(fsrc);
			c=c<<8|fgetc(fsrc);
			c=c<<8|fgetc(fsrc);
			if(c!=('2'<<24|'5'<<16|'5'<<8|'\n'))
			{
				LOG_ERROR("Unsupported PPM file");
				return 1;
			}
		}
		else
		{
			fread(&iw, 1, 4, fsrc);
			fread(&ih, 1, 4, fsrc);
			nread=ftell(fsrc);
			streamsize=srcsize-nread;

		}
		usize=(ptrdiff_t)3*iw*ih;
		imagebuf=(unsigned char*)malloc(usize);
		if(!imagebuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		ptrdiff_t expected=0;
		if(fwd)
		{
			expected=usize;
			nread=fread(imagebuf, 1, usize, fsrc);
			guide_save(imagebuf, iw, ih);
		}
		else
		{
			streambuf=(unsigned char*)malloc(streamsize);
			if(!streambuf)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			expected=streamsize;
			nread=fread(streambuf, 1, streamsize, fsrc);
		}
		if(nread!=expected)
			printf("Truncated  expected %td  read %td", expected, nread);
		fclose(fsrc);
	}
	if(!fwd)
	{
		streamstart=streambuf;
		streamend=streambuf+streamsize;
		streamptr=streambuf;

		dist=*streamptr++;
	}
	int xblocks=(iw+BLOCKX-1)/BLOCKX;
	int yblocks=(ih+BLOCKY-1)/BLOCKY;
	int nblocks=xblocks*yblocks;
	int ncores=query_cpu_cores();
	int nthreads=MINVAR(nblocks, ncores);
	if(maxthreads&&nthreads>maxthreads)
		nthreads=maxthreads;

	int argssize=nthreads*sizeof(ThreadArgs);
	ThreadArgs *args=(ThreadArgs*)malloc(argssize);
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(args, 0, argssize);

	int *blocksizes=0;
	if(fwd)
	{
		blocksizes=(int*)malloc(nblocks*(int)sizeof(int[2]));
		if(!blocksizes)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(blocksizes, 0, nblocks*sizeof(int[2]));
	}
	else
	{
		blocksizes=(int*)streamptr;
		streamptr+=nblocks*sizeof(int[2]);
	}
	int offsetssize=(2*nblocks+1)*(int)sizeof(int);
	int *offsets=(int*)malloc(offsetssize);
	if(!offsets)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	if(fwd)
	{
		offsets[0]=0;
		int asum=0;
		for(int kb2=0;kb2<nblocks*2;++kb2)
		{
			//if(kb==8*2)//
			//	printf("");

			int kb=kb2>>1;
			int bx=kb%xblocks, by=kb/xblocks;
			int x1=bx*BLOCKX, y1=by*BLOCKY;
			int x2=x1+BLOCKX, y2=y1+BLOCKY;
			if(x2>iw)
				x2=iw;
			if(y2>ih)
				y2=ih;
			int area=(x2-x1)*(y2-y1)*4;
			asum+=area;
			//if(asum<0)
			//	LOG_ERROR("");
			offsets[kb2+1]=asum;
		}
		streamsize=asum;
		streambuf=(unsigned char*)malloc(streamsize);
		if(!streambuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		streamstart=streambuf;
		streamend=streambuf+streamsize;
		streamptr=streambuf;
	}
	else
	{
		int sum=0;
		offsets[0]=0;
		for(int kb=0;kb<nblocks*2;++kb)
		{
			sum+=blocksizes[kb];
			offsets[kb+1]=sum;
		}
		if(streamptr+sum!=streamend)
			LOG_ERROR("Corrupt file  expected %d  got %d bytes", sum, (int)(streamend-streamptr));
	}
	unsigned *mixincdfs=0;
	int tlevels, histsize, statssize;
	int maxwidth=iw;
	if(maxwidth>BLOCKX)
		maxwidth=BLOCKX;
	{
		int token=0, bypass=0, nbits=0;

		quantize_pixel(256, &token, &bypass, &nbits);//256 is a valid symbol due to CALIC-like sign packing
		tlevels=token+1;
		statssize=(tlevels+1)*(int)sizeof(int[3*ELEVELS*CLEVELS]);//CDF padding, contains (1<<PROB_BITS)
		histsize=(int)sizeof(int[OCH_COUNT*PRED_COUNT<<8]);
		if(histsize<statssize)
			histsize=statssize;
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
		arg->imagebuf=imagebuf;
		arg->streambuf=streamptr;
		arg->iw=iw;
		arg->ih=ih;
		arg->offsets=offsets;
		arg->chunksizes=blocksizes;
		
		arg->b1=kb;
		kb=(k+1)*nblocks/nthreads;
		arg->b2=kb;
		arg->xblocks=xblocks;
		arg->nblocks=nblocks;

		arg->bufsize=sizeof(short[4*4*4])*(maxwidth+16LL);//4 padded rows * 4 channels max * {pixels, wg_errors}
		arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));
		arg->histsize=histsize;
		arg->hist=(int*)malloc(histsize);
		if(!arg->pixels||!arg->hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		arg->dist=dist;
		
		arg->tlevels=tlevels;
		arg->fwd=fwd;
		arg->mixincdfs=mixincdfs;
#ifdef ENABLE_MT
		arg->loud=0;
#else
		arg->loud=nblocks<MAXPRINTEDBLOCKS;
#endif
	}
#ifdef ENABLE_MT
	if(nthreads>1)
	{
		void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads);
		mt_finish(ctx);
	}
	else
#endif
	{
		for(int k=0;k<nthreads;++k)
			block_thread(args+k);
	}
	FILE *fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
	ptrdiff_t csize=0;
	if(fwd)
	{
		csize+=fwrite("24", 1, 2, fdst);
		csize+=fwrite(&iw, 1, 4, fdst);
		csize+=fwrite(&ih, 1, 4, fdst);
		csize+=fwrite(&dist, 1, 1, fdst);
		csize+=fwrite(blocksizes, 1, nblocks*sizeof(int[2]), fdst);
		for(int kb=0;kb<nblocks*2;++kb)
			csize+=fwrite(streambuf+offsets[kb], 1, blocksizes[kb], fdst);
	}
	else
	{
		csize=srcsize;
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		fwrite(imagebuf, 1, usize, fdst);
	}
	fclose(fdst);
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		_mm_free(arg->pixels);
		free(arg->hist);
	}
	free(args);
	free(imagebuf);
	free(streambuf);
	free(mixincdfs);
	free(offsets);
	if(fwd)
		free(blocksizes);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
		//printf("Mem usage: ");
		//print_size((double)memusage, 8, 4, 0, 0);
		//printf("\n");
	}
	printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t, usize/(t*1024*1024));
#ifdef ENABLE_GUIDE
	if(!fwd&&dist>=4)
	{
		double rmse=sqrt(g_sqe/usize), psnr=20*log10(255/rmse);
		printf("RMSE %12.6lf  PSNR %12.6lf\n", rmse, psnr);
	}
#endif
#else
	(void)csize;
#endif
	return 0;
}