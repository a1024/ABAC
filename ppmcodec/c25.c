#include"codec.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
	#define ENABLE_MT

	#define ENABLE_EDGECASES

	#define ENABLE_AV4	//good with GDCC
	#define ENABLE_AV9	//good
	#define ENABLE_WG	//good with noisy areas

#define ANALYSIS_XSTRIDE 2
#define ANALYSIS_YSTRIDE 2

//	#define AC_VALIDATE
//	#define AC_IMPLEMENTATION
#include"entropy.h"

#define BLOCKSIZE 512
#define MAXPRINTEDBLOCKS 500
#define CLEVELS0 1
#define CLEVELS1 256
#define CLEVELS2 17
#define CLEVELS3 12	//squared
#define CLEVELS4 19
#define CLEVELS5 19
#define CLEVELS6 256
#define CLEVELS7 512
#define NCTX (CLEVELS0+CLEVELS1+CLEVELS2+CLEVELS3*CLEVELS3+CLEVELS4+CLEVELS5+CLEVELS6+CLEVELS7)
#define CTX_COUNT 8

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
static const char *och_names[OCH_COUNT]=
{
#define OCH(LABEL) #LABEL,
	OCHLIST
#undef  OCH
};

#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		3, 3, 3)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		3, 3, 1)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		3, 3, 0)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		3, 3, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		3, 3, 1)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		3, 3, 1)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		3, 3, 0)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		3, 0, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		3, 0, 1)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		3, 0, 1)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		3, 0, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		3, 0, 1)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		3, 0, 1)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		3, 0, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		3, 0, 1)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		3, 0, 1)
typedef enum _RCTIndex
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF, UOFF, VOFF) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][6]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF, UOFF, VOFF) {YIDX, UIDX, VIDX, YOFF, UOFF, VOFF},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF, UOFF, VOFF) #LABEL,
	RCTLIST
#undef  RCT
};

#define PREDLIST\
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
static const char *pred_names[PRED_COUNT]=
{
#define PRED(LABEL) #LABEL,
	PREDLIST
#undef  PRED
};
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
		DST=_mm256_and_si256(_mm256_cvtps_epi32(_a0lo), _mm256_set1_epi32(0xFFFF));\
		DST=_mm256_or_si256(DST, _mm256_slli_epi32(_mm256_cvtps_epi32(_a0hi), 16));\
	}while(0)

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
#define CONFIG_EXP 4
#define CONFIG_MSB 1
#define CONFIG_LSB 0
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

	int fwd, test, loud, x1, x2, y1, y2;
	int bufsize, histsize;
	short *pixels;
	int *hist;

	BList list;
	const unsigned char *decstart, *decend;

	int tlevels;

	//aux
	int blockidx;
	double bestsize;
	int bestrct, predidx[3];
} ThreadArgs;
static void block_thread(void *param)
{
	const int half=128;
	ThreadArgs *args=(ThreadArgs*)param;
	AC3 ec;
	const unsigned char *image=args->fwd?args->src:args->dst;
	unsigned char bestrct=0, predidx[3]={0};
	const unsigned char *combination=0;
	int cdfstride=args->tlevels+1;
	int entropyidx=0;

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
			predidx[0]=PRED_CG;
			predidx[1]=PRED_CG;
			predidx[2]=PRED_CG;
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
		pred=_mm256_slli_epi16(pred, 8);\
		pred=_mm256_srai_epi16(pred, 8);\
		pred=_mm256_sub_epi16(pred, amin);\
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
				//CG = median(N, W, N+W-NW)
				vmin[0]=_mm256_min_epi16(N, W);
				vmax[0]=_mm256_max_epi16(N, W);
				pred=_mm256_sub_epi16(_mm256_add_epi16(N, W), NW);
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
				vmin[1]=_mm256_min_epi16(N3, W3);
				vmax[1]=_mm256_max_epi16(N3, W3);
				pred=_mm256_sub_epi16(_mm256_add_epi16(N3, W3), NW3);
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
				vmin[2]=_mm256_min_epi16(N4, W4);
				vmax[2]=_mm256_max_epi16(N4, W4);
				pred=_mm256_sub_epi16(_mm256_add_epi16(N4, W4), NW4);
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

				vmin[0]=_mm256_min_epi16(vmin[0], NE);
				vmax[0]=_mm256_max_epi16(vmax[0], NE);
				vmin[1]=_mm256_min_epi16(vmin[1], NE3);
				vmax[1]=_mm256_max_epi16(vmax[1], NE3);
				vmin[2]=_mm256_min_epi16(vmin[2], NE4);
				vmax[2]=_mm256_max_epi16(vmax[2], NE4);
				
				//AV4 = clamp((4*(N+W)+NE-NW)>>3, N,W,NE)
#ifdef ENABLE_AV4
				pred=_mm256_add_epi16(N, W);
				pred=_mm256_slli_epi16(pred, 2);
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE, NW));
				pred=_mm256_srai_epi16(pred, 3);
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
				
				//WG
				//gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1
				//gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1
				//pred=(gx*N+gy*W)/(gx+gy)
#ifdef ENABLE_WG
				__m256i gx, gy;
				gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W, WW)), 1);
				gy=_mm256_abs_epi16(_mm256_sub_epi16(W, NW));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N, NW)));
				gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N, NN)), 1));
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
				gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W3, WW3)), 1);
				gy=_mm256_abs_epi16(_mm256_sub_epi16(W3, NW3));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N3, NW3)));
				gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N3, NN3)), 1));
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
				gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W4, WW4)), 1);
				gy=_mm256_abs_epi16(_mm256_sub_epi16(W4, NW4));
				gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N4, NW4)));
				gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N4, NN4)), 1));
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
		args->bestsize=bestsize;
		entropyidx=(int)((
			csizes[combination[0]*PRED_COUNT+predidx[0]]+
			csizes[combination[1]*PRED_COUNT+predidx[1]]+
			csizes[combination[2]*PRED_COUNT+predidx[2]]
		)*gain*(100/3.));//percent
		if(entropyidx>25)
			entropyidx=3;
		else if(entropyidx>17)
			entropyidx=2;
		else if(entropyidx>10)
			entropyidx=1;
		else
			entropyidx=0;
	skip_analysis:
		args->bestrct=bestrct;
		args->predidx[0]=predidx[0];
		args->predidx[1]=predidx[1];
		args->predidx[2]=predidx[2];
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
		blist_init(&args->list);
		ac3_enc_init(&ec, &args->list);
		ac3_enc_bypass_NPOT(&ec, bestrct, RCT_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[0], PRED_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[1], PRED_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[2], PRED_COUNT);
		ac3_enc_bypass(&ec, entropyidx, 2);
	}
	else
	{
		ac3_dec_init(&ec, args->decstart, args->decend);
		bestrct=ac3_dec_bypass_NPOT(&ec, RCT_COUNT);
		predidx[0]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[1]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[2]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		combination=rct_combinations[bestrct];
		entropyidx=ac3_dec_bypass(&ec, 2);
	}
	int yidx=0, uidx=1, vidx=2;
	switch(bestrct)
	{
	case RCT_R_G_B:
	case RCT_R_G_BG:
	case RCT_R_G_BR:
	case RCT_R_GR_BR:
	case RCT_R_GR_BG:
	//case RCT_R_G_B2:
	//case RCT_R_GR_B2:
		yidx=0;
		uidx=1;
		vidx=2;
		break;
	case RCT_G_B_RG:
	case RCT_G_B_RB:
	case RCT_G_BG_RG:
	case RCT_G_BG_RB:
	//case RCT_G_B_R2:
	//case RCT_G_BG_R2:
		yidx=1;
		uidx=2;
		vidx=0;
		break;
	case RCT_B_R_GR:
	case RCT_B_R_GB:
	case RCT_B_RB_GB:
	case RCT_B_RB_GR:
	//case RCT_B_RB_G2:
		yidx=2;
		uidx=0;
		vidx=1;
		break;
	case RCT_G_RG_BR:
	//case RCT_G_RG_B2:
		yidx=1;
		uidx=0;
		vidx=2;
		break;
	case RCT_B_GB_RG:
	//case RCT_B_GB_R2:
		yidx=2;
		uidx=1;
		vidx=0;
		break;
	case RCT_R_BR_GB:
	//case RCT_R_B_G2:
	//case RCT_R_BR_G2:
		yidx=0;
		uidx=2;
		vidx=1;
		break;
	}
	int histlimit=6144<<(3-entropyidx);
	unsigned short *histptr=(unsigned short*)args->hist;
	//{
	//	static const int init_freqs[]={32, 8, 6, 4, 3, 2, 1};
	//	int sum=0;
	//	for(int ks=0;ks<args->tlevels;++ks)
	//	{
	//		int freq=init_freqs[MINVAR(ks, (int)_countof(init_freqs)-1)];
	//		sum+=histptr[ks]=freq;
	//	}
	//	histptr[args->tlevels]=sum;
	//}
	//memfill(histptr+cdfstride, histptr, args->histsize-cdfstride*sizeof(short), cdfstride*sizeof(short));
	memset(histptr, 0, args->histsize);
	memset(args->pixels, 0, args->bufsize);
	int paddedblockwidth=args->x2-args->x1+16;
	ptrdiff_t mixer[3][CTX_COUNT]={0};
	FILLMEM((ptrdiff_t*)mixer, 0x10000/CTX_COUNT, sizeof(mixer), sizeof(ptrdiff_t));
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) short *rows[]=
		{
			args->pixels+(paddedblockwidth*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+(paddedblockwidth*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+(paddedblockwidth*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+(paddedblockwidth*((ky-3LL)&3)+8LL)*4*2,
		};
		int yuv[4]={0};
		int token=0, bypass=0, nbits=0;
		int pred=0, error=0, sym=0;
		ALIGN(16) int grads[12]={0};
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int idx=3*(args->iw*ky+kx);
			short
			//	*NNN	=rows[3]+0*4*2,
				*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NWW	=rows[1]-2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*NEEE	=rows[1]+3*4*2,
			//	*WWW	=rows[0]-3*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
#ifdef ENABLE_EDGECASES
			//if(ky<=args->y1+2)
			//{
				if(ky<=args->y1+1)
				{
					if(ky==args->y1)
						NEEE=NE=NWW=NW=N=W;
					NNW=NW;
					NN=N;
					NNE=NE;
				}
			//	NNN=NN;
			//}
			//if(kx<=args->x1+2)
			//{
				if(kx<=args->x1+1)
				{
					if(kx<=args->x1)
						NWW=NW=W=N;
					WW=W;
				}
			//	WWW=WW;
			//}
			if(kx>=args->x2-3)
			{
				if(kx>=args->x2-1)
				{
					NNE=NN;
					NE=N;
				}
				NEEE=NE;
			}
#endif
			if(args->fwd)
			{
				yuv[0]=args->src[idx+yidx]-128;
				yuv[1]=args->src[idx+uidx]-128;
				yuv[2]=args->src[idx+vidx]-128;
			}
		//	__m128i mNNW	=_mm_loadu_si128((__m128i*)NNW);
			__m128i mNN	=_mm_loadu_si128((__m128i*)NN);
			__m128i mNNE	=_mm_loadu_si128((__m128i*)NNE);
		//	__m128i mNWW	=_mm_loadu_si128((__m128i*)NWW);
			__m128i mNW	=_mm_loadu_si128((__m128i*)NW);
			__m128i mN	=_mm_loadu_si128((__m128i*)N);
			__m128i mNE	=_mm_loadu_si128((__m128i*)NE);
			__m128i mWW	=_mm_loadu_si128((__m128i*)WW);
			__m128i mW	=_mm_loadu_si128((__m128i*)W);
			__m128i mdiff=_mm_sub_epi16(mN, mW);
			__m128i mr=_mm_abs_epi16(mdiff);
			__m128i mx=_mm_add_epi16(_mm_abs_epi16(_mm_sub_epi16(mW, mWW)), mr);
			__m128i my=_mm_add_epi16(_mm_abs_epi16(_mm_sub_epi16(mN, mNN)), mr);
			mr=_mm_mullo_epi16(mr, mr);
			mx=_mm_add_epi16(mx, _mm_abs_epi16(_mm_sub_epi16(mN, mNW)));
			my=_mm_add_epi16(my, _mm_abs_epi16(_mm_sub_epi16(mW, mNW)));
			mx=_mm_add_epi16(mx, _mm_abs_epi16(_mm_sub_epi16(mNE, mN)));
			my=_mm_add_epi16(my, _mm_abs_epi16(_mm_sub_epi16(mNE, mNNE)));
			__m128i meW=_mm_max_epi16(mW, mWW);
			__m128i meN=_mm_max_epi16(mN, mNN);
			meW=_mm_add_epi16(meW, mW);
			meN=_mm_add_epi16(meN, mN);
			meW=_mm_shuffle_epi32(meW, _MM_SHUFFLE(1, 0, 3, 2));
			meN=_mm_shuffle_epi32(meN, _MM_SHUFFLE(1, 0, 3, 2));
			meW=_mm_slli_epi16(meW, 1);
			meN=_mm_slli_epi16(meN, 1);
			mx=_mm_add_epi16(mx, meW);
			my=_mm_add_epi16(my, meN);
			mr=_mm_add_epi16(mr, _mm_set1_epi16(1));
			mx=_mm_add_epi16(mx, _mm_set1_epi16(1));
			my=_mm_add_epi16(my, _mm_set1_epi16(1));
			mdiff=_mm_cvtepi16_epi32(mdiff);
			mr=_mm_cvtepi16_epi32(mr);
			mx=_mm_cvtepi16_epi32(mx);
			my=_mm_cvtepi16_epi32(my);
			mr=_mm_castps_si128(_mm_cvtepi32_ps(mr));
			mx=_mm_castps_si128(_mm_cvtepi32_ps(mx));
			my=_mm_castps_si128(_mm_cvtepi32_ps(my));
			mr=_mm_sub_epi32(mr, _mm_set1_epi32(127<<23));
			mx=_mm_sub_epi32(mx, _mm_set1_epi32(127<<23));
			my=_mm_sub_epi32(my, _mm_set1_epi32(127<<23));

			mr=_mm_and_si128(mr, _mm_set1_epi32(0x7F800000));
			mx=_mm_slli_epi32(mx, 1);
			my=_mm_slli_epi32(my, 1);

			mr=_mm_srli_epi32(mr, 23);
			mx=_mm_srli_epi32(mx, 24);
			my=_mm_srli_epi32(my, 24);
			mr=_mm_or_si128(mr, _mm_and_si128(_mm_srai_epi32(mdiff, 31), _mm_set1_epi32(1)));
			_mm_store_si128((__m128i*)grads+0, mx);
			_mm_store_si128((__m128i*)grads+1, my);
			_mm_store_si128((__m128i*)grads+2, mr);
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
			for(int kc=0;kc<3;++kc)
			{
				int offset=yuv[combination[kc+3]];
				int qeW=grads[kc+0];
				int qeN=grads[kc+4];
				int qeD=grads[kc+8];
				switch(predidx[kc])
				{
				case PRED_CG:
					MEDIAN3_32(pred, N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
					break;
#ifdef ENABLE_AV4
				case PRED_AV4:
					CLAMP3_32(pred,
						(4*(N[kc]+W[kc])+NE[kc]-NW[kc])>>3,
						N[kc], W[kc], NE[kc]
					);
					break;
#endif
#ifdef ENABLE_AV9
				case PRED_AV9:
					CLAMP3_32(pred,
						W[kc]+((10*N[kc]-9*NW[kc]+4*NE[kc]-2*(NN[kc]+WW[kc])+NNW[kc]-(NNE[kc]+NWW[kc]))>>4),
						N[kc], W[kc], NE[kc]
					);
					break;
#endif
#ifdef ENABLE_WG
				case PRED_WG:
					{
						int gx=abs(W[kc]-WW[kc])+abs(N[kc]-NW[kc])+abs(NE[kc]-N[kc])+1;
						int gy=abs(W[kc]-NW[kc])+abs(N[kc]-NN[kc])+abs(NE[kc]-NNE[kc])+1;
						pred=(gx*N[kc]+gy*W[kc])/(gx+gy);
					}
					break;
#endif
				}
				pred+=offset;
				CLAMP2(pred, -128, 127);
				//CLAMP2_32(pred, pred, -128, 127);

				int cdf, freq=0, den;
				unsigned short *curr_hist=histptr+NCTX*cdfstride*kc;
				unsigned short *curr_hist0=curr_hist;
				unsigned short *curr_hist1=curr_hist+cdfstride*(CLEVELS0+((128+(N[kc]+W[kc])/2)&255));
				unsigned short *curr_hist2=curr_hist+cdfstride*(CLEVELS0+CLEVELS1+qeD);
				unsigned short *curr_hist3=curr_hist+cdfstride*(CLEVELS0+CLEVELS1+CLEVELS2+CLEVELS3*qeN+qeW);
				unsigned short *curr_hist4=curr_hist+cdfstride*(CLEVELS0+CLEVELS1+CLEVELS2+CLEVELS3*CLEVELS3+FLOOR_LOG2_P1(NN[kc+4]*NN[kc+4]+N[kc+4]*N[kc+4]+NW[kc+4]*NW[kc+4]+NE[kc+4]*NE[kc+4]));
				unsigned short *curr_hist5=curr_hist+cdfstride*(CLEVELS0+CLEVELS1+CLEVELS2+CLEVELS3*CLEVELS3+CLEVELS4+FLOOR_LOG2_P1(W[kc+4]*W[kc+4]+NW[kc+4]*NW[kc+4]+WW[kc+4]*WW[kc+4]+N[kc+4]*N[kc+4]));
				unsigned short *curr_hist6=curr_hist+cdfstride*(CLEVELS0+CLEVELS1+CLEVELS2+CLEVELS3*CLEVELS3+CLEVELS4+CLEVELS5+((128+pred)&255));
				unsigned short *curr_hist7=curr_hist+cdfstride*(CLEVELS0+CLEVELS1+CLEVELS2+CLEVELS3*CLEVELS3+CLEVELS4+CLEVELS5+CLEVELS6+(((N[kc]-NW[kc])>>5&7)<<6|((W[kc]-NW[kc])>>5&7)<<3|(pred>>5&7)));
				int f0, f1, f2, f3, f4, f5, f6, f7;
				ptrdiff_t
					w0=mixer[kc][0],
					w1=mixer[kc][1],
					w2=mixer[kc][2],
					w3=mixer[kc][3],
					w4=mixer[kc][4],
					w5=mixer[kc][5],
					w6=mixer[kc][6],
					w7=mixer[kc][7];
				ptrdiff_t wsum=w0+w1+w2+w3+w4+w5+w6+w7+1;
				w1=(w1<<16)/wsum;
				w2=(w2<<16)/wsum;
				w3=(w3<<16)/wsum;
				w4=(w4<<16)/wsum;
				w5=(w5<<16)/wsum;
				w6=(w6<<16)/wsum;
				w7=(w7<<16)/wsum;
				w0=0x10000-(w1+w2+w3+w4+w5+w6+w7);
	#define GETCTR(X) (int)((w0*curr_hist0[X]+w1*curr_hist1[X]+w2*curr_hist2[X]+w3*curr_hist3[X]+w4*curr_hist4[X]+w5*curr_hist5[X]+w6*curr_hist6[X]+w7*curr_hist7[X]+0xFFFF)>>16)
//	#define GETCTR(X) (int)((curr_hist0[X]+curr_hist1[X]+curr_hist2[X]+curr_hist3[X]+curr_hist4[X]+curr_hist5[X]+curr_hist6[X]+curr_hist7[X]+1)>>1)
//	#define GETCTR(X) curr_hist0[X]		//fast
				den=GETCTR(args->tlevels)+(args->tlevels<<1);//there are TLEVELS roundings max

				//if(ky==480&&kx==109&&kc==2)//
				//if(ky==147&&kx==317&&kc==1)//
				//if(ky==1&&kx==455&&kc==0)//
				//if(ky==0&&kx==5&&kc==1)//
				//if(ky==0&&kx==0&&kc==0)//
				//	printf("");
				if(args->fwd)
				{
					error=yuv[kc]-pred;
					{
						int upred=half-abs(pred), aval=abs(error);
						if(aval<=upred)
							sym=error<<1^error>>31;//pack sign
						else
							sym=upred+aval;//error sign is known
					}
					quantize_pixel(sym, &token, &bypass, &nbits);
#ifdef _DEBUG
					if(token>=args->tlevels)
						LOG_ERROR("YXC %d %d %d  token %d/%d", ky, kx, kc, token, args->tlevels);
#endif
					cdf=token;
					freq=GETCTR(token)+1;
					switch(token)
					{
					case 25:cdf+=GETCTR(24);
					case 24:cdf+=GETCTR(23);
					case 23:cdf+=GETCTR(22);
					case 22:cdf+=GETCTR(21);
					case 21:cdf+=GETCTR(20);
					case 20:cdf+=GETCTR(19);
					case 19:cdf+=GETCTR(18);
					case 18:cdf+=GETCTR(17);
					case 17:cdf+=GETCTR(16);
					case 16:cdf+=GETCTR(15);
					case 15:cdf+=GETCTR(14);
					case 14:cdf+=GETCTR(13);
					case 13:cdf+=GETCTR(12);
					case 12:cdf+=GETCTR(11);
					case 11:cdf+=GETCTR(10);
					case 10:cdf+=GETCTR( 9);
					case  9:cdf+=GETCTR( 8);
					case  8:cdf+=GETCTR( 7);
					case  7:cdf+=GETCTR( 6);
					case  6:cdf+=GETCTR( 5);
					case  5:cdf+=GETCTR( 4);
					case  4:cdf+=GETCTR( 3);
					case  3:cdf+=GETCTR( 2);
					case  2:cdf+=GETCTR( 1);
					case  1:cdf+=GETCTR( 0);
					default:
						break;
					}
					//for(int ks=0;ks<token;++ks)
					//	cdf+=GETCTR(ks);
					//freq=GETCTR(token);

					ac3_enc_update_NPOT(&ec, cdf, freq, den);
					if(nbits)
						ac3_enc_bypass(&ec, bypass, nbits);
				}
				else
				{
					unsigned code=ac3_dec_getcdf_NPOT(&ec, den);
					cdf=0;
					freq=0;
					token=0;
#if 1
					for(;;)
					{
						unsigned cdf2;
						freq=GETCTR(token)+1;
						cdf2=cdf+freq;
						if(cdf2>code)
							break;
#ifdef _DEBUG
						if(token>=args->tlevels)
							LOG_ERROR("YXC %d %d %d  token %d/%d", ky, kx, kc, token, args->tlevels);
#endif
						cdf=cdf2;
						++token;
					}
#endif
					ac3_dec_update_NPOT(&ec, cdf, freq, den);
					sym=token;
					if(sym>=(1<<CONFIG_EXP))
					{
						int lsb, msb;

						sym-=1<<CONFIG_EXP;
						lsb=sym&((1<<CONFIG_LSB)-1);
						sym>>=CONFIG_LSB;
						msb=sym&((1<<CONFIG_MSB)-1);
						sym>>=CONFIG_MSB;
						nbits=sym+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB);
						bypass=ac3_dec_bypass(&ec, nbits);
						sym=1;
						sym<<=CONFIG_MSB;
						sym|=msb;
						sym<<=nbits;
						sym|=bypass;
						sym<<=CONFIG_LSB;
						sym|=lsb;
					}
					{
						int upred=half-abs(pred);
						if(sym<=(upred<<1))
							error=sym>>1^-(sym&1);
						else
						{
							int negmask=-(pred>0);
							error=sym-upred;
							error^=negmask;
							error-=negmask;
						}
					}
					yuv[kc]=error+pred;
				}
				curr[kc+0]=yuv[kc]-offset;
				curr[kc+4]=(2*W[kc+4]+abs(error)+NEEE[kc+4])>>2;

				f0=curr_hist0[token];
				f1=curr_hist1[token];
				f2=curr_hist2[token];
				f3=curr_hist3[token];
				f4=curr_hist4[token];
				f5=curr_hist5[token];
				f6=curr_hist6[token];
				f7=curr_hist7[token];
				mixer[kc][0]+=f0;
				mixer[kc][1]+=f1;
				mixer[kc][2]+=f2;
				mixer[kc][3]+=f3;
				mixer[kc][4]+=f4;
				mixer[kc][5]+=f5;
				mixer[kc][6]+=f6;
				mixer[kc][7]+=f7;
				w0>>=9;
				w1>>=9;
				w2>>=9;
				w3>>=9;
				w4>>=9;
				w5>>=9;
				w6>>=9;
				w7>>=9;
				curr_hist0[token]=(unsigned short)(f0+w0);
				curr_hist1[token]=(unsigned short)(f1+w1);
				curr_hist2[token]=(unsigned short)(f2+w2);
				curr_hist3[token]=(unsigned short)(f3+w3);
				curr_hist4[token]=(unsigned short)(f4+w4);
				curr_hist5[token]=(unsigned short)(f5+w5);
				curr_hist6[token]=(unsigned short)(f6+w6);
				curr_hist7[token]=(unsigned short)(f7+w7);
				curr_hist0[args->tlevels]+=(unsigned short)w0;
				curr_hist1[args->tlevels]+=(unsigned short)w1;
				curr_hist2[args->tlevels]+=(unsigned short)w2;
				curr_hist3[args->tlevels]+=(unsigned short)w3;
				curr_hist4[args->tlevels]+=(unsigned short)w4;
				curr_hist5[args->tlevels]+=(unsigned short)w5;
				curr_hist6[args->tlevels]+=(unsigned short)w6;
				curr_hist7[args->tlevels]+=(unsigned short)w7;
				if(wsum+f0+f1+f2+f3+f4+f5+f6+f7>0x30000)
				{
					mixer[kc][0]>>=1;
					mixer[kc][1]>>=1;
					mixer[kc][2]>>=1;
					mixer[kc][3]>>=1;
					mixer[kc][4]>>=1;
					mixer[kc][5]>>=1;
					mixer[kc][6]>>=1;
					mixer[kc][7]>>=1;
				}

				if(curr_hist0[args->tlevels]>=histlimit)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist0[ks]-=curr_hist0[ks]>>1;
					curr_hist0[args->tlevels]=sum;
				}
				if(curr_hist1[args->tlevels]>=histlimit)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist1[ks]-=curr_hist1[ks]>>1;
					curr_hist1[args->tlevels]=sum;
				}
				if(curr_hist2[args->tlevels]>=histlimit)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist2[ks]-=curr_hist2[ks]>>1;
					curr_hist2[args->tlevels]=sum;
				}
				if(curr_hist3[args->tlevels]>=histlimit)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist3[ks]-=curr_hist3[ks]>>1;
					curr_hist3[args->tlevels]=sum;
				}
				if(curr_hist4[args->tlevels]>=histlimit)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist4[ks]-=curr_hist4[ks]>>1;
					curr_hist4[args->tlevels]=sum;
				}
				if(curr_hist5[args->tlevels]>=histlimit)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist5[ks]-=curr_hist5[ks]>>1;
					curr_hist5[args->tlevels]=sum;
				}
				if(curr_hist6[args->tlevels]>=histlimit)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist6[ks]-=curr_hist6[ks]>>1;
					curr_hist6[args->tlevels]=sum;
				}
				if(curr_hist7[args->tlevels]>=histlimit)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist7[ks]-=curr_hist7[ks]>>1;
					curr_hist7[args->tlevels]=sum;
				}
			}
			if(!args->fwd)
			{
				args->dst[idx+yidx]=yuv[0]+128;
				args->dst[idx+uidx]=yuv[1]+128;
				args->dst[idx+vidx]=yuv[2]+128;
#ifdef ENABLE_GUIDE
				if(args->test&&memcmp(args->dst+idx, args->src+idx, sizeof(char[3])))
				{
					unsigned char orig[4]={0};
					memcpy(orig, args->src+idx, sizeof(char[3]));
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	if(args->fwd)
		ac3_enc_flush(&ec);
}
int c25_codec(const char *srcfn, const char *dstfn)
{
	const int depth=8;
	double t0;
	ArrayHandle src, dst;
	int headersize, printed;
	int iw, ih;
	const unsigned char *image, *imageend;
	unsigned char *image2;
	CodecID codec;
	int ncores=query_cpu_cores();
	int xblocks, yblocks, nblocks, nthreads, coffset;
	ptrdiff_t start, memusage, argssize;
	ThreadArgs *args;
	int test, fwd;
	int tlevels, histsize, statssize;
	double esize;
	int usize;
	int maxwidth;
	
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
	xblocks=(iw+BLOCKSIZE-1)/BLOCKSIZE;
	yblocks=(ih+BLOCKSIZE-1)/BLOCKSIZE;
	nblocks=xblocks*yblocks, nthreads=MINVAR(nblocks, ncores);
	coffset=(int)sizeof(int)*nblocks;
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
	if(fwd)
	{
		dst=0;
		printed=snprintf(g_buf, G_BUF_SIZE-1, "C01\n%d %d\n", iw, ih);
		array_append(&dst, g_buf, 1, printed, 1, 0, 0);
		start=array_append(&dst, 0, 1, coffset, 1, 0, 0);
		
		image2=0;
	}
	else//integrity check
	{
		dst=0;
		printed=snprintf(g_buf, G_BUF_SIZE-1, "P6\n%d %d\n255\n", iw, ih);
		array_append(&dst, g_buf, 1, printed, 1, 0, 0);
		array_append(&dst, 0, 1, usize, 1, 0, 0);

		//printed=0;
		start=coffset;
		for(int kt=0;kt<nblocks;++kt)
		{
			int size=0;
			memcpy(&size, image+sizeof(int)*kt, sizeof(int));
			start+=size;
		}
		if(image+start!=imageend)
			LOG_ERROR("Corrupt file");
		start=coffset;

		image2=(unsigned char*)malloc(usize);
		if(!image2)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		memset(image2, 0, usize);
	}
	{
		int nlevels=256;
		int token=0, bypass=0, nbits=0;

		quantize_pixel(nlevels, &token, &bypass, &nbits);
		tlevels=token+1;
		statssize=(tlevels+1)*(int)sizeof(short[3][NCTX]);
		histsize=(int)sizeof(int[OCH_COUNT*PRED_COUNT<<8]);
		if(histsize<statssize)
			histsize=statssize;
		maxwidth=iw;
		if(maxwidth>BLOCKSIZE)
			maxwidth=BLOCKSIZE;
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=image;
		arg->dst=fwd?0:dst->data+printed;
		arg->iw=iw;
		arg->ih=ih;
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
		
		arg->tlevels=tlevels;
		arg->fwd=fwd;
		arg->test=test;
#ifdef ENABLE_MT
		arg->loud=0;
#else
		arg->loud=test&&nblocks<MAXPRINTEDBLOCKS;
#endif
	}
	for(int k2=0;k2<=test;++k2)
	{
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
				arg->x1=BLOCKSIZE*kx;
				arg->y1=BLOCKSIZE*ky;
				arg->x2=MINVAR(arg->x1+BLOCKSIZE, iw);
				arg->y2=MINVAR(arg->y1+BLOCKSIZE, ih);
				if(!fwd)
				{
					int size=0;
					memcpy(&size, image+sizeof(int)*((ptrdiff_t)kt+kt2), sizeof(int));
					arg->decstart=image+start;
					start+=size;
					arg->decend=image+start;
				}
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
						int blocksize=((arg->x2-arg->x1)*(arg->y2-arg->y1)*3*depth+7)>>3;
						int kx, ky;

						kx=kt+kt2;
						ky=kx/xblocks;
						kx%=xblocks;
						if(nblocks<MAXPRINTEDBLOCKS)
						{
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
								arg->list.nbytes,
								arg->list.nbytes-arg->bestsize,
								100.*arg->list.nbytes/blocksize,
								(double)blocksize/arg->list.nbytes,
								rct_names[arg->bestrct],
								pred_names[arg->predidx[0]],
								pred_names[arg->predidx[1]],
								pred_names[arg->predidx[2]]
							);
						}
						esize+=arg->bestsize;
#ifdef ABAC_PROFILE_SIZE
						csizes[0]+=arg->csizes[0];
						csizes[1]+=arg->csizes[1];
						csizes[2]+=arg->csizes[2];
#endif
					}
					memcpy(dst->data+start+sizeof(int)*((ptrdiff_t)kt+kt2), &arg->list.nbytes, sizeof(int));
					blist_appendtoarray(&arg->list, &dst);
					blist_clear(&arg->list);
				}
			}
		}
		if(test)
		{
			ptrdiff_t usize=((ptrdiff_t)iw*ih*3*depth+7)>>3;
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
	}
	free(args);
	array_free(&src);
	array_free(&dst);
	return 0;
}