#include"codec.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs, log2
//#include<immintrin.h>//included by "entropy.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
	#define ENABLE_MT

//	#define ENABLE_ENTROPY_ANALYSIS		//best
	#define ENABLE_SAD_ANALYSIS		//fastest
//	#define ENABLE_EXACT_ANALYSIS		//bad & slow

	#define ENABLE_MIXRCT
	#define ENABLE_NONCAUSALRCT

	#define ENABLE_MIX4			//slower
	#define ENABLE_AV4			//faster

	#define ENABLE_EDGECASES		//slower

	#define ENABLE_QUANTIZATION		//good


#define ANALYSIS_XSTRIDE 4
#define ANALYSIS_YSTRIDE 4

#define AC3_PREC
#include"entropy.h"

#define BLOCKSIZE 512
#define MAXPRINTEDBLOCKS 500
#define CLEVELS 13
#define DLEVELS 11

#define MIXBITS 14

#define OCHLIST\
	OCH(R)\
	OCH(G)\
	OCH(B)\
	OCH(RG)\
	OCH(GB)\
	OCH(BR)\
	OCH(R31)\
	OCH(G31)\
	OCH(B31)\
	OCH(R22)\
	OCH(G22)\
	OCH(B22)\
	OCH(R13)\
	OCH(G13)\
	OCH(B13)\
	OCH(RB)\
	OCH(GR)\
	OCH(BG)\
	OCH(RGBY)\
	OCH(RGBU)\
	OCH(RGBV)\
	OCH(RGBU2)\
	OCH(RGBV2)\
	OCH(GBRY)\
	OCH(GBRU)\
	OCH(GBRV)\
	OCH(GBRU2)\
	OCH(GBRV2)\
	OCH(BRGY)\
	OCH(BRGU)\
	OCH(BRGV)\
	OCH(BRGU2)\
	OCH(BRGV2)
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
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		0,	0, 0)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		0,	0, 4)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		0,	4, 0)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		0,	4, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		0,	0, 4)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		0,	0, 4)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		0,	4, 0)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		4,	4, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		4,	0, 4)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		4,	0, 4)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		4,	4, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		4,	0, 4)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		4,	0, 4)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		4,	4, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		4,	0, 4)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		4,	0, 4)\
	RCT(R_G_B31,	OCH_R,		OCH_G,		OCH_B31,	0,	1, 3)\
	RCT(G_B_R31,	OCH_G,		OCH_B,		OCH_R31,	0,	3, 1)\
	RCT(B_R_G31,	OCH_B,		OCH_R,		OCH_G31,	0,	1, 3)\
	RCT(G_BG_R31,	OCH_G,		OCH_BG,		OCH_R31,	4,	3, 1)\
	RCT(G_RG_B31,	OCH_G,		OCH_RG,		OCH_B31,	4,	1, 3)\
	RCT(B_RB_G31,	OCH_B,		OCH_RB,		OCH_G31,	4,	3, 1)\
	RCT(B_GB_R31,	OCH_B,		OCH_GB,		OCH_R31,	4,	1, 3)\
	RCT(R_GR_B31,	OCH_R,		OCH_GR,		OCH_B31,	4,	3, 1)\
	RCT(R_BR_G31,	OCH_R,		OCH_BR,		OCH_G31,	4,	1, 3)\
	RCT(R_G_B22,	OCH_R,		OCH_G,		OCH_B22,	0,	2, 2)\
	RCT(G_B_R22,	OCH_G,		OCH_B,		OCH_R22,	0,	2, 2)\
	RCT(B_R_G22,	OCH_B,		OCH_R,		OCH_G22,	0,	2, 2)\
	RCT(G_BG_R22,	OCH_G,		OCH_BG,		OCH_R22,	4,	2, 2)\
	RCT(G_RG_B22,	OCH_G,		OCH_RG,		OCH_B22,	4,	2, 2)\
	RCT(B_RB_G22,	OCH_B,		OCH_RB,		OCH_G22,	4,	2, 2)\
	RCT(B_GB_R22,	OCH_B,		OCH_GB,		OCH_R22,	4,	2, 2)\
	RCT(R_GR_B22,	OCH_R,		OCH_GR,		OCH_B22,	4,	2, 2)\
	RCT(R_BR_G22,	OCH_R,		OCH_BR,		OCH_G22,	4,	2, 2)\
	RCT(R_G_B13,	OCH_R,		OCH_G,		OCH_B13,	0,	3, 1)\
	RCT(G_B_R13,	OCH_G,		OCH_B,		OCH_R13,	0,	1, 3)\
	RCT(B_R_G13,	OCH_B,		OCH_R,		OCH_G13,	0,	3, 1)\
	RCT(G_BG_R13,	OCH_G,		OCH_BG,		OCH_R13,	4,	1, 3)\
	RCT(G_RG_B13,	OCH_G,		OCH_RG,		OCH_B13,	4,	3, 1)\
	RCT(B_RB_G13,	OCH_B,		OCH_RB,		OCH_G13,	4,	1, 3)\
	RCT(B_GB_R13,	OCH_B,		OCH_GB,		OCH_R13,	4,	3, 1)\
	RCT(R_GR_B13,	OCH_R,		OCH_GR,		OCH_B13,	4,	1, 3)\
	RCT(R_BR_G13,	OCH_R,		OCH_BR,		OCH_G13,	4,	3, 1)\
	RCT(RGB1,	OCH_RGBY,	OCH_RGBU,	OCH_RGBV,	0,	0, 0)\
	RCT(RGB2,	OCH_RGBY,	OCH_RGBU2,	OCH_RGBV,	0,	0, 0)\
	RCT(RGB3,	OCH_RGBY,	OCH_RGBU,	OCH_RGBV2,	0,	0, 0)\
	RCT(GBR1,	OCH_GBRY,	OCH_GBRU,	OCH_GBRV,	0,	0, 0)\
	RCT(GBR2,	OCH_GBRY,	OCH_GBRU2,	OCH_GBRV,	0,	0, 0)\
	RCT(GBR3,	OCH_GBRY,	OCH_GBRU,	OCH_GBRV2,	0,	0, 0)\
	RCT(BRG1,	OCH_BRGY,	OCH_BRGU,	OCH_BRGV,	0,	0, 0)\
	RCT(BRG2,	OCH_BRGY,	OCH_BRGU2,	OCH_BRGV,	0,	0, 0)\
	RCT(BRG3,	OCH_BRGY,	OCH_BRGU,	OCH_BRGV2,	0,	0, 0)
typedef enum _RCTIndex
{
#define RCT(LABEL, YIDX, UIDX, VIDX, UHELP, VH0, VH1) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
typedef enum _RCTInfoIndex
{
	RII_YOCH,
	RII_UOCH,
	RII_VOCH,

	RII_UHELP,
	RII_VH0,
	RII_VH1,

	RII_COUNT,
} RCTInfoIndex;
static const unsigned char rct_combinations[RCT_COUNT][RII_COUNT]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX, UHELP, VH0, VH1) {YIDX, UIDX, VIDX, UHELP, VH0, VH1},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX, UHELP, VH0, VH1) #LABEL,
	RCTLIST
#undef  RCT
};

#define PREDLIST\
	PRED(AV2)\
	PRED(CG)\
	PRED(AV5)\
	PRED(AV9)
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
typedef enum _NBIndex
{
			NB_NNW,		NB_NN,		NB_NNE,
	NB_NWW,		NB_NW,		NB_N,		NB_NE,
	NB_WW,		NB_W,		NB_curr,

	NB_COUNT,
} NBIndex;

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
#define CONFIG_EXP 4
#define CONFIG_MSB 1
#define CONFIG_LSB 0
FORCEINLINE void quantize_pixel(int val, int *token, int *bypass, int *nbits)
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
#if defined ENABLE_MIX4 && !defined ENABLE_AV4
FORCEINLINE int f28_mix4(int v00, int v01, int v10, int v11, int alphax, int alphay)
{
	//v00=v00*((1<<12)-alphax)+v01*alphax;
	v00=((v00<<MIXBITS)+(v01-v00)*alphax)>>(MIXBITS-1);
	v10=((v10<<MIXBITS)+(v11-v10)*alphax)>>(MIXBITS-1);
	v00=((v00<<MIXBITS)+(v10-v00)*alphay)>>(MIXBITS-1);
	return v00;
}
#endif
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
	ThreadArgs *args=(ThreadArgs*)param;
	AC3 ec;
	const unsigned char *image=args->fwd?args->src:args->dst;
	unsigned char bestrct=0, predidx[4]={0};
	const unsigned char *combination=0;
	int ystride=args->iw*3;
	int cdfstride=args->tlevels+1;
	int token=0, bypass=0, nbits=0;
	int halfs[]={128, 128, 128};

	if(args->fwd)
	{
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		unsigned char predsel[OCH_COUNT]={0};
		int dx=(args->x2-args->x1-3)/ANALYSIS_XSTRIDE/5*5, dy=(args->y2-args->y1-2)/ANALYSIS_YSTRIDE;
		if(dx<=0||dy<=0)
		{
			bestrct=RCT_G_BG_RG;
			predidx[0]=PRED_CG;
			predidx[1]=PRED_CG;
			predidx[2]=PRED_CG;
			combination=rct_combinations[bestrct];
			goto skip_analysis;
		}
		
		int count=0;
#ifdef ENABLE_SAD_ANALYSIS
		long long counters[PRED_COUNT*OCH_COUNT]={0};
#elif defined ENABLE_EXACT_ANALYSIS
		int bypasssizes[PRED_COUNT*OCH_COUNT]={0};
		memset(args->hist, 0, args->histsize);
#elif defined ENABLE_ENTROPY_ANALYSIS
		memset(args->hist, 0, args->histsize);
#endif
		for(int ky=args->y1+2;ky<args->y2;ky+=ANALYSIS_YSTRIDE)//analysis loop
		{
			int kx=args->x1+2;
			const unsigned char *ptr=image+3*(args->iw*ky+kx);

			__m256i amin=_mm256_set1_epi16(-128);
			__m256i amax=_mm256_set1_epi16(127);
#ifdef ENABLE_ENTROPY_ANALYSIS
			__m256i amask=_mm256_set1_epi16(255);
#endif
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
			ALIGN(32) short result[16]={0};
			for(;kx<=args->x2-(5*ANALYSIS_XSTRIDE+1+1);kx+=5*ANALYSIS_XSTRIDE, ptr+=15*ANALYSIS_XSTRIDE, count+=5)
			{
				__m256i
					nb0[NB_COUNT],//rgb
					nb1[NB_COUNT],//gbr			FIXME need just curr
					nb2[NB_COUNT],//rgb - gbr
					nb3[NB_COUNT],//gbr - rgb
#ifdef ENABLE_MIXRCT
					nb4[NB_COUNT],//(3*gbr+brg)/4		FIXME need just curr
					nb5[NB_COUNT],//rgb - (3*gbr+brg)/4
					nb6[NB_COUNT],//(gbr+brg)/2		FIXME need just curr
					nb7[NB_COUNT],//rgb - (gbr+brg)/2
					nb8[NB_COUNT],//(gbr+3*brg)/4		FIXME need just curr
					nb9[NB_COUNT],//rgb - (gbr+brg)/4
#endif
					nby[NB_COUNT],//mod256(rgb + 2*gbr + brg)
					nbu[NB_COUNT],//mod256(rgb - gbr)
					nbv[NB_COUNT],//mod256(brg - gbr)
					nbu2[NB_COUNT],//mod256(u - v/4)
					nbv2[NB_COUNT];//mod256(v - u/4)
				__m256i vmin[11], vmax[11], pred;
				{
					__m128i nbb[NB_COUNT]=//8-bit
					{
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-2*ystride-1*3+0)), half8),//NNW
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-2*ystride+0*3+0)), half8),//NN
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-2*ystride+1*3+0)), half8),//NNE
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride-2*3+0)), half8),//NWW
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride-1*3+0)), half8),//NW
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride+0*3+0)), half8),//N
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride+1*3+0)), half8),//NE
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr+0*ystride-2*3+0)), half8),//WW
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr+0*ystride-1*3+0)), half8),//W
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr+0*ystride+0*3+0)), half8),//curr
					};
#ifdef __GNUC__
#pragma GCC unroll 6
#endif
					for(int k=0;k<NB_COUNT;++k)
					{
						__m128i temp;
						__m256i t1;
						nb0[k]=_mm256_cvtepi8_epi16(nbb[k]);
						temp=_mm_shuffle_epi8(nbb[k], shuf);
						nb1[k]=_mm256_cvtepi8_epi16(temp);
						nb2[k]=_mm256_sub_epi16(nb0[k], nb1[k]);
						nb3[k]=_mm256_sub_epi16(nb1[k], nb0[k]);
						
#ifdef ENABLE_MIXRCT
						t1=_mm256_cvtepi8_epi16(_mm_shuffle_epi8(temp, shuf));
						nb4[k]=_mm256_srai_epi16(_mm256_add_epi16(t1, nb1[k]), 1);
						nb6[k]=_mm256_srai_epi16(_mm256_add_epi16(_mm256_add_epi16(nb1[k], t1), _mm256_add_epi16(nb1[k], nb1[k])), 2);
						nb8[k]=_mm256_srai_epi16(_mm256_add_epi16(_mm256_add_epi16(nb1[k], t1), _mm256_add_epi16(t1, t1)), 2);
						nb5[k]=_mm256_sub_epi16(nb0[k], nb4[k]);
						nb7[k]=_mm256_sub_epi16(nb0[k], nb6[k]);
						nb9[k]=_mm256_sub_epi16(nb0[k], nb8[k]);
#endif

						nbu[k]=nb2[k];
					//	nbu[k]=_mm256_srai_epi16(_mm256_slli_epi16(nbu[k], 8), 8);
						temp=_mm_shuffle_epi8(temp, shuf);
						nbv[k]=_mm256_sub_epi16(_mm256_cvtepi8_epi16(temp), nb1[k]);
					//	nbv[k]=_mm256_srai_epi16(_mm256_slli_epi16(nbv[k], 8), 8);
						nby[k]=_mm256_add_epi16(nb1[k], _mm256_srai_epi16(_mm256_add_epi16(nbu[k], nbv[k]), 2));
					//	nby[k]=_mm256_srai_epi16(_mm256_slli_epi16(nby[k], 8), 8);
						nbu2[k]=_mm256_sub_epi16(nbu[k], _mm256_srai_epi16(nbv[k], 2));
						nbv2[k]=_mm256_sub_epi16(nbv[k], _mm256_srai_epi16(nbu[k], 2));
					//	nbu2[k]=_mm256_srai_epi16(_mm256_slli_epi16(nbu2[k], 8), 8);
					//	nbv2[k]=_mm256_srai_epi16(_mm256_slli_epi16(nbv2[k], 8), 8);
					}
				}
#ifdef ENABLE_SAD_ANALYSIS
#define UPDATE(PREDIDX, IDX0, IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7, IDX8, IDX9, IDXA, IDXB, IDXC, IDXD, IDXE)\
	do\
	{\
		pred=_mm256_abs_epi16(pred);\
		_mm256_store_si256((__m256i*)result, pred);\
		counters[IDX0*PRED_COUNT+PREDIDX]+=result[0x0];\
		counters[IDX1*PRED_COUNT+PREDIDX]+=result[0x1];\
		counters[IDX2*PRED_COUNT+PREDIDX]+=result[0x2];\
		counters[IDX3*PRED_COUNT+PREDIDX]+=result[0x3];\
		counters[IDX4*PRED_COUNT+PREDIDX]+=result[0x4];\
		counters[IDX5*PRED_COUNT+PREDIDX]+=result[0x5];\
		counters[IDX6*PRED_COUNT+PREDIDX]+=result[0x6];\
		counters[IDX7*PRED_COUNT+PREDIDX]+=result[0x7];\
		counters[IDX8*PRED_COUNT+PREDIDX]+=result[0x8];\
		counters[IDX9*PRED_COUNT+PREDIDX]+=result[0x9];\
		counters[IDXA*PRED_COUNT+PREDIDX]+=result[0xA];\
		counters[IDXB*PRED_COUNT+PREDIDX]+=result[0xB];\
		counters[IDXC*PRED_COUNT+PREDIDX]+=result[0xC];\
		counters[IDXD*PRED_COUNT+PREDIDX]+=result[0xD];\
		counters[IDXE*PRED_COUNT+PREDIDX]+=result[0xE];\
	}while(0)
#elif defined ENABLE_EXACT_ANALYSIS
#define UPDATE(PREDIDX, IDX0, IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7, IDX8, IDX9, IDXA, IDXB, IDXC, IDXD, IDXE)\
	do\
	{\
		pred=_mm256_xor_si256(_mm256_slli_epi16(pred, 1), _mm256_srai_epi16(pred, 15));\
		_mm256_store_si256((__m256i*)result, pred);\
		quantize_pixel(result[0x0], &token, &bypass, &nbits); bypasssizes[IDX0*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDX0*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0x1], &token, &bypass, &nbits); bypasssizes[IDX1*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDX1*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0x2], &token, &bypass, &nbits); bypasssizes[IDX2*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDX2*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0x3], &token, &bypass, &nbits); bypasssizes[IDX3*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDX3*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0x4], &token, &bypass, &nbits); bypasssizes[IDX4*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDX4*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0x5], &token, &bypass, &nbits); bypasssizes[IDX5*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDX5*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0x6], &token, &bypass, &nbits); bypasssizes[IDX6*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDX6*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0x7], &token, &bypass, &nbits); bypasssizes[IDX7*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDX7*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0x8], &token, &bypass, &nbits); bypasssizes[IDX8*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDX8*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0x9], &token, &bypass, &nbits); bypasssizes[IDX9*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDX9*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0xA], &token, &bypass, &nbits); bypasssizes[IDXA*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDXA*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0xB], &token, &bypass, &nbits); bypasssizes[IDXB*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDXB*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0xC], &token, &bypass, &nbits); bypasssizes[IDXC*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDXC*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0xD], &token, &bypass, &nbits); bypasssizes[IDXD*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDXD*PRED_COUNT+PREDIDX)<<6|token];\
		quantize_pixel(result[0xE], &token, &bypass, &nbits); bypasssizes[IDXE*PRED_COUNT+PREDIDX]+=nbits; ++args->hist[(IDXE*PRED_COUNT+PREDIDX)<<6|token];\
	}while(0)
#elif defined ENABLE_ENTROPY_ANALYSIS
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
#endif
				//AV2 = (N+W)>>1
				pred=_mm256_srai_epi16(_mm256_add_epi16(nb0[NB_N], nb0[NB_W]), 1);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_srai_epi16(_mm256_add_epi16(nb2[NB_N], nb2[NB_W]), 1);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_srai_epi16(_mm256_add_epi16(nb3[NB_N], nb3[NB_W]), 1);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#ifdef ENABLE_MIXRCT
				pred=_mm256_srai_epi16(_mm256_add_epi16(nb5[NB_N], nb5[NB_W]), 1);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31
				);
				pred=_mm256_srai_epi16(_mm256_add_epi16(nb7[NB_N], nb7[NB_W]), 1);

				pred=_mm256_add_epi16(pred, nb6[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22
				);
				pred=_mm256_srai_epi16(_mm256_add_epi16(nb9[NB_N], nb9[NB_W]), 1);

				pred=_mm256_add_epi16(pred, nb8[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13
				);
#endif
#ifdef ENABLE_NONCAUSALRCT
				pred=_mm256_srai_epi16(_mm256_add_epi16(nby[NB_N], nby[NB_W]), 1);

				pred=_mm256_sub_epi16(nby[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY
				);
				pred=_mm256_srai_epi16(_mm256_add_epi16(nbu[NB_N], nbu[NB_W]), 1);

				pred=_mm256_sub_epi16(nbu[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU
				);
				pred=_mm256_srai_epi16(_mm256_add_epi16(nbv[NB_N], nbv[NB_W]), 1);

				pred=_mm256_sub_epi16(nbv[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV
				);
				pred=_mm256_srai_epi16(_mm256_add_epi16(nbu2[NB_N], nbu2[NB_W]), 1);

				pred=_mm256_sub_epi16(nbu2[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2
				);
				pred=_mm256_srai_epi16(_mm256_add_epi16(nbv2[NB_N], nbv2[NB_W]), 1);

				pred=_mm256_sub_epi16(nbv2[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2
				);
#endif

				//CG = median(N, W, N+W-NW)
				vmin[0]=_mm256_min_epi16(nb0[NB_N], nb0[NB_W]);
				vmax[0]=_mm256_max_epi16(nb0[NB_N], nb0[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb0[NB_N], nb0[NB_W]), nb0[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				vmin[1]=_mm256_min_epi16(nb2[NB_N], nb2[NB_W]);
				vmax[1]=_mm256_max_epi16(nb2[NB_N], nb2[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb2[NB_N], nb2[NB_W]), nb2[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				vmin[2]=_mm256_min_epi16(nb3[NB_N], nb3[NB_W]);
				vmax[2]=_mm256_max_epi16(nb3[NB_N], nb3[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb3[NB_N], nb3[NB_W]), nb3[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#ifdef ENABLE_MIXRCT
				vmin[8]=_mm256_min_epi16(nb5[NB_N], nb5[NB_W]);
				vmax[8]=_mm256_max_epi16(nb5[NB_N], nb5[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb5[NB_N], nb5[NB_W]), nb5[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[8]);
				pred=_mm256_min_epi16(pred, vmax[8]);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31
				);
				vmin[9]=_mm256_min_epi16(nb7[NB_N], nb7[NB_W]);
				vmax[9]=_mm256_max_epi16(nb7[NB_N], nb7[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb7[NB_N], nb7[NB_W]), nb7[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[9]);
				pred=_mm256_min_epi16(pred, vmax[9]);

				pred=_mm256_add_epi16(pred, nb6[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22
				);
				vmin[10]=_mm256_min_epi16(nb9[NB_N], nb9[NB_W]);
				vmax[10]=_mm256_max_epi16(nb9[NB_N], nb9[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb9[NB_N], nb9[NB_W]), nb9[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[10]);
				pred=_mm256_min_epi16(pred, vmax[10]);

				pred=_mm256_add_epi16(pred, nb8[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13
				);
#endif
#ifdef ENABLE_NONCAUSALRCT
				vmin[3]=_mm256_min_epi16(nby[NB_N], nby[NB_W]);
				vmax[3]=_mm256_max_epi16(nby[NB_N], nby[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nby[NB_N], nby[NB_W]), nby[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_sub_epi16(nby[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY
				);
				vmin[4]=_mm256_min_epi16(nbu[NB_N], nbu[NB_W]);
				vmax[4]=_mm256_max_epi16(nbu[NB_N], nbu[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb0[NB_N], nbu[NB_W]), nbu[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[4]);
				pred=_mm256_min_epi16(pred, vmax[4]);

				pred=_mm256_sub_epi16(nbu[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU
				);
				vmin[5]=_mm256_min_epi16(nbv[NB_N], nbv[NB_W]);
				vmax[5]=_mm256_max_epi16(nbv[NB_N], nbv[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nbv[NB_N], nbv[NB_W]), nbv[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[5]);
				pred=_mm256_min_epi16(pred, vmax[5]);

				pred=_mm256_sub_epi16(nbv[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV
				);
				vmin[6]=_mm256_min_epi16(nbu2[NB_N], nbu2[NB_W]);
				vmax[6]=_mm256_max_epi16(nbu2[NB_N], nbu2[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nbu2[NB_N], nbu2[NB_W]), nbu2[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[6]);
				pred=_mm256_min_epi16(pred, vmax[6]);

				pred=_mm256_sub_epi16(nbu2[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2
				);
				vmin[7]=_mm256_min_epi16(nbv2[NB_N], nbv2[NB_W]);
				vmax[7]=_mm256_max_epi16(nbv2[NB_N], nbv2[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nbv2[NB_N], nbv2[NB_W]), nbv2[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[7]);
				pred=_mm256_min_epi16(pred, vmax[7]);

				pred=_mm256_sub_epi16(nbv2[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2
				);
#endif
				
				//AV5 = clamp(W+((5*(N-NW)+NE-WW)>>3), N,W,NE)
				vmin[0]=_mm256_min_epi16(vmin[0], nb0[NB_NE]);
				vmax[0]=_mm256_max_epi16(vmin[0], nb0[NB_NE]);
				pred=_mm256_sub_epi16(nb0[NB_N], nb0[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb0[NB_NE], nb0[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb0[NB_W]);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				vmin[1]=_mm256_min_epi16(vmin[1], nb2[NB_NE]);
				vmax[1]=_mm256_max_epi16(vmin[1], nb2[NB_NE]);
				pred=_mm256_sub_epi16(nb2[NB_N], nb2[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb2[NB_NE], nb2[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb2[NB_W]);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				vmin[2]=_mm256_min_epi16(vmin[2], nb3[NB_NE]);
				vmax[2]=_mm256_max_epi16(vmin[2], nb3[NB_NE]);
				pred=_mm256_sub_epi16(nb3[NB_N], nb3[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb3[NB_NE], nb3[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb3[NB_W]);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#ifdef ENABLE_MIXRCT
				vmin[8]=_mm256_min_epi16(vmin[8], nb5[NB_NE]);
				vmax[8]=_mm256_max_epi16(vmin[8], nb5[NB_NE]);
				pred=_mm256_sub_epi16(nb5[NB_N], nb5[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb5[NB_NE], nb5[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb5[NB_W]);
				pred=_mm256_max_epi16(pred, vmin[8]);
				pred=_mm256_min_epi16(pred, vmax[8]);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31
				);
				vmin[9]=_mm256_min_epi16(vmin[9], nb7[NB_NE]);
				vmax[9]=_mm256_max_epi16(vmin[9], nb7[NB_NE]);
				pred=_mm256_sub_epi16(nb7[NB_N], nb7[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb7[NB_NE], nb7[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb7[NB_W]);
				pred=_mm256_max_epi16(pred, vmin[9]);
				pred=_mm256_min_epi16(pred, vmax[9]);

				pred=_mm256_add_epi16(pred, nb6[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22
				);
				vmin[10]=_mm256_min_epi16(vmin[10], nb9[NB_NE]);
				vmax[10]=_mm256_max_epi16(vmin[10], nb9[NB_NE]);
				pred=_mm256_sub_epi16(nb9[NB_N], nb9[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb9[NB_NE], nb9[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb9[NB_W]);
				pred=_mm256_max_epi16(pred, vmin[10]);
				pred=_mm256_min_epi16(pred, vmax[10]);

				pred=_mm256_add_epi16(pred, nb8[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13
				);
#endif
#ifdef ENABLE_NONCAUSALRCT
				vmin[3]=_mm256_min_epi16(vmin[3], nby[NB_NE]);
				vmax[3]=_mm256_max_epi16(vmin[3], nby[NB_NE]);
				pred=_mm256_sub_epi16(nby[NB_N], nby[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nby[NB_NE], nby[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nby[NB_W]);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_sub_epi16(nby[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY
				);
				vmin[4]=_mm256_min_epi16(vmin[4], nbu[NB_NE]);
				vmax[4]=_mm256_max_epi16(vmin[4], nbu[NB_NE]);
				pred=_mm256_sub_epi16(nbu[NB_N], nbu[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nbu[NB_NE], nbu[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nbu[NB_W]);
				pred=_mm256_max_epi16(pred, vmin[4]);
				pred=_mm256_min_epi16(pred, vmax[4]);

				pred=_mm256_sub_epi16(nbu[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU
				);
				vmin[5]=_mm256_min_epi16(vmin[5], nbv[NB_NE]);
				vmax[5]=_mm256_max_epi16(vmin[5], nbv[NB_NE]);
				pred=_mm256_sub_epi16(nbv[NB_N], nbv[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nbv[NB_NE], nbv[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nbv[NB_W]);
				pred=_mm256_max_epi16(pred, vmin[5]);
				pred=_mm256_min_epi16(pred, vmax[5]);

				pred=_mm256_sub_epi16(nbv[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV
				);
				vmin[6]=_mm256_min_epi16(vmin[6], nbu2[NB_NE]);
				vmax[6]=_mm256_max_epi16(vmin[6], nbu2[NB_NE]);
				pred=_mm256_sub_epi16(nbu2[NB_N], nbu2[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nbu2[NB_NE], nbu2[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nbu2[NB_W]);
				pred=_mm256_max_epi16(pred, vmin[6]);
				pred=_mm256_min_epi16(pred, vmax[6]);

				pred=_mm256_sub_epi16(nbu2[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2
				);
				vmin[7]=_mm256_min_epi16(vmin[7], nbv2[NB_NE]);
				vmax[7]=_mm256_max_epi16(vmin[7], nbv2[NB_NE]);
				pred=_mm256_sub_epi16(nbv2[NB_N], nbv2[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nbv2[NB_NE], nbv2[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nbv2[NB_W]);
				pred=_mm256_max_epi16(pred, vmin[7]);
				pred=_mm256_min_epi16(pred, vmax[7]);

				pred=_mm256_sub_epi16(nbv2[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2
				);
#endif
				
				//AV9
				//		1	-2	-1
				//	-1	-9	10	4
				//	-2	16	[?]>>4		clamp(N,W,NE)
				pred=_mm256_add_epi16(nb0[NB_N], _mm256_slli_epi16(nb0[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb0[NB_NN], nb0[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb0[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb0[NB_NW], 3), nb0[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb0[NB_NNW], _mm256_add_epi16(nb0[NB_NNE], nb0[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb0[NB_W], _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(nb2[NB_N], _mm256_slli_epi16(nb2[NB_N], 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb2[NB_NN], nb2[NB_WW]));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb2[NB_NE], 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb2[NB_NW], 3), nb2[NB_NW]));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb2[NB_NNW], _mm256_add_epi16(nb2[NB_NNE], nb2[NB_NWW])));
				pred=_mm256_add_epi16(nb2[NB_W], _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(nb3[NB_N], _mm256_slli_epi16(nb3[NB_N], 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb3[NB_NN], nb3[NB_WW]));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb3[NB_NE], 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb3[NB_NW], 3), nb3[NB_NW]));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb3[NB_NNW], _mm256_add_epi16(nb3[NB_NNE], nb3[NB_NWW])));
				pred=_mm256_add_epi16(nb3[NB_W], _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
#ifdef ENABLE_MIXRCT
				pred=_mm256_add_epi16(nb5[NB_N], _mm256_slli_epi16(nb5[NB_N], 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb5[NB_NN], nb5[NB_WW]));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb5[NB_NE], 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb5[NB_NW], 3), nb5[NB_NW]));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb5[NB_NNW], _mm256_add_epi16(nb5[NB_NNE], nb5[NB_NWW])));
				pred=_mm256_add_epi16(nb5[NB_W], _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin[8]);
				pred=_mm256_min_epi16(pred, vmax[8]);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31,
					OCH_R31, OCH_G31, OCH_B31
				);
				pred=_mm256_add_epi16(nb7[NB_N], _mm256_slli_epi16(nb7[NB_N], 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb7[NB_NN], nb7[NB_WW]));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb7[NB_NE], 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb7[NB_NW], 3), nb7[NB_NW]));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb7[NB_NNW], _mm256_add_epi16(nb7[NB_NNE], nb7[NB_NWW])));
				pred=_mm256_add_epi16(nb7[NB_W], _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin[9]);
				pred=_mm256_min_epi16(pred, vmax[9]);

				pred=_mm256_add_epi16(pred, nb6[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22,
					OCH_R22, OCH_G22, OCH_B22
				);
				pred=_mm256_add_epi16(nb9[NB_N], _mm256_slli_epi16(nb9[NB_N], 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb9[NB_NN], nb9[NB_WW]));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb9[NB_NE], 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb9[NB_NW], 3), nb9[NB_NW]));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb9[NB_NNW], _mm256_add_epi16(nb9[NB_NNE], nb9[NB_NWW])));
				pred=_mm256_add_epi16(nb9[NB_W], _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin[10]);
				pred=_mm256_min_epi16(pred, vmax[10]);

				pred=_mm256_add_epi16(pred, nb8[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13,
					OCH_R13, OCH_G13, OCH_B13
				);
#endif
#ifdef ENABLE_NONCAUSALRCT
				pred=_mm256_add_epi16(nby[NB_N], _mm256_slli_epi16(nby[NB_N], 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nby[NB_NN], nby[NB_WW]));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nby[NB_NE], 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nby[NB_NW], 3), nby[NB_NW]));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nby[NB_NNW], _mm256_add_epi16(nby[NB_NNE], nby[NB_NWW])));
				pred=_mm256_add_epi16(nby[NB_W], _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_sub_epi16(nby[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY,
					OCH_RGBY, OCH_GBRY, OCH_BRGY
				);
				pred=_mm256_add_epi16(nbu[NB_N], _mm256_slli_epi16(nbu[NB_N], 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nbu[NB_NN], nbu[NB_WW]));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nbu[NB_NE], 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nbu[NB_NW], 3), nbu[NB_NW]));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nbu[NB_NNW], _mm256_add_epi16(nbu[NB_NNE], nbu[NB_NWW])));
				pred=_mm256_add_epi16(nbu[NB_W], _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin[4]);
				pred=_mm256_min_epi16(pred, vmax[4]);

				pred=_mm256_sub_epi16(nbu[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU,
					OCH_RGBU, OCH_GBRU, OCH_BRGU
				);
				pred=_mm256_add_epi16(nbv[NB_N], _mm256_slli_epi16(nbv[NB_N], 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nbv[NB_NN], nbv[NB_WW]));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nbv[NB_NE], 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nbv[NB_NW], 3), nbv[NB_NW]));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nbv[NB_NNW], _mm256_add_epi16(nbv[NB_NNE], nbv[NB_NWW])));
				pred=_mm256_add_epi16(nbv[NB_W], _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin[5]);
				pred=_mm256_min_epi16(pred, vmax[5]);

				pred=_mm256_sub_epi16(nbv[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV,
					OCH_RGBV, OCH_GBRV, OCH_BRGV
				);
				pred=_mm256_add_epi16(nbu2[NB_N], _mm256_slli_epi16(nbu2[NB_N], 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nbu2[NB_NN], nbu2[NB_WW]));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nbu2[NB_NE], 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nbu2[NB_NW], 3), nbu2[NB_NW]));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nbu2[NB_NNW], _mm256_add_epi16(nbu2[NB_NNE], nbu2[NB_NWW])));
				pred=_mm256_add_epi16(nbu2[NB_W], _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin[6]);
				pred=_mm256_min_epi16(pred, vmax[6]);

				pred=_mm256_sub_epi16(nbu2[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2,
					OCH_RGBU2, OCH_GBRU2, OCH_BRGU2
				);
				pred=_mm256_add_epi16(nbv2[NB_N], _mm256_slli_epi16(nbv2[NB_N], 2));
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nbv2[NB_NN], nbv2[NB_WW]));
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nbv2[NB_NE], 1));
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nbv2[NB_NW], 3), nbv2[NB_NW]));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nbv2[NB_NNW], _mm256_add_epi16(nbv2[NB_NNE], nbv2[NB_NWW])));
				pred=_mm256_add_epi16(nbv2[NB_W], _mm256_srai_epi16(pred, 4));
				pred=_mm256_max_epi16(pred, vmin[7]);
				pred=_mm256_min_epi16(pred, vmax[7]);

				pred=_mm256_sub_epi16(nbv2[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2,
					OCH_RGBV2, OCH_GBRV2, OCH_BRGV2
				);
#endif
			}
		}
#ifdef ENABLE_SAD_ANALYSIS
		for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)
			csizes[kc]=(double)counters[kc];
#elif defined ENABLE_EXACT_ANALYSIS
		double gain=1./count;
		for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)
		{
			int *curr_hist=args->hist+((size_t)kc<<6);
			for(int ks=0;ks<64;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
					csizes[kc]-=freq*log2((double)freq*gain);
			}
			csizes[kc]+=bypasssizes[kc];
			csizes[kc]/=8;
		}
#elif defined ENABLE_ENTROPY_ANALYSIS
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
#endif
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
		const int nRCTs=RCT_COUNT;
	//	const int nRCTs=16;
		for(int kt=0;kt<nRCTs;++kt)//select best RCT
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
	}
	else
	{
		ac3_dec_init(&ec, args->decstart, args->decend);
		bestrct=ac3_dec_bypass_NPOT(&ec, RCT_COUNT);
		predidx[0]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[1]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[2]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		combination=rct_combinations[bestrct];
	}
	int yidx=0, uidx=1, vidx=2;
	int noncausalrct=0;
	switch(bestrct)
	{
	case RCT_R_G_B:
	case RCT_R_G_BG:
	case RCT_R_G_BR:
	case RCT_R_GR_BR:
	case RCT_R_GR_BG:
	case RCT_R_G_B31:
	case RCT_R_GR_B31:
	case RCT_R_G_B22:
	case RCT_R_GR_B22:
	case RCT_R_G_B13:
	case RCT_R_GR_B13:
		yidx=0;
		uidx=1;
		vidx=2;
		break;
	case RCT_G_B_RG:
	case RCT_G_B_RB:
	case RCT_G_BG_RG:
	case RCT_G_BG_RB:
	case RCT_G_B_R22:
	case RCT_G_BG_R22:
		yidx=1;
		uidx=2;
		vidx=0;
		break;
	case RCT_B_R_GR:
	case RCT_B_R_GB:
	case RCT_B_RB_GB:
	case RCT_B_RB_GR:
	case RCT_B_R_G31:
	case RCT_B_RB_G31:
	case RCT_B_R_G22:
	case RCT_B_RB_G22:
	case RCT_B_R_G13:
	case RCT_B_RB_G13:
		yidx=2;
		uidx=0;
		vidx=1;
		break;
	case RCT_G_RG_BR:
	case RCT_G_RG_B31:
	case RCT_G_RG_B22:
	case RCT_G_RG_B13:
		yidx=1;
		uidx=0;
		vidx=2;
		break;
	case RCT_B_GB_RG:
	case RCT_B_GB_R31:
	case RCT_B_GB_R22:
	case RCT_B_GB_R13:
		yidx=2;
		uidx=1;
		vidx=0;
		break;
	case RCT_R_BR_GB:
	case RCT_R_BR_G31:
	case RCT_R_BR_G22:
	case RCT_R_BR_G13:
		yidx=0;
		uidx=2;
		vidx=1;
		break;
	case RCT_RGB1:	noncausalrct=1; yidx=1; uidx=2; vidx=0; halfs[1]=halfs[2]=256;	break;
	case RCT_RGB2:	noncausalrct=2; yidx=1; uidx=2; vidx=0; halfs[1]=halfs[2]=256;	break;
	case RCT_RGB3:	noncausalrct=3; yidx=1; uidx=2; vidx=0; halfs[1]=halfs[2]=256;	break;
	case RCT_GBR1:	noncausalrct=1; yidx=2; uidx=0; vidx=1; halfs[1]=halfs[2]=256;	break;
	case RCT_GBR2:	noncausalrct=2; yidx=2; uidx=0; vidx=1; halfs[1]=halfs[2]=256;	break;
	case RCT_GBR3:	noncausalrct=3; yidx=2; uidx=0; vidx=1; halfs[1]=halfs[2]=256;	break;
	case RCT_BRG1:	noncausalrct=1; yidx=0; uidx=1; vidx=2; halfs[1]=halfs[2]=256;	break;
	case RCT_BRG2:	noncausalrct=2; yidx=0; uidx=1; vidx=2; halfs[1]=halfs[2]=256;	break;
	case RCT_BRG3:	noncausalrct=3; yidx=0; uidx=1; vidx=2; halfs[1]=halfs[2]=256;	break;
	}
	int yuv[4]={0};
	unsigned char yuvcoeff[]=
	{
		0, 0,
		combination[RII_UHELP], 0,
		combination[RII_VH0], combination[RII_VH1],
	};
	{
		static const int init_freqs[]={32, 8, 6, 4, 3, 2, 1};
		int sum=0;
		for(int ks=0;ks<args->tlevels;++ks)
		{
			int freq=init_freqs[MINVAR(ks, (int)_countof(init_freqs)-1)];
			sum+=args->hist[ks]=freq;
		}
		args->hist[args->tlevels]=sum;
	}
	memfill(args->hist+cdfstride, args->hist, args->histsize-cdfstride*sizeof(int), cdfstride*sizeof(int));
	memset(args->pixels, 0, args->bufsize);
	int paddedblockwidth=args->x2-args->x1+16;
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) short *rows[]=
		{
			args->pixels+(paddedblockwidth*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+(paddedblockwidth*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+(paddedblockwidth*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+(paddedblockwidth*((ky-3LL)&3)+8LL)*4*2,
		};
		int pred=0, error=0, sym=0;
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int idx=3*(args->iw*ky+kx);
			short
				*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NWW	=rows[1]-2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
#ifdef ENABLE_EDGECASES
			if(ky<=args->y1+1)
			{
				if(ky==args->y1)
					NE=NWW=NW=N=W;
				NNW=NW;
				NN=N;
				NNE=NE;
			}
			if(kx<=args->x1+1)
			{
				if(kx<=args->x1)
					NWW=NW=W=N;
				WW=W;
			}
			if(kx>=args->x2-1)
			{
				NNE=NN;
				NE=N;
			}
#endif
			if(args->fwd)
			{
				yuv[0]=args->src[idx+yidx]-128;
				yuv[1]=args->src[idx+uidx]-128;
				yuv[2]=args->src[idx+vidx]-128;
				if(noncausalrct)
				{
					yuv[2]-=yuv[0];
					yuv[1]-=yuv[0];
					yuv[0]+=(yuv[1]+yuv[2])>>2;
					if(noncausalrct==2)
						yuv[1]-=yuv[2]>>2;
					else if(noncausalrct==3)
						yuv[2]-=yuv[1]>>2;
					//yuv[2]=(char)(yuv[2]-yuv[0]);
					//yuv[1]=(char)(yuv[1]-yuv[0]);
					//yuv[0]=(char)(yuv[0]+((yuv[1]+yuv[2])>>2));
					//if(rct8==2)
					//	yuv[1]=(char)(yuv[1]-(yuv[2]>>2));
					//else if(rct8==3)
					//	yuv[2]=(char)(yuv[2]-(yuv[1]>>2));
				}
				//if(abs(yuv[0])>256||abs(yuv[1])>256||abs(yuv[2])>256)
				//	LOG_ERROR("impossible");
			}
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
			for(int kc=0;kc<3;++kc)
			{
				//if(ky==53&&kx==53&&kc==2)//
				//	printf("");
				int kc2=kc<<1;
				int offset=(yuvcoeff[kc<<1|0]*yuv[0]+yuvcoeff[kc<<1|1]*yuv[1])>>2;
				int vd, vx, vy;
				vd=abs(N[kc2]-W[kc2])+1;//1~512
				vx=abs(W[kc2]-WW[kc2])+abs(N[kc2]-NW[kc2])+abs(NE[kc2]-N  [kc2])+MAXVAR(WW[kc2+1], W[kc2+1])+W[kc2+1]+1;//1~2556
				vy=abs(N[kc2]-NN[kc2])+abs(W[kc2]-NW[kc2])+abs(NE[kc2]-NNE[kc2])+MAXVAR(NN[kc2+1], N[kc2+1])+N[kc2+1]+1;
				int qeD=FLOOR_LOG2(vd);//0~9
				int qeN=FLOOR_LOG2(vy);//0~11
				int qeW=FLOOR_LOG2(vx);
				//if((unsigned)qeN>(unsigned)(CLEVELS-2))//
				//	LOG_ERROR("");
				//if((unsigned)qeW>(unsigned)(CLEVELS-2))//
				//	LOG_ERROR("");
				//if((unsigned)qeD>(unsigned)(DLEVELS-2))//
				//	LOG_ERROR("");
				int cdf, freq=0, den;
#ifdef ENABLE_MIX4
				int *curr_hist[4];

				//if(qeN>CLEVELS-2||qeW>CLEVELS-2)//
				//	LOG_ERROR("");
				//UPDATE_MIN(qeN, CLEVELS-2);
				//UPDATE_MIN(qeW, CLEVELS-2);
				curr_hist[0]=args->hist+cdfstride*(CLEVELS*(CLEVELS*(DLEVELS*kc+qeD+0)+qeN+0)+qeW+0);
				curr_hist[1]=args->hist+cdfstride*(CLEVELS*(CLEVELS*(DLEVELS*kc+qeD+0)+qeN+0)+qeW+1);
				curr_hist[2]=args->hist+cdfstride*(CLEVELS*(CLEVELS*(DLEVELS*kc+qeD+0)+qeN+1)+qeW+0);
				curr_hist[3]=args->hist+cdfstride*(CLEVELS*(CLEVELS*(DLEVELS*kc+qeD+1)+qeN+0)+qeW+0);

			//	curr_hist[0]=args->hist+cdfstride*(CLEVELS*(CLEVELS*(DLEVELS*kc+qeD+0)+qeN+0)+qeW+0);
			//	curr_hist[1]=args->hist+cdfstride*(CLEVELS*(CLEVELS*(DLEVELS*kc+qeD+0)+qeN+0)+qeW+1);
			//	curr_hist[2]=args->hist+cdfstride*(CLEVELS*(CLEVELS*(DLEVELS*kc+qeD+0)+qeN+1)+qeW+0);
			//	curr_hist[3]=args->hist+cdfstride*(CLEVELS*(CLEVELS*(DLEVELS*kc+qeD+0)+qeN+1)+qeW+1);
			//	curr_hist[4]=args->hist+cdfstride*(CLEVELS*(CLEVELS*(DLEVELS*kc+qeD+1)+qeN+0)+qeW+0);
			//	curr_hist[5]=args->hist+cdfstride*(CLEVELS*(CLEVELS*(DLEVELS*kc+qeD+1)+qeN+0)+qeW+1);
			//	curr_hist[6]=args->hist+cdfstride*(CLEVELS*(CLEVELS*(DLEVELS*kc+qeD+1)+qeN+1)+qeW+0);
			//	curr_hist[7]=args->hist+cdfstride*(CLEVELS*(CLEVELS*(DLEVELS*kc+qeD+1)+qeN+1)+qeW+1);
#ifdef ENABLE_AV4
#define MIX4(X) (curr_hist[0][X]+curr_hist[1][X]+curr_hist[2][X]+curr_hist[3][X])
//#define MIX4(X) (curr_hist[0][X]+curr_hist[1][X]+curr_hist[2][X]+curr_hist[3][X]+curr_hist[4][X]+curr_hist[5][X]+curr_hist[6][X]+curr_hist[7][X])
#else
				int alphax=(int)((vx-(1LL<<qeW))<<MIXBITS>>qeW);
				int alphay=(int)((vy-(1LL<<qeN))<<MIXBITS>>qeN);
				//if(alphax<0||alphax>(1<<MIXBITS)||alphay<0||alphay>(1<<MIXBITS))
				//	LOG_ERROR("");
				//alphax=(((vx-(1<<qeW))<<MIXBITS)+(1<<qeW>>1))>>qeW;
				//alphay=(((vy-(1<<qeN))<<MIXBITS)+(1<<qeN>>1))>>qeN;
				//CLAMP2(alphax, 0, 1<<MIXBITS);
				//CLAMP2(alphay, 0, 1<<MIXBITS);
				//do
				//{
				//	ALIGN(16) int _dst[4];
				//	__m128i _mx=_mm_set_epi32(0, 0, alphay, alphax);
				//	__m128i _vmax=_mm_set_epi32(0, 0, 1<<14, 1<<14);
				//	_mx=_mm_max_epi32(_mx, _mm_setzero_si128());
				//	_mx=_mm_min_epi32(_mx, _vmax);
				//	_mm_store_si128((__m128i*)_dst, _mx);
				//	alphax=_dst[0];
				//	alphay=_dst[1];
				//}while(0);
				//CLAMP2_32(alphax, alphax, 0, 1<<MIXBITS);
				//CLAMP2_32(alphay, alphay, 0, 1<<MIXBITS);
#define MIX4(X) f28_mix4(curr_hist[0][X], curr_hist[1][X], curr_hist[2][X], curr_hist[3][X], alphax, alphay)
#endif
#else
				int *curr_hist=args->hist+cdfstride*(CLEVELS*CLEVELS*kc+CLEVELS*qeN+qeW);
			//	int *curr_hist=args->hist+cdfstride*(CLEVELS*CLEVELS*kc+CLEVELS*MINVAR(qeN, CLEVELS-1)+MINVAR(qeW, CLEVELS-1));
#define MIX4(X) curr_hist[X]
#endif
				den=MIX4(args->tlevels);
				switch(predidx[kc])
				{
				case PRED_AV2:
					pred=(N[kc2]+W[kc2])>>1;
					break;
				case PRED_CG:
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					break;
				case PRED_AV5:
					CLAMP3_32(pred,
						W[kc2]+((5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])>>3),
						N[kc2], W[kc2], NE[kc2]
					);
					break;
				case PRED_AV9:
					CLAMP3_32(pred,
						W[kc]+((10*N[kc]-9*NW[kc]+4*NE[kc]-2*(NN[kc]+WW[kc])+NNW[kc]-(NNE[kc]+NWW[kc]))>>4),
						N[kc], W[kc], NE[kc]
					);
					break;
				}
				pred+=offset;
				CLAMP2(pred, -halfs[kc], halfs[kc]-1);
				//CLAMP2_32(pred, pred, -128, 127);
				//if(ky==480&&kx==109&&kc==2)//
				//if(ky==147&&kx==317&&kc==1)//
				//	printf("");
				if(args->fwd)
				{
					error=yuv[kc]-pred;
					{
						int upred=halfs[kc]-abs(pred), aval=abs(error);
						if(aval<=upred)
						{
							sym=error;
#ifdef ENABLE_BIASCORR
							{
								int negmask=-((ibias_corr<0)&(sym!=-halfs[ch]));//sign is flipped if SSE correction was negative, to skew the histogram
								sym^=negmask;
								sym-=negmask;
							}
#endif
							sym=sym<<1^sym>>31;//pack sign
						}
						else
							sym=upred+aval;//error sign is known
					}
#ifdef ENABLE_QUANTIZATION
					quantize_pixel(sym, &token, &bypass, &nbits);
#else
					token=sym;
#endif
#ifdef _DEBUG
					if(token>=args->tlevels)
						LOG_ERROR("YXC %d %d %d  token %d/%d", ky, kx, kc, token, args->tlevels);
#endif
					cdf=0;
					for(int ks=0;ks<token;++ks)
						cdf+=MIX4(ks);
					freq=MIX4(token);
					ac3_enc_update_NPOT(&ec, cdf, freq, den);
#ifdef ENABLE_QUANTIZATION
					if(nbits)
						ac3_enc_bypass(&ec, bypass, nbits);
#endif
				}
				else
				{
					unsigned code=ac3_dec_getcdf_NPOT(&ec, den);
					cdf=0;
					token=0;
					for(;;)
					{
						unsigned cdf2;

						freq=MIX4(token);
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
#ifndef ENABLE_QUANTIZATION
					sym=token;
#else
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
#endif
					{
						int upred=halfs[kc]-abs(pred), negmask=0;
						if(sym<=(upred<<1))
						{
							error=sym>>1^-(sym&1);
#ifdef ENABLE_BIASCORR
							negmask=-((ibias_corr<0)&(error!=-half));
#endif
						}
						else
						{
							error=sym-upred;
							negmask=-(pred>0);
						}
						error^=negmask;
						error-=negmask;
					}
					yuv[kc]=error+pred;
					//if((unsigned)(yuv[kc]+256)>=512)//
					//	LOG_ERROR("");
				}
				curr[kc2+0]=yuv[kc]-offset;
				curr[kc2+1]=abs(error);
#ifdef ENABLE_MIX4
#ifdef ENABLE_AV4
				curr_hist[0][token]+=3; curr_hist[0][args->tlevels]+=3;
				++curr_hist[1][token]; ++curr_hist[1][args->tlevels];
				++curr_hist[2][token]; ++curr_hist[2][args->tlevels];
				++curr_hist[3][token]; ++curr_hist[3][args->tlevels];

			//	++curr_hist[0][token]; ++curr_hist[0][args->tlevels];
			//	++curr_hist[1][token]; ++curr_hist[1][args->tlevels];
			//	++curr_hist[2][token]; ++curr_hist[2][args->tlevels];
			//	++curr_hist[3][token]; ++curr_hist[3][args->tlevels];
			//	++curr_hist[4][token]; ++curr_hist[4][args->tlevels];
			//	++curr_hist[5][token]; ++curr_hist[5][args->tlevels];
			//	++curr_hist[6][token]; ++curr_hist[6][args->tlevels];
			//	++curr_hist[7][token]; ++curr_hist[7][args->tlevels];
#ifdef __GNUC__
#pragma GCC unroll 4
#endif
				for(int k=0;k<_countof(curr_hist);++k)
				{
					int *hist2=curr_hist[k];
					if(hist2[args->tlevels]>=8192)
					{
						int sum=0;
						for(int ks=0;ks<args->tlevels;++ks)
							sum+=hist2[ks]=(hist2[ks]+1)>>1;
						hist2[args->tlevels]=sum;
					}
				}
#else
				{
					int inc;
					inc=((1<<MIXBITS)-alphax)*((1<<MIXBITS)-alphay)>>(MIXBITS+MIXBITS-5); curr_hist[0][token]+=inc; curr_hist[0][args->tlevels]+=inc;
					inc=(             alphax)*((1<<MIXBITS)-alphay)>>(MIXBITS+MIXBITS-5); curr_hist[1][token]+=inc; curr_hist[1][args->tlevels]+=inc;
					inc=((1<<MIXBITS)-alphax)*(             alphay)>>(MIXBITS+MIXBITS-5); curr_hist[2][token]+=inc; curr_hist[2][args->tlevels]+=inc;
					inc=(             alphax)*(             alphay)>>(MIXBITS+MIXBITS-5); curr_hist[3][token]+=inc; curr_hist[3][args->tlevels]+=inc;
				}
				for(int kh=0;kh<4;++kh)
				{
					int *hist2=curr_hist[kh];
					if(hist2[args->tlevels]>=10752)//6144	4296	65536
					{
						int sum=0;
						for(int ks=0;ks<args->tlevels;++ks)
							sum+=hist2[ks]=(hist2[ks]+1)>>1;
						hist2[args->tlevels]=sum;
					}
				}
#endif
#else
				++curr_hist[token];
				++curr_hist[args->tlevels];
				if(curr_hist[args->tlevels]>=8192)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist[ks]=(curr_hist[ks]+1)>>1;
					curr_hist[args->tlevels]=sum;
				}
#endif
			}
			if(!args->fwd)
			{
				if(noncausalrct)
				{
					if(noncausalrct==2)
						yuv[1]+=yuv[2]>>2;
					else if(noncausalrct==3)
						yuv[2]+=yuv[1]>>2;
					yuv[0]-=(yuv[1]+yuv[2])>>2;
					yuv[1]+=yuv[0];
					yuv[2]+=yuv[0];
					//if(rct8==2)
					//	yuv[1]=(char)(yuv[1]+((yuv[2])>>2));
					//else if(rct8==3)
					//	yuv[2]=(char)(yuv[2]+((yuv[1])>>2));
					//yuv[0]=(char)(yuv[0]-((yuv[1]+yuv[2])>>2));
					//yuv[1]=(char)(yuv[1]+yuv[0]);
					//yuv[2]=(char)(yuv[2]+yuv[0]);
				}
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
	const int nch=3, depth=8;
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
		int nlevels=256*2;
#ifdef ENABLE_QUANTIZATION
		int token=0, bypass=0, nbits=0;

		quantize_pixel(nlevels, &token, &bypass, &nbits);
		tlevels=token+1;
#else
		tlevels=nlevels;
#endif
		statssize=CLEVELS*CLEVELS*DLEVELS*(tlevels+1)*nch*(int)sizeof(int);
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
				int bx, by;

				arg->blockidx=bx=kt+kt2;
				by=bx/xblocks;
				bx%=xblocks;
				arg->x1=BLOCKSIZE*bx;
				arg->y1=BLOCKSIZE*by;
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
						int blocksize=((arg->x2-arg->x1)*(arg->y2-arg->y1)*nch*depth+7)>>3;
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
			ptrdiff_t usize=((ptrdiff_t)iw*ih*nch*depth+7)>>3;
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
				compare_bufs_8(image2, src->data+headersize, iw, ih, nch, nch, "C01", 0, 1);
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
	if(test)
		pause();
	return 0;
}