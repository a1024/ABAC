#include"codec.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT

//	#define ENABLE_SELPRED
	#define ENABLE_SSE//good
//	#define ESTIMATE_SIZE

	#define ENABLE_WG

#define CODECNAME "C03"
#define AC3_PREC
#include"entropy.h"

#define BLOCKX 512
#define BLOCKY 256
#define MAXPRINTEDBLOCKS 20

#define OLS_STRIDE 31
#define A2_CTXBITS 8
#define A2_NCTX 12	//multiple of 4
#define A2_SSECBITS 6
#define A2_SSEPBITS 4
#define A2_SSECTR 16


#define OCHLIST\
	OCH(R)\
	OCH(G)\
	OCH(B)\
	OCH(RG)\
	OCH(RB)\
	OCH(GB)\
	OCH(GR)\
	OCH(BR)\
	OCH(BG)\
	OCH(R2)\
	OCH(G2)\
	OCH(B2)
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
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		3, 3, 3,	3, 3, 3,	0, 0, 0)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		3, 3, 1,	3, 3, 3,	0, 0, 0)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		3, 3, 0,	3, 3, 3,	0, 0, 0)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		3, 3, 0,	3, 3, 3,	0, 0, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		3, 3, 1,	3, 3, 3,	0, 0, 0)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		3, 3, 1,	3, 3, 3,	0, 0, 0)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		3, 3, 0,	3, 3, 3,	0, 0, 0)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		3, 0, 0,	3, 3, 3,	0, 0, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		3, 0, 0,	3, 3, 3,	0, 0, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		3, 0, 0,	3, 3, 3,	0, 0, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(R_G_B2,	OCH_R,		OCH_G,		OCH_B2,		3, 3, 0,	3, 3, 1,	0, 0, 1)\
	RCT(R_B_G2,	OCH_R,		OCH_B,		OCH_G2,		3, 3, 0,	3, 3, 1,	0, 0, 1)\
	RCT(G_B_R2,	OCH_G,		OCH_B,		OCH_R2,		3, 3, 0,	3, 3, 1,	0, 0, 1)\
	RCT(R_GR_B2,	OCH_R,		OCH_GR,		OCH_B2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(R_BR_G2,	OCH_R,		OCH_BR,		OCH_G2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(G_BG_R2,	OCH_G,		OCH_BG,		OCH_R2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(G_RG_B2,	OCH_G,		OCH_RG,		OCH_B2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(B_RB_G2,	OCH_B,		OCH_RB,		OCH_G2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(B_GB_R2,	OCH_B,		OCH_GB,		OCH_R2,		3, 0, 0,	3, 3, 1,	0, 0, 1)
typedef enum _RCTIndex
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2, YSH2, USH2, VSH2) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][12]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2, YSH2, USH2, VSH2)\
	{YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2, YSH2, USH2, VSH2},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2, YSH2, USH2, VSH2) #LABEL,
	RCTLIST
#undef  RCT
};

#ifdef ENABLE_SELPRED
#define PREDLIST\
	PRED(W)\
	PRED(CG)\
	PRED(AV5)\
	PRED(AV9)\
	PRED(AV12)
#else
#define PREDLIST\
	PRED(CG)
#endif
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
#ifdef ENABLE_SELPRED
	NB_NNWW,	NB_NNW,		NB_NN,		NB_NNE,		NB_NNEE,
	NB_NWW,		NB_NW,		NB_N,		NB_NE,		NB_NEE,
	NB_WW,		NB_W,		NB_curr,
#else
	NB_NW,		NB_N,
	NB_W,		NB_curr,
#endif
	NB_COUNT,
} NBIndex;
static const short av12_icoeffs[12]=
{
	 0x04,	 0x03,	-0x1F,	-0x26,	 0x00,
	 0x07,	-0x9E,	 0xDB,	 0x1E,	 0x13,
	-0x2A,	 0xF3,
};


//WG:
#define WG_DECAY_NUM	493
#define WG_DECAY_SH	9

#define WG_NPREDS	16	//multiple of 4
#define WG_PREDLIST\
	WG_PRED(255, cgrad)\
	WG_PRED(50, wgrad)\
	WG_PRED(50, 3*N-W-NE)\
	WG_PRED(30, 2*NE-W)\
	WG_PRED(132, N+W-NW)\
	WG_PRED(100, N+W-NW+((eN+eW-eNW+16)>>5))\
	WG_PRED(150, N+eN)\
	WG_PRED(100, N)\
	WG_PRED(150, W+eW)\
	WG_PRED(100, W)\
	WG_PRED(165, W+NE-N)\
	WG_PRED(220, N+NE-NNE)\
	WG_PRED(25, N+NE-NNE+((eN+eNE-eNNE)>>2))\
	WG_PRED(165, 3*(N-NN)+NNN)\
	WG_PRED(165, 3*(W-WW)+WWW)\
	WG_PRED(150, (W+NEEE)/2)
static void wg_init(double *weights)
{
	int j=0;
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
	WG_PREDLIST
#undef  WG_PRED
}
FORCEINLINE int wg_predict(
	const double *weights,
	short **rows, const int stride, int kc2,
	int cgrad, const int *perrors, const int *NWerrors, const int *Nerrors, const int *NEerrors, int *preds
)
{
	int j=0;
	short
		NNNWWWW	=rows[3][kc2-4*stride+0],
		NNWWWW	=rows[2][kc2-4*stride+0],
		NNN	=rows[3][kc2+0*stride+0],
		NN	=rows[2][kc2+0*stride+0],
		NNE	=rows[2][kc2+1*stride+0],
		NNEEE	=rows[2][kc2+3*stride+0],
		NW	=rows[1][kc2-1*stride+0],
		N	=rows[1][kc2+0*stride+0],
		NE	=rows[1][kc2+1*stride+0],
		NEE	=rows[1][kc2+2*stride+0],
		NEEE	=rows[1][kc2+3*stride+0],
		WWW	=rows[0][kc2-3*stride+0],
		WW	=rows[0][kc2-2*stride+0],
		W	=rows[0][kc2-1*stride+0],
		eNNE	=rows[2][kc2+1*stride+1],
		eNW	=rows[1][kc2-1*stride+1],
		eN	=rows[1][kc2+0*stride+1],
		eNE	=rows[1][kc2+1*stride+1],
		eW	=rows[0][kc2-1*stride+1];
	int gy=abs(eN)+1, gx=abs(eW)+1;
	int wgrad=(N*gy+W*gx)/(gy+gx);
	(void)NNNWWWW;
	(void)NNWWWW;
	(void)NNEEE;
	(void)NEE;
	MEDIAN3_32(cgrad, N, W, N+W-NW);
	
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST
#undef  WG_PRED
	
#if 1
	__m128i one=_mm_set1_epi32(1);
	__m128i me0=_mm_load_si128((__m128i*)perrors+0);
	__m128i me1=_mm_load_si128((__m128i*)perrors+1);
	__m128i me2=_mm_load_si128((__m128i*)perrors+2);
	__m128i me3=_mm_load_si128((__m128i*)perrors+3);
	__m256d w0=_mm256_load_pd(weights+0*4);
	__m256d w1=_mm256_load_pd(weights+1*4);
	__m256d w2=_mm256_load_pd(weights+2*4);
	__m256d w3=_mm256_load_pd(weights+3*4);
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)NWerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)NWerrors+1));
	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)NWerrors+2));
	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)NWerrors+3));
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)Nerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)Nerrors+1));
	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)Nerrors+2));
	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)Nerrors+3));
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)NEerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)NEerrors+1));
	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)NEerrors+2));
	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)NEerrors+3));
	me0=_mm_add_epi32(me0, one);
	me1=_mm_add_epi32(me1, one);
	me2=_mm_add_epi32(me2, one);
	me3=_mm_add_epi32(me3, one);
	__m256d de0=_mm256_cvtepi32_pd(me0);
	__m256d de1=_mm256_cvtepi32_pd(me1);
	__m256d de2=_mm256_cvtepi32_pd(me2);
	__m256d de3=_mm256_cvtepi32_pd(me3);
	w0=_mm256_div_pd(w0, de0);
	w1=_mm256_div_pd(w1, de1);
	w2=_mm256_div_pd(w2, de2);
	w3=_mm256_div_pd(w3, de3);
	de0=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+0));
	de1=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+1));
	de2=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+2));
	de3=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+3));
	de0=_mm256_mul_pd(de0, w0);
	de1=_mm256_mul_pd(de1, w1);
	de2=_mm256_mul_pd(de2, w2);
	de3=_mm256_mul_pd(de3, w3);
	w0=_mm256_add_pd(w0, w1);
	w0=_mm256_add_pd(w0, w2);
	w0=_mm256_add_pd(w0, w3);
	de0=_mm256_add_pd(de0, de1);
	de0=_mm256_add_pd(de0, de2);
	de0=_mm256_add_pd(de0, de3);
	//[num3 num2 num1 num0]
	//[den3 den2 den1 den0]
	//r = hadd(num, den) = [den3+den2 num3+num2 den1+den0 num1+num0]
	//lo=_mm256_extractf128_pd(r, 0)
	//hi=_mm256_extractf128_pd(r, 1)
	//hi+lo = [den3+den2+den1+den0 num3+num2+num1+num0]
	w0=_mm256_hadd_pd(de0, w0);
	__m128d dp=_mm_add_pd(_mm256_extractf128_pd(w0, 1), _mm256_extractf128_pd(w0, 0));
	dp=_mm_div_pd(dp, _mm_permute_pd(dp, 3));
	__m128i mp=_mm_cvtpd_epi32(dp);
	__m128i mN	=_mm_set_epi32(0, 0, 0, N);
	__m128i mW	=_mm_set_epi32(0, 0, 0, W);
	__m128i mNE	=_mm_set_epi32(0, 0, 0, NE);
	__m128i vmin=_mm_min_epi32(mN, mW);
	__m128i vmax=_mm_max_epi32(mN, mW);
	vmin=_mm_min_epi32(vmin, mNE);
	vmax=_mm_max_epi32(vmax, mNE);
	mp=_mm_max_epi32(mp, vmin);
	mp=_mm_min_epi32(mp, vmax);
	{
		ALIGN(16) int pred[4];
		_mm_store_si128((__m128i*)pred, mp);
		return pred[0];
	}
#else
	{
		double pred2=0, wsum=0;
		int pred;
#ifdef __GNUC__
#pragma GCC unroll 16
#endif
		for(int k=0;k<WG_NPREDS;++k)
		{
			double weight=(double)weights[k]/(perrors[k]+1);
			pred2+=weight*preds[k];
			wsum+=weight;
		}
		pred2/=wsum;
#ifdef _MSC_VER
		pred=_cvt_dtoi_fast(pred2);
#else
		pred=(int)pred2;
#endif
		CLAMP3_32(pred, pred, N, W, NE);
		return pred;
	}
#endif
}
FORCEINLINE void wg_update(int curr, const int *preds, int *perrors, int *Werrors, int *currerrors, int *NEerrors)
{
#ifdef __GNUC__
#pragma GCC unroll 16
#endif
	for(int k=0;k<WG_NPREDS;++k)
	{
		int e2=abs(curr-preds[k]);
		perrors[k]=(perrors[k]+e2)*WG_DECAY_NUM>>WG_DECAY_SH;
		e2<<=3;
		currerrors[k]=(2*Werrors[k]+e2+NEerrors[k])>>2;
		NEerrors[k]+=e2;
	}
}


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
		int lgv=FLOOR_LOG2((unsigned)val);
		int mantissa=val-(1<<lgv);
		*token = (1<<CONFIG_EXP) + (
				(lgv-CONFIG_EXP)<<(CONFIG_MSB+CONFIG_LSB)|
				(mantissa>>(lgv-CONFIG_MSB))<<CONFIG_LSB|
				(mantissa&((1<<CONFIG_LSB)-1))
			);
		*nbits=lgv-(CONFIG_MSB+CONFIG_LSB);
		*bypass=val>>CONFIG_LSB&((1LL<<*nbits)-1);
	}
}
typedef struct _ThreadArgs
{
	const unsigned char *src;
	unsigned char *dst;
	int iw, ih;

	int fwd, test, loud, b1, b2, xblocks, blocksperthread, currentblock, x1, x2, y1, y2;
	int bufsize, ebufsize, histsize;
	short *pixels;
	int *ebuf;
	int *hist;

	BList *lists;
	const int *offsets;
	const unsigned char *decsrc, *decstart, *decend;
	
	//int hist[54<<8];

	unsigned short stats[A2_NCTX*3<<A2_CTXBITS<<8];
#ifdef ENABLE_SSE
	long long sse[3][1<<A2_SSECBITS][1<<A2_SSECBITS][8<<A2_SSEPBITS];
	//long long sse[24<<A2_SSEPBITS];
#endif

	//aux
	int blockidx;
	double bestsize;
	int bestrct, predidx[3];

#ifdef ESTIMATE_SIZE
	int hist2[3][256];
#endif
} ThreadArgs;
static void block_thread(void *param)
{
	const int nch=3;
	ThreadArgs *args=(ThreadArgs*)param;
	AC3 ec;
	const unsigned char *image=args->fwd?args->src:args->dst;
	unsigned char bestrct=0, combination[6]={0}, predidx[4]={0};

	if(args->fwd)
	{
		int ystride=args->iw*3;
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		unsigned char predsel[OCH_COUNT]={0};
		int res=(args->x2-args->x1-3)/5*5*(args->y2-args->y1-2);
		__m256i av12_mcoeffs[12];
		
		memset(args->hist, 0, args->histsize);
		for(int k=0;k<(int)_countof(av12_mcoeffs);++k)
			av12_mcoeffs[k]=_mm256_set1_epi16(av12_icoeffs[k]>>1);
		for(int ky=args->y1+2;ky<args->y2;++ky)//analysis loop
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
			ALIGN(32) short result[16]={0};
			for(;kx<args->x2-5;kx+=5, ptr+=15)
			{
				__m256i
					nb0[NB_COUNT],//rgb
					nb1[NB_COUNT],//gbr
					nb2[NB_COUNT],//rgb - gbr
					nb3[NB_COUNT],//gbr - rgb
					nb4[NB_COUNT],//(gbr+brg)/2
					nb5[NB_COUNT];//rgb - (gbr+brg)/2
				__m256i vmin[4], vmax[4], pred;
				{
					__m128i nb8[NB_COUNT]=//8-bit
					{
#ifdef ESTIMATE_SIZE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride-2*3+0)), half8),//NNWW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride-1*3+0)), half8),//NNW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride+0*3+0)), half8),//NN
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride+1*3+0)), half8),//NNE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride+2*3+0)), half8),//NNEE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride-2*3+0)), half8),//NWW
#endif
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride-1*3+0)), half8),//NW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+0*3+0)), half8),//N
#ifdef ESTIMATE_SIZE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+1*3+0)), half8),//NE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+2*3+0)), half8),//NEE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride-2*3+0)), half8),//WW
#endif
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride-1*3+0)), half8),//W
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride+0*3+0)), half8),//curr
					};
#ifdef __GNUC__
#pragma GCC unroll 4
#endif
					for(int k=0;k<_countof(nb8);++k)
					{
						__m128i temp;
						__m256i t2;
						nb0[k]=_mm256_cvtepi8_epi16(nb8[k]);
						temp=_mm_shuffle_epi8(nb8[k], shuf);
						nb1[k]=_mm256_cvtepi8_epi16(temp);
						t2=_mm256_cvtepi8_epi16(_mm_shuffle_epi8(temp, shuf));
						t2=_mm256_add_epi16(t2, nb1[k]);
						t2=_mm256_srai_epi16(t2, 1);
						nb2[k]=_mm256_sub_epi16(nb0[k], nb1[k]);
						nb3[k]=_mm256_sub_epi16(nb1[k], nb0[k]);
						nb4[k]=t2;
						nb5[k]=_mm256_sub_epi16(nb0[k], t2);
					}
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
#ifdef ENABLE_SELPRED
				//W
				pred=_mm256_sub_epi16(nb0[NB_curr], nb0[NB_W]);
				UPDATE(
					PRED_W,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);

				pred=_mm256_add_epi16(nb2[NB_W], nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_W,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);

				pred=_mm256_add_epi16(nb3[NB_W], nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_W,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);

				pred=_mm256_add_epi16(nb5[NB_W], nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_W,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);
#endif

				//CG
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
				vmin[3]=_mm256_min_epi16(nb5[NB_N], nb5[NB_W]);
				vmax[3]=_mm256_max_epi16(nb5[NB_N], nb5[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb5[NB_N], nb5[NB_W]), nb5[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);
				
#ifdef ENABLE_SELPRED
				//AV5
				//		-5	5	1
				//	-1	8	[?]>>3
				pred=_mm256_sub_epi16(nb0[NB_N], nb0[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb0[NB_NE], nb0[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb0[NB_W]);
				//vmin=_mm256_min_epi16(N, W);
				//vmax=_mm256_max_epi16(N, W);
				vmin[0]=_mm256_min_epi16(vmin[0], nb0[NB_NE]);
				vmax[0]=_mm256_max_epi16(vmax[0], nb0[NB_NE]);
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
				pred=_mm256_sub_epi16(nb2[NB_N], nb2[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb2[NB_NE], nb2[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb2[NB_W]);
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				vmin[1]=_mm256_min_epi16(vmin[1], nb2[NB_NE]);
				vmax[1]=_mm256_max_epi16(vmax[1], nb2[NB_NE]);
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
				pred=_mm256_sub_epi16(nb3[NB_N], nb3[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb3[NB_NE], nb3[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb3[NB_W]);
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				vmin[2]=_mm256_min_epi16(vmin[2], nb3[NB_NE]);
				vmax[2]=_mm256_max_epi16(vmax[2], nb3[NB_NE]);
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
				pred=_mm256_sub_epi16(nb5[NB_N], nb5[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb5[NB_NE], nb5[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb5[NB_W]);
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				vmin[3]=_mm256_min_epi16(vmin[1], nb5[NB_NE]);
				vmax[3]=_mm256_max_epi16(vmax[1], nb5[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//AV9
				//		1	-2	-1
				//	-1	-9	10	4
				//	-2	16	[?]>>4
				pred=_mm256_add_epi16(nb0[NB_N], _mm256_slli_epi16(nb0[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb0[NB_NN], nb0[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb0[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb0[NB_NW], 3), nb0[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb0[NB_NNW], _mm256_add_epi16(nb0[NB_NNE], nb0[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb0[NB_W], _mm256_srai_epi16(pred, 4));
				//vmin=_mm256_min_epi16(N, W);
				//vmax=_mm256_max_epi16(N, W);
				//vmin=_mm256_min_epi16(vmin, NE);
				//vmax=_mm256_max_epi16(vmax, NE);
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
				pred=_mm256_add_epi16(nb2[NB_N], _mm256_slli_epi16(nb2[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb2[NB_NN], nb2[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb2[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb2[NB_NW], 3), nb2[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb2[NB_NNW], _mm256_add_epi16(nb2[NB_NNE], nb2[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb2[NB_W], _mm256_srai_epi16(pred, 4));
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin=_mm256_min_epi16(vmin, NE3);
				//vmax=_mm256_max_epi16(vmax, NE3);
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
				pred=_mm256_add_epi16(nb3[NB_N], _mm256_slli_epi16(nb3[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb3[NB_NN], nb3[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb3[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb3[NB_NW], 3), nb3[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb3[NB_NNW], _mm256_add_epi16(nb3[NB_NNE], nb3[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb3[NB_W], _mm256_srai_epi16(pred, 4));
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin=_mm256_min_epi16(vmin, NE3);
				//vmax=_mm256_max_epi16(vmax, NE3);
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
				pred=_mm256_add_epi16(nb5[NB_N], _mm256_slli_epi16(nb5[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb5[NB_NN], nb5[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb5[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb5[NB_NW], 3), nb5[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb5[NB_NNW], _mm256_add_epi16(nb5[NB_NNE], nb5[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb5[NB_W], _mm256_srai_epi16(pred, 4));
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin=_mm256_min_epi16(vmin, NE3);
				//vmax=_mm256_max_epi16(vmax, NE3);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);
					
				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//AV12
				pred=_mm256_setzero_si256();
				for(int k=0;k<(int)_countof(av12_icoeffs);++k)
					pred=_mm256_add_epi16(pred, _mm256_mullo_epi16(av12_mcoeffs[k], nb0[k]));
				pred=_mm256_srai_epi16(pred, 7);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV12,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_setzero_si256();
				for(int k=0;k<(int)_countof(av12_icoeffs);++k)
					pred=_mm256_add_epi16(pred, _mm256_mullo_epi16(av12_mcoeffs[k], nb2[k]));
				pred=_mm256_srai_epi16(pred, 7);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);
					
				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV12,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_setzero_si256();
				for(int k=0;k<(int)_countof(av12_icoeffs);++k)
					pred=_mm256_add_epi16(pred, _mm256_mullo_epi16(av12_mcoeffs[k], nb3[k]));
				pred=_mm256_srai_epi16(pred, 7);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);
					
				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV12,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_setzero_si256();
				for(int k=0;k<(int)_countof(av12_icoeffs);++k)
					pred=_mm256_add_epi16(pred, _mm256_mullo_epi16(av12_mcoeffs[k], nb5[k]));
				pred=_mm256_srai_epi16(pred, 7);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);
					
				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV12,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);
#endif
			}
		}
		for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)
		{
			int *curr_hist=args->hist+((size_t)kc<<8);
			for(int ks=0;ks<256;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
					csizes[kc]-=freq*log2((double)freq/res);
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
		memcpy(combination, rct_combinations[bestrct], sizeof(combination));
		predidx[0]=predsel[combination[0]];
		predidx[1]=predsel[combination[1]];
		predidx[2]=predsel[combination[2]];
		args->bestsize+=bestsize;
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
		blist_init(args->lists+args->currentblock);
		ac3_enc_init(&ec, args->lists+args->currentblock);
		ac3_enc_bypass_NPOT(&ec, bestrct, RCT_COUNT);
#ifdef ENABLE_SELPRED
		ac3_enc_bypass_NPOT(&ec, predidx[0], PRED_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[1], PRED_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[2], PRED_COUNT);
#endif
	}
	else
	{
		ac3_dec_init(&ec, args->decstart, args->decend);
		bestrct=ac3_dec_bypass_NPOT(&ec, RCT_COUNT);
#ifdef ENABLE_SELPRED
		predidx[0]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[1]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[2]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
#endif
	}
	ALIGN(32) long long mixer[8*3*A2_NCTX]={0};
	FILLMEM(mixer, 0x3000, sizeof(mixer), sizeof(long long));
	FILLMEM(args->stats, 0x8000, sizeof(args->stats), sizeof(short));
#ifdef ENABLE_SSE
	memset(args->sse, 0, sizeof(args->sse));
#endif
	
#ifdef ESTIMATE_SIZE
	memset(args->hist2, 0, sizeof(args->hist2));
#endif
#ifdef ENABLE_WG
	ALIGN(32) double wg_weights[WG_NPREDS*3]={0};
	ALIGN(32) int wg_perrors[WG_NPREDS*3]={0}, wg_preds[WG_NPREDS]={0};
	wg_init(wg_weights+WG_NPREDS*0);
	wg_init(wg_weights+WG_NPREDS*1);
	wg_init(wg_weights+WG_NPREDS*2);
#endif
	memset(args->ebuf, 0, args->ebufsize);
	memset(args->pixels, 0, args->bufsize);
#ifdef ABAC_SIMPLEOVERRIDE
	unsigned short stats2[768]={0};
	FILLMEM(stats2, 0x8000, sizeof(short[768]), sizeof(short));
#endif
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(32) short *rows[]=
		{
			args->pixels+((BLOCKX+16LL)*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKX+16LL)*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKX+16LL)*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKX+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		ALIGN(32) int *erows[]=
		{
			args->ebuf+((BLOCKX+16LL)*((ky-0LL)&3)+8LL)*4*WG_NPREDS,
			args->ebuf+((BLOCKX+16LL)*((ky-1LL)&3)+8LL)*4*WG_NPREDS,
			args->ebuf+((BLOCKX+16LL)*((ky-2LL)&3)+8LL)*4*WG_NPREDS,
			args->ebuf+((BLOCKX+16LL)*((ky-3LL)&3)+8LL)*4*WG_NPREDS,
		};
		int yuv[4]={0};
		int pred=0, error=0;
		const unsigned char *combination=rct_combinations[bestrct];
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int idx=nch*(args->iw*ky+kx);
			short
				*NNN	=rows[3]+0*4*2,
				*NNWW	=rows[2]-2*4*2,
				*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NNEE	=rows[2]+2*4*2,
				*NNEEE	=rows[2]+3*4*2,
				*NWWW	=rows[1]-3*4*2,
				*NWW	=rows[1]-2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*NEE	=rows[1]+2*4*2,
				*NEEE	=rows[1]+3*4*2,
				*WWWW	=rows[0]-4*4*2,
				*WWW	=rows[0]-3*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
			(void)NNN;
			(void)NNWW;
			(void)NNW;
			(void)NN;
			(void)NNE;
			(void)NNEE;
			(void)NNEEE;
			(void)NWWW;
			(void)NWW;
			(void)NW;
			(void)N;
			(void)NE;
			(void)NEE;
			(void)NEEE;
			(void)WWWW;
			(void)WWW;
			(void)WW;
			(void)W;
#if 0
			if(ky<=args->y1+2)
			{
				if(ky<=args->y1+1)
				{
					if(ky==args->y1)
						NEEE=NEE=NE=NWWW=NWW=NW=N=W;
					NNWW=NWW;
					NNW=NW;
					NN=N;
					NNE=NE;
					NNEE=NEE;
					NNEEE=NEEE;
				}
				NNN=NN;
			}
			if(kx<=args->x1+3)
			{
				if(kx<=args->x1+2)
				{
					if(kx<=args->x1+1)
					{
						if(kx<=args->x1)
							NW=W=N;
						WW=W;
						NWW=NW;
					}
					WWW=WW;
				}
				WWWW=WWW;
			}
			if(kx>=args->x2-3)
			{
				if(kx>=args->x2-2)
				{
					if(kx>=args->x2-1)
					{
						NNE=NN;
						NE=N;
					}
					NNEE=NNE;
					NEE=NE;
				}
				NEEE=NEE;
			}
#endif
			if(args->fwd)
			{
				switch(bestrct)
				{
				case RCT_R_G_B:
				case RCT_R_G_BG:
				case RCT_R_G_BR:
				case RCT_R_GR_BR:
				case RCT_R_GR_BG:
				case RCT_R_G_B2:
				case RCT_R_GR_B2:
					yuv[0]=args->src[idx+0]-128;
					yuv[1]=args->src[idx+1]-128;
					yuv[2]=args->src[idx+2]-128;
					break;
				case RCT_G_B_RG:
				case RCT_G_B_RB:
				case RCT_G_BG_RG:
				case RCT_G_BG_RB:
				case RCT_G_B_R2:
				case RCT_G_BG_R2:
					yuv[0]=args->src[idx+1]-128;
					yuv[1]=args->src[idx+2]-128;
					yuv[2]=args->src[idx+0]-128;
					break;
				case RCT_B_R_GR:
				case RCT_B_R_GB:
				case RCT_B_RB_GB:
				case RCT_B_RB_GR:
				case RCT_B_RB_G2:
					yuv[0]=args->src[idx+2]-128;
					yuv[1]=args->src[idx+0]-128;
					yuv[2]=args->src[idx+1]-128;
					break;
				case RCT_G_RG_BR:
				case RCT_G_RG_B2:
					yuv[0]=args->src[idx+1]-128;
					yuv[1]=args->src[idx+0]-128;
					yuv[2]=args->src[idx+2]-128;
					break;
				case RCT_B_GB_RG:
				case RCT_B_GB_R2:
					yuv[0]=args->src[idx+2]-128;
					yuv[1]=args->src[idx+1]-128;
					yuv[2]=args->src[idx+0]-128;
					break;
				case RCT_R_BR_GB:
				case RCT_R_B_G2:
				case RCT_R_BR_G2:
					yuv[0]=args->src[idx+0]-128;
					yuv[1]=args->src[idx+2]-128;
					yuv[2]=args->src[idx+1]-128;
					break;
				}
			}

			for(int kc=0;kc<nch;++kc)
			{
				int kc2=kc<<1, kc3=kc*WG_NPREDS;
				int offset=(yuv[combination[kc+3]]+yuv[combination[kc+6]])>>combination[kc+9];
				int
					*eNW	=erows[1]+kc3-1*4*WG_NPREDS,
					*eN	=erows[1]+kc3+0*4*WG_NPREDS,
					*eNE	=erows[1]+kc3+1*4*WG_NPREDS,
					*eW	=erows[0]+kc3-1*4*WG_NPREDS,
					*ecurr	=erows[0]+kc3+0*4*WG_NPREDS;
				//if(ky==10&&kx==10&&kc==0)//
				//	printf("");
#ifdef ENABLE_SELPRED
				switch(predidx[kc])
				{
				case PRED_W:
					pred=W[kc2];
					break;
				case PRED_CG:
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					break;
				case PRED_AV5:
					CLAMP3_32(pred, W[kc2]+((5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])>>3), N[kc2], W[kc2], NE[kc2]);
					break;
				case PRED_AV9:
					CLAMP3_32(pred,
						W[kc2]+((10*N[kc2]-9*NW[kc2]+4*NE[kc2]-2*(NN[kc2]+WW[kc2])+NNW[kc2]-(NNE[kc2]+NWW[kc2]))>>4),
						N[kc2], W[kc2], NE[kc2]
					);
					break;
				case PRED_AV12:
					pred=(
						av12_icoeffs[ 0]*NNWW[kc2]+
						av12_icoeffs[ 1]*NNW[kc2]+
						av12_icoeffs[ 2]*NN[kc2]+
						av12_icoeffs[ 3]*NNE[kc2]+
						av12_icoeffs[ 4]*NNEE[kc2]+
						av12_icoeffs[ 5]*NWW[kc2]+
						av12_icoeffs[ 6]*NW[kc2]+
						av12_icoeffs[ 7]*N[kc2]+
						av12_icoeffs[ 8]*NE[kc2]+
						av12_icoeffs[ 9]*NEE[kc2]+
						av12_icoeffs[10]*WW[kc2]+
						av12_icoeffs[11]*W[kc2]
					)>>8;
					CLAMP3_32(pred, pred, N[kc2], W[kc2], NE[kc2]);
					break;
				}
#ifdef ENABLE_WG
				pred=wg_predict(wg_weights+WG_NPREDS*kc, rows, 4*2, kc2, pred, wg_perrors+WG_NPREDS*kc, wg_preds);
#endif
#elif defined ENABLE_WG
				pred=wg_predict(wg_weights+WG_NPREDS*kc, rows, 4*2, kc2, 0, wg_perrors+WG_NPREDS*kc, eNW, eN, eNE, wg_preds);
#else
				MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
#endif
				pred+=offset;
				CLAMP2(pred, -128, 127);
				unsigned short *curr_stats[A2_NCTX];
				int grad=abs(NE[kc2+0]-N[kc2+0])+abs(N[kc2+0]-NW[kc2+0])+abs(W[kc2+0]-NW[kc2+0]);
				int ctx[A2_NCTX]=
				{
					pred*3,
					//pred*5+abs(W[kc2+1])/25,
					pred/5,
					//(W[kc2+0]+NE[kc2+0]+NEE[kc2+0]+NEEE[kc2+0])>>1,
					//N[kc2+0],
					grad*2,
					//grad<<2|(ky&1)<<1|(kx&1),//for chroma-subsampled exJPEGs
					N[kc2+0]-W[kc2+0],
					2*FLOOR_LOG2(abs(N[kc2+1])*7)+FLOOR_LOG2(abs(W[kc2+1])*3),
					//N[kc2+0]+W[kc2+0]-NW[kc2+0],
					//N[kc2+0]-NN[kc2+0],
					//W[kc2+0]+WW[kc2+0],
					//NW[kc2+1],
					//N[kc2+1]+W[kc2+1]-NW[kc2+1],
					//abs(W[kc2+1]),
					((NEEE[kc2+0]<0)<<7|(NEE[kc2+0]<0)<<6|(NN[kc2+0]<0)<<5|(WW[kc2+0]<0)<<4|(NE[kc2+0]<0)<<3|(NW[kc2+0]<0)<<2|(W[kc2+0]<0)<<1|(N[kc2+0]<0))*2,
					((NEEE[kc2+1]<0)<<7|(NEE[kc2+1]<0)<<6|(NN[kc2+1]<0)<<5|(WW[kc2+1]<0)<<4|(NE[kc2+1]<0)<<3|(NW[kc2+1]<0)<<2|(W[kc2+1]<0)<<1|(N[kc2+1]<0))*2,
					(THREEWAY(N[kc2+1], 0)+1+(THREEWAY(W[kc2+1], 0)+1+(THREEWAY(NW[kc2+1], 0)+1+(THREEWAY(NE[kc2+1], 0)+1+(THREEWAY(N[kc2+1], W[kc2+1])+1)*3)*3)*3)*3),
					(kc?abs(curr[1]):N[kc2+2]*2),
					(kc>=2?abs(curr[3])*5+abs(curr[1]):W[kc2+2]/4),
					//(kc?abs(curr[0]):0)+(kc>=2?abs(curr[2]):0),
					//(kc?abs(curr[kc2-1]):abs(curr[kc2-2]))+(kc>=2?abs(curr[kc2-2]):abs(curr[kc2-3])),
					(abs(N[kc2+0])+abs(W[kc2+0]))/2-abs(W[kc2/2+1]),
					(abs(N[kc2+1])+abs(W[kc2+1]))*11-abs(NW[kc2+1])/2+abs(NE[kc2+1])*7,
				};
#ifdef __GNUC__
#pragma GCC unroll 12
#endif
				for(int k=0;k<A2_NCTX;++k)
				{
					int idx2=((A2_NCTX*kc+k)<<A2_CTXBITS|(ctx[k]>>(9-A2_CTXBITS)&((1<<A2_CTXBITS)-1)))<<8;
					curr_stats[k]=args->stats+idx2;
				}
				if(args->fwd)
				{
					curr[kc2+0]=yuv[kc];
					curr[kc2+1]=error=yuv[kc]-pred;
				}
				else
					error=0;
				int tidx=1;
				long long *curr_mixer=mixer+kc*8*A2_NCTX;
#ifdef ENABLE_SSE
				long long *curr_sse=args->sse[kc]
					[(N[kc2]-NW[kc2]+256)>>(9-A2_SSECBITS)&((1<<A2_SSECBITS)-1)]
					[(W[kc2]-NW[kc2]+256)>>(9-A2_SSECBITS)&((1<<A2_SSECBITS)-1)];
				//	[(W[kc2]-NW[kc2]+N[kc2]-NN[kc2]+512)>>(10-A2_SSECBITS)&((1<<A2_SSECBITS)-1)]
				//	[(N[kc2]-NW[kc2]+W[kc2]-WW[kc2]+512)>>(10-A2_SSECBITS)&((1<<A2_SSECBITS)-1)];
				//long long *curr_sse=args->sse+(kc*8LL<<A2_SSEPBITS);
#endif
				int shoffset=abs(N[kc2+1])+abs(W[kc2+1])*2+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1])/2-abs(WW[kc2+1])-abs(NN[kc2+1]);
#ifdef __GNUC__
#pragma GCC unroll 8
#endif
				for(int kb=7, e2=0;kb>=0;--kb)
				{
					long long p0=0;
					int wsum=0, bit;
#if 1
					ALIGN(32) long long probs[12];
					ALIGN(32) long long apsum[4], awsum[4];
					__m256i mprob0=_mm256_set_epi64x(curr_stats[ 3][tidx], curr_stats[ 2][tidx], curr_stats[ 1][tidx], curr_stats[ 0][tidx]);
					__m256i mprob1=_mm256_set_epi64x(curr_stats[ 7][tidx], curr_stats[ 6][tidx], curr_stats[ 5][tidx], curr_stats[ 4][tidx]);
					__m256i mprob2=_mm256_set_epi64x(curr_stats[11][tidx], curr_stats[10][tidx], curr_stats[ 9][tidx], curr_stats[ 8][tidx]);
					//ALIGN(32) long long probs[A2_NCTX];
					//for(int k=0;k<A2_NCTX;++k)
					//	probs[k]=curr_stats[k][tidx];
					//__m256i mprob0=_mm256_load_si256((__m256i*)probs+0);
					//__m256i mprob1=_mm256_load_si256((__m256i*)probs+1);
					//__m256i mprob2=_mm256_load_si256((__m256i*)probs+2);
					__m256i mweight0=_mm256_load_si256((__m256i*)curr_mixer+0);
					__m256i mweight1=_mm256_load_si256((__m256i*)curr_mixer+1);
					__m256i mweight2=_mm256_load_si256((__m256i*)curr_mixer+2);
					__m256i mprod0=_mm256_mul_epi32(mprob0, mweight0);
					__m256i mprod1=_mm256_mul_epi32(mprob1, mweight1);
					__m256i mprod2=_mm256_mul_epi32(mprob2, mweight2);
					__m256i mwsum=_mm256_add_epi64(mweight0, mweight1);
					mwsum=_mm256_add_epi64(mwsum, mweight2);
					__m256i mpsum=_mm256_add_epi64(mprod0, mprod1);
					mpsum=_mm256_add_epi64(mpsum, mprod2);
					_mm256_store_si256((__m256i*)awsum, mwsum);
					_mm256_store_si256((__m256i*)apsum, mpsum);
					wsum=(int)(awsum[0]+awsum[1]+awsum[2]+awsum[3]);
					p0=apsum[0]+apsum[1]+apsum[2]+apsum[3];
					wsum+=!wsum;
					p0/=wsum;

					//long long p0_2=0;
					//int wsum2=0;
					//for(int k=0;k<A2_NCTX;++k)
					//{
					//	p0_2+=(long long)curr_mixer[k]*curr_stats[k][tidx];
					//	wsum2+=curr_mixer[k];
					//}
					//if(wsum2)
					//	p0_2/=wsum2;
					//else
					//	p0_2=0x8000, wsum2=1;
					//if(p0!=p0_2)
					//	LOG_ERROR("");
#else
					ALIGN(32) long long probs[A2_NCTX];
					for(int k=0;k<A2_NCTX;++k)
					{
						probs[k]=curr_stats[k][tidx];
						p0+=(long long)curr_mixer[k]*probs[k];
						wsum+=curr_mixer[k];
					}
					if(wsum)
						p0/=wsum;
					else
						p0=0x8000, wsum=1;
#endif
#ifdef ENABLE_SSE
					int sseidx=(int)(p0>>(16-A2_SSEPBITS));
					CLAMP2(sseidx, 0, (1<<A2_SSEPBITS)-1);
					long long ssesum=curr_sse[sseidx]>>A2_SSECTR, ssecount=curr_sse[sseidx]&((1LL<<A2_SSECTR)-1);
					p0+=ssesum/(ssecount+5);
#endif
					CLAMP2(p0, 1, 0xFFFF);
					if(args->fwd)
					{
						bit=error>>kb&1;
						ac3_enc_bin(&ec, bit, (int)p0, 16);
					}
					else
					{
						bit=ac3_dec_bin(&ec, (int)p0, 16);
						error|=bit<<kb;
					}
					e2|=bit<<kb;
					int prob_error=(!bit<<16)-(int)p0;
#ifdef ENABLE_SSE
					++ssecount;
					ssesum+=prob_error;
					if(ssecount>0x4A00)
					{
						ssecount>>=1;
						ssesum>>=1;
					}
					curr_sse[sseidx]=ssesum<<A2_SSECTR|ssecount;
#endif
					int pred2=(char)(e2|1<<kb>>1)-pred;
					int sh=0;
					switch(kb)
					{
					case 7:
						sh=(shoffset+grad+(1<<7))>>7;
						break;
					case 6:
						sh=(shoffset+grad-abs(pred2)+(1<<7))>>7;
						break;
					case 5:
						sh=(shoffset+grad-abs(pred2)+(1<<7))>>7;
						break;
					case 4:
						sh=(shoffset+grad-abs(pred2)/2+(1<<7))>>7;
						break;
					case 3:
						sh=(shoffset+grad-abs(pred2)/4+(1<<7))>>7;
						break;
					case 2:
						sh=(shoffset+grad-abs(pred2)/8+(1<<7))>>7;
						break;
					case 1:
						sh=(shoffset+grad+abs(pred2)/4+(1<<7))>>7;
						break;
					case 0:
						sh=(shoffset-abs(W[kc2+1])+abs(NEE[kc2+1])/2+abs(NN[kc2+1])+7*abs(pred2))>>8;
						break;
					}
					if(sh<0)sh=0;
					sh=FLOOR_LOG2_P1(sh);
					if(abs(prob_error)>16)
					{
						long long dL_dp0=0x7F000000000/((long long)((bit<<16)-(int)p0)*wsum*2);
#if 1
						CLAMP2(dL_dp0, -0x3FFF, 0x3FFF);
						__m256i vmin2=_mm256_set1_epi64x(1);
						__m256i vmax2=_mm256_set1_epi64x(0x40000);
						__m256i lo32=_mm256_set1_epi64x(0xFFFFFFFF);
						__m256i msh=_mm256_set1_epi64x(sh+12LL);
						__m256i mp0=_mm256_set1_epi64x(p0);
						__m256i mdL_dp0=_mm256_set1_epi64x(dL_dp0);
						__m256i update0=_mm256_sub_epi32(mprob0, mp0);
						__m256i update1=_mm256_sub_epi32(mprob1, mp0);
						__m256i update2=_mm256_sub_epi32(mprob2, mp0);
						update0=_mm256_mul_epi32(update0, mdL_dp0);
						update1=_mm256_mul_epi32(update1, mdL_dp0);
						update2=_mm256_mul_epi32(update2, mdL_dp0);
						update0=_mm256_srav_epi32(update0, msh);
						update1=_mm256_srav_epi32(update1, msh);
						update2=_mm256_srav_epi32(update2, msh);
						mweight0=_mm256_sub_epi32(mweight0, update0);
						mweight1=_mm256_sub_epi32(mweight1, update1);
						mweight2=_mm256_sub_epi32(mweight2, update2);
						mweight0=_mm256_max_epi32(mweight0, vmin2);
						mweight1=_mm256_max_epi32(mweight1, vmin2);
						mweight2=_mm256_max_epi32(mweight2, vmin2);
						mweight0=_mm256_min_epi32(mweight0, vmax2);
						mweight1=_mm256_min_epi32(mweight1, vmax2);
						mweight2=_mm256_min_epi32(mweight2, vmax2);
						mweight0=_mm256_and_si256(mweight0, lo32);
						mweight1=_mm256_and_si256(mweight1, lo32);
						mweight2=_mm256_and_si256(mweight2, lo32);
						_mm256_store_si256((__m256i*)curr_mixer+0, mweight0);
						_mm256_store_si256((__m256i*)curr_mixer+1, mweight1);
						_mm256_store_si256((__m256i*)curr_mixer+2, mweight2);
#else
						for(int k=0;k<A2_NCTX;++k)
						{
							int update=(int)((int)dL_dp0*(curr_stats[k][tidx]-(int)p0)>>(sh+17-4-1));
							int mk=curr_mixer[k]-update;
							CLAMP2(mk, 1, 0x40000);
							curr_mixer[k]=mk;
						}
#endif
					}
#if 1
					__m256i msh=_mm256_set1_epi64x(sh);
					msh=_mm256_add_epi64(msh, _mm256_set_epi64x(5, 5, 4, 4));
					__m256i mtruth=_mm256_sllv_epi32(_mm256_set1_epi64x(1), msh);
					mtruth=_mm256_srli_epi32(mtruth, 1);
					mtruth=_mm256_or_si256(mtruth, _mm256_set1_epi64x((long long)!bit<<16));
					__m256i mupdate0=_mm256_sub_epi64(mtruth, mprob0);
					__m256i mupdate1=_mm256_sub_epi64(mtruth, mprob1);
					__m256i mupdate2=_mm256_sub_epi64(mtruth, mprob2);
					mupdate0=_mm256_srav_epi32(mupdate0, msh);
					mupdate1=_mm256_srav_epi32(mupdate1, msh);
					mupdate2=_mm256_srav_epi32(mupdate2, msh);
					mprob0=_mm256_add_epi64(mprob0, mupdate0);
					mprob1=_mm256_add_epi64(mprob1, mupdate1);
					mprob2=_mm256_add_epi64(mprob2, mupdate2);
					_mm256_store_si256((__m256i*)probs+0, mprob0);
					_mm256_store_si256((__m256i*)probs+1, mprob1);
					_mm256_store_si256((__m256i*)probs+2, mprob2);
#ifdef __GNUC__
#pragma GCC unroll 12
#endif
					for(int k=0;k<A2_NCTX;++k)
						curr_stats[k][tidx]=(int)probs[k];
#else
					for(int k=0;k<A2_NCTX;++k)
					{
						int p=curr_stats[k][tidx];
						p+=((!bit<<16)-p)>>(sh+4+(k>A2_NCTX/2));
						CLAMP2(p, 1, 0xFFFF);
						curr_stats[k][tidx]=p;
					//	CLAMP2_32(curr_stats[k][tidx], p, 1, 0xFFFF);
					}
#endif
					curr_mixer+=A2_NCTX;
#ifdef ENABLE_SSE
					curr_sse+=1<<A2_SSEPBITS;
#endif
					tidx+=tidx+bit;
				}
				if(!args->fwd)
				{
					error+=pred;
					error=error<<(32-8)>>(32-8);
					curr[kc2+0]=yuv[kc]=error;
					curr[kc2+1]=error-pred;
				}
#ifdef ESTIMATE_SIZE
				++args->hist2[kc][(curr[kc2+1]+128)&255];
#endif
				curr[kc2+0]-=offset;
#ifdef ENABLE_WG
				wg_update(curr[kc2+0], wg_preds, wg_perrors+WG_NPREDS*kc, eW, ecurr, eNE);
#endif
			}
			if(!args->fwd)
			{
				switch(bestrct)
				{
				case RCT_R_G_B:
				case RCT_R_G_BG:
				case RCT_R_G_BR:
				case RCT_R_GR_BR:
				case RCT_R_GR_BG:
				case RCT_R_G_B2:
				case RCT_R_GR_B2:
					args->dst[idx+0]=yuv[0]+128;
					args->dst[idx+1]=yuv[1]+128;
					args->dst[idx+2]=yuv[2]+128;
					break;
				case RCT_G_B_RG:
				case RCT_G_B_RB:
				case RCT_G_BG_RG:
				case RCT_G_BG_RB:
				case RCT_G_B_R2:
				case RCT_G_BG_R2:
					args->dst[idx+1]=yuv[0]+128;
					args->dst[idx+2]=yuv[1]+128;
					args->dst[idx+0]=yuv[2]+128;
					break;
				case RCT_B_R_GR:
				case RCT_B_R_GB:
				case RCT_B_RB_GB:
				case RCT_B_RB_GR:
				case RCT_B_RB_G2:
					args->dst[idx+2]=yuv[0]+128;
					args->dst[idx+0]=yuv[1]+128;
					args->dst[idx+1]=yuv[2]+128;
					break;
				case RCT_G_RG_BR:
				case RCT_G_RG_B2:
					args->dst[idx+1]=yuv[0]+128;
					args->dst[idx+0]=yuv[1]+128;
					args->dst[idx+2]=yuv[2]+128;
					break;
				case RCT_B_GB_RG:
				case RCT_B_GB_R2:
					args->dst[idx+2]=yuv[0]+128;
					args->dst[idx+1]=yuv[1]+128;
					args->dst[idx+0]=yuv[2]+128;
					break;
				case RCT_R_BR_GB:
				case RCT_R_B_G2:
				case RCT_R_BR_G2:
					args->dst[idx+0]=yuv[0]+128;
					args->dst[idx+2]=yuv[1]+128;
					args->dst[idx+1]=yuv[2]+128;
					break;
				}
#ifdef ENABLE_GUIDE
				if(args->test&&memcmp(args->dst+idx, args->src+idx, sizeof(char)*nch))
				{
					unsigned char orig[4]={0};
					memcpy(orig, args->src+idx, nch*sizeof(char));
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
			}
			__m256i mr=_mm256_load_si256((__m256i*)rows);
			__m256i er=_mm256_load_si256((__m256i*)erows);
			mr=_mm256_add_epi64(mr, _mm256_set1_epi64x(sizeof(short[4*2])));
			er=_mm256_add_epi64(er, _mm256_set1_epi64x(sizeof(int[4*WG_NPREDS])));
			_mm256_store_si256((__m256i*)rows, mr);
			_mm256_store_si256((__m256i*)erows, er);
			//rows[0]+=4*2;
			//rows[1]+=4*2;
			//rows[2]+=4*2;
			//rows[3]+=4*2;
		}
	}
	if(args->fwd)
		ac3_enc_flush(&ec);
}
static void block_manager(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	for(int kb=args->b1;kb<args->b2;++kb)
	{
		int kx, ky;

		args->blockidx=kb;
		args->currentblock=kb-args->b1;
		ky=args->blockidx/args->xblocks;
		kx=args->blockidx%args->xblocks;
		args->x1=BLOCKX*kx;
		args->y1=BLOCKY*ky;
		args->x2=MINVAR(args->x1+BLOCKX, args->iw);
		args->y2=MINVAR(args->y1+BLOCKY, args->ih);
		if(!args->fwd)
		{
			args->decstart=args->decsrc+args->offsets[kb];
			args->decend=args->decsrc+args->offsets[kb+1];
		}
		block_thread(param);
	}
}
int c03_codec(const char *srcfn, const char *dstfn)
{
	const int nch=3;
//	const int depth=8;
	double t0;
	ArrayHandle src, dst;
	int headersize, printed;
	int iw, ih;
	const unsigned char *image, *imageend;
	unsigned char *image2;
	CodecID codec;
	int xblocks, yblocks, nblocks, ncores, nthreads, blocksperthread, coffset;
	ptrdiff_t start, memusage, argssize;
	ThreadArgs *args;
	int test, fwd;
	int histsize;
	double bestsize;
	int usize;
#ifdef ESTIMATE_SIZE
	double esize[3]={0};
#endif
	
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
	
	if(test)
		printf("%s \"%s\"  WH %d*%d\n", CODECNAME, srcfn, iw, ih);
	usize=iw*ih*3;
	xblocks=(iw+BLOCKX-1)/BLOCKX;
	yblocks=(ih+BLOCKY-1)/BLOCKY;
	nblocks=xblocks*yblocks;
	ncores=query_cpu_cores();
	nthreads=MINVAR(nblocks, ncores);
	blocksperthread=(nblocks+nthreads-1)/nthreads;
	coffset=(int)sizeof(int)*nblocks;
	int *offsets=(int*)malloc(coffset+sizeof(int));
	if(!offsets)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	start=0;
	memusage=0;
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
			offsets[kt]=(int)start;
			start+=size;
		}
		offsets[nblocks]=(int)start;
		if(image+start!=imageend)
			LOG_ERROR("Corrupt file");

		image2=(unsigned char*)malloc(usize);
		if(!image2)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(image2, 0, usize);
	}
	argssize=nthreads*sizeof(ThreadArgs);
	args=(ThreadArgs*)malloc(argssize);
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	bestsize=0;
	memusage+=argssize;
	memset(args, 0, argssize);
	histsize=(int)sizeof(int[OCH_COUNT*PRED_COUNT<<8]);
	for(int k=0, bidx=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=image;
		arg->dst=fwd?0:dst->data+printed;
		arg->iw=iw;
		arg->ih=ih;
		
		arg->xblocks=xblocks;
		arg->b1=bidx;
		bidx+=blocksperthread;
		if(bidx>nblocks)
			bidx=nblocks;
		arg->b2=bidx;
		arg->offsets=offsets;
		if(!fwd)
			arg->decsrc=image;

		int listssize=((size_t)arg->b2-arg->b1)*sizeof(BList);
		arg->lists=(BList*)malloc(listssize);

		arg->bufsize=sizeof(short[4*OCH_COUNT*2*(BLOCKX+16LL)]);//4 padded rows * OCH_COUNT * {pixels, wg_errors}
		arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));

		arg->ebufsize=(int)sizeof(int[4][4][WG_NPREDS][BLOCKX+16LL]);
		arg->ebuf=(int*)_mm_malloc(arg->ebufsize, sizeof(__m256i));

		arg->histsize=histsize;
		arg->hist=(int*)malloc(histsize);
		if(!arg->lists||!arg->pixels||!arg->ebuf||!arg->hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=listssize;
		memusage+=arg->ebufsize;
		memusage+=arg->bufsize;
		memusage+=histsize;

		arg->fwd=fwd;
		arg->test=test;
#ifdef DISABLE_MT
		arg->loud=0;
	//	arg->loud=test&&nblocks<MAXPRINTEDBLOCKS;
#else
		arg->loud=0;
#endif
	}
	for(int k2=0;k2<=test;++k2)
	{
#ifdef DISABLE_MT
		for(int k=0;k<nthreads;++k)
			block_manager(args+k);
#else
		void *ctx=mt_exec(block_manager, args, sizeof(ThreadArgs), nthreads);
		mt_finish(ctx);
#endif
		if(fwd)
		{
			for(int kt1=0;kt1<nthreads;++kt1)
			{
				ThreadArgs *arg=args+kt1;
				for(int kt2=0;kt2<blocksperthread;++kt2)
				{
					int bidx=kt1*blocksperthread+kt2;
					if(bidx>=nblocks)
						break;
					memcpy(dst->data+printed+sizeof(int)*bidx, &arg->lists[kt2].nbytes, sizeof(int));
					blist_appendtoarray(arg->lists+kt2, &dst);
					blist_clear(arg->lists+kt2);
				}
				if(fwd)
					bestsize+=arg->bestsize;
			}
		}
		if(test)
		{
			ptrdiff_t csize=dst->count;
			t0=time_sec()-t0;
			if(fwd)
			{
#ifdef ESTIMATE_SIZE
				esize[0]/=8;
				esize[1]/=8;
				esize[2]/=8;
				printf("T %15.2lf\n", esize[0]+esize[1]+esize[2]);
				printf("Y %15.2lf\n", esize[0]);
				printf("U %15.2lf\n", esize[1]);
				printf("V %15.2lf\n", esize[2]);
#endif
				printf("%15.2lf (%+13.2lf) bytes\n", bestsize, csize-bestsize);
				printf("%12td/%12d  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
				printf("Mem usage: ");
				print_size((double)memusage, 8, 4, 0, 0);
				printf("\n");
			}
			printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
			if(!fwd)
				compare_bufs_8(image2, src->data+headersize, iw, ih, nch, nch, CODECNAME, 0, 1);
		}
		if(fwd&&test)//transition to (test) decode
		{
			fwd=0;
			image2=(unsigned char*)malloc(usize);
			if(!image2)
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			memset(image2, 0, usize);
			start=coffset;
			for(int kt=0;kt<nblocks;++kt)
			{
				int size=0;
				memcpy(&size, dst->data+printed+sizeof(int)*kt, sizeof(int));
				offsets[kt]=(int)start;
				start+=size;
			}
			offsets[nblocks]=(int)start;
			for(int kt=0;kt<nthreads;++kt)
			{
				ThreadArgs *arg=args+kt;
				arg->dst=image2;
				arg->fwd=0;
				arg->decsrc=dst->data+printed;
			}
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
		free(arg->lists);
		free(arg->hist);
		_mm_free(arg->pixels);
		_mm_free(arg->ebuf);
	}
	free(args);
	free(offsets);
	array_free(&src);
	array_free(&dst);
	//if(test)
	//	pause();
	return 0;
}