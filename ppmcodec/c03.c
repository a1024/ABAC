#include"codec.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT

	#define ENABLE_SSE//good
//	#define PRINT_SH
	#define ESTIMATE_SIZE

	#define ENABLE_OLS
//	#define ENABLE_WG

#define CODECNAME "C03"
#define AC3_PREC
#include"entropy.h"

#define BLOCKSIZE 1024
#define MAXPRINTEDBLOCKS 20

#define OLS_STRIDE 255
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

#define PREDLIST\
	PRED(W)\
	PRED(CG)\
	PRED(AV5)\
	PRED(AV9)\
	PRED(AV12)
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
	NB_NNWW,	NB_NNW,		NB_NN,		NB_NNE,		NB_NNEE,
	NB_NWW,		NB_NW,		NB_N,		NB_NE,		NB_NEE,
	NB_WW,		NB_W,		NB_curr,
} NBIndex;
static const short av12_icoeffs[12]=
{
	 0x04,	 0x03,	-0x1F,	-0x26,	 0x00,
	 0x07,	-0x9E,	 0xDB,	 0x1E,	 0x13,
	-0x2A,	 0xF3,
};


//OLS4
typedef struct _OLS4Context
{
	int nparams;
	double *nb, *vec, *cov, *cholesky, *params;
} OLS4Context;
static void ols4_init(OLS4Context *ctx, int nparams)
{
	int msize=nparams*nparams;
	memset(ctx, 0, sizeof(*ctx));
	ctx->nparams=nparams;
	ctx->nb=(double*)_mm_malloc(nparams*sizeof(double), sizeof(__m256d));
	ctx->vec=(double*)_mm_malloc(nparams*sizeof(double), sizeof(__m256d));
	ctx->cov=(double*)_mm_malloc(msize*sizeof(double), sizeof(__m256d));
	ctx->cholesky=(double*)malloc(msize*sizeof(double));
	ctx->params=(double*)_mm_malloc(nparams*sizeof(double), sizeof(__m256d));
	if(!ctx->nb||!ctx->vec||!ctx->cov||!ctx->cholesky||!ctx->params)
	{
		LOG_ERROR("Alloc error");
		return;
	}
}
static void ols4_reset(OLS4Context *ctx)
{
	int msize=ctx->nparams*ctx->nparams;
	memset(ctx->nb, 0, ctx->nparams*sizeof(double));
	memset(ctx->vec, 0, ctx->nparams*sizeof(double));
	memset(ctx->cov, 0, msize*sizeof(double));
	memset(ctx->cholesky, 0, msize*sizeof(double));
	memset(ctx->params, 0, ctx->nparams*sizeof(double));
}
static void ols4_free(OLS4Context *ctx)
{
	_mm_free(ctx->nb);
	_mm_free(ctx->vec);
	_mm_free(ctx->cov);
	free(ctx->cholesky);
	_mm_free(ctx->params);
	memset(ctx, 0, sizeof(*ctx));
}
#define OLS4_NPARAMS 16		//multiple of 4,  NPARAMS is unified because of interleaving


//WG:
#define WG_DECAY_NUM	493
#define WG_DECAY_SH	9

#define WG_NPREDS	16
#define WG_PREDLIST\
	WG_PRED(25, spred)\
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
static void wg_init(int *weights)
{
	int j=0;
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
	WG_PREDLIST
#undef  WG_PRED
}
static int wg_predict(
	const int *weights,
	short **rows, const int stride, int kc2,
	int spred, const int *perrors, int *preds
)
{
	double pred2=0, wsum=0;
	int j=0, pred;
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

#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST
#undef  WG_PRED
	
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
static void wg_update(int curr, const int *preds, int *perrors, int *weights)
{
	for(int k=0;k<WG_NPREDS;++k)
	{
		int e2=abs(curr-preds[k]);
		perrors[k]=(perrors[k]+e2)*WG_DECAY_NUM>>WG_DECAY_SH;
	}
}


//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
#define CONFIG_EXP 4
#define CONFIG_MSB 1
#define CONFIG_LSB 0
static void quantize_pixel(int val, int *token, int *bypass, int *nbits)
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
#define ABAC_PROBBITS 16
static int squash(int x)//sigmoid(x) = 1/(1-exp(-x))		logit sum -> prob
{
#ifdef DISABLE_LOGMIX
	x>>=11;
	x+=1<<ABAC_PROBBITS>>1;
	CLAMP2_32(x, x, 1, (1<<ABAC_PROBBITS)-1);
#else
	static const int t[33]=//2^5 table elements, table amplitude 2^12
	{
		   1,    2,    3,    6,   10,   16,   27,   45,   73,  120,  194,
		 310,  488,  747, 1101, 1546, 2047, 2549, 2994, 3348, 3607, 3785,
		3901, 3975, 4022, 4050, 4068, 4079, 4085, 4089, 4092, 4093, 4094,
	};
	int w=x&((1<<(ABAC_PROBBITS-5))-1);
	x=(x>>(ABAC_PROBBITS-5))+16;
	if(x>31)
		return (1<<ABAC_PROBBITS)-1;
	if(x<0)
		return 1;
	x=(t[x]*((1<<(ABAC_PROBBITS-5))-w)+t[x+1]*w+64)>>(12-5);
#endif
	return x;
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
	
	//int hist[54<<8];

	int tlevels;
	unsigned short stats[A2_NCTX*3<<A2_CTXBITS<<8];
#ifdef ENABLE_SSE
	long long sse[3][1<<A2_SSECBITS][1<<A2_SSECBITS][8<<A2_SSEPBITS];
	//long long sse[24<<A2_SSEPBITS];
#endif

	//aux
	int blockidx;
	double bestsize;
	int bestrct, predidx[3];
#ifdef ENABLE_OLS
	OLS4Context ols4[3];
#endif
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
#ifdef ENABLE_OLS
	ALIGN(16) int ols4_ctx0[OLS4_NPARAMS]={0}, ols4_ctx1[OLS4_NPARAMS]={0}, ols4_ctx2[OLS4_NPARAMS]={0};
	static const double lrs[]={0.0018, 0.003, 0.003};
#endif
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
					nb0[13],//rgb
					nb1[13],//gbr
					nb2[13],//rgb - gbr
					nb3[13],//gbr - rgb
					nb4[13],//(gbr+brg)/2
					nb5[13];//rgb - (gbr+brg)/2
				__m256i vmin[4], vmax[4], pred;
				{
					__m128i nb8[13]=//8-bit
					{
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride-2*3+0)), half8),//NNWW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride-1*3+0)), half8),//NNW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride+0*3+0)), half8),//NN
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride+1*3+0)), half8),//NNE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride+2*3+0)), half8),//NNEE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride-2*3+0)), half8),//NWW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride-1*3+0)), half8),//NW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+0*3+0)), half8),//N
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+1*3+0)), half8),//NE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+2*3+0)), half8),//NEE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride-2*3+0)), half8),//WW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride-1*3+0)), half8),//W
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride+0*3+0)), half8),//curr
					};
					for(int k=0;k<13;++k)
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
		args->bestsize=bestsize;
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
#ifdef ENABLE_OLS
	int kols=0, kols2=OLS_STRIDE;
	for(int k=0;k<3;++k)
		ols4_reset(args->ols4+k);
#endif
#ifdef ENABLE_WG
	int wg_weights[WG_NPREDS*3]={0}, wg_perrors[WG_NPREDS*3]={0}, wg_preds[WG_NPREDS]={0};
	wg_init(wg_weights+WG_NPREDS*0);
	wg_init(wg_weights+WG_NPREDS*1);
	wg_init(wg_weights+WG_NPREDS*2);
#endif
	memset(args->pixels, 0, args->bufsize);
#ifdef ABAC_SIMPLEOVERRIDE
	unsigned short stats2[768]={0};
	FILLMEM(stats2, 0x8000, sizeof(short[768]), sizeof(short));
#endif
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) short *rows[]=
		{
			args->pixels+((BLOCKSIZE+16LL)*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKSIZE+16LL)*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKSIZE+16LL)*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKSIZE+16LL)*((ky-3LL)&3)+8LL)*4*2,
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
#ifdef ENABLE_OLS
			ALIGN(16) int preds[4];
			__m128i mNW	=_mm_load_si128((__m128i*)NW);
			__m128i mN	=_mm_load_si128((__m128i*)N);
			__m128i mNE	=_mm_load_si128((__m128i*)NE);
			__m128i mW	=_mm_load_si128((__m128i*)W);
			__m128i vmin=_mm_min_epi16(mN, mW);
			__m128i vmax=_mm_max_epi16(mN, mW);
			__m128i mg=_mm_sub_epi16(_mm_add_epi16(mN, mW), mNW);
			mg=_mm_max_epi16(mg, vmin);
			mg=_mm_min_epi16(mg, vmax);
			ALIGN(16) short cgrads[8];
			_mm_store_si128((__m128i*)cgrads, mg);
			vmin=_mm_min_epi16(vmin, mNE);
			vmax=_mm_max_epi16(vmax, mNE);
			vmin=_mm_slli_epi32(vmin, 16);//need low signed 16 bits -> 32 bit
			vmax=_mm_slli_epi32(vmax, 16);
			vmin=_mm_srai_epi32(vmin, 16);
			vmax=_mm_srai_epi32(vmax, 16);
			int j0=0, j1=0, j2=0;
			ols4_ctx0[j0++]=cgrads	[0];
			ols4_ctx0[j0++]=NNW	[0];
			ols4_ctx0[j0++]=NN	[0];
			ols4_ctx0[j0++]=NNE	[0];
			ols4_ctx0[j0++]=NNEE	[0];
			ols4_ctx0[j0++]=NNEEE	[0];
			ols4_ctx0[j0++]=NWWW	[0];
			ols4_ctx0[j0++]=NWW	[0];
			ols4_ctx0[j0++]=NW	[0];
			ols4_ctx0[j0++]=N	[0];
			ols4_ctx0[j0++]=NE	[0];
			ols4_ctx0[j0++]=NEE	[0];
			ols4_ctx0[j0++]=NEEE	[0];
			ols4_ctx0[j0++]=WWW	[0];
			ols4_ctx0[j0++]=WW	[0];
			ols4_ctx0[j0++]=W	[0];

			ols4_ctx1[j1++]=cgrads	[2];
			ols4_ctx1[j1++]=NNW	[2];
			ols4_ctx1[j1++]=NN	[2];
			ols4_ctx1[j1++]=NNE	[2];
			ols4_ctx1[j1++]=NNEE	[2];
			ols4_ctx1[j1++]=NWW	[2];
			ols4_ctx1[j1++]=NW	[2];
			ols4_ctx1[j1++]=NW	[0];
			ols4_ctx1[j1++]=N	[2];
			ols4_ctx1[j1++]=N	[0];
			ols4_ctx1[j1++]=NE	[2];
			ols4_ctx1[j1++]=NE	[0];
			ols4_ctx1[j1++]=NEE	[2];
			ols4_ctx1[j1++]=WW	[2];
			ols4_ctx1[j1++]=W	[2];
			ols4_ctx1[j1++]=W	[0];

			ols4_ctx2[j2++]=cgrads	[4];
			ols4_ctx2[j2++]=NNW	[4];
			ols4_ctx2[j2++]=NN	[4];
			ols4_ctx2[j2++]=NNE	[4];
			ols4_ctx2[j2++]=NNEE	[4];
			ols4_ctx2[j2++]=NWW	[4];
			ols4_ctx2[j2++]=NW	[4];
			ols4_ctx2[j2++]=NW	[0];
			ols4_ctx2[j2++]=N	[4];
			ols4_ctx2[j2++]=N	[0];
			ols4_ctx2[j2++]=NE	[4];
			ols4_ctx2[j2++]=NE	[0];
			ols4_ctx2[j2++]=NEE	[4];
			ols4_ctx2[j2++]=WW	[4];
			ols4_ctx2[j2++]=W	[4];
			ols4_ctx2[j2++]=W	[0];
			{
				double *pp0=args->ols4[0].params;
				double *pp1=args->ols4[1].params;
				double *pp2=args->ols4[2].params;
				double *pnb0=args->ols4[0].nb;
				double *pnb1=args->ols4[1].nb;
				double *pnb2=args->ols4[2].nb;

				__m128i mnb0=_mm_load_si128((__m128i*)ols4_ctx0+0);
				__m128i mnb1=_mm_load_si128((__m128i*)ols4_ctx1+0);
				__m128i mnb2=_mm_load_si128((__m128i*)ols4_ctx2+0);
				__m256d dnb0=_mm256_cvtepi32_pd(mnb0);
				__m256d dnb1=_mm256_cvtepi32_pd(mnb1);
				__m256d dnb2=_mm256_cvtepi32_pd(mnb2);
				_mm256_store_pd(pnb0+0*4, dnb0);
				_mm256_store_pd(pnb1+0*4, dnb1);
				_mm256_store_pd(pnb2+0*4, dnb2);
				__m256d sum0=_mm256_mul_pd(dnb0, _mm256_load_pd(pp0+0*4));
				__m256d sum1=_mm256_mul_pd(dnb1, _mm256_load_pd(pp1+0*4));
				__m256d sum2=_mm256_mul_pd(dnb2, _mm256_load_pd(pp2+0*4));

				mnb0=_mm_load_si128((__m128i*)ols4_ctx0+1);
				mnb1=_mm_load_si128((__m128i*)ols4_ctx1+1);
				mnb2=_mm_load_si128((__m128i*)ols4_ctx2+1);
				dnb0=_mm256_cvtepi32_pd(mnb0);
				dnb1=_mm256_cvtepi32_pd(mnb1);
				dnb2=_mm256_cvtepi32_pd(mnb2);
				_mm256_store_pd(pnb0+1*4, dnb0);
				_mm256_store_pd(pnb1+1*4, dnb1);
				_mm256_store_pd(pnb2+1*4, dnb2);
				sum0=_mm256_add_pd(sum0, _mm256_mul_pd(dnb0, _mm256_load_pd(pp0+1*4)));
				sum1=_mm256_add_pd(sum1, _mm256_mul_pd(dnb1, _mm256_load_pd(pp1+1*4)));
				sum2=_mm256_add_pd(sum2, _mm256_mul_pd(dnb2, _mm256_load_pd(pp2+1*4)));

				mnb0=_mm_load_si128((__m128i*)ols4_ctx0+2);
				mnb1=_mm_load_si128((__m128i*)ols4_ctx1+2);
				mnb2=_mm_load_si128((__m128i*)ols4_ctx2+2);
				dnb0=_mm256_cvtepi32_pd(mnb0);
				dnb1=_mm256_cvtepi32_pd(mnb1);
				dnb2=_mm256_cvtepi32_pd(mnb2);
				_mm256_store_pd(pnb0+2*4, dnb0);
				_mm256_store_pd(pnb1+2*4, dnb1);
				_mm256_store_pd(pnb2+2*4, dnb2);
				sum0=_mm256_add_pd(sum0, _mm256_mul_pd(dnb0, _mm256_load_pd(pp0+2*4)));
				sum1=_mm256_add_pd(sum1, _mm256_mul_pd(dnb1, _mm256_load_pd(pp1+2*4)));
				sum2=_mm256_add_pd(sum2, _mm256_mul_pd(dnb2, _mm256_load_pd(pp2+2*4)));

				mnb0=_mm_load_si128((__m128i*)ols4_ctx0+3);
				mnb1=_mm_load_si128((__m128i*)ols4_ctx1+3);
				mnb2=_mm_load_si128((__m128i*)ols4_ctx2+3);
				dnb0=_mm256_cvtepi32_pd(mnb0);
				dnb1=_mm256_cvtepi32_pd(mnb1);
				dnb2=_mm256_cvtepi32_pd(mnb2);
				_mm256_store_pd(pnb0+3*4, dnb0);
				_mm256_store_pd(pnb1+3*4, dnb1);
				_mm256_store_pd(pnb2+3*4, dnb2);
				sum0=_mm256_add_pd(sum0, _mm256_mul_pd(dnb0, _mm256_load_pd(pp0+3*4)));
				sum1=_mm256_add_pd(sum1, _mm256_mul_pd(dnb1, _mm256_load_pd(pp1+3*4)));
				sum2=_mm256_add_pd(sum2, _mm256_mul_pd(dnb2, _mm256_load_pd(pp2+3*4)));

				ALIGN(32) double sums[12];
				_mm256_store_pd(sums+0*4, sum0);
				_mm256_store_pd(sums+1*4, sum1);
				_mm256_store_pd(sums+2*4, sum2);

				sums[ 0]+=sums[ 1]+sums[ 2]+sums[ 3];
				sums[ 4]+=sums[ 5]+sums[ 6]+sums[ 7];
				sums[ 8]+=sums[ 9]+sums[10]+sums[11];
#if OLS4_NPARAMS&15
				for(int k=16;k<OLS4_NPARAMS;++k)
				{
					pnb0[k]=(double)ols4_ctx0[k];
					pnb1[k]=(double)ols4_ctx1[k];
					pnb2[k]=(double)ols4_ctx2[k];
					sums[ 0]+=pp0[k]*pnb0[k];
					sums[ 4]+=pp1[k]*pnb1[k];
					sums[ 8]+=pp2[k]*pnb2[k];
				}
#endif
				sums[1]=sums[4];
				sums[2]=sums[8];

				__m128i mp=_mm256_cvtpd_epi32(_mm256_load_pd(sums));//rounded
				mp=_mm_min_epi32(mp, vmax);
				mp=_mm_max_epi32(mp, vmin);
				_mm_store_si128((__m128i*)preds, mp);
			}
#endif
			for(int kc=0;kc<nch;++kc)
			{
				int kc2=kc<<1;
				int offset=(yuv[combination[kc+3]]+yuv[combination[kc+6]])>>combination[kc+9];
#ifdef ENABLE_OLS
				pred=preds[kc];
#else
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
					//p0-=0x8000;
					//p0>>=1;
					//p0=squash((int)p0);
					CLAMP2(p0, 1, 0xFFFF);
					//CLAMP2_32(p0, (int)p0, 1, 0xFFFF);
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
#ifdef PRINT_SH
					if(args->fwd&&ky<<1>=args->ih-5&&ky<<1<=args->ih+5&&kx<<1>=args->iw-5&&kx<<1<=args->iw+5)
					{
						if(kb==7)
							printf("%5d %5d", kx, ky);
						printf(" %d", sh);
						if(kb==0)
						{
							printf(" 0x%02X = %d\n", error&0xFF, error);
							if(kc==2)
								printf("\n");
						}
					//	printf("%5d %5d %d %d\n", kx, ky, kb, sh);
					}
#endif
					if(abs(prob_error)>16)
					{
#if 1
						long long dL_dp0=0x7F000000000/((long long)((bit<<16)-(int)p0)*wsum*2);
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
						//for(int k=0;k<A2_NCTX;++k)
						//{
						//	int update=(int)((int)dL_dp0*(curr_stats[k][tidx]-(int)p0)>>(sh+17-4-1));
						//	int mk=curr_mixer[k]-update;
						//	CLAMP2(mk, 1, 0x40000);
						//	curr_mixer[k]=mk;
						//}
#else
						int dL_dp0=0x7F00000/((bit<<16)-(int)p0);
						for(int k=0;k<A2_NCTX;++k)
						{
							int mk=curr_mixer[k]-(int)((long long)dL_dp0*(curr_stats[k][tidx]-(int)p0)>>sh)/wsum;
							CLAMP2(mk, 1, 0x40000);
							curr_mixer[k]=mk;
						//	CLAMP2_32(curr_mixer[k], mk, 1, 0x40000);
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
				wg_update(curr[kc2+0], wg_preds, wg_perrors+WG_NPREDS*kc, wg_weights+WG_NPREDS*kc);
#endif
			}
#ifdef ENABLE_OLS
			if(kx<5||kols>=kols2)
			{
				double *pnb0=args->ols4[0].nb;
				double *pnb1=args->ols4[1].nb;
				double *pnb2=args->ols4[2].nb;
				double *cov0=args->ols4[0].cov;
				double *cov1=args->ols4[1].cov;
				double *cov2=args->ols4[2].cov;
				double lval0=curr[0]*lrs[0], lr_comp0=1-lrs[0];
				double lval1=curr[2]*lrs[1], lr_comp1=1-lrs[1];
				double lval2=curr[4]*lrs[2], lr_comp2=1-lrs[2];
				__m256d mlr0=_mm256_set1_pd(lrs[0]);
				__m256d mlr1=_mm256_set1_pd(lrs[1]);
				__m256d mlr2=_mm256_set1_pd(lrs[2]);
				for(int ky2=0, midx=0;ky2<OLS4_NPARAMS;++ky2)
				{
					__m256d vy0=_mm256_set1_pd(pnb0[ky2]);
					__m256d vy1=_mm256_set1_pd(pnb1[ky2]);
					__m256d vy2=_mm256_set1_pd(pnb2[ky2]);
					int kx2=0;
					for(;kx2<OLS4_NPARAMS-3;kx2+=4, midx+=4)
					{
						__m256d mcov0=_mm256_load_pd(cov0+midx);
						__m256d mcov1=_mm256_load_pd(cov1+midx);
						__m256d mcov2=_mm256_load_pd(cov2+midx);
						__m256d vx0=_mm256_load_pd(pnb0+kx2);
						__m256d vx1=_mm256_load_pd(pnb1+kx2);
						__m256d vx2=_mm256_load_pd(pnb2+kx2);
						vx0=_mm256_mul_pd(vx0, vy0);
						vx1=_mm256_mul_pd(vx1, vy1);
						vx2=_mm256_mul_pd(vx2, vy2);
						vx0=_mm256_sub_pd(vx0, mcov0);
						vx1=_mm256_sub_pd(vx1, mcov1);
						vx2=_mm256_sub_pd(vx2, mcov2);
						vx0=_mm256_mul_pd(vx0, mlr0);
						vx1=_mm256_mul_pd(vx1, mlr1);
						vx2=_mm256_mul_pd(vx2, mlr2);
						vx0=_mm256_add_pd(vx0, mcov0);
						vx1=_mm256_add_pd(vx1, mcov1);
						vx2=_mm256_add_pd(vx2, mcov2);
						_mm256_store_pd(cov0+midx, vx0);
						_mm256_store_pd(cov1+midx, vx1);
						_mm256_store_pd(cov2+midx, vx2);
					}
#if OLS4_NPARAMS&15
					for(;kx2<OLS4_NPARAMS;++kx2, ++midx)
					{
						cov0[midx]+=(pnb0[kx2]*pnb0[ky2]-cov0[midx])*lrs[0];
						cov1[midx]+=(pnb1[kx2]*pnb1[ky2]-cov1[midx])*lrs[1];
						cov2[midx]+=(pnb2[kx2]*pnb2[ky2]-cov2[midx])*lrs[2];
					}
#endif
				}
				double *vec0=args->ols4[0].vec;
				double *vec1=args->ols4[1].vec;
				double *vec2=args->ols4[2].vec;
				{
					__m256d mlr_comp0=_mm256_set1_pd(lr_comp0);
					__m256d mlr_comp1=_mm256_set1_pd(lr_comp1);
					__m256d mlr_comp2=_mm256_set1_pd(lr_comp2);
					__m256d mval0=_mm256_set1_pd(lval0);
					__m256d mval1=_mm256_set1_pd(lval1);
					__m256d mval2=_mm256_set1_pd(lval2);
					int k=0;
					for(;k<OLS4_NPARAMS-3;k+=4)
					{
						__m256d mvec0=_mm256_load_pd(vec0+k);
						__m256d mvec1=_mm256_load_pd(vec1+k);
						__m256d mvec2=_mm256_load_pd(vec2+k);
						__m256d mtmp0=_mm256_load_pd(pnb0+k);
						__m256d mtmp1=_mm256_load_pd(pnb1+k);
						__m256d mtmp2=_mm256_load_pd(pnb2+k);
						mvec0=_mm256_mul_pd(mvec0, mlr_comp0);
						mvec1=_mm256_mul_pd(mvec1, mlr_comp1);
						mvec2=_mm256_mul_pd(mvec2, mlr_comp2);
						mtmp0=_mm256_mul_pd(mtmp0, mval0);
						mtmp1=_mm256_mul_pd(mtmp1, mval1);
						mtmp2=_mm256_mul_pd(mtmp2, mval2);
						mvec0=_mm256_add_pd(mvec0, mtmp0);
						mvec1=_mm256_add_pd(mvec1, mtmp1);
						mvec2=_mm256_add_pd(mvec2, mtmp2);
						_mm256_store_pd(vec0+k, mvec0);
						_mm256_store_pd(vec1+k, mvec1);
						_mm256_store_pd(vec2+k, mvec2);
					}
#if OLS4_NPARAMS&15
					for(;k<OLS4_NPARAMS;++k)
					{
						vec0[k]=lval0*pnb0[k]+lr_comp0*vec0[k];
						vec1[k]=lval1*pnb1[k]+lr_comp1*vec1[k];
						vec2[k]=lval2*pnb2[k]+lr_comp2*vec2[k];
					}
#endif
				}
				double *chol0=args->ols4[0].cholesky;
				double *chol1=args->ols4[1].cholesky;
				double *chol2=args->ols4[2].cholesky;
				double *params0=args->ols4[0].params;
				double *params1=args->ols4[1].params;
				double *params2=args->ols4[2].params;
				int success0=1;
				int success1=1;
				int success2=1;
				memcpy(chol0, cov0, sizeof(double[OLS4_NPARAMS*OLS4_NPARAMS]));
				memcpy(chol1, cov1, sizeof(double[OLS4_NPARAMS*OLS4_NPARAMS]));
				memcpy(chol2, cov2, sizeof(double[OLS4_NPARAMS*OLS4_NPARAMS]));
				for(int k=0;k<OLS4_NPARAMS*OLS4_NPARAMS;k+=OLS4_NPARAMS+1)
				{
					chol0[k]+=0.0075;
					chol1[k]+=0.0075;
					chol2[k]+=0.0075;
				}
				double sum;
				for(int i=0;i<OLS4_NPARAMS;++i)
				{
					for(int j=0;j<i;++j)
					{
						sum=chol0[i*OLS4_NPARAMS+j];
						for(int k=0;k<j;++k)
							sum-=chol0[i*OLS4_NPARAMS+k]*chol0[j*OLS4_NPARAMS+k];
						chol0[i*OLS4_NPARAMS+j]=sum/chol0[j*OLS4_NPARAMS+j];
					}
					sum=chol0[i*OLS4_NPARAMS+i];
					for(int k=0;k<i;++k)
						sum-=chol0[i*OLS4_NPARAMS+k]*chol0[i*OLS4_NPARAMS+k];
					if(sum<=1e-8)
					{
						success0=0;
						break;
					}
					chol0[i*OLS4_NPARAMS+i]=sqrt(sum);
				}
				for(int i=0;i<OLS4_NPARAMS;++i)
				{
					for(int j=0;j<i;++j)
					{
						sum=chol1[i*OLS4_NPARAMS+j];
						for(int k=0;k<j;++k)
							sum-=chol1[i*OLS4_NPARAMS+k]*chol1[j*OLS4_NPARAMS+k];
						chol1[i*OLS4_NPARAMS+j]=sum/chol1[j*OLS4_NPARAMS+j];
					}
					sum=chol1[i*OLS4_NPARAMS+i];
					for(int k=0;k<i;++k)
						sum-=chol1[i*OLS4_NPARAMS+k]*chol1[i*OLS4_NPARAMS+k];
					if(sum<=1e-8)
					{
						success1=0;
						break;
					}
					chol1[i*OLS4_NPARAMS+i]=sqrt(sum);
				}
				for(int i=0;i<OLS4_NPARAMS;++i)
				{
					for(int j=0;j<i;++j)
					{
						sum=chol2[i*OLS4_NPARAMS+j];
						for(int k=0;k<j;++k)
							sum-=chol2[i*OLS4_NPARAMS+k]*chol2[j*OLS4_NPARAMS+k];
						chol2[i*OLS4_NPARAMS+j]=sum/chol2[j*OLS4_NPARAMS+j];
					}
					sum=chol2[i*OLS4_NPARAMS+i];
					for(int k=0;k<i;++k)
						sum-=chol2[i*OLS4_NPARAMS+k]*chol2[i*OLS4_NPARAMS+k];
					if(sum<=1e-8)
					{
						success2=0;
						break;
					}
					chol2[i*OLS4_NPARAMS+i]=sqrt(sum);
				}
				if(success0)
				{
					for(int i=0;i<OLS4_NPARAMS;++i)
					{
						sum=vec0[i];
						for(int j=0;j<i;++j)
							sum-=chol0[i*OLS4_NPARAMS+j]*params0[j];
						params0[i]=sum/chol0[i*OLS4_NPARAMS+i];
					}
					for(int i=OLS4_NPARAMS-1;i>=0;--i)
					{
						sum=params0[i];
						for(int j=i+1;j<OLS4_NPARAMS;++j)
							sum-=chol0[j*OLS4_NPARAMS+i]*params0[j];
						params0[i]=sum/chol0[i*OLS4_NPARAMS+i];
					}
				}
				if(success1)
				{
					for(int i=0;i<OLS4_NPARAMS;++i)
					{
						sum=vec1[i];
						for(int j=0;j<i;++j)
							sum-=chol1[i*OLS4_NPARAMS+j]*params1[j];
						params1[i]=sum/chol1[i*OLS4_NPARAMS+i];
					}
					for(int i=OLS4_NPARAMS-1;i>=0;--i)
					{
						sum=params1[i];
						for(int j=i+1;j<OLS4_NPARAMS;++j)
							sum-=chol1[j*OLS4_NPARAMS+i]*params1[j];
						params1[i]=sum/chol1[i*OLS4_NPARAMS+i];
					}
				}
				if(success2)
				{
					for(int i=0;i<OLS4_NPARAMS;++i)
					{
						sum=vec2[i];
						for(int j=0;j<i;++j)
							sum-=chol2[i*OLS4_NPARAMS+j]*params2[j];
						params2[i]=sum/chol2[i*OLS4_NPARAMS+i];
					}
					for(int i=OLS4_NPARAMS-1;i>=0;--i)
					{
						sum=params2[i];
						for(int j=i+1;j<OLS4_NPARAMS;++j)
							sum-=chol2[j*OLS4_NPARAMS+i]*params2[j];
						params2[i]=sum/chol2[i*OLS4_NPARAMS+i];
					}
				}
			}
			if(kols>=kols2)
				kols2+=OLS_STRIDE;
#endif
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
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	if(args->fwd)
		ac3_enc_flush(&ec);
}
int c03_codec(const char *srcfn, const char *dstfn)
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
	double bestsize;
	int usize;
#ifdef ABAC_PROFILESIZE
	double abac_csizes[ABAC_TOKEN_BITS*3]={0};
#endif
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
	bestsize=0;
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
		statssize=0;
		histsize=(int)sizeof(int[OCH_COUNT*PRED_COUNT<<8]);
		if(histsize<statssize)
			histsize=statssize;
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=image;
		arg->dst=fwd?0:dst->data+printed;
		arg->iw=iw;
		arg->ih=ih;
		arg->bufsize=sizeof(short[4*OCH_COUNT*2*(BLOCKSIZE+16LL)]);//4 padded rows * OCH_COUNT * {pixels, wg_errors}
		arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));

		arg->histsize=histsize;
		arg->hist=(int*)malloc(histsize);
		//arg->stats=(unsigned*)malloc(statssize);
#ifdef ENABLE_OLS
		for(int k=0;k<3;++k)
			ols4_init(arg->ols4+k, OLS4_NPARAMS);
#endif
		if(!arg->pixels||!arg->hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=arg->bufsize;
		memusage+=histsize;
		//memusage+=statssize;
		
		arg->tlevels=tlevels;
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
#ifdef DISABLE_MT
			for(int k=0;k<nthreads2;++k)
				block_thread(args+k);
#else
			void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
			mt_finish(ctx);
#endif
			if(fwd)
			{
				for(int kt2=0;kt2<nthreads2;++kt2)
				{
					ThreadArgs *arg=args+kt2;
					if(test)
					{
						int res=(arg->x2-arg->x1)*(arg->y2-arg->y1);
						int blocksize=(res*nch*depth+7)>>3;
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
						bestsize+=arg->bestsize;
#ifdef ABAC_PROFILESIZE
						for(int k=0;k<ABAC_TOKEN_BITS*3;++k)
							abac_csizes[k]+=arg->abac_csizes[k];
#endif
#ifdef ESTIMATE_SIZE
						for(int kc=0;kc<3;++kc)
						{
							for(int ks=0;ks<256;++ks)
							{
								int freq=arg->hist2[kc][ks];
								if(freq)
									esize[kc]-=freq*log2((double)freq/res);
							}
						}
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
#ifdef ABAC_PROFILESIZE
				double upsize=(double)iw*ih/8;
				for(int kb=0;kb<ABAC_TOKEN_BITS*3;++kb)
				{
					double cpsize=abac_csizes[kb]/8;
					printf("C%d B%d  %10.2lf/%10.2lf  %8.4lf%%  CR %10.6lf\n",
						kb/ABAC_TOKEN_BITS,
						kb%ABAC_TOKEN_BITS,
						cpsize,
						upsize,
						100.*cpsize/upsize,
						upsize/cpsize
					);
				}
#endif
#ifdef ESTIMATE_SIZE
				esize[0]/=8;
				esize[1]/=8;
				esize[2]/=8;
				printf("T %15.2lf\n", esize[0]+esize[1]+esize[2]);
				printf("Y %15.2lf\n", esize[0]);
				printf("U %15.2lf\n", esize[1]);
				printf("V %15.2lf\n", esize[2]);
#endif
				printf("Best %15.2lf (%+13.2lf) bytes\n", bestsize, csize-bestsize);
				printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
				printf("Mem usage: ");
				print_size((double)memusage, 8, 4, 0, 0);
				printf("\n");
			}
			printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
			if(!fwd)
				compare_bufs_8(image2, src->data+headersize, iw, ih, nch, nch, CODECNAME, 0, 1);
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
		//free(arg->stats);
#ifdef ENABLE_OLS
		for(int k=0;k<3;++k)
			ols4_free(arg->ols4+k);
#endif
	}
	free(args);
	array_free(&src);
	array_free(&dst);
	//if(test)
	//	pause();
	return 0;
}