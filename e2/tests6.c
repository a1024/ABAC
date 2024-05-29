#include"e2.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctype.h>
//#ifdef _MSC_VER
//#include<intrin.h>
//#else
//#include<x86intrin.h>
//#endif
static const char file[]=__FILE__;

//debug
//	#define ENABLE_GUIDE

//experiments
//	#define PRINT_HIST
//	#define PROFILER//SLOW
//	#define TRACK_SSE_RANGES//SLOW

//efficiency
//	#define ENABLE_SSE_4D//3% smaller, but 2x slower
//	#define USE_ABAC//inferior
//	#define USE_ABAC_SSE
//	#define SSE_UPDATE_NB_CELLS//inferior
//	#define COMPENSATE_sRGB//inferior
//	#define ENABLE_LZ//currently inferior
//	#define ENABLE_SSE_PAR//inferior

//memory
//	#define PROBBITS_15//less efficient

//speed
//	#define AVX512
//	#define AVX2

#define AC_IMPLEMENTATION
#include"ac.h"
#if defined AVX2 || defined AVX512
#include<immintrin.h>
#endif
#ifdef PROFILER
#define CHECKPOINTLIST\
	CHECKPOINT(RCT_OPT)\
	CHECKPOINT(OUTSIDE)\
	CHECKPOINT(FETCH_NB)\
	CHECKPOINT(CALC_SUBPREDS)\
	CHECKPOINT(CALC_WEIGHT_AV)\
	CHECKPOINT(CALC_CTX)\
	CHECKPOINT(SSE_LOOP)\
	CHECKPOINT(PRED_TILL_UPDATE)\
	CHECKPOINT(UPDATE_ERRORS)\
	CHECKPOINT(UPDATE_WP)\
	CHECKPOINT(UPDATE_HIST)\
	CHECKPOINT(UPDATE_SSE)
typedef enum ProfilerLabelEnum
{
#define CHECKPOINT(X) PROF_##X,
	CHECKPOINTLIST
#undef  CHECKPOINT
	PROF_COUNT,
} ProfilerLabel;
static const char *prof_labels[]=
{
#define CHECKPOINT(X) #X,
	CHECKPOINTLIST
#undef  CHECKPOINT
};
static long long prof_timestamp=0, prof_cycles[PROF_COUNT]={0};
#define PROF_START() memset(prof_cycles, 0, sizeof(prof_cycles)), prof_timestamp=__rdtsc()
#define PROF(X) prof_cycles[PROF_##X]+=__rdtsc()-prof_timestamp, prof_timestamp=__rdtsc()
static void prof_print()
{
	long long sum=0;
	int maxlen=0;
	for(int k=0;k<_countof(prof_labels);++k)
	{
		int len=(int)strlen(prof_labels[k]);
		UPDATE_MAX(maxlen, len);
		sum+=prof_cycles[k];
	}
	printf("Profiler:\n");
	for(int k=0;k<_countof(prof_labels);++k)
	{
		double percent=100.*prof_cycles[k]/sum;
		printf("%-*s %16lld %6.2lf%%  ", maxlen, prof_labels[k], prof_cycles[k], percent);
		for(int k2=0, npoints=(int)percent;k2<npoints;++k2)
			printf("*");
		printf("\n");
	}
}
#else
#define PROF_START()
#define PROF(...)
#define prof_print()
#endif

#ifdef COMPENSATE_sRGB
static int linear2sRGB(int x)
{
	x+=0x800000;
	x=CLAMP(1, x, 0xFFFFFF);
	if(x<=0xCD2E)
		x=(int)((long long)x*0xCEB851F>>24);
	else
	{
		x=(int)POW_FIX24(x, 0x6AAAAB);
		x=(int)((long long)x*0x10E147B>>24)-0xE147B;
	}
	x=CLAMP(1, x, 0xFFFFFF);
	x-=0x800000;
	return x;
}
static int sRGB2linear(int x)
{
	x+=0x800000;
	x=CLAMP(1, x, 0xFFFFFF);
	if(x<=0xA5AED)
		x=(int)((long long)x*0x13D072>>24);
	else
	{
		x=(int)((((long long)x+0xE147B)<<24)/0x10E147B);
		x=(int)POW_FIX24(x, 0x2666666);
	}
	x=CLAMP(1, x, 0xFFFFFF);
	x-=0x800000;
	return x;
}
#endif

#ifdef COMPENSATE_GAMMA
//these functions expect and return fix24 values
static int nonlinear2linear(int x, int gamma)
{
	x+=0x800000;
	x=CLAMP(1, x, 0xFFFFFF);
	x=log2_fix24(x);
	x=(int)(((long long)x<<24)/gamma);
	x=(int)exp2_fix24(x);
	x=CLAMP(1, x, 0xFFFFFF);
	x-=0x800000;
	return x;
}
static int linear2nonlinear(int x, int gamma)
{
	x+=0x800000;
	x=CLAMP(1, x, 0xFFFFFF);
	x=log2_fix24(x);
	x=(int)((long long)x*gamma>>24);
	x=(int)exp2_fix24(x);
	x=CLAMP(1, x, 0xFFFFFF);
	x-=0x800000;
	return x;
}
static int nonlinear2nonlinear(int x, int srcgamma, int dstgamma)
{
	//dst_nonlinear = src_nonlinear^(dst_gamma/src_gamma)
	x+=0x800000;
	x=CLAMP(1, x, 0xFFFFFF);
	x=log2_fix24(x);
	x=(int)((long long)x*dstgamma/srcgamma);
	x=(int)exp2_fix24(x);
	x=CLAMP(1, x, 0xFFFFFF);
	x-=0x800000;
	return x;
}
#endif

#define SLIC5_CONFIG_EXP 5
#define SLIC5_CONFIG_MSB 2
#define SLIC5_CONFIG_LSB 0

//#define CDF_UPDATE_PERIOD 0x400//number of sub-pixels processed between CDF updates, must be a power-of-two
#define PAD_SIZE 4
#define HIST_EXP 2
#define HIST_MSB 1
#define SSE_X_EXP 1
#define SSE_X_MSB 0
#define SSE_Y_EXP 1
#define SSE_Y_MSB 0
#define SSE_Z_EXP 1
#define SSE_Z_MSB 0
#define SSE_P_EXP 1
#define SSE_P_MSB 0
#define SSE_W 7		//ENABLE SSE INDEX CLAMP WHEN CHANGING CONFIG
#define SSE_H 21
#define SSE_D 21
#define SSE_PREDBITS 5
#define SSE_PRED_LEVELS (1<<SSE_PREDBITS)
#define SSE_FR_SIZE 59049//3^10		separate final round
#define SSE_STAGES 10
#ifdef SSE_UPDATE_NB_CELLS
#define SSE_STEP 5
#else
#define SSE_STEP 1
#endif
#define SSE_LIMIT (640*SSE_STEP)
#ifdef ENABLE_SSE_PAR
#define SSE_SIZE ((SSE_W+SSE_H+SSE_D)<<SSE_PREDBITS)
#elif defined ENABLE_SSE_4D
#define SSE_SIZE (SSE_W*SSE_H*SSE_D<<SSE_PREDBITS)
#endif
//#define PRED_PREC 8
#define PARAM_PREC 8
//#ifdef USE_ABAC
//#define NPAR_CTX 30
//#endif

//	#define ENABLE_CUSTOM1
//	#define ENABLE_CUSTOM1_v2
//	#define ENABLE_CUSTOM1_v3//best		0.1% smaller, but 2x slower
//	#define ENABLE_CUSTOM1_v4//incomplete

//don't forget to update SLIC5_NPREDS in e2.h
#if 0
#define SLIC5_PREDLIST\
	SLIC5_PRED(W+NE-N-((2*(eN+eW)+eNE-eNW+4)>>3))\
	SLIC5_PRED(N-(int)(((long long)eN+eW+eNE)*-0x05C>>PARAM_PREC))\
	SLIC5_PRED(W-(int)(((long long)eN+eW+eNW)*-0x05B>>PARAM_PREC))\
	SLIC5_PRED(N+(int)((-eNN*0x0DFLL-eN*0x0051LL-eNE*0x0BDLL+((long long)N-NN)*0x05C+((long long)NW-W)*0x0102)>>PARAM_PREC))\
	SLIC5_PRED((W+NEE)>>1)\
	SLIC5_PRED((4*N-2*NN+NW+NE)>>2)\
	SLIC5_PRED(N+NE-NNE-eNNE)\
	SLIC5_PRED((N+W)>>1)\
	SLIC5_PRED((4*(N+W+NW+NE)-(NN+WW+NNWW+NNEE)+6)/12)\
	SLIC5_PRED(N+W+NNWW-NNW-NWW+((eN+eW)>>1))\
	SLIC5_PRED(N+W-NW)\
	SLIC5_PRED(2*W-WW+eW-eWW)\
	SLIC5_PRED(paper_GAP)\
	SLIC5_PRED(calic_GAP)\
	SLIC5_PRED(ols)
#else
#define SLIC5_PREDLIST\
	SLIC5_PRED(W+NE-N-((2*(eN+eW)+eNE-eNW+4)>>3))\
	SLIC5_PRED(N-(int)(((long long)eN+eW+eNE)*-0x05C>>PARAM_PREC))\
	SLIC5_PRED(W-(int)(((long long)eN+eW+eNW)*-0x05B>>PARAM_PREC))\
	SLIC5_PRED(N+(int)((-eNN*0x0DFLL-eN*0x051LL-eNE*0x0BDLL+((long long)N-NN)*0x05C+((long long)NW-W)*0x102)>>PARAM_PREC))\
	SLIC5_PRED(3*(N-NN)+NNN)\
	SLIC5_PRED((N+W)>>1)\
	SLIC5_PRED(geomean)\
	SLIC5_PRED(N+W-NW)\
	SLIC5_PRED((W+NEE)>>1)\
	SLIC5_PRED((3*W+NEEE)>>2)\
	SLIC5_PRED((3*(3*W+NE+NEE)-10*N+2)/5)\
	SLIC5_PRED((3*(3*W+NE+NEE)-10*N)/5)\
	SLIC5_PRED((4*N-2*NN+NW+NE)>>2)\
	SLIC5_PRED(N+NE-NNE-eNNE)\
	SLIC5_PRED((4*(N+W+NW+NE)-(NN+WW+NNWW+NNEE)+6)/12)\
	SLIC5_PRED(W+((eW-eWW)>>1))\
	SLIC5_PRED(paper_GAP)\
	SLIC5_PRED(calic_GAP)\
	SLIC5_PRED(N+W-((NW+NN+WW+NE)>>2))\
	SLIC5_PRED(((2*(N+W)-(NW+NN+WW+NE))*9+(WWW+NWW+NNW+NNN+NNE+NEE)*2)/12)\
	SLIC5_PRED(3*(N+W-NW-(NN+WW-NNWW))+NNN+WWW-NNNWWW)\
	SLIC5_PRED(2*(W+NE-N)-(WW+NNEE-NN))\
	SLIC5_PRED((2*W+NEE-N)>>1)\
	SLIC5_PRED(NW+NWW-NNWWW)\
	SLIC5_PRED((14*NE-(NNEE+NNNEE+NNEEE))/11)\
	SLIC5_PRED((NEEE+NEEEE)>>1)\
	SLIC5_PRED((NNNEEEE+NNEEE)>>1)\
	SLIC5_PRED(NNEEEE)\
	SLIC5_PRED((NNWWWW+NNNWWWW)>>1)\
	SLIC5_PRED((WWW+WWWW)>>1)\
	SLIC5_PRED((N+NN)>>1)\
	SLIC5_PRED((NE+NNEE)>>1)\
	SLIC5_PRED((NE+NNE+NEE+NNEE)>>2)
//	SLIC5_PRED(ols)
//	SLIC5_PRED((4*(N+W+NE+NW)-(WW+NWW+NNWW+NNW+NN+NNE+NNEE+NEE))>>3)
//	SLIC5_PRED((9*(NE-NNEE)+NNNNEE+NNNEEE+NNEEEE)/3)
//	SLIC5_PRED(NN+NNW-(NNNN+NNNNW+NNNNWW)/3)
//	SLIC5_PRED(4*(N+NNN)-6*NN-(NNNNW+NNNN+NNNNE)/3)
//	SLIC5_PRED(ols)
//	SLIC5_PRED(N+W+NNWW-NNW-NWW+((eN+eW)>>1))
//	SLIC5_PRED((NNN+WWW+3*(N+W-NN-WW)+1)>>1)
//	SLIC5_PRED((WW+NE)>>1)
//	SLIC5_PRED((W+2*NE-NNE+eNE)>>1)
//	SLIC5_PRED(experimental)
//	SLIC5_PRED(floor_sqrt(floor_sqrt((N+0x800000LL)*(W+0x800000LL))*floor_sqrt((NW+0x800000LL)*(NE+0x800000LL)))-0x800000)
#endif

const char *slic5_prednames[]=
{
#define SLIC5_PRED(X) #X,
	SLIC5_PREDLIST
#undef  SLIC5_PRED
};
#ifdef PROBBITS_15
#define EC_ENC(EC, X, CDF) ac_enc15(EC, X, CDF)
#define EC_DEC(EC, CDF, NLEVELS) ac_dec15(EC, CDF, NLEVELS)
typedef unsigned short CDF_t;
#define CDF_SHIFT 15
#else
#define EC_ENC(EC, X, CDF) ac_enc(EC, X, CDF)
#define EC_DEC(EC, CDF, NLEVELS) ac_dec(EC, CDF, NLEVELS)
typedef unsigned CDF_t;
#define CDF_SHIFT 16
#endif

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
typedef struct HybridUintStruct
{
	unsigned short token, nbits;
	unsigned bypass;
} HybridUint;
static void hybriduint_encode(unsigned val, HybridUint *hu)
{
	int token, bypass, nbits;
	if(val<(1<<SLIC5_CONFIG_EXP))
	{
		token=val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=FLOOR_LOG2((unsigned)val);
		int mantissa=val-(1<<lgv);
		token = (1<<SLIC5_CONFIG_EXP) + (
				(lgv-SLIC5_CONFIG_EXP)<<(SLIC5_CONFIG_MSB+SLIC5_CONFIG_LSB)|
				(mantissa>>(lgv-SLIC5_CONFIG_MSB))<<SLIC5_CONFIG_LSB|
				(mantissa&((1<<SLIC5_CONFIG_LSB)-1))
			);
		nbits=lgv-(SLIC5_CONFIG_MSB+SLIC5_CONFIG_LSB);
		bypass=val>>SLIC5_CONFIG_LSB&((1LL<<nbits)-1);
	}
	hu->token=(unsigned short)token;
	hu->nbits=(unsigned short)nbits;
	hu->bypass=bypass;
}
INLINE int quantize_unsigned(int val, int exp, int msb)
{
	if(val<(1<<exp))
		return val;
	int lgv=FLOOR_LOG2(val);
	int token=(1<<exp)+((lgv-exp)<<msb|(val-(1<<lgv))>>(lgv-msb));
	return token;
}
INLINE int quantize_signed_get_range(int num, int den, int exp, int msb)
{
	int vmax=(int)((1LL<<24)*num/den>>16);
	int token=quantize_unsigned(vmax, exp, msb);
	token<<=1;
	return token;
}
INLINE int quantize_signed(int val, int shift, int exp, int msb, int nlevels)
{
	val>>=shift;
	int negmask=-(val<0);
	int token=quantize_unsigned(abs(val), exp, msb);
	//if((unsigned)token>=(unsigned)(nlevels>>1))//
	//	LOG_ERROR("SSE range error");
	token^=negmask;
	token-=negmask;
	token+=nlevels>>1;
	token=CLAMP(0, token, nlevels-1);
	return token;
}
#define QUANTIZE_HIST(X) (quantize_unsigned(X, HIST_EXP, HIST_MSB)>>1)

INLINE void matmul(double *dst, const double *m1, const double *m2, int h1, int w1h2, int w2)
{
	for(int ky=0;ky<h1;++ky)
	{
		for(int kx=0;kx<w2;++kx)
		{
			double sum=0;
			for(int j=0;j<w1h2;++j)
				sum+=m1[w1h2*ky+j]*m2[w2*j+kx];
			dst[w2*ky+kx]=sum;
		}
	}
}
INLINE void matmul_selftransposed(double *dst, const double *src, int mh, int mw, int dstw)
{
	for(int ky=0;ky<mh;++ky)
	{
		for(int kx=0;kx<mh;++kx)
		{
			double sum=0;
			for(int j=0;j<mw;++j)
				sum+=src[mw*ky+j]*src[mw*kx+j];
			//if((unsigned)(dstw*ky+kx)>=128)
			//	LOG_ERROR("");
			dst[dstw*ky+kx]=sum;
		}
	}
}
#if defined AVX2 || defined AVX512
INLINE void avx2_floor_log2_p1(__m256i *x)//floor_log2()+1
{
	//https://stackoverflow.com/questions/56153183/is-using-avx2-can-implement-a-faster-processing-of-lzcnt-on-a-word-array
#ifdef AVX512
	__m256i thirtytwo=_mm256_set1_epi32(32);
	*x=_mm256_lzcnt_epi32(*x);
	*x=_mm256_sub_epi32(thirtytwo, *x);
#else
	//_mm256_shuffle_epi8 has 4bit range
	__m256i ramplo=_mm256_set_epi8(
		//15 14 13  12  11  10   9   8   7   6   5   4   3   2   1   0
		 4,  4,  4,  4,  4,  4,  4,  4,  3,  3,  3,  3,  2,  2,  1,  0,
		 4,  4,  4,  4,  4,  4,  4,  4,  3,  3,  3,  3,  2,  2,  1,  0
	);
	__m256i ramphi=_mm256_set_epi8(
		//15 14 13  12  11  10   9   8   7   6   5   4   3   2   1   0
		 8,  8,  8,  8,  8,  8,  8,  8,  7,  7,  7,  7,  6,  6,  5,  0,
		 8,  8,  8,  8,  8,  8,  8,  8,  7,  7,  7,  7,  6,  6,  5,  0
	);
	__m256i masklo=_mm256_set1_epi32(15);

	__m256i lo=_mm256_and_si256(*x, masklo);
	__m256i hi=_mm256_srli_epi32(*x, 4);
	lo=_mm256_shuffle_epi8(ramplo, lo);
	hi=_mm256_shuffle_epi8(ramphi, hi);
	*x=_mm256_max_epi32(lo, hi);
#endif
}
#endif
#ifdef USE_ABAC
typedef struct StatNodeStruct
{
	unsigned short n[2], rec[6], weight[7];
#ifdef USE_ABAC_SSE
	unsigned short sse_ctr;
	long long sse_sum;
#endif
} StatNode;
#endif
typedef struct SLIC5CtxStruct
{
	int iw, ih, nch;
	char depths[4];
	int
		maxdepth, maxlevels,
		nlevels[4],
		half[4],
		shift_prec[4],//shift right to bring value to 8-bit, removes predictor-added bits
		rounding_offset[4],
		shift_error[4],
		full_pixels_processed;
	int *pred_errors, *errors, *pixels;
	
#ifdef USE_ABAC
	//int qplevels;
	int xtextures, ytextures, ntextures;
	//int maxtoken, treedepth;//treesize=1<<treedepth, index 0 is unused
	int treesize;
	StatNode *stats;
	unsigned short ctx_weights[SLIC5_NPREDS<<2];
	//int qpreds[30];
	int txid;
#else
	int nhist, cdfsize;
	int *hist, *histsums;
	CDF_t *CDFs;
	int hist_idx;
#endif

#ifdef ENABLE_SSE_PAR
	int sse_w, sse_h, sse_d, sse_size;
	long long *sse, *sse_fr;
	int sse_idx[SSE_STAGES*3+1];
	long long sse_sum[SSE_STAGES*3+1];
	int sse_count[SSE_STAGES*3+1];
#endif
#ifdef ENABLE_SSE_4D
	int sse_w, sse_h, sse_d, sse_p, sse_size;
	long long *sse, *sse_fr, *sse_cfl;
	int sse_idx_x[SSE_STAGES+6], sse_idx_y[SSE_STAGES+6], sse_idx_z[SSE_STAGES+6], sse_idx_p[SSE_STAGES+6], sse_idx[SSE_STAGES+6];
	long long sse_sum[SSE_STAGES+6];
	int sse_count[SSE_STAGES+6];
#endif
#ifdef ENABLE_CUSTOM1
	char *c1_params;
#endif
#ifdef ENABLE_CUSTOM1_v3
	double c1_params[4*4];
	int c1_init[4];
#endif
#ifdef ENABLE_CUSTOM1_v4
	double c1_params[8*4];
	int c1_init[4];
#endif
	int bias_count[4];
	long long bias_sum[4];

	int preds[SLIC5_NPREDS], params[SLIC5_NPREDS<<2];
	long long pred;
	int pred_final;
	int sse_corr;
	int kc, kx, ky, kym[4*4];
	long long pred_error_sums[SLIC5_NPREDS];
#ifdef TRACK_SSE_RANGES
	int sse_ranges[8];
#endif
	ArithmeticCoder *ec;
} SLIC5Ctx;
#define LOAD(BUF, X, Y) BUF[pr->kym[kc<<2|Y]+kx+PAD_SIZE-(X)]
#define LOAD_CH(BUF, C, X, Y) BUF[pr->kym[C<<2|Y]+kx+PAD_SIZE-(X)]
#define LOAD_PRED_ERROR(C, X, Y, P) pr->pred_errors[SLIC5_NPREDS*(pr->kym[C<<2|Y]+kx+PAD_SIZE-(X))+P]
//#define LOAD_PRED_ERROR(C, X, Y, P) pr->pred_errors[(SLIC5_NPREDS*(pr->kym[Y]+kx+PAD_SIZE-(X))+P)<<2|C]
static ptrdiff_t slic5_init(SLIC5Ctx *pr, int iw, int ih, int nch, const char *depths, ArithmeticCoder *ec)//returns memory usage or 0 on failure
{
	if(iw<1||ih<1||nch<1)
	{
		LOG_ERROR("Invalid image");
		return 0;
	}
	memset(pr, 0, sizeof(*pr));
	pr->iw=iw;
	pr->ih=ih;
	pr->nch=nch;
#ifdef SLIC5_OPTIMIZE_RCT
	memcpy(pr->depths, depths, nch);
#else
	if(nch<3)
		memcpy(pr->depths, depths, nch);
	else
	{
#ifdef COMPENSATE_sRGB
		pr->depths[0]=depths[1]+1;	//Y
#else
		pr->depths[0]=depths[1];	//Y
#endif
		pr->depths[1]=depths[2]+1;	//Cb
		pr->depths[2]=depths[0]+1;	//Cr
		if(nch==4)
			pr->depths[3]=depths[3];//a
	}
#endif
	pr->maxdepth=0;
	for(int kc=0;kc<nch;++kc)
	{
		pr->nlevels[kc]=1<<pr->depths[kc];
		//int prec_half=(1<<(pr->depths[kc]+PRED_PREC))>>1;
		UPDATE_MAX(pr->maxdepth, pr->depths[kc]);
	}
	pr->maxlevels=1<<pr->maxdepth;
	for(int kc=0;kc<pr->nch;++kc)
	{
		pr->shift_prec[kc]=24-pr->depths[kc];
		pr->half[kc]=pr->nlevels[kc]>>1;
		pr->rounding_offset[kc]=(1<<pr->shift_prec[kc])>>1;
		//pr->shift_prec[kc]=MAXVAR(8, pr->depths[kc])-8+PRED_PREC;
		pr->shift_error[kc]=1+((pr->depths[kc]<=9)<<1);
	}

	//suppose	curr == -128	pred > 0	rare but happens
	//error = -128-pred
	//abs(error) = 128+pred
	//upred = 128-pred
	//since abs(error)>upred:
	//	sym = abs(error)+upred = 128+pred+128-pred = 256
	HybridUint hu;
	hybriduint_encode(pr->maxlevels, &hu);

	//int extremesym=-(pr->maxlevels>>1);
	//extremesym=extremesym<<1^-(extremesym<0);//pack sign
	//++extremesym;
	//hybriduint_encode(extremesym, &hu);//encode -half

#ifdef USE_ABAC
	pr->treesize=pr->maxdepth*(pr->maxdepth+1)>>1;
	//pr->maxtoken=hu.token;
	//++hu.token;
	//pr->treedepth=ceil_log2(hu.token);
	//pr->treesize=1<<pr->treedepth;

	//pr->qplevels=quantize_signed_get_range(128, 1, 2, 1);
	pr->xtextures=QUANTIZE_HIST(767);
	pr->ytextures=QUANTIZE_HIST(767);
	pr->ntextures=pr->xtextures*pr->ytextures;
#else
	pr->cdfsize=hu.token+1;
#ifdef ENABLE_LZ
	++pr->cdfsize;//make way for LZ escape symbol
#endif
	pr->nhist=QUANTIZE_HIST(767);
#endif
	
#ifdef ENABLE_SSE_PAR
	pr->sse_w=quantize_unsigned(128, SSE_X_EXP, SSE_X_MSB)<<1;
	pr->sse_h=quantize_unsigned(128, SSE_Y_EXP, SSE_Y_MSB)<<1;
	pr->sse_d=quantize_unsigned(128, SSE_Z_EXP, SSE_Z_MSB)<<1;
	pr->sse_size=pr->sse_w*pr->sse_h*pr->sse_d;
#elif defined ENABLE_SSE_4D
	pr->sse_w=quantize_unsigned(128/32, SSE_X_EXP, SSE_X_MSB)<<1;
	pr->sse_h=quantize_unsigned(128*2, SSE_Y_EXP, SSE_Y_MSB)<<1;
	pr->sse_d=quantize_unsigned(128*2, SSE_Z_EXP, SSE_Z_MSB)<<1;
	pr->sse_p=quantize_unsigned((1<<SSE_PREDBITS)>>1, SSE_P_EXP, SSE_P_MSB)<<1;
#endif

	//pr->sse_width=7;
	//pr->sse_height=21;
	//pr->sse_depth=21;
	//pr->sse_width=quantize_signed_get_range(2, 64, SSE_X_EXP, SSE_X_MSB);
	//pr->sse_height=quantize_signed_get_range(2, 1, SSE_Y_EXP, SSE_Y_MSB);
	//pr->sse_depth=quantize_signed_get_range(2, 1, SSE_Z_EXP, SSE_Z_MSB);
	//pr->sse_nplanes=1<<SSE_PREDBITS;
	//pr->sse_planesize=pr->sse_width*pr->sse_height*pr->sse_depth;
	//pr->sse_size=pr->sse_nplanes*pr->sse_planesize;
	
#ifdef TRACK_SSE_RANGES
	pr->sse_ranges[0]=pr->sse_width>>1;
	pr->sse_ranges[1]=pr->sse_width>>1;
	pr->sse_ranges[2]=pr->sse_height>>1;
	pr->sse_ranges[3]=pr->sse_height>>1;
	pr->sse_ranges[4]=pr->sse_depth>>1;
	pr->sse_ranges[5]=pr->sse_depth>>1;
	pr->sse_ranges[6]=pr->sse_nplanes>>1;
	pr->sse_ranges[7]=pr->sse_nplanes>>1;
#endif

	size_t memusage=0, size;
#define ALLOC(PTR, PTR_TYPE, SIZE) size=SIZE, memusage+=size, PTR=(PTR_TYPE)malloc(size)
	ALLOC(pr->pred_errors, int*, (iw+PAD_SIZE*2LL)*sizeof(int[SLIC5_NPREDS*4*4]));
	ALLOC(pr->errors, int*, (iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	ALLOC(pr->pixels, int*, (iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	//pr->pred_errors=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[SLIC5_NPREDS*4*4]));//NPREDS * 4 rows * 4 comps
	//pr->errors=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	//pr->pixels=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
#ifdef USE_ABAC
	ALLOC(pr->stats, StatNode*, (size_t)pr->nch*pr->ntextures*pr->treesize*sizeof(StatNode));
	//ALLOC(pr->stats, StatNode*, (size_t)pr->nch*pr->qplevels*sizeof(StatNode[NPAR_CTX])<<pr->treedepth);
	//ALLOC(pr->stats, StatNode*, (size_t)pr->nch*pr->qplevels*sizeof(StatNode[SLIC5_NPREDS])<<pr->treedepth);
	//ALLOC(pr->stats, StatNode*, (size_t)pr->nch*pr->ntextures*sizeof(StatNode)<<pr->treedepth);
#else
	int total_hist=pr->nhist*pr->nhist;
	ALLOC(pr->hist, int*, (size_t)nch*total_hist*pr->cdfsize*sizeof(int));
	ALLOC(pr->histsums, int*, (size_t)nch*total_hist*sizeof(int));
	ALLOC(pr->CDFs, CDF_t*, (size_t)nch*total_hist*(pr->cdfsize+1LL)*sizeof(CDF_t));
	//pr->hist=(int*)malloc((size_t)nch*total_hist*pr->cdfsize*sizeof(int));//WH: cdfsize * NHIST
	//pr->histsums=(int*)malloc((size_t)nch*total_hist*sizeof(int));
	//pr->CDFs=(CDF_t*)malloc((size_t)nch*total_hist*(pr->cdfsize+1LL)*sizeof(CDF_t));
#endif
#ifdef ENABLE_SSE_PAR
	pr->sse=(long long*)malloc((size_t)nch*pr->sse_size*sizeof(long long[SSE_STAGES]));
	pr->sse_fr=(long long*)malloc(nch*sizeof(long long[SSE_FR_SIZE]));
#elif defined ENABLE_SSE_4D
	ALLOC(pr->sse, long long*, nch*sizeof(long long[SSE_SIZE*SSE_STAGES]));
	ALLOC(pr->sse_fr, long long*, nch*sizeof(long long[SSE_FR_SIZE]));
	//pr->sse=(long long*)malloc(nch*sizeof(long long[SSE_SIZE*SSE_STAGES]));
	//pr->sse_fr=(long long*)malloc(nch*sizeof(long long[SSE_FR_SIZE]));
	if(pr->nch>1)
	{
		ALLOC(pr->sse_cfl, long long*, (nch-1LL)*sizeof(long long[SSE_SIZE<<1]));
		//pr->sse_cfl=(long long*)malloc((nch-1LL)*sizeof(long long[SSE_SIZE<<1]));
		if(!pr->sse_cfl)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
	}
#endif
#ifdef ENABLE_CUSTOM1
	ALLOC(pr->c1_params, char, (size_t)nch*iw*sizeof(char[8]));
#endif
#undef  ALLOC
	if(!pr->pred_errors||!pr->errors||!pr->pixels
#ifdef USE_ABAC
		||!pr->stats
#else
		||!pr->hist||!pr->histsums||!pr->CDFs
#endif
#if defined ENABLE_SSE_4D || defined ENABLE_SSE_PAR
		||!pr->sse||!pr->sse_fr
#endif
#ifdef ENABLE_CUSTOM1
		||!pr->c1_params
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(pr->pred_errors, 0, (iw+PAD_SIZE*2LL)*sizeof(int[SLIC5_NPREDS*4*4]));
	memset(pr->errors, 0, (iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	memset(pr->pixels, 0, (iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	
#ifdef USE_ABAC
	StatNode node=
	{
		{1, 1},
		{0x8000, 0x8000, 0x8000, 0x8000, 0x8000, 0x8000},
		{0x8000, 0x8000, 0x8000, 0x8000, 0x8000, 0x8000, 0x8000},
	};
	memfill(pr->stats, &node, (size_t)pr->nch*pr->ntextures*pr->treesize*sizeof(StatNode), sizeof(StatNode));
	//memfill(pr->stats, &node, (size_t)pr->nch*pr->qplevels*sizeof(StatNode[NPAR_CTX])<<pr->treedepth, sizeof(StatNode));
	//memfill(pr->stats, &node, (size_t)pr->nch*pr->qplevels*sizeof(StatNode[SLIC5_NPREDS])<<pr->treedepth, sizeof(StatNode));
	//memfill(pr->stats, &node, (size_t)pr->nch*pr->ntextures*sizeof(StatNode)<<pr->treedepth, sizeof(StatNode));
	for(int k=0;k<_countof(pr->ctx_weights);++k)
		pr->ctx_weights[k]=0x8000;
#else
	memset(pr->hist, 0, (size_t)nch*total_hist*pr->cdfsize*sizeof(int));
	memset(pr->histsums, 0, (size_t)nch*total_hist*sizeof(int));

#if 0
	for(int qy=0;qy<pr->nhist;++qy)
	{
		for(int qx=0;qx<pr->nhist;++qx)
		{
			int hist_idx=pr->nhist*qy+qx;
			CDF_t *curr_CDF=pr->CDFs+(pr->cdfsize+1)*hist_idx;
			int sum=0, c=0;
			//int num=16+10*MINVAR(qx, qy);
			for(int ks=0;ks<pr->cdfsize;++ks)
			{
				//int freq=num/(ks+1);
				int freq=0x10000/(ks+1);
				curr_CDF[ks]=freq;
				sum+=freq;
			}
			for(int ks=0;ks<pr->cdfsize;++ks)
			{
				int freq=curr_CDF[ks];
				//curr_CDF[ks]=(int)((long long)c*((1LL<<CDF_SHIFT)-pr->cdfsize)/sum)+ks;//guard
				curr_CDF[ks]=(int)((long long)c*(1LL<<CDF_SHIFT)/sum);//no guard
				c+=freq;
			}
			curr_CDF[pr->cdfsize]=1<<CDF_SHIFT;
		}
	}
#endif
	int sum=0;
	//int freq=0xA186;
	for(int ks=0;ks<pr->cdfsize;++ks)//TODO tune initial CDFs
	{
		int freq=0x10000/(ks+1);
		pr->CDFs[ks]=freq;
		sum+=freq;
		//freq=(int)((long long)freq*freq>>16);
	}
	for(int ks=0, c=0;ks<pr->cdfsize;++ks)
	{
		int freq=pr->CDFs[ks];
		pr->CDFs[ks]=(int)((long long)c*((1LL<<CDF_SHIFT)-pr->cdfsize)/sum)+ks;
		c+=freq;
	}
	pr->CDFs[pr->cdfsize]=1<<CDF_SHIFT;
	memfill(pr->CDFs+pr->cdfsize+1, pr->CDFs, ((size_t)nch*total_hist-1LL)*(pr->cdfsize+1LL)*sizeof(CDF_t), (pr->cdfsize+1LL)*sizeof(CDF_t));
#endif
	
#ifdef ENABLE_SSE_PAR
	memset(pr->sse, 0, (size_t)nch*pr->sse_size*sizeof(long long[SSE_STAGES]));
	memset(pr->sse_fr, 0, nch*sizeof(long long[SSE_FR_SIZE]));
#elif defined ENABLE_SSE_4D
	memset(pr->sse, 0, nch*sizeof(long long[SSE_SIZE*SSE_STAGES]));
	memset(pr->sse_fr, 0, nch*sizeof(long long[SSE_FR_SIZE]));
	if(pr->nch>1)
		memset(pr->sse_cfl, 0, (nch-1LL)*sizeof(long long[SSE_SIZE<<1]));
#endif
#ifdef ENABLE_CUSTOM1
	memset(pr->c1_params, 0, (size_t)nch*iw*sizeof(char[8]));
#endif

	for(int k=0;k<(SLIC5_NPREDS<<2);++k)
		pr->params[k]=176;
	//	pr->params[k]=12<<pr->maxdepth;
	//int shift=16-maxdepth;
	//shift=MAXVAR(0, shift);
	//memcpy(pr->params, wp_params, sizeof(pr->params));
	//for(int k=0;k<SLIC5_NPREDS;++k)
	//	pr->params[k]=pr->params[k]>>shift;

	pr->ec=ec;
	return memusage;
}
static void slic5_free(SLIC5Ctx *pr)
{
	free(pr->pred_errors);
	free(pr->errors);
	free(pr->pixels);
#ifdef USE_ABAC
	free(pr->stats);
#else
	free(pr->hist);
	free(pr->histsums);
	free(pr->CDFs);
#endif
#if defined ENABLE_SSE_4D || defined ENABLE_SSE_PAR
	free(pr->sse);
	free(pr->sse_fr);
	if(pr->nch>1)
		free(pr->sse_cfl);
#endif
#ifdef ENABLE_CUSTOM1
	free(pr->c1_params);
#endif
}
#if 0
static void slic5_update_CDFs(SLIC5Ctx *pr)
{
	int nhist=pr->nch*pr->nhist*pr->nhist;
	for(int kt=0;kt<nhist;++kt)//update CDFs (deferred)
	{
		int sum=pr->histsums[kt];
		if(sum>pr->cdfsize)//only when hist was hit enough times
		{
			int *curr_hist=pr->hist+pr->cdfsize*kt;
			CDF_t *curr_CDF=pr->CDFs+(pr->cdfsize+1)*kt;
			long long c=0;
			for(int ks=0;ks<pr->cdfsize;++ks)
			{
				long long freq=curr_hist[ks];
				curr_CDF[ks]=(int)(c*((1LL<<CDF_SHIFT)-pr->cdfsize)/sum)+ks;
				c+=freq;
			}
			curr_CDF[pr->cdfsize]=1<<CDF_SHIFT;
		}
	}
}
#endif
static void slic5_nextrow(SLIC5Ctx *pr, int ky)
{
	pr->ky=ky;
	pr->kym[0]=(pr->iw+PAD_SIZE*2)*(ky&3);
	pr->kym[1]=(pr->iw+PAD_SIZE*2)*((ky-1)&3);
	pr->kym[2]=(pr->iw+PAD_SIZE*2)*((ky-2)&3);
	pr->kym[3]=(pr->iw+PAD_SIZE*2)*((ky-3)&3);
	for(int k=4;k<16;++k)
		pr->kym[k]=pr->kym[k-4]+4*(pr->iw+PAD_SIZE*2);
#if defined ENABLE_CUSTOM1_v3 || defined ENABLE_CUSTOM1_v4
	memset(pr->c1_init, 0, sizeof(pr->c1_init));
#endif

	//int nhist=pr->nch*pr->nhist*pr->nhist;
	//for(int kt=0;kt<nhist;++kt)//force rescale histograms
	//{
	//}
}
#ifdef ENABLE_CUSTOM1
static int custom1_predict(const int *nb, const char *params)
{
	ALIGN(16) long long p[2];
	memcpy(p, params, sizeof(*p));
	p[1]=0;
	__m256i val=_mm256_load_si256((__m256i*)nb);
	__m256i par=_mm256_cvtepi8_epi32(_mm_load_si128((__m128i*)p));
	val=_mm256_srai_epi32(val, 7);
	val=_mm256_mullo_epi32(val, par);
	val=_mm256_hadd_epi32(val, val);
	val=_mm256_hadd_epi32(val, val);
	int sum=(_mm256_extract_epi32(val, 0)+_mm256_extract_epi32(val, 4))>>3;
	return sum;
	//__m128i v2=_mm256_extracti128_si256(_mm256_hadd_epi32(val, val), 0);
	//v2=_mm_hadd_epi32(v2, v2);
	//long long v3=_mm_extract_epi64(v2, 0);
	//v3=(int)v3+(int)(v3>>32);
	//return (int)v3;
}
#endif
#ifdef ENABLE_CUSTOM1_v2
static double g_matrix[16*8];
#endif
INLINE int invert_matrix(double *matrix, int size, double *temprow)
{
	int success=1;
	//double *temp=_malloca(size*sizeof(double[2]));
	for(int it=0;it<size;++it)
	{
		int kp=it;
		for(;kp<size;++kp)
		{
			if(fabs(matrix[((size_t)size<<1)*kp+it])>1e-6)
				break;
		}
		if(kp==size)
		{
			success=0;
			break;
		}
		if(kp!=it)
		{
			memcpy(temprow, matrix+((size_t)size<<1)*it, size*sizeof(double[2]));
			memcpy(matrix+((size_t)size<<1)*it, matrix+((size_t)size<<1)*kp, size*sizeof(double[2]));
			memcpy(matrix+((size_t)size<<1)*kp, temprow, size*sizeof(double[2]));
		}
		double pivot=matrix[((size_t)size<<1)*it+it];
		for(int kx=it;kx<(size<<1);++kx)
			matrix[((size_t)size<<1)*it+kx]/=pivot;
		for(int ky=0;ky<size;++ky)
		{
			if(ky==it)
				continue;
			double factor=matrix[((size_t)size<<1)*ky+it];
			if(fabs(factor)>1e-6)
			{
				for(int kx=it;kx<(size<<1);++kx)
					matrix[((size_t)size<<1)*ky+kx]-=matrix[((size_t)size<<1)*it+kx]*factor;
			}
		}
	}
	//_freea(temp);
	return success;
}
static void slic5_predict(SLIC5Ctx *pr, int kc, int kx)
{
	PROF(OUTSIDE);
	//int idx=(pr->iw*pr->ky+kx)<<2|kc;
	pr->kc=kc;
	pr->kx=kx;
	pr->full_pixels_processed+=!kc;
	//if(!(idx&(CDF_UPDATE_PERIOD-1))||(idx<CDF_UPDATE_PERIOD&&(idx&15)))
	//if(!(kx&(CDF_UPDATE_PERIOD-1)))
	//	slic5_update_CDFs(pr);
	//PROF(UPDATE_CDFs);
	//XY are flipped, no need to check if indices OOB due to padding
	int
#if PAD_SIZE>=4
	//	NNNNWW	=LOAD(pr->pixels,  2, 4),
	//	NNNNW	=LOAD(pr->pixels,  1, 4),
	//	NNNN	=LOAD(pr->pixels,  0, 4),
	//	NNNNE	=LOAD(pr->pixels, -1, 4),
	//	NNNNEE	=LOAD(pr->pixels, -2, 4),
	//	NNNNEEEE=LOAD(pr->pixels, -4, 4),

		NNNWWWW	=LOAD(pr->pixels,  4, 3),
		NNNWWW	=LOAD(pr->pixels,  3, 3),
	//	NNNWW	=LOAD(pr->pixels,  2, 3),
	//	NNNW	=LOAD(pr->pixels,  1, 3),
		NNN	=LOAD(pr->pixels,  0, 3),
	//	NNNE	=LOAD(pr->pixels, -1, 3),
		NNNEE	=LOAD(pr->pixels, -2, 3),
	//	NNNEEE	=LOAD(pr->pixels, -3, 3),
		NNNEEEE	=LOAD(pr->pixels, -4, 3),
		
		NNWWWW	=LOAD(pr->pixels,  4, 2),
		NNWWW	=LOAD(pr->pixels,  3, 2),
#endif
		NNWW	=LOAD(pr->pixels,  2, 2),
		NNW	=LOAD(pr->pixels,  1, 2),
		NN	=LOAD(pr->pixels,  0, 2),
		NNE	=LOAD(pr->pixels, -1, 2),
		NNEE	=LOAD(pr->pixels, -2, 2),
#if PAD_SIZE>=4
		NNEEE	=LOAD(pr->pixels, -3, 2),
		NNEEEE	=LOAD(pr->pixels, -4, 2),
	//	NWWWW	=LOAD(pr->pixels,  4, 1),
	//	NWWW	=LOAD(pr->pixels,  3, 1),
#endif
		NWW	=LOAD(pr->pixels,  2, 1),
		NW	=LOAD(pr->pixels,  1, 1),
		N	=LOAD(pr->pixels,  0, 1),
		NE	=LOAD(pr->pixels, -1, 1),
		NEE	=LOAD(pr->pixels, -2, 1),
#if PAD_SIZE>=4
		NEEE	=LOAD(pr->pixels, -3, 1),
		NEEEE	=LOAD(pr->pixels, -4, 1),
		WWWW	=LOAD(pr->pixels,  4, 0),
		WWW	=LOAD(pr->pixels,  3, 0),
#endif
		WW	=LOAD(pr->pixels,  2, 0),
		W	=LOAD(pr->pixels,  1, 0);
	int
	//	eNNWW	=LOAD(pr->errors,  2, 2),//error = (curr<<8) - pred
	//	eNNW	=LOAD(pr->errors,  1, 2),
		eNN	=LOAD(pr->errors,  0, 2),
		eNNE	=LOAD(pr->errors, -1, 2),
	//	eNNEE	=LOAD(pr->errors, -2, 2),
	//	eNWW	=LOAD(pr->errors,  2, 1),
		eNW	=LOAD(pr->errors,  1, 1),
		eN	=LOAD(pr->errors,  0, 1),
		eNE	=LOAD(pr->errors, -1, 1),
	//	eNEE	=LOAD(pr->errors, -2, 1),
		eWW	=LOAD(pr->errors,  2, 0),
		eW	=LOAD(pr->errors,  1, 0);
	int sh=24-8;
	//int sh=pr->shift_prec[kc];
	int clamp_lo=(int)N, clamp_hi=(int)N;
	clamp_lo=(int)MINVAR(clamp_lo, W);
	clamp_hi=(int)MAXVAR(clamp_hi, W);
	clamp_lo=(int)MINVAR(clamp_lo, NE);
	clamp_hi=(int)MAXVAR(clamp_hi, NE);
	PROF(FETCH_NB);

	int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
	int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
	int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
	int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
	int diff=(dy-dx)>>sh, diff2=(d45-d135)>>sh, diff3=NE-NW;
	int paper_GAP, calic_GAP;
	if(dy+dx>32)
		paper_GAP=(int)(((long long)dx*N+(long long)dy*W)/((long long)dy+dx));
	else if(diff>12)
		paper_GAP=(N+2*W)/3;
	else if(diff<-12)
		paper_GAP=(2*N+W)/3;
	else
		paper_GAP=(N+W)>>1;

	if(diff2>32)
		paper_GAP+=diff3>>2;
	else if(diff2>16)
		paper_GAP+=diff3*3>>4;
	else if(diff2>=-16)
		paper_GAP+=diff3>>3;
	else if(diff2>=-32)
		paper_GAP+=diff3>>4;

	if(diff>80)
		calic_GAP=W;
	else if(diff<-80)
		calic_GAP=N;
	else if(diff>32)
		calic_GAP=(2*N+6*W+NE-NW)>>3;		//c1	[1/4  3/4  1/8  -1/8].[N W NE NW]
	else if(diff>8)
		calic_GAP=(6*N+10*W+3*(NE-NW))>>4;	//c2	[3/8  5/8  3/16  -3/16]
	else if(diff<-32)
		calic_GAP=(6*N+2*W+NE-NW)>>3;		//c3	[3/4  1/4  1/8  -1/8]
	else if(diff<-8)
		calic_GAP=(10*N+6*W+3*(NE-NW))>>4;	//c4	[5/8  3/8  3/16  -3/16]
	else
		calic_GAP=(((N+W)<<1)+NE-NW)>>2;	//c5	[1/2  1/2  1/4  -1/4]

	int geomean;
	long long aN=N+0x800000LL, aW=W+0x800000LL;
	aN=CLAMP(0, aN, 0xFFFFFF);
	aW=CLAMP(0, aW, 0xFFFFFF);
	geomean=(int)floor_sqrt(aN*aW)-0x800000;

#ifdef ENABLE_CUSTOM1_v4
	int ols=0;
	//NNNWWWW NNNWWW NNNWW NNNW NNN  NNNE NNNEE NNNEEE NNNEEEE
	//NNWWWW  NNWWW  NNWW  NNW  NN   NNE  NNEE  NNEEE  NNEEEE
	//NWWWW   NWWW   NWW   NW   N    NE   NEE   NEEE   NEEEE
	//WWWW    WWW    WW    W    curr
#define OLS_NPARAMS 8
#define OLS_NSAMPLES 17
	ALIGN(32) double nb[(OLS_NPARAMS+1)*OLS_NSAMPLES];
	for(int ky2=-2, idx=0;ky2<=0;++ky2)
	{
		for(int kx2=-3;kx2<=3;++kx2)
		{
			if(!ky2&&!kx2)
				break;
			nb[OLS_NSAMPLES*0+idx]=pr->pixels[(pr->kym[-(ky2-1)]+kx+4-kx2-1)<<2|kc];//NW
			nb[OLS_NSAMPLES*1+idx]=pr->pixels[(pr->kym[-(ky2-1)]+kx+4-kx2+0)<<2|kc];//N
			nb[OLS_NSAMPLES*2+idx]=pr->pixels[(pr->kym[-(ky2-1)]+kx+4-kx2+1)<<2|kc];//NE
			nb[OLS_NSAMPLES*3+idx]=pr->pixels[(pr->kym[-(ky2+0)]+kx+4-kx2-1)<<2|kc];//W
			nb[OLS_NSAMPLES*4+idx]=pr->errors[(pr->kym[-(ky2-1)]+kx+4-kx2-1)<<2|kc];//eNW
			nb[OLS_NSAMPLES*5+idx]=pr->errors[(pr->kym[-(ky2-1)]+kx+4-kx2+0)<<2|kc];//eN
			nb[OLS_NSAMPLES*6+idx]=pr->errors[(pr->kym[-(ky2-1)]+kx+4-kx2+1)<<2|kc];//eNE
			nb[OLS_NSAMPLES*7+idx]=pr->errors[(pr->kym[-(ky2+0)]+kx+4-kx2-1)<<2|kc];//eW
			nb[OLS_NSAMPLES*8+idx]=pr->pixels[(pr->kym[-(ky2+0)]+kx+4-kx2+0)<<2|kc];//target
		}
	}
	for(int k=0;k<_countof(nb);++k)
		nb[k]/=1<<24;
	ALIGN(32) double sqm[OLS_NPARAMS*(OLS_NPARAMS<<1)]={0}, temp[OLS_NPARAMS<<1];
	for(int k=0;k<OLS_NPARAMS;++k)//set identity
		sqm[(OLS_NPARAMS<<1)*k+k+OLS_NPARAMS]=1;
	matmul_selftransposed(sqm, nb, OLS_NPARAMS, OLS_NSAMPLES, OLS_NPARAMS<<1);//sqm = NB*NBT
	int success=invert_matrix(sqm, OLS_NPARAMS, temp);//inv(NB*NBT)
	if(success)
	{
		//pred = nb * params,		params += ((inv(NB * NBT) * NB * Y) - params)*lr
		matmul(temp, nb, nb+OLS_NSAMPLES*OLS_NPARAMS, OLS_NPARAMS, OLS_NSAMPLES, 1);
		matmul(nb, sqm, temp, OLS_NPARAMS, OLS_NPARAMS, 1);
		if(pr->c1_init[kc])
		{
			for(int k=0;k<OLS_NPARAMS;++k)
				pr->c1_params[OLS_NPARAMS*kc]+=(nb[k]-pr->c1_params[OLS_NPARAMS*kc])*0.5;
		}
		else
		{
			memcpy(pr->c1_params+OLS_NPARAMS*kc, nb, sizeof(double[OLS_NPARAMS]));
			pr->c1_init[kc]=1;
		}
		int nb3[]=
		{
			NW,	N,	NE,	W,	eNW,	eN,	eNE,	eW,
		};
		double fpred=0;
		for(int k=0;k<OLS_NPARAMS;++k)
			fpred+=nb3[k]*pr->c1_params[kc<<2|k];
		ols=(int)fpred;
	}
#undef  OLS_NSAMPLES
#endif
#ifdef ENABLE_CUSTOM1_v3
	int ols=0;
	//NNNWWWW NNNWWW NNNWW NNNW NNN  NNNE NNNEE NNNEEE NNNEEEE
	//NNWWWW  NNWWW  NNWW  NNW  NN   NNE  NNEE  NNEEE  NNEEEE
	//NWWWW   NWWW   NWW   NW   N    NE   NEE   NEEE   NEEEE
	//WWWW    WWW    WW    W    curr
#define OLS_NPARAMS 4
#define OLS_NSAMPLES 17
	int nb[]=
	{
		NNNWWWW,NNNWWW,	NNNWW,	NNWWWW,	//NNWWW
		NNNWWW,	NNNWW,	NNNW,	NNWWW,	//NNWW
		NNNWW,	NNNW,	NNN,	NNWW,	//NNW
		NNNW,	NNN,	NNNE,	NNW,	//NN
		NNN,	NNNE,	NNNEE,	NN,	//NNE
		NNNE,	NNNEE,	NNNEEE,	NNE,	//NNEE
		NNNEE,	NNNEEE,	NNNEEEE,NNEE,	//NNEEE
		NNWWWW,	NNWWW,	NNWW,	NWWWW,	//NWWW
		NNWWW,	NNWW,	NNW,	NWWW,	//NWW
		NNWW,	NNW,	NN,	NWW,	//NW
		NNW,	NN,	NNE,	NW,	//N
		NN,	NNE,	NNEE,	N,	//NE
		NNE,	NNEE,	NNEEE,	NE,	//NEE
		NNEE,	NNEEE,	NNEEEE,	NEE,	//NEEE
		NWWWW,	NWWW,	NWW,	WWWW,	//WWW
		NWWW,	NWW,	NW,	WWW,	//WW
		NWW,	NW,	N,	WW,	//W
	};
	const int priority[]=
	{
		1, 1, 1, 1, 1, 1, 1,
		1, 1, 2, 4, 2, 1, 1,
		1, 1, 4,
		
		//1, 1, 1, 1, 1, 1, 1,
		//1, 1, 1, 3, 1, 1, 1,
		//1, 1, 3,

		//1, 2, 4,  8, 4, 2, 1,
		//2, 4, 8, 16, 8, 4, 2,
		//4, 8, 16,
	};
	double nb2[_countof(nb)];
	for(int ky=0;ky<4;++ky)
	{
		for(int kx=0;kx<OLS_NSAMPLES;++kx)
			nb2[OLS_NSAMPLES*ky+kx]=(double)nb[kx<<2|ky]/(1<<24);
	}
	double sqm[OLS_NPARAMS*(OLS_NPARAMS<<1)]={0};
	for(int ky=0;ky<OLS_NPARAMS;++ky)
	{
		for(int kx=0;kx<OLS_NPARAMS;++kx)
		{
			double sum=0;
			for(int j=0;j<OLS_NSAMPLES;++j)
				sum+=priority[j]*nb2[OLS_NSAMPLES*ky+j]*nb2[OLS_NSAMPLES*kx+j];
			sqm[(OLS_NPARAMS<<1)*ky+kx]=sum;
		}
		sqm[(OLS_NPARAMS<<1)*ky+ky+OLS_NPARAMS]=1;
	}
	double temp[OLS_NPARAMS<<1];
	int success=invert_matrix(sqm, OLS_NPARAMS, temp);
	if(success||pr->c1_init[kc])
	{
		//pred = nb * params,		params += ((inv(NB * NBT) * NB * Y) - params)*lr
		if(success)
		{
			int targets[]=
			{
				NNWWW,
				NNWW,
				NNW,
				NN,
				NNE,
				NNEE,
				NNEEE,
				NWWW,
				NWW,
				NW,
				N,
				NE,
				NEE,
				NEEE,
				WWW,
				WW,
				W,
			};
			for(int ky=0;ky<OLS_NPARAMS;++ky)
			{
				double sum=0;
				for(int kx=0;kx<OLS_NSAMPLES;++kx)
					sum+=priority[kx]*nb2[OLS_NSAMPLES*ky+kx]*targets[kx]/(1<<24);
				temp[ky]=sum;
			}
			for(int ky=0;ky<OLS_NPARAMS;++ky)
			{
				double sum=0;
				for(int j=0;j<OLS_NPARAMS;++j)
					sum+=sqm[(OLS_NPARAMS<<1)*ky+j+OLS_NPARAMS]*temp[j];
				if(pr->c1_init[kc])
					pr->c1_params[kc<<2|ky]+=(sum-pr->c1_params[kc<<2|ky])*0.2;
				else
					pr->c1_params[kc<<2|ky]=sum;
			}
			pr->c1_init[kc]=1;
		}
		int nb3[]=
		{
			NW,	N,	NE,	W,
		};
		double fpred=0;
		for(int k=0;k<OLS_NPARAMS;++k)
			fpred+=nb3[k]*pr->c1_params[kc<<2|k];
		ols=(int)fpred;
	}
#endif
#ifdef ENABLE_CUSTOM1_v2
	int ols=0;
	ALIGN(32) int nb[]=
	{
		//NNWWW NNWW NNW NN   NNE NNEE NNEEE
		//NWWW  NWW  NW  N    NE  NEE  NEEE
		//WWW   WW   W   curr

		//NW	N	NE	W
		NNWWW,	NNWW,	NNW,	NWWW,
		NNWW,	NNW,	NN,	NWW,
		NNW,	NN,	NNE,	NW,
		NN,	NNE,	NNEE,	N,
		NNE,	NNEE,	NNEEE,	NE,
		NWWW,	NWW,	NW,	WWW,
		NWW,	NW,	N,	WW,
		0x1000000, 0x1000000, 0x1000000, 0x1000000,

		//NNWWW,	NNWW,	NNW,	NN,	NNE,	NWWW,	NWW,	1,//NW
		//NNWW,		NNW,	NN,	NNE,	NNEE,	NWW,	NW,	1,//N
		//NNW,		NN,	NNE,	NNEE,	NNEEE,	NW,	N,	1,//NE
		//NWWW,		NWW,	NW,	N,	NE,	WWW,	WW,	1,//W
	};

	//FIXME adds and divs OVERFLOW, use mullo_epi32, srai_epi32
#if 0
	{
		__m256i shuf0=_mm256_set_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1, 13, 12, 9, 8, 5, 4, 1, 0,
			-1, -1, -1, -1, -1, -1, -1, -1, 13, 12, 9, 8, 5, 4, 1, 0
		);
		__m256i t0=_mm256_load_si256((__m256i*)nb+0);
		__m256i t1=_mm256_load_si256((__m256i*)nb+1);
		__m256i t2=_mm256_load_si256((__m256i*)nb+2);
		__m256i t3=_mm256_load_si256((__m256i*)nb+3);
		t0=_mm256_srai_epi32(t0, 10);//s7.24 -> s1.14 bit	mulhrs: 1.15*1.15>>15=s.15 or 1.14*1.14>>15 = s2.13
		t1=_mm256_srai_epi32(t1, 10);
		t2=_mm256_srai_epi32(t2, 10);
		t3=_mm256_srai_epi32(t3, 10);
		t0=_mm256_shuffle_epi8(t0, shuf0);
		t1=_mm256_shuffle_epi8(t1, shuf0);
		t2=_mm256_shuffle_epi8(t2, shuf0);
		t3=_mm256_shuffle_epi8(t3, shuf0);
		_mm256_store_si256((__m256i*)nb+0, t0);
		_mm256_store_si256((__m256i*)nb+1, t1);
		_mm256_store_si256((__m256i*)nb+2, t2);
		_mm256_store_si256((__m256i*)nb+3, t3);
	}
	__m128i row0=_mm_load_si128((__m128i*)nb+0);//[0, 0, 0, 0, a3, a2, a1, a0]	there is no _mm_madd_pi16 (MMX)
	__m128i row1=_mm_load_si128((__m128i*)nb+1);//[0, 0, 0, 0, b3, b2, b1, b0]
	__m128i row2=_mm_load_si128((__m128i*)nb+2);
	__m128i row3=_mm_load_si128((__m128i*)nb+3);
	__m128i row4=_mm_load_si128((__m128i*)nb+4);
	__m128i row5=_mm_load_si128((__m128i*)nb+5);
	__m128i row6=_mm_load_si128((__m128i*)nb+6);
	__m128i row7=_mm_load_si128((__m128i*)nb+7);
	ALIGN(32) short sqm[8*16]={0};
	__m128i t0, t1, t2, t3, t4, t5, t6, t7;
	//row0
	t0=_mm_madd_epi16(row0, row0);//[0, 0, a3*a3+a2*a2, a1*a1+a0*a0]
	t1=_mm_madd_epi16(row0, row1);//[0, 0, b3*b3+b2*b2, b1*b1+b0*b0]
	t2=_mm_madd_epi16(row0, row2);
	t3=_mm_madd_epi16(row0, row3);
	t4=_mm_madd_epi16(row0, row4);
	t5=_mm_madd_epi16(row0, row5);
	t6=_mm_madd_epi16(row0, row6);
	t7=_mm_madd_epi16(row0, row7);
	sqm[0<<4|0]=_mm_extract_epi32(_mm_hadd_epi32(t0, t0), 0)>>14;//[0, a3*a3+a2*a2+a1*a1+a0*a0, 0, a3*a3+a2*a2+a1*a1+a0*a0]
	sqm[0<<4|1]=_mm_extract_epi32(_mm_hadd_epi32(t1, t1), 0)>>14;
	sqm[0<<4|2]=_mm_extract_epi32(_mm_hadd_epi32(t2, t2), 0)>>14;
	sqm[0<<4|3]=_mm_extract_epi32(_mm_hadd_epi32(t3, t3), 0)>>14;
	sqm[0<<4|4]=_mm_extract_epi32(_mm_hadd_epi32(t4, t4), 0)>>14;
	sqm[0<<4|5]=_mm_extract_epi32(_mm_hadd_epi32(t5, t5), 0)>>14;
	sqm[0<<4|6]=_mm_extract_epi32(_mm_hadd_epi32(t6, t6), 0)>>14;
	sqm[0<<4|7]=_mm_extract_epi32(_mm_hadd_epi32(t7, t7), 0)>>14;

	//row1
	t1=_mm_madd_epi16(row1, row1);
	t2=_mm_madd_epi16(row1, row2);
	t3=_mm_madd_epi16(row1, row3);
	t4=_mm_madd_epi16(row1, row4);
	t5=_mm_madd_epi16(row1, row5);
	t6=_mm_madd_epi16(row1, row6);
	t7=_mm_madd_epi16(row1, row7);
	sqm[1<<4|0]=sqm[0<<4|1];
	sqm[1<<4|1]=_mm_extract_epi32(_mm_hadd_epi32(t1, t1), 0)>>14;
	sqm[1<<4|2]=_mm_extract_epi32(_mm_hadd_epi32(t2, t2), 0)>>14;
	sqm[1<<4|3]=_mm_extract_epi32(_mm_hadd_epi32(t3, t3), 0)>>14;
	sqm[1<<4|4]=_mm_extract_epi32(_mm_hadd_epi32(t4, t4), 0)>>14;
	sqm[1<<4|5]=_mm_extract_epi32(_mm_hadd_epi32(t5, t5), 0)>>14;
	sqm[1<<4|6]=_mm_extract_epi32(_mm_hadd_epi32(t6, t6), 0)>>14;
	sqm[1<<4|7]=_mm_extract_epi32(_mm_hadd_epi32(t7, t7), 0)>>14;

	//row2
	t2=_mm_madd_epi16(row2, row2);
	t3=_mm_madd_epi16(row2, row3);
	t4=_mm_madd_epi16(row2, row4);
	t5=_mm_madd_epi16(row2, row5);
	t6=_mm_madd_epi16(row2, row6);
	t7=_mm_madd_epi16(row2, row7);
	sqm[2<<4|0]=sqm[0<<4|2];
	sqm[2<<4|1]=sqm[1<<4|2];
	sqm[2<<4|2]=_mm_extract_epi32(_mm_hadd_epi32(t2, t2), 0)>>14;
	sqm[2<<4|3]=_mm_extract_epi32(_mm_hadd_epi32(t3, t3), 0)>>14;
	sqm[2<<4|4]=_mm_extract_epi32(_mm_hadd_epi32(t4, t4), 0)>>14;
	sqm[2<<4|5]=_mm_extract_epi32(_mm_hadd_epi32(t5, t5), 0)>>14;
	sqm[2<<4|6]=_mm_extract_epi32(_mm_hadd_epi32(t6, t6), 0)>>14;
	sqm[2<<4|7]=_mm_extract_epi32(_mm_hadd_epi32(t7, t7), 0)>>14;

	//row3
	t3=_mm_madd_epi16(row3, row3);
	t4=_mm_madd_epi16(row3, row4);
	t5=_mm_madd_epi16(row3, row5);
	t6=_mm_madd_epi16(row3, row6);
	t7=_mm_madd_epi16(row3, row7);
	sqm[3<<4|0]=sqm[0<<4|3];
	sqm[3<<4|1]=sqm[1<<4|3];
	sqm[3<<4|2]=sqm[2<<4|3];
	sqm[3<<4|3]=_mm_extract_epi32(_mm_hadd_epi32(t3, t3), 0)>>14;
	sqm[3<<4|4]=_mm_extract_epi32(_mm_hadd_epi32(t4, t4), 0)>>14;
	sqm[3<<4|5]=_mm_extract_epi32(_mm_hadd_epi32(t5, t5), 0)>>14;
	sqm[3<<4|6]=_mm_extract_epi32(_mm_hadd_epi32(t6, t6), 0)>>14;
	sqm[3<<4|7]=_mm_extract_epi32(_mm_hadd_epi32(t7, t7), 0)>>14;

	//row4
	t4=_mm_madd_epi16(row4, row4);
	t5=_mm_madd_epi16(row4, row5);
	t6=_mm_madd_epi16(row4, row6);
	t7=_mm_madd_epi16(row4, row7);
	sqm[4<<4|0]=sqm[0<<4|4];
	sqm[4<<4|1]=sqm[1<<4|4];
	sqm[4<<4|2]=sqm[2<<4|4];
	sqm[4<<4|3]=sqm[3<<4|4];
	sqm[4<<4|4]=_mm_extract_epi32(_mm_hadd_epi32(t4, t4), 0)>>14;
	sqm[4<<4|5]=_mm_extract_epi32(_mm_hadd_epi32(t5, t5), 0)>>14;
	sqm[4<<4|6]=_mm_extract_epi32(_mm_hadd_epi32(t6, t6), 0)>>14;
	sqm[4<<4|7]=_mm_extract_epi32(_mm_hadd_epi32(t7, t7), 0)>>14;

	//rows 5, 6, 7
	t5=_mm_madd_epi16(row5, row5);
	t6=_mm_madd_epi16(row5, row6);
	t7=_mm_madd_epi16(row5, row7);

	t0=_mm_madd_epi16(row6, row6);
	t1=_mm_madd_epi16(row6, row7);

	t2=_mm_madd_epi16(row7, row7);
	sqm[5<<4|0]=sqm[0<<4|5];
	sqm[5<<4|1]=sqm[1<<4|5];
	sqm[5<<4|2]=sqm[2<<4|5];
	sqm[5<<4|3]=sqm[3<<4|5];
	sqm[5<<4|4]=sqm[4<<4|5];
	sqm[5<<4|5]=_mm_extract_epi32(_mm_hadd_epi32(t5, t5), 0)>>14;
	sqm[5<<4|6]=_mm_extract_epi32(_mm_hadd_epi32(t6, t6), 0)>>14;
	sqm[5<<4|7]=_mm_extract_epi32(_mm_hadd_epi32(t7, t7), 0)>>14;

	sqm[6<<4|0]=sqm[0<<4|6];
	sqm[6<<4|1]=sqm[1<<4|6];
	sqm[6<<4|2]=sqm[2<<4|6];
	sqm[6<<4|3]=sqm[3<<4|6];
	sqm[6<<4|4]=sqm[4<<4|6];
	sqm[6<<4|5]=sqm[5<<4|6];
	sqm[6<<4|6]=_mm_extract_epi32(_mm_hadd_epi32(t0, t0), 0)>>14;
	sqm[6<<4|7]=_mm_extract_epi32(_mm_hadd_epi32(t1, t1), 0)>>14;

	sqm[7<<4|0]=sqm[0<<4|7];
	sqm[7<<4|1]=sqm[1<<4|7];
	sqm[7<<4|2]=sqm[2<<4|7];
	sqm[7<<4|3]=sqm[3<<4|7];
	sqm[7<<4|4]=sqm[4<<4|7];
	sqm[7<<4|5]=sqm[5<<4|7];
	sqm[7<<4|6]=sqm[6<<4|7];
	sqm[7<<4|7]=_mm_extract_epi32(_mm_hadd_epi32(t2, t2), 0)>>14;
	int success=1;
	for(int it=0;it<8;++it)//Gaussian elimination
	{
		int kp=it;
		for(;kp<8;++kp)
		{
			if(sqm[kp<<4|it])
				break;
		}
		if(kp==8)
		{
			success=0;
			break;
		}
		if(kp!=it)//swap rows it & kp
		{
			__m256i v0=_mm256_load_si256((__m256i*)sqm+kp);
			__m256i v1=_mm256_load_si256((__m256i*)sqm+it);
			_mm256_store_si256((__m256i*)sqm+it, v0);
			_mm256_store_si256((__m256i*)sqm+kp, v1);
		}
		int pivot=sqm[it<<4|it];
		__m256i mp=_mm256_set1_epi16(pivot);
		__m256i pr=_mm256_load_si256((__m256i*)sqm+it);
		for(int ky=0;ky<8;++ky)
		{
			if(ky==it)
				continue;
			int factor=sqm[ky<<4|it];
			if(factor)
			{
				__m256i mf=_mm256_set1_epi16(factor);
				__m256i row=_mm256_load_si256((__m256i*)sqm+ky);
				__m256i t0=_mm256_mulhrs_epi16(row, mp);//s1.14 bit, mulhrs: 1.14*1.14>>15 = s2.13
				__m256i t1=_mm256_mulhrs_epi16(pr, mf);
				t0=_mm256_sub_epi16(t0, t1);
				t0=_mm256_slli_epi16(t0, 1);//s2.13 -> s.14
				_mm256_store_si256((__m256i*)sqm+ky, t0);
				//memset((__m256i*)sqm+ky, 0, ky*sizeof(short));

				//sqm[ky<<4|it]=0;
				//for(int kx=it+1;kx<16;++kx)
				//	sqm[ky<<4|kx]=((long long)sqm[ky<<4|kx]*pivot-(long long)sqm[it<<4|kx]*factor)>>16;
			}
		}
		sqm[it<<4|it]=1<<14;
		for(int kx=it+1;kx<16;++kx)
			sqm[it<<4|kx]=(sqm[it<<4|kx]<<14)/pivot;//(0x4000<<14)/1 = 0x1000[0000]	OVERFLOW
	}
	if(success)
	{
		//pred = nb * (inv(NB * NBT) * NB * Y)

		//multiply NB * Y
		__m256i shuf1=_mm256_set_epi8(
			13, 12, 9, 8, 5, 4, 1, 0, 13, 12, 9, 8, 5, 4, 1, 0,
			13, 12, 9, 8, 5, 4, 1, 0, 13, 12, 9, 8, 5, 4, 1, 0
		);
		__m256i Y=_mm256_set_epi32(W, NE, N, NW, W, NE, N, NW);
		Y=_mm256_srai_epi32(Y, 10);
		Y=_mm256_shuffle_epi8(Y, shuf1);//s1.14 [W, NE, N, NW] x4
		row0=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(row0), _mm_castsi128_ps(row1), _MM_SHUFFLE(1, 0, 1, 0)));
		row2=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(row2), _mm_castsi128_ps(row3), _MM_SHUFFLE(1, 0, 1, 0)));
		row4=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(row4), _mm_castsi128_ps(row5), _MM_SHUFFLE(1, 0, 1, 0)));
		row6=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(row6), _mm_castsi128_ps(row6), _MM_SHUFFLE(1, 0, 1, 0)));
		__m256i r0=_mm256_inserti128_si256(_mm256_castsi128_si256(row0), row2, 1);
		__m256i r1=_mm256_inserti128_si256(_mm256_castsi128_si256(row4), row6, 1);
		r0=_mm256_madd_epi16(r0, Y);//[d1, d0, c1, c0, b1, b0, a1, a0]
		r1=_mm256_madd_epi16(r1, Y);
		r0=_mm256_hadd_epi32(r0, r1);
		r0=_mm256_srai_epi32(r0, 14);//s1.14

		//multiply inv(NB*NBT) * (NB*Y) two elements at a time
		r0=_mm256_shuffle_epi8(r0, shuf1);
		__m256i col=_mm256_permute4x64_epi64(r0, _MM_SHUFFLE(2, 0, 2, 0));//repeat column twice
		__m256i r2, r3;
		r0=_mm256_inserti128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i*)sqm+ 1)), _mm_load_si128((__m128i*)sqm+ 3), 1);//8+8 int16's
		r1=_mm256_inserti128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i*)sqm+ 5)), _mm_load_si128((__m128i*)sqm+ 7), 1);
		r2=_mm256_inserti128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i*)sqm+ 9)), _mm_load_si128((__m128i*)sqm+11), 1);
		r3=_mm256_inserti128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i*)sqm+13)), _mm_load_si128((__m128i*)sqm+15), 1);
		r0=_mm256_madd_epi16(col, r0);//4+4 int32's
		r1=_mm256_madd_epi16(col, r1);
		r2=_mm256_madd_epi16(col, r2);
		r3=_mm256_madd_epi16(col, r3);
		r0=_mm256_hadd_epi32(r0, r1);//2+2+2+2 int32's
		r2=_mm256_hadd_epi32(r2, r3);
		r0=_mm256_hadd_epi32(r0, r2);//8 int32's	the params
		ALIGN(32) int custom_params[8];
		_mm256_store_si256((__m256i*)custom_params, r0);

		//NNWWW NNWW NNW NN   NNE NNEE NNEEE
		//NWWW  NWW  NW  N    NE  NEE  NEEE
		//WWW   WW   W   curr
		int nb3[]=
		{
			NWW,	NW,	N,	NE,	NEE,	WW,	W,	1,
		};
		long long pred1=0;
		for(int k=0;k<8;++k)
			pred1+=(long long)nb3[k]*custom_params[k];
		pred1+=1LL<<27;
		pred1>>=28;
		custom1=(int)pred1;


		//__m256i s0=_mm256_cvtepi16_epi32(_mm_load_si128((__m128i*)sqm+ 1));
		//__m256i s1=_mm256_cvtepi16_epi32(_mm_load_si128((__m128i*)sqm+ 3));
		//__m256i s2=_mm256_cvtepi16_epi32(_mm_load_si128((__m128i*)sqm+ 5));
		//__m256i s3=_mm256_cvtepi16_epi32(_mm_load_si128((__m128i*)sqm+ 7));
		//__m256i s4=_mm256_cvtepi16_epi32(_mm_load_si128((__m128i*)sqm+ 9));
		//__m256i s5=_mm256_cvtepi16_epi32(_mm_load_si128((__m128i*)sqm+11));
		//__m256i s6=_mm256_cvtepi16_epi32(_mm_load_si128((__m128i*)sqm+13));
		//__m256i s7=_mm256_cvtepi16_epi32(_mm_load_si128((__m128i*)sqm+15));
		//__m128i Y=_mm_set_epi16(0, 0, 0, 0, W, NE, N, NW);//half-SIMD
		//t0=_mm_madd_epi16(row0, Y);//[0, 0, a1, a0]
		//t1=_mm_madd_epi16(row1, Y);//[0, 0, b1, b0]
		//t2=_mm_madd_epi16(row2, Y);
		//t3=_mm_madd_epi16(row3, Y);
		//t4=_mm_madd_epi16(row4, Y);
		//t5=_mm_madd_epi16(row5, Y);
		//t6=_mm_madd_epi16(row6, Y);
		//t7=_mm_madd_epi16(row7, Y);
		//t0=_mm_hadd_epi32(t0, t0);//[0, a1+a0, 0, a1+a0]
		//t1=_mm_hadd_epi32(t1, t1);//[0, b1+b0, 0, b1+b0]
		//t2=_mm_hadd_epi32(t2, t2);
		//t3=_mm_hadd_epi32(t3, t3);
		//t4=_mm_hadd_epi32(t4, t4);
		//t5=_mm_hadd_epi32(t5, t5);
		//t6=_mm_hadd_epi32(t6, t6);
		//t7=_mm_hadd_epi32(t7, t7);
	}
#endif

	//FIXME numerically instable
#if 1
	ALIGN(32) int sqm[8*16]={0};
	__m128i c0=_mm_load_si128((__m128i*)nb+0);
	__m128i c1=_mm_load_si128((__m128i*)nb+1);
	__m128i c2=_mm_load_si128((__m128i*)nb+2);
	__m128i c3=_mm_load_si128((__m128i*)nb+3);
	__m128i c4=_mm_load_si128((__m128i*)nb+4);
	__m128i c5=_mm_load_si128((__m128i*)nb+5);
	__m128i c6=_mm_load_si128((__m128i*)nb+6);
	__m128i c7=_mm_load_si128((__m128i*)nb+7);
	c0=_mm_srai_epi32(c0, 10);//7.24 -> 17.14 bit
	c1=_mm_srai_epi32(c1, 10);
	c2=_mm_srai_epi32(c2, 10);
	c3=_mm_srai_epi32(c3, 10);
	c4=_mm_srai_epi32(c4, 10);
	c5=_mm_srai_epi32(c5, 10);
	c6=_mm_srai_epi32(c6, 10);
	c7=_mm_srai_epi32(c7, 10);
	__m128i val;
#define DOT_COL(Y, X) val=_mm_mullo_epi32(c##Y, c##X), val=_mm_srai_epi32(val, 14), val=_mm_hadd_epi32(val, val), val=_mm_hadd_epi32(val, val), sqm[Y<<4|X]=_mm_extract_epi32(val, 0)
#define CPY_RES(Y, X) sqm[Y<<4|X]=sqm[X<<4|Y]
	DOT_COL(0, 0),	DOT_COL(0, 1),	DOT_COL(0, 2),	DOT_COL(0, 3),	DOT_COL(0, 4),	DOT_COL(0, 5),	DOT_COL(0, 6),	DOT_COL(0, 7);
	CPY_RES(1, 0),	DOT_COL(1, 1),	DOT_COL(1, 2),	DOT_COL(1, 3),	DOT_COL(1, 4),	DOT_COL(1, 5),	DOT_COL(1, 6),	DOT_COL(1, 7);
	CPY_RES(2, 0),	CPY_RES(2, 1),	DOT_COL(2, 2),	DOT_COL(2, 3),	DOT_COL(2, 4),	DOT_COL(2, 5),	DOT_COL(2, 6),	DOT_COL(2, 7);
	CPY_RES(3, 0),	CPY_RES(3, 1),	CPY_RES(3, 2),	DOT_COL(3, 3),	DOT_COL(3, 4),	DOT_COL(3, 5),	DOT_COL(3, 6),	DOT_COL(3, 7);
	CPY_RES(4, 0),	CPY_RES(4, 1),	CPY_RES(4, 2),	CPY_RES(4, 3),	DOT_COL(4, 4),	DOT_COL(4, 5),	DOT_COL(4, 6),	DOT_COL(4, 7);
	CPY_RES(5, 0),	CPY_RES(5, 1),	CPY_RES(5, 2),	CPY_RES(5, 3),	CPY_RES(5, 4),	DOT_COL(5, 5),	DOT_COL(5, 6),	DOT_COL(5, 7);
	CPY_RES(6, 0),	CPY_RES(6, 1),	CPY_RES(6, 2),	CPY_RES(6, 3),	CPY_RES(6, 4),	CPY_RES(6, 5),	DOT_COL(6, 6),	DOT_COL(6, 7);
	CPY_RES(7, 0),	CPY_RES(7, 1),	CPY_RES(7, 2),	CPY_RES(7, 3),	CPY_RES(7, 4),	CPY_RES(7, 5),	CPY_RES(7, 6),	DOT_COL(7, 7);
	for(int k=0;k<8;++k)
		sqm[k<<4|(k+8)]=1<<14;
#undef  DOT_COL
#undef  CPY_RES
	//if(pr->ky==256&&kx==256)//
	//	printf("");

	int success=1;
	for(int k=0;k<_countof(sqm);++k)
		g_matrix[k]=(double)sqm[k]/(1<<14);
	for(int it=0;it<8;++it)
	{
		int kp=it;
		for(;kp<8;++kp)
		{
			if(fabs(g_matrix[kp<<4|it])>1e-4)
				break;
		}
		if(kp==8)
		{
			success=0;
			break;
		}
		if(kp!=it)
		{
			double temp[16];
			memcpy(temp, g_matrix+((size_t)it<<4), sizeof(double[16]));
			memcpy(g_matrix+((size_t)it<<4), g_matrix+((size_t)kp<<4), sizeof(double[16]));
			memcpy(g_matrix+((size_t)kp<<4), temp, sizeof(double[16]));
		}
		double pivot=g_matrix[it<<4|it];
		for(int kx=it;kx<16;++kx)
			g_matrix[it<<4|kx]/=pivot;
		for(int ky=0;ky<8;++ky)
		{
			if(ky==it)
				continue;
			double factor=g_matrix[ky<<4|it];
			if(fabs(factor)>1e-4)
			{
				for(int kx=it;kx<16;++kx)
					g_matrix[ky<<4|kx]-=g_matrix[it<<4|kx]*factor;
			}
		}
	}
	for(int k=0;k<_countof(sqm);++k)
	{
		if(isinf(g_matrix[k]))
		{
			success=0;
			break;
		}
		sqm[k]=(int)(g_matrix[k]*(1<<14));
	}
#if 0
	for(int it=0;it<8;++it)//Gaussian elimination
	{
		int kp=it;
		for(;kp<8;++kp)
		{
			if(sqm[kp<<4|it])
				break;
		}
		if(kp==8)
		{
			success=0;
			break;
		}
		if(kp!=it)//swap rows it & kp
		{
			__m256i v1, v2, v3, v4;
			v1=_mm256_load_si256((__m256i*)sqm+(kp<<1|0));
			v2=_mm256_load_si256((__m256i*)sqm+(kp<<1|1));
			v3=_mm256_load_si256((__m256i*)sqm+(it<<1|0));
			v4=_mm256_load_si256((__m256i*)sqm+(it<<1|1));
			_mm256_store_si256((__m256i*)sqm+(it<<1|0), v1);
			_mm256_store_si256((__m256i*)sqm+(it<<1|1), v2);
			_mm256_store_si256((__m256i*)sqm+(kp<<1|0), v3);
			_mm256_store_si256((__m256i*)sqm+(kp<<1|1), v4);
		}
		int pivot=sqm[it<<4|it];

		sqm[it<<4|it]=1<<14;//normalize pivot row
		for(int kx=it+1;kx<16;++kx)
			sqm[it<<4|kx]=((long long)sqm[it<<4|kx]<<14)/pivot;

		//__m256i mp=_mm256_set1_epi32(pivot);
		__m256i prow0=_mm256_load_si256((__m256i*)sqm+(it<<1|0));
		__m256i prow1=_mm256_load_si256((__m256i*)sqm+(it<<1|1));
		for(int ky=0;ky<8;++ky)
		{
			if(ky==it)
				continue;
			int factor=sqm[ky<<4|it];
			if(factor)
			{
				__m256i mf=_mm256_set1_epi32(factor);
				__m256i r0=_mm256_load_si256((__m256i*)sqm+(ky<<1|0));
				__m256i r1=_mm256_load_si256((__m256i*)sqm+(ky<<1|1));
				__m256i t2=_mm256_mullo_epi32(prow0, mf);
				__m256i t3=_mm256_mullo_epi32(prow1, mf);
				r0=_mm256_sub_epi32(r0, t2);
				r1=_mm256_sub_epi32(r1, t3);
				r0=_mm256_srai_epi32(r0, 14);//s2.13 -> s.14
				r1=_mm256_srai_epi32(r1, 14);
				_mm256_store_si256((__m256i*)sqm+(ky<<1|0), r0);
				_mm256_store_si256((__m256i*)sqm+(ky<<1|1), r1);

				//__m256i mf=_mm256_set1_epi32(factor);
				//__m256i r0=_mm256_load_si256((__m256i*)sqm+(ky<<1|0));
				//__m256i r1=_mm256_load_si256((__m256i*)sqm+(ky<<1|1));
				//__m256i t0=_mm256_mullo_epi32(r0, mp);//s1.14 bit, mulhrs: 1.14*1.14>>15 = s2.13
				//__m256i t1=_mm256_mullo_epi32(r1, mp);
				//__m256i t2=_mm256_mullo_epi32(prow0, mf);
				//__m256i t3=_mm256_mullo_epi32(prow1, mf);
				//t0=_mm256_sub_epi32(t0, t2);
				//t1=_mm256_sub_epi32(t1, t3);
				//t0=_mm256_srai_epi32(t0, 14);//s2.13 -> s.14
				//t1=_mm256_srai_epi32(t1, 14);
				//_mm256_store_si256((__m256i*)sqm+(ky<<1|0), t0);
				//_mm256_store_si256((__m256i*)sqm+(ky<<1|1), t1);

				//sqm[ky<<4|it]=0;
				//for(int kx=it+1;kx<16;++kx)
				//	sqm[ky<<4|kx]=((long long)sqm[ky<<4|kx]*pivot-(long long)sqm[it<<4|kx]*factor)>>14;
			}
		}
	}
#endif
	if(success)
	{
		//pred = nb * (inv(NB * NBT) * NB * Y)

		//multiply NB * Y two elements at a time
		__m128i temp=_mm_set_epi32(W, NE, N, NW);
		temp=_mm_srai_epi32(temp, 10);
		__m256i Y=_mm256_inserti128_si256(_mm256_castsi128_si256(temp), temp, 1);//duplicate column vector
		__m256i r0=_mm256_inserti128_si256(_mm256_castsi128_si256(c0), c1, 1);
		__m256i r1=_mm256_inserti128_si256(_mm256_castsi128_si256(c2), c3, 1);
		__m256i r2=_mm256_inserti128_si256(_mm256_castsi128_si256(c4), c5, 1);
		__m256i r3=_mm256_inserti128_si256(_mm256_castsi128_si256(c6), c7, 1);
		//__m256i r0=_mm256_load_si256((__m256i*)nb+0);//4+4 int32's
		//__m256i r1=_mm256_load_si256((__m256i*)nb+1);
		//__m256i r2=_mm256_load_si256((__m256i*)nb+2);
		//__m256i r3=_mm256_load_si256((__m256i*)nb+3);
		r0=_mm256_mullo_epi32(Y, r0);//4+4 int32's
		r1=_mm256_mullo_epi32(Y, r1);
		r2=_mm256_mullo_epi32(Y, r2);
		r3=_mm256_mullo_epi32(Y, r3);
		r0=_mm256_srai_epi32(r0, 14);
		r1=_mm256_srai_epi32(r1, 14);
		r2=_mm256_srai_epi32(r2, 14);
		r3=_mm256_srai_epi32(r3, 14);
		r0=_mm256_hadd_epi32(r0, r1);//2+2+2+2 int32's
		r2=_mm256_hadd_epi32(r2, r3);
		__m256i col=_mm256_hadd_epi32(r0, r2);//8 int32's

		__m256i r4, r5, r6, r7;
		r0=_mm256_load_si256((__m256i*)sqm+1);
		r1=_mm256_load_si256((__m256i*)sqm+3);
		r2=_mm256_load_si256((__m256i*)sqm+5);
		r3=_mm256_load_si256((__m256i*)sqm+7);
		r4=_mm256_load_si256((__m256i*)sqm+9);
		r5=_mm256_load_si256((__m256i*)sqm+11);
		r6=_mm256_load_si256((__m256i*)sqm+13);
		r7=_mm256_load_si256((__m256i*)sqm+15);
		r0=_mm256_mullo_epi32(col, r0);
		r1=_mm256_mullo_epi32(col, r1);
		r2=_mm256_mullo_epi32(col, r2);
		r3=_mm256_mullo_epi32(col, r3);
		r4=_mm256_mullo_epi32(col, r4);
		r5=_mm256_mullo_epi32(col, r5);
		r6=_mm256_mullo_epi32(col, r6);
		r7=_mm256_mullo_epi32(col, r7);
		r0=_mm256_srai_epi32(r0, 14);
		r1=_mm256_srai_epi32(r1, 14);
		r2=_mm256_srai_epi32(r2, 14);
		r3=_mm256_srai_epi32(r3, 14);
		r4=_mm256_srai_epi32(r4, 14);
		r5=_mm256_srai_epi32(r5, 14);
		r6=_mm256_srai_epi32(r6, 14);
		r7=_mm256_srai_epi32(r7, 14);
		r0=_mm256_hadd_epi32(r0, r1);
		r2=_mm256_hadd_epi32(r2, r3);
		r4=_mm256_hadd_epi32(r4, r5);
		r6=_mm256_hadd_epi32(r6, r7);
		r0=_mm256_hadd_epi32(r0, r2);
		r4=_mm256_hadd_epi32(r4, r6);
		r0=_mm256_hadd_epi32(r0, r4);
		ALIGN(32) int custom_params[8];
		_mm256_store_si256((__m256i*)custom_params, r0);

		//NNWWW NNWW NNW NN   NNE NNEE NNEEE
		//NWWW  NWW  NW  N    NE  NEE  NEEE
		//WWW   WW   W   curr
		int nb3[]=
		{
			NWW,	NW,	N,	NE,	NEE,	WW,	W,	1,
		};
		long long pred1=0;
		for(int k=0;k<8;++k)
			pred1+=(long long)nb3[k]*custom_params[k];
		pred1+=1LL<<13;
		pred1>>=14;
		custom1=(int)pred1;
#if 0
		__m256i Y=_mm256_set_epi64x(W, NE, N, NW);
		//ALIGN(32) int Y[]=
		//{
		//	NW,	0,
		//	N,	0,
		//	NE,	0,
		//	W,	0,
		//};
		//ALIGN(32) int temp[16];
		__m256i row0=_mm256_cvtepi32_epi64(_mm_load_si128((__m128i*)nb+0));
		__m256i row1=_mm256_cvtepi32_epi64(_mm_load_si128((__m128i*)nb+1));
		__m256i row2=_mm256_cvtepi32_epi64(_mm_load_si128((__m128i*)nb+2));
		__m256i row3=_mm256_cvtepi32_epi64(_mm_load_si128((__m128i*)nb+3));
		__m256i row4=_mm256_cvtepi32_epi64(_mm_load_si128((__m128i*)nb+4));
		__m256i row5=_mm256_cvtepi32_epi64(_mm_load_si128((__m128i*)nb+5));
		__m256i row6=_mm256_cvtepi32_epi64(_mm_load_si128((__m128i*)nb+6));
		__m256i row7=_mm256_cvtepi32_epi64(_mm_load_si128((__m128i*)nb+7));
		row0=_mm256_mul_epi32(row0, Y);//[a3, a2, a1, a0]
		row1=_mm256_mul_epi32(row1, Y);//[b3, b2, b1, b0]
		row2=_mm256_mul_epi32(row2, Y);//[...]
		row3=_mm256_mul_epi32(row3, Y);
		row4=_mm256_mul_epi32(row4, Y);
		row5=_mm256_mul_epi32(row5, Y);
		row6=_mm256_mul_epi32(row6, Y);
		row7=_mm256_mul_epi32(row7, Y);
		__m256i t0=_mm256_unpacklo_epi64(row0, row1);//[b2, a2, b0, a0]	there are no _mm256_hadd_epi64, _mm256_srai_epi64
		__m256i t1=_mm256_unpackhi_epi64(row0, row1);//[b3, a3, b1, a1]
		__m256i t2=_mm256_unpacklo_epi64(row2, row3);//[...]
		__m256i t3=_mm256_unpackhi_epi64(row2, row3);
		__m256i t4=_mm256_unpacklo_epi64(row4, row5);
		__m256i t5=_mm256_unpackhi_epi64(row4, row5);
		__m256i t6=_mm256_unpacklo_epi64(row6, row7);
		__m256i t7=_mm256_unpackhi_epi64(row6, row7);
		t0=_mm256_add_epi64(t0, t1);//[b2+b3, a2+a3, b0+b1, a0+a1]
		t2=_mm256_add_epi64(t2, t3);//[d2+d3, c2+c3, d0+d1, c0+c1]
		t4=_mm256_add_epi64(t4, t5);
		t6=_mm256_add_epi64(t6, t7);
		t0=_mm256_permute4x64_epi64(t0, _MM_SHUFFLE(3, 1, 2, 0));//[b2+b3, b0+b1, a2+a3, a0+a1]
		t2=_mm256_permute4x64_epi64(t2, _MM_SHUFFLE(3, 1, 2, 0));//[d2+d3, d0+d1, c2+c3, c0+c1]
		t4=_mm256_permute4x64_epi64(t4, _MM_SHUFFLE(3, 1, 2, 0));
		t6=_mm256_permute4x64_epi64(t6, _MM_SHUFFLE(3, 1, 2, 0));
		row0=_mm256_unpacklo_epi64(t0, t2);//[d0+d1, b0+b1, c2+c3, a0+a1]
		row1=_mm256_unpackhi_epi64(t0, t2);//[d2+d3, b2+b3, c0+c1, a2+a3]
		row2=_mm256_unpacklo_epi64(t4, t6);
		row3=_mm256_unpackhi_epi64(t4, t6);
		row0=_mm256_add_epi64(row0, row1);//[d0+d1+d2+d3, b0+b1+b2+b3, c0+c1+c2+c3, a0+a1+a2+a3] the product
		row2=_mm256_add_epi64(row2, row3);
		ALIGN(32) long long params[8];
		_mm256_store_si256((__m256i*)params+0, row0);
		_mm256_store_si256((__m256i*)params+1, row1);
		for(int k=0;k<8;++k)
			params[k]>>=16;
		//for(int ky=0;ky<8;++ky)//sign-extend int32 to int64
		//{
		//	__m128i v0=_mm_load_si128((__m128i*)sqm+(ky<<2|2));
		//	__m128i v1=_mm_load_si128((__m128i*)sqm+(ky<<2|3));
		//	__m256i v2=_mm256_cvtepi32_epi64(v0);
		//	__m256i v3=_mm256_cvtepi32_epi64(v1);
		//	_mm256_store_si256((__m256i*)sqm+(ky<<1|0), v2);
		//	_mm256_store_si256((__m256i*)sqm+(ky<<1|1), v3);
		//}
#endif
	}
#endif
#endif
#ifdef ENABLE_CUSTOM1
	int ols=0;
	if(kx&&pr->ky&&kx<pr->iw-1)
	{
		ALIGN(32) int nb2[]=
		{
			N, W, NW, NE, eN, eW, eNW, eNE,
		};
		char
			*params_N=pr->c1_params+((size_t)(pr->iw*kc+kx)<<3),
			*params_W=pr->c1_params+((size_t)(pr->iw*kc+kx-1)<<3),
			*params_NE=pr->c1_params+((size_t)(pr->iw*kc+kx+1)<<3);
		
		custom_N =custom1_predict(nb2, params_N);
		custom_W =custom1_predict(nb2, params_W);
		custom_NE=custom1_predict(nb2, params_NE);
		//int pN=custom1_predict(nb2, params_N);
		//int pW=custom1_predict(nb2, params_W);
		//int pNE=custom1_predict(nb2, params_NE);
		//custom1=(3*pW+pN+4*pNE+4)>>3;
#if 0
		ALIGN(16) long long p[6]={0};//{W, 0, N, 0, NE, 0}	*[3, 1, 4]
		memcpy(p+0, params_W , sizeof(*p));
		memcpy(p+2, params_N , sizeof(*p));
		memcpy(p+4, params_NE, sizeof(*p));
		__m256i val=_mm256_load_si256((__m256i*)nb2);
		__m256i pW =_mm256_cvtepi8_epi32(_mm_load_si128((__m128i*)p+0));//W
		__m256i pN =_mm256_cvtepi8_epi32(_mm_load_si128((__m128i*)p+1));//N
		__m256i pNE=_mm256_cvtepi8_epi32(_mm_load_si128((__m128i*)p+2));//NE
		__m256i par=_mm256_slli_epi32(pW, 1);
		par=_mm256_add_epi32(par, pW);
		par=_mm256_add_epi32(par, pN);
		pNE=_mm256_slli_epi32(pNE, 2);
		par=_mm256_add_epi32(par, pNE);
		val=_mm256_mullo_epi32(val, par);
		val=_mm256_hadd_epi32(val, val);
		val=_mm256_hadd_epi32(val, val);
		val=_mm256_hadd_epi32(val, val);
		custom1=_mm256_extract_epi32(val, 0)>>10;//den = 8<<7
#endif
	}
#endif

	int j=-1;
#define SLIC5_PRED(X) pr->preds[++j]=X;
	SLIC5_PREDLIST
#undef  SLIC5_PRED
	PROF(CALC_SUBPREDS);

	long long weights[SLIC5_NPREDS]={0}, wsum=0;
	for(int k=0;k<SLIC5_NPREDS;++k)
	{
		weights[k]=
			(long long)LOAD_PRED_ERROR(kc, -2, 2, k)+	//peNNEE
			(long long)LOAD_PRED_ERROR(kc, -1, 2, k)+	//peNNE
			(long long)LOAD_PRED_ERROR(kc,  0, 2, k)+	//peNN+peNW
			(long long)LOAD_PRED_ERROR(kc,  1, 2, k)+	//peNNW+peNWW
			(long long)LOAD_PRED_ERROR(kc, -2, 1, k)+	//peNEE
			(long long)LOAD_PRED_ERROR(kc, -1, 1, k)+	//peNE
			(long long)LOAD_PRED_ERROR(kc,  0, 1, k)*3+	//peN+peW
			(long long)LOAD_PRED_ERROR(kc,  1, 1, k);	//peNW+peWW
		
		//	(long long)LOAD_PRED_ERROR(kc, -2, 2, k)+	//peNNEE
		//	(long long)LOAD_PRED_ERROR(kc, -1, 2, k)+	//peNNE
		//	(long long)LOAD_PRED_ERROR(kc,  0, 2, k)*2+	//peNN+peNW
		//	(long long)LOAD_PRED_ERROR(kc,  1, 2, k)+	//peNNW+peNWW
		//	(long long)LOAD_PRED_ERROR(kc,  2, 2, k)+	//peNNWW+peNWWW
		//	(long long)LOAD_PRED_ERROR(kc, -2, 1, k)+	//peNEE
		//	(long long)LOAD_PRED_ERROR(kc, -1, 1, k)*3+	//peNE
		//	(long long)LOAD_PRED_ERROR(kc,  0, 1, k)*6+	//peN+peW
		//	(long long)LOAD_PRED_ERROR(kc,  1, 1, k)*2+	//peNW+peWW
		//	(long long)LOAD_PRED_ERROR(kc,  2, 1, k);	//peNWW+peWWW

//#ifdef COMPENSATE_sRGB
		weights[k]=((long long)pr->params[k<<2|kc]<<29)/(weights[k]+1);
//#else
//		weights[k]=((long long)pr->params[k<<2|kc]<<(pr->depths[kc]+PRED_PREC+5))/(weights[k]+1);
//#endif
		wsum+=weights[k];
	}
	if(wsum)
	{
		pr->pred=0;
		for(int k=0;k<SLIC5_NPREDS;++k)
			pr->pred+=weights[k]*pr->preds[k];
		pr->pred+=(wsum>>1)-1;
		pr->pred/=wsum;
	}
	else
		pr->pred=pr->preds[4];
	PROF(CALC_WEIGHT_AV);

#ifdef USE_ABAC
#if 0
	{
		int g3[NPAR_CTX]=
		{
			//X
			(NNE-NN)>>1,			//5	2/64
			(NE-NNE)>>1,			//6	2/64
			eNW,				//9	1/64
			eNE,				//10	1/64
			(3*(NE-NW)+NNWW-NNEE)>>3,	//4	8/512
			(eNE-eNNE)>>1,			//8	2/64
			(eNNE-eNN)>>1,			//7	2/64
			(NW+3*(NE-N)-NEE)>>3,		//1	8/512
			(NWW+3*(N-NW)-NE)>>3,		//2	8/512
			(3*(N-W)+WW-NN)>>3,		//3	8/512

			//Y
			(N-NW)>>1,		//5	2/16
			(N-NN)>>1,		//6	2/16
			eW,			//9	2/1
			eWW,			//10	2/1
			(W-WW)>>1,		//4	2/8
			(eN-eNN)>>1,		//8	2/16
			(eN-eNW)>>1,		//7	2/16
			(NE-NNEE)>>1,		//1	2/8
			(N-NN)>>1,		//2	2/8
			(NW-NNWW)>>1,		//3	2/8

			//Z
			(W-WW)>>1,		//5	2/1
			(W-NW)>>1,		//6	2/1
			eN,			//9	2/1
			eNN,			//10	2/1
			(NW+2*NE+WW-4*W)>>3,	//4	8/4
			(eW-eNW)>>1,		//8	2/1
			(eW-eWW)>>1,		//7	2/1
			(NNE+NEE+2*W-4*NE)>>3,	//1	8/4
			(W+NW+NE+NN-4*N)>>3,	//2	8/4
			(NNW+NWW+W+N-4*NW)>>3,	//3	8/4
		};
		for(int k=0;k<NPAR_CTX;++k)
			pr->qpreds[k]=quantize_signed(g3[k], sh, 2, 1, pr->qplevels);
			//pr->qpreds[k]=quantize_signed(pr->preds[k]-pr->preds[(k+1)%SLIC5_NPREDS], sh, 2, 1, pr->qplevels);
	}
#endif
	int qx=QUANTIZE_HIST((dx+abs(eW+eWW))>>(sh-2)), qy=QUANTIZE_HIST((dy+abs(eN+eNN))>>(sh-2));
	qx=CLAMP(0, qx, pr->xtextures-1);
	qy=CLAMP(0, qy, pr->ytextures-1);
	pr->txid=pr->xtextures*qy+qx;
#else
	int qx=QUANTIZE_HIST((dx+abs(eW+eWW))>>(sh-2)), qy=QUANTIZE_HIST((dy+abs(eN+eNN))>>(sh-2));
	qx=CLAMP(0, qx, pr->nhist-1);
	qy=CLAMP(0, qy, pr->nhist-1);
	pr->hist_idx=pr->nhist*qy+qx;
#endif

#ifdef ENABLE_SSE_PAR
	int g[3][10]=
	{
		{
			//X
			(NNE-NN)>>1,			//5	2/64
			(NE-NNE)>>1,			//6	2/64
			eNW,				//9	1/64
			eNE,				//10	1/64
			(3*(NE-NW)+NNWW-NNEE)>>3,	//4	8/512
			(eNE-eNNE)>>1,			//8	2/64
			(eNNE-eNN)>>1,			//7	2/64
			(NW+3*(NE-N)-NEE)>>3,		//1	8/512
			(NWW+3*(N-NW)-NE)>>3,		//2	8/512
			(3*(N-W)+WW-NN)>>3,		//3	8/512
		},
		{
			//Y
			(N-NW)>>1,		//5	2/16
			(N-NN)>>1,		//6	2/16
			eW,			//9	2/1
			eWW,			//10	2/1
			(W-WW)>>1,		//4	2/8
			(eN-eNN)>>1,		//8	2/16
			(eN-eNW)>>1,		//7	2/16
			(NE-NNEE)>>1,		//1	2/8
			(N-NN)>>1,		//2	2/8
			(NW-NNWW)>>1,		//3	2/8
		},
		{
			//Z
			(W-WW)>>1,		//5	2/1
			(W-NW)>>1,		//6	2/1
			eN,			//9	2/1
			eNN,			//10	2/1
			(NW+2*NE+WW-4*W)>>3,	//4	8/4
			(eW-eNW)>>1,		//8	2/1
			(eW-eWW)>>1,		//7	2/1
			(NNE+NEE+2*W-4*NE)>>3,	//1	8/4
			(W+NW+NE+NN-4*N)>>3,	//2	8/4
			(NNW+NWW+W+N-4*NW)>>3,	//3	8/4
		},
	};
	j=-1;
	for(int ks=0;ks<SSE_STAGES;++ks)
	{
		long long sse_val;
		int qp;

		qp=(int)(pr->pred>>(sh+8-(SSE_PREDBITS-1)));
		qp+=(1<<SSE_PREDBITS)>>1;
		qp=CLAMP(0, qp, (1<<SSE_PREDBITS)-1);

		++j;
		pr->sse_idx[j]=quantize_signed(g[0][ks], sh, SSE_X_EXP, SSE_X_MSB, SSE_W);
		pr->sse_idx[j]=pr->sse_size*ks+(pr->sse_idx[j]<<SSE_PREDBITS|qp);
		sse_val=pr->sse[pr->sse_idx[j]];
		pr->sse_count[j]=(int)(sse_val&0xFFF);
		pr->sse_sum[j]=sse_val>>12;

		++j;
		pr->sse_idx[j]=quantize_signed(g[0][ks], sh, SSE_Y_EXP, SSE_Y_MSB, SSE_H);
		pr->sse_idx[j]=pr->sse_size*ks+((pr->sse_w+pr->sse_idx[j])<<SSE_PREDBITS|qp);
		sse_val=pr->sse[pr->sse_idx[j]];
		pr->sse_count[j]=(int)(sse_val&0xFFF);
		pr->sse_sum[j]=sse_val>>12;

		++j;
		pr->sse_idx[j]=quantize_signed(g[0][ks], sh, SSE_Z_EXP, SSE_Z_MSB, SSE_D);
		pr->sse_idx[j]=pr->sse_size*ks+((pr->sse_w+pr->sse_h+pr->sse_idx[j])<<SSE_PREDBITS|qp);
		sse_val=pr->sse[pr->sse_idx[j]];
		pr->sse_count[j]=(int)(sse_val&0xFFF);
		pr->sse_sum[j]=sse_val>>12;

		//n1/d1 + n2/d2 + n3/d3 = (n1*d2*d3 + d1*n2*d3 + d1*d2*n3)/(d1*d2*d3)
		//sse_val=(
		//	pr->sse_sum[j]*pr->sse_count[j-1]*pr->sse_count[j-2]+
		//	pr->sse_sum[j-1]*pr->sse_count[j]*pr->sse_count[j-2]+
		//	pr->sse_sum[j-2]*pr->sse_count[j]*pr->sse_count[j-1]
		//)/(3LL*pr->sse_count[j]*pr->sse_count[j-1]*pr->sse_count[j-2]+1);
		sse_val=(
			pr->sse_sum[j-2]/(pr->sse_count[j-2]+1LL)+
			pr->sse_sum[j-1]/(pr->sse_count[j-1]+1LL)+
			pr->sse_sum[j-0]/(pr->sse_count[j-0]+1LL)
		)/3;

		pr->pred+=sse_val;
		pr->sse_corr+=(int)sse_val;
	}
	{
		++j;
		int idx=0, temp;
		idx=3*idx+THREEWAY(N, (int)pr->pred)+1;
		idx=3*idx+THREEWAY(W, (int)pr->pred)+1;
		idx=3*idx+THREEWAY(NW, (int)pr->pred)+1;
		idx=3*idx+THREEWAY(NE, (int)pr->pred)+1;
		idx=3*idx+THREEWAY(NN, (int)pr->pred)+1;
		idx=3*idx+THREEWAY(WW, (int)pr->pred)+1;
		temp=N+W-NW, idx=3*idx+THREEWAY(temp, (int)pr->pred)+1;
		temp=(4*(W+N)-(WW+NN)+6)/12, idx=3*idx+THREEWAY(temp, (int)pr->pred)+1;
		temp=2*N-NN, idx=3*idx+THREEWAY(temp, (int)pr->pred)+1;
		temp=2*W-WW, idx=3*idx+THREEWAY(temp, (int)pr->pred)+1;
		if(idx>=SSE_FR_SIZE)
			LOG_ERROR("SSE OOB");
		pr->sse_idx[j]=SSE_FR_SIZE*kc+idx;

		long long sse_val=pr->sse_fr[pr->sse_idx[j]];
		pr->sse_count[j]=(int)(sse_val&0xFFF);
		pr->sse_sum[j]=sse_val>>12;

		sse_val=pr->sse_sum[j]/(pr->sse_count[j]+1LL);

		pr->pred+=sse_val;
		pr->sse_corr+=(int)sse_val;
	}
#endif
#ifdef ENABLE_SSE_4D
#ifdef AVX2
	ALIGN(32) int g[]=
	{
		//X
		(NNE-NN)>>6,			//5	2/64
		(NE-NNE)>>6,			//6	2/64
		eNW>>6,				//9	1/64
		eNE>>6,				//10	1/64
		(3*(NE-NW)+NNWW-NNEE)>>9,	//4	8/512
		(eNE-eNNE)>>6,			//8	2/64
		(eNNE-eNN)>>6,			//7	2/64
		(NW+3*(NE-N)-NEE)>>9,		//1	8/512
		(NWW+3*(N-NW)-NE)>>9,		//2	8/512
		(3*(N-W)+WW-NN)>>9,		//3	8/512
		0,
		0,
		0,
		0,
		0,
		0,

		//Y
		(N-NW)>>4,		//5	2/16
		(N-NN)>>4,		//6	2/16
		eW<<1,			//9	2/1
		eWW<<1,			//10	2/1
		(W-WW)>>3,		//4	2/8
		(eN-eNN)>>4,		//8	2/16
		(eN-eNW)>>4,		//7	2/16
		(NE-NNEE)>>3,		//1	2/8
		(N-NN)>>3,		//2	2/8
		(NW-NNWW)>>3,		//3	2/8
		0,
		0,
		0,
		0,
		0,
		0,

		//Z
		W-WW,			//5	2/1
		W-NW,			//6	2/1
		eN<<1,			//9	2/1
		eNN<<1,			//10	2/1
		(NW+2*NE+WW-4*W)>>2,	//4	8/4
		eW-eNW,			//8	2/1
		eW-eWW,			//7	2/1
		(NNE+NEE+2*W-4*NE)>>2,	//1	8/4
		(W+NW+NE+NN-4*N)>>2,	//2	8/4
		(NNW+NWW+W+N-4*NW)>>2,	//3	8/4
		0,
		0,
		0,
		0,
		0,
		0,
	};
	for(int k=0;k<(kc<<1);k+=2)
	{
		int kc2=k>>1;
		int
			NW	=LOAD_CH(pr->pixels, kc2,  0, 1),
			N	=LOAD_CH(pr->pixels, kc2,  0, 1),
			NE	=LOAD_CH(pr->pixels, kc2, -1, 1),
			W	=LOAD_CH(pr->pixels, kc2,  1, 0),
			curr	=LOAD_CH(pr->pixels, kc2,  0, 0),
			eNW	=LOAD_CH(pr->errors, kc2,  0, 1),
			eN	=LOAD_CH(pr->errors, kc2,  0, 1),
			eNE	=LOAD_CH(pr->errors, kc2, -1, 1),
			eW	=LOAD_CH(pr->errors, kc2,  1, 0),
			ecurr	=LOAD_CH(pr->errors, kc2,  0, 0);
		g[10+k]=(N+W-NW+curr)>>7;	//X 1/32
		g[20+k]=curr<<1;		//Y 2/1
		g[30+k]=N+W;			//Z 2/1

		g[10+k+1]=(W+NE-N+curr)>>7;
		g[20+k+1]=curr+ecurr;
		g[30+k+1]=eN+eW;
	}
	__m128i shift=_mm_set_epi32(0, 0, 0, sh);
	__m256i nlevels[]=
	{
		_mm256_set1_epi32(SSE_W),
		_mm256_set1_epi32(SSE_H),
		_mm256_set1_epi32(SSE_D),
	};
	__m256i half[]=
	{
		_mm256_set1_epi32(SSE_W>>1),
		_mm256_set1_epi32(SSE_H>>1),
		_mm256_set1_epi32(SSE_D>>1),
	};
	for(int k=0;k<_countof(g)-(sizeof(__m256i)/sizeof(int)-1);k+=sizeof(__m256i)/sizeof(int))
	{
		__m256i val=_mm256_load_si256((__m256i*)(g+k));
		__m256i neg=_mm256_srai_epi32(val, 31);
		val=_mm256_sra_epi32(val, shift);
		val=_mm256_abs_epi32(val);
		avx2_floor_log2_p1(&val);
		val=_mm256_xor_si256(val, neg);
		val=_mm256_sub_epi32(val, neg);
		val=_mm256_add_epi32(val, half[k>>4]);//k>>4 == k/(sizeof(__m256i)/sizeof(int))*sizeof(half)/sizeof(g0)
		_mm256_store_si256((__m256i*)(g+k), val);
	}
	__m256i X0=_mm256_load_si256((__m256i*)g+0);
	__m256i X1=_mm256_load_si256((__m256i*)g+1);
	__m256i Y0=_mm256_load_si256((__m256i*)g+2);
	__m256i Y1=_mm256_load_si256((__m256i*)g+3);
	__m256i Z0=_mm256_load_si256((__m256i*)g+4);
	__m256i Z1=_mm256_load_si256((__m256i*)g+5);
	Z0=_mm256_mullo_epi32(Z0, nlevels[1]);
	Z1=_mm256_mullo_epi32(Z1, nlevels[1]);
	Y0=_mm256_add_epi32(Z0, Y0);
	Y1=_mm256_add_epi32(Z1, Y1);
	Y0=_mm256_mullo_epi32(Y0, nlevels[0]);
	Y1=_mm256_mullo_epi32(Y1, nlevels[0]);
	X0=_mm256_add_epi32(Y0, X0);
	X1=_mm256_add_epi32(Y1, X1);
	ALIGN(32) int g2[16];
	_mm256_store_si256((__m256i*)g2+0, X0);
	_mm256_store_si256((__m256i*)g2+1, X1);
#else
#define QUANTIZE(IDX, X, Y, Z)\
	SSE_W*(SSE_H*(pr->sse_idx_z[IDX]=quantize_signed(Z, sh, SSE_Z_EXP, SSE_Z_MSB, SSE_D))+\
	(pr->sse_idx_y[IDX]=quantize_signed(Y, sh, SSE_Y_EXP, SSE_Y_MSB, SSE_H)))+\
	(pr->sse_idx_x[IDX]=quantize_signed(X, sh, SSE_X_EXP, SSE_X_MSB, SSE_W))
	int g2[SSE_STAGES+6]=
	{
		//4D SSE

		QUANTIZE(0, (NNE-NN)>>6,		(N-NW)>>4,		W-WW),//5
		QUANTIZE(1, (NE-NNE)>>6,		(N-NN)>>4,		W-NW),//6
		QUANTIZE(2, eNW>>6,			eW<<1,			eN<<1),//9
		QUANTIZE(3, eNE>>6,			eWW<<1,			eNN<<1),//10
		QUANTIZE(4, (3*(NE-NW)+NNWW-NNEE)>>9,	(W-WW)>>3,		(NW+2*NE+WW-4*W)>>2),//4
		QUANTIZE(5, (eNE-eNNE)>>6,		(eN-eNN)>>4,		eW-eNW),//8
		QUANTIZE(6, (eNNE-eNN)>>6,		(eN-eNW)>>4,		eW-eWW),//7
		QUANTIZE(7, (NW+3*(NE-N)-NEE)>>9,	(NE-NNEE)>>3,		(NNE+NEE+2*W-4*NE)>>2),//1
		QUANTIZE(8, (NWW+3*(N-NW)-NE)>>9,	(N-NN)>>3,		(W+NW+NE+NN-4*N)>>2),//2
		QUANTIZE(9, (3*(N-W)+WW-NN)>>9,		(NW-NNWW)>>3,		(NNW+NWW+W+N-4*NW)>>2),//3
	};
	for(int k=0;k<(kc<<1);k+=2)
	{
		int kc2=k>>1;
		int
			NW	=LOAD_CH(pr->pixels, kc2,  0, 1),
			N	=LOAD_CH(pr->pixels, kc2,  0, 1),
			NE	=LOAD_CH(pr->pixels, kc2, -1, 1),
			W	=LOAD_CH(pr->pixels, kc2,  1, 0),
			curr	=LOAD_CH(pr->pixels, kc2,  0, 0),
			eNW	=LOAD_CH(pr->errors, kc2,  0, 1),
			eN	=LOAD_CH(pr->errors, kc2,  0, 1),
			eNE	=LOAD_CH(pr->errors, kc2, -1, 1),
			eW	=LOAD_CH(pr->errors, kc2,  1, 0),
			ecurr	=LOAD_CH(pr->errors, kc2,  0, 0);

		g2[SSE_STAGES+k  ]=QUANTIZE(SSE_STAGES+k, (N+W-NW+curr)>>7, curr<<1, N+W);
		g2[SSE_STAGES+k+1]=QUANTIZE(SSE_STAGES+k+1, (W+NE-N+curr)>>7, curr+ecurr, eN+eW);
	}
#undef  QUANTIZE
#endif
	PROF(CALC_CTX);

	pr->sse_corr=0;
	for(int k=0;k<SSE_STAGES+1+(kc<<1);++k)
	{
		long long sse_val;
		if(k<SSE_STAGES)
		{
			pr->sse_idx_p[k]=quantize_signed((int)pr->pred, sh+8-(SSE_PREDBITS-1), SSE_P_EXP, SSE_P_MSB, pr->sse_p);//shift is sh+8-(SSE_PREDBITS-1): remove the 8-bits, and insert up to SSE_PREDBITS
			//int qp=(int)(pr->pred>>(sh+8-(SSE_PREDBITS-1)));
			//qp+=(1<<SSE_PREDBITS)>>1;
			//qp=CLAMP(0, qp, (1<<SSE_PREDBITS)-1);
			pr->sse_idx[k]=SSE_SIZE*(SSE_STAGES*kc+k)+(SSE_W*SSE_H*SSE_D)*pr->sse_idx_p[k]+g2[k];
			sse_val=pr->sse[pr->sse_idx[k]];
			
#ifdef TRACK_SSE_RANGES
			int kx=g2[k]%pr->sse_width, ky=g2[k]/pr->sse_width%pr->sse_height, kz=g2[k]/(pr->sse_width*pr->sse_height);
			UPDATE_MIN(pr->sse_ranges[0], kx);
			UPDATE_MAX(pr->sse_ranges[1], kx);
			UPDATE_MIN(pr->sse_ranges[2], ky);
			UPDATE_MAX(pr->sse_ranges[3], ky);
			UPDATE_MIN(pr->sse_ranges[4], kz);
			UPDATE_MAX(pr->sse_ranges[5], kz);
			UPDATE_MIN(pr->sse_ranges[6], qp);
			UPDATE_MAX(pr->sse_ranges[7], qp);
#endif
		}
		else if(k<SSE_STAGES+(kc<<1))
		{
			pr->sse_idx_p[k]=quantize_signed((int)pr->pred, sh+8-(SSE_PREDBITS-1), SSE_P_EXP, SSE_P_MSB, pr->sse_p);
			pr->sse_idx[k]=SSE_SIZE*(k-SSE_STAGES)+(SSE_W*SSE_H*SSE_D)*pr->sse_idx_p[k]+g2[k];
			sse_val=pr->sse_cfl[pr->sse_idx[k]];
		}
		else
		{
			int temp;
			int idx=0;
		//	idx=3*idx+THREEWAY(eN, eW)+1;//X
		//	idx=3*idx+THREEWAY(eN, 0)+1;//X
		//	idx=3*idx+THREEWAY(eW, 0)+1;//X
		//	idx=3*idx+THREEWAY(N, W)+1;//X
			idx=3*idx+THREEWAY(N, (int)pr->pred)+1;
			idx=3*idx+THREEWAY(W, (int)pr->pred)+1;
			idx=3*idx+THREEWAY(NW, (int)pr->pred)+1;
			idx=3*idx+THREEWAY(NE, (int)pr->pred)+1;
			idx=3*idx+THREEWAY(NN, (int)pr->pred)+1;
			idx=3*idx+THREEWAY(WW, (int)pr->pred)+1;
			temp=N+W-NW, idx=3*idx+THREEWAY(temp, (int)pr->pred)+1;
			temp=(4*(W+N)-(WW+NN)+6)/12, idx=3*idx+THREEWAY(temp, (int)pr->pred)+1;
			temp=2*N-NN, idx=3*idx+THREEWAY(temp, (int)pr->pred)+1;
			temp=2*W-WW, idx=3*idx+THREEWAY(temp, (int)pr->pred)+1;
			if(idx>=SSE_FR_SIZE)
				LOG_ERROR("SSE OOB");
			pr->sse_idx[k]=SSE_FR_SIZE*kc+idx;
			//pr->sse_idx[k]=
			//	//(W<NW)<<11|
			//	//(NW<NE)<<10|
			//	(N<W)<<9|
			//	(N+W-NW<(int)pr->pred)<<8|
			//	(W <pr->pred)<<7|
			//	(NW<pr->pred)<<6|
			//	(N <pr->pred)<<5|
			//	(NE<pr->pred)<<4|
			//	(WW<pr->pred)<<3|
			//	(NN<pr->pred)<<2|
			//	(2*N-NN<(int)pr->pred)<<1|
			//	(2*W-WW<(int)pr->pred)
			//;
			sse_val=pr->sse_fr[pr->sse_idx[k]];
		}
		pr->sse_count[k]=(int)(sse_val&0xFFF);
		pr->sse_sum[k]=(int)(sse_val>>12);
		int sse_corr=(int)(pr->sse_sum[k]*SSE_STEP/(pr->sse_count[k]+1LL));
		pr->pred+=sse_corr;
		pr->sse_corr+=sse_corr;
	}
#endif

	int final_corr=(int)(pr->bias_sum[kc]/(pr->bias_count[kc]+1LL));
	pr->pred+=final_corr;
	pr->sse_corr+=final_corr;
	PROF(SSE_LOOP);

	pr->pred=CLAMP(clamp_lo, pr->pred, clamp_hi);
#ifdef COMPENSATE_sRGB
	pr->pred_final=linear2sRGB((int)pr->pred);
	pr->pred_final+=pr->rounding_offset[kc];//for rounding
	pr->pred_final>>=pr->shift_prec[kc];
#else
	pr->pred_final=(int)((pr->pred+pr->rounding_offset[kc])>>pr->shift_prec[kc]);
	//pr->pred_final=(int)((pr->pred+((1<<PRED_PREC)>>1)-1)>>PRED_PREC);
#endif
}
#ifndef USE_ABAC
static void slic5_update_hist(SLIC5Ctx *pr, int token)
{
	int total_hists=pr->nhist*pr->nhist,
		*curr_hist=pr->hist+pr->cdfsize*(total_hists*pr->kc+pr->hist_idx),
		*curr_sum=pr->histsums+total_hists*pr->kc+pr->hist_idx;

	++curr_hist[token];
	++*curr_sum;
	
	if(*curr_sum>pr->cdfsize&&!(*curr_sum&63))
	{
		int sum=*curr_sum;
		CDF_t *curr_CDF=pr->CDFs+(pr->cdfsize+1)*(total_hists*pr->kc+pr->hist_idx);
		long long c=0;
		for(int ks=0;ks<pr->cdfsize;++ks)
		{
			long long freq=curr_hist[ks];
			curr_CDF[ks]=(int)(c*((1LL<<CDF_SHIFT)-pr->cdfsize)/sum)+ks;
			c+=freq;
		}
		curr_CDF[pr->cdfsize]=1<<CDF_SHIFT;
	}

	//if(*curr_sum>(pr->iw*pr->nch<<1))
	if(*curr_sum>=6144)//rescale hist
	{
		int sum=0;
		for(int ks=0;ks<pr->cdfsize;++ks)
		{
			curr_hist[ks]=(curr_hist[ks]+1)>>1;
			sum+=curr_hist[ks];
		}
		*curr_sum=sum;
	}
}
#endif
static void slic5_update(SLIC5Ctx *pr, int curr, int token)
{
	PROF(PRED_TILL_UPDATE);
	int kc=pr->kc, kx=pr->kx;

#ifdef COMPENSATE_sRGB
	curr<<=24-pr->depths[kc];
	curr=sRGB2linear(curr);
#else
	curr<<=pr->shift_prec[kc];
	//curr<<=PRED_PREC;
#endif
	LOAD(pr->pixels, 0, 0)=curr;

	//update WP errors
	int error=curr-(int)pr->pred;
	LOAD(pr->errors, 0, 0)=error;
	int errors[SLIC5_NPREDS]={0}, kbest=0;
	for(int k=0;k<SLIC5_NPREDS;++k)
	{
		errors[k]=abs(curr-pr->preds[k]);
		LOAD_PRED_ERROR(kc, 0, 0, k)=errors[k];
		if(pr->ky&&pr->kx+1<pr->iw)
			LOAD_PRED_ERROR(kc, -1, 1, k)+=errors[k];//eNE += ecurr
		if(errors[kbest]>errors[k])
			kbest=k;

		pr->pred_error_sums[k]+=errors[k];
	}
	PROF(UPDATE_ERRORS);

	//update WP weights
	++pr->params[kbest<<2|kc];
	//if(pr->params[kbest<<2|kc]>(12<<pr->depths[kc]))
	if(pr->params[kbest<<2|kc]>352)
	{
		for(int k=0;k<SLIC5_NPREDS;++k)
			pr->params[k<<2|kc]>>=1;
	}
	PROF(UPDATE_WP);
	
#ifdef ENABLE_CUSTOM1
	ALIGN(32) int nb[]=
	{
		LOAD(pr->pixels,  0, 1),//N
		LOAD(pr->pixels,  1, 0),//W
		LOAD(pr->pixels,  1, 1),//NW
		LOAD(pr->pixels, -1, 1),//NE
		LOAD(pr->errors,  0, 1),
		LOAD(pr->errors,  1, 0),
		LOAD(pr->errors,  1, 1),
		LOAD(pr->errors, -1, 1),
	};
	char *curr_params=pr->c1_params+((size_t)kx<<3);
	
	int vals[3], losses[3];
#define CALC_LOSS(L) L=abs(curr-custom1_predict(nb, curr_params))
#define C1_ONE (1<<7)
	if(kx==256&&pr->ky==256)//
		printf("");
	for(int k=0;k<8;++k)
	{
		const int deviation=C1_ONE/2;
		char *p=curr_params+k;
		vals[0]=CLAMP(-C1_ONE, *p-deviation, C1_ONE-1),
		vals[1]=CLAMP(-C1_ONE, *p+deviation, C1_ONE-1),
		*p=vals[0]; CALC_LOSS(losses[0]);
		*p=vals[1]; CALC_LOSS(losses[1]);
		for(;vals[0]<vals[1];)
		{
			vals[2]=(vals[0]+vals[1])>>1;
			*p=vals[2]; CALC_LOSS(losses[2]);
			if(losses[2]>=losses[0]&&losses[2]>=losses[1])
			{
				if(losses[1]>losses[0])
					*p=vals[0];
				else
					*p=vals[1];
				break;
			}
			if(losses[1]>losses[0])
				losses[1]=losses[2], vals[1]=vals[2];
			else
				losses[0]=losses[2], vals[0]=vals[2];
		}
	}
	if(kx==256&&pr->ky==256)//
		printf("");
#undef  CALC_LOSS
	PROF(UPDATE_CUSTOM1);
#endif

	//update token hist
#ifndef USE_ABAC
	if(token>=0)
		slic5_update_hist(pr, token);
#endif
	PROF(UPDATE_HIST);
	
	//update SSE
#if defined ENABLE_SSE_4D || defined ENABLE_SSE_PAR
	for(int k=0;k<SSE_STAGES+1+(kc<<1);++k)
	{
		long long *sse_ptr;
		if(k<SSE_STAGES)
			sse_ptr=pr->sse;
		else if(k<SSE_STAGES+(kc<<1))
			sse_ptr=pr->sse_cfl;
		else
			sse_ptr=pr->sse_fr;
		pr->sse_count[k]+=SSE_STEP;
		pr->sse_sum[k]+=error;
		if(pr->sse_count[k]>=SSE_LIMIT)
		{
			pr->sse_count[k]>>=1;
			pr->sse_sum[k]>>=1;
		}
		long long sse_val=pr->sse_sum[k]<<12|pr->sse_count[k];
		sse_ptr[pr->sse_idx[k]]=sse_val;
		
#ifdef SSE_UPDATE_NB_CELLS
		int e2=error/SSE_STEP;
		if(e2&&k<SSE_STAGES+(kc<<1))
		{
			int idx2, count;
			long long sum;
			if(pr->sse_idx_y[k]>0)
			{
				idx2=SSE_W*(SSE_H*(SSE_D*pr->sse_idx_p[k]+pr->sse_idx_z[k])+pr->sse_idx_y[k]-1)+pr->sse_idx_x[k];
				if(k<SSE_STAGES)
					idx2+=SSE_SIZE*(SSE_STAGES*kc+k);
				else
					idx2+=SSE_SIZE*(k-SSE_STAGES);
				sse_val=sse_ptr[idx2];
				count=(int)(sse_val&0xFFF);
				sum=sse_val>>12;
				++count;
				sum+=error;
				if(count>=SSE_LIMIT)
				{
					count>>=1;
					sum>>=1;
				}
				sse_ptr[idx2]=sum<<12|count;
			}
			if(pr->sse_idx_y[k]<SSE_H-1)
			{
				idx2=SSE_W*(SSE_H*(SSE_D*pr->sse_idx_p[k]+pr->sse_idx_z[k])+pr->sse_idx_y[k]+1)+pr->sse_idx_x[k];
				if(k<SSE_STAGES)
					idx2+=SSE_SIZE*(SSE_STAGES*kc+k);
				else
					idx2+=SSE_SIZE*(k-SSE_STAGES);
				sse_val=sse_ptr[idx2];
				count=(int)(sse_val&0xFFF);
				sum=sse_val>>12;
				++count;
				sum+=error;
				if(count>=SSE_LIMIT)
				{
					count>>=1;
					sum>>=1;
				}
				sse_ptr[idx2]=sum<<12|count;
			}
			if(pr->sse_idx_z[k]>0)
			{
				idx2=SSE_W*(SSE_H*(SSE_D*pr->sse_idx_p[k]+pr->sse_idx_z[k]-1)+pr->sse_idx_y[k])+pr->sse_idx_x[k];
				if(k<SSE_STAGES)
					idx2+=SSE_SIZE*(SSE_STAGES*kc+k);
				else
					idx2+=SSE_SIZE*(k-SSE_STAGES);
				sse_val=sse_ptr[idx2];
				count=(int)(sse_val&0xFFF);
				sum=sse_val>>12;
				++count;
				sum+=error;
				if(count>=SSE_LIMIT)
				{
					count>>=1;
					sum>>=1;
				}
				sse_ptr[idx2]=sum<<12|count;
			}
			if(pr->sse_idx_z[k]<SSE_D-1)
			{
				idx2=SSE_W*(SSE_H*(SSE_D*pr->sse_idx_p[k]+pr->sse_idx_z[k]+1)+pr->sse_idx_y[k])+pr->sse_idx_x[k];
				if(k<SSE_STAGES)
					idx2+=SSE_SIZE*(SSE_STAGES*kc+k);
				else
					idx2+=SSE_SIZE*(k-SSE_STAGES);
				sse_val=sse_ptr[idx2];
				count=(int)(sse_val&0xFFF);
				sum=sse_val>>12;
				++count;
				sum+=error;
				if(count>=SSE_LIMIT)
				{
					count>>=1;
					sum>>=1;
				}
				sse_ptr[idx2]=sum<<12|count;
			}
		}
#endif
	}
#endif
	++pr->bias_count[kc];
	pr->bias_sum[kc]+=error;
	PROF(UPDATE_SSE);
}
#ifdef ENABLE_LZ
static void slic5_skip_lzpixels(SLIC5Ctx *pr, int *lzpixels, int lzlen)
{
#if 0
	int kx0=pr->kx;
	for(int kx=0;kx<lzlen;++kx)
	{
		for(int kc=0;kc<pr->nch;++kc)
		{
			slic5_predict(pr, kc, kx0+kx, pr->ky);
			slic5_update(pr, lzpixels[kx<<2|kc], -1);
		}
	}
#endif
#if 1
	size_t idx=(size_t)pr->kym0+pr->kx+PAD_SIZE;
	memcpy(pr->pixels+(idx<<2), lzpixels, lzlen*sizeof(int[4]));
	memset(pr->errors+(idx<<2), 0, lzlen*sizeof(int[4]));
	for(int k=0;k<SLIC5_NPREDS;++k)
		memset(pr->pred_errors+((SLIC5_NPREDS*idx+k)<<2), 0, lzlen*sizeof(int[4]));
#endif
}
#endif
#undef  LOAD
#ifdef USE_ABAC
static int slic5_abac_predict(StatNode const *node, int *p0_ctr)
{
	int wsum=node->weight[0]+node->weight[1]+node->weight[2]+node->weight[3]+node->weight[4]+node->weight[5]+node->weight[6]+7;
	*p0_ctr=((node->n[0]+1)<<16)/(node->n[0]+node->n[1]+2);
	long long p0=
		node->rec[0]*(node->weight[0]+1LL)+
		node->rec[1]*(node->weight[1]+1LL)+
		node->rec[2]*(node->weight[2]+1LL)+
		node->rec[3]*(node->weight[3]+1LL)+
		node->rec[4]*(node->weight[4]+1LL)+
		node->rec[5]*(node->weight[5]+1LL)+
		*p0_ctr*(node->weight[6]+1LL);
	p0/=wsum;
#ifdef USE_ABAC_SSE
	p0+=node->sse_sum/(node->sse_ctr+1LL);
#endif
	p0=CLAMP(1, p0, 0xFFFF);
	return (int)p0;
}
static void slic5_abac_update(StatNode *node, int p0_ctr, int p0, int bit)
{
	//small vs large change:  shift right ammount is "resistance", more slows it down, less speeds it up
	static const int shift_p0[]=
	{
		3, 4, 5, 6, 7, 8,
	};
	static const int shift_weight[]=
	{
		5, 5, 9, 10, 11, 12, 12,
	};
	int prob_error=(!bit<<16)-p0;
	node->weight[0]+=(((!bit==(node->rec[0]	>0x8000))<<16)-node->weight[0]+((1<<shift_weight[0])>>1))>>shift_weight[0];
	node->weight[1]+=(((!bit==(node->rec[1]	>0x8000))<<16)-node->weight[1]+((1<<shift_weight[1])>>1))>>shift_weight[1];
	node->weight[2]+=(((!bit==(node->rec[2]	>0x8000))<<16)-node->weight[2]+((1<<shift_weight[2])>>1))>>shift_weight[2];
	node->weight[3]+=(((!bit==(node->rec[3]	>0x8000))<<16)-node->weight[3]+((1<<shift_weight[3])>>1))>>shift_weight[3];
	node->weight[4]+=(((!bit==(node->rec[4]	>0x8000))<<16)-node->weight[4]+((1<<shift_weight[4])>>1))>>shift_weight[4];
	node->weight[5]+=(((!bit==(node->rec[5]	>0x8000))<<16)-node->weight[5]+((1<<shift_weight[5])>>1))>>shift_weight[5];
	node->weight[6]+=(((!bit==(p0_ctr	>0x8000))<<16)-node->weight[6]+((1<<shift_weight[6])>>1))>>shift_weight[6];
	
	node->rec[0]+=((!bit<<16)-node->rec[0]+((1<<shift_p0[0])>>1))>>shift_p0[0];
	node->rec[1]+=((!bit<<16)-node->rec[1]+((1<<shift_p0[1])>>1))>>shift_p0[1];
	node->rec[2]+=((!bit<<16)-node->rec[2]+((1<<shift_p0[2])>>1))>>shift_p0[2];
	node->rec[3]+=((!bit<<16)-node->rec[3]+((1<<shift_p0[3])>>1))>>shift_p0[3];
	node->rec[4]+=((!bit<<16)-node->rec[4]+((1<<shift_p0[4])>>1))>>shift_p0[4];
	node->rec[5]+=((!bit<<16)-node->rec[5]+((1<<shift_p0[5])>>1))>>shift_p0[5];
	if(node->n[bit]+1>0xFFF)
	{
		node->n[0]>>=1;
		node->n[1]>>=1;
	}
	++node->n[bit];
	
#ifdef USE_ABAC_SSE
	if(node->sse_ctr+1>0xFFFF)
	{
		node->sse_ctr>>=1;
		node->sse_sum>>=1;
	}
	++node->sse_ctr;
	node->sse_sum+=prob_error;
#endif
}
#endif
static void slic5_enc(SLIC5Ctx *pr, int curr, int kc, int kx)
{
	slic5_predict(pr, kc, kx);

	int error=curr-pr->pred_final;
#if 1
	int upred=pr->half[kc]-abs(pr->pred_final);
	int abs_error=abs(error), negmask;
	if(abs_error<=upred)
	{
		negmask=-((pr->sse_corr<0)&(error!=-pr->half[kc]));//sign is flipped if SSE correction was negative, to skew the histogram
		error^=negmask;
		error-=negmask;
		error=error<<1^-(error<0);//pack sign
	}
	else
		error=upred+abs_error;//error sign is known
#endif
#if 0
	int upred=pr->pred_final+(pr->nlevels[kc]>>1), abs_error=abs(error), negmask;
	if(upred<(pr->nlevels[kc]>>1))
	{
		if(abs_error<=upred)
		{
			negmask=-(pr->sse_corr<0);//sign is flipped if SSE correction was negative, to skew the histogram
			error^=negmask;
			error-=negmask;
			error=error<<1^-(error<0);//pack sign
		}
		else
			error=upred+abs(error);
	}
	else
	{
		if(abs_error<=pr->nlevels[kc]-upred)
		{
			negmask=-(pr->sse_corr<0);
			error^=negmask;
			error-=negmask;
			error=error<<1^-(error<0);
		}
		else
			error=pr->nlevels[kc]-upred+abs(error);
	}
#endif
	if(error<0)
		LOG_ERROR("Negative symbol %d", error);

#ifdef USE_ABAC
	StatNode *curr_stats=pr->stats+pr->treesize*(pr->ntextures*kc+pr->txid);
	for(int kb=0, MSBidx=0;kb<pr->depths[kc];++kb)
	{
		StatNode *cell=curr_stats+(kb*(kb+1)>>1)+kb-MSBidx;
		int p0_ctr=0, p0=slic5_abac_predict(cell, &p0_ctr);

		int bit=error>>(pr->depths[kc]-1-kb)&1;
		ac_enc_bin(pr->ec, (unsigned short)p0, bit);

		slic5_abac_update(cell, p0_ctr, p0, bit);
		if(MSBidx==kb)
			MSBidx+=!bit;
	}
	slic5_update(pr, curr, 0);
#if 0
	if(hu.token>pr->maxtoken
#ifdef ENABLE_LZ
		-1//LZ escape symbol is OOB here
#endif
	)
		LOG_ERROR("Token OOB %d/%d", hu.token, pr->maxtoken+1);
	StatNode *curr_stats[NPAR_CTX];
	for(int k=0;k<NPAR_CTX;++k)
		curr_stats[k]=pr->stats+(((size_t)pr->qplevels*(NPAR_CTX*kc+k)+pr->qpreds[k])<<pr->treedepth);
		//curr_stats[k]=pr->stats+(((size_t)pr->qplevels*(SLIC5_NPREDS*kc+k)+pr->qpreds[k])<<pr->treedepth);
	//StatNode *curr_stats=pr->stats+(((size_t)pr->ntextures*kc+pr->txid)<<pr->treedepth);
	int abac_idx=1;
	for(int kb=pr->treedepth-1;kb>=0;--kb)
	{
		int wsum=0;
		long long p0_final=0;
		int p0[NPAR_CTX], p0_ctr[NPAR_CTX];
		for(int k=0;k<NPAR_CTX;++k)
		{
			p0[k]=slic5_abac_predict(curr_stats[k]+abac_idx, p0_ctr+k);
			wsum+=pr->ctx_weights[k<<2|kc];
			p0_final+=(long long)pr->ctx_weights[k<<2|kc]*p0[k];
		}
		p0_final/=wsum;
		p0_final=CLAMP(1, p0_final, 0xFFFF);
		//StatNode *node=curr_stats+abac_idx;
		//int p0_ctr=0, p0=slic5_abac_predict(node, &p0_ctr);
		int bit=hu.token>>kb&1;

		ac_enc_bin(pr->ec, (unsigned short)p0_final, bit);
		
		for(int k=0;k<NPAR_CTX;++k)
			slic5_abac_update(curr_stats[k]+abac_idx, p0_ctr[k], p0[k], bit);
		//slic5_abac_update(node, p0_ctr, p0, bit);
		abac_idx<<=1;
		abac_idx|=bit;
	}
#endif
#else
	HybridUint hu;
	hybriduint_encode(error, &hu);
	if(hu.token>=pr->cdfsize
#ifdef ENABLE_LZ
		-1//LZ escape symbol is OOB here
#endif
	)
		LOG_ERROR("Token OOB %d/%d", hu.token, pr->cdfsize);

	EC_ENC(pr->ec, hu.token, pr->CDFs+(pr->cdfsize+1)*(pr->nhist*pr->nhist*kc+pr->hist_idx));
	if(hu.nbits)
	{
		int bypass=hu.bypass, nbits=hu.nbits;
		while(nbits>8)
		{
			ac_enc_bypass(pr->ec, bypass>>(nbits-8)&0xFF, 1<<8);
			nbits-=8;
		}
		ac_enc_bypass(pr->ec, bypass&((1<<nbits)-1), 1<<nbits);
	}

	slic5_update(pr, curr, hu.token);
#endif
}
static int slic5_dec(SLIC5Ctx *pr, int kc, int kx)
{
	slic5_predict(pr, kc, kx);
	
#ifdef USE_ABAC
	int error=0;
	StatNode *curr_stats=pr->stats+pr->treesize*(pr->ntextures*kc+pr->txid);
	for(int kb=0, MSBidx=0;kb<pr->depths[kc];++kb)
	{
		StatNode *cell=curr_stats+(kb*(kb+1)>>1)+kb-MSBidx;
		int p0_ctr=0, p0=slic5_abac_predict(cell, &p0_ctr);

		int bit=ac_dec_bin(pr->ec, (unsigned short)p0);
		error|=bit<<(pr->depths[kc]-1-kb);

		slic5_abac_update(cell, p0_ctr, p0, bit);
		if(MSBidx==kb)
			MSBidx+=!bit;
	}
	int upred=pr->half[kc]-abs(pr->pred_final), negmask;
	if(error<=(upred<<1))
	{
		error=error>>1^-(error&1);
		negmask=-((pr->sse_corr<0)&(error!=-pr->half[kc]));
	}
	else
	{
		error=error-upred;
		negmask=-(pr->pred_final>0);
	}
	error^=negmask;
	error-=negmask;
	int curr=error+pr->pred_final;
	slic5_update(pr, curr, 0);
	return curr;
#if 0
	int token=0;
	StatNode *curr_stats[NPAR_CTX];
	for(int k=0;k<NPAR_CTX;++k)
		curr_stats[k]=pr->stats+(((size_t)pr->qplevels*(NPAR_CTX*kc+k)+pr->qpreds[k])<<pr->treedepth);
	int abac_idx=1;
	for(int kb=pr->treedepth-1;kb>=0;--kb)
	{
		int wsum=0;
		long long p0_final=0;
		int p0[NPAR_CTX], p0_ctr[NPAR_CTX];
		for(int k=0;k<NPAR_CTX;++k)
		{
			p0[k]=slic5_abac_predict(curr_stats[k]+abac_idx, p0_ctr+k);
			wsum+=pr->ctx_weights[k<<2|kc];
			p0_final+=(long long)pr->ctx_weights[k<<2|kc]*p0[k];
		}
		p0_final/=wsum;
		p0_final=CLAMP(1, p0_final, 0xFFFF);

		int bit=ac_dec_bin(pr->ec, (unsigned short)p0_final);
		token|=bit<<kb;
		
		for(int k=0;k<NPAR_CTX;++k)
			slic5_abac_update(curr_stats[k]+abac_idx, p0_ctr[k], p0[k], bit);
		abac_idx<<=1;
		abac_idx|=bit;
	}
#endif
#else
	int token=EC_DEC(pr->ec, pr->CDFs+(pr->cdfsize+1)*(pr->nhist*pr->nhist*kc+pr->hist_idx), pr->cdfsize);
#ifdef ENABLE_LZ
	if(token==pr->cdfsize-1)
		return 0x80000000;
#endif
	int error=token;
	if(error>=(1<<SLIC5_CONFIG_EXP))
	{
		error-=1<<SLIC5_CONFIG_EXP;
		int lsb=error&((1<<SLIC5_CONFIG_LSB)-1);
		error>>=SLIC5_CONFIG_LSB;
		int msb=error&((1<<SLIC5_CONFIG_MSB)-1);
		error>>=SLIC5_CONFIG_MSB;
		int nbits=error+SLIC5_CONFIG_EXP-(SLIC5_CONFIG_MSB+SLIC5_CONFIG_LSB), n=nbits;
		int bypass=0;
		while(n>8)
		{
			n-=8;
			bypass|=ac_dec_bypass(pr->ec, 1<<8)<<n;
		}
		bypass|=ac_dec_bypass(pr->ec, 1<<n);
		error=1;
		error<<=SLIC5_CONFIG_MSB;
		error|=msb;
		error<<=nbits;
		error|=bypass;
		error<<=SLIC5_CONFIG_LSB;
		error|=lsb;
	}
#if 1
	int upred=pr->half[kc]-abs(pr->pred_final), negmask;
	if(error<=(upred<<1))
	{
		error=error>>1^-(error&1);
		negmask=-((pr->sse_corr<0)&(error!=-pr->half[kc]));
	}
	else
	{
		error=error-upred;
		negmask=-(pr->pred_final>0);
	}
	error^=negmask;
	error-=negmask;
#endif
#if 0
	int upred=pr->pred_final+(pr->nlevels[kc]>>1), negmask;
	if(upred<(pr->nlevels[kc]>>1))
	{
		if(error<=(upred<<1))
		{
			negmask=-(pr->sse_corr<0);
			error=error>>1^-(error&1);
			error^=negmask;
			error-=negmask;
		}
		else
			error=error-upred;
	}
	else
	{
		upred=pr->nlevels[kc]-upred;
		if(error<=(upred<<1))
		{
			negmask=-(pr->sse_corr<0);
			error=error>>1^-(error&1);
			error^=negmask;
			error-=negmask;
		}
		else
			error=upred-error;//here error is negative if prediction is positive
	}
#endif
	int curr=error+pr->pred_final;
	slic5_update(pr, curr, token);
	return curr;
#endif
}

#ifdef SLIC5_OPTIMIZE_RCT
#define ORCT_NITER 50
//#define ORCT_NITER 4096
#define ORCT_DELTAGROUP 1
#define ORCT_WATCHDOGTIMEOUT (ORCT_NPARAMS*16)

#define OCRT_PARAMBITS 7
#define ORCT_ONE (1<<OCRT_PARAMBITS)
static void orct_unpack_permutation(char p, char *permutation)
{
	int temp=0;
	switch(p)
	{
	case 0:temp=0x020100;break;//rgb
	case 1:temp=0x020001;break;//grb
	case 2:temp=0x000201;break;//gbr
	case 3:temp=0x010200;break;//rbg
	case 4:temp=0x000102;break;//bgr
	case 5:temp=0x010002;break;//brg
	default:
		LOG_ERROR("Invalid RGB permutation");
		return;
	}
	memcpy(permutation, &temp, sizeof(char[3]));
}
static void orct_fwd(int *comp, const char *params, const char *permutation)
{
	int temp;
	//const int half=(1<<OCRT_PARAMBITS)>>1;
	int c2[3]=
	{
		comp[permutation[0]],
		comp[permutation[1]],
		comp[permutation[2]],
	};
	temp=params[0]*c2[1]+params[1]*c2[2], c2[0]+=(temp>>OCRT_PARAMBITS)+(temp<0);
	temp=params[2]*c2[0]+params[3]*c2[2], c2[1]+=(temp>>OCRT_PARAMBITS)+(temp<0);
	temp=params[4]*c2[0]+params[5]*c2[1], c2[2]+=(temp>>OCRT_PARAMBITS)+(temp<0);
	temp=params[6]*c2[0]+params[7]*c2[2], c2[1]+=(temp>>OCRT_PARAMBITS)+(temp<0);
	//c2[0]+=(params[0]*c2[1]+params[1]*c2[2]+half)>>OCRT_PARAMBITS;
	//c2[1]+=(params[2]*c2[0]+params[3]*c2[2]+half)>>OCRT_PARAMBITS;
	//c2[2]+=(params[4]*c2[0]+params[5]*c2[1]+half)>>OCRT_PARAMBITS;
	//c2[1]+=(params[6]*c2[0]+params[7]*c2[2]+half)>>OCRT_PARAMBITS;
	memcpy(comp, c2, sizeof(c2));
}
static void orct_inv(int *comp, const char *params, const char *permutation)
{
	int temp;
	//const int half=(1<<OCRT_PARAMBITS)>>1;
	int c2[3];
	memcpy(c2, comp, sizeof(c2));
	temp=params[6]*c2[0]+params[7]*c2[2], c2[1]-=(temp>>OCRT_PARAMBITS)+(temp<0);
	temp=params[4]*c2[0]+params[5]*c2[1], c2[2]-=(temp>>OCRT_PARAMBITS)+(temp<0);
	temp=params[2]*c2[0]+params[3]*c2[2], c2[1]-=(temp>>OCRT_PARAMBITS)+(temp<0);
	temp=params[0]*c2[1]+params[1]*c2[2], c2[0]-=(temp>>OCRT_PARAMBITS)+(temp<0);
	//c2[1]-=(params[6]*c2[0]+params[7]*c2[2]+half)>>OCRT_PARAMBITS;
	//c2[2]-=(params[4]*c2[0]+params[5]*c2[1]+half)>>OCRT_PARAMBITS;
	//c2[1]-=(params[2]*c2[0]+params[3]*c2[2]+half)>>OCRT_PARAMBITS;
	//c2[0]-=(params[0]*c2[1]+params[1]*c2[2]+half)>>OCRT_PARAMBITS;
	comp[permutation[0]]=c2[0];
	comp[permutation[1]]=c2[1];
	comp[permutation[2]]=c2[2];
}
static void orct_calc_depths(const char *params, const char *permutation, const char *srcdepths, char *dstdepths)
{
	int extremevals[]=
	{
		//-((1<<srcdepths[0])>>1), -((1<<srcdepths[0])>>2), ((1<<srcdepths[0])>>2)-1, ((1<<srcdepths[0])>>1)-1,
		//-((1<<srcdepths[1])>>1), -((1<<srcdepths[1])>>2), ((1<<srcdepths[1])>>2)-1, ((1<<srcdepths[1])>>1)-1,
		//-((1<<srcdepths[2])>>1), -((1<<srcdepths[2])>>2), ((1<<srcdepths[2])>>2)-1, ((1<<srcdepths[2])>>1)-1,
		-((1<<srcdepths[0])>>1), ((1<<srcdepths[0])>>1)-1,
		-((1<<srcdepths[1])>>1), ((1<<srcdepths[1])>>1)-1,
		-((1<<srcdepths[2])>>1), ((1<<srcdepths[2])>>1)-1,
	};
	int ranges[6]={0};
	//for(int k=0;k<64;++k)
	//for(int k=0;k<27;++k)//zero is useless here
	for(int k=0;k<8;++k)
	{
		//int comp[]={extremevals[k&3], extremevals[4|k>>2&3], extremevals[8|k>>4&3]};
		//int comp[]={extremevals[k%3], extremevals[3+k/3%3], extremevals[6+k/9%3]};
		int comp[]={extremevals[k&1], extremevals[2|k>>1&1], extremevals[4|k>>2&1]};
		orct_fwd(comp, params, permutation);
		UPDATE_MIN(ranges[0], comp[0]);
		UPDATE_MAX(ranges[1], comp[0]);
		UPDATE_MIN(ranges[2], comp[1]);
		UPDATE_MAX(ranges[3], comp[1]);
		UPDATE_MIN(ranges[4], comp[2]);
		UPDATE_MAX(ranges[5], comp[2]);
	}
	ranges[0]=abs(ranges[0]);
	ranges[2]=abs(ranges[2]);
	ranges[4]=abs(ranges[4]);
	++ranges[1];
	++ranges[3];
	++ranges[5];
	UPDATE_MAX(ranges[1], ranges[0]);
	UPDATE_MAX(ranges[3], ranges[2]);
	UPDATE_MAX(ranges[5], ranges[4]);
	dstdepths[0]=ceil_log2(ranges[3])+1;
	dstdepths[1]=ceil_log2(ranges[5])+1;
	dstdepths[2]=ceil_log2(ranges[1])+1;
#if 0
	int gains[3];
	int maxdepth=MAXVAR(srcdepths[0], srcdepths[1]);
	maxdepth=MAXVAR(maxdepth, srcdepths[2]);

	gains[0]=ORCT_ONE+abs(params[0])+abs(params[1]);
	gains[1]=
		abs(params[6]+((params[7]*params[4]+(ORCT_ONE+((params[5]*params[7])>>OCRT_PARAMBITS))*params[2])>>OCRT_PARAMBITS))+
		abs(params[0]*(params[6]+((params[4]*params[7])>>OCRT_PARAMBITS))+(1+((params[5]*params[7])>>OCRT_PARAMBITS))*(1+((params[0]*params[2])>>OCRT_PARAMBITS)))+
		abs(((params[1]*(params[6]+((params[4]*params[7])>>OCRT_PARAMBITS))+(1+((params[5]*params[7])>>OCRT_PARAMBITS))*(((params[1]*params[2])>>OCRT_PARAMBITS)+params[3]))>>OCRT_PARAMBITS)+params[7]);
	gains[2]=
		abs(params[5]+((params[3]*params[6])>>OCRT_PARAMBITS))+
		abs((params[1]*params[5]+params[6]*(ORCT_ONE+((params[1]*params[3])>>OCRT_PARAMBITS)))>>OCRT_PARAMBITS)+
		abs(ORCT_ONE+((params[2]*params[5]+params[6]*(((params[2]*params[3])>>OCRT_PARAMBITS)+params[4]))>>OCRT_PARAMBITS));

	//gains[0]=ORCT_ONE+params[0]+params[1];//Cr	X
	//gains[2]=ORCT_ONE+((gains[0]*params[4]+(ORCT_ONE+params[3]+(gains[0]*params[2]>>OCRT_PARAMBITS))*params[5])>>OCRT_PARAMBITS);//Cb
	//gains[1]=ORCT_ONE+params[3]+((gains[0]*(params[2]+params[6])+gains[2]*params[7])>>OCRT_PARAMBITS);//Y
	dstdepths[0]=(char)ceil_log2((long long)gains[1]<<maxdepth>>OCRT_PARAMBITS);
	dstdepths[1]=(char)ceil_log2((long long)gains[2]<<maxdepth>>OCRT_PARAMBITS);
	dstdepths[2]=(char)ceil_log2((long long)gains[0]<<maxdepth>>OCRT_PARAMBITS);
#endif
}
static double orct_calcloss(Image const *src, int *pixels, int *hist, int maxlevels, const char *params)
{
	//int nlevels[]=
	//{
	//	1<<src->depth[1],
	//	1<<(src->depth[2]+1),
	//	1<<(src->depth[0]+1),
	//};
	char permutation[3]={0}, depths[3];
	orct_unpack_permutation(params[8], permutation);
	orct_calc_depths(params, permutation, src->depth, depths);
	int nlevels[]=
	{
		1<<depths[0],
		1<<depths[1],
		1<<depths[2],
	};
	memset(hist, 0, maxlevels*sizeof(int[3]));
	memset(pixels, 0, (src->iw+2LL)*sizeof(int[2*4]));
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int row0=ky&1, row1=!row0;
		for(int kx=0;kx<src->iw;++kx, ++idx)
		{
			int N, W, NW, pred;
			int comp[]=
			{
				src->data[idx<<2|0],
				src->data[idx<<2|1],
				src->data[idx<<2|2],
			};
			orct_fwd(comp, params, permutation);

#define GET_PIXEL(C, X, Y) pixels[(src->iw*row##Y+kx+1-X)<<2|C]
#define PROCESS_SUBPIXEL(C, PIXEL)\
	GET_PIXEL(C, 0, 0)=PIXEL,\
	N=GET_PIXEL(C, 0, 1),\
	W=GET_PIXEL(C, 1, 0),\
	NW=GET_PIXEL(C, 1, 1),\
	pred=N+W-NW, pred=MEDIAN3(N, W, pred), PIXEL-=pred,\
	++hist[maxlevels*C+((PIXEL+(nlevels[C]>>1))&(nlevels[C]-1))]

			PROCESS_SUBPIXEL(0, comp[1]);
			PROCESS_SUBPIXEL(1, comp[2]);
			PROCESS_SUBPIXEL(2, comp[0]);
#undef  GET_PIXEL
#undef  PROCESS_SUBPIXEL
		}
	}
	int res=src->iw*src->ih;
	int kbest=0;
	double bestsize=0;
	double csizes[3]={0};
	for(int kc=0;kc<3;++kc)
	{
		int *curr_hist=hist+maxlevels*kc;
		for(int ks=0;ks<maxlevels;++ks)
		{
			//if((kc==0||kc==1)&&ks==128)
			//	printf("");
			int freq=curr_hist[ks];
			if(freq)
			{
				double p=(double)freq/res;
				double bitsize=-log2(p);
				csizes[kc]+=freq*bitsize;
			}
		}
	}
	double csize=(csizes[0]+csizes[1]+csizes[2])/8;
	return csize;
}
//typedef struct ORCTInfoStruct
//{
//	double loss;
//	char params[9];
//} ORCTInfo;
static const char *slic5_orct_permutationnames[]=
{
	"rgb",//0
	"grb",//1
	"gbr",//2
	"rbg",//3
	"bgr",//4
	"brg",//5
};
void orct_print_compact(const char *params)
{
	printf("[%s", slic5_orct_permutationnames[params[ORCT_NPARAMS]]);
	for(int k=0;k<ORCT_NPARAMS;++k)
	{
		int val=params[k];
		printf("%c%02X", val<0?'-':'+', abs(val));
	}
	printf("]");
}
static void orct_optimize(Image const *src, char *params, int loud)
{
	static const char paramsets[]=//too lazy to initialize the permutations correctly
	{
		0,	0,//0	RCT_NONE defeats the purpose of optimization
		0,	0,
		0,	0,
		0,	0,
		0,//rgb
		
		-ORCT_ONE,	0,//1	RCT_SubG
		0,		0,
		0,		-ORCT_ONE,
		0,		0,
		0,//rgb

		-ORCT_ONE,	0,//2	RCT_JPEG2000
		0,		0,
		0,		-ORCT_ONE,
		ORCT_ONE>>1,	ORCT_ONE>>1,
		0,//rgb

		-ORCT_ONE,	0,//3	RCT_YCbCr-R_v1
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		0,		ORCT_ONE>>1,
		0,//rgb

		-ORCT_ONE,	0,//4	RCT_A710
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		-(ORCT_ONE>>3),	ORCT_ONE>>4,
		0,//rgb

		-ORCT_ONE,	0,//5	RCT_YCbCr-R_v3
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		ORCT_ONE>>3,	ORCT_ONE>>4,
		0,//rgb

		-ORCT_ONE,	0,//6	RCT_YCbCr-R_v4 "fair luma"
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		0,		(ORCT_ONE+1)/3,
		0,//rgb

		-ORCT_ONE,	0,//7	RCT_YCbCr-R_v5
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		0,		ORCT_ONE*3>>3,
		0,//rgb

		-ORCT_ONE,	0,//8	RCT_YCbCr-R_v6
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		0,		ORCT_ONE*7>>4,
		0,//rgb

		-ORCT_ONE,	0,//9	RCT_YCbCr-R_v7
		ORCT_ONE>>1,	0,
		0,		-ORCT_ONE,
		-(ORCT_ONE>>5),	ORCT_ONE*10>>5,
		0,//rgb

		-(ORCT_ONE*87>>8),	-(ORCT_ONE*169>>8),	//10	RCT_Pei09
		0,			0,
		0,			-ORCT_ONE,
		ORCT_ONE*29>>8,		ORCT_ONE*86>>8,
		0,//rgb
	};
	static const char *rctnames[]=
	{
		"NONE",
		"SubGreen",
		"JPEG2000",
		"YCbCr_R_v1",
		"A710",
		"YCbCr_R_v3",
		"YCbCr_R_v4",
		"YCbCr_R_v5",
		"YCbCr_R_v6",
		"YCbCr_R_v7",
		"Pei09",
	};
	double t_start=time_sec();
	int inflation[]={2, 2, 2, 0};
	int maxdepth=calc_maxdepth(src, inflation), maxlevels=1<<maxdepth;
	int *pixels=(int*)malloc((src->iw+2LL)*sizeof(int[2*4]));//2 rows * 4 channels
	int *hist=(int*)malloc(maxlevels*sizeof(int[4]));//4 channels
	if(!pixels||!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	//ORCTInfo info;
	//memcpy(info.params, params, sizeof(info.params));
#define CALC_LOSS(L) L=orct_calcloss(src, pixels, hist, maxlevels, params2)
	double loss_init, loss_bestsofar, loss_prev, loss_curr;
#if 0
	char params2[ORCT_NPARAMS+1]={0};
	CALC_LOSS(loss_curr);
	loss_init=loss_bestsofar=loss_curr;
	memcpy(params, params2, sizeof(params2));
	for(int kb=7;kb>=0;--kb)
	{
		for(int kc=0;kc<ORCT_NPARAMS;++kc)
		{
			//for(int kc2=0;kc<ORCT_NPARAMS;++kc)
			//{
				for(int inc=1<<kb, end=-inc;inc>=-end;inc-=2<<kb)
				{
					//for(int inc2=1<<kb, end2=-inc;inc2>=-end2;inc2-=2<<kb)
					//{
						params2[kc]+=inc;
						//params2[kc2]+=inc2;
						CALC_LOSS(loss_curr);
						if(loss_bestsofar>loss_curr)
						{
							loss_bestsofar=loss_curr;
							memcpy(params, params2, sizeof(params2));
						}
						else
						{
							params2[kc]-=inc;
							//params2[kc2]-=inc2;
						}
					//}
				}
			//}
		}
		int bestp=params2[ORCT_NPARAMS];
		for(int k=0;k<6;++k)
		{
			if(k==bestp)
				continue;
			params2[ORCT_NPARAMS]=k;
			CALC_LOSS(loss_prev);
			if(loss_bestsofar>loss_prev)
				loss_bestsofar=loss_prev, bestp=k;
		}
		params2[ORCT_NPARAMS]=bestp;
		if(loud)
			printf("%4d/%4d  best %16lf  curr %16lf  elapsed %12lf sec\r",
				(7-kb)+1, 8, loss_bestsofar, loss_curr, time_sec()-t_start
			);
	}
#endif
#if 1
	char params2[]=
	{
		-ORCT_ONE,	0,		//-1	0		JPEG2000 RCT
		0,		0,		//0	0
		0,		-ORCT_ONE,	//0	-1
		ORCT_ONE>>1,	ORCT_ONE>>1,	//1/2	1/2
		0,//rgb
	};

	int best_init=0, bestp=0;
	loss_bestsofar=INFINITY;
	for(int k=0;k<(_countof(paramsets)/(ORCT_NPARAMS+1));++k)//brute force initial RCT
	{
		memcpy(params2, paramsets+(ORCT_NPARAMS+1)*k, sizeof(params2));
		for(int k2=1;k2<6;++k2)//brute force the permutations
		{
			params2[ORCT_NPARAMS]=k2;
			CALC_LOSS(loss_prev);
			if(loss_bestsofar>loss_prev)
				loss_bestsofar=loss_prev, best_init=k, bestp=k2;
		}
		//params2[ORCT_NPARAMS]=bestp;
		//CALC_LOSS(loss_prev);
		//if(loss_bestsofar>loss_prev)
		//	loss_bestsofar=loss_prev, best_init=k;
	}
	memcpy(params2, paramsets+(ORCT_NPARAMS+1)*best_init, sizeof(params2));
	params2[ORCT_NPARAMS]=bestp;
	if(loud)
		printf("Init %d %s  %lf\n", best_init, rctnames[best_init], loss_bestsofar);

	//for(int k=0;k<ORCT_NPARAMS+1;++k)//
	//	params2[k]=rand();
	//params2[ORCT_NPARAMS]%=6;

	//CALC_LOSS(loss_prev);
	//loss_init=loss_bestsofar=loss_prev;
	loss_init=loss_bestsofar;

	//for(int k=1;k<6;++k)//brute force the permutations
	//{
	//	params2[ORCT_NPARAMS]=k;
	//	CALC_LOSS(loss_prev);
	//	if(loss_bestsofar>loss_prev)
	//		loss_bestsofar=loss_prev, bestp=k;
	//}
	//params2[ORCT_NPARAMS]=bestp;
	memcpy(params, params2, sizeof(params2));//save
	if(loud)
		printf("Permutation %d rgb->%s  %lf\n", bestp, slic5_orct_permutationnames[bestp], loss_bestsofar);

	for(int it=0, checkpoint=-1;it<ORCT_NITER;++it)
	{
		char params3[ORCT_NPARAMS+1];
		memcpy(params3, params2, sizeof(params3));

		//slider
		const int deviation2=ORCT_ONE>>4;
		int vals[3]=
		{
			CLAMP(-ORCT_ONE, params2[0]-deviation2, ORCT_ONE-1),
			CLAMP(-ORCT_ONE, params2[0]+deviation2, ORCT_ONE-1),
		};
		//int vals[3]={-(ORCT_ONE-1), 0};
		double losses[3];
		params2[0]=vals[0], params2[1]=-ORCT_ONE-params2[0]; CALC_LOSS(losses[0]);
		params2[0]=vals[1], params2[1]=-ORCT_ONE-params2[0]; CALC_LOSS(losses[1]);
		for(;vals[0]<vals[1]-1;)
		{
			vals[2]=(vals[0]+vals[1])>>1;
			params2[0]=vals[2], params2[1]=-ORCT_ONE-params2[0]; CALC_LOSS(losses[2]);
			if(losses[2]>=losses[0]&&losses[2]>=losses[1])
			{
				if(losses[1]>losses[0])
					loss_curr=losses[0], params2[0]=vals[0];
				else
					loss_curr=losses[1], params2[0]=vals[1];
				params2[1]=-ORCT_ONE-params2[0];
				break;
			}
			loss_curr=losses[2];
			if(losses[1]>losses[0])
				losses[1]=losses[2], vals[1]=vals[2];
			else
				losses[0]=losses[2], vals[0]=vals[2];
		}

		//luma coefficient 1 and 2 (can take any value)
		const int deviation=ORCT_ONE>>4;
		for(int k=2;k<ORCT_NPARAMS;++k)
		{
			char *p=params2+k;
			switch(k>>1)
			{
			case 2://curr+[-1, 1/2]
				vals[0]=*p-deviation;
				vals[1]=*p+(deviation>>1);
				break;
			case 1://curr+[-1/2, 0.99]
			case 3:
				vals[0]=*p-(deviation>>1);
				vals[1]=*p+deviation;
				break;
			}
#if 0
			if((k>>1)==2)
			{
				vals[0]=*p-ORCT_ONE;
				vals[1]=*p+(ORCT_ONE>>1);
				//continue;
			}
				//vals[0]=-(ORCT_ONE-1), vals[1]=(ORCT_ONE-1)>>1;//WORSE
				//continue;
			else if((k>>1)==3)
			{
				vals[0]=*p-(ORCT_ONE>>1);
				vals[1]=*p+ORCT_ONE-1;
			}
				//vals[0]=-((ORCT_ONE-1)>>1), vals[1]=ORCT_ONE-1;
#endif
			vals[0]=CLAMP(-ORCT_ONE, vals[0], ORCT_ONE-1);
			vals[1]=CLAMP(-ORCT_ONE, vals[1], ORCT_ONE-1);
			*p=vals[0]; CALC_LOSS(losses[0]);
			*p=vals[1]; CALC_LOSS(losses[1]);
			for(;vals[0]<vals[1];)
			{
				vals[2]=(vals[0]+vals[1])>>1;
				*p=vals[2]; CALC_LOSS(losses[2]);
				if(losses[2]>=losses[0]&&losses[2]>=losses[1])
				{
					if(losses[1]>losses[0])
						loss_curr=losses[0], *p=vals[0];
					else
						loss_curr=losses[1], *p=vals[1];
					break;
				}
				loss_curr=losses[2];
				if(losses[1]>losses[0])
					losses[1]=losses[2], vals[1]=vals[2];
				else
					losses[0]=losses[2], vals[0]=vals[2];
			}
		}

		//brute force the permutations
		bestp=params2[ORCT_NPARAMS];
		for(int k=0;k<6;++k)
		{
			if(k==bestp)
				continue;
			params2[ORCT_NPARAMS]=k;
			CALC_LOSS(loss_prev);
			if(loss_bestsofar>loss_prev)
				loss_bestsofar=loss_prev, bestp=k;
		}
		params2[ORCT_NPARAMS]=bestp;
		
		if(loss_bestsofar>loss_curr)
		{
			memcpy(params, params2, sizeof(params2));
			loss_bestsofar=loss_curr;
			checkpoint=it;
		}
		if(loud)
			printf("chk %4d  %4d/%4d  curr %16lf  best %16lf  elapsed %12lf sec\r",
				checkpoint+1, it+1, ORCT_NITER, loss_curr, loss_bestsofar, time_sec()-t_start
			);
		if(!memcmp(params2, params3, sizeof(params2))||it-checkpoint>=5)
			break;
	}
#endif
#if 0
	for(int it=0, watchdog=0, nupdates=0;it<ORCT_NITER;++it)
	{
		int idx[ORCT_DELTAGROUP]={0}, stuck=0;
		int params_original_selected[ORCT_DELTAGROUP];
		if(watchdog>=ORCT_WATCHDOGTIMEOUT)//bump if stuck
		{
			if((rand()&7)<1)//try something new (12.5% chance)
				memcpy(params2, paramsets+(ORCT_NPARAMS+1)*(rand()%(_countof(paramsets)/(ORCT_NPARAMS+1))), sizeof(params2)-1);//exclude the permutation
			else
				memcpy(params2, params, sizeof(params2));//reload
			for(int k=0;k<ORCT_NPARAMS;++k)
				params2[k]+=(rand()&7)-4;
				//params2[k]+=((rand()&1)<<1)-1;
			watchdog=0;
			stuck=1;
		}
		else
		{
			for(int k=0;k<ORCT_DELTAGROUP;++k)
			{
				int inc=0;
				idx[k]=rand()%ORCT_NPARAMS;
				while(!(inc=rand()-(RAND_MAX>>1)));//reject zero delta		TODO xoroshiro128+
		
				params_original_selected[k]=params2[idx[k]];
				params2[idx[k]]+=(inc<<2)/RAND_MAX;

				if((rand()&7)<7)//87.5% chance
				{
					if(idx[k]==0)//high chance a1+a2 should be -1
						params[1]=-64-params[0];
					else if(idx[k]==1)
						params[0]=-64-params[1];
				}
			}
		}
		CALC_LOSS(loss_curr);
		
		if(loss_prev<=loss_curr)//revert if worse
		{
			if(stuck)//a bad branch may surpass the local minimum
				loss_prev=loss_curr;
			else
			{
				loss_curr=loss_prev;
				for(int k=0;k<ORCT_DELTAGROUP;++k)
					params2[idx[k]]=params_original_selected[k];
			}
			++watchdog;
		}
		else//save if better
		{
			if(loss_curr<loss_bestsofar)//publish if record best
			{
				memcpy(params, params2, sizeof(params2));
				loss_bestsofar=loss_curr;
				--it;//bis

				//revise permutations
				int p0=params[ORCT_NPARAMS];
				bestp=params[ORCT_NPARAMS];
				for(int k=0;k<6;++k)
				{
					if(k==p0)//skip original permutation
						continue;
					params2[ORCT_NPARAMS]=k;
					CALC_LOSS(loss_curr);
					if(loss_bestsofar>loss_curr)
						loss_bestsofar=loss_curr, bestp=k;
				}
				params2[ORCT_NPARAMS]=bestp;
				if(loud)
					printf("\n");
			}
			loss_prev=loss_curr;
			watchdog=0;
			++nupdates;
		}
		if(loud)
			printf("%4d/%4d  %4d/%4d  %4d updates  prev %16lf  best %16lf  elapsed %12lf sec\r",
				it+1, ORCT_NITER, watchdog, ORCT_WATCHDOGTIMEOUT, nupdates, loss_prev, loss_bestsofar, time_sec()-t_start
			);
	}
#endif
#undef  CALC_LOSS
	if(loud)
	{
		printf("\n");
		printf("CR improvement %lf -> %lf  %+6.2lf%%\n", loss_init, loss_bestsofar, 100.*(1-loss_init/loss_bestsofar));
		printf("Final RCT params:\n");
		//for(int k=0;k<ORCT_NPARAMS;++k)
		//	printf("  %12g,%c", (double)params[k]/ORCT_ONE, k&1?'\n':' ');
		//printf("  p%d  rgb->%s\n", params[ORCT_NPARAMS], slic5_orct_permutationnames[params[ORCT_NPARAMS]]);
		const char chnames[]="rgb";
		char p[3];
		orct_unpack_permutation(params[ORCT_NPARAMS], p);
		printf("  %c += %12g*%c + %12g*%c\n", chnames[p[0]], (double)params[0]/ORCT_ONE, chnames[p[1]], (double)params[1]/ORCT_ONE, chnames[p[2]]);
		printf("  %c += %12g*%c + %12g*%c\n", chnames[p[1]], (double)params[2]/ORCT_ONE, chnames[p[0]], (double)params[3]/ORCT_ONE, chnames[p[2]]);
		printf("  %c += %12g*%c + %12g*%c\n", chnames[p[2]], (double)params[4]/ORCT_ONE, chnames[p[0]], (double)params[5]/ORCT_ONE, chnames[p[1]]);
		printf("  %c += %12g*%c + %12g*%c\n", chnames[p[1]], (double)params[6]/ORCT_ONE, chnames[p[0]], (double)params[7]/ORCT_ONE, chnames[p[2]]);

		printf("Equivalent RCT matrix:\n");
		double m1[9]={0};
		m1[0+p[0]]=1;
		m1[3+p[1]]=1;
		m1[6+p[2]]=1;
		double m2[9]=
		{
			1, (double)params[0]/ORCT_ONE, (double)params[1]/ORCT_ONE,
			0, 1, 0,
			0, 0, 1,
		};
		double m3[9]=
		{
			1, 0, 0,
			(double)params[2]/ORCT_ONE, 1, (double)params[3]/ORCT_ONE,
			0, 0, 1,
		};
		double m4[9]=
		{
			1, 0, 0,
			0, 1, 0,
			(double)params[4]/ORCT_ONE, (double)params[5]/ORCT_ONE, 1,
		};
		double m5[9]=
		{
			1, 0, 0,
			(double)params[6]/ORCT_ONE, 1, (double)params[7]/ORCT_ONE,
			0, 0, 1,
		};
		double mt1[9], mt2[9];
		matmul(mt1, m2, m1, 3, 3, 3);
		matmul(mt2, m3, mt1, 3, 3, 3);
		matmul(mt1, m4, mt2, 3, 3, 3);
		matmul(mt2, m5, mt1, 3, 3, 3);
		for(int k=0;k<9;++k)
		{
			if(!(k%3))
				printf("%.3s", "Cr Y  Cb "+k);
			printf("%15g, ", mt2[k]);
			if(!((k+1)%3))
				printf("%c\n", "RGB"[k/3]);
			//printf(" %15g,%c", mt2[k], (k+1)%3?' ':'\n');
		}

		printf("Raw data: ");
		orct_print_compact(params);
		//for(int k=0;k<ORCT_NPARAMS+1;++k)
		//	printf(" 0x%02X,", params[k]&0xFF);
		printf("\n");
	}

	free(pixels);
	free(hist);
}
#endif

const char *rct_names[]=
{
#define RCT(X) #X,
	RCTLIST
#undef  RCT
};
#ifdef COMPENSATE_sRGB
#define N2L(SRC) sRGB2linear(p[SRC]<<(24-depths[SRC]))
#define L2N(VAL, DST) linear2sRGB(VAL)>>(24-depths[DST])
#endif
#define r p[0]
#define g p[1]
#define b p[2]
static void rct_fwd(int *p, RCTType rct)
{
	switch(rct)//rounding is to prevent luma from getting OOB
	{
	//case RCT_NONE:
	//		[1	0	0]
	//		[0	0	1]
	//		[0	1	0]
	//	break;
	case RCT_SubGreen:
		r-=g;	//[1	-1	0]
		b-=g;	//[0	-1	1]
		break;	//[0	1	0]
	case RCT_JPEG2000:
		r-=g;		//r-g				[1     -1     0  ].RGB
		b-=g;		//b-g				[0     -1     1  ].RGB
		g+=(r+b+2)>>2;	//g+(r-g+b-g)/4 = r/4+g/2+b/4	[1/4    1/2   1/4].RGB
		break;
	case RCT_YCoCg_R:
		r-=b;		//co = r-b			diff(r, b)		[1	0	-1]
		b+=r>>1;	//(r+b)/2
		g-=b;		//cg = g-(r+b)/2		diff(g, av(r, b))	[-1/2	1	-1/2]
		b+=g>>1;	//Y  = (r+b)/2 + (g-(r+b)/2)/2 = r/4+g/2+b/4	av(g, av(r, b))	[1/4	1/2	1/4]
		{
			int temp;
			SWAPVAR(g, b, temp);//move luma from C2 to C1
		}
		break;
	case RCT_YCbCr_R_v1:
		r-=g;		//diff(r, g)            [ 1      -1      0  ].RGB	Cr
		g+=r>>1;
		b-=g;		//diff(b, av(r, g))     [-1/2    -1/2    1  ].RGB	Cb
		g+=b>>1;	//av(b, av(r, g))       [ 1/4     1/4    1/2].RGB	Y
		break;
	case RCT_A710:
		r-=g;		//Cr =	[1	-1	0].RGB
		g+=r>>1;	//	[1/2	1/2	0]
		b-=g;		//Cb =	[-1/2	-1/2	1]
		g+=(2*b-r+4)>>3;//Y  =	[1/4	1/2	1/4]	v2
		break;
	//case RCT_A710://from GFWX
	//	r-=g;			//Cr = [ 1    -1     0]
	//	b-=(2*g+r)>>1;		//Cb = [-1/2  -1/2   1]
	//	g+=(2*b+3*r+4)>>3;	//Y  = [1/4    1/2   1/4]
	//	break;
	case RCT_YCbCr_R_v3:
		r-=g;		//Cr =	[1	-1	0].RGB
		g+=r>>1;	//	[1/2	1/2	0]
		b-=g;		//Cb =	[-1/2	-1/2	1]
		g+=(2*b+r+4)>>3;//Y  =	[1/2	1/4	1/4]	v3
		break;
	case RCT_YCbCr_R_v4:
		r-=g;		//Cr =	[1	-1	0].RGB
		g+=r>>1;	//	[1/2	1/2	0]
		b-=g;		//Cb =	[-1/2	-1/2	1]
		g+=b/3;		//Y  =	[1/3	1/3	1/3]	v4
		break;
	case RCT_YCbCr_R_v5:
		r-=g;		//Cr =	[1	-1	0].RGB
		g+=r>>1;	//	[1/2	1/2	0]
		b-=g;		//Cb =	[-1/2	-1/2	1]
		g+=(3*b+4)>>3;	//Y  =	[5/16	5/16	6/16]	v5
		break;
	case RCT_YCbCr_R_v6:
		r-=g;		//Cr =	[1	-1	0].RGB
		g+=r>>1;	//	[1/2	1/2	0]
		b-=g;		//Cb =	[-1/2	-1/2	1]
		g+=(7*b+8)>>4;	//Y  =	[9/32	9/32	14/32]	v6
		break;
	case RCT_YCbCr_R_v7:
		r-=g;			//Cr =	[1	-1	0].RGB
		g+=r>>1;		//	[1/2	1/2	0]
		b-=g;			//Cb =	[-1/2	-1/2	1]
		g+=(10*b-r+16)>>5;	//Y  =	[5/16	 6/16	5/16]	v7
		break;
	case RCT_Pei09:
		b-=(87*r+169*g+128)>>8;	//Cb = [-87/256  -169/256  1]
		r-=g;			//Cr = [1  -1  0].RGB
		g+=(86*r+29*b+128)>>8;	//Y  = [19493/65536  38619/65536  29/256]	g+86/256*(r-g)+29/256*(b-87/256*r-169/256*g) = 19493/65536*r + 38619/65536*g + 29/256*b
		break;
	default:
		break;
	}
}
static void rct_inv(int *p, RCTType rct)
{
	switch(rct)
	{
	case RCT_SubGreen:
		b+=g;
		r+=g;
		break;
	case RCT_JPEG2000:
		g-=(r+b+2)>>2;
		b+=g;
		r+=g;
		break;
	case RCT_YCoCg_R:
		{
			int temp;
			SWAPVAR(g, b, temp);//move luma from C2 to C1
		}
		b-=g>>1;
		g+=b;
		b-=r>>1;
		r+=b;
		break;
	case RCT_YCbCr_R_v1:
		g-=b>>1;
		b+=g;
		g-=r>>1;
		r+=g;
		break;
	case RCT_A710:
		g-=(2*b-r+4)>>3;
		b+=g;
		g-=r>>1;
		r+=g;
		break;
	//case RCT_A710:
	//	g-=(2*b+3*r+4)>>3;
	//	b+=(2*g+r)>>1;
	//	r+=g;
	//	break;
	case RCT_YCbCr_R_v3:
		g-=(2*b+r+4)>>3;
		b+=g;
		g-=r>>1;
		r+=g;
		break;
	case RCT_YCbCr_R_v4:
		g-=b/3;
		b+=g;
		g-=r>>1;
		r+=g;
		break;
	case RCT_YCbCr_R_v5:
		g-=(3*b+4)>>3;
		b+=g;
		g-=r>>1;
		r+=g;
		break;
	case RCT_YCbCr_R_v6:
		g-=(7*b+8)>>4;
		b+=g;
		g-=r>>1;
		r+=g;
		break;
	case RCT_YCbCr_R_v7:
		g-=(10*b-r+16)>>5;
		b+=g;
		g-=r>>1;
		r+=g;
		break;
	case RCT_Pei09:
		g-=(86*r+29*b+128)>>8;
		r+=g;
		b+=(87*r+169*g+128)>>8;
		break;
	default:
		break;
	}
}
#undef	r
#undef	g
#undef	b
#ifdef COMPENSATE_sRGB
#undef  N2L
#undef  L2N
#endif
RCTType rct_select_best(Image const *src, double *ret_csizes)
{
	int inflation[]={1, 1, 1, 0};
	int maxdepth=calc_maxdepth(src, inflation), maxlevels=1<<maxdepth;
	int *hist=(int*)malloc((size_t)src->nch*maxlevels*sizeof(int[RCT_COUNT]));
	int *pixels=(int*)malloc((src->iw+2LL)*sizeof(int[4*2*RCT_COUNT]));//(iw + 2 padding) * 4 channels * 2 rows * N_RCT	need 2 rows because of NW
	if(!hist||!pixels)
	{
		LOG_ERROR("Alloc error");
		return RCT_NONE;
	}
	memset(hist, 0, (size_t)src->nch*maxlevels*sizeof(int[RCT_COUNT]));
	memset(pixels, 0, src->iw*sizeof(int[4*2*RCT_COUNT]));
	int res=src->iw*src->ih;
	int nlevels[][3]=
	{
		{1<<src->depth[1], 1<<src->depth[2], 1<<src->depth[0]},//RCT_NONE
		{1<<src->depth[1], 1<<(src->depth[2]+1), 1<<(src->depth[0]+1)},
	};
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int row0=ky&1, row1=!row0;
		for(int kx=0;kx<src->iw;++kx, ++idx)
		{
			int v[]=
			{
				src->data[idx<<2|0],
				src->data[idx<<2|1],
				src->data[idx<<2|2],
				//src->data[idx<<2|0]+(nlevels[0][0]>>1),
				//src->data[idx<<2|1]+(nlevels[0][1]>>1),
				//src->data[idx<<2|2]+(nlevels[0][2]>>1),
			};
			int v2[3], N, W, NW, pred;

#define GET_PIXEL(RCT_IDX, C, X, Y) pixels[(src->iw*(RCT_IDX<<1|row##Y)+kx+1-X)<<2|C]
#define PROCESS_SUBPIXEL(RCT_IDX, C, PIXEL, N_IDX)\
	GET_PIXEL(RCT_IDX, C, 0, 0)=PIXEL,\
	N=GET_PIXEL(RCT_IDX, C, 0, 1),\
	W=GET_PIXEL(RCT_IDX, C, 1, 0),\
	NW=GET_PIXEL(RCT_IDX, C, 1, 1),\
	pred=N+W-NW, pred=MEDIAN3(N, W, pred), PIXEL-=pred,\
	++hist[maxlevels*(3*RCT_IDX+C)+((PIXEL+(nlevels[N_IDX][C]>>1))&(nlevels[N_IDX][C]-1))]

			//RCT_NONE
			memcpy(v2, v, sizeof(v2));
			PROCESS_SUBPIXEL(RCT_NONE, 0, v2[1], 0);
			PROCESS_SUBPIXEL(RCT_NONE, 1, v2[2], 0);
			PROCESS_SUBPIXEL(RCT_NONE, 2, v2[0], 0);

			//RCT_SubGreen
			memcpy(v2, v, sizeof(v2));
			rct_fwd(v2, RCT_SubGreen);
			PROCESS_SUBPIXEL(RCT_SubGreen, 0, v2[1], 1);
			PROCESS_SUBPIXEL(RCT_SubGreen, 1, v2[2], 1);
			PROCESS_SUBPIXEL(RCT_SubGreen, 2, v2[0], 1);

			//RCT_JPEG2000
			memcpy(v2, v, sizeof(v2));
			rct_fwd(v2, RCT_JPEG2000);
			PROCESS_SUBPIXEL(RCT_JPEG2000, 0, v2[1], 1);
			PROCESS_SUBPIXEL(RCT_JPEG2000, 1, v2[2], 1);
			PROCESS_SUBPIXEL(RCT_JPEG2000, 2, v2[0], 1);

			//RCT_YCoCg_R
			memcpy(v2, v, sizeof(v2));
			rct_fwd(v2, RCT_YCoCg_R);
			PROCESS_SUBPIXEL(RCT_YCoCg_R, 0, v2[1], 1);
			PROCESS_SUBPIXEL(RCT_YCoCg_R, 1, v2[2], 1);
			PROCESS_SUBPIXEL(RCT_YCoCg_R, 2, v2[0], 1);

			//RCT_YCbCr_R_v1
			memcpy(v2, v, sizeof(v2));
			rct_fwd(v2, RCT_YCbCr_R_v1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v1, 0, v2[1], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v1, 1, v2[2], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v1, 2, v2[0], 1);

			//RCT_A710
			memcpy(v2, v, sizeof(v2));
			rct_fwd(v2, RCT_A710);
			PROCESS_SUBPIXEL(RCT_A710, 0, v2[1], 1);
			PROCESS_SUBPIXEL(RCT_A710, 1, v2[2], 1);
			PROCESS_SUBPIXEL(RCT_A710, 2, v2[0], 1);

			//RCT_YCbCr_R_v3
			memcpy(v2, v, sizeof(v2));
			rct_fwd(v2, RCT_YCbCr_R_v3);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v3, 0, v2[1], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v3, 1, v2[2], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v3, 2, v2[0], 1);

			//RCT_YCbCr_R_v4
			memcpy(v2, v, sizeof(v2));
			rct_fwd(v2, RCT_YCbCr_R_v4);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v4, 0, v2[1], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v4, 1, v2[2], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v4, 2, v2[0], 1);

			//RCT_YCbCr_R_v5
			memcpy(v2, v, sizeof(v2));
			rct_fwd(v2, RCT_YCbCr_R_v5);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v5, 0, v2[1], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v5, 1, v2[2], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v5, 2, v2[0], 1);

			//RCT_YCbCr_R_v6
			memcpy(v2, v, sizeof(v2));
			rct_fwd(v2, RCT_YCbCr_R_v6);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v6, 0, v2[1], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v6, 1, v2[2], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v6, 2, v2[0], 1);

			//RCT_YCbCr_R_v7
			memcpy(v2, v, sizeof(v2));
			rct_fwd(v2, RCT_YCbCr_R_v7);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v7, 0, v2[1], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v7, 1, v2[2], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v7, 2, v2[0], 1);

			//RCT_Pei09
			memcpy(v2, v, sizeof(v2));
			rct_fwd(v2, RCT_Pei09);
			PROCESS_SUBPIXEL(RCT_Pei09, 0, v2[1], 1);
			PROCESS_SUBPIXEL(RCT_Pei09, 1, v2[2], 1);
			PROCESS_SUBPIXEL(RCT_Pei09, 2, v2[0], 1);
#undef  GET_PIXEL
#undef  PROCESS_SUBPIXEL
		}
	}
	int kbest=0;
	double bestsize=0;
	for(int kt=0;kt<RCT_COUNT;++kt)
	{
		double csizes[3]={0};
		for(int kc=0;kc<3;++kc)
		{
			int *curr_hist=hist+maxlevels*(3*kt+kc);
			for(int ks=0;ks<maxlevels;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
				{
					double p=(double)freq/res;
					double bitsize=-log2(p);
					//if(isinf(bitsize))
					//	LOG_ERROR("ZPS");
					csizes[kc]+=freq*bitsize;
				}
			}
		}
		double csize=csizes[0]+csizes[1]+csizes[2];
		if(!kt||bestsize>csize)
			kbest=kt, bestsize=csize;
		if(ret_csizes)
			ret_csizes[kt]=csize/8;
	}
	free(hist);
	free(pixels);
	return kbest;
}

#define LZ_MAP_SIZE 0x100000//must be a power of two
#define LZ_MIN_EMIT 8
#define HASH_CANTOR(A, B) (((A)+(B)+1)*((A)+(B))/2+(B))
//#define LZ_CELLS 16
//typedef struct LZInfoStruct
//{
//	int x0, y0, pixels_saved;
//} LZInfo;

#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
int t47_encode(Image const *src, ArrayHandle *data, SLIC5Curiosity *curiosity, int loud)
{
	PROF_START();
#ifdef ENABLE_GUIDE
	guide=src;
#endif
	double t_start=time_sec();
	int nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);//TODO
	if(nch>src->nch)
		nch=src->nch;
	if(loud)
	{
		int maxdepth=calc_maxdepth(src, 0);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T47 SLIC5-AC  Enc %s  CWHD %d*%d*%d*%d/8\n", g_buf, nch, src->iw, src->ih, maxdepth);
	}
#ifdef SLIC5_OPTIMIZE_RCT
	char rct_params[ORCT_NPARAMS+1]=
	{
		0

		//-ORCT_ONE,	0,//2	RCT_JPEG2000
		//0,		0,
		//0,		-ORCT_ONE,
		//ORCT_ONE>>1,	ORCT_ONE>>1,
		//0,//rgb

		//0x00, 0xC0, 0x0D, 0xF8, 0x02, 0x02, 0x0E, 0x00, 0x00,
	};
	char permutation[3]={0};
	char depths[4]={0};
	memcpy(depths, src->depth, nch*sizeof(char));
	if(nch>=3)
	{
		orct_optimize(src, rct_params, loud);
		orct_unpack_permutation(rct_params[ORCT_NPARAMS], permutation);
		orct_calc_depths(rct_params, permutation, src->depth, depths);
		PROF(RCT_OPT);
	}
	if(curiosity)
		memcpy(curiosity->rct_params, rct_params, sizeof(curiosity->rct_params));
#else
	RCTType rct=RCT_NONE;
	double csizes[RCT_COUNT]={0};
	if(nch>=3)
	{
		double t1=time_sec();
		if(loud)
			printf("Selecting best RCT...\n");
		rct=rct_select_best(src, csizes);
		if(loud)
		{
			for(int kt=0;kt<(int)RCT_COUNT;++kt)
			{
				printf("%20lf %s", csizes[kt], rct_names[kt]);
				if(kt==(int)rct)
					printf(" <- %lf sec", time_sec()-t1);
				printf("\n");
			}
		}
	}
	if(curiosity)
		curiosity->rct=rct;
#endif
	DList list;
	ArithmeticCoder ec;
	SLIC5Ctx pr;
	dlist_init(&list, 1, 1024, 0);
	ac_enc_init(&ec, &list);
#ifdef SLIC5_OPTIMIZE_RCT
	ptrdiff_t memusage=slic5_init(&pr, src->iw, src->ih, nch, depths, &ec);
#else
	ptrdiff_t memusage=slic5_init(&pr, src->iw, src->ih, nch, src->depth, &ec);
#endif
#ifdef ENABLE_LZ
	int *lzmap=(int*)malloc(sizeof(int[LZ_MAP_SIZE]));
#endif
	if(!memusage
#ifdef ENABLE_LZ
		||!lzmap
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
#ifdef ENABLE_LZ
	memset(lzmap, -1, sizeof(int[LZ_MAP_SIZE]));
#endif
	if(nch>=3)//emit best RCT
#ifdef SLIC5_OPTIMIZE_RCT
		dlist_push_back(&list, rct_params, sizeof(rct_params));
#else
		dlist_push_back1(&list, &rct);
#endif
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		slic5_nextrow(&pr, ky);
		for(int kx=0;kx<src->iw;)
		{
			//if(kx==(src->iw>>1)&&ky==(src->ih>>1))//
			//if(kx==648&&ky==354)//
			//if(kx==5&&ky==7)//
			//if(kx==721&&ky==412)//
			//if(kx==252&&ky==1)//
			//if(kx==465&&ky==0)//
			//if(kx==588&&ky==254)//
			//	printf("");

#ifdef ENABLE_LZ
			int hash=-1, lzidx=-1, lzlen=0;
			if(kx<src->iw-LZ_MIN_EMIT)
			{
				const int *srcptr, *dstptr;
				hash=0;
				for(int kx2=0;kx2<4;++kx2)
				{
					for(int kc=0;kc<nch;++kc)
						hash=HASH_CANTOR(hash, src->data[idx+(kx2<<2)+kc]);
				}
				hash&=LZ_MAP_SIZE-1;
				lzidx=lzmap[hash];
				if(lzidx>=0)
				{
					srcptr=src->data+(((size_t)lzidx)<<2);
					dstptr=src->data+idx;
					for(lzlen=0;kx+lzlen<src->iw&&!memcmp(srcptr+((size_t)lzlen<<2), dstptr+((size_t)lzlen<<2), nch*sizeof(int));++lzlen);
				}
				else if(ky)
				{
					srcptr=src->data+idx-((size_t)src->iw<<2);
					dstptr=src->data+idx;
					for(lzlen=0;kx+lzlen<src->iw&&!memcmp(srcptr+((size_t)lzlen<<2), dstptr+((size_t)lzlen<<2), nch*sizeof(int));++lzlen);
				}
			}
			if(lzlen>4)
			{
				int kc=0;
				//int comp[]={src->data[idx|0], src->data[idx|1], src->data[idx|2]};
				//if(nch>=3)
				//	rct_fwd(comp, rct);
				slic5_predict(&pr, 0, kx, ky);
				EC_ENC(&ec, pr.cdfsize-1, pr.CDFs+(pr.cdfsize+1)*(pr.nhist*pr.nhist*0+pr.hist_idx), pr.cdfsize, 0);//escape symbol
				EC_ENC(&ec, lzidx>>24&0xFF, 0, 256, 0);
				if(lzidx>=0)
				{
					EC_ENC(&ec, lzidx>>16&0xFF, 0, 256, 0);
					EC_ENC(&ec, lzidx>>8&0xFF, 0, 256, 0);
					EC_ENC(&ec, lzidx&0xFF, 0, 256, 0);
				}
				EC_ENC(&ec, lzlen&0xFF, 0, 256, 0);
				EC_ENC(&ec, lzlen>>8&0xFF, 0, 256, 0);
				slic5_update_hist(&pr, pr.cdfsize-1);
				slic5_update_CDFs(&pr);
				//slic5_update(&pr, nch>=3?comp[1]:comp[0], pr.cdfsize-1);

				slic5_skip_lzpixels(&pr, src->data+idx, lzlen);
				//lzmap[hash]=idx>>2;
				kx+=lzlen;
				idx+=lzlen<<2;
			}
			else
			{
				if(hash>=0&&lzmap[hash]<0)
					lzmap[hash]=idx>>2;
#endif
			int kc=0;
			if(nch>=3)
			{
				int comp[]={src->data[idx|0], src->data[idx|1], src->data[idx|2]};
#ifdef SLIC5_OPTIMIZE_RCT
				orct_fwd(comp, rct_params, permutation);
#else
				rct_fwd(comp, rct);
#endif
				slic5_enc(&pr, comp[1], 0, kx);//Y		luma is encoded first
				slic5_enc(&pr, comp[2], 1, kx);//Cb
				slic5_enc(&pr, comp[0], 2, kx);//Cr
				kc+=3;
			}
			for(;kc<nch;++kc)//gray/alpha
			{
				int pixel=src->data[idx|kc];
				slic5_enc(&pr, pixel, kc, kx);
			}
			++kx;
			idx+=4;
#ifdef ENABLE_LZ
			}
#endif
		}
	}
	ac_enc_flush(&ec);
	dlist_appendtoarray(&list, data);
	if(curiosity)
	{
#ifndef SLIC5_OPTIMIZE_RCT
		for(int k=0;k<RCT_COUNT;++k)
			curiosity->rct_sizes[k]+=csizes[k];
#endif
		for(int k=0;k<SLIC5_NPREDS;++k)
			curiosity->pred_errors[k]+=pr.pred_error_sums[k];
	}
	if(loud)
	{
		double usize=image_getBMPsize(src);
		printf("\n");
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");

		printf("csize %8d  invCR %10.6lf%%  used %lf MB\n", (int)list.nobj, 100.*list.nobj/usize, (double)memusage/(1024*1024));

#ifdef PRINT_HIST
		const char *ch_labels[]={"Y", "Cb", "Cr", "A"};
		int ch_nhist=pr.nhist*pr.nhist, total_hist=nch*ch_nhist;
		printf("Hist CWH %d*%d*%d:\n", nch, pr.cdfsize, ch_nhist);
		for(int kt=0;kt<total_hist;++kt)
		{
			int *curr_hist=pr.hist+pr.cdfsize*kt;
			if(!(kt%ch_nhist))
				printf("%s  %d-bit\n", ch_labels[kt/ch_nhist], pr.depths[kt/ch_nhist]);
			for(int ks=0;ks<pr.cdfsize;++ks)
				printf("%d%c", curr_hist[ks], ks<pr.cdfsize-1?' ':'\n');
			if(!((kt+1)%pr.nhist)&&kt+1<total_hist)
				printf("\n");
		}
#endif

		printf("Pred errors\n");
		long long sum=0;
		for(int k=0;k<SLIC5_NPREDS;++k)
			sum+=pr.pred_error_sums[k];
		for(int k=0;k<SLIC5_NPREDS;++k)
			printf("  %2d  %17lld  %6.2lf%%  %s\n", k, pr.pred_error_sums[k], 100.*pr.pred_error_sums[k]/sum, slic5_prednames[k]);
		
#ifdef TRACK_SSE_RANGES
		printf("SSE ranges\n");
		for(int k=0;k<_countof(pr.sse_ranges)-1;k+=2)
		{
			int n=(&pr.sse_width)[k>>1];
			printf("  %c %4d ~ %4d / %4d  ", "XYZP"[k>>1], pr.sse_ranges[k], pr.sse_ranges[k+1], n);
			for(int k2=0;k2<n;++k2)
				printf("%c", k2>=pr.sse_ranges[k]&&k2<=pr.sse_ranges[k+1]?'*':'-');
			printf("\n");
		}
#endif

		//printf("Bias\n");
		//for(int kc=0;kc<nch;++kc)
		//	printf("  %lld/%d = %d\n", pr.bias_sum[kc], pr.bias_count[kc], pr.bias_count[kc]?(int)(pr.bias_sum[kc]/pr.bias_count[kc]):0);
		//for(int k=0;k<768;k+=4)
		//{
		//	int q=QUANTIZE_HIST(k);
		//	int q2=floor_log2_32(k)+1;
		//	printf("%3d  %3d  %3d\n", k, q, q2);
		//}

		prof_print();
	}
	slic5_free(&pr);
	dlist_clear(&list);
#ifdef ENABLE_LZ
	free(lzmap);
#endif
	return 1;
}
int t47_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	PROF_START();
	double t_start=time_sec();
#ifdef SLIC5_OPTIMIZE_RCT
	char rct_params[ORCT_NPARAMS+1]={0}, permutation[3]={0};
#else
	RCTType rct=RCT_NONE;
#endif
	int dst_nch=(dst->depth[0]!=0)+(dst->depth[1]!=0)+(dst->depth[2]!=0)+(dst->depth[3]!=0);
	int nch=MINVAR(dst_nch, dst->nch);
#ifdef SLIC5_OPTIMIZE_RCT
	char depths[4]={0};
	memcpy(depths, dst->depth, nch*sizeof(char));
#endif
	if(nch>=3)
	{
#ifdef SLIC5_OPTIMIZE_RCT
		if(srclen<ORCT_NPARAMS+1)
		{
			LOG_ERROR("Corrupt file");
			return 0;
		}
		memcpy(rct_params, data, sizeof(rct_params));
		if(rct_params[ORCT_NPARAMS]>=6)
		{
			LOG_ERROR("Corrupt file");
			return 0;
		}
		orct_unpack_permutation(rct_params[ORCT_NPARAMS], permutation);
		orct_calc_depths(rct_params, permutation, dst->depth, depths);
		data+=sizeof(rct_params);
		srclen-=sizeof(rct_params);
#else
		if(srclen<1)
		{
			LOG_ERROR("Corrupt file");
			return 0;
		}
		rct=*data;
		if((unsigned)rct>=(unsigned)RCT_COUNT)
		{
			LOG_ERROR("Corrupt file");
			return 0;
		}
		++data;
		--srclen;
#endif
	}
	ArithmeticCoder ec;
	SLIC5Ctx pr;
	ac_dec_init(&ec, data, data+srclen);
#ifdef SLIC5_OPTIMIZE_RCT
	ptrdiff_t memusage=slic5_init(&pr, dst->iw, dst->ih, nch, depths, &ec);
#else
	ptrdiff_t memusage=slic5_init(&pr, dst->iw, dst->ih, nch, dst->depth, &ec);
#endif
	if(!memusage)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int ky=0, idx=0;ky<dst->ih;++ky)
	{
		slic5_nextrow(&pr, ky);
		for(int kx=0;kx<dst->iw;)
		{
			//if(kx==252&&ky==16)//
			//if(kx==151&&ky==37)//
			//if(kx==100&&ky==38)//
			//if(kx==5&&ky==7)//
			//if(kx==588&&ky==254)//
			//	printf("");
			
			int comp[3];
			int kc=0;
			int first=slic5_dec(&pr, 0, kx);
#ifdef ENABLE_LZ
			if((unsigned)first==0x80000000)
			{
				int lzidx, lzlen;
				lzidx=(EC_DEC(&ec, 0, 256, 0)&0xFF)<<24;
				if(lzidx<0)
					lzidx=(idx>>2)-dst->iw;
				else
				{
					lzidx|=(EC_DEC(&ec, 0, 256, 0)&0xFF)<<16;
					lzidx|=(EC_DEC(&ec, 0, 256, 0)&0xFF)<<8;
					lzidx|=EC_DEC(&ec, 0, 256, 0)&0xFF;
				}
				lzlen=EC_DEC(&ec, 0, 256, 0)&0xFF;
				lzlen|=(EC_DEC(&ec, 0, 256, 0)&0xFF)<<8;
				slic5_update_hist(&pr, pr.cdfsize-1);
				slic5_update_CDFs(&pr);
				const int *srcptr=dst->data+((size_t)lzidx<<2);
				int *dstptr=dst->data+idx;

				for(int k=0;k<(lzlen<<2);k+=4)//can't use memcpy/memmove here, to support RLE
				{
					dstptr[k|0]=srcptr[k|0];
					dstptr[k|1]=srcptr[k|1];
					dstptr[k|2]=srcptr[k|2];
					dstptr[k|3]=srcptr[k|3];
#ifdef ENABLE_GUIDE
					if(guide&&(
						dstptr[k|0]!=guide->data[(guide->iw*ky+kx)*4+k+0]||
						dstptr[k|1]!=guide->data[(guide->iw*ky+kx)*4+k+1]||
						dstptr[k|2]!=guide->data[(guide->iw*ky+kx)*4+k+2]||
						dstptr[k|3]!=guide->data[(guide->iw*ky+kx)*4+k+3]
					))
						LOG_ERROR("Guide error XY %d %d", kx, ky);
#endif
				}
				slic5_skip_lzpixels(&pr, dst->data+idx, lzlen);
				kx+=lzlen;
				idx+=lzlen<<2;
			}
			else
			{
#endif
				if(nch>=3)
				{
					comp[1]=first;//Y
					comp[2]=slic5_dec(&pr, 1, kx);//Cb
					comp[0]=slic5_dec(&pr, 2, kx);//Cr
#ifdef SLIC5_OPTIMIZE_RCT
					orct_inv(comp, rct_params, permutation);
#else
					rct_inv(comp, rct);
#endif
					dst->data[idx|0]=comp[0];
					dst->data[idx|1]=comp[1];
					dst->data[idx|2]=comp[2];
					kc+=3;
#ifdef ENABLE_GUIDE
					if(guide&&(
						comp[0]!=guide->data[(guide->iw*ky+kx)<<2|0]||
						comp[1]!=guide->data[(guide->iw*ky+kx)<<2|1]||
						comp[2]!=guide->data[(guide->iw*ky+kx)<<2|2]
					))
						LOG_ERROR("Guide error XY %d %d", kx, ky);
#endif
				}
				while(kc<nch)//gray/alpha
				{
					int pixel=kc?slic5_dec(&pr, kc, kx):first;
#ifdef ENABLE_GUIDE
					if(guide&&pixel!=guide->data[(guide->iw*ky+kx)<<2|kc])
						LOG_ERROR("Guide error XY %d %d", kx, ky);
#endif
					dst->data[idx|kc]=pixel;
					++kc;
				}
				++kx;
				idx+=4;
#ifdef ENABLE_LZ
			}
#endif
		}
	}
	if(dst_nch>=3&&nch==1)//grayscale fix
	{
		int res=dst->iw*dst->ih;
		for(int k=0;k<res;++k)
		{
			int val=dst->data[k<<2|0];
			dst->data[k<<2|1]=val;
			dst->data[k<<2|2]=val;
		}
	}
	if(loud)
	{
		printf("\n");
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		prof_print();
	}
	slic5_free(&pr);
	return 1;
}

#if 0
int t47_from_ppm(const char *src, const char *dst)
{
	double t=time_sec();
	FILE *f=fopen(src, "rb");
	if(!f)
	{
		printf("Cannot open \"%s\"\n", src);
		return 0;
	}
	char buf[1024]={0};
	int len=0, idx=0;
	int iw, ih, maxval;
	len=(int)fread(buf, 1, 1024, f);
	if(memcmp(buf, "P6", 2))
	{
		fclose(f);
		printf("Expected a Binary Portable PixMap (P6 header)\n");
		return 0;
	}
	idx=2;
	
	char *end=0;
	for(;idx<1024&&isspace(buf[idx]);++idx);//skip newline after 'P6'
	iw=strtol(buf+idx, &end, 10);
	idx=(int)(end-buf);

	for(;idx<1024&&isspace(buf[idx]);++idx);//skip space after width
	ih=strtol(buf+idx, &end, 10);
	idx=(int)(end-buf);

	for(;idx<1024&&isspace(buf[idx]);++idx);//skip newline after height
	maxval=strtol(buf+idx, &end, 10);
	idx=(int)(end-buf);
	if(maxval!=255)
	{
		fclose(f);
		printf("Expected an 8-bit image\n");
		return 0;
	}
	for(;idx<1024&&isspace(buf[idx]);++idx);//skip newline after maxval
	
	DList list;
	ArithmeticCoder ec;
	SLIC5Ctx pr;
	dlist_init(&list, 1, 65536, 0);
	ac_enc_init(&ec, &list);
	char depths[]={8, 8, 8, 0};
	int success=slic5_init(&pr, iw, ih, 3, depths, &ec);
	if(!success)
	{
		fclose(f);
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int ky=0;ky<ih;++ky)
	{
		slic5_nextrow(&pr, ky);
		for(int kx=0;kx<iw;++kx)
		{
			int r, g, b;
#define CHECK_BUFFER()\
	if(idx>=len)\
	{\
		len=(int)fread(buf, 1, 1024, f);\
		idx=0;\
		if(idx>=len)\
		{\
			LOG_ERROR("Alloc error");\
			return 0;\
		}\
	}
			CHECK_BUFFER()
			r=(unsigned char)buf[idx++]-128;
			CHECK_BUFFER()
			g=(unsigned char)buf[idx++]-128;
			CHECK_BUFFER()
			b=(unsigned char)buf[idx++]-128;

			r-=g;
			b-=g;
			g+=(r+b)>>2;

			slic5_enc(&pr, g, 0, kx);//Y
			slic5_enc(&pr, b, 1, kx);//Cb
			slic5_enc(&pr, r, 2, kx);//Cr
		}
		printf("%6.2lf%%  CR %10lf  Elapsed %12lf\r", 100.*(ky+1)/ih, 3.*iw*(ky+1)/list.nobj, time_sec()-t);
	}
	ac_enc_flush(&ec);

	printf("\n");

	ArrayHandle cdata=0;
	lsim_writeheader(&cdata, iw, ih, 3, depths, 47);
	dlist_appendtoarray(&list, &cdata);
	
	t=time_sec()-t;
	printf("Encoded  %lf sec  CR 3*%d*%d/%d = %lf\n", t, iw, ih, (int)list.nobj, 3.*iw*ih/list.nobj);
	
	slic5_free(&pr);
	dlist_clear(&list);
	fclose(f);

	success=save_file(dst, cdata->data, cdata->count, 1);
	printf("%s \"%s\"\n", success?"Saved":"Failed to save", dst);

	array_free(&cdata);
	return 1;
}
int t47_to_ppm(const char *src, const char *dst)
{
	double t=time_sec();
	ArrayHandle cdata=load_file(src, 1, 0, 0);
	const unsigned char *data=cdata->data;
	size_t srclen=cdata->count;
	LSIMHeader header;
	size_t bytesread=lsim_readheader(data, srclen, &header);
	if(header.iw<1||header.ih<1||header.nch!=3||header.depth[0]!=8||header.depth[1]!=8||header.depth[2]!=8)
	{
		printf("Expected an 8-bit RGB image\n");
		return 0;
	}
	data+=bytesread;
	srclen-=bytesread;
	FILE *f=fopen(dst, "wb");
	if(!f)
	{
		printf("Failed to save \"%s\"\n", dst);
		return 0;
	}
	int written=snprintf(g_buf, G_BUF_SIZE, "P6\n%d %d\n%d\n", header.iw, header.ih, (1<<header.depth[0])-1);
	written=(int)fwrite(g_buf, 1, written, f);
	if(!written)
	{
		printf("File write error");
		return 0;
	}
	ArithmeticCoder ec;
	SLIC5Ctx pr;
	ac_dec_init(&ec, data, data+srclen);
	int success=slic5_init(&pr, header.iw, header.ih, header.nch, header.depth, &ec);
	if(!success)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int ky=0;ky<header.ih;++ky)
	{
		slic5_nextrow(&pr, ky);
		for(int kx=0;kx<header.iw;++kx)
		{
			int g=slic5_dec(&pr, 0, kx);//Y
			int b=slic5_dec(&pr, 1, kx);//Cb
			int r=slic5_dec(&pr, 2, kx);//Cr

			g-=(r+b)>>2;
			b+=g;
			r+=g;

			char triplet[]={r+128, g+128, b+128, 0};
			written=(int)fwrite(triplet, 1, 3, f);
			if(!written)
			{
				printf("File write error");
				return 0;
			}
		}
		printf("%6.2lf%%  Elapsed %12lf\r", 100.*(ky+1)/header.ih, time_sec()-t);
	}

	printf("\n");

	t=time_sec()-t;
	printf("Decoded  %lf sec\n", t);
	
	slic5_free(&pr);
	fclose(f);
	return 1;
}

#define LG_HIST_SIZE 8
void t47_analyze_preds(const char *path)
{
	double t_start=time_sec();
	printf("Predictor Analysis\n");
	static const char *ext[]=
	{
		"png",
		"jpg",
		"jpeg",
		"ppm",
		"pgm",
	};
	ArrayHandle filenames=get_filenames(path, ext, _countof(ext), 1);
	if(!filenames)
	{
		printf("No supported images in \"%s\"\n", path);
		return;
	}
	long long *hist=(long long*)malloc(sizeof(long long[SLIC5_NPREDS<<(LG_HIST_SIZE<<1)]));
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(hist, 0, sizeof(long long[SLIC5_NPREDS<<(LG_HIST_SIZE<<1)]));
	int histcount=sizeof(long long[SLIC5_NPREDS<<(LG_HIST_SIZE<<1)])/sizeof(long long);
	SLIC5Ctx pr;
#if 0
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k);
		printf("%s\n", (char*)fn[0]->data);
		const char *p2="C:/Projects/datasets/dataset-LPCB/";
		if(memcmp(fn[0]->data, p2, strlen(p2)))//
			LOG_ERROR("%s", fn[0]->data);
	}
#endif
	for(int k=0;k<(int)filenames->count;++k)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k);
#if 0
		const char *p2="C:/Projects/datasets/dataset-LPCB/";
		if(memcmp(fn[0]->data, p2, strlen(p2)))//
			LOG_ERROR("%s", fn[0]->data);
#endif
		printf("%5d/%5d  %6.2lf%%  %12lf  \"%s\"\n", k+1, (int)filenames->count, 100.*(k+1)/filenames->count, time_sec()-t_start, (char*)fn[0]->data);
		ptrdiff_t formatsize=get_filesize((char*)fn[0]->data);
		if(formatsize<1)//skip non-images, this check is useless because get_filenames() has already filtered the list
			continue;
		Image *src=image_load((char*)fn[0]->data);
		if(!src)
		{
			printf("\nCannot open \"%s\"\n", (char*)fn[0]->data);
			continue;
		}
		int success=slic5_init(&pr, src->iw, src->ih, 1, src->depth, 0);
		//int ucount=0;
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			slic5_nextrow(&pr, ky);
			for(int kx=0;kx<src->iw;++kx, idx+=4)
			{
				//if(kx==src->iw/2&&ky==src->ih/2)//
				//	printf("");
				int val=src->data[idx];
				slic5_predict(&pr, 0, kx);

				int y0=0;
				for(int kp=0;kp<SLIC5_NPREDS;++kp)
				{
					int idx2=(size_t)kp<<(LG_HIST_SIZE<<1);

					if((idx2|((1<<(LG_HIST_SIZE<<1))-1))>=histcount)//
						LOG_ERROR("IDX %#08X/%#08X", idx2, histcount);
					
					long long *curr_hist=hist+idx2;
					int x=val>>(pr.shift_prec[0]-PRED_PREC+8-LG_HIST_SIZE);
					int y=pr.preds[kp]>>(pr.shift_prec[0]+8-LG_HIST_SIZE);
					x+=(1<<LG_HIST_SIZE)>>1;
					y+=(1<<LG_HIST_SIZE)>>1;
					x=CLAMP(0, x, (1<<LG_HIST_SIZE)-1);
					y=CLAMP(0, y, (1<<LG_HIST_SIZE)-1);
					//if((unsigned)x>=(unsigned)(1<<LG_HIST_SIZE)&&(unsigned)y>=(unsigned)(1<<LG_HIST_SIZE))
					//	LOG_ERROR("XY %#08X %#08X", x, y);
					if((unsigned)(y<<LG_HIST_SIZE|x)>=(1<<(LG_HIST_SIZE<<1)))
						LOG_ERROR("IDX %#08X/%#08X", y<<LG_HIST_SIZE|x, 1<<(LG_HIST_SIZE<<1));
					//if(!kp)
					//	y0=y;
					//else
					//	ucount+=y!=y0;
					++curr_hist[y<<LG_HIST_SIZE|x];
				}
				slic5_update(&pr, val, -1);
			}
		}
		//printf("%d/15\n", ucount);
		slic5_free(&pr);
		free(src);
	}
	printf("\nFinished Processing\n");
#if 1
	printf("Are histograms identical?\n");
	for(int kp=1;kp<SLIC5_NPREDS;++kp)
	{
		long long *curr_hist=hist+((size_t)kp<<(LG_HIST_SIZE<<1));
		int result=memcmp(hist, curr_hist, sizeof(long long[1<<(LG_HIST_SIZE<<1)]));
		printf("memcmp %2d: %d\n", kp, result);
	}
#endif
	printf("Histograms:\n");
	for(int kp=0;kp<SLIC5_NPREDS;++kp)
	{
		long long *curr_hist=hist+((size_t)kp<<(LG_HIST_SIZE<<1));
		printf("%s\n", slic5_prednames[kp]);
		long long vmax=0;
		//for(int k=0;k<(1<<(LG_HIST_SIZE<<1));++k)
		//{
		//	UPDATE_MAX(vmax, curr_hist[k]);
		//}
		for(int ky=0;ky<(1<<LG_HIST_SIZE);++ky)
		{
			vmax=0;
			for(int kx=0;kx<(1<<LG_HIST_SIZE);++kx)
			{
				UPDATE_MAX(vmax, curr_hist[ky<<LG_HIST_SIZE|kx]);
			}
			for(int kx=0;kx<(1<<LG_HIST_SIZE);++kx)
				printf("%X", vmax?(int)(curr_hist[ky<<LG_HIST_SIZE|kx]*15/vmax):0);
			printf("  %16lld\n", vmax);
		}
		printf("\n");
	}
	printf("Done.\n");
	array_free(&filenames);
	free(hist);
	pause();
}
#endif