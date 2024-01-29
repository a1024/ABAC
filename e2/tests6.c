#include"e2.h"
#define AC_IMPLEMENTATION
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctype.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;

//	#define AVX512
	#define AVX2
//	#define PROBBITS_15
//	#define PROFILER//SLOW
//	#define TRACK_SSE_RANGES//SLOW
//	#define DISABLE_SSE

#if defined AVX2 || defined AVX512
#include<immintrin.h>
#endif
#ifdef PROFILER
#define CHECKPOINTLIST\
	CHECKPOINT(OUTSIDE)\
	CHECKPOINT(UPDATE_CDFs)\
	CHECKPOINT(FETCH_NB)\
	CHECKPOINT(CALC_SUBPREDS)\
	CHECKPOINT(CALC_WEIGHT_AV)\
	CHECKPOINT(CALC_CTX)\
	CHECKPOINT(SSE_LOOP)\
	CHECKPOINT(PRED_TILL_UPDATE)\
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
#define PROF_START() memset(prof_cycles, 0, _countof(prof_cycles)), prof_timestamp=__rdtsc()
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

#define SLIC5_CONFIG_EXP 5
#define SLIC5_CONFIG_MSB 2
#define SLIC5_CONFIG_LSB 0

#define CDF_UPDATE_PERIOD 0x400//number of sub-pixels processed between CDF updates, must be a power-of-two
#define PAD_SIZE 2
#define NPREDS 13
#define PRED_PREC 8
#define PARAM_PREC 8
#define HIST_EXP 2
#define HIST_MSB 1
#define SSE_X_EXP 1
#define SSE_X_MSB 0
#define SSE_Y_EXP 1
#define SSE_Y_MSB 0
#define SSE_Z_EXP 1
#define SSE_Z_MSB 0
#define SSE_W 7		//ENABLE SSE INDEX CLAMP WHAN CHANGING CONFIG
#define SSE_H 21
#define SSE_D 21
#define SSE_PREDBITS 5
#define SSE_PRED_LEVELS (1<<SSE_PREDBITS)//(_countof(qlevels_pred)+1)
#define SSE_FR_SIZE 59049//3^10		separate final round
//#define SSE_FR_SIZE (1<<10)
#define SSE_STAGES 10
#define SSE_SIZE (SSE_W*SSE_H*SSE_D*SSE_PRED_LEVELS)
//#define HASH_CANTOR(A, B) (((A)+(B)+1)*((A)+(B))/2+(B))
#ifdef PROBBITS_15
#define EC_ENC(EC, X, CDF, NLEVELS, LG_FMIN) ac_enc15(EC, X, CDF, NLEVELS)
#define EC_DEC(EC, CDF, NLEVELS, LG_FMIN) ac_dec15(EC, CDF, NLEVELS)
typedef unsigned short CDF_t;
#define CDF_SHIFT 15
#else
#define EC_ENC(EC, X, CDF, NLEVELS, LG_FMIN) ac_enc(EC, X, CDF, NLEVELS, LG_FMIN)
#define EC_DEC(EC, CDF, NLEVELS, LG_FMIN) ac_dec(EC, CDF, NLEVELS, LG_FMIN)
typedef unsigned CDF_t;
#define CDF_SHIFT 16
#endif
typedef struct HybridUintStruct
{
	unsigned short token, nbits;
	unsigned bypass;
} HybridUint;

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
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
		int lgv=floor_log2_32((unsigned)val);
		int mantissa=val-(1<<lgv);
		token = (1<<SLIC5_CONFIG_EXP) + (
				(lgv-SLIC5_CONFIG_EXP)<<(SLIC5_CONFIG_MSB+SLIC5_CONFIG_LSB)|
				(mantissa>>(lgv-SLIC5_CONFIG_MSB))<<SLIC5_CONFIG_LSB|
				(mantissa&((1<<SLIC5_CONFIG_LSB)-1))
			);
		nbits=lgv-(SLIC5_CONFIG_MSB+SLIC5_CONFIG_LSB);
		bypass=val>>SLIC5_CONFIG_LSB&((1LL<<nbits)-1);
	}
	hu->token=token;
	hu->bypass=bypass;
	hu->nbits=nbits;
}
static int quantize_unsigned(int val, int exp, int msb)
{
	if(val<(1<<exp))
		return val;
	int lgv=floor_log2_32(val);
	int token=(1<<exp)+((lgv-exp)<<msb|(val-(1<<lgv))>>(lgv-msb));
	return token;
}
static int quantize_signed_get_range(int num, int den, int exp, int msb)
{
	int vmax=(int)((1LL<<24)*num/den>>16);
	int token=quantize_unsigned(vmax, exp, msb);
	token<<=1;
	return token;
}
static int quantize_signed(int val, int shift, int exp, int msb, int nlevels)
{
	val>>=shift;
	int negmask=-(val<0);
	int token=quantize_unsigned(abs(val), exp, msb);
	//if((unsigned)token>=(unsigned)(nlevels>>1))//
	//	LOG_ERROR("SSE range error");
	token^=negmask;
	token-=negmask;
	token+=nlevels>>1;
	//token=CLAMP(0, token, nlevels-1);
	return token;
}
#define QUANTIZE_HIST(X) (quantize_unsigned(X, HIST_EXP, HIST_MSB)>>1)
//static int QUANTIZE_HIST(int x)
//{
//	int token=quantize_unsigned(x, HIST_EXP, HIST_MSB);
//
//	//if(token<(1<<HIST_EXP))
//	//	token>>=1;
//	//else
//	//	token-=(1<<HIST_EXP)>>1;
//
//	token>>=1;
//	//token-=token!=0;
//	
//	//if(token<128)//ORIGINAL
//	//	token>>=1;
//	//else
//	//	token-=64;
//	//
//	return token;
//}
#if defined AVX2 || defined AVX512
static void avx2_floor_log2_p1(__m256i *x)//floor_log2()+1
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
//#define QUANTIZE_HIST(X) quantize_unsigned(X, HIST_EXP, HIST_MSB)
static const int wp_params[]=
{
	0x0004AB65,//wp1 X		bad				2nd worst	worst
	0x0002841C,//wp2						3rd best
	0x00034147,//wp3				best		2nd best	best		2nd best
	0x00016056,//wp4
	0x00016056,//(W+NEE)>>1		best				best		2nd best	best
	0x00016056,//N+NE-NNE		2nd best	2nd best			2nd worst	worst
	//0x00016056,//N+NW-NNW	X	worst
	//0x00016056,//W+NW-NWW X	2nd worst	2nd worst	worst
	//0x00016056,//N+NEE-NNEE X	very bad	worst
	0x00016056,//(N+W)>>1
	0x00016056,//(4*(W+NW+N+NE)-(WW+NNWW+NN+NNEE)+6)/12
	0x00016056,//(W+2*NE-NNE)>>1
	//0x00016056,//(2*(W+N+NE)-(NNW+NNE))>>2 X
	//0x00016056,//NW+NE-NN X
	0x00016056,//(N+NN)>>1
	0x00016056,//grad								3rd best
	0x00028F8A,//paper GAP		bad
	0x00022F58,//CALIC GAP		bad		3rd worst	3rd worst	3rd worst	2nd worst

	-0x005C,
	-0x005B,
	 0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,

	//0x80,//X
	//0x80,
	//0x40, 0x00, 0x30, 0x170, 0x20,

	//0xd, 0xc, 0xc, 0xb,
	//8,
	//8,
	//4, 0, 3, 23, 2,
};
typedef struct SLIC5CtxStruct
{
	int iw, ih, nch;
	char depths[4];
	int
		half[4],
		shift_prec[4],//shift right to bring value to 8-bit, removes predictor-added bits
		shift_error[4];
	int nhist, cdfsize;
	//int sse_width, sse_height, sse_depth, sse_nplanes, sse_planesize, sse_size;
	int *pred_errors, *errors, *pixels;

	int *hist, *histsums;
	CDF_t *CDFs;

#ifndef DISABLE_SSE
	long long *sse, *sse_fr;
	int sse_idx[SSE_STAGES+1];
	long long sse_sum[SSE_STAGES+1];
	int sse_count[SSE_STAGES+1];
#endif
	int bias_count[4];
	long long bias_sum[4];

	int preds[NPREDS], params[_countof(wp_params)];
	long long pred;
	int pred_final;
	int hist_idx;
	int sse_corr;
	int kc, kx, ky, kym0, kym1, kym2;
	long long pred_error_sums[NPREDS];
#ifdef TRACK_SSE_RANGES
	int sse_ranges[8];
#endif
	ArithmeticCoder *ec;
} SLIC5Ctx;
#define LOAD(BUF, X, Y) BUF[(pr->kym##Y+kx-(X)+PAD_SIZE)<<2|kc]
static int slic5_init(SLIC5Ctx *pr, int iw, int ih, int nch, const char *depths, ArithmeticCoder *ec)
{
	PROF_START();
	if(iw<1||ih<1||nch<1)
	{
		LOG_ERROR("Invalid image");
		return 0;
	}
	memset(pr, 0, sizeof(*pr));
	pr->iw=iw;
	pr->ih=ih;
	pr->nch=nch;
	if(nch<3)
		memcpy(pr->depths, depths, nch);
	else
	{
		pr->depths[0]=depths[1];  //Y
		pr->depths[1]=depths[2]+1;//Cb
		pr->depths[2]=depths[0]+1;//Cr
		if(nch==4)
			pr->depths[3]=depths[3];//a
	}
	int maxdepth=0;
	for(int kc=0;kc<nch;++kc)
	{
		int nlevels=1<<pr->depths[kc], prec_half=(1<<(pr->depths[kc]+PRED_PREC))>>1;
		pr->half[kc]=nlevels>>1;
		if(maxdepth<pr->depths[kc])
			maxdepth=pr->depths[kc];
	}
	for(int kc=0;kc<pr->nch;++kc)
	{
		pr->shift_prec[kc]=MAXVAR(8, pr->depths[kc])-8+PRED_PREC;
		pr->shift_error[kc]=1+((pr->depths[kc]<=9)<<1);
	}
	HybridUint hu;
	int extremesym=-((1<<maxdepth)>>1);
	extremesym=extremesym<<1^-(extremesym<0);//pack sign
	hybriduint_encode(extremesym, &hu);//encode -half
	pr->cdfsize=hu.token+1;
	pr->nhist=QUANTIZE_HIST(767);
	//pr->nhist=1<<3;//X

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

	pr->pred_errors=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[NPREDS*4*4]));//NPREDS * 4 rows * 4 comps
	pr->errors=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	pr->pixels=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	pr->hist=(int*)malloc((size_t)nch*pr->nhist*pr->nhist*pr->cdfsize*sizeof(int));//WH: cdfsize * NHIST
	pr->histsums=(int*)malloc((size_t)nch*pr->nhist*pr->nhist*sizeof(int));
	pr->CDFs=(CDF_t*)malloc((size_t)nch*pr->nhist*pr->nhist*(pr->cdfsize+1LL)*sizeof(CDF_t));
#ifndef DISABLE_SSE
	pr->sse=(long long*)malloc(nch*sizeof(long long[SSE_SIZE*SSE_STAGES]));
	pr->sse_fr=(long long*)malloc(nch*sizeof(long long[SSE_FR_SIZE]));
#endif
	if(!pr->pred_errors||!pr->errors||!pr->pixels||!pr->hist||!pr->histsums||!pr->CDFs
#ifndef DISABLE_SSE
		||!pr->sse||!pr->sse_fr
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(pr->pred_errors, 0, (iw+PAD_SIZE*2LL)*sizeof(int[NPREDS*4*4]));
	memset(pr->errors, 0, (iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	memset(pr->pixels, 0, (iw+PAD_SIZE*2LL)*sizeof(int[4*4]));

	memset(pr->hist, 0, (size_t)nch*pr->nhist*pr->nhist*pr->cdfsize*sizeof(int));
	memset(pr->histsums, 0, (size_t)nch*pr->nhist*pr->nhist*sizeof(int));

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
	memfill(pr->CDFs+pr->cdfsize+1, pr->CDFs, ((size_t)nch*pr->nhist*pr->nhist-1LL)*(pr->cdfsize+1LL)*sizeof(CDF_t), (pr->cdfsize+1LL)*sizeof(CDF_t));
	
#ifndef DISABLE_SSE
	memset(pr->sse, 0, nch*sizeof(long long[SSE_SIZE*SSE_STAGES]));
	memset(pr->sse_fr, 0, nch*sizeof(long long[SSE_FR_SIZE]));
#endif

	int shift=16-maxdepth;
	shift=MAXVAR(0, shift);
	memcpy(pr->params, wp_params, sizeof(pr->params));
	for(int k=0;k<NPREDS;++k)
		pr->params[k]=pr->params[k]>>shift;

	pr->ec=ec;
	return 1;
}
static void slic5_free(SLIC5Ctx *pr)
{
	free(pr->pred_errors);
	free(pr->errors);
	free(pr->pixels);
	free(pr->hist);
	free(pr->histsums);
	free(pr->CDFs);
#ifndef DISABLE_SSE
	free(pr->sse);
	free(pr->sse_fr);
#endif
}
static void slic5_update_CDFs(SLIC5Ctx *pr)
{
	int nhist=pr->nch*pr->nhist*pr->nhist;
	for(int kt=0;kt<nhist;++kt)//update CDFs (deferred)
	{
		int sum=pr->histsums[kt];
		if(sum>pr->cdfsize)//only when the CDF is hit enough times
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
static void slic5_nextrow(SLIC5Ctx *pr, int ky)
{
	pr->ky=ky;
	pr->kym0=(pr->iw+PAD_SIZE*2)*(ky&3);
	pr->kym1=(pr->iw+PAD_SIZE*2)*((ky-1)&3);
	pr->kym2=(pr->iw+PAD_SIZE*2)*((ky-2)&3);
}
static void slic5_predict(SLIC5Ctx *pr, int kc, int kx, int ky)
{
	PROF(OUTSIDE);
	int idx=(pr->iw*ky+kx)<<2|kc;
	if(!(idx&(CDF_UPDATE_PERIOD-1))||idx<CDF_UPDATE_PERIOD&&(idx&15))
		slic5_update_CDFs(pr);
	PROF(UPDATE_CDFs);
	//XY are flipped, no need to check if indices OOB due to padding
	int
		NNWW	=LOAD(pr->pixels,  2, 2)<<PRED_PREC,
		NNW	=LOAD(pr->pixels,  1, 2)<<PRED_PREC,
		NN	=LOAD(pr->pixels,  0, 2)<<PRED_PREC,
		NNE	=LOAD(pr->pixels, -1, 2)<<PRED_PREC,
		NNEE	=LOAD(pr->pixels, -2, 2)<<PRED_PREC,
		NWW	=LOAD(pr->pixels,  2, 1)<<PRED_PREC,
		NW	=LOAD(pr->pixels,  1, 1)<<PRED_PREC,
		N	=LOAD(pr->pixels,  0, 1)<<PRED_PREC,
		NE	=LOAD(pr->pixels, -1, 1)<<PRED_PREC,
		NEE	=LOAD(pr->pixels, -2, 1)<<PRED_PREC,
		WW	=LOAD(pr->pixels,  2, 0)<<PRED_PREC,
		W	=LOAD(pr->pixels,  1, 0)<<PRED_PREC;
	int
		eNNWW	=LOAD(pr->errors,  2, 2),//error = (curr<<8) - pred
		eNNW	=LOAD(pr->errors,  1, 2),
		eNN	=LOAD(pr->errors,  0, 2),
		eNNE	=LOAD(pr->errors, -1, 2),
		eNNEE	=LOAD(pr->errors, -2, 2),
		eNWW	=LOAD(pr->errors,  2, 1),
		eNW	=LOAD(pr->errors,  1, 1),
		eN	=LOAD(pr->errors,  0, 1),
		eNE	=LOAD(pr->errors, -1, 1),
		eNEE	=LOAD(pr->errors, -2, 1),
		eWW	=LOAD(pr->errors,  2, 0),
		eW	=LOAD(pr->errors,  1, 0);
	int sh=pr->shift_prec[kc];
	int clamp_lo=(int)N, clamp_hi=(int)N;
	clamp_lo=(int)MINVAR(clamp_lo, W);
	clamp_hi=(int)MAXVAR(clamp_hi, W);
	clamp_lo=(int)MINVAR(clamp_lo, NE);
	clamp_hi=(int)MAXVAR(clamp_hi, NE);
	PROF(FETCH_NB);

	int j=-1;
	pr->preds[++j]=W+NE-N-(int)((2*(eN+eW)+eNE-eNW+4)>>3);
	pr->preds[++j]=N-(int)((eN+eW+eNE)*pr->params[NPREDS+0]>>PARAM_PREC);
	pr->preds[++j]=W-(int)((eN+eW+eNW)*pr->params[NPREDS+1]>>PARAM_PREC);
	pr->preds[++j]=N-(int)((
		(long long)eNW*pr->params[NPREDS+2]+
		(long long)eN*pr->params[NPREDS+3]+
		(long long)eNE*pr->params[NPREDS+4]+
		((long long)NN-N)*pr->params[NPREDS+5]+
		((long long)NW-W)*pr->params[NPREDS+6]
	)>>PARAM_PREC);

	pr->preds[++j]=(W+NEE)>>1;
	pr->preds[++j]=N+NE-NNE-eNNE;
	//pr->preds[++j]=4*NE-(N+NEE+NNEE);
	pr->preds[++j]=(N+W)>>1;
	pr->preds[++j]=(4*(N+W+NW+NE)-(NN+WW+NNWW+NNEE)+6)/12;
	//pr->preds[++j]=(int)((181*(4*((long long)N+W)-((long long)NN+WW))+128*(4*((long long)NW+NE)-((long long)NNWW+NNEE))+1854)/3708);//X
	
	pr->preds[++j]=(W+2*NE-NNE+eNE)>>1;
	//pr->preds[++j]=(2*(W+N+NE)-(NNW+NNE))>>2;//X
	//pr->preds[++j]=(2*(W+NEE)-NNEE)/3;//X
	//pr->preds[++j]=(W+2*N-NNW)>>1;//X
	pr->preds[++j]=N+W-NW;

	pr->preds[++j]=2*W-WW+eW-eWW;

#if 0
	//pr->preds[++j]=(((N+W)<<1)+NW+NE+NN+WW+4)>>3;//X

	pr->preds[++j]=(int)(W -(eW *pr->params[NPREDS+ 7]>>8));
	pr->preds[++j]=(int)(N -(eN *pr->params[NPREDS+ 8]>>8));
	pr->preds[++j]=(int)(NW-(eNW*pr->params[NPREDS+ 9]>>8));
	pr->preds[++j]=(int)(NE-(eNE*pr->params[NPREDS+10]>>8));

	int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
	int vmin2=MINVAR(vmin, NW), vmax2=MAXVAR(vmax, NW);
	vmin2=MINVAR(vmin2, NE), vmax2=MAXVAR(vmax, NE);
	++j, pr->preds[j]=N+W-NW, pr->preds[j]=CLAMP(vmin2, pr->preds[j], vmax2);
	++j, pr->preds[j]=W+NE-N, pr->preds[j]=CLAMP(vmin2, pr->preds[j], vmax2);
	++j, pr->preds[j]=N+NW-NNW, pr->preds[j]=CLAMP(vmin2, pr->preds[j], vmax2);
	vmin2=MINVAR(vmin, NE), vmax2=MAXVAR(vmax, NE);
	vmin2=MINVAR(vmin2, NEE), vmax2=MAXVAR(vmax, NEE);
	++j, pr->preds[j]=N+NE-NNE, pr->preds[j]=CLAMP(vmin2, pr->preds[j], vmax2);
	
	++j, pr->preds[j]=(W+NEE)/2;
	++j, pr->preds[j]=NNNNNN;
	++j, pr->preds[j]=(NEEEE+NEEEEEE)/2;
	++j, pr->preds[j]=(WWWW+WWWWWW)/2;
	++j, pr->preds[j]=(N+W+NEEEEE+NEEEEEEE)/4;
	++j, pr->preds[j]=N*2-NN, pr->preds[j]=CLAMP(vmin2, pr->preds[j], vmax2);
	++j, pr->preds[j]=(N+NNN)/2;
	++j, pr->preds[j]=((N+W)*3-NW*2)/4;
	++j, pr->preds[j]=(int)(N+W-NW), pr->preds[j]=(int)CLAMP(vmin, pr->preds[j], vmax);
#endif

	int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
	int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
	int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
	int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
	int sum2=(dy+dx)>>sh, diff=(dy-dx)>>sh, diff2=(d45-d135)>>sh, diff3=(NE-NW)>>3;
	++j;
	//if(sum2>32)
	//	//pr->preds[j]=(int)(((long long)dx*N+(long long)dy*W)/sum2)+diff2;
	//	pr->preds[j]=(dx*N+dy*W)/sum2+diff2;
	//else if(diff>12)
	//	pr->preds[j]=(N+2*W)/3+diff2;
	//else if(diff<-12)
	//	pr->preds[j]=(2*N+W)/3+diff2;
	//else
		pr->preds[j]=((N+W)>>1)+diff3;
	if(diff2>32)
		pr->preds[j]+=diff3;
	else if(diff2>16)
		pr->preds[j]+=diff3>>1;
	else if(diff2<-32)
		pr->preds[j]-=diff3;
	else if(diff2<-16)
		pr->preds[j]-=diff3>>1;

	++j;
	if(diff>80)
		pr->preds[j]=W;
	else if(diff<-80)
		pr->preds[j]=N;
	else
	{
		pr->preds[j]=(((N+W)<<1)+NE-NW)>>2;
		if(diff>32)
			pr->preds[j]=(pr->preds[j]+W)>>1;
		else if(diff>8)
			pr->preds[j]=(3*pr->preds[j]+W)>>2;
		else if(diff<-32)
			pr->preds[j]=(pr->preds[j]+N)>>1;
		else if(diff<-8)
			pr->preds[j]=(3*pr->preds[j]+N)>>2;
	}
	PROF(CALC_SUBPREDS);

	long long weights[NPREDS]={0}, wsum=0;
	for(int k=0;k<NPREDS;++k)
	{
		weights[k]=
			(long long)pr->pred_errors[(NPREDS*(pr->kym2+kx+2+PAD_SIZE)+k)<<2|kc]+//peNNEE
			(long long)pr->pred_errors[(NPREDS*(pr->kym2+kx+1+PAD_SIZE)+k)<<2|kc]+//peNNE
			(long long)pr->pred_errors[(NPREDS*(pr->kym2+kx  +PAD_SIZE)+k)<<2|kc]+//peNN+peNW
			(long long)pr->pred_errors[(NPREDS*(pr->kym2+kx-1+PAD_SIZE)+k)<<2|kc]+//peNNW+peNWW
			(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx+2+PAD_SIZE)+k)<<2|kc]+//peNEE [not in jxl]
			(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx+1+PAD_SIZE)+k)<<2|kc]+//peNE
			(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx  +PAD_SIZE)+k)<<2|kc]*3+//peN+peW
			(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx-1+PAD_SIZE)+k)<<2|kc];//peNW+peWW
		weights[k]=((long long)pr->params[k]<<16)/(weights[k]+1);
		wsum+=weights[k];
	}
	if(wsum)
	{
		pr->pred=0;
		for(int k=0;k<NPREDS;++k)
			pr->pred+=weights[k]*pr->preds[k];
		pr->pred+=(wsum>>1)-1;
		pr->pred/=wsum;
	}
	else
		pr->pred=pr->preds[0];
	PROF(CALC_WEIGHT_AV);

#if 1
	int qx=QUANTIZE_HIST((dx+abs(eW+eWW))>>(sh-2)), qy=QUANTIZE_HIST((dy+abs(eN+eNN))>>(sh-2));

	//int qx=dx>>(sh+8-3), qy=dy>>(sh+8-3);//X

	//int hist_sh=sh;
	//hist_sh-=pr->depths[kc]==9;
	//hist_sh-=pr->depths[kc]>9;
	//int qx=quantize_unsigned(dx>>hist_sh, HIST_EXP, HIST_MSB)>>(pr->depths[kc]>9);
	//int qy=quantize_unsigned(dy>>hist_sh, HIST_EXP, HIST_MSB)>>(pr->depths[kc]>9);

	qx=CLAMP(0, qx, pr->nhist-1);
	qy=CLAMP(0, qy, pr->nhist-1);
	pr->hist_idx=pr->nhist*qy+qx;
#endif
	//pr->hist_idx=QUANTIZE_HIST((dx+abs(eW+eWW)+dy+abs(eN+eNN))>>(sh-2));
	//pr->hist_idx=CLAMP(0, pr->hist_idx, pr->nhist*pr->nhist-1);

#ifndef DISABLE_SSE
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
		(eW-eNW),		//8	2/1
		(eW-eWW),		//7	2/1
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
	_mm256_store_si256((__m256i*)g+0, X0);
	_mm256_store_si256((__m256i*)g+1, X1);
#else
	int g[]=
	{
		//4D SSE
#define QUANTIZE(X, Y, Z)\
	SSE_W*(SSE_H*quantize_signed(Z, sh, SSE_Z_EXP, SSE_Z_MSB, SSE_D)+\
	quantize_signed(Y, sh, SSE_Y_EXP, SSE_Y_MSB, SSE_H))+\
	quantize_signed(X, sh, SSE_X_EXP, SSE_X_MSB, SSE_W)

		QUANTIZE((NNE-NN)>>6,			(N-NW)>>4,		W-WW),//5
		QUANTIZE((NE-NNE)>>6,			(N-NN)>>4,		W-NW),//6
		QUANTIZE((int)eNW>>6,			(int)eW<<1,		(int)eN<<1),//9
		QUANTIZE((int)eNE>>6,			(int)eWW<<1,		(int)eNN<<1),//10
		QUANTIZE((3*(NE-NW)+NNWW-NNEE)>>9,	(W-WW)>>3,		(NW+2*NE+WW-4*W)>>2),//4
		QUANTIZE((int)(eNE-eNNE)>>6,		(int)(eN-eNN)>>4,	(int)(eW-eNW)),//8
		QUANTIZE((int)(eNNE-eNN)>>6,		(int)(eN-eNW)>>4,	(int)(eW-eWW)),//7
		QUANTIZE((NW+3*(NE-N)-NEE)>>9,		(NE-NNEE)>>3,		(NNE+NEE+2*W-4*NE)>>2),//1
		QUANTIZE((NWW+3*(N-NW)-NE)>>9,		(N-NN)>>3,		(W+NW+NE+NN-4*N)>>2),//2
		QUANTIZE((3*(N-W)+WW-NN)>>9,		(NW-NNWW)>>3,		(NNW+NWW+W+N-4*NW)>>2),//3
#undef  QUANTIZE
	};
#endif
	PROF(CALC_CTX);

	pr->sse_corr=0;
	for(int k=0;k<SSE_STAGES+1;++k)
	{
		long long sse_val;
		if(k<SSE_STAGES)
		{
			int qp=(int)(pr->pred>>(sh+8-(SSE_PREDBITS-1)));
			qp+=(1<<SSE_PREDBITS)>>1;
			qp=CLAMP(0, qp, (1<<SSE_PREDBITS)-1);
			pr->sse_idx[k]=(SSE_W*SSE_H*SSE_D)*qp+g[k];
			sse_val=pr->sse[SSE_SIZE*(SSE_STAGES*kc+k)+pr->sse_idx[k]];
			
#ifdef TRACK_SSE_RANGES
			int kx=g[k]%pr->sse_width, ky=g[k]/pr->sse_width%pr->sse_height, kz=g[k]/(pr->sse_width*pr->sse_height);
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
		else
		{
			int temp;
			int idx=THREEWAY(N, W)+1;
			temp=N+W-NW, idx=3*idx+THREEWAY(temp, (int)pr->pred)+1;
			idx=3*idx+THREEWAY(W, (int)pr->pred)+1;
			idx=3*idx+THREEWAY(NW, (int)pr->pred)+1;
			idx=3*idx+THREEWAY(N, (int)pr->pred)+1;
			idx=3*idx+THREEWAY(NE, (int)pr->pred)+1;
			idx=3*idx+THREEWAY(WW, (int)pr->pred)+1;
			idx=3*idx+THREEWAY(NN, (int)pr->pred)+1;
			temp=2*N-NN, idx=3*idx+THREEWAY(temp, (int)pr->pred)+1;
			temp=2*W-WW, idx=3*idx+THREEWAY(temp, (int)pr->pred)+1;
			if(idx>=SSE_FR_SIZE)
				LOG_ERROR("SSE OOB");
			pr->sse_idx[k]=idx;
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
			sse_val=pr->sse_fr[SSE_FR_SIZE*kc+pr->sse_idx[k]];
		}
		pr->sse_count[k]=(int)(sse_val&0xFFF);
		pr->sse_sum[k]=(int)(sse_val>>12);
		int sse_corr=(int)(pr->sse_sum[k]/(pr->sse_count[k]+1LL));
		pr->pred+=sse_corr;
		pr->sse_corr+=sse_corr;
	}
#endif

	int final_corr=(int)(pr->bias_sum[kc]/(pr->bias_count[kc]+1LL));
	pr->pred+=final_corr;
	pr->sse_corr+=final_corr;
	PROF(SSE_LOOP);

	pr->pred=CLAMP(clamp_lo, pr->pred, clamp_hi);
	pr->pred_final=(int)((pr->pred+((1<<PRED_PREC)>>1)-1)>>PRED_PREC);
	pr->kc=kc;
	pr->kx=kx;
}
static void slic5_update(SLIC5Ctx *pr, int curr, int token)
{
	PROF(PRED_TILL_UPDATE);
	int kc=pr->kc, kx=pr->kx;

	//update WP errors
	int error=(curr<<PRED_PREC)-(int)pr->pred;
	LOAD(pr->errors, 0, 0)=error;
	int errors[NPREDS]={0}, kbest=0;
	for(int k=0;k<NPREDS;++k)
	{
		errors[k]=abs((curr<<PRED_PREC)-pr->preds[k]);
		pr->pred_errors[(NPREDS*(pr->kym0+pr->kx+PAD_SIZE)+k)<<2|kc]=errors[k];
		if(pr->ky&&pr->kx+1<pr->iw)
			pr->pred_errors[(NPREDS*(pr->kym1+pr->kx+1+PAD_SIZE)+k)<<2|kc]+=errors[k];//eNE += ecurr
		if(errors[kbest]>errors[k])
			kbest=k;

		pr->pred_error_sums[k]+=errors[k];
	}

	//update WP weights
	++pr->params[kbest];
	if(pr->params[kbest]>(12<<pr->depths[kc]))
	{
		for(int k=0;k<NPREDS;++k)
			pr->params[k]>>=1;
	}
	PROF(UPDATE_WP);

	LOAD(pr->pixels, 0, 0)=curr;

	//update hist
	int total_hists=pr->nhist*pr->nhist;
	++pr->hist[pr->cdfsize*(total_hists*pr->kc+pr->hist_idx)+token];
	++pr->histsums[total_hists*pr->kc+pr->hist_idx];
	if(pr->histsums[pr->hist_idx]>(pr->iw*pr->nch<<1))//rescale hist
	{
		int *curr_hist=pr->hist+pr->cdfsize*(total_hists*pr->kc+pr->hist_idx);
		int sum=0;
		for(int ks=0;ks<pr->cdfsize;++ks)
		{
			curr_hist[ks]=(curr_hist[ks]+1)>>1;
			sum+=curr_hist[ks];
		}
		pr->histsums[pr->hist_idx]=sum;
	}
	PROF(UPDATE_HIST);
	
	//update SSE
#ifndef DISABLE_SSE
	for(int k=0;k<SSE_STAGES+1;++k)
	{
		++pr->sse_count[k];
		pr->sse_sum[k]+=error;
		if(pr->sse_count[k]>=640)
		{
			pr->sse_count[k]>>=1;
			pr->sse_sum[k]>>=1;
		}
		long long sse_val=pr->sse_sum[k]<<12|pr->sse_count[k];
		if(k<SSE_STAGES)
			pr->sse[SSE_SIZE*(SSE_STAGES*kc+k)+pr->sse_idx[k]]=sse_val;
		else
			pr->sse_fr[SSE_FR_SIZE*kc+pr->sse_idx[k]]=sse_val;
	}
#endif
	++pr->bias_count[kc];
	pr->bias_sum[kc]+=error;
	PROF(UPDATE_SSE);
}
#undef  LOAD
static void slic5_enc(SLIC5Ctx *pr, int curr, int kc, int kx, int ky)
{
	slic5_predict(pr, kc, kx, ky);

	int error=curr-pr->pred_final;
#if 1
	int upred=pr->half[kc]-abs(pr->pred_final);
	int abs_error=abs(error), negmask;
	if(abs_error<=upred)
	{
		negmask=-(pr->sse_corr<0);//sign is flipped if SSE correction was negative, to skew the histogram
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

	HybridUint hu;
	hybriduint_encode(error, &hu);
	if(hu.token>=pr->cdfsize)
		LOG_ERROR("Token OOB %d/%d", hu.token, pr->cdfsize);

	EC_ENC(pr->ec, hu.token, pr->CDFs+(pr->cdfsize+1)*(pr->nhist*pr->nhist*kc+pr->hist_idx), pr->cdfsize, 0);
	if(hu.nbits)
	{
		int bypass=hu.bypass, nbits=hu.nbits;
		while(nbits>8)
		{
			EC_ENC(pr->ec, bypass>>(nbits-8)&0xFF, 0, 1<<8, 16-8);
			nbits-=8;
		}
		EC_ENC(pr->ec, bypass&((1<<nbits)-1), 0, 1<<nbits, 16-nbits);
	}

	slic5_update(pr, curr, hu.token);
}
static int slic5_dec(SLIC5Ctx *pr, int kc, int kx, int ky)
{
	slic5_predict(pr, kc, kx, ky);

	int token=EC_DEC(pr->ec, pr->CDFs+(pr->cdfsize+1)*(pr->nhist*pr->nhist*kc+pr->hist_idx), pr->cdfsize, 0);
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
			bypass|=EC_DEC(pr->ec, 0, 1<<8, 16-8)<<n;
		}
		bypass|=EC_DEC(pr->ec, 0, 1<<n, 16-n);
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
		negmask=-(pr->sse_corr<0);
		error=error>>1^-(error&1);
	}
	else
	{
		negmask=-(pr->pred_final>0);
		error=error-upred;
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
}
typedef enum RCTTypeEnum
{
	RCT_NONE,
	RCT_SubGreen,
	RCT_JPEG2000,
	RCT_YCoCg_R,
	RCT_YCbCr_R_v1,
	RCT_YCbCr_R_v2,
	RCT_YCbCr_R_v3,
	RCT_YCbCr_R_v4,
	RCT_COUNT,
} RCTType;
static RCTType select_best_RCT(Image const *src)
{
	int inflation[]={1, 1, 1, 0};
	int maxdepth=calc_maxdepth(src, inflation), maxlevels=1<<maxdepth;
	int *hist=(int*)malloc(src->nch*maxlevels*sizeof(int[RCT_COUNT]));
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return RCT_NONE;
	}
	memset(hist, 0, src->nch*maxlevels*sizeof(int[RCT_COUNT]));
	int res=src->iw*src->ih;
	int nlevels[]=
	{
		1<<src->depth[0],	//0 r
		1<<src->depth[1],	//1 Y/g
		1<<src->depth[2],	//2 b
		1<<(src->depth[2]+1),	//3 Cb
		1<<(src->depth[0]+1),	//4 Cr
	};
	for(int k=0;k<res;++k)
	{
		int v[]={src->data[k<<2|0], src->data[k<<2|1], src->data[k<<2|2]}, r, g, b;
#define CLAMP_PIXEL(IDX, X, CH) maxlevels*IDX+((X+(nlevels[CH]>>1))&(nlevels[CH]-1))

		//RCT_NONE
		r=v[0], g=v[1], b=v[2];
		++hist[CLAMP_PIXEL( 0, g, 1)];
		++hist[CLAMP_PIXEL( 1, b, 2)];
		++hist[CLAMP_PIXEL( 2, r, 0)];
		
		//RCT_SubGreen
		r=v[0], g=v[1], b=v[2];
		r-=g;
		b-=g;
		++hist[CLAMP_PIXEL( 3, g, 1)];
		++hist[CLAMP_PIXEL( 4, b, 3)];
		++hist[CLAMP_PIXEL( 5, r, 4)];

		//RCT_JPEG2000
		r=v[0], g=v[1], b=v[2];
		r-=g;
		b-=g;
		g+=(r+b)>>1;
		++hist[CLAMP_PIXEL( 6, g, 1)];
		++hist[CLAMP_PIXEL( 7, b, 3)];
		++hist[CLAMP_PIXEL( 8, r, 4)];

		//RCT_YCoCg_R
		r=v[0], g=v[1], b=v[2];
		r-=b;
		b+=r>>1;
		g-=b;
		b+=g>>1;
		++hist[CLAMP_PIXEL( 9, g, 1)];
		++hist[CLAMP_PIXEL(10, b, 3)];
		++hist[CLAMP_PIXEL(11, r, 4)];

		//RCT_YCbCr_R_v1
		r=v[0], g=v[1], b=v[2];
		r-=g;
		g+=r>>1;
		b-=g;
		g+=b>>1;
		++hist[CLAMP_PIXEL(12, g, 1)];
		++hist[CLAMP_PIXEL(13, b, 3)];
		++hist[CLAMP_PIXEL(14, r, 4)];

		//RCT_YCbCr_R_v2
		r=v[0], g=v[1], b=v[2];
		r-=g;
		g+=r>>1;
		b-=g;
		g+=(2*b-r+4)>>3;
		++hist[CLAMP_PIXEL(15, g, 1)];
		++hist[CLAMP_PIXEL(16, b, 3)];
		++hist[CLAMP_PIXEL(17, r, 4)];

		//RCT_YCbCr_R_v3
		r=v[0], g=v[1], b=v[2];
		r-=g;
		g+=r>>1;
		b-=g;
		g+=(2*b+r+4)>>3;
		++hist[CLAMP_PIXEL(18, g, 1)];
		++hist[CLAMP_PIXEL(19, b, 3)];
		++hist[CLAMP_PIXEL(20, r, 4)];

		//RCT_YCbCr_R_v4
		r=v[0], g=v[1], b=v[2];
		r-=g;
		g+=r>>1;
		b-=g;
		g+=b/3;
		++hist[CLAMP_PIXEL(21, g, 1)];
		++hist[CLAMP_PIXEL(22, b, 3)];
		++hist[CLAMP_PIXEL(23, r, 4)];
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
					csizes[kc]+=freq*bitsize;
				}
			}
		}
		double csize=csizes[0]+csizes[1]+csizes[2];
		if(!kt||bestsize>csize)
			kbest=kt, bestsize=csize;
	}
	free(hist);
	return kbest;
}
static const Image *guide=0;
int t47_encode(Image const *src, ArrayHandle *data, int loud)
{
	guide=src;
	double t_start=time_sec();
	if(loud)
	{
		int maxdepth=calc_maxdepth(src, 0);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T47 SLIC5-AC  Enc %s  CWHD %d*%d*%d*%d/8\n", g_buf, src->nch, src->iw, src->ih, maxdepth);
	}
	RCTType best_rct=src->nch>=3?select_best_RCT(src):RCT_NONE;
	double t_selectRCT=time_sec()-t_start;
	DList list;
	ArithmeticCoder ec;
	SLIC5Ctx pr;
	dlist_init(&list, 1, 1024, 0);
	ac_enc_init(&ec, &list);
	int success=slic5_init(&pr, src->iw, src->ih, src->nch, src->depth, &ec);
	if(!success)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	if(src->nch>=3)//emit best RCT
		dlist_push_back1(&list, &best_rct);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		slic5_nextrow(&pr, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			//if(kx==2656&&ky==1118)//
			//	printf("");

			int kc=0;
			if(src->nch>=3)
			{
				int r=src->data[idx|0], g=src->data[idx|1], b=src->data[idx|2];
				switch(best_rct)
				{
				case RCT_SubGreen:
					r-=g;
					b-=g;
					break;
				case RCT_JPEG2000:
					r-=g;
					b-=g;
					g+=(r+b)>>2;
					break;
				case RCT_YCoCg_R:
					r-=b;
					b+=r>>1;
					g-=b;
					b+=g>>1;
					{
						int temp;
						SWAPVAR(g, b, temp);//move luma from C2 to C1
					}
					break;
				case RCT_YCbCr_R_v1:
					r-=g;
					g+=r>>1;
					b-=g;
					g+=b>>1;
					break;
				case RCT_YCbCr_R_v2:
					r-=g;
					g+=r>>1;
					b-=g;
					g+=(2*b-r+4)>>3;
					break;
				case RCT_YCbCr_R_v3:
					r-=g;
					g+=r>>1;
					b-=g;
					g+=(2*b+r+4)>>3;
					break;
				case RCT_YCbCr_R_v4:
					r-=g;
					g+=r>>1;
					b-=g;
					g+=b/3;
					break;
				default:
					break;
				}
				slic5_enc(&pr, g, 0, kx, ky);//Y
				slic5_enc(&pr, b, 1, kx, ky);//Cb
				slic5_enc(&pr, r, 2, kx, ky);//Cr
				kc+=3;
			}
			while(kc<src->nch)//gray/alpha
			{
				int pixel=src->data[idx|kc];
				slic5_enc(&pr, pixel, kc, kx, ky);
				++kc;
			}
		}
	}
	ac_enc_flush(&ec);
	dlist_appendtoarray(&list, data);
	if(loud)
	{
		double usize=image_getBMPsize(src);
		printf("\n");
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		if(src->nch>=3)
			printf(" (select RCT=%d  %lf sec)", best_rct, t_selectRCT);
		printf("\n");
		printf("csize %8d  CR %10.6lf\n", (int)list.nobj, usize/list.nobj);

		const char *ch_labels[]={"Y", "Cb", "Cr", "A"};
		int ch_nhist=pr.nhist*pr.nhist, total_hist=pr.nch*ch_nhist;
		printf("Hist CWH %d*%d*%d:\n", pr.nch, pr.cdfsize, ch_nhist);
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

		printf("Pred errors\n");
		long long sum=0;
		for(int k=0;k<NPREDS;++k)
			sum+=pr.pred_error_sums[k];
		for(int k=0;k<NPREDS;++k)
			printf("  %2d  %12lld  %6.2lf%%\n", k, pr.pred_error_sums[k], 100.*pr.pred_error_sums[k]/sum);
		
#ifdef TRACK_SSE_RANGES
		const char dim_labels[]="XYZP";
		printf("SSE ranges\n");
		for(int k=0;k<_countof(pr.sse_ranges)-1;k+=2)
		{
			int n=(&pr.sse_width)[k>>1];
			printf("  %c %4d ~ %4d / %4d  ", dim_labels[k>>1], pr.sse_ranges[k], pr.sse_ranges[k+1], n);
			for(int k2=0;k2<n;++k2)
				printf("%c", k2>=pr.sse_ranges[k]&&k2<=pr.sse_ranges[k+1]?'*':'-');
			printf("\n");
		}
#endif

		printf("Bias\n");
		for(int kc=0;kc<src->nch;++kc)
			printf("  %lld/%d = %d\n", pr.bias_sum[kc], pr.bias_count[kc], pr.bias_count[kc]?(int)(pr.bias_sum[kc]/pr.bias_count[kc]):0);

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
	return 1;
}
int t47_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	double t_start=time_sec();
	RCTType best_rct=RCT_NONE;
	if(dst->nch>=3)
	{
		if(srclen<1)
		{
			LOG_ERROR("Corrupt file");
			return 0;
		}
		best_rct=*data;
		if((unsigned)best_rct>=(unsigned)RCT_COUNT)
		{
			LOG_ERROR("Corrupt file");
			return 0;
		}
		++data;
		--srclen;
	}
	ArithmeticCoder ec;
	SLIC5Ctx pr;
	ac_dec_init(&ec, data, data+srclen);
	int success=slic5_init(&pr, dst->iw, dst->ih, dst->nch, dst->depth, &ec);
	if(!success)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int ky=0, idx=0;ky<dst->ih;++ky)
	{
		slic5_nextrow(&pr, ky);
		for(int kx=0;kx<dst->iw;++kx, idx+=4)
		{
			//if(kx==2656&&ky==1118)//
			//	printf("");

			int kc=0;
			if(dst->nch>=3)
			{
				int g=slic5_dec(&pr, 0, kx, ky);//Y
				int b=slic5_dec(&pr, 1, kx, ky);//Cb
				int r=slic5_dec(&pr, 2, kx, ky);//Cr
				switch(best_rct)
				{
				case RCT_SubGreen:
					b+=g;
					r+=g;
					break;
				case RCT_JPEG2000:
					g-=(r+b)>>2;
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
				case RCT_YCbCr_R_v2:
					g-=(2*b-r+4)>>3;
					b+=g;
					g-=r>>1;
					r+=g;
					break;
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
				default:
					break;
				}

				//if(guide&&(r!=guide->data[(guide->iw*ky+kx)<<2|0]||g!=guide->data[(guide->iw*ky+kx)<<2|1]||b!=guide->data[(guide->iw*ky+kx)<<2|2]))
				//	LOG_ERROR("Guide error XY %d %d", kx, ky);

				dst->data[idx|0]=r;
				dst->data[idx|1]=g;
				dst->data[idx|2]=b;

				kc+=3;
			}
			while(kc<dst->nch)//gray/alpha
			{
				int pixel=slic5_dec(&pr, kc, kx, ky);
				dst->data[idx|kc]=pixel;
				++kc;
			}
		}
	}
	if(loud)
	{
		printf("\n");
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	slic5_free(&pr);
	return 1;
}

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

			slic5_enc(&pr, g, 0, kx, ky);//Y
			slic5_enc(&pr, b, 1, kx, ky);//Cb
			slic5_enc(&pr, r, 2, kx, ky);//Cr
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
			int g=slic5_dec(&pr, 0, kx, ky);//Y
			int b=slic5_dec(&pr, 1, kx, ky);//Cb
			int r=slic5_dec(&pr, 2, kx, ky);//Cr

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