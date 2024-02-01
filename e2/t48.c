#include"e2.h"
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


//	#define ENABLE_GUIDE

//	#define BACKPROP
	#define ENABLE_PROB_BIAS
//	#define PRINT_AVERAGE_ERROR
//	#define TRACK_SSE_RANGES//SLOW
//	#define DISABLE_SSE//less efficient
//	#define PROBBITS_15//less efficient

//	#define PROFILER//SLOWS DOWN BATCH TEST
//	#define AVX512
	#define AVX2


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

#define LG_SEQ_CTR 4
#define LG_PAR_CTR 1
//#define LG_CTR (LG_SEQ_CTR+LG_PAR_CTR)
//#define NCOUNTERS (1<<LG_CTR)
#define PAD_SIZE 2
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
#define SSE_PRED_LEVELS (1<<SSE_PREDBITS)
#define SSE_FR_SIZE 59049//3^10		separate final round
#define SSE_STAGES 10
#define SSE_SIZE (SSE_W*SSE_H*SSE_D*SSE_PRED_LEVELS)
#define PRED_PREC 8
#define PARAM_PREC 8

#define NPREDS 14
#define PREDLIST\
	PRED(W+NE-N-((2*(eN+eW)+eNE-eNW+4)>>3))\
	PRED(N-(int)(((long long)eN+eW+eNE)*-0x05C>>PARAM_PREC))\
	PRED(W-(int)(((long long)eN+eW+eNW)*-0x05B>>PARAM_PREC))\
	PRED(N+(int)((-eNN*0x0DFLL-eN*0x0051LL-eNE*0x0BDLL+((long long)N-NN)*0x05C+((long long)NW-W)*0x0102)>>PARAM_PREC))\
	PRED((W+NEE)>>1)\
	PRED(N+NE-NNE-eNNE)\
	PRED((N+W)>>1)\
	PRED((4*(N+W+NW+NE)-(NN+WW+NNWW+NNEE)+6)/12)\
	PRED((W+2*NE-NNE+eNE)>>1)\
	PRED(N+W+NNWW-NNW-NWW+((eN+eW)>>1))\
	PRED(N+W-NW)\
	PRED(2*W-WW+eW-eWW)\
	PRED(paper_GAP)\
	PRED(calic_GAP)

static const char *prednames[]=
{
#define PRED(X) #X,
	PREDLIST
#undef  PRED
};
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
typedef struct ProbCellStruct
{
	unsigned short p0[1<<LG_PAR_CTR], weight[1<<LG_PAR_CTR];
	int sse;//signed 24-bit sum | unsigned 8-bit count
} ProbCell;
typedef struct SLIC6CtxStruct
{
	int iw, ih, nch;
	char depths[4];
	int
		half[4],
		shift_prec[4],//shift right to bring value to 8-bit, removes predictor-added bits
		shift_error[4];
	int *pred_errors, *errors, *pixels;

	int nhist, total_hist;
	ProbCell *counters;
	long long *prob_sse;
	//long long prob_sse[1<<LG_SEQ_CTR];
	//int nhist, cdfsize;
	//int *hist, *histsums;
	//CDF_t *CDFs;

#ifndef DISABLE_SSE
	long long *sse, *sse_fr;
	int sse_idx[SSE_STAGES+1];
	long long sse_sum[SSE_STAGES+1];
	int sse_count[SSE_STAGES+1];
#endif
	int bias_count[4];
	long long bias_sum[4];

	int preds[NPREDS], params[NPREDS];
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
} SLIC6Ctx;
#define LOAD(BUF, X, Y) BUF[(pr->kym##Y+kx+PAD_SIZE-(X))<<2|kc]
#define LOAD_PRED_ERROR(C, X, Y, P) pr->pred_errors[(NPREDS*(pr->kym##Y+kx+PAD_SIZE-(X))+P)<<2|C]
#ifdef PRINT_AVERAGE_ERROR
static double av_err[4]={0};
#endif
static int slic6_init(SLIC6Ctx *pr, int iw, int ih, int nch, const char *depths, ArithmeticCoder *ec)
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
	//HybridUint hu;
	//int extremesym=-((1<<maxdepth)>>1);
	//extremesym=extremesym<<1^-(extremesym<0);//pack sign
	//hybriduint_encode(extremesym, &hu);//encode -half
	//pr->cdfsize=hu.token+1;
	pr->nhist=QUANTIZE_HIST(767);
	pr->total_hist=pr->nhist*pr->nhist;

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

	pr->counters=(ProbCell*)malloc((size_t)nch*pr->total_hist*sizeof(ProbCell[1<<LG_SEQ_CTR]));
	pr->prob_sse=(long long*)malloc(sizeof(long long[256<<LG_SEQ_CTR]));
	//pr->counters=(unsigned short(*)[2])malloc((size_t)nch*pr->total_hist*sizeof(short[NCOUNTERS][2]));
	//pr->prob_sse=(long long*)malloc((size_t)nch*pr->total_hist*sizeof(long long));
	//pr->hist=(int*)malloc((size_t)nch*total_hist*pr->cdfsize*sizeof(int));//WH: cdfsize * NHIST
	//pr->histsums=(int*)malloc((size_t)nch*total_hist*sizeof(int));
	//pr->CDFs=(CDF_t*)malloc((size_t)nch*total_hist*(pr->cdfsize+1LL)*sizeof(CDF_t));
#ifndef DISABLE_SSE
	pr->sse=(long long*)malloc(nch*sizeof(long long[SSE_SIZE*SSE_STAGES]));
	pr->sse_fr=(long long*)malloc(nch*sizeof(long long[SSE_FR_SIZE]));
#endif
	if(!pr->pred_errors||!pr->errors||!pr->pixels||!pr->counters||!pr->prob_sse
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

	ProbCell fillval=
	{
		{0x8000, 0x8000},//p0's
		{0x8000, 0x8000},//weights
		0,//SSE
	};
	memfill(pr->counters, &fillval, (size_t)nch*pr->total_hist*sizeof(ProbCell[1<<LG_SEQ_CTR]), sizeof(ProbCell));
	//memfill(pr->counters, &fillval, (size_t)nch*pr->total_hist*sizeof(short[NCOUNTERS][2]), sizeof(short));
	memset(pr->prob_sse, 0, sizeof(long long[256<<LG_SEQ_CTR]));
#if 0
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

#ifndef DISABLE_SSE
	memset(pr->sse, 0, nch*sizeof(long long[SSE_SIZE*SSE_STAGES]));
	memset(pr->sse_fr, 0, nch*sizeof(long long[SSE_FR_SIZE]));
#endif

	for(int k=0;k<NPREDS;++k)
		pr->params[k]=12<<maxdepth;
	//int shift=16-maxdepth;
	//shift=MAXVAR(0, shift);
	//memcpy(pr->params, wp_params, sizeof(pr->params));
	//for(int k=0;k<NPREDS;++k)
	//	pr->params[k]=pr->params[k]>>shift;

	pr->ec=ec;
	return 1;
}
static void slic6_free(SLIC6Ctx *pr)
{
	free(pr->pred_errors);
	free(pr->errors);
	free(pr->pixels);
	free(pr->counters);
	//free(pr->hist);
	//free(pr->histsums);
	//free(pr->CDFs);
#ifndef DISABLE_SSE
	free(pr->sse);
	free(pr->sse_fr);
#endif
}
static void slic6_nextrow(SLIC6Ctx *pr, int ky)
{
	pr->ky=ky;
	pr->kym0=(pr->iw+PAD_SIZE*2)*(ky&3);
	pr->kym1=(pr->iw+PAD_SIZE*2)*((ky-1)&3);
	pr->kym2=(pr->iw+PAD_SIZE*2)*((ky-2)&3);
}
static void slic6_predict(SLIC6Ctx *pr, int kc, int kx)
{
	PROF(OUTSIDE);
	int idx=(pr->iw*pr->ky+kx)<<2|kc;
	pr->kc=kc;
	pr->kx=kx;
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

	int j=-1;
#define PRED(X) pr->preds[++j]=X;
	PREDLIST
#undef  PRED
	PROF(CALC_SUBPREDS);

	long long weights[NPREDS]={0}, wsum=0;
	for(int k=0;k<NPREDS;++k)
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

		//	(long long)pr->pred_errors[(NPREDS*(pr->kym2+kx+PAD_SIZE+2)+k)<<2|kc]+//peNNEE
		//	(long long)pr->pred_errors[(NPREDS*(pr->kym2+kx+PAD_SIZE+1)+k)<<2|kc]+//peNNE
		//	(long long)pr->pred_errors[(NPREDS*(pr->kym2+kx+PAD_SIZE  )+k)<<2|kc]+//peNN+peNW
		//	(long long)pr->pred_errors[(NPREDS*(pr->kym2+kx+PAD_SIZE-1)+k)<<2|kc]+//peNNW+peNWW
		//	(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx+PAD_SIZE+2)+k)<<2|kc]+//peNEE
		//	(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx+PAD_SIZE+1)+k)<<2|kc]+//peNE
		//	(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx+PAD_SIZE  )+k)<<2|kc]*3+//peN+peW
		//	(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx+PAD_SIZE-1)+k)<<2|kc];//peNW+peWW
		weights[k]=((long long)pr->params[k]<<(pr->depths[kc]+PRED_PREC+5))/(weights[k]+1);
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

	int qx=QUANTIZE_HIST((dx+abs(eW+eWW))>>(sh-2)), qy=QUANTIZE_HIST((dy+abs(eN+eNN))>>(sh-2));
	qx=CLAMP(0, qx, pr->nhist-1);
	qy=CLAMP(0, qy, pr->nhist-1);
	pr->hist_idx=pr->nhist*qy+qx;

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
	ALIGN(32) int g2[16];
	_mm256_store_si256((__m256i*)g2+0, X0);
	_mm256_store_si256((__m256i*)g2+1, X1);
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
			pr->sse_idx[k]=(SSE_W*SSE_H*SSE_D)*qp+g2[k];
			sse_val=pr->sse[SSE_SIZE*(SSE_STAGES*kc+k)+pr->sse_idx[k]];
			
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
}
static void slic6_update(SLIC6Ctx *pr, int curr)
{
	PROF(PRED_TILL_UPDATE);
	int kc=pr->kc, kx=pr->kx;

	LOAD(pr->pixels, 0, 0)=curr;

	//update WP errors
	int error=(curr<<PRED_PREC)-(int)pr->pred;
	LOAD(pr->errors, 0, 0)=error;
	int errors[NPREDS]={0}, kbest=0;
	for(int k=0;k<NPREDS;++k)
	{
		errors[k]=abs((curr<<PRED_PREC)-pr->preds[k]);
		LOAD_PRED_ERROR(kc, 0, 0, k)=errors[k];

		//if(LOAD_PRED_ERROR(kc, 0, 0, k)<0)//
		//	LOG_ERROR("");

		//pr->pred_errors[(NPREDS*(pr->kym0+pr->kx+PAD_SIZE)+k)<<2|kc]=errors[k];
		if(pr->ky&&pr->kx+1<pr->iw)
		{
			LOAD_PRED_ERROR(kc, -1, 1, k)+=errors[k];//eNE += ecurr
		
			//if(LOAD_PRED_ERROR(kc, -1, 1, k)<0)//
			//	LOG_ERROR("");
		}
			//pr->pred_errors[(NPREDS*(pr->kym1+pr->kx+PAD_SIZE+1)+k)<<2|kc]+=errors[k];
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
static int p0_calc(ProbCell const *cell, const long long *prob_sse, int *ret_wsum, int *ret_p0)
{
	unsigned wsum=0;
	long long p0_final=0;
	int sse_count, sse_sum;
	for(int k=0;k<(1<<LG_PAR_CTR);++k)
	{
		p0_final+=(long long)cell->p0[k]*cell->weight[k];
		wsum+=cell->weight[k];
	}
	p0_final/=(unsigned long long)wsum+1;
	*ret_wsum=wsum;
	*ret_p0=(int)p0_final;

	sse_count=cell->sse&0xFF;
	sse_sum=cell->sse>>8;
	//p0_final=(p0_final*sse_count + sse_sum*wsum)/(wsum*sse_count+1);//X
	p0_final+=sse_sum/(sse_count+1);
#ifdef ENABLE_PROB_BIAS
	p0_final=CLAMP(1, p0_final, 0xFFFF);
	long long sseval=prob_sse[p0_final>>8];
	sse_count=sseval&0xFF;
	sse_sum=(int)(sseval>>8);
	p0_final+=sse_sum/(sse_count+1);
#endif
	p0_final=CLAMP(1, p0_final, 0xFFFF);
	return (int)p0_final;

#if 0
	int sse_count, sse_sum, sse_corr;
	int p0;
	unsigned wsum=0;
	unsigned long long p0_final=0;
	for(int k=0;k<(1<<LG_PAR_CTR);++k)
	{
		sse_count=cell[k].sse&0xFF;
		sse_sum=cell[k].sse>>8;
		sse_corr=sse_sum/(sse_count+1);
		p0=cell[k].p0+sse_corr;
		p0_final+=(long long)p0*cell[k].weight;
		wsum+=cell[k].weight;
	}
	p0_final/=(unsigned long long)wsum+1;
	p0_final=CLAMP(1, p0_final, 0xFFFF);
	return (int)p0_final;
#endif

#if 0
	unsigned long long p0_final=0;
	unsigned wsum=0;
	p0_final+=(unsigned long long)p0[0][0]*p0[0][1], wsum+=p0[0][1];
	p0_final+=(unsigned long long)p0[1][0]*p0[1][1], wsum+=p0[1][1];
	//p0_final+=(unsigned long long)p0[2][0]*p0[2][1], wsum+=p0[2][1];
	//p0_final+=(unsigned long long)p0[3][0]*p0[3][1], wsum+=p0[3][1];
	//p0_final+=(unsigned long long)p0[4][0]*p0[0][1], wsum+=p0[4][1];
	//p0_final+=(unsigned long long)p0[5][0]*p0[1][1], wsum+=p0[5][1];
	//p0_final+=(unsigned long long)p0[6][0]*p0[2][1], wsum+=p0[6][1];
	//p0_final+=(unsigned long long)p0[7][0]*p0[3][1], wsum+=p0[7][1];
	p0_final/=(unsigned long long)wsum+1;
	return (int)p0_final;
#endif
}
static void p0_update(ProbCell *cell, long long *prob_sse, int wsum, int p0, int p0_final, int bit)
{
	//small vs large change:  shift right ammount is "resistance", more slows it down, less speeds it up
	static const int shift_p0[]=
	{
		6, 7
	};
	static const int shift_weight[]=
	{
		5, 10
	};
	unsigned sse_count;
	int sse_sum;
	int p0_perf=!bit<<16;
	int error=p0_perf-p0_final;

	sse_count=cell->sse&0xFF;
	sse_sum=cell->sse>>8;
	++sse_count;
	sse_sum+=error/3;
	if(sse_count>0xFF)
	{
		sse_count>>=1;
		sse_sum>>=1;
	}
	cell->sse=sse_sum<<8|sse_count;
	
	long long sseval=prob_sse[p0_final>>8];
	sse_count=sseval&0xFF;
	sse_sum=(int)(sseval>>8);
	++sse_count;
	sse_sum+=error>>5;
	if(sse_count>0xFF)
	{
		sse_count>>=1;
		sse_sum>>=1;
	}
	prob_sse[p0_final>>8]=(long long)sse_sum<<8|sse_count;

	int negmask=-(error<0);
	for(int k=0;k<(1<<LG_PAR_CTR);++k)
	{
#ifdef BACKPROP
		int update=(((long long)cell->p0[k]-p0)<<16)/wsum;
		update^=negmask;
		update-=negmask;
		update>>=5;//LR = 0.03125
		update+=cell->weight[k];
		update=CLAMP(1, update, 0xFFFF);
		cell->weight[k]=update;
#else
		cell->weight[k]+=(((!bit==(cell->p0[k]>0x8000))<<16)-cell->weight[k]+((1<<shift_weight[k])>>1))>>shift_weight[k];
#endif
		cell->p0[k]+=((!bit<<16)-cell->p0[k]+((1<<shift_p0[k])>>1))>>shift_p0[k];
	}
#if 0
	const int shift_p0[]=
	{
		6, 7
	};
	const int shift_weight[]=
	{
		5, 10
	};
	unsigned sse_count;
	int sse_sum;
	int error=(!bit<<16)-p0_final;
	for(int k=0;k<(1<<LG_PAR_CTR);++k)
	{
		sse_count=cell[k].sse&0xFF;
		sse_sum=cell[k].sse>>8;
		++sse_count;
		sse_sum+=error;
		if(sse_count>0xFF)
		{
			sse_count>>=1;
			sse_sum>>=1;
		}
		cell[k].sse=sse_sum<<8|sse_count;

		cell[k].weight+=(((!bit==(cell[k].p0>0x8000))<<16)-cell[k].weight+((1<<shift_weight[k])>>1))>>shift_weight[k];

		cell[k].p0+=((!bit<<16)-cell[k].p0+((1<<shift_p0[k])>>1))>>shift_p0[k];
	}
#endif
#if 0
	p0[0][1]+=(((!bit==(p0[0][0]>0x8000))<<16)-p0[0][1]+16)>>5;
	p0[1][1]+=(((!bit==(p0[1][0]>0x8000))<<16)-p0[1][1]+16)>>10;
	//p0[2][1]+=(((!bit==(p0[2][0]>0x8000))<<16)-p0[2][1]+16)>>9;
	//p0[3][1]+=(((!bit==(p0[3][0]>0x8000))<<16)-p0[3][1]+16)>>10;
	//p0[0][1]+=(((!bit==(p0[0][0]>0x8000))<<16)-p0[0][1]+16)>>10;
	//p0[1][1]+=(((!bit==(p0[1][0]>0x8000))<<16)-p0[1][1]+16)>>10;
	//p0[2][1]+=(((!bit==(p0[2][0]>0x8000))<<16)-p0[2][1]+16)>>10;
	//p0[3][1]+=(((!bit==(p0[3][0]>0x8000))<<16)-p0[3][1]+16)>>10;

	p0[0][0]+=((!bit<<16)-p0[0][0]+16)>>6;
	p0[1][0]+=((!bit<<16)-p0[1][0]+16)>>7;
	//p0[2][0]+=((!bit<<16)-p0[2][0]+16)>>8;
	//p0[3][0]+=((!bit<<16)-p0[3][0]+16)>>9;
	//p0[0][0]+=((!bit<<16)-p0[0][0]+16)>>5;
	//p0[1][0]+=((!bit<<16)-p0[1][0]+16)>>6;
	//p0[2][0]+=((!bit<<16)-p0[2][0]+16)>>7;
	//p0[3][0]+=((!bit<<16)-p0[3][0]+16)>>8;
#endif
}
static void slic6_enc(SLIC6Ctx *pr, int curr, int kc, int kx)
{
	slic6_predict(pr, kc, kx);

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
	
#ifdef PRINT_AVERAGE_ERROR
	av_err[kc]+=error;
#endif
	int ctr_idx=0, p0=0, wsum=0;
	ProbCell *cells=pr->counters+(((size_t)pr->total_hist*kc+pr->hist_idx)<<LG_SEQ_CTR);
	//unsigned short *p0=pr->counters;
	for(int k=0;k<error;++k)
	{
		const int bit=1;
		int p0_final=p0_calc(cells+ctr_idx, pr->prob_sse+((size_t)ctr_idx<<8), &p0, &wsum);

		ac_enc_bin(pr->ec, (unsigned short)p0_final, bit);

		p0_update(cells+ctr_idx, pr->prob_sse+((size_t)ctr_idx<<8), p0, wsum, p0_final, bit);
		ctr_idx+=ctr_idx<(1<<LG_SEQ_CTR)-1;

		//if(ctr_idx==30)//
		//	printf("");
	}
	const int bit=0;
	int p0_final=p0_calc(cells+ctr_idx, pr->prob_sse+((size_t)ctr_idx<<8), &p0, &wsum);
	ac_enc_bin(pr->ec, p0_final, bit);
	p0_update(cells+ctr_idx, pr->prob_sse+((size_t)ctr_idx<<8), p0, wsum, p0_final, bit);
	//p0[ctr_idx]+=((!bit<<16)-p0[ctr_idx]+16)>>5;

	slic6_update(pr, curr);
}
static int slic6_dec(SLIC6Ctx *pr, int kc, int kx)
{
	slic6_predict(pr, kc, kx);
	
	int ctr_idx=0;
	ProbCell *cells=pr->counters+(((size_t)pr->total_hist*kc+pr->hist_idx)<<LG_SEQ_CTR);
	//ProbCell *cells=pr->counters+(((size_t)pr->total_hist*kc+pr->hist_idx)<<LG_CTR);
	//unsigned short (*p0)[2]=pr->counters+(((size_t)pr->total_hist*kc+pr->hist_idx)<<LG_CTR);
	//unsigned short *p0=pr->counters+NCOUNTERS*(pr->total_hist*kc+pr->hist_idx);
	//unsigned short *p0=pr->counters;
	int bit=0, error=0, p0=0, wsum=0;
	do
	{
		int p0_final=p0_calc(cells+ctr_idx, pr->prob_sse+((size_t)ctr_idx<<8), &p0, &wsum);
		//int p0_final=p0_calc(cells+ctr_idx, pr->prob_sse[ctr_idx]);

		bit=ac_dec_bin(pr->ec, p0_final);
		error+=bit;
		
		p0_update(cells+ctr_idx, pr->prob_sse+((size_t)ctr_idx<<8), p0, wsum, p0_final, bit);
		//p0_update(cells+ctr_idx, pr->prob_sse+ctr_idx, p0_final, bit);
		ctr_idx+=ctr_idx<(1<<LG_SEQ_CTR)-1;
		//ctr_idx+=(ctr_idx<(((1<<LG_SEQ_CTR)-1)<<LG_PAR_CTR))<<LG_PAR_CTR;

		//p0[ctr_idx]+=((!bit<<16)-p0[ctr_idx]+16)>>5;
		//ctr_idx+=ctr_idx<NCOUNTERS-1;
	}while(bit);

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
	slic6_update(pr, curr);
	return curr;
}
#define r p[0]
#define g p[1]
#define b p[2]
static void rct_fwd(int *p, RCTType rct)
{
	switch(rct)
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
	case RCT_YCbCr_R_v7:
		r-=g;
		g+=r>>1;
		b-=g;
		g+=(10*b-r+16)>>5;
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
	case RCT_YCbCr_R_v7:
		g-=(10*b-r+16)>>5;
		b+=g;
		g-=r>>1;
		r+=g;
		break;
	default:
		break;
	}
}
#undef	r
#undef	g
#undef	b
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
int t48_encode(Image const *src, ArrayHandle *data, int loud)
{
#ifdef ENABLE_GUIDE
	guide=src;
#endif
#ifdef PRINT_AVERAGE_ERROR
	memset(av_err, 0, sizeof(av_err));
#endif
	double t_start=time_sec();
	int nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);//TODO
	if(nch>src->nch)
		nch=src->nch;
	if(loud)
	{
		int maxdepth=calc_maxdepth(src, 0);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T48 CABAC  Enc %s  CWHD %d*%d*%d*%d/8\n", g_buf, nch, src->iw, src->ih, maxdepth);
	}
	RCTType rct=RCT_NONE;
	if(nch>=3)
	{
		double csizes[RCT_COUNT];
		double t1=time_sec();
		if(loud)
			printf("Selecting best RCT...\n");
		rct=rct_select_best(src, csizes);
		if(loud)
		{
			for(int kt=0;kt<RCT_COUNT;++kt)
			{
				printf("%20lf %s", csizes[kt], rct_names[kt]);
				if(kt==rct)
					printf(" <- %lf sec", time_sec()-t1);
				printf("\n");
			}
		}
	}
	DList list;
	ArithmeticCoder ec;
	SLIC6Ctx pr;
	dlist_init(&list, 1, 1024, 0);
	ac_enc_init(&ec, &list);
	int success=slic6_init(&pr, src->iw, src->ih, nch, src->depth, &ec);
	if(!success)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	if(nch>=3)//emit best RCT
		dlist_push_back1(&list, &rct);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		slic6_nextrow(&pr, ky);
		for(int kx=0;kx<src->iw;)
		{
			if(kx==158&&ky==571)//
				printf("");

			int kc=0;
			if(nch>=3)
			{
				int comp[]={src->data[idx|0], src->data[idx|1], src->data[idx|2]};
				rct_fwd(comp, rct);
				slic6_enc(&pr, comp[1], 0, kx);//Y		luma is encoded first
				slic6_enc(&pr, comp[2], 1, kx);//Cb
				slic6_enc(&pr, comp[0], 2, kx);//Cr
				kc+=3;
			}
			for(;kc<nch;++kc)//gray/alpha
			{
				int pixel=src->data[idx|kc];
				slic6_enc(&pr, pixel, kc, kx);
			}
			++kx;
			idx+=4;
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
		printf("\n");

		printf("csize %8d  CR %10.6lf\n", (int)list.nobj, usize/list.nobj);

		printf("Pred errors\n");
		long long sum=0;
		for(int k=0;k<NPREDS;++k)
			sum+=pr.pred_error_sums[k];
		for(int k=0;k<NPREDS;++k)
			printf("  %2d  %17lld  %6.2lf%%  %s\n", k, pr.pred_error_sums[k], 100.*pr.pred_error_sums[k]/sum, prednames[k]);
		
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
		//printf("SSE bias:\n");
		//for(int k=0;k<(1<<LG_SEQ_CTR);++k)
		//{
		//	int sse_count=pr.prob_sse[k]&0xFF, sse_sum=(int)(pr.prob_sse[k]>>8);
		//	printf("  %2d:  %d/%d = %d\n", k, sse_sum, sse_count, sse_sum/(sse_count+1));
		//}
#ifdef PRINT_AVERAGE_ERROR
		printf("Average error\n");
		for(int kc=0;kc<nch;++kc)
			printf("  C%d  %16lf\n", kc, av_err[kc]/((double)pr.nch*pr.iw*pr.ih));
#endif

		//printf("Bias\n");
		//for(int kc=0;kc<nch;++kc)
		//	printf("  %lld/%d = %d\n", pr.bias_sum[kc], pr.bias_count[kc], pr.bias_count[kc]?(int)(pr.bias_sum[kc]/pr.bias_count[kc]):0);

		prof_print();
	}
	slic6_free(&pr);
	dlist_clear(&list);
	return 1;
}
int t48_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	double t_start=time_sec();
	RCTType rct=RCT_NONE;
	int dst_nch=(dst->depth[0]!=0)+(dst->depth[1]!=0)+(dst->depth[2]!=0)+(dst->depth[3]!=0);
	int nch=MINVAR(dst_nch, dst->nch);
	if(nch>=3)
	{
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
	}
	ArithmeticCoder ec;
	SLIC6Ctx pr;
	ac_dec_init(&ec, data, data+srclen);
	int success=slic6_init(&pr, dst->iw, dst->ih, nch, dst->depth, &ec);
	if(!success)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int ky=0, idx=0;ky<dst->ih;++ky)
	{
		slic6_nextrow(&pr, ky);
		for(int kx=0;kx<dst->iw;)
		{
			//if(kx==252&&ky==16)//
			//if(kx==151&&ky==37)//
			//if(kx==100&&ky==38)//
			//if(kx==427&&ky==698)//
			if(kx==158&&ky==571)//
				printf("");

			int kc=0;
			if(nch>=3)
			{
				int comp[3];
				comp[1]=slic6_dec(&pr, 0, kx);//Y
				comp[2]=slic6_dec(&pr, 1, kx);//Cb
				comp[0]=slic6_dec(&pr, 2, kx);//Cr
				rct_inv(comp, rct);
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
				int pixel=slic6_dec(&pr, kc, kx);
				dst->data[idx|kc]=pixel;
				++kc;
			}
			++kx;
			idx+=4;
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
	}
	slic6_free(&pr);
	return 1;
}