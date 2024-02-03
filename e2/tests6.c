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

//debug
//	#define ENABLE_GUIDE

//experiments
//	#define PRINT_HIST
//	#define PROFILER//SLOW
//	#define TRACK_SSE_RANGES//SLOW

//efficiency
	#define ENABLE_SSE_4D
//	#define SSE_UPDATE_NB_CELLS//inferior
//	#define COMPENSATE_GAMMA//crude method, inferior
//	#define ENABLE_LZ//currently inferior
//	#define ENABLE_SSE_PAR//inferior

//memory
//	#define PROBBITS_15//less efficient

//speed
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

#define SLIC5_CONFIG_EXP 5
#define SLIC5_CONFIG_MSB 2
#define SLIC5_CONFIG_LSB 0

#define CDF_UPDATE_PERIOD 0x400//number of sub-pixels processed between CDF updates, must be a power-of-two
#define PAD_SIZE 2
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
#define SSE_W 7		//ENABLE SSE INDEX CLAMP WHAN CHANGING CONFIG
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
#define PRED_PREC 8
#define PARAM_PREC 8

#define NPREDS 14
#define PREDLIST\
	PRED(W+NE-N-((2*(eN+eW)+eNE-eNW+4)>>3))\
	PRED(N-(int)(((long long)eN+eW+eNE)*-0x05C>>PARAM_PREC))\
	PRED(W-(int)(((long long)eN+eW+eNW)*-0x05B>>PARAM_PREC))\
	PRED(N+(int)((-eNN*0x0DFLL-eN*0x0051LL-eNE*0x0BDLL+((long long)N-NN)*0x05C+((long long)NW-W)*0x0102)>>PARAM_PREC))\
	PRED((W+NEE)>>1)\
	PRED((4*N-2*NN+NW+NE)>>2)\
	PRED(N+NE-NNE-eNNE)\
	PRED((N+W)>>1)\
	PRED((4*(N+W+NW+NE)-(NN+WW+NNWW+NNEE)+6)/12)\
	PRED(N+W+NNWW-NNW-NWW+((eN+eW)>>1))\
	PRED(N+W-NW)\
	PRED(2*W-WW+eW-eWW)\
	PRED(paper_GAP)\
	PRED(calic_GAP)
//	PRED((WW+NE)>>1)
//	PRED((W+2*NE-NNE+eNE)>>1)

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
	token=CLAMP(0, token, nlevels-1);
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
typedef struct SLIC5CtxStruct
{
	int iw, ih, nch;
	char depths[4];
	int
		maxdepth, maxlevels,
		nlevels[4],
		half[4],
		shift_prec[4],//shift right to bring value to 8-bit, removes predictor-added bits
		shift_error[4],
		full_pixels_processed;
	int nhist, cdfsize;
	int *pred_errors, *errors, *pixels;

	int *hist, *histsums;
	CDF_t *CDFs;

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
	int bias_count[4];
	long long bias_sum[4];

#ifdef COMPENSATE_GAMMA
	int *eq_hist, *equalizer;
#endif

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
} SLIC5Ctx;
#define LOAD(BUF, X, Y) BUF[(pr->kym##Y+kx+PAD_SIZE-(X))<<2|kc]
#define LOAD_CH(BUF, C, X, Y) BUF[(pr->kym##Y+kx+PAD_SIZE-(X))<<2|C]
#define LOAD_PRED_ERROR(C, X, Y, P) pr->pred_errors[(NPREDS*(pr->kym##Y+kx+PAD_SIZE-(X))+P)<<2|C]
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
	pr->maxdepth=0;
	for(int kc=0;kc<nch;++kc)
	{
		int prec_half=(1<<(pr->depths[kc]+PRED_PREC))>>1;
		pr->nlevels[kc]=1<<pr->depths[kc];
		pr->half[kc]=pr->nlevels[kc]>>1;
		UPDATE_MAX(pr->maxdepth, pr->depths[kc]);
	}
	pr->maxlevels=1<<pr->maxdepth;
	for(int kc=0;kc<pr->nch;++kc)
	{
		pr->shift_prec[kc]=MAXVAR(8, pr->depths[kc])-8+PRED_PREC;
		pr->shift_error[kc]=1+((pr->depths[kc]<=9)<<1);
	}
	HybridUint hu;
	int extremesym=-(pr->maxlevels>>1);
	extremesym=extremesym<<1^-(extremesym<0);//pack sign
	hybriduint_encode(extremesym, &hu);//encode -half
	pr->cdfsize=hu.token+1;
#ifdef ENABLE_LZ
	++pr->cdfsize;//make way for LZ escape symbol
#endif
	pr->nhist=QUANTIZE_HIST(767);
	
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

	pr->pred_errors=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[NPREDS*4*4]));//NPREDS * 4 rows * 4 comps
	pr->errors=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	pr->pixels=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	int total_hist=pr->nhist*pr->nhist;
	pr->hist=(int*)malloc((size_t)nch*total_hist*pr->cdfsize*sizeof(int));//WH: cdfsize * NHIST
	pr->histsums=(int*)malloc((size_t)nch*total_hist*sizeof(int));
	pr->CDFs=(CDF_t*)malloc((size_t)nch*total_hist*(pr->cdfsize+1LL)*sizeof(CDF_t));
#ifdef ENABLE_SSE_PAR
	pr->sse=(long long*)malloc((size_t)nch*pr->sse_size*sizeof(long long[SSE_STAGES]));
	pr->sse_fr=(long long*)malloc(nch*sizeof(long long[SSE_FR_SIZE]));
#elif defined ENABLE_SSE_4D
	pr->sse=(long long*)malloc(nch*sizeof(long long[SSE_SIZE*SSE_STAGES]));
	pr->sse_fr=(long long*)malloc(nch*sizeof(long long[SSE_FR_SIZE]));
	if(pr->nch>1)
	{
		pr->sse_cfl=(long long*)malloc((nch-1LL)*sizeof(long long[SSE_SIZE<<1]));
		if(!pr->sse_cfl)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
	}
#endif
	if(!pr->pred_errors||!pr->errors||!pr->pixels||!pr->hist||!pr->histsums||!pr->CDFs
#if defined ENABLE_SSE_4D || defined ENABLE_SSE_PAR
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
	
#ifdef ENABLE_SSE_PAR
	memset(pr->sse, 0, (size_t)nch*pr->sse_size*sizeof(long long[SSE_STAGES]));
	memset(pr->sse_fr, 0, nch*sizeof(long long[SSE_FR_SIZE]));
#elif defined ENABLE_SSE_4D
	memset(pr->sse, 0, nch*sizeof(long long[SSE_SIZE*SSE_STAGES]));
	memset(pr->sse_fr, 0, nch*sizeof(long long[SSE_FR_SIZE]));
	if(pr->nch>1)
		memset(pr->sse_cfl, 0, (nch-1LL)*sizeof(long long[SSE_SIZE<<1]));
#endif
#ifdef COMPENSATE_GAMMA
	pr->eq_hist=(int*)malloc((size_t)nch*sizeof(int)<<pr->maxdepth);//1MB for 16-bit RGBA
	pr->equalizer=(int*)malloc(nch*((size_t)pr->maxlevels+1)*sizeof(int));
	if(!pr->eq_hist||!pr->equalizer)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(pr->eq_hist, 0, (size_t)nch*sizeof(int)<<pr->maxdepth);
	for(int kc=0;kc<pr->nch;++kc)
	{
		int *eq=pr->equalizer+((size_t)pr->maxlevels+1)*kc;
		for(int ks=0;ks<=pr->nlevels[kc];++ks)
			eq[ks]=ks<<PRED_PREC;
	}
#if 0
	int res=iw*ih;
	for(int kc=0;kc<nch;++kc)
	{
		int *eq=pr->equalizer+((size_t)pr->maxlevels+1)*kc;
		for(int k=0;k<res;++k)
		{
			int val=preview[k<<2|kc]+pr->half[kc];
			val=CLAMP(0, val, pr->nlevels[kc]-1);
			++eq[val];
		}
		int sum=0;
		for(int ks=0;ks<pr->nlevels[kc];++ks)
		{
			int freq=eq[ks];
			eq[ks]=(int)(((long long)sum<<(pr->depths[kc]+PRED_PREC))/res);
			sum+=freq;
		}
		//for(int ks=0;ks<pr->nlevels[kc];++ks)
		//{
		//	int curr=CDF_fwd[ks]>>PRED_PREC, next=ks<pr->nlevels[kc]-1?CDF_fwd[ks+1]>>PRED_PREC:pr->nlevels[kc];
		//	//if(curr<next)
		//	//	memfill(CDF_inv+curr, &ks, ((size_t)next-curr)*sizeof(int), sizeof(ks));
		//	for(int kl=curr;kl<next;++kl)
		//		CDF_inv[kl]=ks;
		//}
	}
#endif
#endif

	for(int k=0;k<NPREDS;++k)
		pr->params[k]=12<<pr->maxdepth;
	//int shift=16-maxdepth;
	//shift=MAXVAR(0, shift);
	//memcpy(pr->params, wp_params, sizeof(pr->params));
	//for(int k=0;k<NPREDS;++k)
	//	pr->params[k]=pr->params[k]>>shift;

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
#if defined ENABLE_SSE_4D || defined ENABLE_SSE_PAR
	free(pr->sse);
	free(pr->sse_fr);
	if(pr->nch>1)
		free(pr->sse_cfl);
#endif
	//memset(pr, 0, sizeof(*pr));//init does this already
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
#ifdef COMPENSATE_GAMMA
	if(pr->full_pixels_processed)
	{
		for(int kc=0;kc<pr->nch;++kc)
		{
			int *hist=pr->eq_hist+((size_t)kc<<pr->maxdepth);
			int *eq=pr->equalizer+((size_t)pr->maxlevels+1)*kc;
			int sum=0;
			for(int ks=0;ks<pr->nlevels[kc];++ks)
			{
				int freq=hist[ks];
				eq[ks]=(int)(((long long)sum<<(pr->depths[kc]+PRED_PREC))/pr->full_pixels_processed);
				sum+=freq;
			}
		}
	}
#endif
}
static void slic5_nextrow(SLIC5Ctx *pr, int ky)
{
	pr->ky=ky;
	pr->kym0=(pr->iw+PAD_SIZE*2)*(ky&3);
	pr->kym1=(pr->iw+PAD_SIZE*2)*((ky-1)&3);
	pr->kym2=(pr->iw+PAD_SIZE*2)*((ky-2)&3);
}
static void slic5_predict(SLIC5Ctx *pr, int kc, int kx)
{
	PROF(OUTSIDE);
	int idx=(pr->iw*pr->ky+kx)<<2|kc;
	pr->kc=kc;
	pr->kx=kx;
	pr->full_pixels_processed+=!kc;
	if(!(idx&(CDF_UPDATE_PERIOD-1))||(idx<CDF_UPDATE_PERIOD&&(idx&15)))
		slic5_update_CDFs(pr);
	PROF(UPDATE_CDFs);
	//XY are flipped, no need to check if indices OOB due to padding
#ifdef COMPENSATE_GAMMA
	const int *eq=pr->equalizer+((size_t)pr->maxlevels+1)*kc;
	int
		NNWW	=LOAD(pr->pixels,  2, 2)+pr->half[kc],
		NNW	=LOAD(pr->pixels,  1, 2)+pr->half[kc],
		NN	=LOAD(pr->pixels,  0, 2)+pr->half[kc],
		NNE	=LOAD(pr->pixels, -1, 2)+pr->half[kc],
		NNEE	=LOAD(pr->pixels, -2, 2)+pr->half[kc],
		NWW	=LOAD(pr->pixels,  2, 1)+pr->half[kc],
		NW	=LOAD(pr->pixels,  1, 1)+pr->half[kc],
		N	=LOAD(pr->pixels,  0, 1)+pr->half[kc],
		NE	=LOAD(pr->pixels, -1, 1)+pr->half[kc],
		NEE	=LOAD(pr->pixels, -2, 1)+pr->half[kc],
		WW	=LOAD(pr->pixels,  2, 0)+pr->half[kc],
		W	=LOAD(pr->pixels,  1, 0)+pr->half[kc];
	NNWW	=eq[CLAMP(0, NNWW	, pr->nlevels[kc]-1)];
	NNW	=eq[CLAMP(0, NNW	, pr->nlevels[kc]-1)];
	NN	=eq[CLAMP(0, NN		, pr->nlevels[kc]-1)];
	NNE	=eq[CLAMP(0, NNE	, pr->nlevels[kc]-1)];
	NNEE	=eq[CLAMP(0, NNEE	, pr->nlevels[kc]-1)];
	NWW	=eq[CLAMP(0, NWW	, pr->nlevels[kc]-1)];
	NW	=eq[CLAMP(0, NW		, pr->nlevels[kc]-1)];
	N	=eq[CLAMP(0, N		, pr->nlevels[kc]-1)];
	NE	=eq[CLAMP(0, NE		, pr->nlevels[kc]-1)];
	NEE	=eq[CLAMP(0, NEE	, pr->nlevels[kc]-1)];
	WW	=eq[CLAMP(0, WW		, pr->nlevels[kc]-1)];
	W	=eq[CLAMP(0, W		, pr->nlevels[kc]-1)];
#else
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
#ifdef COMPENSAGE_GAMMA_v2
	int mean=(2*(N+W+NW+NE)+NNWW+NNW+NN+NNE+NNEE+NWW+NEE+WW+8)>>4;//29-bit max, sign bit included
	long long temp, var=0;
	temp=(long long)NNWW	-mean, var+=temp*temp;
	temp=(long long)NNW	-mean, var+=temp*temp;
	temp=(long long)NN	-mean, var+=temp*temp;
	temp=(long long)NNE	-mean, var+=temp*temp;
	temp=(long long)NNEE	-mean, var+=temp*temp;
	temp=(long long)NWW	-mean, var+=temp*temp;
	temp=(long long)NW	-mean, var+=temp*temp<<1;
	temp=(long long)N	-mean, var+=temp*temp<<1;
	temp=(long long)NE	-mean, var+=temp*temp<<1;
	temp=(long long)NEE	-mean, var+=temp*temp;
	temp=(long long)WW	-mean, var+=temp*temp;
	temp=(long long)W	-mean, var+=temp*temp<<1;
	int sdev=floor_sqrt((var<<PRED_PREC)/15);
#endif
#endif
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
			pr->sse_idx_p[k]=quantize_signed((int)pr->pred, sh+8-(SSE_PREDBITS-1), SSE_P_EXP, SSE_P_MSB, pr->sse_p);
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
#ifdef COMPENSATE_GAMMA
	int L=0, R=pr->nlevels[kc]-1, found=0;
	while(L<=R)//binary search	O(depth)
	{
		int level=eq[pr->pred_final=(L+R)>>1];
		if(pr->pred>level)
			L=pr->pred_final+1;
		else if(pr->pred<level)
			R=pr->pred_final-1;
		else
		{
			found=1;
			break;
		}
	}
	pr->pred_final-=pr->half[kc];
	//pr->pred_final+=(int)((pr->pred+((1<<PRED_PREC)>>1)-1)>>PRED_PREC);
	//pr->pred_final>>=1;
#else
	pr->pred_final=(int)((pr->pred+((1<<PRED_PREC)>>1)-1)>>PRED_PREC);
#endif
}
static void slic5_update_hist(SLIC5Ctx *pr, int token)
{
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
}
static void slic5_update(SLIC5Ctx *pr, int curr, int token)
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
		if(pr->ky&&pr->kx+1<pr->iw)
			LOAD_PRED_ERROR(kc, -1, 1, k)+=errors[k];//eNE += ecurr
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

	//update token hist
	if(token>=0)
		slic5_update_hist(pr, token);

#ifdef COMPENSATE_GAMMA
	//update equalizer hist
	int *hist=pr->eq_hist+((size_t)kc<<pr->maxdepth);
	int ucurr=curr+pr->half[kc];
	ucurr=CLAMP(0, ucurr, pr->nlevels[kc]-1);
	++hist[ucurr];
#endif
	
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
	for(int k=0;k<NPREDS;++k)
		memset(pr->pred_errors+((NPREDS*idx+k)<<2), 0, lzlen*sizeof(int[4]));
#endif
}
#endif
#undef  LOAD
static void slic5_enc(SLIC5Ctx *pr, int curr, int kc, int kx)
{
	slic5_predict(pr, kc, kx);

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
	if(hu.token>=pr->cdfsize
#ifdef ENABLE_LZ
		-1//LZ escape symbol is OOB here
#endif
	)
		LOG_ERROR("Token OOB %d/%d", hu.token, pr->cdfsize-1);

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
static int slic5_dec(SLIC5Ctx *pr, int kc, int kx)
{
	slic5_predict(pr, kc, kx);

	int token=EC_DEC(pr->ec, pr->CDFs+(pr->cdfsize+1)*(pr->nhist*pr->nhist*kc+pr->hist_idx), pr->cdfsize, 0);
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
const char *rct_names[]=
{
#define RCT(X) #X,
	RCTLIST
#undef  RCT
};
#define r p[0]
#define g p[1]
#define b p[2]
static void rct_fwd(int *p, RCTType rct)
{
	switch(rct)//rounding is to prevent luma from getting OOB
	{
	case RCT_SubGreen:
		r-=g;
		b-=g;
		break;
	case RCT_JPEG2000:
		r-=g;		//r-g				[1     -1     0  ].RGB
		b-=g;		//b-g				[0     -1     1  ].RGB
		g+=(r+b+2)>>2;	//g+(r-g+b-g)/4 = r/4+g/2+b/4	[1/4    1/2   1/4].RGB
		break;
	case RCT_YCoCg_R:
		r-=b;		//co = r-b			diff(r, b)
		b+=r>>1;	//(r+b)/2
		g-=b;		//cg = g-(r+b)/2		diff(g, av(r, b))
		b+=g>>1;	//Y  = (r+b)/2 + (g-(r+b)/2)/2 = r/4+g/2+b/4	av(g, av(r, b))
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
	case RCT_YCbCr_R_v2:
		r-=g;		//Cr =	[1	-1	0].RGB
		g+=r>>1;	//	[1/2	1/2	0]
		b-=g;		//Cb =	[-1/2	-1/2	1]
		g+=(2*b-r+4)>>3;//Y  =	[1/4	1/2	1/4]	v2
		break;
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
		g+=(6*b+8)>>4;	//Y  =	[5/16	5/16	6/16]	v5
		break;
	case RCT_YCbCr_R_v6:
		r-=g;		//Cr =	[1	-1	0].RGB
		g+=r>>1;	//	[1/2	1/2	0]
		b-=g;		//Cb =	[-1/2	-1/2	1]
		g+=(14*b+16)>>5;//Y  =	[9/32	9/32	14/32]	v6
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
	case RCT_YCbCr_R_v5:
		g-=(6*b+8)>>4;
		b+=g;
		g-=r>>1;
		r+=g;
		break;
	case RCT_YCbCr_R_v6:
		g-=(14*b+16)>>5;
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
	memset(pixels, 0, src->iw*sizeof(int[4*RCT_COUNT]));
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
			int v[]={src->data[idx<<2|0], src->data[idx<<2|1], src->data[idx<<2|2]}, v2[3], N, W, NW, pred;

#define GET_PIXEL(RCT_IDX, C, X, Y) pixels[(src->iw*(RCT_IDX<<1|row##Y)+kx+1-X)<<2|C]
#define PROCESS_SUBPIXEL(RCT_IDX, C, X, N_IDX)\
	GET_PIXEL(RCT_IDX, C, 0, 0)=X,\
	N=GET_PIXEL(RCT_IDX, C, 0, 1),\
	W=GET_PIXEL(RCT_IDX, C, 1, 0),\
	NW=GET_PIXEL(RCT_IDX, C, 1, 1),\
	pred=N+W-NW, pred=MEDIAN3(N, W, pred), X-=pred,\
	++hist[maxlevels*(3*RCT_IDX+C)+((X+(nlevels[N_IDX][C]>>1))&(nlevels[N_IDX][C]-1))]
			
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

			//RCT_YCbCr_R_v2
			memcpy(v2, v, sizeof(v2));
			rct_fwd(v2, RCT_YCbCr_R_v2);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v2, 0, v2[1], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v2, 1, v2[2], 1);
			PROCESS_SUBPIXEL(RCT_YCbCr_R_v2, 2, v2[0], 1);

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
int t47_encode(Image const *src, ArrayHandle *data, int loud)
{
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
	RCTType rct=RCT_NONE;
	int *CDF=0;
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
	SLIC5Ctx pr;
	dlist_init(&list, 1, 1024, 0);
	ac_enc_init(&ec, &list);
	int success=slic5_init(&pr, src->iw, src->ih, nch, src->depth, &ec);
#ifdef ENABLE_LZ
	int *lzmap=(int*)malloc(sizeof(int[LZ_MAP_SIZE]));
#endif
	if(!success
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
		dlist_push_back1(&list, &rct);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		slic5_nextrow(&pr, ky);
		for(int kx=0;kx<src->iw;)
		{
			if(kx==(src->iw>>1)&&ky==(src->ih>>1))//
				printf("");

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
					rct_fwd(comp, rct);
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
	if(loud)
	{
		double usize=image_getBMPsize(src);
		printf("\n");
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");

		printf("csize %8d  CR %10.6lf\n", (int)list.nobj, usize/list.nobj);

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
	SLIC5Ctx pr;
	ac_dec_init(&ec, data, data+srclen);
	int success=slic5_init(&pr, dst->iw, dst->ih, nch, dst->depth, &ec);
	if(!success)
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
			//	printf("");

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
					int comp[3];
					comp[1]=first;//Y
					comp[2]=slic5_dec(&pr, 1, kx);//Cb
					comp[0]=slic5_dec(&pr, 2, kx);//Cr
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
					int pixel=kc?slic5_dec(&pr, kc, kx):first;
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
	long long *hist=(long long*)malloc(sizeof(long long[NPREDS<<(LG_HIST_SIZE<<1)]));
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(hist, 0, sizeof(long long[NPREDS<<(LG_HIST_SIZE<<1)]));
	int histcount=sizeof(long long[NPREDS<<(LG_HIST_SIZE<<1)])/sizeof(long long);
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
				for(int kp=0;kp<NPREDS;++kp)
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
	for(int kp=1;kp<NPREDS;++kp)
	{
		long long *curr_hist=hist+((size_t)kp<<(LG_HIST_SIZE<<1));
		int result=memcmp(hist, curr_hist, sizeof(long long[1<<(LG_HIST_SIZE<<1)]));
		printf("memcmp %2d: %d\n", kp, result);
	}
#endif
	printf("Histograms:\n");
	for(int kp=0;kp<NPREDS;++kp)
	{
		long long *curr_hist=hist+((size_t)kp<<(LG_HIST_SIZE<<1));
		printf("%s\n", prednames[kp]);
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