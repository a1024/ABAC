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

//debug
//	#define ENABLE_GUIDE
//choose one:
	#define RCT_OVERRIDE rct_curr
//	#define RCT_OVERRIDE RCT_YCbCr_R_v1

//experiments
//	#define PRINT_HIST
//	#define PROFILER//SLOW
//	#define TRACK_SSE_RANGES//SLOW

//efficiency
	#define ENABLE_SSE_4D//3% smaller, but 2x slower
//	#define ENABLE_CfL//must be disabled here
//	#define SSE_UPDATE_NB_CELLS//inferior

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

#define SLIC7_CONFIG_EXP 5
#define SLIC7_CONFIG_MSB 2
#define SLIC7_CONFIG_LSB 0

#define CDF_UPDATE_PERIOD 0x400//number of sub-pixels processed between CDF updates, must be a power-of-two
#define PAD_SIZE 3
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
#ifdef ENABLE_SSE_4D
#define SSE_SIZE (SSE_W*SSE_H*SSE_D<<SSE_PREDBITS)
#endif
//#define PRED_PREC 8
#define PARAM_PREC 8

#define NPREDS 14
#define PREDLIST\
	PRED(W+NE-N-((2*(eN+eW)+eNE-eNW+4)>>3))\
	PRED(N-(int)(((long long)eN+eW+eNE)*-0x05C>>PARAM_PREC))\
	PRED(W-(int)(((long long)eN+eW+eNW)*-0x05B>>PARAM_PREC))\
	PRED(N+(int)((-eNN*0x0DFLL-eN*0x0051LL-eNE*0x0BDLL+((long long)N-NN)*0x05C+((long long)NW-W)*0x0102)>>PARAM_PREC))\
	PRED((3*W+NEEE+2)>>2)\
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
//	PRED(experimental)
//	PRED(floor_sqrt(floor_sqrt((N+0x800000LL)*(W+0x800000LL))*floor_sqrt((NW+0x800000LL)*(NE+0x800000LL)))-0x800000)

static const char *prednames[]=
{
#define PRED(X) #X,
	PREDLIST
#undef  PRED
};
#define EC_ENC(EC, X, CDF, NLEVELS, LG_FMIN) ac_enc(EC, X, CDF, NLEVELS, LG_FMIN)
#define EC_DEC(EC, CDF, NLEVELS, LG_FMIN) ac_dec(EC, CDF, NLEVELS, LG_FMIN)
typedef unsigned CDF_t;
#define CDF_SHIFT 16
typedef struct HybridUintStruct
{
	unsigned short token, nbits;
	unsigned bypass;
} HybridUint;

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
static void hybriduint_encode(unsigned val, HybridUint *hu)
{
	int token, bypass, nbits;
	if(val<(1<<SLIC7_CONFIG_EXP))
	{
		token=val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=floor_log2_32((unsigned)val);
		int mantissa=val-(1<<lgv);
		token = (1<<SLIC7_CONFIG_EXP) + (
				(lgv-SLIC7_CONFIG_EXP)<<(SLIC7_CONFIG_MSB+SLIC7_CONFIG_LSB)|
				(mantissa>>(lgv-SLIC7_CONFIG_MSB))<<SLIC7_CONFIG_LSB|
				(mantissa&((1<<SLIC7_CONFIG_LSB)-1))
			);
		nbits=lgv-(SLIC7_CONFIG_MSB+SLIC7_CONFIG_LSB);
		bypass=val>>SLIC7_CONFIG_LSB&((1LL<<nbits)-1);
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

typedef struct SLIC7CtxStruct
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
	int nhist, cdfsize;
	int *pred_errors, *errors, *pixels;

	int *hist, *histsums;
	CDF_t *CDFs;

#ifdef ENABLE_SSE_4D
	int sse_w, sse_h, sse_d, sse_p, sse_size;
	long long *sse, *sse_fr;
#ifdef ENABLE_CfL
	long long *sse_cfl;
#endif
#ifdef SSE_UPDATE_NB_CELLS
	int sse_idx_x[(SSE_STAGES+6)<<2], sse_idx_y[(SSE_STAGES+6)<<2], sse_idx_z[(SSE_STAGES+6)<<2];
#endif
	int sse_idx_p[(SSE_STAGES+6)<<2];
	int sse_idx[(SSE_STAGES+6)<<2];
	long long sse_sum[(SSE_STAGES+6)<<2];
	int sse_count[(SSE_STAGES+6)<<2];
#endif
	int bias_count[4];
	long long bias_sum[4];

	int preds[NPREDS<<2], params[NPREDS<<2];
	long long pred[4];
	int pred_final[4];
	int hist_idx[4];
	int sse_corr[4];
	int ky, kym0, kym1, kym2;
	long long pred_error_sums[NPREDS];
#ifdef TRACK_SSE_RANGES
	int sse_ranges[8];
#endif
	ArithmeticCoder *ec;
} SLIC7Ctx;
#define LOAD(BUF, X, Y) BUF[(pr->kym##Y+kx+PAD_SIZE-(X))<<2|kc]
#define LOAD_CH(BUF, C, X, Y) BUF[(pr->kym##Y+kx+PAD_SIZE-(X))<<2|C]
#define LOAD_PRED_ERROR(C, X, Y, P) pr->pred_errors[(NPREDS*(pr->kym##Y+kx+PAD_SIZE-(X))+P)<<2|C]
static int slic7_init(SLIC7Ctx *pr, int iw, int ih, int nch, const char *depths, ArithmeticCoder *ec)
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
		pr->depths[0]=depths[1];	//Y
		pr->depths[1]=depths[2]+1;	//Cb
		pr->depths[2]=depths[0]+1;	//Cr
		if(nch==4)
			pr->depths[3]=depths[3];//a
	}
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
	HybridUint hu;
	int extremesym=-(pr->maxlevels>>1);
	extremesym=extremesym<<1^-(extremesym<0);//pack sign
	hybriduint_encode(extremesym, &hu);//encode -half
	pr->cdfsize=hu.token+1;
	pr->nhist=QUANTIZE_HIST(767);
	
#ifdef ENABLE_SSE_4D
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
#ifdef ENABLE_SSE_4D
	pr->sse=(long long*)malloc(nch*sizeof(long long[SSE_SIZE*SSE_STAGES]));
	pr->sse_fr=(long long*)malloc(nch*sizeof(long long[SSE_FR_SIZE]));
#ifdef ENABLE_CfL
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
#endif
	if(!pr->pred_errors||!pr->errors||!pr->pixels||!pr->hist||!pr->histsums||!pr->CDFs
#ifdef ENABLE_SSE_4D
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
			int hist_idx[kc]=pr->nhist*qy+qx;
			CDF_t *curr_CDF=pr->CDFs+(pr->cdfsize+1)*hist_idx[kc];
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
	
#ifdef ENABLE_SSE_4D
	memset(pr->sse, 0, nch*sizeof(long long[SSE_SIZE*SSE_STAGES]));
	memset(pr->sse_fr, 0, nch*sizeof(long long[SSE_FR_SIZE]));
#ifdef ENABLE_CfL
	if(pr->nch>1)
		memset(pr->sse_cfl, 0, (nch-1LL)*sizeof(long long[SSE_SIZE<<1]));
#endif
#endif

	for(int k=0;k<(NPREDS<<2);++k)
		pr->params[k]=12<<pr->maxdepth;
	//memfill(pr->params+NPREDS, pr->params, sizeof(pr->params)-sizeof(int[NPREDS]), sizeof(int[NPREDS]));
	//int shift=16-maxdepth;
	//shift=MAXVAR(0, shift);
	//memcpy(pr->params, wp_params, sizeof(pr->params));
	//for(int k=0;k<NPREDS;++k)
	//	pr->params[k]=pr->params[k]>>shift;

	pr->ec=ec;
	return 1;
}
static void slic7_free(SLIC7Ctx *pr)
{
	free(pr->pred_errors);
	free(pr->errors);
	free(pr->pixels);
	free(pr->hist);
	free(pr->histsums);
	free(pr->CDFs);
#ifdef ENABLE_SSE_4D
	free(pr->sse);
	free(pr->sse_fr);
#ifdef ENABLE_CfL
	if(pr->nch>1)
		free(pr->sse_cfl);
#endif
#endif
	//memset(pr, 0, sizeof(*pr));//init does this already
}
static void slic7_update_CDFs(SLIC7Ctx *pr)
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
static void slic7_nextrow(SLIC7Ctx *pr, int ky)
{
	pr->ky=ky;
	pr->kym0=(pr->iw+PAD_SIZE*2)*(ky&3);
	pr->kym1=(pr->iw+PAD_SIZE*2)*((ky-1)&3);
	pr->kym2=(pr->iw+PAD_SIZE*2)*((ky-2)&3);
}
static void slic7_predict(SLIC7Ctx *pr, int kc, int kx)
{
	PROF(OUTSIDE);
	int idx=(pr->iw*pr->ky+kx)<<2|kc;
	pr->full_pixels_processed+=!kc;
	if(!(idx&(CDF_UPDATE_PERIOD-1))||(idx<CDF_UPDATE_PERIOD&&(idx&15)))
		slic7_update_CDFs(pr);
	PROF(UPDATE_CDFs);
	//XY are flipped, no need to check if indices OOB due to padding
	int
		NNWW	=LOAD(pr->pixels,  2, 2),
		NNW	=LOAD(pr->pixels,  1, 2),
		NN	=LOAD(pr->pixels,  0, 2),
		NNE	=LOAD(pr->pixels, -1, 2),
		NNEE	=LOAD(pr->pixels, -2, 2),
		NWW	=LOAD(pr->pixels,  2, 1),
		NW	=LOAD(pr->pixels,  1, 1),
		N	=LOAD(pr->pixels,  0, 1),
		NE	=LOAD(pr->pixels, -1, 1),
		NEE	=LOAD(pr->pixels, -2, 1),
#if PAD_SIZE>=3
		NEEE	=LOAD(pr->pixels, -3, 1),
#endif
		WW	=LOAD(pr->pixels,  2, 0),
		W	=LOAD(pr->pixels,  1, 0);
	int
		eNNWW	=LOAD(pr->errors,  2, 2),//error = curr - pred
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

	//int experimental;
	//long long aN=N+0x800000LL, aW=W+0x800000LL;
	//aN=CLAMP(0, aN, 0x1000000);
	//aW=CLAMP(0, aW, 0x1000000);
	//experimental=floor_sqrt(aN*aW)-0x800000;

	int j=-1;
#define PRED(X) pr->preds[++j<<2|kc]=X;
	PREDLIST
#undef  PRED
	PROF(CALC_SUBPREDS);
#if 1
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

		weights[k]=((long long)pr->params[k<<2|kc]<<29)/(weights[k]+1);
		wsum+=weights[k];
	}
	if(wsum)
	{
		pr->pred[kc]=0;
		for(int k=0;k<NPREDS;++k)
			pr->pred[kc]+=weights[k]*pr->preds[k<<2|kc];
		pr->pred[kc]+=(wsum>>1)-1;
		pr->pred[kc]/=wsum;
	}
	else
#endif
		pr->pred[kc]=pr->preds[10<<2|kc];//clamped grad
	PROF(CALC_WEIGHT_AV);

	int qx=QUANTIZE_HIST((dx+abs(eW+eWW))>>(sh-2)), qy=QUANTIZE_HIST((dy+abs(eN+eNN))>>(sh-2));
	qx=CLAMP(0, qx, pr->nhist-1);
	qy=CLAMP(0, qy, pr->nhist-1);
	pr->hist_idx[kc]=pr->nhist*qy+qx;

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
#ifdef ENABLE_CfL
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
#endif
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
#ifdef SSE_UPDATE_NB_CELLS
#define QUANTIZE(IDX, X, Y, Z)\
	SSE_W*(SSE_H*(pr->sse_idx_z[IDX<<2|kc]=quantize_signed(Z, sh, SSE_Z_EXP, SSE_Z_MSB, SSE_D))+\
	(pr->sse_idx_y[IDX<<2|kc]=quantize_signed(Y, sh, SSE_Y_EXP, SSE_Y_MSB, SSE_H)))+\
	(pr->sse_idx_x[IDX<<2|kc]=quantize_signed(X, sh, SSE_X_EXP, SSE_X_MSB, SSE_W))
#else
#define QUANTIZE(IDX, X, Y, Z)\
	SSE_W*(SSE_H*quantize_signed(Z, sh, SSE_Z_EXP, SSE_Z_MSB, SSE_D)+\
	quantize_signed(Y, sh, SSE_Y_EXP, SSE_Y_MSB, SSE_H))+\
	quantize_signed(X, sh, SSE_X_EXP, SSE_X_MSB, SSE_W)
#endif
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

	pr->sse_corr[kc]=0;
#ifdef ENABLE_CfL
	for(int k=0;k<SSE_STAGES+1+(kc<<1);++k)
#else
	for(int k=0;k<SSE_STAGES+1;++k)
#endif
	{
		long long sse_val;
		if(k<SSE_STAGES)
		{
			pr->sse_idx_p[k<<2|kc]=quantize_signed((int)pr->pred[kc], sh+8-(SSE_PREDBITS-1), SSE_P_EXP, SSE_P_MSB, pr->sse_p);//shift is sh+8-(SSE_PREDBITS-1): remove the 8-bits, and insert up to SSE_PREDBITS
			pr->sse_idx[k<<2|kc]=SSE_SIZE*(SSE_STAGES*kc+k)+(SSE_W*SSE_H*SSE_D)*pr->sse_idx_p[k<<2|kc]+g2[k];
			sse_val=pr->sse[pr->sse_idx[k<<2|kc]];
			
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
#ifdef ENABLE_CfL
		else if(k<SSE_STAGES+(kc<<1))
		{
			pr->sse_idx_p[k<<2|kc]=quantize_signed((int)pr->pred[kc], sh+8-(SSE_PREDBITS-1), SSE_P_EXP, SSE_P_MSB, pr->sse_p);
			pr->sse_idx[k<<2|kc]=SSE_SIZE*(k-SSE_STAGES)+(SSE_W*SSE_H*SSE_D)*pr->sse_idx_p[k<<2|kc]+g2[k];
			sse_val=pr->sse_cfl[pr->sse_idx[k<<2|kc]];
		}
#endif
		else
		{
			int temp;
			int idx=0;
		//	idx=3*idx+THREEWAY(eN, eW)+1;//X
		//	idx=3*idx+THREEWAY(eN, 0)+1;//X
		//	idx=3*idx+THREEWAY(eW, 0)+1;//X
		//	idx=3*idx+THREEWAY(N, W)+1;//X
			idx=3*idx+THREEWAY(N, (int)pr->pred[kc])+1;
			idx=3*idx+THREEWAY(W, (int)pr->pred[kc])+1;
			idx=3*idx+THREEWAY(NW, (int)pr->pred[kc])+1;
			idx=3*idx+THREEWAY(NE, (int)pr->pred[kc])+1;
			idx=3*idx+THREEWAY(NN, (int)pr->pred[kc])+1;
			idx=3*idx+THREEWAY(WW, (int)pr->pred[kc])+1;
			temp=N+W-NW, idx=3*idx+THREEWAY(temp, (int)pr->pred[kc])+1;
			temp=(4*(W+N)-(WW+NN)+6)/12, idx=3*idx+THREEWAY(temp, (int)pr->pred[kc])+1;
			temp=2*N-NN, idx=3*idx+THREEWAY(temp, (int)pr->pred[kc])+1;
			temp=2*W-WW, idx=3*idx+THREEWAY(temp, (int)pr->pred[kc])+1;
			if(idx>=SSE_FR_SIZE)
				LOG_ERROR("SSE OOB");
			pr->sse_idx[k<<2|kc]=SSE_FR_SIZE*kc+idx;
			//pr->sse_idx[k]=
			//	//(W<NW)<<11|
			//	//(NW<NE)<<10|
			//	(N<W)<<9|
			//	(N+W-NW<(int)pr->pred[kc])<<8|
			//	(W <pr->pred[kc])<<7|
			//	(NW<pr->pred[kc])<<6|
			//	(N <pr->pred[kc])<<5|
			//	(NE<pr->pred[kc])<<4|
			//	(WW<pr->pred[kc])<<3|
			//	(NN<pr->pred[kc])<<2|
			//	(2*N-NN<(int)pr->pred[kc])<<1|
			//	(2*W-WW<(int)pr->pred[kc])
			//;
			sse_val=pr->sse_fr[pr->sse_idx[k<<2|kc]];
		}
		pr->sse_count[k<<2|kc]=(int)(sse_val&0xFFF);
		pr->sse_sum[k<<2|kc]=(int)(sse_val>>12);
		int sse_corr=(int)(pr->sse_sum[k<<2|kc]*SSE_STEP/(pr->sse_count[k<<2|kc]+1LL));
		pr->pred[kc]+=sse_corr;
		pr->sse_corr[kc]+=sse_corr;
	}
#endif

	int final_corr=(int)(pr->bias_sum[kc]/(pr->bias_count[kc]+1LL));
	pr->pred[kc]+=final_corr;
	pr->sse_corr[kc]+=final_corr;
	PROF(SSE_LOOP);

	pr->pred[kc]=CLAMP(clamp_lo, pr->pred[kc], clamp_hi);
	pr->pred_final[kc]=(int)((pr->pred[kc]+pr->rounding_offset[kc])>>pr->shift_prec[kc]);
}
#if 0
static void slic7_update_hist(SLIC7Ctx *pr, int token, int kc)
{
	int total_hists=pr->nhist*pr->nhist;
	++pr->hist[pr->cdfsize*(total_hists*kc+pr->hist_idx[kc])+token];
	++pr->histsums[total_hists*kc+pr->hist_idx[kc]];
	if(pr->histsums[pr->hist_idx[kc]]>(pr->iw*pr->nch<<1))//rescale hist
	{
		int *curr_hist=pr->hist+pr->cdfsize*(total_hists*kc+pr->hist_idx[kc]);
		int sum=0;
		for(int ks=0;ks<pr->cdfsize;++ks)
		{
			curr_hist[ks]=(curr_hist[ks]+1)>>1;
			sum+=curr_hist[ks];
		}
		pr->histsums[pr->hist_idx[kc]]=sum;
	}
	PROF(UPDATE_HIST);
}
#endif
static int slic7_gettoken(SLIC7Ctx *pr, int curr, int kc, HybridUint *hu)
{
	int error=curr-pr->pred_final[kc], e0=error;
#if 1
	int upred=pr->half[kc]-abs(pr->pred_final[kc]);
	int abs_error=abs(error), negmask;
	if(abs_error<=upred)
	{
		negmask=-(pr->sse_corr[kc]<0);//sign is flipped if SSE correction was negative, to skew the histogram
		error^=negmask;
		error-=negmask;
		error=error<<1^-(error<0);//pack sign
	}
	else
		error=upred+abs_error;//error sign is known
#endif
#if 0
	int upred=pr->pred_final[kc]+(pr->nlevels[kc]>>1), abs_error=abs(error), negmask;
	if(upred<(pr->nlevels[kc]>>1))
	{
		if(abs_error<=upred)
		{
			negmask=-(pr->sse_corr[kc]<0);//sign is flipped if SSE correction was negative, to skew the histogram
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
			negmask=-(pr->sse_corr[kc]<0);
			error^=negmask;
			error-=negmask;
			error=error<<1^-(error<0);
		}
		else
			error=pr->nlevels[kc]-upred+abs(error);
	}
#endif

	hybriduint_encode(error, hu);
	if(hu->token>=pr->cdfsize)
		LOG_ERROR("Token OOB %d/%d", hu->token, pr->cdfsize);
	return e0;
}
static int slic7_update(SLIC7Ctx *pr, int curr, int kc, int kx)
{
	PROF(PRED_TILL_UPDATE);

	//update token hist
	HybridUint hu;
	slic7_gettoken(pr, curr, kc, &hu);
	//slic7_update_hist(pr, hu.token, kc);
	int total_hists=pr->nhist*pr->nhist;
	++pr->hist[pr->cdfsize*(total_hists*kc+pr->hist_idx[kc])+hu.token];
	++pr->histsums[total_hists*kc+pr->hist_idx[kc]];
	if(pr->histsums[pr->hist_idx[kc]]>(pr->iw*pr->nch<<1))//rescale hist
	{
		int *curr_hist=pr->hist+pr->cdfsize*(total_hists*kc+pr->hist_idx[kc]);
		int sum=0;
		for(int ks=0;ks<pr->cdfsize;++ks)
		{
			curr_hist[ks]=(curr_hist[ks]+1)>>1;
			sum+=curr_hist[ks];
		}
		pr->histsums[pr->hist_idx[kc]]=sum;
	}
	PROF(UPDATE_HIST);

	//update current pixel
	curr<<=pr->shift_prec[kc];
	LOAD(pr->pixels, 0, 0)=curr;

	//update WP errors
	int error=curr-(int)pr->pred[kc];
	LOAD(pr->errors, 0, 0)=error;
	int errors[NPREDS]={0}, kbest=0;
	for(int k=0;k<NPREDS;++k)
	{
		errors[k]=abs(curr-pr->preds[k<<2|kc]);
		LOAD_PRED_ERROR(kc, 0, 0, k)=errors[k];
		if(pr->ky&&kx+1<pr->iw)
			LOAD_PRED_ERROR(kc, -1, 1, k)+=errors[k];//eNE += ecurr
		if(errors[kbest]>errors[k])
			kbest=k;

		pr->pred_error_sums[k]+=errors[k];
	}

	//update WP weights
	++pr->params[kbest<<2|kc];
	if(pr->params[kbest<<2|kc]>(12<<pr->depths[kc]))
	{
		for(int k=0;k<NPREDS;++k)
			pr->params[k<<2|kc]>>=1;
	}
	PROF(UPDATE_WP);
	
	//update SSE
#ifdef ENABLE_SSE_4D
	for(int k=0;k<SSE_STAGES+1+(kc<<1);++k)
	{
		long long *sse_ptr;
		if(k<SSE_STAGES)
			sse_ptr=pr->sse;
#ifdef ENABLE_CfL
		else if(k<SSE_STAGES+(kc<<1))
			sse_ptr=pr->sse_cfl;
#endif
		else
			sse_ptr=pr->sse_fr;
		pr->sse_count[k<<2|kc]+=SSE_STEP;
		pr->sse_sum[k<<2|kc]+=error;
		if(pr->sse_count[k<<2|kc]>=SSE_LIMIT)
		{
			pr->sse_count[k<<2|kc]>>=1;
			pr->sse_sum[k<<2|kc]>>=1;
		}
		long long sse_val=pr->sse_sum[k<<2|kc]<<12|pr->sse_count[k<<2|kc];
		sse_ptr[pr->sse_idx[k<<2|kc]]=sse_val;
		
#ifdef SSE_UPDATE_NB_CELLS
		int e2=error/SSE_STEP;
		if(e2&&k<SSE_STAGES+(kc<<1))
		{
			int idx2, count;
			long long sum;
			if(pr->sse_idx_y[k<<2|kc]>0)
			{
				idx2=SSE_W*(SSE_H*(SSE_D*pr->sse_idx_p[k<<2|kc]+pr->sse_idx_z[k<<2|kc])+pr->sse_idx_y[k<<2|kc]-1)+pr->sse_idx_x[k<<2|kc];
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
			if(pr->sse_idx_y[k<<2|kc]<SSE_H-1)
			{
				idx2=SSE_W*(SSE_H*(SSE_D*pr->sse_idx_p[k<<2|kc]+pr->sse_idx_z[k<<2|kc])+pr->sse_idx_y[k<<2|kc]+1)+pr->sse_idx_x[k<<2|kc];
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
			if(pr->sse_idx_z[k<<2|kc]>0)
			{
				idx2=SSE_W*(SSE_H*(SSE_D*pr->sse_idx_p[k<<2|kc]+pr->sse_idx_z[k<<2|kc]-1)+pr->sse_idx_y[k<<2|kc])+pr->sse_idx_x[k<<2|kc];
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
			if(pr->sse_idx_z[k<<2|kc]<SSE_D-1)
			{
				idx2=SSE_W*(SSE_H*(SSE_D*pr->sse_idx_p[k<<2|kc]+pr->sse_idx_z[k<<2|kc]+1)+pr->sse_idx_y[k<<2|kc])+pr->sse_idx_x[k<<2|kc];
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
	return error;
}
#undef  LOAD
static int slic7_enc(SLIC7Ctx *prs, int nencoders, RCTType rct, int curr, int kc, int kx)
{
	SLIC7Ctx *pr=prs+rct;
	for(int k=0;k<nencoders;++k)
		slic7_predict(prs+k, kc, kx);

	HybridUint hu;
	int e0=slic7_gettoken(pr, curr, kc, &hu);

	EC_ENC(pr->ec, hu.token, pr->CDFs+(pr->cdfsize+1)*(pr->nhist*pr->nhist*kc+pr->hist_idx[kc]), pr->cdfsize, 0);
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

	if(nencoders==1)
		slic7_update(prs, curr, kc, kx);
	return e0;
}
static int slic7_dec(SLIC7Ctx *prs, int nencoders, RCTType rct, int kc, int kx)
{
	SLIC7Ctx *pr=prs+rct;
	for(int k=0;k<nencoders;++k)
		slic7_predict(prs+k, kc, kx);

	int token=EC_DEC(pr->ec, pr->CDFs+(pr->cdfsize+1)*(pr->nhist*pr->nhist*kc+pr->hist_idx[kc]), pr->cdfsize, 0);
	int error=token;
	if(error>=(1<<SLIC7_CONFIG_EXP))
	{
		error-=1<<SLIC7_CONFIG_EXP;
		int lsb=error&((1<<SLIC7_CONFIG_LSB)-1);
		error>>=SLIC7_CONFIG_LSB;
		int msb=error&((1<<SLIC7_CONFIG_MSB)-1);
		error>>=SLIC7_CONFIG_MSB;
		int nbits=error+SLIC7_CONFIG_EXP-(SLIC7_CONFIG_MSB+SLIC7_CONFIG_LSB), n=nbits;
		int bypass=0;
		while(n>8)
		{
			n-=8;
			bypass|=EC_DEC(pr->ec, 0, 1<<8, 16-8)<<n;
		}
		bypass|=EC_DEC(pr->ec, 0, 1<<n, 16-n);
		error=1;
		error<<=SLIC7_CONFIG_MSB;
		error|=msb;
		error<<=nbits;
		error|=bypass;
		error<<=SLIC7_CONFIG_LSB;
		error|=lsb;
	}
#if 1
	int upred=pr->half[kc]-abs(pr->pred_final[kc]), negmask;
	if(error<=(upred<<1))
	{
		negmask=-(pr->sse_corr[kc]<0);
		error=error>>1^-(error&1);
	}
	else
	{
		negmask=-(pr->pred_final[kc]>0);
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
			negmask=-(pr->sse_corr[kc]<0);
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
			negmask=-(pr->sse_corr[kc]<0);
			error=error>>1^-(error&1);
			error^=negmask;
			error-=negmask;
		}
		else
			error=upred-error;//here error is negative if prediction is positive
	}
#endif
	int curr=error+pr->pred_final[kc];
	if(nencoders==1)
		slic7_update(prs, curr, token, kx);
	return curr;
}
static RCTType slic7_choose_rct(SLIC7Ctx *prs, int nencoders, const int *comp, int *rct_errors, int kx, int ky, int iw, int ih)
{
	static const int weights[]=
	{
		2,  3,  4,  6,  8,  8,  6,  4,  3,  2,//v3
		3,  4,  6,  8, 12, 12,  8,  6,  4,  3,
		4,  6,  8, 12, 16, 16, 12,  8,  6,  4,
		6,  8, 12, 16, 16, 16, 16, 12,  8,  6,
		8, 12, 16, 16,//[?]

		//0, 0, 0, 1, 1, 3, 1, 1, 1, 0,//v2, 3/4 POT
		//1, 2, 2, 2, 3, 4, 3, 2, 2, 1,
		//2, 3, 3, 3, 4, 7, 4, 3, 3, 2,
		//3, 4, 4, 4, 7, 9, 7, 4, 4, 3,
		//4, 4, 4, 7,

		//1,  2,  3,  4,  6,  4,  3,  2,  1,  0,//v1
		//2,  3,  4,  6,  8,  6,  4,  3,  2,  1,
		//3,  4,  6,  8, 12,  8,  6,  4,  3,  2,
		//4,  6,  8, 12, 16, 12,  8,  6,  4,  3,
		//6,  8, 12, 16,

		//1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		//1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		//1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		//1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		//1, 1, 1, 1,
	};
	int v2[3];
	long long ebest=0;
	int best_rct=0;
	//int ky0=ky&3, ky1=(ky-1)&3, ky2=(ky-2)&3;
	for(int k=0;k<nencoders;++k)
	{
		memcpy(v2, comp, sizeof(v2));
		rct_fwd(v2, k);
		int error=0;
		for(int kc2=0;kc2<3;++kc2)
		{
			int e=slic7_update(prs+k, v2[(kc2+1)%3], kc2, kx);
			error+=abs(e);
		}

		long long esum=error*weights[_countof(weights)-1];
		if(ky>=4)
		{
			for(int kx2=0;kx2<=5;++kx2)
			{
				if((unsigned)(kx+kx2)<(unsigned)iw)
					esum+=(long long)rct_errors[RCT_COUNT*(iw*(ky&3)+kx+kx2)+k]*weights[kx2+3];
			}
		}
		for(int ky2=-3;ky2<=0;++ky2)
		{
			for(int kx2=-3;kx2<=5;++kx2)
			{
				if((unsigned)(ky+ky2)<(unsigned)ih&&(unsigned)(kx+kx2)<(unsigned)iw)
					esum+=(long long)rct_errors[RCT_COUNT*(iw*((ky+ky2)&3)+kx+kx2)+k]*weights[10*(ky2+4)+kx2+3];
				if(!ky2&&!kx2)
					break;
			}
		}
		rct_errors[RCT_COUNT*(iw*(ky&3)+kx)+k]=error;
		if(!k||ebest>esum)
			ebest=esum, best_rct=k;

		//int eW=kx?rct_errors[RCT_COUNT*(kx-1)+k]:0;
		//int eNEE=kx<iw-2?rct_errors[RCT_COUNT*(kx+2)+k]:0;
		//int ecurr=(2*eW+2*eNEE+error)/5;
		//rct_errors[RCT_COUNT*kx+k]=ecurr;
		//if(!k||ebest>ecurr)
		//	ebest=ecurr, best_rct=k;

		//if(!k||ebest>rct_errors[RCT_COUNT*kx+k])
		//	ebest=rct_errors[RCT_COUNT*kx+k], best_rct=k;
	}
	return (RCTType)best_rct;
}

#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
int t49_encode(Image const *src, ArrayHandle *data, int loud)
{
#ifdef ENABLE_GUIDE
	guide=src;
#endif
	double t_start=time_sec();
	double usize=image_getBMPsize(src);
	int nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);//TODO
	if(nch>src->nch)
		nch=src->nch;
	if(loud)
	{
		int maxdepth=calc_maxdepth(src, 0);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T49 SLIC7  Enc %s  CWHD %d*%d*%d*%d/8\n", g_buf, nch, src->iw, src->ih, maxdepth);
	}
	DList list;
	ArithmeticCoder ec;
	int nencoders=nch>=3?RCT_COUNT:1;
	SLIC7Ctx *pr=(SLIC7Ctx*)malloc(nencoders*sizeof(SLIC7Ctx));
	int *rct_errors=nch>=3?(int*)malloc(src->iw*sizeof(int[RCT_COUNT<<2])):0;//4 rows * N_RCT
	if(!pr||nch>=3&&!rct_errors)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	dlist_init(&list, 1, 1024, 0);
	ac_enc_init(&ec, &list);
	for(int k=0;k<nencoders;++k)
	{
		int success=slic7_init(pr+k, src->iw, src->ih, nch, src->depth, &ec);
		if(!success)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
	}
	RCTType rct_curr=RCT_SubGreen, rct0=RCT_SubGreen;
	if(nch>=3)
		memset(rct_errors, 0, src->iw*sizeof(int[RCT_COUNT<<2]));
	//long long rct_errors[RCT_COUNT]={0};//crude solution, RCT selection should depend on local-average RCT errors
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		//memset(rct_errors, 0, sizeof(rct_errors));
		for(int k=0;k<nencoders;++k)
			slic7_nextrow(pr+k, ky);
		for(int kx=0;kx<src->iw;)
		{
			int kc=0;
			if(nch>=3)
			{
				int comp[]={src->data[idx|0], src->data[idx|1], src->data[idx|2]}, v2[3];
				if(!kx)
					rct_curr=rct0;
				memcpy(v2, comp, sizeof(v2));
				rct_fwd(v2, RCT_OVERRIDE);
				int eY=slic7_enc(pr, nencoders, RCT_OVERRIDE, v2[1], 0, kx);//Y		luma is encoded first
				int eU=slic7_enc(pr, nencoders, RCT_OVERRIDE, v2[2], 1, kx);//U
				int eV=slic7_enc(pr, nencoders, RCT_OVERRIDE, v2[0], 2, kx);//V
				
				if(ky==256)
					printf("%3d  %-11s  eYUV %d %d %d\n", kx, rct_names[rct_curr], eY, eU, eV);//

				rct_curr=slic7_choose_rct(pr, nencoders, comp, rct_errors, kx, ky, src->iw, src->ih);
				if(!kx)
					rct0=rct_curr;
				kc+=3;
			}
			for(;kc<nch;++kc)//gray/alpha
			{
				int pixel=src->data[idx|kc];
				slic7_enc(pr, 1, RCT_NONE, pixel, kc, kx);
			}
			++kx;
			idx+=4;
		}
		//if(loud&&!(ky&15))
		if(loud)
			printf("%5d  %6.2lf%%  %10lf sec  %14lld  %12lf  %-11s\n",
				ky, 100.*(ky+1)/src->ih, time_sec()-t_start, list.nobj, (ky+1)*usize/(src->ih*list.nobj), rct_names[RCT_OVERRIDE]
			);
	}
	ac_enc_flush(&ec);
	dlist_appendtoarray(&list, data);
	if(loud)
	{
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

		//printf("Pred errors\n");
		//long long sum=0;
		//for(int k=0;k<NPREDS;++k)
		//	sum+=pr.pred_error_sums[k];
		//for(int k=0;k<NPREDS;++k)
		//	printf("  %2d  %17lld  %6.2lf%%  %s\n", k, pr.pred_error_sums[k], 100.*pr.pred_error_sums[k]/sum, prednames[k]);
		
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
	dlist_clear(&list);
	for(int k=0;k<nencoders;++k)
		slic7_free(pr+k);
	free(pr);
	return 1;
}
int t49_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	double t_start=time_sec();
	int dst_nch=(dst->depth[0]!=0)+(dst->depth[1]!=0)+(dst->depth[2]!=0)+(dst->depth[3]!=0);
	int nch=MINVAR(dst_nch, dst->nch);
	ArithmeticCoder ec;
	ac_dec_init(&ec, data, data+srclen);
	int nencoders=nch>=3?RCT_COUNT:1;
	SLIC7Ctx *pr=malloc(nencoders*sizeof(SLIC7Ctx));
	int *rct_errors=nch>=3?(int*)malloc(dst->iw*sizeof(int[RCT_COUNT<<2])):0;
	if(!pr||!rct_errors)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int k=0;k<nencoders;++k)
	{
		int success=slic7_init(pr+k, dst->iw, dst->ih, nch, dst->depth, &ec);
		if(!success)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
	}
	RCTType rct_curr=RCT_SubGreen, rct0=RCT_SubGreen;
	if(nch>=3)
		memset(rct_errors, 0, dst->iw*sizeof(int[RCT_COUNT<<2]));
	//long long rct_errors[RCT_COUNT]={0};
	for(int ky=0, idx=0;ky<dst->ih;++ky)
	{
		//memset(rct_errors, 0, sizeof(rct_errors));
		for(int k=0;k<nencoders;++k)
			slic7_nextrow(pr+k, ky);
		for(int kx=0;kx<dst->iw;)
		{
			int kc=0;
			if(nch>=3)
			{
				int comp[3];
				if(!kx)
					rct_curr=rct0;
				comp[1]=slic7_dec(pr, nencoders, RCT_OVERRIDE, 0, kx);//Y
				comp[2]=slic7_dec(pr, nencoders, RCT_OVERRIDE, 1, kx);//U
				comp[0]=slic7_dec(pr, nencoders, RCT_OVERRIDE, 2, kx);//V
				rct_inv(comp, RCT_OVERRIDE);
				dst->data[idx|0]=comp[0];
				dst->data[idx|1]=comp[1];
				dst->data[idx|2]=comp[2];

				//if(!ky)
				//	printf("%3d  %-11s  RGB %d %d %d\n", kx, rct_names[rct_curr], comp[0], comp[1], comp[2]);//

				rct_curr=slic7_choose_rct(pr, nencoders, comp, rct_errors, kx, ky, dst->iw, dst->ih);
				if(!kx)
					rct0=rct_curr;
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
				int pixel=slic7_dec(pr, 1, RCT_NONE, kc, kx);
				dst->data[idx|kc]=pixel;
				++kc;
			}
			++kx;
			idx+=4;
		}
		if(loud&&!(ky&15))
			printf("%6.2lf%%  %10lf sec\r", 100.*(ky+1)/dst->ih, time_sec()-t_start);
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
	for(int k=0;k<nencoders;++k)
		slic7_free(pr+k);
	free(pr);
	return 1;
}