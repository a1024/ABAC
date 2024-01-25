#include"e2.h"
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;

//	#define PROFILER//SLOW
//	#define TRACK_SSE_RANGES//SLOW
//	#define DISABLE_SSE

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
void prof_print()
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
#define NPREDS 8
#define PRED_PREC 8
#define PARAM_PREC 8
//#define SSE_W 20
//#define SSE_H 20
//#define SSE_D 20
#define HIST_EXP 7
#define HIST_MSB 2
#define SSE_X_EXP 1
#define SSE_X_MSB 0
#define SSE_Y_EXP 1
#define SSE_Y_MSB 0
#define SSE_Z_EXP 1
#define SSE_Z_MSB 0
#define SSE_STAGES 10
#define SSE_FR_SIZE (1<<10)//separate final round
#define SSE_PREDBITS 5
#define SSE_PRED_LEVELS (1<<SSE_PREDBITS)//(_countof(qlevels_pred)+1)
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
int quantize_signed_get_range(int num, int den, int exp, int msb)
{
	int vmax=(int)((1LL<<24)*num/den>>16);
	int token=quantize_unsigned(vmax, exp, msb);
	token<<=1;
	return token;
}
int quantize_signed(int val, int shift, int exp, int msb, int nlevels)
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
#define QUANTIZE_HIST(X) quantize_unsigned(X, HIST_EXP, HIST_MSB)
static const int wp_params[]=
{
	0x0004AB65,//wp1
	0x0002841C,//wp2
	0x00034147,//wp3
	0x00016056,//wp4
	0x00016056,//av
	0x00016056,//grad
	0x00028F8A,//paper GAP
	0x00022F58,//CALIC GAP

	-0x005C,
	-0x005B,
	 0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,

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
	int sse_width, sse_height, sse_depth, sse_nplanes, sse_planesize, sse_size;
	int *pred_errors, *errors, *pixels;
	long long *hist, *histsums;
	unsigned *CDF_bypass, *CDFs;
	long long *sse, *sse_fr;
	int sse_idx[SSE_STAGES+1];
	long long sse_sum[SSE_STAGES+1];
	int sse_count[SSE_STAGES+1];
	int bias_count[4];
	long long bias_sum[4];
	int preds[NPREDS], params[_countof(wp_params)];
	long long pred, pred_sse[SSE_STAGES+1];
	int pred_final;
	int hist_idx;
	int sse_corr;
	int kc, kx, ky, kym0, kym1, kym2;
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
	pr->nhist=QUANTIZE_HIST(255);
	pr->cdfsize=hu.token+1;
	pr->sse_width=7;
	pr->sse_height=21;
	pr->sse_depth=21;
	//pr->sse_width=quantize_signed_get_range(2, 64, SSE_X_EXP, SSE_X_MSB);
	//pr->sse_height=quantize_signed_get_range(2, 1, SSE_Y_EXP, SSE_Y_MSB);
	//pr->sse_depth=quantize_signed_get_range(2, 1, SSE_Z_EXP, SSE_Z_MSB);
	//pr->sse_width=quantize_signed(255, 0, SSE_X_EXP, SSE_X_MSB, SSE_W);
	//pr->sse_height=quantize_signed(255, 0, SSE_Y_EXP, SSE_Y_MSB, SSE_H);
	//pr->sse_depth=quantize_signed(255, 0, SSE_Z_EXP, SSE_Z_MSB, SSE_D);
	pr->sse_nplanes=1<<SSE_PREDBITS;
	pr->sse_planesize=pr->sse_width*pr->sse_height*pr->sse_depth;
	pr->sse_size=pr->sse_nplanes*pr->sse_planesize;
	
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
	pr->hist=(long long*)malloc(pr->nhist*pr->cdfsize*sizeof(long long));//WH: cdfsize * NHIST
	pr->histsums=(long long*)malloc(pr->nhist*sizeof(long long));
	pr->CDF_bypass=(unsigned*)malloc((pr->cdfsize+1LL)*sizeof(int));
	pr->CDFs=(unsigned*)malloc(pr->nhist*(pr->cdfsize+1LL)*sizeof(int));
	pr->sse=(long long*)malloc(pr->sse_size*sizeof(long long[SSE_STAGES]));
	pr->sse_fr=(long long*)malloc(sizeof(long long[SSE_FR_SIZE]));
	if(!pr->pred_errors||!pr->errors||!pr->pixels||!pr->hist||!pr->histsums||!pr->CDFs||!pr->CDF_bypass||!pr->sse||!pr->sse_fr)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(pr->pred_errors, 0, (iw+PAD_SIZE*2LL)*sizeof(int[NPREDS*4*4]));
	memset(pr->errors, 0, (iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	memset(pr->pixels, 0, (iw+PAD_SIZE*2LL)*sizeof(int[4*4]));

	memset(pr->hist, 0, pr->nhist*pr->cdfsize*sizeof(long long));
	memset(pr->histsums, 0, pr->nhist*sizeof(long long));

	int sum=0;
	for(int ks=0;ks<pr->cdfsize;++ks)//TODO tune initial CDFs
		sum+=0x10000/(ks+1);
	for(int ks=0, c=0;ks<pr->cdfsize+1;++ks)
	{
		int freq=0x10000/(ks+1);
		pr->CDF_bypass[ks]=(int)(((long long)c<<16)/sum);
		c+=freq;
	}
	memfill(pr->CDFs, pr->CDF_bypass, pr->nhist*(pr->cdfsize+1LL)*sizeof(int), (pr->cdfsize+1LL)*sizeof(int));

	memset(pr->sse, 0, pr->sse_size*sizeof(long long[SSE_STAGES]));
	memset(pr->sse_fr, 0, sizeof(long long[SSE_FR_SIZE]));

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
	free(pr->sse);
	free(pr->hist);
	free(pr->CDFs);
}
static void slic5_update_CDFs(SLIC5Ctx *pr)
{
	for(int kt=0;kt<pr->nhist;++kt)//update CDFs
	{
		long long sum=pr->histsums[kt];
		if(sum>pr->cdfsize)//switch to histogram when it matures
		{
			long long *curr_hist=pr->hist+pr->cdfsize*kt;
			unsigned *curr_CDF=pr->CDFs+(pr->cdfsize+1)*kt;
			long long c=0;
			for(int ks=0;ks<pr->cdfsize;++ks)
			{
				long long freq=curr_hist[ks];
				curr_CDF[ks]=(int)(c*(0x10000LL-pr->cdfsize)/sum)+ks;
				c+=freq;
			}
			curr_CDF[pr->cdfsize]=0x10000;
		}
		//else//unknown statistics
		//	memcpy(curr_CDF, pr->CDF_bypass, (pr->cdfsize+1)*sizeof(int));
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
		NNWW    =LOAD(pr->pixels,  2, 2)<<PRED_PREC,
		NNW     =LOAD(pr->pixels,  1, 2)<<PRED_PREC,
		NN      =LOAD(pr->pixels,  0, 2)<<PRED_PREC,
		NNE     =LOAD(pr->pixels, -1, 2)<<PRED_PREC,
		NNEE    =LOAD(pr->pixels, -2, 2)<<PRED_PREC,
		NWW     =LOAD(pr->pixels,  2, 1)<<PRED_PREC,
		NW      =LOAD(pr->pixels,  1, 1)<<PRED_PREC,
		N       =LOAD(pr->pixels,  0, 1)<<PRED_PREC,
		NE      =LOAD(pr->pixels, -1, 1)<<PRED_PREC,
		NEE     =LOAD(pr->pixels, -2, 1)<<PRED_PREC,
		WW      =LOAD(pr->pixels,  2, 0)<<PRED_PREC,
		W       =LOAD(pr->pixels,  1, 0)<<PRED_PREC;
	long long
		eNNWW    =LOAD(pr->errors,  2, 2),//error = (curr<<8) - pred
		eNNW     =LOAD(pr->errors,  1, 2),
		eNN      =LOAD(pr->errors,  0, 2),
		eNNE     =LOAD(pr->errors, -1, 2),
		eNNEE    =LOAD(pr->errors, -2, 2),
		eNWW     =LOAD(pr->errors,  2, 1),
		eNW      =LOAD(pr->errors,  1, 1),
		eN       =LOAD(pr->errors,  0, 1),
		eNE      =LOAD(pr->errors, -1, 1),
		eNEE     =LOAD(pr->errors, -2, 1),
		eWW      =LOAD(pr->errors,  2, 0),
		eW       =LOAD(pr->errors,  1, 0);
	int sh=pr->shift_prec[kc];
	int clamp_lo=(int)N, clamp_hi=(int)N;
	clamp_lo=(int)MINVAR(clamp_lo, W);
	clamp_hi=(int)MAXVAR(clamp_hi, W);
	clamp_lo=(int)MINVAR(clamp_lo, NE);
	clamp_hi=(int)MAXVAR(clamp_hi, NE);
	PROF(FETCH_NB);

	int j=-1;
	pr->preds[++j]=(int)(W+NE-N);
	pr->preds[++j]=(int)(N-((eN+eW+eNE)*pr->params[NPREDS+0]>>PARAM_PREC));
	pr->preds[++j]=(int)(W-((eN+eW+eNW)*pr->params[NPREDS+1]>>PARAM_PREC));
	pr->preds[++j]=(int)(N-((
		eNW*pr->params[NPREDS+2]+
		eN*pr->params[NPREDS+3]+
		eNE*pr->params[NPREDS+4]+
		((long long)NN-N)*pr->params[NPREDS+5]+
		((long long)NW-W)*pr->params[NPREDS+6]
	)>>PARAM_PREC));

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
#if 1
	//++j, pr->preds[j]=(int)(N+W+NW+NE+((eN+eW+eNW+eNE)>>2))>>2;//v
	++j, pr->preds[j]=(int)(N+W+NW+NE+((eN+eW+eNW+eNE)>>pr->shift_error[kc]))>>2;
	//++j, pr->preds[j]=(NNWW+NNW+NN+NNE+NNEE+NWW+NW+N+NE+NEE+WW+W)/12;//~
	//++j, pr->preds[j]=(int)(
	//	9*(N+W+((eN+eW)>>(1+(pr->depths[kc]<=9))))+
	//	7*(NW+NE+((eNW+eNE)>>(1+(pr->depths[kc]<=9))))
	//)>>5;
	//++j, pr->preds[j]=(4*(N+W)+2*(NW+NE)+eN+eW+eNW+eNE)>>3;//X
	//++j, pr->preds[j]=2*W-WW+(eW>>1), pr->preds[j]=CLAMP(pr->min_allowed, pr->preds[j], pr->max_allowed);//X

	++j, pr->preds[j]=(int)(N+W-NW);//, pr->preds[j]=(int)MEDIAN3(N, W, pr->preds[j]);

	int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
	int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
	int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
	int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
	++j;
	int sum2=(dy+dx)>>sh, diff=(dy-dx)>>sh, diff2=(NE-NW)>>3;
	if(sum2>(32<<(PRED_PREC+1)))
		pr->preds[j]=(dx*N+dy*W)/sum2+diff2;
	else if(diff>(12<<(PRED_PREC+1)))
		pr->preds[j]=(N+2*W)/3+diff2;
	else if(diff<-(12<<(PRED_PREC+1)))
		pr->preds[j]=(2*N+W)/3+diff2;
	else
		pr->preds[j]=((N+W)>>1)+diff2;
	diff=d45-d135;
	if(diff>(32<<(PRED_PREC+1)))
		pr->preds[j]+=diff2;
	else if(diff>(16<<(PRED_PREC+1)))
		pr->preds[j]+=diff2>>1;
	else if(diff<-(32<<(PRED_PREC+1)))
		pr->preds[j]-=diff2;
	else if(diff<-(16<<(PRED_PREC+1)))
		pr->preds[j]-=diff2>>1;

	++j;
	int disc=(dy-dx)>>sh;
	if(disc>80)
		pr->preds[j]=W;
	else if(disc<-80)
		pr->preds[j]=N;
	else
	{
		pr->preds[j]=(((N+W)<<1)+(NE-NW))>>2;
		if(disc>32)
			pr->preds[j]=(pr->preds[j]+W)>>1;
		else if(disc>8)
			pr->preds[j]=(3*pr->preds[j]+W)>>2;
		else if(disc<-32)
			pr->preds[j]=(pr->preds[j]+N)>>1;
		else if(disc<-8)
			pr->preds[j]=(3*pr->preds[j]+N)>>2;
	}
#endif
	PROF(CALC_SUBPREDS);

	long long weights[NPREDS]={0}, wsum=0;
	for(int k=0;k<NPREDS;++k)
	{
		weights[k]=
			(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx+2+PAD_SIZE)+k)<<2|kc]+//peNEE [not in jxl]
			(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx+1+PAD_SIZE)+k)<<2|kc]+//peNE
			(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx  +PAD_SIZE)+k)<<2|kc]+//peN+peW
			(long long)pr->pred_errors[(NPREDS*(pr->kym1+kx-1+PAD_SIZE)+k)<<2|kc];//peNW+peWW
		weights[k]=((long long)pr->params[k]<<8)/(weights[k]+1);
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
	//if(kx==10&&ky==10)
	//	printf("");
	
	pr->hist_idx=MAXVAR(dx, dy)>>(sh-2);
	pr->hist_idx=QUANTIZE_HIST(pr->hist_idx);
	//if(pr->hist_idx==32)//
	//	printf("");
	//int
	//	qx  =QUANTIZE_HIST(dx  >>(sh-2)),//dx>>shift_pred is in [0 ~ 255*3],  dx>>(shift_pred-2) has quad range
	//	qy  =QUANTIZE_HIST(dy  >>(sh-2));
	//	q45 =QUANTIZE_HIST(d45 >>(sh-2)),
	//	q135=QUANTIZE_HIST(d135>>(sh-2));
	//pr->hist_idx=MAXVAR(qx, qy);
	//if(pr->hist_idx>pr->nhist-1)//
	//	LOG_ERROR("");
	pr->hist_idx=MINVAR(pr->hist_idx, pr->nhist-1);
#ifndef DISABLE_SSE
	int g[]=
	{
		//4D SSE
#define QUANTIZE(X, Y, Z)\
	pr->sse_width*(pr->sse_height*quantize_signed(Z, sh, SSE_Z_EXP, SSE_Z_MSB, pr->sse_depth)+\
	quantize_signed(Y, sh, SSE_Y_EXP, SSE_Y_MSB, pr->sse_height))+\
	quantize_signed(X, sh, SSE_X_EXP, SSE_X_MSB, pr->sse_width)

		QUANTIZE((NNE-NN)>>6,			(N-NW)>>4,		W-WW),//5
		QUANTIZE((NE-NNE)>>6,			(N-NN)>>4,		W-NW),//6
		QUANTIZE((int)eNW>>6,			(int)eW<<1,		(int)eN<<1),//9
		QUANTIZE((int)eNE>>6,			(int)eWW<<1,		(int)eNN<<1),//10
		QUANTIZE((3*(NE-NW)+NNWW-NNEE)>>9,	(W-WW)>>3,		(NW+2*NE+WW-4*W)>>2),//4
		QUANTIZE((int)(eNE-eNNE)>>6,		(int)(eN-eNN)>>4,	(int)(eW-eNW)),//8
		QUANTIZE((int)(eNNE-eNN)>>6,		(int)(eN-eNW)>>4,	(int)(eW-eWW)),//7
		QUANTIZE((NW+3*(NE-N)-NEE)>>9,		(NE-NNEE)>>3,		(NNE+NEE+2*W-4*NE)>>2),//1
		QUANTIZE((NWW+3*(N-NW)-NE)>>9,		(N-NN)>>3,		(W+NW+NE+NN-4*N)>>2),//2
		QUANTIZE((3*(N-W)+WW-NN)>>9,		(NW-NNWW)>>3,		(NNW+NWW+W+N-4*NW)>>2),//3 (weakest)
#if 0
		QUANTIZE((NNE-NN)/32,			(N-NW)/8,		(W-WW)/2),//5
		QUANTIZE((NE-NNE)/32,			(N-NN)/8,		(W-NW)/2),//6
		//	QUANTIZE((int)(eNNE-eNN)/32,	(int)(eN-eNW)/2,	(int)(eW-eWW)/2),//7
		QUANTIZE((int)eNW/16,			(int)eW*4,		(int)eN*4),//9
		QUANTIZE((int)eNE/16,			(int)eWW*2,		(int)eNN*2),//10
		//	QUANTIZE(NW+3*(NE-N)-NEE,	NE-NNEE,		NNE+NEE+2*W-4*NE),//1
		//	QUANTIZE(NWW+3*(N-NW)-NE,	N-NN,			W+NW+NE+NN-4*N),//2
		//	QUANTIZE((3*(N-W)+WW-NN)/128,	(NW-NNWW)/2,		(NNW+NWW+W+N-4*NW)/16),//3 X
		QUANTIZE((3*(NE-NW)+NNWW-NNEE)/128,	(W-WW)/8,		(NW+2*NE+WW-4*W)/16),//4
		QUANTIZE((int)(eNE-eNNE)/32,		(int)(eN-eNN)/8,	(int)(eW-eNW)/2),//8
#endif
#undef  QUANTIZE

#if 0
		pr->sse_width*(pr->sse_height*quantize_signed((int)(W-WW)/2, sh, 2, 1, pr->sse_depth)+quantize_signed((int)(N-NW)*2, sh, 2, 1, pr->sse_height))+quantize_signed((int)(NNE-NN)/64, sh, 2, 1, pr->sse_width),//gx
		pr->sse_width*(pr->sse_height*quantize_signed((int)(W-NW)/2, sh, 2, 1, pr->sse_depth)+quantize_signed((int)(N-NN)*2, sh, 2, 1, pr->sse_height))+quantize_signed((int)(NE-NNE)/64, sh, 2, 1, pr->sse_width),//gy
		pr->sse_width*(pr->sse_height*quantize_signed((int)(eW-eWW)/2, sh, 2, 1, pr->sse_depth)+quantize_signed((int)(eN-eNW)*2, sh, 2, 1, pr->sse_height))+quantize_signed((int)(eNNE-eNN)/64, sh, 2, 1, pr->sse_width),//egx
		pr->sse_width*(pr->sse_height*quantize_signed((int)(eW-eNW)/2, sh, 2, 1, pr->sse_depth)+quantize_signed((int)(eN-eNN)*2, sh, 2, 1, pr->sse_height))+quantize_signed((int)(eNE-eNNE)/64, sh, 2, 1, pr->sse_width),//egy
		pr->sse_width*(pr->sse_height*quantize_signed((int) eN*4 , sh, 2, 1, pr->sse_depth)+quantize_signed((int) eW*2   , sh, 2, 1, pr->sse_height))+quantize_signed((int) eNW/32    , sh, 2, 1, pr->sse_width),//eclose
		pr->sse_width*(pr->sse_height*quantize_signed((int) eNN*4, sh, 2, 1, pr->sse_depth)+quantize_signed((int) eWW*2  , sh, 2, 1, pr->sse_height))+quantize_signed((int) eNE/32    , sh, 2, 1, pr->sse_width),//efar

		//SSE_W*(SSE_H*quantize_signed((int)(W-WW), sh+1, 4, 1, SSE_D)+quantize_signed((int)(N-NW), sh+1, 4, 1, SSE_H))+quantize_signed((int)(NNE-NN), sh+1, 4, 1, SSE_W),//gx		X  bad
		//SSE_W*(SSE_H*quantize_signed((int)(W-NW), sh+1, 4, 1, SSE_D)+quantize_signed((int)(N-NN), sh+1, 4, 1, SSE_H))+quantize_signed((int)(NE-NNE), sh+1, 4, 1, SSE_W),//gy
		//SSE_W*(SSE_H*quantize_signed((int) eN   , sh, 4, 1, SSE_D)+quantize_signed((int) eW   , sh, 4, 1, SSE_H))+quantize_signed((int) eNW    , sh, 4, 1, SSE_W),//eclose
		//SSE_W*(SSE_H*quantize_signed((int) eNN  , sh, 4, 1, SSE_D)+quantize_signed((int) eWW  , sh, 4, 1, SSE_H))+quantize_signed((int) eNE    , sh, 4, 1, SSE_W),//efar

	//	SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)(W-WW)>>shifts[0])+QUANTIZE_SSE((int)(N-NW)>>shifts[1]))+QUANTIZE_SSE((int)(NNE-NN)>>shifts[2]),//gx
	//	SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)(W-NW)>>shifts[0])+QUANTIZE_SSE((int)(N-NN)>>shifts[1]))+QUANTIZE_SSE((int)(NE-NNE)>>shifts[2]),//gy
	//	//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)N>>shifts[3])+QUANTIZE_SSE((int)W>>shifts[4]))+QUANTIZE_SSE((int)NW>>shifts[5]),//close
	//	SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)eN>>shifts[3])+QUANTIZE_SSE((int)eW>>shifts[4]))+QUANTIZE_SSE((int)eNW>>shifts[5]),//eclose
	//	SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)eNN>>shifts[3])+QUANTIZE_SSE((int)eWW>>shifts[4]))+QUANTIZE_SSE((int)eNE>>shifts[5]),//efar
		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)NN>>shifts[3])+QUANTIZE_SSE((int)WW>>shifts[4]))+QUANTIZE_SSE((int)NE>>shifts[5]),//far
		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)(eNW-eNE)>>shifts[3])+QUANTIZE_SSE((int)(eNN-eWW)>>shifts[4]))+QUANTIZE_SSE((int)(eN-eW)>>shifts[5]),//ediff
		
		//+-+-+		+-+-+		+++++
		//-+-+-		+-+--		-----
		//+-?		+-		+-
	//	SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)(NNWW-NNW+NN-NNE+NNEE - NWW+NW-N+NE-NEE + WW-W)/12>>shifts[0])+
	//		QUANTIZE_SSE((int)(NNWW-NNW+NN-NNE+NNEE + NWW-NW+N-NE-NEE + WW-W)/12>>shifts[1]))+
	//		QUANTIZE_SSE((int)(NNWW+NNW+NN+NNE+NNEE - NWW-NW-N-NE-NEE + WW-W)/12>>shifts[2]),

		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)WW>>pr->shift_prec[kc])+QUANTIZE_SSE((int)W>>pr->shift_prec[kc]))+QUANTIZE_SSE((int)NWW>>pr->shift_prec[kc]),
		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)NW>>pr->shift_prec[kc])+QUANTIZE_SSE((int)N>>pr->shift_prec[kc]))+QUANTIZE_SSE((int)NE>>pr->shift_prec[kc]),
		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)NEE>>pr->shift_prec[kc])+QUANTIZE_SSE((int)NNWW>>pr->shift_prec[kc]))+QUANTIZE_SSE((int)NNW>>pr->shift_prec[kc]),
		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)NN>>pr->shift_prec[kc])+QUANTIZE_SSE((int)NNE>>pr->shift_prec[kc]))+QUANTIZE_SSE((int)NNEE>>pr->shift_prec[kc]),
		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)eNN>>pr->shift_prec[kc])+QUANTIZE_SSE((int)eNW>>pr->shift_prec[kc]))+QUANTIZE_SSE((int)eN>>pr->shift_prec[kc]),
		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)eNE>>pr->shift_prec[kc])+QUANTIZE_SSE((int)eWW>>pr->shift_prec[kc]))+QUANTIZE_SSE((int)eW>>pr->shift_prec[kc]),

		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)(W-WW)>>shifts[0])+QUANTIZE_SSE((int)(N-NW)>>shifts[1]))+QUANTIZE_SSE((int)(NNE-NN)>>shifts[2]),//gx
		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)(W-NW)>>shifts[0])+QUANTIZE_SSE((int)(N-NN)>>shifts[1]))+QUANTIZE_SSE((int)(NE-NNE)>>shifts[2]),//gy
		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)eN>>shifts[3])+QUANTIZE_SSE((int)eW>>shifts[4]))+QUANTIZE_SSE((int)eNW>>shifts[5]),//eclose
		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)eNN>>shifts[3])+QUANTIZE_SSE((int)eWW>>shifts[4]))+QUANTIZE_SSE((int)eNE>>shifts[5]),//efar

		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)(W-WW)>>(pr->shift_prec[kc]+4))+QUANTIZE_SSE((int)(N-NW)>>(pr->shift_prec[kc]+4)))+QUANTIZE_SSE((int)(NNE-NN)>>(pr->shift_prec[kc]+4)),
		//SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)(W-NW)>>(pr->shift_prec[kc]+4))+QUANTIZE_SSE((int)(N-NN)>>(pr->shift_prec[kc]+4)))+QUANTIZE_SSE((int)(NE-NNE)>>(pr->shift_prec[kc]+4)),
		//SSE_WIDTH*(QUANTIZE_SSE((int)eN>>(pr->shift_prec[kc]+3))*QUANTIZE_SSE((int)eW>>(pr->shift_prec[kc]+3)))+QUANTIZE_SSE((int)eNW>>(pr->shift_prec[kc]+3)),
		//SSE_WIDTH*(QUANTIZE_SSE((int)eNN>>(pr->shift_prec[kc]+3))*QUANTIZE_SSE((int)eWW>>(pr->shift_prec[kc]+3)))+QUANTIZE_SSE((int)eNE>>(pr->shift_prec[kc]+3)),

		//QUANTIZE_SSE((int)(W-WW)>>(pr->shift_prec[kc]+4))*QUANTIZE_SSE((int)(N-NW)>>(pr->shift_prec[kc]+4))*QUANTIZE_SSE((int)(NNE-NN)>>(pr->shift_prec[kc]+4)),//X  bad
		//QUANTIZE_SSE((int)(W-NW)>>(pr->shift_prec[kc]+4))*QUANTIZE_SSE((int)(N-NN)>>(pr->shift_prec[kc]+4))*QUANTIZE_SSE((int)(NE-NNE)>>(pr->shift_prec[kc]+4)),
		//QUANTIZE_SSE((int)eN>>(pr->shift_prec[kc]+3))*QUANTIZE_SSE((int)eW>>(pr->shift_prec[kc]+3))*QUANTIZE_SSE((int)eNW>>(pr->shift_prec[kc]+3)),
		//QUANTIZE_SSE((int)eNN>>(pr->shift_prec[kc]+3))*QUANTIZE_SSE((int)eWW>>(pr->shift_prec[kc]+3))*QUANTIZE_SSE((int)eNE>>(pr->shift_prec[kc]+3)),
#endif

		//3D SSE
#if 0
		QUANTIZE_SSE((int)(W-WW+NE-NW+NN-NNE)>>(pr->shift_prec[kc]+5))*QUANTIZE_SSE((int)(W-NW+N-NN+NE-NNE)>>(pr->shift_prec[kc]+5)),

		//QUANTIZE_SSE((int)(W-WW)>>(pr->shift_prec[kc]+4))*QUANTIZE_SSE((int)(W-NW)>>(pr->shift_prec[kc]+4)),//X  bad
		//QUANTIZE_SSE((int)(NE-NW)>>(pr->shift_prec[kc]+4))*QUANTIZE_SSE((int)(N-NN)>>(pr->shift_prec[kc]+4)),
		//QUANTIZE_SSE((int)(NN-NNE)>>(pr->shift_prec[kc]+4))*QUANTIZE_SSE((int)(NE-NNE)>>(pr->shift_prec[kc]+4)),
		// 
		//SSE_WIDTH*QUANTIZE_SSE((int)(W-WW+NE-NW+NN-NNE)>>(pr->shift_prec[kc]+5))+QUANTIZE_SSE((int)(W-NW+N-NN+NE-NNE)>>(pr->shift_prec[kc]+5)),//original
		//SSE_WIDTH*QUANTIZE_SSE((int)(W-WW+NE-N)>>(pr->shift_prec[kc]+4))+QUANTIZE_SSE((int)(W-NW+NE-NNE)>>(pr->shift_prec[kc]+4)),//X bad
		//SSE_WIDTH*QUANTIZE_SSE((int)(W-WW+N-NW+NNE-NN)>>(pr->shift_prec[kc]+5))+QUANTIZE_SSE((int)(W-NW+N-NN+NE-NNE)>>(pr->shift_prec[kc]+5)),//X  bad
		//SSE_WIDTH*QUANTIZE_SSE((int)(W-WW+NE-NW+NN-NNE)>>(pr->shift_prec[kc]+5))+QUANTIZE_SSE((int)(W-NNW+N-NN+NNE-NE)>>(pr->shift_prec[kc]+5)),//X  bad
		SSE_WIDTH*QUANTIZE_SSE((int)eN >>(pr->shift_prec[kc]+3))+QUANTIZE_SSE((int)eW >>(pr->shift_prec[kc]+3)),
		SSE_WIDTH*QUANTIZE_SSE((int)eNW>>(pr->shift_prec[kc]+3))+QUANTIZE_SSE((int)eNE>>(pr->shift_prec[kc]+3)),
		SSE_WIDTH*QUANTIZE_SSE((int)eNN>>(pr->shift_prec[kc]+3))+QUANTIZE_SSE((int)eWW>>(pr->shift_prec[kc]+3)),

		//SSE_WIDTH*QUANTIZE_SSE((int)(W-WW+NE-NW+NN-NNE)>>(pr->shift_prec[kc]+5))+QUANTIZE_SSE((int)eW >>(pr->shift_prec[kc]+3)),//X  bad
		//SSE_WIDTH*QUANTIZE_SSE((int)(W-NW+N -NN+NE-NNE)>>(pr->shift_prec[kc]+5))+QUANTIZE_SSE((int)eN >>(pr->shift_prec[kc]+3)),
		//SSE_WIDTH*QUANTIZE_SSE((int)eNW>>(pr->shift_prec[kc]+3))+QUANTIZE_SSE((int)eNE>>(pr->shift_prec[kc]+3)),
		//SSE_WIDTH*QUANTIZE_SSE((int)eNN>>(pr->shift_prec[kc]+3))+QUANTIZE_SSE((int)eWW>>(pr->shift_prec[kc]+3)),
#endif

		//2D SSE
#if 0
		QUANTIZE_SSE((int)(W-WW+NE-NW+NN-NNE)>>(pr->shift_prec[kc]+5)),//gx
		QUANTIZE_SSE((int)(W-NW+N -NN+NE-NNE)>>(pr->shift_prec[kc]+5)),//gy
		QUANTIZE_SSE((int)eW >>(pr->shift_prec[kc]+3)),
		QUANTIZE_SSE((int)eNW>>(pr->shift_prec[kc]+3)),
		QUANTIZE_SSE((int)eN >>(pr->shift_prec[kc]+3)),
		QUANTIZE_SSE((int)eNE>>(pr->shift_prec[kc]+3)),
		QUANTIZE_SSE((int)eWW>>(pr->shift_prec[kc]+3)),
		QUANTIZE_SSE((int)eNN>>(pr->shift_prec[kc]+3)),
#endif
	};
#if 0
	int g[SSE_STAGES];
	g[0]=QUANTIZE_SSE((int)(W-WW+NE-NW+NN-NNE)>>(pr->shift_prec[kc]+5));//gx
	g[1]=QUANTIZE_SSE((int)(W-NW+N -NN+NE-NNE)>>(pr->shift_prec[kc]+5));//gy
	g[2]=QUANTIZE_SSE((int)eW >>(pr->shift_prec[kc]+3));
	g[3]=QUANTIZE_SSE((int)eNW>>(pr->shift_prec[kc]+3));
	g[4]=QUANTIZE_SSE((int)eN >>(pr->shift_prec[kc]+3));
	g[5]=QUANTIZE_SSE((int)eNE>>(pr->shift_prec[kc]+3));
	g[6]=QUANTIZE_SSE((int)eWW>>(pr->shift_prec[kc]+3));
	g[7]=QUANTIZE_SSE((int)eNN>>(pr->shift_prec[kc]+3));
#endif
#if 0
	int g[]=
	{
		QUANTIZE_SSE((int)(W-WW+NE-NW+NN-NNE)>>(pr->shift_prec[kc]+5)),//d>>shift in [-768, 768]
		QUANTIZE_SSE((int)(W-NW+N -NN+NE-NNE)>>(pr->shift_prec[kc]+5)),
	//	qx,
	//	qy,
	//	QUANTIZE_SSE((int)(dx+W-WW+NE-NW+NN-NNE)>>(pr->shift[kc]+3)),
	//	QUANTIZE_SSE((int)(dy+W-NW+N -NN+NE-NNE)>>(pr->shift[kc]+3)),
	//	q45,
	//	q135,
		QUANTIZE_SSE((int)eW >>(pr->shift_prec[kc]+3)),
		QUANTIZE_SSE((int)eNW>>(pr->shift_prec[kc]+3)),
		QUANTIZE_SSE((int)eN >>(pr->shift_prec[kc]+3)),
		QUANTIZE_SSE((int)eNE>>(pr->shift_prec[kc]+3)),
		QUANTIZE_SSE((int)eWW>>(pr->shift_prec[kc]+3)),
		QUANTIZE_SSE((int)eNN>>(pr->shift_prec[kc]+3)),

	//	QUANTIZE_SSE((int)N  >>pr->shift[kc]),
	//	QUANTIZE_SSE((int)W  >>pr->shift[kc]),
	};
#endif
	PROF(CALC_CTX);
	//if(kx==10&&ky==10)//
	//	printf("");
	pr->sse_corr=0;
	for(int k=0;k<SSE_STAGES+1;++k)
	{
		long long sse_val;
		if(k<SSE_STAGES)
		{
			int qp=(int)(pr->pred>>(pr->shift_prec[kc]+8-(SSE_PREDBITS-1)));
			qp+=(1<<SSE_PREDBITS)>>1;
			qp=CLAMP(0, qp, (1<<SSE_PREDBITS)-1);
			pr->sse_idx[k]=pr->sse_planesize*qp+g[k];
			sse_val=pr->sse[pr->sse_size*k+pr->sse_idx[k]];
			
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
			pr->sse_idx[k]=
				//(W<NW)<<11|
				//(NW<NE)<<10|
				(N<W)<<9|
				(N+W-NW<(int)pr->pred)<<8|
				(W <pr->pred)<<7|
				(NW<pr->pred)<<6|
				(N <pr->pred)<<5|
				(NE<pr->pred)<<4|
				(WW<pr->pred)<<3|
				(NN<pr->pred)<<2|
				(2*N-NN<(int)pr->pred)<<1|
				(2*W-WW<(int)pr->pred)
			;
			sse_val=pr->sse_fr[pr->sse_idx[k]];
		}
		pr->sse_count[k]=(int)(sse_val&0xFFF);
		pr->sse_sum[k]=(int)(sse_val>>12);
		//int sse_corr=(int)((pr->sse_sum[k]+(pr->sse_count[k]>>1))/(pr->sse_count[k]+1));
		int sse_corr=(int)(pr->sse_sum[k]/(pr->sse_count[k]+1));
		//int sse_corr=pr->sse_count[k]?(int)(pr->sse_sum[k]/pr->sse_count[k]):0;//X
		pr->pred+=sse_corr;
		pr->sse_corr+=sse_corr;

		pr->pred_sse[k]=pr->pred;
	}
#endif
	//int final_corr=(int)((pr->bias_sum[kc]+(pr->bias_count[kc]>>1))/(pr->bias_count[kc]+1));
	int final_corr=(int)(pr->bias_sum[kc]/(pr->bias_count[kc]+1));
	//int final_corr=pr->bias_count[kc]?(int)(pr->bias_sum[kc]/pr->bias_count[kc]):0;//X
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
	++pr->hist[pr->cdfsize*pr->hist_idx+token];
	++pr->histsums[pr->hist_idx];
	if(pr->histsums[pr->hist_idx]>((long long)pr->iw*pr->nch<<1))//rescale hist
	{
		long long *curr_hist=pr->hist+pr->cdfsize*pr->hist_idx;
		long long sum=0;
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
			pr->sse[pr->sse_size*k+pr->sse_idx[k]]=sse_val;
		else
			pr->sse_fr[pr->sse_idx[k]]=sse_val;
	}
	++pr->bias_count[kc];
	pr->bias_sum[kc]+=error;
	PROF(UPDATE_SSE);
#endif
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
		error=upred+abs_error;
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

	ac_enc(pr->ec, hu.token, pr->CDFs+(pr->cdfsize+1)*pr->hist_idx, pr->cdfsize, 1);
	if(hu.nbits)
	{
		int bypass=hu.bypass, nbits=hu.nbits;
		while(nbits>8)
		{
			ac_enc(pr->ec, bypass>>(nbits-8)&0xFF, 0, 1<<8, 0x10000>>8);
			nbits-=8;
		}
		ac_enc(pr->ec, bypass&((1<<nbits)-1), 0, 1<<nbits, 0x10000>>nbits);
	}

	slic5_update(pr, curr, hu.token);
}
static int slic5_dec(SLIC5Ctx *pr, int kc, int kx, int ky)
{
	slic5_predict(pr, kc, kx, ky);

	int token=ac_dec(pr->ec, pr->CDFs+(pr->cdfsize+1)*pr->hist_idx, pr->cdfsize, 1);
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
			bypass|=ac_dec(pr->ec, 0, 1<<8, 0x10000>>8)<<n;
		}
		bypass|=ac_dec(pr->ec, 0, 1<<n, 0x10000>>n);
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
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		slic5_nextrow(&pr, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			int kc=0;
			if(src->nch>=3)
			{
				int r=src->data[idx|0], g=src->data[idx|1], b=src->data[idx|2];

				r-=g;
				b-=g;
				g+=(r+b)>>2;

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
		printf("\n");
		printf("csize %8d  CR %10.6lf\n", (int)list.nobj, usize/list.nobj);

		printf("Hist WH %d*%d:\n", pr.cdfsize, pr.nhist);
		for(int kt=0;kt<pr.nhist;++kt)
		{
			long long *curr_hist=pr.hist+pr.cdfsize*kt;
			for(int ks=0;ks<pr.cdfsize;++ks)
				printf("%lld%c", curr_hist[ks], ks<pr.cdfsize-1?' ':'\n');
		}
		
#ifdef TRACK_SSE_RANGES
		const char labels[]="XYZP";
		printf("SSE ranges\n");
		for(int k=0;k<_countof(pr.sse_ranges)-1;k+=2)
		{
			int n=(&pr.sse_width)[k>>1];
			printf("  %c %4d ~ %4d / %4d  ", labels[k>>1], pr.sse_ranges[k], pr.sse_ranges[k+1], n);
			for(int k2=0;k2<n;++k2)
				printf("%c", k2>=pr.sse_ranges[k]&&k2<=pr.sse_ranges[k+1]?'*':'-');
			printf("\n");
		}
#endif

		printf("Bias\n");
		for(int kc=0;kc<src->nch;++kc)
			printf("  %lld/%d = %d\n", pr.bias_sum[kc], pr.bias_count[kc], pr.bias_count[kc]?(int)(pr.bias_sum[kc]/pr.bias_count[kc]):0);
		prof_print();
	}
	slic5_free(&pr);
	dlist_clear(&list);
	return 1;
}
int t47_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	double t_start=time_sec();
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
			int kc=0;
			if(dst->nch>=3)
			{
				int g=slic5_dec(&pr, 0, kx, ky);//Y
				int b=slic5_dec(&pr, 1, kx, ky);//Cb
				int r=slic5_dec(&pr, 2, kx, ky);//Cr

				g-=(r+b)>>2;
				b+=g;
				r+=g;
				//if(guide&&(r!=guide->data[(guide->iw*ky+kx)<<2|0]||g!=guide->data[(guide->iw*ky+kx)<<2|1]||b!=guide->data[(guide->iw*ky+kx)<<2|2]))
				//{
				//	r-=g;
				//	b-=g;
				//	g+=(r+b)>>2;
				//	LOG_ERROR("Guide error XY %d %d", kx, ky);
				//}

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
	return 1;
}