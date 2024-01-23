#include"e2.h"
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
static const char file[]=__FILE__;

	#define SLIC5_DISABLE_PALETTE//do NOT enable palette, currently not supported, will leak memory
	#define PRINT_SSE

#define SLIC5_CONFIG_EXP 5
#define SLIC5_CONFIG_MSB 2
#define SLIC5_CONFIG_LSB 0

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

static const int qlevels_hist[]=
{
	0, 1, 2, 3, 5, 7, 11, 13, 15, 23, 27, 31, 47, 55, 63, 95, 127, 191, 255, 392, 511, 767, 1023,

	//1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 15, 19, 23, 27, 31, 39, 47, 55, 63, 79, 95, 111, 127, 159, 191, 223, 255, 323, 392, 446, 511, 639, 767, 895, 1023,
	//1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 15, 19, 23, 27, 31, 39, 47, 55, 63, 79, 95, 111, 127, 159, 191, 223, 255,
//	1, 3, 5, 7, 11, 15, 23, 31, 47, 63, 95, 127, 191, 255, 392, 511, 767, 1023,
	//1, 3, 5, 7, 11, 15, 23, 31, 47, 63, 95, 127, 191, 255, 392,
	//3, 7, 15, 31, 63, 127, 255, 511, 1023,

	//1, 3, 5, 7, 11, 15, 23, 31, 47, 63, 95, 127, 191, 255, 392, 500,
	//1, 3, 5, 7, 11, 15, 23, 31, 47, 63, 95, 127, 255, 512, 1023,
};
static const int qlevels_sse[]=//symmetric to avoid quantized negative zero
{
//	-22, -14, -9, -7, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 7, 9, 14, 22,
	-23, -15, -11, -7, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 7, 11, 15, 23,
//	-31, -23, -15, -11, -7, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 7, 11, 15, 23, 31,
	//-31, -27, -23, -19, -15, -13, -11, -9, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 15, 19, 23, 27, 31,
	//-31, -29, -27, -25, -23, -21, -19, -17, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 21, 23, 25, 27, 29, 31,
};
static const int qlevels_pred[]=
{
	//-1023,-895,-767,-639,-511,-446,-392,-323,-255,-223,
	//-191, -159, -127, -111, -95, -79, -63, -55,
	-47, -39, -31, -27, -23, -19, -15, -13, -11, -9, -7, -6, -5, -4, -3, -2, -1,
	0,
	1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 15, 19, 23, 27, 31, 39, 47,
	//55, 63, 79, 95, 111, 127, 159, 191,
	//223, 255, 323, 392, 446, 511, 639, 767, 895, 1023,
};
#define CDF_UPDATE_PERIOD 0x400//number of sub-pixels processed between CDF updates, must be a power-of-two
#define PAD_SIZE 2
#define NPREDS 10
#define PRED_PREC 8
#define PARAM_PREC 8
#define SSE_STAGES 4
#define SSE_FR_SIZE (1<<10)//separate final round
#define NHIST (_countof(qlevels_hist)+1)
#define SSE_WIDTH (_countof(qlevels_sse)+1)
#define SSE_HEIGHT (_countof(qlevels_pred)+1)
#define SSE_SIZE (SSE_WIDTH*SSE_WIDTH*SSE_WIDTH*SSE_HEIGHT)
//#define HASH_FUNC(A, B) ((A+B+1)*(A+B)/2+B)
static int quantize(int x, const int *qlevels, int nlevels)//experimental, bottleneck, should be turned into an unrolled macro on final release
{
	//if(qlevels==qlevels_hist&&(x==1||x==3))//NOT HIT
	//if(qlevels==qlevels_hist&&(x==0||x==2))
	//	printf("");
	int L=0, R=nlevels-1;
	while(L<=R)
	{
		int mid=(L+R)>>1;
		if(x>qlevels[mid])
			L=mid+1;
		else if(x<qlevels[mid])
			R=mid-1;
		else
			return mid;
	}
	L+=L<nlevels-1&&x>qlevels[L];
	return L;
}
//static int QUANTIZE_SSE(int x)
//{
//	x+=SSE_WIDTH/2;
//	x=CLAMP(0, x, SSE_WIDTH-1);
//	return x;
//}
#define QUANTIZE_HIST(X) quantize(X, qlevels_hist, _countof(qlevels_hist))
#define QUANTIZE_SSE(X) quantize(X, qlevels_sse, _countof(qlevels_sse))
#define QUANTIZE_PRED(X) quantize(X, qlevels_pred, _countof(qlevels_pred))
static const int wp_params[]=
{
	// 0x00EA, 0x01C8, 0x00A2, 0x005E, 0x01F4, 0x0045, 0x0091, 0x0066, 0x003B, 0x0027,-0x0011, 0x001B, 0x00FF, 0x007E, 0x00D1, 0x00F3, 0x008F, 0x0130, 0x018E,-0x00AC, 0x0004,
	// 0x010C, 0x0008,-0x007E, 0x00A2, 0x000E,-0x0069,-0x0073,-0x0125,-0x0092, 0x0000, 0x0078,

	0x0004AB65,//wp1
	0x0002841C,//wp2
	0x00034147,//wp3
	0x00016056,//wp4
	0x00016056,//wp5
	0x00016056,//wp6
	0x00016056,//wp7
	0x00016056,//clamp grad
	0x00028F8A,//journal GAP
	0x00022F58,//CALIC GAP

//	0xD000, 0xC000, 0xC000, 0xB000,//wp weights
//	0xC000,//journal GAP
//	0xD000,//CALIC GAP

//	 0x0DB8,  0x0E22,  0x181F,  0x0BF3,//wp weights
//	//0x181F,//clamp grad
//	 0x181F,//journal GAP
//	 0x1A1F,//CALIC GAP

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
	//	nlevels[4],
		half[4],
	//	min_allowed[4], max_allowed[4],//pixel range
		min_allowed_prec[4], max_allowed_prec[4];//pixel range + prec bits
	int cdfsize;
	int *pred_errors, *errors, *pixels;
	long long *hist, histsums[NHIST];
	unsigned *CDF_bypass, *CDFs;
	long long esum[4];//sum of abs errors so far in current row
	int
	//	shift[4],//shift right to bring value to 8-bit
		shift_prec[4],//shift right to bring value to 8-bit, removes predictor-added bits
		shift_error[4];
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
	ArithmeticCoder *ec;
} SLIC5Ctx;
#define LOAD(BUF, X, Y) BUF[(pr->kym##Y+kx-(X)+PAD_SIZE)<<2|kc]
static int slic5_init(SLIC5Ctx *pr, int iw, int ih, int nch, const char *depths, ArithmeticCoder *ec)
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
		//pr->nlevels[kc]=nlevels;
		pr->half[kc]=nlevels>>1;
		//pr->min_allowed[kc]=-(pr->nlevels[kc]>>1);
		//pr->max_allowed[kc]=(pr->nlevels[kc]>>1)-1;
		pr->min_allowed_prec[kc]=-prec_half;
		pr->max_allowed_prec[kc]=prec_half-1;
		if(maxdepth<pr->depths[kc])
			maxdepth=pr->depths[kc];
	}
	HybridUint hu;
	int extremesym=-((1<<maxdepth)>>1);
	extremesym=extremesym<<1^-(extremesym<0);//pack sign
	hybriduint_encode(extremesym, &hu);//encode -half
	pr->cdfsize=hu.token+1;

	pr->pred_errors=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[NPREDS*4*4]));//NPREDS * 4 rows * 4 comps
	pr->errors=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	pr->pixels=(int*)malloc((iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	pr->hist=(long long*)malloc(pr->cdfsize*sizeof(long long[NHIST]));//WH: cdfsize * NHIST
	pr->CDF_bypass=(unsigned*)malloc((pr->cdfsize+1LL)*sizeof(int));
	pr->CDFs=(unsigned*)malloc((pr->cdfsize+1LL)*sizeof(int[NHIST]));
	pr->sse=(long long*)malloc(sizeof(long long[SSE_STAGES*SSE_SIZE]));
	pr->sse_fr=(long long*)malloc(sizeof(long long[SSE_FR_SIZE]));
	if(!pr->pred_errors||!pr->errors||!pr->pixels||!pr->hist||!pr->CDFs||!pr->CDF_bypass||!pr->sse||!pr->sse_fr)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(pr->pred_errors, 0, (iw+PAD_SIZE*2LL)*sizeof(int[NPREDS*4*4]));
	memset(pr->errors, 0, (iw+PAD_SIZE*2LL)*sizeof(int[4*4]));
	memset(pr->pixels, 0, (iw+PAD_SIZE*2LL)*sizeof(int[4*4]));

	memset(pr->hist, 0, pr->cdfsize*sizeof(long long[NHIST]));
	memset(pr->histsums, 0, sizeof(pr->histsums));
	//long long fillval=1;
	//memfill(pr->hist, &fillval, pr->cdfsize*sizeof(long long[NHIST]), sizeof(fillval));
	//fillval=pr->cdfsize;
	//memfill(pr->histsums, &fillval, sizeof(pr->histsums), sizeof(fillval));

	//memset(pr->CDFs, 0, (pr->cdfsize+1LL)*sizeof(int[NHIST]));
	int sum=0;
	for(int ks=0;ks<pr->cdfsize;++ks)//TODO tune initial CDFs
		sum+=0x10000/(ks+1);
	//	sum+=0x10000/((ks+1)*(ks+1));//X  bad
	for(int ks=0, c=0;ks<pr->cdfsize+1;++ks)
	{
		int freq=0x10000/(ks+1);
	//	int freq=0x10000/((ks+1)*(ks+1));
		pr->CDF_bypass[ks]=(int)(((long long)c<<16)/sum);
		c+=freq;
	}
	//for(int ks=0;ks<pr->cdfsize+1;++ks)//actual bypass
	//	pr->CDF_bypass[ks]=(ks<<16)/pr->cdfsize;
	memfill(pr->CDFs, pr->CDF_bypass, (pr->cdfsize+1)*sizeof(int[NHIST]), (pr->cdfsize+1)*sizeof(int));

	memset(pr->sse, 0, sizeof(long long[SSE_STAGES*SSE_SIZE]));
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
	for(int kt=0;kt<NHIST;++kt)//update CDFs
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
	if(ky)
	{
		//slic5_update_CDFs(pr);
		//for(int kc=0;kc<pr->nch;++kc)
		//{
		//	int x=(int)(pr->esum[kc]/pr->iw);
		//	x=floor_log2(x)+1-PRED_PREC;
		//	pr->shift[kc]=MAXVAR(8, x)-8+PRED_PREC;
		//	pr->esum[kc]=0;
		//}
	}
	else
	{
		for(int kc=0;kc<pr->nch;++kc)//this code is here in case it's dependent on prev row-average error
		{
			pr->shift_prec[kc]=MAXVAR(8, pr->depths[kc])-8+PRED_PREC;
			//pr->shift[kc]=MAXVAR(8, pr->depths[kc])-8;//-((pr->depths[kc]>9)<<1);
			//pr->shift_prec[kc]=pr->shift[kc]+PRED_PREC;
			pr->shift_error[kc]=1+((pr->depths[kc]<=9)<<1);
		}
	}
	pr->ky=ky;
	pr->kym0=(pr->iw+PAD_SIZE*2)*(ky&3);
	pr->kym1=(pr->iw+PAD_SIZE*2)*((ky-1)&3);
	pr->kym2=(pr->iw+PAD_SIZE*2)*((ky-2)&3);
}
static void slic5_predict(SLIC5Ctx *pr, int kc, int kx, int ky)
{
	int idx=(pr->iw*ky+kx)<<2|kc;
	//if(!(idx&idx>>1))
	//if(idx<CDF_UPDATE_PERIOD&&!(idx&idx>>1)||!(idx&(CDF_UPDATE_PERIOD-1)))
	if(!(idx&(CDF_UPDATE_PERIOD-1))||idx<CDF_UPDATE_PERIOD&&(idx&15))
		slic5_update_CDFs(pr);
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
		eN =LOAD(pr->errors,  0, 1),//error = (curr<<8) - pred
		eW =LOAD(pr->errors,  1, 0),
		eNW=LOAD(pr->errors,  1, 1),
		eNE=LOAD(pr->errors, -1, 0),
		eNN=LOAD(pr->errors,  0, 2),
		eWW=LOAD(pr->errors, -2, 0);

	int clamp_lo=(int)N, clamp_hi=(int)N;
	clamp_lo=(int)MINVAR(clamp_lo, W);
	clamp_hi=(int)MAXVAR(clamp_hi, W);
	clamp_lo=(int)MINVAR(clamp_lo, NE);
	clamp_hi=(int)MAXVAR(clamp_hi, NE);

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
	++j, pr->preds[j]=N>>(abs((int)eN)>=abs(N));//v
	++j, pr->preds[j]=W>>(abs((int)eW)>=abs(W));//v
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

	++j;
	int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
	int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
	int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
	int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
	int sum2=(dy+dx)>>pr->shift_prec[kc], diff=(dy-dx)>>pr->shift_prec[kc], diff2=(NE-NW)/8;
	if(sum2>(32<<(PRED_PREC+1)))
		pr->preds[j]=(dx*N+dy*W)/sum2+diff2;
	else if(diff>(12<<(PRED_PREC+1)))
		pr->preds[j]=(N+2*W)/3+diff2;
	else if(diff<-(12<<(PRED_PREC+1)))
		pr->preds[j]=(2*N+W)/3+diff2;
	else
		pr->preds[j]=(N+W)/2+diff2;
	diff=d45-d135;
	if(diff>(32<<(PRED_PREC+1)))
		pr->preds[j]+=diff2;
	else if(diff>(16<<(PRED_PREC+1)))
		pr->preds[j]+=diff2/2;
	else if(diff<-(32<<(PRED_PREC+1)))
		pr->preds[j]-=diff2;
	else if(diff<-(16<<(PRED_PREC+1)))
		pr->preds[j]-=diff2/2;

	++j;
	int disc=(dy-dx)>>pr->shift_prec[kc];
	if(disc>80)
		pr->preds[j]=W;
	else if(disc<-80)
		pr->preds[j]=N;
	else
	{
		pr->preds[j]=(N+W)/2+(NE-NW)/4;
		if(disc>32)
			pr->preds[j]=(pr->preds[j]+W)/2;
		else if(disc>8)
			pr->preds[j]=(3*pr->preds[j]+W)/4;
		else if(disc<-32)
			pr->preds[j]=(pr->preds[j]+N)/2;
		else if(disc<-8)
			pr->preds[j]=(3*pr->preds[j]+N)/4;
	}
#endif

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

	int
		qx  =QUANTIZE_HIST(dx  >>(pr->shift_prec[kc]-2)),//dx>>shift_pred is in [0 ~ 255*3],  dx>>(shift_pred-2) has quad range
		qy  =QUANTIZE_HIST(dy  >>(pr->shift_prec[kc]-2)),
		q45 =QUANTIZE_HIST(d45 >>(pr->shift_prec[kc]-2)),
		q135=QUANTIZE_HIST(d135>>(pr->shift_prec[kc]-2));
	//int
	//	gx=(int)(W-WW+NE-NW+NN-NNE)>>(pr->shift[kc]+2), qx=quantize_signed(gx),//X
	//	gy=(int)(W-NW+N -NN+NE-NNE)>>(pr->shift[kc]+2), qy=quantize_signed(gy),
	//	g45=(int)(((N-W)<<1)+NN-WW)>>(pr->shift[kc]+2), q45=quantize_signed(g45),
	//	g135=(int)(N-NNW+NE-NN)>>(pr->shift[kc]+2), q135=quantize_signed(g135);
	pr->hist_idx=MAXVAR(qx, qy);
	pr->hist_idx=MAXVAR(pr->hist_idx, q45);
	pr->hist_idx=MAXVAR(pr->hist_idx, q135);
	pr->hist_idx=MINVAR(pr->hist_idx, NHIST-1);
	
	//qx=quantize_signed((int)(W-WW+NE-NW+NN-NNE)>>(pr->shift[kc]+2));
	//qy=quantize_signed((int)(W-NW+N -NN+NE-NNE)>>(pr->shift[kc]+2));
	//q45=(int)(((N-W)<<1)+NN-WW)>>(pr->shift[kc]+2), q45=quantize_signed(q45);
	//q135=(int)(N-NNW+NE-NN)>>(pr->shift[kc]+2), q135=quantize_signed(q135);

	int shifts[]=//3		4 is best for 8-bit
	{
		MAXVAR(0, pr->shift_prec[kc]-2),//amplification
		MAXVAR(0, pr->shift_prec[kc]-1),
		pr->shift_prec[kc]+5,//attenuation
		pr->shift_prec[kc]+6,
	};
	int g[]=
	{
		//4D SSE
#if 1
		SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)(W-WW)>>shifts[1])+QUANTIZE_SSE((int)(N-NW)>>shifts[3]))+QUANTIZE_SSE((int)(NNE-NN)>>shifts[1]),
		SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)(W-NW)>>shifts[1])+QUANTIZE_SSE((int)(N-NN)>>shifts[3]))+QUANTIZE_SSE((int)(NE-NNE)>>shifts[1]),
		SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)eN>>shifts[0])+QUANTIZE_SSE((int)eW>>shifts[2]))+QUANTIZE_SSE((int)eNW>>shifts[0]),
		SSE_WIDTH*(SSE_WIDTH*QUANTIZE_SSE((int)eNN>>shifts[0])+QUANTIZE_SSE((int)eWW>>shifts[2]))+QUANTIZE_SSE((int)eNE>>shifts[0]),

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
	pr->sse_corr=0;
	for(int k=0;k<SSE_STAGES+1;++k)
	{
		long long sse_val;
		if(k<SSE_STAGES)
		{
			int qp=QUANTIZE_PRED((int)(pr->pred>>(pr->shift_prec[kc]+4)));
			pr->sse_idx[k]=SSE_WIDTH*SSE_WIDTH*SSE_WIDTH*qp+g[k];
			sse_val=pr->sse[SSE_SIZE*k+pr->sse_idx[k]];
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
		int sse_corr=pr->sse_count[k]?(int)(pr->sse_sum[k]/pr->sse_count[k]):0;
		long long pred0=pr->pred;
		pr->pred+=sse_corr;
		//pr->pred=CLAMP(clamp_lo, pr->pred, clamp_hi);//X  bad
		pr->sse_corr+=sse_corr;
		//pr->sse_corr+=(int)(pr->pred-pred0);//X  bad

		pr->pred_sse[k]=pr->pred;
	}
	int final_corr=pr->bias_count[kc]?(int)(pr->bias_sum[kc]/pr->bias_count[kc]):0;
	pr->pred+=final_corr;
	pr->sse_corr+=final_corr;

	//if(((eN^eW)|(eN^eNW))<0)//clamp only if signs are different		TODO analyze different clamp strategies
	//{
		//clamp_lo=(int)MINVAR(clamp_lo, NEE);//X  bad
		//clamp_hi=(int)MAXVAR(clamp_hi, NEE);
		//clamp_lo=(int)MINVAR(clamp_lo, NW);//X  bad
		//clamp_hi=(int)MAXVAR(clamp_hi, NW);
		pr->pred=CLAMP(clamp_lo, pr->pred, clamp_hi);
	//}
	//pr->pred=MEDIAN3(N, W, pr->pred);//X  bad
	//pr->pred=CLAMP(pr->min_allowed_prec[kc], pr->pred, pr->max_allowed_prec[kc]);//X  bad
	pr->pred_final=(int)((pr->pred+((1<<PRED_PREC)>>1)-1)>>PRED_PREC);
	pr->kc=kc;
	pr->kx=kx;
}
static void slic5_update(SLIC5Ctx *pr, int curr, int token)
{
	int kc=pr->kc, kx=pr->kx;

	int error=(curr<<PRED_PREC)-(int)pr->pred;
	pr->esum[kc]+=abs(error);
	LOAD(pr->errors, 0, 0)=error;
	int errors[NPREDS]={0}, kbest=0;//, kworst=0;
	for(int k=0;k<NPREDS;++k)
	{
		errors[k]=abs((curr<<PRED_PREC)-pr->preds[k]);
		//if(e<0)
		//	LOG_ERROR("");
		pr->pred_errors[(NPREDS*(pr->kym0+pr->kx+PAD_SIZE)+k)<<2|kc]=errors[k];
		if(pr->ky&&pr->kx+1<pr->iw)
			pr->pred_errors[(NPREDS*(pr->kym1+pr->kx+1+PAD_SIZE)+k)<<2|kc]+=errors[k];//eNE += ecurr
		if(errors[kbest]>errors[k])
			kbest=k;
		//if(errors[kworst]<errors[k])
		//	kworst=k;
	}

	//update WP weights
	++pr->params[kbest];
	//pr->params[kbest]+=pr->params[kbest]>>4;//pr->params[kbest]=(int)(pr->params[kbest]*0x14000LL>>16);//X  bad
	//pr->params[kbest]<<=1;//X  bad
	//if(pr->params[kbest]>0x1000000)
	//if(pr->params[kbest]>0xC00)//best for Kodak
	if(pr->params[kbest]>(12<<pr->depths[kc]))
	{
		for(int k=0;k<NPREDS;++k)
			pr->params[k]>>=1;
	}
	//pr->params[kworst]-=pr->params[kworst]>0x10000;//X  bad

	//if(pr->hist_idx==23&&token==32)//
	//	printf("");

	LOAD(pr->pixels, 0, 0)=curr;

	//update hist
	++pr->hist[pr->cdfsize*pr->hist_idx+token];
	++pr->histsums[pr->hist_idx];
	//if(pr->histsums[pr->hist_idx]>((long long)pr->cdfsize*96))
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
	
	//update SSE
	for(int k=0;k<SSE_STAGES+1;++k)
	{
		++pr->sse_count[k];
		pr->sse_sum[k]+=error;
		//pr->sse_sum[k]+=((long long)curr<<PRED_PREC)-pr->pred_sse[k];//X
		if(pr->sse_count[k]>=640)
		{
			pr->sse_count[k]>>=1;
			pr->sse_sum[k]>>=1;
			//pr->sse_sum[k]=(pr->sse_sum[k]+1)>>1;//worse
		}
		long long sse_val=pr->sse_sum[k]<<12|pr->sse_count[k];
		if(k<SSE_STAGES)
			pr->sse[SSE_SIZE*k+pr->sse_idx[k]]=sse_val;
		else
			pr->sse_fr[pr->sse_idx[k]]=sse_val;
	}
	++pr->bias_count[kc];
	pr->bias_sum[kc]+=error;
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
	
#ifndef SLIC5_DISABLE_PALETTE
	unsigned short *palettes[4]={0};
	unsigned short pal_sizes[4]={0};
	//memset(header.palettesizes, 0, sizeof(header.palettesizes));
	int *pal_hist=(int*)malloc(maxlevels*sizeof(int));
	if(!pal_hist)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int kc=0;kc<nch;++kc)
	{
		int depth=im2->depth[kc], nlevels=1<<depth;
		memset(hist, 0, nlevels*sizeof(int));
		for(ptrdiff_t k=0;k<res;++k)
		{
			int val=im2->data[k<<2|kc]+(nlevels>>1);
			val=CLAMP(0, val, nlevels-1);
			++hist[val];
		}
		int palettesize=0;
		for(int k=0;k<nlevels;++k)
			palettesize+=hist[k]!=0;
		if(palettesize>((nlevels+15)>>4))
			pal_sizes[kc]=0;
		else
		{
			pal_sizes[kc]=palettesize;
			im2->depth[kc]=floor_log2_32(palettesize-1)+1;
			palettes[kc]=(unsigned short*)malloc(palettesize*sizeof(short));
			if(!palettes[kc])
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			for(int k=0, idx=0;k<nlevels;++k)
			{
				if(hist[k])
				{
					palettes[kc][idx]=k;
					++idx;
				}
			}
			int half=(palettesize+1)>>1;//ceil_half
			for(int k=0;k<res;++k)
			{
				unsigned short val=im2->data[k<<2|kc];
				int idx=0;
				int L=0, R=palettesize-1;
				while(L<=R)
				{
					idx=(L+R)>>1;
					if(palettes[kc][idx]<val)
						L=idx+1;
					else if(palettes[kc][idx]>val)
						R=idx-1;
					else
						break;
				}
				im2->data[k<<2|kc]=idx-half;
			}
		}
	}
	free(pal_hist);
#endif
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		slic5_nextrow(&pr, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			//if(kx==7&&ky==3)//
			//if(kx==674&&ky==2)//
			//if(kx==2759&&ky==1)//
			//if(kx==1191&&ky==1659)//
			//if(kx==2475&&ky==0)//
			//	printf("");

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

		printf("Hist WH %d*%d:\n", pr.cdfsize, (int)NHIST);
		for(int kt=0;kt<NHIST;++kt)
		{
			long long *curr_hist=pr.hist+pr.cdfsize*kt;
			for(int ks=0;ks<pr.cdfsize;++ks)
				printf("%lld%c", curr_hist[ks], ks<pr.cdfsize-1?' ':'\n');
		}

#if 0
		int kr=3;
		printf("SSE round %d\n", kr+1);
		long long *curr_sse=pr.sse+SSE_SIZE*kr;
		for(int kz=0;kz<SSE_HEIGHT;++kz)
		{
			for(int ky=0;ky<SSE_WIDTH;++ky)
			{
				for(int kx=0;kx<SSE_WIDTH;++kx)
				{
					long long val=curr_sse[SSE_WIDTH*(SSE_WIDTH*kz+ky)+kx];
					long long sum=val>>12;
					int count=(int)(val&0xFFF);
					int corr=count?(int)(sum/count):0;
					printf("%7d%c", corr, kx<SSE_WIDTH-1?' ':'\n');
				}
			}
			printf("\n");
		}
#endif
#if 0
		int xmin=SSE_WIDTH/2, xmax=SSE_WIDTH/2, ymin=SSE_HEIGHT/2, ymax=SSE_HEIGHT/2;
		for(int ks=0;ks<SSE_STAGES;++ks)
		{
			long long *curr_sse=pr.sse+SSE_SIZE*ks;
#ifdef PRINT_SSE
			printf("\nSSE stage %d:\n", ks+1);
#endif
			for(int ky=0;ky<SSE_HEIGHT;++ky)
			{
				for(int kx=0;kx<SSE_WIDTH;++kx)
				{
					long long val=curr_sse[SSE_WIDTH*ky+kx];
					long long sum=val>>12;
					int count=(int)(val&0xFFF);
					int corr=count?(int)(sum/count):0;
#ifdef PRINT_SSE
					printf("%7d%c", corr, kx<SSE_WIDTH-1?' ':'\n');
#endif
					if(count)
					{
						if(xmin>kx)xmin=kx;
						if(xmax<kx)xmax=kx;
						if(ymin>ky)ymin=ky;
						if(ymax<ky)ymax=ky;
					}
				}
			}
		}
		int xmaxreach=abs(xmin-SSE_WIDTH/2), ymaxreach=abs(ymin-SSE_HEIGHT/2), reach;
		reach=abs(xmax-SSE_WIDTH/2), xmaxreach=MAXVAR(xmaxreach, reach);
		reach=abs(ymax-SSE_HEIGHT/2), ymaxreach=MAXVAR(ymaxreach, reach);
		printf("SSE coverage XY %d~%d %d~%d  maxreach %d*%d\n", xmin, xmax, ymin, ymax, xmaxreach, ymaxreach);
#endif
		printf("Bias\n");
		for(int kc=0;kc<src->nch;++kc)
			printf("  %lld/%d = %d\n", pr.bias_sum[kc], pr.bias_count[kc], pr.bias_count[kc]?(int)(pr.bias_sum[kc]/pr.bias_count[kc]):0);
	}
	slic5_free(&pr);
	dlist_clear(&list);
	return 1;
}
int t47_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	double t_start=time_sec();
#ifndef SLIC5_DISABLE_PALETTE
	unsigned short pal_sizes[4];
	//size_t emithistsize=cdfsize*sizeof(short[_countof(qlevels)+1]);
	if(srclen<sizeof(pal_sizes))
	{
		LOG_ERROR("Corrupt file\n");
		//printf("File smaller than header size\n");
		return 0;
	}
	memcpy(pal_sizes, data, sizeof(pal_sizes));
	data+=sizeof(pal_sizes);
	srclen-=sizeof(pal_sizes);
	unsigned short *palettes[4]={0};
	for(int kc=0;kc<dst->nch;++kc)
	{
		if(pal_sizes[kc]>0)
		{
			int bytesize=pal_sizes[kc]*sizeof(short);
			palettes[kc]=(unsigned short*)malloc(bytesize);
			if(!palettes[kc])
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			if(bytesize>srclen)
			{
				LOG_ERROR("Invalid file\n");
				return 0;
			}
			memcpy(palettes[kc], data, bytesize);
			data+=bytesize;
			srclen-=bytesize;
		}
	}
#endif
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
			//if(kx==7&&ky==3)//
			//if(kx==674&&ky==2)//
			//if(kx==1173&&ky==11)//
			//if(kx==0&&ky==0)//
			//if(kx==2759&&ky==1)//
			//if(kx==1191&&ky==1659)//
			//if(kx==2475&&ky==0)//
			//	printf("");

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
#ifndef SLIC5_DISABLE_PALETTE
	for(int kc=0;kc<dst->nch;++kc)
	{
		if(palettes[kc])
		{
			for(int k=0;k<res;++k)
			{
				int val=dst->data[k<<2|kc];
				val+=(pal_sizes[kc]+1)>>1;
				if((unsigned)val>=pal_sizes[kc])
				{
					LOG_ERROR("Palette error");
					return 0;
				}
				dst->data[k<<2|kc]=palettes[kc][val];
			}
			free(palettes[kc]);
		}
	}
#endif
	if(loud)
	{
		printf("\n");
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	return 1;
}