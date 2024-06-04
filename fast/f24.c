#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#ifdef _MSC_VER
#include<intrin.h>
#endif
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT
//	#define PROFILE_BYPASS

	#define USE_PALETTE
	#define USE_AC		//good
//	#define USE_GRABAC	//bad
//	#define CTX_MIN_NW	//bad
//	#define ENABLE_ZERO	//bad	incompatible with quantized entropy coder
	#define ENABLE_WP	//good
	#define ENABLE_WG	//good
//	#define ENABLE_WG2	//2x slower, not worth the penalty
	#define ENABLE_WG3	//good
	#define FOLD_OOB	//good
//	#define USE_T47CTX	//bad
	#define USE_FP64	//4% faster
	#define WG_DECAY_NUM	493	//493	//3	//230
	#define WG_DECAY_SH	9	//9	//2	//8
	#define WG3_DECAY_NUM	493
	#define WG3_DECAY_SH	9

#include"ac.h"
#define BLOCKSIZE 256
#define NCHPOOL 15
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
#if defined ENABLE_GUIDE && defined USE_PALETTE
static Image g_palimage={0};
#endif
#ifdef USE_GRABAC
#define LR_SHIFT_TOKEN 1
#define LR_SHIFT_BYPASS 16
#define UPDATE_P0(P0, BIT, SH)\
	do\
	{\
		int update=(BIT<<16)-P0;\
		P0+=(update+(1<<SH>>1))>>SH;\
	}while(0)
#endif
typedef enum RCTTypeEnum
{
	RCT_NONE,
	RCT_SUBG,
	RCT_JPEG2000_G,
	RCT_JPEG2000_B,
	RCT_JPEG2000_R,
	RCT_1,
	RCT_PEI09,

	RCT_COUNT,
} RCTType;
static const char *rctnames[RCT_COUNT]=
{
	"NONE      ",
	"SubG      ",
	"JPEG2000_G",
	"JPEG2000_B",
	"JPEG2000_R",
	"RCT1      ",
	"Pei09     ",
};
typedef enum PredTypeEnum
{
#ifdef ENABLE_ZERO
	PRED_ZERO,
#endif
//	PRED_N,
//	PRED_W,
//	PRED_AV2,
//	PRED_AV4,
//	PRED_AV5,
	PRED_CG,
//	PRED_CGU,
#ifdef ENABLE_WP
	PRED_WP,
#endif
#ifdef ENABLE_WG
	PRED_WG,
#endif
#ifdef ENABLE_WG2
	PRED_WG2,
#endif
#ifdef ENABLE_WG3
	PRED_WG3,
#endif

	PRED_COUNT,
} PredType;
static const char *prednames[PRED_COUNT]=
{
#ifdef ENABLE_ZERO
	"Zero",
#endif
//	"N   ",
//	"W   ",
//	"Av2 ",
//	"Av4 ",
//	"Av5 ",
	"CG  ",
//	"CGU ",
#ifdef ENABLE_WP
	"WP  ",
#endif
#ifdef ENABLE_WG
	"WG  ",
#endif
#ifdef ENABLE_WG2
	"WG2 ",
#endif
#ifdef ENABLE_WG3
	"WG3 ",
#endif
};
long long rct_gains[RCT_COUNT]={0}, pred_gains[PRED_COUNT]={0};
typedef struct _ThreadArgs
{
	const Image *src;
	Image *dst;
	int fwd, loud, y1, y2;
	int bufsize, histsize;
	int *pixels, *hist;

	DList list;
	const unsigned char *decstart, *decend;

	double bestsize;
	int bestrct, predidx[3];

#ifdef USE_AC
	int tlevels, clevels, statssize;
	unsigned *stats;
#elif defined USE_GRABAC
	unsigned short p0[16*4], p0_bypass[4];
#endif
#ifdef PROFILE_BYPASS
	ptrdiff_t bypasssizes[4];
#endif
} ThreadArgs;
static const int combinations[RCT_COUNT*3]=
{
	//YUV		X/Y orders are irrelevant here
	 1,  2,  0,//RCT_NONE
	 1,  3,  4,//RCT_SUBG
	 5,  6,  7,//RCT_J2KG	default
	 8,  9,  6,//RCT_J2KB
	10,  7,  9,//RCT_J2KR
	11, 12,  7,//RCT_1
	13, 14,  7,//RCT_PEI09
};
static const int depth_inf[NCHPOOL]=
{
	0, 0, 0,	//RCT_NONE
	0, 0,		//RCT_SubG (causal RCT)
	0, 1, 1,	//RCT_J2KG
	0, 1,		//RCT_J2KB
	0,		//RCT_J2KR
	0, 1,		//RCT_1
	0, 1,		//RCT_Pei09
};
#ifdef ENABLE_WP
static const int wp_param_idx[NCHPOOL]=
{
	0, 1, 2,	//RCT_NONE
	2, 1,		//RCT_SubG
	0, 1, 2,	//RCT_J2KG
	0, 2,		//RCT_J2KB
	0,		//RCT_J2KR
	0, 1,		//RCT_1
	0, 1,		//RCT_Pei09
};
#define WP_NPARAMS 11
static const short wp_params[WP_NPARAMS*3]=//signed fixed 7.8 bit
{
	0x0DB8,  0x0E22,  0x181F,  0x0BF3, -0x005C, -0x005B,  0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,//Y
	0x064C,  0x0F31,  0x1040,  0x0BF8, -0x0007, -0x000D, -0x0085, -0x0063, -0x00A2, -0x0017,  0x00F2,//Cb
	0x0B37,  0x110B,  0x121B,  0x0BFC, -0x0001,  0x000E, -0x0188, -0x00E7, -0x00BB, -0x004A,  0x00BA,//Cr
};
#endif
#ifdef USE_AC

#define CDFSTRIDE 64

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
		*nbits=lgv-CONFIG_MSB+CONFIG_LSB;
		*bypass=val>>CONFIG_LSB&((1LL<<*nbits)-1);
	}
}
#ifdef USE_T47CTX
#define quantize_ctx FLOOR_LOG2_P1
#define QUANTIZE_CTX FLOOR_LOG2_P1
#else
static int quantize_ctx(int val)
{
#if 1
	int negmask=val>>31;
	val=abs(val);
	val=FLOOR_LOG2_P1(val);//[0~0x8000] -> [0~16]
	val=FLOOR_LOG2_P1(val);//[0~16] -> [0~5]
	val^=negmask;
	val-=negmask;
#ifdef CHECK_OOB
	if((unsigned)val>=CLEVELS)
		LOG_ERROR("Context OOB");
#endif
	return val;
#else
	return THREEWAY(val, 0)+1;
#endif
}
#define QUANTIZE_CTX(X) quantize_ctx(X)+chalf
#endif
static void update_CDF(const int *hist, unsigned *CDF, int tlevels)
{
	int sum=hist[tlevels], c=0;
	for(int ks=0;ks<tlevels;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*((1LL<<PROB_BITS)-tlevels)/sum)+ks;
		c+=freq;
	}
	CDF[tlevels]=1<<PROB_BITS;
}
static void init_stats(int *hist, unsigned *stats, int tlevels, int clevels, int nch)
{
	int stride=tlevels+1;
#ifdef CTX_MIN_NW
	int chsize=clevels*stride;
#else
	int chsize=clevels*clevels*stride;
#endif
	for(int kc=0;kc<nch;++kc)
	{
		int sum=0;
		int *curr_hist=hist+chsize*kc;
		unsigned *curr_CDF=stats+chsize*kc;

#if 1
		//memset(curr_hist, 0, chsize);
		//memset(curr_CDF, 0, chsize);
		for(int k=0;k<tlevels;++k)
		{
			sum+=curr_hist[k]=1;
			//int val=(tlevels-k)>>1;
			//val+=!val;
			//val*=val;
			//val*=val*val;
			//if(sum+val*val*val<sum)
			//	LOG_ERROR("");
			//sum+=curr_hist[k]=val;
			//sum+=curr_hist[k]=tlevels-(k>>2);
		}
#else
		int val=0x400;
		for(int k=0;k<tlevels;++k)//X
		{
			sum+=curr_hist[k]=val;
			val-=val>>5;
		}
#endif
		curr_hist[tlevels]=sum;
		update_CDF(curr_hist, curr_CDF, tlevels);
		memfill(curr_hist+stride, curr_hist, ((size_t)chsize-stride)*sizeof(int), stride*sizeof(int));
		memfill(curr_CDF+stride, curr_CDF, ((size_t)chsize-stride)*sizeof(int), stride*sizeof(int));
	}
}
#endif
#ifdef ENABLE_WP
static int wp_predict(const short *params, int NN, int NW, int N, int NE, int W, const int *eNW, const int *eN, const int *eNE, const int *eW, int *preds)
{
	int pred;
	long long lpred=0, wsum=0;
	preds[0]=(W+NE-N)<<8;
	preds[1]=(N<<8)-((eN[0]+eW[0]+eNE[0])*params[4]/0x100);
	preds[2]=(W<<8)-((eN[0]+eW[0]+eNW[0])*params[5]/0x100);
	preds[3]=(N<<8)-(((eNW[0]*params[6]+eN[0]*params[7]+eNE[0]*params[8])/0x100)+(NN-N)*params[9]+(NW-W)*params[10]);
	for(int k=0;k<4;++k)
	{
		int weight=(params[k]<<8)/(eNW[k+1]+eN[k+1]+eNE[k+1]+1);
		lpred+=(long long)weight*preds[k];
		wsum+=weight;
	}
	lpred+=wsum/2-1;
	lpred/=wsum+1;
	if(!wsum)
		lpred=preds[0];
	CLAMP3_32(pred, (int)lpred, N<<8, W<<8, NE<<8);
	return pred;
}
static void wp_update(int curr, int pred, const int *preds, int *ecurr, int *eNE)
{
	curr<<=8;
	ecurr[0]=curr-(int)pred;
	for(int k=0;k<4;++k)
	{
		int e2=abs(curr-preds[k]);
		ecurr[k+1]=e2;
		eNE[k+1]+=e2;
	}
}
#endif
#ifdef ENABLE_WG
static int wg_predict(int NNN, int NN, int NNE, int NW, int N, int NE, int NEEE, int WWW, int WW, int W, const int *eW, int *preds)
{
	int pred;
#ifdef USE_FP64
	double pred2=0, wsum=0;
#else
	long long pred2=0, wsum=0;
#endif

	preds[0]=N+W-NW;
	preds[1]=W+NE-N;
	preds[2]=N+NE-NNE;
	preds[3]=N;
	preds[4]=W;
	preds[5]=3*(N-NN)+NNN;
	preds[6]=3*(W-WW)+WWW;
	preds[7]=(W+NEEE)/2;
	
	for(int k=0;k<8;++k)
	{
#ifdef USE_FP64
		double weight=1./(eW[k]+1);
		pred2+=weight*preds[k];
		wsum+=weight;
#else
		int weight=0x1000000/(eW[k]+1);
		pred2+=(long long)weight*preds[k];
		wsum+=weight;
#endif
	}
	pred2/=wsum;
#ifdef USE_FP64
#ifdef _MSC_VER
	pred=_cvt_dtoi_fast(pred2);
#else
	pred=(int)pred2;
#endif
#else
	pred=(int)pred2;
#endif
	//lpred=0;
	//wsum=0;
	//for(int k=0;k<8;++k)
	//{
	//	int weight=0x1000000/(eW[k]+1);
	//	lpred+=(long long)weight*preds[k];
	//	wsum+=weight;
	//}
	//lpred/=wsum+1;
	//pred=(int)lpred;
	CLAMP3_32(pred, pred, N, W, NE);
	return pred;
}
static void wg_update(int curr, const int *preds, int *eW)
{
	for(int k=0;k<8;++k)
		eW[k]=(eW[k]+abs(curr-preds[k]))*WG_DECAY_NUM>>WG_DECAY_SH;
}
#endif
#ifdef ENABLE_WG2
#define STRIDE 7
#define WG2_NPREDS 32
#define WG2_PREDLIST\
	WG2_PRED(W+NE-N-((2*(eN+eW)+eNE-eNW+4)>>3))\
	WG2_PRED(N-(int)(((long long)eN+eW+eNE)*-0x05C>>8))\
	WG2_PRED(W-(int)(((long long)eN+eW+eNW)*-0x05B>>8))\
	WG2_PRED(N+(int)((-eNN*0x0DFLL-eN*0x051LL-eNE*0x0BDLL+((long long)N-NN)*0x05C+((long long)NW-W)*0x102)>>8))\
	WG2_PRED(3*(N-NN)+NNN)\
	WG2_PRED((N+W)>>1)\
	WG2_PRED(N+W-NW)\
	WG2_PRED((W+NEE)>>1)\
	WG2_PRED((3*W+NEEE)>>2)\
	WG2_PRED((3*(3*W+NE+NEE)-10*N+2)/5)\
	WG2_PRED((3*(3*W+NE+NEE)-10*N)/5)\
	WG2_PRED((4*N-2*NN+NW+NE)>>2)\
	WG2_PRED(N+NE-NNE-eNNE)\
	WG2_PRED((4*(N+W+NW+NE)-(NN+WW+NNWW+NNEE)+6)/12)\
	WG2_PRED(W+((eW-eWW)>>1))\
	WG2_PRED(paper_GAP)\
	WG2_PRED(calic_GAP)\
	WG2_PRED(N+W-((NW+NN+WW+NE)>>2))\
	WG2_PRED(((2*(N+W)-(NW+NN+WW+NE))*9+(WWW+NWW+NNW+NNN+NNE+NEE)*2)/12)\
	WG2_PRED(3*(N+W-NW-(NN+WW-NNWW))+NNN+WWW-NNNWWW)\
	WG2_PRED(2*(W+NE-N)-(WW+NNEE-NN))\
	WG2_PRED((2*W+NEE-N)>>1)\
	WG2_PRED(NW+NWW-NNWWW)\
	WG2_PRED((14*NE-(NNEE+NNNEE+NNEEE))/11)\
	WG2_PRED((NEEE+NEEEE)>>1)\
	WG2_PRED((NNNEEEE+NNEEE)>>1)\
	WG2_PRED(NNEEEE)\
	WG2_PRED((NNWWWW+NNNWWWW)>>1)\
	WG2_PRED((WWW+WWWW)>>1)\
	WG2_PRED((N+NN)>>1)\
	WG2_PRED((NE+NNEE)>>1)\
	WG2_PRED((NE+NNE+NEE+NNEE)>>2)
static int wg2_predict(int **rows, int stride, int kc2, int eoffset, int depth, const int *eprev, int *preds)
{
	int sh_up=24-depth, sh_dn=depth-8;
	int
		NNNWWWW	=rows[3][kc2-4*stride]<<sh_up,
		NNNWWW	=rows[3][kc2-3*stride]<<sh_up,
	//	NNNWW	=rows[3][kc2-2*stride]<<sh_up,
	//	NNNW	=rows[3][kc2-1*stride]<<sh_up,
		NNN	=rows[3][kc2+0*stride]<<sh_up,
	//	NNNE	=rows[3][kc2+1*stride]<<sh_up,
		NNNEE	=rows[3][kc2+2*stride]<<sh_up,
	//	NNNEEE	=rows[3][kc2+3*stride]<<sh_up,
		NNNEEEE	=rows[3][kc2+4*stride]<<sh_up,
		NNWWWW	=rows[2][kc2-4*stride]<<sh_up,
		NNWWW	=rows[2][kc2-3*stride]<<sh_up,
		NNWW	=rows[2][kc2-2*stride]<<sh_up,
		NNW	=rows[2][kc2-1*stride]<<sh_up,
		NN	=rows[2][kc2+0*stride]<<sh_up,
		NNE	=rows[2][kc2+1*stride]<<sh_up,
		NNEE	=rows[2][kc2+2*stride]<<sh_up,
		NNEEE	=rows[2][kc2+3*stride]<<sh_up,
		NNEEEE	=rows[2][kc2+4*stride]<<sh_up,
	//	NWWWW	=rows[1][kc2-4*stride]<<sh_up,
	//	NWWW	=rows[1][kc2-3*stride]<<sh_up,
		NWW	=rows[1][kc2-2*stride]<<sh_up,
		NW	=rows[1][kc2-1*stride]<<sh_up,
		N	=rows[1][kc2+0*stride]<<sh_up,
		NE	=rows[1][kc2+1*stride]<<sh_up,
		NEE	=rows[1][kc2+2*stride]<<sh_up,
		NEEE	=rows[1][kc2+3*stride]<<sh_up,
		NEEEE	=rows[1][kc2+4*stride]<<sh_up,
		WWWW	=rows[0][kc2-4*stride]<<sh_up,
		WWW	=rows[0][kc2-3*stride]<<sh_up,
		WW	=rows[0][kc2-2*stride]<<sh_up,
		W	=rows[0][kc2-1*stride]<<sh_up,
	//	eNNWW	=rows[2][kc2-2*stride+eoffset],
	//	eNNW	=rows[2][kc2-1*stride+eoffset],
		eNN	=rows[2][kc2+0*stride+eoffset],
		eNNE	=rows[2][kc2+1*stride+eoffset],
	//	eNNEE	=rows[2][kc2+2*stride+eoffset],
	//	eNWW	=rows[1][kc2-2*stride+eoffset],
		eNW	=rows[1][kc2-1*stride+eoffset],
		eN	=rows[1][kc2+0*stride+eoffset],
		eNE	=rows[1][kc2+1*stride+eoffset],
	//	eNEE	=rows[1][kc2+2*stride+eoffset],
		eWW	=rows[0][kc2-2*stride+eoffset],
		eW	=rows[0][kc2-1*stride+eoffset];
	int pred;
	long long lpred=0, wsum=0;
	int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
	int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
	int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
	int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
	int diff=(dy-dx)>>sh_dn, diff2=(d45-d135)>>sh_dn, diff3=NE-NW;
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
	int j=0;
	
#define WG2_PRED(X) preds[j++]=X;
	WG2_PREDLIST
#undef  WG2_PRED
	
	lpred=0;
	wsum=0;
	for(int k=0;k<WG2_NPREDS;++k)
	{
		int weight=0xFFFFFFFF/(eprev[k]+1);
		lpred+=(long long)weight*preds[k];
		wsum+=weight;
	}
	lpred/=wsum+1;
	pred=(int)lpred;
	CLAMP3_32(pred, (int)lpred, N, W, NE);
	return pred;
}
static void wg2_update(int curr24, const int *preds, int *eprev)
{
	for(int k=0;k<WG2_NPREDS;++k)
	{
		int err=abs(curr24-preds[k]);
		eprev[k]+=err;
		eprev[k]-=(eprev[k]+3)>>2;
	}
}
#else
#define STRIDE 6
#endif
#ifdef ENABLE_WG3
#define WG3_NPREDS 8
static int wg3_predict(int cgrad, int NNN, int NN, int NNE, int NW, int N, int NE, int NEE, int NEEE, int WWW, int WW, int W, const int *eW, int *preds)
{
	int pred;
#ifdef USE_FP64
	double pred2=0, wsum=0;
#else
	long long pred2=0, wsum=0;
#endif
	int j=0;

	preds[j++]=cgrad;
	preds[j++]=(N+W)/2;
	preds[j++]=(W+NEE)/2;
	preds[j++]=W+NE-N;
	preds[j++]=2*N-NN;
	preds[j++]=(4*W+N+NE+NEE+NEEE)/8;
	preds[j++]=(N+NN)/2;
	preds[j++]=(N+NE)/2;

	//preds[j++]=N+W-NW;
	//preds[j++]=W+NE-N;
	//preds[j++]=N+NE-NNE;
	//preds[j++]=N;
	//preds[j++]=W;
	//preds[j++]=3*(N-NN)+NNN;
	//preds[j++]=3*(W-WW)+WWW;
	//preds[j++]=(W+NEEE)/2;
	
	for(int k=0;k<WG3_NPREDS;++k)
	{
#ifdef USE_FP64
		double weight=1./(eW[k]+1);
		pred2+=weight*preds[k];
		wsum+=weight;
#else
		int weight=0x1000000/(eW[k]+1);
		pred2+=(long long)weight*preds[k];
		wsum+=weight;
#endif
	}
	pred2/=wsum;
#ifdef USE_FP64
#ifdef _MSC_VER
	pred=_cvt_dtoi_fast(pred2);
#else
	pred=(int)pred2;
#endif
#else
	pred=(int)pred2;
#endif
	CLAMP3_32(pred, pred, N, W, NE);
	return pred;
}
static void wg3_update(int curr, const int *preds, int *eW)
{
	for(int k=0;k<WG3_NPREDS;++k)
		eW[k]=(eW[k]+abs(curr-preds[k]))*WG3_DECAY_NUM>>WG3_DECAY_SH;
}
#endif
static void block_enc(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
#ifdef USE_AC
	ArithmeticCoder ec;
	int chalf=args->clevels>>1;
#ifdef CTX_MIN_NW
	int nctx=args->clevels;
#else
	int nctx=args->clevels*args->clevels;
#endif
	int cdfstride=args->tlevels+1;
	int token=0, bypass=0, nbits=0;
#elif defined USE_GRABAC
	ArithmeticCoder ec;
#else
	GolombRiceCoder ec;
#endif
	Image const *image=args->src;
	int nlevels[NCHPOOL]={0}, halfs[NCHPOOL]={0};
	double csizes[NCHPOOL*PRED_COUNT]={0};
#ifdef ENABLE_WP
	const short *ch_params[NCHPOOL]={0};
#endif
	int res=image->iw*(args->y2-args->y1);

	for(int k=0;k<NCHPOOL;++k)
	{
		nlevels[k]=1<<(image->depth+depth_inf[k]);
		halfs[k]=nlevels[k]>>1;
#ifdef ENABLE_WP
		ch_params[k]=wp_params+WP_NPARAMS*wp_param_idx[k];
#endif
	}
	memset(args->pixels, 0, args->bufsize);
	memset(args->hist, 0, args->histsize);
	for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
	{
		ALIGN(16) int *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*STRIDE*NCHPOOL,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*STRIDE*NCHPOOL,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*STRIDE*NCHPOOL,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*STRIDE*NCHPOOL,
		};
		short comp[NCHPOOL]={0};
		short
			*comp_none=comp,
			*comp_subg=comp_none+3,
			*comp_j2kg=comp_subg+2,
			*comp_j2kb=comp_j2kg+3,
			*comp_j2kr=comp_j2kb+2,
			*comp_rct1=comp_j2kr+1,
			*comp_pei9=comp_rct1+2;
		int preds[PRED_COUNT]={0};
		int wg_errors[8]={0};
#ifdef ENABLE_WG2
		int wg2_eprev[WG2_NPREDS]={0}, wg2_preds[WG2_NPREDS]={0}, wg2pred=0;
#endif
#ifdef ENABLE_WG3
		int wg3_eprev[WG3_NPREDS]={0}, wg3_preds[WG3_NPREDS]={0};
#endif
		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			int
				*NNN	=rows[3]+0*STRIDE*NCHPOOL,
				*NN	=rows[2]+0*STRIDE*NCHPOOL,
				*NNE	=rows[2]+1*STRIDE*NCHPOOL,
				*NW	=rows[1]-1*STRIDE*NCHPOOL,
				*N	=rows[1]+0*STRIDE*NCHPOOL,
				*NE	=rows[1]+1*STRIDE*NCHPOOL,
				*NEE	=rows[1]+2*STRIDE*NCHPOOL,
				*NEEE	=rows[1]+3*STRIDE*NCHPOOL,
				*WWW	=rows[0]-3*STRIDE*NCHPOOL,
				*WW	=rows[0]-2*STRIDE*NCHPOOL,
				*W	=rows[0]-1*STRIDE*NCHPOOL,
				*curr	=rows[0]+0*STRIDE*NCHPOOL;
			//          NNN
			//          NN  NNE
			//       NW N   NE  .  NEEE
			//WWW WW W  ?
#ifdef FOLD_OOB
			if(ky<=args->y1+2)
			{
				if(ky<=args->y1+1)
				{
					if(ky==args->y1)
						NEEE=NE=NW=N=W;
					NN=N;
					NNE=NE;
				}
				NNN=NN;
			}
			if(kx<=2)
			{
				if(kx<=1)
				{
					if(!kx)
						NW=W=N;
					WW=W;
				}
				WWW=WW;
			}
			if(kx>=image->iw-3)
			{
				if(kx>=image->iw-1)
				{
					NNE=NN;
					NE=N;
				}
				NEEE=NE;
			}
#endif

			//0	r-P(rprev)		NONE		r; g; b;
			//1	g-P(gprev)
			//2	b-P(bprev)
			//
			//3	r-clamp(P(rprev-gprev)+gcurr)		SubG	r-=g; b-=g; g;
			//4	b-clamp(P(bprev-gprev)+gcurr)
			//
			//5	y1-P(y1prev)	y1 = (r+2g+b)/4		JPEG2000-lumaG		r-=g; b-=g; g+=(r+b)>>2;
			//6	u1-P(u1prev)	u1 = b-g
			//7	v1-P(v1prev)	v1 = r-g
			//
			//8	y2-P(y2prev)	y2 = (r+g+2b)/4		JPEG2000-lumaB		r-=b; g-=b; b+=(r+g)>>2;
			//9	u2-P(u2prev)	u2 = r-b
			// X	v2-P(v2prev)	v2 = g-b = -u1
			//
			//10	y3-P(y3prev)	y3 = (2r+g+b)/4		JPEG2000-lumaR		b-=r; g-=r; r+=(b+g)>>2;
			// X	u3-P(u3prev)	u3 = g-r = -v1
			// X	v3-P(v3prev)	v3 = b-r = -u2
			//
			//11	y4-P(y4prev)		RCT1		r-=g; g+=r>>1; b-=g; g+=b>>1;
			//12	u4-P(u4prev)
			// X	v4-P(v4prev)	v4 = v1
			//
			//13	y5-P(y5prev)		Pei09		b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
			//14	u5-P(u5prev)
			// X	v5-P(v5prev)	v5 = v1

			//NONE
			comp_none[0]=image->data[idx+0];//r
			comp_none[1]=image->data[idx+1];//g
			comp_none[2]=image->data[idx+2];//b
			
			//SubG+
			comp_subg[0]=comp_none[0];//cr0
			comp_subg[1]=comp_none[2];//cb0

			//JPEG2000_RCT-lumaG	0 [1] 2
			comp_j2kg[2]=comp_none[0]-comp_none[1];//v1 = r-g		cr1 reused in RCT1 & Pei09
			comp_j2kg[1]=comp_none[2]-comp_none[1];//u1 = b-g
			comp_j2kg[0]=comp_none[1]+((comp_j2kg[2]+comp_j2kg[1])>>2);//y1 = (r+2g+b)/4

			//JPEG2000_RCT-lumaB	0 1 [2]
			{
				int comp_v2=-comp_j2kg[1];//v2 = g-b = -u1
				comp_j2kb[1]=comp_none[0]-comp_none[2];//u2 = r-b
				comp_j2kb[0]=comp_none[2]+((comp_v2+comp_j2kb[1])>>2);//y2 = (r+g+2b)/4
			}

			//JPEG2000_RCT-lumaR	[0] 1 2
			{
				int comp_v3=-comp_j2kb[1];//v3 = b-r = -u2
				int comp_u3=-comp_j2kg[2];//u3 = g-r = -v1
				comp_j2kr[0]=comp_none[1]+((comp_u3+comp_v3)>>2);//y3 = (2r+g+b)/4
			}

			//RCT1
			comp_rct1[0]=comp_none[1]+(comp_j2kg[2]>>1);
			comp_rct1[1]=comp_none[2]-comp_rct1[0];//cb4
			comp_rct1[0]+=comp_rct1[1]>>1;//y4

			//Pei09
			comp_pei9[1]=comp_none[2]-((87*comp_none[0]+169*comp_none[1]+128)>>8);//cb5
			comp_pei9[0]=comp_none[1]+((86*comp_j2kg[2]+29*comp_pei9[1]+128)>>8);//y5
			for(int kc=0;kc<NCHPOOL;++kc)
			{
				int kc2=kc*STRIDE;
#ifdef ENABLE_WP
				int wpred, wp_preds[4];
#endif
#ifdef ENABLE_WG
				int wg_preds[8];
#endif
				int offset=0, val;
				//int
				//	vx=abs(W[kc2+1]-WW[kc2+1])+abs(N[kc2+1]-NW[kc2+1])+abs(NE[kc2+1]-N[kc2+1])+1,
				//	vy=abs(W[kc2+1]-NW[kc2+1])+abs(N[kc2+1]-NN[kc2+1])+abs(NE[kc2+1]-NNE[kc2+1])+1;
				//int update=(vx<vy?N[kc2+1]:W[kc2+1])/4;
#ifdef ENABLE_ZERO
				preds[PRED_ZERO]=0;
#endif
			//	preds[PRED_N]=N[kc2];
			//	preds[PRED_W]=W[kc2];
			//	preds[PRED_AV2]=(N[kc2]+W[kc2])/2;
			//	preds[PRED_AV4]=(4*(N[kc2]+W[kc2])+NE[kc2]-NW[kc2])/8;
			//	CLAMP3_32(preds[PRED_AV5], W[kc2]+(5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])/8, N[kc2], W[kc2], NE[kc2]);
				MEDIAN3_32(preds[PRED_CG], N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
			//	CLAMP2_32(cgrad, cgrad+update, -halfs[kc], halfs[kc]-1);
			//	preds[PRED_CGU]=cgrad;
#ifdef ENABLE_WP
				wpred=wp_predict(ch_params[kc], NN[kc2], NW[kc2], N[kc2], NE[kc2], W[kc2], NW+kc2+1, N+kc2+1, NE+kc2+1, W+kc2+1, wp_preds);
				preds[PRED_WP]=(wpred+127)>>8;
#endif
#ifdef ENABLE_WG
				preds[PRED_WG]=wg_predict(NNN[kc2], NN[kc2], NNE[kc2], NW[kc2], N[kc2], NE[kc2], NEEE[kc2], WWW[kc2], WW[kc2], W[kc2], wg_errors, wg_preds);
#endif
#ifdef ENABLE_WG2
				int depth=image->depth+depth_inf[kc];
				wg2pred=wg2_predict(rows, STRIDE*NCHPOOL, kc2, 6, depth, wg2_eprev, wg2_preds);
				preds[PRED_WG2]=wg2pred>>(24-depth);
#endif
#ifdef ENABLE_WG3
				preds[PRED_WG3]=wg3_predict(preds[PRED_CG], NNN[kc2], NN[kc2], NNE[kc2], NW[kc2], N[kc2], NE[kc2], NEE[kc2], NEEE[kc2], WWW[kc2], WW[kc2], W[kc2], wg3_eprev, wg3_preds);
#endif
				if((unsigned)(kc-3)<2)//{r, g, b, [r-g, b-g], JPEG2000_G, ...}
				{
					offset=comp[1];
					for(int k=0;k<PRED_COUNT;++k)
					{
						preds[k]+=offset;
						CLAMP2_32(preds[k], preds[k], -halfs[kc], halfs[kc]-1);
					}
				}
				for(int k=0;k<PRED_COUNT;++k)
				{
					val=comp[kc]-preds[k];
					val+=halfs[kc];
					val&=nlevels[kc]-1;

					//if(kc*PRED_COUNT+k==3&&val==134)//
					//	printf("");

					++args->hist[(kc*PRED_COUNT+k)<<(image->depth+1)|val];
				}
				curr[kc2+0]=comp[kc]-offset;
#ifdef ENABLE_WP
				wp_update(comp[kc], wpred, wp_preds, curr+kc2+1, NE+kc2+1);
#endif
#ifdef ENABLE_WG
				wg_update(comp[kc], wg_preds, wg_errors);
#endif
#ifdef ENABLE_WG2
				{
					int curr24=comp[kc]<<(24-depth);
					curr[kc2+6]=curr24-wg2pred;
					wg2_update(curr24, wg2_preds, wg2_eprev);
				}
#endif
#ifdef ENABLE_WG3
				wg3_update(comp[kc], wg3_preds, wg3_eprev);
#endif
			}
			rows[0]+=STRIDE*NCHPOOL;
			rows[1]+=STRIDE*NCHPOOL;
			rows[2]+=STRIDE*NCHPOOL;
			rows[3]+=STRIDE*NCHPOOL;
		}
	}
	for(int kc=0;kc<NCHPOOL*PRED_COUNT;++kc)
	{
		int *curr_hist=args->hist+((size_t)kc<<(image->depth+1));
		int nlevels2=nlevels[kc/PRED_COUNT];
		for(int ks=0;ks<nlevels2;++ks)
		{
			int freq=curr_hist[ks];
			if(freq)
				csizes[kc]-=freq*log2((double)freq/res);
		}
		csizes[kc]/=8;
	}
	int bestrct=0;
	double bestsize;
	const int *group;
	int predsel[NCHPOOL]={0};
	for(int k=0;k<NCHPOOL;++k)
	{
		int best=0;
		for(int k2=1;k2<PRED_COUNT;++k2)
		{
			if(csizes[k*PRED_COUNT+best]>csizes[k*PRED_COUNT+k2])
				best=k2;
		}
		predsel[k]=best;
	}
	group=combinations;
	bestrct=0;
	bestsize=
		csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
		csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
		csizes[group[2]*PRED_COUNT+predsel[group[2]]];
#if 1
	for(int k=1;k<RCT_COUNT;++k)
	{
		group=combinations+k*3;
		double csize=
			csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
			csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
			csizes[group[2]*PRED_COUNT+predsel[group[2]]];
		if(bestsize>csize)
			bestsize=csize, bestrct=k;
	}
#endif
	int combination[]=
	{
		combinations[bestrct*3+0],
		combinations[bestrct*3+1],
		combinations[bestrct*3+2],
	};
	int predidx[]=
	{
		predsel[combinations[bestrct*3+0]],
		predsel[combinations[bestrct*3+1]],
		predsel[combinations[bestrct*3+2]],
	};
	int flag=bestrct;
	flag=flag*PRED_COUNT+predidx[2];
	flag=flag*PRED_COUNT+predidx[1];
	flag=flag*PRED_COUNT+predidx[0];

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
			rctnames[bestrct],
			prednames[predidx[0]],
			prednames[predidx[1]],
			prednames[predidx[2]]
		);
		{
			const char *poolnames[]=
			{
				"r-P(rprev)",
				"g-P(gprev)",
				"b-P(bprev)",

				"r-clamp(P(rprev-gprev)+gcurr)",
				"b-clamp(P(bprev-gprev)+gcurr), y=g",

				"J2KG: y1-P(y1prev)",
				"J2KG: u1-P(u1prev)",
				"J2KG: v1-P(v1prev)",

				"J2KB: y2-P(y2prev)",
				"J2KB: u2-P(u2prev), v2=-u1",

				"J2KR: y3-P(y3prev), u3=-v1, v3=-u2",

				"RCT1: y4-P(y4prev)",
				"RCT1: u4-P(u4prev), v4=v1",

				"Pei9: y5-P(y5prev)",
				"Pei9: u5-P(u5prev), v5=v1",
			};
			for(int k2=0;k2<PRED_COUNT;++k2)
				printf("%-14s", prednames[k2]);
			printf("\n");
			for(int k=0;k<NCHPOOL;++k)
			{
				int best=0;
				for(int k2=1;k2<PRED_COUNT;++k2)
				{
					if(csizes[k*PRED_COUNT+best]>csizes[k*PRED_COUNT+k2])
						best=k2;
				}
				for(int k2=0;k2<PRED_COUNT;++k2)
					printf(" %12.2lf%c", csizes[k*PRED_COUNT+k2], k2==best?'*':' ');
				printf("  %s\n", poolnames[k]);
			}
		}
		for(int k=0;k<RCT_COUNT;++k)
		{
			group=combinations+k*3;
			double csize=
				csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
				csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
				csizes[group[2]*PRED_COUNT+predsel[group[2]]];
			printf("%12.2lf %c  %s %s %s %s\n",
				csize,
				k==bestrct?'*':' ',
				prednames[predsel[group[0]]],
				prednames[predsel[group[1]]],
				prednames[predsel[group[2]]],
				rctnames[k]
			);
		}
	}
	dlist_init(&args->list, 1, 1024, 0);
	dlist_push_back(&args->list, &flag, sizeof(char[2]));
#ifdef USE_AC
	ac_enc_init(&ec, &args->list);
	init_stats(args->hist, args->stats, args->tlevels, args->clevels, image->nch);
#elif defined USE_GRABAC
	ac_enc_init(&ec, &args->list);
	for(int k=0;k<_countof(args->p0);++k)
		args->p0[k]=0x8000;
	for(int k=0;k<_countof(args->p0_bypass);++k)
		args->p0_bypass[k]=0x8000;
#else
	gr_enc_init(&ec, &args->list);
#endif
	memset(args->pixels, 0, args->bufsize);
	for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
	{
		//static const int perm[]={1, 2, 0};

		ALIGN(16) int *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*7*4,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*7*4,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*7*4,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*7*4,
		};
		int yuv[4]={0};
		int wp_preds[8]={0}, wp_result=0;//wp_preds are reused in WG
		int wg_errors[8]={0};
#ifdef ENABLE_WG2
		int wg2_eprev[WG2_NPREDS]={0}, wg2_preds[WG2_NPREDS]={0}, wg2pred=0;
		int depths[]=
		{
			image->depth+depth_inf[combination[0]],
			image->depth+depth_inf[combination[1]],
			image->depth+depth_inf[combination[2]],
		};
#endif
#ifdef ENABLE_WG3
		int wg3_eprev[WG3_NPREDS]={0}, wg3_preds[WG3_NPREDS]={0};
#endif
		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			int
				*NNN	=rows[3]+0*7*4,
				*NN	=rows[2]+0*7*4,
				*NNE	=rows[2]+1*7*4,
				*NW	=rows[1]-1*7*4,
				*N	=rows[1]+0*7*4,
				*NE	=rows[1]+1*7*4,
				*NEE	=rows[1]+2*7*4,
				*NEEE	=rows[1]+3*7*4,
				*WWW	=rows[0]-3*7*4,
				*WW	=rows[0]-2*7*4,
				*W	=rows[0]-1*7*4,
				*curr	=rows[0]+0*7*4;
			//          NNN
			//          NN  NNE
			//       NW N   NE  NEE  NEEE
			//WWW WW W  ?
#ifdef FOLD_OOB
			if(ky<=args->y1+2)
			{
				if(ky<=args->y1+1)
				{
					if(ky==args->y1)
						NEEE=NE=NW=N=W;
					NN=N;
					NNE=NE;
				}
				NNN=NN;
			}
			if(kx<=2)
			{
				if(kx<=1)
				{
					if(!kx)
						NW=W=N;
					WW=W;
				}
				WWW=WW;
			}
			if(kx>=image->iw-3)
			{
				if(kx>=image->iw-2)
				{
					if(kx>=image->iw-1)
					{
						NNE=NN;
						NE=N;
					}
					NEE=NE;
				}
				NEEE=NEE;
			}
#endif

			//if(ky==330&&kx==146)//
			//if(ky==28&&kx==624)//
			//if(idx==26147568)//
			//if(idx==2304)//
			//if(ky==0&&kx==64)//
			//if(ky==1634&&kx==433)//
			//	printf("");

			switch(bestrct)//forward RCT
			{
			case RCT_NONE:
				yuv[0]=image->data[idx+0];
				yuv[1]=image->data[idx+1];
				yuv[2]=image->data[idx+2];
				break;
			case RCT_SUBG:
				yuv[0]=image->data[idx+1];
				yuv[1]=image->data[idx+2];
				yuv[2]=image->data[idx+0];
				break;
			case RCT_JPEG2000_G://r-=g; b-=g; g+=(r+b)>>2;
				yuv[0]=image->data[idx+1];
				yuv[1]=image->data[idx+2];
				yuv[2]=image->data[idx+0];
				yuv[1]-=yuv[0];
				yuv[2]-=yuv[0];
				yuv[0]+=(yuv[1]+yuv[2])>>2;
				break;
			case RCT_JPEG2000_B://r-=b; g-=b; b+=(r+g)>>2;
				yuv[0]=image->data[idx+2];
				yuv[1]=image->data[idx+0];
				yuv[2]=image->data[idx+1];
				yuv[1]-=yuv[0];
				yuv[2]-=yuv[0];
				yuv[0]+=(yuv[1]+yuv[2])>>2;
				break;
			case RCT_JPEG2000_R://b-=r; g-=r; r+=(b+g)>>2;
				yuv[0]=image->data[idx+0];
				yuv[1]=image->data[idx+1];
				yuv[2]=image->data[idx+2];
				yuv[1]-=yuv[0];
				yuv[2]-=yuv[0];
				yuv[0]+=(yuv[1]+yuv[2])>>2;
				break;
			case RCT_1://r-=g; g+=r>>1; b-=g; g+=b>>1;
				yuv[0]=image->data[idx+1];
				yuv[1]=image->data[idx+2];
				yuv[2]=image->data[idx+0];
				yuv[2]-=yuv[0];
				yuv[0]+=yuv[2]>>1;
				yuv[1]-=yuv[0];
				yuv[0]+=yuv[1]>>1;
				break;
			case RCT_PEI09://b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
				yuv[0]=image->data[idx+1];
				yuv[1]=image->data[idx+2];
				yuv[2]=image->data[idx+0];
				yuv[1]-=(87*yuv[2]+169*yuv[0]+128)>>8;
				yuv[2]-=yuv[0];
				yuv[0]+=(86*yuv[2]+29*yuv[1]+128)>>8;
				break;
			}
			for(int kc=0;kc<3;++kc)
			{
				int kc2=kc*7, pred=0, offset=0, ch=combination[kc], val, sym;
#ifdef USE_T47CTX
				int
					vx=abs(W[kc2+1]-WW[kc2+1])+abs(N[kc2+1]-NW[kc2+1])+abs(NE[kc2+1]-N[kc2+1])+abs(WW[kc2+1]+W[kc2+1]),
					vy=abs(W[kc2+1]-NW[kc2+1])+abs(N[kc2+1]-NN[kc2+1])+abs(NE[kc2+1]-NNE[kc2+1])+abs(NN[kc2+1]+N[kc2+1]);
#endif
				switch(predidx[kc])
				{
#ifdef ENABLE_ZERO
				case PRED_ZERO:
					pred=0;
					break;
#endif
				//case PRED_N:
				//	pred=N[kc2];
				//	break;
				//case PRED_W:
				//	pred=W[kc2];
				//	break;
				//case PRED_AV2:
				//	pred=(N[kc2]+W[kc2])/2;
				//	break;
				//case PRED_AV4:
				//	pred=(4*(N[kc2]+W[kc2])+NE[kc2]-NW[kc2])/8;
				//	break;
				//case PRED_AV5:
				//	CLAMP3_32(pred, W[kc2]+(5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])/8, N[kc2], W[kc2], NE[kc2]);
				//	break;
				case PRED_CG:
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					break;
				//case PRED_CGU:
				//	{
				//		int
				//			vx=abs(W[kc2+1]-WW[kc2+1])+abs(N[kc2+1]-NW[kc2+1])+abs(NE[kc2+1]-N[kc2+1])+1,
				//			vy=abs(W[kc2+1]-NW[kc2+1])+abs(N[kc2+1]-NN[kc2+1])+abs(NE[kc2+1]-NNE[kc2+1])+1;
				//		int update=(vx<vy?N[kc2+1]:W[kc2+1])/4;
				//		MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
				//		CLAMP2_32(pred, pred+update, -halfs[ch], halfs[ch]-1);
				//	}
				//	break;
#ifdef ENABLE_WP
				case PRED_WP:
					wp_result=wp_predict(ch_params[ch], NN[kc2], NW[kc2], N[kc2], NE[kc2], W[kc2], NW+kc2+2, N+kc2+2, NE+kc2+2, W+kc2+2, wp_preds);
					pred=(wp_result+127)>>8;
					break;
#endif
#ifdef ENABLE_WG
				case PRED_WG:
					pred=wg_predict(NNN[kc2], NN[kc2], NNE[kc2], NW[kc2], N[kc2], NE[kc2], NEEE[kc2], WWW[kc2], WW[kc2], W[kc2], wg_errors, wp_preds);
					break;
#endif
#ifdef ENABLE_WG2
				case PRED_WG2:
					wg2pred=wg2_predict(rows, 7*4, kc2, 2, depths[kc], wg2_eprev, wg2_preds);
					pred=wg2pred>>(24-depths[kc]);
					break;
#endif
#ifdef ENABLE_WG3
				case PRED_WG3:
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					pred=wg3_predict(pred, NNN[kc2], NN[kc2], NNE[kc2], NW[kc2], N[kc2], NE[kc2], NEE[kc2], NEEE[kc2], WWW[kc2], WW[kc2], W[kc2], wg3_eprev, wg3_preds);
					break;
#endif
				}
				if(bestrct==RCT_SUBG&&kc>0)
				{
					offset=rows[0][0];//luma
					pred+=offset;
					CLAMP2_32(pred, pred, -halfs[ch], halfs[ch]-1);
				}
				val=yuv[kc]-pred;
				{
					int upred=halfs[ch]-abs(pred), aval=abs(val);
					if(aval<=upred)
					{
						sym=val;
						//negmask=-((pr->sse_corr<0)&(sym!=-pr->half[kc]));//sign is flipped if SSE correction was negative, to skew the histogram
						//sym^=negmask;
						//sym-=negmask;
						sym=sym<<1^(sym>>31);//pack sign
					}
					else
						sym=upred+aval;//error sign is known
				}
#ifdef USE_AC
#ifdef CTX_MIN_NW
				int cidx=cdfstride*(nctx*kc+MINVAR(N[kc2+1], W[kc2+1]));
#else
#ifdef USE_T47CTX
				int qeN=QUANTIZE_CTX(vy);
				int qeW=QUANTIZE_CTX(vx);
#else
				int qeN=QUANTIZE_CTX(N[kc2+1]);
				int qeW=QUANTIZE_CTX(W[kc2+1]);
#endif
				int cidx=cdfstride*(nctx*kc+args->clevels*qeN+qeW);
#endif
				int *curr_hist=args->hist+cidx;
				unsigned *curr_CDF=args->stats+cidx;
				quantize_pixel(sym, &token, &bypass, &nbits);
#ifdef _DEBUG
				if((unsigned)token>=(unsigned)args->tlevels)
					LOG_ERROR("Token OOB %d/%d", token, args->tlevels);
#endif
				//if(ky==1634&&kx==432&&kc==0)//
				//	printf("");

				ac_enc(&ec, token, curr_CDF);
				if(nbits)
					ac_enc_bypass(&ec, bypass, nbits);//up to 16 bits
				
#ifdef PROFILE_BYPASS
				args->bypasssizes[kc]+=nbits;
#endif
				//if(curr_hist[token]+1>0x1800)//X
				if(curr_hist[args->tlevels]+1>0x1800)
				{
					int sum=0;
					//update_CDF(curr_hist, curr_CDF, args->tlevels);
					for(int k=0;k<args->tlevels;++k)
						sum+=curr_hist[k]>>=1;
					curr_hist[args->tlevels]=sum;
				}
				++curr_hist[token];
				++curr_hist[args->tlevels];
				if((kx&(CDFSTRIDE-1))==CDFSTRIDE-1)		//FIXME this updates just the current CDF at regular intervals
					update_CDF(curr_hist, curr_CDF, args->tlevels);
#ifdef CTX_MIN_NW
				curr[kc2+1]=token;
#else
				curr[kc2+1]=val;
#endif
#else
				int nbits=FLOOR_LOG2(W[kc2+1]+1);
#ifdef PROFILE_BYPASS
				args->bypasssizes[kc]+=nbits;
#endif
				//if(ky<4)
				//{
				//	int nzeros=sym>>nbits, bypass=sym&((1<<nbits)-1);
				//	for(int k=0;k<nzeros;++k)
				//		printf("0");
				//	printf("1\t");
				//	for(int k=nbits-1;k>=0;--k)
				//		printf("%d", bypass>>k&1);
				//	printf("\n");
				//}
#ifdef USE_GRABAC
				int nzeros=sym>>nbits, bypass=sym&((1<<nbits)-1);
				unsigned short
					*p0_nzeros	=args->p0+((size_t)kc<<4),
					*p0_end		=args->p0+((kc+1LL)<<4)-1,
					*p0_bypass	=args->p0_bypass+kc;
				//if(ky==args->y2-1&&!kx)//
				//	printf("");
				for(int k=0;k<nzeros+1;++k)
				{
					int p0=*p0_nzeros;
					int bit=k>=nzeros;
					ac_enc_bin(&ec, p0, bit);
					UPDATE_P0(p0, bit, LR_SHIFT_TOKEN);
					CLAMP2_32(*p0_nzeros, p0, 1, 0xFFFF);
					p0_nzeros+=p0_nzeros<p0_end;
				}
				for(int k=nbits-1;k>=0;--k)
				{
					int p0=*p0_bypass;
					int bit=bypass>>k&1;
					ac_enc_bin(&ec, p0, bit);
					UPDATE_P0(p0, bit, LR_SHIFT_BYPASS);
					CLAMP2_32(*p0_bypass, p0, 1, 0xFFFF);
				}
#else
				gr_enc_POT(&ec, sym, nbits);
#endif
				curr[kc2+1]=(2*W[kc2+1]+sym+NEEE[kc2+1])/4;			//1211852242
				//curr[kc2+1]=(W[kc2+1]+sym+NEE[kc2+1]+NEEE[kc2+1])/4;		//1212314674
				//curr[kc2+1]=(W[kc2+1]+sym+NEEE[kc2+1])/3;			//1213062186
				//curr[kc2+1]=(3*W[kc2+1]+2*sym+N[kc2+1]+2*NEEE[kc2+1])/8;	//1215279290
				//curr[kc2+1]=(W[kc2+1]+2*sym+NEEE[kc2+1])/4;			//1227447530
#endif
				curr[kc2+0]=yuv[kc]-offset;
#ifdef ENABLE_WP
				if(predidx[kc]==PRED_WP)
					wp_update(curr[kc2+0], wp_result, wp_preds, curr+kc2+2, NE+kc2+2);
#endif
#ifdef ENABLE_WG
				if(predidx[kc]==PRED_WG)
					wg_update(curr[kc2+0], wp_preds, wg_errors);
#endif
#ifdef ENABLE_WG2
				if(predidx[kc]==PRED_WG2)
				{
					int curr24=curr[kc2+0]<<(24-depths[kc]);
					curr[kc2+2]=curr24-wg2pred;
					wg2_update(curr24, wg2_preds, wg2_eprev);
				}
#endif
#ifdef ENABLE_WG3
				if(predidx[kc]==PRED_WG3)
					wg3_update(curr[kc2+0], wg3_preds, wg3_eprev);
#endif
			}
			rows[0]+=7*4;
			rows[1]+=7*4;
			rows[2]+=7*4;
			rows[3]+=7*4;
		}
	}
#if defined USE_AC || defined USE_GRABAC
	ac_enc_flush(&ec);
#else
	gr_enc_flush(&ec);
#endif
	if(args->loud)
		printf("Actual %8zd bytes (%+11.2lf)\n\n", args->list.nobj, args->list.nobj-bestsize);
}
static void block_dec(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
#ifdef USE_AC
	ArithmeticCoder ec;
	int chalf=args->clevels>>1;
#ifdef CTX_MIN_NW
	int nctx=args->clevels;
#else
	int nctx=args->clevels*args->clevels;
#endif
	int cdfstride=args->tlevels+1;
	int token=0, bypass=0, nbits=0;
#elif defined USE_GRABAC
	ArithmeticCoder ec;
#else
	GolombRiceCoder ec;
#endif
	Image *image=args->dst;
	const unsigned char *srcstart=args->decstart, *srcend=args->decend;
	int flag=0, bestrct;
	int combination[3]={0}, predidx[3]={0};
	int halfs[NCHPOOL]={0};
#ifdef ENABLE_WP
	const short *ch_params[NCHPOOL]={0};
#endif
	for(int k=0;k<NCHPOOL;++k)
	{
		halfs[k]=1<<(image->depth+depth_inf[k])>>1;
#ifdef ENABLE_WP
		ch_params[k]=wp_params+WP_NPARAMS*wp_param_idx[k];
#endif
	}
	memcpy(&flag, srcstart, sizeof(char[2]));
	srcstart+=sizeof(char[2]);
	{
		int f2=flag;
		predidx[0]=f2%PRED_COUNT;	f2/=PRED_COUNT;
		predidx[1]=f2%PRED_COUNT;	f2/=PRED_COUNT;
		predidx[2]=f2%PRED_COUNT;	f2/=PRED_COUNT;
		bestrct=f2;
		if((unsigned)bestrct>=(unsigned)RCT_COUNT)
		{
			LOG_ERROR("Corrupt file");
			return;
		}
		memcpy(combination, combinations+bestrct*3, sizeof(int[3]));
	}
#ifdef USE_AC
	ac_dec_init(&ec, srcstart, srcend);
	init_stats(args->hist, args->stats, args->tlevels, args->clevels, image->nch);
#elif defined USE_GRABAC
	ac_dec_init(&ec, srcstart, srcend);
	for(int k=0;k<_countof(args->p0);++k)
		args->p0[k]=0x8000;
	for(int k=0;k<_countof(args->p0_bypass);++k)
		args->p0_bypass[k]=0x8000;
#else
	gr_dec_init(&ec, srcstart, srcend);
#endif
	memset(args->pixels, 0, args->bufsize);
	for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
	{
		//static const int perm[]={1, 2, 0};

		ALIGN(16) int *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*7*4,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*7*4,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*7*4,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*7*4,
		};
		int yuv[4]={0};
		int wp_preds[8]={0}, wp_result=0;//wp_preds are reused in WG
		int wg_errors[8]={0};
#ifdef ENABLE_WG2
		int wg2_eprev[WG2_NPREDS]={0}, wg2_preds[WG2_NPREDS]={0}, wg2pred=0;
		int depths[]=
		{
			image->depth+depth_inf[combination[0]],
			image->depth+depth_inf[combination[1]],
			image->depth+depth_inf[combination[2]],
		};
#endif
#ifdef ENABLE_WG3
		int wg3_eprev[WG3_NPREDS]={0}, wg3_preds[WG3_NPREDS]={0};
#endif
		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			int
				*NNN	=rows[3]+0*7*4,
				*NN	=rows[2]+0*7*4,
				*NNE	=rows[2]+1*7*4,
				*NW	=rows[1]-1*7*4,
				*N	=rows[1]+0*7*4,
				*NE	=rows[1]+1*7*4,
				*NEE	=rows[1]+2*7*4,
				*NEEE	=rows[1]+3*7*4,
				*WWW	=rows[0]-3*7*4,
				*WW	=rows[0]-2*7*4,
				*W	=rows[0]-1*7*4,
				*curr	=rows[0]+0*7*4;
			//          NNN
			//          NN  NNE
			//       NW N   NE  NEE  NEEE
			//WWW WW W  ?
#ifdef FOLD_OOB
			if(ky<=args->y1+2)
			{
				if(ky<=args->y1+1)
				{
					if(ky==args->y1)
						NEEE=NE=NW=N=W;
					NN=N;
					NNE=NE;
				}
				NNN=NN;
			}
			if(kx<=2)
			{
				if(kx<=1)
				{
					if(!kx)
						NW=W=N;
					WW=W;
				}
				WWW=WW;
			}
			if(kx>=image->iw-3)
			{
				if(kx>=image->iw-2)
				{
					if(kx>=image->iw-1)
					{
						NNE=NN;
						NE=N;
					}
					NEE=NE;
				}
				NEEE=NEE;
			}
#endif

			//if(ky==330&&kx==146)//
			//if(ky==28&&kx==624)//
			//if(idx==26147568)//
			//if(idx==2304)//
			//if(ky==0&&kx==64)//
			//if(ky==1634&&kx==433)//
			//	printf("");

			for(int kc=0;kc<3;++kc)
			{
				int kc2=kc*7, pred=0, offset=0, ch=combination[kc], val, sym;
#ifdef USE_T47CTX
				int
					vx=abs(W[kc2+1]-WW[kc2+1])+abs(N[kc2+1]-NW[kc2+1])+abs(NE[kc2+1]-N[kc2+1])+abs(WW[kc2+1]+W[kc2+1]),
					vy=abs(W[kc2+1]-NW[kc2+1])+abs(N[kc2+1]-NN[kc2+1])+abs(NE[kc2+1]-NNE[kc2+1])+abs(NN[kc2+1]+N[kc2+1]);
#endif
				switch(predidx[kc])//WP
				{
#ifdef ENABLE_ZERO
				case PRED_ZERO:
					pred=0;
					break;
#endif
				//case PRED_N:
				//	pred=N[kc2];
				//	break;
				//case PRED_W:
				//	pred=W[kc2];
				//	break;
				//case PRED_AV2:
				//	pred=(N[kc2]+W[kc2])/2;
				//	break;
				//case PRED_AV4:
				//	pred=(4*(N[kc2]+W[kc2])+NE[kc2]-NW[kc2])/8;
				//	break;
				//case PRED_AV5:
				//	CLAMP3_32(pred, W[kc2]+(5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])/8, N[kc2], W[kc2], NE[kc2]);
				//	break;
				case PRED_CG:
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					break;
				//case PRED_CGU:
				//	{
				//		int
				//			vx=abs(W[kc2+1]-WW[kc2+1])+abs(N[kc2+1]-NW[kc2+1])+abs(NE[kc2+1]-N[kc2+1])+1,
				//			vy=abs(W[kc2+1]-NW[kc2+1])+abs(N[kc2+1]-NN[kc2+1])+abs(NE[kc2+1]-NNE[kc2+1])+1;
				//		int update=(vx<vy?N[kc2+1]:W[kc2+1])/4;
				//		MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
				//		CLAMP2_32(pred, pred+update, -halfs[ch], halfs[ch]-1);
				//	}
				//	break;
#ifdef ENABLE_WP
				case PRED_WP:
					wp_result=wp_predict(ch_params[ch], NN[kc2], NW[kc2], N[kc2], NE[kc2], W[kc2], NW+kc2+2, N+kc2+2, NE+kc2+2, W+kc2+2, wp_preds);
					pred=(wp_result+127)>>8;
					break;
#endif
#ifdef ENABLE_WG
				case PRED_WG:
					pred=wg_predict(NNN[kc2], NN[kc2], NNE[kc2], NW[kc2], N[kc2], NE[kc2], NEEE[kc2], WWW[kc2], WW[kc2], W[kc2], wg_errors, wp_preds);
					break;
#endif
#ifdef ENABLE_WG2
				case PRED_WG2:
					wg2pred=wg2_predict(rows, 7*4, kc2, 2, image->depth+depth_inf[ch], wg2_eprev, wg2_preds);
					pred=wg2pred>>(24-depths[kc]);
					break;
#endif
#ifdef ENABLE_WG3
				case PRED_WG3:
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					pred=wg3_predict(pred, NNN[kc2], NN[kc2], NNE[kc2], NW[kc2], N[kc2], NE[kc2], NEE[kc2], NEEE[kc2], WWW[kc2], WW[kc2], W[kc2], wg3_eprev, wg3_preds);
					break;
#endif
				}
				if(bestrct==RCT_SUBG&&kc>0)
				{
					offset=rows[0][0];//luma
					pred+=offset;
					CLAMP2_32(pred, pred, -halfs[ch], halfs[ch]-1);
				}
#ifdef USE_AC
#ifdef CTX_MIN_NW
				int cidx=cdfstride*(nctx*kc+MINVAR(N[kc2+1], W[kc2+1]));
#else
#ifdef USE_T47CTX
				int qeN=QUANTIZE_CTX(vy);
				int qeW=QUANTIZE_CTX(vx);
#else
				int qeN=QUANTIZE_CTX(N[kc2+1]);
				int qeW=QUANTIZE_CTX(W[kc2+1]);
#endif
				int cidx=cdfstride*(nctx*kc+args->clevels*qeN+qeW);
#endif
				int *curr_hist=args->hist+cidx;
				unsigned *curr_CDF=args->stats+cidx;
				
				//if(ky==1634&&kx==432&&kc==0)//
				//	printf("");

				token=ac_dec(&ec, curr_CDF, args->tlevels);//try ac_dec_packedsign()
				sym=token;
				if(sym>=(1<<CONFIG_EXP))
				{
					sym-=1<<CONFIG_EXP;
					int lsb=sym&((1<<CONFIG_LSB)-1);
					sym>>=CONFIG_LSB;
					int msb=sym&((1<<CONFIG_MSB)-1);
					sym>>=CONFIG_MSB;
					nbits=sym+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB);
					bypass=ac_dec_bypass(&ec, nbits);
					sym=1;
					sym<<=CONFIG_MSB;
					sym|=msb;
					sym<<=nbits;
					sym|=bypass;
					sym<<=CONFIG_LSB;
					sym|=lsb;
				}
#elif defined USE_GRABAC
				int nbits=FLOOR_LOG2(W[kc2+1]+1), nzeros=0, bypass=0;
				unsigned short
					*p0_nzeros	=args->p0+((size_t)kc<<4),
					*p0_end		=args->p0+((kc+1LL)<<4)-1,
					*p0_bypass	=args->p0_bypass+kc;
				for(;;)
				{
					int p0=*p0_nzeros;
					int bit=ac_dec_bin(&ec, p0);
					UPDATE_P0(p0, bit, LR_SHIFT_TOKEN);
					//p0+=((bit<<16)-p0)>>LR_SHIFT_TOKEN;
					CLAMP2_32(*p0_nzeros, p0, 1, 0xFFFF);
					if(bit)
						break;
					++nzeros;
					p0_nzeros+=p0_nzeros<p0_end;
				}
				for(int k=nbits-1;k>=0;--k)
				{
					int p0=*p0_bypass;
					int bit=ac_dec_bin(&ec, p0);
					bypass|=bit<<k;
					UPDATE_P0(p0, bit, LR_SHIFT_BYPASS);
					//p0+=((bit<<16)-p0)>>LR_SHIFT_BYPASS;
					CLAMP2_32(*p0_bypass, p0, 1, 0xFFFF);
				}
				sym=nzeros<<nbits|bypass;
#else
				sym=gr_dec_POT(&ec, FLOOR_LOG2(W[kc2+1]+1));
#endif
				{
					int upred=halfs[ch]-abs(pred), negmask=0;
					if(sym<=(upred<<1))
					{
						val=sym>>1^-(sym&1);
						//negmask=-((pr->sse_corr<0)&(error!=-pr->half[kc]));
					}
					else
					{
						val=sym-upred;
						negmask=-(pred>0);
					}
					val^=negmask;
					val-=negmask;
				}
#ifdef USE_AC
				//if(curr_hist[token]+1>0x1800)//X
				if(curr_hist[args->tlevels]+1>0x1800)
				{
					int sum=0;
					//update_CDF(curr_hist, curr_CDF, args->tlevels);
					for(int k=0;k<args->tlevels;++k)
						sum+=curr_hist[k]>>=1;
					curr_hist[args->tlevels]=sum;
				}
				++curr_hist[token];
				++curr_hist[args->tlevels];
				if((kx&(CDFSTRIDE-1))==CDFSTRIDE-1)
					update_CDF(curr_hist, curr_CDF, args->tlevels);
#ifdef CTX_MIN_NW
				curr[kc2+1]=token;
#else
				curr[kc2+1]=val;
#endif
#else
				curr[kc2+1]=(2*W[kc2+1]+sym+NEEE[kc2+1])/4;
				//curr[kc2+1]=(W[kc2+1]+sym+NEE[kc2+1]+NEEE[kc2+1])/4;
				//curr[kc2+1]=(W[kc2+1]+sym+NEEE[kc2+1])/3;
				//curr[kc2+1]=(3*W[kc2+1]+2*sym+N[kc2+1]+2*NEEE[kc2+1])/8;
				//curr[kc2+1]=(W[kc2+1]+2*sym+NEEE[kc2+1])/4;
#endif
				val+=pred;
				yuv[kc]=val;
				curr[kc2+0]=yuv[kc]-offset;
#ifdef ENABLE_WP
				if(predidx[kc]==PRED_WP)
					wp_update(curr[kc2+0], wp_result, wp_preds, curr+kc2+2, NE+kc2+2);
#endif
#ifdef ENABLE_WG
				if(predidx[kc]==PRED_WG)
					wg_update(curr[kc2+0], wp_preds, wg_errors);
#endif
#ifdef ENABLE_WG2
				if(predidx[kc]==PRED_WG2)
				{
					int curr24=curr[kc2+0]<<(24-depths[kc]);
					curr[kc2+2]=curr24-wg2pred;
					wg2_update(curr24, wg2_preds, wg2_eprev);
				}
#endif
#ifdef ENABLE_WG3
				if(predidx[kc]==PRED_WG3)
					wg3_update(curr[kc2+0], wg3_preds, wg3_eprev);
#endif
			}
			switch(bestrct)//forward RCT
			{
			case RCT_NONE:
				image->data[idx+0]=yuv[0];
				image->data[idx+1]=yuv[1];
				image->data[idx+2]=yuv[2];
				break;
			case RCT_SUBG:
				image->data[idx+1]=yuv[0];
				image->data[idx+2]=yuv[1];
				image->data[idx+0]=yuv[2];
				break;
			case RCT_JPEG2000_G://r-=g; b-=g; g+=(r+b)>>2;
				yuv[0]-=(yuv[1]+yuv[2])>>2;
				yuv[2]+=yuv[0];
				yuv[1]+=yuv[0];
				image->data[idx+1]=yuv[0];
				image->data[idx+2]=yuv[1];
				image->data[idx+0]=yuv[2];
				break;
			case RCT_JPEG2000_B://r-=b; g-=b; b+=(r+g)>>2;
				yuv[0]-=(yuv[1]+yuv[2])>>2;
				yuv[2]+=yuv[0];
				yuv[1]+=yuv[0];
				image->data[idx+2]=yuv[0];
				image->data[idx+0]=yuv[1];
				image->data[idx+1]=yuv[2];
				break;
			case RCT_JPEG2000_R://b-=r; g-=r; r+=(b+g)>>2;
				yuv[0]-=(yuv[1]+yuv[2])>>2;
				yuv[2]+=yuv[0];
				yuv[1]+=yuv[0];
				image->data[idx+0]=yuv[0];
				image->data[idx+1]=yuv[1];
				image->data[idx+2]=yuv[2];
				break;
			case RCT_1://r-=g; g+=r>>1; b-=g; g+=b>>1;
				yuv[0]-=yuv[1]>>1;
				yuv[1]+=yuv[0];
				yuv[0]-=yuv[2]>>1;
				yuv[2]+=yuv[0];
				image->data[idx+1]=yuv[0];
				image->data[idx+2]=yuv[1];
				image->data[idx+0]=yuv[2];
				break;
			case RCT_PEI09://b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
				yuv[0]-=(86*yuv[2]+29*yuv[1]+128)>>8;
				yuv[2]+=yuv[0];
				yuv[1]+=(87*yuv[2]+169*yuv[0]+128)>>8;
				image->data[idx+1]=yuv[0];
				image->data[idx+2]=yuv[1];
				image->data[idx+0]=yuv[2];
				break;
			}
#ifdef ENABLE_GUIDE
			if(memcmp(image->data+idx, guide->data+idx, sizeof(short)*image->nch))
			{
				short orig[4]={0};
				memcpy(orig, guide->data+idx, image->nch*sizeof(short));
				LOG_ERROR("Guide error XY %d %d", kx, ky);
				printf("");//
			}
#endif
			rows[0]+=7*4;
			rows[1]+=7*4;
			rows[2]+=7*4;
			rows[3]+=7*4;
		}
	}
}
int f24_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	int ncores, nblocks, nthreads;
	int histsize;
#ifdef USE_AC
	int tlevels, clevels, statssize;
#endif

	if(image->nch!=3)
	{
		LOG_ERROR("Expected 3 channels, got %d", image->nch);
		return 1;
	}
#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif
	ncores=query_cpu_cores();
	nblocks=(image->ih+BLOCKSIZE-1)/BLOCKSIZE;
	nthreads=MINVAR(nblocks, ncores);
	ptrdiff_t memusage=0;
	ptrdiff_t argssize=nthreads*sizeof(ThreadArgs);
	ThreadArgs *args=(ThreadArgs*)malloc(argssize);
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memusage+=argssize;
	memset(args, 0, argssize);
#ifdef USE_AC
	{
		int nlevels=2<<image->depth;//chroma-inflated
		int token=0, bypass=0, nbits=0;

		quantize_pixel(nlevels, &token, &bypass, &nbits);
		tlevels=token+1;
#ifdef CTX_MIN_NW
		clevels=token+1;
		statssize=clevels*(tlevels+1)*image->nch*(int)sizeof(int);
#else
#ifdef USE_T47CTX
		clevels=quantize_ctx(nlevels<<1)+1;
#else
		clevels=quantize_ctx(nlevels>>1)<<1|1;
#endif
		statssize=clevels*clevels*(tlevels+1)*image->nch*(int)sizeof(int);
#endif
		if(fwd)
		{
			histsize=sizeof(int[NCHPOOL*PRED_COUNT])<<(image->depth+1);
			UPDATE_MAX(histsize, statssize);
		}
		else
			histsize=statssize;
	}
#else
	histsize=sizeof(int[NCHPOOL*PRED_COUNT])<<(image->depth+1);
#endif
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=src;
		arg->dst=dst;
		arg->bufsize=sizeof(int[4*NCHPOOL*STRIDE])*(image->iw+16LL);//4 padded rows * NCHPOOL * {pixels, hires-errors, 4 pred hires-errors}
		arg->pixels=(int*)_mm_malloc(arg->bufsize, sizeof(__m128i));
#ifdef USE_AC
		arg->tlevels=tlevels;
		arg->clevels=clevels;
		arg->histsize=histsize;
		arg->hist=(int*)malloc(histsize);
		arg->statssize=statssize;
		arg->stats=(unsigned*)malloc(statssize);
#else
		if(fwd)
		{
			arg->histsize=histsize;
			arg->hist=(int*)malloc(arg->histsize);
		}
		else//hist isn't needed in GR decoder
		{
			arg->histsize=0;
			arg->hist=0;
		}
#endif
		if(!arg->pixels||(fwd&&!arg->hist)
#ifdef USE_AC
			||!args->stats
#endif
		)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=(ptrdiff_t)arg->bufsize+arg->histsize;
#ifdef USE_AC
		memusage+=args->statssize;
#endif
		arg->fwd=fwd;
#ifdef DISABLE_MT
		arg->loud=loud;
#else
		arg->loud=0;
#endif
	}
	if(fwd)
	{
		int hist_rct[RCT_COUNT]={0};
		int hist_pred[PRED_COUNT]={0};

#ifdef USE_PALETTE
		int palused;
		Image palimage={0};
		int palnvals[4]={0};
		unsigned short *pal;
#endif
		double bestsize=0;
#ifdef PROFILE_BYPASS
		double bypasstotal=0;
#endif
		ptrdiff_t start_dst, start_blocksizes;
#ifdef USE_PALETTE
		pal=apply_palette_fwd(image, &palimage, palnvals);
		palused=pal!=0&&palimage.data!=0;
		if(palused)
		{
			image=&palimage;
			for(int k=0;k<nthreads;++k)
				args[k].src=&palimage;
#ifdef ENABLE_GUIDE
			image_copy(&g_palimage, &palimage);
			guide=&g_palimage;
#endif
		}

		start_dst=-1;
		for(int kc=0;kc<image->nch;++kc)//data starts with NCH palette sizes
		{
			ptrdiff_t start=array_append(data, palnvals+kc, 1, 1LL<<(image->depth>8), 1, 0, 0);
			if(start_dst==-1)
				start_dst=start;
		}
		if(palused)
		{
			unsigned short *palptr=pal;
			for(int kc=0;kc<image->nch;++kc)//palette can be delta coded
			{
				int count=palnvals[kc];
				if(src->depth>8)
				{
					ptrdiff_t idx=array_append(data, 0, 1, (ptrdiff_t)count<<1, 1, 0, 0);
					unsigned short *dstptr=(unsigned short*)(data[0]->data+idx);
					for(int k=0;k<count;++k)
						memcpy(dstptr+k*sizeof(short), palptr+k, sizeof(short));
				}
				else
				{
					ptrdiff_t idx=array_append(data, 0, 1, count, 1, 0, 0);
					unsigned char *dstptr=data[0]->data+idx;
					for(int k=0;k<count;++k)
						dstptr[k]=(unsigned char)palptr[k];
				}
				palptr+=1LL<<image->depth;
			}
		}

		start_blocksizes=data[0]->count;
		array_append(data, 0, 1, sizeof(int)*nblocks, 1, 0, 0);
#else
		start_blocksizes=start_dst=array_append(data, 0, 1, sizeof(int)*nblocks, 1, 0, 0);
#endif
		for(int kt=0;kt<nblocks;kt+=nthreads)
		{
			int nthreads2=MINVAR(kt+nthreads, nblocks)-kt;
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				arg->y1=BLOCKSIZE*(kt+kt2);
				arg->y2=MINVAR(arg->y1+BLOCKSIZE, image->ih);
			}
#ifdef DISABLE_MT
			for(int k=0;k<nthreads2;++k)
				block_enc(args+k);
#else
			void *ctx=mt_exec(block_enc, args, sizeof(ThreadArgs), nthreads2);
			mt_finish(ctx);
#endif
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				int
					blocksize=(image->iw*(arg->y2-arg->y1)*image->nch*image->depth+7)>>3,
					cbsize=(image->iw*(arg->y2-arg->y1)*image->depth+7)>>3;
				if(loud)
				{
					if(!(kt+kt2))
						printf("block,  nrows,  usize,     best  ->  actual,  (actual-best)\n");
					printf(
						"[%3d]  %4d  %8d  %11.2lf -> %8zd bytes  (%+11.2lf)"
#ifdef PROFILE_BYPASS
						"  bypass %7td+%7td+%7td=%8td (%8.4lf%%)"
#endif
						"  %s %s %s %s\n",
						kt+kt2,
						arg->y2-arg->y1,
						blocksize,
						arg->bestsize,
						arg->list.nobj,
						(double)arg->list.nobj-arg->bestsize,
#ifdef PROFILE_BYPASS
						arg->bypasssizes[0]>>3,
						arg->bypasssizes[1]>>3,
						arg->bypasssizes[2]>>3,
						(arg->bypasssizes[0]+arg->bypasssizes[1]+arg->bypasssizes[2])>>3,
						100.*(arg->bypasssizes[0]+arg->bypasssizes[1]+arg->bypasssizes[2])/(8*arg->list.nobj),
#endif
						rctnames[arg->bestrct],
						prednames[arg->predidx[0]],
						prednames[arg->predidx[1]],
						prednames[arg->predidx[2]]
					);
					bestsize+=arg->bestsize;
					hist_rct[arg->bestrct]+=blocksize;
					hist_pred[arg->predidx[0]]+=cbsize;
					hist_pred[arg->predidx[1]]+=cbsize;
					hist_pred[arg->predidx[2]]+=cbsize;
#ifdef PROFILE_BYPASS
					bypasstotal+=(arg->bypasssizes[0]+arg->bypasssizes[1]+arg->bypasssizes[2])/8.;
#endif
				}
				rct_gains[arg->bestrct]+=blocksize;
				pred_gains[arg->predidx[0]]+=cbsize;
				pred_gains[arg->predidx[1]]+=cbsize;
				pred_gains[arg->predidx[2]]+=cbsize;

				memcpy(data[0]->data+start_blocksizes+sizeof(int)*((ptrdiff_t)kt+kt2), &arg->list.nobj, sizeof(int));
				dlist_appendtoarray(&arg->list, data);
				dlist_clear(&arg->list);
			}
		}
		if(loud)
		{
			ptrdiff_t
				csize=data[0]->count-start_dst,
				usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
			t0=time_sec()-t0;
#ifdef USE_PALETTE
			if(palused)
			{
				int palsize=palnvals[0]+palnvals[1]+palnvals[2]+palnvals[3];
				printf("Palette  %d bytes\n", palsize<<(image->depth>8));
				for(int kc=0;kc<4;++kc)
				{
					if(palnvals[kc])
						printf("  C%d  %5d/%5d levels\n", kc, palnvals[kc], 1<<image->depth);
				}
			}
#endif
			{
				int blocksize=(image->iw*image->nch*image->depth*BLOCKSIZE+7)>>3;
				printf("RCTs:\n");
				for(int k=0;k<RCT_COUNT;++k)
				{
					int nstars=(int)(hist_rct[k]*60LL/usize);
					printf("%s %12d %8.4lf ", rctnames[k], hist_rct[k], (double)hist_rct[k]/blocksize);
					for(int k2=0;k2<nstars;++k2)
						printf("*");
					printf("\n");
				}
				printf("\n");

				printf("Predictors:\n");
				for(int k=0;k<PRED_COUNT;++k)
				{
					int nstars=(int)(hist_pred[k]*60LL/usize);
					printf("%s %12d %8.4lf ", prednames[k], hist_pred[k], (double)hist_pred[k]/blocksize);
					for(int k2=0;k2<nstars;++k2)
						printf("*");
					printf("\n");
				}
				printf("\n");
			}
			printf(
				"best:     %12.2lf\n"
				"actual: [[%9td]]/%9td bytes  (%+11.2lf)  %10lf%% %10lf\n",
				bestsize,
				csize,
				usize,
				csize-bestsize,
				100.*csize/usize,
				(double)usize/csize
			);
#ifdef PROFILE_BYPASS
			printf("bypass %.2lf %lf%%\n", bypasstotal, 100.*bypasstotal/csize);
#endif
			printf("E %16.6lf sec  %16.6lf MB/s  Mem usage: ", t0, usize/(t0*1024*1024));
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
#if 0
			unsigned char *ptr=data[0]->data+start_dst;
			int n=MINVAR((int)data[0]->count, 4000);
			printf("First %d bytes:\n", n);
			for(int k=0;k<n;++k)
				printf("%02X", ptr[k]);
			printf("\n");
#endif
#if 0
			int bins[16]={0};
			for(int k=0;k<(int)data[0]->count;++k)
				++bins[ptr[k]&15];
			printf("Histogram of low 4 bytes:\n");
			for(int k=0;k<16;++k)
			{
				int nstars=bins[k]*240/(int)data[0]->count;
				printf("%2d %8d ", k, bins[k]);
				for(int k2=0;k2<nstars;++k2)
					printf("*");
				printf("\n");
			}
#endif
		}
#ifdef USE_PALETTE
		if(palused)
			image_clear(&palimage);
#endif
	}
	else
	{
		const unsigned char *dstptr;
		int dec_offset=0;
		
#ifdef USE_PALETTE
		//const unsigned char *debugptr=cbuf;
		//ptrdiff_t debugsize=clen;

		int palsize=0;
		int palnvals[4]={0};
		unsigned short *pal=0, *palptr;

		for(int kc=0;kc<dst->nch;++kc)
		{
			int tagsize=1<<(dst->depth>8);
			memcpy(palnvals+kc, cbuf, tagsize);
			cbuf+=tagsize;
			clen-=tagsize;
			palsize+=palnvals[kc];
		}
		if(palsize)
		{
			pal=(unsigned short*)malloc(sizeof(short)*palsize);
			if(!pal)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			palptr=pal;
			for(int kc=0;kc<image->nch;++kc)//palette can be delta coded
			{
				int count=palnvals[kc];
				if(dst->depth>8)
				{
					for(int k=0;k<count;++k)
						memcpy(palptr+k, cbuf+k*sizeof(short), sizeof(short));
					cbuf+=count*sizeof(short);
					clen-=count*sizeof(short);
				}
				else
				{
					for(int k=0;k<count;++k)
						palptr[k]=(unsigned char)cbuf[k];
					cbuf+=count;
					clen-=count;
				}
				palptr+=count;
			}
		}
#endif
		dstptr=cbuf+sizeof(int)*nblocks;
		//integrity check
#if 1
		for(int kt=0;kt<nblocks;++kt)
		{
			int size=0;
			memcpy(&size, cbuf+sizeof(int)*kt, sizeof(int));
			dec_offset+=size;
		}
		if(sizeof(int)*nblocks+dec_offset!=clen)
			LOG_ERROR("Corrupt file");
#endif
		dec_offset=0;
		for(int kt=0;kt<nblocks;kt+=nthreads)
		{
			int nthreads2=MINVAR(kt+nthreads, nblocks)-kt;
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				int size=0;
				memcpy(&size, cbuf+sizeof(int)*((ptrdiff_t)kt+kt2), sizeof(int));
				arg->y1=BLOCKSIZE*(kt+kt2);
				arg->y2=MINVAR(arg->y1+BLOCKSIZE, image->ih);
				arg->decstart=dstptr+dec_offset;
				dec_offset+=size;
				arg->decend=dstptr+dec_offset;
			}
#ifdef DISABLE_MT
			for(int k=0;k<nthreads2;++k)
				block_dec(args+k);
#else
			{
				void *ctx=mt_exec(block_dec, args, sizeof(ThreadArgs), nthreads2);
				mt_finish(ctx);
			}
#endif
		}
#ifdef USE_PALETTE
		if(palsize)
		{
			apply_palette_inv(dst, pal, palnvals);
			free(pal);
#ifdef ENABLE_GUIDE
			image_clear(&g_palimage);
#endif
		}
#endif
		if(loud)
		{
			ptrdiff_t usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
			t0=time_sec()-t0;
			printf("D %16.6lf sec  %16.6lf MB/s  Mem usage: ", t0, usize/(t0*1024*1024));
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
		}
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		_mm_free(arg->pixels);
#ifdef USE_AC
		free(arg->hist);
		free(arg->stats);
#else
		if(fwd)
			free(arg->hist);
#endif
	}
	free(args);
	return 0;
}
void f24_curiosity()
{
	long long sum;

	printf("RCTs:\n");
	sum=0;
	for(int k=0;k<RCT_COUNT;++k)
		sum+=rct_gains[k];
	for(int k=0;k<RCT_COUNT;++k)
	{
		int nstars=(int)(rct_gains[k]*60LL/sum);
		printf("%s %12lld %8.4lf%% ", rctnames[k], rct_gains[k], 100.*rct_gains[k]/sum);
		for(int k2=0;k2<nstars;++k2)
			printf("*");
		printf("\n");
	}
	printf("\n");

	printf("Predictors:\n");
	sum=0;
	for(int k=0;k<PRED_COUNT;++k)
		sum+=pred_gains[k];
	for(int k=0;k<PRED_COUNT;++k)
	{
		int nstars=(int)(pred_gains[k]*60LL/sum);
		printf("%s %12lld %8.4lf%% ", prednames[k], pred_gains[k], 100.*pred_gains[k]/sum);
		for(int k2=0;k2<nstars;++k2)
			printf("*");
		printf("\n");
	}
	printf("\n");
}