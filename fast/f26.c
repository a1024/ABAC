#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
//#ifdef _MSC_VER
//#include<intrin.h>
//#else
//#include<x86intrin.h>
//#endif
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT


#include"ac.h"
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
#define BLOCKSIZE 256
typedef enum _IChannelType
{
	ICH_R,
	ICH_G,
	ICH_B,
	ICH_J2Y,
	ICH_J2U,
	ICH_J2V,
	ICH_R1Y,
	ICH_R1U,
	ICH_P9Y,
	ICH_P9U,

	ICH_ZERO,
	ICH_COUNT,
} IChannelType;
#define ICH_R1V ICH_J2V
#define ICH_P9V ICH_J2V
#define OCHLIST\
	OCH(OCH_R,	0,	ICH_R,		ICH_ZERO)\
	OCH(OCH_G,	0,	ICH_G,		ICH_ZERO)\
	OCH(OCH_B,	0,	ICH_B,		ICH_ZERO)\
	OCH(OCH_RG,	0,	ICH_R,		ICH_G)\
	OCH(OCH_RB,	0,	ICH_R,		ICH_B)\
	OCH(OCH_GR,	0,	ICH_G,		ICH_R)\
	OCH(OCH_GB,	0,	ICH_G,		ICH_B)\
	OCH(OCH_BG,	0,	ICH_B,		ICH_G)\
	OCH(OCH_BR,	0,	ICH_B,		ICH_R)\
	OCH(OCH_J2Y,	0,	ICH_J2Y,	ICH_ZERO)\
	OCH(OCH_J2U,	1,	ICH_J2U,	ICH_ZERO)\
	OCH(OCH_J2V,	1,	ICH_J2V,	ICH_ZERO)\
	OCH(OCH_R1Y,	0,	ICH_R1Y,	ICH_ZERO)\
	OCH(OCH_R1U,	1,	ICH_R1U,	ICH_ZERO)\
	OCH(OCH_P9Y,	0,	ICH_P9Y,	ICH_ZERO)\
	OCH(OCH_P9U,	1,	ICH_P9U,	ICH_ZERO)
typedef enum _OChannelType
{
#define OCH(ONAME, OINF, TARGET, HELPER) ONAME,
	OCHLIST
#undef  OCH
	OCH_COUNT,
} OChannelType;
static const int och_info[OCH_COUNT][3]=
{
#define OCH(ONAME, OINF, TARGET, HELPER) {TARGET, HELPER, OINF},
	OCHLIST
#undef  OCH
};
static const char *och_names[OCH_COUNT]=
{
#define OCH(ONAME, OINF, TARGET, HELPER) #ONAME,
	OCHLIST
#undef  OCH
};
#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		4, 4, 4)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		4, 4, 1)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		4, 4, 0)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		4, 4, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		4, 4, 1)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		4, 4, 1)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		4, 4, 0)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		4, 0, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		4, 0, 1)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		4, 0, 1)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		4, 0, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		4, 0, 1)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		4, 0, 1)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		4, 0, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		4, 0, 1)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		4, 0, 1)\
	RCT(J2K,	OCH_J2Y,	OCH_J2U,	OCH_J2V,	4, 4, 4)\
	RCT(RCT1,	OCH_R1Y,	OCH_R1U,	OCH_J2V,	4, 4, 4)\
	RCT(Pei09,	OCH_P9Y,	OCH_P9U,	OCH_J2V,	4, 4, 4)
typedef enum _RCTType
{
#define RCT(RCTNAME, YIDX, UIDX, VIDX, YOFF, UOFF, VOFF) RCT_##RCTNAME,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTType;
static const int rct_combinations[RCT_COUNT][6]=
{
#define RCT(RCTNAME, YIDX, UIDX, VIDX, YOFF, UOFF, VOFF) {YIDX, UIDX, VIDX, YOFF, UOFF, VOFF},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(RCTNAME, YIDX, UIDX, VIDX, YOFF, UOFF, VOFF) #RCTNAME,
	RCTLIST
#undef  RCT
};

#define PREDLIST\
	PRED(W)\
	PRED(cgrad)\
	PRED(wgrad)
typedef enum _PredType
{
#define PRED(NAME) PRED_##NAME,
	PREDLIST
#undef  PRED
	PRED_COUNT,
} PredType;
static const char *pred_names[PRED_COUNT]=
{
#define PRED(NAME) #NAME,
	PREDLIST
#undef  PRED
};


//WG:

//	#define WG_UPDATE
#define WG_RESCALE_LIMIT 100
#define WG_DECAY_NUM	493
#define WG_DECAY_SH	9

#define WG_NPREDS	8
#define WG_PREDLIST\
	WG_PRED(1.2, N+W-NW)\
	WG_PRED(1.5, N)\
	WG_PRED(1.5, W)\
	WG_PRED(1, W+NE-N)\
	WG_PRED(1, N+NE-NNE)\
	WG_PRED(1, 3*(N-NN)+NNN)\
	WG_PRED(1, 3*(W-WW)+WWW)\
	WG_PRED(1, (W+NEEE)/2)
//	WG_PRED(0.5, NW)
//	WG_PRED(0.5, NE)
static void wg_init(double *weights)
{
	int j=0;
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
	WG_PREDLIST
#undef  WG_PRED
}
static int wg_predict(const double *weights, int NNN, int NN, int NNE, int NW, int N, int NE, int NEE, int NEEE, int WWW, int WW, int W, const int *eW, int *preds)
{
	int pred;
	double pred2=0, wsum=0;
	int j=0;
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST
#undef  WG_PRED
	
	for(int k=0;k<WG_NPREDS;++k)
	{
		double weight=weights[k]/(eW[k]+1);
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
static void wg_update(int curr, const int *preds, int *eW, double *weights)
{
	int kbest=0, ebest=0;
	for(int k=0;k<WG_NPREDS;++k)
	{
		int e2=abs(curr-preds[k]);
		eW[k]=(eW[k]+e2)*WG_DECAY_NUM>>WG_DECAY_SH;
		if(!k||ebest>e2)
			kbest=k, ebest=e2;
	}
#ifdef WG_UPDATE
	++weights[kbest];
	if(weights[kbest]>WG_RESCALE_LIMIT)
	{
		for(int k=0;k<WG_NPREDS;++k)
			weights[k]*=0.5;
	}
#endif
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
		*nbits=lgv-CONFIG_MSB+CONFIG_LSB;
		*bypass=val>>CONFIG_LSB&((1LL<<*nbits)-1);
	}
}

typedef struct _ThreadArgs
{
	const Image *src;
	Image *dst;
	int fwd, loud, y1, y2;
	
	short *pixels;
	int bufsize;

	int *hist, histsize, histindices[OCH_COUNT*PRED_COUNT+1];

	DList list;
	const unsigned char *decstart, *decend;

	double wg_weights[WG_NPREDS*OCH_COUNT];

	//AC
	int tlevels, clevels, statssize;
	unsigned *stats;
} ThreadArgs;
#define CDFSTRIDE 64
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
static void block_thread(void *param)
{
	ArithmeticCoder ec;
	ThreadArgs *args=(ThreadArgs*)param;
	Image const *image=args->fwd?args->src:args->dst;
	int depths[OCH_COUNT]={0}, halfs[OCH_COUNT]={0};
	double bestsize=0;
	int bestrct=0, combination[6]={0}, predidx[3]={0}, flag;
	int res=image->iw*(args->y2-args->y1);
	int nctx=args->clevels*args->clevels, cdfstride=args->tlevels+1, chsize=nctx*cdfstride;

	for(int kc=0;kc<OCH_COUNT;++kc)
	{
		int depth=image->depth+och_info[kc][2];
		UPDATE_MIN(depth, 16);
		depths[kc]=depth;
		halfs[kc]=1<<depth>>1;
	}
	memset(args->pixels, 0, args->bufsize);
		
	if(args->fwd)//encode
	{
		int nlevels[OCH_COUNT]={0};
		double csizes[OCH_COUNT*PRED_COUNT]={0};
		char predsel[OCH_COUNT]={0};

		for(int kc=0;kc<OCH_COUNT;++kc)
			nlevels[kc]=1<<depths[kc];
		for(int kc=0;kc<OCH_COUNT;++kc)
			wg_init(args->wg_weights+WG_NPREDS*kc);
		memset(args->hist, 0, args->histsize);
		for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
		{
			ALIGN(16) short *rows[]=
			{
				args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*OCH_COUNT,
				args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*OCH_COUNT,
				args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*OCH_COUNT,
				args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*OCH_COUNT,
			};
			short input[ICH_COUNT]={0};
			int preds[PRED_COUNT]={0};
			int wg_errors[WG_NPREDS]={0}, wg_preds[WG_NPREDS]={0};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				short
					*NNN	=rows[3]+0*OCH_COUNT,
					*NNW	=rows[2]-1*OCH_COUNT,
					*NN	=rows[2]+0*OCH_COUNT,
					*NNE	=rows[2]+1*OCH_COUNT,
					*NNEE	=rows[2]+2*OCH_COUNT,
					*NW	=rows[1]-1*OCH_COUNT,
					*N	=rows[1]+0*OCH_COUNT,
					*NE	=rows[1]+1*OCH_COUNT,
					*NEE	=rows[1]+2*OCH_COUNT,
					*NEEE	=rows[1]+3*OCH_COUNT,
					*WWW	=rows[0]-3*OCH_COUNT,
					*WW	=rows[0]-2*OCH_COUNT,
					*W	=rows[0]-1*OCH_COUNT,
					*curr	=rows[0]+0*OCH_COUNT;
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
				//NONE
				input[ICH_R]=image->data[idx+0];//r
				input[ICH_G]=image->data[idx+1];//g
				input[ICH_B]=image->data[idx+2];//b

				//JPEG2000	r-=g; b-=g; g+=(r+b)>>2;
				input[ICH_J2V]=input[ICH_R]-input[ICH_G];
				input[ICH_J2U]=input[ICH_B]-input[ICH_G];
				input[ICH_J2Y]=input[ICH_G]+((input[ICH_J2V]+input[ICH_J2U])>>2);

				//RCT1		r-=g; g+=r>>1; b-=g; g+=b>>1;
			//	input[ICH_R1V]=input[ICH_R]-input[ICH_G];
				input[ICH_R1Y]=input[ICH_G]+(input[ICH_R1V]>>1);
				input[ICH_R1U]=input[ICH_B]-input[ICH_G];
				input[ICH_R1Y]+=input[ICH_R1U]>>1;

				//Pei09		b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
				input[ICH_P9U]=input[ICH_B]-((87*input[ICH_R]+169*input[ICH_G]+128)>>8);
			//	input[ICH_P9V]=input[ICH_R]-input[ICH_G];
				input[ICH_P9Y]=input[ICH_G]+((86*input[ICH_P9V]+29*input[ICH_P9U]+128)>>8);

				//if(ky==128&&kx==128)//
				//	printf("");

				for(int kc=0;kc<OCH_COUNT;++kc)
				{
					int target=input[och_info[kc][0]], offset=input[och_info[kc][1]];

					preds[PRED_W]=W[kc];
					MEDIAN3_32(preds[PRED_cgrad], N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
					preds[PRED_wgrad]=wg_predict(
						args->wg_weights+WG_NPREDS*kc,
						NNN[kc],
						NN[kc], NNE[kc],
						NW[kc], N[kc], NE[kc], NEE[kc], NEEE[kc],
						WWW[kc], WW[kc], W[kc],
						wg_errors, wg_preds
					);

					if(offset)
					{
						for(int kp=0;kp<PRED_COUNT;++kp)
						{
							preds[kp]+=offset;
							CLAMP2_32(preds[kp], preds[kp], -halfs[kc], halfs[kc]-1);
						}
					}
					for(int kp=0;kp<PRED_COUNT;++kp)
					{
						int *curr_hist=args->hist+args->histindices[kc*PRED_COUNT+kp];
						int val=target-preds[kp];
						val+=halfs[kc];
						val&=nlevels[kc]-1;
						//if((unsigned)val>=(unsigned)(args->histindices[kc*PRED_COUNT+kp+1]-args->histindices[kc*PRED_COUNT+kp]))//
						//	LOG_ERROR("");
						++curr_hist[val];
					}
					target-=offset;
					wg_update(target, wg_preds, wg_errors, args->wg_weights+WG_NPREDS*kc);
					curr[kc]=target;
				}
				rows[0]+=OCH_COUNT;
				rows[1]+=OCH_COUNT;
				rows[2]+=OCH_COUNT;
				rows[3]+=OCH_COUNT;
			}
		}
		for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)//calculate channel sizes
		{
			int *curr_hist=args->hist+args->histindices[kc];
			int nlevels2=nlevels[kc/PRED_COUNT];
			for(int ks=0;ks<nlevels2;++ks)
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
		for(int kt=0;kt<RCT_COUNT;++kt)//select best R.C.T
		{
			const int *combination=rct_combinations[kt];
			double csize=
				csizes[combination[0]*PRED_COUNT+predsel[combination[0]]]+
				csizes[combination[1]*PRED_COUNT+predsel[combination[1]]]+
				csizes[combination[2]*PRED_COUNT+predsel[combination[2]]];
			if(!kt||bestsize>csize)
				bestsize=csize, bestrct=kt;
		}
		memcpy(combination, rct_combinations[bestrct], sizeof(combination));
		predidx[0]=predsel[combination[0]];
		predidx[1]=predsel[combination[1]];
		predidx[2]=predsel[combination[2]];
		flag=bestrct;
		flag=flag*PRED_COUNT+predidx[2];
		flag=flag*PRED_COUNT+predidx[1];
		flag=flag*PRED_COUNT+predidx[0];//19*3*3*3 = 513
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
				const int *combination=rct_combinations[kt];
				double csize=
					csizes[combination[0]*PRED_COUNT+predsel[combination[0]]]+
					csizes[combination[1]*PRED_COUNT+predsel[combination[1]]]+
					csizes[combination[2]*PRED_COUNT+predsel[combination[2]]];
				printf("%12.2lf %c  %-10s %-10s %-10s %-10s\n",
					csize,
					kt==bestrct?'*':' ',
					rct_names[kt],
					pred_names[predsel[combination[0]]],
					pred_names[predsel[combination[1]]],
					pred_names[predsel[combination[2]]]
				);
			}
		}
		dlist_init(&args->list, 1, 1024, 0);
		dlist_push_back(&args->list, &flag, sizeof(char[2]));
		ac_enc_init(&ec, &args->list);
	}
	else//decode
	{
		const unsigned char *srcstart=args->decstart, *srcend=args->decend;
			
		memcpy(&flag, srcstart, sizeof(char[2]));
		srcstart+=sizeof(char[2]);
		bestrct=flag;
		predidx[0]=bestrct%PRED_COUNT;	bestrct/=PRED_COUNT;
		predidx[1]=bestrct%PRED_COUNT;	bestrct/=PRED_COUNT;
		predidx[2]=bestrct%PRED_COUNT;	bestrct/=PRED_COUNT;
		if((unsigned)bestrct>=(unsigned)RCT_COUNT)
		{
			LOG_ERROR("Corrupt file");
			return;
		}
		memcpy(combination, rct_combinations[bestrct], sizeof(combination));
		ac_dec_init(&ec, srcstart, srcend);
	}
	
	for(int kc=0;kc<image->nch;++kc)
	{
		int *curr_hist=args->hist+chsize*kc;
		unsigned *curr_CDF=args->stats+chsize*kc;
		
		*curr_hist=1;
		memfill(curr_hist+1, curr_hist, sizeof(int)*(args->tlevels-1LL), sizeof(int));
		curr_hist[args->tlevels]=args->tlevels;
		update_CDF(curr_hist, curr_CDF, args->tlevels);
		memfill(curr_hist+cdfstride, curr_hist, ((size_t)chsize-cdfstride)*sizeof(int), cdfstride*sizeof(int));
		memfill(curr_CDF+cdfstride, curr_CDF, ((size_t)chsize-cdfstride)*sizeof(int), cdfstride*sizeof(int));
	}
	for(int kc=0;kc<3;++kc)
	{
		if(predidx[kc]==PRED_wgrad)
			wg_init(args->wg_weights+WG_NPREDS*kc);
	}
	for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) short *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*4,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*4,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*4,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*4,
		};
		short yuv[5]={0};
		int wg_errors[WG_NPREDS]={0}, wg_preds[WG_NPREDS]={0};
		int token=0, bypass=0, nbits=0;

		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			short
				*NNN	=rows[3]+0*4,
				*NNW	=rows[2]-1*4,
				*NN	=rows[2]+0*4,
				*NNE	=rows[2]+1*4,
				*NNEE	=rows[2]+2*4,
				*NW	=rows[1]-1*4,
				*N	=rows[1]+0*4,
				*NE	=rows[1]+1*4,
				*NEE	=rows[1]+2*4,
				*NEEE	=rows[1]+3*4,
				*WWW	=rows[0]-3*4,
				*WW	=rows[0]-2*4,
				*W	=rows[0]-1*4,
				*curr	=rows[0]+0*4;
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
			if(args->fwd)
			{
				int kc=0;
				if(image->nch>=3)
				{
					switch(bestrct)
					{
					case RCT_R_G_B:
					case RCT_R_G_BG:
					case RCT_R_G_BR:
					case RCT_R_GR_BR:
					case RCT_R_GR_BG:
						yuv[0]=image->data[idx+0];
						yuv[1]=image->data[idx+1];
						yuv[2]=image->data[idx+2];
						break;
					case RCT_G_B_RG:
					case RCT_G_B_RB:
					case RCT_G_BG_RG:
					case RCT_G_BG_RB:
						yuv[0]=image->data[idx+1];
						yuv[1]=image->data[idx+2];
						yuv[2]=image->data[idx+0];
						break;
					case RCT_B_R_GR:
					case RCT_B_R_GB:
					case RCT_B_RB_GB:
					case RCT_B_RB_GR:
						yuv[0]=image->data[idx+2];
						yuv[1]=image->data[idx+0];
						yuv[2]=image->data[idx+1];
						break;
					case RCT_G_RG_BR:
						yuv[0]=image->data[idx+1];
						yuv[1]=image->data[idx+0];
						yuv[2]=image->data[idx+2];
						break;
					case RCT_B_GB_RG:
						yuv[0]=image->data[idx+2];
						yuv[1]=image->data[idx+1];
						yuv[2]=image->data[idx+0];
						break;
					case RCT_R_BR_GB:
						yuv[0]=image->data[idx+0];
						yuv[1]=image->data[idx+2];
						yuv[2]=image->data[idx+1];
						break;
					case RCT_J2K:
						yuv[0]=image->data[idx+1];
						yuv[1]=image->data[idx+2];
						yuv[2]=image->data[idx+0];
						yuv[1]-=yuv[0];
						yuv[2]-=yuv[0];
						yuv[0]+=(yuv[1]+yuv[2])>>2;
						break;
					case RCT_RCT1:
						yuv[0]=image->data[idx+1];
						yuv[1]=image->data[idx+2];
						yuv[2]=image->data[idx+0];
						yuv[1]-=yuv[0];
						yuv[0]+=yuv[1]>>1;
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						break;
					case RCT_Pei09:
						yuv[0]=image->data[idx+1];
						yuv[1]=image->data[idx+2];
						yuv[2]=image->data[idx+0];
						yuv[1]-=(87*yuv[2]+169*yuv[0]+128)>>8;
						yuv[2]-=yuv[0];
						yuv[0]+=(86*yuv[2]+29*yuv[1]+128)>>8;
						break;
					}
					kc=3;
				}
				if(kc<image->nch)
					yuv[kc]=image->data[kc];
			}

			//if(!idx)//
			//	printf("");
			
			for(int kc=0;kc<image->nch;++kc)
			{
				int
					vx=(abs(W[kc]-WW[kc])+abs(N[kc]-NW[kc])+abs(NE[kc]-N[kc]))<<9>>depths[kc],
					vy=(abs(W[kc]-NW[kc])+abs(N[kc]-NN[kc])+abs(NE[kc]-NNE[kc]))<<9>>depths[kc];
				int qeN=FLOOR_LOG2_P1(vy);
				int qeW=FLOOR_LOG2_P1(vx);
				qeN=MINVAR(qeN, 8);
				qeW=MINVAR(qeW, 8);
				int cidx=cdfstride*(nctx*kc+args->clevels*qeN+qeW);
				int *curr_hist=args->hist+cidx;
				unsigned *curr_CDF=args->stats+cidx;
				int ch=combination[kc], offset=yuv[combination[kc+3]], pred, error, sym;

				switch(predidx[kc])
				{
				case PRED_W:
					pred=W[kc];
					break;
				case PRED_cgrad:
					MEDIAN3_32(pred, N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
					break;
				case PRED_wgrad:
					pred=wg_predict(
						args->wg_weights+WG_NPREDS*kc,
						NNN[kc],
						NN[kc], NNE[kc],
						NW[kc], N[kc], NE[kc], NEE[kc], NEEE[kc],
						WWW[kc], WW[kc], W[kc],
						wg_errors, wg_preds
					);
					break;
				}
				pred+=offset;
				CLAMP2_32(pred, pred, -halfs[ch], halfs[ch]-1);
				if(args->fwd)
				{
					curr[kc]=yuv[kc];
					error=curr[kc]-pred;
					{
						int upred=halfs[ch]-abs(pred), aval=abs(error);
						if(aval<=upred)
						{
							sym=error;
#ifdef ENABLE_BIASCORR
							{
								int negmask=-((ibias_corr<0)&(sym!=-halfs[ch]));//sign is flipped if SSE correction was negative, to skew the histogram
								sym^=negmask;
								sym-=negmask;
							}
#endif
							sym=sym<<1^(sym>>31);//pack sign
						}
						else
							sym=upred+aval;//error sign is known
					}
					quantize_pixel(sym, &token, &bypass, &nbits);
					ac_enc(&ec, token, curr_CDF);
					if(nbits)
						ac_enc_bypass(&ec, bypass, nbits);//up to 16 bits
				}
				else
				{
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
					{
						int upred=halfs[ch]-abs(pred), negmask=0;
						if(sym<=(upred<<1))
						{
							error=sym>>1^-(sym&1);
#ifdef ENABLE_BIASCORR
							negmask=-((ibias_corr<0)&(error!=-halfs[ch]));
#endif
						}
						else
						{
							error=sym-upred;
							negmask=-(pred>0);
						}
						error^=negmask;
						error-=negmask;
					}
					yuv[kc]=error+pred;
					curr[kc]=yuv[kc];
				}
				++curr_hist[token];
				++curr_hist[args->tlevels];
				if((kx&(CDFSTRIDE-1))==CDFSTRIDE-1)
					update_CDF(curr_hist, curr_CDF, args->tlevels);
				curr[kc]-=offset;
				if(predidx[kc]==PRED_wgrad)
					wg_update(curr[kc], wg_preds, wg_errors, args->wg_weights+WG_NPREDS*kc);
			}
			if(!args->fwd)
			{
				int kc=0;
				if(image->nch>=3)
				{
					switch(bestrct)
					{
					case RCT_R_G_B:
					case RCT_R_G_BG:
					case RCT_R_G_BR:
					case RCT_R_GR_BR:
					case RCT_R_GR_BG:
						image->data[idx+0]=yuv[0];
						image->data[idx+1]=yuv[1];
						image->data[idx+2]=yuv[2];
						break;
					case RCT_G_B_RG:
					case RCT_G_B_RB:
					case RCT_G_BG_RG:
					case RCT_G_BG_RB:
						image->data[idx+1]=yuv[0];
						image->data[idx+2]=yuv[1];
						image->data[idx+0]=yuv[2];
						break;
					case RCT_B_R_GR:
					case RCT_B_R_GB:
					case RCT_B_RB_GB:
					case RCT_B_RB_GR:
						image->data[idx+2]=yuv[0];
						image->data[idx+0]=yuv[1];
						image->data[idx+1]=yuv[2];
						break;
					case RCT_G_RG_BR:
						image->data[idx+1]=yuv[0];
						image->data[idx+0]=yuv[1];
						image->data[idx+2]=yuv[2];
						break;
					case RCT_B_GB_RG:
						image->data[idx+2]=yuv[0];
						image->data[idx+1]=yuv[1];
						image->data[idx+0]=yuv[2];
						break;
					case RCT_R_BR_GB:
						image->data[idx+0]=yuv[0];
						image->data[idx+2]=yuv[1];
						image->data[idx+1]=yuv[2];
						break;
					case RCT_J2K:
						yuv[0]-=(yuv[1]+yuv[2])>>2;
						yuv[2]+=yuv[0];
						yuv[1]+=yuv[0];
						image->data[idx+1]=yuv[0];
						image->data[idx+2]=yuv[1];
						image->data[idx+0]=yuv[2];
						break;
					case RCT_RCT1:
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						yuv[0]-=yuv[1]>>1;
						yuv[1]+=yuv[0];
						image->data[idx+1]=yuv[0];
						image->data[idx+2]=yuv[1];
						image->data[idx+0]=yuv[2];
						break;
					case RCT_Pei09:
						yuv[0]-=(86*yuv[2]+29*yuv[1]+128)>>8;
						yuv[2]+=yuv[0];
						yuv[1]+=(87*yuv[2]+169*yuv[0]+128)>>8;
						image->data[idx+1]=yuv[0];
						image->data[idx+2]=yuv[1];
						image->data[idx+0]=yuv[2];
						break;
					}
					kc=3;
				}
				if(kc<image->nch)
					image->data[kc]=yuv[kc];
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
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	if(args->fwd)
	{
		ac_enc_flush(&ec);
		if(args->loud)
			printf("Actual %8zd bytes (%+11.2lf)\n\n", args->list.nobj, args->list.nobj-bestsize);
	}
}
int f26_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	int ncores=query_cpu_cores();
	int nblocks=(image->ih+BLOCKSIZE-1)/BLOCKSIZE, nthreads=MINVAR(nblocks, ncores);
	int coffset=sizeof(int)*nblocks;
	ptrdiff_t start=0;
	ptrdiff_t memusage=0;
	ptrdiff_t argssize=nthreads*sizeof(ThreadArgs);
	ThreadArgs *args=(ThreadArgs*)malloc(argssize);
	int depths[OCH_COUNT]={0}, histindices[OCH_COUNT*PRED_COUNT+1]={0}, histsize=0;
	int tlevels=0, clevels=0, statssize=0;
	
	if(fwd)
	{
#ifdef ENABLE_GUIDE
		guide=image;
#endif
		histsize=0;
		for(int kc=0;kc<OCH_COUNT;++kc)
		{
			int depth=image->depth+och_info[kc][2];
			UPDATE_MIN(depth, 16);
			depths[kc]=depth;
			for(int kp=0;kp<PRED_COUNT;++kp)
			{
				int nlevels=1<<depth;
				histindices[kc*PRED_COUNT+kp]=histsize;
				histsize+=nlevels;
			}
		}
		histindices[OCH_COUNT*PRED_COUNT]=histsize;
		start=array_append(data, 0, 1, coffset, 1, 0, 0);
	}
	else//integrity check
	{
		start=coffset;
		for(int kt=0;kt<nblocks;++kt)
		{
			int size=0;
			memcpy(&size, cbuf+sizeof(int)*kt, sizeof(int));
			start+=size;
		}
		if(start!=clen)
			LOG_ERROR("Corrupt file");
		start=coffset;
	}
	{
		int nlevels=2<<image->depth;//chroma-inflated
		int token=0, bypass=0, nbits=0;

		quantize_pixel(nlevels, &token, &bypass, &nbits);
		tlevels=token+1;
		clevels=9;
		statssize=clevels*clevels*(tlevels+1)*image->nch*(int)sizeof(int);
		if(fwd)
		{
			UPDATE_MAX(histsize, statssize);
		}
		else
			histsize=statssize;
	}

	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memusage+=argssize;
	memset(args, 0, argssize);
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=src;
		arg->dst=dst;
		arg->bufsize=sizeof(short[4*OCH_COUNT*1])*(image->iw+16LL);//4 padded rows * OCH_COUNT * {pixels}
		arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));

		arg->histsize=histsize*sizeof(int);
		arg->hist=(int*)malloc(arg->histsize);

		arg->statssize=statssize;
		arg->stats=(unsigned*)malloc(statssize);
		if(!arg->pixels||!arg->hist||!arg->stats)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=arg->bufsize;
		memusage+=arg->histsize;
		memusage+=arg->statssize;
		
		arg->tlevels=tlevels;
		arg->clevels=clevels;
		memcpy(arg->histindices, histindices, sizeof(arg->histindices));
		arg->fwd=fwd;
#ifdef DISABLE_MT
		arg->loud=loud;
#else
		arg->loud=0;
#endif
	}
	for(int kt=0;kt<nblocks;kt+=nthreads)
	{
		int nthreads2=MINVAR(kt+nthreads, nblocks)-kt;
		for(int kt2=0;kt2<nthreads2;++kt2)
		{
			ThreadArgs *arg=args+kt2;
			arg->y1=BLOCKSIZE*(kt+kt2);
			arg->y2=MINVAR(arg->y1+BLOCKSIZE, image->ih);
			if(!fwd)
			{
				int size=0;
				memcpy(&size, cbuf+sizeof(int)*((ptrdiff_t)kt+kt2), sizeof(int));
				arg->decstart=cbuf+start;
				start+=size;
				arg->decend=cbuf+start;
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
				if(loud)
					printf("[%d]  %zd\n", kt+kt2, args[kt2].list.nobj);
				memcpy(data[0]->data+start+sizeof(int)*((ptrdiff_t)kt+kt2), &args[kt2].list.nobj, sizeof(int));
				dlist_appendtoarray(&args[kt2].list, data);
				dlist_clear(&args[kt2].list);
			}
		}
	}
	if(loud)
	{
		ptrdiff_t usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
		t0=time_sec()-t0;
		if(fwd)
		{
			ptrdiff_t csize=data[0]->count-start;
			printf("Size %14td/%14td  %16lf%%  %16lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
			printf("Mem usage: ");
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
		}
		printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		_mm_free(arg->pixels);
		free(arg->hist);
		free(arg->stats);
	}
	free(args);
	return 0;
}