#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT


#include"ac.h"
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif

	#define USE_AC
	#define ENABLE_ZERO	//can't use PRED_ZERO with GR coder
//	#define ENABLE_WG

#define BLOCKSIZE 256
#define NCHPOOL 15
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
typedef enum PredTypeEnum
{
#ifdef ENABLE_ZERO
	PRED_ZERO,
#endif
//	PRED_N,
	PRED_W,
//	PRED_AV2,
//	PRED_AV4,
//	PRED_AV5,
	PRED_CG,
//	PRED_CGU,
	PRED_WP,
#ifdef ENABLE_WG
	PRED_WG,
#endif

	PRED_COUNT,
} PredType;
typedef struct _ThreadArgs
{
	const Image *src;
	Image *dst;
	int fwd, loud, y1, y2;
	int bufsize, histsize;
	int *pixels, *hist;

	DList list;
	const unsigned char *decstart, *decend;

#ifdef USE_AC
	int tlevels, clevels, statssize;
	unsigned *stats;
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
static const int wp_param_idx[NCHPOOL]=
{
	2, 0, 1,	//RCT_NONE
	2, 1,		//RCT_SubG
	0, 1, 2,	//RCT_J2KG
	0, 2,		//RCT_J2KB
	0,		//RCT_J2KR
	0, 1,		//RCT_1
	0, 1,		//RCT_Pei09
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
#define WP_NPARAMS 11
static const short wp_params[WP_NPARAMS*3]=//signed fixed 7.8 bit
{
	0x0DB8,  0x0E22,  0x181F,  0x0BF3, -0x005C, -0x005B,  0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,//Y
	0x064C,  0x0F31,  0x1040,  0x0BF8, -0x0007, -0x000D, -0x0085, -0x0063, -0x00A2, -0x0017,  0x00F2,//Cb
	0x0B37,  0x110B,  0x121B,  0x0BFC, -0x0001,  0x000E, -0x0188, -0x00E7, -0x00BB, -0x004A,  0x00BA,//Cr
};
#ifdef USE_AC

#define CDFSTRIDE 64

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
#define CONFIG_EXP 5
#define CONFIG_MSB 2
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
static int quantize_ctx(int val)
{
	int negmask=-(val<0);
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
}
#define QUANTIZE_CTX(X) quantize_ctx(X)+chalf
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
	int stride=tlevels+1, chsize=clevels*clevels*stride;
	for(int kc=0;kc<nch;++kc)
	{
		int sum=0;
		int *curr_hist=hist+chsize*kc;
		unsigned *curr_CDF=stats+chsize*kc;

		//memset(curr_hist, 0, chsize);
		//memset(curr_CDF, 0, chsize);
		for(int k=0;k<tlevels;++k)
		{
			int val=(tlevels-k)>>1;
			val+=!val;
			val*=val;
			if(sum+val*val*val<sum)
				LOG_ERROR("");
			sum+=curr_hist[k]=val*val*val;
		}
		curr_hist[tlevels]=sum;
		update_CDF(curr_hist, curr_CDF, tlevels);
		memfill(curr_hist+stride, curr_hist, ((size_t)chsize-stride)*sizeof(int), stride*sizeof(int));
		memfill(curr_CDF+stride, curr_CDF, ((size_t)chsize-stride)*sizeof(int), stride*sizeof(int));
	}
}
#endif
static int wp_predict(const short *params, int NN, int NW, int N, int NE, int W, const int *eNW, const int *eN, const int *eNE, const int *eW, int *preds)
{
	int pred;
	long long lpred=0, wsum=0;
	preds[0]=(W+NE-N)<<8;
	preds[1]=(N<<8)-((eN[0]+eW[0]+eNE[0])*params[4]>>8);
	preds[2]=(W<<8)-((eN[0]+eW[0]+eNW[0])*params[5]>>8);
	preds[3]=(N<<8)-(((eNW[0]*params[6]+eN[0]*params[7]+eNE[0]*params[8])>>8)+(NN-N)*params[9]+(NW-W)*params[10]);
	for(int k=0;k<4;++k)
	{
		int weight=(params[k]<<8)/(eNW[k+1]+eN[k+1]+eNE[k+1]+1);
		lpred+=(long long)weight*preds[k];
		wsum+=weight;
	}
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
static int wg_predict(int NNN, int NN, int NNE, int NW, int N, int NE, int NEEE, int WWW, int WW, int W, const int *eW, int *preds)
{
	int pred;
	long long lpred=0, wsum=0;

	preds[0]=N+W-NW;
	preds[1]=W+NE-N;
	preds[2]=N+NE-NNE;
	preds[3]=N;
	preds[4]=W;
	preds[5]=3*(N-NN)+NNN;
	preds[6]=3*(W-WW)+WWW;
	preds[7]=(W+NEEE)>>1;
	
	lpred=0;
	wsum=0;
	for(int k=0;k<8;++k)
	{
		int weight=0x1000000/(eW[k]+1);
		lpred+=(long long)weight*preds[k];
		wsum+=weight;
	}
	lpred/=wsum+1;
	pred=(int)lpred;
	CLAMP3_32(pred, (int)lpred, N, W, NE);
	return pred;
}
static void wg_update(int curr, int pred, const int *preds, int *eW)
{
	for(int k=0;k<8;++k)
	{
		int err=abs(curr-preds[k]);
		eW[k]+=err;
		eW[k]-=(eW[k]+3)>>2;
	}
}
static void block_enc(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
#ifdef USE_AC
	ArithmeticCoder ec;
	int chalf=args->clevels>>1, nctx=args->clevels*args->clevels, cdfstride=args->tlevels+1;
	int token=0, bypass=0, nbits=0;
#else
	GolombRiceCoder ec;
#endif
	Image const *image=args->src;
	int nlevels[NCHPOOL]={0}, halfs[NCHPOOL]={0};
	double csizes[NCHPOOL*PRED_COUNT]={0};
	const short *ch_params[NCHPOOL]={0};
	int res=image->iw*(args->y2-args->y1);

	for(int k=0;k<NCHPOOL;++k)
	{
		nlevels[k]=1<<(image->depth+depth_inf[k]);
		halfs[k]=nlevels[k]>>1;
		ch_params[k]=wp_params+WP_NPARAMS*wp_param_idx[k];
	}
	memset(args->pixels, 0, args->bufsize);
	memset(args->hist, 0, args->histsize);
	for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
	{
		ALIGN(16) int *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*6*NCHPOOL,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*6*NCHPOOL,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*6*NCHPOOL,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*6*NCHPOOL,
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
		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			int
				*NNN	=rows[3]+0*6*NCHPOOL,
				*NN	=rows[2]+0*6*NCHPOOL,
				*NNE	=rows[2]+1*6*NCHPOOL,
				*NW	=rows[1]-1*6*NCHPOOL,
				*N	=rows[1]+0*6*NCHPOOL,
				*NE	=rows[1]+1*6*NCHPOOL,
				*NEEE	=rows[1]+3*6*NCHPOOL,
				*WWW	=rows[0]-3*6*NCHPOOL,
				*WW	=rows[0]-2*6*NCHPOOL,
				*W	=rows[0]-1*6*NCHPOOL,
				*curr	=rows[0]+0*6*NCHPOOL;

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
			//
			//13	y5-P(y5prev)		Pei09		b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
			//14	u5-P(u5prev)

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
				int kc2=kc*6;
				int cgrad;
				int wpred, wp_preds[4];
#ifdef ENABLE_WG
				int wg_preds[8];
#endif
				int offset=0, val;
				//int
				//	vx=abs(W[kc2+1]-WW[kc2+1])+abs(N[kc2+1]-NW[kc2+1])+abs(NE[kc2+1]-N[kc2+1])+1,
				//	vy=abs(W[kc2+1]-NW[kc2+1])+abs(N[kc2+1]-NN[kc2+1])+abs(NE[kc2+1]-NNE[kc2+1])+1;
				//int update=(vx<vy?N[kc2+1]:W[kc2+1])/4;
				
				wpred=wp_predict(ch_params[kc], NN[kc2], NW[kc2], N[kc2], NE[kc2], W[kc2], NW+kc2+1, N+kc2+1, NE+kc2+1, W+kc2+1, wp_preds);
				MEDIAN3_32(cgrad, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
				
#ifdef ENABLE_ZERO
				preds[PRED_ZERO]=0;
#endif
			//	preds[PRED_N]=N[kc2];
				preds[PRED_W]=W[kc2];
			//	preds[PRED_AV2]=(N[kc2]+W[kc2])/2;
			//	preds[PRED_AV4]=(4*(N[kc2]+W[kc2])+NE[kc2]-NW[kc2])/8;
			//	preds[PRED_AV5]=W[kc2]+(5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])/8;
				preds[PRED_CG]=cgrad;
			//	CLAMP2_32(cgrad, cgrad+update, -halfs[kc], halfs[kc]-1);
			//	preds[PRED_CGU]=cgrad;
				preds[PRED_WP]=(wpred+127)>>8;
#ifdef ENABLE_WG
				preds[PRED_WG]=wg_predict(NNN[kc2], NN[kc2], NNE[kc2], NW[kc2], N[kc2], NE[kc2], NEEE[kc2], WWW[kc2], WW[kc2], W[kc2], wg_errors, wg_preds);
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
					++args->hist[(kc*PRED_COUNT+k)<<(image->depth+1)|val];
				}
				curr[kc2+0]=comp[kc]-offset;
				wp_update(comp[kc], wpred, wp_preds, curr+kc2+1, NE+kc2+1);
#ifdef ENABLE_WG
				wg_update(comp[kc], preds[PRED_WG], wg_preds, wg_errors);
#endif
			}
			rows[0]+=6*NCHPOOL;
			rows[1]+=6*NCHPOOL;
			rows[2]+=6*NCHPOOL;
			rows[3]+=6*NCHPOOL;
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
	char predsel[NCHPOOL]={0};
	double bestsize=0;
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
	const int *group=combinations;
	bestsize=
		csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
		csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
		csizes[group[2]*PRED_COUNT+predsel[group[2]]];
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
	if(args->loud)
	{
		static const char *rctnames[RCT_COUNT]=
		{
			"NONE",
			"SubG",
			"JPEG2000-G",
			"JPEG2000-B",
			"JPEG2000-R",
			"RCT1",
			"Pei09",
		};
		static const char *prednames[PRED_COUNT]=
		{
#ifdef ENABLE_ZERO
			"Zero",
#endif
		//	"N",
			"W",
		//	"Av2",
		//	"Av4",
		//	"Av5",
			"CG",
		//	"CGU",
			"WP",
#ifdef ENABLE_WG
			"WG",
#endif
		};

		double defsize=
			csizes[combinations[3*2+0]*PRED_COUNT+PRED_WP]+
			csizes[combinations[3*2+1]*PRED_COUNT+PRED_WP]+
			csizes[combinations[3*2+2]*PRED_COUNT+PRED_WP];
		printf("Y %5d~%5d  default %12.2lf bytes  current %12.2lf bytes (%+12.2lf)  %s [YUV: %s %s %s]\n",
			args->y1, args->y2,
			defsize,
			bestsize,
			bestsize-defsize,
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
				printf("  %11s", prednames[k2]);
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
					printf(" %11.2lf%c", csizes[k*PRED_COUNT+k2], k2==best?'*':' ');
				printf("  %s\n", poolnames[k]);
			}
		}
	}
	dlist_init(&args->list, 1, 1024, 0);
	dlist_push_back(&args->list, &flag, sizeof(char[2]));
#ifdef USE_AC
	ac_enc_init(&ec, &args->list);
	init_stats(args->hist, args->stats, args->tlevels, args->clevels, image->nch);
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
		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			int
				*NNN	=rows[3]+0*7*4,
				*NN	=rows[2]+0*7*4,
				*NNE	=rows[2]+1*7*4,
				*NW	=rows[1]-1*7*4,
				*N	=rows[1]+0*7*4,
				*NE	=rows[1]+1*7*4,
				*NEEE	=rows[1]+3*7*4,
				*WWW	=rows[0]-3*7*4,
				*WW	=rows[0]-2*7*4,
				*W	=rows[0]-1*7*4,
				*curr	=rows[0]+0*7*4;

			//if(ky==330&&kx==146)//
			//if(ky==28&&kx==624)//
			//if(idx==26147568)//
			//if(idx==2304)//
			//if(ky==0&&kx==64)//
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
				case PRED_W:
					pred=W[kc2];
					break;
				//case PRED_AV2:
				//	pred=(N[kc2]+W[kc2])/2;
				//	break;
				//case PRED_AV4:
				//	pred=(4*(N[kc2]+W[kc2])+NE[kc2]-NW[kc2])/8;
				//	break;
				//case PRED_AV5:
				//	pred=W[kc2]+(5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])/8;
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
				case PRED_WP:
					wp_result=wp_predict(ch_params[ch], NN[kc2], NW[kc2], N[kc2], NE[kc2], W[kc2], NW+kc2+2, N+kc2+2, NE+kc2+2, W+kc2+2, wp_preds);
					pred=(wp_result+127)>>8;
					break;
#ifdef ENABLE_WG
				case PRED_WG:
					pred=wg_predict(NNN[kc2], NN[kc2], NNE[kc2], NW[kc2], N[kc2], NE[kc2], NEEE[kc2], WWW[kc2], WW[kc2], W[kc2], wg_errors, wp_preds);
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
						sym=sym<<1^-(sym<0);//pack sign
					}
					else
						sym=upred+aval;//error sign is known
				}
#ifdef USE_AC
				int qeN=QUANTIZE_CTX(N[kc2+1]);
				int qeW=QUANTIZE_CTX(W[kc2+1]);
				int cidx=cdfstride*(nctx*kc+args->clevels*qeN+qeW);
				int *curr_hist=args->hist+cidx;
				unsigned *curr_CDF=args->stats+cidx;
				quantize_pixel(sym, &token, &bypass, &nbits);
#ifdef _DEBUG
				if((unsigned)token>=(unsigned)args->tlevels)
					LOG_ERROR("Token OOB %d/%d", token, args->tlevels);
#endif
				ac_enc(&ec, token, curr_CDF);
				if(nbits)
					ac_enc_bypass(&ec, bypass, nbits);//up to 16 bits

				if(curr_hist[args->tlevels]+1>0x1000)
				{
					int sum=0;
					//update_CDF(curr_hist, curr_CDF, args->tlevels);
					for(int k=0;k<args->tlevels;++k)
						sum+=curr_hist[k]=(curr_hist[k]+1)>>1;
					curr_hist[args->tlevels]=sum;
				}
				++curr_hist[token];
				++curr_hist[args->tlevels];
				if((kx&(CDFSTRIDE-1))==CDFSTRIDE-1)		//FIXME this updates just the current CDF at regular intervals
					update_CDF(curr_hist, curr_CDF, args->tlevels);

				curr[kc2+1]=val;
#else
				gr_enc_POT(&ec, sym, FLOOR_LOG2(W[kc2+1]+1));
				curr[kc2+1]=(2*W[kc2+1]+sym+NEEE[kc2+1])>>2;
#endif
				curr[kc2+0]=yuv[kc]-offset;
				if(predidx[kc]==PRED_WP)
					wp_update(curr[kc2+0], wp_result, wp_preds, curr+kc2+2, NE+kc2+2);
#ifdef ENABLE_WG
				else if(predidx[kc]==PRED_WG)
					wg_update(curr[kc2+0], pred, wp_preds, wg_errors);
#endif
			}
			rows[0]+=7*4;
			rows[1]+=7*4;
			rows[2]+=7*4;
			rows[3]+=7*4;
		}
	}
#ifdef USE_AC
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
	int chalf=args->clevels>>1, nctx=args->clevels*args->clevels, cdfstride=args->tlevels+1;
	int token=0, bypass=0, nbits=0;
#else
	GolombRiceCoder ec;
#endif
	Image *image=args->dst;
	const unsigned char *srcstart=args->decstart, *srcend=args->decend;
	int flag=0, bestrct;
	int combination[3]={0}, predidx[3]={0};
	int halfs[NCHPOOL]={0};
	const short *ch_params[NCHPOOL]={0};
	
	for(int k=0;k<NCHPOOL;++k)
	{
		halfs[k]=1<<(image->depth+depth_inf[k])>>1;
		ch_params[k]=wp_params+WP_NPARAMS*wp_param_idx[k];
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
		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			int
				*NNN	=rows[3]+0*7*4,
				*NN	=rows[2]+0*7*4,
				*NNE	=rows[2]+1*7*4,
				*NW	=rows[1]-1*7*4,
				*N	=rows[1]+0*7*4,
				*NE	=rows[1]+1*7*4,
				*NEEE	=rows[1]+3*7*4,
				*WWW	=rows[0]-3*7*4,
				*WW	=rows[0]-2*7*4,
				*W	=rows[0]-1*7*4,
				*curr	=rows[0]+0*7*4;

			//if(ky==330&&kx==146)//
			//if(ky==28&&kx==624)//
			//if(idx==26147568)//
			//if(idx==2304)//
			//if(ky==0&&kx==64)//
			//	printf("");

			for(int kc=0;kc<3;++kc)
			{
				int kc2=kc*7, pred=0, offset=0, ch=combination[kc], val, sym;
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
				case PRED_W:
					pred=W[kc2];
					break;
				//case PRED_AV2:
				//	pred=(N[kc2]+W[kc2])/2;
				//	break;
				//case PRED_AV4:
				//	pred=(4*(N[kc2]+W[kc2])+NE[kc2]-NW[kc2])/8;
				//	break;
				//case PRED_AV5:
				//	pred=W[kc2]+(5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])/8;
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
				case PRED_WP:
					wp_result=wp_predict(ch_params[ch], NN[kc2], NW[kc2], N[kc2], NE[kc2], W[kc2], NW+kc2+2, N+kc2+2, NE+kc2+2, W+kc2+2, wp_preds);
					pred=(wp_result+127)>>8;
					break;
#ifdef ENABLE_WG
				case PRED_WG:
					pred=wg_predict(NNN[kc2], NN[kc2], NNE[kc2], NW[kc2], N[kc2], NE[kc2], NEEE[kc2], WWW[kc2], WW[kc2], W[kc2], wg_errors, wp_preds);
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
				int qeN=QUANTIZE_CTX(N[kc2+1]);
				int qeW=QUANTIZE_CTX(W[kc2+1]);
				int cidx=cdfstride*(nctx*kc+args->clevels*qeN+qeW);
				int *curr_hist=args->hist+cidx;
				unsigned *curr_CDF=args->stats+cidx;

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
				if(curr_hist[args->tlevels]+1>0x1000)
				{
					int sum=0;
					//update_CDF(curr_hist, curr_CDF, args->tlevels);
					for(int k=0;k<args->tlevels;++k)
						sum+=curr_hist[k]=(curr_hist[k]+1)>>1;
					curr_hist[args->tlevels]=sum;
				}
				++curr_hist[token];
				++curr_hist[args->tlevels];
				if((kx&(CDFSTRIDE-1))==CDFSTRIDE-1)
					update_CDF(curr_hist, curr_CDF, args->tlevels);

				curr[kc2+1]=val;
#else
				curr[kc2+1]=(2*W[kc2+1]+sym+NEEE[kc2+1])>>2;
#endif
				val+=pred;
				yuv[kc]=val;
				curr[kc2+0]=yuv[kc]-offset;
				if(predidx[kc]==PRED_WP)
					wp_update(curr[kc2+0], wp_result, wp_preds, curr+kc2+2, NE+kc2+2);
#ifdef ENABLE_WG
				else if(predidx[kc]==PRED_WG)
					wg_update(curr[kc2+0], pred, wp_preds, wg_errors);
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
		clevels=quantize_ctx(nlevels>>1)<<1|1;
		statssize=clevels*clevels*(tlevels+1)*image->nch*(int)sizeof(int);
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
		arg->bufsize=sizeof(int[4*NCHPOOL*6])*(image->iw+16LL);//4 padded rows * NCHPOOL * {pixels, hires-errors, 4 pred hires-errors}
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
		ptrdiff_t dststart=array_append(data, 0, 1, sizeof(int)*nblocks, 1, 0, 0);
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
				if(loud)
					printf("[%3d]  %8zd\n", kt+kt2, args[kt2].list.nobj);
				memcpy(data[0]->data+dststart+sizeof(int)*((ptrdiff_t)kt+kt2), &args[kt2].list.nobj, sizeof(int));
				dlist_appendtoarray(&args[kt2].list, data);
				dlist_clear(&args[kt2].list);
			}
		}
		if(loud)
		{
			ptrdiff_t
				csize=data[0]->count-dststart,
				usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
			t0=time_sec()-t0;
			printf("Size %14td/%14td  %16lf%%  %16lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
			printf("E %16.6lf sec  %16.6lf MB/s  Mem usage: ", t0, usize/(t0*1024*1024));
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
		}
	}
	else
	{
		const unsigned char *dstptr=cbuf+sizeof(int)*nblocks;
		int dec_offset=0;

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