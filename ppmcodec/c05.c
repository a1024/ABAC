#include"codec.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT


//	#define ENABLE_SSE


#define CODECNAME "C05"
#include"entropy.h"

#define BLOCKSIZE 512
#define MAXPRINTEDBLOCKS 50

#define OCHLIST\
	OCH(R)\
	OCH(G)\
	OCH(B)\
	OCH(RG)\
	OCH(RB)\
	OCH(GB)\
	OCH(GR)\
	OCH(BR)\
	OCH(BG)\
	OCH(R2)\
	OCH(G2)\
	OCH(B2)
typedef enum _OCHIndex
{
#define OCH(LABEL) OCH_##LABEL,
	OCHLIST
#undef  OCH
	OCH_COUNT,
} OCHIndex;
static const char *och_names[OCH_COUNT]=
{
#define OCH(LABEL) #LABEL,
	OCHLIST
#undef  OCH
};

#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	3,	3, 3, 0)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	3,	1, 3, 0)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	3,	0, 3, 0)\
	RCT(R_G_B2,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	3,	0, 1, 1)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	0,	0, 3, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	0,	1, 3, 0)\
	RCT(R_GR_B2,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	0,	0, 1, 1)\
	RCT(R_B_G2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	3,	0, 1, 1)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	0,	1, 3, 0)\
	RCT(R_BR_G2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	0,	0, 1, 1)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	3,	0, 3, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	3,	1, 3, 0)\
	RCT(G_B_R2,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	3,	0, 1, 1)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	0,	0, 3, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	0,	1, 3, 0)\
	RCT(G_BG_R2,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	0,	0, 1, 1)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	0,	1, 3, 0)\
	RCT(G_RG_B2,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	0,	0, 1, 1)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	3,	1, 3, 0)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	3,	0, 3, 0)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	0,	0, 3, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	0,	1, 3, 0)\
	RCT(B_RB_G2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	0,	0, 1, 1)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	0,	1, 3, 0)\
	RCT(B_GB_R2,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	0,	0, 1, 1)
typedef enum _RCTIndex
{
#define RCT(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][10]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2)\
	{YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX,  YPERM, UPERM, VPERM,  UOFF1,  VOFF1, VOFF2, VSH2) #LABEL,
	RCTLIST
#undef  RCT
};

#define PREDLIST\
	PRED(CG)
typedef enum _PredIndex
{
#define PRED(LABEL) PRED_##LABEL,
	PREDLIST
#undef  PRED
	PRED_COUNT,
} PredIndex;
static const char *pred_names[PRED_COUNT]=
{
#define PRED(LABEL) #LABEL,
	PREDLIST
#undef  PRED
};
typedef enum _NBIndex
{
	NB_NW,		NB_N,
	NB_W,		NB_curr,
} NBIndex;


typedef struct _ThreadArgs
{
	const unsigned char *src;
	unsigned char *dst;
	int iw, ih;

	int fwd, test, loud, x1, x2, y1, y2;
	int bufsize;
	short pixels[(BLOCKSIZE+16)*4*4*2];//4 padded rows * 4 channels max * {pixel, abs(e)}

	BList list;
	const unsigned char *decstart, *decend;
	
	int hist[OCH_COUNT*PRED_COUNT<<8];

	//aux
	int blockidx;
	double bestsize;
	int bestrct, predidx[3];
	double t_analysis;
} ThreadArgs;
static void block_thread(void *param)
{
	const int half=128, nch=3;
	ThreadArgs *args=(ThreadArgs*)param;
	GolombRiceCoder ec;
	const unsigned char *image=args->fwd?args->src:args->dst;
	unsigned char bestrct=0, combination[6]={0}, predidx[4]={0};
	int ystride=args->iw*3;
	
	if(args->fwd)
	{
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		unsigned char predsel[OCH_COUNT]={0};
		int res=(args->x2-args->x1-1)/5*5*(args->y2-args->y1-1);
		
		args->t_analysis=time_sec();
		memset(args->hist, 0, sizeof(args->hist));
		for(int ky=args->y1+1;ky<args->y2;++ky)//analysis loop
		{
			int kx=args->x1+1;
			const unsigned char *ptr=image+3*(args->iw*ky+kx);

			__m256i amin=_mm256_set1_epi16(-128);
			__m256i amax=_mm256_set1_epi16(127);
			__m128i half8=_mm_set1_epi8(128);
			__m128i shuf=_mm_set_epi8(
				-1,
				12, 14, 13,
				 9, 11, 10,
				 6,  8,  7,
				 3,  5,  4,
				 0,  2,  1
				//15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
			);
			ALIGN(32) short result[16]={0};
			for(;kx<args->x2-4;kx+=5, ptr+=15)
			{
				__m256i
					nb0[4],//rgb
					nb1[4],//gbr
					nb2[4],//rgb - gbr
					nb3[4],//gbr - rgb
					nb4[4],//(gbr+brg)/2
					nb5[4];//rgb - (gbr+brg)/2
				__m256i vmin[4], vmax[4], pred;
				{
					__m128i nb8[4]=//8-bit
					{
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride-1*3+0)), half8),//NW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+0*3+0)), half8),//N
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride-1*3+0)), half8),//W
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride+0*3+0)), half8),//curr
					};
					for(int k=0;k<(int)_countof(nb8);++k)
					{
						__m128i temp;
						__m256i t2;
						nb0[k]=_mm256_cvtepi8_epi16(nb8[k]);
						temp=_mm_shuffle_epi8(nb8[k], shuf);
						nb1[k]=_mm256_cvtepi8_epi16(temp);
						t2=_mm256_cvtepi8_epi16(_mm_shuffle_epi8(temp, shuf));
						t2=_mm256_add_epi16(t2, nb1[k]);
						t2=_mm256_srai_epi16(t2, 1);
						nb2[k]=_mm256_sub_epi16(nb0[k], nb1[k]);
						nb3[k]=_mm256_sub_epi16(nb1[k], nb0[k]);
						nb4[k]=t2;
						nb5[k]=_mm256_sub_epi16(nb0[k], t2);
					}
				}
#define UPDATE(PREDIDX, IDX0, IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7, IDX8, IDX9, IDXA, IDXB, IDXC, IDXD, IDXE)\
	do\
	{\
		pred=_mm256_slli_epi16(pred, 8);\
		pred=_mm256_srai_epi16(pred, 8);\
		pred=_mm256_sub_epi16(pred, amin);\
		_mm256_store_si256((__m256i*)result, pred);\
		++args->hist[(IDX0*PRED_COUNT+PREDIDX)<<8|result[0x0]];\
		++args->hist[(IDX1*PRED_COUNT+PREDIDX)<<8|result[0x1]];\
		++args->hist[(IDX2*PRED_COUNT+PREDIDX)<<8|result[0x2]];\
		++args->hist[(IDX3*PRED_COUNT+PREDIDX)<<8|result[0x3]];\
		++args->hist[(IDX4*PRED_COUNT+PREDIDX)<<8|result[0x4]];\
		++args->hist[(IDX5*PRED_COUNT+PREDIDX)<<8|result[0x5]];\
		++args->hist[(IDX6*PRED_COUNT+PREDIDX)<<8|result[0x6]];\
		++args->hist[(IDX7*PRED_COUNT+PREDIDX)<<8|result[0x7]];\
		++args->hist[(IDX8*PRED_COUNT+PREDIDX)<<8|result[0x8]];\
		++args->hist[(IDX9*PRED_COUNT+PREDIDX)<<8|result[0x9]];\
		++args->hist[(IDXA*PRED_COUNT+PREDIDX)<<8|result[0xA]];\
		++args->hist[(IDXB*PRED_COUNT+PREDIDX)<<8|result[0xB]];\
		++args->hist[(IDXC*PRED_COUNT+PREDIDX)<<8|result[0xC]];\
		++args->hist[(IDXD*PRED_COUNT+PREDIDX)<<8|result[0xD]];\
		++args->hist[(IDXE*PRED_COUNT+PREDIDX)<<8|result[0xE]];\
	}while(0)

				//CG
				vmin[0]=_mm256_min_epi16(nb0[NB_N], nb0[NB_W]);
				vmax[0]=_mm256_max_epi16(nb0[NB_N], nb0[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb0[NB_N], nb0[NB_W]), nb0[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				vmin[1]=_mm256_min_epi16(nb2[NB_N], nb2[NB_W]);
				vmax[1]=_mm256_max_epi16(nb2[NB_N], nb2[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb2[NB_N], nb2[NB_W]), nb2[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				vmin[2]=_mm256_min_epi16(nb3[NB_N], nb3[NB_W]);
				vmax[2]=_mm256_max_epi16(nb3[NB_N], nb3[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb3[NB_N], nb3[NB_W]), nb3[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				vmin[3]=_mm256_min_epi16(nb5[NB_N], nb5[NB_W]);
				vmax[3]=_mm256_max_epi16(nb5[NB_N], nb5[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb5[NB_N], nb5[NB_W]), nb5[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);
			}
		}
		for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)
		{
			int *curr_hist=args->hist+((size_t)kc<<8);
			for(int ks=0;ks<256;++ks)
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
		for(int kt=0;kt<RCT_COUNT;++kt)//select best RCT
		{
			const unsigned char *group=rct_combinations[kt];
			double csize=
				csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
				csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
				csizes[group[2]*PRED_COUNT+predsel[group[2]]];
			if(!kt||bestsize>csize)
				bestsize=csize, bestrct=kt;
		}
		memcpy(combination, rct_combinations[bestrct], sizeof(combination));
		predidx[0]=predsel[combination[0]];
		predidx[1]=predsel[combination[1]];
		predidx[2]=predsel[combination[2]];
		args->bestsize=bestsize;
		args->bestrct=bestrct;
		args->predidx[0]=predidx[0];
		args->predidx[1]=predidx[1];
		args->predidx[2]=predidx[2];
		args->t_analysis=time_sec()-args->t_analysis;
		if(args->loud)
		{
			printf("Y %5d~%5d  best %12.2lf bytes  %s\n",
				args->y1, args->y2,
				bestsize,
				rct_names[bestrct]
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
				const unsigned char *group=rct_combinations[kt];
				double csize=
					csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
					csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
					csizes[group[2]*PRED_COUNT+predsel[group[2]]];
				printf("%12.2lf %c  %-10s\n",
					csize,
					kt==bestrct?'*':' ',
					rct_names[kt]
				);
			}
		}
		blist_init(&args->list);
		gr_enc_init(&ec, &args->list);
		gr_enc_NPOT(&ec, bestrct, RCT_COUNT);
	}
	else
	{
		gr_dec_init(&ec, args->decstart, args->decend);
		bestrct=gr_dec_NPOT(&ec, RCT_COUNT);
	}
#ifdef ENABLE_SSE
	int sse[3][128][2]={0};
#endif
	memset(args->pixels, 0, sizeof(args->pixels));
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) short *rows[]=
		{
			args->pixels+((BLOCKSIZE+16LL)*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKSIZE+16LL)*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKSIZE+16LL)*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKSIZE+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		int yuv[4]={0};
		int pred=0, error=0;
		const unsigned char *combination=rct_combinations[bestrct];
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int idx=nch*(args->iw*ky+kx);
			short
			//	*NNN	=rows[3]+0*4*2,
				*NNNEE	=rows[3]+2*4*2,
			//	*NNWW	=rows[2]-2*4*2,
			//	*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NNEE	=rows[2]+2*4*2,
			//	*NNEEE	=rows[2]+3*4*2,
			//	*NWW	=rows[1]-2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*NEE	=rows[1]+2*4*2,
				*NEEE	=rows[1]+3*4*2,
				*NEEEE	=rows[1]+4*4*2,
			//	*WWWW	=rows[0]-4*4*2,
				*WWW	=rows[0]-3*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
			(void)NNNEE;
			(void)NN;
			(void)NNE;
			(void)NNEE;
			(void)NEEEE;
			if(ky==args->y1)
				NEEEE=NEEE=NEE=NE=NW=N=W;
			else if(kx==args->x1)
				NW=WWW=WW=W=N;
			else if(kx>args->x2-4)
				NEEE-=(kx-(args->x2-4))*4*2;
			int nbypass[]=
			{
				FLOOR_LOG2(W[4+0]+1),
				FLOOR_LOG2(W[4+1]+1),
				FLOOR_LOG2(W[4+2]+1),
			};
			ALIGN(16) short preds[8];
			__m128i mNW	=_mm_loadu_si128((__m128i*)NW);
			__m128i mN	=_mm_loadu_si128((__m128i*)N);
			__m128i mW	=_mm_loadu_si128((__m128i*)W);
			__m128i mmin=_mm_min_epi16(mN, mW);
			__m128i mmax=_mm_max_epi16(mN, mW);
			__m128i mp=_mm_add_epi16(mN, mW);
			mp=_mm_sub_epi16(mp, mNW);
			mp=_mm_max_epi16(mp, mmin);
			mp=_mm_min_epi16(mp, mmax);
			_mm_store_si128((__m128i*)preds, mp);
#ifdef ENABLE_SSE
			int *curr_sse[]=
			{
				sse[0][(preds[0]+half*2)>>2],
				sse[1][(preds[1]+half*2)>>2],
				sse[2][(preds[2]+half*2)>>2],
			};
			//if((unsigned)(preds[0]+half*2)>511||(unsigned)(preds[1]+half*2)>511||(unsigned)(preds[2]+half*2)>511)
			//	LOG_ERROR("");
			preds[0]+=curr_sse[0][1]/(curr_sse[0][0]+5);
			preds[1]+=curr_sse[1][1]/(curr_sse[1][0]+5);
			preds[2]+=curr_sse[2][1]/(curr_sse[2][0]+5);
#endif
			int sym=0, offset;
//	#define UPDATE_FORMULA(IDX) (MAXVAR(N[IDX], W[IDX])+2*(sym+NE[IDX])+NEE[IDX]+NEEE[IDX])/7	//Formula10
//	#define UPDATE_FORMULA(IDX) (MAXVAR(N[IDX], W[IDX])+sym+NE[IDX]+NEEE[IDX])>>2		//Formula11
//	#define UPDATE_FORMULA(IDX) (NE[IDX]>sym?sym+2*NE[IDX]+NEE[IDX]:W[IDX]+2*sym+NE[IDX])>>2//X
	#define UPDATE_FORMULA(IDX) (W[IDX]+sym+NE[IDX]+NEEE[IDX])>>2	//Formula8
//	#define UPDATE_FORMULA(IDX) (sym+NE[IDX]+NEE[IDX]+NEEE[IDX])>>2
//	#define UPDATE_FORMULA(IDX) (2*W[IDX]+sym+NEEE[IDX])>>2		//Formula1
//	#define UPDATE_FORMULA(IDX) (W[IDX]+sym+NEEE[IDX])/3
//	#define UPDATE_FORMULA(IDX) (WWW[IDX]+5*WW[IDX]+10*(W[IDX]+sym)+NE[IDX])/27	//Formula5	X
			if(args->fwd)
			{
				yuv[0]=args->src[idx+combination[3+0]]-128;
				yuv[1]=args->src[idx+combination[3+1]]-128;
				yuv[2]=args->src[idx+combination[3+2]]-128;

				//enc Y
				pred=preds[0];
				error=yuv[0]-pred;
				{
					int upred=half-abs(pred), aval=abs(error);
					if(aval<=upred)
					{
						sym=error;
#ifdef ENABLE_BIASCORR
						{
							int negmask=-((ibias_corr<0)&(sym!=-half));//sign is flipped if SSE correction was negative, to skew the histogram
							sym^=negmask;
							sym-=negmask;
						}
#endif
						sym=sym<<1^sym>>31;//pack sign
					}
					else
						sym=upred+aval;//error sign is known
				}
				gr_enc(&ec, sym, nbypass[0]);
				curr[0]=yuv[0];
				curr[4]=UPDATE_FORMULA(4);
				
				//enc U
				offset=yuv[combination[6]];
				pred=preds[1]+offset;
				CLAMP2_32(pred, pred, -half, half-1);
				error=yuv[1]-pred;
				{
					int upred=half-abs(pred), aval=abs(error);
					if(aval<=upred)
					{
						sym=error;
#ifdef ENABLE_BIASCORR
						{
							int negmask=-((ibias_corr<0)&(sym!=-half));//sign is flipped if SSE correction was negative, to skew the histogram
							sym^=negmask;
							sym-=negmask;
						}
#endif
						sym=sym<<1^sym>>31;//pack sign
					}
					else
						sym=upred+aval;//error sign is known
				}
				gr_enc(&ec, sym, nbypass[1]);
				curr[1]=yuv[1]-offset;
				curr[5]=UPDATE_FORMULA(5);

				//enc V
				offset=(yuv[combination[7]]+yuv[combination[8]])>>combination[9];
				pred=preds[2]+offset;
				CLAMP2_32(pred, pred, -half, half-1);
				error=yuv[2]-pred;
				{
					int upred=half-abs(pred), aval=abs(error);
					if(aval<=upred)
					{
						sym=error;
#ifdef ENABLE_BIASCORR
						{
							int negmask=-((ibias_corr<0)&(sym!=-half));//sign is flipped if SSE correction was negative, to skew the histogram
							sym^=negmask;
							sym-=negmask;
						}
#endif
						sym=sym<<1^sym>>31;//pack sign
					}
					else
						sym=upred+aval;//error sign is known
				}
				gr_enc(&ec, sym, nbypass[2]);
				curr[2]=yuv[2]-offset;
				curr[6]=UPDATE_FORMULA(6);
			}
			else
			{
				//dec Y
				pred=preds[0];
				sym=gr_dec(&ec, nbypass[0]);
				{
					int upred=half-abs(pred), negmask=0;
					if(sym<=(upred<<1))
					{
						error=sym>>1^-(sym&1);
#ifdef ENABLE_BIASCORR
						negmask=-((ibias_corr<0)&(error!=-half));
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
				curr[0]=yuv[0]=error+pred;
				curr[4]=UPDATE_FORMULA(4);

				//decU
				offset=yuv[combination[6]];
				pred=preds[1]+offset;
				CLAMP2_32(pred, pred, -half, half-1);
				sym=gr_dec(&ec, nbypass[1]);
				{
					int upred=half-abs(pred), negmask=0;
					if(sym<=(upred<<1))
					{
						error=sym>>1^-(sym&1);
#ifdef ENABLE_BIASCORR
						negmask=-((ibias_corr<0)&(error!=-half));
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
				yuv[1]=error+pred;
				curr[1]=yuv[1]-offset;
				curr[5]=UPDATE_FORMULA(5);

				//dec V
				offset=(yuv[combination[7]]+yuv[combination[8]])>>combination[9];
				pred=preds[2]+offset;
				CLAMP2_32(pred, pred, -half, half-1);
				sym=gr_dec(&ec, nbypass[2]);
				{
					int upred=half-abs(pred), negmask=0;
					if(sym<=(upred<<1))
					{
						error=sym>>1^-(sym&1);
#ifdef ENABLE_BIASCORR
						negmask=-((ibias_corr<0)&(error!=-half));
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
				yuv[2]=error+pred;
				curr[2]=yuv[2]-offset;
				curr[6]=UPDATE_FORMULA(6);
				
				args->dst[idx+combination[3+0]]=yuv[0]+128;
				args->dst[idx+combination[3+1]]=yuv[1]+128;
				args->dst[idx+combination[3+2]]=yuv[2]+128;
#ifdef ENABLE_GUIDE
				if(args->test&&memcmp(args->dst+idx, args->src+idx, sizeof(char)*nch))
				{
					unsigned char orig[4]={0};
					memcpy(orig, args->src+idx, nch*sizeof(char));
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
			}
#ifdef ENABLE_SSE
			++curr_sse[0][0];
			++curr_sse[1][0];
			++curr_sse[2][0];
			curr_sse[0][1]+=curr[0]-preds[0];
			curr_sse[1][1]+=curr[1]-preds[1];
			curr_sse[2][1]+=curr[2]-preds[2];
#endif
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	if(args->fwd)
		gr_enc_flush(&ec);
}
int c05_codec(const char *srcfn, const char *dstfn)
{
	const int nch=3, depth=8;
	double t0;
	ArrayHandle src, dst;
	int headersize, printed;
	int iw, ih;
	const unsigned char *image, *imageend;
	unsigned char *image2;
	CodecID codec;
	int ncores;
	int xblocks, yblocks, nblocks, nthreads, coffset;
	ptrdiff_t start, memusage, argssize;
	ThreadArgs *args;
	int test, fwd;
	double esize;
	int usize;
	double t_analysis=0;
	
	t0=time_sec();
	src=load_file(srcfn, 1, 3, 1);
	headersize=header_read(src->data, (int)src->count, &iw, &ih, &codec);
	image=src->data+headersize;
	imageend=src->data+src->count;
	if(codec==CODEC_INVALID||codec==CODEC_PGM)
	{
		LOG_ERROR("Unsupported codec %d.\n", codec);
		array_free(&src);
		return 1;
	}
	else if(codec==CODEC_C01&&!dstfn)
	{
		LOG_ERROR(
			"Test mode expects PPM source.\n"
			"Decode mode expects destination filename."
		);
		return 1;
	}
	test=!dstfn;
	fwd=codec==CODEC_PPM;
	
	if(test)
		printf("%s \"%s\"  WH %d*%d\n", CODECNAME, srcfn, iw, ih);
#if 1
	if(test)
	{
		volatile long long sum=0;
		double t=time_sec();
		long long c=__rdtsc();
		unsigned state=0x10000;
		for(int k=1;k<(int)src->count-1;++k)
		{
			int cdf=k%0xFFFD+1, freq=k%(0xFFFF-cdf)+1;
			if(state>>16>(unsigned)freq)
			{
				sum+=2;
				state>>=16;
			}
			state=state/freq<<16|(cdf+state%freq);

			//int den=abs(src->data[k-1])+1;
			//int q=(int)((double)src->data[k]/den);
			//int r=src->data[k]%den;
			//sum+=(long long)q+r*src->data[k+1];
		}
		sum+=4;
		c=__rdtsc()-c;
		t=time_sec()-t;
		printf("DIV  %lf sec %lf MB/s  %lld cycles %lf cycles/byte  result %lld\n", t, src->count/(t*1024*1024), c, (double)c/src->count, sum);
	}
#endif
	ncores=query_cpu_cores();
	usize=iw*ih*3;
	xblocks=(iw+BLOCKSIZE-1)/BLOCKSIZE;
	yblocks=(ih+BLOCKSIZE-1)/BLOCKSIZE;
	nblocks=xblocks*yblocks, nthreads=MINVAR(nblocks, ncores);
	coffset=(int)sizeof(int)*nblocks;
	start=0;
	memusage=0;
	argssize=nthreads*sizeof(ThreadArgs);
	args=(ThreadArgs*)malloc(argssize);
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	esize=0;
	memusage+=argssize;
	memset(args, 0, argssize);
	if(fwd)
	{
		dst=0;
		printed=snprintf(g_buf, G_BUF_SIZE-1, "C01\n%d %d\n", iw, ih);
		array_append(&dst, g_buf, 1, printed, 1, 0, 0);
		start=array_append(&dst, 0, 1, coffset, 1, 0, 0);
		
		image2=0;
	}
	else//integrity check
	{
		dst=0;
		printed=snprintf(g_buf, G_BUF_SIZE-1, "P6\n%d %d\n255\n", iw, ih);
		array_append(&dst, g_buf, 1, printed, 1, 0, 0);
		array_append(&dst, 0, 1, usize, 1, 0, 0);

		//printed=0;
		start=coffset;
		for(int kt=0;kt<nblocks;++kt)
		{
			int size=0;
			memcpy(&size, image+sizeof(int)*kt, sizeof(int));
			start+=size;
		}
		if(image+start!=imageend)
			LOG_ERROR("Corrupt file");
		start=coffset;

		image2=(unsigned char*)malloc(usize);
		if(!image2)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		memset(image2, 0, usize);
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=image;
		arg->dst=fwd?0:dst->data+printed;
		arg->iw=iw;
		arg->ih=ih;

		arg->fwd=fwd;
		arg->test=test;
#ifdef DISABLE_MT
		arg->loud=test&&nblocks<MAXPRINTEDBLOCKS;
#else
		arg->loud=0;
#endif
	}
	for(int k2=0;k2<=test;++k2)
	{
		for(int kt=0;kt<nblocks;kt+=nthreads)
		{
			int nthreads2=MINVAR(kt+nthreads, nblocks)-kt;
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				int kx, ky;

				arg->blockidx=kx=kt+kt2;
				ky=kx/xblocks;
				kx%=xblocks;
				arg->x1=BLOCKSIZE*kx;
				arg->y1=BLOCKSIZE*ky;
				arg->x2=MINVAR(arg->x1+BLOCKSIZE, iw);
				arg->y2=MINVAR(arg->y1+BLOCKSIZE, ih);
				if(!fwd)
				{
					int size=0;
					memcpy(&size, image+sizeof(int)*((ptrdiff_t)kt+kt2), sizeof(int));
					arg->decstart=image+start;
					start+=size;
					arg->decend=image+start;
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
					ThreadArgs *arg=args+kt2;
					if(test)
					{
						int blocksize=((arg->x2-arg->x1)*(arg->y2-arg->y1)*nch*depth+7)>>3;
						int kx, ky;

						t_analysis+=arg->t_analysis;
						kx=kt+kt2;
						ky=kx/xblocks;
						kx%=xblocks;
						if(nblocks<MAXPRINTEDBLOCKS)
						{
							printf(
								"block %4d/%4d  XY %3d %3d  %4d*%4d:  %8d->%16lf->%8zd bytes (%+10.2lf)  %10.6lf%%  CR %10lf  %s\n",
								kt+kt2+1, nblocks,
								kx, ky,
								arg->y2-arg->y1,
								arg->x2-arg->x1,
								blocksize,
								arg->bestsize,
								arg->list.nbytes,
								arg->list.nbytes-arg->bestsize,
								100.*arg->list.nbytes/blocksize,
								(double)blocksize/arg->list.nbytes,
								rct_names[arg->bestrct]
							);
							//printf(
							//	"block %4d/%4d  XY %3d %3d  %4d*%4d:  %8d->%16lf->%8zd bytes (%+10.2lf)  %10.6lf%%  CR %10lf  %s %s %s %s\n",
							//	kt+kt2+1, nblocks,
							//	kx, ky,
							//	arg->y2-arg->y1,
							//	arg->x2-arg->x1,
							//	blocksize,
							//	arg->bestsize,
							//	arg->list.nbytes,
							//	arg->list.nbytes-arg->bestsize,
							//	100.*arg->list.nbytes/blocksize,
							//	(double)blocksize/arg->list.nbytes,
							//	rct_names[arg->bestrct],
							//	pred_names[arg->predidx[0]],
							//	pred_names[arg->predidx[1]],
							//	pred_names[arg->predidx[2]]
							//);
						}
						esize+=arg->bestsize;
#ifdef ABAC_PROFILESIZE
						for(int k=0;k<ABAC_TOKEN_BITS*3;++k)
							abac_csizes[k]+=arg->abac_csizes[k];
#endif
					}
					memcpy(dst->data+start+sizeof(int)*((ptrdiff_t)kt+kt2), &arg->list.nbytes, sizeof(int));
					blist_appendtoarray(&arg->list, &dst);
					blist_clear(&arg->list);
				}
			}
		}
		if(test)
		{
			ptrdiff_t usize=((ptrdiff_t)iw*ih*nch*depth+7)>>3;
			ptrdiff_t csize=dst->count;
			t0=time_sec()-t0;
			if(fwd)
			{
				printf("A %16.6lf sec  %16.6lf MB/s\n", t_analysis, usize/(t_analysis*1024*1024));
				printf("Best %15.2lf (%+13.2lf) bytes\n", esize, csize-esize);
				printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
				printf("Mem usage: ");
				print_size((double)memusage, 8, 4, 0, 0);
				printf("\n");
			}
			printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
			if(!fwd)
				compare_bufs_8(image2, src->data+headersize, iw, ih, nch, nch, CODECNAME, 0, 1);
		}
		if(!k2&&test)
		{
			int usize=iw*ih*3;
			fwd=0;
			image2=(unsigned char*)malloc(usize);
			if(!image2)
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			memset(image2, 0, usize);
			for(int kt=0;kt<nthreads;++kt)
			{
				ThreadArgs *arg=args+kt;
				arg->dst=image2;
				arg->fwd=0;
			}
			image=dst->data+printed;
			start=coffset;
		}
		t0=time_sec();
	}
	if(!test)
		save_file(dstfn, dst->data, dst->count, 1);
	if(image2)
		free(image2);
	free(args);
	array_free(&src);
	array_free(&dst);
	return 0;
}