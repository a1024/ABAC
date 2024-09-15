#include"codec.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT

	#define AC3_PREC

#define AC_IMPLEMENTATION
#include"entropy.h"

#define BLOCKSIZE 512
#define MAXPRINTEDBLOCKS 10
#define CLEVELS 18	//9
#define MIXBITS 14

#define OCHLIST\
	OCH(R)\
	OCH(G)\
	OCH(B)\
	OCH(RG)\
	OCH(RB)\
	OCH(GB)\
	OCH(GR)\
	OCH(BR)\
	OCH(BG)
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
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		3, 3, 3)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		3, 3, 1)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		3, 3, 0)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		3, 3, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		3, 3, 1)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		3, 3, 1)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		3, 3, 0)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		3, 0, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		3, 0, 1)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		3, 0, 1)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		3, 0, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		3, 0, 1)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		3, 0, 1)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		3, 0, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		3, 0, 1)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		3, 0, 1)
typedef enum _RCTIndex
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF, UOFF, VOFF) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][6]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF, UOFF, VOFF) {YIDX, UIDX, VIDX, YOFF, UOFF, VOFF},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF, UOFF, VOFF) #LABEL,
	RCTLIST
#undef  RCT
};

#define PREDLIST\
	PRED(W)\
	PRED(AV2)\
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

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
#define CONFIG_EXP 4
#define CONFIG_MSB 1
#define CONFIG_LSB 0
FORCEINLINE void quantize_pixel(int val, int *token, int *bypass, int *nbits)
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
		*nbits=lgv-(CONFIG_MSB+CONFIG_LSB);
		*bypass=val>>CONFIG_LSB&((1LL<<*nbits)-1);
	}
}
//FORCEINLINE int f28_mix4(int v00, int v01, int v10, int v11, int alphax, int alphay)
//{
//	//v00=v00*((1<<12)-alphax)+v01*alphax;
//	v00=((v00<<MIXBITS)+(v01-v00)*alphax)>>(MIXBITS-1);
//	v10=((v10<<MIXBITS)+(v11-v10)*alphax)>>(MIXBITS-1);
//	v00=((v00<<MIXBITS)+(v10-v00)*alphay)>>(MIXBITS-1);
//	return v00;
//}
typedef struct _ThreadArgs
{
	const unsigned char *src;
	unsigned char *dst;
	int iw, ih;

	int fwd, test, loud, x1, x2, y1, y2;
	int bufsize, histsize;
	short *pixels;
	int *hist;

	BList list;
	const unsigned char *decstart, *decend;

	int clevels, tlevels;

	//aux
	int blockidx;
	double bestsize;
	int bestrct, predidx[3];
} ThreadArgs;
static void block_thread(void *param)
{
	const int nch=3, depth=8, half=128;
	ThreadArgs *args=(ThreadArgs*)param;
	AC3 ec;
	const unsigned char *image=args->fwd?args->src:args->dst;
	unsigned char bestrct=0, combination[6]={0}, predidx[4]={0};
	int ystride=args->iw*3;
//	int nctx=args->clevels*args->clevels;
//	int cdfstride=args->tlevels+1;
//	int chsize=nctx*cdfstride;
	(void)depth;
	if(args->fwd)
	{
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		unsigned char predsel[OCH_COUNT]={0};
		int res=(args->x2-args->x1-1)/5*5*(args->y2-args->y1-1);
		
		memset(args->hist, 0, args->histsize);
		for(int ky=args->y1+1;ky<args->y2;++ky)//analysis loop
		{
			int kx=args->x1+1;
			const unsigned char *ptr=image+3*(args->iw*ky+kx);

			__m256i amin=_mm256_set1_epi16(-128);
			__m256i amax=_mm256_set1_epi16(127);
			__m128i half=_mm_set1_epi8(128);
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
				__m256i NW, N, W, curr, NW2, N2, W2, curr2;
				__m256i vmin, vmax, pred;
				{
					__m128i NW8	=_mm_load_si128((__m128i*)(ptr-1*ystride-1*3+0));
					__m128i N8	=_mm_load_si128((__m128i*)(ptr-1*ystride+0*3+0));
					__m128i W8	=_mm_load_si128((__m128i*)(ptr+0*ystride-1*3+0));
					__m128i curr8	=_mm_load_si128((__m128i*)(ptr+0*ystride+0*3+0));
					__m128i NW82=_mm_shuffle_epi8(NW8, shuf);
					__m128i N82=_mm_shuffle_epi8(N8, shuf);
					__m128i W82=_mm_shuffle_epi8(W8, shuf);
					__m128i curr82=_mm_shuffle_epi8(curr8, shuf);
					NW8	=_mm_xor_si128(NW8	, half);
					N8	=_mm_xor_si128(N8	, half);
					W8	=_mm_xor_si128(W8	, half);
					curr8	=_mm_xor_si128(curr8	, half);
					NW82	=_mm_xor_si128(NW82	, half);
					N82	=_mm_xor_si128(N82	, half);
					W82	=_mm_xor_si128(W82	, half);
					curr82	=_mm_xor_si128(curr82	, half);
					NW	=_mm256_cvtepi8_epi16(NW8);
					N	=_mm256_cvtepi8_epi16(N8);
					W	=_mm256_cvtepi8_epi16(W8);
					curr	=_mm256_cvtepi8_epi16(curr8);
					NW2	=_mm256_cvtepi8_epi16(NW82);
					N2	=_mm256_cvtepi8_epi16(N82);
					W2	=_mm256_cvtepi8_epi16(W82);
					curr2	=_mm256_cvtepi8_epi16(curr82);
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
				//W
				pred=_mm256_sub_epi16(curr, W);
				UPDATE(
					PRED_W,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(_mm256_sub_epi16(W, W2), curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_W,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(_mm256_sub_epi16(W2, W), curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_W,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);

				//(N+W)/2
				pred=_mm256_srai_epi16(_mm256_add_epi16(N, W), 1);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV2,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_srai_epi16(_mm256_add_epi16(_mm256_sub_epi16(N, N2), _mm256_sub_epi16(W, W2)), 1);
				pred=_mm256_add_epi16(pred, curr2);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_AV2,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_srai_epi16(_mm256_add_epi16(_mm256_sub_epi16(N2, N), _mm256_sub_epi16(W2, W)), 1);
				pred=_mm256_add_epi16(pred, curr);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(curr2, pred);
				UPDATE(
					PRED_AV2,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);

				//CG
				vmin=_mm256_min_epi16(N, W);
				vmax=_mm256_max_epi16(N, W);
				pred=_mm256_sub_epi16(_mm256_add_epi16(N, W), NW);
				pred=_mm256_max_epi16(pred, vmin);
				pred=_mm256_min_epi16(pred, vmax);
				pred=_mm256_sub_epi16(curr, pred);
				UPDATE(
					PRED_CG,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				{
					__m256i N3=_mm256_sub_epi16(N, N2);
					__m256i W3=_mm256_sub_epi16(W, W2);
					__m256i NW3=_mm256_sub_epi16(NW, NW2);

					vmin=_mm256_min_epi16(N3, W3);
					vmax=_mm256_max_epi16(N3, W3);
					pred=_mm256_sub_epi16(_mm256_add_epi16(N3, W3), NW3);
					pred=_mm256_max_epi16(pred, vmin);
					pred=_mm256_min_epi16(pred, vmax);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_CG,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
				}
				{
					__m256i N3=_mm256_sub_epi16(N2, N);
					__m256i W3=_mm256_sub_epi16(W2, W);
					__m256i NW3=_mm256_sub_epi16(NW2, NW);

					vmin=_mm256_min_epi16(N3, W3);
					vmax=_mm256_max_epi16(N3, W3);
					pred=_mm256_sub_epi16(_mm256_add_epi16(N3, W3), NW3);
					pred=_mm256_max_epi16(pred, vmin);
					pred=_mm256_min_epi16(pred, vmax);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_CG,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
				}
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
		if(args->loud)
		{
			printf("%3d  X%5d~%5d Y%5d~%5d  best %12.2lf bytes  %s [YUV: %s %s %s]\n",
				args->blockidx,
				args->x1, args->x2,
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
				const unsigned char *group=rct_combinations[kt];
				double csize=
					csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
					csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
					csizes[group[2]*PRED_COUNT+predsel[group[2]]];
				printf("%12.2lf %c  %-10s %-10s %-10s %-10s\n",
					csize,
					kt==bestrct?'*':' ',
					rct_names[kt],
					pred_names[predsel[group[0]]],
					pred_names[predsel[group[1]]],
					pred_names[predsel[group[2]]]
				);
			}
		}
		blist_init(&args->list);
		ac3_enc_init(&ec, &args->list);
		ac3_enc_bypass_NPOT(&ec, bestrct, RCT_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[0], PRED_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[1], PRED_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[2], PRED_COUNT);
	}
	else
	{
		ac3_dec_init(&ec, args->decstart, args->decend);
		bestrct=ac3_dec_bypass_NPOT(&ec, RCT_COUNT);
		predidx[0]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[1]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[2]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
	}
	memset(args->hist, 0, args->histsize);
	//for(int ky=0;ky<args->clevels;++ky)
	//{
	//	for(int kx=0;kx<args->clevels;++kx)
	//	{
	//		static const int init_freqs[]={32, 8, 6, 4, 3, 2, 1};
	//		int *curr_hist=args->hist+cdfstride*(args->clevels*ky+kx);
	//		int sum=0;
	//		for(int ks=0;ks<args->tlevels;++ks)
	//		{
	//			int freq=init_freqs[MINVAR(ks, (int)_countof(init_freqs)-1)];
	//			sum+=curr_hist[ks]=freq;
	//		}
	//		curr_hist[args->tlevels]=sum;
	//	}
	//}
	//memfill(args->hist+chsize, args->hist, sizeof(int)*chsize*(nch-1LL), sizeof(int)*chsize);
	memset(args->pixels, 0, args->bufsize);
	__m128i swap16=_mm_set_epi8(
		13, 12, 15, 14,  9,  8, 11, 10,  5,  4,  7,  6,  1,  0,  3,  2
	//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
	);
	__m128i cmax=_mm_set1_epi16(CLEVELS-2);
	//__m128i one=_mm_set1_epi16(1);
	//__m128i amax=_mm_set1_epi16(1<<MIXBITS);
	__m128i predmaskW=_mm_set_epi16(
		0, 0, 0, 0, 0,
		-(predidx[2]==PRED_W),
		-(predidx[1]==PRED_W),
		-(predidx[0]==PRED_W)
	);
	__m128i predmaskAV2=_mm_set_epi16(
		0, 0, 0, 0, 0,
		-(predidx[2]==PRED_AV2),
		-(predidx[1]==PRED_AV2),
		-(predidx[0]==PRED_AV2)
	);
	__m256i mclevels=_mm256_set1_epi64x(args->clevels);
	__m256i mclevels2=_mm256_set1_epi64x(args->clevels*2);
	__m256i mcdfstride4=_mm256_set1_epi64x((args->tlevels+1LL)*sizeof(int));
	__m256i mgyoffset=_mm256_set_epi64x(1, 1, 0, 0);
	__m256i mgxoffset=_mm256_set_epi64x(1, 0, 1, 0);
	__m256i mhistptr=_mm256_set1_epi64x((long long)args->hist);
	__m256i mtlevels=_mm256_set1_epi64x(args->tlevels);
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
		int token=0, bypass=0, nbits=0;
		int pred=0, error=0, sym=0;
		const unsigned char *combination=rct_combinations[bestrct];
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int idx=nch*(args->iw*ky+kx);
			short
				*NNN	=rows[3]+0*4*2,
			//	*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
			//	*NNEE	=rows[2]+2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*NEE	=rows[1]+2*4*2,
				*NEEE	=rows[1]+3*4*2,
				*WWW	=rows[0]-3*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
			(void)NEEE;
			if(ky<=args->y1+2)
			{
				if(ky<=args->y1+1)
				{
					if(ky==args->y1)
						NEEE=NEE=NE=NW=N=W;
					NN=N;
					NNE=NE;
				}
				NNN=NN;
			}
			if(kx<=args->x1+2)
			{
				if(kx<=args->x1+1)
				{
					if(kx<=args->x1)
						NW=W=N;
					WW=W;
				}
				WWW=WW;
			}
			if(kx>=args->x2-3)
			{
				if(kx>=args->x2-2)
				{
					if(kx>=args->x2-1)
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
				switch(bestrct)
				{
				case RCT_R_G_B:
				case RCT_R_G_BG:
				case RCT_R_G_BR:
				case RCT_R_GR_BR:
				case RCT_R_GR_BG:
				//case RCT_R_G_B2:
				//case RCT_R_GR_B2:
					yuv[0]=args->src[idx+0]-128;
					yuv[1]=args->src[idx+1]-128;
					yuv[2]=args->src[idx+2]-128;
					break;
				case RCT_G_B_RG:
				case RCT_G_B_RB:
				case RCT_G_BG_RG:
				case RCT_G_BG_RB:
				//case RCT_G_B_R2:
				//case RCT_G_BG_R2:
					yuv[0]=args->src[idx+1]-128;
					yuv[1]=args->src[idx+2]-128;
					yuv[2]=args->src[idx+0]-128;
					break;
				case RCT_B_R_GR:
				case RCT_B_R_GB:
				case RCT_B_RB_GB:
				case RCT_B_RB_GR:
				//case RCT_B_RB_G2:
					yuv[0]=args->src[idx+2]-128;
					yuv[1]=args->src[idx+0]-128;
					yuv[2]=args->src[idx+1]-128;
					break;
				case RCT_G_RG_BR:
				//case RCT_G_RG_B2:
					yuv[0]=args->src[idx+1]-128;
					yuv[1]=args->src[idx+0]-128;
					yuv[2]=args->src[idx+2]-128;
					break;
				case RCT_B_GB_RG:
				//case RCT_B_GB_R2:
					yuv[0]=args->src[idx+2]-128;
					yuv[1]=args->src[idx+1]-128;
					yuv[2]=args->src[idx+0]-128;
					break;
				case RCT_R_BR_GB:
				//case RCT_R_B_G2:
				//case RCT_R_BR_G2:
					yuv[0]=args->src[idx+0]-128;
					yuv[1]=args->src[idx+2]-128;
					yuv[2]=args->src[idx+1]-128;
					break;
				}
			}
			__m128i mNNN	=_mm_load_si128((__m128i*)NNN);
			__m128i mNN	=_mm_load_si128((__m128i*)NN);
			__m128i mNNE	=_mm_load_si128((__m128i*)NNE);
			__m128i mNW	=_mm_load_si128((__m128i*)NW);
			__m128i mN	=_mm_load_si128((__m128i*)N);
			__m128i mNE	=_mm_load_si128((__m128i*)NE);
			__m128i mWWW	=_mm_load_si128((__m128i*)WWW);
			__m128i mWW	=_mm_load_si128((__m128i*)WW);
			__m128i mW	=_mm_load_si128((__m128i*)W);
			__m128i t0=_mm_abs_epi16(_mm_sub_epi16(mW, mWW));//abs(W-WW)+abs(N-NW)+abs(NE-N)+abs(eWWW)+abs(eWW)+2*abs(eW)
			__m128i t1=_mm_abs_epi16(_mm_sub_epi16(mN, mNN));//abs(N-NN)+abs(W-NW)+abs(NE-NNE)+abs(eNNN)+abs(eNN)+2*abs(eN)
			t0=_mm_add_epi16(t0, _mm_abs_epi16(_mm_sub_epi16(mN, mNW)));
			t1=_mm_add_epi16(t1, _mm_abs_epi16(_mm_sub_epi16(mW, mNW)));
			t0=_mm_add_epi16(t0, _mm_abs_epi16(_mm_sub_epi16(mNE, mN)));
			t1=_mm_add_epi16(t1, _mm_abs_epi16(_mm_sub_epi16(mNE, mNNE)));
			t0=_mm_add_epi16(t0, _mm_shuffle_epi8(_mm_abs_epi16(mWWW), swap16));
			t1=_mm_add_epi16(t1, _mm_shuffle_epi8(_mm_abs_epi16(mNNN), swap16));
			t0=_mm_add_epi16(t0, _mm_shuffle_epi8(_mm_abs_epi16(mWW), swap16));
			t1=_mm_add_epi16(t1, _mm_shuffle_epi8(_mm_abs_epi16(mNN), swap16));
			t0=_mm_add_epi16(t0, _mm_slli_epi16(_mm_shuffle_epi8(_mm_abs_epi16(mW), swap16), 1));
			t1=_mm_add_epi16(t1, _mm_slli_epi16(_mm_shuffle_epi8(_mm_abs_epi16(mN), swap16), 1));
		//	t0=_mm_slli_epi16(t0, 2);
		//	t1=_mm_slli_epi16(t1, 2);
		//	t0=_mm_add_epi16(t0, one);
		//	t1=_mm_add_epi16(t1, one);
			ALIGN(16) unsigned short grads[16];
			_mm_store_si128((__m128i*)grads+0, t0);//{0, 2, 4, 8, 10, 12}
			_mm_store_si128((__m128i*)grads+1, t1);
#if 1
			ALIGN(16) int g2[]=
			{
				grads[ 0]*grads[ 0]+1,
				grads[ 2]*grads[ 2]+1,
				grads[ 4]*grads[ 4]+1,
				grads[ 8]*grads[ 8]+1,
				grads[10]*grads[10]+1,
				grads[12]*grads[12]+1,
			};
			grads[ 8]=FLOOR_LOG2(g2[0]);
			grads[ 9]=FLOOR_LOG2(g2[1]);
			grads[10]=FLOOR_LOG2(g2[2]);
			grads[11]=FLOOR_LOG2(g2[3]);
			grads[12]=FLOOR_LOG2(g2[4]);
			grads[13]=FLOOR_LOG2(g2[5]);
			//ALIGN(32) long long alphas[]=//((x<<a)+(1<<b>>1))>>b = x<<(a-b)	a>=b
			//{
			//	1<<MIXBITS>>2,
			//	1<<MIXBITS>>2,
			//	1<<MIXBITS>>2,
			//	0,
			//	1<<MIXBITS>>2,
			//	1<<MIXBITS>>2,
			//	1<<MIXBITS>>2,
			//	0,
			//	//(g2[0]-(1<<grads[ 8]))<<(24-grads[ 8])>>(24-MIXBITS),//yax
			//	//(g2[1]-(1<<grads[ 9]))<<(24-grads[ 9])>>(24-MIXBITS),//uax
			//	//(g2[2]-(1<<grads[10]))<<(24-grads[10])>>(24-MIXBITS),//vax
			//	//0,
			//	//(g2[3]-(1<<grads[11]))<<(24-grads[11])>>(24-MIXBITS),//yay
			//	//(g2[4]-(1<<grads[12]))<<(24-grads[12])>>(24-MIXBITS),//uay
			//	//(g2[5]-(1<<grads[13]))<<(24-grads[13])>>(24-MIXBITS),//vay
			//	//0,
			//};
#else
			grads[0]=grads[ 0]*grads[ 0]+1;//vx
			grads[1]=grads[ 2]*grads[ 2]+1;
			grads[2]=grads[ 4]*grads[ 4]+1;
			grads[3]=grads[ 8]*grads[ 8]+1;//vy
			grads[4]=grads[10]*grads[10]+1;
			grads[5]=grads[12]*grads[12]+1;//{0, 1, 2, 3, 4, 5}
		//	grads[0]=grads[ 0];//vx
		//	grads[1]=grads[ 2];
		//	grads[2]=grads[ 4];
		//	grads[3]=grads[ 8];//vy
		//	grads[4]=grads[10];
		//	grads[5]=grads[12];//{0, 1, 2, 3, 4, 5}
			grads[ 8]=FLOOR_LOG2(grads[0]);
			grads[ 9]=FLOOR_LOG2(grads[1]);
			grads[10]=FLOOR_LOG2(grads[2]);
			grads[11]=FLOOR_LOG2(grads[3]);
			grads[12]=FLOOR_LOG2(grads[4]);
			grads[13]=FLOOR_LOG2(grads[5]);
			//if(ky==0&&kx==1)//
			//	printf("");
			ALIGN(32) long long alphas[]=//((x<<a)+(1<<b>>1))>>b = x<<(a-b)	a>=b
			{
				(grads[0]-(1<<grads[ 8]))<<(24-grads[ 8])>>(24-MIXBITS),//yax
				(grads[1]-(1<<grads[ 9]))<<(24-grads[ 9])>>(24-MIXBITS),//uax
				(grads[2]-(1<<grads[10]))<<(24-grads[10])>>(24-MIXBITS),//vax
				0,
				(grads[3]-(1<<grads[11]))<<(24-grads[11])>>(24-MIXBITS),//yay
				(grads[4]-(1<<grads[12]))<<(24-grads[12])>>(24-MIXBITS),//uay
				(grads[5]-(1<<grads[13]))<<(24-grads[13])>>(24-MIXBITS),//vay
				0,
			};
#endif
			//if(
			//	alphas[0]<0||alphas[0]>=(1<<MIXBITS)||
			//	alphas[1]<0||alphas[1]>=(1<<MIXBITS)||
			//	alphas[2]<0||alphas[2]>=(1<<MIXBITS)||
			//	alphas[3]<0||alphas[3]>=(1<<MIXBITS)||
			//	alphas[4]<0||alphas[4]>=(1<<MIXBITS)||
			//	alphas[5]<0||alphas[5]>=(1<<MIXBITS)
			//)
			//	printf("");
			t1=_mm_load_si128((__m128i*)grads+1);
			t1=_mm_min_epi16(t1, cmax);
			_mm_store_si128((__m128i*)grads+1, t1);
#if 1
			ALIGN(32) int *ctxptr[12];
			{
				__m256i gx0=_mm256_set1_epi64x(grads[ 8]);
				__m256i gx1=_mm256_set1_epi64x(grads[ 9]);
				__m256i gx2=_mm256_set1_epi64x(grads[10]);
				__m256i gy0=_mm256_set1_epi64x(grads[11]);
				__m256i gy1=_mm256_set1_epi64x(grads[12]);
				__m256i gy2=_mm256_set1_epi64x(grads[13]);
				gy0=_mm256_add_epi64(gy0, mgyoffset);
				gy1=_mm256_add_epi64(gy1, mgyoffset);
				gy2=_mm256_add_epi64(gy2, mgyoffset);
				gy1=_mm256_add_epi64(gy1, mclevels);
				gy2=_mm256_add_epi64(gy2, mclevels2);
				gy0=_mm256_mul_epi32(gy0, mclevels);
				gy1=_mm256_mul_epi32(gy1, mclevels);
				gy2=_mm256_mul_epi32(gy2, mclevels);
				gx0=_mm256_add_epi64(gx0, mgxoffset);
				gx1=_mm256_add_epi64(gx1, mgxoffset);
				gx2=_mm256_add_epi64(gx2, mgxoffset);
				gy0=_mm256_add_epi64(gy0, gx0);
				gy1=_mm256_add_epi64(gy1, gx1);
				gy2=_mm256_add_epi64(gy2, gx2);
				gy0=_mm256_mul_epi32(gy0, mcdfstride4);
				gy1=_mm256_mul_epi32(gy1, mcdfstride4);
				gy2=_mm256_mul_epi32(gy2, mcdfstride4);
				gy0=_mm256_add_epi64(gy0, mhistptr);
				gy1=_mm256_add_epi64(gy1, mhistptr);
				gy2=_mm256_add_epi64(gy2, mhistptr);
				_mm256_store_si256((__m256i*)ctxptr+0, gy0);
				_mm256_store_si256((__m256i*)ctxptr+1, gy1);
				_mm256_store_si256((__m256i*)ctxptr+2, gy2);
			}
#else
			int *ctxptr[]=
			{
				args->hist+cdfstride*(nctx*0+args->clevels*(grads[11]+0)+grads[ 8]+0),//optimize
				args->hist+cdfstride*(nctx*0+args->clevels*(grads[11]+0)+grads[ 8]+1),
				args->hist+cdfstride*(nctx*0+args->clevels*(grads[11]+1)+grads[ 8]+0),
				args->hist+cdfstride*(nctx*0+args->clevels*(grads[11]+1)+grads[ 8]+1),
				args->hist+cdfstride*(nctx*1+args->clevels*(grads[12]+0)+grads[ 9]+0),
				args->hist+cdfstride*(nctx*1+args->clevels*(grads[12]+0)+grads[ 9]+1),
				args->hist+cdfstride*(nctx*1+args->clevels*(grads[12]+1)+grads[ 9]+0),
				args->hist+cdfstride*(nctx*1+args->clevels*(grads[12]+1)+grads[ 9]+1),
				args->hist+cdfstride*(nctx*2+args->clevels*(grads[13]+0)+grads[10]+0),
				args->hist+cdfstride*(nctx*2+args->clevels*(grads[13]+0)+grads[10]+1),
				args->hist+cdfstride*(nctx*2+args->clevels*(grads[13]+1)+grads[10]+0),
				args->hist+cdfstride*(nctx*2+args->clevels*(grads[13]+1)+grads[10]+1),
			};
#endif
			//t0=_mm_load_si128((__m128i*)alphas);//no need
			//t0=_mm_min_epi16(t0, amax);
			//t0=_mm_max_epi16(t0, _mm_setzero_si128());
			//_mm_store_si128((__m128i*)alphas, t0);
#if 1
			ALIGN(32) long long dens[4];
			{
				ALIGN(32) long long temp[]=
				{
					ctxptr[ 0][args->tlevels],//y00
					ctxptr[ 4][args->tlevels],//u00
					ctxptr[ 8][args->tlevels],//v00
					0,
					ctxptr[ 1][args->tlevels],//y01
					ctxptr[ 5][args->tlevels],//u01
					ctxptr[ 9][args->tlevels],//v01
					0,
					ctxptr[ 2][args->tlevels],//y10
					ctxptr[ 6][args->tlevels],//u10
					ctxptr[10][args->tlevels],//v10
					0,
					ctxptr[ 3][args->tlevels],//y11
					ctxptr[ 7][args->tlevels],//u11
					ctxptr[11][args->tlevels],//v11
					0,
				};
#if 1
				__m256i v00=_mm256_load_si256((__m256i*)temp+0);
				__m256i v01=_mm256_load_si256((__m256i*)temp+1);
				__m256i v10=_mm256_load_si256((__m256i*)temp+2);
				__m256i v11=_mm256_load_si256((__m256i*)temp+3);
				v00=_mm256_add_epi32(v00, v01);
				v00=_mm256_add_epi32(v00, v10);
				v00=_mm256_add_epi32(v00, v11);
			//	v00=_mm256_srli_epi32(v00, 1);
				v00=_mm256_add_epi32(v00, mtlevels);
				_mm256_store_si256((__m256i*)dens, v00);
#else
				__m256i ax=_mm256_load_si256((__m256i*)alphas+0);
				__m256i ay=_mm256_load_si256((__m256i*)alphas+1);
				__m256i v00=_mm256_load_si256((__m256i*)temp+0);
				__m256i v01=_mm256_load_si256((__m256i*)temp+1);
				__m256i v10=_mm256_load_si256((__m256i*)temp+2);
				__m256i v11=_mm256_load_si256((__m256i*)temp+3);
				__m256i tmp0=_mm256_sub_epi32(v01, v00);
				__m256i tmp1=_mm256_sub_epi32(v11, v10);
				tmp0=_mm256_mul_epi32(tmp0, ax);
				tmp1=_mm256_mul_epi32(tmp1, ax);
				v00=_mm256_slli_epi32(v00, MIXBITS);
				v10=_mm256_slli_epi32(v10, MIXBITS);
				tmp0=_mm256_add_epi64(tmp0, v00);
				tmp1=_mm256_add_epi64(tmp1, v10);
				tmp0=_mm256_srli_epi32(tmp0, MIXBITS-1);
				tmp1=_mm256_srli_epi32(tmp1, MIXBITS-1);

				v00=_mm256_sub_epi32(tmp1, tmp0);
				v00=_mm256_mul_epi32(v00, ay);
				tmp0=_mm256_slli_epi32(tmp0, MIXBITS);
				v00=_mm256_add_epi64(v00, tmp0);
				v00=_mm256_srli_epi32(v00, MIXBITS-1);
				v00=_mm256_add_epi32(v00, mtlevels);
				_mm256_store_si256((__m256i*)dens, v00);
#endif
			}
#else
			int dens[]=
			{
				ctxptr[ 0][args->tlevels],
				ctxptr[ 1][args->tlevels],
				ctxptr[ 2][args->tlevels],
				ctxptr[ 3][args->tlevels],
				ctxptr[ 4][args->tlevels],
				ctxptr[ 5][args->tlevels],
				ctxptr[ 6][args->tlevels],
				ctxptr[ 7][args->tlevels],
				ctxptr[ 8][args->tlevels],
				ctxptr[ 9][args->tlevels],
				ctxptr[10][args->tlevels],
				ctxptr[11][args->tlevels],
			};
			dens[0+0*4]=((dens[0+0*4]<<MIXBITS)+(dens[1+0*4]-dens[0+0*4])*alphas[0])>>(MIXBITS-1);//optimize
			dens[2+0*4]=((dens[2+0*4]<<MIXBITS)+(dens[3+0*4]-dens[2+0*4])*alphas[0])>>(MIXBITS-1);
			dens[0+1*4]=((dens[0+1*4]<<MIXBITS)+(dens[1+1*4]-dens[0+1*4])*alphas[1])>>(MIXBITS-1);
			dens[2+1*4]=((dens[2+1*4]<<MIXBITS)+(dens[3+1*4]-dens[2+1*4])*alphas[1])>>(MIXBITS-1);
			dens[0+2*4]=((dens[0+2*4]<<MIXBITS)+(dens[1+2*4]-dens[0+2*4])*alphas[2])>>(MIXBITS-1);
			dens[2+2*4]=((dens[2+2*4]<<MIXBITS)+(dens[3+2*4]-dens[2+2*4])*alphas[2])>>(MIXBITS-1);
			dens[0+0*4]=((dens[0+0*4]<<MIXBITS)+(dens[2+0*4]-dens[0+0*4])*alphas[3])>>(MIXBITS-1);//den0
			dens[0+1*4]=((dens[0+1*4]<<MIXBITS)+(dens[2+1*4]-dens[0+1*4])*alphas[4])>>(MIXBITS-1);//den1
			dens[0+2*4]=((dens[0+2*4]<<MIXBITS)+(dens[2+2*4]-dens[0+2*4])*alphas[5])>>(MIXBITS-1);//den2
			dens[0+0*4]+=args->tlevels;
			dens[0+1*4]+=args->tlevels;
			dens[0+2*4]+=args->tlevels;
#endif
			ALIGN(16) short preds[8];
			__m128i mav2=_mm_add_epi16(mN, mW);
			__m128i mp=_mm_sub_epi16(mav2, mNW);
			mav2=_mm_srai_epi16(mav2, 1);
			t0=_mm_min_epi16(mN, mW);
			t1=_mm_max_epi16(mN, mW);
			mp=_mm_max_epi16(mp, t0);
			mp=_mm_min_epi16(mp, t1);
			mp=_mm_blendv_epi8(mp, mav2, predmaskAV2);//SSE4.1
			mp=_mm_blendv_epi8(mp, mW, predmaskW);
			_mm_store_si128((__m128i*)preds, mp);
			for(int kc=0;kc<nch;++kc)
			{
				int kc2=kc<<1;
				int offset=yuv[combination[kc+3]];
				int cdf=0, freq=0;
				int den=(int)dens[kc];
			//	int den=dens[kc<<2];
				int *curr_hist00=ctxptr[kc<<2|0];
				int *curr_hist01=ctxptr[kc<<2|1];
				int *curr_hist10=ctxptr[kc<<2|2];
				int *curr_hist11=ctxptr[kc<<2|3];
			//	__m256i ax=_mm256_set1_epi64x(alphas[kc+0]);
			//	__m256i ay=_mm256_set1_epi64x(alphas[kc+4]);
			//	int alphax=(int)alphas[kc+0];
			//	int alphay=(int)alphas[kc+4];
				pred=preds[kc2];
#if 0
				int curr_hist[4];
				int
					vx=(abs(W[kc2]-WW[kc2])+abs(N[kc2]-NW[kc2])+abs(NE[kc2]-N  [kc2])+abs(WWW[kc2+1])+abs(WW[kc2+1])+abs(W[kc2+1])*2)<<10>>depth,
					vy=(abs(W[kc2]-NW[kc2])+abs(N[kc2]-NN[kc2])+abs(NE[kc2]-NNE[kc2])+abs(NNN[kc2+1])+abs(NN[kc2+1])+abs(N[kc2+1])*2)<<10>>depth;
				int qeN=FLOOR_LOG2(vy+1);
				int qeW=FLOOR_LOG2(vx+1);
				int alphax2=(((vx+1-(1<<qeW))<<MIXBITS)+(1<<qeW>>1))>>qeW;//qe < MIXBITS, no need to add half
				int alphay2=(((vy+1-(1<<qeN))<<MIXBITS)+(1<<qeN>>1))>>qeN;
				qeN=MINVAR(qeN, CLEVELS-2);
				qeW=MINVAR(qeW, CLEVELS-2);
				curr_hist[0]=args->hist+cdfstride*(nctx*kc+args->clevels*(qeN+0)+qeW+0);
				curr_hist[1]=args->hist+cdfstride*(nctx*kc+args->clevels*(qeN+0)+qeW+1);
				curr_hist[2]=args->hist+cdfstride*(nctx*kc+args->clevels*(qeN+1)+qeW+0);
				curr_hist[3]=args->hist+cdfstride*(nctx*kc+args->clevels*(qeN+1)+qeW+1);
				//if(alphax!=alphax2||alphay!=alphay2)
				//	LOG_ERROR("");
				CLAMP2_32(alphax2, 0, alphax2, 1<<MIXBITS);
				CLAMP2_32(alphay2, 0, alphay2, 1<<MIXBITS);
#endif
#if 0
//#define MIX4(X) f28_mix4(curr_hist00[X], curr_hist01[X], curr_hist10[X], curr_hist11[X], alphax, alphay)
#define MIX4(DST, X)\
	do\
	{\
		int v00=curr_hist00[X];\
		int v01=curr_hist01[X];\
		int v10=curr_hist10[X];\
		int v11=curr_hist11[X];\
		v00=((v00<<MIXBITS)+(v01-v00)*alphax)>>(MIXBITS-1);\
		v10=((v10<<MIXBITS)+(v11-v10)*alphax)>>(MIXBITS-1);\
		v00=((v00<<MIXBITS)+(v10-v00)*alphay)>>(MIXBITS-1);\
		DST=v00+1;\
	}while(0)
#endif
			//	den=MIX4(args->tlevels);
			//	switch(predidx[kc])
			//	{
			//	case PRED_W:
			//		pred=W[kc2];
			//		break;
			//	case PRED_CG:
			//		MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
			//		break;
			//	}
				pred+=offset;
				CLAMP2(pred, -128, 127);
				//if(ky==480&&kx==109&&kc==2)//
				//if(ky==147&&kx==317&&kc==1)//
				//	printf("");
				if(args->fwd)
				{
					curr[kc2+0]=yuv[kc];
					curr[kc2+1]=error=yuv[kc]-pred;
					{
						int upred=half-abs(pred), aval=abs(error);
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
#ifdef _DEBUG
					if(token>=args->tlevels)
						LOG_ERROR("YXC %d %d %d  token %d/%d", ky, kx, kc, token, args->tlevels);
#endif
#if 1
					cdf=0;
					for(int ks=0;;ks+=4)
					{
						ALIGN(32) long long freq2[4];
						__m256i v00=_mm256_cvtepi32_epi64(_mm_loadu_si128((__m128i*)(curr_hist00+ks)));
						__m256i v01=_mm256_cvtepi32_epi64(_mm_loadu_si128((__m128i*)(curr_hist01+ks)));
						__m256i v10=_mm256_cvtepi32_epi64(_mm_loadu_si128((__m128i*)(curr_hist10+ks)));
						__m256i v11=_mm256_cvtepi32_epi64(_mm_loadu_si128((__m128i*)(curr_hist11+ks)));
#if 1
						v00=_mm256_add_epi32(v00, v01);
						v00=_mm256_add_epi32(v00, v10);
						v00=_mm256_add_epi32(v00, v11);
					//	v00=_mm256_srli_epi32(v00, 1);
#else
						__m256i tmp0=_mm256_sub_epi32(v01, v00);
						__m256i tmp1=_mm256_sub_epi32(v11, v10);
						tmp0=_mm256_mul_epi32(tmp0, ax);
						tmp1=_mm256_mul_epi32(tmp1, ax);
						v00=_mm256_slli_epi32(v00, MIXBITS);
						v10=_mm256_slli_epi32(v10, MIXBITS);
						tmp0=_mm256_add_epi64(tmp0, v00);
						tmp1=_mm256_add_epi64(tmp1, v10);
						tmp0=_mm256_srli_epi32(tmp0, MIXBITS-1);
						tmp1=_mm256_srli_epi32(tmp1, MIXBITS-1);

						v00=_mm256_sub_epi32(tmp1, tmp0);
						v00=_mm256_mul_epi32(v00, ay);
						tmp0=_mm256_slli_epi32(tmp0, MIXBITS);
						v00=_mm256_add_epi64(v00, tmp0);
						v00=_mm256_srli_epi32(v00, MIXBITS-1);
#endif
						v00=_mm256_add_epi32(v00, _mm256_set1_epi64x(1));
						_mm256_store_si256((__m256i*)freq2, v00);
					//	MIX4(freq2[0], ks+0);
					//	MIX4(freq2[1], ks+1);
					//	MIX4(freq2[2], ks+2);
					//	MIX4(freq2[3], ks+3);

						freq=(int)freq2[0];
						if(ks+0>=token)
							break;
						cdf+=freq;

						freq=(int)freq2[1];
						if(ks+1>=token)
							break;
						cdf+=freq;

						freq=(int)freq2[2];
						if(ks+2>=token)
							break;
						cdf+=freq;

						freq=(int)freq2[3];
						if(ks+3>=token)
							break;
						cdf+=freq;

					}
#else
					cdf=0;
					for(int ks=0;;++ks)
					{
						MIX4(freq, ks);
						if(ks>=token)
							break;
						cdf+=freq;
					}
#endif
					ac3_enc_update_NPOT(&ec, cdf, freq, den);
					if(nbits)
						ac3_enc_bypass(&ec, bypass, nbits);
				}
				else
				{
					unsigned code=ac3_dec_getcdf_NPOT(&ec, den);
#if 1
					cdf=0;
					token=0;
					for(;;)
					{
						ALIGN(32) long long freq2[4];
						__m256i v00=_mm256_cvtepi32_epi64(_mm_loadu_si128((__m128i*)(curr_hist00+token)));
						__m256i v01=_mm256_cvtepi32_epi64(_mm_loadu_si128((__m128i*)(curr_hist01+token)));
						__m256i v10=_mm256_cvtepi32_epi64(_mm_loadu_si128((__m128i*)(curr_hist10+token)));
						__m256i v11=_mm256_cvtepi32_epi64(_mm_loadu_si128((__m128i*)(curr_hist11+token)));
#if 1
						v00=_mm256_add_epi32(v00, v01);
						v00=_mm256_add_epi32(v00, v10);
						v00=_mm256_add_epi32(v00, v11);
					//	v00=_mm256_srli_epi32(v00, 1);
#else
						__m256i tmp0=_mm256_sub_epi32(v01, v00);
						__m256i tmp1=_mm256_sub_epi32(v11, v10);
						tmp0=_mm256_mul_epi32(tmp0, ax);
						tmp1=_mm256_mul_epi32(tmp1, ax);
						v00=_mm256_slli_epi32(v00, MIXBITS);
						v10=_mm256_slli_epi32(v10, MIXBITS);
						tmp0=_mm256_add_epi64(tmp0, v00);
						tmp1=_mm256_add_epi64(tmp1, v10);
						tmp0=_mm256_srli_epi32(tmp0, MIXBITS-1);
						tmp1=_mm256_srli_epi32(tmp1, MIXBITS-1);

						v00=_mm256_sub_epi32(tmp1, tmp0);
						v00=_mm256_mul_epi32(v00, ay);
						tmp0=_mm256_slli_epi32(tmp0, MIXBITS);
						v00=_mm256_add_epi64(v00, tmp0);
						v00=_mm256_srli_epi32(v00, MIXBITS-1);
#endif
						v00=_mm256_add_epi32(v00, _mm256_set1_epi64x(1));
						_mm256_store_si256((__m256i*)freq2, v00);
					//	MIX4(freq2[0], token+0);
					//	MIX4(freq2[1], token+1);
					//	MIX4(freq2[2], token+2);
					//	MIX4(freq2[3], token+3);

						unsigned cdf2;
						freq=(int)freq2[0];
						cdf2=cdf+freq;
						if(cdf2>code)
							break;
						cdf=cdf2;
						++token;
						
						freq=(int)freq2[1];
						cdf2=cdf+freq;
						if(cdf2>code)
							break;
						cdf=cdf2;
						++token;
						
						freq=(int)freq2[2];
						cdf2=cdf+freq;
						if(cdf2>code)
							break;
						cdf=cdf2;
						++token;
						
						freq=(int)freq2[3];
						cdf2=cdf+freq;
						if(cdf2>code)
							break;
						cdf=cdf2;
						++token;
					}
#else
					cdf=0;
					token=0;
					for(;;)
					{
						unsigned cdf2;

						freq=0;
						MIX4(freq, token);
						cdf2=cdf+freq;
						if(cdf2>code)
							break;
#ifdef _DEBUG
						if(token>=args->tlevels)
							LOG_ERROR("YXC %d %d %d  token %d/%d", ky, kx, kc, token, args->tlevels);
#endif
						cdf=cdf2;
						++token;
					}
#endif
					ac3_dec_update_NPOT(&ec, cdf, freq, den);
					sym=token;
					if(sym>=(1<<CONFIG_EXP))
					{
						int lsb, msb;

						sym-=1<<CONFIG_EXP;
						lsb=sym&((1<<CONFIG_LSB)-1);
						sym>>=CONFIG_LSB;
						msb=sym&((1<<CONFIG_MSB)-1);
						sym>>=CONFIG_MSB;
						nbits=sym+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB);
						bypass=ac3_dec_bypass(&ec, nbits);
						sym=1;
						sym<<=CONFIG_MSB;
						sym|=msb;
						sym<<=nbits;
						sym|=bypass;
						sym<<=CONFIG_LSB;
						sym|=lsb;
					}
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
					curr[kc2+1]=error;
					error+=pred;
					curr[kc2+0]=yuv[kc]=error;
				}
				curr[kc2+0]-=offset;
				{
#if 1
#if 1
					curr_hist00[token]+=8;
					curr_hist01[token]+=8;
					curr_hist10[token]+=8;
					curr_hist11[token]+=8;
					int s00=curr_hist00[args->tlevels]+=8;
					int s01=curr_hist01[args->tlevels]+=8;
					int s10=curr_hist10[args->tlevels]+=8;
					int s11=curr_hist11[args->tlevels]+=8;
#else
					ALIGN(32) long long inc[4];
				//	__m256i axmask=_mm256_set_epi64x(0, ~0ULL, 0, ~0ULL);
				//	__m256i aymask=_mm256_set_epi64x(0, 0, ~0ULL, ~0ULL);
				//	__m256i axoffset=_mm256_set_epi64x(0, 1LL<<MIXBITS|1, 0, 1LL<<MIXBITS|1);//a-b = ~b+(a+1)
				//	__m256i ayoffset=_mm256_set_epi64x(0, 0, 1LL<<MIXBITS|1, 1LL<<MIXBITS|1);
					ax=_mm256_xor_si256(ax, _mm256_set_epi64x(0, ~0ULL, 0, ~0ULL));
					ay=_mm256_xor_si256(ay, _mm256_set_epi64x(0, 0, ~0ULL, ~0ULL));
					ax=_mm256_add_epi64(ax, _mm256_set_epi64x(0, 1LL<<MIXBITS|1, 0, 1LL<<MIXBITS|1));//a-b = ~b+(a+1)
					ay=_mm256_add_epi64(ay, _mm256_set_epi64x(0, 0, 1LL<<MIXBITS|1, 1LL<<MIXBITS|1));
					ax=_mm256_mul_epi32(ax, ay);
					ax=_mm256_srli_epi32(ax, MIXBITS*2-5);
					_mm256_store_si256((__m256i*)inc, ax);
					curr_hist00[token]+=(int)inc[0];
					curr_hist01[token]+=(int)inc[1];
					curr_hist10[token]+=(int)inc[2];
					curr_hist11[token]+=(int)inc[3];
					int s00=curr_hist00[args->tlevels]+=(int)inc[0];
					int s01=curr_hist01[args->tlevels]+=(int)inc[1];
					int s10=curr_hist10[args->tlevels]+=(int)inc[2];
					int s11=curr_hist11[args->tlevels]+=(int)inc[3];
#endif
					if(s00>=10752)//4296	6144	10752	65536
					{
						int sum=0;
						for(int ks=0;ks<args->tlevels;++ks)
							sum+=curr_hist00[ks]>>=1;
						curr_hist00[args->tlevels]=sum;
					}
					if(s01>=10752)
					{
						int sum=0;
						for(int ks=0;ks<args->tlevels;++ks)
							sum+=curr_hist01[ks]>>=1;
						curr_hist01[args->tlevels]=sum;
					}
					if(s10>=10752)
					{
						int sum=0;
						for(int ks=0;ks<args->tlevels;++ks)
							sum+=curr_hist10[ks]>>=1;
						curr_hist10[args->tlevels]=sum;
					}
					if(s11>=10752)
					{
						int sum=0;
						for(int ks=0;ks<args->tlevels;++ks)
							sum+=curr_hist11[ks]>>=1;
						curr_hist11[args->tlevels]=sum;
					}
#else
					int inc;
					inc=((1<<MIXBITS)-alphax)*((1<<MIXBITS)-alphay)>>(MIXBITS+MIXBITS-5); curr_hist00[token]+=inc; curr_hist00[args->tlevels]+=inc;//optimize
					inc=(             alphax)*((1<<MIXBITS)-alphay)>>(MIXBITS+MIXBITS-5); curr_hist01[token]+=inc; curr_hist01[args->tlevels]+=inc;
					inc=((1<<MIXBITS)-alphax)*(             alphay)>>(MIXBITS+MIXBITS-5); curr_hist10[token]+=inc; curr_hist10[args->tlevels]+=inc;
					inc=(             alphax)*(             alphay)>>(MIXBITS+MIXBITS-5); curr_hist11[token]+=inc; curr_hist11[args->tlevels]+=inc;
					if(curr_hist00[args->tlevels]>=10752)//4296	6144	10752	65536
					{
						int sum=0;
						for(int ks=0;ks<args->tlevels;++ks)
							sum+=curr_hist00[ks]>>=1;
						//	sum+=curr_hist00[ks]=(curr_hist00[ks]+1)>>1;
						curr_hist00[args->tlevels]=sum;
					}
					if(curr_hist01[args->tlevels]>=10752)
					{
						int sum=0;
						for(int ks=0;ks<args->tlevels;++ks)
							sum+=curr_hist01[ks]>>=1;
						//	sum+=curr_hist01[ks]=(curr_hist01[ks]+1)>>1;
						curr_hist01[args->tlevels]=sum;
					}
					if(curr_hist10[args->tlevels]>=10752)
					{
						int sum=0;
						for(int ks=0;ks<args->tlevels;++ks)
							sum+=curr_hist10[ks]>>=1;
						//	sum+=curr_hist10[ks]=(curr_hist10[ks]+1)>>1;
						curr_hist10[args->tlevels]=sum;
					}
					if(curr_hist11[args->tlevels]>=10752)
					{
						int sum=0;
						for(int ks=0;ks<args->tlevels;++ks)
							sum+=curr_hist11[ks]>>=1;
						//	sum+=curr_hist11[ks]=(curr_hist11[ks]+1)>>1;
						curr_hist11[args->tlevels]=sum;
					}
#endif
				}
			}
			if(!args->fwd)
			{
				switch(bestrct)
				{
				case RCT_R_G_B:
				case RCT_R_G_BG:
				case RCT_R_G_BR:
				case RCT_R_GR_BR:
				case RCT_R_GR_BG:
				//case RCT_R_G_B2:
				//case RCT_R_GR_B2:
					args->dst[idx+0]=yuv[0]+128;
					args->dst[idx+1]=yuv[1]+128;
					args->dst[idx+2]=yuv[2]+128;
					break;
				case RCT_G_B_RG:
				case RCT_G_B_RB:
				case RCT_G_BG_RG:
				case RCT_G_BG_RB:
				//case RCT_G_B_R2:
				//case RCT_G_BG_R2:
					args->dst[idx+1]=yuv[0]+128;
					args->dst[idx+2]=yuv[1]+128;
					args->dst[idx+0]=yuv[2]+128;
					break;
				case RCT_B_R_GR:
				case RCT_B_R_GB:
				case RCT_B_RB_GB:
				case RCT_B_RB_GR:
				//case RCT_B_RB_G2:
					args->dst[idx+2]=yuv[0]+128;
					args->dst[idx+0]=yuv[1]+128;
					args->dst[idx+1]=yuv[2]+128;
					break;
				case RCT_G_RG_BR:
				//case RCT_G_RG_B2:
					args->dst[idx+1]=yuv[0]+128;
					args->dst[idx+0]=yuv[1]+128;
					args->dst[idx+2]=yuv[2]+128;
					break;
				case RCT_B_GB_RG:
				//case RCT_B_GB_R2:
					args->dst[idx+2]=yuv[0]+128;
					args->dst[idx+1]=yuv[1]+128;
					args->dst[idx+0]=yuv[2]+128;
					break;
				case RCT_R_BR_GB:
				//case RCT_R_B_G2:
				//case RCT_R_BR_G2:
					args->dst[idx+0]=yuv[0]+128;
					args->dst[idx+2]=yuv[1]+128;
					args->dst[idx+1]=yuv[2]+128;
					break;
				}
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
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	if(args->fwd)
		ac3_enc_flush(&ec);
}
int c01_codec(const char *srcfn, const char *dstfn)
{
	const int nch=3, depth=8;
	double t0;
	ArrayHandle src, dst;
	int headersize, printed;
	int iw, ih;
	const unsigned char *image, *imageend;
	unsigned char *image2;
	CodecID codec;
	int ncores=query_cpu_cores();
	int xblocks, yblocks, nblocks, nthreads, coffset;
	ptrdiff_t start, memusage, argssize;
	ThreadArgs *args;
	int test, fwd;
	int tlevels, clevels, histsize, statssize;
	double esize;
	int usize;
	
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
	{
		int nlevels=256;
		int token=0, bypass=0, nbits=0;

		quantize_pixel(nlevels, &token, &bypass, &nbits);
		tlevels=token+1;
		clevels=CLEVELS;
		statssize=(clevels*clevels*(tlevels+1)*nch+4LL)*(int)sizeof(int);//AVX2
	//	statssize=clevels*clevels*(tlevels+1)*nch*(int)sizeof(int);
		histsize=(int)sizeof(int[OCH_COUNT*PRED_COUNT<<8]);
		if(histsize<statssize)
			histsize=statssize;
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=image;
		arg->dst=fwd?0:dst->data+printed;
		arg->iw=iw;
		arg->ih=ih;
		arg->bufsize=sizeof(short[4*OCH_COUNT*2])*(BLOCKSIZE+16LL);//4 padded rows * OCH_COUNT * {pixels, wg_errors}
		arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));

		arg->histsize=histsize;
		arg->hist=(int*)malloc(histsize);

		//arg->statssize=statssize;
		//arg->stats=(unsigned*)malloc(statssize);
		if(!arg->pixels||!arg->hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=arg->bufsize;
		memusage+=arg->histsize;
		//memusage+=arg->statssize;
		
		arg->tlevels=tlevels;
		arg->clevels=clevels;
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

						kx=kt+kt2;
						ky=kx/xblocks;
						kx%=xblocks;
						if(nblocks<MAXPRINTEDBLOCKS)
						{
							//if(!(kt+kt2))
							//	printf("block,  nrows,  usize,     best  ->  actual,  (actual-best)\n");
							printf(
								"block %4d/%4d  XY %3d %3d  %4d*%4d:  %8d->%16lf->%8zd bytes (%+10.2lf)  %10.6lf%%  CR %10lf  %s %s %s %s\n",
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
								rct_names[arg->bestrct],
								pred_names[arg->predidx[0]],
								pred_names[arg->predidx[1]],
								pred_names[arg->predidx[2]]
							);
						}
						esize+=arg->bestsize;
#ifdef ABAC_PROFILE_SIZE
						csizes[0]+=arg->csizes[0];
						csizes[1]+=arg->csizes[1];
						csizes[2]+=arg->csizes[2];
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
				printf("Best %15.2lf (%+13.2lf) bytes\n", esize, csize-esize);
				printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
				printf("Mem usage: ");
				print_size((double)memusage, 8, 4, 0, 0);
				printf("\n");
			}
			printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
			if(!fwd)
				compare_bufs_8(image2, src->data+headersize, iw, ih, nch, nch, "C01", 0, 1);
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
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		_mm_free(arg->pixels);
		free(arg->hist);
		//free(arg->stats);
	}
	free(args);
	array_free(&src);
	array_free(&dst);
	//if(test)
	//	pause();
	return 0;
}