#include"codec.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT
//	#define PROFILE_SIZE


	#define ENCODE_ERROR
//	#define AC_SYMMETRIC
//	#define DISABLE_RCTPROBE
	#define AC3_PREC

#define CODECNAME "C02"
#include"entropy.h"

#define BLOCKSIZE 640
#define MAXPRINTEDBLOCKS 0//200

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
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		3, 3, 3,	3, 3, 3,	0, 0, 0)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		3, 3, 1,	3, 3, 3,	0, 0, 0)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		3, 3, 0,	3, 3, 3,	0, 0, 0)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		3, 3, 0,	3, 3, 3,	0, 0, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		3, 3, 1,	3, 3, 3,	0, 0, 0)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		3, 3, 1,	3, 3, 3,	0, 0, 0)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		3, 3, 0,	3, 3, 3,	0, 0, 0)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		3, 0, 0,	3, 3, 3,	0, 0, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		3, 0, 0,	3, 3, 3,	0, 0, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		3, 0, 0,	3, 3, 3,	0, 0, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(R_G_B2,	OCH_R,		OCH_G,		OCH_B2,		3, 3, 0,	3, 3, 1,	0, 0, 1)\
	RCT(R_B_G2,	OCH_R,		OCH_B,		OCH_G2,		3, 3, 0,	3, 3, 1,	0, 0, 1)\
	RCT(G_B_R2,	OCH_G,		OCH_B,		OCH_R2,		3, 3, 0,	3, 3, 1,	0, 0, 1)\
	RCT(R_GR_B2,	OCH_R,		OCH_GR,		OCH_B2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(R_BR_G2,	OCH_R,		OCH_BR,		OCH_G2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(G_BG_R2,	OCH_G,		OCH_BG,		OCH_R2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(G_RG_B2,	OCH_G,		OCH_RG,		OCH_B2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(B_RB_G2,	OCH_B,		OCH_RB,		OCH_G2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(B_GB_R2,	OCH_B,		OCH_GB,		OCH_R2,		3, 0, 0,	3, 3, 1,	0, 0, 1)
typedef enum _RCTIndex
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2, YSH2, USH2, VSH2) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][12]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2, YSH2, USH2, VSH2)\
	{YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2, YSH2, USH2, VSH2},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2, YSH2, USH2, VSH2) #LABEL,
	RCTLIST
#undef  RCT
};

#define PREDLIST\
	PRED(W)\
	PRED(CG)\
	PRED(AV5)\
	PRED(AV9)\
	PRED(AV12)
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
	NB_NNWW,	NB_NNW,		NB_NN,		NB_NNE,		NB_NNEE,
	NB_NWW,		NB_NW,		NB_N,		NB_NE,		NB_NEE,
	NB_WW,		NB_W,		NB_curr,
} NBIndex;
static const short av12_icoeffs[12]=
{
	 0x04,	 0x03,	-0x1F,	-0x26,	 0x00,
	 0x07,	-0x9E,	 0xDB,	 0x1E,	 0x13,
	-0x2A,	 0xF3,
};

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
		*nbits=lgv-(CONFIG_MSB+CONFIG_LSB);
		*bypass=val>>CONFIG_LSB&((1LL<<*nbits)-1);
	}
}
static int quantize_ctx(int x)
{
	int negmask=x>>31;
	x=abs(x)+1;
	x=FLOOR_LOG2(x);
	x^=negmask;
	x-=negmask;
	return x;

	//x=x<<1^x>>31;
	//if(x<1<<(6-3))
	//	return x;
	//int lgx=FLOOR_LOG2(x);
	//int mantissa=x-(1<<lgx);
	//int token=(1<<(6-3))+((lgx-(6-3))<<3|mantissa>>(lgx-3));
	//return token;
}
#define SCALE_BITS 16
static int squash(int x)//sigmoid(x) = 1/(1-exp(-x))		logit sum -> prob
{
#ifdef DISABLE_LOGMIX
	x>>=11;
	x+=1<<SCALE_BITS>>1;
	CLAMP2_32(x, x, 1, (1<<SCALE_BITS)-1);
#else
	static const int t[33]=//2^5 table elements, table amplitude 2^12
	{
		   1,    2,    3,    6,   10,   16,   27,   45,   73,  120,  194,
		 310,  488,  747, 1101, 1546, 2047, 2549, 2994, 3348, 3607, 3785,
		3901, 3975, 4022, 4050, 4068, 4079, 4085, 4089, 4092, 4093, 4094,
	};
	int w=x&((1<<(SCALE_BITS-5))-1);
	x=(x>>(SCALE_BITS-5))+16;
	if(x>31)
		return (1<<SCALE_BITS)-1;
	if(x<0)
		return 1;
	x=(t[x]*((1<<(SCALE_BITS-5))-w)+t[x+1]*w+64)>>(12-5);
#endif
	return x;
}
static int stretch(int x)//ln(x/(1-x))		probs -> logits
{
#ifndef DISABLE_LOGMIX
	static short t[4096];
	static int initialized=0;
	if(!initialized)
	{
		initialized=1;
		
		int pi=0;
		for(int k=-2047;k<=2047;++k)//invert squash()
		{
			int i=squash(k<<(SCALE_BITS-12))>>(SCALE_BITS-12);
			for(int j=pi;j<=i;++j)
				t[j]=k<<(SCALE_BITS-12);
			pi=i+1;
		}
		t[4095]=(1<<SCALE_BITS>>1)-1;
	}
	x=t[x>>(SCALE_BITS-12)];
#endif
	return x;
}
static void xorshift64(unsigned long long *state)
{
	unsigned long long x=*state;
	x^=x<<13;
	x^=x>>7;
	x^=x<<17;
	*state=x;
}
typedef struct _ThreadArgs
{
	const unsigned char *src;
	unsigned char *dst;
	int iw, ih;

	int fwd, test, loud, x1, x2, y1, y2;
	int bufsize, histsize;
	short *pixels;
	int *hist;

	DList list;
	const unsigned char *decstart, *decend;
	
	int stats0[3][1<<8];
	//int stats1[3][1<<8][1<<8];
	//unsigned char stats2[3][1<<8][1<<8][1<<8];
#ifdef PROFILE_SIZE
	double csizes[24];
#endif

	//aux
	int blockidx;
	double bestsize;
	int bestrct, predidx[3];
} ThreadArgs;
static void block_thread(void *param)
{
	const int nch=3, depth=8, half=128;
	ThreadArgs *args=(ThreadArgs*)param;
	AC2 ec;
	const unsigned char *image=args->fwd?args->src:args->dst;
	unsigned char bestrct=0, combination[6]={0}, predidx[4]={0};
	int ystride=args->iw*3;
	
	if(args->fwd)
	{
		int res;
#ifndef DISABLE_RCTPROBE
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		unsigned char predsel[OCH_COUNT]={0};
		__m256i av12_mcoeffs[12];

		for(int k=0;k<(int)_countof(av12_mcoeffs);++k)
			av12_mcoeffs[k]=_mm256_set1_epi16(av12_icoeffs[k]>>1);
		memset(args->hist, 0, args->histsize);
		for(int ky=args->y1+2;ky<args->y2;++ky)//analysis loop
		{
			int kx=args->x1+2;
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
			for(;kx<args->x2-5;kx+=5, ptr+=15)
			{
				__m256i
					nb0[13],//rgb
					nb1[13],//gbr
					nb2[13],//rgb - gbr
					nb3[13],//gbr - rgb
					nb4[13],//(gbr+brg)/2
					nb5[13];//rgb - (gbr+brg)/2
				__m256i vmin[4], vmax[4], pred;
				{
					__m128i nb8[13]=//8-bit
					{
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride-2*3+0)), half8),//NNWW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride-1*3+0)), half8),//NNW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride+0*3+0)), half8),//NN
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride+1*3+0)), half8),//NNE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride+2*3+0)), half8),//NNEE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride-2*3+0)), half8),//NWW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride-1*3+0)), half8),//NW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+0*3+0)), half8),//N
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+1*3+0)), half8),//NE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+2*3+0)), half8),//NEE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride-2*3+0)), half8),//WW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride-1*3+0)), half8),//W
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride+0*3+0)), half8),//curr
					};
					for(int k=0;k<13;++k)
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
				//W
				pred=_mm256_sub_epi16(nb0[NB_curr], nb0[NB_W]);
				UPDATE(
					PRED_W,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);

				pred=_mm256_add_epi16(nb2[NB_W], nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_W,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);

				pred=_mm256_add_epi16(nb3[NB_W], nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_W,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);

				pred=_mm256_add_epi16(nb5[NB_W], nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_W,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

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

				//AV5
				//		-5	5	1
				//	-1	8	[?]>>3
				pred=_mm256_sub_epi16(nb0[NB_N], nb0[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb0[NB_NE], nb0[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb0[NB_W]);
				//vmin=_mm256_min_epi16(N, W);
				//vmax=_mm256_max_epi16(N, W);
				vmin[0]=_mm256_min_epi16(vmin[0], nb0[NB_NE]);
				vmax[0]=_mm256_max_epi16(vmax[0], nb0[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_sub_epi16(nb2[NB_N], nb2[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb2[NB_NE], nb2[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb2[NB_W]);
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				vmin[1]=_mm256_min_epi16(vmin[1], nb2[NB_NE]);
				vmax[1]=_mm256_max_epi16(vmax[1], nb2[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_sub_epi16(nb3[NB_N], nb3[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb3[NB_NE], nb3[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb3[NB_W]);
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				vmin[2]=_mm256_min_epi16(vmin[2], nb3[NB_NE]);
				vmax[2]=_mm256_max_epi16(vmax[2], nb3[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_sub_epi16(nb5[NB_N], nb5[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb5[NB_NE], nb5[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb5[NB_W]);
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				vmin[3]=_mm256_min_epi16(vmin[1], nb5[NB_NE]);
				vmax[3]=_mm256_max_epi16(vmax[1], nb5[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//AV9
				//		1	-2	-1
				//	-1	-9	10	4
				//	-2	16	[?]>>4
				pred=_mm256_add_epi16(nb0[NB_N], _mm256_slli_epi16(nb0[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb0[NB_NN], nb0[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb0[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb0[NB_NW], 3), nb0[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb0[NB_NNW], _mm256_add_epi16(nb0[NB_NNE], nb0[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb0[NB_W], _mm256_srai_epi16(pred, 4));
				//vmin=_mm256_min_epi16(N, W);
				//vmax=_mm256_max_epi16(N, W);
				//vmin=_mm256_min_epi16(vmin, NE);
				//vmax=_mm256_max_epi16(vmax, NE);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(nb2[NB_N], _mm256_slli_epi16(nb2[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb2[NB_NN], nb2[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb2[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb2[NB_NW], 3), nb2[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb2[NB_NNW], _mm256_add_epi16(nb2[NB_NNE], nb2[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb2[NB_W], _mm256_srai_epi16(pred, 4));
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin=_mm256_min_epi16(vmin, NE3);
				//vmax=_mm256_max_epi16(vmax, NE3);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);
					
				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(nb3[NB_N], _mm256_slli_epi16(nb3[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb3[NB_NN], nb3[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb3[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb3[NB_NW], 3), nb3[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb3[NB_NNW], _mm256_add_epi16(nb3[NB_NNE], nb3[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb3[NB_W], _mm256_srai_epi16(pred, 4));
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin=_mm256_min_epi16(vmin, NE3);
				//vmax=_mm256_max_epi16(vmax, NE3);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);
					
				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_add_epi16(nb5[NB_N], _mm256_slli_epi16(nb5[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb5[NB_NN], nb5[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb5[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb5[NB_NW], 3), nb5[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb5[NB_NNW], _mm256_add_epi16(nb5[NB_NNE], nb5[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb5[NB_W], _mm256_srai_epi16(pred, 4));
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin=_mm256_min_epi16(vmin, NE3);
				//vmax=_mm256_max_epi16(vmax, NE3);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);
					
				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//AV12
				pred=_mm256_setzero_si256();
				for(int k=0;k<(int)_countof(av12_icoeffs);++k)
					pred=_mm256_add_epi16(pred, _mm256_mullo_epi16(av12_mcoeffs[k], nb0[k]));
				pred=_mm256_srai_epi16(pred, 7);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV12,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_setzero_si256();
				for(int k=0;k<(int)_countof(av12_icoeffs);++k)
					pred=_mm256_add_epi16(pred, _mm256_mullo_epi16(av12_mcoeffs[k], nb2[k]));
				pred=_mm256_srai_epi16(pred, 7);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);
					
				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV12,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_setzero_si256();
				for(int k=0;k<(int)_countof(av12_icoeffs);++k)
					pred=_mm256_add_epi16(pred, _mm256_mullo_epi16(av12_mcoeffs[k], nb3[k]));
				pred=_mm256_srai_epi16(pred, 7);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);
					
				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV12,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_setzero_si256();
				for(int k=0;k<(int)_countof(av12_icoeffs);++k)
					pred=_mm256_add_epi16(pred, _mm256_mullo_epi16(av12_mcoeffs[k], nb5[k]));
				pred=_mm256_srai_epi16(pred, 7);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);
					
				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV12,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);
			}
		}
		res=(args->x2-args->x1-3)/5*5*(args->y2-args->y1-2);
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
			printf("X %5d~%5d  Y %5d~%5d  best %12.2lf bytes  %s [YUV: %s %s %s]\n",
				args->x1, args->x2, args->y1, args->y2,
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
#else
		bestrct=RCT_G_RG_B2;
		predidx[0]=PRED_AV12;
		predidx[1]=PRED_AV12;
		predidx[2]=PRED_AV12;
#endif
		res=(args->x2-args->x1)*(args->y2-args->y1)*3;
		dlist_init(&args->list, 1, res, 0);
		ac2_enc_init(&ec, &args->list);
		ac2_enc_bypass_NPOT(&ec, bestrct, RCT_COUNT);
		ac2_enc_bypass_NPOT(&ec, predidx[0], PRED_COUNT);
		ac2_enc_bypass_NPOT(&ec, predidx[1], PRED_COUNT);
		ac2_enc_bypass_NPOT(&ec, predidx[2], PRED_COUNT);
	}
	else
	{
		ac2_dec_init(&ec, args->decstart, args->decend);
		bestrct=ac2_dec_bypass_NPOT(&ec, RCT_COUNT);
		predidx[0]=ac2_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[1]=ac2_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[2]=ac2_dec_bypass_NPOT(&ec, PRED_COUNT);
	}
	
	unsigned long long prngstate=0x3243F6A8885A308D;//pi
	memset(args->stats0, 0, sizeof(args->stats0));
	//memset(args->stats1, 0, sizeof(args->stats1));
	//memset(args->stats2, 0, sizeof(args->stats2));
#ifdef PROFILE_SIZE
	memset(args->csizes, 0, sizeof(args->csizes));
#endif

	memset(args->pixels, 0, args->bufsize);
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
		int pred=0, sym=0;
		const unsigned char *combination=rct_combinations[bestrct];
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int idx=nch*(args->iw*ky+kx);
			short
				*NNN	=rows[3]+0*4*2,
				*NNWW	=rows[2]-2*4*2,
				*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NNEE	=rows[2]+2*4*2,
				*NNEEE	=rows[2]+3*4*2,
				*NWW	=rows[1]-2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*NEE	=rows[1]+2*4*2,
				*NEEE	=rows[1]+3*4*2,
				*WWWW	=rows[0]-4*4*2,
				*WWW	=rows[0]-3*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
			if(ky<=args->y1+2)
			{
				if(ky<=args->y1+1)
				{
					if(ky==args->y1)
						NEEE=NEE=NE=NWW=NW=N=W;
					NNWW=NWW;
					NNW=NW;
					NN=N;
					NNE=NE;
					NNEE=NEE;
					NNEEE=NEEE;
				}
				NNN=NN;
			}
			if(kx<=args->x1+3)
			{
				if(kx<=args->x1+2)
				{
					if(kx<=args->x1+1)
					{
						if(kx<=args->x1)
							NW=W=N;
						WW=W;
						NWW=NW;
					}
					WWW=WW;
				}
				WWWW=WWW;
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
					NNEE=NNE;
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
				case RCT_R_G_B2:
				case RCT_R_GR_B2:
					yuv[0]=args->src[idx+0]-128;
					yuv[1]=args->src[idx+1]-128;
					yuv[2]=args->src[idx+2]-128;
					break;
				case RCT_G_B_RG:
				case RCT_G_B_RB:
				case RCT_G_BG_RG:
				case RCT_G_BG_RB:
				case RCT_G_B_R2:
				case RCT_G_BG_R2:
					yuv[0]=args->src[idx+1]-128;
					yuv[1]=args->src[idx+2]-128;
					yuv[2]=args->src[idx+0]-128;
					break;
				case RCT_B_R_GR:
				case RCT_B_R_GB:
				case RCT_B_RB_GB:
				case RCT_B_RB_GR:
				case RCT_B_RB_G2:
					yuv[0]=args->src[idx+2]-128;
					yuv[1]=args->src[idx+0]-128;
					yuv[2]=args->src[idx+1]-128;
					break;
				case RCT_G_RG_BR:
				case RCT_G_RG_B2:
					yuv[0]=args->src[idx+1]-128;
					yuv[1]=args->src[idx+0]-128;
					yuv[2]=args->src[idx+2]-128;
					break;
				case RCT_B_GB_RG:
				case RCT_B_GB_R2:
					yuv[0]=args->src[idx+2]-128;
					yuv[1]=args->src[idx+1]-128;
					yuv[2]=args->src[idx+0]-128;
					break;
				case RCT_R_BR_GB:
				case RCT_R_B_G2:
				case RCT_R_BR_G2:
					yuv[0]=args->src[idx+0]-128;
					yuv[1]=args->src[idx+2]-128;
					yuv[2]=args->src[idx+1]-128;
					break;
				}
			}
			//if(ky==289&&kx==544)
			//if(ky==128&&kx==128)
			//	printf("");
			for(int kc=0;kc<nch;++kc)
			{
				int kc2=kc<<1;
				int offset=(yuv[combination[kc+3]]+yuv[combination[kc+6]])>>combination[kc+9];

				switch(predidx[kc])
				{
				case PRED_W:
					pred=W[kc2];
					break;
				case PRED_CG:
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					break;
				case PRED_AV5:
					CLAMP3_32(pred, W[kc2]+((5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])>>3), N[kc2], W[kc2], NE[kc2]);
					break;
				case PRED_AV9:
					CLAMP3_32(pred,
						W[kc2]+((10*N[kc2]-9*NW[kc2]+4*NE[kc2]-2*(NN[kc2]+WW[kc2])+NNW[kc2]-(NNE[kc2]+NWW[kc2]))>>4),
						N[kc2],
						W[kc2],
						NE[kc2]
					);
					break;
				case PRED_AV12:
					pred=(
						av12_icoeffs[ 0]*NNWW[kc2]+
						av12_icoeffs[ 1]*NNW[kc2]+
						av12_icoeffs[ 2]*NN[kc2]+
						av12_icoeffs[ 3]*NNE[kc2]+
						av12_icoeffs[ 4]*NNEE[kc2]+
						av12_icoeffs[ 5]*NWW[kc2]+
						av12_icoeffs[ 6]*NW[kc2]+
						av12_icoeffs[ 7]*N[kc2]+
						av12_icoeffs[ 8]*NE[kc2]+
						av12_icoeffs[ 9]*NEE[kc2]+
						av12_icoeffs[10]*WW[kc2]+
						av12_icoeffs[11]*W[kc2]
					)>>8;
					CLAMP3_32(pred, pred, N[kc2], W[kc2], NE[kc2]);
					break;
				}
				pred+=offset;
				CLAMP2_32(pred, pred, -128, 127);

				int *curr_stats0=args->stats0[kc];
				//int *curr_stats1=args->stats1[kc][((N[kc2+0]+W[kc2+0])/2+128)&255];
				//unsigned char *curr_stats2=args->stats2[kc][(WW[kc2+0]+128)&255][(W[kc2+0]+128)&255];
				if(args->fwd)
				{
					//if(ky==3&&kx==128&&!kc)//
					//	printf("");
					curr[kc2+0]=yuv[kc];
					curr[kc2+1]=yuv[kc]-pred;
#ifdef ENCODE_ERROR
					sym=curr[kc2+1];
					int upred=half-abs(pred), aval=abs(sym);
					if(aval<=upred)
					{
#ifdef ENABLE_BIASCORR
						{
							int negmask=ibias_corr>>31&-(sym!=-128);//sign is flipped if SSE correction was negative, to skew the histogram
							sym^=negmask;
							sym-=negmask;
						}
#endif
						sym=sym<<1^sym>>31;//pack sign
					}
					else
						sym=aval+upred;//error sign is known	sym=2^n=256 when pixel=-2^(n-1)=-128 and pred>0
#else
					sym=yuv[kc];
#endif
				}
				else
					sym=0;
				for(int kb=8, tidx=0;kb>=0;--kb)
				{
					long long probs[]=
					{
						curr_stats0[tidx],
					//	curr_stats1[tidx],
					};
					//long long p1=(probs[0]+3LL*probs[1]+2)>>2;
					//long long p1=probs[1];
					long long p1=probs[0];
					int bit;
					p1+=1ULL<<AC2_PROB_BITS>>1;
					//CLAMP2_32(p1, (int)p1, 1, (int)((1ULL<<AC2_PROB_BITS)-1));
					if(args->fwd)
					{
						bit=sym>>kb&1;
						ac2_enc_bin(&ec, (unsigned)p1, bit);
#ifdef PROFILE_SIZE
						if(kb<8)
							args->csizes[kc<<3|kb]-=log2((double)(bit?p1:0x80000000-p1)/0x80000000);
#endif
					}
					else
					{
						bit=ac2_dec_bin(&ec, (unsigned)p1);
						sym|=bit<<kb;
					}
					//if(!kc&&kb==7)
					//	printf("");
					//long long perr=(long long)(bit<<AC2_PROB_BITS)-p1;
					long long collapse=(long long)(bit<<AC2_PROB_BITS)-(int)(1ULL<<AC2_PROB_BITS>>1);
					probs[0]+=((collapse-probs[0])*17+(1<<11>>1))>>11;//30 bit
					//probs[0]+=((collapse-probs[0])*17+(1<<22>>1))>>22;//31 bit
					//probs[1]+=(collapse-probs[1]+(1<<5>>1))>>5;
					const int pmin=-((int)(1ULL<<AC2_PROB_BITS>>1)-1), pmax=(int)(1ULL<<AC2_PROB_BITS>>1)-1;
					CLAMP2(probs[0], pmin, pmax);
					//CLAMP2(probs[1], -((1LL<<AC2_PROB_BITS>>1)-1), (1LL<<AC2_PROB_BITS>>1)-1);
					//long long probs0[]=
					//{
					//	probs[0],
					//	probs[1],
					//};
					//if(probs[0]/(1<<31)||probs[1]/(1<<31))
					//	LOG_ERROR("");
					//probs[0]+=(int)(((long long)(bit<<AC2_PROB_BITS)-probs[0]+(1<<7>>1))>>7);
					//probs[1]+=(int)(((long long)(bit<<AC2_PROB_BITS)-probs[1]+(1<<5>>1))>>5);
					//if(probs[0]/(1<<31)||probs[1]/(1<<31))
					//	LOG_ERROR("");
					curr_stats0[tidx]=(int)probs[0];
					//curr_stats1[tidx]=(int)probs[1];
					if(kb==8&&bit)
						break;
					tidx+=tidx+!bit;
				}
				if(!args->fwd)
				{
					//if(ky==3&&kx==128&&!kc)//
					//	printf("");
#ifdef ENCODE_ERROR
					int upred=half-abs(pred), negmask=0;
					if(sym<=(upred<<1))
					{
						sym=sym>>1^-(sym&1);
#ifdef ENABLE_BIASCORR
						negmask=ibias_corr>>31&-(sym!=-128);
#endif
					}
					else
					{
						sym=sym-upred;
						negmask=(-pred)>>31;
					}
					sym^=negmask;
					sym-=negmask;
					sym+=pred;
#else
					sym=sym<<(32-8)>>(32-8);
#endif
					curr[kc2+0]=yuv[kc]=sym;
					curr[kc2+1]=sym-pred;
				}
				curr[kc2+0]-=offset;
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
				case RCT_R_G_B2:
				case RCT_R_GR_B2:
					args->dst[idx+0]=yuv[0]+128;
					args->dst[idx+1]=yuv[1]+128;
					args->dst[idx+2]=yuv[2]+128;
					break;
				case RCT_G_B_RG:
				case RCT_G_B_RB:
				case RCT_G_BG_RG:
				case RCT_G_BG_RB:
				case RCT_G_B_R2:
				case RCT_G_BG_R2:
					args->dst[idx+1]=yuv[0]+128;
					args->dst[idx+2]=yuv[1]+128;
					args->dst[idx+0]=yuv[2]+128;
					break;
				case RCT_B_R_GR:
				case RCT_B_R_GB:
				case RCT_B_RB_GB:
				case RCT_B_RB_GR:
				case RCT_B_RB_G2:
					args->dst[idx+2]=yuv[0]+128;
					args->dst[idx+0]=yuv[1]+128;
					args->dst[idx+1]=yuv[2]+128;
					break;
				case RCT_G_RG_BR:
				case RCT_G_RG_B2:
					args->dst[idx+1]=yuv[0]+128;
					args->dst[idx+0]=yuv[1]+128;
					args->dst[idx+2]=yuv[2]+128;
					break;
				case RCT_B_GB_RG:
				case RCT_B_GB_R2:
					args->dst[idx+2]=yuv[0]+128;
					args->dst[idx+1]=yuv[1]+128;
					args->dst[idx+0]=yuv[2]+128;
					break;
				case RCT_R_BR_GB:
				case RCT_R_B_G2:
				case RCT_R_BR_G2:
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
		ac2_enc_flush(&ec);
}
int c02_codec(const char *srcfn, const char *dstfn)
{
	//test();//

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
#ifdef ENABLE_MIX4
	int clevels;
#endif
	int histsize;
	double esize;
	int usize;
#ifdef ABAC_PROFILESIZE
	double abac_csizes[ABAC_TOKEN_BITS*3]={0};
#endif
#ifdef ENABLE_ABAC
	stretch(0);//for thread safety
#endif
	
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
	histsize=(int)sizeof(int[OCH_COUNT*PRED_COUNT<<8]);
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
#ifdef ENABLE_ABAC
		arg->statssize=statssize;
#endif
		//arg->stats=(unsigned*)malloc(statssize);
		if(!arg->pixels||!arg->hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=arg->bufsize;
		memusage+=histsize;
		
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
#ifdef PROFILE_SIZE
						double usize=(double)blocksize/nch;
						for(int kp=0;kp<24;++kp)
						{
							double csize=arg->csizes[kp]/8;
							printf("%2d  %12.2lf/%10.2lf  %8.4lf%%  %12.6lf\n", kp, csize, usize, 100.*csize/usize, usize/csize);
						}
#endif
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
								arg->list.nobj,
								arg->list.nobj-arg->bestsize,
								100.*arg->list.nobj/blocksize,
								(double)blocksize/arg->list.nobj,
								rct_names[arg->bestrct],
								pred_names[arg->predidx[0]],
								pred_names[arg->predidx[1]],
								pred_names[arg->predidx[2]]
							);
						}
						esize+=arg->bestsize;
#ifdef ABAC_PROFILESIZE
						for(int k=0;k<ABAC_TOKEN_BITS*3;++k)
							abac_csizes[k]+=arg->abac_csizes[k];
#endif
					}
					memcpy(dst->data+start+sizeof(int)*((ptrdiff_t)kt+kt2), &arg->list.nobj, sizeof(int));
					dlist_appendtoarray(&arg->list, &dst);
					dlist_clear(&arg->list);
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
//#ifdef DISABLE_MT
//				printf("overflows %d/%lld = %lf%%\n", count_overflow, count_total, 100.*count_overflow/count_total);
//#endif
#ifdef ABAC_PROFILESIZE
				double upsize=(double)iw*ih/8;
				for(int kb=0;kb<ABAC_TOKEN_BITS*3;++kb)
				{
					double cpsize=abac_csizes[kb]/8;
					printf("C%d B%d  %10.2lf/%10.2lf  %8.4lf%%  CR %10.6lf\n",
						kb/ABAC_TOKEN_BITS,
						kb%ABAC_TOKEN_BITS,
						cpsize,
						upsize,
						100.*cpsize/upsize,
						upsize/cpsize
					);
				}
#endif
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