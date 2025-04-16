#include"codec.h"
#include"util.h"
#include<stdlib.h>
#include<math.h>//log2
static const char file[]=__FILE__;


	#define BYPASS_ON_INFLATION
//	#define ENABLE_FILEGUARD	//makes using scripts harder

//	#define ESTIMATE_SIZE

//	#define EXPERIMENTAL1
//	#define EXPERIMENTAL2
//	#define ENABLE_O1
//	#define ENABLE_CTR


#define NCTX 4
#define LR_PROB 20
#define LR_MIXER 21
#define SPLITLEVELS 3

#include"entropy.h"
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
static int calc_centroid(int *hist, int x1, int x2, int margin)
{
	long long num=0;
	int den=0;
	for(int k=x1;k<x2;++k)
	{
		num+=(long long)hist[k]*k;
		den+=hist[k];
	}
	int split;
	if(den)
		split=(int)(num/den);
	else
		split=128;
	if(split<x1+margin)
		split=x1+margin;
	if(split>x2-margin)
		split=x2-margin;
	return split;
}
static int decorrelate(unsigned char *buffer, int iw, int ih, unsigned char *ret_params, unsigned char *splits)
{
	int ystride=iw*3;
	double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
	unsigned char predsel[OCH_COUNT]={0};
	int res=(iw-3)/5*5*(ih-2);
	__m256i av12_mcoeffs[12];
	int histsize=(int)sizeof(int[OCH_COUNT*PRED_COUNT<<8]);
	int *hist=(int*)malloc(histsize);
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(hist, 0, histsize);
	for(int k=0;k<(int)_countof(av12_mcoeffs);++k)
		av12_mcoeffs[k]=_mm256_set1_epi16(av12_icoeffs[k]>>1);
	for(int ky=2;ky<ih;++ky)//analysis loop
	{
		int kx=2;
		const unsigned char *ptr=buffer+(iw*ky+kx)*3;

		__m256i amin=_mm256_set1_epi16(-128);
		__m256i amax=_mm256_set1_epi16(127);
		__m256i mask=_mm256_set1_epi16(255);
		__m128i half8=_mm_set1_epi8(-128);
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
		for(;kx<iw-5;kx+=5, ptr+=15)
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
//	pred=_mm256_slli_epi16(pred, 8);
//	pred=_mm256_xor_si256(_mm256_slli_epi16(pred, 7), _mm256_slli_epi16(pred, 15));
#define UPDATE(PREDIDX, IDX0, IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7, IDX8, IDX9, IDXA, IDXB, IDXC, IDXD, IDXE)\
do\
{\
	pred=_mm256_sub_epi16(pred, amin);\
	pred=_mm256_and_si256(pred, mask);\
	_mm256_store_si256((__m256i*)result, pred);\
	++hist[(IDX0*PRED_COUNT+PREDIDX)<<8|result[0x0]];\
	++hist[(IDX1*PRED_COUNT+PREDIDX)<<8|result[0x1]];\
	++hist[(IDX2*PRED_COUNT+PREDIDX)<<8|result[0x2]];\
	++hist[(IDX3*PRED_COUNT+PREDIDX)<<8|result[0x3]];\
	++hist[(IDX4*PRED_COUNT+PREDIDX)<<8|result[0x4]];\
	++hist[(IDX5*PRED_COUNT+PREDIDX)<<8|result[0x5]];\
	++hist[(IDX6*PRED_COUNT+PREDIDX)<<8|result[0x6]];\
	++hist[(IDX7*PRED_COUNT+PREDIDX)<<8|result[0x7]];\
	++hist[(IDX8*PRED_COUNT+PREDIDX)<<8|result[0x8]];\
	++hist[(IDX9*PRED_COUNT+PREDIDX)<<8|result[0x9]];\
	++hist[(IDXA*PRED_COUNT+PREDIDX)<<8|result[0xA]];\
	++hist[(IDXB*PRED_COUNT+PREDIDX)<<8|result[0xB]];\
	++hist[(IDXC*PRED_COUNT+PREDIDX)<<8|result[0xC]];\
	++hist[(IDXD*PRED_COUNT+PREDIDX)<<8|result[0xD]];\
	++hist[(IDXE*PRED_COUNT+PREDIDX)<<8|result[0xE]];\
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
	for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)
	{
		int *curr_hist=hist+((size_t)kc<<8);
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
	unsigned char bestrct=0;
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
	ret_params[0]=bestrct;//best RCT	5 bit
	ret_params[1]=predsel[rct_combinations[bestrct][0]];//pred idx	3 bit each
	ret_params[2]=predsel[rct_combinations[bestrct][1]];
	ret_params[3]=predsel[rct_combinations[bestrct][2]];
	for(int kc=0;kc<3;++kc)
	{
		unsigned char *curr_splits=splits+((size_t)kc<<SPLITLEVELS);
		unsigned char hidx=rct_combinations[bestrct][kc];
		int *curr_hist=hist+(((size_t)hidx*PRED_COUNT+predsel[hidx])<<8);

		curr_splits[1]=calc_centroid(curr_hist, 0, 256, 3);

		curr_splits[2]=calc_centroid(curr_hist, 0, curr_splits[1], 2);
		curr_splits[3]=calc_centroid(curr_hist, curr_splits[1], 256, 2);

		curr_splits[4]=calc_centroid(curr_hist, 0, curr_splits[2], 1);
		curr_splits[5]=calc_centroid(curr_hist, curr_splits[2], curr_splits[1], 1);
		curr_splits[6]=calc_centroid(curr_hist, curr_splits[1], curr_splits[3], 1);
		curr_splits[7]=calc_centroid(curr_hist, curr_splits[3], 256, 1);
	}
	free(hist);
	return 0;
}
int c08_codec(const char *srcfn, const char *dstfn, int nthreads0)
{
#ifdef ESTIMATE_SIZE
	double t=time_sec();
	double csizes[3]={0};
	size_t nqueries=0;
#endif
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec C08 requires both source and destination filenames");
		return 1;
	}
	ptrdiff_t dstsize;
#ifdef ENABLE_FILEGUARD
	dstsize=get_filesize(dstfn);
	if(dstsize>=0)
	{
		LOG_ERROR("Destination file already exists");
		return 1;
	}
#endif
	FILE *fsrc=fopen(srcfn, "rb");
	int tag=0;
	size_t nread=fread(&tag, 1, 2, fsrc), nwritten=0;
	if(nread!=2)
	{
		LOG_ERROR("File is empty");
		return 1;
	}
	int fwd=tag==('P'|'6'<<8);
	int bypass=tag==('B'|'P'<<8);
	if(!fwd&&tag!=('C'|'H'<<8)&&!bypass)
	{
		LOG_ERROR("Unsupported source file");
		return 1;
	}
	int iw=0, ih=0;
	unsigned char dparams[4]={0}, splits[3][1<<SPLITLEVELS]={0};
	int usize=0;
	if(fwd)
	{
		int temp=fgetc(fsrc);
		if(temp!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
		nread=fscanf(fsrc, "%d %d", &iw, &ih);
		if(nread!=2)
		{
			LOG_ERROR("Unsupported PPM file");
			return 1;
		}
		int vmax=0;
		nread=fscanf(fsrc, "%d", &vmax);
		if(nread!=1||vmax!=255)
		{
			LOG_ERROR("Unsupported PPM file");
			return 1;
		}
		temp=fgetc(fsrc);
		if(temp!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
	}
	else
	{
		nread=fread(&iw, 1, 4, fsrc);
		nread+=fread(&ih, 1, 4, fsrc);
		if(nread!=4LL+4)
		{
			LOG_ERROR("Unsupported archive");
			return 1;
		}
		if(!bypass)
		{
			nread=fread(dparams, 1, 4, fsrc);
			if(nread!=4||dparams[0]>=RCT_COUNT||dparams[1]>=PRED_COUNT||dparams[2]>=PRED_COUNT||dparams[3]>=PRED_COUNT)
			{
				LOG_ERROR("Unsupported archive");
				return 1;
			}
			nread =fread(splits[0]+1, 1, (1LL<<SPLITLEVELS)-1, fsrc);
			nread+=fread(splits[1]+1, 1, (1LL<<SPLITLEVELS)-1, fsrc);
			nread+=fread(splits[2]+1, 1, (1LL<<SPLITLEVELS)-1, fsrc);
			if(nread!=((1LL<<SPLITLEVELS)-1)*3)
			{
				LOG_ERROR("Unsupported archive");
				return 1;
			}
		}
	}
	usize=iw*ih*3;
	unsigned char *buffer=(unsigned char*)malloc(usize);
	if(!buffer)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	FILE *fdst=fopen(dstfn, "wb");
	AC5 ec;
	nwritten=0;
	if(fwd)
	{
		nread=fread(buffer, 1, usize, fsrc);
		if(nread!=usize)
		{
			LOG_ERROR("Corrupt PPM");
			return 1;
		}
		decorrelate(buffer, iw, ih, dparams, (unsigned char*)splits);

		nwritten+=fwrite("CH", 1, 2, fdst);
		nwritten+=fwrite(&iw, 1, 4, fdst);
		nwritten+=fwrite(&ih, 1, 4, fdst);
		nwritten+=fwrite(dparams, 1, 4, fdst);
		nwritten+=fwrite(splits[0]+1, 1, (1LL<<SPLITLEVELS)-1, fdst);
		nwritten+=fwrite(splits[1]+1, 1, (1LL<<SPLITLEVELS)-1, fdst);
		nwritten+=fwrite(splits[2]+1, 1, (1LL<<SPLITLEVELS)-1, fdst);
		ac5_enc_init(&ec, fdst);
	}
	else
	{
		nwritten+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		if(bypass)
		{
			printf("Bypass decode\n");
			nread=fread(buffer, 1, usize, fsrc);
			if(nread!=usize)
				printf("Error  read %d/%d bytes\n", (int)nread, usize);
			fwrite(buffer, 1, usize, fdst);
			fclose(fsrc);
			fclose(fdst);
			free(buffer);
			return 0;
		}
		memset(buffer, 0, usize);
		ac5_dec_init(&ec, fsrc);
	}
	static const unsigned char permutations[RCT_COUNT][3]=
	{
		{0, 1, 2},//R_G_B
		{0, 1, 2},//R_G_BG
		{0, 1, 2},//R_G_BR
		{1, 2, 0},//G_B_RG
		{1, 2, 0},//G_B_RB
		{2, 0, 1},//B_R_GR
		{2, 0, 1},//B_R_GB
		{1, 2, 0},//G_BG_RG
		{1, 2, 0},//G_BG_RB
		{1, 0, 2},//G_RG_BR
		{2, 0, 1},//B_RB_GB
		{2, 0, 1},//B_RB_GR
		{2, 1, 0},//B_GB_RG
		{0, 1, 2},//R_GR_BR
		{0, 1, 2},//R_GR_BG
		{0, 2, 1},//R_BR_GB
		{0, 1, 2},//R_G_B2
		{0, 2, 1},//R_B_G2
		{1, 2, 0},//G_B_R2
		{0, 1, 2},//R_GR_B2
		{0, 2, 1},//R_BR_G2
		{1, 2, 0},//G_BG_R2
		{1, 0, 2},//G_RG_B2
		{2, 0, 1},//B_RB_G2
		{2, 1, 0},//B_GB_R2
	};
	const unsigned char *permutation=permutations[dparams[0]];
	const unsigned char *combination=rct_combinations[dparams[0]];
	int psize=(iw+16LL)*sizeof(short[4*4]);//4 padded rows * 4 channels max
	short *pixels=(short*)malloc(psize);
#ifdef ENABLE_O1
	int stats1size=(int)sizeof(short[3*NCTX*256*256]);
	short *stats1=(short*)malloc(stats1size);
#ifdef EXPERIMENTAL1
	int esign1size=(int)sizeof(short[3*NCTX*256]);
	short *esign1=(short*)malloc(esign1size);
#endif
#endif
	if(!pixels
#ifdef ENABLE_O1
		||!stats1
#ifdef EXPERIMENTAL1
		||!esign1
#endif
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
#ifdef ENABLE_O1
	memset(stats1, 0, stats1size);
#ifdef EXPERIMENTAL1
	memset(esign1, 0, esign1size);
#endif
	short mixer1[3][256][NCTX]={0};
	FILLMEM((short*)mixer1, 0x4000, sizeof(mixer1), sizeof(short));
#else
//#ifdef EXPERIMENTAL1
//	unsigned short esign0[]={0x8000, 0x8000, 0x8000};
//#endif
	unsigned short stats0[3][256]={0};
	FILLMEM((unsigned short*)stats0, 0x8000, sizeof(stats0), sizeof(short));
#endif
	memset(pixels, 0, psize);
	for(int ky=0, idx=0;ky<ih;++ky)//codec loop
	{
		ALIGN(32) short *rows[]=
		{
			pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4,
			pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4,
		};
		int yuv[4]={0};
		for(int kx=0;kx<iw;++kx, idx+=3)
		{
			if(fwd)
			{
				yuv[0]=buffer[idx+permutation[0]];
				yuv[1]=buffer[idx+permutation[1]];
				yuv[2]=buffer[idx+permutation[2]];
			}
			for(int kc=0;kc<3;++kc)
			{
				int pred=0;
				int offset=(yuv[combination[kc+3]]+yuv[combination[kc+6]])>>combination[kc+9];
				short
					NNN	=rows[3][kc+0*4],
					NNWW	=rows[2][kc-2*4],
					NNW	=rows[2][kc-1*4],
					NN	=rows[2][kc+0*4],
					NNE	=rows[2][kc+1*4],
					NNEE	=rows[2][kc+2*4],
					NNEEE	=rows[2][kc+3*4],
					NWW	=rows[1][kc-2*4],
					NW	=rows[1][kc-1*4],
					N	=rows[1][kc+0*4],
					NE	=rows[1][kc+1*4],
					NEE	=rows[1][kc+2*4],
					NEEE	=rows[1][kc+3*4],
					WWWW	=rows[0][kc-4*4],
					WWW	=rows[0][kc-3*4],
					WW	=rows[0][kc-2*4],
					W	=rows[0][kc-1*4],
					*curr	=rows[0]+kc;
				(void)NNN;
				(void)NNEEE;
				(void)WWWW;
				if(ky<=2)
				{
					if(ky<=1)
					{
						if(ky==0)
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
				if(kx<=3)
				{
					if(kx<=2)
					{
						if(kx<=1)
						{
							if(kx<=0)
								NW=W=N;
							WW=W;
							NWW=NW;
						}
						WWW=WW;
					}
					WWWW=WWW;
				}
				if(kx>=iw-3)
				{
					if(kx>=iw-2)
					{
						if(kx>=iw-1)
						{
							NNE=NN;
							NE=N;
						}
						NNEE=NNE;
						NEE=NE;
					}
					NEEE=NEE;
				}
				switch(dparams[kc+1])
				{
				case PRED_W:
					pred=W;
					break;
				case PRED_CG:
					MEDIAN3_32(pred, N, W, N+W-NW);
					break;
				case PRED_AV5:
					CLAMP3_32(pred, W+((5*(N-NW)+NE-WW)>>3), N, W, NE);
					break;
				case PRED_AV9:
					CLAMP3_32(pred, W+((10*N-9*NW+4*NE-2*(NN+WW)+NNW-(NNE+NWW))>>4), N, W, NE);
					break;
				case PRED_AV12:
					pred=(
						av12_icoeffs[ 0]*NNWW+
						av12_icoeffs[ 1]*NNW+
						av12_icoeffs[ 2]*NN+
						av12_icoeffs[ 3]*NNE+
						av12_icoeffs[ 4]*NNEE+
						av12_icoeffs[ 5]*NWW+
						av12_icoeffs[ 6]*NW+
						av12_icoeffs[ 7]*N+
						av12_icoeffs[ 8]*NE+
						av12_icoeffs[ 9]*NEE+
						av12_icoeffs[10]*WW+
						av12_icoeffs[11]*W
					)>>8;
					CLAMP3_32(pred, pred, N, W, NE);
					break;
				}
				pred+=offset;
				CLAMP2_32(pred, pred, 0, 255);
				int error=0;
				if(fwd)
				{
					*curr=yuv[kc];
					error=yuv[kc]-pred;
					error+=128;
					error&=255;
					//error<<=24;
					//error>>=24;
				}
#ifdef EXPERIMENTAL1
				ALIGN(16) int ctx[]=
				{
					N+W-NW,
					N,
					W,
					(N+W)>>1,
				};
				{
					__m128i mc=_mm_load_si128((__m128i*)ctx);
					mc=_mm_min_epi32(mc, _mm_set1_epi32(255));
					mc=_mm_max_epi32(mc, _mm_setzero_si128());
					_mm_store_si128((__m128i*)ctx, mc);
				}
				short *curr_esign[]=
				{
					esign1+(((size_t)kc<<8|ctx[0])<<2|0),
					esign1+(((size_t)kc<<8|ctx[1])<<2|1),
					esign1+(((size_t)kc<<8|ctx[2])<<2|2),
					esign1+(((size_t)kc<<8|ctx[3])<<2|3),
				};
				short *curr_stats[]=
				{
					stats1+((((size_t)kc<<8|ctx[0])<<2|0)<<8),
					stats1+((((size_t)kc<<8|ctx[1])<<2|1)<<8),
					stats1+((((size_t)kc<<8|ctx[2])<<2|2)<<8),
					stats1+((((size_t)kc<<8|ctx[3])<<2|3)<<8),
				};
				short *curr_mixer=mixer1[kc][pred];
				short probs[4];
				int temp, p1, proberror;

			//	unsigned short *curr_stats, p1;
				int bit;
				
				probs[0]=curr_esign[0][0];
				probs[1]=curr_esign[1][0];
				probs[2]=curr_esign[2][0];
				probs[3]=curr_esign[3][0];
				p1=(int)((
					+(long long)curr_mixer[0]*probs[0]
					+(long long)curr_mixer[1]*probs[1]
					+(long long)curr_mixer[2]*probs[2]
					+(long long)curr_mixer[3]*probs[3]
				)>>19);
				p1+=0x8000;
				CLAMP2_32(p1, p1, 1, 0xFFFF);
			//	p1=esign1[kc];
				if(fwd)
				{
					bit=error>=128;
					ac5_enc_bin(&ec, p1, bit);
#ifdef ESTIMATE_SIZE
					csizes[kc]-=log2((double)(bit?p1:0x10000-p1)/0x10000);
					++nqueries;
#endif
				}
				else
					bit=ac5_dec_bin(&ec, p1);
				proberror=(bit<<16)-p1;
				temp=probs[0]+(((long long)curr_mixer[0]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_esign[0][0]=temp;
				temp=probs[1]+(((long long)curr_mixer[1]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_esign[1][0]=temp;
				temp=probs[2]+(((long long)curr_mixer[2]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_esign[2][0]=temp;
				temp=probs[3]+(((long long)curr_mixer[3]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_esign[3][0]=temp;
				temp=curr_mixer[0]+(((long long)probs[0]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[0]=temp;
				temp=curr_mixer[1]+(((long long)probs[1]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[1]=temp;
				temp=curr_mixer[2]+(((long long)probs[2]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[2]=temp;
				temp=curr_mixer[3]+(((long long)probs[3]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[3]=temp;
			//	p1+=((bit<<16)-p1+(1<<7>>1))>>7;
			//	esign1[kc]=p1;
				if(bit)
				{
				//	curr_stats=stats0[kc]+128;
					for(int ks=128;ks<256;++ks)
					{
						probs[0]=curr_stats[0][ks];
						probs[1]=curr_stats[1][ks];
						probs[2]=curr_stats[2][ks];
						probs[3]=curr_stats[3][ks];
						p1=(int)((
							+(long long)curr_mixer[0]*probs[0]
							+(long long)curr_mixer[1]*probs[1]
							+(long long)curr_mixer[2]*probs[2]
							+(long long)curr_mixer[3]*probs[3]
						)>>19);
						p1+=0x8000;
						CLAMP2_32(p1, p1, 1, 0xFFFF);
						//if((size_t)(curr_stats-stats0[kc])>255)
						//	LOG_ERROR("");
					//	p1=*curr_stats;
						if(fwd)
						{
							bit=error==ks;
							ac5_enc_bin(&ec, p1, bit);
#ifdef ESTIMATE_SIZE
							csizes[kc]-=log2((double)(bit?p1:0x10000-p1)/0x10000);
							++nqueries;
#endif
						}
						else
							bit=ac5_dec_bin(&ec, p1);
						proberror=(bit<<16)-p1;
						temp=probs[0]+(((long long)curr_mixer[0]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[0][ks]=temp;
						temp=probs[1]+(((long long)curr_mixer[1]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[1][ks]=temp;
						temp=probs[2]+(((long long)curr_mixer[2]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[2][ks]=temp;
						temp=probs[3]+(((long long)curr_mixer[3]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[3][ks]=temp;
						temp=curr_mixer[0]+(((long long)probs[0]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[0]=temp;
						temp=curr_mixer[1]+(((long long)probs[1]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[1]=temp;
						temp=curr_mixer[2]+(((long long)probs[2]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[2]=temp;
						temp=curr_mixer[3]+(((long long)probs[3]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[3]=temp;
					//	p1+=((bit<<16)-p1+(1<<7>>1))>>7;
					//	*curr_stats=p1;
						if(bit)
						{
							error=ks;
							break;
						}
					//	++curr_stats;
					}
				}
				else
				{
				//	curr_stats=stats0[kc]+127;
					for(int ks=127;ks>=0;--ks)
					{
						probs[0]=curr_stats[0][ks];
						probs[1]=curr_stats[1][ks];
						probs[2]=curr_stats[2][ks];
						probs[3]=curr_stats[3][ks];
						p1=(int)((
							+(long long)curr_mixer[0]*probs[0]
							+(long long)curr_mixer[1]*probs[1]
							+(long long)curr_mixer[2]*probs[2]
							+(long long)curr_mixer[3]*probs[3]
						)>>19);
						p1+=0x8000;
						CLAMP2_32(p1, p1, 1, 0xFFFF);
						//if((size_t)(curr_stats-stats0[kc])>255)
						//	LOG_ERROR("");
					//	p1=*curr_stats;
						if(fwd)
						{
							bit=error==ks;
							ac5_enc_bin(&ec, p1, bit);
#ifdef ESTIMATE_SIZE
							csizes[kc]-=log2((double)(bit?p1:0x10000-p1)/0x10000);
							++nqueries;
#endif
						}
						else
							bit=ac5_dec_bin(&ec, p1);
						proberror=(bit<<16)-p1;
						temp=probs[0]+(((long long)curr_mixer[0]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[0][ks]=temp;
						temp=probs[1]+(((long long)curr_mixer[1]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[1][ks]=temp;
						temp=probs[2]+(((long long)curr_mixer[2]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[2][ks]=temp;
						temp=probs[3]+(((long long)curr_mixer[3]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[3][ks]=temp;
						temp=curr_mixer[0]+(((long long)probs[0]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[0]=temp;
						temp=curr_mixer[1]+(((long long)probs[1]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[1]=temp;
						temp=curr_mixer[2]+(((long long)probs[2]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[2]=temp;
						temp=curr_mixer[3]+(((long long)probs[3]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[3]=temp;
					//	p1+=((bit<<16)-p1+(1<<7>>1))>>7;
					//	*curr_stats=p1;
						if(bit)
						{
							error=ks;
							break;
						}
					//	--curr_stats;
					}
				}
#elif defined EXPERIMENTAL2
				ALIGN(16) int ctx[]=
				{
					N+W-NW,
					N,
					W,
					(N+W)>>1,
				};
				{
					__m128i mc=_mm_load_si128((__m128i*)ctx);
					mc=_mm_min_epi32(mc, _mm_set1_epi32(255));
					mc=_mm_max_epi32(mc, _mm_setzero_si128());
					_mm_store_si128((__m128i*)ctx, mc);
				}
				short *curr_stats[]=
				{
					stats1+(((size_t)kc<<8|ctx[0])<<10|0),
					stats1+(((size_t)kc<<8|ctx[1])<<10|1),
					stats1+(((size_t)kc<<8|ctx[2])<<10|2),
					stats1+(((size_t)kc<<8|ctx[3])<<10|3),
				};
				short *curr_mixer=mixer1[kc][pred];
				short probs[4];
				int temp, p1, proberror;

				unsigned char *curr_splits=splits[kc];
				int low=0, range=256, tidx=1;
			//	short *curr_mixer=mixer1+((size_t)kc<<8|((size_t)N+W)>>1);
				int bit;

				while(range)
				{
					int mid=curr_splits[tidx];
					probs[0]=curr_stats[0][mid];
					probs[1]=curr_stats[1][mid];
					probs[2]=curr_stats[2][mid];
					probs[3]=curr_stats[3][mid];
					p1=(int)((
						+(long long)curr_mixer[0]*probs[0]
						+(long long)curr_mixer[1]*probs[1]
						+(long long)curr_mixer[2]*probs[2]
						+(long long)curr_mixer[3]*probs[3]
					)>>19);
					p1+=0x8000;
					CLAMP2_32(p1, p1, 1, 0xFFFF);
					if(fwd)
					{
						bit=error>mid;
						ac5_enc_bin(&ec, p1, bit);
#ifdef ESTIMATE_SIZE
						csizes[kc]-=log2((double)(bit?p1:0x10000-p1)/0x10000);
						++nqueries;
#endif
					}
					else
						bit=ac5_dec_bin(&ec, p1);
					proberror=(bit<<16)-p1;
					temp=probs[0]+(((long long)curr_mixer[0]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[0][mid]=temp;
					temp=probs[1]+(((long long)curr_mixer[1]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[1][mid]=temp;
					temp=probs[2]+(((long long)curr_mixer[2]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[2][mid]=temp;
					temp=probs[3]+(((long long)curr_mixer[3]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[3][mid]=temp;
					temp=curr_mixer[0]+(((long long)probs[0]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[0]=temp;
					temp=curr_mixer[1]+(((long long)probs[1]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[1]=temp;
					temp=curr_mixer[2]+(((long long)probs[2]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[2]=temp;
					temp=curr_mixer[3]+(((long long)probs[3]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[3]=temp;
				//	curr_stats[0][mid]+=((bit<<16)-probs[0]+(1<<7>>1))>>7;
				//	curr_stats[1][mid]+=((bit<<16)-probs[1]+(1<<7>>1))>>7;
				//	curr_stats[2][mid]+=((bit<<16)-probs[2]+(1<<7>>1))>>7;
				//	curr_stats[3][mid]+=((bit<<16)-probs[3]+(1<<7>>1))>>7;
				//	p1+=((bit<<16)-p1+(1<<7>>1))>>7;
				//	curr_stats[mid]=p1;

					//range=bit?range-(mid-low):mid-low;
					//low=bit?mid:low;

					if(bit)
					{
						range=low+range-(mid+1);
						low=mid+1;
					}
					else
						range=mid+1-low;
					
					//int floorhalf=mid-low;//X
					//if(bit)
					//	low+=range-floorhalf;
					//range=floorhalf;
					if(tidx>=(1LL<<SPLITLEVELS>>1)||!range)
					{
						while(range)
						{
							int floorhalf=range>>1;
							mid=low+floorhalf;
							probs[0]=curr_stats[0][mid];
							probs[1]=curr_stats[1][mid];
							probs[2]=curr_stats[2][mid];
							probs[3]=curr_stats[3][mid];
							p1=(int)((
								+(long long)curr_mixer[0]*probs[0]
								+(long long)curr_mixer[1]*probs[1]
								+(long long)curr_mixer[2]*probs[2]
								+(long long)curr_mixer[3]*probs[3]
							)>>20);
							p1+=0x8000;
							CLAMP2_32(p1, p1, 1, 0xFFFF);
						//	p1=curr_stats[mid];
							if(fwd)
							{
								bit=error>mid;
								ac5_enc_bin(&ec, p1, bit);
#ifdef ESTIMATE_SIZE
								csizes[kc]-=log2((double)(bit?p1:0x10000-p1)/0x10000);
								++nqueries;
#endif
							}
							else
								bit=ac5_dec_bin(&ec, p1);
							proberror=(bit<<16)-p1;
							temp=probs[0]+(((long long)curr_mixer[0]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[0][mid]=temp;
							temp=probs[1]+(((long long)curr_mixer[1]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[1][mid]=temp;
							temp=probs[2]+(((long long)curr_mixer[2]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[2][mid]=temp;
							temp=probs[3]+(((long long)curr_mixer[3]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[3][mid]=temp;
							temp=curr_mixer[0]+(((long long)probs[0]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[0]=temp;
							temp=curr_mixer[1]+(((long long)probs[1]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[1]=temp;
							temp=curr_mixer[2]+(((long long)probs[2]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[2]=temp;
							temp=curr_mixer[3]+(((long long)probs[3]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[3]=temp;
						//	p1+=((bit<<16)-p1+(1<<7>>1))>>7;
						//	curr_stats[mid]=p1;
							if(bit)
								low+=range-floorhalf;
							range=floorhalf;
						}
						error=low;
						break;
					}
					tidx=tidx*2+bit;
				}
#else
#ifdef ENABLE_O1
				ALIGN(16) int ctx[]=
				{
					N+W-NW,
					N,
					W,
					(N+W)>>1,
				};
				{
					__m128i mc=_mm_load_si128((__m128i*)ctx);
					mc=_mm_min_epi32(mc, _mm_set1_epi32(255));
					mc=_mm_max_epi32(mc, _mm_setzero_si128());
					_mm_store_si128((__m128i*)ctx, mc);
				}
				short *curr_stats[]=
				{
					stats1+(((size_t)kc<<8|ctx[0])<<10|0),
					stats1+(((size_t)kc<<8|ctx[1])<<10|1),
					stats1+(((size_t)kc<<8|ctx[2])<<10|2),
					stats1+(((size_t)kc<<8|ctx[3])<<10|3),
				};
				short *curr_mixer=mixer1[kc][pred];
				short probs[4];
				int temp, p1, proberror;
#else
				unsigned short *curr_stats=stats0[kc];
#endif

				int bit=0;
				for(int kb=7, tidx=1;kb>=0;--kb)
				{
#ifdef ENABLE_O1
					probs[0]=curr_stats[0][tidx];
					probs[1]=curr_stats[1][tidx];
					probs[2]=curr_stats[2][tidx];
					probs[3]=curr_stats[3][tidx];
					p1=(int)((
						+(long long)curr_mixer[0]*probs[0]
						+(long long)curr_mixer[1]*probs[1]
						+(long long)curr_mixer[2]*probs[2]
						+(long long)curr_mixer[3]*probs[3]
					)>>16);
					p1+=0x8000;
					CLAMP2_32(p1, p1, 1, 0xFFFF);
#elif defined ENABLE_CTR
					unsigned char *ctr=(unsigned char*)(curr_stats+tidx);
					int p1=((ctr[1]*3+1)<<16)/((ctr[0]+ctr[1])*3+2);
#else
					int p1=curr_stats[tidx];
#endif
					if(fwd)
					{
						bit=error>>kb&1;
						ac5_enc_bin(&ec, p1, bit);
#ifdef ESTIMATE_SIZE
						csizes[kc]-=log2((double)(bit?p1:0x10000-p1)/0x10000);
						++nqueries;
#endif
					}
					else
					{
						bit=ac5_dec_bin(&ec, p1);
						error|=bit<<kb;
					}
#ifdef ENABLE_O1
					proberror=(bit<<16)-p1;
					temp=probs[0]+(((long long)curr_mixer[0]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[0][tidx]=temp;
					temp=probs[1]+(((long long)curr_mixer[1]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[1][tidx]=temp;
					temp=probs[2]+(((long long)curr_mixer[2]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[2][tidx]=temp;
					temp=probs[3]+(((long long)curr_mixer[3]*proberror+(1<<LR_PROB>>1))>>LR_PROB); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_stats[3][tidx]=temp;
					temp=curr_mixer[0]+(((long long)probs[0]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[0]=temp;
					temp=curr_mixer[1]+(((long long)probs[1]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[1]=temp;
					temp=curr_mixer[2]+(((long long)probs[2]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[2]=temp;
					temp=curr_mixer[3]+(((long long)probs[3]*proberror+(1LL<<LR_MIXER>>1))>>LR_MIXER); CLAMP2_32(temp, temp, -0x8000, 0x7FFF); curr_mixer[3]=temp;
#elif defined ENABLE_CTR
					++ctr[bit];
					if(ctr[0]+ctr[1]>=255)
					{
						ctr[0]>>=1;
						ctr[1]>>=1;
					//	ctr[0]=(ctr[0]+1)>>1;
					//	ctr[1]=(ctr[1]+1)>>1;
					}
#else
					p1+=((bit<<16)-p1+(1<<7>>1))>>7;
					curr_stats[tidx]=p1;
#endif
					tidx+=tidx+bit;
				}
#endif
				if(!fwd)
				{
					error+=pred;
					error-=128;
					error&=255;
					*curr=yuv[kc]=error;
				}
				*curr-=offset;
			}
			if(!fwd)
			{
				buffer[idx+permutation[0]]=yuv[0];
				buffer[idx+permutation[1]]=yuv[1];
				buffer[idx+permutation[2]]=yuv[2];
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	if(fwd)
		ac5_enc_flush(&ec);
	else
		fwrite(buffer, 1, usize, fdst);
	fclose(fsrc);
	fclose(fdst);
	free(pixels);
#ifdef ENABLE_O1
	free(stats1);
	free(mixer);
#endif
#ifdef BYPASS_ON_INFLATION
	dstsize=get_filesize(dstfn);
#ifdef ESTIMATE_SIZE
	if(fwd)
	{
		t=time_sec()-t;
		csizes[0]/=8;
		csizes[1]/=8;
		csizes[2]/=8;
		printf("%.2lf  ->  %.2lf + %.2lf + %.2lf = %.2lf  ->  %td\n",
			nqueries/8.,
			csizes[0],
			csizes[1],
			csizes[2],
			csizes[0]+csizes[1]+csizes[2],
			dstsize
		);
		printf("%lf sec  %lf MB/s\n", t, usize/(t*1024*1024));
	}
#endif
	if(fwd&&dstsize>usize)
	{
		printf("Bypass encode\n");
		fdst=fopen(dstfn, "wb");
		fwrite("BP", 1, 2, fdst);
		fwrite(&iw, 1, 4, fdst);
		fwrite(&ih, 1, 4, fdst);
		fwrite(buffer, 1, usize, fdst);
		fclose(fdst);
	}
#endif
	free(buffer);
	(void)nwritten;
	(void)och_names;
	(void)rct_names;
	(void)pred_names;
	return 0;
}