static const char file[]=__FILE__;
#include"ppm.h"
#include"util.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"


//	#define ENABLE_GUIDE
#ifndef DISABLE_MT
//	#define ENABLE_MT
#endif


//	#define ENABLE_WP2
	#define ENABLE_WP
	#define ENABLE_SSE
	#define DISABLE_PREDSEL


#define CODECNAME "C04"
#include"entropy.h"

#define BLOCKSIZE 512
#define MAXPRINTEDBLOCKS 50
#ifdef ENABLE_WP
#define NPREDS 3
#endif

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

#ifdef DISABLE_PREDSEL
#define PREDLIST\
	PRED(CG)
#else
#define PREDLIST\
	PRED(W)\
	PRED(CG)\
	PRED(AV5)\
	PRED(AV9)\
	PRED(AV12)
#endif
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

//typedef enum _NBIndex
//{
//	NB_NNWW,	NB_NNW,		NB_NN,		NB_NNE,		NB_NNEE,
//	NB_NWW,		NB_NW,		NB_N,		NB_NE,		NB_NEE,
//	NB_WW,		NB_W,		NB_curr,
//} NBIndex;
typedef enum _NBIndex
{
	NB_NW,		NB_N,
	NB_W,		NB_curr,
} NBIndex;
#if 0
static const short av12_icoeffs[12]=
{
	 0x04,	 0x03,	-0x1F,	-0x26,	 0x00,
	 0x07,	-0x9E,	 0xDB,	 0x1E,	 0x13,
	-0x2A,	 0xF3,
};
ALIGN(32) static const short g_coeffs[]=
{
	//W
	 0x00,		0x00,		0x00,		0x00,		0x00,
	 0x00,		0x00,		0x00,		0x00,		0x00,
	 0x00,		1<<8,		0, 0, 0, 0,

	//CG
	 0x00,		 0x00,		0x00,		0x00,		0x00,
	 0x00,		-(1<<8),	1<<8,		0x00,		0x00,
	 0x00,		 1<<8,		0, 0, 0, 0,

	//AV5
	 0x00,		 0x00,		0x00,		0x00,		0x00,
	 0x00,		-(5<<5),	5<<5,		1<<5,		0x00,
	 -(1<<5),	 8<<5,		0, 0, 0, 0,

	//AV9
	 0x00,		 1<<4,		-(2<<4),	-(1<<4),	0x00,
	 1<<4,		-(9<<4),	10<<4,		4<<4,		0x00,
	 -(2<<4),	 16<<4,		0, 0, 0, 0,

	//AV12
	 0x04,		 0x03,		-0x1F,		-0x26,		 0x00,
	 0x07,		-0x9E,		 0xDB,		 0x1E,		 0x13,
	-0x2A,		 0xF3,		0, 0, 0, 0,
};
#endif


typedef struct _ThreadArgs
{
	const unsigned char *src;
	unsigned char *dst;
	int iw, ih;

	int fwd, test, loud, x1, x2, y1, y2;
	int bufsize;
	short pixels[(BLOCKSIZE+16)*4*4*2];//4 padded rows * 4 channels max * {pixel, abs(e)}
#ifdef ENABLE_WP
	unsigned short wp_errors[(BLOCKSIZE+16)*NPREDS*3];
#endif

	BList list;
	const unsigned char *decstart, *decend;
	
	int hist[OCH_COUNT*PRED_COUNT<<8];

	//aux
	int blockidx;
	double bestsize;
	int bestrct, predidx[3];
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
		//__m256i av12_mcoeffs[12];
		
		memset(args->hist, 0, sizeof(args->hist));
		//for(int k=0;k<(int)_countof(av12_mcoeffs);++k)
		//	av12_mcoeffs[k]=_mm256_set1_epi16(av12_icoeffs[k]>>1);
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
#ifndef DISABLE_PREDSEL
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
#endif

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
				
#ifndef DISABLE_PREDSEL
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
#endif
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
			printf("Y %5d~%5d  best %12.2lf bytes  %s\n",
				args->y1, args->y2,
				bestsize,
				rct_names[bestrct]
			);
			//printf("Y %5d~%5d  best %12.2lf bytes  %s [YUV: %s %s %s]\n",
			//	args->y1, args->y2,
			//	bestsize,
			//	rct_names[bestrct],
			//	pred_names[predidx[0]],
			//	pred_names[predidx[1]],
			//	pred_names[predidx[2]]
			//);

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
				//printf("%12.2lf %c  %-10s %-10s %-10s %-10s\n",
				//	csize,
				//	kt==bestrct?'*':' ',
				//	rct_names[kt],
				//	pred_names[predsel[group[0]]],
				//	pred_names[predsel[group[1]]],
				//	pred_names[predsel[group[2]]]
				//);
			}
		}
		blist_init(&args->list);
		gr_enc_init(&ec, &args->list);
		gr_enc_NPOT(&ec, bestrct, RCT_COUNT);
#ifndef DISABLE_PREDSEL
		gr_enc_NPOT(&ec, predidx[0], PRED_COUNT);
		gr_enc_NPOT(&ec, predidx[1], PRED_COUNT);
		gr_enc_NPOT(&ec, predidx[2], PRED_COUNT);
#endif
	}
	else
	{
		gr_dec_init(&ec, args->decstart, args->decend);
		bestrct=gr_dec_NPOT(&ec, RCT_COUNT);
#ifndef DISABLE_PREDSEL
		predidx[0]=gr_dec_NPOT(&ec, PRED_COUNT);
		predidx[1]=gr_dec_NPOT(&ec, PRED_COUNT);
		predidx[2]=gr_dec_NPOT(&ec, PRED_COUNT);
#endif
	}
#ifdef ENABLE_SSE
	int sse[3][128][2]={0};
#endif
#ifdef ENABLE_WP
	memset(args->wp_errors, 0, sizeof(args->wp_errors));
#endif
	memset(args->pixels, 0, sizeof(args->pixels));
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
#ifdef ENABLE_WP2
		int wp_errors[6]={0};
#elif defined ENABLE_WP
		unsigned short *wp_errors=args->wp_errors+NPREDS*3;
		//int wp_errors[NPREDS*3]={0};
#endif
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
			//	*NNWW	=rows[2]-2*4*2,
			//	*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
			//	*NNE	=rows[2]+1*4*2,
			//	*NNEE	=rows[2]+2*4*2,
			//	*NNEEE	=rows[2]+3*4*2,
			//	*NWW	=rows[1]-2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
			//	*NEE	=rows[1]+2*4*2,
				*NEEE	=rows[1]+3*4*2,
			//	*WWWW	=rows[0]-4*4*2,
			//	*WWW	=rows[0]-3*4*2,
			//	*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
			if(ky==args->y1)
				NEEE=NE=NW=N=W;
			else if(kx==args->x1)
				NW=W=N;
			else if(kx>args->x2-4)
				NEEE-=(kx-(args->x2-4))*4*2;
			int nbypass[]=
			{
				FLOOR_LOG2(W[4+0]+1),
				FLOOR_LOG2(W[4+1]+1),
				FLOOR_LOG2(W[4+2]+1),
			};
#ifdef ENABLE_WP2
			int weights[]=
			{
				wp_errors[0]+1,
				wp_errors[1]+1,
				wp_errors[2]+1,
				wp_errors[3]+1,
				wp_errors[4]+1,
				wp_errors[5]+1,
			};
			int preds[]=
			{
				(weights[0*2+1]*N[0]+weights[0*2+0]*W[0])/(weights[0*2+0]+weights[0*2+1]),
				(weights[1*2+1]*N[1]+weights[1*2+0]*W[1])/(weights[1*2+0]+weights[1*2+1]),
				(weights[2*2+1]*N[2]+weights[2*2+0]*W[2])/(weights[2*2+0]+weights[2*2+1]),
			};
#elif defined ENABLE_WP
			int preds[3];
			int wp_preds[NPREDS*3]=
			{
				N[0],
				W[0],
				N[0]+W[0]-NW[0],
				//NE[0],
				//N[0]+W[0]+NE[0]-NW[0],
				//(N[0]+W[0])>>1,
				//W[0]+NE[0]-N[0],
				
				N[1],
				W[1],
				N[1]+W[1]-NW[1],
				//NE[1],
				//N[1]+W[1]+NE[1]-NW[1],
				//(N[1]+W[1])>>1,
				//W[1]+NE[1]-N[1],
				
				N[2],
				W[2],
				N[2]+W[2]-NW[2],
				//NE[2],
				//N[2]+W[2]+NE[2]-NW[2],
				//(N[2]+W[2])>>1,
				//W[2]+NE[2]-N[2],
			};
			for(int kc=0;kc<3;++kc)
			{
				//int weights[NPREDS]={0};
				//long long pred=0;
				int pred=0;
				int wsum=0;

				int weights[NPREDS]=
				{
					+wp_errors[kc*NPREDS+0-1*NPREDS*3]//W
					+wp_errors[kc*NPREDS+0+0*NPREDS*3]//N
					+wp_errors[kc*NPREDS+0+1*NPREDS*3]//NE
					+wp_errors[kc*NPREDS+0+2*NPREDS*3]//NEE
					,
					+wp_errors[kc*NPREDS+1-1*NPREDS*3]//W
					+wp_errors[kc*NPREDS+1+0*NPREDS*3]//N
					+wp_errors[kc*NPREDS+1+1*NPREDS*3]//NE
					+wp_errors[kc*NPREDS+1+2*NPREDS*3]//NEE
					,
					+wp_errors[kc*NPREDS+2-1*NPREDS*3]//W
					+wp_errors[kc*NPREDS+2+0*NPREDS*3]//N
					+wp_errors[kc*NPREDS+2+1*NPREDS*3]//NE
					+wp_errors[kc*NPREDS+2+2*NPREDS*3],//NEE
				};
				int w2[]=
				{
					weights[1]*weights[2]+1,
					weights[2]*weights[0]+1,
					weights[0]*weights[1]+1,
				};
				pred=
					+w2[0]*wp_preds[kc*NPREDS+0]
					+w2[1]*wp_preds[kc*NPREDS+1]
					+w2[2]*wp_preds[kc*NPREDS+2];
				wsum=w2[0]+w2[1]+w2[2];

				//for(int kp=0;kp<NPREDS;++kp)
				//{
				//	int weight=0x1000000/(
				//		wp_errors[kc*NPREDS+kp-1*NPREDS*3]+//W
				//		wp_errors[kc*NPREDS+kp+0*NPREDS*3]+//N
				//		wp_errors[kc*NPREDS+kp+1*NPREDS*3]+//NE
				//		wp_errors[kc*NPREDS+kp+2*NPREDS*3]+//NEE
				//	1);
				//	//weights[kp]=weight;
				//	wsum+=weight;
				//	pred+=(long long)weight*wp_preds[kc*NPREDS+kp];
				//}
				pred/=wsum;
				CLAMP3_32(pred, (int)pred, N[kc], W[kc], NE[kc]);
				preds[kc]=(int)pred;
			}
#else
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
			//MEDIAN3_32(preds[0], N[0], W[0], N[0]+W[0]-NW[0]);
			//MEDIAN3_32(preds[1], N[1], W[1], N[1]+W[1]-NW[1]);
			//MEDIAN3_32(preds[2], N[2], W[2], N[2]+W[2]-NW[2]);
#endif
#ifdef ENABLE_SSE
			int *curr_sse[]=
			{
				sse[0][(preds[0]+half*2)>>2],
				sse[1][(preds[1]+half*2)>>2],
				sse[2][(preds[2]+half*2)>>2],
			};
			//if((unsigned)(preds[0]+half*2)>511||(unsigned)(preds[1]+half*2)>511||(unsigned)(preds[2]+half*2)>511)
			//	LOG_ERROR("");
			preds[0]+=curr_sse[0][1]/(curr_sse[0][0]+1);
			preds[1]+=curr_sse[1][1]/(curr_sse[1][0]+1);
			preds[2]+=curr_sse[2][1]/(curr_sse[2][0]+1);
#endif
			int sym=0, offset;
#define UPDATE_FORMULA(IDX) (2*W[IDX]+sym+NEEE[IDX])>>2
//#define UPDATE_FORMULA(IDX) (W[IDX]+sym+NN[IDX]+NEEE[IDX])>>2
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
#ifdef ENABLE_WP2
			wp_errors[0*2+0]=(wp_errors[0*2+0]>>1)+abs(curr[0]-N[0]);
			wp_errors[0*2+1]=(wp_errors[0*2+1]>>1)+abs(curr[0]-W[0]);
			wp_errors[1*2+0]=(wp_errors[1*2+0]>>1)+abs(curr[1]-N[1]);
			wp_errors[1*2+1]=(wp_errors[1*2+1]>>1)+abs(curr[1]-W[1]);
			wp_errors[2*2+0]=(wp_errors[2*2+0]>>1)+abs(curr[2]-N[2]);
			wp_errors[2*2+1]=(wp_errors[2*2+1]>>1)+abs(curr[2]-W[2]);
#elif defined ENABLE_WP
			for(int kc=0;kc<3;++kc)
			{
				int c=curr[kc];
				int wsum=0;
				for(int kp=0;kp<NPREDS;++kp)
					wsum+=wp_errors[kc*NPREDS+kp]=((
						+wp_errors[kc*NPREDS+kp-1*NPREDS*3]//W
						+wp_errors[kc*NPREDS+kp+0*NPREDS*3]//N
					)>>2)+abs(c-wp_preds[kc*NPREDS+kp]);

				//for(int kp=0;kp<NPREDS;++kp)
				//	wsum+=wp_errors[kc*NPREDS+kp]+=
				//		abs(c-wp_preds[kc*NPREDS+kp]);
				if(wsum>64)
				{
					for(int kp=0;kp<NPREDS;++kp)
						wp_errors[kc*NPREDS+kp]>>=1;
				}

				//	wp_errors[kc*NPREDS+kp]=((
				//		wp_errors[kc*NPREDS+kp-1*NPREDS*3]+//W
				//		wp_errors[kc*NPREDS+kp+0*NPREDS*3]+//N
				//	abs(c-wp_preds[kc*NPREDS+kp]))>>2);

					//wp_errors[kc*NPREDS+kp]=(wp_errors[kc*NPREDS+kp]>>1)+abs(c-wp_preds[kc*NPREDS+kp]);//good
					//wp_errors[kc*NPREDS+kp]=abs(c-wp_preds[kc*NPREDS+kp]);
					//args->wp_errors[kc*NPREDS+kp]=(wp_errors[kc*NPREDS+kp]+abs(c-wp_preds[kc*NPREDS+kp]))*493>>9;
			}
			wp_errors+=NPREDS*3;
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
int c04_codec(int argc, char **argv)
{
	if(argc!=2&&argc!=3&&argc!=4)
	{
		printf(
			"Usage: \"%s\"  input  output  [maxthreads]    Encode/decode.\n"
			"       \"%s\"  input                          Test without saving.\n"
			"[maxthreads]:\n"
			"  0: nthreads = number of cores (default)\n"
			"  1: Single thread\n"
			, argv[0]
			, argv[0]
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argc>2?argv[2]:0;
	int maxthreads=argc==4?atoi(argv[3]):0;
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
	ncores=query_cpu_cores();
	usize=iw*ih*3;
	xblocks=(iw+BLOCKSIZE-1)/BLOCKSIZE;
	yblocks=(ih+BLOCKSIZE-1)/BLOCKSIZE;
	nblocks=xblocks*yblocks;
	nthreads=MINVAR(nblocks, ncores);
	if(maxthreads&&nthreads>maxthreads)
		nthreads=maxthreads;
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
#ifdef ENABLE_MT
		arg->loud=0;
#else
		arg->loud=test&&nblocks<MAXPRINTEDBLOCKS;
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
#ifdef ENABLE_MT
			void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
			mt_finish(ctx);
#else
			for(int k=0;k<nthreads2;++k)
				block_thread(args+k);
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