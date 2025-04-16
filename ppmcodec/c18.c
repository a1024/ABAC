#include"codec.h"
#include"util.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
#ifndef DISABLE_MT
	#define ENABLE_MT
#endif

	#define ENABLE_AV2	//good
//	#define ENABLE_SSE2	//good for noisy areas		div-free SSE from NBLIC by WangXuan95		2.613108 -> 2.632816
//	#define ENABLE_SSE	//good with MT			obsolete
//	#define USE_ANGULAR_PRED//bad    15% slower
//	#define USE_GR_ESTIM	//bad    Shannon estimate 0.4% better, but 4.9% slower, when both at stride 4
//	#define ENABLE_SSE_SKEW	//bad
//	#define DISABLE_RCTSEL	//causalRCTSEL > RCT,  disabling RCTSEL breaks SSE
//	#define ENABLE_NPOT	//bad

//	#define EXPORT_GR_PARAM
//	#define DISABLE_AC

#define ANALYSIS_XSTRIDE 2	//4	//1
#define ANALYSIS_YSTRIDE 2	//4	//3


#define CODECNAME "C18"
#include"entropy.h"

#define BLOCKSIZE 224	//192	//448
#define MAXPRINTEDBLOCKS 500
//#define MIXBITS 8
#define AC_THRESHOLD 0.21

#ifndef DISABLE_RCTSEL
#define HALF_Y 128
#define HALF_U 128
#define HALF_V 128
#else
#define HALF_Y 128
#define HALF_U 256
#define HALF_V 256
#endif

#ifndef DISABLE_RCTSEL
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
//static const char *och_names[OCH_COUNT]=
//{
//#define OCH(LABEL) #LABEL,
//	OCHLIST
//#undef  OCH
//};

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
	PRED(N)\
	PRED(W)\
	PRED(AV2)\
	PRED(AV3)\
	PRED(AV4)\
	PRED(AV5)\
	PRED(AV6)\
	PRED(AV9)\
	PRED(AV12)\
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
//typedef enum _NBIndex
//{
//	NB_NW,		NB_N,		NB_NE,
//	NB_W,		NB_curr,
//
//	NB_COUNT,
//} NBIndex;
typedef enum _NBIndex
{
	NB_NNWW,	NB_NNW,		NB_NN,		NB_NNE,		NB_NNEE,
	NB_NWW,		NB_NW,		NB_N,		NB_NE,		NB_NEE,
	NB_WW,		NB_W,		NB_curr,

	NB_COUNT,
} NBIndex;
//static const short av5_icoeffs[12]=//X  immediates are faster
//{
//	 0x00,	 0x000,	 0x00,	 0x00,	 0x00,
//	 0x00,	-0x0A0,	 0xA0,	 0x20,	 0x00,
//	-0x20,	 0x100,
//};
//static const short av9_icoeffs[12]=
//{
//	 0x00,	 0x010,	-0x20,	-0x10,	 0x00,
//	-0x10,	-0x090,	 0xA0,	 0x40,	 0x00,
//	-0x20,	 0x100,
//};
static const short av12_icoeffs[12]=
{
	 0x04,	 0x03,	-0x1F,	-0x26,	 0x00,
	 0x07,	-0x9E,	 0xDB,	 0x1E,	 0x13,
	-0x2A,	 0xF3,
};
#endif


typedef struct _ThreadArgs
{
	const unsigned char *src;
	unsigned char *dst;
	int iw, ih;

	int fwd, test, loud, x1, x2, y1, y2;
	int bufsize;
	short pixels[(BLOCKSIZE+16)*4*3*4];//4 padded rows * 3 channels max * {pixel, abs(e1), abs(e2), abs(e3)}	43.5 KB
	
	int bestrct;
	int use_AC[3];
	BList list[3];
	const unsigned char *decstart[3], *decend[3];
	
#if !defined DISABLE_RCTSEL && !defined USE_GR_ESTIM
	int hist[OCH_COUNT*PRED_COUNT<<8];//120 KB
#endif
#ifdef ENABLE_SSE
	int sse1[3][64][64][2];
	int sse2[3][64][64][2];
	int sse3[3][64][64][2];
	int sse4[3][256][2];
#endif
#ifdef ENABLE_SSE2
	int sse1[3][64][64];//147 KB
	int sse2[3][256][16];
	int sse3[3][256][16];
	int sse4[3][256][16];
	int sse5[3][64][64];
	int sse6[3][512];
#endif

	//aux
	int blockidx;
	double bestsize;
	int predidx[3];
#ifndef DISABLE_RCTSEL
	double t_analysis;
#endif
#ifdef EXPORT_GR_PARAM
	unsigned char *grimage;
#endif
} ThreadArgs;
#if 0
AWM_INLINE int max3(int a, int b, int c)
{
	if(a<b)
		a=b;
	if(a<c)
		a=c;
	return a;
}
AWM_INLINE int max4(int v0, int v1, int v2, int v3)
{
	if(v0<v1)
		v0=v1;
	if(v0<v2)
		v0=v2;
	if(v0<v3)
		v0=v3;
	return v0;
}
AWM_INLINE int median3(int a, int b, int c)
{
	MEDIAN3_32(a, a, b, c);
	return a;
}
AWM_INLINE int sum_largest3of7(int v0, int v1, int v2, int v3, int v4, int v5, int v6)
{
	int m0=v0, m1=v0, m2=v0;
	if(m0<v1)m2=m1, m1=m0, m0=v1;
	if(m0<v2)m2=m1, m1=m0, m0=v2;
	if(m0<v3)m2=m1, m1=m0, m0=v3;
	if(m0<v4)m2=m1, m1=m0, m0=v4;
	if(m0<v5)m2=m1, m1=m0, m0=v5;
	if(m0<v6)m2=m1, m1=m0, m0=v6;
	return m0+m1+m2;
}
#endif
static void block_thread(void *param)
{
	const int nch=3;
	ThreadArgs *args=(ThreadArgs*)param;
	GolombRiceCoder ec[3]={0};
	AC3 ec2[3]={0};
#ifndef DISABLE_RCTSEL
//	const unsigned char *image=args->fwd?args->src:args->dst;
	unsigned char bestrct=0, predidx[4]={0};
#endif
//	short mixcoeff[4]={0};
	
	if(args->fwd)
	{
#ifndef DISABLE_RCTSEL
		int ystride=args->iw*3;
#ifdef USE_GR_ESTIM
		__m256i table_floorlog2p1_lo=_mm256_set_epi8(
		//	 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0,
			  4,  4,  4,  4,  4,  4,  4,  4,  3,  3,  3,  3,  2,  2,  1,  0,
			  4,  4,  4,  4,  4,  4,  4,  4,  3,  3,  3,  3,  2,  2,  1,  0
		);
		__m256i table_floorlog2p1_hi=_mm256_set_epi8(
		//	240,224,208,192,176,160,144,128,112, 96, 80, 64, 48, 32, 16,  0,
			  8,  8,  8,  8,  8,  8,  8,  8,  7,  7,  7,  7,  6,  6,  5,  0,
			  8,  8,  8,  8,  8,  8,  8,  8,  7,  7,  7,  7,  6,  6,  5,  0
		);
		int csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
#else
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
#endif
		unsigned char predsel[OCH_COUNT]={0};
#ifndef USE_GR_ESTIM
		int res;
#endif
		__m256i av12_mcoeffs[12];
		
		args->t_analysis=time_sec();
		for(int k=0;k<(int)_countof(av12_mcoeffs);++k)
			av12_mcoeffs[k]=_mm256_set1_epi16(av12_icoeffs[k]>>1);
#if 0
		for(int ky=args->y1+1;ky<args->y2;++ky)
		{
			unsigned char *Nptr=args->src+3*(args->iw*(ky-1)+1);
			unsigned char *currptr=args->src+3*(args->iw*ky+1);
			for(int kx=args->x1+1;kx<args->x2;++kx)
			{
				for(int kc=0;kc<3;++kc)//X  need 12 output channels
				{
					int
						N	=Nptr[kc],
						W	=currptr[kc-3],
						curr	=currptr[kc];
					if(N!=W)
						mixcoeff[kc]+=((curr-W)<<8|1<<8>>1)/(N-W);
				}
			}
		}
		res=(args->x2-args->x1-1)*(args->y2-args->y1-1);
		mixcoeff[0]/=res;
		mixcoeff[1]/=res;
		mixcoeff[2]/=res;
#endif
#ifndef USE_GR_ESTIM
		memset(args->hist, 0, sizeof(args->hist));
		res=(args->x2-args->x1-4)/ANALYSIS_XSTRIDE/5*5*((args->y2-args->y1-2)/ANALYSIS_YSTRIDE);
		if(!res)
			res=(args->x2-args->x1)*(args->y2-args->y1-2);
#endif
		for(int ky=args->y1+2;ky<args->y2-(ANALYSIS_YSTRIDE-1);ky+=ANALYSIS_YSTRIDE)//analysis loop
		{
			int kx=args->x1+2;
			const unsigned char *ptr=args->src+3*(args->iw*ky+kx);

			__m256i amin=_mm256_set1_epi16(-128);
			__m256i amax=_mm256_set1_epi16(127);
#ifndef USE_GR_ESTIM
			__m256i amag=_mm256_set1_epi16(255);
#endif
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
			for(;kx<args->x2-(2+(5*ANALYSIS_XSTRIDE-1));kx+=5*ANALYSIS_XSTRIDE, ptr+=15*ANALYSIS_XSTRIDE)
			{
				__m256i
					nb0[NB_COUNT],//rgb
					nb1[NB_COUNT],//gbr
					nb2[NB_COUNT],//rgb - gbr
					nb3[NB_COUNT],//gbr - rgb
					nb4[NB_COUNT],//(gbr+brg)/2
					nb5[NB_COUNT];//rgb - (gbr+brg)/2
				__m256i vmin[4], vmax[4], pred;
				{
					__m128i nb8[NB_COUNT]=//8-bit
					{
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-2*ystride-2*3+0)), half8),//NNWW
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-2*ystride-1*3+0)), half8),//NNW
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-2*ystride+0*3+0)), half8),//NN
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-2*ystride+1*3+0)), half8),//NNE
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-2*ystride+2*3+0)), half8),//NNEE
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride-2*3+0)), half8),//NWW
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride-1*3+0)), half8),//NW
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride+0*3+0)), half8),//N
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride+1*3+0)), half8),//NE
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride+2*3+0)), half8),//NEE
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr+0*ystride-2*3+0)), half8),//WW
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr+0*ystride-1*3+0)), half8),//W
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr+0*ystride+0*3+0)), half8),//curr
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
#ifdef USE_GR_ESTIM
#define UPDATE(PREDIDX, IDX0, IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7, IDX8, IDX9, IDXA, IDXB, IDXC, IDXD, IDXE)\
	do\
	{\
		pred=_mm256_xor_si256(_mm256_slli_epi16(pred, 1), _mm256_srai_epi16(pred, 15));\
		pred=_mm256_max_epi8(_mm256_shuffle_epi8(table_floorlog2p1_hi, _mm256_srai_epi16(pred, 4)), _mm256_shuffle_epi8(table_floorlog2p1_lo, _mm256_and_si256(pred, _mm256_set1_epi16(15))));\
		_mm256_store_si256((__m256i*)result, pred);\
		csizes[IDX0*PRED_COUNT+PREDIDX]+=result[0x0];\
		csizes[IDX1*PRED_COUNT+PREDIDX]+=result[0x1];\
		csizes[IDX2*PRED_COUNT+PREDIDX]+=result[0x2];\
		csizes[IDX3*PRED_COUNT+PREDIDX]+=result[0x3];\
		csizes[IDX4*PRED_COUNT+PREDIDX]+=result[0x4];\
		csizes[IDX5*PRED_COUNT+PREDIDX]+=result[0x5];\
		csizes[IDX6*PRED_COUNT+PREDIDX]+=result[0x6];\
		csizes[IDX7*PRED_COUNT+PREDIDX]+=result[0x7];\
		csizes[IDX8*PRED_COUNT+PREDIDX]+=result[0x8];\
		csizes[IDX9*PRED_COUNT+PREDIDX]+=result[0x9];\
		csizes[IDXA*PRED_COUNT+PREDIDX]+=result[0xA];\
		csizes[IDXB*PRED_COUNT+PREDIDX]+=result[0xB];\
		csizes[IDXC*PRED_COUNT+PREDIDX]+=result[0xC];\
		csizes[IDXD*PRED_COUNT+PREDIDX]+=result[0xD];\
		csizes[IDXE*PRED_COUNT+PREDIDX]+=result[0xE];\
	}while(0)
#else
#define UPDATE(PREDIDX, IDX0, IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7, IDX8, IDX9, IDXA, IDXB, IDXC, IDXD, IDXE)\
	do\
	{\
		pred=_mm256_sub_epi16(pred, amin);\
		pred=_mm256_and_si256(pred, amag);\
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
#ifdef ENABLE_AV2
				//N
				pred=nb0[NB_N];

				//(N+W+NW+NE)/4
			//	pred=_mm256_add_epi16(nb0[NB_N], nb0[NB_W]);
			//	pred=_mm256_add_epi16(pred, nb0[NB_NW]);
			//	pred=_mm256_add_epi16(pred, nb0[NB_NE]);
			//	pred=_mm256_srai_epi16(pred, 2);

				//(N+W)/2
			//	pred=_mm256_srai_epi16(_mm256_add_epi16(nb0[NB_N], nb0[NB_W]), 1);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_N,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=nb2[NB_N];

			//	pred=_mm256_add_epi16(nb2[NB_N], nb2[NB_W]);
			//	pred=_mm256_add_epi16(pred, nb2[NB_NW]);
			//	pred=_mm256_add_epi16(pred, nb2[NB_NE]);
			//	pred=_mm256_srai_epi16(pred, 2);

			//	pred=_mm256_srai_epi16(_mm256_add_epi16(nb2[NB_N], nb2[NB_W]), 1);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_N,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=nb3[NB_N];

			//	pred=_mm256_add_epi16(nb3[NB_N], nb3[NB_W]);
			//	pred=_mm256_add_epi16(pred, nb3[NB_NW]);
			//	pred=_mm256_add_epi16(pred, nb3[NB_NE]);
			//	pred=_mm256_srai_epi16(pred, 2);

			//	pred=_mm256_srai_epi16(_mm256_add_epi16(nb3[NB_N], nb3[NB_W]), 1);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_N,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=nb5[NB_N];

			//	pred=_mm256_add_epi16(nb5[NB_N], nb5[NB_W]);
			//	pred=_mm256_add_epi16(pred, nb5[NB_NW]);
			//	pred=_mm256_add_epi16(pred, nb5[NB_NE]);
			//	pred=_mm256_srai_epi16(pred, 2);

			//	pred=_mm256_srai_epi16(_mm256_add_epi16(nb5[NB_N], nb5[NB_W]), 1);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_N,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//W
				pred=nb0[NB_W];

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_W,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=nb2[NB_W];

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
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
				pred=nb3[NB_W];

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
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
				pred=nb5[NB_W];

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
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
				
				//N+W-NW
				//pred=_mm256_sub_epi16(_mm256_add_epi16(nb0[NB_N], nb0[NB_W]), nb0[NB_NW]);

				//(W+NE)/2
				//pred=_mm256_srai_epi16(_mm256_add_epi16(nb0[NB_W], nb0[NB_NE]), 1);

				//(N+NE)/2
				//pred=_mm256_srai_epi16(_mm256_add_epi16(nb0[NB_N], nb0[NB_NE]), 1);

				//(4*(N+W)+NE-NW)/8
				//pred=_mm256_add_epi16(nb0[NB_N], nb0[NB_W]);
				//pred=_mm256_slli_epi16(pred, 2);
				//pred=_mm256_add_epi16(pred, nb0[NB_NE]);
				//pred=_mm256_sub_epi16(pred, nb0[NB_NW]);
				//pred=_mm256_srai_epi16(pred, 3);

				//(N+W)/2
				pred=_mm256_srai_epi16(_mm256_add_epi16(nb0[NB_N], nb0[NB_W]), 1);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
			//	pred=_mm256_sub_epi16(_mm256_add_epi16(nb2[NB_N], nb2[NB_W]), nb2[NB_NW]);

			//	pred=_mm256_srai_epi16(_mm256_add_epi16(nb2[NB_W], nb2[NB_NE]), 1);

			//	pred=_mm256_add_epi16(nb2[NB_N], nb2[NB_W]);
			//	pred=_mm256_slli_epi16(pred, 2);
			//	pred=_mm256_add_epi16(pred, nb2[NB_NE]);
			//	pred=_mm256_sub_epi16(pred, nb2[NB_NW]);
			//	pred=_mm256_srai_epi16(pred, 3);

				pred=_mm256_srai_epi16(_mm256_add_epi16(nb2[NB_N], nb2[NB_W]), 1);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
			//	pred=_mm256_sub_epi16(_mm256_add_epi16(nb3[NB_N], nb3[NB_W]), nb3[NB_NW]);

			//	pred=_mm256_srai_epi16(_mm256_add_epi16(nb3[NB_W], nb3[NB_NE]), 1);

			//	pred=_mm256_add_epi16(nb3[NB_N], nb3[NB_W]);
			//	pred=_mm256_slli_epi16(pred, 2);
			//	pred=_mm256_add_epi16(pred, nb3[NB_NE]);
			//	pred=_mm256_sub_epi16(pred, nb3[NB_NW]);
			//	pred=_mm256_srai_epi16(pred, 3);

				pred=_mm256_srai_epi16(_mm256_add_epi16(nb3[NB_N], nb3[NB_W]), 1);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
			//	pred=_mm256_sub_epi16(_mm256_add_epi16(nb5[NB_N], nb5[NB_W]), nb5[NB_NW]);

			//	pred=_mm256_srai_epi16(_mm256_add_epi16(nb5[NB_W], nb5[NB_NE]), 1);

			//	pred=_mm256_add_epi16(nb5[NB_N], nb5[NB_W]);
			//	pred=_mm256_slli_epi16(pred, 2);
			//	pred=_mm256_add_epi16(pred, nb5[NB_NE]);
			//	pred=_mm256_sub_epi16(pred, nb5[NB_NW]);
			//	pred=_mm256_srai_epi16(pred, 3);

				pred=_mm256_srai_epi16(_mm256_add_epi16(nb5[NB_N], nb5[NB_W]), 1);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//AV5
				//		-5	5	1
				//	-1	8	[?]>>3		clamp(N, W, NE)
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

				//AV3
				//	-2	3
				//	3	[?]>>2		clamp(N, W, NE)
				pred=_mm256_add_epi16(nb0[NB_N], nb0[NB_W]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 1));//3*(N+W)
				pred=_mm256_sub_epi16(pred, _mm256_slli_epi16(nb0[NB_NW], 1));
				pred=_mm256_srai_epi16(pred, 2);
				//vmin=_mm256_min_epi16(N, W);
				//vmax=_mm256_max_epi16(N, W);
				//vmin[0]=_mm256_min_epi16(vmin[0], nb0[NB_NE]);
				//vmax[0]=_mm256_max_epi16(vmax[0], nb0[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV3,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(nb2[NB_N], nb2[NB_W]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 1));//3*(N+W)
				pred=_mm256_sub_epi16(pred, _mm256_slli_epi16(nb2[NB_NW], 1));
				pred=_mm256_srai_epi16(pred, 2);
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin[1]=_mm256_min_epi16(vmin[1], nb2[NB_NE]);
				//vmax[1]=_mm256_max_epi16(vmax[1], nb2[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV3,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(nb3[NB_N], nb3[NB_W]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 1));//3*(N+W)
				pred=_mm256_sub_epi16(pred, _mm256_slli_epi16(nb3[NB_NW], 1));
				pred=_mm256_srai_epi16(pred, 2);
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin[2]=_mm256_min_epi16(vmin[2], nb3[NB_NE]);
				//vmax[2]=_mm256_max_epi16(vmax[2], nb3[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV3,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_add_epi16(nb5[NB_N], nb5[NB_W]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 1));//3*(N+W)
				pred=_mm256_sub_epi16(pred, _mm256_slli_epi16(nb5[NB_NW], 1));
				pred=_mm256_srai_epi16(pred, 2);
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin[3]=_mm256_min_epi16(vmin[1], nb5[NB_NE]);
				//vmax[3]=_mm256_max_epi16(vmax[1], nb5[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV3,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//AV4	(4*(N+W)+NE-NW)>>3		clamp(N, W, NE)
				pred=_mm256_slli_epi16(_mm256_add_epi16(nb0[NB_N], nb0[NB_W]), 2);
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb0[NB_NE], nb0[NB_NW]));
				pred=_mm256_srai_epi16(pred, 3);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV4,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_slli_epi16(_mm256_add_epi16(nb2[NB_N], nb2[NB_W]), 2);
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb2[NB_NE], nb2[NB_NW]));
				pred=_mm256_srai_epi16(pred, 3);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV4,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_slli_epi16(_mm256_add_epi16(nb3[NB_N], nb3[NB_W]), 2);
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb3[NB_NE], nb3[NB_NW]));
				pred=_mm256_srai_epi16(pred, 3);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV4,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_slli_epi16(_mm256_add_epi16(nb5[NB_N], nb5[NB_W]), 2);
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb5[NB_NE], nb5[NB_NW]));
				pred=_mm256_srai_epi16(pred, 3);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV4,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//AV6
				//			-1
				//		-5	6	1
				//	-1	8	[?]>>3		clamp(N, W, NE)
				pred=_mm256_sub_epi16(nb0[NB_N], nb0[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));//5*(N-NW)
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb0[NB_NE], nb0[NB_WW]));//5*(N-NW)+NE-WW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb0[NB_N], nb0[NB_NN]));//6*N-5*NW-NN-WW+NE
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb0[NB_W]);//W+(6*N-5*NW-NN-WW+NE)/8
				//vmin=_mm256_min_epi16(N, W);
				//vmax=_mm256_max_epi16(N, W);
				//vmin[0]=_mm256_min_epi16(vmin[0], nb0[NB_NE]);
				//vmax[0]=_mm256_max_epi16(vmax[0], nb0[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV6,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_sub_epi16(nb2[NB_N], nb2[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));//5*(N-NW)
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb2[NB_NE], nb2[NB_WW]));//5*(N-NW)+NE-WW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb2[NB_N], nb2[NB_NN]));//6*N-5*NW-NN-WW+NE
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb2[NB_W]);//W+(6*N-5*NW-NN-WW+NE)/8
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin[1]=_mm256_min_epi16(vmin[1], nb2[NB_NE]);
				//vmax[1]=_mm256_max_epi16(vmax[1], nb2[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV6,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_sub_epi16(nb3[NB_N], nb3[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));//5*(N-NW)
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb3[NB_NE], nb3[NB_WW]));//5*(N-NW)+NE-WW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb3[NB_N], nb3[NB_NN]));//6*N-5*NW-NN-WW+NE
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb3[NB_W]);//W+(6*N-5*NW-NN-WW+NE)/8
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin[2]=_mm256_min_epi16(vmin[2], nb3[NB_NE]);
				//vmax[2]=_mm256_max_epi16(vmax[2], nb3[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV6,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_sub_epi16(nb5[NB_N], nb3[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));//5*(N-NW)
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb5[NB_NE], nb5[NB_WW]));//5*(N-NW)+NE-WW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb5[NB_N], nb5[NB_NN]));//6*N-5*NW-NN-WW+NE
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb5[NB_W]);//W+(6*N-5*NW-NN-WW+NE)/8
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin[3]=_mm256_min_epi16(vmin[1], nb5[NB_NE]);
				//vmax[3]=_mm256_max_epi16(vmax[1], nb5[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV6,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//AV9
				//		1	-2	-1
				//	-1	-9	10	4
				//	-2	16	[?]>>4		clamp(N, W, NE)
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
				//	+4	+3	-31	-38	+0
				//	+7	-158	+219	+30	+19
				//	-42	+243	[?]>>8		clamp(N, W, NE)
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
#ifndef USE_GR_ESTIM
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
#endif
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
		double compsizes[3]={0};
		for(int kt=0;kt<RCT_COUNT;++kt)//select best RCT
		{
			const unsigned char *group=rct_combinations[kt];
			double tcompsizes[]=
			{
				csizes[group[0]*PRED_COUNT+predsel[group[0]]],
				csizes[group[1]*PRED_COUNT+predsel[group[1]]],
				csizes[group[2]*PRED_COUNT+predsel[group[2]]],
			};
#ifdef USE_GR_ESTIM
			int csize=
#else
			double csize=
#endif
				tcompsizes[0]+
				tcompsizes[1]+
				tcompsizes[2];
			if(!kt||bestsize>csize)
			{
				bestsize=csize, bestrct=kt;
				compsizes[0]=tcompsizes[0];
				compsizes[1]=tcompsizes[1];
				compsizes[2]=tcompsizes[2];
			}
		}
		const double threshold=(BLOCKSIZE/ANALYSIS_XSTRIDE)*(BLOCKSIZE/ANALYSIS_YSTRIDE)*AC_THRESHOLD;
		args->use_AC[0]=compsizes[0]<threshold;
		args->use_AC[1]=compsizes[1]<threshold;
		args->use_AC[2]=compsizes[2]<threshold;
#ifdef DISABLE_AC
		args->use_AC[0]=0;
		args->use_AC[1]=0;
		args->use_AC[2]=0;
#endif
		predidx[0]=predsel[rct_combinations[bestrct][0]];
		predidx[1]=predsel[rct_combinations[bestrct][1]];
		predidx[2]=predsel[rct_combinations[bestrct][2]];
		args->bestsize=bestsize;
		args->bestrct=bestrct;
		args->predidx[0]=predidx[0];
		args->predidx[1]=predidx[1];
		args->predidx[2]=predidx[2];
		args->t_analysis=time_sec()-args->t_analysis;
#if 0
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
#endif
#endif
		blist_init(args->list+0);
		blist_init(args->list+1);
		blist_init(args->list+2);
		if(args->use_AC[0])
		{
			ac3_enc_init(ec2+0, args->list+0);
			ac3_enc_bypass_NPOT(ec2+0, predidx[0], PRED_COUNT);
		}
		else
		{
			gr_enc_init(ec+0, args->list+0);
			gr_enc_NPOT(ec+0, predidx[0], PRED_COUNT);
		}
		if(args->use_AC[1])
		{
			ac3_enc_init(ec2+1, args->list+1);
			ac3_enc_bypass_NPOT(ec2+1, predidx[1], PRED_COUNT);
		}
		else
		{
			gr_enc_init(ec+1, args->list+1);
			gr_enc_NPOT(ec+1, predidx[1], PRED_COUNT);
		}
		if(args->use_AC[2])
		{
			ac3_enc_init(ec2+2, args->list+2);
			ac3_enc_bypass_NPOT(ec2+2, predidx[2], PRED_COUNT);
		}
		else
		{
			gr_enc_init(ec+2, args->list+2);
			gr_enc_NPOT(ec+2, predidx[2], PRED_COUNT);
		}
//		gr_enc_init(&ec, &args->list);
//#ifndef DISABLE_RCTSEL
//		gr_enc_NPOT(&ec, bestrct, RCT_COUNT);
//		gr_enc_NPOT(&ec, predidx[0], PRED_COUNT);
//		gr_enc_NPOT(&ec, predidx[1], PRED_COUNT);
//		gr_enc_NPOT(&ec, predidx[2], PRED_COUNT);
//#ifdef ENABLE_SSE2
//		gr_enc(&ec, enable_SSE[0], 1);
//		gr_enc(&ec, enable_SSE[1], 1);
//		gr_enc(&ec, enable_SSE[2], 1);
//#endif
//#endif
	}
	else
	{
		if(args->use_AC[0])
		{
			ac3_dec_init(ec2+0, args->decstart[0], args->decend[0]);
			predidx[0]=ac3_dec_bypass_NPOT(ec2+0, PRED_COUNT);
		}
		else
		{
			gr_dec_init(ec+0, args->decstart[0], args->decend[0]);
			predidx[0]=gr_dec_NPOT(ec+0, PRED_COUNT);
		}
		if(args->use_AC[1])
		{
			ac3_dec_init(ec2+1, args->decstart[1], args->decend[1]);
			predidx[1]=ac3_dec_bypass_NPOT(ec2+1, PRED_COUNT);
		}
		else
		{
			gr_dec_init(ec+1, args->decstart[1], args->decend[1]);
			predidx[1]=gr_dec_NPOT(ec+1, PRED_COUNT);
		}
		if(args->use_AC[2])
		{
			ac3_dec_init(ec2+2, args->decstart[2], args->decend[2]);
			predidx[2]=ac3_dec_bypass_NPOT(ec2+2, PRED_COUNT);
		}
		else
		{
			gr_dec_init(ec+2, args->decstart[2], args->decend[2]);
			predidx[2]=gr_dec_NPOT(ec+2, PRED_COUNT);
		}
//		gr_dec_init(&ec, args->decstart, args->decend);
//#ifndef DISABLE_RCTSEL
//		bestrct=gr_dec_NPOT(&ec, RCT_COUNT);
//		predidx[0]=gr_dec_NPOT(&ec, PRED_COUNT);
//		predidx[1]=gr_dec_NPOT(&ec, PRED_COUNT);
//		predidx[2]=gr_dec_NPOT(&ec, PRED_COUNT);
//#ifdef ENABLE_SSE2
//		enable_SSE[0]=gr_dec(&ec, 1);
//		enable_SSE[1]=gr_dec(&ec, 1);
//		enable_SSE[2]=gr_dec(&ec, 1);
//#endif
//#endif
	}
#if 0
	//__m128i predmask1=_mm_setzero_si128();
	__m128i predmask1=_mm_set_epi8(
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		-(predidx[2]==PRED_W), -(predidx[2]==PRED_W),
		-(predidx[1]==PRED_W), -(predidx[1]==PRED_W),
		-(predidx[0]==PRED_W), -(predidx[0]==PRED_W)
	);
	//__m128i predmask2=_mm_setzero_si128();
	__m128i predmask2=_mm_set_epi8(
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		-(predidx[2]==PRED_AV2), -(predidx[2]==PRED_AV2),
		-(predidx[1]==PRED_AV2), -(predidx[1]==PRED_AV2),
		-(predidx[0]==PRED_AV2), -(predidx[0]==PRED_AV2)
	);
	//__m128i predmask3=_mm_set1_epi32(-1);
	__m128i predmask3=_mm_set_epi8(
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		-(predidx[2]==PRED_CG), -(predidx[2]==PRED_CG),
		-(predidx[1]==PRED_CG), -(predidx[1]==PRED_CG),
		-(predidx[0]==PRED_CG), -(predidx[0]==PRED_CG)
	);
#endif
	//int mixer[3]=
	//{
	//	1<<MIXBITS>>1,
	//	1<<MIXBITS>>1,
	//	1<<MIXBITS>>1,
	//};
#ifdef ENABLE_SSE
	memset(args->sse1, 0, sizeof(args->sse1));
	memset(args->sse2, 0, sizeof(args->sse2));
	memset(args->sse3, 0, sizeof(args->sse3));
	memset(args->sse4, 0, sizeof(args->sse4));
#endif
#ifdef ENABLE_SSE2
	memset(args->sse1, 0, sizeof(args->sse1));
	memset(args->sse2, 0, sizeof(args->sse2));
	memset(args->sse3, 0, sizeof(args->sse3));
	memset(args->sse4, 0, sizeof(args->sse4));
	memset(args->sse5, 0, sizeof(args->sse5));
	memset(args->sse6, 0, sizeof(args->sse6));
#endif
	int *stats0[]=
	{
		args->hist+258*0,//{257 levels (because sym=256 is valid), total count}
		args->hist+258*1,
		args->hist+258*2,
	};
	//const short *coeffs[3]={0};
	//if(predidx[0]==PRED_AV5)
	//	coeffs[0]=av5_icoeffs;
	//else if(predidx[0]==PRED_AV9)
	//	coeffs[0]=av9_icoeffs;
	//else if(predidx[0]==PRED_AV12)
	//	coeffs[0]=av12_icoeffs;
	//if(predidx[1]==PRED_AV5)
	//	coeffs[1]=av5_icoeffs;
	//else if(predidx[1]==PRED_AV9)
	//	coeffs[1]=av9_icoeffs;
	//else if(predidx[1]==PRED_AV12)
	//	coeffs[1]=av12_icoeffs;
	//if(predidx[2]==PRED_AV5)
	//	coeffs[2]=av5_icoeffs;
	//else if(predidx[2]==PRED_AV9)
	//	coeffs[2]=av9_icoeffs;
	//else if(predidx[2]==PRED_AV12)
	//	coeffs[2]=av12_icoeffs;
	memset(args->hist, 0, sizeof(args->hist));
	memset(args->pixels, 0, sizeof(args->pixels));
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(32) short *rows[]=
		{
			args->pixels+((BLOCKSIZE+16LL)*((ky-0LL)&3)+8LL)*3*4,
			args->pixels+((BLOCKSIZE+16LL)*((ky-1LL)&3)+8LL)*3*4,
			args->pixels+((BLOCKSIZE+16LL)*((ky-2LL)&3)+8LL)*3*4,
			args->pixels+((BLOCKSIZE+16LL)*((ky-3LL)&3)+8LL)*3*4,
		};
		int yuv[4]={0};
		int pred=0, error=0;
#ifndef DISABLE_RCTSEL
		const unsigned char *combination=rct_combinations[args->bestrct];
#endif
#ifdef USE_ANGULAR_PRED
	//	int linear[4]={0};
		int angular[4]={0};
#endif
		ALIGN(16) short preds[8]={0};
#ifdef ENABLE_SSE2
		int *curr_sse1[3]={0};
		int *curr_sse2[3]={0};
		int *curr_sse3[3]={0};
		int *curr_sse4[3]={0};
		int *curr_sse5[3]={0};
		int *curr_sse6[3]={0};
#endif
		int nbypass=0;
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int idx=nch*(args->iw*ky+kx);
			short
				*NNN	=rows[3]+0*3*4,
				*NNNE	=rows[3]+1*3*4,
				*NNWW	=rows[2]-2*3*4,
				*NNW	=rows[2]-1*3*4,
				*NN	=rows[2]+0*3*4,
				*NNE	=rows[2]+1*3*4,
				*NNEE	=rows[2]+2*3*4,
				*NWW	=rows[1]-2*3*4,
				*NW	=rows[1]-1*3*4,
				*N	=rows[1]+0*3*4,
				*NE	=rows[1]+1*3*4,
				*NEE	=rows[1]+2*3*4,
				*NEEE	=rows[1]+3*3*4,
				*NEEEE	=rows[1]+4*3*4,
				*WWW	=rows[0]-3*3*4,
				*WW	=rows[0]-2*3*4,
				*W	=rows[0]-1*3*4,
				*curr	=rows[0]+0*3*4;
			(void)NNNE;
			(void)NNN;
			(void)NN;
			(void)NNE;
			(void)NNEE;
			(void)NEEEE;
			if(ky==args->y1)
				NEEEE=NEEE=NEE=NE=NW=N=W;
			else if(kx==args->x1)
				NW=WWW=WW=W=N;
			else if(kx>args->x2-4)
				NEEE-=(kx-(args->x2-4))*3*4;
#if 0
			int nbypass[3];
#ifdef __GNUC__
#pragma GCC unroll 3
#elif defined __clang__
#pragma clang loop unroll_count(3)
#endif
			for(int kc=4;kc<4+3;++kc)
			{
				//if(2*NW[kc]<N[kc]+W[kc]&&2*NNWW[kc]<NNW[kc]+NWW[kc])
				//	nbypass[kc-4]=(NW[kc]+NNWW[kc])>>1;
				//else if(NE[kc]<N[kc]&&NNEE[kc+4]<NNE[kc])
				//	nbypass[kc-4]=(NE[kc]+NNEE[kc])>>1;
				//else if(W[kc]<NW[kc]&&WW[kc]<NWW[kc])
				//	nbypass[kc-4]=(W[kc]+WW[kc])>>1;
				//else if(2*N[kc]<NW[kc]+NE[kc]&&2*NN[kc]<NNW[kc]+NNE[kc])
				//	nbypass[kc-4]=(N[kc]+NN[kc])>>1;
				//else
					nbypass[kc-4]=W[kc];
				//	nbypass[kc-4]=(4*(N[kc]+W[kc])+2*(NW[kc]+NE[kc]+NN[kc]+WW[kc])+WWW[kc]+NWW[kc]+NNW[kc]+NNN[kc]+NNE[kc]+NEE[kc]+NEEE[kc]+NNEE[kc])/24;
				//	nbypass[kc-4]=(2*(N[kc]+W[kc])+NW[kc]+NE[kc]+NN[kc]+WW[kc])>>3;
			}
#if 0
			{
				int xgrads[]=
				{
					abs(W[4+0]-WW[4+0])+abs(WW[4+0]-WWW[4+0])+1,
					abs(W[4+1]-WW[4+1])+abs(WW[4+1]-WWW[4+1])+1,
					abs(W[4+2]-WW[4+2])+abs(WW[4+2]-WWW[4+2])+1,
				};
				int ygrads[]=
				{
					abs(N[4+0]-NN[4+0])+abs(NN[4+0]-NNN[4+0])+1,
					abs(N[4+1]-NN[4+1])+abs(NN[4+1]-NNN[4+1])+1,
					abs(N[4+2]-NN[4+2])+abs(NN[4+2]-NNN[4+2])+1,
				};
				nbypass[0]=((W[4+0]+WW[4+0]+WWW[4+0])*ygrads[0]+(N[4+0]+NN[4+0]+NNN[4+0])*xgrads[0])/(xgrads[0]+ygrads[0]);
				nbypass[1]=((W[4+1]+WW[4+1]+WWW[4+1])*ygrads[1]+(N[4+1]+NN[4+1]+NNN[4+1])*xgrads[1])/(xgrads[1]+ygrads[1]);
				nbypass[2]=((W[4+2]+WW[4+2]+WWW[4+2])*ygrads[2]+(N[4+2]+NN[4+2]+NNN[4+2])*xgrads[2])/(xgrads[2]+ygrads[2]);
			}
#endif
#endif
	#define GR_PREDICT(DST, IDX) DST=(7*W[3*2+IDX]+3*(NE[3*1+IDX]+N[3*3+IDX]+W[3*3+IDX]))>>4, DST+=(DST)<4
	#define GR_UPDATE1(IDX) (MAXVAR(NW[IDX], W[IDX])+sym+NEE[IDX]+MAXVAR(WW[IDX], WWW[IDX]))>>2	//for SW (thru NE)
	#define GR_UPDATE2(IDX) (MAXVAR(WW[IDX], W[IDX])+sym+NE[IDX]+MAXVAR(NEE[IDX], NEEE[IDX]))>>2	//for E (thru W)
	#define GR_UPDATE3(IDX) (W[IDX]+sym+MAXVAR(N[IDX], NE[IDX]))/3
			
//	#define GR_PREDICT(DST, IDX) DST=(N[3*3+IDX]+W[3*3+IDX])>>1, DST+=(DST)<4
//	#define GR_UPDATE1(IDX) 0
//	#define GR_UPDATE2(IDX) 0
//	#define GR_UPDATE3(IDX) (W[IDX]+sym+MAXVAR(N[IDX], NE[IDX]))/3

//	#define GR_PREDICT(DST, IDX) DST=(7*W[3*2+IDX]+3*(NE[3*1+IDX]+N[3*3+IDX]+(IDX==0?W[3*3+IDX]:WW[3*3+IDX])))>>4, DST+=(DST)<4
//	#define GR_UPDATE3(IDX) (MAXVAR(NE[IDX], N[IDX])+(IDX==9?W[IDX]:WW[IDX])+sym)/3
//	#define UPDATE_FORMULA3(IDX) (N[IDX]+(IDX==9?W[IDX]:WW[IDX])+sym)/3
//	#define UPDATE_FORMULA3(IDX) sym
//	#define UPDATE_FORMULA3(IDX) (MAXVAR(N[IDX], W[IDX])+sym)>>1
//	#define UPDATE_FORMULA(IDX) (MAXVAR(WW[IDX], W[IDX])+sym+NE[IDX]+MAXVAR(NEE[IDX], NEEE[IDX]))>>2	//Formula12	best
//	#define UPDATE_FORMULA(IDX) (max3(WW[IDX], W[IDX], NW[IDX])+sym+max3(NNE[IDX], NE[IDX], NEE[IDX])+NEEE[IDX])>>2
//	#define UPDATE_FORMULA(IDX) (W[IDX]+sym+NE[IDX]+NEEE[IDX])>>2	//Formula8
//	#define UPDATE_FORMULA(IDX) (MAXVAR(N[IDX], W[IDX])+sym+NE[IDX]+MAXVAR(NEE[IDX], NEEE[IDX]))>>2
//	#define UPDATE_FORMULA(IDX) (max3(WW[IDX], W[IDX], N[IDX])+sym+NE[IDX]+MAXVAR(NEE[IDX], NEEE[IDX]))>>2
//	#define UPDATE_FORMULA(IDX) (sum_largest3of7(WW[IDX], W[IDX], NW[IDX], N[IDX], NE[IDX], NEE[IDX], NEEE[IDX])+sym)>>2
//	#define UPDATE_FORMULA(IDX) (MAXVAR(WW[IDX], W[IDX])+sym+NE[IDX]+NNE[IDX]+NNNE[IDX]+MAXVAR(NEE[IDX], NEEE[IDX]))/6
//	#define UPDATE_FORMULA(IDX) (MAXVAR(WW[IDX], W[IDX])+sym+NE[IDX]+NNE[IDX]/2+MAXVAR(NEE[IDX], NEEE[IDX]))*2/9
//	#define UPDATE_FORMULA(IDX) (MAXVAR(N[IDX], W[IDX])+sym+NE[IDX]+NEEE[IDX])>>2		//Formula11
//	#define UPDATE_FORMULA(IDX) (MAXVAR(N[IDX], W[IDX])+sym+NE[IDX]+NEE[IDX])>>2		//X
//	#define UPDATE_FORMULA(IDX) (MAXVAR(N[IDX], W[IDX])+2*(sym+NE[IDX])+NEE[IDX]+NEEE[IDX])/7	//Formula10
//	#define UPDATE_FORMULA(IDX) (NE[IDX]>sym?sym+2*NE[IDX]+NEE[IDX]:W[IDX]+2*sym+NE[IDX])>>2//X
//	#define UPDATE_FORMULA(IDX) (sym+NE[IDX]+NEE[IDX]+NEEE[IDX])>>2
//	#define UPDATE_FORMULA(IDX) (2*W[IDX]+sym+NEEE[IDX])>>2		//Formula1
//	#define UPDATE_FORMULA(IDX) (W[IDX]+sym+NEEE[IDX])/3
//	#define UPDATE_FORMULA(IDX) (WWW[IDX]+5*WW[IDX]+10*(W[IDX]+sym)+NE[IDX])/27	//Formula5	X
//	#define UPDATE_FORMULA(IDX) sym
#if 0
			if(!args->use_AC[0])
			{
				nbypass[0]=(7*W[3*2+0]+3*(NE[3*1+0]+N[3*3+0]+W [3*3+0]))>>4;
				nbypass[0]+=nbypass[0]<4;
			//	nbypass[0]+=(nbypass[0]<4)|((W[3*3+0]&N[3*3+0]&1));
#ifndef ENABLE_NPOT
				nbypass[0]=FLOOR_LOG2(nbypass[0]);
#endif
			}
			if(!args->use_AC[1])
			{
				nbypass[1]=(7*W[3*2+1]+3*(NE[3*1+1]+N[3*3+1]+WW[3*3+1]))>>4;
				nbypass[1]+=nbypass[1]<4;
			//	nbypass[1]+=(nbypass[1]<4)|((W[3*3+1]&N[3*3+1]&1));
#ifndef ENABLE_NPOT
				nbypass[1]=FLOOR_LOG2(nbypass[1]);
#endif
			}
			if(!args->use_AC[2])
			{
				nbypass[2]=(7*W[3*2+2]+3*(NE[3*1+2]+N[3*3+2]+WW[3*3+2]))>>4;
				nbypass[2]+=nbypass[2]<4;
			//	nbypass[2]+=(nbypass[2]<4)|((W[3*3+2]&N[3*3+2]&1));
#ifndef ENABLE_NPOT
				nbypass[2]=FLOOR_LOG2(nbypass[2]);
#endif
			}
			//if(nbypass[0]>=7||nbypass[1]>=7||nbypass[2]>=7)//not hit
			//	printf("");
#endif
#if 1
			{
				ALIGN(16) short vmin1[8];
				ALIGN(16) short vmax1[8];
				ALIGN(16) short vmin2[8];
				ALIGN(16) short vmax2[8];
				__m128i mN	=_mm_load_si128((__m128i*)N);
				__m128i mNE	=_mm_load_si128((__m128i*)NE);
				__m128i mW	=_mm_load_si128((__m128i*)W);
				__m128i mmin=_mm_min_epi16(mN, mW);
				__m128i mmax=_mm_max_epi16(mN, mW);
				_mm_store_si128((__m128i*)vmin1, mmin);
				_mm_store_si128((__m128i*)vmax1, mmax);
				mmin=_mm_min_epi16(mmin, mNE);
				mmax=_mm_max_epi16(mmax, mNE);
				_mm_store_si128((__m128i*)vmin2, mmin);
				_mm_store_si128((__m128i*)vmax2, mmax);
#ifdef __GNUC__
#pragma GCC unroll 3
#elif defined __clang__
#pragma clang loop unroll_count(3)
#endif
				for(int kc=0;kc<3;++kc)//linear predictor (better for smooth areas)
				{
					switch(predidx[kc])
					{
					case PRED_N:
						preds[kc]=N[kc];
						break;
					case PRED_W:
						preds[kc]=W[kc];
						break;
					case PRED_AV2:
						preds[kc]=(N[kc]+W[kc]+1)>>1;
						break;
#if 1
					case PRED_AV3:
						preds[kc]=(3*(N[kc]+W[kc])-2*NW[kc]+2)>>2;
						break;
					case PRED_AV4:
						preds[kc]=(4*(N[kc]+W[kc])+NE[kc]-NW[kc]+4)>>3;
						break;
					case PRED_AV5:
						preds[kc]=W[kc]+((5*(N[kc]-NW[kc])+NE[kc]-WW[kc]+4)>>3);
						break;
					case PRED_AV6:
						preds[kc]=W[kc]+((6*N[kc]-5*NW[kc]-NN[kc]-WW[kc]+NE[kc]+4)>>3);
						break;
					case PRED_AV9:
						preds[kc]=W[kc]+((10*N[kc]-9*NW[kc]+4*NE[kc]-2*(NN[kc]+WW[kc])+NNW[kc]-NNE[kc]-NWW[kc]+8)>>4);
						break;
					case PRED_AV12:
						preds[kc]=(
							+av12_icoeffs[ 0]*NNWW[kc]
							+av12_icoeffs[ 1]*NNW[kc]
							+av12_icoeffs[ 2]*NN[kc]
							+av12_icoeffs[ 3]*NNE[kc]
							+av12_icoeffs[ 4]*NNEE[kc]
							+av12_icoeffs[ 5]*NWW[kc]
							+av12_icoeffs[ 6]*NW[kc]
							+av12_icoeffs[ 7]*N[kc]
							+av12_icoeffs[ 8]*NE[kc]
							+av12_icoeffs[ 9]*NEE[kc]
							+av12_icoeffs[10]*WW[kc]
							+av12_icoeffs[11]*W[kc]
							+128
						)>>8;
						break;
					case PRED_CG:
						vmin2[kc]=vmin1[kc];
						vmax2[kc]=vmax1[kc];
						preds[kc]=N[kc]+W[kc]-NW[kc];
						break;
#endif
#if 0
					case PRED_AV5:
					//		-5	5	1
					//	-1	8	[?]>>3
						CLAMP3_32(preds[kc], W[kc]+((5*(N[kc]-NW[kc])+NE[kc]-WW[kc]+4)>>3), N[kc], W[kc], NE[kc]);
						break;
					case PRED_AV9:
					//		1	-2	-1
					//	-1	-9	10	4
					//	-2	16	[?]>>4
						CLAMP3_32(preds[kc], W[kc]+((10*N[kc]-9*NW[kc]+4*NE[kc]-2*(NN[kc]+WW[kc])+NNW[kc]-NNE[kc]-NWW[kc]+8)>>4), N[kc], W[kc], NE[kc]);
						break;
					case PRED_AV12:
						preds[kc]=(
							+av12_icoeffs[ 0]*NNWW[kc]
							+av12_icoeffs[ 1]*NNW[kc]
							+av12_icoeffs[ 2]*NN[kc]
							+av12_icoeffs[ 3]*NNE[kc]
							+av12_icoeffs[ 4]*NNEE[kc]
							+av12_icoeffs[ 5]*NWW[kc]
							+av12_icoeffs[ 6]*NW[kc]
							+av12_icoeffs[ 7]*N[kc]
							+av12_icoeffs[ 8]*NE[kc]
							+av12_icoeffs[ 9]*NEE[kc]
							+av12_icoeffs[10]*WW[kc]
							+av12_icoeffs[11]*W[kc]
							+128
						)>>8;
						CLAMP3_32(preds[kc], preds[kc], N[kc], W[kc], NE[kc]);
						break;
					case PRED_CG:
						MEDIAN3_32(preds[kc], N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
						break;
#endif
					}
				}
				__m128i mp=_mm_load_si128((__m128i*)preds);
				mp=_mm_max_epi16(mp, _mm_load_si128((__m128i*)vmin2));
				mp=_mm_min_epi16(mp, _mm_load_si128((__m128i*)vmax2));

				//__m128i mc=_mm_or_si128(_mm_sub_epi16(mN, mW), _mm_sub_epi16(mN, mNE));//if(N==W&&N==NE)pred=W	does nothing
				//mc=_mm_cmpeq_epi16(mc, _mm_setzero_si128());
				//mp=_mm_blendv_epi8(mp, mW, mc);

				_mm_store_si128((__m128i*)preds, mp);
			}
#endif
#ifdef USE_ANGULAR_PRED
			{//angular predictor (better for high-noise/edges) weak improvement, is it worth 1 IDIV/px? (15% slower)
				int grads[]=
				{
					abs(NE[0]-NNE[0])+abs(N[0]-NN[0])+abs(NW[0]-NNW[0])+abs(W[0]-NW[0])+1,//N
					abs(NE[1]-NNE[1])+abs(N[1]-NN[1])+abs(NW[1]-NNW[1])+abs(W[1]-NW[1])+1,
					abs(NE[2]-NNE[2])+abs(N[2]-NN[2])+abs(NW[2]-NNW[2])+abs(W[2]-NW[2])+1,

					abs(NE[0]-N[0])+abs(N[0]-NW[0])+abs(NW[0]-NWW[0])+abs(W[0]-WW[0])+1,//W
					abs(NE[1]-N[1])+abs(N[1]-NW[1])+abs(NW[1]-NWW[1])+abs(W[1]-WW[1])+1,
					abs(NE[2]-N[2])+abs(N[2]-NW[2])+abs(NW[2]-NWW[2])+abs(W[2]-WW[2])+1,
				};
				//linear[0]=preds[0];
				//linear[1]=preds[1];
				//linear[2]=preds[2];
				angular[0]=(grads[3*1+0]*N[0]+grads[3*0+0]*W[0])/(grads[3*0+0]+grads[3*1+0]);
				angular[1]=(grads[3*1+1]*N[1]+grads[3*0+1]*W[1])/(grads[3*0+1]+grads[3*1+1]);
				angular[2]=(grads[3*1+2]*N[2]+grads[3*0+2]*W[2])/(grads[3*0+2]+grads[3*1+2]);

				//if(mixer[0]!=(1<<MIXBITS>>1)&&abs(linear[0]-angular[0])>8)//
				//	printf("");
				//preds[0]+=(angular[0]-preds[0])*(mixer[0]+(1<<MIXBITS>>1))>>MIXBITS;
				//preds[1]+=(angular[1]-preds[1])*(mixer[1]+(1<<MIXBITS>>1))>>MIXBITS;
				//preds[2]+=(angular[2]-preds[2])*(mixer[2]+(1<<MIXBITS>>1))>>MIXBITS;

				int noise[]=
				{
					(abs(W[0]-NW[0])+abs(N[0]-NW[0])+abs(N[0]-NE[0])+abs(N[0]-NN[0])+abs(W[0]-WW[0])+abs(NE[0]-NEE[0])+abs(NEE[0]-NEEE[0]))*170>>8,
					(abs(W[1]-NW[1])+abs(N[1]-NW[1])+abs(N[1]-NE[1])+abs(N[1]-NN[1])+abs(W[1]-WW[1])+abs(NE[1]-NEE[1])+abs(NEE[1]-NEEE[1]))*170>>8,
					(abs(W[2]-NW[2])+abs(N[2]-NW[2])+abs(N[2]-NE[2])+abs(N[2]-NN[2])+abs(W[2]-WW[2])+abs(NE[2]-NEE[2])+abs(NEE[2]-NEEE[2]))*170>>8,
				};
#define ANGULAR_STRENGTH 8
				if(noise[0]>(1<<ANGULAR_STRENGTH))noise[0]=1<<ANGULAR_STRENGTH;
				if(noise[1]>(1<<ANGULAR_STRENGTH))noise[1]=1<<ANGULAR_STRENGTH;
				if(noise[2]>(1<<ANGULAR_STRENGTH))noise[2]=1<<ANGULAR_STRENGTH;
				preds[0]+=((angular[0]-preds[0])*noise[0]+(1<<ANGULAR_STRENGTH>>1))>>ANGULAR_STRENGTH;
				preds[1]+=((angular[1]-preds[1])*noise[1]+(1<<ANGULAR_STRENGTH>>1))>>ANGULAR_STRENGTH;
				preds[2]+=((angular[2]-preds[2])*noise[2]+(1<<ANGULAR_STRENGTH>>1))>>ANGULAR_STRENGTH;
				//preds[0]=angular[0];
				//preds[1]=angular[1];
				//preds[2]=angular[2];
				//MEDIAN3_32(preds[0], N[0], W[0], N[0]+W[0]-NW[0]);
				//MEDIAN3_32(preds[1], N[1], W[1], N[1]+W[1]-NW[1]);
				//MEDIAN3_32(preds[2], N[2], W[2], N[2]+W[2]-NW[2]);

				//preds[0]=(15*preds[0]+angular[0]+8)>>4;
				//preds[1]=(15*preds[1]+angular[1]+8)>>4;
				//preds[2]=(15*preds[2]+angular[2]+8)>>4;
			}
#endif
#if 0
			{//angular predictor (better for high-noise/edges) bad
				int grads[]=
				{
					abs(NE[0]-NNEE[0])+abs(N[0]-NNE[0])+abs(NW[0]-NN[0])+abs(W[0]-N[0]),//NE
					abs(NE[1]-NNEE[1])+abs(N[1]-NNE[1])+abs(NW[1]-NN[1])+abs(W[1]-N[1]),
					abs(NE[2]-NNEE[2])+abs(N[2]-NNE[2])+abs(NW[2]-NN[2])+abs(W[2]-N[2]),

					abs(NE[0]-NNE[0])+abs(N[0]-NN[0])+abs(NW[0]-NNW[0])+abs(W[0]-NW[0]),//N
					abs(NE[1]-NNE[1])+abs(N[1]-NN[1])+abs(NW[1]-NNW[1])+abs(W[1]-NW[1]),
					abs(NE[2]-NNE[2])+abs(N[2]-NN[2])+abs(NW[2]-NNW[2])+abs(W[2]-NW[2]),
				
					abs(NE[0]-NN[0])+abs(N[0]-NNW[0])+abs(NW[0]-NNWW[0])+abs(W[0]-NWW[0]),//NW
					abs(NE[1]-NN[1])+abs(N[1]-NNW[1])+abs(NW[1]-NNWW[1])+abs(W[1]-NWW[1]),
					abs(NE[2]-NN[2])+abs(N[2]-NNW[2])+abs(NW[2]-NNWW[2])+abs(W[2]-NWW[2]),

					abs(NE[0]-N[0])+abs(N[0]-NW[0])+abs(NW[0]-NWW[0])+abs(W[0]-WW[0]),//W
					abs(NE[1]-N[1])+abs(N[1]-NW[1])+abs(NW[1]-NWW[1])+abs(W[1]-WW[1]),
					abs(NE[2]-N[2])+abs(N[2]-NW[2])+abs(NW[2]-NWW[2])+abs(W[2]-WW[2]),
				};
				long long weights[]=
				{
					(long long)grads[3*1+0]*grads[3*2+0]*grads[3*3+0]+1,
					(long long)grads[3*0+0]*grads[3*2+0]*grads[3*3+0]+1,
					(long long)grads[3*0+0]*grads[3*1+0]*grads[3*3+0]+1,
					(long long)grads[3*0+0]*grads[3*1+0]*grads[3*2+0]+1,

					(long long)grads[3*1+1]*grads[3*2+1]*grads[3*3+1]+1,
					(long long)grads[3*0+1]*grads[3*2+1]*grads[3*3+1]+1,
					(long long)grads[3*0+1]*grads[3*1+1]*grads[3*3+1]+1,
					(long long)grads[3*0+1]*grads[3*1+1]*grads[3*2+1]+1,

					(long long)grads[3*1+2]*grads[3*2+2]*grads[3*3+2]+1,
					(long long)grads[3*0+2]*grads[3*2+2]*grads[3*3+2]+1,
					(long long)grads[3*0+2]*grads[3*1+2]*grads[3*3+2]+1,
					(long long)grads[3*0+2]*grads[3*1+2]*grads[3*2+2]+1,
				};
				//linear[0]=preds[0];
				//linear[1]=preds[1];
				//linear[2]=preds[2];
				angular[0]=(weights[4*0+0]*NE[0]+weights[4*0+1]*N[0]+weights[4*0+2]*NW[0]+weights[4*0+3]*W[0])/(weights[4*0+0]+weights[4*0+1]+weights[4*0+2]+weights[4*0+3]);
				angular[1]=(weights[4*1+0]*NE[1]+weights[4*1+1]*N[0]+weights[4*1+2]*NW[1]+weights[4*1+3]*W[0])/(weights[4*1+0]+weights[4*1+1]+weights[4*1+2]+weights[4*1+3]);
				angular[2]=(weights[4*2+0]*NE[2]+weights[4*2+1]*N[0]+weights[4*2+2]*NW[2]+weights[4*2+3]*W[0])/(weights[4*2+0]+weights[4*2+1]+weights[4*2+2]+weights[4*2+3]);

				//if(mixer[0]!=(1<<MIXBITS>>1)&&abs(linear[0]-angular[0])>8)//
				//	printf("");
				//preds[0]+=(angular[0]-preds[0])*(mixer[0]+(1<<MIXBITS>>1))>>MIXBITS;
				//preds[1]+=(angular[1]-preds[1])*(mixer[1]+(1<<MIXBITS>>1))>>MIXBITS;
				//preds[2]+=(angular[2]-preds[2])*(mixer[2]+(1<<MIXBITS>>1))>>MIXBITS;

				int noise[]=
				{
					abs(N[0]-W[0])+abs(W[0]-NW[0])+abs(N[0]-NW[0])+abs(N[0]-NE[0]),
					abs(N[1]-W[1])+abs(W[1]-NW[1])+abs(N[1]-NW[1])+abs(N[1]-NE[1]),
					abs(N[2]-W[2])+abs(W[2]-NW[2])+abs(N[2]-NW[2])+abs(N[2]-NE[2]),
				};
				//if(noise[0]>256)noise[0]=256;
				//if(noise[1]>256)noise[1]=256;
				//if(noise[2]>256)noise[2]=256;
				preds[0]+=((angular[0]-preds[0])*noise[0]+(1<<10>>1))>>10;
				preds[1]+=((angular[1]-preds[1])*noise[1]+(1<<10>>1))>>10;
				preds[2]+=((angular[2]-preds[2])*noise[2]+(1<<10>>1))>>10;
				//preds[0]=angular[0];
				//preds[1]=angular[1];
				//preds[2]=angular[2];
				//MEDIAN3_32(preds[0], N[0], W[0], N[0]+W[0]-NW[0]);
				//MEDIAN3_32(preds[1], N[1], W[1], N[1]+W[1]-NW[1]);
				//MEDIAN3_32(preds[2], N[2], W[2], N[2]+W[2]-NW[2]);

				//preds[0]=(15*preds[0]+angular[0]+8)>>4;
				//preds[1]=(15*preds[1]+angular[1]+8)>>4;
				//preds[2]=(15*preds[2]+angular[2]+8)>>4;
			}
#endif
#if 0
			ALIGN(16) short preds[8];
			__m128i mNW	=_mm_loadu_si128((__m128i*)NW);
			__m128i mN	=_mm_loadu_si128((__m128i*)N);
			__m128i mNE	=_mm_loadu_si128((__m128i*)NE);
			__m128i mW	=_mm_loadu_si128((__m128i*)W);
			__m128i mmin=_mm_min_epi16(mN, mW);
			__m128i mmax=_mm_max_epi16(mN, mW);
			__m128i mp=_mm_add_epi16(mN, mW);
			mp=_mm_sub_epi16(mp, mNW);

			__m128i mav2=_mm_blendv_epi8(mN, mW, predmask1);
		//	__m128i mav4=_mm_srai_epi16(_mm_add_epi16(mW, mNE), 1);
			mav2=_mm_blendv_epi8(mav2, mp, predmask2);

		//	__m128i mav2=_mm_add_epi16(mp, mNW);
		//	mav2=_mm_add_epi16(mav2, mNE);
		//	mav2=_mm_srai_epi16(mav2, 2);

		//	__m128i mav2=_mm_srai_epi16(mp, 1);
			mp=_mm_max_epi16(mp, mmin);
			mp=_mm_min_epi16(mp, mmax);
			mp=_mm_blendv_epi8(mav2, mp, predmask3);
			_mm_store_si128((__m128i*)preds, mp);
#endif
#ifdef ENABLE_SSE
#if 1
			int den[3];
			int *curr_sse1[]=
			{
				args->sse1[0][((N[0]+W[0]+(NE[0]-NW[0])/2)>>3)&63][(preds[0]>>2)&63],
				args->sse1[1][((N[1]+W[1]+(NE[1]-NW[1])/2)>>3)&63][(preds[1]>>2)&63],
				args->sse1[2][((N[2]+W[2]+(NE[2]-NW[2])/2)>>3)&63][(preds[2]>>2)&63],
			};
			den[0]=curr_sse1[0][0]+5;
			den[1]=curr_sse1[1][0]+5;
			den[2]=curr_sse1[2][0]+5;
			preds[0]+=(curr_sse1[0][1]+(den[0]>>1))/den[0];
			preds[1]+=(curr_sse1[1][1]+(den[1]>>1))/den[1];
			preds[2]+=(curr_sse1[2][1]+(den[2]>>1))/den[2];
#endif
#if 1
			int *curr_sse2[]=
			{
				args->sse2[0][((N[0]-NW[0])>>1)&63][((W[0]-WW[0])>>2)&63],
				args->sse2[1][((N[1]-NW[1])>>1)&63][((W[1]-WW[1])>>2)&63],
				args->sse2[2][((N[2]-NW[2])>>1)&63][((W[2]-WW[2])>>2)&63],
			};
			den[0]=curr_sse2[0][0]+5;
			den[1]=curr_sse2[1][0]+5;
			den[2]=curr_sse2[2][0]+5;
			preds[0]+=(curr_sse2[0][1]+(den[0]>>1))/den[0];
			preds[1]+=(curr_sse2[1][1]+(den[1]>>1))/den[1];
			preds[2]+=(curr_sse2[2][1]+(den[2]>>1))/den[2];
			int *curr_sse3[]=
			{
				args->sse3[0][((W[0]-NW[0])>>1)&63][((N[0]-NN[0])>>2)&63],
				args->sse3[1][((W[1]-NW[1])>>1)&63][((N[1]-NN[1])>>2)&63],
				args->sse3[2][((W[2]-NW[2])>>1)&63][((N[2]-NN[2])>>2)&63],
			};
			den[0]=curr_sse3[0][0]+5;
			den[1]=curr_sse3[1][0]+5;
			den[2]=curr_sse3[2][0]+5;
			preds[0]+=(curr_sse3[0][1]+(den[0]>>1))/den[0];
			preds[1]+=(curr_sse3[1][1]+(den[1]>>1))/den[1];
			preds[2]+=(curr_sse3[2][1]+(den[2]>>1))/den[2];
#endif
#if 1
			int *curr_sse4[]=
			{
				args->sse4[0][preds[0]&255],
				args->sse4[1][preds[1]&255],
				args->sse4[2][preds[2]&255],
			};
			den[0]=curr_sse4[0][0]+5;
			den[1]=curr_sse4[1][0]+5;
			den[2]=curr_sse4[2][0]+5;
			preds[0]+=(curr_sse4[0][1]+(den[0]>>1))/den[0];
			preds[1]+=(curr_sse4[1][1]+(den[1]>>1))/den[1];
			preds[2]+=(curr_sse4[2][1]+(den[2]>>1))/den[2];
#endif
#endif
#ifdef ENABLE_SSE2
#define SSE2_SCALE 4
#define SSE2_DECAY 7
			//if(idx==6099)//
			//	printf("");
			curr_sse1[0]=&args->sse1[0][(N[0]+W[0]+(NE[0]-NW[0])/2+4)>>3&63][(preds[0]+2)>>2&63];
			curr_sse1[1]=&args->sse1[1][(N[1]+W[1]+(NE[1]-NW[1])/2+4)>>3&63][(preds[1]+2)>>2&63];
			curr_sse1[2]=&args->sse1[2][(N[2]+W[2]+(NE[2]-NW[2])/2+4)>>3&63][(preds[2]+2)>>2&63];
			preds[0]+=(*curr_sse1[0]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);//2.613108 -> 2.624025
			preds[1]+=(*curr_sse1[1]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
			preds[2]+=(*curr_sse1[2]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
			curr_sse2[0]=&args->sse2[0][(3*(N[0]-NN[0])+NNN[0]+1)>>1&255][(preds[0]+8)>>4&15];
			curr_sse2[1]=&args->sse2[1][(3*(N[1]-NN[1])+NNN[1]+1)>>1&255][(preds[1]+8)>>4&15];
			curr_sse2[2]=&args->sse2[2][(3*(N[2]-NN[2])+NNN[2]+1)>>1&255][(preds[2]+8)>>4&15];
			curr_sse3[0]=&args->sse3[0][(3*(W[0]-WW[0])+WWW[0]+1)>>1&255][(preds[0]+8)>>4&15];
			curr_sse3[1]=&args->sse3[1][(3*(W[1]-WW[1])+WWW[1]+1)>>1&255][(preds[1]+8)>>4&15];
			curr_sse3[2]=&args->sse3[2][(3*(W[2]-WW[2])+WWW[2]+1)>>1&255][(preds[2]+8)>>4&15];
			preds[0]+=(*curr_sse2[0]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);//2.613108 -> 2.616589
			preds[1]+=(*curr_sse2[1]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
			preds[2]+=(*curr_sse2[2]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
			preds[0]+=(*curr_sse3[0]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);//2.613108 -> 2.618468
			preds[1]+=(*curr_sse3[1]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
			preds[2]+=(*curr_sse3[2]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
			curr_sse4[0]=&args->sse4[0][(N[0]+W[0]-NW[0]+1)>>1&255][(preds[0]+8)>>4&15];
			curr_sse4[1]=&args->sse4[1][(N[1]+W[1]-NW[1]+1)>>1&255][(preds[1]+8)>>4&15];
			curr_sse4[2]=&args->sse4[2][(N[2]+W[2]-NW[2]+1)>>1&255][(preds[2]+8)>>4&15];
			preds[0]+=(*curr_sse4[0]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);//2.613108 -> 2.622749
			preds[1]+=(*curr_sse4[1]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
			preds[2]+=(*curr_sse4[2]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
			curr_sse5[0]=&args->sse5[0][(N[0]+NE[0]+4)>>3&63][(W[0]+NW[0]+4)>>3&63];
			curr_sse5[1]=&args->sse5[1][(N[1]+NE[1]+4)>>3&63][(W[1]+NW[1]+4)>>3&63];
			curr_sse5[2]=&args->sse5[2][(N[2]+NE[2]+4)>>3&63][(W[2]+NW[2]+4)>>3&63];
			preds[0]+=(*curr_sse5[0]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);//2.613108 -> 2.618170
			preds[1]+=(*curr_sse5[1]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
			preds[2]+=(*curr_sse5[2]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
			curr_sse6[0]=&args->sse6[0][preds[0]&511];
			curr_sse6[1]=&args->sse6[1][preds[1]&511];
			curr_sse6[2]=&args->sse6[2][preds[2]&511];
			preds[0]+=(*curr_sse6[0]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);//2.613108 -> 2.624963
			preds[1]+=(*curr_sse6[1]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
			preds[2]+=(*curr_sse6[2]+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
#endif
#if 0
			{
				__m128i mNW	=_mm_load_si128((__m128i*)NW);//flat area correction
				__m128i mN	=_mm_load_si128((__m128i*)N);
				__m128i mW	=_mm_load_si128((__m128i*)W);
				__m128i mp=_mm_load_si128((__m128i*)preds);
				__m128i mc=_mm_or_si128(_mm_sub_epi16(mN, mW), _mm_sub_epi16(mN, mNW));
				mc=_mm_cmpeq_epi16(mc, _mm_setzero_si128());
				mp=_mm_blendv_epi8(mp, mW, mc);
			//	mp=_mm_max_epi16(mp, _mm_set1_epi16(-128));	//worse (?!)
			//	mp=_mm_min_epi16(mp, _mm_set1_epi16(127));
				_mm_store_si128((__m128i*)preds, mp);
			}
			//if(N[0]==W[0]&&N[0]==NW[0])preds[0]=W[0];//flat area correction
			//if(N[1]==W[1]&&N[1]==NW[1])preds[1]=W[1];
			//if(N[2]==W[2]&&N[2]==NW[2])preds[2]=W[2];
			
			//int correction[]=
			//{
			//	preds[0], W[0], N[0], preds[0],
			//	preds[1], W[1], N[1], preds[1],
			//	preds[2], W[2], N[2], preds[2],
			//};
			//int condition[]=
			//{
			//	((N[0]==NN[0])&(W[0]==NW[0]))<<1|((W[0]==WW[0])&(N[0]==NW[0])),//worse
			//	((N[1]==NN[1])&(W[1]==NW[1]))<<1|((W[1]==WW[1])&(N[1]==NW[1])),
			//	((N[2]==NN[2])&(W[2]==NW[2]))<<1|((W[2]==WW[2])&(N[2]==NW[2])),
			//};
			//preds[0]=correction[4*0+condition[0]];
			//preds[1]=correction[4*1+condition[1]];
			//preds[2]=correction[4*2+condition[2]];

			//int pred0[]=
			//{
			//	preds[0],
			//	preds[1],
			//	preds[2],
			//};
			//{
			//	__m128i mp=_mm_load_si128((__m128i*)preds);	//worse
			//	mp=_mm_max_epi16(mp, _mm_set1_epi16(-128));
			//	mp=_mm_min_epi16(mp, _mm_set1_epi16(127));
			//	_mm_store_si128((__m128i*)preds, mp);
			//}
#endif
			int sym=0;
#ifndef DISABLE_RCTSEL
			int offset;
#endif
			if(args->fwd)//enc
			{
				//if(idx==8193027)//
				//if(idx==8193024)//
				//	printf("");
#ifndef DISABLE_RCTSEL
				yuv[0]=args->src[idx+combination[3+0]]-128;
				yuv[1]=args->src[idx+combination[3+1]]-128;
				yuv[2]=args->src[idx+combination[3+2]]-128;
#else
				yuv[0]=args->src[idx+1]-128;
				yuv[1]=args->src[idx+2]-128;
				yuv[2]=args->src[idx+0]-128;

				//Pei09 RCT		b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
				yuv[1]-=(87*yuv[2]+169*yuv[0]+128)>>8;
				yuv[2]-=yuv[0];
				yuv[0]+=(86*yuv[2]+29*yuv[1]+128)>>8;

				//J2K RCT
			//	yuv[1]-=yuv[0];
			//	yuv[2]-=yuv[0];
			//	yuv[0]+=(yuv[1]+yuv[2])>>2;
#endif

				//enc Y
				pred=preds[0];
				error=yuv[0]-pred;
				{
					int upred=HALF_Y-abs(pred), aval=abs(error);
					if(aval<=upred)
					{
						sym=error;
//#if defined ENABLE_SSE && defined ENABLE_SSE_SKEW
//						{
//							int negmask=-((sse_corr[0]<0)&(sym!=-HALF_Y));//sign is flipped if SSE correction was negative, to skew the signal
//							sym^=negmask;
//							sym-=negmask;
//						}
//#endif
						sym=sym<<1^sym>>31;//pack sign
					}
					else
						sym=upred+aval;//error sign is known
				}
				if(args->use_AC[0])
				{
					int den=stats0[0][257]++ + 257, cdf=0, freq=stats0[0][sym]++ + 1;
					for(int k=0;k<sym;++k)
						cdf+=stats0[0][k]+1;
					ac3_enc_update_NPOT(ec2+0, cdf, freq, den);
				}
				else
				{
					GR_PREDICT(nbypass, 0);
#ifdef ENABLE_NPOT
					gr_enc_NPOT(ec+0, sym, nbypass);
#else
					gr_enc(ec+0, sym, FLOOR_LOG2(nbypass));
#endif
					curr[3*1+0]=GR_UPDATE1(3*1+0);
					curr[3*2+0]=GR_UPDATE2(3*2+0);
					curr[3*3+0]=GR_UPDATE3(3*3+0);
#ifdef EXPORT_GR_PARAM
					args->grimage[idx+1]=nbypass-128;
#endif
				}
				curr[0]=yuv[0];
				
				//enc U
#ifndef DISABLE_RCTSEL
				offset=yuv[combination[6]];
				pred=preds[1]+offset;
				CLAMP2(pred, -HALF_U, HALF_U-1);
#else
				pred=preds[1];
#endif
				error=yuv[1]-pred;
				{
					int upred=HALF_U-abs(pred), aval=abs(error);
					if(aval<=upred)
					{
						sym=error;
//#if defined ENABLE_SSE && defined ENABLE_SSE_SKEW
//						{
//							int negmask=-((sse_corr[1]<0)&(sym!=-HALF_U));//sign is flipped if SSE correction was negative, to skew the signal
//							sym^=negmask;
//							sym-=negmask;
//						}
//#endif
						sym=sym<<1^sym>>31;//pack sign
					}
					else
						sym=upred+aval;//error sign is known
				}
				if(args->use_AC[1])
				{
					int den=stats0[1][257]++ + 257, cdf=0, freq=stats0[1][sym]++ + 1;
					for(int k=0;k<sym;++k)
						cdf+=stats0[1][k]+1;
					ac3_enc_update_NPOT(ec2+1, cdf, freq, den);
				}
				else
				{
					GR_PREDICT(nbypass, 1);
#ifdef ENABLE_NPOT
					gr_enc_NPOT(ec+1, sym, nbypass);
#else
					gr_enc(ec+1, sym, FLOOR_LOG2(nbypass));
#endif
					curr[3*1+1]=GR_UPDATE1(3*1+1);
					curr[3*2+1]=GR_UPDATE2(3*2+1);
					curr[3*3+1]=GR_UPDATE3(3*3+1);
#ifdef EXPORT_GR_PARAM
					args->grimage[idx+2]=nbypass-128;
#endif
				}
#ifndef DISABLE_RCTSEL
				curr[1]=yuv[1]-offset;
#else
				curr[1]=yuv[1];
#endif

				//enc V
#ifndef DISABLE_RCTSEL
				offset=(yuv[combination[7]]+yuv[combination[8]])>>combination[9];
				pred=preds[2]+offset;
				CLAMP2(pred, -HALF_V, HALF_V-1);
#else
				pred=preds[2];
#endif
				error=yuv[2]-pred;
				{
					int upred=HALF_V-abs(pred), aval=abs(error);
					if(aval<=upred)
					{
						sym=error;
//#if defined ENABLE_SSE && defined ENABLE_SSE_SKEW
//						{
//							int negmask=-((sse_corr[2]<0)&(sym!=-HALF_V));//sign is flipped if SSE correction was negative, to skew the signal
//							sym^=negmask;
//							sym-=negmask;
//						}
//#endif
						sym=sym<<1^sym>>31;//pack sign
					}
					else
						sym=upred+aval;//error sign is known
				}
				if(args->use_AC[2])
				{
					int den=stats0[2][257]++ + 257, cdf=0, freq=stats0[2][sym]++ + 1;
					for(int k=0;k<sym;++k)
						cdf+=stats0[2][k]+1;
					ac3_enc_update_NPOT(ec2+2, cdf, freq, den);
				}
				else
				{
					GR_PREDICT(nbypass, 2);
#ifdef ENABLE_NPOT
					gr_enc_NPOT(ec+2, sym, nbypass);
#else
					gr_enc(ec+2, sym, FLOOR_LOG2(nbypass));
#endif
					curr[3*1+2]=GR_UPDATE1(3*1+2);
					curr[3*2+2]=GR_UPDATE2(3*2+2);
					curr[3*3+2]=GR_UPDATE3(3*3+2);
#ifdef EXPORT_GR_PARAM
					args->grimage[idx+0]=nbypass-128;
#endif
				}
#ifndef DISABLE_RCTSEL
				curr[2]=yuv[2]-offset;
#else
				curr[2]=yuv[2];
#endif
			}
			else//dec
			{
				//if(idx==8193027)//
				//if(idx==8193024)//
				//	printf("");

				//dec Y
				pred=preds[0];
				if(args->use_AC[0])
				{
					int den=stats0[0][257]++ + 257, cdf=0, freq;
					int code=ac3_dec_getcdf_NPOT(ec2+0, den);
#ifdef _DEBUG
					if(code>=den)
						LOG_ERROR("code >= den at %d", idx);
#endif
					for(sym=0;;)
					{
						int cdf2;

						freq=stats0[0][sym]+1;
						cdf2=cdf+freq;
						if(cdf2>code)
							break;
						cdf=cdf2;
						++sym;
					}
					ac3_dec_update_NPOT(ec2+0, cdf, freq, den);
					++stats0[0][sym];
				}
				else
				{
					GR_PREDICT(nbypass, 0);
#ifdef ENABLE_NPOT
					sym=gr_dec_NPOT(ec+0, nbypass);
#else
					sym=gr_dec(ec+0, FLOOR_LOG2(nbypass));
#endif
					curr[3*1+0]=GR_UPDATE1(3*1+0);
					curr[3*2+0]=GR_UPDATE2(3*2+0);
					curr[3*3+0]=GR_UPDATE3(3*3+0);
				}
				{
					int upred=HALF_Y-abs(pred), negmask=0;
					if(sym<=(upred<<1))
					{
						error=sym>>1^-(sym&1);
//#if defined ENABLE_SSE && defined ENABLE_SSE_SKEW
//						negmask=-((sse_corr[0]<0)&(error!=-HALF_Y));
//#endif
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

				//dec U
#ifndef DISABLE_RCTSEL
				offset=yuv[combination[6]];
				pred=preds[1]+offset;
				CLAMP2(pred, -HALF_U, HALF_U-1);
#else
				pred=preds[1];
#endif
				if(args->use_AC[1])
				{
					int den=stats0[1][257]++ + 257, cdf=0, freq;
					int code=ac3_dec_getcdf_NPOT(ec2+1, den);
#ifdef _DEBUG
					if(code>=den)
						LOG_ERROR("code >= den at %d", idx);
#endif
					for(sym=0;;)
					{
						int cdf2;

						freq=stats0[1][sym]+1;
						cdf2=cdf+freq;
						if(cdf2>code)
							break;
						cdf=cdf2;
						++sym;
					}
					ac3_dec_update_NPOT(ec2+1, cdf, freq, den);
					++stats0[1][sym];
				}
				else
				{
					GR_PREDICT(nbypass, 1);
#ifdef ENABLE_NPOT
					sym=gr_dec_NPOT(ec+1, nbypass);
#else
					sym=gr_dec(ec+1, FLOOR_LOG2(nbypass));
#endif
					curr[3*1+1]=GR_UPDATE1(3*1+1);
					curr[3*2+1]=GR_UPDATE2(3*2+1);
					curr[3*3+1]=GR_UPDATE3(3*3+1);
				}
				{
					int upred=HALF_U-abs(pred), negmask=0;
					if(sym<=(upred<<1))
					{
						error=sym>>1^-(sym&1);
//#if defined ENABLE_SSE && defined ENABLE_SSE_SKEW
//						negmask=-((sse_corr[1]<0)&(error!=-HALF_U));
//#endif
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
#ifndef DISABLE_RCTSEL
				curr[1]=yuv[1]-offset;
#else
				curr[1]=yuv[1];
#endif

				//dec V
#ifndef DISABLE_RCTSEL
				offset=(yuv[combination[7]]+yuv[combination[8]])>>combination[9];
				pred=preds[2]+offset;
				CLAMP2(pred, -HALF_V, HALF_V-1);
#else
				pred=preds[2];
#endif
				if(args->use_AC[2])
				{
					int den=stats0[2][257]++ + 257, cdf=0, freq;
					int code=ac3_dec_getcdf_NPOT(ec2+2, den);
#ifdef _DEBUG
					if(code>=den)
						LOG_ERROR("code >= den at %d", idx);
#endif
					for(sym=0;;)
					{
						int cdf2;

						freq=stats0[2][sym]+1;
						cdf2=cdf+freq;
						if(cdf2>code)
							break;
						cdf=cdf2;
						++sym;
					}
					ac3_dec_update_NPOT(ec2+2, cdf, freq, den);
					++stats0[2][sym];
				}
				else
				{
					GR_PREDICT(nbypass, 2);
#ifdef ENABLE_NPOT
					sym=gr_dec_NPOT(ec+2, nbypass);
#else
					sym=gr_dec(ec+2, FLOOR_LOG2(nbypass));
#endif
					curr[3*1+2]=GR_UPDATE1(3*1+2);
					curr[3*2+2]=GR_UPDATE2(3*2+2);
					curr[3*3+2]=GR_UPDATE3(3*3+2);
				}
				{
					int upred=HALF_V-abs(pred), negmask=0;
					if(sym<=(upred<<1))
					{
						error=sym>>1^-(sym&1);
//#if defined ENABLE_SSE && defined ENABLE_SSE_SKEW
//						negmask=-((sse_corr[2]<0)&(error!=-HALF_V));
//#endif
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
#ifndef DISABLE_RCTSEL
				curr[2]=yuv[2]-offset;
#else
				curr[2]=yuv[2];
#endif
				
#ifndef DISABLE_RCTSEL
				args->dst[idx+combination[3+0]]=yuv[0]+128;
				args->dst[idx+combination[3+1]]=yuv[1]+128;
				args->dst[idx+combination[3+2]]=yuv[2]+128;
#else
				//Pei09 iRCT
				yuv[0]-=(86*yuv[2]+29*yuv[1]+128)>>8;
				yuv[2]+=yuv[0];
				yuv[1]+=(87*yuv[2]+169*yuv[0]+128)>>8;

				//J2K iRCT
			//	yuv[0]-=(yuv[1]+yuv[2])>>2;
			//	yuv[2]+=yuv[0];
			//	yuv[1]+=yuv[0];
				args->dst[idx+1]=yuv[0]+128;
				args->dst[idx+2]=yuv[1]+128;
				args->dst[idx+0]=yuv[2]+128;
#endif
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
			//if(linear[0]!=angular[0])
			//{
			//	int alpha=linear[0]-angular[0];
			//	alpha=(((curr[0]-angular[0])<<MIXBITS)+(alpha>>1))/alpha;
			//	CLAMP2(alpha, 0, 1<<MIXBITS);
			//	mixer[0]+=(alpha-(1<<MIXBITS>>1)-mixer[0])>>6;
			//}
			//if(linear[1]!=angular[1])
			//{
			//	int alpha=linear[1]-angular[1];
			//	alpha=(((curr[1]-angular[1])<<MIXBITS)+(alpha>>1))/alpha;
			//	CLAMP2(alpha, 0, 1<<MIXBITS);
			//	mixer[1]+=(alpha-(1<<MIXBITS>>1)-mixer[1])>>6;
			//}
			//if(linear[2]!=angular[2])
			//{
			//	int alpha=linear[2]-angular[2];
			//	alpha=(((curr[2]-angular[2])<<MIXBITS)+(alpha>>1))/alpha;
			//	CLAMP2(alpha, 0, 1<<MIXBITS);
			//	mixer[2]+=(alpha-(1<<MIXBITS>>1)-mixer[2])>>6;
			//}
#ifdef ENABLE_SSE
			++curr_sse1[0][0];
			++curr_sse1[1][0];
			++curr_sse1[2][0];
			curr_sse1[0][1]+=curr[0]-preds[0];
			curr_sse1[1][1]+=curr[1]-preds[1];
			curr_sse1[2][1]+=curr[2]-preds[2];
#if 1
			++curr_sse2[0][0];
			++curr_sse2[1][0];
			++curr_sse2[2][0];
			curr_sse2[0][1]+=curr[0]-preds[0];
			curr_sse2[1][1]+=curr[1]-preds[1];
			curr_sse2[2][1]+=curr[2]-preds[2];
#endif
#if 1
			++curr_sse3[0][0];
			++curr_sse3[1][0];
			++curr_sse3[2][0];
			curr_sse3[0][1]+=curr[0]-preds[0];
			curr_sse3[1][1]+=curr[1]-preds[1];
			curr_sse3[2][1]+=curr[2]-preds[2];
#endif
#if 1
			++curr_sse4[0][0];
			++curr_sse4[1][0];
			++curr_sse4[2][0];
			curr_sse4[0][1]+=curr[0]-preds[0];
			curr_sse4[1][1]+=curr[1]-preds[1];
			curr_sse4[2][1]+=curr[2]-preds[2];
#endif
#endif
#ifdef ENABLE_SSE2
			{
				int e=(curr[0]-preds[0])<<SSE2_SCALE;
				*curr_sse1[0]=((*curr_sse1[0]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse2[0]=((*curr_sse2[0]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse3[0]=((*curr_sse3[0]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse4[0]=((*curr_sse4[0]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse5[0]=((*curr_sse5[0]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse6[0]=((*curr_sse6[0]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
			}
			{
				int e=(curr[1]-preds[1])<<SSE2_SCALE;
				*curr_sse1[1]=((*curr_sse1[1]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse2[1]=((*curr_sse2[1]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse3[1]=((*curr_sse3[1]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse4[1]=((*curr_sse4[1]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse5[1]=((*curr_sse5[1]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse6[1]=((*curr_sse6[1]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
			}
			{
				int e=(curr[2]-preds[2])<<SSE2_SCALE;
				*curr_sse1[2]=((*curr_sse1[2]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse2[2]=((*curr_sse2[2]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse3[2]=((*curr_sse3[2]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse4[2]=((*curr_sse4[2]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse5[2]=((*curr_sse5[2]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
				*curr_sse6[2]=((*curr_sse6[2]*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
			}
#endif
#if 0
			if(args->use_AC[0]&&stats0[0][257]>=25600)
			{
				int sum=0;
				for(int ks=0;ks<257;++ks)
					sum+=stats0[0][ks]>>=1;
				stats0[0][257]=sum;
			}
			if(args->use_AC[1]&&stats0[1][257]>=25600)
			{
				int sum=0;
				for(int ks=0;ks<257;++ks)
					sum+=stats0[1][ks]>>=1;
				stats0[1][257]=sum;
			}
			if(args->use_AC[2]&&stats0[2][257]>=25600)
			{
				int sum=0;
				for(int ks=0;ks<257;++ks)
					sum+=stats0[2][ks]>>=1;
				stats0[2][257]=sum;
			}
#endif
			__m256i mr=_mm256_load_si256((__m256i*)rows);
			mr=_mm256_add_epi64(mr, _mm256_set1_epi64x(sizeof(short[3*4])));
			_mm256_store_si256((__m256i*)rows, mr);
			//rows[0]+=3*4;
			//rows[1]+=3*4;
			//rows[2]+=3*4;
			//rows[3]+=3*4;
		}
	}
	if(args->fwd)
	{
		if(args->use_AC[0])
			ac3_enc_flush(ec2+0);
		else
			gr_enc_flush(ec+0);
		if(args->use_AC[1])
			ac3_enc_flush(ec2+1);
		else
			gr_enc_flush(ec+1);
		if(args->use_AC[2])
			ac3_enc_flush(ec2+2);
		else
			gr_enc_flush(ec+2);
	}
}
int c18_codec(const char *srcfn, const char *dstfn, int nthreads0)
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
	int usize;
	size_t yuvsizes[3]={0};
#ifndef DISABLE_RCTSEL
	double esize;
	double t_analysis=0;
#endif
#ifdef EXPORT_GR_PARAM
	unsigned char *grimage=0;
#endif
	
	t0=time_sec();
	src=load_file(srcfn, 1, 3, 1);
	headersize=header_read(src->data, (int)src->count, &iw, &ih, &codec);
	image=src->data+headersize;
	imageend=src->data+src->count;
#ifdef EXPORT_GR_PARAM
	grimage=(unsigned char*)malloc(3LL*iw*ih);
	if(!grimage)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
#endif
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
#if 0
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
	nblocks=xblocks*yblocks;
	if(nthreads0)
	{
		int nthreads2=MINVAR(nblocks, ncores);
		nthreads=nthreads0;
		CLAMP2(nthreads, 1, nthreads2);
	}
	else
		nthreads=MINVAR(nblocks, ncores);
	coffset=(int)sizeof(int[3])*nblocks;
	start=0;
	memusage=0;
	argssize=nthreads*sizeof(ThreadArgs);
	args=(ThreadArgs*)malloc(argssize);
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
#ifndef DISABLE_RCTSEL
	esize=0;
#endif
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

			memcpy(&size, image+sizeof(int)*(3LL*kt+0), sizeof(int));
			start+=size/RCT_COUNT>>1;

			memcpy(&size, image+sizeof(int)*(3LL*kt+1), sizeof(int));
			start+=size>>1;

			memcpy(&size, image+sizeof(int)*(3LL*kt+2), sizeof(int));
			start+=size>>1;
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
#ifdef EXPORT_GR_PARAM
		arg->grimage=grimage;
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

					memcpy(&size, image+sizeof(int)*(((ptrdiff_t)kt+kt2)*3+0), sizeof(int));
					arg->bestrct=size%RCT_COUNT;
					arg->use_AC[0]=size/RCT_COUNT&1;
					arg->decstart[0]=image+start;
					start+=size/RCT_COUNT>>1;
					arg->decend[0]=image+start;

					memcpy(&size, image+sizeof(int)*(((ptrdiff_t)kt+kt2)*3+1), sizeof(int));
					arg->use_AC[1]=size&1;
					arg->decstart[1]=image+start;
					start+=size>>1;
					arg->decend[1]=image+start;

					memcpy(&size, image+sizeof(int)*(((ptrdiff_t)kt+kt2)*3+2), sizeof(int));
					arg->use_AC[2]=size&1;
					arg->decstart[2]=image+start;
					start+=size>>1;
					arg->decend[2]=image+start;
				}
			}
#ifdef ENABLE_MT
			if(nthreads>1)
			{
				void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
				mt_finish(ctx);
			}
			else
#endif
			{
				for(int k=0;k<nthreads2;++k)
					block_thread(args+k);
			}
			if(fwd)
			{
				for(int kt2=0;kt2<nthreads2;++kt2)
				{
					ThreadArgs *arg=args+kt2;
#ifndef DISABLE_RCTSEL
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
							size_t csize=arg->list[0].nbytes+arg->list[1].nbytes+arg->list[2].nbytes;
							printf(
								"block %4d/%4d  XY %3d %3d  %4d*%4d:  %8d->%16lf->%7zd+%7zd+%7zd = %8zd bytes %s %s %s (%+10.2lf)  %10.6lf%%  CR %10lf  %s %s %s %s\n",
								kt+kt2+1, nblocks,
								kx, ky,
								arg->y2-arg->y1,
								arg->x2-arg->x1,
								blocksize,
								arg->bestsize,
								arg->list[0].nbytes,
								arg->list[1].nbytes,
								arg->list[2].nbytes,
								csize,
								arg->use_AC[0]?"AC":"GR",
								arg->use_AC[1]?"AC":"GR",
								arg->use_AC[2]?"AC":"GR",
								csize-arg->bestsize,
								100.*csize/blocksize,
								(double)blocksize/csize,
								rct_names[arg->bestrct],
								pred_names[arg->predidx[0]],
								pred_names[arg->predidx[1]],
								pred_names[arg->predidx[2]]
							);
						}
						esize+=arg->bestsize;
						yuvsizes[0]+=arg->list[0].nbytes;
						yuvsizes[1]+=arg->list[1].nbytes;
						yuvsizes[2]+=arg->list[2].nbytes;
#ifdef ABAC_PROFILESIZE
						for(int k=0;k<ABAC_TOKEN_BITS*3;++k)
							abac_csizes[k]+=arg->abac_csizes[k];
#endif
					}
#endif
					int msg[]=
					{
						((int)arg->list[0].nbytes<<1|arg->use_AC[0])*RCT_COUNT+arg->bestrct,
						(int)arg->list[1].nbytes<<1|arg->use_AC[1],
						(int)arg->list[2].nbytes<<1|arg->use_AC[2],
					};
					memcpy(dst->data+start+sizeof(int)*(((ptrdiff_t)kt+kt2)*3+0), msg+0, sizeof(int));
					memcpy(dst->data+start+sizeof(int)*(((ptrdiff_t)kt+kt2)*3+1), msg+1, sizeof(int));
					memcpy(dst->data+start+sizeof(int)*(((ptrdiff_t)kt+kt2)*3+2), msg+2, sizeof(int));
					blist_appendtoarray(arg->list+0, &dst);
					blist_appendtoarray(arg->list+1, &dst);
					blist_appendtoarray(arg->list+2, &dst);
					blist_clear(arg->list+0);
					blist_clear(arg->list+1);
					blist_clear(arg->list+2);
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
#ifndef DISABLE_RCTSEL
				printf("Best %15.2lf (%+13.2lf) bytes\n", esize, csize-esize);
#endif
				printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
				printf("Y%11zd\n", yuvsizes[0]);
				printf("U%11zd\n", yuvsizes[1]);
				printf("V%11zd\n", yuvsizes[2]);
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
#ifdef EXPORT_GR_PARAM
	if(test)
	{
		acme_strftime(g_buf, G_BUF_SIZE-1, "%Y-%m-%d_%H%M%S.PPM");
		printf("Saving \"%s\"...\n", g_buf);
		FILE *fdst2=fopen(g_buf, "wb");
		if(!fdst2)
		{
			LOG_ERROR("Cannot open \"%s\" for writing");
			return 1;
		}
		fprintf(fdst2, "P6\n%d %d\n255\n", iw, ih);
		fwrite(grimage, 1, 3LL*iw*ih, fdst2);
		fclose(fdst2);
	}
	free(grimage);
#endif
	return 0;
}