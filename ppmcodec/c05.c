#include"codec.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT

	#define ENABLE_AV2	//good
	#define ENABLE_SSE	//good with MT
//	#define ENABLE_SSE_SKEW	//bad
//	#define DISABLE_RCTSEL	//causalRCTSEL > RCT,  disabling RCTSEL breaks SSE
//	#define ENABLE_NPOT	//bad


#define CODECNAME "C05"
#include"entropy.h"

#define BLOCKSIZE 448
#define MAXPRINTEDBLOCKS 50

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
	PRED(N)\
	PRED(W)\
	PRED(AV2)\
	PRED(AV5)\
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
	short pixels[(BLOCKSIZE+16)*4*4*3];//4 padded rows * 4 channels max * {pixel, abs(e1), abs(e2)}

	BList list;
	const unsigned char *decstart, *decend;
	
#ifndef DISABLE_RCTSEL
	int hist[OCH_COUNT*PRED_COUNT<<8];
#endif
#ifdef ENABLE_SSE
	int sse1[3][64][64][2];
	int sse2[3][64][64][2];
	//int sse[3][128][2];
#endif

	//aux
	int blockidx;
	double bestsize;
	int bestrct, predidx[3];
#ifndef DISABLE_RCTSEL
	double t_analysis;
#endif
} ThreadArgs;
#if 0
FORCEINLINE int max3(int a, int b, int c)
{
	if(a<b)
		a=b;
	if(a<c)
		a=c;
	return a;
}
FORCEINLINE int max4(int v0, int v1, int v2, int v3)
{
	if(v0<v1)
		v0=v1;
	if(v0<v2)
		v0=v2;
	if(v0<v3)
		v0=v3;
	return v0;
}
FORCEINLINE int median3(int a, int b, int c)
{
	MEDIAN3_32(a, a, b, c);
	return a;
}
FORCEINLINE int sum_largest3of7(int v0, int v1, int v2, int v3, int v4, int v5, int v6)
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
	GolombRiceCoder ec;
#ifndef DISABLE_RCTSEL
//	const unsigned char *image=args->fwd?args->src:args->dst;
	unsigned char bestrct=0, combination[6]={0}, predidx[4]={0};
#endif
//	short mixcoeff[4]={0};
	
	if(args->fwd)
	{
#ifndef DISABLE_RCTSEL
		int ystride=args->iw*3;
	//	size_t csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		unsigned char predsel[OCH_COUNT]={0};
		int res;
		__m256i av12_mcoeffs[12];
	//	unsigned char eW[OCH_COUNT*PRED_COUNT];
		
	//	memset(eW, 1, sizeof(eW));
		args->t_analysis=time_sec();
		for(int k=0;k<(int)_countof(av12_mcoeffs);++k)
			av12_mcoeffs[k]=_mm256_set1_epi16(av12_icoeffs[k]>>1);
		memset(args->hist, 0, sizeof(args->hist));
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
		res=(args->x2-args->x1-4)/5*5*(args->y2-args->y1-2);
		for(int ky=args->y1+2;ky<args->y2;++ky)//analysis loop
		{
			int kx=args->x1+2;
			const unsigned char *ptr=args->src+3*(args->iw*ky+kx);

			__m256i amin=_mm256_set1_epi16(-128);
			__m256i amax=_mm256_set1_epi16(127);
			__m256i amag=_mm256_set1_epi16(255);
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
			for(;kx<args->x2-4-2;kx+=5, ptr+=15)
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
#if 0
#define UPDATE(PREDIDX, IDX0, IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7, IDX8, IDX9, IDXA, IDXB, IDXC, IDXD, IDXE)\
	do\
	{\
		int temp;\
		pred=_mm256_sub_epi16(pred, amin);\
		pred=_mm256_and_si256(pred, amag);\
		_mm256_store_si256((__m256i*)result, pred);\
		temp=FLOOR_LOG2(eW[IDX0*PRED_COUNT+PREDIDX]+1); csizes[IDX0*PRED_COUNT+PREDIDX]+=(result[0x0]>>temp)+1+temp; eW[IDX0*PRED_COUNT+PREDIDX]=result[0x0];\
		temp=FLOOR_LOG2(eW[IDX1*PRED_COUNT+PREDIDX]+1); csizes[IDX1*PRED_COUNT+PREDIDX]+=(result[0x1]>>temp)+1+temp; eW[IDX1*PRED_COUNT+PREDIDX]=result[0x1];\
		temp=FLOOR_LOG2(eW[IDX2*PRED_COUNT+PREDIDX]+1); csizes[IDX2*PRED_COUNT+PREDIDX]+=(result[0x2]>>temp)+1+temp; eW[IDX2*PRED_COUNT+PREDIDX]=result[0x2];\
		temp=FLOOR_LOG2(eW[IDX3*PRED_COUNT+PREDIDX]+1); csizes[IDX3*PRED_COUNT+PREDIDX]+=(result[0x3]>>temp)+1+temp; eW[IDX3*PRED_COUNT+PREDIDX]=result[0x3];\
		temp=FLOOR_LOG2(eW[IDX4*PRED_COUNT+PREDIDX]+1); csizes[IDX4*PRED_COUNT+PREDIDX]+=(result[0x4]>>temp)+1+temp; eW[IDX4*PRED_COUNT+PREDIDX]=result[0x4];\
		temp=FLOOR_LOG2(eW[IDX5*PRED_COUNT+PREDIDX]+1); csizes[IDX5*PRED_COUNT+PREDIDX]+=(result[0x5]>>temp)+1+temp; eW[IDX5*PRED_COUNT+PREDIDX]=result[0x5];\
		temp=FLOOR_LOG2(eW[IDX6*PRED_COUNT+PREDIDX]+1); csizes[IDX6*PRED_COUNT+PREDIDX]+=(result[0x6]>>temp)+1+temp; eW[IDX6*PRED_COUNT+PREDIDX]=result[0x6];\
		temp=FLOOR_LOG2(eW[IDX7*PRED_COUNT+PREDIDX]+1); csizes[IDX7*PRED_COUNT+PREDIDX]+=(result[0x7]>>temp)+1+temp; eW[IDX7*PRED_COUNT+PREDIDX]=result[0x7];\
		temp=FLOOR_LOG2(eW[IDX8*PRED_COUNT+PREDIDX]+1); csizes[IDX8*PRED_COUNT+PREDIDX]+=(result[0x8]>>temp)+1+temp; eW[IDX8*PRED_COUNT+PREDIDX]=result[0x8];\
		temp=FLOOR_LOG2(eW[IDX9*PRED_COUNT+PREDIDX]+1); csizes[IDX9*PRED_COUNT+PREDIDX]+=(result[0x9]>>temp)+1+temp; eW[IDX9*PRED_COUNT+PREDIDX]=result[0x9];\
		temp=FLOOR_LOG2(eW[IDXA*PRED_COUNT+PREDIDX]+1); csizes[IDXA*PRED_COUNT+PREDIDX]+=(result[0xA]>>temp)+1+temp; eW[IDXA*PRED_COUNT+PREDIDX]=result[0xA];\
		temp=FLOOR_LOG2(eW[IDXB*PRED_COUNT+PREDIDX]+1); csizes[IDXB*PRED_COUNT+PREDIDX]+=(result[0xB]>>temp)+1+temp; eW[IDXB*PRED_COUNT+PREDIDX]=result[0xB];\
		temp=FLOOR_LOG2(eW[IDXC*PRED_COUNT+PREDIDX]+1); csizes[IDXC*PRED_COUNT+PREDIDX]+=(result[0xC]>>temp)+1+temp; eW[IDXC*PRED_COUNT+PREDIDX]=result[0xC];\
		temp=FLOOR_LOG2(eW[IDXD*PRED_COUNT+PREDIDX]+1); csizes[IDXD*PRED_COUNT+PREDIDX]+=(result[0xD]>>temp)+1+temp; eW[IDXD*PRED_COUNT+PREDIDX]=result[0xD];\
		temp=FLOOR_LOG2(eW[IDXE*PRED_COUNT+PREDIDX]+1); csizes[IDXE*PRED_COUNT+PREDIDX]+=(result[0xE]>>temp)+1+temp; eW[IDXE*PRED_COUNT+PREDIDX]=result[0xE];\
	}while(0)
#endif
#if 1
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
				pred=_mm256_srai_epi16(_mm256_add_epi16(nb2[NB_N], nb2[NB_W]), 1);

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
#if 1
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
		blist_init(&args->list);
		gr_enc_init(&ec, &args->list);
#ifndef DISABLE_RCTSEL
		gr_enc_NPOT(&ec, bestrct, RCT_COUNT);
		gr_enc_NPOT(&ec, predidx[0], PRED_COUNT);
		gr_enc_NPOT(&ec, predidx[1], PRED_COUNT);
		gr_enc_NPOT(&ec, predidx[2], PRED_COUNT);
#endif
	}
	else
	{
		gr_dec_init(&ec, args->decstart, args->decend);
#ifndef DISABLE_RCTSEL
		bestrct=gr_dec_NPOT(&ec, RCT_COUNT);
		predidx[0]=gr_dec_NPOT(&ec, PRED_COUNT);
		predidx[1]=gr_dec_NPOT(&ec, PRED_COUNT);
		predidx[2]=gr_dec_NPOT(&ec, PRED_COUNT);
#endif
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
#ifdef ENABLE_SSE
	memset(args->sse1, 0, sizeof(args->sse1));
	memset(args->sse2, 0, sizeof(args->sse2));
#endif
	memset(args->pixels, 0, sizeof(args->pixels));
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) short *rows[]=
		{
			args->pixels+((BLOCKSIZE+16LL)*((ky-0LL)&3)+8LL)*4*3,
			args->pixels+((BLOCKSIZE+16LL)*((ky-1LL)&3)+8LL)*4*3,
			args->pixels+((BLOCKSIZE+16LL)*((ky-2LL)&3)+8LL)*4*3,
			args->pixels+((BLOCKSIZE+16LL)*((ky-3LL)&3)+8LL)*4*3,
		};
		int yuv[4]={0};
		int pred=0, error=0;
#ifndef DISABLE_RCTSEL
		const unsigned char *combination=rct_combinations[bestrct];
#endif
		int preds[3]={0};
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int idx=nch*(args->iw*ky+kx);
			short
				*NNN	=rows[3]+0*4*3,
				*NNNE	=rows[3]+1*4*3,
				*NNNEE	=rows[3]+2*4*3,
				*NNWW	=rows[2]-2*4*3,
				*NNW	=rows[2]-1*4*3,
				*NN	=rows[2]+0*4*3,
				*NNE	=rows[2]+1*4*3,
				*NNEE	=rows[2]+2*4*3,
			//	*NNEEE	=rows[2]+3*4*3,
				*NWW	=rows[1]-2*4*3,
				*NW	=rows[1]-1*4*3,
				*N	=rows[1]+0*4*3,
				*NE	=rows[1]+1*4*3,
				*NEE	=rows[1]+2*4*3,
				*NEEE	=rows[1]+3*4*3,
				*NEEEE	=rows[1]+4*4*3,
			//	*WWWW	=rows[0]-4*4*3,
				*WWW	=rows[0]-3*4*3,
				*WW	=rows[0]-2*4*3,
				*W	=rows[0]-1*4*3,
				*curr	=rows[0]+0*4*3;
			(void)NNNEE;
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
				NEEE-=(kx-(args->x2-4))*4*3;
#if 0
			int nbypass[3];
#ifdef __GNUC__
#pragma GCC unroll 3
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
			int nbypass[]=
			{
				(NE[4+0]+W[8+0])>>1,
				(NE[4+1]+W[8+1])>>1,
				(NE[4+2]+W[8+2])>>1,
			};
			//if(nbypass[0]<0)nbypass[0]=0;
			//if(nbypass[1]<0)nbypass[1]=0;
			//if(nbypass[2]<0)nbypass[2]=0;
			nbypass[0]+=nbypass[0]<4;
			nbypass[1]+=nbypass[1]<4;
			nbypass[2]+=nbypass[2]<4;
#ifndef ENABLE_NPOT
			nbypass[0]=FLOOR_LOG2(nbypass[0]);
			nbypass[1]=FLOOR_LOG2(nbypass[1]);
			nbypass[2]=FLOOR_LOG2(nbypass[2]);
#endif
#if 1
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
			for(int kc=0;kc<3;++kc)
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
				}
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
			//int sse_ctx[][2]=
			//{
			//	{abs(N[0]-NW[0]), abs(W[0]-NW[0])},
			//	{abs(N[1]-NW[1]), abs(W[1]-NW[1])},
			//	{abs(N[2]-NW[2]), abs(W[2]-NW[2])},
			//};
			//sse_ctx[0][0]=FLOOR_LOG2_P1(sse_ctx[0][0]);
			//sse_ctx[0][1]=FLOOR_LOG2_P1(sse_ctx[0][1]);
			//sse_ctx[1][0]=FLOOR_LOG2_P1(sse_ctx[1][0]);
			//sse_ctx[1][1]=FLOOR_LOG2_P1(sse_ctx[1][1]);
			//sse_ctx[2][0]=FLOOR_LOG2_P1(sse_ctx[2][0]);
			//sse_ctx[2][1]=FLOOR_LOG2_P1(sse_ctx[2][1]);
#if 1
			int *curr_sse2[]=
			{
				args->sse2[0][((N[0]+W[0]+(NE[0]-NW[0])/2+HALF_Y*2)>>3)&63][(preds[0]>>2)&63],
				args->sse2[1][((N[1]+W[1]+(NE[1]-NW[1])/2+HALF_U*2)>>3)&63][(preds[1]>>2)&63],
				args->sse2[2][((N[2]+W[2]+(NE[2]-NW[2])/2+HALF_V*2)>>3)&63][(preds[2]>>2)&63],

				//args->sse[0][sse_ctx[0][0]][sse_ctx[0][1]],//worse
				//args->sse[1][sse_ctx[1][0]][sse_ctx[1][1]],
				//args->sse[2][sse_ctx[2][0]][sse_ctx[2][1]],

				//args->sse[0][(preds[0]+HALF_Y)>>1&127],//worse
				//args->sse[1][(preds[1]+HALF_U)>>1&127],
				//args->sse[2][(preds[2]+HALF_V)>>1&127],
			};
			preds[0]+=(curr_sse2[0][1]+((curr_sse2[0][0]+5)>>1))/(curr_sse2[0][0]+5);
			preds[1]+=(curr_sse2[1][1]+((curr_sse2[1][0]+5)>>1))/(curr_sse2[1][0]+5);
			preds[2]+=(curr_sse2[2][1]+((curr_sse2[2][0]+5)>>1))/(curr_sse2[2][0]+5);
#endif
#if 1
			int *curr_sse1[]=
			{
				args->sse1[0][((N[0]-NW[0])>>2)&63][((W[0]-NW[0])>>2)&63],
				args->sse1[1][((N[1]-NW[1])>>2)&63][((W[1]-NW[1])>>2)&63],
				args->sse1[2][((N[2]-NW[2])>>2)&63][((W[2]-NW[2])>>2)&63],
			};
			preds[0]+=(curr_sse1[0][1]+((curr_sse1[0][0]+5)>>1))/(curr_sse1[0][0]+5);
			preds[1]+=(curr_sse1[1][1]+((curr_sse1[1][0]+5)>>1))/(curr_sse1[1][0]+5);
			preds[2]+=(curr_sse1[2][1]+((curr_sse1[2][0]+5)>>1))/(curr_sse1[2][0]+5);
#endif
			//int sse_corr[]=
			//{
			//	curr_sse[0][1]/(curr_sse[0][0]+5)+curr_sse[3][1]/(curr_sse[3][0]+5),
			//	curr_sse[1][1]/(curr_sse[1][0]+5)+curr_sse[4][1]/(curr_sse[4][0]+5),
			//	curr_sse[2][1]/(curr_sse[2][0]+5)+curr_sse[5][1]/(curr_sse[5][0]+5),
			//};
			//preds[0]+=sse_corr[0];
			//preds[1]+=sse_corr[1];
			//preds[2]+=sse_corr[2];
#endif
			int sym=0;
#ifndef DISABLE_RCTSEL
			int offset;
#endif
	#define UPDATE_FORMULA1(IDX) (MAXVAR(NW[IDX], W[IDX])+sym+NEE[IDX]+MAXVAR(WW[IDX], WWW[IDX]))>>2	//for SW (using NE)
	#define UPDATE_FORMULA2(IDX) (MAXVAR(WW[IDX], W[IDX])+sym+NE[IDX]+MAXVAR(NEE[IDX], NEEE[IDX]))>>2	//for E (using W)
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
			if(args->fwd)
			{
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
#ifdef ENABLE_NPOT
				gr_enc_NPOT(&ec, sym, nbypass[0]);
#else
				gr_enc(&ec, sym, nbypass[0]);
#endif
				curr[0]=yuv[0];
				curr[4+0]=UPDATE_FORMULA1(4+0);
				curr[8+0]=UPDATE_FORMULA2(8+0);
				
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
#ifdef ENABLE_NPOT
				gr_enc_NPOT(&ec, sym, nbypass[1]);
#else
				gr_enc(&ec, sym, nbypass[1]);
#endif
#ifndef DISABLE_RCTSEL
				curr[1]=yuv[1]-offset;
#else
				curr[1]=yuv[1];
#endif
				curr[4+1]=UPDATE_FORMULA1(4+1);
				curr[8+1]=UPDATE_FORMULA2(8+1);

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
#ifdef ENABLE_NPOT
				gr_enc_NPOT(&ec, sym, nbypass[2]);
#else
				gr_enc(&ec, sym, nbypass[2]);
#endif
#ifndef DISABLE_RCTSEL
				curr[2]=yuv[2]-offset;
#else
				curr[2]=yuv[2];
#endif
				curr[4+2]=UPDATE_FORMULA1(4+2);
				curr[8+2]=UPDATE_FORMULA2(8+2);
			}
			else
			{
				//dec Y
				pred=preds[0];
#ifdef ENABLE_NPOT
				sym=gr_dec_NPOT(&ec, nbypass[0]);
#else
				sym=gr_dec(&ec, nbypass[0]);
#endif
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
				curr[4+0]=UPDATE_FORMULA1(4+0);
				curr[8+0]=UPDATE_FORMULA2(8+0);

				//dec U
#ifndef DISABLE_RCTSEL
				offset=yuv[combination[6]];
				pred=preds[1]+offset;
				CLAMP2(pred, -HALF_U, HALF_U-1);
#else
				pred=preds[1];
#endif
#ifdef ENABLE_NPOT
				sym=gr_dec_NPOT(&ec, nbypass[1]);
#else
				sym=gr_dec(&ec, nbypass[1]);
#endif
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
				curr[4+1]=UPDATE_FORMULA1(4+1);
				curr[8+1]=UPDATE_FORMULA2(8+1);

				//dec V
#ifndef DISABLE_RCTSEL
				offset=(yuv[combination[7]]+yuv[combination[8]])>>combination[9];
				pred=preds[2]+offset;
				CLAMP2(pred, -HALF_V, HALF_V-1);
#else
				pred=preds[2];
#endif
#ifdef ENABLE_NPOT
				sym=gr_dec_NPOT(&ec, nbypass[2]);
#else
				sym=gr_dec(&ec, nbypass[2]);
#endif
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
				curr[4+2]=UPDATE_FORMULA1(4+2);
				curr[8+2]=UPDATE_FORMULA2(8+2);
				
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
#ifdef ENABLE_SSE
#if 1
			++curr_sse2[0][0];
			++curr_sse2[1][0];
			++curr_sse2[2][0];
			curr_sse2[0][1]+=curr[0]-preds[0];
			curr_sse2[1][1]+=curr[1]-preds[1];
			curr_sse2[2][1]+=curr[2]-preds[2];
#endif
#if 1
			++curr_sse1[0][0];
			++curr_sse1[1][0];
			++curr_sse1[2][0];
			curr_sse1[0][1]+=curr[0]-preds[0];
			curr_sse1[1][1]+=curr[1]-preds[1];
			curr_sse1[2][1]+=curr[2]-preds[2];
#endif
#endif
			rows[0]+=4*3;
			rows[1]+=4*3;
			rows[2]+=4*3;
			rows[3]+=4*3;
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
	int usize;
#ifndef DISABLE_RCTSEL
	double esize;
	double t_analysis=0;
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
#endif
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
#ifndef DISABLE_RCTSEL
				printf("Best %15.2lf (%+13.2lf) bytes\n", esize, csize-esize);
#endif
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