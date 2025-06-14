static const char file[]=__FILE__;
#include"ppm.h"
#include"util.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"


//	#define ENABLE_GUIDE
#ifndef DISABLE_MT
	#define ENABLE_MT
#endif

//	#define ENABLE_SELPRED
//	#define ESTIMATE_SIZE
#ifdef _DEBUG
//	#define ENABLE_ZIPF_VIEW
#endif

	#define ENABLE_WG//good
	#define ENABLE_SSE2
	#define ENABLE_SSE3_PROB
//	#define ENABLE_GATHER//bad

//	#define NONBINARY_CODING	//impractical

#define CODECNAME "C15"
#include"entropy.h"

#define BLOCKX 512
#define BLOCKY 512
#define MAXPRINTEDBLOCKS 20

#define A2_NCTX 8	//multiple of 4

#define SSE2_SCALE 11
#define SSE2_DECAY 7

#define SSE3_SCALE 11
#define SSE3_DECAY 7


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

#ifdef ENABLE_SELPRED
#define PREDLIST\
	PRED(W)\
	PRED(CG)\
	PRED(AV5)\
	PRED(AV9)\
	PRED(AV12)
#else
#define PREDLIST\
	PRED(CG)
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

typedef enum _NBIndex
{
#ifdef ENABLE_SELPRED
	NB_NNWW,	NB_NNW,		NB_NN,		NB_NNE,		NB_NNEE,
	NB_NWW,		NB_NW,		NB_N,		NB_NE,		NB_NEE,
	NB_WW,		NB_W,		NB_curr,
#else
	NB_NW,		NB_N,
	NB_W,		NB_curr,
#endif
	NB_COUNT,
} NBIndex;
//static const short av12_icoeffs[12]=
//{
//	 0x04,	 0x03,	-0x1F,	-0x26,	 0x00,
//	 0x07,	-0x9E,	 0xDB,	 0x1E,	 0x13,
//	-0x2A,	 0xF3,
//};


//WG:

#define WG_NPREDS	8	//multiple of 4
#define WG_PREDLIST0\
	WG_PRED(340,	N)\
	WG_PRED(340,	W)\
	WG_PRED(205,	3*(N-NN)+NNN)\
	WG_PRED(205,	3*(W-WW)+WWW)\
	WG_PRED(140,	W+NE-N)\
	WG_PRED(240,	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)\
	WG_PRED(120,	N+W-NW)\
	WG_PRED(120,	N+NE-NNE)
#define WG_PREDLIST1\
	WG_PRED(330,	N)\
	WG_PRED(330,	W)\
	WG_PRED(175,	3*(N-NN)+NNN)\
	WG_PRED(175,	3*(W-WW)+WWW)\
	WG_PRED(180,	W+NE-N)\
	WG_PRED(175,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW)/3)\
	WG_PRED(130,	N+W-NW)\
	WG_PRED(150,	N+NE-NNE)
#define WG_PREDLIST2\
	WG_PRED(330,	N)\
	WG_PRED(330,	W)\
	WG_PRED(200,	3*(N-NN)+NNN)\
	WG_PRED(200,	3*(W-WW)+WWW)\
	WG_PRED(180,	W+NE-N)\
	WG_PRED(175,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW)/3)\
	WG_PRED(140,	N+W-NW)\
	WG_PRED(150,	N+NE-NNE)
#if 0
#define WG_PREDLIST0\
	WG_PRED(340,	N-eN/3)\
	WG_PRED(340,	W-eW/3)\
	WG_PRED(205,	3*(N-NN)+NNN-eN/6-eNN/6+eNNN*2/3)\
	WG_PRED(205,	3*(W-WW)+WWW-eW/6-eWW/6+eWWW*2/3)\
	WG_PRED(140,	W+NE-N+((-13*eN)>>4)+eW/4-(eW>>7))\
	WG_PRED(240,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW-(4*(eN+eW)+eNN+eWW)/2)/3)\
	WG_PRED(120,	N+W-NW+(2*(eN+eW)-eNW)/5)\
	WG_PRED(120,	N+NE-NNE+((eN+eNE+eNNE+4)>>3))
#define WG_PREDLIST1\
	WG_PRED(330,	N+(2*eN+eW)/6)\
	WG_PRED(330,	W+(2*eW+eN)/6)\
	WG_PRED(175,	3*(N-NN)+NNN+eN/6+eNN/6-eWW*2/3)\
	WG_PRED(175,	3*(W-WW)+WWW+eW/6+eWW/6-eNN*2/3)\
	WG_PRED(180,	W+NE-N-((eN+eW+31)>>5))\
	WG_PRED(175,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW+(W-N+NN-NE)/2-(eN+eW+eNN/3+eWW/3))/3)\
	WG_PRED(130,	N+W-NW+(2*(eN+eW)-eNW)/5)\
	WG_PRED(150,	N+NE-NNE+(2*eN+eNE)/10)
#define WG_PREDLIST2\
	WG_PRED(330,	N+(2*eN+eW)/6)\
	WG_PRED(330,	W+(2*eW+eN)/6)\
	WG_PRED(200,	3*(N-NN)+NNN+(eN-eWW)/3)\
	WG_PRED(200,	3*(W-WW)+WWW+(eW-eNN)/3)\
	WG_PRED(180,	W+NE-N-((5*eN+eW+31)>>5))\
	WG_PRED(175,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW+(W-N+NN-NE)/2-(eN+eW+eNN/3+eWW/3))/3)\
	WG_PRED(140,	N+W-NW+(2*(eN+eW)-eNW)/8)\
	WG_PRED(150,	N+NE-NNE+(2*eN+eNE)/10)
#endif
#if 0
#define WG_PREDLIST0\
	WG_PRED(340,	N)\
	WG_PRED(340,	W)\
	WG_PRED(190,	3*(N-NN)+NNN)\
	WG_PRED(190,	3*(W-WW)+WWW)\
	WG_PRED(140,	W+NE-N+eW/3)\
	WG_PRED(240,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW-(4*(eN+eW)+eNN+eWW)/2)/3)\
	WG_PRED(120,	N+W-NW+(2*(eN+eW)-eNW)/7)\
	WG_PRED(120,	N+NE-NNE-eNNE/3)
#define WG_PREDLIST1\
	WG_PRED(330,	N-eN/3)\
	WG_PRED(330,	W-eW/3)\
	WG_PRED(175,	3*(N-NN)+NNN+eNN/6-eWW*2/3)\
	WG_PRED(175,	3*(W-WW)+WWW+eWW/6-eNN*2/3)\
	WG_PRED(180,	W+NE-N-((eN+2*eW+7)>>3))\
	WG_PRED(175,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW+(W-N+NN-NE)/2-(eN+eW+eNN/3+eWW/3))/3)\
	WG_PRED(130,	N+W-NW+(2*(eN+eW)-eNW)/4)\
	WG_PRED(110,	N+NE-NNE)
#define WG_PREDLIST2\
	WG_PRED(330,	N-eN/3)\
	WG_PRED(330,	W-eW/3)\
	WG_PRED(200,	3*(N-NN)+NNN-eWW/3)\
	WG_PRED(200,	3*(W-WW)+WWW-eNN/3)\
	WG_PRED(180,	W+NE-N-((eN+2*eW+7)>>3))\
	WG_PRED(175,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW+(W-N+NN-NE)/2-(eN+eW+eNN/3+eWW/3))/3)\
	WG_PRED(140,	N+W-NW+(2*(eN+eW)-eNW)/8)\
	WG_PRED(110,	N+NE-NNE)
#endif
//	WG_PRED(205,	5*(W-WWWW)+10*(WWW-WW)+WWWWW+eW+eWW+eWWW)
//	WG_PRED(240,	W+((5*(N-NW)+NE-WW)>>3))
//	WG_PRED(240,	W+((10*N-9*NW+4*NE-2*(NN+WW)+NNW-NNE-NWW)>>4))
//	WG_PRED(240,	(0x04*NNWW+0x03*NNW-0x1F*NN-0x26*NNE+0x07*NWW-0x9E*NW+0xDB*N+0x1E*NE+0x13*NEE-0x2A*WW+0xF3*W)>>8)
static void wg_init(double *weights, int kc)
{
	int j=0;
	switch(kc)
	{
	case 0:
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
		WG_PREDLIST0
#undef  WG_PRED
		break;
	case 1:
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
		WG_PREDLIST1
#undef  WG_PRED
		break;
	case 2:
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
		WG_PREDLIST2
#undef  WG_PRED
		break;
	}
}
AWM_INLINE int wg_predict(
	const double *weights, short **rows, const int stride, int kc2,
	int cgrad, const int *perrors, const int *NWerrors, const int *Nerrors, const int *NEerrors, const int *NNEerrors, int *preds
)
{
	int j=0;
	short
		NNNN	=rows[0][kc2+0*stride+0],
		NNNWWWW	=rows[3][kc2-4*stride+0],
		NNNWWW	=rows[3][kc2-3*stride+0],
		NNNW	=rows[3][kc2-1*stride+0],
		NNN	=rows[3][kc2+0*stride+0],
		NNNE	=rows[3][kc2+1*stride+0],
		NNNEEE	=rows[3][kc2+3*stride+0],
		NNWWWW	=rows[2][kc2-4*stride+0],
		NNWW	=rows[2][kc2-2*stride+0],
		NNW	=rows[2][kc2-1*stride+0],
		NN	=rows[2][kc2+0*stride+0],
		NNE	=rows[2][kc2+1*stride+0],
		NNEE	=rows[2][kc2+2*stride+0],
		NNEEE	=rows[2][kc2+3*stride+0],
		NNEEEE	=rows[2][kc2+4*stride+0],
		NWWW	=rows[1][kc2-3*stride+0],
		NWW	=rows[1][kc2-2*stride+0],
		NW	=rows[1][kc2-1*stride+0],
		N	=rows[1][kc2+0*stride+0],
		NE	=rows[1][kc2+1*stride+0],
		NEE	=rows[1][kc2+2*stride+0],
		NEEE	=rows[1][kc2+3*stride+0],
		NEEEE	=rows[1][kc2+4*stride+0],
		WWWWW	=rows[0][kc2-5*stride+0],
		WWWW	=rows[0][kc2-4*stride+0],
		WWW	=rows[0][kc2-3*stride+0],
		WW	=rows[0][kc2-2*stride+0],
		W	=rows[0][kc2-1*stride+0],
		eNNN	=rows[3][kc2+0*stride+1],
		eNN	=rows[2][kc2+0*stride+1],
		eNNE	=rows[2][kc2+1*stride+1],
		eNW	=rows[1][kc2-1*stride+1],
		eN	=rows[1][kc2+0*stride+1],
		eNE	=rows[1][kc2+1*stride+1],
		eNEE	=rows[1][kc2+2*stride+1],
		eNEEE	=rows[1][kc2+3*stride+1],
		eWWWW	=rows[0][kc2-4*stride+1],
		eWWW	=rows[0][kc2-3*stride+1],
		eWW	=rows[0][kc2-2*stride+1],
		eW	=rows[0][kc2-1*stride+1];
	(void)NNNN;
	(void)NNNWWWW;
	(void)NNNWWW;
	(void)NNNW;
	(void)NNN;
	(void)NNNE;
	(void)NNNEEE;
	(void)NNWWWW;
	(void)NNWW;
	(void)NNW;
	(void)NN;
	(void)NNE;
	(void)NNEE;
	(void)NNEEE;
	(void)NNEEEE;
	(void)NWWW;
	(void)NWW;
	(void)NW;
	(void)N;
	(void)NE;
	(void)NEE;
	(void)NEEE;
	(void)NEEEE;
	(void)WWWWW;
	(void)WWWW;
	(void)WWW;
	(void)WW;
	(void)W;
	(void)eNNN;
	(void)eNN;
	(void)eNNE;
	(void)eNW;
	(void)eN;
	(void)eNE;
	(void)eNEE;
	(void)eNEEE;
	(void)eWWWW;
	(void)eWWW;
	(void)eWW;
	(void)eW;
	
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST0
#undef  WG_PRED
//	switch(kc2)
//	{
//	case 0:
//#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
//		WG_PREDLIST0
//#undef  WG_PRED
//		break;
//	case 2:
//#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
//		WG_PREDLIST1
//#undef  WG_PRED
//		break;
//	case 4:
//#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
//		WG_PREDLIST2
//#undef  WG_PRED
//		break;
//	}
	
#if 1
	__m128i one=_mm_set1_epi32(1);
	__m128i me0=_mm_load_si128((__m128i*)perrors+0);
	__m128i me1=_mm_load_si128((__m128i*)perrors+1);
//	__m128i me2=_mm_load_si128((__m128i*)perrors+2);
//	__m128i me3=_mm_load_si128((__m128i*)perrors+3);
	__m256d w0=_mm256_load_pd(weights+0*4);
	__m256d w1=_mm256_load_pd(weights+1*4);
//	__m256d w2=_mm256_load_pd(weights+2*4);
//	__m256d w3=_mm256_load_pd(weights+3*4);
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)NWerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)NWerrors+1));
//	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)NWerrors+2));
//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)NWerrors+3));
	me0=_mm_add_epi32(me0, _mm_slli_epi32(_mm_load_si128((__m128i*)Nerrors+0), 1));
	me1=_mm_add_epi32(me1, _mm_slli_epi32(_mm_load_si128((__m128i*)Nerrors+1), 1));
//	me2=_mm_add_epi32(me2, _mm_slli_epi32(_mm_load_si128((__m128i*)Nerrors+2), 1));
//	me3=_mm_add_epi32(me3, _mm_slli_epi32(_mm_load_si128((__m128i*)Nerrors+3), 1));
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)NEerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)NEerrors+1));
//	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)NEerrors+2));
//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)NEerrors+3));
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)NNEerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)NNEerrors+1));
//	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)NNEerrors+2));
//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)NNEerrors+3));
	me0=_mm_add_epi32(me0, one);
	me1=_mm_add_epi32(me1, one);
//	me2=_mm_add_epi32(me2, one);
//	me3=_mm_add_epi32(me3, one);
	__m256d de0=_mm256_cvtepi32_pd(me0);
	__m256d de1=_mm256_cvtepi32_pd(me1);
//	__m256d de2=_mm256_cvtepi32_pd(me2);
//	__m256d de3=_mm256_cvtepi32_pd(me3);
	w0=_mm256_div_pd(w0, de0);
	w1=_mm256_div_pd(w1, de1);
//	w2=_mm256_div_pd(w2, de2);
//	w3=_mm256_div_pd(w3, de3);
	__m128i pr0=_mm_load_si128((__m128i*)preds+0);
	__m128i pr1=_mm_load_si128((__m128i*)preds+1);
//	pr0=_mm_min_epi32(pr0, _mm_set1_epi32(127));
//	pr1=_mm_min_epi32(pr1, _mm_set1_epi32(127));
//	pr0=_mm_max_epi32(pr0, _mm_set1_epi32(-128));
//	pr1=_mm_max_epi32(pr1, _mm_set1_epi32(-128));
	de0=_mm256_cvtepi32_pd(pr0);
	de1=_mm256_cvtepi32_pd(pr1);
//	de0=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+0));
//	de1=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+1));
//	de2=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+2));
//	de3=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+3));
	de0=_mm256_mul_pd(de0, w0);
	de1=_mm256_mul_pd(de1, w1);
//	de2=_mm256_mul_pd(de2, w2);
//	de3=_mm256_mul_pd(de3, w3);
	w0=_mm256_add_pd(w0, w1);
//	w0=_mm256_add_pd(w0, w2);
//	w0=_mm256_add_pd(w0, w3);
	de0=_mm256_add_pd(de0, de1);
//	de0=_mm256_add_pd(de0, de2);
//	de0=_mm256_add_pd(de0, de3);
	//[num3 num2 num1 num0]
	//[den3 den2 den1 den0]
	//r = hadd(num, den) = [den3+den2 num3+num2 den1+den0 num1+num0]
	//lo=_mm256_extractf128_pd(r, 0)
	//hi=_mm256_extractf128_pd(r, 1)
	//hi+lo = [den3+den2+den1+den0 num3+num2+num1+num0]
	w0=_mm256_hadd_pd(de0, w0);
	__m128d dp=_mm_add_pd(_mm256_extractf128_pd(w0, 1), _mm256_extractf128_pd(w0, 0));
	dp=_mm_div_pd(dp, _mm_permute_pd(dp, 3));
	__m128i mp=_mm_cvtpd_epi32(dp);
	__m128i mN	=_mm_set_epi32(0, 0, 0, N);
	__m128i mW	=_mm_set_epi32(0, 0, 0, W);
	__m128i mNE	=_mm_set_epi32(0, 0, 0, NE);
	__m128i vmin=_mm_min_epi32(mN, mW);
	__m128i vmax=_mm_max_epi32(mN, mW);
	vmin=_mm_min_epi32(vmin, mNE);
	vmax=_mm_max_epi32(vmax, mNE);
	mp=_mm_max_epi32(mp, vmin);
	mp=_mm_min_epi32(mp, vmax);
	{
		ALIGN(16) int pred[4];
		_mm_store_si128((__m128i*)pred, mp);
		return pred[0];
	}
#else
	{
		double pred2=0, wsum=0;
		int pred;
#ifdef __GNUC__
#pragma GCC unroll 8
#endif
		for(int k=0;k<WG_NPREDS;++k)
		{
			double weight=(double)weights[k]/(perrors[k]+1);
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
#endif
}
AWM_INLINE void wg_update(int curr, int kc, const int *preds, int *perrors, int *Nerrors, int *Werrors, int *currerrors, int *NEerrors)
{
	static const int factors[]={97, 99, 99};
	int factor=factors[kc];
#ifdef __GNUC__
#pragma GCC unroll 8
#endif
	for(int k=0;k<WG_NPREDS;++k)
	{
		int pr=preds[k];
	//	CLAMP2(pr, -128, 127);
		int e2=abs(curr-pr)<<1;
		perrors[k]=(perrors[k]+e2)*factor>>7;
	//	currerrors[k]=(e2+NEerrors[k])>>1;
		currerrors[k]=(Werrors[k]+e2+Nerrors[k]+NEerrors[k])>>2;
	//	currerrors[k]=(2*Werrors[k]+e2+NEerrors[k])>>2;
		NEerrors[k]+=e2;
	}
}


typedef struct _ThreadArgs
{
	const unsigned char *src;
	unsigned char *dst;
	int iw, ih;

	int fwd, test, loud, b1, b2, xblocks, blocksperthread, currentblock, x1, x2, y1, y2;
	int bufsize, ebufsize, histsize;
	short *pixels;
	int *ebuf;
	int *hist;

	BList *lists;
	const int *offsets;
	const unsigned char *decsrc, *decstart, *decend;

#if A2_NCTX==4
	unsigned short stats[3][(1+32+32+256)<<8];
#elif A2_NCTX==8
	unsigned short stats[3][(1+32+32+128+512+256+256+256)<<8];
#elif A2_NCTX==12
	unsigned short stats[3][(1+32+32+128+512+256+256+256+256+256+32+32)<<8];
#elif A2_NCTX==16
	unsigned short stats[3][(1+32+32+128+512+256+256+256+256+256+32+32+256+256+256+256)<<8];
#endif
#ifdef ENABLE_SSE2
	int sse1[3][32][64];
//	int sse2[3][8][8][8][8];
//	int sse3[3][32][64];
//	int sse4[3][256];
#endif
#ifdef ENABLE_SSE3_PROB
	long long sse_prob[3][8][64][32];
#endif

	//aux
	int blockidx;
	double bestsize;
	int bestrct, predidx[3];

#ifdef ESTIMATE_SIZE
	int hist2[3][256];
#endif
#ifdef ENABLE_ZIPF_VIEW
	unsigned char *zimage;
#endif
} ThreadArgs;
#if A2_NCTX==4
#define A2CTXLIST\
	A2CTX(2, 15, 0,		0)\
	A2CTX(7, 14, 0+1,	grad)\
	A2CTX(6, 14, 0+1+32,	energy)\
	A2CTX(6, 14, 0+1+32+32,	((kc2?curr[kc2-1]:(N[kc2+1]+W[kc2+1])>>1)+128)&255)
#elif A2_NCTX==8
#define A2CTXLIST\
	A2CTX(2, 15, 0,				0)\
	A2CTX(7, 14, 0+1,			grad)\
	A2CTX(6, 14, 0+1+32,			energy)\
	A2CTX(2, 15, 0+1+32+32,			(kx-args->x1)>>2)\
	A2CTX(3, 15, 0+1+32+32+128,		((N[kc2]-NW[kc2])>>5&7)<<6|((W[kc2]-NW[kc2])>>5&7)<<3|(pred>>5&7))\
	A2CTX(5, 15, 0+1+32+32+128+512,		pred+128)\
	A2CTX(6, 14, 0+1+32+32+128+512+256,	((kc2?curr[kc2-1]:(N[kc2+1]+W[kc2+1])>>1)+128)&255)\
	A2CTX(6, 15, 0+1+32+32+128+512+256+256,	(NN[kc2+1]<0)<<7|(WW[kc2+1]<0)<<6|(NE[kc2+1]<0)<<5|(NW[kc2+1]<0)<<4|(W[kc2+1]>>6&3)<<2|(N[kc2+1]>>6&3))

//	A2CTX(2, 16, 0,				kc?(ky&1)<<1|(kx&1):0)
#elif A2_NCTX==12
#define A2CTXLIST\
	A2CTX(2, 15, 0,						0)\
	A2CTX(7, 14, 0+1,					grad)\
	A2CTX(6, 14, 0+1+32,					energy)\
	A2CTX(2, 15, 0+1+32+32,					(kx-args->x1)>>2)\
	A2CTX(3, 15, 0+1+32+32+128,				((N[kc2]-NW[kc2])>>5&7)<<6|((W[kc2]-NW[kc2])>>5&7)<<3|(pred>>5&7))\
	A2CTX(5, 15, 0+1+32+32+128+512,				pred+128)\
	A2CTX(6, 14, 0+1+32+32+128+512+256,			((kc2?curr[kc2-1]:(N[kc2+1]+W[kc2+1])>>1)+128)&255)\
	A2CTX(6, 15, 0+1+32+32+128+512+256+256,			(NN[kc2+1]<0)<<7|(WW[kc2+1]<0)<<6|(NE[kc2+1]<0)<<5|(NW[kc2+1]<0)<<4|(W[kc2+1]>>6&3)<<2|(N[kc2+1]>>6&3))\
	A2CTX(5, 15, 0+1+32+32+128+512+256+256+256,		(2*N[kc2+1]-NN[kc2+1])>>1&255)\
	A2CTX(5, 15, 0+1+32+32+128+512+256+256+256+256,		(2*W[kc2+1]-WW[kc2+1])>>1&255)\
	A2CTX(5, 15, 0+1+32+32+128+512+256+256+256+256+256,	e090)\
	A2CTX(5, 15, 0+1+32+32+128+512+256+256+256+256+256+32,	e180)
#elif A2_NCTX==16
#define A2CTXLIST\
	A2CTX(2, 15, 0,								0)\
	A2CTX(7, 14, 0+1,							grad)\
	A2CTX(6, 14, 0+1+32,							energy)\
	A2CTX(2, 15, 0+1+32+32,							(kx-args->x1)>>2)\
	A2CTX(3, 15, 0+1+32+32+128,						((N[kc2]-NW[kc2])>>5&7)<<6|((W[kc2]-NW[kc2])>>5&7)<<3|(pred>>5&7))\
	A2CTX(5, 15, 0+1+32+32+128+512,						pred+128)\
	A2CTX(6, 14, 0+1+32+32+128+512+256,					((kc2?curr[kc2-1]:(N[kc2+1]+W[kc2+1])>>1)+128)&255)\
	A2CTX(6, 15, 0+1+32+32+128+512+256+256,					(NN[kc2+1]<0)<<7|(WW[kc2+1]<0)<<6|(NE[kc2+1]<0)<<5|(NW[kc2+1]<0)<<4|(W[kc2+1]>>6&3)<<2|(N[kc2+1]>>6&3))\
	A2CTX(5, 15, 0+1+32+32+128+512+256+256+256,				(2*N[kc2+1]-NN[kc2+1])>>1&255)\
	A2CTX(5, 15, 0+1+32+32+128+512+256+256+256+256,				(2*W[kc2+1]-WW[kc2+1])>>1&255)\
	A2CTX(5, 15, 0+1+32+32+128+512+256+256+256+256+256,			e090)\
	A2CTX(5, 15, 0+1+32+32+128+512+256+256+256+256+256+32,			e180)\
	A2CTX(5, 15, 0+1+32+32+128+512+256+256+256+256+256+32+32,		N[kc2+0]&255)\
	A2CTX(5, 15, 0+1+32+32+128+512+256+256+256+256+256+32+32+256,		W[kc2+0]&255)\
	A2CTX(5, 15, 0+1+32+32+128+512+256+256+256+256+256+32+32+256+256,	N[kc2+1]&255)\
	A2CTX(5, 15, 0+1+32+32+128+512+256+256+256+256+256+32+32+256+256+256,	W[kc2+1]&255)
#endif
#define A2CTX(SH, MIXSH, OFFSET, EXPR) MIXSH,
ALIGN(32) static const long long g_mixsh[]={A2CTXLIST};
#undef  A2CTX
#define A2CTX(SH, MIXSH, OFFSET, EXPR) SH,
ALIGN(32) static const long long g_sh[]={A2CTXLIST};
#undef  A2CTX
#define A2CTX(SH, MIXSH, OFFSET, EXPR) 1LL<<SH>>1,
ALIGN(32) static const long long g_offset[]={A2CTXLIST};
#undef  A2CTX
static void block_thread(void *param)
{
	const int nch=3;
	ThreadArgs *args=(ThreadArgs*)param;
	AC3 ec;
	const unsigned char *image=args->fwd?args->src:args->dst;
	unsigned char bestrct=0;
	const unsigned char *combination=0;
	unsigned char predidx[4]={0};

	if(args->fwd)
	{
		int ystride=args->iw*3;
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		unsigned char predsel[OCH_COUNT]={0};
		int res=(args->x2-args->x1-3)/5*5*(args->y2-args->y1-2);
		
		memset(args->hist, 0, args->histsize);
		for(int ky=args->y1+2;ky<args->y2;++ky)//analysis loop
		{
			int kx=args->x1+2;
			const unsigned char *ptr=image+3*(args->iw*ky+kx);

			__m256i amin=_mm256_set1_epi16(-128);
			__m256i amax=_mm256_set1_epi16(127);
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
			for(;kx<args->x2-5;kx+=5, ptr+=15)
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
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride-1*3+0)), half8),//NW
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride+0*3+0)), half8),//N
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr+0*ystride-1*3+0)), half8),//W
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr+0*ystride+0*3+0)), half8),//curr
					};
#ifdef __GNUC__
#pragma GCC unroll 4
#endif
					for(int k=0;k<NB_COUNT;++k)
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
		combination=rct_combinations[bestrct];
		predidx[0]=predsel[combination[0]];
		predidx[1]=predsel[combination[1]];
		predidx[2]=predsel[combination[2]];
		args->bestsize+=bestsize;
		args->bestrct=bestrct;
		args->predidx[0]=predidx[0];
		args->predidx[1]=predidx[1];
		args->predidx[2]=predidx[2];
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
		blist_init(args->lists+args->currentblock);
		ac3_enc_init(&ec, args->lists+args->currentblock);
		ac3_enc_bypass_NPOT(&ec, bestrct, RCT_COUNT);
	//	ac3_enc_bypass_NPOT(&ec, predidx[0], PRED_COUNT);
	//	ac3_enc_bypass_NPOT(&ec, predidx[1], PRED_COUNT);
	//	ac3_enc_bypass_NPOT(&ec, predidx[2], PRED_COUNT);
	}
	else
	{
		ac3_dec_init(&ec, args->decstart, args->decend);
		bestrct=ac3_dec_bypass_NPOT(&ec, RCT_COUNT);
		combination=rct_combinations[bestrct];
	//	predidx[0]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
	//	predidx[1]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
	//	predidx[2]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
	}
	ALIGN(32) long long curr_sh[A2_NCTX], curr_mixsh[A2_NCTX];
	memcpy(curr_sh, g_sh, sizeof(curr_sh));
	memcpy(curr_mixsh, g_mixsh, sizeof(curr_mixsh));
#if A2_NCTX==16
	ALIGN(32) unsigned short mixer1[3][8*A2_NCTX]={0};
	ALIGN(16) int mixer2[3][8*A2_NCTX/4]={0};
	FILLMEM((unsigned short*)mixer1, 0x4000, sizeof(mixer1), sizeof(short));
	//memset(mixer1, 64, sizeof(mixer1));
	FILLMEM((int*)mixer2, 0x1000, sizeof(mixer2), sizeof(int));
	memset(args->stats, 0, sizeof(args->stats));
#else
	ALIGN(32) long long mixer[8*3*A2_NCTX]={0};
	FILLMEM(mixer, 0x1000, sizeof(mixer), sizeof(long long));
#ifdef NONBINARY_CODING
	memset(args->stats, 0, sizeof(args->stats));
#else
	FILLMEM((unsigned short*)args->stats, 0x8000, sizeof(args->stats), sizeof(short));
#endif
#endif
#ifdef ENABLE_SSE2
	memset(args->sse1, 0, sizeof(args->sse1));
//	memset(args->sse2, 0, sizeof(args->sse2));
//	memset(args->sse3, 0, sizeof(args->sse3));
//	memset(args->sse4, 0, sizeof(args->sse4));
#endif
#ifdef ENABLE_SSE3_PROB
	memset(args->sse_prob, 0, sizeof(args->sse_prob));
#endif
#ifdef ESTIMATE_SIZE
	memset(args->hist2, 0, sizeof(args->hist2));
#endif
#ifdef ENABLE_WG
	ALIGN(32) double wg_weights[WG_NPREDS*3]={0};
	ALIGN(32) int wg_perrors[WG_NPREDS*3]={0}, wg_preds[WG_NPREDS]={0};
	wg_init(wg_weights+WG_NPREDS*0, 0);
	wg_init(wg_weights+WG_NPREDS*1, 1);
	wg_init(wg_weights+WG_NPREDS*2, 2);
#endif
	memset(args->ebuf, 0, args->ebufsize);
	memset(args->pixels, 0, args->bufsize);
#ifdef ABAC_SIMPLEOVERRIDE
	unsigned short stats2[768]={0};
	FILLMEM(stats2, 0x8000, sizeof(short[768]), sizeof(short));
#endif
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(32) short *rows[]=
		{
			args->pixels+((BLOCKX+16LL)*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKX+16LL)*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKX+16LL)*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKX+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		ALIGN(32) int *erows[]=
		{
			args->ebuf+((BLOCKX+16LL)*((ky-0LL)&3)+8LL)*4*WG_NPREDS,
			args->ebuf+((BLOCKX+16LL)*((ky-1LL)&3)+8LL)*4*WG_NPREDS,
			args->ebuf+((BLOCKX+16LL)*((ky-2LL)&3)+8LL)*4*WG_NPREDS,
			args->ebuf+((BLOCKX+16LL)*((ky-3LL)&3)+8LL)*4*WG_NPREDS,
		};
		int yuv[4]={0};
		int pred=0, error=0;
		for(int kx=args->x1;kx<args->x2;++kx)
		{
#ifdef ENABLE_ZIPF_VIEW
			double psize=0;
#endif
			int idx=nch*(args->iw*ky+kx);
			short
				*NNN	=rows[3]+0*4*2,
				*NNWWW	=rows[2]-3*4*2,
				*NNWW	=rows[2]-2*4*2,
				*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NNEE	=rows[2]+2*4*2,
				*NNEEE	=rows[2]+3*4*2,
				*NWWW	=rows[1]-3*4*2,
				*NWW	=rows[1]-2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*NEE	=rows[1]+2*4*2,
				*NEEE	=rows[1]+3*4*2,
				*WWWWW	=rows[0]-5*4*2,
				*WWWW	=rows[0]-4*4*2,
				*WWW	=rows[0]-3*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
			(void)NNN;
			(void)NNWWW;
			(void)NNWW;
			(void)NNW;
			(void)NN;
			(void)NNE;
			(void)NNEE;
			(void)NNEEE;
			(void)NWWW;
			(void)NWW;
			(void)NW;
			(void)N;
			(void)NE;
			(void)NEE;
			(void)NEEE;
			(void)WWWWW;
			(void)WWWW;
			(void)WWW;
			(void)WW;
			(void)W;
#if 0
			if(ky<=args->y1+2)
			{
				if(ky<=args->y1+1)
				{
					if(ky==args->y1)
						NEEE=NEE=NE=NWWW=NWW=NW=N=W;
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
#endif
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
			for(int kc=0;kc<nch;++kc)
			{
				int kc2=kc<<1, kc3=kc*WG_NPREDS;
				int offset=(yuv[combination[kc+3]]+yuv[combination[kc+6]])>>combination[kc+9];
				int
					*eNNE	=erows[2]+kc3+1*4*WG_NPREDS,
					*eNW	=erows[1]+kc3-1*4*WG_NPREDS,
					*eN	=erows[1]+kc3+0*4*WG_NPREDS,
					*eNE	=erows[1]+kc3+1*4*WG_NPREDS,
				//	*eNEE	=erows[1]+kc3+2*4*WG_NPREDS,
				//	*eNEEE	=erows[1]+kc3+3*4*WG_NPREDS,
					*eW	=erows[0]+kc3-1*4*WG_NPREDS,
					*ecurr	=erows[0]+kc3+0*4*WG_NPREDS;
				pred=wg_predict(wg_weights+WG_NPREDS*kc, rows, 4*2, kc2, 0, wg_perrors+WG_NPREDS*kc, eNW, eN, eNE, eNNE, wg_preds);
#ifdef ENABLE_SSE2
				int *curr_sse1=&args->sse1[kc][(N[kc2]-W[kc2]+8)>>4&31][(pred+2)>>2&63];	//2.858631
			//	int *curr_sse1=&args->sse1[kc][(N[kc2]-W[kc2]+16)>>5&15][(pred+1)>>1&127];	//2.858571
			//	int *curr_sse1=&args->sse1[kc][(N[kc2]-W[kc2]+4)>>3&63][(pred+4)>>3&31];	//2.858470
			//	int *curr_sse1=&args->sse1[kc][(N[kc2]+W[kc2]+(NE[kc2]-NW[kc2])/2+4)>>3&63][(pred+4)>>3&31];	//2.858261
			//	int ssectx;
			//	CLAMP3_32(ssectx, (N[kc2]+W[kc2]+(NE[kc2]-NW[kc2])/2)/2, N[kc2], W[kc2], NE[kc2]);	//2.858192
			//	int *curr_sse1=&args->sse1[kc][(N[kc2]+W[kc2]-NW[kc2]+2)>>2&63][(pred+4)>>3&31];	//2.858071
			//	//MEDIAN3_32(ssectx, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);	//2.858031
			//	int *curr_sse1=&args->sse1[kc][(ssectx+2)>>2&63][(pred+4)>>3&31];
			//	int *curr_sse1=&args->sse1[kc][(N[kc2]+W[kc2]+(NE[kc2]-NW[kc2])/2+2)>>2&127][(pred+4)>>4&15];	//2.858019
			//	int *curr_sse1=&args->sse1[kc][(N[kc2]+W[kc2]+(NE[kc2]-NW[kc2])/2+8)>>4&31][(pred+8)>>4&63];	//2.857995
			//	int *curr_sse1=&args->sse1[kc][(N[kc2]-W[kc2]+2)>>2&127][(pred+8)>>4&15];	//2.857902
				int sseval=*curr_sse1;
				pred+=(sseval+(1<<SSE2_SCALE>>1))>>SSE2_SCALE;
			//	pred+=(*curr_sse1+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);

			//	int *curr_sse2=&args->sse2[kc][(N[kc2]-NW[kc2]+32)>>6&7][(W[kc2]-NW[kc2]+32)>>6&7][(N[kc2]-NN[kc2]+32)>>6&7][(W[kc2]-WW[kc2]+32)>>6&7];//2.858805 worse
			//	pred+=(*curr_sse2+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);

			//	int *curr_sse2=&args->sse2[kc][(N[kc2]-NW[kc2]+8)>>4&31][(pred+2)>>2&63];
			//	int *curr_sse3=&args->sse3[kc][(W[kc2]-NW[kc2]+8)>>4&31][(pred+2)>>2&63];
			//	pred+=(*curr_sse2+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
			//	pred+=(*curr_sse3+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
				//int pred1=pred;
				//int pred2=pred;
				//int pred3=pred;
				//int *curr_sse_pred=&args->sse_pred[kc][((N[kc2]-W[kc2]+8)>>4)&63][(pred+2)>>2&63];
				//if(args->fwd&&args->blockidx==0)
				//	printf("%d\t%d\t%d\n", yuv[kc], pred0, pred);

			//	int *curr_sse4=&args->sse4[kc][pred&255];
			//	pred+=(*curr_sse4+(1<<(SSE2_SCALE+SSE2_DECAY)>>1))>>(SSE2_SCALE+SSE2_DECAY);
#endif
				//int pred0=pred;
				pred+=offset;
				CLAMP2(pred, -128, 127);
#if 1
				int grad=
					+2*abs(N[kc2]-W[kc2])
					+(abs(N[kc2]-NE[kc2])+abs(N[kc2]-NW[kc2])+abs(W[kc2]-NW[kc2]))
				//	+abs(N[kc2]-NN[kc2])
				//	+abs(W[kc2]-WW[kc2])
				//	+abs(NEE[kc2]-NEEE[kc2])
				;
				grad*=grad;
				grad=FLOOR_LOG2_P1(grad);//256(2+1+1+1)=1280  FLOOR_LOG2_P1(1280^2)=0~21 -> 22 levels
				int energy=
					+4*(abs(N[kc2+1])+abs(W[kc2+1]))
					+2*(abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1])+abs(NN[kc2+1])+abs(WW[kc2+1]))
					+abs(NWW[kc2+1])+abs(NNW[kc2+1])+abs(NNE[kc2+1])+abs(NEEE[kc2+1]);
				energy*=energy;
				energy=FLOOR_LOG2_P1(energy);//256(4+4+2+2+2+2+2+1+1+1+1)=5632  FLOOR_LOG2_P1(2816^2)=0~25 -> 26 levels
#endif
#if A2_NCTX>=12
				int e090=
					+abs(NE[kc2]-NNE[kc2])
					+abs(N[kc2]-NN[kc2])
					+abs(NW[kc2]-NNW[kc2])
					+2*abs(W[kc2]-NW[kc2])
					+2*abs(N[kc2+1])
					+abs(NN[kc2+1])
				;
				int e180=
					+abs(NE[kc2]-NEE[kc2])
					+abs(N[kc2]-NE[kc2])
					+2*abs(NW[kc2]-N[kc2])
					+abs(W[kc2]-WW[kc2])
					+2*abs(W[kc2+1])
					+abs(WW[kc2+1])
				;
				e090*=e090;
				e180*=e180;
				e090=FLOOR_LOG2_P1(e090);
				e180=FLOOR_LOG2_P1(e180);
#endif
#if A2_NCTX==16
#define A2CTX(SH, MIXSH, OFFSET, EXPR) args->stats[kc]+((OFFSET+(EXPR))<<8),
				short *curr_stats[A2_NCTX]={A2CTXLIST};
#undef  A2CTX
				unsigned short *curr_mixer1=mixer1[kc];
				int *curr_mixer2=mixer2[kc];
#else
#ifdef ENABLE_GATHER
				ALIGN(32) unsigned short *curr_stats[A2_NCTX];
#define A2CTX(SH, MIXSH, OFFSET, EXPR) (OFFSET+(EXPR))<<8,
				ALIGN(32) unsigned long long ctx[A2_NCTX]={A2CTXLIST};
#undef  A2CTX
				unsigned short *statsptr=args->stats[kc];
#ifdef __GNUC__
#pragma GCC unroll 12
#endif
				for(int k=0;k<A2_NCTX;++k)
					curr_stats[k]=statsptr+ctx[k];
				__m256i midx0=_mm256_load_si256((__m256i*)ctx+0);
				__m256i midx1=_mm256_load_si256((__m256i*)ctx+1);
				__m256i midx2=_mm256_load_si256((__m256i*)ctx+2);
				__m256i mtidx=_mm256_set1_epi64x(1);
				__m256i probmask=_mm256_set1_epi64x(0xFFFF);
#else
#define A2CTX(SH, MIXSH, OFFSET, EXPR) args->stats[kc]+((OFFSET+(EXPR))<<8),
				unsigned short *curr_stats[A2_NCTX]={A2CTXLIST};
#undef  A2CTX
#endif
				long long *curr_mixer=mixer+kc*8*A2_NCTX;
#endif
#ifdef NONBINARY_CODING
				int cdf, freq;
				if(args->fwd)
				{
					error=(char)(yuv[kc]-pred);
					error=error<<1^error>>31;
					for(int k=0;;)
					{
						freq=(
							+curr_mixer[0]*curr_stats[0][k]
							+curr_mixer[1]*curr_stats[1][k]
							+curr_mixer[2]*curr_stats[2][k]
							+curr_mixer[3]*curr_stats[3][k]
							+curr_mixer[4]*curr_stats[4][k]
							+curr_mixer[5]*curr_stats[5][k]
							+curr_mixer[6]*curr_stats[6][k]
							+curr_mixer[7]*curr_stats[7][k]
						)+1;
						if(k>=error)
							break;
						++k;
						cdf+=freq;
					}
					ac3_enc_update(&ec, cdf, freq);
				}
				else
				{
					int code=ac3_dec_getcdf(&ec);
					error=0;
					for(;;)
					{
						freq=(
							+curr_mixer[0]*curr_stats[0][error]
							+curr_mixer[1]*curr_stats[1][error]
							+curr_mixer[2]*curr_stats[2][error]
							+curr_mixer[3]*curr_stats[3][error]
							+curr_mixer[4]*curr_stats[4][error]
							+curr_mixer[5]*curr_stats[5][error]
							+curr_mixer[6]*curr_stats[6][error]
							+curr_mixer[7]*curr_stats[7][error]
						)+1;
						int cdf2=cdf+freq;
						if(cdf2>=code)
							break;
						++error;
						cdf=cdf2;
					}
				}
				curr[kc2+0]=yuv[kc];
				curr[kc2+1]=yuv[kc]-pred;
#else
				if(args->fwd)
				{
					curr[kc2+0]=yuv[kc];
					error=yuv[kc]-pred;
					curr[kc2+1]=curr[kc2+0]-pred;
				}
				else
					error=0;
				int tidx=1;
#ifdef __GNUC__
#pragma GCC unroll 8
#endif
				for(int kb=7;kb>=0;--kb)
				{
					int bit;

					//2.860709	baseine
					//2.856320	0	rem 0
					//2.847950	1	rem grad
					//2.852647	2	rem energy
					//2.859236	3	rem (kx-args->x1)>>2
					//2.859232	4	rem ((N[kc2]-NW[kc2])>>5&7)<<6|((W[kc2]-NW[kc2])>>5&7)<<3|(pred>>5&7)
					//2.859212	5	rem pred+128
					//2.852711	6	rem ((kc2?curr[kc2-1]:(N[kc2+1]+W[kc2+1])>>1)+128)&255
					//2.856976	7	rem (NN[kc2+1]<0)<<7|(WW[kc2+1]<0)<<6|(NE[kc2+1]<0)<<5|(NW[kc2+1]<0)<<4|(W[kc2+1]>>6&3)<<2|(N[kc2+1]>>6&3)
					//curr_mixer[0]=0;
#if 1
#if A2_NCTX==4
					long long p0=0;
					int wsum=0;
					ALIGN(32) long long probs[A2_NCTX];
					ALIGN(32) long long apsum[4], awsum[4];
					__m256i mprob0=_mm256_set_epi64x(curr_stats[ 3][tidx], curr_stats[ 2][tidx], curr_stats[ 1][tidx], curr_stats[ 0][tidx]);
					__m256i mweight0=_mm256_load_si256((__m256i*)curr_mixer+0);
					__m256i mprod0=_mm256_mul_epi32(mprob0, mweight0);
					_mm256_store_si256((__m256i*)awsum, mweight0);
					_mm256_store_si256((__m256i*)apsum, mprod0);
					wsum=(int)(awsum[0]+awsum[1]+awsum[2]+awsum[3]);
					p0=apsum[0]+apsum[1]+apsum[2]+apsum[3];
					wsum+=!wsum;
					p0+=wsum>>1;
					p0/=wsum;
#elif A2_NCTX==8
					long long p0=0;
					int wsum=0;
					ALIGN(32) long long probs[A2_NCTX];
					ALIGN(32) long long apsum[4], awsum[4];
					__m256i mprob0=_mm256_set_epi64x(curr_stats[ 3][tidx], curr_stats[ 2][tidx], curr_stats[ 1][tidx], curr_stats[ 0][tidx]);
					__m256i mprob1=_mm256_set_epi64x(curr_stats[ 7][tidx], curr_stats[ 6][tidx], curr_stats[ 5][tidx], curr_stats[ 4][tidx]);
					__m256i mweight0=_mm256_load_si256((__m256i*)curr_mixer+0);
					__m256i mweight1=_mm256_load_si256((__m256i*)curr_mixer+1);
					__m256i mprod0=_mm256_mul_epi32(mprob0, mweight0);
					__m256i mprod1=_mm256_mul_epi32(mprob1, mweight1);
					__m256i mwsum=_mm256_add_epi64(mweight0, mweight1);
					__m256i mpsum=_mm256_add_epi64(mprod0, mprod1);
					_mm256_store_si256((__m256i*)awsum, mwsum);
					_mm256_store_si256((__m256i*)apsum, mpsum);
					wsum=(int)(awsum[0]+awsum[1]+awsum[2]+awsum[3]);
					p0=apsum[0]+apsum[1]+apsum[2]+apsum[3];
					wsum+=!wsum;
					p0+=wsum>>1;
					p0/=wsum;
#elif A2_NCTX==12
					long long p0=0;
					int wsum=0;
					ALIGN(32) long long probs[A2_NCTX];
					ALIGN(32) long long apsum[4], awsum[4];
#ifdef ENABLE_GATHER
					//if(ky==BLOCKY/2&&kx==BLOCKX/2)
					//	printf("");
					__m256i mprob0=_mm256_add_epi64(midx0, mtidx);
					__m256i mprob1=_mm256_add_epi64(midx1, mtidx);
					__m256i mprob2=_mm256_add_epi64(midx2, mtidx);
					mprob0=_mm256_i64gather_epi64((long long*)statsptr, mprob0, 2);
					mprob1=_mm256_i64gather_epi64((long long*)statsptr, mprob1, 2);
					mprob2=_mm256_i64gather_epi64((long long*)statsptr, mprob2, 2);
					mprob0=_mm256_and_si256(mprob0, probmask);
					mprob1=_mm256_and_si256(mprob1, probmask);
					mprob2=_mm256_and_si256(mprob2, probmask);
#else
					__m256i mprob0=_mm256_set_epi64x(curr_stats[ 3][tidx], curr_stats[ 2][tidx], curr_stats[ 1][tidx], curr_stats[ 0][tidx]);
					__m256i mprob1=_mm256_set_epi64x(curr_stats[ 7][tidx], curr_stats[ 6][tidx], curr_stats[ 5][tidx], curr_stats[ 4][tidx]);
					__m256i mprob2=_mm256_set_epi64x(curr_stats[11][tidx], curr_stats[10][tidx], curr_stats[ 9][tidx], curr_stats[ 8][tidx]);
#endif
					__m256i mweight0=_mm256_load_si256((__m256i*)curr_mixer+0);
					__m256i mweight1=_mm256_load_si256((__m256i*)curr_mixer+1);
					__m256i mweight2=_mm256_load_si256((__m256i*)curr_mixer+2);
					//mweight0=_mm256_slli_epi32(mweight0, 1);
					//mweight1=_mm256_slli_epi32(mweight1, 1);
					//mweight2=_mm256_srai_epi32(mweight2, 1);
					__m256i mprod0=_mm256_mul_epi32(mprob0, mweight0);
					__m256i mprod1=_mm256_mul_epi32(mprob1, mweight1);
					__m256i mprod2=_mm256_mul_epi32(mprob2, mweight2);
					__m256i mwsum=_mm256_add_epi64(mweight0, mweight1);
					mwsum=_mm256_add_epi64(mwsum, mweight2);
					__m256i mpsum=_mm256_add_epi64(mprod0, mprod1);
					mpsum=_mm256_add_epi64(mpsum, mprod2);
					_mm256_store_si256((__m256i*)awsum, mwsum);
					_mm256_store_si256((__m256i*)apsum, mpsum);
					wsum=(int)(awsum[0]+awsum[1]+awsum[2]+awsum[3]);
					p0=apsum[0]+apsum[1]+apsum[2]+apsum[3];
					wsum+=!wsum;
					p0+=wsum>>1;
					p0/=wsum;
#elif A2_NCTX==16
					//if(kx==100&&ky==100)//
					//	printf("");

					short probs[16]=
					{
						curr_stats[ 0][tidx], curr_stats[ 1][tidx], curr_stats[ 2][tidx], curr_stats[ 3][tidx],
						curr_stats[ 4][tidx], curr_stats[ 5][tidx], curr_stats[ 6][tidx], curr_stats[ 7][tidx],
						curr_stats[ 8][tidx], curr_stats[ 9][tidx], curr_stats[10][tidx], curr_stats[11][tidx],
						curr_stats[12][tidx], curr_stats[13][tidx], curr_stats[14][tidx], curr_stats[15][tidx],
					};
					int prob32[4]=
					{
						(curr_mixer1[ 0]*probs[ 0]+curr_mixer1[ 1]*probs[ 1]+curr_mixer1[ 2]*probs[ 2]+curr_mixer1[ 3]*probs[ 3])>>12,//(+-128)*(+128)*4>>4 = 0x1000
						(curr_mixer1[ 4]*probs[ 4]+curr_mixer1[ 5]*probs[ 5]+curr_mixer1[ 6]*probs[ 6]+curr_mixer1[ 7]*probs[ 7])>>12,
						(curr_mixer1[ 8]*probs[ 8]+curr_mixer1[ 9]*probs[ 9]+curr_mixer1[10]*probs[10]+curr_mixer1[11]*probs[11])>>12,
						(curr_mixer1[12]*probs[12]+curr_mixer1[13]*probs[13]+curr_mixer1[14]*probs[14]+curr_mixer1[15]*probs[15])>>12,
					};
					//CLAMP2(prob32[0], -0x4000, 0x4000);
					//CLAMP2(prob32[1], -0x4000, 0x4000);
					//CLAMP2(prob32[2], -0x4000, 0x4000);
					//CLAMP2(prob32[3], -0x4000, 0x4000);
					int p0=(int)((
						+(long long)curr_mixer2[0]*prob32[0]
						+(long long)curr_mixer2[1]*prob32[1]
						+(long long)curr_mixer2[2]*prob32[2]
						+(long long)curr_mixer2[3]*prob32[3]
						+(1LL<<10>>1)
					)>>10);//0x10000*0x1000*4>>14 = 0x10000

					//p0=
					//(
					//	+curr_stats[ 0][tidx]+curr_stats[ 1][tidx]+curr_stats[ 2][tidx]+curr_stats[ 3][tidx]
					//	+curr_stats[ 4][tidx]+curr_stats[ 5][tidx]+curr_stats[ 6][tidx]+curr_stats[ 7][tidx]
					//	+curr_stats[ 8][tidx]+curr_stats[ 9][tidx]+curr_stats[10][tidx]+curr_stats[11][tidx]
					//	+curr_stats[12][tidx]+curr_stats[13][tidx]+curr_stats[14][tidx]+curr_stats[15][tidx]
					//)<<4;
					p0+=0x8000;

					//__m128i mmix1=_mm_load_si128((__m128i*)curr_mixer1);
					//__m128i mmix2=_mm_load_si128((__m128i*)curr_mixer2);
					//__m128i mprob8=_mm_set_epi8(
					//	curr_stats[15][tidx], curr_stats[14][tidx], curr_stats[13][tidx], curr_stats[12][tidx],
					//	curr_stats[11][tidx], curr_stats[10][tidx], curr_stats[ 9][tidx], curr_stats[ 8][tidx],
					//	curr_stats[ 7][tidx], curr_stats[ 6][tidx], curr_stats[ 5][tidx], curr_stats[ 4][tidx],
					//	curr_stats[ 3][tidx], curr_stats[ 2][tidx], curr_stats[ 1][tidx], curr_stats[ 0][tidx]
					//);
					//__m128i mprob32=_mm_dpbusd_epi32(mmix1, mprob8, _mm_setzero_si128());//(+-128)*(+128) * 4 = 0x10000
					//mprob32=_mm_srai_epi32(mprob32, 1);
					//mprob32=_mm_min_epi32(mprob32, _mm_set1_epi32(+0x3FFF));
					//mprob32=_mm_max_epi32(mprob32, _mm_set1_epi32(-0x3FFF));
					//mprob32=_mm_and_si128(mprob32, _mm_set1_epi32(0xFFFF));
					//__m128i mp=_mm_madd_epi16(mprob32, mmix2);//0x4000*0x10000 = 0x40000000		2x faster than mullo32
					//mp=_mm_srai_epi32(mp, 14);//0x40000000>>14 = 0x10000
					//ALIGN(16) int probs[4];
					//_mm_store_si128((__m128i*)probs, mp);
					//int p0=(probs[0]+probs[1]+probs[2]+probs[3])>>2;
					//p0+=0x8000;
#endif

					//long long p0_2=0;
					//int wsum2=0;
					//for(int k=0;k<A2_NCTX;++k)
					//{
					//	p0_2+=(long long)curr_mixer[k]*curr_stats[k][tidx];
					//	wsum2+=curr_mixer[k];
					//}
					//if(wsum2)
					//	p0_2/=wsum2;
					//else
					//	p0_2=0x8000, wsum2=1;
					//if(p0!=p0_2)
					//	LOG_ERROR("");
#else
					ALIGN(32) long long probs[A2_NCTX];
					for(int k=0;k<A2_NCTX;++k)
					{
						probs[k]=curr_stats[k][tidx];
						p0+=(long long)curr_mixer[k]*probs[k];
						wsum+=curr_mixer[k];
					}
					if(wsum)
						p0/=wsum;
					else
						p0=0x8000, wsum=1;
#endif
#ifdef ENABLE_SSE3_PROB
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][((error&0x1FE<<kb)+pred)>>3&31][p0>>11&31];	//2.859087
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][(pred+pred2)>>5&15][p0>>11&31];	//2.859014 not worth it
					long long *curr_sse_prob=&args->sse_prob[kc][kb][(((error&0x1FE<<kb)|1<<kb)+pred)>>3&31][p0>>11&31];	//2.859006 CLIC303 record
					long long sseprobval=*curr_sse_prob;

					//long long *curr_sse_prob=&args->sse_prob[kc][kb][pred>>3&31][p0>>11&31];	//2.858704
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][pred>>2&63][p0>>11&31];	//2.858688
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][pred>>4&15][p0>>11&31];	//2.858656
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][(pred+pred2+16)>>5&15][p0>>11&31];	//2.858611
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][pred>>2&63][p0>>12&15];	//2.858482
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][pred>>4&15][p0>>10&63];	//2.858152
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][(partial+pred)>>5&15][p0>>11&31];	//2.858111
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][pred>>3&31][p0>>10&63];	//2.858067
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][pred>>2&63][p0>>10&63];	//2.857870
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][pred>>3&31][p0>>9&127];	//2.857620
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][pred>>2&63][p0>>9&127];	//2.857519
					//long long *curr_sse_prob=&args->sse_prob[kc][kb][(error&0x1FE<<kb)>>3&31][p0>>11&31];	//2.855330
					//int change=(*curr_sse_prob+(1<<(SSE3_DECAY+16)>>1))>>(SSE3_DECAY+16);
					//if(abs(change)>100)
					//	printf("");
					p0+=(sseprobval+(1LL<<SSE3_SCALE>>1))>>SSE3_SCALE;
					//p0+=(int)((*curr_sse_prob+(1LL<<(SSE3_DECAY+SSE3_SCALE)>>1))>>(SSE3_DECAY+SSE3_SCALE));
#endif
					CLAMP2(p0, 0x0001, 0xFFFF);
					if(args->fwd)
					{
						bit=error>>kb&1;
						ac3_enc_bin(&ec, bit, (int)p0, 16);
#ifdef ENABLE_ZIPF_VIEW
						psize-=log2((double)(bit?0x10000-p0:p0)/0x10000);
#endif
					}
					else
					{
						bit=ac3_dec_bin(&ec, (int)p0, 16);
						error|=bit<<kb;
					}
					
					int truth=!bit<<16;
					int prob_error=truth-(int)p0;
#if A2_NCTX==16
					long long dL_dp0=0x171547652/((bit<<16)-p0);
					CLAMP2(dL_dp0, -0x10000, 0x10000);
					int gain=(int)dL_dp0;
					int mixupdate2[]=
					{
						curr_mixer2[0]-(int)(((long long)gain*prob32[0]+(1<<9>>1))>>9),//0x4000*0x1000>>18 = 0x100
						curr_mixer2[1]-(int)(((long long)gain*prob32[1]+(1<<9>>1))>>9),
						curr_mixer2[2]-(int)(((long long)gain*prob32[2]+(1<<9>>1))>>9),
						curr_mixer2[3]-(int)(((long long)gain*prob32[3]+(1<<9>>1))>>9),
					};
					long long mixupdate1_temp_ll[]=
					{
						((long long)gain*curr_mixer2[0]+(1<<14>>1))>>14,
						((long long)gain*curr_mixer2[1]+(1<<14>>1))>>14,
						((long long)gain*curr_mixer2[2]+(1<<14>>1))>>14,
						((long long)gain*curr_mixer2[3]+(1<<14>>1))>>14,
					};
					CLAMP2(mixupdate1_temp_ll[0], -0x80000, 0x80000);
					CLAMP2(mixupdate1_temp_ll[1], -0x80000, 0x80000);
					CLAMP2(mixupdate1_temp_ll[2], -0x80000, 0x80000);
					CLAMP2(mixupdate1_temp_ll[3], -0x80000, 0x80000);
					int mixupdate1_temp[]=
					{
						(int)mixupdate1_temp_ll[0],
						(int)mixupdate1_temp_ll[1],
						(int)mixupdate1_temp_ll[2],
						(int)mixupdate1_temp_ll[3],
					};
					int mixupdate1[]=
					{
						curr_mixer1[ 0]-(int)(((long long)mixupdate1_temp[0]*probs[ 0]+(1<<9>>1))>>9),//(0x4000*0x4000>>8)*0x80>>24 = 0x8
						curr_mixer1[ 1]-(int)(((long long)mixupdate1_temp[0]*probs[ 1]+(1<<9>>1))>>9),
						curr_mixer1[ 2]-(int)(((long long)mixupdate1_temp[0]*probs[ 2]+(1<<9>>1))>>9),
						curr_mixer1[ 3]-(int)(((long long)mixupdate1_temp[0]*probs[ 3]+(1<<9>>1))>>9),
						curr_mixer1[ 4]-(int)(((long long)mixupdate1_temp[1]*probs[ 4]+(1<<9>>1))>>9),
						curr_mixer1[ 5]-(int)(((long long)mixupdate1_temp[1]*probs[ 5]+(1<<9>>1))>>9),
						curr_mixer1[ 6]-(int)(((long long)mixupdate1_temp[1]*probs[ 6]+(1<<9>>1))>>9),
						curr_mixer1[ 7]-(int)(((long long)mixupdate1_temp[1]*probs[ 7]+(1<<9>>1))>>9),
						curr_mixer1[ 8]-(int)(((long long)mixupdate1_temp[2]*probs[ 8]+(1<<9>>1))>>9),
						curr_mixer1[ 9]-(int)(((long long)mixupdate1_temp[2]*probs[ 9]+(1<<9>>1))>>9),
						curr_mixer1[10]-(int)(((long long)mixupdate1_temp[2]*probs[10]+(1<<9>>1))>>9),
						curr_mixer1[11]-(int)(((long long)mixupdate1_temp[2]*probs[11]+(1<<9>>1))>>9),
						curr_mixer1[12]-(int)(((long long)mixupdate1_temp[3]*probs[12]+(1<<9>>1))>>9),
						curr_mixer1[13]-(int)(((long long)mixupdate1_temp[3]*probs[13]+(1<<9>>1))>>9),
						curr_mixer1[14]-(int)(((long long)mixupdate1_temp[3]*probs[14]+(1<<9>>1))>>9),
						curr_mixer1[15]-(int)(((long long)mixupdate1_temp[3]*probs[15]+(1<<9>>1))>>9),
					};
					int probupdate[]=
					{
						probs[ 0]-(int)(((long long)mixupdate1_temp[0]*curr_mixer1[ 0]+(1<<7>>1))>>7),//(0x4000*0x4000>>8)*0x80>>20 = 0x8
						probs[ 1]-(int)(((long long)mixupdate1_temp[0]*curr_mixer1[ 1]+(1<<7>>1))>>7),
						probs[ 2]-(int)(((long long)mixupdate1_temp[0]*curr_mixer1[ 2]+(1<<7>>1))>>7),
						probs[ 3]-(int)(((long long)mixupdate1_temp[0]*curr_mixer1[ 3]+(1<<7>>1))>>7),
						probs[ 4]-(int)(((long long)mixupdate1_temp[1]*curr_mixer1[ 4]+(1<<7>>1))>>7),
						probs[ 5]-(int)(((long long)mixupdate1_temp[1]*curr_mixer1[ 5]+(1<<7>>1))>>7),
						probs[ 6]-(int)(((long long)mixupdate1_temp[1]*curr_mixer1[ 6]+(1<<7>>1))>>7),
						probs[ 7]-(int)(((long long)mixupdate1_temp[1]*curr_mixer1[ 7]+(1<<7>>1))>>7),
						probs[ 8]-(int)(((long long)mixupdate1_temp[2]*curr_mixer1[ 8]+(1<<7>>1))>>7),
						probs[ 9]-(int)(((long long)mixupdate1_temp[2]*curr_mixer1[ 9]+(1<<7>>1))>>7),
						probs[10]-(int)(((long long)mixupdate1_temp[2]*curr_mixer1[10]+(1<<7>>1))>>7),
						probs[11]-(int)(((long long)mixupdate1_temp[2]*curr_mixer1[11]+(1<<7>>1))>>7),
						probs[12]-(int)(((long long)mixupdate1_temp[3]*curr_mixer1[12]+(1<<7>>1))>>7),
						probs[13]-(int)(((long long)mixupdate1_temp[3]*curr_mixer1[13]+(1<<7>>1))>>7),
						probs[14]-(int)(((long long)mixupdate1_temp[3]*curr_mixer1[14]+(1<<7>>1))>>7),
						probs[15]-(int)(((long long)mixupdate1_temp[3]*curr_mixer1[15]+(1<<7>>1))>>7),
					};
					CLAMP2(mixupdate2[0], 1, 0xFFFF);
					CLAMP2(mixupdate2[1], 1, 0xFFFF);
					CLAMP2(mixupdate2[2], 1, 0xFFFF);
					CLAMP2(mixupdate2[3], 1, 0xFFFF);
					CLAMP2(mixupdate1[ 0], 1, 0xFFFF);
					CLAMP2(mixupdate1[ 1], 1, 0xFFFF);
					CLAMP2(mixupdate1[ 2], 1, 0xFFFF);
					CLAMP2(mixupdate1[ 3], 1, 0xFFFF);
					CLAMP2(mixupdate1[ 4], 1, 0xFFFF);
					CLAMP2(mixupdate1[ 5], 1, 0xFFFF);
					CLAMP2(mixupdate1[ 6], 1, 0xFFFF);
					CLAMP2(mixupdate1[ 7], 1, 0xFFFF);
					CLAMP2(mixupdate1[ 8], 1, 0xFFFF);
					CLAMP2(mixupdate1[ 9], 1, 0xFFFF);
					CLAMP2(mixupdate1[10], 1, 0xFFFF);
					CLAMP2(mixupdate1[11], 1, 0xFFFF);
					CLAMP2(mixupdate1[12], 1, 0xFFFF);
					CLAMP2(mixupdate1[13], 1, 0xFFFF);
					CLAMP2(mixupdate1[14], 1, 0xFFFF);
					CLAMP2(mixupdate1[15], 1, 0xFFFF);
					CLAMP2(probupdate[ 0], 1, 0xFFFF);
					CLAMP2(probupdate[ 1], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[ 2], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[ 3], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[ 4], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[ 5], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[ 6], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[ 7], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[ 8], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[ 9], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[10], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[11], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[12], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[13], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[14], -0x7FFF, 0x7FFF);
					CLAMP2(probupdate[15], -0x7FFF, 0x7FFF);
					curr_mixer2[0]=mixupdate2[0];
					curr_mixer2[1]=mixupdate2[1];
					curr_mixer2[2]=mixupdate2[2];
					curr_mixer2[3]=mixupdate2[3];
					curr_mixer1[ 0]=(unsigned char)mixupdate1[ 0];
					curr_mixer1[ 1]=(unsigned char)mixupdate1[ 1];
					curr_mixer1[ 2]=(unsigned char)mixupdate1[ 2];
					curr_mixer1[ 3]=(unsigned char)mixupdate1[ 3];
					curr_mixer1[ 4]=(unsigned char)mixupdate1[ 4];
					curr_mixer1[ 5]=(unsigned char)mixupdate1[ 5];
					curr_mixer1[ 6]=(unsigned char)mixupdate1[ 6];
					curr_mixer1[ 7]=(unsigned char)mixupdate1[ 7];
					curr_mixer1[ 8]=(unsigned char)mixupdate1[ 8];
					curr_mixer1[ 9]=(unsigned char)mixupdate1[ 9];
					curr_mixer1[10]=(unsigned char)mixupdate1[10];
					curr_mixer1[11]=(unsigned char)mixupdate1[11];
					curr_mixer1[12]=(unsigned char)mixupdate1[12];
					curr_mixer1[13]=(unsigned char)mixupdate1[13];
					curr_mixer1[14]=(unsigned char)mixupdate1[14];
					curr_mixer1[15]=(unsigned char)mixupdate1[15];
					curr_stats[ 0][tidx]=(unsigned char)probupdate[ 0];
					curr_stats[ 1][tidx]=(unsigned char)probupdate[ 1];
					curr_stats[ 2][tidx]=(unsigned char)probupdate[ 2];
					curr_stats[ 3][tidx]=(unsigned char)probupdate[ 3];
					curr_stats[ 4][tidx]=(unsigned char)probupdate[ 4];
					curr_stats[ 5][tidx]=(unsigned char)probupdate[ 5];
					curr_stats[ 6][tidx]=(unsigned char)probupdate[ 6];
					curr_stats[ 7][tidx]=(unsigned char)probupdate[ 7];
					curr_stats[ 8][tidx]=(unsigned char)probupdate[ 8];
					curr_stats[ 9][tidx]=(unsigned char)probupdate[ 9];
					curr_stats[10][tidx]=(unsigned char)probupdate[10];
					curr_stats[11][tidx]=(unsigned char)probupdate[11];
					curr_stats[12][tidx]=(unsigned char)probupdate[12];
					curr_stats[13][tidx]=(unsigned char)probupdate[13];
					curr_stats[14][tidx]=(unsigned char)probupdate[14];
					curr_stats[15][tidx]=(unsigned char)probupdate[15];
					curr_mixer1+=A2_NCTX;
					curr_mixer2+=A2_NCTX/4;
#else
					if(abs(prob_error)>512)
					{
#if 1
#if A2_NCTX==4
						long long dL_dp0=0xE0000000000/((long long)((bit<<16)-(int)p0)*wsum);
						CLAMP2(dL_dp0, -0x3FFF, 0x3FFF);
						__m256i vmin2=_mm256_set1_epi64x(1);
						__m256i mp0=_mm256_set1_epi64x(p0);
						__m256i mdL_dp0=_mm256_set1_epi64x(dL_dp0);
						__m256i update0=_mm256_sub_epi32(mprob0, mp0);
						update0=_mm256_mul_epi32(update0, mdL_dp0);
						update0=_mm256_srav_epi32(update0, _mm256_load_si256((__m256i*)curr_mixsh+0));
						mweight0=_mm256_sub_epi32(mweight0, update0);
						mweight0=_mm256_max_epi32(mweight0, vmin2);
						_mm256_store_si256((__m256i*)curr_mixer+0, mweight0);
#elif A2_NCTX==8
						long long dL_dp0=0xE0000000000/((long long)((bit<<16)-(int)p0)*wsum);
						CLAMP2(dL_dp0, -0x3FFF, 0x3FFF);
						__m256i vmin2=_mm256_set1_epi64x(1);
						__m256i mp0=_mm256_set1_epi64x(p0);
						__m256i mdL_dp0=_mm256_set1_epi64x(dL_dp0);
						__m256i update0=_mm256_sub_epi32(mprob0, mp0);
						__m256i update1=_mm256_sub_epi32(mprob1, mp0);
						update0=_mm256_mul_epi32(update0, mdL_dp0);
						update1=_mm256_mul_epi32(update1, mdL_dp0);
						update0=_mm256_srav_epi32(update0, _mm256_load_si256((__m256i*)curr_mixsh+0));
						update1=_mm256_srav_epi32(update1, _mm256_load_si256((__m256i*)curr_mixsh+1));
						mweight0=_mm256_sub_epi32(mweight0, update0);
						mweight1=_mm256_sub_epi32(mweight1, update1);
						mweight0=_mm256_max_epi32(mweight0, vmin2);
						mweight1=_mm256_max_epi32(mweight1, vmin2);
						_mm256_store_si256((__m256i*)curr_mixer+0, mweight0);
						_mm256_store_si256((__m256i*)curr_mixer+1, mweight1);
#elif A2_NCTX==12
						long long dL_dp0=0xE0000000000/((long long)((bit<<16)-(int)p0)*wsum);
						CLAMP2(dL_dp0, -0x3FFF, 0x3FFF);
						//mweight0=_mm256_load_si256((__m256i*)curr_mixer+0);
						//mweight1=_mm256_load_si256((__m256i*)curr_mixer+1);
						//mweight2=_mm256_load_si256((__m256i*)curr_mixer+2);
						__m256i vmin2=_mm256_set1_epi64x(1);
						__m256i mp0=_mm256_set1_epi64x(p0);
						__m256i mdL_dp0=_mm256_set1_epi64x(dL_dp0);
						__m256i update0=_mm256_sub_epi32(mprob0, mp0);
						__m256i update1=_mm256_sub_epi32(mprob1, mp0);
						__m256i update2=_mm256_sub_epi32(mprob2, mp0);
						update0=_mm256_mul_epi32(update0, mdL_dp0);
						update1=_mm256_mul_epi32(update1, mdL_dp0);
						update2=_mm256_mul_epi32(update2, mdL_dp0);
						update0=_mm256_srav_epi32(update0, _mm256_load_si256((__m256i*)curr_mixsh+0));
						update1=_mm256_srav_epi32(update1, _mm256_load_si256((__m256i*)curr_mixsh+1));
						update2=_mm256_srav_epi32(update2, _mm256_load_si256((__m256i*)curr_mixsh+2));
						mweight0=_mm256_sub_epi32(mweight0, update0);
						mweight1=_mm256_sub_epi32(mweight1, update1);
						mweight2=_mm256_sub_epi32(mweight2, update2);
						mweight0=_mm256_max_epi32(mweight0, vmin2);
						mweight1=_mm256_max_epi32(mweight1, vmin2);
						mweight2=_mm256_max_epi32(mweight2, vmin2);
					//	mweight0=_mm256_and_si256(mweight0, lo32);
					//	mweight1=_mm256_and_si256(mweight1, lo32);
					//	mweight2=_mm256_and_si256(mweight2, lo32);
						_mm256_store_si256((__m256i*)curr_mixer+0, mweight0);
						_mm256_store_si256((__m256i*)curr_mixer+1, mweight1);
						_mm256_store_si256((__m256i*)curr_mixer+2, mweight2);
#endif
#else
						long long dL_dp0=0xE0000000000/((long long)((bit<<16)-(int)p0)*wsum);
						for(int k=0;k<A2_NCTX;++k)
						{
							int update=(int)((int)dL_dp0*(curr_stats[k][tidx]-(int)p0)>>(1+17-4-1));
							int mk=curr_mixer[k]-update;
							CLAMP2(mk, 1, 0x40000);
							curr_mixer[k]=mk;
						}
#endif
					}
#if 1
#if A2_NCTX==4
					__m256i mtruth=_mm256_set1_epi64x(truth);
					__m256i mupdate0=_mm256_sub_epi64(_mm256_or_si256(mtruth, _mm256_load_si256((__m256i*)g_offset+0)), mprob0);
					mupdate0=_mm256_srav_epi32(mupdate0, _mm256_load_si256((__m256i*)curr_sh+0));
					mprob0=_mm256_add_epi64(mprob0, mupdate0);
					_mm256_store_si256((__m256i*)probs+0, mprob0);
#elif A2_NCTX==8
					__m256i mtruth=_mm256_set1_epi64x(truth);
					__m256i mupdate0=_mm256_sub_epi64(_mm256_or_si256(mtruth, _mm256_load_si256((__m256i*)g_offset+0)), mprob0);
					__m256i mupdate1=_mm256_sub_epi64(_mm256_or_si256(mtruth, _mm256_load_si256((__m256i*)g_offset+1)), mprob1);
					mupdate0=_mm256_srav_epi32(mupdate0, _mm256_load_si256((__m256i*)curr_sh+0));
					mupdate1=_mm256_srav_epi32(mupdate1, _mm256_load_si256((__m256i*)curr_sh+1));
					mprob0=_mm256_add_epi64(mprob0, mupdate0);
					mprob1=_mm256_add_epi64(mprob1, mupdate1);
					_mm256_store_si256((__m256i*)probs+0, mprob0);
					_mm256_store_si256((__m256i*)probs+1, mprob1);
#elif A2_NCTX==12
					__m256i mtruth=_mm256_set1_epi64x(truth);
					__m256i mupdate0=_mm256_sub_epi64(_mm256_or_si256(mtruth, _mm256_load_si256((__m256i*)g_offset+0)), mprob0);
					__m256i mupdate1=_mm256_sub_epi64(_mm256_or_si256(mtruth, _mm256_load_si256((__m256i*)g_offset+1)), mprob1);
					__m256i mupdate2=_mm256_sub_epi64(_mm256_or_si256(mtruth, _mm256_load_si256((__m256i*)g_offset+2)), mprob2);
					mupdate0=_mm256_srav_epi32(mupdate0, _mm256_load_si256((__m256i*)curr_sh+0));
					mupdate1=_mm256_srav_epi32(mupdate1, _mm256_load_si256((__m256i*)curr_sh+1));
					mupdate2=_mm256_srav_epi32(mupdate2, _mm256_load_si256((__m256i*)curr_sh+2));
					mprob0=_mm256_add_epi64(mprob0, mupdate0);
					mprob1=_mm256_add_epi64(mprob1, mupdate1);
					mprob2=_mm256_add_epi64(mprob2, mupdate2);
					_mm256_store_si256((__m256i*)probs+0, mprob0);
					_mm256_store_si256((__m256i*)probs+1, mprob1);
					_mm256_store_si256((__m256i*)probs+2, mprob2);
#endif
#ifdef __GNUC__
#pragma GCC unroll 12
#endif
					for(int k=0;k<A2_NCTX;++k)//scatter
						curr_stats[k][tidx]=(int)probs[k];
#else
					int truth=!bit<<16;
#ifdef __GNUC__
#pragma GCC unroll 12
#endif
					for(int k=0;k<A2_NCTX;++k)
					{
						int p=curr_stats[k][tidx];
						p+=(truth-p+g_offset[k])>>curr_sh[k];
					//	CLAMP2(p, 1, 0xFFFF);
						curr_stats[k][tidx]=p;
					}
#endif
					curr_mixer+=A2_NCTX;
#endif
#ifdef ENABLE_SSE3_PROB
					long long e=(long long)prob_error<<SSE3_SCALE;

					*curr_sse_prob=sseprobval+((e-sseprobval+(1LL<<SSE3_DECAY>>1))>>SSE3_DECAY);
					//*curr_sse_prob=((sseprobval*((1LL<<SSE3_DECAY)-1)+(1LL<<SSE3_DECAY>>1))>>SSE3_DECAY)+e;
					//*curr_sse_prob=((*curr_sse_prob*((1LL<<SSE3_DECAY)-1)+(1LL<<SSE3_DECAY>>1))>>SSE3_DECAY)+e;
					//if(!kb)
					//	printf("%lld\n", *curr_sse_prob);
					//if(abs(*curr_sse_prob>>2)>0xFFFF)
					//	printf("");
#endif
					tidx+=tidx+bit;
#ifdef ENABLE_GATHER
					mtidx=_mm256_slli_epi32(mtidx, 1);
					if(bit)
						mtidx=_mm256_or_si256(mtidx, _mm256_set1_epi64x(1));
#endif
				}
				if(!args->fwd)
				{
					error=(char)(error+pred);
					yuv[kc]=error;
					curr[kc2+0]=yuv[kc];
					curr[kc2+1]=curr[kc2+0]-pred;
				}
#endif
#ifdef ESTIMATE_SIZE
				++args->hist2[kc][(curr[kc2+1]+128)&255];
#endif
#ifdef ENABLE_SSE2
				//if(curr[kc2+0]-pred0!=curr[kc2+1])//
				//	printf("");
				int e=curr[kc2+1]<<SSE2_SCALE;
				*curr_sse1=sseval+((e-sseval+(1<<SSE2_DECAY>>1))>>SSE2_DECAY);
			//	*curr_sse1=((*curr_sse1*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
			//	*curr_sse2=((*curr_sse2*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
			//	*curr_sse3=((*curr_sse3*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
			//	*curr_sse4=((*curr_sse4*((1<<SSE2_DECAY)-1)+(1<<SSE2_DECAY>>1))>>SSE2_DECAY)+e;
#endif
				curr[kc2+0]-=offset;
#ifdef ENABLE_WG
				wg_update(curr[kc2+0], kc, wg_preds, wg_perrors+WG_NPREDS*kc, eN, eW, ecurr, eNE);
#endif
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
#ifdef ENABLE_ZIPF_VIEW
			else
			{
				psize*=255./24;
				CLAMP2(psize, 0, 255);
				args->zimage[args->iw*ky+kx]=(unsigned char)psize;
			}
#endif
			__m256i mr=_mm256_load_si256((__m256i*)rows);
			__m256i er=_mm256_load_si256((__m256i*)erows);
			mr=_mm256_add_epi64(mr, _mm256_set1_epi64x(sizeof(short[4*2])));
			er=_mm256_add_epi64(er, _mm256_set1_epi64x(sizeof(int[4*WG_NPREDS])));
			_mm256_store_si256((__m256i*)rows, mr);
			_mm256_store_si256((__m256i*)erows, er);
			//rows[0]+=4*2;
			//rows[1]+=4*2;
			//rows[2]+=4*2;
			//rows[3]+=4*2;
		}
		if(ky-args->y1==BLOCKY/7)
		{
			__m256i one=_mm256_set1_epi32(1);
			__m256i vmax=_mm256_set1_epi32(8);
			__m256i sh0=_mm256_load_si256((__m256i*)curr_sh+0);
			__m256i sh1=_mm256_load_si256((__m256i*)curr_sh+1);
		//	__m256i sh2=_mm256_load_si256((__m256i*)curr_sh+2);
			sh0=_mm256_add_epi32(sh0, one);
			sh1=_mm256_add_epi32(sh1, one);
		//	sh2=_mm256_add_epi32(sh2, one);
			sh0=_mm256_min_epi32(sh0, vmax);
			sh1=_mm256_min_epi32(sh1, vmax);
		//	sh2=_mm256_min_epi32(sh2, vmax);
			_mm256_store_si256((__m256i*)curr_sh+0, sh0);
			_mm256_store_si256((__m256i*)curr_sh+1, sh1);
		//	_mm256_store_si256((__m256i*)curr_sh+2, sh2);
			curr_sh[0]=g_sh[0];
			curr_sh[3]=g_sh[3];

			sh0=_mm256_load_si256((__m256i*)curr_mixsh+0);
			sh1=_mm256_load_si256((__m256i*)curr_mixsh+1);
		//	sh2=_mm256_load_si256((__m256i*)curr_mixsh+2);
			sh0=_mm256_sub_epi32(sh0, one);
			sh1=_mm256_sub_epi32(sh1, one);
		//	sh2=_mm256_sub_epi32(sh2, one);
		//	sh0=_mm256_max_epi32(sh0, vmax);
		//	sh1=_mm256_max_epi32(sh1, vmax);
		//	sh2=_mm256_max_epi32(sh2, vmax);
			_mm256_store_si256((__m256i*)curr_mixsh+0, sh0);
			_mm256_store_si256((__m256i*)curr_mixsh+1, sh1);
		//	_mm256_store_si256((__m256i*)curr_mixsh+2, sh2);
		}
	}
	if(args->fwd)
		ac3_enc_flush(&ec);
}
static void block_manager(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	for(int kb=args->b1;kb<args->b2;++kb)
	{
		int kx, ky;

		args->blockidx=kb;
		args->currentblock=kb-args->b1;
		ky=args->blockidx/args->xblocks;
		kx=args->blockidx%args->xblocks;
		args->x1=BLOCKX*kx;
		args->y1=BLOCKY*ky;
		args->x2=MINVAR(args->x1+BLOCKX, args->iw);
		args->y2=MINVAR(args->y1+BLOCKY, args->ih);
		if(!args->fwd)
		{
			args->decstart=args->decsrc+args->offsets[kb];
			args->decend=args->decsrc+args->offsets[kb+1];
		}
		block_thread(param);
	}
}
int c15_codec(int argc, char **argv)
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
	int maxthreads=argc<4?0:atoi(argv[3]);
	const int nch=3;
	double t0;
	ArrayHandle src, dst;
	int headersize, printed;
	int iw, ih;
	const unsigned char *image, *imageend;
	unsigned char *image2;
	CodecID codec;
	int xblocks, yblocks, nblocks, ncores, nthreads, blocksperthread, coffset;
	ptrdiff_t start, memusage, argssize;
	ThreadArgs *args;
	int test, fwd;
	int histsize;
	double bestsize;
	int usize;
#ifdef ESTIMATE_SIZE
	double esize[3]={0};
#endif
#ifdef ENABLE_ZIPF_VIEW
	unsigned char *zimage=0;
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
	{
		printf("%s \"%s\"  WH %d*%d\n", CODECNAME, srcfn, iw, ih);
#ifdef ENABLE_ZIPF_VIEW
		zimage=(unsigned char*)malloc((size_t)iw*ih);
		if(!zimage)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(zimage, 0, (size_t)iw*ih);
#endif
	}
	usize=iw*ih*3;
	xblocks=(iw+BLOCKX-1)/BLOCKX;
	yblocks=(ih+BLOCKY-1)/BLOCKY;
	nblocks=xblocks*yblocks;
	ncores=query_cpu_cores();
	nthreads=MINVAR(nblocks, ncores);
	if(maxthreads&&nthreads>maxthreads)
		nthreads=maxthreads;
	//if(nthreads0)
	//{
	//	int nthreads2=MINVAR(nblocks, ncores);
	//	nthreads=nthreads0;
	//	CLAMP2(nthreads, 1, nthreads2);
	//}
	//else
	//	nthreads=MINVAR(nblocks, ncores);
	blocksperthread=(nblocks+nthreads-1)/nthreads;
	coffset=(int)sizeof(int)*nblocks;
	int *offsets=(int*)malloc(coffset+sizeof(int));
	if(!offsets)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	start=0;
	memusage=0;
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
			offsets[kt]=(int)start;
			start+=size;
		}
		offsets[nblocks]=(int)start;
		if(image+start!=imageend)
			LOG_ERROR("Corrupt file");

		image2=(unsigned char*)malloc(usize);
		if(!image2)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(image2, 0, usize);
	}
	argssize=nthreads*sizeof(ThreadArgs);
	args=(ThreadArgs*)malloc(argssize);
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	bestsize=0;
	memusage+=argssize;
	memset(args, 0, argssize);
	histsize=(int)sizeof(int[OCH_COUNT*PRED_COUNT<<8]);
	for(int k=0, bidx=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=image;
		arg->dst=fwd?0:dst->data+printed;
		arg->iw=iw;
		arg->ih=ih;
		
		arg->xblocks=xblocks;
		arg->b1=bidx;
		bidx+=blocksperthread;
		if(bidx>nblocks)
			bidx=nblocks;
		arg->b2=bidx;
		arg->offsets=offsets;
		if(!fwd)
			arg->decsrc=image;

		int listssize=((size_t)arg->b2-arg->b1)*sizeof(BList);
		arg->lists=(BList*)malloc(listssize);

		arg->bufsize=sizeof(short[4*OCH_COUNT*2*(BLOCKX+16LL)]);//4 padded rows * OCH_COUNT * {pixels, wg_errors}
		arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));

		arg->ebufsize=(int)sizeof(int[4][4][WG_NPREDS][BLOCKX+16LL]);
		arg->ebuf=(int*)_mm_malloc(arg->ebufsize, sizeof(__m256i));

		arg->histsize=histsize;
		arg->hist=(int*)malloc(histsize);
		if(!arg->lists||!arg->pixels||!arg->ebuf||!arg->hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=listssize;
		memusage+=arg->ebufsize;
		memusage+=arg->bufsize;
		memusage+=histsize;

		arg->fwd=fwd;
		arg->test=test;
#ifdef ENABLE_MT
		arg->loud=0;
#else
		arg->loud=0;
	//	arg->loud=test&&nblocks<MAXPRINTEDBLOCKS;
#endif
#ifdef ENABLE_ZIPF_VIEW
		arg->zimage=zimage;
#endif
	}
	for(int k2=0;k2<=test;++k2)
	{
#ifdef ENABLE_MT
		if(nthreads>1)
		{
			void *ctx=mt_exec(block_manager, args, sizeof(ThreadArgs), nthreads);
			mt_finish(ctx);
		}
		else
#endif
		{
			for(int k=0;k<nthreads;++k)
				block_manager(args+k);
		}
		if(fwd)
		{
			for(int kt1=0;kt1<nthreads;++kt1)
			{
				ThreadArgs *arg=args+kt1;
				for(int kt2=0;kt2<blocksperthread;++kt2)
				{
					int bidx=kt1*blocksperthread+kt2;
					if(bidx>=nblocks)
						break;
					memcpy(dst->data+printed+sizeof(int)*bidx, &arg->lists[kt2].nbytes, sizeof(int));
					blist_appendtoarray(arg->lists+kt2, &dst);
					blist_clear(arg->lists+kt2);
				}
				if(fwd)
					bestsize+=arg->bestsize;
			}
		}
		if(test)
		{
			ptrdiff_t csize=dst->count;
			t0=time_sec()-t0;
			if(fwd)
			{
#ifdef ESTIMATE_SIZE
				esize[0]/=8;
				esize[1]/=8;
				esize[2]/=8;
				printf("T %15.2lf\n", esize[0]+esize[1]+esize[2]);
				printf("Y %15.2lf\n", esize[0]);
				printf("U %15.2lf\n", esize[1]);
				printf("V %15.2lf\n", esize[2]);
#endif
				printf("%15.2lf (%+13.2lf) bytes\n", bestsize, csize-bestsize);
				printf("%12td/%12d  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
				printf("Mem usage: ");
				print_size((double)memusage, 8, 4, 0, 0);
				printf("\n");
			}
			printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
			if(!fwd)
				compare_bufs_8(image2, src->data+headersize, iw, ih, nch, nch, CODECNAME, 0, 1);
		}
		if(fwd&&test)//transition to (test) decode
		{
			fwd=0;
			image2=(unsigned char*)malloc(usize);
			if(!image2)
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			memset(image2, 0, usize);
			start=coffset;
			for(int kt=0;kt<nblocks;++kt)
			{
				int size=0;
				memcpy(&size, dst->data+printed+sizeof(int)*kt, sizeof(int));
				offsets[kt]=(int)start;
				start+=size;
			}
			offsets[nblocks]=(int)start;
			for(int kt=0;kt<nthreads;++kt)
			{
				ThreadArgs *arg=args+kt;
				arg->dst=image2;
				arg->fwd=0;
				arg->decsrc=dst->data+printed;
			}
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
		free(arg->lists);
		free(arg->hist);
		_mm_free(arg->pixels);
		_mm_free(arg->ebuf);
	}
	free(args);
	free(offsets);
	array_free(&src);
	array_free(&dst);
#ifdef ENABLE_ZIPF_VIEW
	if(test)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H%M%S.PGM");
		FILE *f=fopen(g_buf, "wb");
		if(!f)
		{
			LOG_ERROR("Cannot open \"%s\" for writing", g_buf);
			return 1;
		}
		fprintf(f, "P5\n%d %d\n255\n", iw, ih);
		fwrite(zimage, 1, (size_t)iw*ih, f);
		fclose(f);
		free(zimage);
	}
#endif
	return 0;
}
