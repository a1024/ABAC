#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<process.h>
#include<immintrin.h>
static const char file[]=__FILE__;

//WG:
#define WG_NBITS 0	//why 0 is best?

#if 1
#define WG_NPREDS	12	//multiple of 4
#define WG_PREDLIST0\
	WG_PRED(210,	N+(23*eN-2*(eNN+eNW)+9*eW)/110)\
	WG_PRED(210,	W+(23*eW-2*(eWW+eNW)+9*eN)/110)\
	WG_PRED(215,	3*(N-NN)+NNN+(eN/2+eNN/2+eNNN)/3)\
	WG_PRED(215,	3*(W-WW)+WWW+(eW/2+eWW/2+eWWW)/3)\
	WG_PRED(140,	W+NE-N+((-13*eN)>>4)+eW*13/32-(eW>>7))\
	WG_PRED(230,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(N-W)/6+(NNN-NE-(5*(eN+eW)+eWW))/2)/3)\
	WG_PRED(120,	N+W-NW+(eN+eW-eNW)/3)\
	WG_PRED(140,	N+NE-NNE+((eN+eNE+eNNE+4)>>3))\
	WG_PRED(45,	(4*(N+NNN)-6*NN+NNNW+NNNE-(NNWW+NNEE)/2+NNE+NNW-NE-NW-eN-eNN+eNNN)/3)\
	WG_PRED(97,	(6*(W+WWW)+20*WW+(eW-eWW+eWWW)/3)/32)\
	WG_PRED(65,	(W+3*NW-NWW-NNWW)/2+eNW/4+eW/6)\
	WG_PRED(40,	(3*NE+NEE+NEEEE-NNEE-NNEEE+(3*eNE+6*eNEE+3*eNEEE)/2)/3)
#define WG_PREDLIST1\
	WG_PRED(250,	N+(3*eN+eNW+eW)/6)\
	WG_PRED(250,	W+(3*eW+eNW+eN)/6)\
	WG_PRED(175,	3*(N-NN)+NNN+(eN/2+eNN/2-eWW)/2)\
	WG_PRED(175,	3*(W-WW)+WWW+(eW/2+eWW/2-eNN)/2)\
	WG_PRED(180,	W+NE-N-((eN+eW+31)>>5))\
	WG_PRED(175,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(W-N+NNN-NE)/2-eN-eW-eWW/3)/3)\
	WG_PRED(130,	N+W-NW+(2*(eN+eW)-eNW)/5)\
	WG_PRED(150,	N+NE-NNE+((eN+eNE+eNNE+8)>>4))\
	WG_PRED(45,	(4*(N+NNN)-6*NN+NNNW+NNNE-(NNWW+NNEE)/2+NNE+NNW-NE-NW-2*eN-eNN+eNNN)/3)\
	WG_PRED(57,	(W+WW+(eW-eWW+eWWW)/4)/2)\
	WG_PRED(35,	(W+3*NW-NWW-NNWW+eNW)/2+eW/3)\
	WG_PRED(40,	(3*NE+NEE+NEEEE-NNEE-NNEEE+(3*eNE+6*eNEE+3*eNEEE)/2)/3)
#define WG_PREDLIST2\
	WG_PRED(270,	N+(3*eN+eW)/6)\
	WG_PRED(270,	W+(3*eW+eN)/6)\
	WG_PRED(200,	3*(N-NN)+NNN+(eN/2-eWW)/2)\
	WG_PRED(200,	3*(W-WW)+WWW+(eW/2-eNN)/2)\
	WG_PRED(180,	W+NE-N-((2*eN+eW+31)>>5))\
	WG_PRED(165,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(W-N+NNN-NE)/2-eN-eW-eWW/3)/3)\
	WG_PRED(140,	N+W-NW+(2*(eN+eW)-eNW)/5)\
	WG_PRED(150,	N+NE-NNE+(eN+eNE+eNNE)/16)\
	WG_PRED(55,	(4*(N+NNN)-6*NN+NNNW+NNNE-(NNWW+NNEE)/2+NNE+NNW-NE-NW)/3-eN-eNN+eNNN)\
	WG_PRED(47,	(W+WW+(eW+eWWW)/3)/2)\
	WG_PRED(22,	(W+3*NW-NWW-NNWW+eNW)/2+eW/3)\
	WG_PRED(40,	(3*NE+NEE+NEEEE-NNEE-NNEEE+(3*eNE+6*eNEE+3*eNEEE)/2)/3)
#endif
#if 0
#define WG_NPREDS	8	//multiple of 4
#define WG_PREDLIST0\
	WG_PRED(210,	N+(23*eN-2*(eNN+eNW)+9*eW)/110)\
	WG_PRED(210,	W+(23*eW-2*(eWW+eNW)+9*eN)/110)\
	WG_PRED(215,	3*(N-NN)+NNN+(eN/2+eNN/2+eNNN)/3)\
	WG_PRED(215,	3*(W-WW)+WWW+(eW/2+eWW/2+eWWW)/3)\
	WG_PRED(140,	W+NE-N+((-13*eN)>>4)+eW*13/32-(eW>>7))\
	WG_PRED(230,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(N-W)/6+(NNN-NE-(5*(eN+eW)+eWW))/2)/3)\
	WG_PRED(120,	N+W-NW+(eN+eW-eNW)/3)\
	WG_PRED(140,	N+NE-NNE+((eN+eNE+eNNE+4)>>3))
#define WG_PREDLIST1\
	WG_PRED(250,	N+(3*eN+eNW+eW)/6)\
	WG_PRED(250,	W+(3*eW+eNW+eN)/6)\
	WG_PRED(175,	3*(N-NN)+NNN+(eN/2+eNN/2-eWW)/2)\
	WG_PRED(175,	3*(W-WW)+WWW+(eW/2+eWW/2-eNN)/2)\
	WG_PRED(180,	W+NE-N-((eN+eW+31)>>5))\
	WG_PRED(175,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(W-N+NNN-NE)/2-eN-eW-eWW/3)/3)\
	WG_PRED(130,	N+W-NW+(2*(eN+eW)-eNW)/5)\
	WG_PRED(150,	N+NE-NNE+((eN+eNE+eNNE+8)>>4))
#define WG_PREDLIST2\
	WG_PRED(270,	N+(3*eN+eW)/6)\
	WG_PRED(270,	W+(3*eW+eN)/6)\
	WG_PRED(200,	3*(N-NN)+NNN+(eN/2-eWW)/2)\
	WG_PRED(200,	3*(W-WW)+WWW+(eW/2-eNN)/2)\
	WG_PRED(180,	W+NE-N-((2*eN+eW+31)>>5))\
	WG_PRED(165,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(W-N+NNN-NE)/2-eN-eW-eWW/3)/3)\
	WG_PRED(140,	N+W-NW+(2*(eN+eW)-eNW)/5)\
	WG_PRED(150,	N+NE-NNE+(eN+eNE+eNNE)/16)
#endif
#if 0
#define WG_NPREDS	20
#define WG_PREDLIST\
	WG_PRED(210,	N+(23*eN-2*(eNN+eNW)+9*eW)/110)\
	WG_PRED(210,	W+(23*eW-2*(eWW+eNW)+9*eN)/110)\
	WG_PRED(215,	3*(N-NN)+NNN+(eN/2+eNN/2+eNNN)/3)\
	WG_PRED(215,	3*(W-WW)+WWW+(eW/2+eWW/2+eWWW)/3)\
	WG_PRED(140,	W+NE-N+((-13*eN)>>4)+eW*13/32-(eW>>7))\
	WG_PRED(230,	(WWW+NN-2*NW+NEE+NEEE+NEEEE+(N-W)/6+(NNN-NE-(5*(eN+eW)+eWW))/2)/3)\
	WG_PRED(120,	N+W-NW+(eN+eW-eNW)/3)\
	WG_PRED(140,	N+NE-NNE+((eN+eNE+eNNE+4)>>3))\
	WG_PRED(50,	3*N-W-NE)\
	WG_PRED(30,	2*NE-W)\
	WG_PRED(132,	N+W-NW)\
	WG_PRED(100,	N+W-NW+((eN+eW-eNW+16)>>5))\
	WG_PRED(100,	N)\
	WG_PRED(100,	W)\
	WG_PRED(165,	W+NE-N)\
	WG_PRED(220,	N+NE-NNE)\
	WG_PRED(25,	N+NE-NNE+((eN+eNE-eNNE)>>2))\
	WG_PRED(165,	3*(N-NN)+NNN)\
	WG_PRED(165,	3*(W-WW)+WWW)\
	WG_PRED(150,	(W+NEEE)/2)
#endif
#if 0
#define WG_PREDLIST\
	WG_PRED(kc==0?230:(kc==1?250:270),	N+(kc2==0?(3*eN+eNW+eW)/16:(kc2==2?(3*eN+eNW+eW)/6:(3*eN+eW)/6)))\
	WG_PRED(kc==0?230:(kc==1?250:270),	W+(kc2==0?(3*eW+eNW+eN)/16:(kc2==2?(3*eW+eNW+eN)/6:(3*eW+eN)/6)))\
	WG_PRED(kc==0?195:(kc==1?185:200),	3*(N-NN)+NNN+(kc2==0?(eN+eNN+eNNN)/3:(kc2==2?(eN+eNN-eWW)/2:(eN-eWW)/2)))\
	WG_PRED(kc==0?195:(kc==1?185:200),	3*(W-WW)+WWW+(kc2==0?(eW+eWW+eWWW)/3:(kc2==2?(eW+eWW-eNN)/2:(eW-eNN)/2)))\
	WG_PRED(185,				W+NE-N-(eW>>5)-(kc2==0?eN:(kc2==2?eN>>5:eN>>4)))\
	WG_PRED(195,				(WWW+(W-N+NN+NNN)/2-NW+NEE+NEEE+NEEEE-eN-eW-eWW)/4)\
	WG_PRED(132,				N+W-NW+(kc2?(2*(eN+eW)-eNW)/5:(eN+eW-eNW)/3))\
	WG_PRED(120,				N+NE-NNE+((eN+eNE+eNNE+(kc2?8:4))>>(kc2?4:3)))\
	WG_PRED(40,				(4*(N+NNN)-6*NN+NNNW+NNNE-(NNWW+NNEE)/2+NNE+NNW-NE-NW-eN)/3)\
	WG_PRED(kc?37:67,			kc2==0 ? (6*(W+WWW)+20*WW+(eW-eWW+eWWW)/3)/32 : (kc2==2 ? (W+WW+(eW-eWW+eWWW)/4)/2 : (W+WW+(eW+eWWW)/3)/2))\
	WG_PRED(40,				(W+3*NW-NWW-NNWW)/2+eW/6+eNW)\
	WG_PRED(30,				(3*NE+NEE+NEEEE-NNEE-NNEEE+(3*eNE+6*eNEE+3*eNEEE)/2)/3)
#endif
//	WG_PRED(0, (W+2*NE-NNE+(eW+2*eNE-eNNE)/4)/2)
//	WG_PRED(0, N+W-NW)
//	WG_PRED(0, (-WW+3*W+N)/3)
//	WG_PRED(0, (4*(N+W)-WW-NN)/6)
//	WG_PRED(0, (10*W+5*(N-WW)-NN+WWW)/10)
//	WG_PRED(0, (15*(N+W)-6*(NN+WW)+NNN+WWW)/20)
//	WG_PRED(0, N+NE-NNE+((eN+eNE-eNNE)>>2))
#if 0
#define WG_NPREDS	16
#define WG_PREDLIST\
	WG_PRED(25, spred)\
	WG_PRED(50, wgrad)\
	WG_PRED(50, 3*N-W-NE)\
	WG_PRED(30, 2*NE-W)\
	WG_PRED(132, N+W-NW)\
	WG_PRED(100, N+W-NW+((eN+eW-eNW+16)>>5))\
	WG_PRED(150, N+eN)\
	WG_PRED(100, N)\
	WG_PRED(150, W+eW)\
	WG_PRED(100, W)\
	WG_PRED(165, W+NE-N)\
	WG_PRED(220, N+NE-NNE)\
	WG_PRED(25, N+NE-NNE+((eN+eNE-eNNE)>>2))\
	WG_PRED(165, 3*(N-NN)+NNN)\
	WG_PRED(165, 3*(W-WW)+WWW)\
	WG_PRED(150, (W+NEEE)/2)
#endif
static void wg_init(double *weights, int kc)
{
	int j=0;
//#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
//	WG_PREDLIST
//#undef  WG_PRED
#if 1
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
#endif
}
FORCE_INLINE int wg_predict(
	const double *weights, const int **rows, const int stride, int kc2,
	int spred, const int *perrors, const int *NWerrors, const int *Nerrors, const int *NEerrors, const int *NNEerrors, int *preds
)
{
	//double fpred=0, wsum=0;
	int j=0;
	ALIGN(16) int pred[4];
	int
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
	//int sh=kc2?1:2;
	//int gy=abs(eN)+1, gx=abs(eW)+1;
	//int wgrad=(N*gy+W*gx)/(gy+gx);

	//int cgrad, cgrad2;
	//MEDIAN3_32(cgrad, N, W, N+W-NW);
	//MEDIAN3_32(cgrad2, N, NE, N+NE-NNE);

//#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
//	WG_PREDLIST
//#undef  WG_PRED
#if 1
	switch(kc2)
	{
	case 0:
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
		WG_PREDLIST0
#undef  WG_PRED
		break;
	case 2:
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
		WG_PREDLIST1
#undef  WG_PRED
		break;
	case 4:
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
		WG_PREDLIST2
#undef  WG_PRED
		break;
	}
#endif
	//if((eW*47>>7)!=(eW*3>>3)-(eW>>7))
	//	LOG_ERROR("%d", eW);
#if 1
	__m128i one=_mm_set1_epi32(1);
	__m128i me0=_mm_load_si128((__m128i*)perrors+0);
	__m128i me1=_mm_load_si128((__m128i*)perrors+1);
	__m128i me2=_mm_load_si128((__m128i*)perrors+2);
//	__m128i me3=_mm_load_si128((__m128i*)perrors+3);
	__m256d w0=_mm256_load_pd(weights+0*4);
	__m256d w1=_mm256_load_pd(weights+1*4);
	__m256d w2=_mm256_load_pd(weights+2*4);
//	__m256d w3=_mm256_load_pd(weights+3*4);
	//me0=_mm_srai_epi32(me0, 1);
	//me1=_mm_srai_epi32(me1, 1);
	//me2=_mm_srai_epi32(me2, 1);
	//me3=_mm_srai_epi32(me3, 1);
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)NWerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)NWerrors+1));
	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)NWerrors+2));
//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)NWerrors+3));
	me0=_mm_add_epi32(me0, _mm_slli_epi32(_mm_load_si128((__m128i*)Nerrors+0), 1));
	me1=_mm_add_epi32(me1, _mm_slli_epi32(_mm_load_si128((__m128i*)Nerrors+1), 1));
	me2=_mm_add_epi32(me2, _mm_slli_epi32(_mm_load_si128((__m128i*)Nerrors+2), 1));
//	me3=_mm_add_epi32(me3, _mm_slli_epi32(_mm_load_si128((__m128i*)Nerrors+3), 1));
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)NEerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)NEerrors+1));
	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)NEerrors+2));
//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)NEerrors+3));
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)NNEerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)NNEerrors+1));
	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)NNEerrors+2));
//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)NNEerrors+3));
#if 1
	me0=_mm_srli_epi32(me0, 1);
	me1=_mm_srli_epi32(me1, 1);
	me2=_mm_srli_epi32(me2, 1);
//	me3=_mm_srli_epi32(me3, 1);
#endif
	me0=_mm_add_epi32(me0, one);
	me1=_mm_add_epi32(me1, one);
	me2=_mm_add_epi32(me2, one);
//	me3=_mm_add_epi32(me3, one);
#if 0
	__m128i emax=_mm_set1_epi32(0x3FFF);
	me0=_mm_min_epi32(me0, emax);
	me1=_mm_min_epi32(me1, emax);
	me2=_mm_min_epi32(me2, emax);
//	me3=_mm_min_epi32(me3, emax);
#endif
	__m256d de0=_mm256_cvtepi32_pd(me0);
	__m256d de1=_mm256_cvtepi32_pd(me1);
	__m256d de2=_mm256_cvtepi32_pd(me2);
//	__m256d de3=_mm256_cvtepi32_pd(me3);
	w0=_mm256_div_pd(w0, de0);
	w1=_mm256_div_pd(w1, de1);
	w2=_mm256_div_pd(w2, de2);
//	w3=_mm256_div_pd(w3, de3);
	de0=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+0));
	de1=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+1));
	de2=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+2));
//	de3=_mm256_cvtepi32_pd(_mm_load_si128((__m128i*)preds+3));
	de0=_mm256_mul_pd(de0, w0);
	de1=_mm256_mul_pd(de1, w1);
	de2=_mm256_mul_pd(de2, w2);
//	de3=_mm256_mul_pd(de3, w3);
	w0=_mm256_add_pd(w0, w1);
	w0=_mm256_add_pd(w0, w2);
//	w0=_mm256_add_pd(w0, w3);
	de0=_mm256_add_pd(de0, de1);
	de0=_mm256_add_pd(de0, de2);
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
	_mm_store_si128((__m128i*)pred, mp);
#else
	for(int k=0;k<WG_NPREDS;++k)
	{
		double weight=(double)weights[k]/(perrors[k]+1);
		fpred+=weight*preds[k];
		wsum+=weight;
	}
	fpred/=wsum;
#ifdef _MSC_VER
	pred=_cvt_dtoi_fast(fpred);
#else
	pred=(int)fpred;
#endif
	CLAMP3_32(pred, pred, N, W, NE);
#endif
	return pred[0];
}
FORCE_INLINE int wg_predict2(
	const double *weights, const int **rows, const int stride, int kc2,
	int spred, const int *perrors, const int *NWerrors, const int *Nerrors, const int *NEerrors, const int *NNEerrors, int *preds
)
{
	int j=0;
	int
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
	
//#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
//	WG_PREDLIST
//#undef  WG_PRED
#if 1
	switch(kc2)
	{
	case 0:
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
		WG_PREDLIST0
#undef  WG_PRED
		break;
	case 2:
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
		WG_PREDLIST1
#undef  WG_PRED
		break;
	case 4:
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
		WG_PREDLIST2
#undef  WG_PRED
		break;
	}
#endif
	if(N==W)
		return N;

	int coeffs[]=
	{
		perrors[ 0]+NWerrors[ 0]+(Nerrors[ 0]<<1)+NEerrors[ 0]+NNEerrors[ 0],
		perrors[ 1]+NWerrors[ 1]+(Nerrors[ 1]<<1)+NEerrors[ 1]+NNEerrors[ 1],
		perrors[ 2]+NWerrors[ 2]+(Nerrors[ 2]<<1)+NEerrors[ 2]+NNEerrors[ 2],
		perrors[ 3]+NWerrors[ 3]+(Nerrors[ 3]<<1)+NEerrors[ 3]+NNEerrors[ 3],
		perrors[ 4]+NWerrors[ 4]+(Nerrors[ 4]<<1)+NEerrors[ 4]+NNEerrors[ 4],
		perrors[ 5]+NWerrors[ 5]+(Nerrors[ 5]<<1)+NEerrors[ 5]+NNEerrors[ 5],
		perrors[ 6]+NWerrors[ 6]+(Nerrors[ 6]<<1)+NEerrors[ 6]+NNEerrors[ 6],
		perrors[ 7]+NWerrors[ 7]+(Nerrors[ 7]<<1)+NEerrors[ 7]+NNEerrors[ 7],
		perrors[ 8]+NWerrors[ 8]+(Nerrors[ 8]<<1)+NEerrors[ 8]+NNEerrors[ 8],
		perrors[ 9]+NWerrors[ 9]+(Nerrors[ 9]<<1)+NEerrors[ 9]+NNEerrors[ 9],
		perrors[10]+NWerrors[10]+(Nerrors[10]<<1)+NEerrors[10]+NNEerrors[10],
		perrors[11]+NWerrors[11]+(Nerrors[11]<<1)+NEerrors[11]+NNEerrors[11],
	//	perrors[12]+NWerrors[12]+(Nerrors[12]<<1)+NEerrors[12]+NNEerrors[12],
	//	perrors[13]+NWerrors[13]+(Nerrors[13]<<1)+NEerrors[13]+NNEerrors[13],
	//	perrors[14]+NWerrors[14]+(Nerrors[14]<<1)+NEerrors[14]+NNEerrors[14],
	//	perrors[15]+NWerrors[15]+(Nerrors[15]<<1)+NEerrors[15]+NNEerrors[15],
	//	perrors[16]+NWerrors[16]+(Nerrors[16]<<1)+NEerrors[16]+NNEerrors[16],
	//	perrors[17]+NWerrors[17]+(Nerrors[17]<<1)+NEerrors[17]+NNEerrors[17],
	//	perrors[18]+NWerrors[18]+(Nerrors[18]<<1)+NEerrors[18]+NNEerrors[18],
	//	perrors[19]+NWerrors[19]+(Nerrors[19]<<1)+NEerrors[19]+NNEerrors[19],
	};
	int cmax=coeffs[0];
	if(cmax<coeffs[ 1])cmax=coeffs[ 1];
	if(cmax<coeffs[ 2])cmax=coeffs[ 2];
	if(cmax<coeffs[ 3])cmax=coeffs[ 3];
	if(cmax<coeffs[ 4])cmax=coeffs[ 4];
	if(cmax<coeffs[ 5])cmax=coeffs[ 5];
	if(cmax<coeffs[ 6])cmax=coeffs[ 6];
	if(cmax<coeffs[ 7])cmax=coeffs[ 7];
	if(cmax<coeffs[ 8])cmax=coeffs[ 8];
	if(cmax<coeffs[ 9])cmax=coeffs[ 9];
	if(cmax<coeffs[10])cmax=coeffs[10];
	if(cmax<coeffs[11])cmax=coeffs[11];
	//if(cmax<coeffs[12])cmax=coeffs[12];
	//if(cmax<coeffs[13])cmax=coeffs[13];
	//if(cmax<coeffs[14])cmax=coeffs[14];
	//if(cmax<coeffs[15])cmax=coeffs[15];
	//if(cmax<coeffs[16])cmax=coeffs[16];
	//if(cmax<coeffs[17])cmax=coeffs[17];
	//if(cmax<coeffs[18])cmax=coeffs[18];
	//if(cmax<coeffs[19])cmax=coeffs[19];
	++cmax;
	coeffs[ 0]=cmax-coeffs[ 0];
	coeffs[ 1]=cmax-coeffs[ 1];
	coeffs[ 2]=cmax-coeffs[ 2];
	coeffs[ 3]=cmax-coeffs[ 3];
	coeffs[ 4]=cmax-coeffs[ 4];
	coeffs[ 5]=cmax-coeffs[ 5];
	coeffs[ 6]=cmax-coeffs[ 6];
	coeffs[ 7]=cmax-coeffs[ 7];
	coeffs[ 8]=cmax-coeffs[ 8];
	coeffs[ 9]=cmax-coeffs[ 9];
	coeffs[10]=cmax-coeffs[10];
	coeffs[11]=cmax-coeffs[11];
	//coeffs[12]=cmax-coeffs[12];
	//coeffs[13]=cmax-coeffs[13];
	//coeffs[14]=cmax-coeffs[14];
	//coeffs[15]=cmax-coeffs[15];
	//coeffs[16]=cmax-coeffs[16];
	//coeffs[17]=cmax-coeffs[17];
	//coeffs[18]=cmax-coeffs[18];
	//coeffs[19]=cmax-coeffs[19];
	coeffs[ 0]*=coeffs[ 0];
	coeffs[ 1]*=coeffs[ 1];
	coeffs[ 2]*=coeffs[ 2];
	coeffs[ 3]*=coeffs[ 3];
	coeffs[ 4]*=coeffs[ 4];
	coeffs[ 5]*=coeffs[ 5];
	coeffs[ 6]*=coeffs[ 6];
	coeffs[ 7]*=coeffs[ 7];
	coeffs[ 8]*=coeffs[ 8];
	coeffs[ 9]*=coeffs[ 9];
	coeffs[10]*=coeffs[10];
	coeffs[11]*=coeffs[11];
	//coeffs[12]*=coeffs[12];
	//coeffs[13]*=coeffs[13];
	//coeffs[14]*=coeffs[14];
	//coeffs[15]*=coeffs[15];
	//coeffs[16]*=coeffs[16];
	//coeffs[17]*=coeffs[17];
	//coeffs[18]*=coeffs[18];
	//coeffs[19]*=coeffs[19];
	int pred=(int)((
		+(long long)coeffs[ 0]*preds[ 0]
		+(long long)coeffs[ 1]*preds[ 1]
		+(long long)coeffs[ 2]*preds[ 2]
		+(long long)coeffs[ 3]*preds[ 3]
		+(long long)coeffs[ 4]*preds[ 4]
		+(long long)coeffs[ 5]*preds[ 5]
		+(long long)coeffs[ 6]*preds[ 6]
		+(long long)coeffs[ 7]*preds[ 7]
		+(long long)coeffs[ 8]*preds[ 8]
		+(long long)coeffs[ 9]*preds[ 9]
		+(long long)coeffs[10]*preds[10]
		+(long long)coeffs[11]*preds[11]
	//	+(long long)coeffs[12]*preds[12]
	//	+(long long)coeffs[13]*preds[13]
	//	+(long long)coeffs[14]*preds[14]
	//	+(long long)coeffs[15]*preds[15]
	//	+(long long)coeffs[16]*preds[16]
	//	+(long long)coeffs[17]*preds[17]
	//	+(long long)coeffs[18]*preds[18]
	//	+(long long)coeffs[19]*preds[19]
	)/(
		+(long long)coeffs[ 0]
		+(long long)coeffs[ 1]
		+(long long)coeffs[ 2]
		+(long long)coeffs[ 3]
		+(long long)coeffs[ 4]
		+(long long)coeffs[ 5]
		+(long long)coeffs[ 6]
		+(long long)coeffs[ 7]
		+(long long)coeffs[ 8]
		+(long long)coeffs[ 9]
		+(long long)coeffs[10]
		+(long long)coeffs[11]
	//	+(long long)coeffs[12]
	//	+(long long)coeffs[13]
	//	+(long long)coeffs[14]
	//	+(long long)coeffs[15]
	//	+(long long)coeffs[16]
	//	+(long long)coeffs[17]
	//	+(long long)coeffs[18]
	//	+(long long)coeffs[19]
	));
#if 0
	__m128i one=_mm_set1_epi32(1);
	__m128i me0=_mm_load_si128((__m128i*)perrors+0);
	__m128i me1=_mm_load_si128((__m128i*)perrors+1);
	__m128i me2=_mm_load_si128((__m128i*)perrors+2);
//	__m128i me3=_mm_load_si128((__m128i*)perrors+3);
//	__m256d w0=_mm256_load_pd(weights+0*4);
//	__m256d w1=_mm256_load_pd(weights+1*4);
//	__m256d w2=_mm256_load_pd(weights+2*4);
//	__m256d w3=_mm256_load_pd(weights+3*4);
	//me0=_mm_srai_epi32(me0, 1);
	//me1=_mm_srai_epi32(me1, 1);
	//me2=_mm_srai_epi32(me2, 1);
	//me3=_mm_srai_epi32(me3, 1);
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)NWerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)NWerrors+1));
	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)NWerrors+2));
//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)NWerrors+3));
	me0=_mm_add_epi32(me0, _mm_slli_epi32(_mm_load_si128((__m128i*)Nerrors+0), 1));
	me1=_mm_add_epi32(me1, _mm_slli_epi32(_mm_load_si128((__m128i*)Nerrors+1), 1));
	me2=_mm_add_epi32(me2, _mm_slli_epi32(_mm_load_si128((__m128i*)Nerrors+2), 1));
//	me3=_mm_add_epi32(me3, _mm_slli_epi32(_mm_load_si128((__m128i*)Nerrors+3), 1));
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)NEerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)NEerrors+1));
	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)NEerrors+2));
//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)NEerrors+3));
	me0=_mm_add_epi32(me0, _mm_load_si128((__m128i*)NNEerrors+0));
	me1=_mm_add_epi32(me1, _mm_load_si128((__m128i*)NNEerrors+1));
	me2=_mm_add_epi32(me2, _mm_load_si128((__m128i*)NNEerrors+2));
//	me3=_mm_add_epi32(me3, _mm_load_si128((__m128i*)NNEerrors+3));
	me0=_mm_add_epi32(me0, one);
	me1=_mm_add_epi32(me1, one);
	me2=_mm_add_epi32(me2, one);
//	me3=_mm_add_epi32(me3, one);

	__m128i fbias=_mm_set1_epi32(127);
	me0=_mm_castps_si128(_mm_cvtepi32_ps(me0));//FLOOR_LOG2
	me1=_mm_castps_si128(_mm_cvtepi32_ps(me1));
	me2=_mm_castps_si128(_mm_cvtepi32_ps(me2));
	me0=_mm_srli_epi32(me0, 23);
	me1=_mm_srli_epi32(me1, 23);
	me2=_mm_srli_epi32(me2, 23);
	me0=_mm_sub_epi32(me0, fbias);
	me1=_mm_sub_epi32(me1, fbias);
	me2=_mm_sub_epi32(me2, fbias);
	ALIGN(16) int esums[WG_NPREDS];
	_mm_store_si128((__m128i*)esums+0, me0);
	_mm_store_si128((__m128i*)esums+1, me1);
	_mm_store_si128((__m128i*)esums+2, me2);
	int t2=esums[0]+esums[1]+esums[2];
	int t3=esums[3]+esums[4]+esums[5];
	int t4=esums[6]+esums[7]+esums[8];
	int t5=esums[9]+esums[10]+esums[11];
	int t0=t2+t3;
	int t1=t4+t5;
	int lweights[WG_NPREDS]=
	{
		esums[1]+esums[2]+(t3+t1),
		esums[0]+esums[2]+(t3+t1),
		esums[0]+esums[1]+(t3+t1),
		(t2+t1)+esums[4]+esums[5],
		(t2+t1)+esums[3]+esums[5],
		(t2+t1)+esums[3]+esums[4],
		esums[7]+esums[8]+(t0+t5),
		esums[6]+esums[8]+(t0+t5),
		esums[6]+esums[7]+(t0+t5),
		(t0+t4)+esums[10]+esums[11],
		(t0+t4)+esums[ 9]+esums[11],
		(t0+t4)+esums[ 9]+esums[10],

		//esums[1]+esums[2]+(t3+t1)+1,
		//esums[0]+esums[2]+(t3+t1)+1,
		//esums[0]+esums[1]+(t3+t1)+1,
		//(t2+t1)+esums[4]+esums[5]+1,
		//(t2+t1)+esums[3]+esums[5]+1,
		//(t2+t1)+esums[3]+esums[4]+1,
		//esums[7]+esums[8]+(t0+t5)+1,
		//esums[6]+esums[8]+(t0+t5)+1,
		//esums[6]+esums[7]+(t0+t5)+1,
		//(t0+t4)+esums[10]+esums[11]+1,
		//(t0+t4)+esums[ 9]+esums[11]+1,
		//(t0+t4)+esums[ 9]+esums[10]+1,

		//esums[1]+esums[2]+(esums[3]+esums[4]+esums[5])+(esums[6]+esums[7]+esums[8]+esums[9]+esums[10]+esums[11]),
		//esums[0]+esums[2]+(esums[3]+esums[4]+esums[5])+(esums[6]+esums[7]+esums[8]+esums[9]+esums[10]+esums[11]),
		//esums[0]+esums[1]+(esums[3]+esums[4]+esums[5])+(esums[6]+esums[7]+esums[8]+esums[9]+esums[10]+esums[11]),
		//(esums[0]+esums[1]+esums[2])+esums[4]+esums[5]+(esums[6]+esums[7]+esums[8]+esums[9]+esums[10]+esums[11]),
		//(esums[0]+esums[1]+esums[2])+esums[3]+esums[5]+(esums[6]+esums[7]+esums[8]+esums[9]+esums[10]+esums[11]),
		//(esums[0]+esums[1]+esums[2])+esums[3]+esums[4]+(esums[6]+esums[7]+esums[8]+esums[9]+esums[10]+esums[11]),
		//(esums[0]+esums[1]+esums[2]+esums[3]+esums[4]+esums[5])+esums[7]+esums[8]+(esums[9]+esums[10]+esums[11]),
		//(esums[0]+esums[1]+esums[2]+esums[3]+esums[4]+esums[5])+esums[6]+esums[8]+(esums[9]+esums[10]+esums[11]),
		//(esums[0]+esums[1]+esums[2]+esums[3]+esums[4]+esums[5])+esums[6]+esums[7]+(esums[9]+esums[10]+esums[11]),
		//(esums[0]+esums[1]+esums[2]+esums[3]+esums[4]+esums[5])+(esums[6]+esums[7]+esums[8])+esums[10]+esums[11],
		//(esums[0]+esums[1]+esums[2]+esums[3]+esums[4]+esums[5])+(esums[6]+esums[7]+esums[8])+esums[ 9]+esums[11],
		//(esums[0]+esums[1]+esums[2]+esums[3]+esums[4]+esums[5])+(esums[6]+esums[7]+esums[8])+esums[ 9]+esums[10],
	};
	int wsum=
		+lweights[ 0]
		+lweights[ 1]
		+lweights[ 2]
		+lweights[ 3]
		+lweights[ 4]
		+lweights[ 5]
		+lweights[ 6]
		+lweights[ 7]
		+lweights[ 8]
		+lweights[ 9]
		+lweights[10]
		+lweights[11]
	;
//	int wmin=lweights[0];
//	int wmax=lweights[0];
//	if(wmin>lweights[ 1])wmin=lweights[ 1];
//	if(wmax<lweights[ 1])wmax=lweights[ 1];
//	if(wmin>lweights[ 2])wmin=lweights[ 2];
//	if(wmax<lweights[ 2])wmax=lweights[ 2];
//	if(wmin>lweights[ 3])wmin=lweights[ 3];
//	if(wmax<lweights[ 3])wmax=lweights[ 3];
//	if(wmin>lweights[ 4])wmin=lweights[ 4];
//	if(wmax<lweights[ 4])wmax=lweights[ 4];
//	if(wmin>lweights[ 5])wmin=lweights[ 5];
//	if(wmax<lweights[ 5])wmax=lweights[ 5];
//	if(wmin>lweights[ 6])wmin=lweights[ 6];
//	if(wmax<lweights[ 6])wmax=lweights[ 6];
//	if(wmin>lweights[ 7])wmin=lweights[ 7];
//	if(wmax<lweights[ 7])wmax=lweights[ 7];
//	if(wmin>lweights[ 8])wmin=lweights[ 8];
//	if(wmax<lweights[ 8])wmax=lweights[ 8];
//	if(wmin>lweights[ 9])wmin=lweights[ 9];
//	if(wmax<lweights[ 9])wmax=lweights[ 9];
//	if(wmin>lweights[10])wmin=lweights[10];
//	if(wmax<lweights[10])wmax=lweights[10];
//	if(wmin>lweights[11])wmin=lweights[11];
//	if(wmax<lweights[11])wmax=lweights[11];
	int sh=FLOOR_LOG2(wsum);
	lweights[ 0]<<=3;
	lweights[ 1]<<=3;
	lweights[ 2]<<=3;
	lweights[ 3]<<=3;
	lweights[ 4]<<=3;
	lweights[ 5]<<=3;
	lweights[ 6]<<=3;
	lweights[ 7]<<=3;
	lweights[ 8]<<=3;
	lweights[ 9]<<=3;
	lweights[10]<<=3;
	lweights[11]<<=3;
	lweights[ 0]>>=sh;
	lweights[ 1]>>=sh;
	lweights[ 2]>>=sh;
	lweights[ 3]>>=sh;
	lweights[ 4]>>=sh;
	lweights[ 5]>>=sh;
	lweights[ 6]>>=sh;
	lweights[ 7]>>=sh;
	lweights[ 8]>>=sh;
	lweights[ 9]>>=sh;
	lweights[10]>>=sh;
	lweights[11]>>=sh;
	if(lweights[ 0]>8)lweights[ 0]=8;
	if(lweights[ 1]>8)lweights[ 1]=8;
	if(lweights[ 2]>8)lweights[ 2]=8;
	if(lweights[ 3]>8)lweights[ 3]=8;
	if(lweights[ 4]>8)lweights[ 4]=8;
	if(lweights[ 5]>8)lweights[ 5]=8;
	if(lweights[ 6]>8)lweights[ 6]=8;
	if(lweights[ 7]>8)lweights[ 7]=8;
	if(lweights[ 8]>8)lweights[ 8]=8;
	if(lweights[ 9]>8)lweights[ 9]=8;
	if(lweights[10]>8)lweights[10]=8;
	if(lweights[11]>8)lweights[11]=8;
	int pred=(int)((
		+((long long)preds[ 0]<<lweights[ 0])
		+((long long)preds[ 1]<<lweights[ 1])
		+((long long)preds[ 2]<<lweights[ 2])
		+((long long)preds[ 3]<<lweights[ 3])
		+((long long)preds[ 4]<<lweights[ 4])
		+((long long)preds[ 5]<<lweights[ 5])
		+((long long)preds[ 6]<<lweights[ 6])
		+((long long)preds[ 7]<<lweights[ 7])
		+((long long)preds[ 8]<<lweights[ 8])
		+((long long)preds[ 9]<<lweights[ 9])
		+((long long)preds[10]<<lweights[10])
		+((long long)preds[11]<<lweights[11])
	)/(
		+(1ULL<<lweights[ 0])
		+(1ULL<<lweights[ 1])
		+(1ULL<<lweights[ 2])
		+(1ULL<<lweights[ 3])
		+(1ULL<<lweights[ 4])
		+(1ULL<<lweights[ 5])
		+(1ULL<<lweights[ 6])
		+(1ULL<<lweights[ 7])
		+(1ULL<<lweights[ 8])
		+(1ULL<<lweights[ 9])
		+(1ULL<<lweights[10])
		+(1ULL<<lweights[11])
	));
	//int pred=(
	//	+lweights[ 0]*preds[ 0]
	//	+lweights[ 1]*preds[ 1]
	//	+lweights[ 2]*preds[ 2]
	//	+lweights[ 3]*preds[ 3]
	//	+lweights[ 4]*preds[ 4]
	//	+lweights[ 5]*preds[ 5]
	//	+lweights[ 6]*preds[ 6]
	//	+lweights[ 7]*preds[ 7]
	//	+lweights[ 8]*preds[ 8]
	//	+lweights[ 9]*preds[ 9]
	//	+lweights[10]*preds[10]
	//	+lweights[11]*preds[11]
	//)/((t0+t1+WG_NPREDS)*(WG_NPREDS-1));
#endif
	int vmax=N, vmin=W;
	if(N<W)vmin=N, vmax=W;
	if(vmin>NE)vmin=NE;
	if(vmax<NE)vmax=NE;
	CLAMP2(pred, vmin, vmax);
	return pred;
}
FORCE_INLINE void wg_update(int curr, int kc, const int *preds, int *perrors, int *Werrors, int *currerrors, int *NEerrors)
{
	static const int factors[]={97, 99, 99};
	int factor=factors[kc];
	for(int k=0;k<WG_NPREDS;++k)
	{
		int e2=abs(curr-preds[k])<<1;
		perrors[k]=(perrors[k]+e2)*factor>>7;
	//	perrors[k]=(e2-perrors[k]+(1<<5>>1))>>5;//slightly worse
		currerrors[k]=(2*Werrors[k]+e2+NEerrors[k])>>2;
		NEerrors[k]+=e2;
	}
}

FORCE_INLINE void wg_init_v3(int *weights)
{
//	int j=0;
//#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
//	WG_PREDLIST0
//#undef  WG_PRED

	for(int k=0;k<WG_NPREDS+1;++k)
		weights[k]=k*(0x10000/WG_NPREDS);
}
FORCE_INLINE int wg_predict_v3(const int *weights, const int *currptr, const int *Nptr, const int *NNptr, const int *NNNptr, int *preds)
{
	int
		NNNWWWW	=NNNptr	[-4*4*2+0],
		NNNWWW	=NNNptr	[-3*4*2+0],
		NNNW	=NNNptr	[-1*4*2+0],
		NNN	=NNNptr	[+0*4*2+0],
		NNNE	=NNNptr	[+1*4*2+0],
		NNNEEE	=NNNptr	[+3*4*2+0],
		NNWWWW	=NNptr	[-4*4*2+0],
		NNWW	=NNptr	[-2*4*2+0],
		NNW	=NNptr	[-1*4*2+0],
		NN	=NNptr	[+0*4*2+0],
		NNE	=NNptr	[+1*4*2+0],
		NNEE	=NNptr	[+2*4*2+0],
		NNEEE	=NNptr	[+3*4*2+0],
		NNEEEE	=NNptr	[+4*4*2+0],
		NWWW	=Nptr	[-3*4*2+0],
		NWW	=Nptr	[-2*4*2+0],
		NW	=Nptr	[-1*4*2+0],
		N	=Nptr	[+0*4*2+0],
		NE	=Nptr	[+1*4*2+0],
		NEE	=Nptr	[+2*4*2+0],
		NEEE	=Nptr	[+3*4*2+0],
		NEEEE	=Nptr	[+4*4*2+0],
		WWWWW	=currptr[-5*4*2+0],
		WWWW	=currptr[-4*4*2+0],
		WWW	=currptr[-3*4*2+0],
		WW	=currptr[-2*4*2+0],
		W	=currptr[-1*4*2+0],
		eNNN	=NNNptr	[+0*4*2+1],
		eNN	=NNptr	[+0*4*2+1],
		eNNE	=NNptr	[+1*4*2+1],
		eNW	=Nptr	[-1*4*2+1],
		eN	=Nptr	[+0*4*2+1],
		eNE	=Nptr	[+1*4*2+1],
		eNEE	=Nptr	[+2*4*2+1],
		eNEEE	=Nptr	[+3*4*2+1],
		eWWWW	=currptr[-4*4*2+1],
		eWWW	=currptr[-3*4*2+1],
		eWW	=currptr[-2*4*2+1],
		eW	=currptr[-1*4*2+1];
	int j=0;

#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST0
#undef  WG_PRED
	if(N==W)
		return N;

	int pred=(
		+(weights[ 1]-weights[ 0])*preds[ 0]
		+(weights[ 2]-weights[ 1])*preds[ 1]
		+(weights[ 3]-weights[ 2])*preds[ 2]
		+(weights[ 4]-weights[ 3])*preds[ 3]
		+(weights[ 5]-weights[ 4])*preds[ 4]
		+(weights[ 6]-weights[ 5])*preds[ 5]
		+(weights[ 7]-weights[ 6])*preds[ 6]
		+(weights[ 8]-weights[ 7])*preds[ 7]
		+(weights[ 9]-weights[ 8])*preds[ 8]
		+(weights[10]-weights[ 9])*preds[ 9]
		+(weights[11]-weights[10])*preds[10]
		+(weights[12]-weights[11])*preds[11]
		+(1<<16>>1)
	)>>16;
	{
		int vmax=N, vmin=W;
		if(N<W)vmin=N, vmin=W;
		if(vmin>NE)vmin=NE;
		if(vmax>NE)vmax=NE;
		CLAMP2(pred, vmin, vmax);
	}
	return pred;
}
FORCE_INLINE void wg_update_v3(int target, const int *preds, int *weights)
{
	int mixin[WG_NPREDS], wsum=0;
	for(int k=0;k<WG_NPREDS;++k)
	{
		int e=abs(target-preds[k])+1;
		e=0x1000000/e;//(NPREDS*2) DIVs per subpixel is too slow
		mixin[k]=e;
		wsum+=e;
	}
	int sum2=0;
	for(int k=0;k<WG_NPREDS;++k)
	{
		int val=mixin[k];
		weights[k]+=((int)(((long long)sum2<<16)/wsum)-weights[k])>>0;
	//	mixin[k]=(int)(((long long)sum2<<16)/wsum);
		sum2+=val;
	}

	//for(int k=0;k<WG_NPREDS;++k)//slow O(N^2)	X  inefficient
	//{
	//	int e=target-preds[k];
	//	e*=e;//instead of abs()
	//	e=FLOOR_LOG2(e+1);
	//	for(int k2=1;k2<WG_NPREDS;++k2)
	//		weights[k2]+=(((k2>k)<<16)-weights[k2])>>e;//FIXME SIMD
	//}
}
void pred_wgrad4(Image *src, int fwd)
{
	int amin[]=
	{
		-(1<<src->depth[0]>>1),
		-(1<<src->depth[1]>>1),
		-(1<<src->depth[2]>>1),
		-(1<<src->depth[3]>>1),
	};
	int amax[]=
	{
		(1<<src->depth[0]>>1)-1,
		(1<<src->depth[1]>>1)-1,
		(1<<src->depth[2]>>1)-1,
		(1<<src->depth[3]>>1)-1,
	};
//	ALIGN(32) int wg_weights_v3[4][(WG_NPREDS+7)&~7];
	ALIGN(32) double wg_weights[4][WG_NPREDS]={0};
	ALIGN(32) int wg_perrors[4][WG_NPREDS]={0}, wg_preds[WG_NPREDS]={0};
	int nch;
	int fwdmask=-fwd;
	int bufsize;
	int *pixels;
	
	int invdist=((1<<16)+g_dist-1)/g_dist;
	int pesize=(src->iw+16)*(int)sizeof(int[4][4][WG_NPREDS]);//4 padded rows * 3 channels * WG_NPREDS
	int *ebuf=(int*)_mm_malloc(pesize, sizeof(__m256i));
	bufsize=(src->iw+16LL)*sizeof(int[4*4*2]);//4 padded rows * 4 channels max
	pixels=(int*)malloc(bufsize);
	if(!pixels||!ebuf)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(ebuf, 0, pesize);
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(int kc=0;kc<nch;++kc)
	//	wg_init_v3(wg_weights_v3[kc]);
		wg_init(wg_weights[kc], kc);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		int *erows[]=
		{
			ebuf+((src->iw+16LL)*((ky-0LL)&3)+8LL)*4*WG_NPREDS,
			ebuf+((src->iw+16LL)*((ky-1LL)&3)+8LL)*4*WG_NPREDS,
			ebuf+((src->iw+16LL)*((ky-2LL)&3)+8LL)*4*WG_NPREDS,
			ebuf+((src->iw+16LL)*((ky-3LL)&3)+8LL)*4*WG_NPREDS,
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
#if 0
			//if(ky==src->ih/2&&kx==src->iw/2)
			//	printf("");
			for(int kc=0;kc<src->nch;++kc)
			{
				int pred=wg_predict_v3(wg_weights_v3[kc], rows[0]+kc, rows[1]+kc, rows[2]+kc, rows[3]+kc, wg_preds);
				{
					int curr=src->data[idx+kc], pred0=pred;
					pred+=1<<WG_NBITS>>1;
					pred>>=WG_NBITS;
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=curr;

					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];

					src->data[idx+kc]=pred;
					rows[0][kc+0]=(fwd?curr:pred)<<WG_NBITS;
					rows[0][kc+4]=rows[0][kc]-pred0;
				}
				wg_update_v3(rows[0][kc+0], wg_preds, wg_weights_v3[kc]);
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
#endif
#if 1
			for(int kc=0;kc<src->nch;++kc)
			{
				int pred;
				int
					kc2=kc<<1,
					kc3=kc*WG_NPREDS,
				//	*eNNW	=erows[2]+kc3-1*4*WG_NPREDS,
				//	*eNN	=erows[2]+kc3+0*4*WG_NPREDS,
					*eNNE	=erows[2]+kc3+1*4*WG_NPREDS,
				//	*eNNEE	=erows[2]+kc3+2*4*WG_NPREDS,
					*eNW	=erows[1]+kc3-1*4*WG_NPREDS,
					*eN	=erows[1]+kc3+0*4*WG_NPREDS,
					*eNE	=erows[1]+kc3+1*4*WG_NPREDS,
				//	*eNEE	=erows[1]+kc3+2*4*WG_NPREDS,
				//	*eNEEE	=erows[1]+kc3+3*4*WG_NPREDS,
					*eW	=erows[0]+kc3-1*4*WG_NPREDS,
					*ecurr	=erows[0]+kc3+0*4*WG_NPREDS;
				pred=wg_predict(wg_weights[kc], rows, 4*2, kc2, 0, wg_perrors[kc], eNW, eN, eNE, eNNE, wg_preds);
				int curr=src->data[idx+kc];
				if(g_dist>1)
				{
					if(fwd)
					{
						curr-=(int)pred;
						curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));//curr/=g_dist
						src->data[idx+kc]=curr;

						curr=g_dist*curr+(int)pred;
						CLAMP2(curr, amin[kc], amax[kc]);
					}
					else
					{
						curr=g_dist*curr+(int)pred;
						CLAMP2(curr, amin[kc], amax[kc]);

						src->data[idx+kc]=curr;
					}
				}
				else
				{
					if(fwd)
					{
						int error=curr-pred;
						error<<=32-src->depth[kc];
						error>>=32-src->depth[kc];
						src->data[idx+kc]=error;
					}
					else
					{
						curr+=pred;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						src->data[idx+kc]=curr;
					}
				}
				//int val;
				//if(fwd)
				//{
				//	val=(curr-(int)pred+g_dist/2)/g_dist;
				//	curr=g_dist*val+(int)pred;
				//}
				//else
				//{
				//	val=g_dist*curr+(int)pred;
				//	curr=val;
				//	CLAMP2(val, amin[kc], amax[kc]);
				//}
				//src->data[idx+kc]=keyboard[KEY_ALT]?pred:val;
				rows[0][kc2+0]=curr;
				rows[0][kc2+1]=curr-pred;
				//{
				//	int curr=src->data[idx+kc], pred0=pred;
				//	pred+=1<<WG_NBITS>>1;
				//	pred>>=WG_NBITS;
				//	int p2=pred;
				//	pred^=fwdmask;
				//	pred-=fwdmask;
				//	pred+=curr;
				//
				//	pred<<=32-src->depth[kc];
				//	pred>>=32-src->depth[kc];
				//
				//	src->data[idx+kc]=keyboard[KEY_ALT]?p2:pred;
				//	rows[0][kc2+0]=(fwd?curr:pred)<<WG_NBITS;
				//	rows[0][kc2+1]=rows[0][kc2]-pred0;
				//}
				wg_update(rows[0][kc2], kc, wg_preds, wg_perrors[kc], eW, ecurr, eNE);
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
			erows[0]+=4*WG_NPREDS;
			erows[1]+=4*WG_NPREDS;
			erows[2]+=4*WG_NPREDS;
			erows[3]+=4*WG_NPREDS;
#endif
		}
	}
	free(pixels);
	_mm_free(ebuf);
}