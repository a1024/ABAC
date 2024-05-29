#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT
	#define DISABLE_WGRAD


#include"ac.h"
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
#define BLOCKSIZE 256
static int cgrad(int N, int W, int NW)
{
	int pred;

	MEDIAN3_32(pred, N, W, N+W-NW);
	return pred;
}
static int clampav(int NW, int N, int NE, int WW, int W)
{
	ALIGN(16) int pred[4];
	__m128i va=_mm_set_epi32(0, 0, 0, N);
	__m128i vb=_mm_set_epi32(0, 0, 0, W);
	__m128i vc=_mm_set_epi32(0, 0, 0, NE);
	__m128i vd=_mm_set_epi32(0, 0, 0, (8*W+5*(N-NW)+NE-WW)>>3);
	__m128i vmin=_mm_min_epi32(va, vb);
	__m128i vmax=_mm_max_epi32(va, vb);
	vmin=_mm_min_epi32(vmin, vc);
	vmax=_mm_max_epi32(vmax, vc);
	vd=_mm_max_epi32(vd, vmin);
	vd=_mm_min_epi32(vd, vmax);
	_mm_store_si128((__m128i*)pred, vd);
	return pred[0];
}
static int clamp(int vmin, int x, int vmax)
{
	int ret;

	MEDIAN3_32(ret, vmin, vmax, x);
	return ret;
}
#define MAXPREDS 20
#define PREDLIST0\
	PRED( 0, curr[0], (N[0]+W[0])>>1)\
	PRED( 0, curr[0], clampav(NW[0], N[0], NE[0], WW[0], W[0]))\
	PRED( 1, curr[0], cgrad(N[0], W[0], NW[0]))
#define PREDLIST1\
	PRED( 2, curr[1], (N[1]+W[1])>>1)\
	PRED( 0, curr[1], clampav(NW[1], N[1], NE[1], WW[1], W[1]))\
	PRED( 3, curr[1], cgrad(N[1], W[1], NW[1]))\
	PRED( 4, curr[1], clamp(-half, ((N[1]-N[0]+W[1]-W[0])>>1)+curr[0], half-1))\
	PRED( 0, curr[1], clamp(-half, clampav(NW[1]-NW[0], N[1]-N[0], NE[1]-NE[0], WW[1]-WW[0], W[1]-W[0]), half-1)+curr[0])\
	PRED( 5, curr[1], clamp(-half, cgrad(N[1]-N[0], W[1]-W[0], NW[1]-NW[0])+curr[0], half-1))\
	PRED( 6, curr[0], clamp(-half, ((N[0]-N[1]+W[0]-W[1])>>1)+curr[1], half-1))\
	PRED( 0, curr[0], clamp(-half, clampav(NW[0]-NW[1], N[0]-N[1], NE[0]-NE[1], WW[0]-WW[1], W[0]-W[1]), half-1)+curr[1])\
	PRED( 7, curr[0], clamp(-half, cgrad(N[0]-N[1], W[0]-W[1], NW[0]-NW[1])+curr[1], half-1))
#define PREDLIST2\
	PRED( 8, curr[2], (N[2]+W[2])>>1)\
	PRED( 0, curr[2], clampav(NW[2], N[2], NE[2], WW[2], W[2]))\
	PRED( 9, curr[2], cgrad(N[2], W[2], NW[2]))\
	PRED(10, curr[2], clamp(-half, ((N[2]-N[0]+W[2]-W[0])>>1)+curr[0], half-1))\
	PRED( 0, curr[2], clamp(-half, clampav(NW[2]-NW[0], N[2]-N[0], NE[2]-NE[0], WW[2]-WW[0], W[2]-W[0]), half-1)+curr[0])\
	PRED(11, curr[2], clamp(-half, cgrad(N[2]-N[0], W[2]-W[0], NW[2]-NW[0])+curr[0], half-1))\
	PRED(12, curr[2], clamp(-half, ((N[2]-N[1]+W[2]-W[1])>>1)+curr[1], half-1))\
	PRED( 0, curr[2], clamp(-half, clampav(NW[2]-NW[1], N[2]-N[1], NE[2]-NE[1], WW[2]-WW[1], W[2]-W[1]), half-1)+curr[1])\
	PRED(13, curr[2], clamp(-half, cgrad(N[2]-N[1], W[2]-W[1], NW[2]-NW[1])+curr[1], half-1))\
	PRED(14, curr[1], clamp(-half, ((N[1]-N[2]+W[1]-W[2])>>1)+curr[2], half-1))\
	PRED( 0, curr[1], clamp(-half, clampav(NW[1]-NW[2], N[1]-N[2], NE[1]-NE[2], WW[1]-WW[2], W[1]-W[2]), half-1)+curr[2])\
	PRED(15, curr[1], clamp(-half, cgrad(N[1]-N[2], W[1]-W[2], NW[1]-NW[2])+curr[2], half-1))\
	PRED(16, curr[0], clamp(-half, ((N[0]-N[2]+W[0]-W[2])>>1)+curr[2], half-1))\
	PRED( 0, curr[0], clamp(-half, clampav(NW[0]-NW[2], N[0]-N[2], NE[0]-NE[2], WW[0]-WW[2], W[0]-W[2]), half-1)+curr[2])\
	PRED(17, curr[0], clamp(-half, cgrad(N[0]-N[2], W[0]-W[2], NW[0]-NW[2])+curr[2], half-1))
#define PREDLIST3\
	PRED(18, curr[3], (N[3]+W[3])>>1)\
	PRED( 0, curr[3], clampav(NW[3], N[3], NE[3], WW[3], W[3]))\
	PRED(19, curr[3], cgrad(N[3], W[3], NW[3]))

static const int combinations_3c[]=
{
	0, 2, 5,
	0, 2, 6,
	0, 2, 4,
	0, 1, 5,
	0, 1, 6,
	0, 1, 4,
	0, 5, 7,//021
//	0, 6, 7,
	0, 4, 7,//021

//	3, 2, 5,
//	3, 2, 6,
//	3, 2, 4,
	1, 3, 5,//102
	1, 3, 6,//102
	1, 3, 4,//102
//	5, 7, 3,
//	6, 7, 3,
	4, 7, 3,//210

//	5, 8, 2,
//	6, 8, 2,
	4, 8, 2,//201
//	1, 5, 8,
	1, 6, 8,//120
	1, 4, 8,//120
//	5, 8, 7,
//	6, 8, 7,
	4, 8, 7,//201
};
