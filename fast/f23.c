#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
//#ifdef _MSC_VER
//#include<intrin.h>
//#else
//#include<x86intrin.h>
//#endif
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT

//	#define UNROLL_DECODER//bad, binary incompatible!
//	#define ENABLE_WGRAD//bad


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
#ifdef ENABLE_WGRAD
static int wgrad(int N, int W, int X, int Y)//(X*N+Y*W)/(X+Y)
{
#ifdef _MSC_VER
	double pred=((double)X*N+(double)Y*W)/((double)X+Y);
	return _cvt_dtoi_fast(pred);
#else
	__m128d mres;
	__m128i mN=_mm_set_epi32(0, 0, 0, N);
	__m128i mW=_mm_set_epi32(0, 0, 0, W);
	__m128i mX=_mm_set_epi32(0, 0, 0, X);
	__m128i mY=_mm_set_epi32(0, 0, 0, Y);
	mN=_mm_mul_epi32(mX, mN);
	mW=_mm_mul_epi32(mY, mW);
	mN=_mm_add_epi32(mN, mW);
	mX=_mm_add_epi32(mX, mY);
	mres=_mm_div_pd(_mm_cvtepi32_pd(mN), _mm_cvtepi32_pd(mX));
	{
		ALIGN(16) int result[4];
		__m128i res2=_mm_cvtpd_epi32(mres);
		_mm_store_si128((__m128i*)result, res2);
		return result[0];
	}
#endif
}
#endif
static int clamp(int vmin, int x, int vmax)
{
	int ret;

	MEDIAN3_32(ret, vmin, vmax, x);
	return ret;
}
#define MAXPREDS 32
#if 1
#define PREDLIST0\
	PRED( 0, 0x00, curr[0], W[0])\
	PRED( 1, 0x00, curr[0], cgrad(N[0], W[0], NW[0]))
#define PREDLIST1\
	PRED( 2, 0x11, curr[1], W[1])\
	PRED( 3, 0x11, curr[1], cgrad(N[1], W[1], NW[1]))\
	PRED( 4, 0x10, curr[1], clamp(-half, W[1]-W[0]+curr[0], half-1))\
	PRED( 5, 0x10, curr[1], clamp(-half, cgrad(N[1]-N[0], W[1]-W[0], NW[1]-NW[0])+curr[0], half-1))\
	PRED( 6, 0x01, curr[0], clamp(-half, W[0]-W[1]+curr[1], half-1))\
	PRED( 7, 0x01, curr[0], clamp(-half, cgrad(N[0]-N[1], W[0]-W[1], NW[0]-NW[1])+curr[1], half-1))
#define PREDLIST2\
	PRED( 8, 0x22, curr[2], W[2])\
	PRED( 9, 0x22, curr[2], cgrad(N[2], W[2], NW[2]))\
	PRED(10, 0x20, curr[2], clamp(-half, W[2]-W[0]+curr[0], half-1))\
	PRED(11, 0x20, curr[2], clamp(-half, cgrad(N[2]-N[0], W[2]-W[0], NW[2]-NW[0])+curr[0], half-1))\
	PRED(12, 0x21, curr[2], clamp(-half, W[2]-W[1]+curr[1], half-1))\
	PRED(13, 0x21, curr[2], clamp(-half, cgrad(N[2]-N[1], W[2]-W[1], NW[2]-NW[1])+curr[1], half-1))\
	PRED(14, 0x12, curr[1], clamp(-half, W[1]-W[2]+curr[2], half-1))\
	PRED(15, 0x12, curr[1], clamp(-half, cgrad(N[1]-N[2], W[1]-W[2], NW[1]-NW[2])+curr[2], half-1))\
	PRED(16, 0x02, curr[0], clamp(-half, W[0]-W[2]+curr[2], half-1))\
	PRED(17, 0x02, curr[0], clamp(-half, cgrad(N[0]-N[2], W[0]-W[2], NW[0]-NW[2])+curr[2], half-1))
#define PREDLIST3\
	PRED(18, 0x33, curr[3], W[3])\
	PRED(19, 0x33, curr[3], cgrad(N[3], W[3], NW[3]))\
	PRED(20, 0x30, curr[3], clamp(-half, W[3]-W[0]+curr[0], half-1))\
	PRED(21, 0x30, curr[3], clamp(-half, cgrad(N[3]-N[0], W[3]-W[0], NW[3]-NW[0])+curr[0], half-1))\
	PRED(22, 0x31, curr[3], clamp(-half, W[3]-W[1]+curr[1], half-1))\
	PRED(23, 0x31, curr[3], clamp(-half, cgrad(N[3]-N[1], W[3]-W[1], NW[3]-NW[1])+curr[1], half-1))\
	PRED(24, 0x32, curr[3], clamp(-half, W[3]-W[2]+curr[2], half-1))\
	PRED(25, 0x32, curr[3], clamp(-half, cgrad(N[3]-N[2], W[3]-W[2], NW[3]-NW[2])+curr[2], half-1))\
	PRED(26, 0x23, curr[2], clamp(-half, W[2]-W[3]+curr[3], half-1))\
	PRED(27, 0x23, curr[2], clamp(-half, cgrad(N[2]-N[3], W[2]-W[3], NW[2]-NW[3])+curr[3], half-1))\
	PRED(28, 0x13, curr[1], clamp(-half, W[1]-W[3]+curr[3], half-1))\
	PRED(29, 0x13, curr[1], clamp(-half, cgrad(N[1]-N[3], W[1]-W[3], NW[1]-NW[3])+curr[3], half-1))\
	PRED(30, 0x03, curr[0], clamp(-half, W[0]-W[3]+curr[3], half-1))\
	PRED(31, 0x03, curr[0], clamp(-half, cgrad(N[0]-N[3], W[0]-W[3], NW[0]-NW[3])+curr[3], half-1))
#endif
#if 0
#define PREDLIST0\
	PRED( 0, 0x00, curr[0], ((N[0]+W[0])>>1))\
	PRED( 1, 0x00, curr[0], cgrad(N[0], W[0], NW[0]))
#define PREDLIST1\
	PRED( 2, 0x11, curr[1], ((N[1]+W[1])>>1))\
	PRED( 3, 0x11, curr[1], cgrad(N[1], W[1], NW[1]))\
	PRED( 4, 0x10, curr[1], clamp(-half, ((N[1]-N[0]+W[1]-W[0])>>1)+curr[0], half-1))\
	PRED( 5, 0x10, curr[1], clamp(-half, cgrad(N[1]-N[0], W[1]-W[0], NW[1]-NW[0])+curr[0], half-1))\
	PRED( 6, 0x01, curr[0], clamp(-half, ((N[0]-N[1]+W[0]-W[1])>>1)+curr[1], half-1))\
	PRED( 7, 0x01, curr[0], clamp(-half, cgrad(N[0]-N[1], W[0]-W[1], NW[0]-NW[1])+curr[1], half-1))
#define PREDLIST2\
	PRED( 8, 0x22, curr[2], ((N[2]+W[2])>>1))\
	PRED( 9, 0x22, curr[2], cgrad(N[2], W[2], NW[2]))\
	PRED(10, 0x20, curr[2], clamp(-half, ((N[2]-N[0]+W[2]-W[0])>>1)+curr[0], half-1))\
	PRED(11, 0x20, curr[2], clamp(-half, cgrad(N[2]-N[0], W[2]-W[0], NW[2]-NW[0])+curr[0], half-1))\
	PRED(12, 0x21, curr[2], clamp(-half, ((N[2]-N[1]+W[2]-W[1])>>1)+curr[1], half-1))\
	PRED(13, 0x21, curr[2], clamp(-half, cgrad(N[2]-N[1], W[2]-W[1], NW[2]-NW[1])+curr[1], half-1))\
	PRED(14, 0x12, curr[1], clamp(-half, ((N[1]-N[2]+W[1]-W[2])>>1)+curr[2], half-1))\
	PRED(15, 0x12, curr[1], clamp(-half, cgrad(N[1]-N[2], W[1]-W[2], NW[1]-NW[2])+curr[2], half-1))\
	PRED(16, 0x02, curr[0], clamp(-half, ((N[0]-N[2]+W[0]-W[2])>>1)+curr[2], half-1))\
	PRED(17, 0x02, curr[0], clamp(-half, cgrad(N[0]-N[2], W[0]-W[2], NW[0]-NW[2])+curr[2], half-1))
#define PREDLIST3\
	PRED(18, 0x33, curr[3], ((N[3]+W[3])>>1))\
	PRED(19, 0x33, curr[3], cgrad(N[3], W[3], NW[3]))\
	PRED(20, 0x30, curr[3], clamp(-half, ((N[3]-N[0]+W[3]-W[0])>>1)+curr[0], half-1))\
	PRED(21, 0x30, curr[3], clamp(-half, cgrad(N[3]-N[0], W[3]-W[0], NW[3]-NW[0])+curr[0], half-1))\
	PRED(22, 0x31, curr[3], clamp(-half, ((N[3]-N[1]+W[3]-W[1])>>1)+curr[1], half-1))\
	PRED(23, 0x31, curr[3], clamp(-half, cgrad(N[3]-N[1], W[3]-W[1], NW[3]-NW[1])+curr[1], half-1))\
	PRED(24, 0x32, curr[3], clamp(-half, ((N[3]-N[2]+W[3]-W[2])>>1)+curr[2], half-1))\
	PRED(25, 0x32, curr[3], clamp(-half, cgrad(N[3]-N[2], W[3]-W[2], NW[3]-NW[2])+curr[2], half-1))\
	PRED(26, 0x23, curr[2], clamp(-half, ((N[2]-N[3]+W[2]-W[3])>>1)+curr[3], half-1))\
	PRED(27, 0x23, curr[2], clamp(-half, cgrad(N[2]-N[3], W[2]-W[3], NW[2]-NW[3])+curr[3], half-1))\
	PRED(28, 0x13, curr[1], clamp(-half, ((N[1]-N[3]+W[1]-W[3])>>1)+curr[3], half-1))\
	PRED(29, 0x13, curr[1], clamp(-half, cgrad(N[1]-N[3], W[1]-W[3], NW[1]-NW[3])+curr[3], half-1))\
	PRED(30, 0x03, curr[0], clamp(-half, ((N[0]-N[3]+W[0]-W[3])>>1)+curr[3], half-1))\
	PRED(31, 0x03, curr[0], clamp(-half, cgrad(N[0]-N[3], W[0]-W[3], NW[0]-NW[3])+curr[3], half-1))
#endif
#if 0
#define PREDLIST0\
	PRED( 0, 0x00, curr[0], W[0])\
	PRED( 1, 0x00, curr[0], N[0])
#define PREDLIST1\
	PRED( 2, 0x11, curr[1], W[1])\
	PRED( 3, 0x11, curr[1], N[1])\
	PRED( 4, 0x10, curr[1], clamp(-half, W[1]-W[0]+curr[0], half-1))\
	PRED( 5, 0x10, curr[1], clamp(-half, N[1]-N[0]+curr[0], half-1))\
	PRED( 6, 0x01, curr[0], clamp(-half, W[0]-W[1]+curr[1], half-1))\
	PRED( 7, 0x01, curr[0], clamp(-half, N[0]-N[1]+curr[1], half-1))
#define PREDLIST2\
	PRED( 8, 0x22, curr[2], W[2])\
	PRED( 9, 0x22, curr[2], N[2])\
	PRED(10, 0x20, curr[2], clamp(-half, W[2]-W[0]+curr[0], half-1))\
	PRED(11, 0x20, curr[2], clamp(-half, N[2]-N[0]+curr[0], half-1))\
	PRED(12, 0x21, curr[2], clamp(-half, W[2]-W[1]+curr[1], half-1))\
	PRED(13, 0x21, curr[2], clamp(-half, N[2]-N[1]+curr[1], half-1))\
	PRED(14, 0x12, curr[1], clamp(-half, W[1]-W[2]+curr[2], half-1))\
	PRED(15, 0x12, curr[1], clamp(-half, N[1]-N[2]+curr[2], half-1))\
	PRED(16, 0x02, curr[0], clamp(-half, W[0]-W[2]+curr[2], half-1))\
	PRED(17, 0x02, curr[0], clamp(-half, N[0]-N[2]+curr[2], half-1))
#define PREDLIST3\
	PRED(18, 0x33, curr[3], W[3])\
	PRED(19, 0x33, curr[3], N[3])\
	PRED(20, 0x30, curr[3], clamp(-half, W[3]-W[0]+curr[0], half-1))\
	PRED(21, 0x30, curr[3], clamp(-half, N[3]-N[0]+curr[0], half-1))\
	PRED(22, 0x31, curr[3], clamp(-half, W[3]-W[1]+curr[1], half-1))\
	PRED(23, 0x31, curr[3], clamp(-half, N[3]-N[1]+curr[1], half-1))\
	PRED(24, 0x32, curr[3], clamp(-half, W[3]-W[2]+curr[2], half-1))\
	PRED(25, 0x32, curr[3], clamp(-half, N[3]-N[2]+curr[2], half-1))\
	PRED(26, 0x23, curr[2], clamp(-half, W[2]-W[3]+curr[3], half-1))\
	PRED(27, 0x23, curr[2], clamp(-half, N[2]-N[3]+curr[3], half-1))\
	PRED(28, 0x13, curr[1], clamp(-half, W[1]-W[3]+curr[3], half-1))\
	PRED(29, 0x13, curr[1], clamp(-half, N[1]-N[3]+curr[3], half-1))\
	PRED(30, 0x03, curr[0], clamp(-half, W[0]-W[3]+curr[3], half-1))\
	PRED(31, 0x03, curr[0], clamp(-half, N[0]-N[3]+curr[3], half-1))
#endif

#if 0
#define PREDLIST0\
	PRED( 0, 0x00, curr[0], wgrad(N[0], W[0], X[0], Y[0]))\
	PRED( 1, 0x00, curr[0], cgrad(N[0], W[0], NW[0]))
#define PREDLIST1\
	PRED( 2, 0x11, curr[1], wgrad(N[1], W[1], X[1], Y[1]))\
	PRED( 3, 0x11, curr[1], cgrad(N[1], W[1], NW[1]))\
	PRED( 4, 0x10, curr[1], clamp(-half, wgrad(N[1]-N[0], W[1]-W[0], X[1]+X[0], Y[1]+Y[0])+curr[0], half-1))\
	PRED( 5, 0x10, curr[1], clamp(-half, cgrad(N[1]-N[0], W[1]-W[0], NW[1]-NW[0])+curr[0], half-1))\
	PRED( 6, 0x01, curr[0], clamp(-half, wgrad(N[0]-N[1], W[0]-W[1], X[0]+X[1], Y[0]+Y[1])+curr[1], half-1))\
	PRED( 7, 0x01, curr[0], clamp(-half, cgrad(N[0]-N[1], W[0]-W[1], NW[0]-NW[1])+curr[1], half-1))
#define PREDLIST2\
	PRED( 8, 0x22, curr[2], wgrad(N[2], W[2], X[2], Y[2]))\
	PRED( 9, 0x22, curr[2], cgrad(N[2], W[2], NW[2]))\
	PRED(10, 0x20, curr[2], clamp(-half, wgrad(N[2]-N[0], W[2]-W[0], X[2]+X[0], Y[2]+Y[0])+curr[0], half-1))\
	PRED(11, 0x20, curr[2], clamp(-half, cgrad(N[2]-N[0], W[2]-W[0], NW[2]-NW[0])+curr[0], half-1))\
	PRED(12, 0x21, curr[2], clamp(-half, wgrad(N[2]-N[1], W[2]-W[1], X[2]+X[1], Y[2]+Y[1])+curr[1], half-1))\
	PRED(13, 0x21, curr[2], clamp(-half, cgrad(N[2]-N[1], W[2]-W[1], NW[2]-NW[1])+curr[1], half-1))\
	PRED(14, 0x12, curr[1], clamp(-half, wgrad(N[1]-N[2], W[1]-W[2], X[1]+X[2], Y[1]+Y[2])+curr[2], half-1))\
	PRED(15, 0x12, curr[1], clamp(-half, cgrad(N[1]-N[2], W[1]-W[2], NW[1]-NW[2])+curr[2], half-1))\
	PRED(16, 0x02, curr[0], clamp(-half, wgrad(N[0]-N[2], W[0]-W[2], X[0]+X[2], Y[0]+Y[2])+curr[2], half-1))\
	PRED(17, 0x02, curr[0], clamp(-half, cgrad(N[0]-N[2], W[0]-W[2], NW[0]-NW[2])+curr[2], half-1))
#define PREDLIST3\
	PRED(18, 0x33, curr[3], wgrad(N[3], W[3], X[3], Y[3]))\
	PRED(19, 0x33, curr[3], cgrad(N[3], W[3], NW[3]))\
	PRED(20, 0x30, curr[3], clamp(-half, wgrad(N[3]-N[0], W[3]-W[0], X[3]+X[0], Y[3]+Y[0])+curr[0], half-1))\
	PRED(21, 0x30, curr[3], clamp(-half, cgrad(N[3]-N[0], W[3]-W[0], NW[3]-NW[0])+curr[0], half-1))\
	PRED(22, 0x31, curr[3], clamp(-half, wgrad(N[3]-N[1], W[3]-W[1], X[3]+X[1], Y[3]+Y[1])+curr[1], half-1))\
	PRED(23, 0x31, curr[3], clamp(-half, cgrad(N[3]-N[1], W[3]-W[1], NW[3]-NW[1])+curr[1], half-1))\
	PRED(24, 0x32, curr[3], clamp(-half, wgrad(N[3]-N[2], W[3]-W[2], X[3]+X[2], Y[3]+Y[2])+curr[2], half-1))\
	PRED(25, 0x32, curr[3], clamp(-half, cgrad(N[3]-N[2], W[3]-W[2], NW[3]-NW[2])+curr[2], half-1))\
	PRED(26, 0x23, curr[2], clamp(-half, wgrad(N[2]-N[3], W[2]-W[3], X[2]+X[3], Y[2]+Y[3])+curr[3], half-1))\
	PRED(27, 0x23, curr[2], clamp(-half, cgrad(N[2]-N[3], W[2]-W[3], NW[2]-NW[3])+curr[3], half-1))\
	PRED(28, 0x13, curr[1], clamp(-half, wgrad(N[1]-N[3], W[1]-W[3], X[1]+X[3], Y[1]+Y[3])+curr[3], half-1))\
	PRED(29, 0x13, curr[1], clamp(-half, cgrad(N[1]-N[3], W[1]-W[3], NW[1]-NW[3])+curr[3], half-1))\
	PRED(30, 0x03, curr[0], clamp(-half, wgrad(N[0]-N[3], W[0]-W[3], X[0]+X[3], Y[0]+Y[3])+curr[3], half-1))\
	PRED(31, 0x03, curr[0], clamp(-half, cgrad(N[0]-N[3], W[0]-W[3], NW[0]-NW[3])+curr[3], half-1))
#endif
#if 0
#define PREDLIST0\
	PRED( 0, 0x00, curr[0], clampav(NW[0], N[0], NE[0], WW[0], W[0]))\
	PRED( 1, 0x00, curr[0], cgrad(N[0], W[0], NW[0]))
#define PREDLIST1\
	PRED( 2, 0x11, curr[1], clampav(NW[1], N[1], NE[1], WW[1], W[1]))\
	PRED( 3, 0x11, curr[1], cgrad(N[1], W[1], NW[1]))\
	PRED( 4, 0x10, curr[1], clamp(-half, clampav(NW[1]-NW[0], N[1]-N[0], NE[1]-NE[0], WW[1]-WW[0], W[1]-W[0])+curr[0], half-1))\
	PRED( 5, 0x10, curr[1], clamp(-half, cgrad(N[1]-N[0], W[1]-W[0], NW[1]-NW[0])+curr[0], half-1))\
	PRED( 6, 0x01, curr[0], clamp(-half, clampav(NW[0]-NW[1], N[0]-N[1], NE[0]-NE[1], WW[0]-WW[1], W[0]-W[1])+curr[1], half-1))\
	PRED( 7, 0x01, curr[0], clamp(-half, cgrad(N[0]-N[1], W[0]-W[1], NW[0]-NW[1])+curr[1], half-1))
#define PREDLIST2\
	PRED( 8, 0x22, curr[2], clampav(NW[2], N[2], NE[2], WW[2], W[2]))\
	PRED( 9, 0x22, curr[2], cgrad(N[2], W[2], NW[2]))\
	PRED(10, 0x20, curr[2], clamp(-half, clampav(NW[2]-NW[0], N[2]-N[0], NE[2]-NE[0], WW[2]-WW[0], W[2]-W[0])+curr[0], half-1))\
	PRED(11, 0x20, curr[2], clamp(-half, cgrad(N[2]-N[0], W[2]-W[0], NW[2]-NW[0])+curr[0], half-1))\
	PRED(12, 0x21, curr[2], clamp(-half, clampav(NW[2]-NW[1], N[2]-N[1], NE[2]-NE[1], WW[2]-WW[1], W[2]-W[1])+curr[1], half-1))\
	PRED(13, 0x21, curr[2], clamp(-half, cgrad(N[2]-N[1], W[2]-W[1], NW[2]-NW[1])+curr[1], half-1))\
	PRED(14, 0x12, curr[1], clamp(-half, clampav(NW[1]-NW[2], N[1]-N[2], NE[1]-NE[2], WW[1]-WW[2], W[1]-W[2])+curr[2], half-1))\
	PRED(15, 0x12, curr[1], clamp(-half, cgrad(N[1]-N[2], W[1]-W[2], NW[1]-NW[2])+curr[2], half-1))\
	PRED(16, 0x02, curr[0], clamp(-half, clampav(NW[0]-NW[2], N[0]-N[2], NE[0]-NE[2], WW[0]-WW[2], W[0]-W[2])+curr[2], half-1))\
	PRED(17, 0x02, curr[0], clamp(-half, cgrad(N[0]-N[2], W[0]-W[2], NW[0]-NW[2])+curr[2], half-1))
#define PREDLIST3\
	PRED(18, 0x33, curr[3], clampav(NW[3], N[3], NE[3], WW[3], W[3]))\
	PRED(19, 0x33, curr[3], cgrad(N[3], W[3], NW[3]))\
	PRED(20, 0x30, curr[3], clamp(-half, clampav(NW[3]-NW[0], N[3]-N[0], NE[3]-NE[0], WW[3]-WW[0], W[3]-W[0])+curr[0], half-1))\
	PRED(21, 0x30, curr[3], clamp(-half, cgrad(N[3]-N[0], W[3]-W[0], NW[3]-NW[0])+curr[0], half-1))\
	PRED(22, 0x31, curr[3], clamp(-half, clampav(NW[3]-NW[1], N[3]-N[1], NE[3]-NE[1], WW[3]-WW[1], W[3]-W[1])+curr[1], half-1))\
	PRED(23, 0x31, curr[3], clamp(-half, cgrad(N[3]-N[1], W[3]-W[1], NW[3]-NW[1])+curr[1], half-1))\
	PRED(24, 0x32, curr[3], clamp(-half, clampav(NW[3]-NW[2], N[3]-N[2], NE[3]-NE[2], WW[3]-WW[2], W[3]-W[2])+curr[2], half-1))\
	PRED(25, 0x32, curr[3], clamp(-half, cgrad(N[3]-N[2], W[3]-W[2], NW[3]-NW[2])+curr[2], half-1))\
	PRED(26, 0x23, curr[2], clamp(-half, clampav(NW[2]-NW[3], N[2]-N[3], NE[2]-NE[3], WW[2]-WW[3], W[2]-W[3])+curr[3], half-1))\
	PRED(27, 0x23, curr[2], clamp(-half, cgrad(N[2]-N[3], W[2]-W[3], NW[2]-NW[3])+curr[3], half-1))\
	PRED(28, 0x13, curr[1], clamp(-half, clampav(NW[1]-NW[3], N[1]-N[3], NE[1]-NE[3], WW[1]-WW[3], W[1]-W[3])+curr[3], half-1))\
	PRED(29, 0x13, curr[1], clamp(-half, cgrad(N[1]-N[3], W[1]-W[3], NW[1]-NW[3])+curr[3], half-1))\
	PRED(30, 0x03, curr[0], clamp(-half, clampav(NW[0]-NW[3], N[0]-N[3], NE[0]-NE[3], WW[0]-WW[3], W[0]-W[3])+curr[3], half-1))\
	PRED(31, 0x03, curr[0], clamp(-half, cgrad(N[0]-N[3], W[0]-W[3], NW[0]-NW[3])+curr[3], half-1))
#endif
typedef struct _ThreadArgs
{
	const Image *src;
	Image *dst;
	int fwd, loud, y1, y2;
	int bufsize, histsize;
	int *pixels, *hist;

	DList list;
	const unsigned char *decstart, *decend;
} ThreadArgs;

static const char *prednames[]=
{
#define PRED(IDX, DFLAG, VAL, EXPR) #VAL " - " #EXPR,
	PREDLIST0
	PREDLIST1
	PREDLIST2
	PREDLIST3
#undef  PRED
};

//indices are permuted for causality, cyclic dependencies are commented out
static const int combinations_4c[]=
{
	//0<-0...
	 0,  4, 10, 20,
	 0,  4, 10, 22,
	 0,  4, 10, 24,
	 0,  4, 10, 18,
	 0,  4, 12, 20,
	 0,  4, 12, 22,
	 0,  4, 12, 24,
	 0,  4, 12, 18,
	 0,  4,  8, 20,
	 0,  4,  8, 22,
	 0,  4,  8, 24,
	 0,  4,  8, 18,
	 0,  4, 20, 26,
	 0,  4, 22, 26,
//	 0,  4, 24, 26,
	 0,  4, 18, 26,

	 0,  2, 10, 20,
	 0,  2, 10, 22,
	 0,  2, 10, 24,
	 0,  2, 10, 18,
	 0,  2, 12, 20,
	 0,  2, 12, 22,
	 0,  2, 12, 24,
	 0,  2, 12, 18,
	 0,  2,  8, 20,
	 0,  2,  8, 22,
	 0,  2,  8, 24,
	 0,  2,  8, 18,
	 0,  2, 20, 26,
	 0,  2, 22, 26,
//	 0,  2, 24, 26,
	 0,  2, 18, 26,

	 0, 10, 14, 20,
	 0, 10, 14, 22,
	 0, 10, 14, 24,
	 0, 10, 14, 18,
//	 0, 12, 14, 20,
//	 0, 12, 14, 22,
//	 0, 12, 14, 24,
//	 0, 12, 14, 18,
	 0,  8, 14, 20,
	 0,  8, 14, 22,
	 0,  8, 14, 24,
	 0,  8, 14, 18,
	 0, 20, 26, 14,
//	 0, 22, 26, 14,
//	 0, 24, 26, 14,
	 0, 18, 26, 14,

	 0, 10, 20, 28,
	 0, 10, 22, 28,
	 0, 10, 24, 28,
	 0, 10, 18, 28,
	 0, 20, 28, 12,
//	 0, 22, 28, 12,
//	 0, 24, 28, 12,
	 0, 18, 28, 12,
	 0,  8, 20, 28,
//	 0,  8, 22, 28,
	 0,  8, 24, 28,
	 0,  8, 18, 28,
	 0, 20, 28, 26,
//	 0, 22, 28, 26,
//	 0, 24, 28, 26,
	 0, 18, 28, 26,

	 
	//0<-1...
//	 6,  4, 10, 20,
//	 6,  4, 10, 22,
//	 6,  4, 10, 24,
//	 6,  4, 10, 18,
//	 6,  4, 12, 20,
//	 6,  4, 12, 22,
//	 6,  4, 12, 24,
//	 6,  4, 12, 18,
//	 6,  4,  8, 20,
//	 6,  4,  8, 22,
//	 6,  4,  8, 24,
//	 6,  4,  8, 18,
//	 6,  4, 20, 26,
//	 6,  4, 22, 26,
//	 6,  4, 24, 26,
//	 6,  4, 18, 26,

	 2,  6, 10, 20,
	 2,  6, 10, 22,
	 2,  6, 10, 24,
	 2,  6, 10, 18,
	 2,  6, 12, 20,
	 2,  6, 12, 22,
	 2,  6, 12, 24,
	 2,  6, 12, 18,
	 2,  6,  8, 20,
	 2,  6,  8, 22,
	 2,  6,  8, 24,
	 2,  6,  8, 18,
	 2,  6, 26, 20,
	 2,  6, 26, 22,
//	 2,  6, 26, 24,
	 2,  6, 26, 18,

//	 6, 10, 14, 20,
//	 6, 10, 14, 22,
//	 6, 10, 14, 24,
//	 6, 10, 14, 18,
//	 6, 12, 14, 20,
//	 6, 12, 14, 22,
//	 6, 12, 14, 24,
//	 6, 12, 14, 18,
	 6,  8, 14, 20,
	 6,  8, 14, 22,
	 6,  8, 14, 24,
	 6,  8, 14, 18,
	 6, 20, 26, 14,
//	 6, 22, 26, 14,
//	 6, 24, 26, 14,
	 6, 18, 26, 14,

	 6, 10, 20, 28,
	 6, 10, 22, 28,
	 6, 10, 24, 28,
	 6, 10, 18, 28,
	 6, 20, 28, 12,
//	 6, 22, 28, 12,
//	 6, 24, 28, 12,
	 6, 18, 28, 12,
	 6,  8, 20, 28,
//	 6,  8, 22, 28,
	 6,  8, 24, 28,
	 6,  8, 18, 28,
	 6, 20, 28, 26,
//	 6, 22, 28, 26,
//	 6, 24, 28, 26,
	 6, 18, 28, 26,

	 
	//0<-2...
//	16,  4, 10, 20,
//	16,  4, 10, 22,
//	16,  4, 10, 24,
//	16,  4, 10, 18,
//	16,  4, 12, 20,
//	16,  4, 12, 22,
//	16,  4, 12, 24,
//	16,  4, 12, 18,
	 8, 16,  4, 20,
	 8, 16,  4, 22,
	 8, 16,  4, 24,
	 8, 16,  4, 18,
//	16,  4, 20, 26,
//	16,  4, 22, 26,
//	16,  4, 24, 26,
//	16,  4, 18, 26,

//	16,  2, 10, 20,
//	16,  2, 10, 22,
//	16,  2, 10, 24,
//	16,  2, 10, 18,
	 2, 12, 16, 20,
	 2, 12, 16, 22,
	 2, 12, 16, 24,
	 2, 12, 16, 18,
	 2,  8, 16, 20,
	 2,  8, 16, 22,
	 2,  8, 16, 24,
	 2,  8, 16, 18,
//	16,  2, 20, 26,
//	16,  2, 22, 26,
//	16,  2, 24, 26,
//	16,  2, 18, 26,

//	16, 14, 10, 20,
//	16, 14, 10, 22,
//	16, 14, 10, 24,
//	16, 14, 10, 18,
//	16, 14, 12, 20,
//	16, 14, 12, 22,
//	16, 14, 12, 24,
//	16, 14, 12, 18,
	 8, 16, 14, 20,
	 8, 16, 14, 22,
	 8, 16, 14, 24,
	 8, 16, 14, 18,
//	16, 14, 26, 20,
//	16, 14, 26, 22,
//	16, 14, 26, 24,
	18, 26, 16, 14,//3201

//	16, 28, 10, 20,
//	16, 28, 10, 22,
//	16, 28, 10, 24,
//	16, 28, 10, 18,
//	16, 28, 12, 20,
//	16, 28, 12, 22,
//	16, 28, 12, 24,
	18, 28, 12, 16,//3120
	 8, 16, 20, 28,//2031
//	 8, 16, 22, 28,
	 8, 16, 24, 28,//2031
	 8, 16, 18, 28,//2031
//	16, 28, 26, 20,
//	16, 28, 26, 22,
//	16, 28, 26, 24,
	18, 28, 26, 16,//3120

	 
	//0<-3...
//	20, 30,  4, 10,
//	22, 30,  4, 10,
//	24, 30,  4, 10,
	18, 30,  4, 10,//3012
//	20, 30,  4, 12,
//	22, 30,  4, 12,
//	24, 30,  4, 12,
	18, 30,  4, 12,//3012
//	30,  4,  8, 20,
//	30,  4,  8, 22,
	 8, 24, 30,  4,//2301
	 8, 18, 30,  4,//2301
//	26, 30,  4, 20,
//	26, 30,  4, 22,
//	26, 30,  4, 24,
	26, 30,  4, 18,//3012

//	30,  2, 10, 20,
	 2, 22, 30, 10,//1302
//	30,  2, 10, 24,
//	30,  2, 10, 18,
//	 2, 12, 20, 30,
	 2, 12, 22, 30,//1230
	 2, 12, 24, 30,//1230
	 2, 12, 18, 30,//1230
//	 2,  8, 20, 30,
	 2,  8, 22, 30,//1230
	 2,  8, 24, 30,//1230
	 2,  8, 18, 30,//1230
//	 2, 26, 30, 20,
	 2, 26, 30, 22,//1302
//	 2, 26, 30, 24,
	 2, 26, 30, 18,//1302

//	30, 14, 10, 20,
//	30, 14, 10, 22,
//	30, 14, 10, 24,
	18, 30, 10, 14,//3021
//	30, 14, 12, 20,
//	30, 14, 12, 22,
//	30, 14, 12, 24,
//	30, 14, 12, 18,
//	30, 14,  8, 20,
	 8, 14, 22, 30,//2130
	 8, 14, 24, 30,//2130
	 8, 14, 18, 30,//2130
//	30, 14, 26, 20,
//	30, 14, 26, 22,
//	30, 14, 26, 24,
	18, 30, 26, 14,//3021

//	20, 30, 28, 10,
//	22, 30, 28, 10,
//	24, 30, 28, 10,
	18, 30, 28, 10,//3012
//	20, 30, 28, 12,
//	22, 30, 28, 12,
//	24, 30, 28, 12,
	18, 30, 28, 12,//3012
//	30, 28,  8, 20,
//	30, 28,  8, 22,
	 8, 24, 30, 28,//2301
	30, 28,  8, 18,
	20, 30, 28, 26,//3012
//	22, 30, 28, 26,
//	24, 30, 28, 26,
	18, 30, 28, 26,//3012
};
static const int combinations_3c[]=
{
	 0,  4, 10,
	 0,  4, 12,
	 0,  4,  8,
	 0,  2, 10,
	 0,  2, 12,
	 0,  2,  8,
	 0, 10, 14,//021
//	 0, 12, 14,
	 0,  8, 14,//021

//	 6,  4, 10,
//	 6,  4, 12,
//	 6,  4,  8,
	 2,  6, 10,//102
	 2,  6, 12,//102
	 2,  6,  8,//102
//	10, 14, 6,
//	12, 14, 6,
	 8, 14, 6,//210

//	10, 16,  4,
//	12, 16,  4,
	 8, 16,  4,//201
//	 2, 10, 16,
	 2, 12, 16,//120
	 2,  8, 16,//120
//	10, 16, 14,
//	12, 16, 14,
	 8, 16, 14,//201
};
static const int combinations_2c[]=
{
	0, 4,//01
	0, 2,//01
//	6, 4,
	2, 6,//10
};
#if 0
static void test4cycles()
{
	static const int dependency_flags[]=
	{
#define PRED(IDX, DFLAG, VAL, EXPR) DFLAG,
		PREDLIST0
		PREDLIST1
		PREDLIST2
		PREDLIST3
#undef  PRED
	};
	for(int kg=0;kg<_countof(combinations_4c);kg+=4)
	{
		int *group=combinations_4c+kg;
		for(int k=0;k<4;++k)
		{
			int kv=k;
			for(int it=0;;++it)
			{
				if(it>=4)
					printf("Cycle at group %d:  %2d, %2d, %2d, %2d", kg>>2, group[0], group[1], group[2], group[3]);
				int val=dependency_flags[kv];
				if((val&15)==(val>>4&15))
					break;
				kv=val&15;
			}
		}
	}
	for(int kg=0;kg<_countof(combinations_3c);kg+=3)
	{
		int *group=combinations_3c+kg;
		for(int k=0;k<3;++k)
		{
			int kv=k;
			for(int it=0;;++it)
			{
				if(it>=3)
					printf("Cycle at group %d:  %2d, %2d, %2d", kg>>2, group[0], group[1], group[2]);
				int val=dependency_flags[kv];
				if((val&15)==(val>>4&15))
					break;
				kv=val&15;
			}
		}
	}
	for(int kg=0;kg<_countof(combinations_2c);kg+=2)
	{
		int *group=combinations_2c+kg;
		for(int k=0;k<2;++k)
		{
			int kv=k;
			for(int it=0;;++it)
			{
				if(it>=2)
					printf("Cycle at group %d:  %2d, %2d", kg>>2, group[0], group[1]);
				int val=dependency_flags[kv];
				if((val&15)==(val>>4&15))
					break;
				kv=val&15;
			}
		}
	}
	printf("Done.\n");
	pause();
}
#endif
static void block_thread(void *param)
{
//	static const int combination_idx[]=
//	{//	 L   L   C   C   C   C   C   C
//		 0,  1,  6,  7, 16, 17, 30, 31,
//		 2,  3,  4,  5, 14, 15, 28, 29,
//		 8,  9, 10, 11, 12, 13, 26, 27,
//		18, 19, 20, 21, 22, 23, 24, 25,
//	};
//	static const int dependency_flags[]=
//	{
//#define PRED(IDX, DFLAG, VAL, EXPR) DFLAG,
//		PREDLIST0
//		PREDLIST1
//		PREDLIST2
//		PREDLIST3
//#undef  PRED
//	};

	GolombRiceCoder ec;
	ThreadArgs *args=(ThreadArgs*)param;
	memset(args->pixels, 0, args->bufsize);
	if(args->fwd)
	{
		Image const *image=args->src;
		int nlevels=1<<image->depth, half=nlevels>>1, mask=nlevels-1;
		int poolNch=0;
		double csizes[MAXPREDS]={0};
		char predsel[MAXPREDS>>1]={0};
		double best_csize=0;
		int comb_idx=0, combination[4]={0}, flag=0;
		//int best_lumas[4]={0}, best_chromas[4]={0}, combination[4]={0};
		//int permutation[4]={0};
		int res=image->iw*(args->y2-args->y1);
		switch(image->nch)
		{
		case 4:poolNch+=7*2;
		case 3:poolNch+=5*2;
		case 2:poolNch+=3*2;
		case 1:poolNch+=1*2;
			break;
		}
		memset(args->hist, 0, args->histsize);
		for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
		{
			ALIGN(16) int *rows[]=
			{
				args->pixels+(((image->iw+16LL)*((ky-0LL)&3)+8)<<3),
				args->pixels+(((image->iw+16LL)*((ky-1LL)&3)+8)<<3),
				args->pixels+(((image->iw+16LL)*((ky-2LL)&3)+8)<<3),
				args->pixels+(((image->iw+16LL)*((ky-3LL)&3)+8)<<3),
			};
#ifdef ENABLE_WGRAD
			ALIGN(16) int X[4]={0}, Y[4]={0};
#endif
			int preds[MAXPREDS]={0}, j;
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				int
#ifdef ENABLE_WGRAD
					*NN	=rows[2]+0*8,
					*NNE	=rows[2]+1*8,
					*NE	=rows[1]+1*8,
					*WW	=rows[0]-2*8,
#endif
					*NW	=rows[1]-1*8,
					*N	=rows[1]+0*8,
					*W	=rows[0]-1*8,
					*curr	=rows[0]+0*8;
				for(int kc=0;kc<image->nch;++kc)
					curr[kc]=image->data[idx+kc];
#ifdef ENABLE_WGRAD
				{
					__m128i mX=_mm_set1_epi32(1);
					__m128i mY=mX;

					__m128i mNN	=_mm_load_si128((__m128i*)NN);
					__m128i mNNE	=_mm_load_si128((__m128i*)NNE);
					__m128i mNW	=_mm_load_si128((__m128i*)NW);
					__m128i mN	=_mm_load_si128((__m128i*)N);
					__m128i mNE	=_mm_load_si128((__m128i*)NE);
					__m128i mWW	=_mm_load_si128((__m128i*)WW);
					__m128i mW	=_mm_load_si128((__m128i*)W);
					mX=_mm_add_epi32(mX, _mm_abs_epi32(_mm_sub_epi32(mW, mWW)));
					mY=_mm_add_epi32(mY, _mm_abs_epi32(_mm_sub_epi32(mW, mNW)));
					mX=_mm_add_epi32(mX, _mm_abs_epi32(_mm_sub_epi32(mN, mNW)));
					mY=_mm_add_epi32(mY, _mm_abs_epi32(_mm_sub_epi32(mN, mNN)));
					mX=_mm_add_epi32(mX, _mm_abs_epi32(_mm_sub_epi32(mNE, mN)));
					mY=_mm_add_epi32(mY, _mm_abs_epi32(_mm_sub_epi32(mNE, mNNE)));
					mX=_mm_add_epi32(mX, _mm_set1_epi32(15));
					mY=_mm_add_epi32(mY, _mm_set1_epi32(15));
					mX=_mm_srli_epi32(mX, 4);
					mY=_mm_srli_epi32(mY, 4);
					_mm_store_si128((__m128i*)X, mX);
					_mm_store_si128((__m128i*)Y, mY);
				}
#endif
#define PRED(IDX, DFLAG, VAL, EXPR) preds[j++]=(VAL)-(EXPR);
				j=0;
				PREDLIST0
				if(image->nch>=2)
				{
					PREDLIST1
					if(image->nch>=3)
					{
						PREDLIST2
						if(image->nch>=4)
						{
							PREDLIST3
						}
					}
				}
#undef  PRED
				for(int k=0;k<j;++k)
				{
					int val=preds[k];
					val+=half;
					val&=mask;
					++args->hist[k<<image->depth|val];
				}
				rows[0]+=8;
				rows[1]+=8;
#ifdef ENABLE_WGRAD
				rows[2]+=8;
				rows[3]+=8;
#endif
			}
		}
		for(int kc=0;kc<poolNch;++kc)
		{
			int *curr_hist=args->hist+((size_t)kc<<image->depth);
			for(int ks=0;ks<nlevels;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
					csizes[kc]-=freq*log2((double)freq/res);
			}
			csizes[kc]/=8;
		}

		//select best combination:
		{
			const int *groups=0, *group;
			int glen=0;

			for(int k=0;k<poolNch;k+=2)
				predsel[k>>1]=csizes[k|1]<csizes[k|0];
			switch(image->nch)
			{
			case 4:groups=combinations_4c;	glen=_countof(combinations_4c);	break;
			case 3:groups=combinations_3c;	glen=_countof(combinations_3c);	break;
			case 2:groups=combinations_2c;	glen=_countof(combinations_2c);	break;
			case 1:groups=0;		glen=0;				break;
			}
			if(groups)
			{
				double csize;

				comb_idx=0;

				group=groups+comb_idx;
				csize=0;
				for(int k=0;k<image->nch;++k)
					csize+=csizes[group[k]|predsel[group[k]>>1]];
				best_csize=csize;
				for(int kg=image->nch;kg<glen;kg+=image->nch)
				{
					group=groups+kg;
					csize=0;
					for(int k=0;k<image->nch;++k)
						csize+=csizes[group[k]|predsel[group[k]>>1]];
					if(best_csize>csize)
					{
						comb_idx=kg;
						best_csize=csize;
					}
				}
				
				group=groups+comb_idx;
				for(int k=0;k<image->nch;++k)
					combination[k]=group[k]|predsel[group[k]>>1];
				comb_idx/=image->nch;
			}
			else
			{
				comb_idx=0;
				combination[0]=predsel[0];
			}
			flag=comb_idx<<image->nch;
			for(int k=0;k<image->nch;++k)
				flag|=predsel[combination[k]>>1]<<k;
		}
		if(args->loud)
		{
			double avsize2=0;
			for(int k=0;k<image->nch;++k)
				avsize2+=csizes[combination[k]];
			printf("Y %5d", args->y1);
			if(image->nch>=3)
			{
				double defsize=csizes[3]+csizes[13]+csizes[7];
				printf("  default %lf (%+lf) bytes", defsize, avsize2-defsize);
			}
			printf("  current %lf bytes\n", avsize2);
			for(int k=0;k<image->nch;++k)
				printf("%s\n", prednames[combination[k]]);
			printf("\n");
		}

		dlist_init(&args->list, 1, 1024, 0);
		dlist_push_back(&args->list, &flag, 1LL+(image->nch==4));
		//{
		//	int flag=combination[3]<<3*5|combination[2]<<2*5|combination[1]<<1*5|combination[0];
		//	dlist_push_back(&args->list, &flag, sizeof(char[3]));
		//}
		gr_enc_init(&ec, &args->list);
		memset(args->pixels, 0, args->bufsize);
		for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
		{
			ALIGN(16) int *rows[]=
			{
				args->pixels+(((image->iw+16LL)*((ky-0LL)&3)+8)<<3),
				args->pixels+(((image->iw+16LL)*((ky-1LL)&3)+8)<<3),
				args->pixels+(((image->iw+16LL)*((ky-2LL)&3)+8)<<3),
				args->pixels+(((image->iw+16LL)*((ky-3LL)&3)+8)<<3),
			};
#ifdef ENABLE_WGRAD
			ALIGN(16) int X[4]={0}, Y[4]={0};
#endif
			ALIGN(16) int val[4]={0};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				int
#ifdef ENABLE_WGRAD
					*NN	=rows[2]+0*8,
					*NNE	=rows[2]+1*8,
					*NE	=rows[1]+1*8,
					*WW	=rows[0]-2*8,
#endif
					*NW	=rows[1]-1*8,
					*N	=rows[1]+0*8,
					*NEEE	=rows[1]+3*8,
					*W	=rows[0]-1*8,
					*curr	=rows[0]+0*8;
#ifdef ENABLE_WGRAD
				{
					__m128i mX=_mm_set1_epi32(1);
					__m128i mY=mX;

					__m128i mNN	=_mm_load_si128((__m128i*)NN);
					__m128i mNNE	=_mm_load_si128((__m128i*)NNE);
					__m128i mNW	=_mm_load_si128((__m128i*)NW);
					__m128i mN	=_mm_load_si128((__m128i*)N);
					__m128i mNE	=_mm_load_si128((__m128i*)NE);
					__m128i mWW	=_mm_load_si128((__m128i*)WW);
					__m128i mW	=_mm_load_si128((__m128i*)W);
					mX=_mm_add_epi32(mX, _mm_abs_epi32(_mm_sub_epi32(mW, mWW)));
					mY=_mm_add_epi32(mY, _mm_abs_epi32(_mm_sub_epi32(mW, mNW)));
					mX=_mm_add_epi32(mX, _mm_abs_epi32(_mm_sub_epi32(mN, mNW)));
					mY=_mm_add_epi32(mY, _mm_abs_epi32(_mm_sub_epi32(mN, mNN)));
					mX=_mm_add_epi32(mX, _mm_abs_epi32(_mm_sub_epi32(mNE, mN)));
					mY=_mm_add_epi32(mY, _mm_abs_epi32(_mm_sub_epi32(mNE, mNNE)));
					mX=_mm_add_epi32(mX, _mm_set1_epi32(15));
					mY=_mm_add_epi32(mY, _mm_set1_epi32(15));
					mX=_mm_srli_epi32(mX, 4);
					mY=_mm_srli_epi32(mY, 4);
					_mm_store_si128((__m128i*)X, mX);
					_mm_store_si128((__m128i*)Y, mY);
				}
#endif
				for(int kc=0;kc<image->nch;++kc)
					curr[kc]=image->data[idx+kc];
				for(int kc=0;kc<image->nch;++kc)
				{
					switch(combination[kc])
					{
#define PRED(IDX, DFLAG, VAL, EXPR) case IDX:val[kc]=(VAL)-(EXPR);break;
					PREDLIST0
					PREDLIST1
					PREDLIST2
					PREDLIST3
#undef  PRED
					}
#ifndef UNROLL_DECODER
					val[kc]=val[kc]<<1^-(val[kc]<0);
					gr_enc_POT(&ec, val[kc], FLOOR_LOG2(W[kc+4]+1));
					curr[kc+4]=(2*W[kc+4]+val[kc]+NEEE[kc+4])>>2;
#endif
				}
#ifdef UNROLL_DECODER
				{
					__m128i mval=_mm_load_si128((__m128i*)val);
					__m128i update=_mm_add_epi32(_mm_slli_epi32(_mm_load_si128((__m128i*)W+1), 1), _mm_load_si128((__m128i*)NEEE+1));
					mval=_mm_xor_si128(_mm_slli_epi32(mval, 1), _mm_srai_epi32(mval, 31));
					update=_mm_srli_epi32(_mm_add_epi32(update, mval), 2);
					_mm_store_si128((__m128i*)val, mval);
					_mm_store_si128((__m128i*)curr+1, update);
				}
				switch(image->nch)
				{
				case 4:gr_enc_POT(&ec, val[3], FLOOR_LOG2(W[3+4]+1));
				case 3:gr_enc_POT(&ec, val[2], FLOOR_LOG2(W[2+4]+1));
				case 2:gr_enc_POT(&ec, val[1], FLOOR_LOG2(W[1+4]+1));
				case 1:gr_enc_POT(&ec, val[0], FLOOR_LOG2(W[0+4]+1));
				}
#endif
				rows[0]+=8;
				rows[1]+=8;
#ifdef ENABLE_WGRAD
				rows[2]+=8;
				rows[3]+=8;
#endif
			}
		}
		gr_enc_flush(&ec);
	}
	else
	{
		Image *image=args->dst;
		int nlevels=1<<image->depth, half=nlevels>>1;
		const unsigned char *srcstart=args->decstart, *srcend=args->decend;
		int combination[4]={0};
		unsigned short flag=0;

		memcpy(&flag, srcstart, 1LL+(image->nch==4));
		srcstart+=1LL+(image->nch==4);
		{
			const int *groups=0;

			int idx=image->nch*(flag>>image->nch);
			switch(image->nch)
			{
			case 4:groups=combinations_4c;	break;
			case 3:groups=combinations_3c;	break;
			case 2:groups=combinations_2c;	break;
			case 1:groups=0;		break;
			}
			if(groups)
			{
				memcpy(combination, groups+idx, sizeof(int)*image->nch);
				for(int k=0;k<image->nch;++k)
					combination[k]|=flag>>k&1;
			}
			else
			{
				combination[0]=flag;
			}
		}

		//combination[0]=flag>>0*5&31;
		//combination[1]=flag>>1*5&31;
		//combination[2]=flag>>2*5&31;
		//combination[3]=flag>>3*5&31;

		gr_dec_init(&ec, srcstart, srcend);
		for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
		{
			ALIGN(16) int *rows[]=
			{
				args->pixels+(((image->iw+16LL)*((ky-0LL)&3)+8)<<3),
				args->pixels+(((image->iw+16LL)*((ky-1LL)&3)+8)<<3),
				args->pixels+(((image->iw+16LL)*((ky-2LL)&3)+8)<<3),
				args->pixels+(((image->iw+16LL)*((ky-3LL)&3)+8)<<3),
			};
#ifdef ENABLE_WGRAD
			ALIGN(16) int X[4]={0}, Y[4]={0};
#endif
			ALIGN(16) int val[4]={0};

			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				int
#ifdef ENABLE_WGRAD
					*NN	=rows[2]+0*8,
					*NNE	=rows[2]+1*8,
					*NE	=rows[1]+1*8,
					*WW	=rows[0]-2*8,
#endif
					*NW	=rows[1]-1*8,
					*N	=rows[1]+0*8,
					*NEEE	=rows[1]+3*8,
					*W	=rows[0]-1*8,
					*curr	=rows[0]+0*8;
#ifdef ENABLE_WGRAD
				{
					__m128i mX=_mm_set1_epi32(1);
					__m128i mY=mX;
				
					__m128i mNN	=_mm_load_si128((__m128i*)NN);
					__m128i mNNE	=_mm_load_si128((__m128i*)NNE);
					__m128i mNW	=_mm_load_si128((__m128i*)NW);
					__m128i mN	=_mm_load_si128((__m128i*)N);
					__m128i mNE	=_mm_load_si128((__m128i*)NE);
					__m128i mWW	=_mm_load_si128((__m128i*)WW);
					__m128i mW	=_mm_load_si128((__m128i*)W);
					mX=_mm_add_epi32(mX, _mm_abs_epi32(_mm_sub_epi32(mW, mWW)));
					mY=_mm_add_epi32(mY, _mm_abs_epi32(_mm_sub_epi32(mW, mNW)));
					mX=_mm_add_epi32(mX, _mm_abs_epi32(_mm_sub_epi32(mN, mNW)));
					mY=_mm_add_epi32(mY, _mm_abs_epi32(_mm_sub_epi32(mN, mNN)));
					mX=_mm_add_epi32(mX, _mm_abs_epi32(_mm_sub_epi32(mNE, mN)));
					mY=_mm_add_epi32(mY, _mm_abs_epi32(_mm_sub_epi32(mNE, mNNE)));
					mX=_mm_add_epi32(mX, _mm_set1_epi32(15));
					mY=_mm_add_epi32(mY, _mm_set1_epi32(15));
					mX=_mm_srli_epi32(mX, 4);
					mY=_mm_srli_epi32(mY, 4);
					_mm_store_si128((__m128i*)X, mX);
					_mm_store_si128((__m128i*)Y, mY);
				}
#endif
#ifdef UNROLL_DECODER
				switch(image->nch)
				{
#if 0
				case 4:
					{
						GolombRiceCoder *pc=&ec;
						int nbypass=FLOOR_LOG2(W[3+4]+1);
						int sym=0, nleadingzeros=0;
						if(!pc->nbits)//cache is empty
							goto read_c3;
						for(;;)//cache reading loop
						{
							nleadingzeros=pc->nbits-FLOOR_LOG2_P1(pc->cache);//count leading zeros
							pc->nbits-=nleadingzeros;//remove accessed zeros
							sym+=nleadingzeros;

							if(pc->nbits)
								break;
						read_c3://cache is empty
							gr_dec_impl_read(pc);
						}
						//now  0 < nbits <= 64
						--pc->nbits;
						//now  0 <= nbits < 64
						pc->cache-=1ULL<<pc->nbits;//remove stop bit

						unsigned bypass=0;
						sym<<=nbypass;
						if(pc->nbits<nbypass)
						{
							nbypass-=pc->nbits;
							bypass|=(int)(pc->cache<<nbypass);
							gr_dec_impl_read(pc);
						}
						if(nbypass)
						{
							pc->nbits-=nbypass;
							bypass|=(int)(pc->cache>>pc->nbits);
							pc->cache&=(1ULL<<pc->nbits)-1;
						}
						val[3]=sym|bypass;
					}
				case 3:
					{
						GolombRiceCoder *pc=&ec;
						int nbypass=FLOOR_LOG2(W[2+4]+1);
						int sym=0, nleadingzeros=0;
						if(!pc->nbits)//cache is empty
							goto read_c2;
						for(;;)//cache reading loop
						{
							nleadingzeros=pc->nbits-FLOOR_LOG2_P1(pc->cache);//count leading zeros
							pc->nbits-=nleadingzeros;//remove accessed zeros
							sym+=nleadingzeros;

							if(pc->nbits)
								break;
						read_c2://cache is empty
							gr_dec_impl_read(pc);
						}
						//now  0 < nbits <= 64
						--pc->nbits;
						//now  0 <= nbits < 64
						pc->cache-=1ULL<<pc->nbits;//remove stop bit

						unsigned bypass=0;
						sym<<=nbypass;
						if(pc->nbits<nbypass)
						{
							nbypass-=pc->nbits;
							bypass|=(int)(pc->cache<<nbypass);
							gr_dec_impl_read(pc);
						}
						if(nbypass)
						{
							pc->nbits-=nbypass;
							bypass|=(int)(pc->cache>>pc->nbits);
							pc->cache&=(1ULL<<pc->nbits)-1;
						}
						val[2]=sym|bypass;
					}
				case 2:
					{
						GolombRiceCoder *pc=&ec;
						int nbypass=FLOOR_LOG2(W[1+4]+1);
						int sym=0, nleadingzeros=0;
						if(!pc->nbits)//cache is empty
							goto read_c1;
						for(;;)//cache reading loop
						{
							nleadingzeros=pc->nbits-FLOOR_LOG2_P1(pc->cache);//count leading zeros
							pc->nbits-=nleadingzeros;//remove accessed zeros
							sym+=nleadingzeros;

							if(pc->nbits)
								break;
						read_c1://cache is empty
							gr_dec_impl_read(pc);
						}
						//now  0 < nbits <= 64
						--pc->nbits;
						//now  0 <= nbits < 64
						pc->cache-=1ULL<<pc->nbits;//remove stop bit

						unsigned bypass=0;
						sym<<=nbypass;
						if(pc->nbits<nbypass)
						{
							nbypass-=pc->nbits;
							bypass|=(int)(pc->cache<<nbypass);
							gr_dec_impl_read(pc);
						}
						if(nbypass)
						{
							pc->nbits-=nbypass;
							bypass|=(int)(pc->cache>>pc->nbits);
							pc->cache&=(1ULL<<pc->nbits)-1;
						}
						val[1]=sym|bypass;
					}
				case 1:
					{
						GolombRiceCoder *pc=&ec;
						int nbypass=FLOOR_LOG2(W[0+4]+1);
						int sym=0, nleadingzeros=0;
						if(!pc->nbits)//cache is empty
							goto read_c0;
						for(;;)//cache reading loop
						{
							nleadingzeros=pc->nbits-FLOOR_LOG2_P1(pc->cache);//count leading zeros
							pc->nbits-=nleadingzeros;//remove accessed zeros
							sym+=nleadingzeros;

							if(pc->nbits)
								break;
						read_c0://cache is empty
							gr_dec_impl_read(pc);
						}
						//now  0 < nbits <= 64
						--pc->nbits;
						//now  0 <= nbits < 64
						pc->cache-=1ULL<<pc->nbits;//remove stop bit

						unsigned bypass=0;
						sym<<=nbypass;
						if(pc->nbits<nbypass)
						{
							nbypass-=pc->nbits;
							bypass|=(int)(pc->cache<<nbypass);
							gr_dec_impl_read(pc);
						}
						if(nbypass)
						{
							pc->nbits-=nbypass;
							bypass|=(int)(pc->cache>>pc->nbits);
							pc->cache&=(1ULL<<pc->nbits)-1;
						}
						val[0]=sym|bypass;
					}
#else
				case 4:val[3]=gr_dec_POT(&ec, FLOOR_LOG2(W[3+4]+1));
				case 3:val[2]=gr_dec_POT(&ec, FLOOR_LOG2(W[2+4]+1));
				case 2:val[1]=gr_dec_POT(&ec, FLOOR_LOG2(W[1+4]+1));
				case 1:val[0]=gr_dec_POT(&ec, FLOOR_LOG2(W[0+4]+1));
#endif
				}
				{
					__m128i mval=_mm_load_si128((__m128i*)val);
					__m128i update=_mm_add_epi32(_mm_slli_epi32(_mm_load_si128((__m128i*)W+1), 1), _mm_load_si128((__m128i*)NEEE+1));
					__m128i mask=_mm_cmpeq_epi32(_mm_and_si128(mval, _mm_set1_epi32(1)), _mm_set1_epi32(1));
					update=_mm_srli_epi32(_mm_add_epi32(update, mval), 2);
					mval=_mm_xor_si128(_mm_srli_epi32(mval, 1), mask);
					_mm_store_si128((__m128i*)curr+1, update);
					_mm_store_si128((__m128i*)val, mval);
				}
#endif
				//for(int kc=0;kc<image->nch;++kc)
				//{
				//	val[kc]=gr_dec_POT(&ec, FLOOR_LOG2(W[kc+4]+1));
				//	curr[kc+4]=(2*W[kc+4]+val[kc]+NEEE[kc+4])>>2;
				//	val[kc]=val[kc]>>1^-(val[kc]&1);
				//}
				for(int kc=0;kc<image->nch;++kc)
				{
#ifndef UNROLL_DECODER
					val[kc]=gr_dec_POT(&ec, FLOOR_LOG2(W[kc+4]+1));
					curr[kc+4]=(2*W[kc+4]+val[kc]+NEEE[kc+4])>>2;
					val[kc]=val[kc]>>1^-(val[kc]&1);
#endif
					switch(combination[kc])
					{
#define PRED(IDX, DFLAG, VAL, EXPR) case IDX:VAL=val[kc]+EXPR;break;
					PREDLIST0
					PREDLIST1
					PREDLIST2
					PREDLIST3
#undef  PRED
					}
				}
				for(int kc=0;kc<image->nch;++kc)
					image->data[idx+kc]=curr[kc];
#ifdef ENABLE_GUIDE
				if(memcmp(image->data+idx, guide->data+idx, sizeof(short)*image->nch))
				{
					short orig[4]={0};
					memcpy(orig, guide->data+idx, image->nch*sizeof(short));
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
				rows[0]+=8;
				rows[1]+=8;
#ifdef ENABLE_WGRAD
				rows[2]+=8;
				rows[3]+=8;
#endif
			}
		}
	}
}
int f23_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	//test4cycles();

	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif
	int ncores=query_cpu_cores();
	int nblocks=(image->ih+BLOCKSIZE-1)/BLOCKSIZE, nthreads=MINVAR(nblocks, ncores);
	ptrdiff_t memusage=0;
	ptrdiff_t argssize=nthreads*sizeof(ThreadArgs);
	ThreadArgs *args=(ThreadArgs*)malloc(argssize);
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memusage+=argssize;
	memset(args, 0, argssize);
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=src;
		arg->dst=dst;
		arg->bufsize=sizeof(int[4*4*2])*(image->iw+16LL);//4 padded rows * 4 channels max * {pixels, errors}
		arg->pixels=(int*)_mm_malloc(arg->bufsize, sizeof(__m128i));
		arg->histsize=sizeof(int[MAXPREDS])<<image->depth;
		arg->hist=(int*)malloc(arg->histsize);
		if(!arg->pixels||!arg->hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=(ptrdiff_t)arg->bufsize+arg->histsize;
		arg->fwd=fwd;
#ifdef DISABLE_MT
		arg->loud=loud;
#else
		arg->loud=0;
#endif
	}
	if(fwd)
	{
		ptrdiff_t dststart=array_append(data, 0, 1, sizeof(int)*nblocks, 1, 0, 0);
		for(int kt=0;kt<nblocks;kt+=nthreads)
		{
			int nthreads2=MINVAR(kt+nthreads, nblocks)-kt;
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				arg->y1=BLOCKSIZE*(kt+kt2);
				arg->y2=MINVAR(arg->y1+BLOCKSIZE, image->ih);
			}
#ifdef DISABLE_MT
			for(int k=0;k<nthreads2;++k)
				block_thread(args+k);
#else
			void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
			mt_finish(ctx);
#endif
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				if(loud)
					printf("[%d]  %zd\n", kt+kt2, args[kt2].list.nobj);
				memcpy(data[0]->data+dststart+sizeof(int)*((ptrdiff_t)kt+kt2), &args[kt2].list.nobj, sizeof(int));
				dlist_appendtoarray(&args[kt2].list, data);
				dlist_clear(&args[kt2].list);
			}
		}
		if(loud)
		{
			ptrdiff_t
				csize=data[0]->count-dststart,
				usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
			t0=time_sec()-t0;
			printf("Size %14td/%14td  %16lf%%  %16lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
			printf("Mem usage: ");
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
			printf("E %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
	}
	else
	{
		const unsigned char *dstptr=cbuf+sizeof(int)*nblocks;
		int dec_offset=0;

		//integrity check
#if 1
		for(int kt=0;kt<nblocks;++kt)
		{
			int size=0;
			memcpy(&size, cbuf+sizeof(int)*kt, sizeof(int));
			dec_offset+=size;
		}
		if(sizeof(int)*nblocks+dec_offset!=clen)
			LOG_ERROR("Corrupt file");
#endif
		dec_offset=0;
		for(int kt=0;kt<nblocks;kt+=nthreads)
		{
			int nthreads2=MINVAR(kt+nthreads, nblocks)-kt;
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				int size=0;
				memcpy(&size, cbuf+sizeof(int)*((ptrdiff_t)kt+kt2), sizeof(int));
				arg->y1=BLOCKSIZE*(kt+kt2);
				arg->y2=MINVAR(arg->y1+BLOCKSIZE, image->ih);
				arg->decstart=dstptr+dec_offset;
				dec_offset+=size;
				arg->decend=dstptr+dec_offset;
			}
#ifdef DISABLE_MT
			for(int k=0;k<nthreads2;++k)
				block_thread(args+k);
#else
			{
				void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
				mt_finish(ctx);
			}
#endif
		}
		if(loud)
		{
			ptrdiff_t usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
			t0=time_sec()-t0;
			printf("D %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		_mm_free(arg->pixels);
		free(arg->hist);
	}
	free(args);
	return 0;
}