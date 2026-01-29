#ifdef _MSC_VER
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif
#elif defined __linux__ && !defined _GNU_SOURCE
#	define _GNU_SOURCE
#	include<stddef.h>//ptrdiff_t
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<math.h>
#include<sys/stat.h>
#if defined _MSC_VER || defined _WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#else
#include<time.h>
#endif
#include<immintrin.h>


#ifdef _MSC_VER
	#define LOUD
	#define PROFILE_SIZE
	#define ENABLE_GUIDE
	#define FIFOVAL
#endif


//	#define USE_DCT8
	#define USE_L1


//c32
#define PREDLIST\
	PRED(N)\
	PRED(W)\
	PRED(3*(N-NN)+NNN)\
	PRED(3*(W-WW)+WWW)\
	PRED(W+NE-N)\
	PRED((WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)\
	PRED(N+W-NW)\
	PRED(N+NE-NNE)\

enum
{
#ifdef USE_DCT8
	BLOCKX=8,
	BLOCKY=8,
	DCT_ROUND_TRIP_GAIN=16,
#else
	BLOCKX=4,
	BLOCKY=4,
	DCT_ROUND_TRIP_GAIN=4,
#endif

	SHIFT=18,
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED

	XPAD=8,
	NCH=3,
	NROWS=4,
	NVAL=2,

	GRBITS=3,

	NCTX_DC=13,
	NCTX_AC=(BLOCKX+BLOCKY)/2,
	NCTX=NCTX_DC+NCTX_AC,
#ifdef USE_DCT8
	NLEVELS=224,
#else
	NLEVELS=72,
#endif
};

//runtime
#if 1
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
#ifdef _MSC_VER
#	define	ALIGN(N) __declspec(align(N))
#	define AWM_INLINE __forceinline static
#else
#	define	ALIGN(N) __attribute__((aligned(N)))
#	define AWM_INLINE __attribute__((always_inline)) inline static
#	ifndef _countof
#		define _countof(A) (sizeof(A)/sizeof(*(A)))
#	endif
#endif
#if defined _M_X64 || defined __x86_64__
#define FLOOR_LOG2(X)\
	(sizeof(X)==8?63-(int32_t)_lzcnt_u64(X):31-_lzcnt_u32((uint32_t)(X)))
#else
AWM_INLINE int floor_log2_64(uint64_t n)
{
	int	logn=-!n;
	int	sh=(n>=1ULL<<32)<<5;	logn+=sh, n>>=sh;
		sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
}
AWM_INLINE int floor_log2_32(uint32_t n)
{
	int	logn=-!n;
	int	sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
}
#define FLOOR_LOG2(X)\
	(sizeof(X)==8?floor_log2_64(X):floor_log2_32((uint32_t)(X)))
#endif
#define CVTFP32_I32(X)  _mm_cvt_ss2si(_mm_set_ss(X))
#define CVTTFP32_I32(X) _mm_cvtt_ss2si(_mm_set_ss(X))
#define CVTFP64_I64(X)  _mm_cvtsd_si64(_mm_set_sd(X))
#define CVTTFP64_I64(X) _mm_cvttsd_si64(_mm_set_sd(X))
static void crash(const char *file, int line, const char *format, ...)
{
	printf("%s(%d):\n", file, line);
	if(format)
	{
		va_list args;
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
	}
	printf("\n");
	exit(1);
}
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT,##__VA_ARGS__)
static double time_sec(void)
{
#ifdef _WIN32
	static long long t0=0;
	LARGE_INTEGER li;
	double t;
	QueryPerformanceCounter(&li);
	if(!t0)
		t0=li.QuadPart;
	t=(double)(li.QuadPart-t0);
	QueryPerformanceFrequency(&li);
	t/=(double)li.QuadPart;
	return t;
#else
	struct timespec t;
	clock_gettime(CLOCK_REALTIME, &t);//<time.h>
	return t.tv_sec+t.tv_nsec*1e-9;
#endif
}
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static uint8_t *g_image=0;
static double g_sqe[3]={0};
static void guide_save(uint8_t *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(uint8_t*)malloc(size);
	if(!g_image)
	{
		CRASH("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
#endif
#endif


#ifdef FIFOVAL
static ptrdiff_t fifoidx=0, fifocap=0, fifoidx2=0;
static uint32_t *fifoval=0;
static void valfifo_enqueue(uint32_t val)
{
	if(fifoidx+1>=fifocap)
	{
		void *p=0;

		if(!fifocap)
			fifocap=1;
		fifocap<<=1;
		p=realloc(fifoval, fifocap*sizeof(uint32_t));
		if(!p)
		{
			CRASH("Alloc error");
			return;
		}
		fifoval=(uint32_t*)p;
	}
	fifoval[fifoidx++]=val;
}
static void valfifo_check(uint32_t val)
{
	uint32_t val0=fifoval[fifoidx2++];
	if(val!=val0)
	{
		--fifoidx2;
		printf(
			"\n"
			"FIFO Error  at %10lld,  remaining %10lld\n"
			"    0x%08X  !=  original 0x%08X\n"
			"\n"
			, fifoidx2
			, fifoidx-fifoidx2//current element was not decoded successfully
			, val, val0
		);
		for(int k=-32;k<32;++k)
		{
			ptrdiff_t idx=fifoidx2+k;
			if((size_t)idx<(size_t)fifoidx)
			{
				printf(
					"%10td  0x%08X"
					, idx
					, fifoval[idx]
				);
				if(idx<fifoidx2)
					printf("  OK");
				if(idx==fifoidx2)
					printf("  !=  corrupt 0x%08X", val);
				printf("\n");
			}
		}
		CRASH("");
	}
}
#endif


#ifdef _MSC_VER
static double sum_before[3][BLOCKX*BLOCKY], sum_after[3][BLOCKX*BLOCKY];
static int cmin[3][BLOCKX*BLOCKY], cmax[3][BLOCKX*BLOCKY];
#endif
static const uint8_t context_table[]=
{
	0, 1, 2, 3,
	1, 2, 3, 4,
	2, 3, 4, 4,
	4, 4, 4, 4,
};
static const float qtable_luma[]=//check NLEVELS
{
#ifdef USE_DCT8
	 5,  3,  3,  5,  7, 12, 15, 18,
	 3,  3,  4,  6,  8, 17, 18, 16,
	 4,  4,  5,  7, 12, 17, 21, 17,
	 4,  5,  7,  9, 15, 26, 24, 19,
	 5,  7, 11, 17, 20, 33, 31, 23,
	 7, 11, 17, 19, 24, 31, 34, 28,
	15, 19, 23, 26, 31, 36, 36, 30,
	22, 28, 29, 31, 34, 30, 31, 30,

	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
#else
	16, 12, 12, 16,//v17
	12, 12, 16, 16,
	12, 16, 16, 16,
	16, 16, 16, 16,

	//16, 12, 12, 12,//v16
	//12, 10, 12, 12,
	//12, 12, 12, 12,
	//12, 12, 12, 12,

	//16, 10, 12, 16,//v15
	//10, 10, 16, 24,
	//12, 16, 24, 24,
	//16, 24, 24, 24,

	//10, 12, 14, 16,//v14
	//12, 20, 32, 32,
	//14, 32, 32, 32,
	//16, 32, 32, 32,

	//10, 24, 28, 32,
	//24, 28, 32, 32,
	//28, 32, 32, 32,
	//32, 32, 32, 32,

	//16, 18, 18, 14,//nice
	//18, 18, 14, 12,
	//18, 14, 12, 12,
	//14, 12, 12, 12,

	//16, 16, 14, 12,//good
	//16, 14, 12, 12,
	//14, 12, 12, 12,
	//12, 12, 12, 12,

	//10, 7, 7, 7,
	// 7, 7, 7, 7,
	// 7, 7, 7, 7,
	// 7, 7, 7, 7,

	//12, 13, 14, 15,
	//13, 14, 15, 16,
	//14, 15, 16, 16,
	//15, 16, 16, 16,

	//12, 12, 14, 16,
	//12, 14, 16, 32,
	//14, 16, 24, 32,
	//16, 32, 32, 32,

	//5, 6, 7, 8,//CRASH overflow
	//6, 7, 8, 9,
	//7, 8, 8, 9,
	//8, 9, 9, 9,

	//5, 3, 4, 5,
	//3, 4, 5, 6,
	//4, 5, 5, 6,
	//5, 6, 6, 6,

	//16, 24, 32, 48,
	//24, 32, 48, 64,
	//32, 48, 48, 64,
	//48, 64, 64, 64,
#endif
};
static const float qtable_chroma[]=
{
#ifdef USE_DCT8
	 5,  5,  7, 14, 30, 30, 30, 30,
	 5,  6,  8, 20, 30, 30, 30, 30,
	 7,  8, 17, 30, 30, 30, 30, 30,
	14, 20, 30, 30, 30, 30, 30, 30,
	30, 30, 30, 30, 30, 30, 30, 30,
	30, 30, 30, 30, 30, 30, 30, 30,
	30, 30, 30, 30, 30, 30, 30, 30,
	30, 30, 30, 30, 30, 30, 30, 30,

	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
	//64, 64, 64, 64, 64, 64, 64, 64,
#else
	16, 10, 12, 12,//v14
	10, 10, 12, 12,
	12, 12, 12, 12,
	12, 12, 12, 12,

	//10, 12, 14, 16,
	//12, 20, 32, 32,
	//14, 32, 32, 32,
	//16, 32, 32, 32,

	//12, 12, 10, 10,
	//12, 10, 10, 10,
	//10, 10, 10, 10,
	//10, 10, 10, 10,

	//12, 12, 13, 13,
	//12, 13, 15, 16,
	//13, 15, 16, 16,
	//13, 16, 16, 16,

	//10, 7, 7, 7,
	// 7, 7, 7, 7,
	// 7, 7, 7, 7,
	// 7, 7, 7, 7,

	//10, 4, 3, 2,
	// 4, 3, 2, 1,
	// 3, 2, 1, 1,
	// 2, 1, 1, 1,

	//12, 13, 14, 15,
	//13, 14, 15, 15,
	//14, 15, 15, 15,
	//15, 15, 15, 15,

	//12, 12, 14, 16,
	//12, 14, 16, 16,
	//14, 16, 16, 16,
	//16, 16, 16, 16,

	//5, 6, 7, 8,//CRASH overflow
	//6, 7, 8, 8,
	//7, 8, 8, 8,
	//8, 8, 8, 8,

	//5, 3, 4, 4,
	//3, 4, 4, 4,
	//4, 4, 4, 4,
	//4, 4, 4, 4,

	//16, 24, 32, 32,
	//24, 32, 32, 32,
	//32, 32, 32, 32,
	//32, 32, 32, 32,
#endif
};
#ifdef USE_DCT8
AWM_INLINE void cvti2f(float *block)
{
	__m256i a[8];

	a[0]=_mm256_load_si256((__m256i*)block+0);
	a[1]=_mm256_load_si256((__m256i*)block+1);
	a[2]=_mm256_load_si256((__m256i*)block+2);
	a[3]=_mm256_load_si256((__m256i*)block+3);
	a[4]=_mm256_load_si256((__m256i*)block+4);
	a[5]=_mm256_load_si256((__m256i*)block+5);
	a[6]=_mm256_load_si256((__m256i*)block+6);
	a[7]=_mm256_load_si256((__m256i*)block+7);

	_mm256_store_ps(block+0*8, _mm256_cvtepi32_ps(a[0]));
	_mm256_store_ps(block+1*8, _mm256_cvtepi32_ps(a[1]));
	_mm256_store_ps(block+2*8, _mm256_cvtepi32_ps(a[2]));
	_mm256_store_ps(block+3*8, _mm256_cvtepi32_ps(a[3]));
	_mm256_store_ps(block+4*8, _mm256_cvtepi32_ps(a[4]));
	_mm256_store_ps(block+5*8, _mm256_cvtepi32_ps(a[5]));
	_mm256_store_ps(block+6*8, _mm256_cvtepi32_ps(a[6]));
	_mm256_store_ps(block+7*8, _mm256_cvtepi32_ps(a[7]));
}
AWM_INLINE void cvtf2i(float *block)
{
	__m256 a[8];

	a[0]=_mm256_load_ps(block+0*8);
	a[1]=_mm256_load_ps(block+1*8);
	a[2]=_mm256_load_ps(block+2*8);
	a[3]=_mm256_load_ps(block+3*8);
	a[4]=_mm256_load_ps(block+4*8);
	a[5]=_mm256_load_ps(block+5*8);
	a[6]=_mm256_load_ps(block+6*8);
	a[7]=_mm256_load_ps(block+7*8);

	_mm256_store_si256((__m256i*)block+0, _mm256_cvtps_epi32(a[0]));
	_mm256_store_si256((__m256i*)block+1, _mm256_cvtps_epi32(a[1]));
	_mm256_store_si256((__m256i*)block+2, _mm256_cvtps_epi32(a[2]));
	_mm256_store_si256((__m256i*)block+3, _mm256_cvtps_epi32(a[3]));
	_mm256_store_si256((__m256i*)block+4, _mm256_cvtps_epi32(a[4]));
	_mm256_store_si256((__m256i*)block+5, _mm256_cvtps_epi32(a[5]));
	_mm256_store_si256((__m256i*)block+6, _mm256_cvtps_epi32(a[6]));
	_mm256_store_si256((__m256i*)block+7, _mm256_cvtps_epi32(a[7]));
}
AWM_INLINE void rgb2yuv(float *c0, float *c1, float *c2)
{
	__m256 half=_mm256_set1_ps(128);
	for(int k=0;k<8;++k)
	{
		__m256 rgb[3], yuv[3];

		rgb[0]=_mm256_load_ps(c0+k*8);
		rgb[1]=_mm256_load_ps(c1+k*8);
		rgb[2]=_mm256_load_ps(c2+k*8);
		
		//https://en.wikipedia.org/wiki/YCbCr
		yuv[0]=_mm256_mul_ps(_mm256_set1_ps(+0.299000f), rgb[0]);
		yuv[1]=_mm256_mul_ps(_mm256_set1_ps(-0.168736f), rgb[0]);
		yuv[2]=_mm256_mul_ps(_mm256_set1_ps(+0.500000f), rgb[0]);
		yuv[0]=_mm256_add_ps(yuv[0], _mm256_mul_ps(_mm256_set1_ps(+0.587000f), rgb[1]));
		yuv[1]=_mm256_add_ps(yuv[1], _mm256_mul_ps(_mm256_set1_ps(-0.331264f), rgb[1]));
		yuv[2]=_mm256_add_ps(yuv[2], _mm256_mul_ps(_mm256_set1_ps(-0.418688f), rgb[1]));
		yuv[0]=_mm256_add_ps(yuv[0], _mm256_mul_ps(_mm256_set1_ps(+0.114000f), rgb[2]));
		yuv[1]=_mm256_add_ps(yuv[1], _mm256_mul_ps(_mm256_set1_ps(+0.500000f), rgb[2]));
		yuv[2]=_mm256_add_ps(yuv[2], _mm256_mul_ps(_mm256_set1_ps(-0.081312f), rgb[2]));
		
		yuv[0]=_mm256_sub_ps(yuv[0], half);

		_mm256_store_ps(c0+k*8, yuv[0]);
		_mm256_store_ps(c1+k*8, yuv[1]);
		_mm256_store_ps(c2+k*8, yuv[2]);
	}
}
AWM_INLINE void yuv2rgb(float *c0, float *c1, float *c2)
{
	__m256 half=_mm256_set1_ps(128);
	__m256 vmin=_mm256_setzero_ps();
	__m256 vmax=_mm256_set1_ps(255);
	for(int k=0;k<8;++k)
	{
		__m256 rgb[3], yuv[3];

		yuv[0]=_mm256_load_ps(c0+k*8);
		yuv[1]=_mm256_load_ps(c1+k*8);
		yuv[2]=_mm256_load_ps(c2+k*8);
		
		yuv[0]=_mm256_add_ps(yuv[0], half);

		rgb[1]=_mm256_add_ps(yuv[0], _mm256_mul_ps(_mm256_set1_ps(-0.344136f), yuv[1]));
		rgb[2]=_mm256_add_ps(yuv[0], _mm256_mul_ps(_mm256_set1_ps(+1.772000f), yuv[1]));
		rgb[0]=_mm256_add_ps(yuv[0], _mm256_mul_ps(_mm256_set1_ps(+1.402000f), yuv[2]));
		rgb[1]=_mm256_add_ps(rgb[1], _mm256_mul_ps(_mm256_set1_ps(-0.714136f), yuv[2]));

		rgb[0]=_mm256_min_ps(rgb[0], vmax);
		rgb[1]=_mm256_min_ps(rgb[1], vmax);
		rgb[2]=_mm256_min_ps(rgb[2], vmax);
		rgb[0]=_mm256_max_ps(rgb[0], vmin);
		rgb[1]=_mm256_max_ps(rgb[1], vmin);
		rgb[2]=_mm256_max_ps(rgb[2], vmin);
		_mm256_store_ps(c0+k*8, rgb[0]);
		_mm256_store_ps(c1+k*8, rgb[1]);
		_mm256_store_ps(c2+k*8, rgb[2]);
	}
}
AWM_INLINE void dct_fwd(float *block)//DCT-II
{
	/*
	DCTii8 =

	DCTii4		0
	0		C * DCTii4
	
	where C = diag(cos(pi*1/16), cos(pi*3/16), cos(pi*5/16), cos(pi*7/16))
	*/

	//Chen factorization
	__m256 c7=_mm256_set1_ps(0.1950903220161282678482848684770f);//cosd(90*7/8) = sind(90*1/8) = s1
	__m256 c6=_mm256_set1_ps(0.3826834323650897717284599840304f);//cosd(90*6/8) = sind(90*2/8) = s2
	__m256 c5=_mm256_set1_ps(0.5555702330196022247428308139485f);//cosd(90*5/8) = sind(90*3/8) = s3
	__m256 c4=_mm256_set1_ps(0.7071067811865475244008443621048f);//cosd(90*4/8) = sind(90*4/8) = s4 = sqrt(2)
	__m256 c3=_mm256_set1_ps(0.8314696123025452370787883776179f);//cosd(90*3/8) = sind(90*5/8) = s5
	__m256 c2=_mm256_set1_ps(0.9238795325112867561281831893968f);//cosd(90*2/8) = sind(90*6/8) = s6
	__m256 c1=_mm256_set1_ps(0.9807852804032304491261822361342f);//cosd(90*1/8) = sind(90*7/8) = s7

	__m256 a[8], b[8];

	a[0]=_mm256_load_ps(block+0*8);
	a[1]=_mm256_load_ps(block+1*8);
	a[2]=_mm256_load_ps(block+2*8);
	a[3]=_mm256_load_ps(block+3*8);
	a[4]=_mm256_load_ps(block+4*8);
	a[5]=_mm256_load_ps(block+5*8);
	a[6]=_mm256_load_ps(block+6*8);
	a[7]=_mm256_load_ps(block+7*8);

	b[0]=_mm256_add_ps(a[0], a[7]);//step 1
	b[1]=_mm256_add_ps(a[1], a[6]);
	b[2]=_mm256_add_ps(a[2], a[5]);
	b[3]=_mm256_add_ps(a[3], a[4]);
	b[4]=_mm256_sub_ps(a[3], a[4]);
	b[5]=_mm256_sub_ps(a[2], a[5]);
	b[6]=_mm256_sub_ps(a[1], a[6]);
	b[7]=_mm256_sub_ps(a[0], a[7]);
	
	a[0]=_mm256_add_ps(b[0], b[3]);//step 2
	a[1]=_mm256_add_ps(b[1], b[2]);
	a[2]=_mm256_sub_ps(b[1], b[2]);
	a[3]=_mm256_sub_ps(b[0], b[3]);
	a[4]=b[4];
	a[5]=_mm256_mul_ps(_mm256_sub_ps(b[6], b[5]), c4);//M0
	a[6]=_mm256_mul_ps(_mm256_add_ps(b[6], b[5]), c4);
	a[7]=b[7];
	
	b[0]=_mm256_add_ps(a[0], a[1]);//step 3
	b[1]=_mm256_mul_ps(_mm256_sub_ps(a[0], a[1]), c4);
	b[2]=_mm256_add_ps(_mm256_mul_ps(a[3], c2), _mm256_mul_ps(a[2], c6));//M1
	b[3]=_mm256_sub_ps(_mm256_mul_ps(a[3], c6), _mm256_mul_ps(a[2], c2));
	b[4]=_mm256_add_ps(a[4], a[5]);
	b[5]=_mm256_sub_ps(a[4], a[5]);
	b[6]=_mm256_sub_ps(a[7], a[6]);
	b[7]=_mm256_add_ps(a[7], a[6]);
	
	a[0]=b[0];//step 4
	a[1]=b[1];
	a[2]=b[2];
	a[3]=b[3];
	a[4]=_mm256_add_ps(_mm256_mul_ps(b[7], c1), _mm256_mul_ps(b[4], c7));//M2
	a[5]=_mm256_add_ps(_mm256_mul_ps(b[6], c5), _mm256_mul_ps(b[5], c3));//M3
	a[6]=_mm256_sub_ps(_mm256_mul_ps(b[6], c3), _mm256_mul_ps(b[5], c5));
	a[7]=_mm256_sub_ps(_mm256_mul_ps(b[7], c7), _mm256_mul_ps(b[4], c1));

	_mm256_store_ps(block+0*8, a[0]);
	_mm256_store_ps(block+4*8, a[1]);
	_mm256_store_ps(block+2*8, a[2]);
	_mm256_store_ps(block+6*8, a[3]);
	_mm256_store_ps(block+1*8, a[4]);
	_mm256_store_ps(block+5*8, a[5]);
	_mm256_store_ps(block+3*8, a[6]);
	_mm256_store_ps(block+7*8, a[7]);
}
AWM_INLINE void dct_inv(float *block)//DCT-III
{
	//Chen factorization
	__m256 c7=_mm256_set1_ps(0.1950903220161282678482848684770f);//cosd(90*7/8) = sind(90*1/8) = s1
	__m256 c6=_mm256_set1_ps(0.3826834323650897717284599840304f);//cosd(90*6/8) = sind(90*2/8) = s2
	__m256 c5=_mm256_set1_ps(0.5555702330196022247428308139485f);//cosd(90*5/8) = sind(90*3/8) = s3
	__m256 c4=_mm256_set1_ps(0.7071067811865475244008443621048f);//cosd(90*4/8) = sind(90*4/8) = s4 = sqrt(2)
	__m256 c3=_mm256_set1_ps(0.8314696123025452370787883776179f);//cosd(90*3/8) = sind(90*5/8) = s5
	__m256 c2=_mm256_set1_ps(0.9238795325112867561281831893968f);//cosd(90*2/8) = sind(90*6/8) = s6
	__m256 c1=_mm256_set1_ps(0.9807852804032304491261822361342f);//cosd(90*1/8) = sind(90*7/8) = s7

	__m256 a[8], b[8];

	b[0]=_mm256_load_ps(block+0*8);
	b[1]=_mm256_load_ps(block+4*8);
	b[2]=_mm256_load_ps(block+2*8);
	b[3]=_mm256_load_ps(block+6*8);
	b[4]=_mm256_load_ps(block+1*8);
	b[5]=_mm256_load_ps(block+5*8);
	b[6]=_mm256_load_ps(block+3*8);
	b[7]=_mm256_load_ps(block+7*8);

	/*
	matrix					inverse
	M0 = [c4 c4;c4 -c4]			itself (scaled)
	M1: [a2;a3] = [c6 c2;-c2 c6][b2;b3]	transpose [b2;b3] = [c6 -c2;c2 c6][a2;a3]
	M2: [a4;a7] = [c7 c1;-c1 c7][b4;b7]	transpose [a4;a7] = [c7 -c1;c1 c7][b4;b7]
	M3: [a5;a6] = [c3 c5;-c5 c3][b5;b6]	transpose [a5;a6] = [c3 -c5;c5 c3][b5;b6]
	*/
	a[0]=b[0];//step 4
	a[1]=b[1];
	a[2]=b[2];
	a[3]=b[3];
	a[4]=_mm256_sub_ps(_mm256_mul_ps(b[4], c7), _mm256_mul_ps(b[7], c1));//M2 (transposed)
	a[5]=_mm256_sub_ps(_mm256_mul_ps(b[5], c3), _mm256_mul_ps(b[6], c5));//M3 (transposed)
	a[6]=_mm256_add_ps(_mm256_mul_ps(b[5], c5), _mm256_mul_ps(b[6], c3));
	a[7]=_mm256_add_ps(_mm256_mul_ps(b[4], c1), _mm256_mul_ps(b[7], c7));

	a[0]=_mm256_mul_ps(a[0], _mm256_set1_ps(0.5f));
	a[1]=_mm256_mul_ps(a[1], c4);
	
	b[0]=_mm256_add_ps(a[0], a[1]);//step 3
	b[1]=_mm256_sub_ps(a[0], a[1]);
	b[2]=_mm256_sub_ps(_mm256_mul_ps(a[2], c6), _mm256_mul_ps(a[3], c2));//M1
	b[3]=_mm256_add_ps(_mm256_mul_ps(a[2], c2), _mm256_mul_ps(a[3], c6));
	b[4]=_mm256_add_ps(a[4], a[5]);
	b[5]=_mm256_sub_ps(a[4], a[5]);
	b[6]=_mm256_sub_ps(a[7], a[6]);
	b[7]=_mm256_add_ps(a[7], a[6]);
	
	a[0]=_mm256_add_ps(b[0], b[3]);//step 2
	a[1]=_mm256_add_ps(b[1], b[2]);
	a[2]=_mm256_sub_ps(b[1], b[2]);
	a[3]=_mm256_sub_ps(b[0], b[3]);
	a[4]=b[4];
	a[5]=_mm256_mul_ps(_mm256_sub_ps(b[6], b[5]), c4);//M0
	a[6]=_mm256_mul_ps(_mm256_add_ps(b[6], b[5]), c4);
	a[7]=b[7];
	
	b[0]=_mm256_add_ps(a[0], a[7]);//step 1
	b[1]=_mm256_add_ps(a[1], a[6]);
	b[2]=_mm256_add_ps(a[2], a[5]);
	b[3]=_mm256_add_ps(a[3], a[4]);
	b[4]=_mm256_sub_ps(a[3], a[4]);
	b[5]=_mm256_sub_ps(a[2], a[5]);
	b[6]=_mm256_sub_ps(a[1], a[6]);
	b[7]=_mm256_sub_ps(a[0], a[7]);

	_mm256_store_ps(block+0*8, b[0]);
	_mm256_store_ps(block+1*8, b[1]);
	_mm256_store_ps(block+2*8, b[2]);
	_mm256_store_ps(block+3*8, b[3]);
	_mm256_store_ps(block+4*8, b[4]);
	_mm256_store_ps(block+5*8, b[5]);
	_mm256_store_ps(block+6*8, b[6]);
	_mm256_store_ps(block+7*8, b[7]);
}
AWM_INLINE void transpose(float *block)
{
	__m128 a[4], b[4];
	
	a[0]=_mm_load_ps(block+0*8);
	a[1]=_mm_load_ps(block+1*8);
	a[2]=_mm_load_ps(block+2*8);
	a[3]=_mm_load_ps(block+3*8);
	b[0]=_mm_load_ps(block+4*8+4);
	b[1]=_mm_load_ps(block+5*8+4);
	b[2]=_mm_load_ps(block+6*8+4);
	b[3]=_mm_load_ps(block+7*8+4);
	
	_MM_TRANSPOSE4_PS(a[0], a[1], a[2], a[3]);
	_MM_TRANSPOSE4_PS(b[0], b[1], b[2], b[3]);

	_mm_store_ps(block+0*8, a[0]);
	_mm_store_ps(block+1*8, a[1]);
	_mm_store_ps(block+2*8, a[2]);
	_mm_store_ps(block+3*8, a[3]);
	_mm_store_ps(block+4*8+4, b[0]);
	_mm_store_ps(block+5*8+4, b[1]);
	_mm_store_ps(block+6*8+4, b[2]);
	_mm_store_ps(block+7*8+4, b[3]);
	
	a[0]=_mm_load_ps(block+0*8+4);
	a[1]=_mm_load_ps(block+1*8+4);
	a[2]=_mm_load_ps(block+2*8+4);
	a[3]=_mm_load_ps(block+3*8+4);
	b[0]=_mm_load_ps(block+4*8);
	b[1]=_mm_load_ps(block+5*8);
	b[2]=_mm_load_ps(block+6*8);
	b[3]=_mm_load_ps(block+7*8);
	
	_MM_TRANSPOSE4_PS(a[0], a[1], a[2], a[3]);
	_MM_TRANSPOSE4_PS(b[0], b[1], b[2], b[3]);

	_mm_store_ps(block+0*8+4, b[0]);
	_mm_store_ps(block+1*8+4, b[1]);
	_mm_store_ps(block+2*8+4, b[2]);
	_mm_store_ps(block+3*8+4, b[3]);
	_mm_store_ps(block+4*8, a[0]);
	_mm_store_ps(block+5*8, a[1]);
	_mm_store_ps(block+6*8, a[2]);
	_mm_store_ps(block+7*8, a[3]);
}
AWM_INLINE void gain(float *block, float *g)
{
	__m256 a[8];

	a[0]=_mm256_load_ps(block+0*8);
	a[1]=_mm256_load_ps(block+1*8);
	a[2]=_mm256_load_ps(block+2*8);
	a[3]=_mm256_load_ps(block+3*8);
	a[4]=_mm256_load_ps(block+4*8);
	a[5]=_mm256_load_ps(block+5*8);
	a[6]=_mm256_load_ps(block+6*8);
	a[7]=_mm256_load_ps(block+7*8);

	a[0]=_mm256_mul_ps(a[0], _mm256_load_ps(g+0*8));
	a[1]=_mm256_mul_ps(a[1], _mm256_load_ps(g+1*8));
	a[2]=_mm256_mul_ps(a[2], _mm256_load_ps(g+2*8));
	a[3]=_mm256_mul_ps(a[3], _mm256_load_ps(g+3*8));
	a[4]=_mm256_mul_ps(a[4], _mm256_load_ps(g+4*8));
	a[5]=_mm256_mul_ps(a[5], _mm256_load_ps(g+5*8));
	a[6]=_mm256_mul_ps(a[6], _mm256_load_ps(g+6*8));
	a[7]=_mm256_mul_ps(a[7], _mm256_load_ps(g+7*8));

	_mm256_store_ps(block+0*8, a[0]);
	_mm256_store_ps(block+1*8, a[1]);
	_mm256_store_ps(block+2*8, a[2]);
	_mm256_store_ps(block+3*8, a[3]);
	_mm256_store_ps(block+4*8, a[4]);
	_mm256_store_ps(block+5*8, a[5]);
	_mm256_store_ps(block+6*8, a[6]);
	_mm256_store_ps(block+7*8, a[7]);
}
#else
AWM_INLINE void cvti2f(float *block)
{
	__m128i a[4];

	a[0]=_mm_load_si128((__m128i*)block+0);
	a[1]=_mm_load_si128((__m128i*)block+1);
	a[2]=_mm_load_si128((__m128i*)block+2);
	a[3]=_mm_load_si128((__m128i*)block+3);

	_mm_store_ps(block+0*4, _mm_cvtepi32_ps(a[0]));
	_mm_store_ps(block+1*4, _mm_cvtepi32_ps(a[1]));
	_mm_store_ps(block+2*4, _mm_cvtepi32_ps(a[2]));
	_mm_store_ps(block+3*4, _mm_cvtepi32_ps(a[3]));
}
AWM_INLINE void cvtf2i(float *block)
{
	__m128 a[4];

	a[0]=_mm_load_ps(block+0*4);
	a[1]=_mm_load_ps(block+1*4);
	a[2]=_mm_load_ps(block+2*4);
	a[3]=_mm_load_ps(block+3*4);

	_mm_store_si128((__m128i*)block+0, _mm_cvtps_epi32(a[0]));
	_mm_store_si128((__m128i*)block+1, _mm_cvtps_epi32(a[1]));
	_mm_store_si128((__m128i*)block+2, _mm_cvtps_epi32(a[2]));
	_mm_store_si128((__m128i*)block+3, _mm_cvtps_epi32(a[3]));
}
AWM_INLINE void rgb2yuv(float *c0, float *c1, float *c2)
{
	__m128 half=_mm_set1_ps(128);
	for(int k=0;k<4;++k)
	{
		__m128 rgb[3], yuv[3];

		rgb[0]=_mm_load_ps(c0+k*4);
		rgb[1]=_mm_load_ps(c1+k*4);
		rgb[2]=_mm_load_ps(c2+k*4);
		
		//https://en.wikipedia.org/wiki/YCbCr
		yuv[0]=_mm_mul_ps(_mm_set1_ps(+0.299000f), rgb[0]);
		yuv[1]=_mm_mul_ps(_mm_set1_ps(-0.168736f), rgb[0]);
		yuv[2]=_mm_mul_ps(_mm_set1_ps(+0.500000f), rgb[0]);
		yuv[0]=_mm_add_ps(yuv[0], _mm_mul_ps(_mm_set1_ps(+0.587000f), rgb[1]));
		yuv[1]=_mm_add_ps(yuv[1], _mm_mul_ps(_mm_set1_ps(-0.331264f), rgb[1]));
		yuv[2]=_mm_add_ps(yuv[2], _mm_mul_ps(_mm_set1_ps(-0.418688f), rgb[1]));
		yuv[0]=_mm_add_ps(yuv[0], _mm_mul_ps(_mm_set1_ps(+0.114000f), rgb[2]));
		yuv[1]=_mm_add_ps(yuv[1], _mm_mul_ps(_mm_set1_ps(+0.500000f), rgb[2]));
		yuv[2]=_mm_add_ps(yuv[2], _mm_mul_ps(_mm_set1_ps(-0.081312f), rgb[2]));
		
		yuv[0]=_mm_sub_ps(yuv[0], half);
	//	yuv[1]=_mm_add_ps(yuv[1], half);
	//	yuv[2]=_mm_add_ps(yuv[2], half);

		_mm_store_ps(c0+k*4, yuv[0]);
		_mm_store_ps(c1+k*4, yuv[1]);
		_mm_store_ps(c2+k*4, yuv[2]);
	}
}
AWM_INLINE void yuv2rgb(float *c0, float *c1, float *c2)
{
	__m128 half=_mm_set1_ps(128);
	__m128 vmin=_mm_setzero_ps();
	__m128 vmax=_mm_set1_ps(255);
	for(int k=0;k<4;++k)
	{
		__m128 rgb[3], yuv[3];

		yuv[0]=_mm_load_ps(c0+k*4);
		yuv[1]=_mm_load_ps(c1+k*4);
		yuv[2]=_mm_load_ps(c2+k*4);
		
		yuv[0]=_mm_add_ps(yuv[0], half);
	//	yuv[1]=_mm_sub_ps(yuv[1], half);
	//	yuv[2]=_mm_sub_ps(yuv[2], half);

		rgb[1]=_mm_add_ps(yuv[0], _mm_mul_ps(_mm_set1_ps(-0.344136f), yuv[1]));
		rgb[2]=_mm_add_ps(yuv[0], _mm_mul_ps(_mm_set1_ps(+1.772f), yuv[1]));
		rgb[0]=_mm_add_ps(yuv[0], _mm_mul_ps(_mm_set1_ps(+1.402f), yuv[2]));
		rgb[1]=_mm_add_ps(rgb[1], _mm_mul_ps(_mm_set1_ps(-0.714136f), yuv[2]));

		rgb[0]=_mm_min_ps(rgb[0], vmax);
		rgb[1]=_mm_min_ps(rgb[1], vmax);
		rgb[2]=_mm_min_ps(rgb[2], vmax);
		rgb[0]=_mm_max_ps(rgb[0], vmin);
		rgb[1]=_mm_max_ps(rgb[1], vmin);
		rgb[2]=_mm_max_ps(rgb[2], vmin);
		_mm_store_ps(c0+k*4, rgb[0]);
		_mm_store_ps(c1+k*4, rgb[1]);
		_mm_store_ps(c2+k*4, rgb[2]);
	}
}
AWM_INLINE void dct_fwd(float *block)//DCT-II
{
	__m128 c3=_mm_set1_ps(0.3826834323650897717284599840304f);//cosd(90*3/4) = sind(90*1/4) = s1
	__m128 c2=_mm_set1_ps(0.7071067811865475244008443621048f);//cosd(90*2/4) = sind(90*2/4) = s2 = sqrt(2)
	__m128 c1=_mm_set1_ps(0.9238795325112867561281831893968f);//cosd(90*1/4) = sind(90*3/4) = s3

	__m128 a[4], b[4];

	a[0]=_mm_load_ps(block+0*4);
	a[1]=_mm_load_ps(block+1*4);
	a[2]=_mm_load_ps(block+2*4);
	a[3]=_mm_load_ps(block+3*4);

	b[0]=_mm_add_ps(a[0], a[3]);//step 1
	b[1]=_mm_add_ps(a[1], a[2]);
	b[2]=_mm_sub_ps(a[1], a[2]);
	b[3]=_mm_sub_ps(a[0], a[3]);

	a[0]=_mm_add_ps(b[0], b[1]);//step 2
	a[1]=_mm_mul_ps(_mm_sub_ps(b[0], b[1]), c2);
	a[2]=_mm_add_ps(_mm_mul_ps(b[3], c1), _mm_mul_ps(b[2], c3));
	a[3]=_mm_sub_ps(_mm_mul_ps(b[3], c3), _mm_mul_ps(b[2], c1));

	_mm_store_ps(block+0*4, a[0]);
	_mm_store_ps(block+2*4, a[1]);
	_mm_store_ps(block+1*4, a[2]);
	_mm_store_ps(block+3*4, a[3]);
}
AWM_INLINE void dct_inv(float *block)//DCT-III
{
	__m128 half=_mm_set1_ps(0.5f);
	__m128 c3=_mm_set1_ps(0.3826834323650897717284599840304f);//cosd(90*3/4) = sind(90*1/4) = s1
	__m128 c2=_mm_set1_ps(0.7071067811865475244008443621048f);//cosd(90*2/4) = sind(90*2/4) = s2 = sqrt(2)
	__m128 c1=_mm_set1_ps(0.9238795325112867561281831893968f);//cosd(90*1/4) = sind(90*3/4) = s3

	__m128 a[4], b[4];

	a[0]=_mm_load_ps(block+0*4);
	a[1]=_mm_load_ps(block+2*4);
	a[2]=_mm_load_ps(block+1*4);
	a[3]=_mm_load_ps(block+3*4);
	
	/*
	matrix					inverse
	M0 = [1 1;1 -1]				same/2
	M1: [a0;a1] = [1 1;c2 -c2][b0;b1]	[b0;b1] = [1/2 c2;1/2 -c2][a0;a1]
	M2: [a2;a3] = [c3 c1;-c1 c3][b2;b3]	transpose [b2;b3] = [c3 -c1;c1 c3][a2;a3]
	*/
	a[0]=_mm_mul_ps(a[0], half);
	a[1]=_mm_mul_ps(a[1], c2);

	b[0]=_mm_add_ps(a[0], a[1]);//step 2
	b[1]=_mm_sub_ps(a[0], a[1]);
	b[2]=_mm_sub_ps(_mm_mul_ps(a[2], c3), _mm_mul_ps(a[3], c1));
	b[3]=_mm_add_ps(_mm_mul_ps(a[2], c1), _mm_mul_ps(a[3], c3));

	a[0]=_mm_add_ps(b[0], b[3]);//step 1
	a[1]=_mm_add_ps(b[1], b[2]);
	a[2]=_mm_sub_ps(b[1], b[2]);
	a[3]=_mm_sub_ps(b[0], b[3]);

	_mm_store_ps(block+0*4, a[0]);
	_mm_store_ps(block+1*4, a[1]);
	_mm_store_ps(block+2*4, a[2]);
	_mm_store_ps(block+3*4, a[3]);
}
AWM_INLINE void transpose(float *block)
{
	__m128 a[4];
	
	a[0]=_mm_load_ps(block+0*4);
	a[1]=_mm_load_ps(block+1*4);
	a[2]=_mm_load_ps(block+2*4);
	a[3]=_mm_load_ps(block+3*4);
	
	_MM_TRANSPOSE4_PS(a[0], a[1], a[2], a[3]);

	_mm_store_ps(block+0*4, a[0]);
	_mm_store_ps(block+1*4, a[1]);
	_mm_store_ps(block+2*4, a[2]);
	_mm_store_ps(block+3*4, a[3]);
}
AWM_INLINE void gain(float *block, float *g)
{
	__m128 a[4];

	a[0]=_mm_load_ps(block+0*4);
	a[1]=_mm_load_ps(block+1*4);
	a[2]=_mm_load_ps(block+2*4);
	a[3]=_mm_load_ps(block+3*4);

	a[0]=_mm_mul_ps(a[0], _mm_load_ps(g+0*4));
	a[1]=_mm_mul_ps(a[1], _mm_load_ps(g+1*4));
	a[2]=_mm_mul_ps(a[2], _mm_load_ps(g+2*4));
	a[3]=_mm_mul_ps(a[3], _mm_load_ps(g+3*4));

	_mm_store_ps(block+0*4, a[0]);
	_mm_store_ps(block+1*4, a[1]);
	_mm_store_ps(block+2*4, a[2]);
	_mm_store_ps(block+3*4, a[3]);
}
#endif
AWM_INLINE void dc_predict(int16_t **rows, int32_t *weights, int32_t *estim, int32_t *ret_p1, int32_t *ret_pred, int32_t *ret_ctx)
{
	int
		NNN	=rows[3][0+0*NCH*NROWS*NVAL],
		NN	=rows[2][0+0*NCH*NROWS*NVAL],
		NNE	=rows[2][0+1*NCH*NROWS*NVAL],
		NW	=rows[1][0-1*NCH*NROWS*NVAL],
		N	=rows[1][0+0*NCH*NROWS*NVAL],
		NE	=rows[1][0+1*NCH*NROWS*NVAL],
		NEE	=rows[1][0+2*NCH*NROWS*NVAL],
		NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
		NEEEE	=rows[1][0+4*NCH*NROWS*NVAL],
		WWWW	=rows[0][0-4*NCH*NROWS*NVAL],
		WWW	=rows[0][0-3*NCH*NROWS*NVAL],
		WW	=rows[0][0-2*NCH*NROWS*NVAL],
		W	=rows[0][0-1*NCH*NROWS*NVAL],
		eW	=rows[0][1-1*NCH*NROWS*NVAL];
	int ctx=FLOOR_LOG2(eW*eW+1);
	int pred=1<<SHIFT>>1, j=0, p1;
	int vmax=N, vmin=W;

	if(ctx>NCTX_DC-1)
		ctx=NCTX_DC-1;
	if(N<W)vmin=N, vmax=W;
#ifdef USE_L1
	if(vmin>NE)vmin=NE;
	if(vmax<NE)vmax=NE;
	if(vmin>NEEE)vmin=NEEE;
	if(vmax<NEEE)vmax=NEEE;
#define PRED(E) estim[j]=E; pred+=weights[j]*estim[j]; ++j;
	j=0;
	PREDLIST
#undef  PRED
	pred>>=SHIFT;
	p1=pred;
#else
	pred=N=W-NW;
#endif
	CLAMP2(pred, vmin, vmax);
	*ret_p1=p1;
	*ret_pred=pred;
	*ret_ctx=ctx;
}
AWM_INLINE void dc_update(int16_t **rows, int32_t *weights, int32_t *estim, int32_t p1, int32_t curr, int32_t sym)
{
	int
		eNEEE	=rows[1][1+3*NCH*NROWS*NVAL],
		eW	=rows[0][1-1*NCH*NROWS*NVAL];
#ifdef USE_L1
	int e=(curr>p1)-(curr<p1), j=0;
#define PRED(...) weights[j]+=e*estim[j]; ++j;
	j=0;
	PREDLIST
#undef  PRED
#endif
	rows[0][1]=(2*eW+(sym<<GRBITS)+eNEEE)>>2;
}

#if 0
#define TEST
static void dcttest()
{
	ALIGN(32) float data[BLOCKX*BLOCKY]=
	{
		1,  3,  3,  1,
		3,  3,  6,  3,
		3,  6, 12,  6,
		1,  3,  6,  3,
	};
	ALIGN(32) float d0[BLOCKX*BLOCKY], diff[BLOCKX*BLOCKY];
	for(int k=0;k<BLOCKX*BLOCKY;++k)
		//data[k]=(float)rand();
		data[k]=(int)k;

	memcpy(d0, data, sizeof(d0));

	dct_fwd(data);
	transpose(data);
	dct_fwd(data);

	dct_inv(data);
	transpose(data);
	dct_inv(data);
	
	for(int k=0;k<BLOCKX*BLOCKY;++k)
		diff[k]=d0[k]-data[k]/16;
	for(int k=0;k<BLOCKX*BLOCKY;++k)
	{
		if(fabsf(diff[k])>1e-2)
			CRASH("");
	}

	printf("SUCCESS\n");
	exit(0);
}
#endif

static void applydct(float *buf, int iw, int ih)
{
	ALIGN(32) float block[BLOCKX*BLOCKY];
	for(int ky=0;ky<ih;ky+=BLOCKY)
	{
		for(int kx=0;kx<iw;kx+=BLOCKX)
		{
			for(int ky2=0;ky2<BLOCKY;++ky2)
				memcpy(block+BLOCKX*ky2, buf+iw*(ky+ky2)+kx, sizeof(float[BLOCKX]));
		}
		dct_fwd(block);
		transpose(block);
		dct_fwd(block);
		for(int kx=0;kx<iw;kx+=BLOCKX)
		{
			for(int ky2=0;ky2<BLOCKY;++ky2)
				memcpy(buf+iw*(ky+ky2)+kx, block+BLOCKX*ky2, sizeof(float[BLOCKX]));
		}
	}
}
static void getqtable(float *buf, int iw, int ih, float *qtable)
{
	memset(qtable, 0, sizeof(float[BLOCKX*BLOCKY]));
	for(int ky=0;ky<ih;ky+=BLOCKY)
	{
		for(int kx=0;kx<iw;kx+=BLOCKX)
		{
			for(int ky2=0;ky2<BLOCKY;++ky2)
			{
				for(int kx2=0;kx2<BLOCKY;++kx2)
					qtable[BLOCKX*ky2+kx2]+=fabsf(buf[iw*(ky+ky2)+kx+kx2]);
			}
		}
	}
}
static void printqtable(float *qtable, const char *label)
{
	printf("%s\n", label);
	for(int ky=0;ky<BLOCKY;++ky)
	{
		for(int kx=0;kx<BLOCKY;++kx)
			printf(" %12.2lf", qtable[BLOCKX*ky+kx]);
		printf("\n");
	}
	printf("\n");
}


#ifdef PROFILE_SIZE
static double csizes[3][BLOCKX*BLOCKY]={0};
#endif

static uint32_t hweight[3*NCTX];
static uint32_t hists[3*NCTX*NLEVELS];
ALIGN(32) static float blocks[3][BLOCKX*BLOCKY]={0}, qtable1[BLOCKX*BLOCKY]={0}, qtable2[BLOCKX*BLOCKY]={0};
int c51_codec(int argc, char **argv)
{
#ifdef TEST
	dcttest();//
#endif

	const uint16_t tag='5'|'1'<<8;

	const char *srcfn=0, *dstfn=0;
	FILE *fsrc=0, *fdst=0;
	uint64_t c=0;
	int fwd=0, iw=0, ih=0;
	int64_t usize=0, ccap=0, csize=0;
	int psize=0;
	int16_t *pixels=0;
	uint8_t *image=0, *imptr=0, *stream=0, *streamptr=0;
#ifdef _MSC_VER
	uint8_t *streamend=0;
#endif
	int32_t weights[NCH][NPREDS]={0}, estim[NPREDS]={0};
	int32_t pred=0, ctx=0, error=0, sym=0, curr=0, cdf=0, freq=0;
	uint64_t low=0, range=0xFFFFFFFFFFFF, code=0;
#ifdef LOUD
	double t=0;
#endif

	if(argc!=3)
	{
		printf(
			"Usage:  \"%s\"  src  dst\n"
			"Only for 24-bit PPM images\n"
			"Built on %s %s\n"
			, argv[0]
			, __DATE__, __TIME__
		);
		return 1;
	}
#ifdef LOUD
	t=time_sec();
#endif
	srcfn=argv[1];
	dstfn=argv[2];
	
	fsrc=fopen(srcfn, "rb");
	if(!fsrc)
	{
		CRASH("Cannot open \"%s\"", srcfn);
		return 1;
	}
	fread(&c, 1, 2, fsrc);
	fwd=c==('P'|'6'<<8);
	if(!fwd&&c!=tag)
	{
		CRASH("Unsupported file \"%s\"", srcfn);
		return 1;
	}
	if(fwd)//parse header
	{
		c=fgetc(fsrc);
		if(c!='\n')
		{
			CRASH("Unsupported PPM file");
			return 1;
		}
		c=fgetc(fsrc);
		while(c=='#')
		{
			c=fgetc(fsrc);
			while(c!='\n')
				c=fgetc(fsrc);
			c=fgetc(fsrc);
		}
		iw=0;
		while((uint32_t)(c-'0')<10)
		{
			iw=10*iw+(int32_t)c-'0';
			c=fgetc(fsrc);
		}
		while(c<=' ')
			c=fgetc(fsrc);
		ih=0;
		while((uint32_t)(c-'0')<10)
		{
			ih=10*ih+(int32_t)c-'0';
			c=fgetc(fsrc);
		}
		while(c=='#')
		{
			c=fgetc(fsrc);
			while(c!='\n')
				c=fgetc(fsrc);
			c=fgetc(fsrc);
		}
		c|=(int64_t)fgetc(fsrc)<<8*1;
		c|=(int64_t)fgetc(fsrc)<<8*2;
		c|=(int64_t)fgetc(fsrc)<<8*3;
		c|=(int64_t)fgetc(fsrc)<<8*4;
		if(c!=(
			(uint64_t)'\n'<<8*0|
			(uint64_t) '2'<<8*1|
			(uint64_t) '5'<<8*2|
			(uint64_t) '5'<<8*3|
			(uint64_t)'\n'<<8*4
		))
		{
			CRASH("Unsupported PPM file");
			return 1;
		}
		ccap=(int64_t)4*iw*ih;
	}
	else
	{
		iw=0;
		ih=0;
		fread(&iw, 1, 3, fsrc);
		fread(&ih, 1, 3, fsrc);
		{
			struct stat info={0};

			stat(srcfn, &info);
			ccap=(int64_t)info.st_size-ftell(fsrc);
		}
	}
	if(iw<1||ih<1)
	{
		CRASH("Unsupported source file");
		return 1;
	}
	usize=(int64_t)3*iw*ih;
	image=(uint8_t*)malloc(usize);
	stream=(uint8_t*)malloc(ccap);
	psize=(iw/BLOCKX+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
	pixels=(int16_t*)malloc(psize);
	if(!image||!stream||!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	streamptr=stream;
#ifdef _MSC_VER
	streamend=stream+ccap;
#endif
	if(fwd)
	{
		fread(image, 1, usize, fsrc);
#ifdef ENABLE_GUIDE
		guide_save(image, iw, ih);
#endif
	}
	else
	{
		c=fread(stream, 1, ccap, fsrc);
	}
	fclose(fsrc);

	ALIGN(32) float qtables[3][BLOCKX*BLOCKY]={0};
	if(fwd)
	{
		int planew=(iw+2*BLOCKX-1)&~(2*BLOCKX-1), planeh=(ih+2*BLOCKY-1)&~(2*BLOCKY-1);
		int cplanew=planew/2, cplaneh=planeh/2;
		int64_t planesize=sizeof(float[3])*planew*planeh;
		float *planes=(float*)malloc(planesize);
		if(!planes)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(planes, 0, planesize);
		int64_t planeres=(int64_t)planew*planeh;
		float *ybuf=planes, *ubuf=ybuf+planeres, *vbuf=ubuf+planeres;
		uint8_t *imptr=image;
		float *yptr=ybuf, *uptr=ubuf, *vptr=vbuf;
		for(int ky=0;ky<ih;++ky)//rgb to yuv
		{
			for(int kx=0;kx<iw;++kx, imptr+=3)
			{
				float r=imptr[0], g=imptr[1], b=imptr[2];
				*yptr++=0.299f*r+0.587f*g+0.114f*b-128;
				*uptr++=-0.168736f*r-0.331264f*g+0.5f*b;
				*vptr++=0.5f*r-0.418688f*g-0.081312f*b;
			}
		}
		//yptr=ybuf;
		//uptr=ubuf;
		//vptr=vbuf;
		for(int kx=iw;kx<planew;++kx)//right remainder extension
		{
			for(int ky=0;ky<ih;++ky)
			{
				ybuf[planew*ky+kx]=ybuf[planew*ky+kx-1];
				ubuf[planew*ky+kx]=ubuf[planew*ky+kx-1];
				vbuf[planew*ky+kx]=vbuf[planew*ky+kx-1];
			}
		}
		for(int ky=ih;ky<planeh;++ky)//bottom remainder extension
		{
			memcpy(ybuf+planew*ky, ybuf+planew*(ky-1), sizeof(float)*planew);
			memcpy(ubuf+planew*ky, ubuf+planew*(ky-1), sizeof(float)*planew);
			memcpy(vbuf+planew*ky, vbuf+planew*(ky-1), sizeof(float)*planew);
		}
		uptr=ubuf;
		vptr=vbuf;
		for(int ky=0;ky<ih;ky+=2)//chroma subsampling
		{
			for(int kx=0;kx<iw;kx+=2)
			{
				*uptr++=0.25f*(
					+ubuf[planew*(ky+0)+(kx+0)]
					+ubuf[planew*(ky+0)+(kx+1)]
					+ubuf[planew*(ky+1)+(kx+0)]
					+ubuf[planew*(ky+1)+(kx+1)]
				);
				*vptr++=0.25f*(
					+vbuf[planew*(ky+0)+(kx+0)]
					+vbuf[planew*(ky+0)+(kx+1)]
					+vbuf[planew*(ky+1)+(kx+0)]
					+vbuf[planew*(ky+1)+(kx+1)]
				);
			}
		}
		applydct(ybuf, planew, planeh);
		applydct(ubuf, cplanew, cplaneh);
		applydct(vbuf, cplanew, cplaneh);

		getqtable(ybuf, planew, planeh, qtables[0]);
		getqtable(ubuf, cplanew, cplaneh, qtables[1]);
		getqtable(vbuf, cplanew, cplaneh, qtables[2]);
		float gain=(float)BLOCKX*BLOCKY/(planew*planeh);
		for(int k=0;k<BLOCKX*BLOCKY;++k)
			qtables[0][k]*=gain;
		gain=(float)BLOCKX*BLOCKY/(cplanew*cplaneh);
		for(int k=0;k<BLOCKX*BLOCKY;++k)
		{
			qtables[1][k]*=gain;
			qtables[2][k]*=gain;
		}
#ifdef LOUD
		printqtable(qtables[0], "Y");
		printqtable(qtables[1], "U");
		printqtable(qtables[2], "V");
#endif

#if 0
		for(int y1=0;y1<ih;y1+=2*BLOCKY)
		{
			int y2=y1+2*BLOCKY;
			if(y2>ih)
				y2=ih;
			int dy=y2-y1;

			for(int x1=0;x1<iw;x1+=2*BLOCKX)
			{
				int x2=x1+2*BLOCKY;
				if(x2>iw)
					x2=iw;
				int dx=x2-x1;

				for(int ky2=0;ky2<dy;ky2+=2)
				{
					imptr=image+3*(iw*(y1+ky2)+x1);
					for(int kx2=0;kx2<dx;kx2+=2, imptr+=3)
					{
						((int*)blocks[0])[BLOCKX*ky2+kx2]=imptr[0];
						((int*)blocks[1])[BLOCKX*ky2+kx2]=imptr[1];
						((int*)blocks[2])[BLOCKX*ky2+kx2]=imptr[2];
					}
					for(int kx2=dx;kx2<BLOCKX;++kx2)
					{
						blocks[0][BLOCKX*ky2+kx2]=blocks[0][BLOCKX*ky2+dx-1];
						blocks[1][BLOCKX*ky2+kx2]=blocks[1][BLOCKX*ky2+dx-1];
						blocks[2][BLOCKX*ky2+kx2]=blocks[2][BLOCKX*ky2+dx-1];
					}
				}
				for(int ky2=dy;ky2<BLOCKY;++ky2)
				{
					memcpy(blocks[0]+BLOCKX*ky2, blocks[0]+BLOCKX*(dy-1), sizeof(float[BLOCKX]));
					memcpy(blocks[1]+BLOCKX*ky2, blocks[1]+BLOCKX*(dy-1), sizeof(float[BLOCKX]));
					memcpy(blocks[2]+BLOCKX*ky2, blocks[2]+BLOCKX*(dy-1), sizeof(float[BLOCKX]));
				}
			}
		}
#endif
		free(planes);
		exit(0);
	}
	else
	{
	}
#if 0
	for(int k=0;k<NCH*NPREDS;++k)
		((int32_t*)weights)[k]=(1<<SHIFT)/NPREDS;
	memset(pixels, 0, psize);
	memset(hists, 0, sizeof(hists));
	memset(hweight, 0, sizeof(hweight));
	if(fwd)
	{
		for(int k=0;k<BLOCKX*BLOCKY;++k)
		{
			qtable1[k]=1/(DCT_ROUND_TRIP_GAIN*qtable_luma[k]);
			qtable2[k]=1/(DCT_ROUND_TRIP_GAIN*qtable_chroma[k]);
		}
	}
	else
	{
		for(int k=0;k<BLOCKX*BLOCKY;++k)
		{
			qtable1[k]=qtable_luma[k];
			qtable2[k]=qtable_chroma[k];
		}
		code=0;
		code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);//load
		code=code<<32|*(uint32_t*)streamptr; streamptr+=sizeof(uint32_t);
	}
	for(int ky=0;ky<ih;ky+=BLOCKY)
	{
		int16_t *rows[]=
		{
			pixels+(XPAD*NCH*NROWS+(ky/BLOCKY-0LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky/BLOCKY-1LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky/BLOCKY-2LL+NROWS)%NROWS)*NVAL,
			pixels+(XPAD*NCH*NROWS+(ky/BLOCKY-3LL+NROWS)%NROWS)*NVAL,
		};
		int y2=ky+BLOCKY, dy;

		if(y2>ih)
			y2=ih;
		dy=y2-ky;
		for(int kx=0;kx<iw;kx+=BLOCKX)
		{
			int x2=kx+BLOCKX, dx;
			int den=0;

			if(x2>iw)
				x2=iw;
			dx=x2-kx;
			if(fwd)
			{
				for(int ky2=0;ky2<dy;++ky2)
				{
					imptr=image+3*(iw*(ky+ky2)+kx);
					for(int kx2=0;kx2<dx;++kx2, imptr+=3)
					{
						((int*)blocks[0])[BLOCKX*ky2+kx2]=imptr[0];
						((int*)blocks[1])[BLOCKX*ky2+kx2]=imptr[1];
						((int*)blocks[2])[BLOCKX*ky2+kx2]=imptr[2];
					}
					for(int kx2=dx;kx2<BLOCKX;++kx2)
					{
						blocks[0][BLOCKX*ky2+kx2]=blocks[0][BLOCKX*ky2+dx-1];
						blocks[1][BLOCKX*ky2+kx2]=blocks[1][BLOCKX*ky2+dx-1];
						blocks[2][BLOCKX*ky2+kx2]=blocks[2][BLOCKX*ky2+dx-1];
					}
				}
				for(int ky2=dy;ky2<BLOCKY;++ky2)
				{
					memcpy(blocks[0]+BLOCKX*ky2, blocks[0]+BLOCKX*(dy-1), sizeof(float[BLOCKX]));
					memcpy(blocks[1]+BLOCKX*ky2, blocks[1]+BLOCKX*(dy-1), sizeof(float[BLOCKX]));
					memcpy(blocks[2]+BLOCKX*ky2, blocks[2]+BLOCKX*(dy-1), sizeof(float[BLOCKX]));
				}
				cvti2f(blocks[0]);
				cvti2f(blocks[1]);
				cvti2f(blocks[2]);
				rgb2yuv(blocks[0], blocks[1], blocks[2]);
				dct_fwd(blocks[0]);
				dct_fwd(blocks[1]);
				dct_fwd(blocks[2]);
				transpose(blocks[0]);
				transpose(blocks[1]);
				transpose(blocks[2]);
				dct_fwd(blocks[0]);
				dct_fwd(blocks[1]);
				dct_fwd(blocks[2]);
#ifdef _MSC_VER
				for(int kc=0;kc<3;++kc)for(int k=0;k<BLOCKX*BLOCKY;++k)sum_before[kc][k]+=fabsf(blocks[kc][k]);
#endif
				gain(blocks[0], qtable1);
				gain(blocks[1], qtable2);
				gain(blocks[2], qtable2);
				cvtf2i(blocks[0]);
				cvtf2i(blocks[1]);
				cvtf2i(blocks[2]);
#ifdef _MSC_VER
				for(int kc=0;kc<3;++kc)
				{
					for(int k=0;k<BLOCKX*BLOCKY;++k)
					{
						int val=((int*)blocks[kc])[k];
						sum_after[kc][k]+=abs(val);
						if(cmin[kc][k]>val)cmin[kc][k]=val;
						if(cmax[kc][k]<val)cmax[kc][k]=val;
					}
				}
#endif
			}
			for(int kc=0;kc<3;++kc)
			{
				for(int k=0;k<BLOCKX*BLOCKY;++k)
				{
					int32_t p1=0;
				//	int ctx2;

				//	if(ky==2852&&kx==4288)//
				//	if(ky==2852&&kx==4284)//
				//	if(ky==0&&kx==8)//
				//	if(ky==0&&kx==12)//
				//	if(ky==0&&kx==2804)//
				//	if(ky==0&&kx==8&&!kc&&!k)//
				//	if(ky==2852&&kx==4284&&kc==1&&k==8)//
				//	if(ky==12&&kx==1172&&!kc&&!k)//
				//	if(ky==8&&kx==1676&&!kc&&!k)//
				//	if(ky==0&&kx==0&&kc==1&&!k)//
				//	if(ky==0&&kx==0&&kc==0&&k==1)//
				//	if(ky==0&&kx==4&&kc==0&&k==1)//
				//	if(ky==0&&kx==256&&kc==1&&k==0)//
				//	if(ky==0&&kx==256&&kc==0&&k==0)//
				//	if(ky==0&&kx==256&&kc==0&&k==0)//
				//	if(ky==2652&&kx==3132&&kc==1&&k==7)//
				//	if(ky==0&&kx==0&&kc==0&&k==4)//
				//	if(ky==872&&kx==1544&&kc==0&&k==0)//
				//		printf("");

					if(!k)
					{
						//if(ky==0&&kx==0&&!kc)//
						//	printf("");

						dc_predict(rows, weights[kc], estim, &p1, &pred, &ctx);
					//	ctx2=ctx+(ctx<NCTX_DC-1);
					}
					else
					{
						pred=0;
						ctx=k/BLOCKX+k%BLOCKX-1;
						if(ctx>NCTX_AC-1)
							ctx=NCTX_AC-1;
					//	ctx2=ctx+(ctx<NCTX_AC-1);
						ctx+=NCTX_DC;
					//	ctx2+=NCTX_DC;
					}
					ctx+=NCTX*kc;
				//	ctx2+=NCTX*kc;
					uint32_t *currhist=hists+NLEVELS*ctx;
				//	uint32_t *currhist2=hists+NLEVELS*ctx2;
				//	den=hweight[ctx]+hweight[ctx2]+NLEVELS;
					den=hweight[ctx]+NLEVELS;
					if(fwd)
					{
						int t;
#ifdef FIFOVAL
						valfifo_enqueue(ctx<<24^p1<<12^pred);
#endif
						curr=((int*)blocks[kc])[k];
						CLAMP2(curr, -(NLEVELS>>1), (NLEVELS>>1)-1);
						error=curr-pred;
						if(error>(NLEVELS>>1)-1)
							error-=NLEVELS;
						if(error<-(NLEVELS>>1))
							error+=NLEVELS;
					//	error=((curr-pred+(NLEVELS>>1))&(NLEVELS-1))-(NLEVELS>>1);//POT
						sym=error<<1^error>>31;

						if(range<=0xFFFF)
						{
							*(uint32_t*)streamptr=(uint32_t)(low>>32);
							streamptr+=4;
							low<<=32;
							range=range<<32|0xFFFFFFFF;
							if(range>~low)
								range=~low;
						}
						for(t=0, cdf=0;;++t)
						{
						//	freq=currhist[t]+currhist2[t]+1;
							freq=currhist[t]+1;
							if(t>=sym)
								break;
							cdf+=freq;
						}
#ifdef _MSC_VER
						if(!freq||(uint32_t)freq>(uint32_t)den||(uint32_t)cdf>(uint32_t)den)
							CRASH("");
#endif
#ifdef FIFOVAL
						valfifo_enqueue(freq<<16|cdf);
						valfifo_enqueue(curr);
#endif
						low+=range*cdf/den;
						range=range*freq/den-1;
#ifdef PROFILE_SIZE
						csizes[kc][k]+=-log2((double)freq/den);
#endif
					}
					else
					{
#ifdef FIFOVAL
						valfifo_check(ctx<<24^p1<<12^pred);
#endif
						if(range<=0xFFFF)//stall: unpredictable branch
						{
#ifdef _MSC_VER
							if(streamptr>streamend)
								CRASH("");
#endif
							code=code<<32|*(uint32_t*)streamptr;
							streamptr+=sizeof(uint32_t);
							low<<=32;
							range=range<<32|0xFFFFFFFF;
							if(range>~low)
								range=~low;
						}
						int c2=(int)(((code-low+1)*den-1)/range);
						sym=0;
						cdf=0;
						for(;;)
						{
						//	freq=currhist[sym]+currhist2[sym]+1;
							freq=currhist[sym]+1;
							if(cdf+freq>c2)
								break;
							cdf+=freq;
							++sym;
						}
#ifdef FIFOVAL
						valfifo_check(freq<<16|cdf);
#endif
						low+=range*cdf/den;
						range=range*freq/den-1;

						error=sym>>1^-(sym&1);
						curr=error+pred;
						if(curr>(NLEVELS>>1)-1)
							curr-=NLEVELS;
						if(curr<-(NLEVELS>>1))
							curr+=NLEVELS;
					//	curr=((error+pred+(NLEVELS>>1))&(NLEVELS-1))-(NLEVELS>>1);//POT
#ifdef FIFOVAL
						valfifo_check(curr);
#endif
						((int*)blocks[kc])[k]=curr;
					}
					
					++currhist[sym];
					++hweight[ctx];
					if(hweight[ctx]>=0xFFFF-NLEVELS)
					{
						den=0;
						for(int k=0;k<NLEVELS;++k)
							den+=currhist[k]>>=1;
						hweight[ctx]=den;
					}

					//++currhist2[sym];
					//++hweight[ctx2];
					//if(hweight[ctx2]>=0x5000)
					//{
					//	den=0;
					//	for(int k=0;k<NLEVELS;++k)
					//		den+=currhist2[k]>>=1;
					//	hweight[ctx2]=den;
					//}
					if(!k)
					{
						dc_update(rows, weights[kc], estim, p1, curr, sym);
						rows[0][0]=curr;
						rows[0]+=NROWS*NVAL;
						rows[1]+=NROWS*NVAL;
						rows[2]+=NROWS*NVAL;
						rows[3]+=NROWS*NVAL;
					}
				}
			}
			if(!fwd)
			{
				cvti2f(blocks[0]);
				cvti2f(blocks[1]);
				cvti2f(blocks[2]);
				gain(blocks[0], qtable1);
				gain(blocks[1], qtable2);
				gain(blocks[2], qtable2);
				dct_inv(blocks[0]);
				dct_inv(blocks[1]);
				dct_inv(blocks[2]);
				transpose(blocks[0]);
				transpose(blocks[1]);
				transpose(blocks[2]);
				dct_inv(blocks[0]);
				dct_inv(blocks[1]);
				dct_inv(blocks[2]);
				yuv2rgb(blocks[0], blocks[1], blocks[2]);
				cvtf2i(blocks[0]);
				cvtf2i(blocks[1]);
				cvtf2i(blocks[2]);
				for(int ky2=0;ky2<dy;++ky2)
				{
					imptr=image+3*(iw*(ky+ky2)+kx);
					for(int kx2=0;kx2<dx;++kx2, imptr+=3)
					{
						imptr[0]=((int*)blocks[0])[BLOCKX*ky2+kx2];
						imptr[1]=((int*)blocks[1])[BLOCKX*ky2+kx2];
						imptr[2]=((int*)blocks[2])[BLOCKX*ky2+kx2];
					}
				}
#ifdef ENABLE_GUIDE
				for(int ky2=0;ky2<dy;++ky2)
				{
					imptr=image+3*(iw*(ky+ky2)+kx);
					for(int kx2=0;kx2<dx;++kx2, imptr+=3)
					{
						int diff;

						diff=g_image[imptr-image+0]-imptr[0]; g_sqe[0]+=diff*diff;
						diff=g_image[imptr-image+1]-imptr[1]; g_sqe[1]+=diff*diff;
						diff=g_image[imptr-image+2]-imptr[2]; g_sqe[2]+=diff*diff;
					}
				}
#endif
			}
		}
	}
#endif
	free(pixels);

	fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
	if(fwd)
	{
		*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;//flush
		*(uint32_t*)streamptr=(uint32_t)(low>>32); streamptr+=sizeof(uint32_t); low<<=32;
		
		csize+=fwrite(&tag, 1, 2, fdst);
		csize+=fwrite(&iw, 1, 3, fdst);
		csize+=fwrite(&ih, 1, 3, fdst);
		csize+=fwrite(stream, 1, streamptr-stream, fdst);
	}
	else
	{
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		fwrite(image, 1, usize, fdst);
	}
	fclose(fdst);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
#ifdef PROFILE_SIZE
		double sum=0;
		for(int kc=0;kc<3;++kc)
		{
			for(int k=0;k<BLOCKX*BLOCKY;++k)
				sum+=csizes[kc][k];
		}
		if(!sum)
			sum=1;
		sum=100/sum;
		for(int kc=0;kc<3;++kc)
		{
			double sum2=0;
			for(int k=0;k<BLOCKX*BLOCKY;++k)
				sum2+=csizes[kc][k];
			printf("%c: %12.2lf\n", "YUV"[kc], sum2);
			for(int ky=0;ky<BLOCKY;++ky)
			{
				for(int kx=0;kx<BLOCKX;++kx)
					printf(" %12.2lf %8.4lf%%"
						, csizes[kc][BLOCKX*ky+kx]
						, csizes[kc][BLOCKX*ky+kx]*sum
					);
				printf("\n");
			}
		}
#endif
#ifdef _MSC_VER
		for(int kc=0;kc<3;++kc)
		{
			printf("%c before quant:\n", "YUV"[kc]);
			for(int ky=0;ky<BLOCKY;++ky)
			{
				for(int kx=0;kx<BLOCKY;++kx)
					printf(" %16.0lf", sum_before[kc][BLOCKX*ky+kx]);
				printf("\n");
			}
			printf("\n");

			printf("%c after quant:\n", "YUV"[kc]);
			for(int ky=0;ky<BLOCKY;++ky)
			{
				for(int kx=0;kx<BLOCKY;++kx)
					printf(" %16.0lf", sum_after[kc][BLOCKX*ky+kx]);
				printf("\n");
			}
			printf("\n");

			printf("%c range after quant:\n", "YUV"[kc]);
			for(int ky=0;ky<BLOCKY;++ky)
			{
				for(int kx=0;kx<BLOCKY;++kx)
					printf(" %7d ~ %7d", cmin[kc][BLOCKX*ky+kx], cmax[kc][BLOCKX*ky+kx]);
				printf("\n");
			}
			printf("\n");
		}
#endif
		printf("WH %5d*%5d  \"%s\"\n", iw, ih, srcfn);
		printf("%9td->%9td  %8.4lf%%  %12.6lf:1  BPD %12.6lf\n"
			, usize
			, csize
			, 100.*csize/usize
			, (double)usize/csize
			, 8.*csize/usize
		);
	}
	printf("%c  %12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB\n"
		, 'D'+fwd
		, t
		, usize/(t*1024*1024)
		, t*1024*1024*1000/usize
	);
	if(!fwd)
	{
		double invres=1/((double)iw*ih);
		double rmse[]=
		{
			sqrt(g_sqe[0]*invres),
			sqrt(g_sqe[1]*invres),
			sqrt(g_sqe[2]*invres),
			0,
		};
		rmse[3]=(rmse[0]+rmse[1]+rmse[2])/3;
		double psnr[]=
		{
			-20*log10(rmse[0]/255),
			-20*log10(rmse[1]/255),
			-20*log10(rmse[2]/255),
			-20*log10(rmse[3]/255),
		};
		printf(
			"RMSE  PSNR\n"
			"T %12.6lf  %12.6lf\n"
			"Y %12.6lf  %12.6lf\n"
			"U %12.6lf  %12.6lf\n"
			"V %12.6lf  %12.6lf\n"
			, rmse[3], psnr[3]
			, rmse[0], psnr[0]
			, rmse[1], psnr[1]
			, rmse[2], psnr[2]
		);
	}
#endif
	(void)&time_sec;
	(void)csize;
	return 0;
}
