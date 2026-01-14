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


	#define USE_L1
//	#define USE_RLE
	#define USE_DIVAC


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
	BLOCKX=4,
	BLOCKY=4,

	SHIFT=18,
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED

	XPAD=8,
	NCH=3,
	NROWS=4,
	NVAL=2,

	GRBITS=3,
#ifdef USE_RLE
	NCTX_DC=18,
	NCTX_RUN=1,
	NCTX_AC=1,
	NCTX=NCTX_DC+NCTX_RUN+NCTX_AC,
	NLEVELS=256,
#else
	NCTX_DC=18,
	NCTX_AC=1,
	NCTX=NCTX_DC+(BLOCKX*BLOCKY-1)*NCTX_AC,
	NLEVELS=512,
#endif

	PROBBITS=12,	//do NOT change this
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
#endif
static const float qtable_luma[]=
{
	16, 10, 12, 16,//v15
	10, 10, 16, 24,
	12, 16, 24, 24,
	16, 24, 24, 24,

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
};
static const float qtable_chroma[]=
{
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
};
#ifdef USE_RLE
static int zigzag[]=
{
	//0x0, 0x1, 0x2, 0x3,
	//0x4, 0x5, 0x6, 0x7,
	//0x8, 0x9, 0xA, 0xB,
	//0xC, 0xD, 0xE, 0xF,

	0x0,
	0x1, 0x4,
	0x8, 0x5, 0x2,
	0x3, 0x6, 0x9, 0xC,
	0xD, 0xA, 0x7,
	0xB, 0xE,
	0xF,
};
#endif
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
	__m128 a[8];

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

AWM_INLINE void dct4y4_fwd(float *block)//DCT-II
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
AWM_INLINE void dct4y4_inv(float *block)//DCT-III
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
AWM_INLINE void dct8y8_fwd(float *block)//DCT-II
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
AWM_INLINE void dct8y8_inv(float *block)//DCT-III
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

AWM_INLINE void dctii_4(float *data, const int stride)
{
	__m128 a[4], b[4];

	a[0]=_mm_load_ps(data+0*stride);
	a[1]=_mm_load_ps(data+1*stride);
	a[2]=_mm_load_ps(data+2*stride);
	a[3]=_mm_load_ps(data+3*stride);

	//Hbar
	b[0]=_mm_add_ps(a[0], a[3]);
	b[1]=_mm_add_ps(a[1], a[2]);
	b[2]=_mm_sub_ps(a[0], a[3]);
	b[3]=_mm_sub_ps(a[1], a[2]);

	//Wc
	b[2]=_mm_mul_ps(b[2], _mm_set1_ps(0.5411961001461969843997232053664f));//0.5 sec( 1 pi/(2*4))
	b[3]=_mm_mul_ps(b[3], _mm_set1_ps(1.3065629648763765278566431734272f));//0.5 sec( 3 pi/(2*4))

	//two half-DCTs
	a[0]=_mm_add_ps(b[0], b[1]);
	a[1]=_mm_sub_ps(b[0], b[1]);
	a[2]=_mm_add_ps(b[2], b[3]);
	a[3]=_mm_sub_ps(b[2], b[3]);

	//Bc
	a[2]=_mm_add_ps(_mm_mul_ps(a[2], _mm_set1_ps(1.4142135623730950488016887242097f)), a[3]);//sqrt2

	//PeoT
	_mm_store_ps(data+0*stride, a[0]);
	_mm_store_ps(data+1*stride, a[2]);
	_mm_store_ps(data+2*stride, a[1]);
	_mm_store_ps(data+3*stride, a[3]);
}
AWM_INLINE void dctiii_4(float *data, const int stride)
{
	__m128 a[4], b[4];

	//Peo
	a[0]=_mm_load_ps(data+0*stride);
	a[2]=_mm_load_ps(data+1*stride);
	a[1]=_mm_load_ps(data+2*stride);
	a[3]=_mm_load_ps(data+3*stride);

	//BcT
	a[3]=_mm_add_ps(a[3], a[2]);
	a[2]=_mm_mul_ps(a[2], _mm_set1_ps(1.4142135623730950488016887242097f));
	
	//two half-DCTs
	b[0]=_mm_add_ps(a[0], a[3]);
	b[1]=_mm_add_ps(a[1], a[2]);
	b[2]=_mm_sub_ps(a[0], a[3]);
	b[3]=_mm_sub_ps(a[1], a[2]);

	//Wc
	b[2]=_mm_mul_ps(b[2], _mm_set1_ps(0.5411961001461969843997232053664f));//0.5 sec( 1 pi/(2*4))
	b[3]=_mm_mul_ps(b[3], _mm_set1_ps(1.3065629648763765278566431734272f));//0.5 sec( 3 pi/(2*4))

	//HbarT
	a[0]=_mm_add_ps(b[0], b[2]);
	a[1]=_mm_add_ps(b[1], b[3]);
	a[2]=_mm_sub_ps(b[0], b[2]);
	a[3]=_mm_sub_ps(b[1], b[3]);

	_mm_store_ps(data+0*stride, a[0]);
	_mm_store_ps(data+1*stride, a[1]);
	_mm_store_ps(data+2*stride, a[2]);
	_mm_store_ps(data+3*stride, a[3]);
}
AWM_INLINE void dctii_8(float *data, const int stride)
{
	__m128 w78=_mm_set1_ps(2.5629154477415061787960862961777f);//Wij = sec(i*Pi/(2j))
	__m128 w58=_mm_set1_ps(0.8999762231364157046385095409419f);
	__m128 w38=_mm_set1_ps(0.6013448869350452805437218239092f);
	__m128 w18=_mm_set1_ps(0.5097955791041591689419398039878f);
	__m128 w34=_mm_set1_ps(1.3065629648763765278566431734272f);
	__m128 w14=_mm_set1_ps(0.5411961001461969843997232053664f);
	__m128 sqrt2=_mm_set1_ps(1.4142135623730950488016887242097f);

	__m128 a[8], b[8];

	a[0]=_mm_load_ps(data+0*stride);
	a[1]=_mm_load_ps(data+1*stride);
	a[2]=_mm_load_ps(data+2*stride);
	a[3]=_mm_load_ps(data+3*stride);
	a[4]=_mm_load_ps(data+4*stride);
	a[5]=_mm_load_ps(data+5*stride);
	a[6]=_mm_load_ps(data+6*stride);
	a[7]=_mm_load_ps(data+7*stride);

	b[0]=_mm_add_ps(a[0], a[7]);
	b[1]=_mm_add_ps(a[1], a[6]);
	b[2]=_mm_add_ps(a[2], a[5]);
	b[3]=_mm_add_ps(a[3], a[4]);
	b[4]=_mm_sub_ps(a[3], a[4]);
	b[5]=_mm_sub_ps(a[2], a[5]);
	b[6]=_mm_sub_ps(a[1], a[6]);
	b[7]=_mm_sub_ps(a[0], a[7]);

	b[4]=_mm_mul_ps(b[4], w78);
	b[5]=_mm_mul_ps(b[5], w58);
	b[6]=_mm_mul_ps(b[6], w38);
	b[7]=_mm_mul_ps(b[7], w18);

	a[0]=_mm_add_ps(b[0], b[3]);
	a[1]=_mm_add_ps(b[1], b[2]);
	a[2]=_mm_sub_ps(b[1], b[2]);
	a[3]=_mm_sub_ps(b[0], b[3]);
	a[4]=_mm_sub_ps(b[7], b[4]);
	a[5]=_mm_sub_ps(b[6], b[5]);
	a[6]=_mm_add_ps(b[6], b[5]);
	a[7]=_mm_add_ps(b[7], b[4]);

	a[2]=_mm_mul_ps(a[2], w34);
	a[3]=_mm_mul_ps(a[3], w14);
	a[4]=_mm_mul_ps(a[4], w14);
	a[5]=_mm_mul_ps(a[5], w34);

	b[0]=_mm_add_ps(a[0], a[1]);
	b[1]=_mm_sub_ps(a[0], a[1]);
	b[2]=_mm_add_ps(a[3], a[2]);
	b[3]=_mm_sub_ps(a[3], a[2]);
	b[4]=_mm_add_ps(a[4], a[5]);
	b[5]=_mm_sub_ps(a[4], a[5]);
	b[6]=_mm_add_ps(a[7], a[6]);
	b[7]=_mm_sub_ps(a[7], a[6]);

	b[3]=_mm_add_ps(_mm_mul_ps(b[3], sqrt2), b[2]);
	b[4]=_mm_add_ps(_mm_mul_ps(b[4], sqrt2), b[5]);

	b[7]=_mm_add_ps(_mm_mul_ps(b[7], sqrt2), b[4]);

	b[4]=_mm_add_ps(b[4], b[6]);
	b[6]=_mm_add_ps(b[6], b[5]);

	_mm_store_ps(data+0*stride, b[0]);
	_mm_store_ps(data+1*stride, b[7]);
	_mm_store_ps(data+2*stride, b[3]);
	_mm_store_ps(data+3*stride, b[4]);
	_mm_store_ps(data+4*stride, b[1]);
	_mm_store_ps(data+5*stride, b[6]);
	_mm_store_ps(data+6*stride, b[2]);
	_mm_store_ps(data+7*stride, b[5]);
}
AWM_INLINE void dctiii_8(float *data, const int stride)
{
	__m128 w78=_mm_set1_ps(2.5629154477415061787960862961777f);//Wij = sec(i*Pi/(2j))
	__m128 w58=_mm_set1_ps(0.8999762231364157046385095409419f);
	__m128 w38=_mm_set1_ps(0.6013448869350452805437218239092f);
	__m128 w18=_mm_set1_ps(0.5097955791041591689419398039878f);
	__m128 w34=_mm_set1_ps(1.3065629648763765278566431734272f);
	__m128 w14=_mm_set1_ps(0.5411961001461969843997232053664f);
	__m128 sqrt2=_mm_set1_ps(1.4142135623730950488016887242097f);

	__m128 a[8], b[8];

	a[0]=_mm_load_ps(data+0*stride);
	a[1]=_mm_load_ps(data+1*stride);
	a[2]=_mm_load_ps(data+2*stride);
	a[3]=_mm_load_ps(data+3*stride);
	a[4]=_mm_load_ps(data+4*stride);
	a[5]=_mm_load_ps(data+5*stride);
	a[6]=_mm_load_ps(data+6*stride);
	a[7]=_mm_load_ps(data+7*stride);

	b[0]=a[0];
	b[1]=_mm_mul_ps(sqrt2, a[1]);
	b[2]=a[2];
	b[3]=_mm_add_ps(a[1], a[3]);
	b[4]=a[4];
	b[5]=_mm_add_ps(a[3], a[5]);
	b[6]=a[6];
	b[7]=_mm_add_ps(a[5], a[7]);

	a[0]=b[0];
	a[1]=b[1];
	a[2]=_mm_mul_ps(sqrt2, b[2]);
	a[3]=_mm_mul_ps(sqrt2, b[3]);
	a[4]=b[4];
	a[5]=b[5];
	a[6]=_mm_add_ps(b[2], b[6]);
	a[7]=_mm_add_ps(b[3], b[7]);

	b[0]=_mm_add_ps(a[0], a[4]);
	b[1]=_mm_add_ps(a[1], a[5]);
	b[2]=_mm_add_ps(a[2], a[6]);
	b[3]=_mm_add_ps(a[3], a[7]);
	b[4]=_mm_sub_ps(a[0], a[4]);
	b[5]=_mm_sub_ps(a[1], a[5]);
	b[6]=_mm_sub_ps(a[2], a[6]);
	b[7]=_mm_sub_ps(a[3], a[7]);

	b[2]=_mm_mul_ps(b[2], w14);
	b[3]=_mm_mul_ps(b[3], w14);
	b[6]=_mm_mul_ps(b[6], w34);
	b[7]=_mm_mul_ps(b[7], w34);

	a[0]=_mm_add_ps(b[0], b[2]);
	a[1]=_mm_add_ps(b[1], b[3]);
	a[2]=_mm_sub_ps(b[0], b[2]);
	a[3]=_mm_sub_ps(b[1], b[3]);
	a[4]=_mm_add_ps(b[4], b[6]);
	a[5]=_mm_add_ps(b[5], b[7]);
	a[6]=_mm_sub_ps(b[4], b[6]);
	a[7]=_mm_sub_ps(b[5], b[7]);

	a[1]=_mm_mul_ps(a[1], w18);
	a[3]=_mm_mul_ps(a[3], w38);
	a[5]=_mm_mul_ps(a[5], w58);
	a[7]=_mm_mul_ps(a[7], w78);

	b[0]=_mm_add_ps(a[0], a[1]);
	b[1]=_mm_sub_ps(a[0], a[1]);
	b[2]=_mm_add_ps(a[2], a[3]);
	b[3]=_mm_sub_ps(a[2], a[3]);
	b[4]=_mm_add_ps(a[4], a[5]);
	b[5]=_mm_sub_ps(a[4], a[5]);
	b[6]=_mm_add_ps(a[6], a[7]);
	b[7]=_mm_sub_ps(a[6], a[7]);

	_mm_store_ps(data+0*stride, b[0]);
	_mm_store_ps(data+1*stride, b[4]);
	_mm_store_ps(data+2*stride, b[6]);
	_mm_store_ps(data+3*stride, b[2]);
	_mm_store_ps(data+4*stride, b[3]);
	_mm_store_ps(data+5*stride, b[7]);
	_mm_store_ps(data+6*stride, b[5]);
	_mm_store_ps(data+7*stride, b[1]);
}
AWM_INLINE void dctii_16(float *data, const int stride)
{
	__m128 a[16];

	//Hbar_16
	float *p1=data, *p2=data+15*stride;
	for(int k=0;k<8;++k)
	{
		__m128 t0=_mm_load_ps(p1);
		__m128 t1=_mm_load_ps(p2);
		p1+=stride;
		p2-=stride;
		a[k+0*8]=_mm_add_ps(t0, t1);
		a[k+1*8]=_mm_sub_ps(t0, t1);
	}

	//Wc_16
	a[0x8]=_mm_mul_ps(a[0x8], _mm_set1_ps(0.5024192861881557055116701192801f));//0.5 sec( 1 pi/(2*16))
	a[0x9]=_mm_mul_ps(a[0x9], _mm_set1_ps(0.5224986149396888806285753190567f));//0.5 sec( 3 pi/(2*16))
	a[0xA]=_mm_mul_ps(a[0xA], _mm_set1_ps(0.5669440348163577036805379151549f));//0.5 sec( 5 pi/(2*16))
	a[0xB]=_mm_mul_ps(a[0xB], _mm_set1_ps(0.6468217833599901295483601116520f));//0.5 sec( 7 pi/(2*16))
	a[0xC]=_mm_mul_ps(a[0xC], _mm_set1_ps(0.7881546234512502247339824871974f));//0.5 sec( 9 pi/(2*16))
	a[0xD]=_mm_mul_ps(a[0xD], _mm_set1_ps(1.0606776859903474713404517472331f));//0.5 sec(11 pi/(2*16))
	a[0xE]=_mm_mul_ps(a[0xE], _mm_set1_ps(1.7224470982383339278159153641566f));//0.5 sec(13 pi/(2*16))
	a[0xF]=_mm_mul_ps(a[0xF], _mm_set1_ps(5.1011486186891638581062454923454f));//0.5 sec(15 pi/(2*16))

	dctii_8((float*)(a+0*8), sizeof(__m128)/sizeof(float));
	dctii_8((float*)(a+1*8), sizeof(__m128)/sizeof(float));

	//Bc_16
	a[0x8]=_mm_add_ps(_mm_mul_ps(a[0x8], _mm_set1_ps(1.4142135623730950488016887242097f)), a[0x9]);//sqrt2
	a[0x9]=_mm_add_ps(a[0x9], a[0xA]);
	a[0xA]=_mm_add_ps(a[0xA], a[0xB]);
	a[0xB]=_mm_add_ps(a[0xB], a[0xC]);
	a[0xC]=_mm_add_ps(a[0xC], a[0xD]);
	a[0xD]=_mm_add_ps(a[0xD], a[0xE]);
	a[0xE]=_mm_add_ps(a[0xE], a[0xF]);

	//Peo16T
	for(int k=0;k<8;++k)
		_mm_store_ps(data+(2*k)*stride, a[k]);
	for(int k=0;k<8;++k)
		_mm_store_ps(data+(2*k+1)*stride, a[k+8]);
}
AWM_INLINE void dctiii_16(float *data, const int stride)
{
	__m128 a[16];

	//Peo16
	for(int k=0;k<8;++k)
		a[k]=_mm_load_ps(data+(2*k)*stride);
	for(int k=0;k<8;++k)
		a[k+8]=_mm_load_ps(data+(2*k+1)*stride);

	//Bc_16T
	a[0xF]=_mm_add_ps(a[0xF], a[0xE]);
	a[0xE]=_mm_add_ps(a[0xE], a[0xD]);
	a[0xD]=_mm_add_ps(a[0xD], a[0xC]);
	a[0xC]=_mm_add_ps(a[0xC], a[0xB]);
	a[0xB]=_mm_add_ps(a[0xB], a[0xA]);
	a[0xA]=_mm_add_ps(a[0xA], a[0x9]);
	a[0x9]=_mm_add_ps(a[0x9], a[0x8]);
	a[0x8]=_mm_mul_ps(a[0x8], _mm_set1_ps(1.4142135623730950488016887242097f));

	dctiii_8((float*)(a+0*8), sizeof(__m128)/sizeof(float));
	dctiii_8((float*)(a+1*8), sizeof(__m128)/sizeof(float));

	//Wc_16
	a[0x8]=_mm_mul_ps(a[0x8], _mm_set1_ps(0.5024192861881557055116701192801f));//0.5 sec( 1 pi/(2*16))
	a[0x9]=_mm_mul_ps(a[0x9], _mm_set1_ps(0.5224986149396888806285753190567f));//0.5 sec( 3 pi/(2*16))
	a[0xA]=_mm_mul_ps(a[0xA], _mm_set1_ps(0.5669440348163577036805379151549f));//0.5 sec( 5 pi/(2*16))
	a[0xB]=_mm_mul_ps(a[0xB], _mm_set1_ps(0.6468217833599901295483601116520f));//0.5 sec( 7 pi/(2*16))
	a[0xC]=_mm_mul_ps(a[0xC], _mm_set1_ps(0.7881546234512502247339824871974f));//0.5 sec( 9 pi/(2*16))
	a[0xD]=_mm_mul_ps(a[0xD], _mm_set1_ps(1.0606776859903474713404517472331f));//0.5 sec(11 pi/(2*16))
	a[0xE]=_mm_mul_ps(a[0xE], _mm_set1_ps(1.7224470982383339278159153641566f));//0.5 sec(13 pi/(2*16))
	a[0xF]=_mm_mul_ps(a[0xF], _mm_set1_ps(5.1011486186891638581062454923454f));//0.5 sec(15 pi/(2*16))

	//Hbar16T
	for(int k=0;k<8;++k)
		_mm_store_ps(data+k*stride, _mm_add_ps(a[k], a[k+8]));
	for(int k=0;k<8;++k)
		_mm_store_ps(data+(k+8)*stride, _mm_sub_ps(a[7-k], a[7-k+8]));
}
AWM_INLINE void dctii_32(float *data, const int stride)
{
	__m128 a[32];

	//Hbar_32
	float *p1=data, *p2=data+31*stride;
	for(int k=0;k<16;++k)
	{
		__m128 t0=_mm_load_ps(p1);
		__m128 t1=_mm_load_ps(p2);
		p1+=stride;
		p2-=stride;
		a[k+0*16]=_mm_add_ps(t0, t1);
		a[k+1*16]=_mm_sub_ps(t0, t1);
	}

	//Wc_32
	a[0x10]=_mm_mul_ps(a[0x10], _mm_set1_ps( 0.5006029982351963013455041067664f));//0.5 sec( 1 pi/(2*32))
	a[0x11]=_mm_mul_ps(a[0x11], _mm_set1_ps( 0.5054709598975436599844445856070f));//0.5 sec( 3 pi/(2*32))
	a[0x12]=_mm_mul_ps(a[0x12], _mm_set1_ps( 0.5154473099226245469749513056493f));//0.5 sec( 5 pi/(2*32))
	a[0x13]=_mm_mul_ps(a[0x13], _mm_set1_ps( 0.5310425910897841744757339323572f));//0.5 sec( 7 pi/(2*32))
	a[0x14]=_mm_mul_ps(a[0x14], _mm_set1_ps( 0.5531038960344445278293808381371f));//0.5 sec( 9 pi/(2*32))
	a[0x15]=_mm_mul_ps(a[0x15], _mm_set1_ps( 0.5829349682061338736738307012526f));//0.5 sec(11 pi/(2*32))
	a[0x16]=_mm_mul_ps(a[0x16], _mm_set1_ps( 0.6225041230356648161572561567628f));//0.5 sec(13 pi/(2*32))
	a[0x17]=_mm_mul_ps(a[0x17], _mm_set1_ps( 0.6748083414550057460259687110410f));//0.5 sec(15 pi/(2*32))
	a[0x18]=_mm_mul_ps(a[0x18], _mm_set1_ps( 0.7445362710022984497769811919729f));//0.5 sec(17 pi/(2*32))
	a[0x19]=_mm_mul_ps(a[0x19], _mm_set1_ps( 0.8393496454155270387392637466254f));//0.5 sec(19 pi/(2*32))
	a[0x1A]=_mm_mul_ps(a[0x1A], _mm_set1_ps( 0.9725682378619606936976894140525f));//0.5 sec(21 pi/(2*32))
	a[0x1B]=_mm_mul_ps(a[0x1B], _mm_set1_ps( 1.1694399334328849551557702840422f));//0.5 sec(23 pi/(2*32))
	a[0x1C]=_mm_mul_ps(a[0x1C], _mm_set1_ps( 1.4841646163141662772433269374281f));//0.5 sec(25 pi/(2*32))
	a[0x1D]=_mm_mul_ps(a[0x1D], _mm_set1_ps( 2.0577810099534115508565544797104f));//0.5 sec(27 pi/(2*32))
	a[0x1E]=_mm_mul_ps(a[0x1E], _mm_set1_ps( 3.4076084184687187857011913334591f));//0.5 sec(29 pi/(2*32))
	a[0x1F]=_mm_mul_ps(a[0x1F], _mm_set1_ps(10.1900081235480568112121092010356f));//0.5 sec(31 pi/(2*32))

	dctii_16((float*)(a+0*16), sizeof(__m128)/sizeof(float));
	dctii_16((float*)(a+1*16), sizeof(__m128)/sizeof(float));

	//Bc_16
	a[0x10]=_mm_add_ps(_mm_mul_ps(a[0x10], _mm_set1_ps(1.4142135623730950488016887242097f)), a[0x11]);//sqrt2
	a[0x11]=_mm_add_ps(a[0x11], a[0x12]);
	a[0x12]=_mm_add_ps(a[0x12], a[0x13]);
	a[0x13]=_mm_add_ps(a[0x13], a[0x14]);
	a[0x14]=_mm_add_ps(a[0x14], a[0x15]);
	a[0x15]=_mm_add_ps(a[0x15], a[0x16]);
	a[0x16]=_mm_add_ps(a[0x16], a[0x17]);
	a[0x17]=_mm_add_ps(a[0x17], a[0x18]);
	a[0x18]=_mm_add_ps(a[0x18], a[0x19]);
	a[0x19]=_mm_add_ps(a[0x19], a[0x1A]);
	a[0x1A]=_mm_add_ps(a[0x1A], a[0x1B]);
	a[0x1B]=_mm_add_ps(a[0x1B], a[0x1C]);
	a[0x1C]=_mm_add_ps(a[0x1C], a[0x1D]);
	a[0x1D]=_mm_add_ps(a[0x1D], a[0x1E]);
	a[0x1E]=_mm_add_ps(a[0x1E], a[0x1F]);

	//Peo16T
	for(int k=0;k<16;++k)
		_mm_store_ps(data+(2*k)*stride, a[k]);
	for(int k=0;k<16;++k)
		_mm_store_ps(data+(2*k+1)*stride, a[k+16]);
}
AWM_INLINE void dctiii_32(float *data, const int stride)
{
	__m128 a[32];

	//Peo16
	for(int k=0;k<16;++k)
		a[k]=_mm_load_ps(data+(2*k)*stride);
	for(int k=0;k<16;++k)
		a[k+16]=_mm_load_ps(data+(2*k+1)*stride);

	//Bc_16T
	a[0x1F]=_mm_add_ps(a[0x1F], a[0x1E]);
	a[0x1E]=_mm_add_ps(a[0x1E], a[0x1D]);
	a[0x1D]=_mm_add_ps(a[0x1D], a[0x1C]);
	a[0x1C]=_mm_add_ps(a[0x1C], a[0x1B]);
	a[0x1B]=_mm_add_ps(a[0x1B], a[0x1A]);
	a[0x1A]=_mm_add_ps(a[0x1A], a[0x19]);
	a[0x19]=_mm_add_ps(a[0x19], a[0x18]);
	a[0x18]=_mm_add_ps(a[0x18], a[0x17]);
	a[0x17]=_mm_add_ps(a[0x17], a[0x16]);
	a[0x16]=_mm_add_ps(a[0x16], a[0x15]);
	a[0x15]=_mm_add_ps(a[0x15], a[0x14]);
	a[0x14]=_mm_add_ps(a[0x14], a[0x13]);
	a[0x13]=_mm_add_ps(a[0x13], a[0x12]);
	a[0x12]=_mm_add_ps(a[0x12], a[0x11]);
	a[0x11]=_mm_add_ps(a[0x11], a[0x10]);
	a[0x10]=_mm_mul_ps(a[0x10], _mm_set1_ps(1.4142135623730950488016887242097f));

	dctiii_16((float*)(a+0*16), sizeof(__m128)/sizeof(float));
	dctiii_16((float*)(a+1*16), sizeof(__m128)/sizeof(float));

	//Wc_16
	a[0x10]=_mm_mul_ps(a[0x10], _mm_set1_ps( 0.5006029982351963013455041067664f));//0.5 sec( 1 pi/(2*32))
	a[0x11]=_mm_mul_ps(a[0x11], _mm_set1_ps( 0.5054709598975436599844445856070f));//0.5 sec( 3 pi/(2*32))
	a[0x12]=_mm_mul_ps(a[0x12], _mm_set1_ps( 0.5154473099226245469749513056493f));//0.5 sec( 5 pi/(2*32))
	a[0x13]=_mm_mul_ps(a[0x13], _mm_set1_ps( 0.5310425910897841744757339323572f));//0.5 sec( 7 pi/(2*32))
	a[0x14]=_mm_mul_ps(a[0x14], _mm_set1_ps( 0.5531038960344445278293808381371f));//0.5 sec( 9 pi/(2*32))
	a[0x15]=_mm_mul_ps(a[0x15], _mm_set1_ps( 0.5829349682061338736738307012526f));//0.5 sec(11 pi/(2*32))
	a[0x16]=_mm_mul_ps(a[0x16], _mm_set1_ps( 0.6225041230356648161572561567628f));//0.5 sec(13 pi/(2*32))
	a[0x17]=_mm_mul_ps(a[0x17], _mm_set1_ps( 0.6748083414550057460259687110410f));//0.5 sec(15 pi/(2*32))
	a[0x18]=_mm_mul_ps(a[0x18], _mm_set1_ps( 0.7445362710022984497769811919729f));//0.5 sec(17 pi/(2*32))
	a[0x19]=_mm_mul_ps(a[0x19], _mm_set1_ps( 0.8393496454155270387392637466254f));//0.5 sec(19 pi/(2*32))
	a[0x1A]=_mm_mul_ps(a[0x1A], _mm_set1_ps( 0.9725682378619606936976894140525f));//0.5 sec(21 pi/(2*32))
	a[0x1B]=_mm_mul_ps(a[0x1B], _mm_set1_ps( 1.1694399334328849551557702840422f));//0.5 sec(23 pi/(2*32))
	a[0x1C]=_mm_mul_ps(a[0x1C], _mm_set1_ps( 1.4841646163141662772433269374281f));//0.5 sec(25 pi/(2*32))
	a[0x1D]=_mm_mul_ps(a[0x1D], _mm_set1_ps( 2.0577810099534115508565544797104f));//0.5 sec(27 pi/(2*32))
	a[0x1E]=_mm_mul_ps(a[0x1E], _mm_set1_ps( 3.4076084184687187857011913334591f));//0.5 sec(29 pi/(2*32))
	a[0x1F]=_mm_mul_ps(a[0x1F], _mm_set1_ps(10.1900081235480568112121092010356f));//0.5 sec(31 pi/(2*32))

	//Hbar16T
	for(int k=0;k<16;++k)
		_mm_store_ps(data+k*stride, _mm_add_ps(a[k], a[k+16]));
	for(int k=0;k<16;++k)
		_mm_store_ps(data+(k+16)*stride, _mm_sub_ps(a[15-k], a[15-k+16]));
}

AWM_INLINE void transpose4x4(float *block)
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
AWM_INLINE void transpose8x8(float *block)
{
	__m128 a[4], b[4];

	//	A	B
	//	C	D
	
	//transpose A
	a[0]=_mm_load_ps(block+0*8);
	a[1]=_mm_load_ps(block+1*8);
	a[2]=_mm_load_ps(block+2*8);
	a[3]=_mm_load_ps(block+3*8);
	
	_MM_TRANSPOSE4_PS(a[0], a[1], a[2], a[3]);

	_mm_store_ps(block+0*8, a[0]);
	_mm_store_ps(block+1*8, a[1]);
	_mm_store_ps(block+2*8, a[2]);
	_mm_store_ps(block+3*8, a[3]);

	//transpose B and C
	a[0]=_mm_load_ps(block+0*8+4);
	a[1]=_mm_load_ps(block+1*8+4);
	a[2]=_mm_load_ps(block+2*8+4);
	a[3]=_mm_load_ps(block+3*8+4);
	
	_MM_TRANSPOSE4_PS(a[0], a[1], a[2], a[3]);

	b[0]=_mm_load_ps(block+4*8);
	b[1]=_mm_load_ps(block+5*8);
	b[2]=_mm_load_ps(block+6*8);
	b[3]=_mm_load_ps(block+7*8);
	_mm_store_ps(block+4*8, a[0]);
	_mm_store_ps(block+5*8, a[1]);
	_mm_store_ps(block+6*8, a[2]);
	_mm_store_ps(block+7*8, a[3]);

	_MM_TRANSPOSE4_PS(b[0], b[1], b[2], b[3]);

	_mm_store_ps(block+0*8+4, a[0]);
	_mm_store_ps(block+1*8+4, a[1]);
	_mm_store_ps(block+2*8+4, a[2]);
	_mm_store_ps(block+3*8+4, a[3]);

	//transpose D
	a[0]=_mm_load_ps(block+4*8+4);
	a[1]=_mm_load_ps(block+5*8+4);
	a[2]=_mm_load_ps(block+6*8+4);
	a[3]=_mm_load_ps(block+7*8+4);
	
	_MM_TRANSPOSE4_PS(a[0], a[1], a[2], a[3]);

	_mm_store_ps(block+4*8+4, a[0]);
	_mm_store_ps(block+5*8+4, a[1]);
	_mm_store_ps(block+6*8+4, a[2]);
	_mm_store_ps(block+7*8+4, a[3]);
}
AWM_INLINE void transpose16x16(float *block)
{
	__m128 a[4], b[4];
	
#define TRANSPOSE_DIAG(IDX)\
	do\
	{\
		a[0]=_mm_load_ps(block+(0+4*(IDX))*16+4*(IDX));\
		a[1]=_mm_load_ps(block+(1+4*(IDX))*16+4*(IDX));\
		a[2]=_mm_load_ps(block+(2+4*(IDX))*16+4*(IDX));\
		a[3]=_mm_load_ps(block+(3+4*(IDX))*16+4*(IDX));\
		_MM_TRANSPOSE4_PS(a[0], a[1], a[2], a[3]);\
		_mm_store_ps(block+(0+4*(IDX))*16+4*(IDX), a[0]);\
		_mm_store_ps(block+(1+4*(IDX))*16+4*(IDX), a[1]);\
		_mm_store_ps(block+(2+4*(IDX))*16+4*(IDX), a[2]);\
		_mm_store_ps(block+(3+4*(IDX))*16+4*(IDX), a[3]);\
	}while(0)
#define TRANSPOSE_NOND(Y, X)\
	do\
	{\
		a[0]=_mm_load_ps(block+(0+4*(Y))*16+4*(X));\
		a[1]=_mm_load_ps(block+(1+4*(Y))*16+4*(X));\
		a[2]=_mm_load_ps(block+(2+4*(Y))*16+4*(X));\
		a[3]=_mm_load_ps(block+(3+4*(Y))*16+4*(X));\
		b[0]=_mm_load_ps(block+(0+4*(X))*16+4*(Y));\
		b[1]=_mm_load_ps(block+(1+4*(X))*16+4*(Y));\
		b[2]=_mm_load_ps(block+(2+4*(X))*16+4*(Y));\
		b[3]=_mm_load_ps(block+(3+4*(X))*16+4*(Y));\
		_MM_TRANSPOSE4_PS(a[0], a[1], a[2], a[3]);\
		_MM_TRANSPOSE4_PS(b[0], b[1], b[2], b[3]);\
		_mm_store_ps(block+(0+4*(Y))*16+4*(X), b[0]);\
		_mm_store_ps(block+(1+4*(Y))*16+4*(X), b[1]);\
		_mm_store_ps(block+(2+4*(Y))*16+4*(X), b[2]);\
		_mm_store_ps(block+(3+4*(Y))*16+4*(X), b[3]);\
		_mm_store_ps(block+(0+4*(X))*16+4*(Y), a[0]);\
		_mm_store_ps(block+(1+4*(X))*16+4*(Y), a[1]);\
		_mm_store_ps(block+(2+4*(X))*16+4*(Y), a[2]);\
		_mm_store_ps(block+(3+4*(X))*16+4*(Y), a[3]);\
	}while(0)
	
	//	0	1	2	3
	//	4	5	6	7
	//	8	9	A	B
	//	C	D	E	F
	TRANSPOSE_DIAG(0);
	TRANSPOSE_DIAG(1);
	TRANSPOSE_DIAG(2);
	TRANSPOSE_DIAG(3);
	TRANSPOSE_NOND(0, 1);
	TRANSPOSE_NOND(0, 2);
	TRANSPOSE_NOND(0, 3);
	TRANSPOSE_NOND(1, 2);
	TRANSPOSE_NOND(1, 3);
	TRANSPOSE_NOND(2, 3);

#undef  TRANSPOSE_DIAG
#undef  TRANSPOSE_NOND
}
AWM_INLINE void transpose32x32(float *block)
{
	__m128 a[4], b[4];
	
#define TRANSPOSE_DIAG(IDX)\
	do\
	{\
		a[0]=_mm_load_ps(block+(0+4*(IDX))*32+4*(IDX));\
		a[1]=_mm_load_ps(block+(1+4*(IDX))*32+4*(IDX));\
		a[2]=_mm_load_ps(block+(2+4*(IDX))*32+4*(IDX));\
		a[3]=_mm_load_ps(block+(3+4*(IDX))*32+4*(IDX));\
		_MM_TRANSPOSE4_PS(a[0], a[1], a[2], a[3]);\
		_mm_store_ps(block+(0+4*(IDX))*32+4*(IDX), a[0]);\
		_mm_store_ps(block+(1+4*(IDX))*32+4*(IDX), a[1]);\
		_mm_store_ps(block+(2+4*(IDX))*32+4*(IDX), a[2]);\
		_mm_store_ps(block+(3+4*(IDX))*32+4*(IDX), a[3]);\
	}while(0)
#define TRANSPOSE_NOND(Y, X)\
	do\
	{\
		a[0]=_mm_load_ps(block+(0+4*(Y))*32+4*(X));\
		a[1]=_mm_load_ps(block+(1+4*(Y))*32+4*(X));\
		a[2]=_mm_load_ps(block+(2+4*(Y))*32+4*(X));\
		a[3]=_mm_load_ps(block+(3+4*(Y))*32+4*(X));\
		b[0]=_mm_load_ps(block+(0+4*(X))*32+4*(Y));\
		b[1]=_mm_load_ps(block+(1+4*(X))*32+4*(Y));\
		b[2]=_mm_load_ps(block+(2+4*(X))*32+4*(Y));\
		b[3]=_mm_load_ps(block+(3+4*(X))*32+4*(Y));\
		_MM_TRANSPOSE4_PS(a[0], a[1], a[2], a[3]);\
		_MM_TRANSPOSE4_PS(b[0], b[1], b[2], b[3]);\
		_mm_store_ps(block+(0+4*(Y))*32+4*(X), b[0]);\
		_mm_store_ps(block+(1+4*(Y))*32+4*(X), b[1]);\
		_mm_store_ps(block+(2+4*(Y))*32+4*(X), b[2]);\
		_mm_store_ps(block+(3+4*(Y))*32+4*(X), b[3]);\
		_mm_store_ps(block+(0+4*(X))*32+4*(Y), a[0]);\
		_mm_store_ps(block+(1+4*(X))*32+4*(Y), a[1]);\
		_mm_store_ps(block+(2+4*(X))*32+4*(Y), a[2]);\
		_mm_store_ps(block+(3+4*(X))*32+4*(Y), a[3]);\
	}while(0)
	
	//	00	01	02	03	04	05	06	07
	//	08	09	0A	0B	0C	0D	0E	0F
	//	10	11	12	13	14	15	16	17
	//	18	19	1A	1B	1C	1D	1E	1F
	//	20	21	22	23	24	25	26	27
	//	28	29	2A	2B	2C	2D	2E	2F
	//	30	31	32	33	34	35	36	37
	//	38	39	3A	3B	3C	3D	3E	3F
	TRANSPOSE_DIAG(0);
	TRANSPOSE_DIAG(1);
	TRANSPOSE_DIAG(2);
	TRANSPOSE_DIAG(3);
	TRANSPOSE_DIAG(4);
	TRANSPOSE_DIAG(5);
	TRANSPOSE_DIAG(6);
	TRANSPOSE_DIAG(7);
	TRANSPOSE_NOND(0, 1);
	TRANSPOSE_NOND(0, 2);
	TRANSPOSE_NOND(0, 3);
	TRANSPOSE_NOND(0, 4);
	TRANSPOSE_NOND(0, 5);
	TRANSPOSE_NOND(0, 6);
	TRANSPOSE_NOND(0, 7);

	TRANSPOSE_NOND(1, 2);
	TRANSPOSE_NOND(1, 3);
	TRANSPOSE_NOND(1, 4);
	TRANSPOSE_NOND(1, 5);
	TRANSPOSE_NOND(1, 6);
	TRANSPOSE_NOND(1, 7);

	TRANSPOSE_NOND(2, 3);
	TRANSPOSE_NOND(2, 4);
	TRANSPOSE_NOND(2, 5);
	TRANSPOSE_NOND(2, 6);
	TRANSPOSE_NOND(2, 7);

	TRANSPOSE_NOND(3, 4);
	TRANSPOSE_NOND(3, 5);
	TRANSPOSE_NOND(3, 6);
	TRANSPOSE_NOND(3, 7);

	TRANSPOSE_NOND(4, 5);
	TRANSPOSE_NOND(4, 6);
	TRANSPOSE_NOND(4, 7);

	TRANSPOSE_NOND(5, 6);
	TRANSPOSE_NOND(5, 7);

	TRANSPOSE_NOND(6, 7);

#undef  TRANSPOSE_DIAG
#undef  TRANSPOSE_NOND
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
AWM_INLINE void dc_update(int16_t **rows, int32_t *weights, int32_t *estim, int32_t p1, int32_t curr, int32_t error)
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
	rows[0][1]=(2*eW+((error<<1^error>>31)<<GRBITS)+eNEEE)>>2;
}

#if 0
#define TEST
static void dcttest()
{
	ALIGN(16) float data[16]=
	{
		1,  3,  3,  1,
		3,  3,  6,  3,
		3,  6, 12,  6,
		1,  3,  6,  3,
	};
	ALIGN(16) float d0[16], diff[16];
	for(int k=0;k<16;++k)
		data[k]=(float)rand();
		//data[k]=(int)k;

	memcpy(d0, data, sizeof(d0));

	dct4y4_fwd(data);
	transpose4x4(data);
	dct4y4_fwd(data);

	dct4y4_inv(data);
	transpose4x4(data);
	dct4y4_inv(data);
	
	for(int k=0;k<16;++k)
		diff[k]=d0[k]-data[k]/4;
	for(int k=0;k<16;++k)
	{
		if(fabsf(diff[k])>1e-2)
			CRASH("");
	}

	printf("SUCCESS\n");
	exit(0);
}
#endif


#ifdef PROFILE_SIZE
static double csizes[3][BLOCKX*BLOCKY*2]={0};
#endif

static uint32_t hweight[3*NCTX];
static uint32_t hists[3*NCTX*NLEVELS];
ALIGN(32) static float blocks[3][BLOCKX*BLOCKY]={0}, qtable1[BLOCKX*BLOCKY]={0}, qtable2[BLOCKX*BLOCKY]={0};

#ifndef USE_DIVAC
uint16_t enccdf[3*NCTX*NLEVELS];//enc
uint32_t CDF2sym[3*NCTX<<PROBBITS];//dec
static void normalize_enc(int kc, int rescale)
{
	uint64_t invsum=0;
	uint32_t cdf=0, sum2=0, bypass=0;
	uint32_t *currh=hists+256*kc;
	uint16_t *currs=enccdf+256*kc;
	int ks=0;

	invsum=hweight[kc];
	if(invsum<256)
	{
		invsum+=256;
		bypass=1;
	}
	invsum=(0x100000000ULL+invsum-1)/invsum*((1ULL<<PROBBITS)-256);

	for(ks=0, cdf=0, sum2=0;ks<=255;++ks)
	{
		uint32_t freq=currh[ks];
		currs[ks]=(uint32_t)(cdf*invsum>>32)+ks;
		cdf+=freq+bypass;
		if(rescale)
		{
			freq>>=1;
			currh[ks]=freq;
			sum2+=freq;
		}
	}
	if(rescale)
		hweight[kc]=sum2;
}
static void normalize_dec(int kc, int rescale)
{
	uint64_t invsum=0;
	uint32_t cdf=0, sum2=0, bypass=0;
	uint32_t *currh=hists+256*kc, *currCDF2sym=CDF2sym+((ptrdiff_t)kc<<PROBBITS);
	uint16_t CDF[257];
	int ks=0;
	
	invsum=hweight[kc];
	if(invsum<256)
	{
		invsum+=256;
		bypass=1;
	}
	invsum=(0x100000000ULL+invsum-1)/invsum*((1ULL<<PROBBITS)-256);

	for(ks=0, cdf=0, sum2=0;ks<=255;++ks)
	{
		uint32_t freq=currh[ks];
		CDF[ks]=(uint32_t)(cdf*invsum>>32)+ks;
#ifdef _MSC_VER
		if(CDF[ks]>(1<<PROBBITS))
			CRASH("");
#endif
		cdf+=freq+bypass;
		if(rescale)
		{
			freq>>=1;
			currh[ks]=freq;
			sum2+=freq;
		}
	}
	CDF[256]=1<<PROBBITS;
	if(rescale)
		hweight[kc]=sum2;
	for(ks=0, cdf=0;ks<256;++ks)
	{
		uint32_t val=0;
		int freq=CDF[ks], end=CDF[ks+1];
		for(val=(end-freq)<<(PROBBITS+8)|freq<<8|ks;freq<end;)//for AC
			currCDF2sym[freq++]=val;
	}
}
#endif
int c47_codec(int argc, char **argv)
{
#ifdef TEST
	dcttest();//
#endif

	const uint16_t tag='4'|'7'<<8;

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
#ifdef USE_RLE
	char rle[3][BLOCKX*BLOCKY*2+1]={0};
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
		//for(int ky=0;ky<ih;++ky)
		//{
		//	for(int kx=0;kx<iw;++kx)
		//	{
		//		image[3*(iw*ky+kx)+0]=64*(kx^ky);
		//		image[3*(iw*ky+kx)+1]=64*(kx^ky);
		//		image[3*(iw*ky+kx)+2]=64*(kx^ky);
		//	}
		//}
		//for(ptrdiff_t k=0;k<usize;++k)//
		//	image[k]=rand()>>8&15;

		guide_save(image, iw, ih);
#endif
	}
	else
	{
		c=fread(stream, 1, ccap, fsrc);
	}
	fclose(fsrc);

	for(int k=0;k<NCH*NPREDS;++k)
		((int32_t*)weights)[k]=(1<<SHIFT)/NPREDS;
	memset(pixels, 0, psize);
	memset(hists, 0, sizeof(hists));
	memset(hweight, 0, sizeof(hweight));
	if(fwd)
	{
		for(int k=0;k<BLOCKX*BLOCKY;++k)
		{
			qtable1[k]=1/(4*qtable_luma[k]);
			qtable2[k]=1/(4*qtable_chroma[k]);
		}
#ifndef USE_DIVAC
		for(int kc=0;kc<3*NCTX;++kc)
			normalize_enc(kc, 0);
#endif
	}
	else
	{
		for(int k=0;k<BLOCKX*BLOCKY;++k)
		{
			qtable1[k]=qtable_luma[k];
			qtable2[k]=qtable_chroma[k];
		}
#ifndef USE_DIVAC
		for(int kc=0;kc<3*NCTX;++kc)
			normalize_dec(kc, 0);
#endif
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
#ifdef USE_RLE
			int remaining=0;
#endif
#ifdef USE_DIVAC
			int den=0;
#ifdef USE_RLE
			int nlevels=0;
#endif
#endif

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
				dct4y4_fwd(blocks[0]);
				dct4y4_fwd(blocks[1]);
				dct4y4_fwd(blocks[2]);
				transpose4x4(blocks[0]);
				transpose4x4(blocks[1]);
				transpose4x4(blocks[2]);
				dct4y4_fwd(blocks[0]);
				dct4y4_fwd(blocks[1]);
				dct4y4_fwd(blocks[2]);
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
				for(int kc=0;kc<3;++kc)for(int k=0;k<BLOCKX*BLOCKY;++k)sum_after[kc][k]+=abs(((int*)blocks[kc])[k]);
#endif

				//RLE: {DC, run, AC, run, AC, run, ...AC/run}
#ifdef USE_RLE
				for(int kc=0;kc<3;++kc)
				{
					int kd=1, val;
					{
						float *idiotgcc1=blocks[kc];
						size_t idiotgcc2=(size_t)idiotgcc1;
						val=((int*)idiotgcc2)[0];//DC can be zero
					}
					CLAMP2(val, -128, 127);
					rle[kc][0]=val;
					rle[kc][1]=0;
					for(int ks=1;ks<BLOCKX*BLOCKY;++ks)
					{
						val=((int*)blocks[kc])[zigzag[ks]];
						CLAMP2(val, -128, 127);
						if(!val)
						{
							++rle[kc][kd];
#ifdef _MSC_VER
							if(rle[kc][kd]>BLOCKX*BLOCKY-1)
								CRASH("");
#endif
						}
						else
						{
							++kd;
							rle[kc][kd++]=(val<<1^val>>31)-1;//AC cannot be zero outsize of a run
							rle[kc][kd]=0;
						}
					}
				}
#endif
			}
			for(int kc=0;kc<3;++kc)
			{
#ifdef USE_RLE
				remaining=BLOCKX*BLOCKY;
				for(int k=0;remaining>0;++k)
#else
				for(int k=0;k<BLOCKX*BLOCKY;++k)
#endif
				{
					int32_t p1=0;
#ifdef USE_DIVAC
#ifdef USE_RLE
					nlevels=k&1?remaining+1:NLEVELS;
#endif
#endif

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
				//		printf("");

					if(!k)
					{
						//if(ky==0&&kx==0&&!kc)//
						//	printf("");

						dc_predict(rows, weights[kc], estim, &p1, &pred, &ctx);
					}
					else
					{
						pred=0;
#ifdef USE_RLE
						ctx=(k&1)+NCTX_DC;
#else
						ctx=NCTX_DC-1+k;
#endif
						//ctx=FLOOR_LOG2(rows[0][k+1]+1);
						//if(ctx>NCTX_AC-1)
						//	ctx=NCTX_AC-1;
						//ctx=0;
						//ctx+=NCTX*kc+NCTX_AC*k+(NCTX_DC-NCTX_AC);
					}
					//if(ctx>=NCTX)//
					//	CRASH("");
					ctx+=NCTX*kc;
#ifdef USE_DIVAC
					uint32_t *currhist=hists+NLEVELS*ctx;
					den=hweight[ctx]+NLEVELS;
#ifdef USE_RLE
					if(k&1)
					{
						den=nlevels;
						for(int k=0;k<nlevels;++k)
							den+=currhist[k];
					}
#endif
#endif
					if(fwd)
					{
#ifdef FIFOVAL
						valfifo_enqueue(ctx<<24^p1<<12^pred);
#endif
#ifdef USE_RLE
						curr=rle[kc][k];
						//curr=((int*)blocks[kc])[k];
//#ifdef _MSC_VER
//						if(abs(curr)>128)//
//							CRASH("");
//#endif
//						CLAMP2(curr, -128, 127);

						error=curr-pred;
						sym=(uint8_t)error;
#else
						curr=((int*)blocks[kc])[k];
						CLAMP2(curr, -(NLEVELS>>1), (NLEVELS>>1)-1);
						error=((curr-pred+(NLEVELS>>1))&(NLEVELS-1))-(NLEVELS>>1);
						sym=error<<1^error>>31;
#endif
#ifndef USE_DIVAC
						cdf=enccdf[ctx<<8|sym];
						freq=1<<PROBBITS;
						if(sym<255)
							freq=enccdf[ctx<<8|(sym+1)];
						freq-=cdf;
#ifdef FIFOVAL
						valfifo_enqueue(ctx<<24^p1<<12^pred);
						valfifo_enqueue(freq<<(PROBBITS+8)|cdf<<8|sym);
#endif
#endif

						if(range<=0xFFFF)
						{
							*(uint32_t*)streamptr=(uint32_t)(low>>32);
							streamptr+=4;
							low<<=32;
							range=range<<32|0xFFFFFFFF;
							if(range>~low)
								range=~low;
						}
#ifdef USE_DIVAC
						int t=0;
						cdf=0;
						for(;;)
						{
							freq=currhist[t]+1;
							if(t>=sym)
								break;
							cdf+=freq;
							++t;
						}
#ifdef _MSC_VER
						if(!freq||(uint32_t)freq>(uint32_t)den||(uint32_t)cdf>(uint32_t)den)
							CRASH("");
#endif
#ifdef FIFOVAL
						valfifo_enqueue(freq<<16|cdf);
#endif
						low+=range*cdf/den;
						range=range*freq/den-1;
#ifdef PROFILE_SIZE
						csizes[kc][k]+=-log2((double)freq/den);
#endif
#else
						low+=range*cdf>>PROBBITS;
						range=(range*freq>>PROBBITS)-1;
#ifdef PROFILE_SIZE
						csizes[kc][k]+=-log2(freq*(1./(1<<PROBBITS)));
#endif
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
#ifdef USE_DIVAC
						int c2=(int)(((code-low+1)*den-1)/range);
						sym=0;
						cdf=0;
						for(;;)
						{
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
#else
						int c2=(int)(((code-low)<<PROBBITS|((1ULL<<PROBBITS)-1))/range);//stall: DIV64
#ifdef _MSC_VER
						if(c2>>PROBBITS)
							CRASH("Dec error X %d  Y %d  C%d", kx, ky, kc);
#endif
						uint32_t info=CDF2sym[ctx<<PROBBITS|c2];//stall: cache miss
						sym=(uint8_t)info;
						cdf=info<<PROBBITS>>(32-PROBBITS);
						freq=info>>(PROBBITS+8);
#ifdef FIFOVAL
						valfifo_check(ctx<<24^p1<<12^pred);
						valfifo_check(info);
#endif
						low+=range*cdf>>PROBBITS;
						range=(range*freq>>PROBBITS)-1;
#endif
#ifdef USE_RLE
						rle[kc][k]=curr=(int8_t)sym+pred;
#else
						error=sym>>1^-(sym&1);
						curr=((error+pred+(NLEVELS>>1))&(NLEVELS-1))-(NLEVELS>>1);
						((int*)blocks[kc])[k]=curr;
#endif
						error=curr-pred;
					}
					
					++currhist[sym];
					++hweight[ctx];
#ifdef USE_DIVAC
					if(hweight[ctx]>=0x4000)
					{
#ifdef USE_RLE
						nlevels=k&1?BLOCKX*BLOCKY+1:NLEVELS;
						den=0;
						for(int k=0;k<nlevels;++k)
							den+=currhist[k]>>=1;
						hweight[ctx]=den;
#else
						den=0;
						for(int k=0;k<NLEVELS;++k)
							den+=currhist[k]>>=1;
						hweight[ctx]=den;
#endif
					}
#else
					int count=hweight[ctx];
					if(count>=NLEVELS&&!(count&(count-1)))
					{
						if(fwd)
							normalize_enc(ctx, hweight[ctx]>=0x2000);
						else
							normalize_dec(ctx, hweight[ctx]>=0x2000);
					}
#endif
					if(!k)
					{
						dc_update(rows, weights[kc], estim, p1, curr, error);
						rows[0][0]=curr;
						rows[0]+=NROWS*NVAL;
						rows[1]+=NROWS*NVAL;
						rows[2]+=NROWS*NVAL;
						rows[3]+=NROWS*NVAL;
					}
					//else
					//	rows[0][k+1]=(2*rows[0][k+1-1*NCH*NROWS*NVAL]+sym+rows[1][k+1+3*NCH*NROWS*NVAL])>>3;
					//if(k==BLOCKX*BLOCKY-1)
					//{
					//	rows[0]+=NROWS*NVAL;
					//	rows[1]+=NROWS*NVAL;
					//	rows[2]+=NROWS*NVAL;
					//	rows[3]+=NROWS*NVAL;
					//}
#ifdef USE_RLE
					if(k&1)
						remaining-=curr;
					else
						--remaining;
#endif
				}
#ifdef USE_RLE
#ifdef _MSC_VER
				if(remaining)
					CRASH("");
#endif
#endif
			}
			if(!fwd)
			{
#ifdef USE_RLE
				for(int kc=0;kc<3;++kc)
				{
					{
						float *idiotgcc1=blocks[kc];
						size_t idiotgcc2=(size_t)idiotgcc1;
						((int*)idiotgcc2)[0]=rle[kc][0];
					}
					//((int*)(size_t)blocks[kc])[0]=rle[kc][0];
					for(int ks=1, kd=1;kd<BLOCKX*BLOCKY;++kd)
					{
						if(rle[kc][ks])
						{
							((int*)blocks[kc])[zigzag[kd]]=0;
							--rle[kc][ks];
						}
						else
						{
							int val;

							++ks;
							val=(uint8_t)rle[kc][ks++]+1;
							((int*)blocks[kc])[zigzag[kd]]=val>>1^-(val&1);
						}
					}
				}
#endif
				cvti2f(blocks[0]);
				cvti2f(blocks[1]);
				cvti2f(blocks[2]);
				gain(blocks[0], qtable1);
				gain(blocks[1], qtable2);
				gain(blocks[2], qtable2);
				dct4y4_inv(blocks[0]);
				dct4y4_inv(blocks[1]);
				dct4y4_inv(blocks[2]);
				transpose4x4(blocks[0]);
				transpose4x4(blocks[1]);
				transpose4x4(blocks[2]);
				dct4y4_inv(blocks[0]);
				dct4y4_inv(blocks[1]);
				dct4y4_inv(blocks[2]);
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
#if 0
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
	exit(0);
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
			for(int k=0;k<BLOCKX*BLOCKY*2;++k)
				sum+=csizes[kc][k];
		}
		if(!sum)
			sum=1;
		sum=100/sum;
		for(int kc=0;kc<3;++kc)
		{
			double sum2=0;
			for(int k=0;k<BLOCKX*BLOCKY*2;++k)
				sum2+=csizes[kc][k];
			printf("%c: %12.2lf\n", "YUV"[kc], sum2);
			for(int k=0;k<BLOCKX*BLOCKY*2;++k)
				printf("  %12.2lf %8.4lf%% %c\n"
					, csizes[kc][k]
					, csizes[kc][k]*sum
					, k&1?'R':'C'
				);
			//for(int ky=0;ky<BLOCKY;++ky)
			//{
			//	for(int kx=0;kx<BLOCKY;++kx)
			//		printf("  %12.2lf %8.4lf%%"
			//			, csizes[kc][ky*BLOCKX+kx]
			//			, csizes[kc][ky*BLOCKX+kx]*sum
			//		);
			//	printf("\n");
			//}
			//printf("\n");
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
