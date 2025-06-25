#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<math.h>
#include<sys/stat.h>
#include<immintrin.h>
#if defined _MSC_VER || defined _WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#include<process.h>
#define THREAD_CALL __stdcall
typedef unsigned THREAD_RET;
#else
#include<pthread.h>
#define THREAD_CALL
typedef void *THREAD_RET;
#endif


//	#define RELEASE
//	#define PRING_BLOCKSIZES
	#define PROFILE_TIME

#ifdef _MSC_VER
	#define LOUD			//size & time
#endif
#if defined _MSC_VER && !defined RELEASE
	#define DISABLE_MT

//	#define INTERCEPT_FWRITE
	#define PROFILE_SIZE

//	#define ESTIMATE_SIZE		//DEBUG		checks for zero frequency, visualizes context usage
	#define ENABLE_GUIDE		//DEBUG		checks interleaved pixels
	#define ANS_VAL			//DEBUG

//	#define PRINT_L1_BOUNDS
#endif


#define XMAXBLOCK 1024
#define YMAXBLOCK 1024
#define XLANES 4
#define YLANES 4
#define NLANES 16

#define GRBITS 3
#define NCTX 18		//NCTX*3+3  total

#define PROBBITS 12	//12 bit max		CDF2sym {freq<<20 | bias<<8 | sym}
#define RANS_STATE_BITS 31
#define RANS_RENORM_BITS 16

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
#ifdef _MSC_VER
#define	ALIGN(N) __declspec(align(N))
#define AWM_INLINE __forceinline static
#else
#define	ALIGN(N) __attribute__((aligned(N)))
#define AWM_INLINE __attribute__((always_inline)) inline static
#ifndef _countof
#define _countof(A) (sizeof(A)/sizeof(*(A)))
#endif
#endif
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
static void memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
{
	size_t copied;
	char *d=(char*)dst;
	const char *s=(const char*)src;
#ifdef _DEBUG
	if(!dstbytes||!srcbytes)
		return;
#endif
	if(dstbytes<srcbytes)
	{
		memcpy(dst, src, dstbytes);
		return;
	}
	copied=srcbytes;
	memcpy(d, s, srcbytes);
	while((copied<<1)<=dstbytes)
	{
		memcpy(d+copied, d, copied);
		copied<<=1;
	}
	if(copied<dstbytes)
		memcpy(d+copied, d, dstbytes-copied);
}
#define FILLMEM(PTR, DATA, ASIZE, ESIZE)\
	do\
	{\
		*(PTR)=(DATA);\
		memfill((PTR)+1, PTR, (ASIZE)-(ESIZE), ESIZE);\
	}while(0)
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
static int query_cpu_cores(void)
{
#ifdef _WIN32
	SYSTEM_INFO info;
	GetNativeSystemInfo(&info);
	return info.dwNumberOfProcessors;
#else
	return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}
#ifdef PROFILE_TIME
#define PROFLIST\
	PROFLABEL(read)\
	PROFLABEL(main)\
	PROFLABEL(rem)\
	PROFLABEL(write)\

typedef enum _ProfLabel
{
#define PROFLABEL(LABEL) PROF_##LABEL,
	PROFLIST
#undef  PROFLABEL
	PROF_COUNT,
} ProfLabel;
static const char *profnames[]=
{
#define PROFLABEL(LABEL) #LABEL,
	PROFLIST
#undef  PROFLABEL
};
static double g_profinfo[128]={0};
static ptrdiff_t g_volume[128]={0};
static double g_t=0;
static int g_profstart=0, g_profend=0;
static void prof_start()
{
	g_t=time_sec();
	memset(g_profinfo, 0, sizeof(g_profinfo));
	memset(g_volume, 0, sizeof(g_volume));
}
static void prof_checkpoint(int idx, ptrdiff_t volume)
{
	if(idx>=_countof(g_profinfo))
	{
		CRASH("Profiler OOB");
		return;
	}
	{
		double t2=time_sec();
		g_profinfo[idx]+=t2-g_t;
		g_volume[idx]+=volume;
		g_t=t2;
		if(g_profend<idx+1)
			g_profend=idx+1;
	}
}
#define PROF(LABEL, VOLUME) prof_checkpoint(PROF_##LABEL, VOLUME)
static void prof_print(void)
{
	double sum=0;
	for(int k=g_profstart;k<g_profend;++k)
		sum+=g_profinfo[k];
	const int scale=4;//ms
	printf("1 char = %d ms\n", scale);
	printf("|");
	double csum=0;
	int prev=0;
	for(int k=g_profstart;k<g_profend;++k)
	{
		double val=g_profinfo[k];
		csum+=val;
		int next=(int)(csum*1000/scale);
		int nstars=next-prev;
		prev=next;
		for(int k2=0;k2<nstars/2;++k2)
			printf("-");
		printf("%d", k-g_profstart+1);
		for(int k2=nstars/2;k2<nstars;++k2)
			printf("-");
		printf("|");
	}
	printf("\n");
	ptrdiff_t smax=0;
	for(int k=g_profstart;k<g_profend;++k)
	{
		double elapsed=g_profinfo[k];
		printf("%8.4lf%%  %12.6lf sec  %10td bytes  %12.6lf MB/s  %2d  %s\n",
			100.*elapsed/sum,
			elapsed,
			g_volume[k],
			g_volume[k]/(elapsed*1024*1024),
			k-g_profstart+1,
			profnames[k]
		);
		if(smax<g_volume[k])
			smax=g_volume[k];
	}
	printf("\n");
	printf("%lf sec  %12.6lf MB/s\n"
		, sum
		, smax/(sum*1024*1024)
	);
	printf("\n");
}
#else
#define prof_start()
#define PROF(...)
#define prof_print()
#endif

#ifdef INTERCEPT_FWRITE
static ptrdiff_t my_fwrite(const void *buf, ptrdiff_t esize, ptrdiff_t count, FILE *f)
{
	ptrdiff_t ret;

	double t=time_sec();
	ret=fwrite(buf, esize, count, f);
	t=time_sec()-t;
	printf("%12lld bytes  %12.6lf sec  %12.6lf MB/s\n", (int64_t)esize*count, t, (double)esize*count/(1024*1024*t));
	return ret;
}
#define fwrite my_fwrite
#endif

#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static double g_sqe[3]={0};
static void guide_save(unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		CRASH("");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		CRASH("");
		printf("");
	}
}
static void guide_update_sqe(const unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	double diff;
	diff=g_image[idx+0]-image[idx+0]; g_sqe[0]+=diff*diff;
	diff=g_image[idx+1]-image[idx+1]; g_sqe[1]+=diff*diff;
	diff=g_image[idx+2]-image[idx+2]; g_sqe[2]+=diff*diff;
}
#else
#define guide_save(...)
#define guide_check(...)
#define guide_update_sqe(...)
#endif

#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
	OCH(CX40) OCH(C0X4) OCH(C40X)\
	OCH(CX31) OCH(C3X1) OCH(C31X)\
	OCH(CX13) OCH(C1X3) OCH(C13X)\
	OCH(CX22) OCH(C2X2) OCH(C22X)
#if 0
#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
	OCH(Y310) OCH(Y031) OCH(Y103)\
	OCH(Y301) OCH(Y130) OCH(Y013)\
	OCH(Y211) OCH(Y121) OCH(Y112)\
	OCH(CX40) OCH(C0X4) OCH(C40X)\
	OCH(CX31) OCH(C3X1) OCH(C31X)\
	OCH(CX13) OCH(C1X3) OCH(C13X)\
	OCH(CX22) OCH(C2X2) OCH(C22X)
#endif
typedef enum _OCHIndex
{
#define OCH(X) OCH_##X,
	OCHLIST
#undef  OCH
	OCH_COUNT,
	OCH_C4X0=OCH_CX40,
	OCH_C04X=OCH_C0X4,
	OCH_CX04=OCH_C40X,
	OCH_R=OCH_Y400,
	OCH_G=OCH_Y040,
	OCH_B=OCH_Y004,
	OCH_BG=OCH_C04X,
	OCH_BR=OCH_C40X,
	OCH_RG=OCH_CX40,
	OCH_RB=OCH_CX04,
	OCH_GB=OCH_C0X4,
	OCH_GR=OCH_C4X0,
	OCH_R1=OCH_CX13,
	OCH_G1=OCH_C3X1,
	OCH_B1=OCH_C13X,
	OCH_R2=OCH_CX22,
	OCH_G2=OCH_C2X2,
	OCH_B2=OCH_C22X,
	OCH_R3=OCH_CX31,
	OCH_G3=OCH_C1X3,
	OCH_B3=OCH_C31X,
} OCHIndex;
static const char *och_names[]=
{
#define OCH(X) #X,
	OCHLIST
#undef  OCH
};
typedef enum _RCTInfoIdx
{
	II_OCH_Y,
	II_OCH_U,
	II_OCH_V,

	II_PERM_Y,
	II_PERM_U,
	II_PERM_V,

	II_COEFF_U_SUB_Y,
	II_COEFF_V_SUB_Y,
	II_COEFF_V_SUB_U,

//	II_COEFF_Y_SUB_U,
//	II_COEFF_Y_SUB_V,
//	II_COEFF_U_SUB_V_NBLI,
//	II_COEFF_V_SUB_U_NBLI,

	II_COUNT,
} RCTInfoIdx;
//YUV = RCT * RGB	watch out for permutation in last row
//luma: averaging	chroma: subtraction
//example: _400_40X_3X1 == [1 0 0; -1 0 1; -3/4 1 -1/4]
#define RCTLIST\
	RCT(_400_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_400_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4)\
	RCT(_400_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0)\
	RCT(_040_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0)\
	RCT(_040_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4)\
	RCT(_004_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4)\
	RCT(_004_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0)\
	RCT(_040_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0)\
	RCT(_040_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4)\
	RCT(_040_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4)\
	RCT(_004_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0)\
	RCT(_004_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4)\
	RCT(_004_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4)\
	RCT(_400_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0)\
	RCT(_400_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4)\
	RCT(_400_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)\
	RCT(_400_0X0_13X,	OCH_R,		OCH_G,		OCH_B1,		0, 1, 2,	0,  1, 3)\
	RCT(_400_4X0_13X,	OCH_R,		OCH_GR,		OCH_B1,		0, 1, 2,	4,  1, 3)\
	RCT(_400_00X_3X1,	OCH_R,		OCH_B,		OCH_G1,		0, 2, 1,	0,  3, 1)\
	RCT(_400_40X_3X1,	OCH_R,		OCH_BR,		OCH_G1,		0, 2, 1,	4,  3, 1)\
	RCT(_040_00X_X13,	OCH_G,		OCH_B,		OCH_R1,		1, 2, 0,	0,  1, 3)\
	RCT(_040_04X_X13,	OCH_G,		OCH_BG,		OCH_R1,		1, 2, 0,	4,  1, 3)\
	RCT(_040_X40_13X,	OCH_G,		OCH_RG,		OCH_B1,		1, 0, 2,	4,  3, 1)\
	RCT(_004_X04_3X1,	OCH_B,		OCH_RB,		OCH_G1,		2, 0, 1,	4,  1, 3)\
	RCT(_004_04X_X13,	OCH_B,		OCH_GB,		OCH_R1,		2, 1, 0,	4,  3, 1)\
	RCT(_400_0X0_22X,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,  2, 2)\
	RCT(_400_4X0_22X,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	4,  2, 2)\
	RCT(_400_00X_2X2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,  2, 2)\
	RCT(_400_40X_2X2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	4,  2, 2)\
	RCT(_040_00X_X22,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,  2, 2)\
	RCT(_040_04X_X22,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	4,  2, 2)\
	RCT(_040_X40_22X,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	4,  2, 2)\
	RCT(_004_X04_2X2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	4,  2, 2)\
	RCT(_004_0X4_X22,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	4,  2, 2)\
	RCT(_400_0X0_31X,	OCH_R,		OCH_G,		OCH_B3,		0, 1, 2,	0,  3, 1)\
	RCT(_400_4X0_31X,	OCH_R,		OCH_GR,		OCH_B3,		0, 1, 2,	4,  3, 1)\
	RCT(_400_00X_1X3,	OCH_R,		OCH_B,		OCH_G3,		0, 2, 1,	0,  1, 3)\
	RCT(_400_40X_1X3,	OCH_R,		OCH_BR,		OCH_G3,		0, 2, 1,	4,  1, 3)\
	RCT(_040_00X_X31,	OCH_G,		OCH_B,		OCH_R3,		1, 2, 0,	0,  3, 1)\
	RCT(_040_04X_X31,	OCH_G,		OCH_BG,		OCH_R3,		1, 2, 0,	4,  3, 1)\
	RCT(_040_X40_31X,	OCH_G,		OCH_RG,		OCH_B3,		1, 0, 2,	4,  1, 3)\
	RCT(_004_X04_1X3,	OCH_B,		OCH_RB,		OCH_G3,		2, 0, 1,	4,  3, 1)\
	RCT(_004_0X4_X31,	OCH_B,		OCH_GB,		OCH_R3,		2, 1, 0,	4,  1, 3)
#if 0
	RCT(_211_4X0_40X,	OCH_Y211,	OCH_C4X0,	OCH_C40X,	0, 1, 2,	4,  4, 0,	1, 1, 0, 0)\
	RCT(_211_4X0_31X,	OCH_Y211,	OCH_C4X0,	OCH_C31X,	0, 1, 2,	4,  4, 0,	1, 1, 0, 1)\
	RCT(_211_3X1_40X,	OCH_Y211,	OCH_C3X1,	OCH_C40X,	0, 1, 2,	4,  4, 0,	1, 1, 1, 0)\
	RCT(_310_4X0_40X,	OCH_Y310,	OCH_C4X0,	OCH_C40X,	0, 1, 2,	4,  4, 0,	1, 0, 0, 0)\
	RCT(_310_4X0_31X,	OCH_Y310,	OCH_C4X0,	OCH_C31X,	0, 1, 2,	4,  4, 0,	1, 0, 0, 1)\
	RCT(_310_3X1_40X,	OCH_Y310,	OCH_C3X1,	OCH_C40X,	0, 1, 2,	4,  4, 0,	1, 0, 1, 0)\
	RCT(_301_4X0_40X,	OCH_Y301,	OCH_C4X0,	OCH_C40X,	0, 1, 2,	4,  4, 0,	0, 1, 0, 0)\
	RCT(_301_4X0_31X,	OCH_Y301,	OCH_C4X0,	OCH_C31X,	0, 1, 2,	4,  4, 0,	0, 1, 0, 1)\
	RCT(_301_3X1_40X,	OCH_Y301,	OCH_C3X1,	OCH_C40X,	0, 1, 2,	4,  4, 0,	0, 1, 1, 0)\
	RCT(_121_04X_X40,	OCH_Y121,	OCH_C04X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	1, 1, 0, 0)\
	RCT(_121_04X_X31,	OCH_Y121,	OCH_C04X,	OCH_CX31,	1, 2, 0,	4,  4, 0,	1, 1, 0, 1)\
	RCT(_121_13X_X40,	OCH_Y121,	OCH_C13X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	1, 1, 1, 0)\
	RCT(_031_04X_X40,	OCH_Y031,	OCH_C04X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	0, 1, 0, 0)\
	RCT(_031_04X_X31,	OCH_Y031,	OCH_C04X,	OCH_CX31,	1, 2, 0,	4,  4, 0,	0, 1, 0, 1)\
	RCT(_031_13X_X40,	OCH_Y031,	OCH_C13X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	0, 1, 1, 0)\
	RCT(_130_40X_X40,	OCH_Y130,	OCH_C04X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	1, 0, 0, 0)\
	RCT(_130_40X_X31,	OCH_Y130,	OCH_C04X,	OCH_CX31,	1, 2, 0,	4,  4, 0,	1, 0, 0, 1)\
	RCT(_130_31X_X40,	OCH_Y130,	OCH_C13X,	OCH_CX40,	1, 2, 0,	4,  4, 0,	1, 0, 1, 0)\
	RCT(_112_X04_0X4,	OCH_Y112,	OCH_CX04,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	1, 1, 0, 0)\
	RCT(_112_X04_1X3,	OCH_Y112,	OCH_CX04,	OCH_C1X3,	2, 0, 1,	4,  4, 0,	1, 1, 0, 1)\
	RCT(_112_X13_0X4,	OCH_Y112,	OCH_CX13,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	1, 1, 1, 0)\
	RCT(_013_X04_0X4,	OCH_Y013,	OCH_CX04,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	0, 1, 0, 0)\
	RCT(_013_X04_1X3,	OCH_Y013,	OCH_CX04,	OCH_C1X3,	2, 0, 1,	4,  4, 0,	0, 1, 0, 1)\
	RCT(_013_X13_0X4,	OCH_Y013,	OCH_CX13,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	0, 1, 1, 0)\
	RCT(_103_X40_0X4,	OCH_Y103,	OCH_CX04,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	1, 0, 0, 0)\
	RCT(_103_X40_1X3,	OCH_Y103,	OCH_CX04,	OCH_C1X3,	2, 0, 1,	4,  4, 0,	1, 0, 0, 1)\
	RCT(_103_X31_0X4,	OCH_Y103,	OCH_CX13,	OCH_C0X4,	2, 0, 1,	4,  4, 0,	1, 0, 1, 0)
#endif
typedef enum _RCTIndex
{
#define RCT(LABEL, ...) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][II_COUNT]=
{
#define RCT(LABEL, ...) {__VA_ARGS__},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, ...) #LABEL,
	RCTLIST
#undef  RCT
};

static int16_t nthreads, effort, dist, profiler, fwd;
static int32_t iw, ih, rowstride;
static ptrdiff_t res, usize;
static uint16_t xrem, yrem;
static int32_t qw, qh;
static uint16_t xnblocks, ynblocks, nblocks, blockxstops[512], blockystops[512];
static uint8_t *image, *stream;
static ptrdiff_t srcsize;
static int bsizes[512];
ALIGN(32) static int ans_permute[256][8];

#ifdef ANS_VAL
typedef struct _ANSVALHeader
{
	uint16_t esize, count;
	unsigned idx;
	struct _ANSVALHeader *above, *below;
	unsigned char data[];
} ANSVALNode;
static ANSVALNode *debugstack=0;
static int ansvalidx=0, ansvalcount=0;
static void ansval_push(const void *data, int esize, int count)
{
	int size=sizeof(ANSVALNode)+esize*count;
	ANSVALNode *node=(ANSVALNode*)malloc(size);
	if(!node)
	{
		printf("Alloc error\n");
		exit(1);
	}
	memset(node, 0, size);
	node->esize=esize;
	node->count=count;
	node->idx=ansvalidx++;
	node->above=0;
	node->below=debugstack;
	if(debugstack)
		debugstack->above=node;
	memcpy(node->data, data, size-sizeof(ANSVALNode));
	debugstack=node;
	ansvalcount=ansvalidx;
}
static void ansval_printr(const void *data, int esize, int count, const void *xdata)//print elements in reverse because little-endian
{
	const unsigned char *p=(const unsigned char*)data, *p2=(const unsigned char*)xdata;
	int size=count*esize, k;
	for(k=0;k<size;k+=esize)
	{
		int k2=esize-1;
		printf(" ");
		for(;k2>=0;--k2)
		{
			int val=p[k+k2];
			if(p2)
				val^=p2[k+k2];
			if(p2&&!val)
				printf("--");
			else
				printf("%02X", val);
		}
	}
	printf("\n");
}
static void* ansval_ptrguard(const void *start, const void *end, const void *ptr, ptrdiff_t nbytes)
{
	size_t istart=(size_t)start, iend=(size_t)end;
	ptrdiff_t size=iend-istart;
	size_t ip1=(size_t)ptr, ip2=ip1+nbytes;
	int problems[]=
	{
		size<0,
		(size_t)(ip1-istart)>=(size_t)size,
		(size_t)(ip2-istart)>=(size_t)size,
	};
	if(problems[0]||problems[1]||problems[2])
	{
		printf("\nOOB\n");
		printf("  inc     %+16td bytes\n", nbytes);
		printf("  start   %016zd  %16d\n", istart, 0);
		if(nbytes<0)
		{
			printf("  after   %016zd  %16td%s\n", ip2, ip2-istart, problems[2]?"  <-":"");
			printf("  before  %016zd  %16td%s\n", ip1, ip1-istart, problems[1]?"  <-":"");
		}
		else
		{
			printf("  before  %016zd  %16td%s\n", ip1, ip1-istart, problems[1]?"  <-":"");
			printf("  after   %016zd  %16td%s\n", ip2, ip2-istart, problems[2]?"  <-":"");
		}
		printf("  end     %016zd  %16td%s\n", iend, size, problems[0]?"  <-":"");
		CRASH("\n");
		return 0;
	}
	return (void*)(nbytes<0?ip2:ip1);
}
static void ansval_check(const void *data, int esize, int count)
{
	--ansvalidx;
	if(!debugstack)
	{
		printf("Debug stack is empty\n");
		ansval_printr(data, esize, count, 0);
		CRASH("");
	}
	else if(debugstack->esize!=esize||debugstack->count!=count||memcmp(data, debugstack->data, esize*count))
	{
		printf("\n\nValidation Error  [enc ^ | v dec]  total %d\n", ansvalcount);
		if(debugstack->above)
		{
			ANSVALNode *node=debugstack->above;
			if(node->above)
			{
				ANSVALNode *node2=node->above;
				printf("[%10d] Verified:   ", node2->idx);
				ansval_printr(node2->data, node2->esize, node2->count, 0);
			}
			printf("[%10d] Verified:   ", node->idx);
			ansval_printr(node->data, node->esize, node->count, 0);
			printf("\n");
		}

		printf("[%10d] Original:   ", debugstack->idx);
		ansval_printr(debugstack->data, esize, count, 0);

		printf("[%10d] Corrupt:    ", debugstack->idx);
		ansval_printr(data, esize, count, 0);
		
		if(debugstack->esize==esize&&debugstack->count==count)
		{
			printf("[%10d] XOR:        ", debugstack->idx);
			ansval_printr(debugstack->data, esize, count, data);
		}
		if(debugstack->below)
		{
			ANSVALNode *node=debugstack->below;
			printf("\n");
			printf("[%10d] Below:      ", node->idx);
			ansval_printr(node->data, node->esize, node->count, 0);
			if(node->below)
			{
				node=node->below;
				printf("[%10d] Below:      ", node->idx);
				ansval_printr(node->data, node->esize, node->count, 0);
			}
		}
		printf("\n\n");
		CRASH("");
	}
	if(debugstack->below)
		debugstack=debugstack->below;
}
#endif
//LIFO Bypass Coder
typedef struct _BitPackerLIFO//bwd enc / fwd dec
{
	unsigned long long state;
	int enc_nwritten, dec_navailable;//bitcounts, only for tracking renorms
	unsigned char *dstbwdptr;
	const unsigned char *srcfwdptr, *streamend;
} BitPackerLIFO;
AWM_INLINE void bitpacker_enc_init(BitPackerLIFO *ec, const unsigned char *bufstart, unsigned char *bufptr0_OOB)
{
	memset(ec, 0, sizeof(*ec));
	ec->state=1ULL<<32;
	ec->enc_nwritten=33;
	ec->streamend=bufstart;
	ec->dstbwdptr=bufptr0_OOB;
}
AWM_INLINE void bitpacker_dec_init(BitPackerLIFO *ec, const unsigned char *bufptr0_start, const unsigned char *bufend)
{
	memset(ec, 0, sizeof(*ec));
	ec->srcfwdptr=bufptr0_start+8;
	ec->streamend=bufend;
	ec->state=*(const unsigned long long*)bufptr0_start;
	ec->dec_navailable=64-(int)_lzcnt_u64(ec->state);
}
AWM_INLINE void bitpacker_enc_flush(BitPackerLIFO *ec)
{
	ec->dstbwdptr-=8;
#ifdef _DEBUG
	if(ec->dstbwdptr<ec->streamend)
		CRASH("IntPacker Encoder OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
	*(unsigned long long*)ec->dstbwdptr=ec->state;
}
AWM_INLINE void bitpacker_enc(BitPackerLIFO *ec, int inbits, int sym)
{
#ifdef _DEBUG
	if(!inbits)
		CRASH("BitPacker inbits=0");
#endif
	//renorm then push inbits
	ec->enc_nwritten+=inbits;
	if(ec->enc_nwritten>64)//renorm on overflow
	{
		ec->enc_nwritten-=32;
		ec->dstbwdptr-=4;
#ifdef _DEBUG
		if(ec->dstbwdptr<ec->streamend)
			CRASH("IntPacker OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
		*(unsigned*)ec->dstbwdptr=(unsigned)ec->state;
		ec->state>>=32;
#ifdef ANS_VAL
		ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
	}
	ec->state=ec->state<<inbits|sym;
#ifdef ANS_VAL
	ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
}
AWM_INLINE int bitpacker_dec(BitPackerLIFO *ec, int outbits)
{
#ifdef _DEBUG
	if(!outbits)
		CRASH("BitPacker outbits=0");
#endif
	int sym=ec->state&((1ULL<<outbits)-1);

	//pop outbits then renorm
#ifdef ANS_VAL
	ansval_check(&ec->state, sizeof(ec->state), 1);
#endif
	ec->dec_navailable-=outbits;
	ec->state>>=outbits;
	if(ec->dec_navailable<=32)
	{
#ifdef ANS_VAL
		ansval_check(&ec->state, sizeof(ec->state), 1);
#endif
		ec->dec_navailable+=32;
#ifdef _DEBUG
		if(ec->srcfwdptr+4>ec->streamend)
			CRASH("IntPacker OOB:  srcfwdptr = 0x%016zX >= 0x%016zX", ec->srcfwdptr, ec->streamend);
#endif
		ec->state=ec->state<<32|*(const unsigned*)ec->srcfwdptr;
		ec->srcfwdptr+=4;
	}
	return sym;
}
//https://github.com/rygorous/ryg_rans
//https://github.com/samtools/htscodecs
typedef struct _rANS_SIMD_SymInfo	//16 bytes/level	4KB/ctx = 1<<12 bytes
{
	unsigned smax, invf, cdf;
	uint16_t negf, sh;
} rANS_SIMD_SymInfo;
static void enc_hist2stats(int *hist, rANS_SIMD_SymInfo *syminfo, unsigned long long *ctxmask, int ctxidx)
{
	int sum=0, count=0;
	for(int ks=0;ks<256;++ks)
	{
		int freq=hist[ks];
		sum+=freq;
		count+=freq!=0;
	}
	int rare=sum<12*256/8;
	*ctxmask|=(unsigned long long)rare<<ctxidx;
#ifdef ESTIMATE_SIZE
	int count0=count, sum0=sum;
#endif
	if(rare)
	{
		for(int ks=0;ks<256;++ks)//bypass
			hist[ks]=1;
		sum=256;
		count=256;
	}
	else if(count==1)//disallow degenerate distribution
	{
		for(int ks=0;ks<256;++ks)
		{
			int freq=hist[ks];
			if(freq==(1<<PROBBITS))
			{
				--freq;
				if(!ks)
					++hist[ks+1];
				else
					++hist[ks-1];
				break;
			}
		}
		count=2;
	}
	int sum2=0;
	for(int ks=0, ks2=0;ks<256;++ks)//absent symbols get zero freqs
	{
		int freq=hist[ks];
		hist[ks]=(int)(sum2*((1ULL<<PROBBITS)-count)/sum)+ks2;
		ks2+=freq!=0;
		sum2+=freq;
	}
	//for(int ks=0;ks<256;++ks)//never allows zero freqs	INEFFICIENT
	//{
	//	int freq=hist[ks];
	//	hist[ks]=(int)(sum2*((1ULL<<PROBBITS)-256)/sum)+ks;
	//	sum2+=freq;
	//}
#ifdef ESTIMATE_SIZE
	double e=sum0;
	if(count==count0)
	{
		double norm=1./0x1000;
		e=0;
		for(int ks=0;ks<256;++ks)//estimate
		{
			int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
			if(freq)
			{
				double p=freq*norm;
				e-=p*log2(p);
			}
			if(e!=e)
				CRASH("");
		}
		e*=sum/8.;
	}
	if(ctxidx&&!(ctxidx%NCTX))
		printf("\n");
	printf("%c  ctx %3d  %12.2lf / %9d bytes%8.2lf%%  %3d %s",
		ctxidx<3*NCTX?"YUV"[ctxidx/NCTX]:"yuv"[ctxidx-3*NCTX],
		ctxidx, e, sum0, 100.*e/sum0, count0, count==count0?"levels":"bypass"
	);
	if(count==count0&&count<256)
	{
		printf(" %3d", count);
		int fmax=0;
		for(int ks=0;ks<256;++ks)
		{
			int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
			if(fmax<freq)
				fmax=freq;
		}
		for(int ks=0;ks<256;++ks)
		{
			int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks];
			if(!(ks&15))
				printf(" ");

			int shade=48+freq*(255-48)/fmax;
			colorprintf(shade<<16|shade<<8|shade, freq?0x808080:COLORPRINTF_BK_DEFAULT, "%c", "0123456789ABCDEF"[ks&15]);
			//int shade=freq*255/fmax;
			//colorprintf(freq?0xFFFFFF:0x808080, shade<<16|0<<8|shade, "%c", "0123456789ABCDEF"[ks&15]);

			//printf("%c", freq?"0123456789ABCDEF"[ks&15]:'-');
		}
#if 0
		int printmissing=count>128, printcount=printmissing?256-count:count;
		if(printmissing)
			printf(" MISSING %3d: ", printcount);
		else
			printf("            : ");
		//printf(" %3d %-7s: ", printcount, printmissing?"MISSING":"       ");
		for(int ks=0, printed=0;ks<256;++ks)
		{
			int ks2=((ks>>1^-(ks&1))+128)&255;
			int freq=(ks2<256-1?hist[ks2+1]:1<<PROBBITS)-hist[ks2];
			if(printmissing!=(freq!=0))
			{
				printf(" %02X", ks2);
				//++printed;
				//if(printed&1)
				//	printf("%02x", ks2);
				//else
				//	printf("%02X", ks2);
				//if(printed>=90)
				//{
				//	printf("...%+4d more", printcount-printed);
				//	break;
				//}
			}
		}
#endif
	}
	printf("\n");
#if 0
	if(ctxidx==26)
	{
		const int amplitude=512;
		printf("Context %d: (1 star = %d steps)\n", ctxidx, 4096/amplitude);
		for(int ks=0;ks<256;++ks)
		{
			int freq=(ks<256-1?hist[ks+1]:1<<PROBBITS)-hist[ks], nstars=freq*amplitude>>PROBBITS;
			printf("%3d %4d ", ks, freq);
			for(int k=0;k<nstars;++k)
				printf("*");
			printf("\n");
		}
	}
#endif
#endif
	int next=1<<PROBBITS;
	for(int ks=255;ks>=0;--ks)
	{
		rANS_SIMD_SymInfo *info=syminfo+ks;
		int curr=hist[ks];
		int freq=next-curr;
		next=curr;
		hist[ks]=freq;
		info->smax=(freq<<(RANS_STATE_BITS-PROBBITS))-1;//rescale freq to match the rANS state, and decrement to use _mm256_cmpgt_epi32 instead of '>='
		info->cdf=curr;
		info->negf=(1<<PROBBITS)-freq;
		//encoding:  state  =  q<<16|(cdf+r)
		//div-free:  state  =  q*M+cdf+state-q*freq  =  state+q*(M-freq)+cdf  =  state+(state*invf>>sh)*(M-freq)+cdf
		//sh = FLOOR_LOG2(freq)+32
		//invf = ceil(2^sh/freq)		state is 31 bits
		if(freq<2)
		{
			//freq=1
			//ideally  q = state*inv(1)>>sh(1) = state*2^32>>32
			//here  q' = state*(2^32-1)>>32 = floor(state-state/2^32) = state-1  if  1 <= x < 2^32
			//enc  state = (state/1)*M+cdf+state%1  =  state+q*(M-1)+cdf
			//but  q' = state-1
			//so  state = state+(state-1+1)*(M-1)+cdf  =  state+q'*(M-1)+(cdf+M-1)
			info->sh=0;
			info->invf=0xFFFFFFFF;
			info->cdf+=(1<<PROBBITS)-1;
		}
		else
		{
			info->sh=31-_lzcnt_u32(freq);//eg: x/2 = x*0x80000000>>32>>0
			unsigned long long inv=((0x100000000ULL<<info->sh)+freq-1)/freq;
			info->invf=(unsigned)inv;
			if(inv>0xFFFFFFFF)
			{
				--info->sh;
				info->invf=(unsigned)(inv>>1);
			}
		}
	}
}
static void enc_packhist(BitPackerLIFO *ec, const int *hist, unsigned long long ctxmask, int ctxidx)//histogram must be normalized to PROBBITS, with spike at 128
{
	if(ctxmask>>ctxidx&1)
		return;
	int sum=0;
	uint16_t CDF[257];
	for(int ks=0;ks<256;++ks)//integrage to zigzag CDF to be packed backwards
	{
		int sym=((ks>>1^-(ks&1))+128)&255;
		int freq=hist[sym];
		CDF[ks]=sum;//separate buffer for faster access in 2nd loop
		sum+=freq;
	}
	CDF[256]=1<<PROBBITS;
	
	int cdfW=CDF[0];
	int CDFlevels=1<<PROBBITS;
	int startsym=0;
	for(int ks=1;ks<=256;++ks)//push GR.k
	{
		int next=CDF[ks], freq=next-cdfW;
		int nbypass=31-_lzcnt_u32(CDFlevels);
		if(ks>1)
			nbypass-=7;
		if(nbypass<0)
			nbypass=0;
		CDF[ks]=nbypass<<PROBBITS|freq;
		cdfW=next;
		CDFlevels-=freq;
		startsym=ks;
		if(!CDFlevels)
			break;
	}
	for(int ks=startsym;ks>0;--ks)//encode GR
	{
		int freq=CDF[ks], nbypass=freq>>PROBBITS;
		freq&=(1<<PROBBITS)-1;
		int nzeros=freq>>nbypass, bypass=freq&((1<<nbypass)-1);
		if(nbypass)
			bitpacker_enc(ec, nbypass, bypass);
		bitpacker_enc(ec, 1, 1);
		while(nzeros)
		{
			bitpacker_enc(ec, 1, 0);
			--nzeros;
		}
#ifdef ANS_VAL
		ansval_push(&ks, sizeof(ks), 1);
#endif
	}
}
static void dec_unpackhist(BitPackerLIFO *ec, unsigned *CDF2sym, unsigned long long ctxmask, int ctxidx)
{
	uint16_t hist[257];
	if(ctxmask>>ctxidx&1)//rare context
	{
		for(int ks=0;ks<256;++ks)//bypass
			hist[ks]=(1<<PROBBITS)/256;
	}
	else
	{
		uint16_t CDF[257]={0};
		int CDFlevels=1<<PROBBITS;
		CDF[0]=0;
		for(int ks=0;ks<256;++ks)//decode GR
		{
			int freq=-1;//stop bit doesn't count
			int nbypass=31-_lzcnt_u32(CDFlevels);
			int ks2=ks+1;
			if(ks2>1)
				nbypass-=7;
			if(nbypass<0)
				nbypass=0;
#ifdef ANS_VAL
			ansval_check(&ks2, sizeof(ks2), 1);
#endif
			int bit=0;
			do
			{
				bit=bitpacker_dec(ec, 1);
				++freq;
			}while(!bit);
			if(nbypass)
				freq=freq<<nbypass|bitpacker_dec(ec, nbypass);

			CDF[ks]=freq;
			CDFlevels-=freq;
			if(CDFlevels<=0)
			{
#ifdef _DEBUG
				if(CDFlevels<0)
					CRASH("CDF unpack error");
#endif
				break;
			}
		}
		if(CDFlevels)
			CRASH("CDF unpack error");
		for(int ks=0;ks<256;++ks)//undo zigzag
		{
			int sym=((ks>>1^-(ks&1))+128)&255;
			hist[sym]=CDF[ks];
		}
	}
	int sum=0;
	for(int ks=0;ks<256;++ks)//integrate
	{
		int freq=hist[ks];
		hist[ks]=sum;
		sum+=freq;
	}
	hist[256]=1<<PROBBITS;
	for(int ks=0;ks<256;++ks)//CDF2sym contains {freq, (state&0xFFF)-cdf, sym}
	{
		int cdf=hist[ks], next=hist[ks+1], freq=next-cdf;
		int val=(freq<<PROBBITS|0)<<8|ks;
		for(int ks2=cdf;ks2<next;++ks2, val+=1<<8)
			CDF2sym[ks2]=val;
	}
}
AWM_INLINE void gather32(int *dst, const int *src, const int *offsets)
{
#ifdef EMULATE_GATHER
	volatile int *ptr=dst;
	ptr[0]=src[offsets[0]];
	ptr[1]=src[offsets[1]];
	ptr[2]=src[offsets[2]];
	ptr[3]=src[offsets[3]];
	ptr[4]=src[offsets[4]];
	ptr[5]=src[offsets[5]];
	ptr[6]=src[offsets[6]];
	ptr[7]=src[offsets[7]];
#else
	_mm256_store_si256((__m256i*)dst, _mm256_i32gather_epi32(src, _mm256_load_si256((__m256i*)offsets), sizeof(int)));
#endif
}
AWM_INLINE void dec_yuv(__m256i *mstate, int kc, const __m256i *ctx0, const int *CDF2syms, unsigned char **pstreamptr, const unsigned char *streamend, __m256i *syms)
{
	const unsigned char *streamptr=*pstreamptr;
	__m256i decctx[2];
	{
		decctx[1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(*ctx0, 1));
		decctx[0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(*ctx0));
	}
	decctx[0]=_mm256_slli_epi32(decctx[0], PROBBITS);
	decctx[1]=_mm256_slli_epi32(decctx[1], PROBBITS);
#ifdef ANS_VAL
	ansval_check(mstate, sizeof(int), NLANES);
#endif
	{
		__m256i mprobmask=_mm256_set1_epi32((1<<PROBBITS)-1);
		__m256i rem0=_mm256_and_si256(mstate[0], mprobmask);
		__m256i rem1=_mm256_and_si256(mstate[1], mprobmask);
		decctx[0]=_mm256_or_si256(decctx[0], rem0);
		decctx[1]=_mm256_or_si256(decctx[1], rem1);
	}
#ifdef ANS_VAL
	ALIGN(32) int debugctx[NLANES];
	memcpy(debugctx, decctx, sizeof(int[NLANES]));
#endif
	const int *statsptr=CDF2syms+((ptrdiff_t)NCTX*kc<<PROBBITS);
	gather32((int*)(decctx+0), statsptr, (int*)(decctx+0));
	gather32((int*)(decctx+1), statsptr, (int*)(decctx+1));

	//update		state = (state>>12)*freq+(rem-cdf)	rem-cdf is prebaked
	{
		__m256i mfreq0=_mm256_srli_epi32(decctx[0], PROBBITS+8);//1 <= freq <= 0xF01
		__m256i mfreq1=_mm256_srli_epi32(decctx[1], PROBBITS+8);
#ifdef ANS_VAL
		__m256i mdebugfreq[1];
		mdebugfreq[0]=_mm256_packus_epi32(mfreq0, mfreq1);
		mdebugfreq[0]=_mm256_permute4x64_epi64(mdebugfreq[0], _MM_SHUFFLE(3, 1, 2, 0));
		ansval_check(mdebugfreq, sizeof(int16_t), NLANES);
#endif
		mstate[0]=_mm256_srli_epi32(mstate[0], PROBBITS);
		mstate[1]=_mm256_srli_epi32(mstate[1], PROBBITS);
		mstate[0]=_mm256_mullo_epi32(mstate[0], mfreq0);//10 cycles
		mstate[1]=_mm256_mullo_epi32(mstate[1], mfreq1);
	}
	{
		__m256i mbias0=_mm256_slli_epi32(decctx[0], PROBBITS);
		__m256i mbias1=_mm256_slli_epi32(decctx[1], PROBBITS);
		mbias0=_mm256_srli_epi32(mbias0, 32-PROBBITS);
		mbias1=_mm256_srli_epi32(mbias1, 32-PROBBITS);
		mstate[0]=_mm256_add_epi32(mstate[0], mbias0);
		mstate[1]=_mm256_add_epi32(mstate[1], mbias1);
	}
	__m256i symmask=_mm256_set1_epi32(255);
	decctx[0]=_mm256_and_si256(decctx[0], symmask);
	decctx[1]=_mm256_and_si256(decctx[1], symmask);
	decctx[0]=_mm256_packus_epi16(decctx[0], decctx[1]);
	syms[0]=_mm256_permute4x64_epi64(decctx[0], _MM_SHUFFLE(3, 1, 2, 0));
#ifdef ANS_VAL
	ansval_check(mstate, sizeof(int), NLANES);
#endif
	//renorm
	{
		__m256i cond0, idx0, lo0, renorm0;
		__m256i cond1, idx1, lo1, renorm1;
		int mask0, mask1;
		__m256i smin=_mm256_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
		cond0=_mm256_cmpgt_epi32(smin, mstate[0]);//signed comparison
		cond1=_mm256_cmpgt_epi32(smin, mstate[1]);
		mask0=_mm256_movemask_ps(_mm256_castsi256_ps(cond0));
		mask1=_mm256_movemask_ps(_mm256_castsi256_ps(cond1));
		idx0=_mm256_load_si256((const __m256i*)ans_permute[mask0]);
		idx1=_mm256_load_si256((const __m256i*)ans_permute[mask1]);
		mask0=_mm_popcnt_u32(mask0);
		mask1=_mm_popcnt_u32(mask1);
#ifdef _MSC_VER
		if(streamptr+sizeof(int16_t)*((ptrdiff_t)mask0+mask1)>streamend)
			CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		lo0=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr)); streamptr+=mask0*sizeof(int16_t);
		lo1=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr)); streamptr+=mask1*sizeof(int16_t);
		renorm0=_mm256_slli_epi32(mstate[0], 16);
		renorm1=_mm256_slli_epi32(mstate[1], 16);
		lo0=_mm256_permutevar8x32_epi32(lo0, idx0);
		lo1=_mm256_permutevar8x32_epi32(lo1, idx1);
		renorm0=_mm256_or_si256(renorm0, lo0);
		renorm1=_mm256_or_si256(renorm1, lo1);
		mstate[0]=_mm256_blendv_epi8(mstate[0], renorm0, cond0);
		mstate[1]=_mm256_blendv_epi8(mstate[1], renorm1, cond1);
	}
	*pstreamptr=(unsigned char*)(size_t)streamptr;
}
static void decorr1d(unsigned char *data, int count, int bytestride, int bestrct, int *rhist)
{
	const unsigned char *combination=rct_combinations[bestrct];
	int yidx=combination[II_PERM_Y];
	int uidx=combination[II_PERM_U];
	int vidx=combination[II_PERM_V];
	int ufromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	int vc0=combination[II_COEFF_V_SUB_Y];
	int vc1=combination[II_COEFF_V_SUB_U];

	unsigned char *ptr=data;
	int prevy=0, prevu=0, prevv=0, offset=0;
	for(int k=0;k<count;++k)
	{
		//if(k==560)//
		//	printf("");

		int y=ptr[yidx]-128;
		int u=ptr[uidx]-128;
		int v=ptr[vidx]-128;
		int sym;
		ptr[0]=sym=(unsigned char)(y-prevy+128);
		++rhist[256*0+sym];
		prevy=y;

		offset=y&ufromy;
		prevu+=offset;
		CLAMP2(prevu, -128, 127);
		ptr[1]=sym=(unsigned char)(u-prevu+128);
		++rhist[256*1+sym];
		prevu=u-offset;

		offset=vc0*y+vc1*u;
		int vpred=(prevv+offset)>>2;
		CLAMP2(vpred, -128, 127);
		ptr[2]=sym=(unsigned char)(v-vpred+128);
		++rhist[256*2+sym];
		prevv=4*v-offset;
		ptr+=bytestride;
	}
}
static void encode1d(unsigned char *data, int count, int bytestride, unsigned *pstate, unsigned char **pstreamptr, const unsigned char *streamend, const rANS_SIMD_SymInfo *rsyminfo)
{
	unsigned char *streamptr=*pstreamptr;
	unsigned state=*pstate;
	unsigned char *ptr=data+(count-(ptrdiff_t)1)*bytestride;
	const rANS_SIMD_SymInfo *info=0;
	for(int k=0;k<count;++k)
	{
		//if(k==count-1)//
		//	printf("");

		info=rsyminfo+ptr[2]+256*2;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)//"streamend" is buffer start
				CRASH("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(unsigned short*)streamptr=(unsigned short)state;
			state>>=RANS_RENORM_BITS;
		}
		state+=((unsigned long long)state*info->invf>>32>>info->sh)*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif

		info=rsyminfo+ptr[1]+256*1;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)
				CRASH("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(unsigned short*)streamptr=(unsigned short)state;
			state>>=RANS_RENORM_BITS;
		}
		state+=((unsigned long long)state*info->invf>>32>>info->sh)*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif

		info=rsyminfo+ptr[0]+256*0;
		if(state>info->smax)
		{
			streamptr-=2;
#ifdef _DEBUG
			if(streamptr<=streamend)
				CRASH("OOB ptr %016zX <= %016zX", streamptr, streamend);
#endif
			*(unsigned short*)streamptr=(unsigned short)state;
			state>>=RANS_RENORM_BITS;
		}
		state+=((unsigned long long)state*info->invf>>32>>info->sh)*info->negf+info->cdf;
#ifdef ANS_VAL
		ansval_push(&state, sizeof(state), 1);
#endif
		ptr-=bytestride;
	}
	*pstreamptr=streamptr;
	*pstate=state;
}
static void decode1d(unsigned char *data, int count, int bytestride, int bestrct, unsigned *pstate, const unsigned char **pstreamptr, const unsigned char *streamend, unsigned *rCDF2syms)
{
	const unsigned char *combination=rct_combinations[bestrct];
	int yidx=combination[II_PERM_Y];
	int uidx=combination[II_PERM_U];
	int vidx=combination[II_PERM_V];
	int ufromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	int vc0=combination[II_COEFF_V_SUB_Y];
	int vc1=combination[II_COEFF_V_SUB_U];

	const unsigned char *streamptr=*pstreamptr;
	unsigned state=*pstate;
	unsigned char *ptr=data;
	int prevy=0, prevu=0, prevv=0, offset=0;
	int y=0, u=0, v=0;
	for(int k=0;k<count;++k)
	{
		unsigned info;

		//yuv = (char)(error+N-128)
		info=rCDF2syms[0<<PROBBITS|(state&((1<<PROBBITS)-1))];
		y=(char)(info+prevy-128);
		prevy=y;
#ifdef ANS_VAL
		ansval_check(&state, sizeof(state), 1);
#endif
		state=(state>>PROBBITS)*(info>>(PROBBITS+8))+(info<<PROBBITS>>(32-PROBBITS));
		if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
		{
#ifdef _DEBUG
			if(streamptr+2>streamend)
				CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(unsigned short*)streamptr;
			streamptr+=2;
		}

		//if(k==560)//
		//	printf("");//

		offset=y&ufromy;
		prevu+=offset;
		CLAMP2(prevu, -128, 127);
		info=rCDF2syms[1<<PROBBITS|(state&((1<<PROBBITS)-1))];
		u=(char)(info+prevu-128);
		prevu=u-offset;
#ifdef ANS_VAL
		ansval_check(&state, sizeof(state), 1);
#endif
		state=(state>>PROBBITS)*(info>>(PROBBITS+8))+(info<<PROBBITS>>(32-PROBBITS));
		if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
		{
#ifdef _DEBUG
			if(streamptr+2>streamend)
				CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(unsigned short*)streamptr;
			streamptr+=2;
		}

		offset=vc0*y+vc1*u;
		int vpred=(prevv+offset)>>2;
		CLAMP2(vpred, -128, 127);
		info=rCDF2syms[2<<PROBBITS|(state&((1<<PROBBITS)-1))];
		v=(char)(info+vpred-128);
		prevv=4*v-offset;
#ifdef ANS_VAL
		ansval_check(&state, sizeof(state), 1);
#endif
		state=(state>>PROBBITS)*(info>>(PROBBITS+8))+(info<<PROBBITS>>(32-PROBBITS));
		if(state<(1<<(RANS_STATE_BITS-RANS_RENORM_BITS)))
		{
#ifdef _DEBUG
			if(streamptr+2>streamend)
				CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(unsigned short*)streamptr;
			streamptr+=2;
		}
		ptr[yidx]=y+128;
		ptr[uidx]=u+128;
		ptr[vidx]=v+128;
		ptr+=bytestride;
	}
	*pstreamptr=streamptr;
	*pstate=state;
}

typedef struct _ThreadArgs
{
	int threadidx, b1, b2;
	uint8_t *chunk;
	ptrdiff_t clen;
	void *buffer;
} ThreadArgs;
static THREAD_RET THREAD_CALL c34_thread(void *param)
{
	ThreadArgs *arg;
	int blockwidth;
	ptrdiff_t cap, blocksize;
	uint16_t *ctxbuf;
	uint8_t *block;
	int cbufsize;
	int16_t *cbuf;
	int kb;
	int paddedwidth;
	__m256i dist_rcp, mdist;
	const int hsize=(int)sizeof(int[3*NCTX<<8]);
	const int encstatsize=(int)sizeof(rANS_SIMD_SymInfo[3*NCTX<<8]);
	const int CDF2syms_size=(int)sizeof(int[3*NCTX<<PROBBITS]);
	int *hists;
	rANS_SIMD_SymInfo *encstats;
	uint32_t *CDF2syms;
	uint8_t *streamstart=0, *streamend=0, *streamptr=0;

	(void)streamstart;
	arg=(ThreadArgs*)param;
	blocksize=(ptrdiff_t)3*YMAXBLOCK*XMAXBLOCK;
	if(blocksize>usize)
		blocksize=usize;
	cap=sizeof(uint8_t[5*YMAXBLOCK*XMAXBLOCK])*((ptrdiff_t)arg->b2-arg->b1);
	if(cap>(int)sizeof(uint8_t[5])*res)
		cap=sizeof(uint8_t[5])*res;
	if(cap<sizeof(uint8_t[8*YMAXBLOCK*XMAXBLOCK]))
		cap=sizeof(uint8_t[8*YMAXBLOCK*XMAXBLOCK]);
	if(!fwd)
		cap=blocksize;
	arg->buffer=malloc(cap);
	blockwidth=XMAXBLOCK;
	if(blockwidth>iw)
		blockwidth=iw;
	paddedwidth=blockwidth+16;
	cbufsize=(int)sizeof(int16_t[4*6*NLANES])*paddedwidth;//4 padded rows  *  {Y*NLANES, U*NLANES, V*NLANES,  eY*NLANES, eU*NLANES, eV*NLANES}
	cbuf=(int16_t*)_mm_malloc(cbufsize, sizeof(__m256i));

	if(fwd)
	{
		hists=(int*)malloc(hsize);//enc-only
		encstats=(rANS_SIMD_SymInfo*)_mm_malloc(encstatsize, sizeof(__m256i));
		CDF2syms=0;

		ctxbuf=(uint16_t*)arg->buffer;
		block=(uint8_t*)arg->buffer+blocksize;
		streamstart=(uint8_t*)arg->buffer;
		streamptr=streamend=(uint8_t*)arg->buffer+cap;
	}
	else
	{
		hists=0;
		encstats=0;
		CDF2syms=(uint32_t*)_mm_malloc(CDF2syms_size, sizeof(__m256i));//dec-only
		block=(uint8_t*)arg->buffer;
		ctxbuf=0;
	}
	if(!arg->buffer||!cbuf||(fwd?!hists||!encstats:!CDF2syms))
	{
		CRASH("Alloc error");
		return 1;
	}
	dist_rcp=_mm256_set1_epi16(0x7FFF);
	mdist=_mm256_set1_epi16(1);
	if(dist>1)
	{
		dist_rcp=_mm256_set1_epi16(((1<<16)+dist-1)/dist);//x/dist  ->  {x*=inv; x=(x>>16)+((unsigned)x>>31);}
		mdist=_mm256_set1_epi16(dist);
	}
	for(kb=fwd?arg->b2-1:arg->b1;fwd?kb>=arg->b1:kb<arg->b2;kb+=fwd?-1:1)//block loop
	{
		int bx=kb%xnblocks, by=kb/xnblocks;
		int x1=blockxstops[bx], x2=blockxstops[bx+1], blockw=x2-x1, iblockw=blockw/XLANES;
		int y1=blockystops[by], y2=blockystops[by+1], blockh=y2-y1, iblockh=blockh/YLANES;
		int bestrct=0;
		uint8_t *streamptr0=streamptr;
		uint64_t ctxmask=0;
		__m256i mstate[2];

		if(blockw%XLANES||blockh%YLANES)
		{
			CRASH("Invalid block size  WH %d*%d  (expected multiples of 4)", blockw, blockh);
			return 1;
		}

		//interleave
		{
			const uint8_t *slowptrs0[NLANES], *slowptrs[NLANES];
			uint8_t *fastptr;
			int kx, ky, k;

			fastptr=block;
			for(ky=0;ky<YLANES;++ky)//spread slow pointers
			{
				for(kx=0;kx<XLANES;++kx)
					slowptrs0[XLANES*ky+kx]=image+3*(iw*(y1+iblockh*ky)+x1+iblockw*kx);
			}
			for(ky=0;ky<iblockh;++ky)//interleave
			{
				memcpy((void*)slowptrs, slowptrs0, sizeof(slowptrs));
				for(kx=0;kx<iblockw;++kx)
				{
					for(k=0;k<NLANES;++k)
						*fastptr++=*slowptrs[k]++;
					for(k=0;k<NLANES;++k)
						*fastptr++=*slowptrs[k]++;
					for(k=0;k<NLANES;++k)
						*fastptr++=*slowptrs[k]++;
				}
				for(k=0;k<NLANES;++k)
					slowptrs0[k]+=rowstride;
			}
		}

		//analysis
		{
			__m256i mcounters[OCH_COUNT];
			__m128i half8=_mm_set1_epi8(-128);
			__m256i wordmask=_mm256_set1_epi64x(0xFFFF);
			int kx, ky, k, kt;
			uint8_t *imptr;
			long long minerr;
			int64_t counters[OCH_COUNT];

			imptr=block;
			memset(mcounters, 0, sizeof(mcounters));
			for(ky=0;ky<iblockh;++ky)
			{
				__m256i prev[OCH_COUNT];//16-bit
				memset(prev, 0, sizeof(prev));
				for(kx=0;kx<iblockw;++kx)
				{
					__m256i r=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+0), half8));
					__m256i g=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+1), half8));
					__m256i b=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+2), half8));
					imptr+=3*NLANES;
					r=_mm256_slli_epi16(r, 2);
					g=_mm256_slli_epi16(g, 2);
					b=_mm256_slli_epi16(b, 2);
					__m256i rg=_mm256_sub_epi16(r, g);
					__m256i gb=_mm256_sub_epi16(g, b);
					__m256i br=_mm256_sub_epi16(b, r);
					__m256i t0, t1, t2;
#define UPDATE(IDXA, IDXB, IDXC, A0, B0, C0)\
	do\
	{\
		__m256i ta=_mm256_sub_epi16(A0, prev[IDXA]);\
		__m256i tb=_mm256_sub_epi16(B0, prev[IDXB]);\
		__m256i tc=_mm256_sub_epi16(C0, prev[IDXC]);\
		prev[IDXA]=A0;\
		prev[IDXB]=B0;\
		prev[IDXC]=C0;\
		ta=_mm256_abs_epi16(ta);\
		tb=_mm256_abs_epi16(tb);\
		tc=_mm256_abs_epi16(tc);\
		ta=_mm256_add_epi16(ta, _mm256_srli_epi64(ta, 32));\
		tb=_mm256_add_epi16(tb, _mm256_srli_epi64(tb, 32));\
		tc=_mm256_add_epi16(tc, _mm256_srli_epi64(tc, 32));\
		ta=_mm256_add_epi16(ta, _mm256_srli_epi64(ta, 16));\
		tb=_mm256_add_epi16(tb, _mm256_srli_epi64(tb, 16));\
		tc=_mm256_add_epi16(tc, _mm256_srli_epi64(tc, 16));\
		mcounters[IDXA]=_mm256_add_epi64(mcounters[IDXA], _mm256_and_si256(ta, wordmask));\
		mcounters[IDXB]=_mm256_add_epi64(mcounters[IDXB], _mm256_and_si256(tb, wordmask));\
		mcounters[IDXC]=_mm256_add_epi64(mcounters[IDXC], _mm256_and_si256(tc, wordmask));\
	}while(0)
					UPDATE(OCH_Y400, OCH_Y040, OCH_Y004, r, g, b);
					UPDATE(OCH_CX40, OCH_C0X4, OCH_C40X, rg, gb, br);
					t0=_mm256_add_epi16(rg, _mm256_srai_epi16(gb, 2));//r-(3*g+b)/4 = r-g-(b-g)/4
					t1=_mm256_add_epi16(rg, _mm256_srai_epi16(br, 2));//g-(3*r+b)/4 = g-r-(b-r)/4
					t2=_mm256_add_epi16(br, _mm256_srai_epi16(rg, 2));//b-(3*r+g)/4 = b-r-(g-r)/4
					UPDATE(OCH_CX31, OCH_C3X1, OCH_C31X, t0, t1, t2);
					t0=_mm256_add_epi16(br, _mm256_srai_epi16(gb, 2));//r-(g+3*b)/4 = r-b-(g-b)/4
					t1=_mm256_add_epi16(gb, _mm256_srai_epi16(br, 2));//g-(r+3*b)/4 = g-b-(r-b)/4
					t2=_mm256_add_epi16(gb, _mm256_srai_epi16(rg, 2));//b-(r+3*g)/4 = b-g-(r-g)/4
					UPDATE(OCH_CX13, OCH_C1X3, OCH_C13X, t0, t1, t2);
					t0=_mm256_srai_epi16(_mm256_sub_epi16(rg, br), 1);//r-(g+b)/2 = (r-g + r-b)/2
					t1=_mm256_srai_epi16(_mm256_sub_epi16(gb, rg), 1);//g-(r+b)/2 = (g-r + g-b)/2
					t2=_mm256_srai_epi16(_mm256_sub_epi16(br, gb), 1);//b-(r+g)/2 = (b-r + b-g)/2
					UPDATE(OCH_CX22, OCH_C2X2, OCH_C22X, t0, t1, t2);
				}
			}
			for(k=0;k<OCH_COUNT;++k)
			{
				ALIGN(32) long long temp[4]={0};
				_mm256_store_si256((__m256i*)temp, mcounters[k]);
				counters[k]=temp[0]+temp[1]+temp[2]+temp[3];
			}
			for(kt=0, minerr=0;kt<RCT_COUNT;++kt)
			{
				const unsigned char *rct=rct_combinations[kt];
				long long currerr=
					+counters[rct[0]]
					+counters[rct[1]]
					+counters[rct[2]]
				;
//#ifdef LOUD
//				printf("%-14s %12lld + %12lld + %12lld = %12lld%s\n",
//					rct_names[kt],
//					counters[rct[0]],
//					counters[rct[1]],
//					counters[rct[2]],
//					currerr,
//					!kt||minerr>currerr?" <-":""
//				);
//#endif
				if(!kt||minerr>currerr)
				{
					minerr=currerr;
					bestrct=kt;
				}
			}
//#ifdef LOUD
//			printf("%s  WH %d*%d = %td bytes\n", rct_names[bestrct], iw, ih, usize);
//#endif
		}

		//main
		{
			const uint8_t *combination;
			int yidx, uidx, vidx;
			__m256i uhelpmask;
			__m256i vc0;
			__m256i vc1;
			__m256i mctxmax=_mm256_set1_epi16(NCTX-1);
			__m256i mctxuoffset=_mm256_set1_epi16(NCTX);
			__m256i mctxvoffset=_mm256_set1_epi16(NCTX*2);
			__m256i amin=_mm256_set1_epi16(-128);
			__m256i amax=_mm256_set1_epi16(127);
			__m128i half8=_mm_set1_epi8(-128);
			__m256i bytemask=_mm256_set1_epi16(255);
			__m256i wordmask=_mm256_set1_epi32(0xFFFF);
			__m256i myuv[3];
			uint16_t *ctxptr=ctxbuf;
			uint8_t *imptr=block;

			memset(cbuf, 0, cbufsize);
			memset(myuv, 0, sizeof(myuv));
			if(fwd)
				memset(hists, 0, hsize);
			else
			{
				streamptr=stream+sizeof(int)*nblocks;
				for(int kb0=0;kb0<kb;++kb0)
					streamptr+=bsizes[kb0];
				streamstart=streamptr;
				streamend=streamptr+bsizes[kb];

				bestrct=*streamptr++;
				{
					BitPackerLIFO ec;

					ctxmask=*(uint64_t*)streamptr;
					streamptr+=8;
					bitpacker_dec_init(&ec, streamptr, streamend);
					for(int kc=0;kc<3*NCTX;++kc)
						dec_unpackhist(&ec, CDF2syms+((ptrdiff_t)kc<<PROBBITS), ctxmask, kc);
					streamptr=(unsigned char*)(size_t)ec.srcfwdptr;
				}
#ifdef _DEBUG
				if(streamptr>streamend)
					CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
				memcpy(mstate, streamptr, sizeof(mstate));
				streamptr+=sizeof(mstate);
			}
			combination=rct_combinations[bestrct];
			yidx=combination[II_PERM_Y]*NLANES;
			uidx=combination[II_PERM_U]*NLANES;
			vidx=combination[II_PERM_V]*NLANES;
			uhelpmask=_mm256_set1_epi16(-(combination[II_COEFF_U_SUB_Y]!=0));
			vc0=_mm256_set1_epi16(combination[II_COEFF_V_SUB_Y]);
			vc1=_mm256_set1_epi16(combination[II_COEFF_V_SUB_U]);
			for(int ky=0;ky<iblockh;++ky)//main loop
			{
				ALIGN(32) int16_t *rows[]=
				{
					cbuf+(paddedwidth*((ky-0LL)&3)+8LL)*6*NLANES,
					cbuf+(paddedwidth*((ky-1LL)&3)+8LL)*6*NLANES,
					cbuf+(paddedwidth*((ky-2LL)&3)+8LL)*6*NLANES,
					cbuf+(paddedwidth*((ky-3LL)&3)+8LL)*6*NLANES,
				};
				ALIGN(32) uint16_t syms[3*NLANES]={0};
				__m256i NW[3], N[3], W[3];
				__m256i eW[3], ecurr[3], eNEE[3], eNEEE[3];
				memset(NW, 0, sizeof(NW));
				memset(N, 0, sizeof(N));
				memset(W, 0, sizeof(W));
				memset(eW, 0, sizeof(eW));
				memset(ecurr, 0, sizeof(ecurr));
				memset(eNEE, 0, sizeof(eNEE));
				memset(eNEEE, 0, sizeof(eNEEE));
				//                       (__m256i*)rows[-Y]+E+C+X*6
				eNEE[0]=_mm256_load_si256((__m256i*)rows[1]+3+0+2*6);
				eNEE[1]=_mm256_load_si256((__m256i*)rows[1]+3+1+2*6);
				eNEE[2]=_mm256_load_si256((__m256i*)rows[1]+3+2+2*6);
				for(int kx=0;kx<iblockw;++kx)
				{
					__m256i
						predY, ctxY,
						predU, ctxU,
						predV, ctxV;
					//__m256i predYUV0[3];
					__m256i msyms, moffset;

					//                     (__m256i*)rows[1]+E+C+X*6
					N[0]=_mm256_load_si256((__m256i*)rows[1]+0+0+0*6);//y neighbors
					N[1]=_mm256_load_si256((__m256i*)rows[1]+0+1+0*6);//u
					N[2]=_mm256_load_si256((__m256i*)rows[1]+0+2+0*6);//v

					//if(kb==1&&ky==35&&kx==84)//
					//if(kb==1&&ky==0&&kx==1)//
					//if(kb==1&&ky==35&&kx==83)//
					//if(kb==0&&ky==0&&kx==2)//
					//	printf("");

					//context
					{
						//context = FLOOR_LOG2(eW*eW+1)
						__m256i one=_mm256_set1_epi32(1);
						__m256i cy0=_mm256_and_si256(eW[0], wordmask), cy1=_mm256_srli_epi32(eW[0], 16);
						__m256i cu0=_mm256_and_si256(eW[1], wordmask), cu1=_mm256_srli_epi32(eW[1], 16);
						__m256i cv0=_mm256_and_si256(eW[2], wordmask), cv1=_mm256_srli_epi32(eW[2], 16);
						cy0=_mm256_mullo_epi32(cy0, cy0);
						cy1=_mm256_mullo_epi32(cy1, cy1);
						cu0=_mm256_mullo_epi32(cu0, cu0);
						cu1=_mm256_mullo_epi32(cu1, cu1);
						cv0=_mm256_mullo_epi32(cv0, cv0);
						cv1=_mm256_mullo_epi32(cv1, cv1);
						cy0=_mm256_add_epi32(cy0, one);
						cy1=_mm256_add_epi32(cy1, one);
						cu0=_mm256_add_epi32(cu0, one);
						cu1=_mm256_add_epi32(cu1, one);
						cv0=_mm256_add_epi32(cv0, one);
						cv1=_mm256_add_epi32(cv1, one);
						//FLOOR_LOG2_32x8(X) = _mm256_sub_epi32(_mm256_srli_epi32(_mm256_castps_si256(_mm256_cvtepi32_ps(X)), 23), _mm256_set1_epi32(127))
						cy0=_mm256_castps_si256(_mm256_cvtepi32_ps(cy0));
						cy1=_mm256_castps_si256(_mm256_cvtepi32_ps(cy1));
						cu0=_mm256_castps_si256(_mm256_cvtepi32_ps(cu0));
						cu1=_mm256_castps_si256(_mm256_cvtepi32_ps(cu1));
						cv0=_mm256_castps_si256(_mm256_cvtepi32_ps(cv0));
						cv1=_mm256_castps_si256(_mm256_cvtepi32_ps(cv1));
						cy0=_mm256_srli_epi32(cy0, 23);
						cy1=_mm256_srli_epi32(cy1, 23);
						cu0=_mm256_srli_epi32(cu0, 23);
						cu1=_mm256_srli_epi32(cu1, 23);
						cv0=_mm256_srli_epi32(cv0, 23);
						cv1=_mm256_srli_epi32(cv1, 23);
						__m256i expbias=_mm256_set1_epi32(127);
						cy0=_mm256_sub_epi32(cy0, expbias);
						cy1=_mm256_sub_epi32(cy1, expbias);
						cu0=_mm256_sub_epi32(cu0, expbias);
						cu1=_mm256_sub_epi32(cu1, expbias);
						cv0=_mm256_sub_epi32(cv0, expbias);
						cv1=_mm256_sub_epi32(cv1, expbias);
						cy1=_mm256_slli_epi32(cy1, 16);
						cu1=_mm256_slli_epi32(cu1, 16);
						cv1=_mm256_slli_epi32(cv1, 16);
						ctxY=_mm256_or_si256(cy0, cy1);
						ctxU=_mm256_or_si256(cu0, cu1);
						ctxV=_mm256_or_si256(cv0, cv1);
						ctxY=_mm256_min_epi16(ctxY, mctxmax);
						ctxU=_mm256_min_epi16(ctxU, mctxmax);
						ctxV=_mm256_min_epi16(ctxV, mctxmax);
					}

					//predict
					{
#if 0
						const int borderW=3;
						const int borderN=3;
						const int borderE=3;
						int cond_cg=(unsigned)(kx-3*NLANES*borderW)>=(unsigned)((iblockw-1)*3*NLANES*(borderW+borderE))
							||(unsigned)(ky-borderN)>=(unsigned)(iblockh-borderN);
						__m256i mcg[3];
#endif
						__m256i ymin=_mm256_min_epi16(N[0], W[0]);
						__m256i ymax=_mm256_max_epi16(N[0], W[0]);
						__m256i umin=_mm256_min_epi16(N[1], W[1]);
						__m256i umax=_mm256_max_epi16(N[1], W[1]);
						__m256i vmin=_mm256_min_epi16(N[2], W[2]);
						__m256i vmax=_mm256_max_epi16(N[2], W[2]);
						predY=_mm256_add_epi16(N[0], W[0]);//N+W-NW
						predU=_mm256_add_epi16(N[1], W[1]);
						predV=_mm256_add_epi16(N[2], W[2]);
						predY=_mm256_sub_epi16(predY, NW[0]);
						predU=_mm256_sub_epi16(predU, NW[1]);
						predV=_mm256_sub_epi16(predV, NW[2]);
#if 0
						mcg[0]=predY;
						mcg[1]=predU;
						mcg[2]=predV;
						if(effort==1)//predict
						{
							/*
							effort 1
							0	N+W-NW
							1	N
							2	NE
							3	W
							*/

							//N+W-NW
							L1preds[0*3+0]=predY;
							L1preds[0*3+1]=predU;
							L1preds[0*3+2]=predV;

							//N
							L1preds[1*3+0]=N[0];
							L1preds[1*3+1]=N[1];
							L1preds[1*3+2]=N[2];

							//NE
							L1preds[2*3+0]=_mm256_load_si256((__m256i*)rows[1]+0+0+1*6);
							L1preds[2*3+1]=_mm256_load_si256((__m256i*)rows[1]+0+1+1*6);
							L1preds[2*3+2]=_mm256_load_si256((__m256i*)rows[1]+0+2+1*6);

							//W
							L1preds[3*3+0]=W[0];
							L1preds[3*3+1]=W[1];
							L1preds[3*3+2]=W[2];
					

							//mix
							__m256i mp[6], t[6];
							mp[0]=_mm256_setzero_si256();
							mp[1]=_mm256_setzero_si256();
							mp[2]=_mm256_setzero_si256();
							mp[3]=_mm256_setzero_si256();
							mp[4]=_mm256_setzero_si256();
							mp[5]=_mm256_setzero_si256();
		#if defined __GNUC__ && !defined PROFILER
		#pragma GCC unroll 4
		#endif
							for(int k=0;k<L1_NPREDS1;++k)
							{
								//16 -> 32		3 lo 3 hi registers
								t[0]=_mm256_slli_epi32(L1preds[k*3+0], 16);
								t[1]=_mm256_slli_epi32(L1preds[k*3+1], 16);
								t[2]=_mm256_slli_epi32(L1preds[k*3+2], 16);
								t[3]=_mm256_srai_epi32(L1preds[k*3+0], 16);
								t[4]=_mm256_srai_epi32(L1preds[k*3+1], 16);
								t[5]=_mm256_srai_epi32(L1preds[k*3+2], 16);
								t[0]=_mm256_srai_epi32(t[0], 16);
								t[1]=_mm256_srai_epi32(t[1], 16);
								t[2]=_mm256_srai_epi32(t[2], 16);
								t[0]=_mm256_mullo_epi32(t[0], _mm256_load_si256((__m256i*)L1weights+k*6+0));
								t[1]=_mm256_mullo_epi32(t[1], _mm256_load_si256((__m256i*)L1weights+k*6+1));
								t[2]=_mm256_mullo_epi32(t[2], _mm256_load_si256((__m256i*)L1weights+k*6+2));
								t[3]=_mm256_mullo_epi32(t[3], _mm256_load_si256((__m256i*)L1weights+k*6+3));
								t[4]=_mm256_mullo_epi32(t[4], _mm256_load_si256((__m256i*)L1weights+k*6+4));
								t[5]=_mm256_mullo_epi32(t[5], _mm256_load_si256((__m256i*)L1weights+k*6+5));
								mp[0]=_mm256_add_epi32(mp[0], t[0]);
								mp[1]=_mm256_add_epi32(mp[1], t[1]);
								mp[2]=_mm256_add_epi32(mp[2], t[2]);
								mp[3]=_mm256_add_epi32(mp[3], t[3]);
								mp[4]=_mm256_add_epi32(mp[4], t[4]);
								mp[5]=_mm256_add_epi32(mp[5], t[5]);
							}
							//__m256i rcon=_mm256_set1_epi32((1<<L1_SH1)-1);
							//mp[0]=_mm256_add_epi32(mp[0], _mm256_and_si256(rcon, _mm256_srli_epi32(mp[0], 31)));//rounding to zero
							//mp[1]=_mm256_add_epi32(mp[1], _mm256_and_si256(rcon, _mm256_srli_epi32(mp[1], 31)));
							//mp[2]=_mm256_add_epi32(mp[2], _mm256_and_si256(rcon, _mm256_srli_epi32(mp[2], 31)));
							//mp[3]=_mm256_add_epi32(mp[3], _mm256_and_si256(rcon, _mm256_srli_epi32(mp[3], 31)));
							//mp[4]=_mm256_add_epi32(mp[4], _mm256_and_si256(rcon, _mm256_srli_epi32(mp[4], 31)));
							//mp[5]=_mm256_add_epi32(mp[5], _mm256_and_si256(rcon, _mm256_srli_epi32(mp[5], 31)));

							__m256i rcon=_mm256_set1_epi32(1<<L1_SH1>>1);
							mp[0]=_mm256_add_epi32(mp[0], rcon);//rounding to nearest
							mp[1]=_mm256_add_epi32(mp[1], rcon);
							mp[2]=_mm256_add_epi32(mp[2], rcon);
							mp[3]=_mm256_add_epi32(mp[3], rcon);
							mp[4]=_mm256_add_epi32(mp[4], rcon);
							mp[5]=_mm256_add_epi32(mp[5], rcon);

							mp[0]=_mm256_srai_epi32(mp[0], L1_SH1);
							mp[1]=_mm256_srai_epi32(mp[1], L1_SH1);
							mp[2]=_mm256_srai_epi32(mp[2], L1_SH1);
							mp[3]=_mm256_slli_epi32(mp[3], 16-L1_SH1);
							mp[4]=_mm256_slli_epi32(mp[4], 16-L1_SH1);
							mp[5]=_mm256_slli_epi32(mp[5], 16-L1_SH1);
							//32 -> 16
							predY=_mm256_blend_epi16(mp[0], mp[3], 0xAA);
							predU=_mm256_blend_epi16(mp[1], mp[4], 0xAA);
							predV=_mm256_blend_epi16(mp[2], mp[5], 0xAA);


							//loosen pred range
							if(!cond_cg)
							{
								t[0]=_mm256_load_si256((__m256i*)rows[1]+0+0+1*6);//NE
								t[1]=_mm256_load_si256((__m256i*)rows[1]+0+1+1*6);
								t[2]=_mm256_load_si256((__m256i*)rows[1]+0+2+1*6);
								ymin=_mm256_min_epi16(ymin, t[0]);
								ymax=_mm256_max_epi16(ymax, t[0]);
								umin=_mm256_min_epi16(umin, t[1]);
								umax=_mm256_max_epi16(umax, t[1]);
								vmin=_mm256_min_epi16(vmin, t[2]);
								vmax=_mm256_max_epi16(vmax, t[2]);
								t[0]=_mm256_load_si256((__m256i*)rows[1]+0+0+3*6);//NEEE
								t[1]=_mm256_load_si256((__m256i*)rows[1]+0+1+3*6);
								t[2]=_mm256_load_si256((__m256i*)rows[1]+0+2+3*6);
								ymin=_mm256_min_epi16(ymin, t[0]);
								ymax=_mm256_max_epi16(ymax, t[0]);
								umin=_mm256_min_epi16(umin, t[1]);
								umax=_mm256_max_epi16(umax, t[1]);
								vmin=_mm256_min_epi16(vmin, t[2]);
								vmax=_mm256_max_epi16(vmax, t[2]);
							}
						}
						else if(effort==2)
						{
							__m256i cache[3];
							/*
							effort 2
							0	N
							1	W
							2	3*(N-NN)+NNN
							3	3*(W-WW)+WWW
							4	W+NE-N
							5	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4
							6	N+W-NW
							7	N+NE-NNE
							*/

							//N
							L1preds[0*3+0]=N[0];
							L1preds[0*3+1]=N[1];
							L1preds[0*3+2]=N[2];

							//W
							L1preds[1*3+0]=W[0];
							L1preds[1*3+1]=W[1];
							L1preds[1*3+2]=W[2];

							//3*(N-NN)+NNN
							cache[0]=_mm256_sub_epi16(N[0], _mm256_load_si256((__m256i*)rows[2]+0+0+0*6));//N-NN
							cache[1]=_mm256_sub_epi16(N[1], _mm256_load_si256((__m256i*)rows[2]+0+1+0*6));
							cache[2]=_mm256_sub_epi16(N[2], _mm256_load_si256((__m256i*)rows[2]+0+2+0*6));
							cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
							cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
							L1preds[2*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[3]+0+0+0*6));//+NNN
							L1preds[2*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[3]+0+1+0*6));
							L1preds[2*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[3]+0+2+0*6));

							//3*(W-WW)+WWW
							cache[0]=_mm256_sub_epi16(W[0], _mm256_load_si256((__m256i*)rows[0]+0+0-2*6));//W-WW
							cache[1]=_mm256_sub_epi16(W[1], _mm256_load_si256((__m256i*)rows[0]+0+1-2*6));
							cache[2]=_mm256_sub_epi16(W[2], _mm256_load_si256((__m256i*)rows[0]+0+2-2*6));
							cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
							cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
							L1preds[3*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[0]+0+0-3*6));//+WWW
							L1preds[3*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[0]+0+1-3*6));
							L1preds[3*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[0]+0+2-3*6));

							//W+NE-N
							cache[0]=_mm256_sub_epi16(W[0], N[0]);
							cache[1]=_mm256_sub_epi16(W[1], N[1]);
							cache[2]=_mm256_sub_epi16(W[2], N[2]);
							L1preds[4*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+0+1*6));//+NE
							L1preds[4*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+1+1*6));
							L1preds[4*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+2+1*6));

							//N+W-NW
							L1preds[5*3+0]=predY;
							L1preds[5*3+1]=predU;
							L1preds[5*3+2]=predV;

							//N+NE-NNE
							cache[0]=_mm256_add_epi16(N[0], _mm256_load_si256((__m256i*)rows[1]+0+0+1*6));//N+NE
							cache[1]=_mm256_add_epi16(N[1], _mm256_load_si256((__m256i*)rows[1]+0+1+1*6));
							cache[2]=_mm256_add_epi16(N[2], _mm256_load_si256((__m256i*)rows[1]+0+2+1*6));
							L1preds[6*3+0]=_mm256_sub_epi16(cache[0], _mm256_load_si256((__m256i*)rows[2]+0+0+1*6));//NNE
							L1preds[6*3+1]=_mm256_sub_epi16(cache[1], _mm256_load_si256((__m256i*)rows[2]+0+1+1*6));
							L1preds[6*3+2]=_mm256_sub_epi16(cache[2], _mm256_load_si256((__m256i*)rows[2]+0+2+1*6));

							//(WWWW+WWW+NNN+NEE+NEEE+NEEEE-(NW+N))>>2
							cache[0]=_mm256_load_si256((__m256i*)rows[0]+0+0-4*6);//WWWW
							cache[1]=_mm256_load_si256((__m256i*)rows[0]+0+1-4*6);
							cache[2]=_mm256_load_si256((__m256i*)rows[0]+0+2-4*6);
							cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[0]+0+0-3*6));//+WWW
							cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[0]+0+1-3*6));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[0]+0+2-3*6));
							cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[3]+0+0+0*6));//+NNN
							cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[3]+0+1+0*6));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[3]+0+2+0*6));
							cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+0+2*6));//+NEE
							cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+1+2*6));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+2+2*6));
							cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+0+3*6));//+NEEE
							cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+1+3*6));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+2+3*6));
							cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+0+4*6));//+NEEEE
							cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+1+4*6));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+2+4*6));
							cache[0]=_mm256_sub_epi16(cache[0], _mm256_add_epi16(N[0], NW[0]));
							cache[1]=_mm256_sub_epi16(cache[1], _mm256_add_epi16(N[1], NW[1]));
							cache[2]=_mm256_sub_epi16(cache[2], _mm256_add_epi16(N[2], NW[2]));
							L1preds[7*3+0]=_mm256_srai_epi16(cache[0], 2);
							L1preds[7*3+1]=_mm256_srai_epi16(cache[1], 2);
							L1preds[7*3+2]=_mm256_srai_epi16(cache[2], 2);


							//mix
							__m256i mp[6], t[6];
							mp[0]=_mm256_setzero_si256();
							mp[1]=_mm256_setzero_si256();
							mp[2]=_mm256_setzero_si256();
							mp[3]=_mm256_setzero_si256();
							mp[4]=_mm256_setzero_si256();
							mp[5]=_mm256_setzero_si256();

							//mp[0]=_mm256_load_si256((__m256i*)L1coeffs+L1_NPREDS3*6+0);//bias
							//mp[1]=_mm256_load_si256((__m256i*)L1coeffs+L1_NPREDS3*6+1);
							//mp[2]=_mm256_load_si256((__m256i*)L1coeffs+L1_NPREDS3*6+2);
							//mp[3]=_mm256_load_si256((__m256i*)L1coeffs+L1_NPREDS3*6+3);
							//mp[4]=_mm256_load_si256((__m256i*)L1coeffs+L1_NPREDS3*6+4);
							//mp[5]=_mm256_load_si256((__m256i*)L1coeffs+L1_NPREDS3*6+5);
		#ifdef PRINT_L1_BOUNDS
							for(int k2=0;k2<6*8;++k2)
							{
								int bias=L1coeffs[L1_NPREDS3*6*8+k2];
								if(bmin>bias)bmin=bias;
								if(bmax<bias)bmax=bias;
							}
		#endif
		#if defined __GNUC__ && !defined PROFILER
		#pragma GCC unroll 8
		#endif
							for(int k=0;k<L1_NPREDS2;++k)
							{
								//16 -> 32		3 lo 3 hi registers
								t[0]=_mm256_slli_epi32(L1preds[k*3+0], 16);
								t[1]=_mm256_slli_epi32(L1preds[k*3+1], 16);
								t[2]=_mm256_slli_epi32(L1preds[k*3+2], 16);
								t[3]=_mm256_srai_epi32(L1preds[k*3+0], 16);
								t[4]=_mm256_srai_epi32(L1preds[k*3+1], 16);
								t[5]=_mm256_srai_epi32(L1preds[k*3+2], 16);
								t[0]=_mm256_srai_epi32(t[0], 16);
								t[1]=_mm256_srai_epi32(t[1], 16);
								t[2]=_mm256_srai_epi32(t[2], 16);
								t[0]=_mm256_mullo_epi32(t[0], _mm256_load_si256((__m256i*)L1weights+k*6+0));
								t[1]=_mm256_mullo_epi32(t[1], _mm256_load_si256((__m256i*)L1weights+k*6+1));
								t[2]=_mm256_mullo_epi32(t[2], _mm256_load_si256((__m256i*)L1weights+k*6+2));
								t[3]=_mm256_mullo_epi32(t[3], _mm256_load_si256((__m256i*)L1weights+k*6+3));
								t[4]=_mm256_mullo_epi32(t[4], _mm256_load_si256((__m256i*)L1weights+k*6+4));
								t[5]=_mm256_mullo_epi32(t[5], _mm256_load_si256((__m256i*)L1weights+k*6+5));
								mp[0]=_mm256_add_epi32(mp[0], t[0]);
								mp[1]=_mm256_add_epi32(mp[1], t[1]);
								mp[2]=_mm256_add_epi32(mp[2], t[2]);
								mp[3]=_mm256_add_epi32(mp[3], t[3]);
								mp[4]=_mm256_add_epi32(mp[4], t[4]);
								mp[5]=_mm256_add_epi32(mp[5], t[5]);
							}
							__m256i rcon=_mm256_set1_epi32(1<<L1_SH2>>1);
							mp[0]=_mm256_add_epi32(mp[0], rcon);//rounding to nearest
							mp[1]=_mm256_add_epi32(mp[1], rcon);
							mp[2]=_mm256_add_epi32(mp[2], rcon);
							mp[3]=_mm256_add_epi32(mp[3], rcon);
							mp[4]=_mm256_add_epi32(mp[4], rcon);
							mp[5]=_mm256_add_epi32(mp[5], rcon);

							mp[0]=_mm256_srai_epi32(mp[0], L1_SH2);
							mp[1]=_mm256_srai_epi32(mp[1], L1_SH2);
							mp[2]=_mm256_srai_epi32(mp[2], L1_SH2);
							mp[3]=_mm256_srai_epi32(mp[3], L1_SH2-16);
							mp[4]=_mm256_srai_epi32(mp[4], L1_SH2-16);
							mp[5]=_mm256_srai_epi32(mp[5], L1_SH2-16);
							//32 -> 16
							predY=_mm256_blend_epi16(mp[0], mp[3], 0xAA);
							predU=_mm256_blend_epi16(mp[1], mp[4], 0xAA);
							predV=_mm256_blend_epi16(mp[2], mp[5], 0xAA);


							//loosen pred range
							if(!cond_cg)
							{
								cache[0]=_mm256_load_si256((__m256i*)rows[1]+0+0+1*6);//NE
								cache[1]=_mm256_load_si256((__m256i*)rows[1]+0+1+1*6);
								cache[2]=_mm256_load_si256((__m256i*)rows[1]+0+2+1*6);
								ymin=_mm256_min_epi16(ymin, cache[0]);
								ymax=_mm256_max_epi16(ymax, cache[0]);
								umin=_mm256_min_epi16(umin, cache[1]);
								umax=_mm256_max_epi16(umax, cache[1]);
								vmin=_mm256_min_epi16(vmin, cache[2]);
								vmax=_mm256_max_epi16(vmax, cache[2]);
								cache[0]=_mm256_load_si256((__m256i*)rows[1]+0+0+3*6);//NEEE	crisp DIV2K -0.18%, noisy GDCC +0.13%
								cache[1]=_mm256_load_si256((__m256i*)rows[1]+0+1+3*6);
								cache[2]=_mm256_load_si256((__m256i*)rows[1]+0+2+3*6);
								ymin=_mm256_min_epi16(ymin, cache[0]);
								ymax=_mm256_max_epi16(ymax, cache[0]);
								umin=_mm256_min_epi16(umin, cache[1]);
								umax=_mm256_max_epi16(umax, cache[1]);
								vmin=_mm256_min_epi16(vmin, cache[2]);
								vmax=_mm256_max_epi16(vmax, cache[2]);
							}
						}
						else if(effort==3)
						{
							__m256i cache[3];
							/*
							effort 3
							0	N
							1	W
							2	NNN
							3	WWW
							4	NEEE
							5	3*(N-NN)+NNN
							6	3*(W-WW)+WWW
							7	W+NE-N
							8	N+W-NW
							9	N+NE-NNE
							10	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-(NW+N))>>2
							*/

							//N
							L1preds[0*3+0]=N[0];
							L1preds[0*3+1]=N[1];
							L1preds[0*3+2]=N[2];

							//W
							L1preds[1*3+0]=W[0];
							L1preds[1*3+1]=W[1];
							L1preds[1*3+2]=W[2];

							//NNN
							L1preds[2*3+0]=_mm256_load_si256((__m256i*)rows[3]+0+0+0*6);
							L1preds[2*3+1]=_mm256_load_si256((__m256i*)rows[3]+0+1+0*6);
							L1preds[2*3+2]=_mm256_load_si256((__m256i*)rows[3]+0+2+0*6);

							//WWW
							L1preds[3*3+0]=_mm256_load_si256((__m256i*)rows[0]+0+0-3*6);
							L1preds[3*3+1]=_mm256_load_si256((__m256i*)rows[0]+0+1-3*6);
							L1preds[3*3+2]=_mm256_load_si256((__m256i*)rows[0]+0+2-3*6);

							//NEEE
							L1preds[4*3+0]=_mm256_load_si256((__m256i*)rows[1]+0+0+3*6);
							L1preds[4*3+1]=_mm256_load_si256((__m256i*)rows[1]+0+1+3*6);
							L1preds[4*3+2]=_mm256_load_si256((__m256i*)rows[1]+0+2+3*6);

							//3*(N-NN)+NNN
							cache[0]=_mm256_sub_epi16(N[0], _mm256_load_si256((__m256i*)rows[2]+0+0+0*6));//N-NN
							cache[1]=_mm256_sub_epi16(N[1], _mm256_load_si256((__m256i*)rows[2]+0+1+0*6));
							cache[2]=_mm256_sub_epi16(N[2], _mm256_load_si256((__m256i*)rows[2]+0+2+0*6));
							cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
							cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
							L1preds[5*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[3]+0+0+0*6));//+NNN
							L1preds[5*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[3]+0+1+0*6));
							L1preds[5*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[3]+0+2+0*6));

							//3*(W-WW)+WWW
							cache[0]=_mm256_sub_epi16(W[0], _mm256_load_si256((__m256i*)rows[0]+0+0-2*6));//W-WW
							cache[1]=_mm256_sub_epi16(W[1], _mm256_load_si256((__m256i*)rows[0]+0+1-2*6));
							cache[2]=_mm256_sub_epi16(W[2], _mm256_load_si256((__m256i*)rows[0]+0+2-2*6));
							cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
							cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
							L1preds[6*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[0]+0+0-3*6));//+WWW
							L1preds[6*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[0]+0+1-3*6));
							L1preds[6*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[0]+0+2-3*6));

							//W+NE-N
							cache[0]=_mm256_sub_epi16(W[0], N[0]);
							cache[1]=_mm256_sub_epi16(W[1], N[1]);
							cache[2]=_mm256_sub_epi16(W[2], N[2]);
							L1preds[7*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+0+1*6));//+NE
							L1preds[7*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+1+1*6));
							L1preds[7*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+2+1*6));

							//N+W-NW
							L1preds[8*3+0]=predY;
							L1preds[8*3+1]=predU;
							L1preds[8*3+2]=predV;

							//N+NE-NNE
							cache[0]=_mm256_add_epi16(N[0], _mm256_load_si256((__m256i*)rows[1]+0+0+1*6));//N+NE
							cache[1]=_mm256_add_epi16(N[1], _mm256_load_si256((__m256i*)rows[1]+0+1+1*6));
							cache[2]=_mm256_add_epi16(N[2], _mm256_load_si256((__m256i*)rows[1]+0+2+1*6));
							L1preds[9*3+0]=_mm256_sub_epi16(cache[0], _mm256_load_si256((__m256i*)rows[2]+0+0+1*6));//NNE
							L1preds[9*3+1]=_mm256_sub_epi16(cache[1], _mm256_load_si256((__m256i*)rows[2]+0+1+1*6));
							L1preds[9*3+2]=_mm256_sub_epi16(cache[2], _mm256_load_si256((__m256i*)rows[2]+0+2+1*6));
					
							//(WWWW+WWW+NNN+NNEE+NEEE+NEEEE-(N+W))>>2
		#if 1
							cache[0]=_mm256_load_si256((__m256i*)rows[0]+0+0-4*6);//WWWW
							cache[1]=_mm256_load_si256((__m256i*)rows[0]+0+1-4*6);
							cache[2]=_mm256_load_si256((__m256i*)rows[0]+0+2-4*6);
							cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[0]+0+0-3*6));//+WWW
							cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[0]+0+1-3*6));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[0]+0+2-3*6));
							cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[3]+0+0+0*6));//+NNN
							cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[3]+0+1+0*6));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[3]+0+2+0*6));
							cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[2]+0+0+2*6));//+NNEE
							cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[2]+0+1+2*6));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[2]+0+2+2*6));
							cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+0+3*6));//+NEEE
							cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+1+3*6));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+2+3*6));
							cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+0+4*6));//+NEEEE
							cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+1+4*6));
							cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+2+4*6));
							cache[0]=_mm256_sub_epi16(cache[0], _mm256_add_epi16(N[0], W[0]));//-(N+W)
							cache[1]=_mm256_sub_epi16(cache[1], _mm256_add_epi16(N[1], W[1]));
							cache[2]=_mm256_sub_epi16(cache[2], _mm256_add_epi16(N[2], W[2]));
							L1preds[10*3+0]=_mm256_srai_epi16(cache[0], 2);
							L1preds[10*3+1]=_mm256_srai_epi16(cache[1], 2);
							L1preds[10*3+2]=_mm256_srai_epi16(cache[2], 2);
		#endif


							//mix
							__m256i mp[6], t[6];
							mp[0]=_mm256_setzero_si256();
							mp[1]=_mm256_setzero_si256();
							mp[2]=_mm256_setzero_si256();
							mp[3]=_mm256_setzero_si256();
							mp[4]=_mm256_setzero_si256();
							mp[5]=_mm256_setzero_si256();

							//mp[0]=_mm256_load_si256((__m256i*)L1coeffs+L1_NPREDS3*6+0);//bias
							//mp[1]=_mm256_load_si256((__m256i*)L1coeffs+L1_NPREDS3*6+1);
							//mp[2]=_mm256_load_si256((__m256i*)L1coeffs+L1_NPREDS3*6+2);
							//mp[3]=_mm256_load_si256((__m256i*)L1coeffs+L1_NPREDS3*6+3);
							//mp[4]=_mm256_load_si256((__m256i*)L1coeffs+L1_NPREDS3*6+4);
							//mp[5]=_mm256_load_si256((__m256i*)L1coeffs+L1_NPREDS3*6+5);
		#ifdef PRINT_L1_BOUNDS
							for(int k2=0;k2<6*8;++k2)
							{
								int bias=L1coeffs[L1_NPREDS3*6*8+k2];
								if(bmin>bias)bmin=bias;
								if(bmax<bias)bmax=bias;
							}
		#endif
		//#if defined __GNUC__ && !defined PROFILER
		//#pragma GCC unroll 11
		//#endif
							for(int k=0;k<L1_NPREDS3;++k)
							{
								//16 -> 32		3 lo 3 hi registers
								t[0]=_mm256_slli_epi32(L1preds[k*3+0], 16);
								t[1]=_mm256_slli_epi32(L1preds[k*3+1], 16);
								t[2]=_mm256_slli_epi32(L1preds[k*3+2], 16);
								t[3]=_mm256_srai_epi32(L1preds[k*3+0], 16);
								t[4]=_mm256_srai_epi32(L1preds[k*3+1], 16);
								t[5]=_mm256_srai_epi32(L1preds[k*3+2], 16);
								t[0]=_mm256_srai_epi32(t[0], 16);
								t[1]=_mm256_srai_epi32(t[1], 16);
								t[2]=_mm256_srai_epi32(t[2], 16);
								t[0]=_mm256_mullo_epi32(t[0], _mm256_load_si256((__m256i*)L1weights+k*6+0));
								t[1]=_mm256_mullo_epi32(t[1], _mm256_load_si256((__m256i*)L1weights+k*6+1));
								t[2]=_mm256_mullo_epi32(t[2], _mm256_load_si256((__m256i*)L1weights+k*6+2));
								t[3]=_mm256_mullo_epi32(t[3], _mm256_load_si256((__m256i*)L1weights+k*6+3));
								t[4]=_mm256_mullo_epi32(t[4], _mm256_load_si256((__m256i*)L1weights+k*6+4));
								t[5]=_mm256_mullo_epi32(t[5], _mm256_load_si256((__m256i*)L1weights+k*6+5));
								mp[0]=_mm256_add_epi32(mp[0], t[0]);
								mp[1]=_mm256_add_epi32(mp[1], t[1]);
								mp[2]=_mm256_add_epi32(mp[2], t[2]);
								mp[3]=_mm256_add_epi32(mp[3], t[3]);
								mp[4]=_mm256_add_epi32(mp[4], t[4]);
								mp[5]=_mm256_add_epi32(mp[5], t[5]);
							}
							__m256i rcon=_mm256_set1_epi32(1<<L1_SH3>>1);
							mp[0]=_mm256_add_epi32(mp[0], rcon);//rounding to nearest
							mp[1]=_mm256_add_epi32(mp[1], rcon);
							mp[2]=_mm256_add_epi32(mp[2], rcon);
							mp[3]=_mm256_add_epi32(mp[3], rcon);
							mp[4]=_mm256_add_epi32(mp[4], rcon);
							mp[5]=_mm256_add_epi32(mp[5], rcon);

							mp[0]=_mm256_srai_epi32(mp[0], L1_SH3);
							mp[1]=_mm256_srai_epi32(mp[1], L1_SH3);
							mp[2]=_mm256_srai_epi32(mp[2], L1_SH3);
							mp[3]=_mm256_srai_epi32(mp[3], L1_SH3-16);
							mp[4]=_mm256_srai_epi32(mp[4], L1_SH3-16);
							mp[5]=_mm256_srai_epi32(mp[5], L1_SH3-16);
							//32 -> 16
							predY=_mm256_blend_epi16(mp[0], mp[3], 0xAA);
							predU=_mm256_blend_epi16(mp[1], mp[4], 0xAA);
							predV=_mm256_blend_epi16(mp[2], mp[5], 0xAA);


							//loosen pred range
							if(!cond_cg)
							{
								cache[0]=_mm256_load_si256((__m256i*)rows[1]+0+0+1*6);//NE
								cache[1]=_mm256_load_si256((__m256i*)rows[1]+0+1+1*6);
								cache[2]=_mm256_load_si256((__m256i*)rows[1]+0+2+1*6);
								ymin=_mm256_min_epi16(ymin, cache[0]);
								ymax=_mm256_max_epi16(ymax, cache[0]);
								umin=_mm256_min_epi16(umin, cache[1]);
								umax=_mm256_max_epi16(umax, cache[1]);
								vmin=_mm256_min_epi16(vmin, cache[2]);
								vmax=_mm256_max_epi16(vmax, cache[2]);
								cache[0]=_mm256_load_si256((__m256i*)rows[1]+0+0+3*6);//NEEE	crisp DIV2K -0.18%, noisy GDCC +0.13%
								cache[1]=_mm256_load_si256((__m256i*)rows[1]+0+1+3*6);
								cache[2]=_mm256_load_si256((__m256i*)rows[1]+0+2+3*6);
								ymin=_mm256_min_epi16(ymin, cache[0]);
								ymax=_mm256_max_epi16(ymax, cache[0]);
								umin=_mm256_min_epi16(umin, cache[1]);
								umax=_mm256_max_epi16(umax, cache[1]);
								vmin=_mm256_min_epi16(vmin, cache[2]);
								vmax=_mm256_max_epi16(vmax, cache[2]);
							}
						}
						predYUV0[0]=predY;
						predYUV0[1]=predU;
						predYUV0[2]=predV;
						if(cond_cg)
						{
							predY=mcg[0];
							predU=mcg[1];
							predV=mcg[2];
						}
#endif

						predY=_mm256_max_epi16(predY, ymin);
						predU=_mm256_max_epi16(predU, umin);
						predV=_mm256_max_epi16(predV, vmin);
						predY=_mm256_min_epi16(predY, ymax);
						predU=_mm256_min_epi16(predU, umax);
						predV=_mm256_min_epi16(predV, vmax);

						//predY=_mm256_setzero_si256();//
						//predU=_mm256_setzero_si256();//
						//predV=_mm256_setzero_si256();//
					}
					
					//if(kb==1&&ky==35&&kx==83)//
					//if(kb==0&&ky==0&&kx==0)//
					//if(kb==0&&ky==0&&kx==2)//
					//	printf("");

					if(dist>1)
					{
						__m256i val[3], tmp;
						if(fwd)
						{
							myuv[0]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+0*NLANES)), half8));//load rgb
							myuv[1]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+1*NLANES)), half8));
							myuv[2]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+2*NLANES)), half8));

							//val=(curr-(int)pred)/dist
							val[0]=_mm256_sub_epi16(myuv[0], predY);
							val[1]=_mm256_sub_epi16(myuv[1], predU);
							val[2]=_mm256_sub_epi16(myuv[2], predV);
							__m256i cond[3];
							cond[0]=_mm256_srai_epi16(val[0], 15);
							cond[1]=_mm256_srai_epi16(val[1], 15);
							cond[2]=_mm256_srai_epi16(val[2], 15);
							val[0]=_mm256_mulhi_epi16(val[0], dist_rcp);
							val[1]=_mm256_mulhi_epi16(val[1], dist_rcp);
							val[2]=_mm256_mulhi_epi16(val[2], dist_rcp);
							val[0]=_mm256_sub_epi16(val[0], cond[0]);//(x*inv>>16)-(x>>31)
							val[1]=_mm256_sub_epi16(val[1], cond[1]);
							val[2]=_mm256_sub_epi16(val[2], cond[2]);

							//curr=dist*val+pred
							myuv[0]=_mm256_mullo_epi16(val[0], mdist);
							myuv[1]=_mm256_mullo_epi16(val[1], mdist);
							myuv[2]=_mm256_mullo_epi16(val[2], mdist);

							myuv[0]=_mm256_add_epi16(myuv[0], predY);
							myuv[1]=_mm256_add_epi16(myuv[1], predU);
							myuv[2]=_mm256_add_epi16(myuv[2], predV);

							tmp=_mm256_set1_epi16(-128);
							myuv[0]=_mm256_max_epi16(myuv[0], tmp);
							myuv[1]=_mm256_max_epi16(myuv[1], tmp);
							myuv[2]=_mm256_max_epi16(myuv[2], tmp);
							tmp=_mm256_set1_epi16(127);
							myuv[0]=_mm256_min_epi16(myuv[0], tmp);
							myuv[1]=_mm256_min_epi16(myuv[1], tmp);
							myuv[2]=_mm256_min_epi16(myuv[2], tmp);
					
							//RCT
							val[0]=_mm256_sub_epi16(val[0], val[1]);
							val[2]=_mm256_sub_epi16(val[2], val[1]);
							tmp=_mm256_add_epi16(val[0], val[2]);
							tmp=_mm256_srai_epi16(tmp, 2);
							val[1]=_mm256_add_epi16(val[1], tmp);

							ecurr[0]=val[0];
							ecurr[1]=val[1];
							ecurr[2]=val[2];
#ifdef _DEBUG
							for(int k=0;k<48;++k)
							{
								int v=((int16_t*)val)[k];
								if((unsigned)(v+128)>=256)
									CRASH("");
							}
#endif
							tmp=_mm256_set1_epi16(128);
							val[0]=_mm256_add_epi16(val[0], tmp);
							val[1]=_mm256_add_epi16(val[1], tmp);
							val[2]=_mm256_add_epi16(val[2], tmp);
							tmp=_mm256_set1_epi16(255);
							val[0]=_mm256_and_si256(val[0], tmp);
							val[1]=_mm256_and_si256(val[1], tmp);
							val[2]=_mm256_and_si256(val[2], tmp);
					
							ctxU=_mm256_add_epi16(ctxU, mctxuoffset);
							ctxV=_mm256_add_epi16(ctxV, mctxvoffset);
							ctxY=_mm256_slli_epi16(ctxY, 8);
							ctxU=_mm256_slli_epi16(ctxU, 8);
							ctxV=_mm256_slli_epi16(ctxV, 8);
							ctxY=_mm256_or_si256(ctxY, val[0]);
							ctxU=_mm256_or_si256(ctxU, val[1]);
							ctxV=_mm256_or_si256(ctxV, val[2]);
							_mm256_store_si256((__m256i*)syms+0, ctxY);
							_mm256_store_si256((__m256i*)syms+1, ctxU);
							_mm256_store_si256((__m256i*)syms+2, ctxV);
							_mm256_store_si256((__m256i*)ctxptr+0, ctxY);
							_mm256_store_si256((__m256i*)ctxptr+1, ctxU);
							_mm256_store_si256((__m256i*)ctxptr+2, ctxV);

#if 0
							for(int k=0;k<48;++k)
							{
								if((unsigned)syms[k]>=(unsigned)(18*3<<8))
									CRASH("");
							}
#endif
							int *pa, *pb, *pc, va, vb, vc;
							pa=hists+syms[0*16+0x0]; pb=hists+syms[1*16+0x0]; pc=hists+syms[2*16+0x0]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0x1]; pb=hists+syms[1*16+0x1]; pc=hists+syms[2*16+0x1]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0x2]; pb=hists+syms[1*16+0x2]; pc=hists+syms[2*16+0x2]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0x3]; pb=hists+syms[1*16+0x3]; pc=hists+syms[2*16+0x3]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0x4]; pb=hists+syms[1*16+0x4]; pc=hists+syms[2*16+0x4]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0x5]; pb=hists+syms[1*16+0x5]; pc=hists+syms[2*16+0x5]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0x6]; pb=hists+syms[1*16+0x6]; pc=hists+syms[2*16+0x6]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0x7]; pb=hists+syms[1*16+0x7]; pc=hists+syms[2*16+0x7]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0x8]; pb=hists+syms[1*16+0x8]; pc=hists+syms[2*16+0x8]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0x9]; pb=hists+syms[1*16+0x9]; pc=hists+syms[2*16+0x9]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0xA]; pb=hists+syms[1*16+0xA]; pc=hists+syms[2*16+0xA]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0xB]; pb=hists+syms[1*16+0xB]; pc=hists+syms[2*16+0xB]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0xC]; pb=hists+syms[1*16+0xC]; pc=hists+syms[2*16+0xC]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0xD]; pb=hists+syms[1*16+0xD]; pc=hists+syms[2*16+0xD]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0xE]; pb=hists+syms[1*16+0xE]; pc=hists+syms[2*16+0xE]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							pa=hists+syms[0*16+0xF]; pb=hists+syms[1*16+0xF]; pc=hists+syms[2*16+0xF]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
							ctxptr+=3*NLANES;
						}
						else
						{
							__m128i msyms24[3];

							dec_yuv(mstate, 0, &ctxY, (int*)CDF2syms, &streamptr, streamend, myuv+0);
							dec_yuv(mstate, 1, &ctxU, (int*)CDF2syms, &streamptr, streamend, myuv+1);
							dec_yuv(mstate, 2, &ctxV, (int*)CDF2syms, &streamptr, streamend, myuv+2);
					
							tmp=_mm256_set1_epi16(128);
							val[0]=_mm256_sub_epi16(myuv[0], tmp);
							val[1]=_mm256_sub_epi16(myuv[1], tmp);
							val[2]=_mm256_sub_epi16(myuv[2], tmp);
							ecurr[0]=val[0];
							ecurr[1]=val[1];
							ecurr[2]=val[2];
					
							//RCT
							tmp=_mm256_add_epi16(val[0], val[2]);
							tmp=_mm256_srai_epi16(tmp, 2);
							val[1]=_mm256_sub_epi16(val[1], tmp);
							val[0]=_mm256_add_epi16(val[0], val[1]);
							val[2]=_mm256_add_epi16(val[2], val[1]);

							//val=dist*curr+pred
							val[0]=_mm256_mullo_epi16(val[0], mdist);
							val[1]=_mm256_mullo_epi16(val[1], mdist);
							val[2]=_mm256_mullo_epi16(val[2], mdist);
							val[0]=_mm256_add_epi16(val[0], predY);
							val[1]=_mm256_add_epi16(val[1], predU);
							val[2]=_mm256_add_epi16(val[2], predV);

							tmp=_mm256_set1_epi16(-128);
							val[0]=_mm256_max_epi16(val[0], tmp);
							val[1]=_mm256_max_epi16(val[1], tmp);
							val[2]=_mm256_max_epi16(val[2], tmp);
							tmp=_mm256_set1_epi16(127);
							val[0]=_mm256_min_epi16(val[0], tmp);
							val[1]=_mm256_min_epi16(val[1], tmp);
							val[2]=_mm256_min_epi16(val[2], tmp);
							myuv[0]=val[0];
							myuv[1]=val[1];
							myuv[2]=val[2];

							msyms24[0]=_mm_packs_epi16(_mm256_extracti128_si256(val[0], 0), _mm256_extracti128_si256(val[0], 1));
							msyms24[1]=_mm_packs_epi16(_mm256_extracti128_si256(val[1], 0), _mm256_extracti128_si256(val[1], 1));
							msyms24[2]=_mm_packs_epi16(_mm256_extracti128_si256(val[2], 0), _mm256_extracti128_si256(val[2], 1));
							msyms24[0]=_mm_xor_si128(msyms24[0], _mm_set1_epi8(-128));
							msyms24[1]=_mm_xor_si128(msyms24[1], _mm_set1_epi8(-128));
							msyms24[2]=_mm_xor_si128(msyms24[2], _mm_set1_epi8(-128));
							_mm_store_si128((__m128i*)(imptr+0*NLANES), msyms24[0]);//store rgb
							_mm_store_si128((__m128i*)(imptr+1*NLANES), msyms24[1]);
							_mm_store_si128((__m128i*)(imptr+2*NLANES), msyms24[2]);
#ifdef ENABLE_GUIDE
							for(int kc=0;kc<3;++kc)
							{
								for(int k=0;k<NLANES;++k)
								{
									double diff=imptr[kc*NLANES+k]-g_image[3*(iw*(y1+YLANES*ky+k/XLANES)+x1+XLANES*kx+k%XLANES)+combination[II_PERM_Y+kc]];
									g_sqe[kc]+=diff*diff;
								}
							}
#endif
						}
						_mm256_store_si256((__m256i*)rows[0]+0+0+0*6, myuv[0]);
						_mm256_store_si256((__m256i*)rows[0]+0+1+0*6, myuv[1]);
						_mm256_store_si256((__m256i*)rows[0]+0+2+0*6, myuv[2]);

						W[0]=myuv[0];
						W[1]=myuv[1];
						W[2]=myuv[2];
						ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(ecurr[0], 1), _mm256_srai_epi16(ecurr[0], 15));
						ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(ecurr[1], 1), _mm256_srai_epi16(ecurr[1], 15));
						ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(ecurr[2], 1), _mm256_srai_epi16(ecurr[2], 15));
					}
					else if(fwd)
					{
						//if(ky==1&&kx==1296)//
						//if(ky==0&&kx==48)//
						//if(ky==1&&kx==1296)//
						//	printf("");

						__m256i ctxblendmask=_mm256_set1_epi16(255);
						myuv[0]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+yidx)), half8));//load yuv
						myuv[1]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+uidx)), half8));
						myuv[2]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+vidx)), half8));

						//encode Y
						W[0]=myuv[0];
						_mm256_store_si256((__m256i*)rows[0]+0+0+0*6, myuv[0]);//store Y neighbors
						msyms=_mm256_sub_epi16(myuv[0], predY);//sub pred
#ifdef SAVE_RESIDUALS
						{
							ptrdiff_t idx=imptr-interleaved-isize+(ptrdiff_t)0*NLANES;
							int16_t syms2[16];
							_mm256_storeu_si256((__m256i*)syms2, msyms);
							for(int k=0;k<16;++k)
								residuals[idx+k]=(unsigned char)(syms2[k]+128);
						}
#endif
						ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));//ecurr = pack_sign(yuv-pred)
						msyms=_mm256_sub_epi16(msyms, amin);
						ctxY=_mm256_slli_epi16(ctxY, 8);
						ctxY=_mm256_blendv_epi8(ctxY, msyms, ctxblendmask);
						_mm256_store_si256((__m256i*)syms, ctxY);
						_mm256_store_si256((__m256i*)ctxptr+0, ctxY);//store Y  ctx|residuals

						//encode U
						moffset=_mm256_and_si256(myuv[0], uhelpmask);
						predU=_mm256_add_epi16(predU, moffset);
						predU=_mm256_max_epi16(predU, amin);
						predU=_mm256_min_epi16(predU, amax);

						msyms=_mm256_sub_epi16(myuv[1], predU);
#ifdef SAVE_RESIDUALS
						{
							ptrdiff_t idx=imptr-interleaved-isize+(ptrdiff_t)1*NLANES;
							int16_t syms2[16];
							_mm256_storeu_si256((__m256i*)syms2, msyms);
							for(int k=0;k<16;++k)
								residuals[idx+k]=(unsigned char)(syms2[k]+128);
						}
#endif
						ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
						msyms=_mm256_sub_epi16(msyms, amin);
						ctxU=_mm256_add_epi16(ctxU, mctxuoffset);
						ctxU=_mm256_slli_epi16(ctxU, 8);
						ctxU=_mm256_blendv_epi8(ctxU, msyms, ctxblendmask);
						_mm256_store_si256((__m256i*)syms+1, ctxU);
						_mm256_store_si256((__m256i*)ctxptr+1, ctxU);//store U  ctx|residuals
						W[1]=_mm256_sub_epi16(myuv[1], moffset);
						_mm256_store_si256((__m256i*)rows[0]+0+1+0*6, W[1]);//store U neighbors

						//encode V
						moffset=_mm256_mullo_epi16(vc0, myuv[0]);
						moffset=_mm256_add_epi16(moffset, _mm256_mullo_epi16(vc1, myuv[1]));
						moffset=_mm256_srai_epi16(moffset, 2);
						predV=_mm256_add_epi16(predV, moffset);
						predV=_mm256_max_epi16(predV, amin);
						predV=_mm256_min_epi16(predV, amax);

						msyms=_mm256_sub_epi16(myuv[2], predV);
#ifdef SAVE_RESIDUALS
						{
							ptrdiff_t idx=imptr-interleaved-isize+(ptrdiff_t)2*NLANES;
							int16_t syms2[16];
							_mm256_storeu_si256((__m256i*)syms2, msyms);
							for(int k=0;k<16;++k)
								residuals[idx+k]=(unsigned char)(syms2[k]+128);
						}
#endif
						ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
						msyms=_mm256_sub_epi16(msyms, amin);
						ctxV=_mm256_add_epi16(ctxV, mctxvoffset);
						ctxV=_mm256_slli_epi16(ctxV, 8);
						ctxV=_mm256_blendv_epi8(ctxV, msyms, ctxblendmask);
						_mm256_store_si256((__m256i*)syms+2, ctxV);
						_mm256_store_si256((__m256i*)ctxptr+2, ctxV);//store V  ctx|residuals		ctxptr+NLANES*(C*2+R)
						W[2]=_mm256_sub_epi16(myuv[2], moffset);
						_mm256_store_si256((__m256i*)rows[0]+0+2+0*6, W[2]);//store V neighbors
#if 1
						int *pa, *pb, *pc, va, vb, vc;
						pa=hists+syms[0*16+0x0]; pb=hists+syms[1*16+0x0]; pc=hists+syms[2*16+0x0]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0x1]; pb=hists+syms[1*16+0x1]; pc=hists+syms[2*16+0x1]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0x2]; pb=hists+syms[1*16+0x2]; pc=hists+syms[2*16+0x2]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0x3]; pb=hists+syms[1*16+0x3]; pc=hists+syms[2*16+0x3]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0x4]; pb=hists+syms[1*16+0x4]; pc=hists+syms[2*16+0x4]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0x5]; pb=hists+syms[1*16+0x5]; pc=hists+syms[2*16+0x5]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0x6]; pb=hists+syms[1*16+0x6]; pc=hists+syms[2*16+0x6]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0x7]; pb=hists+syms[1*16+0x7]; pc=hists+syms[2*16+0x7]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0x8]; pb=hists+syms[1*16+0x8]; pc=hists+syms[2*16+0x8]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0x9]; pb=hists+syms[1*16+0x9]; pc=hists+syms[2*16+0x9]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0xA]; pb=hists+syms[1*16+0xA]; pc=hists+syms[2*16+0xA]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0xB]; pb=hists+syms[1*16+0xB]; pc=hists+syms[2*16+0xB]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0xC]; pb=hists+syms[1*16+0xC]; pc=hists+syms[2*16+0xC]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0xD]; pb=hists+syms[1*16+0xD]; pc=hists+syms[2*16+0xD]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0xE]; pb=hists+syms[1*16+0xE]; pc=hists+syms[2*16+0xE]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
						pa=hists+syms[0*16+0xF]; pb=hists+syms[1*16+0xF]; pc=hists+syms[2*16+0xF]; va=*pa+1; vb=*pb+1; vc=*pc+1; *pa=va; *pb=vb; *pc=vc;
#endif
#if 0
						++hists[syms[0*16+0x0]];
						++hists[syms[1*16+0x0]];
						++hists[syms[2*16+0x0]];
						++hists[syms[0*16+0x1]];
						++hists[syms[1*16+0x1]];
						++hists[syms[2*16+0x1]];
						++hists[syms[0*16+0x2]];
						++hists[syms[1*16+0x2]];
						++hists[syms[2*16+0x2]];
						++hists[syms[0*16+0x3]];
						++hists[syms[1*16+0x3]];
						++hists[syms[2*16+0x3]];
						++hists[syms[0*16+0x4]];
						++hists[syms[1*16+0x4]];
						++hists[syms[2*16+0x4]];
						++hists[syms[0*16+0x5]];
						++hists[syms[1*16+0x5]];
						++hists[syms[2*16+0x5]];
						++hists[syms[0*16+0x6]];
						++hists[syms[1*16+0x6]];
						++hists[syms[2*16+0x6]];
						++hists[syms[0*16+0x7]];
						++hists[syms[1*16+0x7]];
						++hists[syms[2*16+0x7]];
						++hists[syms[0*16+0x8]];
						++hists[syms[1*16+0x8]];
						++hists[syms[2*16+0x8]];
						++hists[syms[0*16+0x9]];
						++hists[syms[1*16+0x9]];
						++hists[syms[2*16+0x9]];
						++hists[syms[0*16+0xA]];
						++hists[syms[1*16+0xA]];
						++hists[syms[2*16+0xA]];
						++hists[syms[0*16+0xB]];
						++hists[syms[1*16+0xB]];
						++hists[syms[2*16+0xB]];
						++hists[syms[0*16+0xC]];
						++hists[syms[1*16+0xC]];
						++hists[syms[2*16+0xC]];
						++hists[syms[0*16+0xD]];
						++hists[syms[1*16+0xD]];
						++hists[syms[2*16+0xD]];
						++hists[syms[0*16+0xE]];
						++hists[syms[1*16+0xE]];
						++hists[syms[2*16+0xE]];
						++hists[syms[0*16+0xF]];
						++hists[syms[1*16+0xF]];
						++hists[syms[2*16+0xF]];
#endif
						ctxptr+=3*NLANES;
					}
					else
					{
						//if(ky==1&&kx==1296)//
						//if(ky==0&&kx==48)//
						//if(kb==0&&ky==0&&kx==2)//
						//	printf("");

						//decode main
						__m128i msyms8;

						//decode Y
						dec_yuv(mstate, 0, &ctxY, (int*)CDF2syms, &streamptr, streamend, myuv+0);//residuals from [0 ~ 255]
						dec_yuv(mstate, 1, &ctxU, (int*)CDF2syms, &streamptr, streamend, myuv+1);
						dec_yuv(mstate, 2, &ctxV, (int*)CDF2syms, &streamptr, streamend, myuv+2);

						//yuv = (char)(sym+pred-128)	= (unsigned char)(sym+pred)-128
						myuv[0]=_mm256_add_epi16(myuv[0], predY);
						myuv[0]=_mm256_and_si256(myuv[0], bytemask);
						myuv[0]=_mm256_add_epi16(myuv[0], amin);
						msyms=_mm256_sub_epi16(myuv[0], predY);//sub pred
						ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
						msyms8=_mm_packs_epi16(_mm256_extracti128_si256(myuv[0], 0), _mm256_extracti128_si256(myuv[0], 1));
						msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
						_mm_store_si128((__m128i*)(imptr+yidx), msyms8);//store Y bytes
						W[0]=myuv[0];
						_mm256_store_si256((__m256i*)rows[0]+0+0+0*6, myuv[0]);//store Y neighbors


						//decode U
						moffset=_mm256_and_si256(myuv[0], uhelpmask);
						predU=_mm256_add_epi16(predU, moffset);
						predU=_mm256_max_epi16(predU, amin);
						predU=_mm256_min_epi16(predU, amax);

						myuv[1]=_mm256_add_epi16(myuv[1], predU);
						myuv[1]=_mm256_and_si256(myuv[1], bytemask);
						myuv[1]=_mm256_add_epi16(myuv[1], amin);
						msyms=_mm256_sub_epi16(myuv[1], predU);//sub pred
						ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
						msyms8=_mm_packs_epi16(_mm256_extracti128_si256(myuv[1], 0), _mm256_extracti128_si256(myuv[1], 1));
						msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
						_mm_store_si128((__m128i*)(imptr+uidx), msyms8);//store U bytes
						W[1]=_mm256_sub_epi16(myuv[1], moffset);//subtract Uoffset from U
						_mm256_store_si256((__m256i*)rows[0]+0+1+0*6, W[1]);//store U neighbors
				

						//decode V
						moffset=_mm256_mullo_epi16(vc0, myuv[0]);
						moffset=_mm256_add_epi16(moffset, _mm256_mullo_epi16(vc1, myuv[1]));
						moffset=_mm256_srai_epi16(moffset, 2);
						predV=_mm256_add_epi16(predV, moffset);
						predV=_mm256_max_epi16(predV, amin);
						predV=_mm256_min_epi16(predV, amax);
				
						myuv[2]=_mm256_add_epi16(myuv[2], predV);
						myuv[2]=_mm256_and_si256(myuv[2], bytemask);
						myuv[2]=_mm256_add_epi16(myuv[2], amin);
						msyms=_mm256_sub_epi16(myuv[2], predV);
						ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
						msyms8=_mm_packs_epi16(_mm256_extracti128_si256(myuv[2], 0), _mm256_extracti128_si256(myuv[2], 1));
						msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
						_mm_store_si128((__m128i*)(imptr+vidx), msyms8);//store V bytes
						W[2]=_mm256_sub_epi16(myuv[2], moffset);//subtract Voffset from V
						_mm256_store_si256((__m256i*)rows[0]+0+2+0*6, W[2]);//store V neighbors
#ifdef ENABLE_GUIDE
						{
							int bad=0;
							for(int kc=0;kc<3;++kc)
							{
								for(int k=0;k<NLANES;++k)
								{
									int v0=g_image[3*(iw*(y1+k/XLANES*iblockh+ky)+x1+k%XLANES*iblockw+kx)+kc];
									int v1=imptr[NLANES*kc+k];
									bad|=v1!=v0;
								}
							}
							if(bad)
							{
								int invperm[3]={0};
								invperm[combination[II_PERM_Y]]=0;
								invperm[combination[II_PERM_U]]=1;
								invperm[combination[II_PERM_V]]=2;
								printf("Block %d RCT %d %s\n", kb, bestrct, rct_names[bestrct]);
								for(int kc=0;kc<3;++kc)
								{
									for(int k=0;k<NLANES;++k)
									{
										int kx2=x1+k%XLANES*iblockw+kx;
										int ky2=y1+k/XLANES*iblockh+ky;
										int v0=g_image[3*(iw*ky2+kx2)+kc];
										int v1=imptr[NLANES*kc+k];
										printf("CXYL \"%c\" %d %d %3d  \"%c\" %5d %5d  OG 0x%02X %4d  %c=  %c 0x%02X %4d  diff %4d\n"
											, "YUV"[invperm[kc]], kx, ky, k
											, "RGB"[kc], kx2, ky2
											, v0, v0-128
											, v1==v0?'=':'!', v1==v0?' ':'X'
											, v1, v1-128
											, v0-v1
										);
									}
									printf("\n");
								}
								CRASH("");
								return 1;
							}
						}
#endif
					}
#if 0
					if(effort==1)//update
					{
						__m256i mu[3];

						mu[0]=_mm256_sub_epi16(W[0], predYUV0[0]);
						mu[1]=_mm256_sub_epi16(W[1], predYUV0[1]);
						mu[2]=_mm256_sub_epi16(W[2], predYUV0[2]);
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 4
#endif
						for(int k=0;k<L1_NPREDS1;++k)//update
						{
							__m256i mc[6];
							mc[0]=_mm256_sign_epi16(L1preds[k*3+0], mu[0]);//L1
							mc[1]=_mm256_sign_epi16(L1preds[k*3+1], mu[1]);
							mc[2]=_mm256_sign_epi16(L1preds[k*3+2], mu[2]);
							//mc[0]=_mm256_mullo_epi16(L1preds[k*3+0], mu[0]);//L2
							//mc[1]=_mm256_mullo_epi16(L1preds[k*3+1], mu[1]);
							//mc[2]=_mm256_mullo_epi16(L1preds[k*3+2], mu[2]);

							//16 -> 32	3 lo 3 hi registers
							mc[3]=_mm256_srai_epi32(mc[0], 16);
							mc[4]=_mm256_srai_epi32(mc[1], 16);
							mc[5]=_mm256_srai_epi32(mc[2], 16);
							mc[0]=_mm256_slli_epi32(mc[0], 16);
							mc[1]=_mm256_slli_epi32(mc[1], 16);
							mc[2]=_mm256_slli_epi32(mc[2], 16);
							mc[0]=_mm256_srai_epi32(mc[0], 16);
							mc[1]=_mm256_srai_epi32(mc[1], 16);
							mc[2]=_mm256_srai_epi32(mc[2], 16);
							mc[0]=_mm256_add_epi32(mc[0], _mm256_load_si256((__m256i*)L1weights+k*6+0));//update coeffs
							mc[1]=_mm256_add_epi32(mc[1], _mm256_load_si256((__m256i*)L1weights+k*6+1));
							mc[2]=_mm256_add_epi32(mc[2], _mm256_load_si256((__m256i*)L1weights+k*6+2));
							mc[3]=_mm256_add_epi32(mc[3], _mm256_load_si256((__m256i*)L1weights+k*6+3));
							mc[4]=_mm256_add_epi32(mc[4], _mm256_load_si256((__m256i*)L1weights+k*6+4));
							mc[5]=_mm256_add_epi32(mc[5], _mm256_load_si256((__m256i*)L1weights+k*6+5));
							_mm256_store_si256((__m256i*)L1weights+k*6+0, mc[0]);
							_mm256_store_si256((__m256i*)L1weights+k*6+1, mc[1]);
							_mm256_store_si256((__m256i*)L1weights+k*6+2, mc[2]);
							_mm256_store_si256((__m256i*)L1weights+k*6+3, mc[3]);
							_mm256_store_si256((__m256i*)L1weights+k*6+4, mc[4]);
							_mm256_store_si256((__m256i*)L1weights+k*6+5, mc[5]);
						}
					}
					else if(effort==2)//update
					{
						__m256i mu[3];

						mu[0]=_mm256_sub_epi16(W[0], predYUV0[0]);
						mu[1]=_mm256_sub_epi16(W[1], predYUV0[1]);
						mu[2]=_mm256_sub_epi16(W[2], predYUV0[2]);
		#if defined __GNUC__ && !defined PROFILER
		#pragma GCC unroll 8
		#endif
						for(int k=0;k<L1_NPREDS2;++k)//update
						{
							__m256i mc[6];
							mc[0]=_mm256_sign_epi16(L1preds[k*3+0], mu[0]);
							mc[1]=_mm256_sign_epi16(L1preds[k*3+1], mu[1]);
							mc[2]=_mm256_sign_epi16(L1preds[k*3+2], mu[2]);
							//16 -> 32	3 lo 3 hi registers
							mc[3]=_mm256_srai_epi32(mc[0], 16);
							mc[4]=_mm256_srai_epi32(mc[1], 16);
							mc[5]=_mm256_srai_epi32(mc[2], 16);
							mc[0]=_mm256_slli_epi32(mc[0], 16);
							mc[1]=_mm256_slli_epi32(mc[1], 16);
							mc[2]=_mm256_slli_epi32(mc[2], 16);
							mc[0]=_mm256_srai_epi32(mc[0], 16);
							mc[1]=_mm256_srai_epi32(mc[1], 16);
							mc[2]=_mm256_srai_epi32(mc[2], 16);
							mc[0]=_mm256_add_epi32(mc[0], _mm256_load_si256((__m256i*)L1weights+k*6+0));//update coeffs
							mc[1]=_mm256_add_epi32(mc[1], _mm256_load_si256((__m256i*)L1weights+k*6+1));
							mc[2]=_mm256_add_epi32(mc[2], _mm256_load_si256((__m256i*)L1weights+k*6+2));
							mc[3]=_mm256_add_epi32(mc[3], _mm256_load_si256((__m256i*)L1weights+k*6+3));
							mc[4]=_mm256_add_epi32(mc[4], _mm256_load_si256((__m256i*)L1weights+k*6+4));
							mc[5]=_mm256_add_epi32(mc[5], _mm256_load_si256((__m256i*)L1weights+k*6+5));
							_mm256_store_si256((__m256i*)L1weights+k*6+0, mc[0]);
							_mm256_store_si256((__m256i*)L1weights+k*6+1, mc[1]);
							_mm256_store_si256((__m256i*)L1weights+k*6+2, mc[2]);
							_mm256_store_si256((__m256i*)L1weights+k*6+3, mc[3]);
							_mm256_store_si256((__m256i*)L1weights+k*6+4, mc[4]);
							_mm256_store_si256((__m256i*)L1weights+k*6+5, mc[5]);
						}
					}
					else if(effort==3)//update
					{
						__m256i mu[3];

						mu[0]=_mm256_sub_epi16(W[0], predYUV0[0]);
						mu[1]=_mm256_sub_epi16(W[1], predYUV0[1]);
						mu[2]=_mm256_sub_epi16(W[2], predYUV0[2]);
		//#if defined __GNUC__ && !defined PROFILER
		//#pragma GCC unroll 11
		//#endif
						for(int k=0;k<L1_NPREDS3;++k)//update
						{
							__m256i mc[6];
							mc[0]=_mm256_sign_epi16(L1preds[k*3+0], mu[0]);
							mc[1]=_mm256_sign_epi16(L1preds[k*3+1], mu[1]);
							mc[2]=_mm256_sign_epi16(L1preds[k*3+2], mu[2]);
							//16 -> 32	3 lo 3 hi registers
							mc[3]=_mm256_srai_epi32(mc[0], 16);
							mc[4]=_mm256_srai_epi32(mc[1], 16);
							mc[5]=_mm256_srai_epi32(mc[2], 16);
							mc[0]=_mm256_slli_epi32(mc[0], 16);
							mc[1]=_mm256_slli_epi32(mc[1], 16);
							mc[2]=_mm256_slli_epi32(mc[2], 16);
							mc[0]=_mm256_srai_epi32(mc[0], 16);
							mc[1]=_mm256_srai_epi32(mc[1], 16);
							mc[2]=_mm256_srai_epi32(mc[2], 16);
							mc[0]=_mm256_add_epi32(mc[0], _mm256_load_si256((__m256i*)L1weights+k*6+0));//update coeffs
							mc[1]=_mm256_add_epi32(mc[1], _mm256_load_si256((__m256i*)L1weights+k*6+1));
							mc[2]=_mm256_add_epi32(mc[2], _mm256_load_si256((__m256i*)L1weights+k*6+2));
							mc[3]=_mm256_add_epi32(mc[3], _mm256_load_si256((__m256i*)L1weights+k*6+3));
							mc[4]=_mm256_add_epi32(mc[4], _mm256_load_si256((__m256i*)L1weights+k*6+4));
							mc[5]=_mm256_add_epi32(mc[5], _mm256_load_si256((__m256i*)L1weights+k*6+5));
							_mm256_store_si256((__m256i*)L1weights+k*6+0, mc[0]);
							_mm256_store_si256((__m256i*)L1weights+k*6+1, mc[1]);
							_mm256_store_si256((__m256i*)L1weights+k*6+2, mc[2]);
							_mm256_store_si256((__m256i*)L1weights+k*6+3, mc[3]);
							_mm256_store_si256((__m256i*)L1weights+k*6+4, mc[4]);
							_mm256_store_si256((__m256i*)L1weights+k*6+5, mc[5]);
						}
					}
#endif

					//if(kb==1&&ky==35&&kx==83)//
					//	printf("");

					//context update = (2*eW+(e<<3)+max(eNEE, eNEEE))>>2
					eNEEE[0]=_mm256_load_si256((__m256i*)rows[1]+3+0+3*6);
					eNEEE[1]=_mm256_load_si256((__m256i*)rows[1]+3+1+3*6);
					eNEEE[2]=_mm256_load_si256((__m256i*)rows[1]+3+2+3*6);

					ecurr[0]=_mm256_slli_epi16(ecurr[0], GRBITS);
					ecurr[1]=_mm256_slli_epi16(ecurr[1], GRBITS);
					ecurr[2]=_mm256_slli_epi16(ecurr[2], GRBITS);
					ecurr[0]=_mm256_avg_epu16(ecurr[0], _mm256_max_epi16(eNEE[0], eNEEE[0]));
					ecurr[1]=_mm256_avg_epu16(ecurr[1], _mm256_max_epi16(eNEE[1], eNEEE[1]));
					ecurr[2]=_mm256_avg_epu16(ecurr[2], _mm256_max_epi16(eNEE[2], eNEEE[2]));
					eW[0]=_mm256_avg_epu16(eW[0], ecurr[0]);
					eW[1]=_mm256_avg_epu16(eW[1], ecurr[1]);
					eW[2]=_mm256_avg_epu16(eW[2], ecurr[2]);

					_mm256_store_si256((__m256i*)rows[0]+3+0+0*6, eW[0]);//store current contexts
					_mm256_store_si256((__m256i*)rows[0]+3+1+0*6, eW[1]);
					_mm256_store_si256((__m256i*)rows[0]+3+2+0*6, eW[2]);
					eNEE[0]=eNEEE[0];
					eNEE[1]=eNEEE[1];
					eNEE[2]=eNEEE[2];
					NW[0]=N[0];
					NW[1]=N[1];
					NW[2]=N[2];
					rows[0]+=6*NLANES;
					rows[1]+=6*NLANES;
					rows[2]+=6*NLANES;
					rows[3]+=6*NLANES;
					imptr+=3*NLANES;
				}
			}
		}

		//encode/deinterleave
		if(fwd)//encode main
		{
			for(int kc=0;kc<3*NCTX;++kc)
				enc_hist2stats(hists+(ptrdiff_t)256*kc, encstats+(ptrdiff_t)256*kc, &ctxmask, kc);
			mstate[1]=mstate[0]=_mm256_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
			uint16_t *ctxptr2=ctxbuf+(ptrdiff_t)3*NLANES*iblockw*iblockh-NLANES;
			for(int ky=iblockh-1;ky>=0;--ky)
			{
#ifdef ESTIMATE_SIZE
				int kc=2;
#endif
				for(int kx=3*iblockw-1;kx>=0;--kx)
				{
					__m256i mmax[2], minvf[2], mcdf[2], mnegf_sh[2];
#ifdef _DEBUG
					if(ctxptr2<ctxbuf)
						CRASH("OOB");
#endif
					//if(!kb&&!ky&&!kx)//
					//if(!kb&&!ky&&kx==2*3+1)//
					//if(!kb&&!ky&&kx==2*3+2)//
					//if(!kb&&!ky&&kx==2*3+0)//
					//	printf("");

					{
						__m256i s0, s1, s2, s3;
						__m256i t0, t1, t2, t3;
#define SHUFFLE_PS(LO, HI, IMM8_HHLL) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(LO), _mm256_castsi256_ps(HI), IMM8_HHLL))

						s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(encstats+ctxptr2[0*8+0])));
						s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(encstats+ctxptr2[0*8+1])));
						s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(encstats+ctxptr2[0*8+2])));
						s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(encstats+ctxptr2[0*8+3])));
						s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)(encstats+ctxptr2[0*8+4])), 1);
						s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)(encstats+ctxptr2[0*8+5])), 1);
						s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)(encstats+ctxptr2[0*8+6])), 1);
						s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)(encstats+ctxptr2[0*8+7])), 1);
						t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));//_MM_TRANSPOSE4_PS
						t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
						t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
						t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
						mmax	[0]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
						minvf	[0]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
						mcdf	[0]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
						mnegf_sh[0]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));

						s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(encstats+ctxptr2[1*8+0])));
						s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(encstats+ctxptr2[1*8+1])));
						s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(encstats+ctxptr2[1*8+2])));
						s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(encstats+ctxptr2[1*8+3])));
						s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)(encstats+ctxptr2[1*8+4])), 1);
						s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)(encstats+ctxptr2[1*8+5])), 1);
						s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)(encstats+ctxptr2[1*8+6])), 1);
						s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)(encstats+ctxptr2[1*8+7])), 1);
						t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));
						t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
						t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
						t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
						mmax	[1]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
						minvf	[1]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
						mcdf	[1]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
						mnegf_sh[1]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));

						ctxptr2-=NLANES;
					}
#ifdef ESTIMATE_SIZE
					{
						ALIGN(32) int anegf[NLANES]={0};
						memcpy(anegf, mnegf_sh, sizeof(anegf));
						const double norm=1./(1<<PROBBITS);
						for(int k=0;k<NLANES;++k)
						{
							int freq=(1<<PROBBITS)-(anegf[k]&0xFFFF);
							if((unsigned)(freq-1)>=(unsigned)((1<<PROBBITS)-1))
								CRASH("freq = %d", freq);
							esize[kc*NLANES+k]-=log2(freq*norm)*0.125;
						}
					}
					--kc;
					if(kc<0)
						kc=2;
#endif

					//enc renorm
					__m256i cond0=_mm256_cmpgt_epi32(mstate[0], mmax[0]);
					__m256i cond1=_mm256_cmpgt_epi32(mstate[1], mmax[1]);
					int mask0=_mm256_movemask_ps(_mm256_castsi256_ps(cond0));
					int mask1=_mm256_movemask_ps(_mm256_castsi256_ps(cond1));
					__m256i idx0=_mm256_load_si256((__m256i*)ans_permute[mask0]);
					__m256i idx1=_mm256_load_si256((__m256i*)ans_permute[mask1]);
					__m256i emit0=_mm256_permutevar8x32_epi32(_mm256_and_si256(mstate[0], cond0), idx0);
					__m256i emit1=_mm256_permutevar8x32_epi32(_mm256_and_si256(mstate[1], cond1), idx1);
					emit0=_mm256_and_si256(emit0, _mm256_set1_epi32(0xFFFF));
					emit1=_mm256_and_si256(emit1, _mm256_set1_epi32(0xFFFF));
					emit0=_mm256_packus_epi32(emit0, emit1);
					emit0=_mm256_permute4x64_epi64(emit0, _MM_SHUFFLE(3, 1, 2, 0));
					__m128i e1=_mm256_extractf128_si256(emit0, 1);
					__m128i e0=_mm256_castsi256_si128(emit0);
					mask1=_mm_popcnt_u32(mask1);
					mask0=_mm_popcnt_u32(mask0);
#ifdef _DEBUG
					if(streamptr-(2*((ptrdiff_t)mask0+mask1)+sizeof(__m128i))<=streamstart)
						CRASH("OOB ptr %016zX <= %016zX", streamptr, streamstart);
#endif
					_mm_storeu_si128((__m128i*)streamptr-1, e1); streamptr-=mask1*sizeof(short);
					_mm_storeu_si128((__m128i*)streamptr-1, e0); streamptr-=mask0*sizeof(short);
					{
						__m256i state0=_mm256_srli_epi32(mstate[0], 16);
						__m256i state1=_mm256_srli_epi32(mstate[1], 16);
						mstate[0]=_mm256_blendv_epi8(mstate[0], state0, cond0);
						mstate[1]=_mm256_blendv_epi8(mstate[1], state1, cond1);
					}
#ifdef ANS_VAL
					ansval_push(mstate, sizeof(int), NLANES);
#endif
					//enc update		state += (state*invf>>32>>sh)*negf+cdf		state = state/freq<<12|(cdf+state%freq)
					{
						__m256i lo0=_mm256_mul_epu32(mstate[0], minvf[0]);//q = mulhi32(state, invf)
						__m256i lo1=_mm256_mul_epu32(mstate[1], minvf[1]);
						__m256i hi0=_mm256_mul_epu32(_mm256_srli_epi64(mstate[0], 32), _mm256_srli_epi64(minvf[0], 32));
						__m256i hi1=_mm256_mul_epu32(_mm256_srli_epi64(mstate[1], 32), _mm256_srli_epi64(minvf[1], 32));
						minvf[0]=_mm256_blend_epi32(_mm256_srli_epi64(lo0, 32), hi0, 0xAA);
						minvf[1]=_mm256_blend_epi32(_mm256_srli_epi64(lo1, 32), hi1, 0xAA);
					}
					{
						__m256i sh0=_mm256_srli_epi32(mnegf_sh[0], 16);
						__m256i sh1=_mm256_srli_epi32(mnegf_sh[1], 16);
						minvf[0]=_mm256_srlv_epi32(minvf[0], sh0);
						minvf[1]=_mm256_srlv_epi32(minvf[1], sh1);
					}
					mstate[0]=_mm256_add_epi32(mstate[0], mcdf[0]);
					mstate[1]=_mm256_add_epi32(mstate[1], mcdf[1]);
					{
						__m256i lomask=_mm256_set1_epi32(0xFFFF);
						__m256i negf0=_mm256_and_si256(mnegf_sh[0], lomask);
						__m256i negf1=_mm256_and_si256(mnegf_sh[1], lomask);
						minvf[0]=_mm256_mullo_epi32(minvf[0], negf0);
						minvf[1]=_mm256_mullo_epi32(minvf[1], negf1);
#ifdef ANS_VAL
						__m256i mM=_mm256_set1_epi32(1<<PROBBITS);
						negf0=_mm256_sub_epi16(mM, negf0);
						negf1=_mm256_sub_epi16(mM, negf1);
						negf0=_mm256_packus_epi32(negf0, negf1);
						negf0=_mm256_permute4x64_epi64(negf0, _MM_SHUFFLE(3, 1, 2, 0));
						ALIGN(32) unsigned short freqs[NLANES];
						_mm256_store_si256((__m256i*)freqs, negf0);
						ansval_push(freqs, sizeof(short), NLANES);
#endif
					}
					mstate[0]=_mm256_add_epi32(mstate[0], minvf[0]);
					mstate[1]=_mm256_add_epi32(mstate[1], minvf[1]);
#ifdef ANS_VAL
					ansval_push(mstate, sizeof(int), NLANES);
#endif
				}
			}
			//flush
			streamptr-=sizeof(mstate);
#ifdef _DEBUG
			if(streamptr<=streamstart)
				CRASH("OOB ptr %016zX <= %016zX", streamptr, streamstart);
#endif
			memcpy(streamptr, mstate, sizeof(mstate));
			//pack hists
			{
				BitPackerLIFO ec;

				bitpacker_enc_init(&ec, image, streamptr);
				for(int kc=3*NCTX-1;kc>=0;--kc)
					enc_packhist(&ec, hists+(ptrdiff_t)256*kc, ctxmask, kc);
				bitpacker_enc_flush(&ec);
				streamptr=ec.dstbwdptr;

				streamptr-=8;
				*(unsigned long long*)streamptr=ctxmask;

				streamptr--;
				*streamptr=bestrct;
			}
		}
		else//deinterleave
		{
			uint8_t *slowptrs0[NLANES], *slowptrs[NLANES];
			const uint8_t *fastptr;
			int kx, ky, k;

			fastptr=block;
			for(ky=0;ky<YLANES;++ky)//spread slow pointers
			{
				for(kx=0;kx<XLANES;++kx)
					slowptrs0[XLANES*ky+kx]=image+3*(iw*(y1+iblockh*ky)+x1+iblockw*kx);
			}
			for(ky=0;ky<iblockh;++ky)//interleave
			{
				memcpy(slowptrs, slowptrs0, sizeof(slowptrs));
				for(kx=0;kx<iblockw;++kx)
				{
					for(k=0;k<NLANES;++k)
						*slowptrs[k]++=*fastptr++;
					for(k=0;k<NLANES;++k)
						*slowptrs[k]++=*fastptr++;
					for(k=0;k<NLANES;++k)
						*slowptrs[k]++=*fastptr++;
				}
				for(k=0;k<NLANES;++k)
					slowptrs0[k]+=rowstride;
			}
		}
		if(fwd)
			bsizes[kb]=(int)(streamptr0-streamptr);
	}
	if(fwd)
	{
		arg->chunk=streamptr;
		arg->clen=streamend-streamptr;
	}
	_mm_free(cbuf);
	if(fwd)
	{
		free(hists);
		_mm_free(encstats);
	}
	else
		_mm_free(CDF2syms);
	return 0;
}
static int partition_dim(int qdim, int maxblocksize, uint16_t *stops)
{
	int nblocks, prev;
	
	nblocks=0;
	prev=0;
	stops[nblocks++]=0;
	while((prev+=maxblocksize)<qdim)stops[nblocks++]=prev;
	stops[nblocks]=qdim;
	return nblocks;
}
int c34_codec(int argc, char **argv)
{
#ifdef LOUD
	double t=time_sec();
	double t_main;
#endif
	const char *srcfn=0, *dstfn=0;
	int ncores;

	ptrdiff_t streamsize=0;
	unsigned char *streamstart=0, *streamend=0;

	int argssize;
	ThreadArgs *args=0;

	ptrdiff_t csize=0;
	uint8_t *rbuf=0, *rstream=0;
	ptrdiff_t rsize;

	prof_start();

	(void)streamstart;
	(void)streamend;
	(void)memfill;
	(void)time_sec;
	(void)och_names;
	(void)rct_names;
	if(argc<3)
	{
		printf(
			"Usage:  \"%s\"  input  output  [-d distortion  -m nthreads  -p]\n"
		//	"Usage:  \"%s\"  input  output  [-e effort  -d distortion  -m nthreads  -p]\n"
		//	"  -e effort      0: No prediction    1: CG    2: L1A    3: L1B (default)\n"
			"  -d distortion  0: Lossless (default)    otherwise 4 <= dist. <= 16\n"
			"  -m nthreads    0: MT (default)    1: ST    2...: specify number of threads\n"
			"  -p             Performance profiler\n"
			, argv[0]
		);
		return 1;
	}
	srcfn=argv[1];
	dstfn=argv[2];
	nthreads=0, effort=3, dist=0, profiler=0;

	//parse command args
	{
		int k;

		for(k=3;k<argc;++k)
		{
			const char *arg=argv[k];
			if(arg[0]!='-'||!arg[0]||(arg[1]&&!arg[2]&&k+1>=argc&&(arg[1]&0xDF)!='P'))
			{
				CRASH("Invalid argument  \"%s\"", argv[k]);
				return 1;
			}
			switch(arg[1]&0xDF)
			{
			case 'E':
				if(arg[2])
					effort=atoi(arg+2);
				else
					effort=atoi(argv[++k]);
				if(effort<0||effort>3)
				{
					CRASH("Invalid argument  effort %d  expected from {0, ...3}", effort);
					return 1;
				}
				break;
			case 'D':
				if(arg[2])
					dist=atoi(arg+2);
				else
					dist=atoi(argv[++k]);
				if(dist&&(dist<4||dist>16))
				{
					CRASH("Invalid argument  dist. %d  expected 0 or from {4, ...16}", dist);
					return 1;
				}
				break;
			case 'M':
				if(arg[2])
					nthreads=atoi(arg+2);
				else
					nthreads=atoi(argv[++k]);
				if(nthreads<0)
				{
					CRASH("Invalid argument  nthreads %d  expected from {0=MT, 1=ST, 2...=MT}", nthreads);
					return 1;
				}
				break;
			case 'P':
				profiler=1;
				break;
			}
		}
	}
	ncores=query_cpu_cores();
	if(!nthreads||nthreads>ncores)
		nthreads=ncores;
	{
		struct stat info={0};
		int error=stat(srcfn, &info);
		if(error)
		{
			CRASH("Cannot stat \"%s\"", srcfn);
			return 1;
		}
		srcsize=info.st_size;
	}

	//read source
	{
		ptrdiff_t nread;
		int c;
		FILE *fsrc;
		
		fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			CRASH("Cannot open \"%s\"", srcfn);
			return 1;
		}
		c=0;
		fread(&c, 1, 2, fsrc);
		fwd=c==('P'|'6'<<8);
		if(!fwd&&c!=('3'|'4'<<8))
		{
			CRASH("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
		{
			c=fgetc(fsrc);
			if(c!='\n')
			{
				CRASH("Invalid PPM file");
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
			while((unsigned)(c-'0')<10)
			{
				iw=10*iw+c-'0';
				c=fgetc(fsrc);
			}
			while(c<=' ')
				c=fgetc(fsrc);
			ih=0;
			while((unsigned)(c-'0')<10)
			{
				ih=10*ih+c-'0';
				c=fgetc(fsrc);
			}
			while(c<=' ')
				c=fgetc(fsrc);
			while(c=='#')
			{
				c=fgetc(fsrc);
				while(c!='\n')
					c=fgetc(fsrc);
				c=fgetc(fsrc);
			}
			c=c<<8|fgetc(fsrc);
			c=c<<8|fgetc(fsrc);
			c=c<<8|fgetc(fsrc);
			if(c!=('2'<<24|'5'<<16|'5'<<8|'\n'))
			{
				CRASH("Unsupported PPM file");
				return 1;
			}
		}
		else
		{
			fread(&iw, 1, 4, fsrc);
			fread(&ih, 1, 4, fsrc);
			dist=0;
			fread(&dist, 1, 1, fsrc);

			nread=ftell(fsrc);
			streamsize=srcsize-nread;

		}
		if(iw<1||ih<1||iw>0x8000||ih>0x8000)
		{
			CRASH("Unsupported image dimensions  WH %d*%d", iw, ih);
			return 1;
		}
		rowstride=3*iw;
		res=(ptrdiff_t)iw*ih;
		usize=3*res;
		image=(unsigned char*)malloc(usize);
		if(!image)
		{
			CRASH("Alloc error");
			return 1;
		}
		{
			ptrdiff_t expected=0;
			if(fwd)
			{
				expected=usize;
				nread=fread(image, 1, usize, fsrc);
				guide_save(image, iw, ih);
			}
			else
			{
				stream=(unsigned char*)malloc(streamsize+sizeof(__m256i));
				if(!stream)
				{
					CRASH("Alloc error");
					return 1;
				}
				expected=streamsize;
				nread=fread(stream, 1, streamsize, fsrc);
			}
			if(nread!=expected)
				printf("Truncated  expected %td  read %td", expected, nread);
		}
		fclose(fsrc);
	}
	memset(ans_permute, 0, sizeof(ans_permute));
	if(fwd)
	{
		//generate encode permutations		eg: mask = MSB 0b11000101 LSB  ->  LO {x, x, x, x, 0, 2, 6, 7} HI
		for(int km=0;km<256;++km)
		{
			int *curr=ans_permute[km];
			int kb2=7;
			for(int kb=7;kb>=0;--kb)
			{
				int bit=km>>kb&1;
				if(bit)
					curr[kb2--]=kb;
			}
		}
	}
	else
	{
		//generate decode permutations		eg: mask = MSB 0b11000101 LSB  ->  LO {0, x, 1, x, x, x, 2, 3} HI
		for(int km=0;km<256;++km)
		{
			int *curr=ans_permute[km];
			int idx=0;
			for(int kb=0;kb<8;++kb)
			{
				int bit=km>>kb&1;
				if(bit)
					curr[kb]=idx++;
			}
		}
		streamstart=stream;
		streamend=stream+streamsize;
	}
	xrem=iw&3;
	yrem=ih&3;
	qw=iw&~3;
	qh=ih&~3;
	xnblocks=partition_dim(qw, XMAXBLOCK, blockxstops);
	ynblocks=partition_dim(qh, YMAXBLOCK, blockystops);
	nblocks=xnblocks*ynblocks;
	if(!fwd)
	{
		memcpy(bsizes, stream, sizeof(int)*nblocks);
		//streamptr+=sizeof(int)*nblocks;

#ifdef PRING_BLOCKSIZES
		for(int k=0;k<nblocks;++k)
			printf("%3d %7d\n", k, bsizes[k]);
#endif
	}
	if(nthreads>nblocks)
		nthreads=nblocks;
#ifdef DISABLE_MT
	nthreads=1;
#endif
	argssize=sizeof(ThreadArgs)*nthreads;
	args=(ThreadArgs*)malloc(argssize);
	if(!args)
	{
		CRASH("Alloc error");
		return 1;
	}
	memset(args, 0, argssize);
	{
		int k, kb;

		for(k=0, kb=0;k<nthreads;++k)
		{
			ThreadArgs *arg=args+k;
			arg->threadidx=k;
			arg->b1=kb;
			kb=(k+1)*nblocks/nthreads;
			arg->b2=kb;
		}
	}
	PROF(read, srcsize);

	//decode remainders
	if(!fwd&&(xrem||yrem))
	{
		const int CDF2symsize=(int)sizeof(int[3<<PROBBITS]);
		uint32_t *CDF2sym;
		int rowstride=iw*3;
		int rct=7;
		uint64_t ctxmask=-1;
		uint8_t *streamptr;

		CDF2sym=(uint32_t*)malloc(CDF2symsize);
		if(!CDF2sym)
		{
			CRASH("Alloc error");
			return 1;
		}
		streamptr=stream+nblocks*sizeof(int);
		for(int k=0;k<nblocks;++k)
			streamptr+=bsizes[k];
		{
			BitPackerLIFO ec;

			ctxmask=*streamptr++;
			bitpacker_dec_init(&ec, streamptr, streamend);
			dec_unpackhist(&ec, CDF2sym+((ptrdiff_t)0<<PROBBITS), ctxmask, 3*NCTX+0);
			dec_unpackhist(&ec, CDF2sym+((ptrdiff_t)1<<PROBBITS), ctxmask, 3*NCTX+1);
			dec_unpackhist(&ec, CDF2sym+((ptrdiff_t)2<<PROBBITS), ctxmask, 3*NCTX+2);
			streamptr=(unsigned char*)(size_t)ec.srcfwdptr;
		}

#ifdef _DEBUG
		if(streamptr>streamend)
			CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		{
			unsigned state=*(unsigned*)streamptr;
			streamptr+=4;
			for(int ky=0;ky<yrem;++ky)
				decode1d(image+rowstride*(qh+ky), iw, 3, rct, &state, (const unsigned char**)&streamptr, streamend, CDF2sym);
			for(int kx=0;kx<xrem;++kx)
				decode1d(image+3*(qw+kx), qh, rowstride, rct, &state, (const unsigned char**)&streamptr, streamend, CDF2sym);
		}
		free(CDF2sym);
		PROF(rem, 3*(xrem*ih+yrem*iw));
	}
#ifdef LOUD
	t_main=time_sec();
#endif
#ifndef DISABLE_MT
	if(nthreads>1)
	{
		void *handles[128];
		int k;

		for(k=0;k<nthreads;++k)
		{
			int error;

#ifdef _WIN32
			handles[k]=(void*)_beginthreadex(0, 0, c34_thread, args+k, 0, 0);
			error=!handles[k];
#else
			error=pthread_create((pthread_t*)&handle[k], 0, c34_thread, args+k);
#endif
			if(error)
			{
				CRASH("Alloc error");
				return 1;
			}
		}
#ifdef _WIN32
		WaitForMultipleObjects(nthreads, (HANDLE*)handles, TRUE, INFINITE);
		for(int k=0;k<nthreads;++k)
			CloseHandle(handles[k]);
#else
		for(int k=0;k<nthreads;++k)
			pthread_join((pthread_t)handles[k], 0);
#endif
	}
	else
#endif
	{
		c34_thread(args);
	}
#ifdef LOUD
	t_main=time_sec()-t_main;
#endif
	PROF(main, usize);

	//encode remainders
	rsize=0;
	if(fwd&&(xrem||yrem))
	{
		const int hsize=(int)sizeof(int[3*256]);
		const int encstatsize=(int)sizeof(rANS_SIMD_SymInfo[3<<8]);
		int32_t *hist;
		rANS_SIMD_SymInfo *encstats;
		int rowstride=iw*3;
		int rct=7;
		uint64_t ctxmask=0;
		int rcap;
		uint8_t *streamptr;

		rcap=4*(xrem*ih+yrem*iw);
		rbuf=(uint8_t*)malloc(rcap);
		hist=(int*)malloc(hsize);
		encstats=(rANS_SIMD_SymInfo*)malloc(encstatsize);
		if(!rbuf||!hist||!encstats)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(hist, 0, hsize);

		for(int ky=0;ky<yrem;++ky)//remainder rows
			decorr1d(image+rowstride*(qh+ky), iw, 3, rct, hist);
		for(int kx=0;kx<xrem;++kx)//remainder columns
			decorr1d(image+3*(qw+kx), qh, rowstride, rct, hist);
			
		enc_hist2stats(hist+(ptrdiff_t)256*0, encstats+(ptrdiff_t)256*0, &ctxmask, 0);
		enc_hist2stats(hist+(ptrdiff_t)256*1, encstats+(ptrdiff_t)256*1, &ctxmask, 1);
		enc_hist2stats(hist+(ptrdiff_t)256*2, encstats+(ptrdiff_t)256*2, &ctxmask, 2);

		streamstart=rbuf;
		streamptr=streamend=rbuf+rcap;
		{
			unsigned state=1<<(RANS_STATE_BITS-RANS_RENORM_BITS);
			for(int kx=xrem-1;kx>=0;--kx)
				encode1d(image+3*(qw+kx), qh, rowstride, &state, &streamptr, streamstart, encstats);
			for(int ky=yrem-1;ky>=0;--ky)
				encode1d(image+rowstride*(qh+ky), iw, 3, &state, &streamptr, streamstart, encstats);
			//flush
			streamptr-=4;
#ifdef _DEBUG
			if(streamptr<=rbuf)
				CRASH("OOB ptr %016zX <= %016zX", streamptr, streamstart);
#endif
			*(unsigned*)streamptr=state;
		}
		{
			BitPackerLIFO ec;
			bitpacker_enc_init(&ec, streamstart, streamptr);
			enc_packhist(&ec, hist+256*2, ctxmask, 2);
			enc_packhist(&ec, hist+256*1, ctxmask, 1);
			enc_packhist(&ec, hist+256*0, ctxmask, 0);
			bitpacker_enc_flush(&ec);
			streamptr=ec.dstbwdptr;
		}
		--streamptr;
		*streamptr=(uint8_t)ctxmask;

		rstream=streamptr;
		rsize=streamend-streamptr;

		free(hist);
		free(encstats);
		PROF(rem, 3*(xrem*ih+yrem*iw));
	}

	//write output
	{
		ptrdiff_t dstsize=0;
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing", dstfn);
			return 1;
		}
		if(fwd)
		{
			csize+=fwrite("34", 1, 2, fdst);
			csize+=fwrite(&iw, 1, 4, fdst);
			csize+=fwrite(&ih, 1, 4, fdst);
			csize+=fwrite(&dist, 1, 1, fdst);
			csize+=fwrite(bsizes, 1, nblocks*sizeof(int), fdst);
			for(int kt=0;kt<nthreads;++kt)
			{
				ThreadArgs *arg=args+kt;
				csize+=fwrite(arg->chunk, 1, arg->clen, fdst);
				//uint8_t *streamptr=arg->chunk+arg->clen;
				//for(int kb=arg->b1;kb<arg->b2;++kb)//reverse block order
				//{
				//	streamptr-=bsizes[kb];
				//	csize+=fwrite(streamptr, 1, bsizes[kb], fdst);
				//}
			}
			csize+=fwrite(rstream, 1, rsize, fdst);

			dstsize=csize;

			for(int kt=0;kt<nthreads;++kt)
				free(args[kt].buffer);
			free(rbuf);
		}
		else
		{
			csize=srcsize;
			dstsize+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			dstsize+=fwrite(image, 1, usize, fdst);
		}
		fclose(fdst);
		PROF(write, dstsize);
		(void)dstsize;
	}
#if 0
	if(fwd)
	{
		uint8_t header[128]={0};
		FILE_SEGMENT_ELEMENT chunks[128]=
		{
			header, 0,
		};
		HANDLE hdst=CreateFileA(dstfn, GENERIC_WRITE, 0, 0, CREATE_ALWAYS, FILE_FLAG_OVERLAPPED|FILE_FLAG_NO_BUFFERING, 0);
		if(hdst==INVALID_HANDLE_VALUE)
		{
			CRASH("Cannot open \"%s\" for writing", dstfn);
			return 1;
		}
	}
	else
	{
	}
#endif
	free(args);
	free(image);
	if(!fwd)
		free(stream);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		//printf("Mem usage: ");
		//print_size((double)memusage, 8, 4, 0, 0);
		//printf("\n");
		printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
	}
	printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t_main, usize/(t_main*1024*1024));//
	printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t, usize/(t*1024*1024));
#ifdef ENABLE_GUIDE
	if(!fwd&&dist>=4)
	{
		double rmse=sqrt((g_sqe[0]+g_sqe[1]+g_sqe[2])/usize), psnr=20*log10(255/rmse);
		printf("RMSE %12.6lf  PSNR %12.6lf\n", rmse, psnr);
	}
	if(!fwd)
		free(g_image);
#endif
#else
	(void)csize;
#endif
#ifndef _MSC_VER
	if(profiler)
#endif
		prof_print();
	return 0;
}
