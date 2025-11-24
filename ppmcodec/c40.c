#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<math.h>
#include<time.h>
#include<immintrin.h>
#include<sys/stat.h>
#if defined _WIN32 || defined WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#endif


	#define PROFILE_TIME		//should be on

#ifdef _MSC_VER
	#define PROFILE_SIZE
	#define LOUD			//size & time

//	#define ESTIMATE_SIZE2
//	#define ESTIMATE_SIZE		//DEBUG		checks for zero frequency, visualizes context usage
	#define ENABLE_GUIDE		//DEBUG		checks interleaved pixels
//	#define ANS_VAL			//DEBUG

//	#define SAVE_RESIDUALS
//	#define PRINT_L1_BOUNDS
//	#define TEST_INTERLEAVE
#endif

//	#define USE_FSE
	#define ENABLE_RCT_EXTENSION
//	#define EMULATE_GATHER		//gather is a little faster
	#define INTERLEAVESIMD		//2.5x faster interleave


#define ANALYSIS_XSTRIDE 2
#define ANALYSIS_YSTRIDE 2

#define DEFAULT_EFFORT_LEVEL 2
#define L1_NPREDS1 4
#define L1_NPREDS2 8
#define L1_NPREDS3 20
#define L1_SH1 15	//L1_SH1 <= 16
#define L1_SH2 17
#define L1_SH3 19

//3*17+3=54 contexts
#define GRBITS 3
#define NCTX 18		//18*3+3 = 57 total

#define XCODERS 4	//xrem 1~3 cols
#define YCODERS 4	//yrem 1~3 rows

#define NCODERS 16

#define PROBBITS 12	//12 bit max	James Bonfield's CDF2sym: {freq<<20 | bias<<8 | sym}

#define RANS_STATE_BITS 31
#define RANS_RENORM_BITS 16

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
#if defined _WIN32 || defined WIN32
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
static unsigned char *g_image=0;
static double g_sqe[3]={0};
static void guide_save(const unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		CRASH("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(const unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		CRASH("");
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif
#ifdef PROFILE_SIZE
#include<stdarg.h>
static void profile_size(const unsigned char *dstbwdptr, const char *msg, ...)
{
	static ptrdiff_t size=0;
	static const unsigned char *prev=0;
	if(prev)
	{
		ptrdiff_t diff=prev-dstbwdptr;
		size+=diff;
		printf("%10td (%+10td) bytes", size, diff);
		if(msg)
		{
			va_list args;
			va_start(args, msg);
			printf("  ");
			vprintf(msg, args);
			va_end(args);
		}
		printf("\n");
	}
	prev=dstbwdptr;
}
#else
#define profile_size(...)
#endif
#ifdef PROFILE_TIME
typedef struct _SpeedProfilerInfo
{
	double t;
	ptrdiff_t size;
	const char *msg;
} SpeedProfilerInfo;
#define PROF_CAP 128
static double prof_timestamp=0;
static SpeedProfilerInfo prof_data[PROF_CAP]={0};
static int prof_count=0;
static void prof_checkpoint(ptrdiff_t size, const char *msg)
{
	double t2=time_sec();
	if(prof_timestamp)
	{
		SpeedProfilerInfo *info=prof_data+prof_count++;
		if(prof_count>PROF_CAP)
		{
			CRASH("Profiler OOB");
			return;
		}
		info->t=t2-prof_timestamp;
		info->size=size;
		info->msg=msg;
		//double delta=t2-t;
		//printf("%16.12lf sec", delta);
		//if(size)
		//	printf(" %12.6lf MB/s %10td bytes", size/(delta*1024*1024), size);
		//if(msg)
		//	printf(" %s", msg);
		//printf("\n");
	}
	prof_timestamp=t2;
}
static void prof_print(ptrdiff_t usize)
{
	static char buf[2048]={0};
	double timesum=0, tmax=0;
	for(int k=0;k<prof_count;++k)
	{
		double t=prof_data[k].t;
		timesum+=t;
		if(tmax<t)
			tmax=t;
	}
	int prev=0;
	double csum=0;
	//int colors[128]={0};
	//srand((unsigned)__rdtsc());
	//colorgen(colors, prof_count, 64, 300, 100);
	const int scale=5;
	printf("1 char = %d ms\n", scale);
	printf("|");
	for(int k=0;k<prof_count;++k)
	{
		SpeedProfilerInfo *info=prof_data+k;
		csum+=info->t;
		int curr=(int)(csum*1000/scale);//fixed scale
		int space=curr-prev;
		int len=0;
		if(info->msg)
			len=(int)strlen(info->msg);
		if(space>2047)//printf("%*s", HUGE, ""); CRASHES
			space=2047;
		if(info->msg&&space>=len)
		{
			int labelstart=(space-len)>>1;
			int labelend=labelstart+len;

			memset(buf, '-', labelstart);
			buf[labelstart]=0;
			printf("%s%s", buf, info->msg);
			//colorprintf(colors[k], colors[k], buf);
			//colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%s", info->msg);
			memset(buf, '-', (ptrdiff_t)space-labelend);
			buf[space-labelend]=0;
			printf("%s", buf);
			//colorprintf(colors[k], colors[k], buf);
		}
		else
		{
			memset(buf, '-', space);
			buf[space]=0;
			printf("%s", buf);
			//colorprintf(colors[k], colors[k], buf);
		}
		printf("|");
		prev=curr;
	}
	printf("\n");
	for(int k=0;k<prof_count;++k)
	{
		SpeedProfilerInfo *info=prof_data+k;
		printf("%16.7lf ms %8.4lf%% ", info->t*1000, 100.*info->t/timesum);
		if(info->size)
			printf(" %16.6lf MB/s %16.6lf ms/MB %10td bytes "
				, info->size/(info->t*1024*1024)
				, info->t*1024*1024*1000/info->size
				, info->size
			);
		if(info->msg)
			printf("%s", info->msg);
		else
			printf("%d", k);
		//if(info->msg)
		//	colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%-20s", info->msg);
		//else
		//	colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%-20s", "");
		printf("\n");
	}
	printf("\n");
	printf("%16.7lf ms %12.6lf MB/s Total\n", timesum*1000, usize/(timesum*1024*1024));
	printf("\n");
	prof_count=0;
	prof_timestamp=0;
}
#else
#define prof_checkpoint(...)
#define prof_print(...)
#endif

//cRCT
#if 1
#ifndef ENABLE_RCT_EXTENSION
#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
	OCH(CX40) OCH(C0X4) OCH(C40X)
#endif
#ifdef ENABLE_RCT_EXTENSION
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
#endif
typedef enum _OCHIndex
{
#define OCH(X) OCH_##X,
	OCHLIST
#undef  OCH
	OCH_COUNT,

	OCH_R=OCH_Y400,
	OCH_G=OCH_Y040,
	OCH_B=OCH_Y004,
	OCH_C4X0=OCH_CX40,
	OCH_C04X=OCH_C0X4,
	OCH_CX04=OCH_C40X,
	OCH_BG=OCH_C04X,
	OCH_BR=OCH_C40X,
	OCH_RG=OCH_CX40,
	OCH_RB=OCH_CX04,
	OCH_GB=OCH_C0X4,
	OCH_GR=OCH_C4X0,
#ifdef ENABLE_RCT_EXTENSION
	OCH_R1=OCH_CX13,
	OCH_G1=OCH_C3X1,
	OCH_B1=OCH_C13X,
	OCH_R2=OCH_CX22,
	OCH_G2=OCH_C2X2,
	OCH_B2=OCH_C22X,
	OCH_R3=OCH_CX31,
	OCH_G3=OCH_C1X3,
	OCH_B3=OCH_C31X,
#endif
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
//	II_COEFF_U_SUB_V2,
//	II_COEFF_V_SUB_U2,

	II_COUNT,
} RCTInfoIdx;
//YUV = RCT * RGB	watch out for permutation in last row
//luma: averaging	chroma: subtraction
//example: _400_40X_3X1 == [1 0 0; -1 0 1; -3/4 1 -1/4]
#ifndef ENABLE_RCT_EXTENSION
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
	RCT(_400_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)
#endif
#ifdef ENABLE_RCT_EXTENSION
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
#endif

#ifdef ANS_VAL
#define ANS_VAL_HISTSIZE 128
typedef struct _ANSVALHeader
{
	unsigned short esize, count;
	unsigned idx;
	struct _ANSVALHeader *above, *below;
	unsigned char data[];
} ANSVALNode;
static ANSVALNode *debugstack=0;
static int ansvalidx=0, ansvalmax=0;
static void ansval_push(const void *data, int esize, int count)
{
	int size=count*esize;
	ANSVALNode *node=(ANSVALNode*)malloc(sizeof(ANSVALNode)+size);
	if(!node)
	{
		printf("Alloc error\n");
		exit(1);
	}
	memset(node, 0, sizeof(ANSVALNode)+size);
	node->esize=esize;
	node->count=count;
	node->idx=ansvalidx++;
	node->above=0;
	node->below=debugstack;
	if(debugstack)
		debugstack->above=node;
	memcpy(node->data, data, size);
	debugstack=node;
	++ansvalmax;
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
		printf("\n\nValidation Error  [enc ^ | v dec]\n");
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
#else
#define ansval_push(...)
#define ansval_check(...)
#endif


//LIFO Bypass Coder
#define BITPACKERMAX 32
typedef struct _BitPackerLIFO//bwd enc / fwd dec
{
	uint64_t state;
	int32_t enc_nwritten, dec_navailable;//bitcounts, only for tracking renorms
	uint8_t *dstbwdptr;
	const uint8_t *srcfwdptr, *streamend;
} BitPackerLIFO;
AWM_INLINE void bitpacker_enc_init(BitPackerLIFO *ec, const uint8_t *bufstart, uint8_t *bufptr0_OOB)
{
	memset(ec, 0, sizeof(*ec));
	ec->state=1ULL<<32;
	ec->enc_nwritten=33;
	ec->streamend=bufstart;
	ec->dstbwdptr=bufptr0_OOB;
}
AWM_INLINE void bitpacker_dec_init(BitPackerLIFO *ec, const uint8_t *bufptr0_start, const uint8_t *bufend)
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
	if(inbits>BITPACKERMAX)
		CRASH("BitPacker inbits %d", inbits);
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
	ansval_push(&inbits, sizeof(inbits), 1);
	ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
}
AWM_INLINE int bitpacker_dec(BitPackerLIFO *ec, int outbits)
{
#ifdef _DEBUG
	if(outbits>BITPACKERMAX)
		CRASH("BitPacker outbits %d", outbits);
#endif
	int sym=ec->state&((1ULL<<outbits)-1);

	//pop outbits then renorm
#ifdef ANS_VAL
	ansval_check(&ec->state, sizeof(ec->state), 1);
	ansval_check(&outbits, sizeof(outbits), 1);
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
		if(ec->srcfwdptr>=ec->streamend)
			CRASH("IntPacker OOB:  srcfwdptr = 0x%016zX >= 0x%016zX", ec->srcfwdptr, ec->streamend);
#endif
		ec->state=ec->state<<32|*(const unsigned*)ec->srcfwdptr;
		ec->srcfwdptr+=4;
	}
	return sym;
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


//FSE		the following code was copied from https://github.com/Cyan4973/FiniteStateEntropy
#if 1
static const uint32_t SHIFT_mask[]=
{
	0x0,		0x1,		0x3,		0x7,
	0xF,		0x1F,		0x3F,		0x7F,
	0xFF,		0x1FF,		0x3FF,		0x7FF,
	0xFFF,		0x1FFF,		0x3FFF,		0x7FFF,
	0xFFFF,		0x1FFFF,	0x3FFFF,	0x7FFFF,
	0xFFFFF,	0x1FFFFF,	0x3FFFFF,	0x7FFFFF,
	0xFFFFFF,	0x1FFFFFF,	0x3FFFFFF,	0x7FFFFFF,
	0xFFFFFFF,	0x1FFFFFFF,	0x3FFFFFFF,	0x7FFFFFFF,
};
typedef struct _tANS_CState_t
{
	ptrdiff_t value;
	const void *stateTable, *symbolTT;
	unsigned stateLog;
} tANS_CState_t;
typedef struct _tANS_DState_t
{
	size_t state;
	const void *table;
} tANS_DState_t;
typedef struct _tANS_symbolCompressionTransform
{
	int32_t deltaFindState;
	uint32_t deltaNbBits;
} tANS_symbolCompressionTransform;
typedef unsigned tANS_CTable;
typedef unsigned tANS_DTable;
typedef struct _tANS_DTableHeader
{
	uint16_t tableLog, fastMode;
} tANS_DTableHeader;
typedef struct _tANS_decode_t
{
	uint16_t newState;
	uint8_t symbol, nbBits;
} tANS_decode_t;
#define tANS_MAX_SYMBOL_VALUE 255
#define tANS_FUNCTION_TYPE uint8_t
#define tANS_TABLELOG_ABSOLUTE_MAX 15
#define tANS_MAX_MEMORY_USAGE 14
#define tANS_DEFAULT_MEMORY_USAGE 13
#define tANS_MIN_TABLELOG 5
#define tANS_MAX_TABLELOG (tANS_MAX_MEMORY_USAGE-2)
#define tANS_DEFAULT_TABLELOG (tANS_DEFAULT_MEMORY_USAGE-2)
#define tANS_TABLESTEP(tableSize) (((tableSize)>>1) + ((tableSize)>>3) + 3)
#define tANS_CTABLE_SIZE_U32(maxTableLog, maxSymbolValue)	(1+(1<<((maxTableLog)-1))+(((maxSymbolValue)+1)*2))
#define tANS_DTABLE_SIZE_U32(maxTableLog)			(1+(1<<(maxTableLog)))
#define tANS_DECODE_TYPE tANS_decode_t

AWM_INLINE void tANS_initCState(tANS_CState_t *statePtr, const tANS_CTable *ct)
{
	const void *ptr=ct;
	const uint16_t *u16ptr=(const uint16_t*)ptr;
	const uint32_t tableLog=*(const uint16_t*)ptr;
	statePtr->value=(ptrdiff_t)1<<tableLog;
	statePtr->stateTable=u16ptr+2;
	statePtr->symbolTT=ct+1+(tableLog?(1<<(tableLog-1)):1);
	statePtr->stateLog=tableLog;
}
AWM_INLINE void tANS_encodeSymbol(BitPackerLIFO *ec, tANS_CState_t *statePtr, unsigned symbol)
{
	tANS_symbolCompressionTransform const symbolTT=((const tANS_symbolCompressionTransform*)(statePtr->symbolTT))[symbol];
	const uint16_t *const stateTable=(const uint16_t*)(statePtr->stateTable);
	uint32_t const nbBitsOut=(uint32_t)((statePtr->value+symbolTT.deltaNbBits)>>16);
	bitpacker_enc(ec, nbBitsOut, SHIFT_mask[nbBitsOut]&(int)statePtr->value);
	//BIT_addBits(bitC, statePtr->value, nbBitsOut);
	statePtr->value=stateTable[(statePtr->value>>nbBitsOut)+symbolTT.deltaFindState];
}
AWM_INLINE void tANS_flushCState(BitPackerLIFO *ec, const tANS_CState_t *statePtr)
{
	bitpacker_enc(ec, statePtr->stateLog, SHIFT_mask[statePtr->stateLog]&(int)statePtr->value);
	//BIT_addBits(bitC, statePtr->value, statePtr->stateLog);
	//bitpacker_enc_flush(ec);
	//BIT_flushBits(bitC);
}

AWM_INLINE void tANS_initDState(tANS_DState_t *DStatePtr, BitPackerLIFO *ec, const tANS_DTable *dt)
{
	const void *ptr=dt;
	const tANS_DTableHeader *const DTableH=(const tANS_DTableHeader*)ptr;
	DStatePtr->state=bitpacker_dec(ec, DTableH->tableLog);
	//DStatePtr->state=BIT_readBits(bitD, DTableH->tableLog);
	//BIT_reloadDStream(bitD);
	DStatePtr->table=dt+1;
}
AWM_INLINE BYTE tANS_decodeSymbol(tANS_DState_t *DStatePtr, BitPackerLIFO *ec)
{
	tANS_decode_t const DInfo=((const tANS_decode_t*)(DStatePtr->table))[DStatePtr->state];
	uint32_t const nbBits=DInfo.nbBits;
	BYTE const symbol=DInfo.symbol;
	size_t const lowBits=bitpacker_dec(ec, nbBits);
	//size_t const lowBits=BIT_readBits(bitD, nbBits);

	DStatePtr->state=DInfo.newState+lowBits;
	return symbol;
}

static void normalizehist(const uint32_t *hist, uint16_t *nhist)
{
	int hsum=0, nusedlevels=0;
	for(int ks=0;ks<256;++ks)//faster than maintaining hist sum
	{
		int freq=hist[ks];
		hsum+=freq;
		nusedlevels+=freq!=0;
	}
	long long rsum=(((1LL<<PROBBITS)-nusedlevels)<<24)/hsum;//adaptive: allow all symbols
	uint16_t CDF[257]={0};
	for(int ks=0, c=0, c2=0;ks<256;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*rsum>>24)+c2;
		c+=freq;
		c2+=freq!=0;
	}
	CDF[256]=1<<PROBBITS;
	for(int ks=0;ks<256;++ks)
		nhist[ks]=CDF[ks+1]-CDF[ks];
}
static tANS_CTable *tANS_createCTable(unsigned maxSymbolValue, unsigned tableLog)
{
	size_t size;
	if(tableLog>tANS_TABLELOG_ABSOLUTE_MAX)
		tableLog=tANS_TABLELOG_ABSOLUTE_MAX;
	size=tANS_CTABLE_SIZE_U32(tableLog, maxSymbolValue)*sizeof(uint32_t);
	return (tANS_CTable*)malloc(size);
}
static tANS_DTable *tANS_createDTable(unsigned tableLog)
{
	if(tableLog>tANS_TABLELOG_ABSOLUTE_MAX)
		tableLog=tANS_TABLELOG_ABSOLUTE_MAX;
	return (tANS_DTable*)malloc(tANS_DTABLE_SIZE_U32(tableLog)*sizeof(uint32_t));
}
//tANS_buildCTable_wksp() :
//Same as tANS_buildCTable(), but using an externally allocated scratch buffer (`workSpace`).
//wkspSize should be sized to handle worst case situation, which is `1<<max_tableLog*sizeof(tANS_FUNCTION_TYPE)`
//workSpace must also be properly aligned with tANS_FUNCTION_TYPE requirements
static size_t tANS_buildCTable_wksp(tANS_CTable *ct, const short *normalizedCounter, unsigned maxSymbolValue, unsigned tableLog, void *workSpace, size_t wkspSize)
{
	uint32_t const tableSize=1<<tableLog;
	uint32_t const tableMask=tableSize-1;
	void *const ptr=ct;
	uint16_t *const tableU16=(uint16_t*)ptr+2;
	void *const FSCT = ((uint32_t*)ptr)+1/*header*/+(tableLog ? tableSize>>1:1);
	tANS_symbolCompressionTransform *const symbolTT=(tANS_symbolCompressionTransform*)(FSCT);
	uint32_t const step=tANS_TABLESTEP(tableSize);
	uint32_t cumul[tANS_MAX_SYMBOL_VALUE+2];

	tANS_FUNCTION_TYPE *const tableSymbol=(tANS_FUNCTION_TYPE*)workSpace;
	uint32_t highThreshold=tableSize-1;

	//CTable header
	if(((size_t)1<<tableLog)*sizeof(tANS_FUNCTION_TYPE)>wkspSize)
	{
		CRASH("tableLog too large");
		return 1;
	}
	tableU16[-2]=(uint16_t)tableLog;
	tableU16[-1]=(uint16_t)maxSymbolValue;
	if(tableLog>=16)//required for threshold strategy to work
	{
		CRASH("tableLog too large");
		return 1;
	}

	//For explanations on how to distribute symbol values over the table
	//http://fastcompression.blogspot.fr/2014/02/fse-distributing-symbol-values.html

#ifdef __clang_analyzer__
	memset(tableSymbol, 0, sizeof(*tableSymbol)*tableSize);//useless initialization, just to keep scan-build happy
#endif

	//symbol start positions
	{
		uint32_t u;
		cumul[0]=0;
		for(u=1;u<=maxSymbolValue+1;++u)
		{
			if(normalizedCounter[u-1]==-1)//Low probability symbol
			{
				cumul[u]=cumul[u-1]+1;
				tableSymbol[highThreshold--]=(tANS_FUNCTION_TYPE)(u-1);
			}
			else
				cumul[u]=cumul[u-1]+normalizedCounter[u-1];
		}
		cumul[maxSymbolValue+1]=tableSize+1;
	}

	//Spread symbols
	{
		uint32_t position=0, symbol;
		for(symbol=0;symbol<=maxSymbolValue;++symbol)
		{
			int freq=normalizedCounter[symbol], nbOccurrences;
			for(nbOccurrences=0;nbOccurrences<freq;++nbOccurrences)
			{
				tableSymbol[position]=(tANS_FUNCTION_TYPE)symbol;
				position=(position+step)&tableMask;
				while(position>highThreshold)
					position=(position+step)&tableMask;//Low probability area
			}
		}
		if(position)//Must have initialized all positions
		{
			CRASH("Must have initialized all positions");
			return 1;
		}
	}

	//Build table
	{
		uint32_t u;
		for(u=0;u<tableSize;++u)
		{
			tANS_FUNCTION_TYPE s=tableSymbol[u];//note : static analyzer may not understand tableSymbol is properly initialized
			tableU16[cumul[s]++]=(uint16_t)(tableSize+u);//TableU16 : sorted by symbol order; gives next state value
		}
	}

	//Build Symbol Transformation Table
	{
		unsigned total=0, s;
		for(s=0;s<=maxSymbolValue;++s)
		{
			switch(normalizedCounter[s])
			{
			case  0://filling nonetheless, for compatibility with tANS_getMaxNbBits()
				symbolTT[s].deltaNbBits=((tableLog+1)<<16)-(1<<tableLog);
//#ifdef _DEBUG
//				symbolTT[s].deltaFindState=0xCDCDCDCD;//AWM-20251122
//#endif
				break;

			case -1:
			case  1:
				symbolTT[s].deltaNbBits=(tableLog<<16)-(1<<tableLog);
				symbolTT[s].deltaFindState=total-1;
				++total;
				break;
			default:
				{
					uint32_t const maxBitsOut=tableLog-(31-_lzcnt_u32(normalizedCounter[s]-1));
					uint32_t const minStatePlus=normalizedCounter[s]<<maxBitsOut;
					symbolTT[s].deltaNbBits=(maxBitsOut<<16)-minStatePlus;
					symbolTT[s].deltaFindState=total-normalizedCounter[s];
					total+=normalizedCounter[s];
				}
				break;
			}
		}
	}

#if defined _DEBUG && 0		//debug : symbol costs
	printf("table statistics\n");
	{
		uint32_t symbol;
		for(symbol=0;symbol<=maxSymbolValue;++symbol)
		{
			printf("%3u: w=%3i,   maxBits=%u, fracBits=%.2f"
				, symbol
				, normalizedCounter[symbol]
				, tANS_getMaxNbBits(symbolTT, symbol)
				, (double)tANS_bitCost(symbolTT, tableLog, symbol, 8)/256
			);
		}
	}
#endif
	return 0;
}
static size_t tANS_buildCTable(tANS_CTable *ct, const short *normalizedCounter, unsigned maxSymbolValue, unsigned tableLog)
{
	uint8_t tableSymbol[1<<14>>2];//memset() is not necessary, even if static analyzer complains about it
	return tANS_buildCTable_wksp(ct, normalizedCounter, maxSymbolValue, tableLog, tableSymbol, sizeof(tableSymbol));
}
static size_t tANS_buildDTable(tANS_DTable *dt, const short *normalizedCounter, unsigned maxSymbolValue, unsigned tableLog)
{
	void *const tdPtr=dt+1;//because *dt is unsigned, 32-bits aligned on 32-bits
	tANS_DECODE_TYPE *const tableDecode=(tANS_DECODE_TYPE*)tdPtr;
	uint16_t symbolNext[tANS_MAX_SYMBOL_VALUE+1];

	uint32_t const maxSV1 = maxSymbolValue+1;
	uint32_t const tableSize=1<<tableLog;
	uint32_t highThreshold=tableSize-1;

	//Sanity Checks
	if(maxSymbolValue>tANS_MAX_SYMBOL_VALUE)
	{
		CRASH("maxSymbolValue too large");
		return 1;
	}
	if(tableLog>tANS_MAX_TABLELOG)
	{
		CRASH("tableLog too large");
		return 1;
	}

	//Init, lay down lowprob symbols
	{
		tANS_DTableHeader DTableH;
		DTableH.tableLog=(uint16_t)tableLog;
		DTableH.fastMode=1;
		{
			int16_t const largeLimit=(int16_t)(1<<(tableLog-1));
			uint32_t s;
			for(s=0;s<maxSV1;++s)
			{
				if(normalizedCounter[s]==-1)
				{
					tableDecode[highThreshold--].symbol=(tANS_FUNCTION_TYPE)s;
					symbolNext[s]=1;
				}
				else
				{
					if(normalizedCounter[s]>=largeLimit)
						DTableH.fastMode=0;
					symbolNext[s]=normalizedCounter[s];
				}
			}
		}
		memcpy(dt, &DTableH, sizeof(DTableH));
	}

	//Spread symbols
	{
		const uint32_t tableMask=tableSize-1, step=tANS_TABLESTEP(tableSize);
		uint32_t s, position=0;
		for(s=0;s<maxSV1;++s)
		{
			int i;
			for(i=0;i<normalizedCounter[s];++i)
			{
				tableDecode[position].symbol=(tANS_FUNCTION_TYPE)s;
				position=(position+step)&tableMask;
				while(position>highThreshold)
					position=(position+step)&tableMask;//lowprob area
			}
		}
		if(position!=0)//position must reach all cells once, otherwise normalizedCounter is incorrect
		{
			CRASH("position must reach all cells once, otherwise normalizedCounter is incorrect");
			return 1;
		}
	}

	//Build Decoding table
	{
		uint32_t u;
		for(u=0;u<tableSize;++u)
		{
			tANS_FUNCTION_TYPE const symbol=(tANS_FUNCTION_TYPE)tableDecode[u].symbol;
			uint32_t const nextState=symbolNext[symbol]++;
			tableDecode[u].nbBits=(BYTE)(tableLog-(31-_lzcnt_u32(nextState)));
			tableDecode[u].newState=(uint16_t)((nextState<<tableDecode[u].nbBits)-tableSize);
		}
	}
	return 0;
}
static size_t tANS_buildCTable_rle(tANS_CTable *ct, BYTE symbolValue)
{
	void *ptr=ct;
	uint16_t *tableU16=(uint16_t*)ptr+2;
	void *FSCTptr=(uint32_t*)ptr+2;
	tANS_symbolCompressionTransform *symbolTT=(tANS_symbolCompressionTransform*)FSCTptr;

	//header
	tableU16[-2]=(uint16_t)0;
	tableU16[-1]=(uint16_t)symbolValue;

	//Build table
	tableU16[0]=0;
	tableU16[1]=0;//just in case

	//Build Symbol Transformation Table
	symbolTT[symbolValue].deltaNbBits=0;
	symbolTT[symbolValue].deltaFindState=0;
	return 0;
}
static size_t tANS_buildDTable_rle(tANS_DTable *dt, BYTE symbolValue)
{
	void *ptr=dt;
	tANS_DTableHeader *const DTableH=(tANS_DTableHeader*)ptr;
	void *dPtr=dt+1;
	tANS_decode_t *const cell=(tANS_decode_t*)dPtr;

	DTableH->tableLog=0;
	DTableH->fastMode=0;

	cell->newState=0;
	cell->symbol=symbolValue;
	cell->nbBits=0;
	return 0;
}
static size_t tANS_buildCTable_raw(tANS_CTable *ct, unsigned nbBits)
{
	const unsigned tableSize=1<<nbBits;
	const unsigned tableMask=tableSize-1;
	const unsigned maxSymbolValue=tableMask;
	void *const ptr=ct;
	uint16_t *const tableU16=(uint16_t*)ptr+2;
	void *const FSCT=(uint32_t*)ptr+1/*header*/+(tableSize>>1);//assumption: tableLog >= 1
	tANS_symbolCompressionTransform *const symbolTT = (tANS_symbolCompressionTransform*) (FSCT);
	unsigned s;

	//Sanity checks
	if(nbBits<1)
	{
		CRASH("nbBits");//min size
		return 1;
	}

	//header
	tableU16[-2]=(uint16_t)nbBits;
	tableU16[-1]=(uint16_t)maxSymbolValue;

	//Build table
	for(s=0;s<tableSize;++s)
		tableU16[s]=(uint16_t)(tableSize+s);

	//Build Symbol Transformation Table
	{
		const uint32_t deltaNbBits=(nbBits<<16)-(1<<nbBits);
		for(s=0;s<=maxSymbolValue;++s)
		{
			symbolTT[s].deltaNbBits=deltaNbBits;
			symbolTT[s].deltaFindState=s-1;
		}
	}
	return 0;
}
static size_t tANS_buildDTable_raw(tANS_DTable *dt, unsigned nbBits)
{
	void *ptr=dt;
	tANS_DTableHeader *const DTableH=(tANS_DTableHeader*)ptr;
	void *dPtr=dt+1;
	tANS_decode_t *const dinfo=(tANS_decode_t*)dPtr;
	const unsigned tableSize=1<<nbBits;
	const unsigned tableMask=tableSize-1;
	const unsigned maxSV1=tableMask+1;
	unsigned s;

	//Sanity checks
	if(nbBits<1)//min size
	{
		CRASH("nbBits");
		return 1;
	}

	//Build Decoding Table
	DTableH->tableLog=(uint16_t)nbBits;
	DTableH->fastMode=1;
	for(s=0;s<maxSV1;++s)
	{
		dinfo[s].newState=0;
		dinfo[s].symbol=(BYTE)s;
		dinfo[s].nbBits=(BYTE)nbBits;
	}
	return 0;
}
#endif
//end of copied FSE code

//rANS	https://github.com/rygorous/ryg_rans	https://github.com/samtools/htscodecs
typedef struct _rANS_SIMD_SymInfo	//16 bytes/level	4KB/ctx = 1<<12 bytes
{
	uint32_t smax, invf, cdf;
	uint16_t negf, sh;
} rANS_SIMD_SymInfo;
static void enc_hist2stats(int *hist, rANS_SIMD_SymInfo *syminfo, unsigned long long *bypassmask, int ctxidx)
{
	int sum=0, count=0;
	for(int ks=0;ks<256;++ks)
	{
		int freq=hist[ks];
		sum+=freq;
		count+=freq!=0;
	}
	int rare=sum<12*256/8;
	*bypassmask|=(unsigned long long)rare<<ctxidx;
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

			printf("%c", "0123456789ABCDEF"[ks&15]);
			//int shade=48+freq*(255-48)/fmax;
			//colorprintf(shade<<16|shade<<8|shade, freq?0x808080:COLORPRINTF_BK_DEFAULT, "%c", "0123456789ABCDEF"[ks&15]);
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
			info->sh=31-(int)_lzcnt_u32(freq);//eg: x/2 = x*0x80000000>>32>>0
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
static void enc_packhist(BitPackerLIFO *ec, const int *hist, unsigned long long bypassmask, int ctxidx, int userans)//histogram must be normalized to PROBBITS, with spike at 128
{
	if(bypassmask>>ctxidx&1)
		return;
	int sum=0;
	uint16_t hist2[256];
	if(userans)//32-bit hist
	{
		for(int ks=0;ks<256;++ks)
			hist2[ks]=hist[ks];
	}
	else//16-bit hist
		memcpy(hist2, hist, sizeof(hist2));
	unsigned short CDF[257];
	for(int ks=0;ks<256;++ks)//integrage to zigzag CDF to be packed backwards
	{
		int sym=((ks>>1^-(ks&1))+128)&255;
		int freq=hist2[sym];
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
#ifdef ANS_VAL
		ansval_push(&freq, sizeof(freq), 1);
#endif
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
static void dec_unpackhist(BitPackerLIFO *ec, unsigned long long bypassmask, int ctxidx, uint16_t *hist)
{
	if(bypassmask>>ctxidx&1)//rare context
	{
		for(int ks=0;ks<256;++ks)//bypass
			hist[ks]=(1<<PROBBITS)/256;
	}
	else
	{
		unsigned short CDF[257]={0};
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
#ifdef ANS_VAL
			ansval_check(&freq, sizeof(freq), 1);
#endif

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
	//int sum=0;
	//for(int ks=0;ks<256;++ks)//integrate
	//{
	//	int freq=hist[ks];
	//	hist[ks]=sum;
	//	sum+=freq;
	//}
	//hist[256]=1<<PROBBITS;
	//for(int ks=0;ks<256;++ks)//CDF2sym contains {freq, (state&0xFFF)-cdf, sym}
	//{
	//	int cdf=hist[ks], next=hist[ks+1], freq=next-cdf;
	//	int val=(freq<<PROBBITS|0)<<8|ks;
	//	for(int ks2=cdf;ks2<next;++ks2, val+=1<<8)
	//		CDF2sym[ks2]=val;
	//}
}
static void dec_hist2CDF2sym(const uint16_t *hist, unsigned *CDF2sym)
{
	int cdf=0;
	for(int ks=0;ks<256;++ks)//CDF2sym contains {freq, (state&0xFFF)-cdf, sym}
	{
		int freq=hist[ks], next=cdf+freq;
		int val=(freq<<PROBBITS|0)<<8|ks;
		for(int ks2=cdf;ks2<next;++ks2, val+=1<<8)
			CDF2sym[ks2]=val;
		cdf=next;
	}
}

AWM_INLINE void dec_yuv(__m256i *mstate, int kc, const __m256i *ctx0, const int *CDF2syms, const int *ans_permute, unsigned char **pstreamptr, const unsigned char *streamend, __m256i *syms)
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
	ansval_check(mstate, sizeof(int), NCODERS);
#endif
	{
		__m256i mprobmask=_mm256_set1_epi32((1<<PROBBITS)-1);
		__m256i rem0=_mm256_and_si256(mstate[0], mprobmask);
		__m256i rem1=_mm256_and_si256(mstate[1], mprobmask);
		decctx[0]=_mm256_or_si256(decctx[0], rem0);
		decctx[1]=_mm256_or_si256(decctx[1], rem1);
	}
#ifdef ANS_VAL
	ALIGN(32) int debugctx[NCODERS];
	memcpy(debugctx, decctx, sizeof(int[NCODERS]));
#endif
	const int *statsptr=CDF2syms+((ptrdiff_t)NCTX*kc<<PROBBITS);
	gather32((int*)(decctx+0), statsptr, (int*)(decctx+0));
	gather32((int*)(decctx+1), statsptr, (int*)(decctx+1));
	//decctx[0]=_mm256_i32gather_epi32(statsptr, decctx[0], sizeof(int));
	//decctx[1]=_mm256_i32gather_epi32(statsptr, decctx[1], sizeof(int));

	//update		state = (state>>12)*freq+(rem-cdf)	rem-cdf is prebaked
	{
		__m256i mfreq0=_mm256_srli_epi32(decctx[0], PROBBITS+8);//1 <= freq <= 0xF01
		__m256i mfreq1=_mm256_srli_epi32(decctx[1], PROBBITS+8);
#ifdef ANS_VAL
		__m256i mdebugfreq[1];
		mdebugfreq[0]=_mm256_packus_epi32(mfreq0, mfreq1);
		mdebugfreq[0]=_mm256_permute4x64_epi64(mdebugfreq[0], _MM_SHUFFLE(3, 1, 2, 0));
		ALIGN(32) unsigned short freqs[NCODERS];
		memcpy(freqs, mdebugfreq, sizeof(freqs));
		ansval_check(freqs, sizeof(short), NCODERS);
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
	ansval_check(mstate, sizeof(int), NCODERS);
#endif
	//renorm
	{
		__m256i smin=_mm256_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
#ifdef _DEBUG
		if(streamptr>streamend)
			CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__m256i lo0=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr));
		__m256i cond0=_mm256_cmpgt_epi32(smin, mstate[0]);//FIXME this is signed comparison
		__m256i cond1=_mm256_cmpgt_epi32(smin, mstate[1]);
		int mask0=_mm256_movemask_ps(_mm256_castsi256_ps(cond0));
		int mask1=_mm256_movemask_ps(_mm256_castsi256_ps(cond1));
		__m256i idx0=_mm256_load_si256((const __m256i*)ans_permute+mask0);
		__m256i idx1=_mm256_load_si256((const __m256i*)ans_permute+mask1);
		mask0=_mm_popcnt_u32(mask0);
		mask1=_mm_popcnt_u32(mask1);
		__m256i renorm0=_mm256_slli_epi32(mstate[0], 16);
		__m256i renorm1=_mm256_slli_epi32(mstate[1], 16);
		streamptr+=mask0*sizeof(short);
		lo0=_mm256_permutevar8x32_epi32(lo0, idx0);
#ifdef _DEBUG
		if(streamptr>streamend)
			CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__m256i lo1=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr));
		streamptr+=mask1*sizeof(short);
		lo1=_mm256_permutevar8x32_epi32(lo1, idx1);
		renorm0=_mm256_or_si256(renorm0, lo0);
		renorm1=_mm256_or_si256(renorm1, lo1);

		mstate[0]=_mm256_blendv_epi8(mstate[0], renorm0, cond0);
		mstate[1]=_mm256_blendv_epi8(mstate[1], renorm1, cond1);
	}
	*pstreamptr=(unsigned char*)(size_t)streamptr;
}
AWM_INLINE void transpose16(__m128i *data)
{
#if 1
	__m128i a[16], b[16];
	a[0x0]=_mm_load_si128((__m128i*)data+0x0);
	a[0x1]=_mm_load_si128((__m128i*)data+0x1);
	a[0x2]=_mm_load_si128((__m128i*)data+0x2);
	a[0x3]=_mm_load_si128((__m128i*)data+0x3);
	a[0x4]=_mm_load_si128((__m128i*)data+0x4);
	a[0x5]=_mm_load_si128((__m128i*)data+0x5);
	a[0x6]=_mm_load_si128((__m128i*)data+0x6);
	a[0x7]=_mm_load_si128((__m128i*)data+0x7);
	a[0x8]=_mm_load_si128((__m128i*)data+0x8);
	a[0x9]=_mm_load_si128((__m128i*)data+0x9);
	a[0xA]=_mm_load_si128((__m128i*)data+0xA);
	a[0xB]=_mm_load_si128((__m128i*)data+0xB);
	a[0xC]=_mm_load_si128((__m128i*)data+0xC);
	a[0xD]=_mm_load_si128((__m128i*)data+0xD);
	a[0xE]=_mm_load_si128((__m128i*)data+0xE);
	a[0xF]=_mm_load_si128((__m128i*)data+0xF);

	b[0x0]=_mm_unpacklo_epi8(a[0x0], a[0x1]);
	b[0x1]=_mm_unpackhi_epi8(a[0x0], a[0x1]);
	b[0x2]=_mm_unpacklo_epi8(a[0x2], a[0x3]);
	b[0x3]=_mm_unpackhi_epi8(a[0x2], a[0x3]);
	b[0x4]=_mm_unpacklo_epi8(a[0x4], a[0x5]);
	b[0x5]=_mm_unpackhi_epi8(a[0x4], a[0x5]);
	b[0x6]=_mm_unpacklo_epi8(a[0x6], a[0x7]);
	b[0x7]=_mm_unpackhi_epi8(a[0x6], a[0x7]);
	b[0x8]=_mm_unpacklo_epi8(a[0x8], a[0x9]);
	b[0x9]=_mm_unpackhi_epi8(a[0x8], a[0x9]);
	b[0xA]=_mm_unpacklo_epi8(a[0xA], a[0xB]);
	b[0xB]=_mm_unpackhi_epi8(a[0xA], a[0xB]);
	b[0xC]=_mm_unpacklo_epi8(a[0xC], a[0xD]);
	b[0xD]=_mm_unpackhi_epi8(a[0xC], a[0xD]);
	b[0xE]=_mm_unpacklo_epi8(a[0xE], a[0xF]);
	b[0xF]=_mm_unpackhi_epi8(a[0xE], a[0xF]);

	a[0x0]=_mm_unpacklo_epi16(b[0x0], b[0x2]);
	a[0x1]=_mm_unpackhi_epi16(b[0x0], b[0x2]);
	a[0x2]=_mm_unpacklo_epi16(b[0x1], b[0x3]);
	a[0x3]=_mm_unpackhi_epi16(b[0x1], b[0x3]);
	a[0x4]=_mm_unpacklo_epi16(b[0x4], b[0x6]);
	a[0x5]=_mm_unpackhi_epi16(b[0x4], b[0x6]);
	a[0x6]=_mm_unpacklo_epi16(b[0x5], b[0x7]);
	a[0x7]=_mm_unpackhi_epi16(b[0x5], b[0x7]);
	a[0x8]=_mm_unpacklo_epi16(b[0x8], b[0xA]);
	a[0x9]=_mm_unpackhi_epi16(b[0x8], b[0xA]);
	a[0xA]=_mm_unpacklo_epi16(b[0x9], b[0xB]);
	a[0xB]=_mm_unpackhi_epi16(b[0x9], b[0xB]);
	a[0xC]=_mm_unpacklo_epi16(b[0xC], b[0xE]);
	a[0xD]=_mm_unpackhi_epi16(b[0xC], b[0xE]);
	a[0xE]=_mm_unpacklo_epi16(b[0xD], b[0xF]);
	a[0xF]=_mm_unpackhi_epi16(b[0xD], b[0xF]);

	b[0x0]=_mm_unpacklo_epi32(a[0x0], a[0x4]);
	b[0x1]=_mm_unpackhi_epi32(a[0x0], a[0x4]);
	b[0x2]=_mm_unpacklo_epi32(a[0x1], a[0x5]);
	b[0x3]=_mm_unpackhi_epi32(a[0x1], a[0x5]);
	b[0x4]=_mm_unpacklo_epi32(a[0x2], a[0x6]);
	b[0x5]=_mm_unpackhi_epi32(a[0x2], a[0x6]);
	b[0x6]=_mm_unpacklo_epi32(a[0x3], a[0x7]);
	b[0x7]=_mm_unpackhi_epi32(a[0x3], a[0x7]);
	b[0x8]=_mm_unpacklo_epi32(a[0x8], a[0xC]);
	b[0x9]=_mm_unpackhi_epi32(a[0x8], a[0xC]);
	b[0xA]=_mm_unpacklo_epi32(a[0x9], a[0xD]);
	b[0xB]=_mm_unpackhi_epi32(a[0x9], a[0xD]);
	b[0xC]=_mm_unpacklo_epi32(a[0xA], a[0xE]);
	b[0xD]=_mm_unpackhi_epi32(a[0xA], a[0xE]);
	b[0xE]=_mm_unpacklo_epi32(a[0xB], a[0xF]);
	b[0xF]=_mm_unpackhi_epi32(a[0xB], a[0xF]);

	a[0x0]=_mm_unpacklo_epi64(b[0x0], b[0x8]);
	a[0x1]=_mm_unpackhi_epi64(b[0x0], b[0x8]);
	a[0x2]=_mm_unpacklo_epi64(b[0x1], b[0x9]);
	a[0x3]=_mm_unpackhi_epi64(b[0x1], b[0x9]);
	a[0x4]=_mm_unpacklo_epi64(b[0x2], b[0xA]);
	a[0x5]=_mm_unpackhi_epi64(b[0x2], b[0xA]);
	a[0x6]=_mm_unpacklo_epi64(b[0x3], b[0xB]);
	a[0x7]=_mm_unpackhi_epi64(b[0x3], b[0xB]);
	a[0x8]=_mm_unpacklo_epi64(b[0x4], b[0xC]);
	a[0x9]=_mm_unpackhi_epi64(b[0x4], b[0xC]);
	a[0xA]=_mm_unpacklo_epi64(b[0x5], b[0xD]);
	a[0xB]=_mm_unpackhi_epi64(b[0x5], b[0xD]);
	a[0xC]=_mm_unpacklo_epi64(b[0x6], b[0xE]);
	a[0xD]=_mm_unpackhi_epi64(b[0x6], b[0xE]);
	a[0xE]=_mm_unpacklo_epi64(b[0x7], b[0xF]);
	a[0xF]=_mm_unpackhi_epi64(b[0x7], b[0xF]);

	_mm_store_si128((__m128i*)data+0x0, a[0x0]);
	_mm_store_si128((__m128i*)data+0x1, a[0x1]);
	_mm_store_si128((__m128i*)data+0x2, a[0x2]);
	_mm_store_si128((__m128i*)data+0x3, a[0x3]);
	_mm_store_si128((__m128i*)data+0x4, a[0x4]);
	_mm_store_si128((__m128i*)data+0x5, a[0x5]);
	_mm_store_si128((__m128i*)data+0x6, a[0x6]);
	_mm_store_si128((__m128i*)data+0x7, a[0x7]);
	_mm_store_si128((__m128i*)data+0x8, a[0x8]);
	_mm_store_si128((__m128i*)data+0x9, a[0x9]);
	_mm_store_si128((__m128i*)data+0xA, a[0xA]);
	_mm_store_si128((__m128i*)data+0xB, a[0xB]);
	_mm_store_si128((__m128i*)data+0xC, a[0xC]);
	_mm_store_si128((__m128i*)data+0xD, a[0xD]);
	_mm_store_si128((__m128i*)data+0xE, a[0xE]);
	_mm_store_si128((__m128i*)data+0xF, a[0xF]);
#endif

#if 0
	//https://pzemtsov.github.io/2014/10/01/how-to-transpose-a-16x16-matrix.html
	__m128i t0, t1, t2, t3;
#define SHUFFLE(LO, HI, IMM8_HHLL) _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(LO), _mm_castsi128_ps(HI), IMM8_HHLL))
#define TRANSPOSE4(S0, S1, S2, S3,  T0, T1, T2, T3,  D0, D1, D2, D3)\
	do\
	{\
		T0=SHUFFLE(S0, S1, _MM_SHUFFLE(1, 0, 1, 0));\
		T2=SHUFFLE(S0, S1, _MM_SHUFFLE(3, 2, 3, 2));\
		T1=SHUFFLE(S2, S3, _MM_SHUFFLE(1, 0, 1, 0));\
		T3=SHUFFLE(S2, S3, _MM_SHUFFLE(3, 2, 3, 2));\
		D0=SHUFFLE(T0, T1, _MM_SHUFFLE(2, 0, 2, 0));\
		D1=SHUFFLE(T0, T1, _MM_SHUFFLE(3, 1, 3, 1));\
		D2=SHUFFLE(T2, T3, _MM_SHUFFLE(2, 0, 2, 0));\
		D3=SHUFFLE(T2, T3, _MM_SHUFFLE(3, 1, 3, 1));\
	}while(0)
	TRANSPOSE4(data[0x0], data[0x1], data[0x2], data[0x3],  t0, t1, t2, t3,  data[0x0], data[0x1], data[0x2], data[0x3]);
	TRANSPOSE4(data[0x4], data[0x5], data[0x6], data[0x7],  t0, t1, t2, t3,  data[0x4], data[0x5], data[0x6], data[0x7]);
	TRANSPOSE4(data[0x8], data[0x9], data[0xA], data[0xB],  t0, t1, t2, t3,  data[0x8], data[0x9], data[0xA], data[0xB]);
	TRANSPOSE4(data[0xC], data[0xD], data[0xE], data[0xF],  t0, t1, t2, t3,  data[0xC], data[0xD], data[0xE], data[0xF]);
	t0=_mm_set_epi8(
		15, 11,  7,  3,
		14, 10,  6,  2,
		13,  9,  5,  1,
		12,  8,  4,  0
	);
	/*
	transpose the 4*4 registers themselves while shuffling
	exchange:
	1 <-> 4
	2 <-> 8
	3 <-> C
	6 <-> 9
	7 <-> D
	B <-> E
	*/
	//diagonals
	data[0x0]=_mm_shuffle_epi8(data[0x0], t0);
	data[0x5]=_mm_shuffle_epi8(data[0x5], t0);
	data[0xA]=_mm_shuffle_epi8(data[0xA], t0);
	data[0xF]=_mm_shuffle_epi8(data[0xF], t0);
	//nondiagonals
	t1=data[0x4]; data[0x4]=_mm_shuffle_epi8(data[0x1], t0); data[0x1]=_mm_shuffle_epi8(t1, t0);
	t1=data[0x8]; data[0x8]=_mm_shuffle_epi8(data[0x2], t0); data[0x2]=_mm_shuffle_epi8(t1, t0);
	t1=data[0xC]; data[0xC]=_mm_shuffle_epi8(data[0x3], t0); data[0x3]=_mm_shuffle_epi8(t1, t0);
	t1=data[0x9]; data[0x9]=_mm_shuffle_epi8(data[0x6], t0); data[0x6]=_mm_shuffle_epi8(t1, t0);
	t1=data[0xD]; data[0xD]=_mm_shuffle_epi8(data[0x7], t0); data[0x7]=_mm_shuffle_epi8(t1, t0);
	t1=data[0xE]; data[0xE]=_mm_shuffle_epi8(data[0xB], t0); data[0xB]=_mm_shuffle_epi8(t1, t0);
	TRANSPOSE4(data[0x0], data[0x1], data[0x2], data[0x3],  t0, t1, t2, t3,  data[0x0], data[0x1], data[0x2], data[0x3]);
	TRANSPOSE4(data[0x4], data[0x5], data[0x6], data[0x7],  t0, t1, t2, t3,  data[0x4], data[0x5], data[0x6], data[0x7]);
	TRANSPOSE4(data[0x8], data[0x9], data[0xA], data[0xB],  t0, t1, t2, t3,  data[0x8], data[0x9], data[0xA], data[0xB]);
	TRANSPOSE4(data[0xC], data[0xD], data[0xE], data[0xF],  t0, t1, t2, t3,  data[0xC], data[0xD], data[0xE], data[0xF]);
#undef  SHUFFLE
#undef  TRANSPOSE4
#endif
}
static void interleave_blocks_fwd(const unsigned char *original, int iw, int ih, unsigned char *interleaved)
{
	//original[ih][iw][3]
	//interleaved[ih/YCODERS][iw/XCODERS][3][NCODERS]	contiguous & aligned
	//xrem[ih%YCODERS][iw][3]
	//yrem[ih-iw%YCODERS][iw%XCODERS][3]

	//only difference between fwd and inv:		swap assignments (const slow->fast)
	int rowstride=3*iw;
	int ixyblockw=iw/XCODERS;
	int ixyblockh=ih/YCODERS;
	int blockxbytes=3*NCODERS*ixyblockw;
#ifdef INTERLEAVESIMD
	int SIMDxcount=blockxbytes&~((int)sizeof(__m128i[NCODERS])-1);
	__m256i slowinc=_mm256_set1_epi64x(sizeof(__m128i));
#endif
	unsigned char *fastptr=interleaved;
	ALIGN(32) const unsigned char *slowptrs[NCODERS]={0}, *slowptrs0[NCODERS]={0};
	for(int ky=0;ky<YCODERS;++ky)//spread slow pointers
	{
		for(int kx=0;kx<XCODERS;++kx)
			slowptrs0[XCODERS*ky+kx]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);
	}
	for(int ky=0;ky<ixyblockh;++ky)//interleave
	{
		int kx=0;
		memcpy((void*)slowptrs, slowptrs0, sizeof(slowptrs));
#ifdef INTERLEAVESIMD
		for(;kx<SIMDxcount;kx+=(int)sizeof(__m128i[NCODERS]))
		{
			/*
			SIMD interleave
			|000111222333444555666777888999AAABBBCCCDDDEEEFFFGGGHHHIIIJJJKKKLLLMMMNNN	24 pixels * 3 channels
			|^        ^        ^        ^        ^        ^        ^        ^		8-lane serial interleave
			|^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^		8-lane 4xSIMD interleave step 1
			|    ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^	8-lane 4xSIMD interleave step 2
			|        ^        ^        ^        ^        ^        ^        ^        ^	8-lane serial remainder
			|0369CFIL 0369CFIL 0369CFIL 147ADGJM 147ADGJM 147ADGJM 258BEHKN 258BEHKN 258BEHKN	serial behavior
			|0369CFIL0369CFIL0369CFIL147ADGJM  147ADGJM147ADGJM258BEHKN258BEHKN  2 5 8 B E H K N	SIMD quotients & serial remainder

			slowptrs[NCODERS]		fastptr  (aligned because NCODERS == sizeof(__m256i))
			A16x16 B16x16		->	A16x16T C16x16T B16x16T D16x16T
			C16x16 D16x16

			load, transpose, store blocks A, B
			load, transpose, store blocks C, D

			speed:
			SIMD maximum estimate	(3/load + 3/store)*32+80*0.5*2 = 272 cycles for 32*32 bytes  ->  3GHz*32*32/272 ~= 10770 MB/s without remainders, cache misses
			actual SIMD		MSVC 4900 MB/s		GCC 3900/3000 MB/s
			actual serial		1900 MB/s
			*/
			__m128i block[16];
			block[0x0]=_mm_loadu_si128((__m128i*)slowptrs[0x00]);
			block[0x1]=_mm_loadu_si128((__m128i*)slowptrs[0x01]);
			block[0x2]=_mm_loadu_si128((__m128i*)slowptrs[0x02]);
			block[0x3]=_mm_loadu_si128((__m128i*)slowptrs[0x03]);
			block[0x4]=_mm_loadu_si128((__m128i*)slowptrs[0x04]);
			block[0x5]=_mm_loadu_si128((__m128i*)slowptrs[0x05]);
			block[0x6]=_mm_loadu_si128((__m128i*)slowptrs[0x06]);
			block[0x7]=_mm_loadu_si128((__m128i*)slowptrs[0x07]);
			block[0x8]=_mm_loadu_si128((__m128i*)slowptrs[0x08]);
			block[0x9]=_mm_loadu_si128((__m128i*)slowptrs[0x09]);
			block[0xA]=_mm_loadu_si128((__m128i*)slowptrs[0x0A]);
			block[0xB]=_mm_loadu_si128((__m128i*)slowptrs[0x0B]);
			block[0xC]=_mm_loadu_si128((__m128i*)slowptrs[0x0C]);
			block[0xD]=_mm_loadu_si128((__m128i*)slowptrs[0x0D]);
			block[0xE]=_mm_loadu_si128((__m128i*)slowptrs[0x0E]);
			block[0xF]=_mm_loadu_si128((__m128i*)slowptrs[0x0F]);
			transpose16(block);
			_mm_store_si128((__m128i*)fastptr+0x0, block[0x0]);
			_mm_store_si128((__m128i*)fastptr+0x1, block[0x1]);
			_mm_store_si128((__m128i*)fastptr+0x2, block[0x2]);
			_mm_store_si128((__m128i*)fastptr+0x3, block[0x3]);
			_mm_store_si128((__m128i*)fastptr+0x4, block[0x4]);
			_mm_store_si128((__m128i*)fastptr+0x5, block[0x5]);
			_mm_store_si128((__m128i*)fastptr+0x6, block[0x6]);
			_mm_store_si128((__m128i*)fastptr+0x7, block[0x7]);
			_mm_store_si128((__m128i*)fastptr+0x8, block[0x8]);
			_mm_store_si128((__m128i*)fastptr+0x9, block[0x9]);
			_mm_store_si128((__m128i*)fastptr+0xA, block[0xA]);
			_mm_store_si128((__m128i*)fastptr+0xB, block[0xB]);
			_mm_store_si128((__m128i*)fastptr+0xC, block[0xC]);
			_mm_store_si128((__m128i*)fastptr+0xD, block[0xD]);
			_mm_store_si128((__m128i*)fastptr+0xE, block[0xE]);
			_mm_store_si128((__m128i*)fastptr+0xF, block[0xF]);

			fastptr+=sizeof(__m128i[NCODERS]);
			{
				__m256i mp[4];
				mp[0]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+0));
				mp[1]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+1));
				mp[2]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+2));
				mp[3]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+3));
				_mm256_store_si256((__m256i*)slowptrs+0, mp[0]);
				_mm256_store_si256((__m256i*)slowptrs+1, mp[1]);
				_mm256_store_si256((__m256i*)slowptrs+2, mp[2]);
				_mm256_store_si256((__m256i*)slowptrs+3, mp[3]);
			}
		}
#endif
#if 1
		for(;kx<blockxbytes;kx+=NCODERS)
		{
			/*
			toy example
			CWH 3*4*5 original
			4 coders
			2*2 blocks
			4 Y's then 4 U's then 4 V's

			original	channel[block][pixel]
			r00 g00 b00 r01 g01 b01  r10 g10 b10 r11 g11 b11
			r02 g02 b02 r03 g03 b03  r12 g12 b12 r13 g13 b13
			r20 g20 b20 r21 g21 b21  r30 g30 b30 r31 g31 b31
			r22 g22 b22 r23 g23 b23  r32 g32 b32 r33 g33 b33

			interleaved
			r00 r10 r20 r30 g00 g10 g20 g30 b00 b10 b20 b30  r01 r11 r21 r31 g01 g11 g21 g31 b01 b11 b21 b31
			r02 r12 r22 r32 g02 g12 g22 g32 b02 b12 b22 b32  r03 r13 r23 r33 g03 g13 g23 g33 b03 b13 b23 b33
			*/
			fastptr[0x00]=*slowptrs[0x00]++;
			fastptr[0x01]=*slowptrs[0x01]++;
			fastptr[0x02]=*slowptrs[0x02]++;
			fastptr[0x03]=*slowptrs[0x03]++;
			fastptr[0x04]=*slowptrs[0x04]++;
			fastptr[0x05]=*slowptrs[0x05]++;
			fastptr[0x06]=*slowptrs[0x06]++;
			fastptr[0x07]=*slowptrs[0x07]++;
			fastptr[0x08]=*slowptrs[0x08]++;
			fastptr[0x09]=*slowptrs[0x09]++;
			fastptr[0x0A]=*slowptrs[0x0A]++;
			fastptr[0x0B]=*slowptrs[0x0B]++;
			fastptr[0x0C]=*slowptrs[0x0C]++;
			fastptr[0x0D]=*slowptrs[0x0D]++;
			fastptr[0x0E]=*slowptrs[0x0E]++;
			fastptr[0x0F]=*slowptrs[0x0F]++;
			fastptr+=NCODERS;
//#if defined __GNUC__ && !defined PROFILER
//#pragma GCC unroll 32
//#endif
//			for(int k=0;k<NCODERS;++k)
//				*fastptr++=*slowptrs[k]++;
		}
#endif
		for(int k=0;k<NCODERS;++k)
			slowptrs0[k]+=rowstride;
	}
}
static void interleave_blocks_inv(const unsigned char *interleaved, int iw, int ih, unsigned char *original)
{
	//original[ih][iw][3]
	//interleaved[ih/YCODERS][iw/XCODERS][3][NCODERS]	contiguous & aligned
	//xrem[ih%YCODERS][iw][3]
	//yrem[ih-iw%YCODERS][iw%XCODERS][3]

	//only difference between fwd and inv:		swap assignments (const slow->fast)
	int rowstride=3*iw;
	int ixyblockw=iw/XCODERS;
	int ixyblockh=ih/YCODERS;
	int blockxbytes=3*NCODERS*ixyblockw;
#ifdef INTERLEAVESIMD
	int SIMDxcount=blockxbytes&~((int)sizeof(__m128i[NCODERS])-1);
	__m256i slowinc=_mm256_set1_epi64x(sizeof(__m128i));
#endif
	const unsigned char *fastptr=interleaved;
	ALIGN(32) unsigned char *slowptrs[NCODERS]={0}, *slowptrs0[NCODERS]={0};
	for(int ky=0;ky<YCODERS;++ky)//spread slow pointers
	{
		for(int kx=0;kx<XCODERS;++kx)
			slowptrs0[XCODERS*ky+kx]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);
	}
	for(int ky=0;ky<ixyblockh;++ky)//interleave
	{
		int kx=0;
		memcpy((void*)slowptrs, slowptrs0, sizeof(slowptrs));
#ifdef INTERLEAVESIMD
		for(;kx<SIMDxcount;kx+=(int)sizeof(__m128i[NCODERS]))
		{
			/*
			SIMD interleave
			|000111222333444555666777888999AAABBBCCCDDDEEEFFFGGGHHHIIIJJJKKKLLLMMMNNN	24 pixels * 3 channels
			|^        ^        ^        ^        ^        ^        ^        ^		8-lane serial interleave
			|^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^		8-lane 4xSIMD interleave step 1
			|    ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^     ^^^^	8-lane 4xSIMD interleave step 2
			|        ^        ^        ^        ^        ^        ^        ^        ^	8-lane serial remainder
			|0369CFIL 0369CFIL 0369CFIL 147ADGJM 147ADGJM 147ADGJM 258BEHKN 258BEHKN 258BEHKN	serial behavior
			|0369CFIL0369CFIL0369CFIL147ADGJM  147ADGJM147ADGJM258BEHKN258BEHKN  2 5 8 B E H K N	SIMD quotients & serial remainder

			slowptrs[NCODERS]		fastptr  (aligned because NCODERS == sizeof(__m256i))
			A16x16 B16x16		->	A16x16T C16x16T B16x16T D16x16T
			C16x16 D16x16

			load, transpose, store blocks A, B
			load, transpose, store blocks C, D

			speed:
			SIMD maximum estimate	(3/load + 3/store)*32+80*0.5*2 = 272 cycles for 32*32 bytes  ->  3GHz*32*32/272 ~= 10770 MB/s without remainders, cache misses
			actual SIMD		MSVC 4900 MB/s		GCC 3900/3000 MB/s
			actual serial		1900 MB/s
			*/
			__m128i block[16];
			block[0x0]=_mm_load_si128((__m128i*)fastptr+0x0);
			block[0x1]=_mm_load_si128((__m128i*)fastptr+0x1);
			block[0x2]=_mm_load_si128((__m128i*)fastptr+0x2);
			block[0x3]=_mm_load_si128((__m128i*)fastptr+0x3);
			block[0x4]=_mm_load_si128((__m128i*)fastptr+0x4);
			block[0x5]=_mm_load_si128((__m128i*)fastptr+0x5);
			block[0x6]=_mm_load_si128((__m128i*)fastptr+0x6);
			block[0x7]=_mm_load_si128((__m128i*)fastptr+0x7);
			block[0x8]=_mm_load_si128((__m128i*)fastptr+0x8);
			block[0x9]=_mm_load_si128((__m128i*)fastptr+0x9);
			block[0xA]=_mm_load_si128((__m128i*)fastptr+0xA);
			block[0xB]=_mm_load_si128((__m128i*)fastptr+0xB);
			block[0xC]=_mm_load_si128((__m128i*)fastptr+0xC);
			block[0xD]=_mm_load_si128((__m128i*)fastptr+0xD);
			block[0xE]=_mm_load_si128((__m128i*)fastptr+0xE);
			block[0xF]=_mm_load_si128((__m128i*)fastptr+0xF);
			transpose16(block);
			_mm_storeu_si128((__m128i*)slowptrs[0x00], block[0x0]);
			_mm_storeu_si128((__m128i*)slowptrs[0x01], block[0x1]);
			_mm_storeu_si128((__m128i*)slowptrs[0x02], block[0x2]);
			_mm_storeu_si128((__m128i*)slowptrs[0x03], block[0x3]);
			_mm_storeu_si128((__m128i*)slowptrs[0x04], block[0x4]);
			_mm_storeu_si128((__m128i*)slowptrs[0x05], block[0x5]);
			_mm_storeu_si128((__m128i*)slowptrs[0x06], block[0x6]);
			_mm_storeu_si128((__m128i*)slowptrs[0x07], block[0x7]);
			_mm_storeu_si128((__m128i*)slowptrs[0x08], block[0x8]);
			_mm_storeu_si128((__m128i*)slowptrs[0x09], block[0x9]);
			_mm_storeu_si128((__m128i*)slowptrs[0x0A], block[0xA]);
			_mm_storeu_si128((__m128i*)slowptrs[0x0B], block[0xB]);
			_mm_storeu_si128((__m128i*)slowptrs[0x0C], block[0xC]);
			_mm_storeu_si128((__m128i*)slowptrs[0x0D], block[0xD]);
			_mm_storeu_si128((__m128i*)slowptrs[0x0E], block[0xE]);
			_mm_storeu_si128((__m128i*)slowptrs[0x0F], block[0xF]);

			fastptr+=sizeof(__m128i[NCODERS]);
			{
				__m256i mp[4];
				mp[0]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+0));
				mp[1]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+1));
				mp[2]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+2));
				mp[3]=_mm256_add_epi64(slowinc, _mm256_load_si256((__m256i*)slowptrs+3));
				_mm256_store_si256((__m256i*)slowptrs+0, mp[0]);
				_mm256_store_si256((__m256i*)slowptrs+1, mp[1]);
				_mm256_store_si256((__m256i*)slowptrs+2, mp[2]);
				_mm256_store_si256((__m256i*)slowptrs+3, mp[3]);
			}
		}
#endif
#if 1
		for(;kx<blockxbytes;kx+=NCODERS)
		{
			/*
			toy example
			CWH 3*4*5 original
			4 coders
			2*2 blocks
			4 Y's then 4 U's then 4 V's

			original	channel[block][pixel]
			r00 g00 b00 r01 g01 b01  r10 g10 b10 r11 g11 b11
			r02 g02 b02 r03 g03 b03  r12 g12 b12 r13 g13 b13
			r20 g20 b20 r21 g21 b21  r30 g30 b30 r31 g31 b31
			r22 g22 b22 r23 g23 b23  r32 g32 b32 r33 g33 b33

			interleaved
			r00 r10 r20 r30 g00 g10 g20 g30 b00 b10 b20 b30  r01 r11 r21 r31 g01 g11 g21 g31 b01 b11 b21 b31
			r02 r12 r22 r32 g02 g12 g22 g32 b02 b12 b22 b32  r03 r13 r23 r33 g03 g13 g23 g33 b03 b13 b23 b33
			*/
			*slowptrs[0x00]++=fastptr[0x00];
			*slowptrs[0x01]++=fastptr[0x01];
			*slowptrs[0x02]++=fastptr[0x02];
			*slowptrs[0x03]++=fastptr[0x03];
			*slowptrs[0x04]++=fastptr[0x04];
			*slowptrs[0x05]++=fastptr[0x05];
			*slowptrs[0x06]++=fastptr[0x06];
			*slowptrs[0x07]++=fastptr[0x07];
			*slowptrs[0x08]++=fastptr[0x08];
			*slowptrs[0x09]++=fastptr[0x09];
			*slowptrs[0x0A]++=fastptr[0x0A];
			*slowptrs[0x0B]++=fastptr[0x0B];
			*slowptrs[0x0C]++=fastptr[0x0C];
			*slowptrs[0x0D]++=fastptr[0x0D];
			*slowptrs[0x0E]++=fastptr[0x0E];
			*slowptrs[0x0F]++=fastptr[0x0F];
			fastptr+=NCODERS;
//#if defined __GNUC__ && !defined PROFILER
//#pragma GCC unroll 32
//#endif
//			for(int k=0;k<NCODERS;++k)
//				*slowptrs[k]++=*fastptr++;
		}
#endif
		for(int k=0;k<NCODERS;++k)
			slowptrs0[k]+=rowstride;
	}
}
static void save_ppm(const char *fn, const unsigned char *image, int iw, int ih)
{
	FILE *fdst=fopen(fn, "wb");
	if(!fdst)
	{
		CRASH("Cannot open \"%s\" for writing", fn);
		return;
	}
	fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	fwrite(image, 1, (ptrdiff_t)3*iw*ih, fdst);
	fclose(fdst);
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
			if(streamptr>streamend)
				CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			state=state<<16|*(unsigned short*)streamptr;
			streamptr+=2;
		}

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
			if(streamptr>streamend)
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
			if(streamptr>streamend)
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
static void encode1d_fse(unsigned char *data, int count, int bytestride, BitPackerLIFO *ec, tANS_CState_t *states)
{
	unsigned char *ptr=data+(count-(ptrdiff_t)1)*bytestride;
	for(int k=0;k<count;++k)
	{
		ansval_push(ptr+2, 1, 1); tANS_encodeSymbol(ec, states+2, ptr[2]);
		ansval_push(ptr+1, 1, 1); tANS_encodeSymbol(ec, states+1, ptr[1]);
		ansval_push(ptr+0, 1, 1); tANS_encodeSymbol(ec, states+0, ptr[0]);
		ptr-=bytestride;
	}
}
static void decode1d_fse(unsigned char *data, int count, int bytestride, int bestrct, BitPackerLIFO *ec, tANS_DState_t *states)
{
	const unsigned char *combination=rct_combinations[bestrct];
	int yidx=combination[II_PERM_Y];
	int uidx=combination[II_PERM_U];
	int vidx=combination[II_PERM_V];
	int ufromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	int vc0=combination[II_COEFF_V_SUB_Y];
	int vc1=combination[II_COEFF_V_SUB_U];

	unsigned char *ptr=data;
	int prevy=0, prevu=0, prevv=0, offset1=0, offset2=0;
	int y=0, u=0, v=0;
	for(int k=0;k<count;++k)
	{
		y=tANS_decodeSymbol(states+0, ec); ansval_check(&y, 1, 1);
		u=tANS_decodeSymbol(states+1, ec); ansval_check(&u, 1, 1);
		v=tANS_decodeSymbol(states+2, ec); ansval_check(&v, 1, 1);
		
		y=(int8_t)(y+prevy-128);

		offset1=y&ufromy;
		prevu+=offset1;
		CLAMP2(prevu, -128, 127);
		u=(int8_t)(u+prevu-128);

		offset2=vc0*y+vc1*u;
		int vpred=(prevv+offset2)>>2;
		CLAMP2(vpred, -128, 127);
		v=(int8_t)(v+vpred-128);

		prevy=y;
		prevu=u-offset1;
		prevv=4*v-offset2;

		ptr[yidx]=y+128;
		ptr[uidx]=u+128;
		ptr[vidx]=v+128;
		ptr+=bytestride;
	}
}
int c40_codec(int argc, char **argv)
{
	if(argc!=3&&argc!=4&&argc!=5&&argc!=6)
	{
		printf(
			"Usage: \"%s\"  input  output  [Effort]  [Dist]    Encode/decode.\n"
			"  Effort  =  0 CG / 1~3 L1  +  10 Profiler  +  100 FSE (otherwise, rANS).\n"
			"  Dist    =  lossy distortion. 5 <= Dist <= 17.\n"
			"Built on %s %s\n"
			, argv[0]
			, __DATE__, __TIME__
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];
	int param1=argc<4?DEFAULT_EFFORT_LEVEL:atoi(argv[3]), dist=argc<5?1:atoi(argv[4]);
	int effort=param1%10, profile=param1/10%10, userans=!(param1/100);
	if(dist>1)
		CLAMP2(dist, 5, 17);
#ifdef ESTIMATE_SIZE
	double esize[3*NCODERS]={0};
#endif
#ifdef LOUD
	ptrdiff_t usize2=0;
	{
		struct stat info={0};
		stat(srcfn, &info);
		usize2=info.st_size;
	}
	double t=time_sec();
#endif
	prof_checkpoint(0, 0);
	if(!srcfn||!dstfn)
	{
		CRASH("Codec requires both source and destination filenames");
		return 1;
	}
	int fwd=0, iw=0, ih=0, rowstride=0;
	ptrdiff_t usize=0, cap=0;
	unsigned char *image=0, *imptr=0, *streamptr=0, *streamstart=0, *streamend=0;
	int psize=0;
	short *pixels=0;
	ptrdiff_t cheadersize=0, csize=0;
	uint64_t degenmask=0;
	uint64_t bypassmask=0;//0: rare context (bypass)  1: emit stats		3*NCTX+3 = 57 flags
	int bestrct=0, npreds=0, sh=0;
	size_t fseerror=0;
	uint32_t *fsetables[3*NCTX+3]={0};
	tANS_CState_t fsecstates[3*NCTX+3]={0};
	tANS_DState_t fsedstates[3*NCTX+3]={0};
	BitPackerLIFO ec={0};
	{
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			CRASH("Cannot open \"%s\"", srcfn);
			return 1;
		}
		int tag=0;
		fread(&tag, 1, 2, fsrc);
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('F'|'2'<<8))
		{
			CRASH("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
		{
#ifdef LOUD
			{
				char buf[1024]={0};
				time_t tstamp=time(0);
				struct tm *tformat=localtime(&tstamp);
				strftime(buf, sizeof(buf)-1, "%Y-%m-%d_%H%M%S\n", tformat);
			}
			//print_timestamp("%Y-%m-%d_%H%M%S\n");
#endif
			int temp=fgetc(fsrc);
			if(temp!='\n')
			{
				CRASH("Invalid PPM file");
				return 1;
			}
			int nread=fscanf(fsrc, "%d %d", &iw, &ih);
			if(nread!=2)
			{
				CRASH("Unsupported PPM file");
				return 1;
			}
			int vmax=0;
			nread=fscanf(fsrc, "%d", &vmax);
			if(nread!=1||vmax!=255)
			{
				CRASH("Unsupported PPM file");
				return 1;
			}
			temp=fgetc(fsrc);
			if(temp!='\n')
			{
				CRASH("Invalid PPM file");
				return 1;
			}
		}
		else
		{
			int flags=0;

			iw=0;
			ih=0;
			dist=0;
			bypassmask=0;
			degenmask=0;
			fread(&iw, 1, 3, fsrc);
			fread(&ih, 1, 3, fsrc);
			fread(&flags, 1, 2, fsrc);
			bestrct=flags>>3;
			userans=flags>>2&1;
			effort=flags&3;
			fread(&dist, 1, 1, fsrc);
			if(!userans)
				fread(&degenmask, 1, 8, fsrc);
			fread(&bypassmask, 1, 8, fsrc);
			cheadersize=ftell(fsrc);
		}
		if(iw<1||ih<1)
		{
			CRASH("Unsupported source file");
			return 1;
		}
		rowstride=3*iw;
		usize=(ptrdiff_t)3*iw*ih;
		cap=(ptrdiff_t)4*iw*ih;
	//	image=(unsigned char*)_mm_malloc(cap+sizeof(__m256i), 0x1000);
		image=(unsigned char*)malloc(cap+sizeof(__m256i));
		if(!image)
		{
			CRASH("Alloc error");
			return 1;
		}
		if(fwd)
		{
			fread(image, 1, usize, fsrc);//read image
			streamstart=image;
			streamend=image+cap;
			streamptr=streamend;//bwd-bwd ANS encoding
			profile_size(streamptr, "start");
		}
		else
		{
			struct stat info={0};
			stat(srcfn, &info);
			csize=info.st_size;
			streamstart=image+cap-(csize-cheadersize)-sizeof(__m256i);
			streamend=image+cap-sizeof(__m256i);
			streamptr=streamstart;
			fread(streamstart, 1, csize-cheadersize, fsrc);//read stream
		}
		fclose(fsrc);
	}
	prof_checkpoint(fwd?usize:csize, "fread");
	int blockw=iw/XCODERS;
	int blockh=ih/YCODERS;
	int qxbytes=blockw*XCODERS*3;//iw/XCODERS*XCODERS*3
	int ixcount=blockw*NCODERS, ixbytes=3*ixcount;//ix = interleaved circular buffer width		iw/XCODERS*NCODERS
	int xremw=iw-blockw*XCODERS, yremh=ih-blockh*YCODERS;
	int xrembytes=3*xremw;
	int nctx=3*NCTX+3*(xremw||yremh);
	ptrdiff_t isize=(ptrdiff_t)ixbytes*blockh;
	ptrdiff_t interleavedsize=isize<<fwd;//fwd ? interleave residuals & context : pack residuals
	unsigned char *interleaved=(unsigned char*)_mm_malloc(interleavedsize, sizeof(__m256i));
	if(!interleaved)
	{
		CRASH("Alloc error");
		return 1;
	}
	(void)xrembytes;
	const int hsize=nctx*(int)sizeof(int[256]);//3 channels
	int *hists=fwd?(int*)malloc(hsize):0;//fwd-only

	int CDF2syms_size=nctx*(int)sizeof(int[1<<PROBBITS]);
	if(fwd)//DIV-free rANS encoder reuses this as SIMD symbol info
		CDF2syms_size=nctx*(int)sizeof(rANS_SIMD_SymInfo[256]);
	unsigned *CDF2syms=(unsigned*)_mm_malloc(CDF2syms_size, sizeof(__m256i));
	
	__m256i mstate[2];
	int ans_permute_size=sizeof(__m256i[256]);
	int *ans_permute=userans?(int*)_mm_malloc(ans_permute_size, sizeof(__m256i)):0;

	psize=(int)sizeof(short[4*6*NCODERS])*(blockw+16);//4 padded rows  *  {Y*NCODERS, U*NCODERS, V*NCODERS,  eY*NCODERS, eU*NCODERS, eV*NCODERS} = 2*3*32 = 192 channels  ~48*iw bytes
	pixels=(short*)_mm_malloc(psize, sizeof(__m256i));//~188 KB for 4K/12MP
	if((fwd&&!hists)||(userans&&(!CDF2syms||!ans_permute))||!pixels)
	{
		CRASH("Alloc error");
		return 1;
	}
	if(userans)
		memset(ans_permute, 0, ans_permute_size);//_mm256_permutevar8x32_epi32 can't clear elements like _mm256_shuffle_epi8 does
	else
	{
		for(int k=0;k<nctx;++k)
		{
			if(fwd)
				fsetables[k]=tANS_createCTable(255, PROBBITS);
			else
				fsetables[k]=tANS_createDTable(PROBBITS);
			if(!fsetables[k])
			{
				CRASH("Alloc error");
				return 1;
			}
		}
	}
	if(fwd)
	{
		memset(hists, 0, hsize);
#ifdef TEST_INTERLEAVE
		guide_save(image, iw, ih);
		save_ppm("20250227_1225AM_original.PPM", image, iw, ih);
		interleave_blocks_fwd(image, iw, ih, interleaved);
		save_ppm("20250226_1153PM_interleaved.PPM", interleaved, ixcount, blockh);
		interleave_blocks_inv(interleaved, iw, ih, image);
		save_ppm("20250227_1244AM_deinterleaved.PPM", image, iw, ih);
		if(memcmp(image, g_image, usize))
			CRASH("ERROR");
		printf("SUCCESS\n");
		exit(0);
#endif
		interleave_blocks_fwd(image, iw, ih, interleaved+isize);//reuse memory: read 8-bit pixel, write 16-bit context<<8|residual
		guide_save(interleaved+isize, ixcount, blockh);
		prof_checkpoint(usize, "interleave");

		//analysis
		ALIGN(32) long long counters[OCH_COUNT]={0};
		__m256i mcounters[OCH_COUNT];//64-bit
		__m128i half8=_mm_set1_epi8(-128);
		__m256i wordmask=_mm256_set1_epi64x(0xFFFF);
		memset(mcounters, 0, sizeof(mcounters));
		imptr=interleaved+isize;
		for(int ky=0;ky<blockh;ky+=ANALYSIS_YSTRIDE)//analysis
		{
			__m256i prev[OCH_COUNT];//16-bit
			memset(prev, 0, sizeof(prev));
			for(int kx=0;kx<ixbytes-3*NCODERS;kx+=3*NCODERS*ANALYSIS_XSTRIDE)
			{
				__m256i r=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+0), half8));
				__m256i g=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+1), half8));
				__m256i b=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+2), half8));
				imptr+=3*NCODERS*ANALYSIS_XSTRIDE;
				r=_mm256_slli_epi16(r, 2);
				g=_mm256_slli_epi16(g, 2);
				b=_mm256_slli_epi16(b, 2);
				__m256i rg=_mm256_sub_epi16(r, g);
				__m256i gb=_mm256_sub_epi16(g, b);
				__m256i br=_mm256_sub_epi16(b, r);
#ifdef ENABLE_RCT_EXTENSION
				__m256i t0, t1, t2;
#endif
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
#ifdef ENABLE_RCT_EXTENSION
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
#endif
			}
			imptr+=ixbytes*(ANALYSIS_YSTRIDE-1);
		}
		for(int k=0;k<OCH_COUNT;++k)
		{
			ALIGN(32) long long temp[4]={0};
			_mm256_store_si256((__m256i*)temp, mcounters[k]);
			counters[k]=temp[0]+temp[1]+temp[2]+temp[3];
		}
		long long minerr=0;
		for(int kt=0;kt<RCT_COUNT;++kt)
		{
			const unsigned char *rct=rct_combinations[kt];
			long long currerr=
				+counters[rct[0]]
				+counters[rct[1]]
				+counters[rct[2]]
			;
#ifdef LOUD
			printf("%-14s %12lld + %12lld + %12lld = %12lld%s\n",
				rct_names[kt],
				counters[rct[0]],
				counters[rct[1]],
				counters[rct[2]],
				currerr,
				!kt||minerr>currerr?" <-":""
			);
#endif
			if(!kt||minerr>currerr)
			{
				minerr=currerr;
				bestrct=kt;
			}
		}
		//printf("%2d ", bestrct);
		prof_checkpoint(usize, "analysis");

		if(userans)
		{
			//generate encode permutations		eg: mask = MSB 0b11000101 LSB  ->  LO {x, x, x, x, 0, 2, 6, 7} HI
			for(int km=0;km<256;++km)
			{
				int *curr=ans_permute+((ptrdiff_t)km<<3);
				int kb2=7;
				for(int kb=7;kb>=0;--kb)
				{
					int bit=km>>kb&1;
					if(bit)
						curr[kb2--]=kb;
				}
			}
			prof_checkpoint(ans_permute_size, "gen permutation");
		}
	}
	else
	{
		if(userans)
		{
			//BitPackerLIFO ec;
			bitpacker_dec_init(&ec, streamptr, streamend);
			for(int kc=0;kc<nctx;++kc)
			{
				uint16_t hist[256];
				dec_unpackhist(&ec, bypassmask, kc, hist);
				dec_hist2CDF2sym(hist, CDF2syms+((ptrdiff_t)kc<<PROBBITS));
			}
			streamptr=(unsigned char*)(size_t)ec.srcfwdptr;
			prof_checkpoint((ptrdiff_t)CDF2syms_size, "unpack histograms");

			//generate decode permutations		eg: mask = MSB 0b11000101 LSB  ->  LO {0, x, 1, x, x, x, 2, 3} HI
			for(int km=0;km<256;++km)
			{
				int *curr=ans_permute+((ptrdiff_t)km<<3);
				int idx=0;
				for(int kb=0;kb<8;++kb)
				{
					int bit=km>>kb&1;
					if(bit)
						curr[kb]=idx++;
				}
			}
			prof_checkpoint(ans_permute_size, "gen permutation");

#ifdef _DEBUG
			if(streamptr>streamend)
				CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			memcpy(mstate, streamptr, sizeof(mstate));
			streamptr+=sizeof(mstate);
		}
		else
		{
			bitpacker_dec_init(&ec, streamptr, streamend);
			for(int kc=0;kc<nctx;++kc)
			{
				uint16_t nhist[256]={0};
				if(degenmask>>kc&1)//degenerate distribution
				{
					int sym=bitpacker_dec(&ec, 8);
					fseerror=tANS_buildDTable_rle(fsetables[kc], sym);
					if(fseerror)
					{
						CRASH("FSE Error");
						return 1;
					}
				}
				else if(bypassmask>>kc&1)//bypass distribution
				{
					fseerror=tANS_buildDTable_raw(fsetables[kc], 8);
					if(fseerror)
					{
						CRASH("FSE Error");
						return 1;
					}
				}
				else//ordinary distribution
				{
					uint32_t vmax=255, probbits=PROBBITS;
					dec_unpackhist(&ec, bypassmask, kc, nhist);
					
					fseerror=tANS_buildDTable(fsetables[kc], (short*)nhist, vmax, probbits);
					if(fseerror)
					{
						CRASH("FSE Error");
						return 1;
					}
				}
#if defined _DEBUG && 0
				if(kc==0)
				{
					printf("dec.ctx = 0:\n");
					for(int ks=0;ks<256;++ks)
						printf("%3d  %04X\n", ks, nhist[ks]);
					printf("\n");
				}
#endif
			}
#if defined _DEBUG && 0
			printf("probbits:\n");
#endif
			for(int kc=0;kc<nctx;++kc)
			{
				tANS_initDState(fsedstates+kc, &ec, fsetables[kc]);
#if defined _DEBUG && 0
				const tANS_DTableHeader *const DTableH=(const tANS_DTableHeader*)fsetables[kc];
				printf("%c %2d %2d load 0x%04tX %s\n"
					, "YUVR"[kc/NCTX]
					, kc%NCTX
					, DTableH->tableLog
					, fsedstates[kc].state
					, degenmask>>kc&1?"degen":(bypassmask>>kc&1?"bypass":"normal")
				);
#endif
			}
		}
	}
	switch(effort)
	{
	case 0://use CG
		npreds=0;
		break;
	case 1://use L1
		npreds=L1_NPREDS1;
		sh=L1_SH1;
		break;
	case 2:
		npreds=L1_NPREDS2;
		sh=L1_SH2;
		break;
	case 3:
		npreds=L1_NPREDS3;
		sh=L1_SH3;
		break;
	}
#ifdef LOUD
	if(fwd)
#else
	if(profile)
#endif
		printf("%s  NPREDS=%d  %td bytes\n", rct_names[bestrct], npreds, usize);
	int L1statesize=0;
	int *L1state=0;
	if(effort)
	{
		L1statesize=(int)sizeof(int[2*NCODERS*3*(L1_NPREDS3+1)]);//{preds, coeffs} * (NPREDS+{bias}) * 3 channels * NCODERS
		L1state=(int*)_mm_malloc(L1statesize, sizeof(__m256i));
		if(!L1state)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(L1state, 0, L1statesize);
	}
#ifdef PRINT_L1_BOUNDS
	int cmin=0, cmax=0;
	int bmin=0, bmax=0;
#endif
	const unsigned char *combination=rct_combinations[bestrct];
	int
		yidx=combination[II_PERM_Y]*NCODERS,
		uidx=combination[II_PERM_U]*NCODERS,
		vidx=combination[II_PERM_V]*NCODERS;
	__m256i uhelpmask=_mm256_set1_epi16(-(combination[II_COEFF_U_SUB_Y]!=0));
	__m256i vc0=_mm256_set1_epi16(combination[II_COEFF_V_SUB_Y]);
	__m256i vc1=_mm256_set1_epi16(combination[II_COEFF_V_SUB_U]);
	int paddedwidth=blockw+16;
	memset(pixels, 0, psize);
	__m256i mctxmax=_mm256_set1_epi16(NCTX-1);
	__m256i mctxuoffset=_mm256_set1_epi16(NCTX);
	__m256i mctxvoffset=_mm256_set1_epi16(NCTX*2);
	__m256i amin=_mm256_set1_epi16(-128);
	__m256i amax=_mm256_set1_epi16(127);
	__m128i half8=_mm_set1_epi8(-128);
	__m256i bytemask=_mm256_set1_epi16(255);
	__m256i wordmask=_mm256_set1_epi32(0xFFFF);
	__m256i myuv[3];
	__m256i dist_rcp=_mm256_set1_epi16(0x7FFF), mdist=_mm256_set1_epi16(1);
#ifdef SAVE_RESIDUALS
	unsigned char *residuals=0;
	if(fwd)
	{
		residuals=(unsigned char*)malloc(isize);
		if(!residuals)
		{
			CRASH("Alloc error");
			return 1;
		}
		memset(residuals, 0, isize);
	}
#endif
	if(dist>1)
	{
		dist_rcp=_mm256_set1_epi16(((1<<16)+dist-1)/dist);//x/dist  ->  {x*=inv; x=(x>>16)+((unsigned)x>>31);}
		mdist=_mm256_set1_epi16(dist);
	}
	memset(myuv, 0, sizeof(myuv));
	unsigned char *ctxptr=interleaved;
	imptr=interleaved+(fwd?isize:0);
	__m256i *L1preds=effort?(__m256i*)L1state:0;
	int *L1weights=effort?(int*)(L1state+1*(ptrdiff_t)NCODERS*3*(L1_NPREDS3+1)):0;
	if(effort)
	{
		for(int k=0, count=3*NCODERS*(npreds+1), init=(1<<sh)/npreds;k<count;++k)
			L1weights[k]=init;
	}
	for(int ky=0;ky<blockh;++ky)//main coding loop
	{
		ALIGN(32) short *rows[]=
		{
			pixels+(paddedwidth*((ky-0LL)&3)+8LL)*6*NCODERS,
			pixels+(paddedwidth*((ky-1LL)&3)+8LL)*6*NCODERS,
			pixels+(paddedwidth*((ky-2LL)&3)+8LL)*6*NCODERS,
			pixels+(paddedwidth*((ky-3LL)&3)+8LL)*6*NCODERS,
		};
		ALIGN(32) unsigned short syms[3*NCODERS]={0};
		__m256i NW[3], N[3], W[3];
		__m256i eW[3], ecurr[3], eNEE[3], eNEEE[3];
		memset(NW, 0, sizeof(NW));
		memset(N, 0, sizeof(N));
		memset(W, 0, sizeof(W));
		memset(eW, 0, sizeof(eW));
		memset(ecurr, 0, sizeof(ecurr));
		memset(eNEE, 0, sizeof(eNEE));
		memset(eNEEE, 0, sizeof(eNEEE));
		//                     (__m256i*)rows[-Y]+E*3+C+X*6
		eNEE[0]=_mm256_load_si256((__m256i*)rows[1]+1*3+0+2*6);
		eNEE[1]=_mm256_load_si256((__m256i*)rows[1]+1*3+1+2*6);
		eNEE[2]=_mm256_load_si256((__m256i*)rows[1]+1*3+2+2*6);
		for(int kx=0;kx<blockw;++kx)
		{
			//                   (__m256i*)rows[1]+E*3+C+X*6
			N[0]=_mm256_load_si256((__m256i*)rows[1]+0+0+0*6);//y neighbors
			N[1]=_mm256_load_si256((__m256i*)rows[1]+0+1+0*6);//u
			N[2]=_mm256_load_si256((__m256i*)rows[1]+0+2+0*6);//v
			__m256i
				predY, ctxY,
				predU, ctxU,
				predV, ctxV;
			__m256i predYUV0[3];
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
			{
				const int borderW=3;
				const int borderN=3;
				const int borderE=3;
				int cond_cg=(unsigned)(kx-borderW)>=(unsigned)(blockw-(borderW+borderE))
					||(unsigned)(ky-borderN)>=(unsigned)(blockh-borderN);
				__m256i mcg[3];
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
				mcg[0]=predY;
				mcg[1]=predU;
				mcg[2]=predV;

				if(effort==1)//predict
				{
					/*
					effort 1
					0	N+W-NW
					1	N
					2	W
					3	W+NE-N
					*/

					//N+W-NW
					L1preds[0*3+0]=predY;
					L1preds[0*3+1]=predU;
					L1preds[0*3+2]=predV;

					//N
					L1preds[1*3+0]=N[0];
					L1preds[1*3+1]=N[1];
					L1preds[1*3+2]=N[2];

					//W
					L1preds[2*3+0]=W[0];
					L1preds[2*3+1]=W[1];
					L1preds[2*3+2]=W[2];

					//W+NE-N
					L1preds[3*3+0]=_mm256_sub_epi16(_mm256_add_epi16(W[0], _mm256_load_si256((__m256i*)rows[1]+0+0+1*6)), N[0]);
					L1preds[3*3+1]=_mm256_sub_epi16(_mm256_add_epi16(W[1], _mm256_load_si256((__m256i*)rows[1]+0+1+1*6)), N[1]);
					L1preds[3*3+2]=_mm256_sub_epi16(_mm256_add_epi16(W[2], _mm256_load_si256((__m256i*)rows[1]+0+2+1*6)), N[2]);
					

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
					0	N+W-NW
					1	N
					2	W
					3	W+NE-N
					4	3*(N-NN)+NNN
					5	3*(W-WW)+WWW
					6	N+NE-NNE
					7	NEE
					8	NN
					9	WW
					10	2*N-NN
					11	2*W-WW
					12	NEEE
					13	NEEEE
					14	NNWW
					15	NNEE
					16	N+NW-NNW
					17	W+NW-NWW
					18	(WWWW+NEEEE)>>1
					19	(WWW+NNN+NEEE-NW)>>1

					(__m256i*)rows[-Y]+E*3+C+X*6
					*/

					//N+W-NW
					L1preds[0*3+0]=predY;
					L1preds[0*3+1]=predU;
					L1preds[0*3+2]=predV;

					//N
					L1preds[1*3+0]=N[0];
					L1preds[1*3+1]=N[1];
					L1preds[1*3+2]=N[2];

					//W
					L1preds[2*3+0]=W[0];
					L1preds[2*3+1]=W[1];
					L1preds[2*3+2]=W[2];

					//W+NE-N
					cache[0]=_mm256_add_epi16(W[0], _mm256_load_si256((__m256i*)rows[1]+0+0+1*6));
					cache[1]=_mm256_add_epi16(W[1], _mm256_load_si256((__m256i*)rows[1]+0+1+1*6));
					cache[2]=_mm256_add_epi16(W[2], _mm256_load_si256((__m256i*)rows[1]+0+2+1*6));
					L1preds[3*3+0]=_mm256_sub_epi16(cache[0], N[0]);
					L1preds[3*3+1]=_mm256_sub_epi16(cache[1], N[1]);
					L1preds[3*3+2]=_mm256_sub_epi16(cache[2], N[2]);

					//3*(N-NN)+NNN
					cache[0]=_mm256_sub_epi16(N[0], _mm256_load_si256((__m256i*)rows[2]+0+0+0*6));//N-NN
					cache[1]=_mm256_sub_epi16(N[1], _mm256_load_si256((__m256i*)rows[2]+0+1+0*6));
					cache[2]=_mm256_sub_epi16(N[2], _mm256_load_si256((__m256i*)rows[2]+0+2+0*6));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
					L1preds[4*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[3]+0+0+0*6));//+NNN
					L1preds[4*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[3]+0+1+0*6));
					L1preds[4*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[3]+0+2+0*6));

					//3*(W-WW)+WWW
					cache[0]=_mm256_sub_epi16(W[0], _mm256_load_si256((__m256i*)rows[0]+0+0-2*6));//W-WW
					cache[1]=_mm256_sub_epi16(W[1], _mm256_load_si256((__m256i*)rows[0]+0+1-2*6));
					cache[2]=_mm256_sub_epi16(W[2], _mm256_load_si256((__m256i*)rows[0]+0+2-2*6));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
					L1preds[5*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[0]+0+0-3*6));//+WWW
					L1preds[5*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[0]+0+1-3*6));
					L1preds[5*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[0]+0+2-3*6));

					//N+NE-NNE
					cache[0]=_mm256_add_epi16(N[0], _mm256_load_si256((__m256i*)rows[1]+0+0+1*6));//N+NE
					cache[1]=_mm256_add_epi16(N[1], _mm256_load_si256((__m256i*)rows[1]+0+1+1*6));
					cache[2]=_mm256_add_epi16(N[2], _mm256_load_si256((__m256i*)rows[1]+0+2+1*6));
					L1preds[6*3+0]=_mm256_sub_epi16(cache[0], _mm256_load_si256((__m256i*)rows[2]+0+0+1*6));//NNE
					L1preds[6*3+1]=_mm256_sub_epi16(cache[1], _mm256_load_si256((__m256i*)rows[2]+0+1+1*6));
					L1preds[6*3+2]=_mm256_sub_epi16(cache[2], _mm256_load_si256((__m256i*)rows[2]+0+2+1*6));

					//NEE
					L1preds[7*3+0]=_mm256_load_si256((__m256i*)rows[1]+0+0+2*6);
					L1preds[7*3+1]=_mm256_load_si256((__m256i*)rows[1]+0+1+2*6);
					L1preds[7*3+2]=_mm256_load_si256((__m256i*)rows[1]+0+2+2*6);

					//NN
					L1preds[8*3+0]=_mm256_load_si256((__m256i*)rows[2]+0+0+0*6);
					L1preds[8*3+1]=_mm256_load_si256((__m256i*)rows[2]+0+1+0*6);
					L1preds[8*3+2]=_mm256_load_si256((__m256i*)rows[2]+0+2+0*6);

					//WW
					L1preds[9*3+0]=_mm256_load_si256((__m256i*)rows[0]+0+0-2*6);
					L1preds[9*3+1]=_mm256_load_si256((__m256i*)rows[0]+0+1-2*6);
					L1preds[9*3+2]=_mm256_load_si256((__m256i*)rows[0]+0+2-2*6);

					//2*N-NN
					L1preds[10*3+0]=_mm256_sub_epi16(_mm256_slli_epi16(N[0], 1), _mm256_load_si256((__m256i*)rows[2]+0+0+0*6));
					L1preds[10*3+1]=_mm256_sub_epi16(_mm256_slli_epi16(N[1], 1), _mm256_load_si256((__m256i*)rows[2]+0+1+0*6));
					L1preds[10*3+2]=_mm256_sub_epi16(_mm256_slli_epi16(N[2], 1), _mm256_load_si256((__m256i*)rows[2]+0+2+0*6));

					//2*W-WW
					L1preds[11*3+0]=_mm256_sub_epi16(_mm256_slli_epi16(W[0], 1), _mm256_load_si256((__m256i*)rows[0]+0+0-2*6));
					L1preds[11*3+1]=_mm256_sub_epi16(_mm256_slli_epi16(W[1], 1), _mm256_load_si256((__m256i*)rows[0]+0+1-2*6));
					L1preds[11*3+2]=_mm256_sub_epi16(_mm256_slli_epi16(W[2], 1), _mm256_load_si256((__m256i*)rows[0]+0+2-2*6));

					//NEEE
					L1preds[12*3+0]=_mm256_load_si256((__m256i*)rows[1]+0+0+3*6);
					L1preds[12*3+1]=_mm256_load_si256((__m256i*)rows[1]+0+1+3*6);
					L1preds[12*3+2]=_mm256_load_si256((__m256i*)rows[1]+0+2+3*6);

					//NEEEE
					L1preds[13*3+0]=_mm256_load_si256((__m256i*)rows[1]+0+0+4*6);
					L1preds[13*3+1]=_mm256_load_si256((__m256i*)rows[1]+0+1+4*6);
					L1preds[13*3+2]=_mm256_load_si256((__m256i*)rows[1]+0+2+4*6);

					//NNWW
					L1preds[14*3+0]=_mm256_load_si256((__m256i*)rows[2]+0+0-2*6);
					L1preds[14*3+1]=_mm256_load_si256((__m256i*)rows[2]+0+1-2*6);
					L1preds[14*3+2]=_mm256_load_si256((__m256i*)rows[2]+0+2-2*6);

					//NNEE
					L1preds[15*3+0]=_mm256_load_si256((__m256i*)rows[2]+0+0+2*6);
					L1preds[15*3+1]=_mm256_load_si256((__m256i*)rows[2]+0+1+2*6);
					L1preds[15*3+2]=_mm256_load_si256((__m256i*)rows[2]+0+2+2*6);

					//N+NW-NNW
					cache[0]=_mm256_add_epi16(N[0], NW[0]);//N+NW
					cache[1]=_mm256_add_epi16(N[1], NW[1]);
					cache[2]=_mm256_add_epi16(N[2], NW[2]);
					L1preds[16*3+0]=_mm256_sub_epi16(cache[0], _mm256_load_si256((__m256i*)rows[2]+0+0-1*6));//NNW
					L1preds[16*3+1]=_mm256_sub_epi16(cache[1], _mm256_load_si256((__m256i*)rows[2]+0+1-1*6));
					L1preds[16*3+2]=_mm256_sub_epi16(cache[2], _mm256_load_si256((__m256i*)rows[2]+0+2-1*6));

					//W+NW-NWW
					cache[0]=_mm256_add_epi16(W[0], NW[0]);//W+NW
					cache[1]=_mm256_add_epi16(W[1], NW[1]);
					cache[2]=_mm256_add_epi16(W[2], NW[2]);
					L1preds[17*3+0]=_mm256_sub_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+0-2*6));//NWW
					L1preds[17*3+1]=_mm256_sub_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+1-2*6));
					L1preds[17*3+2]=_mm256_sub_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+2-2*6));

					//(WWWW+NEEEE)>>1
					cache[0]=_mm256_add_epi16(_mm256_load_si256((__m256i*)rows[0]+0+0-4*6), _mm256_load_si256((__m256i*)rows[1]+0+0+4*6));//WWWW+NEEEE
					cache[1]=_mm256_add_epi16(_mm256_load_si256((__m256i*)rows[0]+0+1-4*6), _mm256_load_si256((__m256i*)rows[1]+0+1+4*6));
					cache[2]=_mm256_add_epi16(_mm256_load_si256((__m256i*)rows[0]+0+2-4*6), _mm256_load_si256((__m256i*)rows[1]+0+2+4*6));
					L1preds[18*3+0]=_mm256_srai_epi16(cache[0], 1);
					L1preds[18*3+1]=_mm256_srai_epi16(cache[1], 1);
					L1preds[18*3+2]=_mm256_srai_epi16(cache[2], 1);

					//(WWW+NNN+NEEE-NW)>>1
					cache[0]=_mm256_add_epi16(_mm256_load_si256((__m256i*)rows[0]+0+0-3*6), _mm256_load_si256((__m256i*)rows[3]+0+0+4*6));//WWW+NNN
					cache[1]=_mm256_add_epi16(_mm256_load_si256((__m256i*)rows[0]+0+1-3*6), _mm256_load_si256((__m256i*)rows[3]+0+1+4*6));
					cache[2]=_mm256_add_epi16(_mm256_load_si256((__m256i*)rows[0]+0+2-3*6), _mm256_load_si256((__m256i*)rows[3]+0+2+4*6));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+0+3*6));//+NEEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+1+3*6));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+2+3*6));
					cache[0]=_mm256_sub_epi16(cache[0], NW[0]);//-NW
					cache[1]=_mm256_sub_epi16(cache[1], NW[1]);
					cache[2]=_mm256_sub_epi16(cache[2], NW[2]);
					L1preds[19*3+0]=_mm256_srai_epi16(cache[0], 1);
					L1preds[19*3+1]=_mm256_srai_epi16(cache[1], 1);
					L1preds[19*3+2]=_mm256_srai_epi16(cache[2], 1);


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

				predY=_mm256_max_epi16(predY, ymin);
				predU=_mm256_max_epi16(predU, umin);
				predV=_mm256_max_epi16(predV, vmin);
				predY=_mm256_min_epi16(predY, ymax);
				predU=_mm256_min_epi16(predU, umax);
				predV=_mm256_min_epi16(predV, vmax);
			}
			__m256i msyms, moffset;
			if(fwd)
			{
				__m256i ctxblendmask=_mm256_set1_epi16(255);
				myuv[0]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+yidx)), half8));//load yuv
				myuv[1]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+uidx)), half8));
				myuv[2]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+vidx)), half8));

				//encode Y
				msyms=_mm256_sub_epi16(myuv[0], predY);//sub pred
				if(dist>1)
				{
					__m256i t0=_mm256_srai_epi16(msyms, 15);
					msyms=_mm256_mulhi_epi16(msyms, dist_rcp);
					msyms=_mm256_sub_epi16(msyms, t0);
					myuv[0]=_mm256_mullo_epi16(msyms, mdist);
					myuv[0]=_mm256_add_epi16(myuv[0], predY);
					myuv[0]=_mm256_max_epi16(myuv[0], amin);
					myuv[0]=_mm256_min_epi16(myuv[0], amax);
#ifdef ENABLE_GUIDE
					__m128i modified=_mm_packs_epi16(_mm256_extracti128_si256(myuv[0], 0), _mm256_extracti128_si256(myuv[0], 1));
					modified=_mm_xor_si128(modified, _mm_set1_epi8(-128));
					_mm_storeu_si128((__m128i*)(g_image+(imptr-interleaved-isize)+yidx), modified);
#endif
				}
				W[0]=myuv[0];
				_mm256_store_si256((__m256i*)rows[0]+0+0+0*6, myuv[0]);//store Y neighbors
#ifdef SAVE_RESIDUALS
				{
					ptrdiff_t idx=imptr-interleaved-isize+(ptrdiff_t)0*NCODERS;
					short syms2[16];
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
				if(dist>1)
				{
					__m256i t0=_mm256_srai_epi16(msyms, 15);
					msyms=_mm256_mulhi_epi16(msyms, dist_rcp);
					msyms=_mm256_sub_epi16(msyms, t0);
					myuv[1]=_mm256_mullo_epi16(msyms, mdist);
					myuv[1]=_mm256_add_epi16(myuv[1], predU);
					myuv[1]=_mm256_max_epi16(myuv[1], amin);
					myuv[1]=_mm256_min_epi16(myuv[1], amax);
#ifdef ENABLE_GUIDE
					__m128i modified=_mm_packs_epi16(_mm256_extracti128_si256(myuv[1], 0), _mm256_extracti128_si256(myuv[1], 1));
					modified=_mm_xor_si128(modified, _mm_set1_epi8(-128));
					_mm_storeu_si128((__m128i*)(g_image+(imptr-interleaved-isize)+uidx), modified);
#endif
				}
#ifdef SAVE_RESIDUALS
				{
					ptrdiff_t idx=imptr-interleaved-isize+(ptrdiff_t)1*NCODERS;
					short syms2[16];
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
				if(dist>1)
				{
					__m256i t0=_mm256_srai_epi16(msyms, 15);
					msyms=_mm256_mulhi_epi16(msyms, dist_rcp);
					msyms=_mm256_sub_epi16(msyms, t0);
					myuv[2]=_mm256_mullo_epi16(msyms, mdist);
					myuv[2]=_mm256_add_epi16(myuv[2], predV);
					myuv[2]=_mm256_max_epi16(myuv[2], amin);
					myuv[2]=_mm256_min_epi16(myuv[2], amax);
#ifdef ENABLE_GUIDE
					__m128i modified=_mm_packs_epi16(_mm256_extracti128_si256(myuv[2], 0), _mm256_extracti128_si256(myuv[2], 1));
					modified=_mm_xor_si128(modified, _mm_set1_epi8(-128));
					_mm_storeu_si128((__m128i*)(g_image+(imptr-interleaved-isize)+vidx), modified);
#endif
				}
#ifdef SAVE_RESIDUALS
				{
					ptrdiff_t idx=imptr-interleaved-isize+(ptrdiff_t)2*NCODERS;
					short syms2[16];
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
				_mm256_store_si256((__m256i*)ctxptr+2, ctxV);//store V  ctx|residuals		ctxptr+NCODERS*(C*2+R)
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
				ctxptr+=sizeof(short[3][NCODERS]);
			}
			else
			{
				//decode main
				__m128i msyms8;
				
				//decode Y
				if(userans)
				{
					dec_yuv(mstate, 0, &ctxY, (int*)CDF2syms, ans_permute, &streamptr, streamend, myuv+0);//residuals from [0 ~ 255]
					dec_yuv(mstate, 1, &ctxU, (int*)CDF2syms, ans_permute, &streamptr, streamend, myuv+1);
					dec_yuv(mstate, 2, &ctxV, (int*)CDF2syms, ans_permute, &streamptr, streamend, myuv+2);
				}
				else
				{
					ALIGN(32) uint16_t ctx[3*16], yuv[3*16];
					ctxU=_mm256_add_epi16(ctxU, mctxuoffset);
					ctxV=_mm256_add_epi16(ctxV, mctxvoffset);
					_mm256_store_si256((__m256i*)ctx+0, ctxY);
					_mm256_store_si256((__m256i*)ctx+1, ctxU);
					_mm256_store_si256((__m256i*)ctx+2, ctxV);

					yuv[0x00]=tANS_decodeSymbol(fsedstates+ctx[0x00], &ec);
					yuv[0x01]=tANS_decodeSymbol(fsedstates+ctx[0x01], &ec);
					yuv[0x02]=tANS_decodeSymbol(fsedstates+ctx[0x02], &ec);
					yuv[0x03]=tANS_decodeSymbol(fsedstates+ctx[0x03], &ec);
					yuv[0x04]=tANS_decodeSymbol(fsedstates+ctx[0x04], &ec);
					yuv[0x05]=tANS_decodeSymbol(fsedstates+ctx[0x05], &ec);
					yuv[0x06]=tANS_decodeSymbol(fsedstates+ctx[0x06], &ec);
					yuv[0x07]=tANS_decodeSymbol(fsedstates+ctx[0x07], &ec);
					yuv[0x08]=tANS_decodeSymbol(fsedstates+ctx[0x08], &ec);
					yuv[0x09]=tANS_decodeSymbol(fsedstates+ctx[0x09], &ec);
					yuv[0x0A]=tANS_decodeSymbol(fsedstates+ctx[0x0A], &ec);
					yuv[0x0B]=tANS_decodeSymbol(fsedstates+ctx[0x0B], &ec);
					yuv[0x0C]=tANS_decodeSymbol(fsedstates+ctx[0x0C], &ec);
					yuv[0x0D]=tANS_decodeSymbol(fsedstates+ctx[0x0D], &ec);
					yuv[0x0E]=tANS_decodeSymbol(fsedstates+ctx[0x0E], &ec);
					yuv[0x0F]=tANS_decodeSymbol(fsedstates+ctx[0x0F], &ec);
					yuv[0x10]=tANS_decodeSymbol(fsedstates+ctx[0x10], &ec);
					yuv[0x11]=tANS_decodeSymbol(fsedstates+ctx[0x11], &ec);
					yuv[0x12]=tANS_decodeSymbol(fsedstates+ctx[0x12], &ec);
					yuv[0x13]=tANS_decodeSymbol(fsedstates+ctx[0x13], &ec);
					yuv[0x14]=tANS_decodeSymbol(fsedstates+ctx[0x14], &ec);
					yuv[0x15]=tANS_decodeSymbol(fsedstates+ctx[0x15], &ec);
					yuv[0x16]=tANS_decodeSymbol(fsedstates+ctx[0x16], &ec);
					yuv[0x17]=tANS_decodeSymbol(fsedstates+ctx[0x17], &ec);
					yuv[0x18]=tANS_decodeSymbol(fsedstates+ctx[0x18], &ec);
					yuv[0x19]=tANS_decodeSymbol(fsedstates+ctx[0x19], &ec);
					yuv[0x1A]=tANS_decodeSymbol(fsedstates+ctx[0x1A], &ec);
					yuv[0x1B]=tANS_decodeSymbol(fsedstates+ctx[0x1B], &ec);
					yuv[0x1C]=tANS_decodeSymbol(fsedstates+ctx[0x1C], &ec);
					yuv[0x1D]=tANS_decodeSymbol(fsedstates+ctx[0x1D], &ec);
					yuv[0x1E]=tANS_decodeSymbol(fsedstates+ctx[0x1E], &ec);
					yuv[0x1F]=tANS_decodeSymbol(fsedstates+ctx[0x1F], &ec);
					yuv[0x20]=tANS_decodeSymbol(fsedstates+ctx[0x20], &ec);
					yuv[0x21]=tANS_decodeSymbol(fsedstates+ctx[0x21], &ec);
					yuv[0x22]=tANS_decodeSymbol(fsedstates+ctx[0x22], &ec);
					yuv[0x23]=tANS_decodeSymbol(fsedstates+ctx[0x23], &ec);
					yuv[0x24]=tANS_decodeSymbol(fsedstates+ctx[0x24], &ec);
					yuv[0x25]=tANS_decodeSymbol(fsedstates+ctx[0x25], &ec);
					yuv[0x26]=tANS_decodeSymbol(fsedstates+ctx[0x26], &ec);
					yuv[0x27]=tANS_decodeSymbol(fsedstates+ctx[0x27], &ec);
					yuv[0x28]=tANS_decodeSymbol(fsedstates+ctx[0x28], &ec);
					yuv[0x29]=tANS_decodeSymbol(fsedstates+ctx[0x29], &ec);
					yuv[0x2A]=tANS_decodeSymbol(fsedstates+ctx[0x2A], &ec);
					yuv[0x2B]=tANS_decodeSymbol(fsedstates+ctx[0x2B], &ec);
					yuv[0x2C]=tANS_decodeSymbol(fsedstates+ctx[0x2C], &ec);
					yuv[0x2D]=tANS_decodeSymbol(fsedstates+ctx[0x2D], &ec);
					yuv[0x2E]=tANS_decodeSymbol(fsedstates+ctx[0x2E], &ec);
					yuv[0x2F]=tANS_decodeSymbol(fsedstates+ctx[0x2F], &ec);

					myuv[0]=_mm256_load_si256((__m256i*)yuv+0);
					myuv[1]=_mm256_load_si256((__m256i*)yuv+1);
					myuv[2]=_mm256_load_si256((__m256i*)yuv+2);
				}

//#ifdef ANS_VAL
//				ALIGN(32) unsigned char debugvals[6*NCODERS];
//				msyms8=_mm256_packus_epi16(ctxY0, ctxY1);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_store_si256((__m256i*)debugvals+0, msyms8);
//				msyms8=_mm256_packus_epi16(myuv[0], myuv[1]);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_store_si256((__m256i*)debugvals+1, msyms8);
//				ansval_check(debugvals+0*NCODERS, 1, NCODERS);//contexts
//				ansval_check(debugvals+1*NCODERS, 1, NCODERS);//residuals
//#endif

				//yuv = (char)(sym+pred-128)	= (unsigned char)(sym+pred)-128
				if(dist>1)
				{
					msyms=_mm256_add_epi16(myuv[0], amin);
					myuv[0]=_mm256_mullo_epi16(msyms, mdist);
					myuv[0]=_mm256_add_epi16(myuv[0], predY);
					myuv[0]=_mm256_max_epi16(myuv[0], amin);
					myuv[0]=_mm256_min_epi16(myuv[0], amax);
				}
				else
				{
					myuv[0]=_mm256_add_epi16(myuv[0], predY);
					myuv[0]=_mm256_and_si256(myuv[0], bytemask);
					myuv[0]=_mm256_add_epi16(myuv[0], amin);
					msyms=_mm256_sub_epi16(myuv[0], predY);//sub pred
				}
				ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
				msyms8=_mm_packs_epi16(_mm256_extracti128_si256(myuv[0], 0), _mm256_extracti128_si256(myuv[0], 1));
				msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
				_mm_store_si128((__m128i*)(imptr+yidx), msyms8);//store Y bytes
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+yidx, g_image+(imptr-interleaved)+yidx, NCODERS))
				{
					printf("original  decoded  original-decoded  XYC0 %d %d %d\n", kx, ky, yidx);
					for(int k=0;k<NCODERS;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+yidx+k],
							imptr[yidx+k],
							g_image[imptr-interleaved+yidx+k]-imptr[yidx+k]
						);
					CRASH("guide error XYC0 %d %d %d/%d", kx, ky, yidx, NCODERS);
				}
#endif
//#ifdef ANS_VAL
//				ALIGN(32) unsigned char actx[3*NCODERS];
//				msyms8=_mm256_packus_epi16(ctxY0, ctxY1);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_storeu_si256((__m256i*)actx+0, msyms8);
//				ansval_check(actx+0*NCODERS, 1, NCODERS);
//				ansval_check(imptr+kx+yidx, 1, NCODERS);
//#endif
				W[0]=myuv[0];
				_mm256_store_si256((__m256i*)rows[0]+0+0+0*6, myuv[0]);//store Y neighbors


				//decode U
				moffset=_mm256_and_si256(myuv[0], uhelpmask);
				predU=_mm256_add_epi16(predU, moffset);
				predU=_mm256_max_epi16(predU, amin);
				predU=_mm256_min_epi16(predU, amax);
//#ifdef ANS_VAL
//				msyms8=_mm256_packus_epi16(ctxU0, ctxU1);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_store_si256((__m256i*)debugvals+2, msyms8);
//				msyms8=_mm256_packus_epi16(myuv[2], myuv[3]);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_store_si256((__m256i*)debugvals+3, msyms8);
//				ansval_check(debugvals+2*NCODERS, 1, NCODERS);//contexts
//				ansval_check(debugvals+3*NCODERS, 1, NCODERS);//residuals
//#endif
				
				if(dist>1)
				{
					msyms=_mm256_add_epi16(myuv[1], amin);
					myuv[1]=_mm256_mullo_epi16(msyms, mdist);
					myuv[1]=_mm256_add_epi16(myuv[1], predU);
					myuv[1]=_mm256_max_epi16(myuv[1], amin);
					myuv[1]=_mm256_min_epi16(myuv[1], amax);
				}
				else
				{
					myuv[1]=_mm256_add_epi16(myuv[1], predU);
					myuv[1]=_mm256_and_si256(myuv[1], bytemask);
					myuv[1]=_mm256_add_epi16(myuv[1], amin);
					msyms=_mm256_sub_epi16(myuv[1], predU);//sub pred
				}
				ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
				msyms8=_mm_packs_epi16(_mm256_extracti128_si256(myuv[1], 0), _mm256_extracti128_si256(myuv[1], 1));
				msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
				_mm_store_si128((__m128i*)(imptr+uidx), msyms8);//store U bytes
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+uidx, g_image+(imptr-interleaved)+uidx, NCODERS))
				{
					printf("original  decoded  original-decoded  XYC1 %d %d %d\n", kx, ky, uidx);
					for(int k=0;k<NCODERS;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+uidx+k],
							imptr[uidx+k],
							g_image[imptr-interleaved+uidx+k]-imptr[uidx+k]
						);
					CRASH("guide error XYC1 %d %d %d/%d", kx, ky, uidx, NCODERS);
				}
#endif
//#ifdef ANS_VAL
//				msyms8=_mm256_packus_epi16(ctxU0, ctxU1);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_storeu_si256((__m256i*)actx+1, msyms8);
//				ansval_check(actx+1*NCODERS, 1, NCODERS);
//				ansval_check(imptr+kx+uidx, 1, NCODERS);
//#endif
				W[1]=_mm256_sub_epi16(myuv[1], moffset);//subtract Uoffset from U
				_mm256_store_si256((__m256i*)rows[0]+0+1+0*6, W[1]);//store U neighbors
				

				//decode V
				moffset=_mm256_mullo_epi16(vc0, myuv[0]);
				moffset=_mm256_add_epi16(moffset, _mm256_mullo_epi16(vc1, myuv[1]));
				moffset=_mm256_srai_epi16(moffset, 2);
				predV=_mm256_add_epi16(predV, moffset);
				predV=_mm256_max_epi16(predV, amin);
				predV=_mm256_min_epi16(predV, amax);
				
				if(dist>1)
				{
					msyms=_mm256_add_epi16(myuv[2], amin);
					myuv[2]=_mm256_mullo_epi16(msyms, mdist);
					myuv[2]=_mm256_add_epi16(myuv[2], predV);
					myuv[2]=_mm256_max_epi16(myuv[2], amin);
					myuv[2]=_mm256_min_epi16(myuv[2], amax);
				}
				else
				{
					myuv[2]=_mm256_add_epi16(myuv[2], predV);
					myuv[2]=_mm256_and_si256(myuv[2], bytemask);
					myuv[2]=_mm256_add_epi16(myuv[2], amin);
					msyms=_mm256_sub_epi16(myuv[2], predV);
				}
				ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
				msyms8=_mm_packs_epi16(_mm256_extracti128_si256(myuv[2], 0), _mm256_extracti128_si256(myuv[2], 1));
				msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
				_mm_store_si128((__m128i*)(imptr+vidx), msyms8);//store V bytes
#ifdef ENABLE_GUIDE
				if(memcmp(imptr+vidx, g_image+(imptr-interleaved)+vidx, NCODERS))
				{
					printf("original  decoded  original-decoded  XYC2 %d %d %d\n", kx, ky, vidx);
					for(int k=0;k<NCODERS;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+vidx+k],
							imptr[vidx+k],
							g_image[imptr-interleaved+vidx+k]-imptr[vidx+k]
						);
					CRASH("guide error XYC2 %d %d %d/%d", kx, ky, vidx, NCODERS);
				}
#endif
				W[2]=_mm256_sub_epi16(myuv[2], moffset);//subtract Voffset from V
				_mm256_store_si256((__m256i*)rows[0]+0+2+0*6, W[2]);//store V neighbors
			}
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
			rows[0]+=6*NCODERS;
			rows[1]+=6*NCODERS;
			rows[2]+=6*NCODERS;
			rows[3]+=6*NCODERS;
			imptr+=3*NCODERS;
		}
		//printf("%8d/%8d\r", ky+1, blockh);
	}
	prof_checkpoint(isize, "main");
#ifdef ENABLE_GUIDE
	if(dist>1&&!fwd)
	{
		double rmse[]=
		{
			sqrt(g_sqe[0]*3/isize),
			sqrt(g_sqe[1]*3/isize),
			sqrt(g_sqe[2]*3/isize),
			sqrt((g_sqe[0]+g_sqe[1]+g_sqe[2])/isize),
		};
		double psnr[]=
		{
			20*log10(255/rmse[0]),
			20*log10(255/rmse[1]),
			20*log10(255/rmse[2]),
			20*log10(255/rmse[3]),
		};
		printf("RMSE PSNR\n");
		printf("T %12.6lf %12.6lf\n", rmse[3], psnr[3]);
		printf("Y %12.6lf %12.6lf\n", rmse[0], psnr[0]);
		printf("U %12.6lf %12.6lf\n", rmse[1], psnr[1]);
		printf("V %12.6lf %12.6lf\n", rmse[2], psnr[2]);
	}
#endif

#ifdef PRINT_L1_BOUNDS
	printf("Coeff %8d ~ %8d\n", cmin, cmax);
	printf("Bias  %8d ~ %8d\n", bmin, bmax);
#endif
	if(effort)
		_mm_free(L1state);
	if(fwd)
	{
		//All ANS encoding here is  bwd image - bwd stream

		//(original FSE encoding was  bwd image - fwd stream)
		rANS_SIMD_SymInfo *syminfo=0;
		int nhistsize=nctx*(int)sizeof(int16_t[256]);
		int16_t *nhist=0;
		uint8_t probbits[3*NCTX+3]={0};

		if(userans)//normalize/integrate hists
		{
			syminfo=(rANS_SIMD_SymInfo*)CDF2syms;
			for(int kc=0;kc<nctx;++kc)
				enc_hist2stats(hists+(ptrdiff_t)256*kc, syminfo+(ptrdiff_t)256*kc, &bypassmask, kc);
		}
		else
		{
#ifdef ESTIMATE_SIZE2
			double esize=0;
#endif
			nhist=(int16_t*)malloc(nhistsize);
			if(!nhist)
			{
				CRASH("Alloc error");
				return 1;
			}
			memset(nhist, 0, nhistsize);
			for(int kc=0;kc<nctx;++kc)
			{
				int16_t *currnhist=nhist+kc*256;
				uint32_t *hcurr=(uint32_t*)hists+256*kc;
				int sum=0, nlevels=0, vmax=0, fmax=0;
				for(int ks=0;ks<256;++ks)
				{
					int freq=hcurr[ks];
					sum+=freq;
					nlevels+=freq!=0;
					if(fmax<freq)
					{
						fmax=freq;
						vmax=ks;
					}
				}
				if(nlevels==1)//degenerate distribution
				{
					currnhist[0]=vmax;
					fseerror=tANS_buildCTable_rle(fsetables[kc], vmax);
					if(fseerror)
					{
						CRASH("FSE Error");
						return 1;
					}
					degenmask|=1ULL<<kc;
				}
				else if(sum>256)//ordinary distribution
				{
					normalizehist(hcurr, (uint16_t*)currnhist);
					probbits[kc]=PROBBITS;
					//fseerror=tANS_normalizeCount((short*)currnhist, PROBBITS, hcurr, sum, 255);
					//if(fseerror)
					//{
					//	CRASH("Normalize failed");
					//	return 1;
					//}
					//probbits[kc]=(uint8_t)fseerror;
					fseerror=tANS_buildCTable(fsetables[kc], (short*)currnhist, 255, probbits[kc]);
					if(fseerror)
					{
						CRASH("FSE Error");
						return 1;
					}
				}
				else//bypass distribution
				{
					fseerror=tANS_buildCTable_raw(fsetables[kc], 8);
					if(fseerror)
					{
						CRASH("FSE Error");
						return 1;
					}
					bypassmask|=1ULL<<kc;
				}
#ifdef ESTIMATE_SIZE2
				double e=0, invsum=1./sum;
				for(int ks=0;ks<256;++ks)
				{
					int freq=hcurr[ks];
					if(freq)
						e-=freq*log2(freq*invsum);
				}
				e/=8;
				esize+=e;
				printf("%c  %2d  hsum=%12d  %12.2lf bytes  %s\n"
					, "YUV"[kc/NCTX]
					, kc%NCTX
					, sum
					, e
					, degenmask>>kc&1?"-----DEGEN":(bypassmask>>kc&1?"BYPASS----":"--normal--")
				);
#endif
			}
#ifdef ESTIMATE_SIZE2
			printf("\n%12.2lf\n\n", esize);
#endif
			bitpacker_enc_init(&ec, streamstart, streamptr);
			for(int kc=0;kc<nctx;++kc)
				tANS_initCState(fsecstates+kc, fsetables[kc]);
		}
		
		//encode remainder
		if(xremw||yremh)
		{
			int32_t *rhist=hists+3*NCTX*256;
			for(int ky=0;ky<yremh;++ky)
				decorr1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, bestrct, rhist);
			for(int kx=0;kx<xremw;++kx)
				decorr1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, bestrct, rhist);
			if(userans)
			{
				rANS_SIMD_SymInfo *rsyminfo=(rANS_SIMD_SymInfo*)CDF2syms+3*NCTX*256;

				unsigned state=1<<(RANS_STATE_BITS-RANS_RENORM_BITS);
				for(int kx=xremw-1;kx>=0;--kx)
					encode1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, &state, &streamptr, image, rsyminfo);
				for(int ky=yremh-1;ky>=0;--ky)
					encode1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, &state, &streamptr, image, rsyminfo);
				//flush
				streamptr-=4;
#ifdef _DEBUG
				if(streamptr<=image)
					CRASH("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
				*(unsigned*)streamptr=state;
			}
			else
			{
				for(int kx=xremw-1;kx>=0;--kx)
					encode1d_fse(image+qxbytes+3*kx, blockh*YCODERS, rowstride, &ec, fsecstates+3*NCTX);
				for(int ky=yremh-1;ky>=0;--ky)
					encode1d_fse(image+rowstride*(blockh*YCODERS+ky), iw, 3, &ec, fsecstates+3*NCTX);
				streamptr=ec.dstbwdptr;
			}
			prof_checkpoint(usize-isize, "encode remainder");
			profile_size(streamptr, "/ %9td bytes remainders", usize-isize);
		}

		//encode main
		unsigned short *ctxptr2=(unsigned short*)(interleaved+(isize<<1));
		if(userans)
		{
			mstate[1]=mstate[0]=_mm256_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
			for(int ky=blockh-1;ky>=0;--ky)
			{
#ifdef ESTIMATE_SIZE
				int kc=2;
#endif
				for(int kx=ixbytes-NCODERS;kx>=0;kx-=NCODERS)//ixbytes = iw/XCODERS*NCODERS*3
				{
					__m256i mmax[2], minvf[2], mcdf[2], mnegf_sh[2];
					{
						__m256i s0, s1, s2, s3;
						__m256i t0, t1, t2, t3;
#define SHUFFLE_PS(LO, HI, IMM8_HHLL) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(LO), _mm256_castsi256_ps(HI), IMM8_HHLL))
						
						ctxptr2-=NCODERS;

						s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+0])));
						s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+1])));
						s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+2])));
						s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+3])));
						s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+4])), 1);
						s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+5])), 1);
						s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+6])), 1);
						s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[0*8+7])), 1);
						t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));//_MM_TRANSPOSE4_PS
						t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
						t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
						t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
						mmax	[0]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
						minvf	[0]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
						mcdf	[0]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
						mnegf_sh[0]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));

						s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+0])));
						s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+1])));
						s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+2])));
						s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+3])));
						s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+4])), 1);
						s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+5])), 1);
						s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+6])), 1);
						s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[1*8+7])), 1);
						t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));
						t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
						t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
						t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
						mmax	[1]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
						minvf	[1]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
						mcdf	[1]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
						mnegf_sh[1]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));
					}
#ifdef ESTIMATE_SIZE
					{
						ALIGN(32) int anegf[NCODERS]={0};
						memcpy(anegf, mnegf_sh, sizeof(anegf));
						const double norm=1./(1<<PROBBITS);
						for(int k=0;k<NCODERS;++k)
						{
							int freq=(1<<PROBBITS)-(anegf[k]&0xFFFF);
							if((unsigned)(freq-1)>=(unsigned)((1<<PROBBITS)-1))
								CRASH("freq = %d", freq);
							esize[kc*NCODERS+k]-=log2(freq*norm)*0.125;
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
					__m256i idx0=_mm256_load_si256((__m256i*)ans_permute+mask0);
					__m256i idx1=_mm256_load_si256((__m256i*)ans_permute+mask1);
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
					if(streamptr-(2*((ptrdiff_t)mask0+mask1)+sizeof(__m128i))<=image)
						CRASH("OOB ptr %016zX <= %016zX", streamptr, image);
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
					ansval_push(mstate, sizeof(int), NCODERS);
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
						ALIGN(32) unsigned short freqs[NCODERS];
						_mm256_store_si256((__m256i*)freqs, negf0);
						ansval_push(freqs, sizeof(short), NCODERS);
#endif
					}
					mstate[0]=_mm256_add_epi32(mstate[0], minvf[0]);
					mstate[1]=_mm256_add_epi32(mstate[1], minvf[1]);
#ifdef ANS_VAL
					ansval_push(mstate, sizeof(int), NCODERS);
#endif
				}
				//printf("%8d/%8d\r", blockh-1-ky, blockh);
			}
			//flush
			streamptr-=sizeof(mstate);
#ifdef _DEBUG
			if(streamptr<=image)
				CRASH("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
			memcpy(streamptr, mstate, sizeof(mstate));
		}
		else
		{
			int ctr=0;
			while((ptrdiff_t)ctxptr2>(ptrdiff_t)interleaved)
			{
				int sym;

				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
#if defined _DEBUG && 0
				if(ctr==blockw*blockh-1)//
					printf("");
#endif
				sym=*--ctxptr2; tANS_encodeSymbol(&ec, fsecstates+(sym>>8), (uint8_t)sym);
				++ctr;
			}
			if(ctr!=blockw*blockh)
				CRASH("");
			for(int kc=nctx-1;kc>=0;--kc)
			{
#if defined _DEBUG && 0
				printf("flush 0x%04tX\n", fsecstates[kc].value);
#endif
				tANS_flushCState(&ec, fsecstates+kc);
			}
			//fseerror=BIT_closeCStream(&fsecstream);
			//if(!fseerror)
			//{
			//	CRASH("INFLATION: Stream did not fit into allocated buffer");
			//	return 1;
			//}
			//streamptr+=fseerror;
			streamptr=ec.dstbwdptr;
		}
		prof_checkpoint(isize, "encode main");
		profile_size(streamptr, "/ %9td bytes main", isize);

		//pack hists
		if(userans)
		{
			//BitPackerLIFO ec;
			bitpacker_enc_init(&ec, image, streamptr);
			for(int kc=nctx-1;kc>=0;--kc)
				enc_packhist(&ec, hists+(ptrdiff_t)256*kc, bypassmask, kc, userans);
			bitpacker_enc_flush(&ec);
			streamptr=ec.dstbwdptr;
		}
		else
		{
			for(int kc=nctx-1;kc>=0;--kc)
			{
				//if(!kc)
				//	printf("");
				int16_t *currnhist=nhist+kc*256;
				if(degenmask>>kc&1)//degenerate distribution
					bitpacker_enc(&ec, 8, (uint8_t)currnhist[0]);
					//*streamptr++=(uint8_t)currnhist[0];
				else if(!(bypassmask>>kc&1))//ordinary distribution
					enc_packhist(&ec, (int*)currnhist, bypassmask, kc, userans);
				//{
				//	fseerror=tANS_writeNCount(streamptr, streamend-streamptr, (short*)currnhist, 255, probbits[kc]);
				//	if(fseerror)
				//	{
				//		CRASH("FSE Error");
				//		return 1;
				//	}
				//}
#if defined _DEBUG && 0
				if(kc==0)
				{
					printf("enc.ctx = 0:\n");
					for(int ks=0;ks<256;++ks)
						printf("%3d  %04X\n", ks, currnhist[ks]);
					printf("\n");
				}
#endif
			}
			bitpacker_enc_flush(&ec);
			streamptr=ec.dstbwdptr;
			free(nhist);
		}
		prof_checkpoint((ptrdiff_t)nctx*256, "pack histograms");
		profile_size(streamptr, "/ %9d bytes overhead", ((nctx*12<<8)+7)>>3);

		//save compressed file
		{
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				CRASH("Cannot open \"%s\" for writing", fdst);
				return 1;
			}
			ptrdiff_t csize2=0;
			csize2+=fwrite("F2", 1, 2, fdst);
			csize2+=fwrite(&iw, 1, 3, fdst);
			csize2+=fwrite(&ih, 1, 3, fdst);
			int flags=bestrct<<3|userans<<2|effort;
			csize2+=fwrite(&flags, 1, 2, fdst);
			csize2+=fwrite(&dist, 1, 1, fdst);
			if(!userans)
				csize2+=fwrite(&degenmask, 1, 8, fdst);
			csize2+=fwrite(&bypassmask, 1, 8, fdst);
#ifdef _DEBUG
			if(streamptr>streamend)
				CRASH("OOB ptr %016zX > %016zX", streamptr, streamend);
			if(streamptr<streamstart)
				CRASH("OOB ptr %016zX < %016zX", streamptr, image);
#endif
			csize2+=fwrite(streamptr, 1, streamend-streamptr, fdst);
			fclose(fdst);
			
#ifdef ESTIMATE_SIZE
			double etotal=0;
			for(int k=0;k<3*NCODERS;++k)
			{
				etotal+=esize[k];
				printf("E %12.2lf\n", esize[k]);//
			}
			printf("Total estimate  %12.2lf bytes\n", etotal);
#endif
#ifdef LOUD
			printf("%s  WH %d*%d\n", srcfn, iw, ih);
			printf("%8td/%8td bytes\n", csize2, usize2);
#endif
			(void)csize2;
			prof_checkpoint(csize2, "fwrite");
		}
		free(hists);
#ifdef SAVE_RESIDUALS
		{
			unsigned char *result=(unsigned char*)malloc(usize);
			if(!result)
			{
				CRASH("Alloc error");
				return 1;
			}
			memset(result, -128, usize);
			interleave_blocks_inv(residuals, iw, ih, result);
			{
				const char fn[]="20250603_0628PM.PPM";
				FILE *fdst2=fopen(fn, "wb");
				if(!fdst2)
				{
					CRASH("Cannot open \"%s\" for writing", fn);
					return 1;
				}
				fprintf(fdst2, "P6\n%d %d\n255\n", iw, ih);
				fwrite(result, 1, usize, fdst2);
				fclose(fdst2);
			}
			free(residuals);
			free(result);
		}
#endif
	}
	else
	{
		//deinterleave
		interleave_blocks_inv(interleaved, iw, ih, image);
		prof_checkpoint(usize, "deinterleave");

		if(xremw||yremh)
		{
#ifdef _DEBUG
			if(streamptr>streamend)
				CRASH("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			if(userans)
			{
				uint32_t *rCDF2syms=CDF2syms+(3*NCTX<<PROBBITS);
				unsigned state=*(unsigned*)streamptr;
				streamptr+=4;
				for(int ky=0;ky<yremh;++ky)
					decode1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, bestrct, &state, (const unsigned char**)&streamptr, streamend, rCDF2syms);
				for(int kx=0;kx<xremw;++kx)
					decode1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, bestrct, &state, (const unsigned char**)&streamptr, streamend, rCDF2syms);
			}
			else
			{
				for(int ky=0;ky<yremh;++ky)
					decode1d_fse(image+rowstride*(blockh*YCODERS+ky), iw, 3, bestrct, &ec, fsedstates+3*NCTX);
				for(int kx=0;kx<xremw;++kx)
					decode1d_fse(image+qxbytes+3*kx, blockh*YCODERS, rowstride, bestrct, &ec, fsedstates+3*NCTX);
			}
			prof_checkpoint(usize-isize, "remainder");
		}

		//save PPM file
		save_ppm(dstfn, image, iw, ih);
		prof_checkpoint(usize, "fwrite");
	}
	_mm_free(pixels);
	_mm_free(interleaved);
	free(image);
	if(userans)
		_mm_free(CDF2syms);
	else
	{
		for(int k=0;k<nctx;++k)
			free(fsetables[k]);
	}

#ifdef LOUD
	t=time_sec()-t;
	printf("%c  %12.6lf sec  %12.6lf MB/s\n", 'D'+fwd, t, usize/(t*1024*1024));
#endif
#ifdef PROFILE_TIME
#ifdef __GNUC__
	if(profile)
#endif
		prof_print(usize);
#endif
	(void)och_names;
	(void)rct_names;
	return 0;
}
