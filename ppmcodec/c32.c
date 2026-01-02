#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


	#define PROFILE_TIME		//should be on

#ifdef _MSC_VER
	#define PROFILE_SIZE
	#define LOUD			//size & time

//	#define ESTIMATE_SIZE		//DEBUG		checks for zero frequency, visualizes context usage
	#define ENABLE_GUIDE		//DEBUG		checks interleaved pixels
	#define ANS_VAL			//DEBUG

//	#define PRINT_L1_BOUNDS
//	#define WG4_PRINTMAXERR
//	#define WG4_SERIALDEBUG
//	#define TEST_INTERLEAVE
#endif

	#define USE_YCBCR
//	#define ONE_LOSSYMEAN
//	#define ONE_LOSSYCOV
	#define ENABLE_RCT_EXTENSION
//	#define EMULATE_GATHER		//gather is a little faster
	#define VNATIVE			//disables useless v-chroma scaling by 4
//	#define NAIVEMIX		//3.5 3.8x slower	55 58 MB/s vs 196 224 MB/s
//	#define WG_COMMONMIX		//bad
	#define WG_ENABLE_eW		//eW is bad with blocks unlike in eBench
	#define INTERLEAVESIMD		//2.5x faster interleave


//#define NEAR_DISTORTION 7		//d >= 4

//3*17+3=54 contexts
#define GRBITS 3
#define NCTX 18		//18*3+3 = 57 total

#define XCODERS 4	//xrem 1~3 cols
#define YCODERS 4	//yrem 1~3 rows

#define NCODERS 16

#define PROBBITS 12	//12 bit max	James Bonfield's CDF2sym: {freq<<20 | bias<<8 | sym}

#define RANS_STATE_BITS 31
#define RANS_RENORM_BITS 16

#include"entropy.h"
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
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(const unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		LOG_ERROR("");
		printf("");
	}
}
//static void guide_update_error(const unsigned char *image, int kx, int ky)
//{
//	int idx=3*(g_iw*ky+kx);
//	double diff;
//	diff=g_image[idx+0]-image[idx+0], g_sqe[0]+=diff*diff;
//	diff=g_image[idx+1]-image[idx+1], g_sqe[1]+=diff*diff;
//	diff=g_image[idx+2]-image[idx+2], g_sqe[2]+=diff*diff;
//}
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
			LOG_ERROR("Profiler OOB");
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
	//printf("| ");
	int colors[128]={0};
	srand((unsigned)__rdtsc());
	colorgen(colors, prof_count, 64, 300, 100);
	//colorgen0(colors, prof_count, 0xC0C0C0);
	printf("1 char = 4 ms\n");
	printf("|");
	for(int k=0;k<prof_count;++k)
	{
		SpeedProfilerInfo *info=prof_data+k;
		csum+=info->t;
		int curr=(int)(csum*250);//fixed scale
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
			colorprintf(colors[k], colors[k], buf);
			colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%s", info->msg);
			memset(buf, '-', (ptrdiff_t)space-labelend);
			buf[space-labelend]=0;
			colorprintf(colors[k], colors[k], buf);
		//	colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%*s%s%*s", labelstart, "", info->msg, space-labelend, "");
		}
		else
		{
			memset(buf, '-', space);
			buf[space]=0;
			colorprintf(colors[k], colors[k], buf);
		}
		//	colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%*s", space, "");
		printf("|");
		prev=curr;
	}
	printf("\n");
	for(int k=0;k<prof_count;++k)
	{
		SpeedProfilerInfo *info=prof_data+k;
		//int nstars=(int)(info->t/tmax*64+0.5);
	//	colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "  %c", k+'A');
		printf("%16.7lf ms %8.4lf%% ", info->t*1000, 100.*info->t/timesum);
	//	printf("  %c: %16.12lf sec %8.4lf%% ", k+'A', info->t, info->t/timesum);
		if(info->size)
			printf(" %16.6lf MB/s %10td bytes ", info->size/(info->t*1024*1024), info->size);
		if(info->msg)
			colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%-20s", info->msg);
		else// if(nstars)
			colorprintf(COLORPRINTF_TXT_DEFAULT, colors[k], "%-20s", "");
		//for(int k2=0;k2<nstars;++k2)
		//	printf("*");
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
//	II_COEFF_U_SUB_V_NBLI,
//	II_COEFF_V_SUB_U_NBLI,

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


#define WG_NPREDS 8

#ifdef WG4_PRINTMAXERR
static unsigned maxerror[16*8]={0};
#endif
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
AWM_INLINE void wg_mix_pt2(const __m256i *wgpreds, const __m256i *prederrors, __m256i *result)//process halfwave
{
	//mix(const short wgpreds[8, stride 6][16], unsigned short prederrors[8, stride 2][16]) -> float[8][16]
	__m256i coeff[16];

	//convert 16 shorts * 8 regs -> 8 floats * 16 regs		each 2 float regs make half a wave
	__m256i wordmask=_mm256_set1_epi32(0xFFFF);
	coeff[0*2+0]=_mm256_and_si256(prederrors[0], wordmask);
	coeff[1*2+0]=_mm256_and_si256(prederrors[1], wordmask);
	coeff[2*2+0]=_mm256_and_si256(prederrors[2], wordmask);
	coeff[3*2+0]=_mm256_and_si256(prederrors[3], wordmask);
	coeff[4*2+0]=_mm256_and_si256(prederrors[4], wordmask);
	coeff[5*2+0]=_mm256_and_si256(prederrors[5], wordmask);
	coeff[6*2+0]=_mm256_and_si256(prederrors[6], wordmask);
	coeff[7*2+0]=_mm256_and_si256(prederrors[7], wordmask);
	coeff[0*2+1]=_mm256_srli_epi32(prederrors[0], 16);
	coeff[1*2+1]=_mm256_srli_epi32(prederrors[1], 16);
	coeff[2*2+1]=_mm256_srli_epi32(prederrors[2], 16);
	coeff[3*2+1]=_mm256_srli_epi32(prederrors[3], 16);
	coeff[4*2+1]=_mm256_srli_epi32(prederrors[4], 16);
	coeff[5*2+1]=_mm256_srli_epi32(prederrors[5], 16);
	coeff[6*2+1]=_mm256_srli_epi32(prederrors[6], 16);
	coeff[7*2+1]=_mm256_srli_epi32(prederrors[7], 16);
	
	//convert 16 shorts * 8 regs -> 8 floats * 16 regs
#if 1
	__m256i estim[16];
	estim[0*2+0]=_mm256_slli_epi32(wgpreds[0*3+0], 16);
	estim[1*2+0]=_mm256_slli_epi32(wgpreds[1*3+0], 16);
	estim[2*2+0]=_mm256_slli_epi32(wgpreds[2*3+0], 16);
	estim[3*2+0]=_mm256_slli_epi32(wgpreds[3*3+0], 16);
	estim[4*2+0]=_mm256_slli_epi32(wgpreds[4*3+0], 16);
	estim[5*2+0]=_mm256_slli_epi32(wgpreds[5*3+0], 16);
	estim[6*2+0]=_mm256_slli_epi32(wgpreds[6*3+0], 16);
	estim[7*2+0]=_mm256_slli_epi32(wgpreds[7*3+0], 16);
	estim[0*2+0]=_mm256_srai_epi32(estim[0*2+0], 16);//preds are signed
	estim[1*2+0]=_mm256_srai_epi32(estim[1*2+0], 16);
	estim[2*2+0]=_mm256_srai_epi32(estim[2*2+0], 16);
	estim[3*2+0]=_mm256_srai_epi32(estim[3*2+0], 16);
	estim[4*2+0]=_mm256_srai_epi32(estim[4*2+0], 16);
	estim[5*2+0]=_mm256_srai_epi32(estim[5*2+0], 16);
	estim[6*2+0]=_mm256_srai_epi32(estim[6*2+0], 16);
	estim[7*2+0]=_mm256_srai_epi32(estim[7*2+0], 16);
	estim[0*2+1]=_mm256_srai_epi32(wgpreds[0*3+0], 16);
	estim[1*2+1]=_mm256_srai_epi32(wgpreds[1*3+0], 16);
	estim[2*2+1]=_mm256_srai_epi32(wgpreds[2*3+0], 16);
	estim[3*2+1]=_mm256_srai_epi32(wgpreds[3*3+0], 16);
	estim[4*2+1]=_mm256_srai_epi32(wgpreds[4*3+0], 16);
	estim[5*2+1]=_mm256_srai_epi32(wgpreds[5*3+0], 16);
	estim[6*2+1]=_mm256_srai_epi32(wgpreds[6*3+0], 16);
	estim[7*2+1]=_mm256_srai_epi32(wgpreds[7*3+0], 16);
#endif
#ifdef WG4_PRINTMAXERR
	for(int k=0;k<16*8;++k)
	{
		unsigned e=((unsigned*)coeff)[k];
		if(maxerror[k]<e)
			maxerror[k]=e;
	}
#endif
	
	/*
	|	WG5	DIV2K -0.13%
	|		GDCC -0.08%
	0	210	240 ^	N
	1	210	240 ^	W
	2	215	180 v	3*(N-NN)+NNN
	3	215	180 v	3*(W-WW)+WWW		4*(W+WWW)-6*WW-WWWW
	4	140	140	W+NE-N
	5	230	160 v	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4
	6	120	120	N+W-NW
	7	140	120 v	N+NE-NNE
	*/
#define NUMBITS 15
#define DENBITS 7
#define DIVLUTSIZE (1<<DENBITS)
	static int divlookup[WG_NPREDS+1][DIVLUTSIZE]={0};
	if(!*(int*)divlookup)
	{
		static const int weights[WG_NPREDS]=
		{
			240,
			240,
			180,
			180,
			160,
			160,
			120,
			120,
		};
		for(int kp=0;kp<WG_NPREDS;++kp)
		{
			for(int k=0;k<DIVLUTSIZE;++k)
				divlookup[kp][k]=(weights[kp]<<NUMBITS)/(k+1);
		}
		for(int k=0;k<DIVLUTSIZE;++k)
			divlookup[WG_NPREDS][k]=(1<<NUMBITS)/(k+1);
	}
	__m256i sh[16];
	__m256i one=_mm256_set1_epi32(1);
	sh[0x0]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0x0], one));
	sh[0x1]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0x1], one));
	sh[0x2]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0x2], one));
	sh[0x3]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0x3], one));
	sh[0x4]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0x4], one));
	sh[0x5]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0x5], one));
	sh[0x6]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0x6], one));
	sh[0x7]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0x7], one));
	sh[0x8]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0x8], one));
	sh[0x9]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0x9], one));
	sh[0xA]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0xA], one));
	sh[0xB]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0xB], one));
	sh[0xC]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0xC], one));
	sh[0xD]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0xD], one));
	sh[0xE]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0xE], one));
	sh[0xF]=FLOOR_LOG2_32x8(_mm256_add_epi32(coeff[0xF], one));
	__m256i offset=_mm256_set1_epi32(DENBITS-1);
	sh[0x0]=_mm256_sub_epi32(sh[0x0], offset);
	sh[0x1]=_mm256_sub_epi32(sh[0x1], offset);
	sh[0x2]=_mm256_sub_epi32(sh[0x2], offset);
	sh[0x3]=_mm256_sub_epi32(sh[0x3], offset);
	sh[0x4]=_mm256_sub_epi32(sh[0x4], offset);
	sh[0x5]=_mm256_sub_epi32(sh[0x5], offset);
	sh[0x6]=_mm256_sub_epi32(sh[0x6], offset);
	sh[0x7]=_mm256_sub_epi32(sh[0x7], offset);
	sh[0x8]=_mm256_sub_epi32(sh[0x8], offset);
	sh[0x9]=_mm256_sub_epi32(sh[0x9], offset);
	sh[0xA]=_mm256_sub_epi32(sh[0xA], offset);
	sh[0xB]=_mm256_sub_epi32(sh[0xB], offset);
	sh[0xC]=_mm256_sub_epi32(sh[0xC], offset);
	sh[0xD]=_mm256_sub_epi32(sh[0xD], offset);
	sh[0xE]=_mm256_sub_epi32(sh[0xE], offset);
	sh[0xF]=_mm256_sub_epi32(sh[0xF], offset);
	offset=_mm256_setzero_si256();
	sh[0x0]=_mm256_max_epi32(sh[0x0], offset);
	sh[0x1]=_mm256_max_epi32(sh[0x1], offset);
	sh[0x2]=_mm256_max_epi32(sh[0x2], offset);
	sh[0x3]=_mm256_max_epi32(sh[0x3], offset);
	sh[0x4]=_mm256_max_epi32(sh[0x4], offset);
	sh[0x5]=_mm256_max_epi32(sh[0x5], offset);
	sh[0x6]=_mm256_max_epi32(sh[0x6], offset);
	sh[0x7]=_mm256_max_epi32(sh[0x7], offset);
	sh[0x8]=_mm256_max_epi32(sh[0x8], offset);
	sh[0x9]=_mm256_max_epi32(sh[0x9], offset);
	sh[0xA]=_mm256_max_epi32(sh[0xA], offset);
	sh[0xB]=_mm256_max_epi32(sh[0xB], offset);
	sh[0xC]=_mm256_max_epi32(sh[0xC], offset);
	sh[0xD]=_mm256_max_epi32(sh[0xD], offset);
	sh[0xE]=_mm256_max_epi32(sh[0xE], offset);
	sh[0xF]=_mm256_max_epi32(sh[0xF], offset);
	coeff[0x0]=_mm256_srlv_epi32(coeff[0x0], sh[0x0]);
	coeff[0x1]=_mm256_srlv_epi32(coeff[0x1], sh[0x1]);
	coeff[0x2]=_mm256_srlv_epi32(coeff[0x2], sh[0x2]);
	coeff[0x3]=_mm256_srlv_epi32(coeff[0x3], sh[0x3]);
	coeff[0x4]=_mm256_srlv_epi32(coeff[0x4], sh[0x4]);
	coeff[0x5]=_mm256_srlv_epi32(coeff[0x5], sh[0x5]);
	coeff[0x6]=_mm256_srlv_epi32(coeff[0x6], sh[0x6]);
	coeff[0x7]=_mm256_srlv_epi32(coeff[0x7], sh[0x7]);
	coeff[0x8]=_mm256_srlv_epi32(coeff[0x8], sh[0x8]);
	coeff[0x9]=_mm256_srlv_epi32(coeff[0x9], sh[0x9]);
	coeff[0xA]=_mm256_srlv_epi32(coeff[0xA], sh[0xA]);
	coeff[0xB]=_mm256_srlv_epi32(coeff[0xB], sh[0xB]);
	coeff[0xC]=_mm256_srlv_epi32(coeff[0xC], sh[0xC]);
	coeff[0xD]=_mm256_srlv_epi32(coeff[0xD], sh[0xD]);
	coeff[0xE]=_mm256_srlv_epi32(coeff[0xE], sh[0xE]);
	coeff[0xF]=_mm256_srlv_epi32(coeff[0xF], sh[0xF]);
	gather32((int*)(coeff+0x0), divlookup[0], (int*)(coeff+0x0));
	gather32((int*)(coeff+0x1), divlookup[0], (int*)(coeff+0x1));
	gather32((int*)(coeff+0x2), divlookup[1], (int*)(coeff+0x2));
	gather32((int*)(coeff+0x3), divlookup[1], (int*)(coeff+0x3));
	gather32((int*)(coeff+0x4), divlookup[2], (int*)(coeff+0x4));
	gather32((int*)(coeff+0x5), divlookup[2], (int*)(coeff+0x5));
	gather32((int*)(coeff+0x6), divlookup[3], (int*)(coeff+0x6));
	gather32((int*)(coeff+0x7), divlookup[3], (int*)(coeff+0x7));
	gather32((int*)(coeff+0x8), divlookup[4], (int*)(coeff+0x8));
	gather32((int*)(coeff+0x9), divlookup[4], (int*)(coeff+0x9));
	gather32((int*)(coeff+0xA), divlookup[5], (int*)(coeff+0xA));
	gather32((int*)(coeff+0xB), divlookup[5], (int*)(coeff+0xB));
	gather32((int*)(coeff+0xC), divlookup[6], (int*)(coeff+0xC));
	gather32((int*)(coeff+0xD), divlookup[6], (int*)(coeff+0xD));
	gather32((int*)(coeff+0xE), divlookup[7], (int*)(coeff+0xE));
	gather32((int*)(coeff+0xF), divlookup[7], (int*)(coeff+0xF));
	coeff[0x0]=_mm256_srlv_epi32(coeff[0x0], sh[0x0]);
	coeff[0x1]=_mm256_srlv_epi32(coeff[0x1], sh[0x1]);
	coeff[0x2]=_mm256_srlv_epi32(coeff[0x2], sh[0x2]);
	coeff[0x3]=_mm256_srlv_epi32(coeff[0x3], sh[0x3]);
	coeff[0x4]=_mm256_srlv_epi32(coeff[0x4], sh[0x4]);
	coeff[0x5]=_mm256_srlv_epi32(coeff[0x5], sh[0x5]);
	coeff[0x6]=_mm256_srlv_epi32(coeff[0x6], sh[0x6]);
	coeff[0x7]=_mm256_srlv_epi32(coeff[0x7], sh[0x7]);
	coeff[0x8]=_mm256_srlv_epi32(coeff[0x8], sh[0x8]);
	coeff[0x9]=_mm256_srlv_epi32(coeff[0x9], sh[0x9]);
	coeff[0xA]=_mm256_srlv_epi32(coeff[0xA], sh[0xA]);
	coeff[0xB]=_mm256_srlv_epi32(coeff[0xB], sh[0xB]);
	coeff[0xC]=_mm256_srlv_epi32(coeff[0xC], sh[0xC]);
	coeff[0xD]=_mm256_srlv_epi32(coeff[0xD], sh[0xD]);
	coeff[0xE]=_mm256_srlv_epi32(coeff[0xE], sh[0xE]);
	coeff[0xF]=_mm256_srlv_epi32(coeff[0xF], sh[0xF]);
	offset=_mm256_set1_epi32((1<<DENBITS>>2)/WG_NPREDS);
	coeff[0x0]=_mm256_add_epi32(coeff[0x0], offset);
	coeff[0x1]=_mm256_add_epi32(coeff[0x1], offset);
	coeff[0x2]=_mm256_add_epi32(coeff[0x2], offset);
	coeff[0x3]=_mm256_add_epi32(coeff[0x3], offset);
	coeff[0x4]=_mm256_add_epi32(coeff[0x4], offset);
	coeff[0x5]=_mm256_add_epi32(coeff[0x5], offset);
	coeff[0x6]=_mm256_add_epi32(coeff[0x6], offset);
	coeff[0x7]=_mm256_add_epi32(coeff[0x7], offset);
	coeff[0x8]=_mm256_add_epi32(coeff[0x8], offset);
	coeff[0x9]=_mm256_add_epi32(coeff[0x9], offset);
	coeff[0xA]=_mm256_add_epi32(coeff[0xA], offset);
	coeff[0xB]=_mm256_add_epi32(coeff[0xB], offset);
	coeff[0xC]=_mm256_add_epi32(coeff[0xC], offset);
	coeff[0xD]=_mm256_add_epi32(coeff[0xD], offset);
	coeff[0xE]=_mm256_add_epi32(coeff[0xE], offset);
	coeff[0xF]=_mm256_add_epi32(coeff[0xF], offset);
	__m256i sum0, sum1;
	sum0=_mm256_add_epi32(coeff[0*2+0], coeff[1*2+0]); sum1=_mm256_add_epi32(coeff[0*2+1], coeff[1*2+1]);
	sum0=_mm256_add_epi32(sum0, coeff[2*2+0]); sum1=_mm256_add_epi32(sum1, coeff[2*2+1]);
	sum0=_mm256_add_epi32(sum0, coeff[3*2+0]); sum1=_mm256_add_epi32(sum1, coeff[3*2+1]);
	sum0=_mm256_add_epi32(sum0, coeff[4*2+0]); sum1=_mm256_add_epi32(sum1, coeff[4*2+1]);
	sum0=_mm256_add_epi32(sum0, coeff[5*2+0]); sum1=_mm256_add_epi32(sum1, coeff[5*2+1]);
	sum0=_mm256_add_epi32(sum0, coeff[6*2+0]); sum1=_mm256_add_epi32(sum1, coeff[6*2+1]);
	sum0=_mm256_add_epi32(sum0, coeff[7*2+0]); sum1=_mm256_add_epi32(sum1, coeff[7*2+1]);
	__m256i sh0, sh1;
	sh0=FLOOR_LOG2_32x8(sum0);
	sh1=FLOOR_LOG2_32x8(sum1);
	offset=_mm256_set1_epi32(DENBITS-2);
	sh0=_mm256_sub_epi32(sh0, offset);
	sh1=_mm256_sub_epi32(sh1, offset);
	coeff[0*2+0]=_mm256_srlv_epi32(coeff[0*2+0], sh0); coeff[0*2+1]=_mm256_srlv_epi32(coeff[0*2+1], sh1);
	coeff[1*2+0]=_mm256_srlv_epi32(coeff[1*2+0], sh0); coeff[1*2+1]=_mm256_srlv_epi32(coeff[1*2+1], sh1);
	coeff[2*2+0]=_mm256_srlv_epi32(coeff[2*2+0], sh0); coeff[2*2+1]=_mm256_srlv_epi32(coeff[2*2+1], sh1);
	coeff[3*2+0]=_mm256_srlv_epi32(coeff[3*2+0], sh0); coeff[3*2+1]=_mm256_srlv_epi32(coeff[3*2+1], sh1);
	coeff[4*2+0]=_mm256_srlv_epi32(coeff[4*2+0], sh0); coeff[4*2+1]=_mm256_srlv_epi32(coeff[4*2+1], sh1);
	coeff[5*2+0]=_mm256_srlv_epi32(coeff[5*2+0], sh0); coeff[5*2+1]=_mm256_srlv_epi32(coeff[5*2+1], sh1);
	coeff[6*2+0]=_mm256_srlv_epi32(coeff[6*2+0], sh0); coeff[6*2+1]=_mm256_srlv_epi32(coeff[6*2+1], sh1);
	coeff[7*2+0]=_mm256_srlv_epi32(coeff[7*2+0], sh0); coeff[7*2+1]=_mm256_srlv_epi32(coeff[7*2+1], sh1);
	sum0=_mm256_add_epi32(coeff[0*2+0], coeff[1*2+0]); sum1=_mm256_add_epi32(coeff[0*2+1], coeff[1*2+1]);
	sum0=_mm256_add_epi32(sum0, coeff[2*2+0]); sum1=_mm256_add_epi32(sum1, coeff[2*2+1]);
	sum0=_mm256_add_epi32(sum0, coeff[3*2+0]); sum1=_mm256_add_epi32(sum1, coeff[3*2+1]);
	sum0=_mm256_add_epi32(sum0, coeff[4*2+0]); sum1=_mm256_add_epi32(sum1, coeff[4*2+1]);
	sum0=_mm256_add_epi32(sum0, coeff[5*2+0]); sum1=_mm256_add_epi32(sum1, coeff[5*2+1]);
	sum0=_mm256_add_epi32(sum0, coeff[6*2+0]); sum1=_mm256_add_epi32(sum1, coeff[6*2+1]);
	sum0=_mm256_add_epi32(sum0, coeff[7*2+0]); sum1=_mm256_add_epi32(sum1, coeff[7*2+1]);
#if 1
	__m256i ipred[2];
	ipred[0]=_mm256_srli_epi32(sum0, 1);
	ipred[1]=_mm256_srli_epi32(sum1, 1);
	sum0=_mm256_sub_epi32(sum0, one);
	sum1=_mm256_sub_epi32(sum1, one);
	ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[0*2+0], estim[0*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[0*2+1], estim[0*2+1]));
	ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[1*2+0], estim[1*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[1*2+1], estim[1*2+1]));
	ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[2*2+0], estim[2*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[2*2+1], estim[2*2+1]));
	ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[3*2+0], estim[3*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[3*2+1], estim[3*2+1]));
	ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[4*2+0], estim[4*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[4*2+1], estim[4*2+1]));
	ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[5*2+0], estim[5*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[5*2+1], estim[5*2+1]));
	ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[6*2+0], estim[6*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[6*2+1], estim[6*2+1]));
	ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[7*2+0], estim[7*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[7*2+1], estim[7*2+1]));
	sum0=_mm256_i32gather_epi32(divlookup[WG_NPREDS], sum0, sizeof(int));
	sum1=_mm256_i32gather_epi32(divlookup[WG_NPREDS], sum1, sizeof(int));
	ipred[0]=_mm256_mullo_epi32(ipred[0], sum0);
	ipred[1]=_mm256_mullo_epi32(ipred[1], sum1);
	ipred[0]=_mm256_srai_epi32(ipred[0], NUMBITS);
	ipred[1]=_mm256_srai_epi32(ipred[1], NUMBITS);
	ipred[0]=_mm256_slli_epi32(ipred[0], 16);
	ipred[1]=_mm256_slli_epi32(ipred[1], 16);
	ipred[0]=_mm256_srli_epi32(ipred[0], 16);
	result[0]=_mm256_or_si256(ipred[0], ipred[1]);
#endif
#if 0
	coeff[0*2+1]=_mm256_slli_epi32(coeff[0*2+1], 16);
	coeff[1*2+1]=_mm256_slli_epi32(coeff[1*2+1], 16);
	coeff[2*2+1]=_mm256_slli_epi32(coeff[2*2+1], 16);
	coeff[3*2+1]=_mm256_slli_epi32(coeff[3*2+1], 16);
	coeff[4*2+1]=_mm256_slli_epi32(coeff[4*2+1], 16);
	coeff[5*2+1]=_mm256_slli_epi32(coeff[5*2+1], 16);
	coeff[6*2+1]=_mm256_slli_epi32(coeff[6*2+1], 16);
	coeff[7*2+1]=_mm256_slli_epi32(coeff[7*2+1], 16);
	coeff[0*2+0]=_mm256_or_si256(coeff[0*2+0], coeff[0*2+1]);
	coeff[1*2+0]=_mm256_or_si256(coeff[1*2+0], coeff[1*2+1]);
	coeff[2*2+0]=_mm256_or_si256(coeff[2*2+0], coeff[2*2+1]);
	coeff[3*2+0]=_mm256_or_si256(coeff[3*2+0], coeff[3*2+1]);
	coeff[4*2+0]=_mm256_or_si256(coeff[4*2+0], coeff[4*2+1]);
	coeff[5*2+0]=_mm256_or_si256(coeff[5*2+0], coeff[5*2+1]);
	coeff[6*2+0]=_mm256_or_si256(coeff[6*2+0], coeff[6*2+1]);
	coeff[7*2+0]=_mm256_or_si256(coeff[7*2+0], coeff[7*2+1]);
#if 0
	__m256i tmp[8];
	tmp[0]=_mm256_mulhi_epi16(coeff[0*2+0], wgpreds[0*3+0]);
	tmp[1]=_mm256_mulhi_epi16(coeff[1*2+0], wgpreds[1*3+0]);
	tmp[2]=_mm256_mulhi_epi16(coeff[2*2+0], wgpreds[2*3+0]);
	tmp[3]=_mm256_mulhi_epi16(coeff[3*2+0], wgpreds[3*3+0]);
	tmp[4]=_mm256_mulhi_epi16(coeff[4*2+0], wgpreds[4*3+0]);
	tmp[5]=_mm256_mulhi_epi16(coeff[5*2+0], wgpreds[5*3+0]);
	tmp[6]=_mm256_mulhi_epi16(coeff[6*2+0], wgpreds[6*3+0]);
	tmp[7]=_mm256_mulhi_epi16(coeff[7*2+0], wgpreds[7*3+0]);
	for(int k=0;k<16*8;++k)
	{
		int val=((short*)tmp)[k];
		if(val!=0&&val!=-1)
			LOG_ERROR("");
	}
#endif
	coeff[0*2+0]=_mm256_mullo_epi16(coeff[0*2+0], wgpreds[0*3+0]);
	coeff[1*2+0]=_mm256_mullo_epi16(coeff[1*2+0], wgpreds[1*3+0]);
	coeff[2*2+0]=_mm256_mullo_epi16(coeff[2*2+0], wgpreds[2*3+0]);
	coeff[3*2+0]=_mm256_mullo_epi16(coeff[3*2+0], wgpreds[3*3+0]);
	coeff[4*2+0]=_mm256_mullo_epi16(coeff[4*2+0], wgpreds[4*3+0]);
	coeff[5*2+0]=_mm256_mullo_epi16(coeff[5*2+0], wgpreds[5*3+0]);
	coeff[6*2+0]=_mm256_mullo_epi16(coeff[6*2+0], wgpreds[6*3+0]);
	coeff[7*2+0]=_mm256_mullo_epi16(coeff[7*2+0], wgpreds[7*3+0]);
	
	__m256i bias0=_mm256_srli_epi32(sum0, 1);
	__m256i bias1=_mm256_srli_epi32(sum1, 1);
	sum0=_mm256_sub_epi32(sum0, one);
	sum1=_mm256_sub_epi32(sum1, one);
	sum0=_mm256_i32gather_epi32(divlookup[WG_NPREDS], sum0, sizeof(int));
	sum1=_mm256_i32gather_epi32(divlookup[WG_NPREDS], sum1, sizeof(int));

	//0: hi,  1: lo
	coeff[0*2+1]=_mm256_slli_epi32(coeff[0*2+0], 16);
	coeff[1*2+1]=_mm256_slli_epi32(coeff[1*2+0], 16);
	coeff[2*2+1]=_mm256_slli_epi32(coeff[2*2+0], 16);
	coeff[3*2+1]=_mm256_slli_epi32(coeff[3*2+0], 16);
	coeff[4*2+1]=_mm256_slli_epi32(coeff[4*2+0], 16);
	coeff[5*2+1]=_mm256_slli_epi32(coeff[5*2+0], 16);
	coeff[6*2+1]=_mm256_slli_epi32(coeff[6*2+0], 16);
	coeff[7*2+1]=_mm256_slli_epi32(coeff[7*2+0], 16);
	coeff[0*2+1]=_mm256_srai_epi32(coeff[0*2+1], 16);
	coeff[1*2+1]=_mm256_srai_epi32(coeff[1*2+1], 16);
	coeff[2*2+1]=_mm256_srai_epi32(coeff[2*2+1], 16);
	coeff[3*2+1]=_mm256_srai_epi32(coeff[3*2+1], 16);
	coeff[4*2+1]=_mm256_srai_epi32(coeff[4*2+1], 16);
	coeff[5*2+1]=_mm256_srai_epi32(coeff[5*2+1], 16);
	coeff[6*2+1]=_mm256_srai_epi32(coeff[6*2+1], 16);
	coeff[7*2+1]=_mm256_srai_epi32(coeff[7*2+1], 16);
	coeff[0*2+0]=_mm256_srai_epi32(coeff[0*2+0], 16);
	coeff[1*2+0]=_mm256_srai_epi32(coeff[1*2+0], 16);
	coeff[2*2+0]=_mm256_srai_epi32(coeff[2*2+0], 16);
	coeff[3*2+0]=_mm256_srai_epi32(coeff[3*2+0], 16);
	coeff[4*2+0]=_mm256_srai_epi32(coeff[4*2+0], 16);
	coeff[5*2+0]=_mm256_srai_epi32(coeff[5*2+0], 16);
	coeff[6*2+0]=_mm256_srai_epi32(coeff[6*2+0], 16);
	coeff[7*2+0]=_mm256_srai_epi32(coeff[7*2+0], 16);

	coeff[0*2+1]=_mm256_add_epi32(coeff[0*2+1], coeff[1*2+1]); coeff[0*2+0]=_mm256_add_epi32(coeff[0*2+0], coeff[1*2+0]);
	coeff[0*2+1]=_mm256_add_epi32(coeff[0*2+1], coeff[2*2+1]); coeff[0*2+0]=_mm256_add_epi32(coeff[0*2+0], coeff[2*2+0]);
	coeff[0*2+1]=_mm256_add_epi32(coeff[0*2+1], coeff[3*2+1]); coeff[0*2+0]=_mm256_add_epi32(coeff[0*2+0], coeff[3*2+0]);
	coeff[0*2+1]=_mm256_add_epi32(coeff[0*2+1], coeff[4*2+1]); coeff[0*2+0]=_mm256_add_epi32(coeff[0*2+0], coeff[4*2+0]);
	coeff[0*2+1]=_mm256_add_epi32(coeff[0*2+1], coeff[5*2+1]); coeff[0*2+0]=_mm256_add_epi32(coeff[0*2+0], coeff[5*2+0]);
	coeff[0*2+1]=_mm256_add_epi32(coeff[0*2+1], coeff[6*2+1]); coeff[0*2+0]=_mm256_add_epi32(coeff[0*2+0], coeff[6*2+0]);
	coeff[0*2+1]=_mm256_add_epi32(coeff[0*2+1], coeff[7*2+1]); coeff[0*2+0]=_mm256_add_epi32(coeff[0*2+0], coeff[7*2+0]);
	coeff[0*2+1]=_mm256_add_epi32(coeff[0*2+1], bias0);
	coeff[0*2+0]=_mm256_add_epi32(coeff[0*2+0], bias1);
	coeff[0*2+1]=_mm256_mullo_epi32(coeff[0*2+1], sum0);
	coeff[0*2+0]=_mm256_mullo_epi32(coeff[0*2+0], sum1);
	coeff[0*2+1]=_mm256_srai_epi32(coeff[0*2+1], NUMBITS);
	coeff[0*2+0]=_mm256_srai_epi32(coeff[0*2+0], NUMBITS);
	coeff[0*2+1]=_mm256_slli_epi32(coeff[0*2+1], 16);
	coeff[0*2+0]=_mm256_slli_epi32(coeff[0*2+0], 16);
	coeff[0*2+1]=_mm256_srli_epi32(coeff[0*2+1], 16);
	result[0]=_mm256_or_si256(coeff[0*2+1], coeff[0*2+0]);
#endif

#if 0
#define NUMBITS 20
#define DENBITS 7
#define DIVLUTSIZE (1<<DENBITS)
	static int divlookup[DIVLUTSIZE]={0};
	if(!*divlookup)
	{
		for(int k=0;k<DIVLUTSIZE;++k)
			divlookup[k]=(1<<NUMBITS)/(k+1);
	}
	for(int kl=0;kl<NCODERS/2;++kl)
	{
		static const int weights[]=
		{
			240,
			240,
			180,
			180,
			140,
			160,
			120,
			120,
		};
		long long coeff2[WG_NPREDS], ipred=0, wsum=0;
		int sh=0;
		for(int kp=0;kp<WG_NPREDS;++kp)
		{
			int e=((int*)coeff)[NCODERS/2*kp+kl];
			sh=FLOOR_LOG2(e+1);
			coeff2[kp]=((long long)weights[kp]*divlookup[e<<(DENBITS-1)>>sh]<<(DENBITS-1)>>sh)+(1<<DENBITS>>1)/WG_NPREDS;
			wsum+=coeff2[kp];
		}
		sh=FLOOR_LOG2(wsum)-(DENBITS-2);
		wsum=0;
		for(int kp=0;kp<WG_NPREDS;++kp)
		{
			long long c=coeff2[kp]>>=sh;
			wsum+=c;
			coeff2[kp]=c;
		}
		ipred=(wsum>>1)-1LL;
		for(int kp=0;kp<WG_NPREDS;++kp)
			ipred+=(long long)coeff2[kp]*((int*)estim)[NCODERS/2*kp+kl];
		((int*)result)[kl]=(int)(ipred*divlookup[wsum-1]>>NUMBITS);
	}
#endif
}
AWM_INLINE void wg_mix(int kc, const short *wgWerrors, const short *rowscurrNNptr, const short *rowsNptr, const __m256i *wgpreds, __m256i *preds)//process channel wave (NCODERS lanes)
{
	//pred/errors layout:  {Y[NCODERS], U[NCODERS], V[NCODERS]} * NPREDS
	const short *chWerrors=wgWerrors+(ptrdiff_t)NCODERS*kc;
	const short *Nptr=rowsNptr+(ptrdiff_t)NCODERS*kc;
	const short *currNNptr=rowscurrNNptr+(ptrdiff_t)NCODERS*kc;
	__m256i prederrors[sizeof(short[WG_NPREDS*NCODERS])/(sizeof(__m256i))];//2 regs per lane  *  NPREDS

	//errors = pI+pNW+2*pN+pNE+pNNE		(__m256i*)erows[-Y]+(P+X*NPREDS)*3+C
	prederrors[0]=_mm256_load_si256((const __m256i*)currNNptr+(0-1*WG_NPREDS)*3+0);//pW
	prederrors[1]=_mm256_load_si256((const __m256i*)currNNptr+(1-1*WG_NPREDS)*3+0);
	prederrors[2]=_mm256_load_si256((const __m256i*)currNNptr+(2-1*WG_NPREDS)*3+0);
	prederrors[3]=_mm256_load_si256((const __m256i*)currNNptr+(3-1*WG_NPREDS)*3+0);
	prederrors[4]=_mm256_load_si256((const __m256i*)currNNptr+(4-1*WG_NPREDS)*3+0);
	prederrors[5]=_mm256_load_si256((const __m256i*)currNNptr+(5-1*WG_NPREDS)*3+0);
	prederrors[6]=_mm256_load_si256((const __m256i*)currNNptr+(6-1*WG_NPREDS)*3+0);
	prederrors[7]=_mm256_load_si256((const __m256i*)currNNptr+(7-1*WG_NPREDS)*3+0);
	
	prederrors[0]=_mm256_add_epi16(prederrors[0], _mm256_load_si256((const __m256i*)Nptr+(0+0*WG_NPREDS)*3+0));//pN
	prederrors[1]=_mm256_add_epi16(prederrors[1], _mm256_load_si256((const __m256i*)Nptr+(1+0*WG_NPREDS)*3+0));
	prederrors[2]=_mm256_add_epi16(prederrors[2], _mm256_load_si256((const __m256i*)Nptr+(2+0*WG_NPREDS)*3+0));
	prederrors[3]=_mm256_add_epi16(prederrors[3], _mm256_load_si256((const __m256i*)Nptr+(3+0*WG_NPREDS)*3+0));
	prederrors[4]=_mm256_add_epi16(prederrors[4], _mm256_load_si256((const __m256i*)Nptr+(4+0*WG_NPREDS)*3+0));
	prederrors[5]=_mm256_add_epi16(prederrors[5], _mm256_load_si256((const __m256i*)Nptr+(5+0*WG_NPREDS)*3+0));
	prederrors[6]=_mm256_add_epi16(prederrors[6], _mm256_load_si256((const __m256i*)Nptr+(6+0*WG_NPREDS)*3+0));
	prederrors[7]=_mm256_add_epi16(prederrors[7], _mm256_load_si256((const __m256i*)Nptr+(7+0*WG_NPREDS)*3+0));
	
	prederrors[0]=_mm256_add_epi16(prederrors[0], _mm256_load_si256((const __m256i*)Nptr+(0-1*WG_NPREDS)*3+0));//pNW
	prederrors[1]=_mm256_add_epi16(prederrors[1], _mm256_load_si256((const __m256i*)Nptr+(1-1*WG_NPREDS)*3+0));
	prederrors[2]=_mm256_add_epi16(prederrors[2], _mm256_load_si256((const __m256i*)Nptr+(2-1*WG_NPREDS)*3+0));
	prederrors[3]=_mm256_add_epi16(prederrors[3], _mm256_load_si256((const __m256i*)Nptr+(3-1*WG_NPREDS)*3+0));
	prederrors[4]=_mm256_add_epi16(prederrors[4], _mm256_load_si256((const __m256i*)Nptr+(4-1*WG_NPREDS)*3+0));
	prederrors[5]=_mm256_add_epi16(prederrors[5], _mm256_load_si256((const __m256i*)Nptr+(5-1*WG_NPREDS)*3+0));
	prederrors[6]=_mm256_add_epi16(prederrors[6], _mm256_load_si256((const __m256i*)Nptr+(6-1*WG_NPREDS)*3+0));
	prederrors[7]=_mm256_add_epi16(prederrors[7], _mm256_load_si256((const __m256i*)Nptr+(7-1*WG_NPREDS)*3+0));

	prederrors[0]=_mm256_slli_epi16(prederrors[0], 1);
	prederrors[1]=_mm256_slli_epi16(prederrors[1], 1);
	prederrors[2]=_mm256_slli_epi16(prederrors[2], 1);
	prederrors[3]=_mm256_slli_epi16(prederrors[3], 1);
	prederrors[4]=_mm256_slli_epi16(prederrors[4], 1);
	prederrors[5]=_mm256_slli_epi16(prederrors[5], 1);
	prederrors[6]=_mm256_slli_epi16(prederrors[6], 1);
	prederrors[7]=_mm256_slli_epi16(prederrors[7], 1);
	
	prederrors[0]=_mm256_add_epi16(prederrors[0], _mm256_load_si256((const __m256i*)currNNptr+(0+0*WG_NPREDS)*3+0));//pNN
	prederrors[1]=_mm256_add_epi16(prederrors[1], _mm256_load_si256((const __m256i*)currNNptr+(1+0*WG_NPREDS)*3+0));
	prederrors[2]=_mm256_add_epi16(prederrors[2], _mm256_load_si256((const __m256i*)currNNptr+(2+0*WG_NPREDS)*3+0));
	prederrors[3]=_mm256_add_epi16(prederrors[3], _mm256_load_si256((const __m256i*)currNNptr+(3+0*WG_NPREDS)*3+0));
	prederrors[4]=_mm256_add_epi16(prederrors[4], _mm256_load_si256((const __m256i*)currNNptr+(4+0*WG_NPREDS)*3+0));
	prederrors[5]=_mm256_add_epi16(prederrors[5], _mm256_load_si256((const __m256i*)currNNptr+(5+0*WG_NPREDS)*3+0));
	prederrors[6]=_mm256_add_epi16(prederrors[6], _mm256_load_si256((const __m256i*)currNNptr+(6+0*WG_NPREDS)*3+0));
	prederrors[7]=_mm256_add_epi16(prederrors[7], _mm256_load_si256((const __m256i*)currNNptr+(7+0*WG_NPREDS)*3+0));
	
	prederrors[0]=_mm256_add_epi16(prederrors[0], _mm256_load_si256((const __m256i*)currNNptr+(0+1*WG_NPREDS)*3+0));//pNNE
	prederrors[1]=_mm256_add_epi16(prederrors[1], _mm256_load_si256((const __m256i*)currNNptr+(1+1*WG_NPREDS)*3+0));
	prederrors[2]=_mm256_add_epi16(prederrors[2], _mm256_load_si256((const __m256i*)currNNptr+(2+1*WG_NPREDS)*3+0));
	prederrors[3]=_mm256_add_epi16(prederrors[3], _mm256_load_si256((const __m256i*)currNNptr+(3+1*WG_NPREDS)*3+0));
	prederrors[4]=_mm256_add_epi16(prederrors[4], _mm256_load_si256((const __m256i*)currNNptr+(4+1*WG_NPREDS)*3+0));
	prederrors[5]=_mm256_add_epi16(prederrors[5], _mm256_load_si256((const __m256i*)currNNptr+(5+1*WG_NPREDS)*3+0));
	prederrors[6]=_mm256_add_epi16(prederrors[6], _mm256_load_si256((const __m256i*)currNNptr+(6+1*WG_NPREDS)*3+0));
	prederrors[7]=_mm256_add_epi16(prederrors[7], _mm256_load_si256((const __m256i*)currNNptr+(7+1*WG_NPREDS)*3+0));
	
	prederrors[0]=_mm256_add_epi16(prederrors[0], _mm256_load_si256((const __m256i*)Nptr+(0+1*WG_NPREDS)*3+0));//pNE
	prederrors[1]=_mm256_add_epi16(prederrors[1], _mm256_load_si256((const __m256i*)Nptr+(1+1*WG_NPREDS)*3+0));
	prederrors[2]=_mm256_add_epi16(prederrors[2], _mm256_load_si256((const __m256i*)Nptr+(2+1*WG_NPREDS)*3+0));
	prederrors[3]=_mm256_add_epi16(prederrors[3], _mm256_load_si256((const __m256i*)Nptr+(3+1*WG_NPREDS)*3+0));
	prederrors[4]=_mm256_add_epi16(prederrors[4], _mm256_load_si256((const __m256i*)Nptr+(4+1*WG_NPREDS)*3+0));
	prederrors[5]=_mm256_add_epi16(prederrors[5], _mm256_load_si256((const __m256i*)Nptr+(5+1*WG_NPREDS)*3+0));
	prederrors[6]=_mm256_add_epi16(prederrors[6], _mm256_load_si256((const __m256i*)Nptr+(6+1*WG_NPREDS)*3+0));
	prederrors[7]=_mm256_add_epi16(prederrors[7], _mm256_load_si256((const __m256i*)Nptr+(7+1*WG_NPREDS)*3+0));
	
	prederrors[0]=_mm256_add_epi16(prederrors[0], _mm256_load_si256((const __m256i*)currNNptr+(0-2*WG_NPREDS)*3+0));//pWW
	prederrors[1]=_mm256_add_epi16(prederrors[1], _mm256_load_si256((const __m256i*)currNNptr+(1-2*WG_NPREDS)*3+0));
	prederrors[2]=_mm256_add_epi16(prederrors[2], _mm256_load_si256((const __m256i*)currNNptr+(2-2*WG_NPREDS)*3+0));
	prederrors[3]=_mm256_add_epi16(prederrors[3], _mm256_load_si256((const __m256i*)currNNptr+(3-2*WG_NPREDS)*3+0));
	prederrors[4]=_mm256_add_epi16(prederrors[4], _mm256_load_si256((const __m256i*)currNNptr+(4-2*WG_NPREDS)*3+0));
	prederrors[5]=_mm256_add_epi16(prederrors[5], _mm256_load_si256((const __m256i*)currNNptr+(5-2*WG_NPREDS)*3+0));
	prederrors[6]=_mm256_add_epi16(prederrors[6], _mm256_load_si256((const __m256i*)currNNptr+(6-2*WG_NPREDS)*3+0));
	prederrors[7]=_mm256_add_epi16(prederrors[7], _mm256_load_si256((const __m256i*)currNNptr+(7-2*WG_NPREDS)*3+0));
#ifndef WG_ENABLE_eW
	(void)chWerrors;
#else
	prederrors[0]=_mm256_add_epi16(prederrors[0], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)chWerrors+(0*3+0)), 2));//pI	(__m256i*)(wgWerrors+(P*3+C)*NCODERS)+R
	prederrors[1]=_mm256_add_epi16(prederrors[1], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)chWerrors+(1*3+0)), 2));
	prederrors[2]=_mm256_add_epi16(prederrors[2], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)chWerrors+(2*3+0)), 2));
	prederrors[3]=_mm256_add_epi16(prederrors[3], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)chWerrors+(3*3+0)), 2));
	prederrors[4]=_mm256_add_epi16(prederrors[4], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)chWerrors+(4*3+0)), 2));
	prederrors[5]=_mm256_add_epi16(prederrors[5], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)chWerrors+(5*3+0)), 2));
	prederrors[6]=_mm256_add_epi16(prederrors[6], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)chWerrors+(6*3+0)), 2));
	prederrors[7]=_mm256_add_epi16(prederrors[7], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)chWerrors+(7*3+0)), 2));
#endif

	wg_mix_pt2(wgpreds+kc, prederrors, preds);

//	__m256i p2;
//	int mask;
//again:
//	wg_mix_pt2v0(wgpreds+kc, prederrors, preds);
//
//	wg_mix_pt2(wgpreds+kc, prederrors, &p2);
//	p2=_mm256_xor_si256(p2, *preds);
//	mask=_mm256_movemask_epi8(p2);
//	if(mask)
//		goto again;

	//__m256i ires[2];
	//wg_mix_pt2(wgpreds+kc, prederrors, ires);//P*3+C
	//ires[0]=_mm256_slli_epi32(ires[0], 16);
	//ires[1]=_mm256_slli_epi32(ires[1], 16);
	//ires[0]=_mm256_srli_epi32(ires[0], 16);
	//preds[0]=_mm256_or_si256(ires[0], ires[1]);
}


//https://github.com/rygorous/ryg_rans
//https://github.com/samtools/htscodecs
typedef struct _rANS_SIMD_SymInfo	//16 bytes/level	4KB/ctx = 1<<12 bytes
{
	unsigned smax, invf, cdf;
	unsigned short negf, sh;
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
				LOG_ERROR("");
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
			info->sh=FLOOR_LOG2(freq);//eg: x/2 = x*0x80000000>>32>>0
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
	unsigned short CDF[257];
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
		int nbypass=FLOOR_LOG2(CDFlevels);
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
	unsigned short hist[257];
	if(ctxmask>>ctxidx&1)//rare context
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
			int nbypass=FLOOR_LOG2(CDFlevels);
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
					LOG_ERROR("CDF unpack error");
#endif
				break;
			}
		}
		if(CDFlevels)
			LOG_ERROR("CDF unpack error");
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
AWM_INLINE void dec_yuv(
	__m256i *mstate,
	const __m256i *ctx0,
	const int *CDF2syms,
	const int *ans_permute,
	unsigned char **pstreamptr,
	const unsigned char *streamend,
	__m256i *syms
)
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
	gather32((int*)(decctx+0), (int*)CDF2syms, (int*)(decctx+0));
	gather32((int*)(decctx+1), (int*)CDF2syms, (int*)(decctx+1));
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
		ansval_check(mdebugfreq, sizeof(short), NCODERS);
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
		__m256i idx0, idx1, lo0, lo1, renorm0, renorm1;
		__m256i smin=_mm256_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
		__m256i cond0=_mm256_cmpgt_epi32(smin, mstate[0]);//FIXME this is signed comparison
		__m256i cond1=_mm256_cmpgt_epi32(smin, mstate[1]);
		int mask0=_mm256_movemask_ps(_mm256_castsi256_ps(cond0));
		int mask1=_mm256_movemask_ps(_mm256_castsi256_ps(cond1));
		idx0=_mm256_load_si256((const __m256i*)ans_permute+mask0);
		idx1=_mm256_load_si256((const __m256i*)ans_permute+mask1);
		mask0=_mm_popcnt_u32(mask0);
		mask1=_mm_popcnt_u32(mask1);
#ifdef _DEBUG
		if(streamptr+((ptrdiff_t)mask0+mask1)*sizeof(short)>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		lo0=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr)); streamptr+=mask0*sizeof(short);
		lo1=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr)); streamptr+=mask1*sizeof(short);
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
		LOG_ERROR("Cannot open \"%s\" for writing", fn);
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
				LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, streamend);
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
				LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, streamend);
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
				LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, streamend);
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
				LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
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
				LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
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
				LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
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
#if 0
static void matmul3(double *res, const double *mleft, const double *mright)
{
	for(int ky=0;ky<3;++ky)
	{
		for(int kx=0;kx<3;++kx)
		{
			double r=0;
			for(int k=0;k<3;++k)
				r+=mleft[3*ky+k]*mright[3*k+kx];
			res[3*ky+kx]=r;
		}
	}
}
static void pca(const double *cov, int16_t *lossytransform)
{
	double cov0[]=
	{
		/*
		0  1  2  3  4  5
		rr gg bb rg gb br

		rr rg br
		rg gg gb
		br gb bb

		0 3 5
		3 1 4
		5 4 2
		*/
		cov[0], cov[3], cov[5],
		cov[3], cov[1], cov[4],
		cov[5], cov[4], cov[2],
	};
	double trans[]={1, 0, 0, 0, 1, 0, 0, 0, 1};

	for(int k=0;k<100;++k)
	{
		static const int pairs[][2]=
		{//	{i, j}
			{0, 1},
			{0, 2},
			{1, 2},
		};

		for(int k2=0;k2<3;++k2)
		{
			double a, b, d;
			double tau, t, s, c;
			int i=pairs[k2][0], j=pairs[k2][1];
			double temp[9];

			a=cov0[3*i+i];
			b=cov0[3*j+j];
			d=cov0[3*i+j];
			tau=(b-a)/(2*d);
			t=((tau>0)-(tau<0))/(fabs(tau)+sqrt(1+tau*tau));
			c=1/sqrt(1+t*t);
			s=t*c;

			/*
			Jacobi rotations
			(0,1)
			c	s	0
			-s	c	0
			0	0	1

			(0,2)
			c	0	s
			0	1	0
			-s	0	c

			(1,2)
			1	0	0
			0	c	s
			0	-s	c

			T' = T * J
			C' = JT * C * J
			*/
			double jacobi[]={1, 0, 0, 0, 1, 0, 0, 0, 1};
			jacobi[3*i+i]=c;
			jacobi[3*i+j]=s;
			jacobi[3*j+i]=-s;
			jacobi[3*j+j]=c;
			matmul3(temp, trans, jacobi);
			memcpy(trans, temp, sizeof(temp));

			matmul3(temp, cov0, jacobi);
			memcpy(cov0, temp, sizeof(temp));
			jacobi[3*i+j]=-s;//transpose the Jacobi
			jacobi[3*j+i]=s;
			matmul3(temp, jacobi, cov0);
			memcpy(cov0, temp, sizeof(temp));
		}
		if(fabs(cov0[1])+fabs(cov0[2])+fabs(cov0[5])<1e-6)
			break;
#if 0
#define GETPARAMS(A, B, D)\
tau=((B)-(A))/(2*(D)), t=((tau>0)-(tau<0))/(fabs(tau)+sqrt(1+tau*tau)), c=1/sqrt(1+t*t), s=t*c
		GETPARAMS(cov0[0], cov0[4], cov0[1]);
#endif
	}
				
	CLAMP2(trans[0], -15, 15);
	CLAMP2(trans[1], -15, 15);
	CLAMP2(trans[2], -15, 15);
	CLAMP2(trans[3], -15, 15);
	CLAMP2(trans[4], -15, 15);
	CLAMP2(trans[5], -15, 15);
	CLAMP2(trans[6], -15, 15);
	CLAMP2(trans[7], -15, 15);
	CLAMP2(trans[8], -15, 15);
	lossytransform[0]=(int16_t)(trans[0]*0x4000);//inverse (transpose) matrix
	lossytransform[1]=(int16_t)(trans[1]*0x4000);
	lossytransform[2]=(int16_t)(trans[2]*0x4000);
	lossytransform[3]=(int16_t)(trans[3]*0x4000);
	lossytransform[4]=(int16_t)(trans[4]*0x4000);
	lossytransform[5]=(int16_t)(trans[5]*0x4000);
	lossytransform[6]=(int16_t)(trans[6]*0x4000);
	lossytransform[7]=(int16_t)(trans[7]*0x4000);
	lossytransform[8]=(int16_t)(trans[8]*0x4000);
}
#endif
int c32_codec(int argc, char **argv)
{
	if(argc!=3&&argc!=4&&argc!=5)
	{
		printf(
			"Usage: \"%s\"  input  output  [P]  [D]    Encode/decode.\n"
			"P  =  1 Force CG / 2 Force WG4 | 4 Profile\n"
			"D  =  lossy distortion. 4 <= D <= 16\n"
			, argv[0]
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];
	int nthreads0=argc<4?0:atoi(argv[3]), dist=argc<5?1:atoi(argv[4]);
	if(dist>1)
		CLAMP2(dist, 3, 31);
#ifdef ESTIMATE_SIZE
	double esize[3*NCODERS]={0};
#endif
#ifdef LOUD
	ptrdiff_t usize2=get_filesize(srcfn);
	double t=time_sec();
#endif
	prof_checkpoint(0, 0);
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	int fwd=0, iw=0, ih=0, rowstride=0;
	ptrdiff_t usize=0, cap=0;
	unsigned char *image=0, *imptr=0, *streamptr=0, *streamstart=0, *streamend=0;
	int psize=0;
	short *pixels=0;
	int wgsize=0, wgstatesize=0;
	short *wgerrors=0, *wgstate=0;
	ptrdiff_t cheadersize=0, csize=0;
	{
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		int tag=0;
		fread(&tag, 1, 2, fsrc);
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('3'|'2'<<8))
		{
			LOG_ERROR("Unsupported file \"%s\"", srcfn);
			return 1;
		}
		if(fwd)
		{
#ifdef LOUD
			print_timestamp("%Y-%m-%d_%H%M%S\n");
#endif
			int temp=fgetc(fsrc);
			if(temp!='\n')
			{
				LOG_ERROR("Invalid PPM file");
				return 1;
			}
			int nread=fscanf(fsrc, "%d %d", &iw, &ih);
			if(nread!=2)
			{
				LOG_ERROR("Unsupported PPM file");
				return 1;
			}
			int vmax=0;
			nread=fscanf(fsrc, "%d", &vmax);
			if(nread!=1||vmax!=255)
			{
				LOG_ERROR("Unsupported PPM file");
				return 1;
			}
			temp=fgetc(fsrc);
			if(temp!='\n')
			{
				LOG_ERROR("Invalid PPM file");
				return 1;
			}
		}
		else
		{
			fread(&iw, 1, 4, fsrc);
			fread(&ih, 1, 4, fsrc);
			cheadersize=ftell(fsrc);
		}
		if(iw<1||ih<1)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		rowstride=3*iw;
		usize=(ptrdiff_t)3*iw*ih;
		cap=(ptrdiff_t)4*iw*ih;
	//	image=(unsigned char*)_mm_malloc(cap+sizeof(__m256i), 0x1000);
		image=(unsigned char*)malloc(cap+sizeof(__m256i));
		if(!image)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		if(fwd)
		{
			fread(image, 1, usize, fsrc);//read image
			streamptr=streamstart=image+cap;//bwd-bwd ANS encoding
			profile_size(streamptr, "start");
		}
		else
		{
			csize=get_filesize(srcfn);
			streamptr=streamstart=image+cap-(csize-cheadersize)-sizeof(__m256i);
			streamend=image+cap-sizeof(__m256i);
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
	ptrdiff_t isize=(ptrdiff_t)ixbytes*blockh;
	ptrdiff_t interleavedsize=isize<<fwd;//fwd ? interleave residuals & context : pack residuals
	unsigned char *interleaved=(unsigned char*)_mm_malloc(interleavedsize, sizeof(__m256i));
	if(!interleaved)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	(void)xrembytes;
	int bestrct=0, use_wg4=0;
	//ALIGN(32) int16_t lossymean[3*NCODERS]={0}, lossytransform[9*NCODERS]={0};
	//int L1sh=0;
	unsigned long long ctxmask=0;//3*NCTX+3 = 54 flags	0: rare context (bypass)  1: emit stats
	const int hsize=(int)sizeof(int[3*NCTX<<8]);//3 channels
	int *hists=fwd?(int*)malloc(hsize):0;//fwd-only
	const int rhsize=(int)sizeof(int[3*256]);
	int *rhist=fwd?(int*)malloc(rhsize):0;

	int CDF2syms_size=(int)sizeof(int[3*NCTX<<PROBBITS]);
	if(fwd)//DIV-free rANS encoder reuses these as SIMD symbol info
		CDF2syms_size=(int)sizeof(rANS_SIMD_SymInfo[3*NCTX<<8]);
	unsigned *CDF2syms=(unsigned*)_mm_malloc(CDF2syms_size, sizeof(__m256i));

	int rCDF2syms_size=(int)sizeof(int[3<<PROBBITS]);
	if(fwd)
		rCDF2syms_size=(int)sizeof(rANS_SIMD_SymInfo[3<<8]);
	unsigned *rCDF2syms=(unsigned*)_mm_malloc(rCDF2syms_size, sizeof(__m256i));

	int ans_permute_size=sizeof(__m256i[256]);
	int *ans_permute=(int*)_mm_malloc(ans_permute_size, sizeof(__m256i));

	psize=(int)sizeof(short[4*6*NCODERS])*(blockw+16);//4 padded rows  *  {Y*NCODERS, U*NCODERS, V*NCODERS,  eY*NCODERS, eU*NCODERS, eV*NCODERS} = 2*3*32 = 192 channels  ~48*iw bytes
	pixels=(short*)_mm_malloc(psize, sizeof(__m256i));//~188 KB for 4K/12MP
	if((fwd&&(!hists||!rhist))||!CDF2syms||!rCDF2syms||!ans_permute||!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(ans_permute, 0, ans_permute_size);//_mm256_permutevar8x32_epi32 can't clear elements like _mm256_shuffle_epi8 does
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
			LOG_ERROR("ERROR");
		printf("SUCCESS\n");
		exit(0);
#endif
		interleave_blocks_fwd(image, iw, ih, interleaved+isize);//reuse memory: read 8-bit pixel, write 16-bit context<<8|residual
		guide_save(interleaved+isize, ixcount, blockh);
		prof_checkpoint(usize, "interleave");
		if(dist<=1)//lossless analysis
		{
			ALIGN(32) long long counters[OCH_COUNT]={0};
			__m256i mcounters[OCH_COUNT];//64-bit
			__m128i half8=_mm_set1_epi8(-128);
			__m256i wordmask=_mm256_set1_epi64x(0xFFFF);
			memset(mcounters, 0, sizeof(mcounters));
			imptr=interleaved+isize;
			for(int ky=0;ky<blockh;++ky)
			{
				__m256i prev[OCH_COUNT];//16-bit
				memset(prev, 0, sizeof(prev));
				for(int kx=0;kx<ixbytes-3*NCODERS;kx+=3*NCODERS)
				{
					__m256i r=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+0), half8));
					__m256i g=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+1), half8));
					__m256i b=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+2), half8));
					imptr+=3*NCODERS;
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

//			bestrct=0;//
//#ifdef __GNUC__
//#error remove above
//#endif
			//printf("%2d ", bestrct);
		}
#if 0
		else//lossy analysis
		{
			ALIGN(32) int64_t sums[6*NCODERS]={0};
			__m128i half8=_mm_set1_epi8(-128);
#if 1
			ALIGN(32) int64_t means[3*NCODERS]={0};
			imptr=interleaved+isize;
			for(int ky=0;ky<blockh;++ky)
			{
				__m256i mmean[3*2];//64-bit
				memset(mmean, 0, sizeof(mmean));
				for(int kx=0;kx<ixbytes-3*NCODERS;kx+=3*NCODERS)
				{
					__m256i r=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+0), half8));
					__m256i g=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+1), half8));
					__m256i b=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+2), half8));
					imptr+=3*NCODERS;

					mmean[0]=_mm256_add_epi32(mmean[0], _mm256_srai_epi32(_mm256_slli_epi32(r, 16), 16));
					mmean[1]=_mm256_add_epi32(mmean[1], _mm256_srai_epi32(_mm256_slli_epi32(g, 16), 16));
					mmean[2]=_mm256_add_epi32(mmean[2], _mm256_srai_epi32(_mm256_slli_epi32(b, 16), 16));
					mmean[3]=_mm256_add_epi32(mmean[3], _mm256_srai_epi32(r, 16));
					mmean[4]=_mm256_add_epi32(mmean[4], _mm256_srai_epi32(g, 16));
					mmean[5]=_mm256_add_epi32(mmean[5], _mm256_srai_epi32(b, 16));
#if 0
#define ACC16_32(X) _mm256_add_epi32(_mm256_srai_epi32(X, 16), _mm256_srai_epi32(_mm256_slli_epi32(X, 16), 16))
					r=ACC16_32(r);
					g=ACC16_32(g);
					b=ACC16_32(b);
					mmean[0]=_mm256_add_epi32(mmean[0], r);
					mmean[1]=_mm256_add_epi32(mmean[1], g);
					mmean[2]=_mm256_add_epi32(mmean[2], b);
#endif
				}
				{
					ALIGN(32) int tmp[3*NCODERS];
					memcpy(tmp, mmean, sizeof(tmp));
					for(int k=0;k<3*NCODERS/2;++k)//lo
						means[2*k+0]+=tmp[k];
					for(int k=0;k<3*NCODERS/2;++k)//hi
						means[2*k+1]+=tmp[k+3*NCODERS/2];
				}
			}
#ifdef ONE_LOSSYMEAN
			{
				int64_t sum[3]={0};
				for(int k=0;k<NCODERS;++k)
				{
					sum[0]+=means[NCODERS*0+k];
					sum[1]+=means[NCODERS*1+k];
					sum[2]+=means[NCODERS*2+k];
				}
				for(int k=0;k<NCODERS;++k)
				{
					means[NCODERS*0+k]=sum[0];
					means[NCODERS*1+k]=sum[1];
					means[NCODERS*2+k]=sum[2];
				}
			}
#endif
			{
#ifdef ONE_LOSSYMEAN
				double invcount=1./((double)blockh*ixcount);
#else
				double invcount=1./((double)blockh*blockw);
#endif
				for(int k=0;k<3*NCODERS;++k)
					lossymean[k]=(int16_t)round((double)means[k]*invcount);
			}
			__m256i mrmean=_mm256_load_si256((__m256i*)lossymean+0);
			__m256i mgmean=_mm256_load_si256((__m256i*)lossymean+1);
			__m256i mbmean=_mm256_load_si256((__m256i*)lossymean+2);
			//__m256i mrmean=_mm256_set1_epi16(lossymean[0]);
			//__m256i mgmean=_mm256_set1_epi16(lossymean[1]);
			//__m256i mbmean=_mm256_set1_epi16(lossymean[2]);
#endif
			imptr=interleaved+isize;
			for(int ky=0;ky<blockh;++ky)
			{
				__m256i mcov[6*2];//64-bit
				memset(mcov, 0, sizeof(mcov));
				for(int kx=0;kx<ixbytes-3*NCODERS;kx+=3*NCODERS)
				{
					__m256i r=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+0), half8));
					__m256i g=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+1), half8));
					__m256i b=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+2), half8));
					r=_mm256_sub_epi16(r, mrmean);
					g=_mm256_sub_epi16(g, mgmean);
					b=_mm256_sub_epi16(b, mbmean);
					imptr+=3*NCODERS;

					__m256i rlo=_mm256_srai_epi32(_mm256_slli_epi32(r, 16), 16);
					__m256i glo=_mm256_srai_epi32(_mm256_slli_epi32(g, 16), 16);
					__m256i blo=_mm256_srai_epi32(_mm256_slli_epi32(b, 16), 16);
					__m256i rhi=_mm256_srai_epi32(r, 16);
					__m256i ghi=_mm256_srai_epi32(g, 16);
					__m256i bhi=_mm256_srai_epi32(b, 16);
					mcov[0+0]=_mm256_add_epi32(mcov[0+0], _mm256_mullo_epi32(rlo, rlo));
					mcov[0+1]=_mm256_add_epi32(mcov[0+1], _mm256_mullo_epi32(glo, glo));
					mcov[0+2]=_mm256_add_epi32(mcov[0+2], _mm256_mullo_epi32(blo, blo));
					mcov[0+3]=_mm256_add_epi32(mcov[0+3], _mm256_mullo_epi32(rlo, glo));
					mcov[0+4]=_mm256_add_epi32(mcov[0+4], _mm256_mullo_epi32(glo, blo));
					mcov[0+5]=_mm256_add_epi32(mcov[0+5], _mm256_mullo_epi32(blo, rlo));
					mcov[6+0]=_mm256_add_epi32(mcov[6+0], _mm256_mullo_epi32(rhi, rhi));
					mcov[6+1]=_mm256_add_epi32(mcov[6+1], _mm256_mullo_epi32(ghi, ghi));
					mcov[6+2]=_mm256_add_epi32(mcov[6+2], _mm256_mullo_epi32(bhi, bhi));
					mcov[6+3]=_mm256_add_epi32(mcov[6+3], _mm256_mullo_epi32(rhi, ghi));
					mcov[6+4]=_mm256_add_epi32(mcov[6+4], _mm256_mullo_epi32(ghi, bhi));
					mcov[6+5]=_mm256_add_epi32(mcov[6+5], _mm256_mullo_epi32(bhi, rhi));
#if 0
					__m256i rr=_mm256_mullo_epi16(r, r);
					__m256i gg=_mm256_mullo_epi16(g, g);
					__m256i bb=_mm256_mullo_epi16(b, b);
					__m256i rg=_mm256_mullo_epi16(r, g);
					__m256i gb=_mm256_mullo_epi16(g, b);
					__m256i br=_mm256_mullo_epi16(b, r);
#define ACC16_32(X) _mm256_add_epi32(_mm256_srai_epi32(X, 16), _mm256_srai_epi32(_mm256_slli_epi32(X, 16), 16))
					rr=ACC16_32(rr);
					gg=ACC16_32(gg);
					bb=ACC16_32(bb);
					rg=ACC16_32(rg);
					gb=ACC16_32(gb);
					br=ACC16_32(br);
					mcov[0]=_mm256_add_epi32(mcov[0], rr);
					mcov[1]=_mm256_add_epi32(mcov[1], gg);
					mcov[2]=_mm256_add_epi32(mcov[2], bb);
					mcov[3]=_mm256_add_epi32(mcov[3], rg);
					mcov[4]=_mm256_add_epi32(mcov[4], gb);
					mcov[5]=_mm256_add_epi32(mcov[5], br);
#endif
				}
				{
					ALIGN(32) int tmp[6*NCODERS];
					memcpy(tmp, mcov, sizeof(tmp));
					for(int k=0;k<6*NCODERS/2;++k)//lo
						sums[2*k+0]+=tmp[k];
					for(int k=0;k<6*NCODERS/2;++k)//hi
						sums[2*k+1]+=tmp[k+6*NCODERS/2];
				}
			}
#ifdef ONE_LOSSYCOV
			{
				int64_t cov[6]={0};
				for(int k=0;k<NCODERS;++k)
				{
					for(int k2=0;k2<6;++k2)
						cov[k2]+=sums[NCODERS*k2+k];
				}
				for(int k=0;k<NCODERS;++k)
				{
					for(int k2=0;k2<6;++k2)
						sums[NCODERS*k2+k]=cov[k2];
				}
			}
#endif
#ifdef _MSC_VER
			for(int k=0;k<NCODERS;++k)
				printf(
					"%10d %10d %10d\n"
					, lossymean[NCODERS*0+k]
					, lossymean[NCODERS*1+k]
					, lossymean[NCODERS*2+k]
				);
#endif
			for(int k=0;k<NCODERS;++k)
			{
#ifdef ONE_LOSSYCOV
				double invcount=1./((double)blockh*ixcount);
#else
				double invcount=1./((double)blockh*blockw);
#endif
				double cov[]=
				{
					sums[NCODERS*0+k]*invcount,
					sums[NCODERS*1+k]*invcount,
					sums[NCODERS*2+k]*invcount,
					sums[NCODERS*3+k]*invcount,
					sums[NCODERS*4+k]*invcount,
					sums[NCODERS*5+k]*invcount,
				};
				int16_t ltt[9];
				pca(cov, ltt);
#ifdef _MSC_VER
				printf(
					"%10.5lf %10.5lf %10.5lf\n"
					"%10.5lf %10.5lf %10.5lf\n"
					"%10.5lf %10.5lf %10.5lf\n\n"
					, (double)ltt[0]/0x4000
					, (double)ltt[1]/0x4000
					, (double)ltt[2]/0x4000
					, (double)ltt[3]/0x4000
					, (double)ltt[4]/0x4000
					, (double)ltt[5]/0x4000
					, (double)ltt[6]/0x4000
					, (double)ltt[7]/0x4000
					, (double)ltt[8]/0x4000
				);
#endif
				for(int k2=0;k2<9;++k2)
					lossytransform[NCODERS*k2+k]=ltt[k2];
			}
		}
#endif
		prof_checkpoint(usize, "analysis");
		switch(nthreads0&3)
		{
		case 1://use CG
			use_wg4=0;
			break;
		case 2://use WG4
			use_wg4=1;
			break;
		default://use L1 loss
			use_wg4=2;
			break;
		}
#ifndef LOUD
		if(nthreads0&4)
#endif
			printf("%s  WG4=%d  %td bytes  dist=%d\n", rct_names[bestrct], use_wg4, usize, dist);

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
	else
	{
		//decode flags, stats
		int flags=*streamptr++;
		use_wg4=flags&3;
		flags>>=2;
		bestrct=flags%RCT_COUNT;
		dist=*streamptr++;

		//if(dist>1)
		//{
		//	memcpy(lossymean	, streamptr, sizeof(lossymean		)); streamptr+=sizeof(lossymean		);
		//	memcpy(lossytransform	, streamptr, sizeof(lossytransform	)); streamptr+=sizeof(lossytransform	);
		//}

		ctxmask=*(unsigned long long*)streamptr;
		streamptr+=8;
		BitPackerLIFO ec;
		bitpacker_dec_init(&ec, streamptr, streamend);
		for(int kc=0;kc<3*NCTX;++kc)
			dec_unpackhist(&ec, CDF2syms+((ptrdiff_t)kc<<PROBBITS), ctxmask, kc);
		if(xremw||yremh)
		{
			dec_unpackhist(&ec, rCDF2syms+((ptrdiff_t)0<<PROBBITS), ctxmask, 3*NCTX+0);
			dec_unpackhist(&ec, rCDF2syms+((ptrdiff_t)1<<PROBBITS), ctxmask, 3*NCTX+1);
			dec_unpackhist(&ec, rCDF2syms+((ptrdiff_t)2<<PROBBITS), ctxmask, 3*NCTX+2);
		}
		streamptr=(unsigned char*)(size_t)ec.srcfwdptr;
		prof_checkpoint((ptrdiff_t)CDF2syms_size+rCDF2syms_size, "unpack histograms");

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
	}
	int L1statesize=0;
	int *L1state=0;
	if(use_wg4==1)
	{
		wgsize=(int)sizeof(short[2*NCODERS*3*WG_NPREDS])*(blockw+16);//2 padded rows  *  {WGY*NCODERS, WGU*NCODERS, WGV*NCODERS} * NPREDS = 3*32*8 = 768 channels  ~96*iw bytes
		wgerrors=(short*)_mm_malloc(wgsize, sizeof(__m256i));//~375 KB for 4K/12MP		NNEerrors = currerrors
		wgstatesize=(int)sizeof(short[2*NCODERS*3*WG_NPREDS]);//{preds, Wprederrors} * NPREDS * 3 channels * NCODERS
		wgstate=(short*)_mm_malloc(wgstatesize, sizeof(__m256i));
		if(!wgerrors||!wgstate)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(wgerrors, 0, wgsize);
		memset(wgstate, 0, wgstatesize);
	}
	else if(use_wg4==2)
	{
		L1statesize=(int)sizeof(int[2*NCODERS*3*(WG_NPREDS+1)]);//{preds, coeffs} * (NPREDS+{bias}) * 3 channels * NCODERS
		L1state=(int*)_mm_malloc(L1statesize, sizeof(__m256i));
		if(!L1state)
		{
			LOG_ERROR("Alloc error");
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
	//__m256i mlossymean[3];
	__m256i mlossytrans[9];
	if(dist>1)
	{
		dist_rcp=_mm256_set1_epi16(((1<<16)+dist-1)/dist);//x/dist  ->  {x*=inv; x=(x>>16)+((unsigned)x>>31);}
		mdist=_mm256_set1_epi16(dist);
		//for(int k=0;k<3;++k)
		//	mlossymean[k]=_mm256_set1_epi16(lossymean[k]);
#ifdef USE_YCBCR
		if(fwd)
		{
			mlossytrans[0]=_mm256_set1_epi16((int16_t)(+0.299	*0x4000+0.5));
			mlossytrans[1]=_mm256_set1_epi16((int16_t)(+0.587	*0x4000+0.5));
			mlossytrans[2]=_mm256_set1_epi16((int16_t)(+0.114	*0x4000+0.5));
			mlossytrans[3]=_mm256_set1_epi16((int16_t)(-0.168736	*0x4000+0.5));
			mlossytrans[4]=_mm256_set1_epi16((int16_t)(-0.331264	*0x4000+0.5));
			mlossytrans[5]=_mm256_set1_epi16((int16_t)(+0.5		*0x4000+0.5));
			mlossytrans[6]=_mm256_set1_epi16((int16_t)(+0.5		*0x4000+0.5));
			mlossytrans[7]=_mm256_set1_epi16((int16_t)(-0.418688	*0x4000+0.5));
			mlossytrans[8]=_mm256_set1_epi16((int16_t)(-0.081312	*0x4000+0.5));
		}
		else
		{
			mlossytrans[0]=_mm256_set1_epi16((int16_t)(+1		*0x4000+0.5));
			mlossytrans[1]=_mm256_set1_epi16((int16_t)(+0		*0x4000+0.5));
			mlossytrans[2]=_mm256_set1_epi16((int16_t)(+1.402	*0x4000+0.5));
			mlossytrans[3]=_mm256_set1_epi16((int16_t)(+1		*0x4000+0.5));
			mlossytrans[4]=_mm256_set1_epi16((int16_t)(-0.344136	*0x4000+0.5));
			mlossytrans[5]=_mm256_set1_epi16((int16_t)(-0.714136	*0x4000+0.5));
			mlossytrans[6]=_mm256_set1_epi16((int16_t)(+1		*0x4000+0.5));
			mlossytrans[7]=_mm256_set1_epi16((int16_t)(+1.772	*0x4000+0.5));
			mlossytrans[8]=_mm256_set1_epi16((int16_t)(+0		*0x4000+0.5));
		}
#else
		if(fwd)
		{
			for(int ky=0;ky<3;++ky)
				for(int kx=0;kx<3;++kx)
					mlossytrans[3*kx+ky]=_mm256_load_si256((__m256i*)lossytransform+3*ky+kx);
				//	mlossytrans[3*kx+ky]=_mm256_set1_epi16(lossytransform[3*ky+kx]);
		}
		else//transpose the transform
		{
			for(int k=0;k<9;++k)
				mlossytrans[k]=_mm256_load_si256((__m256i*)lossytransform+k);
			//	mlossytrans[k]=_mm256_set1_epi16(lossytransform[k]);
		}
#endif
	}
	memset(myuv, 0, sizeof(myuv));
	unsigned char *ctxptr=interleaved;
	imptr=interleaved+(fwd?isize:0);
	__m256i mstate[2];
	__m256i *wgpreds=use_wg4==2?(__m256i*)L1state:(__m256i*)wgstate;
	short *wgWerrors=use_wg4==2?(short*)(L1state+1*(ptrdiff_t)NCODERS*3*(WG_NPREDS+1)):wgstate+1*(ptrdiff_t)NCODERS*3*WG_NPREDS;
	if(use_wg4==2)
	{
		static const int weights0[]=
		{
			100000,//0	N
			100000,//1	W
			 80000,//2	3*(N-NN)+NNN
			 80000,//3	3*(W-WW)+WWW
			 50000,//4	W+NE-N
			 50000,//5	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4
			150000,//6	N+W-NW
			 50000,//7	N+NE-NNE
			0,
		};
		int *L1coeffs=(int*)wgWerrors;
		for(int k=0;k<WG_NPREDS;++k)
			FILLMEM(L1coeffs+6*8*k, weights0[k], sizeof(int[6*8]), sizeof(int));
	}
	//short *wgWerrors=wgstate+0*(ptrdiff_t)NCODERS*3*WG_NPREDS;
	//__m256i *wgpreds=(__m256i*)(wgstate+1*(ptrdiff_t)NCODERS*3*WG_NPREDS);
	//__m256i rcon=_mm256_set1_epi32(1<<L1sh>>1);
	//__m256i mL1sh=_mm256_set1_epi32(L1sh);
	if(!fwd)
	{
#ifdef _DEBUG
		if(streamptr>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		memcpy(mstate, streamptr, sizeof(mstate));
		streamptr+=sizeof(mstate);
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
		ALIGN(32) short *erows[]=
		{
			wgerrors+(paddedwidth*((ky-0LL)&1)+8LL)*NCODERS*3*WG_NPREDS,
			wgerrors+(paddedwidth*((ky-1LL)&1)+8LL)*NCODERS*3*WG_NPREDS,
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
		//memset(wgpreds, 0, sizeof(short[NCODERS*3*WG_NPREDS]));
		//                       (__m256i*)rows[-Y]+E+C+X*6
		eNEE[0]=_mm256_load_si256((__m256i*)rows[1]+3+0+2*6);
		eNEE[1]=_mm256_load_si256((__m256i*)rows[1]+3+1+2*6);
		eNEE[2]=_mm256_load_si256((__m256i*)rows[1]+3+2+2*6);
		for(int kx=0;kx<ixbytes;kx+=3*NCODERS)
		{
			//                     (__m256i*)rows[1]+E+C+X*6
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
				if(use_wg4)//predict
				{
					__m256i cache[3];
					/*
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
					wgpreds[0*3+0]=N[0];
					wgpreds[0*3+1]=N[1];
					wgpreds[0*3+2]=N[2];

					//W
					wgpreds[1*3+0]=W[0];
					wgpreds[1*3+1]=W[1];
					wgpreds[1*3+2]=W[2];

					//3*(N-NN)+NNN
					cache[0]=_mm256_sub_epi16(N[0], _mm256_load_si256((__m256i*)rows[2]+0+0+0*6));//N-NN
					cache[1]=_mm256_sub_epi16(N[1], _mm256_load_si256((__m256i*)rows[2]+0+1+0*6));
					cache[2]=_mm256_sub_epi16(N[2], _mm256_load_si256((__m256i*)rows[2]+0+2+0*6));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
					wgpreds[2*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[3]+0+0+0*6));//+NNN
					wgpreds[2*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[3]+0+1+0*6));
					wgpreds[2*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[3]+0+2+0*6));

					//3*(W-WW)+WWW
					cache[0]=_mm256_sub_epi16(W[0], _mm256_load_si256((__m256i*)rows[0]+0+0-2*6));//W-WW
					cache[1]=_mm256_sub_epi16(W[1], _mm256_load_si256((__m256i*)rows[0]+0+1-2*6));
					cache[2]=_mm256_sub_epi16(W[2], _mm256_load_si256((__m256i*)rows[0]+0+2-2*6));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
					wgpreds[3*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[0]+0+0-3*6));//+WWW
					wgpreds[3*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[0]+0+1-3*6));
					wgpreds[3*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[0]+0+2-3*6));

					//W+NE-N
					cache[0]=_mm256_sub_epi16(W[0], N[0]);
					cache[1]=_mm256_sub_epi16(W[1], N[1]);
					cache[2]=_mm256_sub_epi16(W[2], N[2]);
					wgpreds[4*3+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)rows[1]+0+0+1*6));//+NE
					wgpreds[4*3+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)rows[1]+0+1+1*6));
					wgpreds[4*3+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)rows[1]+0+2+1*6));

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
					wgpreds[5*3+0]=_mm256_srai_epi16(cache[0], 2);
					wgpreds[5*3+1]=_mm256_srai_epi16(cache[1], 2);
					wgpreds[5*3+2]=_mm256_srai_epi16(cache[2], 2);

					//N+W-NW
					wgpreds[6*3+0]=predY;
					wgpreds[6*3+1]=predU;
					wgpreds[6*3+2]=predV;

					//N+NE-NNE
					cache[0]=_mm256_add_epi16(N[0], _mm256_load_si256((__m256i*)rows[1]+0+0+1*6));//N+NE
					cache[1]=_mm256_add_epi16(N[1], _mm256_load_si256((__m256i*)rows[1]+0+1+1*6));
					cache[2]=_mm256_add_epi16(N[2], _mm256_load_si256((__m256i*)rows[1]+0+2+1*6));
					wgpreds[7*3+0]=_mm256_sub_epi16(cache[0], _mm256_load_si256((__m256i*)rows[2]+0+0+1*6));//NNE
					wgpreds[7*3+1]=_mm256_sub_epi16(cache[1], _mm256_load_si256((__m256i*)rows[2]+0+1+1*6));
					wgpreds[7*3+2]=_mm256_sub_epi16(cache[2], _mm256_load_si256((__m256i*)rows[2]+0+2+1*6));


					//mix
					if(use_wg4==2)
					{
						__m256i mp[6], t[6];
						int *L1coeffs=(int*)wgWerrors;
						mp[0]=_mm256_setzero_si256();
						mp[1]=_mm256_setzero_si256();
						mp[2]=_mm256_setzero_si256();
						mp[3]=_mm256_setzero_si256();
						mp[4]=_mm256_setzero_si256();
						mp[5]=_mm256_setzero_si256();

						//mp[0]=_mm256_load_si256((__m256i*)L1coeffs+WG_NPREDS*6+0);//bias
						//mp[1]=_mm256_load_si256((__m256i*)L1coeffs+WG_NPREDS*6+1);
						//mp[2]=_mm256_load_si256((__m256i*)L1coeffs+WG_NPREDS*6+2);
						//mp[3]=_mm256_load_si256((__m256i*)L1coeffs+WG_NPREDS*6+3);
						//mp[4]=_mm256_load_si256((__m256i*)L1coeffs+WG_NPREDS*6+4);
						//mp[5]=_mm256_load_si256((__m256i*)L1coeffs+WG_NPREDS*6+5);
#ifdef PRINT_L1_BOUNDS
						for(int k2=0;k2<6*8;++k2)
						{
							int bias=L1coeffs[WG_NPREDS*6*8+k2];
							if(bmin>bias)bmin=bias;
							if(bmax<bias)bmax=bias;
						}
#endif
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
						for(int k=0;k<WG_NPREDS;++k)
						{
							__m256i c[6];
							c[0]=_mm256_load_si256((__m256i*)L1coeffs+k*6+0);
							c[1]=_mm256_load_si256((__m256i*)L1coeffs+k*6+1);
							c[2]=_mm256_load_si256((__m256i*)L1coeffs+k*6+2);
							c[3]=_mm256_load_si256((__m256i*)L1coeffs+k*6+3);
							c[4]=_mm256_load_si256((__m256i*)L1coeffs+k*6+4);
							c[5]=_mm256_load_si256((__m256i*)L1coeffs+k*6+5);
							//16 -> 32		3 lo 3 hi registers
							t[0]=_mm256_slli_epi32(wgpreds[k*3+0], 16);
							t[1]=_mm256_slli_epi32(wgpreds[k*3+1], 16);
							t[2]=_mm256_slli_epi32(wgpreds[k*3+2], 16);
							t[3]=_mm256_srai_epi32(wgpreds[k*3+0], 16);
							t[4]=_mm256_srai_epi32(wgpreds[k*3+1], 16);
							t[5]=_mm256_srai_epi32(wgpreds[k*3+2], 16);
							t[0]=_mm256_srai_epi32(t[0], 16);
							t[1]=_mm256_srai_epi32(t[1], 16);
							t[2]=_mm256_srai_epi32(t[2], 16);
							t[0]=_mm256_mullo_epi32(t[0], c[0]);
							t[1]=_mm256_mullo_epi32(t[1], c[1]);
							t[2]=_mm256_mullo_epi32(t[2], c[2]);
							t[3]=_mm256_mullo_epi32(t[3], c[3]);
							t[4]=_mm256_mullo_epi32(t[4], c[4]);
							t[5]=_mm256_mullo_epi32(t[5], c[5]);
							mp[0]=_mm256_add_epi32(mp[0], t[0]);
							mp[1]=_mm256_add_epi32(mp[1], t[1]);
							mp[2]=_mm256_add_epi32(mp[2], t[2]);
							mp[3]=_mm256_add_epi32(mp[3], t[3]);
							mp[4]=_mm256_add_epi32(mp[4], t[4]);
							mp[5]=_mm256_add_epi32(mp[5], t[5]);
						}
#define L1SH 19
						__m256i rcon=_mm256_set1_epi32(1<<L1SH>>1);
						mp[0]=_mm256_add_epi32(mp[0], rcon);
						mp[3]=_mm256_add_epi32(mp[3], rcon);
					//	rcon=_mm256_slli_epi32(rcon, 1);
						mp[1]=_mm256_add_epi32(mp[1], rcon);
						mp[4]=_mm256_add_epi32(mp[4], rcon);
						mp[2]=_mm256_add_epi32(mp[2], rcon);
						mp[5]=_mm256_add_epi32(mp[5], rcon);

						mp[0]=_mm256_srai_epi32(mp[0], L1SH);
						mp[1]=_mm256_srai_epi32(mp[1], L1SH);
						mp[2]=_mm256_srai_epi32(mp[2], L1SH);
						mp[3]=_mm256_srai_epi32(mp[3], L1SH);
						mp[4]=_mm256_srai_epi32(mp[4], L1SH);
						mp[5]=_mm256_srai_epi32(mp[5], L1SH);
						//32 -> 16
						mp[0]=_mm256_slli_epi32(mp[0], 16);
						mp[1]=_mm256_slli_epi32(mp[1], 16);
						mp[2]=_mm256_slli_epi32(mp[2], 16);
						mp[3]=_mm256_slli_epi32(mp[3], 16);
						mp[4]=_mm256_slli_epi32(mp[4], 16);
						mp[5]=_mm256_slli_epi32(mp[5], 16);
						mp[0]=_mm256_srli_epi32(mp[0], 16);
						mp[1]=_mm256_srli_epi32(mp[1], 16);
						mp[2]=_mm256_srli_epi32(mp[2], 16);
						predY=_mm256_or_si256(mp[0], mp[3]);
						predU=_mm256_or_si256(mp[1], mp[4]);
						predV=_mm256_or_si256(mp[2], mp[5]);
					}
					else
					{
						wg_mix(0, wgWerrors, erows[0], erows[1], wgpreds, &predY);
						wg_mix(1, wgWerrors, erows[0], erows[1], wgpreds, &predU);
						wg_mix(2, wgWerrors, erows[0], erows[1], wgpreds, &predV);
					}


					//loosen pred range
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
				predYUV0[0]=predY;
				predYUV0[1]=predU;
				predYUV0[2]=predV;

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
			__m256i msyms, moffset;
			//if(ky==261&&kx==76800)//lane 3
			//	printf("");
			if(dist>1)
			{
				__m256i msyms0[3];
				if(fwd)
				{
					__m256i msyms1[3];
					__m256i ctxblendmask=_mm256_set1_epi16(255);
					__m256i mrgb[3];
					mrgb[0]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+0*NCODERS)), half8));//load yuv
					mrgb[1]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+1*NCODERS)), half8));
					mrgb[2]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+2*NCODERS)), half8));

					//forward transform
#ifndef USE_YCBCR
					mrgb[0]=_mm256_sub_epi16(mrgb[0], mlossymean[0]);
					mrgb[1]=_mm256_sub_epi16(mrgb[1], mlossymean[1]);
					mrgb[2]=_mm256_sub_epi16(mrgb[2], mlossymean[2]);
#endif
					mrgb[0]=_mm256_slli_epi16(mrgb[0], 4);
					mrgb[1]=_mm256_slli_epi16(mrgb[1], 4);
					mrgb[2]=_mm256_slli_epi16(mrgb[2], 4);
					myuv[0]=_mm256_mulhi_epi16(mlossytrans[0+0], mrgb[0]);
					myuv[1]=_mm256_mulhi_epi16(mlossytrans[3+0], mrgb[0]);
					myuv[2]=_mm256_mulhi_epi16(mlossytrans[6+0], mrgb[0]);
					myuv[0]=_mm256_add_epi16(myuv[0], _mm256_mulhi_epi16(mlossytrans[0+1], mrgb[1]));
					myuv[1]=_mm256_add_epi16(myuv[1], _mm256_mulhi_epi16(mlossytrans[3+1], mrgb[1]));
					myuv[2]=_mm256_add_epi16(myuv[2], _mm256_mulhi_epi16(mlossytrans[6+1], mrgb[1]));
					myuv[0]=_mm256_add_epi16(myuv[0], _mm256_mulhi_epi16(mlossytrans[0+2], mrgb[2]));
					myuv[1]=_mm256_add_epi16(myuv[1], _mm256_mulhi_epi16(mlossytrans[3+2], mrgb[2]));
					myuv[2]=_mm256_add_epi16(myuv[2], _mm256_mulhi_epi16(mlossytrans[6+2], mrgb[2]));
					myuv[0]=_mm256_srai_epi16(myuv[0], 2);
					myuv[1]=_mm256_srai_epi16(myuv[1], 2);
					myuv[2]=_mm256_srai_epi16(myuv[2], 2);

					msyms0[0]=_mm256_sub_epi16(myuv[0], predY);
					msyms0[1]=_mm256_sub_epi16(myuv[1], predU);
					msyms0[2]=_mm256_sub_epi16(myuv[2], predV);

					__m256i cond0=_mm256_srai_epi16(msyms0[0], 15);
					__m256i cond1=_mm256_srai_epi16(msyms0[1], 15);
					__m256i cond2=_mm256_srai_epi16(msyms0[2], 15);
					msyms0[0]=_mm256_mulhi_epi16(msyms0[0], dist_rcp);
					msyms0[1]=_mm256_mulhi_epi16(msyms0[1], dist_rcp);
					msyms0[2]=_mm256_mulhi_epi16(msyms0[2], dist_rcp);
					msyms0[0]=_mm256_sub_epi16(msyms0[0], cond0);//(x*inv>>16)-(x>>31)
					msyms0[1]=_mm256_sub_epi16(msyms0[1], cond1);
					msyms0[2]=_mm256_sub_epi16(msyms0[2], cond2);
					myuv[0]=_mm256_mullo_epi16(msyms0[0], mdist);
					myuv[1]=_mm256_mullo_epi16(msyms0[1], mdist);
					myuv[2]=_mm256_mullo_epi16(msyms0[2], mdist);
					myuv[0]=_mm256_add_epi16(myuv[0], predY);
					myuv[1]=_mm256_add_epi16(myuv[1], predU);
					myuv[2]=_mm256_add_epi16(myuv[2], predV);
					myuv[0]=_mm256_max_epi16(myuv[0], _mm256_set1_epi16(-128));
					myuv[1]=_mm256_max_epi16(myuv[1], _mm256_set1_epi16(-128));
					myuv[2]=_mm256_max_epi16(myuv[2], _mm256_set1_epi16(-128));
					myuv[0]=_mm256_min_epi16(myuv[0], _mm256_set1_epi16(127));
					myuv[1]=_mm256_min_epi16(myuv[1], _mm256_set1_epi16(127));
					myuv[2]=_mm256_min_epi16(myuv[2], _mm256_set1_epi16(127));

					msyms1[0]=_mm256_sub_epi16(msyms0[0], amin);
					msyms1[1]=_mm256_sub_epi16(msyms0[1], amin);
					msyms1[2]=_mm256_sub_epi16(msyms0[2], amin);
#ifdef _MSC_VER
					int mask=0;
					mask|=_mm256_movemask_epi8(_mm256_cmpgt_epi16(msyms1[0], _mm256_set1_epi16(255)));
					mask|=_mm256_movemask_epi8(_mm256_cmpgt_epi16(msyms1[1], _mm256_set1_epi16(255)));
					mask|=_mm256_movemask_epi8(_mm256_cmpgt_epi16(msyms1[2], _mm256_set1_epi16(255)));
					mask|=_mm256_movemask_epi8(_mm256_cmpgt_epi16(_mm256_setzero_si256(), msyms1[0]));
					mask|=_mm256_movemask_epi8(_mm256_cmpgt_epi16(_mm256_setzero_si256(), msyms1[1]));
					mask|=_mm256_movemask_epi8(_mm256_cmpgt_epi16(_mm256_setzero_si256(), msyms1[2]));
					if(mask)
						LOG_ERROR("");
#endif
					ctxU=_mm256_add_epi16(ctxU, mctxuoffset);
					ctxV=_mm256_add_epi16(ctxV, mctxvoffset);
					ctxY=_mm256_slli_epi16(ctxY, 8);
					ctxU=_mm256_slli_epi16(ctxU, 8);
					ctxV=_mm256_slli_epi16(ctxV, 8);
					ctxY=_mm256_blendv_epi8(ctxY, msyms1[0], ctxblendmask);
					ctxU=_mm256_blendv_epi8(ctxU, msyms1[1], ctxblendmask);
					ctxV=_mm256_blendv_epi8(ctxV, msyms1[2], ctxblendmask);
					_mm256_store_si256((__m256i*)syms+0, ctxY);
					_mm256_store_si256((__m256i*)syms+1, ctxU);
					_mm256_store_si256((__m256i*)syms+2, ctxV);
					_mm256_store_si256((__m256i*)ctxptr+0, ctxY);
					_mm256_store_si256((__m256i*)ctxptr+1, ctxU);
					_mm256_store_si256((__m256i*)ctxptr+2, ctxV);

#if defined _MSC_VER && 0
					for(int k=0;k<_countof(syms);++k)
					{
						if(syms[k]>=3*NCTX*256)
							LOG_ERROR("");
					}
#endif
					ctxptr+=sizeof(short[3][NCODERS]);

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
				}
				else
				{
					dec_yuv(mstate, &ctxY, (int*)CDF2syms+((ptrdiff_t)NCTX*0<<PROBBITS), ans_permute, &streamptr, streamend, myuv+0);//residuals from [0 ~ 255]
					dec_yuv(mstate, &ctxU, (int*)CDF2syms+((ptrdiff_t)NCTX*1<<PROBBITS), ans_permute, &streamptr, streamend, myuv+1);
					dec_yuv(mstate, &ctxV, (int*)CDF2syms+((ptrdiff_t)NCTX*2<<PROBBITS), ans_permute, &streamptr, streamend, myuv+2);

					myuv[0]=_mm256_sub_epi16(myuv[0], _mm256_set1_epi16(128));
					myuv[1]=_mm256_sub_epi16(myuv[1], _mm256_set1_epi16(128));
					myuv[2]=_mm256_sub_epi16(myuv[2], _mm256_set1_epi16(128));
					msyms0[0]=myuv[0];
					msyms0[1]=myuv[1];
					msyms0[2]=myuv[2];
					myuv[0]=_mm256_mullo_epi16(myuv[0], mdist);
					myuv[1]=_mm256_mullo_epi16(myuv[1], mdist);
					myuv[2]=_mm256_mullo_epi16(myuv[2], mdist);

					myuv[0]=_mm256_add_epi16(myuv[0], predY);
					myuv[1]=_mm256_add_epi16(myuv[1], predU);
					myuv[2]=_mm256_add_epi16(myuv[2], predV);
					myuv[0]=_mm256_max_epi16(myuv[0], _mm256_set1_epi16(-128));
					myuv[1]=_mm256_max_epi16(myuv[1], _mm256_set1_epi16(-128));
					myuv[2]=_mm256_max_epi16(myuv[2], _mm256_set1_epi16(-128));
					myuv[0]=_mm256_min_epi16(myuv[0], _mm256_set1_epi16(127));
					myuv[1]=_mm256_min_epi16(myuv[1], _mm256_set1_epi16(127));
					myuv[2]=_mm256_min_epi16(myuv[2], _mm256_set1_epi16(127));

					//inverse transform
					__m256i mrgb[3], tmp[3];
					tmp[0]=_mm256_slli_epi16(myuv[0], 4);
					tmp[1]=_mm256_slli_epi16(myuv[1], 4);
					tmp[2]=_mm256_slli_epi16(myuv[2], 4);
					mrgb[0]=_mm256_mulhi_epi16(mlossytrans[0+0], tmp[0]);
					mrgb[1]=_mm256_mulhi_epi16(mlossytrans[3+0], tmp[0]);
					mrgb[2]=_mm256_mulhi_epi16(mlossytrans[6+0], tmp[0]);
					mrgb[0]=_mm256_add_epi16(mrgb[0], _mm256_mulhi_epi16(mlossytrans[0+1], tmp[1]));
					mrgb[1]=_mm256_add_epi16(mrgb[1], _mm256_mulhi_epi16(mlossytrans[3+1], tmp[1]));
					mrgb[2]=_mm256_add_epi16(mrgb[2], _mm256_mulhi_epi16(mlossytrans[6+1], tmp[1]));
					mrgb[0]=_mm256_add_epi16(mrgb[0], _mm256_mulhi_epi16(mlossytrans[0+2], tmp[2]));
					mrgb[1]=_mm256_add_epi16(mrgb[1], _mm256_mulhi_epi16(mlossytrans[3+2], tmp[2]));
					mrgb[2]=_mm256_add_epi16(mrgb[2], _mm256_mulhi_epi16(mlossytrans[6+2], tmp[2]));
					mrgb[0]=_mm256_srai_epi16(mrgb[0], 2);
					mrgb[1]=_mm256_srai_epi16(mrgb[1], 2);
					mrgb[2]=_mm256_srai_epi16(mrgb[2], 2);
#ifndef USE_YCBCR
					mrgb[0]=_mm256_add_epi16(mrgb[0], mlossymean[0]);
					mrgb[1]=_mm256_add_epi16(mrgb[1], mlossymean[1]);
					mrgb[2]=_mm256_add_epi16(mrgb[2], mlossymean[2]);
#endif

					__m128i msyms8[3];
					msyms8[0]=_mm_packs_epi16(_mm256_extracti128_si256(mrgb[0], 0), _mm256_extracti128_si256(mrgb[0], 1));
					msyms8[1]=_mm_packs_epi16(_mm256_extracti128_si256(mrgb[1], 0), _mm256_extracti128_si256(mrgb[1], 1));
					msyms8[2]=_mm_packs_epi16(_mm256_extracti128_si256(mrgb[2], 0), _mm256_extracti128_si256(mrgb[2], 1));
					msyms8[0]=_mm_xor_si128(msyms8[0], _mm_set1_epi8(-128));
					msyms8[1]=_mm_xor_si128(msyms8[1], _mm_set1_epi8(-128));
					msyms8[2]=_mm_xor_si128(msyms8[2], _mm_set1_epi8(-128));
					_mm_store_si128((__m128i*)(imptr+0*NCODERS), msyms8[0]);//store RGB bytes
					_mm_store_si128((__m128i*)(imptr+1*NCODERS), msyms8[1]);
					_mm_store_si128((__m128i*)(imptr+2*NCODERS), msyms8[2]);
#if defined ENABLE_GUIDE && 1
					for(int kc=0;kc<3;++kc)
					{
						for(int k=0;k<NCODERS;++k)
						{
							double diff=imptr[kc*NCODERS+k]-g_image[imptr-interleaved+kc*NCODERS+k];
							g_sqe[kc]+=diff*diff;
						}
					}
#endif
				}
				W[0]=myuv[0];
				W[1]=myuv[1];
				W[2]=myuv[2];
				_mm256_store_si256((__m256i*)rows[0]+0+0+0*6, myuv[0]);//store neighbors
				_mm256_store_si256((__m256i*)rows[0]+0+1+0*6, myuv[1]);
				_mm256_store_si256((__m256i*)rows[0]+0+2+0*6, myuv[2]);
				ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(msyms0[0], 1), _mm256_srai_epi16(msyms0[0], 15));
				ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms0[1], 1), _mm256_srai_epi16(msyms0[1], 15));
				ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms0[2], 1), _mm256_srai_epi16(msyms0[2], 15));
			}
			else if(fwd)
			{
				__m256i ctxblendmask=_mm256_set1_epi16(255);
				myuv[0]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+yidx)), half8));//load yuv
				myuv[1]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+uidx)), half8));
				myuv[2]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+vidx)), half8));

				//encode Y
				msyms=_mm256_sub_epi16(myuv[0], predY);//sub pred
#if 0
				if(dist>1)
				{
					__m256i cond=_mm256_srai_epi16(msyms, 15);
					msyms=_mm256_mulhi_epi16(msyms, dist_rcp);
					msyms=_mm256_sub_epi16(msyms, cond);//(x*inv>>16)-(x>>31)
					myuv[0]=_mm256_mullo_epi16(msyms, mdist);
					myuv[0]=_mm256_add_epi16(myuv[0], predY);
					myuv[0]=_mm256_max_epi16(myuv[0], _mm256_set1_epi16(-128));
					myuv[0]=_mm256_min_epi16(myuv[0], _mm256_set1_epi16(127));
				}
#endif
				W[0]=myuv[0];
				_mm256_store_si256((__m256i*)rows[0]+0+0+0*6, myuv[0]);//store Y neighbors
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
#if 0
				if(dist>1)
				{
					__m256i cond=_mm256_srai_epi16(msyms, 15);
					msyms=_mm256_mulhi_epi16(msyms, dist_rcp);
					msyms=_mm256_sub_epi16(msyms, cond);//(x*inv>>16)-(x>>31)
					myuv[1]=_mm256_mullo_epi16(msyms, mdist);
					myuv[1]=_mm256_add_epi16(myuv[1], predU);
					myuv[1]=_mm256_max_epi16(myuv[1], _mm256_set1_epi16(-128));
					myuv[1]=_mm256_min_epi16(myuv[1], _mm256_set1_epi16(127));
				}
#endif
				ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
				msyms=_mm256_sub_epi16(msyms, amin);
				ctxU=_mm256_add_epi16(ctxU, mctxuoffset);
				ctxU=_mm256_slli_epi16(ctxU, 8);
				ctxU=_mm256_blendv_epi8(ctxU, msyms, ctxblendmask);
				_mm256_store_si256((__m256i*)syms+1, ctxU);
				_mm256_store_si256((__m256i*)ctxptr+1, ctxU);//store U  ctx|residuals
#ifdef VNATIVE
				W[1]=_mm256_sub_epi16(myuv[1], moffset);
				_mm256_store_si256((__m256i*)rows[0]+0+1+0*6, W[1]);//store U neighbors
#else
				msyms=_mm256_sub_epi16(myuv[1], moffset);
				W[1]=msyms;
				_mm256_store_si256((__m256i*)rows[0]+0+1+0*6, msyms);//store U neighbors
#endif
				
				//encode V
				moffset=_mm256_mullo_epi16(vc0, myuv[0]);
				moffset=_mm256_add_epi16(moffset, _mm256_mullo_epi16(vc1, myuv[1]));
#ifdef VNATIVE
				moffset=_mm256_srai_epi16(moffset, 2);
#endif
				predV=_mm256_add_epi16(predV, moffset);
#ifndef VNATIVE
				predV=_mm256_srai_epi16(predV, 2);
#endif
				predV=_mm256_max_epi16(predV, amin);
				predV=_mm256_min_epi16(predV, amax);

				msyms=_mm256_sub_epi16(myuv[2], predV);
#if 0
				if(dist>1)
				{
					__m256i cond=_mm256_srai_epi16(msyms, 15);
					msyms=_mm256_mulhi_epi16(msyms, dist_rcp);
					msyms=_mm256_sub_epi16(msyms, cond);//(x*inv>>16)-(x>>31)
					myuv[2]=_mm256_mullo_epi16(msyms, mdist);
					myuv[2]=_mm256_add_epi16(myuv[2], predV);
					myuv[2]=_mm256_max_epi16(myuv[2], _mm256_set1_epi16(-128));
					myuv[2]=_mm256_min_epi16(myuv[2], _mm256_set1_epi16(127));
#ifdef _DEBUG
					for(int k=0;k<48;++k)
					{
						int v=((short*)myuv)[k];
						if((unsigned)(v+128)>=256)
							LOG_ERROR("");
					}
#endif
				}
#endif
				ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms, 1), _mm256_srai_epi16(msyms, 15));
				msyms=_mm256_sub_epi16(msyms, amin);
				ctxV=_mm256_add_epi16(ctxV, mctxvoffset);
				ctxV=_mm256_slli_epi16(ctxV, 8);
				ctxV=_mm256_blendv_epi8(ctxV, msyms, ctxblendmask);
				_mm256_store_si256((__m256i*)syms+2, ctxV);
				_mm256_store_si256((__m256i*)ctxptr+2, ctxV);//store V  ctx|residuals		ctxptr+NCODERS*(C*2+R)
#ifdef VNATIVE
				W[2]=_mm256_sub_epi16(myuv[2], moffset);
				_mm256_store_si256((__m256i*)rows[0]+0+2+0*6, W[2]);//store V neighbors
#else
				msyms=_mm256_slli_epi16(myuv[2], 2);
				msyms=_mm256_sub_epi16(msyms, moffset);
				W[2]=msyms;
				_mm256_store_si256((__m256i*)rows[0]+0+2+0*6, msyms);//store V neighbors
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
				//if(ky==1&&kx==1296)//
				//if(ky==0&&kx==48)//
				//	printf("");

				//decode main
				__m256i msyms0[3];
				__m128i msyms8;
				
				//decode Y
				dec_yuv(mstate, &ctxY, (int*)CDF2syms+((ptrdiff_t)NCTX*0<<PROBBITS), ans_permute, &streamptr, streamend, myuv+0);//residuals from [0 ~ 255]
				dec_yuv(mstate, &ctxU, (int*)CDF2syms+((ptrdiff_t)NCTX*1<<PROBBITS), ans_permute, &streamptr, streamend, myuv+1);
				dec_yuv(mstate, &ctxV, (int*)CDF2syms+((ptrdiff_t)NCTX*2<<PROBBITS), ans_permute, &streamptr, streamend, myuv+2);
#ifdef _DEBUG
				if(streamptr>streamend)
					LOG_ERROR("OOB stream %8.4lf%%  image %8.4lf%%"
						, 100.*(streamptr-streamstart)/(streamend-streamstart)
						, 100.*(ixbytes*ky+kx)/(ixbytes*blockh)
					);
#endif
#if 0
				if(dist>1)
				{
					__m256i tmp=_mm256_set1_epi16(128);
					myuv[0]=_mm256_sub_epi16(myuv[0], tmp);
					myuv[1]=_mm256_sub_epi16(myuv[1], tmp);
					myuv[2]=_mm256_sub_epi16(myuv[2], tmp);
					msyms0[0]=myuv[0];
					msyms0[1]=myuv[1];
					msyms0[2]=myuv[2];
					myuv[0]=_mm256_mullo_epi16(myuv[0], mdist);
					myuv[1]=_mm256_mullo_epi16(myuv[1], mdist);
					myuv[2]=_mm256_mullo_epi16(myuv[2], mdist);
					myuv[0]=_mm256_add_epi16(myuv[0], tmp);
					myuv[1]=_mm256_add_epi16(myuv[1], tmp);
					myuv[2]=_mm256_add_epi16(myuv[2], tmp);
				}
#endif

				//yuv = (char)(sym+pred-128)	= (unsigned char)(sym+pred)-128
				myuv[0]=_mm256_add_epi16(myuv[0], predY);
				myuv[0]=_mm256_and_si256(myuv[0], bytemask);
				myuv[0]=_mm256_add_epi16(myuv[0], amin);
				//msyms=_mm256_sub_epi16(myuv[0], predY);//sub pred
				msyms8=_mm_packs_epi16(_mm256_extracti128_si256(myuv[0], 0), _mm256_extracti128_si256(myuv[0], 1));
				msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
				_mm_store_si128((__m128i*)(imptr+yidx), msyms8);//store Y bytes
#ifdef ENABLE_GUIDE
				if(dist<=1&&memcmp(imptr+yidx, g_image+(imptr-interleaved)+yidx, NCODERS))
				{
					printf("original  decoded  original-decoded  XYC0 %d %d %d\n", kx, ky, yidx);
					for(int k=0;k<NCODERS;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+yidx+k],
							imptr[yidx+k],
							g_image[imptr-interleaved+yidx+k]-imptr[yidx+k]
						);
					LOG_ERROR("guide error XYC0 %d %d %d/%d", kx, ky, yidx, NCODERS);
				}
#endif
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
				//msyms=_mm256_sub_epi16(myuv[1], predU);//sub pred
				msyms8=_mm_packs_epi16(_mm256_extracti128_si256(myuv[1], 0), _mm256_extracti128_si256(myuv[1], 1));
				msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
				_mm_store_si128((__m128i*)(imptr+uidx), msyms8);//store U bytes
#ifdef ENABLE_GUIDE
				if(dist<=1&&memcmp(imptr+uidx, g_image+(imptr-interleaved)+uidx, NCODERS))
				{
					printf("original  decoded  original-decoded  XYC1 %d %d %d\n", kx, ky, uidx);
					for(int k=0;k<NCODERS;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+uidx+k],
							imptr[uidx+k],
							g_image[imptr-interleaved+uidx+k]-imptr[uidx+k]
						);
					LOG_ERROR("guide error XYC1 %d %d %d/%d", kx, ky, uidx, NCODERS);
				}
#endif
				W[1]=_mm256_sub_epi16(myuv[1], moffset);//subtract Uoffset from U
				_mm256_store_si256((__m256i*)rows[0]+0+1+0*6, W[1]);//store U neighbors
				

				//decode V
				moffset=_mm256_mullo_epi16(vc0, myuv[0]);
				moffset=_mm256_add_epi16(moffset, _mm256_mullo_epi16(vc1, myuv[1]));
#ifdef VNATIVE
				moffset=_mm256_srai_epi16(moffset, 2);
				predV=_mm256_add_epi16(predV, moffset);
#else
				predV=_mm256_add_epi16(predV, moffset);
				predV=_mm256_srai_epi16(predV, 2);
#endif
				predV=_mm256_max_epi16(predV, amin);
				predV=_mm256_min_epi16(predV, amax);
				
				myuv[2]=_mm256_add_epi16(myuv[2], predV);
				myuv[2]=_mm256_and_si256(myuv[2], bytemask);
				myuv[2]=_mm256_add_epi16(myuv[2], amin);
				//msyms=_mm256_sub_epi16(myuv[2], predV);
				msyms8=_mm_packs_epi16(_mm256_extracti128_si256(myuv[2], 0), _mm256_extracti128_si256(myuv[2], 1));
				msyms8=_mm_xor_si128(msyms8, _mm_set1_epi8(-128));
				_mm_store_si128((__m128i*)(imptr+vidx), msyms8);//store V bytes
#ifdef ENABLE_GUIDE
				if(dist<=1&&memcmp(imptr+vidx, g_image+(imptr-interleaved)+vidx, NCODERS))
				{
					printf("original  decoded  original-decoded  XYC2 %d %d %d\n", kx, ky, vidx);
					for(int k=0;k<NCODERS;++k)
						printf("0x%02X  0x%02X  %4d\n",
							g_image[imptr-interleaved+vidx+k],
							imptr[vidx+k],
							g_image[imptr-interleaved+vidx+k]-imptr[vidx+k]
						);
					LOG_ERROR("guide error XYC2 %d %d %d/%d", kx, ky, vidx, NCODERS);
				}
#endif
#ifdef VNATIVE
				W[2]=_mm256_sub_epi16(myuv[2], moffset);//subtract Voffset from V
#else
				W[2]=_mm256_slli_epi16(myuv[2], 2);
				W[2]=_mm256_sub_epi16(W[2], moffset);//subtract Voffset from V
#endif
				_mm256_store_si256((__m256i*)rows[0]+0+2+0*6, W[2]);//store V neighbors
				
#if defined ENABLE_GUIDE && 0
				if(dist>1)
				{
					for(int kc=0;kc<3;++kc)
					{
						for(int k=0;k<NCODERS;++k)
						{
							double diff=imptr[kc*NCODERS+k]-g_image[imptr-interleaved+kc*NCODERS+k];
							g_sqe[kc]+=diff*diff;
						}
					}
				}
#endif
//#ifdef ENABLE_GUIDE
//				for(int kx2=0;kx2<3*NCODERS;kx2+=3)
//					guide_check(image, kx+kx2, ky);
//#endif
				if(dist<=1)
				{
					msyms0[0]=_mm256_sub_epi16(myuv[0], predY);
					msyms0[1]=_mm256_sub_epi16(myuv[1], predU);
					msyms0[2]=_mm256_sub_epi16(myuv[2], predV);
				}
				ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(msyms0[0], 1), _mm256_srai_epi16(msyms0[0], 15));
				ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms0[1], 1), _mm256_srai_epi16(msyms0[1], 15));
				ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms0[2], 1), _mm256_srai_epi16(msyms0[2], 15));
			}
			if(use_wg4==1)//update
			{
				for(int kr=0;kr<3;++kr)
				{
					__m256i curr=W[kr];
					__m256i me0=_mm256_sub_epi16(curr, wgpreds[0*3+kr]);
					__m256i me1=_mm256_sub_epi16(curr, wgpreds[1*3+kr]);
					__m256i me2=_mm256_sub_epi16(curr, wgpreds[2*3+kr]);
					__m256i me3=_mm256_sub_epi16(curr, wgpreds[3*3+kr]);
					__m256i me4=_mm256_sub_epi16(curr, wgpreds[4*3+kr]);
					__m256i me5=_mm256_sub_epi16(curr, wgpreds[5*3+kr]);
					__m256i me6=_mm256_sub_epi16(curr, wgpreds[6*3+kr]);
					__m256i me7=_mm256_sub_epi16(curr, wgpreds[7*3+kr]);
					me0=_mm256_abs_epi16(me0);
					me1=_mm256_abs_epi16(me1);
					me2=_mm256_abs_epi16(me2);
					me3=_mm256_abs_epi16(me3);
					me4=_mm256_abs_epi16(me4);
					me5=_mm256_abs_epi16(me5);
					me6=_mm256_abs_epi16(me6);
					me7=_mm256_abs_epi16(me7);
					//me0=_mm256_xor_si256(_mm256_slli_epi16(me0, 1), _mm256_srai_epi16(me0, 15));
					//me1=_mm256_xor_si256(_mm256_slli_epi16(me1, 1), _mm256_srai_epi16(me1, 15));
					//me2=_mm256_xor_si256(_mm256_slli_epi16(me2, 1), _mm256_srai_epi16(me2, 15));
					//me3=_mm256_xor_si256(_mm256_slli_epi16(me3, 1), _mm256_srai_epi16(me3, 15));
					//me4=_mm256_xor_si256(_mm256_slli_epi16(me4, 1), _mm256_srai_epi16(me4, 15));
					//me5=_mm256_xor_si256(_mm256_slli_epi16(me5, 1), _mm256_srai_epi16(me5, 15));
					//me6=_mm256_xor_si256(_mm256_slli_epi16(me6, 1), _mm256_srai_epi16(me6, 15));
					//me7=_mm256_xor_si256(_mm256_slli_epi16(me7, 1), _mm256_srai_epi16(me7, 15));
					if(kx)
					{
						__m256i best=_mm256_min_epi16(me0, me1);
						best=_mm256_min_epi16(best, me2);
						best=_mm256_min_epi16(best, me3);
						best=_mm256_min_epi16(best, me4);
						best=_mm256_min_epi16(best, me5);
						best=_mm256_min_epi16(best, me6);
						best=_mm256_min_epi16(best, me7);
						me0=_mm256_sub_epi16(me0, best);
						me1=_mm256_sub_epi16(me1, best);
						me2=_mm256_sub_epi16(me2, best);
						me3=_mm256_sub_epi16(me3, best);
						me4=_mm256_sub_epi16(me4, best);
						me5=_mm256_sub_epi16(me5, best);
						me6=_mm256_sub_epi16(me6, best);
						me7=_mm256_sub_epi16(me7, best);
					}
					
					//                 (__m256i*)erows[-Y]+(P+X*NPREDS)*3+C
					_mm256_store_si256((__m256i*)erows[0]+(0+0*WG_NPREDS)*3+kr, me0);//ecurr
					_mm256_store_si256((__m256i*)erows[0]+(1+0*WG_NPREDS)*3+kr, me1);
					_mm256_store_si256((__m256i*)erows[0]+(2+0*WG_NPREDS)*3+kr, me2);
					_mm256_store_si256((__m256i*)erows[0]+(3+0*WG_NPREDS)*3+kr, me3);
					_mm256_store_si256((__m256i*)erows[0]+(4+0*WG_NPREDS)*3+kr, me4);
					_mm256_store_si256((__m256i*)erows[0]+(5+0*WG_NPREDS)*3+kr, me5);
					_mm256_store_si256((__m256i*)erows[0]+(6+0*WG_NPREDS)*3+kr, me6);
					_mm256_store_si256((__m256i*)erows[0]+(7+0*WG_NPREDS)*3+kr, me7);
					//__m256i t[8];
					//t[0]=_mm256_add_epi16(me0, _mm256_load_si256((__m256i*)erows[1]+(0+1*WG_NPREDS)*3+kr));//eNE+=e
					//t[1]=_mm256_add_epi16(me1, _mm256_load_si256((__m256i*)erows[1]+(1+1*WG_NPREDS)*3+kr));
					//t[2]=_mm256_add_epi16(me2, _mm256_load_si256((__m256i*)erows[1]+(2+1*WG_NPREDS)*3+kr));
					//t[3]=_mm256_add_epi16(me3, _mm256_load_si256((__m256i*)erows[1]+(3+1*WG_NPREDS)*3+kr));
					//t[4]=_mm256_add_epi16(me4, _mm256_load_si256((__m256i*)erows[1]+(4+1*WG_NPREDS)*3+kr));
					//t[5]=_mm256_add_epi16(me5, _mm256_load_si256((__m256i*)erows[1]+(5+1*WG_NPREDS)*3+kr));
					//t[6]=_mm256_add_epi16(me6, _mm256_load_si256((__m256i*)erows[1]+(6+1*WG_NPREDS)*3+kr));
					//t[7]=_mm256_add_epi16(me7, _mm256_load_si256((__m256i*)erows[1]+(7+1*WG_NPREDS)*3+kr));
					//_mm256_store_si256((__m256i*)erows[1]+(0+1*WG_NPREDS)*3+kr, t[0]);
					//_mm256_store_si256((__m256i*)erows[1]+(1+1*WG_NPREDS)*3+kr, t[1]);
					//_mm256_store_si256((__m256i*)erows[1]+(2+1*WG_NPREDS)*3+kr, t[2]);
					//_mm256_store_si256((__m256i*)erows[1]+(3+1*WG_NPREDS)*3+kr, t[3]);
					//_mm256_store_si256((__m256i*)erows[1]+(4+1*WG_NPREDS)*3+kr, t[4]);
					//_mm256_store_si256((__m256i*)erows[1]+(5+1*WG_NPREDS)*3+kr, t[5]);
					//_mm256_store_si256((__m256i*)erows[1]+(6+1*WG_NPREDS)*3+kr, t[6]);
					//_mm256_store_si256((__m256i*)erows[1]+(7+1*WG_NPREDS)*3+kr, t[7]);

#ifdef WG_ENABLE_eW
					__m256i t2[8];
					t2[0]=_mm256_load_si256((__m256i*)wgWerrors+0*3+kr);
					t2[1]=_mm256_load_si256((__m256i*)wgWerrors+1*3+kr);
					t2[2]=_mm256_load_si256((__m256i*)wgWerrors+2*3+kr);
					t2[3]=_mm256_load_si256((__m256i*)wgWerrors+3*3+kr);
					t2[4]=_mm256_load_si256((__m256i*)wgWerrors+4*3+kr);
					t2[5]=_mm256_load_si256((__m256i*)wgWerrors+5*3+kr);
					t2[6]=_mm256_load_si256((__m256i*)wgWerrors+6*3+kr);
					t2[7]=_mm256_load_si256((__m256i*)wgWerrors+7*3+kr);
					me0=_mm256_slli_epi16(me0, 5);
					me1=_mm256_slli_epi16(me1, 5);
					me2=_mm256_slli_epi16(me2, 5);
					me3=_mm256_slli_epi16(me3, 5);
					me4=_mm256_slli_epi16(me4, 5);
					me5=_mm256_slli_epi16(me5, 5);
					me6=_mm256_slli_epi16(me6, 5);
					me7=_mm256_slli_epi16(me7, 5);
					me0=_mm256_sub_epi16(me0, t2[0]);
					me1=_mm256_sub_epi16(me1, t2[1]);
					me2=_mm256_sub_epi16(me2, t2[2]);
					me3=_mm256_sub_epi16(me3, t2[3]);
					me4=_mm256_sub_epi16(me4, t2[4]);
					me5=_mm256_sub_epi16(me5, t2[5]);
					me6=_mm256_sub_epi16(me6, t2[6]);
					me7=_mm256_sub_epi16(me7, t2[7]);
					__m256i rcon=_mm256_set1_epi16(1<<3>>1);
					me0=_mm256_add_epi16(me0, rcon);
					me1=_mm256_add_epi16(me1, rcon);
					me2=_mm256_add_epi16(me2, rcon);
					me3=_mm256_add_epi16(me3, rcon);
					me4=_mm256_add_epi16(me4, rcon);
					me5=_mm256_add_epi16(me5, rcon);
					me6=_mm256_add_epi16(me6, rcon);
					me7=_mm256_add_epi16(me7, rcon);
					me0=_mm256_srai_epi16(me0, 3);
					me1=_mm256_srai_epi16(me1, 3);
					me2=_mm256_srai_epi16(me2, 3);
					me3=_mm256_srai_epi16(me3, 3);
					me4=_mm256_srai_epi16(me4, 3);
					me5=_mm256_srai_epi16(me5, 3);
					me6=_mm256_srai_epi16(me6, 3);
					me7=_mm256_srai_epi16(me7, 3);
					t2[0]=_mm256_add_epi16(t2[0], me0);
					t2[1]=_mm256_add_epi16(t2[1], me1);
					t2[2]=_mm256_add_epi16(t2[2], me2);
					t2[3]=_mm256_add_epi16(t2[3], me3);
					t2[4]=_mm256_add_epi16(t2[4], me4);
					t2[5]=_mm256_add_epi16(t2[5], me5);
					t2[6]=_mm256_add_epi16(t2[6], me6);
					t2[7]=_mm256_add_epi16(t2[7], me7);
					_mm256_store_si256((__m256i*)wgWerrors+0*3+kr, t2[0]);
					_mm256_store_si256((__m256i*)wgWerrors+1*3+kr, t2[1]);
					_mm256_store_si256((__m256i*)wgWerrors+2*3+kr, t2[2]);
					_mm256_store_si256((__m256i*)wgWerrors+3*3+kr, t2[3]);
					_mm256_store_si256((__m256i*)wgWerrors+4*3+kr, t2[4]);
					_mm256_store_si256((__m256i*)wgWerrors+5*3+kr, t2[5]);
					_mm256_store_si256((__m256i*)wgWerrors+6*3+kr, t2[6]);
					_mm256_store_si256((__m256i*)wgWerrors+7*3+kr, t2[7]);
#endif
				}
				erows[0]+=NCODERS*3*WG_NPREDS;
				erows[1]+=NCODERS*3*WG_NPREDS;
			}
			else if(use_wg4==2)//update
			{
				__m256i mu[3];

				mu[0]=_mm256_sub_epi16(W[0], predYUV0[0]);
				mu[1]=_mm256_sub_epi16(W[1], predYUV0[1]);
				mu[2]=_mm256_sub_epi16(W[2], predYUV0[2]);
#if 0
				{
					__m256i t0, t1, t2;

					t0=_mm256_cmpgt_epi16(W[0], predYUV0[0]);
					t1=_mm256_cmpgt_epi16(W[1], predYUV0[1]);
					t2=_mm256_cmpgt_epi16(W[2], predYUV0[2]);
					t0=_mm256_srli_epi16(t0, 15);
					t1=_mm256_srli_epi16(t1, 15);
					t2=_mm256_srli_epi16(t2, 15);
					mu[0]=_mm256_cmpgt_epi16(predYUV0[0], W[0]);
					mu[1]=_mm256_cmpgt_epi16(predYUV0[1], W[1]);
					mu[2]=_mm256_cmpgt_epi16(predYUV0[2], W[2]);
					mu[0]=_mm256_add_epi16(mu[0], t0);
					mu[1]=_mm256_add_epi16(mu[1], t1);
					mu[2]=_mm256_add_epi16(mu[2], t2);
				}
#endif
				int *L1coeffs=(int*)wgWerrors;

				//if(wgpreds[0*3+0].m256i_i16[0])//
				//	printf("");

#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 8
#endif
				for(int k=0;k<WG_NPREDS;++k)//update
				{
#if 0
					for(int k2=0;k2<16*3;++k2)//
					{
						int prod=((short*)wgpreds)[k2]*((short*)mu)[k2];
						int hi=prod>>16;
						if(hi!=0&&hi!=-1)
							LOG_ERROR("");
					}
#endif
					__m256i mc[6];
					mc[0]=_mm256_sign_epi16(wgpreds[k*3+0], mu[0]);
					mc[1]=_mm256_sign_epi16(wgpreds[k*3+1], mu[1]);
					mc[2]=_mm256_sign_epi16(wgpreds[k*3+2], mu[2]);
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
					mc[0]=_mm256_add_epi32(mc[0], _mm256_load_si256((__m256i*)L1coeffs+k*6+0));//update coeffs
					mc[1]=_mm256_add_epi32(mc[1], _mm256_load_si256((__m256i*)L1coeffs+k*6+1));
					mc[2]=_mm256_add_epi32(mc[2], _mm256_load_si256((__m256i*)L1coeffs+k*6+2));
					mc[3]=_mm256_add_epi32(mc[3], _mm256_load_si256((__m256i*)L1coeffs+k*6+3));
					mc[4]=_mm256_add_epi32(mc[4], _mm256_load_si256((__m256i*)L1coeffs+k*6+4));
					mc[5]=_mm256_add_epi32(mc[5], _mm256_load_si256((__m256i*)L1coeffs+k*6+5));
					_mm256_store_si256((__m256i*)L1coeffs+k*6+0, mc[0]);
					_mm256_store_si256((__m256i*)L1coeffs+k*6+1, mc[1]);
					_mm256_store_si256((__m256i*)L1coeffs+k*6+2, mc[2]);
					_mm256_store_si256((__m256i*)L1coeffs+k*6+3, mc[3]);
					_mm256_store_si256((__m256i*)L1coeffs+k*6+4, mc[4]);
					_mm256_store_si256((__m256i*)L1coeffs+k*6+5, mc[5]);
				}
				//update bias	16 -> 32	3 lo 3 hi registers
			//	mu[3]=_mm256_srai_epi32(mu[0], 16);
			//	mu[4]=_mm256_srai_epi32(mu[1], 16);
			//	mu[5]=_mm256_srai_epi32(mu[2], 16);
			//	mu[0]=_mm256_slli_epi32(mu[0], 16);
			//	mu[1]=_mm256_slli_epi32(mu[1], 16);
			//	mu[2]=_mm256_slli_epi32(mu[2], 16);
			//	mu[0]=_mm256_srai_epi32(mu[0], 16);
			//	mu[1]=_mm256_srai_epi32(mu[1], 16);
			//	mu[2]=_mm256_srai_epi32(mu[2], 16);
			//	mu[0]=_mm256_add_epi32(mu[0], _mm256_load_si256((__m256i*)L1coeffs+WG_NPREDS*6+0));
			//	mu[1]=_mm256_add_epi32(mu[1], _mm256_load_si256((__m256i*)L1coeffs+WG_NPREDS*6+1));
			//	mu[2]=_mm256_add_epi32(mu[2], _mm256_load_si256((__m256i*)L1coeffs+WG_NPREDS*6+2));
			//	mu[3]=_mm256_add_epi32(mu[3], _mm256_load_si256((__m256i*)L1coeffs+WG_NPREDS*6+3));
			//	mu[4]=_mm256_add_epi32(mu[4], _mm256_load_si256((__m256i*)L1coeffs+WG_NPREDS*6+4));
			//	mu[5]=_mm256_add_epi32(mu[5], _mm256_load_si256((__m256i*)L1coeffs+WG_NPREDS*6+5));
			//	_mm256_store_si256((__m256i*)L1coeffs+WG_NPREDS*6+0, mu[0]);
			//	_mm256_store_si256((__m256i*)L1coeffs+WG_NPREDS*6+1, mu[1]);
			//	_mm256_store_si256((__m256i*)L1coeffs+WG_NPREDS*6+2, mu[2]);
			//	_mm256_store_si256((__m256i*)L1coeffs+WG_NPREDS*6+3, mu[3]);
			//	_mm256_store_si256((__m256i*)L1coeffs+WG_NPREDS*6+4, mu[4]);
			//	_mm256_store_si256((__m256i*)L1coeffs+WG_NPREDS*6+5, mu[5]);
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
	if(use_wg4==1)
	{
		_mm_free(wgerrors);
		_mm_free(wgstate);
	}
	else if(use_wg4==2)
	{
		_mm_free(L1state);
	}
	if(fwd)//all rANS encoding is bwd-bwd
	{
		rANS_SIMD_SymInfo *syminfo=(rANS_SIMD_SymInfo*)CDF2syms;
		rANS_SIMD_SymInfo *rsyminfo=(rANS_SIMD_SymInfo*)rCDF2syms;

		//normalize/integrate hists
		for(int kc=0;kc<3*NCTX;++kc)
			enc_hist2stats(hists+(ptrdiff_t)256*kc, syminfo+(ptrdiff_t)256*kc, &ctxmask, kc);
		
		//encode remainder
		if(xremw||yremh)
		{
			memset(rhist, 0, rhsize);
			for(int ky=0;ky<yremh;++ky)
				decorr1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, bestrct, rhist);
			for(int kx=0;kx<xremw;++kx)
				decorr1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, bestrct, rhist);
			enc_hist2stats(rhist+(ptrdiff_t)256*0, rsyminfo+(ptrdiff_t)256*0, &ctxmask, 3*NCTX+0);
			enc_hist2stats(rhist+(ptrdiff_t)256*1, rsyminfo+(ptrdiff_t)256*1, &ctxmask, 3*NCTX+1);
			enc_hist2stats(rhist+(ptrdiff_t)256*2, rsyminfo+(ptrdiff_t)256*2, &ctxmask, 3*NCTX+2);
			
			unsigned state=1<<(RANS_STATE_BITS-RANS_RENORM_BITS);
			for(int kx=xremw-1;kx>=0;--kx)
				encode1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, &state, &streamptr, image, rsyminfo);
			for(int ky=yremh-1;ky>=0;--ky)
				encode1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, &state, &streamptr, image, rsyminfo);
			//flush
			streamptr-=4;
#ifdef _DEBUG
			if(streamptr<=image)
				LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
			*(unsigned*)streamptr=state;
			prof_checkpoint(usize-isize, "encode remainder");
			profile_size(streamptr, "/ %9td bytes remainders", usize-isize);
		}

		//encode main
		mstate[1]=mstate[0]=_mm256_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
		unsigned short *ctxptr2=(unsigned short*)(interleaved+(isize<<1)-sizeof(short[NCODERS]));
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

					ctxptr2-=NCODERS;
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
							LOG_ERROR("freq = %d", freq);
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
					LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
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
		}
		//flush
		streamptr-=sizeof(mstate);
#ifdef _DEBUG
		if(streamptr<=image)
			LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
		memcpy(streamptr, mstate, sizeof(mstate));
		prof_checkpoint(isize, "encode main");
		profile_size(streamptr, "/ %9td bytes main", isize);

		//pack hists
		{
			BitPackerLIFO ec;
			bitpacker_enc_init(&ec, image, streamptr);
			if(xremw||yremh)
			{
				enc_packhist(&ec, rhist+256*2, ctxmask, 3*NCTX+2);
				enc_packhist(&ec, rhist+256*1, ctxmask, 3*NCTX+1);
				enc_packhist(&ec, rhist+256*0, ctxmask, 3*NCTX+0);
			}
			for(int kc=3*NCTX-1;kc>=0;--kc)
				enc_packhist(&ec, hists+(ptrdiff_t)256*kc, ctxmask, kc);
			bitpacker_enc_flush(&ec);
			streamptr=ec.dstbwdptr;

			streamptr-=8;
			*(unsigned long long*)streamptr=ctxmask;
		}
		prof_checkpoint(((ptrdiff_t)3*NCTX+((xremw||yremh)?3:0))*256, "pack histograms");
		profile_size(streamptr, "/ %9d bytes overhead", (3*NCTX+3)*12<<8>>3);

		//save compressed file
		{
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", fdst);
				return 1;
			}
			ptrdiff_t csize2=0;
			csize2+=fwrite("32", 1, 2, fdst);
			csize2+=fwrite(&iw, 1, 4, fdst);
			csize2+=fwrite(&ih, 1, 4, fdst);
			int flags=bestrct<<2|(use_wg4&3);
			csize2+=fwrite(&flags, 1, 1, fdst);
			csize2+=fwrite(&dist, 1, 1, fdst);
			//if(dist>1)
			//{
			//	csize2+=fwrite(lossymean, 1, sizeof(lossymean), fdst);
			//	csize2+=fwrite(lossytransform, 1, sizeof(lossytransform), fdst);
			//}
#ifdef _DEBUG
			if(streamptr>streamstart)
				LOG_ERROR("OOB ptr %016zX > %016zX", streamptr, streamstart);
			if(streamptr<image)
				LOG_ERROR("OOB ptr %016zX < %016zX", streamptr, image);
#endif
			csize2+=fwrite(streamptr, 1, streamstart-streamptr, fdst);
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
		free(rhist);
#ifdef WG4_PRINTMAXERR
		for(int k=0;k<16*8;++k)
		{
			printf("maxerror[%3d] %7d\n", k, maxerror[k]);
			if((k&15)==15)
				printf("\n");
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
				LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
			unsigned state=*(unsigned*)streamptr;
			streamptr+=4;
			for(int ky=0;ky<yremh;++ky)
				decode1d(image+rowstride*(blockh*YCODERS+ky), iw, 3, bestrct, &state, (const unsigned char**)&streamptr, streamend, rCDF2syms);
			for(int kx=0;kx<xremw;++kx)
				decode1d(image+qxbytes+3*kx, blockh*YCODERS, rowstride, bestrct, &state, (const unsigned char**)&streamptr, streamend, rCDF2syms);
			prof_checkpoint(usize-isize, "remainder");
		}

		//save PPM file
		save_ppm(dstfn, image, iw, ih);
		prof_checkpoint(usize, "fwrite");
	}
	_mm_free(pixels);
	_mm_free(CDF2syms);
	_mm_free(rCDF2syms);
	_mm_free(interleaved);
	free(image);

#ifdef LOUD
	t=time_sec()-t;
	printf("%c  %12.6lf sec  %12.6lf MB/s\n", 'D'+fwd, t, usize/(t*1024*1024));
#endif
#ifdef PROFILE_TIME
#ifdef __GNUC__
	if(nthreads0&4)
#endif
		prof_print(usize);
#endif
	(void)och_names;
	(void)rct_names;
	return 0;
}
