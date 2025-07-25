#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
//#ifdef _WIN32
//#include<Windows.h>
//#endif
static const char file[]=__FILE__;


	#define PROFILE_TIME		//should be on
	#define INTMIX

#ifdef _MSC_VER
	#define PROFILE_SIZE
	#define LOUD			//size & time

//	#define ESTIMATE_SIZE		//DEBUG		checks for zero frequency, visualizes context usage
	#define ENABLE_GUIDE		//DEBUG		checks interleaved pixels
//	#define ANS_VAL			//DEBUG

//	#define WG4_SERIALDEBUG
//	#define TEST_INTERLEAVE
#endif

//	#define ENABLE_MA		//bad		DIV2K +2.15 MB +0.7%
//	#define EMULATE_GATHER
//	#define CHECK_FLAT		//slower encode
//	#define DISABLE_ANALYSIS	//HORRIBLE	analysis is crucial
//	#define SERIAL_MAIN		//2x slower
//	#define DISABLE_WG

//	#define NAIVEMIX		//3.5 3.8x slower	55 58 MB/s vs 196 224 MB/s
//	#define WG_COMMONMIX		//bad
	#define WG_ENABLE_eW		//eW is bad with blocks unlike in eBench
	#define INTERLEAVESIMD		//2.5x faster interleave


//3*17+3=54 contexts
#define GRBITS 3
#define NCTX 17

#define XCODERS 4	//xrem 1~3 cols
#define YCODERS 8	//yrem 1~7 rows

#define NCODERS 32

#define PROBBITS 12	//12 bit max	James Bonfield's CDF2sym: {freq<<20 | bias<<8 | sym}

#define RANS_STATE_BITS 31
#define RANS_RENORM_BITS 16

#include"entropy.h"
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
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
	int printed=0;
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
#if 1
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
#endif
#if 0
		if(info->msg&&space>=len)
		{
			int labelstart=(space-len)>>1;
			int labelend=labelstart+len;
			if(labelstart>2&&space-labelend>2)
			{
				for(int k2=0;k2<labelstart;++k2)
					printf("%c", (printed+k2)&1||k2==labelstart-1?' ':k+'A');
			}
			else
				printf("%*s", labelstart, "");
			for(int k2=labelstart;k2<labelend;++k2)
				printf("%c", info->msg[k2-labelstart]);
			if(labelstart>2&&space-labelend>2)
			{
				for(int k2=labelend;k2<space;++k2)
					printf("%c", (printed+k2)&1||k2==labelend?' ':k+'A');
			}
			else
				printf("%*s", space-labelend, "");
		}
		else
		{
			for(int k2=0;k2<space;++k2)
				printf("%c", (printed+k2)&1?' ':k+'A');
		}
		printf(" |");
		if(k<prof_count-1)
			printf(" ");
#endif
		prev=curr;
		printed+=space+2+(k<prof_count-1);
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
			printf(" %12.6lf MB/s %10td bytes ", info->size/(info->t*1024*1024), info->size);
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


#define WG_NPREDS	8	//multiple of 4

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
	//wgpreds	short{Y[32], U[32], V[32]}{NPREDS} = __m256i[NPREDS][6]		6 = 3*NCODERS/16	16 = sizeof(__m256i)/sizeof(short)
	//	stride = sizeof(__m256i[6]) = 192 bytes

	//prederrors	uint16{32 lanes}{NPREDS} = __m256i[NPREDS][2]
	//	stride = sizeof(__m256i[2]) = 64 bytes

	//__m256[16] pr/coeffs		float{8 lanes}{NPREDS}{even/odd} = float[2 (even/odd)][8 PREDS][8 "strided" lanes]
	//	faststride = sizeof(__m256) = 32 bytes		slowstride = sizeof(__m256[8]) = 256 bytes


	//mix(const short wgpreds[8, stride 6][16], unsigned short prederrors[8, stride 2][16]) -> float[8][16]
	__m256i estim[16], coeff[16];

	//convert 16 shorts * 8 regs -> 8 floats * 16 regs		each 2 float regs make half a wave
	__m256i wordmask=_mm256_set1_epi32(0xFFFF);
	coeff[0*2+0]=_mm256_and_si256(prederrors[0*2+0], wordmask);
	coeff[1*2+0]=_mm256_and_si256(prederrors[1*2+0], wordmask);
	coeff[2*2+0]=_mm256_and_si256(prederrors[2*2+0], wordmask);
	coeff[3*2+0]=_mm256_and_si256(prederrors[3*2+0], wordmask);
	coeff[4*2+0]=_mm256_and_si256(prederrors[4*2+0], wordmask);
	coeff[5*2+0]=_mm256_and_si256(prederrors[5*2+0], wordmask);
	coeff[6*2+0]=_mm256_and_si256(prederrors[6*2+0], wordmask);
	coeff[7*2+0]=_mm256_and_si256(prederrors[7*2+0], wordmask);
	coeff[0*2+1]=_mm256_srli_epi32(prederrors[0*2+0], 16);
	coeff[1*2+1]=_mm256_srli_epi32(prederrors[1*2+0], 16);
	coeff[2*2+1]=_mm256_srli_epi32(prederrors[2*2+0], 16);
	coeff[3*2+1]=_mm256_srli_epi32(prederrors[3*2+0], 16);
	coeff[4*2+1]=_mm256_srli_epi32(prederrors[4*2+0], 16);
	coeff[5*2+1]=_mm256_srli_epi32(prederrors[5*2+0], 16);
	coeff[6*2+1]=_mm256_srli_epi32(prederrors[6*2+0], 16);
	coeff[7*2+1]=_mm256_srli_epi32(prederrors[7*2+0], 16);
	
	//convert 16 shorts * 8 regs -> 8 floats * 16 regs
	estim[0*2+0]=_mm256_slli_epi32(wgpreds[0*6+0], 16);
	estim[1*2+0]=_mm256_slli_epi32(wgpreds[1*6+0], 16);
	estim[2*2+0]=_mm256_slli_epi32(wgpreds[2*6+0], 16);
	estim[3*2+0]=_mm256_slli_epi32(wgpreds[3*6+0], 16);
	estim[4*2+0]=_mm256_slli_epi32(wgpreds[4*6+0], 16);
	estim[5*2+0]=_mm256_slli_epi32(wgpreds[5*6+0], 16);
	estim[6*2+0]=_mm256_slli_epi32(wgpreds[6*6+0], 16);
	estim[7*2+0]=_mm256_slli_epi32(wgpreds[7*6+0], 16);
	estim[0*2+0]=_mm256_srai_epi32(estim[0*2+0], 16);//preds are signed
	estim[1*2+0]=_mm256_srai_epi32(estim[1*2+0], 16);
	estim[2*2+0]=_mm256_srai_epi32(estim[2*2+0], 16);
	estim[3*2+0]=_mm256_srai_epi32(estim[3*2+0], 16);
	estim[4*2+0]=_mm256_srai_epi32(estim[4*2+0], 16);
	estim[5*2+0]=_mm256_srai_epi32(estim[5*2+0], 16);
	estim[6*2+0]=_mm256_srai_epi32(estim[6*2+0], 16);
	estim[7*2+0]=_mm256_srai_epi32(estim[7*2+0], 16);
	estim[0*2+1]=_mm256_srai_epi32(wgpreds[0*6+0], 16);
	estim[1*2+1]=_mm256_srai_epi32(wgpreds[1*6+0], 16);
	estim[2*2+1]=_mm256_srai_epi32(wgpreds[2*6+0], 16);
	estim[3*2+1]=_mm256_srai_epi32(wgpreds[3*6+0], 16);
	estim[4*2+1]=_mm256_srai_epi32(wgpreds[4*6+0], 16);
	estim[5*2+1]=_mm256_srai_epi32(wgpreds[5*6+0], 16);
	estim[6*2+1]=_mm256_srai_epi32(wgpreds[6*6+0], 16);
	estim[7*2+1]=_mm256_srai_epi32(wgpreds[7*6+0], 16);
#ifdef WG4_PRINTMAXERR
	for(int k=0;k<128;++k)
	{
		float e=((float*)coeff)[k];
		if(maxerror<e)
			maxerror=e;
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
#if 1
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
			140,
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
	//coeff[0x0]=_mm256_i32gather_epi32(divlookup[0], coeff[0x0], sizeof(int));
	//coeff[0x1]=_mm256_i32gather_epi32(divlookup[0], coeff[0x1], sizeof(int));
	//coeff[0x2]=_mm256_i32gather_epi32(divlookup[1], coeff[0x2], sizeof(int));
	//coeff[0x3]=_mm256_i32gather_epi32(divlookup[1], coeff[0x3], sizeof(int));
	//coeff[0x4]=_mm256_i32gather_epi32(divlookup[2], coeff[0x4], sizeof(int));
	//coeff[0x5]=_mm256_i32gather_epi32(divlookup[2], coeff[0x5], sizeof(int));
	//coeff[0x6]=_mm256_i32gather_epi32(divlookup[3], coeff[0x6], sizeof(int));
	//coeff[0x7]=_mm256_i32gather_epi32(divlookup[3], coeff[0x7], sizeof(int));
	//coeff[0x8]=_mm256_i32gather_epi32(divlookup[4], coeff[0x8], sizeof(int));
	//coeff[0x9]=_mm256_i32gather_epi32(divlookup[4], coeff[0x9], sizeof(int));
	//coeff[0xA]=_mm256_i32gather_epi32(divlookup[5], coeff[0xA], sizeof(int));
	//coeff[0xB]=_mm256_i32gather_epi32(divlookup[5], coeff[0xB], sizeof(int));
	//coeff[0xC]=_mm256_i32gather_epi32(divlookup[6], coeff[0xC], sizeof(int));
	//coeff[0xD]=_mm256_i32gather_epi32(divlookup[6], coeff[0xD], sizeof(int));
	//coeff[0xE]=_mm256_i32gather_epi32(divlookup[7], coeff[0xE], sizeof(int));
	//coeff[0xF]=_mm256_i32gather_epi32(divlookup[7], coeff[0xF], sizeof(int));
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
	__m256i ipred[2];
	ipred[0]=_mm256_srli_epi32(sum0, 1);
	ipred[1]=_mm256_srli_epi32(sum1, 1);
	coeff[0*2+0]=_mm256_mullo_epi32(coeff[0*2+0], estim[0*2+0]);
	coeff[1*2+0]=_mm256_mullo_epi32(coeff[1*2+0], estim[1*2+0]);
	coeff[2*2+0]=_mm256_mullo_epi32(coeff[2*2+0], estim[2*2+0]);
	coeff[3*2+0]=_mm256_mullo_epi32(coeff[3*2+0], estim[3*2+0]);
	coeff[4*2+0]=_mm256_mullo_epi32(coeff[4*2+0], estim[4*2+0]);
	coeff[5*2+0]=_mm256_mullo_epi32(coeff[5*2+0], estim[5*2+0]);
	coeff[6*2+0]=_mm256_mullo_epi32(coeff[6*2+0], estim[6*2+0]);
	coeff[7*2+0]=_mm256_mullo_epi32(coeff[7*2+0], estim[7*2+0]);
	coeff[0*2+1]=_mm256_mullo_epi32(coeff[0*2+1], estim[0*2+1]);
	coeff[1*2+1]=_mm256_mullo_epi32(coeff[1*2+1], estim[1*2+1]);
	coeff[2*2+1]=_mm256_mullo_epi32(coeff[2*2+1], estim[2*2+1]);
	coeff[3*2+1]=_mm256_mullo_epi32(coeff[3*2+1], estim[3*2+1]);
	coeff[4*2+1]=_mm256_mullo_epi32(coeff[4*2+1], estim[4*2+1]);
	coeff[5*2+1]=_mm256_mullo_epi32(coeff[5*2+1], estim[5*2+1]);
	coeff[6*2+1]=_mm256_mullo_epi32(coeff[6*2+1], estim[6*2+1]);
	coeff[7*2+1]=_mm256_mullo_epi32(coeff[7*2+1], estim[7*2+1]);
	ipred[0]=_mm256_add_epi32(ipred[0], coeff[0*2+0]);
	ipred[0]=_mm256_add_epi32(ipred[0], coeff[1*2+0]);
	ipred[0]=_mm256_add_epi32(ipred[0], coeff[2*2+0]);
	ipred[0]=_mm256_add_epi32(ipred[0], coeff[3*2+0]);
	ipred[0]=_mm256_add_epi32(ipred[0], coeff[4*2+0]);
	ipred[0]=_mm256_add_epi32(ipred[0], coeff[5*2+0]);
	ipred[0]=_mm256_add_epi32(ipred[0], coeff[6*2+0]);
	ipred[0]=_mm256_add_epi32(ipred[0], coeff[7*2+0]);
	ipred[1]=_mm256_add_epi32(ipred[1], coeff[0*2+1]);
	ipred[1]=_mm256_add_epi32(ipred[1], coeff[1*2+1]);
	ipred[1]=_mm256_add_epi32(ipred[1], coeff[2*2+1]);
	ipred[1]=_mm256_add_epi32(ipred[1], coeff[3*2+1]);
	ipred[1]=_mm256_add_epi32(ipred[1], coeff[4*2+1]);
	ipred[1]=_mm256_add_epi32(ipred[1], coeff[5*2+1]);
	ipred[1]=_mm256_add_epi32(ipred[1], coeff[6*2+1]);
	ipred[1]=_mm256_add_epi32(ipred[1], coeff[7*2+1]);
	//ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[0*2+0], estim[0*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[0*2+1], estim[0*2+1]));
	//ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[1*2+0], estim[1*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[1*2+1], estim[1*2+1]));
	//ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[2*2+0], estim[2*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[2*2+1], estim[2*2+1]));
	//ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[3*2+0], estim[3*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[3*2+1], estim[3*2+1]));
	//ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[4*2+0], estim[4*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[4*2+1], estim[4*2+1]));
	//ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[5*2+0], estim[5*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[5*2+1], estim[5*2+1]));
	//ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[6*2+0], estim[6*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[6*2+1], estim[6*2+1]));
	//ipred[0]=_mm256_add_epi32(ipred[0], _mm256_mullo_epi32(coeff[7*2+0], estim[7*2+0])); ipred[1]=_mm256_add_epi32(ipred[1], _mm256_mullo_epi32(coeff[7*2+1], estim[7*2+1]));
	gather32((int*)&sum0, divlookup[WG_NPREDS], (int*)&sum0);
	gather32((int*)&sum1, divlookup[WG_NPREDS], (int*)&sum1);
	//sum0=_mm256_i32gather_epi32(divlookup[WG_NPREDS], sum0, sizeof(int));
	//sum1=_mm256_i32gather_epi32(divlookup[WG_NPREDS], sum1, sizeof(int));
	ipred[0]=_mm256_mullo_epi32(ipred[0], sum0);
	ipred[1]=_mm256_mullo_epi32(ipred[1], sum1);
	result[0]=_mm256_srai_epi32(ipred[0], NUMBITS);
	result[1]=_mm256_srai_epi32(ipred[1], NUMBITS);
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
#if 0
	for(int kl=0;kl<NCODERS/2;++kl)
	{
		static const int weights[]=
		{
			240<<10,
			240<<10,
			180<<10,
			180<<10,
			140<<10,
			160<<10,
			120<<10,
			120<<10,
		};
#if 1
		long long ipred=0;
		unsigned isum=0;
		for(int kp=0;kp<WG_NPREDS;++kp)
		{
			int weight=((int*)coeff)[NCODERS/2*kp+kl];
			weight=weights[kp]/weight;
			int pred=((int*)estim)[NCODERS/2*kp+kl];
			ipred+=(long long)weight*pred;
			isum+=weight;
		}
		((int*)result)[kl]=(int)(ipred/isum);
#endif
#if 0
		//int it=0;
		long long ipred=0;
		unsigned isum=0;
		int powers[WG_NPREDS], maxpower;
	//again:
		maxpower=0;
		for(int kp=0;kp<WG_NPREDS;++kp)
		{
			int ie=((int*)coeff)[NCODERS/2*kp+kl];
			powers[kp]=FLOOR_LOG2_P1(ie);
			if(maxpower<powers[kp])
				maxpower=powers[kp];
		}
		++maxpower;
		for(int kp=0;kp<WG_NPREDS;++kp)
		{
			int pred=((int*)estim)[NCODERS/2*kp+kl];
			int coeff=weights[kp]*(maxpower-powers[kp]);
			ipred+=(long long)coeff*pred;
			isum+=coeff;
		}
		//ipred/=isum;
		((int*)result)[kl]=(int)ipred/isum;
#endif
	}
#endif
#if 0
#ifdef WG4_DEBUG
	ALIGN(32) float apr[sizeof(__m256[16])/(sizeof(float))], acf[sizeof(__m256[16])/(sizeof(float))];
	memcpy(acf, coeff, sizeof(acf));
	memcpy(apr, pr, sizeof(apr));
#endif
	//coeff = 1/error
	coeff[0x0]=_mm256_rcp_ps(coeff[0x0]);//4/0.5	~11 bit prec, max error < (1.5 / 0x1000)
	coeff[0x1]=_mm256_rcp_ps(coeff[0x1]);
	coeff[0x2]=_mm256_rcp_ps(coeff[0x2]);
	coeff[0x3]=_mm256_rcp_ps(coeff[0x3]);
	coeff[0x4]=_mm256_rcp_ps(coeff[0x4]);
	coeff[0x5]=_mm256_rcp_ps(coeff[0x5]);
	coeff[0x6]=_mm256_rcp_ps(coeff[0x6]);
	coeff[0x7]=_mm256_rcp_ps(coeff[0x7]);
	coeff[0x8]=_mm256_rcp_ps(coeff[0x8]);
	coeff[0x9]=_mm256_rcp_ps(coeff[0x9]);
	coeff[0xA]=_mm256_rcp_ps(coeff[0xA]);
	coeff[0xB]=_mm256_rcp_ps(coeff[0xB]);
	coeff[0xC]=_mm256_rcp_ps(coeff[0xC]);
	coeff[0xD]=_mm256_rcp_ps(coeff[0xD]);
	coeff[0xE]=_mm256_rcp_ps(coeff[0xE]);
	coeff[0xF]=_mm256_rcp_ps(coeff[0xF]);
	__m256 rcpmask=_mm256_castsi256_ps(_mm256_set1_epi32(0xFFFFE000));//use 11 bits
	coeff[0x0]=_mm256_and_ps(coeff[0x0], rcpmask);
	coeff[0x1]=_mm256_and_ps(coeff[0x1], rcpmask);
	coeff[0x2]=_mm256_and_ps(coeff[0x2], rcpmask);
	coeff[0x3]=_mm256_and_ps(coeff[0x3], rcpmask);
	coeff[0x4]=_mm256_and_ps(coeff[0x4], rcpmask);
	coeff[0x5]=_mm256_and_ps(coeff[0x5], rcpmask);
	coeff[0x6]=_mm256_and_ps(coeff[0x6], rcpmask);
	coeff[0x7]=_mm256_and_ps(coeff[0x7], rcpmask);
	coeff[0x8]=_mm256_and_ps(coeff[0x8], rcpmask);
	coeff[0x9]=_mm256_and_ps(coeff[0x9], rcpmask);
	coeff[0xA]=_mm256_and_ps(coeff[0xA], rcpmask);
	coeff[0xB]=_mm256_and_ps(coeff[0xB], rcpmask);
	coeff[0xC]=_mm256_and_ps(coeff[0xC], rcpmask);
	coeff[0xD]=_mm256_and_ps(coeff[0xD], rcpmask);
	coeff[0xE]=_mm256_and_ps(coeff[0xE], rcpmask);
	coeff[0xF]=_mm256_and_ps(coeff[0xF], rcpmask);
//#ifdef EXPERIMENT_RCPPS
//	unsigned *ptr=(unsigned*)coeff;
//	for(int k=0;k<sizeof(coeff)/sizeof(int);++k)
//	{
//		if(ptr[k]&4095)
//			LOG_ERROR("%08X", ptr[k]);
//	}
//#endif
	//__m256 num=_mm256_set1_ps(1);
	//coeff[0x0]=_mm256_div_ps(num, coeff[0x0]);//11/5
	//coeff[0x1]=_mm256_div_ps(num, coeff[0x1]);
	//coeff[0x2]=_mm256_div_ps(num, coeff[0x2]);
	//coeff[0x3]=_mm256_div_ps(num, coeff[0x3]);
	//coeff[0x4]=_mm256_div_ps(num, coeff[0x4]);
	//coeff[0x5]=_mm256_div_ps(num, coeff[0x5]);
	//coeff[0x6]=_mm256_div_ps(num, coeff[0x6]);
	//coeff[0x7]=_mm256_div_ps(num, coeff[0x7]);
	//coeff[0x8]=_mm256_div_ps(num, coeff[0x8]);
	//coeff[0x9]=_mm256_div_ps(num, coeff[0x9]);
	//coeff[0xA]=_mm256_div_ps(num, coeff[0xA]);
	//coeff[0xB]=_mm256_div_ps(num, coeff[0xB]);
	//coeff[0xC]=_mm256_div_ps(num, coeff[0xC]);
	//coeff[0xD]=_mm256_div_ps(num, coeff[0xD]);
	//coeff[0xE]=_mm256_div_ps(num, coeff[0xE]);
	//coeff[0xF]=_mm256_div_ps(num, coeff[0xF]);

//	__m256 g0=_mm256_set1_ps(210);
// 	__m256 g1=g0;
//	__m256 g2=_mm256_set1_ps(215);
// 	__m256 g3=g2;
//	__m256 g4=_mm256_set1_ps(140);
//	__m256 g5=_mm256_set1_ps(230);
//	__m256 g6=_mm256_set1_ps(120);
//	__m256 g7=_mm256_set1_ps(140);
	coeff[0x0]=_mm256_mul_ps(coeff[0x0], g0);
	coeff[0x1]=_mm256_mul_ps(coeff[0x1], g0);
	coeff[0x2]=_mm256_mul_ps(coeff[0x2], g1);
	coeff[0x3]=_mm256_mul_ps(coeff[0x3], g1);
	coeff[0x4]=_mm256_mul_ps(coeff[0x4], g2);
	coeff[0x5]=_mm256_mul_ps(coeff[0x5], g2);
	coeff[0x6]=_mm256_mul_ps(coeff[0x6], g3);
	coeff[0x7]=_mm256_mul_ps(coeff[0x7], g3);
	coeff[0x8]=_mm256_mul_ps(coeff[0x8], g4);
	coeff[0x9]=_mm256_mul_ps(coeff[0x9], g4);
	coeff[0xA]=_mm256_mul_ps(coeff[0xA], g5);
	coeff[0xB]=_mm256_mul_ps(coeff[0xB], g5);
	coeff[0xC]=_mm256_mul_ps(coeff[0xC], g6);
	coeff[0xD]=_mm256_mul_ps(coeff[0xD], g6);
	coeff[0xE]=_mm256_mul_ps(coeff[0xE], g7);
	coeff[0xF]=_mm256_mul_ps(coeff[0xF], g7);

	//pr = (sum: coeff[i]*wg[i])/(sum: coeff[i])	NPREDS=8  (*16 interleaved floats)
	pr[0x0]=_mm256_mul_ps(pr[0x0], coeff[0x0]);
	pr[0x1]=_mm256_mul_ps(pr[0x1], coeff[0x1]);
	pr[0x2]=_mm256_mul_ps(pr[0x2], coeff[0x2]);
	pr[0x3]=_mm256_mul_ps(pr[0x3], coeff[0x3]);
	pr[0x4]=_mm256_mul_ps(pr[0x4], coeff[0x4]);
	pr[0x5]=_mm256_mul_ps(pr[0x5], coeff[0x5]);
	pr[0x6]=_mm256_mul_ps(pr[0x6], coeff[0x6]);
	pr[0x7]=_mm256_mul_ps(pr[0x7], coeff[0x7]);
	pr[0x8]=_mm256_mul_ps(pr[0x8], coeff[0x8]);
	pr[0x9]=_mm256_mul_ps(pr[0x9], coeff[0x9]);
	pr[0xA]=_mm256_mul_ps(pr[0xA], coeff[0xA]);
	pr[0xB]=_mm256_mul_ps(pr[0xB], coeff[0xB]);
	pr[0xC]=_mm256_mul_ps(pr[0xC], coeff[0xC]);
	pr[0xD]=_mm256_mul_ps(pr[0xD], coeff[0xD]);
	pr[0xE]=_mm256_mul_ps(pr[0xE], coeff[0xE]);
	pr[0xF]=_mm256_mul_ps(pr[0xF], coeff[0xF]);

	coeff[0x0]=_mm256_add_ps(coeff[0x0], coeff[0x2]);
	coeff[0x1]=_mm256_add_ps(coeff[0x1], coeff[0x3]);
	coeff[0x4]=_mm256_add_ps(coeff[0x4], coeff[0x6]);
	coeff[0x5]=_mm256_add_ps(coeff[0x5], coeff[0x7]);
	coeff[0x8]=_mm256_add_ps(coeff[0x8], coeff[0xA]);
	coeff[0x9]=_mm256_add_ps(coeff[0x9], coeff[0xB]);
	coeff[0xC]=_mm256_add_ps(coeff[0xC], coeff[0xE]);
	coeff[0xD]=_mm256_add_ps(coeff[0xD], coeff[0xF]);

	coeff[0x0]=_mm256_add_ps(coeff[0x0], coeff[0x4]);
	coeff[0x1]=_mm256_add_ps(coeff[0x1], coeff[0x5]);
	coeff[0x8]=_mm256_add_ps(coeff[0x8], coeff[0xC]);
	coeff[0x9]=_mm256_add_ps(coeff[0x9], coeff[0xD]);
	
	coeff[0x0]=_mm256_add_ps(coeff[0x0], coeff[0x8]);
	coeff[0x1]=_mm256_add_ps(coeff[0x1], coeff[0x9]);
	
	pr[0x0]=_mm256_add_ps(pr[0x0], pr[0x2]);
	pr[0x1]=_mm256_add_ps(pr[0x1], pr[0x3]);
	pr[0x4]=_mm256_add_ps(pr[0x4], pr[0x6]);
	pr[0x5]=_mm256_add_ps(pr[0x5], pr[0x7]);
	pr[0x8]=_mm256_add_ps(pr[0x8], pr[0xA]);
	pr[0x9]=_mm256_add_ps(pr[0x9], pr[0xB]);
	pr[0xC]=_mm256_add_ps(pr[0xC], pr[0xE]);
	pr[0xD]=_mm256_add_ps(pr[0xD], pr[0xF]);

	pr[0x0]=_mm256_add_ps(pr[0x0], pr[0x4]);
	pr[0x1]=_mm256_add_ps(pr[0x1], pr[0x5]);
	pr[0x8]=_mm256_add_ps(pr[0x8], pr[0xC]);
	pr[0x9]=_mm256_add_ps(pr[0x9], pr[0xD]);

	pr[0x0]=_mm256_add_ps(pr[0x0], pr[0x8]);
	pr[0x1]=_mm256_add_ps(pr[0x1], pr[0x9]);

	result[0]=_mm256_div_ps(pr[0], coeff[0]);//11/5
	result[1]=_mm256_div_ps(pr[1], coeff[1]);
#endif
#ifdef WG4_DEBUG
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
		//int it=0;
		//float simdpred=((float*)result)[kl];
		long long ipred=0;
		unsigned isum=0;
		int powers[WG_NPREDS], maxpower;
		//float f1pred, w1sum, pred1;
	//again:
		maxpower=0;
		for(int kp=0;kp<WG_NPREDS;++kp)
		{
			int ie=(int)acf[NCODERS/2*kp+kl];
			powers[kp]=FLOOR_LOG2_P1(ie);
			if(maxpower<powers[kp])
				maxpower=powers[kp];
		}
		++maxpower;
		for(int kp=0;kp<WG_NPREDS;++kp)
		{
			int pred=(int)apr[NCODERS/2*kp+kl];
			int coeff=weights[kp]*(maxpower-powers[kp]);
			ipred+=(long long)coeff*pred;
			isum+=coeff;
		}
		//ipred/=isum;
		((float*)result)[kl]=(float)ipred/isum;
#if 0
		f1pred=0, w1sum=0;
		for(int kp=0;kp<WG_NPREDS;++kp)
		{
			float fe=acf[NCODERS/2*kp+kl], fsp=apr[NCODERS/2*kp+kl];
			float w1=1.f/fe;
			*(unsigned*)&w1&=0xFFFFE000;
			w1*=weights[kp];
			w1sum+=w1;
			f1pred+=w1*fsp;

			int ie=(int)fe, ip=(int)fsp;
			int lg=FLOOR_LOG2(ie);
		}
		pred1=f1pred/w1sum;
		if(fabsf(pred1-simdpred)>5)
		{
			if(it++<2)
				goto again;//
			LOG_ERROR("");
		}
		int a=CVTFP32_I32(simdpred);
		int b=(int)round(pred1);
		if(a!=b)
		{
			printf("\n");
			printf("SIMD %16.10lf %8d\n", simdpred, a);
			printf("EMUL %16.10lf %8d\n", pred1, b);
			printf("Error %16.10lf\n", simdpred-pred1);
			if(it++<2)
				goto again;//
			LOG_ERROR("");
		}
#endif
	}
#endif
}
AWM_INLINE void wg_mix_pt2f(const __m256i *wgpreds, const __m256i *prederrors, __m256 *result)//process halfwave
{
	//wgpreds	short{Y[32], U[32], V[32]}{NPREDS} = __m256i[NPREDS][6]		6 = 3*NCODERS/16	16 = sizeof(__m256i)/sizeof(short)
	//	stride = sizeof(__m256i[6]) = 192 bytes

	//prederrors	uint16{32 lanes}{NPREDS} = __m256i[NPREDS][2]
	//	stride = sizeof(__m256i[2]) = 64 bytes

	//__m256[16] pr/coeffs		float{8 lanes}{NPREDS}{even/odd} = float[2 (even/odd)][8 PREDS][8 "strided" lanes]
	//	faststride = sizeof(__m256) = 32 bytes		slowstride = sizeof(__m256[8]) = 256 bytes


	//mix(const short wgpreds[8, stride 6][16], unsigned short prederrors[8, stride 2][16]) -> float[8][16]
	__m256 pr[16], coeff[16];

	//convert 16 shorts * 8 regs -> 8 floats * 16 regs		each 2 float regs make half a wave
	{
		__m256i t0[8];
		t0[0]=_mm256_slli_epi32(prederrors[0*2+0], 16);
		t0[1]=_mm256_slli_epi32(prederrors[1*2+0], 16);
		t0[2]=_mm256_slli_epi32(prederrors[2*2+0], 16);
		t0[3]=_mm256_slli_epi32(prederrors[3*2+0], 16);
		t0[4]=_mm256_slli_epi32(prederrors[4*2+0], 16);
		t0[5]=_mm256_slli_epi32(prederrors[5*2+0], 16);
		t0[6]=_mm256_slli_epi32(prederrors[6*2+0], 16);
		t0[7]=_mm256_slli_epi32(prederrors[7*2+0], 16);
		t0[0]=_mm256_srli_epi32(t0[0], 16);//errors are unsigned
		t0[1]=_mm256_srli_epi32(t0[1], 16);
		t0[2]=_mm256_srli_epi32(t0[2], 16);
		t0[3]=_mm256_srli_epi32(t0[3], 16);
		t0[4]=_mm256_srli_epi32(t0[4], 16);
		t0[5]=_mm256_srli_epi32(t0[5], 16);
		t0[6]=_mm256_srli_epi32(t0[6], 16);
		t0[7]=_mm256_srli_epi32(t0[7], 16);
		coeff[0*2+0]=_mm256_cvtepi32_ps(t0[0]);
		coeff[1*2+0]=_mm256_cvtepi32_ps(t0[1]);
		coeff[2*2+0]=_mm256_cvtepi32_ps(t0[2]);
		coeff[3*2+0]=_mm256_cvtepi32_ps(t0[3]);
		coeff[4*2+0]=_mm256_cvtepi32_ps(t0[4]);
		coeff[5*2+0]=_mm256_cvtepi32_ps(t0[5]);
		coeff[6*2+0]=_mm256_cvtepi32_ps(t0[6]);
		coeff[7*2+0]=_mm256_cvtepi32_ps(t0[7]);
		t0[0]=_mm256_srli_epi32(prederrors[0*2+0], 16);
		t0[1]=_mm256_srli_epi32(prederrors[1*2+0], 16);
		t0[2]=_mm256_srli_epi32(prederrors[2*2+0], 16);
		t0[3]=_mm256_srli_epi32(prederrors[3*2+0], 16);
		t0[4]=_mm256_srli_epi32(prederrors[4*2+0], 16);
		t0[5]=_mm256_srli_epi32(prederrors[5*2+0], 16);
		t0[6]=_mm256_srli_epi32(prederrors[6*2+0], 16);
		t0[7]=_mm256_srli_epi32(prederrors[7*2+0], 16);
		coeff[0*2+1]=_mm256_cvtepi32_ps(t0[0]);
		coeff[1*2+1]=_mm256_cvtepi32_ps(t0[1]);
		coeff[2*2+1]=_mm256_cvtepi32_ps(t0[2]);
		coeff[3*2+1]=_mm256_cvtepi32_ps(t0[3]);
		coeff[4*2+1]=_mm256_cvtepi32_ps(t0[4]);
		coeff[5*2+1]=_mm256_cvtepi32_ps(t0[5]);
		coeff[6*2+1]=_mm256_cvtepi32_ps(t0[6]);
		coeff[7*2+1]=_mm256_cvtepi32_ps(t0[7]);
#if 0
		__m128i hi[8];
		__m256i i32[8];
		hi[0]=_mm256_extracti128_si256(prederrors[0], 1);
		hi[1]=_mm256_extracti128_si256(prederrors[1], 1);
		hi[2]=_mm256_extracti128_si256(prederrors[2], 1);
		hi[3]=_mm256_extracti128_si256(prederrors[3], 1);
		hi[4]=_mm256_extracti128_si256(prederrors[4], 1);
		hi[5]=_mm256_extracti128_si256(prederrors[5], 1);
		hi[6]=_mm256_extracti128_si256(prederrors[6], 1);
		hi[7]=_mm256_extracti128_si256(prederrors[7], 1);
		i32[0]=_mm256_cvtepi16_epi32(hi[0]);
		i32[1]=_mm256_cvtepi16_epi32(hi[1]);
		i32[2]=_mm256_cvtepi16_epi32(hi[2]);
		i32[3]=_mm256_cvtepi16_epi32(hi[3]);
		i32[4]=_mm256_cvtepi16_epi32(hi[4]);
		i32[5]=_mm256_cvtepi16_epi32(hi[5]);
		i32[6]=_mm256_cvtepi16_epi32(hi[6]);
		i32[7]=_mm256_cvtepi16_epi32(hi[7]);
		coeff[0x1]=_mm256_cvtepi32_ps(i32[0]);
		coeff[0x3]=_mm256_cvtepi32_ps(i32[1]);
		coeff[0x5]=_mm256_cvtepi32_ps(i32[2]);
		coeff[0x7]=_mm256_cvtepi32_ps(i32[3]);
		coeff[0x9]=_mm256_cvtepi32_ps(i32[4]);
		coeff[0xB]=_mm256_cvtepi32_ps(i32[5]);
		coeff[0xD]=_mm256_cvtepi32_ps(i32[6]);
		coeff[0xF]=_mm256_cvtepi32_ps(i32[7]);
		i32[0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[0]));
		i32[1]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[1]));
		i32[2]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[2]));
		i32[3]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[3]));
		i32[4]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[4]));
		i32[5]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[5]));
		i32[6]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[6]));
		i32[7]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(prederrors[7]));
		coeff[0x0]=_mm256_cvtepi32_ps(i32[0]);
		coeff[0x2]=_mm256_cvtepi32_ps(i32[1]);
		coeff[0x4]=_mm256_cvtepi32_ps(i32[2]);
		coeff[0x6]=_mm256_cvtepi32_ps(i32[3]);
		coeff[0x8]=_mm256_cvtepi32_ps(i32[4]);
		coeff[0xA]=_mm256_cvtepi32_ps(i32[5]);
		coeff[0xC]=_mm256_cvtepi32_ps(i32[6]);
		coeff[0xE]=_mm256_cvtepi32_ps(i32[7]);
#endif
	}
	
	//convert 16 shorts * 8 regs -> 8 floats * 16 regs
	{
		__m256i t0[8];
		t0[0]=_mm256_slli_epi32(wgpreds[0*6+0], 16);
		t0[1]=_mm256_slli_epi32(wgpreds[1*6+0], 16);
		t0[2]=_mm256_slli_epi32(wgpreds[2*6+0], 16);
		t0[3]=_mm256_slli_epi32(wgpreds[3*6+0], 16);
		t0[4]=_mm256_slli_epi32(wgpreds[4*6+0], 16);
		t0[5]=_mm256_slli_epi32(wgpreds[5*6+0], 16);
		t0[6]=_mm256_slli_epi32(wgpreds[6*6+0], 16);
		t0[7]=_mm256_slli_epi32(wgpreds[7*6+0], 16);
		t0[0]=_mm256_srai_epi32(t0[0], 16);//preds are signed
		t0[1]=_mm256_srai_epi32(t0[1], 16);
		t0[2]=_mm256_srai_epi32(t0[2], 16);
		t0[3]=_mm256_srai_epi32(t0[3], 16);
		t0[4]=_mm256_srai_epi32(t0[4], 16);
		t0[5]=_mm256_srai_epi32(t0[5], 16);
		t0[6]=_mm256_srai_epi32(t0[6], 16);
		t0[7]=_mm256_srai_epi32(t0[7], 16);
		pr[0*2+0]=_mm256_cvtepi32_ps(t0[0]);
		pr[1*2+0]=_mm256_cvtepi32_ps(t0[1]);
		pr[2*2+0]=_mm256_cvtepi32_ps(t0[2]);
		pr[3*2+0]=_mm256_cvtepi32_ps(t0[3]);
		pr[4*2+0]=_mm256_cvtepi32_ps(t0[4]);
		pr[5*2+0]=_mm256_cvtepi32_ps(t0[5]);
		pr[6*2+0]=_mm256_cvtepi32_ps(t0[6]);
		pr[7*2+0]=_mm256_cvtepi32_ps(t0[7]);
		t0[0]=_mm256_srai_epi32(wgpreds[0*6+0], 16);
		t0[1]=_mm256_srai_epi32(wgpreds[1*6+0], 16);
		t0[2]=_mm256_srai_epi32(wgpreds[2*6+0], 16);
		t0[3]=_mm256_srai_epi32(wgpreds[3*6+0], 16);
		t0[4]=_mm256_srai_epi32(wgpreds[4*6+0], 16);
		t0[5]=_mm256_srai_epi32(wgpreds[5*6+0], 16);
		t0[6]=_mm256_srai_epi32(wgpreds[6*6+0], 16);
		t0[7]=_mm256_srai_epi32(wgpreds[7*6+0], 16);
		pr[0*2+1]=_mm256_cvtepi32_ps(t0[0]);
		pr[1*2+1]=_mm256_cvtepi32_ps(t0[1]);
		pr[2*2+1]=_mm256_cvtepi32_ps(t0[2]);
		pr[3*2+1]=_mm256_cvtepi32_ps(t0[3]);
		pr[4*2+1]=_mm256_cvtepi32_ps(t0[4]);
		pr[5*2+1]=_mm256_cvtepi32_ps(t0[5]);
		pr[6*2+1]=_mm256_cvtepi32_ps(t0[6]);
		pr[7*2+1]=_mm256_cvtepi32_ps(t0[7]);
#if 0
		__m256i t0[16];
		t0[0*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[0*6+0], 1));
		t0[1*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[1*6+0], 1));
		t0[2*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[2*6+0], 1));
		t0[3*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[3*6+0], 1));
		t0[4*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[4*6+0], 1));
		t0[5*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[5*6+0], 1));
		t0[6*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[6*6+0], 1));
		t0[7*2+1]=_mm256_cvtepi16_epi32(_mm256_extracti128_si256(wgpreds[7*6+0], 1));
		t0[0*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[0*6+0]));
		t0[1*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[1*6+0]));
		t0[2*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[2*6+0]));
		t0[3*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[3*6+0]));
		t0[4*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[4*6+0]));
		t0[5*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[5*6+0]));
		t0[6*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[6*6+0]));
		t0[7*2+0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(wgpreds[7*6+0]));

		pr[0x0]=_mm256_cvtepi32_ps(t0[0x0]);
		pr[0x1]=_mm256_cvtepi32_ps(t0[0x1]);
		pr[0x2]=_mm256_cvtepi32_ps(t0[0x2]);
		pr[0x3]=_mm256_cvtepi32_ps(t0[0x3]);
		pr[0x4]=_mm256_cvtepi32_ps(t0[0x4]);
		pr[0x5]=_mm256_cvtepi32_ps(t0[0x5]);
		pr[0x6]=_mm256_cvtepi32_ps(t0[0x6]);
		pr[0x7]=_mm256_cvtepi32_ps(t0[0x7]);
		pr[0x8]=_mm256_cvtepi32_ps(t0[0x8]);
		pr[0x9]=_mm256_cvtepi32_ps(t0[0x9]);
		pr[0xA]=_mm256_cvtepi32_ps(t0[0xA]);
		pr[0xB]=_mm256_cvtepi32_ps(t0[0xB]);
		pr[0xC]=_mm256_cvtepi32_ps(t0[0xC]);
		pr[0xD]=_mm256_cvtepi32_ps(t0[0xD]);
		pr[0xE]=_mm256_cvtepi32_ps(t0[0xE]);
		pr[0xF]=_mm256_cvtepi32_ps(t0[0xF]);
#endif
	}
#ifdef WG4_SERIALDEBUG
	ALIGN(32) float cvtpe[sizeof(__m256[16])/sizeof(float)], cvtsp[sizeof(__m256[16])/sizeof(float)];
	memcpy(cvtpe, coeff, sizeof(cvtpe));
	memcpy(cvtsp, pr, sizeof(cvtsp));
#endif

	//coeff = 1/error
	coeff[0x0]=_mm256_rcp_ps(coeff[0x0]);//4/0.5	~11 bit prec, max error < (1.5 / 0x1000)
	coeff[0x1]=_mm256_rcp_ps(coeff[0x1]);
	coeff[0x2]=_mm256_rcp_ps(coeff[0x2]);
	coeff[0x3]=_mm256_rcp_ps(coeff[0x3]);
	coeff[0x4]=_mm256_rcp_ps(coeff[0x4]);
	coeff[0x5]=_mm256_rcp_ps(coeff[0x5]);
	coeff[0x6]=_mm256_rcp_ps(coeff[0x6]);
	coeff[0x7]=_mm256_rcp_ps(coeff[0x7]);
	coeff[0x8]=_mm256_rcp_ps(coeff[0x8]);
	coeff[0x9]=_mm256_rcp_ps(coeff[0x9]);
	coeff[0xA]=_mm256_rcp_ps(coeff[0xA]);
	coeff[0xB]=_mm256_rcp_ps(coeff[0xB]);
	coeff[0xC]=_mm256_rcp_ps(coeff[0xC]);
	coeff[0xD]=_mm256_rcp_ps(coeff[0xD]);
	coeff[0xE]=_mm256_rcp_ps(coeff[0xE]);
	coeff[0xF]=_mm256_rcp_ps(coeff[0xF]);
	//__m256 num=_mm256_set1_ps(1);
	//coeff[0x0]=_mm256_div_ps(num, coeff[0x0]);//11/5
	//coeff[0x1]=_mm256_div_ps(num, coeff[0x1]);
	//coeff[0x2]=_mm256_div_ps(num, coeff[0x2]);
	//coeff[0x3]=_mm256_div_ps(num, coeff[0x3]);
	//coeff[0x4]=_mm256_div_ps(num, coeff[0x4]);
	//coeff[0x5]=_mm256_div_ps(num, coeff[0x5]);
	//coeff[0x6]=_mm256_div_ps(num, coeff[0x6]);
	//coeff[0x7]=_mm256_div_ps(num, coeff[0x7]);
	//coeff[0x8]=_mm256_div_ps(num, coeff[0x8]);
	//coeff[0x9]=_mm256_div_ps(num, coeff[0x9]);
	//coeff[0xA]=_mm256_div_ps(num, coeff[0xA]);
	//coeff[0xB]=_mm256_div_ps(num, coeff[0xB]);
	//coeff[0xC]=_mm256_div_ps(num, coeff[0xC]);
	//coeff[0xD]=_mm256_div_ps(num, coeff[0xD]);
	//coeff[0xE]=_mm256_div_ps(num, coeff[0xE]);
	//coeff[0xF]=_mm256_div_ps(num, coeff[0xF]);
	
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
	__m256 g0=_mm256_set1_ps(240);
	__m256 g1=_mm256_set1_ps(240);
	__m256 g2=_mm256_set1_ps(180);
	__m256 g3=_mm256_set1_ps(180);
	__m256 g4=_mm256_set1_ps(140);
	__m256 g5=_mm256_set1_ps(160);
	__m256 g6=_mm256_set1_ps(120);
	__m256 g7=_mm256_set1_ps(120);

//	__m256 g0=_mm256_set1_ps(210);
// 	__m256 g1=g0;
//	__m256 g2=_mm256_set1_ps(215);
// 	__m256 g3=g2;
//	__m256 g4=_mm256_set1_ps(140);
//	__m256 g5=_mm256_set1_ps(230);
//	__m256 g6=_mm256_set1_ps(120);
//	__m256 g7=_mm256_set1_ps(140);
	coeff[0x0]=_mm256_mul_ps(coeff[0x0], g0);
	coeff[0x1]=_mm256_mul_ps(coeff[0x1], g0);
	coeff[0x2]=_mm256_mul_ps(coeff[0x2], g1);
	coeff[0x3]=_mm256_mul_ps(coeff[0x3], g1);
	coeff[0x4]=_mm256_mul_ps(coeff[0x4], g2);
	coeff[0x5]=_mm256_mul_ps(coeff[0x5], g2);
	coeff[0x6]=_mm256_mul_ps(coeff[0x6], g3);
	coeff[0x7]=_mm256_mul_ps(coeff[0x7], g3);
	coeff[0x8]=_mm256_mul_ps(coeff[0x8], g4);
	coeff[0x9]=_mm256_mul_ps(coeff[0x9], g4);
	coeff[0xA]=_mm256_mul_ps(coeff[0xA], g5);
	coeff[0xB]=_mm256_mul_ps(coeff[0xB], g5);
	coeff[0xC]=_mm256_mul_ps(coeff[0xC], g6);
	coeff[0xD]=_mm256_mul_ps(coeff[0xD], g6);
	coeff[0xE]=_mm256_mul_ps(coeff[0xE], g7);
	coeff[0xF]=_mm256_mul_ps(coeff[0xF], g7);

	//pr = (sum: coeff[i]*wg[i])/(sum: coeff[i])	NPREDS=8  (*16 interleaved floats)
	pr[0x0]=_mm256_mul_ps(pr[0x0], coeff[0x0]);
	pr[0x1]=_mm256_mul_ps(pr[0x1], coeff[0x1]);
	pr[0x2]=_mm256_mul_ps(pr[0x2], coeff[0x2]);
	pr[0x3]=_mm256_mul_ps(pr[0x3], coeff[0x3]);
	pr[0x4]=_mm256_mul_ps(pr[0x4], coeff[0x4]);
	pr[0x5]=_mm256_mul_ps(pr[0x5], coeff[0x5]);
	pr[0x6]=_mm256_mul_ps(pr[0x6], coeff[0x6]);
	pr[0x7]=_mm256_mul_ps(pr[0x7], coeff[0x7]);
	pr[0x8]=_mm256_mul_ps(pr[0x8], coeff[0x8]);
	pr[0x9]=_mm256_mul_ps(pr[0x9], coeff[0x9]);
	pr[0xA]=_mm256_mul_ps(pr[0xA], coeff[0xA]);
	pr[0xB]=_mm256_mul_ps(pr[0xB], coeff[0xB]);
	pr[0xC]=_mm256_mul_ps(pr[0xC], coeff[0xC]);
	pr[0xD]=_mm256_mul_ps(pr[0xD], coeff[0xD]);
	pr[0xE]=_mm256_mul_ps(pr[0xE], coeff[0xE]);
	pr[0xF]=_mm256_mul_ps(pr[0xF], coeff[0xF]);

	coeff[0x0]=_mm256_add_ps(coeff[0x0], coeff[0x2]);
	coeff[0x1]=_mm256_add_ps(coeff[0x1], coeff[0x3]);
	coeff[0x4]=_mm256_add_ps(coeff[0x4], coeff[0x6]);
	coeff[0x5]=_mm256_add_ps(coeff[0x5], coeff[0x7]);
	coeff[0x8]=_mm256_add_ps(coeff[0x8], coeff[0xA]);
	coeff[0x9]=_mm256_add_ps(coeff[0x9], coeff[0xB]);
	coeff[0xC]=_mm256_add_ps(coeff[0xC], coeff[0xE]);
	coeff[0xD]=_mm256_add_ps(coeff[0xD], coeff[0xF]);

	coeff[0x0]=_mm256_add_ps(coeff[0x0], coeff[0x4]);
	coeff[0x1]=_mm256_add_ps(coeff[0x1], coeff[0x5]);
	coeff[0x8]=_mm256_add_ps(coeff[0x8], coeff[0xC]);
	coeff[0x9]=_mm256_add_ps(coeff[0x9], coeff[0xD]);
	
	coeff[0x0]=_mm256_add_ps(coeff[0x0], coeff[0x8]);
	coeff[0x1]=_mm256_add_ps(coeff[0x1], coeff[0x9]);
	
	pr[0x0]=_mm256_add_ps(pr[0x0], pr[0x2]);
	pr[0x1]=_mm256_add_ps(pr[0x1], pr[0x3]);
	pr[0x4]=_mm256_add_ps(pr[0x4], pr[0x6]);
	pr[0x5]=_mm256_add_ps(pr[0x5], pr[0x7]);
	pr[0x8]=_mm256_add_ps(pr[0x8], pr[0xA]);
	pr[0x9]=_mm256_add_ps(pr[0x9], pr[0xB]);
	pr[0xC]=_mm256_add_ps(pr[0xC], pr[0xE]);
	pr[0xD]=_mm256_add_ps(pr[0xD], pr[0xF]);

	pr[0x0]=_mm256_add_ps(pr[0x0], pr[0x4]);
	pr[0x1]=_mm256_add_ps(pr[0x1], pr[0x5]);
	pr[0x8]=_mm256_add_ps(pr[0x8], pr[0xC]);
	pr[0x9]=_mm256_add_ps(pr[0x9], pr[0xD]);

	pr[0x0]=_mm256_add_ps(pr[0x0], pr[0x8]);
	pr[0x1]=_mm256_add_ps(pr[0x1], pr[0x9]);

	result[0]=_mm256_div_ps(pr[0], coeff[0]);//11/5
	result[1]=_mm256_div_ps(pr[1], coeff[1]);
#ifdef WG4_SERIALDEBUG
	for(int k=0;k<NCODERS/2;++k)
	{
		int it=0;
		float simdpred=((float*)result)[8*(k&1)+(k>>1)];
		double f2pred, w2sum, pred2;
		float f1pred, w1sum, pred1;
	again:
		f2pred=0, w2sum=0;
		f1pred=0, w1sum=0;
		for(int kp=0;kp<WG_NPREDS;++kp)
		{
			int subpred=((short*)wgpreds)[3*NCODERS*kp+k], err=((unsigned short*)prederrors)[NCODERS*kp+k];
			double w2=1./err;
			w2sum+=w2;
			f2pred+=w2*subpred;

			int k2=k>>1, odd=k&1;
			float fe=cvtpe[2*NCODERS*odd+8*kp+k2], fsp=cvtsp[2*NCODERS*odd+8*kp+k2];
			float w1=1.f/fe;
			w1sum+=w1;
			f1pred+=w1*fsp;
		}
		pred1=f1pred/w1sum;
		pred2=f2pred/w2sum;
		if(fabsf(pred1-simdpred)>5)
		{
			if(it++<2)
				goto again;//
			LOG_ERROR("");
		}
		if(fabs(pred2-simdpred)>5)
		{
			if(it++<2)
				goto again;//
			LOG_ERROR("");
		}
		//pred=(int)(fpred/wsum+0.5f);
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
	
	//errors = pI+pNW+2*pN+pNE+pNNE		(__m256i*)(erows[-Y]+((P+X*NPREDS)*3+C)*NCODERS)+R
	prederrors[0x0]=_mm256_load_si256((const __m256i*)(currNNptr+((0-1*WG_NPREDS)*3+0)*NCODERS)+0);//pW
	prederrors[0x1]=_mm256_load_si256((const __m256i*)(currNNptr+((0-1*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0x2]=_mm256_load_si256((const __m256i*)(currNNptr+((1-1*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0x3]=_mm256_load_si256((const __m256i*)(currNNptr+((1-1*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0x4]=_mm256_load_si256((const __m256i*)(currNNptr+((2-1*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0x5]=_mm256_load_si256((const __m256i*)(currNNptr+((2-1*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0x6]=_mm256_load_si256((const __m256i*)(currNNptr+((3-1*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0x7]=_mm256_load_si256((const __m256i*)(currNNptr+((3-1*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0x8]=_mm256_load_si256((const __m256i*)(currNNptr+((4-1*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0x9]=_mm256_load_si256((const __m256i*)(currNNptr+((4-1*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0xA]=_mm256_load_si256((const __m256i*)(currNNptr+((5-1*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0xB]=_mm256_load_si256((const __m256i*)(currNNptr+((5-1*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0xC]=_mm256_load_si256((const __m256i*)(currNNptr+((6-1*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0xD]=_mm256_load_si256((const __m256i*)(currNNptr+((6-1*WG_NPREDS)*3+0)*NCODERS)+1);
	prederrors[0xE]=_mm256_load_si256((const __m256i*)(currNNptr+((7-1*WG_NPREDS)*3+0)*NCODERS)+0);
	prederrors[0xF]=_mm256_load_si256((const __m256i*)(currNNptr+((7-1*WG_NPREDS)*3+0)*NCODERS)+1);

	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], _mm256_load_si256((const __m256i*)(Nptr+((0+0*WG_NPREDS)*3+0)*NCODERS)+0));//pN
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], _mm256_load_si256((const __m256i*)(Nptr+((0+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], _mm256_load_si256((const __m256i*)(Nptr+((1+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], _mm256_load_si256((const __m256i*)(Nptr+((1+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], _mm256_load_si256((const __m256i*)(Nptr+((2+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], _mm256_load_si256((const __m256i*)(Nptr+((2+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], _mm256_load_si256((const __m256i*)(Nptr+((3+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], _mm256_load_si256((const __m256i*)(Nptr+((3+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], _mm256_load_si256((const __m256i*)(Nptr+((4+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], _mm256_load_si256((const __m256i*)(Nptr+((4+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], _mm256_load_si256((const __m256i*)(Nptr+((5+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], _mm256_load_si256((const __m256i*)(Nptr+((5+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], _mm256_load_si256((const __m256i*)(Nptr+((6+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], _mm256_load_si256((const __m256i*)(Nptr+((6+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], _mm256_load_si256((const __m256i*)(Nptr+((7+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], _mm256_load_si256((const __m256i*)(Nptr+((7+0*WG_NPREDS)*3+0)*NCODERS)+1));

	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], _mm256_load_si256((const __m256i*)(Nptr+((0-1*WG_NPREDS)*3+0)*NCODERS)+0));//pNW
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], _mm256_load_si256((const __m256i*)(Nptr+((0-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], _mm256_load_si256((const __m256i*)(Nptr+((1-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], _mm256_load_si256((const __m256i*)(Nptr+((1-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], _mm256_load_si256((const __m256i*)(Nptr+((2-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], _mm256_load_si256((const __m256i*)(Nptr+((2-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], _mm256_load_si256((const __m256i*)(Nptr+((3-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], _mm256_load_si256((const __m256i*)(Nptr+((3-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], _mm256_load_si256((const __m256i*)(Nptr+((4-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], _mm256_load_si256((const __m256i*)(Nptr+((4-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], _mm256_load_si256((const __m256i*)(Nptr+((5-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], _mm256_load_si256((const __m256i*)(Nptr+((5-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], _mm256_load_si256((const __m256i*)(Nptr+((6-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], _mm256_load_si256((const __m256i*)(Nptr+((6-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], _mm256_load_si256((const __m256i*)(Nptr+((7-1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], _mm256_load_si256((const __m256i*)(Nptr+((7-1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x0]=_mm256_slli_epi16(prederrors[0x0], 1);
	prederrors[0x1]=_mm256_slli_epi16(prederrors[0x1], 1);
	prederrors[0x2]=_mm256_slli_epi16(prederrors[0x2], 1);
	prederrors[0x3]=_mm256_slli_epi16(prederrors[0x3], 1);
	prederrors[0x4]=_mm256_slli_epi16(prederrors[0x4], 1);
	prederrors[0x5]=_mm256_slli_epi16(prederrors[0x5], 1);
	prederrors[0x6]=_mm256_slli_epi16(prederrors[0x6], 1);
	prederrors[0x7]=_mm256_slli_epi16(prederrors[0x7], 1);
	prederrors[0x8]=_mm256_slli_epi16(prederrors[0x8], 1);
	prederrors[0x9]=_mm256_slli_epi16(prederrors[0x9], 1);
	prederrors[0xA]=_mm256_slli_epi16(prederrors[0xA], 1);
	prederrors[0xB]=_mm256_slli_epi16(prederrors[0xB], 1);
	prederrors[0xC]=_mm256_slli_epi16(prederrors[0xC], 1);
	prederrors[0xD]=_mm256_slli_epi16(prederrors[0xD], 1);
	prederrors[0xE]=_mm256_slli_epi16(prederrors[0xE], 1);
	prederrors[0xF]=_mm256_slli_epi16(prederrors[0xF], 1);

	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], _mm256_load_si256((__m256i*)(currNNptr+((0+0*WG_NPREDS)*3+0)*NCODERS)+0));//pNN
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], _mm256_load_si256((__m256i*)(currNNptr+((0+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], _mm256_load_si256((__m256i*)(currNNptr+((1+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], _mm256_load_si256((__m256i*)(currNNptr+((1+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], _mm256_load_si256((__m256i*)(currNNptr+((2+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], _mm256_load_si256((__m256i*)(currNNptr+((2+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], _mm256_load_si256((__m256i*)(currNNptr+((3+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], _mm256_load_si256((__m256i*)(currNNptr+((3+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], _mm256_load_si256((__m256i*)(currNNptr+((4+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], _mm256_load_si256((__m256i*)(currNNptr+((4+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], _mm256_load_si256((__m256i*)(currNNptr+((5+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], _mm256_load_si256((__m256i*)(currNNptr+((5+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], _mm256_load_si256((__m256i*)(currNNptr+((6+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], _mm256_load_si256((__m256i*)(currNNptr+((6+0*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], _mm256_load_si256((__m256i*)(currNNptr+((7+0*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], _mm256_load_si256((__m256i*)(currNNptr+((7+0*WG_NPREDS)*3+0)*NCODERS)+1));

	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], _mm256_load_si256((__m256i*)(currNNptr+((0+1*WG_NPREDS)*3+0)*NCODERS)+0));//pNNE
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], _mm256_load_si256((__m256i*)(currNNptr+((0+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], _mm256_load_si256((__m256i*)(currNNptr+((1+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], _mm256_load_si256((__m256i*)(currNNptr+((1+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], _mm256_load_si256((__m256i*)(currNNptr+((2+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], _mm256_load_si256((__m256i*)(currNNptr+((2+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], _mm256_load_si256((__m256i*)(currNNptr+((3+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], _mm256_load_si256((__m256i*)(currNNptr+((3+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], _mm256_load_si256((__m256i*)(currNNptr+((4+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], _mm256_load_si256((__m256i*)(currNNptr+((4+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], _mm256_load_si256((__m256i*)(currNNptr+((5+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], _mm256_load_si256((__m256i*)(currNNptr+((5+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], _mm256_load_si256((__m256i*)(currNNptr+((6+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], _mm256_load_si256((__m256i*)(currNNptr+((6+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], _mm256_load_si256((__m256i*)(currNNptr+((7+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], _mm256_load_si256((__m256i*)(currNNptr+((7+1*WG_NPREDS)*3+0)*NCODERS)+1));

	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], _mm256_load_si256((__m256i*)(Nptr+((0+1*WG_NPREDS)*3+0)*NCODERS)+0));//pNE
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], _mm256_load_si256((__m256i*)(Nptr+((0+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], _mm256_load_si256((__m256i*)(Nptr+((1+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], _mm256_load_si256((__m256i*)(Nptr+((1+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], _mm256_load_si256((__m256i*)(Nptr+((2+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], _mm256_load_si256((__m256i*)(Nptr+((2+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], _mm256_load_si256((__m256i*)(Nptr+((3+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], _mm256_load_si256((__m256i*)(Nptr+((3+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], _mm256_load_si256((__m256i*)(Nptr+((4+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], _mm256_load_si256((__m256i*)(Nptr+((4+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], _mm256_load_si256((__m256i*)(Nptr+((5+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], _mm256_load_si256((__m256i*)(Nptr+((5+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], _mm256_load_si256((__m256i*)(Nptr+((6+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], _mm256_load_si256((__m256i*)(Nptr+((6+1*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], _mm256_load_si256((__m256i*)(Nptr+((7+1*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], _mm256_load_si256((__m256i*)(Nptr+((7+1*WG_NPREDS)*3+0)*NCODERS)+1));

	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], _mm256_load_si256((__m256i*)(currNNptr+((0-2*WG_NPREDS)*3+0)*NCODERS)+0));//pWW
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], _mm256_load_si256((__m256i*)(currNNptr+((0-2*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], _mm256_load_si256((__m256i*)(currNNptr+((1-2*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], _mm256_load_si256((__m256i*)(currNNptr+((1-2*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], _mm256_load_si256((__m256i*)(currNNptr+((2-2*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], _mm256_load_si256((__m256i*)(currNNptr+((2-2*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], _mm256_load_si256((__m256i*)(currNNptr+((3-2*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], _mm256_load_si256((__m256i*)(currNNptr+((3-2*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], _mm256_load_si256((__m256i*)(currNNptr+((4-2*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], _mm256_load_si256((__m256i*)(currNNptr+((4-2*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], _mm256_load_si256((__m256i*)(currNNptr+((5-2*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], _mm256_load_si256((__m256i*)(currNNptr+((5-2*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], _mm256_load_si256((__m256i*)(currNNptr+((6-2*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], _mm256_load_si256((__m256i*)(currNNptr+((6-2*WG_NPREDS)*3+0)*NCODERS)+1));
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], _mm256_load_si256((__m256i*)(currNNptr+((7-2*WG_NPREDS)*3+0)*NCODERS)+0));
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], _mm256_load_si256((__m256i*)(currNNptr+((7-2*WG_NPREDS)*3+0)*NCODERS)+1));
#ifndef WG_ENABLE_eW
	(void)chWerrors;
#else
	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(0*3+0)*NCODERS)+0), 2));//pI	(__m256i*)(wgWerrors+(P*3+C)*NCODERS)+R
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(0*3+0)*NCODERS)+1), 2));
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(1*3+0)*NCODERS)+0), 2));
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(1*3+0)*NCODERS)+1), 2));
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(2*3+0)*NCODERS)+0), 2));
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(2*3+0)*NCODERS)+1), 2));
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(3*3+0)*NCODERS)+0), 2));
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(3*3+0)*NCODERS)+1), 2));
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(4*3+0)*NCODERS)+0), 2));
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(4*3+0)*NCODERS)+1), 2));
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(5*3+0)*NCODERS)+0), 2));
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(5*3+0)*NCODERS)+1), 2));
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(6*3+0)*NCODERS)+0), 2));
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(6*3+0)*NCODERS)+1), 2));
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(7*3+0)*NCODERS)+0), 2));
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], _mm256_srli_epi16(_mm256_load_si256((const __m256i*)(chWerrors+(7*3+0)*NCODERS)+1), 2));
#endif

#ifdef INTMIX
	__m256i ires[4];
	wg_mix_pt2(wgpreds+(ptrdiff_t)2*kc+0, prederrors+0, ires+0);//P*6+2*C+R
	wg_mix_pt2(wgpreds+(ptrdiff_t)2*kc+1, prederrors+1, ires+2);
	ires[0]=_mm256_slli_epi32(ires[0], 16);
	ires[2]=_mm256_slli_epi32(ires[2], 16);
	ires[1]=_mm256_slli_epi32(ires[1], 16);
	ires[3]=_mm256_slli_epi32(ires[3], 16);
	ires[0]=_mm256_srli_epi32(ires[0], 16);
	ires[2]=_mm256_srli_epi32(ires[2], 16);
	preds[0]=_mm256_or_si256(ires[0], ires[1]);
	preds[1]=_mm256_or_si256(ires[2], ires[3]);
#else
	__m256i one=_mm256_set1_epi16(1);
	prederrors[0x0]=_mm256_add_epi16(prederrors[0x0], one);
	prederrors[0x1]=_mm256_add_epi16(prederrors[0x1], one);
	prederrors[0x2]=_mm256_add_epi16(prederrors[0x2], one);
	prederrors[0x3]=_mm256_add_epi16(prederrors[0x3], one);
	prederrors[0x4]=_mm256_add_epi16(prederrors[0x4], one);
	prederrors[0x5]=_mm256_add_epi16(prederrors[0x5], one);
	prederrors[0x6]=_mm256_add_epi16(prederrors[0x6], one);
	prederrors[0x7]=_mm256_add_epi16(prederrors[0x7], one);
	prederrors[0x8]=_mm256_add_epi16(prederrors[0x8], one);
	prederrors[0x9]=_mm256_add_epi16(prederrors[0x9], one);
	prederrors[0xA]=_mm256_add_epi16(prederrors[0xA], one);
	prederrors[0xB]=_mm256_add_epi16(prederrors[0xB], one);
	prederrors[0xC]=_mm256_add_epi16(prederrors[0xC], one);
	prederrors[0xD]=_mm256_add_epi16(prederrors[0xD], one);
	prederrors[0xE]=_mm256_add_epi16(prederrors[0xE], one);
	prederrors[0xF]=_mm256_add_epi16(prederrors[0xF], one);
	__m256 result[4];
	wg_mix_pt2f(wgpreds+(ptrdiff_t)2*kc+0, prederrors+0, result+0);//P*6+2*C+R
	wg_mix_pt2f(wgpreds+(ptrdiff_t)2*kc+1, prederrors+1, result+2);

	__m256i ires[4];
	ires[0]=_mm256_cvtps_epi32(result[0]);
	ires[1]=_mm256_cvtps_epi32(result[1]);
	ires[2]=_mm256_cvtps_epi32(result[2]);
	ires[3]=_mm256_cvtps_epi32(result[3]);
	ires[0]=_mm256_slli_epi32(ires[0], 16);
	ires[2]=_mm256_slli_epi32(ires[2], 16);
	ires[1]=_mm256_slli_epi32(ires[1], 16);
	ires[3]=_mm256_slli_epi32(ires[3], 16);
	ires[0]=_mm256_srli_epi32(ires[0], 16);
	ires[2]=_mm256_srli_epi32(ires[2], 16);
	preds[0]=_mm256_or_si256(ires[0], ires[1]);
	preds[1]=_mm256_or_si256(ires[2], ires[3]);
#endif
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
		int cdfW=0;
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

			cdfW+=freq;
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
AWM_INLINE void dec_yuv(__m256i *mstate, int kc, const __m256i *ctx0, const __m256i *ctx1, const int *CDF2syms, const int *ans_permute, unsigned char **pstreamptr, const unsigned char *streamend, __m256i *syms)
{
	const unsigned char *streamptr=*pstreamptr;
	__m256i decctx[4];
	{
		__m128i t0=_mm256_extracti128_si256(*ctx0, 1);
		__m128i t1=_mm256_extracti128_si256(*ctx1, 1);
		decctx[0]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(*ctx0));
		decctx[2]=_mm256_cvtepi16_epi32(_mm256_castsi256_si128(*ctx1));
		decctx[1]=_mm256_cvtepi16_epi32(t0);
		decctx[3]=_mm256_cvtepi16_epi32(t1);
	}
	decctx[0]=_mm256_slli_epi32(decctx[0], PROBBITS);
	decctx[1]=_mm256_slli_epi32(decctx[1], PROBBITS);
	decctx[2]=_mm256_slli_epi32(decctx[2], PROBBITS);
	decctx[3]=_mm256_slli_epi32(decctx[3], PROBBITS);
#ifdef ANS_VAL
	ansval_check(mstate, sizeof(int), NCODERS);
#endif
	{
		__m256i mprobmask=_mm256_set1_epi32((1<<PROBBITS)-1);
		__m256i rem0=_mm256_and_si256(mstate[0], mprobmask);
		__m256i rem1=_mm256_and_si256(mstate[1], mprobmask);
		__m256i rem2=_mm256_and_si256(mstate[2], mprobmask);
		__m256i rem3=_mm256_and_si256(mstate[3], mprobmask);
		decctx[0]=_mm256_or_si256(decctx[0], rem0);
		decctx[1]=_mm256_or_si256(decctx[1], rem1);
		decctx[2]=_mm256_or_si256(decctx[2], rem2);
		decctx[3]=_mm256_or_si256(decctx[3], rem3);
	}
#ifdef ANS_VAL
	ALIGN(32) int debugctx[NCODERS];
	memcpy(debugctx, decctx, sizeof(int[NCODERS]));
#endif
	const int *statsptr=CDF2syms+((ptrdiff_t)NCTX*kc<<PROBBITS);
	gather32((int*)(decctx+0), statsptr, (int*)(decctx+0));
	gather32((int*)(decctx+1), statsptr, (int*)(decctx+1));
	gather32((int*)(decctx+2), statsptr, (int*)(decctx+2));
	gather32((int*)(decctx+3), statsptr, (int*)(decctx+3));
	//decctx[0]=_mm256_i32gather_epi32(statsptr, decctx[0], sizeof(int));//FIXME try not using gather
	//decctx[1]=_mm256_i32gather_epi32(statsptr, decctx[1], sizeof(int));
	//decctx[2]=_mm256_i32gather_epi32(statsptr, decctx[2], sizeof(int));
	//decctx[3]=_mm256_i32gather_epi32(statsptr, decctx[3], sizeof(int));

	//update		state = (state>>12)*freq+(rem-cdf)	rem-cdf is prebaked
	{
		__m256i mfreq0=_mm256_srli_epi32(decctx[0], PROBBITS+8);//1 <= freq <= 0xF01
		__m256i mfreq1=_mm256_srli_epi32(decctx[1], PROBBITS+8);
		__m256i mfreq2=_mm256_srli_epi32(decctx[2], PROBBITS+8);
		__m256i mfreq3=_mm256_srli_epi32(decctx[3], PROBBITS+8);
#ifdef ANS_VAL
		__m256i mdebugfreq[2];
		mdebugfreq[0]=_mm256_packus_epi32(mfreq0, mfreq1);
		mdebugfreq[1]=_mm256_packus_epi32(mfreq2, mfreq3);
		mdebugfreq[0]=_mm256_permute4x64_epi64(mdebugfreq[0], _MM_SHUFFLE(3, 1, 2, 0));
		mdebugfreq[1]=_mm256_permute4x64_epi64(mdebugfreq[1], _MM_SHUFFLE(3, 1, 2, 0));
		ALIGN(32) unsigned short freqs[NCODERS];
		memcpy(freqs, mdebugfreq, sizeof(freqs));
		ansval_check(freqs, sizeof(short), NCODERS);
#endif
		mstate[0]=_mm256_srli_epi32(mstate[0], PROBBITS);
		mstate[1]=_mm256_srli_epi32(mstate[1], PROBBITS);
		mstate[2]=_mm256_srli_epi32(mstate[2], PROBBITS);
		mstate[3]=_mm256_srli_epi32(mstate[3], PROBBITS);
		mstate[0]=_mm256_mullo_epi32(mstate[0], mfreq0);//10 cycles
		mstate[1]=_mm256_mullo_epi32(mstate[1], mfreq1);
		mstate[2]=_mm256_mullo_epi32(mstate[2], mfreq2);
		mstate[3]=_mm256_mullo_epi32(mstate[3], mfreq3);
	}
	//__m256i lobyte=_mm256_set1_epi32(255);
	//__m256i msym0=_mm256_and_si256(decctx[0], lobyte);
	//__m256i msym1=_mm256_and_si256(decctx[1], lobyte);
	//__m256i msym2=_mm256_and_si256(decctx[2], lobyte);
	//__m256i msym3=_mm256_and_si256(decctx[3], lobyte);
	//msym0=_mm256_packus_epi32(msym0, msym1);
	//msym1=_mm256_packus_epi32(msym2, msym3);
	//msym0=_mm256_permute4x64_epi64(msym0, _MM_SHUFFLE(3, 1, 2, 0));
	//msym1=_mm256_permute4x64_epi64(msym1, _MM_SHUFFLE(3, 1, 2, 0));
	{
		__m256i mbias0=_mm256_slli_epi32(decctx[0], PROBBITS);
		__m256i mbias1=_mm256_slli_epi32(decctx[1], PROBBITS);
		__m256i mbias2=_mm256_slli_epi32(decctx[2], PROBBITS);
		__m256i mbias3=_mm256_slli_epi32(decctx[3], PROBBITS);
		mbias0=_mm256_srli_epi32(mbias0, 32-PROBBITS);
		mbias1=_mm256_srli_epi32(mbias1, 32-PROBBITS);
		mbias2=_mm256_srli_epi32(mbias2, 32-PROBBITS);
		mbias3=_mm256_srli_epi32(mbias3, 32-PROBBITS);
		mstate[0]=_mm256_add_epi32(mstate[0], mbias0);
		mstate[1]=_mm256_add_epi32(mstate[1], mbias1);
		mstate[2]=_mm256_add_epi32(mstate[2], mbias2);
		mstate[3]=_mm256_add_epi32(mstate[3], mbias3);
	}
	//decctx[0]=_mm256_slli_epi32(decctx[0], 24);
	//decctx[1]=_mm256_slli_epi32(decctx[1], 24);
	//decctx[2]=_mm256_slli_epi32(decctx[2], 24);
	//decctx[3]=_mm256_slli_epi32(decctx[3], 24);
	//decctx[0]=_mm256_srai_epi32(decctx[0], 24);
	//decctx[1]=_mm256_srai_epi32(decctx[1], 24);
	//decctx[2]=_mm256_srai_epi32(decctx[2], 24);
	//decctx[3]=_mm256_srai_epi32(decctx[3], 24);
	//decctx[0]=_mm256_packs_epi16(decctx[0], decctx[1]);
	//decctx[2]=_mm256_packs_epi16(decctx[2], decctx[3]);
	__m256i symmask=_mm256_set1_epi32(255);
	decctx[0]=_mm256_and_si256(decctx[0], symmask);
	decctx[1]=_mm256_and_si256(decctx[1], symmask);
	decctx[2]=_mm256_and_si256(decctx[2], symmask);
	decctx[3]=_mm256_and_si256(decctx[3], symmask);
	decctx[0]=_mm256_packus_epi16(decctx[0], decctx[1]);
	decctx[2]=_mm256_packus_epi16(decctx[2], decctx[3]);
	syms[0]=_mm256_permute4x64_epi64(decctx[0], _MM_SHUFFLE(3, 1, 2, 0));
	syms[1]=_mm256_permute4x64_epi64(decctx[2], _MM_SHUFFLE(3, 1, 2, 0));
#ifdef ANS_VAL
	ansval_check(mstate, sizeof(int), NCODERS);
#endif
	//renorm
	{
#ifdef _DEBUG
		if(streamptr>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__m256i lo0=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr));
		__m256i smin=_mm256_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
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
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__m256i lo1=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr));
		streamptr+=mask1*sizeof(short);
		lo1=_mm256_permutevar8x32_epi32(lo1, idx1);
		renorm0=_mm256_or_si256(renorm0, lo0);
		renorm1=_mm256_or_si256(renorm1, lo1);
		
#ifdef _DEBUG
		if(streamptr>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__m256i lo2=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr));
		__m256i cond2=_mm256_cmpgt_epi32(smin, mstate[2]);
		__m256i cond3=_mm256_cmpgt_epi32(smin, mstate[3]);
		int mask2=_mm256_movemask_ps(_mm256_castsi256_ps(cond2));
		int mask3=_mm256_movemask_ps(_mm256_castsi256_ps(cond3));
		__m256i idx2=_mm256_load_si256((const __m256i*)ans_permute+mask2);
		__m256i idx3=_mm256_load_si256((const __m256i*)ans_permute+mask3);
		mask2=_mm_popcnt_u32(mask2);
		mask3=_mm_popcnt_u32(mask3);
		__m256i renorm2=_mm256_slli_epi32(mstate[2], 16);
		__m256i renorm3=_mm256_slli_epi32(mstate[3], 16);
		streamptr+=mask2*sizeof(short);
		lo2=_mm256_permutevar8x32_epi32(lo2, idx2);
#ifdef _DEBUG
		if(streamptr>streamend)
			LOG_ERROR("OOB ptr %016zX >= %016zX", streamptr, streamend);
#endif
		__m256i lo3=_mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i*)streamptr));
		streamptr+=mask3*sizeof(short);
		lo3=_mm256_permutevar8x32_epi32(lo3, idx3);
		renorm2=_mm256_or_si256(renorm2, lo2);
		renorm3=_mm256_or_si256(renorm3, lo3);

		mstate[0]=_mm256_blendv_epi8(mstate[0], renorm0, cond0);
		mstate[1]=_mm256_blendv_epi8(mstate[1], renorm1, cond1);
		mstate[2]=_mm256_blendv_epi8(mstate[2], renorm2, cond2);
		mstate[3]=_mm256_blendv_epi8(mstate[3], renorm3, cond3);
	}
	*pstreamptr=(unsigned char*)(size_t)streamptr;
}
AWM_INLINE void transpose16x2(__m256i *data)
{
	//https://pzemtsov.github.io/2014/10/01/how-to-transpose-a-16x16-matrix.html
	__m256i t0, t1, t2, t3;
#define SHUFFLE(LO, HI, IMM8_HHLL) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(LO), _mm256_castsi256_ps(HI), IMM8_HHLL))
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
	t0=_mm256_set_epi8(
		15, 11,  7,  3,
		14, 10,  6,  2,
		13,  9,  5,  1,
		12,  8,  4,  0,
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
	data[0x0]=_mm256_shuffle_epi8(data[0x0], t0);
	data[0x5]=_mm256_shuffle_epi8(data[0x5], t0);
	data[0xA]=_mm256_shuffle_epi8(data[0xA], t0);
	data[0xF]=_mm256_shuffle_epi8(data[0xF], t0);
	//nondiagonals
	t1=data[0x4]; data[0x4]=_mm256_shuffle_epi8(data[0x1], t0); data[0x1]=_mm256_shuffle_epi8(t1, t0);
	t1=data[0x8]; data[0x8]=_mm256_shuffle_epi8(data[0x2], t0); data[0x2]=_mm256_shuffle_epi8(t1, t0);
	t1=data[0xC]; data[0xC]=_mm256_shuffle_epi8(data[0x3], t0); data[0x3]=_mm256_shuffle_epi8(t1, t0);
	t1=data[0x9]; data[0x9]=_mm256_shuffle_epi8(data[0x6], t0); data[0x6]=_mm256_shuffle_epi8(t1, t0);
	t1=data[0xD]; data[0xD]=_mm256_shuffle_epi8(data[0x7], t0); data[0x7]=_mm256_shuffle_epi8(t1, t0);
	t1=data[0xE]; data[0xE]=_mm256_shuffle_epi8(data[0xB], t0); data[0xB]=_mm256_shuffle_epi8(t1, t0);
	TRANSPOSE4(data[0x0], data[0x1], data[0x2], data[0x3],  t0, t1, t2, t3,  data[0x0], data[0x1], data[0x2], data[0x3]);
	TRANSPOSE4(data[0x4], data[0x5], data[0x6], data[0x7],  t0, t1, t2, t3,  data[0x4], data[0x5], data[0x6], data[0x7]);
	TRANSPOSE4(data[0x8], data[0x9], data[0xA], data[0xB],  t0, t1, t2, t3,  data[0x8], data[0x9], data[0xA], data[0xB]);
	TRANSPOSE4(data[0xC], data[0xD], data[0xE], data[0xF],  t0, t1, t2, t3,  data[0xC], data[0xD], data[0xE], data[0xF]);
#undef  SHUFFLE
#undef  TRANSPOSE4
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
	int SIMDxcount=blockxbytes&~((int)sizeof(__m256i[NCODERS])-1);
	__m256i slowinc=_mm256_set1_epi64x(sizeof(__m256i));
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
		for(;kx<SIMDxcount;kx+=(int)sizeof(__m256i[NCODERS]))
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
			__m256i block[16];
			block[0x0]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x10]+0, (__m128i*)slowptrs[0x00]+0);
			block[0x1]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x11]+0, (__m128i*)slowptrs[0x01]+0);
			block[0x2]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x12]+0, (__m128i*)slowptrs[0x02]+0);
			block[0x3]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x13]+0, (__m128i*)slowptrs[0x03]+0);
			block[0x4]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x14]+0, (__m128i*)slowptrs[0x04]+0);
			block[0x5]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x15]+0, (__m128i*)slowptrs[0x05]+0);
			block[0x6]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x16]+0, (__m128i*)slowptrs[0x06]+0);
			block[0x7]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x17]+0, (__m128i*)slowptrs[0x07]+0);
			block[0x8]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x18]+0, (__m128i*)slowptrs[0x08]+0);
			block[0x9]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x19]+0, (__m128i*)slowptrs[0x09]+0);
			block[0xA]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x1A]+0, (__m128i*)slowptrs[0x0A]+0);
			block[0xB]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x1B]+0, (__m128i*)slowptrs[0x0B]+0);
			block[0xC]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x1C]+0, (__m128i*)slowptrs[0x0C]+0);
			block[0xD]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x1D]+0, (__m128i*)slowptrs[0x0D]+0);
			block[0xE]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x1E]+0, (__m128i*)slowptrs[0x0E]+0);
			block[0xF]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x1F]+0, (__m128i*)slowptrs[0x0F]+0);
			transpose16x2(block);
			_mm256_store_si256((__m256i*)fastptr+0x00, block[0x0]);
			_mm256_store_si256((__m256i*)fastptr+0x01, block[0x1]);
			_mm256_store_si256((__m256i*)fastptr+0x02, block[0x2]);
			_mm256_store_si256((__m256i*)fastptr+0x03, block[0x3]);
			_mm256_store_si256((__m256i*)fastptr+0x04, block[0x4]);
			_mm256_store_si256((__m256i*)fastptr+0x05, block[0x5]);
			_mm256_store_si256((__m256i*)fastptr+0x06, block[0x6]);
			_mm256_store_si256((__m256i*)fastptr+0x07, block[0x7]);
			_mm256_store_si256((__m256i*)fastptr+0x08, block[0x8]);
			_mm256_store_si256((__m256i*)fastptr+0x09, block[0x9]);
			_mm256_store_si256((__m256i*)fastptr+0x0A, block[0xA]);
			_mm256_store_si256((__m256i*)fastptr+0x0B, block[0xB]);
			_mm256_store_si256((__m256i*)fastptr+0x0C, block[0xC]);
			_mm256_store_si256((__m256i*)fastptr+0x0D, block[0xD]);
			_mm256_store_si256((__m256i*)fastptr+0x0E, block[0xE]);
			_mm256_store_si256((__m256i*)fastptr+0x0F, block[0xF]);
			block[0x0]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x10]+1, (__m128i*)slowptrs[0x00]+1);
			block[0x1]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x11]+1, (__m128i*)slowptrs[0x01]+1);
			block[0x2]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x12]+1, (__m128i*)slowptrs[0x02]+1);
			block[0x3]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x13]+1, (__m128i*)slowptrs[0x03]+1);
			block[0x4]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x14]+1, (__m128i*)slowptrs[0x04]+1);
			block[0x5]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x15]+1, (__m128i*)slowptrs[0x05]+1);
			block[0x6]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x16]+1, (__m128i*)slowptrs[0x06]+1);
			block[0x7]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x17]+1, (__m128i*)slowptrs[0x07]+1);
			block[0x8]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x18]+1, (__m128i*)slowptrs[0x08]+1);
			block[0x9]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x19]+1, (__m128i*)slowptrs[0x09]+1);
			block[0xA]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x1A]+1, (__m128i*)slowptrs[0x0A]+1);
			block[0xB]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x1B]+1, (__m128i*)slowptrs[0x0B]+1);
			block[0xC]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x1C]+1, (__m128i*)slowptrs[0x0C]+1);
			block[0xD]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x1D]+1, (__m128i*)slowptrs[0x0D]+1);
			block[0xE]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x1E]+1, (__m128i*)slowptrs[0x0E]+1);
			block[0xF]=_mm256_loadu2_m128i((__m128i*)slowptrs[0x1F]+1, (__m128i*)slowptrs[0x0F]+1);
			transpose16x2(block);
			_mm256_store_si256((__m256i*)fastptr+0x10, block[0x0]);
			_mm256_store_si256((__m256i*)fastptr+0x11, block[0x1]);
			_mm256_store_si256((__m256i*)fastptr+0x12, block[0x2]);
			_mm256_store_si256((__m256i*)fastptr+0x13, block[0x3]);
			_mm256_store_si256((__m256i*)fastptr+0x14, block[0x4]);
			_mm256_store_si256((__m256i*)fastptr+0x15, block[0x5]);
			_mm256_store_si256((__m256i*)fastptr+0x16, block[0x6]);
			_mm256_store_si256((__m256i*)fastptr+0x17, block[0x7]);
			_mm256_store_si256((__m256i*)fastptr+0x18, block[0x8]);
			_mm256_store_si256((__m256i*)fastptr+0x19, block[0x9]);
			_mm256_store_si256((__m256i*)fastptr+0x1A, block[0xA]);
			_mm256_store_si256((__m256i*)fastptr+0x1B, block[0xB]);
			_mm256_store_si256((__m256i*)fastptr+0x1C, block[0xC]);
			_mm256_store_si256((__m256i*)fastptr+0x1D, block[0xD]);
			_mm256_store_si256((__m256i*)fastptr+0x1E, block[0xE]);
			_mm256_store_si256((__m256i*)fastptr+0x1F, block[0xF]);

			fastptr+=sizeof(__m256i[NCODERS]);
			block[0]=_mm256_load_si256((__m256i*)slowptrs+0);
			block[1]=_mm256_load_si256((__m256i*)slowptrs+1);
			block[2]=_mm256_load_si256((__m256i*)slowptrs+2);
			block[3]=_mm256_load_si256((__m256i*)slowptrs+3);
			block[4]=_mm256_load_si256((__m256i*)slowptrs+4);
			block[5]=_mm256_load_si256((__m256i*)slowptrs+5);
			block[6]=_mm256_load_si256((__m256i*)slowptrs+6);
			block[7]=_mm256_load_si256((__m256i*)slowptrs+7);
			block[0]=_mm256_add_epi64(block[0], slowinc);
			block[1]=_mm256_add_epi64(block[1], slowinc);
			block[2]=_mm256_add_epi64(block[2], slowinc);
			block[3]=_mm256_add_epi64(block[3], slowinc);
			block[4]=_mm256_add_epi64(block[4], slowinc);
			block[5]=_mm256_add_epi64(block[5], slowinc);
			block[6]=_mm256_add_epi64(block[6], slowinc);
			block[7]=_mm256_add_epi64(block[7], slowinc);
			_mm256_store_si256((__m256i*)slowptrs+0, block[0]);
			_mm256_store_si256((__m256i*)slowptrs+1, block[1]);
			_mm256_store_si256((__m256i*)slowptrs+2, block[2]);
			_mm256_store_si256((__m256i*)slowptrs+3, block[3]);
			_mm256_store_si256((__m256i*)slowptrs+4, block[4]);
			_mm256_store_si256((__m256i*)slowptrs+5, block[5]);
			_mm256_store_si256((__m256i*)slowptrs+6, block[6]);
			_mm256_store_si256((__m256i*)slowptrs+7, block[7]);
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
			fastptr[0x10]=*slowptrs[0x10]++;
			fastptr[0x11]=*slowptrs[0x11]++;
			fastptr[0x12]=*slowptrs[0x12]++;
			fastptr[0x13]=*slowptrs[0x13]++;
			fastptr[0x14]=*slowptrs[0x14]++;
			fastptr[0x15]=*slowptrs[0x15]++;
			fastptr[0x16]=*slowptrs[0x16]++;
			fastptr[0x17]=*slowptrs[0x17]++;
			fastptr[0x18]=*slowptrs[0x18]++;
			fastptr[0x19]=*slowptrs[0x19]++;
			fastptr[0x1A]=*slowptrs[0x1A]++;
			fastptr[0x1B]=*slowptrs[0x1B]++;
			fastptr[0x1C]=*slowptrs[0x1C]++;
			fastptr[0x1D]=*slowptrs[0x1D]++;
			fastptr[0x1E]=*slowptrs[0x1E]++;
			fastptr[0x1F]=*slowptrs[0x1F]++;
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

	//only difference between inv and fwd:		swap assignments (slow<-const fast)
	int rowstride=3*iw;
	int ixyblockw=iw/XCODERS;
	int ixyblockh=ih/YCODERS;
	int blockxbytes=3*NCODERS*ixyblockw;
#ifdef INTERLEAVESIMD
	int SIMDxcount=blockxbytes&~((int)sizeof(__m256i[NCODERS])-1);
	__m256i slowinc=_mm256_set1_epi64x(sizeof(__m256i));
#endif
	const unsigned char *fastptr=interleaved;
	unsigned char *slowptrs[NCODERS]={0}, *slowptrs0[NCODERS]={0};
	for(int ky=0;ky<YCODERS;++ky)//spread slow pointers
	{
		for(int kx=0;kx<XCODERS;++kx)
			slowptrs0[XCODERS*ky+kx]=original+3*(iw*ixyblockh*ky+ixyblockw*kx);
	}
	for(int ky=0;ky<ixyblockh;++ky)//deinterleave
	{
		int kx=0;
		memcpy(slowptrs, slowptrs0, sizeof(slowptrs));
#ifdef INTERLEAVESIMD
		for(;kx<SIMDxcount;kx+=(int)sizeof(__m256i[NCODERS]))
		{
			__m256i block[16];
			block[0x0]=_mm256_load_si256((__m256i*)fastptr+0x00);
			block[0x1]=_mm256_load_si256((__m256i*)fastptr+0x01);
			block[0x2]=_mm256_load_si256((__m256i*)fastptr+0x02);
			block[0x3]=_mm256_load_si256((__m256i*)fastptr+0x03);
			block[0x4]=_mm256_load_si256((__m256i*)fastptr+0x04);
			block[0x5]=_mm256_load_si256((__m256i*)fastptr+0x05);
			block[0x6]=_mm256_load_si256((__m256i*)fastptr+0x06);
			block[0x7]=_mm256_load_si256((__m256i*)fastptr+0x07);
			block[0x8]=_mm256_load_si256((__m256i*)fastptr+0x08);
			block[0x9]=_mm256_load_si256((__m256i*)fastptr+0x09);
			block[0xA]=_mm256_load_si256((__m256i*)fastptr+0x0A);
			block[0xB]=_mm256_load_si256((__m256i*)fastptr+0x0B);
			block[0xC]=_mm256_load_si256((__m256i*)fastptr+0x0C);
			block[0xD]=_mm256_load_si256((__m256i*)fastptr+0x0D);
			block[0xE]=_mm256_load_si256((__m256i*)fastptr+0x0E);
			block[0xF]=_mm256_load_si256((__m256i*)fastptr+0x0F);
			transpose16x2(block);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x10]+0, (__m128i*)slowptrs[0x00]+0, block[0x0]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x11]+0, (__m128i*)slowptrs[0x01]+0, block[0x1]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x12]+0, (__m128i*)slowptrs[0x02]+0, block[0x2]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x13]+0, (__m128i*)slowptrs[0x03]+0, block[0x3]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x14]+0, (__m128i*)slowptrs[0x04]+0, block[0x4]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x15]+0, (__m128i*)slowptrs[0x05]+0, block[0x5]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x16]+0, (__m128i*)slowptrs[0x06]+0, block[0x6]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x17]+0, (__m128i*)slowptrs[0x07]+0, block[0x7]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x18]+0, (__m128i*)slowptrs[0x08]+0, block[0x8]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x19]+0, (__m128i*)slowptrs[0x09]+0, block[0x9]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x1A]+0, (__m128i*)slowptrs[0x0A]+0, block[0xA]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x1B]+0, (__m128i*)slowptrs[0x0B]+0, block[0xB]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x1C]+0, (__m128i*)slowptrs[0x0C]+0, block[0xC]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x1D]+0, (__m128i*)slowptrs[0x0D]+0, block[0xD]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x1E]+0, (__m128i*)slowptrs[0x0E]+0, block[0xE]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x1F]+0, (__m128i*)slowptrs[0x0F]+0, block[0xF]);
			block[0x0]=_mm256_load_si256((__m256i*)fastptr+0x10);
			block[0x1]=_mm256_load_si256((__m256i*)fastptr+0x11);
			block[0x2]=_mm256_load_si256((__m256i*)fastptr+0x12);
			block[0x3]=_mm256_load_si256((__m256i*)fastptr+0x13);
			block[0x4]=_mm256_load_si256((__m256i*)fastptr+0x14);
			block[0x5]=_mm256_load_si256((__m256i*)fastptr+0x15);
			block[0x6]=_mm256_load_si256((__m256i*)fastptr+0x16);
			block[0x7]=_mm256_load_si256((__m256i*)fastptr+0x17);
			block[0x8]=_mm256_load_si256((__m256i*)fastptr+0x18);
			block[0x9]=_mm256_load_si256((__m256i*)fastptr+0x19);
			block[0xA]=_mm256_load_si256((__m256i*)fastptr+0x1A);
			block[0xB]=_mm256_load_si256((__m256i*)fastptr+0x1B);
			block[0xC]=_mm256_load_si256((__m256i*)fastptr+0x1C);
			block[0xD]=_mm256_load_si256((__m256i*)fastptr+0x1D);
			block[0xE]=_mm256_load_si256((__m256i*)fastptr+0x1E);
			block[0xF]=_mm256_load_si256((__m256i*)fastptr+0x1F);
			transpose16x2(block);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x10]+1, (__m128i*)slowptrs[0x00]+1, block[0x0]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x11]+1, (__m128i*)slowptrs[0x01]+1, block[0x1]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x12]+1, (__m128i*)slowptrs[0x02]+1, block[0x2]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x13]+1, (__m128i*)slowptrs[0x03]+1, block[0x3]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x14]+1, (__m128i*)slowptrs[0x04]+1, block[0x4]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x15]+1, (__m128i*)slowptrs[0x05]+1, block[0x5]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x16]+1, (__m128i*)slowptrs[0x06]+1, block[0x6]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x17]+1, (__m128i*)slowptrs[0x07]+1, block[0x7]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x18]+1, (__m128i*)slowptrs[0x08]+1, block[0x8]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x19]+1, (__m128i*)slowptrs[0x09]+1, block[0x9]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x1A]+1, (__m128i*)slowptrs[0x0A]+1, block[0xA]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x1B]+1, (__m128i*)slowptrs[0x0B]+1, block[0xB]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x1C]+1, (__m128i*)slowptrs[0x0C]+1, block[0xC]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x1D]+1, (__m128i*)slowptrs[0x0D]+1, block[0xD]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x1E]+1, (__m128i*)slowptrs[0x0E]+1, block[0xE]);
			_mm256_storeu2_m128i((__m128i*)slowptrs[0x1F]+1, (__m128i*)slowptrs[0x0F]+1, block[0xF]);
			fastptr+=sizeof(__m256i[NCODERS]);
			block[0]=_mm256_load_si256((__m256i*)slowptrs+0);
			block[1]=_mm256_load_si256((__m256i*)slowptrs+1);
			block[2]=_mm256_load_si256((__m256i*)slowptrs+2);
			block[3]=_mm256_load_si256((__m256i*)slowptrs+3);
			block[4]=_mm256_load_si256((__m256i*)slowptrs+4);
			block[5]=_mm256_load_si256((__m256i*)slowptrs+5);
			block[6]=_mm256_load_si256((__m256i*)slowptrs+6);
			block[7]=_mm256_load_si256((__m256i*)slowptrs+7);
			block[0]=_mm256_add_epi64(block[0], slowinc);
			block[1]=_mm256_add_epi64(block[1], slowinc);
			block[2]=_mm256_add_epi64(block[2], slowinc);
			block[3]=_mm256_add_epi64(block[3], slowinc);
			block[4]=_mm256_add_epi64(block[4], slowinc);
			block[5]=_mm256_add_epi64(block[5], slowinc);
			block[6]=_mm256_add_epi64(block[6], slowinc);
			block[7]=_mm256_add_epi64(block[7], slowinc);
			_mm256_store_si256((__m256i*)slowptrs+0, block[0]);
			_mm256_store_si256((__m256i*)slowptrs+1, block[1]);
			_mm256_store_si256((__m256i*)slowptrs+2, block[2]);
			_mm256_store_si256((__m256i*)slowptrs+3, block[3]);
			_mm256_store_si256((__m256i*)slowptrs+4, block[4]);
			_mm256_store_si256((__m256i*)slowptrs+5, block[5]);
			_mm256_store_si256((__m256i*)slowptrs+6, block[6]);
			_mm256_store_si256((__m256i*)slowptrs+7, block[7]);
		}
#endif
#if 1
		for(;kx<blockxbytes;kx+=NCODERS)
		{
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
			*slowptrs[0x10]++=fastptr[0x10];
			*slowptrs[0x11]++=fastptr[0x11];
			*slowptrs[0x12]++=fastptr[0x12];
			*slowptrs[0x13]++=fastptr[0x13];
			*slowptrs[0x14]++=fastptr[0x14];
			*slowptrs[0x15]++=fastptr[0x15];
			*slowptrs[0x16]++=fastptr[0x16];
			*slowptrs[0x17]++=fastptr[0x17];
			*slowptrs[0x18]++=fastptr[0x18];
			*slowptrs[0x19]++=fastptr[0x19];
			*slowptrs[0x1A]++=fastptr[0x1A];
			*slowptrs[0x1B]++=fastptr[0x1B];
			*slowptrs[0x1C]++=fastptr[0x1C];
			*slowptrs[0x1D]++=fastptr[0x1D];
			*slowptrs[0x1E]++=fastptr[0x1E];
			*slowptrs[0x1F]++=fastptr[0x1F];
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
#if 0
	//X  where is the PPM header?
	HANDLE fdst=CreateFileA(fn, GENERIC_WRITE, 0, 0, CREATE_ALWAYS, FILE_FLAG_NO_BUFFERING|FILE_FLAG_WRITE_THROUGH|FILE_FLAG_SEQUENTIAL_SCAN, 0);
//	HANDLE fdst=CreateFileA(fn, GENERIC_WRITE, 0, 0, CREATE_ALWAYS, FILE_FLAG_SEQUENTIAL_SCAN, 0);
	if(fdst==INVALID_HANDLE_VALUE)
	{
		LOG_ERROR("Cannot open \"%s\" for writing  GetLastError %d", fn, (int)GetLastError());
		return;
	}
	DWORD nwritten=0;
	int success=WriteFile(fdst, image, 3*iw*ih, &nwritten, 0);
	if(!success)
	{
		LOG_ERROR("Cannot save \"%s\"  GetLastError %d", fn, (int)GetLastError());
		return;
	}
	success=FlushFileBuffers(fdst);
	if(!success)
	{
		LOG_ERROR("FlushFileBuffers  GetLastError %d", fn, (int)GetLastError());
		return;
	}
	success=CloseHandle(fdst);
	if(!success)
	{
		LOG_ERROR("CloseHandle  GetLastError %d", fn, (int)GetLastError());
		return;
	}
#else
	FILE *fdst=fopen(fn, "wb");
	if(!fdst)
	{
		LOG_ERROR("Cannot open \"%s\" for writing", fn);
		return;
	}
	fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
	fwrite(image, 1, (ptrdiff_t)3*iw*ih, fdst);
	fclose(fdst);
#endif
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
int c29_codec(int argc, char **argv)
{
	if(argc!=3&&argc!=4)
	{
		printf(
			"Usage: \"%s\"  input  output  [N]    Encode/decode.\n"
			"N  =  1 Force CG / 2 Force WG4 / 3 Lossy | 4 Profile\n"
			, argv[0]
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2], nthreads0=argc<4?0:atoi(argv[3]);
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
	int wgsize=0;
	short *wgerrors=0;
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
		if(!fwd&&tag!=('2'|'9'<<8))
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
	unsigned long long ctxmask=0;//3*NCTX+3 = 54 flags	0: rare context (bypass)  1: emit stats
	const int hsize=(int)sizeof(int[3*NCTX<<8]);//3 channels
	int *hists=fwd?(int*)malloc(hsize):0;//fwd-only
	const int rhsize=(int)sizeof(int[3*256]);
	int *rhist=fwd?(int*)malloc(rhsize):0;

	int CDF2syms_size=(int)sizeof(int[3*NCTX<<PROBBITS]);
	if(fwd)//DIV-free rANS encoder reuses these as SIMD symbol info
		CDF2syms_size=(int)sizeof(rANS_SIMD_SymInfo[3*NCTX<<8]);
	unsigned *CDF2syms=(unsigned*)_mm_malloc(CDF2syms_size, sizeof(__m128i*));

	int rCDF2syms_size=(int)sizeof(int[3<<PROBBITS]);
	if(fwd)
		rCDF2syms_size=(int)sizeof(rANS_SIMD_SymInfo[3<<8]);
	unsigned *rCDF2syms=(unsigned*)_mm_malloc(rCDF2syms_size, sizeof(__m128i*));

	int ans_permute_size=sizeof(__m256i[256]);
	int *ans_permute=(int*)_mm_malloc(ans_permute_size, sizeof(__m256i));

	psize=(int)sizeof(short[4*6*NCODERS])*(blockw+16);//4 padded rows  *  {Y*NCODERS, U*NCODERS, V*NCODERS,  eY*NCODERS, eU*NCODERS, eV*NCODERS} = 2*3*32 = 192 channels  ~48*iw bytes
	pixels=(short*)_mm_malloc(psize, sizeof(__m256i));//~188 KB for 4K/12MP
	wgsize=(int)sizeof(short[2*NCODERS*3*WG_NPREDS])*(blockw+16);//2 padded rows  *  {WGY*NCODERS, WGU*NCODERS, WGV*NCODERS} * NPREDS = 3*32*8 = 768 channels  ~96*iw bytes
	wgerrors=(short*)_mm_malloc(wgsize, sizeof(__m256i));//~375 KB for 4K/12MP		NNEerrors = currerrors
	const int wgstatesize=(int)sizeof(short[2*NCODERS*3*WG_NPREDS]);//{preds, Wprederrors} * NPREDS * 3 channels * NCODERS
	short *wgstate=(short*)_mm_malloc(wgstatesize, sizeof(__m256i));
	if((fwd&&(!hists||!rhist))||!CDF2syms||!rCDF2syms||!ans_permute||!pixels||!wgerrors||!wgstate)
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
#if 0
		{
			__m256i debugblock1[16]={0}, debugblock2[16]={0};
			unsigned char *ptr=(unsigned char*)debugblock1;
			for(int kr=0;kr<2;++kr)
			{
				for(int ky=0;ky<16;++ky)
				{
					for(int kx=0;kx<16;++kx)//Y R X
						ptr[ky<<5|kr<<4|kx]=(ky<<4|kx)*kr;
					//	ptr[ky<<5|kr<<4|kx]=(kr?15-ky:ky)<<4|(kr?15-kx:kx);
				}
			}
			//for(int k=0;k<512;++k)
			//	ptr[k]=(unsigned char)k;
			memcpy(debugblock2, debugblock1, sizeof(debugblock2));
			transpose16x2(debugblock2);
			transpose16x2(debugblock2);
			if(memcmp(debugblock2, debugblock1, sizeof(debugblock1)))
				LOG_ERROR("ERROR");
			else
				LOG_ERROR("SUCCESS");
		}
#endif
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
#ifdef DISABLE_ANALYSIS
		bestrct=RCT__040_04X_X13;
		use_wg4=1;
		if(nthreads0&1)//force CG
			use_wg4=0;
		else if(nthreads0&2)//force WG4
			use_wg4=1;
#else
		prof_checkpoint(usize, "interleave");
		{
#ifdef CHECK_FLAT
			long long flatctr[3]={0};
			__m256i mflat[3];
			memset(mflat, 0, sizeof(mflat));
#endif
			ALIGN(32) long long counters[OCH_COUNT]={0};
			__m256i mcounters[OCH_COUNT];//64-bit
			__m128i half8=_mm_set1_epi8(-128);
			__m256i wordmask=_mm256_set1_epi64x(0xFFFF);
			memset(mcounters, 0, sizeof(mcounters));
			imptr=interleaved+isize;
			for(int ky=0;ky<blockh;++ky)//analysis
			{
				__m256i prev[OCH_COUNT*2];//16-bit
				memset(prev, 0, sizeof(prev));
				for(int kx=0;kx<ixbytes-3*NCODERS;kx+=3*NCODERS)
				{
					__m256i r0=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+0), half8));
					__m256i r1=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+1), half8));
					__m256i g0=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+2), half8));
					__m256i g1=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+3), half8));
					__m256i b0=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+4), half8));
					__m256i b1=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_load_si128((__m128i*)imptr+5), half8));
					imptr+=3*NCODERS;
					r0=_mm256_slli_epi16(r0, 2);
					r1=_mm256_slli_epi16(r1, 2);
					g0=_mm256_slli_epi16(g0, 2);
					g1=_mm256_slli_epi16(g1, 2);
					b0=_mm256_slli_epi16(b0, 2);
					b1=_mm256_slli_epi16(b1, 2);
#ifdef CHECK_FLAT
					__m256i cmp0=_mm256_abs_epi16(_mm256_cmpeq_epi16(r0, prev[OCH_Y400*2+0]));
					__m256i cmp1=_mm256_abs_epi16(_mm256_cmpeq_epi16(g0, prev[OCH_Y040*2+0]));
					__m256i cmp2=_mm256_abs_epi16(_mm256_cmpeq_epi16(b0, prev[OCH_Y004*2+0]));
					cmp0=_mm256_add_epi16(cmp0, _mm256_abs_epi16(_mm256_cmpeq_epi16(r1, prev[OCH_Y400*2+1])));
					cmp1=_mm256_add_epi16(cmp1, _mm256_abs_epi16(_mm256_cmpeq_epi16(g1, prev[OCH_Y040*2+1])));
					cmp2=_mm256_add_epi16(cmp2, _mm256_abs_epi16(_mm256_cmpeq_epi16(b1, prev[OCH_Y004*2+1])));
					cmp0=_mm256_add_epi16(cmp0, _mm256_srli_epi64(cmp0, 32));
					cmp1=_mm256_add_epi16(cmp1, _mm256_srli_epi64(cmp1, 32));
					cmp2=_mm256_add_epi16(cmp2, _mm256_srli_epi64(cmp2, 32));
					cmp0=_mm256_add_epi16(cmp0, _mm256_srli_epi64(cmp0, 16));
					cmp1=_mm256_add_epi16(cmp1, _mm256_srli_epi64(cmp1, 16));
					cmp2=_mm256_add_epi16(cmp2, _mm256_srli_epi64(cmp2, 16));
					cmp0=_mm256_and_si256(cmp0, wordmask);
					cmp1=_mm256_and_si256(cmp1, wordmask);
					cmp2=_mm256_and_si256(cmp2, wordmask);
					mflat[0]=_mm256_add_epi16(mflat[0], cmp0);
					mflat[1]=_mm256_add_epi16(mflat[1], cmp1);
					mflat[2]=_mm256_add_epi16(mflat[2], cmp2);
#endif
					__m256i rg0=_mm256_sub_epi16(r0, g0);
					__m256i rg1=_mm256_sub_epi16(r1, g1);
					__m256i gb0=_mm256_sub_epi16(g0, b0);
					__m256i gb1=_mm256_sub_epi16(g1, b1);
					__m256i br0=_mm256_sub_epi16(b0, r0);
					__m256i br1=_mm256_sub_epi16(b1, r1);
					__m256i t0, t1, t2, t3, t4, t5;
#define UPDATE(IDXA, IDXB, IDXC, A0, A1, B0, B1, C0, C1)\
	do\
	{\
		__m256i ta0=_mm256_sub_epi16(A0, prev[IDXA*2+0]);\
		__m256i ta1=_mm256_sub_epi16(A1, prev[IDXA*2+1]);\
		__m256i tb0=_mm256_sub_epi16(B0, prev[IDXB*2+0]);\
		__m256i tb1=_mm256_sub_epi16(B1, prev[IDXB*2+1]);\
		__m256i tc0=_mm256_sub_epi16(C0, prev[IDXC*2+0]);\
		__m256i tc1=_mm256_sub_epi16(C1, prev[IDXC*2+1]);\
		prev[IDXA*2+0]=A0;\
		prev[IDXA*2+1]=A1;\
		prev[IDXB*2+0]=B0;\
		prev[IDXB*2+1]=B1;\
		prev[IDXC*2+0]=C0;\
		prev[IDXC*2+1]=C1;\
		ta0=_mm256_abs_epi16(ta0);\
		ta1=_mm256_abs_epi16(ta1);\
		tb0=_mm256_abs_epi16(tb0);\
		tb1=_mm256_abs_epi16(tb1);\
		tc0=_mm256_abs_epi16(tc0);\
		tc1=_mm256_abs_epi16(tc1);\
		ta0=_mm256_add_epi16(ta0, ta1);\
		tb0=_mm256_add_epi16(tb0, tb1);\
		tc0=_mm256_add_epi16(tc0, tc1);\
		ta0=_mm256_add_epi16(ta0, _mm256_srli_epi64(ta0, 32));\
		tb0=_mm256_add_epi16(tb0, _mm256_srli_epi64(tb0, 32));\
		tc0=_mm256_add_epi16(tc0, _mm256_srli_epi64(tc0, 32));\
		ta0=_mm256_add_epi16(ta0, _mm256_srli_epi64(ta0, 16));\
		tb0=_mm256_add_epi16(tb0, _mm256_srli_epi64(tb0, 16));\
		tc0=_mm256_add_epi16(tc0, _mm256_srli_epi64(tc0, 16));\
		mcounters[IDXA]=_mm256_add_epi64(mcounters[IDXA], _mm256_and_si256(ta0, wordmask));\
		mcounters[IDXB]=_mm256_add_epi64(mcounters[IDXB], _mm256_and_si256(tb0, wordmask));\
		mcounters[IDXC]=_mm256_add_epi64(mcounters[IDXC], _mm256_and_si256(tc0, wordmask));\
	}while(0)
					UPDATE(OCH_Y400, OCH_Y040, OCH_Y004, r0, r1, g0, g1, b0, b1);
					UPDATE(OCH_CX40, OCH_C0X4, OCH_C40X, rg0, rg1, gb0, gb1, br0, br1);
					t0=_mm256_add_epi16(rg0, _mm256_srai_epi16(gb0, 2));//r-(3*g+b)/4 = r-g-(b-g)/4
					t1=_mm256_add_epi16(rg1, _mm256_srai_epi16(gb1, 2));
					t2=_mm256_add_epi16(rg0, _mm256_srai_epi16(br0, 2));//g-(3*r+b)/4 = g-r-(b-r)/4
					t3=_mm256_add_epi16(rg1, _mm256_srai_epi16(br1, 2));
					t4=_mm256_add_epi16(br0, _mm256_srai_epi16(rg0, 2));//b-(3*r+g)/4 = b-r-(g-r)/4
					t5=_mm256_add_epi16(br1, _mm256_srai_epi16(rg1, 2));
					UPDATE(OCH_CX31, OCH_C3X1, OCH_C31X, t0, t1, t2, t3, t4, t5);
					t0=_mm256_add_epi16(br0, _mm256_srai_epi16(gb0, 2));//r-(g+3*b)/4 = r-b-(g-b)/4
					t1=_mm256_add_epi16(br1, _mm256_srai_epi16(gb1, 2));
					t2=_mm256_add_epi16(gb0, _mm256_srai_epi16(br0, 2));//g-(r+3*b)/4 = g-b-(r-b)/4
					t3=_mm256_add_epi16(gb1, _mm256_srai_epi16(br1, 2));
					t4=_mm256_add_epi16(gb0, _mm256_srai_epi16(rg0, 2));//b-(r+3*g)/4 = b-g-(r-g)/4
					t5=_mm256_add_epi16(gb1, _mm256_srai_epi16(rg1, 2));
					UPDATE(OCH_CX13, OCH_C1X3, OCH_C13X, t0, t1, t2, t3, t4, t5);
					t0=_mm256_srai_epi16(_mm256_sub_epi16(rg0, br0), 1);//r-(g+b)/2 = (r-g + r-b)/2
					t1=_mm256_srai_epi16(_mm256_sub_epi16(rg1, br1), 1);
					t2=_mm256_srai_epi16(_mm256_sub_epi16(gb0, rg0), 1);//g-(r+b)/2 = (g-r + g-b)/2
					t3=_mm256_srai_epi16(_mm256_sub_epi16(gb1, rg1), 1);
					t4=_mm256_srai_epi16(_mm256_sub_epi16(br0, gb0), 1);//b-(r+g)/2 = (b-r + b-g)/2
					t5=_mm256_srai_epi16(_mm256_sub_epi16(br1, gb1), 1);
					UPDATE(OCH_CX22, OCH_C2X2, OCH_C22X, t0, t1, t2, t3, t4, t5);
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
				//printf("%-14s %12lld + %12lld + %12lld = %12lld\n", rct_names[kt], counters[rct[0]], counters[rct[1]], counters[rct[2]], currerr);
#endif
				if(!kt||minerr>currerr)
				{
					minerr=currerr;
					bestrct=kt;
				}
			}
#ifdef DISABLE_WG
			use_wg4=0;
#else
			const unsigned char *rct=rct_combinations[bestrct];
#ifdef CHECK_FLAT
			{
				ALIGN(32) long long temp[4]={0};
				_mm256_store_si256((__m256i*)temp, mflat[0]);
				flatctr[0]=temp[0]+temp[1]+temp[2]+temp[3];
				_mm256_store_si256((__m256i*)temp, mflat[1]);
				flatctr[1]=temp[0]+temp[1]+temp[2]+temp[3];
				_mm256_store_si256((__m256i*)temp, mflat[2]);
				flatctr[2]=temp[0]+temp[1]+temp[2]+temp[3];

			}
			use_wg4=flatctr[rct[0]]<isize/6;
#else
			use_wg4=counters[rct[0]]+counters[rct[1]]+counters[rct[2]] > isize*2;//FIXME tune
#endif
			if(nthreads0&1)//force CG
				use_wg4=0;
			else if(nthreads0&2)//force WG4
				use_wg4=1;
#endif
		}
		prof_checkpoint(usize, "analysis");
#ifdef LOUD
		printf("%s  WG4=%d  %td bytes\n", rct_names[bestrct], use_wg4, usize);
#endif
#endif

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
		use_wg4=flags&1;
		bestrct=flags>>1;

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
	memset(wgerrors, 0, wgsize);
	__m256i mctxmax=_mm256_set1_epi16(NCTX-1);
	__m256i mctxuoffset=_mm256_set1_epi16(NCTX);
	__m256i mctxvoffset=_mm256_set1_epi16(NCTX*2);
	__m256i amin=_mm256_set1_epi16(-128);
#ifndef ENABLE_MA
	__m256i amax=_mm256_set1_epi16(127);
#endif
	__m128i half8=_mm_set1_epi8(-128);
	__m256i bytemask=_mm256_set1_epi16(255);
	__m256i wordmask=_mm256_set1_epi32(0xFFFF);
	__m256i myuv[6];
	memset(myuv, 0, sizeof(myuv));
	unsigned char *ctxptr=interleaved;
	imptr=interleaved+(fwd?isize:0);
	__m256i mstate[4];
	short *wgWerrors=wgstate+0*(ptrdiff_t)NCODERS*3*WG_NPREDS;
	__m256i *wgpreds=(__m256i*)(wgstate+1*(ptrdiff_t)NCODERS*3*WG_NPREDS);
	memset(wgstate, 0, wgstatesize);
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
#ifdef SERIAL_MAIN
		int ufromy=-(combination[II_COEFF_U_SUB_Y]!=0);
		int vcoeff0=combination[II_COEFF_V_SUB_Y];
		int vcoeff1=combination[II_COEFF_V_SUB_U];
		//short prev[2*NCODERS]={0};
		int yuv[3]={0}, uoffset=0, voffset=0, errors[3]={0}, syms[3]={0};
		for(int kx=0;kx<ixbytes;kx+=3*NCODERS)
		{
			short
				*NW	=rows[1]-1*6*NCODERS,
				*N	=rows[1]+0*6*NCODERS,
				*NEE	=rows[1]+2*6*NCODERS,
				*NEEE	=rows[1]+3*6*NCODERS,
				*W	=rows[0]-1*6*NCODERS,
				*curr	=rows[0]+0*6*NCODERS;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
			for(int k=0;k<NCODERS;++k)
			{
				int ctx0=FLOOR_LOG2(W[k+0*NCODERS]*W[k+0*NCODERS]+1);
				int ctx1=FLOOR_LOG2(W[k+1*NCODERS]*W[k+1*NCODERS]+1);
				int ctx2=FLOOR_LOG2(W[k+2*NCODERS]*W[k+2*NCODERS]+1);
				if(ctx0>NCTX-1)ctx0=NCTX-1;
				if(ctx1>NCTX-1)ctx1=NCTX-1;
				if(ctx2>NCTX-1)ctx2=NCTX-1;
				ctx1+=NCTX;
				ctx2+=NCTX*2;

				int vmax0=W[k+0*NCODERS], vmin0=W[k+0*NCODERS];
				int vmax1=W[k+1*NCODERS], vmin1=W[k+1*NCODERS];
				int vmax2=W[k+2*NCODERS], vmin2=W[k+2*NCODERS];
				if(vmax0<vmin0)vmin0=N[k+0*NCODERS], vmax0=W[k+0*NCODERS];
				if(vmax1<vmin1)vmin1=N[k+1*NCODERS], vmax1=W[k+1*NCODERS];
				if(vmax2<vmin2)vmin2=N[k+2*NCODERS], vmax2=W[k+2*NCODERS];
				int pred0=vmin0+vmax0-NW[k+0*NCODERS];
				int pred1=vmin1+vmax1-NW[k+1*NCODERS];
				int pred2=vmin2+vmax2-NW[k+2*NCODERS];
				CLAMP2(pred0, vmin0, vmax0);
				CLAMP2(pred1, vmin1, vmax1);
				CLAMP2(pred2, vmin2, vmax2);
				if(fwd)
				{
					yuv[0]=imptr[k+yidx]-128;
					yuv[1]=imptr[k+uidx]-128;
					yuv[2]=imptr[k+vidx]-128;
					uoffset=yuv[0]&ufromy;
					voffset=vcoeff0*yuv[0]+vcoeff1*yuv[1];
					pred1+=uoffset;
					pred2+=voffset;
					pred2>>=2;
					CLAMP2(pred1, -128, 127);
					CLAMP2(pred2, -128, 127);
					errors[0]=(unsigned char)(yuv[0]-pred0+128);
					errors[1]=(unsigned char)(yuv[1]-pred1+128);
					errors[2]=(unsigned char)(yuv[2]-pred2+128);
					((unsigned short*)ctxptr)[k+0*NCODERS]=syms[0]=ctx0<<8|errors[0];
					((unsigned short*)ctxptr)[k+1*NCODERS]=syms[1]=ctx1<<8|errors[1];
					((unsigned short*)ctxptr)[k+2*NCODERS]=syms[2]=ctx2<<8|errors[2];
					++hists[syms[0]];
					++hists[syms[1]];
					++hists[syms[2]];
				}
				else
				{
					LOG_ERROR("");
				}
				curr[k+(0+0)*NCODERS]=yuv[0];
				curr[k+(0+1)*NCODERS]=yuv[1]-uoffset;
				curr[k+(0+2)*NCODERS]=(yuv[2]<<2)-voffset;
				curr[k+(3+0)*NCODERS]=(2*W[k+(3+0)*NCODERS]+((errors[0]<<1^errors[0]>>31)<<GRBITS)+MAXVAR(NEE[k+(3+0)*NCODERS], NEEE[k+(3+0)*NCODERS]))>>2;
				curr[k+(3+1)*NCODERS]=(2*W[k+(3+1)*NCODERS]+((errors[1]<<1^errors[1]>>31)<<GRBITS)+MAXVAR(NEE[k+(3+1)*NCODERS], NEEE[k+(3+1)*NCODERS]))>>2;
				curr[k+(3+2)*NCODERS]=(2*W[k+(3+2)*NCODERS]+((errors[2]<<1^errors[2]>>31)<<GRBITS)+MAXVAR(NEE[k+(3+2)*NCODERS], NEEE[k+(3+2)*NCODERS]))>>2;
			}
			if(fwd)
				ctxptr+=sizeof(short[3][NCODERS]);
			rows[0]+=6*NCODERS;
			rows[1]+=6*NCODERS;
			rows[2]+=6*NCODERS;
			rows[3]+=6*NCODERS;
			imptr+=3*NCODERS;
		}
		(void)mctxuoffset;
		(void)mctxvoffset;
		(void)mctxmax;
		(void)uhelpmask;
		(void)vc0;
		(void)vc1;
		(void)amin;
		(void)amax;
		(void)bytemask;
		(void)half8;
		(void)wordmask;

		(void)erows;
		(void)wgpreds;
		(void)wgWerrors;
#else
		ALIGN(32) unsigned short syms[3*NCODERS]={0};
		__m256i NW[6], N[6], W[6];
		__m256i eW[6], ecurr[6], eNEE[6], eNEEE[6];
		memset(NW, 0, sizeof(NW));
		memset(N, 0, sizeof(N));
		memset(W, 0, sizeof(W));
		memset(eW, 0, sizeof(eW));
		memset(ecurr, 0, sizeof(ecurr));
		memset(eNEE, 0, sizeof(eNEE));
		memset(eNEEE, 0, sizeof(eNEEE));
		memset(wgpreds, 0, sizeof(short[NCODERS*3*WG_NPREDS]));
		//                       (__m256i*)(rows[-Y]+(E+C+X*6)*NCODERS)+R
		eNEE[0]=_mm256_load_si256((__m256i*)(rows[1]+(3+0+2*6)*NCODERS)+0);
		eNEE[1]=_mm256_load_si256((__m256i*)(rows[1]+(3+0+2*6)*NCODERS)+1);
		eNEE[2]=_mm256_load_si256((__m256i*)(rows[1]+(3+1+2*6)*NCODERS)+0);
		eNEE[3]=_mm256_load_si256((__m256i*)(rows[1]+(3+1+2*6)*NCODERS)+1);
		eNEE[4]=_mm256_load_si256((__m256i*)(rows[1]+(3+2+2*6)*NCODERS)+0);
		eNEE[5]=_mm256_load_si256((__m256i*)(rows[1]+(3+2+2*6)*NCODERS)+1);
		for(int kx=0;kx<ixbytes;kx+=3*NCODERS)
		{
			//						(E+C+X*6)*NCODERS
			N	[0]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+0*6)*NCODERS)+0);//y0 neighbors
			N	[1]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+0*6)*NCODERS)+1);//y1
			N	[2]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+0*6)*NCODERS)+0);//u0
			N	[3]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+0*6)*NCODERS)+1);//u1
			N	[4]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+0*6)*NCODERS)+0);//v0
			N	[5]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+0*6)*NCODERS)+1);//v1
			__m256i
				predY0, predY1, ctxY0, ctxY1,
				predU0, predU1, ctxU0, ctxU1,
				predV0, predV1, ctxV0, ctxV1;
			{
				//context = FLOOR_LOG2(eW*eW+1)
				__m256i one=_mm256_set1_epi32(1);
				__m256i cy0=_mm256_and_si256(eW[0], wordmask), cy1=_mm256_srli_epi32(eW[0], 16);
				__m256i cy2=_mm256_and_si256(eW[1], wordmask), cy3=_mm256_srli_epi32(eW[1], 16);
				__m256i cu0=_mm256_and_si256(eW[2], wordmask), cu1=_mm256_srli_epi32(eW[2], 16);
				__m256i cu2=_mm256_and_si256(eW[3], wordmask), cu3=_mm256_srli_epi32(eW[3], 16);
				__m256i cv0=_mm256_and_si256(eW[4], wordmask), cv1=_mm256_srli_epi32(eW[4], 16);
				__m256i cv2=_mm256_and_si256(eW[5], wordmask), cv3=_mm256_srli_epi32(eW[5], 16);
				cy0=_mm256_mullo_epi32(cy0, cy0);
				cy1=_mm256_mullo_epi32(cy1, cy1);
				cy2=_mm256_mullo_epi32(cy2, cy2);
				cy3=_mm256_mullo_epi32(cy3, cy3);
				cu0=_mm256_mullo_epi32(cu0, cu0);
				cu1=_mm256_mullo_epi32(cu1, cu1);
				cu2=_mm256_mullo_epi32(cu2, cu2);
				cu3=_mm256_mullo_epi32(cu3, cu3);
				cv0=_mm256_mullo_epi32(cv0, cv0);
				cv1=_mm256_mullo_epi32(cv1, cv1);
				cv2=_mm256_mullo_epi32(cv2, cv2);
				cv3=_mm256_mullo_epi32(cv3, cv3);
				cy0=_mm256_add_epi32(cy0, one);
				cy1=_mm256_add_epi32(cy1, one);
				cy2=_mm256_add_epi32(cy2, one);
				cy3=_mm256_add_epi32(cy3, one);
				cu0=_mm256_add_epi32(cu0, one);
				cu1=_mm256_add_epi32(cu1, one);
				cu2=_mm256_add_epi32(cu2, one);
				cu3=_mm256_add_epi32(cu3, one);
				cv0=_mm256_add_epi32(cv0, one);
				cv1=_mm256_add_epi32(cv1, one);
				cv2=_mm256_add_epi32(cv2, one);
				cv3=_mm256_add_epi32(cv3, one);
				//FLOOR_LOG2_32x8(X) = _mm256_sub_epi32(_mm256_srli_epi32(_mm256_castps_si256(_mm256_cvtepi32_ps(X)), 23), _mm256_set1_epi32(127))
				cy0=_mm256_castps_si256(_mm256_cvtepi32_ps(cy0));
				cy1=_mm256_castps_si256(_mm256_cvtepi32_ps(cy1));
				cy2=_mm256_castps_si256(_mm256_cvtepi32_ps(cy2));
				cy3=_mm256_castps_si256(_mm256_cvtepi32_ps(cy3));
				cu0=_mm256_castps_si256(_mm256_cvtepi32_ps(cu0));
				cu1=_mm256_castps_si256(_mm256_cvtepi32_ps(cu1));
				cu2=_mm256_castps_si256(_mm256_cvtepi32_ps(cu2));
				cu3=_mm256_castps_si256(_mm256_cvtepi32_ps(cu3));
				cv0=_mm256_castps_si256(_mm256_cvtepi32_ps(cv0));
				cv1=_mm256_castps_si256(_mm256_cvtepi32_ps(cv1));
				cv2=_mm256_castps_si256(_mm256_cvtepi32_ps(cv2));
				cv3=_mm256_castps_si256(_mm256_cvtepi32_ps(cv3));
				cy0=_mm256_srli_epi32(cy0, 23);
				cy1=_mm256_srli_epi32(cy1, 23);
				cy2=_mm256_srli_epi32(cy2, 23);
				cy3=_mm256_srli_epi32(cy3, 23);
				cu0=_mm256_srli_epi32(cu0, 23);
				cu1=_mm256_srli_epi32(cu1, 23);
				cu2=_mm256_srli_epi32(cu2, 23);
				cu3=_mm256_srli_epi32(cu3, 23);
				cv0=_mm256_srli_epi32(cv0, 23);
				cv1=_mm256_srli_epi32(cv1, 23);
				cv2=_mm256_srli_epi32(cv2, 23);
				cv3=_mm256_srli_epi32(cv3, 23);
				__m256i expbias=_mm256_set1_epi32(127);
				cy0=_mm256_sub_epi32(cy0, expbias);
				cy1=_mm256_sub_epi32(cy1, expbias);
				cy2=_mm256_sub_epi32(cy2, expbias);
				cy3=_mm256_sub_epi32(cy3, expbias);
				cu0=_mm256_sub_epi32(cu0, expbias);
				cu1=_mm256_sub_epi32(cu1, expbias);
				cu2=_mm256_sub_epi32(cu2, expbias);
				cu3=_mm256_sub_epi32(cu3, expbias);
				cv0=_mm256_sub_epi32(cv0, expbias);
				cv1=_mm256_sub_epi32(cv1, expbias);
				cv2=_mm256_sub_epi32(cv2, expbias);
				cv3=_mm256_sub_epi32(cv3, expbias);
				cy1=_mm256_slli_epi32(cy1, 16);
				cy3=_mm256_slli_epi32(cy3, 16);
				cu1=_mm256_slli_epi32(cu1, 16);
				cu3=_mm256_slli_epi32(cu3, 16);
				cv1=_mm256_slli_epi32(cv1, 16);
				cv3=_mm256_slli_epi32(cv3, 16);
				ctxY0=_mm256_or_si256(cy0, cy1);
				ctxY1=_mm256_or_si256(cy2, cy3);
				ctxU0=_mm256_or_si256(cu0, cu1);
				ctxU1=_mm256_or_si256(cu2, cu3);
				ctxV0=_mm256_or_si256(cv0, cv1);
				ctxV1=_mm256_or_si256(cv2, cv3);
				ctxY0=_mm256_min_epi16(ctxY0, mctxmax);
				ctxY1=_mm256_min_epi16(ctxY1, mctxmax);
				ctxU0=_mm256_min_epi16(ctxU0, mctxmax);
				ctxU1=_mm256_min_epi16(ctxU1, mctxmax);
				ctxV0=_mm256_min_epi16(ctxV0, mctxmax);
				ctxV1=_mm256_min_epi16(ctxV1, mctxmax);
			}
			{
				__m256i ymin0=_mm256_min_epi16(N[0], W[0]);
				__m256i ymax0=_mm256_max_epi16(N[0], W[0]);
				__m256i ymin1=_mm256_min_epi16(N[1], W[1]);
				__m256i ymax1=_mm256_max_epi16(N[1], W[1]);
				__m256i umin0=_mm256_min_epi16(N[2], W[2]);
				__m256i umax0=_mm256_max_epi16(N[2], W[2]);
				__m256i umin1=_mm256_min_epi16(N[3], W[3]);
				__m256i umax1=_mm256_max_epi16(N[3], W[3]);
				__m256i vmin0=_mm256_min_epi16(N[4], W[4]);
				__m256i vmax0=_mm256_max_epi16(N[4], W[4]);
				__m256i vmin1=_mm256_min_epi16(N[5], W[5]);
				__m256i vmax1=_mm256_max_epi16(N[5], W[5]);
				predY0=_mm256_add_epi16(N[0], W[0]);
				predY1=_mm256_add_epi16(N[1], W[1]);
				predU0=_mm256_add_epi16(N[2], W[2]);
				predU1=_mm256_add_epi16(N[3], W[3]);
				predV0=_mm256_add_epi16(N[4], W[4]);
				predV1=_mm256_add_epi16(N[5], W[5]);
				predY0=_mm256_sub_epi16(predY0, NW[0]);
				predY1=_mm256_sub_epi16(predY1, NW[1]);
				predU0=_mm256_sub_epi16(predU0, NW[2]);
				predU1=_mm256_sub_epi16(predU1, NW[3]);
				predV0=_mm256_sub_epi16(predV0, NW[4]);
				predV1=_mm256_sub_epi16(predV1, NW[5]);
				if(use_wg4)//predict
				{
					__m256i cache[6];
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
					wgpreds[0*6+0]=N[0];
					wgpreds[0*6+1]=N[1];
					wgpreds[0*6+2]=N[2];
					wgpreds[0*6+3]=N[3];
					wgpreds[0*6+4]=N[4];
					wgpreds[0*6+5]=N[5];

					//W
					wgpreds[1*6+0]=W[0];
					wgpreds[1*6+1]=W[1];
					wgpreds[1*6+2]=W[2];
					wgpreds[1*6+3]=W[3];
					wgpreds[1*6+4]=W[4];
					wgpreds[1*6+5]=W[5];

					//3*(N-NN)+NNN
					cache[0]=_mm256_sub_epi16(N[0], _mm256_load_si256((__m256i*)(rows[2]+(0+0+0*6)*NCODERS)+0));//N-NN
					cache[1]=_mm256_sub_epi16(N[1], _mm256_load_si256((__m256i*)(rows[2]+(0+0+0*6)*NCODERS)+1));
					cache[2]=_mm256_sub_epi16(N[2], _mm256_load_si256((__m256i*)(rows[2]+(0+1+0*6)*NCODERS)+0));
					cache[3]=_mm256_sub_epi16(N[3], _mm256_load_si256((__m256i*)(rows[2]+(0+1+0*6)*NCODERS)+1));
					cache[4]=_mm256_sub_epi16(N[4], _mm256_load_si256((__m256i*)(rows[2]+(0+2+0*6)*NCODERS)+0));
					cache[5]=_mm256_sub_epi16(N[5], _mm256_load_si256((__m256i*)(rows[2]+(0+2+0*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_slli_epi16(cache[3], 1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_slli_epi16(cache[4], 1));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_slli_epi16(cache[5], 1));
					wgpreds[2*6+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[3]+(0+0+0*6)*NCODERS)+0));//+NNN
					wgpreds[2*6+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[3]+(0+0+0*6)*NCODERS)+1));
					wgpreds[2*6+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[3]+(0+1+0*6)*NCODERS)+0));
					wgpreds[2*6+3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[3]+(0+1+0*6)*NCODERS)+1));
					wgpreds[2*6+4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[3]+(0+2+0*6)*NCODERS)+0));
					wgpreds[2*6+5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[3]+(0+2+0*6)*NCODERS)+1));

					//4*(W+WWW)-6*WW-WWWW		X  HORRIBLE
#if 0
					cache[0]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-1*6)*NCODERS)+0);//W
					cache[1]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-1*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-1*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-1*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-1*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-1*6)*NCODERS)+1);
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[0]+(0+0-3*6)*NCODERS)+0));//W+WWW
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[0]+(0+0-3*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[0]+(0+1-3*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[0]+(0+1-3*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[0]+(0+2-3*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[0]+(0+2-3*6)*NCODERS)+1));
					wgpreds[3*6+0]=_mm256_slli_epi16(cache[0], 2);//4*(W+WWW)
					wgpreds[3*6+1]=_mm256_slli_epi16(cache[1], 2);
					wgpreds[3*6+2]=_mm256_slli_epi16(cache[2], 2);
					wgpreds[3*6+3]=_mm256_slli_epi16(cache[3], 2);
					wgpreds[3*6+4]=_mm256_slli_epi16(cache[4], 2);
					wgpreds[3*6+5]=_mm256_slli_epi16(cache[5], 2);
					cache[0]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-2*6)*NCODERS)+0);//WW
					cache[1]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-2*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-2*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-2*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-2*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-2*6)*NCODERS)+1);
					cache[0]=_mm256_add_epi16(_mm256_slli_epi16(cache[0], 2), _mm256_slli_epi16(cache[0], 1));//6*WW
					cache[1]=_mm256_add_epi16(_mm256_slli_epi16(cache[1], 2), _mm256_slli_epi16(cache[1], 1));
					cache[2]=_mm256_add_epi16(_mm256_slli_epi16(cache[2], 2), _mm256_slli_epi16(cache[2], 1));
					cache[3]=_mm256_add_epi16(_mm256_slli_epi16(cache[3], 2), _mm256_slli_epi16(cache[3], 1));
					cache[4]=_mm256_add_epi16(_mm256_slli_epi16(cache[4], 2), _mm256_slli_epi16(cache[4], 1));
					cache[5]=_mm256_add_epi16(_mm256_slli_epi16(cache[5], 2), _mm256_slli_epi16(cache[5], 1));
					wgpreds[3*6+0]=_mm256_sub_epi16(wgpreds[3*6+0], cache[0]);//4*(W+WWW)-6*WW
					wgpreds[3*6+1]=_mm256_sub_epi16(wgpreds[3*6+1], cache[1]);
					wgpreds[3*6+2]=_mm256_sub_epi16(wgpreds[3*6+2], cache[2]);
					wgpreds[3*6+3]=_mm256_sub_epi16(wgpreds[3*6+3], cache[3]);
					wgpreds[3*6+4]=_mm256_sub_epi16(wgpreds[3*6+4], cache[4]);
					wgpreds[3*6+5]=_mm256_sub_epi16(wgpreds[3*6+5], cache[5]);
					wgpreds[3*6+0]=_mm256_sub_epi16(wgpreds[3*6+0], _mm256_load_si256((__m256i*)(rows[0]+(0+0-4*6)*NCODERS)+0));//4*(W+WWW)-6*WW-WWWW
					wgpreds[3*6+1]=_mm256_sub_epi16(wgpreds[3*6+1], _mm256_load_si256((__m256i*)(rows[0]+(0+0-4*6)*NCODERS)+1));
					wgpreds[3*6+2]=_mm256_sub_epi16(wgpreds[3*6+2], _mm256_load_si256((__m256i*)(rows[0]+(0+1-4*6)*NCODERS)+0));
					wgpreds[3*6+3]=_mm256_sub_epi16(wgpreds[3*6+3], _mm256_load_si256((__m256i*)(rows[0]+(0+1-4*6)*NCODERS)+1));
					wgpreds[3*6+4]=_mm256_sub_epi16(wgpreds[3*6+4], _mm256_load_si256((__m256i*)(rows[0]+(0+2-4*6)*NCODERS)+0));
					wgpreds[3*6+5]=_mm256_sub_epi16(wgpreds[3*6+5], _mm256_load_si256((__m256i*)(rows[0]+(0+2-4*6)*NCODERS)+1));
#endif
					//3*(W-WW)+WWW
#if 1
					cache[0]=_mm256_sub_epi16(W[0], _mm256_load_si256((__m256i*)(rows[0]+(0+0-2*6)*NCODERS)+0));//W-WW
					cache[1]=_mm256_sub_epi16(W[1], _mm256_load_si256((__m256i*)(rows[0]+(0+0-2*6)*NCODERS)+1));
					cache[2]=_mm256_sub_epi16(W[2], _mm256_load_si256((__m256i*)(rows[0]+(0+1-2*6)*NCODERS)+0));
					cache[3]=_mm256_sub_epi16(W[3], _mm256_load_si256((__m256i*)(rows[0]+(0+1-2*6)*NCODERS)+1));
					cache[4]=_mm256_sub_epi16(W[4], _mm256_load_si256((__m256i*)(rows[0]+(0+2-2*6)*NCODERS)+0));
					cache[5]=_mm256_sub_epi16(W[5], _mm256_load_si256((__m256i*)(rows[0]+(0+2-2*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_slli_epi16(cache[0], 1));//*3
					cache[1]=_mm256_add_epi16(cache[1], _mm256_slli_epi16(cache[1], 1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_slli_epi16(cache[2], 1));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_slli_epi16(cache[3], 1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_slli_epi16(cache[4], 1));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_slli_epi16(cache[5], 1));
					wgpreds[3*6+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[0]+(0+0-3*6)*NCODERS)+0));//+WWW
					wgpreds[3*6+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[0]+(0+0-3*6)*NCODERS)+1));
					wgpreds[3*6+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[0]+(0+1-3*6)*NCODERS)+0));
					wgpreds[3*6+3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[0]+(0+1-3*6)*NCODERS)+1));
					wgpreds[3*6+4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[0]+(0+2-3*6)*NCODERS)+0));
					wgpreds[3*6+5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[0]+(0+2-3*6)*NCODERS)+1));
#endif
					//2*W-WW	X  HORRIBLE
#if 0
					cache[0]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-2*6)*NCODERS)+0);//WW
					cache[1]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-2*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-2*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-2*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-2*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-2*6)*NCODERS)+1);
					wgpreds[3*6+0]=_mm256_sub_epi16(_mm256_slli_epi16(W[0], 1), cache[0]);
					wgpreds[3*6+1]=_mm256_sub_epi16(_mm256_slli_epi16(W[1], 1), cache[1]);
					wgpreds[3*6+2]=_mm256_sub_epi16(_mm256_slli_epi16(W[2], 1), cache[2]);
					wgpreds[3*6+3]=_mm256_sub_epi16(_mm256_slli_epi16(W[3], 1), cache[3]);
					wgpreds[3*6+4]=_mm256_sub_epi16(_mm256_slli_epi16(W[4], 1), cache[4]);
					wgpreds[3*6+5]=_mm256_sub_epi16(_mm256_slli_epi16(W[5], 1), cache[5]);
#endif

					//W+NE-N
					cache[0]=_mm256_sub_epi16(W[0], N[0]);
					cache[1]=_mm256_sub_epi16(W[1], N[1]);
					cache[2]=_mm256_sub_epi16(W[2], N[2]);
					cache[3]=_mm256_sub_epi16(W[3], N[3]);
					cache[4]=_mm256_sub_epi16(W[4], N[4]);
					cache[5]=_mm256_sub_epi16(W[5], N[5]);
					wgpreds[4*6+0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+0));//+NE
					wgpreds[4*6+1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+1));
					wgpreds[4*6+2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+0));
					wgpreds[4*6+3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+1));
					wgpreds[4*6+4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+0));
					wgpreds[4*6+5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+1));

					//(WWWW+WWW+NNN+NEE+NEEE+NEEEE-(NW+N))>>2
#if 1
					cache[0]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-4*6)*NCODERS)+0);//WWWW
					cache[1]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-4*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-4*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-4*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-4*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-4*6)*NCODERS)+1);
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[0]+(0+0-3*6)*NCODERS)+0));//+WWW
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[0]+(0+0-3*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[0]+(0+1-3*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[0]+(0+1-3*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[0]+(0+2-3*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[0]+(0+2-3*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[3]+(0+0+0*6)*NCODERS)+0));//+NNN
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[3]+(0+0+0*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[3]+(0+1+0*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[3]+(0+1+0*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[3]+(0+2+0*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[3]+(0+2+0*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[1]+(0+0+2*6)*NCODERS)+0));//+NEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[1]+(0+0+2*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[1]+(0+1+2*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[1]+(0+1+2*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[1]+(0+2+2*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[1]+(0+2+2*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[1]+(0+0+3*6)*NCODERS)+0));//+NEEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[1]+(0+0+3*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[1]+(0+1+3*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[1]+(0+1+3*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[1]+(0+2+3*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[1]+(0+2+3*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[1]+(0+0+4*6)*NCODERS)+0));//+NEEEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[1]+(0+0+4*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[1]+(0+1+4*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[1]+(0+1+4*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[1]+(0+2+4*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[1]+(0+2+4*6)*NCODERS)+1));
					cache[0]=_mm256_sub_epi16(cache[0], _mm256_add_epi16(N[0], NW[0]));
					cache[1]=_mm256_sub_epi16(cache[1], _mm256_add_epi16(N[1], NW[1]));
					cache[2]=_mm256_sub_epi16(cache[2], _mm256_add_epi16(N[2], NW[2]));
					cache[3]=_mm256_sub_epi16(cache[3], _mm256_add_epi16(N[3], NW[3]));
					cache[4]=_mm256_sub_epi16(cache[4], _mm256_add_epi16(N[4], NW[4]));
					cache[5]=_mm256_sub_epi16(cache[5], _mm256_add_epi16(N[5], NW[5]));
				//	cache[0]=_mm256_sub_epi16(cache[0], _mm256_slli_epi16(N[0], 1));
				//	cache[1]=_mm256_sub_epi16(cache[1], _mm256_slli_epi16(N[1], 1));
				//	cache[2]=_mm256_sub_epi16(cache[2], _mm256_slli_epi16(N[2], 1));
				//	cache[3]=_mm256_sub_epi16(cache[3], _mm256_slli_epi16(N[3], 1));
				//	cache[4]=_mm256_sub_epi16(cache[4], _mm256_slli_epi16(N[4], 1));
				//	cache[5]=_mm256_sub_epi16(cache[5], _mm256_slli_epi16(N[5], 1));
					wgpreds[5*6+0]=_mm256_srai_epi16(cache[0], 2);
					wgpreds[5*6+1]=_mm256_srai_epi16(cache[1], 2);
					wgpreds[5*6+2]=_mm256_srai_epi16(cache[2], 2);
					wgpreds[5*6+3]=_mm256_srai_epi16(cache[3], 2);
					wgpreds[5*6+4]=_mm256_srai_epi16(cache[4], 2);
					wgpreds[5*6+5]=_mm256_srai_epi16(cache[5], 2);
#endif
					//(WWWW+NNN+NEEE+NEEEE)>>2
#if 0
					cache[0]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-4*6)*NCODERS)+0);//WWWW
					cache[1]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-4*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-4*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-4*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-4*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-4*6)*NCODERS)+1);
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[3]+(0+0+0*6)*NCODERS)+0));//NNN	vs NNNE?
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[3]+(0+0+0*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[3]+(0+1+0*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[3]+(0+1+0*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[3]+(0+2+0*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[3]+(0+2+0*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[1]+(0+0+3*6)*NCODERS)+0));//NEEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[1]+(0+0+3*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[1]+(0+1+3*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[1]+(0+1+3*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[1]+(0+2+3*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[1]+(0+2+3*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[1]+(0+0+4*6)*NCODERS)+0));//NEEEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[1]+(0+0+4*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[1]+(0+1+4*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[1]+(0+1+4*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[1]+(0+2+4*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[1]+(0+2+4*6)*NCODERS)+1));
					wgpreds[5*6+0]=_mm256_srai_epi16(cache[0], 2);
					wgpreds[5*6+1]=_mm256_srai_epi16(cache[1], 2);
					wgpreds[5*6+2]=_mm256_srai_epi16(cache[2], 2);
					wgpreds[5*6+3]=_mm256_srai_epi16(cache[3], 2);
					wgpreds[5*6+4]=_mm256_srai_epi16(cache[4], 2);
					wgpreds[5*6+5]=_mm256_srai_epi16(cache[5], 2);
#endif
					//(WWWW+WWW+NNNN+NNNNE+NNN+NNNE+NEEE+NEEEE+N+W-2*NW)>>3		X
					//(WWWWW+WWWW+WWW+NNNN+NNNNE+NNN+NNNE+NEEE+NEEEE+W-2*NW)>>3	X
					//(WWWWW+WWWW+WWW+NNNN+NNNNE+NNN+NNNE+NEEE+NEEEE-W)>>3		X
					//(WWWW+WWW+NNNN+NNNNE+NNN+NNNE+NEEE+NEEEE-4*W)>>2		X
					//(WWWW+WWW+NNNN+NNNNE+NNN+NNNE+NEEE+NEEEE-(W+NW+N+NE))>>2	X
					//(WWWW+WWW+NNNN+NNNNE+NNN+NNNE+NEEE+NEEEE-2*(N+W))>>2		X
#if 0
					cache[0]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-4*6)*NCODERS)+0);//WWWW
					cache[1]=_mm256_load_si256((__m256i*)(rows[0]+(0+0-4*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-4*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[0]+(0+1-4*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-4*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[0]+(0+2-4*6)*NCODERS)+1);
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[0]+(0+0-3*6)*NCODERS)+0));//WWW
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[0]+(0+0-3*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[0]+(0+1-3*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[0]+(0+1-3*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[0]+(0+2-3*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[0]+(0+2-3*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[0]+(0+0+0*6)*NCODERS)+0));//NNNN
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[0]+(0+0+0*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[0]+(0+0+1*6)*NCODERS)+0));//NNNNE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[0]+(0+0+1*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[0]+(0+1+1*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[0]+(0+1+1*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[0]+(0+2+1*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[0]+(0+2+1*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[3]+(0+0+0*6)*NCODERS)+0));//NNN
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[3]+(0+0+0*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[3]+(0+1+0*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[3]+(0+1+0*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[3]+(0+2+0*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[3]+(0+2+0*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[3]+(0+0+1*6)*NCODERS)+0));//NNNE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[3]+(0+0+1*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[3]+(0+1+1*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[3]+(0+1+1*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[3]+(0+2+1*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[3]+(0+2+1*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[1]+(0+0+3*6)*NCODERS)+0));//NEEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[1]+(0+0+3*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[1]+(0+1+3*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[1]+(0+1+3*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[1]+(0+2+3*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[1]+(0+2+3*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[1]+(0+0+4*6)*NCODERS)+0));//NEEEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[1]+(0+0+4*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[1]+(0+1+4*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[1]+(0+1+4*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[1]+(0+2+4*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[1]+(0+2+4*6)*NCODERS)+1));
					cache[0]=_mm256_add_epi16(cache[0], N[0]);
					cache[1]=_mm256_add_epi16(cache[1], N[1]);
					cache[2]=_mm256_add_epi16(cache[2], N[2]);
					cache[3]=_mm256_add_epi16(cache[3], N[3]);
					cache[4]=_mm256_add_epi16(cache[4], N[4]);
					cache[5]=_mm256_add_epi16(cache[5], N[5]);
					cache[0]=_mm256_add_epi16(cache[0], W[0]);
					cache[1]=_mm256_add_epi16(cache[1], W[1]);
					cache[2]=_mm256_add_epi16(cache[2], W[2]);
					cache[3]=_mm256_add_epi16(cache[3], W[3]);
					cache[4]=_mm256_add_epi16(cache[4], W[4]);
					cache[5]=_mm256_add_epi16(cache[5], W[5]);
					cache[0]=_mm256_sub_epi16(cache[0], _mm256_slli_epi16(NW[0], 1));//-2*NW
					cache[1]=_mm256_sub_epi16(cache[1], _mm256_slli_epi16(NW[1], 1));
					cache[2]=_mm256_sub_epi16(cache[2], _mm256_slli_epi16(NW[2], 1));
					cache[3]=_mm256_sub_epi16(cache[3], _mm256_slli_epi16(NW[3], 1));
					cache[4]=_mm256_sub_epi16(cache[4], _mm256_slli_epi16(NW[4], 1));
					cache[5]=_mm256_sub_epi16(cache[5], _mm256_slli_epi16(NW[5], 1));
					wgpreds[5*6+0]=_mm256_srai_epi16(cache[0], 3);
					wgpreds[5*6+1]=_mm256_srai_epi16(cache[1], 3);
					wgpreds[5*6+2]=_mm256_srai_epi16(cache[2], 3);
					wgpreds[5*6+3]=_mm256_srai_epi16(cache[3], 3);
					wgpreds[5*6+4]=_mm256_srai_epi16(cache[4], 3);
					wgpreds[5*6+5]=_mm256_srai_epi16(cache[5], 3);
#endif

					//N+W-NW
					wgpreds[6*6+0]=predY0;
					wgpreds[6*6+1]=predY1;
					wgpreds[6*6+2]=predU0;
					wgpreds[6*6+3]=predU1;
					wgpreds[6*6+4]=predV0;
					wgpreds[6*6+5]=predV1;

					//N+NE-NNE
#if 1
					cache[0]=_mm256_add_epi16(N[0], _mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+0));//N+NE
					cache[1]=_mm256_add_epi16(N[1], _mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(N[2], _mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(N[3], _mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(N[4], _mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(N[5], _mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+1));
					wgpreds[7*6+0]=_mm256_sub_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[2]+(0+0+1*6)*NCODERS)+0));//NNE
					wgpreds[7*6+1]=_mm256_sub_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[2]+(0+0+1*6)*NCODERS)+1));
					wgpreds[7*6+2]=_mm256_sub_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[2]+(0+1+1*6)*NCODERS)+0));
					wgpreds[7*6+3]=_mm256_sub_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[2]+(0+1+1*6)*NCODERS)+1));
					wgpreds[7*6+4]=_mm256_sub_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[2]+(0+2+1*6)*NCODERS)+0));
					wgpreds[7*6+5]=_mm256_sub_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[2]+(0+2+1*6)*NCODERS)+1));
#endif
					//NE+NEE-NNEEE		X  bad
#if 0
					cache[0]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+0);//NE
					cache[1]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+1);
					cache[0]=_mm256_add_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[1]+(0+0+2*6)*NCODERS)+0));//+NEE
					cache[1]=_mm256_add_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[1]+(0+0+2*6)*NCODERS)+1));
					cache[2]=_mm256_add_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[1]+(0+1+2*6)*NCODERS)+0));
					cache[3]=_mm256_add_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[1]+(0+1+2*6)*NCODERS)+1));
					cache[4]=_mm256_add_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[1]+(0+2+2*6)*NCODERS)+0));
					cache[5]=_mm256_add_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[1]+(0+2+2*6)*NCODERS)+1));
					wgpreds[7*6+0]=_mm256_sub_epi16(cache[0], _mm256_load_si256((__m256i*)(rows[2]+(0+0+3*6)*NCODERS)+0));//-NNEEE
					wgpreds[7*6+1]=_mm256_sub_epi16(cache[1], _mm256_load_si256((__m256i*)(rows[2]+(0+0+3*6)*NCODERS)+1));
					wgpreds[7*6+2]=_mm256_sub_epi16(cache[2], _mm256_load_si256((__m256i*)(rows[2]+(0+1+3*6)*NCODERS)+0));
					wgpreds[7*6+3]=_mm256_sub_epi16(cache[3], _mm256_load_si256((__m256i*)(rows[2]+(0+1+3*6)*NCODERS)+1));
					wgpreds[7*6+4]=_mm256_sub_epi16(cache[4], _mm256_load_si256((__m256i*)(rows[2]+(0+2+3*6)*NCODERS)+0));
					wgpreds[7*6+5]=_mm256_sub_epi16(cache[5], _mm256_load_si256((__m256i*)(rows[2]+(0+2+3*6)*NCODERS)+1));
#endif

					//mix
					__m256i result[6];
					wg_mix(0, wgWerrors, erows[0], erows[1], wgpreds, result+0);
					wg_mix(1, wgWerrors, erows[0], erows[1], wgpreds, result+2);
					wg_mix(2, wgWerrors, erows[0], erows[1], wgpreds, result+4);
					predY0=result[0];
					predY1=result[1];
					predU0=result[2];
					predU1=result[3];
					predV0=result[4];
					predV1=result[5];

#if 0
					__m256d de[16];//de[P*2+highhalf]		prederrors[P*2+halflaneidx]
					de[0*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[0*2+0], 1));
					de[1*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[1*2+0], 1));
					de[2*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[2*2+0], 1));
					de[3*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[3*2+0], 1));
					de[4*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[4*2+0], 1));
					de[5*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[5*2+0], 1));
					de[6*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[6*2+0], 1));
					de[7*2+0]=_mm256_cvtepi32_pd(_mm256_extracti128_si256(prederrors[7*2+0], 1));
					de[0*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[0*2+0]));
					de[1*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[1*2+0]));
					de[2*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[2*2+0]));
					de[3*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[3*2+0]));
					de[4*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[4*2+0]));
					de[5*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[5*2+0]));
					de[6*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[6*2+0]));
					de[7*2+1]=_mm256_cvtepi32_pd(_mm256_castsi256_si128(prederrors[7*2+0]));
#endif

					//loosen pred range
					ymin0=_mm256_min_epi16(ymin0, NW[0]);
					ymax0=_mm256_max_epi16(ymax0, NW[0]);
					ymin1=_mm256_min_epi16(ymin1, NW[1]);
					ymax1=_mm256_max_epi16(ymax1, NW[1]);
					umin0=_mm256_min_epi16(umin0, NW[2]);
					umax0=_mm256_max_epi16(umax0, NW[2]);
					umin1=_mm256_min_epi16(umin1, NW[3]);
					umax1=_mm256_max_epi16(umax1, NW[3]);
					vmin0=_mm256_min_epi16(vmin0, NW[4]);
					vmax0=_mm256_max_epi16(vmax0, NW[4]);
					vmin1=_mm256_min_epi16(vmin1, NW[5]);
					vmax1=_mm256_max_epi16(vmax1, NW[5]);
					cache[0]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+0);//NE
					cache[1]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+1*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+1*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+1*6)*NCODERS)+1);
					ymin0=_mm256_min_epi16(ymin0, cache[0]);
					ymax0=_mm256_max_epi16(ymax0, cache[0]);
					ymin1=_mm256_min_epi16(ymin1, cache[1]);
					ymax1=_mm256_max_epi16(ymax1, cache[1]);
					umin0=_mm256_min_epi16(umin0, cache[2]);
					umax0=_mm256_max_epi16(umax0, cache[2]);
					umin1=_mm256_min_epi16(umin1, cache[3]);
					umax1=_mm256_max_epi16(umax1, cache[3]);
					vmin0=_mm256_min_epi16(vmin0, cache[4]);
					vmax0=_mm256_max_epi16(vmax0, cache[4]);
					vmin1=_mm256_min_epi16(vmin1, cache[5]);
					vmax1=_mm256_max_epi16(vmax1, cache[5]);
#if 1
					cache[0]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+3*6)*NCODERS)+0);//NEEE	crisp DIV2K -0.18%, noisy GDCC +0.13%
					cache[1]=_mm256_load_si256((__m256i*)(rows[1]+(0+0+3*6)*NCODERS)+1);
					cache[2]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+3*6)*NCODERS)+0);
					cache[3]=_mm256_load_si256((__m256i*)(rows[1]+(0+1+3*6)*NCODERS)+1);
					cache[4]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+3*6)*NCODERS)+0);
					cache[5]=_mm256_load_si256((__m256i*)(rows[1]+(0+2+3*6)*NCODERS)+1);
					//cache[0]=_mm256_srai_epi16(_mm256_add_epi16(NE[0], _mm256_load_si256((__m256i*)(rows[1]+(0+0+2*6)*NCODERS)+0)), 1);//(NE+NEE)>>1	X  bad
					//cache[1]=_mm256_srai_epi16(_mm256_add_epi16(NE[1], _mm256_load_si256((__m256i*)(rows[1]+(0+0+2*6)*NCODERS)+1)), 1);
					//cache[2]=_mm256_srai_epi16(_mm256_add_epi16(NE[2], _mm256_load_si256((__m256i*)(rows[1]+(0+1+2*6)*NCODERS)+0)), 1);
					//cache[3]=_mm256_srai_epi16(_mm256_add_epi16(NE[3], _mm256_load_si256((__m256i*)(rows[1]+(0+1+2*6)*NCODERS)+1)), 1);
					//cache[4]=_mm256_srai_epi16(_mm256_add_epi16(NE[4], _mm256_load_si256((__m256i*)(rows[1]+(0+2+2*6)*NCODERS)+0)), 1);
					//cache[5]=_mm256_srai_epi16(_mm256_add_epi16(NE[5], _mm256_load_si256((__m256i*)(rows[1]+(0+2+2*6)*NCODERS)+1)), 1);
					//cache[0]=_mm256_add_epi16(NE[0], _mm256_srai_epi16(_mm256_sub_epi16(_mm256_load_si256((__m256i*)(rows[1]+(0+0+2*6)*NCODERS)+0), NE[0]), 2));//NE+((NEE-NE)>>2)	X  worse
					//cache[1]=_mm256_add_epi16(NE[1], _mm256_srai_epi16(_mm256_sub_epi16(_mm256_load_si256((__m256i*)(rows[1]+(0+0+2*6)*NCODERS)+1), NE[1]), 2));
					//cache[2]=_mm256_add_epi16(NE[2], _mm256_srai_epi16(_mm256_sub_epi16(_mm256_load_si256((__m256i*)(rows[1]+(0+1+2*6)*NCODERS)+0), NE[2]), 2));
					//cache[3]=_mm256_add_epi16(NE[3], _mm256_srai_epi16(_mm256_sub_epi16(_mm256_load_si256((__m256i*)(rows[1]+(0+1+2*6)*NCODERS)+1), NE[3]), 2));
					//cache[4]=_mm256_add_epi16(NE[4], _mm256_srai_epi16(_mm256_sub_epi16(_mm256_load_si256((__m256i*)(rows[1]+(0+2+2*6)*NCODERS)+0), NE[4]), 2));
					//cache[5]=_mm256_add_epi16(NE[5], _mm256_srai_epi16(_mm256_sub_epi16(_mm256_load_si256((__m256i*)(rows[1]+(0+2+2*6)*NCODERS)+1), NE[5]), 2));
					ymin0=_mm256_min_epi16(ymin0, cache[0]);
					ymax0=_mm256_max_epi16(ymax0, cache[0]);
					ymin1=_mm256_min_epi16(ymin1, cache[1]);
					ymax1=_mm256_max_epi16(ymax1, cache[1]);
					umin0=_mm256_min_epi16(umin0, cache[2]);
					umax0=_mm256_max_epi16(umax0, cache[2]);
					umin1=_mm256_min_epi16(umin1, cache[3]);
					umax1=_mm256_max_epi16(umax1, cache[3]);
					vmin0=_mm256_min_epi16(vmin0, cache[4]);
					vmax0=_mm256_max_epi16(vmax0, cache[4]);
					vmin1=_mm256_min_epi16(vmin1, cache[5]);
					vmax1=_mm256_max_epi16(vmax1, cache[5]);
#endif
				}
				predY0=_mm256_max_epi16(predY0, ymin0);
				predY1=_mm256_max_epi16(predY1, ymin1);
				predU0=_mm256_max_epi16(predU0, umin0);
				predU1=_mm256_max_epi16(predU1, umin1);
				predV0=_mm256_max_epi16(predV0, vmin0);
				predV1=_mm256_max_epi16(predV1, vmin1);
				predY0=_mm256_min_epi16(predY0, ymax0);
				predY1=_mm256_min_epi16(predY1, ymax1);
				predU0=_mm256_min_epi16(predU0, umax0);
				predU1=_mm256_min_epi16(predU1, umax1);
				predV0=_mm256_min_epi16(predV0, vmax0);
				predV1=_mm256_min_epi16(predV1, vmax1);
			}

			__m256i msyms0, msyms1, moffset0, moffset1;
			if(fwd)
			{
				__m256i ctxblendmask=_mm256_set1_epi16(255);
#if 0
				{
					__m128i y0=_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+yidx)+0), half8);
					__m128i y1=_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+yidx)+1), half8);
					__m128i u0=_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+uidx)+0), half8);
					__m128i u1=_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+uidx)+1), half8);
					__m128i v0=_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+vidx)+0), half8);
					__m128i v1=_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+vidx)+1), half8);
					if(combination[II_COEFF_V_SUB_Y]+combination[II_COEFF_V_SUB_U])
					{
						__m128i t0, t1;
						switch(combination[II_COEFF_V_SUB_U])
						{
						case 0:
							v0=_mm_sub_epi8(v0, y0);
							v1=_mm_sub_epi8(v1, y1);
							break;
						case 1:
							break;
						case 2:
							v0=_mm_sub_epi8(v0, y0);
							v1=_mm_sub_epi8(v1, y1);
							break;
						case 3:
							break;
						case 4:
							v0=_mm_sub_epi8(v0, u0);
							v1=_mm_sub_epi8(v1, u1);
							break;
						}
					}
					if(combination[II_COEFF_U_SUB_Y])
					{
						u0=_mm_sub_epi8(u0, y0);
						u1=_mm_sub_epi8(u1, y1);
					}
				}
#endif
				myuv[0]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+yidx)+0), half8));//load yuv
				myuv[1]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+yidx)+1), half8));
				myuv[2]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+uidx)+0), half8));
				myuv[3]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+uidx)+1), half8));
				myuv[4]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+vidx)+0), half8));
				myuv[5]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+vidx)+1), half8));
#ifdef ENABLE_MA
				//v-=(combination[II_COEFF_V_SUB_Y]*y+combination[II_COEFF_V_SUB_U]*u)>>2;
				moffset0=_mm256_mullo_epi16(vc0, myuv[0]);
				moffset1=_mm256_mullo_epi16(vc0, myuv[1]);
				moffset0=_mm256_add_epi16(moffset0, _mm256_mullo_epi16(vc1, myuv[2]));
				moffset1=_mm256_add_epi16(moffset1, _mm256_mullo_epi16(vc1, myuv[3]));
				moffset0=_mm256_srai_epi16(moffset0, 2);
				moffset1=_mm256_srai_epi16(moffset1, 2);
				myuv[4]=_mm256_sub_epi16(myuv[4], moffset0);
				myuv[5]=_mm256_sub_epi16(myuv[5], moffset1);
				myuv[4]=_mm256_slli_epi16(myuv[4], 8);
				myuv[5]=_mm256_slli_epi16(myuv[5], 8);
				myuv[4]=_mm256_srai_epi16(myuv[4], 8);
				myuv[5]=_mm256_srai_epi16(myuv[5], 8);

				//u-=y&-(combination[II_COEFF_U_SUB_Y]!=0);
				moffset0=_mm256_and_si256(myuv[0], uhelpmask);
				moffset1=_mm256_and_si256(myuv[1], uhelpmask);
				myuv[2]=_mm256_sub_epi16(myuv[2], moffset0);
				myuv[3]=_mm256_sub_epi16(myuv[3], moffset1);
				myuv[2]=_mm256_slli_epi16(myuv[2], 8);
				myuv[3]=_mm256_slli_epi16(myuv[3], 8);
				myuv[2]=_mm256_srai_epi16(myuv[2], 8);
				myuv[3]=_mm256_srai_epi16(myuv[3], 8);
				moffset1=moffset0=_mm256_setzero_si256();
#endif

				W[0]=myuv[0];
				W[1]=myuv[1];
				_mm256_store_si256((__m256i*)(rows[0]+(0+0+0*6)*NCODERS)+0, myuv[0]);//store Y neighbors
				_mm256_store_si256((__m256i*)(rows[0]+(0+0+0*6)*NCODERS)+1, myuv[1]);
				msyms0=_mm256_sub_epi16(myuv[0], predY0);//sub pred
				msyms1=_mm256_sub_epi16(myuv[1], predY1);
				ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));//ecurr = pack_sign(yuv-pred)		FIXME try abs
				ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms0=_mm256_sub_epi16(msyms0, amin);
				msyms1=_mm256_sub_epi16(msyms1, amin);
				ctxY0=_mm256_slli_epi16(ctxY0, 8);
				ctxY1=_mm256_slli_epi16(ctxY1, 8);
				ctxY0=_mm256_blendv_epi8(ctxY0, msyms0, ctxblendmask);
				ctxY1=_mm256_blendv_epi8(ctxY1, msyms1, ctxblendmask);
				_mm256_store_si256((__m256i*)syms+0, ctxY0);
				_mm256_store_si256((__m256i*)syms+1, ctxY1);
				_mm256_store_si256((__m256i*)(ctxptr+NCODERS*(0*2+0)), ctxY0);//store Y  ctx|residuals		ctxptr+NCODERS*(C*2+R)
				_mm256_store_si256((__m256i*)(ctxptr+NCODERS*(0*2+1)), ctxY1);
				
#ifndef ENABLE_MA
				moffset0=_mm256_and_si256(myuv[0], uhelpmask);
				moffset1=_mm256_and_si256(myuv[1], uhelpmask);
				predU0=_mm256_add_epi16(predU0, moffset0);
				predU1=_mm256_add_epi16(predU1, moffset1);
				predU0=_mm256_max_epi16(predU0, amin);
				predU1=_mm256_max_epi16(predU1, amin);
				predU0=_mm256_min_epi16(predU0, amax);
				predU1=_mm256_min_epi16(predU1, amax);
#endif
				msyms0=_mm256_sub_epi16(myuv[2], predU0);
				msyms1=_mm256_sub_epi16(myuv[3], predU1);
				ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));
				ecurr[3]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms0=_mm256_sub_epi16(msyms0, amin);
				msyms1=_mm256_sub_epi16(msyms1, amin);
				ctxU0=_mm256_add_epi16(ctxU0, mctxuoffset);
				ctxU1=_mm256_add_epi16(ctxU1, mctxuoffset);
				ctxU0=_mm256_slli_epi16(ctxU0, 8);
				ctxU1=_mm256_slli_epi16(ctxU1, 8);
				ctxU0=_mm256_blendv_epi8(ctxU0, msyms0, ctxblendmask);
				ctxU1=_mm256_blendv_epi8(ctxU1, msyms1, ctxblendmask);
				_mm256_store_si256((__m256i*)syms+2, ctxU0);
				_mm256_store_si256((__m256i*)syms+3, ctxU1);
				_mm256_store_si256((__m256i*)(ctxptr+NCODERS*(1*2+0)), ctxU0);//store U  ctx|residuals		ctxptr+NCODERS*(C*2+R)
				_mm256_store_si256((__m256i*)(ctxptr+NCODERS*(1*2+1)), ctxU1);
#ifdef ENABLE_MA
				W[2]=myuv[2];
				W[3]=myuv[3];
				_mm256_store_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+0, myuv[2]);
				_mm256_store_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+1, myuv[3]);
#else
				msyms0=_mm256_sub_epi16(myuv[2], moffset0);
				msyms1=_mm256_sub_epi16(myuv[3], moffset1);
				W[2]=msyms0;
				W[3]=msyms1;
				_mm256_store_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+0, msyms0);//store U neighbors
				_mm256_store_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+1, msyms1);
#endif
				
#ifndef ENABLE_MA
				moffset0=_mm256_mullo_epi16(vc0, myuv[0]);
				moffset1=_mm256_mullo_epi16(vc0, myuv[1]);
				moffset0=_mm256_add_epi16(moffset0, _mm256_mullo_epi16(vc1, myuv[2]));
				moffset1=_mm256_add_epi16(moffset1, _mm256_mullo_epi16(vc1, myuv[3]));
				predV0=_mm256_add_epi16(predV0, moffset0);
				predV1=_mm256_add_epi16(predV1, moffset1);
				predV0=_mm256_srai_epi16(predV0, 2);
				predV1=_mm256_srai_epi16(predV1, 2);
				predV0=_mm256_max_epi16(predV0, amin);
				predV1=_mm256_max_epi16(predV1, amin);
				predV0=_mm256_min_epi16(predV0, amax);
				predV1=_mm256_min_epi16(predV1, amax);
#endif
				msyms0=_mm256_sub_epi16(myuv[4], predV0);
				msyms1=_mm256_sub_epi16(myuv[5], predV1);
				ecurr[4]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));
				ecurr[5]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms0=_mm256_sub_epi16(msyms0, amin);
				msyms1=_mm256_sub_epi16(msyms1, amin);
				ctxV0=_mm256_add_epi16(ctxV0, mctxvoffset);
				ctxV1=_mm256_add_epi16(ctxV1, mctxvoffset);
				ctxV0=_mm256_slli_epi16(ctxV0, 8);
				ctxV1=_mm256_slli_epi16(ctxV1, 8);
				ctxV0=_mm256_blendv_epi8(ctxV0, msyms0, ctxblendmask);
				ctxV1=_mm256_blendv_epi8(ctxV1, msyms1, ctxblendmask);
				_mm256_store_si256((__m256i*)syms+4, ctxV0);
				_mm256_store_si256((__m256i*)syms+5, ctxV1);
				_mm256_store_si256((__m256i*)(ctxptr+NCODERS*(2*2+0)), ctxV0);//store V  ctx|residuals		ctxptr+NCODERS*(C*2+R)
				_mm256_store_si256((__m256i*)(ctxptr+NCODERS*(2*2+1)), ctxV1);
#ifdef ENABLE_MA
				W[4]=myuv[4];
				W[5]=myuv[5];
				_mm256_store_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+0, myuv[4]);
				_mm256_store_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+1, myuv[5]);
#else
				msyms0=_mm256_slli_epi16(myuv[4], 2);
				msyms1=_mm256_slli_epi16(myuv[5], 2);
				msyms0=_mm256_sub_epi16(msyms0, moffset0);
				msyms1=_mm256_sub_epi16(msyms1, moffset1);
				W[4]=msyms0;
				W[5]=msyms1;
				_mm256_store_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+0, msyms0);//store V neighbors
				_mm256_store_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+1, msyms1);
#endif
				++hists[syms[0*32+0x00]];
				++hists[syms[1*32+0x00]];
				++hists[syms[2*32+0x00]];
				++hists[syms[0*32+0x01]];
				++hists[syms[1*32+0x01]];
				++hists[syms[2*32+0x01]];
				++hists[syms[0*32+0x02]];
				++hists[syms[1*32+0x02]];
				++hists[syms[2*32+0x02]];
				++hists[syms[0*32+0x03]];
				++hists[syms[1*32+0x03]];
				++hists[syms[2*32+0x03]];
				++hists[syms[0*32+0x04]];
				++hists[syms[1*32+0x04]];
				++hists[syms[2*32+0x04]];
				++hists[syms[0*32+0x05]];
				++hists[syms[1*32+0x05]];
				++hists[syms[2*32+0x05]];
				++hists[syms[0*32+0x06]];
				++hists[syms[1*32+0x06]];
				++hists[syms[2*32+0x06]];
				++hists[syms[0*32+0x07]];
				++hists[syms[1*32+0x07]];
				++hists[syms[2*32+0x07]];
				++hists[syms[0*32+0x08]];
				++hists[syms[1*32+0x08]];
				++hists[syms[2*32+0x08]];
				++hists[syms[0*32+0x09]];
				++hists[syms[1*32+0x09]];
				++hists[syms[2*32+0x09]];
				++hists[syms[0*32+0x0A]];
				++hists[syms[1*32+0x0A]];
				++hists[syms[2*32+0x0A]];
				++hists[syms[0*32+0x0B]];
				++hists[syms[1*32+0x0B]];
				++hists[syms[2*32+0x0B]];
				++hists[syms[0*32+0x0C]];
				++hists[syms[1*32+0x0C]];
				++hists[syms[2*32+0x0C]];
				++hists[syms[0*32+0x0D]];
				++hists[syms[1*32+0x0D]];
				++hists[syms[2*32+0x0D]];
				++hists[syms[0*32+0x0E]];
				++hists[syms[1*32+0x0E]];
				++hists[syms[2*32+0x0E]];
				++hists[syms[0*32+0x0F]];
				++hists[syms[1*32+0x0F]];
				++hists[syms[2*32+0x0F]];
				++hists[syms[0*32+0x10]];
				++hists[syms[1*32+0x10]];
				++hists[syms[2*32+0x10]];
				++hists[syms[0*32+0x11]];
				++hists[syms[1*32+0x11]];
				++hists[syms[2*32+0x11]];
				++hists[syms[0*32+0x12]];
				++hists[syms[1*32+0x12]];
				++hists[syms[2*32+0x12]];
				++hists[syms[0*32+0x13]];
				++hists[syms[1*32+0x13]];
				++hists[syms[2*32+0x13]];
				++hists[syms[0*32+0x14]];
				++hists[syms[1*32+0x14]];
				++hists[syms[2*32+0x14]];
				++hists[syms[0*32+0x15]];
				++hists[syms[1*32+0x15]];
				++hists[syms[2*32+0x15]];
				++hists[syms[0*32+0x16]];
				++hists[syms[1*32+0x16]];
				++hists[syms[2*32+0x16]];
				++hists[syms[0*32+0x17]];
				++hists[syms[1*32+0x17]];
				++hists[syms[2*32+0x17]];
				++hists[syms[0*32+0x18]];
				++hists[syms[1*32+0x18]];
				++hists[syms[2*32+0x18]];
				++hists[syms[0*32+0x19]];
				++hists[syms[1*32+0x19]];
				++hists[syms[2*32+0x19]];
				++hists[syms[0*32+0x1A]];
				++hists[syms[1*32+0x1A]];
				++hists[syms[2*32+0x1A]];
				++hists[syms[0*32+0x1B]];
				++hists[syms[1*32+0x1B]];
				++hists[syms[2*32+0x1B]];
				++hists[syms[0*32+0x1C]];
				++hists[syms[1*32+0x1C]];
				++hists[syms[2*32+0x1C]];
				++hists[syms[0*32+0x1D]];
				++hists[syms[1*32+0x1D]];
				++hists[syms[2*32+0x1D]];
				++hists[syms[0*32+0x1E]];
				++hists[syms[1*32+0x1E]];
				++hists[syms[2*32+0x1E]];
				++hists[syms[0*32+0x1F]];
				++hists[syms[1*32+0x1F]];
				++hists[syms[2*32+0x1F]];
				ctxptr+=sizeof(short[3][NCODERS]);
			}
			else
			{
				//decode main
				__m256i msyms8;
				
				//decode Y
				dec_yuv(mstate, 0, &ctxY0, &ctxY1, (int*)CDF2syms, ans_permute, &streamptr, streamend, myuv+0*2);//residuals from [0 ~ 255]

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
				myuv[0]=_mm256_add_epi16(myuv[0], predY0);
				myuv[1]=_mm256_add_epi16(myuv[1], predY1);
				myuv[0]=_mm256_and_si256(myuv[0], bytemask);
				myuv[1]=_mm256_and_si256(myuv[1], bytemask);
				myuv[0]=_mm256_add_epi16(myuv[0], amin);
				myuv[1]=_mm256_add_epi16(myuv[1], amin);
				msyms0=_mm256_sub_epi16(myuv[0], predY0);//sub pred
				msyms1=_mm256_sub_epi16(myuv[1], predY1);
				ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));
				ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms8=_mm256_packs_epi16(myuv[0], myuv[1]);
				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
				msyms8=_mm256_xor_si256(msyms8, _mm256_set1_epi8(-128));
				_mm256_store_si256((__m256i*)(imptr+yidx), msyms8);//store Y bytes
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
					LOG_ERROR("guide error XYC0 %d %d %d/%d", kx, ky, yidx, NCODERS);
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
				W[1]=myuv[1];
				_mm256_store_si256((__m256i*)(rows[0]+(0+0+0*6)*NCODERS)+0, myuv[0]);//store Y neighbors
				_mm256_store_si256((__m256i*)(rows[0]+(0+0+0*6)*NCODERS)+1, myuv[1]);


				//decode U
#ifndef ENABLE_MA
				moffset0=_mm256_and_si256(myuv[0], uhelpmask);
				moffset1=_mm256_and_si256(myuv[1], uhelpmask);
				predU0=_mm256_add_epi16(predU0, moffset0);
				predU1=_mm256_add_epi16(predU1, moffset1);
				predU0=_mm256_max_epi16(predU0, amin);
				predU1=_mm256_max_epi16(predU1, amin);
				predU0=_mm256_min_epi16(predU0, amax);
				predU1=_mm256_min_epi16(predU1, amax);
#endif
				dec_yuv(mstate, 1, &ctxU0, &ctxU1, (int*)CDF2syms, ans_permute, &streamptr, streamend, myuv+1*2);
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
				
				myuv[2]=_mm256_add_epi16(myuv[2], predU0);
				myuv[3]=_mm256_add_epi16(myuv[3], predU1);
				myuv[2]=_mm256_and_si256(myuv[2], bytemask);
				myuv[3]=_mm256_and_si256(myuv[3], bytemask);
				myuv[2]=_mm256_add_epi16(myuv[2], amin);
				myuv[3]=_mm256_add_epi16(myuv[3], amin);
				msyms0=_mm256_sub_epi16(myuv[2], predU0);//sub pred
				msyms1=_mm256_sub_epi16(myuv[3], predU1);
				ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));
				ecurr[3]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
#ifdef ENABLE_MA
				__m256i mrgb[6];
				//u-=y&-(combination[II_COEFF_U_SUB_Y]!=0);
				moffset0=_mm256_and_si256(myuv[0], uhelpmask);
				moffset1=_mm256_and_si256(myuv[1], uhelpmask);
				mrgb[2]=_mm256_add_epi16(myuv[2], moffset0);
				mrgb[3]=_mm256_add_epi16(myuv[3], moffset1);
				mrgb[2]=_mm256_slli_epi16(mrgb[2], 8);
				mrgb[3]=_mm256_slli_epi16(mrgb[3], 8);
				mrgb[2]=_mm256_srai_epi16(mrgb[2], 8);
				mrgb[3]=_mm256_srai_epi16(mrgb[3], 8);
				moffset1=moffset0=_mm256_setzero_si256();
				msyms8=_mm256_packs_epi16(mrgb[2], mrgb[3]);
#else
				msyms8=_mm256_packs_epi16(myuv[2], myuv[3]);
#endif
				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
				msyms8=_mm256_xor_si256(msyms8, _mm256_set1_epi8(-128));
				_mm256_store_si256((__m256i*)(imptr+uidx), msyms8);//store U bytes
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
					LOG_ERROR("guide error XYC1 %d %d %d/%d", kx, ky, uidx, NCODERS);
				}
#endif
//#ifdef ANS_VAL
//				msyms8=_mm256_packus_epi16(ctxU0, ctxU1);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_storeu_si256((__m256i*)actx+1, msyms8);
//				ansval_check(actx+1*NCODERS, 1, NCODERS);
//				ansval_check(imptr+kx+uidx, 1, NCODERS);
//#endif
#ifdef ENABLE_MA
				W[2]=myuv[2];
				W[3]=myuv[3];
				_mm256_store_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+0, myuv[2]);
				_mm256_store_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+1, myuv[3]);
#else
				W[2]=_mm256_sub_epi16(myuv[2], moffset0);//subtract Uoffset from U
				W[3]=_mm256_sub_epi16(myuv[3], moffset1);
				_mm256_store_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+0, W[2]);//store U neighbors
				_mm256_store_si256((__m256i*)(rows[0]+(0+1+0*6)*NCODERS)+1, W[3]);
#endif
				

				//decode V
#ifndef ENABLE_MA
				moffset0=_mm256_mullo_epi16(vc0, myuv[0]);
				moffset1=_mm256_mullo_epi16(vc0, myuv[1]);
				moffset0=_mm256_add_epi16(moffset0, _mm256_mullo_epi16(vc1, myuv[2]));
				moffset1=_mm256_add_epi16(moffset1, _mm256_mullo_epi16(vc1, myuv[3]));
				predV0=_mm256_add_epi16(predV0, moffset0);
				predV1=_mm256_add_epi16(predV1, moffset1);
				predV0=_mm256_srai_epi16(predV0, 2);
				predV1=_mm256_srai_epi16(predV1, 2);
				predV0=_mm256_max_epi16(predV0, amin);
				predV1=_mm256_max_epi16(predV1, amin);
				predV0=_mm256_min_epi16(predV0, amax);
				predV1=_mm256_min_epi16(predV1, amax);
#endif
				dec_yuv(mstate, 2, &ctxV0, &ctxV1, (int*)CDF2syms, ans_permute, &streamptr, streamend, myuv+2*2);
//#ifdef ANS_VAL
//				msyms8=_mm256_packus_epi16(ctxV0, ctxV1);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_store_si256((__m256i*)debugvals+4, msyms8);
//				msyms8=_mm256_packus_epi16(myuv[4], myuv[5]);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_store_si256((__m256i*)debugvals+5, msyms8);
//				ansval_check(debugvals+4*NCODERS, 1, NCODERS);//contexts
//				ansval_check(debugvals+5*NCODERS, 1, NCODERS);//residuals
//#endif
				
				myuv[4]=_mm256_add_epi16(myuv[4], predV0);
				myuv[5]=_mm256_add_epi16(myuv[5], predV1);
				myuv[4]=_mm256_and_si256(myuv[4], bytemask);
				myuv[5]=_mm256_and_si256(myuv[5], bytemask);
				myuv[4]=_mm256_add_epi16(myuv[4], amin);
				myuv[5]=_mm256_add_epi16(myuv[5], amin);
				msyms0=_mm256_sub_epi16(myuv[4], predV0);
				msyms1=_mm256_sub_epi16(myuv[5], predV1);
				ecurr[4]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));
				ecurr[5]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
#ifdef ENABLE_MA
				//v-=(combination[II_COEFF_V_SUB_Y]*y+combination[II_COEFF_V_SUB_U]*u)>>2;
				moffset0=_mm256_mullo_epi16(vc0, myuv[0]);
				moffset1=_mm256_mullo_epi16(vc0, myuv[1]);
				moffset0=_mm256_add_epi16(moffset0, _mm256_mullo_epi16(vc1, mrgb[2]));
				moffset1=_mm256_add_epi16(moffset1, _mm256_mullo_epi16(vc1, mrgb[3]));
				moffset0=_mm256_srai_epi16(moffset0, 2);
				moffset1=_mm256_srai_epi16(moffset1, 2);
				mrgb[4]=_mm256_add_epi16(myuv[4], moffset0);
				mrgb[5]=_mm256_add_epi16(myuv[5], moffset1);
				mrgb[4]=_mm256_slli_epi16(mrgb[4], 8);
				mrgb[5]=_mm256_slli_epi16(mrgb[5], 8);
				mrgb[4]=_mm256_srai_epi16(mrgb[4], 8);
				mrgb[5]=_mm256_srai_epi16(mrgb[5], 8);
				moffset1=moffset0=_mm256_setzero_si256();
				msyms8=_mm256_packs_epi16(mrgb[4], mrgb[5]);
#else
				msyms8=_mm256_packs_epi16(myuv[4], myuv[5]);
#endif
				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
				msyms8=_mm256_xor_si256(msyms8, _mm256_set1_epi8(-128));
				_mm256_store_si256((__m256i*)(imptr+vidx), msyms8);//store V bytes
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
					LOG_ERROR("guide error XYC2 %d %d %d/%d", kx, ky, vidx, NCODERS);
				}
#endif
//#ifdef ANS_VAL
//				msyms8=_mm256_packus_epi16(ctxV0, ctxV1);
//				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
//				_mm256_storeu_si256((__m256i*)actx, msyms8);
//				ansval_check(actx+2*NCODERS, 1, NCODERS);
//				ansval_check(imptr+kx+vidx, 1, NCODERS);
//#endif
#ifdef ENABLE_MA
				W[4]=myuv[4];
				W[5]=myuv[5];
				_mm256_store_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+0, myuv[4]);
				_mm256_store_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+1, myuv[5]);
#else
				W[4]=_mm256_slli_epi16(myuv[4], 2);
				W[5]=_mm256_slli_epi16(myuv[5], 2);
				W[4]=_mm256_sub_epi16(W[4], moffset0);//subtract Voffset from V
				W[5]=_mm256_sub_epi16(W[5], moffset1);
				_mm256_store_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+0, W[4]);//store V neighbors
				_mm256_store_si256((__m256i*)(rows[0]+(0+2+0*6)*NCODERS)+1, W[5]);
#endif
				
//#ifdef ENABLE_GUIDE
//				for(int kx2=0;kx2<3*NCODERS;kx2+=3)
//					guide_check(image, kx+kx2, ky);
//#endif
			}
			if(use_wg4)//update
			{
#if 0
				for(int kp=0;kp<8;++kp)
				{
					__m256i t[8];

					//e = abs(curr-wgpred[i])<<1
					__m256i me0=_mm256_sub_epi16(W[0], wgpreds[kp*6+0]);//W is now (YUV-offset)
					__m256i me1=_mm256_sub_epi16(W[1], wgpreds[kp*6+1]);
					__m256i me2=_mm256_sub_epi16(W[2], wgpreds[kp*6+2]);
					__m256i me3=_mm256_sub_epi16(W[3], wgpreds[kp*6+3]);
					__m256i me4=_mm256_sub_epi16(W[4], wgpreds[kp*6+4]);
					__m256i me5=_mm256_sub_epi16(W[5], wgpreds[kp*6+5]);
					me0=_mm256_abs_epi16(me0);
					me1=_mm256_abs_epi16(me1);
					me2=_mm256_abs_epi16(me2);
					me3=_mm256_abs_epi16(me3);
					me4=_mm256_abs_epi16(me4);
					me5=_mm256_abs_epi16(me5);
					me0=_mm256_slli_epi16(me0, 1);
					me1=_mm256_slli_epi16(me1, 1);
					me2=_mm256_slli_epi16(me2, 1);
					me3=_mm256_slli_epi16(me3, 1);
					me4=_mm256_slli_epi16(me4, 1);
					me5=_mm256_slli_epi16(me5, 1);
				//	me0=_mm256_xor_si256(_mm256_slli_epi16(me0, 1), _mm256_srai_epi16(me0, 15));//X
				//	me1=_mm256_xor_si256(_mm256_slli_epi16(me1, 1), _mm256_srai_epi16(me1, 15));
				//	me2=_mm256_xor_si256(_mm256_slli_epi16(me2, 1), _mm256_srai_epi16(me2, 15));
				//	me3=_mm256_xor_si256(_mm256_slli_epi16(me3, 1), _mm256_srai_epi16(me3, 15));
				//	me4=_mm256_xor_si256(_mm256_slli_epi16(me4, 1), _mm256_srai_epi16(me4, 15));
				//	me5=_mm256_xor_si256(_mm256_slli_epi16(me5, 1), _mm256_srai_epi16(me5, 15));

					//eW += (e-eW+(1<<7>>1))>>7
#ifndef WG_ENABLE_eW
					(void)wgWerrors;
#else
					{
						__m256i bias=_mm256_set1_epi16(1<<4>>1);
						__m256i peW[6];
						peW[0]=_mm256_load_si256((__m256i*)wgWerrors+kp*6+0);
						peW[1]=_mm256_load_si256((__m256i*)wgWerrors+kp*6+1);
						peW[2]=_mm256_load_si256((__m256i*)wgWerrors+kp*6+2);
						peW[3]=_mm256_load_si256((__m256i*)wgWerrors+kp*6+3);
						peW[4]=_mm256_load_si256((__m256i*)wgWerrors+kp*6+4);
						peW[5]=_mm256_load_si256((__m256i*)wgWerrors+kp*6+5);
						t[0]=_mm256_sub_epi16(me0, peW[0]);//eW += (e-eW+(1<<7>>1))>>7
						t[1]=_mm256_sub_epi16(me1, peW[1]);
						t[2]=_mm256_sub_epi16(me2, peW[2]);
						t[3]=_mm256_sub_epi16(me3, peW[3]);
						t[4]=_mm256_sub_epi16(me4, peW[4]);
						t[5]=_mm256_sub_epi16(me5, peW[5]);
						t[0]=_mm256_add_epi16(t[0], bias);
						t[1]=_mm256_add_epi16(t[1], bias);
						t[2]=_mm256_add_epi16(t[2], bias);
						t[3]=_mm256_add_epi16(t[3], bias);
						t[4]=_mm256_add_epi16(t[4], bias);
						t[5]=_mm256_add_epi16(t[5], bias);
						t[0]=_mm256_srai_epi16(t[0], 4);
						t[1]=_mm256_srai_epi16(t[1], 4);
						t[2]=_mm256_srai_epi16(t[2], 4);
						t[3]=_mm256_srai_epi16(t[3], 4);
						t[4]=_mm256_srai_epi16(t[4], 4);
						t[5]=_mm256_srai_epi16(t[5], 4);
						peW[0]=_mm256_add_epi16(peW[0], t[0]);
						peW[1]=_mm256_add_epi16(peW[1], t[1]);
						peW[2]=_mm256_add_epi16(peW[2], t[2]);
						peW[3]=_mm256_add_epi16(peW[3], t[3]);
						peW[4]=_mm256_add_epi16(peW[4], t[4]);
						peW[5]=_mm256_add_epi16(peW[5], t[5]);
						_mm256_store_si256((__m256i*)wgWerrors+kp*6+0, peW[0]);
						_mm256_store_si256((__m256i*)wgWerrors+kp*6+1, peW[1]);
						_mm256_store_si256((__m256i*)wgWerrors+kp*6+2, peW[2]);
						_mm256_store_si256((__m256i*)wgWerrors+kp*6+3, peW[3]);
						_mm256_store_si256((__m256i*)wgWerrors+kp*6+4, peW[4]);
						_mm256_store_si256((__m256i*)wgWerrors+kp*6+5, peW[5]);
					}
#endif
					//				(__m256i*)(erows[-Y]+((P+X*NPREDS)*3+C)*NCODERS)+R
					__m256i eNW0	=_mm256_load_si256((__m256i*)(erows[1]+((kp-1*WG_NPREDS)*3+0)*NCODERS)+0);//eNW
					__m256i eNW1	=_mm256_load_si256((__m256i*)(erows[1]+((kp-1*WG_NPREDS)*3+0)*NCODERS)+1);
					__m256i eNW2	=_mm256_load_si256((__m256i*)(erows[1]+((kp-1*WG_NPREDS)*3+1)*NCODERS)+0);
					__m256i eNW3	=_mm256_load_si256((__m256i*)(erows[1]+((kp-1*WG_NPREDS)*3+1)*NCODERS)+1);
					__m256i eNW4	=_mm256_load_si256((__m256i*)(erows[1]+((kp-1*WG_NPREDS)*3+2)*NCODERS)+0);
					__m256i eNW5	=_mm256_load_si256((__m256i*)(erows[1]+((kp-1*WG_NPREDS)*3+2)*NCODERS)+1);
					__m256i eN0	=_mm256_load_si256((__m256i*)(erows[1]+((kp+0*WG_NPREDS)*3+0)*NCODERS)+0);//eN
					__m256i eN1	=_mm256_load_si256((__m256i*)(erows[1]+((kp+0*WG_NPREDS)*3+0)*NCODERS)+1);
					__m256i eN2	=_mm256_load_si256((__m256i*)(erows[1]+((kp+0*WG_NPREDS)*3+1)*NCODERS)+0);
					__m256i eN3	=_mm256_load_si256((__m256i*)(erows[1]+((kp+0*WG_NPREDS)*3+1)*NCODERS)+1);
					__m256i eN4	=_mm256_load_si256((__m256i*)(erows[1]+((kp+0*WG_NPREDS)*3+2)*NCODERS)+0);
					__m256i eN5	=_mm256_load_si256((__m256i*)(erows[1]+((kp+0*WG_NPREDS)*3+2)*NCODERS)+1);
					__m256i eNE0	=_mm256_load_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+0)*NCODERS)+0);//eNE
					__m256i eNE1	=_mm256_load_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+0)*NCODERS)+1);
					__m256i eNE2	=_mm256_load_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+1)*NCODERS)+0);
					__m256i eNE3	=_mm256_load_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+1)*NCODERS)+1);
					__m256i eNE4	=_mm256_load_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+2)*NCODERS)+0);
					__m256i eNE5	=_mm256_load_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+2)*NCODERS)+1);
					__m256i eNEE0	=_mm256_load_si256((__m256i*)(erows[1]+((kp+2*WG_NPREDS)*3+0)*NCODERS)+0);//eNEE
					__m256i eNEE1	=_mm256_load_si256((__m256i*)(erows[1]+((kp+2*WG_NPREDS)*3+0)*NCODERS)+1);
					__m256i eNEE2	=_mm256_load_si256((__m256i*)(erows[1]+((kp+2*WG_NPREDS)*3+1)*NCODERS)+0);
					__m256i eNEE3	=_mm256_load_si256((__m256i*)(erows[1]+((kp+2*WG_NPREDS)*3+1)*NCODERS)+1);
					__m256i eNEE4	=_mm256_load_si256((__m256i*)(erows[1]+((kp+2*WG_NPREDS)*3+2)*NCODERS)+0);
					__m256i eNEE5	=_mm256_load_si256((__m256i*)(erows[1]+((kp+2*WG_NPREDS)*3+2)*NCODERS)+1);

					//eNE += e
					eNE0=_mm256_add_epi16(eNE0, me0);
					eNE1=_mm256_add_epi16(eNE1, me1);
					eNE2=_mm256_add_epi16(eNE2, me2);
					eNE3=_mm256_add_epi16(eNE3, me3);
					eNE4=_mm256_add_epi16(eNE4, me4);
					eNE5=_mm256_add_epi16(eNE5, me5);
				//	eNE0=_mm256_srli_epi16(eNE0, 1);
				//	eNE1=_mm256_srli_epi16(eNE1, 1);
				//	eNE2=_mm256_srli_epi16(eNE2, 1);
				//	eNE3=_mm256_srli_epi16(eNE3, 1);
				//	eNE4=_mm256_srli_epi16(eNE4, 1);
				//	eNE5=_mm256_srli_epi16(eNE5, 1);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+0)*NCODERS)+0, eNE0);//eNE
					_mm256_store_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+0)*NCODERS)+1, eNE1);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+1)*NCODERS)+0, eNE2);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+1)*NCODERS)+1, eNE3);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+2)*NCODERS)+0, eNE4);
					_mm256_store_si256((__m256i*)(erows[1]+((kp+1*WG_NPREDS)*3+2)*NCODERS)+1, eNE5);

					t[0]=_mm256_load_si256((__m256i*)(erows[0]+((kp-1*WG_NPREDS)*3+0)*NCODERS)+0);//eW
					t[1]=_mm256_load_si256((__m256i*)(erows[0]+((kp-1*WG_NPREDS)*3+0)*NCODERS)+1);
					t[2]=_mm256_load_si256((__m256i*)(erows[0]+((kp-1*WG_NPREDS)*3+1)*NCODERS)+0);
					t[3]=_mm256_load_si256((__m256i*)(erows[0]+((kp-1*WG_NPREDS)*3+1)*NCODERS)+1);
					t[4]=_mm256_load_si256((__m256i*)(erows[0]+((kp-1*WG_NPREDS)*3+2)*NCODERS)+0);
					t[5]=_mm256_load_si256((__m256i*)(erows[0]+((kp-1*WG_NPREDS)*3+2)*NCODERS)+1);
#if 1
					//curr = (eW + min(eN, eW) + e + min(eN, eNEE))>>2
					//curr = (2*eW+e+eNE)>>2	X
					t[0]=_mm256_add_epi32(t[0], _mm256_min_epi16(eNW0, eN0));
					t[1]=_mm256_add_epi32(t[1], _mm256_min_epi16(eNW1, eN1));
					t[2]=_mm256_add_epi32(t[2], _mm256_min_epi16(eNW2, eN2));
					t[3]=_mm256_add_epi32(t[3], _mm256_min_epi16(eNW3, eN3));
					t[4]=_mm256_add_epi32(t[4], _mm256_min_epi16(eNW4, eN4));
					t[5]=_mm256_add_epi32(t[5], _mm256_min_epi16(eNW5, eN5));
					t[0]=_mm256_add_epi16(t[0], me0);
					t[1]=_mm256_add_epi16(t[1], me1);
					t[2]=_mm256_add_epi16(t[2], me2);
					t[3]=_mm256_add_epi16(t[3], me3);
					t[4]=_mm256_add_epi16(t[4], me4);
					t[5]=_mm256_add_epi16(t[5], me5);
					t[0]=_mm256_add_epi16(t[0], _mm256_min_epi16(eN0, eNEE0));
					t[1]=_mm256_add_epi16(t[1], _mm256_min_epi16(eN1, eNEE1));
					t[2]=_mm256_add_epi16(t[2], _mm256_min_epi16(eN2, eNEE2));
					t[3]=_mm256_add_epi16(t[3], _mm256_min_epi16(eN3, eNEE3));
					t[4]=_mm256_add_epi16(t[4], _mm256_min_epi16(eN4, eNEE4));
					t[5]=_mm256_add_epi16(t[5], _mm256_min_epi16(eN5, eNEE5));
					t[0]=_mm256_srli_epi16(t[0], 2);
					t[1]=_mm256_srli_epi16(t[1], 2);
					t[2]=_mm256_srli_epi16(t[2], 2);
					t[3]=_mm256_srli_epi16(t[3], 2);
					t[4]=_mm256_srli_epi16(t[4], 2);
					t[5]=_mm256_srli_epi16(t[5], 2);
					_mm256_store_si256((__m256i*)(erows[0]+((kp+0*WG_NPREDS)*3+0)*NCODERS)+0, t[0]);//ecurr
					_mm256_store_si256((__m256i*)(erows[0]+((kp+0*WG_NPREDS)*3+0)*NCODERS)+1, t[1]);
					_mm256_store_si256((__m256i*)(erows[0]+((kp+0*WG_NPREDS)*3+1)*NCODERS)+0, t[2]);
					_mm256_store_si256((__m256i*)(erows[0]+((kp+0*WG_NPREDS)*3+1)*NCODERS)+1, t[3]);
					_mm256_store_si256((__m256i*)(erows[0]+((kp+0*WG_NPREDS)*3+2)*NCODERS)+0, t[4]);
					_mm256_store_si256((__m256i*)(erows[0]+((kp+0*WG_NPREDS)*3+2)*NCODERS)+1, t[5]);
#else
					//curr = eW + ((e-eW+(1<<3>>1))>>3)	X
					__m256i s0=_mm256_sub_epi16(me0, t[0]);
					__m256i s1=_mm256_sub_epi16(me1, t[1]);
					__m256i s2=_mm256_sub_epi16(me2, t[2]);
					__m256i s3=_mm256_sub_epi16(me3, t[3]);
					__m256i s4=_mm256_sub_epi16(me4, t[4]);
					__m256i s5=_mm256_sub_epi16(me5, t[5]);
					__m256i bias=_mm256_set1_epi16(1<<2>>1);
					s0=_mm256_add_epi16(s0, bias);
					s1=_mm256_add_epi16(s1, bias);
					s2=_mm256_add_epi16(s2, bias);
					s3=_mm256_add_epi16(s3, bias);
					s4=_mm256_add_epi16(s4, bias);
					s5=_mm256_add_epi16(s5, bias);
					s0=_mm256_srai_epi16(s0, 2);
					s1=_mm256_srai_epi16(s1, 2);
					s2=_mm256_srai_epi16(s2, 2);
					s3=_mm256_srai_epi16(s3, 2);
					s4=_mm256_srai_epi16(s4, 2);
					s5=_mm256_srai_epi16(s5, 2);
					s0=_mm256_add_epi16(s0, t[0]);
					s1=_mm256_add_epi16(s1, t[1]);
					s2=_mm256_add_epi16(s2, t[2]);
					s3=_mm256_add_epi16(s3, t[3]);
					s4=_mm256_add_epi16(s4, t[4]);
					s5=_mm256_add_epi16(s5, t[5]);
					_mm256_store_si256((__m256i*)(erows[0]+((kp+0*WG_NPREDS)*3+0)*NCODERS)+0, s0);//ecurr
					_mm256_store_si256((__m256i*)(erows[0]+((kp+0*WG_NPREDS)*3+0)*NCODERS)+1, s1);
					_mm256_store_si256((__m256i*)(erows[0]+((kp+0*WG_NPREDS)*3+1)*NCODERS)+0, s2);
					_mm256_store_si256((__m256i*)(erows[0]+((kp+0*WG_NPREDS)*3+1)*NCODERS)+1, s3);
					_mm256_store_si256((__m256i*)(erows[0]+((kp+0*WG_NPREDS)*3+2)*NCODERS)+0, s4);
					_mm256_store_si256((__m256i*)(erows[0]+((kp+0*WG_NPREDS)*3+2)*NCODERS)+1, s5);
#endif
				}
#endif
				for(int kr=0;kr<6;++kr)//reuse register but fragmented access
				{
					__m256i curr=W[kr];
					__m256i me0=_mm256_sub_epi16(curr, wgpreds[0*6+kr]);
					__m256i me1=_mm256_sub_epi16(curr, wgpreds[1*6+kr]);
					__m256i me2=_mm256_sub_epi16(curr, wgpreds[2*6+kr]);
					__m256i me3=_mm256_sub_epi16(curr, wgpreds[3*6+kr]);
					__m256i me4=_mm256_sub_epi16(curr, wgpreds[4*6+kr]);
					__m256i me5=_mm256_sub_epi16(curr, wgpreds[5*6+kr]);
					__m256i me6=_mm256_sub_epi16(curr, wgpreds[6*6+kr]);
					__m256i me7=_mm256_sub_epi16(curr, wgpreds[7*6+kr]);
					me0=_mm256_abs_epi16(me0);
					me1=_mm256_abs_epi16(me1);
					me2=_mm256_abs_epi16(me2);
					me3=_mm256_abs_epi16(me3);
					me4=_mm256_abs_epi16(me4);
					me5=_mm256_abs_epi16(me5);
					me6=_mm256_abs_epi16(me6);
					me7=_mm256_abs_epi16(me7);
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
					
					//                 (__m256i*)(erows[-Y]+((P+X*NPREDS)*3+C)*NCODERS)+R
					_mm256_store_si256((__m256i*)erows[0]+(0+0*WG_NPREDS)*6+kr, me0);//ecurr
					_mm256_store_si256((__m256i*)erows[0]+(1+0*WG_NPREDS)*6+kr, me1);
					_mm256_store_si256((__m256i*)erows[0]+(2+0*WG_NPREDS)*6+kr, me2);
					_mm256_store_si256((__m256i*)erows[0]+(3+0*WG_NPREDS)*6+kr, me3);
					_mm256_store_si256((__m256i*)erows[0]+(4+0*WG_NPREDS)*6+kr, me4);
					_mm256_store_si256((__m256i*)erows[0]+(5+0*WG_NPREDS)*6+kr, me5);
					_mm256_store_si256((__m256i*)erows[0]+(6+0*WG_NPREDS)*6+kr, me6);
					_mm256_store_si256((__m256i*)erows[0]+(7+0*WG_NPREDS)*6+kr, me7);
					//__m256i t[8];
					//t[0]=_mm256_add_epi16(me0, _mm256_load_si256((__m256i*)erows[1]+(0+1*WG_NPREDS)*6+kr));//eNE+=e
					//t[1]=_mm256_add_epi16(me1, _mm256_load_si256((__m256i*)erows[1]+(1+1*WG_NPREDS)*6+kr));
					//t[2]=_mm256_add_epi16(me2, _mm256_load_si256((__m256i*)erows[1]+(2+1*WG_NPREDS)*6+kr));
					//t[3]=_mm256_add_epi16(me3, _mm256_load_si256((__m256i*)erows[1]+(3+1*WG_NPREDS)*6+kr));
					//t[4]=_mm256_add_epi16(me4, _mm256_load_si256((__m256i*)erows[1]+(4+1*WG_NPREDS)*6+kr));
					//t[5]=_mm256_add_epi16(me5, _mm256_load_si256((__m256i*)erows[1]+(5+1*WG_NPREDS)*6+kr));
					//t[6]=_mm256_add_epi16(me6, _mm256_load_si256((__m256i*)erows[1]+(6+1*WG_NPREDS)*6+kr));
					//t[7]=_mm256_add_epi16(me7, _mm256_load_si256((__m256i*)erows[1]+(7+1*WG_NPREDS)*6+kr));
					//_mm256_store_si256((__m256i*)erows[1]+(0+1*WG_NPREDS)*6+kr, t[0]);
					//_mm256_store_si256((__m256i*)erows[1]+(1+1*WG_NPREDS)*6+kr, t[1]);
					//_mm256_store_si256((__m256i*)erows[1]+(2+1*WG_NPREDS)*6+kr, t[2]);
					//_mm256_store_si256((__m256i*)erows[1]+(3+1*WG_NPREDS)*6+kr, t[3]);
					//_mm256_store_si256((__m256i*)erows[1]+(4+1*WG_NPREDS)*6+kr, t[4]);
					//_mm256_store_si256((__m256i*)erows[1]+(5+1*WG_NPREDS)*6+kr, t[5]);
					//_mm256_store_si256((__m256i*)erows[1]+(6+1*WG_NPREDS)*6+kr, t[6]);
					//_mm256_store_si256((__m256i*)erows[1]+(7+1*WG_NPREDS)*6+kr, t[7]);

#ifdef WG_ENABLE_eW
					__m256i t2[8];
					t2[0]=_mm256_load_si256((__m256i*)wgWerrors+0*6+kr);
					t2[1]=_mm256_load_si256((__m256i*)wgWerrors+1*6+kr);
					t2[2]=_mm256_load_si256((__m256i*)wgWerrors+2*6+kr);
					t2[3]=_mm256_load_si256((__m256i*)wgWerrors+3*6+kr);
					t2[4]=_mm256_load_si256((__m256i*)wgWerrors+4*6+kr);
					t2[5]=_mm256_load_si256((__m256i*)wgWerrors+5*6+kr);
					t2[6]=_mm256_load_si256((__m256i*)wgWerrors+6*6+kr);
					t2[7]=_mm256_load_si256((__m256i*)wgWerrors+7*6+kr);
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
					_mm256_store_si256((__m256i*)wgWerrors+0*6+kr, t2[0]);
					_mm256_store_si256((__m256i*)wgWerrors+1*6+kr, t2[1]);
					_mm256_store_si256((__m256i*)wgWerrors+2*6+kr, t2[2]);
					_mm256_store_si256((__m256i*)wgWerrors+3*6+kr, t2[3]);
					_mm256_store_si256((__m256i*)wgWerrors+4*6+kr, t2[4]);
					_mm256_store_si256((__m256i*)wgWerrors+5*6+kr, t2[5]);
					_mm256_store_si256((__m256i*)wgWerrors+6*6+kr, t2[6]);
					_mm256_store_si256((__m256i*)wgWerrors+7*6+kr, t2[7]);
#endif
				}
				erows[0]+=NCODERS*3*WG_NPREDS;
				erows[1]+=NCODERS*3*WG_NPREDS;
			}
			//context update = (2*eW+(e<<3)+max(eNEE, eNEEE))>>2
			eNEEE[0]=_mm256_load_si256((__m256i*)(rows[1]+(3+0+3*6)*NCODERS)+0);
			eNEEE[1]=_mm256_load_si256((__m256i*)(rows[1]+(3+0+3*6)*NCODERS)+1);
			eNEEE[2]=_mm256_load_si256((__m256i*)(rows[1]+(3+1+3*6)*NCODERS)+0);
			eNEEE[3]=_mm256_load_si256((__m256i*)(rows[1]+(3+1+3*6)*NCODERS)+1);
			eNEEE[4]=_mm256_load_si256((__m256i*)(rows[1]+(3+2+3*6)*NCODERS)+0);
			eNEEE[5]=_mm256_load_si256((__m256i*)(rows[1]+(3+2+3*6)*NCODERS)+1);
			eW[0]=_mm256_slli_epi16(eW[0], 1);
			eW[1]=_mm256_slli_epi16(eW[1], 1);
			eW[2]=_mm256_slli_epi16(eW[2], 1);
			eW[3]=_mm256_slli_epi16(eW[3], 1);
			eW[4]=_mm256_slli_epi16(eW[4], 1);
			eW[5]=_mm256_slli_epi16(eW[5], 1);
			eW[0]=_mm256_add_epi16(eW[0], _mm256_max_epi16(eNEE[0], eNEEE[0]));
			eW[1]=_mm256_add_epi16(eW[1], _mm256_max_epi16(eNEE[1], eNEEE[1]));
			eW[2]=_mm256_add_epi16(eW[2], _mm256_max_epi16(eNEE[2], eNEEE[2]));
			eW[3]=_mm256_add_epi16(eW[3], _mm256_max_epi16(eNEE[3], eNEEE[3]));
			eW[4]=_mm256_add_epi16(eW[4], _mm256_max_epi16(eNEE[4], eNEEE[4]));
			eW[5]=_mm256_add_epi16(eW[5], _mm256_max_epi16(eNEE[5], eNEEE[5]));
			ecurr[0]=_mm256_slli_epi16(ecurr[0], GRBITS);
			ecurr[1]=_mm256_slli_epi16(ecurr[1], GRBITS);
			ecurr[2]=_mm256_slli_epi16(ecurr[2], GRBITS);
			ecurr[3]=_mm256_slli_epi16(ecurr[3], GRBITS);
			ecurr[4]=_mm256_slli_epi16(ecurr[4], GRBITS);
			ecurr[5]=_mm256_slli_epi16(ecurr[5], GRBITS);
			eW[0]=_mm256_add_epi16(eW[0], ecurr[0]);
			eW[1]=_mm256_add_epi16(eW[1], ecurr[1]);
			eW[2]=_mm256_add_epi16(eW[2], ecurr[2]);
			eW[3]=_mm256_add_epi16(eW[3], ecurr[3]);
			eW[4]=_mm256_add_epi16(eW[4], ecurr[4]);
			eW[5]=_mm256_add_epi16(eW[5], ecurr[5]);
			eW[0]=_mm256_srli_epi16(eW[0], 2);
			eW[1]=_mm256_srli_epi16(eW[1], 2);
			eW[2]=_mm256_srli_epi16(eW[2], 2);
			eW[3]=_mm256_srli_epi16(eW[3], 2);
			eW[4]=_mm256_srli_epi16(eW[4], 2);
			eW[5]=_mm256_srli_epi16(eW[5], 2);
			_mm256_store_si256((__m256i*)(rows[0]+(3+0+0*6)*NCODERS)+0, eW[0]);//store current contexts
			_mm256_store_si256((__m256i*)(rows[0]+(3+0+0*6)*NCODERS)+1, eW[1]);
			_mm256_store_si256((__m256i*)(rows[0]+(3+1+0*6)*NCODERS)+0, eW[2]);
			_mm256_store_si256((__m256i*)(rows[0]+(3+1+0*6)*NCODERS)+1, eW[3]);
			_mm256_store_si256((__m256i*)(rows[0]+(3+2+0*6)*NCODERS)+0, eW[4]);
			_mm256_store_si256((__m256i*)(rows[0]+(3+2+0*6)*NCODERS)+1, eW[5]);
			eNEE[0]=eNEEE[0];
			eNEE[1]=eNEEE[1];
			eNEE[2]=eNEEE[2];
			eNEE[3]=eNEEE[3];
			eNEE[4]=eNEEE[4];
			eNEE[5]=eNEEE[5];
			NW[0]=N[0];
			NW[1]=N[1];
			NW[2]=N[2];
			NW[3]=N[3];
			NW[4]=N[4];
			NW[5]=N[5];
			rows[0]+=6*NCODERS;
			rows[1]+=6*NCODERS;
			rows[2]+=6*NCODERS;
			rows[3]+=6*NCODERS;
			imptr+=3*NCODERS;
		}
#endif
	}
	prof_checkpoint(isize, "main");
	_mm_free(wgerrors);
	_mm_free(wgstate);
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
			profile_size(streamptr, "/ %9td bytes remainders _|rem %d %d", usize-isize, yremh, xremw);
		}

		//encode main
		mstate[0]=_mm256_set1_epi32(1<<(RANS_STATE_BITS-RANS_RENORM_BITS));
		mstate[1]=mstate[0];
		mstate[2]=mstate[0];
		mstate[3]=mstate[0];
		unsigned short *ctxptr2=(unsigned short*)(interleaved+(isize<<1)-sizeof(short[NCODERS]));
		for(int ky=blockh-1;ky>=0;--ky)
		{
#ifdef ESTIMATE_SIZE
			int kc=2;
#endif
			for(int kx=ixbytes-NCODERS;kx>=0;kx-=NCODERS)//ixbytes = iw/XCODERS*NCODERS*3
			{
				__m256i mmax[4], minvf[4], mcdf[4], mnegf_sh[4];
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

					s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[2*8+0])));
					s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[2*8+1])));
					s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[2*8+2])));
					s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[2*8+3])));
					s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[2*8+4])), 1);
					s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[2*8+5])), 1);
					s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[2*8+6])), 1);
					s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[2*8+7])), 1);
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[2]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[2]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[2]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[2]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));

					s0=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[3*8+0])));
					s1=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[3*8+1])));
					s2=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[3*8+2])));
					s3=_mm256_castsi128_si256(	_mm_load_si128((__m128i*)(syminfo+ctxptr2[3*8+3])));
					s0=_mm256_inserti128_si256(s0,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[3*8+4])), 1);
					s1=_mm256_inserti128_si256(s1,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[3*8+5])), 1);
					s2=_mm256_inserti128_si256(s2,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[3*8+6])), 1);
					s3=_mm256_inserti128_si256(s3,	_mm_load_si128((__m128i*)(syminfo+ctxptr2[3*8+7])), 1);
					t0=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(1, 0, 1, 0));
					t2=SHUFFLE_PS(s0, s1, _MM_SHUFFLE(3, 2, 3, 2));
					t1=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(1, 0, 1, 0));
					t3=SHUFFLE_PS(s2, s3, _MM_SHUFFLE(3, 2, 3, 2));
					mmax	[3]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
					minvf	[3]=SHUFFLE_PS(t0, t1, _MM_SHUFFLE(3, 1, 3, 1));
					mcdf	[3]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(2, 0, 2, 0));
					mnegf_sh[3]=SHUFFLE_PS(t2, t3, _MM_SHUFFLE(3, 1, 3, 1));

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
				__m256i cond2=_mm256_cmpgt_epi32(mstate[2], mmax[2]);
				__m256i cond3=_mm256_cmpgt_epi32(mstate[3], mmax[3]);
				int mask0=_mm256_movemask_ps(_mm256_castsi256_ps(cond0));
				int mask1=_mm256_movemask_ps(_mm256_castsi256_ps(cond1));
				int mask2=_mm256_movemask_ps(_mm256_castsi256_ps(cond2));
				int mask3=_mm256_movemask_ps(_mm256_castsi256_ps(cond3));
				__m256i idx0=_mm256_load_si256((__m256i*)ans_permute+mask0);
				__m256i idx1=_mm256_load_si256((__m256i*)ans_permute+mask1);
				__m256i idx2=_mm256_load_si256((__m256i*)ans_permute+mask2);
				__m256i idx3=_mm256_load_si256((__m256i*)ans_permute+mask3);
				__m256i emit0=_mm256_permutevar8x32_epi32(_mm256_and_si256(mstate[0], cond0), idx0);
				__m256i emit1=_mm256_permutevar8x32_epi32(_mm256_and_si256(mstate[1], cond1), idx1);
				__m256i emit2=_mm256_permutevar8x32_epi32(_mm256_and_si256(mstate[2], cond2), idx2);
				__m256i emit3=_mm256_permutevar8x32_epi32(_mm256_and_si256(mstate[3], cond3), idx3);
				emit0=_mm256_and_si256(emit0, _mm256_set1_epi32(0xFFFF));
				emit1=_mm256_and_si256(emit1, _mm256_set1_epi32(0xFFFF));
				emit2=_mm256_and_si256(emit2, _mm256_set1_epi32(0xFFFF));
				emit3=_mm256_and_si256(emit3, _mm256_set1_epi32(0xFFFF));
				emit0=_mm256_packus_epi32(emit0, emit1);
				emit2=_mm256_packus_epi32(emit2, emit3);
				emit0=_mm256_permute4x64_epi64(emit0, _MM_SHUFFLE(3, 1, 2, 0));
				emit2=_mm256_permute4x64_epi64(emit2, _MM_SHUFFLE(3, 1, 2, 0));
				__m128i e3=_mm256_extractf128_si256(emit2, 1);
				__m128i e2=_mm256_castsi256_si128(emit2);
				__m128i e1=_mm256_extractf128_si256(emit0, 1);
				__m128i e0=_mm256_castsi256_si128(emit0);
				mask3=_mm_popcnt_u32(mask3);
				mask2=_mm_popcnt_u32(mask2);
				mask1=_mm_popcnt_u32(mask1);
				mask0=_mm_popcnt_u32(mask0);
#ifdef _DEBUG
				if(streamptr-(2*((ptrdiff_t)mask0+mask1+mask2+mask3)+sizeof(__m128i))<=image)
					LOG_ERROR("OOB ptr %016zX <= %016zX", streamptr, image);
#endif
				_mm_storeu_si128((__m128i*)streamptr-1, e3); streamptr-=mask3*sizeof(short);
				_mm_storeu_si128((__m128i*)streamptr-1, e2); streamptr-=mask2*sizeof(short);
				_mm_storeu_si128((__m128i*)streamptr-1, e1); streamptr-=mask1*sizeof(short);
				_mm_storeu_si128((__m128i*)streamptr-1, e0); streamptr-=mask0*sizeof(short);
				{
					__m256i state0=_mm256_srli_epi32(mstate[0], 16);
					__m256i state1=_mm256_srli_epi32(mstate[1], 16);
					__m256i state2=_mm256_srli_epi32(mstate[2], 16);
					__m256i state3=_mm256_srli_epi32(mstate[3], 16);
					mstate[0]=_mm256_blendv_epi8(mstate[0], state0, cond0);
					mstate[1]=_mm256_blendv_epi8(mstate[1], state1, cond1);
					mstate[2]=_mm256_blendv_epi8(mstate[2], state2, cond2);
					mstate[3]=_mm256_blendv_epi8(mstate[3], state3, cond3);
				}
#ifdef ANS_VAL
				ansval_push(mstate, sizeof(int), NCODERS);
#endif
				//enc update		state += (state*invf>>32>>sh)*negf+cdf		state = state/freq<<12|(cdf+state%freq)
				{
					__m256i lo0=_mm256_mul_epu32(mstate[0], minvf[0]);//q = mulhi32(state, invf)
					__m256i lo1=_mm256_mul_epu32(mstate[1], minvf[1]);
					__m256i lo2=_mm256_mul_epu32(mstate[2], minvf[2]);
					__m256i lo3=_mm256_mul_epu32(mstate[3], minvf[3]);
					__m256i hi0=_mm256_mul_epu32(_mm256_srli_epi64(mstate[0], 32), _mm256_srli_epi64(minvf[0], 32));
					__m256i hi1=_mm256_mul_epu32(_mm256_srli_epi64(mstate[1], 32), _mm256_srli_epi64(minvf[1], 32));
					__m256i hi2=_mm256_mul_epu32(_mm256_srli_epi64(mstate[2], 32), _mm256_srli_epi64(minvf[2], 32));
					__m256i hi3=_mm256_mul_epu32(_mm256_srli_epi64(mstate[3], 32), _mm256_srli_epi64(minvf[3], 32));
					minvf[0]=_mm256_blend_epi32(_mm256_srli_epi64(lo0, 32), hi0, 0xAA);
					minvf[1]=_mm256_blend_epi32(_mm256_srli_epi64(lo1, 32), hi1, 0xAA);
					minvf[2]=_mm256_blend_epi32(_mm256_srli_epi64(lo2, 32), hi2, 0xAA);
					minvf[3]=_mm256_blend_epi32(_mm256_srli_epi64(lo3, 32), hi3, 0xAA);
				}
				{
					__m256i sh0=_mm256_srli_epi32(mnegf_sh[0], 16);
					__m256i sh1=_mm256_srli_epi32(mnegf_sh[1], 16);
					__m256i sh2=_mm256_srli_epi32(mnegf_sh[2], 16);
					__m256i sh3=_mm256_srli_epi32(mnegf_sh[3], 16);
					minvf[0]=_mm256_srlv_epi32(minvf[0], sh0);
					minvf[1]=_mm256_srlv_epi32(minvf[1], sh1);
					minvf[2]=_mm256_srlv_epi32(minvf[2], sh2);
					minvf[3]=_mm256_srlv_epi32(minvf[3], sh3);
				}
				mstate[0]=_mm256_add_epi32(mstate[0], mcdf[0]);
				mstate[1]=_mm256_add_epi32(mstate[1], mcdf[1]);
				mstate[2]=_mm256_add_epi32(mstate[2], mcdf[2]);
				mstate[3]=_mm256_add_epi32(mstate[3], mcdf[3]);
				{
					__m256i lomask=_mm256_set1_epi32(0xFFFF);
					__m256i negf0=_mm256_and_si256(mnegf_sh[0], lomask);
					__m256i negf1=_mm256_and_si256(mnegf_sh[1], lomask);
					__m256i negf2=_mm256_and_si256(mnegf_sh[2], lomask);
					__m256i negf3=_mm256_and_si256(mnegf_sh[3], lomask);
					minvf[0]=_mm256_mullo_epi32(minvf[0], negf0);
					minvf[1]=_mm256_mullo_epi32(minvf[1], negf1);
					minvf[2]=_mm256_mullo_epi32(minvf[2], negf2);
					minvf[3]=_mm256_mullo_epi32(minvf[3], negf3);
#ifdef ANS_VAL
					__m256i mM=_mm256_set1_epi32(1<<PROBBITS);
					negf0=_mm256_sub_epi16(mM, negf0);
					negf1=_mm256_sub_epi16(mM, negf1);
					negf2=_mm256_sub_epi16(mM, negf2);
					negf3=_mm256_sub_epi16(mM, negf3);
					negf0=_mm256_packus_epi32(negf0, negf1);
					negf2=_mm256_packus_epi32(negf2, negf3);
					negf0=_mm256_permute4x64_epi64(negf0, _MM_SHUFFLE(3, 1, 2, 0));
					negf2=_mm256_permute4x64_epi64(negf2, _MM_SHUFFLE(3, 1, 2, 0));
					ALIGN(32) unsigned short freqs[NCODERS];
					_mm256_store_si256((__m256i*)freqs+0, negf0);
					_mm256_store_si256((__m256i*)freqs+1, negf2);
					ansval_push(freqs, sizeof(short), NCODERS);
#endif
				}
				mstate[0]=_mm256_add_epi32(mstate[0], minvf[0]);
				mstate[1]=_mm256_add_epi32(mstate[1], minvf[1]);
				mstate[2]=_mm256_add_epi32(mstate[2], minvf[2]);
				mstate[3]=_mm256_add_epi32(mstate[3], minvf[3]);
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
			csize2+=fwrite("29", 1, 2, fdst);
			csize2+=fwrite(&iw, 1, 4, fdst);
			csize2+=fwrite(&ih, 1, 4, fdst);
			int flags=bestrct<<1|use_wg4;
			csize2+=fwrite(&flags, 1, 1, fdst);
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