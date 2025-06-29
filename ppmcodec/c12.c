#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<sys/stat.h>
#ifdef _MSC_VER
#include<intrin.h>
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>//QueryPerformanceCounter
#elif defined __GNUC__
#include<x86intrin.h>
#include<time.h>
#endif
#include<immintrin.h>


#ifdef _MSC_VER
	#define LOUD

//	#define PRINT_RCT
	#define ESTIMATE_SIZE
	#define ENABLE_GUIDE
//	#define PRINTBITS
//	#define PRINTGR
#endif

//	#define USE_NONLINEARITY


#define L1SH 18
#define PREDLIST\
	PRED(N)\
	PRED(W)\
	PRED(3*(W-WW)+WWW)\
	PRED(3*(N-NN)+NNN)\
	PRED(N+W-NW)\
	PRED(W+NE-N)\
	PRED(N+NE-NNE)\
	PRED((WWWW+WWW+NNN+NEE+NEEE+NEEEE-N-W)/4)\

enum
{
#define PRED(EXPR) +1
	L1NPREDS=PREDLIST,
#undef  PRED
};


#define GRBITS 6
#define GRLIMIT 24
#define PROBBITS_STORE 24
#define PROBBITS_USE 15


#ifdef _MSC_VER
#define INLINE __forceinline static
#else
#define INLINE __attribute__((always_inline)) inline static
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
#ifdef _MSC_VER
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
static void crash(const char *file, int line, const char *format, ...)
{
	printf("%s(%d):\n", file, line);
	if(format)
	{
		va_list args;
		printf("\n");
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
		printf("\n");
	}
	printf("CRASH\n");
	exit(1);
}
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT, ##__VA_ARGS__)
#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static double g_sqe=0;
static void guide_save(unsigned char *image, int iw, int ih)
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
static void guide_check(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
	{
		CRASH("Guide error  XY %d %d", kx, ky);
		printf("");
	}
}
static void guide_update(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx), diff;
	diff=g_image[idx+0]-image[idx+0]; g_sqe+=diff*diff;
	diff=g_image[idx+1]-image[idx+1]; g_sqe+=diff*diff;
	diff=g_image[idx+2]-image[idx+2]; g_sqe+=diff*diff;
}
#else
#define guide_save(...)
#define guide_check(...)
#define guide_update(...)
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
//	II_COEFF_U_SUB_V2,
//	II_COEFF_V_SUB_U2,

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

#ifdef ESTIMATE_SIZE
static int32_t hist[3][256]={0};
#endif
static uint32_t stats[3][256][8][GRLIMIT+8];
typedef struct _ACState
{
	uint64_t low, range, code;
	unsigned char *ptr, *end;
	uint32_t symidx, totalsyms;
} ACState;
#ifdef USE_NONLINEARITY
static int32_t squashtable[1<<PROBBITS_USE];
#endif
INLINE void codebit(ACState *ac, uint32_t *pp0, int32_t *bit, const int fwd)
{
	int rbit;
	uint64_t r2, mid;

	//uint32_t cell=*pp0;
	//int32_t p0=cell&((1<<PROBBITS_USE)-1);
	//int32_t ctr=cell>>PROBBITS_USE;

	int32_t p00=*pp0;
	int32_t p0=p00>>(PROBBITS_STORE-PROBBITS_USE);

	//int32_t p001=*pp01, p002=*pp02, p003=*pp03;
	//int32_t p0=(p002-((p001+p003-2*p002)>>7))>>(PROBBITS_STORE-PROBBITS_USE);
	//CLAMP2(p0, 1, (1<<PROBBITS_USE)-1);

#if 1
	if(ac->range<=0xFFFF)
	{
		if(ac->ptr>=ac->end)
		{
			CRASH("ERROR at %d/%d  inflation %8.4lf%%\n"
				, (int32_t)ac->symidx
				, (int32_t)ac->totalsyms
				, 100.*ac->totalsyms/ac->symidx
			);
			exit(1);
		}
		if(fwd)
			*(uint32_t*)ac->ptr=(uint32_t)(ac->low>>32);
		else
			ac->code=ac->code<<32|*(uint32_t*)ac->ptr;
		ac->ptr+=4;
		ac->low<<=32;
		ac->range=ac->range<<32|0xFFFFFFFF;
		if(ac->range>~ac->low)
			ac->range=~ac->low;
	}
#endif
#if 0
	int renorm=ac->range<=0xFFFF;
	uint32_t str=(uint32_t)(ac->low>>32);
	if(fwd)
	{
		if(renorm)
			*(uint32_t*)ac->ptr=str;
	}
	else
	{
		uint64_t code2=ac->code<<32|*(uint32_t*)ac->ptr;
		if(renorm)
			ac->code=code2;
	}
	ac->ptr+=4*renorm;
	ac->low<<=32*renorm;
	ac->range<<=32*renorm;
	ac->range|=0xFFFFFFFF&-renorm;
	if(ac->range>~ac->low)
		ac->range=~ac->low;
#endif
	r2=ac->range*((uint64_t)p0+(p0<(1<<PROBBITS_USE>>1)))>>PROBBITS_USE;
	mid=ac->low+r2;
	ac->range-=r2;
	--r2;
	if(fwd)
		rbit=*bit;
	else
		*bit=rbit=ac->code>=mid;
	if(rbit)
		ac->low=mid;
	else
		ac->range=r2;

#if 0
	++ctr;
	if(ctr>(1<<(32-PROBBITS_USE))-1)
		ctr=(1<<(32-PROBBITS_USE))-1;
	{
		static const int gaintable[]=
		{
			((1<<17)+ 1-1) /  1,
			((1<<17)+ 2-1) /  2,
			((1<<17)+ 3-1) /  3,
			((1<<17)+ 4-1) /  4,
			((1<<17)+ 5-1) /  5,
			((1<<17)+ 6-1) /  6,
			((1<<17)+ 7-1) /  7,
			((1<<17)+ 8-1) /  8,
			((1<<17)+ 9-1) /  9,
			((1<<17)+10-1) / 10,
			((1<<17)+11-1) / 11,
			((1<<17)+12-1) / 12,
			((1<<17)+13-1) / 13,
			((1<<17)+14-1) / 14,
			((1<<17)+15-1) / 15,
			((1<<17)+16-1) / 16,
			((1<<17)+17-1) / 17,
		};
		int gain=32-_lzcnt_u32(ctr);
		gain=gaintable[gain];
		p0+=(((int64_t)!*bit<<PROBBITS_USE)-p0)*gain>>(17+3);
	//	p0+=((!*bit<<PROBBITS_USE)-p0)/gain>>3;
	}
	*pp0=ctr<<PROBBITS_USE|p0;
#endif

	int32_t truth=!rbit<<PROBBITS_STORE;
	*pp0=p00+((truth-p00)>>6);

	//*pp02=p002+(((!*bit<<PROBBITS_STORE)-p002)>>6);
}
INLINE void mainloop(int iw, int ih, int bestrct, int dist, uint8_t *image, uint8_t *stream, ACState *ac, const int lossy, const int fwd)
{
	int
		yidx=rct_combinations[bestrct][II_PERM_Y],
		uidx=rct_combinations[bestrct][II_PERM_U],
		vidx=rct_combinations[bestrct][II_PERM_V],
		cu0=rct_combinations[bestrct][II_COEFF_U_SUB_Y],
		cv0=rct_combinations[bestrct][II_COEFF_V_SUB_Y],
		cv1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
	int32_t ky, kx, idx;
	int32_t psize=0;
	int16_t *pixels=0;
	int32_t paddedwidth=iw+16;
	int32_t coeffs[3][L1NPREDS+1]={0};
	int32_t invdist=((1<<16)+dist-1)/dist;
	uint8_t *imptr=image;

	psize=(int32_t)sizeof(int16_t[4*3*2])*(iw+16);//4 padded rows * 3 channels * {pixels, nbypass}
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error\n");
		free(image);
		free(stream);
		return;
	}
	memset(pixels, 0, psize);
//	FILLMEM((uint32_t*)stats, 0xFFFFFFFF<<PROBBITS_USE|1<<PROBBITS_USE>>1, sizeof(stats), sizeof(int32_t));
//	FILLMEM((uint32_t*)stats, 1<<PROBBITS_USE>>1, sizeof(stats), sizeof(int32_t));
	FILLMEM((uint32_t*)stats, 1<<PROBBITS_STORE>>1, sizeof(stats), sizeof(int32_t));
//	memset(stats, 0, sizeof(stats));
	FILLMEM((int32_t*)coeffs, (1<<L1SH)/L1NPREDS, sizeof(coeffs), sizeof(int32_t));
	for(ky=0, idx=0;ky<ih;++ky)
	{
		char yuv[4]={0};
		int16_t *rows[]=
		{
			pixels+(paddedwidth*((ky-0)&3)+8)*3*2,
			pixels+(paddedwidth*((ky-1)&3)+8)*3*2,
			pixels+(paddedwidth*((ky-2)&3)+8)*3*2,
			pixels+(paddedwidth*((ky-3)&3)+8)*3*2,
		};
		for(kx=0;kx<iw;++kx, ++idx)
		{
			int kc;
			int offset;

			if(fwd)
			{
				yuv[0]=imptr[yidx]-128;
				yuv[1]=imptr[uidx]-128;
				yuv[2]=imptr[vidx]-128;
			}
			offset=0;
			for(kc=0;kc<3;++kc)
			{
				int32_t
					NNN	=rows[3][0+0*3*2],
					NNE	=rows[2][0+1*3*2],
					NN	=rows[2][0+0*3*2],
					NW	=rows[1][0-1*3*2],
					N	=rows[1][0+0*3*2],
					NE	=rows[1][0+1*3*2],
					NEE	=rows[1][0+2*3*2],
					NEEE	=rows[1][0+3*3*2],
					NEEEE	=rows[1][0+4*3*2],
					WWWW	=rows[0][0-4*3*2],
					WWW	=rows[0][0-3*3*2],
					WW	=rows[0][0-2*3*2],
					W	=rows[0][0-1*3*2],
					eNEE	=rows[1][1+2*3*2],
					eNEEE	=rows[1][1+3*3*2],
					eW	=rows[0][1-1*3*2];
				int32_t pred, vmax, vmin, pred0;
				int32_t error;
				int32_t nbypass, nzeros=-1, bypass=0;
				int32_t tidx=0;
				uint32_t *statsptr;
				int32_t bit;

				int32_t preds[]=
				{
#define PRED(EXPR) EXPR,
					PREDLIST
#undef  PRED
				};
				pred=coeffs[kc][L1NPREDS];
				{
					int j=0;
					for(;j<L1NPREDS;++j)
						pred+=coeffs[kc][j]*preds[j];
				}
				pred+=1<<L1SH>>1;
				pred>>=L1SH;
				pred0=pred;
				vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);

				pred+=offset;
				CLAMP2(pred, -128, 127);

				//if(ky==10&&kx==2019)//
				//	printf("");

				nbypass=31-GRBITS-_lzcnt_u32(eW+1);
				CLAMP2(nbypass, 0, 7);
				//{
				//	float fval=(float)(eW+1);
				//	size_t addr=(size_t)&fval;
				//	int32_t bits=*(int32_t*)addr;
				//	nbypass=(bits>>23)-127-GRBITS;
				//	if(nbypass<0)
				//		nbypass=0;
				//}
				statsptr=stats[kc][(pred+128)&255][nbypass];
				//{
				//	int ctx=(pred+128)&255;
				//	statsptr1=stats[kc][ctx+0][nbypass];
				//	statsptr2=stats[kc][ctx+1][nbypass];
				//	statsptr3=stats[kc][ctx+2][nbypass];
				//}
				int upred=128-abs(pred);
					
				//if(ky==193&&kx==975&&!kc)//
				//if(ky==415&&kx==996&&!kc)//
				//	printf("");

				if(fwd)
				{
					if(lossy)
					{
						int e2=yuv[kc]-pred;
						if(e2<0)
							e2+=dist-1;
						e2=e2*invdist>>16;
						error=e2<<1^e2>>31;
						yuv[kc]=e2*dist+pred;
						CLAMP2(yuv[kc], -128, 127);
					}
					else
					{
						error=yuv[kc]-pred;
						{
							int negmask=error>>31;
							int abserr=(error^negmask)-negmask;
							error=error<<1^negmask;
							if(upred<abserr)
								error=upred+abserr;
							if(error==256)
							{
								error=(char)(yuv[kc]-pred);
								error=error<<1^error>>31;
							}
						}
#if 0
						if(pred>0&&yuv[kc]==-128)
						{
							error=(char)(yuv[kc]-pred);
							error=error<<1^error>>31;
						}
						else
						{
							error=yuv[kc]-pred;
							{
								int negmask=error>>31;
								int abserr=(error^negmask)-negmask;
								error=error<<1^negmask;
								if(upred<abserr)
									error=upred+abserr;
							}
						}
#endif
						//error=(char)(yuv[kc]-pred);
						//error=error<<1^error>>31;
					}
					nzeros=error>>nbypass;
					bypass=error&((1<<nbypass)-1);
#ifdef ESTIMATE_SIZE
					++hist[kc][error];
#endif
#ifdef PRINTGR
					{
						float fval=(float)(error+1);
						size_t addr=(size_t)&fval;
						int32_t bits=*(int32_t*)addr;
						bits=(bits>>23)-127;
						gr_bypsum+=nbypass;
						gr_symlen+=bits;
					}
#endif
				}
				else
					error=0;
				do
				{
					bit=nzeros--<=0;
					codebit(ac, statsptr+tidx, &bit, fwd);
#ifdef PRINTBITS
					if(fwd&&(unsigned)(idx-(usize>>2))<1000)printf("%c", '0'+bit);//
#endif
					++tidx;
					if(tidx==GRLIMIT)
					{
						tidx=1;
						nbypass=8;
						if(fwd)
							bypass=error;
						break;
					}
				}while(!bit);
				{
					int32_t kb=nbypass-1;

					for(;kb>=0;--kb)
					{
						bit=bypass>>kb&1;
						codebit(ac, statsptr+GRLIMIT+8-nbypass+kb, &bit, fwd);
						bypass|=bit<<kb;
#ifdef PRINTBITS
						if(fwd&&(unsigned)(idx-(usize>>2))<1000)printf("%c", '0'+bit);//
#endif
					}
				}
				if(!fwd)
				{
					error=(tidx-1)<<nbypass|bypass;
					if(lossy)
					{
						error=error>>1^-(error&1);
						yuv[kc]=error*dist+pred;
						CLAMP2(yuv[kc], -128, 127);
					}
					else
					{
						if(2*pred+error==256)
						{
							error=error>>1^-(error&1);
							yuv[kc]=(char)(error+pred);
						}
						else
						{
							int negmask=pred>>31;
							int sym=error;
							int e2=upred-sym;
							error=sym>>1^-(sym&1);
							e2=(e2^negmask)-negmask;
							if((upred<<1)<sym)
								error=e2;
							yuv[kc]=error+pred;
						}
						//error=error>>1^-(error&1);
						//yuv[kc]=(char)(error+pred);
					}
				}

				{
					int32_t curr=yuv[kc]-offset;
					int32_t k, e;

					e=(curr>pred0)-(curr<pred0);
					coeffs[kc][L1NPREDS]+=e;
					for(k=0;k<L1NPREDS;++k)
						coeffs[kc][k]+=e*preds[k];

					error=yuv[kc]-pred;
					error=error<<1^error>>31;
					rows[0][0]=curr;
					rows[0][1]=(4*eW+2*(error<<GRBITS)+eNEE+eNEEE)>>3;
				//	rows[0][1]=(2*eW+(error<<GRBITS)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
				}
				offset=(kc ? cv0*yuv[0]+cv1*yuv[1] : cu0*yuv[0])>>2;
				//offset=kc ? yuv[vhelpidx] : yuv[uhelpidx];
				rows[0]+=2;
				rows[1]+=2;
				rows[2]+=2;
				rows[3]+=2;
			}
			if(!fwd)
			{
				imptr[yidx]=yuv[0]+128;
				imptr[uidx]=yuv[1]+128;
				imptr[vidx]=yuv[2]+128;
#ifdef ENABLE_GUIDE
				if(lossy)
					guide_update(image, kx, ky);
				else
					guide_check(image, kx, ky);
#endif
			}
			imptr+=3;
		}
	}
	free(pixels);
}
int c12_codec(int argc, char **argv)
{
#if 0
	{
		printf("pred\\target  naive modular arithmetic sign packing\n");
		printf("\t");
		for(int k=0;k<8;++k)
			printf(" %3d", k-4);
		printf("\n\n");
		for(int kp=-4;kp<4;++kp)
		{
			printf(" %3d\t", kp);
			for(int kt=-4;kt<4;++kt)
			{
				int e=kt-kp, e0;

				e0=e<<(32-3)>>(32-3);
				e0=e0<<1^e0>>31;
				printf(" %3d", e0);
			}
			printf("\n");
		}
		printf("\n");

		printf("pred\\target  CALIC sign deduction\n");
		printf("\t");
		for(int k=0;k<8;++k)
			printf(" %3d", k-4);
		printf("\n\n");
		for(int kp=-4;kp<4;++kp)
		{
			printf(" %3d\t", kp);
			for(int kt=-4;kt<4;++kt)
			{
				int e=kt-kp, e1;

				if(kt==-4&&kp>0)
				{
					e1=e<<(32-3)>>(32-3);
					e1=e1<<1^e1>>31;
				}
				else
				{
					int upred=4-abs(kp);
					int negmask=e>>31;
					int abse=(e^negmask)-negmask;
					e1=e<<1^negmask;
					if(upred<abse)
						e1=upred+abse;
				}
				printf(" %3d", e1);

				//deduce kt from e1 and kp
				{
					int kt2=e1>>1^-(e1&1);
					kt2+=kp;
					kt2=kt2<<(32-3)>>(32-3);
					if(!(kp>0&&kt2==-4))
					{
						int upred=4-abs(kp);
						int negmask=kp>>31;
						int e2=upred-e1;
						kt2=e1>>1^-(e1&1);
						e2=(e2^negmask)-negmask;
						if(2*upred<e1)
							kt2=e2;
						kt2+=kp;
					}
					if(kt2!=kt)
						CRASH("ERROR");
				}
			}
			printf("\n");
		}
		printf("\n");
		exit(0);
	}
#endif
	if(argc!=3&&argc!=4)
	{
		printf(
			"Usage:  \"%s\"  input  output  [dist]\n"
			"  dist=1 for lossless (default).  Or 4 <= dist <= 16 for lossy.\n"
			, argv[0]
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];
	int dist=argc<4?1:atoi(argv[3]);
	if(dist!=1)
		CLAMP2(dist, 4, 16);
	ptrdiff_t srcsize=0, dstsize=0;
	int fwd=0;
	int32_t iw=0, ih=0;
	ptrdiff_t res=0, usize=0, csize=0;
	unsigned char *image=0, *imptr=0, *imend=0, *stream=0, *streamptr=0, *streamend=0;
	ptrdiff_t streamsize=0;
	int bestrct=0;
	ACState ac=
	{
		0, 0xFFFFFFFFFFFF, 0,
		0, 0,
	};
#ifdef LOUD
	double t=time_sec();
#endif
#ifdef PRINTGR
	long long gr_symlen=0, gr_bypsum=0;
#endif
	
	//read source
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
		if(!fwd&&c!=('1'|'2'<<8))
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
			bestrct=0;
			fread(&bestrct, 1, 1, fsrc);
			dist=1;
			fread(&dist, 1, 1, fsrc);

			nread=ftell(fsrc);
			streamsize=srcsize-nread;

		}
		if(iw<1||ih<1)
		{
			CRASH("Unsupported image dimensions  WH %d*%d", iw, ih);
			return 1;
		}
		res=(ptrdiff_t)iw*ih;
		usize=3*res;
		if(fwd)
			streamsize=usize;
		image=(unsigned char*)malloc(usize);
		stream=(unsigned char*)malloc(streamsize+sizeof(char[32]));
		if(!image||!stream)
		{
			CRASH("Alloc error");
			return 1;
		}
		imend=image+usize;
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
				expected=streamsize;
				nread=fread(stream, 1, streamsize, fsrc);
			}
			if(nread!=expected)
				printf("Truncated  expected %td  read %td", expected, nread);
		}
		fclose(fsrc);
	}
	if(fwd)
	{
		//analysis
		int64_t counters[OCH_COUNT]={0}, minerr=0;
		int prev[OCH_COUNT]={0};

		imptr=image;
		while(imptr<imend)
		{
			int r, g, b, rg, gb, br;

			r=imptr[0]<<2;
			g=imptr[1]<<2;
			b=imptr[2]<<2;
			imptr+=3;
			rg=r-g;
			gb=g-b;
			br=b-r;
#define UPDATE(I0, E0, I1, E1, I2, E2)\
	do\
	{\
		int t0=E0;\
		int t1=E1;\
		int t2=E2;\
		counters[I0]+=abs(t0-prev[I0]);\
		counters[I1]+=abs(t1-prev[I1]);\
		counters[I2]+=abs(t2-prev[I2]);\
		prev[I0]=t0;\
		prev[I1]=t1;\
		prev[I2]=t2;\
	}while(0)

			UPDATE(
				OCH_Y400, r,
				OCH_Y040, g,
				OCH_Y004, b
			);
			UPDATE(
				OCH_CX40, rg,
				OCH_C0X4, gb,
				OCH_C40X, br
			);
			UPDATE(
				OCH_CX31, rg+(gb>>2),//r-(3*g+b)/4 = r-g-(b-g)/4
				OCH_C3X1, rg+(br>>2),//g-(3*r+b)/4 = g-r-(b-r)/4
				OCH_C31X, br+(rg>>2) //b-(3*r+g)/4 = b-r-(g-r)/4
			);
			UPDATE(
				OCH_CX13, br+(gb>>2),//r-(g+3*b)/4 = r-b-(g-b)/4
				OCH_C1X3, gb+(br>>2),//g-(r+3*b)/4 = g-b-(r-b)/4
				OCH_C13X, gb+(rg>>2) //b-(r+3*g)/4 = b-g-(r-g)/4
			);
			UPDATE(
				OCH_CX22, (rg-br)>>1,//r-(g+b)/2 = (r-g + r-b)/2
				OCH_C2X2, (rg-gb)>>1,//g-(r+b)/2 = (g-r + g-b)/2
				OCH_C22X, (br-gb)>>1 //b-(r+g)/2 = (b-r + b-g)/2
			);
#undef  UPDATE
			//counters[0]+=abs(r	-prev[0]);
			//counters[1]+=abs(g	-prev[1]);
			//counters[2]+=abs(b	-prev[2]);
			//counters[3]+=abs(rg	-prev[3]);
			//counters[4]+=abs(gb	-prev[4]);
			//counters[5]+=abs(br	-prev[5]);
			//prev[0]=r;
			//prev[1]=g;
			//prev[2]=b;
			//prev[3]=rg;
			//prev[4]=gb;
			//prev[5]=br;
		}
		imptr=image;
		{
			int kt;

#ifdef PRINT_RCT
			for(kt=0;kt<OCH_COUNT;++kt)
				printf("%d %16lld\n", kt, counters[kt]);
			printf("\n");
#endif
			for(kt=0;kt<RCT_COUNT;++kt)
			{
				const unsigned char *combination=rct_combinations[kt];
				long long currerr=
					+counters[combination[0]]
					+counters[combination[1]]
					+counters[combination[2]]
				;
				if(!kt||minerr>currerr)
				{
					minerr=currerr;
					bestrct=kt;
				}
#ifdef PRINT_RCT
				printf("RCT%02d %16lld%s\n", kt, currerr, kt==bestrct?" <-":"");
#endif
			}
#ifdef LOUD
			printf("\"%s\"  WH %d*%d  %lld bytes  RCT %d %s\n"
				, srcfn
				, iw, ih, usize
				, bestrct, rct_names[bestrct]
			);
#else
			(void)och_names;
			(void)rct_names;
#endif
		}
		streamptr=stream;
		streamend=stream+usize;

		ac.ptr=streamptr;
		ac.end=streamend;
	}
	else
	{
		imptr=image;
		streamptr=stream;
		streamend=stream+srcsize;

		ac.code=*(uint32_t*)streamptr;//load
		streamptr+=4;
		ac.code=ac.code<<32|*(uint32_t*)streamptr;
		streamptr+=4;
		ac.ptr=streamptr;
		ac.end=streamend;

		csize=srcsize;
	}
	if(dist>1)
	{
		if(fwd)
			mainloop(iw, ih, bestrct, dist, image, stream, &ac, 1, 1);
		else
			mainloop(iw, ih, bestrct, dist, image, stream, &ac, 1, 0);
	}
	else
	{
		if(fwd)
			mainloop(iw, ih, bestrct, dist, image, stream, &ac, 0, 1);
		else
			mainloop(iw, ih, bestrct, dist, image, stream, &ac, 0, 0);
	}
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing\n", dstfn);
			free(image);
			free(stream);
			return 1;
		}
		if(fwd)
		{
			*(uint32_t*)ac.ptr=(uint32_t)(ac.low>>32);//flush
			ac.ptr+=4;
			*(uint32_t*)ac.ptr=(uint32_t)ac.low;
			ac.ptr+=4;

			csize=ac.ptr-stream;

			dstsize+=fwrite("12", 1, 2, fdst);
			dstsize+=fwrite(&iw, 1, 4, fdst);
			dstsize+=fwrite(&ih, 1, 4, fdst);
			dstsize+=fwrite(&bestrct, 1, 1, fdst);
			dstsize+=fwrite(&dist, 1, 1, fdst);
			dstsize+=fwrite(stream, 1, csize, fdst);
			csize=dstsize;
		}
		else
		{
			dstsize+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			dstsize+=fwrite(image, 1, usize, fdst);
			//usize=dstsize;
		}
		fclose(fdst);
	}
	free(image);
	free(stream);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
#ifdef PRINTBITS
		printf("\n");
#endif
#ifdef ESTIMATE_SIZE
		{
			double csizes[3]={0};
			int kc, ks;
			for(kc=0;kc<3;++kc)
			{
				double norm;
				int32_t sum=0;
				for(ks=0;ks<256;++ks)
					sum+=hist[kc][ks];
				norm=1./sum;
				for(ks=0;ks<256;++ks)
				{
					int32_t freq=hist[kc][ks];
					if(freq)
						csizes[kc]-=freq*log(freq*norm);
				}
				csizes[kc]*=1./(M_LN2*8);//convert ln->log2
			}
			printf("TYUV %12.2lf %12.2lf %12.2lf %12.2lf\n"
				, csizes[0]+csizes[1]+csizes[2]
				, csizes[0]
				, csizes[1]
				, csizes[2]
			);
		}
#endif
#ifdef PRINTGR
		printf("%12lld bypass %12.6lf bits\n", gr_bypsum, (double)gr_bypsum/usize);
		printf("%12lld symlen %12.6lf bits\n", gr_symlen, (double)gr_symlen/usize);
#endif
		printf("%9td->%9td  %8.4lf%%  %12.6lf\n"
			, usize
			, csize
			, 100.*csize/usize
			, (double)usize/csize
		);
	}
	printf("%c  %12.6lf sec  %12.6lf MB/s\n"
		, 'D'+fwd
		, t
		, usize/(t*1024*1024)
	);
#ifdef ENABLE_GUIDE
	if(!fwd&&dist>1)
	{
		double rmse=sqrt(g_sqe/((double)3*iw*ih)), psnr=20*log10(255/rmse);
		printf("RMSE %12.6lf PSNR %12.6lf\n", rmse, psnr);
	}
#endif
#endif
	(void)time_sec;
	return 0;
}
