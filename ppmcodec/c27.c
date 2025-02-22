#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define BYPASS_SSE
//	#define DISABLE_RCT
//	#define ENTROPY0
	#define CALICSIGNPRED
#ifdef _MSC_VER
	#define LOUD
//	#define ENABLE_MEMCHECKS
//	#define ENABLE_GUIDE
#endif

#define GRBITS 6
#define ANALYSIS_STRIDE 4

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

#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
	OCH(Y310) OCH(Y031) OCH(Y103)\
	OCH(Y301) OCH(Y130) OCH(Y013)\
	OCH(Y211) OCH(Y121) OCH(Y112)\
	OCH(CX40) OCH(C0X4) OCH(C40X)\
	OCH(CX31) OCH(C3X1) OCH(C31X)\
	OCH(CX13) OCH(C1X3) OCH(C13X)\
	OCH(CX22) OCH(C2X2) OCH(C22X)
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

	II_COEFF_Y_SUB_U,
	II_COEFF_Y_SUB_V,
	II_COEFF_U_SUB_V_NBLI,
	II_COEFF_V_SUB_U_NBLI,

	II_COUNT,
} RCTInfoIdx;
//YUV = RCT * RGB	watch out for permutation in last row
//luma: averaging	chroma: subtraction
//example: _400_40X_3X1 == [1 0 0; -1 0 1; -3/4 1 -1/4]
#define RCTLIST\
	RCT(_400_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0,	0, 0, 0, 0)\
	RCT(_400_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4,	0, 0, 0, 0)\
	RCT(_400_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0,	0, 0, 0, 0)\
	RCT(_040_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0,	0, 0, 0, 0)\
	RCT(_040_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4,	0, 0, 0, 0)\
	RCT(_004_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4,	0, 0, 0, 0)\
	RCT(_004_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0,	0, 0, 0, 0)\
	RCT(_040_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0,	0, 0, 0, 0)\
	RCT(_040_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4,	0, 0, 0, 0)\
	RCT(_040_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4,	0, 0, 0, 0)\
	RCT(_004_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0,	0, 0, 0, 0)\
	RCT(_004_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4,	0, 0, 0, 0)\
	RCT(_004_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4,	0, 0, 0, 0)\
	RCT(_400_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0,	0, 0, 0, 0)\
	RCT(_400_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4,	0, 0, 0, 0)\
	RCT(_400_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4,	0, 0, 0, 0)\
	RCT(_400_0X0_13X,	OCH_R,		OCH_G,		OCH_B1,		0, 1, 2,	0,  1, 3,	0, 0, 0, 0)\
	RCT(_400_4X0_13X,	OCH_R,		OCH_GR,		OCH_B1,		0, 1, 2,	4,  1, 3,	0, 0, 0, 0)\
	RCT(_400_00X_3X1,	OCH_R,		OCH_B,		OCH_G1,		0, 2, 1,	0,  3, 1,	0, 0, 0, 0)\
	RCT(_400_40X_3X1,	OCH_R,		OCH_BR,		OCH_G1,		0, 2, 1,	4,  3, 1,	0, 0, 0, 0)\
	RCT(_040_00X_X13,	OCH_G,		OCH_B,		OCH_R1,		1, 2, 0,	0,  1, 3,	0, 0, 0, 0)\
	RCT(_040_04X_X13,	OCH_G,		OCH_BG,		OCH_R1,		1, 2, 0,	4,  1, 3,	0, 0, 0, 0)\
	RCT(_040_X40_13X,	OCH_G,		OCH_RG,		OCH_B1,		1, 0, 2,	4,  3, 1,	0, 0, 0, 0)\
	RCT(_004_X04_3X1,	OCH_B,		OCH_RB,		OCH_G1,		2, 0, 1,	4,  1, 3,	0, 0, 0, 0)\
	RCT(_004_04X_X13,	OCH_B,		OCH_GB,		OCH_R1,		2, 1, 0,	4,  3, 1,	0, 0, 0, 0)\
	RCT(_400_0X0_22X,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,  2, 2,	0, 0, 0, 0)\
	RCT(_400_4X0_22X,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	4,  2, 2,	0, 0, 0, 0)\
	RCT(_400_00X_2X2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,  2, 2,	0, 0, 0, 0)\
	RCT(_400_40X_2X2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	4,  2, 2,	0, 0, 0, 0)\
	RCT(_040_00X_X22,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,  2, 2,	0, 0, 0, 0)\
	RCT(_040_04X_X22,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	4,  2, 2,	0, 0, 0, 0)\
	RCT(_040_X40_22X,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	4,  2, 2,	0, 0, 0, 0)\
	RCT(_004_X04_2X2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	4,  2, 2,	0, 0, 0, 0)\
	RCT(_004_0X4_X22,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	4,  2, 2,	0, 0, 0, 0)\
	RCT(_400_0X0_31X,	OCH_R,		OCH_G,		OCH_B3,		0, 1, 2,	0,  3, 1,	0, 0, 0, 0)\
	RCT(_400_4X0_31X,	OCH_R,		OCH_GR,		OCH_B3,		0, 1, 2,	4,  3, 1,	0, 0, 0, 0)\
	RCT(_400_00X_1X3,	OCH_R,		OCH_B,		OCH_G3,		0, 2, 1,	0,  1, 3,	0, 0, 0, 0)\
	RCT(_400_40X_1X3,	OCH_R,		OCH_BR,		OCH_G3,		0, 2, 1,	4,  1, 3,	0, 0, 0, 0)\
	RCT(_040_00X_X31,	OCH_G,		OCH_B,		OCH_R3,		1, 2, 0,	0,  3, 1,	0, 0, 0, 0)\
	RCT(_040_04X_X31,	OCH_G,		OCH_BG,		OCH_R3,		1, 2, 0,	4,  3, 1,	0, 0, 0, 0)\
	RCT(_040_X40_31X,	OCH_G,		OCH_RG,		OCH_B3,		1, 0, 2,	4,  1, 3,	0, 0, 0, 0)\
	RCT(_004_X04_1X3,	OCH_B,		OCH_RB,		OCH_G3,		2, 0, 1,	4,  3, 1,	0, 0, 0, 0)\
	RCT(_004_0X4_X31,	OCH_B,		OCH_GB,		OCH_R3,		2, 1, 0,	4,  1, 3,	0, 0, 0, 0)\
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
#ifdef BYPASS_SSE
ALIGN(16) static int bypass_sse[64][4]={0};
#endif
#ifdef ENTROPY0
int hist[OCH_COUNT][256]={0};
#endif
int c27_codec(const char *srcfn, const char *dstfn, int nthreads0)
{
#ifdef LOUD
	double t=time_sec();
	long long bitcount[3][3]={0};
#endif
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	ptrdiff_t srcsize;
	unsigned char *srcbuf=0, *srcptr=0, *srcend=0;
	{
		srcsize=get_filesize(srcfn);
		if(srcsize<1)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		srcbuf=(unsigned char*)malloc(srcsize+sizeof(__m256i[2]));
		if(!srcbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		fread(srcbuf, 1, srcsize, fsrc);
		fclose(fsrc);
		srcbuf[srcsize]=0;
	}
	srcptr=srcbuf;
	srcend=srcbuf+srcsize;
	int tag=*(unsigned short*)srcptr;
	srcptr+=2;
	int fwd=tag==('P'|'6'<<8);
	if(!fwd&&tag!=('2'|'7'<<8))
	{
		LOG_ERROR("Unsupported source file");
		return 1;
	}
	int iw=0, ih=0;
	if(fwd)//encode
	{
		if(*srcptr++ != '\n')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		while((unsigned)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++ - '0';
		while(*srcptr==' ')
			++srcptr;
		while((unsigned)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++ - '0';
		if(memcmp(srcptr, "\n255\n", 5))
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		srcptr+=5;
	}
	else//decode
	{
		if(srcptr+4*3>=srcend)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		memcpy(&iw, srcptr, 4); srcptr+=4;
		memcpy(&ih, srcptr, 4); srcptr+=4;
	}
	if(iw<1||ih<1)
	{
		LOG_ERROR("Unsupported source file");
		return 1;
	}
	ptrdiff_t usize=(ptrdiff_t)3*iw*ih;

	int nbits=fwd?64:0;
	unsigned long long cache=0;
	
	int bestrct=0;
	ptrdiff_t dstbufsize=0;
	unsigned char *dstbuf=0, *image=0, *streamptr=0;
#ifdef _MSC_VER
	unsigned char *streamend=0;
#endif
	if(fwd)
	{
		dstbufsize=(ptrdiff_t)4*iw*ih;
		dstbuf=malloc(dstbufsize);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		streamptr=dstbuf;
#ifdef _MSC_VER
		streamend=dstbuf+dstbufsize;
#endif
		image=srcptr;
		guide_save(image, iw, ih);
		
		int prev[OCH_COUNT]={0};
#ifdef ENTROPY0
		memset(hist, 0, sizeof(hist));
#else
		long long counters[OCH_COUNT]={0};
#endif
		while(srcptr<srcend)
		{
			int
#ifdef ENTROPY0
				r=srcptr[0]-128,
				g=srcptr[1]-128,
				b=srcptr[2]-128,
#else
				r=(srcptr[0]-128)<<2,
				g=(srcptr[1]-128)<<2,
				b=(srcptr[2]-128)<<2,
#endif
				rg=r-g,
				gb=g-b,
				br=b-r;
			int temp;
			srcptr+=3*ANALYSIS_STRIDE;
#ifdef ENTROPY0
#define UPDATE(IDX, EXPR) temp=EXPR, ++hist[IDX][(temp-prev[IDX])&255], prev[IDX]=temp
#else
#define UPDATE(IDX, EXPR) temp=EXPR, counters[IDX]+=abs(temp-prev[IDX]), prev[IDX]=temp
#endif
			UPDATE(OCH_Y400, r);
			UPDATE(OCH_Y040, g);
			UPDATE(OCH_Y004, b);
			UPDATE(OCH_Y301, r+(br>>2));
			UPDATE(OCH_Y031, g-(gb>>2));
			UPDATE(OCH_Y013, b+(gb>>2));
			UPDATE(OCH_Y310, r-(rg>>2));
			UPDATE(OCH_Y130, g+(rg>>2));
			UPDATE(OCH_Y103, b-(br>>2));
			UPDATE(OCH_Y211, r+((br-rg)>>2));	//(2*r+g+b)/4 = r + (g-r + b-r)/4
			UPDATE(OCH_Y121, g+((rg-gb)>>2));	//(r+2*g+b)/4 = g + (r-g + b-g)/4
			UPDATE(OCH_Y112, b+((gb-br)>>2));	//(r+g+2*b)/4 = b + (r-b + g-b)/4
			UPDATE(OCH_CX40, rg);
			UPDATE(OCH_C40X, br);
			UPDATE(OCH_C0X4, gb);
			UPDATE(OCH_CX31, rg+(gb>>2));	//r-(3*g+b)/4 = r-g-(b-g)/4
			UPDATE(OCH_C3X1, -rg-(br>>2));	//g-(3*r+b)/4 = g-r-(b-r)/4
			UPDATE(OCH_C31X, br+(rg>>2));	//b-(3*r+g)/4 = b-r-(g-r)/4
			UPDATE(OCH_CX13, -br-(gb>>2));	//r-(g+3*b)/4 = r-b-(g-b)/4
			UPDATE(OCH_C1X3, gb+(br>>2));	//g-(r+3*b)/4 = g-b-(r-b)/4
			UPDATE(OCH_C13X, -gb-(rg>>2));	//b-(r+3*g)/4 = b-g-(r-g)/4
			UPDATE(OCH_CX22, (rg-br)>>1);	//r-(g+b)/2 = (r-g + r-b)/2
			UPDATE(OCH_C2X2, (gb-rg)>>1);	//g-(r+b)/2 = (g-r + g-b)/2
			UPDATE(OCH_C22X, (br-gb)>>1);	//b-(r+g)/2 = (b-r + b-g)/2
#undef  UPDATE
		}
#ifdef ENTROPY0
		double csizes[OCH_COUNT]={0};
		double norm=ANALYSIS_STRIDE/((double)iw*ih);
		for(int kc=0;kc<OCH_COUNT;++kc)
		{
			int *curr_hist=hist[kc];
			double e=0;
			for(int ks=0;ks<256;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
					e-=freq*log2(freq*norm);
			}
			csizes[kc]=e/8;
		}
		double bestsize=0;
		for(int kt=0;kt<RCT_COUNT;++kt)
		{
			const unsigned char *rct=rct_combinations[kt];
			double size=
				+csizes[rct[0]]
				+csizes[rct[1]]
				+csizes[rct[2]]
			;
#ifdef LOUD
			printf("%-14s %12.2lf + %12.2lf + %12.2lf = %12.2lf\n", rct_names[kt], csizes[rct[0]], csizes[rct[1]], csizes[rct[2]], size);
#endif
			if(!kt||bestsize>size)
			{
				bestsize=size;
				bestrct=kt;
			}
		}
#else
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
			printf("%-14s %12lld + %12lld + %12lld = %12lld\n", rct_names[kt], counters[rct[0]], counters[rct[1]], counters[rct[2]], currerr);
#endif
			if(!kt||minerr>currerr)
			{
				minerr=currerr;
				bestrct=kt;
			}
		}
#endif
		//no renorm here
		nbits-=7;
		cache|=(unsigned long long)bestrct<<nbits;
#ifdef LOUD
		double t2=time_sec()-t;
		for(int k=0;k<OCH_COUNT;++k)
#ifdef ENTROPY0
			printf("%s %12.2lf\n", och_names[k], csizes[k]);
#else
			printf("%s %12lld\n", och_names[k], counters[k]);
#endif
		printf("A  %12.6lf sec  %12.6lf MB/s  %s\n", t2, srcsize/(t2*1024*1024), rct_names[bestrct]);
#else
		(void)och_names;
		(void)rct_names;
#endif
	}
	else
	{
		dstbufsize=usize;
		dstbuf=malloc(dstbufsize);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		streamptr=srcptr;
#ifdef _MSC_VER
		streamend=srcend+sizeof(__m256i[2]);
#endif
		image=dstbuf;
		
		cache=*(unsigned long long*)streamptr;
		streamptr+=8;
		nbits+=64-7;
		bestrct=cache>>nbits&((1ULL<<7)-1);
		cache&=(1ULL<<nbits)-1;
		if(bestrct>=RCT_COUNT)
		{
			LOG_ERROR("Invalid file");
			return 1;
		}
	}
	const unsigned char *combination=rct_combinations[bestrct];
	int
		yidx=combination[II_PERM_Y],
		uidx=combination[II_PERM_U],
		vidx=combination[II_PERM_V];
	int vfromy=-(combination[II_COEFF_U_SUB_Y]!=0);
	int yfromu=-(combination[II_COEFF_Y_SUB_U]!=0);
	int yfromv=-(combination[II_COEFF_Y_SUB_V]!=0);
	int vfromu=-(combination[II_COEFF_V_SUB_U_NBLI]!=0);
	int ufromv=-(combination[II_COEFF_U_SUB_V_NBLI]!=0);
	int psize=(int)sizeof(short[4*4*2])*(iw+16);//4 padded rows * 3 channels max * {pixels, errors}
	short *pixels=(short*)malloc(psize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(pixels, 0, psize);
#ifdef BYPASS_SSE
	memset(bypass_sse, 0, sizeof(bypass_sse));
#endif
	int paddedwidth=iw+16;
	unsigned char *ptr=image;
	for(int ky=0;ky<ih;++ky)
	{
		ALIGN(32) short *rows[]=
		{
			pixels+(paddedwidth*((ky-0LL)&3)+8LL)*4*2,
			pixels+(paddedwidth*((ky-1LL)&3)+8LL)*4*2,
			pixels+(paddedwidth*((ky-2LL)&3)+8LL)*4*2,
			pixels+(paddedwidth*((ky-3LL)&3)+8LL)*4*2,
		};
		ALIGN(16) short yuv[8]={0};
#ifdef BYPASS_SSE
		ALIGN(16) short preds[8]={0};
#endif
		ALIGN(16) int nbypass[4]={0};
		ALIGN(16) int nzeros[4]={0}, bypass[4]={0};
		__m128i msym=_mm_setzero_si128();
		__m128i mW=_mm_setzero_si128();
		__m128i meW=_mm_setzero_si128();
		for(int kx=0;kx<iw;++kx, ptr+=3)
		{
			__m128i mNW	=_mm_loadu_si128((__m128i*)(rows[1]-1*4*2));
			__m128i mN	=_mm_loadu_si128((__m128i*)(rows[1]+0*4*2));
		//	__m128i mW	=_mm_loadu_si128((__m128i*)(rows[0]-1*4*2));//15.42% 16.38%
			__m128i mmin=_mm_min_epi16(mN, mW);
			__m128i mmax=_mm_max_epi16(mN, mW);
			__m128i mp=_mm_add_epi16(mN, mW);
			mp=_mm_sub_epi16(mp, mNW);
			mp=_mm_max_epi16(mp, mmin);
			mp=_mm_min_epi16(mp, mmax);
#ifdef BYPASS_SSE
			_mm_store_si128((__m128i*)preds, mp);
#endif

			__m128i mnb=_mm_cvtepi16_epi32(meW);
#ifdef BYPASS_SSE
			int *ssectx=bypass_sse[(preds[0]>>2&0x30)|(preds[1]>>4&0x0C)|(preds[2]>>6&0x03)];
		//	int *ssectx=bypass_sse[(2*preds[0]+preds[1]+preds[2])>>4&63];
			__m128i msse=_mm_load_si128((__m128i*)ssectx);
			mnb=_mm_add_epi32(mnb, _mm_srai_epi32(msse, GRBITS+4));
#endif
			mnb=_mm_add_epi32(mnb, _mm_set1_epi32(1));
			mnb=_mm_sub_epi32(_mm_srli_epi32(_mm_castps_si128(_mm_cvtepi32_ps(mnb)), 23), _mm_set1_epi32(127+GRBITS));
			mnb=_mm_max_epi32(mnb, _mm_setzero_si128());
#ifdef BYPASS_SSE
			mnb=_mm_min_epi32(mnb, _mm_set_epi32(0, 9-1, 9-1, 8-1));
#endif
			//mnb=_mm_min_epi32(mnb, _mm_set_epi32(0, 9-1, 9-1, 8-1));//0<=nbypass<=depth-1		redundant
			//mnb=_mm_min_epi32(mnb, _mm_set_epi32(0, 9-3, 9-3, 8-3));//X  bad
			_mm_store_si128((__m128i*)nbypass, mnb);
			
#ifdef CALICSIGNPRED
			__m128i mhp=_mm_abs_epi16(mp);
#ifdef DISABLE_RCT
			mhp=_mm_sub_epi16(_mm_set1_epi16(128), mhp);
#else
			mhp=_mm_sub_epi16(_mm_set_epi16(0, 0, 0, 0, 0, 256, 256, 128), mhp);
#endif
#endif
#if 0
			short
			//	yeNEE	=rows[1][0+2*4*2+4],
			//	ueNEE	=rows[1][1+2*4*2+4],
			//	veNEE	=rows[1][2+2*4*2+4],
			//	yeNEEE	=rows[1][0+3*4*2+4],
			//	ueNEEE	=rows[1][1+3*4*2+4],
			//	veNEEE	=rows[1][2+3*4*2+4],
				yeW	=rows[0][0-1*4*2+4],
				ueW	=rows[0][1-1*4*2+4],
				veW	=rows[0][2-1*4*2+4];
			int ynbypass=FLOOR_LOG2(yeW+1)-GRBITS;
			int unbypass=FLOOR_LOG2(ueW+1)-GRBITS;
			int vnbypass=FLOOR_LOG2(veW+1)-GRBITS;
			if(ynbypass<0)ynbypass=0;
			if(unbypass<0)unbypass=0;
			if(vnbypass<0)vnbypass=0;
#ifdef CALICSIGNPRED
			int yhpred=128-abs(preds[0]);
			int uhpred=256-abs(preds[1]);
			int vhpred=256-abs(preds[2]);
#endif
#endif
			int ysym=0;
			int usym=0;
			int vsym=0;
			if(fwd)
			{
				//__m128i myuv=_mm_set_epi32(0, 0, 0, *(unsigned*)ptr);
				//myuv=_mm_cvtepi8_epi32(myuv);
				//myuv=_mm_sub_epi32(myuv, _mm_set1_epi32(128));

				yuv[0]=ptr[yidx]-128;
				yuv[1]=ptr[uidx]-128;
				yuv[2]=ptr[vidx]-128;
#ifndef DISABLE_RCT
				yuv[2]-=(combination[II_COEFF_V_SUB_Y]*yuv[0]+combination[II_COEFF_V_SUB_U]*yuv[1])>>2;//16.31%
				yuv[1]-=yuv[0]&vfromy;
				yuv[0]+=((yuv[1]&yfromu)+(yuv[2]&yfromv))>>2;
				yuv[2]-=(yuv[1]&vfromu)>>2;
				yuv[1]-=(yuv[2]&ufromv)>>2;
#endif
				__m128i myuv=_mm_load_si128((__m128i*)yuv);
				_mm_storeu_si128((__m128i*)rows[0], myuv);
				mW=myuv;
				//rows[0][0+0]=yuv[0];
				//rows[0][1+0]=yuv[1];
				//rows[0][2+0]=yuv[2];
				//mW=_mm_loadu_si128((__m128i*)yuv);
				
				__m128i me=_mm_sub_epi16(myuv, mp);
#ifdef CALICSIGNPRED
				__m128i ae=_mm_abs_epi16(me);
				__m128i pe=_mm_xor_si128(_mm_slli_epi16(me, 1), _mm_srai_epi16(me, 15));
				__m128i msum=_mm_add_epi16(ae, mhp);
				ae=_mm_cmpgt_epi16(ae, mhp);
				pe=_mm_blendv_epi8(pe, msum, ae);
#else
				__m128i pe=_mm_xor_si128(_mm_slli_epi16(me, 1), _mm_srai_epi16(me, 15));
#endif
				msym=pe;
				pe=_mm_cvtepi16_epi32(pe);
				__m128i mnz=_mm_srlv_epi32(pe, mnb);//AVX2
				_mm_store_si128((__m128i*)nzeros, mnz);
				__m128i bmask=_mm_sllv_epi32(_mm_set1_epi32(1), mnb);//AVX2
				bmask=_mm_sub_epi32(bmask, _mm_set1_epi32(1));
				pe=_mm_and_si128(pe, bmask);
				_mm_store_si128((__m128i*)bypass, pe);
#if 0
				int yerror=yuv[0]-preds[0];
				int uerror=yuv[1]-preds[1];
				int verror=yuv[2]-preds[2];
#ifdef CALICSIGNPRED
				int ynegmask=yerror>>31;
				int unegmask=uerror>>31;
				int vnegmask=verror>>31;
				int yabserr=(yerror^ynegmask)-ynegmask;
				int uabserr=(uerror^unegmask)-unegmask;
				int vabserr=(verror^vnegmask)-vnegmask;
				ysym=yerror<<1^ynegmask;
				usym=uerror<<1^unegmask;
				vsym=verror<<1^vnegmask;
				if(yhpred<yabserr)ysym=yhpred+yabserr;
				if(uhpred<uabserr)usym=uhpred+uabserr;
				if(vhpred<vabserr)vsym=vhpred+vabserr;
#else
				ysym=yerror<<1^yerror>>31;
				usym=uerror<<1^uerror>>31;
				vsym=verror<<1^verror>>31;
#endif
				int ynzeros=ysym>>ynbypass;
				int unzeros=usym>>unbypass;
				int vnzeros=vsym>>vnbypass;
				int ybypass=ysym&((1<<ynbypass)-1);
				int ubypass=usym&((1<<unbypass)-1);
				int vbypass=vsym&((1<<vnbypass)-1);
#endif
#ifdef LOUD
				bitcount[0][0]+=nzeros[0];
				bitcount[1][0]+=nzeros[1];
				bitcount[2][0]+=nzeros[2];
				++bitcount[0][1];
				++bitcount[1][1];
				++bitcount[2][1];
				bitcount[0][2]+=nbypass[0];
				bitcount[1][2]+=nbypass[1];
				bitcount[2][2]+=nbypass[2];
#endif
				//buffer: {c,c,c,b,b,a,a,a, f,f,f,e,e,e,d,c}, cache: MSB gg[hhh]000 LSB	nbits 6->3, code h is about to be emitted
				//written 64-bit words are byte-reversed because the CPU is little-endian
				if(nzeros[0]>=nbits)//fill the rest of cache with zeros, and flush
				{
					*(unsigned long long*)streamptr=cache;
					streamptr+=8;
					nzeros[0]-=nbits;
					cache=0;
					nbits=64;

					switch(nzeros[0]>>6&7)//just flush zeros
					{
					case 7:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 6:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 5:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 4:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 3:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 2:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 1:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 0:
						break;
					}
					nzeros[0]&=63;
					//while(nzeros[0]>=64)//just flush zeros
					//{
					//	nzeros[0]-=64;
					//	*(unsigned long long*)streamptr=0;
					//	streamptr+=8;
					//}
#if defined _MSC_VER && defined ENABLE_MEMCHECKS
					if(streamptr>streamend)
						LOG_ERROR("Out of memory XY %d %d", kx, ky);
#endif
				}
				//here  0 <= nzeros < nbits <= 64
				nbits-=nzeros[0]+1;//emit remaining zeros to cache
				//if(nbits<0)
				//	LOG_ERROR("");
				cache|=1ULL<<nbits;//append stop bit


				if(nzeros[1]>=nbits)//fill the rest of cache with zeros, and flush
				{
					*(unsigned long long*)streamptr=cache;
					streamptr+=8;
					nzeros[1]-=nbits;
					cache=0;
					nbits=64;
					switch(nzeros[1]>>6&7)//just flush zeros
					{
					case 7:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 6:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 5:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 4:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 3:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 2:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 1:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 0:
						break;
					}
					nzeros[1]&=63;
					//while(nzeros[1]>=64)//just flush zeros
					//{
					//	nzeros[1]-=64;
					//	*(unsigned long long*)streamptr=0;
					//	streamptr+=8;
					//}
				}
				//here  0 <= nzeros < nbits <= 64
				nbits-=nzeros[1]+1;//emit remaining zeros to cache
				//if(nbits<0)
				//	LOG_ERROR("");
				cache|=1ULL<<nbits;//append stop bit


				if(nzeros[2]>=nbits)//fill the rest of cache with zeros, and flush
				{
					*(unsigned long long*)streamptr=cache;
					streamptr+=8;
					nzeros[2]-=nbits;
					cache=0;
					nbits=64;

					switch(nzeros[2]>>6&7)//just flush zeros
					{
					case 7:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 6:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 5:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 4:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 3:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 2:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 1:*(unsigned long long*)streamptr=0; streamptr+=8;
					case 0:
						break;
					}
					nzeros[2]&=63;
					//while(nzeros[2]>=64)//just flush zeros
					//{
					//	nzeros[2]-=64;
					//	*(unsigned long long*)streamptr=0;
					//	streamptr+=8;
					//}
#if defined _MSC_VER && defined ENABLE_MEMCHECKS
					if(streamptr>streamend)
						LOG_ERROR("Out of memory XY %d %d", kx, ky);
#endif
				}
				//here  0 <= nzeros < nbits <= 64
				nbits-=nzeros[2]+1;//emit remaining zeros to cache
				//if(nbits<0)
				//	LOG_ERROR("");
				cache|=1ULL<<nbits;//append stop bit

				int totalnbypass=nbypass[0]+nbypass[1]+nbypass[2];
				int totalbypass=bypass[0];
				totalbypass=totalbypass<<nbypass[1]|bypass[1];
				totalbypass=totalbypass<<nbypass[2]|bypass[2];
				if(totalnbypass>=nbits)//not enough free bits in cache:  fill cache, write to list, and repeat
				{
					totalnbypass-=nbits;
					cache|=(unsigned long long)totalbypass>>totalnbypass;
					totalbypass&=(1<<totalnbypass)-1;
					*(unsigned long long*)streamptr=cache;
					streamptr+=8;
					cache=0;
					nbits=64;
#if defined _MSC_VER && defined ENABLE_MEMCHECKS
					if(streamptr>streamend)
						LOG_ERROR("Out of memory XY %d %d", kx, ky);
#endif
				}
				//now there is room for bypass:  0 <= nbypass < nbits <= 64
				nbits-=totalnbypass;//emit remaining bypass to cache
				cache|=(unsigned long long)totalbypass<<nbits;
			}
			else
			{
				//cache: MSB 00[hhh]ijj LSB		nbits 6->3, h is about to be read (past codes must be cleared from cache)

				ysym=nbits;
				nbits=(int)_lzcnt_u64(cache);
				//nbits=(int)_lzcnt_u64(cache);//dec: nbits = number of read bits	X  don't shift data. how to tell apart {0b0000|00000, 0100...} from {0b000001|000, 00...}?
				ysym-=64-nbits;
				nbits=64-nbits;//number of unread bits
				if(!cache)//must encounter stop bit
				{
					unsigned char *p2=streamptr+8;
					
#if defined _MSC_VER && defined ENABLE_MEMCHECKS
					if(streamptr+sizeof(__m256i[2])>streamend)
						LOG_ERROR("Out of memory XY %d %d", kx, ky);
#endif
					__m256i v0=_mm256_loadu_si256((__m256i*)streamptr+0);
					__m256i v1=_mm256_loadu_si256((__m256i*)streamptr+1);
					v0=_mm256_cmpeq_epi64(v0, _mm256_setzero_si256());
					v1=_mm256_cmpeq_epi64(v1, _mm256_setzero_si256());
					int m0=_mm256_movemask_pd(_mm256_castsi256_pd(v0));
					int m1=_mm256_movemask_pd(_mm256_castsi256_pd(v1));
					streamptr+=(_tzcnt_u16(~(m1<<4|m0))<<3)+8;//advance by [8~72] bytes
					cache=*(unsigned long long*)(streamptr-8);
					//do
					//{
					//	cache=*(unsigned long long*)streamptr;
					//	streamptr+=8;
					//}
					//while(!cache);
					ysym+=((int)(streamptr-p2)<<(6-3));//increments of 64 instead of 8
					nbits=(int)_lzcnt_u64(cache);
					ysym+=nbits;
					nbits=64-nbits;
				}
				--nbits;
				cache-=1ULL<<nbits;//remove stop bit

				usym=nbits;
				nbits=(int)_lzcnt_u64(cache);
				usym-=64-nbits;
				nbits=64-nbits;//number of unread bits
				if(!cache)//must encounter stop bit
				{
					unsigned char *p2=streamptr+8;
					
#if defined _MSC_VER && defined ENABLE_MEMCHECKS
					if(streamptr+sizeof(__m256i[2])>streamend)
						LOG_ERROR("Out of memory XY %d %d", kx, ky);
#endif
					__m256i v0=_mm256_loadu_si256((__m256i*)streamptr+0);
					__m256i v1=_mm256_loadu_si256((__m256i*)streamptr+1);
					v0=_mm256_cmpeq_epi64(v0, _mm256_setzero_si256());
					v1=_mm256_cmpeq_epi64(v1, _mm256_setzero_si256());
					int m0=_mm256_movemask_pd(_mm256_castsi256_pd(v0));
					int m1=_mm256_movemask_pd(_mm256_castsi256_pd(v1));
					streamptr+=(_tzcnt_u16(~(m1<<4|m0))<<3)+8;//advance by [8~72] bytes
					cache=*(unsigned long long*)(streamptr-8);
					//do
					//{
					//	cache=*(unsigned long long*)streamptr;
					//	streamptr+=8;
					//}
					//while(!cache);
					usym+=((int)(streamptr-p2)<<(6-3));//increments of 64 instead of 8
					nbits=(int)_lzcnt_u64(cache);
					usym+=nbits;
					nbits=64-nbits;
				}
				--nbits;
				cache-=1ULL<<nbits;//remove stop bit

				vsym=nbits;
				nbits=(int)_lzcnt_u64(cache);
				vsym-=64-nbits;
				nbits=64-nbits;//number of unread bits
				if(!cache)//must encounter stop bit
				{
					unsigned char *p2=streamptr+8;
					
#if defined _MSC_VER && defined ENABLE_MEMCHECKS
					if(streamptr+sizeof(__m256i[2])>streamend)
						LOG_ERROR("Out of memory XY %d %d", kx, ky);
#endif
					__m256i v0=_mm256_loadu_si256((__m256i*)streamptr+0);
					__m256i v1=_mm256_loadu_si256((__m256i*)streamptr+1);
					v0=_mm256_cmpeq_epi64(v0, _mm256_setzero_si256());
					v1=_mm256_cmpeq_epi64(v1, _mm256_setzero_si256());
					int m0=_mm256_movemask_pd(_mm256_castsi256_pd(v0));
					int m1=_mm256_movemask_pd(_mm256_castsi256_pd(v1));
					streamptr+=(_tzcnt_u16(~(m1<<4|m0))<<3)+8;//advance by [8~72] bytes
					cache=*(unsigned long long*)(streamptr-8);
					//do
					//{
					//	cache=*(unsigned long long*)streamptr;
					//	streamptr+=8;
					//}
					//while(!cache);
					vsym+=((int)(streamptr-p2)<<(6-3));//increments of 64 instead of 8
					nbits=(int)_lzcnt_u64(cache);
					vsym+=nbits;
					nbits=64-nbits;
				}
				--nbits;
				cache-=1ULL<<nbits;//remove stop bit

				int totalbypass=0;
				int totalnbypass=nbypass[0]+nbypass[1]+nbypass[2];
				if(nbits<totalnbypass)
				{
					totalnbypass-=nbits;
					totalbypass|=(int)(cache<<totalnbypass);
					cache=*(unsigned long long*)streamptr;
					streamptr+=8;
					nbits=64;
#if defined _MSC_VER && defined ENABLE_MEMCHECKS
					if(streamptr>streamend)
						LOG_ERROR("Out of memory XY %d %d", kx, ky);
#endif
				}
				if(totalnbypass)
				{
					nbits-=totalnbypass;
					totalbypass|=(int)(cache>>nbits);
					cache&=(1ULL<<nbits)-1;
				}
				vsym=vsym<<nbypass[2]|(totalbypass&((1<<nbypass[2])-1));
				totalbypass>>=nbypass[2];
				usym=usym<<nbypass[1]|(totalbypass&((1<<nbypass[1])-1));
				totalbypass>>=nbypass[1];
				ysym=ysym<<nbypass[0]|(totalbypass&((1<<nbypass[0])-1));

				__m128i ps=_mm_set_epi16(0, 0, 0, 0, 0, vsym, usym, ysym);
				msym=ps;
#ifdef CALICSIGNPRED
				__m128i negmask=_mm_srai_epi16(mp, 15);
				__m128i pe2=_mm_sub_epi16(mhp, ps);
				__m128i pe=_mm_xor_si128(_mm_srli_epi16(ps, 1), _mm_sub_epi16(_mm_setzero_si128(), _mm_and_si128(ps, _mm_set1_epi16(1))));
				pe2=_mm_xor_si128(pe2, negmask);
				pe2=_mm_sub_epi16(pe2, negmask);
				ps=_mm_cmpgt_epi16(ps, _mm_slli_epi16(mhp, 1));
				pe=_mm_blendv_epi8(pe, pe2, ps);
				pe=_mm_add_epi16(pe, mp);
#else
				__m128i pe=_mm_xor_si128(_mm_srli_epi16(ps, 1), _mm_sub_epi16(_mm_setzero_si128(), _mm_and_si128(ps, _mm_set1_epi16(1))));
#endif
				mW=pe;
				_mm_storeu_si128((__m128i*)rows[0], pe);
				_mm_store_si128((__m128i*)yuv, pe);
#if 0
#ifdef CALICSIGNPRED
				int ynegmask=preds[0]>>31;
				int unegmask=preds[1]>>31;
				int vnegmask=preds[2]>>31;
				int ye2=yhpred-ysym;
				int ue2=uhpred-usym;
				int ve2=vhpred-vsym;
				int yerror=ysym>>1^-(ysym&1);
				int uerror=usym>>1^-(usym&1);
				int verror=vsym>>1^-(vsym&1);
				ye2=(ye2^ynegmask)-ynegmask;
				ue2=(ue2^unegmask)-unegmask;
				ve2=(ve2^vnegmask)-vnegmask;
				if((yhpred<<1)<ysym)yerror=ye2;
				if((uhpred<<1)<usym)uerror=ue2;
				if((vhpred<<1)<vsym)verror=ve2;
#else
				int yerror=ysym>>1^-(ysym&1);
				int uerror=usym>>1^-(usym&1);
				int verror=vsym>>1^-(vsym&1);
#endif
				yuv[0]=yerror+preds[0];
				yuv[1]=uerror+preds[1];
				yuv[2]=verror+preds[2];
				rows[0][0+0]=yuv[0];
				rows[0][1+0]=yuv[1];
				rows[0][2+0]=yuv[2];
				mW=_mm_loadu_si128((__m128i*)yuv);
#endif
#ifndef DISABLE_RCT
				yuv[1]+=(yuv[2]&ufromv)>>2;
				yuv[2]+=(yuv[1]&vfromu)>>2;
				yuv[0]-=((yuv[1]&yfromu)+(yuv[2]&yfromv))>>2;
				yuv[1]+=yuv[0]&vfromy;
				yuv[2]+=(combination[II_COEFF_V_SUB_Y]*yuv[0]+combination[II_COEFF_V_SUB_U]*yuv[1])>>2;
#endif
				ptr[yidx]=yuv[0]+128;
				ptr[uidx]=yuv[1]+128;
				ptr[vidx]=yuv[2]+128;
				guide_check(image, kx, ky);
			}
			meW=_mm_slli_epi16(meW, 1);
			meW=_mm_add_epi16(meW, _mm_slli_epi16(msym, GRBITS));
			__m128i mbias=_mm_max_epi16(_mm_loadu_si128((__m128i*)(rows[1]+0+2*4*2+4)), _mm_loadu_si128((__m128i*)(rows[1]+0+3*4*2+4)));
			//__m128i mbias=_mm_loadu_si128((__m128i*)(rows[1]+0+3*4*2+4));//X
			meW=_mm_srli_epi16(_mm_add_epi16(meW, mbias), 2);
			_mm_storeu_si128((__m128i*)(rows[0]+0+4), meW);
			//rows[0][0+4]=(2*yeW+(ysym<<GRBITS)+MAXVAR(rows[1][0+2*4*2+4], rows[1][0+3*4*2+4]))>>2;
			//rows[0][1+4]=(2*ueW+(usym<<GRBITS)+MAXVAR(rows[1][1+2*4*2+4], rows[1][1+3*4*2+4]))>>2;
			//rows[0][2+4]=(2*veW+(vsym<<GRBITS)+MAXVAR(rows[1][2+2*4*2+4], rows[1][2+3*4*2+4]))>>2;
			//rows[0][0+4]=(2*yeW+(ysym<<GRBITS)+MAXVAR(yeNEE, yeNEEE))>>2;
			//rows[0][1+4]=(2*ueW+(usym<<GRBITS)+MAXVAR(ueNEE, ueNEEE))>>2;
			//rows[0][2+4]=(2*veW+(vsym<<GRBITS)+MAXVAR(veNEE, veNEEE))>>2;
#ifdef BYPASS_SSE
			msym=_mm_sub_epi16(_mm_slli_epi16(msym, GRBITS), meW);
			msym=_mm_cvtepi16_epi32(msym);
			msym=_mm_slli_epi32(msym, 4);
			msym=_mm_sub_epi32(msym, msse);
			msym=_mm_add_epi32(msym, _mm_set1_epi32(1<<5>>1));
			msym=_mm_srai_epi32(msym, 5);
			msse=_mm_add_epi32(msse, msym);
			_mm_store_si128((__m128i*)ssectx, msse);
#endif
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			LOG_ERROR("Cannot open \"%s\" for writing", fdst);
			return 1;
		}
		if(fwd)
		{
			*(unsigned long long*)streamptr=(unsigned long long)cache;
			streamptr+=8;

			fwrite("27", 1, 2, fdst);
			fwrite(&iw, 1, 4, fdst);
			fwrite(&ih, 1, 4, fdst);
			fwrite(dstbuf, 1, streamptr-dstbuf, fdst);
		}
		else
		{
			srcsize=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(image, 1, usize, fdst);
			srcsize+=usize;
		}
		fclose(fdst);
	}
	free(pixels);
	free(dstbuf);
	free(srcbuf);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		printf("%s\n", srcfn);
		printf("Flags: (usize %td)\n", srcsize);
		printf("Y %12.2lf %12.2lf %12.2lf %12.2lf\n", bitcount[0][0]/8., bitcount[0][1]/8., bitcount[0][2]/8., (bitcount[0][0]+bitcount[0][1]+bitcount[0][2])/8.);
		printf("U %12.2lf %12.2lf %12.2lf %12.2lf\n", bitcount[1][0]/8., bitcount[1][1]/8., bitcount[1][2]/8., (bitcount[1][0]+bitcount[1][1]+bitcount[1][2])/8.);
		printf("V %12.2lf %12.2lf %12.2lf %12.2lf\n", bitcount[2][0]/8., bitcount[2][1]/8., bitcount[2][2]/8., (bitcount[2][0]+bitcount[2][1]+bitcount[2][2])/8.);
		printf("%8td/%8td bytes\n", 2+4+4+streamptr-dstbuf, srcsize);
	}
	printf("%c  %12.6lf sec  %12.6lf MB/s\n", 'D'+fwd, t, srcsize/(t*1024*1024));
#endif

	(void)nthreads0;
#ifndef LOUD
	(void)rct_names;
#endif
	return 0;
}