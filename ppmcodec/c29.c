#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


#ifdef _MSC_VER
	#define LOUD
	#define ENABLE_MEMCHECKS
	#define ENABLE_GUIDE
#endif

#define NCODERS 32
#define GRBITS 3

//CDF2sym with 3*2*(8+3)+1=67 contexts:  16bit -> 4+ MB  12 bit -> 268 KB
#define PROB_BITS 12	//15 bit max

//DIV-free state is 31 bits
#define RANS_BYTE_L (1<<15)	//rANS decoder renorm lower bound


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


//WG:

#define WG_NPREDS	8	//multiple of 4
#if 0
#define WG_PREDLIST0\
	WG_PRED(340,	N)\
	WG_PRED(340,	W)\
	WG_PRED(205,	3*(N-NN)+NNN)\
	WG_PRED(205,	3*(W-WW)+WWW)\
	WG_PRED(140,	W+NE-N)\
	WG_PRED(240,	(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4)\
	WG_PRED(120,	N+W-NW)\
	WG_PRED(120,	N+NE-NNE)
#define WG_PREDLIST1\
	WG_PRED(330,	N)\
	WG_PRED(330,	W)\
	WG_PRED(175,	3*(N-NN)+NNN)\
	WG_PRED(175,	3*(W-WW)+WWW)\
	WG_PRED(180,	W+NE-N)\
	WG_PRED(175,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW)/3)\
	WG_PRED(130,	N+W-NW)\
	WG_PRED(150,	N+NE-NNE)
#define WG_PREDLIST2\
	WG_PRED(330,	N)\
	WG_PRED(330,	W)\
	WG_PRED(200,	3*(N-NN)+NNN)\
	WG_PRED(200,	3*(W-WW)+WWW)\
	WG_PRED(180,	W+NE-N)\
	WG_PRED(175,	(WWW+NNN+NEE+NEEE+NEEEE-2*NW)/3)\
	WG_PRED(140,	N+W-NW)\
	WG_PRED(150,	N+NE-NNE)
#endif

//https://github.com/rygorous/ryg_rans
//https://github.com/samtools/htscodecs
typedef struct _rANS_SIMD_SymInfo	//16 bytes/level	4KB/ctx = 1<<12 bytes
{
	unsigned vmax, invf, cdf;
	unsigned short negf, sh;
} rANS_SIMD_SymInfo;
static void enc_hist2stats(const int *hist, rANS_SIMD_SymInfo *syminfo)
{
	int sum=0;
	for(int ks=0;ks<256;++ks)
		sum+=hist[ks];
	int sum2=0;
	for(int ks=0;ks<256;++ks)
	{
		int freq=hist[ks];
		syminfo[ks].cdf=(int)(sum2*((1ULL<<PROB_BITS)-256)/sum)+ks;
		sum2+=freq;
	}
	int next=1<<PROB_BITS;
	for(int ks=255;ks>=0;--ks)
	{
		rANS_SIMD_SymInfo *info=syminfo+ks;
		int curr=info->cdf;
		int freq=next-curr;
		next=curr;
		info->vmax=(RANS_BYTE_L>>PROB_BITS<<16)*freq-1;
		info->negf=(1<<PROB_BITS)-freq;
		//encoding:  state  =  q<<16|(cdf+r)
		//div-free:  state  =  q*M+cdf+state-q*freq  =  state+q*(M-freq)+cdf  =  state+(state*invf>>sh)*(M-freq)+cdf
		//sh = FLOOR_LOG2(freq)+32
		//invf = ceil(2^sh/freq)		state is 31 bits
		if(freq<2)
		{
			//ideally  q = state*inv(1)>>sh(1) = state*2^32>>32
			//here  q' = state*(2^32-1)>>32 = floor(state-state/2^32) = state-1  if  1 <= x < 2^32
			//enc  state = (state/1)*M+cdf+state%1  =  state+q*(M-1)+cdf
			//but  q' = state-1
			//so  state = state+(state-1+1)*(M-1)+cdf  =  state+q'*(M-1)+(cdf+M-1)
			info->sh=0;
			info->invf=0xFFFFFFFF;
			info->cdf+=(1<<PROB_BITS)-1;
		}
		else
		{
			info->sh=63-_lzcnt_u32(freq);//FLOOR_LOG2(freq)+32
			info->invf=((1ULL<<info->sh)+freq-1)/freq;
		}
	}
}
int c29_codec(const char *srcfn, const char *dstfn, int nthreads0)
{
#ifdef LOUD
	double t=time_sec();
#endif
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	int fwd=0, iw=0, ih=0, rowstride=0, blockw=0, rwidth=0;
	ptrdiff_t usize=0, cap=0;
	unsigned char *image=0, *imptr=0, *streamptr=0, *streamstart=0;
	unsigned char *context=0;
	int psize=0;
	short *pixels=0;
	int *wgsize=0;
	short *wgerrors=0;
	ptrdiff_t uheadersize=0, cheadersize=0, csize=0;
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
			uheadersize=ftell(fsrc);
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
		blockw=iw/NCODERS;
		rwidth=iw%NCODERS;
		rowstride=3*iw;
		usize=(ptrdiff_t)3*iw*ih;
		cap=(ptrdiff_t)4*iw*ih;
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
		}
		else
		{
			csize=get_filesize(fsrc);
			streamptr=streamstart=image+cap-sizeof(__m256i)-(csize-cheadersize);
			fread(streamstart, 1, csize-cheadersize, fsrc);//read stream
		}
		fclose(fsrc);
	}
	int bestrct=0, use_wg4=0;
	int hsize=(int)sizeof(int[3*2*(8+GRBITS)<<8]);//3 channels  *  2*(8+GRBITS) contexts
	int *hists=(int*)malloc(hsize);
	int rhsize=(int)sizeof(int[3*256]);
	int *rhist=(int*)malloc(rhsize);

	int tempsize;
	int CDF2syms_size=(int)sizeof(char[3*2*(8+GRBITS)<<PROB_BITS]);//DIV-free ANS encoder reuses these as SIMD symbol info
	tempsize=(int)sizeof(rANS_SIMD_SymInfo[3*2*(8+GRBITS)<<8]);
	if(CDF2syms_size<tempsize)CDF2syms_size=tempsize;
	unsigned char *CDF2syms=(unsigned char*)_mm_malloc(CDF2syms_size, sizeof(__m128i*));

	int rCDF2syms_size=(int)sizeof(char[3<<PROB_BITS]);
	tempsize=(int)sizeof(rANS_SIMD_SymInfo[3<<8]);
	if(rCDF2syms_size<tempsize)rCDF2syms_size=tempsize;
	unsigned char *rCDF2syms=(unsigned char*)_mm_malloc(rCDF2syms_size, sizeof(__m128i*));

	psize=(int)sizeof(short[4*3*NCODERS*2])*(blockw+16);//4 padded rows  *  {Y*NCODERS, eY*NCODERS,  U*NCODERS, eU*NCODERS,  V*NCODERS, eV*NCODERS} = 3*2*32 = 192 channels  ~48*iw bytes
	pixels=(short*)_mm_malloc(psize, sizeof(__m256i));
	wgsize=(int)sizeof(short[2*3*WG_NPREDS*NCODERS])*(blockw+16);//2 padded rows  *  {WGY*NCODERS, WGU*NCODERS, WGV*NCODERS} = 3*8*32 = 768 channels  ~96*iw bytes	total ~570 KB for 4K/12MP
	wgerrors=(short*)_mm_malloc(wgsize, sizeof(__m256i));
	if(!hists||!rhist||!pixels||!wgerrors||!CDF2syms||!rCDF2syms)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(hists, 0, hsize);
	memset(rhist, 0, rhsize);
	if(fwd)
	{
		ALIGN(32) long long counters[OCH_COUNT]={0};
		{
			unsigned char *temprow=(unsigned char*)pixels;
			unsigned char *tptrs[NCODERS]={0}, *tptrs0[NCODERS]={0};
			{
				unsigned char *temp=temprow;
				for(int k=0;k<NCODERS;++k)
				{
					tptrs[k]=temp;
					temp+=blockw;
				}
			}
			memcpy(tptrs0, tptrs, sizeof(tptrs0));
			imptr=image;
			for(int ky=0;ky<ih;++ky)//interleave coder lanes
			{
				unsigned char *dstptr=imptr;
				memcpy(temprow, imptr, rowstride);
				memcpy(tptrs, tptrs0, sizeof(tptrs));
				int kx=0;
				for(;kx<=rowstride-3*NCODERS;kx+=3*NCODERS)
				{
					//toy example with 2 coders:
					// p0  p2  p4                                      p1  p3  p5
					//|r00 g00 b00|r01 g01 b01|r02 g02 b02|r03 g03 b03|r04 g04 b04|r05 g05 b05|r06 g06 b06|r07 g07 b07|
					//|r00 r04|g00 g04|b00 b04|r01 r05|g01 g05|b01 b05|r02 r06|g02 g06|b02 b06|r03 r07|g03 g07|b03 b07|
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
					for(int k=0;k<NCODERS;++k)
						*dstptr++=*tptrs[k]++;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
					for(int k=0;k<NCODERS;++k)
						*dstptr++=*tptrs[k]++;
#if defined __GNUC__ && !defined PROFILER
#pragma GCC unroll 32
#endif
					for(int k=0;k<NCODERS;++k)
						*dstptr++=*tptrs[k]++;
				}
				//if(kx<rowstride)
				//{
				//	do
				//	{
				//		kx+=3;
				//	}
				//	while(kx<rowstride);
				//}
				imptr+=rowstride;
			}
		}
		guide_save(image, iw, ih);
		{
			__m256i mcounters[OCH_COUNT];//64-bit
			__m128i half8=_mm_set1_epi8(-128);
			__m256i wordmask=_mm256_set1_epi64x(0xFFFF);
			imptr=image;
			for(int ky=0;ky<ih;++ky)//analysis
			{
				__m256i prev[OCH_COUNT*2];//16-bit
				memset(prev, 0, sizeof(prev));
				for(int kx=0;kx<rowstride-3*NCODERS;kx+=3*NCODERS)
				{
					__m256i r0=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx)+0), half8));
					__m256i r1=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx)+1), half8));
					__m256i g0=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx)+2), half8));
					__m256i g1=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx)+3), half8));
					__m256i b0=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx)+4), half8));
					__m256i b1=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx)+5), half8));
					r0=_mm256_slli_epi16(r0, 2);
					r1=_mm256_slli_epi16(r1, 2);
					g0=_mm256_slli_epi16(g0, 2);
					g1=_mm256_slli_epi16(g1, 2);
					b0=_mm256_slli_epi16(b0, 2);
					b1=_mm256_slli_epi16(b1, 2);
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
		__m256i a0=_mm256_sub_epi16(A0, prev[IDXA*2+0]);\
		__m256i a1=_mm256_sub_epi16(A1, prev[IDXA*2+1]);\
		__m256i b0=_mm256_sub_epi16(B0, prev[IDXB*2+0]);\
		__m256i b1=_mm256_sub_epi16(B1, prev[IDXB*2+1]);\
		__m256i c0=_mm256_sub_epi16(C0, prev[IDXC*2+0]);\
		__m256i c1=_mm256_sub_epi16(C1, prev[IDXC*2+1]);\
		prev[IDXA*2+0]=A0;\
		prev[IDXA*2+1]=A1;\
		prev[IDXB*2+0]=B0;\
		prev[IDXB*2+1]=B1;\
		prev[IDXC*2+0]=C0;\
		prev[IDXC*2+1]=C1;\
		a0=_mm256_abs_epi16(a0);\
		a1=_mm256_abs_epi16(a1);\
		b0=_mm256_abs_epi16(b0);\
		b1=_mm256_abs_epi16(b1);\
		c0=_mm256_abs_epi16(c0);\
		c1=_mm256_abs_epi16(c1);\
		a0=_mm256_add_epi16(a0, a1);\
		b0=_mm256_add_epi16(b0, b1);\
		c0=_mm256_add_epi16(c0, c1);\
		a0=_mm256_add_epi16(a0, _mm256_srli_epi64(a0, 32));\
		b0=_mm256_add_epi16(b0, _mm256_srli_epi64(b0, 32));\
		c0=_mm256_add_epi16(c0, _mm256_srli_epi64(c0, 32));\
		a0=_mm256_add_epi16(a0, _mm256_srli_epi64(a0, 16));\
		b0=_mm256_add_epi16(b0, _mm256_srli_epi64(b0, 16));\
		c0=_mm256_add_epi16(c0, _mm256_srli_epi64(c0, 16));\
		mcounters[IDXA]=_mm256_add_epi64(mcounters[IDXA], _mm256_and_si256(a0, wordmask));\
		mcounters[IDXB]=_mm256_add_epi64(mcounters[IDXB], _mm256_and_si256(b0, wordmask));\
		mcounters[IDXC]=_mm256_add_epi64(mcounters[IDXC], _mm256_and_si256(c0, wordmask));\
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
				imptr+=rowstride;
			}
			for(int k=0;k<OCH_COUNT;++k)
			{
				ALIGN(32) long long temp[8]={0};
				_mm256_store_si256((__m256i*)temp+0, mcounters[k*2+0]);
				_mm256_store_si256((__m256i*)temp+1, mcounters[k*2+1]);
				counters[k]=temp[0]+temp[1]+temp[2]+temp[3]+temp[4]+temp[5]+temp[6]+temp[7];
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
				printf("%-14s %12lld + %12lld + %12lld = %12lld\n", rct_names[kt], counters[rct[0]], counters[rct[1]], counters[rct[2]], currerr);
#endif
				if(!kt||minerr>currerr)
				{
					minerr=currerr;
					bestrct=kt;
				}
			}
		}
		{
			const unsigned char *rct=rct_combinations[bestrct];
			long long count=3*((long long)iw-rwidth)*ih * 2;//FIXME tune
			use_wg4=counters[rct[0]]+counters[rct[1]]+counters[rct[2]]>count;
		}
		context=(unsigned char*)malloc(usize+sizeof(__m256i));
		if(!context)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
	}
	else
	{
		//FIXME decode flags, stats
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
	__m256i mctxmax=_mm256_set1_epi16(2*(GRBITS+8)-1);
	__m256i mctxuoffset=_mm256_set1_epi16(2*(GRBITS+8));
	__m256i mctxvoffset=_mm256_set1_epi16(2*(GRBITS+8)*2);
	__m256i amin=_mm256_set1_epi16(-128);
	__m256i amax=_mm256_set1_epi16(127);
	__m256i half16=_mm256_set1_epi16(128);
	__m128i half8=_mm_set1_epi8(-128);
	__m256i bytemask=_mm256_set1_epi16(255);
	__m256i wordmask=_mm256_set1_epi64x(0xFFFF);
	__m256i myuv[6];
	memset(myuv, 0, sizeof(myuv));
	imptr=image;
	unsigned char *ctxptr=context;
	for(int ky=0;ky<ih;++ky)//main coding loop
	{
		ALIGN(32) short *rows[]=
		{
			pixels+(paddedwidth*((ky-0LL)&3)+8LL)*3*NCODERS*2,
			pixels+(paddedwidth*((ky-1LL)&3)+8LL)*3*NCODERS*2,
			pixels+(paddedwidth*((ky-2LL)&3)+8LL)*3*NCODERS*2,
			pixels+(paddedwidth*((ky-3LL)&3)+8LL)*3*NCODERS*2,
		};
		ALIGN(32) short *erows[]=
		{
			wgerrors+(paddedwidth*((ky-0LL)&1)+8LL)*3*WG_NPREDS*NCODERS,
			wgerrors+(paddedwidth*((ky-1LL)&1)+8LL)*3*WG_NPREDS*NCODERS,
		};
		ALIGN(32) unsigned short syms[32]={0};
		__m256i NW[6], N[6], W[6];
		__m256i eW[6], ecurr[6], eNEE[6], eNEEE[6];
		memset(NW, 0, sizeof(NW));
		memset(W, 0, sizeof(W));
		memset(eW, 0, sizeof(eW));
		memset(ecurr, 0, sizeof(ecurr));
		eNEE[0]=_mm256_load_si256((__m256i*)(rows[1]+2*3*NCODERS*2)+0+1*2);
		eNEE[1]=_mm256_load_si256((__m256i*)(rows[1]+2*3*NCODERS*2)+1+1*2);
		eNEE[2]=_mm256_load_si256((__m256i*)(rows[1]+2*3*NCODERS*2)+0+3*2);
		eNEE[3]=_mm256_load_si256((__m256i*)(rows[1]+2*3*NCODERS*2)+1+3*2);
		eNEE[4]=_mm256_load_si256((__m256i*)(rows[1]+2*3*NCODERS*2)+0+5*2);
		eNEE[5]=_mm256_load_si256((__m256i*)(rows[1]+2*3*NCODERS*2)+1+5*2);
		memset(eNEEE, 0, sizeof(eNEEE));
		for(int kx=0;kx<rowstride-3*NCODERS;kx+=3*NCODERS)
		{
			N	[0]=_mm256_load_si256((__m256i*)(rows[1]+0*3*NCODERS*2)+0+0*2);//y0
			N	[1]=_mm256_load_si256((__m256i*)(rows[1]+0*3*NCODERS*2)+1+0*2);//y1
			N	[2]=_mm256_load_si256((__m256i*)(rows[1]+0*3*NCODERS*2)+0+2*2);//u0
			N	[3]=_mm256_load_si256((__m256i*)(rows[1]+0*3*NCODERS*2)+1+2*2);//u1
			N	[4]=_mm256_load_si256((__m256i*)(rows[1]+0*3*NCODERS*2)+0+4*2);//v0
			N	[5]=_mm256_load_si256((__m256i*)(rows[1]+0*3*NCODERS*2)+1+4*2);//v1
			__m256i
				predY0, predY1, ctxY0, ctxY1,
				predU0, predU1, ctxU0, ctxU1,
				predV0, predV1, ctxV0, ctxV1;
			{
				__m256i one=_mm256_set1_epi16(1);
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
				cy0=_mm256_add_epi16(cy0, one);
				cy1=_mm256_add_epi16(cy1, one);
				cy2=_mm256_add_epi16(cy2, one);
				cy3=_mm256_add_epi16(cy3, one);
				cu0=_mm256_add_epi16(cu0, one);
				cu1=_mm256_add_epi16(cu1, one);
				cu2=_mm256_add_epi16(cu2, one);
				cu3=_mm256_add_epi16(cu3, one);
				cv0=_mm256_add_epi16(cv0, one);
				cv1=_mm256_add_epi16(cv1, one);
				cv2=_mm256_add_epi16(cv2, one);
				cv3=_mm256_add_epi16(cv3, one);
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
				//if(use_wg4)//FIXME
				if(0)
				{
					/*
					N
					W
					3*(N-NN)+NNN
					3*(W-WW)+WWW
					W+NE-N
					(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW)/4
					N+W-NW
					N+NE-NNE
					*/
					predY0=_mm256_setzero_si256();
					predY1=_mm256_setzero_si256();
					predU0=_mm256_setzero_si256();
					predU1=_mm256_setzero_si256();
					predV0=_mm256_setzero_si256();
					predV1=_mm256_setzero_si256();
					//ymin0=_mm256_min_epi16(ymin0, NE[0]);
					//ymax0=_mm256_max_epi16(ymax0, NE[0]);
					//ymin1=_mm256_min_epi16(ymin1, NE[1]);
					//ymax1=_mm256_max_epi16(ymax1, NE[1]);
					//umin0=_mm256_min_epi16(umin0, NE[2]);
					//umax0=_mm256_max_epi16(umax0, NE[2]);
					//umin1=_mm256_min_epi16(umin1, NE[3]);
					//umax1=_mm256_max_epi16(umax1, NE[3]);
					//vmin0=_mm256_min_epi16(vmin0, NE[4]);
					//vmax0=_mm256_max_epi16(vmax0, NE[4]);
					//vmin1=_mm256_min_epi16(vmin1, NE[5]);
					//vmax1=_mm256_max_epi16(vmax1, NE[5]);
				}
				else
				{
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
			__m256i msyms0, msyms1, msyms8, mctx, moffset0, moffset1;
			if(fwd)
			{
				myuv[0]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx+yidx)+0), half8));//load yuv
				myuv[1]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx+yidx)+1), half8));
				myuv[2]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx+uidx)+0), half8));
				myuv[3]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx+uidx)+1), half8));
				myuv[4]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx+vidx)+0), half8));
				myuv[5]=_mm256_cvtepi8_epi16(_mm_add_epi8(_mm_loadu_si128((__m128i*)(imptr+kx+vidx)+1), half8));
				
				msyms0=_mm256_sub_epi16(myuv[0], predY0);//sub pred
				msyms1=_mm256_sub_epi16(myuv[1], predY1);
				ecurr[0]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));//ecurr = pack_sign(yuv-pred)		FIXME try abs
				ecurr[1]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms0=_mm256_add_epi16(msyms0, half16);
				msyms1=_mm256_add_epi16(msyms1, half16);
				msyms0=_mm256_and_si256(msyms0, bytemask);//sym = (yuv-pred+128)&255
				msyms1=_mm256_and_si256(msyms1, bytemask);
				msyms8=_mm256_packus_epi16(msyms0, msyms1);//FIXME shift-or
				mctx=_mm256_packus_epi16(ctxY0, ctxY1);
				ctxY0=_mm256_slli_epi16(ctxY0, 8);
				ctxY1=_mm256_slli_epi16(ctxY1, 8);
				ctxY0=_mm256_or_si256(ctxY0, msyms0);
				ctxY1=_mm256_or_si256(ctxY1, msyms1);
				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
				mctx=_mm256_permute4x64_epi64(mctx, _MM_SHUFFLE(3, 1, 2, 0));
				//msyms8=_mm256_add_epi8(msyms8, half16);
				_mm256_storeu_si256((__m256i*)(imptr+kx+yidx), msyms8);//store residuals
				_mm256_storeu_si256((__m256i*)(ctxptr+kx+yidx), mctx);//store contexts
				_mm256_store_si256((__m256i*)syms+0, ctxY0);
				_mm256_store_si256((__m256i*)syms+1, ctxY1);
				W[0]=myuv[0];
				W[1]=myuv[1];
				_mm256_store_si256((__m256i*)(rows[0]+0*3*NCODERS*2)+0+0*2, myuv[0]);//store pixel
				_mm256_store_si256((__m256i*)(rows[0]+0*3*NCODERS*2)+1+0*2, myuv[1]);
				++hists[syms[0x00]];
				++hists[syms[0x01]];
				++hists[syms[0x02]];
				++hists[syms[0x03]];
				++hists[syms[0x04]];
				++hists[syms[0x05]];
				++hists[syms[0x06]];
				++hists[syms[0x07]];
				++hists[syms[0x08]];
				++hists[syms[0x09]];
				++hists[syms[0x0A]];
				++hists[syms[0x0B]];
				++hists[syms[0x0C]];
				++hists[syms[0x0D]];
				++hists[syms[0x0E]];
				++hists[syms[0x0F]];
				++hists[syms[0x10]];
				++hists[syms[0x11]];
				++hists[syms[0x12]];
				++hists[syms[0x13]];
				++hists[syms[0x14]];
				++hists[syms[0x15]];
				++hists[syms[0x16]];
				++hists[syms[0x17]];
				++hists[syms[0x18]];
				++hists[syms[0x19]];
				++hists[syms[0x1A]];
				++hists[syms[0x1B]];
				++hists[syms[0x1C]];
				++hists[syms[0x1D]];
				++hists[syms[0x1E]];
				++hists[syms[0x1F]];

				moffset0=_mm256_and_si256(myuv[0], uhelpmask);
				moffset1=_mm256_and_si256(myuv[1], uhelpmask);
				predU0=_mm256_add_epi16(predU0, moffset0);
				predU1=_mm256_add_epi16(predU1, moffset1);
				predU0=_mm256_max_epi16(predU0, amin);
				predU1=_mm256_max_epi16(predU1, amin);
				predU0=_mm256_min_epi16(predU0, amax);
				predU1=_mm256_min_epi16(predU1, amax);
				
				msyms0=_mm256_sub_epi16(myuv[2], predU0);
				msyms1=_mm256_sub_epi16(myuv[3], predU1);
				ecurr[2]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));
				ecurr[3]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms0=_mm256_add_epi16(msyms0, half16);
				msyms1=_mm256_add_epi16(msyms1, half16);
				msyms0=_mm256_and_si256(msyms0, bytemask);
				msyms1=_mm256_and_si256(msyms1, bytemask);
				msyms8=_mm256_packus_epi16(msyms0, msyms1);
				mctx=_mm256_packus_epi16(ctxU0, ctxU1);
				ctxU0=_mm256_add_epi16(ctxU0, mctxuoffset);
				ctxU1=_mm256_add_epi16(ctxU1, mctxuoffset);
				ctxU0=_mm256_slli_epi16(ctxU0, 8);
				ctxU1=_mm256_slli_epi16(ctxU1, 8);
				ctxU0=_mm256_or_si256(ctxU0, msyms0);
				ctxU1=_mm256_or_si256(ctxU1, msyms1);
				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
				mctx=_mm256_permute4x64_epi64(mctx, _MM_SHUFFLE(3, 1, 2, 0));
				_mm256_storeu_si256((__m256i*)(imptr+kx+yidx), msyms8);
				_mm256_storeu_si256((__m256i*)(ctxptr+kx+yidx), mctx);
				msyms0=_mm256_sub_epi16(myuv[2], moffset0);
				msyms1=_mm256_sub_epi16(myuv[3], moffset1);
				_mm256_store_si256((__m256i*)syms+0, ctxU0);
				_mm256_store_si256((__m256i*)syms+1, ctxU1);
				W[2]=msyms0;
				W[3]=msyms1;
				_mm256_store_si256((__m256i*)(rows[0]+0*3*NCODERS*2)+0+2*2, msyms0);
				_mm256_store_si256((__m256i*)(rows[0]+0*3*NCODERS*2)+1+2*2, msyms1);
				++hists[syms[0x00]];
				++hists[syms[0x01]];
				++hists[syms[0x02]];
				++hists[syms[0x03]];
				++hists[syms[0x04]];
				++hists[syms[0x05]];
				++hists[syms[0x06]];
				++hists[syms[0x07]];
				++hists[syms[0x08]];
				++hists[syms[0x09]];
				++hists[syms[0x0A]];
				++hists[syms[0x0B]];
				++hists[syms[0x0C]];
				++hists[syms[0x0D]];
				++hists[syms[0x0E]];
				++hists[syms[0x0F]];
				++hists[syms[0x10]];
				++hists[syms[0x11]];
				++hists[syms[0x12]];
				++hists[syms[0x13]];
				++hists[syms[0x14]];
				++hists[syms[0x15]];
				++hists[syms[0x16]];
				++hists[syms[0x17]];
				++hists[syms[0x18]];
				++hists[syms[0x19]];
				++hists[syms[0x1A]];
				++hists[syms[0x1B]];
				++hists[syms[0x1C]];
				++hists[syms[0x1D]];
				++hists[syms[0x1E]];
				++hists[syms[0x1F]];

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
				
				msyms0=_mm256_sub_epi16(myuv[4], predV0);
				msyms1=_mm256_sub_epi16(myuv[5], predV1);
				ecurr[4]=_mm256_xor_si256(_mm256_slli_epi16(msyms0, 1), _mm256_srai_epi16(msyms0, 15));
				ecurr[5]=_mm256_xor_si256(_mm256_slli_epi16(msyms1, 1), _mm256_srai_epi16(msyms1, 15));
				msyms0=_mm256_add_epi16(msyms0, half16);
				msyms1=_mm256_add_epi16(msyms1, half16);
				msyms0=_mm256_and_si256(msyms0, bytemask);
				msyms1=_mm256_and_si256(msyms1, bytemask);
				msyms8=_mm256_packus_epi16(msyms0, msyms1);
				mctx=_mm256_packus_epi16(ctxV0, ctxV1);
				ctxV0=_mm256_add_epi16(ctxV0, mctxvoffset);
				ctxV1=_mm256_add_epi16(ctxV1, mctxvoffset);
				ctxV0=_mm256_slli_epi16(ctxV0, 8);
				ctxV1=_mm256_slli_epi16(ctxV1, 8);
				ctxV0=_mm256_or_si256(ctxV0, msyms0);
				ctxV1=_mm256_or_si256(ctxV1, msyms1);
				msyms8=_mm256_permute4x64_epi64(msyms8, _MM_SHUFFLE(3, 1, 2, 0));
				mctx=_mm256_permute4x64_epi64(mctx, _MM_SHUFFLE(3, 1, 2, 0));
				_mm256_storeu_si256((__m256i*)(imptr+kx+yidx), msyms8);
				_mm256_storeu_si256((__m256i*)(ctxptr+kx+yidx), mctx);
				_mm256_store_si256((__m256i*)syms+0, ctxV0);
				_mm256_store_si256((__m256i*)syms+1, ctxV1);
				msyms0=_mm256_sub_epi16(myuv[4], moffset0);
				msyms1=_mm256_sub_epi16(myuv[5], moffset1);
				W[4]=msyms0;
				W[5]=msyms1;
				_mm256_store_si256((__m256i*)(rows[0]+0*3*NCODERS*2)+0+4*2, msyms0);
				_mm256_store_si256((__m256i*)(rows[0]+0*3*NCODERS*2)+1+4*2, msyms1);
				++hists[syms[0x00]];
				++hists[syms[0x01]];
				++hists[syms[0x02]];
				++hists[syms[0x03]];
				++hists[syms[0x04]];
				++hists[syms[0x05]];
				++hists[syms[0x06]];
				++hists[syms[0x07]];
				++hists[syms[0x08]];
				++hists[syms[0x09]];
				++hists[syms[0x0A]];
				++hists[syms[0x0B]];
				++hists[syms[0x0C]];
				++hists[syms[0x0D]];
				++hists[syms[0x0E]];
				++hists[syms[0x0F]];
				++hists[syms[0x10]];
				++hists[syms[0x11]];
				++hists[syms[0x12]];
				++hists[syms[0x13]];
				++hists[syms[0x14]];
				++hists[syms[0x15]];
				++hists[syms[0x16]];
				++hists[syms[0x17]];
				++hists[syms[0x18]];
				++hists[syms[0x19]];
				++hists[syms[0x1A]];
				++hists[syms[0x1B]];
				++hists[syms[0x1C]];
				++hists[syms[0x1D]];
				++hists[syms[0x1E]];
				++hists[syms[0x1F]];
			}
			else
			{
				//FIXME maintain W[]
			}
			eNEEE[0]=_mm256_load_si256((__m256i*)(rows[1]+3*3*NCODERS*2)+0+1*2);
			eNEEE[1]=_mm256_load_si256((__m256i*)(rows[1]+3*3*NCODERS*2)+1+1*2);
			eNEEE[2]=_mm256_load_si256((__m256i*)(rows[1]+3*3*NCODERS*2)+0+3*2);
			eNEEE[3]=_mm256_load_si256((__m256i*)(rows[1]+3*3*NCODERS*2)+1+3*2);
			eNEEE[4]=_mm256_load_si256((__m256i*)(rows[1]+3*3*NCODERS*2)+0+5*2);
			eNEEE[5]=_mm256_load_si256((__m256i*)(rows[1]+3*3*NCODERS*2)+1+5*2);
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
			_mm256_store_si256((__m256i*)(rows[0]+0*3*NCODERS*2)+0+1*2, eW[0]);
			_mm256_store_si256((__m256i*)(rows[0]+0*3*NCODERS*2)+1+1*2, eW[1]);
			_mm256_store_si256((__m256i*)(rows[0]+0*3*NCODERS*2)+0+3*2, eW[2]);
			_mm256_store_si256((__m256i*)(rows[0]+0*3*NCODERS*2)+1+3*2, eW[3]);
			_mm256_store_si256((__m256i*)(rows[0]+0*3*NCODERS*2)+0+5*2, eW[4]);
			_mm256_store_si256((__m256i*)(rows[0]+0*3*NCODERS*2)+1+5*2, eW[5]);
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
			rows[0]+=3*NCODERS*2;
			rows[1]+=3*NCODERS*2;
			rows[2]+=3*NCODERS*2;
			rows[3]+=3*NCODERS*2;
			erows[0]+=3*WG_NPREDS*NCODERS;
			erows[1]+=3*WG_NPREDS*NCODERS;
		}
		imptr+=rowstride;
		ctxptr+=rowstride;
	}
	_mm_free(pixels);
	_mm_free(wgerrors);
	if(blockw*NCODERS<iw)
	{
		int sym=0;
		int vfromy=-(combination[II_COEFF_U_SUB_Y]!=0);
		int vc0=combination[II_COEFF_V_SUB_Y];
		int vc1=combination[II_COEFF_V_SUB_U];
		int offset=0;
		unsigned short N[32*3]={0}, yuv[32*3]={0};
		int rembytes=rowstride-3*NCODERS*blockw;
		imptr=image+3*NCODERS*blockw;
		yidx=combination[II_PERM_Y],
		uidx=combination[II_PERM_U],
		vidx=combination[II_PERM_V];
		for(int ky=0;ky<ih;++ky)//remainder coding loop		simple differentiation, static-o0 coding
		{
			for(int kx=0;kx<rembytes;kx+=3)
			{
				//sym = (yuv-N+128)&255
				if(fwd)
				{
					yuv[kx+0]=imptr[kx+yidx];
					yuv[kx+1]=imptr[kx+uidx];
					yuv[kx+2]=imptr[kx+vidx];
					imptr[kx+0]=sym=(unsigned char)(yuv[kx+0]-N[kx+0]+128);
					N[kx+0]=yuv[kx+0];
					++rhist[256*0+sym];

					offset=yuv[kx+0]&vfromy;
					N[kx+1]+=offset;
					CLAMP2(N[kx+1], -128, 127);
					imptr[kx+1]=sym=(unsigned char)(yuv[kx+1]-N[kx+1]+128);
					N[kx+1]=yuv[kx+1]-offset;
					++rhist[256*1+sym];

					offset=vc0*yuv[kx+0]+vc1*yuv[kx+1];
					int vpred=(N[kx+2]+offset)>>2;
					CLAMP2(vpred, -128, 127);
					imptr[kx+2]=sym=(unsigned char)(yuv[kx+2]-vpred+128);
					N[kx+2]=4*yuv[kx+2]-offset;
					++rhist[256*2+sym];
				}
				else
				{
					//FIXME
				}
			}
			imptr+=rowstride;
		}
	}
	if(fwd)
	{
		//normalize/integrate hists
		rANS_SIMD_SymInfo *syminfo=(rANS_SIMD_SymInfo*)CDF2syms;
		rANS_SIMD_SymInfo *rsyminfo=(rANS_SIMD_SymInfo*)rCDF2syms;
		if(blockw*NCODERS<iw)
		{
			enc_hist2stats(rhist+(ptrdiff_t)256*0, rsyminfo+(ptrdiff_t)256*0);
			enc_hist2stats(rhist+(ptrdiff_t)256*1, rsyminfo+(ptrdiff_t)256*1);
			enc_hist2stats(rhist+(ptrdiff_t)256*2, rsyminfo+(ptrdiff_t)256*2);
		}
		for(int kc=0;kc<3*2*(8+GRBITS);++kc)
			enc_hist2stats(hists+(ptrdiff_t)256*kc, syminfo+(ptrdiff_t)256*kc);

		//encode bwd-bwd	remainder (aux lane)
		unsigned state_aux=RANS_BYTE_L;
		{
			rANS_SIMD_SymInfo *info=0;
			int rembytes=rowstride-3*NCODERS*blockw;
			imptr=image+usize-rembytes;
			for(int ky=ih-1;ky>=0;--ky)
			{
				for(int kx=rembytes-3;kx>=0;kx-=3)
				{
					info=rsyminfo+imptr[kx+2]+256*2;
					if(state_aux>info->vmax)
					{
						*(unsigned short*)streamptr=(unsigned short)state_aux;
						streamptr+=2;
						state_aux>>=16;
					}
					state_aux+=((unsigned long long)state_aux*info->invf>>info->sh)*info->negf+info->cdf;

					info=rsyminfo+imptr[kx+1]+256*1;
					if(state_aux>info->vmax)
					{
						*(unsigned short*)streamptr=(unsigned short)state_aux;
						streamptr+=2;
						state_aux>>=16;
					}
					state_aux+=((unsigned long long)state_aux*info->invf>>info->sh)*info->negf+info->cdf;

					info=rsyminfo+imptr[kx+0]+256*0;
					if(state_aux>info->vmax)
					{
						*(unsigned short*)streamptr=(unsigned short)state_aux;
						streamptr+=2;
						state_aux>>=16;
					}
					state_aux+=((unsigned long long)state_aux*info->invf>>info->sh)*info->negf+info->cdf;
				}
				imptr-=rowstride;
			}
		}

		//encode bwd-bwd	quotient
		__m256i mstate[]=
		{
			_mm256_set1_epi32(RANS_BYTE_L),
			_mm256_set1_epi32(RANS_BYTE_L),
			_mm256_set1_epi32(RANS_BYTE_L),
			_mm256_set1_epi32(RANS_BYTE_L),
		};
		imptr=image+usize-rowstride;
		ctxptr=context+usize-rowstride;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=3*blockw-3;kx>=0;kx-=3)
			{
				for(int k=0;k<32;k+=8)
				{
				}
			}
			imptr-=rowstride;
			ctxptr-=rowstride;
		}
		//encode hists zigzag (aux lane)	skip null contexts
		//flush ANS states
		//save compressed file
	}
	else
	{
		//save PPM file
	}
	free(hists);
	free(rhist);
	_mm_free(CDF2syms);
	_mm_free(rCDF2syms);
#if 0
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			LOG_ERROR("Cannot open \"%s\" for writing", fdst);
			return 1;
		}
		if(fwd)
		{
			//FIXME enc bestrct, flush streamptr

			fwrite("29", 1, 2, fdst);
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
	free(image);
	if(fwd)
		free(context);
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		printf("%s    %td bytes\n", srcfn, srcsize);
		printf("%8td/%8td bytes\n", 2+4+4+streamptr-dstbuf, srcsize);
	}
	printf("%c  %12.6lf sec  %12.6lf MB/s\n", 'D'+fwd, t, srcsize/(t*1024*1024));
#endif
#endif
	return 0;
}