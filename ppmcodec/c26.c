#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


	#define LOUD
//	#define ENABLE_GUIDE

	#define ENABLE_ZERO	//good with astro, requires static-o0 to be fast
	#define ENABLE_W	//good with GDCC
	#define ENABLE_NW	//good
	#define ENABLE_AV2	//good, but redundant with WG
	#define ENABLE_CG
	#define ENABLE_SELECT	//good with synth
	#define ENABLE_IZ	//weak with CG/AV4
	#define ENABLE_AV3	//weak with CG/AV4
	#define ENABLE_AV4	//good with GDCC
	#define ENABLE_AV6	//good
	#define ENABLE_AV8	//?		LPF
	#define ENABLE_AV9	//good				GDCC 0.31% slower, 0.37% smaller
//	#define ENABLE_WG	//for noisy areas		GDCC 5.69% slower, 0.07% smaller


#define CBITS 12
#define USTREAM_SIZE 0x2000
#define BLOCKX 768
#define BLOCKY 768

#define ANALYSIS_XSTRIDE 2
#define ANALYSIS_YSTRIDE 2

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
#if 0
#define OCHLIST\
	OCH(R)\
	OCH(G)\
	OCH(B)\
	OCH(RG)\
	OCH(GB)\
	OCH(BR)\
	OCH(R2)\
	OCH(G2)\
	OCH(B2)\
	OCH(RB)\
	OCH(GR)\
	OCH(BG)
typedef enum _OCHIndex
{
#define OCH(LABEL) OCH_##LABEL,
	OCHLIST
#undef  OCH
	OCH_COUNT,
} OCHIndex;
//static const char *och_names[OCH_COUNT]=
//{
//#define OCH(LABEL) #LABEL,
//	OCHLIST
//#undef  OCH
//};
#endif

typedef enum _OCHIndex
{
	OCH_R,
	OCH_G,
	OCH_B,
	OCH_RG,
	OCH_GB,
	OCH_BR,
	OCH_R2,
	OCH_G2,
	OCH_B2,
	OCH_GR=OCH_RG,
	OCH_BG=OCH_GB,
	OCH_RB=OCH_BR,

	OCH_COUNT=9,
} OCHIndex;
typedef enum _RCTInfoIdx
{
	II_OCH_Y,
	II_OCH_U,
	II_OCH_V,

	II_PERM_Y,
	II_PERM_U,
	II_PERM_V,

	II_HELP_U,
	II_HELP_V0,
	II_HELP_V1,

	II_COUNT,
} RCTInfoIdx;
#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 2)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  2, 0)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  2, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 2)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 2)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  2, 0)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	2,  2, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	2,  0, 2)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	2,  0, 2)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	2,  2, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	2,  0, 2)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	2,  0, 2)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	2,  2, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	2,  0, 2)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	2,  0, 2)\
	RCT(R_G_B2,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,  1, 1)\
	RCT(R_GR_B2,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	2,  1, 1)\
	RCT(R_B_G2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,  1, 1)\
	RCT(R_BR_G2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	2,  1, 1)\
	RCT(G_B_R2,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,  1, 1)\
	RCT(G_BG_R2,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	2,  1, 1)\
	RCT(G_RG_B2,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	2,  1, 1)\
	RCT(B_RB_G2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	2,  1, 1)\
	RCT(B_GB_R2,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	2,  1, 1)
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

#define PREDLIST\
	PRED(ZERO)\
	PRED(W)\
	PRED(NW)\
	PRED(AV2)\
	PRED(SELECT)\
	PRED(CG)\
	PRED(IZ)\
	PRED(AV3)\
	PRED(AV4)\
	PRED(AV6)\
	PRED(AV8)\
	PRED(AV9)
typedef enum _PredIndex
{
#define PRED(LABEL) PRED_##LABEL,
	PREDLIST
#undef  PRED
	PRED_COUNT,
} PredIndex;
static const char *pred_names[PRED_COUNT]=
{
#define PRED(LABEL) #LABEL,
	PREDLIST
#undef  PRED
};

//https://github.com/samtools/htscodecs
unsigned char *rans_compress_O0_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_O0_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int out_sz);

int c26_codec(const char *srcfn, const char *dstfn, int nthreads0)
{
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
		srcbuf=(unsigned char*)malloc(srcsize+16);
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
	if(!fwd&&tag!=('2'|'6'<<8))
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
	int xblocks=(iw+BLOCKX-1)/BLOCKX;
	int yblocks=(ih+BLOCKY-1)/BLOCKY;
	int nblocks=xblocks*yblocks;

	if(fwd)
	{
#if 0
		double esizes[3]={0};
#endif
		unsigned char *image=srcptr;
		guide_save(image, iw, ih);

		//analysis
		int *flags=(int*)malloc(sizeof(int)*nblocks);
		int ystride=3*iw;
		int histssize=(int)sizeof(int[PRED_COUNT*OCH_COUNT*256]);
		int *hists=(int*)malloc(histssize);
		//int *stats=(int*)malloc(histssize);
		if(!flags||!hists)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		//memset(stats, 0, histssize);
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		unsigned char predsel[OCH_COUNT]={0}, bestrct=0;
		for(int kb=0;kb<nblocks;++kb)
		{
			int bx=kb%xblocks, by=kb/xblocks;
			int count=0;
			int x1=BLOCKX*bx, x2=x1+BLOCKX;
			int y1=BLOCKY*by, y2=y1+BLOCKY;
			if(x2>iw)
				x2=iw;
			if(y2>ih)
				y2=ih;
			memset(hists, 0, histssize);
			for(int ky=y1+2;ky<y2;ky+=ANALYSIS_YSTRIDE)//analysis loop
			{
				int kx=x1+2;
				const unsigned char *ptr=image+3*(iw*ky+kx);

				__m256i amin=_mm256_set1_epi16(-128);
				__m256i amax=_mm256_set1_epi16(127);
				__m256i amask=_mm256_set1_epi16(255);
				__m128i half8=_mm_set1_epi8(128);
				__m128i shuf=_mm_set_epi8(
					-1,
					12, 14, 13,
					 9, 11, 10,
					 6,  8,  7,
					 3,  5,  4,
					 0,  2,  1
					//15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
				);
				ALIGN(32) short result[16]={0};
				for(;kx<=x2-(5*ANALYSIS_XSTRIDE+1+1);kx+=5*ANALYSIS_XSTRIDE, ptr+=15*ANALYSIS_XSTRIDE, count+=5)
				{
					//		NNW	NN	NNE
					//	NWW	NW	N	NE
					//	WW	W	?
					//1		rgb
					//2		gbr			unused exept for curr2
					//3		rgb - gbr
					//4		gbr - rgb
					//5 curr	(gbr+brg)/2
					//5 ...		rgb - (gbr+brg)/2
					__m256i
					//	NNWW,	NNWW3,	NNWW4,	NNWW5,
						NNW,	NNW3,	NNW4,	NNW5,
						NN,	NN3,	NN4,	NN5,
						NNE,	NNE3,	NNE4,	NNE5,
						NWW,	NWW3,	NWW4,	NWW5,
						NW,	NW3,	NW4,	NW5,
						N,	N3,	N4,	N5,
						NE,	NE3,	NE4,	NE5,
						NEE,	NEE3,	NEE4,	NEE5,
						WW,	WW3,	WW4,	WW5,
						W,	W3,	W4,	W5,
						curr,	curr2,		curr5;
					__m256i vmin[4], vmax[4], pred;
					{
					//	__m128i NNWW8	=_mm_loadu_si128((__m128i*)(ptr-2*ystride-2*3+0));
						__m128i NNW8	=_mm_loadu_si128((__m128i*)(ptr-2*ystride-1*3+0));
						__m128i NN8	=_mm_loadu_si128((__m128i*)(ptr-2*ystride+0*3+0));
						__m128i NNE8	=_mm_loadu_si128((__m128i*)(ptr-2*ystride+1*3+0));
						__m128i NWW8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride-2*3+0));
						__m128i NW8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride-1*3+0));
						__m128i N8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride+0*3+0));
						__m128i NE8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride+1*3+0));
						__m128i NEE8	=_mm_loadu_si128((__m128i*)(ptr-1*ystride+2*3+0));
						__m128i WW8	=_mm_loadu_si128((__m128i*)(ptr+0*ystride-2*3+0));
						__m128i W8	=_mm_loadu_si128((__m128i*)(ptr+0*ystride-1*3+0));
						__m128i curr8	=_mm_loadu_si128((__m128i*)(ptr+0*ystride+0*3+0));
					//	NNWW8	=_mm_xor_si128(NNWW8	, half8);
						NNW8	=_mm_xor_si128(NNW8	, half8);
						NN8	=_mm_xor_si128(NN8	, half8);
						NNE8	=_mm_xor_si128(NNE8	, half8);
						NWW8	=_mm_xor_si128(NWW8	, half8);
						NW8	=_mm_xor_si128(NW8	, half8);
						N8	=_mm_xor_si128(N8	, half8);
						NE8	=_mm_xor_si128(NE8	, half8);
						NEE8	=_mm_xor_si128(NEE8	, half8);
						WW8	=_mm_xor_si128(WW8	, half8);
						W8	=_mm_xor_si128(W8	, half8);
						curr8	=_mm_xor_si128(curr8	, half8);
					//	__m128i NNWW82	=_mm_shuffle_epi8(NNWW8	, shuf);
						__m128i NNW82	=_mm_shuffle_epi8(NNW8	, shuf);
						__m128i NN82	=_mm_shuffle_epi8(NN8	, shuf);
						__m128i NNE82	=_mm_shuffle_epi8(NNE8	, shuf);
						__m128i NWW82	=_mm_shuffle_epi8(NWW8	, shuf);
						__m128i NW82	=_mm_shuffle_epi8(NW8	, shuf);
						__m128i N82	=_mm_shuffle_epi8(N8	, shuf);
						__m128i NE82	=_mm_shuffle_epi8(NE8	, shuf);
						__m128i NEE82	=_mm_shuffle_epi8(NEE8	, shuf);
						__m128i WW82	=_mm_shuffle_epi8(WW8	, shuf);
						__m128i W82	=_mm_shuffle_epi8(W8	, shuf);
						__m128i curr82	=_mm_shuffle_epi8(curr8	, shuf);
					//	NNWW	=_mm256_cvtepi8_epi16(NNWW8);
						NNW	=_mm256_cvtepi8_epi16(NNW8);
						NN	=_mm256_cvtepi8_epi16(NN8);
						NNE	=_mm256_cvtepi8_epi16(NNE8);
						NWW	=_mm256_cvtepi8_epi16(NWW8);
						NW	=_mm256_cvtepi8_epi16(NW8);
						N	=_mm256_cvtepi8_epi16(N8);
						NE	=_mm256_cvtepi8_epi16(NE8);
						NEE	=_mm256_cvtepi8_epi16(NEE8);
						WW	=_mm256_cvtepi8_epi16(WW8);
						W	=_mm256_cvtepi8_epi16(W8);
						curr	=_mm256_cvtepi8_epi16(curr8);
					//	__m256i NNWW2	=_mm256_cvtepi8_epi16(NNWW82);
						__m256i NNW2	=_mm256_cvtepi8_epi16(NNW82);
						__m256i NN2	=_mm256_cvtepi8_epi16(NN82);
						__m256i NNE2	=_mm256_cvtepi8_epi16(NNE82);
						__m256i NWW2	=_mm256_cvtepi8_epi16(NWW82);
						__m256i NW2	=_mm256_cvtepi8_epi16(NW82);
						__m256i N2	=_mm256_cvtepi8_epi16(N82);
						__m256i NE2	=_mm256_cvtepi8_epi16(NE82);
						__m256i NEE2	=_mm256_cvtepi8_epi16(NEE82);
						__m256i WW2	=_mm256_cvtepi8_epi16(WW82);
						__m256i W2	=_mm256_cvtepi8_epi16(W82);
						curr2	=_mm256_cvtepi8_epi16(curr82);
					//	NNWW3	=_mm256_sub_epi16(NNWW	, NNWW2	);
						NNW3	=_mm256_sub_epi16(NNW	, NNW2	);
						NN3	=_mm256_sub_epi16(NN	, NN2	);
						NNE3	=_mm256_sub_epi16(NNE	, NNE2	);
						NWW3	=_mm256_sub_epi16(NWW	, NWW2	);
						NW3	=_mm256_sub_epi16(NW	, NW2	);
						N3	=_mm256_sub_epi16(N	, N2	);
						NE3	=_mm256_sub_epi16(NE	, NE2	);
						NEE3	=_mm256_sub_epi16(NEE	, NEE2	);
						WW3	=_mm256_sub_epi16(WW	, WW2	);
						W3	=_mm256_sub_epi16(W	, W2	);
					//	curr3	=_mm256_sub_epi16(curr	, curr2	);
					//	NNWW4	=_mm256_sub_epi16(NNWW2	, NNWW	);
						NNW4	=_mm256_sub_epi16(NNW2	, NNW	);
						NN4	=_mm256_sub_epi16(NN2	, NN	);
						NNE4	=_mm256_sub_epi16(NNE2	, NNE	);
						NWW4	=_mm256_sub_epi16(NWW2	, NWW	);
						NW4	=_mm256_sub_epi16(NW2	, NW	);
						N4	=_mm256_sub_epi16(N2	, N	);
						NE4	=_mm256_sub_epi16(NE2	, NE	);
						NEE4	=_mm256_sub_epi16(NEE2	, NEE	);
						WW4	=_mm256_sub_epi16(WW2	, WW	);
						W4	=_mm256_sub_epi16(W2	, W	);
					//	curr4	=_mm256_sub_epi16(curr2	, curr	);
					
					//	NNWW5	=_mm256_add_epi16(NNWW2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NNWW82	, shuf)));
						NNW5	=_mm256_add_epi16(NNW2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NNW82	, shuf)));
						NN5	=_mm256_add_epi16(NN2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NN82	, shuf)));
						NNE5	=_mm256_add_epi16(NNE2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NNE82	, shuf)));
						NWW5	=_mm256_add_epi16(NWW2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NWW82	, shuf)));
						NW5	=_mm256_add_epi16(NW2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NW82	, shuf)));
						N5	=_mm256_add_epi16(N2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(N82	, shuf)));
						NE5	=_mm256_add_epi16(NE2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NE82	, shuf)));
						NEE5	=_mm256_add_epi16(NEE2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(NEE82	, shuf)));
						WW5	=_mm256_add_epi16(WW2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(WW82	, shuf)));
						W5	=_mm256_add_epi16(W2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(W82	, shuf)));
						curr5	=_mm256_add_epi16(curr2	, _mm256_cvtepi8_epi16(_mm_shuffle_epi8(curr82	, shuf)));
					//	NNWW5	=_mm256_srai_epi16(NNWW5, 1);
						NNW5	=_mm256_srai_epi16(NNW5	, 1);
						NN5	=_mm256_srai_epi16(NN5	, 1);
						NNE5	=_mm256_srai_epi16(NNE5	, 1);
						NWW5	=_mm256_srai_epi16(NWW5	, 1);
						NW5	=_mm256_srai_epi16(NW5	, 1);
						N5	=_mm256_srai_epi16(N5	, 1);
						NE5	=_mm256_srai_epi16(NE5	, 1);
						NEE5	=_mm256_srai_epi16(NEE5	, 1);
						WW5	=_mm256_srai_epi16(WW5	, 1);
						W5	=_mm256_srai_epi16(W5	, 1);
						curr5	=_mm256_srai_epi16(curr5, 1);
					//	NNWW5	=_mm256_sub_epi16(NNWW	, NNWW5);
						NNW5	=_mm256_sub_epi16(NNW	, NNW5);
						NN5	=_mm256_sub_epi16(NN	, NN5);
						NNE5	=_mm256_sub_epi16(NNE	, NNE5);
						NWW5	=_mm256_sub_epi16(NWW	, NWW5);
						NW5	=_mm256_sub_epi16(NW	, NW5);
						N5	=_mm256_sub_epi16(N	, N5);
						NE5	=_mm256_sub_epi16(NE	, NE5);
						NEE5	=_mm256_sub_epi16(NEE	, NEE5);
						WW5	=_mm256_sub_epi16(WW	, WW5);
						W5	=_mm256_sub_epi16(W	, W5);
					}
	#define UPDATE(PREDIDX, IDX0, IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7, IDX8, IDX9, IDXA, IDXB, IDXC, IDXD, IDXE)\
		do\
		{\
			pred=_mm256_sub_epi16(pred, amin);\
			pred=_mm256_and_si256(pred, amask);\
			_mm256_store_si256((__m256i*)result, pred);\
			++hists[(IDX0*PRED_COUNT+PREDIDX)<<8|result[0x0]];\
			++hists[(IDX1*PRED_COUNT+PREDIDX)<<8|result[0x1]];\
			++hists[(IDX2*PRED_COUNT+PREDIDX)<<8|result[0x2]];\
			++hists[(IDX3*PRED_COUNT+PREDIDX)<<8|result[0x3]];\
			++hists[(IDX4*PRED_COUNT+PREDIDX)<<8|result[0x4]];\
			++hists[(IDX5*PRED_COUNT+PREDIDX)<<8|result[0x5]];\
			++hists[(IDX6*PRED_COUNT+PREDIDX)<<8|result[0x6]];\
			++hists[(IDX7*PRED_COUNT+PREDIDX)<<8|result[0x7]];\
			++hists[(IDX8*PRED_COUNT+PREDIDX)<<8|result[0x8]];\
			++hists[(IDX9*PRED_COUNT+PREDIDX)<<8|result[0x9]];\
			++hists[(IDXA*PRED_COUNT+PREDIDX)<<8|result[0xA]];\
			++hists[(IDXB*PRED_COUNT+PREDIDX)<<8|result[0xB]];\
			++hists[(IDXC*PRED_COUNT+PREDIDX)<<8|result[0xC]];\
			++hists[(IDXD*PRED_COUNT+PREDIDX)<<8|result[0xD]];\
			++hists[(IDXE*PRED_COUNT+PREDIDX)<<8|result[0xE]];\
		}while(0)
					vmin[0]=_mm256_min_epi16(N, W);
					vmax[0]=_mm256_max_epi16(N, W);
					vmin[1]=_mm256_min_epi16(N3, W3);
					vmax[1]=_mm256_max_epi16(N3, W3);
					vmin[2]=_mm256_min_epi16(N4, W4);
					vmax[2]=_mm256_max_epi16(N4, W4);
					vmin[3]=_mm256_min_epi16(N5, W5);
					vmax[3]=_mm256_max_epi16(N5, W5);
				
					//ZERO
#ifdef ENABLE_ZERO
					pred=curr;
					UPDATE(
						PRED_ZERO,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					pred=_mm256_sub_epi16(curr, curr2);
					UPDATE(
						PRED_ZERO,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					pred=_mm256_sub_epi16(curr2, curr);
					UPDATE(
						PRED_ZERO,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					pred=_mm256_sub_epi16(curr, curr5);
					UPDATE(
						PRED_ZERO,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
#endif
				
					//W
#ifdef ENABLE_W
					pred=_mm256_sub_epi16(curr, W);
					UPDATE(
						PRED_W,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					pred=_mm256_add_epi16(W3, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_W,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					pred=_mm256_add_epi16(W4, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_W,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					pred=_mm256_add_epi16(W5, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_W,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
#endif

					//NW
#ifdef ENABLE_NW
					pred=NW;

					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_NW,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					pred=NW3;

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_NW,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					pred=NW4;

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_NW,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					pred=NW5;

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_NW,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
#endif

					//AV2 = (N+W)>>1
#ifdef ENABLE_AV2
					pred=_mm256_srai_epi16(_mm256_add_epi16(N, W), 1);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV2,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					pred=_mm256_srai_epi16(_mm256_add_epi16(N3, W3), 1);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV2,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					pred=_mm256_srai_epi16(_mm256_add_epi16(N4, W4), 1);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_AV2,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					pred=_mm256_srai_epi16(_mm256_add_epi16(N5, W5), 1);

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
#endif

					//SELECT = abs(N-NW)>abs(W-NW) ? N : W
#ifdef ENABLE_SELECT
					pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N, NW)), _mm256_abs_epi16(_mm256_sub_epi16(W, NW)));
					pred=_mm256_blendv_epi8(W, N, pred);

					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_SELECT,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N3, NW3)), _mm256_abs_epi16(_mm256_sub_epi16(W3, NW3)));
					pred=_mm256_blendv_epi8(W3, N3, pred);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_SELECT,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N4, NW4)), _mm256_abs_epi16(_mm256_sub_epi16(W4, NW4)));
					pred=_mm256_blendv_epi8(W4, N4, pred);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_SELECT,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					pred=_mm256_cmpgt_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N5, NW5)), _mm256_abs_epi16(_mm256_sub_epi16(W5, NW5)));
					pred=_mm256_blendv_epi8(W5, N5, pred);

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_SELECT,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
#endif

					//CG = median(N, W, N+W-NW)
#ifdef ENABLE_CG
					pred=_mm256_sub_epi16(_mm256_add_epi16(N, W), NW);
					pred=_mm256_max_epi16(pred, vmin[0]);
					pred=_mm256_min_epi16(pred, vmax[0]);

					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_CG,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					pred=_mm256_sub_epi16(_mm256_add_epi16(N3, W3), NW3);
					pred=_mm256_max_epi16(pred, vmin[1]);
					pred=_mm256_min_epi16(pred, vmax[1]);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_CG,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					pred=_mm256_sub_epi16(_mm256_add_epi16(N4, W4), NW4);
					pred=_mm256_max_epi16(pred, vmin[2]);
					pred=_mm256_min_epi16(pred, vmax[2]);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_CG,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					pred=_mm256_sub_epi16(_mm256_add_epi16(N5, W5), NW5);
					pred=_mm256_max_epi16(pred, vmin[3]);
					pred=_mm256_min_epi16(pred, vmax[3]);

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_CG,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
#endif
					vmin[0]=_mm256_min_epi16(vmin[0], NE);
					vmax[0]=_mm256_max_epi16(vmax[0], NE);
					vmin[1]=_mm256_min_epi16(vmin[1], NE3);
					vmax[1]=_mm256_max_epi16(vmax[1], NE3);
					vmin[2]=_mm256_min_epi16(vmin[2], NE4);
					vmax[2]=_mm256_max_epi16(vmax[2], NE4);
					vmin[3]=_mm256_min_epi16(vmin[3], NE5);
					vmax[3]=_mm256_max_epi16(vmax[3], NE5);

					//IZ = clamp((3*(N+W)-2*NW)>>2, N,W,NE)
#ifdef ENABLE_IZ
					pred=_mm256_add_epi16(N, W);
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(_mm256_sub_epi16(pred, NW), 1));
					pred=_mm256_srai_epi16(pred, 2);
					pred=_mm256_max_epi16(pred, vmin[0]);
					pred=_mm256_min_epi16(pred, vmax[0]);

					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_IZ,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					pred=_mm256_add_epi16(N3, W3);
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(_mm256_sub_epi16(pred, NW3), 1));
					pred=_mm256_srai_epi16(pred, 2);
					pred=_mm256_max_epi16(pred, vmin[1]);
					pred=_mm256_min_epi16(pred, vmax[1]);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_IZ,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					pred=_mm256_add_epi16(N4, W4);
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(_mm256_sub_epi16(pred, NW4), 1));
					pred=_mm256_srai_epi16(pred, 2);
					pred=_mm256_max_epi16(pred, vmin[2]);
					pred=_mm256_min_epi16(pred, vmax[2]);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_IZ,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					pred=_mm256_add_epi16(N5, W5);
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(_mm256_sub_epi16(pred, NW5), 1));
					pred=_mm256_srai_epi16(pred, 2);
					pred=_mm256_max_epi16(pred, vmin[3]);
					pred=_mm256_min_epi16(pred, vmax[3]);

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_IZ,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
#endif

					//AV3 = clamp((5*(N+W)-2*NW)>>3, N,W,NE)
#ifdef ENABLE_AV3
					pred=_mm256_add_epi16(N, W);
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
					pred=_mm256_sub_epi16(pred, _mm256_slli_epi16(NW, 1));
					pred=_mm256_srai_epi16(pred, 3);
					pred=_mm256_max_epi16(pred, vmin[0]);
					pred=_mm256_min_epi16(pred, vmax[0]);

					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV3,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					pred=_mm256_add_epi16(N3, W3);
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
					pred=_mm256_sub_epi16(pred, _mm256_slli_epi16(NW3, 1));
					pred=_mm256_srai_epi16(pred, 3);
					pred=_mm256_max_epi16(pred, vmin[1]);
					pred=_mm256_min_epi16(pred, vmax[1]);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV3,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					pred=_mm256_add_epi16(N4, W4);
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
					pred=_mm256_sub_epi16(pred, _mm256_slli_epi16(NW4, 1));
					pred=_mm256_srai_epi16(pred, 3);
					pred=_mm256_max_epi16(pred, vmin[2]);
					pred=_mm256_min_epi16(pred, vmax[2]);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_AV3,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					pred=_mm256_add_epi16(N5, W5);
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
					pred=_mm256_sub_epi16(pred, _mm256_slli_epi16(NW5, 1));
					pred=_mm256_srai_epi16(pred, 3);
					pred=_mm256_max_epi16(pred, vmin[3]);
					pred=_mm256_min_epi16(pred, vmax[3]);

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV3,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
#endif

					//AV4 = clamp((4*(N+W)+NE-NW)>>3, N,W,NE)
#ifdef ENABLE_AV4
					pred=_mm256_add_epi16(N, W);
					pred=_mm256_slli_epi16(pred, 2);
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE, NW));
					pred=_mm256_srai_epi16(pred, 3);
					pred=_mm256_max_epi16(pred, vmin[0]);
					pred=_mm256_min_epi16(pred, vmax[0]);

					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV4,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					pred=_mm256_add_epi16(N3, W3);
					pred=_mm256_slli_epi16(pred, 2);
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE3, NW3));
					pred=_mm256_srai_epi16(pred, 3);
					pred=_mm256_max_epi16(pred, vmin[1]);
					pred=_mm256_min_epi16(pred, vmax[1]);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV4,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					pred=_mm256_add_epi16(N4, W4);
					pred=_mm256_slli_epi16(pred, 2);
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE4, NW4));
					pred=_mm256_srai_epi16(pred, 3);
					pred=_mm256_max_epi16(pred, vmin[2]);
					pred=_mm256_min_epi16(pred, vmax[2]);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_AV4,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					pred=_mm256_add_epi16(N5, W5);
					pred=_mm256_slli_epi16(pred, 2);
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE5, NW5));
					pred=_mm256_srai_epi16(pred, 3);
					pred=_mm256_max_epi16(pred, vmin[3]);
					pred=_mm256_min_epi16(pred, vmax[3]);

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV4,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
#endif
				
					//AV6
					//			-1
					//		-5	6	1
					//	-1	8	[?]>>3		clamp(N,W,NE)
#ifdef ENABLE_AV6
					pred=_mm256_sub_epi16(N, NW);
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));//5*(N-NW)
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE, WW));//5*(N-NW)+NE-WW
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(N, NN));//6*N-5*NW-NN-WW+NE
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), W);//W+(6*N-5*NW-NN-WW+NE)/8
					pred=_mm256_max_epi16(pred, vmin[0]);
					pred=_mm256_min_epi16(pred, vmax[0]);

					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV6,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					pred=_mm256_sub_epi16(N3, NW3);
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));//5*(N-NW)
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE3, WW3));//5*(N-NW)+NE-WW
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(N3, NN3));//6*N-5*NW-NN-WW+NE
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), W3);//W+(6*N-5*NW-NN-WW+NE)/8
					pred=_mm256_max_epi16(pred, vmin[1]);
					pred=_mm256_min_epi16(pred, vmax[1]);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV6,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					pred=_mm256_sub_epi16(N4, NW4);
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));//5*(N-NW)
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE4, WW4));//5*(N-NW)+NE-WW
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(N4, NN4));//6*N-5*NW-NN-WW+NE
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), W4);//W+(6*N-5*NW-NN-WW+NE)/8
					pred=_mm256_max_epi16(pred, vmin[2]);
					pred=_mm256_min_epi16(pred, vmax[2]);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_AV6,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					pred=_mm256_sub_epi16(N5, NW5);
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));//5*(N-NW)
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NE5, WW5));//5*(N-NW)+NE-WW
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(N5, NN5));//6*N-5*NW-NN-WW+NE
					pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), W5);//W+(6*N-5*NW-NN-WW+NE)/8
					pred=_mm256_max_epi16(pred, vmin[3]);
					pred=_mm256_min_epi16(pred, vmax[3]);

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV6,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
#endif
				
					//AV8
					//			1	1
					//		1	1	1	1
					//	1	1	[?]>>3
#ifdef ENABLE_AV8
					pred=_mm256_add_epi16(N, W);
					pred=_mm256_add_epi16(pred, _mm256_add_epi16(NN, WW));
					pred=_mm256_add_epi16(pred, _mm256_add_epi16(NW, NE));
					pred=_mm256_add_epi16(pred, _mm256_add_epi16(NNE, NEE));
					pred=_mm256_srai_epi16(pred, 3);

					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV8,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					pred=_mm256_add_epi16(N3, W3);
					pred=_mm256_add_epi16(pred, _mm256_add_epi16(NN3, WW3));
					pred=_mm256_add_epi16(pred, _mm256_add_epi16(NW3, NE3));
					pred=_mm256_add_epi16(pred, _mm256_add_epi16(NNE3, NEE3));
					pred=_mm256_srai_epi16(pred, 3);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV8,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					pred=_mm256_add_epi16(N4, W4);
					pred=_mm256_add_epi16(pred, _mm256_add_epi16(NN4, WW4));
					pred=_mm256_add_epi16(pred, _mm256_add_epi16(NW4, NE4));
					pred=_mm256_add_epi16(pred, _mm256_add_epi16(NNE4, NEE4));
					pred=_mm256_srai_epi16(pred, 3);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_AV8,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					pred=_mm256_add_epi16(N5, W5);
					pred=_mm256_add_epi16(pred, _mm256_add_epi16(NN5, WW5));
					pred=_mm256_add_epi16(pred, _mm256_add_epi16(NW5, NE5));
					pred=_mm256_add_epi16(pred, _mm256_add_epi16(NNE5, NEE5));
					pred=_mm256_srai_epi16(pred, 3);

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV8,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
#endif
				
					//AV9
					//		1	-2	-1
					//	-1	-9	10	4
					//	-2	16	[?]>>4		clamp(N,W,NE)
#ifdef ENABLE_AV9
					pred=_mm256_add_epi16(N, _mm256_slli_epi16(N, 2));//5*N
					pred=_mm256_sub_epi16(pred, _mm256_add_epi16(NN, WW));//5*N - (NN+WW)
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(NE, 1));//5*N-NN-WW + 2*NE
					pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(NW, 3), NW));//2*(5*N-NN-WW+2*NE) - 9*NW
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NNW, _mm256_add_epi16(NNE, NWW)));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
					pred=_mm256_add_epi16(W, _mm256_srai_epi16(pred, 4));
					pred=_mm256_max_epi16(pred, vmin[0]);
					pred=_mm256_min_epi16(pred, vmax[0]);

					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV9,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B,
						OCH_R, OCH_G, OCH_B
					);
					pred=_mm256_add_epi16(N3, _mm256_slli_epi16(N3, 2));
					pred=_mm256_sub_epi16(pred, _mm256_add_epi16(NN3, WW3));
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(NE3, 1));
					pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(NW3, 3), NW3));
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NNW3, _mm256_add_epi16(NNE3, NWW3)));
					pred=_mm256_add_epi16(W3, _mm256_srai_epi16(pred, 4));
					pred=_mm256_max_epi16(pred, vmin[1]);
					pred=_mm256_min_epi16(pred, vmax[1]);

					pred=_mm256_add_epi16(pred, curr2);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV9,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR,
						OCH_RG, OCH_GB, OCH_BR
					);
					pred=_mm256_add_epi16(N4, _mm256_slli_epi16(N4, 2));
					pred=_mm256_sub_epi16(pred, _mm256_add_epi16(NN4, WW4));
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(NE4, 1));
					pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(NW4, 3), NW4));
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NNW4, _mm256_add_epi16(NNE4, NWW4)));
					pred=_mm256_add_epi16(W4, _mm256_srai_epi16(pred, 4));
					pred=_mm256_max_epi16(pred, vmin[2]);
					pred=_mm256_min_epi16(pred, vmax[2]);

					pred=_mm256_add_epi16(pred, curr);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr2, pred);
					UPDATE(
						PRED_AV9,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB,
						OCH_GR, OCH_BG, OCH_RB
					);
					pred=_mm256_add_epi16(N5, _mm256_slli_epi16(N5, 2));
					pred=_mm256_sub_epi16(pred, _mm256_add_epi16(NN5, WW5));
					pred=_mm256_add_epi16(pred, _mm256_slli_epi16(NE5, 1));
					pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(NW5, 3), NW5));
					pred=_mm256_add_epi16(pred, _mm256_sub_epi16(NNW5, _mm256_add_epi16(NNE5, NWW5)));
					pred=_mm256_add_epi16(W5, _mm256_srai_epi16(pred, 4));
					pred=_mm256_max_epi16(pred, vmin[3]);
					pred=_mm256_min_epi16(pred, vmax[3]);

					pred=_mm256_add_epi16(pred, curr5);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(curr, pred);
					UPDATE(
						PRED_AV9,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2,
						OCH_R2, OCH_G2, OCH_B2
					);
#endif
				
					//WG
					//gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1
					//gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1
					//pred=(gx*N+gy*W)/(gx+gy)
#ifdef ENABLE_WG
					{
						__m256i gx, gy;
					//	gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W, WW)), 1);
						gx=_mm256_abs_epi16(_mm256_sub_epi16(W, WW));
						gy=_mm256_abs_epi16(_mm256_sub_epi16(W, NW));
						gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N, NW)));
					//	gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N, NN)), 1));
						gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(N, NN)));
						gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE, N)));
						gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE, NNE)));
						gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
						gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
						MIX2_16x16(pred, N, W, gx, gy);

						pred=_mm256_sub_epi16(curr, pred);
						UPDATE(
							PRED_WG,
							OCH_R, OCH_G, OCH_B,
							OCH_R, OCH_G, OCH_B,
							OCH_R, OCH_G, OCH_B,
							OCH_R, OCH_G, OCH_B,
							OCH_R, OCH_G, OCH_B
						);
					//	gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W3, WW3)), 1);
						gx=_mm256_abs_epi16(_mm256_sub_epi16(W3, WW3));
						gy=_mm256_abs_epi16(_mm256_sub_epi16(W3, NW3));
						gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N3, NW3)));
					//	gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N3, NN3)), 1));
						gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(N3, NN3)));
						gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE3, N3)));
						gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE3, NNE3)));
						gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
						gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
						MIX2_16x16(pred, N3, W3, gx, gy);

						pred=_mm256_add_epi16(pred, curr2);
						pred=_mm256_max_epi16(pred, amin);
						pred=_mm256_min_epi16(pred, amax);
						pred=_mm256_sub_epi16(curr, pred);
						UPDATE(
							PRED_WG,
							OCH_RG, OCH_GB, OCH_BR,
							OCH_RG, OCH_GB, OCH_BR,
							OCH_RG, OCH_GB, OCH_BR,
							OCH_RG, OCH_GB, OCH_BR,
							OCH_RG, OCH_GB, OCH_BR
						);
					//	gx=_mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(W4, WW4)), 1);
						gx=_mm256_abs_epi16(_mm256_sub_epi16(W4, WW4));
						gy=_mm256_abs_epi16(_mm256_sub_epi16(W4, NW4));
						gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N4, NW4)));
					//	gy=_mm256_add_epi16(gy, _mm256_slli_epi16(_mm256_abs_epi16(_mm256_sub_epi16(N4, NN4)), 1));
						gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(N4, NN4)));
						gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE4, N4)));
						gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE4, NNE4)));
						gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
						gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
						MIX2_16x16(pred, N4, W4, gx, gy);

						pred=_mm256_add_epi16(pred, curr);
						pred=_mm256_max_epi16(pred, amin);
						pred=_mm256_min_epi16(pred, amax);
						pred=_mm256_sub_epi16(curr2, pred);
						UPDATE(
							PRED_WG,
							OCH_GR, OCH_BG, OCH_RB,
							OCH_GR, OCH_BG, OCH_RB,
							OCH_GR, OCH_BG, OCH_RB,
							OCH_GR, OCH_BG, OCH_RB,
							OCH_GR, OCH_BG, OCH_RB
						);
						gx=_mm256_abs_epi16(_mm256_sub_epi16(W5, WW5));
						gy=_mm256_abs_epi16(_mm256_sub_epi16(W5, NW5));
						gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(N5, NW5)));
						gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(N5, NN5)));
						gx=_mm256_add_epi16(gx, _mm256_abs_epi16(_mm256_sub_epi16(NE5, N5)));
						gy=_mm256_add_epi16(gy, _mm256_abs_epi16(_mm256_sub_epi16(NE5, NNE5)));
						gx=_mm256_add_epi16(gx, _mm256_set1_epi16(1));
						gy=_mm256_add_epi16(gy, _mm256_set1_epi16(1));
						MIX2_16x16(pred, N5, W5, gx, gy);

						pred=_mm256_add_epi16(pred, curr5);
						pred=_mm256_max_epi16(pred, amin);
						pred=_mm256_min_epi16(pred, amax);
						pred=_mm256_sub_epi16(curr, pred);
						UPDATE(
							PRED_WG,
							OCH_R2, OCH_G2, OCH_B2,
							OCH_R2, OCH_G2, OCH_B2,
							OCH_R2, OCH_G2, OCH_B2,
							OCH_R2, OCH_G2, OCH_B2,
							OCH_R2, OCH_G2, OCH_B2
						);
					}
#endif
				}
			}
			if(!count)
			{
				if(bx==xblocks-1&&bx>=1)
					flags[kb]=flags[kb-1];
				else if(by==yblocks-1&&kb>=xblocks)
					flags[kb]=flags[kb-xblocks];
				else//too small to analyze
					flags[kb]=RCT_R_G_B<<24|
						PRED_CG<<0*8|
						PRED_CG<<1*8|
						PRED_CG<<2*8;
				continue;
			}
			double gain=1./count;
			for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)
			{
				int *curr_hist=hists+((size_t)kc<<8);
				double e=0;
				for(int ks=0;ks<256;++ks)
				{
					int freq=curr_hist[ks];
					if(freq)
						e-=freq*log2((double)freq*gain);
				}
				csizes[kc]=e/8;
			}
			for(int kc=0;kc<OCH_COUNT;++kc)//select best predictors
			{
				int bestpred=0;
				for(int kp=1;kp<PRED_COUNT;++kp)
				{
					if(csizes[kc*PRED_COUNT+bestpred]>csizes[kc*PRED_COUNT+kp])
						bestpred=kp;
				}
				predsel[kc]=bestpred;
			}
			for(int kt=0;kt<RCT_COUNT;++kt)//select best RCT
			{
				const unsigned char *group=rct_combinations[kt];
				double csize=
					csizes[group[II_OCH_Y]*PRED_COUNT+predsel[group[II_OCH_Y]]]+
					csizes[group[II_OCH_U]*PRED_COUNT+predsel[group[II_OCH_U]]]+
					csizes[group[II_OCH_V]*PRED_COUNT+predsel[group[II_OCH_V]]];
				if(!kt||bestsize>csize)
					bestsize=csize, bestrct=kt;
			}
			const unsigned char *group=rct_combinations[bestrct];
			unsigned char predidx[]=
			{
				predsel[group[II_OCH_Y]],
				predsel[group[II_OCH_U]],
				predsel[group[II_OCH_V]],
			};
			flags[kb]=bestrct<<24|
				predidx[0]<<0*8|
				predidx[1]<<1*8|
				predidx[2]<<2*8;
#if 0
			esizes[0]+=csizes[group[II_OCH_Y]*PRED_COUNT+predidx[0]];
			esizes[1]+=csizes[group[II_OCH_U]*PRED_COUNT+predidx[1]];
			esizes[2]+=csizes[group[II_OCH_V]*PRED_COUNT+predidx[2]];

			//int ochidx[]=
			//{
			//	group[II_OCH_Y]*PRED_COUNT+predidx[0],
			//	group[II_OCH_U]*PRED_COUNT+predidx[1],
			//	group[II_OCH_V]*PRED_COUNT+predidx[2],
			//};
			//for(int kc=0;kc<3;++kc)
			//{
			//	int *hist=hists+ochidx[kc]*256LL;
			//	int *curr_stats=stats+ochidx[kc]*256LL;
			//	for(int ks=0;ks<256;++ks)
			//		curr_stats[ks]+=hist[ks];
			//}
#endif
		}
#if 0
		//double esize=0;
		//for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)
		//{
		//	int *curr_hist=stats+((size_t)kc<<8);
		//	double e=0;
		//	int sum=0;
		//	for(int ks=0;ks<256;++ks)
		//		sum+=curr_hist[ks];
		//	double gain=1./sum;
		//	for(int ks=0;ks<256;++ks)
		//	{
		//		int freq=curr_hist[ks];
		//		if(freq)
		//			e-=freq*log2((double)freq*gain);
		//	}
		//	esize+=e/8;
		//}
		for(int kb=0;kb<nblocks;++kb)
		{
			int curr_flags=flags[kb];
			printf("block %5d  %-8s %-10s %-10s %-10s\n",
				kb,
				rct_names[curr_flags>>24&255],
				pred_names[curr_flags>>0*8&255],
				pred_names[curr_flags>>1*8&255],
				pred_names[curr_flags>>2*8&255]
			);
		}
		printf("T %12.2lf\n", esizes[0]+esizes[1]+esizes[2]);
		printf("Y %12.2lf\n", esizes[0]);
		printf("U %12.2lf\n", esizes[1]);
		printf("V %12.2lf\n", esizes[2]);
		//printf("E %12.2lf\n", esize);
#endif
		free(hists);

		int nhistograms=0;
		int predhist[3][RCT_COUNT][PRED_COUNT]={0};
		short *blockmap=(short*)malloc(sizeof(int[3])*nblocks);
		if(!blockmap)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(predhist, -1, sizeof(predhist));
		for(int kb=0;kb<nblocks;++kb)
		{
			int flag=flags[kb];
			int rct=flag>>24;
			int preds[]=
			{
				flag>>0*8&255,
				flag>>1*8&255,
				flag>>2*8&255,
			};
			int idx;
			idx=predhist[0][rct][preds[0]];
			if(idx==-1)
				predhist[0][rct][preds[0]]=idx=nhistograms++;
			blockmap[3*kb+0]=idx;

			idx=predhist[1][rct][preds[1]];
			if(idx==-1)
				predhist[1][rct][preds[1]]=idx=nhistograms++;
			blockmap[3*kb+1]=idx;

			idx=predhist[2][rct][preds[2]];
			if(idx==-1)
				predhist[2][rct][preds[2]]=idx=nhistograms++;
			blockmap[3*kb+2]=idx;
		}
		int nstreams=nhistograms<<CBITS;
		unsigned char *ustreams=(unsigned char*)malloc(USTREAM_SIZE*nstreams);
		int *ustreamsizes=(int*)malloc(nstreams*sizeof(int));
		int pbufsize=(int)sizeof(short[4*3])*(iw+16);
		short *pixels=(short*)malloc(pbufsize);
		unsigned cstreamcap=(unsigned)sizeof(char[USTREAM_SIZE*2]);
		unsigned char *cstream=(unsigned char*)malloc(cstreamcap);
		hists=(int*)malloc(nstreams*sizeof(int[256]));//
		if(!ustreams||!ustreamsizes||!pixels||!hists)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(hists, 0, nstreams*sizeof(int[256]));//
		memset(ustreamsizes, 0, nstreams*sizeof(int));
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			LOG_ERROR("Cannot open \"%s\"\n", dstfn);
			return 1;
		}
#if 0
		int nfullstreams=0;
#endif
		fwrite(flags, 4, nblocks, fdst);//FIXME compress the flags
		srcptr=image;
		memset(pixels, 0, pbufsize);
		int paddedwidth=iw+16;
		for(int ky=0;ky<ih;++ky)
		{
			short *rows[]=
			{
				pixels+(paddedwidth*((ky-0LL)&3)+8LL)*3,
				pixels+(paddedwidth*((ky-1LL)&3)+8LL)*3,
				pixels+(paddedwidth*((ky-2LL)&3)+8LL)*3,
				pixels+(paddedwidth*((ky-3LL)&3)+8LL)*3,
			};
			//short regW[3]={0};
			int kb=(iw/BLOCKX)*(ky/BLOCKX);
			int *flagptr=flags+kb;
			short *streamidx=blockmap+3*kb;
			short yuv[3]={0};
			for(int kx=0;kx<iw;++kx, srcptr+=3)
			{
				int flag=*flagptr;
				const unsigned char *combination=rct_combinations[flag>>24];
				char help[]=
				{
					0,
					combination[II_HELP_U],
					combination[II_HELP_V0],

					0,
					0,
					combination[II_HELP_V1],
				};

				yuv[0]=srcptr[combination[II_PERM_Y]];
				yuv[1]=srcptr[combination[II_PERM_U]];
				yuv[2]=srcptr[combination[II_PERM_V]];
#if defined _GNUC__ && !defined PROFILER
#pragma GCC unroll 3
#endif
				for(int kc=0;kc<3;++kc)
				{
					int offset=(help[kc]*yuv[0]+help[kc+3]*yuv[1])>>1;
					int
						NNW	=rows[2][-1*3],
						NN	=rows[2][+0*3],
						NNE	=rows[2][+1*3],
						NWW	=rows[1][-2*3],
						NW	=rows[1][-1*3],
						N	=rows[1][+0*3],
						NE	=rows[1][+1*3],
						NEE	=rows[1][+2*3],
						WW	=rows[0][-2*3],
						W	=rows[0][-1*3];
					int pred0=0;
					int pred=0;
					switch(flag>>kc*8&255)
					{
					case PRED_ZERO:
						pred=0;
						pred0=pred+offset;//
						break;
					case PRED_W:
						pred=W;
						pred0=pred+offset;//
						break;
					case PRED_NW:
						pred=NW;
						pred0=pred+offset;//
						break;
					case PRED_AV2:
						pred0=N+W+2*offset;//
						pred=(N+W)>>1;
						break;
					case PRED_SELECT:
						pred=abs(W-NW)<abs(N-NW)?N:W;
						pred0=pred+offset;//
						break;
					case PRED_CG:
						{
							int vmax=N, vmin=W;
							if(N<W)vmin=N, vmax=W;
							pred=N+W-NW;
							pred0=pred+offset;//
							CLAMP2(pred, vmin, vmax);
						}
						break;
					case PRED_IZ:
						{
							int vmax=N, vmin=W;
							if(N<W)vmin=N, vmax=W;
							if(vmin>NE)vmin=NE;
							if(vmax<NE)vmax=NE;
							pred0=3*(N+W)-2*NW+4*offset;//
							pred=(3*(N+W)-2*NW)>>2;
							CLAMP2(pred, vmin, vmax);
						}
						break;
					case PRED_AV3:
						{
							int vmax=N, vmin=W;
							if(N<W)vmin=N, vmax=W;
							if(vmin>NE)vmin=NE;
							if(vmax<NE)vmax=NE;
							pred0=5*(N+W)-2*NW+8*offset;//
							pred=(5*(N+W)-2*NW)>>3;
							CLAMP2(pred, vmin, vmax);
						}
						break;
					case PRED_AV4:
						{
							int vmax=N, vmin=W;
							if(N<W)vmin=N, vmax=W;
							if(vmin>NE)vmin=NE;
							if(vmax<NE)vmax=NE;
							pred0=4*(N+W)+NE-NW+8*offset;//
							pred=(4*(N+W)+NE-NW)>>3;
							CLAMP2(pred, vmin, vmax);
						}
						break;
					case PRED_AV6:
						{
							int vmax=N, vmin=W;
							if(N<W)vmin=N, vmax=W;
							if(vmin>NE)vmin=NE;
							if(vmax<NE)vmax=NE;
							pred0=16*W+6*N-5*NW-NN-WW+NE+16*offset;//
							pred=W+((6*N-5*NW-NN-WW+NE)>>3);
							CLAMP2(pred, vmin, vmax);
						}
						break;
					case PRED_AV8:
						pred=(NN+NNE+NW+N+NE+NEE+WW+W)>>3;
						break;
					case PRED_AV9:
						{
							int vmax=N, vmin=W;
							if(N<W)vmin=N, vmax=W;
							if(vmin>NE)vmin=NE;
							if(vmax<NE)vmax=NE;
							pred0=16*W+10*N-9*NW+4*NE-2*(NN+WW)+NNW-(NNE+NWW)+16*offset;//
							pred=W+((10*N-9*NW+4*NE-2*(NN+WW)+NNW-(NNE+NWW))>>4);
							CLAMP2(pred, vmin, vmax);
						}
						break;
					}
					pred+=offset;
					CLAMP2(pred, -128, 127);
					int curr=yuv[kc]-128;
					rows[0][0]=curr-offset;
					int error=(curr-pred)&255;
					int ctx=pred0&((1<<CBITS)-1);
					//int ctx=pred>>(8-CBITS)&((1<<CBITS)-1);
					int curridx=streamidx[kc]<<CBITS|ctx;
					int *size=ustreamsizes+curridx;
					unsigned char *ustream=ustreams+USTREAM_SIZE*curridx;
					ustream[(*size)++]=error;
					++hists[curridx<<8|error];//
#if 1
					if(*size>=USTREAM_SIZE)
					{
						cstreamcap=(unsigned)sizeof(char[USTREAM_SIZE*2]);
						rans_compress_O0_32x16_avx2(ustream, USTREAM_SIZE, cstream, &cstreamcap);
						fwrite(&curridx, 1, 2, fdst);		//stream index
						fwrite(&cstreamcap, 1, 4, fdst);	//compressed size
						fwrite(cstream, 1, cstreamcap, fdst);	//compressed stream
						*size=0;
#if 0
						++nfullstreams;
#endif
					}
#endif
					++rows[0];
					++rows[1];
					++rows[2];
					++rows[3];
				}
				if(kx%BLOCKX==BLOCKX-1)
				{
					++kb;
					++flagptr;
					streamidx+=3;
				}
			}
		}
		int nonzerostreams=0;
		double csize=0;
		for(int k=0;k<nstreams;++k)
		{
			int *hist=hists+k*256LL;
			double e=0;
			int sum=0;
			for(int ks=0;ks<256;++ks)
				sum+=hist[ks];
			if(!sum)
				continue;
			++nonzerostreams;
			double norm=1./sum;
			for(int ks=0;ks<256;++ks)
			{
				int freq=hist[ks];
				if(freq)
					e-=freq*log2(freq*norm);
			}
			e/=8;
		//	printf("stream %5d  %12.2lf/%8d\n", k, e, sum);
			csize+=e;
#if 1
			int *size=ustreamsizes+k;
			if(*size)
			{
				unsigned char *ustream=ustreams+USTREAM_SIZE*k;
				cstreamcap=(unsigned)sizeof(char[USTREAM_SIZE*2]);
				rans_compress_O0_32x16_avx2(ustream, *size, cstream, &cstreamcap);
				fwrite(&k, 1, 2, fdst);			//stream index
				fwrite(&cstreamcap, 1, 4, fdst);	//compressed size
				fwrite(cstream, 1, cstreamcap, fdst);	//compressed stream
			}
#endif
		}
	//	printf("csize %12.2lf/%8d  nstreams %8d/%8d\n", csize, 3*iw*ih, nonzerostreams, nstreams);
		printf("%12.2lf\t", csize);
//#if 0
//		printf("Number of full streams  %d\n", nfullstreams);
//#endif
		fclose(fdst);
		free(pixels);
		free(ustreams);
		free(ustreamsizes);
		free(cstream);
		free(hists);
	}
	else
	{
		LOG_ERROR("WIP");
		return 1;
	}

	
	//LOG_ERROR("WIP");

	(void)nthreads0;
	(void)rct_names;
	(void)pred_names;
	return 0;
}