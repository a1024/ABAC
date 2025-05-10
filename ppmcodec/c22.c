#include"util.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
#include<immintrin.h>
static const char file[]=__FILE__;


#ifndef DISABLE_MT
	#define ENABLE_MT
#endif
//	#define ENABLE_GUIDE
//	#define LOUD

	#define ENABLE_FASTANALYSIS	//LPCB 3~14% faster, 0.04% larger than entropy	good
//	#define ENABLE_ENTROPYANALYSIS
	#define USE_RCT8		//LPCB 3~5% faster, 0.2% larger			good

//	#define USE_ROWPRED565	//(6*N+5*(NW+NE)+8)>>4			A70 37.24%
	#define USE_ROWPRED484	//(2*N+NE+NW+2)>>2	LPCB 34.45%	A70 36.98%*	synth 8.14%			exJPEGs 35.91%	best
//	#define USE_ROWPRED3A3	//(10*N+3*(NW+NE)+8)>>4	LPCB 34.38%	A70 37.23%
//	#define USE_ROWPREDN	//N			LPCB 35.04%			synth 6.71%, 3.6~23.3% faster	exJPEGs 31.85%
//	#define USE_PRED31	//(3*N-NN+1)>>1
//	#define USE_PRED21	//2*N-NN	worse than N

	#define USE_O1		//good


#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static void guide_save(unsigned char *image, int iw, int ih)
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
static void guide_check(unsigned char *image, int kx, int ky)
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

#if defined ENABLE_FASTANALYSIS || defined ENABLE_ENTROPYANALYSIS
typedef enum _OCHType
{
	OCH_R,
	OCH_G,
	OCH_B,
	OCH_RG,
	OCH_GB,
	OCH_BR,
	OCH_RB=OCH_BR,
	OCH_GR=OCH_RG,
	OCH_BG=OCH_GB,

	OCH_COUNT=6,
} OCHType;
static const unsigned char rct_indices[][8]=
{//	output channels			permutation	helper index
	{OCH_R,	OCH_G,	OCH_B,		0, 1, 2,	3, 3},
	{OCH_R,	OCH_G,	OCH_BG,		0, 1, 2,	3, 1},
	{OCH_R,	OCH_G,	OCH_BR,		0, 1, 2,	3, 0},
	{OCH_R,	OCH_GR,	OCH_BR,		0, 1, 2,	0, 0},
	{OCH_R,	OCH_GR,	OCH_BG,		0, 1, 2,	0, 1},
	{OCH_R,	OCH_BR,	OCH_GB,		0, 2, 1,	0, 1},
	{OCH_G,	OCH_B,	OCH_RG,		1, 2, 0,	3, 0},
	{OCH_G,	OCH_B,	OCH_RB,		1, 2, 0,	3, 1},
	{OCH_G,	OCH_BG,	OCH_RG,		1, 2, 0,	0, 0},
	{OCH_G,	OCH_BG,	OCH_RB,		1, 2, 0,	0, 1},
	{OCH_G,	OCH_RG,	OCH_BR,		1, 0, 2,	0, 1},
	{OCH_B,	OCH_R,	OCH_GR,		2, 0, 1,	3, 1},
	{OCH_B,	OCH_R,	OCH_GB,		2, 0, 1,	3, 0},
	{OCH_B,	OCH_RB,	OCH_GB,		2, 0, 1,	0, 0},
	{OCH_B,	OCH_RB,	OCH_GR,		2, 0, 1,	0, 1},
	{OCH_B,	OCH_GB,	OCH_RG,		2, 1, 0,	0, 1},
};
#endif
#ifdef ENABLE_ENTROPYANALYSIS
static int hist[6][256]={0};
#endif

//https://github.com/samtools/htscodecs
unsigned char *rans_compress_O0_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_O0_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int out_sz);
unsigned char *rans_compress_O1_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_O1_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int out_sz);
#ifdef USE_O1
#define ENTROPY_ENC rans_compress_O1_32x16_avx2
#define ENTROPY_DEC rans_uncompress_O1_32x16_avx2
#else
#define ENTROPY_ENC rans_compress_O0_32x16_avx2
#define ENTROPY_DEC rans_uncompress_O0_32x16_avx2
#endif

#ifdef ENABLE_MT
typedef struct _ThreadArgs
{
	unsigned char *in, *out, *ret;
	unsigned insize, outsize;
} ThreadArgs;
static void c20_enc(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	args->ret=ENTROPY_ENC(args->in, args->insize, args->out, &args->outsize);
}
static void c20_dec(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	args->ret=ENTROPY_DEC(args->in, args->insize, args->out, args->outsize);
}
#endif
int c22_codec(int argc, char **argv)
{
	if(argc!=3&&argc!=4)
	{
		printf(
			"Usage: \"%s\"  input  output  [maxthreads]    Encode/decode.\n"
			"[maxthreads]:\n"
			"  0: nthreads = number of cores (default)\n"
			"  1: Single thread\n"
			, argv[0]
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];
	int maxthreads=argc<4?0:atoi(argv[3]);
#if defined _MSC_VER && defined LOUD
	double ptime=0, etime=0;
	double elapsed=time_sec();
	unsigned long long cycles=__rdtsc();
	ptrdiff_t csize_actual=0;
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
		if(srcsize<0)
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
	if(srcsize<=2)
	{
		LOG_ERROR("File is empty");
		return 1;
	}
	int tag=*(unsigned short*)srcptr;
	srcptr+=2;
	int fwd=tag==('P'|'6'<<8);
	if(!fwd&&tag!=('C'|'H'<<8))
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
		if(*srcptr++ != ' ')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
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
	
	int psize=(iw+32LL)*sizeof(short[4*3]);//4 padded rows * 3 channels
	short *pixels=(short*)_mm_malloc(psize, sizeof(__m256i));
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(pixels, 0x80, psize);
	//memset(pixels, 0, psize);
	
	if(fwd)//encode
	{
		ptrdiff_t dstbufsize=4LL*iw*ih;
		unsigned short *dstbuf=0, *dstptr=0;
		dstbuf=(unsigned short*)malloc(dstbufsize);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		dstptr=dstbuf;
		unsigned char *image=srcptr;
		guide_save(image, iw, ih);
		unsigned char *dptr=(unsigned char*)dstptr;
		int res=iw*ih;
#if defined _MSC_VER && defined LOUD
		ptime=time_sec();
#endif
#ifdef ENABLE_FASTANALYSIS
		unsigned char bestrct=0;
		{
			__m256i mcounters[OCH_COUNT];//r g b rg gb br rb gr bg
			memset(mcounters, 0, sizeof(mcounters));
			for(int ky=1;ky<ih;++ky)
			{
				int rowstride=iw*3, idx=rowstride*ky;
				__m256i blend0=_mm256_set_epi8(
				//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
					 0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0,
					-1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1
				);
				__m256i blend1=_mm256_set_epi8(
				//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
					-1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,
					 0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0
				);
				__m256i blend2=_mm256_set_epi8(
				//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
					 0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,
					 0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0
				);
				//LSB
				//r00 g00 b00 r01 g01 b01 r02 g02 b02 r03 g03 b03 r04 g04 b04 r05	g05 b05 r06 g06 b06 r07 g07 b07 r08 g08 b08 r09 g09 b09 r10 g10		MSB
				//b10 r11 g11 b11 r12 g12 b12 r13 g13 b13 r14 g14 b14 r15 g15 b15	r16 g16 b16 r17 g17 b17 r18 g18 b18 r19 g19 b19 r20 g20 b20 r21
				//g21 b21 r22 g22 b22 r23 g23 b23 r24 g24 b24 r25 g25 b25 r26 g26	b26 r27 g27 b27 r28 g28 b28 r29 g29 b29 r30 g30 b30 r31 g31 b31
				//
				//  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15	  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
				//r00 r11 r22 r01 r12 r23 r02 r13 r24 r03 r14 r25 r04 r15 r26 r05	r16 r27 r06 r17 r28 r07 r18 r29 r08 r19 r30 r09 r20 r31 r10 r21
				//g21 g00 g11 g22 g01 g12 g23 g02 g13 g24 g03 g14 g25 g04 g15 g26	g05 g16 g27 g06 g17 g28 g07 g18 g29 g08 g19 g30 g09 g20 g31 g10
				//b10 b21 b00 b11 b22 b01 b12 b23 b02 b13 b24 b03 b14 b25 b04 b15	b26 b05 b16 b27 b06 b17 b28 b07 b18 b29 b08 b19 b30 b09 b20 b31
				__m256i maskr=_mm256_set_epi64x(0x0000FFFFFFFFFFFF, -1, 0x0000FFFFFFFFFFFF, -1);
				__m256i shufg=_mm256_set_epi8(
				//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
					-1, -1, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
					-1, -1, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1
				);
				__m256i shufb=_mm256_set_epi8(
				//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
					-1, -1, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,
					-1, -1, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2
				);
				__m256i half=_mm256_set1_epi8(-128);
				//__m256i mask_sh1=_mm256_set1_epi8(127);
				//__m256i mask_sh2=_mm256_set1_epi8(63);
				for(int kx=0;kx<iw-31;kx+=32, idx+=96)
				{
				//	__m256i mNW0	=_mm256_loadu_si256((__m256i*)(image+idx-rowstride+0*32+3));
				//	__m256i mNW1	=_mm256_loadu_si256((__m256i*)(image+idx-rowstride+1*32+3));
				//	__m256i mNW2	=_mm256_loadu_si256((__m256i*)(image+idx-rowstride+2*32+3));
					__m256i mN0	=_mm256_loadu_si256((__m256i*)(image+idx-rowstride+0*32+3));
					__m256i mN1	=_mm256_loadu_si256((__m256i*)(image+idx-rowstride+1*32+3));
					__m256i mN2	=_mm256_loadu_si256((__m256i*)(image+idx-rowstride+2*32+3));
				//	__m256i mNE0	=_mm256_loadu_si256((__m256i*)(image+idx-rowstride+0*32+3));
				//	__m256i mNE1	=_mm256_loadu_si256((__m256i*)(image+idx-rowstride+1*32+3));
				//	__m256i mNE2	=_mm256_loadu_si256((__m256i*)(image+idx-rowstride+2*32+3));
					__m256i mcurr0	=_mm256_loadu_si256((__m256i*)(image+idx+0*32));
					__m256i mcurr1	=_mm256_loadu_si256((__m256i*)(image+idx+1*32));
					__m256i mcurr2	=_mm256_loadu_si256((__m256i*)(image+idx+2*32));
					
				//	__m256i rNW	=_mm256_blendv_epi8(mNW0, mNW1, blend1);
				//	__m256i gNW	=_mm256_blendv_epi8(mNW0, mNW1, blend2);
				//	__m256i bNW	=_mm256_blendv_epi8(mNW0, mNW1, blend0);
					__m256i rN	=_mm256_blendv_epi8(mN0,  mN1,  blend1);
					__m256i gN	=_mm256_blendv_epi8(mN0,  mN1,  blend2);
					__m256i bN	=_mm256_blendv_epi8(mN0,  mN1,  blend0);
				//	__m256i rNE	=_mm256_blendv_epi8(mNE0, mNE1, blend1);
				//	__m256i gNE	=_mm256_blendv_epi8(mNE0, mNE1, blend2);
				//	__m256i bNE	=_mm256_blendv_epi8(mNE0, mNE1, blend0);
					__m256i rcurr	=_mm256_blendv_epi8(mcurr0, mcurr1, blend1);
					__m256i gcurr	=_mm256_blendv_epi8(mcurr0, mcurr1, blend2);
					__m256i bcurr	=_mm256_blendv_epi8(mcurr0, mcurr1, blend0);
				//	rNW	=_mm256_blendv_epi8(rNW, mNW2, blend2);
				//	gNW	=_mm256_blendv_epi8(gNW, mNW2, blend0);
				//	bNW	=_mm256_blendv_epi8(bNW, mNW2, blend1);
					rN	=_mm256_blendv_epi8(rN,  mN2,  blend2);
					gN	=_mm256_blendv_epi8(gN,  mN2,  blend0);
					bN	=_mm256_blendv_epi8(bN,  mN2,  blend1);
				//	rNE	=_mm256_blendv_epi8(rNE, mNE2, blend2);
				//	gNE	=_mm256_blendv_epi8(gNE, mNE2, blend0);
				//	bNE	=_mm256_blendv_epi8(bNE, mNE2, blend1);
					rcurr	=_mm256_blendv_epi8(rcurr, mcurr2, blend2);
					gcurr	=_mm256_blendv_epi8(gcurr, mcurr2, blend0);
					bcurr	=_mm256_blendv_epi8(bcurr, mcurr2, blend1);
					
				//	rNW	=_mm256_and_si256(rNW, maskr);
				//	gNW	=_mm256_shuffle_epi8(gNW, shufg);
				//	bNW	=_mm256_shuffle_epi8(bNW, shufb);
					rN	=_mm256_and_si256(rN, maskr);
					gN	=_mm256_shuffle_epi8(gN, shufg);
					bN	=_mm256_shuffle_epi8(bN, shufb);
				//	rNE	=_mm256_and_si256(rNE, maskr);
				//	gNE	=_mm256_shuffle_epi8(gNE, shufg);
				//	bNE	=_mm256_shuffle_epi8(bNE, shufb);
					rcurr	=_mm256_and_si256(rcurr, maskr);
					gcurr	=_mm256_shuffle_epi8(gcurr, shufg);
					bcurr	=_mm256_shuffle_epi8(bcurr, shufb);

				//	__m256i pr=_mm256_srli_epi16(rN, 1);
				//	__m256i pg=_mm256_srli_epi16(gN, 1);
				//	__m256i pb=_mm256_srli_epi16(bN, 1);
				//	pr=_mm256_and_si256(pr, mask_sh1);
				//	pg=_mm256_and_si256(pg, mask_sh1);
				//	pb=_mm256_and_si256(pb, mask_sh1);
				//	pr=_mm256_add_epi8(pr, _mm256_and_si256(_mm256_srli_epi16(rNW, 2), mask_sh2));
				//	pg=_mm256_add_epi8(pg, _mm256_and_si256(_mm256_srli_epi16(gNW, 2), mask_sh2));
				//	pb=_mm256_add_epi8(pb, _mm256_and_si256(_mm256_srli_epi16(bNW, 2), mask_sh2));
				//	pr=_mm256_add_epi8(pr, _mm256_and_si256(_mm256_srli_epi16(rNE, 2), mask_sh2));
				//	pg=_mm256_add_epi8(pg, _mm256_and_si256(_mm256_srli_epi16(gNE, 2), mask_sh2));
				//	pb=_mm256_add_epi8(pb, _mm256_and_si256(_mm256_srli_epi16(bNE, 2), mask_sh2));

					__m256i dr1=_mm256_add_epi8(_mm256_sub_epi8(rcurr, rN), half);
					__m256i dg1=_mm256_add_epi8(_mm256_sub_epi8(gcurr, gN), half);
					__m256i db1=_mm256_add_epi8(_mm256_sub_epi8(bcurr, bN), half);
#if 1
					mcounters[0]=_mm256_add_epi64(mcounters[0], _mm256_sad_epu8(rcurr, rN));
					mcounters[1]=_mm256_add_epi64(mcounters[1], _mm256_sad_epu8(gcurr, gN));
					mcounters[2]=_mm256_add_epi64(mcounters[2], _mm256_sad_epu8(bcurr, bN));
					mcounters[3]=_mm256_add_epi64(mcounters[3], _mm256_sad_epu8(dr1, dg1));
					mcounters[4]=_mm256_add_epi64(mcounters[4], _mm256_sad_epu8(dg1, db1));
					mcounters[5]=_mm256_add_epi64(mcounters[5], _mm256_sad_epu8(db1, dr1));
#else
					__m256i dr2=_mm256_add_epi8(_mm256_sub_epi8(rcurr, pr), half);
					__m256i dg2=_mm256_add_epi8(_mm256_sub_epi8(gcurr, pg), half);
					__m256i db2=_mm256_add_epi8(_mm256_sub_epi8(bcurr, pb), half);
					mcounters[ 0]=_mm256_add_epi64(mcounters[ 0], _mm256_sad_epu8(rcurr, rN));
					mcounters[ 1]=_mm256_add_epi64(mcounters[ 1], _mm256_sad_epu8(rcurr, pr));
					mcounters[ 2]=_mm256_add_epi64(mcounters[ 2], _mm256_sad_epu8(gcurr, gN));
					mcounters[ 3]=_mm256_add_epi64(mcounters[ 3], _mm256_sad_epu8(gcurr, pg));
					mcounters[ 4]=_mm256_add_epi64(mcounters[ 4], _mm256_sad_epu8(bcurr, bN));
					mcounters[ 5]=_mm256_add_epi64(mcounters[ 5], _mm256_sad_epu8(bcurr, pb));
					mcounters[ 6]=_mm256_add_epi64(mcounters[ 6], _mm256_sad_epu8(dr1, dg1));
					mcounters[ 7]=_mm256_add_epi64(mcounters[ 7], _mm256_sad_epu8(dr2, dg2));
					mcounters[ 8]=_mm256_add_epi64(mcounters[ 8], _mm256_sad_epu8(dg1, db1));
					mcounters[ 9]=_mm256_add_epi64(mcounters[ 9], _mm256_sad_epu8(dg2, db2));
					mcounters[10]=_mm256_add_epi64(mcounters[10], _mm256_sad_epu8(db1, dr1));
					mcounters[11]=_mm256_add_epi64(mcounters[11], _mm256_sad_epu8(db2, dr2));
#endif
				}
			}
			ALIGN(32) unsigned long long counters[OCH_COUNT][4]={0};
			_mm256_store_si256((__m256i*)counters+0, mcounters[0]);
			_mm256_store_si256((__m256i*)counters+1, mcounters[1]);
			_mm256_store_si256((__m256i*)counters+2, mcounters[2]);
			_mm256_store_si256((__m256i*)counters+3, mcounters[3]);
			_mm256_store_si256((__m256i*)counters+4, mcounters[4]);
			_mm256_store_si256((__m256i*)counters+5, mcounters[5]);
			counters[0][0]+=counters[0][1]+counters[0][2]+counters[0][3];
			counters[1][0]+=counters[1][1]+counters[1][2]+counters[1][3];
			counters[2][0]+=counters[2][1]+counters[2][2]+counters[2][3];
			counters[3][0]+=counters[3][1]+counters[3][2]+counters[3][3];
			counters[4][0]+=counters[4][1]+counters[4][2]+counters[4][3];
			counters[5][0]+=counters[5][1]+counters[5][2]+counters[5][3];

			unsigned long long minerr=0;
			for(int kt=0;kt<_countof(rct_indices);++kt)
			{
				const unsigned char *rct=rct_indices[kt];
				unsigned long long currerr=
					+counters[rct[0]][0]
					+counters[rct[1]][0]
					+counters[rct[2]][0]
				;
				if(!kt||minerr>currerr)
				{
					minerr=currerr;
					bestrct=kt;
				}
			}
		}
		const unsigned char *rct=rct_indices[bestrct];
		int yidx=rct[3+0];
		int uidx=rct[3+1];
		int vidx=rct[3+2];
		int uhelpidx=rct[6+0];
		int vhelpidx=rct[6+1];

#elif defined ENABLE_ENTROPYANALYSIS
		unsigned char bestrct=0;
		{
			memset(hist, 0, sizeof(hist));
			for(int ky=1;ky<ih;++ky)
			{
				int rowstride=iw*3, idx=rowstride*ky;
				__m256i blend0=_mm256_set_epi8(
				//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
					 0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0,
					-1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1
				);
				__m256i blend1=_mm256_set_epi8(
				//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
					-1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,
					 0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0
				);
				__m256i blend2=_mm256_set_epi8(
				//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
					 0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,
					 0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0
				);
				//LSB
				//r00 g00 b00 r01 g01 b01 r02 g02 b02 r03 g03 b03 r04 g04 b04 r05	g05 b05 r06 g06 b06 r07 g07 b07 r08 g08 b08 r09 g09 b09 r10 g10		MSB
				//b10 r11 g11 b11 r12 g12 b12 r13 g13 b13 r14 g14 b14 r15 g15 b15	r16 g16 b16 r17 g17 b17 r18 g18 b18 r19 g19 b19 r20 g20 b20 r21
				//g21 b21 r22 g22 b22 r23 g23 b23 r24 g24 b24 r25 g25 b25 r26 g26	b26 r27 g27 b27 r28 g28 b28 r29 g29 b29 r30 g30 b30 r31 g31 b31
				//
				//  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15	  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
				//r00 r11 r22 r01 r12 r23 r02 r13 r24 r03 r14 r25 r04 r15 r26 r05	r16 r27 r06 r17 r28 r07 r18 r29 r08 r19 r30 r09 r20 r31 r10 r21
				//g21 g00 g11 g22 g01 g12 g23 g02 g13 g24 g03 g14 g25 g04 g15 g26	g05 g16 g27 g06 g17 g28 g07 g18 g29 g08 g19 g30 g09 g20 g31 g10
				//b10 b21 b00 b11 b22 b01 b12 b23 b02 b13 b24 b03 b14 b25 b04 b15	b26 b05 b16 b27 b06 b17 b28 b07 b18 b29 b08 b19 b30 b09 b20 b31
				__m256i maskr=_mm256_set_epi64x(0x0000FFFFFFFFFFFF, -1, 0x0000FFFFFFFFFFFF, -1);
				__m256i shufg=_mm256_set_epi8(
				//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
					-1, -1, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
					-1, -1, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1
				);
				__m256i shufb=_mm256_set_epi8(
				//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
					-1, -1, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,
					-1, -1, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2
				);
				ALIGN(32) unsigned char errors[192]={0};
				for(int kx=0;kx<iw-31;kx+=32, idx+=96)
				{
					__m256i mN0	=_mm256_loadu_si256((__m256i*)(image+idx-rowstride+0*32));
					__m256i mN1	=_mm256_loadu_si256((__m256i*)(image+idx-rowstride+1*32));
					__m256i mN2	=_mm256_loadu_si256((__m256i*)(image+idx-rowstride+2*32));
					__m256i mcurr0	=_mm256_loadu_si256((__m256i*)(image+idx+0*32));
					__m256i mcurr1	=_mm256_loadu_si256((__m256i*)(image+idx+1*32));
					__m256i mcurr2	=_mm256_loadu_si256((__m256i*)(image+idx+2*32));

					__m256i rN	=_mm256_blendv_epi8(mN0, mN1, blend1);
					__m256i gN	=_mm256_blendv_epi8(mN0, mN1, blend2);
					__m256i bN	=_mm256_blendv_epi8(mN0, mN1, blend0);
					__m256i rcurr	=_mm256_blendv_epi8(mcurr0, mcurr1, blend1);
					__m256i gcurr	=_mm256_blendv_epi8(mcurr0, mcurr1, blend2);
					__m256i bcurr	=_mm256_blendv_epi8(mcurr0, mcurr1, blend0);
					rN	=_mm256_blendv_epi8(rN, mN2, blend2);
					gN	=_mm256_blendv_epi8(gN, mN2, blend0);
					bN	=_mm256_blendv_epi8(bN, mN2, blend1);
					rcurr	=_mm256_blendv_epi8(rcurr, mcurr2, blend2);
					gcurr	=_mm256_blendv_epi8(gcurr, mcurr2, blend0);
					bcurr	=_mm256_blendv_epi8(bcurr, mcurr2, blend1);

					rN	=_mm256_and_si256(rN, maskr);
					gN	=_mm256_shuffle_epi8(gN, shufg);
					bN	=_mm256_shuffle_epi8(bN, shufb);
					rcurr	=_mm256_and_si256(rcurr, maskr);
					gcurr	=_mm256_shuffle_epi8(gcurr, shufg);
					bcurr	=_mm256_shuffle_epi8(bcurr, shufb);

					__m256i dr=_mm256_sub_epi8(rcurr, rN);
					__m256i dg=_mm256_sub_epi8(gcurr, gN);
					__m256i db=_mm256_sub_epi8(bcurr, bN);
					__m256i drg=_mm256_sub_epi8(dr, dg);
					__m256i dgb=_mm256_sub_epi8(dg, db);
					__m256i dbr=_mm256_sub_epi8(db, dr);

					_mm256_store_si256((__m256i*)errors+0, dr);
					_mm256_store_si256((__m256i*)errors+1, dg);
					_mm256_store_si256((__m256i*)errors+2, db);
					_mm256_store_si256((__m256i*)errors+3, drg);
					_mm256_store_si256((__m256i*)errors+4, dgb);
					_mm256_store_si256((__m256i*)errors+5, dbr);

					++hist[0][errors[32*0+ 0]]; ++hist[1][errors[32*1+ 0]]; ++hist[2][errors[32*2+ 0]]; ++hist[3][errors[32*3+ 0]]; ++hist[4][errors[32*4+ 0]]; ++hist[5][errors[32*5+ 0]];
					++hist[0][errors[32*0+ 1]]; ++hist[1][errors[32*1+ 1]]; ++hist[2][errors[32*2+ 1]]; ++hist[3][errors[32*3+ 1]]; ++hist[4][errors[32*4+ 1]]; ++hist[5][errors[32*5+ 1]];
					++hist[0][errors[32*0+ 2]]; ++hist[1][errors[32*1+ 2]]; ++hist[2][errors[32*2+ 2]]; ++hist[3][errors[32*3+ 2]]; ++hist[4][errors[32*4+ 2]]; ++hist[5][errors[32*5+ 2]];
					++hist[0][errors[32*0+ 3]]; ++hist[1][errors[32*1+ 3]]; ++hist[2][errors[32*2+ 3]]; ++hist[3][errors[32*3+ 3]]; ++hist[4][errors[32*4+ 3]]; ++hist[5][errors[32*5+ 3]];
					++hist[0][errors[32*0+ 4]]; ++hist[1][errors[32*1+ 4]]; ++hist[2][errors[32*2+ 4]]; ++hist[3][errors[32*3+ 4]]; ++hist[4][errors[32*4+ 4]]; ++hist[5][errors[32*5+ 4]];
					++hist[0][errors[32*0+ 5]]; ++hist[1][errors[32*1+ 5]]; ++hist[2][errors[32*2+ 5]]; ++hist[3][errors[32*3+ 5]]; ++hist[4][errors[32*4+ 5]]; ++hist[5][errors[32*5+ 5]];
					++hist[0][errors[32*0+ 6]]; ++hist[1][errors[32*1+ 6]]; ++hist[2][errors[32*2+ 6]]; ++hist[3][errors[32*3+ 6]]; ++hist[4][errors[32*4+ 6]]; ++hist[5][errors[32*5+ 6]];
					++hist[0][errors[32*0+ 7]]; ++hist[1][errors[32*1+ 7]]; ++hist[2][errors[32*2+ 7]]; ++hist[3][errors[32*3+ 7]]; ++hist[4][errors[32*4+ 7]]; ++hist[5][errors[32*5+ 7]];
					++hist[0][errors[32*0+ 8]]; ++hist[1][errors[32*1+ 8]]; ++hist[2][errors[32*2+ 8]]; ++hist[3][errors[32*3+ 8]]; ++hist[4][errors[32*4+ 8]]; ++hist[5][errors[32*5+ 8]];
					++hist[0][errors[32*0+ 9]]; ++hist[1][errors[32*1+ 9]]; ++hist[2][errors[32*2+ 9]]; ++hist[3][errors[32*3+ 9]]; ++hist[4][errors[32*4+ 9]]; ++hist[5][errors[32*5+ 9]];
					++hist[0][errors[32*0+10]]; ++hist[1][errors[32*1+10]]; ++hist[2][errors[32*2+10]]; ++hist[3][errors[32*3+10]]; ++hist[4][errors[32*4+10]]; ++hist[5][errors[32*5+10]];
					++hist[0][errors[32*0+11]]; ++hist[1][errors[32*1+11]]; ++hist[2][errors[32*2+11]]; ++hist[3][errors[32*3+11]]; ++hist[4][errors[32*4+11]]; ++hist[5][errors[32*5+11]];
					++hist[0][errors[32*0+12]]; ++hist[1][errors[32*1+12]]; ++hist[2][errors[32*2+12]]; ++hist[3][errors[32*3+12]]; ++hist[4][errors[32*4+12]]; ++hist[5][errors[32*5+12]];
					++hist[0][errors[32*0+13]]; ++hist[1][errors[32*1+13]]; ++hist[2][errors[32*2+13]]; ++hist[3][errors[32*3+13]]; ++hist[4][errors[32*4+13]]; ++hist[5][errors[32*5+13]];
					++hist[0][errors[32*0+14]]; ++hist[1][errors[32*1+14]]; ++hist[2][errors[32*2+14]]; ++hist[3][errors[32*3+14]]; ++hist[4][errors[32*4+14]]; ++hist[5][errors[32*5+14]];
					++hist[0][errors[32*0+15]]; ++hist[1][errors[32*1+15]]; ++hist[2][errors[32*2+15]]; ++hist[3][errors[32*3+15]]; ++hist[4][errors[32*4+15]]; ++hist[5][errors[32*5+15]];
					++hist[0][errors[32*0+16]]; ++hist[1][errors[32*1+16]]; ++hist[2][errors[32*2+16]]; ++hist[3][errors[32*3+16]]; ++hist[4][errors[32*4+16]]; ++hist[5][errors[32*5+16]];
					++hist[0][errors[32*0+17]]; ++hist[1][errors[32*1+17]]; ++hist[2][errors[32*2+17]]; ++hist[3][errors[32*3+17]]; ++hist[4][errors[32*4+17]]; ++hist[5][errors[32*5+17]];
					++hist[0][errors[32*0+18]]; ++hist[1][errors[32*1+18]]; ++hist[2][errors[32*2+18]]; ++hist[3][errors[32*3+18]]; ++hist[4][errors[32*4+18]]; ++hist[5][errors[32*5+18]];
					++hist[0][errors[32*0+19]]; ++hist[1][errors[32*1+19]]; ++hist[2][errors[32*2+19]]; ++hist[3][errors[32*3+19]]; ++hist[4][errors[32*4+19]]; ++hist[5][errors[32*5+19]];
					++hist[0][errors[32*0+20]]; ++hist[1][errors[32*1+20]]; ++hist[2][errors[32*2+20]]; ++hist[3][errors[32*3+20]]; ++hist[4][errors[32*4+20]]; ++hist[5][errors[32*5+20]];
					++hist[0][errors[32*0+21]]; ++hist[1][errors[32*1+21]]; ++hist[2][errors[32*2+21]]; ++hist[3][errors[32*3+21]]; ++hist[4][errors[32*4+21]]; ++hist[5][errors[32*5+21]];
					++hist[0][errors[32*0+22]]; ++hist[1][errors[32*1+22]]; ++hist[2][errors[32*2+22]]; ++hist[3][errors[32*3+22]]; ++hist[4][errors[32*4+22]]; ++hist[5][errors[32*5+22]];
					++hist[0][errors[32*0+23]]; ++hist[1][errors[32*1+23]]; ++hist[2][errors[32*2+23]]; ++hist[3][errors[32*3+23]]; ++hist[4][errors[32*4+23]]; ++hist[5][errors[32*5+23]];
					++hist[0][errors[32*0+24]]; ++hist[1][errors[32*1+24]]; ++hist[2][errors[32*2+24]]; ++hist[3][errors[32*3+24]]; ++hist[4][errors[32*4+24]]; ++hist[5][errors[32*5+24]];
					++hist[0][errors[32*0+25]]; ++hist[1][errors[32*1+25]]; ++hist[2][errors[32*2+25]]; ++hist[3][errors[32*3+25]]; ++hist[4][errors[32*4+25]]; ++hist[5][errors[32*5+25]];
					++hist[0][errors[32*0+26]]; ++hist[1][errors[32*1+26]]; ++hist[2][errors[32*2+26]]; ++hist[3][errors[32*3+26]]; ++hist[4][errors[32*4+26]]; ++hist[5][errors[32*5+26]];
					++hist[0][errors[32*0+27]]; ++hist[1][errors[32*1+27]]; ++hist[2][errors[32*2+27]]; ++hist[3][errors[32*3+27]]; ++hist[4][errors[32*4+27]]; ++hist[5][errors[32*5+27]];
					++hist[0][errors[32*0+28]]; ++hist[1][errors[32*1+28]]; ++hist[2][errors[32*2+28]]; ++hist[3][errors[32*3+28]]; ++hist[4][errors[32*4+28]]; ++hist[5][errors[32*5+28]];
					++hist[0][errors[32*0+29]]; ++hist[1][errors[32*1+29]]; ++hist[2][errors[32*2+29]]; ++hist[3][errors[32*3+29]]; ++hist[4][errors[32*4+29]]; ++hist[5][errors[32*5+29]];
					++hist[0][errors[32*0+30]]; ++hist[1][errors[32*1+30]]; ++hist[2][errors[32*2+30]]; ++hist[3][errors[32*3+30]]; ++hist[4][errors[32*4+30]]; ++hist[5][errors[32*5+30]];
					++hist[0][errors[32*0+31]]; ++hist[1][errors[32*1+31]]; ++hist[2][errors[32*2+31]]; ++hist[3][errors[32*3+31]]; ++hist[4][errors[32*4+31]]; ++hist[5][errors[32*5+31]];
				}
			}
			double csizes[6]={0};
			for(int kc=0;kc<6;++kc)
			{
				unsigned *curr_hist=(unsigned*)hist[kc];
				double e=0;
				unsigned hsum=0;
				for(int ks=0;ks<256;++ks)
					hsum+=curr_hist[ks];
				for(int ks=0;ks<256;++ks)
				{
					unsigned freq=curr_hist[ks];
					if(freq)
						e-=freq*log2((double)freq/hsum);
				}
				csizes[kc]=e/8;
			}
			double bestsize=0;
			for(int kt=0;kt<_countof(rct_indices);++kt)
			{
				const unsigned char *rct=rct_indices[kt];
				double size=csizes[rct[0]]+csizes[rct[1]]+csizes[rct[2]];
				if(!kt||bestsize>size)
				{
					bestsize=size;
					bestrct=kt;
				}
			}
		}
		const unsigned char *rct=rct_indices[bestrct];
		int yidx=rct[3+0];
		int uidx=rct[3+1];
		int vidx=rct[3+2];
		int uhelpidx=rct[6+0];
		int vhelpidx=rct[6+1];
#endif
		//deinterleave |rgbr gbrg brgb rgbr|gbrg brgb rgbr gbrg|brgb rgbr gbrg brgb| -> |yyyy yyyy yyyy yyyy|uuuu uuuu uuuu uuuu|vvvv vvvv vvvv vvvv| (in order):
		__m128i getr0=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 15, 12,  9,  6,  3,  0
		);
		__m128i getr1=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
			-1, -1, -1, -1, -1, 14, 11,  8,  5,  2, -1, -1, -1, -1, -1, -1
		);
		__m128i getr2=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
			13, 10,  7,  4,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
		);
		__m128i getg0=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 13, 10,  7,  4,  1
		);
		__m128i getg1=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
			-1, -1, -1, -1, -1, 15, 12,  9,  6,  3,  0, -1, -1, -1, -1, -1
		);
		__m128i getg2=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
			14, 11,  8,  5,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
		);
		__m128i getb0=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
			-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 14, 11,  8,  5,  2
		);
		__m128i getb1=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
			-1, -1, -1, -1, -1, -1, 13, 10,  7,  4,  1, -1, -1, -1, -1, -1
		);
		__m128i getb2=_mm_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
			15, 12,  9,  6,  3,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
		);
		__m128i
			gety0, gety1, gety2,
			getu0, getu1, getu2,
			getv0, getv1, getv2;
		switch(bestrct)
		{
		default:
		case  0:gety0=getr0, gety1=getr1, gety2=getr2,  getu0=getg0, getu1=getg1, getu2=getg2,  getv0=getb0, getv1=getb1, getv2=getb2;break;//0, 1, 2
		case  1:gety0=getr0, gety1=getr1, gety2=getr2,  getu0=getg0, getu1=getg1, getu2=getg2,  getv0=getb0, getv1=getb1, getv2=getb2;break;//0, 1, 2
		case  2:gety0=getr0, gety1=getr1, gety2=getr2,  getu0=getg0, getu1=getg1, getu2=getg2,  getv0=getb0, getv1=getb1, getv2=getb2;break;//0, 1, 2
		case  3:gety0=getr0, gety1=getr1, gety2=getr2,  getu0=getg0, getu1=getg1, getu2=getg2,  getv0=getb0, getv1=getb1, getv2=getb2;break;//0, 1, 2
		case  4:gety0=getr0, gety1=getr1, gety2=getr2,  getu0=getg0, getu1=getg1, getu2=getg2,  getv0=getb0, getv1=getb1, getv2=getb2;break;//0, 1, 2
		case  5:gety0=getr0, gety1=getr1, gety2=getr2,  getu0=getb0, getu1=getb1, getu2=getb2,  getv0=getg0, getv1=getg1, getv2=getg2;break;//0, 2, 1
		case  6:gety0=getg0, gety1=getg1, gety2=getg2,  getu0=getb0, getu1=getb1, getu2=getb2,  getv0=getr0, getv1=getr1, getv2=getr2;break;//1, 2, 0
		case  7:gety0=getg0, gety1=getg1, gety2=getg2,  getu0=getb0, getu1=getb1, getu2=getb2,  getv0=getr0, getv1=getr1, getv2=getr2;break;//1, 2, 0
		case  8:gety0=getg0, gety1=getg1, gety2=getg2,  getu0=getb0, getu1=getb1, getu2=getb2,  getv0=getr0, getv1=getr1, getv2=getr2;break;//1, 2, 0
		case  9:gety0=getg0, gety1=getg1, gety2=getg2,  getu0=getb0, getu1=getb1, getu2=getb2,  getv0=getr0, getv1=getr1, getv2=getr2;break;//1, 2, 0
		case 10:gety0=getg0, gety1=getg1, gety2=getg2,  getu0=getr0, getu1=getr1, getu2=getr2,  getv0=getb0, getv1=getb1, getv2=getb2;break;//1, 0, 2
		case 11:gety0=getb0, gety1=getb1, gety2=getb2,  getu0=getr0, getu1=getr1, getu2=getr2,  getv0=getg0, getv1=getg1, getv2=getg2;break;//2, 0, 1
		case 12:gety0=getb0, gety1=getb1, gety2=getb2,  getu0=getr0, getu1=getr1, getu2=getr2,  getv0=getg0, getv1=getg1, getv2=getg2;break;//2, 0, 1
		case 13:gety0=getb0, gety1=getb1, gety2=getb2,  getu0=getr0, getu1=getr1, getu2=getr2,  getv0=getg0, getv1=getg1, getv2=getg2;break;//2, 0, 1
		case 14:gety0=getb0, gety1=getb1, gety2=getb2,  getu0=getr0, getu1=getr1, getu2=getr2,  getv0=getg0, getv1=getg1, getv2=getg2;break;//2, 0, 1
		case 15:gety0=getb0, gety1=getb1, gety2=getb2,  getu0=getg0, getu1=getg1, getu2=getg2,  getv0=getr0, getv1=getr1, getv2=getr2;break;//2, 1, 0
		}
		__m256i helperUmask=_mm256_set1_epi8(-(uhelpidx!=3));
		__m256i helperVsel=_mm256_set1_epi8(-(vhelpidx==1));
		__m256i helperVmask=_mm256_set1_epi8(-(vhelpidx!=3));
		char helperUoffsetscalar=(char)(0x80&-(uhelpidx!=3));
		char helperVoffsetscalar=(char)(0x80&-(vhelpidx!=3));
		__m256i helperUoffset=_mm256_set1_epi8(helperUoffsetscalar);
		__m256i helperVoffset=_mm256_set1_epi8(helperVoffsetscalar);
		__m256i one=_mm256_set1_epi16(1);
		__m256i two=_mm256_set1_epi16(2);
		__m256i eight=_mm256_set1_epi16(8);
		//__m128i mhalf=_mm_set1_epi8(-128);
		__m256i mmin=_mm256_set1_epi16(-128);
		__m256i mmax=_mm256_set1_epi16(127);
		(void)one;
		(void)two;
		(void)eight;
		(void)mmin;
		(void)mmax;
		unsigned char *planes[]=
		{
			dptr+0*res,
			dptr+1*res,
			dptr+2*res,
		};
		for(int ky=0, idx=0, idx2=0;ky<ih;++ky)
		{
			ALIGN(32) unsigned char *rows[]=
			{
				(unsigned char*)(pixels+((iw+32LL)*(((ky-0LL)&3)+0*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-1LL)&3)+0*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-2LL)&3)+0*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-3LL)&3)+0*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-0LL)&3)+1*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-1LL)&3)+1*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-2LL)&3)+1*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-3LL)&3)+1*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-0LL)&3)+2*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-1LL)&3)+2*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-2LL)&3)+2*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-3LL)&3)+2*4)+16LL)),
			};
			int kx=0;
#if 1
			//for(;kx<=iw-16;kx+=16, idx+=3*16, idx2+=16)
			for(;kx<=iw-32;kx+=32, idx+=3*32, idx2+=32)
			{
				//process 32 pixels (96 SPs) at a time
#if defined USE_ROWPRED3A3
				__m256i NW0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]-1));
				__m256i NW1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]-1));
				__m256i NW2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]-1));
				__m256i N0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+0));
				__m256i N1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+0));
				__m256i N2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+0));
				__m256i NE0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+1));
				__m256i NE1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+1));
				__m256i NE2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+1));
				__m256i mp0=_mm256_add_epi16(NW0, NE0);//(10*N+3*(NW+NE)+8)>>4
				__m256i mp1=_mm256_add_epi16(NW1, NE1);
				__m256i mp2=_mm256_add_epi16(NW2, NE2);
				__m256i t0=_mm256_add_epi16(mp0, N0);
				__m256i t1=_mm256_add_epi16(mp1, N1);
				__m256i t2=_mm256_add_epi16(mp2, N2);
				t0	=_mm256_slli_epi16(t0, 1);
				t1	=_mm256_slli_epi16(t1, 1);
				t2	=_mm256_slli_epi16(t2, 1);
				mp0	=_mm256_add_epi16(mp0, t0);
				mp1	=_mm256_add_epi16(mp1, t1);
				mp2	=_mm256_add_epi16(mp2, t2);
				t0	=_mm256_slli_epi16(N0, 3);
				t1	=_mm256_slli_epi16(N1, 3);
				t2	=_mm256_slli_epi16(N2, 3);
				mp0	=_mm256_add_epi16(mp0, t0);
				mp1	=_mm256_add_epi16(mp1, t1);
				mp2	=_mm256_add_epi16(mp2, t2);
				mp0	=_mm256_add_epi16(mp0, eight);
				mp1	=_mm256_add_epi16(mp1, eight);
				mp2	=_mm256_add_epi16(mp2, eight);
				mp0	=_mm256_srai_epi16(mp0, 4);
				mp1	=_mm256_srai_epi16(mp1, 4);
				mp2	=_mm256_srai_epi16(mp2, 4);
				__m128i ypred=_mm_packs_epi16(_mm256_castsi256_si128(mp0), _mm256_extracti128_si256(mp0, 1));
				__m128i upred=_mm_packs_epi16(_mm256_castsi256_si128(mp1), _mm256_extracti128_si256(mp1, 1));
				__m128i vpred=_mm_packs_epi16(_mm256_castsi256_si128(mp2), _mm256_extracti128_si256(mp2, 1));
#elif defined USE_ROWPRED484
				__m256i NW0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]-1));
				__m256i NW1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]-1));
				__m256i NW2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]-1));
				__m256i N0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+0));
				__m256i N1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+0));
				__m256i N2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+0));
				__m256i NE0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+1));
				__m256i NE1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+1));
				__m256i NE2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+1));
				NW0=_mm256_avg_epu8(NW0, NE0);
				NW1=_mm256_avg_epu8(NW1, NE1);
				NW2=_mm256_avg_epu8(NW2, NE2);
				N0=_mm256_avg_epu8(N0, NW0);
				N1=_mm256_avg_epu8(N1, NW1);
				N2=_mm256_avg_epu8(N2, NW2);
				//__m256i mp0=_mm256_slli_epi16(N0, 1);//(2*N+NW+NE+2)>>2
				//__m256i mp1=_mm256_slli_epi16(N1, 1);
				//__m256i mp2=_mm256_slli_epi16(N2, 1);
				//mp0=_mm256_add_epi16(mp0, NW0);
				//mp1=_mm256_add_epi16(mp1, NW1);
				//mp2=_mm256_add_epi16(mp2, NW2);
				//mp0=_mm256_add_epi16(mp0, NE0);
				//mp1=_mm256_add_epi16(mp1, NE1);
				//mp2=_mm256_add_epi16(mp2, NE2);
				//mp0=_mm256_add_epi16(mp0, two);
				//mp1=_mm256_add_epi16(mp1, two);
				//mp2=_mm256_add_epi16(mp2, two);
				//mp0=_mm256_srai_epi16(mp0, 2);
				//mp1=_mm256_srai_epi16(mp1, 2);
				//mp2=_mm256_srai_epi16(mp2, 2);
				//__m128i ypred=_mm_packs_epi16(_mm256_castsi256_si128(mp0), _mm256_extracti128_si256(mp0, 1));
				//__m128i upred=_mm_packs_epi16(_mm256_castsi256_si128(mp1), _mm256_extracti128_si256(mp1, 1));
				//__m128i vpred=_mm_packs_epi16(_mm256_castsi256_si128(mp2), _mm256_extracti128_si256(mp2, 1));
#elif defined USE_ROWPRED565
				__m256i NW0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]-1));
				__m256i NW1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]-1));
				__m256i NW2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]-1));
				__m256i N0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+0));
				__m256i N1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+0));
				__m256i N2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+0));
				__m256i NE0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+1));
				__m256i NE1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+1));
				__m256i NE2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+1));
				__m256i mp0=_mm256_add_epi16(NW0, NE0);//(6*N+5*(NW+NE)+8)>>4
				__m256i mp1=_mm256_add_epi16(NW1, NE1);
				__m256i mp2=_mm256_add_epi16(NW2, NE2);
				__m256i t0=_mm256_add_epi16(mp0, N0);
				__m256i t1=_mm256_add_epi16(mp1, N1);
				__m256i t2=_mm256_add_epi16(mp2, N2);
				t0	=_mm256_slli_epi16(t0, 2);
				t1	=_mm256_slli_epi16(t1, 2);
				t2	=_mm256_slli_epi16(t2, 2);
				mp0	=_mm256_add_epi16(mp0, t0);
				mp1	=_mm256_add_epi16(mp1, t1);
				mp2	=_mm256_add_epi16(mp2, t2);
				t0	=_mm256_slli_epi16(N0, 1);
				t1	=_mm256_slli_epi16(N1, 1);
				t2	=_mm256_slli_epi16(N2, 1);
				mp0	=_mm256_add_epi16(mp0, t0);
				mp1	=_mm256_add_epi16(mp1, t1);
				mp2	=_mm256_add_epi16(mp2, t2);
				mp0	=_mm256_add_epi16(mp0, eight);
				mp1	=_mm256_add_epi16(mp1, eight);
				mp2	=_mm256_add_epi16(mp2, eight);
				mp0	=_mm256_srai_epi16(mp0, 4);
				mp1	=_mm256_srai_epi16(mp1, 4);
				mp2	=_mm256_srai_epi16(mp2, 4);
				__m128i ypred=_mm_packs_epi16(_mm256_castsi256_si128(mp0), _mm256_extracti128_si256(mp0, 1));
				__m128i upred=_mm_packs_epi16(_mm256_castsi256_si128(mp1), _mm256_extracti128_si256(mp1, 1));
				__m128i vpred=_mm_packs_epi16(_mm256_castsi256_si128(mp2), _mm256_extracti128_si256(mp2, 1));
#elif defined USE_PRED31
				__m256i NN0	=_mm256_loadu_si256((__m256i*)(rows[2+0*4]+0));
				__m256i NN1	=_mm256_loadu_si256((__m256i*)(rows[2+1*4]+0));
				__m256i NN2	=_mm256_loadu_si256((__m256i*)(rows[2+2*4]+0));
				__m256i N0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+0));
				__m256i N1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+0));
				__m256i N2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+0));
				__m256i mp0=_mm256_slli_epi16(N0, 1);//clamp(3*N-NN+1)>>1
				__m256i mp1=_mm256_slli_epi16(N1, 1);
				__m256i mp2=_mm256_slli_epi16(N2, 1);
				mp0=_mm256_add_epi16(mp0, N0);
				mp1=_mm256_add_epi16(mp1, N1);
				mp2=_mm256_add_epi16(mp2, N2);
				mp0=_mm256_sub_epi16(mp0, NN0);
				mp1=_mm256_sub_epi16(mp1, NN1);
				mp2=_mm256_sub_epi16(mp2, NN2);
				mp0=_mm256_add_epi16(mp0, one);
				mp1=_mm256_add_epi16(mp1, one);
				mp2=_mm256_add_epi16(mp2, one);
				mp0=_mm256_srai_epi16(mp0, 1);
				mp1=_mm256_srai_epi16(mp1, 1);
				mp2=_mm256_srai_epi16(mp2, 1);
				mp0=_mm256_max_epi16(mp0, mmin);
				mp1=_mm256_max_epi16(mp1, mmin);
				mp2=_mm256_max_epi16(mp2, mmin);
				mp0=_mm256_min_epi16(mp0, mmax);
				mp1=_mm256_min_epi16(mp1, mmax);
				mp2=_mm256_min_epi16(mp2, mmax);
				__m128i ypred=_mm_packs_epi16(_mm256_castsi256_si128(mp0), _mm256_extracti128_si256(mp0, 1));
				__m128i upred=_mm_packs_epi16(_mm256_castsi256_si128(mp1), _mm256_extracti128_si256(mp1, 1));
				__m128i vpred=_mm_packs_epi16(_mm256_castsi256_si128(mp2), _mm256_extracti128_si256(mp2, 1));
#elif defined USE_PRED21
				__m256i NN0	=_mm256_loadu_si256((__m256i*)(rows[2+0*4]+0));
				__m256i NN1	=_mm256_loadu_si256((__m256i*)(rows[2+1*4]+0));
				__m256i NN2	=_mm256_loadu_si256((__m256i*)(rows[2+2*4]+0));
				__m256i N0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+0));
				__m256i N1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+0));
				__m256i N2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+0));
				__m256i mp0=_mm256_slli_epi16(N0, 1);//clamp(2*N-NN)
				__m256i mp1=_mm256_slli_epi16(N1, 1);
				__m256i mp2=_mm256_slli_epi16(N2, 1);
				mp0=_mm256_sub_epi16(mp0, NN0);
				mp1=_mm256_sub_epi16(mp1, NN1);
				mp2=_mm256_sub_epi16(mp2, NN2);
				mp0=_mm256_max_epi16(mp0, mmin);
				mp1=_mm256_max_epi16(mp1, mmin);
				mp2=_mm256_max_epi16(mp2, mmin);
				mp0=_mm256_min_epi16(mp0, mmax);
				mp1=_mm256_min_epi16(mp1, mmax);
				mp2=_mm256_min_epi16(mp2, mmax);
				__m128i ypred=_mm_packs_epi16(_mm256_castsi256_si128(mp0), _mm256_extracti128_si256(mp0, 1));
				__m128i upred=_mm_packs_epi16(_mm256_castsi256_si128(mp1), _mm256_extracti128_si256(mp1, 1));
				__m128i vpred=_mm_packs_epi16(_mm256_castsi256_si128(mp2), _mm256_extracti128_si256(mp2, 1));
#elif defined USE_ROWPREDN
				__m128i ypred=_mm_loadu_si128((__m128i*)(rows[1+0*4]+0));
				__m128i upred=_mm_loadu_si128((__m128i*)(rows[1+1*4]+0));
				__m128i vpred=_mm_loadu_si128((__m128i*)(rows[1+2*4]+0));
#endif

				__m128i p0=_mm_loadu_si128((__m128i*)(image+idx+0*16));
				__m128i p1=_mm_loadu_si128((__m128i*)(image+idx+1*16));
				__m128i p2=_mm_loadu_si128((__m128i*)(image+idx+2*16));
				__m128i p3=_mm_loadu_si128((__m128i*)(image+idx+3*16));
				__m128i p4=_mm_loadu_si128((__m128i*)(image+idx+4*16));
				__m128i p5=_mm_loadu_si128((__m128i*)(image+idx+5*16));

				//p0=_mm_sub_epi8(p0, mhalf);
				//p1=_mm_sub_epi8(p1, mhalf);
				//p2=_mm_sub_epi8(p2, mhalf);

				//deinterleave |rgbr gbrg brgb rgbr|gbrg brgb rgbr gbrg|brgb rgbr gbrg brgb| -> |yyyy yyyy yyyy yyyy|uuuu uuuu uuuu uuuu|vvvv vvvv vvvv vvvv| (in order):
				__m128i y0=_mm_shuffle_epi8(p0, gety0);
				__m128i u0=_mm_shuffle_epi8(p0, getu0);
				__m128i v0=_mm_shuffle_epi8(p0, getv0);
				__m128i y3=_mm_shuffle_epi8(p3, gety0);
				__m128i u3=_mm_shuffle_epi8(p3, getu0);
				__m128i v3=_mm_shuffle_epi8(p3, getv0);
				y0=_mm_or_si128(y0, _mm_shuffle_epi8(p1, gety1));
				u0=_mm_or_si128(u0, _mm_shuffle_epi8(p1, getu1));
				v0=_mm_or_si128(v0, _mm_shuffle_epi8(p1, getv1));
				y3=_mm_or_si128(y3, _mm_shuffle_epi8(p4, gety1));
				u3=_mm_or_si128(u3, _mm_shuffle_epi8(p4, getu1));
				v3=_mm_or_si128(v3, _mm_shuffle_epi8(p4, getv1));
				y0=_mm_or_si128(y0, _mm_shuffle_epi8(p2, gety2));
				u0=_mm_or_si128(u0, _mm_shuffle_epi8(p2, getu2));
				v0=_mm_or_si128(v0, _mm_shuffle_epi8(p2, getv2));
				y3=_mm_or_si128(y3, _mm_shuffle_epi8(p5, gety2));
				u3=_mm_or_si128(u3, _mm_shuffle_epi8(p5, getu2));
				v3=_mm_or_si128(v3, _mm_shuffle_epi8(p5, getv2));
				__m256i my=_mm256_inserti128_si256(_mm256_castsi128_si256(y0), y3, 1);
				__m256i mu=_mm256_inserti128_si256(_mm256_castsi128_si256(u0), u3, 1);
				__m256i mv=_mm256_inserti128_si256(_mm256_castsi128_si256(v0), v3, 1);
				__m256i helperU0=_mm256_and_si256(my, helperUmask);
				__m256i helperV0=_mm256_blendv_epi8(my, mu, helperVsel);
				helperV0=_mm256_and_si256(helperV0, helperVmask);
				mu=_mm256_sub_epi8(mu, helperU0);
				mv=_mm256_sub_epi8(mv, helperV0);
				mu=_mm256_sub_epi8(mu, helperUoffset);
				mv=_mm256_sub_epi8(mv, helperVoffset);
				_mm256_storeu_si256((__m256i*)rows[0+0*4], my);
				_mm256_storeu_si256((__m256i*)rows[0+1*4], mu);
				_mm256_storeu_si256((__m256i*)rows[0+2*4], mv);
//#ifdef USE_ROWPREDN
//				_mm_storeu_si128((__m128i*)rows[0+0*4], my);
//				_mm_storeu_si128((__m128i*)rows[0+1*4], mu);
//				_mm_storeu_si128((__m128i*)rows[0+2*4], mv);
//#else
//				_mm256_storeu_si256((__m256i*)rows[0+0*4], _mm256_cvtepi8_epi16(my));
//				_mm256_storeu_si256((__m256i*)rows[0+1*4], _mm256_cvtepi8_epi16(mu));
//				_mm256_storeu_si256((__m256i*)rows[0+2*4], _mm256_cvtepi8_epi16(mv));
//#endif
				
				my=_mm256_sub_epi8(my, N0);
				mu=_mm256_sub_epi8(mu, N1);
				mv=_mm256_sub_epi8(mv, N2);
				//my=_mm_sub_epi8(my, ypred);
				//mu=_mm_sub_epi8(mu, upred);
				//mv=_mm_sub_epi8(mv, vpred);
				
				_mm256_storeu_si256((__m256i*)(planes[0]+idx2), my);
				_mm256_storeu_si256((__m256i*)(planes[1]+idx2), mu);
				_mm256_storeu_si256((__m256i*)(planes[2]+idx2), mv);
				//_mm_storeu_si128((__m128i*)(planes[0]+idx2), my);
				//_mm_storeu_si128((__m128i*)(planes[1]+idx2), mu);
				//_mm_storeu_si128((__m128i*)(planes[2]+idx2), mv);

				__m256i rows0=_mm256_load_si256((__m256i*)rows+0);
				__m256i rows1=_mm256_load_si256((__m256i*)rows+1);
				__m256i rows2=_mm256_load_si256((__m256i*)rows+2);
				rows0=_mm256_add_epi64(rows0, _mm256_set1_epi64x(sizeof(char[32])));
				rows1=_mm256_add_epi64(rows1, _mm256_set1_epi64x(sizeof(char[32])));
				rows2=_mm256_add_epi64(rows2, _mm256_set1_epi64x(sizeof(char[32])));
				_mm256_store_si256((__m256i*)rows+0, rows0);
				_mm256_store_si256((__m256i*)rows+1, rows1);
				_mm256_store_si256((__m256i*)rows+2, rows2);
			}
#endif
//#ifdef __GNUC__
//#pragma GCC unroll 2
//#endif
			for(;kx<iw;++kx, idx+=3, ++idx2)
			{
				//if(ky==0&&kx==6&&kc==0)//
				//if(ky==399&&kx==714)//
				//if(kx==4&&ky==1)//
				//if(ky==2&&kx==4)//
				//if(ky==355&&kx==647)//
				//if(ky==3&&kx==104)//
				//if(ky==3&&kx==103)//
				//if(ky==4&&kx==138)//
				//	printf("");
				
				//short *pcurr=rows[0]+0*4;
				unsigned char yuv[]=
				{
					image[idx+yidx],
					image[idx+uidx],
					image[idx+vidx],
					//image[idx+yidx]-128,
					//image[idx+uidx]-128,
					//image[idx+vidx]-128,
					0,
				};
#if defined USE_ROWPRED3A3
				int ypred=(10*rows[1+0*4][+0]+3*(rows[1+0*4][-1]+rows[1+0*4][+1])+8)>>4;
				int upred=(10*rows[1+1*4][+0]+3*(rows[1+1*4][-1]+rows[1+1*4][+1])+8)>>4;
				int vpred=(10*rows[1+2*4][+0]+3*(rows[1+2*4][-1]+rows[1+2*4][+1])+8)>>4;
#elif defined USE_ROWPRED484
				int ypred=(rows[1+0*4][+0]+((rows[1+0*4][-1]+rows[1+0*4][+1]+1)>>1)+1)>>1;//6 steps (2 roundings)
				int upred=(rows[1+1*4][+0]+((rows[1+1*4][-1]+rows[1+1*4][+1]+1)>>1)+1)>>1;
				int vpred=(rows[1+2*4][+0]+((rows[1+2*4][-1]+rows[1+2*4][+1]+1)>>1)+1)>>1;
				//int ypred=(2*rows[1+0*4][+0]+rows[1+0*4][-1]+rows[1+0*4][+1]+2)>>2;//5 steps
				//int upred=(2*rows[1+1*4][+0]+rows[1+1*4][-1]+rows[1+1*4][+1]+2)>>2;
				//int vpred=(2*rows[1+2*4][+0]+rows[1+2*4][-1]+rows[1+2*4][+1]+2)>>2;
#elif defined USE_ROWPRED565
				int ypred=(6*rows[1+0*4][+0]+5*(rows[1+0*4][-1]+rows[1+0*4][+1])+8)>>4;
				int upred=(6*rows[1+1*4][+0]+5*(rows[1+1*4][-1]+rows[1+1*4][+1])+8)>>4;
				int vpred=(6*rows[1+2*4][+0]+5*(rows[1+2*4][-1]+rows[1+2*4][+1])+8)>>4;
#elif defined USE_PRED31
				int ypred=(3*rows[1+0*4][+0]-rows[2+0*4][+0]+1)>>1;
				int upred=(3*rows[1+1*4][+0]-rows[2+1*4][+0]+1)>>1;
				int vpred=(3*rows[1+2*4][+0]-rows[2+2*4][+0]+1)>>1;
				CLAMP2(ypred, -128, 127);
				CLAMP2(upred, -128, 127);
				CLAMP2(vpred, -128, 127);
#elif defined USE_PRED21
				int ypred=2*rows[1+0*4][+0]-rows[2+0*4][+0];
				int upred=2*rows[1+1*4][+0]-rows[2+1*4][+0];
				int vpred=2*rows[1+2*4][+0]-rows[2+2*4][+0];
				CLAMP2(ypred, -128, 127);
				CLAMP2(upred, -128, 127);
				CLAMP2(vpred, -128, 127);
#elif defined USE_ROWPREDN
				int ypred=rows[1+0*4][+0];
				int upred=rows[1+1*4][+0];
				int vpred=rows[1+2*4][+0];
#endif

				//if(ky==1&&kx==192+12)//
				//if(ky==2&&kx==4)//
				//	printf("");
				yuv[2]-=yuv[vhelpidx]-helperVoffsetscalar;
				yuv[1]-=yuv[uhelpidx]-helperUoffsetscalar;
				rows[0+0*4][+0]=yuv[0];
				rows[0+1*4][+0]=yuv[1];
				rows[0+2*4][+0]=yuv[2];

				dptr[res*0+idx2]=yuv[0]-ypred;
				dptr[res*1+idx2]=yuv[1]-upred;
				dptr[res*2+idx2]=yuv[2]-vpred;

				__m256i rows0=_mm256_load_si256((__m256i*)rows+0);
				__m256i rows1=_mm256_load_si256((__m256i*)rows+1);
				__m256i rows2=_mm256_load_si256((__m256i*)rows+2);
				rows0=_mm256_add_epi64(rows0, _mm256_set1_epi64x(sizeof(char)));
				rows1=_mm256_add_epi64(rows1, _mm256_set1_epi64x(sizeof(char)));
				rows2=_mm256_add_epi64(rows2, _mm256_set1_epi64x(sizeof(char)));
				_mm256_store_si256((__m256i*)rows+0, rows0);
				_mm256_store_si256((__m256i*)rows+1, rows1);
				_mm256_store_si256((__m256i*)rows+2, rows2);
			}
		}
		unsigned char *cdata[3]={0};
		unsigned cdatasize[3]={0};
#if defined _MSC_VER && defined LOUD
		ptime=time_sec()-ptime;
		etime=time_sec();
#endif
#ifdef ENABLE_MT
		if(maxthreads>1)
		{
			ThreadArgs args[]=
			{
				{dptr+res*0, 0, 0, res, 0},
				{dptr+res*1, 0, 0, res, 0},
				{dptr+res*2, 0, 0, res, 0},
			};
			void *mt=mt_exec(c20_enc, args, sizeof(*args), _countof(args));
			mt_finish(mt);
			cdata[0]=args[0].ret;	cdatasize[0]=args[0].outsize;
			cdata[1]=args[1].ret;	cdatasize[1]=args[1].outsize;
			cdata[2]=args[2].ret;	cdatasize[2]=args[2].outsize;
		}
		else
#endif
		{
			cdata[0]=ENTROPY_ENC(dptr+res*0, res, cdata[0], cdatasize+0);
			cdata[1]=ENTROPY_ENC(dptr+res*1, res, cdata[1], cdatasize+1);
			cdata[2]=ENTROPY_ENC(dptr+res*2, res, cdata[2], cdatasize+2);
		}
		if(!cdata[0]||!cdata[1]||!cdata[2])
		{
			LOG_ERROR("HTS Encode Error");
			return 1;
		}
#if defined _MSC_VER && defined LOUD
		etime=time_sec()-etime;
#endif
		{
			ptrdiff_t csize=0;
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
				return 1;
			}
			csize+=fwrite("CH", 1, 2, fdst);
			csize+=fwrite(&iw, 1, 4, fdst);
			csize+=fwrite(&ih, 1, 4, fdst);
			
#if defined ENABLE_FASTANALYSIS || defined ENABLE_ENTROPYANALYSIS
			csize+=fwrite(&bestrct, 1, 1, fdst);
#endif
			csize+=fwrite(cdatasize+0, 1, 4, fdst);
			csize+=fwrite(cdatasize+1, 1, 4, fdst);
			csize+=fwrite(cdata[0], 1, cdatasize[0], fdst);
			csize+=fwrite(cdata[1], 1, cdatasize[1], fdst);
			csize+=fwrite(cdata[2], 1, cdatasize[2], fdst);
			
#if defined _MSC_VER && defined LOUD
			csize_actual=csize;
#else
			(void)csize;
#endif
			fclose(fdst);
		}
		free(cdata[0]);
		free(cdata[1]);
		free(cdata[2]);
		//free(cdata);

		free(dstbuf);
	}
	else//decode
	{
		unsigned char yuv[4]={0};
#if defined ENABLE_FASTANALYSIS || defined ENABLE_ENTROPYANALYSIS
		int flag=*srcptr++;
		const unsigned char *rct=rct_indices[flag];
		int yidx=rct[3+0];
		int uidx=rct[3+1];
		int vidx=rct[3+2];
		int uhelpidx=rct[6+0];
		int vhelpidx=rct[6+1];
		//char *puhelper=yuv+uhelpidx;
		//char *pvhelper=yuv+vhelpidx;
#endif
		int imsize=3*iw*ih;
		unsigned char *image=(unsigned char*)malloc(imsize);
		if(!image)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		unsigned char *im2=(unsigned char*)malloc(imsize);
		if(!im2)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		int res=iw*ih;
		unsigned char *planes[]=
		{
			im2+res*0,
			im2+res*1,
			im2+res*2,
		};
		unsigned char *ret[3]={0};
		unsigned cdatasize[2]={0};
		memcpy(cdatasize+0, srcptr, 4); srcptr+=4;
		memcpy(cdatasize+1, srcptr, 4); srcptr+=4;
		unsigned char *cdata[3]={0};
		cdata[0]=srcptr;
		cdata[1]=cdata[0]+cdatasize[0];
		cdata[2]=cdata[1]+cdatasize[1];
#if defined _MSC_VER && defined LOUD
		etime=time_sec();
#endif
#ifdef ENABLE_MT
		if(maxthreads>1)
		{
			ThreadArgs args[]=
			{
				{cdata[0], planes[0], 0, cdatasize[0],			res},
				{cdata[1], planes[1], 0, cdatasize[1],			res},
				{cdata[2], planes[2], 0, (unsigned)(srcend-cdata[2]),	res},
			};
			void *mt=mt_exec(c20_dec, args, sizeof(*args), _countof(args));
			mt_finish(mt);
			ret[0]=args[0].ret;
			ret[1]=args[1].ret;
			ret[2]=args[2].ret;
		}
#endif
		{
			ret[0]=ENTROPY_DEC(cdata[0], cdatasize[0], planes[0], res);
			ret[1]=ENTROPY_DEC(cdata[1], cdatasize[1], planes[1], res);
			ret[2]=ENTROPY_DEC(cdata[2], (unsigned)(srcend-cdata[2]), planes[2], res);
		}
		if(!ret[0]||!ret[1]||!ret[2])
		{
			LOG_ERROR("HTS Decode Error");
			return 1;
		}
#if defined _MSC_VER && defined LOUD
		etime=time_sec()-etime;
		ptime=time_sec();
#endif
		//interleave |yyyy yyyy yyyy yyyy|uuuu uuuu uuuu uuuu|vvvv vvvv vvvv vvvv| -> |rgbr gbrg brgb rgbr|gbrg brgb rgbr gbrg|brgb rgbr gbrg brgb| (in order):
		__m256i setr0=_mm256_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
		//	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 15, 12,  9,  6,  3,  0
			 5, -1, -1,  4, -1, -1,  3, -1, -1,  2, -1, -1,  1, -1, -1,  0,
			 5, -1, -1,  4, -1, -1,  3, -1, -1,  2, -1, -1,  1, -1, -1,  0
		);
		__m256i setr1=_mm256_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
		//	-1, -1, -1, -1, -1, 14, 11,  8,  5,  2, -1, -1, -1, -1, -1, -1
			-1, 10, -1, -1,  9, -1, -1,  8, -1, -1,  7, -1, -1,  6, -1, -1,
			-1, 10, -1, -1,  9, -1, -1,  8, -1, -1,  7, -1, -1,  6, -1, -1
		);
		__m256i setr2=_mm256_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
		//	13, 10,  7,  4,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
			-1, -1, 15, -1, -1, 14, -1, -1, 13, -1, -1, 12, -1, -1, 11, -1,
			-1, -1, 15, -1, -1, 14, -1, -1, 13, -1, -1, 12, -1, -1, 11, -1
		);
		__m256i setg0=_mm256_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
		//	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 13, 10,  7,  4,  1
			-1, -1,  4, -1, -1,  3, -1, -1,  2, -1, -1,  1, -1, -1,  0, -1,
			-1, -1,  4, -1, -1,  3, -1, -1,  2, -1, -1,  1, -1, -1,  0, -1
		);
		__m256i setg1=_mm256_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
		//	-1, -1, -1, -1, -1, 15, 12,  9,  6,  3,  0, -1, -1, -1, -1, -1
			10, -1, -1,  9, -1, -1,  8, -1, -1,  7, -1, -1,  6, -1, -1,  5,
			10, -1, -1,  9, -1, -1,  8, -1, -1,  7, -1, -1,  6, -1, -1,  5
		);
		__m256i setg2=_mm256_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
		//	14, 11,  8,  5,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
			-1, 15, -1, -1, 14, -1, -1, 13, -1, -1, 12, -1, -1, 11, -1, -1,
			-1, 15, -1, -1, 14, -1, -1, 13, -1, -1, 12, -1, -1, 11, -1, -1
		);
		__m256i setb0=_mm256_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
		//	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 14, 11,  8,  5,  2
			-1,  4, -1, -1,  3, -1, -1,  2, -1, -1,  1, -1, -1,  0, -1, -1,
			-1,  4, -1, -1,  3, -1, -1,  2, -1, -1,  1, -1, -1,  0, -1, -1
		);
		__m256i setb1=_mm256_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
		//	-1, -1, -1, -1, -1, -1, 13, 10,  7,  4,  1, -1, -1, -1, -1, -1
			-1, -1,  9, -1, -1,  8, -1, -1,  7, -1, -1,  6, -1, -1,  5, -1,
			-1, -1,  9, -1, -1,  8, -1, -1,  7, -1, -1,  6, -1, -1,  5, -1
		);
		__m256i setb2=_mm256_set_epi8(
		//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
		//	15, 12,  9,  6,  3,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
			15, -1, -1, 14, -1, -1, 13, -1, -1, 12, -1, -1, 11, -1, -1, 10,
			15, -1, -1, 14, -1, -1, 13, -1, -1, 12, -1, -1, 11, -1, -1, 10
		);
		__m256i
			getr0, getr1, getr2,
			getg0, getg1, getg2,
			getb0, getb1, getb2;
		switch(flag)
		{
		default:
		case  0:getr0=setr0, getr1=setr1, getr2=setr2,  getg0=setg0, getg1=setg1, getg2=setg2,  getb0=setb0, getb1=setb1, getb2=setb2;break;//0, 1, 2
		case  1:getr0=setr0, getr1=setr1, getr2=setr2,  getg0=setg0, getg1=setg1, getg2=setg2,  getb0=setb0, getb1=setb1, getb2=setb2;break;//0, 1, 2
		case  2:getr0=setr0, getr1=setr1, getr2=setr2,  getg0=setg0, getg1=setg1, getg2=setg2,  getb0=setb0, getb1=setb1, getb2=setb2;break;//0, 1, 2
		case  3:getr0=setr0, getr1=setr1, getr2=setr2,  getg0=setg0, getg1=setg1, getg2=setg2,  getb0=setb0, getb1=setb1, getb2=setb2;break;//0, 1, 2
		case  4:getr0=setr0, getr1=setr1, getr2=setr2,  getg0=setg0, getg1=setg1, getg2=setg2,  getb0=setb0, getb1=setb1, getb2=setb2;break;//0, 1, 2
		case  5:getr0=setr0, getr1=setr1, getr2=setr2,  getg0=setb0, getg1=setb1, getg2=setb2,  getb0=setg0, getb1=setg1, getb2=setg2;break;//0, 2, 1
		case  6:getr0=setg0, getr1=setg1, getr2=setg2,  getg0=setb0, getg1=setb1, getg2=setb2,  getb0=setr0, getb1=setr1, getb2=setr2;break;//1, 2, 0
		case  7:getr0=setg0, getr1=setg1, getr2=setg2,  getg0=setb0, getg1=setb1, getg2=setb2,  getb0=setr0, getb1=setr1, getb2=setr2;break;//1, 2, 0
		case  8:getr0=setg0, getr1=setg1, getr2=setg2,  getg0=setb0, getg1=setb1, getg2=setb2,  getb0=setr0, getb1=setr1, getb2=setr2;break;//1, 2, 0
		case  9:getr0=setg0, getr1=setg1, getr2=setg2,  getg0=setb0, getg1=setb1, getg2=setb2,  getb0=setr0, getb1=setr1, getb2=setr2;break;//1, 2, 0
		case 10:getr0=setg0, getr1=setg1, getr2=setg2,  getg0=setr0, getg1=setr1, getg2=setr2,  getb0=setb0, getb1=setb1, getb2=setb2;break;//1, 0, 2
		case 11:getr0=setb0, getr1=setb1, getr2=setb2,  getg0=setr0, getg1=setr1, getg2=setr2,  getb0=setg0, getb1=setg1, getb2=setg2;break;//2, 0, 1
		case 12:getr0=setb0, getr1=setb1, getr2=setb2,  getg0=setr0, getg1=setr1, getg2=setr2,  getb0=setg0, getb1=setg1, getb2=setg2;break;//2, 0, 1
		case 13:getr0=setb0, getr1=setb1, getr2=setb2,  getg0=setr0, getg1=setr1, getg2=setr2,  getb0=setg0, getb1=setg1, getb2=setg2;break;//2, 0, 1
		case 14:getr0=setb0, getr1=setb1, getr2=setb2,  getg0=setr0, getg1=setr1, getg2=setr2,  getb0=setg0, getb1=setg1, getb2=setg2;break;//2, 0, 1
		case 15:getr0=setb0, getr1=setb1, getr2=setb2,  getg0=setg0, getg1=setg1, getg2=setg2,  getb0=setr0, getb1=setr1, getb2=setr2;break;//2, 1, 0
		}
		__m256i helperUmask	=_mm256_set1_epi8(-(uhelpidx!=3));
		__m256i helperVsel	=_mm256_set1_epi8(-(vhelpidx==1));
		__m256i helperVmask	=_mm256_set1_epi8(-(vhelpidx!=3));
		char helperUoffsetscalar=(char)(0x80&-(uhelpidx!=3));
		char helperVoffsetscalar=(char)(0x80&-(vhelpidx!=3));
		__m256i helperUoffset=_mm256_set1_epi8(helperUoffsetscalar);
		__m256i helperVoffset=_mm256_set1_epi8(helperVoffsetscalar);
		__m256i one=_mm256_set1_epi16(1);
		__m256i two=_mm256_set1_epi16(2);
		__m256i eight=_mm256_set1_epi16(8);
		//__m128i mhalf=_mm_set1_epi8(-128);
		__m256i mmin=_mm256_set1_epi16(-128);
		__m256i mmax=_mm256_set1_epi16(127);
		(void)one;
		(void)two;
		(void)eight;
		(void)mmin;
		(void)mmax;
		for(int ky=0, idx=0, idx2=0;ky<ih;++ky)
		{
			ALIGN(32) unsigned char *rows[]=
			{
				(unsigned char*)(pixels+((iw+32LL)*(((ky-0LL)&3)+0*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-1LL)&3)+0*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-2LL)&3)+0*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-3LL)&3)+0*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-0LL)&3)+1*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-1LL)&3)+1*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-2LL)&3)+1*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-3LL)&3)+1*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-0LL)&3)+2*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-1LL)&3)+2*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-2LL)&3)+2*4)+16LL)),
				(unsigned char*)(pixels+((iw+32LL)*(((ky-3LL)&3)+2*4)+16LL)),
			};
			int kx=0;
#if 1
			//for(;kx<=iw-16;kx+=16, idx+=3*16, idx2+=16)
			for(;kx<=iw-32;kx+=32, idx+=3*32, idx2+=32)
			{
				//process 32 pixels (96 SPs) at a time
#if defined USE_ROWPRED3A3
				__m256i NW0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]-1));
				__m256i NW1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]-1));
				__m256i NW2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]-1));
				__m256i N0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+0));
				__m256i N1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+0));
				__m256i N2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+0));
				__m256i NE0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+1));
				__m256i NE1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+1));
				__m256i NE2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+1));
				__m256i mp0=_mm256_add_epi16(NW0, NE0);//(10*N+3*(NW+NE)+8)>>4
				__m256i mp1=_mm256_add_epi16(NW1, NE1);
				__m256i mp2=_mm256_add_epi16(NW2, NE2);
				__m256i t0=_mm256_add_epi16(mp0, N0);
				__m256i t1=_mm256_add_epi16(mp1, N1);
				__m256i t2=_mm256_add_epi16(mp2, N2);
				t0	=_mm256_slli_epi16(t0, 1);
				t1	=_mm256_slli_epi16(t1, 1);
				t2	=_mm256_slli_epi16(t2, 1);
				mp0	=_mm256_add_epi16(mp0, t0);
				mp1	=_mm256_add_epi16(mp1, t1);
				mp2	=_mm256_add_epi16(mp2, t2);
				t0	=_mm256_slli_epi16(N0, 3);
				t1	=_mm256_slli_epi16(N1, 3);
				t2	=_mm256_slli_epi16(N2, 3);
				mp0	=_mm256_add_epi16(mp0, t0);
				mp1	=_mm256_add_epi16(mp1, t1);
				mp2	=_mm256_add_epi16(mp2, t2);
				mp0	=_mm256_add_epi16(mp0, eight);
				mp1	=_mm256_add_epi16(mp1, eight);
				mp2	=_mm256_add_epi16(mp2, eight);
				mp0	=_mm256_srai_epi16(mp0, 4);
				mp1	=_mm256_srai_epi16(mp1, 4);
				mp2	=_mm256_srai_epi16(mp2, 4);
				__m128i ypred=_mm_packs_epi16(_mm256_castsi256_si128(mp0), _mm256_extracti128_si256(mp0, 1));
				__m128i upred=_mm_packs_epi16(_mm256_castsi256_si128(mp1), _mm256_extracti128_si256(mp1, 1));
				__m128i vpred=_mm_packs_epi16(_mm256_castsi256_si128(mp2), _mm256_extracti128_si256(mp2, 1));
#elif defined USE_ROWPRED484
				__m256i NW0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]-1));
				__m256i NW1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]-1));
				__m256i NW2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]-1));
				__m256i N0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+0));
				__m256i N1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+0));
				__m256i N2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+0));
				__m256i NE0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+1));
				__m256i NE1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+1));
				__m256i NE2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+1));
				NW0=_mm256_avg_epu8(NW0, NE0);
				NW1=_mm256_avg_epu8(NW1, NE1);
				NW2=_mm256_avg_epu8(NW2, NE2);
				N0=_mm256_avg_epu8(N0, NW0);
				N1=_mm256_avg_epu8(N1, NW1);
				N2=_mm256_avg_epu8(N2, NW2);
				//__m256i mp0=_mm256_slli_epi16(N0, 1);//(2*N+NW+NE+2)>>2
				//__m256i mp1=_mm256_slli_epi16(N1, 1);
				//__m256i mp2=_mm256_slli_epi16(N2, 1);
				//mp0=_mm256_add_epi16(mp0, NW0);
				//mp1=_mm256_add_epi16(mp1, NW1);
				//mp2=_mm256_add_epi16(mp2, NW2);
				//mp0=_mm256_add_epi16(mp0, NE0);
				//mp1=_mm256_add_epi16(mp1, NE1);
				//mp2=_mm256_add_epi16(mp2, NE2);
				//mp0=_mm256_add_epi16(mp0, two);
				//mp1=_mm256_add_epi16(mp1, two);
				//mp2=_mm256_add_epi16(mp2, two);
				//mp0=_mm256_srai_epi16(mp0, 2);
				//mp1=_mm256_srai_epi16(mp1, 2);
				//mp2=_mm256_srai_epi16(mp2, 2);
				//__m128i ypred=_mm_packs_epi16(_mm256_castsi256_si128(mp0), _mm256_extracti128_si256(mp0, 1));
				//__m128i upred=_mm_packs_epi16(_mm256_castsi256_si128(mp1), _mm256_extracti128_si256(mp1, 1));
				//__m128i vpred=_mm_packs_epi16(_mm256_castsi256_si128(mp2), _mm256_extracti128_si256(mp2, 1));
#elif defined USE_ROWPRED565
				__m256i NW0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]-1));
				__m256i NW1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]-1));
				__m256i NW2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]-1));
				__m256i N0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+0));
				__m256i N1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+0));
				__m256i N2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+0));
				__m256i NE0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+1));
				__m256i NE1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+1));
				__m256i NE2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+1));
				__m256i mp0=_mm256_add_epi16(NW0, NE0);//(6*N+5*(NW+NE)+8)>>4
				__m256i mp1=_mm256_add_epi16(NW1, NE1);
				__m256i mp2=_mm256_add_epi16(NW2, NE2);
				__m256i t0=_mm256_add_epi16(mp0, N0);
				__m256i t1=_mm256_add_epi16(mp1, N1);
				__m256i t2=_mm256_add_epi16(mp2, N2);
				t0	=_mm256_slli_epi16(t0, 2);
				t1	=_mm256_slli_epi16(t1, 2);
				t2	=_mm256_slli_epi16(t2, 2);
				mp0	=_mm256_add_epi16(mp0, t0);
				mp1	=_mm256_add_epi16(mp1, t1);
				mp2	=_mm256_add_epi16(mp2, t2);
				t0	=_mm256_slli_epi16(N0, 1);
				t1	=_mm256_slli_epi16(N1, 1);
				t2	=_mm256_slli_epi16(N2, 1);
				mp0	=_mm256_add_epi16(mp0, t0);
				mp1	=_mm256_add_epi16(mp1, t1);
				mp2	=_mm256_add_epi16(mp2, t2);
				mp0	=_mm256_add_epi16(mp0, eight);
				mp1	=_mm256_add_epi16(mp1, eight);
				mp2	=_mm256_add_epi16(mp2, eight);
				mp0	=_mm256_srai_epi16(mp0, 4);
				mp1	=_mm256_srai_epi16(mp1, 4);
				mp2	=_mm256_srai_epi16(mp2, 4);
				__m128i ypred=_mm_packs_epi16(_mm256_castsi256_si128(mp0), _mm256_extracti128_si256(mp0, 1));
				__m128i upred=_mm_packs_epi16(_mm256_castsi256_si128(mp1), _mm256_extracti128_si256(mp1, 1));
				__m128i vpred=_mm_packs_epi16(_mm256_castsi256_si128(mp2), _mm256_extracti128_si256(mp2, 1));
#elif defined USE_PRED31
				__m256i NN0	=_mm256_loadu_si256((__m256i*)(rows[2+0*4]+0));
				__m256i NN1	=_mm256_loadu_si256((__m256i*)(rows[2+1*4]+0));
				__m256i NN2	=_mm256_loadu_si256((__m256i*)(rows[2+2*4]+0));
				__m256i N0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+0));
				__m256i N1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+0));
				__m256i N2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+0));
				__m256i mp0=_mm256_slli_epi16(N0, 1);//clamp((3*N-NN+1)>>1)
				__m256i mp1=_mm256_slli_epi16(N1, 1);
				__m256i mp2=_mm256_slli_epi16(N2, 1);
				mp0=_mm256_add_epi16(mp0, N0);
				mp1=_mm256_add_epi16(mp1, N1);
				mp2=_mm256_add_epi16(mp2, N2);
				mp0=_mm256_sub_epi16(mp0, NN0);
				mp1=_mm256_sub_epi16(mp1, NN1);
				mp2=_mm256_sub_epi16(mp2, NN2);
				mp0=_mm256_add_epi16(mp0, one);
				mp1=_mm256_add_epi16(mp1, one);
				mp2=_mm256_add_epi16(mp2, one);
				mp0=_mm256_srai_epi16(mp0, 1);
				mp1=_mm256_srai_epi16(mp1, 1);
				mp2=_mm256_srai_epi16(mp2, 1);
				mp0=_mm256_max_epi16(mp0, mmin);
				mp1=_mm256_max_epi16(mp1, mmin);
				mp2=_mm256_max_epi16(mp2, mmin);
				mp0=_mm256_min_epi16(mp0, mmax);
				mp1=_mm256_min_epi16(mp1, mmax);
				mp2=_mm256_min_epi16(mp2, mmax);
				__m128i ypred=_mm_packs_epi16(_mm256_castsi256_si128(mp0), _mm256_extracti128_si256(mp0, 1));
				__m128i upred=_mm_packs_epi16(_mm256_castsi256_si128(mp1), _mm256_extracti128_si256(mp1, 1));
				__m128i vpred=_mm_packs_epi16(_mm256_castsi256_si128(mp2), _mm256_extracti128_si256(mp2, 1));
#elif defined USE_PRED21
				__m256i NN0	=_mm256_loadu_si256((__m256i*)(rows[2+0*4]+0));
				__m256i NN1	=_mm256_loadu_si256((__m256i*)(rows[2+1*4]+0));
				__m256i NN2	=_mm256_loadu_si256((__m256i*)(rows[2+2*4]+0));
				__m256i N0	=_mm256_loadu_si256((__m256i*)(rows[1+0*4]+0));
				__m256i N1	=_mm256_loadu_si256((__m256i*)(rows[1+1*4]+0));
				__m256i N2	=_mm256_loadu_si256((__m256i*)(rows[1+2*4]+0));
				__m256i mp0=_mm256_slli_epi16(N0, 1);//clamp(2*N-NN)
				__m256i mp1=_mm256_slli_epi16(N1, 1);
				__m256i mp2=_mm256_slli_epi16(N2, 1);
				mp0=_mm256_sub_epi16(mp0, NN0);
				mp1=_mm256_sub_epi16(mp1, NN1);
				mp2=_mm256_sub_epi16(mp2, NN2);
				mp0=_mm256_max_epi16(mp0, mmin);
				mp1=_mm256_max_epi16(mp1, mmin);
				mp2=_mm256_max_epi16(mp2, mmin);
				mp0=_mm256_min_epi16(mp0, mmax);
				mp1=_mm256_min_epi16(mp1, mmax);
				mp2=_mm256_min_epi16(mp2, mmax);
				__m128i ypred=_mm_packs_epi16(_mm256_castsi256_si128(mp0), _mm256_extracti128_si256(mp0, 1));
				__m128i upred=_mm_packs_epi16(_mm256_castsi256_si128(mp1), _mm256_extracti128_si256(mp1, 1));
				__m128i vpred=_mm_packs_epi16(_mm256_castsi256_si128(mp2), _mm256_extracti128_si256(mp2, 1));
#elif defined USE_ROWPREDN
				__m128i ypred=_mm_loadu_si128((__m128i*)(rows[1+0*4]+0));
				__m128i upred=_mm_loadu_si128((__m128i*)(rows[1+1*4]+0));
				__m128i vpred=_mm_loadu_si128((__m128i*)(rows[1+2*4]+0));
#endif

				//if(ky==1&&kx==192)//
				//if(ky==2&&kx==0)//
				//	printf("");

				__m256i p0=_mm256_loadu_si256((__m256i*)(planes[0]+idx2));
				__m256i p1=_mm256_loadu_si256((__m256i*)(planes[1]+idx2));
				__m256i p2=_mm256_loadu_si256((__m256i*)(planes[2]+idx2));
				p0=_mm256_add_epi8(p0, N0);
				p1=_mm256_add_epi8(p1, N1);
				p2=_mm256_add_epi8(p2, N2);
				_mm256_storeu_si256((__m256i*)rows[0+0*4], p0);
				_mm256_storeu_si256((__m256i*)rows[0+1*4], p1);
				_mm256_storeu_si256((__m256i*)rows[0+2*4], p2);
//#ifdef USE_ROWPREDN
//				_mm_storeu_si128((__m128i*)rows[0+0*4], p0);
//				_mm_storeu_si128((__m128i*)rows[0+1*4], p1);
//				_mm_storeu_si128((__m128i*)rows[0+2*4], p2);
//#else
//				_mm256_storeu_si256((__m256i*)rows[0+0*4], _mm256_cvtepi8_epi16(p0));
//				_mm256_storeu_si256((__m256i*)rows[0+1*4], _mm256_cvtepi8_epi16(p1));
//				_mm256_storeu_si256((__m256i*)rows[0+2*4], _mm256_cvtepi8_epi16(p2));
//#endif

				__m256i helperU0=_mm256_and_si256(p0, helperUmask);
				p1=_mm256_add_epi8(p1, helperU0);
				p1=_mm256_sub_epi8(p1, helperUoffset);
				__m256i helperV0=_mm256_blendv_epi8(p0, p1, helperVsel);
				helperV0=_mm256_and_si256(helperV0, helperVmask);
				p2=_mm256_add_epi8(p2, helperV0);
				p2=_mm256_sub_epi8(p2, helperVoffset);

				//interleave |yyyy yyyy yyyy yyyy|uuuu uuuu uuuu uuuu|vvvv vvvv vvvv vvvv| -> |rgbr gbrg brgb rgbr|gbrg brgb rgbr gbrg|brgb rgbr gbrg brgb| (in order):
				__m256i r0=_mm256_shuffle_epi8(p0, getr0);
				__m256i r1=_mm256_shuffle_epi8(p0, getr1);
				__m256i r2=_mm256_shuffle_epi8(p0, getr2);
				__m256i g0=_mm256_shuffle_epi8(p1, getg0);
				__m256i g1=_mm256_shuffle_epi8(p1, getg1);
				__m256i g2=_mm256_shuffle_epi8(p1, getg2);
				__m256i b0=_mm256_shuffle_epi8(p2, getb0);
				__m256i b1=_mm256_shuffle_epi8(p2, getb1);
				__m256i b2=_mm256_shuffle_epi8(p2, getb2);
				r0=_mm256_or_si256(r0, g0);
				r1=_mm256_or_si256(r1, g1);
				r2=_mm256_or_si256(r2, g2);
				r0=_mm256_or_si256(r0, b0);
				r1=_mm256_or_si256(r1, b1);
				r2=_mm256_or_si256(r2, b2);

				//r0=_mm256_sub_epi8(r0, mhalf);
				//r1=_mm256_sub_epi8(r1, mhalf);
				//r2=_mm256_sub_epi8(r2, mhalf);
				_mm_storeu_si128((__m128i*)(image+idx+0*16), _mm256_castsi256_si128(r0));
				_mm_storeu_si128((__m128i*)(image+idx+1*16), _mm256_castsi256_si128(r1));
				_mm_storeu_si128((__m128i*)(image+idx+2*16), _mm256_castsi256_si128(r2));
				_mm_storeu_si128((__m128i*)(image+idx+3*16), _mm256_extracti128_si256(r0, 1));
				_mm_storeu_si128((__m128i*)(image+idx+4*16), _mm256_extracti128_si256(r1, 1));
				_mm_storeu_si128((__m128i*)(image+idx+5*16), _mm256_extracti128_si256(r2, 1));
#ifdef ENABLE_GUIDE
				for(int k=0;k<16;++k)
					guide_check(image, kx+k, ky);
#endif
				
				__m256i rows0=_mm256_load_si256((__m256i*)rows+0);
				__m256i rows1=_mm256_load_si256((__m256i*)rows+1);
				__m256i rows2=_mm256_load_si256((__m256i*)rows+2);
				rows0=_mm256_add_epi64(rows0, _mm256_set1_epi64x(sizeof(char[32])));
				rows1=_mm256_add_epi64(rows1, _mm256_set1_epi64x(sizeof(char[32])));
				rows2=_mm256_add_epi64(rows2, _mm256_set1_epi64x(sizeof(char[32])));
				_mm256_store_si256((__m256i*)rows+0, rows0);
				_mm256_store_si256((__m256i*)rows+1, rows1);
				_mm256_store_si256((__m256i*)rows+2, rows2);
			}
#endif
//#ifdef __GNUC__
//#pragma GCC unroll 2
//#endif
			for(;kx<iw;++kx, idx+=3, ++idx2)
			{
				//if(ky==0&&kx==6&&kc==0)//
				//if(ky==399&&kx==714)//
				//if(ky==2&&kx==4)//
				//if(ky==355&&kx==647)//
				//if(ky==3&&kx==104)//
				//if(ky==3&&kx==103)//
				//if(ky==4&&kx==138)//
				//	printf("");

				//short *pcurr=rows[0]+0*4;
				
#if defined USE_ROWPRED3A3
				int ypred=(10*rows[1+0*4][+0]+3*(rows[1+0*4][-1]+rows[1+0*4][+1])+8)>>4;
				int upred=(10*rows[1+1*4][+0]+3*(rows[1+1*4][-1]+rows[1+1*4][+1])+8)>>4;
				int vpred=(10*rows[1+2*4][+0]+3*(rows[1+2*4][-1]+rows[1+2*4][+1])+8)>>4;
#elif defined USE_ROWPRED484
				int ypred=(rows[1+0*4][+0]+((rows[1+0*4][-1]+rows[1+0*4][+1]+1)>>1)+1)>>1;
				int upred=(rows[1+1*4][+0]+((rows[1+1*4][-1]+rows[1+1*4][+1]+1)>>1)+1)>>1;
				int vpred=(rows[1+2*4][+0]+((rows[1+2*4][-1]+rows[1+2*4][+1]+1)>>1)+1)>>1;
				//int ypred=(2*rows[1+0*4][+0]+rows[1+0*4][-1]+rows[1+0*4][+1]+2)>>2;
				//int upred=(2*rows[1+1*4][+0]+rows[1+1*4][-1]+rows[1+1*4][+1]+2)>>2;
				//int vpred=(2*rows[1+2*4][+0]+rows[1+2*4][-1]+rows[1+2*4][+1]+2)>>2;
#elif defined USE_ROWPRED565
				int ypred=(6*rows[1+0*4][+0]+5*(rows[1+0*4][-1]+rows[1+0*4][+1])+8)>>4;
				int upred=(6*rows[1+1*4][+0]+5*(rows[1+1*4][-1]+rows[1+1*4][+1])+8)>>4;
				int vpred=(6*rows[1+2*4][+0]+5*(rows[1+2*4][-1]+rows[1+2*4][+1])+8)>>4;
#elif defined USE_PRED31
				int ypred=(3*rows[1+0*4][+0]-rows[2+0*4][+0]+1)>>1;
				int upred=(3*rows[1+1*4][+0]-rows[2+1*4][+0]+1)>>1;
				int vpred=(3*rows[1+2*4][+0]-rows[2+2*4][+0]+1)>>1;
				CLAMP2(ypred, -128, 127);
				CLAMP2(upred, -128, 127);
				CLAMP2(vpred, -128, 127);
#elif defined USE_PRED21
				int ypred=2*rows[1+0*4][+0]-rows[2+0*4][+0];
				int upred=2*rows[1+1*4][+0]-rows[2+1*4][+0];
				int vpred=2*rows[1+2*4][+0]-rows[2+2*4][+0];
				CLAMP2(ypred, -128, 127);
				CLAMP2(upred, -128, 127);
				CLAMP2(vpred, -128, 127);
#elif defined USE_ROWPREDN
				int ypred=rows[1+0*4][+0];
				int upred=rows[1+1*4][+0];
				int vpred=rows[1+2*4][+0];
#endif
				
				rows[0+0*4][+0]=yuv[0]=(unsigned char)(planes[0][idx2]+ypred);
				rows[0+1*4][+0]=yuv[1]=(unsigned char)(planes[1][idx2]+upred);
				rows[0+2*4][+0]=yuv[2]=(unsigned char)(planes[2][idx2]+vpred);
				yuv[1]+=yuv[uhelpidx]-helperUoffsetscalar;
				yuv[2]+=yuv[vhelpidx]-helperVoffsetscalar;
				
				image[idx+yidx]=yuv[0];
				image[idx+uidx]=yuv[1];
				image[idx+vidx]=yuv[2];
				//image[idx+yidx]=yuv[0]+128;
				//image[idx+uidx]=yuv[1]+128;
				//image[idx+vidx]=yuv[2]+128;

				guide_check(image, kx, ky);
				__m256i rows0=_mm256_load_si256((__m256i*)rows+0);
				__m256i rows1=_mm256_load_si256((__m256i*)rows+1);
				__m256i rows2=_mm256_load_si256((__m256i*)rows+2);
				rows0=_mm256_add_epi64(rows0, _mm256_set1_epi64x(sizeof(char)));
				rows1=_mm256_add_epi64(rows1, _mm256_set1_epi64x(sizeof(char)));
				rows2=_mm256_add_epi64(rows2, _mm256_set1_epi64x(sizeof(char)));
				_mm256_store_si256((__m256i*)rows+0, rows0);
				_mm256_store_si256((__m256i*)rows+1, rows1);
				_mm256_store_si256((__m256i*)rows+2, rows2);
			}
		}
#if defined _MSC_VER && defined LOUD
		ptime=time_sec()-ptime;
#endif
		free(im2);
		{
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
				return 1;
			}
			fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(image, 1, imsize, fdst);
			fclose(fdst);
		}
		free(image);
	}
	_mm_free(pixels);
	free(srcbuf);
#if defined _MSC_VER && defined LOUD
	cycles=__rdtsc()-cycles;
	elapsed=time_sec()-elapsed;
	ptrdiff_t usize=3LL*iw*ih;
	if(fwd)
		printf(" %9td   /%9td %12.6lf:1\n", csize_actual, usize, (double)usize/csize_actual);
	printf("%c%12lf sec %12lf MB/s  %12lld cycles %12lf C/B  pred %lf sec  EC %lf sec\n",
		'D'+fwd,
		elapsed,
		usize/(elapsed*1024*1024),
		cycles,
		(double)cycles/usize,
		ptime,
		etime
	);
#endif
	return 0;
}