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
#include<emmintrin.h>


#ifdef _MSC_VER
	#define LOUD
	#define PRINT_PREDS
//	#define PRINT_RCT
	#define ENABLE_GUIDE
//	#define ENABLE_VAL
//	#define SAVE_RESIDUALS
#endif

	#define ENABLE_L1

	#define BLOCKX 1024
	#define BLOCKY 1024
//	#define BLOCKX 448
//	#define BLOCKY 448
//	#define BLOCKX 64
//	#define BLOCKY 64

//from libjxl
//sym = packsign(delta) = 0b00001MMBB...BBL
//token = offset + 0bGGGGMML,  where G = bits of FLOOR_LOG2(sym)-E,  bypass = 0bBB...BB  FLOOR_LOG2(sym)-(M+L) bits
#define CONFIG_EXP 4	//revise AVX2 CDF search to change config
#define CONFIG_MSB 1	//410->1+24+1	411->1+32+1	421->1+48+1
#define CONFIG_LSB 1	//		511->1+44+1
#define NLEVELS 33

#define NCTX 16
#define PBITS 2

#define NPREDS 13

#define L1PREDS 10
#define L1SH 19

#define PROBBITS_STORE	24
#define PROBBITS_USE	16

#define RATE 9


#ifdef _MSC_VER
#define AWM_INLINE __forceinline static
#else
#define AWM_INLINE __attribute__((always_inline)) inline static
#endif
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
AWM_INLINE int32_t floor_log2(uint32_t n)
{
#ifdef _MSC_VER
	unsigned long index=0;
	_BitScanReverse(&index, n);//3 cycles
	return n?index:-1;
#elif defined __GNUC__
	return n?31-__builtin_clz(n):-1;
#else
	if(!n)
		return -1;
	if(n>0x7FFFFFFF)
		return 31;
	{
		//9 cycles excluding CMOVs above
		__m128i t0=_mm_castpd_si128(_mm_cvtsi32_sd(_mm_setzero_pd(), n));
		t0=_mm_srli_epi64(t0, 52);
		return _mm_cvtsi128_si32(t0)-1023;
	}
#endif
}
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
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
		printf("\n");
	}
	printf("CRASH\n");
#ifdef _MSC_VER
	system("pause");
#endif
	exit(1);
}
#define CRASH(FORMAT, ...) crash(__FILE__, __LINE__, FORMAT, ##__VA_ARGS__)
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
#else
#define guide_save(...)
#define guide_check(...)
#endif
#ifdef ENABLE_VAL
typedef struct _ValInfo
{
	int cdf, freq;
	unsigned long long low0, range0, low1, range1;
} ValInfo;
static uint32_t valcap=1, valcount=0, validx=0;
static ValInfo *valdata=0;
static void val_enc(int cdf, int freq, unsigned long long low0, unsigned long long range0, unsigned long long low1, unsigned long long range1)
{
	{
		void *ptr;

		++valcount;
		if(valcount>=valcap)
		{
			valcap<<=1;
			ptr=realloc(valdata, sizeof(ValInfo)*valcap);
			if(!ptr)
			{
				CRASH("Alloc error");
				return;
			}
			valdata=(ValInfo*)ptr;
		}
	}
	{
		ValInfo info={cdf, freq, low0, range0, low1, range1};
		valdata[valcount-1]=info;
	}
}
static int val_dec(int cdf, int freq, unsigned long long low0, unsigned long long range0, unsigned long long code0, unsigned long long low1, unsigned long long range1, unsigned long long code1)
{
	ValInfo *info=valdata+validx;
	if(info->cdf!=cdf||info->freq!=freq||info->low0!=low0||info->range0!=range0||info->low1!=low1||info->range1!=range1)
	{
		printf("OG  %d %d\n", info->cdf, info->freq);
		printf("  0: %016llX %016llX\n", info->low0, info->range0);
		printf("  1: %016llX %016llX\n", info->low1, info->range1);
		printf("X   %d %d\n", cdf, freq);
		printf("  0: %016llX %016llX %016llX\n", low0, range0, code0);
		printf("  1: %016llX %016llX %016llX\n", low1, range1, code1);
		return 1;
	}
	++validx;
	return 0;
}
#define VAL_DEC(CDF, FREQ, LOW0, RANGE0, CODE0, LOW1, RANGE1, CODE1, ...) if(val_dec(CDF, FREQ, LOW0, RANGE0, CODE0, LOW1, RANGE1, CODE1))CRASH(__VA_ARGS__)
#else
#define val_enc(...)
#define val_dec(...)
#define VAL_DEC(...)
#endif
static unsigned char* load_file(const char *fn, ptrdiff_t *ret_size)
{
	struct stat info={0};
	int error;
	ptrdiff_t size, nread;
	FILE *fsrc;
	unsigned char *buf;

	error=stat(fn, &info);
	if(error)
	{
		printf("Cannot stat \"%s\"\n", fn);
		return 0;
	}
	size=info.st_size;
	fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		printf("Cannot open \"%s\"\n", fn);
		return 0;
	}
	buf=(unsigned char*)malloc(size+16);
	if(!buf)
	{
		CRASH("Alloc error\n");
		return 0;
	}
	nread=fread(buf, 1, size, fsrc);
	if(nread!=size)
		printf("Error reading \"%s\"\n", fn);
	fclose(fsrc);
	*ret_size=size;
	return buf;
}
#define OCHLIST\
	OCH(Y400) OCH(Y040) OCH(Y004)\
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
#if 0
typedef enum _OCHType
{
	OCH_R,
	OCH_G,
	OCH_B,
	OCH_RG,
	OCH_GB,
	OCH_BR,
//	OCH_R2,
//	OCH_G2,
//	OCH_B2,

	OCH_COUNT,

	OCH_RB=OCH_BR,
	OCH_GR=OCH_RG,
	OCH_BG=OCH_GB,

} OCHType;
static const unsigned char rct_indices[][8]=
{//	output channels			permutation	helper index
	{OCH_R,	OCH_G,	OCH_B,		0, 1, 2,	3, 3},// 0
	{OCH_R,	OCH_G,	OCH_BG,		0, 1, 2,	3, 1},// 1
	{OCH_R,	OCH_G,	OCH_BR,		0, 1, 2,	3, 0},// 2
	{OCH_R,	OCH_GR,	OCH_BR,		0, 1, 2,	0, 0},// 3
	{OCH_R,	OCH_GR,	OCH_BG,		0, 1, 2,	0, 1},// 4
	{OCH_R,	OCH_BR,	OCH_GB,		0, 2, 1,	0, 1},// 5
	{OCH_G,	OCH_B,	OCH_RG,		1, 2, 0,	3, 0},// 6
	{OCH_G,	OCH_B,	OCH_RB,		1, 2, 0,	3, 1},// 7
	{OCH_G,	OCH_BG,	OCH_RG,		1, 2, 0,	0, 0},// 8
	{OCH_G,	OCH_BG,	OCH_RB,		1, 2, 0,	0, 1},// 9
	{OCH_G,	OCH_RG,	OCH_BR,		1, 0, 2,	0, 1},//10
	{OCH_B,	OCH_R,	OCH_GR,		2, 0, 1,	3, 1},//11
	{OCH_B,	OCH_R,	OCH_GB,		2, 0, 1,	3, 0},//12
	{OCH_B,	OCH_RB,	OCH_GB,		2, 0, 1,	0, 0},//13
	{OCH_B,	OCH_RB,	OCH_GR,		2, 0, 1,	0, 1},//14
	{OCH_B,	OCH_GB,	OCH_RG,		2, 1, 0,	0, 1},//15
};
#endif
#define PREDLIST\
	PRED(W)\
	PRED(AV2)\
	PRED(SELECT)\
	PRED(NW)\
	PRED(CG)\
	PRED(AV3)\
	PRED(IZ)\
	PRED(AV4)\
	PRED(AV6)\
	PRED(AV8)\
	PRED(AV9)\
	PRED(WG)\
	PRED(L1)
typedef enum _PredType
{
#define PRED(LABEL) PRED_##LABEL,
	PREDLIST
#undef  PRED

	PRED_COUNT,
} PredType;
static const char *pred_names[]=
{
#define PRED(LABEL) #LABEL,
	PREDLIST
#undef  PRED
};
typedef struct _QuantInfoFwd
{
	uint8_t token, nbypass;
	uint16_t bypass;
} QuantInfoFwd;
typedef struct _QuantInfoInv
{
	uint8_t sym, nbits;
} QuantInfoInv;
static uint32_t stats[OCH_COUNT][1<<PBITS][NCTX][NLEVELS+1], mixin[NLEVELS*2];
int s07_codec(const char *command, const char *srcfn, const char *dstfn)
{
	//const int mem=sizeof(stats)+sizeof(mixer11)+sizeof(mixer12)+sizeof(mixer13)+sizeof(mixer14)+sizeof(mixer2);
	unsigned char *srcbuf=0, *srcptr=0;
	ptrdiff_t srcsize=0, dstsize=0;
	int fwd=0;
	int32_t iw=0, ih=0;
	ptrdiff_t res=0, usize=0, csize=0;
	unsigned char *image=0;
	unsigned char *dstbuf=0;
	unsigned char *streamptr=0, *streamend=0;
	unsigned char *imptr=0;
	uint8_t *decinfo=0;
	//int bestrct=0;
	//int NWratio=0;
	//uint8_t *predidx=0;
#ifdef LOUD
	double t=time_sec();
#endif
	uint64_t low=0, range=0xFFFFFFFFFFFF, code=0;
#ifdef SAVE_RESIDUALS
	uint8_t *residuals;
#endif

	int xblocks, yblocks, nblocks;

	(void)memfill;

	srcbuf=load_file(srcfn, &srcsize);
	if(!srcbuf)
	{
		CRASH("Cannot open \"%s\"", srcfn);
		return 1;
	}
	srcptr=srcbuf;
	{
		int tag=*(uint16_t*)srcptr;
		fwd=tag==('P'|'6'<<8);
		if(!fwd&&tag!=('0'|'7'<<8))
		{
			CRASH("Unsupported file  \"%s\"\n", srcfn);
			free(srcbuf);
			return 1;
		}
		srcptr+=2;
	}
	if(fwd)
	{
		if(*srcptr++!='\n')
		{
			CRASH("Unsupported file\n");
			free(srcbuf);
			return 1;
		}
		while(*srcptr=='#')//skip comment(s)
		{
			++srcptr;
			while(*srcptr++!='\n');
		}
		iw=0;
		while((uint32_t)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++-'0';
		while(*srcptr<=' ')++srcptr;//skip space
		ih=0;
		while((uint32_t)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++-'0';
		if(*srcptr++!='\n')
		{
			CRASH("Unsupported file\n");
			free(srcbuf);
			return 1;
		}
		while(*srcptr=='#')//skip comment(s)
		{
			++srcptr;
			while(*srcptr++!='\n');
		}
		if(memcmp(srcptr, "255\n", 4))
		{
			CRASH("Unsupported file\n");
			free(srcbuf);
			return 1;
		}
		srcptr+=4;
		image=srcptr;
		guide_save(image, iw, ih);
	}
	else
	{
		iw=((uint32_t*)srcptr)[0];
		ih=((uint32_t*)srcptr)[1];
		srcptr+=sizeof(uint32_t[2]);
	}
	if(iw<1||ih<1)
	{
		CRASH("Invalid file\n");
		free(srcbuf);
		return 1;
	}
	res=(ptrdiff_t)iw*ih;
	usize=3*res;
	dstbuf=(unsigned char*)malloc(usize);
	xblocks=(iw+BLOCKX-1)/BLOCKX;
	yblocks=(ih+BLOCKY-1)/BLOCKY;
	nblocks=xblocks*yblocks;
	decinfo=(uint8_t*)malloc(nblocks*sizeof(uint8_t[4]));
#ifdef SAVE_RESIDUALS
	residuals=0;
	if(fwd)
	{
		residuals=(uint8_t*)malloc(usize);
		if(!residuals)
		{
			CRASH("Alloc error");
			free(srcbuf);
			return 1;
		}
	}
#endif
	if(!dstbuf||!decinfo)
	{
		CRASH("Alloc error");
		free(srcbuf);
		return 1;
	}
	if(fwd)
	{
		//analysis
		int ky, kx, kb;
	//	int64_t counters[OCH_COUNT*NPREDS]={0}, minerr=0;
	//	int rowstride=3*iw;
		int32_t paddedwidth=iw+16;
		int32_t psize=(int32_t)sizeof(int16_t[4*OCH_COUNT])*paddedwidth;//4 padded rows * 6 channels
		int16_t *pixels=(int16_t*)malloc(psize);
		int32_t hsize=(int32_t)sizeof(int32_t[OCH_COUNT*NPREDS*512])*nblocks;
		int32_t *hist=(int32_t*)malloc(hsize);
		int32_t coeffs[OCH_COUNT][L1PREDS]={0};

		if(!pixels||!hist)
		{
			CRASH("Alloc error\n");
			free(srcbuf);
			free(dstbuf);
			return 1;
		}
		memset(pixels, 0, psize);
		memset(hist, 0, hsize);

		imptr=image;
		for(ky=0;ky<ih;++ky)
		{
			int by=ky/BLOCKY;
			int16_t *rows[]=
			{
				pixels+(paddedwidth*((ky-0)&3)+8)*OCH_COUNT,
				pixels+(paddedwidth*((ky-1)&3)+8)*OCH_COUNT,
				pixels+(paddedwidth*((ky-2)&3)+8)*OCH_COUNT,
				pixels+(paddedwidth*((ky-3)&3)+8)*OCH_COUNT,
			};
			for(kx=0;kx<iw;++kx, imptr+=3)
			{
				int bx=kx/BLOCKX;
				int32_t *currhist=hist+(xblocks*by+bx)*(OCH_COUNT*NPREDS*512);
				int16_t
				//	*NNN	=rows[3]+0*OCH_COUNT,
				//	*NN	=rows[2]+0*OCH_COUNT,
				//	*NW	=rows[1]-1*OCH_COUNT,
				//	*N	=rows[1]+0*OCH_COUNT,
				//	*NE	=rows[1]+1*OCH_COUNT,
				//	*WWW	=rows[0]-3*OCH_COUNT,
				//	*WW	=rows[0]-2*OCH_COUNT,
				//	*W	=rows[0]-1*OCH_COUNT,
					*curr	=rows[0]+0*OCH_COUNT;
				int32_t kc;
				
				curr[0x0]=imptr[0];		//r
				curr[0x1]=imptr[1];		//g
				curr[0x2]=imptr[2];		//b
				curr[0x3]=curr[0]-curr[1];	//rg
				curr[0x4]=curr[1]-curr[2];	//gb
				curr[0x5]=curr[2]-curr[0];	//br
				curr[0x6]=+curr[3]+curr[4]/4;	//r-(3*g+b)/4 = r-g-(b-g)/4
				curr[0x7]=-curr[3]-curr[5]/4;	//g-(3*r+b)/4 = g-r-(b-r)/4
				curr[0x8]=+curr[5]+curr[3]/4;	//b-(3*r+g)/4 = b-r-(g-r)/4
				curr[0x9]=-curr[5]-curr[4]/4;	//r-(g+3*b)/4 = r-b-(g-b)/4
				curr[0xA]=+curr[4]+curr[5]/4;	//g-(r+3*b)/4 = g-b-(r-b)/4
				curr[0xB]=-curr[4]-curr[3]/4;	//b-(r+3*g)/4 = b-g-(r-g)/4
				curr[0xC]=(+curr[3]-curr[5])/2;//r-(g+b)/2 = (r-g + r-b)/2
				curr[0xD]=(-curr[3]+curr[4])/2;//g-(r+b)/2 = (g-r + g-b)/2
				curr[0xE]=(+curr[5]-curr[4])/2;//b-(r+g)/2 = (b-r + b-g)/2

				for(kc=0;kc<OCH_COUNT;++kc)
				{
					int16_t
					//	NNNN	=rows[0][kc+0*OCH_COUNT],//X
						NNN	=rows[3][kc+0*OCH_COUNT],
						NNW	=rows[2][kc-1*OCH_COUNT],
						NN	=rows[2][kc+0*OCH_COUNT],
						NNE	=rows[2][kc+1*OCH_COUNT],
						NWW	=rows[1][kc-2*OCH_COUNT],
						NW	=rows[1][kc-1*OCH_COUNT],
						N	=rows[1][kc+0*OCH_COUNT],
						NE	=rows[1][kc+1*OCH_COUNT],
						NEE	=rows[1][kc+2*OCH_COUNT],
						NEEE	=rows[1][kc+3*OCH_COUNT],
						NEEEE	=rows[1][kc+4*OCH_COUNT],
						WWWW	=rows[0][kc-4*OCH_COUNT],
						WWW	=rows[0][kc-3*OCH_COUNT],
						WW	=rows[0][kc-2*OCH_COUNT],
						W	=rows[0][kc-1*OCH_COUNT],
						curr2	=rows[0][kc+0*OCH_COUNT];
					int32_t kp;
					int16_t preds[NPREDS], l1preds[L1PREDS];
					int32_t tmp, vmin, vmax;
					int16_t gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1;
					int16_t gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1;
					int j=0, i=0;
					
					(void)NNN	;
					(void)NNW	;
					(void)NN	;
					(void)NNE	;
					(void)NWW	;
					(void)NW	;
					(void)N		;
					(void)NE	;
					(void)NEE	;
					(void)NEEE	;
					(void)NEEEE	;
					(void)WWWW	;
					(void)WWW	;
					(void)WW	;
					(void)W		;
					(void)curr2	;

					preds[j++]=W;
					preds[j++]=(N+W)/2;
					preds[j++]=abs(N-NW)>abs(W-NW)?N:W;
					preds[j++]=NW;

					vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					tmp=N+W-NW;
					CLAMP2(tmp, vmin, vmax);
					preds[j++]=tmp;

					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					
					tmp=(5*(N+W)-2*NW)/8;
					CLAMP2(tmp, vmin, vmax);
					preds[j++]=tmp;

					tmp=(3*(N+W)-2*NW)/4;
					CLAMP2(tmp, vmin, vmax);
					preds[j++]=tmp;
					
					tmp=(4*(N+W)+NE-NW)/8;
					CLAMP2(tmp, vmin, vmax);
					preds[j++]=tmp;

					tmp=W+(6*N-5*NW-NN-WW+NE)/8;
					CLAMP2(tmp, vmin, vmax);
					preds[j++]=tmp;

					preds[j++]=(NN+NNE+NW+N+NE+NEE+WW+W)/8;

					tmp=W+(10*N+9*NW+4*NE-2*(NN+WW)+NNW-(NNE+NWW))/16;
					CLAMP2(tmp, vmin, vmax);
					preds[j++]=tmp;

					preds[j++]=(gx*N+gy*W)/(gx+gy);
					
					//if(ky==ih/2&&kx==iw/2)//
					//	printf("");

					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					l1preds[i++]=N;
					l1preds[i++]=W;
					l1preds[i++]=3*(N-NN)+NNN;
					l1preds[i++]=3*(W-WW)+WWW;
					l1preds[i++]=4*(W+WWW)-6*WW-WWWW;
				//	l1preds[i++]=4*(N+NNN)-6*NN-NNNN;
					l1preds[i++]=N+W-NW;
					l1preds[i++]=W+NE-N;
					l1preds[i++]=N+NE-NNE;
					l1preds[i++]=(WWWW+NNN+NEEE+NEEEE)/4;
					l1preds[i++]=NEEE;
					//tmp=(
					//	+coeffs[kc][0]*NW
					//	+coeffs[kc][1]*N
					//	+coeffs[kc][2]*NE
					//	+coeffs[kc][3]*W
					//	+(1<<16>>1)
					//)>>16;
					tmp=1<<L1SH>>1;
					for(i=0;i<L1PREDS;++i)
						tmp+=coeffs[kc][i]*l1preds[i];
					tmp>>=L1SH;
					CLAMP2(tmp, vmin, vmax);
					preds[j++]=tmp;//tmp = L1
#if 0
					int16_t preds[]=
					{
						W,
						(N+W)/2,
						abs(N-NW)>abs(W-NW)?N:W,
						NW,
						N+W-NW,
						(5*(N+W)-2*NW)/8,
						(3*(N+W)-2*NW)/4,
						(4*(N+W)+NE-NW)/8,
						W+(6*N-5*NW-NN-WW+NE)/8,
						(NN+NNE+NW+N+NE+NEE+WW+W)/8,
						W+(10*N+9*NW+4*NE-2*(NN+WW)+NNW-(NNE+NWW))/16,
						(gx*N+gy*W)/(gx+gy),

						//N,
						//W,
						//(N+W)/2,
						//3*(N-NN)+NNN,
						//3*(W-WW)+WWW,
						//N+W-NW,
						//W+NE-N,
						//W+NEE-NE,
						//N+NE-NNE,
						//(WWWW+NNN+NEEE+NEEEE)/4,
						//NWratio>128?W:(NWratio<128?N:(N+W)/2),
						//NEEE,
					};
#endif
					for(kp=0;kp<NPREDS;++kp)
						++currhist[512*(OCH_COUNT*kp+kc)+((curr2-preds[kp]+256)&511)];
					//	counters[kc+kp*OCH_COUNT]+=abs(curr[kc]-preds[kp]);

					{
						int e=(curr2>tmp)-(curr2<tmp);
						for(i=0;i<L1PREDS;++i)
							coeffs[kc][i]+=e*l1preds[i];
					}
					//{
					//	int e=preds[NPREDS-1];
					//	e=(curr2>e)-(curr2<e);
					//	coeffs[kc][0]+=e*NW;
					//	coeffs[kc][1]+=e*N;
					//	coeffs[kc][2]+=e*NE;
					//	coeffs[kc][3]+=e*W;
					//}
				}
				rows[0]+=OCH_COUNT;
				rows[1]+=OCH_COUNT;
				rows[2]+=OCH_COUNT;
				rows[3]+=OCH_COUNT;
			}
		}
		imptr=image;
		for(kb=0;kb<nblocks;++kb)
		{
			int kc, kp, ks, kt;
			double csizes[OCH_COUNT*NPREDS], bestsize;
			uint8_t predsel[OCH_COUNT];
			int bestrct=0;

			for(kc=0;kc<OCH_COUNT;++kc)
			{
				for(kp=0;kp<NPREDS;++kp)
				{
					double e=0, norm;
					int32_t *currhist=hist+((ptrdiff_t)OCH_COUNT*(NPREDS*kb+kp)+kc)*512, sum=0;
					for(ks=0;ks<512;++ks)
						sum+=currhist[ks];
					norm=1./sum;
					for(ks=0;ks<512;++ks)
					{
						int32_t freq=currhist[ks];
						if(freq)
							e-=freq*log(freq*norm);
					}
					csizes[OCH_COUNT*kp+kc]=e*(1./(M_LN2*8));
				}
			}
			for(kc=0;kc<OCH_COUNT;++kc)//select preds
			{
				double *currsizes=csizes+kc;
				int bestpred=0;
				for(kp=bestpred+1;kp<NPREDS;++kp)
				{
					if(currsizes[OCH_COUNT*bestpred]>currsizes[OCH_COUNT*kp])
						bestpred=kp;
				}
				predsel[kc]=bestpred;
			}
			for(kt=0;kt<RCT_COUNT;++kt)//select RCT
			{
				const unsigned char *rctinfo=rct_combinations[kt];
				double csize=
					csizes[rctinfo[0]+predsel[rctinfo[0]]*OCH_COUNT]+
					csizes[rctinfo[1]+predsel[rctinfo[1]]*OCH_COUNT]+
					csizes[rctinfo[2]+predsel[rctinfo[2]]*OCH_COUNT];
				if(!kt||bestsize>csize)
					bestsize=csize, bestrct=kt;
			}
			{
				const unsigned char *rctinfo=rct_combinations[bestrct];
				decinfo[4*kb+0]=bestrct;
				decinfo[4*kb+1]=predsel[rctinfo[0]];
				decinfo[4*kb+2]=predsel[rctinfo[1]];
				decinfo[4*kb+3]=predsel[rctinfo[2]];
			}
		}
#ifdef PRINT_PREDS
		{
			int kc, bx, by;
			
			for(by=0;by<yblocks;++by)
			{
				for(bx=0;bx<xblocks;++bx)
					printf(" %-6s", rct_names[decinfo[4*(xblocks*by+bx)+0]]);
				printf("\n");
			}
			printf("\n");
			for(kc=0;kc<3;++kc)
			{
				for(by=0;by<yblocks;++by)
				{
					for(bx=0;bx<xblocks;++bx)
						printf(" %-6s", pred_names[decinfo[4*(xblocks*by+bx)+kc+1]]);
					printf("\n");
				}
				printf("\n");
			}
		}
#else
		(void)och_names;
		(void)rct_names;
		(void)pred_names;
#endif
		free(pixels);
		free(hist);
#if 0
		{
			int kt;

#ifdef PRINT_RCT
			for(kt=0;kt<OCH_COUNT;++kt)
				printf("%d %16lld\n", kt, counters[kt]);
			printf("\n");
#endif
			for(kt=0;kt<_countof(rct_indices);++kt)
			{
				const unsigned char *rctinfo=rct_indices[kt];
				int64_t currerr=
					+counters[rctinfo[0]]+counters[rctinfo[0]+OCH_COUNT]
					+counters[rctinfo[1]]+counters[rctinfo[1]+OCH_COUNT]
					+counters[rctinfo[2]]+counters[rctinfo[2]+OCH_COUNT]
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
		}
		{
			const unsigned char *rctinfo=rct_indices[bestrct];
			int64_t Nerror=
				+counters[rctinfo[0]]
				+counters[rctinfo[1]]
				+counters[rctinfo[2]]
			;
			int64_t Werror=
				+counters[rctinfo[0]+OCH_COUNT]
				+counters[rctinfo[1]+OCH_COUNT]
				+counters[rctinfo[2]+OCH_COUNT]
			;
			NWratio=(int32_t)(((Werror+1)<<8)/(Nerror+Werror+2));
			CLAMP2(NWratio, 1, 255);
		}
#endif

		streamend=dstbuf+usize;
		streamptr=dstbuf;

	//	*streamptr++=bestrct;
	//	*streamptr++=NWratio;
	//	*streamptr++=predidx[0];
	//	*streamptr++=predidx[1];
	//	*streamptr++=predidx[2];
	}
	else
	{
		imptr=dstbuf;
		streamend=srcptr+srcsize;
		streamptr=srcptr;

	//	bestrct=*streamptr++;
	//	NWratio=*streamptr++;
		memcpy(decinfo, streamptr, 4*nblocks);
		streamptr+=4*nblocks;
	//	predidx[0]=*streamptr++;
	//	predidx[1]=*streamptr++;
	//	predidx[2]=*streamptr++;

		code=0;
		code=code<<32|*(uint32_t*)streamptr;
		streamptr+=4;
		code=code<<32|*(uint32_t*)streamptr;
		streamptr+=4;

		csize=srcsize;
	}
	{
		//int
		//	yidx=rct_indices[bestrct][3+0],
		//	uidx=rct_indices[bestrct][3+1],
		//	vidx=rct_indices[bestrct][3+2],
		//	uhelpidx=rct_indices[bestrct][6+0],
		//	vhelpidx=rct_indices[bestrct][6+1];
		int32_t ky, kx, idx;
		int32_t psize=0;
		int16_t *pixels=0;
		int32_t paddedwidth=iw+16;
#ifdef ENABLE_L1
		int32_t coeffs[3][L1PREDS+1]={0};
#endif
		QuantInfoFwd qtablefwd[256];
		QuantInfoInv qtableinv[33];
		//int32_t coeffs2[3][4]={0};
		
		if(fwd)//init qtables
		{
			int ks;

			for(ks=0;ks<256;++ks)
			{
				QuantInfoFwd *info=qtablefwd+ks;

				if(ks<(1<<CONFIG_EXP))
				{
					info->nbypass=0;
					info->token=ks;
					info->bypass=0;
				}
				else
				{
					info->nbypass=floor_log2(ks)-(CONFIG_LSB+CONFIG_MSB);
					info->token=(1<<CONFIG_EXP)-((CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB))<<(CONFIG_MSB+CONFIG_LSB)) + (
						info->nbypass<<(CONFIG_MSB+CONFIG_LSB)|
						(ks>>info->nbypass&((1<<CONFIG_MSB)-1)<<CONFIG_LSB)|
						(ks&((1<<CONFIG_LSB)-1))
					);
					info->bypass=ks>>CONFIG_LSB&((1LL<<info->nbypass)-1);
				}
			}
		}
		else
		{
			int kt;

			for(kt=0;kt<33;++kt)
			{
				QuantInfoInv *info=qtableinv+kt;
				if(kt<(1<<CONFIG_EXP))
				{
					info->sym=kt;
					info->nbits=0;
				}
				else
				{
					info->nbits=(kt>>(CONFIG_MSB+CONFIG_LSB))-((1<<(CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB)))-(CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB)));
					info->sym=
						+(((1<<(CONFIG_MSB+CONFIG_LSB))+(kt&((1<<CONFIG_MSB)-1)<<CONFIG_LSB))<<info->nbits)//7 instructions
					//	+(bypass<<CONFIG_LSB)
						+(kt&((1<<CONFIG_LSB)-1))
					;
				}
			}
		}
		psize=(int32_t)sizeof(int16_t[4*3*2])*paddedwidth;//4 padded rows * 3 channels * {pixels, nbypass}
		pixels=(int16_t*)malloc(psize);
		if(!pixels)
		{
			CRASH("Alloc error\n");
			free(srcbuf);
			free(dstbuf);
			return 1;
		}
		memset(pixels, 0, psize);
		memset(stats, 0, sizeof(stats));
		{
			int kctx, ks;

			for(kctx=0;kctx<NCTX;++kctx)
			{
				for(ks=0;ks<NLEVELS;++ks)//init bypass
					stats[0][0][kctx][ks]=ks*(((1<<PROBBITS_STORE)-(NLEVELS<<(PROBBITS_STORE-PROBBITS_USE)))/NLEVELS);
				stats[0][0][kctx][NLEVELS]=1<<PROBBITS_STORE;
			}
		}
		//const int LOL_1=sizeof(stats[0][0]);
		//const int LOL_2=sizeof(stats);
		memfill(stats[0][1], stats[0][0], sizeof(stats)-sizeof(stats[0][0]), sizeof(stats[0][0]));
		//memcpy(stats[1], stats[0], sizeof(stats[0]));
		//memcpy(stats[2], stats[0], sizeof(stats[0]));
		memset(mixin, 0, sizeof(mixin));//
		FILLMEM(
			(uint32_t*)mixin,
			(1<<RATE>>1),
			sizeof(int32_t[NLEVELS]),
			sizeof(int32_t)
		);
		FILLMEM(
			(uint32_t*)mixin+NLEVELS,
			(1<<PROBBITS_STORE)-(NLEVELS<<(PROBBITS_STORE-PROBBITS_USE))+(1<<RATE>>1),
			sizeof(int32_t[NLEVELS]),
			sizeof(int32_t)
		);
		FILLMEM((int32_t*)coeffs, (1<<L1SH)/L1PREDS, sizeof(coeffs), sizeof(int32_t));
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
			int by=ky/BLOCKY;
			for(kx=0;kx<iw;++kx, ++idx)
			{
				int bx=kx/BLOCKX;
				uint8_t *pidx=decinfo+4*(xblocks*by+bx);
				const uint8_t *combination=rct_combinations[pidx[0]];
				int
					yidx	=combination[II_PERM_Y],
					uidx	=combination[II_PERM_U],
					vidx	=combination[II_PERM_V],
					ufromy	=combination[II_COEFF_U_SUB_Y],
					vc0	=combination[II_COEFF_V_SUB_Y],
					vc1	=combination[II_COEFF_V_SUB_U];
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
					int och=combination[kc];
					int32_t
						NNNN	=rows[0][0+0*3*2],
						NNN	=rows[3][0+0*3*2],
						NNNE	=rows[3][0+1*3*2],
						NNW	=rows[2][0-1*3*2],
						NN	=rows[2][0+0*3*2],
						NNE	=rows[2][0+1*3*2],
						NWW	=rows[1][0-2*3*2],
						NW	=rows[1][0-1*3*2],
						N	=rows[1][0+0*3*2],
						NE	=rows[1][0+1*3*2],
						NEE	=rows[1][0+2*3*2],
						NEEE	=rows[1][0+3*3*2],
						NEEEE	=rows[1][0+4*3*2],
						WWWWWW	=rows[0][0-6*3*2],
						WWWWW	=rows[0][0-5*3*2],
						WWWW	=rows[0][0-4*3*2],
						WWW	=rows[0][0-3*3*2],
						WW	=rows[0][0-2*3*2],
						W	=rows[0][0-1*3*2],
						eNEE	=rows[1][1+2*3*2],
						eNEEE	=rows[1][1+3*3*2],
						eW	=rows[0][1-1*3*2];
					//int32_t p1, p2;
					int32_t vmin, vmax, vmin0, vmax0, vmin1, vmax1, p0;
					int32_t pred, ctx, pctx, error, token, cdf, freq, nbypass, bypass=0;
					uint32_t *currstats, *currstats2, *currstats3, *currstats4;
					int32_t preds[L1PREDS+1];
#if 0
					int32_t preds[]=
					{
#if 1
						N,
						W,
						3*(N-NN)+NNN,
						3*(W-WW)+WWW,
						4*(W+WWW)-6*WW-WWWW,
						4*(N+NNN)-6*NN-NNNN,
					//	5*(W+WWWW)-10*(WW-WWW)+WWWWW,//bad
					//	(21*W-35*WW+35*WWW-21*WWWW+7*WWWWW-WWWWWW+3*(NE-NNE)+NNNE)/7,
						N+W-NW,
						W+NE-N,
						N+NE-NNE,
						(WWWW+NNN+NEEE+NEEEE)/4,
						NEEE,
					//	0,
#endif
#if 0
						NNW-NN,
						NN-N,
						NNE-NN,
						NWW-NW,
						NW-N,
						N,
						NE-N,
						NEE-NE,
						NEEE-NEE,
						WW-W,
						W,
#endif
#if 0
					//	0,
						N,
						W,
						3*(N-NN)+NNN,
						3*(W-WW)+WWW,
						N+W-NW,
						W+NE-N,
						W+NEE-NE,
						N+NE-NNE,
						(WWWW+NNN+NEEE+NEEEE)/4,
						NWratio>128?W:(NWratio<128?N:(N+W)/2),
						NEEE,
#endif
#if 0
						0,
						W,
						(N+W)/2,
						abs(N-NW)>abs(W-NW)?N:W,
						NW,
						N+W-NW,
						(5*(N+W)-2*NW)/8,
						(3*(N+W)-2*NW)/4,
						W+(6*N-5*NW-NN-WW+NE)/8,
						(NN+NNE+NW+N+NE+NEE+WW+W)/8,
						W+(10*N+9*NW+4*NE-2*(NN+WW)+NNW-(NNE+NWW))/16,
						0,
						3*(N-NN)+NNN,
						3*(W-WW)+WWW,
#endif
					};
#endif
#ifdef ENABLE_VAL
					uint64_t val_low, val_range, val_code;
#endif
					(void)NNNN	;
					(void)NNN	;
					(void)NNNE	;
					(void)NNW	;
					(void)NN	;
					(void)NNE	;
					(void)NWW	;
					(void)NW	;
					(void)N		;
					(void)NE	;
					(void)NEE	;
					(void)NEEE	;
					(void)NEEEE	;
					(void)WWWWWW	;
					(void)WWWWW	;
					(void)WWWW	;
					(void)WWW	;
					(void)WW	;
					(void)W		;
					(void)eNEE	;
					(void)eNEEE	;
					(void)eW	;

					pred=0;
					vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					vmin0=vmin;
					vmax0=vmax;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					vmin1=vmin;
					vmax1=vmax;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					
					//if(ky==ih/2&&kx==iw/2)//
					//	printf("");

					{
						int j=0;
						preds[j++]=N;
						preds[j++]=W;
						preds[j++]=3*(N-NN)+NNN;
						preds[j++]=3*(W-WW)+WWW;
						preds[j++]=4*(W+WWW)-6*WW-WWWW;
						preds[j++]=4*(N+NNN)-6*NN-NNNN;
						preds[j++]=N+W-NW;
						preds[j++]=W+NE-N;
						preds[j++]=N+NE-NNE;
						preds[j++]=(WWWW+NNN+NEEE+NEEEE)/4;
						preds[j++]=NEEE;
						p0=1<<L1SH>>1;
						for(j=0;j<L1PREDS+1;++j)
							p0+=coeffs[kc][j]*preds[j];
						p0>>=L1SH;
					}
					switch(pidx[kc+1])
					{
					case 0:pred=W;break;
					case 1:pred=(N+W)/2;break;
					case 2:pred=abs(N-NW)>abs(W-NW)?N:W;break;
					case 3:pred=NW;break;
					case 4:
						pred=N+W-NW;
						CLAMP2(pred, vmin0, vmax0);
						break;
					case 5:
						pred=(5*(N+W)-2*NW)/8;
						CLAMP2(pred, vmin1, vmax1);
						break;
					case 6:
						pred=(3*(N+W)-2*NW)/4;
						CLAMP2(pred, vmin1, vmax1);
						break;
					case 7:
						pred=(4*(N+W)+NE-NW)/8;
						CLAMP2(pred, vmin1, vmax1);
						break;
					case 8:
						pred=W+(6*N-5*NW-NN-WW+NE)/8;
						CLAMP2(pred, vmin1, vmax1);
						break;
					case 9:pred=(NN+NNE+NW+N+NE+NEE+WW+W)/8;break;
					case 10:
						pred=W+(10*N-9*NW+4*NE-2*(NN+WW)+NNW-NNE-NWW)/16;
						CLAMP2(pred, vmin1, vmax1);
						break;
					case 11:
						{
							int gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1;
							int gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1;
							pred=(gx*N+gy*W)/(gx+gy);
						}
						break;
					case 12:
						pred=p0;
						CLAMP2(pred, vmin, vmax);
						//pred=(
						//	+coeffs2[kc][0]*NW
						//	+coeffs2[kc][1]*N
						//	+coeffs2[kc][2]*NE
						//	+coeffs2[kc][3]*W
						//	+(1<<16>>1)
						//)>>16;
						//CLAMP2(pred, vmin, vmax);
						break;
					}
				//	preds[L1PREDS-1]=pred;
					//p1=pred;
					//pred=(
					//	+coeffs[kc][0]*NW
					//	+coeffs[kc][1]*N
					//	+coeffs[kc][2]*NE
					//	+coeffs[kc][3]*W
					//	+coeffs[kc][4]*p1
					//	+(1<<L1SH>>1)
					//)>>L1SH;
					//p2=pred;
					//CLAMP2(pred, vmin, vmax);
					if(kc)
					{
						pred+=offset;
						CLAMP2(pred, -128, 127);
					}
					ctx=floor_log2(eW*eW+1);
					if(ctx>NCTX-1)
						ctx=NCTX-1;
					pctx=pred>>(8-PBITS)&((1<<PBITS)-1);
					{
						int ctx2=ctx+(ctx<NCTX-1);
						int pctx2=pctx+(pctx<(1<<PBITS)-1);
						currstats=stats[och][pctx][ctx];
						currstats2=stats[och][pctx][ctx2];
						currstats3=stats[och][pctx2][ctx];
						currstats4=stats[och][pctx2][ctx2];
					}
	#define GETCDF(X) (((currstats[X]+currstats2[X]+currstats3[X]+currstats4[X])>>(PROBBITS_STORE+2-PROBBITS_USE))+(X))
//	#define GETCDF(X) (((2*currstats[X]+currstats2[X]+currstats3[X])>>(PROBBITS_STORE+2-PROBBITS_USE))+(X))
//	#define GETCDF(X) (((currstats[X]+currstats2[X])>>(PROBBITS_STORE+1-PROBBITS_USE))+(X))
					if(fwd)
					{
						QuantInfoFwd *info;

						error=(char)(yuv[kc]-pred);
#ifdef SAVE_RESIDUALS
						residuals[3*(iw*ky+kx)+kc]=(uint8_t)(error+128);
#endif
						error=error<<1^error>>31;

						info=qtablefwd+error;
						token=info->token;
						nbypass=info->nbypass;
						bypass=info->bypass;
#ifdef _DEBUG
						if(bypass>=(1<<nbypass))
							CRASH("");
#endif
					}
					else
					{
						error=0;
						token=0;
					}
					//if(ky==13&&kx==785)//
					//if(ky==10&&kx==iw/2)//
					//if(ky==0&&kx==2)//
					//if(ky==0&&kx==17&&kc==1)//
					//if(ky==13&&kx==2039&&kc==2)//
					//if(ky==31&&kx==2037&&kc==0)//
					//	printf("");

					//token
					if(range<0x10000)
					{
						if(streamptr>=streamend)
						{
							int symidx=3*idx+kc, totalsyms=3*iw*ih;
							
							CRASH("inflation %d/%d  %8.4lf%%\n"
								, (int32_t)symidx
								, (int32_t)totalsyms
								, 100.*totalsyms/symidx
							);
							return 1;
						}
						if(fwd)
							*(uint32_t*)streamptr=(uint32_t)(low>>32);
						else
							code=code<<32|*(uint32_t*)streamptr;
						streamptr+=4;
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
					}
					if(fwd)
					{
						cdf =GETCDF(token);
						freq=GETCDF(token+1);
					}
					else
					{
						int32_t c=(int32_t)(((code-low)<<PROBBITS_USE|((1<<PROBBITS_USE)-1))/range);
						token=0;
						cdf=0;
						for(;;)
						{
							freq=GETCDF(token+1);
							if(freq>c)
								break;
#ifdef _DEBUG
							if(token>=NLEVELS)
								CRASH("");
#endif
							++token;
							cdf=freq;
						}
						{
							QuantInfoInv *info=qtableinv+token;
							nbypass=info->nbits;
							error=info->sym;
						}
					}
					freq-=cdf;
#ifdef _DEBUG
					if(
						freq<0
						||(uint32_t)cdf>=(unsigned)(1<<PROBBITS_USE)
						||(uint32_t)(cdf+freq)>(uint32_t)(1<<PROBBITS_USE)
					)
						CRASH("Invalid stats");
					{
						uint64_t low0=low, range0=range;
#endif
#ifdef ENABLE_VAL
					val_low=low;
					val_range=range;
					val_code=code;
#endif
					low+=range*cdf>>PROBBITS_USE;
					range=(range*freq>>PROBBITS_USE)-1;
#ifdef ENABLE_VAL
					if(fwd)
						val_enc(cdf, freq, val_low, val_range, low, range);//
					else
						VAL_DEC(cdf, freq, val_low, val_range, val_code, low, range, code, "XY  %d %d", kx, ky);
#endif
#ifdef _DEBUG
					if(!fwd&&(code<low||code>low+range||low<low0||low+range>low0+range0||range>=range0))
						CRASH("");
					}
#endif
					{
						int ks;
						uint32_t *currmixin=mixin+NLEVELS-1-token;

						for(ks=1;ks<NLEVELS;++ks)
						{
							uint32_t cell;

							cell=currstats[ks];
							cell+=(int32_t)(currmixin[ks]-cell)>>RATE;
							currstats[ks]=cell;

							cell=currstats2[ks];
							cell+=(int32_t)(currmixin[ks]-cell)>>RATE;
							currstats2[ks]=cell;

							cell=currstats3[ks];
							cell+=(int32_t)(currmixin[ks]-cell)>>RATE;
							currstats3[ks]=cell;

							cell=currstats4[ks];
							cell+=(int32_t)(currmixin[ks]-cell)>>RATE;
							currstats4[ks]=cell;
//#ifdef _DEBUG
//							if(cell>(1<<PROBBITS_STORE)-NLEVELS||(currstats!=currstats2&&ks&&cell<currstats4[ks-1]))//
//								CRASH("");
//#endif
						}
					}

					if(nbypass)//bypass
					{
						if(range<0x10000)
						{
							if(streamptr>=streamend)
							{
								int symidx=3*idx+kc, totalsyms=3*iw*ih;

								CRASH("inflation %d/%d  %8.4lf%%\n"
									, (int32_t)symidx
									, (int32_t)totalsyms
									, 100.*totalsyms/symidx
								);
								return 1;
							}
							if(fwd)
								*(uint32_t*)streamptr=(uint32_t)(low>>32);
							else
								code=code<<32|*(uint32_t*)streamptr;
							streamptr+=4;
							low<<=32;
							range=range<<32|0xFFFFFFFF;
							if(range>~low)
								range=~low;
						}
						if(!fwd)
						{
							bypass=(int32_t)((((code-low+1)<<nbypass)-1)/range);
#ifdef _DEBUG
							if(bypass>=(1<<nbypass))
								CRASH("");
#endif
						}
#ifdef _DEBUG
						{
							uint64_t low0=low, range0=range;
#endif
#ifdef ENABLE_VAL
						val_low=low;
						val_range=range;
						val_code=code;
#endif
						low+=range*bypass>>nbypass;
						range=(range>>nbypass)-1;
#ifdef ENABLE_VAL
						if(fwd)
							val_enc(bypass, 1, val_low, val_range, low, range);//
						else
							VAL_DEC(bypass, 1, val_low, val_range, val_code, low, range, code, "XY  %d %d", kx, ky);
#endif
#ifdef _DEBUG
							if(!fwd&&(code<low||code>low+range||low<low0||low+range>low0+range0||range>=range0))
								CRASH("");
						}
#endif
					}
					else
						bypass=0;

					if(!fwd)
					{
						error+=bypass<<CONFIG_LSB;
						error=error>>1^-(error&1);
						yuv[kc]=(char)(error+pred);
					}
					rows[0][0]=yuv[kc]-offset;
					error=yuv[kc]-pred;
					rows[0][1]=(2*eW+((error<<1^error>>31)<<3)+(eNEE>eNEEE?eNEE:eNEEE))>>2;
					//{
					//	int32_t e=rows[0][0];
					//	e=(e>p2)-(e<p2);
					//	coeffs[kc][0]+=e*NW;
					//	coeffs[kc][1]+=e*N;
					//	coeffs[kc][2]+=e*NE;
					//	coeffs[kc][3]+=e*W;
					//	coeffs[kc][4]+=e*p1;
					//}
#ifdef ENABLE_L1
				//	if(pidx[kc+1]==NPREDS-1)
					{
						int32_t e=rows[0][0], k;
						e=(e>p0)-(e<p0);
						for(k=0;k<L1PREDS+1;++k)
							coeffs[kc][k]+=e*preds[k];
					}
#endif
					offset=(kc?vc0*yuv[0]+vc1*yuv[1]:ufromy*yuv[0])/4;
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
					guide_check(dstbuf, kx, ky);
				}
				imptr+=3;
			}
		}
		free(pixels);
	}
	{
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing\n", dstfn);
			free(srcbuf);
			free(dstbuf);
			return 1;
		}
		if(fwd)
		{
			*(uint32_t*)streamptr=(uint32_t)(low>>32);
			streamptr+=4;
			*(uint32_t*)streamptr=(uint32_t)low;
			streamptr+=4;

			csize=streamptr-dstbuf;

			dstsize+=fwrite("07", 1, 2, fdst);
			dstsize+=fwrite(&iw, 1, 4, fdst);
			dstsize+=fwrite(&ih, 1, 4, fdst);
			dstsize+=fwrite(decinfo, 1, nblocks*4, fdst);
			dstsize+=fwrite(dstbuf, 1, csize, fdst);
		}
		else
		{
			dstsize+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			dstsize+=fwrite(dstbuf, 1, usize, fdst);
		}
		fclose(fdst);
	}
	free(srcbuf);
	free(dstbuf);
#ifdef SAVE_RESIDUALS
	if(fwd)
	{
		const char fn[]="20250530_1220am.ppm";
		FILE *fdst=fopen(fn, "wb");
		if(!fdst)
		{
			CRASH("Cannot open \"%s\" for writing", fn);
			return 1;
		}
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		fwrite(residuals, 1, (ptrdiff_t)3*iw*ih, fdst);
		fclose(fdst);
		free(residuals);
	}
#endif
#ifdef LOUD
	t=time_sec()-t;
	if(fwd)
	{
		printf("%9d->%9d  %8.4lf%%  %12.6lf\n"
			, srcsize
			, dstsize
			, 100.*dstsize/srcsize
			, (double)srcsize/dstsize
		);
	}
	printf("%c  %12.6lf  %12.6lf MB/s\n"
		, 'D'+fwd
		, t
		, usize/(t*1024*1024)
	);
#endif
	(void)time_sec;
	return 0;
}