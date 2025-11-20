#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#endif
#include<sys/stat.h>
#if defined _WIN32 || defined WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>
#else
#include<time.h>
#endif


#ifdef _MSC_VER
	#define ENABLE_GUIDE

	#define ESTIMATE_SIZE
#endif


#define NPREDS 8
#define SHIFT 20

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
#ifdef _MSC_VER
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
#endif
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
		CRASH("XY %d %d  RGB 0x%08X 0x%08X OG", kx, ky, *(unsigned*)(image+idx), *(unsigned*)(g_image+idx));
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif

	#define ENABLE_EXTENDED_RCT
#if 1
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

	II_COUNT,
} RCTInfoIdx;
#ifdef ENABLE_EXTENDED_RCT
#define OCHLIST\
	OCH(YX00) OCH(Y0X0) OCH(Y00X)\
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
	OCH_R=OCH_YX00,
	OCH_G=OCH_Y0X0,
	OCH_B=OCH_Y00X,
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
#endif
#ifndef ENABLE_EXTENDED_RCT
typedef enum _OCHIndex
{
	OCH_R,
	OCH_G,
	OCH_B,
	OCH_RG,
	OCH_GB,
	OCH_BR,
	OCH_COUNT,
	OCH_GR=OCH_RG,
	OCH_BG=OCH_GB,
	OCH_RB=OCH_BR,
} OCHIndex;
#endif
#ifdef ENABLE_EXTENDED_RCT
#define RCTLIST\
	RCT(_X00_0X0_00X,	OCH_R,		OCH_G,		OCH_B,		0, 1, 2,	0,  0, 0)\
	RCT(_X00_0X0_04X,	OCH_R,		OCH_G,		OCH_BG,		0, 1, 2,	0,  0, 4)\
	RCT(_X00_0X0_40X,	OCH_R,		OCH_G,		OCH_BR,		0, 1, 2,	0,  4, 0)\
	RCT(_0X0_00X_X40,	OCH_G,		OCH_B,		OCH_RG,		1, 2, 0,	0,  4, 0)\
	RCT(_0X0_00X_X04,	OCH_G,		OCH_B,		OCH_RB,		1, 2, 0,	0,  0, 4)\
	RCT(_00X_X00_4X0,	OCH_B,		OCH_R,		OCH_GR,		2, 0, 1,	0,  0, 4)\
	RCT(_00X_X00_0X4,	OCH_B,		OCH_R,		OCH_GB,		2, 0, 1,	0,  4, 0)\
	RCT(_0X0_04X_X40,	OCH_G,		OCH_BG,		OCH_RG,		1, 2, 0,	4,  4, 0)\
	RCT(_0X0_04X_X04,	OCH_G,		OCH_BG,		OCH_RB,		1, 2, 0,	4,  0, 4)\
	RCT(_0X0_X40_40X,	OCH_G,		OCH_RG,		OCH_BR,		1, 0, 2,	4,  0, 4)\
	RCT(_00X_X04_0X4,	OCH_B,		OCH_RB,		OCH_GB,		2, 0, 1,	4,  4, 0)\
	RCT(_00X_X04_4X0,	OCH_B,		OCH_RB,		OCH_GR,		2, 0, 1,	4,  0, 4)\
	RCT(_00X_0X4_X40,	OCH_B,		OCH_GB,		OCH_RG,		2, 1, 0,	4,  0, 4)\
	RCT(_X00_4X0_40X,	OCH_R,		OCH_GR,		OCH_BR,		0, 1, 2,	4,  4, 0)\
	RCT(_X00_4X0_04X,	OCH_R,		OCH_GR,		OCH_BG,		0, 1, 2,	4,  0, 4)\
	RCT(_X00_40X_0X4,	OCH_R,		OCH_BR,		OCH_GB,		0, 2, 1,	4,  0, 4)\
	RCT(_X00_0X0_13X,	OCH_R,		OCH_G,		OCH_B1,		0, 1, 2,	0,  1, 3)\
	RCT(_X00_4X0_13X,	OCH_R,		OCH_GR,		OCH_B1,		0, 1, 2,	4,  1, 3)\
	RCT(_X00_00X_3X1,	OCH_R,		OCH_B,		OCH_G1,		0, 2, 1,	0,  3, 1)\
	RCT(_X00_40X_3X1,	OCH_R,		OCH_BR,		OCH_G1,		0, 2, 1,	4,  3, 1)\
	RCT(_0X0_00X_X13,	OCH_G,		OCH_B,		OCH_R1,		1, 2, 0,	0,  1, 3)\
	RCT(_0X0_04X_X13,	OCH_G,		OCH_BG,		OCH_R1,		1, 2, 0,	4,  1, 3)\
	RCT(_0X0_X40_13X,	OCH_G,		OCH_RG,		OCH_B1,		1, 0, 2,	4,  3, 1)\
	RCT(_00X_X04_3X1,	OCH_B,		OCH_RB,		OCH_G1,		2, 0, 1,	4,  1, 3)\
	RCT(_00X_04X_X13,	OCH_B,		OCH_GB,		OCH_R1,		2, 1, 0,	4,  3, 1)\
	RCT(_X00_0X0_22X,	OCH_R,		OCH_G,		OCH_B2,		0, 1, 2,	0,  2, 2)\
	RCT(_X00_4X0_22X,	OCH_R,		OCH_GR,		OCH_B2,		0, 1, 2,	4,  2, 2)\
	RCT(_X00_00X_2X2,	OCH_R,		OCH_B,		OCH_G2,		0, 2, 1,	0,  2, 2)\
	RCT(_X00_40X_2X2,	OCH_R,		OCH_BR,		OCH_G2,		0, 2, 1,	4,  2, 2)\
	RCT(_0X0_00X_X22,	OCH_G,		OCH_B,		OCH_R2,		1, 2, 0,	0,  2, 2)\
	RCT(_0X0_04X_X22,	OCH_G,		OCH_BG,		OCH_R2,		1, 2, 0,	4,  2, 2)\
	RCT(_0X0_X40_22X,	OCH_G,		OCH_RG,		OCH_B2,		1, 0, 2,	4,  2, 2)\
	RCT(_00X_X04_2X2,	OCH_B,		OCH_RB,		OCH_G2,		2, 0, 1,	4,  2, 2)\
	RCT(_00X_0X4_X22,	OCH_B,		OCH_GB,		OCH_R2,		2, 1, 0,	4,  2, 2)\
	RCT(_X00_0X0_31X,	OCH_R,		OCH_G,		OCH_B3,		0, 1, 2,	0,  3, 1)\
	RCT(_X00_4X0_31X,	OCH_R,		OCH_GR,		OCH_B3,		0, 1, 2,	4,  3, 1)\
	RCT(_X00_00X_1X3,	OCH_R,		OCH_B,		OCH_G3,		0, 2, 1,	0,  1, 3)\
	RCT(_X00_40X_1X3,	OCH_R,		OCH_BR,		OCH_G3,		0, 2, 1,	4,  1, 3)\
	RCT(_0X0_00X_X31,	OCH_G,		OCH_B,		OCH_R3,		1, 2, 0,	0,  3, 1)\
	RCT(_0X0_04X_X31,	OCH_G,		OCH_BG,		OCH_R3,		1, 2, 0,	4,  3, 1)\
	RCT(_0X0_X40_31X,	OCH_G,		OCH_RG,		OCH_B3,		1, 0, 2,	4,  1, 3)\
	RCT(_00X_X04_1X3,	OCH_B,		OCH_RB,		OCH_G3,		2, 0, 1,	4,  3, 1)\
	RCT(_00X_0X4_X31,	OCH_B,		OCH_GB,		OCH_R3,		2, 1, 0,	4,  1, 3)
#endif
#ifndef ENABLE_EXTENDED_RCT
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
#endif
static int crct_analysis(unsigned char *image, int iw, int ih)
{
	long long counters[OCH_COUNT]={0};
	int prev[OCH_COUNT]={0};
	for(ptrdiff_t k=0, len=(ptrdiff_t)3*iw*ih;k<len;k+=3)
	{
		int
			r=image[k+0]<<2,
			g=image[k+1]<<2,
			b=image[k+2]<<2,
			rg=r-g,
			gb=g-b,
			br=b-r;
		counters[0]+=abs(r -prev[0]);
		counters[1]+=abs(g -prev[1]);
		counters[2]+=abs(b -prev[2]);
		counters[3]+=abs(rg-prev[3]);
		counters[4]+=abs(gb-prev[4]);
		counters[5]+=abs(br-prev[5]);
		prev[0]=r;
		prev[1]=g;
		prev[2]=b;
		prev[3]=rg;
		prev[4]=gb;
		prev[5]=br;
#ifdef ENABLE_EXTENDED_RCT
#define UPDATE(IDXA, IDXB, IDXC, A0, B0, C0)\
	do\
	{\
		int a0=A0, b0=B0, c0=C0;\
		counters[IDXA]+=abs(a0-prev[IDXA]);\
		counters[IDXB]+=abs(b0-prev[IDXB]);\
		counters[IDXC]+=abs(c0-prev[IDXC]);\
		prev[IDXA]=a0;\
		prev[IDXB]=b0;\
		prev[IDXC]=c0;\
	}while(0)
		//r-(3*g+b)/4 = r-g-(b-g)/4
		//g-(3*r+b)/4 = g-r-(b-r)/4
		//b-(3*r+g)/4 = b-r-(g-r)/4
		UPDATE(OCH_CX31, OCH_C3X1, OCH_C31X, rg+(gb>>2), rg+(br>>2), br+(rg>>2));

		//r-(g+3*b)/4 = r-b-(g-b)/4
		//g-(r+3*b)/4 = g-b-(r-b)/4
		//b-(r+3*g)/4 = b-g-(r-g)/4
		UPDATE(OCH_CX13, OCH_C1X3, OCH_C13X, br+(gb>>2), gb+(br>>2), gb+(rg>>2));

		//r-(g+b)/2 = (r-g + r-b)/2
		//g-(r+b)/2 = (g-r + g-b)/2
		//b-(r+g)/2 = (b-r + b-g)/2
		UPDATE(OCH_CX22, OCH_C2X2, OCH_C22X, (rg-br)>>1, (gb-rg)>>1, (br-gb)>>1);
#undef  UPDATE
#endif
	}
	int bestrct=0;
	long long minerr=0;
	for(int kt=0;kt<RCT_COUNT;++kt)
	{
		const unsigned char *rct=rct_combinations[kt];
		long long currerr=
			+counters[rct[0]]
			+counters[rct[1]]
			+counters[rct[2]]
		;
		if(!kt||minerr>currerr)
		{
			minerr=currerr;
			bestrct=kt;
		}
	}
	return bestrct;
}

#define NCTX 18
static unsigned hist[3*NCTX*257], CDF[3*NCTX*257];
static void update_CDF(const unsigned *hist, unsigned *CDF, const int nlevels)
{
	int hsum=0;
	for(int ks=0;ks<nlevels;++ks)//faster than maintaining hist sum
		hsum+=hist[ks];
	long long rsum=(((1LL<<16)-nlevels)<<24)/hsum;//adaptive: allow all symbols
	for(int ks=0, c=0;ks<nlevels;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*rsum>>24)+ks;
		c+=freq;
	}
}
int c39_codec(int argc, char **argv)
{
	if(argc!=3)
	{
		printf(
			"Usage: \"%s\"  input  output    Encode/decode.\n"
			, argv[0]
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];
#ifdef ESTIMATE_SIZE
	double t=time_sec();
	double csizes[3]={0};
#endif
	if(!srcfn||!dstfn)
	{
		CRASH("Codec requires both source and destination filenames");
		return 1;
	}
	FILE *fsrc=fopen(srcfn, "rb");
	int tag=0;
	size_t nread=fread(&tag, 1, 2, fsrc), nwritten=0;
	if(nread!=2)
	{
		CRASH("File is empty");
		return 1;
	}
	int fwd=tag==('P'|'6'<<8);
	if(!fwd&&tag!=('3'|'9'<<8))
	{
		CRASH("Unsupported source file");
		return 1;
	}
	int iw=0, ih=0;
	int usize=0;
	int bestrct=0;
	if(fwd)
	{
		int temp=fgetc(fsrc);
		if(temp!='\n')
		{
			CRASH("Invalid PPM file");
			return 1;
		}
		nread=fscanf(fsrc, "%d %d", &iw, &ih);
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
		iw=0;
		ih=0;
		nread=fread(&iw, 1, 3, fsrc);
		nread+=fread(&ih, 1, 3, fsrc);
		nread+=fread(&bestrct, 1, 1, fsrc);
		if(nread!=3LL+3+1)
		{
			CRASH("Unsupported archive");
			return 1;
		}
	}
	usize=iw*ih*3;
	unsigned char *buffer=(unsigned char*)malloc(usize);
	if(!buffer)
	{
		CRASH("Alloc error");
		return 1;
	}
	FILE *fdst=fopen(dstfn, "wb");
	nwritten=0;
	if(fwd)
	{
		nread=fread(buffer, 1, usize, fsrc);
		if(nread!=usize)
		{
			CRASH("Corrupt PPM");
			return 1;
		}
		guide_save(buffer, iw, ih);
		bestrct=crct_analysis(buffer, iw, ih);

		nwritten+=fwrite("39", 1, 2, fdst);
		nwritten+=fwrite(&iw, 1, 3, fdst);
		nwritten+=fwrite(&ih, 1, 3, fdst);
		nwritten+=fwrite(&bestrct, 1, 1, fdst);
	}
	else
	{
		nwritten+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		memset(buffer, 0, usize);
	}
	int padw=iw+8*2;
	int psize=padw*(int)sizeof(int[4*4*2]);//4 padded rows * 4 channels max * {pixel, error}
	int *pixels=(int*)_mm_malloc(psize, sizeof(__m128i));
	int sbufsize=iw*16;//rowsize = iw*3, allocate iw*16 just in case
	uint16_t *sbuf=(uint16_t*)malloc(sbufsize);//stream buffer
	if(!pixels||!sbuf)
	{
		CRASH("Alloc error");
		return 1;
	}
	memset(sbuf, 0, sbufsize);
	memset(hist, 0, sizeof(hist));
	for(int k=0;k<257;++k)
		CDF[k]=k<<8;
	for(int k=1;k<3*NCTX;++k)
		memcpy(CDF+k*257, CDF, sizeof(int[257]));
	memset(pixels, 0, psize);
	int rowsize=iw*3;
	int yidx=rct_combinations[bestrct][II_PERM_Y];
	int uidx=rct_combinations[bestrct][II_PERM_U];
	int vidx=rct_combinations[bestrct][II_PERM_V];
	int umask=-(rct_combinations[bestrct][II_COEFF_U_SUB_Y]!=0);
	int vc0=rct_combinations[bestrct][II_COEFF_V_SUB_Y];
	int vc1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
	__m128i weights[NPREDS];
	for(int k=0;k<4*NPREDS;++k)
		((int*)weights)[k]=(1<<SHIFT)/NPREDS;
	char modified[3*NCTX]={0};
	for(int ky=0, idx=0;ky<ih;++ky)
	{
		ALIGN(32) int *rows[]=
		{
			pixels+(padw*((ky-0LL)&3)+8LL)*4*2,
			pixels+(padw*((ky-1LL)&3)+8LL)*4*2,
			pixels+(padw*((ky-2LL)&3)+8LL)*4*2,
			pixels+(padw*((ky-3LL)&3)+8LL)*4*2,
		};
		ALIGN(16) int ctx[4], preds[4]={0};
		unsigned *curr_CDF=0;
		unsigned cdfs[3]={0}, freqs[3]={0};
		int yuv[3]={0};
		ALIGN(16) int errors[4]={0};
		int offset1, offset2;
		unsigned long long low=0, range=0xFFFFFFFFFFFF, code=0;
#ifdef AC_VALIDATE
		unsigned long long lo0, r0;
#endif
		uint16_t *sptr=sbuf;
		int streamsize=0;
		if(!fwd)
		{
			nread=fread(&streamsize, 1, 4, fsrc);
			if(streamsize)
			{
				nread=fread(sbuf, sizeof(uint16_t), streamsize, fsrc);
				code=*sptr++;
				code=*sptr++|code<<16;
				code=*sptr++|code<<16;
			}
			else
				nread=fread(buffer+rowsize*ky, 1, rowsize, fsrc);
		}
		for(int kx=0;kx<iw;++kx, idx+=3)
		{
			int *curr=rows[0];

			//ctx = FLOOR_LOG2(eW*eW+1)
			__m128i eW	=_mm_load_si128((__m128i*)rows[0]-1*2+1);	//rows[-Y]+X*2+E
			__m128i mc=_mm_add_epi32(_mm_mullo_epi32(eW, eW), _mm_set1_epi32(1));
			mc=_mm_castps_si128(_mm_cvtepi32_ps(mc));
			mc=_mm_srli_epi32(mc, 23);
			mc=_mm_sub_epi32(mc, _mm_set1_epi32(127));
			mc=_mm_min_epi32(mc, _mm_set1_epi32(NCTX-1));
			mc=_mm_add_epi32(mc, _mm_set_epi32(0, NCTX*2, NCTX, 0));
			_mm_store_si128((__m128i*)ctx, mc);

			modified[ctx[0]]=1;
			modified[ctx[1]]=1;
			modified[ctx[2]]=1;

			ctx[0]*=257;
			ctx[1]*=257;
			ctx[2]*=257;

			__m128i mNW	=_mm_load_si128((__m128i*)rows[1]-1*2+0);	//rows[-Y]+X*2+E
			__m128i mN	=_mm_load_si128((__m128i*)rows[1]+0*2+0);
			__m128i mNE	=_mm_load_si128((__m128i*)rows[1]+1*2+0);
			__m128i mNEEE	=_mm_load_si128((__m128i*)rows[1]+3*2+0);
			__m128i mW	=_mm_load_si128((__m128i*)rows[0]-1*2+0);
			__m128i e0=_mm_sub_epi32(_mm_add_epi32(mN, mW), mNW);
			__m128i e1=mN;
			__m128i e2=mW;
			__m128i e3=_mm_sub_epi32(_mm_add_epi32(mW, mNE), mN);
			__m128i e4=_mm_sub_epi32(mN, _mm_load_si128((__m128i*)rows[2]+0*2+0));
			__m128i e5=_mm_sub_epi32(mW, _mm_load_si128((__m128i*)rows[0]-2*2+0));
			e4=_mm_add_epi32(e4, _mm_slli_epi32(e4, 1));
			e5=_mm_add_epi32(e5, _mm_slli_epi32(e5, 1));
			e4=_mm_add_epi32(e4, _mm_load_si128((__m128i*)rows[3]+0*2+0));
			e5=_mm_add_epi32(e5, _mm_load_si128((__m128i*)rows[0]-3*2+0));
			__m128i e6=_mm_sub_epi32(_mm_add_epi32(mN, mNE), _mm_load_si128((__m128i*)rows[2]+1*2+0));
			__m128i e7=_mm_load_si128((__m128i*)rows[1]+2*2+0);
			__m128i t0=_mm_mullo_epi32(weights[0], e0);
			__m128i t1=_mm_mullo_epi32(weights[1], e1);
			__m128i t2=_mm_mullo_epi32(weights[2], e2);
			__m128i t3=_mm_mullo_epi32(weights[3], e3);
			__m128i t4=_mm_mullo_epi32(weights[4], e4);
			__m128i t5=_mm_mullo_epi32(weights[5], e5);
			__m128i t6=_mm_mullo_epi32(weights[6], e6);
			__m128i t7=_mm_mullo_epi32(weights[7], e7);
			__m128i mmin=_mm_min_epi32(mN, mW);
			__m128i mmax=_mm_max_epi32(mN, mW);
			mmin=_mm_min_epi32(mmin, mNE);
			mmax=_mm_max_epi32(mmax, mNE);
			mmin=_mm_min_epi32(mmin, mNEEE);
			mmax=_mm_max_epi32(mmax, mNEEE);
			t0=_mm_add_epi32(t0, t1);
			t2=_mm_add_epi32(t2, t3);
			t4=_mm_add_epi32(t4, t5);
			t6=_mm_add_epi32(t6, t7);
			t0=_mm_add_epi32(t0, t2);
			t4=_mm_add_epi32(t4, t6);
			t0=_mm_add_epi32(t0, t4);
			t0=_mm_add_epi32(t0, _mm_set1_epi32(1<<SHIFT>>1));
			t0=_mm_srai_epi32(t0, SHIFT);
			__m128i mp0=t0;
			t0=_mm_max_epi32(t0, mmin);
			t0=_mm_min_epi32(t0, mmax);
			_mm_store_si128((__m128i*)preds, t0);
			if(fwd||!streamsize)
			{
				yuv[0]=buffer[idx+yidx]-128;//g
				yuv[1]=buffer[idx+uidx]-128;//b
				yuv[2]=buffer[idx+vidx]-128;//r

				offset1=yuv[0]&umask;
				offset2=(vc0*yuv[0]+vc1*yuv[1])>>2;
				preds[1]+=offset1;
				preds[2]+=offset2;
				CLAMP2(preds[1], -128, 127);
				CLAMP2(preds[2], -128, 127);

				curr[0]=yuv[0];
				curr[1]=yuv[1]-offset1;
				curr[2]=yuv[2]-offset2;

				errors[0]=yuv[0]-preds[0];
				errors[1]=yuv[1]-preds[1];
				errors[2]=yuv[2]-preds[2];
				errors[0]=(char)errors[0];
				errors[1]=(char)errors[1];
				errors[2]=(char)errors[2];
				errors[0]=errors[0]<<1^errors[0]>>31;//pack sign
				errors[1]=errors[1]<<1^errors[1]>>31;
				errors[2]=errors[2]<<1^errors[2]>>31;

				if(fwd)
				{
					curr_CDF=CDF+ctx[0];
					cdfs[0]=curr_CDF[errors[0]];
					freqs[0]=curr_CDF[errors[0]+1]-cdfs[0];
					
					curr_CDF=CDF+ctx[1];
					cdfs[1]=curr_CDF[errors[1]];
					freqs[1]=curr_CDF[errors[1]+1]-cdfs[1];
			
					curr_CDF=CDF+ctx[2];
					cdfs[2]=curr_CDF[errors[2]];
					freqs[2]=curr_CDF[errors[2]+1]-cdfs[2];
#ifdef ESTIMATE_SIZE
					csizes[0]-=log2((double)freqs[0]/0x10000);
					csizes[1]-=log2((double)freqs[1]/0x10000);
					csizes[2]-=log2((double)freqs[2]/0x10000);
#endif
					while(range<0x10000)
					{
						*sptr++=(uint16_t)(low>>32);
						low=low<<16&0xFFFFFFFFFFFF;
						range=range<<16|0xFFFF;
						if(range>(low^0xFFFFFFFFFFFF))
							range=low^0xFFFFFFFFFFFF;
					}
#ifdef AC_VALIDATE
					lo0=low; r0=range;
#endif
					low+=range*cdfs[0]>>16;
					range=(range*freqs[0]>>16)-1;
#ifdef AC_VALIDATE
					acval_enc(0, cdfs[0], freqs[0], lo0, lo0+r0, low, low+range, 0, 0);//
#endif

					while(range<0x10000)
					{
						*sptr++=(uint16_t)(low>>32);
						low=low<<16&0xFFFFFFFFFFFF;
						range=range<<16|0xFFFF;
						if(range>(low^0xFFFFFFFFFFFF))
							range=low^0xFFFFFFFFFFFF;
					}
#ifdef AC_VALIDATE
					lo0=low; r0=range;
#endif
					low+=range*cdfs[1]>>16;
					range=(range*freqs[1]>>16)-1;
#ifdef AC_VALIDATE
					acval_enc(0, cdfs[1], freqs[1], lo0, lo0+r0, low, low+range, 0, 0);//
#endif

					while(range<0x10000)
					{
						*sptr++=(uint16_t)(low>>32);
						low=low<<16&0xFFFFFFFFFFFF;
						range=range<<16|0xFFFF;
						if(range>(low^0xFFFFFFFFFFFF))
							range=low^0xFFFFFFFFFFFF;
					}
#ifdef AC_VALIDATE
					lo0=low; r0=range;
#endif
					low+=range*cdfs[2]>>16;
					range=(range*freqs[2]>>16)-1;
#ifdef AC_VALIDATE
					acval_enc(0, cdfs[2], freqs[2], lo0, lo0+r0, low, low+range, 0, 0);//
#endif
				}
			}
			else
			{
				unsigned cdf, freq, sym;
				unsigned long long code2;
				while(!(range>>16))
				{
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
				//c2=(unsigned)(((code-low)<<16|0xFFFF)/range);
				code2=(code-low)<<16|0xFFFF;
				curr_CDF=CDF+ctx[0];
				sym=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym+2];
					//if((unsigned)freq>c2)
					if(range*freq>code2)
					{
						sym+=range*curr_CDF[sym+1]<=code2;
						break;
					}
					sym+=2;
				}
				cdf=curr_CDF[sym];
				freq=curr_CDF[sym+1];
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				freq-=cdf;
				errors[0]=sym;
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
#ifdef AC_VALIDATE
				acval_dec(0, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				
				while(!(range>>16))
				{
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
				//c2=(unsigned)(((code-low)<<16|0xFFFF)/range);
				code2=(code-low)<<16|0xFFFF;
				curr_CDF=CDF+ctx[1];
				sym=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym+2];
					//if((unsigned)freq>c2)
					if(range*freq>code2)
					{
						sym+=range*curr_CDF[sym+1]<=code2;
						break;
					}
					sym+=2;
				}
				cdf=curr_CDF[sym];
				freq=curr_CDF[sym+1];
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				freq-=cdf;
				errors[1]=sym;
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
#ifdef AC_VALIDATE
				acval_dec(0, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				
				while(!(range>>16))
				{
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
				//c2=(unsigned)(((code-low)<<16|0xFFFF)/range);
				code2=(code-low)<<16|0xFFFF;
				curr_CDF=CDF+ctx[2];
				sym=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym+2];
					//if((unsigned)freq>c2)
					if(range*freq>code2)
					{
						sym+=range*curr_CDF[sym+1]<=code2;
						break;
					}
					sym+=2;
				}
				cdf=curr_CDF[sym];
				freq=curr_CDF[sym+1];
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				freq-=cdf;
				errors[2]=sym;
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
#ifdef AC_VALIDATE
				acval_dec(0, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				yuv[0]=errors[0]>>1^-(errors[0]&1);//unpack sign
				yuv[1]=errors[1]>>1^-(errors[1]&1);
				yuv[2]=errors[2]>>1^-(errors[2]&1);
				yuv[0]+=preds[0];
				yuv[0]=(char)yuv[0];
				offset1=yuv[0]&umask;
				preds[1]+=offset1;
				CLAMP2(preds[1], -128, 127);
				yuv[1]+=preds[1];
				yuv[1]=(char)yuv[1];
				offset2=(vc0*yuv[0]+vc1*yuv[1])>>2;
				preds[2]+=offset2;
				CLAMP2(preds[2], -128, 127);
				yuv[2]+=preds[2];
				yuv[2]=(char)yuv[2];
				curr[0]=yuv[0];
				curr[1]=yuv[1]-offset1;
				curr[2]=yuv[2]-offset2;
				buffer[idx+yidx]=(unsigned char)(yuv[0]+128);
				buffer[idx+uidx]=(unsigned char)(yuv[1]+128);
				buffer[idx+vidx]=(unsigned char)(yuv[2]+128);
				guide_check(buffer, kx, ky);
			}

			__m128i me=_mm_sub_epi32(_mm_load_si128((__m128i*)rows[0]+0*2+0), mp0);
			weights[0]=_mm_add_epi32(weights[0], _mm_sign_epi32(e0, me));
			weights[1]=_mm_add_epi32(weights[1], _mm_sign_epi32(e1, me));
			weights[2]=_mm_add_epi32(weights[2], _mm_sign_epi32(e2, me));
			weights[3]=_mm_add_epi32(weights[3], _mm_sign_epi32(e3, me));
			weights[4]=_mm_add_epi32(weights[4], _mm_sign_epi32(e4, me));
			weights[5]=_mm_add_epi32(weights[5], _mm_sign_epi32(e5, me));
			weights[6]=_mm_add_epi32(weights[6], _mm_sign_epi32(e6, me));
			weights[7]=_mm_add_epi32(weights[7], _mm_sign_epi32(e7, me));

			//ecurr = (2*eW + (error<<3) + max(eNEE, eNEEE))>>2
			me=_mm_max_epi32(_mm_load_si128((__m128i*)rows[1]+2*2+1), _mm_load_si128((__m128i*)rows[1]+3*2+1));
			me=_mm_avg_epu16(me, _mm_slli_epi32(_mm_load_si128((__m128i*)errors), 3));
			me=_mm_avg_epu16(me, eW);
			_mm_store_si128((__m128i*)rows[0]+0*2+1, me);

			++hist[errors[0]+ctx[0]];
			++hist[errors[1]+ctx[1]];
			++hist[errors[2]+ctx[2]];
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
		if(fwd)
		{
			*sptr++=(uint16_t)(low>>32);
			*sptr++=(uint16_t)(low>>16);
			*sptr++=(uint16_t)low;

			streamsize=(int)(sptr-sbuf);//number of emits
			if(streamsize<<1>rowsize+4&&0)//bypass row
			{
				streamsize=0;
				fwrite(&streamsize, 1, 4, fdst);
				fwrite(buffer+rowsize*ky, 1, rowsize, fdst);
			}
			else
			{
				fwrite(&streamsize, 1, 4, fdst);
				fwrite(sbuf, sizeof(uint16_t), streamsize, fdst);
			}
		}

		for(int k=0;k<3*NCTX;++k)
		{
			if(modified[k])
			{
				update_CDF(hist+k*257, CDF+k*257, 256);
				modified[k]=0;
			}
		}
	}
	if(!fwd)
		fwrite(buffer, 1, usize, fdst);
	fclose(fsrc);
	fclose(fdst);
	_mm_free(pixels);
	free(sbuf);
#ifdef ESTIMATE_SIZE
	ptrdiff_t dstsize;
	{
		struct stat info={0};
		stat(dstfn, &info);
		dstsize=info.st_size;
	}
	t=time_sec()-t;
	if(fwd)
	{
		csizes[0]/=8;
		csizes[1]/=8;
		csizes[2]/=8;
		printf("%.2lf + %.2lf + %.2lf = %.2lf  ->  %td\n",
			csizes[0],
			csizes[1],
			csizes[2],
			csizes[0]+csizes[1]+csizes[2],
			dstsize
		);
	}
	printf("%c  %lf sec  %lf MB/s\n", 'D'+fwd, t, usize/(t*1024*1024));
#endif
	free(buffer);
	(void)nwritten;
	return 0;
}