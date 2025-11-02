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
static const char file[]=__FILE__;

#define YPAD 4
#define XPAD 4

#define NCTX 18

#define WPSH 16
#define L1SH 20

#define PREDLIST\
	PRED(N)\
	PRED(W)\
	PRED(NW)\
	PRED(NE)\
	PRED(N+W-NW)\
	PRED(W+NE-N)\
	PRED(N+NE-NNE)\
	PRED(N+NW-NNW)\
	PRED(W+NW-NWW)\
	PRED((WWWW+NNNN)>>1)\
	PRED(3*(N-NN)+NNN)\
	PRED(3*(W-WW)+WWW)\
	PRED((WWW+NEEE)>>1)\

enum
{
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED
};


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

static int hist[3][NCTX][256], hist_wp[3][NCTX][256];

#ifndef FLOOR_LOG2
#define FLOOR_LOG2(X)		(sizeof(X)==8?63-(int)_lzcnt_u64((unsigned long long)(X)):31-(int)_lzcnt_u32((unsigned)(X)))
#endif
#ifndef CLAMP2
#define CLAMP2(X, LO, HI)\
	do\
	{\
		if((X)<(LO))X=LO;\
		if((X)>(HI))X=HI;\
	}while(0)
#endif
#ifndef LOG_ERROR
static int log_error(const char *fn, int line, int quit, const char *format, ...)
{
	ptrdiff_t size=strlen(fn), start=size-1;
	for(;start>=0&&fn[start]!='/'&&fn[start]!='\\';--start);
	start+=start==-1||fn[start]=='/'||fn[start]=='\\';

	printf("\n%s(%d): ", fn+start, line);
	va_list args;
	va_start(args, format);
	vprintf(format, args);
	va_end(args);

	if(quit)
		exit(0);
	return 0;
}
#define LOG_ERROR(MSG, ...) log_error(file, __LINE__, 1, MSG, ##__VA_ARGS__)
#endif
static unsigned char* load_ppm(const char *fn, int *ret_iw, int *ret_ih)
{
	FILE *fsrc=fopen(fn, "rb");
	if(!fsrc)
	{
		LOG_ERROR("Cannot open \"%s\"", fn);
		return 0;
	}
	int tag=0;
	fread(&tag, 1, 2, fsrc);
	if(tag!=('P'|'6'<<8))
	{
		LOG_ERROR("Unsupported file \"%s\"", fn);
		return 0;
	}
#ifdef LOUD
	print_timestamp("%Y-%m-%d_%H%M%S\n");
#endif
	int temp=fgetc(fsrc);
	if(temp!='\n')
	{
		LOG_ERROR("Invalid PPM file");
		return 0;
	}
	int iw=0, ih=0;
	int nread=fscanf(fsrc, "%d %d", &iw, &ih);
	if(nread!=2)
	{
		LOG_ERROR("Unsupported PPM file");
		return 0;
	}
	int vmax=0;
	nread=fscanf(fsrc, "%d", &vmax);
	if(nread!=1||vmax!=255)
	{
		LOG_ERROR("Unsupported PPM file");
		return 0;
	}
	temp=fgetc(fsrc);
	if(temp!='\n')
	{
		LOG_ERROR("Invalid PPM file");
		return 0;
	}
	if(iw<1||ih<1)
	{
		LOG_ERROR("Unsupported source file");
		return 0;
	}
	ptrdiff_t size=(ptrdiff_t)3*iw*ih;
	unsigned char *image=(unsigned char*)malloc(size+sizeof(__m256i));
	fread(image, 1, size, fsrc);//read image
	fclose(fsrc);
	if(ret_iw)*ret_iw=iw;
	if(ret_ih)*ret_ih=ih;
	return image;
}
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
static void calc_csizes(const int *stats, double *csizes)
{
	for(int kc=0;kc<3;++kc)//estimate entropy
	{
		double e=0;
		for(int kctx=0;kctx<NCTX;++kctx)
		{
			const int *currhist=stats+256*(NCTX*kc+kctx);
			//int *currhist=hist[kc][kctx];
			int sum=0;
			for(int ks=0;ks<256;++ks)
				sum+=currhist[ks];
			if(!sum)
				continue;
			for(int ks=0;ks<256;++ks)
			{
				int freq=currhist[ks];
				if(freq)
					e-=freq*log2((double)freq/sum);
			}
		}
		csizes[kc]=e/8;
	}
}
int c37_codec(int argc, char **argv)
{
	if(argc!=2)
	{
		printf(
			"Usage:  \"%s\"  image.ppm\n"
			, argv[0]
		);
		return 1;
	}
	const char *fn=argv[1];

	ptrdiff_t size=0;
	int psize=0;
	short *pixels=0;
	int wpsize=0;
	short *wpixels=0;
	int iw=0, ih=0;
	unsigned char *frame=load_ppm(fn, &iw, &ih);
	if(!frame)
	{
		LOG_ERROR("No frame");
		return 1;
	}
	int bestrct=crct_analysis(frame, iw, ih);
	size=3LL*iw*ih;
	psize=(iw+2*XPAD)*(int)sizeof(short[YPAD*3*3]);//YPAD=4 padded rows * 3 channels * {pixel, l1error, wperror}
	pixels=(short*)malloc(psize);
	wpsize=(iw+2*XPAD)*(int)sizeof(short[YPAD*3*NPREDS]);
	wpixels=(short*)malloc(wpsize);
	if(!pixels||!wpixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(pixels, 0, psize);
	memset(wpixels, 0, wpsize);
	ptrdiff_t idx=0;
	int weights[NPREDS]={0};
	int corrections[NPREDS]={0};
	for(int kp=0;kp<NPREDS;++kp)
		weights[kp]=(1<<L1SH)/NPREDS;
	int yidx=rct_combinations[bestrct][II_PERM_Y];
	int uidx=rct_combinations[bestrct][II_PERM_U];
	int vidx=rct_combinations[bestrct][II_PERM_V];
	int umask=rct_combinations[bestrct][II_COEFF_U_SUB_Y];
	int vc0=rct_combinations[bestrct][II_COEFF_V_SUB_Y];
	int vc1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
	memset(hist, 0, sizeof(hist));
	memset(hist_wp, 0, sizeof(hist_wp));
	for(int ky=0;ky<ih;++ky)
	{
		short *rows[YPAD]={0};
		for(int ky2=0;ky2<YPAD;++ky2)
			rows[ky2]=pixels+(ky-ky2+YPAD)%YPAD*3*3*(iw+2*XPAD)+3*3*XPAD;
		short *wrows[YPAD]={0};
		for(int ky2=0;ky2<YPAD;++ky2)
			wrows[ky2]=wpixels+(ky-ky2+YPAD)%YPAD*3*NPREDS*(iw+2*XPAD)+3*NPREDS*XPAD;
		for(int kx=0;kx<iw;++kx, idx+=3)
		{
			int yuv[3]={0};
			yuv[0]=frame[idx+yidx]-128;
			yuv[1]=frame[idx+uidx]-128;
			yuv[2]=frame[idx+vidx]-128;
			int offset=0;
			for(int kc=0;kc<3;++kc)
			{
				int
					NNNN	=rows[0][+0*3*3],
					NNNNE	=rows[0][+1*3*3],
					NNN	=rows[3][+0*3*3],
					NNW	=rows[2][-1*3*3],
					NNWWW	=rows[2][-3*3*3],
					NN	=rows[2][+0*3*3],
					NNE	=rows[2][+1*3*3],
					NWW	=rows[1][-2*3*3],
					NW	=rows[1][-1*3*3],
					N	=rows[1][+0*3*3],
					NE	=rows[1][+1*3*3],
					NEE	=rows[1][+2*3*3],
					NEEE	=rows[1][+3*3*3],
					WWWW	=rows[0][-4*3*3],
					WWW	=rows[0][-3*3*3],
					WW	=rows[0][-2*3*3],
					W	=rows[0][-1*3*3];
				int
					eNEE	=rows[1][+2*3*3+1],
					eNEEE	=rows[1][+3*3*3+1],
					eW	=rows[0][-0*3*3+1];
				int
					fNEE	=rows[1][+2*3*3+2],
					fNEEE	=rows[1][+3*3*3+2],
					fW	=rows[0][-0*3*3+2];
				int preds[]=
				{
#define PRED(...) __VA_ARGS__,
					PREDLIST
#undef  PRED
				};

				//WP predict
				int wperrors[NPREDS]={0};
				int wpweights[NPREDS]={0};
				for(int kp=0;kp<NPREDS;++kp)
				{
					//		1
					//		2
					//		5	1
					//	5	?
					wpweights[kp]=
						+wrows	[3][kp+0*3*NPREDS]	//+eNNN
						+wrows	[2][kp+0*3*NPREDS]	//+eNN	*2
						+wrows	[1][kp+0*3*NPREDS]*3	//+eN	*5
						+wrows	[1][kp+1*3*NPREDS]	//+eNE
						+wrows	[0][kp-1*3*NPREDS]*3	//+eW	*5

					//	+wrows	[2][kp-1*3*NPREDS]	//+eNNW		X
					//	+wrows	[1][kp-1*3*NPREDS]*2	//+eNW	*2	X
					//	+wrows	[0][kp-1*3*NPREDS]*5	//+eW	*5	X
					//	+wg_perrors[kc][kp]		//+I
						+1
					;
				}
				memcpy(wperrors, wpweights, sizeof(wperrors));
				{
					double rsum=0;
					double rweights[NPREDS]={0};
					for(int kp=0;kp<NPREDS;++kp)
						rsum+=rweights[kp]=1./wpweights[kp];
					if(rsum)
					{
						rsum=1./rsum;
						for(int kp=0;kp<NPREDS;++kp)
							rweights[kp]*=rsum;
						for(int kp=0;kp<NPREDS;++kp)
							wpweights[kp]=(int)_mm_cvttsd_si64(_mm_set_sd(rweights[kp]*(1<<WPSH)));
					}
					else
						for(int kp=0;kp<NPREDS;++kp)
							wpweights[kp]=(1<<WPSH)/NPREDS;
				}
				int wppred=1<<WPSH>>1;
				for(int kp=0;kp<NPREDS;++kp)
					wppred+=wpweights[kp]*preds[kp];
				wppred>>=WPSH;

				//L1 predict
#if 0
				for(int kp=0;kp<NPREDS;++kp)
					weights[kp]-=wperrors[kp]>>5;
#endif
				int pred=1<<L1SH>>1;
				for(int kp=0;kp<NPREDS;++kp)
				//	pred+=weights[kp]*(preds[kp]+(corrections[kp]>>(L1SH-5)));
					pred+=weights[kp]*preds[kp];
				pred>>=L1SH;

				if(ky==ih-7&&kx>iw/4&&kx<iw*3/4&&!kc)
				{
					static int ctr=0;
					++ctr;
#define LGWIDTH 6
					int nprinted=0;

					//printf("\n");
					//printf("%4d\n", yuv[kc]-offset);

					//print WP weights
#if 1
					nprinted+=printf("%5d %4d %4d WP ", ctr, yuv[kc]-offset, wppred);
					for(int kp=0, wsum=0;kp<NPREDS;++kp)
					{
						int nsum=wsum+wpweights[kp];
						int nstars=abs((nsum>>(WPSH-LGWIDTH))-(wsum>>(WPSH-LGWIDTH)));
						if(kp)
							nprinted+=printf(" ");
						for(int k2=0;k2<nstars;++k2)
							nprinted+=printf("%c", (wpweights[kp]<0?'\b':'0'+kp));
						wsum=nsum;
					}
					//	printf("%*s|", wpweights[kp]<<LGWIDTH>>WPSH, "");
					//printf("\n");
					nprinted+=printf("%*s", 100-nprinted, "");
#endif

#if 1
					int sum=0;
					for(int kp=0;kp<NPREDS;++kp)
						sum+=wperrors[kp];
					for(int kp=0, wsum=0;kp<NPREDS;++kp)
					{
						int nsum=wsum+wperrors[kp];
						int nstars=sum?abs((nsum<<LGWIDTH)/sum-(wsum<<LGWIDTH)/sum):(1<<LGWIDTH)/NPREDS;
						if(kp)
							nprinted+=printf(" ");
						for(int k2=0;k2<nstars;++k2)
							nprinted+=printf("%c", (wperrors[kp]<0?'\b':'0'+kp));
						wsum=nsum;
					}
					nprinted+=printf("%*s", 200-nprinted, "");

#endif

					//print L1 weights
#if 1
					nprinted+=printf("%5d %4d %4d L1 ", ctr, yuv[kc]-offset, pred);
					for(int kp=0, wsum=0;kp<NPREDS;++kp)
					{
						int nsum=wsum+weights[kp];
						int nstars=(nsum>>(L1SH-LGWIDTH))-(wsum>>(L1SH-LGWIDTH));
						if(kp)
							nprinted+=printf(" ");
						for(int k2=0;k2<nstars;++k2)
							nprinted+=printf("%c", (weights[kp]<0?'\b':'0'+kp));
						wsum=nsum;

						//printf("%*s", abs(weights[kp])<<6>>L1SH, weights[kp]<0?"<":"|");
					}
#endif
					printf("\n");
				}

				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;

				//WP clamp
				int wp0=wppred;
				CLAMP2(wppred, vmin, vmax);
				wppred+=offset;
				CLAMP2(wppred, -128, 127);
				
				//L1 clamp
				int p0=pred;
				CLAMP2(pred, vmin, vmax);
				pred+=offset;
				CLAMP2(pred, -128, 127);

				//L1 context
				int ctx=FLOOR_LOG2(eW*eW+1);
				if(ctx>NCTX-1)
					ctx=NCTX-1;

				//WP context
				int wpctx=FLOOR_LOG2(fW*fW+1);
				if(wpctx>NCTX-1)
					wpctx=NCTX-1;
					
				int curr=yuv[kc];//reveal current

				int delta=curr-pred;
				++hist[kc][ctx][(delta+128)&255];//accumulate stats

				int wpdelta=curr-wppred;
				++hist_wp[kc][wpctx][(wpdelta+128)&255];//accumulate stats
				
				curr-=offset;
				rows[0][0]=curr;

				//WP update
				for(int kp=0;kp<NPREDS;++kp)
					wrows[0][kp]=abs(curr-preds[kp]);

				//L1 update
				int esign=(curr>p0)-(curr<p0);//update weights
				//int esign=curr-p0;
				for(int kp=0;kp<NPREDS;++kp)
					corrections[kp]+=esign*weights[kp]>>(L1SH-4);
				for(int kp=0;kp<NPREDS;++kp)
					weights[kp]+=esign*preds[kp];

				delta=delta<<1^delta>>31;//update L1 context
				rows[0][1]=(2*eW+(delta<<3)+(eNEE>eNEEE?eNEE:eNEEE))>>2;

				wpdelta=wpdelta<<1^wpdelta>>31;//update WP context
				rows[0][2]=(2*fW+(wpdelta<<3)+(fNEE>fNEEE?fNEE:fNEEE))>>2;
					
				offset=(kc ? vc0*yuv[0]+vc1*yuv[1] : umask*yuv[0])>>2;

				for(int ky2=0;ky2<YPAD;++ky2)
					rows[ky2]+=3;
				for(int ky2=0;ky2<YPAD;++ky2)
					wrows[ky2]+=NPREDS;
			}
		}
	}

	double csizes[3]={0};
	double csizes_wp[3]={0};
	calc_csizes((int*)hist, csizes);
	calc_csizes((int*)hist_wp, csizes_wp);
	printf(
		"TYUV %12.2lf %12.2lf %12.2lf %12.2lf  %6.2lf %6.2lf %6.2lf %6.2lf%%  L1  sh=%d\n"
		, csizes[0]+csizes[1]+csizes[2]
		, csizes[0]
		, csizes[1]
		, csizes[2]
		, (csizes[0]+csizes[1]+csizes[2])*100/(3*iw*ih)
		, csizes[0]*100/(iw*ih)
		, csizes[1]*100/(iw*ih)
		, csizes[2]*100/(iw*ih)
		, L1SH
	);
	printf(
		"TYUV %12.2lf %12.2lf %12.2lf %12.2lf  %6.2lf %6.2lf %6.2lf %6.2lf%%  WP\n"
		, csizes_wp[0]+csizes_wp[1]+csizes_wp[2]
		, csizes_wp[0]
		, csizes_wp[1]
		, csizes_wp[2]
		, (csizes_wp[0]+csizes_wp[1]+csizes_wp[2])*100/(3*iw*ih)
		, csizes_wp[0]*100/(iw*ih)
		, csizes_wp[1]*100/(iw*ih)
		, csizes_wp[2]*100/(iw*ih)
	);
	printf("\"%s\"  CWH 3*%d*%d  %d bytes\n", fn, iw, ih, 3*iw*ih);
	free(pixels);
	free(wpixels);
	free(frame);

	exit(0);//this is not a codec
	return 0;
}
