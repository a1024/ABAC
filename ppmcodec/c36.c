#include"util.h"
//#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
//#define _CRT_SECURE_NO_WARNINGS
//#endif
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


	#define ERROR_BOUNDS
	#define LOSSY


#define PREDLIST\
	PRED(p)\
	PRED(3*(p-pp)+ppp)\
	PRED(N+p-pN)\
	PRED(W+p-pW)\
	PRED(N)\
	PRED(W)\
	PRED(NW)\
	PRED(NE)\
	PRED(N+W-NW)\
	PRED(W+NE-N)\
	PRED(3*(N-NN)+NNN)\
	PRED(3*(W-WW)+WWW)\

#if 0
	PRED(2*p-pp)\
	PRED(p+pN-ppN)\
	PRED(p+pW-ppW)\
	PRED(p+pS-ppS)\
	PRED(p+pE-ppE)\

#endif
enum
{
	NCTX=18,
	SH=18,
#define PRED(...) +1
	NPREDS=PREDLIST,
#undef  PRED

	NFRAMES=4,
	XPAD=4,
	YPAD=4,
	NCH=3,
	NROWS=2*YPAD+1,
	NVAL=NFRAMES+1,
#ifdef LOSSY
	DIST=8,
	INVDIST=((1<<16)+DIST-1)/DIST,
#endif
};


//cRCT
#if 1
	#define ENABLE_EXTENDED_RCT
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

#ifdef LOSSY
static void calc_psnr(const int64_t *esum, int64_t count, double *psnr, int print)
{
	double mse[4];

	mse[0]=(double)esum[0]/count;
	mse[1]=(double)esum[1]/count;
	mse[2]=(double)esum[2]/count;
	mse[3]=(mse[0]+mse[1]+mse[2])/3;
	psnr[0]=mse[0]?-10*log10(mse[0]*(1./(255*255))):INFINITY;
	psnr[1]=mse[1]?-10*log10(mse[1]*(1./(255*255))):INFINITY;
	psnr[2]=mse[2]?-10*log10(mse[2]*(1./(255*255))):INFINITY;
	psnr[3]=mse[3]?-10*log10(mse[3]*(1./(255*255))):INFINITY;
	if(print)
		printf(
			"\n"
			"RMSE %12.6lf %12.6lf %12.6lf %12.6lf\n"
			"PSNR %12.6lf %12.6lf %12.6lf %12.6lf\n"
			"\n"
			, sqrt(mse[3]), sqrt(mse[0]), sqrt(mse[1]), sqrt(mse[2])
			, psnr[3], psnr[0], psnr[1], psnr[2]
		);
}
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
static int hist[3][NCTX][256];
int c36_codec(int argc, char **argv)
{
	if(argc!=2)
	{
		printf(
			"Usage:  \"%s\"  folder\n"
			"  folder must contain only PPM frames of same resolution\n"
			, argv[0]
		);
		return 1;
	}
	const char *ext[]=
	{
		"PPM",
	};
	const char *path=argv[1];

	ptrdiff_t size=0, offset=0;
	int iw0=0, ih0=0;
	unsigned char *frames[NFRAMES]={0};
	ArrayHandle filenames=get_filenames(path, (const char**)ext, _countof(ext), 1);
	int psize=0;
	int16_t *pixels=0;
	//int cpsize=0;
	//short *cpixels=0;
#ifdef LOSSY
	int64_t etotal[3]={0};
#endif
	double total_csize[3]={0};
	double t=time_sec();
#ifdef ERROR_BOUNDS
	int emin[3]={0}, emax[3]={0};
	int64_t esum0=0, ecount0=0;
	int64_t ehist[8]={0};
#endif

	if(!filenames)
	{
		LOG_ERROR("No frames in \"%s\"", path);
		return 1;
	}
	printf("Frame  TYUV (bytes)  TYUV (%%)  %d frames\n", (int)filenames->count);
	for(int k0=0;k0<(int)filenames->count;++k0)
	{
#ifdef LOSSY
		int64_t esum[3]={0};
#endif
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, k0);
		int iw=0, ih=0;
		unsigned char *curr=load_ppm((char*)fn[0]->data, &iw, &ih);
		if(!curr)
		{
			LOG_ERROR("No frame");
			return 1;
		}
		int bestrct=crct_analysis(curr, iw, ih);
		if(!frames[0])
		{
			iw0=iw;
			ih0=ih;
			size=3LL*iw*ih;
			offset=3LL*YPAD*iw;
			for(int kf=0;kf<NFRAMES;++kf)
			{
				frames[kf]=(unsigned char*)malloc(size+offset*2);
				if(!frames[kf])
				{
					LOG_ERROR("Alloc error");
					return 1;
				}
				memset(frames[kf], 0, size+offset*2);
			}
			psize=(iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);
			//psize=(int)sizeof(short[NFRAMES*(2*YPAD+1)*3])*(iw+2*XPAD);
			pixels=(int16_t*)malloc(psize);
			//cpsize=(int)sizeof(short[(2*YPAD+1)*3])*(iw+2*XPAD);
			//cpixels=(short*)malloc(cpsize);
			if(!pixels)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
		}
		else if(iw!=iw0||ih!=ih0)
		{
			LOG_ERROR("Frame resolution changed");
			return 1;
		}
		memcpy(frames[0]+offset, curr, size);
		free(curr);
		memset(pixels, 0, psize);
		//memset(cpixels, 0, cpsize);
		ptrdiff_t idx=0;
		unsigned weights[NPREDS]={0};
		for(int kp=0;kp<NPREDS;++kp)
			weights[kp]=(1<<SH)/NPREDS;
		int yidx=rct_combinations[bestrct][II_PERM_Y];
		int uidx=rct_combinations[bestrct][II_PERM_U];
		int vidx=rct_combinations[bestrct][II_PERM_V];
		int uc0=rct_combinations[bestrct][II_COEFF_U_SUB_Y];
		int vc0=rct_combinations[bestrct][II_COEFF_V_SUB_Y];
		int vc1=rct_combinations[bestrct][II_COEFF_V_SUB_U];
		memset(hist, 0, sizeof(hist));
		for(int ky=0;ky<ih;++ky)
		{
			int16_t *rows[NROWS]={0};
			for(int k=0;k<NROWS;++k)
				rows[k]=pixels+(XPAD*NCH*NROWS+(ky-k+NROWS)%NROWS)*NVAL;
			//short *rows[NFRAMES][2*YPAD+1]={0};
			//for(int kf=0;kf<NFRAMES;++kf)
			//{
			//	for(int ky2=0;ky2<2*YPAD+1;++ky2)
			//		rows[kf][ky2]=pixels+3*(iw+2*XPAD)*((2*YPAD+1)*kf+(ky-ky2+2*YPAD+1)%(2*YPAD+1));
			//}
			//short *crows[2*YPAD+1]={0};
			//for(int ky2=0;ky2<2*YPAD+1;++ky2)
			//	crows[ky2]=cpixels+(ky-ky2+2*YPAD+1)%(2*YPAD+1)*3*(iw+2*XPAD);
			for(int kx=0;kx<iw;++kx)
			{
				int yuv[NFRAMES][3]={0};
				yuv[0][0]=frames[0][idx+yidx];
				yuv[0][1]=frames[0][idx+uidx];
				yuv[0][2]=frames[0][idx+vidx];
				for(int kf=1;kf<NFRAMES;++kf)
				{
					yuv[kf][0]=frames[kf][idx+yidx+3LL*iw];
					yuv[kf][1]=frames[kf][idx+uidx+3LL*iw];
					yuv[kf][2]=frames[kf][idx+vidx+3LL*iw];
				}
				int crct[NFRAMES]={0};
				for(int kc=0;kc<3;++kc, ++idx)
				{
					int
						pppNNN	=rows[YPAD+3][3+0*NCH*NROWS*NVAL],
						pppNN	=rows[YPAD+2][3+0*NCH*NROWS*NVAL],
						pppNW	=rows[YPAD+1][3-1*NCH*NROWS*NVAL],
						pppN	=rows[YPAD+1][3+0*NCH*NROWS*NVAL],
						pppNE	=rows[YPAD+1][3+1*NCH*NROWS*NVAL],
						pppNEE	=rows[YPAD+1][3+2*NCH*NROWS*NVAL],
						pppNEEE	=rows[YPAD+1][3+3*NCH*NROWS*NVAL],
						pppWWW	=rows[YPAD+0][3-3*NCH*NROWS*NVAL],
						pppWW	=rows[YPAD+0][3-2*NCH*NROWS*NVAL],
						pppW	=rows[YPAD+0][3-1*NCH*NROWS*NVAL],
						ppp	=rows[YPAD+0][3+0*NCH*NROWS*NVAL],
						pppE	=rows[YPAD+0][3+1*NCH*NROWS*NVAL],
						pppEE	=rows[YPAD+0][3+2*NCH*NROWS*NVAL],
						pppSW	=rows[YPAD-1][3-1*NCH*NROWS*NVAL],
						pppS	=rows[YPAD-1][3+0*NCH*NROWS*NVAL],
						pppSE	=rows[YPAD-1][3+1*NCH*NROWS*NVAL],
						pppSS	=rows[YPAD-2][3+0*NCH*NROWS*NVAL];
					int
						ppNNN	=rows[YPAD+3][2+0*NCH*NROWS*NVAL],
						ppNN	=rows[YPAD+2][2+0*NCH*NROWS*NVAL],
						ppNW	=rows[YPAD+1][2-1*NCH*NROWS*NVAL],
						ppN	=rows[YPAD+1][2+0*NCH*NROWS*NVAL],
						ppNE	=rows[YPAD+1][2+1*NCH*NROWS*NVAL],
						ppNEE	=rows[YPAD+1][2+2*NCH*NROWS*NVAL],
						ppNEEE	=rows[YPAD+1][2+3*NCH*NROWS*NVAL],
						ppWWW	=rows[YPAD+0][2-3*NCH*NROWS*NVAL],
						ppWW	=rows[YPAD+0][2-2*NCH*NROWS*NVAL],
						ppW	=rows[YPAD+0][2-1*NCH*NROWS*NVAL],
						pp	=rows[YPAD+0][2+0*NCH*NROWS*NVAL],
						ppE	=rows[YPAD+0][2+1*NCH*NROWS*NVAL],
						ppEE	=rows[YPAD+0][2+2*NCH*NROWS*NVAL],
						ppSW	=rows[YPAD-1][2-1*NCH*NROWS*NVAL],
						ppS	=rows[YPAD-1][2+0*NCH*NROWS*NVAL],
						ppSE	=rows[YPAD-1][2+1*NCH*NROWS*NVAL],
						ppSS	=rows[YPAD-2][2+0*NCH*NROWS*NVAL];
					int
						pNNN	=rows[YPAD+3][1+0*NCH*NROWS*NVAL],
						pNN	=rows[YPAD+2][1+0*NCH*NROWS*NVAL],
						pNW	=rows[YPAD+1][1-1*NCH*NROWS*NVAL],
						pN	=rows[YPAD+1][1+0*NCH*NROWS*NVAL],
						pNE	=rows[YPAD+1][1+1*NCH*NROWS*NVAL],
						pNEE	=rows[YPAD+1][1+2*NCH*NROWS*NVAL],
						pNEEE	=rows[YPAD+1][1+3*NCH*NROWS*NVAL],
						pWWW	=rows[YPAD+0][1-3*NCH*NROWS*NVAL],
						pWW	=rows[YPAD+0][1-2*NCH*NROWS*NVAL],
						pW	=rows[YPAD+0][1-1*NCH*NROWS*NVAL],
						p	=rows[YPAD+0][1+0*NCH*NROWS*NVAL],
						pE	=rows[YPAD+0][1+1*NCH*NROWS*NVAL],
						pEE	=rows[YPAD+0][1+2*NCH*NROWS*NVAL],
						pSW	=rows[YPAD-1][1-1*NCH*NROWS*NVAL],
						pS	=rows[YPAD-1][1+0*NCH*NROWS*NVAL],
						pSE	=rows[YPAD-1][1+1*NCH*NROWS*NVAL],
						pSS	=rows[YPAD-2][1+0*NCH*NROWS*NVAL];
					int
						NNN	=rows[YPAD+3][0+0*NCH*NROWS*NVAL],
						NN	=rows[YPAD+2][0+0*NCH*NROWS*NVAL],
						NW	=rows[YPAD+1][0-1*NCH*NROWS*NVAL],
						N	=rows[YPAD+1][0+0*NCH*NROWS*NVAL],
						NE	=rows[YPAD+1][0+1*NCH*NROWS*NVAL],
						NEE	=rows[YPAD+1][0+2*NCH*NROWS*NVAL],
						NEEE	=rows[YPAD+1][0+3*NCH*NROWS*NVAL],
						WWW	=rows[YPAD+0][0-3*NCH*NROWS*NVAL],
						WW	=rows[YPAD+0][0-2*NCH*NROWS*NVAL],
						W	=rows[YPAD+0][0-1*NCH*NROWS*NVAL];
					int
						eNEE	=rows[YPAD+1][NFRAMES+2*NCH*NROWS*NVAL],
						eNEEE	=rows[YPAD+1][NFRAMES+3*NCH*NROWS*NVAL],
						eW	=rows[YPAD+0][NFRAMES-1*NCH*NROWS*NVAL];
#if 0
					int
						pppNN	=rows[3][YPAD+2][+0*3],
						pppNW	=rows[3][YPAD+1][-1*3],
						pppN	=rows[3][YPAD+1][+0*3],
						pppNE	=rows[3][YPAD+1][+1*3],
						pppWW	=rows[3][YPAD+0][-2*3],
						pppW	=rows[3][YPAD+0][-1*3],
						ppp	=rows[3][YPAD+0][+0*3],
						pppE	=rows[3][YPAD+0][+1*3],
						pppEE	=rows[3][YPAD+0][+2*3],
						pppSW	=rows[3][YPAD-1][-1*3],
						pppS	=rows[3][YPAD-1][+0*3],
						pppSE	=rows[3][YPAD-1][+1*3],
						pppSS	=rows[3][YPAD-2][+0*3];
					int
						ppNN	=rows[2][YPAD+2][+0*3],
						ppNW	=rows[2][YPAD+1][-1*3],
						ppN	=rows[2][YPAD+1][+0*3],
						ppNE	=rows[2][YPAD+1][+1*3],
						ppWW	=rows[2][YPAD+0][-2*3],
						ppW	=rows[2][YPAD+0][-1*3],
						pp	=rows[2][YPAD+0][+0*3],
						ppE	=rows[2][YPAD+0][+1*3],
						ppEE	=rows[2][YPAD+0][+2*3],
						ppSW	=rows[2][YPAD-1][-1*3],
						ppS	=rows[2][YPAD-1][+0*3],
						ppSE	=rows[2][YPAD-1][+1*3],
						ppSS	=rows[2][YPAD-2][+0*3];
					int
						pNN	=rows[1][YPAD+2][+0*3],
						pNW	=rows[1][YPAD+1][-1*3],
						pN	=rows[1][YPAD+1][+0*3],
						pNE	=rows[1][YPAD+1][+1*3],
						pWW	=rows[1][YPAD+0][-2*3],
						pW	=rows[1][YPAD+0][-1*3],
						p	=rows[1][YPAD+0][+0*3],
						pE	=rows[1][YPAD+0][+1*3],
						pEE	=rows[1][YPAD+0][+2*3],
						pSW	=rows[1][YPAD-1][-1*3],
						pS	=rows[1][YPAD-1][+0*3],
						pSE	=rows[1][YPAD-1][+1*3],
						pSS	=rows[1][YPAD-2][+0*3];
					int
						NNN	=rows[0][YPAD+3][+0*3],
						NN	=rows[0][YPAD+2][+0*3],
						NW	=rows[0][YPAD+1][-1*3],
						N	=rows[0][YPAD+1][+0*3],
						NE	=rows[0][YPAD+1][+1*3],
						NEE	=rows[0][YPAD+1][+2*3],
						NEEE	=rows[0][YPAD+1][+3*3],
						WWW	=rows[0][YPAD+0][-3*3],
						WW	=rows[0][YPAD+0][-2*3],
						W	=rows[0][YPAD+0][-1*3];
					int
						eNEE	=crows[YPAD+1][+2*3+1],
						eNEEE	=crows[YPAD+1][+3*3+1],
						eW	=crows[YPAD+0][+0*3+1];
#endif
					int preds[]=
					{
#define PRED(...) __VA_ARGS__,
						PREDLIST
#undef  PRED
					};
					int pred=1<<SH>>1;
					for(int kp=0;kp<NPREDS;++kp)
						pred+=weights[kp]*preds[kp];
					pred>>=SH;
					int p0=pred;
					int vmax=N, vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(pred, vmin, vmax);
					pred+=crct[0];
					CLAMP2(pred, 0, 255);

					int ctx=FLOOR_LOG2(eW*eW+1);
					if(ctx>NCTX-1)
						ctx=NCTX-1;
					
					int curr=yuv[0][kc];
					int delta=curr-pred;
					//if(delta)
					//	printf("");
#ifdef LOSSY
					//delta=(int8_t)delta;
					delta=(delta*INVDIST>>16)-(DIST>1?delta>>31:0);
#endif

					++hist[kc][ctx][(delta+128)&255];//accumulate stats
#ifdef LOSSY
					int c0=curr;
					curr=delta*DIST+pred;
					CLAMP2(curr, 0, 255);
					//if(curr!=c0)
					//	LOG_ERROR("");
					c0-=curr;
					esum[kc]+=(int64_t)c0*c0;
#ifdef ERROR_BOUNDS
					esum0+=abs(c0); ++ecount0; ++ehist[abs(c0)];
					if(emin[kc]>c0)emin[kc]=c0;
					if(emax[kc]<c0)emax[kc]=c0;
#endif
#endif

					rows[YPAD+0][0+0*NCH*NROWS*NVAL]=yuv[0][kc]-crct[0];
					for(int kf=1;kf<NFRAMES;++kf)
						rows[YPAD+YPAD][0+0*NCH*NROWS*NVAL]=yuv[kf][kc]-crct[kf];
					//for(int kf=0;kf<NFRAMES;++kf)//update circular buffers
					//	rows[kf][YPAD+YPAD*(kf!=0)][0]=yuv[kf][kc]-crct[kf];

					curr-=crct[0];
					int esign=(curr>p0)-(curr<p0);//update weights
					for(int kp=0;kp<NPREDS;++kp)
						weights[kp]+=esign*preds[kp];

					delta=delta<<1^delta>>31;//update context
					rows[YPAD+0][NFRAMES+0*NCH*NROWS*NVAL]=(2*eW+(delta<<3)+MAXVAR(eNEE, eNEEE))>>2;
					//crows[YPAD+0][0]=(2*eW+(delta<<3)+MAXVAR(eNEE, eNEEE))>>2;
					
					for(int kf=0;kf<NFRAMES;++kf)//update crct
						crct[kf]=(kc ? vc0*yuv[kf][0]+vc1*yuv[kf][1] : uc0*yuv[kf][0])>>2;

					for(int k=0;k<NROWS;++k)
						rows[k]+=NROWS*NVAL;
					//for(int kf=0;kf<NFRAMES;++kf)//increment pointers
					//{
					//	for(int ky2=0;ky2<2*YPAD+1;++ky2)
					//		rows[kf][ky2]+=2;
					//}
					//for(int ky2=0;ky2<2*YPAD+1;++ky2)
					//	++crows[ky2];
				}
			}
		}

		double csize[3]={0};
		for(int kc=0;kc<3;++kc)//estimate entropy
		{
			double e=0;
			for(int kctx=0;kctx<NCTX;++kctx)
			{
				int *currhist=hist[kc][kctx];
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
			total_csize[kc]+=csize[kc]=e/8;
		}
		printf(
			"%5d    %12.2lf %12.2lf %12.2lf %12.2lf    %6.2lf %6.2lf %6.2lf %6.2lf%%"
			, k0
			, csize[0]+csize[1]+csize[2]
			, csize[0]
			, csize[1]
			, csize[2]
			, (csize[0]+csize[1]+csize[2])*100/(3*iw*ih)
			, csize[0]*100/(iw*ih)
			, csize[1]*100/(iw*ih)
			, csize[2]*100/(iw*ih)
		);
#ifdef LOSSY
		double psnr[4];
		calc_psnr(esum, (int64_t)iw*ih, psnr, 0);
		printf("  %8.4lf %8.4lf %8.4lf %8.4lf"
			, psnr[0]
			, psnr[1]
			, psnr[2]
			, psnr[3]
		);
		etotal[0]+=esum[0];
		etotal[1]+=esum[1];
		etotal[2]+=esum[2];
#endif
		printf("\n");

		unsigned char *temp=frames[0];
		for(int kf=1;kf<NFRAMES;++kf)
			frames[kf-1]=frames[kf];
		frames[NFRAMES-1]=temp;
	}
	t=time_sec()-t;
	printf(
		"\n"
		"%5d    %12.2lf %12.2lf %12.2lf %12.2lf    %6.2lf %6.2lf %6.2lf %6.2lf%%    %12.2lf"
		, (int)filenames->count
		, total_csize[0]+total_csize[1]+total_csize[2]
		, total_csize[0]
		, total_csize[1]
		, total_csize[2]
		, (total_csize[0]+total_csize[1]+total_csize[2])*100/(filenames->count*3*iw0*ih0)
		, total_csize[0]*100/(filenames->count*iw0*ih0)
		, total_csize[1]*100/(filenames->count*iw0*ih0)
		, total_csize[2]*100/(filenames->count*iw0*ih0)
		, (double)filenames->count*3*iw0*ih0
	);
#ifdef LOSSY
	double psnr[4];
	calc_psnr(etotal, (int64_t)filenames->count*iw0*ih0, psnr, 1);
	printf("  %8.4lf %8.4lf %8.4lf %8.4lf"
		, psnr[0]
		, psnr[1]
		, psnr[2]
		, psnr[3]
	);
#endif
	printf("\n");
#ifdef ERROR_BOUNDS
	printf("%d ~ %d\n", emin[0], emax[0]);
	printf("%d ~ %d\n", emin[1], emax[1]);
	printf("%d ~ %d\n", emin[2], emax[2]);
	printf("%12.6lf\n", (double)esum0/ecount0);
	printf("0  %16lld\n", ehist[0]);
	printf("1  %16lld\n", ehist[1]);
	printf("2  %16lld\n", ehist[2]);
	printf("3  %16lld\n", ehist[3]);
	printf("4  %16lld\n", ehist[4]);
	printf("5  %16lld\n", ehist[5]);
	printf("6  %16lld\n", ehist[6]);
	printf("7  %16lld\n", ehist[7]);
#endif
	double usize=(double)filenames->count*3*iw0*ih0;
	printf("%12.6lf sec  %12.6lf MB/s  %12.6lf ms/MB\n"
		, t
		, usize/(t*1024*1024)
		, t*1024*1024*1000/usize
	);
	for(int k=0;k<NFRAMES;++k)
		free(frames[k]);
	array_free(&filenames);
	free(pixels);
	//free(cpixels);

	exit(0);//this is not a codec
	return 0;
}
