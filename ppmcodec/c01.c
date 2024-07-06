#include"codec.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT


//select one:
//	#define ENABLE_MIX4//good
//	#define ENABLE_MIX8//slightly worse
	#define ENABLE_ABAC2
//	#define ENABLE_ABAC//bad
//	#define ENABLE_CALICCTX//bad
//	#define ENABLE_HASH//bad

	#define AC3_PREC

#define CODECNAME "C01"
#define AC_IMPLEMENTATION
#include"entropy.h"

#define BLOCKSIZE 384
#define MAXPRINTEDBLOCKS 200
#ifdef ENABLE_CALICCTX
#define CLEVELS (3*3*3*3*3*3*3*3)
#elif defined ENABLE_HASH
#define CLEVELS 0x1000
#define HASH_CANTOR(A, B) (((A)+(B))*((A)+(B)+1)/2+(B))
#elif defined ENABLE_MIX8
#define MIXBITSX 5
#define MIXBITSY 3
#define MIXBITSZ 3
#define CXLEVELS ((1<<(8-MIXBITSX))+1)
#define CYLEVELS (8+1)
#define CZLEVELS (8+1)
#elif defined ENABLE_ABAC2
#define A2_CTXBITS 7
#define A2_NCTX 8
#elif defined ENABLE_ABAC
#define ABAC_TLEVELS 256
#define ABAC_TOKEN_BITS 8
//#define ABAC_TLEVELS 25
//#define ABAC_TOKEN_BITS 5
#define ABAC_TREESIZE (1<<ABAC_TOKEN_BITS)	//25/32 is used
#define ABAC_PCTX 23//1//18
#define ABAC_ECTX 13//2//16
#define ABAC_NCTX (ABAC_PCTX+ABAC_ECTX)
#define ABAC_NCTRS 5
#define ABAC_CLEVELS 32
#define ABAC_PROBBITS 16	//minimum 12
#define ABAC_MIXBITS 7
	#define ABAC_SIMPLEOVERRIDE
//	#define DISABLE_LOGMIX

//	#define ABAC_PROFILESIZE
#else
#define CLEVELS 9
#define MIXBITS 8
#endif

#define OCHLIST\
	OCH(R)\
	OCH(G)\
	OCH(B)\
	OCH(RG)\
	OCH(RB)\
	OCH(GB)\
	OCH(GR)\
	OCH(BR)\
	OCH(BG)\
	OCH(R2)\
	OCH(G2)\
	OCH(B2)
typedef enum _OCHIndex
{
#define OCH(LABEL) OCH_##LABEL,
	OCHLIST
#undef  OCH
	OCH_COUNT,
} OCHIndex;
static const char *och_names[OCH_COUNT]=
{
#define OCH(LABEL) #LABEL,
	OCHLIST
#undef  OCH
};

#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B,		3, 3, 3,	3, 3, 3,	0, 0, 0)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG,		3, 3, 1,	3, 3, 3,	0, 0, 0)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR,		3, 3, 0,	3, 3, 3,	0, 0, 0)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG,		3, 3, 0,	3, 3, 3,	0, 0, 0)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB,		3, 3, 1,	3, 3, 3,	0, 0, 0)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR,		3, 3, 1,	3, 3, 3,	0, 0, 0)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB,		3, 3, 0,	3, 3, 3,	0, 0, 0)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG,		3, 0, 0,	3, 3, 3,	0, 0, 0)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB,		3, 0, 0,	3, 3, 3,	0, 0, 0)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(R_GR_BR,	OCH_R,		OCH_GR,		OCH_BR,		3, 0, 0,	3, 3, 3,	0, 0, 0)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB,		3, 0, 1,	3, 3, 3,	0, 0, 0)\
	RCT(R_G_B2,	OCH_R,		OCH_G,		OCH_B2,		3, 3, 0,	3, 3, 1,	0, 0, 1)\
	RCT(R_B_G2,	OCH_R,		OCH_B,		OCH_G2,		3, 3, 0,	3, 3, 1,	0, 0, 1)\
	RCT(G_B_R2,	OCH_G,		OCH_B,		OCH_R2,		3, 3, 0,	3, 3, 1,	0, 0, 1)\
	RCT(R_GR_B2,	OCH_R,		OCH_GR,		OCH_B2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(R_BR_G2,	OCH_R,		OCH_BR,		OCH_G2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(G_BG_R2,	OCH_G,		OCH_BG,		OCH_R2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(G_RG_B2,	OCH_G,		OCH_RG,		OCH_B2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(B_RB_G2,	OCH_B,		OCH_RB,		OCH_G2,		3, 0, 0,	3, 3, 1,	0, 0, 1)\
	RCT(B_GB_R2,	OCH_B,		OCH_GB,		OCH_R2,		3, 0, 0,	3, 3, 1,	0, 0, 1)
typedef enum _RCTIndex
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2, YSH2, USH2, VSH2) RCT_##LABEL,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTIndex;
static const unsigned char rct_combinations[RCT_COUNT][12]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2, YSH2, USH2, VSH2)\
	{YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2, YSH2, USH2, VSH2},
	RCTLIST
#undef  RCT
};
static const char *rct_names[RCT_COUNT]=
{
#define RCT(LABEL, YIDX, UIDX, VIDX, YOFF1, UOFF1, VOFF1, YOFF2, UOFF2, VOFF2, YSH2, USH2, VSH2) #LABEL,
	RCTLIST
#undef  RCT
};

#define PREDLIST\
	PRED(W)\
	PRED(CG)\
	PRED(AV5)\
	PRED(AV9)\
	PRED(AV12)
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
typedef enum _NBIndex
{
	NB_NNWW,	NB_NNW,		NB_NN,		NB_NNE,		NB_NNEE,
	NB_NWW,		NB_NW,		NB_N,		NB_NE,		NB_NEE,
	NB_WW,		NB_W,		NB_curr,
} NBIndex;
static const short av12_icoeffs[12]=
{
	 0x04,	 0x03,	-0x1F,	-0x26,	 0x00,
	 0x07,	-0x9E,	 0xDB,	 0x1E,	 0x13,
	-0x2A,	 0xF3,
};

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
#define CONFIG_EXP 4
#define CONFIG_MSB 1
#define CONFIG_LSB 0
static void quantize_pixel(int val, int *token, int *bypass, int *nbits)
{
	if(val<(1<<CONFIG_EXP))
	{
		*token=val;//token
		*nbits=0;
		*bypass=0;
	}
	else
	{
		int lgv=FLOOR_LOG2((unsigned)val);
		int mantissa=val-(1<<lgv);
		*token = (1<<CONFIG_EXP) + (
				(lgv-CONFIG_EXP)<<(CONFIG_MSB+CONFIG_LSB)|
				(mantissa>>(lgv-CONFIG_MSB))<<CONFIG_LSB|
				(mantissa&((1<<CONFIG_LSB)-1))
			);
		*nbits=lgv-(CONFIG_MSB+CONFIG_LSB);
		*bypass=val>>CONFIG_LSB&((1LL<<*nbits)-1);
	}
}
#ifdef ENABLE_MIX4
static int f28_mix4(int v00, int v01, int v10, int v11, int alphax, int alphay)
{
	//v00=v00*((1<<12)-alphax)+v01*alphax;
	v00+=(v01-v00)*alphax>>MIXBITS;
	v10+=(v11-v10)*alphax>>MIXBITS;
	v00+=(v10-v00)*alphay>>MIXBITS;
	//v00=((v00<<MIXBITS)+(v01-v00)*alphax+(1<<(MIXBITS-1)>>1))>>MIXBITS;//X
	//v10=((v10<<MIXBITS)+(v11-v10)*alphax+(1<<(MIXBITS-1)>>1))>>MIXBITS;
	//v00=((v00<<MIXBITS)+(v10-v00)*alphay+(1<<(MIXBITS-1)>>1))>>MIXBITS;
	return v00;
}
#elif defined ENABLE_MIX8
static int f28_mix8(int v000, int v001, int v010, int v011, int v100, int v101, int v110, int v111, int alphax, int alphay, int alphaz)
{
	v000=((v000<<MIXBITSX)+(v001-v000)*alphax)>>MIXBITSX;
	v010=((v010<<MIXBITSX)+(v011-v010)*alphax)>>MIXBITSX;
	v100=((v100<<MIXBITSX)+(v101-v100)*alphax)>>MIXBITSX;
	v110=((v110<<MIXBITSX)+(v111-v110)*alphax)>>MIXBITSX;
	v000=((v000<<MIXBITSY)+(v010-v000)*alphay)>>MIXBITSY;
	v100=((v100<<MIXBITSY)+(v110-v100)*alphay)>>MIXBITSY;
	v000=((v000<<MIXBITSZ)+(v100-v000)*alphaz)>>MIXBITSZ;
	//v000=((v000<<MIXBITSX)+(v001-v000)*alphax+(1<<MIXBITSX>>1))>>MIXBITSX;//X
	//v010=((v010<<MIXBITSX)+(v011-v010)*alphax+(1<<MIXBITSX>>1))>>MIXBITSX;
	//v100=((v100<<MIXBITSX)+(v101-v100)*alphax+(1<<MIXBITSX>>1))>>MIXBITSX;
	//v110=((v110<<MIXBITSX)+(v111-v110)*alphax+(1<<MIXBITSX>>1))>>MIXBITSX;
	//v000=((v000<<MIXBITSY)+(v010-v000)*alphay+(1<<MIXBITSY>>1))>>MIXBITSY;
	//v100=((v100<<MIXBITSY)+(v110-v100)*alphay+(1<<MIXBITSY>>1))>>MIXBITSY;
	//v000=((v000<<MIXBITSZ)+(v100-v000)*alphaz+(1<<MIXBITSZ>>1))>>MIXBITSZ;
	return v000;
}
#elif defined ENABLE_ABAC
static int squash(int x)//sigmoid(x) = 1/(1-exp(-x))		logit sum -> prob
{
#ifdef DISABLE_LOGMIX
	x>>=11;
	x+=1<<ABAC_PROBBITS>>1;
	CLAMP2_32(x, x, 1, (1<<ABAC_PROBBITS)-1);
#else
	static const int t[33]=//2^5 table elements, table amplitude 2^12
	{
		   1,    2,    3,    6,   10,   16,   27,   45,   73,  120,  194,
		 310,  488,  747, 1101, 1546, 2047, 2549, 2994, 3348, 3607, 3785,
		3901, 3975, 4022, 4050, 4068, 4079, 4085, 4089, 4092, 4093, 4094,
	};
	int w=x&((1<<(ABAC_PROBBITS-5))-1);
	x=(x>>(ABAC_PROBBITS-5))+16;
	if(x>31)
		return (1<<ABAC_PROBBITS)-1;
	if(x<0)
		return 1;
	x=(t[x]*((1<<(ABAC_PROBBITS-5))-w)+t[x+1]*w+64)>>(12-5);
#endif
	return x;
}
static int stretch(int x)//ln(x/(1-x))		probs -> logits
{
#ifndef DISABLE_LOGMIX
	static short t[4096];
	static int initialized=0;
	if(!initialized)
	{
		initialized=1;
		
		int pi=0;
		for(int k=-2047;k<=2047;++k)//invert squash()
		{
			int i=squash(k<<(ABAC_PROBBITS-12))>>(ABAC_PROBBITS-12);
			for(int j=pi;j<=i;++j)
				t[j]=k<<(ABAC_PROBBITS-12);
			pi=i+1;
		}
		t[4095]=(1<<ABAC_PROBBITS>>1)-1;
	}
	x=t[x>>(ABAC_PROBBITS-12)];
#endif
	return x;
}
//static int g_compare=0;
//static int
//	g_probs1[180]={0},
//	g_probs2[180]={0};
static unsigned abac_predict(const unsigned short **stats, int tidx, const int **mixer, int *alphas, int *weights, int *logits)
{
	//const int prob_sh=22+16-ABAC_PROBBITS;//17 bit

	const int prob_sh=11+16-ABAC_PROBBITS;//16 bit
	const int coeff=96;//192

	long long p0=0;
	for(int k=0;k<ABAC_NCTX;++k)
	{
		for(int k2=0;k2<ABAC_NCTRS;++k2)
		{
			int k3=ABAC_NCTRS*k+k2;
			int w0=mixer[0][k3]+((mixer[1][k3]-mixer[0][k3])*alphas[0]>>(ABAC_MIXBITS+2));
			int w1=mixer[2][k3]+((mixer[3][k3]-mixer[2][k3])*alphas[0]>>(ABAC_MIXBITS+2));
			weights[k3]=w0+((w1-w0)*alphas[1]>>(ABAC_MIXBITS+2));
			//int w0=mixer[0][k3]+((mixer[1][k3]-mixer[0][k3])*coeff>>(ABAC_MIXBITS+2));
			//int w1=mixer[2][k3]+((mixer[3][k3]-mixer[2][k3])*coeff>>(ABAC_MIXBITS+2));
			//weights[k3]=w0+((w1-w0)*coeff>>(ABAC_MIXBITS+2));
			logits[k3]=stretch(stats[k][ABAC_NCTRS*tidx+k2]);
			p0+=(long long)weights[k3]*logits[k3];

			//if(g_compare==1)//
			//	g_probs1[k3]=logits[k3], g_probs2[k2]=weights[k3];
			//else if(g_compare==2)//
			//	g_probs1[k3]-=logits[k3], g_probs2[k2]-=weights[k3];
		}
	}
	p0+=1LL<<prob_sh>>1;
	p0>>=prob_sh;
	p0/=ABAC_NCTX*ABAC_NCTRS;
	//p0=((((p0<<ABAC_PROBBITS)+0x8000)>>16)+ABAC_NCTX/2)/ABAC_NCTX;
	//p0+=1LL<<ABAC_PROBBITS>>1;
	//CLAMP2_32(p0, (int)p0, 1, (1<<ABAC_PROBBITS)-1);
	p0=squash((int)p0);
	return (unsigned)p0;
}
static void abac_update(const int *weights, const int *logits, unsigned short **stats, int tidx, int **mixer, int *alphas, unsigned p0final, int bit)
{
	//const int pupdate_sh=31;//17 bit
	//const int mupdate_sh=2*ABAC_MIXBITS+10;

	const int pupdate_sh=27;//16 bit
	const int mupdate_sh=2*ABAC_MIXBITS+10;
	int pupdate_offset=1LL<<pupdate_sh;

	int offset=16<<mupdate_sh;
	int err=(!bit<<ABAC_PROBBITS)-p0final;
	for(int k=0;k<ABAC_NCTX;++k)
	{
		for(int k2=0;k2<ABAC_NCTRS;++k2)
		{
			int k3=ABAC_NCTRS*k+k2;
			int p0=stats[k][ABAC_NCTRS*tidx+k2];
			int m0=mixer[0][k3];
			int m1=mixer[1][k3];
			int m2=mixer[2][k3];
			int m3=mixer[3][k3];
			int mk=weights[k3];
			long long update=(long long)err*logits[k3];
			p0+=(int)(((long long)err*mk+((long long)pupdate_offset<<k2))<<k2>>pupdate_sh);
			m0+=(int)((update*((1<<ABAC_MIXBITS)-alphas[0])*((1<<ABAC_MIXBITS)-alphas[1])-offset)>>mupdate_sh);
			m1+=(int)((update*((1<<ABAC_MIXBITS)-alphas[0])*(                  alphas[1])-offset)>>mupdate_sh);
			m2+=(int)((update*(                  alphas[0])*((1<<ABAC_MIXBITS)-alphas[1])-offset)>>mupdate_sh);
			m3+=(int)((update*(                  alphas[0])*(                  alphas[1])-offset)>>mupdate_sh);
			CLAMP2_32(p0, p0, -(1<<ABAC_PROBBITS>>1)+1, (1<<ABAC_PROBBITS>>1)-1);
			//CLAMP2_32(m0, m0, 0, 0x7FFFFF);
			stats[k][ABAC_NCTRS*tidx+k2]=p0;
			mixer[0][k3]=m0;
			mixer[1][k3]=m1;
			mixer[2][k3]=m2;
			mixer[2][k3]=m3;
		}
	}
#if 0
	const int sh=4;
	//int err=(!bit<<ABAC_PROBBITS)-p0;
	for(int k=0;k<ABAC_NCTX;++k)
	{
		int p0=stats[k][tidx];
		p0+=((!bit<<16)-p0+(1<<sh>>1))>>sh;
	//	CLAMP2_32(p0, p0, 1, 0xFFFF);
		stats[k][tidx]=p0;
	}
#endif
}
#endif
typedef struct _ThreadArgs
{
	const unsigned char *src;
	unsigned char *dst;
	int iw, ih;

	int fwd, test, loud, x1, x2, y1, y2;
	int bufsize, histsize;
	short *pixels;
	int *hist;

	DList list;
	const unsigned char *decstart, *decend;
	
#ifdef ENABLE_MIX4
	int clevels;
#endif
	int tlevels;
#ifdef ENABLE_ABAC2
	unsigned short stats[A2_NCTX*3<<A2_CTXBITS<<8];
#elif defined ENABLE_ABAC
	int statssize;
#ifdef ABAC_PROFILESIZE
	double abac_csizes[ABAC_TOKEN_BITS*3];
#endif
//	int mixer[ABAC_NCTX*3<<ABAC_TOKEN_BITS];
	int mixer[ABAC_NCTRS*ABAC_NCTX*3*((1<<(8-ABAC_MIXBITS))+1)*(1<<(8-ABAC_MIXBITS)|1)];
//	int mixer[ABAC_NCTX*ABAC_TOKEN_BITS*3*2];
//	int mixer[ABAC_NCTX*ABAC_TOKEN_BITS*3];
#endif

	//aux
	int blockidx;
	double bestsize;
	int bestrct, predidx[3];
} ThreadArgs;
//static void check_result(const short *result)
//{
//	for(int k=0;k<15;++k)
//	{
//		if((unsigned)result[k]>255)
//			LOG_ERROR("");
//	}
//}
static void block_thread(void *param)
{
	const int nch=3, depth=8, half=128;
	ThreadArgs *args=(ThreadArgs*)param;
	AC3 ec;
	const unsigned char *image=args->fwd?args->src:args->dst;
	unsigned char bestrct=0, combination[6]={0}, predidx[4]={0};
	int ystride=args->iw*3;
	int cdfstride=args->tlevels+1;
#ifdef ENABLE_MIX4
	int nctx=args->clevels*args->clevels;
	int chsize=nctx*cdfstride;
#elif defined ENABLE_ABAC
	unsigned short *stats=(unsigned short*)args->hist;
#endif
	
//#ifdef AC_VALIDATE
//	acval_disable=1;
//#endif
	if(args->fwd)
	{
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		unsigned char predsel[OCH_COUNT]={0};
		int res=(args->x2-args->x1-3)/5*5*(args->y2-args->y1-2);
		__m256i av12_mcoeffs[12];

		for(int k=0;k<(int)_countof(av12_mcoeffs);++k)
			av12_mcoeffs[k]=_mm256_set1_epi16(av12_icoeffs[k]>>1);
		memset(args->hist, 0, args->histsize);
		for(int ky=args->y1+2;ky<args->y2;++ky)//analysis loop
		{
			int kx=args->x1+2;
			const unsigned char *ptr=image+3*(args->iw*ky+kx);

			__m256i amin=_mm256_set1_epi16(-128);
			__m256i amax=_mm256_set1_epi16(127);
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
			for(;kx<args->x2-5;kx+=5, ptr+=15)
			{
				__m256i
					nb0[13],//rgb
					nb1[13],//gbr
					nb2[13],//rgb - gbr
					nb3[13],//gbr - rgb
					nb4[13],//(gbr+brg)/2
					nb5[13];//rgb - (gbr+brg)/2
				__m256i vmin[4], vmax[4], pred;
				{
					__m128i nb8[13]=//8-bit
					{
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride-2*3+0)), half8),//NNWW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride-1*3+0)), half8),//NNW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride+0*3+0)), half8),//NN
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride+1*3+0)), half8),//NNE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-2*ystride+2*3+0)), half8),//NNEE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride-2*3+0)), half8),//NWW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride-1*3+0)), half8),//NW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+0*3+0)), half8),//N
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+1*3+0)), half8),//NE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr-1*ystride+2*3+0)), half8),//NEE
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride-2*3+0)), half8),//WW
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride-1*3+0)), half8),//W
						_mm_xor_si128(_mm_load_si128((__m128i*)(ptr+0*ystride+0*3+0)), half8),//curr
					};
					for(int k=0;k<13;++k)
					{
						__m128i temp;
						__m256i t2;
						nb0[k]=_mm256_cvtepi8_epi16(nb8[k]);
						temp=_mm_shuffle_epi8(nb8[k], shuf);
						nb1[k]=_mm256_cvtepi8_epi16(temp);
						t2=_mm256_cvtepi8_epi16(_mm_shuffle_epi8(temp, shuf));
						t2=_mm256_add_epi16(t2, nb1[k]);
						t2=_mm256_srai_epi16(t2, 1);
						nb2[k]=_mm256_sub_epi16(nb0[k], nb1[k]);
						nb3[k]=_mm256_sub_epi16(nb1[k], nb0[k]);
						nb4[k]=t2;
						nb5[k]=_mm256_sub_epi16(nb0[k], t2);
					}
				}
#define UPDATE(PREDIDX, IDX0, IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7, IDX8, IDX9, IDXA, IDXB, IDXC, IDXD, IDXE)\
	do\
	{\
		pred=_mm256_slli_epi16(pred, 8);\
		pred=_mm256_srai_epi16(pred, 8);\
		pred=_mm256_sub_epi16(pred, amin);\
		_mm256_store_si256((__m256i*)result, pred);\
		++args->hist[(IDX0*PRED_COUNT+PREDIDX)<<8|result[0x0]];\
		++args->hist[(IDX1*PRED_COUNT+PREDIDX)<<8|result[0x1]];\
		++args->hist[(IDX2*PRED_COUNT+PREDIDX)<<8|result[0x2]];\
		++args->hist[(IDX3*PRED_COUNT+PREDIDX)<<8|result[0x3]];\
		++args->hist[(IDX4*PRED_COUNT+PREDIDX)<<8|result[0x4]];\
		++args->hist[(IDX5*PRED_COUNT+PREDIDX)<<8|result[0x5]];\
		++args->hist[(IDX6*PRED_COUNT+PREDIDX)<<8|result[0x6]];\
		++args->hist[(IDX7*PRED_COUNT+PREDIDX)<<8|result[0x7]];\
		++args->hist[(IDX8*PRED_COUNT+PREDIDX)<<8|result[0x8]];\
		++args->hist[(IDX9*PRED_COUNT+PREDIDX)<<8|result[0x9]];\
		++args->hist[(IDXA*PRED_COUNT+PREDIDX)<<8|result[0xA]];\
		++args->hist[(IDXB*PRED_COUNT+PREDIDX)<<8|result[0xB]];\
		++args->hist[(IDXC*PRED_COUNT+PREDIDX)<<8|result[0xC]];\
		++args->hist[(IDXD*PRED_COUNT+PREDIDX)<<8|result[0xD]];\
		++args->hist[(IDXE*PRED_COUNT+PREDIDX)<<8|result[0xE]];\
	}while(0)
				//W
				pred=_mm256_sub_epi16(nb0[NB_curr], nb0[NB_W]);
				UPDATE(
					PRED_W,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);

				pred=_mm256_add_epi16(nb2[NB_W], nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_W,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);

				pred=_mm256_add_epi16(nb3[NB_W], nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_W,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);

				pred=_mm256_add_epi16(nb5[NB_W], nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_W,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//CG
				vmin[0]=_mm256_min_epi16(nb0[NB_N], nb0[NB_W]);
				vmax[0]=_mm256_max_epi16(nb0[NB_N], nb0[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb0[NB_N], nb0[NB_W]), nb0[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				vmin[1]=_mm256_min_epi16(nb2[NB_N], nb2[NB_W]);
				vmax[1]=_mm256_max_epi16(nb2[NB_N], nb2[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb2[NB_N], nb2[NB_W]), nb2[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				vmin[2]=_mm256_min_epi16(nb3[NB_N], nb3[NB_W]);
				vmax[2]=_mm256_max_epi16(nb3[NB_N], nb3[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb3[NB_N], nb3[NB_W]), nb3[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				vmin[3]=_mm256_min_epi16(nb5[NB_N], nb5[NB_W]);
				vmax[3]=_mm256_max_epi16(nb5[NB_N], nb5[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(nb5[NB_N], nb5[NB_W]), nb5[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_CG,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//AV5
				//		-5	5	1
				//	-1	8	[?]>>3
				pred=_mm256_sub_epi16(nb0[NB_N], nb0[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb0[NB_NE], nb0[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb0[NB_W]);
				//vmin=_mm256_min_epi16(N, W);
				//vmax=_mm256_max_epi16(N, W);
				vmin[0]=_mm256_min_epi16(vmin[0], nb0[NB_NE]);
				vmax[0]=_mm256_max_epi16(vmax[0], nb0[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_sub_epi16(nb2[NB_N], nb2[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb2[NB_NE], nb2[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb2[NB_W]);
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				vmin[1]=_mm256_min_epi16(vmin[1], nb2[NB_NE]);
				vmax[1]=_mm256_max_epi16(vmax[1], nb2[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);

				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_sub_epi16(nb3[NB_N], nb3[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb3[NB_NE], nb3[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb3[NB_W]);
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				vmin[2]=_mm256_min_epi16(vmin[2], nb3[NB_NE]);
				vmax[2]=_mm256_max_epi16(vmax[2], nb3[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);

				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_sub_epi16(nb5[NB_N], nb5[NB_NW]);
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(pred, 2));
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb5[NB_NE], nb5[NB_WW]));
				pred=_mm256_add_epi16(_mm256_srai_epi16(pred, 3), nb5[NB_W]);
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				vmin[3]=_mm256_min_epi16(vmin[1], nb5[NB_NE]);
				vmax[3]=_mm256_max_epi16(vmax[1], nb5[NB_NE]);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);

				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV5,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//AV9
				//		1	-2	-1
				//	-1	-9	10	4
				//	-2	16	[?]>>4
				pred=_mm256_add_epi16(nb0[NB_N], _mm256_slli_epi16(nb0[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb0[NB_NN], nb0[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb0[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb0[NB_NW], 3), nb0[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb0[NB_NNW], _mm256_add_epi16(nb0[NB_NNE], nb0[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb0[NB_W], _mm256_srai_epi16(pred, 4));
				//vmin=_mm256_min_epi16(N, W);
				//vmax=_mm256_max_epi16(N, W);
				//vmin=_mm256_min_epi16(vmin, NE);
				//vmax=_mm256_max_epi16(vmax, NE);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_add_epi16(nb2[NB_N], _mm256_slli_epi16(nb2[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb2[NB_NN], nb2[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb2[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb2[NB_NW], 3), nb2[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb2[NB_NNW], _mm256_add_epi16(nb2[NB_NNE], nb2[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb2[NB_W], _mm256_srai_epi16(pred, 4));
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin=_mm256_min_epi16(vmin, NE3);
				//vmax=_mm256_max_epi16(vmax, NE3);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);
					
				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_add_epi16(nb3[NB_N], _mm256_slli_epi16(nb3[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb3[NB_NN], nb3[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb3[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb3[NB_NW], 3), nb3[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb3[NB_NNW], _mm256_add_epi16(nb3[NB_NNE], nb3[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb3[NB_W], _mm256_srai_epi16(pred, 4));
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin=_mm256_min_epi16(vmin, NE3);
				//vmax=_mm256_max_epi16(vmax, NE3);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);
					
				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_add_epi16(nb5[NB_N], _mm256_slli_epi16(nb5[NB_N], 2));//5*N
				pred=_mm256_sub_epi16(pred, _mm256_add_epi16(nb5[NB_NN], nb5[NB_WW]));//5*N - (NN+WW)
				pred=_mm256_add_epi16(pred, _mm256_slli_epi16(nb5[NB_NE], 1));//5*N-NN-WW + 2*NE
				pred=_mm256_sub_epi16(_mm256_slli_epi16(pred, 1), _mm256_add_epi16(_mm256_slli_epi16(nb5[NB_NW], 3), nb5[NB_NW]));//2*(5*N-NN-WW+2*NE) - 9*NW
				pred=_mm256_add_epi16(pred, _mm256_sub_epi16(nb5[NB_NNW], _mm256_add_epi16(nb5[NB_NNE], nb5[NB_NWW])));//2*(5*N-NN-WW+2*NE)-9*NW + NNW-NNE-NWW
				pred=_mm256_add_epi16(nb5[NB_W], _mm256_srai_epi16(pred, 4));
				//vmin=_mm256_min_epi16(N3, W3);
				//vmax=_mm256_max_epi16(N3, W3);
				//vmin=_mm256_min_epi16(vmin, NE3);
				//vmax=_mm256_max_epi16(vmax, NE3);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);
					
				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV9,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);

				//AV12
				pred=_mm256_setzero_si256();
				for(int k=0;k<(int)_countof(av12_icoeffs);++k)
					pred=_mm256_add_epi16(pred, _mm256_mullo_epi16(av12_mcoeffs[k], nb0[k]));
				pred=_mm256_srai_epi16(pred, 7);
				pred=_mm256_max_epi16(pred, vmin[0]);
				pred=_mm256_min_epi16(pred, vmax[0]);

				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV12,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B,
					OCH_R, OCH_G, OCH_B
				);
				pred=_mm256_setzero_si256();
				for(int k=0;k<(int)_countof(av12_icoeffs);++k)
					pred=_mm256_add_epi16(pred, _mm256_mullo_epi16(av12_mcoeffs[k], nb2[k]));
				pred=_mm256_srai_epi16(pred, 7);
				pred=_mm256_max_epi16(pred, vmin[1]);
				pred=_mm256_min_epi16(pred, vmax[1]);
					
				pred=_mm256_add_epi16(pred, nb1[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV12,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR,
					OCH_RG, OCH_GB, OCH_BR
				);
				pred=_mm256_setzero_si256();
				for(int k=0;k<(int)_countof(av12_icoeffs);++k)
					pred=_mm256_add_epi16(pred, _mm256_mullo_epi16(av12_mcoeffs[k], nb3[k]));
				pred=_mm256_srai_epi16(pred, 7);
				pred=_mm256_max_epi16(pred, vmin[2]);
				pred=_mm256_min_epi16(pred, vmax[2]);
					
				pred=_mm256_add_epi16(pred, nb0[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb1[NB_curr], pred);
				UPDATE(
					PRED_AV12,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB,
					OCH_GR, OCH_BG, OCH_RB
				);
				pred=_mm256_setzero_si256();
				for(int k=0;k<(int)_countof(av12_icoeffs);++k)
					pred=_mm256_add_epi16(pred, _mm256_mullo_epi16(av12_mcoeffs[k], nb5[k]));
				pred=_mm256_srai_epi16(pred, 7);
				pred=_mm256_max_epi16(pred, vmin[3]);
				pred=_mm256_min_epi16(pred, vmax[3]);
					
				pred=_mm256_add_epi16(pred, nb4[NB_curr]);
				pred=_mm256_max_epi16(pred, amin);
				pred=_mm256_min_epi16(pred, amax);
				pred=_mm256_sub_epi16(nb0[NB_curr], pred);
				UPDATE(
					PRED_AV12,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2,
					OCH_R2, OCH_G2, OCH_B2
				);
			}
		}
		for(int kc=0;kc<OCH_COUNT*PRED_COUNT;++kc)
		{
			int *curr_hist=args->hist+((size_t)kc<<8);
			for(int ks=0;ks<256;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
					csizes[kc]-=freq*log2((double)freq/res);
			}
			csizes[kc]/=8;
			//if(kc>=OCH_R2*PRED_COUNT)
			//	csizes[kc]=INFINITY;
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
				csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
				csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
				csizes[group[2]*PRED_COUNT+predsel[group[2]]];
			if(!kt||bestsize>csize)
				bestsize=csize, bestrct=kt;
		}
		memcpy(combination, rct_combinations[bestrct], sizeof(combination));
		predidx[0]=predsel[combination[0]];
		predidx[1]=predsel[combination[1]];
		predidx[2]=predsel[combination[2]];
		args->bestsize=bestsize;
		args->bestrct=bestrct;
		args->predidx[0]=predidx[0];
		args->predidx[1]=predidx[1];
		args->predidx[2]=predidx[2];
		if(args->loud)
		{
			printf("Y %5d~%5d  best %12.2lf bytes  %s [YUV: %s %s %s]\n",
				args->y1, args->y2,
				bestsize,
				rct_names[bestrct],
				pred_names[predidx[0]],
				pred_names[predidx[1]],
				pred_names[predidx[2]]
			);

			for(int kp=0;kp<PRED_COUNT;++kp)
				printf("%14s", pred_names[kp]);
			printf("\n");
			for(int kc=0;kc<OCH_COUNT;++kc)
			{
				for(int kp=0;kp<PRED_COUNT;++kp)
					printf(" %12.2lf%c", csizes[kc*PRED_COUNT+kp], kp==predsel[kc]?'*':' ');
				printf("  %s\n", och_names[kc]);
			}

			for(int kt=0;kt<RCT_COUNT;++kt)
			{
				const unsigned char *group=rct_combinations[kt];
				double csize=
					csizes[group[0]*PRED_COUNT+predsel[group[0]]]+
					csizes[group[1]*PRED_COUNT+predsel[group[1]]]+
					csizes[group[2]*PRED_COUNT+predsel[group[2]]];
				printf("%12.2lf %c  %-10s %-10s %-10s %-10s\n",
					csize,
					kt==bestrct?'*':' ',
					rct_names[kt],
					pred_names[predsel[group[0]]],
					pred_names[predsel[group[1]]],
					pred_names[predsel[group[2]]]
				);
			}
		}
		dlist_init(&args->list, 1, BLOCKSIZE*BLOCKSIZE*3, 0);
		ac3_enc_init(&ec, &args->list);
		ac3_enc_bypass_NPOT(&ec, bestrct, RCT_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[0], PRED_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[1], PRED_COUNT);
		ac3_enc_bypass_NPOT(&ec, predidx[2], PRED_COUNT);
	}
	else
	{
		ac3_dec_init(&ec, args->decstart, args->decend);
		bestrct=ac3_dec_bypass_NPOT(&ec, RCT_COUNT);
		predidx[0]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[1]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predidx[2]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
	}
#if defined ENABLE_CALICCTX || defined ENABLE_HASH || defined ENABLE_MIX8
	{
		static const int init_freqs[]={32, 8, 6, 4, 3, 2, 1};
		int sum=0;
		for(int ks=0;ks<args->tlevels;++ks)
			sum+=args->hist[ks]=init_freqs[MINVAR(ks, (int)_countof(init_freqs)-1)];
		args->hist[args->tlevels]=sum;
		memfill(args->hist+cdfstride, args->hist, args->histsize-cdfstride*sizeof(int), cdfstride*sizeof(int));
	}
#elif defined ENABLE_ABAC2
	int mixer[8*3*A2_NCTX]={0};
	FILLMEM(mixer, 0x8000, sizeof(mixer), sizeof(int));
	FILLMEM(args->stats, 0x8000, sizeof(args->stats), sizeof(short));
#elif defined ENABLE_ABAC
	{
		static const unsigned short prob0[]=//NPOT tree initialized to bypass
		{
			//0,
			//(16<<16)/25-0x8000,
			//0, (8<<16)/9-0x8000,
			//0, 0, 0, 0,
			//0, 0, 0, 0, 0, 0, 0, 0,
			//0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			
			(1<<16)/2,
			(16<<16)/25,
			(1<<16)/2, (8<<16)/9,
			(1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2,
			(1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2,
			(1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2, (1<<16)/2,
		};
		memfill(stats, prob0, args->statssize, sizeof(prob0));
	}
	//FILLMEM(stats, 0x8000, args->statssize, sizeof(short));
	memset(args->mixer, 0, sizeof(args->mixer));
#ifdef ABAC_PROFILESIZE
	memset(args->abac_csizes, 0, sizeof(args->abac_csizes));
#endif
#else
	for(int ky=0;ky<args->clevels;++ky)
	{
		for(int kx=0;kx<args->clevels;++kx)
		{
			static const int init_freqs[]={32, 8, 6, 4, 3, 2, 1};
			int *curr_hist=args->hist+cdfstride*(args->clevels*ky+kx);
			int sum=0;
			for(int ks=0;ks<args->tlevels;++ks)
			{
				int freq=init_freqs[MINVAR(ks, (int)_countof(init_freqs)-1)];
				sum+=curr_hist[ks]=freq;
			}
			curr_hist[args->tlevels]=sum;
		}
	}
	memfill(args->hist+chsize, args->hist, sizeof(int)*chsize*(nch-1LL), sizeof(int)*chsize);
#endif
	memset(args->pixels, 0, args->bufsize);
#ifdef ABAC_SIMPLEOVERRIDE
	unsigned short stats2[768]={0};
	FILLMEM(stats2, 0x8000, sizeof(short[768]), sizeof(short));
#endif
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) short *rows[]=
		{
			args->pixels+((BLOCKSIZE+16LL)*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKSIZE+16LL)*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKSIZE+16LL)*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKSIZE+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		int yuv[4]={0};
		int token=0, bypass=0, nbits=0;
		int pred=0, error=0, sym=0;
		const unsigned char *combination=rct_combinations[bestrct];
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int idx=nch*(args->iw*ky+kx);
			short
				*NNN	=rows[3]+0*4*2,
				*NNWW	=rows[2]-2*4*2,
				*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NNEE	=rows[2]+2*4*2,
				*NNEEE	=rows[2]+3*4*2,
				*NWW	=rows[1]-2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*NEE	=rows[1]+2*4*2,
				*NEEE	=rows[1]+3*4*2,
				*WWWW	=rows[0]-4*4*2,
				*WWW	=rows[0]-3*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
			if(ky<=args->y1+2)
			{
				if(ky<=args->y1+1)
				{
					if(ky==args->y1)
						NEEE=NEE=NE=NWW=NW=N=W;
					NNWW=NWW;
					NNW=NW;
					NN=N;
					NNE=NE;
					NNEE=NEE;
					NNEEE=NEEE;
				}
				NNN=NN;
			}
			if(kx<=args->x1+3)
			{
				if(kx<=args->x1+2)
				{
					if(kx<=args->x1+1)
					{
						if(kx<=args->x1)
							NW=W=N;
						WW=W;
						NWW=NW;
					}
					WWW=WW;
				}
				WWWW=WWW;
			}
			if(kx>=args->x2-3)
			{
				if(kx>=args->x2-2)
				{
					if(kx>=args->x2-1)
					{
						NNE=NN;
						NE=N;
					}
					NNEE=NNE;
					NEE=NE;
				}
				NEEE=NEE;
			}
			if(args->fwd)
			{
				switch(bestrct)
				{
				case RCT_R_G_B:
				case RCT_R_G_BG:
				case RCT_R_G_BR:
				case RCT_R_GR_BR:
				case RCT_R_GR_BG:
				case RCT_R_G_B2:
				case RCT_R_GR_B2:
					yuv[0]=args->src[idx+0]-128;
					yuv[1]=args->src[idx+1]-128;
					yuv[2]=args->src[idx+2]-128;
					break;
				case RCT_G_B_RG:
				case RCT_G_B_RB:
				case RCT_G_BG_RG:
				case RCT_G_BG_RB:
				case RCT_G_B_R2:
				case RCT_G_BG_R2:
					yuv[0]=args->src[idx+1]-128;
					yuv[1]=args->src[idx+2]-128;
					yuv[2]=args->src[idx+0]-128;
					break;
				case RCT_B_R_GR:
				case RCT_B_R_GB:
				case RCT_B_RB_GB:
				case RCT_B_RB_GR:
				case RCT_B_RB_G2:
					yuv[0]=args->src[idx+2]-128;
					yuv[1]=args->src[idx+0]-128;
					yuv[2]=args->src[idx+1]-128;
					break;
				case RCT_G_RG_BR:
				case RCT_G_RG_B2:
					yuv[0]=args->src[idx+1]-128;
					yuv[1]=args->src[idx+0]-128;
					yuv[2]=args->src[idx+2]-128;
					break;
				case RCT_B_GB_RG:
				case RCT_B_GB_R2:
					yuv[0]=args->src[idx+2]-128;
					yuv[1]=args->src[idx+1]-128;
					yuv[2]=args->src[idx+0]-128;
					break;
				case RCT_R_BR_GB:
				case RCT_R_B_G2:
				case RCT_R_BR_G2:
					yuv[0]=args->src[idx+0]-128;
					yuv[1]=args->src[idx+2]-128;
					yuv[2]=args->src[idx+1]-128;
					break;
				}
			}
			for(int kc=0;kc<nch;++kc)
			{
				int kc2=kc<<1;
				int offset=(yuv[combination[kc+3]]+yuv[combination[kc+6]])>>combination[kc+9];

				switch(predidx[kc])
				{
				case PRED_W:
					pred=W[kc2];
					break;
				case PRED_CG:
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					break;
				case PRED_AV5:
					CLAMP3_32(pred, W[kc2]+((5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])>>3), N[kc2], W[kc2], NE[kc2]);
					break;
				case PRED_AV9:
					CLAMP3_32(pred,
						W[kc2]+((10*N[kc2]-9*NW[kc2]+4*NE[kc2]-2*(NN[kc2]+WW[kc2])+NNW[kc2]-(NNE[kc2]+NWW[kc2]))>>4),
						N[kc2],
						W[kc2],
						NE[kc2]
					);
					break;
				case PRED_AV12:
					pred=(
						av12_icoeffs[ 0]*NNWW[kc2]+
						av12_icoeffs[ 1]*NNW[kc2]+
						av12_icoeffs[ 2]*NN[kc2]+
						av12_icoeffs[ 3]*NNE[kc2]+
						av12_icoeffs[ 4]*NNEE[kc2]+
						av12_icoeffs[ 5]*NWW[kc2]+
						av12_icoeffs[ 6]*NW[kc2]+
						av12_icoeffs[ 7]*N[kc2]+
						av12_icoeffs[ 8]*NE[kc2]+
						av12_icoeffs[ 9]*NEE[kc2]+
						av12_icoeffs[10]*WW[kc2]+
						av12_icoeffs[11]*W[kc2]
					)>>8;
					CLAMP3_32(pred, pred, N[kc2], W[kc2], NE[kc2]);
					break;
				}
				pred+=offset;
				CLAMP2_32(pred, pred, -128, 127);
#ifdef ENABLE_ABAC2
				unsigned short *curr_stats[A2_NCTX];
				int ctx[A2_NCTX]=
				{
					0,
					pred,
					//N[kc2+0],
					W[kc2+0]/32,
					N[kc2+0]-NW[kc2+0],
					NE[kc2+0]-N[kc2+0],
					W[kc2+0]-NW[kc2+0],
					N[kc2+0]-W[kc2+0],
					//N[kc2+0]+W[kc2+0]-NW[kc2+0],
					//N[kc+0]+NN[kc2+0],
					//W[kc+0]+WW[kc2+0],
					//NW[kc2+1],
					N[kc2+0]+W[kc2+0]-NW[kc2+0],
				};
				for(int k=0;k<A2_NCTX;++k)
				{
					int idx2=((A2_NCTX*kc+k)<<A2_CTXBITS|((ctx[k]+256)>>(9-A2_CTXBITS)&((1<<A2_CTXBITS)-1)))<<8;
					curr_stats[k]=args->stats+idx2;
				}
				if(args->fwd)
				{
					curr[kc2+0]=yuv[kc];
					curr[kc2+1]=error=yuv[kc]-pred;
				}
				else
					error=0;
				int tidx=1;
				int *curr_mixer=mixer+kc*8*A2_NCTX;
				for(int kb=7, e2=0;kb>=0;--kb)
				{
					long long p0=0;
					int wsum=0, bit;
					for(int k=0;k<A2_NCTX;++k)
					{
						p0+=(long long)curr_mixer[k]*curr_stats[k][tidx];
						wsum+=curr_mixer[k];
					}
					if(wsum)
						p0/=wsum;
					else
						p0=0x8000, wsum=1;
					int p00=(int)p0;
					CLAMP2_32(p0, (int)p0, 1, 0xFFFF);
					if(args->fwd)
					{
						bit=error>>kb&1;
						ac3_enc_bin(&ec, bit, (int)p0, 16);
					}
					else
					{
						bit=ac3_dec_bin(&ec, (int)p0, 16);
						error|=bit<<kb;
					}
					e2|=bit<<kb;
					int pred2=(char)(e2|1<<kb>>1)-pred;
					//if((unsigned)(p00-1)<0xFFFE)
					if(abs((!bit<<16)-(int)p0)>256)
					{
						int pbit=bit?0x10000-(int)p0:(int)p0;
						long long dL_dp0=-(1LL<<32)/pbit;//fixed 47.16 bit
						dL_dp0^=-bit;
						dL_dp0+=bit;
						for(int k=0;k<A2_NCTX;++k)
						{
							int diff=curr_stats[k][tidx]-(int)p0;
							long long grad=dL_dp0*diff/wsum;
							long long wnew=grad*2750>>16;
							wnew=curr_mixer[k]-wnew;
							CLAMP2_32(curr_mixer[k], (int)wnew, 1, 0xFFFF);
						}
					}
					{
						int sh;
						sh=(abs(N[kc2+1])+abs(W[kc2+1])+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1]));//{76543}
						switch(kb)
						{
						default:
							sh=sh*9>>8;
							break;
						case 3:
						case 2:
							sh=(9*sh+abs(pred2))>>8;
							//sh=(15*sh/28+abs(pred2)/16)>>4;
							break;
						case 1:
							sh=(8*sh+4*abs(pred2))>>8;
							//sh=(7*sh/28+abs(pred2)/8)>>3;
							break;
						case 0:
							sh=abs(pred2)>>5;//{0}
							break;
						}
						sh=FLOOR_LOG2(sh+1);
						for(int k=0;k<A2_NCTX;++k)
						{
							int p=curr_stats[k][tidx];
							p+=((!bit<<16)-p)>>(5+sh);	//5
							CLAMP2_32(curr_stats[k][tidx], p, 1, 0xFFFF);
						}
					}
					curr_mixer+=A2_NCTX;
					tidx+=tidx+bit;
				}
				if(!args->fwd)
				{
					error+=pred;
					error=error<<(32-8)>>(32-8);
					curr[kc2+0]=yuv[kc]=error;
					curr[kc2+1]=error-pred;
				}
#elif defined ENABLE_ABAC
				int tidx, token2;
				int logits[ABAC_NCTRS*ABAC_NCTX];
				unsigned short *curr_stats[ABAC_NCTX];
				int ctx[]=
				{
					//pixel contexts
#if 0
					//kx>>5,
					pred,
					//N[kc2]+W[kc2]-NW[kc2],
					//NNWW	[kc2],
					//NNW	[kc2],
					//NN	[kc2],
					//NNE	[kc2],
					//NNEE	[kc2],
					//NNEEE	[kc2],
					//NWW	[kc2],
					//NW	[kc2],
					//N	[kc2],
					//NE	[kc2],
					//NEE	[kc2],
					//NEEE	[kc2],
					//WWWW	[kc2],
					//WWW	[kc2],
					//WW	[kc2],
					//W	[kc2],
#endif
#if 1
					kx>>5,
					pred,
					N[kc2+0]*5,
					W[kc2+0],
					N[kc2+0]+W[kc2+0]-NW[kc2+0]-NE[kc2+0],
					(N[kc2+0]+W[kc2+0]-NW[kc2+0])*2,
					(W[kc2+0]+NE[kc2+0]-N[kc2+0])/4,
					(N[kc2+0]+NE[kc2+0]-NNE[kc2+0])*3,
					2*N[kc2+0]+NE[kc2+0]-2*NNE[kc2+0],
					N[kc2+0]-NN[kc2+0]-NE[kc2+0]+NNE[kc2+0],
					(NE[kc2+0]+NEE[kc2+0]-NNEEE[kc2+0])/2,
					W[kc2+0]+NW[kc2+0]-NWW[kc2+0],
					W[kc2+0]-WW[kc2+0],
					N[kc2+0]-NN[kc2+0],
					3*(N[kc2+0]-NN[kc2+0])+NNN[kc2+0],
					3*(W[kc2+0]-WW[kc2+0])+WWW[kc2+0],
					W[kc2+0]+N[kc2+0]+NE[kc2+0]+NEE[kc2+0]+NEEE[kc2+0],
					(WWWW[kc2+0]-NEEE[kc2+0])>7,
					(WWW[kc2+0]-NEE[kc2+0])>6,
					NW[kc2+0]-W[kc2+0],
					NW[kc2+0]-N[kc2+0],
					NW[kc2+0]-NNWW[kc2+0],
					NNE[kc2+0]-NN[kc2+0],
#endif

					//error contexts
#if 0
					//N[kc2+1]+W[kc2+1]-NW[kc2+1],
					//NNWW	[kc2+1],
					//NNW	[kc2+1],
					//NN	[kc2+1],
					//NNE	[kc2+1],
					//NNEE	[kc2+1],
					//NNEEE	[kc2+1],
					//NWW	[kc2+1],
					//NW	[kc2+1],
					N	[kc2+1],
					//NE	[kc2+1],
					//NEE	[kc2+1],
					//NEEE	[kc2+1],
					//WWWW	[kc2+1],
					//WWW	[kc2+1],
					//WW	[kc2+1],
					W	[kc2+1],
#endif
#if 1
					N[kc2+1],
					W[kc2+1]*3,
					N[kc2+1]*17-NN[kc2+1],
					W[kc2+1]*9-WW[kc2+1]/4,
					((N[kc2+1]-NN[kc2+1])*3+NNN[kc2+1])/6,
					((W[kc2+1]-WW[kc2+1])*3+WWW[kc2+1])*6,
					(2*abs(N[kc2+1])+abs(W[kc2+1])-abs(NW[kc2+1])+abs(NE[kc2+1]))/32,
					(2*abs(N[kc2+1])+abs(W[kc2+1])+abs(NW[kc2+1])+abs(NE[kc2+1]))/2,
					2*N[kc2+1]+NN[kc2+1],
					(W[kc2+1]+WW[kc2+1]+WWW[kc2+1])/32,
					N[kc2+1]+NN[kc2+1]+NE[kc2+1]+NNE[kc2+1],
					W[kc2+1]+WW[kc2+1]+NW[kc2+1]+NWW[kc2+1],
					6*N[kc2+1]-NNW[kc2+1]-NN[kc2+1]-NNE[kc2+1],
#endif
				};
				for(int k=0;k<ABAC_PCTX;++k)
				{
					int x=ctx[k];
					x>>=3;
					CLAMP2_32(x, x, -ABAC_CLEVELS/2, ABAC_CLEVELS/2-1);
					ctx[k]=x&(ABAC_CLEVELS-1);
				}
				for(int k=ABAC_PCTX;k<ABAC_NCTX;++k)
				{
					int x=ctx[k], neg=x<0;
					x=abs(x);
					x=FLOOR_LOG2(x+1);
					if(!neg)
						x+=ABAC_CLEVELS/2;
					//x>>=3;//X
					ctx[k]=x&(ABAC_CLEVELS-1);
				}
				memset(ctx, 0, sizeof(ctx));//
				for(int k=0;k<ABAC_NCTX;++k)
					curr_stats[k]=stats+(k*ABAC_CLEVELS+ctx[k])*ABAC_TREESIZE;
				
				int mixalphas[]=
				{
					(W[kc2+1]+128)&((1<<ABAC_MIXBITS)-1),
					(N[kc2+1]+128)&((1<<ABAC_MIXBITS)-1),
				};
				int mixidx[]=
				{
					(W[kc2+1]+128)&255>>ABAC_MIXBITS,
					(N[kc2+1]+128)&255>>ABAC_MIXBITS,
				};
				int *curr_mixer[]=
				{
					args->mixer+ABAC_NCTRS*ABAC_NCTX*((1LL<<(8-ABAC_MIXBITS)|1)*((1LL<<(8-ABAC_MIXBITS)|1)*kc+mixidx[1]+0)+mixidx[0]+0),
					args->mixer+ABAC_NCTRS*ABAC_NCTX*((1LL<<(8-ABAC_MIXBITS)|1)*((1LL<<(8-ABAC_MIXBITS)|1)*kc+mixidx[1]+0)+mixidx[0]+1),
					args->mixer+ABAC_NCTRS*ABAC_NCTX*((1LL<<(8-ABAC_MIXBITS)|1)*((1LL<<(8-ABAC_MIXBITS)|1)*kc+mixidx[1]+1)+mixidx[0]+0),
					args->mixer+ABAC_NCTRS*ABAC_NCTX*((1LL<<(8-ABAC_MIXBITS)|1)*((1LL<<(8-ABAC_MIXBITS)|1)*kc+mixidx[1]+1)+mixidx[0]+1),
				};
				int weights[ABAC_NCTRS*ABAC_NCTX];

			//	int *curr_mixer=args->mixer+(ABAC_NCTX*(size_t)kc<<ABAC_TOKEN_BITS);
			//	int *curr_mixer=args->mixer+ABAC_NCTX*(size_t)kc*ABAC_TOKEN_BITS*2;
			//	int *curr_mixer=args->mixer+ABAC_NCTX*(size_t)kc*ABAC_TOKEN_BITS;

#elif defined ENABLE_MIX4
				int cdf, freq=0, den;
				int
					vx=(abs(W[kc2]-WW[kc2])+abs(N[kc2]-NW[kc2])+abs(NE[kc2]-N  [kc2])+abs(WWW[kc2+1])+abs(WW[kc2+1])+abs(W[kc2+1])*2)<<10>>depth,
					vy=(abs(N[kc2]-NN[kc2])+abs(W[kc2]-NW[kc2])+abs(NE[kc2]-NNE[kc2])+abs(NNN[kc2+1])+abs(NN[kc2+1])+abs(N[kc2+1])*2)<<10>>depth;
				int qeN=FLOOR_LOG2(vy+1);
				int qeW=FLOOR_LOG2(vx+1);
				int *curr_hist[4];
				int alphax, alphay;

				qeN=MINVAR(qeN, CLEVELS-2);
				qeW=MINVAR(qeW, CLEVELS-2);
				curr_hist[0]=args->hist+cdfstride*(args->clevels*(args->clevels*kc+qeN+0)+qeW+0);
				curr_hist[1]=args->hist+cdfstride*(args->clevels*(args->clevels*kc+qeN+0)+qeW+1);
				curr_hist[2]=args->hist+cdfstride*(args->clevels*(args->clevels*kc+qeN+1)+qeW+0);
				curr_hist[3]=args->hist+cdfstride*(args->clevels*(args->clevels*kc+qeN+1)+qeW+1);
				alphax=(((vx+1-(1<<qeW))<<MIXBITS)+(1<<qeW>>1))>>qeW;
				alphay=(((vy+1-(1<<qeN))<<MIXBITS)+(1<<qeN>>1))>>qeN;
				CLAMP2_32(alphax, 0, alphax, 1<<MIXBITS);
				CLAMP2_32(alphay, 0, alphay, 1<<MIXBITS);
#define MIXCDF(X) f28_mix4(curr_hist[0][X], curr_hist[1][X], curr_hist[2][X], curr_hist[3][X], alphax, alphay)
				den=MIXCDF(args->tlevels);
#elif defined ENABLE_MIX8
				int cdf, freq=0, den;
				int ctx[3]=
				{
				//	W[kc2+0],
				//	(N[kc2+0]+W[kc2+0])>>1,
					N[kc2+0]+W[kc2+0]-NW[kc2+0],
					abs(N[kc2+0]-W[kc2+0]),
					(abs(N[kc2+1])+abs(W[kc2+1])+abs(NW[kc2+1])+abs(NE[kc2+1]))>>2,
				};
				int qctx[3]=
				{
					ctx[0]>>MIXBITSX&(CXLEVELS-2),
					FLOOR_LOG2(ctx[1]+1),
					FLOOR_LOG2(ctx[2]+1),
				};
				int alphas[3]=
				{
					ctx[0]&((1<<MIXBITSX)-1),
					(((ctx[1]+1-(1<<qctx[1]))<<MIXBITSY)+(1<<qctx[1]>>1))>>qctx[1],
					(((ctx[2]+1-(1<<qctx[2]))<<MIXBITSZ)+(1<<qctx[2]>>1))>>qctx[2],
				};
				UPDATE_MIN(qctx[1], CYLEVELS-2);
				UPDATE_MIN(qctx[2], CZLEVELS-2);
				int *curr_hist[8]=
				{
					args->hist+cdfstride*(CXLEVELS*(CYLEVELS*(CZLEVELS*kc+qctx[2]+0)+qctx[1]+0)+qctx[0]+0),
					args->hist+cdfstride*(CXLEVELS*(CYLEVELS*(CZLEVELS*kc+qctx[2]+0)+qctx[1]+0)+qctx[0]+1),
					args->hist+cdfstride*(CXLEVELS*(CYLEVELS*(CZLEVELS*kc+qctx[2]+0)+qctx[1]+1)+qctx[0]+0),
					args->hist+cdfstride*(CXLEVELS*(CYLEVELS*(CZLEVELS*kc+qctx[2]+0)+qctx[1]+1)+qctx[0]+1),
					args->hist+cdfstride*(CXLEVELS*(CYLEVELS*(CZLEVELS*kc+qctx[2]+1)+qctx[1]+0)+qctx[0]+0),
					args->hist+cdfstride*(CXLEVELS*(CYLEVELS*(CZLEVELS*kc+qctx[2]+1)+qctx[1]+0)+qctx[0]+1),
					args->hist+cdfstride*(CXLEVELS*(CYLEVELS*(CZLEVELS*kc+qctx[2]+1)+qctx[1]+1)+qctx[0]+0),
					args->hist+cdfstride*(CXLEVELS*(CYLEVELS*(CZLEVELS*kc+qctx[2]+1)+qctx[1]+1)+qctx[0]+1),
				};
				int dens[8]=
				{
					(1<<24)/curr_hist[0][args->tlevels],
					(1<<24)/curr_hist[1][args->tlevels],
					(1<<24)/curr_hist[2][args->tlevels],
					(1<<24)/curr_hist[3][args->tlevels],
					(1<<24)/curr_hist[4][args->tlevels],
					(1<<24)/curr_hist[5][args->tlevels],
					(1<<24)/curr_hist[6][args->tlevels],
					(1<<24)/curr_hist[7][args->tlevels],
				};
#define MIXCDF(X)\
	f28_mix8(\
		(int)((long long)curr_hist[0][X]*dens[0]>>(24-14)),\
		(int)((long long)curr_hist[1][X]*dens[1]>>(24-14)),\
		(int)((long long)curr_hist[2][X]*dens[2]>>(24-14)),\
		(int)((long long)curr_hist[3][X]*dens[3]>>(24-14)),\
		(int)((long long)curr_hist[4][X]*dens[4]>>(24-14)),\
		(int)((long long)curr_hist[5][X]*dens[5]>>(24-14)),\
		(int)((long long)curr_hist[6][X]*dens[6]>>(24-14)),\
		(int)((long long)curr_hist[7][X]*dens[7]>>(24-14)),\
		alphas[0],\
		alphas[1],\
		alphas[2]\
	)

//#define MIXCDF(X)\
//	f28_mix8(\
//		(curr_hist[0][X]<<14)/curr_hist[0][args->tlevels],\
//		(curr_hist[1][X]<<14)/curr_hist[1][args->tlevels],\
//		(curr_hist[2][X]<<14)/curr_hist[2][args->tlevels],\
//		(curr_hist[3][X]<<14)/curr_hist[3][args->tlevels],\
//		(curr_hist[4][X]<<14)/curr_hist[4][args->tlevels],\
//		(curr_hist[5][X]<<14)/curr_hist[5][args->tlevels],\
//		(curr_hist[6][X]<<14)/curr_hist[6][args->tlevels],\
//		(curr_hist[7][X]<<14)/curr_hist[7][args->tlevels],\
//		alphas[0],\
//		alphas[1],\
//		alphas[2]\
//	)

//#define MIXCDF(X)\
//	f28_mix8(\
//		curr_hist[0][X],\
//		curr_hist[1][X],\
//		curr_hist[2][X],\
//		curr_hist[3][X],\
//		curr_hist[4][X],\
//		curr_hist[5][X],\
//		curr_hist[6][X],\
//		curr_hist[7][X],\
//		alphas[0],\
//		alphas[1],\
//		alphas[2]\
//	)
				//if((unsigned)qctx[0]>=CXLEVELS||(unsigned)qctx[1]>=CYLEVELS||(unsigned)qctx[2]>=CZLEVELS)
				//	LOG_ERROR("");
				//for(int k=0;k<8;++k)
				//{
				//	int *hist2=curr_hist[k];
				//	if(hist2<args->hist||hist2+args->tlevels>=args->hist+args->histsize)
				//		LOG_ERROR("");
				//}

				//den=1<<14;
				den=MIXCDF(args->tlevels);
#elif defined ENABLE_CALICCTX
				int cdf, freq=0, den;
				int *curr_hist;
				int ctx=0;
			//	ctx=ctx*3+THREEWAY(N[kc2+1], 0)+1;
			//	ctx=ctx*3+THREEWAY(W[kc2+1], 0)+1;
				ctx=ctx*3+THREEWAY((N[kc2+1]+W[kc2+1])/8, 0)+1;
				ctx=ctx*3+THREEWAY(NW[kc2]/4, (N[kc2]+W[kc2])/(2*4))+1;
				ctx=ctx*3+THREEWAY(N[kc2]/8, W[kc2]/8)+1;
			//	ctx=ctx*3+THREEWAY(N[kc2], NE[kc2])+1;
			//	ctx=ctx*3+THREEWAY(N[kc2], NW[kc2])+1;
			//	ctx=ctx*3+THREEWAY(W[kc2], NW[kc2])+1;
				curr_hist=args->hist+cdfstride*(CLEVELS*kc+ctx);
#define MIXCDF(X) curr_hist[X]
				//int
				//	alphax=abs(W[kc2]-WW[kc2])+abs(N[kc2]-NW[kc2])+abs(NE[kc2]-N  [kc2])+abs(WW[kc2+1])+abs(W[kc2+1]),
				//	alphay=abs(N[kc2]-NN[kc2])+abs(W[kc2]-NW[kc2])+abs(NE[kc2]-NNE[kc2])+abs(NN[kc2+1])+abs(N[kc2+1]);
				//int *curr_hist[4];
				//
				//curr_hist[0]=args->hist+cdfstride*(nctx*kc+args->clevels*(0)+0);
				//curr_hist[1]=args->hist+cdfstride*(nctx*kc+args->clevels*(0)+1);
				//curr_hist[2]=args->hist+cdfstride*(nctx*kc+args->clevels*(1)+0);
				//curr_hist[3]=args->hist+cdfstride*(nctx*kc+args->clevels*(1)+1);
//#define MIXCDF(X) f28_mix4(curr_hist[0][X], curr_hist[1][X], curr_hist[2][X], curr_hist[3][X], alphax, alphay)
				den=MIXCDF(args->tlevels);
#elif defined ENABLE_HASH
				int cdf, freq=0, den;
				int *curr_hist;
				int ctx=0;//int32_t!

			//	ctx=HASH_CANTOR(ctx, NNWW[kc2])>>8;
			//	ctx=HASH_CANTOR(ctx, NNW[kc2])>>8;
			//	ctx=HASH_CANTOR(ctx, NN[kc2])>>8;
			//	ctx=HASH_CANTOR(ctx, NNE[kc2])>>8;
			//	ctx=HASH_CANTOR(ctx, NNEE[kc2])>>8;
			//	ctx=HASH_CANTOR(ctx, NWW[kc2])>>8;
			//	ctx=HASH_CANTOR(ctx, NW[kc2])>>8;
			//	ctx=HASH_CANTOR(ctx, NE[kc2])>>8;
			//	ctx=HASH_CANTOR(ctx, NEE[kc2])>>8;
			//	ctx=HASH_CANTOR(ctx, WW[kc2])>>8;
			//	ctx=HASH_CANTOR(ctx, N[kc2])>>8;
			//	ctx=HASH_CANTOR(ctx, W[kc2])>>8;
			//	ctx=HASH_CANTOR(N[kc2], W[kc2])>>8;

				int v0=(3*(W[kc2+0]-WW[kc2+0])+WWW[kc2+0])>>9;
				int v1=(3*(N[kc2+0]-NN[kc2+0])+NNN[kc2+0])>>9;
				int v2=(3*(W[kc2+1]-WW[kc2+1])+WWW[kc2+1])>>9;
				int v3=(3*(N[kc2+1]-NN[kc2+1])+NNN[kc2+1])>>9;
				v0=HASH_CANTOR(v0, v1);
				v2=HASH_CANTOR(v2, v3);
				ctx=HASH_CANTOR(v0, v2);

				ctx&=CLEVELS-1;
				curr_hist=args->hist+cdfstride*(CLEVELS*kc+ctx);
#define MIXCDF(X) curr_hist[X]
#else
				int cdf, freq=0, den;
				int
					vx=(abs(W[kc2]-WW[kc2])+abs(N[kc2]-NW[kc2])+abs(NE[kc2]-N  [kc2])+abs(WWW[kc2+1])+abs(WW[kc2+1])+abs(W[kc2+1])*2)<<10>>depth,
					vy=(abs(N[kc2]-NN[kc2])+abs(W[kc2]-NW[kc2])+abs(NE[kc2]-NNE[kc2])+abs(NNN[kc2+1])+abs(NN[kc2+1])+abs(N[kc2+1])*2)<<10>>depth;
				int qeN=FLOOR_LOG2(vy+1);
				int qeW=FLOOR_LOG2(vx+1);
				qeN=MINVAR(qeN, CLEVELS-1);
				qeW=MINVAR(qeW, CLEVELS-1);
				int *curr_hist=args->hist+cdfstride*(nctx*kc+args->clevels*MINVAR(qeN, CLEVELS-1)+MINVAR(qeW, CLEVELS-1));
#define MIXCDF(X) curr_hist[X]
				den=MIXCDF(args->tlevels);
#endif
#ifndef ENABLE_ABAC2
				//if(den<=0)
				//{
				//	for(int k=0;k<8;++k)
				//	{
				//		for(int ks=0;ks<=args->tlevels;++ks)
				//			printf(" %d", curr_hist[k][ks]);
				//		printf("\n\n");
				//	}
				//}
				//if(ky==480&&kx==109&&kc==2)//
				//if(ky==147&&kx==317&&kc==1)//
				//if(ky==2562&&kx==1038&&kc==1)//
				//if(ky==208&&kx==6339&&kc==2)//
				//if(ky==208&&kx==6339&&kc==1)//
				//if(ky==208&&kx==6340&&kc==0)//
				//if(ky==261&&kx==6700&&kc==0)//
				//if(ky==5420&&kx==3267&&kc==1)//
				//if(ky==0&&kx==0&&kc==0)//
				//if(ky==1&&kx==41&&kc==2)//
				//if(ky==1&&kx==41&&kc==1)//
				//if(ky==5420&&kx==3267&&kc==1)//
				//if(ky==5420&&kx==3267&&kc==1)//
				//if(ky==128&&kx==128&&kc==0)//
				//if(ky==216&&kx==507&&kc==0)//
				//if(ky==0&&kx==86&&kc==0)//
				//if(ky==128&&kx==128&&kc==0)//
				//if(ky==1&&kx==32&&kc==0)//
				//	printf("");
//#ifdef AC_VALIDATE
//				//if(ky==5420&&kx==3072)//
//				if(args->blockidx==260)//
//					acval_disable=0;
//#endif
				if(args->fwd)
				{
					curr[kc2+0]=yuv[kc];
					curr[kc2+1]=error=yuv[kc]-pred;
#ifndef ENABLE_ABAC
					{
						int upred=half-abs(pred), aval=abs(error);
						if(aval<=upred)
						{
							sym=error;
#ifdef ENABLE_BIASCORR
							{
								int negmask=-((ibias_corr<0)&(sym!=-halfs[ch]));//sign is flipped if SSE correction was negative, to skew the histogram
								sym^=negmask;
								sym-=negmask;
							}
#endif
							sym=sym<<1^sym>>31;//pack sign
						}
						else
							sym=upred+aval;//error sign is known
					}
					quantize_pixel(sym, &token, &bypass, &nbits);
#ifdef _DEBUG
					if(token>=args->tlevels)
						LOG_ERROR("YXC %d %d %d  token %d/%d", ky, kx, kc, token, args->tlevels);
#endif
#endif
#ifdef ENABLE_ABAC
					tidx=1;
					token2=0;
					//error<<=32-ABAC_TOKEN_BITS;
					//error>>=32-ABAC_TOKEN_BITS;
					for(int kb=ABAC_TOKEN_BITS-1, bit=0;kb>=0;--kb)
					{
						//int *mixercell=curr_mixer+ABAC_NCTX*kb*2+bit;
						//if(ky==216&&kx==508&&kc==2&&kb==4)//
						//if(ky==0&&kx==86&&kc==2&&kb==2)//
						//if(ky==128&&kx==128&&kc==0&&kb==0)//
						//if(ky==1&&kx==32&&kc==1&&kb==0)//
						//if(ky==1&&kx==33&&kc==0&&kb==7)//
						//{
						//	g_compare=1;
						//	printf("");
						//}
#ifdef ABAC_SIMPLEOVERRIDE
						int p0=stats2[kc<<8|tidx];
#else
						int p0=abac_predict((const unsigned short**)curr_stats, tidx, curr_mixer, mixalphas, weights, logits);
#endif
						
						//if(ky==1&&kx==33&&kc==0&&kb==7)//
						//	g_compare=0;
						bit=error>>kb&1;
						ac3_enc_bin(&ec, bit, p0, ABAC_PROBBITS);
#ifdef ABAC_PROFILESIZE
						args->abac_csizes[kc*ABAC_TOKEN_BITS+kb]-=log2((double)(bit?(1<<ABAC_PROBBITS)-p0:p0)/(1<<ABAC_PROBBITS));
#endif
						token2|=bit<<kb;

						//if(ky==1&&kx==32)//
						//	printf("0x%04X %d\n", p0, bit);
						
#ifdef ABAC_SIMPLEOVERRIDE
						p0+=((!bit<<ABAC_PROBBITS)-p0)>>7;
						CLAMP2_32(stats2[kc<<8|tidx], p0, 1, 0xFFFF);
#else
						abac_update(weights, logits, curr_stats, tidx, curr_mixer, mixalphas, p0, bit);
#endif
						tidx+=tidx+bit;
						//if(tidx==7)//NPOT tree
						//	break;
					}
#else
					cdf=0;
					for(int ks=0;ks<token;++ks)
						cdf+=MIXCDF(ks);
					freq=MIXCDF(token);
					ac3_enc_update_NPOT(&ec, cdf, freq, den);
#endif
#ifndef ENABLE_ABAC
					if(nbits)
						ac3_enc_bypass(&ec, bypass, nbits);
#endif
				}
				else
				{
#ifdef ENABLE_ABAC
					tidx=1;
					token2=0;
					for(int kb=ABAC_TOKEN_BITS-1, bit=0;kb>=0;--kb)
					{
						//if(ky==1&&kx==33&&kc==0&&kb==7)//
						//	g_compare=2;
						//int *mixercell=curr_mixer+ABAC_NCTX*kb*2+bit;
#ifdef ABAC_SIMPLEOVERRIDE
						int p0=stats2[kc<<8|tidx];
#else
						int p0=abac_predict((const unsigned short**)curr_stats, tidx, curr_mixer, mixalphas, weights, logits);
#endif
						
						//if(ky==1&&kx==33&&kc==0&&kb==7)//
						//{
						//	g_compare=0;
						//	printf("");
						//}
						bit=ac3_dec_bin(&ec, p0, ABAC_PROBBITS);
						token2|=bit<<kb;

						//if(ky==1&&kx==32)//
						//	printf("0x%04X %d\n", p0, bit);
						
#ifdef ABAC_SIMPLEOVERRIDE
						p0+=((!bit<<ABAC_PROBBITS)-p0)>>7;
						CLAMP2_32(stats2[kc<<8|tidx], p0, 1, 0xFFFF);
#else
						abac_update(weights, logits, curr_stats, tidx, curr_mixer, mixalphas, p0, bit);
#endif
						tidx+=tidx+bit;
						//if(tidx==7)//NPOT tree
						//	break;
					}
					error=token2;
#else
					unsigned code=ac3_dec_getcdf_NPOT(&ec, den);
					cdf=0;
					token=0;
					for(;;)
					{
						unsigned cdf2;

						freq=MIXCDF(token);
						cdf2=cdf+freq;
						if(cdf2>code)
							break;
#ifdef _DEBUG
						if(token>=args->tlevels)
							LOG_ERROR("YXC %d %d %d  token %d/%d", ky, kx, kc, token, args->tlevels);
#endif
						cdf=cdf2;
						++token;
					}
					ac3_dec_update_NPOT(&ec, cdf, freq, den);
#endif
#ifndef ENABLE_ABAC
					sym=token;
					if(sym>=(1<<CONFIG_EXP))
					{
						int lsb, msb;

						sym-=1<<CONFIG_EXP;
						lsb=sym&((1<<CONFIG_LSB)-1);
						sym>>=CONFIG_LSB;
						msb=sym&((1<<CONFIG_MSB)-1);
						sym>>=CONFIG_MSB;
						nbits=sym+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB);
						bypass=ac3_dec_bypass(&ec, nbits);
						sym=1;
						sym<<=CONFIG_MSB;
						sym|=msb;
						sym<<=nbits;
						sym|=bypass;
						sym<<=CONFIG_LSB;
						sym|=lsb;
					}
					{
						int upred=half-abs(pred), negmask=0;
						if(sym<=(upred<<1))
						{
							error=sym>>1^-(sym&1);
#ifdef ENABLE_BIASCORR
							negmask=-((ibias_corr<0)&(error!=-half));
#endif
						}
						else
						{
							error=sym-upred;
							negmask=-(pred>0);
						}
						error^=negmask;
						error-=negmask;
					}
#endif
					curr[kc2+1]=error;
					error+=pred;
#ifdef ENABLE_ABAC
					error=error<<(32-ABAC_TOKEN_BITS)>>(32-ABAC_TOKEN_BITS);
					curr[kc2+1]=error-pred;
#endif
					curr[kc2+0]=yuv[kc]=error;
					//if((unsigned)(yuv[kc]+128)>255)
					//	LOG_ERROR("");
				}
				//X
#if 0
				{
					int c=curr[kc2+0];
					int kbest=0, kworst=0;
					int ebest=0, eworst=0;
					for(int kp=0;kp<PRED_COUNT;++kp)
					{
						int e=abs(c-preds[kp]);
						if(!kp||ebest>e)
							ebest=e, kbest=kp;
						if(!kp||eworst>e)
							eworst=e, kworst=kp;
					}
					{
						int inc=predmix[kc][kworst]>1&&predmix[kc][kbest]<255;
						predmix[kc][kbest]+=inc;
						predmix[kc][kworst]-=inc;
					}
				}
#endif
#if defined ENABLE_MIX4
				{
					int inc;
					inc=((1<<MIXBITS)-alphax)*((1<<MIXBITS)-alphay)>>(MIXBITS+MIXBITS-5); curr_hist[0][token]+=inc; curr_hist[0][args->tlevels]+=inc;
					inc=(             alphax)*((1<<MIXBITS)-alphay)>>(MIXBITS+MIXBITS-5); curr_hist[1][token]+=inc; curr_hist[1][args->tlevels]+=inc;
					inc=((1<<MIXBITS)-alphax)*(             alphay)>>(MIXBITS+MIXBITS-5); curr_hist[2][token]+=inc; curr_hist[2][args->tlevels]+=inc;
					inc=(             alphax)*(             alphay)>>(MIXBITS+MIXBITS-5); curr_hist[3][token]+=inc; curr_hist[3][args->tlevels]+=inc;
				}
				for(int kh=0;kh<4;++kh)
				{
					int *hist2=curr_hist[kh];
					if(hist2[args->tlevels]>=0xD400)//4296	6144	10752	65536
					{
						int sum=0;
						for(int ks=0;ks<args->tlevels;++ks)
							sum+=hist2[ks]=(hist2[ks]+1)>>1;
						hist2[args->tlevels]=sum;
					}
				}
#elif defined ENABLE_MIX8
				{
					int inc[8]=
					{
						(int)((long long)((1<<MIXBITSX)-alphas[0])*((1<<MIXBITSY)-alphas[1])*((1<<MIXBITSZ)-alphas[2])>>(MIXBITSX+MIXBITSY+MIXBITSZ-5)),
						(int)((long long)(              alphas[0])*((1<<MIXBITSY)-alphas[1])*((1<<MIXBITSZ)-alphas[2])>>(MIXBITSX+MIXBITSY+MIXBITSZ-5)),
						(int)((long long)((1<<MIXBITSX)-alphas[0])*(              alphas[1])*((1<<MIXBITSZ)-alphas[2])>>(MIXBITSX+MIXBITSY+MIXBITSZ-5)),
						(int)((long long)(              alphas[0])*(              alphas[1])*((1<<MIXBITSZ)-alphas[2])>>(MIXBITSX+MIXBITSY+MIXBITSZ-5)),
						(int)((long long)((1<<MIXBITSX)-alphas[0])*((1<<MIXBITSY)-alphas[1])*(              alphas[2])>>(MIXBITSX+MIXBITSY+MIXBITSZ-5)),
						(int)((long long)(              alphas[0])*((1<<MIXBITSY)-alphas[1])*(              alphas[2])>>(MIXBITSX+MIXBITSY+MIXBITSZ-5)),
						(int)((long long)((1<<MIXBITSX)-alphas[0])*(              alphas[1])*(              alphas[2])>>(MIXBITSX+MIXBITSY+MIXBITSZ-5)),
						(int)((long long)(              alphas[0])*(              alphas[1])*(              alphas[2])>>(MIXBITSX+MIXBITSY+MIXBITSZ-5)),
					};
					for(int kh=0;kh<8;++kh)
					{
						int *hist2=curr_hist[kh];
						//if(inc[kh]<0)
						//	LOG_ERROR("");
						hist2[token]+=inc[kh];
						hist2[args->tlevels]+=inc[kh];
						if(hist2[args->tlevels]>=(1<<11))//4296	6144	10752	65536	0xD400
						{
							int sum=0;
							for(int ks=0;ks<args->tlevels;++ks)
								sum+=hist2[ks]=(hist2[ks]+1)>>1;
							hist2[args->tlevels]=sum;
						}
					}
				}
#elif !defined ENABLE_ABAC
				++curr_hist[token];
				++curr_hist[args->tlevels];
				if(curr_hist[args->tlevels]>=6144)
				{
					int sum=0;
					for(int ks=0;ks<args->tlevels;++ks)
						sum+=curr_hist[ks]=(curr_hist[ks]+1)>>1;
					curr_hist[args->tlevels]=sum;
				}
#endif
#endif
				curr[kc2+0]-=offset;
			}
			if(!args->fwd)
			{
				switch(bestrct)
				{
				case RCT_R_G_B:
				case RCT_R_G_BG:
				case RCT_R_G_BR:
				case RCT_R_GR_BR:
				case RCT_R_GR_BG:
				case RCT_R_G_B2:
				case RCT_R_GR_B2:
					args->dst[idx+0]=yuv[0]+128;
					args->dst[idx+1]=yuv[1]+128;
					args->dst[idx+2]=yuv[2]+128;
					break;
				case RCT_G_B_RG:
				case RCT_G_B_RB:
				case RCT_G_BG_RG:
				case RCT_G_BG_RB:
				case RCT_G_B_R2:
				case RCT_G_BG_R2:
					args->dst[idx+1]=yuv[0]+128;
					args->dst[idx+2]=yuv[1]+128;
					args->dst[idx+0]=yuv[2]+128;
					break;
				case RCT_B_R_GR:
				case RCT_B_R_GB:
				case RCT_B_RB_GB:
				case RCT_B_RB_GR:
				case RCT_B_RB_G2:
					args->dst[idx+2]=yuv[0]+128;
					args->dst[idx+0]=yuv[1]+128;
					args->dst[idx+1]=yuv[2]+128;
					break;
				case RCT_G_RG_BR:
				case RCT_G_RG_B2:
					args->dst[idx+1]=yuv[0]+128;
					args->dst[idx+0]=yuv[1]+128;
					args->dst[idx+2]=yuv[2]+128;
					break;
				case RCT_B_GB_RG:
				case RCT_B_GB_R2:
					args->dst[idx+2]=yuv[0]+128;
					args->dst[idx+1]=yuv[1]+128;
					args->dst[idx+0]=yuv[2]+128;
					break;
				case RCT_R_BR_GB:
				case RCT_R_B_G2:
				case RCT_R_BR_G2:
					args->dst[idx+0]=yuv[0]+128;
					args->dst[idx+2]=yuv[1]+128;
					args->dst[idx+1]=yuv[2]+128;
					break;
				}
#ifdef ENABLE_GUIDE
				if(args->test&&memcmp(args->dst+idx, args->src+idx, sizeof(char)*nch))
				{
					unsigned char orig[4]={0};
					memcpy(orig, args->src+idx, nch*sizeof(char));
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	if(args->fwd)
		ac3_enc_flush(&ec);
}
int c01_codec(const char *srcfn, const char *dstfn)
{
	const int nch=3, depth=8;
	double t0;
	ArrayHandle src, dst;
	int headersize, printed;
	int iw, ih;
	const unsigned char *image, *imageend;
	unsigned char *image2;
	CodecID codec;
	int ncores=query_cpu_cores();
	int xblocks, yblocks, nblocks, nthreads, coffset;
	ptrdiff_t start, memusage, argssize;
	ThreadArgs *args;
	int test, fwd;
#ifdef ENABLE_MIX4
	int clevels;
#endif
	int tlevels, histsize, statssize;
	double esize;
	int usize;
#ifdef ABAC_PROFILESIZE
	double abac_csizes[ABAC_TOKEN_BITS*3]={0};
#endif
#ifdef ENABLE_ABAC
	stretch(0);//for thread safety
#endif
	
	t0=time_sec();
	src=load_file(srcfn, 1, 3, 1);
	headersize=header_read(src->data, (int)src->count, &iw, &ih, &codec);
	image=src->data+headersize;
	imageend=src->data+src->count;
	if(codec==CODEC_INVALID||codec==CODEC_PGM)
	{
		LOG_ERROR("Unsupported codec %d.\n", codec);
		array_free(&src);
		return 1;
	}
	else if(codec==CODEC_C01&&!dstfn)
	{
		LOG_ERROR(
			"Test mode expects PPM source.\n"
			"Decode mode expects destination filename."
		);
		return 1;
	}
	test=!dstfn;
	fwd=codec==CODEC_PPM;
	
	if(test)
		printf("%s \"%s\"  WH %d*%d\n", CODECNAME, srcfn, iw, ih);
	usize=iw*ih*3;
	xblocks=(iw+BLOCKSIZE-1)/BLOCKSIZE;
	yblocks=(ih+BLOCKSIZE-1)/BLOCKSIZE;
	nblocks=xblocks*yblocks, nthreads=MINVAR(nblocks, ncores);
	coffset=(int)sizeof(int)*nblocks;
	start=0;
	memusage=0;
	argssize=nthreads*sizeof(ThreadArgs);
	args=(ThreadArgs*)malloc(argssize);
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	esize=0;
	memusage+=argssize;
	memset(args, 0, argssize);
	if(fwd)
	{
		dst=0;
		printed=snprintf(g_buf, G_BUF_SIZE-1, "C01\n%d %d\n", iw, ih);
		array_append(&dst, g_buf, 1, printed, 1, 0, 0);
		start=array_append(&dst, 0, 1, coffset, 1, 0, 0);
		
		image2=0;
	}
	else//integrity check
	{
		dst=0;
		printed=snprintf(g_buf, G_BUF_SIZE-1, "P6\n%d %d\n255\n", iw, ih);
		array_append(&dst, g_buf, 1, printed, 1, 0, 0);
		array_append(&dst, 0, 1, usize, 1, 0, 0);

		//printed=0;
		start=coffset;
		for(int kt=0;kt<nblocks;++kt)
		{
			int size=0;
			memcpy(&size, image+sizeof(int)*kt, sizeof(int));
			start+=size;
		}
		if(image+start!=imageend)
			LOG_ERROR("Corrupt file");
		start=coffset;

		image2=(unsigned char*)malloc(usize);
		if(!image2)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		memset(image2, 0, usize);
	}
	{
		int nlevels=256;
		int token=0, bypass=0, nbits=0;

		quantize_pixel(nlevels, &token, &bypass, &nbits);
		tlevels=token+1;
#if defined ENABLE_CALICCTX || defined ENABLE_HASH
		statssize=(tlevels+1)*nch*(int)sizeof(int[CLEVELS]);
#elif defined ENABLE_MIX8
		statssize=(tlevels+1)*nch*(int)sizeof(int[CXLEVELS*CYLEVELS*CZLEVELS]);
#elif defined ENABLE_ABAC2
		statssize=0;
#elif defined ENABLE_ABAC
		statssize=nch*(int)sizeof(short[ABAC_NCTRS*ABAC_NCTX*ABAC_CLEVELS*ABAC_TREESIZE]);
#else
		clevels=CLEVELS;
		statssize=clevels*clevels*(tlevels+1)*nch*(int)sizeof(int);
#endif
		histsize=(int)sizeof(int[OCH_COUNT*PRED_COUNT<<8]);
		if(histsize<statssize)
			histsize=statssize;
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=image;
		arg->dst=fwd?0:dst->data+printed;
		arg->iw=iw;
		arg->ih=ih;
		arg->bufsize=sizeof(short[4*OCH_COUNT*2])*(BLOCKSIZE+16LL);//4 padded rows * OCH_COUNT * {pixels, wg_errors}
		arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));

		arg->histsize=histsize;
		arg->hist=(int*)malloc(histsize);
#ifdef ENABLE_ABAC
		arg->statssize=statssize;
#endif
		//arg->stats=(unsigned*)malloc(statssize);
		if(!arg->pixels||!arg->hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=arg->bufsize;
		memusage+=histsize;
		//memusage+=statssize;
		
		arg->tlevels=tlevels;
#ifdef ENABLE_MIX4
		arg->clevels=clevels;
#endif
		arg->fwd=fwd;
		arg->test=test;
#ifdef DISABLE_MT
		arg->loud=test&&nblocks<MAXPRINTEDBLOCKS;
#else
		arg->loud=0;
#endif
	}
	for(int k2=0;k2<=test;++k2)
	{
		for(int kt=0;kt<nblocks;kt+=nthreads)
		{
			int nthreads2=MINVAR(kt+nthreads, nblocks)-kt;
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				int kx, ky;

				arg->blockidx=kx=kt+kt2;
				ky=kx/xblocks;
				kx%=xblocks;
				arg->x1=BLOCKSIZE*kx;
				arg->y1=BLOCKSIZE*ky;
				arg->x2=MINVAR(arg->x1+BLOCKSIZE, iw);
				arg->y2=MINVAR(arg->y1+BLOCKSIZE, ih);
				if(!fwd)
				{
					int size=0;
					memcpy(&size, image+sizeof(int)*((ptrdiff_t)kt+kt2), sizeof(int));
					arg->decstart=image+start;
					start+=size;
					arg->decend=image+start;
				}
			}
#ifdef DISABLE_MT
			for(int k=0;k<nthreads2;++k)
				block_thread(args+k);
#else
			void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
			mt_finish(ctx);
#endif
			if(fwd)
			{
				for(int kt2=0;kt2<nthreads2;++kt2)
				{
					ThreadArgs *arg=args+kt2;
					if(test)
					{
						int blocksize=((arg->x2-arg->x1)*(arg->y2-arg->y1)*nch*depth+7)>>3;
						int kx, ky;

						kx=kt+kt2;
						ky=kx/xblocks;
						kx%=xblocks;
						if(nblocks<MAXPRINTEDBLOCKS)
						{
							//if(!(kt+kt2))
							//	printf("block,  nrows,  usize,     best  ->  actual,  (actual-best)\n");
							printf(
								"block %4d/%4d  XY %3d %3d  %4d*%4d:  %8d->%16lf->%8zd bytes (%+10.2lf)  %10.6lf%%  CR %10lf  %s %s %s %s\n",
								kt+kt2+1, nblocks,
								kx, ky,
								arg->y2-arg->y1,
								arg->x2-arg->x1,
								blocksize,
								arg->bestsize,
								arg->list.nobj,
								arg->list.nobj-arg->bestsize,
								100.*arg->list.nobj/blocksize,
								(double)blocksize/arg->list.nobj,
								rct_names[arg->bestrct],
								pred_names[arg->predidx[0]],
								pred_names[arg->predidx[1]],
								pred_names[arg->predidx[2]]
							);
						}
						esize+=arg->bestsize;
#ifdef ABAC_PROFILESIZE
						for(int k=0;k<ABAC_TOKEN_BITS*3;++k)
							abac_csizes[k]+=arg->abac_csizes[k];
#endif
					}
					memcpy(dst->data+start+sizeof(int)*((ptrdiff_t)kt+kt2), &arg->list.nobj, sizeof(int));
					dlist_appendtoarray(&arg->list, &dst);
					dlist_clear(&arg->list);
				}
			}
		}
		if(test)
		{
			ptrdiff_t usize=((ptrdiff_t)iw*ih*nch*depth+7)>>3;
			ptrdiff_t csize=dst->count;
			t0=time_sec()-t0;
			if(fwd)
			{
#ifdef ABAC_PROFILESIZE
				double upsize=(double)iw*ih/8;
				for(int kb=0;kb<ABAC_TOKEN_BITS*3;++kb)
				{
					double cpsize=abac_csizes[kb]/8;
					printf("C%d B%d  %10.2lf/%10.2lf  %8.4lf%%  CR %10.6lf\n",
						kb/ABAC_TOKEN_BITS,
						kb%ABAC_TOKEN_BITS,
						cpsize,
						upsize,
						100.*cpsize/upsize,
						upsize/cpsize
					);
				}
#endif
				printf("Best %15.2lf (%+13.2lf) bytes\n", esize, csize-esize);
				printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
				printf("Mem usage: ");
				print_size((double)memusage, 8, 4, 0, 0);
				printf("\n");
			}
			printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
			if(!fwd)
				compare_bufs_8(image2, src->data+headersize, iw, ih, nch, nch, CODECNAME, 0, 1);
		}
		if(!k2&&test)
		{
			int usize=iw*ih*3;
			fwd=0;
			image2=(unsigned char*)malloc(usize);
			if(!image2)
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			memset(image2, 0, usize);
			for(int kt=0;kt<nthreads;++kt)
			{
				ThreadArgs *arg=args+kt;
				arg->dst=image2;
				arg->fwd=0;
			}
			image=dst->data+printed;
			start=coffset;
		}
		t0=time_sec();
	}
	if(!test)
		save_file(dstfn, dst->data, dst->count, 1);
	if(image2)
		free(image2);
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		_mm_free(arg->pixels);
		free(arg->hist);
		//free(arg->stats);
	}
	free(args);
	array_free(&src);
	array_free(&dst);
	//if(test)
	//	pause();
	return 0;
}