static const char file[]=__FILE__;
#include"ppm.h"
#include"util.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"


//	#define ENABLE_GUIDE
#ifndef DISABLE_MT
	#define ENABLE_MT
#endif


//select one:
//	#define USE_GRCODER
	#define ENABLE_MIX4//good
//	#define ENABLE_MIX8//slightly worse
//	#define ENABLE_ABAC2
//	#define ENABLE_ABAC//bad
//	#define ENABLE_CALICCTX//bad
//	#define ENABLE_HASH//bad

//	#define MIXERLAYERS
//	#define PREDICT_SIGN
//	#define ENABLE_OLS
	#define DISABLE_PREDSEL
	#define ENABLE_WG

#define CODECNAME "C02"
#include"entropy.h"

#define BLOCKSIZE 512
#define MAXPRINTEDBLOCKS 200
#ifdef MIXERLAYERS
#define MIXERINIT 0x100000
#define MIXERCLAMP 0x1000000
#else
#define MIXERINIT 0x3000
#define MIXERCLAMP 0xC000
#endif
#define EBITS 0		//up to 8
#define NBITS 0
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
#define A2_NCTX 12
#define A2_CTXBITS 8
#define A2_NCTX2 3
#define A2_CTXBITS2 6
#ifdef MIXERLAYERS
#define A2_CTXGROUPSIZE 2
#define A2_NCTXGROUPS (A2_NCTX/A2_CTXGROUPSIZE)
#endif
#define A2_SSEBITS 6
#define A2_SSECTR 16
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

#ifdef ENABLE_OLS
#define OLS_REACH 1
#define OLS_NPARAMS0 (2*(OLS_REACH+1)*OLS_REACH*3+0)
#define OLS_NPARAMS1 (2*(OLS_REACH+1)*OLS_REACH*3+1)
#define OLS_NPARAMS2 (2*(OLS_REACH+1)*OLS_REACH*3+2)
#define OLS_NPARAMS (2*(OLS_REACH+1)*OLS_REACH*3*3+3)
#endif

#define PERMUTATIONLIST\
	PERM(RGB, 0, 1, 2)\
	PERM(GBR, 1, 2, 0)\
	PERM(BRG, 2, 0, 1)\
	PERM(BGR, 2, 1, 0)\
	PERM(RBG, 0, 2, 1)\
	PERM(GRB, 1, 0, 2)
typedef enum _CRCTPermutation
{
#define PERM(L, A, B, C) PERM_##L,
	PERMUTATIONLIST
#undef  PERM
	PERM_COUNT,
} CRCTPermutation;
static int permutations[6][3]=
{
#define PERM(L, A, B, C) {A, B, C},
	PERMUTATIONLIST
#undef  PERM
};
const char *permnames[]=
{
#define PERM(L, A, B, C) #L,
	PERMUTATIONLIST
#undef  PERM
};

#if 1
#define OCHLIST\
	OCH(Rgb00)\
	OCH(Rgb0G)\
	OCH(Rgb1F)\
	OCH(Rgb2E)\
	OCH(Rgb3D)\
	OCH(Rgb4C)\
	OCH(Rgb5B)\
	OCH(Rgb6A)\
	OCH(Rgb79)\
	OCH(Rgb88)\
	OCH(Rgb97)\
	OCH(RgbA6)\
	OCH(RgbB5)\
	OCH(RgbC4)\
	OCH(RgbD3)\
	OCH(RgbE2)\
	OCH(RgbF1)\
	OCH(RgbG0)\
	OCH(Gbr00)\
	OCH(Gbr0G)\
	OCH(Gbr1F)\
	OCH(Gbr2E)\
	OCH(Gbr3D)\
	OCH(Gbr4C)\
	OCH(Gbr5B)\
	OCH(Gbr6A)\
	OCH(Gbr79)\
	OCH(Gbr88)\
	OCH(Gbr97)\
	OCH(GbrA6)\
	OCH(GbrB5)\
	OCH(GbrC4)\
	OCH(GbrD3)\
	OCH(GbrE2)\
	OCH(GbrF1)\
	OCH(GbrG0)\
	OCH(Brg00)\
	OCH(Brg0G)\
	OCH(Brg1F)\
	OCH(Brg2E)\
	OCH(Brg3D)\
	OCH(Brg4C)\
	OCH(Brg5B)\
	OCH(Brg6A)\
	OCH(Brg79)\
	OCH(Brg88)\
	OCH(Brg97)\
	OCH(BrgA6)\
	OCH(BrgB5)\
	OCH(BrgC4)\
	OCH(BrgD3)\
	OCH(BrgE2)\
	OCH(BrgF1)\
	OCH(BrgG0)
typedef enum _OCHIndex
{
#define OCH(X) OCH_##X,
	OCHLIST
#undef  OCH
	OCH_COUNT,
} OCHIndex;
const char *ochnames[]=
{
#define OCH(X) #X,
	OCHLIST
#undef  OCH
};

const int lumapool[]=
{
	OCH_Rgb00,
	OCH_Gbr00,
	OCH_Brg00,
};
//const int chroma1pool[]=
//{//			R	G	B
//	OCH_Rgb00,//	X
//	OCH_Rgb0G,//	X	X
//	OCH_RgbG0,//	X		X
//	OCH_Gbr00,//		X
//	OCH_Gbr0G,//		X	X
//	OCH_GbrG0,//	X	X
//	OCH_Brg00,//			X
//	OCH_Brg0G,//	X		X
//	OCH_BrgG0,//		X	X
//};//chroma2 unfiltered pool includes all OCH_COUNT=54 output channels
const int chroma1pool[3][4]=
{
	{OCH_Gbr00, OCH_Gbr0G, OCH_Brg00, OCH_BrgG0},//R
	{OCH_Rgb00, OCH_RgbG0, OCH_Brg00, OCH_Brg0G},//G
	{OCH_Rgb00, OCH_Rgb0G, OCH_Gbr00, OCH_GbrG0},//B
};
#endif
#if 0
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

#ifdef DISABLE_PREDSEL
#define PREDLIST\
	PRED(CG)
#else
#define PREDLIST\
	PRED(W)\
	PRED(CG)\
	PRED(AV5)\
	PRED(AV9)\
	PRED(AV12)
#endif
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
#endif

typedef enum _NBIndex
{
	NB_NW,		NB_N,		NB_NE,
	NB_W,		NB_curr,
} NBIndex;

#if 0
#define PREDLIST\
	PRED(N)\
	PRED(W)\
	PRED(CG)\
	PRED(AV5)\
	PRED(AV9)\
	PRED(AV12)\
	PRED(WG)
#endif
#define PREDLIST\
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
typedef enum _N2Index
{
	N2_NNWW,	N2_NNW,		N2_NN,		N2_NNE,		N2_NNEE,
	N2_NWW,		N2_NW,		N2_N,		N2_NE,		N2_NEE,
	N2_WW,		N2_W,		N2_curr,
} N2Index;
static const short av12_icoeffs[12]=
{
	 0x04,	 0x03,	-0x1F,	-0x26,	 0x00,
	 0x07,	-0x9E,	 0xDB,	 0x1E,	 0x13,
	-0x2A,	 0xF3,
};


//WG
#ifdef ENABLE_WG
#define WG_DECAY_NUM	493
#define WG_DECAY_SH	9

#define WG_NPREDS	16
#define WG_PREDLIST\
	WG_PRED(100, spred)\
	WG_PRED(50, wgrad)\
	WG_PRED(50, 3*N-W-NE)\
	WG_PRED(50, cgrad)\
	WG_PRED(132, N+W-NW)\
	WG_PRED(100, N+W-NW+((eN+eW-eNW+16)>>5))\
	WG_PRED(150, N+eN)\
	WG_PRED(100, N)\
	WG_PRED(150, W+eW)\
	WG_PRED(100, W)\
	WG_PRED(165, W+NE-N)\
	WG_PRED(220, N+NE-NNE)\
	WG_PRED(25, N+NE-NNE+((eN+eNE-eNNE)>>2))\
	WG_PRED(165, 3*(N-NN)+NNN)\
	WG_PRED(165, 3*(W-WW)+WWW)\
	WG_PRED(150, (W+NEEE)/2)
//	WG_PRED(30, 2*NE-W)
static void wg_init(int *weights)
{
	int j=0;
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
	WG_PREDLIST
#undef  WG_PRED
}
static int wg_predict(
	const int *weights,
	short **rows, const int stride, int kc2,
	int spred, int ols, const long long *perrors, int *preds
)
{
	double pred2=0, wsum=0;
	int j=0, pred;
	short
		NNNWWWW	=rows[3][kc2-4*stride+0],
		NNWWWW	=rows[2][kc2-4*stride+0],
		NNN	=rows[3][kc2+0*stride+0],
		NN	=rows[2][kc2+0*stride+0],
		NNE	=rows[2][kc2+1*stride+0],
		NNEEE	=rows[2][kc2+3*stride+0],
		NW	=rows[1][kc2-1*stride+0],
		N	=rows[1][kc2+0*stride+0],
		NE	=rows[1][kc2+1*stride+0],
		NEE	=rows[1][kc2+2*stride+0],
		NEEE	=rows[1][kc2+3*stride+0],
		WWW	=rows[0][kc2-3*stride+0],
		WW	=rows[0][kc2-2*stride+0],
		W	=rows[0][kc2-1*stride+0],
		eNNE	=rows[2][kc2+1*stride+1],
		eNW	=rows[1][kc2-1*stride+1],
		eN	=rows[1][kc2+0*stride+1],
		eNE	=rows[1][kc2+1*stride+1],
		eW	=rows[0][kc2-1*stride+1];

	int gy=abs(eN)+1, gx=abs(eW)+1;
	//int
	//	gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1,
	//	gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1;
	int wgrad=(N*gy+W*gx)/(gy+gx);
	int cgrad;
	MEDIAN3_32(cgrad, N, W, N+W-NW);
	(void)NNNWWWW;
	(void)NNWWWW;
	(void)NNEEE;
	(void)NEE;
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST
#undef  WG_PRED
	
	for(int k=0;k<WG_NPREDS;++k)
	{
		double weight=(double)weights[k]/(perrors[k]+1);
		pred2+=weight*preds[k];
		wsum+=weight;
	}
	pred2/=wsum;
#ifdef _MSC_VER
	pred=_cvt_dtoi_fast(pred2);
#else
	pred=(int)pred2;
#endif
	CLAMP3_32(pred, pred, N, W, NE);
	return pred;
}
static void wg_update(int curr, const int *preds, long long *perrors, int *weights)
{
	for(int k=0;k<WG_NPREDS;++k)
	{
		int e2=abs(curr-preds[k]);
		perrors[k]=(perrors[k]+e2)*WG_DECAY_NUM>>WG_DECAY_SH;
	}
}
#endif


//OLS
#if 0
/*
example with size 2:
matrix buffer before call:
	Fwd  Fwd  1  0
	Fwd  Fwd  0  1
matrix buffer after call, on success:
	1  0  Inv  Inv
	0  1  Inv  Inv
temprow buffer [size*2]
*/
static int invert_matrix(double *matrix, int size, double *temprow)
{
	int success;

	success=1;
	for(int it=0;it<size;++it)
	{
		int kp;
		double pivot;

		kp=it;
		for(;kp<size;++kp)
		{
			if(fabs(matrix[((size_t)size<<1)*kp+it])>1e-6)
				break;
		}
		if(kp==size)
		{
			success=0;
			break;
		}
		if(kp!=it)
		{
			memcpy(temprow, matrix+((size_t)size<<1)*it, size*sizeof(double[2]));
			memcpy(matrix+((size_t)size<<1)*it, matrix+((size_t)size<<1)*kp, size*sizeof(double[2]));
			memcpy(matrix+((size_t)size<<1)*kp, temprow, size*sizeof(double[2]));
		}
		pivot=matrix[((size_t)size<<1)*it+it];
		for(int kx=it;kx<(size<<1);++kx)
			matrix[((size_t)size<<1)*it+kx]/=pivot;
		for(int ky=0;ky<size;++ky)
		{
			double factor;

			if(ky==it)
				continue;
			factor=matrix[((size_t)size<<1)*ky+it];
			if(fabs(factor)>1e-6)
			{
				for(int kx=it;kx<(size<<1);++kx)
					matrix[((size_t)size<<1)*ky+kx]-=matrix[((size_t)size<<1)*it+kx]*factor;
			}
		}
	}
	//_freea(temp);
	return success;
}
#endif


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
#if 0
static void xorshift64(unsigned long long *state)
{
	unsigned long long x=*state;
	x^=x<<13;
	x^=x>>7;
	x^=x<<17;
	*state=x;
}
#endif
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
	int bufsize;
	//int histsize;
	short *pixels;
	//int *hist;

	BList list;
	const unsigned char *decstart, *decend;
	
	int hist[54<<8];
//#ifdef ENABLE_OLS
//	long long ols_imat[3][OLS_NPARAMS2*OLS_NPARAMS2];
//	long long ols_ivec[3][OLS_NPARAMS2];
//	double ols_mat[3][OLS_NPARAMS2*OLS_NPARAMS2<<1];
//	double ols_vec[3][OLS_NPARAMS2];
//	int ols_params[3][OLS_NPARAMS2];
//#endif

//#ifdef ENABLE_MIX4
	int clevels;
//#endif
	int tlevels;
//#ifdef ENABLE_ABAC2
//	unsigned short stats2[3][A2_NCTX2][1<<A2_CTXBITS2][1<<A2_CTXBITS2][1<<8];
//	unsigned stats[A2_NCTX*3<<A2_CTXBITS<<8];
//	unsigned long long stats0[3][1<<8];
//#ifdef PREDICT_SIGN
//	long long sse[9*3<<A2_SSEBITS];
//#else
//	long long sse[8*3<<A2_SSEBITS];
//#endif
//#elif defined ENABLE_ABAC
//	int statssize;
//#ifdef ABAC_PROFILESIZE
//	double abac_csizes[ABAC_TOKEN_BITS*3];
//#endif
////	int mixer[ABAC_NCTX*3<<ABAC_TOKEN_BITS];
//	int mixer[ABAC_NCTRS*ABAC_NCTX*3*((1<<(8-ABAC_MIXBITS))+1)*(1<<(8-ABAC_MIXBITS)|1)];
////	int mixer[ABAC_NCTX*ABAC_TOKEN_BITS*3*2];
////	int mixer[ABAC_NCTX*ABAC_TOKEN_BITS*3];
//#endif

	//aux
	int blockidx;
	double bestsize;
	int permutation, helper1, alpha1, alpha2;
	int predsel[3];
	//int bestrct, predidx[3];
} ThreadArgs;
static void block_thread(void *param)
{
	const int nch=3;
	ThreadArgs *args=(ThreadArgs*)param;
#ifdef USE_GRCODER
	GolombRiceCoder ec;
#else
	AC3 ec;
#endif
	const unsigned char *image=args->fwd?args->src:args->dst;
	//unsigned char bestrct=0, combination[6]={0}, predidx[4]={0};
#ifdef ENABLE_MIX4
	int cdfstride=args->tlevels+1;
	int nctx=args->clevels*args->clevels;
	int chsize=nctx*cdfstride;
#elif defined ENABLE_ABAC
	unsigned short *stats=(unsigned short*)args->hist;
#endif
	
	int permutation=0, helper1=0, alpha1=0, alpha2;
	int rgbidx[3]={0};//derived from permutation
	int predsel[3]={0};
#ifdef ENABLE_OLS
	int ols_success[3]={0};
#endif
	if(args->fwd)
	{
		int ystride=args->iw*3;
		int res=(args->x2-args->x1-2)/5*5*(args->y2-args->y1-1);
		__m256i ramp[]=
		{
			_mm256_set1_epi16(0),
			_mm256_set1_epi16(1),
			_mm256_set1_epi16(2),
			_mm256_set1_epi16(3),
			_mm256_set1_epi16(4),
			_mm256_set1_epi16(5),
			_mm256_set1_epi16(6),
			_mm256_set1_epi16(7),
			_mm256_set1_epi16(8),
			_mm256_set1_epi16(9),
			_mm256_set1_epi16(10),
			_mm256_set1_epi16(11),
			_mm256_set1_epi16(12),
			_mm256_set1_epi16(13),
			_mm256_set1_epi16(14),
			_mm256_set1_epi16(15),
			_mm256_set1_epi16(16),
		};
#ifdef ENABLE_OLS
		long long ols_ivec[OLS_NPARAMS0+3]={0};
		memset(args->ols_imat, 0, sizeof(args->ols_imat));
		memset(args->ols_ivec, 0, sizeof(args->ols_ivec));
#endif
		memset(args->hist, 0, sizeof(args->hist));
		for(int ky=args->y1+1;ky<args->y2;++ky)//analysis loop
		{
			int kx=args->x1+1;
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
					pred, vmin, vmax,
					rgb[5],
					gbr[5],
					brg[5];
				{
					__m128i nb8[]=//8-bit
					{
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride-1*3+0)), half8),//NW
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride+0*3+0)), half8),//N
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr-1*ystride+1*3+0)), half8),//NE
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr+0*ystride-1*3+0)), half8),//W
						_mm_xor_si128(_mm_loadu_si128((__m128i*)(ptr+0*ystride+0*3+0)), half8),//curr
					};
					for(int k=0;k<(int)_countof(rgb);++k)
					{
						__m128i temp;
						rgb[k]=_mm256_cvtepi8_epi16(nb8[k]);
						temp=_mm_shuffle_epi8(nb8[k], shuf);
						gbr[k]=_mm256_cvtepi8_epi16(temp);
						brg[k]=_mm256_cvtepi8_epi16(_mm_shuffle_epi8(temp, shuf));
					}
				}
#define UPDATE(IDX0, IDX1, IDX2, IDX3, IDX4, IDX5, IDX6, IDX7, IDX8, IDX9, IDXA, IDXB, IDXC, IDXD, IDXE)\
	do\
	{\
		pred=_mm256_slli_epi16(pred, 8);\
		pred=_mm256_srai_epi16(pred, 8);\
		pred=_mm256_sub_epi16(pred, amin);\
		_mm256_store_si256((__m256i*)result, pred);\
		++args->hist[(IDX0)<<8|result[0x0]];\
		++args->hist[(IDX1)<<8|result[0x1]];\
		++args->hist[(IDX2)<<8|result[0x2]];\
		++args->hist[(IDX3)<<8|result[0x3]];\
		++args->hist[(IDX4)<<8|result[0x4]];\
		++args->hist[(IDX5)<<8|result[0x5]];\
		++args->hist[(IDX6)<<8|result[0x6]];\
		++args->hist[(IDX7)<<8|result[0x7]];\
		++args->hist[(IDX8)<<8|result[0x8]];\
		++args->hist[(IDX9)<<8|result[0x9]];\
		++args->hist[(IDXA)<<8|result[0xA]];\
		++args->hist[(IDXB)<<8|result[0xB]];\
		++args->hist[(IDXC)<<8|result[0xC]];\
		++args->hist[(IDXD)<<8|result[0xD]];\
		++args->hist[(IDXE)<<8|result[0xE]];\
	}while(0)
				//CG
				vmin=_mm256_min_epi16(rgb[NB_N], rgb[NB_W]);
				vmax=_mm256_max_epi16(rgb[NB_N], rgb[NB_W]);
				pred=_mm256_sub_epi16(_mm256_add_epi16(rgb[NB_N], rgb[NB_W]), rgb[NB_NW]);
				pred=_mm256_max_epi16(pred, vmin);
				pred=_mm256_min_epi16(pred, vmax);

				pred=_mm256_sub_epi16(rgb[NB_curr], pred);
				UPDATE(
					OCH_Rgb00, OCH_Gbr00, OCH_Brg00,
					OCH_Rgb00, OCH_Gbr00, OCH_Brg00,
					OCH_Rgb00, OCH_Gbr00, OCH_Brg00,
					OCH_Rgb00, OCH_Gbr00, OCH_Brg00,
					OCH_Rgb00, OCH_Gbr00, OCH_Brg00
				);
				for(int k=0;k<=16;++k)
				{
					__m256i helpers[]=
					{
						_mm256_srai_epi16(_mm256_add_epi16(_mm256_mullo_epi16(gbr[0], ramp[k]), _mm256_mullo_epi16(brg[0], ramp[16-k])), 4),//NW
						_mm256_srai_epi16(_mm256_add_epi16(_mm256_mullo_epi16(gbr[1], ramp[k]), _mm256_mullo_epi16(brg[1], ramp[16-k])), 4),//N
						_mm256_srai_epi16(_mm256_add_epi16(_mm256_mullo_epi16(gbr[2], ramp[k]), _mm256_mullo_epi16(brg[2], ramp[16-k])), 4),//NE (unused!)
						_mm256_srai_epi16(_mm256_add_epi16(_mm256_mullo_epi16(gbr[3], ramp[k]), _mm256_mullo_epi16(brg[3], ramp[16-k])), 4),//W
						_mm256_srai_epi16(_mm256_add_epi16(_mm256_mullo_epi16(gbr[4], ramp[k]), _mm256_mullo_epi16(brg[4], ramp[16-k])), 4),//curr
					};
					__m256i nb2[]=
					{
						_mm256_sub_epi16(rgb[0], helpers[0]),//NW
						_mm256_sub_epi16(rgb[1], helpers[1]),//N
						_mm256_sub_epi16(rgb[2], helpers[2]),//NE (unused!)
						_mm256_sub_epi16(rgb[3], helpers[3]),//W
					};
					vmin=_mm256_min_epi16(nb2[NB_N], nb2[NB_W]);
					vmax=_mm256_max_epi16(nb2[NB_N], nb2[NB_W]);
					pred=_mm256_sub_epi16(_mm256_add_epi16(nb2[NB_N], nb2[NB_W]), nb2[NB_NW]);
					pred=_mm256_max_epi16(pred, vmin);
					pred=_mm256_min_epi16(pred, vmax);

					pred=_mm256_add_epi16(pred, helpers[NB_curr]);
					pred=_mm256_max_epi16(pred, amin);
					pred=_mm256_min_epi16(pred, amax);
					pred=_mm256_sub_epi16(rgb[NB_curr], pred);
					UPDATE(
						OCH_Rgb0G+k, OCH_Gbr0G+k, OCH_Brg0G+k,
						OCH_Rgb0G+k, OCH_Gbr0G+k, OCH_Brg0G+k,
						OCH_Rgb0G+k, OCH_Gbr0G+k, OCH_Brg0G+k,
						OCH_Rgb0G+k, OCH_Gbr0G+k, OCH_Brg0G+k,
						OCH_Rgb0G+k, OCH_Gbr0G+k, OCH_Brg0G+k
					);
				}
#ifdef ENABLE_OLS
				ALIGN(32) short nbNW[16], nbN[16], nbNE[16], nbW[16], nbcurr[16];
				_mm256_store_si256((__m256i*)nbNW	, rgb[0]);
				_mm256_store_si256((__m256i*)nbN	, rgb[1]);
				_mm256_store_si256((__m256i*)nbNE	, rgb[2]);
				_mm256_store_si256((__m256i*)nbW	, rgb[3]);
				_mm256_store_si256((__m256i*)nbcurr	, rgb[4]);
				for(int kx2=0;kx2<5;++kx2)
				{
					int j=0;
					ols_ivec[j++]=nbNW	[kx2*3+0];
					ols_ivec[j++]=nbNW	[kx2*3+1];
					ols_ivec[j++]=nbNW	[kx2*3+2];
					ols_ivec[j++]=nbN	[kx2*3+0];
					ols_ivec[j++]=nbN	[kx2*3+1];
					ols_ivec[j++]=nbN	[kx2*3+2];
					ols_ivec[j++]=nbNE	[kx2*3+0];
					ols_ivec[j++]=nbNE	[kx2*3+1];
					ols_ivec[j++]=nbNE	[kx2*3+2];
					ols_ivec[j++]=nbW	[kx2*3+0];
					ols_ivec[j++]=nbW	[kx2*3+1];
					ols_ivec[j++]=nbW	[kx2*3+2];
					ols_ivec[j++]=nbcurr	[kx2*3+0];
					ols_ivec[j++]=nbcurr	[kx2*3+1];
					ols_ivec[j++]=nbcurr	[kx2*3+2];
					for(int kc=0;kc<3;++kc)
					{
						int nparams=OLS_NPARAMS0+kc;
						for(int y=0;y<nparams;++y)
						{
							int vy=ols_ivec[y];
							for(int x=0;x<nparams;++x)
							{
								int vx=ols_ivec[x];
								args->ols_imat[kc][nparams*y+x]+=(long long)vy*vx;
							}
						}
						for(int k=0;k<nparams;++k)
							args->ols_ivec[kc][k]+=ols_ivec[k]*ols_ivec[nparams];
					}
				}
#endif
			}
		}
		double csizes[OCH_COUNT]={0};
		for(int kc=0;kc<OCH_COUNT;++kc)
		{
			int *curr_hist=args->hist+((size_t)kc<<8);
			for(int ks=0;ks<256;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
					csizes[kc]-=freq*log2((double)freq/res);
			}
			csizes[kc]/=8;
		}
		const int poolsize=OCH_COUNT/3;
		int cidx[3]={0};
		double bestsize=0;
		for(int k0=0;k0<(int)_countof(lumapool);++k0)//215 RCTs
		{
			int kl=lumapool[k0];
			for(int k1=0;k1<(int)_countof(chroma1pool[0]);++k1)
			{
				int kc1=chroma1pool[k0][k1];
				//if(kc1/poolsize==kl/poolsize)
				//	continue;
				for(int kc2=0;kc2<OCH_COUNT;++kc2)
				{
					if(kc2/poolsize==kl/poolsize||kc2/poolsize==kc1/poolsize)
						continue;
					double size=csizes[kl]+csizes[kc1]+csizes[kc2];
					//printf("%16.6lf %d %d %d\n", size, kl, kc1, kc2);
					if(!bestsize||bestsize>size)
						bestsize=size, cidx[0]=kl, cidx[1]=kc1, cidx[2]=kc2;
				}
			}
		}
		if(cidx[0]/poolsize==cidx[1]/poolsize||cidx[0]/poolsize==cidx[2]/poolsize)
			LOG_ERROR("RCT selection error %d %d %d", cidx[0], cidx[1], cidx[2]);
		rgbidx[0]=cidx[0]/poolsize;
		rgbidx[1]=cidx[1]/poolsize;
		rgbidx[2]=cidx[2]/poolsize;
		helper1=cidx[1]%poolsize>0;
		alpha1=cidx[2]%poolsize;
		if(alpha1)
		{
			--alpha1;
			alpha2=16-alpha1;
		}
		else
		{
			alpha1=0;
			alpha2=0;
		}
		switch(rgbidx[2]<<8|rgbidx[1]<<4|rgbidx[0])
		{
		case 0x210:permutation=0;break;
		case 0x021:permutation=1;break;
		case 0x102:permutation=2;break;
		case 0x012:permutation=3;break;
		case 0x120:permutation=4;break;
		case 0x201:permutation=5;break;
		default:
			LOG_ERROR("Invalid permutation %d %d %d", rgbidx[0], rgbidx[1], rgbidx[2]);
			break;
		}
		memset(args->hist, 0, sizeof(args->hist));
		memset(args->pixels, 0, args->bufsize);
		for(int ky=args->y1+2;ky<args->y2;++ky)//analysis loop #2
		{
			int kx=args->x1+2;
			const unsigned char *ptr=image+3*(args->iw*ky+kx);
			ALIGN(16) short *rows[]=
			{
				args->pixels+((BLOCKSIZE+16LL)*((ky-0LL)&3)+8LL)*4*2,
				args->pixels+((BLOCKSIZE+16LL)*((ky-1LL)&3)+8LL)*4*2,
				args->pixels+((BLOCKSIZE+16LL)*((ky-2LL)&3)+8LL)*4*2,
				args->pixels+((BLOCKSIZE+16LL)*((ky-3LL)&3)+8LL)*4*2,
			};
			int yuv[4]={0};
			for(;kx<args->x2-2;++kx, ptr+=3)
			{
				int offset=0;
				int idx=nch*(args->iw*ky+kx);
				short
					*NNWW	=rows[2]-2*4*2,
					*NNW	=rows[2]-1*4*2,
					*NN	=rows[2]+0*4*2,
					*NNE	=rows[2]+1*4*2,
					*NNEE	=rows[2]+2*4*2,
					*NWW	=rows[1]-2*4*2,
					*NW	=rows[1]-1*4*2,
					*N	=rows[1]+0*4*2,
					*NE	=rows[1]+1*4*2,
					*NEE	=rows[1]+2*4*2,
					*WW	=rows[0]-2*4*2,
					*W	=rows[0]-1*4*2,
					*curr	=rows[0]+0*4*2;
				if(args->fwd)
				{
					yuv[0]=args->src[idx+rgbidx[0]]-128;
					yuv[1]=args->src[idx+rgbidx[1]]-128;
					yuv[2]=args->src[idx+rgbidx[2]]-128;
				}
				for(int kc=0;kc<3;++kc)
				{
					int preds[PRED_COUNT];
					int kc2=kc<<1;
					switch(kc)
					{
					case 0:
						offset=0;
						break;
					case 1:
						offset=yuv[0]&-helper1;
						break;
					case 2:
						offset=(alpha1*yuv[0]+alpha2*yuv[1])>>4;
						break;
					}
					//preds[PRED_W]=W[kc2];
					//preds[PRED_N]=N[kc2];
					MEDIAN3_32(preds[PRED_CG], N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					CLAMP3_32(preds[PRED_AV5], W[kc2]+((5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])>>3), N[kc2], W[kc2], NE[kc2]);
					CLAMP3_32(preds[PRED_AV9],
						W[kc2]+((10*N[kc2]-9*NW[kc2]+4*NE[kc2]-2*(NN[kc2]+WW[kc2])+NNW[kc2]-(NNE[kc2]+NWW[kc2]))>>4),
						N[kc2],
						W[kc2],
						NE[kc2]
					);

					preds[PRED_AV12]=(
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
					CLAMP3_32(preds[PRED_AV12], preds[PRED_AV12], N[kc2], W[kc2], NE[kc2]);

					//int
					//	gy=abs(N[kc2]-NN[kc2])+abs(W[kc2]-NW[kc2])+1,
					//	gx=abs(W[kc2]-WW[kc2])+abs(N[kc2]-NW[kc2])+1;
					//preds[PRED_WG]=(N[kc2]*gy+W[kc2]*gx)/(gy+gx);

					for(int kp=0;kp<PRED_COUNT;++kp)
					{
						int pred=preds[kp]+offset;
						CLAMP2(pred, -128, 127);
						int error=yuv[kc]-pred;
						error+=128;
						error&=255;
						++args->hist[(kc*PRED_COUNT+kp)<<8|error];
					}
					curr[kc2]=yuv[kc]-offset;
				}
				rows[0]+=4*2;
				rows[1]+=4*2;
				rows[2]+=4*2;
				rows[3]+=4*2;
			}
		}
		memset(csizes, 0, sizeof(csizes));
		res=(args->x2-args->x1-4)/5*5*(args->y2-args->y1-2);
		for(int kc=0;kc<PRED_COUNT*3;++kc)
		{
			int *curr_hist=args->hist+((size_t)kc<<8);
			for(int ks=0;ks<256;++ks)
			{
				int freq=curr_hist[ks];
				if(freq)
					csizes[kc]-=freq*log2((double)freq/res);
			}
			csizes[kc]/=8;
		}
		double bestsizes[3]={0};
		for(int kc=0;kc<3;++kc)
		{
			bestsizes[kc]=0;
			for(int kp=0;kp<PRED_COUNT;++kp)
			{
				double size=csizes[kc*PRED_COUNT+kp];
				if(!kp||bestsizes[kc]>size)
					bestsizes[kc]=size, predsel[kc]=kp;
			}
		}

#ifdef ENABLE_OLS
		//naive OLS solver		params = inv(NB * NBT) * NB * Targets	->	pred = NBT * params
		memset(args->ols_mat, 0, sizeof(args->ols_mat));
		for(int kc=0;kc<3;++kc)
		{
			double temprow[(OLS_NPARAMS0+2)*2];
			int nparams=OLS_NPARAMS0+kc;
			int matsize=nparams*nparams;
			long long *curr_imat=args->ols_imat[kc];
			//long long imax=0;
			//for(int k=0;k<matsize;++k)//X
			//{
			//	long long mag=llabs(curr_imat[k]);
			//	if(imax<mag)
			//		imax=mag;
			//}
			//double gain=1/(double)imax;
			double gain=1./(double)((long long)res<<16);
			double *curr_mat=args->ols_mat[kc];
			for(int y=0;y<nparams;++y)
			{
				for(int x=0;x<nparams;++x)
				{
					double *ptr=curr_mat+nparams*2*y+x;
					*ptr=(double)curr_imat[nparams*y+x]*gain;
					if(x==y)
						*ptr+=0.00005;//regularization
				}
			}
			for(int k=0;k<nparams;++k)
				curr_mat[nparams*2*k+k+nparams]=1;
			//double LOL_1[]=//
			//{
			//	1, 2, 2, 1, 0, 0,
			//	3, 4, 3, 0, 1, 0,
			//	9, 7, 1, 0, 0, 1
			//};
			//memcpy(curr_mat, LOL_1, sizeof(LOL_1));
			//nparams=3;
#if 0
			printf("Forward:\n");
			for(int y=0;y<nparams;++y)
			{
				for(int x=0;x<nparams;++x)
					printf("%14lf ", curr_mat[nparams*2*y+x]);
				printf("\n");
			}
			printf("\n");
#endif
			ols_success[kc]=invert_matrix(curr_mat, nparams, temprow);
			if(!ols_success[kc])
				continue;
#if 0
			printf("Inverse:  %s\n", ols_success[kc]?"SUCCESS":"FAILED");
			for(int y=0;y<nparams;++y)
			{
				for(int x=0;x<nparams;++x)
					printf("%14lf ", curr_mat[nparams*2*y+x+nparams]);
				printf("\n");
			}
			printf("\n");
#endif
			long long *curr_ivec=args->ols_ivec[kc];
			double psum=0;
			for(int kp=0;kp<nparams;++kp)
			{
				double param=0;
				for(int j=0;j<nparams;++j)
					param+=curr_mat[nparams*2*kp+j+nparams]*(double)curr_ivec[j];
				//param*=256;
				//CLAMP2(param, -0x7FFF, 0x7FFF);
				temprow[kp]=param;
				psum+=param;
			}
			ols_success[kc]=fabs(psum)>1e-9;
			if(!ols_success[kc])
				continue;
			for(int kp=0;kp<nparams;++kp)
				args->ols_params[kc][kp]=(int)(temprow[kp]*256/psum);
		}
#endif

		blist_init(&args->list);
#ifdef USE_GRCODER
		gr_enc_init(&ec, &args->list);
		gr_enc(&ec, permutation, PERM_COUNT);
		gr_enc(&ec, helper1, 2);
		gr_enc(&ec, alpha1, 17);
		gr_enc(&ec, alpha2, 17);
		gr_enc(&ec, predsel[0], PRED_COUNT);
		gr_enc(&ec, predsel[1], PRED_COUNT);
		gr_enc(&ec, predsel[2], PRED_COUNT);
#else
		ac3_enc_init(&ec, &args->list);
		ac3_enc_bypass_NPOT(&ec, permutation, PERM_COUNT);
		ac3_enc_bypass_NPOT(&ec, helper1, 2);
		ac3_enc_bypass_NPOT(&ec, alpha1, 17);
		ac3_enc_bypass_NPOT(&ec, alpha2, 17);

		ac3_enc_bypass_NPOT(&ec, predsel[0], PRED_COUNT);
		ac3_enc_bypass_NPOT(&ec, predsel[1], PRED_COUNT);
		ac3_enc_bypass_NPOT(&ec, predsel[2], PRED_COUNT);
#endif
#ifdef ENABLE_OLS
		for(int kc=0;kc<3;++kc)
		{
			int nparams=OLS_NPARAMS0+kc;
			ac3_enc_bypass(&ec, ols_success[kc], 1);
			if(ols_success[kc])
			{
				for(int kp=0;kp<nparams;++kp)
				{
					int param=args->ols_params[kc][kp]+0x8000;
					ac3_enc_bypass(&ec, param, 16);
				}
			}
		}
#endif

		args->bestsize=bestsizes[0]+bestsizes[1]+bestsizes[2];
		//args->bestsize=bestsize;
		args->permutation=permutation;
		args->helper1=helper1;
		args->alpha1=alpha1;
		args->alpha2=alpha2;
		memcpy(args->predsel, predsel, sizeof(args->predsel));
	}
	else
	{
#ifdef USE_GRCODER
		gr_dec_init(&ec, args->decstart, args->decend);
		permutation=gr_dec(&ec, PERM_COUNT);
		helper1=gr_dec(&ec, 2);
		alpha1=gr_dec(&ec, 17);
		alpha2=gr_dec(&ec, 17);
		memcpy(rgbidx, permutations[permutation], sizeof(int[3]));
		
		predsel[0]=gr_dec(&ec, PRED_COUNT);
		predsel[1]=gr_dec(&ec, PRED_COUNT);
		predsel[2]=gr_dec(&ec, PRED_COUNT);
#else
		ac3_dec_init(&ec, args->decstart, args->decend);
		permutation=ac3_dec_bypass_NPOT(&ec, PERM_COUNT);
		helper1=ac3_dec_bypass_NPOT(&ec, 2);
		alpha1=ac3_dec_bypass_NPOT(&ec, 17);
		alpha2=ac3_dec_bypass_NPOT(&ec, 17);
		memcpy(rgbidx, permutations[permutation], sizeof(int[3]));
		
		predsel[0]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predsel[1]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
		predsel[2]=ac3_dec_bypass_NPOT(&ec, PRED_COUNT);
#endif
#ifdef ENABLE_OLS
		for(int kc=0;kc<3;++kc)
		{
			int nparams=OLS_NPARAMS0+kc;
			ols_success[kc]=ac3_dec_bypass(&ec, 2);		ac3_dec_update_NPOT(&ec, ols_success[kc], ols_success[kc]+1, 2);
			if(ols_success[kc])
			{
				for(int kp=0;kp<nparams;++kp)
				{
					args->ols_params[kc][kp]=ac3_dec_bypass(&ec, 16);	ac3_dec_update_NPOT(&ec, args->ols_params[kc][kp], args->ols_params[kc][kp]+1, 0x10000);
					args->ols_params[kc][kp]-=0x8000;
				}
			}
		}
#endif
	}
#if 0
	if(args->fwd)
	{
		double csizes[OCH_COUNT*PRED_COUNT]={0}, bestsize=0;
		unsigned char predsel[OCH_COUNT]={0};
		int res=(args->x2-args->x1-3)/5*5*(args->y2-args->y1-2);
		__m256i av12_mcoeffs[12];
		
		memset(args->hist, 0, args->histsize);
		for(int k=0;k<(int)_countof(av12_mcoeffs);++k)
			av12_mcoeffs[k]=_mm256_set1_epi16(av12_icoeffs[k]>>1);
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
#ifndef DISABLE_PREDSEL
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
#endif

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
				
#ifndef DISABLE_PREDSEL
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
#endif
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
		ac3_enc_update_NPOT(&ec, bestrct, bestrct+1, RCT_COUNT);
		ac3_enc_update_NPOT(&ec, predidx[0], predidx[0]+1, PRED_COUNT);
		ac3_enc_update_NPOT(&ec, predidx[1], predidx[1]+1, PRED_COUNT);
		ac3_enc_update_NPOT(&ec, predidx[2], predidx[2]+1, PRED_COUNT);
	}
	else
	{
		ac3_dec_init(&ec, args->decstart, args->decend);
		bestrct=ac3_dec_getcdf_NPOT(&ec, RCT_COUNT);		ac3_dec_update_NPOT(&ec, bestrct, bestrct+1, RCT_COUNT);
		predidx[0]=ac3_dec_getcdf_NPOT(&ec, PRED_COUNT);	ac3_dec_update_NPOT(&ec, predidx[0], predidx[0]+1, PRED_COUNT);
		predidx[1]=ac3_dec_getcdf_NPOT(&ec, PRED_COUNT);	ac3_dec_update_NPOT(&ec, predidx[1], predidx[1]+1, PRED_COUNT);
		predidx[2]=ac3_dec_getcdf_NPOT(&ec, PRED_COUNT);	ac3_dec_update_NPOT(&ec, predidx[2], predidx[2]+1, PRED_COUNT);
	}
#endif
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
#ifdef PREDICT_SIGN
	int mixer[9*3*A2_NCTX]={0};
#else
	int mixer[8*3*(1+A2_NCTX+A2_NCTX2)]={0};
#ifdef MIXERLAYERS
	int mixer2[8*3*A2_NCTXGROUPS]={0};
	FILLMEM(mixer2, MIXERINIT, sizeof(mixer2), sizeof(int));
#endif
#endif
	FILLMEM(mixer, MIXERINIT, sizeof(mixer), sizeof(int));
	memset(args->stats0, 0, sizeof(args->stats0));
	memset(args->stats, 0, sizeof(args->stats));
	memset(args->stats2, 0, sizeof(args->stats2));
	//FILLMEM(args->stats, 0x0000, sizeof(args->stats), sizeof(short));
	//FILLMEM((short*)args->stats2, 0x0000, sizeof(args->stats2), sizeof(short));
	memset(args->sse, 0, sizeof(args->sse));
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

#ifdef ENABLE_WG
	//unsigned long long prngstate=0x3243F6A8885A308D;//pi
	int wg_weights[WG_NPREDS*3]={0};
	int wg_preds[WG_NPREDS]={0};
	long long wg_perrors[WG_NPREDS*3]={0};
	wg_init(wg_weights+WG_NPREDS*0);
	wg_init(wg_weights+WG_NPREDS*1);
	wg_init(wg_weights+WG_NPREDS*2);
#endif
#ifdef USE_GRCODER
	FILLMEM(args->pixels, 2, args->bufsize, sizeof(short));
#else
	memset(args->pixels, 0, args->bufsize);
#endif
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
		int pred=0, error=0;
		//const unsigned char *combination=rct_combinations[bestrct];
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int offset=0;
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
			(void)NNWW;
			(void)NNW;
			(void)NNEE;
			(void)NNEEE;
			(void)WWWW;
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
				yuv[0]=args->src[idx+rgbidx[0]]-128;
				yuv[1]=args->src[idx+rgbidx[1]]-128;
				yuv[2]=args->src[idx+rgbidx[2]]-128;
			}
			for(int kc=0;kc<nch;++kc)
			{
				int kc2=kc<<1;
				int ols=0;
				switch(kc)
				{
				case 0:
					offset=0;
#ifdef ENABLE_OLS
					ols=(
						args->ols_params[0][ 0]*NW	[0]+
						args->ols_params[0][ 1]*NW	[2]+
						args->ols_params[0][ 2]*NW	[6]+
						args->ols_params[0][ 3]*N	[0]+
						args->ols_params[0][ 4]*N	[2]+
						args->ols_params[0][ 5]*N	[6]+
						args->ols_params[0][ 6]*NE	[0]+
						args->ols_params[0][ 7]*NE	[2]+
						args->ols_params[0][ 8]*NE	[6]+
						args->ols_params[0][ 9]*W	[0]+
						args->ols_params[0][10]*W	[2]+
						args->ols_params[0][11]*W	[6]
					)>>8;
#endif
					break;
				case 1:
					offset=helper1?yuv[0]<<EBITS>>NBITS:0;
#ifdef ENABLE_OLS
					ols=(
						args->ols_params[1][ 0]*NW	[0]+
						args->ols_params[1][ 1]*NW	[2]+
						args->ols_params[1][ 2]*NW	[6]+
						args->ols_params[1][ 3]*N	[0]+
						args->ols_params[1][ 4]*N	[2]+
						args->ols_params[1][ 5]*N	[6]+
						args->ols_params[1][ 6]*NE	[0]+
						args->ols_params[1][ 7]*NE	[2]+
						args->ols_params[1][ 8]*NE	[6]+
						args->ols_params[1][ 9]*W	[0]+
						args->ols_params[1][10]*W	[2]+
						args->ols_params[1][11]*W	[6]+
						args->ols_params[1][12]*curr	[0]
					)>>8;
#endif
					break;
				case 2:
					offset=(alpha1*yuv[0]+alpha2*yuv[1])<<EBITS>>NBITS>>4;
#ifdef ENABLE_OLS
					ols=(
						args->ols_params[2][ 0]*NW	[0]+
						args->ols_params[2][ 1]*NW	[2]+
						args->ols_params[2][ 2]*NW	[6]+
						args->ols_params[2][ 3]*N	[0]+
						args->ols_params[2][ 4]*N	[2]+
						args->ols_params[2][ 5]*N	[6]+
						args->ols_params[2][ 6]*NE	[0]+
						args->ols_params[2][ 7]*NE	[2]+
						args->ols_params[2][ 8]*NE	[6]+
						args->ols_params[2][ 9]*W	[0]+
						args->ols_params[2][10]*W	[2]+
						args->ols_params[2][11]*W	[6]+
						args->ols_params[2][12]*curr	[0]+
						args->ols_params[2][13]*curr	[2]
					)>>8;
#endif
					break;
				}
				//int offset=(yuv[combination[kc+3]]+yuv[combination[kc+6]])>>combination[kc+9];
#if 1
				switch(predsel[kc])
				{
				//case PRED_N:
				//	pred=N[kc2];
				//	break;
				//case PRED_W:
				//	pred=W[kc2];
				//	break;
				case PRED_CG:
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
					break;
				case PRED_AV5:
					CLAMP3_32(pred, W[kc2]+((5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])>>3), N[kc2], W[kc2], NE[kc2]);
					break;
				case PRED_AV9:
					CLAMP3_32(pred,
						W[kc2]+((10*N[kc2]-9*NW[kc2]+4*NE[kc2]-2*(NN[kc2]+WW[kc2])+NNW[kc2]-(NNE[kc2]+NWW[kc2]))>>4),
						N[kc2], W[kc2], NE[kc2]
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
				//case PRED_WG:
				//	{
				//		int
				//			gy=abs(N[kc2]-NN[kc2])+abs(W[kc2]-NW[kc2])+1,
				//			gx=abs(W[kc2]-WW[kc2])+abs(N[kc2]-NW[kc2])+1;
				//		pred=(N[kc2]*gy+W[kc2]*gx)/(gy+gx);
				//	}
				//	break;
				}
#else
				CLAMP3_32(pred, W[kc2]+((5*(N[kc2]-NW[kc2])+NE[kc2]-WW[kc2])>>3), N[kc2], W[kc2], NE[kc2]);
#endif
#if 0
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
#endif
#ifdef ENABLE_WG
				pred=wg_predict(wg_weights+WG_NPREDS*kc, rows, 4*2, kc2, pred, ols, wg_perrors+WG_NPREDS*kc, wg_preds);
#endif
				pred+=offset;
				CLAMP2(pred, -(128<<EBITS>>NBITS), 127<<EBITS>>NBITS);
#ifdef USE_GRCODER
				if(args->fwd)
				{
					curr[kc2+0]=error=yuv[kc];
					error-=pred;
					error=error<<1^error>>31;
					gr_enc_POT(&ec, error, FLOOR_LOG2(W[kc2+1]+1));
					curr[kc2+1]=(2*W[kc2+1]+error+NEEE[kc2+1])>>2;
				}
				else
				{
					error=gr_dec_POT(&ec, FLOOR_LOG2(W[kc2+1]+1));
					curr[kc2+1]=(2*W[kc2+1]+error+NEEE[kc2+1])>>2;
					error=error>>1^-(error&1);
					error+=pred;
					yuv[kc]=curr[kc2+0]=error;
				}
#elif defined ENABLE_ABAC2
				unsigned long long *curr_stats0=args->stats0[kc];
				unsigned *curr_stats[A2_NCTX];
				int grads=abs(NE[kc2+0]-N[kc2+0])+abs(N[kc2+0]-NW[kc2+0])+abs(W[kc2+0]-NW[kc2+0]);
				int ctx[]=
				{
					pred*3,
					//pred*5+abs(W[kc2+1])/25,
					pred/5,
					//(W[kc2+0]+NE[kc2+0]+NEE[kc2+0]+NEEE[kc2+0])>>1,
					//N[kc2+0],
					grads*2,
					N[kc2+0]-W[kc2+0],

					(2*FLOOR_LOG2(abs(N[kc2+1]<<NBITS>>EBITS)*7)+FLOOR_LOG2(abs(W[kc2+1]<<NBITS>>EBITS)*3))<<EBITS>>NBITS,
					//N[kc2+0]+W[kc2+0]-NW[kc2+0],
					//N[kc2+0]-NN[kc2+0],
					//W[kc2+0]+WW[kc2+0],
					//NW[kc2+1],
					//N[kc2+1]+W[kc2+1]-NW[kc2+1],
					//abs(W[kc2+1]),
					((NEEE[kc2+0]<0)<<7|(NEE[kc2+0]<0)<<6|(NN[kc2+0]<0)<<5|(WW[kc2+0]<0)<<4|(NE[kc2+0]<0)<<3|(NW[kc2+0]<0)<<2|(W[kc2+0]<0)<<1|(N[kc2+0]<0))*2,
					((NEEE[kc2+1]<0)<<7|(NEE[kc2+1]<0)<<6|(NN[kc2+1]<0)<<5|(WW[kc2+1]<0)<<4|(NE[kc2+1]<0)<<3|(NW[kc2+1]<0)<<2|(W[kc2+1]<0)<<1|(N[kc2+1]<0))*2,
					(THREEWAY(N[kc2+1], 0)+1+(THREEWAY(W[kc2+1], 0)+1+(THREEWAY(NW[kc2+1], 0)+1+(THREEWAY(NE[kc2+1], 0)+1+(THREEWAY(N[kc2+1], W[kc2+1])+1)*3)*3)*3)*3),

					(kc?abs(curr[1]):N[kc2+2]*2),
					(kc>=2?abs(curr[3])*5+abs(curr[1]):W[kc2+2]/4),
					//(kc?abs(curr[0]):0)+(kc>=2?abs(curr[2]):0),
					//(kc?abs(curr[kc2-1]):abs(curr[kc2-2]))+(kc>=2?abs(curr[kc2-2]):abs(curr[kc2-3])),
					(abs(N[kc2+0])+abs(W[kc2+0]))/2-abs(W[kc2/2+1]),
					(abs(N[kc2+1])+abs(W[kc2+1]))*11-abs(NW[kc2+1])/2+abs(NE[kc2+1])*7,

					//0,
					//(WWWW[kc2+0]+WWW[kc2+0]+NEEE[kc2+0]+NNN[kc2+0])>>2,
				};
				unsigned short *curr_stats2[A2_NCTX2];
				int ctx2[][2]=
				{
					{N[kc2+0], NN[kc2+0]},
					{W[kc2+0], WW[kc2+0]},
					{N[kc2+1], W[kc2+1]},
				};
				for(int k=0;k<A2_NCTX;++k)
				{
					int idx2=((A2_NCTX*kc+k)<<A2_CTXBITS|(ctx[k]>>(9+EBITS-NBITS-A2_CTXBITS)&((1<<A2_CTXBITS)-1)))<<8;
					curr_stats[k]=args->stats+idx2;
				}
				for(int k=0;k<A2_NCTX2;++k)
				{
					int a=(ctx2[k][0]+128)>>(8-A2_CTXBITS2)&((1<<A2_CTXBITS2)-1);
					int b=(ctx2[k][1]+128)>>(8-A2_CTXBITS2)&((1<<A2_CTXBITS2)-1);
					curr_stats2[k]=args->stats2[kc][k][b][a];
				}
				//if(ky==330&&kx==183)//
				//	printf("");
				if(args->fwd)
				{
					curr[kc2+0]=yuv[kc]<<EBITS>>NBITS;
					curr[kc2+1]=curr[kc2+0]-pred;
					error=yuv[kc]-(((pred<<NBITS)+(1<<EBITS>>1))>>EBITS);
#ifdef PREDICT_SIGN
					const int half=128;
					int upred=half-abs(pred), aval=abs(error);
					if(aval<=upred)
					{
#ifdef ENABLE_BIASCORR
						{
							int negmask=ibias_corr>>31&-(error!=-half);//sign is flipped if SSE correction was negative, to skew the histogram
							error^=negmask;
							error-=negmask;
						}
#endif
						error=error<<1^error>>31;//pack sign
					}
					else
						error=aval+upred;//error sign is known	sym=2^n=256 when pixel=-2^(n-1)=-128 and pred>0
#endif
				}
				else
					error=0;
#ifdef PREDICT_SIGN
				int *curr_mixer=mixer+kc*9*A2_NCTX;
				long long *curr_sse=args->sse+(kc*9LL<<A2_SSEBITS);
				for(int kb=8, tidx=0, e2=0;kb>=0;--kb)
#else
				int *curr_mixer=mixer+kc*8*(1+A2_NCTX+A2_NCTX2);
#ifdef MIXERLAYERS
				int *curr_mixer2=mixer2+kc*8*A2_NCTXGROUPS;
#endif
				long long *curr_sse=args->sse+(kc*8LL<<A2_SSEBITS);
				for(int kb=7, tidx=1, e2=0;kb>=0;--kb)
#endif
				{
					long long p1=0;
					int wsum=0, bit;
#ifdef MIXERLAYERS
					long long probs[A2_NCTXGROUPS];
					int wsumX[A2_NCTXGROUPS];
					for(int k2=0, j=0;k2<A2_NCTXGROUPS;++k2)
					{
						probs[k2]=0;
						wsumX[k2]=0;
						for(int k=0;k<A2_CTXGROUPSIZE;++k, ++j)
						{
							int m1=curr_mixer[j];
							probs[k2]+=(long long)m1*curr_stats[j][tidx];
							wsumX[k2]+=m1;
						}
						//probs[k2]+=wsumX[k2];
						//probs[k2]/=wsumX[k2];
						probs[k2]/=MIXERINIT*A2_CTXGROUPSIZE;
						p0+=curr_mixer2[k2]*probs[k2];
						wsum+=curr_mixer2[k2];
					}
					//p0+=wsum;
					//p0/=wsum;
					p0=p0*6>>25;
					//p0/=MIXERINIT*A2_NCTXGROUPS;
#else
					int probs[1+A2_NCTX+A2_NCTX2];
					int j=0;
					{
						unsigned long long cell=curr_stats0[tidx];
						int prevseq=(cell>>62);
						const int c0=1, c1=29;
						int n0=0, n1=0;
						switch(prevseq)
						{
						case 0:
							n0=(cell>>42&0x7F)*c0+(cell>>28&0x7F)*c0+(cell>>14&0x7F)*c0+(cell>>0&0x7F)*c1;
							n1=(cell>>49&0x7F)*c0+(cell>>35&0x7F)*c0+(cell>>21&0x7F)*c0+(cell>>7&0x7F)*c1;
							break;
						case 1:
							n0=(cell>>42&0x7F)*c0+(cell>>28&0x7F)*c0+(cell>>14&0x7F)*c1+(cell>>0&0x7F)*c0;
							n1=(cell>>49&0x7F)*c0+(cell>>35&0x7F)*c0+(cell>>21&0x7F)*c1+(cell>>7&0x7F)*c0;
							break;
						case 2:
							n0=(cell>>42&0x7F)*c0+(cell>>28&0x7F)*5+(cell>>14&0x7F)*c0+(cell>>0&0x7F)*c0;
							n1=(cell>>49&0x7F)*c0+(cell>>35&0x7F)*5+(cell>>21&0x7F)*c0+(cell>>7&0x7F)*c0;
							break;
						case 3:
							n0=(cell>>42&0x7F)*c1+(cell>>28&0x7F)*c0+(cell>>14&0x7F)*c0+(cell>>0&0x7F)*c0;
							n1=(cell>>49&0x7F)*c1+(cell>>35&0x7F)*c0+(cell>>21&0x7F)*c0+(cell>>7&0x7F)*c0;
							break;
						}
						probs[j]=(int)(((n1*2LL+1)<<16)/((n0+n1)*2+2));
						p1+=(long long)curr_mixer[j]*probs[j];
						wsum+=curr_mixer[j];
						++j;
					}
					for(int k=0;k<A2_NCTX;++k)
					{
						unsigned cell=curr_stats[k][tidx];
						int prevbit=cell>>31;
						int n0, n1;
						if(prevbit)
							n0=(cell>>14&0x7F)*7+(cell>>0&0x7F), n1=(cell>>21&0x7F)*7+(cell>>7&0x7F);
						else
							n0=(cell>>14&0x7F)+(cell>>0&0x7F)*7, n1=(cell>>21&0x7F)+(cell>>7&0x7F)*7;
						//int n0=curr_stats[k][tidx]&0xFF, n1=curr_stats[k][tidx]>>8&0xFF;
						probs[j]=(int)(((n1*2LL+1)<<16)/((n0+n1)*2+2));
						p1+=(long long)curr_mixer[j]*probs[j];
						wsum+=curr_mixer[j];
						++j;
					}
					for(int k=0;k<A2_NCTX2;++k)
					{
						int n0=curr_stats2[k][tidx]&0xFF, n1=curr_stats2[k][tidx]>>8&0xFF;
						probs[j]=(int)(((n1*2LL+1)<<16)/((n0+n1)*2+2));
						//probs[j]=curr_stats2[k][tidx]*64;
						p1+=(long long)probs[j]*curr_mixer[j];
						wsum+=curr_mixer[j];
						++j;
					}
					p1/=wsum;
#endif
					int sseidx=(int)(p1>>(16-A2_SSEBITS));
					CLAMP2_32(sseidx, sseidx, 0, (1<<A2_SSEBITS)-1);
					long long ssesum=curr_sse[sseidx]>>A2_SSECTR, ssecount=curr_sse[sseidx]&((1LL<<A2_SSECTR)-1);
					p1+=ssesum/(ssecount+32);
					CLAMP2_32(p1, (int)p1, 1, 0xFFFF);
					if(args->fwd)
					{
						bit=error>>kb&1;
						ac3_enc_bin(&ec, !bit, (int)p1, 16);
					}
					else
					{
						bit=!ac3_dec_bin(&ec, (int)p1, 16);
						error|=bit<<kb;
					}
					e2|=bit<<kb;
					int prob_error=(bit<<16)-(int)p1;
					++ssecount;
					ssesum+=prob_error;
					if(ssecount>0x4A00)
					{
						ssecount>>=1;
						ssesum>>=1;
					}
					curr_sse[sseidx]=ssesum<<A2_SSECTR|ssecount;
					{
						//int grads=ctx[2]/2;
						//grads=abs(ctx[1]/2);
						//grads=pred/36;
						grads=0;
						int pred2=((char)(e2|1<<kb>>1)<<EBITS>>NBITS)-pred;
						int sh=0;
						switch(kb)
						{
						case 7:
							sh=(grads+abs(N[kc2+1])+abs(W[kc2+1])*2+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1])/2-abs(WW[kc2+1])-abs(NN[kc2+1]))>>7;
							if(sh<-1)sh=-1;
							sh=FLOOR_LOG2(sh+1)+1;
							break;
						case 6:
							sh=(grads-abs(pred2)+abs(N[kc2+1])+abs(W[kc2+1])*2+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1])/2-abs(WW[kc2+1])-abs(NN[kc2+1]))>>7;
							if(sh<-1)sh=-1;
							sh=FLOOR_LOG2(sh+1)+1;
							break;
						case 5:
							sh=(grads-abs(pred2)+abs(N[kc2+1])+abs(W[kc2+1])*2+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1])/2-abs(WW[kc2+1])-abs(NN[kc2+1]))>>7;
							if(sh<-1)sh=-1;
							sh=FLOOR_LOG2(sh+1)+1;
						//	sh=(9*(abs(N[kc2+1])+abs(W[kc2+1])+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1]))+abs(pred2))>>8;
						//	sh=FLOOR_LOG2(sh);
							break;
						case 4:
							sh=(grads-abs(pred2)/2+abs(N[kc2+1])+abs(W[kc2+1])*2+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1])/2-abs(WW[kc2+1])-abs(NN[kc2+1]))>>7;
							if(sh<-1)sh=-1;
							sh=FLOOR_LOG2(sh+1)+1;
						//	sh=(9*(abs(N[kc2+1])+abs(W[kc2+1])+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1]))+abs(pred2))>>8;
						//	sh=FLOOR_LOG2(sh);
							break;
						case 3:
							sh=(grads-abs(pred2)/4+abs(N[kc2+1])+abs(W[kc2+1])*2+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1])/2-abs(WW[kc2+1])-abs(NN[kc2+1]))>>7;
							if(sh<-1)sh=-1;
							sh=FLOOR_LOG2(sh+1)+1;
						//	sh=(9*(abs(N[kc2+1])+abs(W[kc2+1])+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1]))+abs(pred2))>>8;
						//	sh=FLOOR_LOG2(sh);
							break;
						case 2:
							sh=(grads-abs(pred2)/8+abs(N[kc2+1])+abs(W[kc2+1])*2+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1])/2-abs(WW[kc2+1])-abs(NN[kc2+1]))>>7;
							if(sh<-1)sh=-1;
							sh=FLOOR_LOG2(sh+1)+1;
						//	sh=(9*(abs(N[kc2+1])+abs(W[kc2+1])+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1]))+abs(pred2))>>8;
						//	sh=FLOOR_LOG2(sh);
							break;
						case 1:
							sh=(grads+abs(pred2)/4+abs(N[kc2+1])+abs(W[kc2+1])*2+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1])/2-abs(WW[kc2+1])-abs(NN[kc2+1]))>>7;
							if(sh<-1)sh=-1;
							sh=FLOOR_LOG2(sh+1)+1;
						//	sh=(8*(abs(N[kc2+1])+abs(W[kc2+1])+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1]))+2*abs(pred2))>>8;
						//	sh=FLOOR_LOG2(sh);
							break;
						case 0:
							sh=(abs(N[kc2+1])+abs(W[kc2+1])+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1])-abs(WW[kc2+1])+7*abs(pred2))>>8;
							if(sh<0)sh=0;
							sh=FLOOR_LOG2(sh)+1;
						//	sh=(ctx[2]/2+abs(pred2)/4+abs(N[kc2+1])+abs(W[kc2+1])*2+abs(NW[kc2+1])+abs(NE[kc2+1])+abs(NEE[kc2+1])/2-abs(WW[kc2+1])-abs(NN[kc2+1]))>>7;
						//	if(sh<-1)sh=-1;
						//	sh=FLOOR_LOG2(sh+1);
							break;
						}
						sh+=NBITS-EBITS;
						if(sh<0)sh=0;
						if(abs(prob_error)>64)
						{
#ifdef MIXERLAYERS
#if 0
							//p[i] = (mA*pa+mB*pb+mC*pc+mD*pd)/(mA+mB+mC+mD)
							//p = (m1*p1+m2*p2+m3*p3)/(m1+m2+m3)
							//L = -log2(pbit)
							//dL/dm[i] = 1/(bit-p) * (p[i]-p)/wsum				'A2_NCTXGROUPS' values from mixer2
							//dL/dm[X] = 1/(bit-p) * m[i]/wsum * (p[X]-p[i])/wsumX[i]	'A2_NCTX' values from mixer
							long long dL_dp0=0x7C000000000/((bit<<16)-(int)p0);
							for(int k2=0, j=0;k2<A2_NCTXGROUPS;++k2)
							{
								long long dL_dpi=((long long)dL_dp0*curr_mixer2[k2]/wsum);
								for(int k=0;k<A2_CTXGROUPSIZE;++k, ++j)
								{
									int mk=curr_mixer2[j]-(int)((long long)dL_dpi*(curr_stats[j][tidx]-probs[k2])/wsumX[k2]>>9);
									//CLAMP2_32(mk, mk, 1, MIXERCLAMP);
									curr_mixer[j]=mk;
								}
								int mk=curr_mixer2[k2]-(int)(dL_dp0*(probs[k2]-(int)p0)/wsum>>9);
								//CLAMP2_32(mk, mk, 1, MIXERCLAMP);
								curr_mixer2[k2]=mk;
							}
#else
							int error=(!bit<<16)-(int)p0;
							for(int k2=0, j=0;k2<A2_NCTXGROUPS;++k2)
							{
								for(int k=0;k<A2_CTXGROUPSIZE;++k, ++j)
								{
									int mk=curr_mixer2[j];
									mk+=(int)((long long)error*curr_stats[j][tidx]>>24);
									//CLAMP2_32(mk, mk, -MIXERCLAMP, MIXERCLAMP);
									curr_mixer[j]=mk;
								}
								int mk=curr_mixer2[k2];
								mk+=(int)((long long)error*probs[k2]>>24);
								//CLAMP2_32(mk, mk, -MIXERCLAMP, MIXERCLAMP);
								curr_mixer2[k2]=mk;
							}
#endif
#else
							long long dL_dp1=0x1600000000/((!bit<<16)-(int)p1);
							for(int k=0;k<1+A2_NCTX+A2_NCTX2;++k)
							{
								int mk=curr_mixer[k]-(int)(dL_dp1*(probs[k]-(int)p1)/wsum>>9);
								CLAMP2_32(mk, mk, 1, 0xC000);
								curr_mixer[k]=mk;
							}
#endif
						}
						{
							unsigned long long cell=curr_stats0[tidx];
							int prevseq=(cell>>62);
							int ctrs[]=
							{
								(int)(cell>> 0&0x7F),
								(int)(cell>> 7&0x7F),
								(int)(cell>>14&0x7F),
								(int)(cell>>21&0x7F),
								(int)(cell>>28&0x7F),
								(int)(cell>>35&0x7F),
								(int)(cell>>42&0x7F),
								(int)(cell>>49&0x7F),
							};
							int *c=ctrs+prevseq*2;
							//switch(prevseq)
							//{
							//case 0:
							//	c[0]=cell>>0&0x7F;
							//	c[1]=cell>>7&0x7F;
							//	break;
							//case 1:
							//	c[0]=cell>>14&0x7F;
							//	c[1]=cell>>21&0x7F;
							//	break;
							//case 2:
							//	c[0]=cell>>28&0x7F;
							//	c[1]=cell>>35&0x7F;
							//	break;
							//case 3:
							//	c[0]=cell>>42&0x7F;
							//	c[1]=cell>>49&0x7F;
							//	break;
							//}
							if(c[bit]>=127)
							{
								c[0]>>=1;
								c[1]>>=1;
							}
							++c[bit];
							prevseq=prevseq<<1|bit;
							curr_stats0[tidx]=
								(unsigned long long)prevseq<<62|
								(unsigned long long)c[7]<<49|
								(unsigned long long)c[6]<<42|
								(unsigned long long)c[5]<<35|
								(unsigned long long)c[4]<<28|
								(unsigned long long)c[3]<<21|
								(unsigned long long)c[2]<<14|
								(unsigned long long)c[1]<<7|
								(unsigned long long)c[0];

						}
						for(int k=0;k<A2_NCTX;++k)
						{
							unsigned cell=curr_stats[k][tidx];
							int prevbit=cell>>31;
							int c[2];
							if(prevbit)
								c[0]=cell>>14&0x7F, c[1]=cell>>21&0x7F;
							else
								c[0]=cell>>0&0x7F, c[1]=cell>>7&0x7F;
							if(c[bit]>=127)
							{
								c[0]>>=1;
								c[1]>>=1;
							}
							++c[bit];
							if(prevbit)
								cell=bit<<31|c[1]<<21|c[0]<<14|cell&0x3FFF;
							else
								cell=bit<<31|cell&0xFFFC000|c[1]<<7|c[0];
							curr_stats[k][tidx]=cell;

							//int c[]={curr_stats[k][tidx]&0xFF, curr_stats[k][tidx]>>8&0xFF};
							//if(c[bit]>=255)
							//{
							//	c[0]>>=1;
							//	c[1]>>=1;
							//}
							//++c[bit];
							//curr_stats[k][tidx]=c[1]<<8|c[0];

							//p+=((!bit<<16)-p)>>(sh+4+(k>A2_NCTX/2));
							//CLAMP2_32(p, p, 1, 0xFFFF);
							//curr_stats[k][tidx]=p;
						}
						for(int k=0;k<A2_NCTX2;++k)
						{
							int c[]={curr_stats2[k][tidx]&0xFF, curr_stats2[k][tidx]>>8&0xFF};
							if(c[bit]>=127)
							{
								c[0]>>=1;
								c[1]>>=1;
							}
							++c[bit];
							curr_stats2[k][tidx]=c[1]<<8|c[0];
						}
						//j=A2_NCTX;
						//for(int k=0;k<A2_NCTX2;++k)
						//{
						//	unsigned char run=probs[j];
						//	int bit2=run>0, count=abs(run);
						//	if(bit==bit2)
						//	{
						//		xorshift64(&prngstate);
						//		count+=!(prngstate&((1ULL<<(count+1))-1));
						//		if(count>(1<<8>>1)-1)
						//			count=(1<<8>>1)-1;
						//	}
						//	else
						//		bit2=bit, count=1;
						//	run=bit2?count:-count;
						//	curr_stats2[k][tidx]=run;
						//	++j;
						//}
					}
					curr_mixer+=1+A2_NCTX+A2_NCTX2;
#ifdef MIXERLAYERS
					//curr_mixer2+=A2_NCTXGROUPS;//X
#endif
					curr_sse+=1<<A2_SSEBITS;
#ifdef PREDICT_SIGN
					tidx+=tidx+!bit;
#else
					tidx+=tidx+bit;
#endif
				}
				if(!args->fwd)
				{
#ifdef PREDICT_SIGN
					const int half=128;
					int upred=half-abs(pred), negmask=0;
					if(error<=(upred<<1))
					{
						error=error>>1^-(error&1);
#ifdef ENABLE_BIASCORR
						negmask=ibias_corr>>31&-(error!=-half);
#endif
					}
					else
					{
						error=error-upred;
						negmask=(-pred)>>31;
					}
					error^=negmask;
					error-=negmask;
					error+=pred;
#else
					error+=((pred<<NBITS)+(1<<EBITS>>1))>>EBITS;
					error=error<<(32-8)>>(32-8);
#endif
					yuv[kc]=error;
					curr[kc2+0]=yuv[kc]<<EBITS>>NBITS;
					curr[kc2+1]=curr[kc2+0]-pred;
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
				int token=0, bypass=0, nbits=0;
				const int depth=8;
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
#if 0
#define MIXCDF(X)\
	f28_mix8(\
		(curr_hist[0][X]<<14)/curr_hist[0][args->tlevels],\
		(curr_hist[1][X]<<14)/curr_hist[1][args->tlevels],\
		(curr_hist[2][X]<<14)/curr_hist[2][args->tlevels],\
		(curr_hist[3][X]<<14)/curr_hist[3][args->tlevels],\
		(curr_hist[4][X]<<14)/curr_hist[4][args->tlevels],\
		(curr_hist[5][X]<<14)/curr_hist[5][args->tlevels],\
		(curr_hist[6][X]<<14)/curr_hist[6][args->tlevels],\
		(curr_hist[7][X]<<14)/curr_hist[7][args->tlevels],\
		alphas[0],\
		alphas[1],\
		alphas[2]\
	)
#endif
#if 0
#define MIXCDF(X)\
	f28_mix8(\
		curr_hist[0][X],\
		curr_hist[1][X],\
		curr_hist[2][X],\
		curr_hist[3][X],\
		curr_hist[4][X],\
		curr_hist[5][X],\
		curr_hist[6][X],\
		curr_hist[7][X],\
		alphas[0],\
		alphas[1],\
		alphas[2]\
	)
#endif
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
				const int depth=8;
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
				int sym=0;
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
				if(args->fwd)
				{
					curr[kc2+0]=yuv[kc];
					curr[kc2+1]=error=yuv[kc]-pred;
#ifndef ENABLE_ABAC
					{
						const int half=128;
						int upred=half-abs(pred), aval=abs(error);
						if(aval<=upred)
						{
							sym=error;
#ifdef ENABLE_BIASCORR
							{
								int negmask=-((ibias_corr<0)&(sym!=-half));//sign is flipped if SSE correction was negative, to skew the histogram
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
						const int half=128;
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
#ifndef USE_GRCODER
				//NWW[kc2+1]+=curr[kc2+1]>>1;//bad
				//NW[kc2+1]+=curr[kc2+1]<<1;
				//N[kc2+1]+=curr[kc2+1]>>4;
				//NE[kc2+1]+=curr[kc2+1]<<1;
				//NEE[kc2+1]+=curr[kc2+1]>>1;
				//NEEE[kc2+1]+=curr[kc2+1]>>1;
				//W[kc2+1]-=N[kc2+1];//bad
				//curr[kc2+1]=(2*W[kc2+1]+curr[kc2+1]+NEEE[kc2+1])>>2;//bad
#ifdef ENABLE_WG
				wg_update(curr[kc2+0], wg_preds, wg_perrors+WG_NPREDS*kc, wg_weights+WG_NPREDS*kc);
#endif
#endif
			}
			if(!args->fwd)
			{
				args->dst[idx+rgbidx[0]]=yuv[0]+128;
				args->dst[idx+rgbidx[1]]=yuv[1]+128;
				args->dst[idx+rgbidx[2]]=yuv[2]+128;
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
#ifdef USE_GRCODER
		gr_enc_flush(&ec);
#else
		ac3_enc_flush(&ec);
#endif
}
int c02_codec(int argc, char **argv)
{
	if(argc!=2&&argc!=3&&argc!=4)
	{
		printf(
			"Usage: \"%s\"  input  output  [maxthreads]    Encode/decode.\n"
			"       \"%s\"  input                          Test without saving.\n"
			"[maxthreads]:\n"
			"  0: nthreads = number of cores (default)\n"
			"  1: Single thread\n"
			, argv[0]
			, argv[0]
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argc>2?argv[2]:0;
	int maxthreads=argc==4?atoi(argv[3]):0;
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
	int tlevels;
//	int histsize;
//	int statssize;
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
	nblocks=xblocks*yblocks;
	nthreads=MINVAR(nblocks, ncores);
	if(maxthreads&&nthreads>maxthreads)
		nthreads=maxthreads;
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
	else
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
		//statssize=0;
#elif defined ENABLE_ABAC
		statssize=nch*(int)sizeof(short[ABAC_NCTRS*ABAC_NCTX*ABAC_CLEVELS*ABAC_TREESIZE]);
#else
		clevels=CLEVELS;
		//statssize=clevels*clevels*(tlevels+1)*nch*(int)sizeof(int);
#endif
		//histsize=(int)sizeof(int[OCH_COUNT*PRED_COUNT<<8]);
		//if(histsize<statssize)
		//	histsize=statssize;
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=image;
		arg->dst=fwd?0:dst->data+printed;
		arg->iw=iw;
		arg->ih=ih;
		arg->bufsize=sizeof(short[4*4*2*(BLOCKSIZE+16LL)]);//4 padded rows * 4 channels max * {pixels, wg_errors}
		arg->pixels=(short*)_mm_malloc(arg->bufsize, sizeof(__m128i));

		//arg->histsize=histsize;
		//arg->hist=(int*)malloc(histsize);
#ifdef ENABLE_ABAC
		arg->statssize=statssize;
#endif
		//arg->stats=(unsigned*)malloc(statssize);
		if(!arg->pixels)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=arg->bufsize;
		//memusage+=histsize;
		//memusage+=statssize;
		
		arg->tlevels=tlevels;
#ifdef ENABLE_MIX4
		arg->clevels=clevels;
#endif
		arg->fwd=fwd;
		arg->test=test;
#ifdef ENABLE_MT
		arg->loud=0;
#else
		arg->loud=test&&nblocks<MAXPRINTEDBLOCKS;
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
#ifdef ENABLE_MT
			void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
			mt_finish(ctx);
#else
			for(int k=0;k<nthreads2;++k)
				block_thread(args+k);
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
								"block %4d/%4d  XY %3d %3d  %4d*%4d:  %8d->%16lf->%8zd bytes (%+10.2lf)  %10.6lf%%  CR %10lf  %s %d (%2d+%2d)/16  %s %s %s\n",
								kt+kt2+1, nblocks,
								kx, ky,
								arg->y2-arg->y1,
								arg->x2-arg->x1,
								blocksize,
								arg->bestsize,
								arg->list.nbytes,
								arg->list.nbytes-arg->bestsize,
								100.*arg->list.nbytes/blocksize,
								(double)blocksize/arg->list.nbytes,
								permnames[arg->permutation],
								arg->helper1,
								arg->alpha1,
								arg->alpha2,
								pred_names[arg->predsel[0]],
								pred_names[arg->predsel[1]],
								pred_names[arg->predsel[2]]
							);
						}
						esize+=arg->bestsize;
#ifdef ABAC_PROFILESIZE
						for(int k=0;k<ABAC_TOKEN_BITS*3;++k)
							abac_csizes[k]+=arg->abac_csizes[k];
#endif
					}
					memcpy(dst->data+start+sizeof(int)*((ptrdiff_t)kt+kt2), &arg->list.nbytes, sizeof(int));
					blist_appendtoarray(&arg->list, &dst);
					blist_clear(&arg->list);
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
		//free(arg->hist);
		//free(arg->stats);
	}
	free(args);
	array_free(&src);
	array_free(&dst);
	//if(test)
	//	pause();
	return 0;
}