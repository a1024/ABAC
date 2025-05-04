#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT

//	#define DISABLE_CRCT
//	#define DISABLE_RCT2
	#define ENABLE_RLE
	#define ENABLE_NEG


#include"ac.h"
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
#define BLOCKDX 512
#define BLOCKDY 512
#define MAXPRINTEDBLOCKS 500

#define PREC_BITS 8

#define RLE_BITS 7
#define RLE_ESCAPESYM 0


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
const int chroma1pool[3][4]=
{
	{OCH_Gbr00, OCH_Gbr0G, OCH_Brg00, OCH_BrgG0},//R
	{OCH_Rgb00, OCH_RgbG0, OCH_Brg00, OCH_Brg0G},//G
	{OCH_Rgb00, OCH_Rgb0G, OCH_Gbr00, OCH_GbrG0},//B
};

#define RCT2LIST\
	RCT2(JPEG2000)\
	RCT2(Pei09)\
	RCT2(RCT1)
typedef enum _RCT2Type
{
#define RCT2(X) RCT2_##X,
	RCT2LIST
#undef  RCT2
	RCT2_COUNT
} RCT2Type;
const char *rct2names[]=
{
#define RCT2(X) #X,
	RCT2LIST
#undef  RCT2
};

#define PREDLIST\
	PRED(N)\
	PRED(W)\
	PRED(CG)\
	PRED(AV5)\
	PRED(AV9)\
	PRED(AV12)\
	PRED(WG)
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
//typedef enum _NBIndex
//{
//	NB_NW,		NB_N,
//	NB_W,		NB_curr,
//} NBIndex;
//typedef enum _N2Index
//{
//	N2_NNWW,	N2_NNW,		N2_NN,		N2_NNE,		N2_NNEE,
//	N2_NWW,		N2_NW,		N2_N,		N2_NE,		N2_NEE,
//	N2_WW,		N2_W,		N2_curr,
//} N2Index;
static const short av12_icoeffs[12]=
{
	 0x04,	 0x03,	-0x1F,	-0x26,	 0x00,
	 0x07,	-0x9E,	 0xDB,	 0x1E,	 0x13,
	-0x2A,	 0xF3,
};


typedef struct _ThreadArgs
{
	const Image *src;
	Image *dst;
	
	int fwd, loud, x1, x2, y1, y2;
	int pixels[(BLOCKDX+16)*4*(OCH_COUNT+PERM_COUNT*RCT2_COUNT*3)*2];//4 padded rows  *  EXTENDED_OCH_COUNT

	DList list;
	const unsigned char *decstart, *decend;

	
	//OCH_COUNT				* {+, -} * PRED_COUNT
	//PERM_COUNT * RCT_COUNT * {YUV899}	* {+, -} * PRED_COUNT
	int hist[(OCH_COUNT + RCT2_COUNT*PERM_COUNT*(1+2+2)) * PRED_COUNT*2 << 8];
	
	int blockidx;
	double bestsize;
	int use_proper_rct;
	int rct2_idx;
	int perm_idx;
	int helper1, alpha1, alpha2;
	int predidx[3];
	int pneg[3];
#ifdef ENABLE_RLE
	int use_rle;
#endif
	size_t nbypass;
} ThreadArgs;
static double calc_csize_fromhist(int *hist, int nlevels, int res)
{
	//int sum=0;
	//for(int k=0;k<nlevels;++k)//
	//	sum+=hist[k];
	//if(sum!=res)//
	//	LOG_ERROR("");

	double csize=0;
	for(int ks=0;ks<nlevels;++ks)
	{
		int freq=hist[ks];
		if(freq)
			csize-=freq*log2((double)freq/res);
	}
	return csize/8;
}
static void update(int *hist, int **rows, int kc, int offset, int neg, int pixel, int srcdepth, int dstdepth, int ehalf)
{
	int
		NNWW	=rows[2][kc-2*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2],
		NNW	=rows[2][kc-1*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2],
		NN	=rows[2][kc+0*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2],
		NNE	=rows[2][kc+1*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2],
		NNEE	=rows[2][kc+2*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2],
		NWW	=rows[1][kc-2*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2],
		NW	=rows[1][kc-1*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2],
		N	=rows[1][kc+0*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2],
		NE	=rows[1][kc+1*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2],
		NEE	=rows[1][kc+2*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2],
		WW	=rows[0][kc-2*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2],
		W	=rows[0][kc-1*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2],
		*curr	=rows[0]+kc+0*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2;
	int preds[PRED_COUNT];

	preds[PRED_W]=W;
	preds[PRED_N]=N;
	MEDIAN3_32(preds[PRED_CG], N, W, N+W-NW);
	CLAMP3_32(preds[PRED_AV5], W+((5*(N-NW)+NE-WW)>>3), N, W, NE);
	CLAMP3_32(preds[PRED_AV9], W+((10*N-9*NW+4*NE-2*(NN+WW)+NNW-(NNE+NWW))>>4), N, W, NE);
	
	preds[PRED_AV12]=(
		av12_icoeffs[ 0]*NNWW+
		av12_icoeffs[ 1]*NNW+
		av12_icoeffs[ 2]*NN+
		av12_icoeffs[ 3]*NNE+
		av12_icoeffs[ 4]*NNEE+
		av12_icoeffs[ 5]*NWW+
		av12_icoeffs[ 6]*NW+
		av12_icoeffs[ 7]*N+
		av12_icoeffs[ 8]*NE+
		av12_icoeffs[ 9]*NEE+
		av12_icoeffs[10]*WW+
		av12_icoeffs[11]*W
	)>>8;
	CLAMP3_32(preds[PRED_AV12], preds[PRED_AV12], N, W, NE);
	
	int
		gy=abs(N-NN)+abs(W-NW)+1,
		gx=abs(W-WW)+abs(N-NW)+1;
	preds[PRED_WG]=(N*gy+W*gx)/(gy+gx);
	
	int dstlevels=1<<dstdepth;
	for(int kp=0;kp<PRED_COUNT;++kp)
	{
		int pred=preds[kp];
		if(neg)
			pred=-pred;
		pred+=offset;
		CLAMP2_32(pred, pred, -ehalf, ehalf-1);
		int error=pixel-((pred+(1<<PREC_BITS>>1))>>PREC_BITS);
		error>>=srcdepth-dstdepth;
		error+=dstlevels>>1;
		error&=dstlevels-1;
		++hist[error];
		hist+=dstlevels;
	}
	*curr=(pixel<<PREC_BITS)-offset;
	if(neg)
		*curr=-*curr;
}
static void block_thread(void *param)
{
	GolombRiceCoder ec;
	ThreadArgs *args=(ThreadArgs*)param;
	const Image *image=args->fwd?args->src:args->dst;
	int nlevels=1<<image->depth, half=nlevels>>1, ehalf=half<<PREC_BITS, mask=nlevels-1;
#ifdef ENABLE_RLE
	int use_rle=0;
#endif
	int use_proper_rct=0, mod_sh[3]={0};
	int perm_idx=0;
	int rct2_idx=0;
	int helper1=0, alpha1=0, alpha2;
	int rgbidx[3]={0};//derived from permutation
	int pneg[3]={0}, predidx[3]={0};

	if(args->fwd)
	{
		int res=(args->x2-args->x1)*(args->y2-args->y1);
		memset(args->hist, 0, sizeof(args->hist));
		memset(args->pixels, 0, sizeof(args->pixels));
		//unsigned char flags[(OCH_COUNT + RCT2_COUNT*PERM_COUNT*3)*2]={0};
		for(int ky=args->y1;ky<args->y2;++ky)//analysis loop
		{
			int kx=args->x1;
			const short *ptr=image->data+3*(image->iw*ky+kx);
			ALIGN(16) int *rows[]=
			{
				args->pixels+((BLOCKDX+16LL)*((ky-0LL)&3)+8LL)*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2,
				args->pixels+((BLOCKDX+16LL)*((ky-1LL)&3)+8LL)*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2,
				args->pixels+((BLOCKDX+16LL)*((ky-2LL)&3)+8LL)*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2,
				args->pixels+((BLOCKDX+16LL)*((ky-3LL)&3)+8LL)*(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2,
			};
			int rgb[3]={0};
			int yuv0[3], yuv[3];
			for(;kx<args->x2;++kx, ptr+=3)
			{
				int idx=image->nch*(image->iw*ky+kx), hidx, offset=0;
				rgb[0]=image->data[idx+0];
				rgb[1]=image->data[idx+1];
				rgb[2]=image->data[idx+2];

#ifndef DISABLE_CRCT
				hidx=OCH_Rgb00<<1|0; update(args->hist+((size_t)hidx*PRED_COUNT<<8), rows, hidx, 0, 0, rgb[0], image->depth, 8, ehalf);	//++flags[hidx];
				hidx=OCH_Rgb00<<1|1; update(args->hist+((size_t)hidx*PRED_COUNT<<8), rows, hidx, 0, 1, rgb[0], image->depth, 8, ehalf);	//++flags[hidx];
				hidx=OCH_Gbr00<<1|0; update(args->hist+((size_t)hidx*PRED_COUNT<<8), rows, hidx, 0, 0, rgb[1], image->depth, 8, ehalf);	//++flags[hidx];
				hidx=OCH_Gbr00<<1|1; update(args->hist+((size_t)hidx*PRED_COUNT<<8), rows, hidx, 0, 1, rgb[1], image->depth, 8, ehalf);	//++flags[hidx];
				hidx=OCH_Brg00<<1|0; update(args->hist+((size_t)hidx*PRED_COUNT<<8), rows, hidx, 0, 0, rgb[2], image->depth, 8, ehalf);	//++flags[hidx];
				hidx=OCH_Brg00<<1|1; update(args->hist+((size_t)hidx*PRED_COUNT<<8), rows, hidx, 0, 1, rgb[2], image->depth, 8, ehalf);	//++flags[hidx];
				for(int k=0;k<=16;++k)
				{
					int offset[]=
					{
#define MIX4(A, B) ((B)+(((A)-(B))*k>>4))
						MIX4(rgb[1]<<PREC_BITS, rgb[2]<<PREC_BITS),
						MIX4(rgb[2]<<PREC_BITS, rgb[0]<<PREC_BITS),
						MIX4(rgb[0]<<PREC_BITS, rgb[1]<<PREC_BITS),
#undef  MIX4
					};
					hidx=(OCH_Rgb0G+k)<<1|0; update(args->hist+((size_t)hidx*PRED_COUNT<<8), rows, hidx, offset[0], 0, rgb[0], image->depth, 8, ehalf);	//++flags[hidx];
					hidx=(OCH_Rgb0G+k)<<1|1; update(args->hist+((size_t)hidx*PRED_COUNT<<8), rows, hidx, offset[0], 1, rgb[0], image->depth, 8, ehalf);	//++flags[hidx];
					hidx=(OCH_Gbr0G+k)<<1|0; update(args->hist+((size_t)hidx*PRED_COUNT<<8), rows, hidx, offset[1], 0, rgb[1], image->depth, 8, ehalf);	//++flags[hidx];
					hidx=(OCH_Gbr0G+k)<<1|1; update(args->hist+((size_t)hidx*PRED_COUNT<<8), rows, hidx, offset[1], 1, rgb[1], image->depth, 8, ehalf);	//++flags[hidx];
					hidx=(OCH_Brg0G+k)<<1|0; update(args->hist+((size_t)hidx*PRED_COUNT<<8), rows, hidx, offset[2], 0, rgb[2], image->depth, 8, ehalf);	//++flags[hidx];
					hidx=(OCH_Brg0G+k)<<1|1; update(args->hist+((size_t)hidx*PRED_COUNT<<8), rows, hidx, offset[2], 1, rgb[2], image->depth, 8, ehalf);	//++flags[hidx];
				}
#endif
#ifndef DISABLE_RCT2
				int *phist=args->hist+OCH_COUNT*PRED_COUNT*2*256;
				for(int kp=0;kp<PERM_COUNT;++kp)
				{
					int *perm=permutations[kp];
					yuv0[0]=rgb[perm[0]];
					yuv0[1]=rgb[perm[1]];
					yuv0[2]=rgb[perm[2]];

					memcpy(yuv, yuv0, sizeof(yuv));
					yuv[1]-=yuv[0];
					yuv[2]-=yuv[0];
					yuv[0]+=(yuv[1]+yuv[2])>>2;
					hidx=((RCT2_JPEG2000*PERM_COUNT+kp)*3+0)<<1|0; update(phist, rows, OCH_COUNT*2+hidx, 0, 0, yuv[0], image->depth+0, 8, ehalf); phist+=PRED_COUNT*1<<8;	//++flags[hidx];
					hidx=((RCT2_JPEG2000*PERM_COUNT+kp)*3+0)<<1|1; update(phist, rows, OCH_COUNT*2+hidx, 0, 1, yuv[0], image->depth+0, 8, ehalf); phist+=PRED_COUNT*1<<8;	//++flags[hidx];
					hidx=((RCT2_JPEG2000*PERM_COUNT+kp)*3+1)<<1|0; update(phist, rows, OCH_COUNT*2+hidx, 0, 0, yuv[1], image->depth+1, 9, ehalf); phist+=PRED_COUNT*2<<8;	//++flags[hidx];
					hidx=((RCT2_JPEG2000*PERM_COUNT+kp)*3+1)<<1|1; update(phist, rows, OCH_COUNT*2+hidx, 0, 1, yuv[1], image->depth+1, 9, ehalf); phist+=PRED_COUNT*2<<8;	//++flags[hidx];
					hidx=((RCT2_JPEG2000*PERM_COUNT+kp)*3+2)<<1|0; update(phist, rows, OCH_COUNT*2+hidx, 0, 0, yuv[2], image->depth+1, 9, ehalf); phist+=PRED_COUNT*2<<8;	//++flags[hidx];
					hidx=((RCT2_JPEG2000*PERM_COUNT+kp)*3+2)<<1|1; update(phist, rows, OCH_COUNT*2+hidx, 0, 1, yuv[2], image->depth+1, 9, ehalf); phist+=PRED_COUNT*2<<8;	//++flags[hidx];
					
					memcpy(yuv, yuv0, sizeof(yuv));
					yuv[1]-=(87*yuv[2]+169*yuv[0]+128)>>8;
					yuv[2]-=yuv[0];
					yuv[0]+=(86*yuv[2]+29*yuv[1]+128)>>8;
					hidx=((RCT2_Pei09*PERM_COUNT+kp)*3+0)<<1|0; update(phist, rows, OCH_COUNT*2+hidx, 0, 0, yuv[0], image->depth+0, 8, ehalf); phist+=PRED_COUNT*1<<8;	//++flags[hidx];
					hidx=((RCT2_Pei09*PERM_COUNT+kp)*3+0)<<1|1; update(phist, rows, OCH_COUNT*2+hidx, 0, 1, yuv[0], image->depth+0, 8, ehalf); phist+=PRED_COUNT*1<<8;	//++flags[hidx];
					hidx=((RCT2_Pei09*PERM_COUNT+kp)*3+1)<<1|0; update(phist, rows, OCH_COUNT*2+hidx, 0, 0, yuv[1], image->depth+1, 9, ehalf); phist+=PRED_COUNT*2<<8;	//++flags[hidx];
					hidx=((RCT2_Pei09*PERM_COUNT+kp)*3+1)<<1|1; update(phist, rows, OCH_COUNT*2+hidx, 0, 1, yuv[1], image->depth+1, 9, ehalf); phist+=PRED_COUNT*2<<8;	//++flags[hidx];
					hidx=((RCT2_Pei09*PERM_COUNT+kp)*3+2)<<1|0; update(phist, rows, OCH_COUNT*2+hidx, 0, 0, yuv[2], image->depth+1, 9, ehalf); phist+=PRED_COUNT*2<<8;	//++flags[hidx];
					hidx=((RCT2_Pei09*PERM_COUNT+kp)*3+2)<<1|1; update(phist, rows, OCH_COUNT*2+hidx, 0, 1, yuv[2], image->depth+1, 9, ehalf); phist+=PRED_COUNT*2<<8;	//++flags[hidx];
					
					memcpy(yuv, yuv0, sizeof(yuv));
					yuv[2]-=yuv[0];
					yuv[0]+=yuv[2]>>1;
					yuv[1]-=yuv[0];
					yuv[0]+=yuv[1]>>1;
					hidx=((RCT2_RCT1*PERM_COUNT+kp)*3+0)<<1|0; update(phist, rows, OCH_COUNT*2+hidx, 0, 0, yuv[0], image->depth+0, 8, ehalf); phist+=PRED_COUNT*1<<8;	//++flags[hidx];
					hidx=((RCT2_RCT1*PERM_COUNT+kp)*3+0)<<1|1; update(phist, rows, OCH_COUNT*2+hidx, 0, 1, yuv[0], image->depth+0, 8, ehalf); phist+=PRED_COUNT*1<<8;	//++flags[hidx];
					hidx=((RCT2_RCT1*PERM_COUNT+kp)*3+1)<<1|0; update(phist, rows, OCH_COUNT*2+hidx, 0, 0, yuv[1], image->depth+1, 9, ehalf); phist+=PRED_COUNT*2<<8;	//++flags[hidx];
					hidx=((RCT2_RCT1*PERM_COUNT+kp)*3+1)<<1|1; update(phist, rows, OCH_COUNT*2+hidx, 0, 1, yuv[1], image->depth+1, 9, ehalf); phist+=PRED_COUNT*2<<8;	//++flags[hidx];
					hidx=((RCT2_RCT1*PERM_COUNT+kp)*3+2)<<1|0; update(phist, rows, OCH_COUNT*2+hidx, 0, 0, yuv[2], image->depth+1, 9, ehalf); phist+=PRED_COUNT*2<<8;	//++flags[hidx];
					hidx=((RCT2_RCT1*PERM_COUNT+kp)*3+2)<<1|1; update(phist, rows, OCH_COUNT*2+hidx, 0, 1, yuv[2], image->depth+1, 9, ehalf); phist+=PRED_COUNT*2<<8;	//++flags[hidx];
				}
				//if(phist!=args->hist+_countof(args->hist))
				//	LOG_ERROR("");
#endif
				rows[0]+=(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2;
				rows[1]+=(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2;
				rows[2]+=(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2;
				rows[3]+=(OCH_COUNT+RCT2_COUNT*3*PERM_COUNT)*2;
			}
		}
		//for(int k=1;k<(int)_countof(flags);++k)
		//{
		//	if(flags[k]!=*flags)
		//		LOG_ERROR("");
		//}
		double csizes[(OCH_COUNT + RCT2_COUNT*PERM_COUNT*3) * PRED_COUNT*2]={0};
		{
			int *phist=args->hist+OCH_COUNT*PRED_COUNT*2*256;
#ifndef DISABLE_CRCT
			for(int kc=0;kc<OCH_COUNT*PRED_COUNT*2;++kc)
				csizes[kc]=calc_csize_fromhist(args->hist+((size_t)kc<<8), 1<<8, res);
#endif
#ifndef DISABLE_RCT2
			for(int kc=0;kc<PERM_COUNT*RCT2_COUNT*3;kc+=3)
			{
				for(int kc2=0;kc2<PRED_COUNT*2;++kc2)
					csizes[(kc+0)*PRED_COUNT*2+kc2+OCH_COUNT*PRED_COUNT*2]=calc_csize_fromhist(phist, 1<<8, res), phist+=1<<8;
				for(int kc2=0;kc2<PRED_COUNT*2;++kc2)
					csizes[(kc+1)*PRED_COUNT*2+kc2+OCH_COUNT*PRED_COUNT*2]=calc_csize_fromhist(phist, 2<<8, res), phist+=2<<8;
				for(int kc2=0;kc2<PRED_COUNT*2;++kc2)
					csizes[(kc+2)*PRED_COUNT*2+kc2+OCH_COUNT*PRED_COUNT*2]=calc_csize_fromhist(phist, 2<<8, res), phist+=2<<8;
			}
			//if(phist!=args->hist+_countof(args->hist))
			//	LOG_ERROR("");
#endif
		}

		const int poolsize=OCH_COUNT/3;
		int cidx[3]={0};
		double bestsize=0;
		int predsel[OCH_COUNT+RCT2_COUNT*PERM_COUNT*3]={0};
		for(int kc=0;kc<(int)_countof(predsel);++kc)
		{
			int argmin=0;
			bestsize=0;
			for(int kp=0;kp<PRED_COUNT*2;++kp)
			{
				double size=csizes[kc*PRED_COUNT*2+kp];
				if(!bestsize||bestsize>size)
					bestsize=size, argmin=kp;
			}
			predsel[kc]=argmin;
		}
		bestsize=0;
#ifndef DISABLE_CRCT
		for(int k0=0;k0<(int)_countof(lumapool);++k0)
		{
			int kl=lumapool[k0];
			for(int k1=0;k1<(int)_countof(chroma1pool[0]);++k1)
			{
				int kc1=chroma1pool[k0][k1];
				for(int kc2=0;kc2<OCH_COUNT;++kc2)
				{
					if(kc2/poolsize==kl/poolsize||kc2/poolsize==kc1/poolsize)
						continue;

					double size=
						csizes[kl*PRED_COUNT*2+predsel[kl]]+
						csizes[kc1*PRED_COUNT*2+predsel[kc1]]+
						csizes[kc2*PRED_COUNT*2+predsel[kc2]];
					if(!bestsize||bestsize>size)
						bestsize=size, cidx[0]=kl, cidx[1]=kc1, cidx[2]=kc2;
				}
			}
		}
		if(cidx[0]/poolsize==cidx[1]/poolsize||cidx[0]/poolsize==cidx[2]/poolsize)
			LOG_ERROR("RCT selection error %d %d %d", cidx[0], cidx[1], cidx[2]);
#endif
#ifndef DISABLE_RCT2
		for(int kp=0;kp<PERM_COUNT;++kp)
		{
			for(int kt=0;kt<RCT2_COUNT;++kt)
			{
				int krct=(kt*PERM_COUNT+kp)*3+OCH_COUNT;
				double size=
					csizes[(krct+0)*PRED_COUNT*2+predsel[krct+0]]+
					csizes[(krct+1)*PRED_COUNT*2+predsel[krct+1]]+
					csizes[(krct+2)*PRED_COUNT*2+predsel[krct+2]];
				if(!bestsize||bestsize>size)
					bestsize=size, perm_idx=kp, rct2_idx=kt, use_proper_rct=1;
			}
		}
#endif
		dlist_init(&args->list, 1, 1024, 0);
		gr_enc_init(&ec, &args->list);
		args->bestsize=bestsize;
#ifdef ENABLE_RLE
		args->use_rle=use_rle=args->bestsize<BLOCKDX*BLOCKDY/64;
		gr_enc(&ec, use_rle, 2);
#endif
		gr_enc(&ec, use_proper_rct, 2);
		args->use_proper_rct=use_proper_rct;
		if(use_proper_rct)
		{
			predidx[0]=predsel[(perm_idx*RCT2_COUNT+rct2_idx)*3+0+OCH_COUNT];
			predidx[1]=predsel[(perm_idx*RCT2_COUNT+rct2_idx)*3+1+OCH_COUNT];
			predidx[2]=predsel[(perm_idx*RCT2_COUNT+rct2_idx)*3+2+OCH_COUNT];
			gr_enc(&ec, perm_idx, PERM_COUNT);
			gr_enc(&ec, rct2_idx, RCT2_COUNT);
			memcpy(rgbidx, permutations[perm_idx], sizeof(int[3]));
			args->perm_idx=perm_idx;
			args->rct2_idx=rct2_idx;
		}
		else
		{
			predidx[0]=predsel[cidx[0]];
			predidx[1]=predsel[cidx[1]];
			predidx[2]=predsel[cidx[2]];
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
			case 0x210:perm_idx=0;break;
			case 0x021:perm_idx=1;break;
			case 0x102:perm_idx=2;break;
			case 0x012:perm_idx=3;break;
			case 0x120:perm_idx=4;break;
			case 0x201:perm_idx=5;break;
			default:
				LOG_ERROR("Invalid permutation %d %d %d", rgbidx[0], rgbidx[1], rgbidx[2]);
				break;
			}
			gr_enc(&ec, perm_idx, PERM_COUNT);
			gr_enc(&ec, helper1, 2);
			gr_enc(&ec, alpha1, 17);
			gr_enc(&ec, alpha2, 17);
			args->perm_idx=perm_idx;
			args->helper1=helper1;
			args->alpha1=alpha1;
			args->alpha2=alpha2;
		}
		gr_enc(&ec, predidx[0], PRED_COUNT*2);
		gr_enc(&ec, predidx[1], PRED_COUNT*2);
		gr_enc(&ec, predidx[2], PRED_COUNT*2);
#if 0
		int ystride=image->iw*3;
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
		memset(args->hist, 0, sizeof(args->hist));
		for(int ky=args->y1+1;ky<args->y2;++ky)//analysis loop
		{
			int kx=args->x1+1;
			const short *ptr=image->data+3*(image->iw*ky+kx);
			__m256i amin=_mm256_set1_epi16(-half);
			__m256i amax=_mm256_set1_epi16(half-1);
			__m256i mmask=_mm256_set1_epi16(mask);
			ALIGN(32) short result[16]={0};
			for(;kx<args->x2-5;kx+=5, ptr+=15)
			{
				__m256i pred;
				__m256i comp[3][5]=
				{
					{
						_mm256_loadu_si256((__m256i*)(ptr-1*ystride-1*3+0)),//NW
						_mm256_loadu_si256((__m256i*)(ptr-1*ystride+0*3+0)),//N
						_mm256_loadu_si256((__m256i*)(ptr+0*ystride-1*3+0)),//W
						_mm256_loadu_si256((__m256i*)(ptr+0*ystride+0*3+0)),//curr
					},
					{
						_mm256_loadu_si256((__m256i*)(ptr-1*ystride-1*3+1)),//NW
						_mm256_loadu_si256((__m256i*)(ptr-1*ystride+0*3+1)),//N
						_mm256_loadu_si256((__m256i*)(ptr+0*ystride-1*3+1)),//W
						_mm256_loadu_si256((__m256i*)(ptr+0*ystride+0*3+1)),//curr
					},
					{
						_mm256_loadu_si256((__m256i*)(ptr-1*ystride-1*3+2)),//NW
						_mm256_loadu_si256((__m256i*)(ptr-1*ystride+0*3+2)),//N
						_mm256_loadu_si256((__m256i*)(ptr+0*ystride-1*3+2)),//W
						_mm256_loadu_si256((__m256i*)(ptr+0*ystride+0*3+2)),//curr
					},
				};
#define UPDATE(IDX0, IDX3, IDX6, IDX9, IDXC)\
	do\
	{\
		pred=_mm256_sub_epi16(pred, amin);\
		pred=_mm256_and_si256(pred, mmask);\
		_mm256_store_si256((__m256i*)result, pred);\
		++args->hist[(IDX0)<<8|result[0x0]];\
		++args->hist[(IDX3)<<8|result[0x3]];\
		++args->hist[(IDX6)<<8|result[0x6]];\
		++args->hist[(IDX9)<<8|result[0x9]];\
		++args->hist[(IDXC)<<8|result[0xC]];\
	}while(0)
				for(int kc=0;kc<3;++kc)
				{
					int coffset=kc*(OCH_COUNT/3);
					int kc1=kc+1, kc2=kc+2;
					kc1-=3&-(kc1>=3);
					kc2-=3&-(kc2>=3);
					pred=_mm256_sub_epi16(_mm256_add_epi16(comp[kc][NB_N], comp[kc][NB_W]), comp[kc][NB_NW]);
					pred=_mm256_max_epi16(pred, _mm256_min_epi16(comp[kc][NB_N], comp[kc][NB_W]));
					pred=_mm256_min_epi16(pred, _mm256_max_epi16(comp[kc][NB_N], comp[kc][NB_W]));

					pred=_mm256_sub_epi16(comp[kc][NB_curr], pred);
					UPDATE(
						OCH_Rgb00+coffset,
						OCH_Rgb00+coffset,
						OCH_Rgb00+coffset,
						OCH_Rgb00+coffset,
						OCH_Rgb00+coffset
					);
					for(int k=0;k<=16;++k)
					{
						__m256i helpers[]=
						{
							_mm256_srai_epi16(_mm256_add_epi16(_mm256_mullo_epi16(comp[kc1][0], ramp[k]), _mm256_mullo_epi16(comp[kc2][0], ramp[16-k])), 4),//NW
							_mm256_srai_epi16(_mm256_add_epi16(_mm256_mullo_epi16(comp[kc1][1], ramp[k]), _mm256_mullo_epi16(comp[kc2][1], ramp[16-k])), 4),//N
							_mm256_srai_epi16(_mm256_add_epi16(_mm256_mullo_epi16(comp[kc1][2], ramp[k]), _mm256_mullo_epi16(comp[kc2][2], ramp[16-k])), 4),//W
							_mm256_srai_epi16(_mm256_add_epi16(_mm256_mullo_epi16(comp[kc1][3], ramp[k]), _mm256_mullo_epi16(comp[kc2][3], ramp[16-k])), 4),//curr
						};
						__m256i nb2[]=
						{
							_mm256_sub_epi16(comp[kc][0], helpers[0]),//NW
							_mm256_sub_epi16(comp[kc][1], helpers[1]),//N
							_mm256_sub_epi16(comp[kc][2], helpers[2]),//W
						};
						pred=_mm256_sub_epi16(_mm256_add_epi16(nb2[NB_N], nb2[NB_W]), nb2[NB_NW]);
						pred=_mm256_max_epi16(pred, _mm256_min_epi16(nb2[NB_N], nb2[NB_W]));
						pred=_mm256_min_epi16(pred, _mm256_max_epi16(nb2[NB_N], nb2[NB_W]));

						pred=_mm256_add_epi16(pred, helpers[NB_curr]);
						pred=_mm256_max_epi16(pred, amin);
						pred=_mm256_min_epi16(pred, amax);
						pred=_mm256_sub_epi16(comp[kc][NB_curr], pred);
						UPDATE(
							OCH_Rgb0G+coffset+k,
							OCH_Rgb0G+coffset+k,
							OCH_Rgb0G+coffset+k,
							OCH_Rgb0G+coffset+k,
							OCH_Rgb0G+coffset+k
						);
					}
				}
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
				for(int kc2=0;kc2<OCH_COUNT;++kc2)
				{
					if(kc2/poolsize==kl/poolsize||kc2/poolsize==kc1/poolsize)
						continue;
					double size=csizes[kl]+csizes[kc1]+csizes[kc2];
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
		memset(args->pixels, 0, sizeof(args->pixels));
		res=(args->x2-args->x1)*(args->y2-args->y1);
		for(int ky=args->y1;ky<args->y2;++ky)//analysis loop #2
		{
			int kx=args->x1;
			const short *ptr=image->data+3*(image->iw*ky+kx);
			ALIGN(16) int *rows[]=
			{
				args->pixels+((BLOCKDX+16LL)*((ky-0LL)&3)+8LL)*4*2,
				args->pixels+((BLOCKDX+16LL)*((ky-1LL)&3)+8LL)*4*2,
				args->pixels+((BLOCKDX+16LL)*((ky-2LL)&3)+8LL)*4*2,
				args->pixels+((BLOCKDX+16LL)*((ky-3LL)&3)+8LL)*4*2,
			};
			int yuv[4]={0};
			for(;kx<args->x2;++kx, ptr+=3)
			{
				int offset=0;
				int idx=image->nch*(image->iw*ky+kx);
				int
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
				yuv[0]=image->data[idx+rgbidx[0]];
				yuv[1]=image->data[idx+rgbidx[1]];
				yuv[2]=image->data[idx+rgbidx[2]];
				for(int kc=0;kc<3*2;++kc)
				{
					int preds[PRED_COUNT];
					int kc0=kc>>1, neg=kc&1;
					switch(kc0)
					{
					case 0:
						offset=0;
						break;
					case 1:
						offset=yuv[0]<<PREC_BITS&-helper1;
						break;
					case 2:
						offset=(alpha1*yuv[0]+alpha2*yuv[1])<<PREC_BITS>>4;
						break;
					}
					preds[PRED_W]=W[kc];
					preds[PRED_N]=N[kc];
					MEDIAN3_32(preds[PRED_CG], N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
					CLAMP3_32(preds[PRED_AV5], W[kc]+((5*(N[kc]-NW[kc])+NE[kc]-WW[kc])>>3), N[kc], W[kc], NE[kc]);
					CLAMP3_32(preds[PRED_AV9],
						W[kc]+((10*N[kc]-9*NW[kc]+4*NE[kc]-2*(NN[kc]+WW[kc])+NNW[kc]-(NNE[kc]+NWW[kc]))>>4),
						N[kc],
						W[kc],
						NE[kc]
					);

					preds[PRED_AV12]=(
						av12_icoeffs[ 0]*NNWW[kc]+
						av12_icoeffs[ 1]*NNW[kc]+
						av12_icoeffs[ 2]*NN[kc]+
						av12_icoeffs[ 3]*NNE[kc]+
						av12_icoeffs[ 4]*NNEE[kc]+
						av12_icoeffs[ 5]*NWW[kc]+
						av12_icoeffs[ 6]*NW[kc]+
						av12_icoeffs[ 7]*N[kc]+
						av12_icoeffs[ 8]*NE[kc]+
						av12_icoeffs[ 9]*NEE[kc]+
						av12_icoeffs[10]*WW[kc]+
						av12_icoeffs[11]*W[kc]
					)>>8;
					CLAMP3_32(preds[PRED_AV12], preds[PRED_AV12], N[kc], W[kc], NE[kc]);

					int
						gy=abs(N[kc]-NN[kc])+abs(W[kc]-NW[kc])+1,
						gx=abs(W[kc]-WW[kc])+abs(N[kc]-NW[kc])+1;
					preds[PRED_WG]=(N[kc]*gy+W[kc]*gx)/(gy+gx);
					
					int pixel=yuv[kc0];
					for(int kp=0;kp<PRED_COUNT;++kp)
					{
						int pred=preds[kp];
						if(neg)
							pred=-pred;
						pred+=offset;
						CLAMP2_32(pred, pred, -ehalf, ehalf-1);
						int error=yuv[kc0]-((pred+(1<<PREC_BITS>>1))>>PREC_BITS);
						error>>=image->depth-8;
						error+=128;
						error&=255;
						++args->hist[(kc0*PRED_COUNT+kp)<<9|neg<<8|error];
					}
					curr[kc]=(yuv[kc0]<<PREC_BITS)-offset;
					if(neg)
						curr[kc]=-curr[kc];
				}
				rows[0]+=4*2;
				rows[1]+=4*2;
				rows[2]+=4*2;
				rows[3]+=4*2;
			}
		}
		memset(csizes, 0, sizeof(csizes));
		for(int kc=0;kc<PRED_COUNT*3*2;++kc)
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
			for(int kp=0;kp<PRED_COUNT*2;++kp)
			{
				double size=csizes[kc*PRED_COUNT+kp];
				if(!kp||bestsizes[kc]>size)
					bestsizes[kc]=size, predsel[kc]=kp;
			}
		}
		dlist_init(&args->list, 1, 1024, 0);
		gr_enc_init(&ec, &args->list);
		gr_enc(&ec, permutation, PERM_COUNT);
		gr_enc(&ec, helper1, 2);
		gr_enc(&ec, alpha1, 17);
		gr_enc(&ec, alpha2, 17);
		gr_enc(&ec, predsel[0], PRED_COUNT);
		gr_enc(&ec, predsel[1], PRED_COUNT);
		gr_enc(&ec, predsel[2], PRED_COUNT);
		args->bestsize=bestsizes[0]+bestsizes[1]+bestsizes[2];
#ifdef ENABLE_RLE
		args->use_rle=use_rle=args->bestsize<BLOCKDX*BLOCKDY/64;
		gr_enc(&ec, use_rle, 2);
#endif
		args->permutation=permutation;
		args->helper1=helper1;
		args->alpha1=alpha1;
		args->alpha2=alpha2;
		//memcpy(args->predsel, predsel, sizeof(args->predsel));
#endif
	}
	else
	{
		gr_dec_init(&ec, args->decstart, args->decend);
#ifdef ENABLE_RLE
		use_rle=gr_dec(&ec, 2);
#endif
		use_proper_rct=gr_dec(&ec, 2);
		perm_idx=gr_dec(&ec, PERM_COUNT);
		memcpy(rgbidx, permutations[perm_idx], sizeof(int[3]));
		if(use_proper_rct)
			rct2_idx=gr_dec(&ec, RCT2_COUNT);
		else
		{
			helper1=gr_dec(&ec, 2);
			alpha1=gr_dec(&ec, 17);
			alpha2=gr_dec(&ec, 17);
		}
		predidx[0]=gr_dec(&ec, PRED_COUNT*2);
		predidx[1]=gr_dec(&ec, PRED_COUNT*2);
		predidx[2]=gr_dec(&ec, PRED_COUNT*2);

#if 0
		gr_dec_init(&ec, args->decstart, args->decend);
		permutation=gr_dec(&ec, PERM_COUNT);
		helper1=gr_dec(&ec, 2);
		alpha1=gr_dec(&ec, 17);
		alpha2=gr_dec(&ec, 17);
		memcpy(rgbidx, permutations[permutation], sizeof(int[3]));
		
		predsel[0]=gr_dec(&ec, PRED_COUNT);
		predsel[1]=gr_dec(&ec, PRED_COUNT);
		predsel[2]=gr_dec(&ec, PRED_COUNT);
#ifdef ENABLE_RLE
		use_rle=gr_dec(&ec, 2);
#endif
#endif
	}
	mod_sh[0]=32-image->depth;
	mod_sh[1]=32-image->depth;
	mod_sh[2]=32-image->depth;
	if(use_proper_rct)
	{
		helper1=0;
		alpha1=0;
		alpha2=0;
		--mod_sh[1];
		--mod_sh[2];
	}
	args->nbypass=0;
	args->pneg[0]=pneg[0]=predidx[0]>=PRED_COUNT, args->predidx[0]=predidx[0]%=PRED_COUNT;
	args->pneg[1]=pneg[1]=predidx[1]>=PRED_COUNT, args->predidx[1]=predidx[1]%=PRED_COUNT;
	args->pneg[2]=pneg[2]=predidx[2]>=PRED_COUNT, args->predidx[2]=predidx[2]%=PRED_COUNT;
#ifndef ENABLE_NEG
	pneg[0]=0;
	pneg[1]=0;
	pneg[2]=0;
#endif
	memset(args->pixels, 0, sizeof(args->pixels));
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) int *rows[]=
		{
			args->pixels+((BLOCKDX+16LL)*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKDX+16LL)*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKDX+16LL)*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKDX+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		int yuv[4]={0};
		int pred=0, error=0;
		int eprev[3]={0};
#ifdef ENABLE_RLE
		int run=0;
		int nbypass0=7;
#endif
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int offset=0;
			int idx=image->nch*(image->iw*ky+kx);
			int
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
				yuv[0]=image->data[idx+rgbidx[0]];
				yuv[1]=image->data[idx+rgbidx[1]];
				yuv[2]=image->data[idx+rgbidx[2]];
				if(use_proper_rct)
				{
					switch(rct2_idx)
					{
					case RCT2_JPEG2000:
						yuv[1]-=yuv[0];
						yuv[2]-=yuv[0];
						yuv[0]+=(yuv[1]+yuv[2])>>2;
						break;
					case RCT2_Pei09:
						yuv[1]-=(87*yuv[2]+169*yuv[0]+128)>>8;
						yuv[2]-=yuv[0];
						yuv[0]+=(86*yuv[2]+29*yuv[1]+128)>>8;
						break;
					case RCT2_RCT1:
						yuv[2]-=yuv[0];
						yuv[0]+=yuv[2]>>1;
						yuv[1]-=yuv[0];
						yuv[0]+=yuv[1]>>1;
						break;
					}
				}
			}
			for(int kc=0;kc<image->nch;++kc)
			{
				int ols=0;
				switch(kc)
				{
				case 0:
					offset=0;
					break;
				case 1:
					offset=helper1?yuv[0]<<PREC_BITS:0;
					break;
				case 2:
					offset=(alpha1*yuv[0]+alpha2*yuv[1])<<PREC_BITS>>4;
					break;
				}
				switch(predidx[kc])
				{
				case PRED_N:
					pred=N[kc];
					break;
				case PRED_W:
					pred=W[kc];
					break;
				case PRED_CG:
					MEDIAN3_32(pred, N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
					break;
				case PRED_AV5:
					CLAMP3_32(pred,
						W[kc]+((5*(N[kc]-NW[kc])+NE[kc]-WW[kc])>>3),
						N[kc], W[kc], NE[kc]
					);
					break;
				case PRED_AV9:
					CLAMP3_32(pred,
						W[kc]+((10*N[kc]-9*NW[kc]+4*NE[kc]-2*(NN[kc]+WW[kc])+NNW[kc]-(NNE[kc]+NWW[kc]))>>4),
						N[kc], W[kc], NE[kc]
					);
					break;
				case PRED_AV12:
					pred=(
						av12_icoeffs[ 0]*NNWW[kc]+
						av12_icoeffs[ 1]*NNW[kc]+
						av12_icoeffs[ 2]*NN[kc]+
						av12_icoeffs[ 3]*NNE[kc]+
						av12_icoeffs[ 4]*NNEE[kc]+
						av12_icoeffs[ 5]*NWW[kc]+
						av12_icoeffs[ 6]*NW[kc]+
						av12_icoeffs[ 7]*N[kc]+
						av12_icoeffs[ 8]*NE[kc]+
						av12_icoeffs[ 9]*NEE[kc]+
						av12_icoeffs[10]*WW[kc]+
						av12_icoeffs[11]*W[kc]
					)>>8;
					CLAMP3_32(pred, pred, N[kc], W[kc], NE[kc]);
					break;
				case PRED_WG:
					{
						int
							gy=abs(N[kc]-NN[kc])+abs(W[kc]-NW[kc])+1,
							gx=abs(W[kc]-WW[kc])+abs(N[kc]-NW[kc])+1;
						pred=(N[kc]*gy+W[kc]*gx)/(gy+gx);
					}
					break;
				}
				if(pneg[kc])
					pred=-pred;
				pred+=offset;
				CLAMP2_32(pred, pred, -ehalf, ehalf-1);
				int nbypass=FLOOR_LOG2(W[kc+4]+1)-PREC_BITS;
				if(nbypass<0)
					nbypass=0;
				args->nbypass+=nbypass;
				//if(ky==257&&kx==256&&kc==1)//
				//if(ky==256&&kx==257&&kc==0)//
				//if(ky==511&&kx==123&&kc==0)//
				//if(ky==101&&kx==174&&kc==0)//
				//	printf("");
				if(args->fwd)
				{
					//if(args->blockidx==1)
					//	printf("");
					error=yuv[kc];
					curr[kc]=error<<PREC_BITS;
					error-=(pred+(1<<PREC_BITS>>1))>>PREC_BITS;
					error<<=mod_sh[kc];
					error>>=mod_sh[kc];
					error=error<<1^error>>31;
					//if(error<0)
					//	LOG_ERROR("");
#ifdef ENABLE_RLE
					if(use_rle)
					{
						//if(ky==256&&kx==257&&kc==0)//
						//if(ky==256&&kx==256&&kc==1)//
						//	printf("");
						if(error==eprev[kc])
						{
							++run;
							if(run==1)
								nbypass0=nbypass;
						}
						else
						{
							if(run>0)
							{
								gr_enc_POT(&ec, RLE_ESCAPESYM, nbypass0);//escape symbol
								gr_enc_POT(&ec, run-1, RLE_BITS);
								run=0;
							}
							gr_enc_POT(&ec, error+(error>=RLE_ESCAPESYM), nbypass);
							nbypass0=nbypass;
						}
					}
					else
#endif
						gr_enc_POT(&ec, error, nbypass);
				}
				else
				{
					//if(args->blockidx==1)
					//	printf("");
#ifdef ENABLE_RLE
					if(use_rle)
					{
						//if(ky==256&&kx==257&&kc==0)//
						//if(ky==256&&kx==256&&kc==1)//
						//	printf("");
						if(run>0)
						{
							--run;
							error=eprev[kc];
						}
						else
						{
							error=gr_dec_POT(&ec, nbypass);
							if(error==RLE_ESCAPESYM)
							{
								run=gr_dec_POT(&ec, RLE_BITS)+1;
								--run;
								error=eprev[kc];
							}
							else
								error-=error>=RLE_ESCAPESYM;
						}
					}
					else
#endif
						error=gr_dec_POT(&ec, nbypass);
					int pixel=(error>>1^-(error&1))+((pred+(1<<PREC_BITS>>1))>>PREC_BITS);
					pixel<<=mod_sh[kc];
					pixel>>=mod_sh[kc];
					yuv[kc]=pixel;
					curr[kc]=pixel<<PREC_BITS;
				}
				eprev[kc]=error;
				error=curr[kc]-pred;
				error=error<<1^error>>31;
				curr[kc+4]=(2*W[kc+4]+error+NEEE[kc+4])>>2;
				curr[kc]-=offset;
				if(pneg[kc])
					curr[kc]=-curr[kc];
			}
			if(!args->fwd)
			{
				if(use_proper_rct)
				{
					switch(rct2_idx)
					{
					case RCT2_JPEG2000:
						yuv[0]-=(yuv[1]+yuv[2])>>2;
						yuv[2]+=yuv[0];
						yuv[1]+=yuv[0];
						break;
					case RCT2_Pei09:
						yuv[0]-=(86*yuv[2]+29*yuv[1]+128)>>8;
						yuv[2]+=yuv[0];
						yuv[1]+=(87*yuv[2]+169*yuv[0]+128)>>8;
						break;
					case RCT2_RCT1:
						yuv[0]-=yuv[1]>>1;
						yuv[1]+=yuv[0];
						yuv[0]-=yuv[2]>>1;
						yuv[2]+=yuv[0];
						break;
					}
				}
				args->dst->data[idx+rgbidx[0]]=yuv[0];
				args->dst->data[idx+rgbidx[1]]=yuv[1];
				args->dst->data[idx+rgbidx[2]]=yuv[2];
#ifdef ENABLE_GUIDE
				if(memcmp(image->data+idx, guide->data+idx, sizeof(short)*image->nch))
				{
					short orig[4]={0};
					memcpy(orig, guide->data+idx, image->nch*sizeof(short));
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
#ifdef ENABLE_RLE
		if(args->fwd)
		{
			if(run>0)
			{
				gr_enc_POT(&ec, RLE_ESCAPESYM, nbypass0);//escape symbol
				gr_enc_POT(&ec, run-1, RLE_BITS);
				run=0;
			}
		}
#endif
	}
	if(args->fwd)
		gr_enc_flush(&ec);
}
int f29_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	if(image->nch!=3||image->depth!=8)
	{
		LOG_ERROR("Unsupported  nch %d", image->nch);
		return 1;
	}
	int ncores=query_cpu_cores();
	int
		xblocks=(image->iw+BLOCKDX-1)/BLOCKDX,
		yblocks=(image->ih+BLOCKDY-1)/BLOCKDY,
		nblocks=xblocks*yblocks, nthreads=MINVAR(nblocks, ncores);
	int coffset=sizeof(int)*nblocks;
	ptrdiff_t start=0;
	ptrdiff_t memusage=0;
	ptrdiff_t argssize=nthreads*sizeof(ThreadArgs);
	ThreadArgs *args=(ThreadArgs*)malloc(argssize);
	double esize=0;
	
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	if(fwd)
	{
#ifdef ENABLE_GUIDE
		guide=image;
#endif
		start=array_append(data, 0, 1, coffset, 1, 0, 0);
	}
	else//integrity check
	{
		start=coffset;
		for(int kt=0;kt<nblocks;++kt)
		{
			int size=0;
			memcpy(&size, cbuf+sizeof(int)*kt, sizeof(int));
			start+=size;
		}
		if(start!=(ptrdiff_t)clen)
			LOG_ERROR("Corrupt file");
		start=coffset;
	}

	memusage+=argssize;
	memset(args, 0, argssize);
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=src;
		arg->dst=dst;

		arg->fwd=fwd;
#ifdef DISABLE_MT
		arg->loud=loud;
#else
		arg->loud=0;
#endif
	}
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
			arg->x1=BLOCKDX*kx;
			arg->y1=BLOCKDY*ky;
			arg->x2=MINVAR(arg->x1+BLOCKDX, image->iw);
			arg->y2=MINVAR(arg->y1+BLOCKDY, image->ih);
			if(!fwd)
			{
				int size=0;
				memcpy(&size, cbuf+sizeof(int)*((ptrdiff_t)kt+kt2), sizeof(int));
				arg->decstart=cbuf+start;
				start+=size;
				arg->decend=cbuf+start;
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
				if(loud)
				{
					int blocksize=((arg->x2-arg->x1)*(arg->y2-arg->y1)*image->nch*image->depth+7)>>3;
					int kx, ky;

					kx=kt+kt2;
					ky=kx/xblocks;
					kx%=xblocks;
					if(nblocks<MAXPRINTEDBLOCKS)
					{
						//printf("NB %10.2lf ", arg->nbypass/8.);
						//if(!(kt+kt2))
						//	printf("block,  nrows,  usize,     best  ->  actual,  (actual-best)\n");
						printf(
							"block %4d/%4d  XY %3d %3d  %4d*%4d:  %8d->%16lf->%8zd bytes (%+10.2lf)  %10.6lf%%  CR %10lf  %s",
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
							permnames[arg->perm_idx]
						);
						if(arg->use_proper_rct)
							printf(" %s", rct2names[arg->rct2_idx]);
						else
							printf(" %d (%2d+%2d)/16", arg->helper1, arg->alpha1, arg->alpha2);
						printf(
							"  %c%s %c%s %c%s\n",
							arg->pneg[0]?'-':'+', pred_names[arg->predidx[0]],
							arg->pneg[1]?'-':'+', pred_names[arg->predidx[1]],
							arg->pneg[2]?'-':'+', pred_names[arg->predidx[2]]
						);
					}
					esize+=arg->bestsize;
				}
				memcpy(data[0]->data+start+sizeof(int)*((ptrdiff_t)kt+kt2), &arg->list.nobj, sizeof(int));
				dlist_appendtoarray(&arg->list, data);
				dlist_clear(&arg->list);
			}
		}
	}
	if(loud)
	{
		ptrdiff_t usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
		t0=time_sec()-t0;
		if(fwd)
		{
			ptrdiff_t csize=data[0]->count-start;
			printf("Best %15.2lf (%+13.2lf) bytes\n", esize, csize-esize);
			printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
			printf("Mem usage: ");
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
		}
		printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
	}
	free(args);
	return 0;
}