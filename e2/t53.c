#include"e2.h"
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctype.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;


	#define ENABLE_GUIDE
	#define ENABLE_PREVIEW


//T53 Adaptive block
#define BLOCKSIZE 48
//#define CTXSIZE 64

static int blocks[BLOCKSIZE*BLOCKSIZE*6*4], buf2[BLOCKSIZE*BLOCKSIZE*6*4], buf3[BLOCKSIZE*BLOCKSIZE*6*4];
//static int buf4[CTXSIZE*CTXSIZE*6*4];
typedef enum T53RCTTypeEnum
{
	T53_NONE,
	T53_DUAL_RGP,
	T53_DUAL_RGN,
	T53_DUAL_RBP,
	T53_DUAL_RBN,
	T53_DUAL_GBP,
	T53_DUAL_GBN,
	T53_TRIPLE_RP_BP,
	T53_TRIPLE_RP_BN,
	T53_TRIPLE_RN_BP,
	T53_TRIPLE_RN_BN,
	T53_RCT_COUNT,
} T53RCTType;
typedef enum T53PREDTypeEnum
{
	T53_PRED_ZERO,
	T53_PRED_W,
	T53_PRED_NWW,
	T53_PRED_NW,
	T53_PRED_NNW,
	T53_PRED_N,
	T53_PRED_NNE,
	T53_PRED_NE,
	T53_PRED_GRAD,	//median3(N, W, N+W-NW)
	T53_PRED_AV,	//(N+W)/2
	T53_PRED_AV2,	//(4*(N+W)+NE-NW)/8
	T53_PRED_COUNT,
} T53PREDType;
#define RCT_TRIPLE_FWD(R, G, B) R-=G, G+=R>>1, B-=G, G+=B>>1
#define RCT_TRIPLE_INV(R, G, B) G-=B>>1, B+=G, G-=R>>1, R+=G
#define RCT_DUAL_FWD(A, B) A-=B, B+=A>>1
#define RCT_DUAL_INV(A, B) B-=A>>1, A+=B
static void apply_rct_fwd(int *comp, int rct)
{
	switch(rct)
	{
	case T53_TRIPLE_RP_BP:
		RCT_TRIPLE_FWD(comp[0], comp[1], comp[2]);
		break;
	case T53_TRIPLE_RP_BN:
		comp[2]=-comp[2];
		RCT_TRIPLE_FWD(comp[0], comp[1], comp[2]);
		break;
	case T53_TRIPLE_RN_BP:
		comp[0]=-comp[0];
		RCT_TRIPLE_FWD(comp[0], comp[1], comp[2]);
		break;
	case T53_TRIPLE_RN_BN:
		comp[0]=-comp[0];
		comp[2]=-comp[2];
		RCT_TRIPLE_FWD(comp[0], comp[1], comp[2]);
		break;
	case T53_DUAL_RGP:
		RCT_DUAL_FWD(comp[0], comp[1]);
		break;
	case T53_DUAL_RGN:
		comp[0]=-comp[0];
		RCT_DUAL_FWD(comp[0], comp[1]);
		break;
	case T53_DUAL_RBP:
		RCT_DUAL_FWD(comp[0], comp[2]);
		break;
	case T53_DUAL_RBN:
		comp[0]=-comp[0];
		RCT_DUAL_FWD(comp[0], comp[2]);
		break;
	case T53_DUAL_GBP:
		RCT_DUAL_FWD(comp[1], comp[2]);
		break;
	case T53_DUAL_GBN:
		comp[2]=-comp[2];
		RCT_DUAL_FWD(comp[1], comp[2]);
		break;
	}
}
static void apply_rct_inv(int *comp, int rct)
{
	switch(rct)
	{
	case T53_TRIPLE_RP_BP:
		RCT_TRIPLE_INV(comp[0], comp[1], comp[2]);
		break;
	case T53_TRIPLE_RP_BN:
		RCT_TRIPLE_INV(comp[0], comp[1], comp[2]);
		comp[2]=-comp[2];
		break;
	case T53_TRIPLE_RN_BP:
		RCT_TRIPLE_INV(comp[0], comp[1], comp[2]);
		comp[0]=-comp[0];
		break;
	case T53_TRIPLE_RN_BN:
		RCT_TRIPLE_INV(comp[0], comp[1], comp[2]);
		comp[0]=-comp[0];
		comp[2]=-comp[2];
		break;
	case T53_DUAL_RGP:
		RCT_DUAL_INV(comp[0], comp[1]);
		break;
	case T53_DUAL_RGN:
		RCT_DUAL_INV(comp[0], comp[1]);
		comp[0]=-comp[0];
		break;
	case T53_DUAL_RBP:
		RCT_DUAL_INV(comp[0], comp[2]);
		break;
	case T53_DUAL_RBN:
		RCT_DUAL_INV(comp[0], comp[2]);
		comp[0]=-comp[0];
		break;
	case T53_DUAL_GBP:
		RCT_DUAL_INV(comp[1], comp[2]);
		break;
	case T53_DUAL_GBN:
		RCT_DUAL_INV(comp[1], comp[2]);
		comp[2]=-comp[2];
		break;
	}
}
static int apply_pred(int NW, int N, int NE, int W, int predtype)
{
	switch(predtype)
	{
	case T53_PRED_W:
		return W;
	case T53_PRED_NWW:
		return (W+NW)>>1;
	case T53_PRED_NW:
		return NW;
	case T53_PRED_NNW:
		return (NW+N)>>1;
	case T53_PRED_N:
		return N;
	case T53_PRED_NNE:
		return (N+NE)>>1;
	case T53_PRED_NE:
		return NE;
	case T53_PRED_GRAD:	//median3(N, W, N+W-NW)
		{
			int pred=N+W-NW;
			pred=MEDIAN3(N, W, pred);
			return pred;
		}
	case T53_PRED_AV:	//(N+W)/2
		return (N+W)>>1;
		break;
	case T53_PRED_AV2:	//(4*(N+W)+NE-NW)/8
		return (4*(N+W)+NE-NW)>>3;
	}
	return 0;
}

#define CDFSIZE 127	//power-of-two minus one

#define T53_SYM_EXP 5
#define T53_SYM_MSB 2
#define T53_SYM_LSB 0
typedef struct HybridUintStruct
{
	unsigned short token, nbits;
	unsigned bypass;
} HybridUint;
static void hybriduint_encode(int val, int exp, int msb, int lsb, HybridUint *hu)
{
	int token, bypass, nbits;
	val=val<<1^-(val<0);//pack sign
	if(val<(1<<exp))
	{
		token=val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=floor_log2_32((unsigned)val);
		int mantissa=val-(1<<lgv);
		token = (1<<exp) + (
				(lgv-exp)<<(msb+lsb)|
				(mantissa>>(lgv-msb))<<lsb|
				(mantissa&((1<<lsb)-1))
			);
		nbits=lgv-(msb+lsb);
		bypass=val>>lsb&((1LL<<nbits)-1);
	}
	hu->token=token;
	hu->bypass=bypass;
	hu->nbits=nbits;
}
#define LOAD(X, Y) ((unsigned)(ky2+(Y))<BLOCKSIZE*2&&(unsigned)(kx2+(X))<BLOCKSIZE*3?buf3[(BLOCKSIZE*3*(ky2+(Y))+kx2+(X))<<2|kc]:0)
//#define LOAD2(X, Y) ((unsigned)(ky2+(Y))<CTXSIZE*2&&(unsigned)(kx2+(X))<CTXSIZE*3?buf4[(CTXSIZE*3*(ky2+(Y))+kx2+(X))<<2|kc]:0)
static void prep_block(Image const *src, int kx, int ky, T53RCTType *ret_rct, T53PREDType *ret_pred, unsigned short *hist, unsigned *CDF)
{
	int ytop=MAXVAR(ky-BLOCKSIZE, 0), ybottom=MINVAR(ky+BLOCKSIZE, src->ih);
	int xleft=MAXVAR(kx-BLOCKSIZE, 0), xright=MINVAR(kx+BLOCKSIZE*2, src->iw);
	memset(blocks, 0, sizeof(blocks));
	for(int ky2=ytop;ky2<ybottom;++ky2)
	{
		int xright2=ky2<ky?xright:kx;
		if(xleft<xright2)
			memcpy(blocks+BLOCKSIZE*3*4*(ky2-ky+BLOCKSIZE), src->data+(((size_t)src->iw*ky2+xleft)<<2), ((size_t)xright2-xleft)*sizeof(int[4]));
	}
	ytop+=BLOCKSIZE-ky;
	ybottom+=BLOCKSIZE-ky;
	xleft+=BLOCKSIZE-kx;
	xright+=BLOCKSIZE-kx;
	//int errors[T53_RCT_COUNT*T53_PRED_COUNT]={0};
	int best_rct=0, best_pred=0, besterror=0;
	for(int kt=0;kt<T53_RCT_COUNT;++kt)
	{
		memcpy(buf2, blocks, sizeof(buf2));
		for(int k=0;k<BLOCKSIZE*BLOCKSIZE;++k)
			apply_rct_fwd(buf2+((size_t)k<<2), kt);
		for(int kp=0;kp<T53_PRED_COUNT;++kp)
		{
			int error=0;
			memcpy(buf3, buf2, sizeof(buf3));
			for(int kc=0;kc<3;++kc)
			{
				for(int ky2=ytop;ky2<ybottom;++ky2)
				{
					for(int kx2=xleft;kx2<xright;++kx2)
					{
						if(ky2>=BLOCKSIZE&&kx2>=BLOCKSIZE)
							break;
						int
							NW	=LOAD(-1, -1),
							N	=LOAD( 0, -1),
							NE	=LOAD( 1, -1),
							W	=LOAD(-1,  0);
						int pred=apply_pred(NW, N, NE, W, kp);
						int curr=buf3[(BLOCKSIZE*3*ky2+kx2)<<2|kc];
						error+=abs(curr-pred);
					}
				}
			}
			if(!kt&&!kp||besterror>error)
				besterror=error, best_rct=kt, best_pred=kp;
			//errors[T53_PRED_COUNT*kt+kp]=error;
		}
	}
	*ret_rct=best_rct;
	*ret_pred=best_pred;

	memset(hist, 0, sizeof(short[(CDFSIZE+1)*4]));
	memcpy(buf3, blocks, sizeof(buf2));
	for(int k=0;k<BLOCKSIZE*BLOCKSIZE;++k)
		apply_rct_fwd(buf2+((size_t)k<<2), best_rct);
	for(int kc=0;kc<3;++kc)
	{
		unsigned short *curr_hist=hist+(CDFSIZE+1)*kc;
		unsigned *curr_CDF=CDF+(CDFSIZE+1)*kc;
		for(int ky2=ytop;ky2<ybottom;++ky2)
		{
			for(int kx2=xleft;kx2<xright;++kx2)
			{
				if(ky2>=BLOCKSIZE&&kx2>=BLOCKSIZE)
					break;
				int
					NW	=LOAD(-1, -1),
					N	=LOAD( 0, -1),
					NE	=LOAD( 1, -1),
					W	=LOAD(-1,  0);
				int pred=apply_pred(NW, N, NE, W, best_pred);
				int curr=buf3[(BLOCKSIZE*3*ky2+kx2)<<2|kc];
				HybridUint hu;
				hybriduint_encode(curr-pred, T53_SYM_EXP, T53_SYM_MSB, T53_SYM_LSB, &hu);
				if(hu.token>=CDFSIZE)
					LOG_ERROR("");
				++curr_hist[hu.token];
				++curr_hist[CDFSIZE];
			}
		}
		int sum=curr_hist[CDFSIZE];
		if(sum)
		{
			long long c=0;
			for(int ks=0;ks<CDFSIZE;++ks)
			{
				long long freq=curr_hist[ks];
				curr_CDF[ks]=(int)(c*((1LL<<16)-CDFSIZE)/sum)+ks;
				if(curr_CDF[ks]>0x10000)
					LOG_ERROR("");
				c+=freq;
			}
			curr_CDF[CDFSIZE]=1<<16;
		}
		else
		//if(!curr_CDF[CDFSIZE])
		{
			for(int ks=0;ks<CDFSIZE;++ks)
			{
				int freq=0x10000/(ks+1);
				curr_CDF[ks]=freq;
				sum+=freq;
			}
			for(int ks=0, c=0;ks<CDFSIZE;++ks)
			{
				int freq=curr_CDF[ks];
				curr_CDF[ks]=(int)((long long)c*((1LL<<16)-CDFSIZE)/sum)+ks;
				if(curr_CDF[ks]>0x10000)
					LOG_ERROR("");
				c+=freq;
			}
			curr_CDF[CDFSIZE]=1<<16;
		}
	}
}
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
int t53_encode(Image const *src, ArrayHandle *data, int loud)
{
#ifdef ENABLE_GUIDE
	guide=src;
#endif
	double t_start=time_sec();
	int nblocksx=(src->iw+BLOCKSIZE-1)/BLOCKSIZE;
	int nblocksy=(src->ih+BLOCKSIZE-1)/BLOCKSIZE;
	int nblocks=nblocksx*nblocksy;
	unsigned short *hist=(unsigned short*)malloc(sizeof(short[(CDFSIZE+1)*4]));//cdfsize=127, last element is histsum
	unsigned *CDF=(unsigned*)malloc(sizeof(int[(CDFSIZE+1)*4]));
	//char *RCTs=(char*)malloc(nblocks*sizeof(char));
#ifdef ENABLE_PREVIEW
	char *preview_rct=(char*)malloc(nblocks*sizeof(char));//
	char *preview_pred=(char*)malloc(nblocks*sizeof(char));//
#endif
	if(!hist||!CDF
#ifdef ENABLE_PREVIEW
		||!preview_rct||!preview_pred
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(CDF, 0, sizeof(int[(CDFSIZE+1)*4]));
	ArithmeticCoder ec;
	DList list;
	dlist_init(&list, 1, 0x10000, 0);
	ac_enc_init(&ec, &list);
	for(int ky=0;ky<src->ih;ky+=BLOCKSIZE)
	{
		int yend=MINVAR(ky+BLOCKSIZE, src->ih);
		for(int kx=0;kx<src->iw;kx+=BLOCKSIZE)
		{
			//if(kx==0&&ky==BLOCKSIZE)//
			//	printf("");

			int xend=MINVAR(kx+BLOCKSIZE, src->iw);
			T53RCTType rct=0;
			T53PREDType predtype=0;
			prep_block(src, kx, ky, &rct, &predtype, hist, CDF);
			
#ifdef ENABLE_PREVIEW
			preview_rct[nblocksx*(ky/BLOCKSIZE)+kx/BLOCKSIZE]=rct<<4;//
			preview_pred[nblocksx*(ky/BLOCKSIZE)+kx/BLOCKSIZE]=predtype<<4;//
#endif

			for(int ky2=BLOCKSIZE, ky3=ky;ky3<yend;++ky2, ++ky3)
			{
				for(int kx2=BLOCKSIZE, kx3=kx;kx3<xend;++kx2, ++kx3)
				{
					int disable_NE=kx3==xend-1&&ky3>ky;
					int *curr=buf3+(((size_t)BLOCKSIZE*3*ky2+kx2)<<2);
					memcpy(curr, src->data+(((size_t)src->iw*ky3+kx3)<<2), sizeof(int[4]));
					apply_rct_fwd(curr, rct);
					for(int kc=0;kc<src->nch;++kc)
					{
						int
							NW	=LOAD(-1, -1),
							N	=LOAD( 0, -1),
							NE	=disable_NE?0:LOAD( 1, -1),
							W	=LOAD(-1,  0);
						int pred=apply_pred(NW, N, NE, W, predtype);

						HybridUint hu;
						hybriduint_encode(curr[kc]-pred, T53_SYM_EXP, T53_SYM_MSB, T53_SYM_LSB, &hu);
						ac_enc(&ec, hu.token, CDF+(CDFSIZE+1)*kc, CDFSIZE, 0);
						if(hu.nbits)
						{
							int bypass=hu.bypass, nbits=hu.nbits;
							while(nbits>8)
							{
								ac_enc(&ec, bypass>>(nbits-8)&0xFF, 0, 1<<8, 16-8);
								nbits-=8;
							}
							ac_enc(&ec, bypass&((1<<nbits)-1), 0, 1<<nbits, 16-nbits);
						}
					}
				}
			}
		}
	}
#ifdef ENABLE_PREVIEW
	save_mono8("preview_rct.PNG", preview_rct, nblocksx, nblocksy, 1);
	save_mono8("preview_pred.PNG", preview_pred, nblocksx, nblocksy, 1);
	free(preview_rct);
	free(preview_pred);
#endif
	ac_enc_flush(&ec);
	dlist_appendtoarray(&list, data);
	if(loud)
	{
		double usize=image_getBMPsize(src);
		printf("\n");
		printf("Enc ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");

		printf("csize %8d  invCR %10.6lf%%  usize %.0lf\n", (int)list.nobj, 100.*list.nobj/usize, usize);
	}
	dlist_clear(&list);
	free(hist);
	free(CDF);
	return 1;
}
int t53_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	double t_start=time_sec();
	//int nblocksx=(dst->iw+BLOCKSIZE-1)/BLOCKSIZE;
	//int nblocksy=(dst->ih+BLOCKSIZE-1)/BLOCKSIZE;
	//int nblocks=nblocksx*nblocksy;
	unsigned short *hist=(unsigned short*)malloc(sizeof(short[(CDFSIZE+1)*4]));//cdfsize=127, last element is histsum
	unsigned *CDF=(unsigned*)malloc(sizeof(int[(CDFSIZE+1)*4]));
	//char *RCTs=(char*)malloc(nblocks*sizeof(char));
	if(!hist||!CDF)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(CDF, 0, sizeof(int[(CDFSIZE+1)*4]));
	ArithmeticCoder ec;
	ac_dec_init(&ec, data, data+srclen);
	for(int ky=0;ky<dst->ih;ky+=BLOCKSIZE)
	{
		int yend=MINVAR(ky+BLOCKSIZE, dst->ih);
		for(int kx=0;kx<dst->iw;kx+=BLOCKSIZE)
		{
			//if(kx==0&&ky==BLOCKSIZE)//
			//	printf("");

			int xend=MINVAR(kx+BLOCKSIZE, dst->iw);
			T53RCTType rct=0;
			T53PREDType predtype=0;
			prep_block(dst, kx, ky, &rct, &predtype, hist, CDF);
			for(int ky2=BLOCKSIZE, ky3=ky;ky3<yend;++ky2, ++ky3)
			{
				for(int kx2=BLOCKSIZE, kx3=kx;kx3<xend;++kx2, ++kx3)
				{
					int disable_NE=kx3==xend-1&&ky3>ky;
					int *curr=buf3+(((size_t)BLOCKSIZE*3*ky2+kx2)<<2);
					for(int kc=0;kc<dst->nch;++kc)
					{
						int
							NW	=LOAD(-1, -1),
							N	=LOAD( 0, -1),
							NE	=disable_NE?0:LOAD( 1, -1),
							W	=LOAD(-1,  0);
						int pred=apply_pred(NW, N, NE, W, predtype);
						
						int token=ac_dec(&ec, CDF+(CDFSIZE+1)*kc, CDFSIZE, 0);
						int error=token;
						if(error>=(1<<T53_SYM_EXP))
						{
							error-=1<<T53_SYM_EXP;
							int lsb=error&((1<<T53_SYM_LSB)-1);
							error>>=T53_SYM_LSB;
							int msb=error&((1<<T53_SYM_MSB)-1);
							error>>=T53_SYM_MSB;
							int nbits=error+T53_SYM_EXP-(T53_SYM_MSB+T53_SYM_LSB), n=nbits;
							int bypass=0;
							while(n>8)
							{
								n-=8;
								bypass|=ac_dec(&ec, 0, 1<<8, 16-8)<<n;
							}
							bypass|=ac_dec(&ec, 0, 1<<n, 16-n);
							error=1;
							error<<=T53_SYM_MSB;
							error|=msb;
							error<<=nbits;
							error|=bypass;
							error<<=T53_SYM_LSB;
							error|=lsb;
						}
						error=error>>1^-(error&1);
						curr[kc]=error+pred;
					}
					int *currdst=dst->data+(((size_t)dst->iw*ky3+kx3)<<2);
					memcpy(currdst, curr, sizeof(int[4]));
					apply_rct_inv(currdst, rct);
#ifdef ENABLE_GUIDE
					if(guide&&memcmp(currdst, guide->data+(((size_t)dst->iw*ky3+kx3)<<2), dst->nch*sizeof(int)))
						LOG_ERROR("Guide error  XY %d %d", kx3, ky3);
#endif
				}
			}
		}
	}
	if(loud)
	{
		printf("\n");
		printf("Dec ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	free(hist);
	free(CDF);
	return 1;
}
//#undef  LOAD