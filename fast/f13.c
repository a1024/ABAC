#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


	#define NCTX 2
	#define RESCALE_PERIOD	0x1000
//	#define DEBUG_CHECKS
//	#define PROFILER 1


#ifdef PROFILER
#define CHECKPOINTLIST\
	CHECKPOINT(PREP)\
	CHECKPOINT(CTX)\
	CHECKPOINT(TOKEN)\
	CHECKPOINT(CODE)\
	CHECKPOINT(STATS)
#endif
#include"ac.h"
#include"profiler.h"

#ifdef DEBUG_CHECKS
static short gbuf[65536]={0};
static int global_idx=0;
#endif

static int CG(int N, int W, int NW)
{
	int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
	int pred=N+W-NW;
	pred=CLAMP(vmin, pred, vmax);
	return pred;
}
static int clamp(int x, int half)
{
	return CLAMP(-half, x, half-1);
}
#define NPREDS 3	//preds {0, (N+W)/2, CG} * {0, subtract previous channel}
#define NRCH (NPREDS*3)
#define NGCH NPREDS
#define NBCH (NPREDS*2)
#define NCH (NRCH+NGCH+NBCH)
#define CHLIST\
	CHANNEL(g)\
	CHANNEL(g-((gN+gW)>>1))\
	CHANNEL(g-CG(gN, gW, gNW))\
	CHANNEL(b)\
	CHANNEL(b-((bN+bW)>>1))\
	CHANNEL(b-CG(bN, bW, bNW))\
	CHANNEL(b-g)\
	CHANNEL(b-clamp(((bN-gN+bW-gW)>>1)+g, half))\
	CHANNEL(b-clamp(CG(bN-gN, bW-gW, bNW-gNW)+g, half))\
	CHANNEL(r)\
	CHANNEL(r-((rN+rW)>>1))\
	CHANNEL(r-CG(rN, rW, rNW))\
	CHANNEL(r-g)\
	CHANNEL(r-clamp(((rN-gN+rW-gW)>>1)+g, half))\
	CHANNEL(r-clamp(CG(rN-gN, rW-gW, rNW-gNW)+g, half))\
	CHANNEL(r-((g+b)>>1))\
	CHANNEL(r-clamp(((2*(rN+rW)-(gN+gW+bN+bW))>>2)+((g+b)>>1), half))\
	CHANNEL(r-clamp((CG(2*rN-(gN+bN), 2*rW-(gW+bW), 2*rNW-(gNW+bNW))+g+b)>>1, half))
static const char *chnames[]=
{
#define CHANNEL(X) #X,
	CHLIST
#undef  CHANNEL
};
//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
#define CONFIG_EXP 5
#define CONFIG_MSB 2
#define CONFIG_LSB 0
static int quantize_pixel(int pixel, int *bypass, int *nbits)
{
	int token;
	if(pixel<(1<<CONFIG_EXP))
	{
		token=pixel;
		*nbits=0;
		*bypass=0;
	}
	else
	{
		int lgv=floor_log2_32((unsigned)pixel);
		int mantissa=pixel-(1<<lgv);
		token = (1<<CONFIG_EXP) + (
				(lgv-CONFIG_EXP)<<(CONFIG_MSB+CONFIG_LSB)|
				(mantissa>>(lgv-CONFIG_MSB))<<CONFIG_LSB|
				(mantissa&((1<<CONFIG_LSB)-1))
			);
		*nbits=lgv-CONFIG_MSB+CONFIG_LSB;
		*bypass=pixel>>CONFIG_LSB&((1LL<<*nbits)-1);
	}
	return token;
}
static void update_stats(const int *hist, unsigned *CDF, int nqlevels)
{
	int sum=hist[nqlevels];
	for(int ks=0, c=0;ks<nqlevels;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*(0x10000LL-nqlevels)/sum)+ks;
		//CDF[ks]=c*(0x10000-nqlevels)/sum+ks;//X  init_stats passes large hist
		c+=freq;
	}
}
static void init_stats(int *hist, unsigned *CDF, int nclevels, int nqlevels)
{
	int sum=0;
	for(int ks=0;ks<nqlevels;++ks)
	{
		int val=nqlevels-ks;
		sum+=hist[ks]=val*val;
	}
	hist[nqlevels]=sum;
	update_stats(hist, CDF, nqlevels);
	CDF[nqlevels]=0x10000;
	memfill(hist+nqlevels+1, hist, sizeof(int)*((size_t)nclevels*NCTX-1LL)*(nqlevels+1LL), sizeof(int)*(nqlevels+1LL));
	memfill(CDF +nqlevels+1, CDF , sizeof(int)*((size_t)nclevels*NCTX-1LL)*(nqlevels+1LL), sizeof(int)*(nqlevels+1LL));
	//*alpha=0x4000;
	//memfill(alpha+1, alpha, (nclevels-1LL)*sizeof(short), sizeof(short));

	//*hist=1;
	//memfill(hist+1, hist, sizeof(int)*(nqlevels-1LL), sizeof(int));
	//hist[nqlevels]=nqlevels;
	//for(int k=0;k<nqlevels+1;++k)
	//	CDF[k]=(k<<16)/nqlevels;
}
static void restale_hist(int *hist, int nqlevels)
{
	int sum=0;
	for(int ks=0;ks<nqlevels;++ks)
		sum+=hist[ks]=(hist[ks]+1)>>1;
	hist[nqlevels]=sum;
}
#if 1
static int quantize_ctx(int x)
{
	int negmask=-(x<0);
	x=abs(x);
	x=floor_log2_32(x)+1;
	//x=floor_log2_32(x)+1;
	//x=floor_log2_32(x)+1;
	x^=negmask;
	x-=negmask;
	return x;
}
#define QUANTIZE(X) quantize_ctx(X)+(nclevels>>1)
//#define QUANTIZE(X) (((X)>4)-((X)<-4)+1)
//#define QUANTIZE(X) (THREEWAY(X, 0)+1)	//X  inefficient
#endif
static void init_rows(short **rows, short *errors, int iw, int ky)
{
	for(int k=0;k<4;++k)
		rows[k]=errors+((iw+2LL)*(((size_t)ky-k)&(4-1LL))+1);
}
#if 0
static void enc2(int pixel, int pred, int nlevels, short **rows, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nqlevels)
{
#ifdef DEBUG_CHECKS
	static ctr=0;//
	++ctr;
#endif
	int
		eNW	=rows[1][-1],
		eN	=rows[1][+0],
		eNE	=rows[1][+1],
		eW	=rows[0][-1];

	int ctx=QUANTIZE(eN);
	ctx=3*ctx+QUANTIZE(eW);
	ctx=3*ctx+QUANTIZE(eNW);
	ctx=3*ctx+QUANTIZE(eNE);
#ifdef DEBUG_CHECKS
	//if(ctr==15451)
	//	printf("ctx %d  nb %d %d %d %d\n", ctx, eN, eW, eNW, eNE);
	if((unsigned)ctx>=3*3*3*3)
		LOG_ERROR("ctx %d", ctx);
#endif
	ctx*=nqlevels+1;
	unsigned *curr_CDF=CDF+ctx;
	int *curr_hist=hist+ctx;
	PROF(CTX);

	int delta=pixel-pred;
	rows[0][0]=(short)delta;
	delta+=nlevels>>1;
	delta&=nlevels-1;
	delta-=nlevels>>1;
	delta=delta<<1^-(delta<0);

	int token, bypass, nbits;
	token=quantize_pixel(delta, &bypass, &nbits);
#if 0
	if(delta<(1<<CONFIG_EXP))
	{
		token=delta;
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=floor_log2_32((unsigned)delta);
		int mantissa=delta-(1<<lgv);
		token = (1<<CONFIG_EXP) + (
				(lgv-CONFIG_EXP)<<(CONFIG_MSB+CONFIG_LSB)|
				(mantissa>>(lgv-CONFIG_MSB))<<CONFIG_LSB|
				(mantissa&((1<<CONFIG_LSB)-1))
			);
		nbits=lgv-CONFIG_MSB+CONFIG_LSB;
		bypass=delta>>CONFIG_LSB&((1LL<<nbits)-1);
	}
#endif
	PROF(TOKEN);

	ac_enc(ec, token, curr_CDF);
	if(nbits)
		ac_enc_bypass(ec, bypass, 1<<nbits);
	PROF(CODE);
#ifdef DEBUG_CHECKS
	if(token>=nqlevels)
		LOG_ERROR("enc2: token %d/%d", token, nqlevels);
#endif
	++curr_hist[token];
	++curr_hist[nqlevels];
	if(curr_hist[nqlevels]>=RESCALE_PERIOD)
	{
		update_stats(curr_hist, curr_CDF, nqlevels);
		restale_hist(curr_hist, nqlevels);
	}
	PROF(STATS);
}
#endif
#define F13_ENC2
static void enc_self0(Image const *src, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels)
{
	short *rows[4]={0};
	for(int ky=0, idx=kc;ky<src->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, src->iw, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=3)
		{
			const int pred=0;
#include"f13.h"
			//enc2(src->data[idx], 0, nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void enc_self2(Image const *src, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels)
{
	short *rows[4]={0};
	int rowstride=src->iw*3;
	for(int ky=0, idx=kc;ky<src->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, src->iw, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=3)
		{
			//if(ky==0&&kx==766)//
			//if(ky==5&&kx==671)//
			//if(ky==4&&kx==671)//
			//if(ky==10&&kx==164)//
			//if(ky==1&&kx==77)//
			//	printf("");

			int
				N	=ky?src->data[idx-rowstride]:0,
				W	=kx?src->data[idx-3]:0;
			int pred=(N+W)>>1;
#include"f13.h"
			//enc2(src->data[idx], (N+W)>>1, nlevels, rows, ec, hist, CDF, nqlevels);
#ifdef DEBUG_CHECKS
			if(global_idx<0x10000)
				gbuf[global_idx++]=src->data[idx];
#endif

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void enc_self3(Image const *src, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels)
{
	short *rows[4]={0};
	int rowstride=src->iw*3;
	for(int ky=0, idx=kc;ky<src->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, src->iw, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=3)
		{
			int
				NW	=kx&ky?src->data[idx-rowstride-3]:0,
				N	=ky?src->data[idx-rowstride]:0,
				W	=kx?src->data[idx-3]:0;
			int pred=CG(N, W, NW);
#include"f13.h"
			//enc2(src->data[idx], CG(N, W, NW), nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void enc_aux0 (Image const *src, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels, int kc_offset)
{
	short *rows[4]={0};
	for(int ky=0, idx=kc;ky<src->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, src->iw, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=3)
		{
			int pred=src->data[idx+kc_offset];
#include"f13.h"
			//enc2(src->data[idx], src->data[idx+kc_offset], nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void enc_aux2 (Image const *src, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels, int kc_offset)
{
	short *rows[4]={0};
	int rowstride=src->iw*3, half=nlevels>>1;
	for(int ky=0, idx=kc;ky<src->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, src->iw, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=3)
		{
			int
				cN	=ky?src->data[idx-rowstride]:0,
				pN	=ky?src->data[idx-rowstride+kc_offset]:0,
				cW	=kx?src->data[idx-3]:0,
				pW	=kx?src->data[idx-3+kc_offset]:0,
				p	=src->data[idx+kc_offset];
			int pred=clamp(((cN-pN+cW-pW)>>1)+p, half);
#include"f13.h"
			//enc2(src->data[idx], clamp(((cN-pN+cW-pW)>>1)+p, half), nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void enc_aux3 (Image const *src, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels, int kc_offset)
{
	short *rows[4]={0};
	int rowstride=src->iw*3, half=nlevels>>1;
	for(int ky=0, idx=kc;ky<src->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, src->iw, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=3)
		{
			int
				cNW	=ky&kx?src->data[idx-rowstride-3]:0,
				pNW	=ky&kx?src->data[idx-rowstride-3+kc_offset]:0,
				cN	=ky?src->data[idx-rowstride]:0,
				pN	=ky?src->data[idx-rowstride+kc_offset]:0,
				cW	=kx?src->data[idx-3]:0,
				pW	=kx?src->data[idx-3+kc_offset]:0,
				p	=src->data[idx+kc_offset];
			int pred=clamp(CG(cN-pN, cW-pW, cNW-pNW)+p, half);
#include"f13.h"
			//enc2(src->data[idx], clamp(CG(cN-pN, cW-pW, cNW-pNW)+p, half), nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void enc_red0 (Image const *src, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels, int g_offset, int b_offset)
{
	short *rows[4]={0};
	for(int ky=0, idx=kc;ky<src->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, src->iw, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=3)
		{
			int pred=(src->data[idx+g_offset]+src->data[idx+b_offset])>>1;
#include"f13.h"
			//enc2(src->data[idx], src->data[idx+kc_offset], nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void enc_red2 (Image const *src, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels, int g_offset, int b_offset)
{
	short *rows[4]={0};
	int rowstride=src->iw*3, half=nlevels>>1;
	for(int ky=0, idx=kc;ky<src->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, src->iw, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=3)
		{
			int
				rN	=ky?src->data[idx-rowstride]:0,
				gN	=ky?src->data[idx-rowstride+g_offset]:0,
				bN	=ky?src->data[idx-rowstride+b_offset]:0,
				rW	=kx?src->data[idx-3]:0,
				gW	=kx?src->data[idx-3+g_offset]:0,
				bW	=kx?src->data[idx-3+b_offset]:0,
				g	=src->data[idx+g_offset],
				b	=src->data[idx+b_offset];
			int pred=clamp(((2*(rN+rW)-(gN+gW+bN+bW))>>2)+((g+b)>>1), half);
#include"f13.h"
			//enc2(src->data[idx], clamp(((cN-pN+cW-pW)>>1)+p, half), nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void enc_red3 (Image const *src, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels, int g_offset, int b_offset)
{
	short *rows[4]={0};
	int rowstride=src->iw*3, half=nlevels>>1;
	for(int ky=0, idx=kc;ky<src->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, src->iw, ky);
		for(int kx=0;kx<src->iw;++kx, idx+=3)
		{
			int
				rNW	=ky&kx?src->data[idx-rowstride-3]:0,
				gNW	=ky&kx?src->data[idx-rowstride-3+g_offset]:0,
				bNW	=ky&kx?src->data[idx-rowstride-3+b_offset]:0,
				rN	=ky?src->data[idx-rowstride]:0,
				gN	=ky?src->data[idx-rowstride+g_offset]:0,
				bN	=ky?src->data[idx-rowstride+b_offset]:0,
				rW	=kx?src->data[idx-3]:0,
				gW	=kx?src->data[idx-3+g_offset]:0,
				bW	=kx?src->data[idx-3+b_offset]:0,
				g	=src->data[idx+g_offset],
				b	=src->data[idx+b_offset];
			int pred=clamp((CG(2*rN-(gN+bN), 2*rW-(gW+bW), 2*rNW-(gNW+bNW))+g+b)>>1, half);
#include"f13.h"
			//enc2(src->data[idx], clamp(CG(cN-pN, cW-pW, cNW-pNW)+p, half), nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
#undef  F13_ENC2
#if 0
static int dec2(int pred, int nlevels, short **rows, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nqlevels)
{
#ifdef DEBUG_CHECKS
	static ctr=0;//
	++ctr;
#endif
	int
		eNW	=rows[1][-1],
		eN	=rows[1][+0],
		eNE	=rows[1][+1],
		eW	=rows[0][-1];
	
	int ctx=QUANTIZE(eN);
	ctx=3*ctx+QUANTIZE(eW);
	ctx=3*ctx+QUANTIZE(eNW);
	ctx=3*ctx+QUANTIZE(eNE);
#ifdef DEBUG_CHECKS
	//if(ctr==15451)
	//	printf("ctx %d  nb %d %d %d %d\n", ctx, eN, eW, eNW, eNE);
#endif
	ctx*=nqlevels+1;
	unsigned *curr_CDF=CDF+ctx;
	int *curr_hist=hist+ctx;
	PROF(CTX);

	int token=ac_dec_packedsign(ec, curr_CDF, nqlevels), bypass, nbits;
	PROF(CODE);
	int delta=token;
	if(delta>=(1<<CONFIG_EXP))
	{
		delta-=1<<CONFIG_EXP;
		int lsb=delta&((1<<CONFIG_LSB)-1);
		delta>>=CONFIG_LSB;
		int msb=delta&((1<<CONFIG_MSB)-1);
		delta>>=CONFIG_MSB;
		nbits=delta+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB);
		bypass=ac_dec_bypass(ec, 1<<nbits);
		delta=1;
		delta<<=CONFIG_MSB;
		delta|=msb;
		delta<<=nbits;
		delta|=bypass;
		delta<<=CONFIG_LSB;
		delta|=lsb;
	}
	delta=delta>>1^-(delta&1);
	delta+=pred;
	delta+=nlevels>>1;
	delta&=nlevels-1;
	delta-=nlevels>>1;
	rows[0][0]=(short)(delta-pred);
	PROF(TOKEN);
#ifdef DEBUG_CHECKS
	if(token>=nqlevels)
		LOG_ERROR("dec2: token %d/%d  ctr %d", token, nqlevels, ctr);
#endif
	++curr_hist[token];
	++curr_hist[nqlevels];
	if(curr_hist[nqlevels]>=RESCALE_PERIOD)
	{
		update_stats(curr_hist, curr_CDF, nqlevels);
		restale_hist(curr_hist, nqlevels);
	}
	PROF(STATS);
	return delta;
}
#endif
#define F13_DEC2
static void dec_self0(Image *dst, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels)
{
	short *rows[4]={0};
	for(int ky=0, idx=kc;ky<dst->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, dst->iw, ky);
		for(int kx=0;kx<dst->iw;++kx, idx+=3)
		{
			const int pred=0;
#include"f13.h"
			dst->data[idx]=(short)delta;
			//dst->data[idx]=(short)dec2(0, nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void dec_self2(Image *dst, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels)
{
	short *rows[4]={0};
	int rowstride=dst->iw*3;
	for(int ky=0, idx=kc;ky<dst->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, dst->iw, ky);
		for(int kx=0;kx<dst->iw;++kx, idx+=3)
		{
			//if(ky==0&&kx==766)//
			//if(ky==5&&kx==671)//
			//if(ky==4&&kx==671)//
			//if(ky==10&&kx==164)//
			//if(ky==1&&kx==77)//
			//	printf("");

			int
				N	=ky?dst->data[idx-rowstride]:0,
				W	=kx?dst->data[idx-3]:0;
			int pred=(N+W)>>1;
#include"f13.h"
			dst->data[idx]=(short)delta;
			//dst->data[idx]=(short)dec2((N+W)>>1, nlevels, rows, ec, hist, CDF, nqlevels);
#ifdef DEBUG_CHECKS
	if(global_idx<0x10000)
	{
		if(dst->data[idx]!=gbuf[global_idx])
			LOG_ERROR("Guide error idx %d", global_idx);
		++global_idx;
	}
#endif

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void dec_self3(Image *dst, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels)
{
	short *rows[4]={0};
	int rowstride=dst->iw*3;
	for(int ky=0, idx=kc;ky<dst->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, dst->iw, ky);
		for(int kx=0;kx<dst->iw;++kx, idx+=3)
		{
			int
				NW	=kx&ky?dst->data[idx-rowstride-3]:0,
				N	=ky?dst->data[idx-rowstride]:0,
				W	=kx?dst->data[idx-3]:0;
			int pred=CG(N, W, NW);
#include"f13.h"
			dst->data[idx]=(short)delta;
			//dst->data[idx]=(short)dec2(CG(N, W, NW), nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void dec_aux0 (Image *dst, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels, int kc_offset)
{
	short *rows[4]={0};
	for(int ky=0, idx=kc;ky<dst->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, dst->iw, ky);
		for(int kx=0;kx<dst->iw;++kx, idx+=3)
		{
			int pred=dst->data[idx+kc_offset];
#include"f13.h"
			dst->data[idx]=(short)delta;
			//dst->data[idx]=(short)dec2(dst->data[idx+kc_offset], nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void dec_aux2 (Image *dst, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels, int kc_offset)
{
	short *rows[4]={0};
	int rowstride=dst->iw*3, half=nlevels>>1;
	for(int ky=0, idx=kc;ky<dst->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, dst->iw, ky);
		for(int kx=0;kx<dst->iw;++kx, idx+=3)
		{
			int
				cN	=ky?dst->data[idx-rowstride]:0,
				pN	=ky?dst->data[idx-rowstride+kc_offset]:0,
				cW	=kx?dst->data[idx-3]:0,
				pW	=kx?dst->data[idx-3+kc_offset]:0,
				p	=dst->data[idx+kc_offset];
			int pred=clamp(((cN-pN+cW-pW)>>1)+p, half);
#include"f13.h"
			dst->data[idx]=(short)delta;
			//dst->data[idx]=(short)dec2(clamp(((cN-pN+cW-pW)>>1)+p, half), nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void dec_aux3 (Image *dst, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels, int kc_offset)
{
	short *rows[4]={0};
	int rowstride=dst->iw*3, half=nlevels>>1;
	for(int ky=0, idx=kc;ky<dst->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, dst->iw, ky);
		for(int kx=0;kx<dst->iw;++kx, idx+=3)
		{
			int
				cNW	=ky&kx?dst->data[idx-rowstride-3]:0,
				pNW	=ky&kx?dst->data[idx-rowstride-3+kc_offset]:0,
				cN	=ky?dst->data[idx-rowstride]:0,
				pN	=ky?dst->data[idx-rowstride+kc_offset]:0,
				cW	=kx?dst->data[idx-3]:0,
				pW	=kx?dst->data[idx-3+kc_offset]:0,
				p	=dst->data[idx+kc_offset];
			int pred=clamp(CG(cN-pN, cW-pW, cNW-pNW)+p, half);
#include"f13.h"
			dst->data[idx]=(short)delta;
			//dst->data[idx]=(short)dec2(clamp(CG(cN-pN, cW-pW, cNW-pNW)+p, half), nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void dec_red0 (Image *dst, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels, int g_offset, int b_offset)
{
	short *rows[4]={0};
	for(int ky=0, idx=kc;ky<dst->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, dst->iw, ky);
		for(int kx=0;kx<dst->iw;++kx, idx+=3)
		{
			int pred=(dst->data[idx+g_offset]+dst->data[idx+b_offset])>>1;
#include"f13.h"
			dst->data[idx]=(short)delta;
			//dst->data[idx]=(short)dec2(dst->data[idx+kc_offset], nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void dec_red2 (Image *dst, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels, int g_offset, int b_offset)
{
	short *rows[4]={0};
	int rowstride=dst->iw*3, half=nlevels>>1;
	for(int ky=0, idx=kc;ky<dst->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, dst->iw, ky);
		for(int kx=0;kx<dst->iw;++kx, idx+=3)
		{
			int
				rN	=ky?dst->data[idx-rowstride]:0,
				gN	=ky?dst->data[idx-rowstride+g_offset]:0,
				bN	=ky?dst->data[idx-rowstride+b_offset]:0,
				rW	=kx?dst->data[idx-3]:0,
				gW	=kx?dst->data[idx-3+g_offset]:0,
				bW	=kx?dst->data[idx-3+b_offset]:0,
				g	=dst->data[idx+g_offset],
				b	=dst->data[idx+b_offset];
			int pred=clamp(((2*(rN+rW)-(gN+gW+bN+bW))>>2)+((g+b)>>1), half);
#include"f13.h"
			dst->data[idx]=(short)delta;
			//dst->data[idx]=(short)dec2(clamp(((cN-pN+cW-pW)>>1)+p, half), nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
static void dec_red3 (Image *dst, int kc, short *errors, ArithmeticCoder *ec, int *hist, unsigned *CDF, int nlevels, int nclevels, int nqlevels, int g_offset, int b_offset)
{
	short *rows[4]={0};
	int rowstride=dst->iw*3, half=nlevels>>1;
	for(int ky=0, idx=kc;ky<dst->ih;++ky)
	{
		int alpha=0x4000;
		init_rows(rows, errors, dst->iw, ky);
		for(int kx=0;kx<dst->iw;++kx, idx+=3)
		{
			int
				rNW	=ky&kx?dst->data[idx-rowstride-3]:0,
				gNW	=ky&kx?dst->data[idx-rowstride-3+g_offset]:0,
				bNW	=ky&kx?dst->data[idx-rowstride-3+b_offset]:0,
				rN	=ky?dst->data[idx-rowstride]:0,
				gN	=ky?dst->data[idx-rowstride+g_offset]:0,
				bN	=ky?dst->data[idx-rowstride+b_offset]:0,
				rW	=kx?dst->data[idx-3]:0,
				gW	=kx?dst->data[idx-3+g_offset]:0,
				bW	=kx?dst->data[idx-3+b_offset]:0,
				g	=dst->data[idx+g_offset],
				b	=dst->data[idx+b_offset];
			int pred=clamp((CG(2*rN-(gN+bN), 2*rW-(gW+bW), 2*rNW-(gNW+bNW))+g+b)>>1, half);
#include"f13.h"
			dst->data[idx]=(short)delta;
			//dst->data[idx]=(short)dec2(clamp(CG(cN-pN, cW-pW, cNW-pNW)+p, half), nlevels, rows, ec, hist, CDF, nqlevels);

			++rows[0];
			++rows[1];
			++rows[2];
			++rows[3];
		}
	}
}
#undef  F13_DEC2
int f13_encode(Image const *src, ArrayHandle *data, int loud)
{
	PROF_START();
	if(src->nch!=3)
		LOG_ERROR("Expected 3 channels, got %d", src->nch);
	int depth=src->depth, nlevels=1<<depth, half=nlevels>>1;
	int rowstride=src->iw*3;
	size_t histsize=sizeof(int[NCH])<<src->depth;
	int *hist=(int*)malloc(histsize);
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(hist, 0, histsize);
	double t0=time_sec();
	for(int ky=1;ky<src->ih;++ky)
	{
		for(int kx=1;kx<src->iw;++kx)
		{
#define LOAD(C, X, Y) src->data[idx-rowstride*Y-3*X+C]
			int
				idx=(src->iw*ky+kx)*3,
				r	=LOAD(0, 0, 0),
				rN	=LOAD(0, 0, 1),
				rW	=LOAD(0, 1, 0),
				rNW	=LOAD(0, 1, 1),
				g	=LOAD(1, 0, 0),
				gN	=LOAD(1, 0, 1),
				gW	=LOAD(1, 1, 0),
				gNW	=LOAD(1, 1, 1),
				b	=LOAD(2, 0, 0),
				bN	=LOAD(2, 0, 1),
				bW	=LOAD(2, 1, 0),
				bNW	=LOAD(2, 1, 1);
			int preds[]=
			{
#define CHANNEL(X) X,
				CHLIST
#undef  CHANNEL
			};
			for(int kt=0;kt<NCH;++kt)
			{
				preds[kt]+=half;
				preds[kt]&=nlevels-1;
				++hist[kt<<depth|preds[kt]];
			}
		}
	}
	double t1=time_sec();
	ptrdiff_t field=(src->iw-1LL)*(src->ih-1LL);
	double esizes[NRCH+NGCH+NBCH]={0};
	for(int kt=0;kt<NRCH+NGCH+NBCH;++kt)
	{
		int *curr_hist=hist+((size_t)kt<<depth);
		for(int ks=0;ks<nlevels;++ks)
		{
			int freq=curr_hist[ks];
			if(freq)
				esizes[kt]-=freq*log2((double)freq/field);
		}
		esizes[kt]/=8;
	}
	double t2=time_sec();
	int ridx=0, gidx=0, bidx=0;
	for(int kt=1;kt<NGCH;++kt)
	{
		if(esizes[gidx]>esizes[kt])
			gidx=kt;
	}
	bidx=NGCH;
	for(int kt=NGCH+1;kt<NGCH+NBCH;++kt)
	{
		if(esizes[bidx]>esizes[kt])
			bidx=kt;
	}
	ridx=NGCH+NBCH;
	for(int kt=NGCH+NBCH+1;kt<NCH;++kt)
	{
		if(esizes[ridx]>esizes[kt])
			ridx=kt;
	}
	//gidx=0, bidx=3, ridx=9;//
	if(loud)
	{
		ptrdiff_t usize=field*depth>>3;
		printf("U   %11td\n", usize);
		for(int kt=0;kt<NRCH+NGCH+NBCH;++kt)
		{
			char c=' ';
			if(kt==ridx)
				c='r';
			if(kt==gidx)
				c='g';
			if(kt==bidx)
				c='b';
			printf("%2d  %14.2lf  %6.2lf%%  %14.6lf  %c  %s\n", kt, esizes[kt], 100.*esizes[kt]/usize, usize/esizes[kt], c, chnames[kt]);
		}
		double esize=esizes[ridx]+esizes[gidx]+esizes[bidx];
		printf("Estimated size %.2lf bytes  CR %lf\n", esize, 3*usize/esize);
		printf("Sampling\t%16.6lf sec  %lf B/s\n", t1-t0, usize/(t1-t0));
		printf("Analysis\t%16.6lf sec\n", t2-t1);
	}
	free(hist);
	PROF(PREP);

	double t3=time_sec();
	int token, bypass, nbits;
	token=quantize_pixel(nlevels, &bypass, &nbits);
	int nqlevels=token+1;
	//int nclevels=1;
	int nclevels=quantize_ctx(nlevels)<<1|1;
	histsize=sizeof(int[NCTX])*(nqlevels+1LL)*nclevels;//nclevels ^ NCTX * ntokenlevels
	hist=(int*)malloc(histsize);
	unsigned *CDF=(unsigned*)malloc(histsize);
	//unsigned short *alpha=(unsigned short*)malloc(nclevels*sizeof(short));

	size_t bufsize=sizeof(short[4])*(src->iw+2LL);//4 padded rows * {errors}
	short *errors=(short*)malloc(bufsize);
	if(!hist||!CDF||!errors)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(errors, 0, bufsize);

	ArithmeticCoder ec;
	DList list;
	dlist_init(&list, 1, 256, 0);
	int predidx=gidx;			//Y
	predidx=NBCH*predidx+bidx-NGCH;		//Cb
	predidx=NRCH*predidx+ridx-(NGCH+NBCH);	//Cr
	dlist_push_back1(&list, &predidx);
	ac_enc_init(&ec, &list);
	
	ptrdiff_t csizes[3]={0};//YUV
	init_stats(hist, CDF, nclevels, nqlevels);
	switch(gidx)
	{
	case 0:enc_self0(src, 1, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);break;
	case 1:enc_self2(src, 1, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);break;
	case 2:enc_self3(src, 1, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);break;
	}
	csizes[0]=list.nobj;
	init_stats(hist, CDF, nclevels, nqlevels);
	switch(bidx)
	{
	case 3:enc_self0(src, 2, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);		break;
	case 4:enc_self2(src, 2, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);		break;
	case 5:enc_self3(src, 2, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);		break;
	case 6:enc_aux0 (src, 2, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, -1);	break;
	case 7:enc_aux2 (src, 2, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, -1);	break;
	case 8:enc_aux3 (src, 2, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, -1);	break;
	}
	csizes[1]=list.nobj-csizes[0];
	init_stats(hist, CDF, nclevels, nqlevels);
	switch(ridx)
	{
	case  9:enc_self0(src, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);		break;
	case 10:enc_self2(src, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);		break;
	case 11:enc_self3(src, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);		break;
	case 12:enc_aux0 (src, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, 1);	break;
	case 13:enc_aux2 (src, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, 1);	break;
	case 14:enc_aux3 (src, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, 1);	break;
	case 15:enc_red0 (src, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, 1, 2);	break;
	case 16:enc_red2 (src, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, 1, 2);	break;
	case 17:enc_red3 (src, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, 1, 2);	break;
	}
	ac_enc_flush(&ec);
	dlist_appendtoarray(&list, data);
	csizes[2]=list.nobj-(csizes[0]+csizes[1]);
	if(loud)
	{
		ptrdiff_t usize=((ptrdiff_t)src->iw*src->ih*src->nch*src->depth+7)>>3;
		ptrdiff_t csize=list.nobj;
		t3=time_sec()-t3;
		printf("YUV  %14td  %14td  %14td\n", csizes[0], csizes[1], csizes[2]);
		printf("%14td/%14td = %10.6lf%%  CR %lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
		printf("Encoding\t%16.6lf sec\n", t3);
		printf("E       \t%16.6lf sec\n", time_sec()-t0);
		prof_print();
	}
	dlist_clear(&list);
	free(hist);
	free(CDF);
	free(errors);
	return 0;
}
int f13_decode(const unsigned char *data, size_t len, Image *dst, int loud)
{
	PROF_START();
	double t=time_sec();
	if(dst->nch!=3)
		LOG_ERROR("Expected 3 channels, got %d", dst->nch);
	int depth=dst->depth, nlevels=1<<depth;
	int predidx=*data;
	int ridx, gidx, bidx;
	ridx=predidx%NRCH+NGCH+NBCH,	predidx/=NRCH;//Cr
	bidx=predidx%NBCH+NGCH,		predidx/=NBCH;//Cb
	gidx=predidx%NGCH,		predidx/=NGCH;//Y
	if(predidx)
	{
		LOG_ERROR("Corrupt file");
		return 1;
	}
	++data;
	--len;
	ArithmeticCoder ec;
	ac_dec_init(&ec, data, data+len);

	int token, bypass, nbits;
	token=quantize_pixel(nlevels, &bypass, &nbits);
	int nqlevels=token+1;
	//int nclevels=1;
	int nclevels=quantize_ctx(nlevels)<<1|1;
	size_t histsize=sizeof(int[NCTX])*(nqlevels+1LL)*nclevels;//nclevels ^ NCTX * ntokenlevels
	int *hist=(int*)malloc(histsize);
	unsigned *CDF=(unsigned*)malloc(histsize);
	//unsigned short *alpha=(unsigned short*)malloc(nclevels*sizeof(short));

	size_t bufsize=sizeof(short[4])*(dst->iw+2LL);//4 padded rows * {errors}
	short *errors=(short*)malloc(bufsize);
	if(!hist||!CDF||!errors)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(errors, 0, bufsize);
	
#ifdef DEBUG_CHECKS
	global_idx=0;//
#endif

	init_stats(hist, CDF, nclevels, nqlevels);
	switch(gidx)
	{
	case 0:dec_self0(dst, 1, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);break;
	case 1:dec_self2(dst, 1, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);break;
	case 2:dec_self3(dst, 1, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);break;
	}
	init_stats(hist, CDF, nclevels, nqlevels);
	switch(bidx)
	{
	case 3:dec_self0(dst, 2, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);		break;
	case 4:dec_self2(dst, 2, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);		break;
	case 5:dec_self3(dst, 2, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);		break;
	case 6:dec_aux0 (dst, 2, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, -1);	break;
	case 7:dec_aux2 (dst, 2, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, -1);	break;
	case 8:dec_aux3 (dst, 2, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, -1);	break;
	}
	init_stats(hist, CDF, nclevels, nqlevels);
	switch(ridx)
	{
	case  9:dec_self0(dst, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);		break;
	case 10:dec_self2(dst, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);		break;
	case 11:dec_self3(dst, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels);		break;
	case 12:dec_aux0 (dst, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, 1);	break;
	case 13:dec_aux2 (dst, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, 1);	break;
	case 14:dec_aux3 (dst, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, 1);	break;
	case 15:dec_red0 (dst, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, 1, 2);	break;
	case 16:dec_red2 (dst, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, 1, 2);	break;
	case 17:dec_red3 (dst, 0, errors, &ec, hist, CDF, nlevels, nclevels, nqlevels, 1, 2);	break;
	}
	if(loud)
	{
		t=time_sec()-t;
		printf("D       \t%16.6lf sec\n", t);
		prof_print();
	}
	free(hist);
	free(CDF);
	free(errors);
	return 0;
}