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


#define NPREDS 1
#define NCTX 8
#define PADSIZE 3

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
#define CONFIG_EXP 5
#define CONFIG_MSB 2
#define CONFIG_LSB 0
static void quantize(int val, int *token, int *nbits, int *bypass)
{
	if(val<(1<<CONFIG_EXP))
	{
		*token=val;//token
		if(nbits)
			*nbits=0;
		if(bypass)
			*bypass=0;
	}
	else
	{
		int lgv=floor_log2_32((unsigned)val);
		int mantissa=val-(1<<lgv);
		*token = (1<<CONFIG_EXP) + (
				(lgv-CONFIG_EXP)<<(CONFIG_MSB+CONFIG_LSB)|
				(mantissa>>(lgv-CONFIG_MSB))<<CONFIG_LSB|
				(mantissa&((1<<CONFIG_LSB)-1))
			);
		lgv-=CONFIG_MSB+CONFIG_LSB;
		if(nbits)
			*nbits=lgv;
		if(bypass)
			*bypass=val>>CONFIG_LSB&((1LL<<lgv)-1);
	}
}
#if 0
static int dequantize(int token)
{
	if(token<(1<<CONFIG_EXP))
	{
		//if(ret_nbits)
		//	*ret_nbits=0;
		return token;
	}
	
	token-=1<<CONFIG_EXP;
	int lsb=token&((1<<CONFIG_LSB)-1);
	token>>=CONFIG_LSB;
	int msb=token&((1<<CONFIG_MSB)-1);
	token>>=CONFIG_MSB;
	int nbits=token+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB);
	
	token=1;
	token<<=CONFIG_MSB;
	token|=msb;
	token<<=nbits;
	//token|=bypass;
	token<<=CONFIG_LSB;
	token|=lsb;
	//if(ret_nbits)
	//	*ret_nbits=nbits;
	return token;
}
#endif
static int ctx_quantize(int ctx, int clevels)
{
	int negmask=-(ctx<0);
	ctx=ctx<<1^negmask;
	ctx=floor_log2_32(ctx);
	ctx^=negmask;
	ctx-=negmask;
	ctx+=clevels>>1;
	return ctx;
}
static int ctx_dequantize_floor(int ctx, int clevels)
{
	ctx-=clevels>>1;
	int negmask=-(ctx<0);
	ctx^=negmask;
	ctx-=negmask;
	ctx=1<<ctx;
	ctx^=negmask;
	ctx-=negmask;
	return ctx;
}
static int get_qlevels(int depth)
{
	switch(depth)
	{
	case 0:		return 1;
	case 8:		return 45;
	case 9:		return 49;
	case 16:	return 77;
	case 17:	return 81;
	}
	LOG_ERROR("Unsupported bit depth");
	return 0;
}
#define LOAD_P(X, Y) buf[(kym[-(Y)]+kx+PADSIZE+(X))<<3|0|kc]
#define LOAD_E(X, Y) buf[(kym[-(Y)]+kx+PADSIZE+(X))<<3|4|kc]
static void pt2(ArithmeticCoder *ec, int *curr, int *buf, char *stats, int *kym, int kc, int kx, int ky, int nlevels, int qlevels, int clevels, int fwd)
{
	int
		NNE	=LOAD_P( 1, -2),
		NW	=LOAD_P(-1, -1),
		N	=LOAD_P( 0, -1),
		NE	=LOAD_P( 1, -1),
		NEE	=LOAD_P( 2, -1),
		NEEE	=LOAD_P( 3, -1),
		W	=LOAD_P(-1,  0),

		eNW	=LOAD_E(-1, -1),
		eN	=LOAD_E( 0, -1),
		eNE	=LOAD_E( 1, -1),
		eW	=LOAD_E(-1,  0);

	int pred=N+W-NW;
	pred=MEDIAN3(N, W, pred);

	int ctx[]=
	{
		(N+W)>>1,
		//(eN+eW)>>1,
		//(N+eN)>>1,
		//(W+eW)>>1,
		W+NE-N,
		(W+NEE)>>1,
		//(3*W+NEEE)>>2,
		//N+NE-NNE,
	};
	//int ctx=eN+eW-eNW;
	//ctx=MEDIAN3(eN, eW, eNW);
	int qctx[_countof(ctx)], cstart[_countof(ctx)], cend[_countof(ctx)], alpha[_countof(ctx)];
	for(int k=0;k<_countof(ctx);++k)
	{
		qctx[k]=ctx_quantize(ctx[k], clevels);
		cstart[k]=ctx_dequantize_floor(qctx[k]-1, clevels);
		cend[k]=ctx_dequantize_floor(qctx[k]+1, clevels);
		if((unsigned)(ctx[k]-cstart[k])>=(unsigned)(cend[k]-cstart[k]))
			LOG_ERROR("Alg error");
		alpha[k]=(int)(((long long)(ctx[k]-cstart[k])<<16)/(cend[k]-cstart[k]));
		if(alpha[k]<0||alpha[k]>0x10000)
			LOG_ERROR("Bad alpha");
	}
	//int cstart=ctx_dequantize_floor(qctx-1, clevels), cend=ctx_dequantize_floor(qctx+1, clevels);
	//int negctx=ctx<0;
	//int actx=ctx<<1^-(ctx<0);
	//int qctx=0;
	//quantize(actx, &qctx, 0, 0);
	//int cstart=dequantize(qctx), cend=dequantize(qctx+(negctx?-1:1));
	//if((unsigned)(ctx-cstart)>=(unsigned)(cend-cstart))
	//	LOG_ERROR("Alg error");
	//int alpha=(int)(((long long)(ctx-cstart)<<16)/(cend-cstart));

	//if(!kc&&kx==256&&ky==5)//
	//	printf("");

	int error=0, token=0, bypass=0, nbits=0;
	if(fwd)
	{
		error=*curr-pred;
		error+=nlevels>>1;
		error&=nlevels-1;
		error-=nlevels>>1;
		quantize(error<<1^-(error<0), &token, &nbits, &bypass);
	}
	for(int kb=6, abac_idx=1;kb>=0;--kb)
	{
		char *curr_stats=stats+37*NCTX*(4*abac_idx+kc);
		int p0=0;
		for(int kp=0;kp<_countof(ctx);++kp)
		{
			p0+=curr_stats[qctx[kp]]<<8;
			p0+=((curr_stats[qctx[kp]-1]<<16)+(curr_stats[qctx[kp]+1]-curr_stats[qctx[kp]])*alpha[kp])>>8;
		}
		p0/=_countof(ctx)<<1;
		p0+=0x8000;
		//p0>>=1;
		p0=CLAMP(1, p0, 0xFFFF);

		int bit;
		if(fwd)
		{
			bit=token>>kb&1;
			ac_enc_bin(ec, p0, bit);
		}
		else
		{
			bit=ac_dec_bin(ec, p0);
			token|=bit<<kb;
		}

		int p0_perf=(!bit<<8)-128;
		for(int kp=0;kp<_countof(ctx);++kp)
		{
			int update=curr_stats[qctx[kp]-1]+((p0_perf-curr_stats[qctx[kp]-1])*(0x10000-alpha[kp])>>20);
			curr_stats[qctx[kp]-1]=CLAMP(-128, update, 127);

			update=curr_stats[qctx[kp]]+((p0_perf-curr_stats[qctx[kp]])>>3);
			curr_stats[qctx[kp]]=CLAMP(-128, update, 127);

			update=curr_stats[qctx[kp]+1]+((p0_perf-curr_stats[qctx[kp]+1])*alpha[kp]>>20);
			curr_stats[qctx[kp]+1]=CLAMP(-128, update, 127);
		}

		abac_idx<<=1;
		abac_idx|=bit;
	}
	if(fwd)
	{
		if(nbits)
		{
			while(nbits>8)
			{
				ac_enc(ec, bypass>>(nbits-8)&0xFF, 0, 1<<8, 16-8);
				nbits-=8;
			}
			ac_enc(ec, bypass&((1<<nbits)-1), 0, 1<<nbits, 16-nbits);
		}
	}
	else
	{
		error=token;
		if(error>=(1<<CONFIG_EXP))
		{
			error-=1<<CONFIG_EXP;
			int lsb=error&((1<<CONFIG_LSB)-1);
			error>>=CONFIG_LSB;
			int msb=error&((1<<CONFIG_MSB)-1);
			error>>=CONFIG_MSB;
			int nbits=error+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB), n=nbits;
			int bypass=0;
			while(n>8)
			{
				n-=8;
				bypass|=ac_dec(ec, 0, 1<<8, 16-8)<<n;
			}
			bypass|=ac_dec(ec, 0, 1<<n, 16-n);
			error=1;
			error<<=CONFIG_MSB;
			error|=msb;
			error<<=nbits;
			error|=bypass;
			error<<=CONFIG_LSB;
			error|=lsb;
		}
		error=error>>1^-(error&1);
		*curr=error+pred;
		*curr+=nlevels>>1;
		*curr&=nlevels-1;
		*curr-=nlevels>>1;
	}
	LOAD_P(0, 0)=*curr;
	LOAD_E(0, 0)=error;
}
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
int t55_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t_start=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif
	int nch=(image->depth[0]!=0)+(image->depth[1]!=0)+(image->depth[2]!=0)+(image->depth[3]!=0);
	UPDATE_MIN(nch, image->nch);
	int depths[]=
	{
		image->depth[1],
		image->depth[2]+1,
		image->depth[0]+1,
		image->depth[3],
	};
	int qlevels[]=
	{
		get_qlevels(depths[0]),
		get_qlevels(depths[1]),
		get_qlevels(depths[2]),
		get_qlevels(depths[3]),
	};
	char maxdepth=depths[0];
	UPDATE_MAX(maxdepth, depths[1]);
	UPDATE_MAX(maxdepth, depths[2]);
	UPDATE_MAX(maxdepth, depths[3]);
	int maxlevels=1<<maxdepth;
	int *buf=(int*)malloc((image->iw+PADSIZE*2LL)*sizeof(int[4*PADSIZE*2*2]));//4 channels * (PARSIZE*2) rows * (pixels + errors)
	char *stats=(char*)malloc(sizeof(char[4*128*37*NCTX]));//4 channels * 128 max tree size * 82 max qlevels in context * NCTX
	if(!buf||!stats)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(buf, 0, (image->iw+PADSIZE*2LL)*sizeof(int[4*PADSIZE*2*2]));
	memset(stats, 0, sizeof(char[4*128*37*NCTX]));
	//char fillval=0x80;
	//memfill(stats, &fillval, sizeof(char[4*128*37*NCTX]), sizeof(fillval));
	DList list;
	ArithmeticCoder ec;
	dlist_init(&list, 1, 256, 0);
	if(fwd)
		ac_enc_init(&ec, &list);
	else
		ac_dec_init(&ec, cbuf, cbuf+clen);
	for(int ky=0, idx=0;ky<image->ih;++ky)
	{
		int kym[PADSIZE*2]={0};
		for(int k=0;k<_countof(kym);++k)
		{
			int y=ky-k;
			MODVAR(y, y, PADSIZE*2);
			kym[k]=(image->iw+PADSIZE*2)*y;
		}
		for(int kx=0;kx<image->iw;++kx, ++idx)
		{
			int kc=0;
			int comp[4]={0};
			if(fwd)
				memcpy(comp, image->data+((size_t)idx<<2), sizeof(int[4]));
			if(nch>=3)
			{
				if(fwd)
				{
					comp[0]-=comp[1];
					comp[2]-=comp[1];
					comp[1]+=(comp[0]+comp[2])>>2;
				}
				pt2(&ec, comp+1, buf, stats, kym, 0, kx, ky, 1<<depths[0], qlevels[0], (depths[0]<<1)+3, fwd);
				pt2(&ec, comp+2, buf, stats, kym, 1, kx, ky, 1<<depths[1], qlevels[1], (depths[1]<<1)+3, fwd);
				pt2(&ec, comp+0, buf, stats, kym, 2, kx, ky, 1<<depths[2], qlevels[2], (depths[2]<<1)+3, fwd);
				if(!fwd)
				{
					comp[1]-=(comp[0]+comp[2])>>2;
					comp[2]+=comp[1];
					comp[0]+=comp[1];
				}
				kc+=3;
			}
			for(;kc<nch;++kc)
				pt2(&ec, comp+kc, buf, stats, kym, kc, kx, ky, 1<<depths[3], qlevels[3], (depths[3]<<1)+3, fwd);
#ifdef ENABLE_GUIDE
			if(!fwd&&memcmp(comp, guide->data+((size_t)idx<<2), sizeof(int[4])))
			{
				int comp2[4]=
				{
					guide->data[idx<<2|0],
					guide->data[idx<<2|1],
					guide->data[idx<<2|2],
					guide->data[idx<<2|3],
				};
				comp[0]-=comp[1];
				comp[2]-=comp[1];
				comp[1]+=(comp[0]+comp[2])>>2;

				comp2[0]-=comp2[1];
				comp2[2]-=comp2[1];
				comp2[1]+=(comp2[0]+comp2[2])>>2;
				printf("YUV0 %d %d %d\n", comp[0], comp[1], comp[2]);
				printf("YUV1 %d %d %d\n", comp2[0], comp2[1], comp2[2]);
				LOG_ERROR("Guide error XY %d %d", kx, ky);
			}
#endif
			if(!fwd)
				memcpy(dst->data+((size_t)idx<<2), comp, sizeof(int[4]));
		}
	}
	if(fwd)
	{
		ac_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
	}
	if(loud)
	{
		t_start=time_sec()-t_start;
		if(fwd)
		{
			double usize=image_getBMPsize(image);
			printf("csize %12lld  %7.3lf%%\n", list.nobj, 100.*list.nobj/usize);
		}
		printf("%c %12.3lf sec\n", 'D'+fwd, t_start);
	}
	dlist_clear(&list);
	free(buf);
	free(stats);
	return 1;
}