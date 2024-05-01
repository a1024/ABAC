#include"fast.h"
#include<stdlib.h>
#include<math.h>
#define EC_USE_ARRAY
#define AC_IMPLEMENTATION
#include"ac.h"
static const char file[]=__FILE__;
#define CODECTAG 'a'


//	#define ENABLE_GUIDE
//	#define ENABLE_SSE
//	#define ESTIMATE_CHSIZES

//activate only one:
	#define USE_GRCODER
//	#define USE_ABAC

	#define GRCODER_USE_NPOT
//	#define DUALCTRS//X


#define PADSIZE 2	//power of two

#ifdef USE_GRCODER
#define EC_TYPE GolombRiceCoder
#define EC_ENC_INIT	gr_enc_init
#define EC_DEC_INIT	gr_dec_init
#define EC_FLUSH	gr_enc_flush
#elif defined USE_ABAC
#define EC_TYPE ArithmeticCoder
#define EC_ENC_INIT	ac_enc_init
#define EC_DEC_INIT	ac_dec_init
#define EC_FLUSH	ac_enc_flush
#else
#define EC_TYPE ArithmeticCoder
#define EC_ENC_INIT	ac_enc_init
#define EC_DEC_INIT	ac_dec_init
#define EC_FLUSH	ac_enc_flush
#endif
#ifndef USE_GRCODER
static int quantize(int val, int clevels)
{
	int negmask=-(val<0);
	val=abs(val);
	val=floor_log2_32(val)+1;
	val>>=1;
	val^=negmask;
	val-=negmask;
	val+=clevels>>1;
	return val;
}
#endif
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
int f01_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t_start=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif
	int depth=image->depth;
	if(image->nch>1)
	{
		++depth;
		if(depth>16)
			depth=16;
	}
#if 0
	{
		int *hist=(int*)malloc(sizeof(int[]));
		if(!hist)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		free(hist);
	}
#endif
	int nlevels=1<<depth, half=nlevels>>1;
#ifdef USE_ABAC
	int clevels=quantize(half, 0)<<1|1;
	int nctx=clevels*clevels;
	unsigned short *stats=(unsigned short*)malloc(sizeof(short)*image->nch*nctx<<depth);
#endif
	size_t rowssize=sizeof(short[2*2*PADSIZE])*image->nch*(image->iw+PADSIZE*2LL);
	short *rows=(short*)malloc(rowssize);//N padded rows * {pixels, errors} * nch
#ifdef ENABLE_SSE
	long long *sse=(long long*)malloc(sizeof(long long)*image->nch*nctx);
#endif
	if(!rows
#ifdef USE_ABAC
		||!stats
#endif
#ifdef ENABLE_SSE
		||!sse
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
#ifdef USE_ABAC
#ifdef DUALCTRS
	unsigned short fillval=0x0101;
#else
	unsigned short fillval=0x8000;
#endif
	memfill(stats, &fillval, sizeof(short)*image->nch*nctx<<depth, sizeof(short));
#endif
	for(int k=0;k<(image->iw+PADSIZE*2)*2*PADSIZE;++k)//initialize errors far from zero, for maximum bypass when using Golomb-Rice code
	{
		for(int kc=0;kc<image->nch;++kc)
		{
			//if((size_t)(image->nch*(k<<1|0)+kc)>=rowssize)
			//	LOG_ERROR("");
			//if((size_t)(image->nch*(k<<1|1)+kc)>=rowssize)
			//	LOG_ERROR("");
			rows[image->nch*(k<<1|0)+kc]=0;
			rows[image->nch*(k<<1|1)+kc]=(short)(half>>1);
		}
	}
	//memset(rows, 0, sizeof(short[2*2])*image->nch*(image->iw+2LL));
#ifdef ENABLE_SSE
	memset(sse, 0, sizeof(long long)*image->nch*nctx);
#endif
#ifdef ESTIMATE_CHSIZES
	size_t chsizes[4]={0};
	int compressidx=0;
#endif
	//DList list;
	EC_TYPE ec;
	//dlist_init(&list, 1, 256, 0);
	ptrdiff_t nvals=(ptrdiff_t)image->iw*image->ih*image->nch;
	size_t maxsize=((size_t)nvals*image->depth+7)/8;
	size_t startidx=0;
	if(fwd)
	{
#if 0
		if(loud)
		{
			volatile double t1=time_sec();
			for(int k=0;k<nvals;++k)//iteration speed test
			{
				int val=image->data[k];
				int negmask=-(val<0);
				val=abs(val);
				val=floor_log2_32(val);
				val^=negmask;
				val-=negmask;
				image->data[k]=val;
			}
			//	image->data[k]&=half-1;
			t1=time_sec()-t1;
			printf("iterating through the image takes %lf sec\n", t1);
		}
#endif
		startidx=array_append(data, 0, 1, maxsize+1, 1, 0, 0);
		data[0]->data[startidx]=CODECTAG;
		EC_ENC_INIT(&ec, data[0]->data+startidx+1, data[0]->data+startidx+maxsize);
	}
		//EC_ENC_INIT(&ec, &list);
	else
	{
		if(!clen)
		{
			memset(image->data, 0, sizeof(short)*image->iw*image->ih*image->nch);
			goto finish_dec_bypass;
		}
		if(!cbuf[0])
		{
			if(image->depth<=8)
			{
				for(ptrdiff_t k=0;k<nvals;++k)
					dst->data[k]=(char)cbuf[k+1];
			}
			else
				memcpy(dst->data, cbuf+1, nvals*sizeof(short));
			goto finish_dec_bypass;
		}
		if(cbuf[0]!=CODECTAG)
		{
			LOG_ERROR("Invalid tag");
			return 0;
		}
		++cbuf;//skip bypass tag
		--clen;
		EC_DEC_INIT(&ec, cbuf, cbuf+clen);
	}
	for(int ky=0, idx=0;ky<image->ih;++ky)
	{
		int kym[]=
		{
			(image->iw+PADSIZE*2)*((ky-0)&(2*PADSIZE-1)),
			(image->iw+PADSIZE*2)*((ky-1)&(2*PADSIZE-1)),
			(image->iw+PADSIZE*2)*((ky-2)&(2*PADSIZE-1)),
			(image->iw+PADSIZE*2)*((ky-3)&(2*PADSIZE-1)),
		};
		for(int kx=0;kx<image->iw;++kx, ++idx)
		{
			//if(kx==256&&ky==256)
			//	printf("");

//#define LOADIDX(C, X, Y, E) (image->nch*((kym[-(Y)]+PADSIZE+kx+(X))<<1|E)+C>=rowssize?LOG_ERROR(""):image->nch*((kym[-(Y)]+PADSIZE+kx+(X))<<1|E)+C)
#define LOADIDX(C, X, Y, E) (image->nch*((kym[-(Y)]+PADSIZE+kx+(X))<<1|E)+C)
			short *pixels=image->data+(image->nch*idx+0);
			short *comp=rows+LOADIDX(0, 0, 0, 0);
			if(fwd)
			{
				memcpy(comp, pixels, image->nch*sizeof(short));
				for(int k=1;k<image->nch;++k)
					comp[k]-=comp[0];
			}
			for(int kc=0;kc<image->nch;++kc)
			{
				int
					NW	=rows[LOADIDX(kc, -1, -1, 0)],
					N	=rows[LOADIDX(kc,  0, -1, 0)],
					W	=rows[LOADIDX(kc, -1,  0, 0)],
				//	eNNWW	=rows[LOADIDX(kc, -2, -2, 1)],
				//	eNNW	=rows[LOADIDX(kc, -1, -2, 1)],
				//	eNN	=rows[LOADIDX(kc,  0, -2, 1)],
				//	eNNE	=rows[LOADIDX(kc,  1, -2, 1)],
				//	eNNEE	=rows[LOADIDX(kc,  2, -2, 1)],
				//	eNWW	=rows[LOADIDX(kc, -2, -1, 1)],
					eNW	=rows[LOADIDX(kc, -1, -1, 1)],
					eN	=rows[LOADIDX(kc,  0, -1, 1)],
					eNE	=rows[LOADIDX(kc,  1, -1, 1)],
				//	eNEE	=rows[LOADIDX(kc,  2, -1, 1)],
				//	eWW	=rows[LOADIDX(kc, -2,  0, 1)],
					eW	=rows[LOADIDX(kc, -1,  0, 1)];
				int pred=N+W-NW;
				pred=MEDIAN3(N, W, pred);
#ifdef USE_GRCODER
				int ctx;
				//ctx=abs(N-W)+((eN+eW)>>8);//best so far
				ctx=abs(N-W)+((eN+eW+eNW+eNE)>>9);
				//ctx=abs(N-W)+((eN+eW+eNW+eNE+eNN+eWW+eNEE)>>9);
				//ctx=abs(N-W);
				//if(kc)
				//	ctx=(eN+eW+eNW+eNE+eNN+eWW+eNEE)>>8;
				//else
				//	ctx=abs(N-W);
				//int ctx=(eNNWW+eNNW+eNN+eNNE+eNNEE+eNWW+eNW+eN+eNE+eNEE+eWW+eW)>>(depth+1);
				//if(ctx<0)
				//	LOG_ERROR("");
				
				//int ctx=(abs(eN)+abs(eW)+abs(eNW)+abs(eNE))>>10;
				//int ctx=(abs(eN)+abs(eW)+abs(eNW)+abs(eNE)+abs(eNN)+abs(eWW))>>10;
				//int ctx=(abs(eN)+abs(eW))>>9;
				//int ctx=(abs(eN)+abs(eW)+abs(eNW)+abs(eNE))>>3;
				//int ctx=(3*(abs(eN)+abs(eW))+abs(eNW)+abs(eNE))>>3;
				//int ctx=((!ky?clevels-1:abs(eN))+(!kx?clevels-1:abs(eW)))>>1;
				//int ctx=floor_log2_32((!ky?clevels-1:abs(eN))+(!kx?clevels-1:abs(eW)))+1;
				//int ctx=((!ky?clevels-1:floor_log2_32(abs(eN))+1)+(!kx?clevels-1:floor_log2_32(abs(eN))+1))>>1;
#else
				int ctx=clevels*(!ky?clevels-1:quantize(eN, clevels))+(!kx?clevels-1:quantize(eW, clevels));
#endif
#ifdef ENABLE_SSE
				long long *curr_sse=sse+(size_t)image->nch*ctx+kc;
				long long sum=*curr_sse>>12;
				int count=(int)(*curr_sse&0xFFF);
				pred+=(int)(sum/(count+1LL));
				pred=CLAMP(-half, pred, half-1);
#endif
				
#ifdef USE_GRCODER
#ifndef GRCODER_USE_NPOT
				ctx=floor_log2_32(ctx)+1;
#endif
				int delta, sym, curr;
				//if(ky==1&&kx==102&&kc==0)//
				//if(!idx)//
				//if(kx==5&&ky==1&&kc==2)//
				//	printf("");
				if(fwd)
				{
					delta=comp[kc]-pred;
					delta+=half;
					delta&=nlevels-1;
					delta-=half;
					sym=delta<<1^-(delta<0);
#ifdef GRCODER_USE_NPOT
					if(!gr_enc(&ec, sym, ctx))
					{
						if(loud)
						{
							printf(
								"Compression FAILED\n"
								"    CXY/CWH  %d %d %d / %d %d %d\n"
								"    IDX %d/%d = %lf%%\n"
								"Bypassing...\n",
								kc, kx, ky, image->nch, image->iw, image->ih,
								idx, (int)nvals, 100.*(image->nch*idx+kc)/nvals
							);
#ifdef ESTIMATE_CHSIZES
							compressidx=image->nch*idx+kc;
#endif
						}
						goto finish_enc_bypass;
					}
#ifdef ESTIMATE_CHSIZES
					static ptrdiff_t bitidx=0;
					if(!idx&&!kc)
						bitidx=0;
					ptrdiff_t bitsize=((ec.dstptr-ec.dststart)<<3)+64-ec.nbits;
					//if(bitsize<bitidx)
					//	LOG_ERROR("");
					chsizes[kc]+=bitsize-bitidx;
					bitidx=bitsize;
#endif
#else
					gr_enc_bin(&ec, sym, ctx);
#endif
				}
				else
				{
#ifdef GRCODER_USE_NPOT
					sym=gr_dec(&ec, ctx);
#else
					sym=gr_dec_bin(&ec, ctx);
#endif
					delta=sym>>1^-(sym&1);
					curr=delta+pred;
					curr+=half;
					curr&=nlevels-1;
					curr-=half;
					comp[kc]=(short)curr;
				}
#else
				int delta=0;
				if(fwd)
				{
					delta=comp[kc]-pred;
					delta+=half;
					delta&=nlevels-1;
				}
				unsigned short *curr_stats=stats+(((size_t)image->nch*ctx+kc)<<depth);
				int idx=1;
				for(int kb=depth-1;kb>=0;--kb)
				{
#ifdef DUALCTRS
					unsigned char *p=(unsigned char*)(curr_stats+idx);
					int p0=(p[0]<<16)/(p[0]+p[1]);
					p0=CLAMP(1, p0, 0xFFFF);
#else
					unsigned short *p=curr_stats+idx;
					int p0=*p;
#endif
					int bit;
					if(fwd)
					{
						bit=delta>>kb&1;
						ac_enc_bin(&ec, p0, bit);
					}
					else
					{
						bit=ac_dec_bin(&ec, p0);
						delta|=bit<<kb;
					}
#ifdef DUALCTRS
					if(p[bit]+1>255)
					{
						p[0]>>=1;
						p[1]>>=1;
					}
					++p[bit];
#else
					p0+=((!bit<<16)-p0)>>5;
					*p=CLAMP(1, p0, 0xFFFF);
#endif
					idx=idx<<1|bit;
				}
				if(!fwd)
				{
					int curr=delta;
					curr+=pred;
					curr&=nlevels-1;
					curr-=half;
					comp[kc]=curr;
				}
				delta-=half;
#endif
				rows[LOADIDX(0, 0, 0, 1)]=(short)abs(delta);
#ifdef ENABLE_SSE
				sum+=delta;
				++count;
				if(count>640)
				{
					count>>=1;
					sum>>=1;
				}
				*curr_sse=sum<<12|count;
#endif
			}
			if(!fwd)
			{
				memcpy(pixels, comp, image->nch*sizeof(short));
				for(int k=1;k<image->nch;++k)
					pixels[k]+=pixels[0];
#ifdef ENABLE_GUIDE
				if(guide&&memcmp(pixels, guide->data+(image->nch*idx+0), image->nch*sizeof(short)))
				{
					short pixels0[4]={0};
					memcpy(pixels0, guide->data+(image->nch*idx+0), image->nch*sizeof(short));
					for(int k=1;k<image->nch;++k)
					{
						pixels[k]-=pixels[0];
						pixels0[k]-=pixels0[0];
					}
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
			}
		}
	}
	if(fwd)
	{
		if(EC_FLUSH(&ec))
			data[0]->count=ec.dstptr+1-ec.dststart;
		else
		{
			if(loud)
				printf(
					"Compression FAILED\n"
					"    at final flush  (100%%)\n"
					"Bypassing...\n"
				);
		finish_enc_bypass:
			ec.dststart[-1]=0;
			if(image->depth<=8)
			{
				for(ptrdiff_t k=0;k<nvals;++k)
					ec.dststart[k]=(char)image->data[k];
			}
			else
				memcpy(ec.dststart, image->data, nvals*sizeof(short));
			//ec.dstptr=ec.dstend;
		}
		//dlist_appendtoarray(&list, data);
	}
finish_dec_bypass:
	if(loud)
	{
		t_start=time_sec()-t_start;
		if(fwd)
		{
#ifdef ESTIMATE_CHSIZES
			ptrdiff_t uchsize=compressidx?(ptrdiff_t)compressidx*image->depth/(8LL*image->nch):((ptrdiff_t)image->iw*image->ih*image->depth+7)/8;
			size_t ctotal=0;
			if(compressidx)
				printf("Gave up at  %.3lf / %lld bytes\n", compressidx*image->depth/8., maxsize);
			for(int kc=0;kc<image->nch;++kc)
			{
				double csize=chsizes[kc]/8.;
				printf("C%d  %12.3lf  %15.6lf%%  %8.6lf\n", kc, csize, 100.*csize/uchsize, uchsize/csize);
				ctotal+=chsizes[kc];
			}
			if(compressidx)
				printf("T   %12.3lf  %15.6lf%%  %8.6lf\n",
					ctotal/8., 100.*ctotal/(compressidx*image->depth), compressidx/(ctotal/8.)
				);
#endif
			ptrdiff_t usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)/8;
			ptrdiff_t csize=data[0]->count-startidx;
			printf("csize %12td  %10.6lf%%  %8.6lf\n", csize, 100.*csize/usize, (double)usize/csize);
		}
		printf("%c %15.6lf sec\n", 'D'+fwd, t_start);
		//prof_print();
	}
	//int step=!fwd-fwd;
	//for(size_t ky=fwd?(size_t)image->ih-1:0;ky<(size_t)image->ih;ky+=step)//ANS loop
	//{
	//	for(size_t kx=fwd?(size_t)image->iw-1:0;kx<(size_t)image->iw;kx+=step)
	//	{
	//		for(int kc=fwd?(size_t)image->nch-1;kc<image->nch;kc+=step)
	//		{
	//		}
	//	}
	//}
	//dlist_clear(&list);
#ifdef USE_ABAC
	free(stats);
#endif
	free(rows);
	return 1;
}