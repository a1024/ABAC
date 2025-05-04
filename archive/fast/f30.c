#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT


#include"ac.h"
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif

#define BLOCKDX 768
#define BLOCKDY 768
#define MAXPRINTEDBLOCKS 500


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

typedef struct _ThreadArgs
{
	const Image *src;
	Image *dst;
	
	int fwd, loud, x1, x2, y1, y2;
	int pixels[(BLOCKDX+16)*4*4*2];//4 padded rows  *  EXTENDED_OCH_COUNT

	DList list;
	const unsigned char *decstart, *decend;

	int *hist;
	int histsize, clevels, tlevels;
	
	int blockidx;
} ThreadArgs;
static void block_thread(void *param)
{
	AC3 ec;
	ThreadArgs *args=(ThreadArgs*)param;
	const Image *image=args->fwd?args->src:args->dst;
	//int nlevels=1<<image->depth, half=nlevels>>1, mask=nlevels-1;
	int nlevels[]=
	{
		2<<image->depth,
		1<<image->depth,
		2<<image->depth,
		1<<image->depth,
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};

	if(args->fwd)
	{
		dlist_init(&args->list, 1, 1024, 0);
		ac3_enc_init(&ec, &args->list);
	}
	else
	{
		ac3_dec_init(&ec, args->decstart, args->decend);
	}
	int res=(args->x2-args->x1)*(args->y2-args->y1);
	memset(args->hist, 0, args->histsize);
	memset(args->pixels, 0, sizeof(args->pixels));
	for(int ky=args->y1;ky<args->y2;++ky)//enc loop
	{
		int kx=args->x1;
		ALIGN(16) int *rows[]=
		{
			args->pixels+((BLOCKDX+16LL)*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKDX+16LL)*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKDX+16LL)*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+((BLOCKDX+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		short rgba[4]={0};
		for(;kx<args->x2;++kx)
		{
			int idx=image->nch*(image->iw*ky+kx);
			if(args->fwd)
			{
				memcpy(rgba, image->data+idx, sizeof(short)*image->nch);
				if(image->nch>=3)
				{
					rgba[0]-=rgba[1];
					rgba[2]-=rgba[1];
					rgba[1]+=(rgba[0]+rgba[2])>>2;
				}
			}
			for(int kc=0;kc<image->nch;++kc)
			{
				int
					NW	=rows[1][kc-1*4*2+0],
					N	=rows[1][kc+0*4*2+0],
					W	=rows[0][kc-1*4*2+0],
					eNW	=rows[1][kc-1*4*2+4],
					eN	=rows[1][kc+0*4*2+4],
					eNE	=rows[1][kc+1*4*2+4],
				//	eNEEE	=rows[1][kc+3*4*2+4],
					eW	=rows[0][kc+0*4*2+4],
					*curr	=rows[0]+kc;
				int pred, ctx;
				MEDIAN3_32(pred, N, W, N+W-NW);
				ctx=(eN+eW)*3+eNW+eNE+1;
				ctx=FLOOR_LOG2(ctx);
				int *curr_hist=args->hist+(args->tlevels+1)*(args->clevels*kc+ctx);

				int error, sym, token, bypass, nbits;

				//if(ky==1&&kx==5)
				//	printf("");

				int den=curr_hist[args->tlevels]+args->tlevels, cdf=0, freq;
				if(args->fwd)
				{
					error=rgba[kc]-pred;
					{
						int upred=halfs[kc]-abs(pred), aval=abs(error);
						if(aval<=upred)
						{
							sym=error;
#ifdef ENABLE_BIASCORR
							{
								int negmask=-((ibias_corr<0)&(sym!=-halfs[kc]));//sign is flipped if SSE correction was negative, to skew the histogram
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
					for(int k=0;k<token;++k)
						cdf+=curr_hist[k]+1;
					freq=curr_hist[token]+1;

					ac3_enc_update_NPOT(&ec, cdf, freq, den);
					if(nbits)
						ac3_enc_bypass(&ec, bypass, nbits);
				}
				else
				{
					int code=ac3_dec_getcdf_NPOT(&ec, den);
					token=0;
					cdf=0;
					for(;;)
					{
						int cdf2;

						freq=curr_hist[token]+1;
						cdf2=cdf+freq;
						if(cdf2>code)
							break;
						cdf=cdf2;
						++token;
					}
					ac3_dec_update_NPOT(&ec, cdf, freq, den);
					
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
						int upred=halfs[kc]-abs(pred), negmask=0;
						if(sym<=(upred<<1))
						{
							error=sym>>1^-(sym&1);
#ifdef ENABLE_BIASCORR
							negmask=-((ibias_corr<0)&(error!=-halfs[ch]));
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
					rgba[kc]=error+pred;
				}
				curr[0]=rgba[kc];
				curr[4]=abs(error);
				//curr[4]=(2*eW+abs(error)+eNEEE+2)>>2;
				//if(curr[4]<0)
				//	LOG_ERROR("");
				
				if(curr_hist[args->tlevels]>=2048)
				{
					int sum=0;
					for(int k=0;k<args->tlevels;++k)
						sum+=curr_hist[k]>>=1;
					curr_hist[args->tlevels]=sum;
				}
				++curr_hist[token];
				++curr_hist[args->tlevels];
			}
			if(!args->fwd)
			{
				if(image->nch>=3)
				{
					rgba[1]-=(rgba[0]+rgba[2])>>2;
					rgba[2]+=rgba[1];
					rgba[0]+=rgba[1];
				}
				memcpy(args->dst->data+idx, rgba, sizeof(short)*image->nch);
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
	}
	if(args->fwd)
		ac3_enc_flush(&ec);
}
int f30_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
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
	int histsize, clevels, tlevels;
	
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	{
		int token, bypass, nbits;
		quantize_pixel(2<<image->depth, &token, &bypass, &nbits);
		tlevels=token+1;
		clevels=FLOOR_LOG2(8<<image->depth);
		histsize=(size_t)image->nch*clevels*(tlevels+1LL)*sizeof(int);
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
		
		arg->hist=(int*)malloc(histsize);
		if(!arg->hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		arg->histsize=histsize;
		arg->clevels=clevels;
		arg->tlevels=tlevels;

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
						printf(
							"block %4d/%4d  XY %3d %3d  %4d*%4d:  %8d->%8zd bytes  %10.6lf%%  CR %10lf\n",
							kt+kt2+1, nblocks,
							kx, ky,
							arg->y2-arg->y1,
							arg->x2-arg->x1,
							blocksize,
							arg->list.nobj,
							100.*arg->list.nobj/blocksize,
							(double)blocksize/arg->list.nobj
						);
					}
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
			printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
			printf("Mem usage: ");
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
		}
		printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		free(arg->hist);
	}
	free(args);
	return 0;
}