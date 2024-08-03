#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#ifdef _MSC_VER
#include<intrin.h>
#endif
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


typedef struct _ThreadArgs
{
	const Image *src;
	Image *dst;
	int fwd, loud, x1, x2, y1, y2;
	
	int pixels[(BLOCKDX+16)*2*4*2];//2 padded rows * 4 channels * {pixel, error}

	int sse[4][2048];

	//AC
	DList list;
	const unsigned char *decstart, *decend;

	//aux
	int blockidx;
} ThreadArgs;
static void block_thread(void *param)
{
	GolombRiceCoder ec;
	ThreadArgs *args=(ThreadArgs*)param;
	Image const *image=args->fwd?args->src:args->dst;

	if(args->fwd)
	{
		dlist_init(&args->list, 1, 1024, 0);
		gr_enc_init(&ec, &args->list);
	}
	else//decode
		gr_dec_init(&ec, args->decstart, args->decend);
	memset(args->sse, 0, sizeof(args->sse));
	args->pixels[0]=0;
	args->pixels[1]=0;
	args->pixels[2]=0;
	args->pixels[3]=0;
	args->pixels[4]=1<<(image->depth>>1);
	args->pixels[5]=1<<(image->depth>>1);
	args->pixels[6]=1<<(image->depth>>1);
	args->pixels[7]=1<<(image->depth>>1);
	memfill(args->pixels+8, args->pixels, sizeof(args->pixels)-sizeof(int[8]), sizeof(int[8]));
	//memset(args->pixels, 0, sizeof(args->pixels));
	for(int ky=args->y1;ky<args->y2-1;ky+=2)//codec loop
	{
		int *rows[]=
		{
			args->pixels+((BLOCKDX+16LL)*((ky-0LL)&1)+8)*4*2,
			args->pixels+((BLOCKDX+16LL)*((ky-1LL)&1)+8)*4*2,
		};
		int perm[]={1, 2, 3, 0};
		int comp[4]={0};
		int pred=0, nbypass=0, val=0, error=0;
		for(int kx=args->x1;kx<args->x2-1;kx+=2)
		{
			int
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*NEE	=rows[1]+2*4*2,
				*NEEE	=rows[1]+3*4*2,
				*NEEEE	=rows[1]+4*4*2,
				*WWWW	=rows[0]-4*4*2,
				*WWW	=rows[0]-3*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
			//if(ky==222&&kx==468)//
			//if(ky==2&&kx==2)//
			//if(ky==44&&kx==4820)//
			//	printf("");
			if(args->fwd)
			{
				comp[0]=image->data[image->iw*(ky+0)+kx+0];//r
				comp[1]=image->data[image->iw*(ky+0)+kx+1];//g
				comp[2]=image->data[image->iw*(ky+1)+kx+0];//g
				comp[3]=image->data[image->iw*(ky+1)+kx+1];//b

				comp[0]-=comp[1];	//RCT3
				comp[1]+=comp[0]>>1;
				comp[2]-=comp[1];
				comp[1]+=comp[2]>>1;
				comp[3]-=comp[1];
				comp[1]+=comp[3]>>1;
				
				//comp[0]-=comp[1];	//RCT2
				//comp[2]-=comp[3];
				//comp[1]+=comp[0]>>1;
				//comp[3]+=comp[2]>>1;
				//comp[0]-=comp[2];
				//comp[1]-=comp[3];
				//comp[2]+=comp[0]>>1;
				//comp[3]+=comp[1]>>1;

				//comp[0]-=comp[1];	//RCT1
				//comp[2]-=comp[1];
				//comp[3]-=comp[1];
				//comp[1]+=(comp[0]+comp[2]+comp[3])/3;
			}
			for(int kc=0;kc<4;++kc)
			{
				//if(W[kc+4]<0)//
				//	LOG_ERROR("");
				int kc2=perm[kc];
				MEDIAN3_32(pred, N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
				nbypass=FLOOR_LOG2(W[kc+4]+1);
				if(args->fwd)
				{
					val=comp[kc2];
					error=val-pred;
					error=error<<1^error>>31;
					gr_enc_POT(&ec, error, nbypass);
				}
				else
				{
					error=gr_dec_POT(&ec, nbypass);
					val=error;
					val=val>>1^-(val&1);
					val+=pred;
					comp[kc2]=val;
				}
				curr[kc]=comp[kc2];
				curr[kc+4]=(WWW[kc+4]+5*WW[kc+4]+10*W[kc+4]+10*error+NE[kc+4])/27;	//Formula5
			//	curr[kc+4]=(WW[kc+4]+3*(W[kc+4]+error)+N[kc+4]+NEEE[kc+4])/9;
			//	curr[kc+4]=(3*WWW[kc+4]+11*WW[kc+4]+23*W[kc+4]+23*error+2*NE[kc+4]+2*NEEEE[kc+4])>>6;
			//	curr[kc+4]=(WW[kc+4]+7*(W[kc+4]+error)+NEEE[kc+4])>>4;	//Formula4
			//	curr[kc+4]=(WWW[kc+4]+WW[kc+4]+14*W[kc+4]+13*error+NEEE[kc+4]+2*NEEEE[kc+4])>>5;
			//	curr[kc+4]=(8*W[kc+4]+7*error+NEEE[kc+4])>>4;	//Formula3
			//	curr[kc+4]=(WW[kc+4]+7*W[kc+4]+6*error+NEEE[kc+4]+NEEEE[kc+4])>>4;
			//	curr[kc+4]=(4*W[kc+4]+3*error+NEEE[kc+4])>>3;	//Formula2
			//	curr[kc+4]=(2*W[kc+4]+error+NEEE[kc+4])>>2;	//Formula1
			}
			if(!args->fwd)
			{
				comp[1]-=comp[3]>>1;
				comp[3]+=comp[1];
				comp[1]-=comp[2]>>1;
				comp[2]+=comp[1];
				comp[1]-=comp[0]>>1;
				comp[0]+=comp[1];

				//comp[3]-=comp[1]>>1;
				//comp[2]-=comp[0]>>1;
				//comp[1]+=comp[3];
				//comp[0]+=comp[2];
				//comp[3]-=comp[2]>>1;
				//comp[1]-=comp[0]>>1;
				//comp[2]+=comp[3];
				//comp[0]+=comp[1];
				
				//comp[1]-=(comp[0]+comp[2]+comp[3])/3;
				//comp[3]+=comp[1];
				//comp[2]+=comp[1];
				//comp[0]+=comp[1];
				
				args->dst->data[image->iw*(ky+0)+kx+0]=comp[0];//r
				args->dst->data[image->iw*(ky+0)+kx+1]=comp[1];//g
				args->dst->data[image->iw*(ky+1)+kx+0]=comp[2];//g
				args->dst->data[image->iw*(ky+1)+kx+1]=comp[3];//b
#ifdef ENABLE_GUIDE
				if(
					memcmp(image->data+image->iw*(ky+0)+kx, guide->data+image->iw*(ky+0)+kx, sizeof(short[2]))||
					memcmp(image->data+image->iw*(ky+1)+kx, guide->data+image->iw*(ky+1)+kx, sizeof(short[2]))
				)
				{
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
			}
			

			rows[0]+=4*2;
			rows[1]+=4*2;
		}
	}
	if(args->fwd)
		gr_enc_flush(&ec);
}
int f32_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	int xblocks, yblocks, nblocks, ncores, nthreads, coffset;
	ptrdiff_t start=0;
	ptrdiff_t memusage=0;
	ptrdiff_t argssize;
	ThreadArgs *args;
	size_t usize;
	
	//if(fwd)
	//	image_snapshot8(image);
	ncores=query_cpu_cores();
	usize=((size_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
	xblocks=(image->iw+BLOCKDX-1)/BLOCKDX;
	yblocks=(image->ih+BLOCKDY-1)/BLOCKDY;
	nblocks=xblocks*yblocks, nthreads=MINVAR(nblocks, ncores);
	coffset=sizeof(int)*nblocks;
	argssize=sizeof(ThreadArgs)*nthreads;
	args=(ThreadArgs*)malloc(argssize);
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
				int blocksize=(arg->y2-arg->y1)*image->iw*image->nch*image->depth/8;
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
			printf("Mem usage: ");
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
			printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
		}
		printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
	}
	free(args);
	return 0;
}