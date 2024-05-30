#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
//#ifdef _MSC_VER
//#include<intrin.h>
//#else
//#include<x86intrin.h>
//#endif
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define DISABLE_MT


#include"ac.h"
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif

#define BLOCKSIZE 256

typedef enum RCTTypeEnum
{
	RCT_NONE,
	RCT_SUBG,
	RCT_JPEG2000,
	RCT_1,
	RCT_PEI09,
} RCTType;
typedef struct _ThreadArgs
{
	const Image *src;
	Image *dst;
	int fwd, loud, y1, y2;
	int bufsize, histsize;
	int *pixels, *hist;

	DList list;
	const unsigned char *decstart, *decend;
} ThreadArgs;
static const int combinations[]=
{
	 0,  1,  2,
	 1,  3,  4,
	 5,  6,  7,
	 8,  9,  7,
	10, 11,  7,
};
static const short wp_params[]=//signed fixed 7.8 bit
{
	0x0DB8,  0x0E22,  0x181F,  0x0BF3, -0x005C, -0x005B,  0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,//Y
	0x064C,  0x0F31,  0x1040,  0x0BF8, -0x0007, -0x000D, -0x0085, -0x0063, -0x00A2, -0x0017,  0x00F2,//Cb
	0x0B37,  0x110B,  0x121B,  0x0BFC, -0x0001,  0x000E, -0x0188, -0x00E7, -0x00BB, -0x004A,  0x00BA,//Cr
};
static int wp_predict(const short *params, int NN, int NW, int N, int NE, int W, const int *eNW, const int *eN, const int *eNE, const int *eW, int *preds)
{
	int pred;
	long long lpred=0, wsum=0;
	preds[0]=(W+NE-N)<<8;
	preds[1]=(N<<8)-((eN[0]+eW[0]+eNE[0])*params[4]>>8);
	preds[2]=(W<<8)-((eN[0]+eW[0]+eNW[0])*params[5]>>8);
	preds[3]=(N<<8)-(((eNW[0]*params[6]+eN[0]*params[7]+eNE[0]*params[8])>>8)+(NN-N)*params[9]+(NW-W)*params[10]);
	for(int k=0;k<4;++k)
	{
		int weight=(params[k]<<8)/(eNW[k+1]+eN[k+1]+eNE[k+1]+1);
		lpred+=(long long)weight*preds[k];
		wsum+=weight;
	}
	lpred/=wsum+1;
	if(!wsum)
		lpred=preds[0];
	CLAMP3_32(pred, (int)lpred, N<<8, W<<8, NE<<8);
	return pred;
}
static void wp_update(int curr, int pred, const int *preds, int *ecurr, int *eNE)
{
	curr<<=8;
	ecurr[0]=curr-(int)pred;
	for(int k=0;k<4;++k)
	{
		int e2=abs(curr-preds[k]);
		ecurr[k+1]=e2;
		eNE[k+1]+=e2;
	}
}
static void block_enc(void *param)
{
	const short *ch_params[]=
	{
		wp_params+0*11,
		wp_params+1*11,
		wp_params+2*11,

		wp_params+0*11,
		wp_params+2*11,

		wp_params+1*11,
		wp_params+2*11,
		wp_params+0*11,

		wp_params+1*11,
		wp_params+2*11,

		wp_params+1*11,
		wp_params+2*11,
	};
	GolombRiceCoder ec;
	ThreadArgs *args=(ThreadArgs*)param;
	Image const *image=args->src;
	//int nlevels=1<<image->depth, half=nlevels>>1, mask=nlevels-1;
	int depths[]=
	{
		image->depth,
		image->depth,
		image->depth,

		image->depth,
		image->depth,
		
		image->depth,
		image->depth+1,
		image->depth+1,
		
		image->depth,
		image->depth+1,
		
		image->depth,
		image->depth+1,
	};
	int nlevels[12]={0}, halfs[12]={0};
	double csizes[24]={0};
	int res=image->iw*(args->y2-args->y1);

	for(int k=0;k<12;++k)
	{
		nlevels[k]=1<<depths[k];
		halfs[k]=nlevels[k]>>1;
	}
	memset(args->pixels, 0, args->bufsize);
	memset(args->hist, 0, args->histsize);
	for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
	{
		ALIGN(16) int *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*72,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*72,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*72,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*72,
		};
		int comp[12]={0};
		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			int
				*NN	=rows[2]+0*72,
				*NW	=rows[1]-1*72,
				*N	=rows[1]+0*72,
				*NE	=rows[1]+1*72,
				*W	=rows[0]-1*72,
				*curr	=rows[0]+0*72;

			//0	r-P(rprev)		NONE		r; g; b;
			//1	g-P(gprev)
			//2	b-P(bprev)
			//
			//3	r-clamp(P(rprev-gprev)+gcurr)	SubG	r-=g; b-=g; g;
			//4	b-clamp(P(bprev-gprev)+gcurr)
			//
			//5	y1-P(y1prev)		JPEG2000	r-=g; b-=g; g+=(r+g)>>2;
			//6	u1-P(u1prev)
			//7	v1-P(v1prev)
			//
			//8	y2-P(y2prev)		RCT1		r-=g; g+=r>>1; b-=g; g+=b>>1;
			//9	u2-P(u2prev)
			//
			//10	y3-P(y3prev)		Pei09		b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
			//11	u3-P(u3prev)

			//NONE
			comp[0]=image->data[idx+0];//r
			comp[1]=image->data[idx+1];//g
			comp[2]=image->data[idx+2];//b
			
			//SubG+
			comp[3]=comp[0];//cr0
			comp[4]=comp[2];//cb0

			//JPEG2000_RCT
			comp[7]=comp[0]-comp[1];//cr1		cr1 reused in RCT1 & Pei09
			comp[6]=comp[2]-comp[1];//cb1
			comp[5]=comp[1]+((comp[7]+comp[6])>>2);//y1

			//RCT1
			comp[8]=comp[1]+(comp[7]>>1);
			comp[9]=comp[2]-comp[8];//cb2
			comp[8]+=comp[9]>>1;//y2

			//Pei09
			comp[11]=comp[2]-((87*comp[0]+169*comp[1]+128)>>8);//cb3
			comp[10]=comp[1]+((86*comp[7]+29*comp[11]+128)>>8);//y3
			
			//if(ky==10&&kx==10)//
			//	printf("");

			for(int kc=0;kc<12;++kc)
			{
				int kc2=kc*6;
				int pred, preds[4], cgrad, offset=0, val;
				
				pred=wp_predict(ch_params[kc], NN[kc2], NW[kc2], N[kc2], NE[kc2], W[kc2], NW+kc2+1, N+kc2+1, NE+kc2+1, W+kc2+1, preds);
				//const short *params=ch_params[kc];
				//int preds[]=
				//{
				//	(W[kc2]+NE[kc2]-N[kc2])<<8,
				//	(N[kc2]<<8)-((N[kc2+1]+W[kc2+1]+NE[kc2+1])*params[4]>>8),
				//	(W[kc2]<<8)-((N[kc2+1]+W[kc2+1]+NW[kc2+1])*params[5]>>8),
				//	(N[kc2]<<8)-(((NW[kc2+1]*params[6]+N[kc2+1]*params[7]+NE[kc2+1]*params[8])>>8)+(NN[kc2]-N[kc2])*params[9]+(NW[kc2]-W[kc2])*params[10]),
				//};
				//int pred, cgrad, offset=0;
				//long long lpred=0, wsum=0;
				//
				//for(int k=0;k<4;++k)
				//{
				//	int weight=(params[k]<<8)/(NW[kc2+1+k]+N[kc2+1+k]+NE[kc2+1+k]+1);
				//	lpred+=(long long)weight*preds[k];
				//	wsum+=weight;
				//}
				//lpred/=wsum+1;
				//if(!wsum)
				//	lpred=preds[0];
				//CLAMP3_32(pred, (int)lpred, N[kc2]<<8, W[kc2]<<8, NE[kc2]<<8);
				MEDIAN3_32(cgrad, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);

				if((unsigned)(kc-3)<2)
				{
					offset=comp[1];
					pred+=offset;
					cgrad+=offset;
					CLAMP2_32(pred, pred, -halfs[0], halfs[0]-1);
					CLAMP2_32(cgrad, cgrad, -halfs[0], halfs[0]-1);
				}
				val=comp[kc]-cgrad;
				val+=halfs[kc];
				val&=nlevels[kc]-1;
				++args->hist[(kc<<1|0)<<(image->depth+1)|val];

				val=comp[kc]-((pred+127)>>8);
				val+=halfs[kc];
				val&=nlevels[kc]-1;
				++args->hist[(kc<<1|1)<<(image->depth+1)|val];

				curr[kc2+0]=comp[kc]-offset;
				wp_update(comp[kc], pred, preds, curr+kc2+1, NE+kc2+1);
				//curr[kc2+1]=(comp[kc]<<8)-(int)lpred;
				//for(int k=0;k<4;++k)
				//{
				//	int e2=abs((comp[kc]<<8)-preds[k]);
				//	curr[kc2+1+k]=e2;
				//	NE[kc2+1+k]+=e2;
				//}
			}
			rows[0]+=72;
			rows[1]+=72;
			rows[2]+=72;
			rows[3]+=72;
		}
	}
	for(int kc=0;kc<24;++kc)
	{
		int *curr_hist=args->hist+((size_t)kc<<(image->depth+1));
		for(int ks=0;ks<nlevels[kc>>1];++ks)
		{
			int freq=curr_hist[ks];
			if(freq)
				csizes[kc]-=freq*log2((double)freq/res);
		}
		csizes[kc]/=8;
	}
	int bestrct=0;
	char predsel[12]={0};
	double bestsize=0;
	for(int k=0;k<24;k+=2)
		predsel[k>>1]=csizes[k|1]<csizes[k|0];
	const int *group=combinations;
	bestsize=csizes[group[0]<<1|predsel[group[0]]]+csizes[group[1]<<1|predsel[group[1]]]+csizes[group[2]<<1|predsel[group[2]]];
	for(int k=0;k<(int)_countof(combinations)/3;++k)
	{
		group=combinations+k*3;
		double csize=csizes[group[0]<<1|predsel[group[0]]]+csizes[group[1]<<1|predsel[group[1]]]+csizes[group[2]<<1|predsel[group[2]]];
		if(bestsize>csize)
			bestsize=csize, bestrct=k;
	}
	int combination[]=
	{
		combinations[bestrct*3+0]<<1|predsel[combinations[bestrct*3+0]],
		combinations[bestrct*3+1]<<1|predsel[combinations[bestrct*3+1]],
		combinations[bestrct*3+2]<<1|predsel[combinations[bestrct*3+2]],
	};
	int flag=bestrct<<3|(combination[2]&1)<<2|(combination[1]&1)<<1|(combination[0]&1);
	if(args->loud)
	{
		static const char *rctnames[]=
		{
			"NONE",
			"SubG",
			"JPEG2000",
			"RCT1",
			"Pei09",
		};
		static const char *prednames[]=
		{
			"CG",
			"WP",
		};

		double defsize=csizes[11]+csizes[13]+csizes[15];
		printf("Y %5d~%5d  default %lf (%+lf) bytes  current %lf bytes  %s [%s %s %s]\n",
			args->y1, args->y2,
			defsize, bestsize-defsize,
			bestsize,
			rctnames[bestrct],
			prednames[combination[0]&1],
			prednames[combination[1]&1],
			prednames[combination[2]&1]
		);
	}
	dlist_init(&args->list, 1, 1024, 0);
	dlist_push_back(&args->list, &flag, 1);
	gr_enc_init(&ec, &args->list);
	memset(args->pixels, 0, args->bufsize);
	for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
	{
		static const int perm[]={1, 2, 0};

		ALIGN(16) int *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*28,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*28,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*28,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*28,
		};
		int comp[4]={0};
		int preds[4]={0};
		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			int
				*NN	=rows[2]+0*28,
				*NW	=rows[1]-1*28,
				*N	=rows[1]+0*28,
				*NE	=rows[1]+1*28,
				*NEEE	=rows[1]+3*28,
				*W	=rows[0]-1*28,
				*curr	=rows[0]+0*28;
			comp[0]=image->data[idx+0];
			comp[1]=image->data[idx+1];
			comp[2]=image->data[idx+2];
			switch(bestrct)
			{
			case RCT_NONE:
			case RCT_SUBG:
				break;
			case RCT_JPEG2000:
				comp[0]-=comp[1];
				comp[2]-=comp[1];
				comp[1]+=(comp[0]+comp[2])>>2;
				break;
			case RCT_1:
				comp[0]-=comp[1];
				comp[1]+=comp[0]>>1;
				comp[2]-=comp[1];
				comp[1]+=comp[2]>>1;
				break;
			case RCT_PEI09:
				comp[2]-=(87*comp[0]+169*comp[1]+128)>>8;
				comp[0]-=comp[1];
				comp[1]+=(86*comp[0]+29*comp[2]+128)>>8;
				break;
			}

			//if(ky==330&&kx==146)//
			//	printf("");

			for(int kc0=0;kc0<3;++kc0)
			{
				int kc=perm[kc0], kc2=kc*7, p0=0, pred, offset=0, ch=combination[kc0]>>1, val, sym;
				if(combination[kc0]&1)//WP
				{
					p0=wp_predict(ch_params[ch], NN[kc2], NW[kc2], N[kc2], NE[kc2], W[kc2], NW+kc2+2, N+kc2+2, NE+kc2+2, W+kc2+2, preds);
					pred=(p0+127)>>8;
				}
				else//CG
				{
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
				}
				if(bestrct==RCT_SUBG&&kc0>0)
				{
					offset=rows[0][1];
					pred+=offset;
					CLAMP2_32(pred, pred, -halfs[ch], halfs[ch]-1);
				}
				val=comp[kc]-pred;
				val<<=32-depths[ch];
				val>>=32-depths[ch];
				sym=val<<1^(val>>31);
				gr_enc_POT(&ec, sym, FLOOR_LOG2(W[kc2+1]+1));
				curr[kc2+0]=comp[kc]-offset;
				curr[kc2+1]=(2*W[kc2+1]+sym+NEEE[kc2+1])>>2;
				if(combination[kc0]&1)//WP
					wp_update(curr[kc2+0], p0, preds, curr+kc2+2, NE+kc2+2);
			}
			rows[0]+=28;
			rows[1]+=28;
			rows[2]+=28;
			rows[3]+=28;
		}
	}
	gr_enc_flush(&ec);
}
static void block_dec(void *param)
{
	const short *ch_params[]=
	{
		wp_params+0*11,
		wp_params+1*11,
		wp_params+2*11,

		wp_params+0*11,
		wp_params+2*11,

		wp_params+1*11,
		wp_params+2*11,
		wp_params+0*11,

		wp_params+1*11,
		wp_params+2*11,

		wp_params+1*11,
		wp_params+2*11,
	};
	GolombRiceCoder ec;
	ThreadArgs *args=(ThreadArgs*)param;
	Image *image=args->dst;
	const unsigned char *srcstart=args->decstart, *srcend=args->decend;
	int flag=0, bestrct;
	int combination[3]={0};
	int depths[]=
	{
		image->depth,
		image->depth,
		image->depth,

		image->depth,
		image->depth,
		
		image->depth,
		image->depth+1,
		image->depth+1,
		
		image->depth,
		image->depth+1,
		
		image->depth,
		image->depth+1,
	};
	int halfs[12]={0};
	
	for(int k=0;k<12;++k)
		halfs[k]=1<<depths[k]>>1;
	memset(args->pixels, 0, args->bufsize);
	memcpy(&flag, srcstart, sizeof(char));
	srcstart+=sizeof(char);
	bestrct=flag>>3;
	memcpy(combination, combinations+3*bestrct, sizeof(int[3]));
	for(int k=0;k<3;++k)
		combination[k]=combination[k]<<1|(flag>>k&1);
	gr_dec_init(&ec, srcstart, srcend);
	for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
	{
		static const int perm[]={1, 2, 0};

		ALIGN(16) int *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*28,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*28,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*28,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*28,
		};
		int comp[4]={0};
		int preds[4]={0};
		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			int
				*NN	=rows[2]+0*28,
				*NW	=rows[1]-1*28,
				*N	=rows[1]+0*28,
				*NE	=rows[1]+1*28,
				*NEEE	=rows[1]+3*28,
				*W	=rows[0]-1*28,
				*curr	=rows[0]+0*28;

			//if(ky==330&&kx==146)//
			//	printf("");

			for(int kc0=0;kc0<3;++kc0)
			{
				int kc=perm[kc0], kc2=kc*7, p0=0, pred, offset=0, ch=combination[kc0]>>1, val, sym;
				if(combination[kc0]&1)//WP
				{
					p0=wp_predict(ch_params[ch], NN[kc2], NW[kc2], N[kc2], NE[kc2], W[kc2], NW+kc2+2, N+kc2+2, NE+kc2+2, W+kc2+2, preds);
					pred=(p0+127)>>8;
				}
				else//CG
				{
					MEDIAN3_32(pred, N[kc2], W[kc2], N[kc2]+W[kc2]-NW[kc2]);
				}
				if(bestrct==RCT_SUBG&&kc0>0)
				{
					offset=rows[0][1];
					pred+=offset;
					CLAMP2_32(pred, pred, -halfs[ch], halfs[ch]-1);
				}
				sym=gr_dec_POT(&ec, FLOOR_LOG2(W[kc2+1]+1));
				val=sym>>1^-(sym&1);
				val+=pred;
				val<<=32-depths[ch];
				val>>=32-depths[ch];
				comp[kc]=val;
				curr[kc2+0]=comp[kc]-offset;
				curr[kc2+1]=(2*W[kc2+1]+sym+NEEE[kc2+1])>>2;
				if(combination[kc0]&1)//WP
					wp_update(curr[kc2+0], p0, preds, curr+kc2+2, NE+kc2+2);
			}
			switch(bestrct)
			{
			case RCT_NONE:
			case RCT_SUBG:
				break;
			case RCT_JPEG2000:
				comp[1]-=(comp[0]+comp[2])>>2;
				comp[2]+=comp[1];
				comp[0]+=comp[1];
				break;
			case RCT_1:
				comp[1]-=comp[2]>>1;
				comp[2]+=comp[1];
				comp[1]-=comp[0]>>1;
				comp[0]+=comp[1];
				break;
			case RCT_PEI09:
				comp[1]-=(86*comp[0]+29*comp[2]+128)>>8;
				comp[0]+=comp[1];
				comp[2]+=(87*comp[0]+169*comp[1]+128)>>8;
				break;
			}
			image->data[idx+0]=comp[0];
			image->data[idx+1]=comp[1];
			image->data[idx+2]=comp[2];
#ifdef ENABLE_GUIDE
			if(memcmp(image->data+idx, guide->data+idx, sizeof(short)*image->nch))
			{
				short orig[4]={0};
				memcpy(orig, guide->data+idx, image->nch*sizeof(short));
				LOG_ERROR("Guide error XY %d %d", kx, ky);
				printf("");//
			}
#endif
			rows[0]+=28;
			rows[1]+=28;
			rows[2]+=28;
			rows[3]+=28;
		}
	}
}
int f24_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	if(image->nch!=3)
	{
		LOG_ERROR("Expected 3 channels, got %d", image->nch);
		return 1;
	}
#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif
	int ncores=query_cpu_cores();
	int nblocks=(image->ih+BLOCKSIZE-1)/BLOCKSIZE, nthreads=MINVAR(nblocks, ncores);
	ptrdiff_t memusage=0;
	ptrdiff_t argssize=nthreads*sizeof(ThreadArgs);
	ThreadArgs *args=(ThreadArgs*)malloc(argssize);
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memusage+=argssize;
	memset(args, 0, argssize);
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=src;
		arg->dst=dst;
		arg->bufsize=sizeof(int[4*12*6])*(image->iw+16LL);//4 padded rows * 4 channels max * {pixels, errors}
		arg->pixels=(int*)_mm_malloc(arg->bufsize, sizeof(__m128i));
		arg->histsize=sizeof(int[24])<<(image->depth+1);
		arg->hist=(int*)malloc(arg->histsize);
		if(!arg->pixels||!arg->hist)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=(ptrdiff_t)arg->bufsize+arg->histsize;
		arg->fwd=fwd;
#ifdef DISABLE_MT
		arg->loud=loud;
#else
		arg->loud=0;
#endif
	}
	if(fwd)
	{
		ptrdiff_t dststart=array_append(data, 0, 1, sizeof(int)*nblocks, 1, 0, 0);
		for(int kt=0;kt<nblocks;kt+=nthreads)
		{
			int nthreads2=MINVAR(kt+nthreads, nblocks)-kt;
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				arg->y1=BLOCKSIZE*(kt+kt2);
				arg->y2=MINVAR(arg->y1+BLOCKSIZE, image->ih);
			}
#ifdef DISABLE_MT
			for(int k=0;k<nthreads2;++k)
				block_enc(args+k);
#else
			void *ctx=mt_exec(block_enc, args, sizeof(ThreadArgs), nthreads2);
			mt_finish(ctx);
#endif
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				if(loud)
					printf("[%d]  %zd\n", kt+kt2, args[kt2].list.nobj);
				memcpy(data[0]->data+dststart+sizeof(int)*((ptrdiff_t)kt+kt2), &args[kt2].list.nobj, sizeof(int));
				dlist_appendtoarray(&args[kt2].list, data);
				dlist_clear(&args[kt2].list);
			}
		}
		if(loud)
		{
			ptrdiff_t
				csize=data[0]->count-dststart,
				usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
			t0=time_sec()-t0;
			printf("Size %14td/%14td  %16lf%%  %16lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
			printf("Mem usage: ");
			print_size((double)memusage, 8, 4, 0, 0);
			printf("\n");
			printf("E %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
	}
	else
	{
		const unsigned char *dstptr=cbuf+sizeof(int)*nblocks;
		int dec_offset=0;

		//integrity check
#if 1
		for(int kt=0;kt<nblocks;++kt)
		{
			int size=0;
			memcpy(&size, cbuf+sizeof(int)*kt, sizeof(int));
			dec_offset+=size;
		}
		if(sizeof(int)*nblocks+dec_offset!=clen)
			LOG_ERROR("Corrupt file");
#endif
		dec_offset=0;
		for(int kt=0;kt<nblocks;kt+=nthreads)
		{
			int nthreads2=MINVAR(kt+nthreads, nblocks)-kt;
			for(int kt2=0;kt2<nthreads2;++kt2)
			{
				ThreadArgs *arg=args+kt2;
				int size=0;
				memcpy(&size, cbuf+sizeof(int)*((ptrdiff_t)kt+kt2), sizeof(int));
				arg->y1=BLOCKSIZE*(kt+kt2);
				arg->y2=MINVAR(arg->y1+BLOCKSIZE, image->ih);
				arg->decstart=dstptr+dec_offset;
				dec_offset+=size;
				arg->decend=dstptr+dec_offset;
			}
#ifdef DISABLE_MT
			for(int k=0;k<nthreads2;++k)
				block_dec(args+k);
#else
			{
				void *ctx=mt_exec(block_dec, args, sizeof(ThreadArgs), nthreads2);
				mt_finish(ctx);
			}
#endif
		}
		if(loud)
		{
			ptrdiff_t usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
			t0=time_sec()-t0;
			printf("D %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		_mm_free(arg->pixels);
		free(arg->hist);
	}
	free(args);
	return 0;
}