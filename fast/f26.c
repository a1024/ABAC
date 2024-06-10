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
#if 0
static int cgrad(int N, int W, int NW)
{
	int pred;

	MEDIAN3_32(pred, N, W, N+W-NW);
	return pred;
}
static int clampav(int NW, int N, int NE, int WW, int W)
{
	ALIGN(16) int pred[4];
	__m128i va=_mm_set_epi32(0, 0, 0, N);
	__m128i vb=_mm_set_epi32(0, 0, 0, W);
	__m128i vc=_mm_set_epi32(0, 0, 0, NE);
	__m128i vd=_mm_set_epi32(0, 0, 0, (8*W+5*(N-NW)+NE-WW)>>3);
	__m128i vmin=_mm_min_epi32(va, vb);
	__m128i vmax=_mm_max_epi32(va, vb);
	vmin=_mm_min_epi32(vmin, vc);
	vmax=_mm_max_epi32(vmax, vc);
	vd=_mm_max_epi32(vd, vmin);
	vd=_mm_min_epi32(vd, vmax);
	_mm_store_si128((__m128i*)pred, vd);
	return pred[0];
}
static int clamp(int vmin, int x, int vmax)
{
	int ret;

	MEDIAN3_32(ret, vmin, vmax, x);
	return ret;
}
#endif

typedef enum _IChannelType
{
	ICH_ZERO,
	ICH_R,
	ICH_G,
	ICH_B,
	ICH_J2Y,
	ICH_J2U,
	ICH_J2V,
	ICH_R1Y,
	ICH_R1U,
	ICH_P9Y,
	ICH_P9U,

	ICH_COUNT,
} IChannelType;
#define ICH_R1V ICH_J2V
#define ICH_P9V ICH_J2V
#define OCHLIST\
	OCH(OCH_R,	0,	ICH_R,		ICH_ZERO)\
	OCH(OCH_G,	0,	ICH_G,		ICH_ZERO)\
	OCH(OCH_B,	0,	ICH_B,		ICH_ZERO)\
	OCH(OCH_RG,	0,	ICH_R,		ICH_G)\
	OCH(OCH_RB,	0,	ICH_R,		ICH_B)\
	OCH(OCH_GR,	0,	ICH_G,		ICH_R)\
	OCH(OCH_GB,	0,	ICH_G,		ICH_B)\
	OCH(OCH_BG,	0,	ICH_B,		ICH_G)\
	OCH(OCH_BR,	0,	ICH_B,		ICH_R)\
	OCH(OCH_J2Y,	0,	ICH_J2Y,	ICH_ZERO)\
	OCH(OCH_J2U,	1,	ICH_J2U,	ICH_ZERO)\
	OCH(OCH_J2V,	1,	ICH_J2V,	ICH_ZERO)\
	OCH(OCH_R1Y,	0,	ICH_R1Y,	ICH_ZERO)\
	OCH(OCH_R1U,	1,	ICH_R1U,	ICH_ZERO)\
	OCH(OCH_P9Y,	0,	ICH_P9Y,	ICH_ZERO)\
	OCH(OCH_P9U,	1,	ICH_P9U,	ICH_ZERO)
typedef enum _OChannelType
{
#define OCH(ONAME, OINF, TARGET, HELPER) ONAME,
	OCHLIST
#undef  OCH
	OCH_COUNT,
} OChannelType;
static const char och_inflation[OCH_COUNT]=
{
#define OCH(ONAME, OINF, TARGET, HELPER) OINF,
	OCHLIST
#undef  OCH
};
static const int och_dependencies[OCH_COUNT][2]=
{
#define OCH(ONAME, OINF, TARGET, HELPER) {TARGET, HELPER},
	OCHLIST
#undef  OCH
};
#define RCTLIST\
	RCT(R_G_B,	OCH_R,		OCH_G,		OCH_B)\
	RCT(R_G_BG,	OCH_R,		OCH_G,		OCH_BG)\
	RCT(R_G_BR,	OCH_R,		OCH_G,		OCH_BR)\
	RCT(G_B_RG,	OCH_G,		OCH_B,		OCH_RG)\
	RCT(G_B_RB,	OCH_G,		OCH_B,		OCH_RB)\
	RCT(B_R_GR,	OCH_B,		OCH_R,		OCH_GR)\
	RCT(B_R_GB,	OCH_B,		OCH_R,		OCH_GB)\
	RCT(G_BG_RG,	OCH_G,		OCH_BG,		OCH_RG)\
	RCT(G_BG_RB,	OCH_G,		OCH_BG,		OCH_RB)\
	RCT(G_RG_BR,	OCH_G,		OCH_RG,		OCH_BR)\
	RCT(B_RB_GB,	OCH_B,		OCH_RB,		OCH_GB)\
	RCT(B_RB_GR,	OCH_B,		OCH_RB,		OCH_GR)\
	RCT(B_GB_RG,	OCH_B,		OCH_GB,		OCH_RG)\
	RCT(G_GR_BR,	OCH_R,		OCH_GR,		OCH_BR)\
	RCT(R_GR_BG,	OCH_R,		OCH_GR,		OCH_BG)\
	RCT(R_BR_GB,	OCH_R,		OCH_BR,		OCH_GB)\
	RCT(J2K,	OCH_J2Y,	OCH_J2U,	OCH_J2V)\
	RCT(RCT1,	OCH_R1Y,	OCH_R1U,	OCH_J2V)\
	RCT(Pei09,	OCH_P9Y,	OCH_P9U,	OCH_J2V)
typedef enum _RCTType
{
#define RCT(RCTNAME, YIDX, UIDX, VIDX) RCT_##RCTNAME,
	RCTLIST
#undef  RCT
	RCT_COUNT,
} RCTType;
static const int rctcombinations[RCT_COUNT][3]=
{
#define RCT(RCTNAME, YIDX, UIDX, VIDX) {YIDX, UIDX, VIDX},
	RCTLIST
#undef  RCT
};

#define PREDLIST\
	PRED(W)\
	PRED(cgrad)\
	PRED(wgrad)
typedef enum _PredType
{
#define PRED(NAME) PRED_##NAME,
	PREDLIST
#undef  PRED
	PRED_COUNT,
} PredType;


//WG:

//	#define WG_UPDATE
#define WG_RESCALE_LIMIT 100
#define WG_DECAY_NUM	493
#define WG_DECAY_SH	9

#define WG_NPREDS	8
#define WG_PREDLIST\
	WG_PRED(1.2, N+W-NW)\
	WG_PRED(1.5, N)\
	WG_PRED(1.5, W)\
	WG_PRED(1, W+NE-N)\
	WG_PRED(1, N+NE-NNE)\
	WG_PRED(1, 3*(N-NN)+NNN)\
	WG_PRED(1, 3*(W-WW)+WWW)\
	WG_PRED(1, (W+NEEE)/2)
//	WG_PRED(0.5, NW)
//	WG_PRED(0.5, NE)
static void wg_init(double *weights)
{
	int j=0;
#define WG_PRED(WEIGHT, EXPR) weights[j++]=WEIGHT;
	WG_PREDLIST
#undef  WG_PRED
}
static int wg_predict(const double *weights, int NNN, int NN, int NNE, int NW, int N, int NE, int NEE, int NEEE, int WWW, int WW, int W, const int *eW, int *preds)
{
	int pred;
	double pred2=0, wsum=0;
	int j=0;
#define WG_PRED(WEIGHT, EXPR) preds[j++]=EXPR;
	WG_PREDLIST
#undef  WG_PRED
	
	for(int k=0;k<WG_NPREDS;++k)
	{
		double weight=weights[k]/(eW[k]+1);
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
static void wg_update(int curr, const int *preds, int *eW, double *weights)
{
	int kbest=0, ebest=0;
	for(int k=0;k<WG_NPREDS;++k)
	{
		int e2=abs(curr-preds[k]);
		eW[k]=(eW[k]+e2)*WG_DECAY_NUM>>WG_DECAY_SH;
		if(!k||ebest>e2)
			kbest=k, ebest=e2;
	}
#ifdef WG_UPDATE
	++weights[kbest];
	if(weights[kbest]>WG_RESCALE_LIMIT)
	{
		for(int k=0;k<WG_NPREDS;++k)
			weights[k]*=0.5;
	}
#endif
}


typedef struct _ThreadArgs
{
	const Image *src;
	Image *dst;
	int fwd, loud, y1, y2;
	int bufsize, histsize;
	int *pixels, *hist;

	DList list;
	const unsigned char *decstart, *decend;

	double wg_weights[WG_NPREDS*OCH_COUNT];
} ThreadArgs;
static void block_thread(void *param)
{
	GolombRiceCoder ec;
	ThreadArgs *args=(ThreadArgs*)param;
	Image const *image=args->fwd?args->src:args->dst;
	memset(args->pixels, 0, args->bufsize);

	if(image->nch<3)
	{
		double csizes[OCH_COUNT*PRED_COUNT]={0};
		char predsel[OCH_COUNT]={0};
		memset(args->hist, 0, args->histsize);
		for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
		{
			ALIGN(16) int *rows[]=
			{
				args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*OCH_COUNT,
				args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*OCH_COUNT,
				args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*OCH_COUNT,
				args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*OCH_COUNT,
			};
			short input[ICH_COUNT]={0};
			int preds[PRED_COUNT]={0};
			int wg_errors[WG_NPREDS]={0}, wg_preds[WG_NPREDS]={0};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				int
					*NNN	=rows[3]+0*OCH_COUNT,
					*NNW	=rows[2]-1*OCH_COUNT,
					*NN	=rows[2]+0*OCH_COUNT,
					*NNE	=rows[2]+1*OCH_COUNT,
					*NNEE	=rows[2]+2*OCH_COUNT,
					*NW	=rows[1]-1*OCH_COUNT,
					*N	=rows[1]+0*OCH_COUNT,
					*NE	=rows[1]+1*OCH_COUNT,
					*NEE	=rows[1]+2*OCH_COUNT,
					*NEEE	=rows[1]+3*OCH_COUNT,
					*WWW	=rows[0]-3*OCH_COUNT,
					*WW	=rows[0]-2*OCH_COUNT,
					*W	=rows[0]-1*OCH_COUNT,
					*curr	=rows[0]+0*OCH_COUNT;
				if(ky<=args->y1+2)
				{
					if(ky<=args->y1+1)
					{
						if(ky==args->y1)
							NEEE=NE=NW=N=W;
						NN=N;
						NNE=NE;
					}
					NNN=NN;
				}
				if(kx<=2)
				{
					if(kx<=1)
					{
						if(!kx)
							NW=W=N;
						WW=W;
					}
					WWW=WW;
				}
				if(kx>=image->iw-3)
				{
					if(kx>=image->iw-1)
					{
						NNE=NN;
						NE=N;
					}
					NEEE=NE;
				}
				//NONE
				input[ICH_R]=image->data[idx+0];//r
				input[ICH_G]=image->data[idx+1];//g
				input[ICH_B]=image->data[idx+2];//b

				//JPEG2000	r-=g; b-=g; g+=(r+b)>>2;
				input[ICH_J2V]=input[ICH_R]-input[ICH_G];
				input[ICH_J2U]=input[ICH_B]-input[ICH_G];
				input[ICH_J2Y]=input[ICH_G]+((input[ICH_J2V]+input[ICH_J2U])>>2);

				//RCT1		r-=g; g+=r>>1; b-=g; g+=b>>1;
			//	input[ICH_R1V]=input[ICH_R]-input[ICH_G];
				input[ICH_R1Y]=input[ICH_G]+(input[ICH_R1V]>>1);
				input[ICH_R1U]=input[ICH_B]-input[ICH_G];
				input[ICH_R1Y]+=input[ICH_R1U]>>1;

				//Pei09		b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
				input[ICH_P9U]=(87*input[ICH_R]+169*input[ICH_G]+128)>>8;
			//	input[ICH_P9V]=input[ICH_R]-input[ICH_G];
				input[ICH_P9Y]=(86*input[ICH_P9V]+29*input[ICH_P9U]+128)>>8;
				for(int kc=0;kc<OCH_COUNT;++kc)
				{
					int offset=curr[dependencies[kc][1]];

					preds[PRED_W]=W[kc];
					MEDIAN3_32(preds[PRED_cgrad], N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
					preds[PRED_wgrad]=wg_predict(
						args->wg_weights+WG_NPREDS*kc,
						NNN[kc],
						NN[kc], NNE[kc],
						NW[kc], N[kc], NE[kc], NEE[kc], NEEE[kc],
						WWW[kc], WW[kc], W[kc],
						wg_errors, wg_preds
					);

					if(offset)
					{
						for(int kp=0;kp<PRED_COUNT;++kp)
						{
							preds[kp]+=offset;
							CLAMP2_32(preds[kp], preds[kp], -halfs[kc], halfs[kc]-1);
						}
					}
					curr[dependencies[kc][0]];
				}

				rows[0]+=OCH_COUNT;
				rows[1]+=OCH_COUNT;
				rows[2]+=OCH_COUNT;
				rows[3]+=OCH_COUNT;
			}
		}
	}
	else
	{
	}
}
int f26_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
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
		arg->bufsize=sizeof(int[4*4*2])*(image->iw+16LL);//4 padded rows * 4 channels max * {pixels, errors}
		arg->pixels=(int*)_mm_malloc(arg->bufsize, sizeof(__m128i));
		arg->histsize=sizeof(int[MAXPREDS])<<image->depth;
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
				block_thread(args+k);
#else
			void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
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
				block_thread(args+k);
#else
			{
				void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
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