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

//	#define USE_NONLINEARITY	//bad with backpropagation
	#define USE_BACKPROP		//good, keep enabled


#define BLOCKSIZE 256	//65536
#define CTXBITS 8
#define CTXSIZE (1<<CTXBITS)
//#define NCTRS 2
#define NPREDS 10
#define PREDLIST\
	PRED(1, pred)\
	PRED(1, (N[kc+4]+W[kc+4]+NW[kc+4]+NE[kc+4])/4)\
	PRED(1, 2*N[kc+0]-NN[kc+0])\
	PRED(1, 2*W[kc+0]-WW[kc+0])\
	PRED(1, N[kc+0]+W[kc+0])\
	PRED(1, N[kc+4]+NN[kc+4])\
	PRED(1, W[kc+4]+WW[kc+4])\
	PRED(1, N[kc+4]+W[kc+4])\
	PRED(1, kc?curr[0]:W[0])\
	PRED(1, kc?curr[4]:W[4])
//	PRED(1, 3*(N[kc]-NN[kc])+NNN[kc])\
//	PRED(1, N[kc]+NW[kc]-NNW[kc])\
//	PRED(1, NN[kc])\
//	PRED(1, N[kc]+NE[kc]-NNE[kc])\
//	PRED(1, 2*NE[kc]-NNEE[kc])\
//	PRED(1, NW[kc])\
//	PRED(1, N[kc])\
//	PRED(1, NE[kc])\
//	PRED(1, (W[kc]+NEE[kc])/2)\
//	PRED(1, (3*W[kc]+NEEE[kc])/4)\
//	PRED(1, 3*(W[kc]-WW[kc])+WWW[kc])\
//	PRED(1, WW[kc])\
//	PRED(1, W[kc])
//	PRED(1, NNN	[kc+4])\
//	PRED(1, NNW	[kc+4])\
//	PRED(1, NN	[kc+4])\
//	PRED(1, NNE	[kc+4])\
//	PRED(1, NNEE	[kc+4])\
//	PRED(1, NW	[kc+4])\
//	PRED(1, N	[kc+4])\
//	PRED(1, NE	[kc+4])\
//	PRED(1, NEE	[kc+4])\
//	PRED(1, NEEE	[kc+4])\
//	PRED(1, WWW	[kc+4])\
//	PRED(1, WW	[kc+4])\
//	PRED(1, W	[kc+4])
//	PRED(1, N[kc+0])
//	PRED(1, N[kc+4])
//	PRED(1, W[kc+0])
//	PRED(1, W[kc+4])
typedef unsigned short Counter_t;


#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
#include"ac.h"
typedef struct _ThreadArgs
{
	const Image *src;
	Image *dst;
	char depths[4];
	int maxdepth;
	int fwd, loud, y1, y2;
	int *pixels, bufsize;
	DList list;
	const unsigned char *decstart, *decend;

	Counter_t *stats;
	int statssize;
	int *mixer, mixersize;
#ifdef USE_NONLINEARITY
	const int *stretch, *squash;
#endif
} ThreadArgs;

static void print_i32(int *data, int count, int stride, int mag)
{
	for(int k=0;k<count;k+=stride)
	{
		int val=data[k];
		int nstars=val*64/mag;
		printf("%8d %8d ", k, val);
		for(int k2=0;k2<nstars;++k2)
			printf("*");
		printf("\n");
	}
}
#ifdef USE_NONLINEARITY
//stats -> stretch -> mix -> squash -> AC	probability has nonlinear resolution,  each table has one extra element
//prob->mix:	stretch(x) = ln(x/(1-x))
//mix->prob:	squash(x) = 1/(1+exp(-x))
static void stretchsquash_inittables(int *stretch, int *squash, int probbits, int mixbits)
{
	int plevels=1<<probbits, mlevels=1<<mixbits, mhalf=mlevels>>1;
	memset(stretch, 0, plevels*sizeof(int));
	memset(squash, 0, mlevels*sizeof(int));
	for(int km=-mhalf;km<=mhalf;++km)//prob = squash(mix)
	{
		int p;
		double x, rx;
		x=1/(1+exp(-7.*km/mhalf));//FIXME tune the coefficient
		rx=x*plevels;
		if(x<0.5)
			rx=ceil(rx);
		else if(x>0.5)
			rx=floor(rx);
#ifdef _MSC_VER
		p=_cvt_dtoi_fast(rx);
#else
		p=(int)rx;
#endif
		squash[km+mhalf]=p;
	}
	int c=0;
	for(int km=0;km<=mlevels;++km)//mix = stretch(prob)		just invert squash
	{
		int prob=squash[km];
		for(int kp=c;kp<prob;++kp)
			stretch[kp]=km;
		c=prob;
	}
	for(int kp=c;kp<=mlevels;++kp)
		stretch[kp]=plevels-1;

#if 0
	printf("mix = stretch(prob):\n");
	print_i32(stretch, plevels, 256, mlevels);
	printf("\n");

	printf("prob = squash(mix):\n");
	print_i32(squash, mlevels, 256, plevels);
	printf("\n");

	printf("squash(stretch(p)):\n");
	for(int kp=0;kp<plevels;kp+=256)
	{
		int kp2=squash[stretch[kp]];
		int nstars=kp2*64/plevels;
		printf("%8d %8d ", kp, kp2);
		for(int k2=0;k2<nstars;++k2)
			printf("*");
		printf("\n");
	}
	printf("\n");
#endif
}
#if 0
static int squash(int d)//return p = 1/(1+exp(-d)), d scaled by 8 bits, p scaled by 12 bits
{
	static int initialized=0;
	static short table[4096];
	if(!initialized)
	{
		static const int t[33]=
		{
			1, 2, 3, 6, 10, 16, 27, 45, 73, 120,
			194, 310, 488, 747, 1101, 1546, 2047, 2549, 2994, 3348,
			3607, 3785, 3901, 3975, 4022, 4050, 4068, 4079, 4085, 4089,
			4092, 4093, 4094,
		};
		for(int i=-2048;i<2048;++i)
		{
			int d=(i>>7)+16, w=i&127;
			//table[i+2048]=((t[d]<<7)+(t[d+1]-t[d])*w)>>7;
			table[i+2048]=(t[d]*(128-w)+t[d+1]*w+64) >> 7;
		}
		initialized=1;
	}
	d+=2048;
	if(d<0)
		return 0;
	if(d>4095)
		return 4095;
	return table[d];
}
//Inverse of squash.
// d = stretch(p) = ln(p/(1-p)), d scaled by 8 bits, p by 12 bits.
// d has range -2047 to 2047 representing -8 to 8.  
// p has range 0 to 4095 representing 0 to 1.
static int stretch(int p)
{
	static int initialized=0;
	static short table[4096];
	if(!initialized)
	{
		int pi=0;
		for(int x=-2047;x<=2047;++x)//invert squash()
		{
			int i=squash(x);
			for(int j=pi;j<=i;++j)
				table[j]=x;
			pi=i+1;
		}
		table[4095]=2047;
		initialized=1;
	}
	p=CLAMP(0, p, 4095);//
	return table[p];
}
#endif
#endif
static void block_thread(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	int fwd=args->src!=0;
	const Image *image=fwd?args->src:args->dst;
	ArithmeticCoder ec;
	//int mixer[4*NPREDS]={0};

	*args->mixer=0x8000;
	memfill(args->mixer+1, args->mixer, args->mixersize-sizeof(int), sizeof(int));
	if(fwd)
	{
		dlist_init(&args->list, 1, 1024, 0);
		ac_enc_init(&ec, &args->list);
	}
	else
		ac_dec_init(&ec, args->decstart, args->decend);
	*args->stats=0x8000;
	memfill(args->stats+1, args->stats, args->statssize-sizeof(Counter_t), sizeof(Counter_t));
#if 0
	for(int kc=0;kc<image->nch;++kc)
	{
		int nlevels=1<<args->depths[kc], half=nlevels>>1;
		Counter_t *curr_stats=args->stats+((ptrdiff_t)kc*CTXSIZE<<args->maxdepth);
		*curr_stats=0x8000;//unused
		for(int ks=0;ks<nlevels;++ks)
		{
			int idx=1;
			int inc=half-ks;
			inc=half-abs(inc);
			//inc*=inc;
			//inc=(int)((long long)nlevels*nlevels>>2)-inc;
			inc>>=args->depths[kc]-8;
			//int den=abs(ks-half)+16;
			for(int kb=args->depths[kc]-1;kb>=0;--kb)
			{
				int bit=ks>>kb&1;
				int p0=curr_stats[idx];
				p0+=bit?-inc:inc;
				//p0+=((!bit<<16)-p0)/den>>3;
				CLAMP2_32(curr_stats[idx], p0, 1, 0xFFFF);
			}
		}
		//memfill(curr_stats+(1LL<<args->maxdepth), curr_stats, sizeof(Counter_t[CTXSIZE*NPREDS*NCTRS-1LL])<<args->maxdepth, sizeof(Counter_t)<<args->maxdepth);
		memfill(curr_stats+(1LL<<args->maxdepth), curr_stats, sizeof(Counter_t[CTXSIZE*NPREDS-1LL])<<args->maxdepth, sizeof(Counter_t)<<args->maxdepth);
	}
	//for(int k=0;k<args->statssize/sizeof(Counter_t);++k)
	//{
	//	if(!args->stats[k])
	//		LOG_ERROR("");
	//}
#endif
	//memset(args->stats, 0, args->statssize);
	memset(args->pixels, 0, args->bufsize);
	for(int ky=args->y1, idx=image->nch*image->iw*args->y1;ky<args->y2;++ky)
	{
		ALIGN(16) int *rows[]=
		{
			args->pixels+((image->iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			args->pixels+((image->iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			args->pixels+((image->iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			args->pixels+((image->iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		short yuv[4]={0};
		int pred=0;
		//Counter_t *curr_stats[NPREDS*NCTRS]={0};
		Counter_t *curr_stats[NPREDS]={0};
		int idx3[NPREDS]={0};
		int *curr_mixer=0;
		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			int
				*NNN	=rows[3]+0*4*2,
				*NNW	=rows[2]-1*4*2,
				*NN	=rows[2]+0*4*2,
				*NNE	=rows[2]+1*4*2,
				*NNEE	=rows[2]+2*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*NEE	=rows[1]+2*4*2,
				*NEEE	=rows[1]+3*4*2,
				*WWW	=rows[0]-3*4*2,
				*WW	=rows[0]-2*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
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
				if(kx>=image->iw-2)
				{
					if(kx>=image->iw-1)
					{
						NNE=NN;
						NE=N;
					}
					NEE=NE;
				}
				NEEE=NEE;
			}
			//if(ky==1&&kx==20)//
			//if(ky==1&&kx==6)//
			//	printf("");

			if(fwd)
			{
				memcpy(yuv, image->data+idx, sizeof(short)*image->nch);
				if(image->nch>=3)
				{
					short temp;

					//Pei09 RCT	b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
					yuv[2]-=(87*yuv[0]+169*yuv[1]+128)>>8;
					yuv[0]-=yuv[1];
					yuv[1]+=(86*yuv[0]+29*yuv[2]+128)>>8;

					ROTATE3(yuv[0], yuv[1], yuv[2], temp);
				}
			}
			for(int kc=0;kc<image->nch;++kc)
			{
				static const int fixed_weights[NPREDS]=
				{
#define PRED(WEIGHT, EXPR) WEIGHT,
					PREDLIST
#undef  PRED
				};
				//int
				//	vx=(abs(W[kc]-WW[kc])+abs(N[kc]-NW[kc])+abs(NE[kc]-N[kc])+abs(W[kc+1]+WW[kc+1]))>>2,
				//	vy=(abs(W[kc]-NW[kc])+abs(N[kc]-NN[kc])+abs(NE[kc]-NNE[kc])+abs(N[kc+1]+NN[kc+1]))>>2;
				MEDIAN3_32(pred, N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
				int idx2[NPREDS]=
				{
#define PRED(WEIGHT, EXPR) EXPR,
					PREDLIST
#undef  PRED
#if 0
					//0,
					//0,
					curr[0],
					curr[4],
					(N[kc+0]+W[kc+0])/2,
					pred,
					pred,
					//N[kc+0]+W[kc+0]-NW[kc+0],
					//N[kc+4]+W[kc+4]-NW[kc+4],
					//N[kc+0],
					//N[kc+4],
					//W[kc+0],
					//W[kc+4],
					(N[kc+4]+W[kc+4]+NW[kc+4]+NE[kc+4])/4,
					//(3*(N[kc+4]+W[kc+4])+NW[kc+4]+NE[kc+4])/8,
					(N[kc+4]+W[kc+4])/2,
					//MINVAR(abs(N[kc+4]), abs(W[kc+4])),
					//vx,
					//vy,
					//N[kc+4]*3,
					//W[kc+4]*3,
#endif
				};
				int val=0;
				int treeidx=1;

				//if(ky==1&&kx==5&&kc==1)//
				//	printf("");

				for(int kp=0;kp<NPREDS;++kp)
				{
					int t=idx2[kp];
					t=abs(t);
					t>>=args->depths[kc]-CTXBITS;
					//t=FLOOR_LOG2_P1(t);
					if(t>CTXSIZE-1)
						t=CTXSIZE-1;
					//if(t>args->maxdepth-1)
					//	t=args->maxdepth-1;
					idx3[kp]=t;
					//t+=CTXSIZE*NCTRS*(NPREDS*kc+kp);
					t+=CTXSIZE*(NPREDS*kc+kp);
					//t+=NPREDS*kc;
					//t+=args->maxdepth*kc;
					t<<=args->maxdepth;
					curr_stats[kp]=args->stats+t;
				}
				//curr_mixer=args->mixer+NCTRS*NPREDS*args->maxdepth*kc;
				curr_mixer=args->mixer+NPREDS*args->maxdepth*kc;

				if(fwd)
				{
					val=yuv[kc];
					curr[kc+0]=val;
					val-=pred;
					val<<=32-args->depths[kc];
					val>>=32-args->depths[kc];
					curr[kc+4]=abs(val);
				}

				for(int kb=args->depths[kc]-1;kb>=0;--kb)
				{
					int bit;
#ifdef USE_BACKPROP
					long long psum=0;
					int wsum=0;
					int p0_unclamped;
#endif
					long long p0=0;
					//if(ky==1&&kx==6&&kc==1&&kb==8)//
					//	printf("");
					for(int kp=0;kp<NPREDS;++kp)
					{
						int weight=curr_mixer[kp]*fixed_weights[kp];
						int temp=curr_stats[kp][treeidx];
#ifdef USE_NONLINEARITY
						temp=args->stretch[temp];
#endif
						p0+=(long long)weight*temp;
#ifdef USE_BACKPROP
						wsum+=weight;
#endif
					}
#ifdef USE_BACKPROP
					psum=p0;
					p0/=wsum;
					p0_unclamped=(int)p0;
#else
					p0/=NPREDS<<15;
#endif
					CLAMP2_32(p0, (int)p0, 1, 0xFFFF);
#ifdef USE_NONLINEARITY
					p0=args->squash[p0];
#endif
					if(fwd)
					{
						bit=val>>kb&1;
						ac_enc_bin(&ec, (int)p0, bit);
					}
					else
					{
						bit=ac_dec_bin(&ec, (int)p0);
						val|=bit<<kb;
					}
#ifdef USE_BACKPROP
					int prob=(int)(bit?0x10000-p0:p0);
					if((unsigned)(p0_unclamped-1)<0xFFFF-1)
					{
						long long dL_dp0=-(1LL<<32)/prob;
						dL_dp0^=-bit;
						dL_dp0+=bit;
						for(int kp=0;kp<NPREDS;++kp)
						{
							int p2=curr_stats[kp][treeidx];
							int diff=p2-(int)p0;
							int update=(int)(dL_dp0*diff/wsum)*2500>>16;//0.037 = 2424	0.07 = 4588
							update=curr_mixer[kp]-update;
							CLAMP2_32(curr_mixer[kp], (int)update, 1, 0xFFFF);

							p2+=((!bit<<16)-p2)>>5;
							CLAMP2_32(curr_stats[kp][treeidx], p2, 1, 0xFFFF);
						}
					}
#else
					int proberror=(!bit<<16)-(int)p0;
					for(int kp=0;kp<NPREDS;++kp)
					{
						int prob=curr_stats[kp][treeidx];
						curr_mixer[kp]+=(int)((long long)proberror*prob>>24);//FIXME tune this
						prob+=((!bit<<16)-prob)>>5;//FIXME tune this
						CLAMP2_32(curr_stats[kp][treeidx], prob, 1, 0xFFFF);
					}
#endif
					treeidx=treeidx<<1|bit;
					curr_mixer+=NPREDS;
				}
				if(!fwd)
				{
					val<<=32-args->depths[kc];//extend sign
					val>>=32-args->depths[kc];
					curr[kc+4]=abs(val);
					val+=pred;
					val<<=32-args->depths[kc];
					val>>=32-args->depths[kc];
					curr[kc+0]=val;
					yuv[kc]=val;
				}
			}
			if(!fwd)
			{
				if(image->nch>=3)
				{
					short temp;

					ROTATE3(yuv[2], yuv[1], yuv[0], temp);

					yuv[1]-=(86*yuv[0]+29*yuv[2]+128)>>8;
					yuv[0]+=yuv[1];
					yuv[2]+=(87*yuv[0]+169*yuv[1]+128)>>8;
				}
				memcpy(args->dst->data+idx, yuv, sizeof(short)*image->nch);
			}
#ifdef ENABLE_GUIDE
			if(memcmp(image->data+idx, guide->data+idx, sizeof(short)*image->nch))
			{
				short temp;
				short orig[4]={0};
				memcpy(orig, guide->data+idx, image->nch*sizeof(short));

				//Pei09 RCT	b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
				orig[2]-=(87*orig[0]+169*orig[1]+128)>>8;
				orig[0]-=orig[1];
				orig[1]+=(86*orig[0]+29*orig[2]+128)>>8;
				ROTATE3(orig[0], orig[1], orig[2], temp);

				yuv[2]-=(87*yuv[0]+169*yuv[1]+128)>>8;
				yuv[0]-=yuv[1];
				yuv[1]+=(86*yuv[0]+29*yuv[2]+128)>>8;
				ROTATE3(yuv[0], yuv[1], yuv[2], temp);
				LOG_ERROR("Guide error XY %d %d", kx, ky);
				printf("");//
			}
#endif
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	if(fwd)
		ac_enc_flush(&ec);
}
int f25_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	int ncores, nblocks, nthreads, coffset;
	char depths[4]={0};
	int maxdepth;
	ptrdiff_t memusage, argssize;
	ThreadArgs *args;
#ifdef USE_NONLINEARITY
	int *stretch, *squash;
#endif
	ptrdiff_t start;
	
#ifdef ENABLE_GUIDE
	if(fwd)
		guide=image;
#endif
	memfill(depths, &image->depth, sizeof(char)*image->nch, sizeof(char));
	if(image->nch>=3)
	{
		depths[1]+=depths[1]<16;
		depths[2]+=depths[2]<16;
	}
	maxdepth=depths[0];
	for(int kc=1;kc<image->nch;++kc)
	{
		UPDATE_MAX(maxdepth, depths[kc]);
	}

	ncores=query_cpu_cores();
	nblocks=(image->ih+BLOCKSIZE-1)/BLOCKSIZE;
	nthreads=MINVAR(nblocks, ncores);
	coffset=sizeof(int)*nblocks;

	memusage=0;
	argssize=nthreads*sizeof(ThreadArgs);
	args=(ThreadArgs*)malloc(argssize);
#ifdef USE_NONLINEARITY
	stretch=(int*)malloc(sizeof(int[0x10001]));
	squash=(int*)malloc(sizeof(int[0x10001]));
#endif
	if(!args
#ifdef USE_NONLINEARITY
		||!stretch||!squash
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memusage+=argssize;
	memset(args, 0, argssize);
#ifdef USE_NONLINEARITY
	stretchsquash_inittables(stretch, squash, 16, 16);
#endif
	if(fwd)
		start=array_append(data, 0, 1, coffset, 1, 0, 0);
	else//integrity check
	{
		start=coffset;
		for(int kt=0;kt<nblocks;++kt)
		{
			int size=0;
			memcpy(&size, cbuf+sizeof(int)*kt, sizeof(int));
			start+=size;
		}
		if(start!=clen)
			LOG_ERROR("Corrupt file");
		start=coffset;
	}
	for(int kt=0;kt<nthreads;++kt)
	{
		ThreadArgs *arg=args+kt;
		arg->src=src;
		arg->dst=dst;
		memcpy(arg->depths, depths, sizeof(char[4]));
		arg->maxdepth=maxdepth;
		arg->fwd=fwd;
#ifdef DISABLE_MT
		arg->loud=loud;
#else
		arg->loud=0;
#endif
		arg->bufsize=sizeof(int[4*4*2])*(image->iw+16LL);//4 padded rows * NCHPOOL * {pixels, hires-errors, 4 pred hires-errors}
		arg->pixels=(int*)_mm_malloc(arg->bufsize, sizeof(__m128i));

		//arg->statssize=(int)sizeof(Counter_t[CTXSIZE*NPREDS*NCTRS])*image->nch<<maxdepth;
		arg->statssize=(int)sizeof(Counter_t[CTXSIZE*NPREDS])*image->nch<<maxdepth;
		arg->stats=(Counter_t*)malloc(arg->statssize);
		
		//arg->mixersize=sizeof(int[NPREDS*NCTRS])*image->nch*maxdepth;
		arg->mixersize=sizeof(int[NPREDS])*image->nch*maxdepth;
		arg->mixer=(int*)malloc(arg->mixersize);
		if(!arg->pixels||!arg->stats||!arg->mixer)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memusage+=arg->bufsize;
		memusage+=arg->statssize;
		
#ifdef USE_NONLINEARITY
		arg->stretch=stretch;
		arg->squash=squash;
#endif
	}
	for(int kb=0;kb<nblocks;kb+=nthreads)
	{
		int nthreads2=MINVAR(kb+nthreads, nblocks)-kb;
		for(int kt=0;kt<nthreads2;++kt)
		{
			ThreadArgs *arg=args+kt;
			arg->y1=BLOCKSIZE*(kb+kt);
			arg->y2=MINVAR(arg->y1+BLOCKSIZE, image->ih);
			if(!fwd)
			{
				int size=0;
				memcpy(&size, cbuf+sizeof(int)*((ptrdiff_t)kb+kt), sizeof(int));
				arg->decstart=cbuf+start;
				start+=size;
				arg->decend=cbuf+start;
			}
		}
#ifdef DISABLE_MT
		for(int kt=0;kt<nthreads2;++kt)
			block_thread(args+kt);
#else
		void *ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
		mt_finish(ctx);
#endif
		if(fwd)
		{
			for(int kt=0;kt<nthreads2;++kt)
			{
				ThreadArgs *arg=args+kt;
				int
					blocksize=(image->iw*(arg->y2-arg->y1)*image->nch*image->depth+7)>>3,
					cbsize=(image->iw*(arg->y2-arg->y1)*image->depth+7)>>3;

				memcpy(data[0]->data+start+sizeof(int)*((ptrdiff_t)kb+kt), &arg->list.nobj, sizeof(int));
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
			printf(
				"%9td/%9td bytes  %10lf%% %10lf\n",
				csize,
				usize,
				100.*csize/usize,
				(double)usize/csize
			);
		}
		printf("%c %16.6lf sec  %16.6lf MB/s  Mem usage: ", 'D'+fwd, t0, usize/(t0*1024*1024));
		print_size((double)memusage, 8, 4, 0, 0);
		printf("\n");
	}
	for(int kt=0;kt<nthreads;++kt)
	{
		ThreadArgs *arg=args+kt;
		_mm_free(arg->pixels);
		free(arg->stats);
	}
	free(args);
	return 0;
}