#include"codec.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>//abs
//#include<immintrin.h>//included by "entropy.h"
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
#ifndef DISABLE_MT
//	#define ENABLE_MT
#endif

	#define EXPONENTIAL_WINDOW
//	#define PROFILE_SIZE

#define CODECNAME "C06"
#include"entropy.h"

#define BLOCKSIZE 512
#define MAXPRINTEDBLOCKS 200
#define NCTX 8
#define CTXBITS 8
//static unsigned short g_stats0[3][1<<8][2];
static unsigned short g_stats1[3][NCTX][1<<CTXBITS][1<<8][2];
//static int g_mixer[3][1<<CTXBITS][NCTX+1];
typedef struct _ThreadArgs
{
	const unsigned char *src;
	unsigned char *dst;
	int iw, ih;

	int fwd, test, loud, x1, x2, y1, y2;
	short pixels[(BLOCKSIZE+16)*4*4];//4 padded rows * 4 channels max

	BList list;
	const unsigned char *decstart, *decend;

	unsigned short stats0[3][1<<8][2];
	unsigned short stats1[3][NCTX][1<<CTXBITS][1<<8][2];
	int mixer[3][1<<CTXBITS][NCTX+1];

	//aux
	int blockidx;
	double bestsize;
#ifdef PROFILE_SIZE
	double csizes[24];
#endif
} ThreadArgs;
#define SCALE_BITS 16
static int squash(int x)//sigmoid(x) = 1/(1+exp(-x))		logit sum -> prob
{
#ifdef DISABLE_LOGMIX
	x>>=11;
	x+=1<<SCALE_BITS>>1;
	CLAMP2_32(x, x, 1, (1<<SCALE_BITS)-1);
#else
	static const int t[33]=//2^5 table elements, table amplitude 2^12
	{
		   1,    2,    3,    6,   10,   16,   27,   45,   73,  120,  194,
		 310,  488,  747, 1101, 1546, 2047, 2549, 2994, 3348, 3607, 3785,
		3901, 3975, 4022, 4050, 4068, 4079, 4085, 4089, 4092, 4093, 4094,
	};
	int w=x&((1<<(SCALE_BITS-5))-1);
	x=(x>>(SCALE_BITS-5))+16;
	if(x>31)
		return (1<<SCALE_BITS)-1;
	if(x<0)
		return 1;
	x=(t[x]*((1<<(SCALE_BITS-5))-w)+t[x+1]*w+64)>>(12-5);
#endif
	return x;
}
static int pred_cgrad(int near1, int near2, int far)
{
	int pred;
	MEDIAN3_32(pred, near1, near2, near1+near2-far);
	return pred;
}
static void block_thread(void *param)
{
	const int nch=3;
	ThreadArgs *args=(ThreadArgs*)param;
	AC4 ec;
	const unsigned char *image=args->fwd?args->src:args->dst;
	
	if(args->fwd)
	{
		blist_init(&args->list);
		ac4_enc_init(&ec, &args->list);
	}
	else
		ac4_dec_init(&ec, args->decstart, args->decend);
#ifdef PROFILE_SIZE
	memset(args->csizes, 0, sizeof(args->csizes));
#endif
#ifdef EXPONENTIAL_WINDOW
	FILLMEM((unsigned short*)args->stats0, 0x8000, sizeof(args->stats0), sizeof(short));
	if(args->blockidx)//load stats
		memcpy(args->stats1, g_stats1, sizeof(args->stats1));
	else
		FILLMEM((unsigned short*)args->stats1, 0x8000, sizeof(args->stats1), sizeof(short));
#else
	memset(args->stats0, 0, sizeof(args->stats0));
	memset(args->stats1, 0, sizeof(args->stats1));
#endif
	FILLMEM((int*)args->mixer, 0x4000, sizeof(args->mixer), sizeof(int));
	memset(args->pixels, 0, sizeof(args->pixels));
	for(int ky=args->y1;ky<args->y2;++ky)//codec loop
	{
		ALIGN(16) short *rows[]=
		{
			args->pixels+((BLOCKSIZE+16LL)*((ky-0LL)&3)+8LL)*4,
			args->pixels+((BLOCKSIZE+16LL)*((ky-1LL)&3)+8LL)*4,
			args->pixels+((BLOCKSIZE+16LL)*((ky-2LL)&3)+8LL)*4,
			args->pixels+((BLOCKSIZE+16LL)*((ky-3LL)&3)+8LL)*4,
		};
		int yuv[4]={0};
		unsigned short (*curr_ctx[NCTX+1])[2]={0};
		int probs[NCTX+1]={0};
		int *curr_mixer=0;
		int val=0;
		for(int kx=args->x1;kx<args->x2;++kx)
		{
			int idx=nch*(args->iw*ky+kx);
			short
				*NNN	=rows[3]+0*4,
				*NNWW	=rows[2]-2*4,
				*NNW	=rows[2]-1*4,
				*NN	=rows[2]+0*4,
				*NNE	=rows[2]+1*4,
				*NNEE	=rows[2]+2*4,
				*NNEEE	=rows[2]+3*4,
				*NWW	=rows[1]-2*4,
				*NW	=rows[1]-1*4,
				*N	=rows[1]+0*4,
				*NE	=rows[1]+1*4,
				*NEE	=rows[1]+2*4,
				*NEEE	=rows[1]+3*4,
				*WWWW	=rows[0]-4*4,
				*WWW	=rows[0]-3*4,
				*WW	=rows[0]-2*4,
				*W	=rows[0]-1*4,
				*curr	=rows[0]+0*4;
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
				yuv[0]=args->src[idx+1];
				yuv[1]=args->src[idx+2];
				yuv[2]=args->src[idx+0];
			}
			for(int kc=0;kc<nch;++kc)
			{
				//if(idx==2319)//
				//if(idx==7695)//
				//	printf("");
#define GETCTX(X) ((X)>>(8-CTXBITS)&((1<<CTXBITS)-1))
				switch(kc)
				{
				case 0:
					{
						int x;
						int gx=abs(W[0]-WW[0])+abs(N[0]-NW[0])+abs(NE[0]-N[0])+1;
						int gy=abs(N[0]-NN[0])+abs(W[0]-NW[0])+abs(NE[0]-NNE[0])+1;
						curr_mixer=args->mixer[0][GETCTX((N[0]+W[0])>>1)];
						//curr_mixer=args->mixer[0][GETCTX((gx*N[0]+gy*W[0])/(gx+gy))];
						curr_ctx[0]=args->stats1[0][0][GETCTX(pred_cgrad(N[0], W[0], NW[0]))];
						curr_ctx[1]=args->stats1[0][1][GETCTX(pred_cgrad(N[0], NE[0], NNE[0]))];
						x=(6*W[0]-4*WW[0]+WWW[0])/3; CLAMP2_32(x, x, 0, 255); curr_ctx[2]=args->stats1[0][2][GETCTX(x)];
						x=(6*N[0]-4*NN[0]+NNN[0])/3; CLAMP2_32(x, x, 0, 255); curr_ctx[3]=args->stats1[0][3][GETCTX(x)];
						curr_ctx[4]=args->stats1[0][4][GETCTX(N[0])];
						curr_ctx[5]=args->stats1[0][5][GETCTX(W[0])];
						curr_ctx[6]=args->stats1[0][6][GETCTX(pred_cgrad(W[0], NE[0], N[0]))];
						curr_ctx[7]=args->stats1[0][7][GETCTX((gx*N[0]+gy*W[0])/(gx+gy))];
					//	curr_ctx[7]=args->stats1[0][7][GETCTX((N[0]+W[0]+NW[0]+NE[0])>>2)];
						curr_ctx[8]=args->stats0[0];
						//curr_ctx[3]=args->stats1[0][2][GETCTX(pred_cgrad(N[0], NW[0], NNW[0]))];
						//curr_ctx[3]=args->stats1[0][2][GETCTX(pred_cgrad(W[0], WW[0], WWW[0]))];
						//curr_ctx[4]=args->stats1[0][3][GETCTX(pred_cgrad(N[0], NN[0], NNN[0]))];
						//curr_ctx[1]=args->stats1[0][0][GETCTX(N[0]+W[0]-NW[0])];
						//curr_ctx[2]=args->stats1[0][1][GETCTX(N[0]+NE[0]-NNE[0])];
						//curr_ctx[3]=args->stats1[0][2][GETCTX(3*(W[0]-WW[0])+WWW[0])];
						//curr_ctx[4]=args->stats1[0][3][GETCTX(3*(N[0]-NN[0])+NNN[0])];
					}
					break;
				case 1:
					{
						int gx=abs(W[1]-WW[1]-(W[0]-WW[0]))+abs(N[1]-NW[1]-(N[0]-NW[0]))+abs(NE[1]-N[1]-(NE[0]-N[0]))+1;
						int gy=abs(N[1]-NN[1]-(N[0]-NN[0]))+abs(W[1]-NW[1]-(W[0]-NW[0]))+abs(NE[1]-NNE[1]-(NE[0]-NNE[0]))+1;
						int x;
						curr_mixer=args->mixer[1][GETCTX(curr[0])];
						//x=((N[1]-N[0]+W[1]-W[0])>>1)+curr[0]; CLAMP2_32(x, x, 0, 255); curr_mixer=args->mixer[1][GETCTX(x)];
						x=pred_cgrad(N[1]-N[0], W [1]-W [0], NW [1]-NW [0])+curr[0];	CLAMP2_32(x, x, 0, 255); curr_ctx[0]=args->stats1[1][0][GETCTX(x)];
						x=pred_cgrad(N[1]-N[0], NE[1]-NE[0], NNE[1]-NNE[0])+curr[0];	CLAMP2_32(x, x, 0, 255); curr_ctx[1]=args->stats1[1][1][GETCTX(x)];
						x=pred_cgrad(W[1]-W[0], WW[1]-WW[0], WWW[1]-WWW[0])+curr[0];	CLAMP2_32(x, x, 0, 255); curr_ctx[2]=args->stats1[1][2][GETCTX(x)];
						x=pred_cgrad(N[1]-N[0], NN[1]-NN[0], NNN[1]-NNN[0])+curr[0];	CLAMP2_32(x, x, 0, 255); curr_ctx[3]=args->stats1[1][3][GETCTX(x)];
						x=pred_cgrad(N[1], curr[0], N[0]);				CLAMP2_32(x, x, 0, 255); curr_ctx[4]=args->stats1[1][4][GETCTX(x)];
						x=pred_cgrad(W[1], curr[0], W[0]);				CLAMP2_32(x, x, 0, 255); curr_ctx[5]=args->stats1[1][5][GETCTX(x)];
						x=pred_cgrad(W[1]-W[0], NE[1]-NE[0], N[1]-N[0])+curr[0];	CLAMP2_32(x, x, 0, 255); curr_ctx[6]=args->stats1[1][6][GETCTX(x)];
						x=(gx*(N[1]-N[0])+gy*(W[1]-W[0]))/(gx+gy)+curr[0];		CLAMP2_32(x, x, 0, 255); curr_ctx[7]=args->stats1[1][7][GETCTX(x)];
					//	x=((N[1]-N[0]+W[1]-W[0]+NW[1]-NW[0]+NE[1]-NE[0])>>2)+curr[0];	CLAMP2_32(x, x, 0, 255); curr_ctx[7]=args->stats1[1][7][GETCTX(x)];
						curr_ctx[8]=args->stats0[1];
						//x=(6*W[1]-W[0]-4*(WW[1]-WW[0])+WWW[1]-WWW[0])/3+curr[0]; CLAMP2_32(x, x, 0, 255); curr_ctx[3]=args->stats1[1][2][GETCTX(x)];
						//x=(6*N[1]-N[0]-4*(NN[1]-NN[0])+NNN[1]-NNN[0])/3+curr[0]; CLAMP2_32(x, x, 0, 255); curr_ctx[4]=args->stats1[1][3][GETCTX(x)];
						//curr_ctx[1]=args->stats1[1][0][GETCTX(N[1]-N[0]+W[1]-W[0]-NW[1]+NW[0]+curr[0])];
						//curr_ctx[2]=args->stats1[1][1][GETCTX(N[1]-N[0]+NE[1]-NE[0]-NNE[1]+NNE[0]+curr[0])];
						//curr_ctx[3]=args->stats1[1][2][GETCTX(2*(W[1]-W[0])-WW[1]+WW[0]+curr[0])];
						//curr_ctx[4]=args->stats1[1][3][GETCTX(2*(N[1]-N[0])-NN[1]+NN[0]+curr[0])];
					}
					break;
				case 2:
					{
						int gx=abs(W[2]-WW[2]-(W[0]-WW[0]))+abs(N[2]-NW[2]-(N[0]-NW[0]))+abs(NE[2]-N[2]-(NE[0]-N[0]))+1;
						int gy=abs(N[2]-NN[2]-(N[0]-NN[0]))+abs(W[2]-NW[2]-(W[0]-NW[0]))+abs(NE[2]-NNE[2]-(NE[0]-NNE[0]))+1;
						int x;
						curr_mixer=args->mixer[2][GETCTX(curr[0])];
						//x=((N[2]-N[0]+W[2]-W[0])>>1)+curr[0]; CLAMP2_32(x, x, 0, 255); curr_mixer=args->mixer[2][GETCTX(x)];
						x=pred_cgrad(N[2]-N[0], W [2]-W [0], NW [2]-NW [0])+curr[0];	CLAMP2_32(x, x, 0, 255); curr_ctx[0]=args->stats1[2][0][GETCTX(x)];
						x=pred_cgrad(N[2]-N[1], NE[2]-NE[1], NNE[2]-NNE[1])+curr[1];	CLAMP2_32(x, x, 0, 255); curr_ctx[1]=args->stats1[2][1][GETCTX(x)];
						x=pred_cgrad(W[2]-W[0], WW[2]-WW[0], WWW[2]-WWW[0])+curr[0];	CLAMP2_32(x, x, 0, 255); curr_ctx[2]=args->stats1[2][2][GETCTX(x)];
						x=pred_cgrad(N[2]-N[1], NN[2]-NN[1], NNN[2]-NNN[1])+curr[1];	CLAMP2_32(x, x, 0, 255); curr_ctx[3]=args->stats1[2][3][GETCTX(x)];
						x=pred_cgrad(N[1], curr[0], N[0]);				CLAMP2_32(x, x, 0, 255); curr_ctx[4]=args->stats1[2][4][GETCTX(x)];
						x=pred_cgrad(W[1], curr[0], W[0]);				CLAMP2_32(x, x, 0, 255); curr_ctx[5]=args->stats1[2][5][GETCTX(x)];
						x=pred_cgrad(W[1]-W[0], NE[1]-NE[0], N[1]-N[0])+curr[0];	CLAMP2_32(x, x, 0, 255); curr_ctx[6]=args->stats1[2][6][GETCTX(x)];
						x=(gx*(N[2]-N[0])+gy*(W[2]-W[0]))/(gx+gy)+curr[0];		CLAMP2_32(x, x, 0, 255); curr_ctx[7]=args->stats1[1][7][GETCTX(x)];
					//	x=((N[1]-N[0]+W[1]-W[0]+NW[1]-NW[0]+NE[1]-NE[0])>>2)+curr[0];	CLAMP2_32(x, x, 0, 255); curr_ctx[7]=args->stats1[2][7][GETCTX(x)];
						curr_ctx[8]=args->stats0[2];
						//x=(6*W[2]-W[1]-4*(WW[2]-WW[1])+WWW[2]-WWW[1])/3+curr[1]; CLAMP2_32(x, x, 0, 255); curr_ctx[3]=args->stats1[2][2][GETCTX(x)];
						//x=(6*N[2]-N[0]-4*(NN[2]-NN[0])+NNN[2]-NNN[0])/3+curr[0]; CLAMP2_32(x, x, 0, 255); curr_ctx[4]=args->stats1[2][3][GETCTX(x)];
						//curr_ctx[1]=args->stats1[2][0][GETCTX(N[2]-N[0]+W[2]-W[0]-NW[2]+NW[0]+curr[0])];
						//curr_ctx[2]=args->stats1[2][1][GETCTX(N[2]-N[0]+NE[2]-NE[0]-NNE[2]+NNE[0]+curr[0])];
						//curr_ctx[3]=args->stats1[2][2][GETCTX(2*(W[2]-W[0])-WW[2]+WW[0]+curr[0])];
						//curr_ctx[4]=args->stats1[2][3][GETCTX(2*(N[2]-N[0])-NN[2]+NN[0]+curr[0])];
					}
					break;
				}
				val=0;
				if(args->fwd)
					val=yuv[kc];
				for(int kb=7, tidx=1;kb>=0;--kb)
				{
					int bit=0;
					long long p1=0;
					//if(ky==3&&kx==232&&kc==2&&kb==7)
					//	printf("");
					for(int k=0;k<NCTX+1;++k)
					{
						unsigned short *cell=curr_ctx[k][tidx];
						//if((size_t)(curr_ctx[k][tidx]-args->stats0)>sizeof(args->stats0)+sizeof(args->stats1))
						//	LOG_ERROR("");
						//if((size_t)(curr_mixer-args->mixer)>sizeof(args->mixer))
						//	LOG_ERROR("");
#ifdef EXPONENTIAL_WINDOW
						probs[k]=((cell[0]+cell[1]-0x10000+1)>>1);
#else
						probs[k]=(int)(((cell[1]+1)<<16)/((cell[0]+cell[1])+2))-0x8000;
#endif
						p1+=(long long)curr_mixer[k]*probs[k];
					//	p1+=(long long)probs[k];
					}
					p1=squash((int)((p1*110+(1<<23>>1))>>23));
					//p1=p1*26>>22;
					//p1/=NCTX+1;
					//p1=squash((int)(p1*9>>4));
					//p1+=0x8000;
					//CLAMP2_32(p1, (int)p1, 1, 0xFFFF);

					if(args->fwd)
					{
						bit=val>>kb&1;
						ac4_enc_bin(&ec, (int)p1, bit);
#ifdef PROFILE_SIZE
						args->csizes[kc<<3|kb]-=log2((double)(bit?p1:0x10000-p1)/0x10000);
#endif
					}
					else
					{
						bit=ac4_dec_bin(&ec, (int)p1);
						val|=bit<<kb;
					}
					//if(idx==7695)//
					//if(idx==2319)//
					//	printf("%d %d\n", bit, (int)p1);

					int proberror=(bit<<16)-(int)p1;
					int cellupdate0=proberror*14>>12;
				//	int cellupdate1=proberror*8>>9;
					for(int k=0;k<NCTX+1;++k)
					{
						unsigned short *cell=curr_ctx[k][tidx];
#ifdef EXPONENTIAL_WINDOW
						CLAMP2_32(cell[0], cell[0]+(int)((long long)curr_mixer[k]*cellupdate0>>11), 1, 0xFFFF);
						//CLAMP2_32(cell[1], cell[1]+cellupdate1, 1, 0xFFFF);
						//cell[0]+=((bit<<15)-cell[0])>>5;
						cell[1]+=((bit<<16)-cell[1]+(1<<6>>1))>>6;
#else
						++cell[bit];
						if(cell[bit]>=(k?128:4096))
						{
							cell[0]=(cell[0]+1)>>1;
							cell[1]=(cell[1]+1)>>1;
						}
#endif

						CLAMP2_32(curr_mixer[k], curr_mixer[k]+(int)((long long)probs[k]*proberror>>20), -0x8000, 0x8000);
					}
					tidx+=tidx+bit;
				}
				if(!args->fwd)
					yuv[kc]=val;
				curr[kc]=yuv[kc];
			}
			if(!args->fwd)
			{
				args->dst[idx+1]=yuv[0];
				args->dst[idx+2]=yuv[1];
				args->dst[idx+0]=yuv[2];
#ifdef ENABLE_GUIDE
				if(args->test&&memcmp(args->dst+idx, args->src+idx, sizeof(char)*nch))
				{
					unsigned char orig[4]={0};
					memcpy(orig, args->src+idx, nch*sizeof(char));
					LOG_ERROR("Guide error XY %d %d", kx, ky);
					printf("");//
				}
#endif
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	if(args->fwd)
		ac4_enc_flush(&ec);
	if(!args->blockidx)//save stats
		memcpy(g_stats1, args->stats1, sizeof(args->stats1));
}
int c06_codec(const char *srcfn, const char *dstfn)
{
	const int nch=3, depth=8;
	double t0;
	//double ttable;
	ArrayHandle src, dst;
	int headersize, printed;
	int iw, ih;
	const unsigned char *image, *imageend;
	unsigned char *image2;
	CodecID codec;
	int ncores;
	int xblocks, yblocks, nblocks, nthreads, coffset;
	ptrdiff_t start, memusage, argssize;
	ThreadArgs *args;
	int test, fwd;
	//unsigned short *divtable1, *divtable2, *divtable3;
	//unsigned *statetable;
	double esize;
	int usize;
#ifdef PROFILE_SIZE
	double csizes[24]={0};
#endif
	
	t0=time_sec();
	src=load_file(srcfn, 1, 3, 1);
	headersize=header_read(src->data, (int)src->count, &iw, &ih, &codec);
	image=src->data+headersize;
	imageend=src->data+src->count;
	if(codec==CODEC_INVALID||codec==CODEC_PGM)
	{
		LOG_ERROR("Unsupported codec %d.\n", codec);
		array_free(&src);
		return 1;
	}
	else if(codec==CODEC_C01&&!dstfn)
	{
		LOG_ERROR(
			"Test mode expects PPM source.\n"
			"Decode mode expects destination filename."
		);
		return 1;
	}
	test=!dstfn;
	fwd=codec==CODEC_PPM;
	
	if(test)
		printf("%s \"%s\"  WH %d*%d\n", CODECNAME, srcfn, iw, ih);
	usize=iw*ih*3;
	ncores=query_cpu_cores();
	xblocks=(iw+BLOCKSIZE-1)/BLOCKSIZE;
	yblocks=(ih+BLOCKSIZE-1)/BLOCKSIZE;
	nblocks=xblocks*yblocks, nthreads=MINVAR(nblocks, ncores);
	coffset=(int)sizeof(int)*nblocks;
	start=0;
	memusage=0;
	argssize=nthreads*sizeof(ThreadArgs);
	args=(ThreadArgs*)malloc(argssize);
	//divtable1=(unsigned short*)malloc(sizeof(short[1<<(A2_CTRBITS1<<1)]));
	//divtable2=(unsigned short*)malloc(sizeof(short[1<<(A2_CTRBITS2<<1)]));
	//divtable3=(unsigned short*)malloc(sizeof(short[1<<(A2_CTRBITS3<<1)]));
	//statetable=(unsigned*)malloc(sizeof(int[256]));
	if(!args)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
#if 0
	ttable=time_sec();
	for(int k=0;k<1<<(A2_CTRBITS1<<1);++k)
	{
		long long mag=8;
		int n0=k&((1<<A2_CTRBITS1)-1), n1=k>>A2_CTRBITS1&((1<<A2_CTRBITS1)-1);
		int total=(int)((n0+n1)*mag+2);
		int p1=(int)((((n1*mag+1)<<16)+total/2)/total);
		divtable1[k]=p1;
	}
	for(int k=0;k<1<<(A2_CTRBITS2<<1);++k)
	{
		long long mag=8;
		int n0=k&((1<<A2_CTRBITS2)-1), n1=k>>A2_CTRBITS2&((1<<A2_CTRBITS2)-1);
		int total=(int)((n0+n1)*mag+2);
		int p1=(int)((((n1*mag+1)<<16)+total-1)/total);
		divtable2[k]=p1;
	}
	for(int k=0;k<1<<(A2_CTRBITS3<<1);++k)
	{
		long long mag=8;
		int n0=k&((1<<A2_CTRBITS3)-1), n1=k>>A2_CTRBITS3&((1<<A2_CTRBITS3)-1);
		int total=(int)((n0+n1)*mag+2);
		int p1=(int)((((n1*mag+1)<<16)+0)/total);
		divtable3[k]=p1;
	}
	//for(int k=0;k<256;++k)
	//{
	//	const unsigned char *cell=g_state_table[k];
	//	int total=(cell[2]+cell[3])*8+2;
	//	int p1=(((cell[3]*8+1)<<16)+total/2)/total;
	//	statetable[k]=p1<<16|cell[1]<<8|cell[0];
	//}
	ttable=time_sec()-ttable;
#endif
	esize=0;
	memusage+=argssize;
	memset(args, 0, argssize);
	if(fwd)
	{
		dst=0;
		printed=snprintf(g_buf, G_BUF_SIZE-1, "C01\n%d %d\n", iw, ih);
		array_append(&dst, g_buf, 1, printed, 1, 0, 0);
		start=array_append(&dst, 0, 1, coffset, 1, 0, 0);
		
		image2=0;
	}
	else//integrity check
	{
		dst=0;
		printed=snprintf(g_buf, G_BUF_SIZE-1, "P6\n%d %d\n255\n", iw, ih);
		array_append(&dst, g_buf, 1, printed, 1, 0, 0);
		array_append(&dst, 0, 1, usize, 1, 0, 0);

		//printed=0;
		start=coffset;
		for(int kt=0;kt<nblocks;++kt)
		{
			int size=0;
			memcpy(&size, image+sizeof(int)*kt, sizeof(int));
			start+=size;
		}
		if(image+start!=imageend)
			LOG_ERROR("Corrupt file");
		start=coffset;

		image2=(unsigned char*)malloc(usize);
		if(!image2)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
		memset(image2, 0, usize);
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadArgs *arg=args+k;
		arg->src=image;
		arg->dst=fwd?0:dst->data+printed;
		arg->iw=iw;
		arg->ih=ih;
		
		arg->fwd=fwd;
		arg->test=test;
#ifdef ENABLE_MT
		arg->loud=0;
#else
		arg->loud=test&&nblocks<MAXPRINTEDBLOCKS;
#endif
	}
	for(int k2=0;k2<=test;++k2)
	{
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
				arg->x1=BLOCKSIZE*kx;
				arg->y1=BLOCKSIZE*ky;
				arg->x2=MINVAR(arg->x1+BLOCKSIZE, iw);
				arg->y2=MINVAR(arg->y1+BLOCKSIZE, ih);
				if(!fwd)
				{
					int size=0;
					memcpy(&size, image+sizeof(int)*((ptrdiff_t)kt+kt2), sizeof(int));
					arg->decstart=image+start;
					start+=size;
					arg->decend=image+start;
				}
			}
#ifdef ENABLE_MT
			void *ctx;
			if(kt)
			{
				ctx=mt_exec(block_thread, args, sizeof(ThreadArgs), nthreads2);
				mt_finish(ctx);
			}
			else
			{
				block_thread(args);
				ctx=mt_exec(block_thread, args+1, sizeof(ThreadArgs), nthreads2-1);
				mt_finish(ctx);
			}
#else
			for(int k=0;k<nthreads2;++k)
				block_thread(args+k);
#endif
			if(fwd)
			{
				for(int kt2=0;kt2<nthreads2;++kt2)
				{
					ThreadArgs *arg=args+kt2;
					if(test)
					{
						int blocksize=((arg->x2-arg->x1)*(arg->y2-arg->y1)*nch*depth+7)>>3;
						int kx, ky;

						kx=kt+kt2;
						ky=kx/xblocks;
						kx%=xblocks;
						if(nblocks<MAXPRINTEDBLOCKS)
						{
							printf(
								"block %4d/%4d  XY %3d %3d  %4d*%4d:  %8d->%16lf->%8zd bytes (%+10.2lf)  %10.6lf%%  CR %10lf\n",
								kt+kt2+1, nblocks,
								kx, ky,
								arg->y2-arg->y1,
								arg->x2-arg->x1,
								blocksize,
								arg->bestsize,
								arg->list.nbytes,
								arg->list.nbytes-arg->bestsize,
								100.*arg->list.nbytes/blocksize,
								(double)blocksize/arg->list.nbytes
							);
						}
						esize+=arg->bestsize;
#ifdef ABAC_PROFILESIZE
						for(int k=0;k<ABAC_TOKEN_BITS*3;++k)
							abac_csizes[k]+=arg->abac_csizes[k];
#endif
					}
					memcpy(dst->data+start+sizeof(int)*((ptrdiff_t)kt+kt2), &arg->list.nbytes, sizeof(int));
					blist_appendtoarray(&arg->list, &dst);
					blist_clear(&arg->list);
#ifdef PROFILE_SIZE
					for(int k=0;k<24;++k)
						csizes[k]+=arg->csizes[k];
#endif
				}
			}
		}
		if(test)
		{
			ptrdiff_t usize=((ptrdiff_t)iw*ih*nch*depth+7)>>3;
			ptrdiff_t csize=dst->count;
			t0=time_sec()-t0;
			if(fwd)
			{
#ifdef PROFILE_SIZE
				for(int k=0;k<24;++k)
				{
					double size=csizes[k]/8;
					printf("C%d B%d  %12.2lf\n", k>>3, k&7, size);
					if(!((k+1)&7))
						printf("\n");
				}
#endif
				printf("Best %15.2lf (%+13.2lf) bytes\n", esize, csize-esize);
				printf("%12td/%12td  %10.6lf%%  %10lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
				printf("Mem usage: ");
				print_size((double)memusage, 8, 4, 0, 0);
				printf("\n");
			}
			printf("%c %16.6lf sec  %16.6lf MB/s\n", 'D'+fwd, t0, usize/(t0*1024*1024));
			if(!fwd)
				compare_bufs_8(image2, src->data+headersize, iw, ih, nch, nch, CODECNAME, 0, 1);
		}
		if(!k2&&test)
		{
			int usize=iw*ih*3;
			fwd=0;
			image2=(unsigned char*)malloc(usize);
			if(!image2)
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			memset(image2, 0, usize);
			for(int kt=0;kt<nthreads;++kt)
			{
				ThreadArgs *arg=args+kt;
				arg->dst=image2;
				arg->fwd=0;
			}
			image=dst->data+printed;
			start=coffset;
		}
		t0=time_sec();
	}
	if(!test)
		save_file(dstfn, dst->data, dst->count, 1);
	if(image2)
		free(image2);
	free(args);
	array_free(&src);
	array_free(&dst);
	//if(test)
	//	pause();
	return 0;
}