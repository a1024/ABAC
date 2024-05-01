#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
static const char file[]=__FILE__;


//	#define PROFILER 1
	#define TRACK_MIXER
//	#define FULL_TREE
//	#define ENABLE_SSE
//	#define LEARNED_DENOMINATOR
//	#define LINEAR_QUANTIZER_BITS 8
//	#define DISABLE_RCT
//	#define PERMUTE_YUV


#ifdef PROFILER
#define CHECKPOINTLIST\
	CHECKPOINT(INIT)\
	CHECKPOINT(CTX)\
	CHECKPOINT(P0)\
	CHECKPOINT(EC)\
	CHECKPOINT(UPDATE)\
	CHECKPOINT(FINISH)
#endif
#include"ac.h"
#include"profiler.h"
#define MIXER_SHIFT 19
#define STATS_SHIFT 24
#define PADX 8
#define PADY 8
int f09_disable_ctx=-1;
#if 1
//21 predictors		first 4 are overwritten at chroma
#define PREDLIST\
	PRED(eNNN)\
	PRED(eWWW)\
	PRED(eWWW*3)\
	PRED(eWWWWW-7*eWWWW+21*eWWW+35*(-eWW+eW))\
	PRED(eN+eW)\
	PRED(4*(eNE+eNNNE)-6*eNNE+eNNNNE)\
	PRED(eWWW+eWW)\
	PRED(eNN)\
	PRED(eWW)\
	PRED(3*(eW-eWW)+eWWW)\
	PRED(pred*7)\
	PRED(pred*5)\
	PRED(pred*3)\
	PRED(pred)\
	PRED(rows[1][4+0*8])\
	PRED(rows[0][4-1*8])\
	PRED((rows[0][0-1*8]+rows[1][0+1*8])>>1)\
	PRED(eN*3)\
	PRED(eW*3)\
	PRED((eN+eW)*3)\
	PRED(eN-eW)\
	PRED(eNW)\
	PRED(eNW*7)\
	PRED(7*(eW-eWWWWWW)-21*(eWW-eWWWWW)+35*(eWWW-eWWWW)+eWWWWWWW)\
	PRED(eNE)\
	PRED(eNE*3)\
	PRED(eW+eNE-eNN)\
	PRED(eNE+eNEE-eNNEEE)\
	PRED((eWWW+eWW+eW+eN+eNE+eNEE+eNEEE)>>2)\
	PRED((3*(N+W+NW+NE)+eN+eW+eNW+eNE)>>4)\
	PRED(3*(eN-eNN)+eNNN)\
	PRED((eW*N+eN*W)>>(depth-2))\
	PRED(4*(N+NNN)-6*NN-NNNN)\
	PRED(4*(W+WWW)-6*WW-WWWW)\
	PRED(7*(eNE-eNEEEEEE)-21*(eNEE-eNEEEEE)+35*(eNEEE-eNEEEE)+eNEEEEEEE)\
	PRED(eNEEE)\
	PRED(eN+eW-(eNW+eNE))
#endif
#if 0
//8 predictors
#define PREDLIST\
	PRED((eW+kx*nlevels/image->iw-half)>>1)\
	PRED(pred)\
	PRED(eNW)\
	PRED((eN+eNE)>>1)\
	PRED(eNEEE)\
	PRED(eN+eW-((eNW+eNE)>>1))\
	PRED(eNE)\
	PRED((eN+eW)>>1)
#endif
const char *f09_prednames[]=
{
#define PRED(X) #X,
	PREDLIST
#undef  PRED
};
#ifdef LINEAR_QUANTIZER_BITS
#define QUANTIZE(X) (((X)>>(depth-LINEAR_QUANTIZER_BITS))+(1<<LINEAR_QUANTIZER_BITS>>1))
#else
static int quantize_ctx(int x)//signed
{
	int negmask=x>>31;
	x=floor_log2_32(abs(x))+1;
	x^=negmask;
	x-=negmask;
	return x;
}
#define QUANTIZE(X) quantize_ctx(X)+qhalf
#endif
int f09_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	PROF_START();
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	int depth=image->depth, nlevels=1<<depth, half=nlevels>>1, mask=nlevels-1;
#ifdef LINEAR_QUANTIZER_BITS
	int qhalf=1<<LINEAR_QUANTIZER_BITS>>1;
#else
	int qhalf=depth;
#endif
	int qlevels=qhalf<<1|1;

	size_t bufsize=(image->iw+PADX*2LL)*sizeof(short[PADY*4*2]);//PADY padded rows * 4 channels max * {pixels, errors}
	short *pixels=(short*)malloc(bufsize);

#ifdef FULL_TREE
	//8 bit -> 256*(8*2+1)*4*8*2 = 278528 bytes
	//16 bit -> 65536*(16*2+1)*4*8*2 = 132 MB
	int treesize=1<<depth;
#else
	//8 bit -> 8*(8+1)/2*(8*2+1)*4*8*2 = 39168 bytes
	//16 bit -> 16*(16+1)/2*(16*2+1)*4*8*2 = 280.5 KB
	int treesize=depth*(depth+1)>>1;
#endif
	size_t statssize=sizeof(short[4*F09_NCTX])*qlevels*treesize;
	short *stats=(short*)malloc(statssize);

	size_t mixersize=depth*sizeof(int[4*(F09_NCTX+1LL)]);
	int *mixer=(int*)malloc(mixersize);
#ifdef ENABLE_SSE
	int ssebits=7, sselevels=1<<ssebits;
	size_t ssesize=sizeof(long long[4])*depth*sselevels;
	long long *sse=(long long*)malloc(ssesize);
#endif
	if(!pixels||!stats||!mixer
#ifdef ENABLE_SSE
		||!sse
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(pixels, 0, bufsize);
	memset(stats, 0, statssize);
	for(int ks=0;ks<nlevels;++ks)
	{
		int sh=depth+depth+7;
		int weight=nlevels-1-ks;
		weight*=weight;
#ifdef FULL_TREE
		for(int kb=depth-1, idx2=1;kb>=0;--kb)
		{
			int bit=ks>>kb&1;
			short *ctr=stats+idx2;
			*ctr+=(0x8000-*ctr)*(nlevels-1-ks)>>(depth+7);
			idx2=idx2<<1|bit;
		}
#else
		for(int kb=0, MSBidx=0;kb<depth;++kb)
		{
			int bit=ks>>(depth-1-kb)&1;
			short *ctr=stats+(kb*(kb+1LL)>>1)+kb-MSBidx;

			int val=*ctr+0x8000;
			val+=(int)((0x10000LL-val)*weight>>sh);
			*ctr=(short)(CLAMP(1, val, 0xFFFF)-0x8000);

			MSBidx+=(!bit)&-(MSBidx==kb);
		}
#endif
	}
	memfill(stats+treesize, stats, statssize-treesize*sizeof(short), treesize*sizeof(short));
	*mixer=0x8000;
	memfill(mixer+1, mixer, mixersize-sizeof(*mixer), sizeof(*mixer));
#ifdef ENABLE_SSE
	memset(sse, 0, ssesize);
#endif
	PROF(INIT);
#ifdef TRACK_MIXER
	long long av_mixer[4*(F09_NCTX+1LL)]={0};
#endif
	DList list={0};
	ArithmeticCoder ec={0};
	dlist_init(&list, 1, 1024, 0);
	if(fwd)
		ac_enc_init(&ec, &list);
	else
		ac_dec_init(&ec, cbuf, cbuf+clen);
	for(int ky=0, idx=0;ky<image->ih;++ky)
	{
		short *rows[PADY];
		for(int k=0;k<PADY;++k)
			rows[k]=pixels+(((image->iw+PADX*2LL)*(((size_t)ky-k)&(PADY-1))+PADX)<<3);
		for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
		{
			short *curr=rows[0];
			if(fwd)
			{
				if(image->nch>=3)
				{
#ifdef PERMUTE_YUV
					curr[0]=image->data[idx+1];//Y
					curr[1]=image->data[idx+2];//U/Cb
					curr[2]=image->data[idx+0];//V/Cr
					curr[2]-=curr[0];		curr[2]=(short)(((curr[2]+half)&mask)-half);
					curr[1]-=curr[0];		curr[1]=(short)(((curr[1]+half)&mask)-half);
					curr[0]+=(curr[1]+curr[2])>>2;	curr[0]=(short)(((curr[0]+half)&mask)-half);
#else
					curr[0]=image->data[idx+0];
					curr[1]=image->data[idx+1];
					curr[2]=image->data[idx+2];
					curr[0]-=curr[1];		curr[0]=(short)(((curr[0]+half)&mask)-half);
					curr[2]-=curr[1];		curr[2]=(short)(((curr[2]+half)&mask)-half);
					curr[1]+=(curr[0]+curr[2])>>2;	curr[1]=(short)(((curr[1]+half)&mask)-half);
#endif
				}
				else
					memcpy(curr, image->data+idx, image->nch*sizeof(short));
			}
			for(int kc=0;kc<image->nch;++kc)
			{
				//if(idx==9555&&kc==1)//
				//if(idx==2427&&kc==2)//
				//if(idx==6&&kc==1)//
				//if(kx==image->iw/2&&ky==image->ih/2)//
				//	printf("");

				int
					NNNN		=rows[4][kc+0*8+0],
					NNN		=rows[3][kc+0*8+0],
					NNNE		=rows[3][kc+1*8+0],
					NNWWW		=rows[2][kc-3*8+0],
					NNWW		=rows[2][kc-2*8+0],
					NNW		=rows[2][kc-1*8+0],
					NN		=rows[2][kc+0*8+0],
					NNE		=rows[2][kc+1*8+0],
					NNEE		=rows[2][kc+2*8+0],
					NWW		=rows[1][kc-2*8+0],
					NW		=rows[1][kc-1*8+0],
					N		=rows[1][kc+0*8+0],
					NE		=rows[1][kc+1*8+0],
					NEE		=rows[1][kc+2*8+0],
					NEEE		=rows[1][kc+3*8+0],
					WWWWWW		=rows[0][kc-6*8+0],
					WWWWW		=rows[0][kc-5*8+0],
					WWWW		=rows[0][kc-4*8+0],
					WWW		=rows[0][kc-3*8+0],
					WW		=rows[0][kc-2*8+0],
					W		=rows[0][kc-1*8+0],
					eNNNNNNN	=rows[7][kc+0*8+4],
					eNNNNNN		=rows[6][kc+0*8+4],
					eNNNNN		=rows[5][kc+0*8+4],
					eNNNN		=rows[4][kc+0*8+4],
					eNNNNE		=rows[4][kc+1*8+4],
					eNNN		=rows[3][kc+0*8+4],
					eNNNE		=rows[3][kc+1*8+4],
					eNNNEE		=rows[3][kc+2*8+4],
					eNNNEEE		=rows[3][kc+3*8+4],
					eNNNEEEEEE	=rows[3][kc+6*8+4],
					eNNWWW		=rows[2][kc-3*8+4],
					eNNWW		=rows[2][kc-2*8+4],
					eNNW		=rows[2][kc-1*8+4],
					eNN		=rows[2][kc+0*8+4],
					eNNE		=rows[2][kc+1*8+4],
					eNNEE		=rows[2][kc+2*8+4],
					eNNEEE		=rows[2][kc+3*8+4],
					eNNEEEE		=rows[2][kc+4*8+4],
					eNNEEEEE	=rows[2][kc+5*8+4],
					eNNEEEEEE	=rows[2][kc+6*8+4],
					eNNEEEEEEE	=rows[2][kc+7*8+4],
					eNWWW		=rows[1][kc-3*8+4],
					eNWW		=rows[1][kc-2*8+4],
					eNW		=rows[1][kc-1*8+4],
					eN		=rows[1][kc+0*8+4],
					eNE		=rows[1][kc+1*8+4],
					eNEE		=rows[1][kc+2*8+4],
					eNEEE		=rows[1][kc+3*8+4],
					eNEEEE		=rows[1][kc+4*8+4],
					eNEEEEE		=rows[1][kc+5*8+4],
					eNEEEEEE	=rows[1][kc+6*8+4],
					eNEEEEEEE	=rows[1][kc+7*8+4],
					eWWWWWWW	=rows[0][kc-7*8+4],
					eWWWWWW		=rows[0][kc-6*8+4],
					eWWWWW		=rows[0][kc-5*8+4],
					eWWWW		=rows[0][kc-4*8+4],
					eWWW		=rows[0][kc-3*8+4],
					eWW		=rows[0][kc-2*8+4],
					eW		=rows[0][kc-1*8+4];
				(void)NNNN;
				(void)NNN;
				(void)NNNE;
				(void)NNWWW;
				(void)NNWW;
				(void)NNW;
				(void)NN;
				(void)NNE;
				(void)NNEE;
				(void)NWW;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)WWWWWW;
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
				(void)eNNNNNNN;
				(void)eNNNNNN;
				(void)eNNNNN;
				(void)eNNNN;
				(void)eNNNNE;
				(void)eNNN;
				(void)eNNNE;
				(void)eNNNEE;
				(void)eNNNEEE;
				(void)eNNNEEEEEE;
				(void)eNNWWW;
				(void)eNNWW;
				(void)eNNW;
				(void)eNN;
				(void)eNNE;
				(void)eNNEE;
				(void)eNNEEE;
				(void)eNNEEEE;
				(void)eNNEEEEE;
				(void)eNNEEEEEE;
				(void)eNNEEEEEEE;
				(void)eNWWW;
				(void)eNWW;
				(void)eNW;
				(void)eN;
				(void)eNE;
				(void)eNEE;
				(void)eNEEE;
				(void)eNEEEE;
				(void)eNEEEEE;
				(void)eNEEEEEE;
				(void)eNEEEEEEE;
				(void)eWWWWWWW;
				(void)eWWWWWW;
				(void)eWWWWW;
				(void)eWWWW;
				(void)eWWW;
				(void)eWW;
				(void)eW;

				int pred=N+W-NW;
				int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
				pred=CLAMP(vmin, pred, vmax);

				int ctx[]=
				{
#define PRED(X) X,
					PREDLIST
#undef  PRED

					//(eW+kx*nlevels/image->iw-half)>>1,//0
					//pred,//1
					//eNW,//2
					//(eN+eNE)>>1,//3
					//eNEEE,//4
					//eN+eW-((eNW+eNE)>>1),//5
					//eNE,//6
					//(eN+eW)>>1,//7
				};
				if(kc)
				{
					ctx[2]=rows[0][0];//luma pixel
					ctx[0]=rows[0][4];//luma error (strong helper)
					ctx[3]=rows[0][4];
					ctx[1]=rows[0][4]*3;
					ctx[6]=rows[0][4]*7;
					ctx[4]=rows[0][4]+eW-rows[0][4-1*8];
					ctx[5]=rows[0][4]+eN-rows[1][4];
					ctx[6]=rows[0][4]+eN+eW-(rows[1][4]+rows[0][4-1*8]+eNW)+rows[1][4-1*8];
					if(kc==2)
					{
						ctx[7]=rows[0][4+1];//cb error
						ctx[8]=rows[0][0+1];//cb pixel
						ctx[9]=rows[0][4+1]+eN+eW-(rows[1][4+1]+rows[0][4-1*8+1]+eNW)+rows[1][4-1*8+1];
					}
				}
				for(int k=0;k<F09_NCTX;++k)
				{
					int x=ctx[k];
					x=CLAMP(-half, x, half-1);
					ctx[k]=QUANTIZE(x);
				}
				if((unsigned)f09_disable_ctx<F09_NCTX)
					ctx[f09_disable_ctx]=0;
#if 0
				int c0=quantize_ctx(N)+qhalf;
				int c1=quantize_ctx(W)+qhalf;
				int c2=quantize_ctx((N+W)>>1)+qhalf;
				int c3=W+NE-N;
				c3=quantize_ctx(CLAMP(-half, c3, half-1))+qhalf;
				int c4=quantize_ctx(eN)+qhalf;
				int c5=quantize_ctx(eW)+qhalf;
				int c6=quantize_ctx((eN+eW)>>1)+qhalf;
				int c7=eW+eNE-eN;
				c7=quantize_ctx(CLAMP(-half, c7, half-1))+qhalf;
#endif

#if 0
				if(
					(unsigned)ctx[0]>=qlevels||
					(unsigned)ctx[1]>=qlevels||
					(unsigned)ctx[2]>=qlevels||
					(unsigned)ctx[3]>=qlevels||
					(unsigned)ctx[4]>=qlevels||
					(unsigned)ctx[5]>=qlevels||
					(unsigned)ctx[6]>=qlevels||
					(unsigned)ctx[7]>=qlevels
				)//
					LOG_ERROR("");
#endif

				short *curr_stats[F09_NCTX];
				for(int k=0;k<F09_NCTX;++k)
					curr_stats[k]=stats+treesize*(qlevels*(F09_NCTX*kc+k)+ctx[k]);
				int *curr_mixer=mixer+depth*(F09_NCTX+1LL)*kc;
				PROF(CTX);

				int sym=0;
				if(fwd)
				{
					//if((unsigned)(curr[kc]+half)>=nlevels)
					//	LOG_ERROR("");
					sym=curr[kc]-pred;
					sym+=half;
					sym&=mask;
					sym-=half;
					curr[kc+4]=(short)sym;
					sym=sym<<1^(sym>>31);
				}
				int tidx=0;
#ifdef FULL_TREE
				++tidx;
#endif
#ifdef ENABLE_SSE
				int sse_idx=((kc?rows[0][4]:(eN+eW)>>1)+half)>>(depth-ssebits);
				//int sse_idx=QUANTIZE(kc?rows[0][4]:(eN+eW)>>1);
				//if(((size_t)kc*depth*sselevels+sse_idx)*sizeof(long long)>=ssesize)
				//	LOG_ERROR("");
				long long *curr_sse=sse+((size_t)kc*depth*sselevels+sse_idx);
#endif
				for(int kb=0;kb<depth;++kb)
				{
#ifdef FULL_TREE
					int idx2=tidx;
#else
					int idx2=(kb*(kb+1LL)>>1)+kb-tidx;
#endif
#ifdef LEARNED_DENOMINATOR
					long long p0=0, wsum=0;
					for(int k=0;k<F09_NCTX;++k)
					{
						wsum+=(long long)curr_mixer[k];
						p0+=(long long)curr_mixer[k]*curr_stats[k][idx2];
					}
					p0+=(wsum>>1)-1;
					p0/=wsum;
#else
					long long p0=curr_mixer[F09_NCTX];
					for(int k=0;k<F09_NCTX;++k)
						p0+=(long long)curr_mixer[k]*curr_stats[k][idx2];
					p0+=(1LL<<20)-1;
					p0>>=21;

					//int lgsum=floor_log2(wsum);
					//p0+=(1LL<<lgsum>>1)-1;
					//p0>>=lgsum;

					//p0+=wsum>>1;
					//p0/=wsum;
#endif
					p0+=0x8000;
#ifdef ENABLE_SSE
					long long sse_sum=*curr_sse>>12;
					int sse_count=*curr_sse&0xFFF;
					p0+=sse_sum/(sse_count+1LL);
#endif
					p0=CLAMP(1, p0, 0xFFFF);
					PROF(P0);
					
					//if(idx==6132&&kc==1&&kb==5)//
					//if(idx==3&&kc==0&&kb==0)//
					//	printf("");

					int bit=0;
					if(fwd)
					{
						bit=sym>>(depth-1-kb)&1;
						ac_enc_bin(&ec, (unsigned short)p0, bit);
					}
					else
					{
						bit=ac_dec_bin(&ec, (unsigned short)p0);
						sym|=bit<<(depth-1-kb);
					}
					PROF(EC);
					
#ifdef LEARNED_DENOMINATOR
					//p0 = (w1*a+w2*b)/(w1+w2)
					//L = (1/2)sq(p0_correct-p0)
					//dL/da = (p0-p0_correct)*w1/wsum
					//dL/w1 = (p0-p0_correct)*a/wsum
					long long error=p0-((long long)!bit<<16);
					//wsum*=7;
					for(int k=0;k<F09_NCTX;++k)
					{
						int m=curr_mixer[k];
						int s=curr_stats[k][idx2];
						long long mupdate=error*s/wsum;
						long long supdate=error*m/wsum;
						m-=(int)((mupdate+3)>>3);
						s-=(int)((supdate+3)>>3);
						curr_mixer[k]=m;
						curr_stats[k][idx2]=CLAMP(-0x7FFF, s, 0x7FFF);
					}
#else
					//L = (1/2)sq(p0_correct-p0)		p0 = sum: m[k]*s[k]
					//dL/dbias = p0-p0_correct
					//dL/ds[k] = (p0-p0_correct)*m[k]
					//dL/dm[k] = (p0-p0_correct)*s[k]
					int error=(int)p0-(!bit<<16);
					int bupdate=(error>>6)+(error>>31);
					bupdate=CLAMP(-128, bupdate, 128);
					bupdate=curr_mixer[F09_NCTX]-bupdate;
					bupdate=CLAMP(-0x2000, bupdate, 0x2000);
					curr_mixer[F09_NCTX]=bupdate;
					for(int k=0;k<F09_NCTX;++k)
					{
						int m=curr_mixer[k];
						int s=curr_stats[k][idx2];
						long long mupdate=(long long)error*s;
						long long supdate=(long long)error*m;
						mupdate=(mupdate+(1LL<<MIXER_SHIFT>>1)-1)>>MIXER_SHIFT;
						supdate=(supdate+(1LL<<STATS_SHIFT>>1)-1)>>STATS_SHIFT;
						m-=(int)mupdate;
						s-=(int)supdate;
						curr_mixer[k]=m;
						curr_stats[k][idx2]=(short)CLAMP(-0x7FFF, s, 0x7FFF);
					}
#endif
					curr_mixer+=F09_NCTX+1LL;
#ifdef ENABLE_SSE
					++sse_count;
					sse_sum-=error;
					if(sse_count>640)
					{
						sse_count>>=1;
						sse_sum=(sse_sum+1LL)>>1;
					}
					*curr_sse=sse_sum<<12|sse_count;
					++curr_sse;
#endif
#ifdef FULL_TREE
					tidx=tidx<<1|bit;
#else
					tidx+=(!bit)&-(tidx==kb);
#endif
					PROF(UPDATE);
				}
				if(!fwd)
				{
					sym=sym>>1^-(sym&1);
					curr[kc+4]=(short)sym;
					sym+=pred;
					sym+=half;
					sym&=mask;
					sym-=half;
					curr[kc]=(short)sym;
				}
			}
			if(!fwd)
			{
				short *rgb=dst->data+idx;
				if(image->nch>=3)
				{
#ifdef PERMUTE_YUV
					rgb[0]=curr[2];//Cr
					rgb[1]=curr[0];//Y
					rgb[2]=curr[1];//Cb
#else
					rgb[0]=curr[0];
					rgb[1]=curr[1];
					rgb[2]=curr[2];
#endif
					rgb[1]-=(rgb[0]+rgb[2])>>2;	rgb[1]=(short)(((rgb[1]+half)&mask)-half);
					rgb[2]+=rgb[1];			rgb[2]=(short)(((rgb[2]+half)&mask)-half);
					rgb[0]+=rgb[1];			rgb[0]=(short)(((rgb[0]+half)&mask)-half);
				}
				else
					memcpy(rgb, curr, image->nch*sizeof(short));
			}
			for(int k=0;k<PADY;++k)
				rows[k]+=8;
		}
#ifdef TRACK_MIXER
		for(int k=0;k<(int)_countof(av_mixer);++k)
			av_mixer[k]+=(long long)mixer[k];
#endif
	}
	if(fwd)
	{
		ac_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
		PROF(FINISH);
	}
	if(loud)
	{
		t0=time_sec()-t0;
		if(fwd)
		{
			ptrdiff_t usize=image_getBMPsize(image);
			ptrdiff_t csize=list.nobj;
			ptrdiff_t usedsize=bufsize+statssize+mixersize;
#ifdef ENABLE_SSE
			usedsize+=ssesize;
#endif
			printf("csize %12td  %10.6lf%%  CR %8.6lf  used %.2lf KB\n",
				csize,
				100.*csize/usize,
				(double)usize/csize,
				(double)usedsize/1024
			);
#ifdef TRACK_MIXER
			printf("mixer\n");
			for(int kc=0;kc<image->nch;++kc)
			{
				for(int k=0;k<F09_NCTX;++k)
					printf("%2d  %12lld  %s\n", k, av_mixer[(F09_NCTX+1LL)*kc+k], f09_prednames[k]);
				printf("\n");
			}
#endif
		}
		printf("F09  %c %15.6lf sec\n", 'D'+fwd, t0);
		prof_print();
	}
	dlist_clear(&list);
	free(pixels);
	free(stats);
	free(mixer);
#ifdef ENABLE_SSE
	free(sse);
#endif
	return 1;
}