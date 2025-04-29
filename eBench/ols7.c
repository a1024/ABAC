#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;

#define NPREDS 9
void pred_ols7(Image *src, int fwd)
{
	//interleaved  inefficient and slow
#if 0
	int weights[8][NPREDS+1]={0};
	int w2f=src->iw>>1, w2c=(src->iw+1)>>1;
	int bufsize=(w2c+8*2)*(int)sizeof(int[8*4]);//8 padded half-rows * 4 channels max
	int *pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	for(int ky=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((w2c+16LL)*(((ky-0LL)&3)+0)+8)*4-1,
			pixels+((w2c+16LL)*(((ky-1LL)&3)+0)+8)*4-1,
			pixels+((w2c+16LL)*(((ky-2LL)&3)+0)+8)*4-1,
			pixels+((w2c+16LL)*(((ky-3LL)&3)+0)+8)*4-1,
			pixels+((w2c+16LL)*(((ky-0LL)&3)+4)+8)*4-1,
			pixels+((w2c+16LL)*(((ky-1LL)&3)+4)+8)*4-1,
			pixels+((w2c+16LL)*(((ky-2LL)&3)+4)+8)*4-1,
			pixels+((w2c+16LL)*(((ky-3LL)&3)+4)+8)*4-1,
		};
		int *ptrA=src->data+4LL*src->iw*ky+0*4LL*w2f;
		int *ptrB=src->data+4LL*src->iw*ky+1*4LL*w2f;
		int idx=0;
		for(int kx=0;kx<w2f;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				++rows[0];
				++rows[1];
				++rows[2];
				++rows[3];
				++rows[4];
				++rows[5];
				++rows[6];
				++rows[7];
				if(!src->depth[kc])
					continue;
				int
					aNNNN	=rows[0+0][+0*4],
					aNNNW	=rows[3+0][-1*4],
					aNNN	=rows[3+0][+0*4],
					aNNNE	=rows[3+0][+1*4],
					aNN	=rows[2+0][+0*4],
					aNNE	=rows[2+0][+1*4],
					aNW	=rows[1+0][-1*4],
					aN	=rows[1+0][+0*4],
					aNE	=rows[1+0][+1*4],
					aNEE	=rows[1+0][+2*4],
					aNEEE	=rows[1+0][+3*4],
					aNEEEE	=rows[1+0][+4*4],
					aNEEEEE	=rows[1+0][+5*4],
					aWWWWWW	=rows[0+0][-6*4],
					aWWWWW	=rows[0+0][-5*4],
					aWWWW	=rows[0+0][-4*4],
					aWWW	=rows[0+0][-3*4],
					aWW	=rows[0+0][-2*4],
					aW	=rows[0+0][-1*4];
				int
					bNNNN	=rows[0+4][+0*4],
					bNNNW	=rows[3+4][-1*4],
					bNNN	=rows[3+4][+0*4],
					bNNNE	=rows[3+4][+1*4],
					bNN	=rows[2+4][+0*4],
					bNNE	=rows[2+4][+1*4],
					bNW	=rows[1+4][-1*4],
					bN	=rows[1+4][+0*4],
					bNE	=rows[1+4][+1*4],
					bNEE	=rows[1+4][+2*4],
					bNEEE	=rows[1+4][+3*4],
					bNEEEE	=rows[1+4][+4*4],
					bNEEEEE	=rows[1+4][+5*4],
					bWWWWWW	=rows[0+4][-6*4],
					bWWWWW	=rows[0+4][-5*4],
					bWWWW	=rows[0+4][-4*4],
					bWWW	=rows[0+4][-3*4],
					bWW	=rows[0+4][-2*4],
					bW	=rows[0+4][-1*4];
				int preds[]=
				{
					aN,
					aW,
					3*(aN-aNN)+aNNN,
					3*(aW-aWW)+aWWW,
					aW+aNE-aN,
					aN+aW-aNW,
					(aWWWWW+aWW-aW+aNNN+aN+aNEEEEE)>>2,
					aN+aNE-aNNE,
					aNEEE,

					bN,
					bW,
					3*(bN-bNN)+bNNN,
					3*(bW-bWW)+bWWW,
					bW+bNE-bN,
					bN+bW-bNW,
					(bWWWWW+bWW-bW+bNNN+bN+bNEEEEE)>>2,
					bN+bNE-bNNE,
					bNEEE,
				};
				int *currwa=weights[kc+0*4];
				int *currwb=weights[kc+1*4];
				long long preda=currwa[NPREDS];
				long long predb=currwb[NPREDS];
				for(int k=0;k<NPREDS;++k)
				{
					preda+=(long long)currwa[k]*preds[k+0*NPREDS];
					predb+=(long long)currwb[k]*preds[k+1*NPREDS];
				}
				preda+=1LL<<(src->depth[kc]*2+5)>>1;
				predb+=1LL<<(src->depth[kc]*2+5)>>1;
				preda>>=src->depth[kc]*2+5;
				predb>>=src->depth[kc]*2+5;
				int amax=aN, amin=aW;
				int bmax=bN, bmin=bW;
				if(aN<aW)amin=aN, amax=aW;
				if(bN<bW)bmin=bN, bmax=bW;
				if(amin>aNW)amin=aNW;
				if(amax<aNW)amax=aNW;
				if(bmin>bNW)bmin=bNW;
				if(bmax<bNW)bmax=bNW;
				if(amin>aNE)amin=aNE;
				if(amax<aNE)amax=aNE;
				if(bmin>bNE)bmin=bNE;
				if(bmax<bNE)bmax=bNE;
				if(amin>aNEEE)amin=aNEEE;
				if(amax<aNEEE)amax=aNEEE;
				if(bmin>bNEEE)bmin=bNEEE;
				if(bmax<bNEEE)bmax=bNEEE;
				CLAMP2(preda, amin, amax);
				CLAMP2(predb, bmin, bmax);

				int curra=ptrA[idx];
				int currb=ptrB[idx];
				int pred2a=(int)preda;
				int pred2b=(int)predb;
				if(!keyboard[KEY_SHIFT])
				{
					pred2a^=-fwd;
					pred2b^=-fwd;
					pred2a+=fwd;
					pred2b+=fwd;
					pred2a+=curra;
					pred2b+=currb;
					pred2a<<=32-src->depth[kc];
					pred2b<<=32-src->depth[kc];
					pred2a>>=32-src->depth[kc];
					pred2b>>=32-src->depth[kc];
				}
				ptrA[idx]=(int)pred2a;
				ptrB[idx]=(int)pred2b;
				if(!fwd)
				{
					curra=pred2a;
					currb=pred2b;
				}
				rows[0][0]=curra;
				rows[4][0]=currb;

				//update
				int ea=curra-(int)preda;
				int eb=currb-(int)predb;
				ea=(((ea>0)-(ea<0))<<(src->depth[kc]-5))-(((ea>4)-(ea<-4))<<(src->depth[kc]-6));
				eb=(((eb>0)-(eb<0))<<(src->depth[kc]-5))-(((eb>4)-(eb<-4))<<(src->depth[kc]-6));
				currwa[NPREDS]+=ea;
				currwb[NPREDS]+=eb;
				for(int k=0;k<NPREDS;++k)
				{
					currwa[k]+=ea*preds[k+0*NPREDS];
					currwb[k]+=eb*preds[k+1*NPREDS];
				}
			}
		}
		if(w2f<w2c)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				int
					aW	=rows[0+0][-1*4],
					bW	=rows[0+4][-1*4];
				int preda=aW;
				int predb=bW;
				int curra=ptrA[idx];
				int currb=ptrB[idx];
				if(!keyboard[KEY_SHIFT])
				{
					preda^=-fwd;
					predb^=-fwd;
					preda+=fwd;
					predb+=fwd;
					preda+=curra;
					predb+=currb;
					preda<<=32-src->depth[kc];
					predb<<=32-src->depth[kc];
					preda>>=32-src->depth[kc];
					predb>>=32-src->depth[kc];
				}
				ptrA[idx]=(int)preda;
				ptrB[idx]=(int)predb;
				if(!fwd)
				{
					curra=preda;
					currb=predb;
				}
				rows[0][0]=curra;
				rows[4][0]=currb;
			}
		}
	}
#endif

	//scalar
#if 1
	int weights[4][NPREDS+1]={0};
	int bufsize=(src->iw+8*2)*(int)sizeof(int[4*4]);//4 padded rows * 4 channels max
	int *pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	//if(loud_transforms)
	//	DisableProcessWindowsGhosting();
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8)*4-1,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8)*4-1,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8)*4-1,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8)*4-1,
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<4;++kc, ++idx)
			{
				++rows[0];
				++rows[1];
				++rows[2];
				++rows[3];
				if(!src->depth[kc])
					continue;
				int
					NNNN	=rows[0][+0*4],
					NNNW	=rows[3][-1*4],
					NNN	=rows[3][+0*4],
					NNNE	=rows[3][+1*4],
					NN	=rows[2][+0*4],
					NNE	=rows[2][+1*4],
					NW	=rows[1][-1*4],
					N	=rows[1][+0*4],
					NE	=rows[1][+1*4],
					NEE	=rows[1][+2*4],
					NEEE	=rows[1][+3*4],
					NEEEE	=rows[1][+4*4],
					NEEEEE	=rows[1][+5*4],
					WWWWWW	=rows[0][-6*4],
					WWWWW	=rows[0][-5*4],
					WWWW	=rows[0][-4*4],
					WWW	=rows[0][-3*4],
					WW	=rows[0][-2*4],
					W	=rows[0][-1*4];
				int preds[]=
				{
					N,
					W,
					3*(N-NN)+NNN,
					3*(W-WW)+WWW,
					W+NE-N,
					N+W-NW,
					(WWWWW+WW-W+NNN+N+NEEEEE)>>2,
					N+NE-NNE,
					NEEE,
				};
				//if(ky==src->ih/2&&kx==src->iw/2)//
				//	printf("");
				int *currw=weights[kc];
				int predW=currw[NPREDS];
				for(int k=0;k<NPREDS;++k)
					predW+=currw[k]*preds[k];
				predW+=1LL<<(src->depth[kc]+11)>>1;
				predW>>=src->depth[kc]+11;
			//	predW+=1LL<<(src->depth[kc]*2+5)>>1;
			//	predW>>=src->depth[kc]*2+5;
				int pred=predW;
				int vmax=N, vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NW)vmin=NW;
				if(vmax<NW)vmax=NW;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				if(vmin>NEEE)vmin=NEEE;
				if(vmax<NEEE)vmax=NEEE;
				CLAMP2(pred, vmin, vmax);

				int curr=src->data[idx];
				int pred2=(int)pred;
				if(!keyboard[KEY_SHIFT])
				{
					pred2^=-fwd;
					pred2+=fwd;
					pred2+=curr;
					pred2<<=32-src->depth[kc];
					pred2>>=32-src->depth[kc];
				}
				src->data[idx]=(int)pred2;
				if(!fwd)
					curr=pred2;
				rows[0][0]=curr;

				//update
				int e=curr-(int)pred;
				e=(((e>0)-(e<0))<<1)-((e>4)-(e<-4));
			//	e=(((e>0)-(e<0))<<(src->depth[kc]-5))-(((e>4)-(e<-4))<<(src->depth[kc]-6));
				currw[NPREDS]+=e;
				for(int k=0;k<NPREDS;++k)
					currw[k]+=e*preds[k];
			}
		}
	}
#endif
}
