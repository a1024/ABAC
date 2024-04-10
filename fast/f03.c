#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"ac.h"
static const char file[]=__FILE__;

	#define PACK_SIGNBIT
	#define ESTIMATE_CSIZES

int f03_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	Image im2={0};
	image_copy_nodata(&im2, image);
	if(!im2.data)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	char depths[4]={0};
	memfill(depths, &im2.depth, sizeof(char[4]), sizeof(char));
	if(im2.nch>=3)
	{
		depths[1]+=depths[1]<16;
		depths[2]+=depths[2]<16;
	}
	unsigned short *stats=(unsigned short*)malloc(sizeof(short)<<MAXVAR(depths[0], depths[1]));
	if(!stats)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	if(fwd)
	{
		memcpy(im2.data, image->data, sizeof(short)*im2.iw*im2.ih*im2.nch);
		rct_JPEG2000_32(&im2, 1);
		pred_clampgrad(&im2, 1, depths);
#ifdef PACK_SIGNBIT
		for(ptrdiff_t k=0, n=(ptrdiff_t)im2.iw*im2.ih*im2.nch;k<n;++k)
		{
			//if(k==383188)
			//if(k==365032)
			//	printf("");
			int val=im2.data[k];
			val=val<<1^-(val<0);
			im2.data[k]=val;
		}
#endif
	}
	else
		memset(im2.data, 0, sizeof(short)*im2.iw*im2.ih*im2.nch);
#ifdef ESTIMATE_CSIZES
	double csizes[64]={0};
#endif
	DList list;
	dlist_init(&list, 1, 256, 0);
	ArithmeticCoder ec;
	if(fwd)
		ac_enc_init(&ec, &list);
	else
		ac_dec_init(&ec, cbuf, cbuf+clen);
	int dx=im2.nch;
	int dy=im2.nch*im2.iw;
	const int sh[]=
	{
		7,	// 0
		7,	// 1
		7,	// 2
		7,	// 3
		7,	// 4
		6,	// 5
		5,	// 6
		5,	// 7
		5,	// 8
		5,	// 9
		5,	//10
		5,	//11
		5,	//12
		5,	//13
		5,	//14
		5,	//15
	};
	for(int kc=0;kc<image->nch;++kc)
	{
		int depth=depths[kc], nlevels=1<<depth;
		//int half=nlevels>>1;
		int fillval=0x8000;
		memfill(stats, &fillval, sizeof(short)<<depth, sizeof(short));
		//for(int kb=0;kb<depth;++kb)
		for(int kb=depth-1;kb>=0;--kb)
		{
			//int MSBoffset=half&-(kb==depth-1);
			//unsigned
			//	pmaska=(1<<(kb+1))-1, pmaskb=0,
			//	fmaska=(1<< kb   )-1, fmaskb=0;
			//unsigned pmask=~0U<<kb, fmask=~0U<<(kb+1);
			unsigned
				pmaska=~0U<< kb   , pmaskb=~pmaska>>1,
				fmaska=~0U<<(kb+1), fmaskb=~fmaska>>1;
			for(int ky=0, idx=kc;ky<im2.ih;++ky)
			{
				for(int kx=0;kx<im2.iw;++kx, idx+=im2.nch)
				{
#ifdef PACK_SIGNBIT
#define LOAD(MA, MB, X, Y) ((unsigned)(ky+(Y))<(unsigned)im2.ih&&(unsigned)(kx+(X))<(unsigned)im2.iw?(unsigned short)im2.data[idx+(Y)*dy+(X)*dx]&MA:0)|MB
#else
#define LOAD(MA, MB, X, Y) ((unsigned)(ky+(Y))<(unsigned)im2.ih&&(unsigned)(kx+(X))<(unsigned)im2.iw?(im2.data[idx+(Y)*dy+(X)*dx]+half)&MA:0)|MB
#endif
					int
						NW	=LOAD(pmaska, pmaskb, -1, -1),
						N	=LOAD(pmaska, pmaskb,  0, -1),
						NE	=LOAD(pmaska, pmaskb,  1, -1),
						W	=LOAD(pmaska, pmaskb, -1,  0),
						curr	=LOAD(fmaska, fmaskb,  0,  0),
						E	=LOAD(fmaska, fmaskb,  1,  0),
						SW	=LOAD(fmaska, fmaskb, -1,  1),
						S	=LOAD(fmaska, fmaskb,  0,  1),
						SE	=LOAD(fmaska, fmaskb,  1,  1);
#undef  LOAD
#if 1
					int predNW=N+W-NW;	predNW=MEDIAN3(N, W, predNW);
					int predSW=S+W-SW;	predSW=MEDIAN3(S, W, predSW);
					int predSE=S+E-SE;	predSE=MEDIAN3(S, E, predSE);
					int predNE=N+E-NE;	predNE=MEDIAN3(N, E, predNE);
					int pred=(8*curr-(predNW+predSW+predSE+predNE))>>2;
					pred=CLAMP(0, pred, nlevels-1);
#endif
#if 0
					int pred=4*curr-(N+W+S+E)+((NW+NE+SW+SE)>>2);
					//pred+=half;
					pred=CLAMP(0, pred, nlevels-1);
					//pred=CLAMP(-half, pred, half-1);
#endif
#if 0
					int dist=abs(NW-curr), dmin=dist, pred=(NW+curr)>>1;
#define UPDATE_PRED(X) dist=abs(X-curr); if(dmin>dist)dmin=dist, pred=(X+curr)>>1
					UPDATE_PRED(N);
					UPDATE_PRED(NE);
					UPDATE_PRED(W);
					UPDATE_PRED(E);
					UPDATE_PRED(SW);
					UPDATE_PRED(S);
					UPDATE_PRED(SE);
#undef  UPDATE_PRED
#endif
#if 0
					int pred=0;
					int grads[]=
					{
						abs(NW	-curr),
						abs(N	-curr),
						abs(NE	-curr),
						abs(W	-curr),
						abs(E	-curr),
						abs(SW	-curr),
						abs(S	-curr),
						abs(SE	-curr),
					};
					int vmin=grads[0], gidx=0;
					for(int k=1;k<_countof(grads);++k)
					{
						if(vmin>grads[k])
							vmin=grads[k], gidx=k;
					}
					switch(gidx)
					{
					case 0:pred=(NW	+curr)>>1;break;
					case 1:pred=(N	+curr)>>1;break;
					case 2:pred=(NE	+curr)>>1;break;
					case 3:pred=(W	+curr)>>1;break;
					case 4:pred=(E	+curr)>>1;break;
					case 5:pred=(SW	+curr)>>1;break;
					case 6:pred=(S	+curr)>>1;break;
					case 7:pred=(SE	+curr)>>1;break;
					}
#endif
					//int pred=curr;
					//int pred=(NW+N+NE+W+curr+E+SW+S+SE+4)/9;//X

					int p0=stats[pred];
					
					//if(!kx&&!ky)//
					//if(kc==0&&kb==15&kx==297&&ky==190)//
					//if(kc==0&&kb==15&kx==258&&ky==182)//
					//{
					//	if(fwd)
					//		printf("");
					//	else
					//		printf("");
					//}

					short *val=im2.data+idx;
					int bit;
					if(fwd)
					{
#ifdef PACK_SIGNBIT
						bit=*val>>kb&1;
#else
						bit=(*val+half)>>kb&1;
#endif
						ac_enc_bin(&ec, p0, bit);
#ifdef ESTIMATE_CSIZES
						int p=bit?0x10000-p0:p0;
						csizes[kc<<4|kb]-=log2((double)p/0x10000);
#endif
					}
					else
					{
						bit=ac_dec_bin(&ec, p0);
						*val|=bit<<kb;
#ifndef PACK_SIGNBIT
						*val-=MSBoffset;
#endif
					}
					p0+=((!bit<<16)-p0)>>sh[kb];
					p0=CLAMP(1, p0, 0xFFFF);
					stats[pred]=p0;
				}
			}
		}
	}
	if(fwd)
	{
		ac_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
	}
	else
	{
		//int half[]=
		//{
		//	1<<depths[0]>>1,
		//	1<<depths[1]>>1,
		//	1<<depths[2]>>1,
		//	1<<depths[3]>>1,
		//};
		//for(ptrdiff_t k=0, nvals=(ptrdiff_t)im2.iw*im2.ih*im2.nch;k<nvals;k+=im2.nch)
		//{
		//	for(int kc=0;kc<im2.nch;++kc)
		//		dst->data[k+kc]=im2.data[k+kc];//-half[kc];
		//}
#ifdef PACK_SIGNBIT
		for(ptrdiff_t k=0, n=(ptrdiff_t)im2.iw*im2.ih*im2.nch;k<n;++k)
		{
			int val=(unsigned short)im2.data[k];
			val=val>>1^-(val&1);
			//if((val&0xFFFF)==0x4000)//
			//if(k==365032)//
			//	printf("");
			dst->data[k]=(short)val;
		}
#else
		memcpy(dst->data, im2.data, sizeof(short)*im2.iw*im2.ih*im2.nch);
#endif
		pred_clampgrad(dst, 0, depths);
		rct_JPEG2000_32(dst, 0);
	}
	if(loud)
	{
		t0=time_sec()-t0;
		if(fwd)
		{
			ptrdiff_t usize=image_getBMPsize(image);
			ptrdiff_t csize=list.nobj;
#ifdef ESTIMATE_CSIZES
			double ctotal=0;
			ptrdiff_t res=(ptrdiff_t)image->iw*image->ih;
			int maxdepth=MAXVAR(depths[0], depths[1]);
			printf("Bx   %14.2lf\n", res/8.);
			for(int kb=0;kb<maxdepth;++kb)
			{
				printf("B%2d ", kb);
				for(int kc=0;kc<image->nch;++kc)
				{
					double csize=csizes[kc<<4|kb];
					printf(" %14.2lf %6.2lf%% %12.6lf%c",
						csize/8,
						csize/res*100,
						res/csize,
						kc==image->nch-1?'\n':' '
					);
					ctotal+=csize;
				}
			}
			ctotal/=8;
			printf("T    %14.2lf %6.2lf%% %12.6lf\n",
				ctotal,
				ctotal/usize*100,
				usize/ctotal
			);
#endif
			printf("csize %14lld  %10.6lf%%  %8.6lf\n", csize, 100.*csize/usize, (double)usize/csize);
		}
		printf("%c %15.6lf sec\n", 'D'+fwd, t0);
	}
	dlist_clear(&list);
	image_clear(&im2);
	free(stats);
	return 1;
}