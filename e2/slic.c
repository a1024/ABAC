//SLIC: A simple lossless image codec

#define SLIC_IMPLEMENTATION
//#pragma once

#ifndef _INC_SLIC_H
#define _INC_SLIC_H
#ifdef __cplusplus
extern "C"
{
#endif
#include"util.h"//DList, array


//nch:   must be from {1, 2, 3, 4}.
//depth: must be from [1~16]. If depth<=8, data must be in bytes, otherwise data must be in little-endian uint16's (shorts).
//data:  must be unsigned integers shifted leftmost. For example:
//	A 14-bit subpixel must be stored like this: 0bXXXX_XXXX_XXXX_XX00
//	A 5-bit integer must be stored like this: 0bXXXX_X000
int slic_encode(int iw, int ih, int nch, int depth, const void *pixels, ArrayHandle *data);
void* slic_decode(const void *data, int len, int *iw, int *ih, int *nch, int *depth);


#ifdef SLIC_IMPLEMENTATION
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
typedef struct SLICHeaderStruct
{
	char tag[4];//"SLIC"
	int iw, ih;
	unsigned char nch, depth;//nch: 1, 2, 3 or 4, depth: [1~16]
	short alpha;
	int bookmarks[4];
	unsigned char data[];
} SLICHeader;
#define SLIC_LOAD(BUF, X, Y) ((unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?BUF[iw*(ky+(Y))+kx+(X)]:0)
static void slic_rct(short *buf, int iw, int ih, int depth, int fwd)//reversible color transform: YCoCb
{
	size_t res=(size_t)iw*ih;
	int mask=(short)(0xFFFF0000>>depth);
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
			short
				r=buf[idx],
				g=buf[idx+res],
				b=buf[idx+(res<<1)];

			if(fwd)
			{
				r-=g;
				g+=r>>1&mask;
				b-=g;
				g+=b>>1&mask;
			}
			else
			{
				g-=b>>1&mask;
				b+=g;
				g-=r>>1&mask;
				r+=g;
			}

			buf[idx]=r;
			buf[idx+res]=g;
			buf[idx+(res<<1)]=b;
		}
	}
}
#if 0
static void slic_predict(const short *src, short *dst, int iw, int ih, int fwd)
{
	const short *pixels=fwd?src:dst;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			short
				N =SLIC_LOAD(pixels,  0, -1),
				W =SLIC_LOAD(pixels, -1,  0),
				NW=SLIC_LOAD(pixels, -1, -1);
			short vmin, vmax;
			if(N<W)
				vmin=N, vmax=W;
			else
				vmin=W, vmax=N;
			int pred=N+W-NW;
			if(pred<vmin)pred=vmin;
			if(pred>vmax)pred=vmax;

			int idx=iw*ky+kx;
			if(fwd)
				dst[idx]=src[idx]-pred;
			else
				dst[idx]=src[idx]+pred;
		}
	}
}
#endif
//static size_t slic_calc_ctx_size(int depth)
//{
//	size_t size=0;
//	int np=depth;
//	if(depth>8)
//		np=8;
//	for(int k=0;k<depth;++k)
//		size+=np<<k;
//	return size;
//}
int slic_encode(int iw, int ih, int nch, int depth, const void *pixels, ArrayHandle *data)
{
	if(iw<1||ih<1 || nch<1||nch>4 || depth<1||depth>16 || !pixels||!data)
		return 0;
	size_t res=(size_t)iw*ih;
	//size_t ctrsize=slic_calc_ctx_size(depth);
	unsigned short *ctr=(unsigned short*)malloc(0x1000000*sizeof(short));//32MB
	short *src=(short*)malloc(res*(nch+1)*sizeof(short));
	unsigned short *prob=(unsigned short*)malloc(depth*sizeof(short));
	unsigned char *dst=0;
	size_t dstlen=0;
	int has_alpha=0;
	SLICHeader header=
	{
		{'S', 'L', 'I', 'C'},
		iw, ih,
		nch, depth,
		0,
		{0, 0, 0, 0},
	};
	if(!src||!ctr||!prob)
		return -1;
	if(depth>8)
	{
		const unsigned short *src0=(const unsigned short*)pixels;
		for(int kc=0;kc<nch;++kc)
		{
			for(int k=0;k<res;++k)
				src[res*kc+k]=src0[nch*k+kc]-0x8000;
		}
		//for(int k=0;k<n;++k)
		//	src[k]=src0[k]-0x8000;
	}
	else
	{
		const unsigned char *src0=(const unsigned char*)pixels;
		for(int kc=0;kc<nch;++kc)
		{
			for(int k=0;k<res;++k)
				src[res*kc+k]=(src0[nch*k+kc]<<8)-0x8000;
		}
		//for(int k=0;k<n;++k)
		//	src[k]=(src0[k]<<8)-0x8000;
	}
	if(nch==2||nch==4)
	{
		short alpha=src[nch-1];
		for(int k=1;k<res;++k)
		{
			if(src[k*nch+nch-1]!=alpha)
			{
				has_alpha=1;
				break;
			}
		}
		if(!has_alpha)
			header.alpha=alpha;
	}
	if(nch==3||nch==4)
		slic_rct(src, iw, ih, depth, 1);
	//short *buf=src+res*nch;
	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, &header, sizeof(header));
	for(int kc=0;kc<nch;++kc)//for each channel
	{
		if((nch==2||nch==4)&&!has_alpha&&kc==nch-1)
			continue;
		short *buf=src+res*kc;
		//slic_predict(src+res*kc, buf, iw, ih, 1);
		for(int k=0;k<0x1000000;++k)
			ctr[k]=0x8000;
		unsigned state=0x10000;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx)
			{
				short
					N =SLIC_LOAD(buf,  0, -1),
					W =SLIC_LOAD(buf, -1,  0);
				unsigned short pred=((N+W+1)>>1)+0x8000;
				if(depth<8)
					pred>>=16-depth;
				else
					pred>>=8;

				//if(kx==(iw>>1)&&ky==(ih>>1))//
				//	printf("");

				int idx=iw*ky+kx;
				int ctridx=0;
				int bitidx=0;
				for(int kb=16-1;kb>=16-depth;--kb, ++bitidx)
				{
					unsigned short *p0=ctr+(ctridx<<8|pred);

					prob[bitidx]=*p0;
					int bit=buf[idx]>>kb&1;

					int p0_new=*p0+(((!bit<<16)-*p0)>>4);
					if(p0_new<     1)p0_new=     1;
					if(p0_new>0xFFFF)p0_new=0xFFFF;
					*p0=p0_new;

					ctridx<<=1;
					++ctridx;
					ctridx+=bit;
				}//bit loop

				for(int kb=16-depth;kb<16;++kb)
				{
					--bitidx;
					int p0=prob[bitidx];
					int bit=buf[idx]>>kb&1;
					int cdf=bit?p0:0, freq=bit?0x10000-p0:p0;

					if(state>=(unsigned)(freq<<16))//renorm
					{
						dlist_push_back(&list, &state, 2);
						state>>=16;
					}
					state=state/freq<<16|(cdf+state%freq);//update
				}//bit loop
			}//x-loop
		}//y-loop
		dlist_push_back(&list, &state, 4);
		header.bookmarks[kc]=(int)list.nobj;
	}//channel loop

	size_t startidx=dlist_appendtoarray(&list, data);
	SLICHeader *header2=(SLICHeader*)(data[0]->data+startidx);
	memcpy(header2->bookmarks, header.bookmarks, sizeof(header2->bookmarks));

	free(src);
	free(ctr);
	return 0;
}
void* slic_decode(const void *data, int len, int *iw, int *ih, int *nch, int *depth)
{
	return 0;
}
#undef  SLIC_LOAD
#endif//SLIC_IMPLEMENTATION

#ifdef __cplusplus
}
#endif
#endif//_INC_SLIC_H
