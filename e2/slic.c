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


//nch:    Must be from 1 to 4. The channels are interleaved and packed.
//depth:  Must be from [1~16]. If depth<=8, data must be in bytes, otherwise data must be in little-endian uint16's (shorts).
//pixels: Must be unsigned integers shifted leftmost. For example:
//	A 5-bit subpixel must be stored like this: 0bXXXX_X000
//	A 14-bit subpixel must be stored like this: 0bXXXX_XXXX_XXXX_XX00
//loud: enables printing compression & timing details to standard output
int slic_encode(int iw, int ih, int nch, int depth, const void *pixels, ArrayHandle *data, int loud);
void* slic_decode(const void *data, int len, int *ret_iw, int *ret_ih, int *ret_nch, int *ret_depth);


#ifdef SLIC_IMPLEMENTATION
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>//log2
#include"ac.h"

	#define DISABLE_LOUD
//	#define USE_ANS

typedef struct SLICHeaderStruct
{
	char tag[4];//"SLIC"
	int iw, ih;
	unsigned char nch, depth;//nch: 1, 2, 3 or 4, depth: [1~16]
	short alpha;
	//unsigned char data[];
} SLICHeader;
#define FLOOR_LOG2(LOGN, N, SH, TEMP)\
	TEMP=N,\
	SH=(TEMP>=1<<8)<<3,	LOGN =SH, TEMP>>=SH,\
	SH=(TEMP>=1<<4)<<2,	LOGN+=SH, TEMP>>=SH,\
	SH=(TEMP>=1<<2)<<1,	LOGN+=SH, TEMP>>=SH,\
	SH= TEMP>=1<<1,		LOGN+=SH;
//#define SLIC_LOAD(BUF, X, Y) ((unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?BUF[iw*(ky+(Y))+kx+(X)]:0)
static void slic_rct(short *buf, int iw, int ih, int depth, int fwd)//reversible color transform: YCmCb-R
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
static void slic_predict(const short *src, short *dst, int iw, int ih, int fwd)//one channel pixels are packed
{
	const short *pixels=fwd?src:dst;
	for(int ky=0, idx=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, ++idx)
		{
			short
				N =ky?pixels[idx-iw]:0,			//SLIC_LOAD(pixels,  0, -1),
				W =kx?pixels[idx-1]:0,			//SLIC_LOAD(pixels, -1,  0),
				NW=kx&&ky?pixels[idx-iw-1]:0;	//SLIC_LOAD(pixels, -1, -1);
			short vmin, vmax, pred;
			if(N<W)
				vmin=N, vmax=W;
			else
				vmin=W, vmax=N;

			if(NW<vmin)
				pred=vmax;
			else if(NW>vmax)
				pred=vmin;
			else
				pred=(short)(N+W-NW);//shouldn't overflow, this is clamped gradient

			pred^=-fwd;//negate pred if fwd
			pred+=fwd;
			dst[idx]=src[idx]+pred;

			//if(fwd)
			//	dst[idx]=src[idx]-pred;
			//else
			//	dst[idx]=src[idx]+pred;
		}
	}
}
typedef struct CounterStruct
{
	unsigned short n[2];
} Counter;
static Counter hist[256][64];
#ifdef USE_ANS
typedef struct ProbSymStruct
{
	unsigned short p0, bit;
} ProbSym;
#endif
int slic_encode(int iw, int ih, int nch, int depth, const void *pixels, ArrayHandle *data, int loud)
{
	if(iw<1||ih<1 || nch<1||nch>4 || depth<1||depth>16 || !pixels||!data)
		return 0;
#ifndef DISABLE_LOUD
	double t_start=time_sec();
#endif
	size_t res=(size_t)iw*ih;
	short *src=(short*)malloc(res*(nch+1)*sizeof(short));
#ifdef USE_ANS
	DList prob;
	dlist_init(&prob, sizeof(ProbSym), 1024, 0);//(4*bitcount) bytes
#endif
	SLICHeader header=
	{
		{'S', 'L', 'I', 'C'},
		iw, ih,
		nch, depth,
		0,//alpha
	};
	if(!src)
	{
		LOG_ERROR2("Allocation error");
		return -1;
	}
	if(depth<=8)
	{
		const unsigned char *src0=(const unsigned char*)pixels;
		for(int kc=0;kc<nch;++kc)
		{
			for(int k=0;k<(int)res;++k)
				src[res*kc+k]=(src0[nch*k+kc]<<8)-0x8000;
		}
	}
	else
	{
		const unsigned short *src0=(const unsigned short*)pixels;
		for(int kc=0;kc<nch;++kc)
		{
			for(int k=0;k<(int)res;++k)
				src[res*kc+k]=src0[nch*k+kc]-0x8000;
		}
	}
	short *buf;
	int has_alpha=0;
	if(nch==2||nch==4)//if there is alpha: check if alpha is redundant to store it in the header
	{
		buf=src+res*(nch-1);
		short alpha=buf[0];
		for(int k=1;k<(int)res;++k)
		{
			if(buf[k]!=alpha)
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
	
#ifndef DISABLE_LOUD
	double csizes[4][7]={0};
	int usizes[4][7]={0};
#endif

	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, &header, sizeof(header));
	buf=src+res*nch;
	const short half=1<<(depth-1);
#ifdef USE_ANS
	BANSEncoder ec;
	bans_enc_init(&ec, &list);
#else
	ArithmeticCoder ec;
	ac_enc_init(&ec, &list);
#endif
	for(int kc=0;kc<nch;++kc)//for each channel
	{
		if((nch==2||nch==4)&&!has_alpha&&kc==nch-1)
			continue;
		slic_predict(src+res*kc, buf, iw, ih, 1);
		
		int bypassctr[]={1, 1};
		{
			Counter fillval={1, 1};
			memfill(hist, &fillval, sizeof(hist), sizeof(fillval));
		}
		for(int ky=0, idx=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, ++idx)
			{
				//if(kc==0&&ky==2&&kx==72)//
				//	printf("");

				short
					W=kx?buf[idx-1]:0,	//SLIC_LOAD(buf, -1, 0),
					N=ky?buf[idx-iw]:0;	//SLIC_LOAD(buf, 0, -1);
				short ctx=((N+W)/2+0x8000)>>8;
				Counter *histk=hist[ctx];
				int hidx=0;
				
				short pixel=buf[idx]>>(16-depth);

				unsigned short sym=(pixel<<1)^-((unsigned short)pixel>=half);//L-permutation
				sym&=(1<<depth)-1;
				
				int token, bypass;//the scheme for high bitdepth from JPEG XL
				int nbits;
				if(sym)
				{
					int sh, temp;
					FLOOR_LOG2(nbits, sym, sh, temp);
					++nbits;
				}
				else
					nbits=0;
				//int nbits=sym?floor_log2(sym)+1:0;
				if(nbits<=4)
					token=sym, bypass=0, nbits=0;
				else
				{
					nbits-=2;
					token=16+((nbits-3)<<1)+(sym>>nbits)-2;
					bypass=sym&((1<<nbits)-1);
				}
				for(int kb=6-1;kb>=0;--kb)//encode token (6 bits max)
				{
					Counter *ctr=histk+hidx;
					int p0=((unsigned)ctr->n[0]<<16)/(ctr->n[0]+ctr->n[1]);
					p0=CLAMP(1, p0, 0xFFFF);
					
					int bit=token>>kb&1;
#ifdef USE_ANS
					ProbSym ps={p0, bit};
					dlist_push_back1(&prob, &ps);
					//bans_enc(&ec, p0, bit);
#else
					ac_enc_bin(&ec, p0, bit);
#endif

#ifndef DISABLE_LOUD
					if(loud)
					{
						double p=(double)(bit?0x10000-p0:p0)/0x10000;
						p=-log2(p);
						int kb2=kb+1;
						csizes[kc][kb2]+=p;
						++usizes[kc][kb2];
					}
#endif

					int cnew=ctr->n[bit];
					++cnew;
					if(cnew<0x10000)
						ctr->n[bit]=cnew;
					else
					{
						ctr->n[0]>>=1;
						ctr->n[1]>>=1;
						ctr->n[0]+=!ctr->n[0];
						ctr->n[1]+=!ctr->n[1];
						++ctr->n[bit];
					}

					hidx<<=1;
					hidx+=bit+1;
				}
				for(int kb=nbits-1;kb>=0;--kb)
				{
					int p0=((unsigned)bypassctr[0]<<16)/(bypassctr[0]+bypassctr[1]);
					p0=CLAMP(1, p0, 0xFFFF);
					
					int bit=bypass>>kb&1;
#ifdef USE_ANS
					ProbSym ps={p0, bit};
					dlist_push_back1(&prob, &ps);
					//bans_enc(&ec, p0, bit);
#else
					ac_enc_bin(&ec, p0, bit);
#endif
					
#ifndef DISABLE_LOUD
					if(loud)
					{
						double p=(double)(bit?0x10000-p0:p0)/0x10000;
						p=-log2(p);
						csizes[kc][0]+=p;
						++usizes[kc][0];
					}
#endif
					
					int cnew=bypassctr[bit];//update
					++cnew;
					if(cnew<0x10000)
						bypassctr[bit]=cnew;
					else
					{
						bypassctr[0]>>=1;
						bypassctr[1]>>=1;
						bypassctr[0]+=!bypassctr[0];
						bypassctr[1]+=!bypassctr[1];
						++bypassctr[bit];
					}
					//++bypassctr[bit];
				}
			}//x loop
		}//y loop
#ifndef DISABLE_LOUD
		if(loud)
			printf("C%d Bypass p0 %lf\n", kc, (double)bypassctr[0]/(bypassctr[0]+bypassctr[1]));
#endif
	}//channel loop
#ifdef USE_ANS
	ArrayHandle probarr=0;
	dlist_appendtoarray(&prob, &probarr);
	dlist_clear(&prob);
	for(ptrdiff_t kb=probarr->count-1;kb>=0;--kb)
	{
		ProbSym *ps=(ProbSym*)array_at(&probarr, kb);
		bans_enc(&ec, ps->p0, ps->bit);
	}
	bans_enc_flush(&ec);
	array_free(&probarr);
#else
	ac_enc_flush(&ec);
#endif
	
#ifndef DISABLE_LOUD
	if(loud)
	{
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");

		double chsizes[5]={0};
		printf("\tC0\t\tC1\t\tC2\t\tC3\n\n");
		for(int kb=7-1;kb>=0;--kb)//6 token + any bypass = 7 bits
		{
			printf("B%2d  ", kb);
			for(int kc=0;kc<nch;++kc)
			{
				double csize=csizes[kc][kb];
				int usize=usizes[kc][kb];
				if(usize)
					printf(" %15.6lf", usize/csize);
				else
					printf(" %15s", "N/A");
				chsizes[kc]+=csize;
			}
			printf("\n");
		}
		printf("\n");
		chsizes[4]=chsizes[0]+chsizes[1]+chsizes[2]+chsizes[3];
		int real_nch=nch-(!has_alpha&&(nch==2||nch==4));
		printf("Total %15.6lf %15.6lf %15.6lf %15.6lf %15.6lf\n", res*depth/chsizes[0], res*depth/chsizes[1], res*depth/chsizes[2], res*depth/chsizes[3], res*depth*real_nch/chsizes[4]);
		printf("Total size\t%8d\t\t\t     %15.6lf\n", (int)list.nobj, (double)res*real_nch/list.nobj);
	}
#endif

	dlist_appendtoarray(&list, data);
	dlist_clear(&list);
	free(src);
	return 0;
}
void* slic_decode(const void *data, int len, int *ret_iw, int *ret_ih, int *ret_nch, int *ret_depth)
{
	SLICHeader header;
	if(!data||len<sizeof(header)||!ret_iw||!ret_ih||!ret_nch||!ret_depth)
	{
		LOG_ERROR2("Invalid file/args");
		return 0;
	}
	memcpy(&header, data, sizeof(header));
	int iw=header.iw, ih=header.ih, nch=header.nch, depth=header.depth, res=iw*ih;
	if(memcmp(header.tag, "SLIC", 4)||(unsigned)nch-1>4-1||(unsigned)depth-1>16-1)
	{
		LOG_ERROR2("Invalid file");
		return 0;
	}
	short *dst=(short*)malloc(res*(nch+1)*sizeof(short));
	if(!dst)
	{
		LOG_ERROR2("Allocation error");
		return 0;
	}
	short *buf=dst+res*nch;
	const unsigned char *src=(const unsigned char*)data;
#ifdef USE_ANS
	BANSDecoder ec;
	bans_dec_init(&ec, src+sizeof(header), src+len);
#else
	ArithmeticCoder ec;
	ac_dec_init(&ec, src+sizeof(header), src+len);
#endif
	for(int kc=0;kc<nch;++kc)
	{
		if((nch==2||nch==4)&&kc==nch-1&&
#ifdef USE_ANS
			ec.srcptr==ec.srcstart
#else
			ec.srcptr==ec.srcend
#endif
		)
		{
			memfill(dst+res*kc, &header.alpha, res*sizeof(short), sizeof(short));
			continue;
		}
		
		int bypassctr[]={1, 1};
		{
			Counter fillval={1, 1};
			memfill(hist, &fillval, sizeof(hist), sizeof(fillval));
		}
		for(int ky=0, idx=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, ++idx)
			{
				//if(kc==0&&ky==2&&kx==72)//
				//	printf("");
				
				short
					W=kx?buf[idx-1]:0,	//SLIC_LOAD(buf, -1, 0),
					N=ky?buf[idx-iw]:0;	//SLIC_LOAD(buf, 0, -1);
				short ctx=((N+W)/2+0x8000)>>8;//division by 2 and shift right 8
				Counter *histk=hist[ctx];
				int hidx=0;
				
				int token=0;
				for(int kb=6-1;kb>=0;--kb)//decode token (6 bits max)
				{
					Counter *ctr=histk+hidx;
					int p0=((unsigned)ctr->n[0]<<16)/(ctr->n[0]+ctr->n[1]);
					p0=CLAMP(1, p0, 0xFFFF);
					
#ifdef USE_ANS
					int bit=bans_dec(&ec, p0);
#else
					int bit=ac_dec_bin(&ec, p0);
#endif
					token|=bit<<kb;

					int cnew=ctr->n[bit];//update
					++cnew;
					if(cnew<0x10000)
						ctr->n[bit]=cnew;
					else
					{
						ctr->n[0]>>=1;
						ctr->n[1]>>=1;
						ctr->n[0]+=!ctr->n[0];
						ctr->n[1]+=!ctr->n[1];
						++ctr->n[bit];
					}

					hidx<<=1;
					hidx+=bit+1;
				}
				unsigned short sym;
				if(token<16)
					sym=token;
				else
				{
					int bypass=0, nbits=((token-16)>>1)+3;
					for(int kb=nbits-1;kb>=0;--kb)//decode bypass
					{
						int p0=((unsigned)bypassctr[0]<<16)/(bypassctr[0]+bypassctr[1]);
						p0=CLAMP(1, p0, 0xFFFF);
						
#ifdef USE_ANS
						int bit=bans_dec(&ec, p0);
#else
						int bit=ac_dec_bin(&ec, p0);
#endif
						bypass|=bit<<kb;
						
						int cnew=bypassctr[bit];//update
						++cnew;
						if(cnew<0x10000)
							bypassctr[bit]=cnew;
						else
						{
							bypassctr[0]>>=1;
							bypassctr[1]>>=1;
							bypassctr[0]+=!bypassctr[0];
							bypassctr[1]+=!bypassctr[1];
							++bypassctr[bit];
						}
						//++bypassctr[bit];
					}
					sym=1<<(nbits+1)|(token&1)<<nbits|bypass;
				}

				short pixel=(sym>>1)^-(sym&1);//inverse L-permutation
				buf[idx]=pixel<<(16-depth);
			}//x-loop
		}//y-loop
		slic_predict(buf, dst+res*kc, iw, ih, 0);
	}//channel-loop
	
	if(nch==3||nch==4)
		slic_rct(dst, iw, ih, depth, 0);

	int pxsize=depth<=8?sizeof(char):sizeof(short);
	void *ret=malloc(res*nch*pxsize);
	if(!ret)
	{
		LOG_ERROR2("Allocation error");
		return 0;
	}

	//pack pixels, interleave channels, and add half
	if(depth<=8)
	{
		unsigned char *dst2=(unsigned char*)ret;
		for(int kc=0;kc<nch;++kc)
		{
			for(int k=0;k<res;++k)
				dst2[nch*k+kc]=(dst[res*kc+k]>>8)+0x80;
		}
	}
	else
	{
		unsigned short *dst2=(unsigned short*)ret;
		for(int kc=0;kc<nch;++kc)
		{
			for(int k=0;k<res;++k)
				dst2[nch*k+kc]=dst[res*kc+k]+0x8000;
		}
	}
	free(dst);
	*ret_iw=iw;
	*ret_ih=ih;
	*ret_nch=nch;
	*ret_depth=depth;
	return ret;
}
#undef  SLIC_LOAD
#endif//SLIC_IMPLEMENTATION

#ifdef __cplusplus
}
#endif
#endif//_INC_SLIC_H
