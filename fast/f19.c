#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;


	#define JPEG2000_RCT
	#define DWP


#include"ac.h"

typedef struct SymbolStruct
{
	//unsigned short cdf, freq;
	char token, nbits;
	short bypass, ctx;
} Symbol;
//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
#define CONFIG_EXP 5
#define CONFIG_MSB 2
#define CONFIG_LSB 0
static void quantize_pixel(int x, int *token, int *bypass, int *nbits)
{
	if(x<(1<<CONFIG_EXP))
	{
		*token=x;//token
		*nbits=0;
		*bypass=0;
	}
	else
	{
		int lgv=floor_log2_32((unsigned)x);
		int mantissa=x-(1<<lgv);
		*token = (1<<CONFIG_EXP) + (
				(lgv-CONFIG_EXP)<<(CONFIG_MSB+CONFIG_LSB)|
				(mantissa>>(lgv-CONFIG_MSB))<<CONFIG_LSB|
				(mantissa&((1<<CONFIG_LSB)-1))
			);
		*nbits=lgv-CONFIG_MSB+CONFIG_LSB;
		*bypass=x>>CONFIG_LSB&((1LL<<*nbits)-1);
	}
}
static int quantize_ctx(int x)
{
	int negmask=-(x<0);
	x=abs(x);
	x=floor_log2_32(x)+1;
	//x=floor_log2_32(x)+1;
	x^=negmask;
	x-=negmask;
	return x;
}
#ifdef DWP
#define DWP_NBITS 12
#define DWP_LGBLOCKSIZE 4
#define DWP_BLOCKSIZE (1<<DWP_LGBLOCKSIZE)
#define DWP_NPREDS 8
#define DWP_PREDLIST\
	DWP_PRED(N+W-NW)\
	DWP_PRED(W+((5*(N-NW)+NE-WW)>>3))\
	DWP_PRED((N+W+NE-NW)>>1)\
	DWP_PRED((7*W+5*(N-NW)+NE)>>3)\
	DWP_PRED((N+W)>>1)\
	DWP_PRED((W+NEE)>>1)\
	DWP_PRED((8*(N+W)+3*(NE-NW))>>4)\
	DWP_PRED((3*(N+W)-2*NW)>>2)

//	DWP_PRED((N+W+NE-NW)>>1)
//	DWP_PRED(W+(eW>>1))
//	DWP_PRED(N+(eN>>1))
//	DWP_PRED(W+NW-NWW)
//	DWP_PRED(W+NE-N)
//	DWP_PRED((3*N+W-(NNW+NNE))>>1)
//	DWP_PRED(N+W-NW+((eN+eW-eNW)>>1))
//	DWP_PRED(N+NE-NNE

//	DWP_PRED((7*W+5*(N-NW)+NE)>>3)
static void update_wp(int *weights, int *errors)
{
	for(int kc=0;kc<3;++kc)
	{
		long long
			w0=0x400000LL/((long long)errors[kc+4*0]+1LL),
			w1=0x400000LL/((long long)errors[kc+4*1]+1LL),
			w2=0x400000LL/((long long)errors[kc+4*2]+1LL),
			w3=0x400000LL/((long long)errors[kc+4*3]+1LL),
			w4=0x400000LL/((long long)errors[kc+4*4]+1LL),
			w5=0x400000LL/((long long)errors[kc+4*5]+1LL),
			w6=0x400000LL/((long long)errors[kc+4*6]+1LL),
			w7=0x400000LL/((long long)errors[kc+4*7]+1LL),
			sum=w0+w1+w2+w3+w4+w5+w6+w7;
		sum+=!sum;
		w0=((w0<<DWP_NBITS)+(sum>>1))/sum;
		w1=((w1<<DWP_NBITS)+(sum>>1))/sum;
		w2=((w2<<DWP_NBITS)+(sum>>1))/sum;
		w3=((w3<<DWP_NBITS)+(sum>>1))/sum;
		w4=((w4<<DWP_NBITS)+(sum>>1))/sum;
		w5=((w5<<DWP_NBITS)+(sum>>1))/sum;
		w6=((w6<<DWP_NBITS)+(sum>>1))/sum;
		w7=((w7<<DWP_NBITS)+(sum>>1))/sum;
		weights[kc+4*0]=(int)w0;
		weights[kc+4*1]=(int)w1;
		weights[kc+4*2]=(int)w2;
		weights[kc+4*3]=(int)w3;
		weights[kc+4*4]=(int)w4;
		weights[kc+4*5]=(int)w5;
		weights[kc+4*6]=(int)w6;
		weights[kc+4*7]=(int)w7;
	}
	memset(errors, 0, sizeof(int[4*DWP_NPREDS]));
}
#endif
#define QUANTIZE(X) (quantize_ctx(X)+(clevels>>1))
#define LOAD(BUF, C, X, Y) ((unsigned)(ky+(Y))<(unsigned)image->ih&&(unsigned)(kx+(X))<(unsigned)image->iw?BUF[(image->iw*(ky+(Y))+kx+(X))*image->nch+C]:0)
int f19_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;
	//Image *im2=0, tempimage={0};
	ANSCoder ec;
	DList list;
	dlist_init(&list, 1, 0x10000, 0);
	int depth=image->depth+2;
	UPDATE_MIN(depth, 16);
	int nlevels=1<<depth, half=nlevels>>1;
	int clevels=quantize_ctx(half)<<1|1;
	int token, bypass, nbits;
	quantize_pixel(nlevels, &token, &bypass, &nbits);
	int qlevels=token+1;
	int ncdfs=image->nch*clevels;
	size_t cdfsize=sizeof(int)*ncdfs*(qlevels+1LL);
	int *hist=(int*)malloc(cdfsize);
	size_t ebufsize=sizeof(short[4*8])*(image->iw+8LL);//4 padded rows * 4 channels max * {pixels, errors}
	short *pixels=(short*)malloc(ebufsize);
#ifdef DWP
	int prederrors[4*DWP_NPREDS]={0};
	int nblocks=(image->iw+DWP_BLOCKSIZE-1)>>DWP_LGBLOCKSIZE;
	int *weights=(int*)malloc(sizeof(int[4*DWP_NPREDS])*nblocks);
	if(!weights)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	*weights=(1<<DWP_NBITS)/DWP_NPREDS;
	memfill(weights+1, weights, sizeof(int[4*DWP_NPREDS])*nblocks-sizeof(int), sizeof(int));
#endif
	if(!hist||!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(hist, 0, cdfsize);
	memset(pixels, 0, ebufsize);
	unsigned *CDFs=(unsigned*)hist;
	int perm[]={1, 2, 0, 3};
	short curr[4]={0};
	int ctx;
	unsigned *curr_CDF;
	unsigned cdf;
	int freq;
	if(fwd)
	{
		Symbol *sym;
		int nvals=image->nch*image->iw*image->ih;
		size_t bufsize=sizeof(Symbol)*nvals;
		Symbol *buf=(Symbol*)malloc(bufsize);
		//im2=&tempimage;
		//image_copy(im2, image);
		if(!buf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(buf, 0, bufsize);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
#ifdef DWP
			int *wk=weights;
#endif
			short *rows[]=
			{
				pixels+(((image->iw+8LL)*((ky-0LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-1LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-2LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-3LL)&3)+4)<<3),
			};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				memcpy(curr, image->data+idx, sizeof(short)*image->nch);
#ifdef JPEG2000_RCT
				if(image->nch>=3)
				{
					curr[0]-=curr[1];		//curr[0]=((curr[0]+half)&(nlevels-1))-half;
					curr[2]-=curr[1];		//curr[2]=((curr[2]+half)&(nlevels-1))-half;
					curr[1]+=(curr[0]+curr[2])>>2;	//curr[1]=((curr[1]+half)&(nlevels-1))-half;
				}
#endif
				for(int kc0=0;kc0<image->nch;++kc0)
				{
					int
						kc	=perm[kc0],
						NNW	=rows[2][kc-2*8+0],
						NNE	=rows[2][kc+2*8+0],
						NWW	=rows[1][kc-2*8+0],
						NW	=rows[1][kc-1*8+0],
						N	=rows[1][kc+0*8+0],
						NE	=rows[1][kc+1*8+0],
						NEE	=rows[1][kc+2*8+0],
						WW	=rows[0][kc-2*8+0],
						W	=rows[0][kc-1*8+0],
						eNW	=rows[1][kc-1*8+4],
						eN	=rows[1][kc+0*8+4],
						eNE	=rows[1][kc+1*8+4],
						eW	=rows[0][kc-1*8+4],
						offset	=0;
#ifdef DWP
					int preds[]=
					{
#define DWP_PRED(X) X,
						DWP_PREDLIST
#undef  DWP_PRED
					};
					int pred=0;
					for(int k=0;k<DWP_NPREDS;++k)
						pred+=wk[k<<2|kc]*preds[k];
					pred+=1<<DWP_NBITS>>1;
					pred>>=DWP_NBITS;
					int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
					UPDATE_MIN(vmin, NE);
					UPDATE_MAX(vmax, NE);
					pred=CLAMP(vmin, pred, vmax);
#else
					int pred=N+W-NW;
					int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
					pred=CLAMP(vmin, pred, vmax);
					pred=((pred<<2)-pred+((N+W)>>1))>>2;
					//pred=(6*pred+N+W)>>3;
#endif
#ifndef JPEG2000_RCT
					if(kc0)
					{
						offset=rows[0][1];
						if(kc0>1)
							offset=(2*offset+rows[0][2])>>1;
						pred+=offset;
						pred=CLAMP(-half, pred, half-1);
					}
#endif
					ctx=abs(eN)>abs(eW)?eN:eW;	//52283528
					//ctx=MAXVAR(abs(eN), abs(eW));	//52358926
					//ctx=(abs(eN)+abs(eW))>>1;	//52391324
					ctx=QUANTIZE(ctx);

					int val=curr[kc];
#ifdef DWP
					for(int k=0;k<DWP_NPREDS;++k)
						prederrors[k]+=abs(preds[k]-val);
#endif
					//int val=image->data[idx+kc];
					rows[0][kc]=val-offset;
					val-=pred;
					//val+=half;
					//val&=nlevels-1;
					//val-=half;
					rows[0][kc+4]=val;
					val=val<<1^-(val<0);
					quantize_pixel(val, &token, &bypass, &nbits);
					//if(token>=qlevels)//
					//	LOG_ERROR("Token OOB %d/%d", token, qlevels);

					++hist[(clevels*kc0+ctx)*(qlevels+1)+token];
					buf[idx+kc0].token=token;
					buf[idx+kc0].nbits=nbits;
					buf[idx+kc0].bypass=bypass;
					buf[idx+kc0].ctx=ctx;
				}
				rows[0]+=8;
				rows[1]+=8;
				rows[2]+=8;
				rows[3]+=8;
#ifdef DWP
				if(kx&&!(kx&(DWP_BLOCKSIZE-1)))
				{
					update_wp(wk, prederrors);
					wk+=4*DWP_NPREDS;
				}
#endif
			}
		}
		for(int kt=0;kt<ncdfs;++kt)
		{
			int *curr_hist=hist+(size_t)kt*(qlevels+1LL);
			curr_CDF=CDFs+(size_t)kt*(qlevels+1LL);
			int sum=0;
			for(int ks=0;ks<qlevels;++ks)
				sum+=curr_hist[ks];
			if(sum)
			{
				int c=0;
				for(int ks=0;ks<qlevels;++ks)
				{
					int freq=curr_hist[ks];
					curr_CDF[ks]=c*(0x10000LL-qlevels)/sum+ks;
					c+=freq;
					dlist_push_back(&list, curr_CDF+ks, sizeof(short));
				}
				curr_CDF[qlevels]=0x10000;
				//dlist_push_back(&list, curr_CDF, sizeof(short)*qlevels);
			}
			else
			{
				//memset(curr_CDF, 0, sizeof(short)*(qlevels+1LL));
				curr_CDF[0]=0xFFFF;
				dlist_push_back(&list, curr_CDF, sizeof(short));
			}
		}
		ans_enc_init(&ec, &list);
		for(int idx=nvals-1, kc0=image->nch-1;idx>=0;--idx, --kc0)
		{
			sym=buf+idx;
			kc0+=image->nch&-(kc0<0);
			if(sym->nbits)
			{
				cdf=(sym->bypass<<16)>>sym->nbits, freq=0x10000>>sym->nbits;
				ans_enc_update(&ec, cdf, freq);
			}
			curr_CDF=CDFs+(clevels*kc0+sym->ctx)*(qlevels+1)+sym->token;
			cdf=curr_CDF[0];
			freq=curr_CDF[1]-cdf;
			ans_enc_update(&ec, cdf, freq);
		}
		ans_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
		if(loud)
		{
			t0=time_sec()-t0;
			ptrdiff_t usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
			ptrdiff_t csize=list.nobj;
			printf("Overhead %d CDFs = %zd bytes\n", ncdfs, cdfsize);
			printf("Memory usage: %17.2lf MB\n", (cdfsize+ebufsize+bufsize)/(1024.*1024.));

			printf("%14td/%14td = %10.6lf%%  CR %lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
			printf("E %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
			printf("\n");
		}
		dlist_clear(&list);
		free(buf);
	}
	else//decode
	{
		unsigned short c;
		unsigned char *curr_CDF2sym;
		size_t CDF2symsize=sizeof(char[0x10000])*ncdfs;
		unsigned char *CDF2sym=(unsigned char*)malloc(CDF2symsize);
		for(int kt=0;kt<ncdfs;++kt)
		{
			unsigned short start=0;
			memcpy(&start, cbuf, sizeof(start));
			if(start==0xFFFF)
			{
				cbuf+=sizeof(short);
				clen-=sizeof(short);
				continue;
			}
			if(start)
			{
				LOG_ERROR("Invalid bitstream");
				return 1;
			}
			curr_CDF=CDFs+(size_t)kt*(qlevels+1LL);
			for(int ks=0;ks<qlevels;++ks)//won't overflow thanks to adaptive guard at encoder
			{
				curr_CDF[ks]=0;
				memcpy(curr_CDF+ks, cbuf, sizeof(short));
				cbuf+=sizeof(short);
				clen-=sizeof(short);
			}
			curr_CDF[qlevels]=0x10000;
			curr_CDF2sym=CDF2sym+((size_t)kt<<16);
			int ks=0;
			for(unsigned c=0;c<0x10000;++c)
			{
				ks+=c>=curr_CDF[ks+1];
				curr_CDF2sym[c]=(unsigned char)ks;
			}
		}
		ans_dec_init(&ec, cbuf, cbuf+clen);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
#ifdef DWP
			int *wk=weights;
#endif
			short *rows[]=
			{
				pixels+(((image->iw+8LL)*((ky-0LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-1LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-2LL)&3)+4)<<3),
				pixels+(((image->iw+8LL)*((ky-3LL)&3)+4)<<3),
			};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				for(int kc0=0;kc0<image->nch;++kc0)
				{
					int
						kc	=perm[kc0],
						NNW	=rows[2][kc-2*8+0],
						NNE	=rows[2][kc+2*8+0],
						NWW	=rows[1][kc-2*8+0],
						NW	=rows[1][kc-1*8+0],
						N	=rows[1][kc+0*8+0],
						NE	=rows[1][kc+1*8+0],
						NEE	=rows[1][kc+2*8+0],
						WW	=rows[0][kc-2*8+0],
						W	=rows[0][kc-1*8+0],
						eNW	=rows[1][kc-1*8+4],
						eN	=rows[1][kc+0*8+4],
						eNE	=rows[1][kc+1*8+4],
						eW	=rows[0][kc-1*8+4],
						offset	=0;
#ifdef DWP
					int preds[]=
					{
#define DWP_PRED(X) X,
						DWP_PREDLIST
#undef  DWP_PRED
					};
					int pred=0;
					for(int k=0;k<DWP_NPREDS;++k)
						pred+=wk[k<<2|kc]*preds[k];
					pred+=1<<DWP_NBITS>>1;
					pred>>=DWP_NBITS;
					int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
					UPDATE_MIN(vmin, NE);
					UPDATE_MAX(vmax, NE);
					pred=CLAMP(vmin, pred, vmax);
#else
					int pred=N+W-NW;
					int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
					pred=CLAMP(vmin, pred, vmax);
					pred=((pred<<2)-pred+((N+W)>>1))>>2;
					//pred=(6*pred+N+W)>>3;
#endif
#ifndef JPEG2000_RCT
					if(kc0)
					{
						offset=rows[0][1];
						if(kc0>1)
							offset=(2*offset+rows[0][2])>>1;
						pred+=offset;
						pred=CLAMP(-half, pred, half-1);
					}
#endif
					ctx=abs(eN)>abs(eW)?eN:eW;
					//ctx=MAXVAR(abs(eN), abs(eW));
					//ctx=(abs(eN)+abs(eW))>>1;
					ctx=QUANTIZE(ctx);

					c=(unsigned short)ec.state;
					curr_CDF2sym=CDF2sym+(((size_t)clevels*kc0+ctx)<<16);
					token=curr_CDF2sym[c];

					curr_CDF=CDFs+(clevels*kc0+ctx)*(qlevels+1)+token;
					cdf=curr_CDF[0];
					freq=curr_CDF[1]-cdf;
					ans_dec_update(&ec, cdf, freq);
					
					int delta=token, bypass;
					if(delta>=(1<<CONFIG_EXP))
					{
						delta-=1<<CONFIG_EXP;
						int lsb=delta&((1<<CONFIG_LSB)-1);
						delta>>=CONFIG_LSB;
						int msb=delta&((1<<CONFIG_MSB)-1);
						delta>>=CONFIG_MSB;
						int nbits=delta+CONFIG_EXP-(CONFIG_MSB+CONFIG_LSB), n=nbits;
						c=(unsigned short)ec.state;
						bypass=c>>(16-nbits);
						cdf=(bypass<<16)>>nbits;
						freq=0x10000>>nbits;
						ans_dec_update(&ec, cdf, freq);
						delta=1;
						delta<<=CONFIG_MSB;
						delta|=msb;
						delta<<=nbits;
						delta|=bypass;
						delta<<=CONFIG_LSB;
						delta|=lsb;
					}
					delta=delta>>1^-(delta&1);
					rows[0][kc+4]=delta;
					delta+=pred;
					//delta+=half;
					//delta&=nlevels-1;
					//delta-=half;
					curr[kc]=delta;
					//image->data[idx+kc]=delta;
					rows[0][kc]=delta-offset;
#ifdef DWP
					for(int k=0;k<DWP_NPREDS;++k)
						prederrors[k]+=abs(preds[k]-delta);
#endif
				}
#ifdef JPEG2000_RCT
				if(image->nch>=3)
				{
					curr[1]-=(curr[0]+curr[2])>>2;	//curr[1]=((curr[1]+half)&(nlevels-1))-half;
					curr[2]+=curr[1];		//curr[2]=((curr[2]+half)&(nlevels-1))-half;
					curr[0]+=curr[1];		//curr[0]=((curr[0]+half)&(nlevels-1))-half;
				}
#endif
				memcpy(image->data+idx, curr, sizeof(short)*image->nch);
				rows[0]+=8;
				rows[1]+=8;
				rows[2]+=8;
				rows[3]+=8;
#ifdef DWP
				if(kx&&!(kx&(DWP_BLOCKSIZE-1)))
				{
					update_wp(wk, prederrors);
					wk+=4*DWP_NPREDS;
				}
#endif
			}
		}
		if(loud)
		{
			t0=time_sec()-t0;
			ptrdiff_t usize=((ptrdiff_t)image->iw*image->ih*image->nch*image->depth+7)>>3;
			printf("Memory usage: %17.2lf MB\n", (cdfsize+ebufsize+CDF2symsize)/(1024.*1024.));
			printf("D %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
		free(CDF2sym);
	}
	free(pixels);
	free(hist);
	return 0;
}