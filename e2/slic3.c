#include"ac.h"
#include<stdio.h>   //FILE/fopen/fread/fwrite/fclose
#include<stdlib.h>  //malloc/free		(optional) exit() in slic2_error()
#include<string.h>  //memset/memcpy
#include<stdarg.h>  //just for for slic2_error(): va_list
#include<sys/stat.h>//stat: (standard) to get file size fast
static const char file[]=__FILE__;


//	#define DISABLE_PALETTE
//	#define DISABLE_RCT
//	#define DISABLE_PREDICTOR


static void slic2_rct(short *buf, int iw, int ih, int depth, int fwd)//reversible color transform: YCmCb-R
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
				N =ky?pixels[idx-iw]:0,
				W =kx?pixels[idx-1]:0,
				NW=kx&&ky?pixels[idx-iw-1]:0;
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
		}
	}
}


#define SLIC2_HISTLEN 40//TODO this should depend on given bitdepth
typedef struct SLI2HeaderStruct
{
	char tag[4];//"SLI2"
	int iw, ih;

	unsigned char nch, depth;//nch: 1, 2, 3 or 4, depth: [1~16]
	short alpha;		//TODO (generalize) if palette is degenerate, skip encoding channel

	unsigned short palettesizes[4];//per channel, nonzero means enabled (depth must be > 8)
	unsigned short hist[4][SLIC2_HISTLEN];
} SLI2Header;
static SLI2Header slic2_header;
#define FLOOR_LOG2_16BIT(LOGN, N, SH, TEMP)\
	TEMP=N,\
	SH=(TEMP>=1<<8)<<3,	LOGN =SH, TEMP>>=SH,\
	SH=(TEMP>=1<<4)<<2,	LOGN+=SH, TEMP>>=SH,\
	SH=(TEMP>=1<<2)<<1,	LOGN+=SH, TEMP>>=SH,\
	SH= TEMP>=1<<1,		LOGN+=SH;
int slic3_encode(int iw, int ih, int nch, int depth, const void *src, ArrayHandle *ret)
{
	if(iw<1||ih<1 || nch<1||nch>4 || depth<1||depth>16 || !src)
		return 0;
	size_t res=(size_t)iw*ih;
	short *buf=(short*)malloc(res*(nch+2LL)*sizeof(short));//the image (separate channels), followed by current residues, then the residues are replaced with {bypass}, followed by {nbits, token} channel
	unsigned short *palettes[4]={0};
	if(!buf)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(slic2_header.tag, "SLI2", 4);
	slic2_header.iw=iw;
	slic2_header.ih=ih;
	slic2_header.nch=nch;
	slic2_header.depth=depth;
	unsigned char truedepth[]={depth, depth, depth, depth};
	if(depth<=8)
	{
		const unsigned char *src0=(const unsigned char*)src;
		memset(slic2_header.palettesizes, 0, sizeof(slic2_header.palettesizes));
		for(int kc=0;kc<nch;++kc)
		{
			for(int k=0;k<res;++k)
				buf[res*kc+k]=(src0[nch*k+kc]<<8)-0x8000;
		}
	}
	else
	{
#ifndef DISABLE_PALETTE
		int histlen=1<<depth;
		int *hist=(int*)malloc(histlen*sizeof(int));
		if(!hist)
		{
			LOG_ERROR("Allocation error");
			return 0;
		}
#endif
		const unsigned short *src0=(const unsigned short*)src;
		for(int kc=0;kc<nch;++kc)
		{
#ifndef DISABLE_PALETTE
			memset(hist, 0, histlen*sizeof(int));
			for(int k=0;k<res;++k)
			{
				unsigned short val=src0[nch*k+kc]>>(16-depth);
				++hist[val];
			}
			int palettesize=0;
			for(int k=0;k<histlen;++k)
				palettesize+=hist[k]!=0;
			if(palettesize<=(1<<(depth-8)))
			{
				if(palettesize>1)
				{
					int sh, temp;
					FLOOR_LOG2_16BIT(truedepth[kc], palettesize-1, sh, temp);//truedepth = number of bits in maximum unsigned symbol
					++truedepth[kc];
				}
				else
					truedepth[kc]=1;
				slic2_header.palettesizes[kc]=palettesize;
				palettes[kc]=(unsigned short*)malloc(palettesize*sizeof(short));
				for(int k=0, idx=0;k<histlen;++k)
				{
					if(hist[k])
					{
						palettes[kc][idx]=k;
						++idx;
					}
				}
				for(int k=0;k<res;++k)
				{
					unsigned short val=src0[nch*k+kc];
					int idx=0;
					int L=0, R=palettesize;
					while(L<=R)
					{
						idx=(L+R)>>1;
						if(palettes[kc][idx]<val)
							L=idx+1;
						else if(palettes[kc][idx]>val)
							R=idx-1;
						else
							break;
					}
					buf[res*kc+k]=(short)((idx-((palettesize+1)>>1))<<(16-truedepth[kc]));
				}
			}
			else
#endif
			{
				for(int k=0;k<res;++k)
					buf[res*kc+k]=src0[nch*k+kc]-0x8000;
			}
		}
#ifndef DISABLE_PALETTE
		free(hist);
#endif
	}
	short *pixels, *tokens;
	int has_alpha=0;
	if(nch==2||nch==4)//if there is alpha: check if alpha is redundant to store it in the header
	{
		pixels=buf+res*(nch-1LL);//last channel
		short alpha=pixels[0];
		for(int k=1;k<res;++k)
		{
			if(pixels[k]!=alpha)
			{
				has_alpha=1;
				break;
			}
		}
		if(!has_alpha)
			slic2_header.alpha=alpha;
	}
#ifndef DISABLE_RCT
	if(nch==3||nch==4)
	{
		int mindepth=truedepth[0];//don't let RCT increase true depth in case of using palette
		if(mindepth>truedepth[1])mindepth=truedepth[1];
		if(mindepth>truedepth[2])mindepth=truedepth[2];
		slic2_rct(buf, iw, ih, mindepth, 1);
	}
#endif
	
	unsigned CDF2[SLIC2_HISTLEN+1];
	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, &slic2_header, sizeof(slic2_header));
	for(int kc=0;kc<4;++kc)//insert palettes, if any
	{
		if(palettes[kc])
			dlist_push_back(&list, palettes[kc], slic2_header.palettesizes[kc]*sizeof(short));
	}
	pixels=buf+res*nch;
	tokens=pixels+res;
	ANSEncoder ec;
	ans_enc_init(&ec, &list);
	for(int kc=nch-1;kc>=0;--kc)//for each channel
	{
		if((nch==2||nch==4)&&!has_alpha&&kc==nch-1)
			continue;
		const unsigned short half=1<<(truedepth[kc]-1);
#ifdef DISABLE_PREDICTOR
		pixels=buf+res*kc;
#else
		slic_predict(buf+res*kc, pixels, iw, ih, 1);
#endif
		memset(CDF2, 0, SLIC2_HISTLEN*sizeof(int));
		for(int ky=0, idx=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, ++idx)
			{
				//if(kc==1&&idx==iw*2+506)
				//if(kc==0&&idx==2908)//
				//	printf("");

				short pixel=pixels[idx]>>(16-truedepth[kc]);

				unsigned short sym=(pixel<<1)^-((unsigned short)pixel>=half);//L-permutation
				sym&=(1<<truedepth[kc])-1;
				
				int token, bypass;//the scheme for high bitdepth from JPEG XL
				int nbits;
				if(sym)
				{
					int sh, temp;
					FLOOR_LOG2_16BIT(nbits, sym, sh, temp);
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
				pixels[idx]=bypass;
				tokens[idx]=nbits<<8|token;
				//if(token>=SLIC2_HISTLEN)//
				//{
				//	LOG_ERROR("Token error");
				//	return 0;
				//}
				++CDF2[token];
			}//x-loop
		}//y-loop

		for(int kt=0, sum=0;kt<SLIC2_HISTLEN;++kt)
		{
			unsigned freq=CDF2[kt];
			CDF2[kt]=(unsigned)((unsigned long long)sum*(0x10000LL-SLIC2_HISTLEN)/res)+kt;
			sum+=freq;

			slic2_header.hist[kc][kt]=(unsigned short)CDF2[kt];
		}
		CDF2[SLIC2_HISTLEN]=0x10000;

		for(int k=(int)res-1;k>=0;--k)
		{
			int token=tokens[k]&0xFF, nbits=tokens[k]>>8, bypass=pixels[k];
			if(nbits)
				ans_enc(&ec, bypass, 0, 1<<nbits);
			ans_enc(&ec, token, CDF2, SLIC2_HISTLEN);
		}
	}//channel-loop
	ans_enc_flush(&ec);
	size_t start=dlist_appendtoarray(&list, ret);
	dlist_clear(&list);
	free(buf);
	for(int kc=0;kc<4;++kc)
	{
		if(palettes[kc])
			free(palettes[kc]);
	}
	SLI2Header *header2=(SLI2Header*)(ret[0]->data+start);
	memcpy(header2->hist, slic2_header.hist, sizeof(header2->hist));
	return 1;
}
void* slic3_decode(const unsigned char *src, int len, int *ret_iw, int *ret_ih, int *ret_nch, int *ret_depth, int *ret_dummy_alpha, int force_alpha)
{
	if(!src||len<sizeof(SLI2Header)||!ret_iw||!ret_ih||!ret_nch||!ret_depth)
	{
		LOG_ERROR("Invalid file/args");
		return 0;
	}
	memcpy(&slic2_header, src, sizeof(slic2_header));
	int iw=slic2_header.iw, ih=slic2_header.ih, nch=slic2_header.nch, depth=slic2_header.depth, res=iw*ih;
	if(memcmp(slic2_header.tag, "SLI2", 4)||(unsigned)nch-1>4-1||(unsigned)depth-1>16-1)
	{
		LOG_ERROR("Invalid file");
		return 0;
	}
	unsigned char truedepth[4];
	const unsigned char *srcstart=src+sizeof(slic2_header);
	unsigned short *palettes[4]={0};
	for(int kc=0;kc<4;++kc)
	{
		if(slic2_header.palettesizes[kc]>0)
		{
			int bytesize=slic2_header.palettesizes[kc]*sizeof(short);
			palettes[kc]=(unsigned short*)malloc(bytesize);
			if(!palettes[kc])
			{
				LOG_ERROR("Allocation error");
				return 0;
			}
			memcpy(palettes[kc], srcstart, bytesize);
			srcstart+=bytesize;

			{
				int sh, temp;
				FLOOR_LOG2_16BIT(truedepth[kc], slic2_header.palettesizes[kc]-1, sh, temp);
				++truedepth[kc];
			}
		}
		else
			truedepth[kc]=slic2_header.depth;
	}
	short *dst=(short*)malloc(res*(nch+1LL)*sizeof(short));
	if(!dst)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	short *buf=dst+res*nch;
	unsigned CDF2[SLIC2_HISTLEN+1]={0};
	ANSDecoder ec;
	ans_dec_init(&ec, srcstart, src+len);
	for(int kc=0;kc<nch;++kc)
	{
		if((nch==2||nch==4)&&kc==nch-1&&ec.srcptr==ec.srcstart)
		{
			memfill(dst+res*kc, &slic2_header.alpha, res*sizeof(short), sizeof(short));
			continue;
		}
#ifdef DISABLE_PREDICTOR
		buf=dst+res*kc;
#endif

		for(int kt=0, overflow=0;kt<SLIC2_HISTLEN;++kt)
		{
			if(overflow)
				CDF2[kt]=0x10000-SLIC2_HISTLEN+kt;
			else
			{
				unsigned cdf=slic2_header.hist[kc][kt];
				CDF2[kt]=cdf;
				if(kt<SLIC2_HISTLEN-1)
					overflow|=cdf>slic2_header.hist[kc][kt+1];
			}
		}
		CDF2[SLIC2_HISTLEN]=0x10000;

		for(int k=0;k<res;++k)
		{
			//if(kc==1&&k==iw*2+506)
			//if(kc==0&&k==2908)//
			//if(kc==0&&k==2633)//
			//	printf("");

			int token=ans_dec(&ec, CDF2, SLIC2_HISTLEN);
			unsigned short sym;
			if(token<16)
				sym=token;
			else
			{
				int nbits=((token-16)>>1)+3;
				int bypass=ans_dec(&ec, 0, 1<<nbits);
				sym=1<<(nbits+1)|(token&1)<<nbits|bypass;
			}

			short pixel=(sym>>1)^-(sym&1);//inverse L-permutation
			buf[k]=pixel<<(16-truedepth[kc]);
		}
#ifndef DISABLE_PREDICTOR
		slic_predict(buf, dst+res*kc, iw, ih, 0);
#endif
	}//channel-loop
	
#ifndef DISABLE_RCT
	if(nch==3||nch==4)
	{
		int mindepth=truedepth[0];//don't let RCT increase true depth in case of using palette
		if(mindepth>truedepth[1])mindepth=truedepth[1];
		if(mindepth>truedepth[2])mindepth=truedepth[2];
		slic2_rct(dst, iw, ih, mindepth, 0);
	}
#endif

	for(int kc=0;kc<4;++kc)
	{
		if(palettes[kc])
		{
			short *pixels=dst+kc*res;
			for(int k=0;k<res;++k)
			{
				short val=pixels[k];
				val>>=16-truedepth[kc];
				val+=(slic2_header.palettesizes[kc]+1)>>1;
				if((unsigned)val>=slic2_header.palettesizes[kc])
				{
					LOG_ERROR("Palette error");
					return 0;
				}
				pixels[k]=(short)(palettes[kc][val]-0x8000);
			}
			free(palettes[kc]);
		}
	}

	int pxsize=depth<=8?sizeof(char):sizeof(short);
	int add_alpha=nch==3&&force_alpha;
	void *ret=malloc((size_t)res*(nch+add_alpha)*pxsize);
	if(!ret)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	//pack pixels, interleave channels, and add half
	if(depth<=8)
	{
		unsigned char *dst2=(unsigned char*)ret;
		if(add_alpha)
		{
			int black=0xFF000000;
			memfill(dst2, &black, (size_t)res*4*sizeof(char), sizeof(black));
		}
		for(int kc=0;kc<nch;++kc)
		{
			for(int k=0;k<res;++k)
				dst2[nch*k+kc]=(dst[res*kc+k]>>8)+0x80;
		}
	}
	else
	{
		unsigned short *dst2=(unsigned short*)ret;
		if(add_alpha)
		{
			long long black=0xFFFF000000000000;
			memfill(dst2, &black, (size_t)res*4*sizeof(short), sizeof(black));
		}
		for(int kc=0;kc<nch;++kc)
		{
			for(int k=0;k<res;++k)
				dst2[nch*k+kc]=dst[res*kc+k]+0x8000;
		}
	}
	free(dst);
	if(ret_iw           )*ret_iw=iw;
	if(ret_ih           )*ret_ih=ih;
	if(ret_nch          )*ret_nch=add_alpha?4:nch;
	if(ret_depth        )*ret_depth=depth;
	if(ret_dummy_alpha  )*ret_dummy_alpha=add_alpha;
	return ret;
}
