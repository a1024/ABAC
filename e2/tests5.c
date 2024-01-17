#include"e2.h"
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
static const char file[]=__FILE__;

	#define SLIC4_DISABLE_PALETTE

#define SLIC4_CONFIG_EXP 6
#define SLIC4_CONFIG_MSB 3
#define SLIC4_CONFIG_LSB 1

#define SLIC4_HISTLEN 240
static const int qlevels[]=
{
	1, 3, 5, 7, 11, 15, 23, 31, 47, 63, 95, 127, 191, 255, 392, 500
};
typedef struct SLIC4HeaderStruct
{
	unsigned short palettesizes[4];
	unsigned short hist[_countof(qlevels)+1][SLIC4_HISTLEN];
} SLIC4Header;
static SLIC4Header header;
static unsigned g_hist[_countof(qlevels)+1][SLIC4_HISTLEN+1];

typedef struct TempHybridStruct
{
	unsigned short token, nbits;
	unsigned bypass;
} TempHybrid;

//from libjxl		packsign(pixel) = 0b00001MMBB...BBL	token = offset + 0bGGGGMML,  where G = bits of lg(packsign(pixel)),  bypass = 0bBB...BB
void hybriduint_encode(int val, TempHybrid *hu)
{
	int token, bypass, nbits;
	val=(val<<1)^-(val<0);//pack sign
	if(val<(1<<SLIC4_CONFIG_EXP))
	{
		token=val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=floor_log2_32((unsigned)val);
		int mantissa=val-(1<<lgv);
		token = (1<<SLIC4_CONFIG_EXP) + (
				(lgv-SLIC4_CONFIG_EXP)<<(SLIC4_CONFIG_MSB+SLIC4_CONFIG_LSB)|
				(mantissa>>(lgv-SLIC4_CONFIG_MSB))<<SLIC4_CONFIG_LSB|
				(mantissa&((1<<SLIC4_CONFIG_LSB)-1))
			);
		nbits=lgv-(SLIC4_CONFIG_MSB+SLIC4_CONFIG_LSB);
		bypass=val>>SLIC4_CONFIG_LSB&((1LL<<nbits)-1);
	}
	hu->token=token;
	hu->bypass=bypass;
	hu->nbits=nbits;
}
static int quantize(int x)
{
	int q;
	if(x<qlevels[8])
	{
		if(x<qlevels[4])
		{
			if(x<qlevels[2])
				q=x<qlevels[1]?x>=qlevels[0]:2;
			else
				q=x<qlevels[3]?3:4;
		}
		else
		{
			if(x<qlevels[6])
				q=x<qlevels[5]?5:6;
			else
				q=x<qlevels[7]?7:8;
		}
	}
	else
	{
		if(x<qlevels[12])
		{
			if(x<qlevels[10])
				q=x<qlevels[9]?9:10;
			else
				q=x<qlevels[11]?11:12;
		}
		else
		{
			if(x<qlevels[14])
				q=x<qlevels[13]?13:14;
			else
				q=x<qlevels[15]?15:16;
		}
	}
	return q;
}
static void slic4_pred(Image const *im, int kc, int kx, int ky, int scale, long long *sse, int *ret, long long **cell)
{
	int
		idx=im->iw*ky+kx,
		NN=ky>=2 ?im->data[(idx-im->iw*2)<<2|kc]:0,
		WW=kx>=2 ?im->data[(idx       -2)<<2|kc]:0,
		N =ky    ?im->data[(idx-im->iw  )<<2|kc]:0,
		W =kx    ?im->data[(idx       -1)<<2|kc]:0,
		NW=kx&&ky?im->data[(idx-im->iw-1)<<2|kc]:0;
	int vmax=MAXVAR(N, W), vmin=MINVAR(N, W), energy=vmax-vmin;
	int pred=N+W-NW;
	pred=CLAMP(vmin, pred, vmax);
	int
		hist_ctx=quantize(energy>>scale),
		sse_ctx=hist_ctx;
		//sse_ctx=nqlevels*quantize(abs(W-WW)>>shift)+quantize(abs(N-NN)>>shift),
		//sse_ctx=quantize(abs(W-WW)>>shift)*quantize(abs(N-NN)>>shift),
	*cell=sse+sse_ctx;
	int count=(int)(**cell&0xFFF);
	long long sum=**cell>>12;
	int corr=count?(int)(sum/count):0;
	pred+=corr;
	ret[0]=pred;
	ret[1]=hist_ctx;
}
static void slic4_update_sse(long long *cell, int error)
{
	int count=(int)(*cell&0xFFF);
	long long sum=*cell>>12;
	++count;
	sum+=error;
	if(count>640)
	{
		count>>=1;
		sum>>=1;
	}
	*cell=sum<<12|count;
}
Image const *guide=0;
int t46_encode(Image const *src, ArrayHandle *data, int loud)
{
	guide=src;
	double t_start=time_sec();
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int maxdepth=calc_maxdepth(src, 0), maxlevels=1<<maxdepth;
	double bpp=
		((double)src->src_depth[0]+src->src_depth[1]+src->src_depth[2]+src->src_depth[3])/
		(8*((src->src_depth[0]!=0)+(src->src_depth[1]!=0)+(src->src_depth[2]!=0)+(src->src_depth[3]!=0)));
	double usize=image_getBMPsize(src);
	int nch=src->nch;
	//int nch=get_nch32(src->data, res);//FIXME: differentiate between just gray and just alpha
	if(loud)
	{
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T46 SLIC4  Enc %s  CWHD %d*%d*%d*%d/8\n", g_buf, nch, src->iw, src->ih, maxdepth);
	}
	Image *im2=0, *im3=0, *im4=0;
	image_copy(&im2, src);
	image_copy(&im3, src);
	image_copy(&im4, src);
	int *hist=(int*)malloc(maxlevels*sizeof(int));
	const int nqlevels=_countof(qlevels)+1;
	long long *sse=(long long*)malloc(nqlevels*nqlevels*sizeof(long long));
	if(!im2||!im3||!im4||!hist||!sse)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	
	unsigned short *palettes[4]={0};
#ifdef SLIC4_DISABLE_PALETTE
	memset(header.palettesizes, 0, sizeof(header.palettesizes));
#else
	for(int kc=0;kc<nch;++kc)
	{
		int depth=im2->depth[kc], nlevels=1<<depth;
		memset(hist, 0, nlevels*sizeof(int));
		for(ptrdiff_t k=0;k<res;++k)
		{
			int val=im2->data[k<<2|kc]+(nlevels>>1);
			val=CLAMP(0, val, nlevels-1);
			++hist[val];
		}
		int palettesize=0;
		for(int k=0;k<nlevels;++k)
			palettesize+=hist[k]!=0;
		if(palettesize>((nlevels+15)>>4))
			header.palettesizes[kc]=0;
		else
		{
			header.palettesizes[kc]=palettesize;
			im2->depth[kc]=floor_log2_32(palettesize-1)+1;
			palettes[kc]=(unsigned short*)malloc(palettesize*sizeof(short));
			if(!palettes[kc])
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			for(int k=0, idx=0;k<nlevels;++k)
			{
				if(hist[k])
				{
					palettes[kc][idx]=k;
					++idx;
				}
			}
			int half=(palettesize+1)>>1;//ceil_half
			for(int k=0;k<res;++k)
			{
				unsigned short val=im2->data[k<<2|kc];
				int idx=0;
				int L=0, R=palettesize-1;
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
				im2->data[k<<2|kc]=idx-half;
			}
		}
	}
#endif
	free(hist);
	rct_JPEG2000_32(im2, 1);
	//memset(header.hist, 0, sizeof(header.hist));

	memset(g_hist, 0, sizeof(g_hist));
	for(int kc=0;kc<nch;++kc)//for each channel
	{
		if(header.palettesizes[kc]==1)
			continue;
		int shift=(MAXVAR(8, src->depth[kc])-8)>>2;
		memset(sse, 0, nqlevels*nqlevels*sizeof(long long));
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				//if(kc==0&&kx==114&&ky==2)//
				//if(kc==0&&kx==3&&ky==0)//
				//if(kc==0&&kx==4&&ky==0)//
				//if(kc==0&&kx==509&&ky==0)//
				//if(kc==0&&kx==0&&ky==0)//
				//if(kc==0&&kx==23&&ky==0)//
				//	printf("");

				int ctx[2];
				long long *cell;
				slic4_pred(im2, kc, kx, ky, shift, sse, ctx, &cell);
				int curr=im2->data[idx<<2|kc];
				int error=curr-ctx[0];
				TempHybrid hu;
				hybriduint_encode(error, &hu);
				if(hu.token>=SLIC4_HISTLEN)
					LOG_ERROR("Token OOB %d/%d", hu.token, SLIC4_HISTLEN);
				++g_hist[ctx[1]][hu.token];
				im3->data[idx<<2|kc]=hu.token<<16|ctx[1];
				im4->data[idx<<2|kc]=hu.bypass;
				slic4_update_sse(cell, error);
			}
		}
	}
	free(im2);
	if(loud)
	{
		printf("Histograms:\n");
		for(int kt=0;kt<_countof(qlevels)+1;++kt)
		{
			for(int ks=0;ks<SLIC4_HISTLEN;++ks)
				printf("%d ", g_hist[kt][ks]);
			printf("\n");
		}
		printf("SSE:\n");
		for(int k1=0;k1<_countof(qlevels)+1;++k1)
		{
			for(int k2=0;k2<_countof(qlevels)+1;++k2)
			{
				long long *cell=sse+17*k1+k2;
				int count=(int)(*cell&0xFFF);
				long long sum=*cell>>12;
				int corr=count?(int)(sum/count):0;
				printf("%3d ", corr);
			}
			printf("\n");
		}
	}
	free(sse);
	for(int kt=0;kt<_countof(qlevels)+1;++kt)
	{
		int weight=0, bins_used=0;
		for(int ks=0;ks<SLIC4_HISTLEN;++ks)
		{
			unsigned freq=g_hist[kt][ks];
			weight+=freq;
			bins_used+=freq!=0;
		}
		if(weight)
		{
			for(int ks=0, ks2=0, c=0;ks<SLIC4_HISTLEN;++ks)
			{
				unsigned freq=g_hist[kt][ks];
				header.hist[kt][ks]=g_hist[kt][ks]=(unsigned)((unsigned long long)c*(0x10000LL-bins_used)/weight)+ks2;
				ks2+=freq!=0;
				c+=freq;
			}
		}
		else
			memset(header.hist[kt], 0, sizeof(short[SLIC4_HISTLEN]));
		g_hist[kt][SLIC4_HISTLEN]=0x10000;
	}
	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, &header, sizeof(header));
	for(int kc=0;kc<4;++kc)//insert palettes, if any
	{
		if(palettes[kc])
		{
			dlist_push_back(&list, palettes[kc], header.palettesizes[kc]*sizeof(short));
			free(palettes[kc]);
			palettes[kc]=0;
		}
	}
	ANSCoder ec;
	ans_enc_init(&ec, &list);
	for(int kc=nch-1;kc>=0;--kc)//for each channel
	{
		if(header.palettesizes[kc]==1)
			continue;
		for(ptrdiff_t k=res-1;k>=0;--k)
		{
			//if((k<<2|kc)==6600)//
			//if((k<<2|kc)==4)//
			//if((k<<2|kc)==16)//
			//if((k<<2|kc)==2306)//
			//if((k<<2|kc)==0)//
			//if((k<<2|kc)==92)//
			//	printf("");

			int token=im3->data[k<<2|kc]>>16, ctx=im3->data[k<<2|kc]&0xFFFF, bypass=im4->data[k<<2|kc];
			if(token>=(1<<SLIC4_CONFIG_EXP))
			{
				int nbits=((token-(1<<SLIC4_CONFIG_EXP))>>(SLIC4_CONFIG_MSB+SLIC4_CONFIG_LSB))+SLIC4_CONFIG_EXP-(SLIC4_CONFIG_MSB+SLIC4_CONFIG_LSB);
				if(!nbits)
					LOG_ERROR("Bypass nbits 0");
				//unsigned char *ptr=(unsigned char*)&bypass;
				while(nbits>8)//encode bypass LSB-first
				{
					ans_enc(&ec, bypass&0xFF, 0, 256);
					nbits-=8;
					bypass>>=8;
				}
				//while(nbits>8)//encode bypass MSB-first
				//{
				//	nbits-=8;
				//	ans_enc(&ec, bypass>>nbits&0xFF, 0, 256);
				//}
				ans_enc(&ec, bypass&((1<<nbits)-1), 0, 1<<nbits);
			}
			ans_enc(&ec, token, g_hist[ctx], SLIC4_HISTLEN);
		}
	}
	ans_enc_flush(&ec);
	dlist_appendtoarray(&list, data);
	
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Encode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
		printf("csize %8d  CR %10.6lf\n", (int)list.nobj, usize/list.nobj);
	}
	dlist_clear(&list);
	free(im3);
	free(im4);
	return 1;
}
int t46_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	double t_start=time_sec();
	int res=dst->iw*dst->ih;
	if(srclen<sizeof(header))
		return 0;
	memcpy(&header, data, sizeof(header));
	data+=sizeof(header);
	srclen-=sizeof(header);

	unsigned short *palettes[4]={0};
	for(int kc=0;kc<dst->nch;++kc)
	{
		if(header.palettesizes[kc]>0)
		{
			int bytesize=header.palettesizes[kc]*sizeof(short);
			palettes[kc]=(unsigned short*)malloc(bytesize);
			if(!palettes[kc])
			{
				LOG_ERROR("Alloc error");
				return 0;
			}
			if(bytesize>srclen)
			{
				printf("Invalid file\n");
				return 0;
			}
			memcpy(palettes[kc], data, bytesize);
			data+=bytesize;
			srclen-=bytesize;
		}
	}
	const int nqlevels=_countof(qlevels)+1;
	long long *sse=(long long*)malloc(nqlevels*nqlevels*sizeof(long long));
	if(!sse)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int kt=0;kt<_countof(qlevels)+1;++kt)
	{
		for(int ks=0, overflow=0;ks<SLIC4_HISTLEN;++ks)
		{
			if(overflow)
				g_hist[kt][ks]=0x10000;//-SLIC4_HISTLEN+ks;
			else
			{
				unsigned cdf=header.hist[kt][ks];
				g_hist[kt][ks]=cdf;
				if(ks<SLIC4_HISTLEN-1)
					overflow|=cdf>header.hist[kt][ks+1];
			}
		}
		g_hist[kt][SLIC4_HISTLEN]=0x10000;
	}
	ANSCoder ec;
	ans_dec_init(&ec, data, data+srclen);
	for(int kc=0;kc<dst->nch;++kc)
	{
		if(header.palettesizes[kc]==1)
		{
			for(int k=0;k<res;++k)
				dst->data[k<<2|kc]=palettes[kc][0];
			continue;
		}
		int shift=(MAXVAR(8, dst->depth[kc])-8)>>2;
		memset(sse, 0, nqlevels*nqlevels*sizeof(long long));
		for(int ky=0, idx=0;ky<dst->ih;++ky)
		{
			for(int kx=0;kx<dst->iw;++kx, ++idx)
			{
				//if(kc==0&&kx==114&&ky==2)//
				//if(kc==0&&kx==1&&ky==0)//
				//if(kc==0&&kx==3&&ky==0)//
				//if(kc==0&&kx==4&&ky==0)//
				//if(kc==0&&kx==509&&ky==0)//
				//if(kc==0&&kx==0&&ky==0)//
				//if(kc==0&&kx==23&&ky==0)//
				//	printf("");

				int ctx[2];
				long long *cell;
				slic4_pred(dst, kc, kx, ky, shift, sse, ctx, &cell);
				unsigned *CDF=g_hist[ctx[1]];
				int sym=ans_dec(&ec, CDF, SLIC4_HISTLEN);
				if(sym>=(1<<SLIC4_CONFIG_EXP))
				{
					sym-=1<<SLIC4_CONFIG_EXP;
					int lsb=sym&((1<<SLIC4_CONFIG_LSB)-1);
					sym>>=SLIC4_CONFIG_LSB;
					int msb=sym&((1<<SLIC4_CONFIG_MSB)-1);
					sym>>=SLIC4_CONFIG_MSB;
					int nbits=sym+SLIC4_CONFIG_EXP-(SLIC4_CONFIG_MSB+SLIC4_CONFIG_LSB);
					int bypass=nbits&7?ans_dec(&ec, 0, 1<<(nbits&7)):0;
					for(int k=0, nb=nbits>>3;k<nb;++k)
						bypass=bypass<<8|ans_dec(&ec, 0, 256);
					sym=1;
					sym<<=SLIC4_CONFIG_MSB;
					sym|=msb;
					sym<<=nbits;
					sym|=bypass;
					sym<<=SLIC4_CONFIG_LSB;
					sym|=lsb;
				}
				int error=(sym>>1)^-(sym&1), curr=error+ctx[0];
				//if(guide&&curr!=guide->data[idx<<2|kc])
				//	LOG_ERROR("Guide error");
				dst->data[idx<<2|kc]=curr;
				slic4_update_sse(cell, error);
			}
		}
	}
	free(sse);
	rct_JPEG2000_32(dst, 0);
	for(int kc=0;kc<dst->nch;++kc)
	{
		if(palettes[kc])
		{
			for(int k=0;k<res;++k)
			{
				int val=dst->data[k<<2|kc];
				val+=(header.palettesizes[kc]+1)>>1;
				if((unsigned)val>=header.palettesizes[kc])
				{
					LOG_ERROR("Palette error");
					return 0;
				}
				dst->data[k<<2|kc]=palettes[kc][val];
			}
			free(palettes[kc]);
		}
	}
	if(loud)
	{
		printf("\n");//skip progress line
		printf("Decode elapsed ");
		timedelta2str(0, 0, time_sec()-t_start);
		printf("\n");
	}
	return 1;
}