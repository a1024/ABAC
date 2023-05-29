#include"battle.h"
#include"lodepng.h"
#include<stdio.h>//for debugging
#include<stdlib.h>
#include<string.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#include<math.h>
#include"rans_common.h"
static const char file[]=__FILE__;

	#define PROF(...)//

const int tag_arans='A'|'N'<<8|'0'<<16|'5'<<24;

typedef struct SortedHistInfoStruct
{
	int	sym,  //symbol
		freq, //original freq
		qfreq;//quantized freq
} SortedHistInfo;
int histinfo_byfreq(const void *left, const void *right)
{
	SortedHistInfo const *a, *b;

	a=(SortedHistInfo const*)left;
	b=(SortedHistInfo const*)right;
	return (a->freq>b->freq)-(a->freq<b->freq);
}
int histinfo_bysym(const void *left, const void *right)
{
	SortedHistInfo const *a, *b;

	a=(SortedHistInfo const*)left;
	b=(SortedHistInfo const*)right;
	return (a->sym>b->sym)-(a->sym<b->sym);
}
int rans_calc_histogram(const unsigned char *buffer, int nsymbols, int bytestride, unsigned char *hist, int prob_bits, int integrate)//hist is unsigned char due to alignment issues, but it's 16bit
{
	int prob_sum=1<<prob_bits;
	//MY_ASSERT(ANS_NLEVELS<prob_sum, "Channel depth %d >= PROB_BITS %d", ANS_NLEVELS, prob_sum);//what if ANS_NLEVELS = 2^N-1 ?
	if(!nsymbols)
	{
		memset(hist, 0, ANS_NLEVELS*sizeof(short));//16bit
		LOG_ERROR("Symbol count is zero");
	}
	SortedHistInfo h[ANS_NLEVELS];
	for(int k=0;k<ANS_NLEVELS;++k)
	{
		h[k].sym=k;
		h[k].freq=0;
	}
	int bytesize=nsymbols*bytestride;
	PROF(HISTOGRAM_INIT);
	for(int k=0;k<bytesize;k+=bytestride)//this loop takes 73% of encode time
		++h[buffer[k]].freq;
	PROF(HISTOGRAM_LOOKUP);
	for(int k=0;k<ANS_NLEVELS;++k)
		h[k].qfreq=((long long)h[k].freq<<ANS_PROB_BITS)/nsymbols;

	//print_histogram(h, nsymbols);//
	
	if(nsymbols!=prob_sum)
	{
		const int prob_max=prob_sum-1;

		isort(h, ANS_NLEVELS, sizeof(SortedHistInfo), histinfo_byfreq);
		int idx=0;
		for(;idx<ANS_NLEVELS&&!h[idx].freq;++idx);
		for(;idx<ANS_NLEVELS&&!h[idx].qfreq;++idx)
			++h[idx].qfreq;
		for(idx=ANS_NLEVELS-1;idx>=0&&h[idx].qfreq>=prob_max;--idx);
		for(++idx;idx<ANS_NLEVELS;++idx)
			h[idx].qfreq=prob_max;

		int error=-prob_sum;//too much -> +ve error & vice versa
		for(int k=0;k<ANS_NLEVELS;++k)
			error+=h[k].qfreq;
		if(error>0)
		{
			while(error)
			{
				for(int k=0;k<ANS_NLEVELS&&error;++k)
				{
					int dec=h[k].qfreq>1;
					h[k].qfreq-=dec, error-=dec;
				}
			}
		}
		else
		{
			while(error)
			{
				for(int k=ANS_NLEVELS-1;k>=0&&error;--k)
				{
					int inc=h[k].qfreq<prob_max;
					h[k].qfreq+=inc, error+=inc;
				}
			}
		}
		if(error)
			LOG_ERROR("Internal error: histogram adds up to %d != %d", prob_sum+error, prob_sum);
		isort(h, ANS_NLEVELS, sizeof(SortedHistInfo), histinfo_bysym);
	}
	int sum=0;
	for(int k=0;k<ANS_NLEVELS;++k)
	{
		if(h[k].qfreq>0xFFFF)
			LOG_ERROR("Internal error: symbol %d has probability %d", k, h[k].qfreq);
		memcpy(hist+((size_t)k<<1), integrate?&sum:&h[k].qfreq, 2);//2-byte alignment
		sum+=h[k].qfreq;
	}
	if(sum!=ANS_L)
		LOG_ERROR("Internal error: CDF ends with 0x%08X, should end with 0x%08X", sum, ANS_L);
	return 1;
}
#if 0
void print_histogram(SymbolInfo *hist, int nsymbols)
{
	printf("s\tf\tCDF,\timsize %d\n", nsymbols);
	for(int k=0;k<nsymbols;++k)
	{
		SymbolInfo *sym=hist+k;
		if(sym->freq)
		{
			if(!sym->qfreq)
				printf("[%3d] s %02X q %04X f %04X UNDERFLOW\n", k, sym->idx, sym->qfreq, sym->freq);
			else if(sym->qfreq==0xFFFF||sym->qfreq==0x10000)
				printf("[%3d] s %02X q %04X f %04X OVERFLOW\n", k, sym->idx, sym->qfreq, sym->freq);
			else
				printf("[%3d] s %02X q %04X f %04X\n", k, sym->idx, sym->qfreq, sym->freq);
		}
	}
}
#endif

int rans_prep(const void *hist_ptr, int bytespersymbol, SymbolInfo **info, unsigned char **CDF2sym, int loud)
{
	int tempsize=bytespersymbol*(ANS_NLEVELS*sizeof(SymbolInfo)+(ANS_L&-(CDF2sym!=0)));
	*info=(SymbolInfo*)malloc(tempsize);
	if(!*info)
		LOG_ERROR("Failed to allocate temp buffer");
	if(CDF2sym)
		*CDF2sym=(unsigned char*)*info+bytespersymbol*ANS_NLEVELS*sizeof(SymbolInfo);
	for(int kc=0;kc<bytespersymbol;++kc)
	{
		const unsigned short *c_histogram=(const unsigned short*)hist_ptr+(kc<<ANS_DEPTH);
		SymbolInfo *c_info=*info+(kc<<ANS_DEPTH);
		unsigned char *c_CDF2sym=CDF2sym?*CDF2sym+(kc<<ANS_PROB_BITS):0;
		int sum=0;
		for(int sym=0;sym<ANS_NLEVELS;++sym)
		{
			SymbolInfo *p=c_info+sym;
			memcpy(&p->freq, c_histogram+sym, 2);//alignment
			p->neg_freq=-p->freq;
			p->CDF=sum;
			p->reserved0=0;

			if(p->freq<2)//0 freq: don't care, 1 freq:		//Ryg's fast rANS encoder
			{
			//	p->shift=0;
			//	p->inv_freq=0xFFFFFFFF;
				p->bias=sum+ANS_L-1;
			}
			else
			{
			//	p->shift=ceil_log2(p->freq)-1;
			//	p->inv_freq=(unsigned)(((0x100000000<<p->shift)+p->freq-1)/p->freq);
				p->bias=sum;
			}
			if(p->freq)
				p->invf=(0x0001000000000000+p->freq-1)/p->freq;
			else
				p->invf=0;

			p->renorm_limit=p->freq<<(32-ANS_PROB_BITS);

			if(CDF2sym&&sym)
			{
				for(int k2=c_info[sym-1].CDF;k2<(int)p->CDF;++k2)
					c_CDF2sym[k2]=sym-1;
			}
			sum+=p->freq;
		}
		if(CDF2sym)
		{
			for(int k2=c_info[ANS_NLEVELS-1].CDF;k2<ANS_L;++k2)
				c_CDF2sym[k2]=ANS_NLEVELS-1;
		}
		if(sum!=ANS_L)
			LOG_ERROR("histogram sum = %d != %d", sum, ANS_L);
		if(loud)
		{
#ifdef ANS_PRINT_HISTOGRAM
			static int printed=0;
			if(printed<1)
			{
				printf("s\tf\tCDF\n");
				for(int k=0;k<ANS_NLEVELS;++k)
				{
					auto &si=c_info[k];
					if(c_histogram[k])
						printf("%3d\t%5d = %04X\t%04X\n", k, c_histogram[k], c_histogram[k], si.CDF);
				}
				++printed;
			}
#endif
		}
	}
	//	if(!calc_hist_derivaties((const unsigned short*)hist_ptr+kc*ANS_NLEVELS, info+(kc<<ANS_DEPTH), CDF2sym+ANS_L*kc, loud))
	//		return false;
	return 1;
}

static unsigned short g_hist[256], g_freq[256];
int arans_encode(const void *src, ptrdiff_t nbytes, int bytestride, ArrayHandle *out)
{
	DList list;
	const unsigned char *buf=(const unsigned char*)src;
	
	long long cycles=__rdtsc();
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, 12);//header = {tag, csize}

	unsigned state=0x10000;
	unsigned char hist_idx=0;
	memset(g_hist, -1, sizeof(g_hist));
	for(int k=0;k<256;++k)
		g_freq[k]=256;
	for(ptrdiff_t ks=0;ks<nbytes;ks+=bytestride)
	{
		unsigned char sym=buf[ks];
		unsigned short CDF=0;
		for(int k=0;k<sym;++k)
			CDF+=g_freq[k];
		unsigned short freq=g_freq[sym];
		if(state>=(unsigned)(freq<<16))
		{
			dlist_push_back(&list, &state, 2);
			state>>=16;
		}
		state=state/freq<<16|(state%freq+CDF);

		//update
		++hist_idx;
		if(g_hist[hist_idx]!=0xFFFF&&g_freq[g_hist[hist_idx]]>0)
		{
			g_freq[g_hist[hist_idx]]-=4;
			g_freq[sym]+=4;
		}
		g_hist[hist_idx]=sym;
	}

	size_t dstidx;
	if(*out)
	{
		if(out[0]->esize!=1)
			return RANS_INVALID_DST;
		dstidx=out[0]->count;
		ARRAY_APPEND(*out, 0, 0, 1, 0);
	}
	else
	{
		dstidx=0;
		ARRAY_ALLOC(char, *out, 0, 0, 0, 0);
	}
	dlist_appendtoarray(&list, out);
	memcpy(out[0]->data+dstidx, &tag_arans, 4);
	memcpy(out[0]->data+dstidx+4, &list.nobj, 8);

	cycles=__rdtsc()-cycles;

	dlist_clear(&list);
	return RANS_SUCCESS;
}

static int decode_error(size_t p, size_t srcstart, ptrdiff_t ks)
{
	printf("Buffer underrun p<srcstart, %p < %p at %d\n", (void*)p, (void*)srcstart, (int)ks);
	return 0;
}
#define READ_GUARD(NBYTES, IDX)		if((srcptr-=NBYTES)<srcstart)return decode_error((size_t)srcptr, (size_t)srcstart, IDX)
int arans_decode(const unsigned char *srcdata, ptrdiff_t srclen, ptrdiff_t nbytes, int symbytes, int is_signed, void *dstbuf, int loud)
{
	const int
		histsize=ANS_NLEVELS*sizeof(short), lghistsize=9,
		infosize=ANS_NLEVELS*sizeof(SymbolInfo), lginfosize=13;
	const unsigned char
		*data=(const unsigned char*)srcdata,
		*srcptr=data,
		*srcstart,
		*hist;
	unsigned char *pixels=(unsigned char*)dstbuf;
	SymbolInfo *info;
	unsigned char *CDF2sym;
	int internalheadersize=4+8+symbytes*ANS_NLEVELS;
	int chmask=symbytes-1;
	size_t csize;

	if(srclen<internalheadersize)
		return RANS_BUFFER_OVERRUN;
	if(memcmp(srcptr, &tag_arans, 4))
		return RANS_INVALID_TAG;
	srcptr+=4;

	memcpy(&csize, srcptr, 8);
	srcptr+=8;

	hist=srcptr;
	rans_prep(hist, symbytes, &info, &CDF2sym, 0);
	srcstart=srcptr+(symbytes<<9);
	srcptr=srcstart+csize;

	long long t1=__rdtsc();

	unsigned state[16]={0};
	READ_GUARD(symbytes*sizeof(unsigned), 0);
	memcpy(state, srcptr, symbytes*sizeof(unsigned));
	for(ptrdiff_t ks=0, idx=0;ks<nbytes;ks+=symbytes)
	{
		for(int kc=0;kc<symbytes;++kc, ++idx)
		{
			unsigned short c=(unsigned short)state[kc];
			unsigned char val=CDF2sym[kc<<ANS_PROB_BITS|c];
			SymbolInfo *p=info+(kc<<ANS_DEPTH|val);
			//if(!p->freq)
			//	LOG_ERROR("Symbol 0x%02X has zero frequency", s);
			
			if(is_signed)
			{
				int neg=val&1;
				val>>=1;
				val^=-neg;
				val+=neg;
				val|=neg<<7&-!val;
			}
			pixels[idx]=val;
			PROF(FETCH);

#ifdef ANS_PRINT_STATE2
			printf("dec: 0x%08X = 0x%04X*(0x%08X>>%d)+0x%04X-0x%08X\n", si.freq*(state>>ANS_PROB_BITS)+c-si.CDF, (int)si.freq, state, ANS_PROB_BITS, c, si.CDF);
#endif
			state[kc]=p->freq*(state[kc]>>ANS_PROB_BITS)+c-p->CDF;
			PROF(UPDATE);

			if(state[kc]<ANS_L)
			{
				if(idx>=nbytes-1)//shouldn't need this
					break;//
				READ_GUARD(2, idx);
				state[kc]<<=16;
				memcpy(state+kc, srcptr, 2);
			}
			PROF(RENORM);
		}
	}
	free(info);
	return RANS_SUCCESS;
}


//	#define T21_DISABLE_TRANSFORMS
//	#define T21_DISABLE_SHUFFLE

//test 21:

//void normalize_histogram(unsigned *srchist, int nlevels, int nsymbols, unsigned short *CDF);
//void t16_prepblock(const unsigned char *b2, const unsigned short *CDF, int bw, int bh, int kc, int bx, int by, int alpha, int blockw, int blockh, int margin, unsigned *CDF2, int *xend, int *yend);
double t21_calcentropy_u16(unsigned short *hist, int histsum)
{
	double entropy=0;
	for(int k=0;k<256;++k)
	{
		int freq=hist[k];
		if(freq)
		{
			double p=(double)freq/histsum;
			entropy-=p*log2(p);
		}
	}
	return entropy;
}
void t21_normalize_histogram(unsigned *srchist, int nlevels, int nsymbols, unsigned short *CDF)//hist is unsigned char due to alignment issues, but it's 16bit
{
	SortedHistInfo h[512];
	if(!nsymbols)//bypass
	{
		for(int k=0;k<nlevels;++k)
			CDF[k]=(unsigned short)(k<<8);
		return;
	}
	for(int k=0;k<nlevels;++k)
	{
		h[k].sym=k;
		h[k].freq=srchist[k];
	}
	for(int k=0;k<nlevels;++k)
		h[k].qfreq=((long long)h[k].freq<<16)/nsymbols;
#if 0
	if(nsymbols!=0x10000)
	{
		const int prob_max=0x10000-(nlevels-1);

		isort(h, nlevels, sizeof(SortedHistInfo), histinfo_byfreq);
		int idx=0;
		for(;idx<nlevels&&!h[idx].freq;++idx);
		for(;idx<nlevels&&!h[idx].qfreq;++idx)
			++h[idx].qfreq;
		for(idx=nlevels-1;idx>=0&&h[idx].qfreq>=prob_max;--idx);
		for(++idx;idx<nlevels;++idx)
			h[idx].qfreq=prob_max;

		int error=-0x10000;//too much -> +ve error & vice versa
		for(int k=0;k<nlevels;++k)
			error+=h[k].qfreq;
		if(error>0)
		{
			while(error)
			{
				for(int k=0;k<nlevels&&error;++k)
				{
					int dec=h[k].qfreq>1;
					h[k].qfreq-=dec, error-=dec;
				}
			}
		}
		else
		{
			while(error)
			{
				for(int k=nlevels-1;k>=0&&error;--k)
				{
					int inc=h[k].qfreq<prob_max;
					h[k].qfreq+=inc, error+=inc;
				}
			}
		}
		isort(h, nlevels, sizeof(SortedHistInfo), histinfo_bysym);
	}
#endif
	int sum=0;
	for(int k=0;k<nlevels;++k)
	{
		CDF[k]=sum;
		sum+=h[k].qfreq;
	}
}
#if 0
typedef struct BlockIdxStruct
{
	unsigned char x, y;
} BlockIdx;
typedef struct SwapInfoStruct
{
	BlockIdx b1, b2;
} SwapInfo;
#endif
void t21_calchist_u16(const unsigned char *buf2, int bw, int kc, int x1, int x2, int y1, int y2, unsigned short *hist)
{
	memset(hist, 0, 256*sizeof(short));
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			unsigned char sym=buf2[(bw*ky+kx)<<2|kc];
			++hist[sym];
		}
	}
}
void t21_calchist_u32(const unsigned char *buf2, int bw, int kc, int x1, int x2, int y1, int y2, unsigned *hist)
{
	memset(hist, 0, 256*sizeof(int));
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			unsigned char sym=buf2[(bw*ky+kx)<<2|kc];
			++hist[sym];
		}
	}
}
int t21_blockdist(unsigned short *h1, unsigned short *h2)
{
	int sum=0;
	for(int k=0;k<256;++k)
	{
		int dist=abs(h1[k]-h2[k]);
		sum+=dist;
	}
	return sum;
}
void t21_blockswap(unsigned char *buf2, int bw, int b1, int b2, int blocksize)
{
	if(b1==b2)
		return;
	int cx=bw/blocksize;
	int x1=b1%cx*blocksize, y1=b1/cx*blocksize, x2=b2%cx*blocksize, y2=b2/cx*blocksize;
	for(int ky=0;ky<blocksize;++ky)
	{
		for(int kx=0;kx<blocksize;++kx)
		{
			int idx1=bw*(y1+ky)+x1+kx, idx2=bw*(y2+ky)+x2+kx;

			//unsigned char temp;
			//SWAPVAR(buf2[idx1<<2  ], buf2[idx2<<2  ], temp);
			//SWAPVAR(buf2[idx1<<2|1], buf2[idx2<<2|1], temp);
			//SWAPVAR(buf2[idx1<<2|2], buf2[idx2<<2|2], temp);
			//buf2[idx1<<2|1]+=16;//
			//buf2[idx2<<2|1]+=16;//

			int temp=((int*)buf2)[idx1];
			((int*)buf2)[idx1]=((int*)buf2)[idx2];
			((int*)buf2)[idx2]=temp;
		}
	}
}

static double *g_entropy=0;
int t21_cmp_entropy(const void *p1, const void *p2)
{
	unsigned short idx1=*(const unsigned short*)p1, idx2=*(const unsigned short*)p2;

	return (g_entropy[idx1]<g_entropy[idx2])-(g_entropy[idx1]>g_entropy[idx2]);//descending order (hard->easy)

	//return (g_entropy[idx1]>g_entropy[idx2])-(g_entropy[idx1]<g_entropy[idx2]);//ascending order (easy->hard)
}
ArrayHandle t21_getblockorder(unsigned char *buf2, int bw, int bh, int blocksize)
{
	int cx=bw/blocksize, cy=bh/blocksize, nblocks=cx*cy;
	unsigned short *hist=(unsigned short*)malloc(nblocks*256LL*sizeof(short));
	char *visited=(char*)malloc(nblocks);
	double *entropy=(double*)malloc(nblocks*sizeof(double));
	if(!hist||!visited||!entropy)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	for(int by=0;by<cy;++by)
	{
		int ky=by*blocksize;
		for(int bx=0;bx<cx;++bx)
		{
			int kx=bx*blocksize;
			int blockidx=cx*by+bx;
			t21_calchist_u16(buf2, bw, 1, kx, kx+blocksize, ky, ky+blocksize, hist+blockidx);
			entropy[blockidx]=t21_calcentropy_u16(hist+blockidx, blocksize*blocksize);
		}
	}
	ArrayHandle order;
	ARRAY_ALLOC(unsigned short, order, 0, nblocks, 0, 0);
	memset(visited, 0, nblocks);
	
	unsigned short *ptr=(unsigned short*)array_at(&order, 0);
	
	for(int kb=0;kb<order->count;++kb)
		ptr[kb]=kb;

	g_entropy=entropy;
	isort(ptr, order->count, order->esize, t21_cmp_entropy);
	g_entropy=0;

#if 0
	int b1=0, b2;
	*ptr=b1;
	++ptr;
	visited[b1]=1;
	for(;;)
	{
		int done=0;
		b2=-1;
		for(int kb=0;kb<nblocks;++kb)
		{
			if(visited[kb])
				++done;
			else if(b2==-1)
				b2=kb;
		}
		if(done==nblocks)//visited all
			break;
		
		int dist=t21_blockdist(hist+b1, hist+b2);
		for(int kb=0;kb<nblocks;++kb)
		{
			if(!visited[kb]&&kb!=b2)
			{
				int d2=t21_blockdist(hist+b1, hist+kb);
				if(dist>d2)
				{
					dist=d2;
					b2=kb;
				}
			}
		}
		*ptr=b2;
		++ptr;
		visited[b2]=1;
		b1=b2;
	}
#endif

#if 0
	for(int b1=0;b1<nblocks-1;++b1)
	{
		int b2=b1+1;
		int d12=t21_blockdist(hist+b1, hist+b2);
		for(int b3=0;b3<nblocks;++b3)
		{
			if(b3!=b1&&b3!=b2&&!visited[b3])
			{
				int d13=t21_blockdist(hist+b1, hist+b3);
				if(d12>d13)
				{
					d12=d13;
					b2=b3;
				}
			}
		}
		ptr=(unsigned short*)array_at(&order, b1+1);
		*ptr=b2;
		visited[b2]=1;
	}
#endif

	free(hist);
	free(visited);
	free(entropy);
	return order;
}
#if 0
ArrayHandle t21_sortblocks(unsigned char *buf2, int bw, int bh, int blocksize)
{
	ArrayHandle swaps;
	int cx=bw/blocksize, cy=bh/blocksize, nblocks=cx*cy;
	unsigned short *hist=(unsigned short*)malloc(nblocks*256LL*sizeof(short));
	//int *distances=(int*)malloc(nblocks*sizeof(int));
	//BlockIdx *order=(BlockIdx*)malloc(nblocks*sizeof(BlockIdx));
	if(!hist)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	for(int by=0;by<cy;++by)
	{
		int ky=by*blocksize;
		for(int bx=0;bx<cx;++bx)
		{
			int kx=bx*blocksize;
			int blockidx=cx*by+bx;
			t21_calchist(buf2, bw, 1, kx, kx+blocksize, ky, ky+blocksize, hist+blockidx);
		}
	}
	ARRAY_ALLOC(SwapInfo, swaps, 0, 0, 0, 0);
	
	//for(int b1=0;b1<nblocks-1;++b1)
	for(int by=0;by<cy-1;++by)
	{
		for(int bx=0;bx<cx-1;++bx)
		{
			int b1=cx*by+bx;
			int swapidx1=b1+1, swapidx2=b1+cx;
			int dright=t21_blockdist(hist+b1, hist+swapidx1);
			int dbottom=t21_blockdist(hist+b1, hist+swapidx2);
			for(int b2=swapidx1+1;b2<nblocks;++b2)
			{
				if(b2!=swapidx2)
				{
					int d2=t21_blockdist(hist+b1, hist+b2);
					if(dright>d2)
					{
						dright=d2;
						swapidx1=b2;
					}
					if(dbottom>d2)
					{
						dbottom=d2;
						swapidx2=b2;
					}
				}
			}
#ifndef T21_DISABLE_SHUFFLE
			int idx0=b1+1;
			if(swapidx1!=idx0)
			{
				int x1=idx0%cx, y1=idx0/cx, x2=swapidx1%cx, y2=swapidx1/cx;
				//if(swaps->count==42)
				//	printf("");
				SwapInfo *ptr=(SwapInfo*)ARRAY_APPEND(swaps, 0, 1, 1, 0);
				ptr->b1.x=x1;
				ptr->b1.y=y1;
				ptr->b2.x=x2;
				ptr->b2.y=y2;

				//printf("(%d %d) -> (%d %d)\n", x1, y1, x2, y2);

				t21_blockswap(buf2, bw, idx0, swapidx1, blocksize);
			}
			idx0=b1+cx;
			if(swapidx2!=idx0)
			{
				int x1=idx0%cx, y1=idx0/cx, x2=swapidx2%cx, y2=swapidx2/cx;
				//if(swaps->count==42)
				//	printf("");
				SwapInfo *ptr=(SwapInfo*)ARRAY_APPEND(swaps, 0, 1, 1, 0);
				ptr->b1.x=x1;
				ptr->b1.y=y1;
				ptr->b2.x=x2;
				ptr->b2.y=y2;

				//printf("(%d %d) -> (%d %d)\n", x1, y1, x2, y2);

				t21_blockswap(buf2, bw, idx0, swapidx2, blocksize);
			}
#endif
		}
	}

	free(hist);
	//free(distances);
	return swaps;
}
#endif
void t21_prepblock(const unsigned char *buf2, const unsigned short *CDF0, int bw, int bh, int kc, int bx1, int by1, int bx2, int by2, int alpha, int blocksize, unsigned *CDF2, int *xend, int *yend)
{
	int overflow=0;//CDF0 overflow can happen only once
	if(bx1!=-1&&by1!=-1)
	{
		int x1=bx1*blocksize, x2=x1+blocksize,
			y1=by1*blocksize, y2=y1+blocksize;
		if(x2>bw)
			x2=bw;
		if(y2>bh)
			y2=bh;
		//memset(CDF2, 0, 256*sizeof(unsigned));
		t21_calchist_u32(buf2, bw, kc, x1, x2, y1, y2, CDF2);
		int sum=0, count=(x2-x1)*(y2-y1),
			cdf1, f1, f2, freq;
		for(int sym=0;sym<256;++sym)
		{
			cdf1=!overflow?CDF0[sym]:0x10000;
			if(sym<255)
				overflow|=cdf1>CDF0[sym+1];
			f1=(sym<255&&!overflow?CDF0[sym+1]:0x10000)-cdf1;

			f2=(int)(((long long)CDF2[sym]<<16)/count);//normalize

			freq=f1+(int)(((long long)f2-f1)*alpha>>16);//blend

			freq=((unsigned)(freq*0xFF00)>>16)+1;//guard
			if(freq<0||freq>0xFF01)
				LOG_ERROR("Impossible freq 0x%04X / 0x10000", freq);
			CDF2[sym]=sum;
			sum+=freq;
			if(sum>0x10000)
				LOG_ERROR("ANS CDF sum 0x%04X, freq 0x%04X", sum, freq);
		}
		CDF2[256]=0x10000;
	}
	else
	{
		for(int sym=0;sym<256;++sym)
		{
			if(overflow)
				CDF2[sym]=0xFF00|sym;
			else
			{
				int cdf=CDF0[sym];
				CDF2[sym]=((unsigned)(cdf*0xFF00)>>16)+sym;
				if(sym<255)
					overflow|=cdf>CDF0[sym+1];
			}
		}
		CDF2[256]=0x10000;
	}
	*yend=(by2+1)*blocksize;
	*xend=(bx2+1)*blocksize;
	if(*yend>bh)
		*yend=bh;
	if(*xend>bw)
		*xend=bw;
}
double t21_encblock(const unsigned char *buf2, int bw, int kc, int x1, int x2, int y1, int y2, const unsigned *CDF2, unsigned *state, DList *list)
{
	double csize=0;
	for(int ky=y2-1;ky>=y1;--ky)
	{
		for(int kx=x2-1;kx>=x1;--kx)
		{
			unsigned char sym=buf2[(bw*ky+kx)<<2|kc];

			int cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;

			//if(kc==0&&ky==0&&kx==4)//
			//	printf("CXY %d %d %d  sym 0x%02X cdf 0x%04X freq 0x%04X  state 0x%08X\n", kc, kx, ky, sym, cdf, freq, *state);

			if(!freq)
				LOG_ERROR("ZPS");

			double p=freq/65536.;
			csize-=log2(p);
						
			if(*state>=(unsigned)(freq<<16))//renorm
			{
				dlist_push_back(list, state, 2);
				*state>>=16;
			}
			debug_enc_update(*state, cdf, freq, kx, ky, 0, kc, sym);
			*state=*state/freq<<16|(cdf+*state%freq);//update
		}
	}
	csize/=8;
	return csize;
}
void t21_decblock(int bw, int kc, int x1, int x2, int y1, int y2, unsigned *state, const unsigned *CDF2, unsigned char *buf, const unsigned char **srcptr, const unsigned char *srcstart)
{
	for(int ky2=y1;ky2<y2;++ky2)
	{
		for(int kx2=x1;kx2<x2;++kx2)//for each pixel
		{
			unsigned c=(unsigned short)*state;
			int sym=0;

			//if(kc==0&&ky2==0&&kx2==4)//
			//	kx2=4;

			int L=0, R=256, found=0;
			while(L<=R)
			{
				sym=(L+R)>>1;
				if(CDF2[sym]<c)
					L=sym+1;
				else if(CDF2[sym]>c)
					R=sym-1;
				else
				{
					found=1;
					break;
				}
			}
			if(!found)
				sym=L+(L<256&&CDF2[L]<c)-1;
			else
				for(;sym<256-1&&CDF2[sym+1]==c;++sym);

			//if(sym!=debug_ptr[(bw*ky2+kx2)<<2|kc])//
			//	LOG_ERROR("");

			buf[(bw*ky2+kx2)<<2|kc]=(unsigned char)sym;

			unsigned cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;

			debug_dec_update(*state, cdf, freq, kx2, ky2, 0, kc, sym);
			*state=freq*(*state>>16)+c-cdf;//update
			if(*state<0x10000)//renorm
			{
				*state<<=16;
				if(*srcptr-2>=srcstart)
				{
					*srcptr-=2;
					memcpy(state, *srcptr, 2);
				}
			}
		}
	}
}
size_t test21_encode(const unsigned char *src, int bw, int bh, int alpha, int blocksize, ArrayHandle *data, int loud)
{
	int res=bw*bh;
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	unsigned *hist=(unsigned*)malloc(768LL*sizeof(unsigned));
	unsigned short *CDF=(unsigned short*)malloc(768LL*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!buf2||!hist||!CDF||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
#ifndef T21_DISABLE_TRANSFORMS
	apply_transforms_fwd(buf2, bw, bh);
#endif
	
	int cx=(bw+blocksize-1)/blocksize,
		cy=(bh+blocksize-1)/blocksize;
	ArrayHandle order=t21_getblockorder(buf2, bw, bh, blocksize);
	//ArrayHandle swaps=t21_sortblocks(buf2, bw, bh, blocksize);
	
	//lodepng_encode_file("t21.PNG", buf2, bw, bh, LCT_RGBA, 8);//
	
#if 1
	short *porder=(short*)order->data;
	int repeated_indices=0;
	for(int k=0;k<order->count-1;++k)
	{
		for(int k2=k+1;k2<order->count;++k2)
		{
			if(porder[k]==porder[k2])
			{
				printf("%d %d\n", k, k2);
				++repeated_indices;
			}
		}
	}
	printf("%d repeated indices\n\n", repeated_indices);
	
	//for(int ky=0;ky<cy;++ky)
	//{
	//	for(int kx=0;kx<cx;++kx)
	//		printf(" %3d", porder[cx*ky+kx]);
	//	printf("\n");
	//}
	//printf("\n");

	//for(int k=0;k<order->count;++k)
	//	printf("%3d %3d\n", k, porder[k]);
#if 0
	for(int kx=0;kx<cx;++kx)
		printf("     %2d", kx);
	printf("\n\n");
	for(int ky=0;ky<cy;++ky)
	{
		for(int kx=0;kx<cx;++kx)
		{
			short idx=*(short*)array_at(&order, cx*ky+kx);

			if(idx==1)//
				idx=1;

			int bx=idx%cx, by=idx/cx;
			printf("  %2d %2d", bx, by);
		}
		printf("   %2d\n", ky);
	}
	printf("\n");
#endif
#endif
	
	memset(hist, 0, 768LL*sizeof(unsigned));
	for(int kc=0;kc<3;++kc)
	{
		for(int k=0;k<res;++k)
		{
			unsigned char sym=buf2[k<<2|kc];
			++hist[kc<<8|sym];
		}
	}
	t21_normalize_histogram(hist, 256, res, CDF);//this is just to pack the histogram, CDF is renormalized again with ramp guard
	t21_normalize_histogram(hist+256, 256, res, CDF+256);
	t21_normalize_histogram(hist+512, 256, res, CDF+512);
	
	DList list;
	dlist_init(&list, 1, 1024, 0);

	int ansbookmarks[3]={0};
	dlist_push_back(&list, 0, 12);
	dlist_push_back(&list, &order->count, 4);
	dlist_push_back(&list, order->data, order->count*order->esize);
	//dlist_push_back(&list, &swaps->count, 4);
	//dlist_push_back(&list, swaps->data, swaps->count*swaps->esize);
	dlist_push_back(&list, CDF, 768*sizeof(short));
	
	for(int kc=0;kc<3;++kc)
	{
		unsigned state=0x10000;
		
		double csize=0;
		int blockcount=0;

		int bx2, by2, bx1, by1;
		int xend=0, yend=0;
		if(cy*blocksize<bh)
		{
			if(cx*blocksize<bw)
			{
				bx1=(cx-1)*blocksize, by1=(cy-1)*blocksize, bx2=cx*blocksize, by2=cy*blocksize;
				t21_prepblock(buf2, CDF+((size_t)kc<<8), bw, bh, kc, bx1, by1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
				t21_encblock(buf2, bw, kc, bx2*blocksize, xend, by2*blocksize, yend, CDF2, &state, &list);
			}
			for(int kb=cx-1;kb>0;--kb)
			{
				bx1=(kb-1)*blocksize, by1=cy*blocksize, bx2=kb*blocksize, by2=cy*blocksize;
				t21_prepblock(buf2, CDF+((size_t)kc<<8), bw, bh, kc, bx1, by1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
				t21_encblock(buf2, bw, kc, bx2*blocksize, xend, by2*blocksize, yend, CDF2, &state, &list);
			}
			bx2=0, by2=cy*blocksize;
			t21_prepblock(buf2, CDF+((size_t)kc<<8), bw, bh, kc, -1, -1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
			t21_encblock(buf2, bw, kc, bx2*blocksize, xend, by2*blocksize, yend, CDF2, &state, &list);
		}
		if(cx*blocksize<bw)
		{
			for(int kb=cy-1;kb>0;--kb)
			{
				bx1=cx*blocksize, by1=(kb-1)*blocksize, bx2=cx*blocksize, by2=kb*blocksize;
				t21_prepblock(buf2, CDF+((size_t)kc<<8), bw, bh, kc, bx1, by1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
				t21_encblock(buf2, bw, kc, bx2*blocksize, xend, by2*blocksize, yend, CDF2, &state, &list);
			}
			bx2=cx*blocksize, by2=0;
			t21_prepblock(buf2, CDF+((size_t)kc<<8), bw, bh, kc, -1, -1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
			t21_encblock(buf2, bw, kc, bx2*blocksize, xend, by2*blocksize, yend, CDF2, &state, &list);
		}
		unsigned short b1, b2;
		for(int kb=(int)order->count-1;kb>0;--kb)
		{
			b2=*(unsigned short*)array_at(&order, kb  );
			b1=*(unsigned short*)array_at(&order, kb-1);
			bx2=b2%cx, by2=b2/cx;
			bx1=b1%cx, by1=b1/cx;
			t21_prepblock(buf2, CDF+((size_t)kc<<8), bw, bh, kc, bx1, by1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
			csize+=t21_encblock(buf2, bw, kc, bx2*blocksize, xend, by2*blocksize, yend, CDF2, &state, &list);
			++blockcount;
			//printf("%3d %lf\n", kb, csize);//
		}
		b2=*(unsigned short*)array_at(&order, 0);
		bx2=b2%cx, by2=b2/cx;
		t21_prepblock(buf2, CDF+((size_t)kc<<8), bw, bh, kc, -1, -1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
		t21_encblock(buf2, bw, kc, bx2*blocksize, xend, by2*blocksize, yend, CDF2, &state, &list);

#if 0
		for(int by=bycount-1;by>=0;--by)
		{
			int ky=by*blocksize;
			for(int bx=bxcount-1;bx>=0;--bx)//for each block
			{
				//if(kc==0&&bx==0&&by==1)
				//	kc=0;

				int kx=bx*blocksize;
				int xend=0, yend=0;
				t16_prepblock(buf2, CDF+((size_t)kc<<8), bw, bh, kc, bx, by, alpha, blocksize, margin, CDF2, &xend, &yend);

				//if(kc==0&&bx==0&&by==0)
				//	printf("");
				//	print_CDF(CDF2, buf2, bw, bh, kc, kx, xend, ky, yend);

				//encode block
				for(int ky2=yend-1;ky2>=ky;--ky2)
				{
					for(int kx2=xend-1;kx2>=kx;--kx2)//for each pixel
					{
						unsigned char sym=buf2[(bw*ky2+kx2)<<2|kc];

						int cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;

						//if(kc==0&&ky2==0&&kx2==0)//
						//	printf("sym 0x%02X cdf 0x%04X freq 0x%04X\n", sym, cdf, freq);

						if(!freq)
							LOG_ERROR("ZPS");
						
						if(state>=(unsigned)(freq<<16))//renorm
						{
							dlist_push_back(&list, &state, 2);
							state>>=16;
						}
						debug_enc_update(state, cdf, freq, kx2, ky2, 0, kc, sym);
						state=state/freq<<16|(cdf+state%freq);//update
					}
				}
			}
		}
#endif
		dlist_push_back(&list, &state, 4);
		ansbookmarks[kc]=(int)list.nobj;

		if(loud)
			printf("C%d %lf\n", kc, csize/(blockcount*blocksize*blocksize));
	}
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, ansbookmarks, 12);
	memcpy(data[0]->data+dststart+12, &order->count, 4);
	//memcpy(data[0]->data+dststart+12, &swaps->count, 4);

	if(loud)
	{
		int overhead=16+(int)order->count*2+768*2;
		printf("alpha 0x%04X block %d\n", alpha, blocksize);
		printf("Overhead %7d\n", overhead);
		printf("Red      %7d\n", ansbookmarks[0]-overhead);
		printf("Green    %7d\n", ansbookmarks[1]-ansbookmarks[0]);
		printf("Blue     %7d\n", ansbookmarks[2]-ansbookmarks[1]);
	}

	dlist_clear(&list);
	array_free(&order);
	//array_free(&swaps);
	free(CDF2);
	free(CDF);
	free(hist);
	free(buf2);
	return 1;
}
int    test21_decode(const unsigned char *data, size_t srclen, int bw, int bh, int alpha, int blocksize, unsigned char *buf)
{
	int cdflen=768LL*sizeof(short), overhead=12LL+cdflen;
	int res=bw*bh;
	
	const unsigned char *srcptr, *srcstart, *srcend=data+srclen;
	if(data+overhead>=srcend)
	{
		LOG_ERROR("Invalid file");
		return 0;
	}

	int ansbookmarks[3], blockcount;
	memcpy(ansbookmarks, data, 12);
	memcpy(&blockcount, data+12, 4);
	if((unsigned)blockcount>(unsigned)(bw*bh))
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	overhead+=blockcount*sizeof(short);
	if((unsigned)ansbookmarks[2]<(unsigned)overhead)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	if((unsigned)ansbookmarks[2]>(unsigned)srclen)
	{
		LOG_ERROR("Incomplete file");
		return 0;
	}
	
	unsigned short *CDF=(unsigned short*)malloc(cdflen);
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	const unsigned char *orderptr=data+16;
	memcpy(CDF, data+16+blockcount*sizeof(short), cdflen);
	
	int cx=(bw+blocksize-1)/blocksize,
		cy=(bh+blocksize-1)/blocksize;
	for(int kc=0;kc<3;++kc)
	{
		unsigned state;
		srcptr=data+ansbookmarks[kc];
		srcstart=kc?data+ansbookmarks[kc-1]:data+overhead;
		srcptr-=4;
		if(srcptr<srcstart)
			LOG_ERROR("ANS buffer overflow");
		memcpy(&state, srcptr, 4);
		
		unsigned short b1, b2;
		int bx2, by2, bx1, by1;
		int xend=0, yend=0;

		memcpy(&b2, orderptr, 2);
		bx2=b2%cx, by2=b2/cx;
		t21_prepblock(buf, CDF+((size_t)kc<<8), bw, bh, kc, -1, -1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
		t21_decblock(bw, kc, bx2*blocksize, xend, by2*blocksize, yend, &state, CDF2, buf, &srcptr, srcstart);
		for(int kb=1;kb<blockcount;++kb)
		{
			memcpy(&b1, orderptr+(kb-1)*sizeof(short), 2);
			memcpy(&b2, orderptr+kb*sizeof(short), 2);
			bx1=b1%cx, by1=b1/cx;
			bx2=b2%cx, by2=b2/cx;
			t21_prepblock(buf, CDF+((size_t)kc<<8), bw, bh, kc, bx1, by1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
			t21_decblock(bw, kc, bx2*blocksize, xend, by2*blocksize, yend, &state, CDF2, buf, &srcptr, srcstart);
		}
		if(cx*blocksize<bw)
		{
			bx2=cx*blocksize, by2=0;
			t21_prepblock(buf, CDF+((size_t)kc<<8), bw, bh, kc, -1, -1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
			t21_decblock(bw, kc, bx2*blocksize, xend, by2*blocksize, yend, &state, CDF2, buf, &srcptr, srcstart);
			for(int kb=1;kb<cy;++kb)
			{
				bx1=cx*blocksize, by1=(kb-1)*blocksize, bx2=cx*blocksize, by2=kb*blocksize;
				t21_prepblock(buf, CDF+((size_t)kc<<8), bw, bh, kc, bx1, by1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
				t21_decblock(bw, kc, bx2*blocksize, xend, by2*blocksize, yend, &state, CDF2, buf, &srcptr, srcstart);
			}
		}
		if(cy*blocksize<bh)
		{
			bx2=0, by2=cy*blocksize;
			t21_prepblock(buf, CDF+((size_t)kc<<8), bw, bh, kc, -1, -1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
			t21_decblock(bw, kc, bx2*blocksize, xend, by2*blocksize, yend, &state, CDF2, buf, &srcptr, srcstart);
			for(int kb=1;kb<cx;++kb)
			{
				bx1=(kb-1)*blocksize, by1=cy*blocksize, bx2=kb*blocksize, by2=cy*blocksize;
				t21_prepblock(buf, CDF+((size_t)kc<<8), bw, bh, kc, bx1, by1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
				t21_decblock(bw, kc, bx2*blocksize, xend, by2*blocksize, yend, &state, CDF2, buf, &srcptr, srcstart);
			}
			if(cx*blocksize<bw)
			{
				bx1=(cx-1)*blocksize, by1=(cy-1)*blocksize, bx2=cx*blocksize, by2=cy*blocksize;
				t21_prepblock(buf, CDF+((size_t)kc<<8), bw, bh, kc, bx1, by1, bx2, by2, alpha, blocksize, CDF2, &xend, &yend);
				t21_decblock(bw, kc, bx2*blocksize, xend, by2*blocksize, yend, &state, CDF2, buf, &srcptr, srcstart);
			}
		}
#if 0
		for(int by=0;by<bycount;++by)
		{
			int ky=by*blocksize;
			for(int bx=0;bx<bxcount;++bx)//for each block
			{
				//if(kc==0&&bx==0&&by==1)
				//	kc=0;

				int kx=bx*blocksize;
				int xend=0, yend=0;
				t16_prepblock(buf, CDF+((size_t)kc<<8), bw, bh, kc, bx, by, alpha, blocksize, margin, CDF2, &xend, &yend);
				for(int ky2=ky;ky2<yend;++ky2)
				{
					for(int kx2=kx;kx2<xend;++kx2)//for each pixel
					{
						unsigned c=(unsigned short)state;
						int sym=0;

						int L=0, R=256, found=0;
						while(L<=R)
						{
							sym=(L+R)>>1;
							if(CDF2[sym]<c)
								L=sym+1;
							else if(CDF2[sym]>c)
								R=sym-1;
							else
							{
								found=1;
								break;
							}
						}
						if(!found)
							sym=L+(L<256&&CDF2[L]<c)-1;
						else
							for(;sym<256-1&&CDF2[sym+1]==c;++sym);

						//if(sym!=debug_ptr[(bw*ky2+kx2)<<2|kc])//
						//	LOG_ERROR("");

						buf[(bw*ky2+kx2)<<2|kc]=(unsigned char)sym;

						unsigned cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;

						debug_dec_update(state, cdf, freq, kx2, ky2, 0, kc, sym);
						state=freq*(state>>16)+c-cdf;//update
						if(state<0x10000)//renorm
						{
							state<<=16;
							if(srcptr-2>=srcstart)
							{
								srcptr-=2;
								memcpy(&state, srcptr, 2);
							}
						}
					}
				}
			}
		}
#endif
	}
//#ifndef T21_DISABLE_SHUFFLE
//	int cx=bw/blocksize;
//	for(int ks=swapcount-1;ks>=0;--ks)
//	{
//		SwapInfo *p=(SwapInfo*)swapptr+ks;
//		int b1=bxcount*p->b1.y+p->b1.x,
//			b2=bxcount*p->b2.y+p->b2.x;
//		t21_blockswap(buf, bw, b1, b2, blocksize);
//	}
//#endif
	free(CDF);
	free(CDF2);
	for(int k=0;k<res;++k)//set alpha
		buf[k<<2|3]=0xFF;
#ifndef T21_DISABLE_TRANSFORMS
	apply_transforms_inv(buf, bw, bh);
#endif
	return 1;
}



//	#define T22_VISITED
//	#define T22_BLOCKTEST

//test 22: interlaced parallelogram blocks
void t22_accumulaterow(const unsigned char *buf, int bw, int bh, int kc, int ky, int x1, int x2, unsigned *hist)
{
	if((unsigned)ky<(unsigned)bh)
	{
		if(x1<0)
			x1=0;
		if(x2>bw)
			x2=bw;
		//const unsigned char *row=buf+(bw*ky<<2);
		for(int kx=x1;kx<x2;++kx)
		{
			unsigned char sym=buf[(bw*ky+kx)<<2|kc];
			++hist[sym];
			++hist[256];
		}
	}
}
void t22_prepblock(const unsigned char *buf2, const unsigned short *CDF0, int bw, int bh, int kc, int bx, int by, int alpha, int blockw, int blockh, unsigned *CDF2, int *ybounds, int *xbounds)//block has even dimensions, ybounds[2], xbounds[blockh]
{
	int overflow=0;

	//convention: even block is above odd block
	int oddflag=bx&1;
	ybounds[0]=by*blockh+oddflag;
	ybounds[1]=ybounds[0]+blockh-1;
	xbounds[0]=bx*(blockw>>1)-(blockw>>1)-oddflag;
	xbounds[1]=xbounds[0]+blockw;
	for(int ky=2;ky<blockh;ky+=2)
	{
		xbounds[ky  ]=xbounds[ky-2]-2;
		xbounds[ky+1]=xbounds[ky  ]+blockw;
	}

	memset(CDF2, 0, 257*sizeof(int));
	int ystart=by*blockh;
	t22_accumulaterow(buf2, bw, bh, kc, ystart-2, xbounds[0], xbounds[1], CDF2);
	t22_accumulaterow(buf2, bw, bh, kc, ystart-1, xbounds[0]-1, xbounds[1]+1, CDF2);
	for(int ky=0;ky<blockh;ky+=2)
	{
		int xstart, xend1, xend2;

		if(oddflag)//odd block
			xstart=xbounds[ky]-1, xend1=xbounds[ky]+(blockw>>1)+1, xend2=xbounds[ky];
		else//even block
			xstart=xbounds[ky]-2, xend1=xbounds[ky], xend2=xbounds[ky]+(blockw>>1)-1;

		t22_accumulaterow(buf2, bw, bh, kc, ystart+ky, xstart, xend1, CDF2);
		t22_accumulaterow(buf2, bw, bh, kc, ystart+ky+1, xstart-1, xend2, CDF2);
	}

	int ylimit=bh-(oddflag==(bh&1));//if odd block && odd bh || even block && even bh: subtract 1

	//int ylimit=bh&~oddflag;//if odd block && odd bh: subtract 1
	//if(!oddflag)//if even block && even bh: subtract 1
	//	ylimit-=!(ylimit&1);
	
	if(ybounds[1]>ylimit)//clamp bounds
		ybounds[1]=ylimit;
	for(int ky=0;ky<blockh;ky+=2)
	{
		xbounds[ky  ]=CLAMP(0, xbounds[ky  ], bw);
		xbounds[ky|1]=CLAMP(0, xbounds[ky|1], bw);
	}

	int cdf1, f1, f2, freq, sum;
	if(CDF2[256])//histogram exists
	{
		sum=0;
		for(int sym=0;sym<256;++sym)
		{
			cdf1=!overflow?CDF0[sym]:0x10000;
			if(sym<255)
				overflow|=cdf1>CDF0[sym+1];
			f1=(sym<255&&!overflow?CDF0[sym+1]:0x10000)-cdf1;

			f2=(int)(((long long)CDF2[sym]<<16)/CDF2[256]);//normalize

			freq=f1+(int)(((long long)f2-f1)*alpha>>16);//blend

			freq=((unsigned)(freq*0xFF00)>>16)+1;//guard
			if(freq<0||freq>0xFF01)
				LOG_ERROR("Impossible freq 0x%04X / 0x10000", freq);
			CDF2[sym]=sum;
			sum+=freq;
			if(sum>0x10000)
				LOG_ERROR("ANS CDF sum 0x%04X, freq 0x%04X", sum, freq);
		}
		CDF2[256]=0x10000;
	}
	else
	{
		for(int sym=0;sym<256;++sym)
		{
			if(overflow)
				CDF2[sym]=0xFF00|sym;
			else
			{
				int cdf=CDF0[sym];
				CDF2[sym]=((unsigned)(cdf*0xFF00)>>16)+sym;
				if(sym<255)
					overflow|=cdf>CDF0[sym+1];
			}
		}
		CDF2[256]=0x10000;
	}
}
size_t test22_encode(const unsigned char *src, int bw, int bh, int alpha, int *blockw, int *blockh, ArrayHandle *data, int loud, int *csizes)
{
	int res=bw*bh;
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	unsigned short *CDF0=(unsigned short*)malloc(768LL*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	
#ifndef T22_BLOCKTEST
#ifdef T22_VISITED
	char *mask=(char*)malloc(res);//REMOVEME
#endif
#endif

	int hmax=blockh[0];
	if(hmax<blockh[1])
		hmax=blockh[1];
	if(hmax<blockh[2])
		hmax=blockh[2];
	int *xbounds=(int*)malloc(hmax*sizeof(int));

	if(!buf2||!CDF0||!CDF2||!xbounds)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
#ifndef T21_DISABLE_TRANSFORMS
	apply_transforms_fwd(buf2, bw, bh);
#endif
	for(int kc=0;kc<3;++kc)
	{
		memset(CDF2, 0, 256LL*sizeof(unsigned));
		for(int k=0;k<res;++k)
		{
			unsigned char sym=buf2[k<<2|kc];
			++CDF2[sym];
		}
		t21_normalize_histogram(CDF2, 256, res, CDF0+((size_t)kc<<8));
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);

	int ansbookmarks[3]={0};
	dlist_push_back(&list, 0, 12);
	dlist_push_back(&list, CDF0, 768*sizeof(short));
	
#ifdef T22_BLOCKTEST
	memset(buf2, 0, (size_t)res<<2);
#endif

	int ybounds[2];
	for(int kc=0;kc<3;++kc)
	{
#ifndef T22_BLOCKTEST
#ifdef T22_VISITED
		memset(mask, 0, res);//
#endif
#endif

		unsigned state=0x10000;
		int cx=(bw+blockh[kc]+(blockw[kc]>>1)+(blockw[kc]>>1)-1)/(blockw[kc]>>1),
			cy=(bh+blockh[kc]-1)/blockh[kc];
		for(int by=cy-1;by>=0;--by)
		{
			//if(kc==1&&by==1)
			//	printf("");
			int kblock=0;
			for(int bx=cx-1;bx>=-5;--bx)
			{
				//if(by==1&&bx==0)
				//	printf("");
				//if(kc==1&&by==1&&bx==6)
				//	printf("");
				t22_prepblock(buf2, CDF0+((size_t)kc<<8), bw, bh, kc, bx, by, alpha, blockw[kc], blockh[kc], CDF2, ybounds, xbounds);
				for(int ky=ybounds[1]-1, yidx=(ybounds[1]-ybounds[0])>>1;ky>=ybounds[0];ky-=2, --yidx)
				{
					int x1=xbounds[yidx<<1], x2=xbounds[yidx<<1|1];
					for(int kx=x2-1;kx>=x1;--kx)
					{
						//if(kc==1&&kx==767&&ky==1)
						//	printf("");
#ifdef T22_BLOCKTEST
						if(buf2[(bw*ky+kx)<<2|kc])
							LOG_ERROR("Already visited CXY %d %d %d", kc, kx, ky);
						buf2[(bw*ky+kx)<<2|kc]+=((unsigned char)kblock*64)%(255-32)+32;
#else
#ifdef T22_VISITED
						if(mask[bw*ky+kx])//
							LOG_ERROR("Pixel already encoded");
						mask[bw*ky+kx]=1+kblock%255;
#endif

						unsigned char sym=buf2[(bw*ky+kx)<<2|kc];

						int cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;

						//if(kc==0&&ky==511&&kx==512)//
						//	printf("CXY %d %d %d  sym 0x%02X cdf 0x%04X freq 0x%04X  state 0x%08X\n", kc, kx, ky, sym, cdf, freq, state);

						if(!freq)
							LOG_ERROR("ZPS");

						//double p=freq/65536.;
						//csize-=log2(p);
						//unsigned s0=state;
						
						if(state>=(unsigned)(freq<<16))//renorm
						{
							dlist_push_back(&list, &state, 2);
							state>>=16;
						}
						debug_enc_update(state, cdf, freq, kx, ky, 0, kc, sym);
						state=state/freq<<16|(cdf+state%freq);//update
						//if(!state)
						//	printf("");
#endif
					}
				}
				++kblock;
			}
		}
#ifdef T22_BLOCKTEST
		for(int k=0;k<res;++k)
		{
			if(!buf2[k<<2|kc])
				LOG_ERROR("Missed CXY %d %d %d", kc, k%bw, k/bw);
		}
#else
		dlist_push_back(&list, &state, 4);
		ansbookmarks[kc]=(int)list.nobj;
#ifdef T22_VISITED
		for(int k=0;k<res;++k)
		{
			if(!mask[k])
				LOG_ERROR("Missed %d %d", k%bw, k/bw);
		}
#endif
#endif
	}
#ifdef T22_BLOCKTEST
	for(int k=0;k<res;++k)
		buf2[k<<2|3]=0xFF;
	lodepng_encode_file("T22.PNG", buf2, bw, bh, LCT_RGBA, 8);
#else
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, ansbookmarks, 12);
#endif
	
	int overhead=12+768*2;
	int ch[]=
	{
		ansbookmarks[0]-overhead,
		ansbookmarks[1]-ansbookmarks[0],
		ansbookmarks[2]-ansbookmarks[1],
	};
	if(csizes)
	{
		csizes[0]=ch[0];
		csizes[1]=ch[1];
		csizes[2]=ch[2];
	}
	if(loud)
	{
		printf("alpha 0x%04X block %dx%d %dx%d %dx%d\n", alpha, blockw[0], blockh[0], blockw[1], blockh[1], blockw[2], blockh[2]);
		printf("Overhead %7d\n", overhead);
		printf("Red      %7d  %lf\n", ch[0], (double)res/ch[0]);
		printf("Green    %7d  %lf\n", ch[1], (double)res/ch[1]);
		printf("Blue     %7d  %lf\n", ch[2], (double)res/ch[2]);
	}

	dlist_clear(&list);
	free(buf2);
	free(CDF0);
	free(CDF2);
	free(xbounds);
#ifndef T22_BLOCKTEST
#ifdef T22_VISITED
	free(mask);//REMOVEME
#endif
#endif
	return 1;
}
int    test22_decode(const unsigned char *data, size_t srclen, int bw, int bh, int alpha, int *blockw, int *blockh, unsigned char *buf)
{
	int cdflen=768LL*sizeof(short), overhead=12LL+cdflen;
	int res=bw*bh;
	
	const unsigned char *srcptr, *srcstart, *srcend=data+srclen;
	if(data+overhead>=srcend)
	{
		LOG_ERROR("Invalid file");
		return 0;
	}

	int ansbookmarks[3];
	memcpy(ansbookmarks, data, 12);
	if((unsigned)ansbookmarks[2]<(unsigned)overhead)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	if((unsigned)ansbookmarks[2]>(unsigned)srclen)
	{
		LOG_ERROR("Incomplete file");
		return 0;
	}
	
	unsigned short *CDF0=(unsigned short*)malloc(cdflen);
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	
#ifdef T22_VISITED
	char *mask=(char*)malloc(res);//REMOVEME
#endif

	int hmax=blockh[0];
	if(hmax<blockh[1])
		hmax=blockh[1];
	if(hmax<blockh[2])
		hmax=blockh[2];
	int *xbounds=(int*)malloc(hmax*sizeof(int));

	if(!CDF0||!CDF2||!xbounds)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	const unsigned char *orderptr=data+12;
	memcpy(CDF0, data+12, cdflen);
	
	int ybounds[2];
	for(int kc=0;kc<3;++kc)
	{
		unsigned state;
		srcptr=data+ansbookmarks[kc];
		srcstart=kc?data+ansbookmarks[kc-1]:data+overhead;
		srcptr-=4;
		if(srcptr<srcstart)
			LOG_ERROR("ANS buffer overflow");
		memcpy(&state, srcptr, 4);
		
#ifdef T22_VISITED
		memset(mask, 0, res);//
#endif

		int cx=(bw+blockh[kc]+(blockw[kc]>>1)+(blockw[kc]>>1)-1)/(blockw[kc]>>1)+1,
			cy=(bh+blockh[kc]-1)/blockh[kc];
		for(int by=0;by<cy;++by)
		{
			int kblock=0;
			for(int bx=-5;bx<cx;++bx)
			{
				//if(kc==1&&by==0&&bx==6)//
				//	printf("");

				t22_prepblock(buf, CDF0+((size_t)kc<<8), bw, bh, kc, bx, by, alpha, blockw[kc], blockh[kc], CDF2, ybounds, xbounds);
				for(int ky=ybounds[0], yidx=0;ky<ybounds[1];ky+=2, ++yidx)
				{
					int x1=xbounds[yidx<<1], x2=xbounds[yidx<<1|1];
					for(int kx=x1;kx<x2;++kx)
					{
#ifdef T22_VISITED
						if(mask[bw*ky+kx])//
							LOG_ERROR("Pixel already encoded CXY %d %d %d", kc, kx, ky);
						mask[bw*ky+kx]=1+kblock%255;
#endif

						unsigned c=(unsigned short)state;
						int sym=0;

						//if(kc==1&&ky==1&&kx==767)//
						//	printf("");

						int L=0, R=256, found=0;
						while(L<=R)
						{
							sym=(L+R)>>1;
							if(CDF2[sym]<c)
								L=sym+1;
							else if(CDF2[sym]>c)
								R=sym-1;
							else
							{
								found=1;
								break;
							}
						}
						if(!found)
							sym=L+(L<256&&CDF2[L]<c)-1;
						else
							for(;sym<256-1&&CDF2[sym+1]==c;++sym);

						//if(sym!=debug_ptr[(bw*ky2+kx2)<<2|kc])//
						//	LOG_ERROR("");

						buf[(bw*ky+kx)<<2|kc]=(unsigned char)sym;

						unsigned cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;

						debug_dec_update(state, cdf, freq, kx, ky, 0, kc, sym);
						state=freq*(state>>16)+c-cdf;//update
						if(state<0x10000)//renorm
						{
							state<<=16;
							if(srcptr-2>=srcstart)
							{
								srcptr-=2;
								memcpy(&state, srcptr, 2);
							}
						}
					}
				}
				++kblock;
			}
		}
#ifdef T22_VISITED
		//lodepng_encode_file("T22mask.PNG", mask, bw, bh, LCT_GREY, 8);
		for(int k=0;k<res;++k)
		{
			if(!mask[k])
				LOG_ERROR("Missed CXY %d %d %d", kc, k%bw, k/bw);
		}
#endif
	}
//#ifndef T21_DISABLE_SHUFFLE
//	int cx=bw/blocksize;
//	for(int ks=swapcount-1;ks>=0;--ks)
//	{
//		SwapInfo *p=(SwapInfo*)swapptr+ks;
//		int b1=bxcount*p->b1.y+p->b1.x,
//			b2=bxcount*p->b2.y+p->b2.x;
//		t21_blockswap(buf, bw, b1, b2, blocksize);
//	}
//#endif
	free(CDF0);
	free(CDF2);
#ifdef T22_VISITED
	free(mask);//REMOVEME
#endif
	for(int k=0;k<res;++k)//set alpha
		buf[k<<2|3]=0xFF;
#ifndef T21_DISABLE_TRANSFORMS
	apply_transforms_inv(buf, bw, bh);
#endif
	return 1;
}


//test 23: test 16 optimizer
#define T23_WMIN 8
#define T23_WMAX 768
void t23_prepblock(const unsigned char *b2, const unsigned short *CDF, int bw, int bh, int kc, int kx, int ky, int alpha, int blockw, int blockh, int margin, unsigned *CDF2, int *xend, int *yend)
{
	*yend=ky+blockh<=bh?ky+blockh:bh;
	*xend=kx+blockw<=bw?kx+blockw:bw;
	
	memset(CDF2, 0, 256*sizeof(unsigned));

	int count2=0;

	int left=kx-margin;
	if(left<0)
		left=0;
	int right=kx+blockw+margin;
	if(right>bw)
		right=bw;
	int top=ky-margin;
	if(top<0)
		top=0;

	if(left<kx)//if left block is available
	{
		for(int ky2=ky;ky2<*yend;++ky2)
		{
			for(int kx2=left;kx2<kx;++kx2)//for each pixel
			{
				int sym=b2[(bw*ky2+kx2)<<2|kc];
				int dist=kx-kx2;
				if(dist<0||dist>margin)
					LOG_ERROR("Wrong distance");
					
				int inc=(margin<<1|1)-dist;

				if(!inc)
					LOG_ERROR("Zero inc");

				CDF2[sym]+=inc;
				count2+=inc;
			}
		}
	}
	if(top<ky)//if top block is available
	{
		for(int ky2=top;ky2<ky;++ky2)
		{
			for(int kx2=kx;kx2<*xend;++kx2)//for each pixel
			{
				unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];
				int dist=ky-ky2;
				if(dist<0||dist>margin)
					LOG_ERROR("Wrong distance");
					
				int inc=(margin<<1|1)-dist;
					
				if(!inc)
					LOG_ERROR("Zero inc");

				CDF2[sym]+=inc;
				count2+=inc;
					
				if(count2<inc)
					LOG_ERROR("OVERFLOW");
			}
		}
	}
	if(left<kx&&top<ky)//if topleft block is available
	{
		for(int ky2=top;ky2<ky;++ky2)
		{
			for(int kx2=left;kx2<kx;++kx2)//for each pixel
			{
				unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];
				int dist=kx-kx2+ky-ky2;
					
				if(dist<0||dist>(margin<<1))
					LOG_ERROR("Wrong distance");
					
				int inc=(margin<<1|1)-dist;
					
				if(!inc)
					LOG_ERROR("Zero inc");

				CDF2[sym]+=inc;
				count2+=inc;
			}
		}
	}
	if(right>kx+blockw&&top<ky)//if topright block is available
	{
		for(int ky2=top;ky2<ky;++ky2)
		{
			for(int kx2=kx+blockw;kx2<right;++kx2)//for each pixel
			{
				unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];
				int dist=kx2-(kx+blockw)+ky-ky2;
				//int dist=MAXVAR(kx2-(kx+blockw), ky-ky2);

				if(dist<0||dist>(margin<<1))
					LOG_ERROR("Wrong distance");
					
				int inc=(margin<<1|1)-dist;
					
				if(!inc)
					LOG_ERROR("Zero inc");

				CDF2[sym]+=inc;
				count2+=inc;
			}
		}
	}
		
	int overflow=0;//CDF overflow can happen only once
	if(count2)
	{
		int sum=0;
		for(int sym=0;sym<256;++sym)
		{
			int cdf1=!overflow?CDF[sym]:0x10000;
			if(sym<255)
				overflow|=cdf1>CDF[sym+1];
			int f1=(sym<255&&!overflow?CDF[sym+1]:0x10000)-cdf1;

			int f2=(int)(((long long)CDF2[sym]<<16)/count2);//normalize

			int freq=f1+(int)(((long long)f2-f1)*alpha>>16);//blend

			freq=((unsigned)(freq*0xFF00)>>16)+1;//guard
			if(freq<0||freq>0xFF01)
				LOG_ERROR("Impossible freq 0x%04X / 0x10000", freq);
			CDF2[sym]=sum;
			sum+=freq;
			if(sum>0x10000)
				LOG_ERROR("ANS CDF sum 0x%04X, freq 0x%04X", sum, freq);
		}
		CDF2[256]=0x10000;
	}
	else
	{
		for(int sym=0;sym<256;++sym)
		{
			if(overflow)
				CDF2[sym]=0xFF00|sym;
			else
			{
				int cdf=CDF[sym];
				CDF2[sym]=((unsigned)(cdf*0xFF00)>>16)+sym;
				if(sym<255)
					overflow|=cdf>CDF[sym+1];
			}
		}
		CDF2[256]=0x10000;
	}
}
double t23_calcloss(const unsigned char *buf, unsigned short *CDF0, int iw, int ih, int kc, int kx, int ky, int alpha, int width, int margin, int wmax, unsigned *CDF2)
{
	double csize=0;
	int xend=0, yend=0;
	int xstart=MAXVAR(kx-wmax, 0), w2;
	for(int kx2=MAXVAR(kx-width, 0);;kx2-=width)
	{
		w2=width;
		if(kx2<xstart)
			w2=width-(xstart-kx2), kx2=xstart;
		if(w2>kx-kx2)
			w2=kx-kx2;
		if(!w2)
			break;
		t23_prepblock(buf, CDF0, iw, ih, kc, kx2, ky, alpha, w2, 1, margin, CDF2, &xend, &yend);
		for(int kx3=kx2;kx3<xend;++kx3)
		{
			unsigned char sym=buf[(iw*ky+kx3)<<2|kc];
			int freq=CDF2[sym+1]-CDF2[sym];
			double p=freq/65536.;
			csize-=log2(p);
		}
		if(kx2<=xstart)
			break;
	}
	return csize/8;
}
int t23_findbestwidth(const unsigned char *buf2, int iw, int ih, int kc, int kx, int ky, int alpha, int margin, int width, unsigned short *CDF0, unsigned *CDF2)
{
	double csize0, csize;
	csize=t23_calcloss(buf2, CDF0, iw, ih, kc, kx, ky, alpha, width, margin, T23_WMAX, CDF2);
	for(;width>T23_WMIN;)
	{
		csize0=csize;
		--width;
		csize=t23_calcloss(buf2, CDF0, iw, ih, kc, kx, ky, alpha, width, margin, T23_WMAX, CDF2);
		if(width<=T23_WMIN||csize>=csize0)//cancel last change and break
		{
			if(width<T23_WMIN||csize>=csize0)
			{
				++width;
				csize=csize0;
			}
			break;
		}
	}
	for(;width<T23_WMAX;)
	{
		csize0=csize;
		++width;
		csize=t23_calcloss(buf2, CDF0, iw, ih, kc, kx, ky, alpha, width, margin, T23_WMAX, CDF2);
		if(width>=T23_WMAX||csize>=csize0)
		{
			if(width>T23_WMAX||csize>=csize0)
			{
				--width;
				csize=csize0;
			}
			break;
		}
	}
	return width;
}
int whist[3][T23_WMAX-T23_WMIN+1]={0};
size_t test23_encode(const unsigned char *src, int iw, int ih, ArrayHandle *data, int loud, int *csizes)
{
	int res=iw*ih;
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	unsigned short *CDF0=(unsigned short*)malloc(768LL*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	int *gbounds=(int*)malloc(((size_t)iw+1)*sizeof(int));
	if(!buf2||!CDF0||!CDF2||!gbounds)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);

	addbuf(buf2, iw, ih, 3, 4, 128);
	colortransform_ycocb_fwd((char*)buf2, iw, ih);
	pred_jxl_opt_v2((char*)buf2, iw, ih, jxlparams_i16, loud);
	//pred_jxl_optimizeall(buf2, iw, ih, loud);
#if 1
	if(loud)
	{
		for(int kc=0;kc<3;++kc)//
		{
			for(int kp=0;kp<11;++kp)
			{
				short val=jxlparams_i16[11*kc+kp];
				printf(" %c0x%04X,", val<0?'-':' ', abs(val));
			}
			printf("\n");
		}
	}
#endif
	pred_jxl_apply((char*)buf2, iw, ih, jxlparams_i16, 1);
	addbuf(buf2, iw, ih, 3, 4, 128);
	
	for(int kc=0;kc<3;++kc)
	{
		memset(CDF2, 0, 256LL*sizeof(unsigned));
		for(int k=0;k<res;++k)
		{
			unsigned char sym=buf2[k<<2|kc];
			++CDF2[sym];
		}
		t21_normalize_histogram(CDF2, 256, res, CDF0+((size_t)kc<<8));
	}
	
	DList list;
	dlist_init(&list, 1, 1024, 0);

	int ansbookmarks[3]={0};
	dlist_push_back(&list, 0, 12);
	dlist_push_back(&list, jxlparams_i16, 33*sizeof(short));
	dlist_push_back(&list, CDF0, 768*sizeof(short));
	
	int alpha=0xD3E7, margin=37, width=23;
	for(int kc=0;kc<3;++kc)
	{
		switch(kc)
		{
		case 0:margin=26;break;
		case 1:margin=37;break;
		case 2:margin=26;break;
		}
		unsigned state=0x10000;
		for(int ky=ih-1;ky>=0;--ky)
		{
			if(loud)
				printf("CY %d %5d  %5.2lf%%\r", kc, ih-ky, 100.*(ih*kc+ih-1-ky)/(ih*3));
			switch(kc)
			{
			case 0:width= 8;break;
			case 1:width=23;break;
			case 2:width= 8;break;
			}
			int ngroups=0;
			gbounds[0]=0;
			for(int kx=0;kx<iw;++ngroups)
			{
				//if(kc==0&&ky==2&&kx==512)//
				//	printf("");

				width=t23_findbestwidth(buf2, iw, ih, kc, kx, ky, alpha, margin, width, CDF0, CDF2);

				++whist[kc][width-T23_WMIN];//

				gbounds[ngroups+1]=gbounds[ngroups]+width;
				if(gbounds[ngroups+1]>iw)
					gbounds[ngroups+1]=iw;
				kx=gbounds[ngroups+1];
			}
			for(int kg=ngroups-1;kg>=0;--kg)
			{
				int xend=0, yend=0;
				t23_prepblock(buf2, CDF0, iw, ih, kc, gbounds[kg], ky, alpha, gbounds[kg+1]-gbounds[kg], 1, margin, CDF2, &xend, &yend);
				for(int kx=gbounds[kg+1]-1;kx>=gbounds[kg];--kx)
				{
					//if(kc==0&&ky==2&&kx==512)//
					//	printf("");

					unsigned char sym=buf2[(iw*ky+kx)<<2|kc];

					int cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;

					if(!freq)
						LOG_ERROR("ZPS");

					if(state>=(unsigned)(freq<<16))//renorm
					{
						dlist_push_back(&list, &state, 2);
						state>>=16;
					}
					debug_enc_update(state, cdf, freq, kx, ky, 0, kc, sym);
					state=state/freq<<16|(cdf+state%freq);//update
				}
			}
		}
		dlist_push_back(&list, &state, 4);
		ansbookmarks[kc]=(int)list.nobj;
	}
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, ansbookmarks, 12);
	
	int overhead=12+33*2+768*2;
	int ch[]=
	{
		ansbookmarks[0]-overhead,
		ansbookmarks[1]-ansbookmarks[0],
		ansbookmarks[2]-ansbookmarks[1],
	};
	if(csizes)
	{
		csizes[0]=ch[0];
		csizes[1]=ch[1];
		csizes[2]=ch[2];
	}
	if(loud)
	{
		printf("\n");
		printf("alpha 0x%04X  margin %d\n", alpha, margin);
		printf("Total    %7d  %lf\n", ansbookmarks[2], 3.*res/ansbookmarks[2]);
		printf("Overhead %7d\n", overhead);
		printf("Red      %7d  %lf\n", ch[0], (double)res/ch[0]);
		printf("Green    %7d  %lf\n", ch[1], (double)res/ch[1]);
		printf("Blue     %7d  %lf\n", ch[2], (double)res/ch[2]);

#if 1
		printf("\tC0\tC1\tC2\n");
		for(int k=T23_WMIN;k<=T23_WMAX;++k)
		{
			int f1=whist[0][k-T23_WMIN], f2=whist[1][k-T23_WMIN], f3=whist[2][k-T23_WMIN];
			if(f1||f2||f3)
				printf("W %3d  %5d  %5d  %5d\n", k, f1, f2, f3);
		}
#endif
	}

	dlist_clear(&list);
	free(gbounds);
	free(buf2);
	free(CDF0);
	free(CDF2);
	return 1;
}
int    test23_decode(const unsigned char *data, size_t srclen, int iw, int ih, unsigned char *buf, int loud)
{
	const int cdflen=768LL*sizeof(short), overhead=12LL+33*sizeof(short)+cdflen;
	int res=iw*ih;
	
	const unsigned char *srcptr, *srcstart, *srcend=data+srclen;
	if(data+overhead>=srcend)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}

	unsigned ansbookmarks[3];
	memcpy(ansbookmarks, data, 12);
	if(ansbookmarks[2]<(unsigned)overhead||ansbookmarks[2]>srclen)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	
	unsigned short *CDF0=(unsigned short*)malloc(cdflen);
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	short jxlparams[33];
	memcpy(jxlparams, data+12, 33*sizeof(short));
	memcpy(CDF0, data+12+33*sizeof(short), cdflen);
	
	int alpha=0xD3E7, margin=37, width=23;
	for(int kc=0;kc<3;++kc)
	{
		switch(kc)
		{
		case 0:margin=26;break;
		case 1:margin=37;break;
		case 2:margin=26;break;
		}

		unsigned state;
		srcptr=data+ansbookmarks[kc];
		srcstart=kc?data+ansbookmarks[kc-1]:data+overhead;
		srcptr-=4;
		if(srcptr<srcstart)
			LOG_ERROR("ANS buffer overflow");
		memcpy(&state, srcptr, 4);
		
		for(int ky=0;ky<ih;++ky)
		{
			if(loud)
				printf("CY %d %5d  %5.2lf%%\r", kc, ky, 100.*(ih*kc+ky)/(ih*3));
			switch(kc)
			{
			case 0:width= 8;break;
			case 1:width=23;break;
			case 2:width= 8;break;
			}
			for(int kx=0;kx<iw;)//for each block
			{
				//if(kc==0&&kx>=512&&ky==2)
				//	printf("");

				width=t23_findbestwidth(buf, iw, ih, kc, kx, ky, alpha, margin, width, CDF0, CDF2);

				int xend=0, yend=0;
				t23_prepblock(buf, CDF0, iw, ih, kc, kx, ky, alpha, width, 1, margin, CDF2, &xend, &yend);
				for(;kx<xend;++kx)//for each pixel
				{
					unsigned c=(unsigned short)state;
					int sym=0;
					
					//if(kc==0&&ky==2&&kx==512)//
					//	printf("");

					int L=0, R=256, found=0;
					while(L<=R)
					{
						sym=(L+R)>>1;
						if(CDF2[sym]<c)
							L=sym+1;
						else if(CDF2[sym]>c)
							R=sym-1;
						else
						{
							found=1;
							break;
						}
					}
					if(!found)
						sym=L+(L<256&&CDF2[L]<c)-1;
					else
						for(;sym<256-1&&CDF2[sym+1]==c;++sym);

					buf[(iw*ky+kx)<<2|kc]=(unsigned char)sym;

					unsigned cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;
						
					debug_dec_update(state, cdf, freq, kx, ky, 0, kc, sym);
					state=freq*(state>>16)+c-cdf;//update
					if(state<0x10000)//renorm
					{
						state<<=16;
						if(srcptr-2>=srcstart)
						{
							srcptr-=2;
							memcpy(&state, srcptr, 2);
						}
					}
				}
			}
		}
	}
	free(CDF0);
	free(CDF2);

	addbuf(buf, iw, ih, 3, 4, 128);
	pred_jxl_apply((char*)buf, iw, ih, jxlparams, 0);
	colortransform_ycocb_inv((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);

	for(int k=0;k<res;++k)//set alpha
		buf[k<<2|3]=0xFF;

	return 1;
}


//experiment 24: test16 with a mix of adaptive and pre-calculated: alpha, group width and margin		group height is fixed at 1
void e24_addhist_unchecked(const unsigned char *buf2, int iw, int ih, int kc, int x1, int x2, int y1, int y2, int x0, int y0, int maxinc, unsigned *CDF2)
{
	if(x1<0)
		x1=0;
	if(x2>iw)
		x2=iw;
	if(y1<0)
		y1=0;
	if(y2>ih)
		y2=ih;
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			unsigned char sym=buf2[(iw*ky+kx)<<2|kc];
			int inc=maxinc-(abs(ky-y0)+abs(kx-x0));//X
			if(inc>0)
			{
				CDF2[sym]+=inc;
				CDF2[256]+=inc;
			}
		}
	}
}
void e24_prepblock(const unsigned char *buf2, const unsigned short *CDF0, int iw, int ih, int kc, int x, int y, int maxinc, T24Params const *p, unsigned *CDF2)
{
	int overflow=0;
	int sum, cdf1, f1, f2, freq;
	if(p->alpha&&p->mleft&&p->mtop)
	{
		memset(CDF2, 0, 257*sizeof(unsigned));
		int x0=x+(p->gwidth>>1);
		if(p->mtop)
			e24_addhist_unchecked(buf2, iw, ih, kc, x-p->mleft, x+p->gwidth+p->mright, y-p->mtop, y, x0, y, maxinc, CDF2);
		if(p->mleft)
			e24_addhist_unchecked(buf2, iw, ih, kc, x-p->mleft, x, y, y+1, x0, y, maxinc, CDF2);
		if(!CDF2[256])
			goto just_static;

		sum=0;
		for(int sym=0;sym<256;++sym)
		{
			cdf1=!overflow?CDF0[sym]:0x10000;
			if(sym<255)
				overflow|=cdf1>CDF0[sym+1];
			f1=(sym<255&&!overflow?CDF0[sym+1]:0x10000)-cdf1;

			f2=(int)(((long long)CDF2[sym]<<16)/CDF2[256]);//normalize

			freq=f1+(int)(((long long)f2-f1)*p->alpha/0xFF);//blend

			freq=((unsigned)(freq*0xFF00)>>16)+1;//guard
			if(freq<0||freq>0xFF01)
				LOG_ERROR("Impossible freq 0x%04X / 0x10000", freq);
			CDF2[sym]=sum;
			sum+=freq;
			if(sum>0x10000)
				LOG_ERROR("ANS CDF sum 0x%04X, freq 0x%04X", sum, freq);
		}
	}
	else
	{
	just_static:
		for(int sym=0;sym<256;++sym)
		{
			if(overflow)
				CDF2[sym]=0xFF00|sym;
			else
			{
				int cdf=CDF0[sym];
				CDF2[sym]=((unsigned)(cdf*0xFF00)>>16)+sym;
				if(sym<255)
					overflow|=cdf>CDF0[sym+1];
			}
		}
	}
	CDF2[256]=0x10000;
}
double e24_estimate(const unsigned char *src, int iw, int ih, int cstart, int cend, unsigned char *gw0, unsigned char *maxinc, unsigned char *encounter_threshold, T24Params const *params, double *ret_csizes, int loud)//{16, 64, 0xBF}
{
	int res=iw*ih;
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	unsigned short *CDF0=(unsigned short*)malloc(256LL*6*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(512LL*sizeof(unsigned));
	T24Params *gwidths=(T24Params*)malloc(3LL*iw*ih*sizeof(T24Params));
	int *ngroups=(int*)malloc(3LL*ih*sizeof(int));
	if(!buf2||!CDF0||!CDF2||!gwidths||!ngroups)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
	apply_transforms_fwd(buf2, iw, ih);
	
	for(int kc=cstart;kc<cend;++kc)
	{
		memset(CDF2, 0, 256LL*sizeof(unsigned));
		for(int k=0;k<res;++k)
		{
			unsigned char sym=buf2[k<<2|kc];
			++CDF2[sym];
		}
		t21_normalize_histogram(CDF2, 256, res, CDF0+((size_t)kc<<8));
	}

	unsigned *CDF3=CDF2+256;
	if(params)
		goto encode;
	for(int kc=cstart;kc<cend;++kc)
	{
		for(int ky=0;ky<ih;++ky)
		{
			//if(ky==(ih>>1))//
			//	ky=ih>>1;

			int ng=0;
			for(int kx=0;kx<iw;++ng)
			{
				//if(ky==(ih>>1)&&kx>(iw>>1))//
				//	ky=ih>>1;

				memset(CDF2, 0, 256*sizeof(unsigned));
				int gend=kx+gw0[kc],//16
					unique=1;
				if(gend>iw)
					gend=iw;
				for(int kx2=kx;kx2<gend;++kx2)
				{
					unsigned char sym=buf2[(iw*ky+kx2)<<2|kc];
					++CDF2[sym];
					unique&=CDF2[sym]==1;
				}
				T24Params *p=gwidths+res*kc+iw*ky+ng;
				//if(unique)
				{
				//	for(;gend<iw;++gend)//expand group while next sym is unique
				//	{
				//		unsigned char sym=buf2[(iw*ky+kx+gend)<<2|kc];
				//		if(CDF2[sym])
				//			break;
				//		++CDF2[sym];
				//	}
				//	p->gwidth=gend-kx;
				//
				//	p->mleft=0;
				//	p->mtop=0;
				//	p->mright=0;
				//
				//	p->alpha=0;
				//}
				//else
				//{
					for(;gend<iw;++gend)//expand group while next sym is encountered before
					{
						unsigned char sym=buf2[(iw*ky+kx+gend)<<2|kc];
						if(!CDF2[sym])
							break;
						++CDF2[sym];
					}
					int gwidth=gend-kx;
					int mleft=0, mtop=0, mright=0, mcount=0, alpha=0;
					char exleft=1, extop=1, exright=1;
					int inc=maxinc[kc];//64
					int nencountered, ntotal;
					memset(CDF3, 0, 256*sizeof(unsigned));
					for(;inc>0&&(exleft||extop||exright);--inc)//expand margin left/top/right-ward as long as the percentage of values encountered in current block is above 75%
					{
						//if(exleft)
						//{
							exleft=0;
							if(kx>mleft)
							{
								nencountered=0, ntotal=0;
								for(int k=MAXVAR(ky-mtop, 0);k<=ky;++k)
								{
									unsigned char sym=buf2[(iw*k+kx-1)<<2|kc];
									int freq=CDF2[sym];
									nencountered+=freq!=0;
									++ntotal;
								}
								if(nencountered*0xFF>=encounter_threshold[kc]*ntotal)
								{
									for(int k=MAXVAR(ky-mtop, 0);k<=ky;++k)
									{
										unsigned char sym=buf2[(iw*k+kx-1)<<2|kc];
										CDF3[sym]+=inc;
										mcount+=inc;
									}
									++mleft;
									exleft=1;
								}
							}
						//}
						//if(extop)
						//{
							extop=0;
							if(ky>mtop)
							{
								nencountered=0, ntotal=0;
								for(int k=MAXVAR(kx-mleft, 0), kend=MINVAR(kx+gwidth+mright, iw);k<kend;++k)
								{
									unsigned char sym=buf2[(iw*(ky-mtop-1)+k)<<2|kc];
									int freq=CDF2[sym];
									nencountered+=freq!=0;
									++ntotal;
								}
								if(nencountered*0xFF>=encounter_threshold[kc]*ntotal)
								{
									for(int k=MAXVAR(kx-mleft, 0), kend=MINVAR(kx+gwidth+mright, iw);k<kend;++k)
									{
										unsigned char sym=buf2[(iw*(ky-mtop-1)+k)<<2|kc];
										CDF3[sym]+=inc;
										mcount+=inc;
									}
									++mtop;
									extop=1;
								}
							}
						//}
						//if(exright)
						//{
							exright=0;
							if(mtop&&kx+gwidth<iw-mright)
							{
								nencountered=0, ntotal=0;
								for(int k=MAXVAR(ky-mtop, 0);k<ky;++k)
								{
									unsigned char sym=buf2[(iw*k+kx+gwidth+mright-1)<<2|kc];
									int freq=CDF2[sym];
									nencountered+=freq!=0;
									++ntotal;
								}
								if(nencountered*0xFF>=encounter_threshold[kc]*ntotal)
								{
									for(int k=MAXVAR(ky-mtop, 0);k<ky;++k)
									{
										unsigned char sym=buf2[(iw*k+kx+gwidth+mright-1)<<2|kc];
										CDF3[sym]+=inc;
										mcount+=inc;
									}
									++mright;
									exright=1;
								}
							}
						//}
					}
					
					//determine best alpha
					if(mcount)
					{
						double asum=0;
						int acount=0;
						for(int kx2=kx;kx2<kx+gwidth;++kx2)
						{
							unsigned char sym=buf2[(iw*ky+kx2)<<2|kc];
							double
								f_static=(double)(CDF0[kc<<8|sym]-(sym<255&&CDF0[kc<<8|sym]<CDF0[kc<<8|(sym+1)]?CDF0[kc<<8|(sym+1)]:0x10000))/65536.,
								f_causal=(double)CDF3[sym]/mcount,
								f_group=(double)CDF2[sym]/gwidth;
							if(f_causal!=f_static)
							{
								double ak=(f_group-f_static)/(f_causal-f_static);
								ak=CLAMP(0, ak, 1);
								asum+=ak;
								++acount;
							}
						}
						if(acount)
							alpha=(unsigned char)(0xFF*asum/acount);
						else
							alpha=0xFF;
					}
					else
						alpha=0;

					p->gwidth=gend-kx;

					p->mleft=mleft;
					p->mtop=mtop;
					p->mright=mright;

					p->alpha=mleft||mtop?alpha:0;//alpha is useless without the margin
				}
				kx=gend;
			}
			ngroups[ih*kc+ky]=ng;
		}

		//calculate 2nd pass CDF
		memset(CDF2, 0, 256*sizeof(unsigned));
		int bcount=0;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kg=0, ng=ngroups[ih*kc+ky], kx=0;kg<ng;++kg)
			{
				T24Params *p=gwidths+res*kc+iw*ky+kg;
				if(!p->alpha)
				{
					for(int kx2=kx;kx2<kx+p->gwidth;++kx2)
					{
						unsigned char sym=buf2[(iw*ky+kx2)<<2|kc];
						++CDF2[sym];
						++bcount;
					}
				}
				kx+=p->gwidth;
			}
		}
		t21_normalize_histogram(CDF2, 256, bcount, CDF0+((3LL+kc)<<8));
	}
	
encode:
	{
		double csizes[3]={0};
		int gcount[3]={0}, zero_alpha[3]={0};
		double av_gw[3]={0}, av_margin[9]={0}, av_alpha[3]={0};
		for(int kc=cstart;kc<cend;++kc)
		{
			double chsize=0;
			for(int ky=0;ky<ih;++ky)
			{
				for(int kg=0, ng=params?(iw+params[kc].gwidth-1)/params[kc].gwidth:ngroups[ih*kc+ky], kx=0;kg<ng;++kg)
				{
					T24Params const *p=params?params+kc:gwidths+res*kc+iw*ky+kg;
					
					int gwidth=p->gwidth;
					if(gwidth>iw-kx)
						gwidth=iw-kx;

					av_gw[kc]+=gwidth;
					av_margin[3*kc  ]+=p->mleft;
					av_margin[3*kc+1]+=p->mtop;
					av_margin[3*kc+2]+=p->mright;
					av_alpha[kc]+=p->alpha;
					++gcount[kc];
					zero_alpha[kc]+=!p->alpha;

					e24_prepblock(buf2, CDF0+(((size_t)kc+(!params&&p->alpha?3:0))<<8), iw, ih, kc, kx, ky, maxinc[kc], p, CDF2);
					int kx2=kx;
					for(;kx2<kx+gwidth;++kx2)
					{
						unsigned char sym=buf2[(iw*ky+kx2)<<2|kc];
						int freq=CDF2[sym+1]-CDF2[sym];
						double p=(double)freq/0x10000, bitsize=-log2(p);//Zipf's law
						chsize+=bitsize;
					}
					kx+=gwidth;
				}
			}
			csizes[kc]=chsize/8;
			if(ret_csizes)
				ret_csizes[kc]=csizes[kc];
		}
		double csize=csizes[0]+csizes[1]+csizes[2];

		int usize=res*3;
		if(loud)
		{
			printf("T %14lf  CR %lf\n", csize, usize/csize);
			printf("R %14lf  CR %lf  W %lf  M %10lf %10lf %10lf  A %lf%% Z %lf%%\n", csizes[0], res/csizes[0], av_gw[0]/gcount[0], av_margin[0]/gcount[0], av_margin[1]/gcount[0], av_margin[2]/gcount[0], 100.*av_alpha[0]/(0xFF*gcount[0]), 100.*zero_alpha[0]/gcount[0]);
			printf("G %14lf  CR %lf  W %lf  M %10lf %10lf %10lf  A %lf%% Z %lf%%\n", csizes[1], res/csizes[1], av_gw[1]/gcount[1], av_margin[3]/gcount[1], av_margin[4]/gcount[1], av_margin[5]/gcount[1], 100.*av_alpha[1]/(0xFF*gcount[1]), 100.*zero_alpha[1]/gcount[1]);
			printf("B %14lf  CR %lf  W %lf  M %10lf %10lf %10lf  A %lf%% Z %lf%%\n", csizes[2], res/csizes[2], av_gw[2]/gcount[2], av_margin[6]/gcount[2], av_margin[7]/gcount[2], av_margin[8]/gcount[2], 100.*av_alpha[2]/(0xFF*gcount[2]), 100.*zero_alpha[2]/gcount[2]);
		}

		free(buf2);
		free(CDF0);
		free(CDF2);
		free(gwidths);
		free(ngroups);
		return csize;
	}
}


//test 25 (the 4th optimizer for test16): image is divided into large blocks, each lblock has an sblock selected (least compressible part) on which the T16 params are optimized
typedef struct T25ParamsStruct
{
	short
		gwidth,//>=1
		mleft,
		mtop,
		mright,
		alpha,//0~0xFF
		maxinc;//>=1;
} T25Params;
typedef struct T25ParamsPackedStruct
{
	unsigned char
		gwidth,//>=1
		mleft,
		mtop,
		mright,
		alpha,//0~0xFF
		maxinc;//>=1;
} T25ParamsPacked;
typedef struct RectStruct
{
	int x1, x2, y1, y2;
} Rect;
#define T25_PARAM(P, IDX) ((short*)(P))[IDX]
static T25Params t25_limits={32, 40, 40, 40, 255, 96};
static int t25_ctr=0;
int t25_incparam(T25Params *param, int pidx, int step)
{
	int prevval=0;
	switch(pidx)
	{
	case 0:prevval=param->gwidth, param->gwidth+=step; if(param->gwidth<8)param->gwidth=8; break;
	case 1:prevval=param->mleft , param->mleft +=step; if(param->mleft <0)param->mleft =0; break;
	case 2:prevval=param->mtop  , param->mtop  +=step; if(param->mtop  <0)param->mtop  =0; break;
	case 3:prevval=param->mright, param->mright+=step; if(param->mright<0)param->mright=0; break;
	case 4:prevval=param->alpha , param->alpha +=step; if(param->alpha <0)param->alpha =0; else if(param->alpha>0xFF)param->alpha=0xFF; break;
	case 5:prevval=param->maxinc, param->maxinc+=step; if(param->maxinc<1)param->maxinc=1; break;
	}
	if(T25_PARAM(param, pidx)>T25_PARAM(&t25_limits, pidx))
		T25_PARAM(param, pidx)=T25_PARAM(&t25_limits, pidx);

	if(!BETWEEN(1, param->gwidth, t25_limits.gwidth)||!BETWEEN(0, param->mleft, t25_limits.mleft)||!BETWEEN(0, param->mtop, t25_limits.mtop)||!BETWEEN(0, param->mright, t25_limits.mright)||!BETWEEN(0, param->alpha, t25_limits.alpha)||!BETWEEN(0, param->maxinc, t25_limits.maxinc))
		LOG_ERROR("Invalid params INC  W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d", param->gwidth, param->mleft, param->mtop, param->mright, param->alpha, param->maxinc);
	
	return prevval;
}
void t25_normalize_histogram(unsigned *srchist, int nlevels, int nsymbols, unsigned short *CDF)//hist is unsigned char due to alignment issues, but it's 16bit
{
	if(!nsymbols)//bypass
	{
		for(int k=0;k<nlevels;++k)
			CDF[k]=(unsigned short)(k<<8);
		return;
	}
	unsigned sum=0, qfreq;
	for(int sym=0;sym<nlevels;++sym)
	{
		qfreq=((long long)srchist[sym]<<16)/nsymbols;
		CDF[sym]=sum;
		sum+=qfreq;
	}
}
void t25_addhist(const unsigned char *buf2, int iw, int ih, int kc, int x1, int x2, int y1, int y2, int x0a, int x0b, int y0, int maxinc, unsigned *CDF2)
{
	if(x1<0)
		x1=0;
	if(x2>iw)
		x2=iw;
	if(y1<0)
		y1=0;
	if(y2>ih)
		y2=ih;
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			unsigned char sym=buf2[(iw*ky+kx)<<2|kc];
			int dist=abs(ky-y0);
			if(kx<x0a)
				dist+=abs(kx-x0a);
			else if(kx>x0b)
				dist+=abs(kx-x0b);
			int inc=maxinc-dist;
			if(inc>0)
			{
				CDF2[sym]+=inc;
				CDF2[256]+=inc;
			}
		}
	}
}
int t25_prepblock(const unsigned char *buf2, const unsigned short *CDF0, int iw, int ih, int kc, int x1, int x2, int y, T25Params const *p, unsigned *CDF2)
{
	int overflow=0;
	int sum, cdf1, f1, f2, freq;
	memset(CDF2, 0, 257*sizeof(unsigned));
	if(p->mtop)
		t25_addhist(buf2, iw, ih, kc, x1-p->mleft, x2+p->mright, y-p->mtop, y, x1, x2, y, p->maxinc, CDF2);
	if(p->mleft)
		t25_addhist(buf2, iw, ih, kc, x1-p->mleft, x1, y, y+1, x1, x2, y, p->maxinc, CDF2);

	if(CDF2[256])
	{
		sum=0;
		for(int sym=0;sym<256;++sym)
		{
			cdf1=!overflow?CDF0[sym]:0x10000;
			if(sym<255)
				overflow|=cdf1>CDF0[sym+1];
			f1=(sym<255&&!overflow?CDF0[sym+1]:0x10000)-cdf1;

			f2=(int)(((long long)CDF2[sym]<<16)/CDF2[256]);//normalize

			freq=f1+(int)(((long long)f2-f1)*p->alpha/0xFF);//blend

			freq=(int)((long long)freq*0xFF00>>16)+1;//guard
			//freq=CLAMP(0, freq, 0xFF01);
			if(freq<0||freq>0xFF01)
			{
				printf("Impossible freq 0x%04X  f1 0x%04X  f2 0x%04X  W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d\n", freq, f1, f2, p->gwidth, p->mleft, p->mtop, p->mright, p->alpha, p->maxinc);
				return 0;
			}
				//LOG_ERROR("Impossible freq 0x%04X / 0x10000", freq);
			CDF2[sym]=sum;
			sum+=freq;
			if(sum>0x10000)
			{
				printf("ANS CDF sum 0x%04X, freq 0x%04X", sum, freq);
				return 0;
			}
		}
	}
	else
	{
		for(int sym=0;sym<256;++sym)
		{
			if(overflow)
				CDF2[sym]=0xFF00|sym;
			else
			{
				int cdf=CDF0[sym];
				CDF2[sym]=((unsigned)(cdf*0xFF00)>>16)+sym;
				if(sym<255)
					overflow|=cdf>CDF0[sym+1];
			}
		}
	}
	CDF2[256]=0x10000;
	return 1;
}
double t25_calcloss(const unsigned char *buf, int iw, int ih, int kc, Rect const *r, T25Params const *param, const unsigned short *CDF0, unsigned *CDF2, int loud)
{
	double chsize=0;
	if(!BETWEEN(1, param->gwidth, t25_limits.gwidth)||!BETWEEN(0, param->mleft, t25_limits.mleft)||!BETWEEN(0, param->mtop, t25_limits.mtop)||!BETWEEN(0, param->mright, t25_limits.mright)||!BETWEEN(0, param->alpha, t25_limits.alpha)||!BETWEEN(0, param->maxinc, t25_limits.maxinc))
		LOG_ERROR("Invalid params LOSS W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d", param->gwidth, param->mleft, param->mtop, param->mright, param->alpha, param->maxinc);
	//if(loud)
	//	printf("W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d\r",
	//		param->gwidth,
	//		param->mleft,
	//		param->mtop,
	//		param->mright,
	//		param->alpha,
	//		param->maxinc);
	for(int ky=r->y1;ky<r->y2;++ky)
	{
		for(int kx=r->x1;kx<r->x2;)
		{
			int xend=MINVAR(kx+param->gwidth, r->x2);
			int success=t25_prepblock(buf, CDF0, iw, ih, kc, kx, xend, ky, param, CDF2);
			if(!success)
				return 0;
				
			int kx2=kx;
			for(;kx2<kx+param->gwidth;++kx2)
			{
				unsigned char sym=buf[(iw*ky+kx2)<<2|kc];
				int freq=CDF2[sym+1]-CDF2[sym];
				double prob=(double)freq/0x10000, bitsize=-log2(prob);//Zipf's law
				chsize+=bitsize;
			}
			kx+=param->gwidth;
		}
	}
	++t25_ctr;
	return chsize;
}
#if 0
int t25_opt2(const unsigned char *buf, int iw, int ih, int kc, Rect const *r, const unsigned short *CDF0, unsigned *CDF2, T25Params *param, double *csize, int pidx, int step, int loud)
{
	int went_fwd=0;
	for(int subit=0;subit<20;++subit)
	{
		double csize0=*csize;
		short prevval=t25_incparam(param, pidx, step);
		if(prevval==T25_PARAM(param, pidx))//out of range
			break;
		if(loud)
			printf("W%c%3d  MLTR %c%3d %c%3d %c%3d  A%c0x%02X I%c%3d\r",
				pidx==0?'>':' ', param->gwidth,
				pidx==1?'>':' ', param->mleft,
				pidx==2?'>':' ', param->mtop,
				pidx==3?'>':' ', param->mright,
				pidx==4?'>':' ', param->alpha,
				pidx==5?'>':' ', param->maxinc);
		*csize=t25_calcloss(buf, iw, ih, kc, r, param, CDF0, CDF2);
		if(*csize>csize0)//cancel last change and break
		{
			T25_PARAM(param, pidx)=prevval;
			*csize=csize0;
			break;
		}
		went_fwd=1;
	}
	return went_fwd;
}
double t25_optimize(const unsigned char *buf, int iw, int ih, int kc, Rect const *r, T25Params *param, const unsigned short *CDF0, unsigned *CDF2, int loud)
{
	double csize00, csize;
	int prevval0;
	int steps[]={16, 8, 4, 2, 1};

	for(int ks=0;ks<COUNTOF(steps);++ks)
	{
		int step=steps[ks];
		for(int it=0, improve=1;it<64&&improve;++it)
		{
			improve=0;
			for(int pidx=0;pidx<sizeof(T25Params)/sizeof(short);++pidx)
			{
				prevval0=T25_PARAM(param, pidx);
				csize00=csize=t25_calcloss(buf, iw, ih, kc, r, param, CDF0, CDF2);

				int went_fwd=t25_opt2(buf, iw, ih, kc, r, CDF0, CDF2, param, &csize, pidx, step, loud);
				if(!went_fwd)
					t25_opt2(buf, iw, ih, kc, r, CDF0, CDF2, param, &csize, pidx, -step, loud);

				if(csize>csize00)//prevent CR from worsening
				{
					T25_PARAM(param, pidx)=prevval0;
					csize=csize00;
				}
			}
		}
	}
	return csize;
}
#endif
double t25_optimize_v2(const unsigned char *buf, int iw, int ih, int kc, Rect const *r, T25Params *param, const unsigned short *CDF0, unsigned *CDF2, int loud)
{
	double csize0;
	int steps[]={16, 8, 4, 2, 1};
	//	limit[]={ 1,  1, 1, 1, 1, 1};
	//int steps[]={64, 16, 4, 1},
	//	limit[]={4, 4, 4, 4};
	
	csize0=t25_calcloss(buf, iw, ih, kc, r, param, CDF0, CDF2, loud);
	if(!csize0)
		LOG_ERROR("Start W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d", param->gwidth, param->mleft, param->mtop, param->mright, param->alpha, param->maxinc);
	for(int ks=0;ks<COUNTOF(steps);++ks)
	{
		int step=steps[ks];
		int bestpidx=0, beststep=0;
		double bestcsize=csize0;
			
		for(int pidx=0;pidx<sizeof(T25Params)/sizeof(short);++pidx)
		{
			double csize;
			short prevval;

			prevval=t25_incparam(param, pidx, step);
			if(T25_PARAM(param, pidx)!=prevval)
			{
				csize=t25_calcloss(buf, iw, ih, kc, r, param, CDF0, CDF2, loud);
				if(!csize0)
					LOG_ERROR("Plus W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d", param->gwidth, param->mleft, param->mtop, param->mright, param->alpha, param->maxinc);
				T25_PARAM(param, pidx)=prevval;
				if(bestcsize>csize)
					bestcsize=csize, bestpidx=pidx, beststep=step;
			}

			prevval=t25_incparam(param, pidx, -step);
			if(T25_PARAM(param, pidx)!=prevval)
			{
				csize=t25_calcloss(buf, iw, ih, kc, r, param, CDF0, CDF2, loud);
				if(!csize0)
					LOG_ERROR("Minus W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d", param->gwidth, param->mleft, param->mtop, param->mright, param->alpha, param->maxinc);
				T25_PARAM(param, pidx)=prevval;
				if(bestcsize>csize)
					bestcsize=csize, bestpidx=pidx, beststep=-step;
			}
		}
		if(bestcsize<csize0)
		{
			t25_incparam(param, bestpidx, beststep);
			csize0=bestcsize;
		}
	}
	return csize0;
}
#if 0
int t25_optimizeall(const unsigned char *buf, int iw, int ih, int x1, int x2, int y1, int y2, int loud)
{
	int res=iw*ih;
	unsigned short *CDF0=(unsigned short*)malloc(256LL*3*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	for(int kc=0;kc<3;++kc)
	{
		memset(CDF2, 0, 256LL*sizeof(unsigned));
		for(int k=0;k<res;++k)
		{
			unsigned char sym=buf[k<<2|kc];
			++CDF2[sym];
		}
		t25_normalize_histogram(CDF2, 256, res, CDF0+((size_t)kc<<8));
	}
	
	if(x1<0)
		x1=0;
	if(x2>iw)
		x2=iw;
	if(y1<0)
		y1=0;
	if(y2>ih)
		y2=ih;
	
	double csizes[3]={0};
	int steps[]={4, 2, 1};
	int usize=(x2-x1)*(y2-y1);
	for(int kc=0;kc<3;++kc)
	{
		for(int ks=0;ks<COUNTOF(steps);++ks)
		{
			for(int it=0, improve=1;it<64&&improve;++it)
			{
				improve=0;
				for(int pidx=0;pidx<sizeof(t25_params)/sizeof(t25_params->gwidth);++pidx)
				{
					csizes[kc]=t25_optimize(buf, iw, ih, kc, x1, x2, y1, y2, t25_params+kc, pidx, steps[ks], CDF0+((size_t)kc<<8), CDF2);
					if(loud)
						io_render();
				}
			}
		}
		csizes[kc]/=8;
		t25_cr[kc]=usize>0&&csizes[kc]?(x2-x1)*(y2-y1)/csizes[kc]:0;
	}

	free(CDF0);
	free(CDF2);
	return 1;
}
#endif

void t25_calchist(const unsigned char *buf, int iw, int ih, int kc, int x1, int x2, int y1, int y2, unsigned *hist)
{
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			int idx=iw*ky+kx;
			unsigned char sym=buf[idx<<2|kc];
			++hist[sym];
		}
	}
}
double t25_calccsize(const unsigned char *buf, int iw, int ih, int kc, Rect const *r, unsigned *hist)
{
	memset(hist, 0, 256*sizeof(int));
	t25_calchist(buf, iw, ih, kc, r->x1, r->x2, r->y1, r->y2, hist);
	//for(int ky=r->y1;ky<r->y2;++ky)
	//{
	//	for(int kx=r->x1;kx<r->x2;++kx)
	//	{
	//		int idx=iw*ky+kx;
	//		unsigned char sym=buf[idx<<2|kc];
	//		++hist[sym];
	//	}
	//}
	int count=(r->x2-r->x1)*(r->y2-r->y1);
	double bitsize=0;
	for(int ky=r->y1;ky<r->y2;++ky)
	{
		for(int kx=r->x1;kx<r->x2;++kx)
		{
			int idx=iw*ky+kx;
			unsigned char sym=buf[idx<<2|kc];
			unsigned freq=hist[sym];
			double prob=(double)freq/count;
			bitsize-=log2(prob);//Zipf's law
		}
	}
	double csize=bitsize/8;
	return csize;
}
void t25_selectsmallblock(const unsigned char *buf, int iw, int ih, int kc, Rect const *lb, Rect *sb, int sbw, int sbh, unsigned *hist)
{
	int bw=lb->x2-lb->x1,
		bh=lb->y2-lb->y1,
		nbx=bw/sbw, nby=bh/sbh;
	if(nbx<=1||nby<=1)
	{
		*sb=*lb;
		return;
	}
	//memset(hist, 0, 256*sizeof(int));
	//for(int ky=lb->y1;ky<lb->y2;++ky)
	//{
	//	for(int kx=lb->x1;kx<lb->x2;++kx)
	//	{
	//		int idx=iw*ky+kx;
	//		unsigned char sym=buf[idx<<2|kc];
	//		++hist[sym];
	//	}
	//}
	double bestcsize=0;
	int bestx=0, besty=0;
	int it=0;
	for(int by=0;by<nby;++by)
	{
		int ky=lb->y1+by*sbh;
		for(int bx=0;bx<nbx;++bx, ++it)
		{
			int kx=lb->x1+bx*sbw;
			Rect r_cand={kx, kx+sbw, ky, ky+sbh};
			double csize=t25_calccsize(buf, iw, ih, kc, &r_cand, hist);
			if(!it||bestcsize>csize)
				bestcsize=csize, bestx=kx, besty=ky;
		}
	}
	sb->x1=bestx;
	sb->x2=bestx+sbw;
	sb->y1=besty;
	sb->y2=besty+sbh;
}

static T25Params t25_params[3]=
{
	{ 8, 26, 26, 26, 0xD4, 32},
	{23, 32, 32, 32, 0xD4, 32},
	{ 8, 26, 26, 26, 0xD4, 32},
};
int t25_encode(const unsigned char *src, int iw, int ih, int *blockw, int *blockh, ArrayHandle *data, int loud)
{
	int res=iw*ih;
	unsigned char *buf2=(unsigned char*)malloc((size_t)res<<2);
	unsigned short *CDF0=(unsigned short*)malloc(768LL*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!buf2||!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(buf2, src, (size_t)res<<2);
#if 0
	apply_transforms_fwd(buf2, iw, ih);
#else
	addbuf(buf2, iw, ih, 3, 4, 128);
	colortransform_ycocb_fwd((char*)buf2, iw, ih);
	pred_jxl_opt_v2((char*)buf2, iw, ih, jxlparams_i16, loud);
#if 1
	if(loud)
	{
		for(int kc=0;kc<3;++kc)//
		{
			for(int kp=0;kp<11;++kp)
			{
				short val=jxlparams_i16[11*kc+kp];
				printf(" %c0x%04X,", val<0?'-':' ', abs(val));
			}
			printf("\n");
		}
	}
#endif
	pred_jxl_apply((char*)buf2, iw, ih, jxlparams_i16, 1);
	addbuf(buf2, iw, ih, 3, 4, 128);
#endif

	for(int kc=0;kc<3;++kc)
	{
		memset(CDF2, 0, 256LL*sizeof(unsigned));
		t25_calchist(buf2, iw, ih, kc, 0, iw, 0, ih, CDF2);
		t21_normalize_histogram(CDF2, 256, res, CDF0+((size_t)kc<<8));
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);

	int ansbookmarks[3]={0};
	dlist_push_back(&list, 0, 12);
	dlist_push_back(&list, jxlparams_i16, 33*sizeof(short));
	dlist_push_back(&list, CDF0, 768*sizeof(short));

	int overhead[4]={0, 0, 0, (int)list.nobj};//

	for(int kc=0;kc<3;++kc)
	{
		int lbw=blockw[kc], lbh=blockh[kc],
		//	sbw=sbsizes[kc<<1], sbh=sbsizes[kc<<1|1],
			lbx=(iw+lbw-1)/lbw, lby=(ih+lbh-1)/lbh;
		Rect lblock;
		//Rect sblock;
		ArrayHandle params;
		ARRAY_ALLOC(T25ParamsPacked, params, 0, 0, (size_t)lbx*lby, 0);
		for(int by=0;by<lby;++by)
		{
			lblock.y1=by*lbh;
			lblock.y2=MINVAR(lblock.y1+lbh, ih);
			for(int bx=0;bx<lbx;++bx)
			{
				lblock.x1=bx*lbw;
				lblock.x2=MINVAR(lblock.x1+lbw, iw);

				t25_ctr=0;

#if 0
				lblock.x1=640;//
				lblock.x2=lblock.x1+lbw;
				lblock.y1=384;
				lblock.y2=lblock.y1+lbh;
				t25_params[kc].gwidth=1;
				t25_params[kc].mleft=30;
				t25_params[kc].mtop=6;
				t25_params[kc].mright=32;
				t25_params[kc].alpha=0xD7;
				t25_params[kc].maxinc=32;
#endif

				t25_optimize_v2(buf2, iw, ih, kc, &lblock, t25_params+kc, CDF0, CDF2, loud);//optimize for whole large block

				//t25_selectsmallblock(buf2, iw, ih, kc, &lblock, &sblock, sbw, sbh, CDF2);
				//t25_optimize(buf2, iw, ih, kc, &sblock, t25_params+kc, CDF0, CDF2);//optimize for small block

				T25ParamsPacked *pp=(T25ParamsPacked*)ARRAY_APPEND(params, 0, 1, 1, 0);
				for(int k=0;k<sizeof(T25ParamsPacked);++k)
					((unsigned char*)pp)[k]=(unsigned char)T25_PARAM(t25_params+kc, k);
				//pp->gwidth=(unsigned char)t25_params[kc].gwidth;
				//pp->mleft =(unsigned char)t25_params[kc].mleft;
				//pp->mtop  =(unsigned char)t25_params[kc].mtop;
				//pp->mright=(unsigned char)t25_params[kc].mright;
				//pp->alpha =(unsigned char)t25_params[kc].alpha;
				//pp->maxinc=(unsigned char)t25_params[kc].maxinc;

				if(loud)
					printf("CXY %d %4d %4d  W %3d  MLTR %3d %3d %3d  A 0x%02X I %3d  %d iters\r", kc, lblock.x1, lblock.y1, t25_params[kc].gwidth, t25_params[kc].mleft, t25_params[kc].mtop, t25_params[kc].mright, t25_params[kc].alpha, t25_params[kc].maxinc, t25_ctr);
			}
		}
		if(loud)
			printf("\n");
		dlist_push_back(&list, params->data, params->count*params->esize);

		overhead[kc]=(int)(params->count*params->esize);//

		//if(kc==2)//
		//	kc=2;

		unsigned state=0x10000;
		for(int ky=ih-1;ky>=0;--ky)
		{
			int by=ky/lbh;
			lblock.y1=by*lbh;
			lblock.y2=MINVAR(lblock.y1+lbh, ih);
			for(int bx=lbx-1;bx>=0;--bx)
			{
				T25ParamsPacked *pp=(T25ParamsPacked*)array_at(&params, lbx*by+bx);
				T25Params param={pp->gwidth, pp->mleft, pp->mtop, pp->mright, pp->alpha, pp->maxinc};
				int ng=(lbw+param.gwidth-1)/param.gwidth;
				lblock.x1=bx*lbw;
				lblock.x2=MINVAR(lblock.x1+lbw, iw);
				for(int kg=ng-1;kg>=0;--kg)
				{
					int x1=lblock.x1+kg*param.gwidth, x2=MINVAR(x1+param.gwidth, lblock.x2);
					int success=t25_prepblock(buf2, CDF0, iw, ih, kc, x1, x2, ky, &param, CDF2);
					if(!success)
						LOG_ERROR("t25_prepblock error");
					for(int kx=x2-1;kx>=x1;--kx)
					{
						unsigned char sym=buf2[(iw*ky+kx)<<2|kc];

						int cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;

						//if(kc==0&&ky==511&&kx==512)//
						//	printf("CXY %d %d %d  sym 0x%02X cdf 0x%04X freq 0x%04X  state 0x%08X\n", kc, kx, ky, sym, cdf, freq, state);

						if(!freq)
							LOG_ERROR("ZPS");

						//double prob=freq/65536.;
						//csize-=log2(prob);
						//unsigned s0=state;
						
						if(state>=(unsigned)(freq<<16))//renorm
						{
							dlist_push_back(&list, &state, 2);
							state>>=16;
						}
						debug_enc_update(state, cdf, freq, kx, ky, 0, kc, sym);
						state=state/freq<<16|(cdf+state%freq);//update
					}
				}
			}
		}
		dlist_push_back(&list, &state, 4);
		ansbookmarks[kc]=(int)list.nobj;

		array_free(&params);
	}
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, ansbookmarks, 12);
	
	int chsizes[]=
	{
		ansbookmarks[0]-overhead[3]    -overhead[0],
		ansbookmarks[1]-ansbookmarks[0]-overhead[1],
		ansbookmarks[2]-ansbookmarks[1]-overhead[2],
	};
	if(loud)
	{
		int totaloverhead=overhead[0]+overhead[1]+overhead[2]+overhead[3], totalch=chsizes[0]+chsizes[1]+chsizes[2];
		printf("Total    %7d  %lf\n", totaloverhead+totalch, 3.*res/list.nobj);
		printf("Overhead %7d\n", totaloverhead);
		printf("Red      %7d  %lf\n", chsizes[0], (double)res/chsizes[0]);
		printf("Green    %7d  %lf\n", chsizes[1], (double)res/chsizes[1]);
		printf("Blue     %7d  %lf\n", chsizes[2], (double)res/chsizes[2]);
	}

	dlist_clear(&list);
	free(buf2);
	free(CDF0);
	free(CDF2);
	return 1;
}
int t25_decode(const unsigned char *data, size_t srclen, int iw, int ih, int *blockw, int *blockh, unsigned char *buf, int loud)
{
	const int cdflen=768LL*sizeof(short), overhead=12LL+33*sizeof(short)+cdflen;
	int res=iw*ih;
	
	const unsigned char *srcptr, *srcstart, *srcend=data+srclen;
	if(data+overhead>=srcend)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}

	unsigned ansbookmarks[3];
	memcpy(ansbookmarks, data, 12);
	if(ansbookmarks[2]<(unsigned)overhead||ansbookmarks[2]>srclen)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	
	unsigned short *CDF0=(unsigned short*)malloc(cdflen);
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	short jxlparams[33];
	memcpy(jxlparams, data+12, 33*sizeof(short));
	memcpy(CDF0, data+12+33*sizeof(short), cdflen);
	int lbx[]=
	{
		(iw+blockw[0]-1)/blockw[0],
		(iw+blockw[1]-1)/blockw[1],
		(iw+blockw[2]-1)/blockw[2],
	};
	int lby[]=
	{
		(ih+blockh[0]-1)/blockh[0],
		(ih+blockh[1]-1)/blockh[1],
		(ih+blockh[2]-1)/blockh[2],
	};
	T25ParamsPacked *params[]=
	{
		(T25ParamsPacked*)(data+overhead),
		(T25ParamsPacked*)(data+ansbookmarks[0]),
		(T25ParamsPacked*)(data+ansbookmarks[1]),
	};
	int pcount[]=
	{
		lbx[0]*lby[0],
		lbx[1]*lby[1],
		lbx[2]*lby[2],
	};
	Rect block;
	for(int kc=0;kc<3;++kc)
	{
		unsigned state;
		srcptr=data+ansbookmarks[kc];
		srcstart=kc?data+ansbookmarks[kc-1]:data+overhead;
		srcptr-=4;
		if(srcptr<srcstart)
			LOG_ERROR("ANS buffer overflow");
		memcpy(&state, srcptr, 4);

		for(int ky=0;ky<ih;++ky)
		{
			int by=ky/blockh[kc];
			block.y1=by*blockh[kc];
			block.y2=MINVAR(block.y1+blockh[kc], ih);
			for(int bx=0;bx<lbx[kc];++bx)
			{
				T25ParamsPacked *pp=params[kc]+lbx[kc]*by+bx;
				T25Params param={pp->gwidth, pp->mleft, pp->mtop, pp->mright, pp->alpha, pp->maxinc};
				int ng=(blockw[kc]+param.gwidth-1)/param.gwidth;
				block.x1=bx*blockw[kc];
				block.x2=MINVAR(block.x1+blockw[kc], iw);
				for(int kg=0;kg<ng;++kg)
				{
					int x1=block.x1+kg*param.gwidth, x2=MINVAR(x1+param.gwidth, block.x2);
					int success=t25_prepblock(buf, CDF0, iw, ih, kc, x1, x2, ky, &param, CDF2);
					if(!success)
						LOG_ERROR("t25_prepblock error");
					for(int kx=x1;kx<x2;++kx)
					{
						unsigned c=(unsigned short)state;
						int sym=0;
					
						//if(kc==0&&ky==2&&kx==512)//
						//	printf("");

						int L=0, R=256, found=0;
						while(L<=R)
						{
							sym=(L+R)>>1;
							if(CDF2[sym]<c)
								L=sym+1;
							else if(CDF2[sym]>c)
								R=sym-1;
							else
							{
								found=1;
								break;
							}
						}
						if(!found)
							sym=L+(L<256&&CDF2[L]<c)-1;
						else
							for(;sym<256-1&&CDF2[sym+1]==c;++sym);

						buf[(iw*ky+kx)<<2|kc]=(unsigned char)sym;

						unsigned cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;
						
						debug_dec_update(state, cdf, freq, kx, ky, 0, kc, sym);
						state=freq*(state>>16)+c-cdf;//update
						if(state<0x10000)//renorm
						{
							state<<=16;
							if(srcptr-2>=srcstart)
							{
								srcptr-=2;
								memcpy(&state, srcptr, 2);
							}
						}
					}
				}
			}
		}
	}
	free(CDF0);
	free(CDF2);

	addbuf(buf, iw, ih, 3, 4, 128);
	pred_jxl_apply((char*)buf, iw, ih, jxlparams, 0);
	colortransform_ycocb_inv((char*)buf, iw, ih);
	addbuf(buf, iw, ih, 3, 4, 128);

	for(int k=0;k<res;++k)//set alpha
		buf[k<<2|3]=0xFF;
	return 1;
}