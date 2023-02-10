#include"battle.h"
#include<stdio.h>//for debugging
#include<string.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#include"rans_common.h"
static const char file[]=__FILE__;

	#define PROF(...)//

const int tag_rans4='A'|'N'<<8|'0'<<16|'4'<<24;

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
		memcpy(hist+(k<<1), integrate?&sum:&h[k].qfreq, 2);//2-byte alignment
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
				p->shift=0;
				p->inv_freq=0xFFFFFFFF;
				p->bias=sum+ANS_L-1;
			}
			else
			{
				p->shift=ceil_log2(p->freq)-1;
				p->inv_freq=(unsigned)(((0x100000000<<p->shift)+p->freq-1)/p->freq);
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

int rans4_encode(const void *src, ptrdiff_t nbytes, int symbytes, int is_signed, ArrayHandle *out, int loud)//symbytes: up to 16
{
	const int infosize=ANS_NLEVELS*sizeof(SymbolInfo), lginfosize=13;
	DList list;
	const unsigned char *buf=(const unsigned char*)src;
	size_t dstidx;
	SymbolInfo *info;
	int internalheadersize=4+8+symbytes*ANS_NLEVELS*sizeof(short);
	int chmask=symbytes-1;

	if(*out)
	{
		if(out[0]->esize!=1)
			return RANS_INVALID_DST;
		dstidx=out[0]->count;
		ARRAY_APPEND(*out, 0, internalheadersize, 1, 0);
	}
	else
	{
		dstidx=0;
		ARRAY_ALLOC(char, *out, 0, internalheadersize, 0, 0);
	}
	memcpy(out[0]->data+dstidx, &tag_rans4, 4);
	dstidx+=4+8;
	for(int kc=0;kc<symbytes;++kc)
		rans_calc_histogram(buf+kc, (int)(nbytes/symbytes), symbytes, out[0]->data+dstidx+kc*ANS_NLEVELS*sizeof(short), ANS_PROB_BITS, 0);
	rans_prep(out[0]->data+dstidx, symbytes, &info, 0, 0);
	dstidx-=8;
	dlist_init(&list, 1, 1024, 0);

	long long t1=__rdtsc();

	unsigned state[16]={0};
	for(int kc=0;kc<symbytes;++kc)
		state[kc]=ANS_L;
	for(ptrdiff_t ks=nbytes-symbytes, idx=nbytes-1;ks>=0;ks-=symbytes)
	{
		for(int kc=symbytes-1;kc>=0;--kc, --idx)
		{
			unsigned char val=buf[idx];
			if(is_signed)
			{
				int neg=val<0;
				val^=-neg;
				val+=neg;
				val<<=1;
				val|=neg;
			}
			SymbolInfo *p=info+(kc<<ANS_DEPTH|val);

			//renormalize
			if(state[kc]>=p->renorm_limit)
		//	if(state>=(unsigned)(si.freq<<(32-ANS_PROB_BITS)))
			{
				dlist_push_back(&list, state+kc, 2);
				state[kc]>>=16;
			}
			PROF(RENORM);

			//update
#ifdef ANS_PRINT_STATE2
			printf("enc: 0x%08X = 0x%08X+(0x%08X*0x%08X>>(32+%d))*0x%04X+0x%08X\n", state+((unsigned)((long long)state*si.inv_freq>>32)>>si.shift)*si.neg_freq+si.bias, state, state, si.inv_freq, si.shift, si.neg_freq, si.bias);
#endif
			state[kc]+=(((long long)state[kc]*p->inv_freq>>32)>>p->shift)*p->neg_freq+p->bias;//Ryg's division-free rANS encoder	https://github.com/rygorous/ryg_rans/blob/master/rans_byte.h
			//state=(state/p->freq<<ANS_PROB_BITS)+state%p->freq+p->CDF;
			PROF(UPDATE);
		}
	}
	dlist_push_back(&list, state, symbytes*sizeof(unsigned));

	memcpy(out[0]->data+dstidx, &list.nobj, sizeof(size_t));
	dlist_appendtoarray(&list, out);

	t1=__rdtsc()-t1;

	free(info);
	dlist_clear(&list);
	return RANS_SUCCESS;
}

static int decode_error(size_t p, size_t srcstart, ptrdiff_t ks)
{
	printf("Buffer underrun p<srcstart, %p < %p at %d\n", (void*)p, (void*)srcstart, (int)ks);
	return 0;
}
#define READ_GUARD(NBYTES, IDX)		if((srcptr-=NBYTES)<srcstart)return decode_error((size_t)srcptr, (size_t)srcstart, IDX)
int rans4_decode(const unsigned char *srcdata, ptrdiff_t srclen, ptrdiff_t nbytes, int symbytes, int is_signed, void *dstbuf, int loud)
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
	if(memcmp(srcptr, &tag_rans4, 4))
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