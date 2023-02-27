#include"battle.h"
#include<stdio.h>//for debugging
#include<string.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#include"rans_common.h"
#define _USE_MATH_DEFINES
#include<math.h>//experimental double precision
static const char file[]=__FILE__;

	#define PROF(...)//

typedef struct SortedHistInfoStruct
{
	int	sym,  //symbol
		freq, //original freq
		qfreq;//quantized freq
} SortedHistInfo;
static int histinfo_byfreq(const void *left, const void *right)
{
	SortedHistInfo const *a, *b;

	a=(SortedHistInfo const*)left;
	b=(SortedHistInfo const*)right;
	return (a->freq>b->freq)-(a->freq<b->freq);
}
static int histinfo_bysym(const void *left, const void *right)
{
	SortedHistInfo const *a, *b;

	a=(SortedHistInfo const*)left;
	b=(SortedHistInfo const*)right;
	return (a->sym>b->sym)-(a->sym<b->sym);
}
static int rans_calc_histogram(const unsigned char *buffer, int nsymbols, int bytestride, unsigned char *hist, int prob_bits, int integrate)//hist is unsigned char due to alignment issues, but it's 16bit
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
	if(buffer)
	{
		for(int k=0;k<bytesize;k+=bytestride)//this loop takes 73% of encode time
			++h[buffer[k]].freq;
	}
	else
	{
		nsymbols=0;
		for(int k=0;k<ANS_NLEVELS;++k)
		{
			h[k].freq=((unsigned short*)hist)[k];
			nsymbols+=h[k].freq;
		}
		if(!nsymbols)
			LOG_ERROR("Invalid predictor, sum = %d", nsymbols);
	}
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
	if(!buffer)
	{
		for(int k=0;k<ANS_NLEVELS;++k)
			((unsigned short*)hist)[k]=h[k].qfreq;
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

static int rans_prep(const void *hist_ptr, int bytespersymbol, SymbolInfo **info, unsigned char **CDF2sym, int loud)
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

static const int tag_rans4='A'|'N'<<8|'0'<<16|'4'<<24;
int rans4_encode(const void *src, ptrdiff_t nbytes, int symbytes, int is_signed, ArrayHandle *out, int loud, unsigned short *custom_pred)//symbytes: up to 16
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
	dstidx+=(size_t)4+8;
	unsigned char *histptr;
	if(custom_pred)
		histptr=(unsigned char*)custom_pred;
	else
		histptr=out[0]->data+dstidx;
	for(int kc=0;kc<symbytes;++kc)
		rans_calc_histogram(custom_pred?0:buf+kc, (int)(nbytes/symbytes), symbytes, histptr+kc*ANS_NLEVELS*sizeof(short), ANS_PROB_BITS, 0);
	rans_prep(histptr, symbytes, &info, 0, 0);
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
int rans4_decode(const unsigned char *srcdata, ptrdiff_t srclen, ptrdiff_t nbytes, int symbytes, int is_signed, void *dstbuf, int loud, unsigned short *custom_pred)
{
	const int
		histsize=ANS_NLEVELS*sizeof(short), lghistsize=9,
		infosize=ANS_NLEVELS*sizeof(SymbolInfo), lginfosize=13;
	const unsigned char
		*data=(const unsigned char*)srcdata,
		*srcptr=data,
		*srcstart;
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

	const void *histptr;
	if(custom_pred)
		histptr=custom_pred;
	else
		histptr=srcptr;
	rans_prep(histptr, symbytes, &info, &CDF2sym, 0);
	srcstart=srcptr+((size_t)symbytes<<9);
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


const int tag_rans8='A'|'N'<<8|'0'<<16|'8'<<24;
//binary rANS
int rans8_testencode(const void *src, int bw, int bh, int bitdepth, int bytestride, unsigned short prob_MPS, ArrayHandle *out)
{
	size_t dstidx;
	if(*out)
	{
		if(out[0]->esize!=1)
			return RANS_INVALID_DST;
		dstidx=out[0]->count;
	}
	else
	{
		dstidx=0;
		ARRAY_ALLOC(char, *out, 0, 0, 0, 0);
	}
	long long t1=__rdtsc();

	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, 12);//{tag, csize}

	const unsigned char *buf=(const unsigned char*)src;
	unsigned state=0x10000;
	if(!prob_MPS)
		LOG_ERROR("Binary rANS: prob_MPS 0x%04X", prob_MPS);
	//unsigned short prob_MPS=0xC000;//0.75
	//unsigned short *prob=(unsigned short*)malloc(bitdepth*sizeof(short));
	//for(int k=0;k<bitdepth;++k)
	//	prob[k]=0x8000;
	int bitstride=bytestride<<3;
	ptrdiff_t bitidx=((ptrdiff_t)(bh*bw-1)*bitstride)+bitdepth-1, rowbits=(ptrdiff_t)bw*bitstride;
	for(int ky=bh-1;ky>=0;--ky)
	{
		for(int kx=bw-1;kx>=0;--kx, bitidx-=bitstride)
		{
			ptrdiff_t bitidx2=bitidx;
			for(int kp=bitdepth-1;kp>=0;--kp, --bitidx2)
			{
				unsigned char bit=buf[bitidx2>>3]>>(bitidx2&7)&1, b0=0;
				ptrdiff_t bitidx0=bitidx2-rowbits;
				if(bitidx0>=0)
					b0=buf[bitidx0>>3]>>(bitidx0&7)&1;

				unsigned short p0=prob_MPS;
				p0^=-b0;//flip p0 to LPS if b0 is true
				p0+=b0;

				unsigned short CDF=p0&-bit, freq=(p0^-bit)+bit;
				if(state>=(unsigned)(freq<<16))
				{
					dlist_push_back(&list, &state, 2);
					state>>=16;
				}
				state=state/freq<<16|(CDF+state%freq);
			}
		}
	}
	dlist_appendtoarray(&list, out);
	memcpy(out[0]->data+dstidx, &tag_rans8, 4);
	memcpy(out[0]->data+dstidx+4, &list.nobj, 8);
	t1=__rdtsc()-t1;

	dlist_clear(&list);
	return RANS_SUCCESS;
}

static const int tag_rans0b='A'|'N'<<8|'0'<<16|'A'<<24;
//divides bitplanes into blocks and stores hamming weight of each block, prob_bits is in {2, 4, 8}
int rans0b_testencode(const void *src, int bw, int bh, int bitdepth, int bytestride, int blockdim, int prob_bits, ArrayHandle *out)
{
	size_t dstidx;
	if(*out)
	{
		if(out[0]->esize!=1)
			return RANS_INVALID_DST;
		dstidx=out[0]->count;
	}
	else
	{
		dstidx=0;
		ARRAY_ALLOC(char, *out, 0, 0, 0, 0);
	}
	long long t1=__rdtsc();

	if(bw%blockdim||bh%blockdim)
		LOG_ERROR("Invalid dimensions");

	int nbx=bw/blockdim, nby=bh/blockdim;
	int auxsize=(nbx*nby*bitdepth*prob_bits+7)>>3;
	unsigned char *probbuf=(unsigned char*)malloc(auxsize);
	memset(probbuf, 0, auxsize);

	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, 12);//{tag, csize}
	const unsigned char *buf=(const unsigned char*)src;
	unsigned state=0x10000;
	int blocksize=blockdim*blockdim;
	for(int kp=bitdepth-1;kp>=0;--kp)//for each bitplane
	{
		//for each block
		for(int ky=nby-1;ky>=0;--ky)
		{
			for(int kx=nbx-1;kx>=0;--kx)
			{
				int freq=blocksize;//count of zero bits

				//for each bit
				for(int ky2=blockdim-1;ky2>=0;--ky2)
				{
					for(int kx2=blockdim-1;kx2>=0;--kx2)
					{
						int bitidx=bytestride*(bw*(blockdim*ky+ky2)+blockdim*kx+kx2)+kp;//unoptimized
						unsigned char bit=buf[bitidx>>3]>>(bitidx&7)&1;
						freq-=bit;
					}
				}

				int qfreq=((freq<<prob_bits)+(blocksize>>1))/blocksize;
				qfreq+=!qfreq&&freq!=0;							//allow 0% and 100% prob only if the unquantized freq really is so
				qfreq-=qfreq==1<<prob_bits&&freq<blocksize;
				int blockidx=nbx*ky+kx;
				int bitidx2=blockidx*prob_bits;
				if(bitidx2>=auxsize<<3)//
					LOG_ERROR("OOB");//
				probbuf[bitidx2>>3]|=qfreq<<(bitidx2&7);

				if(qfreq&&qfreq<1<<prob_bits)//otherwise the state won't be changed
				{
					unsigned short p0=qfreq<<(16-prob_bits);

					//for each bit
					for(int ky2=blockdim-1;ky2>=0;--ky2)
					{
						for(int kx2=blockdim-1;kx2>=0;--kx2)
						{
							int bitidx=bytestride*(bw*(blockdim*ky+ky2)+blockdim*kx+kx2)+kp;
							if(bitidx>=bw*bh*bytestride<<3)//
								LOG_ERROR("OOB");//
							unsigned char bit=buf[bitidx>>3]>>(bitidx&7)&1;

							unsigned short CDF=p0&-bit, freq=(p0^-bit)+bit;
							if(state>=(unsigned)(freq<<16))
							{
								dlist_push_back(&list, &state, 2);
								state>>=16;
							}
							state=state/freq<<16|(CDF+state%freq);
						}
					}
				}
			}
		}
	}
	dlist_push_back(&list, probbuf, auxsize);
	dlist_appendtoarray(&list, out);
	memcpy(out[0]->data+dstidx, &tag_rans0b, 4);
	memcpy(out[0]->data+dstidx+4, &list.nobj, 8);
	t1=__rdtsc()-t1;

	dlist_clear(&list);
	return RANS_SUCCESS;
}

//2D adaptive rANS
size_t rans0c_testencode(const void *src, int bw, int bh, int symbytes, int bytestride, ArrayHandle *out)
{
	double invsqrt2=1/sqrt(2.), invsqrt2pi=1/sqrt(2*M_PI);
	size_t csize=4;
	const unsigned char *buf=(const unsigned char*)src;
	unsigned state=0x10000;
	for(int ky=bh-1;ky>=0;--ky)
	{
		for(int kx=bw-1;kx>=0;--kx)
		{
			for(int kc=symbytes-1;kc>=0;--kc)
			{
				//printf("%d %d %d\r", kx, ky, kc);
				unsigned char
					s00=buf[bytestride*(bw*ky+kx)+kc],
					s01=kx>0?buf[bytestride*(bw*ky+kx-1)+kc]:0,
					s10=ky>0?buf[bytestride*(bw*(ky-1)+kx)+kc]:0;
				double
					mean=(s01+s10)*0.5,
					conf=1./(abs(s01-s10)+3),
					a=(0-mean)*conf,		//truncated normal distribution
					b=(255-mean)*conf,
					x=(s00-mean)*conf,
					x1=(s00+1-mean)*conf,
					Phi_a=0.5*(1+erf(a*invsqrt2)),
					Phi_b=0.5*(1+erf(b*invsqrt2)),
					Phi_x=0.5*(1+erf(x*invsqrt2)),
					Phi_x1=0.5*(1+erf(x1*invsqrt2)),
					den=Phi_b-Phi_a,
					d_CDF=(Phi_x-Phi_a)/den,
					d_freq=(Phi_x1-Phi_x)/den;
				unsigned short freq=(int)(0xFFFF*d_freq)+1, CDF=(int)(0x10000*d_CDF);

				if(state>=(unsigned)(freq<<16))//renorm
				{
					csize+=2;
					state>>=16;
				}

				//if(!freq)
				//	LOG_ERROR("%d %d %d freq = %d", kx, ky, kc, freq);
				state=state/freq<<16|(CDF+state%freq);
			}
		}
	}
	csize+=4;
	return csize;
}

static const int tag_rans05='A'|'N'<<8|'0'<<16|'5'<<24;
typedef struct FreqInfoStruct
{
	double freq;//original freq
	unsigned
		sym,    //symbol
		qfreq;  //quantized freq
} FreqInfo;
static int freqinfo_byfreq(const void *left, const void *right)
{
	FreqInfo const *a, *b;

	a=(FreqInfo const*)left;
	b=(FreqInfo const*)right;
	return (a->freq>b->freq)-(a->freq<b->freq);
}
static int freqinfo_bysym(const void *left, const void *right)
{
	FreqInfo const *a, *b;

	a=(FreqInfo const*)left;
	b=(FreqInfo const*)right;
	return (a->sym>b->sym)-(a->sym<b->sym);
}
static FreqInfo g_temp[256];
Prob* rans5_preptable()
{
	const int nmean=256, nconf=256, nlevels=256;
	const State prob_max=ONE-255;//ONE-(nlevels-1)
	size_t len=(size_t)nmean*nconf*nlevels;
	Prob *table=(Prob*)malloc(len*sizeof(Prob));
	if(!table)
	{
		LOG_ERROR("Couldn't allocate arANS lookup table");
		return 0;
	}
	//memset(table, 0, len*sizeof(int));
	double invsqrt2=1/sqrt(2.), invsqrt2pi=1/sqrt(2*M_PI);
	for(int diff=0;diff<nconf;++diff)//abs(left-top)
	{
		//printf("\r%3d / 256", diff);
		double conf=1./((diff<<8)/nconf+3);
		for(int sum=0;sum<nmean;++sum)//left+top
		{
			int mean=(sum<<8)/nmean;
			double dsum=0;
			for(int ks=0;ks<nlevels;++ks)
			{
				double
					a=(0-mean)*conf,
					b=(255-mean)*conf,
					x=(ks-mean)*conf,
					Phi_a=0.5*(1+erf(a*invsqrt2)),
					Phi_b=0.5*(1+erf(b*invsqrt2)),
					Phi_x=0.5*(1+erf(x*invsqrt2)),
					den=Phi_b-Phi_a,
					d_freq=invsqrt2pi*exp(-0.5*x*x)*conf/den;//column width=1
					//d_CDF=(Phi_x-Phi_a)/den;
				//if(d_freq>1)
				//	LOG_ERROR("freq %lf", d_freq);
				g_temp[ks].freq=d_freq;
				g_temp[ks].sym=ks;
				g_temp[ks].qfreq=(unsigned)(ONE*d_freq);
				dsum+=d_freq;
				//if(!sum&&!diff)
				//	printf("%3d: %lf\n", ks, d_freq);
			}

			isort(g_temp, nlevels, sizeof(FreqInfo), freqinfo_byfreq);

			int idx;
			for(idx=0;idx<nlevels;++idx)//zero freq is not allowed for any symbol in lookup table
				g_temp[idx].qfreq+=!g_temp[idx].qfreq;
			//for(idx=nlevels-1;idx>=0&&g_temp[idx].qfreq>prob_max;--idx)
			//	g_temp[idx].qfreq=prob_max;

			long long error=-(long long)ONE;//too much -> +ve error & vice versa
			for(int k=0;k<nlevels;++k)
				error+=g_temp[k].qfreq;
			//printf("error = %lld\r", (long long)error);//
			if(error>0)
			{
				while(error)//decrement starting from smaller frequencies
				{
					for(int k=0;k<nlevels&&error;++k)
					{
						int dec=g_temp[k].qfreq>1;
						g_temp[k].qfreq-=dec, error-=dec;
					}
				}
			}
			else
			{
				while(error)//increment starting from high frequencies
				{
					for(int k=nlevels-1;k>=0&&error;--k)
					{
						int inc=g_temp[k].qfreq<prob_max;
						g_temp[k].qfreq+=inc, error+=inc;
					}
				}
			}

			isort(g_temp, nlevels, sizeof(FreqInfo), freqinfo_bysym);

			State fsum=0;
			idx=nlevels*(nmean*diff+sum);
			for(int ks=0;ks<nlevels;++ks)
			{
				Prob f=g_temp[ks].qfreq;
				table[idx|ks]=(Prob)fsum;
				fsum+=f;
				//if(!sum&&!diff)
				//	printf("%3d: 0x%04X\n", ks, f);
			}
			if(fsum!=ONE)
				LOG_ERROR("diff=%d, sum=%d: fsum = 0x%04llX", diff, sum, (long long)fsum);
		}
	}
	return table;
}
long long rans5_encode(const void *src, int bw, int bh, int symbytes, int bytestride, ArrayHandle *out, Prob *table)
{
	const int nmean=256, nconf=256, nlevels=256;
	
	long long cycles=__rdtsc();
	size_t dstidx;
	if(*out)
	{
		if(out[0]->esize!=1)
			return RANS_INVALID_DST;
		dstidx=out[0]->count;
	}
	else
		dstidx=0;

	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, 12);

	const unsigned char *buf=(const unsigned char*)src;
	State state=ONE;
	for(int ky=bh-1;ky>=0;--ky)
	{
		for(int kx=bw-1;kx>=0;--kx)
		{
			for(int kc=symbytes-1;kc>=0;--kc)
			{
				unsigned char
					s00=buf[bytestride*(bw*ky+kx)+kc],
					s01=kx>0?buf[bytestride*(bw*ky+kx-1)+kc]:0,
					s10=ky>0?buf[bytestride*(bw*(ky-1)+kx)+kc]:0;
				int sum=(s01+s10)>>1, diff=abs(s01-s10);
				int idx=nlevels*(nmean*diff+sum)+s00;
				Prob CDF=table[idx], freq=(Prob)((s00<255?(State)table[idx+1]:ONE)-CDF);
				
				if(state>=(State)((State)freq<<PROBBITS))//renorm
				{
					dlist_push_back(&list, &state, sizeof(Prob));
					state>>=PROBBITS;
				}

				state=state/freq<<PROBBITS|(state%freq+CDF);
			}
		}
	}
	dlist_push_back(&list, &state, sizeof(state));
	dlist_appendtoarray(&list, out);

	memcpy(out[0]->data+dstidx, &tag_rans05, 4);
	memcpy(out[0]->data+dstidx+4, &list.nobj, 8);
	dlist_clear(&list);
	cycles=__rdtsc()-cycles;
	return cycles;
}