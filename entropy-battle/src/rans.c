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
#include"lodepng.h"//for debugging
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

static int rans_prep(const void *hist_ptr, int symbytes, SymbolInfo **info, unsigned char **CDF2sym, int loud)
{
	int tempsize=symbytes*(ANS_NLEVELS*sizeof(SymbolInfo)+(ANS_L&-(CDF2sym!=0)));
	*info=(SymbolInfo*)malloc(tempsize);
	if(!*info)
		LOG_ERROR("Failed to allocate temp buffer");
	if(CDF2sym)
		*CDF2sym=(unsigned char*)*info+symbytes*ANS_NLEVELS*sizeof(SymbolInfo);
	for(int kc=0;kc<symbytes;++kc)
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
int rans4_encode(const void *src, ptrdiff_t nsymbols, int symbytes, int bytestride, ArrayHandle *out, unsigned short *custom_pred)//symbytes: up to 16
{
	const int infosize=ANS_NLEVELS*sizeof(SymbolInfo), lginfosize=13;
	DList list;
	const unsigned char *buf=(const unsigned char*)src;
	size_t dstidx;
	SymbolInfo *info;
	int internalheadersize=4+8+(int)((size_t)symbytes*ANS_NLEVELS*sizeof(short));

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
		rans_calc_histogram(custom_pred?0:buf+kc, (int)nsymbols, bytestride, histptr+kc*ANS_NLEVELS*sizeof(short), ANS_PROB_BITS, 0);
	rans_prep(histptr, symbytes, &info, 0, 0);
	dstidx-=8;
	dlist_init(&list, 1, 1024, 0);

	long long t1=__rdtsc();

	unsigned state=ANS_L;
	//unsigned state[16]={0};
	//for(int kc=0;kc<symbytes;++kc)
	//	state[kc]=ANS_L;
	for(ptrdiff_t ks=(nsymbols-1)*bytestride;ks>=0;ks-=bytestride)
	{
		for(int kc=symbytes-1;kc>=0;--kc)
		{
			unsigned char val=buf[ks+kc];
			//if(is_signed)
			//{
			//	int neg=val<0;
			//	val^=-neg;
			//	val+=neg;
			//	val<<=1;
			//	val|=neg;
			//}
			SymbolInfo *p=info+(kc<<ANS_DEPTH|val);

			//renormalize
			if(state>=p->renorm_limit)
		//	if(state>=(unsigned)(si.freq<<(32-ANS_PROB_BITS)))
			{
				dlist_push_back(&list, &state, 2);
				state>>=16;
			}
			PROF(RENORM);

			//update
#ifdef ANS_PRINT_STATE2
			printf("enc: 0x%08X = 0x%08X+(0x%08X*0x%08X>>(32+%d))*0x%04X+0x%08X\n", state+((unsigned)((long long)state*si.inv_freq>>32)>>si.shift)*si.neg_freq+si.bias, state, state, si.inv_freq, si.shift, si.neg_freq, si.bias);
#endif
			state+=(((long long)state*p->inv_freq>>32)>>p->shift)*p->neg_freq+p->bias;//Ryg's division-free rANS encoder	https://github.com/rygorous/ryg_rans/blob/master/rans_byte.h
			//state=(state/p->freq<<ANS_PROB_BITS)+state%p->freq+p->CDF;
			PROF(UPDATE);
		}
	}
	dlist_push_back(&list, &state, 4);

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
int rans4_decode(const unsigned char *srcdata, ptrdiff_t srclen, ptrdiff_t nsymbols, int symbytes, int bytestride, void *dstbuf, unsigned short *custom_pred)
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
	int internalheadersize=4+8+(int)((size_t)symbytes*ANS_NLEVELS*sizeof(short));
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

	unsigned state=0;
	READ_GUARD(4, 0);
	memcpy(&state, srcptr, sizeof(unsigned));
	//unsigned state[16]={0};
	//READ_GUARD(symbytes*sizeof(unsigned), 0);
	//memcpy(state, srcptr, symbytes*sizeof(unsigned));
	ptrdiff_t nbytes=nsymbols*bytestride;
	for(ptrdiff_t ks=0;ks<nbytes;ks+=bytestride)
	{
		for(int kc=0;kc<symbytes;++kc)
		{
			unsigned short c=(unsigned short)state;
			unsigned char val=CDF2sym[kc<<ANS_PROB_BITS|c];
			SymbolInfo *p=info+(kc<<ANS_DEPTH|val);
			//if(!p->freq)
			//	LOG_ERROR("Symbol 0x%02X has zero frequency", s);
			
			//if(is_signed)
			//{
			//	int neg=val&1;
			//	val>>=1;
			//	val^=-neg;
			//	val+=neg;
			//	val|=neg<<7&-!val;
			//}
			int idx=(int)ks+kc;
			pixels[ks+kc]=val;
			PROF(FETCH);

#ifdef ANS_PRINT_STATE2
			printf("dec: 0x%08X = 0x%04X*(0x%08X>>%d)+0x%04X-0x%08X\n", si.freq*(state>>ANS_PROB_BITS)+c-si.CDF, (int)si.freq, state, ANS_PROB_BITS, c, si.CDF);
#endif
			state=p->freq*(state>>ANS_PROB_BITS)+c-p->CDF;
			PROF(UPDATE);

			if(state<ANS_L)
			{
				if(idx>=nbytes-1)//shouldn't need this
					break;//
				READ_GUARD(2, idx);
				state<<=16;
				memcpy(&state, srcptr, 2);
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


#if 0
	#define ALLOC_BLOCK

//test4: 8 bit  3 channels  no alpha
static const int tag_test04='T'|'0'<<8|'0'<<16|'4'<<24;
#define BLOCKDIM 512
#ifndef ALLOC_BLOCK
static short tblock[BLOCKDIM*BLOCKDIM], trow[BLOCKDIM];
#endif
unsigned tCDF[2048];
static void debug_get_minmax(short *block, int bw, int bh, int *ret_minmax)
{
	int k, res=bw*bh;

	int vmin=0, vmax=0;
	for(k=0;k<res;++k)
	{
		if(!k||vmin>block[k])
			vmin=block[k];
		if(!k||vmax<block[k])
			vmax=block[k];
	}
	if(ret_minmax)
	{
		ret_minmax[0]=vmin;
		ret_minmax[1]=vmax;
	}
}
static void debug_print_CDF(int nlevels, unsigned normal)
{
	int mask=nlevels-1, vmax=0;
	for(int sym=0;sym<nlevels;++sym)
	{
		unsigned short cdf=tCDF[sym], freq;
		{
			int next;
			if(sym==mask)
				next=normal;
			else
			{
				next=tCDF[sym+1];
				//if(!allowZF&&next<cdf)//CDF will overflow with trailing zero frequency symbols  (uint16)0x10000==0
				//	next=normal;
			}
			freq=(unsigned short)(next-cdf);
//#ifdef _DEBUG
//			if(normal==0x10000&&!freq)
//				LOG_ERROR("rANS zero freq sym 0x%02X", sym);
//#endif
		}
		if(vmax<freq)
			vmax=freq;
	}
	for(int sym=0;sym<nlevels;++sym)
	{
		unsigned short cdf=tCDF[sym], freq;
		{
			int next;
			if(sym==mask)
				next=normal;
			else
			{
				next=tCDF[sym+1];
				//if(!allowZF&&next<cdf)//CDF will overflow with trailing zero frequency symbols  (uint16)0x10000==0
				//	next=normal;
			}
			freq=(unsigned short)(next-cdf);
//#ifdef _DEBUG
//			if(normal==0x10000&&!freq)
//				LOG_ERROR("rANS zero freq sym 0x%02X", sym);
//#endif
		}
		printf("%4d 0x%04X ", sym, freq);
		for(int k2=0, end=freq*64/vmax;k2<end;++k2)
			printf("*");
		printf("\n");
	}
	printf("\n");
}
static void debug_save_block(short *block, int bw, int bh, int offset)
{
	int k, res=bw*bh;

	for(k=0;k<res;++k)//swap bytes
	{
		unsigned short val=block[k]+offset;
		block[k]=val>>8|val<<8;
	}

	lodepng_encode_file("out.PNG", (unsigned char*)block, bw, bh, LCT_GREY, 16);//

	for(k=0;k<res;++k)//swap bytes back
	{
		unsigned short val=block[k];
		block[k]=(val>>8|val<<8)-offset;
	}
}
static int test4_encode4_prep(unsigned *CDF, int nlevels)
{
	unsigned prob_max=0x10000-(nlevels-1);
#if 0
	int half=nlevels>>1;
	double sdev=(double)ssdev*half/0x10000;
	double in_gain=1./(sqrt(2)*sdev), x1=0.5*(1+erf(-half*in_gain)), x2=0.5*(1+erf(half*in_gain)), invZ=1/(x2-x1);
	invZ*=0x10000;
	int error=0;
	for(int k=0;k<nlevels;++k)
	{
		double x=0.5*(1+erf((k+1-half)*in_gain));
		double f=(x-x1)*invZ;
		if(f<1)
			f=1;
		if(f>prob_max)
			f=prob_max;
		CDF[k]=(unsigned)f;
		error+=CDF[k];
		x1=x;
	}
#endif
	int error=0;
	for(int k=0;k<nlevels;++k)
		error+=CDF[k];

	error-=0x10000;
	if(error>0)
	{
		while(error)//decrement starting from smaller frequencies
		{
			for(int k=0;k<nlevels&&error;++k)
			{
				int dec=CDF[k]>1;
				CDF[k]-=dec, error-=dec;
			}
		}
	}
	else
	{
		while(error)//increment starting from high frequencies
		{
			for(int k=nlevels-1;k>=0&&error;--k)
			{
				int inc=CDF[k]<prob_max;
				CDF[k]+=inc, error+=inc;
			}
		}
	}
	int sum=0;
	for(int k=0;k<nlevels;++k)
	{
		unsigned freq=CDF[k];
		CDF[k]=sum;
		sum+=freq;
	}
	return sum;
}
static int test4_encode4_ans(const short *block, int dx, int dy, int rowstride, unsigned short nbits, DList *list)
{
	int nlevels=1<<nbits, vmax=nlevels-1, half=nlevels>>1;

	//static int oneshot=0;//DEBUG
#if 0
	if(!oneshot&&dx==128&&dy==128)//DEBUG
	{
		memset(tCDF, 0, nlevels*sizeof(unsigned));
		for(int ky=0;ky<dy;++ky)
		{
			for(int kx=0;kx<dx;++kx)
			{
				unsigned short sym=block[rowstride*ky+kx]-mean+half;
				if(sym>=nlevels)
					LOG_ERROR("Out of range  sym 0x%04X  at XY %d %d", sym, kx, ky);
				if(abs(sym-half)>256)
					printf("sym 0x%04X=%d\n", sym, sym);
				++tCDF[sym];
			}
		}
		int sum=0;
		for(int k=0;k<nlevels;++k)
		{
			unsigned freq=tCDF[k];
			tCDF[k]=sum;
			sum+=freq;
		}
		debug_print_CDF(nlevels, sum);
	}
#endif

	//lossless CDF
#if 1
	memset(tCDF, 0, nlevels*sizeof(unsigned));
	for(int ky=0;ky<dy;++ky)
	{
		for(int kx=0;kx<dx;++kx)
		{
			unsigned short sym=block[rowstride*ky+kx]+half;
			if(sym>=nlevels)
				LOG_ERROR("Out of range  sym 0x%04X  at XY %d %d", sym, kx, ky);
			++tCDF[sym];
		}
	}
	int start=0, end=nlevels;
	for(;!tCDF[start];++start);
	for(;!tCDF[end-1];--end);
	
	int nsym=end-start;
	size_t csize=list->nobj;
	dlist_push_back(list, "A", 1);//ANS marker
	dlist_push_back(list, &start, 2);
	dlist_push_back(list, &nsym, 2);
	for(int k=start;k<end;++k)
	{
		unsigned freq=tCDF[k];
		if(freq<0x80)
			dlist_push_back(list, &freq, 1);
		else
		{
			unsigned char *p=(unsigned char*)&freq, p2[]={p[0]|0x80, p[1]};
			dlist_push_back(list, p2, 2);
		}
	}

	int count=dx*dy, count_npot=(count-1)&count, sum;
	if(count_npot)//for boundary blocks
		sum=test4_encode4_prep(tCDF, nlevels);
	else
	{
		int lg_count=floor_log2(count);
		sum=0;
		for(int k=0;k<nlevels;++k)//works only for POT count
		{
			unsigned freq=tCDF[k];
			freq<<=16-lg_count;
			tCDF[k]=sum;
			sum+=freq;
		}
	}
	if(sum!=0x10000)
		LOG_ERROR("CDF sum 0x%04X", sum);

	size_t cdfsize=list->nobj-csize;
	if(cdfsize==7)
		cdfsize=7;
#endif
	//double precision
#if 0
	double sdev=0;
	for(int ky=0;ky<dy;++ky)
	{
		for(int kx=0;kx<dx;++kx)
		{
			short val=block[rowstride*ky+kx]-mean;
			sdev+=val*val;
		}
	}
	sdev/=dx*dy;
	sdev=sqrt(sdev);

	unsigned short ssdev=(unsigned short)(sdev*0x10000/half*2/5);
	dlist_push_back(list, &ssdev, 2);

	test4_encode4_prep(ssdev, tCDF, nlevels);
#endif

	//DEBUG
#if 0
	if(!oneshot&&dx==128&&dy==128)//DEBUG
	{
		//debug_print_CDF(nlevels, 0x10000);
		oneshot=1;
	}
#endif
	
	unsigned state=0x10000;
	for(int ky=dy-1;ky>=0;--ky)
	{
		for(int kx=dx-1;kx>=0;--kx)
		{
			unsigned short sym=block[rowstride*ky+kx]+half;
			if(sym>nlevels)
				LOG_ERROR("nlevels %d sym %d=0x%04X", nlevels, sym, sym);

			unsigned cdf=tCDF[sym], freq;
			{
				unsigned next;
				if(sym==vmax)
					next=0x10000;
				else
				{
					next=tCDF[sym+1];
					if(next<cdf)//CDF will overflow with trailing zero frequency symbols  (uint16)0x10000==0
						next=0x10000;
				}
				freq=(unsigned)(next-cdf);
#ifdef _DEBUG
				if(!freq)
					LOG_ERROR("rANS zero freq sym 0x%02X at XY %5d,%5d", sym, kx, ky);
#endif
			}

			if(state>=(unsigned)(freq<<16))//renorm
			{
				dlist_push_back(list, &state, 2);
				state>>=16;
			}

			state=state/freq<<16|(cdf+state%freq);//update
		}
	}
	dlist_push_back(list, &state, 4);
	csize=list->nobj-csize;

	int res=dx*dy;
	if((csize<<3)>(size_t)res*nbits)//bypass is better
	{
		printf("Bypassing %d*%d %d bit = %g bytes  was %lld + %lld = %lld bytes\n", dx, dy, nbits, res*nbits/8., cdfsize, csize-cdfsize, csize);

		dlist_pop_back(list, csize);
		dlist_push_back(list, "B", 1);//bypass marker
		unsigned char *dst=(unsigned char*)block;
		int ks=0, kd=0;
		//if(res&7)
		//	memset(block+res, 0, (BLOCKDIM*BLOCKDIM-res)*sizeof(short));
		//	LOG_ERROR("Invalid block res %d = %d*%d", res, dx, dy);
		switch(nbits)
		{
		case 8:
			for(ks=1;ks<res;++ks)
				dlist_push_back(list, block+ks, 1);
			break;
		case 9:
			for(ks=0;ks<res;ks+=8)
			{
				unsigned char p3[]=
				{
					(unsigned char)block[ks  ],
					(unsigned char)block[ks+1],
					(unsigned char)block[ks+2],
					(unsigned char)block[ks+3],
					(unsigned char)block[ks+4],
					(unsigned char)block[ks+5],
					(unsigned char)block[ks+6],
					(unsigned char)block[ks+7],

					(block[ks+7]>>8)<<7|
					(block[ks+6]>>8)<<6|
					(block[ks+5]>>8)<<5|
					(block[ks+4]>>8)<<4|
					(block[ks+3]>>8)<<3|
					(block[ks+2]>>8)<<2|
					(block[ks+1]>>8)<<1|
					(block[ks  ]>>8)
				};
				dlist_push_back(list, p3, 9);
			}
			break;
		case 10:
			for(ks=0;ks<res;ks+=8)
			{
				unsigned char p3[]=
				{
					(unsigned char)block[ks  ],
					(unsigned char)block[ks+1],
					(unsigned char)block[ks+2],
					(unsigned char)block[ks+3],
					(unsigned char)block[ks+4],
					(unsigned char)block[ks+5],
					(unsigned char)block[ks+6],
					(unsigned char)block[ks+7],

					(block[ks+7]>>8&1)<<7|
					(block[ks+6]>>8&1)<<6|
					(block[ks+5]>>8&1)<<5|
					(block[ks+4]>>8&1)<<4|
					(block[ks+3]>>8&1)<<3|
					(block[ks+2]>>8&1)<<2|
					(block[ks+1]>>8&1)<<1|
					(block[ks  ]>>8&1),

					(block[ks+7]>>9)<<7|
					(block[ks+6]>>9)<<6|
					(block[ks+5]>>9)<<5|
					(block[ks+4]>>9)<<4|
					(block[ks+3]>>9)<<3|
					(block[ks+2]>>9)<<2|
					(block[ks+1]>>9)<<1|
					(block[ks  ]>>9)
				};
				dlist_push_back(list, p3, 10);
			}
			break;
		case 11:
			for(ks=0;ks<res;ks+=8)
			{
				unsigned char p3[]=
				{
					(unsigned char)block[ks  ],
					(unsigned char)block[ks+1],
					(unsigned char)block[ks+2],
					(unsigned char)block[ks+3],
					(unsigned char)block[ks+4],
					(unsigned char)block[ks+5],
					(unsigned char)block[ks+6],
					(unsigned char)block[ks+7],

					(block[ks+7]>>8&1)<<7|
					(block[ks+6]>>8&1)<<6|
					(block[ks+5]>>8&1)<<5|
					(block[ks+4]>>8&1)<<4|
					(block[ks+3]>>8&1)<<3|
					(block[ks+2]>>8&1)<<2|
					(block[ks+1]>>8&1)<<1|
					(block[ks  ]>>8&1),

					(block[ks+7]>>9&1)<<7|
					(block[ks+6]>>9&1)<<6|
					(block[ks+5]>>9&1)<<5|
					(block[ks+4]>>9&1)<<4|
					(block[ks+3]>>9&1)<<3|
					(block[ks+2]>>9&1)<<2|
					(block[ks+1]>>9&1)<<1|
					(block[ks  ]>>9&1),

					(block[ks+7]>>10)<<7|
					(block[ks+6]>>10)<<6|
					(block[ks+5]>>10)<<5|
					(block[ks+4]>>10)<<4|
					(block[ks+3]>>10)<<3|
					(block[ks+2]>>10)<<2|
					(block[ks+1]>>10)<<1|
					(block[ks  ]>>10)
				};
				dlist_push_back(list, p3, 11);
			}
			break;
		}
		dlist_push_back(list, block, kd);
		return 1;
	}
	return 0;
}
static int test4_encode3_sq(short *block, int bw, int bh, int nbits, DList *list, int *ret_ntotal, int *ret_nbypass, short *temp)
{
	int ntotal=1, nbypass=0, vmax=(1<<nbits)-1, mindim=32;//3
	ArrayHandle sizes=dwt2d_gensizes(bw, bh, mindim, mindim, 0);
	squeeze_2d_fwd(block, (DWTSize*)sizes->data, (int)sizes->count, 0, temp);
	//squeeze_2d_fwd(block, (DWTSize*)sizes->data, (int)sizes->count, vmax, temp);

	//int minmax[2];
	//debug_get_minmax(block, bw, bh, minmax);
	//printf("block [%d, %d]\n", minmax[0], minmax[1]);

	//debug_save_block(block, bw, bh, 1<<(nbits+2-1));

	//test4_encode4_ans(block, bw, bh, bw, nbits+2, list);
	DWTSize *p=(DWTSize*)array_at(&sizes, sizes->count-1), *p2;
	nbypass+=test4_encode4_ans(block, p->w, p->h, bw, nbits, list);
	for(int k=(int)sizes->count-1;k>0;--k)
	{
		p=(DWTSize*)array_at(&sizes, k);
		p2=(DWTSize*)array_at(&sizes, k-1);
		nbypass+=test4_encode4_ans(block        +p->w, p2->w-p->w,       p->h, bw, nbits+1, list);
		nbypass+=test4_encode4_ans(block+bw*p->h     ,       p->w, p2->h-p->h, bw, nbits+1, list);
		nbypass+=test4_encode4_ans(block+bw*p->h+p->w, p2->w-p->w, p2->h-p->h, bw, nbits+2, list);
		ntotal+=3;
	}

	array_free(&sizes);
	if(ret_ntotal)
		*ret_ntotal+=ntotal;
	if(ret_nbypass)
		*ret_nbypass+=nbypass;
	return 1;
}
static int test4_encode2_xyb(const unsigned char *src, int bw, int bh, int x, int y, int dx, int dy, DList *list, int *ret_ntotal, int *ret_nbypass)
{
	int ntotal=0, nbypass=0;
	int kx, ky;
	
#ifdef ALLOC_BLOCK
	short *trow=(short*)malloc((dx>dy?dx:dy)*sizeof(short));
	short *tblock=(short*)malloc((size_t)dx*dy*sizeof(short));
	if(!tblock)
		return 0;
#endif

	for(ky=0;ky<dy;++ky)//X = R-G
	{
		for(kx=0;kx<dx;++kx)
		{
			int srcidx=(bw*(y+ky)+x+kx)<<2;
			tblock[dx*ky+kx]=src[srcidx]-src[srcidx|1]+255;
		}
	}
	test4_encode3_sq(tblock, dx, dy, 9, list, ret_ntotal, ret_nbypass, trow);

	for(ky=0;ky<dy;++ky)//Y = R+G
	{
		for(kx=0;kx<dx;++kx)
		{
			int srcidx=(bw*(y+ky)+x+kx)<<2;
			tblock[dx*ky+kx]=src[srcidx]+src[srcidx|1];
		}
	}
	test4_encode3_sq(tblock, dx, dy, 9, list, ret_ntotal, ret_nbypass, trow);

	for(ky=0;ky<dy;++ky)//B
	{
		for(kx=0;kx<dx;++kx)
		{
			int srcidx=(bw*(y+ky)+x+kx)<<2;
			tblock[dx*ky+kx]=src[srcidx|2];
		}
	}
	test4_encode3_sq(tblock, dx, dy, 8, list, ret_ntotal, ret_nbypass, trow);

#ifdef ALLOC_BLOCK
	free(trow);
	free(tblock);
#endif
	return 1;
}
size_t test4_encode(const unsigned char *src, int bw, int bh, ArrayHandle *data)
{
	DList list;
	size_t dststart;
	int ntotal=0, nbypass=0;
#ifdef ALLOC_BLOCK
	int blockw=bw, blockh=bh;
#else
	int blockw=BLOCKDIM, blockh=BLOCKDIM;
#endif
	if(*data)
	{
		if(data[0]->esize!=1)
			return 0;
		dststart=data[0]->count;
	}
	else
		dststart=0;
	dlist_init(&list, 1, 1024, 0);
	{
		int kx, xend, ky, yend;
		for(ky=0, yend=bh-(blockh-1);ky<yend;ky+=blockh)
		{
			for(kx=0, xend=bw-(blockw-1);kx<xend;kx+=blockw)
				test4_encode2_xyb(src, bw, bh, kx, ky, blockw, blockh, &list, &ntotal, &nbypass);
			if(kx<bw)
				test4_encode2_xyb(src, bw, bh, kx, ky, bw-kx, blockh, &list, &ntotal, &nbypass);
		}
		if(ky<bh)
		{
			for(kx=0, xend=bw-(blockw-1);kx<xend;kx+=blockw)
				test4_encode2_xyb(src, bw, bh, kx, ky, blockw, bh-ky, &list, &ntotal, &nbypass);
			if(kx<bw)
				test4_encode2_xyb(src, bw, bh, kx, ky, bw-kx, bh-ky, &list, &ntotal, &nbypass);
		}
	}

	printf("Bypass %d / %d blocks\n", nbypass, ntotal);

	size_t csize=list.nobj;
	dlist_appendtoarray(&list, data);
	dlist_clear(&list);
	return csize;
}
int test4_decode(const unsigned char *data, size_t srclen, int bw, int bh, unsigned char *buf)
{
	return 0;
}
#endif


//test4:  8 bit  3 channels  no alpha
static const int tag_test04='T'|'0'<<8|'0'<<16|'4'<<24;
typedef struct CDFInfoStruct
{
	unsigned sym, freq, qfreq;
} CDFInfo;
static int cdfinfo_byfreq(const void *left, const void *right)
{
	CDFInfo const *a, *b;

	a=(CDFInfo const*)left;
	b=(CDFInfo const*)right;
	return (a->freq>b->freq)-(a->freq<b->freq);
}
static int cdfinfo_bysym(const void *left, const void *right)
{
	CDFInfo const *a, *b;

	a=(CDFInfo const*)left;
	b=(CDFInfo const*)right;
	return (a->sym>b->sym)-(a->sym<b->sym);
}
static void balance_CDF(CDFInfo *info, int nlevels, unsigned short *CDF)
{
	unsigned prob_max=0x10000-(nlevels-1);
	int sum=0;
	for(int k=0;k<nlevels;++k)
	{
		CDFInfo *c=info+k;
		sum+=c->freq;
	}
	for(int k=0;k<nlevels;++k)
	{
		CDFInfo *c=info+k;
		c->qfreq=(unsigned)(((unsigned long long)c->freq<<16)/sum);
	}
	isort(info, nlevels, sizeof(CDFInfo), cdfinfo_byfreq);
	int idx=0;
	for(;idx<nlevels&&!info[idx].freq;++idx);
	for(;idx<nlevels&&!info[idx].qfreq;++idx)
		++info[idx].qfreq;
	for(idx=nlevels-1;idx>=0&&info[idx].qfreq>=prob_max;--idx);
	for(++idx;idx<nlevels;++idx)
		info[idx].qfreq=prob_max;

	int error=-0x10000;//too much -> +ve error & vice versa
	for(int k=0;k<nlevels;++k)
		error+=info[k].qfreq;
	if(error>0)
	{
		while(error)
		{
			for(int k=0;k<nlevels&&error;++k)
			{
				int dec=info[k].qfreq>1;
				info[k].qfreq-=dec, error-=dec;
			}
		}
	}
	else
	{
		while(error)
		{
			for(int k=nlevels-1;k>=0&&error;--k)
			{
				int inc=info[k].qfreq<prob_max;
				info[k].qfreq+=inc, error+=inc;
			}
		}
	}
	isort(info, nlevels, sizeof(CDFInfo), cdfinfo_bysym);
	sum=0;
	for(int k=0;k<nlevels;++k)
	{
		unsigned freq=info[k].qfreq;
		CDF[k]=(unsigned short)sum;
		sum+=freq;
	}
}
static void CDF_init(CDFInfo *CDF, int nlevels)
{
	for(int k=0;k<nlevels;++k)
	{
		CDFInfo *c=CDF+k;
		c->sym=k;
		c->freq=0;
		c->qfreq=0;
	}
}
static void CDF_fill2(const short *buf, int rowstride, int dx, int dy, int nbits, CDFInfo *CDF)
{
	int half=1<<(nbits-1);
	for(int ky=0;ky<dy;++ky)
	{
		for(int kx=0;kx<dx;++kx)
		{
			unsigned short val=buf[rowstride*ky+kx]+half;
			if(val>=(half<<1))
				LOG_ERROR("fill_CDF2()  val %d >= nlevels %d", val, half<<1);
			++CDF[val].freq;
		}
	}
}
static void CDF_fill(const short *buf, DWTSize *sizes, int nsizes, int basenbits, CDFInfo *CDF1, CDFInfo *CDF2)
{
	int bw=sizes->w, bh=sizes->h;
	for(int ks=nsizes-2;ks>=0;--ks)
	{
		DWTSize *p2=sizes+ks, *p=sizes+ks+1;
		CDF_fill2(buf        +p->w, bw, p2->w-p->w,       p->h, basenbits  , CDF1);
		CDF_fill2(buf+bw*p->h     , bw,       p->w, p2->h-p->h, basenbits  , CDF1);
		CDF_fill2(buf+bw*p->h+p->w, bw, p2->w-p->w, p2->h-p->h, basenbits+1, CDF2);
	}
}
static void test4_encode3_ans(const short *buf, int rowstride, int dx, int dy, int nbits, unsigned *state, unsigned short *CDF, DList *list)
{
	int half=1<<(nbits-1), vmax=(1<<nbits)-1;
	size_t csize=list->nobj;
	for(int ky=dy-1;ky>=0;--ky)
	{
		for(int kx=dx-1;kx>=0;--kx)
		{
			unsigned short sym=buf[rowstride*ky+kx]+half;
			int cdf, freq, next;
			if(CDF)
			{
				cdf=CDF[sym];
				if(sym==vmax)
					next=0x10000;
				else
				{
					next=CDF[sym+1];
					if(next<cdf)
						next=0x10000;
				}
				freq=next-cdf;
#ifdef _DEBUG
				if(!freq)
					LOG_ERROR("ANS zero freq sym 0x%02X at XY %5d,%5d", sym, kx, ky);
#endif
			}
			else
			{
				cdf=sym<<(16-nbits);
				freq=1<<(16-nbits);
			}

			if(*state>=(unsigned)(freq<<16))//renorm
			{
				dlist_push_back(list, state, 2);
				*state>>=16;
			}

			*state=*state/freq<<16|(cdf+*state%freq);//update
		}
	}
	//printf("%d*%d\t%lld\n", dx, dy, list->nobj-csize);
}
static void test4_encode2_ch(const short *buf, DWTSize *sizes, int nsizes, int basenbits, unsigned short *CDF1, unsigned short *CDF2, DList *list)
{
	int bw=sizes->w, bh=sizes->h;
	unsigned state=0x10000;
	dlist_push_back(list, &tag_test04, 4);
	test4_encode3_ans(buf, bw, sizes[nsizes-1].w, sizes[nsizes-1].h, basenbits-1, &state, 0, list);
	for(int ks=nsizes-2;ks>=0;--ks)
	{
		DWTSize *p2=sizes+ks, *p=sizes+ks+1;
		test4_encode3_ans(buf        +p->w, bw, p2->w-p->w,       p->h, basenbits  , &state, CDF1, list);
		test4_encode3_ans(buf+bw*p->h     , bw,       p->w, p2->h-p->h, basenbits  , &state, CDF1, list);
		test4_encode3_ans(buf+bw*p->h+p->w, bw, p2->w-p->w, p2->h-p->h, basenbits+1, &state, CDF2, list);
	}
	dlist_push_back(list, &state, 4);
}
size_t test4_encode(const unsigned char *src, int bw, int bh, ArrayHandle *data)
{
	DList list;
	int res=bw*bh, maxdim=bw>bh?bw:bh;
	short *bufX=(short*)malloc(((size_t)res*3+maxdim)*sizeof(short));
	CDFInfo *CDFinfo09=(CDFInfo*)malloc(3584*sizeof(CDFInfo));
	if(!bufX||!CDFinfo09)
		return 0;
	short *bufY=bufX+res, *bufB=bufY+res, *temprow=bufB+res;
	CDFInfo *CDFinfo10=CDFinfo09+512, *CDFinfo11=CDFinfo10+1024;

	for(int k=0;k<res;++k)//XYB from JPEG XL
	{
		int idx=k<<2;
		bufX[k]=src[idx]-src[idx|1];
		bufY[k]=src[idx]+src[idx|1];
		bufB[k]=src[idx|2];
	}

	int dstop=3;
	ArrayHandle sizes=dwt2d_gensizes(bw, bh, dstop, dstop, 0);
	DWTSize *psizes=(DWTSize*)sizes->data;
	int nsizes=(int)sizes->count;
	squeeze_2d_fwd(bufX, psizes, 0, nsizes, 1, 0, temprow);//squeeze transform from JPEG XL
	squeeze_2d_fwd(bufY, psizes, 0, nsizes, 1, 0, temprow);
	squeeze_2d_fwd(bufB, psizes, 0, nsizes, 1, 0, temprow);

	CDF_init(CDFinfo09,  512);
	CDF_init(CDFinfo10, 1024);
	CDF_init(CDFinfo11, 2048);
	CDF_fill(bufX, psizes, nsizes, 10, CDFinfo10, CDFinfo11);
	CDF_fill(bufY, psizes, nsizes, 10, CDFinfo10, CDFinfo11);
	CDF_fill(bufB, psizes, nsizes,  9, CDFinfo09, CDFinfo10);
	unsigned short *CDF09=(unsigned short*)malloc(3584*sizeof(short));
	if(!CDF09)
	{
		free(bufX);
		free(CDFinfo09);
		return 0;
	}
	unsigned short *CDF10=CDF09+512, *CDF11=CDF10+1024;
	balance_CDF(CDFinfo09,  512, CDF09);
	balance_CDF(CDFinfo10, 1024, CDF10);
	balance_CDF(CDFinfo11, 2048, CDF11);
	free(CDFinfo09);
	
	dlist_init(&list, 1, 1024, 0);
	//size_t csize=list.nobj;
	dlist_push_back(&list, CDF09, 3584*sizeof(short));
	//printf("CDF overhead: %lld\n", list.nobj-csize);
	//csize=list.nobj;
	test4_encode2_ch(bufX, psizes, nsizes, 10, CDF10, CDF11, &list);
	//printf("X: %lld\n", list.nobj-csize);
	//csize=list.nobj;
	test4_encode2_ch(bufY, psizes, nsizes, 10, CDF10, CDF11, &list);
	//printf("Y: %lld\n", list.nobj-csize);
	//csize=list.nobj;
	test4_encode2_ch(bufB, psizes, nsizes,  9, CDF09, CDF10, &list);
	//printf("B: %lld\n", list.nobj-csize);
	//csize=list.nobj;

	free(bufX);
	free(CDF09);
	array_free(&sizes);
	dlist_appendtoarray(&list, data);
	dlist_clear(&list);
	return 3584*sizeof(short);//return CDF overhead
}
int    test4_decode(const unsigned char *data, size_t srclen, int bw, int bh, unsigned char *buf)
{
	return 0;
}


#if 0
	#define T5_XYB_SQUEEZE
	#define T5_BLOCKDIM 256

//test5:  8 bit, 3 channels, no alpha (for prototype simplicity)		HEVC-like square blocks (useless for lossless)
static const int tag_test05='T'|'0'<<8|'0'<<16|'5'<<24;
typedef struct BlockInfoStruct
{
	unsigned short dim;//block is always square
	unsigned char lgdim;
	long long mean[3], sdev[3];
} BlockInfo;
BlockInfo t5_info[256];
int    test5_calcvar(const short *buf, int bw, int bh, int x, int y, BlockInfo *b)
{
	int xend=MINVAR(x+b->dim, bw), yend=MINVAR(y+b->dim, bh);
	b->sdev[0]=0;
	b->sdev[1]=0;
	b->sdev[2]=0;
	for(int ky=y;ky<yend;++ky)
	{
		for(int kx=x;kx<xend;++kx)
		{
			long long val;

			val=(buf[3*(bw*ky+kx)  ]<<16)-b->mean[0], b->sdev[0]+=val*val>>16;
			val=(buf[3*(bw*ky+kx)+1]<<16)-b->mean[1], b->sdev[1]+=val*val>>16;
			val=(buf[3*(bw*ky+kx)+2]<<16)-b->mean[2], b->sdev[2]+=val*val>>16;
		}
	}
	int count=b->w*b->h;
	b->sdev[0]/=count;
	b->sdev[1]/=count;
	b->sdev[2]/=count;
	//sqrt omitted
	return 0;
}
int    test5_calcmeansdev(const short *buf, int bw, int bh, int x, int y, BlockInfo *b)
{
	int xend=MINVAR(x+b->dim, bw), yend=MINVAR(y+b->dim, bh);
	int count=0;
	b->mean[0]=0;
	b->mean[1]=0;
	b->mean[2]=0;
	for(int ky=y;ky<yend;++ky)
	{
		for(int kx=x;kx<xend;++kx)
		{
			b->mean[0]+=buf[3*(bw*ky+kx)  ];
			b->mean[1]+=buf[3*(bw*ky+kx)+1];
			b->mean[2]+=buf[3*(bw*ky+kx)+2];
			++count;
		}
	}
	b->mean[0]=(b->mean[0]<<16)/count;
	b->mean[1]=(b->mean[1]<<16)/count;
	b->mean[2]=(b->mean[2]<<16)/count;

	test5_calcvar(buf, bw, bh, x, y, b);
	return count;
}
void   test5_emit_pbits(int x, int y, int dim, ArrayHandle *pbits)
{
	if(dim>4)
	{
		unsigned char temp=t5_info[y<<4|x].dim<dim;
		ARRAY_APPEND(*pbits, &temp, 1, 1, 0);
		if(temp)
		{
			dim>>=1;
			test5_emit_pbits(x    , y    , dim, pbits);
			test5_emit_pbits(x+dim, y    , dim, pbits);
			test5_emit_pbits(x    , y+dim, dim, pbits);
			test5_emit_pbits(x+dim, y+dim, dim, pbits);
		}
	}
}
int    test5_encode2_block(const short *buf, int bw, int bh, int x, int y, DList *list)
{
	ArrayHandle pbits;
	BlockInfo *p;
	memset(t5_info, 0, sizeof(t5_info));
	//for(int k=0;k<256;++k)
	//{
	//	BlockInfo *p=t5_info+k;
	//	p->lgw=2;
	//	p->lgh=2;
	//	p->w=4;
	//	p->h=4;
	//	memset(p->mean, 0, sizeof(long long[3]));
	//	memset(p->sdev, 0, sizeof(long long[3]));
	//}
	for(int ky=y, yend=MINVAR(y+T5_BLOCKDIM, bh);ky<yend;ky+=4)
	{
		for(int kx=0, xend=MINVAR(x+T5_BLOCKDIM, bw);kx<xend;kx+=4)
		{
			p=t5_info+(ky<<4|kx);
			p->dim=4;
			p->lgdim=2;
			test5_calcmeansdev(buf, bw, bh, kx, ky, 4, 4, p);
		}
	}
	int threshold_mean=8<<16, threshold_sdev=16*16<<16;
	for(int it=1;it<4;++it)
	{
		int step=1<<it, halfstep=step>>1;
		for(int ky=0;ky<16;ky+=step)
		{
			for(int kx=0;kx<16;kx+=step)//for each subblock
			{
				int merge=1;
				p=t5_info+(ky<<4|kx);
				for(int k=0;k<3;++k)
				{
					merge&=llabs(p->mean[k]-p[            halfstep].mean[k])<threshold_mean&&llabs(p->sdev[k]-p[            halfstep].sdev[k])<threshold_sdev;
					merge&=llabs(p->mean[k]-p[halfstep<<4         ].mean[k])<threshold_mean&&llabs(p->sdev[k]-p[halfstep<<4         ].sdev[k])<threshold_sdev;
					merge&=llabs(p->mean[k]-p[halfstep<<4|halfstep].mean[k])<threshold_mean&&llabs(p->sdev[k]-p[halfstep<<4|halfstep].sdev[k])<threshold_sdev;
				}
				if(merge)
				{
					p->dim<<=1;
					p->lgdim+=1;
					for(int k=0;k<3;++k)
						p->mean[k]=(p->mean[k]+p[halfstep].mean[k]+p[halfstep<<4].mean[k]+p[halfstep<<4|halfstep].mean[k])>>2;
					test5_calcvar(buf, bw, bh, x+(kx<<2), y+(ky<<2), p);

					int dim=p->dim>>2;
					memset(t5_info+(ky<<4|(kx+1)), 0, (dim-1)*sizeof(BlockInfo));
					for(int ky2=1;ky2<dim;++ky2)
						memset(t5_info+((ky+ky2)<<4|kx), 0, dim*sizeof(BlockInfo));
				}
#if 0
				int idx=ky<<4|kx;
				BlockInfo *p=t5_info+idx, *p2;
				if(p->lgw&&p->lgh)
				{
					int w1=(1<<(p->lgw-2));
					if(kx+w1<16)
					{
						p2=t5_info+idx+w1;
						if(p2->lgw==p->lgw&&p2->lgh==p->lgh)
						{
							int merge=1;
							for(int k=0;k<3;++k)
								merge&=llabs(p->mean[k]-p2->mean[k])<threshold_mean&&llabs(p->sdev[k]-p2->sdev[k])<threshold_sdev;
							if(merge)//merge blocks horizontally
							{
								int lgw=p->lgw;
								t5_info[idx  ].lgw=3;
								t5_info[idx|1].lgw=0;
							}
						}
					}
				}
#endif
			}
		}
	}
	ARRAY_ALLOC(char, pbits, 0, 0, 0, 0);
	test5_emit_pbits(0, 0, 256, &pbits);
	//unsigned char temp=t5_info->dim<256;//dim < max
	//ARRAY_APPEND(pbits, &temp, 1, 1, 0);
	//if(temp)
	//{
	//	temp=t5_info->dim<128;
	//	ARRAY_APPEND(pbits, &temp, 1, 1, 0);
	//
	//	temp=t5_info[8].dim<128;
	//	ARRAY_APPEND(pbits, &temp, 1, 1, 0);
	//
	//	temp=t5_info[8<<4].dim<128;
	//	ARRAY_APPEND(pbits, &temp, 1, 1, 0);
	//
	//	temp=t5_info[8<<4|8].dim<128;
	//	ARRAY_APPEND(pbits, &temp, 1, 1, 0);
	//}

	//for(int ky=0;ky<16;++ky)
	//{
	//	for(int kx=0;kx<16;)//for each subblock
	//	{
	//		p=t5_info+(ky<<4|kx);
	//
	//		kx+=p->dim;
	//	}
	//}
}
size_t test5_encode(const unsigned char *src, int bw, int bh, ArrayHandle *data)
{
	size_t res=(size_t)bw*bh, srclen=res<<2;
	short *buf=(short*)malloc(res*3*sizeof(short));
	if(!buf)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	int maxdim=MAXVAR(bw, bh);
	short *temprow=(short*)malloc(maxdim*sizeof(short));
	if(!temprow)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

#ifdef T5_XYB_SQUEEZE
	int nbits=11;
	for(int ks=0, kd=0;ks<srclen-3;ks+=4, kd+=3)//XYB from JPEG XL
	{
		buf[kd]=src[ks]-src[ks|1];
		buf[kd+1]=src[ks]+src[ks|1];
		buf[kd+2]=src[ks|2];
	}
	{
		int dstop=3;
		ArrayHandle sizes=dwt2d_gensizes(bw, bh, dstop, dstop, 0);
		DWTSize *psizes=(DWTSize*)sizes->data;
		int nsizes=(int)sizes->count;
		squeeze_2d_fwd(buf  , psizes, nsizes, 3, 0, temprow);//squeeze transform from JPEG XL
		squeeze_2d_fwd(buf+1, psizes, nsizes, 3, 0, temprow);
		squeeze_2d_fwd(buf+2, psizes, nsizes, 3, 0, temprow);
		array_free(&sizes);
	}
#else
	int nbits=8;
	for(int ks=0, kd=0;ks<srclen-3;ks+=4, kd+=3)
	{
		buf[kd]=src[ks];
		buf[kd+1]=src[ks|1];
		buf[kd+2]=src[ks|2];
	}
#endif
	for(int ky=0;ky<bh-255;ky+=T5_BLOCKDIM)
	{
		for(int kx=0;kx<bw-255;kx+=T5_BLOCKDIM)//for each full block		TODO handle boundary
			test5_encode2_block(buf, bw, bh, kx, ky);
	}
}
#endif


	#define T5_ANS_64BIT

#ifdef T5_ANS_64BIT
typedef unsigned long long T5ANSState;
#else
typedef unsigned T5ANSState;//doesn't work
#endif
//test5 v5:  8 bit, 3 channels, no alpha			subband prediction and adaptive binary ANS
static const int tag_test05='T'|'0'<<8|'0'<<16|'5'<<24;
typedef short (*T5_CALCMEAN)(const short *buf, int bw, int bh, int x, int y, int c);
static short t5_calcmean_diffx(const short *buf, int bw, int bh, int x, int y, int c)
{
	int idx=3*(bw*y+x)+c;

	short val;
	if(x>0&&x+1<bw)
		val=(buf[idx-3]-buf[idx+3])>>2;
	else		//TODO better boundary handling
		val=buf[idx];
	
	//short val=buf[idx];
	//if(x+1<bw)
	//	val-=buf[idx+3];

	return val;
}
static short t5_calcmean_diffy(const short *buf, int bw, int bh, int x, int y, int c)
{
	int idx=3*(bw*y+x)+c;

	short val;
	if(y>0&&y+1<bh)
		val=(buf[idx-3*bw]-buf[idx+3*bw])>>2;
	else
		val=buf[idx];
	//short val=buf[idx];
	//if(y+1<bh)
	//	val-=buf[idx+3*bw];

	return val;
}
static short t5_calcmean_diffxy(const short *buf, int bw, int bh, int x, int y, int c)
{
	int idx=3*(bw*y+x)+c;
	
	short val;
	if(x>0&&x+1<bw&&y>0&&y+1<bh)
		val=-(buf[idx-3*bw-3]-buf[idx-3*bw+3]-buf[idx+3*bw-3]+buf[idx+3*bw+3])>>2;
	else
		val=buf[idx];
	//if(x>0)
	//{
	//	if(y>0)
	//		val=-(buf[idx]-buf[idx-3]-buf[idx-3*bw]+buf[idx-3*bw-3]);
	//	else
	//		val=-(buf[idx]-buf[idx-3]);
	//}
	//else
	//{
	//	if(y>0)
	//		val=-(buf[idx]-buf[idx-3*bw]);
	//	else
	//		val=buf[idx];
	//}
	//short val=buf[idx];
	//if(x+1<bw)
	//{
	//	val-=buf[idx+3];
	//	if(y+1<bh)
	//		val-=buf[idx+3*bw]-buf[idx+3*bw+3];
	//}
	//else
	//{
	//	if(y+1<bh)
	//		val-=buf[idx+3*bw];
	//}

	return val;
}
//static void t5_predict(short *buf, int bw, int bh, int srcx, int srcy, int srcdx, int srcdy, int dstx, int dsty, int dstdx, int dstdy, CALCMEANCONF calcmeanconf, unsigned short *CDF, unsigned short *freq)
//static void t5_predict(const short *buf, int bw, int quadrant, int dx, int dy, unsigned short *CDF, unsigned short *freq)
//{
//}
long long error_func_p32(long long x)
{
	//approximation 3 from Wikipedia
	const unsigned c[]=
	{
		0x120DCCEB,//0.0705230784
		0x0AD2FE74,//0.0422820123
		0x025F8DA3,//0.0092705272
		0x0009F660,//0.0001520143
		0x00122007,//0.0002765672
		0x0002D27E,//0.0000430638
	};
	int neg=x<0;
	unsigned long long x0=(unsigned long long)llabs(x), lo;//32.32 bits
	
	if(x0>0x48E7A0CE1)//erf(4.5565498399) = 0x0.FFFFFFFF800003 ~= 1 in 32.32bit
		lo=0x100000000;
	else
	{
		long long hi;
	#define MUL_P32(DST, A, B)	(lo=_mul128(A, B, &hi), DST=hi<<32|lo>>32)
		MUL_P32(lo, x0, c[5]), lo+=c[4];
		MUL_P32(lo, lo, x0), lo+=c[3];
		MUL_P32(lo, lo, x0), lo+=c[2];
		MUL_P32(lo, lo, x0), lo+=c[1];
		MUL_P32(lo, lo, x0), lo+=c[0];
		MUL_P32(lo, lo, x0), lo+=0x100000000;

		lo=_div128(1, lo>>1, lo, &hi);//round((0x100000000<<32)/lo) = (1<<64|lo>>1)/lo
		MUL_P32(lo, lo, lo);
		MUL_P32(lo, lo, lo);
		MUL_P32(lo, lo, lo);
		MUL_P32(lo, lo, lo);

		lo=0x100000000-lo;
	}
	lo^=-neg;
	lo+=neg;
	lo>>=16;
	return (int)lo;
}
int error_func_p16(int x)
{
	//approximation 1
#if 0
	const int coeff[]=//1 - (1+(c[0]+(c[1]+(c[2]+c[3]*x)*x)*x)*x)^-4, max error=0.0005 = 0x20*2^-16
	{
		0x4745,//0.278393
		0x3AFB,//0.230389
		0x0040,//0.000972
		0x13FF,//0.078108
	};
	int x2=(int)((long long)coeff[3]*x>>16)+coeff[2];
	x2=(int)((long long)x2*x>>16)+coeff[1];
	x2=(int)((long long)x2*x>>16)+coeff[0];
	x2=(int)((long long)x2*x>>16)+0x10000;
	x2=0x10000/x2;
	x2=(int)((long long)x2*x2>>16);
	x2=(int)((long long)x2*x2>>16);
	x2=0x10000-x2;
	return x2;
#endif

	//approximation 3
#if 1
	const unsigned c[]=
	{
		0x120E,//0x120DCCEB,//0.0705230784
		0x0AD3,//0x0AD2FE74,//0.0422820123
		0x0260,//0x025F8DA3,//0.0092705272
		0x000A,//0x0009F660,//0.0001520143
		0x0012,//0x00122007,//0.0002765672
		0x0003,//0x0002D27E,//0.0000430638
	};
	int neg=x<0;
	unsigned long long x0=(unsigned long long)abs(x), res;//16.16 bits

	if(x0>0x32A1F)//erf(3.16453508) = 0x0.FFFF8000003 ~= 1 in 16.16 bit
		res=0x10000;
	else
	{
		res=(x0*c[5]>>16)+c[4];
		res=(res*x0>>16)+c[3];
		res=(res*x0>>16)+c[2];
		res=(res*x0>>16)+c[1];
		res=(res*x0>>16)+c[0];
		res=(res*x0>>16)+0x10000;

		res=(0x100000000|res>>1)/res;
		res=res*res>>16;
		res=res*res>>16;
		res=res*res>>16;
		res=res*res>>16;

		res=0x10000-res;
	}
	res^=-neg;
	res+=neg;
	return (int)res;
#endif
}
static int calc_phi(int sym, int mean, int nbits, int conf)
{
	//const int invsqrt2=0xB505;
	int x;

	x=sym-mean;
	x<<=16-nbits;//normalize to fixed 16 bit [-0x0.FFFF, 0x0.FFFF] exclusive
	//x-=1<<(16-nbits-2);//subtract half

	//++conf;//conf from [0x0.0001, 0x1.0000]
	//conf=(int)((long long)conf*invsqrt2>>16);//divide by sqrt2
	//conf+=!conf;

	//x*=conf;//multiply by conf without shift
	x=(int)((long long)x*conf>>16);//x is in [-1, 1] * conf is 16.16, tuned (pre-divided by sqrt2)
	//x<<=nbits;
	x=error_func_p16(x);
	x=x<<nbits|sym;//guard against zero change which causes division by zero CRASH
	return x;
}
static int cmin=0x10000, cmax=0, ccount=0;
static long long csum=0;
unsigned g_conf=0;
size_t ncorrect=0, ntotal=0;//DEBUG
//ArrayHandle LOL_debug=0;//
static void t5_enc2(int bw, int x, int y, int c, short nbits, const short *buf_mean, const unsigned short *buf_conf, const short *buf, T5ANSState *state, DList *list)
{
	static int ncalls=0;
	++ncalls;
	int idx=3*(bw*y+x)+c;
	int mean=buf_mean[idx]+(1<<(nbits-1)),
		conf=g_conf,
		//conf=buf_conf[idx],
		sym=buf[idx]+(1<<(nbits-1));
	
	//if(LOL_debug->data[idx]!=1)//
	//	LOG_ERROR("Overlap in enc2");//
	//LOL_debug->data[idx]=2;//
	
	if(sym<0||sym>=(1<<nbits))
		LOG_ERROR("Range error sym 0x%04X", sym);

	//conf=0xFFFF;
//	conf=(int)((0x193200193>>(16-nbits))/(abs(mean-(1<<(nbits-1)))+1));//(16bit+nbits)/[1bit, nbits] = conf[nbits.16bit, 0.16bit] * x[-0.16, 0.16] -> param[-nbits.16bit, nbits.16bit]
	conf=0x516666;//not adaptive
//	conf=150<<16;//not adaptive
//	conf=(0x64C63F>>(16-nbits))/(abs(mean-(1<<(nbits-1)))+1);
//	conf=(0x5AA7F0>>(16-nbits))/(abs(mean-(1<<(nbits-1)))+1);
//	conf=(0x5AA993>>(16-nbits))/(abs(mean-(1<<(nbits-1)))+1);
//	conf=(0x64C63F>>(16-nbits))/(abs(mean-(1<<(nbits-1)))+1);
//	conf=(0x10000<<nbits)/(abs(mean-(1<<(nbits-1)))+1);
	//conf=0x10000/(abs(sym-mean)+1);//CHEATING
	//conf=0x200000/((abs(mean-(1<<(nbits-1)))+1)<<(16-nbits))-1;
	//conf=(0x10000<<nbits)/(abs(mean-(1<<(nbits-1)))+1)-1;
	//conf*=conf;
	//if(conf>0xFFFF)
	//	conf=0xFFFF;
	//if(conf<0)
	//	conf=0;

#if 1
	if(cmin>conf)
		cmin=conf;
	if(cmax<conf)
		cmax=conf;
	csum+=conf;
	++ccount;
#endif

	//adaptive ANS
#if 1
	//mean=1<<(nbits-1);//half
	//mean=sym;//CHEAT
	int erf_start=calc_phi(0       , mean, nbits, conf),//always -65536
		erf_end  =calc_phi(1<<nbits, mean, nbits, conf),//always 65536 (proof required)
		erf_curr =calc_phi(sym     , mean, nbits, conf),
		erf_next =calc_phi(sym+1   , mean, nbits, conf);
	
	//if(erf_start!=-65536||erf_end!=65536)
	//	LOG_ERROR("CDF range %d ~ %d", erf_start, erf_end);

#ifdef T5_ANS_64BIT
	long long den=(long long)erf_end-erf_start, CDF, freq;
	//den=den<<nbits|((1LL<<nbits)-1);
	CDF =((long long)(erf_curr-erf_start)<<32|den>>1)/den;
	freq=((long long)(erf_next-erf_curr )<<32|den>>1)/den;

	//CDF|=(long long)sym<<(16-nbits);
	//freq|=1LL<<(16-nbits);

	static size_t sym_count=0;
	++sym_count;
	unsigned f_bypass=(unsigned)(1LL<<(32-nbits));
	if(freq<f_bypass)
		f_bypass=(unsigned)(1LL<<(32-nbits));


	if(CDF+freq>0x100000000)
		LOG_ERROR("CDF is not normalized CDF 0x%08llX freq 0x%08llX", CDF, freq);
	if(!freq)
		LOG_ERROR("Zero freq");

	//CDF=0, freq=0xFFFFFFFF;//CHEAT

	if(*state>=((T5ANSState)freq<<32))//renorm
	{
		dlist_push_back(list, state, 4);
		*state>>=32;
	}

	*state=*state/freq<<32|(CDF+*state%freq);//update
#else
	int den=erf_end-erf_start,
		CDF =(int)(((long long)(erf_curr-erf_start)<<(16+nbits)|den>>1)/den),
		freq=(int)(((long long)(erf_next-erf_curr )<<(16+nbits)|den>>1)/den);
	CDF =(CDF |1<<(nbits-1))>>nbits;//doesn't work
	freq=(freq|1<<(nbits-1))>>nbits;
	CDF =CDF &((1<<nbits)-1)|sym;
	freq=freq&((1<<nbits)-1)|1;

	if(*state>=((T5ANSState)freq<<16))//renorm
	{
		dlist_push_back(list, state, 2);
		*state>>=16;
	}

	*state=*state/freq<<16|(CDF+*state%freq);//update
#endif
#endif
	
	//adaptive binary ANS
#if 0
	static unsigned short p0[11]={0};
	int revealed=0,
		erf_start=calc_phi(0       , mean, nbits, conf),
		erf_end  =calc_phi(1<<nbits, mean, nbits, conf);
	for(int kb=nbits-1;kb>=0;--kb)//predict bits in order of decode (MSB -> LSB)
	{
		int truebit=1<<kb,
			erf_mid=calc_phi(revealed|truebit, mean, nbits, conf);
		//if(abs(erf_curr_half-erf_revealed)>0xFFFF)
		//	LOG_ERROR("Overflow in p0 expression");
		int den=erf_end-erf_start;
		p0[kb]=(unsigned short)(((long long)(erf_mid-erf_start)<<16|den>>1)/den);

		unsigned correction=(unsigned short)((1<<(15-kb))^-(mean>>kb&1));
		//unsigned correction=1<<(nbits-1-kb);
		//correction^=-(mean>>kb&1);
		//if(mean>>kb&1)
		//	correction=0xFFFF-correction;
		//correction=correction+((0xFFFF-correction)&-(mean>>kb&1));
		int temp=p0[kb]-0x8000;
		temp=(int)((long long)temp*correction>>16);
		p0[kb]=(unsigned short)(temp+0x8000);

		if(p0[kb]<1)//clamp p0
			p0[kb]=1;
		if(p0[kb]>=0xFFFF)
			p0[kb]=0xFFFF;
		//p0[kb]=mean>>kb&1?0xFFFF:0x0001;

		int bit=sym>>kb&1;

		int correct=bit!=(p0[kb]<0x8000);//p0 can't be changed at this point but conf can be updated for next time
		conf=correct<<15|conf>>1;

		if(ncalls==1)//
		{
			double den2=(double)(1LL<<(16+nbits));
			printf("%2d %lf\t%lf\t%lf\t0x%04X\tmean %d\tsym %d\tcorrect %d\n", kb, erf_start/den2, erf_mid/den2, erf_end/den2, p0[kb], mean>>kb&1, bit, correct);
		}

		ncorrect+=correct;//DEBUG
		++ntotal;
		
		if(bit)
		{
			revealed|=truebit;
			erf_start=erf_mid;
		}
		else
			erf_end=erf_mid;
	}
	for(int kb=0;kb<nbits;++kb)//encode bits in reverse (MSB <- LSB)
	{
		//p0[kb]=0x8000;//bypass all DEBUG

		int bit=sym>>kb&1;
		int CDF=p0[kb]&-bit, freq=(unsigned short)((p0[kb]^-bit)+bit);//fetch
		
		if(*state>=(unsigned)(freq<<16))//renorm
		{
			dlist_push_back(list, state, 2);
			*state>>=16;
		}
		*state=*state/freq<<16|(CDF+*state%freq);//update
	}
#endif
}
static void t5_update_conf(unsigned short *conf, short error, int nbits)
{
	unsigned c2=*conf+1;
	error=abs(error)+1;
	//c2=(c2<<(nbits>>2))/error;
	c2=(c2<<(nbits>>1))/error;
	//c2=(c2<<nbits)/error;
	if(c2<1)
		c2=1;
	if(c2>0x10000)
		c2=0x10000;
	*conf=(unsigned short)(c2-1);
}
static void t5_predict(int bw, int x, int y, int c, short nbits, short mean, unsigned short *conf, const short *buf, short *buf_mean, unsigned short *buf_conf)
{
	int idx=3*(bw*y+x)+c;
	buf_mean[idx]=mean;
	buf_conf[idx]=*conf;
	t5_update_conf(conf, buf[idx]-mean, nbits);
	
	//if(LOL_debug->data[idx])//
	//	LOG_ERROR("Overlap in predict");//
	//LOL_debug->data[idx]=1;//
}
#if 0
static void t5_enc2(int bw, int x, int y, int c, short nbits, short mean, unsigned short *conf, const short *buf, unsigned *state, DList *list)
{
	int idx=3*(bw*x+y)+c;
	//buf_mean[idx]=mean;
	//buf_conf[idx]=*conf;

	//binary adaptive ANS
	unsigned short p0[11];
	int sym=buf[idx]+(1<<(nbits-1));
	if(sym<0||sym>=(1<<nbits))
		LOG_ERROR("Range error sym 0x%04X", sym);
	int revealed=0,
		erf_revealed=calc_erf(revealed, mean, nbits, *conf);
	for(int kb=nbits-1;kb>=0;--kb)//predict bits in reverse
	{
		int erf_curr_full=calc_erf(revealed|1<<(kb+1), mean, nbits, *conf),
			erf_curr_half=calc_erf(revealed|1<<kb, mean, nbits, *conf),
			p0=((erf_curr_half-erf_revealed)<<16)/(erf_curr_full-erf_revealed);

		int bit=sym>>kb&1;
		int CDF=p0&-bit, freq=p0^-bit;//fetch
		
		if(*state>=(unsigned)(freq<<16))//renorm
		{
			dlist_push_back(list, state, 2);
			*state>>=16;
		}
		*state=*state/freq<<16|(CDF+*state%freq);//update

		if(bit)
		{
			revealed|=1<<kb;
			erf_revealed=erf_curr_full;
		}
	}

	t5_update_conf(conf, buf[idx]-mean);
}
#endif
size_t test5_encode(const unsigned char *src, int bw, int bh, ArrayHandle *data)
{
	size_t res=(size_t)bw*bh, srclen=res<<2;
	int maxdim=MAXVAR(bw, bh);
	short
		*buf=(short*)malloc(res*3*sizeof(short)),
		*temprow=(short*)malloc(maxdim*sizeof(short)),
		*buf_mean=(short*)malloc(res*3*sizeof(short));
	unsigned short *buf_conf=(unsigned short*)malloc(res*3*sizeof(short));
	DList list;

	if(!buf||!temprow||!buf_mean)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	for(int ks=0, kd=0;ks<srclen-3;ks+=4, kd+=3)//XYB from JPEG XL
	{
		//buf[kd  ]=src[ks]-src[ks|1];
		//buf[kd+1]=src[ks]+src[ks|1]-256;
		//buf[kd+2]=src[ks|2]-128;

		buf[kd  ]=src[ks  ]-128;//
		buf[kd+1]=src[ks|1]-128;//
		buf[kd+2]=src[ks|2]-128;//
	}

	int dstop=3;
	ArrayHandle sizes=dwt2d_gensizes(bw, bh, dstop, dstop, 0);
	DWTSize *psizes=(DWTSize*)sizes->data;
	int nsizes=(int)sizes->count;

	const int m0=0, c0=0x4000;//default values for mean & conf
	unsigned short conf;

	ncorrect=0;
	//ARRAY_ALLOC(char, LOL_debug, 0, 3*bw*bh, 0, 0);//
	
	for(int ks=0;ks<nsizes-1;++ks)
	//for(int ks=0;ks<1;++ks)
	{
		int dx=psizes[ks].w, dy=psizes[ks].h, px=dx>>1, py=dy>>1;
		for(int kc=0;kc<3;++kc)
		{
			squeeze_2d_fwd(buf+kc, psizes, ks, ks+2, 3, 0, temprow);//squeeze transform from JPEG XL

#if 0
			if(kc==2)
			{
				snprintf(g_buf, G_BUF_SIZE, "squeeze%02d.PNG", ks);
				save_16bit(g_buf, buf, bw, bh, 3, 1024, 16-11, 1);
			}
#endif

			for(int kq=1;kq<4;++kq)
			{
				T5_CALCMEAN calcmean=0;
				int ox, oy;
				short nbits;
				switch(kq)
				{
				case 1:calcmean=t5_calcmean_diffx , nbits= 9+(kc<2), ox=px, oy= 0;break;
				case 2:calcmean=t5_calcmean_diffy , nbits= 9+(kc<2), ox= 0, oy=py;break;
				case 3:calcmean=t5_calcmean_diffxy, nbits=10+(kc<2), ox=px, oy=py;break;
				default:LOG_ERROR("");return 0;
				}
				conf=c0;

				for(int ky=0;ky<py;++ky)//predict
				{
					for(int kx=0;kx<px;++kx)
					{
						short mean=calcmean(buf, bw, bh, kx, ky, kc);
						t5_predict(bw, ox+kx, oy+ky, kc, nbits, mean, &conf, buf, buf_mean, buf_conf);
					}
					if(dx&1)
						t5_predict(bw, dx-1, oy+ky, kc, nbits, m0, &conf, buf, buf_mean, buf_conf);
				}
				if(dy&1)
				{
					for(int kx=0;kx<px;++kx)
						t5_predict(bw, ox+kx, dy-1, kc, nbits, m0, &conf, buf, buf_mean, buf_conf);
					if(dx&1)
						t5_predict(bw, dx-1, dy-1, kc, nbits, m0, &conf, buf, buf_mean, buf_conf);
				}
			}
		}
	}
	conf=c0;
	for(int ky=psizes[nsizes-1].h-1;ky>=0;--ky)//bypass top-left subband
	{
		for(int kx=psizes[nsizes-1].w-1;kx>=0;--kx)
		{
			int idx=3*(bw*ky+kx);
			for(int kc=0;kc<3;++kc)
				t5_predict(bw, kx, ky, kc, 8+(kc<2), m0, &conf, buf, buf_mean, buf_conf);
		}
	}
	
	//calculate sdev
#if 0
	size_t len=3*res;
	long long sum=0;
	for(int k=0;k<(int)len;k+=3)
		sum+=buf[k]*buf[k];
	double sdev=sqrt((double)sum/len);
	printf("sdev %lf\n", sdev);
#endif

#ifdef T5_ANS_64BIT
	T5ANSState state=0x100000000;
#else
	T5ANSState state=0x10000;
#endif
	dlist_init(&list, 1, 1024, 0);
	for(int ks=0;ks<nsizes-1;++ks)
	//for(int ks=0;ks<1;++ks)
	{
		int dx=psizes[ks].w, dy=psizes[ks].h, px=dx>>1, py=dy>>1;
		for(int kc=2;kc>=0;--kc)
		{
			for(int kq=3;kq>=1;--kq)
			{
				int ox, oy;
				short nbits;
				switch(kq)
				{
				case 1:nbits= 9+(kc<2), ox=px, oy= 0;break;
				case 2:nbits= 9+(kc<2), ox= 0, oy=py;break;
				case 3:nbits=10+(kc<2), ox=px, oy=py;break;
				default:LOG_ERROR("");return 0;
				}

				if(dy&1)
				{
					if(dx&1)
						t5_enc2(bw, dx-1, dy-1, kc, nbits, buf_mean, buf_conf, buf, &state, &list);
					for(int kx=px-1;kx>=0;--kx)
						t5_enc2(bw, ox+kx, dy-1, kc, nbits, buf_mean, buf_conf, buf, &state, &list);
				}
				for(int ky=py-1;ky>=0;--ky)
				{
					if(dx&1)
						t5_enc2(bw, dx-1, oy+ky, kc, nbits, buf_mean, buf_conf, buf, &state, &list);
					for(int kx=px-1;kx>=0;--kx)
						t5_enc2(bw, ox+kx, oy+ky, kc, nbits, buf_mean, buf_conf, buf, &state, &list);
				}
			}
		}
	}
	for(int ky=psizes[nsizes-1].h-1;ky>=0;--ky)//top-left subband
	{
		for(int kx=psizes[nsizes-1].w-1;kx>=0;--kx)
		{
			int idx=3*(bw*ky+kx);
			for(int kc=0;kc<3;++kc)
				t5_enc2(bw, kx, ky, kc, 8+(kc<2), buf_mean, buf_conf, buf, &state, &list);
		}
	}
#ifdef T5_ANS_64BIT
	dlist_push_back(&list, &state, 8);
#else
	dlist_push_back(&list, &state, 4);
#endif
	
	printf("conf [0x%04X, 0x%04X] av 0x%04X %lf\n", cmin, cmax, (int)(csum/ccount), (double)csum/ccount);//
	//printf("size %lld / %d CR %lf\n", list.nobj, bw*bh*3, (double)bw*bh*3/list.nobj);//
	//printf("Predicted %lld / %lld = %.2lf%% bits  size %lld / %d CR %lf\n", ncorrect, ntotal, 100.*ncorrect/ntotal, list.nobj, bw*bh*3, (double)bw*bh*3/list.nobj);//
	//for(int k=0, end=3*bw*bh;k<end;++k)//
	//{
	//	if(LOL_debug->data[k]!=2)
	//		LOG_ERROR("Missing pixels");
	//}
	//array_free(&LOL_debug);

	dlist_appendtoarray(&list, data);
	dlist_clear(&list);

	//DEBUG save mean & conf
#if 0
	for(int k=0, end=3*bw*bh;k<end;++k)
	{
		short val=buf[k]+128;
		if(val<0)
			val=0;
		if(val>255)
			val=255;
		buf[k]=val;
	}
	save_16bit("buf.PNG", buf, 0, bw, bh, 3, 0, 8, 1);
#endif
#if 0
	save_16bit("buf.PNG", buf, 0, bw, bh, 3, 128, 16-8, 1);
	save_16bit("mean.PNG", buf_mean, 0, bw, bh, 3, 128, 16-8, 1);
	save_16bit("diff.PNG", buf, buf_mean, bw, bh, 3, 128, 16-8, 1);
	for(int k=0, len=res*3;k<len;k+=3)
	{
		short orig=buf[k], mean=buf_mean[k], diff=buf[k]-buf_mean[k];
		buf[k  ]=diff;
		buf[k+1]=mean;
		buf[k+2]=orig;
	}
	save_16bit("diff-mean-orig.PNG", buf, 0, bw, bh, 3, 128, 16-8, 1);

	//save_16bit("mean.PNG", buf_mean, 0, bw, bh, 3, 1024, 16-11, 1);
	//save_16bit("conf.PNG", buf_conf, 0, bw, bh, 3, 0, 0, 1);
	//for(int k=0, end=3*bw*bh;k<end;++k)//little->big endian
	//{
	//	buf_mean[k]+=1024;
	//	buf_mean[k]<<=16-11;
	//	buf_mean[k]=_byteswap_ushort(buf_mean[k]);
	//	buf_conf[k]=_byteswap_ushort(buf_conf[k]);
	//}
	//lodepng_encode_file("mean.PNG", (unsigned char*)buf_mean, bw, bh, LCT_RGB, 16);
	//lodepng_encode_file("conf.PNG", (unsigned char*)buf_conf, bw, bh, LCT_RGB, 16);
#endif

	array_free(&sizes);
	free(temprow);
	free(buf);
	free(buf_mean);
	free(buf_conf);
	return 1;
}
int    test5_decode(const unsigned char *data, size_t srclen, int bw, int bh, unsigned char *buf)
{
	return 0;
}


int g_offset=0x12000;
static int t6_calc_phi(int sym, int mean, int sdev)
{
	const int sqrt2=0x16A0A;
	int x;

	x=sym-mean;
	x<<=16;

	//sdev<<=16;
	sdev=(int)((long long)sdev*sqrt2>>16);
	sdev+=!sdev;

	x=((long long)x<<16)/sdev;

	x=error_func_p16(x);
	x=x<<8|sym;
	return x;
}
size_t test6_encode(const unsigned char *src, int bw, int bh, ArrayHandle *data)
{
	size_t res=(size_t)bw*bh, srclen=res<<2;
	unsigned char *buf=(unsigned char*)malloc(srclen);
	if(!buf)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	for(int k=0;k<srclen;k+=4)
	{
		buf[k  ]=src[k  ]-src[k|1]+128;//XGZ
		buf[k|1]=src[k|1];
		buf[k|2]=src[k|2]-src[k|1]+128;
		buf[k|3]=0xFF;

		//buf[k  ]=src[k];//RYZ
		//buf[k|1]=src[k|1]-src[k]+128;
		//buf[k|2]=src[k|2]-src[k]+128;
		//buf[k|3]=0xFF;
	}
	for(int ky=bh-1;ky>=0;--ky)//diff-x
	{
		for(int kx=bw-1;kx>0;--kx)
		{
			int idx=(bw*ky+kx)<<2;
			for(int kc=2;kc>=0;--kc)
				buf[idx|kc]+=128-buf[(idx-4)|kc];
		}
		if(ky>0)
		{
			for(int kc=2;kc>=0;--kc)//diff-y for first pixel in each row
				buf[(bw*ky)<<2|kc]+=128-buf[(bw*(ky-1))<<2|kc];
		}
	}

#if 0
	lodepng_encode_file("out.PNG", buf, bw, bh, LCT_RGBA, 8);
	save_mono8("ch_r_g.PNG", buf  , bw, bh, 4);
	save_mono8("ch_g.PNG", buf+1, bw, bh, 4);
	save_mono8("ch_b_g.PNG", buf+2, bw, bh, 4);
#endif

	DList list;
	dlist_init(&list, 1, 1024, 0);

	unsigned long long state=0x100000000;
	int rowlen=bw<<2;
	for(int ky=bh-1;ky>=0;--ky)
	{
		for(int kx=bw-1;kx>=0;--kx)
		{
			int idx=(bw*ky+kx)<<2;
			for(int kc=2;kc>=0;--kc)
			{
				unsigned char
					//t2=ky>=2?buf[(idx-(rowlen<<1))|kc]:128,
					t1=ky>=1?buf[(idx-rowlen)|kc]:0,
					sym=buf[idx|kc];
				//int sdev=(abs(t1-128)<<16)+g_offset;
				int sdev=(abs(t1-128)<<16)+0x12900;
				//int sdev=abs(t2-t1)+1;

				//sdev-=64;
				//sdev=sdev*sdev>>7;
				//sdev+=64;

				int erf_start=t6_calc_phi(0    , 128, sdev),
					erf_end  =t6_calc_phi(256  , 128, sdev),
					erf_curr =t6_calc_phi(sym  , 128, sdev),
					erf_next =t6_calc_phi(sym+1, 128, sdev);
				long long den=(long long)erf_end-erf_start, CDF, freq;

				CDF =((long long)(erf_curr-erf_start)<<32|den>>1)/den;
				freq=((long long)(erf_next-erf_curr )<<32|den>>1)/den;

				if(state>=((T5ANSState)freq<<32))//renorm
				{
					dlist_push_back(&list, &state, 4);
					state>>=32;
				}

				state=state/freq<<32|(CDF+state%freq);//update
			}
		}
	}
	dlist_push_back(&list, &state, 8);

	dlist_appendtoarray(&list, data);
	dlist_clear(&list);
	free(buf);
	return 1;
}

#if 0
void print_matrix_fixed(long long *matrix, int bw, int bh, int nbits)
{
	double nlevels=1<<nbits;
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
			printf("\t%lf", (double)matrix[bw*ky+kx]/nlevels);
		printf("\n");
	}
	printf("\n");
}
void impl_ref(long long *m, short dx, short dy)
{
#ifdef _DEBUG
	long long pivot;
#endif
	long long coeff;
	int mindim=dx<dy?dx:dy, it, ky, kx, npivots, kpivot;
	for(it=0, npivots=0;it<mindim;++it)//iteration
	{
		for(ky=npivots;ky<dy;++ky)//find pivot
		{
			if(m[dx*ky+it])
			{
#ifdef _DEBUG
				pivot=m[dx*ky+it];
#endif
				kpivot=ky;
				++npivots;
				break;
			}
		}
		if(ky<dy)
		{
			if(ky>it)
				for(kx=0;kx<dx;++kx)//swap rows
					coeff=m[dx*it+kx], m[dx*it+kx]=m[dx*ky+kx], m[dx*ky+kx]=coeff;
			for(++ky;ky<dy;++ky)//subtract pivot row
			{
				coeff=(m[dx*ky+it]<<16)/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*ky+kx]-=coeff*m[dx*kpivot+kx]>>16;
			}
		}
	}
}
void impl_rref(long long *m, short dx, short dy)
{
#ifdef _DEBUG
	long long pivot;
#endif
	long long coeff;
	int mindim=dx<dy?dx:dy, it, ky, kx, npivots, kpivot;
	for(it=0, npivots=0;it<mindim;++it)//iteration
	{
		kpivot=-1;
		for(ky=npivots;ky<dy;++ky)//find pivot
		{
			if(m[dx*ky+it])
			{
#ifdef _DEBUG
				pivot=m[dx*ky+it];
#endif
				kpivot=ky;
				++npivots;
				break;
			}
		}
		if(kpivot==-1)
			continue;
		if(kpivot>npivots-1)
		{
			for(kx=0;kx<dx;++kx)//swap rows
				coeff=m[dx*kpivot+kx], m[dx*kpivot+kx]=m[dx*(npivots-1)+kx], m[dx*(npivots-1)+kx]=coeff;
			kpivot=npivots-1;
		}
		for(ky=0;ky<dy;++ky)
		{
			if(ky==kpivot)//normalize pivot row
			{
				coeff=0x100000000/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
				{
					int idx=dx*kpivot+kx;
					m[idx]*=coeff;
					m[idx]>>=16;
				}
			}
			else//subtract pivot row from all other rows
			{
				coeff=(m[dx*ky+it]<<16)/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*ky+kx]-=coeff*m[dx*kpivot+kx]>>16;
			}
			//print_matrix_fixed(m, dx, dy, 16);
		}
	}
}
long long impl_det(long long *m, int dx)//m is destroyed
{
	int k, dxplus1=dx+1;
	long long result;

	//print_matrix_debug(m, dx, dx);//
	impl_ref(m, dx, dx);
	//print_matrix_debug(m, dx, dx);//

	result=m[0];//accumulate diagonal
	for(k=1;k<dx;++k)
		result=result*m[dxplus1*k]>>16;
	return result;
}
void impl_matinv(long long *m, short dx)//resize m to (dy * 2dx) temporarily,		dx==dy always
{
	int k, dy=dx, size=dx*dy;
			//print_matrix_fixed(m, dx<<1, dy, 16);
	for(k=size-dx;k>=0;k-=dx)//expand M into [M, 0]
	{
		memcpy(m+((size_t)k<<1), m+k, dx*sizeof(long long));
		memset(m+((size_t)k<<1)+dx, 0, dx*sizeof(long long));
	}
			//print_matrix_fixed(m, dx<<1, dy, 16);

	for(k=0;k<dx;++k)//add identity: [M, I]
		m[(dx<<1)*k+dx+k]=0x10000;
			//print_matrix_fixed(m, dx<<1, dy, 16);

	impl_rref(m, dx<<1, dy);//[I, M^-1]
			//print_matrix_fixed(m, dx<<1, dy, 16);

	for(k=0;k<size;k+=dx)//pack M^-1
		memcpy(m+k, m+((size_t)k<<1)+dx, dx*sizeof(long long));
			//print_matrix_fixed(m, dx<<1, dy, 16);
}
int test7_enc2(const unsigned char *buf, int bw, int bh, int x1, int y1, int x2, int y2, unsigned *cube, long long *plane, long long *line, DList *list)
{
	const int meanbits=16;//8
	int count=(x2-x1)*(y2-y1);
	long long mean[3]={0};
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			int idx=bw*ky+kx;
			mean[0]+=buf[idx<<2  ];
			mean[1]+=buf[idx<<2|1];
			mean[2]+=buf[idx<<2|2];
		}
	}
	mean[0]=(mean[0]<<meanbits)/count;
	mean[1]=(mean[1]<<meanbits)/count;
	mean[2]=(mean[2]<<meanbits)/count;

	long long cov[9]={0}, invcov[18]={0};
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			int idx=bw*ky+kx;
			int r=((int)buf[idx<<2  ]<<meanbits)-(int)mean[0],
				g=((int)buf[idx<<2|1]<<meanbits)-(int)mean[1],
				b=((int)buf[idx<<2|2]<<meanbits)-(int)mean[2];
			// r*r    r*g    r*b
			// r*g    g*g    g*b
			// r*b    g*b    b*b
			int rr=(int)((long long)r*r>>meanbits),
				gg=(int)((long long)g*g>>meanbits),
				bb=(int)((long long)b*b>>meanbits),
				rg=(int)((long long)r*g>>meanbits),
				rb=(int)((long long)r*b>>meanbits),
				gb=(int)((long long)g*b>>meanbits);
			cov[0]+=rr, cov[1]+=rg, cov[2]+=rb;
			cov[3]+=rg, cov[4]+=gg, cov[5]+=gb;
			cov[6]+=rb, cov[7]+=gb, cov[8]+=bb;
		}
	}
	for(int k=0;k<9;++k)
		cov[k]/=count;

	memcpy(invcov, cov, 9*sizeof(long long));
	impl_matinv(invcov, 3);//invert covariance

	memcpy(invcov+9, cov, 9*sizeof(long long));
	long long det=impl_det(invcov+9, 3);

	const double twopicubed=248.0502134423985614;
	long long gain=(long long)(65536./sqrt(llabs(det)*(twopicubed/65536.)));//TODO fixed sqrt

#if 0
	static int call=0;
	++call;
	printf("Covariance %d, det %lf\n", call, det/65536.);
	print_matrix_fixed(cov, 3, 3, 16);
	printf("Mean %d\n", call);
	print_matrix_fixed(mean, 1, 3, meanbits);
	printf("Invcov %d\n", call);
	print_matrix_fixed(invcov, 3, 3, 16);
#endif

	memset(cube , 0, 0x1000000*sizeof(unsigned));
	memset(plane, 0, 0x10000*sizeof(unsigned));
	memset(line , 0, 0x100*sizeof(unsigned));
	for(int color=0;color<0x1000000;++color)
	{
		if(color==0x808080)
			color=0x808080;
		int rgb[]=
		{
			((color    &255)<<meanbits)-mean[0],
			((color>> 8&255)<<meanbits)-mean[1],
			((color>>16&255)<<meanbits)-mean[2],
		};
		//f(b,g,r) = 1/sqrt((2pi)^3 * |Cov|) * exp (-(1/2) (x-mu)^T inv(Cov) (x-mu))
		long long t[]=
		{
			(invcov[0]*rgb[0]+invcov[1]*rgb[1]+invcov[2]*rgb[2])>>meanbits,
			(invcov[3]*rgb[0]+invcov[4]*rgb[1]+invcov[5]*rgb[2])>>meanbits,
			(invcov[6]*rgb[0]+invcov[7]*rgb[1]+invcov[8]*rgb[2])>>meanbits,
		};
		long long f=-(t[0]*rgb[0]+t[1]*rgb[1]+t[2]*rgb[2])>>(meanbits+1);//(...)*-0.5
		f=(long long)(65536.*exp(f/65536.));
		cube[color]=(unsigned)(f*gain+1);//no sr16 for 32 bit precision, add 1 to guard against zero frequency
	}

	//normalization
	long long rowsum=0, cdfsum=0;
	for(int color=0, kb=0;kb<256;++kb)
	{
		for(int kg=0;kg<256;++kg)
		{
			int c0=color;
			rowsum=0;
			for(int kr=0;kr<256;++kr, ++color)
				rowsum+=cube[color];

			color=c0;
			cdfsum=0;
			for(int kr=0;kr<256;++kr, ++color)
			{
				unsigned freq=(unsigned)(((unsigned long long)cube[color]<<16)/rowsum);
				cube[color]=(unsigned)cdfsum;
				cdfsum+=freq;
			}
			if(cdfsum>0x10000)
				LOG_ERROR("GB[%d, %d] CDF sum 0x%08llX > 0x10000", kg, kb, cdfsum);
			plane[kb<<8|kg]=rowsum;
		}
	}
	for(int kb=0;kb<256;++kb)
	{
		rowsum=0;
		for(int kg=0;kg<256;++kg)
			rowsum+=plane[kb<<8|kg];
		
		cdfsum=0;
		for(int kg=0;kg<256;++kg)
		{
			long long freq=((unsigned long long)plane[kb<<8|kg]<<16)/rowsum;
			plane[kb<<8|kg]=cdfsum;
			cdfsum+=freq;
		}
		if(cdfsum>0x10000)
			LOG_ERROR("B[%d] CDF sum 0x%08llX > 0x10000", kb, cdfsum);
		line[kb]=rowsum;
	}
	rowsum=0;
	for(int kb=0;kb<256;++kb)
		rowsum+=line[kb];
	cdfsum=0;
	for(int kb=0;kb<256;++kb)
	{
		long long freq=((unsigned long long)line[kb]<<16)/rowsum;
		line[kb]=cdfsum;
		cdfsum+=freq;
	}
	if(cdfsum>0x10000)
		LOG_ERROR("Master CDF sum 0x%08llX > 0x10000", cdfsum);

	for(int k=0;k<3;++k)
		dlist_push_back(list, mean, 4);
	for(int k=0;k<9;++k)
		dlist_push_back(list, cov, 4);

	unsigned state=0x10000;
	for(int ky=y2-1;ky>=y1;--ky)
	{
		for(int kx=x2-1;kx>=x1;--kx)
		{
			unsigned CDF[3], freq[3];
			int v1, v2;
			int idx=(bw*ky+kx)<<2;
			//unsigned char rgb[]={buf[idx], buf[idx|1], buf[idx|2]};

			v1=buf[idx];
			CDF[2]=(unsigned)line[v1];
			freq[2]=(unsigned)((v1<255?line[v1+1]:0x10000)-CDF[2]);
			
			v2=buf[idx|1];
			v1=v1<<8|v2;
			CDF[1]=(unsigned)plane[v1];
			freq[1]=(unsigned)((v2<255?plane[v1+1]:0x10000)-CDF[1]);
			
			v2=buf[idx|1];
			v1=v1<<8|v2;
			CDF[0]=cube[v1];
			freq[0]=(v2<255?cube[v1+1]:0x10000)-CDF[0];

			//if(!freq[0]||!freq[1]||!freq[2])
			//	LOG_ERROR("RGB freq 0x%04X 0x%04X 0x%04X", freq[0], freq[1], freq[2]);

			for(int kc=2;kc>=0;--kc)
			{
				freq[kc]+=!freq[kc];//CHEAT

				if(state>=(freq[kc]<<16))//renorm
				{
					dlist_push_back(list, &state, 2);
					state>>=16;
				}

				state=state/freq[kc]<<16|(CDF[kc]+state%freq[kc]);//update
			}
		}
	}
	dlist_push_back(list, &state, 4);
	return 1;
}
size_t test7_encode(const unsigned char *src, int bw, int bh, int is_unsigned, ArrayHandle *data)
{
	size_t res=(size_t)bw*bh, len=res<<2;
	unsigned char *buf=(unsigned char*)malloc(len);
	if(!buf)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}

	memcpy(buf, src, len);
	if(is_unsigned)
		addhalf(buf, bw, bh, 3, 4, 128);
	colortransform_ycocg_fwd((char*)buf, bw, bh);
	ArrayHandle sizes=dwt2d_gensizes(bw, bh, 3, 3, 0);
	char *temp=(char*)malloc(MAXVAR(bw, bh));
	for(int kc=0;kc<3;++kc)
		dwt2d_cdf53_fwd((char*)buf+kc, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
	addhalf(buf, bw, bh, 3, 4, 128);//mean must be in [0, 255]

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	int fw2=bw>>1, fh2=bh>>1;//floor half

	unsigned  *cube =(unsigned *)malloc((size_t)0x1000000*sizeof(unsigned));
	long long *plane=(long long*)malloc((size_t)0x10000*sizeof(long long));
	long long *line =(long long*)malloc((size_t)0x100*sizeof(long long));
	test7_enc2(buf, bw, bh, 0  , 0  , fw2, fh2, cube, plane, line, &list);
	test7_enc2(buf, bw, bh, fw2, 0  , bw , fh2, cube, plane, line, &list);
	test7_enc2(buf, bw, bh, 0  , fh2, fw2, bh , cube, plane, line, &list);
	test7_enc2(buf, bw, bh, fw2, fh2, bw , bh , cube, plane, line, &list);

	//test7_enc2(buf, bw, bh, 0  , 0  , fw2, fh2, 0, 0, 0, &list);
	//test7_enc2(buf, bw, bh, fw2, 0  , bw , fh2, 0, 0, 0, &list);
	//test7_enc2(buf, bw, bh, 0  , fh2, fw2, bh , 0, 0, 0, &list);
	//test7_enc2(buf, bw, bh, fw2, fh2, bw , bh , 0, 0, 0, &list);

	free(cube);
	free(plane);
	free(line);
	free(buf);

	dlist_appendtoarray(&list, data);

	dlist_clear(&list);
	return 1;
}
#endif


static void normalize_histogram(unsigned *srchist, int nlevels, int nsymbols, unsigned short *CDF)//hist is unsigned char due to alignment issues, but it's 16bit
{
	SortedHistInfo h[512];
	for(int k=0;k<nlevels;++k)
	{
		h[k].sym=k;
		h[k].freq=srchist[k];
	}
	for(int k=0;k<nlevels;++k)
		h[k].qfreq=((long long)h[k].freq<<16)/nsymbols;
	
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
	int sum=0;
	for(int k=0;k<nlevels;++k)
	{
		CDF[k]=sum;
		sum+=h[k].qfreq;
	}
}
	static const int point=5, offset=112;
size_t test8_encode(const unsigned char *src, int bw, int bh, int transform, ArrayHandle *data)
{
	int res=bw*bh;
	unsigned char *b2=(unsigned char*)malloc((size_t)res<<2);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(b2, src, (size_t)res<<2);
	if(transform)
		apply_transforms_fwd(b2, bw, bh);
	

	const int headsize=1<<(8-point)*3, remsize=1<<point,
		masklo=remsize-1, maskhi=(1<<(8-point))-1,
		totalhsize=headsize+remsize*3;
	
	unsigned *hist=(unsigned*)malloc((size_t)totalhsize*sizeof(unsigned));
	unsigned short *CDFs=(unsigned short*)malloc((size_t)totalhsize*sizeof(short));
	if(!hist||!CDFs)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(hist, 0, (size_t)totalhsize*sizeof(unsigned));
	unsigned *histhead=hist, *histr=histhead+headsize, *histg=histr+remsize, *histb=histg+remsize;
	for(int k=0;k<res;++k)
	{
		unsigned char r=b2[k<<2]-offset, g=b2[k<<2|1]-offset, b=b2[k<<2|2]-offset;
		int head=(b>>point&maskhi)<<((8-point)<<1)|(g>>point&maskhi)<<(8-point)|r>>point&maskhi;
		r&=masklo;
		g&=masklo;
		b&=masklo;
		++histhead[head];
		if(!head)
		{
			++histr[head<<point|r];
			++histg[head<<point|g];
			++histb[head<<point|b];
		}
	}

	unsigned short *CDFh=CDFs, *CDFr=CDFh+headsize, *CDFg=CDFr+headsize*remsize, *CDFb=CDFg+headsize*remsize;
	normalize_histogram(histhead, headsize, res, CDFh);
	normalize_histogram(histr   , remsize , res, CDFr);
	normalize_histogram(histg   , remsize , res, CDFg);
	normalize_histogram(histb   , remsize , res, CDFb);
	free(hist);

	DList list;
	dlist_init(&list, 1, 1024, 0);

	dlist_push_back(&list, 0, 4);
	dlist_push_back(&list, CDFs, (size_t)totalhsize*sizeof(short));

	double csize[4]={0};//

	unsigned state=0x10000;
	for(int k=res-1;k>=0;--k)
	{
		int idx=k<<2;
		unsigned char r=b2[idx]-offset, g=b2[idx|1]-offset, b=b2[idx|2]-offset;
		int head=(b>>point&maskhi)<<((8-point)<<1)|(g>>point&maskhi)<<(8-point)|r>>point&maskhi;
		r&=masklo;
		g&=masklo;
		b&=masklo;

		unsigned CDF[4], freq[4];
		CDF[0]=CDFh[head], freq[0]=(head<headsize-1?CDFh[head+1]:0x10000)-CDF[0];
		if(head)
		{
			CDF[1]=CDFr[r], freq[1]=(r<masklo   ?CDFr[r+1]:0x10000)-CDF[1];
			CDF[2]=CDFg[g], freq[2]=(g<masklo   ?CDFg[g+1]:0x10000)-CDF[2];
			CDF[3]=CDFb[b], freq[3]=(b<masklo   ?CDFb[b+1]:0x10000)-CDF[3];
		}
		else//bypass rare pixels
		{
			CDF[1]=r<<(16-point), freq[1]=1<<(16-point);
			CDF[2]=g<<(16-point), freq[2]=1<<(16-point);
			CDF[3]=b<<(16-point), freq[3]=1<<(16-point);
		}

		//if(!k)
		//{
		//	printf("k %d\n", k);
		//	printf("head %3d CDF 0x%04X freq 0x%04X\n", head, CDF[0], freq[0]);
		//	printf("r    %3d CDF 0x%04X freq 0x%04X\n", r,    CDF[1], freq[1]);
		//	printf("g    %3d CDF 0x%04X freq 0x%04X\n", g,    CDF[2], freq[2]);
		//	printf("b    %3d CDF 0x%04X freq 0x%04X\n", b,    CDF[3], freq[3]);
		//}

		for(int k=3;k>=0;--k)
		{
			if(!freq[k])
				LOG_ERROR("ZPS XY[%d %d]", k%bw, k/bw);

			double p=freq[k]/65536.;//
			csize[k]-=log2(p);//

			if(state>=(freq[k]<<16))//renorm
			{
				dlist_push_back(&list, &state, 2);
				state>>=16;
			}

			state=state/freq[k]<<16|(CDF[k]+state%freq[k]);//update
		}
	}
	dlist_push_back(&list, &state, 4);

	printf("head ~ %lf\n", csize[0]/((8-point)*3));//
	printf("r    ~ %lf\n", csize[1]/point);//
	printf("g    ~ %lf\n", csize[2]/point);//
	printf("b    ~ %lf\n", csize[3]/point);//

	size_t dststart=0;
	if(*data)
	{
		if(data[0]->esize!=1)
			LOG_ERROR("Invalid destination array");
		dststart=data[0]->count;
	}
	dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, &list.nobj, 4);

	dlist_clear(&list);
	free(CDFs);
	free(b2);
	return 1;
}
static int t8_dec2(unsigned *state, const unsigned char *srcstart, const unsigned char **srcptr, unsigned short *CDF, int nbits)
{
	unsigned short c;
	unsigned cdf, freq;
	int sym;
	int nlevels=1<<nbits;
	c=(unsigned short)*state;
	if(CDF)
	{
		for(sym=0;sym<nlevels-1;++sym)//find first sym where c>CDF[sym]		//linear search		FIXME use binary search
		{
			unsigned next=CDF[sym+1];
			if(next<CDF[sym])
				next=0x10000;
			if(c<next)
				break;
		}

		cdf=CDF[sym];
		freq=(sym<nlevels-1?CDF[sym+1]:0x10000)-cdf;
	}
	else
	{
		sym=c>>(16-nbits);
		cdf=sym<<(16-nbits), freq=1<<(16-nbits);
	}

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
	return sym;
}
int test8_decode(const unsigned char *data, size_t srclen, int bw, int bh, int detransform, unsigned char *buf)
{
	const int point=5, offset=112;

	const int headsize=1<<(8-point)*3, remsize=1<<point,
		masklo=remsize-1, maskhi=(1<<(8-point))-1,
		totalhsize=headsize+headsize*remsize*3,
		overhead=4+totalhsize*sizeof(short);
	
	const unsigned char *srcptr=data, *srcstart, *srcend=data+srclen;
	if(srcptr+overhead>=srcend)
	{
		LOG_ERROR("Invalid file");
		return 0;
	}

	size_t anslen=0;
	memcpy(&anslen, srcptr, 4);
	srcptr+=4;
	if(anslen<overhead)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	anslen-=overhead;
	
	unsigned short *CDFs=(unsigned short*)malloc((size_t)totalhsize*sizeof(short));
	if(!CDFs)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(CDFs, srcptr, (size_t)totalhsize*sizeof(short));
	srcptr+=(size_t)totalhsize*sizeof(short);
	unsigned short *CDFh=CDFs, *CDFr=CDFh+headsize, *CDFg=CDFr+headsize*remsize, *CDFb=CDFg+headsize*remsize;
	
	srcstart=srcptr;

	if(srcptr+anslen>data+srclen)
		LOG_ERROR("Inconsistent file");
	srcptr+=anslen;
	unsigned state;
	if(srcptr-4<srcstart)
		LOG_ERROR("ANS buffer overflow: ptr %#016llX <= start %#016llX", srcptr-data, srcstart-data);
	srcptr-=4;
	memcpy(&state, srcptr, 4);
	for(int k=0, res=bw*bh;k<res;++k)
	{
		int head=t8_dec2(&state, srcstart, &srcptr, CDFh, 8-point);
		int r   =t8_dec2(&state, srcstart, &srcptr, head?CDFr+remsize*head:0, point);
		int g   =t8_dec2(&state, srcstart, &srcptr, head?CDFg+remsize*head:0, point);
		int b   =t8_dec2(&state, srcstart, &srcptr, head?CDFb+remsize*head:0, point);

		r|=(head                &maskhi)<<point;
		g|=(head>> (8-point)    &maskhi)<<point;
		b|=(head>>((8-point)<<1)&maskhi)<<point;

		buf[k<<2  ]=r+offset;
		buf[k<<2|1]=g+offset;
		buf[k<<2|2]=b+offset;
		buf[k<<2|3]=0xFF;
	}
	if(detransform)
		apply_transforms_inv(buf, bw, bh);

	free(CDFs);
	return 1;
}

//test10
void t10_enc2(unsigned *state, unsigned cdf, unsigned freq, DList *list)
{
	if(*state>=(freq<<16))//renorm
	{
		dlist_push_back(list, state, 2);
		*state>>=16;
	}

	*state=*state/freq<<16|(cdf+*state%freq);//update
}
//int debug_px=205087;//
size_t test10_encode(const unsigned char *src, int bw, int bh, ArrayHandle *data)
{
	int res=bw*bh;
	unsigned char *b2=(unsigned char*)malloc((size_t)res<<2);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(b2, src, (size_t)res<<2);
	apply_transforms_fwd(b2, bw, bh);

	int totalhsize=256*3;
	unsigned *hist=(unsigned*)malloc((size_t)totalhsize*sizeof(unsigned));
	unsigned short *CDFs=(unsigned short*)malloc((size_t)totalhsize*sizeof(short));
	if(!hist||!CDFs)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(hist, 0, (size_t)totalhsize*sizeof(unsigned));
	unsigned *h1=hist, *h2=h1+256, *h3=h2+256;
	unsigned short *CDF1=CDFs, *CDF2=CDF1+256, *CDF3=CDF2+256;

	int count_j2=0, count_j3=0;
	for(int k=0;k<res;++k)
	{
		int r=b2[k<<2]-112, g=b2[k<<2|1]-112, b=b2[k<<2|2]-96;
		int joint1=(b>>6&3)<<6|(g>>5&7)<<3|(r>>5&7);
		++h1[joint1];
		if(!joint1)
		{
			r-=14, r&=31;
			g-=14, g&=31;
			b-=24, b&=63;
			int joint2=(b>>4&3)<<6|(g>>2&7)<<3|(r>>2&7);
			++h2[joint2];
			++count_j2;
			if(!joint2)
			{
				--r, r&=3;
				--g, r&=3;
				b-=6, b&=15;
				int joint3=(b&15)<<4|(g&3)<<2|(r&3);
				++h3[joint3];
				++count_j3;
			}
		}
	}
	normalize_histogram(h1, 256, res     , CDF1);
	normalize_histogram(h2, 256, count_j2, CDF2);
	normalize_histogram(h3, 256, count_j3, CDF3);

	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	dlist_push_back(&list, 0, 4);
	dlist_push_back(&list, CDFs, (size_t)totalhsize*sizeof(short));

#if 0
	int nc1=0, nc2=0, nc3=0, nbypass=0;//
	double cs1=0, cs2=0, cs3=0;
	long long fsum=0;
#endif

	unsigned state=0x10000;
	for(int k=res-1;k>=0;--k)
	{
		//if(k==debug_px)//
		//	k=debug_px;

		unsigned cdf[4], freq[4], nsym;
		int r=b2[k<<2]-112, g=b2[k<<2|1]-112, b=b2[k<<2|2]-96;
		int joint1=(b>>6&3)<<6|(g>>5&7)<<3|(r>>5&7), joint2=0, joint3=0;
		r-=14, r&=31;
		g-=14, g&=31;
		b-=24, b&=63;

		//predict in order of decode
		cdf[0]=CDF1[joint1], freq[0]=(joint1<255&&cdf[0]<CDF1[joint1+1]?CDF1[joint1+1]:0x10000)-cdf[0];//, nc1+=8, cs1-=log2(freq[0]/65536.), fsum+=freq[0];
		if(joint1)//bypass
		{
			cdf[1]=r<<(16-5), freq[1]=1<<(16-5);
			cdf[2]=g<<(16-5), freq[2]=1<<(16-5);
			cdf[3]=b<<(16-6), freq[3]=1<<(16-6);
			nsym=4;

			//nbypass+=5+5+6;
		}
		else
		{
			joint2=(b>>4&3)<<6|(g>>2&7)<<3|(r>>2&7);
			--r, r&=3;
			--g, r&=3;
			b-=6, b&=15;
			joint3=(b&15)<<4|(g&3)<<2|(r&3);
			cdf[1]=CDF2[joint2], freq[1]=(joint2<255&&cdf[1]<CDF2[joint2+1]?CDF2[joint2+1]:0x10000)-cdf[1];//, nc2+=8, cs2-=log2(freq[1]/65536.), fsum+=freq[1];
			if(joint2)//bypass
				cdf[2]=joint3<<(16-8), freq[2]=1<<(16-8);//, nbypass+=8;
			else
				cdf[2]=CDF3[joint3], freq[2]=(joint3<255&&cdf[2]<CDF3[joint3+1]?CDF3[joint3+1]:0x10000)-cdf[2];//, nc3+=8, cs3-=log2(freq[2]/65536.), fsum+=freq[2];
			nsym=3;
		}
		//if(k==debug_px)//
		//	printf("k %d nsym=%d J[%02X %02X %02X] RGB[%02X %02X %02X] RGB0[%02X %02X %02X]\n", k, nsym, joint1, joint2, joint3, r, g, b, b2[k<<2], b2[k<<2|1], b2[k<<2|2]);

		//encode
		for(int k2=nsym-1;k2>=0;--k2)
		{
			//unsigned s0=state;
			t10_enc2(&state, cdf[k2], freq[k2], &list);
			//if(k==debug_px)
			//	printf("sym [%d] CDF %04X freq %04X state 0x%08X -> 0x%08X\n", k2, cdf[k2], freq[k2], s0, state);
		}
	}
	dlist_push_back(&list, &state, 4);

#if 0
	printf("joint1   %14lf <- %7d CR %8lf\n", cs1/8, nc1/8, nc1/cs1);
	printf("joint2   %14lf <- %7d CR %8lf\n", cs2/8, nc2/8, nc2/cs2);
	printf("joint3   %14lf <- %7d CR %8lf\n", cs3/8, nc3/8, nc3/cs3);
	printf("bypass   %7d\n", nbypass/8);
	printf("overhead %7d\n", 256*3*2);
	printf("csize    %14lf\n", (cs1+cs2+cs3+nbypass+256*3*2)/8);
	printf("usize    %7d\n", (nc1+nc2+nc3+nbypass)/8);
	printf("CR       %14lf\n", (nc1+nc2+nc3+nbypass)/(cs1+cs2+cs3+nbypass+256*3*2));
	//printf("c1 %d c2 %d c3 %d bypass %d total %d j1 %lf j2 %lf j3 %lf fav %lf\n", nc1/8, nc2/8, nc3/8, nbypass/8, (nc1+nc2+nc3+nbypass)/8, cs1/8., cs2/8., cs3/8., fsum*8./(nc1+nc2+nc3));//
#endif

	size_t dststart=0;
	if(*data)
	{
		if(data[0]->esize!=1)
			LOG_ERROR("Invalid destination array");
		dststart=data[0]->count;
	}
	dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, &list.nobj, 4);

	dlist_clear(&list);
	free(CDFs);
	free(b2);
	return 1;
}
void t10_fill_CDF2sym(const unsigned short *CDF, unsigned char *CDF2sym)
{
	int k2=0, end=0;
	for(int sym=0;sym<256;++sym, k2=end)
	{
		if(sym==255)
			end=0x10000;
		else
		{
			end=CDF[sym+1];
			if(end<k2)
				end=0x10000;
		}

		for(;k2<end;++k2)
			CDF2sym[k2]=sym;

		if(end==0x10000)//now redundant because	'k2' is recycled 'end'		//to avoid overwriting whole CDF2sym with 0xFF with trailing zero-freq symbols
			break;
	}
}
static int t10_dec2(unsigned *state, const unsigned char *srcstart, const unsigned char **srcptr, unsigned short *CDF, unsigned char *CDF2sym, int nbits)
{
	unsigned short c;
	unsigned cdf, freq;
	int sym;
	int nlevels=1<<nbits;
	c=(unsigned short)*state;
	if(CDF)
	{
		sym=CDF2sym[c];

		cdf=CDF[sym];
		freq=(sym<nlevels-1&&cdf<CDF[sym+1]?CDF[sym+1]:0x10000)-cdf;
	}
	else
	{
		sym=c>>(16-nbits);
		cdf=sym<<(16-nbits), freq=1<<(16-nbits);
	}

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
	return sym;
}
//const unsigned char *debug_ptr=0;//
int test10_decode(const unsigned char *data, size_t srclen, int bw, int bh, unsigned char *buf)
{
	const int cdflen=3*256*sizeof(short), overhead=4+cdflen;
	int res=bw*bh;
	
	const unsigned char *srcptr, *srcstart, *srcend=data+srclen;
	if(data+overhead>=srcend)
	{
		LOG_ERROR("Invalid file");
		return 0;
	}

	size_t csize=0;
	memcpy(&csize, data, 4);
	if(csize<overhead)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	if(csize>srclen)
	{
		LOG_ERROR("Incomplete file");
		return 0;
	}
	
	unsigned short *CDFs=(unsigned short*)malloc(cdflen);
	unsigned char *CDF2sym=(unsigned char*)malloc((size_t)65536*3);
	if(!CDFs||!CDF2sym)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(CDFs, data+4, cdflen);
	unsigned short *CDF1=CDFs, *CDF2=CDF1+256, *CDF3=CDF2+256;
	unsigned char *CDF2sym1=CDF2sym, *CDF2sym2=CDF2sym1+0x10000, *CDF2sym3=CDF2sym2+0x10000;

	t10_fill_CDF2sym(CDF1, CDF2sym1);
	t10_fill_CDF2sym(CDF2, CDF2sym2);
	t10_fill_CDF2sym(CDF3, CDF2sym3);
	
	srcstart=data+overhead;
	srcptr=data+csize;

	unsigned state;
	if(srcptr-4<srcstart)
		LOG_ERROR("ANS buffer overflow: ptr %#016llX <= start %#016llX", srcptr-data, srcstart-data);
	srcptr-=4;
	memcpy(&state, srcptr, 4);
	for(int k=0, res=bw*bh;k<res;++k)
	{
		//if(k==debug_px)//
		//	k=debug_px;

		int r, g, b;
		int joint1=t10_dec2(&state, srcstart, &srcptr, CDF1, CDF2sym1, 8), joint2=0, joint3=0;
		if(joint1)//bypass
		{
			r=t10_dec2(&state, srcstart, &srcptr, 0, 0, 5);
			g=t10_dec2(&state, srcstart, &srcptr, 0, 0, 5);
			b=t10_dec2(&state, srcstart, &srcptr, 0, 0, 6);
		}
		else
		{
			joint2=t10_dec2(&state, srcstart, &srcptr, CDF2, CDF2sym2, 8);
			if(joint2)//bypass
				joint3=t10_dec2(&state, srcstart, &srcptr, 0, 0, 8);
			else
				joint3=t10_dec2(&state, srcstart, &srcptr, CDF3, CDF2sym3, 8);

			r=joint3   & 3;
			g=joint3>>2& 3;
			b=joint3>>4&15;

			r=(joint2   &7)<<2|(r+1)& 3;
			g=(joint2>>3&7)<<2|(g+1)& 3;
			b=(joint2>>6&3)<<4|(b+6)&15;
		}
		buf[k<<2  ]=((joint1   &7)<<5|(r+14)&31)+112;
		buf[k<<2|1]=((joint1>>3&7)<<5|(g+14)&31)+112;
		buf[k<<2|2]=((joint1>>6&3)<<6|(b+24)&63)+ 96;
		buf[k<<2|3]=0xFF;

		//if(debug_ptr&&memcmp(buf+((size_t)k<<2), debug_ptr+((size_t)k<<2), 4))//
		//	LOG_ERROR("Error at pixel number %d", k);//
	}
	apply_transforms_inv(buf, bw, bh);

	free(CDFs);
	return 1;
}