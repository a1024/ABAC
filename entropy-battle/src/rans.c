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
				p->inv_freq=0xFFFFFFFFFFFFFFFF;
			//	p->shift=0;
			//	p->inv_freq=0xFFFFFFFF;
				p->bias=sum+ANS_L-1;
			}
			else
			{
				unsigned long long dummy=0;
				p->inv_freq=_udiv128(1, 0, p->freq, &dummy);
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

static const int tag_rans4='A'|'N'<<8|'0'<<16|'4'<<24;
void umul128_64_trunc128(unsigned long long *res, const unsigned long long *a, const unsigned long long b)
{
	unsigned long long t=0, t2;
	res[0]=_umul128(a[0], b, &t);
	res[1]=_umul128(a[1], b, &t2)+t;
}
void umul128_64(unsigned long long *res, const unsigned long long *a, const unsigned long long b)
{
	unsigned long long t=0;
	res[0]=_umul128(a[0], b, &t);
	res[1]=_umul128(a[1], b, res+2)+t;
	res[2]+=res[1]<t;

	//unsigned long long lo[2]={0}, hi[2]={0};
	//lo[0]=_umul128(a[0], b, lo+1);
	//hi[0]=_umul128(a[1], b, hi+1);
	//res[0]=lo[0];
	//res[1]=lo[1]+hi[0];
	//res[2]=hi[1]+(res[1]<lo[1]);
}
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
			//state+=(((long long)state*p->inv_freq>>32)>>p->shift)*p->neg_freq+p->bias;//Ryg's division-free rANS encoder	https://github.com/rygorous/ryg_rans/blob/master/rans_byte.h
#if 0
			unsigned s1=state, s2;
			unsigned long long q[2]={0}, xfif[2];
			q[0]=_umul128(s1, p->inv_freq, q+1);
			umul128_64_trunc128(xfif, q, p->freq);
			q[1]+=xfif[1]<s1;
			//qhi+=qlo>>63;
			s2=s1+(unsigned)(q[1]*p->neg_freq)+p->bias;
#endif
			//unsigned s1=state, s2=s1+(((long long)s1*p->inv_freq>>32)>>p->shift)*p->neg_freq+p->bias;
			state=state/p->freq<<ANS_PROB_BITS|(p->CDF+state%p->freq);
			//if(s2!=state)
			//	LOG_ERROR("Debug break");
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


void normalize_histogram(unsigned *srchist, int nlevels, int nsymbols, unsigned short *CDF)//hist is unsigned char due to alignment issues, but it's 16bit
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


void t14_calc_sdev(const unsigned char *b2, int iw, int ih, int lgblockdim, int x1, int y1, unsigned char *sdev)
{
	int blockdim=1<<lgblockdim;
	int xend=MINVAR(x1+blockdim, iw),
		yend=MINVAR(y1+blockdim, ih), count=(xend-x1)*(yend-y1);
	unsigned long long variance[3]={0};
	//memset(sdev, 0, 3*sizeof(long long));
	for(int ky=y1;ky<yend;++ky)
	{
		for(int kx=x1;kx<xend;++kx)
		{
			int idx=(iw*ky+kx)<<2;
			for(int kc=0;kc<3;++kc)
			{
				int val=b2[idx|kc]-128;
				variance[kc]+=val*val;
			}
		}
	}
	unsigned temp[]=
	{
		(unsigned)ceil(sqrt((double)(variance[0]<<2)/count)),//8 bit
		(unsigned)ceil(sqrt((double)(variance[1]<<2)/count)),
		(unsigned)ceil(sqrt((double)(variance[2]<<2)/count)),
	};
	if(!temp[0])//
		temp[0]=0;
	if(!temp[1])//
		temp[1]=0;
	if(!temp[2])//
		temp[2]=0;

	for(int kc=0;kc<3;++kc)
	{
		if(temp[kc]>255)
			sdev[kc]=255;
		else
			sdev[kc]=temp[kc];
	}
	sdev[3]=0;

	//sdev[0]=(sdev[0]<<16)/count;
	//sdev[1]=(sdev[1]<<16)/count;
	//sdev[2]=(sdev[2]<<16)/count;
	//sdev[0]=(unsigned long long)(sqrt(sdev[0]/65536.)*65536.);//7.16 bit
	//sdev[1]=(unsigned long long)(sqrt(sdev[1]/65536.)*65536.);
	//sdev[2]=(unsigned long long)(sqrt(sdev[2]/65536.)*65536.);

	//int sh=16-(lgblockdim<<1);
	//if(sh>0)//divide by block pixel count, 16 bit fixed precision
	//{
	//	variance[0]<<=sh;
	//	variance[1]<<=sh;
	//	variance[2]<<=sh;
	//}
	//else
	//{
	//	variance[0]>>=-sh;
	//	variance[1]>>=-sh;
	//	variance[2]>>=-sh;
	//}
}
int t14_sdev_bilinear(int xtail, int ytail, unsigned char TL, unsigned char TR, unsigned char BL, unsigned char BR)
{
	int top, bot;
	top=(TL<<16)+(TR-TL)*xtail;
	bot=(BL<<16)+(BR-BL)*xtail;
	top=top+(int)((long long)(bot-top)*ytail>>16);
	return top;
}
int t14_enc2(unsigned char *b2, int bw, int bh, int lgblockdim, ArrayHandle *data)
{
	apply_transforms_fwd(b2, bw, bh);
	DList list;
	dlist_init(&list, 1, 1024, 0);
	
	int blockdim=1<<lgblockdim;
	int w2=(bw+blockdim-1)/blockdim,
		h2=(bh+blockdim-1)/blockdim, res2=w2*h2;
	unsigned char *sdev=(unsigned char*)malloc((size_t)res2<<2);
	if(!sdev)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	
	memset(sdev, 0xFF, (size_t)res2<<2);
	for(int ky=0;ky<h2;++ky)
	{
		for(int kx=0;kx<w2;++kx)//for each block
		{
			t14_calc_sdev(b2, bw, bh, lgblockdim, blockdim*kx, blockdim*ky, sdev+(((size_t)w2*ky+kx)<<2));
			
			//if(!sdev[(w2*ky+kx)<<2])//
			//	sdev[(w2*ky+kx)<<2]=0;
			//if(!sdev[(w2*ky+kx)<<2|1])//
			//	sdev[(w2*ky+kx)<<2|1]=0;
			//if(!sdev[(w2*ky+kx)<<2|2])//
			//	sdev[(w2*ky+kx)<<2|2]=0;

			//if(sdev[(w2*ky+kx)<<2]==0xFF||sdev[(w2*ky+kx)<<2|1]==0xFF||sdev[(w2*ky+kx)<<2|2]==0xFF)//
			//if(!sdev[(w2*ky+kx)<<2]||!sdev[(w2*ky+kx)<<2|1]||!sdev[(w2*ky+kx)<<2|2])//
			//	LOG_ERROR("Debug break");
		}
	}
	unsigned state=0x10000;
	unsigned mask=(1<<lgblockdim)-1, half=1<<(lgblockdim-1);
	int sh=16-lgblockdim;
	for(int ky=bh-1;ky>=0;--ky)
	{
		for(int kx=bw-1;kx>=0;--kx)
		{
			int kx2=kx>>lgblockdim, ky2=ky>>lgblockdim;
			int idx=(bw*ky+kx)<<2;
			
			kx2-=kx<(int)(kx&~mask|half);//get current quartet coords
			ky2-=ky<(int)(ky&~mask|half);
			int x1=kx2>0?kx2-1:0;
			int y1=ky2>0?ky2-1:0;
			int x2=kx2<w2-1?kx2+1:w2-1;
			int y2=ky2<h2-1?ky2+1:h2-1;
			unsigned quad[]=
			{
				((unsigned*)sdev)[w2*y1+x1], ((unsigned*)sdev)[w2*y1+x2],
				((unsigned*)sdev)[w2*y2+x2], ((unsigned*)sdev)[w2*y2+x2],
			};
			//kx2+=kx>=(int)(kx&~mask|half);//get current quartet coords
			//ky2+=ky>=(int)(ky&~mask|half);
			//kx2-=kx2>=w2;
			//ky2-=ky2>=h2;
			//decx=kx2>0;
			//decy=ky2>0;
			//unsigned quad[]=
			//{
			//	((unsigned*)sdev)[w2*(ky2-decy)+kx2-decx], ((unsigned*)sdev)[w2*(ky2-decy)+kx2],
			//	((unsigned*)sdev)[w2* ky2      +kx2-decx], ((unsigned*)sdev)[w2* ky2      +kx2],
			//};
			const unsigned char *comp=(const unsigned char*)quad;
			unsigned xtail=(kx+half)&mask, ytail=(ky+half)&mask;
			if(sh>0)
			{
				xtail<<=sh;
				ytail<<=sh;
			}
			else
			{
				xtail>>=-sh;
				ytail>>=-sh;
			}
			for(int kc=2;kc>=0;--kc)
			{
				int sdev=t14_sdev_bilinear(xtail, ytail, comp[kc], comp[kc+4], comp[kc+8], comp[kc+12]);//calculate sdev with bilinear interpolation

				double v2[]={comp[kc], comp[kc+4], comp[kc+8], comp[kc+12]};
				double top=v2[0]+(v2[1]-v2[0])*(xtail/65536.);
				double bot=v2[2]+(v2[3]-v2[2])*(xtail/65536.);
				double mid=top+(bot-top)*(ytail/65536.);
				if(fabs(mid-sdev/65536.)>0x2p-16)
					LOG_ERROR("Debug break");

				unsigned char sym=b2[idx|kc];
				if(sdev)
				{
					long long start=-0x8000000000/sdev, end=0x7F00000000/sdev, curr=((long long)(sym-128)<<32)/sdev, next=((long long)(sym+1-128)<<32)/sdev;
					start=error_func_p32(start)<<8;
					end=error_func_p32(end)<<8|0xFF;
					curr=error_func_p32(curr)<<8|sym;
					next=(error_func_p32(next)<<8)+((long long)sym+1);
					long long den=end-start;
					if(!den)
						LOG_ERROR("Division by zero");
					unsigned cdf=(unsigned)(((curr-start)<<16)/den), freq=(unsigned)(((next-curr)<<16)/den);
					//unsigned cdf=(unsigned)(((long long)(curr-start)<<16)/den), freq=(unsigned)(((long long)(next-curr)<<16)/den);
					
					if(!freq)
						LOG_ERROR("Division by zero");

					if(state>=(freq<<16))//renorm
					{
						dlist_push_back(&list, &state, 2);
						state>>=16;
					}
					state=state/freq<<16|(cdf+state%freq);//update
				}
				//if not sdev (const 128) don't encode block
			}
		}
	}

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
	//if(w2>1||h2>1)
	//	t14_enc2(sdev, w2, h2, lgblockdim, data);
	//else
	//	ARRAY_APPEND(*data, sdev, (size_t)res2<<2, 1, 0);
	free(sdev);
	return 1;
}
size_t test14_encode(const unsigned char *src, int bw, int bh, int lgblockdim, ArrayHandle *data)
{
	int res=bw*bh;
	unsigned char *b2=(unsigned char*)malloc((size_t)res<<2);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(b2, src, (size_t)res<<2);
	t14_enc2(b2, bw, bh, lgblockdim, data);
	free(b2);
	return 1;
}


//	#define T16_ENABLE_CONF

int get_mean_p16(const unsigned char *buf, int bw, int kc, int x1, int x2, int y1, int y2)
{
	long long sum=0;
	int res=(x2-x1)*(y2-y1);
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
			sum+=buf[(bw*ky+kx)<<2|kc];
	}
	int mean=(int)(((long long)sum<<16)/res);
	return mean;
}
int get_conf_p16(const unsigned char *buf, int bw, int kc, int x1, int x2, int y1, int y2, int mean, int *trivial)
{
	long long sum=0;
	int res=(x2-x1)*(y2-y1);
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			long long x=((long long)buf[(bw*ky+kx)<<2|kc]<<16)-mean;
			sum+=x*x>>16;
		}
	}
	if(!sum)
	{
		*trivial=1;
		return 0;
	}
	sum=(sum<<1)/res;//res is an integer, need variance*2 inside sqrt
	int conf=(int)(65536/sqrt(sum/65536.));//TODO don't use floating point
	return conf;
}
void print_rep(int count, char c)
{
	for(int k2=0;k2<count;++k2)
		printf("%c", c);
}
void print_CDF(const unsigned *CDF, const unsigned char *b2, int bw, int bh, int kc, int x1, int x2, int y1, int y2)
{
	const int graphwidth=32;
	int res=bw*bh;
#if 1
	int *CDF0=(int*)malloc(256*sizeof(int));
	if(!CDF0)
	{
		LOG_ERROR("Allocation error");
		return;
	}
	memset(CDF0, 0, 256*sizeof(int));
	int count=(x2-x1)*(y2-y1);
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			unsigned char sym=b2[(bw*ky+kx)<<2];
			++CDF0[sym];
		}
	}
	int sum=0;
	for(int sym=0;sym<256;++sym)
	{
		int freq=CDF0[sym];
		freq=(int)(((long long)freq*0xFF00)/count)+1;
		CDF0[sym]=sum;
		sum+=freq;
	}
#endif

	int mean=get_mean_p16(b2, bw, kc, x1, x2, y1, y2);
	//int mean=128<<16;
	int bypass=0;
	int conf=get_conf_p16(b2, bw, kc, x1, x2, y1, y2, mean, &bypass);

	//conf*=2.5;//
	//conf<<=1;//
	
	//printf("ACTUAL CDF,  STATIC ALPHA CDF,  CONF CDF\n");
	double entropy[3]={0};

	int start=(int)((long long)(0-mean)*conf>>16),
		end=(int)((long long)((256<<16)-mean)*conf>>16);
	int
		estart=error_func_p16(start),
		eend=error_func_p16(end);
	for(int k=0;k<256;++k)
	{
		int x=(int)((long long)((k<<16)-mean)*conf>>16);
		int e=error_func_p16(x);
		int den=eend-estart;
		int cdf2=(int)(((long long)(e-estart)*0xFF00)/den)+k;//guard

		int x2=(int)((long long)(((k+1)<<16)-mean)*conf>>16);
		int e2=error_func_p16(x2);
		int f2=(int)(((long long)(e-estart)*0xFF00)/den)+1;
		
#if 1
		printf("%3d", k);

		int c0=CDF0[k]*graphwidth>>16;
		printf(" %04X ", CDF0[k]);
		print_rep(c0, '*');
		print_rep(graphwidth-c0, ' ');

		int c1=CDF[k]*graphwidth>>16;
		printf(" %04X ", CDF[k]);
		print_rep(c1, '*');
		print_rep(graphwidth-c1, ' ');

		int c2=cdf2*graphwidth>>16;
		printf(" %04X ", cdf2);
		print_rep(c2, '*');
		print_rep(graphwidth-c2, ' ');

		printf("\n");
#endif

		double p0=((k<255?CDF0[k+1]:0x10000)-CDF0[k])/65536.;
		if(p0)
			entropy[0]-=p0*log2(p0);
		double p1=((k<255?CDF[k+1]:0x10000)-CDF[k])/65536.;
		if(p1)
			entropy[1]-=p1*log2(p1);
		double p2=f2/65536.;
		if(p2)
			entropy[2]-=p2*log2(p2);
	}
	printf("BPP %lf\t%lf\t%lf\n", entropy[0], entropy[1], entropy[2]);
	//printf("CR %lf\t%lf\t%lf\n", 8/entropy[0], 8/entropy[1], 8/entropy[2]);
	//printf("mean 0x%08X conf 0x%08X\n", mean, conf);
	//printf("Actual         %lf\n", 8/entropy[0]);
	//printf("Static Alpha   %lf\n", 8/entropy[1]);
	//printf("Conf           %lf\n", 8/entropy[2]);
	//printf("\n");
	free(CDF0);
}
void t16_prepblock(const unsigned char *b2, const unsigned short *CDF, int bw, int bh, int kc, int bx, int by, int alpha, int blockw, int blockh, int margin, unsigned *CDF2, int *xend, int *yend)
{
	int kx=bx*blockw, ky=by*blockh;

	*yend=ky+blockh<=bh?ky+blockh:bh;
	*xend=kx+blockw<=bw?kx+blockw:bw;
				
	int overflow=0;//CDF overflow can happen only once
	if(!bx&&!by)//first block has zero alpha
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
	else
	{
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

#ifdef T16_ENABLE_CONF
		double sdev=0;
#endif

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
					//int inc=(kx2>>1)+1;
					//int inc=(kx2*3>>2)+1;
					//int inc=kx2<<1|1;
					//int inc=1;
					//int inc=kx2*kx2+1;
					//int inc=0x10000/(kx2+1);

					if(!inc)
						LOG_ERROR("Zero inc");

					CDF2[sym]+=inc;
					count2+=inc;
					
#ifdef T16_ENABLE_CONF
					sym-=128;
					sdev+=inc*sym*sym;
#endif

					//if(count2<inc)
					//	LOG_ERROR("OVERFLOW");
				}
			}
			//count2+=(blocksize*(blocksize+1)>>1)**ysize;
			//count2+=blocksize**ysize;
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
					//int inc=(ky2>>1)+1;
					//int inc=(ky2*3>>2)+1;
					//int inc=ky2<<1|1;
					//int inc=1;
					//int inc=ky2*ky2+1;
					//int inc=0x10000/(ky2+1);
					
					if(!inc)
						LOG_ERROR("Zero inc");

					CDF2[sym]+=inc;
					count2+=inc;
					
#ifdef T16_ENABLE_CONF
					sym-=128;
					sdev+=inc*sym*sym;
#endif

					if(count2<inc)
						LOG_ERROR("OVERFLOW");
				}
			}
			//count2+=(blocksize*(blocksize+1)>>1)**xsize;
			//count2+=blocksize**xsize;
		}
		if(left<kx&&top<ky)//if topleft block is available
		{
			for(int ky2=top;ky2<ky;++ky2)
			{
				for(int kx2=left;kx2<kx;++kx2)//for each pixel
				{
					unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];
					int dist=kx-kx2+ky-ky2;
					//int dist=MAXVAR(kx-kx2, ky-ky2);
					
					if(dist<0||dist>(margin<<1))
						LOG_ERROR("Wrong distance");
					
					int inc=(margin<<1|1)-dist;
					
					if(!inc)
						LOG_ERROR("Zero inc");

					CDF2[sym]+=inc;
					count2+=inc;
					
#ifdef T16_ENABLE_CONF
					sym-=128;
					sdev+=inc*sym*sym;
#endif
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
					
#ifdef T16_ENABLE_CONF
					sym-=128;
					sdev+=inc*sym*sym;
#endif
				}
			}
		}
#if 0
		int xoffset, yoffset;
		int xsize2, ysize2, count2;
		if(!bx)
		{
			if(!by)
				xoffset=0, yoffset=0;
			else
				xoffset=0, yoffset=blocksize;
		}
		else
			xoffset=blocksize, yoffset=0;
		ysize2=ky-yoffset+blocksize<=bh?blocksize:bh-(ky-yoffset);
		xsize2=kx-xoffset+blocksize<=bw?blocksize:bw-(kx-xoffset);
		count2=xsize2*ysize2;

		for(int ky2=0;ky2<ysize2;++ky2)
		{
			for(int kx2=0;kx2<xsize2;++kx2)//for each pixel
			{
				unsigned char sym=b2[(bw*(ky-yoffset+ky2)+kx-xoffset+kx2)<<2|kc];
				++h2[sym];
			}
		}
#endif
#ifdef T16_ENABLE_CONF
		sdev/=count2;
		sdev=sqrt(sdev);
		double conf=1/(sqrt(2)*sdev);
#endif

		int sum=0;
		for(int sym=0;sym<256;++sym)
		{
			int cdf1=!overflow?CDF[sym]:0x10000;
			if(sym<255)
				overflow|=cdf1>CDF[sym+1];
			int f1=(sym<255&&!overflow?CDF[sym+1]:0x10000)-cdf1;

			int f2=(int)(((long long)CDF2[sym]<<16)/count2);//normalize
			
#ifdef T16_ENABLE_CONF
			double start=erf((0-128)*conf), end=erf((256-128)*conf), x=erf((sym-128)*conf), x2=erf((sym+1-128)*conf);
			int f3=(int)((x2-x)*0x10000/(end-start));
#endif

			//if(f1||f2)//
			//	printf("");

			//int freq=f1+(int)(((long long)f3-f1)*alpha>>16);//blend
			//freq=freq+(int)(((long long)f2-freq)*alpha>>16);//blend

			//int freq=f2;
			int freq=f1+(int)(((long long)f2-f1)*alpha>>16);//blend

#ifdef T16_ENABLE_CONF
			freq=freq+(int)(((long long)f3-freq)*alpha>>16);//blend2
#endif

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
}
size_t test16_encode(const unsigned char *src, int bw, int bh, int alpha, int *blockw, int *blockh, int *margin, ArrayHandle *data, int loud, int *csizes)
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
	unsigned short *CDF=(unsigned short*)malloc((size_t)totalhsize*sizeof(short));
	unsigned *h2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!hist||!CDF||!h2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(hist, 0, (size_t)totalhsize*sizeof(unsigned));
	for(int kc=0;kc<3;++kc)
	{
		for(int k=0;k<res;++k)
		{
			unsigned char sym=b2[k<<2|kc];
			++hist[kc<<8|sym];
		}
	}
	normalize_histogram(hist, 256, res, CDF);//this is just to pack the histogram, CDF is renormalized again with ramp guard
	normalize_histogram(hist+256, 256, res, CDF+256);
	normalize_histogram(hist+512, 256, res, CDF+512);

	DList list;
	dlist_init(&list, 1, 1024, 0);

	int ansbookmarks[3]={0};
	dlist_push_back(&list, 0, 12);
	dlist_push_back(&list, CDF, 768*sizeof(short));

	for(int kc=0;kc<3;++kc)
	{
		int bxcount=(bw+blockw[kc]-1)/blockw[kc],
			bycount=(bh+blockh[kc]-1)/blockh[kc];

		unsigned state=0x10000;
		for(int by=bycount-1;by>=0;--by)
		{
			int ky=by*blockh[kc];
			for(int bx=bxcount-1;bx>=0;--bx)//for each block
			{
				//if(kc==0&&bx==0&&by==0)
				//	kc=0;

				int kx=bx*blockw[kc];
				int xend=0, yend=0;
				t16_prepblock(b2, CDF+((size_t)kc<<8), bw, bh, kc, bx, by, alpha, blockw[kc], blockh[kc], margin[kc], h2, &xend, &yend);

				//if(kc==0&&bx==0&&by==0)
				//	print_CDF(h2, b2, bw, bh, kc, kx, xend, ky, yend);

				//encode block
				for(int ky2=yend-1;ky2>=ky;--ky2)
				{
					for(int kx2=xend-1;kx2>=kx;--kx2)//for each pixel
					{
						unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];

						int cdf=h2[sym], freq=h2[sym+1]-cdf;
#if 0
						int cdf=CDF[kc<<8|sym], freq=(sym<255&&cdf<CDF[kc<<8|(sym+1)]?CDF[kc<<8|(sym+1)]:0x10000)-cdf, cdf2, f2;
						if(bx||by)
						{
							cdf2=h2[sym], f2=h2[sym+1]-cdf2;
							cdf2=(cdf2<<16)/h2[256];
							f2=(f2<<16)/h2[256];
							cdf=cdf+(int)(((long long)cdf2-cdf)*alpha>>16);
							freq=freq+(int)(((long long)f2-freq)*alpha>>16);
							cdf=(cdf*0xFF00>>16)+sym;
							freq=(freq*0xFF00>>16)+1;
						}
#endif

						//if(kc==2&&ky2==0&&kx2==0)//
						//	printf("sym 0x%02X cdf 0x%04X freq 0x%04X\n", sym, cdf, freq);

						if(!freq)
							LOG_ERROR("ZPS");
						
						if(state>=(unsigned)(freq<<16))//renorm
						{
							dlist_push_back(&list, &state, 2);
							state>>=16;
						}
						state=state/freq<<16|(cdf+state%freq);//update
					}
				}
			}
		}
		dlist_push_back(&list, &state, 4);
		ansbookmarks[kc]=(int)list.nobj;
	}
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, ansbookmarks, 12);
	
	int overhead=12+(int)(totalhsize*sizeof(short));
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
		printf("alpha 0x%04X\n", alpha);
		printf("Total    %7d  %lf\n", ansbookmarks[2], 3.*res/ansbookmarks[2]);
		printf("Overhead %7d\n", overhead);
		printf("Red      %7d  %lf  %dx%d  M %d\n", ch[0], (double)res/ch[0], blockw[0], blockh[0], margin[0]);
		printf("Green    %7d  %lf  %dx%d  M %d\n", ch[1], (double)res/ch[1], blockw[1], blockh[1], margin[1]);
		printf("Blue     %7d  %lf  %dx%d  M %d\n", ch[2], (double)res/ch[2], blockw[2], blockh[2], margin[2]);
	}

	dlist_clear(&list);
	free(b2);
	free(hist);
	free(CDF);
	free(h2);
	return 1;
}
int threeway_uint32(const void *p1, const void *p2)
{
	unsigned v1=*(const unsigned*)p1, v2=*(const unsigned*)p2;
	return (v1>v2)-(v1<v2);
}
unsigned char *debug_ptr=0;
int test16_decode(const unsigned char *data, size_t srclen, int bw, int bh, int alpha, int *blockw, int *blockh, int *margin, unsigned char *buf)
{
	const int cdflen=768LL*sizeof(short), overhead=12LL+cdflen;
	int res=bw*bh;
	
	const unsigned char *srcptr, *srcstart, *srcend=data+srclen;
	if(data+overhead>=srcend)
	{
		LOG_ERROR("Invalid file");
		return 0;
	}

	unsigned ansbookmarks[3];
	memcpy(ansbookmarks, data, 12);
	if(ansbookmarks[2]<(unsigned)overhead)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	if(ansbookmarks[2]>srclen)
	{
		LOG_ERROR("Incomplete file");
		return 0;
	}
	
	unsigned short *CDF=(unsigned short*)malloc(cdflen);
	unsigned *h2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF||!h2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(CDF, data+12, cdflen);
	
	for(int kc=0;kc<3;++kc)
	{
		int bxcount=(bw+blockw[kc]-1)/blockw[kc],
			bycount=(bh+blockh[kc]-1)/blockh[kc];

		unsigned state;
		srcptr=data+ansbookmarks[kc];
		srcstart=kc?data+ansbookmarks[kc-1]:data+overhead;
		srcptr-=4;
		if(srcptr<srcstart)
			LOG_ERROR("ANS buffer overflow");
		memcpy(&state, srcptr, 4);
		
		for(int by=0;by<bycount;++by)
		{
			int ky=by*blockh[kc];
			for(int bx=0;bx<bxcount;++bx)//for each block
			{
				//if(kc==0&&bx==0&&by==0)
				//	kc=0;

				int kx=bx*blockw[kc];
				int xend=0, yend=0;
				t16_prepblock(buf, CDF+((size_t)kc<<8), bw, bh, kc, bx, by, alpha, blockw[kc], blockh[kc], margin[kc], h2, &xend, &yend);
				for(int ky2=ky;ky2<yend;++ky2)
				{
					for(int kx2=kx;kx2<xend;++kx2)//for each pixel
					{
						unsigned c=(unsigned short)state;
						size_t sym=0;
						int found=binary_search(h2, 257, sizeof(int), threeway_uint32, &c, &sym);
						sym-=!found;//binary_search gives insertion index

						//if(sym!=debug_ptr[(bw*ky2+kx2)<<2|kc])//
						//	LOG_ERROR("");

						buf[(bw*ky2+kx2)<<2|kc]=(unsigned char)sym;

						unsigned cdf=h2[sym], freq=h2[sym+1]-cdf;

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
	free(CDF);
	free(h2);
	for(int k=0;k<res;++k)//set alpha
		buf[k<<2|3]=0xFF;
	apply_transforms_inv(buf, bw, bh);
	return 1;
}

#if 0
double t17_estimatecsize(const unsigned char *src, int bw, int kc, int x1, int x2, int y1, int y2, int mean, int conf)
{
	double csize=0;
	int start=(int)((long long)(0-mean)*conf>>16),
		end=(int)((long long)((256<<16)-mean)*conf>>16);
	int
		estart=error_func_p16(start),
		eend=error_func_p16(end);
	int den=eend-estart;
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			int sym=src[(bw*ky+kx)<<2|kc];
			
			int curr=(int)((long long)(( sym   <<16)-mean)*conf>>16),
				next=(int)((long long)(((sym+1)<<16)-mean)*conf>>16);
			int e1=error_func_p16(curr),
				e2=error_func_p16(next);
			int freq=(int)(((long long)(e2-e1)<<16)/den);
			freq=((unsigned)(freq*0xFF00)>>16)+1;//guard

			double p=freq/65536., bitsize=-log2(p);
			csize+=bitsize;
		}
	}
	csize/=8;
	return csize;
}
void test17_saveconf(const unsigned char *src, int bw, int bh, int blocksize)
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

	int bxcount=(bw+blocksize-1)/blocksize,
		bycount=(bh+blocksize-1)/blocksize;
	int confsize=bxcount*bycount;
	unsigned char *bufsize=(unsigned char*)malloc((size_t)confsize<<2);
	//unsigned short *bufconf=(unsigned short*)malloc(confsize*sizeof(short));
	if(!bufsize)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	double csize=0;
	for(int kc=0;kc<3;++kc)
	{
		double ch_size=0;
		for(int by=0;by<bycount;++by)
		{
			int ky=by*blocksize;
			int yend=(by+1)*blocksize;
			if(yend>bh)
				yend=bh;
			for(int bx=0;bx<bxcount;++bx)//for each block
			{
				int kx=bx*blocksize;
				int xend=(bx+1)*blocksize;
				if(xend>bw)
					xend=bw;
				int trivial=0;
				unsigned conf=get_conf_p16(b2, bw, kc, kx, xend, ky, yend, 128<<16, &trivial);

				if(conf<0||conf>0xFFFF)
					LOG_ERROR("Strange confidence 0x%04X", conf);
				double bsize=t17_estimatecsize(b2, bw, kc, kx, xend, ky, yend, 128<<16, conf);

				bsize+=2;//overhead

				int bypasssize=(xend-kx)*(yend-ky);
				int size=(int)(bsize*255/bypasssize);
				if(size>255)
					size=255;
				bufsize[(bxcount*by+bx)<<2|kc]=(unsigned char)size;

				printf(" %.1lf", bsize);
				ch_size+=bsize;

				//if(conf>0xFFFE)
				//	conf=0xFFFE;
				//if(trivial)
				//	conf=0xFFFF;
				//bufconf[bxcount*by+bx]=conf;
			}
			printf("\n");
		}
		printf("csize %lf CR %lf\n\n", ch_size, res/ch_size);
		csize+=ch_size;
		//char c='X';
		//switch(kc)
		//{
		//case 0:c='R';break;
		//case 1:c='G';break;
		//case 2:c='B';break;
		//}
		//snprintf(g_buf, G_BUF_SIZE, "test17-%d%c.PNG", kc, c);
		//save_16bit(g_buf, bufconf, 0, bxcount, bycount, 1, 0, 0, 0);
	}
	for(int k=0;k<confsize;++k)
		bufsize[k<<2|3]=0xFF;
	lodepng_encode_file("test17-size.PNG", bufsize, bxcount, bycount, LCT_RGBA, 8);
	free(b2);
	free(bufsize);
	//free(bufconf);

	printf("csize %lf CR %lf\n\n", csize, res*3/csize);
}
#endif
#if 1
void t17_prepblock(const unsigned char *b2, unsigned short *CDF, int conf, int bw, int bh, int kc, int bx, int by, int alpha, int blocksize, unsigned *h2, int *xsize, int *ysize)
{
	int kx=bx*blocksize, ky=by*blocksize;

	*ysize=ky+blocksize<=bh?blocksize:bh-ky;
	*xsize=kx+blocksize<=bw?blocksize:bw-kx;
				
	int overflow=0;//CDF overflow can happen only once
	if(!bx&&!by)//first block has zero alpha
	{
		if(CDF)
		{
			for(int sym=0;sym<256;++sym)
			{
				if(overflow)
					h2[sym]=0xFF00|sym;
				else
				{
					int cdf=CDF[sym];
					h2[sym]=((unsigned)(cdf*0xFF00)>>16)+sym;
					if(sym<255)
						overflow|=cdf>CDF[sym+1];
				}
			}
			h2[256]=0x10000;
		}
		else
		{
			int mean=128<<16;
			int start=(int)((long long)(0-mean)*conf>>16),
				end=(int)((long long)((256<<16)-mean)*conf>>16);
			int
				estart=error_func_p16(start),
				eend=error_func_p16(end);
			int den=eend-estart;
			int sum=0;
			for(int sym=0;sym<256;++sym)
			{
				int curr=(int)((long long)(( sym   <<16)-mean)*conf>>16),
					next=(int)((long long)(((sym+1)<<16)-mean)*conf>>16);
				int e1=error_func_p16(curr),
					e2=error_func_p16(next);
				int freq=(int)(((long long)(e2-e1)<<16)/den);
				freq=((unsigned)(freq*0xFF00)>>16)+1;//guard
				h2[sym]=sum;
				sum+=freq;
			}
			h2[256]=0x10000;
		}
	}
	else
	{
		memset(h2, 0, 256*sizeof(unsigned));

		int count2=0;
		if(bx)//if left block available
		{
			for(int ky2=0;ky2<*ysize;++ky2)
			{
				for(int kx2=0;kx2<blocksize;++kx2)//for each pixel
				{
					unsigned char sym=b2[(bw*(ky+ky2)+kx-blocksize+kx2)<<2|kc];
					++h2[sym];
				}
			}
			count2+=blocksize**ysize;
		}
		if(by)//if top block available
		{
			for(int ky2=0;ky2<blocksize;++ky2)
			{
				for(int kx2=0;kx2<*xsize;++kx2)//for each pixel
				{
					unsigned char sym=b2[(bw*(ky-blocksize+ky2)+kx+kx2)<<2|kc];
					++h2[sym];
				}
			}
			count2+=blocksize**xsize;
		}

		int sum=0;
		if(CDF)
		{
			for(int sym=0;sym<256;++sym)
			{
				int cdf1=!overflow?CDF[sym]:0x10000;
				if(sym<255)
					overflow|=cdf1>CDF[sym+1];
				int f1=(sym<255&&!overflow?CDF[sym+1]:0x10000)-cdf1;

				int f2=h2[sym];
				int freq=(f2<<16)/count2;//normalize

				//if(f1||f2)//
				//	printf("");

				freq=f1+(int)(((long long)freq-f1)*alpha>>16);//blend
				freq=((unsigned)(freq*0xFF00)>>16)+1;//guard
				if(freq<0||freq>0xFF01)
					LOG_ERROR("Impossible freq 0x%04X / 0x10000", freq);
				h2[sym]=sum;
				sum+=freq;
				if(sum>0x10000)
					LOG_ERROR("ANS CDF sum %d, freq %d", sum, freq);
			}
			h2[256]=0x10000;
		}
		else
		{
			int mean=128<<16;
			int start=(int)((long long)(0-mean)*conf>>16),
				end=(int)((long long)((256<<16)-mean)*conf>>16);
			int
				estart=error_func_p16(start),
				eend=error_func_p16(end);
			int den=eend-estart;
			int sum=0;
			for(int sym=0;sym<256;++sym)
			{
				int curr=(int)((long long)(( sym   <<16)-mean)*conf>>16),
					next=(int)((long long)(((sym+1)<<16)-mean)*conf>>16);
				int e1=error_func_p16(curr),
					e2=error_func_p16(next);
				int f_static=(int)(((long long)(e2-e1)<<16)/den);
				
				int f_dynamic=h2[sym];
				f_dynamic=(f_dynamic<<16)/count2;//normalize

				int freq=f_static+(int)(((long long)f_dynamic-f_static)*alpha>>16);//blend
				freq=((unsigned)(freq*0xFF00)>>16)+1;//guard

				h2[sym]=sum;
				sum+=freq;
			}
			h2[256]=0x10000;
		}
	}
}
size_t test17_encode(const unsigned char *src, int bw, int bh, int blocksize, int alpha, ArrayHandle *data, int loud)
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

	int totalhsize=256;
	unsigned *hist=(unsigned*)malloc((size_t)totalhsize*sizeof(unsigned));
	unsigned short *CDF=(unsigned short*)malloc((size_t)totalhsize*sizeof(short));
	unsigned *h2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!hist||!CDF||!h2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memset(hist, 0, (size_t)totalhsize*sizeof(unsigned));

	int trivial[2]={0};
	int conf[2]=
	{
		get_conf_p16(b2, bw, 0, 0, bw, 0, bh, 128<<16, trivial)<<2,
		get_conf_p16(b2, bw, 2, 0, bw, 0, bh, 128<<16, trivial+1)<<2,
	};
	for(int kc=1;kc<2;++kc)
	{
		for(int k=0;k<res;++k)
		{
			unsigned char sym=b2[k<<2|kc];
			++hist[sym];
		}
	}
	normalize_histogram(hist, 256, res, CDF);//this is just to pack the histogram, CDF is renormalized again with ramp guard
	//normalize_histogram(hist+256, 256, res, CDF+256);
	//normalize_histogram(hist+512, 256, res, CDF+512);

	DList list;
	dlist_init(&list, 1, 1024, 0);

	int ansbookmarks[3]={0};
	dlist_push_back(&list, 0, 12);
	dlist_push_back(&list, conf, 2);
	dlist_push_back(&list, conf+1, 2);
	dlist_push_back(&list, CDF, totalhsize*sizeof(short));

	int bxcount=(bw+blocksize-1)/blocksize,
		bycount=(bh+blocksize-1)/blocksize;
	for(int kc=0;kc<3;++kc)
	{
		unsigned state=0x10000;
		for(int by=bycount-1;by>=0;--by)
		{
			int ky=by*blocksize;
			for(int bx=bxcount-1;bx>=0;--bx)//for each block
			{
				int kx=bx*blocksize;
				int xsize=0, ysize=0;
				t17_prepblock(b2, kc==1?CDF:0, conf[kc>>1], bw, bh, kc, bx, by, alpha, blocksize, h2, &xsize, &ysize);

#if 0
				printf("kc %d XY %3d %3d \t", kc, bx, by);//
				//if(kc==0&&bx==1&&by==1)//
					print_CDF(h2, b2, bw, bh, kc, kx, kx+xsize, ky, ky+ysize);
#endif

				for(int ky2=ysize-1;ky2>=0;--ky2)
				{
					for(int kx2=xsize-1;kx2>=0;--kx2)//for each pixel
					{
						unsigned char sym=b2[(bw*(ky+ky2)+kx+kx2)<<2|kc];

						int cdf=h2[sym], freq=h2[sym+1]-cdf;

						//if(kc==2&&ky+ky2==0&&kx+kx2==0)//
						//	printf("sym 0x%02X cdf 0x%04X freq 0x%04X\n", sym, cdf, freq);

						if(!freq)
							LOG_ERROR("ZPS");
						
						if(state>=(unsigned)(freq<<16))//renorm
						{
							dlist_push_back(&list, &state, 2);
							state>>=16;
						}
						state=state/freq<<16|(cdf+state%freq);//update
					}
				}
			}
		}
		dlist_push_back(&list, &state, 4);
		ansbookmarks[kc]=(int)list.nobj;
	}
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, ansbookmarks, 12);

	if(loud)
	{
		printf("Overhead %7d\n", (int)(totalhsize*sizeof(short)));
		printf("Red      %7d\n", ansbookmarks[0]-(int)(totalhsize*sizeof(short)));
		printf("Green    %7d\n", ansbookmarks[1]-ansbookmarks[0]);
		printf("Blue     %7d\n", ansbookmarks[2]-ansbookmarks[1]);
	}

	dlist_clear(&list);
	free(b2);
	free(hist);
	free(CDF);
	free(h2);
	return 1;
}
int test17_decode(const unsigned char *data, size_t srclen, int bw, int bh, int blocksize, int alpha, unsigned char *buf)
{
	const int cdflen=256LL*sizeof(short), overhead=12LL+2+2+cdflen;
	int res=bw*bh;
	
	const unsigned char *srcptr, *srcstart, *srcend=data+srclen;
	if(data+overhead>=srcend)
	{
		LOG_ERROR("Invalid file");
		return 0;
	}

	unsigned ansbookmarks[3];
	memcpy(ansbookmarks, data, 12);
	if(ansbookmarks[2]<(unsigned)overhead)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	if(ansbookmarks[2]>srclen)
	{
		LOG_ERROR("Incomplete file");
		return 0;
	}

	unsigned conf[2]={0};
	memcpy(conf, data+12, 2);
	memcpy(conf+1, data+14, 2);
	
	unsigned short *CDF=(unsigned short*)malloc(cdflen);
	unsigned *h2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!CDF||!h2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(CDF, data+16, cdflen);
	
	int bxcount=(bw+blocksize-1)/blocksize,
		bycount=(bh+blocksize-1)/blocksize;
	for(int kc=0;kc<3;++kc)
	{
		unsigned state;
		srcptr=data+ansbookmarks[kc];
		srcstart=kc?data+ansbookmarks[kc-1]:data+overhead;
		srcptr-=4;
		if(srcptr<srcstart)
			LOG_ERROR("ANS buffer overflow");
		memcpy(&state, srcptr, 4);
		
		for(int by=0;by<bycount;++by)
		{
			int ky=by*blocksize;
			for(int bx=0;bx<bxcount;++bx)//for each block
			{
				int kx=bx*blocksize;
				int xsize, ysize;
				t17_prepblock(buf, kc==1?CDF:0, conf[kc>>1], bw, bh, kc, bx, by, alpha, blocksize, h2, &xsize, &ysize);
				for(int ky2=0;ky2<ysize;++ky2)
				{
					for(int kx2=0;kx2<xsize;++kx2)//for each pixel
					{
						unsigned c=(unsigned short)state;
						size_t sym=0;
						int found=binary_search(h2, 257, sizeof(int), threeway_uint32, &c, &sym);
						sym-=!found;//binary_search gives insertion index

						//if(sym!=debug_ptr[(bw*(ky+ky2)+kx+kx2)<<2|kc])//
						//	LOG_ERROR("");

						buf[(bw*(ky+ky2)+kx+kx2)<<2|kc]=(unsigned char)sym;

						unsigned cdf=h2[sym], freq=h2[sym+1]-cdf;

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
	free(CDF);
	free(h2);
	for(int k=0;k<res;++k)//set alpha
		buf[k<<2|3]=0xFF;
	apply_transforms_inv(buf, bw, bh);
	return 1;
}
#endif

#if 0
size_t test18_encode(const unsigned char *src, int bw, int bh, ArrayHandle *data)
{
	int res=bw*bh;
	unsigned char *b2=(unsigned char*)malloc((size_t)res<<2);
	if(!b2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(b2, src, (size_t)res<<2);

	addbuf(b2, bw, bh, 3, 4, 128);//unsigned char -> signed char
	colortransform_ycocgt_fwd((char*)b2, bw, bh);
	ArrayHandle sizes=dwt2d_gensizes(bw, bh, 3, 3, 0);
	unsigned char *temp=(unsigned char*)malloc(MAXVAR(bw, bh));
	DWTSize *psizes=(DWTSize*)sizes->data;
	int nsizes=(int)sizes->count;
	for(int kc=0;kc<3;++kc)
		//dwt2d_haar_fwd ((char*)b2+kc, psizes, 0, nsizes, 4, (char*)temp);
		dwt2d_squeeze_fwd((char*)b2+kc, psizes, 0, nsizes, 4, (char*)temp);
		//dwt2d_cdf53_fwd((char*)b2+kc, psizes, 0, nsizes, 4, (char*)temp);
		//dwt2d_cdf97_fwd((char*)b2+kc, psizes, 0, nsizes, 4, (char*)temp);
	array_free(&sizes);
	free(temp);
	addbuf(b2, bw, bh, 3, 4, 128);//unsigned char -> signed char
	
	const int totalhsize=256*3*4;
	unsigned *hist=(unsigned*)malloc((size_t)totalhsize*sizeof(unsigned));
	unsigned short *CDF=(unsigned short*)malloc((size_t)totalhsize*sizeof(short));
	if(!hist||!CDF)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	
	int px=bw>>1, py=bh>>1;
	memset(hist, 0, (size_t)totalhsize*sizeof(unsigned));
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<bh;++ky)
		{
			for(int kx=0;kx<bw;++kx)
			{
				unsigned char sym=b2[(bw*ky+kx)<<2|kc];
				int quad=(ky>=py)<<1|(kx>=px);
				++hist[kc<<10|quad<<8|sym];
			}
		}
	}
#if 0
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<py;++ky)
		{
			for(int kx=0;kx<px;++kx)
			{
				unsigned char sym=b2[(bw*ky+kx)<<2|kc];
				++hist[kc<<10|sym];
			}
			for(int kx=px;kx<bw;++kx)
			{
				unsigned char sym=b2[(bw*ky+kx)<<2|kc];
				++hist[kc<<10|1<<8|sym];
			}
		}
		for(int ky=py;ky<bh;++ky)
		{
			for(int kx=0;kx<px;++kx)
			{
				unsigned char sym=b2[(bw*ky+kx)<<2|kc];
				++hist[kc<<10|2<<8|sym];
			}
			for(int kx=px;kx<bw;++kx)
			{
				unsigned char sym=b2[(bw*ky+kx)<<2|kc];
				++hist[kc<<10|3<<8|sym];
			}
		}
	}
#endif
	int usizes[]={px*py, (bw-px)*py, px*(bh-py), (bw-px)*(bh-py)};
	normalize_histogram(hist       , 256, usizes[0], CDF       );//512*12 = 6 KB overhead
	normalize_histogram(hist+256   , 256, usizes[1], CDF+256   );
	normalize_histogram(hist+256* 2, 256, usizes[2], CDF+256* 2);
	normalize_histogram(hist+256* 3, 256, usizes[3], CDF+256* 3);
	normalize_histogram(hist+256* 4, 256, usizes[0], CDF+256* 4);
	normalize_histogram(hist+256* 5, 256, usizes[1], CDF+256* 5);
	normalize_histogram(hist+256* 6, 256, usizes[2], CDF+256* 6);
	normalize_histogram(hist+256* 7, 256, usizes[3], CDF+256* 7);
	normalize_histogram(hist+256* 8, 256, usizes[0], CDF+256* 8);
	normalize_histogram(hist+256* 9, 256, usizes[1], CDF+256* 9);
	normalize_histogram(hist+256*10, 256, usizes[2], CDF+256*10);
	normalize_histogram(hist+256*11, 256, usizes[3], CDF+256*11);
	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, 3*sizeof(int));
	dlist_push_back(&list, CDF, totalhsize*sizeof(short));
	int ansbookmarks[3]={0};
	int sizedetails[12]={0};
	for(int kc=0;kc<3;++kc)
	{
		unsigned state=0x10000;
		for(int ky=bh-1;ky>=0;--ky)
		{
			for(int kx=bw-1;kx>=0;--kx)
			{
				unsigned char sym=b2[(bw*ky+kx)<<2|kc];
				int quad=(ky>=py)<<1|(kx>=px);
				unsigned short *CDFptr=CDF+(kc<<10|quad<<8);

				unsigned cdf=CDFptr[sym], freq=(sym<255&&cdf<CDFptr[sym+1]?CDFptr[sym+1]:0x10000)-cdf;

				if(state>=(unsigned)(freq<<16))//renorm
				{
					sizedetails[kc<<2|quad]+=2;//

					dlist_push_back(&list, &state, 2);
					state>>=16;
				}
				state=state/freq<<16|(cdf+state%freq);//update
			}
		}
		++sizedetails[kc<<2  ];//
		++sizedetails[kc<<2|1];//
		++sizedetails[kc<<2|2];//
		++sizedetails[kc<<2|3];//

		dlist_push_back(&list, &state, 4);
		ansbookmarks[kc]=(int)list.nobj;
	}
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, ansbookmarks, 12);
	
	printf("Overhead %7d\n", 6144);
	printf("red\n");
	printf(" %7d %7d %9lf %9lf\n", sizedetails[ 0], sizedetails[ 1], (double)px*    py /sizedetails[ 0], (double)(bw-px)*    py /sizedetails[ 1]);
	printf(" %7d %7d %9lf %9lf\n", sizedetails[ 2], sizedetails[ 3], (double)px*(bh-py)/sizedetails[ 2], (double)(bw-px)*(bh-py)/sizedetails[ 3]);
	printf("green\n");
	printf(" %7d %7d %9lf %9lf\n", sizedetails[ 4], sizedetails[ 5], (double)px*    py /sizedetails[ 4], (double)(bw-px)*    py /sizedetails[ 5]);
	printf(" %7d %7d %9lf %9lf\n", sizedetails[ 6], sizedetails[ 7], (double)px*(bh-py)/sizedetails[ 6], (double)(bw-px)*(bh-py)/sizedetails[ 7]);
	printf("blue\n");
	printf(" %7d %7d %9lf %9lf\n", sizedetails[ 8], sizedetails[ 9], (double)px*    py /sizedetails[ 8], (double)(bw-px)*    py /sizedetails[ 9]);
	printf(" %7d %7d %9lf %9lf\n", sizedetails[10], sizedetails[11], (double)px*(bh-py)/sizedetails[10], (double)(bw-px)*(bh-py)/sizedetails[11]);
	printf("TOTAL %7lld %9lf\n", list.nobj, (double)bw*bh*3/list.nobj);
#if 0
	printf("Overhead      %7d\n", 6144);
	printf("R topleft     %7d\n", sizedetails[ 0]);
	printf("R topright    %7d\n", sizedetails[ 1]);
	printf("R bottomleft  %7d\n", sizedetails[ 2]);
	printf("R bottomright %7d\n", sizedetails[ 3]);
	printf("G topleft     %7d\n", sizedetails[ 4]);
	printf("G topright    %7d\n", sizedetails[ 5]);
	printf("G bottomleft  %7d\n", sizedetails[ 6]);
	printf("G bottomright %7d\n", sizedetails[ 7]);
	printf("B topleft     %7d\n", sizedetails[ 8]);
	printf("B topright    %7d\n", sizedetails[ 9]);
	printf("B bottomleft  %7d\n", sizedetails[10]);
	printf("B bottomright %7d\n", sizedetails[11]);
	printf("TOTAL         %7lld\n", list.nobj);
#endif

	dlist_clear(&list);
	free(b2);
	free(hist);
	free(CDF);
	return 1;
}
#endif


#ifdef DEBUG_ANS
SList states={0};
int debug_channel=0;
void debug_enc_update(unsigned state, unsigned cdf, unsigned freq, int kx, int ky, int kq, int kc, unsigned char sym)
{
	if(kc==debug_channel)
	{
		if(!states.count)
			slist_init(&states, sizeof(DebugANSInfo), 0);

		state=state/freq<<16|(cdf+state%freq);//update

		DebugANSInfo info={state, cdf, freq, (unsigned)states.count, kx, ky, kq, kc, sym};
		STACK_PUSH(&states, &info);
	}
}
void debug_dec_update(unsigned state, unsigned cdf, unsigned freq, int kx, int ky, int kq, int kc, unsigned char sym)
{
	if(kc==debug_channel)
	{
		if(!states.count)
			LOG_ERROR("Nothing to decode");
		DebugANSInfo *i0=(DebugANSInfo*)STACK_TOP(&states), info;
		memcpy(&info, i0, sizeof(info));

		if(info.state!=state||info.cdf!=cdf||info.freq!=freq||kx!=info.kx||ky!=info.ky||kq!=info.kq||kc!=info.kc||info.sym!=sym)
			LOG_ERROR("Decode error  (%d decodes remaining)", info.id);

		state=freq*(state>>16)+(unsigned short)state-cdf;//update

		STACK_POP(&states);
	}
}
#endif

#if 0
void t19_calchist(const unsigned char *b2, int bw, int kc, int x1, int x2, int y1, int y2, int xslope, int yslope, int constant, int maxinc, unsigned *CDF2, int *hcount)
{
	//int constant=maxinc+(y2-y1)*yslope+(x2-x1)*xslope;
	for(int ky2=y1;ky2<y2;++ky2)
	{
		int yterm=(ky2-y1)*yslope;
		for(int kx2=x1;kx2<x2;++kx2)//for each pixel
		{
			unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];
			int inc=(kx2-x1)*xslope+yterm+constant;
			if(inc<=0||inc>maxinc)
				LOG_ERROR("Wrong increment");
			CDF2[sym]+=inc;
			*hcount+=inc;
			if(*hcount<inc)
				LOG_ERROR("Histogram OVERFLOW");
		}
	}
}
void t19_prepblock(const unsigned short *CDF, const unsigned char *b2, int kc, DWTSize *sizes, int nsizes, int it, int *bandbounds, int bx, int by, int alpha, int blocksize, int margin, unsigned *CDF2, int *blockbounds)
{
	int hcount=0;
	if(b2&&bandbounds&&blockbounds)
	{
		blockbounds[0]=bandbounds[0]+bx*blocksize;
		blockbounds[1]=blockbounds[0]+blocksize;
		blockbounds[2]=bandbounds[2]+by*blocksize;
		blockbounds[3]=blockbounds[2]+blocksize;
		if(blockbounds[1]>bandbounds[1])
			blockbounds[1]=bandbounds[1];
		if(blockbounds[3]>bandbounds[3])
			blockbounds[3]=bandbounds[3];

		int left=blockbounds[0]-margin, right=blockbounds[1]+margin, top=blockbounds[2]-margin;
		if(left<bandbounds[0])
			left=bandbounds[0];
		if(right>bandbounds[1])
			right=bandbounds[1];
		if(top<bandbounds[2])
			top=bandbounds[2];

		int marginxleft=blockbounds[0]-left,
			marginxright=right-blockbounds[1],
			marginy=blockbounds[2]-top,
			maxmarginx;
		if(marginxleft<0)
			marginxleft=0;
		if(marginxright<0)
			marginxright=0;
		if(marginy<0)
			marginy=0;
		maxmarginx=MAXVAR(marginxleft, marginxright);
		int maxinc=maxmarginx+marginy+1;
		memset(CDF2, 0, 256*sizeof(unsigned));
		if(maxinc>0)
		{
			int dwtblockbounds[4];
			memcpy(dwtblockbounds, blockbounds, sizeof(dwtblockbounds));
			for(;;)
			{
				dwtblockbounds[0]>>=1;
				dwtblockbounds[1]>>=1;
				dwtblockbounds[2]>>=1;
				dwtblockbounds[3]>>=1;
				if(dwtblockbounds[0]>=dwtblockbounds[1]||dwtblockbounds[2]>=dwtblockbounds[3])
					break;
				t19_calchist(b2, sizes->w, kc, blockbounds[0], blockbounds[1], top, blockbounds[2], 0, 0, maxinc, maxinc, CDF2, &hcount);//DWT block(s)
			}

			//LPF
#if 0
			//for(int k=0;k<255;++k)
			//	CDF2[k]+=CDF[k+1]>>1;
			//for(int k=255;k>0;--k)
			//	CDF2[k]+=CDF[k-1]>>1;
			for(int k=0;k<255;++k)
				CDF2[k]=(CDF2[k]+CDF[k+1])>>1;
			for(int k=255;k>0;--k)
				CDF2[k]=(CDF2[k]+CDF[k-1])>>1;
#endif

			if(left<blockbounds[0])
			{
				t19_calchist(b2, sizes->w, kc, left, blockbounds[0], blockbounds[2], blockbounds[3], 1, 0, maxinc-(blockbounds[0]-left), maxinc, CDF2, &hcount);//left
				if(top<blockbounds[2])
					t19_calchist(b2, sizes->w, kc, left, blockbounds[0], top, blockbounds[2], 1, 1, maxinc-((blockbounds[0]-left)+(blockbounds[2]-top)-1), maxinc, CDF2, &hcount);//topleft
			}
			if(top<blockbounds[2])
			{
				t19_calchist(b2, sizes->w, kc, blockbounds[0], blockbounds[1], top, blockbounds[2], 0, 1, maxinc-(blockbounds[2]-top), maxinc, CDF2, &hcount);//top
				if(blockbounds[1]<right)
					t19_calchist(b2, sizes->w, kc, blockbounds[1], right, top, blockbounds[2], -1, 1, maxinc-(blockbounds[2]-top), maxinc, CDF2, &hcount);//topright
			}
#if 0
			hcount=0;
			for(int k=0;k<256;++k)
				hcount+=CDF[k];
#endif
		}
	}

	int overflow=0;
	if(hcount)
	{
		int sum=0;
		for(int sym=0;sym<256;++sym)
		{
			int cdf1=!overflow?CDF[sym]:0x10000;
			if(sym<255)
				overflow|=cdf1>CDF[sym+1];
			int f1=(sym<255&&!overflow?CDF[sym+1]:0x10000)-cdf1;

			int f2=CDF2[sym];
			int freq=(int)(((long long)f2<<16)/hcount);//normalize

			//if(f1||f2)//
			//	printf("");

			freq=f1+(int)(((long long)freq-f1)*alpha>>16);//blend
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
size_t test19_encode(const unsigned char *src, int bw, int bh, int alpha, int blocksize, int margin, ArrayHandle *data, int loud)
{
	int res=bw*bh, totalhsize=768;
	unsigned char *b2=(unsigned char*)malloc((size_t)res<<2);
	unsigned *hist=(unsigned*)malloc((size_t)totalhsize*sizeof(unsigned));
	unsigned short *CDF=(unsigned short*)malloc((size_t)totalhsize*sizeof(short));
	unsigned *CDF2=(unsigned*)malloc(257LL*sizeof(unsigned));
	if(!b2||!hist||!CDF||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(b2, src, (size_t)res<<2);

	ArrayHandle sizes=dwt2d_gensizes(bw, bh, 3, 3, 0);
	DWTSize *psizes=(DWTSize*)sizes->data;
	int nsizes=(int)sizes->count;

#if 1
	addbuf(b2, bw, bh, 3, 4, 128);//unsigned char -> signed char
	colortransform_ycocgt_fwd((char*)b2, bw, bh);
	pred_grad_fwd((char*)b2, bw, bh, 3, 4);

	unsigned char *temp=(unsigned char*)malloc(MAXVAR(bw, bh));
	for(int kc=0;kc<3;++kc)
		dwt2d_lazy_fwd((char*)b2+kc, psizes, 0, nsizes, 4, (char*)temp);
		//dwt2d_haar_fwd   ((char*)b2+kc, psizes, 0, nsizes, 4, (char*)temp);
		//dwt2d_squeeze_fwd((char*)b2+kc, psizes, 0, nsizes, 4, (char*)temp);
		//dwt2d_cdf53_fwd  ((char*)b2+kc, psizes, 0, nsizes, 4, (char*)temp);
		//dwt2d_cdf97_fwd  ((char*)b2+kc, psizes, 0, nsizes, 4, (char*)temp);
	addbuf(b2, bw, bh, 3, 4, 128);//back to unsigned char
	free(temp);
#endif
	
	memset(hist, 0, (size_t)totalhsize*sizeof(unsigned));
	for(int kc=0;kc<3;++kc)
	{
		for(int k=0;k<res;++k)
		{
			unsigned char sym=b2[k<<2|kc];
			++hist[kc<<8|sym];
		}
	}
	normalize_histogram(hist, 256, res, CDF);//this is just to pack the histogram, CDF is renormalized again with ramp guard
	normalize_histogram(hist+256, 256, res, CDF+256);
	normalize_histogram(hist+512, 256, res, CDF+512);
	free(hist);

	DList list;
	int ansbookmarks[3]={0};
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, 12);//ansbookmarks
	dlist_push_back(&list, CDF, 768*sizeof(short));

	for(int kc=0;kc<3;++kc)
	{
		unsigned state=0x10000;
		for(int ks=0;ks<nsizes-1;++ks)//for each DWT size
		{
			int w2=psizes[ks].w, h2=psizes[ks].h, px=w2>>1, py=h2>>1;
			for(int kq=3;kq>=1;--kq)//for each high quadrant
			{
				int bounds[]={px&-(kq&1), kq&1?w2:px, py&-(kq>>1), kq>>1?h2:py};//x1, x2, y1, y2
				int blockbounds[4];
				//int x1=px&-(kq&1), x2=x1?w2:px, y1=py&-(kq>>1), y2=y1?h2:py;

				int bxcount=(bounds[1]-bounds[0]+blocksize-1)/blocksize,
					bycount=(bounds[3]-bounds[2]+blocksize-1)/blocksize;
				for(int by=bycount-1;by>=0;--by)
				{
					for(int bx=bxcount-1;bx>=0;--bx)//for each block
					{
						t19_prepblock(CDF+(kc<<8), b2, kc, psizes, nsizes, ks, bounds, bx, by, alpha, blocksize, margin, CDF2, blockbounds);

						//encode block
						for(int ky=blockbounds[3]-1;ky>=blockbounds[2];--ky)
						{
							for(int kx=blockbounds[1]-1;kx>=blockbounds[0];--kx)//for each pixel
							{
								unsigned char sym=b2[(bw*ky+kx)<<2|kc];

								int cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;

								//unsigned s0=state;

								if(!freq)
									LOG_ERROR("ZPS");
						
								if(state>=(unsigned)(freq<<16))//renorm
								{
									dlist_push_back(&list, &state, 2);
									state>>=16;
								}
								debug_enc_update(state, cdf, freq, kx, ky, kq, kc, sym);
								state=state/freq<<16|(cdf+state%freq);//update
								
								//if(kc==0&&kx==3&&ky==0)//
								//	printf("\nsym 0x%02X cdf 0x%04X freq 0x%04X, state 0x%08X -> 0x%08X\n", sym, cdf, freq, s0, state);
							}
						}
					}
				}
			}
		}

		//encode DC block
		int w0=psizes[nsizes-1].w, h0=psizes[nsizes-1].h;
		t19_prepblock(CDF+(kc<<8), 0, kc, psizes, nsizes, nsizes-1, 0, 0, 0, alpha, blocksize, margin, CDF2, 0);
		for(int ky=h0-1;ky>=0;--ky)
		{
			for(int kx=w0-1;kx>=0;--kx)
			{
				unsigned char sym=b2[(bw*ky+kx)<<2|kc];
				int cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;
				//int cdf=CDF[kc<<8|sym], freq=(sym<255&&cdf<CDF[kc<<8|(sym+1)]?CDF[kc<<8|(sym+1)]:0x10000)-cdf;

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
		dlist_push_back(&list, &state, 4);
		ansbookmarks[kc]=(int)list.nobj;
	}
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, ansbookmarks, 12);
	if(loud)
	{
		int osize=(int)(totalhsize*sizeof(short)), rsize=ansbookmarks[0]-osize, bsize=ansbookmarks[1]-ansbookmarks[0], gsize=ansbookmarks[2]-ansbookmarks[1], s0=bw*bh;
		printf("alpha 0x%04X block %d margin %d\n", alpha, blocksize, margin);
		printf("Overhead %7d\n", osize);
		printf("Red      %7d  %lf\n", rsize, (double)s0/rsize);
		printf("Green    %7d  %lf\n", bsize, (double)s0/bsize);
		printf("Blue     %7d  %lf\n", gsize, (double)s0/gsize);
		printf("Uncompr. %7d\n", s0);
	}
	dlist_clear(&list);
	array_free(&sizes);
	free(b2);
	free(CDF);
	free(CDF2);
	return 1;
}
int test19_decode(const unsigned char *data, size_t srclen, int bw, int bh, int alpha, int blocksize, int margin, unsigned char *buf)
{
	const int cdflen=768LL*sizeof(short), overhead=12LL+cdflen;
	int res=bw*bh;
	
	const unsigned char *srcptr, *srcstart, *srcend=data+srclen;
	if(data+overhead>=srcend)
	{
		LOG_ERROR("Invalid file");
		return 0;
	}

	unsigned ansbookmarks[3];
	memcpy(ansbookmarks, data, 12);
	if(ansbookmarks[2]<(unsigned)overhead)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	if(ansbookmarks[2]>srclen)
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
	memcpy(CDF, data+12, cdflen);
	
	ArrayHandle sizes=dwt2d_gensizes(bw, bh, 3, 3, 0);
	DWTSize *psizes=(DWTSize*)sizes->data;
	int nsizes=(int)sizes->count;

	for(int kc=0;kc<3;++kc)
	{
		unsigned state;
		srcptr=data+ansbookmarks[kc];
		srcstart=data+(kc?ansbookmarks[kc-1]:overhead);
		srcptr-=4;
		if(srcptr<srcstart)
			LOG_ERROR("ANS buffer overflow");
		memcpy(&state, srcptr, 4);
		
		t19_prepblock(CDF+(kc<<8), 0, kc, psizes, nsizes, nsizes-1, 0, 0, 0, alpha, blocksize, margin, CDF2, 0);
		for(int ky=0;ky<psizes[nsizes-1].h;++ky)
		{
			for(int kx=0;kx<psizes[nsizes-1].w;++kx)
			{
				unsigned c=(unsigned short)state;
				size_t sym=0;
				int found=binary_search(CDF2, 257, sizeof(int), threeway_uint32, &c, &sym);
				sym-=!found;//binary_search gives insertion index

				buf[(bw*ky+kx)<<2|kc]=(unsigned char)sym;
				int cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;
				//int cdf=CDF[kc<<8|sym], freq=(sym<255&&cdf<CDF[kc<<8|(sym+1)]?CDF[kc<<8|(sym+1)]:0x10000)-cdf;

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
		
		for(int ks=nsizes-2;ks>=0;--ks)//for each DWT size		smallest to largest
		{
			int w2=psizes[ks].w, h2=psizes[ks].h, px=w2>>1, py=h2>>1;
			for(int kq=1;kq<4;++kq)//for each high quadrant
			{
				int bounds[]={px&-(kq&1), kq&1?w2:px, py&-(kq>>1), kq>>1?h2:py};//x1, x2, y1, y2
				int bxcount=(bounds[1]-bounds[0]+blocksize-1)/blocksize,
					bycount=(bounds[3]-bounds[2]+blocksize-1)/blocksize;
				for(int by=0;by<bycount;++by)
				{
					for(int bx=0;bx<bxcount;++bx)//for each block
					{
						//if(kc==0&&bx==0&&by==0)
						//	kc=0;
						int blockbounds[4];
						t19_prepblock(CDF+(kc<<8), buf, kc, psizes, nsizes, ks, bounds, bx, by, alpha, blocksize, margin, CDF2, blockbounds);
						for(int ky=blockbounds[2];ky<blockbounds[3];++ky)
						{
							for(int kx=blockbounds[0];kx<blockbounds[1];++kx)//for each pixel
							{
								unsigned c=(unsigned short)state;
								size_t sym=0;
								int found=binary_search(CDF2, 257, sizeof(int), threeway_uint32, &c, &sym);
								sym-=!found;//binary_search gives insertion index
								
								//if(debug_ptr&&sym!=debug_ptr[(bw*ky+kx)<<2|kc])//
								//	LOG_ERROR("");
								buf[(bw*ky+kx)<<2|kc]=(unsigned char)sym;

								unsigned cdf=CDF2[sym], freq=CDF2[sym+1]-cdf;
								
								debug_dec_update(state, cdf, freq, kx, ky, kq, kc, sym);
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
		}
	}
	free(CDF);
	free(CDF2);
	for(int k=0;k<res;++k)//set alpha
		buf[k<<2|3]=0xFF;
#if 1
	unsigned char *temp=(unsigned char*)malloc(MAXVAR(bw, bh));
	addbuf(buf, bw, bh, 3, 4, 128);//unsigned char -> signed char
	for(int kc=0;kc<3;++kc)
		dwt2d_lazy_inv((char*)buf+kc, psizes, 0, nsizes, 4, (char*)temp);
		//dwt2d_haar_inv   ((char*)buf+kc, psizes, 0, nsizes, 4, (char*)temp);
		//dwt2d_squeeze_inv((char*)buf+kc, psizes, 0, nsizes, 4, (char*)temp);
		//dwt2d_cdf53_inv  ((char*)buf+kc, psizes, 0, nsizes, 4, (char*)temp);
		//dwt2d_cdf97_inv  ((char*)buf+kc, psizes, 0, nsizes, 4, (char*)temp);
	
	pred_grad_inv((char*)buf, bw, bh, 3, 4);
	colortransform_ycocgt_inv((char*)buf, bw, bh);
	addbuf(buf, bw, bh, 3, 4, 128);//back to unsigned char
	free(temp);
#endif
	return 1;
}
#endif

#if 0
void t20_prepblock(const unsigned char *b2, unsigned short *CDF, int bw, int bh, int kc, int bx, int by, int alpha, int blocksize, int margin, unsigned *CDF2, int *xend, int *yend)
{
	int kx=bx*blocksize, ky=by*blocksize;

	*yend=ky+blocksize<=bh?ky+blocksize:bh;
	*xend=kx+blocksize<=bw?kx+blocksize:bw;
				
	int overflow=0;//CDF overflow can happen only once
	if(!bx&&!by||!alpha)//first block has zero alpha
	{
		//if(bx==1&&!by)
		//	bx=1;
		for(int sym=0;sym<256;++sym)
		{
			if(overflow)
				CDF2[sym]=0x10000;
				//CDF2[sym]=0xFF00|sym;
			else
			{
				int cdf0=CDF[sym], cdf=cdf0;
				//cdf=((unsigned)(cdf*0xFF00)>>16)+sym;
				CDF2[sym]=cdf;
				if(sym<255&&cdf0>CDF[sym+1])
					overflow=1;
			}
		}
		CDF2[256]=0x10000;
	}
	else
	{
		memset(CDF2, 0, 256*sizeof(unsigned));

		int count2=0;

		int left=kx-margin;
		if(left<0)
			left=0;
		int right=kx+blocksize+margin;
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
					unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];
					int dist=kx-kx2;
					if(dist<0||dist>margin)
						LOG_ERROR("Wrong distance");
					
					int inc=(margin<<1|1)-dist;
					//int inc=(kx2>>1)+1;
					//int inc=(kx2*3>>2)+1;
					//int inc=kx2<<1|1;
					//int inc=1;
					//int inc=kx2*kx2+1;
					//int inc=0x10000/(kx2+1);

					if(!inc)
						LOG_ERROR("Zero inc");

					CDF2[sym]+=inc;
					count2+=inc;
					//if(count2<inc)
					//	LOG_ERROR("OVERFLOW");
				}
			}
			//count2+=(blocksize*(blocksize+1)>>1)**ysize;
			//count2+=blocksize**ysize;
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
					//int inc=(ky2>>1)+1;
					//int inc=(ky2*3>>2)+1;
					//int inc=ky2<<1|1;
					//int inc=1;
					//int inc=ky2*ky2+1;
					//int inc=0x10000/(ky2+1);
					
					if(!inc)
						LOG_ERROR("Zero inc");

					CDF2[sym]+=inc;
					count2+=inc;
					if(count2<inc)
						LOG_ERROR("OVERFLOW");
				}
			}
			//count2+=(blocksize*(blocksize+1)>>1)**xsize;
			//count2+=blocksize**xsize;
		}
		if(left<kx&&top<ky)//if topleft block is available
		{
			for(int ky2=top;ky2<ky;++ky2)
			{
				for(int kx2=left;kx2<kx;++kx2)//for each pixel
				{
					unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];
					int dist=kx-kx2+ky-ky2;
					//int dist=MAXVAR(kx-kx2, ky-ky2);
					
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
		if(right>kx+blocksize&&top<ky)//if topright block is available
		{
			for(int ky2=top;ky2<ky;++ky2)
			{
				for(int kx2=kx+blocksize;kx2<right;++kx2)//for each pixel
				{
					unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];
					int dist=kx2-(kx+blocksize)+ky-ky2;
					//int dist=MAXVAR(kx2-(kx+blocksize), ky-ky2);

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

		int sum=0;
		for(int sym=0;sym<256;++sym)
		{
			int cdf1=!overflow?CDF[sym]:0x10000;
			if(sym<255)
				overflow|=cdf1>CDF[sym+1];
			int f1=(sym<255&&!overflow?CDF[sym+1]:0x10000)-cdf1;

			int f2=CDF2[sym];
			int freq=(int)(((long long)f2<<16)/count2);//normalize

			//if(f1||f2)//
			//	printf("");

			freq=f1+(int)(((long long)freq-f1)*alpha>>16);//blend
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
}
size_t test20_encode(const unsigned char *src, int bw, int bh, int blocksize, int margin, int alpha, int blockcount, ArrayHandle *data, int loud)
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

	int macroblocksize=blocksize*blockcount,
		macroblockcx=(bw+macroblocksize-1)/macroblocksize,
		macroblockcy=(bh+macroblocksize-1)/macroblocksize,
		macroblockcount=macroblockcx*macroblockcy;
	unsigned *hist=(unsigned*)malloc(257LL*sizeof(int));
	unsigned short *CDF=(unsigned short*)malloc(macroblockcount*768LL*sizeof(short));
	if(!hist||!CDF)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	for(int kc=0;kc<3;++kc)
	{
		for(int mby=0;mby<macroblockcy;++mby)
		{
			int yend=macroblocksize*(mby+1);
			if(yend>bh)
				yend=bh;
			int ystart=macroblocksize*mby, dy=yend-ystart;
			for(int mbx=0;mbx<macroblockcx;++mbx)
			{
				int xend=macroblocksize*(mbx+1);
				if(xend>bw)
					xend=bw;
				int xstart=macroblocksize*mbx, dx=xend-xstart;

				memset(hist, 0, 256LL*sizeof(int));
				for(int ky=ystart;ky<yend;++ky)
				{
					for(int kx=xstart;kx<xend;++kx)
					{
						unsigned char sym=b2[(bw*ky+kx)<<2|kc];
						++hist[sym];
					}
				}
				unsigned short *CDFk=CDF+((3*(macroblockcx*mby+mbx)+kc)<<8);
				int sum=0, count=dy*dx;
				for(int sym=0;sym<256;++sym)//normalize & accumulate histogram to fit in an unsigned short & quickly fetch CDF and freq
				{
					int freq=hist[sym];
					//freq=(int)(((long long)freq<<16)/count);
					freq=(int)(((long long)freq*0xFFFF)/count);
					//freq=(int)(((long long)freq*0xFF00)/count+1);
					CDFk[sym]=sum;
					sum+=freq;
					if(sum>0x10000)
						LOG_ERROR("CDF Overflow");
				}
				//if(kc==0&&mby==255&&mbx==356)
				//	mbx=356;
				//normalize_histogram(hist, 256, dx*dy, CDF+((3*(macroblockcx*mby+mbx)+kc)<<8));
			}
		}
	}

	DList list;
	dlist_init(&list, 1, 1024, 0);

	int ansbookmarks[3]={0};
	dlist_push_back(&list, 0, 12);
	dlist_push_back(&list, CDF, macroblockcount*768LL*sizeof(short));

	int bxcount=(bw+blocksize-1)/blocksize,
		bycount=(bh+blocksize-1)/blocksize;
	for(int kc=0;kc<3;++kc)
	{
		unsigned state=0x10000;
		for(int by=bycount-1;by>=0;--by)
		{
			int mby=by/blockcount;
			int ky=by*blocksize;
			for(int bx=bxcount-1;bx>=0;--bx)//for each block
			{
				//if(kc==0&&bx==0&&by==0)
				//	kc=0;
				
				int mbx=bx/blockcount;
				int kx=bx*blocksize;
				int xend=0, yend=0;
				t20_prepblock(b2, CDF+((3*(macroblockcx*mby+mbx)+kc)<<8), bw, bh, kc, bx, by, alpha, blocksize, margin, hist, &xend, &yend);

				//if(kc==0&&bx==0&&by==0)
				//	print_CDF(hist, b2, bw, bh, kc, kx, xend, ky, yend);

				//encode block
				for(int ky2=yend-1;ky2>=ky;--ky2)
				{
					for(int kx2=xend-1;kx2>=kx;--kx2)//for each pixel
					{
						unsigned char sym=b2[(bw*ky2+kx2)<<2|kc];

						int cdf=hist[sym], freq=hist[sym+1]-cdf;

						//if(kc==0&&ky2==0&&kx2==15)//
						//	printf("\n\nsym 0x%02X cdf 0x%04X freq 0x%04X\n\n", sym, cdf, freq);

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
		}
		dlist_push_back(&list, &state, 4);
		ansbookmarks[kc]=(int)list.nobj;
	}
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, ansbookmarks, 12);
	
	if(loud)
	{
		int osize=(int)(macroblockcount*768LL*sizeof(short)), rsize=ansbookmarks[0]-osize, bsize=ansbookmarks[1]-ansbookmarks[0], gsize=ansbookmarks[2]-ansbookmarks[1], s0=bw*bh;
		printf("alpha 0x%04X block %d margin %d\n", alpha, blocksize, margin);
		printf("Overhead %7d\n", osize);
		printf("Red      %7d  %lf\n", rsize, (double)s0/rsize);
		printf("Green    %7d  %lf\n", bsize, (double)s0/bsize);
		printf("Blue     %7d  %lf\n", gsize, (double)s0/gsize);
		printf("Uncompr. %7d\n", s0);

		//printf("alpha 0x%04X block %d margin %d\n", alpha, blocksize, margin);
		//printf("Overhead %7d\n", (int)(macroblockcount*768LL*sizeof(short)));
		//printf("Red      %7d\n", ansbookmarks[0]-(int)(macroblockcount*768LL*sizeof(short)));
		//printf("Green    %7d\n", ansbookmarks[1]-ansbookmarks[0]);
		//printf("Blue     %7d\n", ansbookmarks[2]-ansbookmarks[1]);
	}

	dlist_clear(&list);
	free(b2);
	free(hist);
	free(CDF);
	return 1;
}
int test20_decode(const unsigned char *data, size_t srclen, int bw, int bh, int blocksize, int margin, int alpha, int blockcount, unsigned char *buf)
{
	int macroblocksize=blocksize*blockcount,
		macroblockcx=(bw+macroblocksize-1)/macroblocksize,
		macroblockcy=(bh+macroblocksize-1)/macroblocksize,
		macroblockcount=macroblockcx*macroblockcy;
	int cdflen=macroblockcount*768*(int)sizeof(short), overhead=12+cdflen;
	int res=bw*bh;
	
	const unsigned char *srcptr, *srcstart, *srcend=data+srclen;
	if(data+overhead>=srcend)
	{
		LOG_ERROR("Invalid file");
		return 0;
	}

	unsigned ansbookmarks[3];
	memcpy(ansbookmarks, data, 12);
	if(ansbookmarks[2]<(unsigned)overhead)
	{
		LOG_ERROR("Corrupt file");
		return 0;
	}
	if(ansbookmarks[2]>srclen)
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
	memcpy(CDF, data+12, cdflen);
	
	int bxcount=(bw+blocksize-1)/blocksize,
		bycount=(bh+blocksize-1)/blocksize;
	for(int kc=0;kc<3;++kc)
	{
		unsigned state;
		srcptr=data+ansbookmarks[kc];
		srcstart=kc?data+ansbookmarks[kc-1]:data+overhead;
		srcptr-=4;
		if(srcptr<srcstart)
			LOG_ERROR("ANS buffer overflow");
		memcpy(&state, srcptr, 4);
		
		for(int by=0;by<bycount;++by)
		{
			int mby=by/blockcount;
			int ky=by*blocksize;
			for(int bx=0;bx<bxcount;++bx)//for each block
			{
				//if(kc==0&&bx==0&&by==0)
				//	kc=0;
				
				int mbx=bx/blockcount;
				int kx=bx*blocksize;
				int xend=0, yend=0;
				t20_prepblock(buf, CDF+((3*(macroblockcx*mby+mbx)+kc)<<8), bw, bh, kc, bx, by, alpha, blocksize, margin, CDF2, &xend, &yend);
				for(int ky2=ky;ky2<yend;++ky2)
				{
					for(int kx2=kx;kx2<xend;++kx2)//for each pixel
					{
						unsigned c=(unsigned short)state;
						int sym=0;
						int L=0, R=256, found=0, c0;
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
							sym=L+(L<257&&CDF2[L]<c)-1;
						else
							for(c0=CDF2[sym];sym<255&&CDF2[sym+1]==c0;++sym);

						//int found=binary_search(CDF2, 257, sizeof(int), threeway_uint32, &c, &sym);
						//if(found)
						//	for(int cdf=CDF2[sym];sym<255&&CDF2[sym+1]==cdf;++sym);
						//else//binary_search gives insertion index
						//	--sym;

						//if(kc==0&&ky2==14&&kx2==14)
						//	kx2=14;
						//if(sym!=debug_ptr[(bw*ky2+kx2)<<2|kc])//
						//	LOG_ERROR("");

						buf[(bw*ky2+kx2)<<2|kc]=(unsigned char)sym;

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
	free(CDF);
	free(CDF2);
	for(int k=0;k<res;++k)//set alpha
		buf[k<<2|3]=0xFF;
	apply_transforms_inv(buf, bw, bh);
	return 1;
}
#endif

#if 0
void e10_iter(unsigned char *buf, int iw, int ih, int kc, long long *hist);
int  e10_cmp_overhead(const void *p1, const void *p2)
{
	short idx1=*(const short*)p1, idx2=*(const short*)p2;
	return (idx1>idx2)-(idx1<idx2);
}
void e10_getminiCDF(int *nb, unsigned short *overhead, int overhead_count, int *miniCDF)
{
	int permutation=0;
	char temp;

#define SORT_STEP(A, B)\
	if(nb[1+A]<nb[1+B])\
		permutation+=0;\
	else if(nb[1+A]>nb[1+B])\
		permutation+=1, temp=nb[1+A], nb[1+A]=nb[1+B], nb[1+B]=temp;\
	else\
		permutation+=2;

	SORT_STEP(0, 1);
	permutation*=3;
	SORT_STEP(0, 2);
	permutation*=3;
	SORT_STEP(0, 3);

	permutation*=3;
	SORT_STEP(1, 2);
	permutation*=3;
	SORT_STEP(1, 3);

	permutation*=3;
	SORT_STEP(2, 3);
#undef  SORT_STEP

	size_t idx;
	int found=binary_search(overhead, overhead_count, 10*sizeof(short), e10_cmp_overhead, &permutation, &idx);
	if(!found)
		LOG_ERROR("Permutation %d not found in overhead", permutation);

	unsigned short *ok=overhead+10*idx+1;
	
	miniCDF[0]=0;                //case 0: CDF[1]-CDF[0]: -129<x<nb[1]
	miniCDF[1]=ok[1]+nb[1]  +128;//case 1: CDF[2]-CDF[1]: x==nb[1]
	miniCDF[2]=ok[2]+nb[1]+1+128;//case 2: CDF[3]-CDF[2]: nb[1]<x<nb[2]
	miniCDF[3]=ok[3]+nb[2]  +128;//case 3: CDF[4]-CDF[3]: x==nb[2]
	miniCDF[4]=ok[4]+nb[2]+1+128;//case 4: CDF[5]-CDF[4]: nb[2]<x<nb[3]
	miniCDF[5]=ok[5]+nb[3]  +128;//case 5: CDF[6]-CDF[5]: x==nb[3]
	miniCDF[6]=ok[6]+nb[3]+1+128;//case 6: CDF[7]-CDF[6]: nb[3]<x<nb[4]
	miniCDF[7]=ok[7]+nb[4]  +128;//case 7: CDF[8]-CDF[7]: x==nb[4]
	miniCDF[8]=ok[8]+nb[4]+1+128;//case 8: CDF[9]-CDF[8]: nb[4]<x<128
	miniCDF[9]=0x10000;

	int sum=0x10000;
	if(nb[4]==127)//case 8 is off
		miniCDF[8]=sum;
	else
		sum=miniCDF[8];

	if(nb[3]==nb[4])//cases 7 & 6 are off
		miniCDF[6]=miniCDF[7]=sum;
	else if(nb[3]+1==nb[4])//case 6 is off
		sum=miniCDF[6]=miniCDF[7];
	else
		sum=miniCDF[6];

	if(nb[2]==nb[3])//cases 5 & 4 are off
		miniCDF[4]=miniCDF[5]=sum;
	else if(nb[2]+1==nb[3])//case 4 is off
		sum=miniCDF[4]=miniCDF[5];
	else
		sum=miniCDF[4];

	if(nb[1]==nb[2])//cases 3 & 2 are off
		miniCDF[2]=miniCDF[3]=sum;
	else if(nb[1]+1==nb[2])//case 2 is off
		sum=miniCDF[2]=miniCDF[3];
	else
		sum=miniCDF[2];

	if(-128==nb[1])//case 0 is off, case 1 cannot be off
		miniCDF[1]=0;
	
	//miniCDF[0]=0;//CDF[1]-CDF[0]: -128<x<nb[1]
	//miniCDF[1]=              ok[1]+nb[1]  +128;                //CDF[2]-CDF[1]: x==nb[1]
	//miniCDF[2]=nb[1]+1<nb[2]?ok[2]+nb[1]+1+128:ok[3]+nb[2]+128;//CDF[3]-CDF[2]: nb[1]<x<nb[2]
	//miniCDF[3]=              ok[3]+nb[2]  +128;                //CDF[4]-CDF[3]: x==nb[2]
	//miniCDF[4]=nb[2]+1<nb[3]?ok[4]+nb[2]+1+128:ok[5]+nb[3]+128;//CDF[5]-CDF[4]: nb[2]<x<nb[3]
	//miniCDF[5]=              ok[5]+nb[3]  +128;                //CDF[6]-CDF[5]: x==nb[3]
	//miniCDF[6]=nb[3]+1<nb[4]?ok[6]+nb[3]+1+128:ok[7]+nb[4]+128;//CDF[7]-CDF[6]: nb[3]<x<nb[4]
	//miniCDF[7]=              ok[7]+nb[4]  +128;                //CDF[8]-CDF[7]: x==nb[4]
	//miniCDF[8]=nb[4]+1<128  ?ok[8]+nb[4]+1+128:0x10000;        //CDF[9]-CDF[8]: nb[4]<x<128
	//miniCDF[9]=0x10000;

	//miniCDF[9]=0x10000;
	//miniCDF[8]=nb[4]+1<128  ?ok[8]+nb[4]+1+128:miniCDF[9];//CDF[9]-CDF[8]: nb[4]<x<128
	//miniCDF[7]=              ok[7]+nb[4]  +128;           //CDF[8]-CDF[7]: x==nb[4]
	//miniCDF[6]=nb[3]+1<nb[4]?ok[6]+nb[3]+1+128:miniCDF[7];//CDF[7]-CDF[6]: nb[3]<x<nb[4]
	//miniCDF[5]=              ok[5]+nb[3]  +128;           //CDF[6]-CDF[5]: x==nb[3]
	//miniCDF[4]=nb[2]+1<nb[3]?ok[4]+nb[2]+1+128:miniCDF[5];//CDF[5]-CDF[4]: nb[2]<x<nb[3]
	//miniCDF[3]=              ok[3]+nb[2]  +128;           //CDF[4]-CDF[3]: x==nb[2]
	//miniCDF[2]=nb[1]+1<nb[2]?ok[2]+nb[1]+1+128:miniCDF[3];//CDF[3]-CDF[2]: nb[1]<x<nb[2]
	//miniCDF[1]=              ok[1]+nb[1]  +128;           //CDF[2]-CDF[1]: x==nb[1]
	//miniCDF[0]=0;//CDF[1]-CDF[0]: -128<x<nb[1]

	//int miniCDF[]=
	//{
	//	(int)((long long)ok[0]*0xFF00>>16),            //0
	//	(int)((long long)ok[1]*0xFF00>>16)+nb[1]  +128,//1
	//	(int)((long long)ok[2]*0xFF00>>16)+nb[1]+1+128,//2
	//	(int)((long long)ok[3]*0xFF00>>16)+nb[2]  +128,//3
	//	(int)((long long)ok[4]*0xFF00>>16)+nb[2]+1+128,//4
	//	(int)((long long)ok[5]*0xFF00>>16)+nb[3]  +128,//5
	//	(int)((long long)ok[6]*0xFF00>>16)+nb[3]+1+128,//6
	//	(int)((long long)ok[7]*0xFF00>>16)+nb[4]  +128,//7
	//	(int)((long long)ok[8]*0xFF00>>16)+nb[4]+1+128,//8
	//	0x10000,//9
	//};
}
int    e10_encode_ch(const unsigned char *src, int bw, int bh, int kc, ArrayHandle *data, int loud)
{
	int res=bw*bh;
	int histlen=(729LL*9+729LL*3)*sizeof(long long),//9 histograms with 729 cases + 729 sets of means for 3 gradients
		overheadlen=75LL*10*sizeof(short);

	unsigned char *b2=(unsigned char*)malloc((size_t)res<<2);
	long long *hist=(long long*)malloc(histlen);
	unsigned short *overhead=(unsigned short*)malloc(overheadlen);
	if(!b2||!hist||!overhead)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(b2, src, (size_t)res<<2);

#if 1
	addbuf(b2, bw, bh, 3, 4, 128);
	colortransform_ycocgt_fwd((char*)b2, bw, bh);
	//pred_grad_fwd((char*)b2, bw, bh, 3, 4);
	addbuf(b2, bw, bh, 3, 4, 128);
#endif
	if(loud)//
		printf("Channel %d:\n", kc);

	memset(hist, 0, histlen);

	e10_iter(b2, bw, bh, kc, hist);
	long long *grad=hist+729LL*9;
	int overheadcount=0;
	for(int ks=0;ks<729;++ks)
	{
		long long *hk=hist+9*ks, *gk=grad+3*ks;
		int sum=0;
		for(int k=0;k<9;++k)
			sum+=(int)hk[k];
		if(sum)
		{
			if(loud)//
			{
				int rem=256-(int)((gk[0]+gk[1]+gk[2])/sum);
				int ranges[]=
				{
					rem>>1,
					1,
					(int)(gk[0]/sum),
					1,
					(int)(gk[1]/sum),
					1,
					(int)(gk[2]/sum),
					1,
					rem-(rem>>1),
				};
				printf("[%3d] %7d ", ks, sum);
				for(int k=0;k<9;++k)
				{
					if(hk[k])
						printf(" %7d/%2d", (int)hk[k], ranges[k]);
					else
						printf(" -------/%2d", ranges[k]);
				}
				printf("\n");
			}
			if(overheadcount>=75)
			{
				LOG_ERROR("Unexpected number of cases %d", ks);
				return 0;
			}
			unsigned short *ok=overhead+10*overheadcount;
			ok[0]=ks;
			int sum2=0;
			for(int k=0;k<9;++k)
			{
				int freq=(int)(hk[k]*0xFF00/sum);
				ok[1+k]=sum2;
				sum2+=freq;
			}
			++overheadcount;
		}
	}
	free(hist);

	int casehist[9]={0};

	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, 5);//int csize, char overheadcount
	dlist_push_back(&list, overhead, overheadlen);
	unsigned state=0x10000;
	for(int ky=bh-1;ky>=0;--ky)
	{
		for(int kx=bw-1;kx>=0;--kx)
		{
			int idx=bw*ky+kx;
			char
				topleft =kx-1>=0&&ky-1>=0?b2[(idx-bw-1)<<2|kc]-128:0,
				top     =         ky-1>=0?b2[(idx-bw  )<<2|kc]-128:0,
				topright=kx+1<bw&&ky-1>=0?b2[(idx-bw+1)<<2|kc]-128:0,
				left    =kx-1>=0         ?b2[(idx   -1)<<2|kc]-128:0,
				curr    =                 b2[ idx      <<2|kc]-128  ;
			int nb[]={-128, topleft, top, topright, left, 128};
			
			//if(kc==0&&ky==1&&kx==20)//
			//	kx=20;

			int miniCDF[10];
			e10_getminiCDF(nb, overhead, overheadcount, miniCDF);

			int caseidx;
				 if(curr<nb[1])		caseidx=0;
			else if(curr==nb[1])	caseidx=1;
			else if(curr<nb[2])		caseidx=2;
			else if(curr==nb[2])	caseidx=3;
			else if(curr<nb[3])		caseidx=4;
			else if(curr==nb[3])	caseidx=5;
			else if(curr<nb[4])		caseidx=6;
			else if(curr==nb[4])	caseidx=7;
			else					caseidx=8;

			++casehist[caseidx];

			int cdf, freq, base;
			if(caseidx&1)//curr == nb[1+(caseidx>>1)]
			{
				cdf=miniCDF[caseidx];
				freq=(caseidx<9-1?miniCDF[caseidx+1]:0x10000)-cdf;
			}
			else//nb[caseidx>>1] <[=] curr < nb[1+(caseidx>>1)]		if caseidx==0 [<=] else [<]
			{
				base=nb[caseidx>>1]+(caseidx!=0);//open range if case > 0
				freq=(miniCDF[caseidx+1]-miniCDF[caseidx])/(nb[(caseidx>>1)+1]-base);
				cdf=miniCDF[caseidx]+(curr-base)*freq;
				//cdf=miniCDF[caseidx]+(curr-nb[caseidx>>1])*(miniCDF[caseidx+1]-miniCDF[caseidx])/den;
			}
			//cdf=(int)((long long)cdf*0xFF00>>16)+curr;//guard
			//freq=(int)((long long)freq*0xFF00>>16)+1;

			if(!freq)
				LOG_ERROR("ZPS");
						
			if(state>=(unsigned)(freq<<16))//renorm
			{
				dlist_push_back(&list, &state, 2);
				state>>=16;
			}
			debug_enc_update(state, cdf, freq, kx, ky, 0, kc, curr);
			state=state/freq<<16|(cdf+state%freq);//update
		}
	}
	dlist_push_back(&list, &state, 4);
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, &list.nobj, 4);
	memcpy(data[0]->data+dststart+4, &overheadcount, 1);

	if(loud)
	{
		printf("ch %d  usize %7d  csize %7d  CR %lf\n", kc, res, (int)list.nobj, (double)res/list.nobj);
		for(int k=0;k<9;++k)
			printf("\tcase %d %7d\n", k, casehist[k]);
		printf("\n\n");
	}

	dlist_clear(&list);
	free(overhead);
	free(b2);
	return 1;
}
size_t e10_decode_ch(const unsigned char *data, size_t datastart, size_t datalen, int bw, int bh, int kc, unsigned char *buf)
{
	int res=bw*bh;
	int overheadlen=75LL*10*sizeof(short);

	unsigned short *overhead=(unsigned short*)malloc(overheadlen);
	if(!overhead)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	const unsigned char *srcstart=data+datastart, *srcptr;

	int csize=0, overheadcount=0;
	memcpy(&csize, srcstart, 4);
	memcpy(&overheadcount, srcstart+4, 1);
	if(datastart+5+overheadcount*10*sizeof(short)>datalen)
	{
		LOG_ERROR("Unexpected EOF");
		return 0;
	}
	memcpy(overhead, srcstart+5, overheadcount*10*sizeof(short));

	unsigned state;
	srcstart+=5+overheadcount*10*sizeof(short);
	srcptr=data+datastart+csize-4;
	if(srcptr<srcstart)
	{
		LOG_ERROR("ANS overflow");
		return 0;
	}
	memcpy(&state, srcptr, 4);
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
		{
			int idx=bw*ky+kx;
			char
				topleft =kx-1>=0&&ky-1>=0?buf[(idx-bw-1)<<2|kc]-128:0,
				top     =         ky-1>=0?buf[(idx-bw  )<<2|kc]-128:0,
				topright=kx+1<bw&&ky-1>=0?buf[(idx-bw+1)<<2|kc]-128:0,
				left    =kx-1>=0         ?buf[(idx   -1)<<2|kc]-128:0;
			int nb[]={-128, topleft, top, topright, left, 128};

			//if(kc==0&&ky==1&&kx==20)//
			//	kx=20;
			
			int miniCDF[10];
			e10_getminiCDF(nb, overhead, overheadcount, miniCDF);

			int c=(unsigned short)state;
			int caseidx=0;
			int L=0, R=10, found=0;
			while(L<=R)
			{
				caseidx=(L+R)>>1;
				if(miniCDF[caseidx]<c)
					L=caseidx+1;
				else if(miniCDF[caseidx]>c)
					R=caseidx-1;
				else
				{
					found=1;
					break;
				}
			}
			if(!found)
				caseidx=L+(L<10&&miniCDF[L]<c)-1;
			else
				for(;caseidx<10-1&&miniCDF[caseidx+1]==c;++caseidx);

			int curr, cdf, freq, base;
			if(caseidx&1)
			{
				curr=nb[(caseidx>>1)+1];
				cdf=miniCDF[caseidx];
				freq=miniCDF[caseidx+1]-cdf;
			}
			else
			{
				base=nb[caseidx>>1]+(caseidx!=0);//open range if case > 0
				freq=(miniCDF[caseidx+1]-miniCDF[caseidx])/(nb[(caseidx>>1)+1]-base);
				curr=base+(c-miniCDF[caseidx])/freq;
				if(curr>nb[(caseidx>>1)+1]-1)
					curr=nb[(caseidx>>1)+1]-1;
				cdf=miniCDF[caseidx]+(curr-base)*freq;
			}

			buf[idx<<2|kc]=curr+128;

			debug_dec_update(state, cdf, freq, kx, ky, 0, kc, curr);
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
	return datastart+csize;
}

static void   ans_enc_update(unsigned *state, int cdf, int freq, DList *list
#ifdef DEBUG_ANS
	, int kc, int kx, int ky, int curr
#endif
)
{
	if(!freq)
		LOG_ERROR("ZPS");
	if(*state>=(unsigned)(freq<<16))//renorm
	{
		dlist_push_back(list, state, 2);
		*state>>=16;
	}
	debug_enc_update(*state, cdf, freq, kx, ky, 0, kc, curr);
	*state=*state/freq<<16|(cdf+*state%freq);//update
}
static int    ans_CDF2sym(int c, int *CDF, int start, int end, int norm)
{
	int sym=0;
	int L=start, R=end, found=0;
	int vmin=CDF[start], range=CDF[end]-vmin, test;
	//vmin=(int)((long long)vmin*norm/range)+start;
	while(L<=R)
	{
		sym=(L+R)>>1;
		test=CDF[sym];
		if(norm)
			test=(int)((long long)(test-vmin)*norm/range)+sym-start;
		if(test<c)
			L=sym+1;
		else if(test>c)
			R=sym-1;
		else
		{
			found=1;
			break;
		}
	}
	if(!found)
		sym=L+(L<end&&(norm?(int)((long long)(CDF[L]-vmin)*norm/range)+L-start:CDF[L])<c)-1;
	else
		for(;sym<end-1&&(norm?(int)((long long)(CDF[sym+1]-vmin)*norm/range)+sym+1-start:CDF[sym+1])==c;++sym);
	return sym;
}
static void ans_dec_update(unsigned *state, int cdf, int freq, const unsigned char *srcstart, const unsigned char **srcptr)
{
	*state=freq*(*state>>16)+(unsigned short)*state-cdf;//update
	if(*state<0x10000)//renorm
	{
		*state<<=16;
		if(*srcptr-2>=srcstart)
		{
			*srcptr-=2;
			memcpy(&*state, *srcptr, 2);
		}
	}
}
int    e10dash_encode_ch(const unsigned char *src, int bw, int bh, int kc, int alpha, int blocksize, int margin, ArrayHandle *data, int loud)
{
	int res=bw*bh;
	int histlen=(729LL*9+729LL*3)*sizeof(long long),//9 histograms with 729 cases + 729 sets of means for 3 gradients
		overheadlen=75LL*10*sizeof(short);

	unsigned char *b2=(unsigned char*)malloc((size_t)res<<2);
	long long *hist=(long long*)malloc(histlen);
	unsigned short *overhead=(unsigned short*)malloc(overheadlen);
	unsigned short *CDF0=(unsigned short*)malloc(256LL*sizeof(short));
	int *CDF2=(int*)malloc(257*sizeof(int));
	if(!b2||!hist||!overhead||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	memcpy(b2, src, (size_t)res<<2);

#if 1
	addbuf(b2, bw, bh, 3, 4, 128);
	colortransform_ycocgt_fwd((char*)b2, bw, bh);
	pred_grad_fwd((char*)b2, bw, bh, 3, 4);
	addbuf(b2, bw, bh, 3, 4, 128);
#endif

	int *hist2=(int*)hist;
	memset(hist2, 0, 256LL*sizeof(unsigned));
	for(int k=0;k<res;++k)
	{
		unsigned char sym=b2[k<<2|kc];
		++hist2[kc<<8|sym];
	}
	normalize_histogram(hist2, 256, res, CDF0);

	double csize_cases=0, csize_remainder=0, p;//
	int remcount=0;//

	//if(loud)//
	//	printf("Channel %d:\n", kc);

	memset(hist, 0, histlen);

	e10_iter(b2, bw, bh, kc, hist);
	long long *grad=hist+729LL*9;
	int overheadcount=0;
	for(int ks=0;ks<729;++ks)
	{
		long long *hk=hist+9*ks, *gk=grad+3*ks;
		int sum=0;
		for(int k=0;k<9;++k)
			sum+=(int)hk[k];
		if(sum)
		{
#if 0
			if(loud)//
			{
				int rem=256-(int)((gk[0]+gk[1]+gk[2])/sum);
				int ranges[]=
				{
					rem>>1,
					1,
					(int)(gk[0]/sum),
					1,
					(int)(gk[1]/sum),
					1,
					(int)(gk[2]/sum),
					1,
					rem-(rem>>1),
				};
				printf("[%3d] %7d ", ks, sum);
				for(int k=0;k<9;++k)
				{
					if(hk[k])
						printf(" %7d/%2d", (int)hk[k], ranges[k]);
					else
						printf(" -------/%2d", ranges[k]);
				}
				printf("\n");
			}
#endif
			if(overheadcount>=75)
			{
				LOG_ERROR("Unexpected number of cases %d", ks);
				return 0;
			}
			unsigned short *ok=overhead+10*overheadcount;
			ok[0]=ks;
			int sum2=0;
			for(int k=0;k<9;++k)
			{
				int freq=(int)(hk[k]*0xFF00/sum);
				ok[1+k]=sum2;
				sum2+=freq;
			}
			++overheadcount;
		}
	}
	overheadlen=overheadcount*10LL*sizeof(short);
	free(hist);

	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, 5);//int csize, char overheadcount
	dlist_push_back(&list, CDF0, 256LL*sizeof(short));
	dlist_push_back(&list, overhead, overheadlen);

	int bxcount=(bw+blocksize-1)/blocksize,
		bycount=(bh+blocksize-1)/blocksize;
	unsigned state=0x10000;
	for(int by=bycount-1;by>=0;--by)
	{
		int ky=by*blocksize;
		for(int bx=bxcount-1;bx>=0;--bx)
		{
			int kx=bx*blocksize;
			int xend=0, yend=0;
			t16_prepblock(b2, CDF0, bw, bh, kc, bx, by, alpha, blocksize, margin, CDF2, &xend, &yend);
			for(int ky2=yend-1;ky2>=ky;--ky2)
			{
				for(int kx2=xend-1;kx2>=kx;--kx2)
				{
					int idx=bw*ky2+kx2;
					int
						topleft =kx2-1>=0  &&ky2-1>=0?b2[(idx-bw-1)<<2|kc]-128:0,
						top     =            ky2-1>=0?b2[(idx-bw  )<<2|kc]-128:0,
						topright=kx2+1<xend&&ky2-1>=0?b2[(idx-bw+1)<<2|kc]-128:0,
						left    =kx2-1>=0            ?b2[(idx   -1)<<2|kc]-128:0,
						curr    =                     b2[ idx      <<2|kc]-128  ;
					int nb[]={-129, topleft, top, topright, left, 128};
					
					//if(ky2==(bh>>1)&&kx2==(bw>>1))//
					//	kx2=bw>>1;
					//if(kc==1&&ky2==0&&kx2==0)//
					//	kx2=0;

					int miniCDF[10];
					e10_getminiCDF(nb, overhead, overheadcount, miniCDF);

					int caseidx, cdf, freq, start, end, range, norm;
						 if(curr<nb[1])		caseidx=0;
					else if(curr==nb[1])	caseidx=1;
					else if(curr<nb[2])		caseidx=2;
					else if(curr==nb[2])	caseidx=3;
					else if(curr<nb[3])		caseidx=4;
					else if(curr==nb[3])	caseidx=5;
					else if(curr<nb[4])		caseidx=6;
					else if(curr==nb[4])	caseidx=7;
					else					caseidx=8;
						 
					if(!(caseidx&1))
					{
						start=nb[caseidx>>1]+1+128;
						end=nb[(caseidx>>1)+1]+128;
						if(start+1<end)
						{
							norm=0x10000-(end-start);
							range=CDF2[end]-CDF2[start];
							curr+=128;

							cdf=(int)((long long)(CDF2[curr]-CDF2[start])*norm/range)+curr-start;
							freq=(int)((long long)(CDF2[curr+1]-CDF2[curr])*norm/range)+1;

							//cdf=(int)((long long)CDF2[curr]*norm/range)+curr-start;
							//freq=(int)((long long)CDF2[curr+1]*norm/range)+curr+1-start-cdf;
							//cdf-=((int)((long long)CDF2[start]*norm/range)+start);

							//cdf=(int)((long long)(CDF2[curr+128]-CDF2[base+128])*norm/range)+curr-base;
							//freq=(int)((long long)(CDF2[curr+128+1]-CDF2[curr+128])*norm/range)+1;
							
							p=(double)freq/0x10000;//
							csize_remainder-=log2(p);
							++remcount;

							ans_enc_update(&state, cdf, freq, &list
#ifdef DEBUG_ANS
								, kc, kx2, ky2, curr
#endif
							);
						}
					}

					cdf=miniCDF[caseidx], freq=miniCDF[caseidx+1]-cdf;

					p=(double)freq/0x10000;//
					csize_cases-=log2(p);

					ans_enc_update(&state, cdf, freq, &list
#ifdef DEBUG_ANS
						, kc, kx2, ky2, caseidx
#endif
					);

				}
			}
		}
	}
	dlist_push_back(&list, &state, 4);
	size_t dststart=dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, &list.nobj, 4);
	memcpy(data[0]->data+dststart+4, &overheadcount, 1);

	if(loud)
	{
		csize_cases/=8;
		csize_remainder/=8;
		printf("ch %d  usize %7d  csize %7d  CR %lf  cases %.2lf  rem %.2lf / %d  total %.2lf\n\n", kc, res, (int)list.nobj, (double)res/list.nobj, csize_cases, csize_remainder, remcount, csize_cases+csize_remainder);
	}

	dlist_clear(&list);
	free(overhead);
	free(b2);
	return 1;
}
size_t e10dash_decode_ch(const unsigned char *data, size_t datastart, size_t datalen, int bw, int bh, int kc, int alpha, int blocksize, int margin, unsigned char *buf)
{
	int res=bw*bh;
	int overheadlen=75LL*10*sizeof(short);

	if(datastart+5+256*sizeof(short)>datalen)
	{
		LOG_ERROR("Unexpected EOF");
		return 0;
	}

	unsigned short *overhead=(unsigned short*)malloc(overheadlen);
	unsigned short *CDF0=(unsigned short*)malloc(256LL*sizeof(short));
	int *CDF2=(int*)malloc(257*sizeof(int));
	if(!overhead||!CDF0||!CDF2)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	const unsigned char *srcstart=data+datastart, *srcptr;

	int csize=0, overheadcount=0;
	memcpy(&csize, srcstart, 4);
	memcpy(&overheadcount, srcstart+4, 1);
	if(datastart+5+overheadcount*10*sizeof(short)>datalen)
	{
		LOG_ERROR("Unexpected EOF");
		return 0;
	}

	memcpy(CDF0, srcstart+5, 256LL*sizeof(short));
	memcpy(overhead, srcstart+5+256LL*sizeof(short), overheadcount*10LL*sizeof(short));

	unsigned state;
	srcstart+=5+256LL*sizeof(short)+overheadcount*10LL*sizeof(short);
	srcptr=data+datastart+csize-4;
	if(srcptr<srcstart)
	{
		LOG_ERROR("ANS overflow");
		return 0;
	}
	memcpy(&state, srcptr, 4);
	int bxcount=(bw+blocksize-1)/blocksize,
		bycount=(bh+blocksize-1)/blocksize;
	for(int by=0;by<bycount;++by)
	{
		int ky=by*blocksize;
		for(int bx=0;bx<bxcount;++bx)//for each block
		{
			int kx=bx*blocksize;
			int xend=0, yend=0;
			t16_prepblock(buf, CDF0, bw, bh, kc, bx, by, alpha, blocksize, margin, CDF2, &xend, &yend);
			for(int ky2=ky;ky2<yend;++ky2)
			{
				for(int kx2=kx;kx2<xend;++kx2)//for each pixel
				{
					int idx=bw*ky2+kx2;
					int
						topleft =kx2-1>=0  &&ky2-1>=0?buf[(idx-bw-1)<<2|kc]-128:0,
						top     =            ky2-1>=0?buf[(idx-bw  )<<2|kc]-128:0,
						topright=kx2+1<xend&&ky2-1>=0?buf[(idx-bw+1)<<2|kc]-128:0,
						left    =kx2-1>=0            ?buf[(idx   -1)<<2|kc]-128:0;
					int nb[]={-129, topleft, top, topright, left, 128};
					
					//if(kc==1&&ky2==0&&kx2==0)//
					//	kx2=0;
			
					int miniCDF[10], cdf, freq, curr;
					e10_getminiCDF(nb, overhead, overheadcount, miniCDF);

					int caseidx=ans_CDF2sym((unsigned short)state, miniCDF, 0, 9, 0);

					cdf=miniCDF[caseidx];
					freq=miniCDF[caseidx+1]-cdf;
					debug_dec_update(state, cdf, freq, kx2, ky2, 0, kc, caseidx);
					ans_dec_update(&state, cdf, freq, srcstart, &srcptr);

					int start, end, range, norm;
					if(caseidx&1)
						curr=nb[(caseidx>>1)+1]+128;
					else
					{
						start=nb[caseidx>>1]+1+128;
						end=nb[(caseidx>>1)+1]+128;
						if(start==end-1)
							curr=start;
						else
						{
							norm=0x10000-(end-start);

							curr=ans_CDF2sym((unsigned short)state, CDF2, start, end, norm);

							range=CDF2[end]-CDF2[start];
							cdf=(int)((long long)(CDF2[curr]-CDF2[start])*norm/range)+curr-start;
							freq=(int)((long long)(CDF2[curr+1]-CDF2[curr])*norm/range)+1;

							//cdf=(int)((long long)CDF2[curr]*norm/range)+curr - ((int)((long long)CDF2[start]*norm/range)-start);
							//freq=(int)((long long)CDF2[curr+1]*norm/range) - (int)((long long)CDF2[curr]*norm/range) + 1;

							//cdf=(int)((long long)(CDF2[curr]-CDF2[base+128])*norm/range)+curr-base;
							//freq=(int)((long long)(CDF2[curr+1]-CDF2[curr])*norm/range)+1;

							debug_dec_update(state, cdf, freq, kx2, ky2, 0, kc, curr);
							ans_dec_update(&state, cdf, freq, srcstart, &srcptr);
						}
					}
					buf[idx<<2|kc]=curr;

					//if(ky2==(bh>>1)&&kx2==(bw>>1))//
					//	kx2=bw>>1;
				}
			}
		}
	}
	free(overhead);
	free(CDF0);
	free(CDF2);
	return datastart+csize;
}
#endif