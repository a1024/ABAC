#include"battle.h"
#include<stdio.h>//for debugging
#include<stdlib.h>
#include<string.h>
#include<smmintrin.h>//SSE4.1 _mm_mullo_epi32
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#include"rans_common.h"
static const char file[]=__FILE__;

//	#define ENABLE_CHECK_DECODE//debug
	#define PROF(...)//

static const int tag_rans0a='A'|'N'<<8|'0'<<16|'A'<<24;

//#define WRITE_GUARD(N)	if(p+N>max_size)bytes=(unsigned char*)realloc(bytes, max_size<<=1)
int rans_sse2_encode(const void *src, size_t nbytes, int symbytes, int is_signed, ArrayHandle *out)
{
	const int
		infosize=ANS_NLEVELS*sizeof(SymbolInfo), lginfosize=13,//log2(256*32)
		lghistsize=9;
	DList list;
	const unsigned char *buf=(const unsigned char*)src;
	size_t dstidx;
	SymbolInfo *info;
	int internalheadersize=4+8+symbytes*ANS_NLEVELS*sizeof(short);
	int chmask=symbytes-1;

	if(nbytes&15)
		return RANS_INVALID_NBYTES;
	if(symbytes&chmask)
		return RANS_INVALID_SYMBYTES;
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
	memcpy(out[0]->data+dstidx, &tag_rans0a, 4);
	dstidx+=(size_t)(4+8);
	for(int kc=0;kc<symbytes;++kc)
		rans_calc_histogram(buf+kc, (int)(nbytes/symbytes), symbytes, out[0]->data+dstidx+((size_t)kc<<lghistsize), ANS_PROB_BITS, 0);
	
#ifndef ENABLE_CHECK_DECODE
	rans_prep(out[0]->data+dstidx, symbytes, &info, 0, 0);
#else
	unsigned char *CDF2sym=0;
	rans_prep(out[0]->data+dstidx, symbytes, &info, &CDF2sym, 0);
#endif

	dstidx-=8;
	dlist_init(&list, 1, 1024, 0);

	long long t1=__rdtsc();
	
	__m128i negtwo=_mm_set1_epi8(0xFE), one=_mm_set1_epi8(1);
	__m128i sr32=_mm_set_epi8(
		-1, -1, -1, -1,
		15, 14, 13, 12,
		-1, -1, -1, -1,
		 7,  6,  5,  4
	);
	__m128i sr16=_mm_set_epi8(
		-1, -1, -1, -1,
		13, 12, 11, 10,
		-1, -1, -1, -1,
		 5,  4,  3,  2
	);
	__m128i sr16_hi=_mm_set_epi8(
		13, 12, 11, 10,
		-1, -1, -1, -1,
		 5,  4,  3,  2,
		-1, -1, -1, -1
	);
	__m128i state=_mm_set1_epi32(0x10000);
	ALIGN(16) unsigned a_state[4];
	ALIGN(16) unsigned char a_val[16];
	for(ptrdiff_t ks=nbytes-16;ks>=0;ks-=16)
	{
		const unsigned char *ptr;
		if(is_signed)
		{
			__m128i val=_mm_loadu_si128((__m128i*)(buf+ks));

			__m128i neg=_mm_cmplt_epi8(val, _mm_setzero_si128());//neg=v<0
			val=_mm_xor_si128(val, neg);//v^=-neg
			neg=_mm_and_si128(neg, one);
			val=_mm_add_epi8(val, neg); //v+=neg
			
			val=_mm_slli_epi16(val, 1); //v<<=1;
			val=_mm_and_si128(val, negtwo);

			val=_mm_or_si128(val, neg); //v|=neg;

			_mm_store_si128((__m128i*)a_val, val);
			ptr=a_val+16-4;
		}
		else
			ptr=buf+ks+16-4;

		for(int kc=16-4;kc>=0;kc-=4, ptr-=4)
		{
			size_t abs_idx=ks|kc;
			SymbolInfo *ip[4]=
			{
				info+(((abs_idx  )&chmask)<<ANS_DEPTH)+ptr[0],
				info+(((abs_idx|1)&chmask)<<ANS_DEPTH)+ptr[1],
				info+(((abs_idx|2)&chmask)<<ANS_DEPTH)+ptr[2],
				info+(((abs_idx|3)&chmask)<<ANS_DEPTH)+ptr[3],
			};
			__m128i
				invf_lo=_mm_set_epi32(ip[3]->invfcomp[0], ip[2]->invfcomp[0], ip[1]->invfcomp[0], ip[0]->invfcomp[0]),
				invf_hi=_mm_set_epi32(ip[3]->invfcomp[1], ip[2]->invfcomp[1], ip[1]->invfcomp[1], ip[0]->invfcomp[1]),
				negf=_mm_set_epi32(ip[3]->neg_freq, ip[2]->neg_freq, ip[1]->neg_freq, ip[0]->neg_freq),
				cdf=_mm_set_epi32(ip[3]->CDF, ip[2]->CDF, ip[1]->CDF, ip[0]->CDF);

			_mm_store_si128((__m128i*)a_state, state);
			for(int k2=3;k2>=0;--k2)
			{
				if(a_state[k2]>=ip[k2]->renorm_limit)
				{
					dlist_push_back(&list, a_state+k2, 2);
					//WRITE_GUARD(2);
					//memcpy(bytes+p, a_state+k2, 2);
					//p+=2;

					a_state[k2]>>=16;
				}
			}
			state=_mm_load_si128((__m128i*)a_state);
#ifdef ENABLE_CHECK_DECODE
			ALIGN(16) unsigned debug_state1[4];
			memcpy(debug_state1, a_state, sizeof(a_state));
#endif
			
			//mul_sr48(state, invf)		bits [32+48:48] of full 96bit result
			__m128i prod0=_mm_mul_epu32(state, invf_lo);
			__m128i prod1=_mm_mul_epu32(state, invf_hi);
			__m128i u0=_mm_alignr_epi8(state, state, 4);
			__m128i u1=_mm_alignr_epi8(invf_lo, invf_lo, 4);
			__m128i u2=_mm_alignr_epi8(invf_hi, invf_hi, 4);
			__m128i prod2=_mm_mul_epu32(u0, u1);
			__m128i prod3=_mm_mul_epu32(u0, u2);

			prod0=_mm_shuffle_epi8(prod0, sr32);//shift-in with zeros
			prod0=_mm_add_epi64(prod0, prod1);
			prod0=_mm_shuffle_epi8(prod0, sr16);//also clear high 32 bits

			prod2=_mm_shuffle_epi8(prod2, sr32);
			prod2=_mm_add_epi64(prod2, prod3);
			prod2=_mm_shuffle_epi8(prod2, sr16_hi);

			prod0=_mm_or_si128(prod0, prod2);

			//state += prod0*negf + CDF
			prod0=_mm_mullo_epi32(prod0, negf);//SSE4.1

			state=_mm_add_epi32(state, prod0);
			state=_mm_add_epi32(state, cdf);
			
			//check decode
#ifdef ENABLE_CHECK_DECODE
			ALIGN(16) unsigned debug_state2[4];
			_mm_store_si128((__m128i*)debug_state2, state);
			for(int k2=0;k2<4;++k2)
			{
				int ch_idx=(int)((abs_idx|k2)&chmask);
				unsigned short c=(unsigned short)debug_state2[k2];
				unsigned char val2=CDF2sym[ch_idx<<ANS_PROB_BITS|c];
				SymbolInfo *p=info+(ch_idx<<ANS_DEPTH|val2);

				unsigned s1=p->freq*(debug_state2[k2]>>ANS_PROB_BITS)+c-p->CDF;

				unsigned char val1=buf[abs_idx|k2];
				if(is_signed)
				{
					int neg=(val1&0x80)!=0;
					val1^=-neg;
					val1+=neg;
					val1<<=1;
					val1|=neg;
				}

				if(s1!=debug_state1[k2]||val2!=val1)
					LOG_ERROR("rANS encode error at %d, state 0x%08X != 0x%08X, val 0x%02X != 0x%02X", abs_idx|k2, s1, debug_state1[k2], val2, val1);
			}
#endif
		}
	}
	dlist_push_back(&list, &state, sizeof(state));
	//WRITE_GUARD(sizeof(state));
	//memcpy(bytes+p, &state, sizeof(state));
	//p+=sizeof(state);

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
#define READ_GUARD(NBYTES, IDX)	if((srcptr-=NBYTES)<srcstart)return decode_error((size_t)srcptr, (size_t)srcstart, IDX)
int rans_sse2_decode(const void *srcdata, size_t srclen, void *dstbuf, size_t nbytes, int symbytes, int is_signed
#ifdef ENABLE_GUIDE
	, unsigned char *guide
#endif
)
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

	if(nbytes&15)
		return RANS_INVALID_NBYTES;
	if(symbytes&chmask)
		return RANS_INVALID_SYMBYTES;
	if(srclen<internalheadersize)
		return RANS_BUFFER_OVERRUN;
	if(memcmp(srcptr, &tag_rans0a, 4))
		return RANS_INVALID_TAG;
	srcptr+=4;

	memcpy(&csize, srcptr, 8);
	srcptr+=8;

	hist=srcptr;
	rans_prep(hist, symbytes, &info, &CDF2sym, 0);
	srcstart=srcptr+((size_t)symbytes<<lghistsize);
	srcptr=srcstart+csize;

	long long t1=__rdtsc();
	
	__m128i
		one=_mm_set1_epi8(1),
		hex7F=_mm_set1_epi8(0x7F), hex80=_mm_set1_epi8(0x80);
	__m128i state;
	ALIGN(16) unsigned sk[4], slo[4], freq[4], cdf[4];
	READ_GUARD(sizeof(sk), 0);
	memcpy(sk, srcptr, sizeof(sk));
	state=_mm_load_si128((__m128i*)sk);
	for(ptrdiff_t ks=0;ks<(ptrdiff_t)nbytes;ks+=16)
	{
		__m128i mc, mfreq, mcdf, lo, hi;

		for(int kc=0;kc<16;kc+=4)
		{
#ifdef ENABLE_CHECK_DECODE
			ALIGN(16) unsigned a_state1[4];
			memcpy(a_state1, sk, sizeof(sk));
#endif
			for(int kd=0;kd<4;++kd)//decode			c=(uint16)state, val=CDF2sym[ch][c]
			{
				ptrdiff_t abs_idx=ks|kc|kd;
				int ch_idx=(int)(abs_idx&chmask);

				unsigned short c=((unsigned short*)sk)[kd<<1];
				unsigned char val=CDF2sym[ch_idx<<ANS_PROB_BITS|c];
				SymbolInfo *p=info+(ch_idx<<ANS_DEPTH|val);

				slo[kd]=c;
				sk[kd]=((unsigned short*)sk)[kd<<1|1];
				freq[kd]=p->freq;
				cdf[kd]=p->CDF;
				pixels[abs_idx]=val;
			}

			state=_mm_load_si128((__m128i*)sk);
			mfreq=_mm_load_si128((__m128i*)freq);
			mc   =_mm_load_si128((__m128i*)slo);
			mcdf =_mm_load_si128((__m128i*)cdf);

			lo=_mm_mullo_epi16(state, mfreq);//update	state = freq*(state>>16) + c - CDF
			hi=_mm_mulhi_epu16(state, mfreq);
			hi=_mm_slli_epi32(hi, 16);
			state=_mm_or_si128(lo, hi);
			state=_mm_add_epi32(state, mc);
			state=_mm_sub_epi32(state, mcdf);

			_mm_store_si128((__m128i*)sk, state);

#ifdef ENABLE_CHECK_DECODE
			for(int kd=0;kd<4;++kd)
			{
				ptrdiff_t abs_idx=ks|kc|kd;
				int ch_idx=(int)(abs_idx&chmask);

				unsigned short c=(unsigned short)a_state1[kd];
				unsigned char val=CDF2sym[ch_idx<<ANS_PROB_BITS|c];
				SymbolInfo *p=info+(ch_idx<<ANS_DEPTH|val);

				unsigned s2=p->freq*(a_state1[kd]>>ANS_PROB_BITS)+c-p->CDF;
				if(is_signed)
				{
					int neg=val&1;
					val>>=1;
					val^=-neg;
					val+=neg;
					val|=neg<<7&-!val;
				}

				if(s2!=sk[kd]||val!=pixels[abs_idx]
#ifdef ENABLE_GUIDE
					||guide&&val!=guide[abs_idx]
#endif
				)
					LOG_ERROR("Decode error at (0x%X|0x%X|0x%X)=0x%X", ks, kc, kd, abs_idx);
			}
#endif

			for(int kd=0;kd<4;++kd)//renormalize
			{
				if(sk[kd]<ANS_L)
				{
					READ_GUARD(2, ks|kc|kd);

					//((unsigned short*)sk)[kd<<1|1]=((unsigned short*)sk)[kd<<1];
					sk[kd]<<=16;

					memcpy(sk+kd, srcptr, 2);
				}
			}
		}
#if 1
		if(is_signed)
		{
			__m128i v=_mm_loadu_si128((__m128i*)(pixels+ks));

			__m128i neg=_mm_and_si128(v, one);//neg=v&1
			neg=_mm_cmpeq_epi8(neg, one);

			v=_mm_srli_epi16(v, 1);//v>>=1
			v=_mm_and_si128(v, hex7F);

			v=_mm_xor_si128(v, neg);//v^=-neg
			v=_mm_sub_epi8(v, neg); //v+=neg

			neg=_mm_and_si128(neg, hex80);//v|=neg<<7&-!val
			neg=_mm_and_si128(neg, _mm_cmpeq_epi8(v, _mm_setzero_si128()));
			v=_mm_or_si128(v, neg);

			_mm_storeu_si128((__m128i*)(pixels+ks), v);
		}
#endif
	}
	free(info);

	t1=__rdtsc()-t1;
	return RANS_SUCCESS;
}
