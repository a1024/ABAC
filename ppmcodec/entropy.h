#pragma once
#ifndef INC_AC_H
#define INC_AC_H
#include"util.h"
#include<stdio.h>
#include<string.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>//_udiv128, _umul128
#endif
#ifdef __cplusplus
extern "C"
{
#endif

	
//	#define AC_VALIDATE

#ifdef AC_VALIDATE
void acval_enc(int sym, int cdf, int freq, unsigned long long lo1, unsigned long long hi1, unsigned long long lo2, unsigned long long hi2, unsigned long long cache, int nbits);
void acval_dec(int sym, int cdf, int freq, unsigned long long lo1, unsigned long long hi1, unsigned long long lo2, unsigned long long hi2, unsigned long long cache, int nbits, unsigned long long code);
void acval_dump();
extern int acval_disable;
#ifdef AC_IMPLEMENTATION
typedef struct ACVALStruct
{
	unsigned long long
		lo1, hi1,
		lo2, hi2,
		cache, code;
	int nbits;
	int sym, cdf, freq;
} ACVAL;
#define ACVAL_CBUFSIZE 128
ArrayHandle acval=0;
ACVAL acval_cbuf[ACVAL_CBUFSIZE]={0};
int acval_idx=0;
int acval_disable=0;
void acval_enc(int sym, int cdf, int freq, unsigned long long lo1, unsigned long long hi1, unsigned long long lo2, unsigned long long hi2, unsigned long long cache, int nbits)
{
	ACVAL val=
	{
		lo1, hi1,
		lo2, hi2,
		cache, 0,
		nbits,
		sym, cdf, freq,
	};
	if(acval_disable)
		return;
	//if(lo1>>48||hi1>>48||lo2>>48||hi2>>48)
	//	LOG_ERROR2("");
	if(!acval)
		ARRAY_ALLOC(ACVAL, acval, 0, 0, 0, 0);
	ARRAY_APPEND(acval, &val, 1, 1, 0);
}
void acval_print(int idx, ACVAL const *p, int dec)
{
	printf("%9d sym %04X cdf %08X freq %08X range %016llX~%016llX:%016llX -> %016llX~%016llX:%016llX",
		idx, p->sym, p->cdf, p->freq, p->lo1, p->hi1, p->hi1-p->lo1, p->lo2, p->hi2, p->hi2-p->lo2
	);
	//printf("%9d sym %03X cdf %04X freq %04X range %012llX~%012llX:%012llX->%012llX~%012llX:%012llX cache %016llX nbits %2d",
	//	idx, p->sym, p->cdf, p->freq, p->lo1, p->hi1, p->hi1-p->lo1, p->lo2, p->hi2, p->hi2-p->lo2, p->cache, p->nbits
	////	idx, p->sym, p->cdf, p->freq, p->lo1, p->hi1, p->lo2, p->hi2, p->cache, p->nbits
	//);
	if(dec)
	{
		printf(" code %012llX ", p->code);
		//print_binn(p->code, 48);
	}
	printf("\n");
}
void acval_dump()
{
	printf("\nEnc:\n");
	for(int k=MAXVAR(acval_idx-ACVAL_CBUFSIZE, 0), end=MINVAR(acval_idx+ACVAL_CBUFSIZE, (int)acval->count);k<end;++k)
	{
		ACVAL *p=(ACVAL*)array_at(&acval, k);
		if(k==acval_idx)
			printf("\n");
		acval_print(k, p, 0);
		if(k==acval_idx)
			printf("\n");
	}
	printf("\nDec:\n");
	for(int k=(acval_idx+1)%ACVAL_CBUFSIZE, k2=acval_idx-(ACVAL_CBUFSIZE-1), end=acval_idx%ACVAL_CBUFSIZE;;k=(k+1)%ACVAL_CBUFSIZE, ++k2)
	{
		if(k2<0)
			continue;
		if(k==end)
			printf("\n");
		acval_print(k2, acval_cbuf+k, 1);
		if(k==end)
			break;
	}
}
void acval_dec(int sym, int cdf, int freq, unsigned long long lo1, unsigned long long hi1, unsigned long long lo2, unsigned long long hi2, unsigned long long cache, int nbits, unsigned long long code)
{
	ACVAL *val, val2=
	{
		lo1, hi1,
		lo2, hi2,
		cache, code,
		nbits,
		sym, cdf, freq,
	};
	
	if(acval_disable)
		return;

	acval_cbuf[acval_idx%ACVAL_CBUFSIZE]=val2;

	if(acval_idx>=(int)acval->count)
		LOG_ERROR2("AC validation index error");

	val=(ACVAL*)array_at(&acval, acval_idx);
	//if(lo1>>48||hi1>>48||lo2>>48||hi2>>48)
	//	LOG_ERROR2("");
	if(sym!=val->sym||cdf!=val->cdf||freq!=val->freq||lo1!=val->lo1||hi1!=val->hi1||lo2!=val->lo2||hi2!=val->hi2)
	{
		acval_dump();
		LOG_ERROR2("AC validation error");
	}
	++acval_idx;
}
#endif
#else
#define acval_enc(...)
#define acval_dec(...)
#endif


//arithmetic coder with up to 31-bit stats/renorms, uses 64-bit registers
#define AC3_PROB_BITS 16			//16-bit stats are best
#define AC3_RENORM 32	//multiple of 8!	//32-bit renorm is best
#if AC3_RENORM<=AC3_PROB_BITS
#define AC3_RENORM_STATEMENT while
#else
#define AC3_RENORM_STATEMENT if
#endif
typedef struct _AC3
{
	unsigned long long low, range;
	DList *dst;
	const unsigned char *srcptr, *srcend;
	unsigned long long code;
} AC3;
INLINE void ac3_enc_init(AC3 *ec, DList *dst)
{
	memset(ec, 0, sizeof(*ec));
	ec->low=0;
	ec->range=~0LLU;
	ec->dst=dst;
}
INLINE void ac3_dec_init(AC3 *ec, const unsigned char *start, unsigned const char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->low=0;
	ec->range=~0LLU;
	ec->srcptr=start;
	ec->srcend=end;
	
	if(ec->srcptr+8>ec->srcend)
	{
		LOG_ERROR2("buffer overflow");
		return;
	}
	for(int k=0;k<8;k+=AC3_RENORM/8)
		memcpy((unsigned char*)&ec->code+8-AC3_RENORM/8-k, ec->srcptr+k, AC3_RENORM/8);
	ec->srcptr+=8;
}
INLINE void ac3_enc_renorm(AC3 *ec)//fast renorm by F. Rubin 1979
{
	unsigned long long rmax;

	dlist_push_back(ec->dst, (unsigned char*)&ec->low+8-AC3_RENORM/8, AC3_RENORM/8);
	ec->range=ec->range<<AC3_RENORM|((1LL<<AC3_RENORM)-1);
	ec->low<<=AC3_RENORM;

	rmax=~ec->low;
	if(ec->range>rmax)//clamp hi to register size after renorm
		ec->range=rmax;
}
INLINE void ac3_dec_renorm(AC3 *ec)//fast renorm by F. Rubin 1979
{
	unsigned long long rmax;

	ec->range=ec->range<<AC3_RENORM|((1LL<<AC3_RENORM)-1);
	ec->code<<=AC3_RENORM;
	ec->low<<=AC3_RENORM;
	if(ec->srcptr+AC3_RENORM/8<=ec->srcend)
		memcpy(&ec->code, ec->srcptr, AC3_RENORM/8);
	ec->srcptr+=AC3_RENORM/8;

	rmax=~ec->low;
	if(ec->range>rmax)
		ec->range=rmax;
}
INLINE void ac3_enc_flush(AC3 *ec)
{
	unsigned long long code=ec->low+ec->range;
	int n=FLOOR_LOG2(ec->low^code);
	if(n>0)
		code&=~((1LL<<n)-1);	//minimize final code
#if 0
	int flushbits=get_lsb_index(code);//FIXME tail-chaining parallel decoders
#endif
	for(int k=8-AC3_RENORM/8;k>=0&&code;k-=AC3_RENORM/8)
	{
		dlist_push_back(ec->dst, (unsigned char*)&code+8-(AC3_RENORM/8), AC3_RENORM/8);
		code<<=AC3_RENORM;
	}
}
INLINE void ac3_enc_update(AC3 *ec, unsigned cdf, unsigned freq)
{
	unsigned long long r2;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!freq)
		LOG_ERROR2("ZPS");
	if(cdf+freq<cdf)
		LOG_ERROR2("Invalid CDF");
#endif
	r2=ec->range>>AC3_PROB_BITS;
	AC3_RENORM_STATEMENT(!r2)//this loop runs twice when freq=1 -> range=0
	{
		ac3_enc_renorm(ec);
		r2=ec->range>>AC3_PROB_BITS;
	}
#ifdef AC3_PREC
	ec->low+=r2*cdf+((ec->range&(~0ULL>>(64-AC3_PROB_BITS)))*cdf>>AC3_PROB_BITS);
	ec->range=r2*freq+((ec->range&(~0ULL>>(64-AC3_PROB_BITS)))*freq>>AC3_PROB_BITS)-1;
#else
	ec->low+=r2*cdf;
	ec->range=r2*freq-1;//must decrement hi because decoder fails when code == hi2
#endif
	acval_enc(0, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
INLINE unsigned ac3_dec_getcdf(AC3 *ec)
{
	unsigned long long r2;
	r2=ec->range>>AC3_PROB_BITS;
	AC3_RENORM_STATEMENT(!r2)
	{
		ac3_dec_renorm(ec);
		r2=ec->range>>AC3_PROB_BITS;
	}
#ifdef AC3_PREC
	{
		unsigned long long diff=ec->code-ec->low;
		return (unsigned)_udiv128(diff>>(64-AC3_PROB_BITS), diff<<AC3_PROB_BITS|((1LL<<AC3_PROB_BITS)-1), ec->range, &diff);
	}
#else
	return (unsigned)((ec->code-ec->low)/r2);
#endif
}
INLINE void ac3_dec_update(AC3 *ec, unsigned cdf, unsigned freq)
{
	unsigned long long r2=ec->range>>AC3_PROB_BITS;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!freq)
		LOG_ERROR2("ZPS");
	if(cdf+freq<cdf)
		LOG_ERROR2("Invalid CDF");
#endif
#ifdef AC3_PREC
	ec->low+=r2*cdf+((ec->range&(~0ULL>>(64-AC3_PROB_BITS)))*cdf>>AC3_PROB_BITS);
	ec->range=r2*freq+((ec->range&(~0ULL>>(64-AC3_PROB_BITS)))*freq>>AC3_PROB_BITS)-1;
#else
	ec->low+=r2*cdf;
	ec->range=r2*freq-1;//must decrement hi because decoder fails when code == hi2
#endif
	acval_dec(0, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
}
INLINE void ac3_enc_update_NPOT(AC3 *ec, unsigned cdf, unsigned freq, unsigned den)
{
	unsigned long long q, r;
#ifdef AC_VALIDATE
	unsigned long long lo0, r0;
	if(!freq)
		LOG_ERROR2("ZPS");
	if(cdf+freq<cdf)
		LOG_ERROR2("Invalid CDF");
#endif
	AC3_RENORM_STATEMENT(ec->range<(unsigned)den)//only when freq=1 -> range=0, this loop runs twice
		ac3_enc_renorm(ec);
#ifdef AC_VALIDATE
	lo0=ec->low, r0=ec->range;
#endif
	q=ec->range/den;
	r=ec->range%den;
#ifdef AC3_PREC
	ec->low+=q*cdf+r*cdf/den;
	ec->range=q*freq+r*freq/den-1;
#else
	ec->low+=r2*cdf;
	ec->range=r2*freq-1;//must decrement hi because decoder fails when code == hi2
#endif
	acval_enc(0, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
INLINE unsigned ac3_dec_getcdf_NPOT(AC3 *ec, unsigned den)
{
	AC3_RENORM_STATEMENT(ec->range<(unsigned)den)
		ac3_dec_renorm(ec);
#ifdef AC3_PREC
	{
		unsigned long long lo, hi, lo0;
		lo0=lo=_umul128(ec->code-ec->low, den, &hi);
		lo+=den-1LL;
	//	lo+=den>>1;
	//	lo+=ec->range>>1;//X
		hi+=lo<lo0;
		return (unsigned)_udiv128(hi, lo, ec->range, &hi);
	}
#else
	return (unsigned)((ec->code-ec->low)/r2);
#endif
}
INLINE void ac3_dec_update_NPOT(AC3 *ec, unsigned cdf, unsigned freq, unsigned den)
{
	unsigned long long q=ec->range/den, r=ec->range%den;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!freq)
		LOG_ERROR2("ZPS");
	if(cdf+freq<cdf)
		LOG_ERROR2("Invalid CDF");
#endif
#ifdef AC3_PREC
	ec->low+=q*cdf+r*cdf/den;
	ec->range=q*freq+r*freq/den-1;
#else
	ec->low+=r2*cdf;
	ec->range=r2*freq-1;//must decrement hi because decoder fails when code == hi2
#endif
	acval_dec(0, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
}
INLINE void ac3_enc_bypass(AC3 *ec, int bypass, int nbits)
{
	unsigned long long r2;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
#endif
	r2=ec->range>>nbits;
	AC3_RENORM_STATEMENT(!r2)
	{
		ac3_enc_renorm(ec);
		r2=ec->range>>nbits;
	}
#ifdef AC3_PREC
	ec->low+=r2*bypass+((ec->range&(~0ULL>>(64-nbits)))*bypass>>nbits);
	ec->range=r2-1;
	//ec->low+=r2*bypass+((ec->range&(~0ULL>>(64-nbits)))*bypass>>nbits);
	//ec->range=r2*1+((ec->range&(~0ULL>>(64-nbits)))*1>>nbits)-1;
#else
	ec->low+=r2*bypass;
	ec->range=r2-1;
#endif
	acval_enc(bypass, bypass, 1, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
INLINE int ac3_dec_bypass(AC3 *ec, int nbits)
{
	unsigned long long r2;
	int bypass;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
#endif
	r2=ec->range>>nbits;
	AC3_RENORM_STATEMENT(!r2)
	{
		ac3_dec_renorm(ec);
		r2=ec->range>>nbits;
	}
#ifdef AC3_PREC
	{
		unsigned long long diff=ec->code-ec->low;
		bypass=(int)_udiv128(diff>>(64-nbits), diff<<nbits|((1LL<<nbits)-1), ec->range, &diff);
	}
	ec->low+=r2*bypass+((ec->range&(~0ULL>>(64-nbits)))*bypass>>nbits);
	ec->range=r2-1;
#else
	bypass=(int)((ec->code-ec->low)/r2);
	ec->low+=r2*bypass;
	ec->range=r2-1;
#endif
	acval_dec(bypass, bypass, 1, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return bypass;
}
INLINE void ac3_enc_bypass_NPOT(AC3 *ec, int bypass, int nlevels)
{
	unsigned long long q, r;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
#endif
	
	AC3_RENORM_STATEMENT(ec->range<(unsigned)nlevels)
		ac3_enc_renorm(ec);
	q=ec->range/nlevels;
	r=ec->range%nlevels;
#ifdef AC3_PREC
	ec->low+=q*bypass+r*bypass/nlevels;
	ec->range=q-1;
#else
	ec->low+=r2*bypass;
	ec->range=r2-1;
#endif
	acval_enc(bypass, bypass, 1, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
INLINE int ac3_dec_bypass_NPOT(AC3 *ec, int nlevels)
{
	unsigned long long q, r;
	int bypass;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
#endif
	
	AC3_RENORM_STATEMENT(ec->range<(unsigned)nlevels)
		ac3_dec_renorm(ec);
#ifdef AC3_PREC
	{
		unsigned long long lo, hi;
		lo=_umul128(ec->code-ec->low, nlevels, &hi);
		hi+=lo+nlevels-1<lo;
		lo+=nlevels-1;
		bypass=(int)_udiv128(hi, lo, ec->range, &hi);
	}
	q=ec->range/nlevels;
	r=ec->range%nlevels;
	ec->low+=q*bypass+r*bypass/nlevels;
	ec->range=q-1;
#else
	bypass=(int)((ec->code-ec->low)*nlevels/ec->range);
	ec->low+=r2*bypass;
	ec->range=r2-1;
#endif
	acval_dec(bypass, bypass, 1, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return bypass;
}
INLINE void ac3_enc(AC3 *ec, int sym, const unsigned *CDF)
{
	unsigned long long r2;
	unsigned cdf, freq;

	r2=ec->range>>AC3_PROB_BITS;
	AC3_RENORM_STATEMENT(!r2)//only when freq=1 -> range=0, this loop runs twice
	{
		ac3_enc_renorm(ec);
		r2=ec->range>>AC3_PROB_BITS;
	}
	cdf=CDF[sym];
	freq=CDF[sym+1]-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!freq)
		LOG_ERROR2("ZPS");
	if(cdf+freq<cdf)
		LOG_ERROR2("Invalid CDF");
#endif
#ifdef AC3_PREC
	ec->low+=r2*cdf+((ec->range&(~0ULL>>(64-AC3_PROB_BITS)))*cdf>>AC3_PROB_BITS);
	ec->range=r2*freq+((ec->range&(~0ULL>>(64-AC3_PROB_BITS)))*freq>>AC3_PROB_BITS)-1;
#else
	ec->low+=r2*cdf;
	ec->range=r2*freq-1;//must decrement hi because decoder fails when code == hi2
#endif
	acval_enc(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
INLINE int ac3_dec(AC3 *ec, const unsigned *CDF, int nlevels)
{
	unsigned long long r2;
	unsigned cdf, freq;
	int range, sym;

	r2=ec->range>>AC3_PROB_BITS;
	AC3_RENORM_STATEMENT(!r2)
	{
		ac3_dec_renorm(ec);
		r2=ec->range>>AC3_PROB_BITS;
	}
#ifdef AC3_PREC
	{
		unsigned long long diff=ec->code-ec->low;
		cdf=(unsigned)_udiv128(diff>>(64-AC3_PROB_BITS), diff<<AC3_PROB_BITS|((1LL<<AC3_PROB_BITS)-1), ec->range, &diff);
	}
#else
	cdf=(unsigned)((ec->code-ec->low)/r2);
#endif
	sym=0;
	range=nlevels;
	while(range)
	{
		int floorhalf=range>>1;
		if(cdf>=CDF[sym+floorhalf+1])
			sym+=range-floorhalf;
		range=floorhalf;
	}
	cdf=CDF[sym];
	freq=CDF[sym+1]-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!freq)
		LOG_ERROR2("ZPS");
	if(cdf+freq<cdf)
		LOG_ERROR2("Invalid CDF");
#endif
#ifdef AC3_PREC
	ec->low+=r2*cdf+((ec->range&(~0ULL>>(64-AC3_PROB_BITS)))*cdf>>AC3_PROB_BITS);
	ec->range=r2*freq+((ec->range&(~0ULL>>(64-AC3_PROB_BITS)))*freq>>AC3_PROB_BITS)-1;
#else
	ec->low+=r2*cdf;
	ec->range=r2*freq-1;//must decrement hi because decoder fails when code == hi2
#endif
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return sym;
}
INLINE void ac3_enc_bin(AC3 *ec, int bit, unsigned p0, int probbits)
{
	unsigned long long mid;

#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!p0||p0>=1<<probbits)
		LOG_ERROR2("ZPS");
#endif
	mid=ec->range>>probbits;
	AC3_RENORM_STATEMENT(!mid)
	{
		ac3_enc_renorm(ec);
		mid=ec->range>>probbits;
	}
	mid*=p0;
#ifdef AC3_PREC
	mid+=(ec->range&((1LL<<probbits)-1))*p0>>probbits;
#endif
	if(bit)
	{
		ec->low+=mid;
		ec->range-=mid;
	}
	else
		ec->range=mid-1;
	acval_enc(bit, bit?p0:0, bit?(1<<probbits)-p0:p0, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
INLINE int ac3_dec_bin(AC3 *ec, unsigned p0, int probbits)
{
	unsigned long long mid, t2;
	int bit;
	
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!p0||p0>=1<<probbits)
		LOG_ERROR2("ZPS");
#endif
	mid=ec->range>>probbits;
	AC3_RENORM_STATEMENT(!mid)
	{
		ac3_dec_renorm(ec);
		mid=ec->range>>probbits;
	}
	mid*=p0;
#ifdef AC3_PREC
	mid+=(ec->range&((1LL<<probbits)-1))*p0>>probbits;
#endif
	t2=ec->low+mid;
	bit=ec->code>=t2;
	ec->range-=mid;
	if(bit)
		ec->low=t2;
	else
		ec->range=mid-1;
	acval_dec(bit, bit?p0:0, bit?(1<<probbits)-p0:p0, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return bit;
}


//Golomb-Rice Coder
typedef struct GolombRiceCoderStruct
{
	unsigned long long cache;
	int nbits;//enc: number of free bits in cache, dec: number of unread bits in cache
	int is_enc;//for padding
	const unsigned char *srcptr, *srcend, *srcstart;
	unsigned char *dstptr, *dstend, *dststart;
	DList *dst;
} GolombRiceCoder;
INLINE size_t gr_enc_flush(GolombRiceCoder *ec)
{
#ifdef EC_USE_ARRAY
	if(ec->dstptr+sizeof(ec->cache)>ec->dstend)//compression failed
		return 0;
	memcpy(ec->dstptr, &ec->cache, sizeof(ec->cache));
	ec->dstptr+=sizeof(ec->cache);
	return 1;
#else
	dlist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));//size is qword-aligned
	return 1;
#endif
}
INLINE int gr_dec_impl_read(GolombRiceCoder *ec)
{
	if(ec->srcptr+sizeof(ec->cache)>ec->srcend)
	{
		LOG_ERROR2("buffer overflow");
		return 1;
	}
	ec->nbits=sizeof(ec->cache)<<3;
	memcpy(&ec->cache, ec->srcptr, sizeof(ec->cache));
	ec->srcptr+=sizeof(ec->cache);
	return 0;
}
INLINE void gr_enc_init(GolombRiceCoder *ec,
#ifdef EC_USE_ARRAY
	unsigned char *start, unsigned char *end
#else
	DList *dst
#endif
)
{
	memset(ec, 0, sizeof(*ec));
	ec->cache=0;
	ec->nbits=sizeof(ec->cache)<<3;
	ec->is_enc=1;
#ifdef EC_USE_ARRAY
	ec->dstptr=start;
	ec->dstend=end;
	ec->dststart=start;
#else
	ec->dst=dst;
#endif
}
INLINE void gr_dec_init(GolombRiceCoder *ec, const unsigned char *start, const unsigned char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->cache=0;
	ec->nbits=0;
	ec->is_enc=0;
	ec->srcptr=start;
	ec->srcend=end;
	ec->srcstart=start;
}

INLINE int gr_enc(GolombRiceCoder *ec, unsigned sym, unsigned magnitude)
{
	//buffer: {c,c,c,b,b,a,a,a, f,f,f,e,e,e,d,c}, cache: MSB gg[hhh]000 LSB	nbits 6->3, code h is about to be emitted
	//written 64-bit words are byte-reversed because the CPU is little-endian

	//magnitude+=!magnitude;
	int nbypass=32-_lzcnt_u32(magnitude);
	int nzeros=sym/magnitude, bypass=sym%magnitude;
	if(nzeros>=ec->nbits)//fill the rest of cache with zeros, and flush
	{
		nzeros-=ec->nbits;
		if(!gr_enc_flush(ec))
			return 0;
		//dlist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		if(nzeros>=(int)(sizeof(ec->cache)<<3))//just flush zeros
		{
			ec->cache=0;
			do
			{
				nzeros-=(sizeof(ec->cache)<<3);
				if(!gr_enc_flush(ec))
					return 0;
				//dlist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
			}
			while(nzeros>(int)(sizeof(ec->cache)<<3));
		}
		ec->cache=0;
		ec->nbits=(sizeof(ec->cache)<<3);
	}
	//now there is room for zeros:  0 <= nzeros < nbits <= 64
	ec->nbits-=nzeros;//emit remaining zeros to cache

	int nunused=(1<<nbypass)-magnitude;//truncated binary code
	if(bypass<nunused)	//emit(bypass, nbypass-1)
		--nbypass;
	else			//emit(bypass+nunused, nbypass)
		bypass+=nunused;

	bypass|=1<<nbypass;//append 1 stop bit
	++nbypass;
	if(nbypass>=ec->nbits)//not enough free bits in cache:  fill cache, write to list, and repeat
	{
		nbypass-=ec->nbits;
		ec->cache|=(unsigned long long)bypass>>nbypass;
		bypass&=(1<<nbypass)-1;
		if(!gr_enc_flush(ec))
			return 0;
		//dlist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		ec->cache=0;
		ec->nbits=sizeof(ec->cache)<<3;
	}
	//now there is room for bypass:  0 <= nbypass < nbits <= 64
	ec->nbits-=nbypass;//emit remaining bypass to cache
	ec->cache|=(unsigned long long)bypass<<ec->nbits;
	return 1;
}
INLINE unsigned gr_dec(GolombRiceCoder *ec, unsigned magnitude)
{
	//cache: MSB 00[hhh]ijj LSB		nbits 6->3, h is about to be read (past codes must be cleared from cache)
	
	//magnitude+=!magnitude;
	int nbypass=FLOOR_LOG2(magnitude);
	int sym=0, nleadingzeros=0;
	if(!ec->nbits)//cache is empty
		goto read;
	for(;;)//cache reading loop
	{
		nleadingzeros=ec->nbits-FLOOR_LOG2_P1(ec->cache);//count leading zeros
		ec->nbits-=nleadingzeros;//remove accessed zeros
		sym+=nleadingzeros;

		if(ec->nbits)
			break;
	read://cache is empty
		if(gr_dec_impl_read(ec))
			return 0;
	}
	//now  0 < nbits <= 64
	--ec->nbits;
	//now  0 <= nbits < 64
	ec->cache-=1ULL<<ec->nbits;//remove stop bit
	//ec->cache&=(1ULL<<ec->nbits)-1;

	//nleadingzeros=(sizeof(ec->cache)<<3)-ec->nbits;//X  won't work when removed bit is LSB
	//ec->cache<<=nleadingzeros;
	//ec->cache>>=nleadingzeros;

	sym*=magnitude;
	unsigned bypass=0, nunused=(1<<(nbypass+1))-magnitude;
	if(ec->nbits<nbypass)
	{
		nbypass-=ec->nbits;
		bypass|=(int)(ec->cache<<nbypass);
		if(gr_dec_impl_read(ec))
			return 0;
	}
	if(nbypass)
	{
		ec->nbits-=nbypass;
		bypass|=(int)(ec->cache>>ec->nbits);
		ec->cache&=(1ULL<<ec->nbits)-1;
	}
	if(bypass>=nunused)
	{
		if(ec->nbits<1)
		{
			if(gr_dec_impl_read(ec))
				return 0;
		}
		--ec->nbits;
		bypass=(bypass<<1|(int)(ec->cache>>ec->nbits))-nunused;
		ec->cache&=(1ULL<<ec->nbits)-1;
	}
	sym+=bypass;
	return sym;
}

INLINE int gr_enc_POT(GolombRiceCoder *ec, int sym, int nbypass)
{
	//buffer: {c,c,c,b,b,a,a,a, f,f,f,e,e,e,d,c}, cache: MSB gg[hhh]000 LSB	nbits 6->3, code h is about to be emitted
	//written 64-bit words are byte-reversed because the CPU is little-endian

	int nzeros=sym>>nbypass, bypass=sym&((1<<nbypass)-1);
	if(nzeros>=ec->nbits)//fill the rest of cache with zeros, and flush
	{
		nzeros-=ec->nbits;
		if(!gr_enc_flush(ec))
			return 0;
		//dlist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		if(nzeros>=(int)(sizeof(ec->cache)<<3))//just flush zeros
		{
			ec->cache=0;
			do
			{
				nzeros-=(sizeof(ec->cache)<<3);
				if(!gr_enc_flush(ec))
					return 0;
				//dlist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
			}
			while(nzeros>(int)(sizeof(ec->cache)<<3));
		}
		ec->cache=0;
		ec->nbits=(sizeof(ec->cache)<<3);
	}
	//now there is room for zeros:  0 <= nzeros < nbits <= 64
	ec->nbits-=nzeros;//emit remaining zeros to cache

	bypass|=1<<nbypass;//append 1 stop bit
	++nbypass;
	if(nbypass>=ec->nbits)//not enough free bits in cache:  fill cache, write to list, and repeat
	{
		nbypass-=ec->nbits;
		ec->cache|=(unsigned long long)bypass>>nbypass;
		bypass&=(1<<nbypass)-1;
		if(!gr_enc_flush(ec))
			return 0;
		//dlist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		ec->cache=0;
		ec->nbits=sizeof(ec->cache)<<3;
	}
	//now there is room for bypass:  0 <= nbypass < nbits <= 64
	ec->nbits-=nbypass;//emit remaining bypass to cache
	ec->cache|=(unsigned long long)bypass<<ec->nbits;
	return 1;
}
INLINE unsigned gr_dec_POT(GolombRiceCoder *ec, int nbypass)
{
	//cache: MSB 00[hhh]ijj LSB		nbits 6->3, h is about to be read (past codes must be cleared from cache)
	
	int sym=0, nleadingzeros=0;
	if(!ec->nbits)//cache is empty
		goto read;
	for(;;)//cache reading loop
	{
		nleadingzeros=ec->nbits-FLOOR_LOG2_P1(ec->cache);//count leading zeros
		ec->nbits-=nleadingzeros;//remove accessed zeros
		sym+=nleadingzeros;

		if(ec->nbits)
			break;
	read://cache is empty
		if(gr_dec_impl_read(ec))
			return 0;
	}
	//now  0 < nbits <= 64
	--ec->nbits;
	//now  0 <= nbits < 64
	ec->cache-=1ULL<<ec->nbits;//remove stop bit

	unsigned bypass=0;
	sym<<=nbypass;
	if(ec->nbits<nbypass)
	{
		nbypass-=ec->nbits;
		bypass|=(int)(ec->cache<<nbypass);
		if(gr_dec_impl_read(ec))
			return 0;
	}
	if(nbypass)
	{
		ec->nbits-=nbypass;
		bypass|=(int)(ec->cache>>ec->nbits);
		ec->cache&=(1ULL<<ec->nbits)-1;
	}
	return sym|bypass;
}


#ifdef __cplusplus
}
#endif
#endif
