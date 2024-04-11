#pragma once
#ifndef INC_AC_H
#define INC_AC_H
#include"util.h"
#include<stdio.h>
#include<string.h>
#ifdef __cplusplus
extern "C"
{
#endif


//	#define AC_VALIDATE

#ifdef AC_VALIDATE
void acval_enc(int sym, int cdf, int freq, unsigned long long lo1, unsigned long long hi1, unsigned long long lo2, unsigned long long hi2, unsigned long long cache, int nbits);
void acval_dec(int sym, int cdf, int freq, unsigned long long lo1, unsigned long long hi1, unsigned long long lo2, unsigned long long hi2, unsigned long long cache, int nbits, unsigned long long code);
void acval_dump();
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
	if(lo1>>48||hi1>>48||lo2>>48||hi2>>48)
		LOG_ERROR2("");
	if(!acval)
		ARRAY_ALLOC(ACVAL, acval, 0, 0, 0, 0);
	ARRAY_APPEND(acval, &val, 1, 1, 0);
}
void acval_print(int idx, ACVAL const *p, int dec)
{
	printf("%9d sym %03X cdf %04X freq %04X range %012llX~%012llX:%012llX->%012llX~%012llX:%012llX",
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

	acval_cbuf[acval_idx%ACVAL_CBUFSIZE]=val2;

	if(acval_idx>=acval->count)
		LOG_ERROR2("AC validation index error");

	val=(ACVAL*)array_at(&acval, acval_idx);
	if(lo1>>48||hi1>>48||lo2>>48||hi2>>48)
		LOG_ERROR2("");
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

	
//	#define DEBUG_ANS

#ifdef DEBUG_ANS
typedef struct DebugANSInfoStruct
{
	unsigned s0, state, cdf, freq, id, kx, ky;
	unsigned char kq, kc;
	unsigned short sym;
} DebugANSInfo;
extern SList states;
extern int debug_channel;
void debug_enc_update(unsigned state, unsigned cdf, unsigned freq, int kx, int ky, int kq, int kc, unsigned char sym);
void debug_dec_update(unsigned state, unsigned cdf, unsigned freq, int kx, int ky, int kq, int kc, unsigned char sym);
#ifdef AC_IMPLEMENTATION
SList states={0};
int debug_channel=0;
void debug_ans_print(DebugANSInfo *info, int dec)
{
	printf("%6d state 0x%08X%s0x%08X cdf 0x%04X freq 0x%04X sym 0x%02X\n", info->id, info->s0, dec?"<-":"->", info->state, info->cdf, info->freq, info->sym);
}
void debug_enc_dump(DebugANSInfo *i2)
{
	debug_ans_print(i2, 1);
	printf("\n");
	for(int k=0;k<20&&states.count;++k)
	{
		DebugANSInfo *i0=(DebugANSInfo*)STACK_TOP(&states);
		debug_ans_print(i0, 0);
		STACK_POP(&states);
	}
}
void debug_enc_update(unsigned state, unsigned cdf, unsigned freq, int kx, int ky, int kq, int kc, unsigned char sym)
{
	if(kc==debug_channel)
	{
		unsigned s0=state;

		if(!states.count)
			slist_init(&states, sizeof(DebugANSInfo), 0);

		state=state/freq<<16|(cdf+state%freq);//enc update

		DebugANSInfo info={s0, state, cdf, freq, (unsigned)states.count, kx, ky, kq, kc, sym};
		STACK_PUSH(&states, &info);
	}
}
void debug_dec_update(unsigned state, unsigned cdf, unsigned freq, int kx, int ky, int kq, int kc, unsigned char sym)
{
	if(kc==debug_channel)
	{
		if(!states.count)
			LOG_ERROR2("Nothing to decode");
		DebugANSInfo *i0=(DebugANSInfo*)STACK_TOP(&states), info;
		memcpy(&info, i0, sizeof(info));

		unsigned s0=freq*(state>>16)+(unsigned short)state-cdf;//dec update

		if(info.s0!=s0||info.state!=state||info.cdf!=cdf||info.freq!=freq||kx!=info.kx||ky!=info.ky||kq!=info.kq||kc!=info.kc||info.sym!=sym)
		{
			DebugANSInfo i2={s0, state, cdf, freq, (unsigned)states.count-1, kx, ky, kq, kc, sym};
			debug_enc_dump(&i2);
			LOG_ERROR2("Decode error  (%d decodes remaining)", info.id);
		}

		STACK_POP(&states);
	}
}
#endif
#else
#define debug_enc_update(...)
#define debug_dec_update(...)
#endif


//arithmetic coder
#define PROB_BITS 16
typedef struct ArithmeticCoderStruct
{
	unsigned long long low, range;//FIXME range is 32-bit after initial renorm
//#ifndef DISABLE_CACHE
//	unsigned long long cache;//cache is read/written MSB->LSB
//	int nbits;//enc: number of free bits in cache,  dec: number of unread bits in cache
//	//int cidx;//number of renorms remaining in cache before it gets updated
//	int is_enc;
//#endif
	ArrayHandle *arr;
	DList *list;
	const unsigned char *srcptr, *srcend;
	unsigned long long code;
} ArithmeticCoder;
INLINE void ac_enc_init(ArithmeticCoder *ec,
#ifdef EC_USE_ARRAY
	ArrayHandle *arr
#else
	DList *list
#endif
)
{
	memset(ec, 0, sizeof(*ec));
	ec->low=0;
	//ec->range=~0LLU>>PROB_BITS;
	ec->range=~0LLU>>(PROB_BITS<<1);
//#ifndef DISABLE_CACHE
//	ec->cache=0;
//	ec->nbits=sizeof(ec->low)<<3;
//	//ec->cidx=4;
//	//ec->is_enc=1;
//#endif
#ifdef EC_USE_ARRAY
	array_append(arr, 0, 1, 0, 0, 0, 0);
	ec->arr=arr;
#else
	ec->list=list;
#endif
}
INLINE void ac_dec_init(ArithmeticCoder *ec, const unsigned char *start, unsigned const char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->low=0;
	//ec->range=~0LLU>>PROB_BITS;
	ec->range=~0LLU>>(PROB_BITS<<1);
//#ifndef DISABLE_CACHE
//	ec->cache=0;
//	ec->nbits=PROB_BITS;
//	//ec->cidx=1;
//	//ec->is_enc=0;
//#endif
	ec->srcptr=start;
	ec->srcend=end;
	
//#ifndef DISABLE_CACHE
//	if(ec->srcptr+8>ec->srcend)
//	{
//		LOG_ERROR2("buffer overflow");
//		return;
//	}
//	memcpy(&ec->cache, ec->srcptr, 8);
//	ec->srcptr+=8;
//	ec->code=ec->cache>>PROB_BITS;//leave PROB_BITS bits in cache
//#else
	if(ec->srcptr+6>ec->srcend)
	{
		LOG_ERROR2("buffer overflow");
		return;
	}
	memcpy((unsigned char*)&ec->code+4, ec->srcptr+0, 2);
	memcpy((unsigned char*)&ec->code+2, ec->srcptr+2, 2);
	memcpy((unsigned char*)&ec->code+0, ec->srcptr+4, 2);
	ec->srcptr+=6;
//#endif
}
INLINE void ac_enc_renorm(ArithmeticCoder *ec)//fast renorm by F. Rubin 1979
{
//#ifndef DISABLE_CACHE
//	if(ec->nbits<PROB_BITS)
//	//if(ec->cidx)
//	//	--ec->cidx;
//	//else
//	{
//#ifdef EC_USE_ARRAY
//		array_append(ec->arr, &ec->cache, 1, 8, 1, 0, 0);
//#else
//		dlist_push_back(ec->list, &ec->cache, 8);
//#endif
//		ec->nbits=sizeof(ec->low)<<3;
//		//ec->cidx=3;
//		ec->cache=0;
//	}
//	//((unsigned short*)&ec->cache)[ec->cidx]=((unsigned short*)&ec->low)[2];
//	//((unsigned short*)&ec->range)[1]=((unsigned short*)&ec->range)[0];
//	//((unsigned short*)&ec->low)[2]=((unsigned short*)&ec->low)[1];
//	//((unsigned short*)&ec->range)[0]=~0;
//	//((unsigned short*)&ec->low)[1]=0;
//#if 1
//	ec->nbits-=PROB_BITS;
//	ec->range<<=PROB_BITS;
//	ec->cache|=(ec->low>>((sizeof(ec->low)<<3)-PROB_BITS*2))<<ec->nbits;//append top PROB_BITS bits of 'low' to cache
//	ec->low<<=PROB_BITS;
//	ec->range|=(1LL<<PROB_BITS)-1;
//	//ec->range&=~0LLU>>PROB_BITS;//X
//	ec->low&=~0LLU>>PROB_BITS;
//#endif
//#else
#ifdef EC_USE_ARRAY
	array_append(ec->arr, (unsigned char*)&ec->low+4, 1, 2, 1, 0, 0);
#else
	dlist_push_back(ec->list, (unsigned char*)&ec->low+4, 2);
#endif
	ec->range<<=PROB_BITS;
	ec->low<<=PROB_BITS;
	ec->range|=(1LL<<PROB_BITS)-1;
	ec->low&=~0LLU>>PROB_BITS;
//#endif

	if(ec->low+ec->range>(~0LLU>>PROB_BITS))//clamp hi to register size after renorm
		ec->range=(~0LLU>>PROB_BITS)-ec->low;
}
INLINE void ac_dec_renorm(ArithmeticCoder *ec)//fast renorm by F. Rubin 1979
{
//#ifndef DISABLE_CACHE
//	if(ec->nbits<PROB_BITS)
//	//if(ec->cidx)
//	//	--ec->cidx;
//	//else
//	{
//		if(ec->srcptr+8>ec->srcend)
//		{
//#ifdef AC_VALIDATE
//			printf("buffer overflow\n");
//			acval_dump();
//#endif
//			LOG_ERROR2("buffer overflow");
//			return;
//		}
//		ec->nbits=sizeof(ec->low)<<3;
//		//ec->cidx=3;
//		memcpy(&ec->cache, ec->srcptr, 8);
//		ec->srcptr+=8;
//	}
//	//((unsigned short*)&ec->range)[1]=((unsigned short*)&ec->range)[0];
//	//((unsigned short*)&ec->low)[2]=((unsigned short*)&ec->low)[1];
//	//((unsigned short*)&ec->range)[0]=~0;
//	//((unsigned short*)&ec->code)[0]=((unsigned short*)&ec->cache)[ec->cidx];
//	//((unsigned short*)&ec->low)[1]=0;
//#if 1
//	ec->nbits-=PROB_BITS;
//	ec->range<<=PROB_BITS;
//	ec->low<<=PROB_BITS;
//	ec->code<<=PROB_BITS;
//	ec->range|=(1LL<<PROB_BITS)-1;
//	ec->code|=ec->cache>>ec->nbits&((1LL<<PROB_BITS)-1);
//	ec->low&=~0LLU>>PROB_BITS;
//	//ec->range&=~0LLU>>PROB_BITS;//X
//	ec->code&=~0LLU>>PROB_BITS;
//#endif
//#else
	if(ec->srcptr+2>ec->srcend)
	{
#ifdef AC_VALIDATE
		printf("buffer overflow\n");
		acval_dump();
#endif
		LOG_ERROR2("buffer overflow");
		return;
	}
	ec->code<<=PROB_BITS;
	ec->range<<=PROB_BITS;
	ec->low<<=PROB_BITS;
	ec->range|=(1LL<<PROB_BITS)-1;
	memcpy(&ec->code, ec->srcptr, 2);
	ec->srcptr+=2;
	ec->low&=~0LLU>>PROB_BITS;
	ec->code&=~0LLU>>PROB_BITS;
//#endif

	if(ec->low+ec->range>(~0LLU>>PROB_BITS))//clamp hi to register size after renorm
		ec->range=(~0LLU>>PROB_BITS)-ec->low;
}
INLINE void ac_enc_flush(ArithmeticCoder *ec)
{
//#ifndef DISABLE_CACHE
//	ac_enc_renorm(ec);
//	ac_enc_renorm(ec);
//	ac_enc_renorm(ec);
//#ifdef EC_USE_ARRAY
//	array_append(ec->arr, &ec->cache, 1, 8, 1, 0, 0);
//#else
//	dlist_push_back(ec->list, &ec->cache, 8);
//#endif
//#else
#ifdef EC_USE_ARRAY
	array_append(ec->arr, (unsigned char*)&ec->low+4, 1, 2, 1, 0, 0);
	array_append(ec->arr, (unsigned char*)&ec->low+2, 1, 2, 1, 0, 0);
	array_append(ec->arr, (unsigned char*)&ec->low+0, 1, 2, 1, 0, 0);
#else
	dlist_push_back(ec->list, (unsigned char*)&ec->low+4, 2);
	dlist_push_back(ec->list, (unsigned char*)&ec->low+2, 2);
	dlist_push_back(ec->list, (unsigned char*)&ec->low+0, 2);
#endif
//#endif
}

INLINE void ac_enc(ArithmeticCoder *ec, int sym, const unsigned *CDF)//CDF is 16 bit
{
	unsigned cdf=CDF[sym];
	int freq=CDF[sym+1]-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(freq<=0)
		LOG_ERROR2("ZPS");
	if(cdf>0x10000||cdf+freq>0x10000)
		LOG_ERROR2("");
	if(ec->range==~0LLU)
		LOG_ERROR2("");
#endif
	ec->low+=ec->range*cdf>>16;
	ec->range*=freq;
	ec->range>>=16;
	--ec->range;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))//only when freq=1 -> range=0, this loop runs twice
		ac_enc_renorm(ec);
	acval_enc(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);//
	//acval_enc(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, ec->cache, ec->cidx);//
}
INLINE int ac_dec(ArithmeticCoder *ec, const unsigned *CDF, int nlevels)
{
	unsigned cdf;
	int freq;
	int sym=0;
	
#ifdef LINEAR_SEARCH
	unsigned long long range=ec->code-ec->low;
	//range=((range<<16)+(ec->range>>1))/ec->range;
	//for(;sym<nlevels&&CDF[sym]<=range;++sym);
	sym=1;
	for(;sym<nlevels&&(ec->range*CDF[sym]>>16)<=range;++sym);
	--sym;
#else
	int L=0, R=nlevels;
	unsigned long long range=ec->code-ec->low;
	while(R)//binary search		lg(nlevels) multiplications per symbol
	{
		int floorhalf=R>>1;
		sym=L+floorhalf;
		L+=(R-floorhalf)&-(ec->range*CDF[sym]>>16<=range);
		R=floorhalf;
	}
#endif

	cdf=CDF[sym];
	freq=CDF[sym+1]-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(freq<=0||cdf>0x10000||cdf+freq>0x10000)
		LOG_ERROR2("ZPS");
#endif
	ec->low+=ec->range*cdf>>16;
	ec->range*=freq;
	ec->range>>=16;
	--ec->range;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))
		ac_dec_renorm(ec);
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);//
	//acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, ec->cache, ec->cidx, ec->code);//
	return sym;
}
INLINE int ac_dec_POT(ArithmeticCoder *ec, const unsigned *CDF, int nbits)
{
	unsigned cdf;
	int freq;
	int sym=0;
	
	unsigned long long range=ec->code-ec->low;
	switch(nbits)
	{
	default:
		LOG_ERROR2("Unsupported bit depth");
		break;
	case 9:sym|=(ec->range*CDF[sym|256]>>16<=range)<<8;
	case 8:sym|=(ec->range*CDF[sym|128]>>16<=range)<<7;
	case 7:sym|=(ec->range*CDF[sym| 64]>>16<=range)<<6;
	case 6:sym|=(ec->range*CDF[sym| 32]>>16<=range)<<5;
	case 5:sym|=(ec->range*CDF[sym| 16]>>16<=range)<<4;
	case 4:sym|=(ec->range*CDF[sym|  8]>>16<=range)<<3;
	case 3:sym|=(ec->range*CDF[sym|  4]>>16<=range)<<2;
	case 2:sym|=(ec->range*CDF[sym|  2]>>16<=range)<<1;
	case 1:sym|=(ec->range*CDF[sym|  1]>>16<=range);
		break;
	}
	//for(int mask=nlevels>>1;mask;mask>>=1)
	//	sym|=mask&-(ec->range*CDF[sym|mask]>>16<=range);

	cdf=CDF[sym];
	freq=CDF[sym+1]-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(freq<=0||cdf>0x10000||cdf+freq>0x10000)
		LOG_ERROR2("ZPS");
#endif
	ec->low+=ec->range*cdf>>16;
	ec->range*=freq;
	ec->range>>=16;
	--ec->range;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))
		ac_dec_renorm(ec);
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);//
	//acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, ec->cache, ec->cidx, ec->code);//
	return sym;
}
#ifdef DEDICATED_DECODER
INLINE int ac_dec_4bit(ArithmeticCoder *ec, const unsigned *CDF)
{
	unsigned cdf;
	int freq;
	unsigned short sym=0;
	
	unsigned long long range=ec->code-ec->low;
	sym|=8&-(ec->range*CDF[sym|8]>>16<=range);
	sym|=4&-(ec->range*CDF[sym|4]>>16<=range);
	sym|=2&-(ec->range*CDF[sym|2]>>16<=range);
	sym|=1&-(ec->range*CDF[sym|1]>>16<=range);
#if 0
	int L=0, R=16;
	int floorhalf=8;
	sym=8;
	L+=8&-(ec->range*CDF[8]>>16<=range);
	R=8;

	floorhalf=4;
	sym=L+4;
	L+=4&-(ec->range*CDF[sym]>>16<=range);
	R=4;

	floorhalf=2;
	sym=L+2;
	L+=2&-(ec->range*CDF[sym]>>16<=range);
	R=2;

	floorhalf=1;
	sym=L+1;
	L+=1&-(ec->range*CDF[sym]>>16<=range);
	R=1;

	floorhalf=0;
	sym=L+0;
	L+=1&-(ec->range*CDF[sym]>>16<=range);
	R=0;
#endif
#if 0
	__m256i permutation=_mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);
	__m256i srl32=_mm256_set_epi8(
		-1, -1, -1, -1, 15, 14, 13, 12, -1, -1, -1, -1,  7,  6,  5,  4,
		-1, -1, -1, -1, 15, 14, 13, 12, -1, -1, -1, -1,  7,  6,  5,  4
	);
	__m256i srl16=_mm256_set_epi8(
		-1, -1, 15, 14, 13, 12, 11, 10, -1, -1,  7,  6,  5,  4,  3,  2,
		-1, -1, 15, 14, 13, 12, 11, 10, -1, -1,  7,  6,  5,  4,  3,  2
	);
	__m256i range0=_mm256_set1_epi64x(ec->range);
	__m256i range1=_mm256_set1_epi64x(ec->code-ec->low);
	__m256i c0=_mm256_loadu_si256((__m256i const*)(CDF+0));
	__m256i c1=_mm256_loadu_si256((__m256i const*)(CDF+8));
	c0=_mm256_permutevar8x32_epi32(c0, permutation);
	c1=_mm256_permutevar8x32_epi32(c1, permutation);
	//c0=_mm256_permute4x64_epi64(c0, _MM_SHUFFLE(3, 1, 2, 0));//slower but across lanes
	//c1=_mm256_permute4x64_epi64(c1, _MM_SHUFFLE(3, 1, 2, 0));
	__m256i p0=_mm256_mul_epu32(range0, c0);
	__m256i p2=_mm256_mul_epu32(range0, c1);
	c0=_mm256_shuffle_epi8(c0, srl32);
	c1=_mm256_shuffle_epi8(c1, srl32);
	__m256i p1=_mm256_mul_epu32(range0, c0);
	__m256i p3=_mm256_mul_epu32(range0, c1);
	p0=_mm256_shuffle_epi8(p0, srl16);
	p1=_mm256_shuffle_epi8(p1, srl16);
	p2=_mm256_shuffle_epi8(p2, srl16);
	p3=_mm256_shuffle_epi8(p3, srl16);
	p0=_mm256_cmpgt_epi64(p0, range1);
	p1=_mm256_cmpgt_epi64(p1, range1);
	p2=_mm256_cmpgt_epi64(p2, range1);
	p3=_mm256_cmpgt_epi64(p3, range1);
	sym|=_mm256_movemask_pd(_mm256_castsi256_pd(p3)), sym<<=4;
	sym|=_mm256_movemask_pd(_mm256_castsi256_pd(p2)), sym<<=4;
	sym|=_mm256_movemask_pd(_mm256_castsi256_pd(p1)), sym<<=4;
	sym|=_mm256_movemask_pd(_mm256_castsi256_pd(p0));
	sym-=sym>>1&0x5555;//sym = Hamming weight of 16-bit mask
	sym=(sym&0x3333)+(sym>>2&0x3333);
	sym*=0x1111;
	sym>>=12;
	sym=15-sym;
#endif

	cdf=CDF[sym];
	freq=CDF[sym+1]-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(freq<=0||cdf>0x10000||cdf+freq>0x10000)
		LOG_ERROR2("ZPS");
#endif
	ec->low+=ec->range*cdf>>16;
	ec->range*=freq;
	ec->range>>=16;
	--ec->range;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))
		ac_dec_renorm(ec);
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);//
	//acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, ec->cache, ec->cidx, ec->code);//
	return sym;
}
INLINE int ac_dec_5bit(ArithmeticCoder *ec, const unsigned *CDF)
{
	unsigned cdf;
	int freq;
	unsigned sym=0;
	
	unsigned long long range=ec->code-ec->low;
	sym|=16&-(ec->range*CDF[sym|16]>>16<=range);
	sym|= 8&-(ec->range*CDF[sym| 8]>>16<=range);
	sym|= 4&-(ec->range*CDF[sym| 4]>>16<=range);
	sym|= 2&-(ec->range*CDF[sym| 2]>>16<=range);
	sym|= 1&-(ec->range*CDF[sym| 1]>>16<=range);
#if 0
	__m256i permutation=_mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);
	//__m256i swapmid32=_mm256_set_epi8(
	//	15, 14, 13, 12,		7,  6,  5,  4,		11, 10,  9,  8,		3,  2,  1,  0,
	//	15, 14, 13, 12,		7,  6,  5,  4,		11, 10,  9,  8,		3,  2,  1,  0
	//);
	__m256i srl32=_mm256_set_epi8(
		-1, -1, -1, -1, 15, 14, 13, 12, -1, -1, -1, -1,  7,  6,  5,  4,
		-1, -1, -1, -1, 15, 14, 13, 12, -1, -1, -1, -1,  7,  6,  5,  4
	);
	__m256i srl16=_mm256_set_epi8(
		-1, -1, 15, 14, 13, 12, 11, 10, -1, -1,  7,  6,  5,  4,  3,  2,
		-1, -1, 15, 14, 13, 12, 11, 10, -1, -1,  7,  6,  5,  4,  3,  2
	);
	__m256i range0=_mm256_set1_epi64x(ec->range);
	__m256i range1=_mm256_set1_epi64x(ec->code-ec->low);
	__m256i c0=_mm256_loadu_si256((__m256i const*)(CDF+ 0));
	__m256i c1=_mm256_loadu_si256((__m256i const*)(CDF+ 8));
	__m256i c2=_mm256_loadu_si256((__m256i const*)(CDF+16));
	__m256i c3=_mm256_loadu_si256((__m256i const*)(CDF+24));
	c0=_mm256_permutevar8x32_epi32(c0, permutation);
	c1=_mm256_permutevar8x32_epi32(c1, permutation);
	c2=_mm256_permutevar8x32_epi32(c2, permutation);
	c3=_mm256_permutevar8x32_epi32(c3, permutation);
	//c0=_mm256_shuffle_epi8(c0, swapmid32);//X
	//c1=_mm256_shuffle_epi8(c1, swapmid32);
	//c2=_mm256_shuffle_epi8(c2, swapmid32);
	//c3=_mm256_shuffle_epi8(c3, swapmid32);
	//c0=_mm256_permute4x64_epi64(c0, _MM_SHUFFLE(3, 1, 2, 0));//slower but across lanes	X
	//c1=_mm256_permute4x64_epi64(c1, _MM_SHUFFLE(3, 1, 2, 0));
	//c2=_mm256_permute4x64_epi64(c2, _MM_SHUFFLE(3, 1, 2, 0));
	//c3=_mm256_permute4x64_epi64(c3, _MM_SHUFFLE(3, 1, 2, 0));
	__m256i p0=_mm256_mul_epu32(range0, c0);
	__m256i p2=_mm256_mul_epu32(range0, c1);
	__m256i p4=_mm256_mul_epu32(range0, c2);
	__m256i p6=_mm256_mul_epu32(range0, c3);
	c0=_mm256_shuffle_epi8(c0, srl32);
	c1=_mm256_shuffle_epi8(c1, srl32);
	c2=_mm256_shuffle_epi8(c2, srl32);
	c3=_mm256_shuffle_epi8(c3, srl32);
	__m256i p1=_mm256_mul_epu32(range0, c0);
	__m256i p3=_mm256_mul_epu32(range0, c1);
	__m256i p5=_mm256_mul_epu32(range0, c2);
	__m256i p7=_mm256_mul_epu32(range0, c3);
	p0=_mm256_shuffle_epi8(p0, srl16);
	p1=_mm256_shuffle_epi8(p1, srl16);
	p2=_mm256_shuffle_epi8(p2, srl16);
	p3=_mm256_shuffle_epi8(p3, srl16);
	p4=_mm256_shuffle_epi8(p4, srl16);
	p5=_mm256_shuffle_epi8(p5, srl16);
	p6=_mm256_shuffle_epi8(p6, srl16);
	p7=_mm256_shuffle_epi8(p7, srl16);
	p0=_mm256_cmpgt_epi64(p0, range1);
	p1=_mm256_cmpgt_epi64(p1, range1);
	p2=_mm256_cmpgt_epi64(p2, range1);
	p3=_mm256_cmpgt_epi64(p3, range1);
	p4=_mm256_cmpgt_epi64(p4, range1);
	p5=_mm256_cmpgt_epi64(p5, range1);
	p6=_mm256_cmpgt_epi64(p6, range1);
	p7=_mm256_cmpgt_epi64(p7, range1);
	sym|=_mm256_movemask_pd(_mm256_castsi256_pd(p7)), sym<<=4;
	sym|=_mm256_movemask_pd(_mm256_castsi256_pd(p6)), sym<<=4;
	sym|=_mm256_movemask_pd(_mm256_castsi256_pd(p5)), sym<<=4;
	sym|=_mm256_movemask_pd(_mm256_castsi256_pd(p4)), sym<<=4;
	sym|=_mm256_movemask_pd(_mm256_castsi256_pd(p3)), sym<<=4;
	sym|=_mm256_movemask_pd(_mm256_castsi256_pd(p2)), sym<<=4;
	sym|=_mm256_movemask_pd(_mm256_castsi256_pd(p1)), sym<<=4;
	sym|=_mm256_movemask_pd(_mm256_castsi256_pd(p0));
	sym-=sym>>1&0x55555555;//sym = Hamming weight of 32-bit mask
	sym=(sym&0x33333333)+(sym>>2&0x33333333);
	sym=(sym&0x0F0F0F0F)+(sym>>4&0x0F0F0F0F);
	sym*=0x01010101;
	sym>>=24;
	sym=31-sym;
#endif
	cdf=CDF[sym];
	freq=CDF[sym+1]-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(freq<=0||cdf>0x10000||cdf+freq>0x10000)
		LOG_ERROR2("ZPS");
#endif
	ec->low+=ec->range*cdf>>16;
	ec->range*=freq;
	ec->range>>=16;
	--ec->range;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))
		ac_dec_renorm(ec);
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);//
	//acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, ec->cache, ec->cidx, ec->code);//
	return sym;
}
#endif

INLINE void ac_enc_p16(ArithmeticCoder *ec, int sym, const unsigned short *CDF, int nlevels)//CDF is 16 bit
{
	unsigned cdf=CDF[sym];
	int freq=(sym>=nlevels-1?0x8000:CDF[sym+1])-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(freq<=0)
		LOG_ERROR2("ZPS");
	if(cdf>0x10000||cdf+freq>0x10000)
		LOG_ERROR2("");
	if(ec->range==~0LLU)
		LOG_ERROR2("");
#endif
	//if(cdf>>15||freq>>15)
	//	LOG_ERROR2("LOL_1");
	ec->low+=ec->range*cdf>>16;//FIXME range is 32-bit, use __umul()
	ec->range=ec->range*freq>>16;
	--ec->range;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))//only when freq=1 -> range=0, this loop runs twice
		ac_enc_renorm(ec);
	acval_enc(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
INLINE int ac_dec_p16_POT(ArithmeticCoder *ec, const unsigned short *CDF, int nbits)
{
	unsigned cdf;
	int freq;
	int sym=0;
	
	unsigned long long range=ec->code-ec->low;
	switch(nbits)
	{
	default:
		LOG_ERROR2("Unsupported bit depth");
		break;
	case 9:sym|=(ec->range*CDF[sym|256]>>16<=range)<<8;
	case 8:sym|=(ec->range*CDF[sym|128]>>16<=range)<<7;
	case 7:sym|=(ec->range*CDF[sym| 64]>>16<=range)<<6;
	case 6:sym|=(ec->range*CDF[sym| 32]>>16<=range)<<5;
	case 5:sym|=(ec->range*CDF[sym| 16]>>16<=range)<<4;
	case 4:sym|=(ec->range*CDF[sym|  8]>>16<=range)<<3;
	case 3:sym|=(ec->range*CDF[sym|  4]>>16<=range)<<2;
	case 2:sym|=(ec->range*CDF[sym|  2]>>16<=range)<<1;
	case 1:sym|=(ec->range*CDF[sym|  1]>>16<=range);
		break;
	}

	cdf=CDF[sym];
	freq=(sym>=(1<<nbits)-1?0x8000:CDF[sym+1])-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(freq<=0||cdf>0x10000||cdf+freq>0x10000)
		LOG_ERROR2("ZPS");
#endif
	ec->low+=ec->range*cdf>>16;
	ec->range=ec->range*freq>>16;
	--ec->range;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))
		ac_dec_renorm(ec);
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return sym;
}

INLINE void ac_enc_bypass(ArithmeticCoder *ec, int sym, int nlevels)//CDF is 16 bit
{
	unsigned cdf=(sym<<16)/nlevels;
	int freq=((sym+1)<<16)/nlevels-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(freq<=0)
		LOG_ERROR2("ZPS");
#endif
	ec->low+=ec->range*cdf>>16;
	ec->range*=freq;
	ec->range>>=16;
	--ec->range;//must decrement hi because decoder fails when code == hi2
	if(ec->range<(1LL<<PROB_BITS))
		ac_enc_renorm(ec);
	acval_enc(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);//
	//acval_enc(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, ec->cache, ec->cidx);//
}
INLINE int ac_dec_bypass(ArithmeticCoder *ec, int nlevels)
{
	unsigned cdf;
	int freq;
	int sym;
	
	cdf=(int)(((ec->code-ec->low)<<16)/ec->range);
	sym=(cdf*nlevels>>16)+1;//this is to handle the case when code == lo2
	cdf=(sym<<16)/nlevels;
	if(ec->low+(ec->range*cdf>>16)>ec->code)
	{
		--sym;
		cdf=(sym<<16)/nlevels;
	}

	freq=((sym+1)<<16)/nlevels-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(freq<=0)
		LOG_ERROR2("ZPS");
#endif
	ec->low+=ec->range*cdf>>16;
	ec->range*=freq;
	ec->range>>=16;
	--ec->range;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))
		ac_dec_renorm(ec);
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);//
	//acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, ec->cache, ec->cidx, ec->code);//
	return sym;
}

INLINE void ac_enc_bin(ArithmeticCoder *ec, unsigned short p0, int bit)
{
	unsigned long long r2=ec->range*p0>>16;
#ifdef AC_VALIDATE
	if(!p0)//reject degenerate distribution
		LOG_ERROR2("ZPS");
#endif
	if(bit)
	{
		ec->low+=r2;
		ec->range-=r2;
	}
	else
		ec->range=r2-1;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))
		ac_enc_renorm(ec);
	acval_enc(bit, bit?p0:0, bit?0x10000-p0:p0, ec->low, ec->low+ec->range, bit?ec->low+r2:ec->low, bit?ec->low+ec->range:ec->low+r2-1, 0, 0);
	//acval_enc(bit, bit?p0:0, bit?0x10000-p0:p0, ec->low, ec->low+ec->range, bit?ec->low+r2:ec->low, bit?ec->low+ec->range:ec->low+r2-1, ec->cache, ec->cidx);
}
INLINE int ac_dec_bin(ArithmeticCoder *ec, unsigned short p0)//binary AC decoder doesn't do binary search
{
	unsigned long long r2=ec->range*p0>>16;
#ifdef AC_VALIDATE
	if(!p0)//reject degenerate distribution
		LOG_ERROR2("ZPS");
#endif
	int bit=ec->code>=ec->low+r2;
	if(bit)
	{
		ec->low+=r2;
		ec->range-=r2;
	}
	else
		ec->range=r2-1;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))
		ac_dec_renorm(ec);
	acval_dec(bit, bit?p0:0, bit?0x10000-p0:p0, ec->low, ec->low+ec->range, bit?ec->low+r2:ec->low, bit?ec->low+ec->range:ec->low+r2-1, 0, 0, ec->code);
	//acval_dec(bit, bit?p0:0, bit?0x10000-p0:p0, ec->low, ec->low+ec->range, bit?ec->low+r2:ec->low, bit?ec->low+ec->range:ec->low+r2-1, ec->cache, ec->cidx, ec->code);
	return bit;
}


//asymmetric numeral systems coder
typedef struct ANSCoderStruct
{
	unsigned state;
	ArrayHandle *arr;
	DList *list;
	const unsigned char *srcptr, *srcstart;
} ANSCoder;
INLINE void ans_enc_init(ANSCoder *ec,
#ifdef EC_USE_ARRAY
	ArrayHandle *arr
#else
	DList *list
#endif
)
{
	memset(ec, 0, sizeof(*ec));
	ec->state=0x10000;
#ifdef EC_USE_ARRAY
	ec->arr=arr;
#else
	ec->list=list;
#endif
}
INLINE void ans_dec_init(ANSCoder *ec, const unsigned char *start, const unsigned char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->srcptr=end;
	ec->srcstart=start;

	ec->srcptr-=4;
	if(ec->srcptr<ec->srcstart)
		LOG_ERROR2("ANS buffer overflow");
	memcpy(&ec->state, ec->srcptr, 4);
}
INLINE void ans_enc_flush(ANSCoder *ec)
{
#ifdef EC_USE_ARRAY
	ARRAY_APPEND(*ec->arr, &ec->state, 4, 1, 0);
#else
	dlist_push_back(ec->list, &ec->state, 4);
#endif
}

INLINE void ans_enc(ANSCoder *ec, int sym, const unsigned *CDF, int nlevels)
{
	int cdf, freq;
	if(CDF)
		cdf=CDF[sym], freq=CDF[sym+1]-cdf;
	else//bypass
		cdf=(sym<<16)/nlevels, freq=((sym+1)<<16)/nlevels-cdf;
	if(!freq)
		LOG_ERROR2("ZPS");
	if((ec->state>>16)>=(unsigned)freq)//renorm
	{
#ifdef EC_USE_ARRAY
		ARRAY_APPEND(*ec->arr, &ec->state, 2, 1, 0);
#else
		dlist_push_back(ec->list, &ec->state, 2);
#endif
		ec->state>>=16;
	}
	debug_enc_update(ec->state, cdf, freq, 0, 0, 0, 0, sym);
	ec->state=ec->state/freq<<16|(cdf+ec->state%freq);//update
}
INLINE int ans_dec(ANSCoder *ec, const unsigned *CDF, int nlevels)
{
	unsigned c=(unsigned short)ec->state;
	int sym=0;

	unsigned cdf, freq;
	if(CDF)
	{
		int L=0, R=nlevels;
		while(R)//binary search		lg(nlevels) multiplications per symbol
		{
			int floorhalf=R>>1;
			sym=L+floorhalf;
			L+=(R-floorhalf)&-(CDF[sym]<=c);
			R=floorhalf;
		}
		//int L=0, R=nlevels, found=0;
		//while(L<=R)
		//{
		//	sym=(L+R)>>1;
		//	if(CDF[sym]<c)
		//		L=sym+1;
		//	else if(CDF[sym]>c)
		//		R=sym-1;
		//	else
		//	{
		//		found=1;
		//		break;
		//	}
		//}
		//if(found)
		//	for(;sym<nlevels-1&&CDF[sym+1]==c;++sym);
		//else
		//	sym=L+(L<nlevels&&CDF[L]<c)-(L!=0);

		cdf=CDF[sym], freq=CDF[sym+1]-cdf;
	}
	else//bypass
	{
		sym=c*nlevels>>16;
		cdf=(sym<<16)/nlevels, freq=((sym+1)<<16)/nlevels-cdf;
	}
	if(!freq)
		LOG_ERROR2("ZPS");
	debug_dec_update(ec->state, cdf, freq, 0, 0, 0, 0, sym);
	ec->state=freq*(ec->state>>16)+c-cdf;//update
	if(ec->state<0x10000)//renorm
	{
		ec->state<<=16;
		if(ec->srcptr-2>=ec->srcstart)
		{
			ec->srcptr-=2;
			memcpy(&ec->state, ec->srcptr, 2);
		}
	}
	return sym;
}
INLINE int ans_dec_POT(ANSCoder *ec, const unsigned *CDF, int nbits)
{
	unsigned c=(unsigned short)ec->state;
	int sym=0;

	unsigned cdf, freq;
	switch(nbits)
	{
	default:
		LOG_ERROR2("Unsupported bit depth");
		break;
	case 9:sym|=(CDF[sym|256]<=c)<<8;
	case 8:sym|=(CDF[sym|128]<=c)<<7;
	case 7:sym|=(CDF[sym| 64]<=c)<<6;
	case 6:sym|=(CDF[sym| 32]<=c)<<5;
	case 5:sym|=(CDF[sym| 16]<=c)<<4;
	case 4:sym|=(CDF[sym|  8]<=c)<<3;
	case 3:sym|=(CDF[sym|  4]<=c)<<2;
	case 2:sym|=(CDF[sym|  2]<=c)<<1;
	case 1:sym|=(CDF[sym|  1]<=c);
		break;
	}

	cdf=CDF[sym], freq=CDF[sym+1]-cdf;
	if(!freq)
		LOG_ERROR2("ZPS");
	debug_dec_update(ec->state, cdf, freq, 0, 0, 0, 0, sym);
	ec->state=freq*(ec->state>>16)+c-cdf;//update
	if(ec->state<0x10000)//renorm
	{
		ec->state<<=16;
		if(ec->srcptr-2>=ec->srcstart)
		{
			ec->srcptr-=2;
			memcpy(&ec->state, ec->srcptr, 2);
		}
	}
	return sym;
}

INLINE void ans_enc15(ANSCoder *ec, int sym, const unsigned short *CDF, int nlevels)
{
	int cdf, freq;
	if(CDF)
		cdf=CDF[sym]<<1, freq=(CDF[sym+1]<<1)-cdf;
	else//bypass
		cdf=(sym<<16)/nlevels, freq=((sym+1)<<16)/nlevels-cdf;
	if(!freq)
		LOG_ERROR2("ZPS");
	if((ec->state>>16)>=(unsigned)freq)//renorm
	{
#ifdef EC_USE_ARRAY
		ARRAY_APPEND(*ec->arr, &ec->state, 2, 1, 0);
#else
		dlist_push_back(ec->list, &ec->state, 2);
#endif
		ec->state>>=16;
	}
	debug_enc_update(ec->state, cdf, freq, 0, 0, 0, 0, sym);
	ec->state=ec->state/freq<<16|(cdf+ec->state%freq);//update
}
INLINE int ans_dec15(ANSCoder *ec, const unsigned short *CDF, int nlevels)
{
	unsigned c=(unsigned short)ec->state;
	int sym=0;

	unsigned cdf, freq;
	if(CDF)
	{
		int L=0, R=nlevels, found=0;
		while(L<=R)
		{
			sym=(L+R)>>1;
			if(CDF[sym]<c)
				L=sym+1;
			else if(CDF[sym]>c)
				R=sym-1;
			else
			{
				found=1;
				break;
			}
		}
		if(found)
			for(;sym<nlevels-1&&CDF[sym+1]==c;++sym);
		else
			sym=L+(L<nlevels&&CDF[L]<c)-(L!=0);

		cdf=CDF[sym]<<1, freq=(CDF[sym+1]<<1)-cdf;
	}
	else//bypass
	{
		sym=c*nlevels>>16;
		cdf=(sym<<16)/nlevels, freq=((sym+1)<<16)/nlevels-cdf;
	}
	if(!freq)
		LOG_ERROR2("ZPS");
	debug_dec_update(ec->state, cdf, freq, 0, 0, 0, 0, sym);
	ec->state=freq*(ec->state>>16)+c-cdf;//update
	if(ec->state<0x10000)//renorm
	{
		ec->state<<=16;
		if(ec->srcptr-2>=ec->srcstart)
		{
			ec->srcptr-=2;
			memcpy(&ec->state, ec->srcptr, 2);
		}
	}
	return sym;
}

INLINE void ans_enc_bin(ANSCoder *ec, unsigned short p0, int bit)
{
	int cdf=bit?p0:0, freq=bit?0x10000-p0:p0;
	if((ec->state>>16)>=(unsigned)freq)//renorm
	{
#ifdef EC_USE_ARRAY
		ARRAY_APPEND(*ec->arr, &ec->state, 2, 1, 0);
#else
		dlist_push_back(ec->list, &ec->state, 2);
#endif
		ec->state>>=16;
	}
	debug_enc_update(ec->state, cdf, freq, 0, 0, 0, 0, bit);
	ec->state=ec->state/freq<<16|(cdf+ec->state%freq);//update
}
INLINE int ans_dec_bin(ANSCoder *ec, unsigned short p0)
{
	unsigned c=(unsigned short)ec->state;
	int bit=c>=p0;
	
	int cdf=bit?p0:0, freq=bit?0x10000-p0:p0;

	debug_dec_update(ec->state, cdf, freq, 0, 0, 0, 0, bit);
	ec->state=freq*(ec->state>>16)+c-cdf;//update
	if(ec->state<0x10000)//renorm
	{
		ec->state<<=16;
		if(ec->srcptr-2>=ec->srcstart)
		{
			ec->srcptr-=2;
			memcpy(&ec->state, ec->srcptr, 2);
		}
	}
	return bit;
}


//Golomb-Rice Coder
typedef struct GolombRiceCoderStruct
{
	unsigned long long cache;
	unsigned nbits;//enc: number of free bits in cache, dec: number of unread bits in cache
	int is_enc;
	const unsigned char *srcptr, *srcend, *srcstart;
	unsigned char *dstptr, *dstend, *dststart;
	DList *list;
} GolombRiceCoder;
INLINE void gr_enc_init(GolombRiceCoder *ec,
#ifdef EC_USE_ARRAY
	unsigned char *start, unsigned char *end
#else
	DList *list
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
	ec->list=list;
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
INLINE size_t gr_enc_flush(GolombRiceCoder *ec)
{
#ifdef EC_USE_ARRAY
	if(ec->dstptr+sizeof(ec->cache)>ec->dstend)//compression failed
		return 0;
	memcpy(ec->dstptr, &ec->cache, sizeof(ec->cache));
	ec->dstptr+=sizeof(ec->cache);
	return 1;
#else
	dlist_push_back(ec->list, &ec->cache, sizeof(ec->cache));//size is qword-aligned
	return 1;
#endif
	//while(ec->nbits<(sizeof(ec->cache)<<3))//loops up to 2 times
	//{
	//	dlist_push_back(ec->list, (unsigned*)&ec->cache+1, 4);
	//	ec->nbits+=32;
	//	ec->cache<<=32;
	//}
}

INLINE int gr_enc(GolombRiceCoder *ec, unsigned sym, unsigned magnitude)
{
	//buffer: {c,c,c,b,b,a,a,a, f,f,f,e,e,e,d,c}, cache: MSB gg[hhh]000 LSB	nbits 6->3, code h is about to be emitted
	//written 64-bit words are byte-reversed because the CPU is big-endian
	magnitude+=!magnitude;
	unsigned nbypass=floor_log2_32(magnitude)+1;
	unsigned nzeros=sym/magnitude, bypass=sym%magnitude;
	if(nzeros>=ec->nbits)//fill the rest of cache with zeros, and flush
	{
		nzeros-=ec->nbits;
		if(!gr_enc_flush(ec))
			return 0;
		//dlist_push_back(ec->list, &ec->cache, sizeof(ec->cache));
		if(nzeros>=(sizeof(ec->cache)<<3))//just flush zeros
		{
			ec->cache=0;
			do
			{
				nzeros-=(sizeof(ec->cache)<<3);
				if(!gr_enc_flush(ec))
					return 0;
				//dlist_push_back(ec->list, &ec->cache, sizeof(ec->cache));
			}
			while(nzeros>(sizeof(ec->cache)<<3));
		}
		ec->cache=0;
		ec->nbits=(sizeof(ec->cache)<<3);
	}
	//now there is room for zeros:  0 <= nzeros < nbits <= 64
	if(nzeros)//emit remaining zeros to cache
		ec->nbits-=nzeros;

	unsigned nunused=(1<<nbypass)-magnitude;//truncated binary code
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
		//dlist_push_back(ec->list, &ec->cache, sizeof(ec->cache));
		ec->cache=0;
		ec->nbits=(sizeof(ec->cache)<<3);
	}
	//now there is room for bypass:  0 <= nbypass < nbits <= 64
	if(nbypass)//emit remaining bypass to cache
	{
		ec->nbits-=nbypass;
		ec->cache|=(unsigned long long)bypass<<ec->nbits;
	}
	return 1;
}
INLINE unsigned gr_dec(GolombRiceCoder *ec, unsigned magnitude)
{
	//cache: MSB 00[hhh]ijj LSB		nbits 6->3, h is about to be read (past codes must be cleared from cache)
	
	magnitude+=!magnitude;
	unsigned nbypass=floor_log2_32(magnitude);
	int sym=0, nleadingzeros=0;
	if(!ec->nbits)//cache is empty
		goto read;
	for(;;)//cache reading loop
	{
		nleadingzeros=ec->nbits-1-floor_log2(ec->cache);//count leading zeros
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
	ec->cache&=(1ULL<<ec->nbits)-1;//remove stop bit

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


#ifdef __cplusplus
}
#endif
#endif
