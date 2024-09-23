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
#include"blist.h"
#ifdef __cplusplus
extern "C"
{
#endif
#ifdef __GNUC__
#define FORCEINLINE __attribute__((always_inline)) inline static
#else
#define FORCEINLINE __forceinline static
#endif

	
//	#define AC_VALIDATE

#ifdef AC_VALIDATE
void acval_enc(int sym, int cdf, int freq, unsigned long long lo1, unsigned long long hi1, unsigned long long lo2, unsigned long long hi2, unsigned long long cache, int nbits);
void acval_dec(int sym, int cdf, int freq, unsigned long long lo1, unsigned long long hi1, unsigned long long lo2, unsigned long long hi2, unsigned long long cache, int nbits, unsigned long long code);
void acval_dump();
//extern int acval_disable;
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
extern ArrayHandle acval;
extern ACVAL acval_cbuf[ACVAL_CBUFSIZE];
extern int acval_idx;
#ifdef AC_IMPLEMENTATION
ArrayHandle acval=0;
ACVAL acval_cbuf[ACVAL_CBUFSIZE]={0};
int acval_idx=0;
//int acval_disable=0;
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
	//if(acval_disable)
	//	return;
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
	
	//if(acval_disable)
	//	return;

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
void debug_enc_update(unsigned state, unsigned cdf, int freq, int kx, int ky, int kq, int kc, unsigned char sym);
void debug_dec_update(unsigned state, unsigned cdf, int freq, int kx, int ky, int kq, int kc, unsigned char sym);
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
void debug_enc_update(unsigned state, unsigned cdf, int freq, int kx, int ky, int kq, int kc, unsigned char sym)
{
	if(freq<=0)
		LOG_ERROR2("ANS: invalid frequency %d", freq);
	//if(kc==debug_channel)
	{
		unsigned s0=state;

		if(!states.count)
			slist_init(&states, sizeof(DebugANSInfo), 0);

		//if(state>>16>(unsigned)freq)//enc renorm	X
		//	state>>=16;
		state=state/freq<<16|(cdf+state%freq);//enc update

		DebugANSInfo info={s0, state, cdf, freq, (unsigned)states.count, kx, ky, kq, kc, sym};
		STACK_PUSH(&states, &info);
	}
}
void debug_dec_update(unsigned state, unsigned cdf, int freq, int kx, int ky, int kq, int kc, unsigned char sym)
{
	if(freq<=0)
		LOG_ERROR2("ANS: invalid frequency %d", freq);
	//if(kc==debug_channel)
	{
		if(!states.count)
			LOG_ERROR2("Nothing to decode");
		DebugANSInfo *i0=(DebugANSInfo*)STACK_TOP(&states), info;
		memcpy(&info, i0, sizeof(info));

		//unsigned s0=state;
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


//arithmetic coder (paq8px)
#define AC2_PROB_BITS 16
typedef struct _AC2
{
	unsigned x1, x2, pending_bits, code;
	unsigned char cache, nbits;//enc: number of written bits	dec: number of remaining bits
	BList *dst;
	const unsigned char *srcptr, *srcend;
} AC2;
FORCEINLINE void ac2_bitwrite(AC2 *ec, int bit)
{
	ec->cache=ec->cache<<1|bit;
	++ec->nbits;
	if(ec->nbits==8)
	{
		blist_push_back1(ec->dst, &ec->cache);
		ec->cache=0;
		ec->nbits=0;
	}
}
FORCEINLINE void ac2_bitwrite_withpending(AC2 *ec, int bit)
{
	ac2_bitwrite(ec, bit);
	bit^=1;
	while(ec->pending_bits)
	{
		ac2_bitwrite(ec, bit);
		--ec->pending_bits;
	}
}
FORCEINLINE int ac2_bitread(AC2 *ec)
{
	if(!ec->nbits)
	{
		ec->cache=ec->srcptr<ec->srcend?*ec->srcptr++:0;
		ec->nbits=8;
	}
	--ec->nbits;
	return ec->cache>>ec->nbits&1;
}
FORCEINLINE void ac2_enc_init(AC2 *ec, BList *dst)
{
	memset(ec, 0, sizeof(*ec));
	ec->x1=0;
	ec->x2=~0;
	ec->dst=dst;
}
FORCEINLINE void ac2_dec_init(AC2 *ec, const unsigned char *start, unsigned const char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->x1=0;
	ec->x2=~0;
	ec->srcptr=start;
	ec->srcend=end;
	for(int k=0;k<32;++k)
		ec->code=ec->code<<1|ac2_bitread(ec);
}
FORCEINLINE void ac2_enc_flush(AC2 *ec)
{
	do
	{
		ac2_bitwrite_withpending(ec, ec->x1>>31);
		ec->x1<<=1;
	}while(ec->nbits||ec->x1);//while(ec->nbits);
}
FORCEINLINE void ac2_enc_renorm(AC2 *ec)
{
	while(!((ec->x1^ec->x2)>>31))
	{
		ac2_bitwrite_withpending(ec, ec->x1>>31);
		ec->x1<<=1;
		ec->x2=ec->x2<<1|1;
	}
	while(ec->x1>=0x40000000&&ec->x2<0xC0000000)
	{
		++ec->pending_bits;
		ec->x1=ec->x1<<1&0x7FFFFFFF;
		ec->x2=ec->x2<<1|0x80000001;
	}
#ifdef AC_VALIDATE
	if(ec->x1>=ec->x2)
		LOG_ERROR("Invalid AC range %08X~%08X", ec->x1, ec->x2);
#endif
}
FORCEINLINE void ac2_dec_renorm(AC2 *ec)
{
	while(!((ec->x1^ec->x2)>>31))
	{
		ec->x1<<=1;
		ec->x2=ec->x2<<1|1;
		ec->code=(ec->code<<1)|ac2_bitread(ec);
	}
	while(ec->x1>=0x40000000&&ec->x2<0xC0000000)
	{
		ec->x1=ec->x1<<1&0x7FFFFFFF;
		ec->x2=ec->x2<<1|0x80000001;
		ec->code=(ec->code<<1^0x80000000)+ac2_bitread(ec);
	}
#ifdef AC_VALIDATE
	if(ec->x1>=ec->x2)
		LOG_ERROR("Invalid AC range %08X~%08X", ec->x1, ec->x2);
#endif
}

FORCEINLINE void ac2_enc_update(AC2 *ec, unsigned cdf_curr, unsigned cdf_next)
{
	unsigned range, x1, x2;

	ac2_enc_renorm(ec);
	range=ec->x2-ec->x1;
	x1=ec->x1+(unsigned)((unsigned long long)range*cdf_curr>>AC2_PROB_BITS);
	x2=ec->x1+(unsigned)((unsigned long long)range*cdf_next>>AC2_PROB_BITS)-1;//must decrement hi because decoder fails when code == hi2
	acval_enc(0, cdf_curr, cdf_next-cdf_curr, ec->x1, ec->x2, x1, x2, 0, 0);
	ec->x1=x1;
	ec->x2=x2;
}
FORCEINLINE unsigned ac2_dec_getcdf(AC2 *ec)
{
	unsigned cdf;

	ac2_dec_renorm(ec);
	cdf=(unsigned)(((unsigned long long)(ec->code-ec->x1)<<AC2_PROB_BITS|((1LL<<AC2_PROB_BITS)-1))/(ec->x2-ec->x1));
	return cdf;
}
FORCEINLINE void ac2_dec_update(AC2 *ec, unsigned cdf_curr, unsigned cdf_next)
{
	unsigned range, x1, x2;

	range=ec->x2-ec->x1;
	x1=ec->x1+(unsigned)((unsigned long long)range*cdf_curr>>AC2_PROB_BITS);
	x2=ec->x1+(unsigned)((unsigned long long)range*cdf_next>>AC2_PROB_BITS)-1;//must decrement hi because decoder fails when code == hi2
	acval_dec(0, cdf_curr, cdf_next-cdf_curr, ec->x1, ec->x2, x1, x2, 0, 0, ec->code);
	ec->x1=x1;
	ec->x2=x2;
}
FORCEINLINE void ac2_enc_bypass(AC2 *ec, unsigned sym, int nbits)
{
	unsigned range, x1, x2;
	ac2_enc_renorm(ec);
	range=ec->x2-ec->x1;
	x1=ec->x1+(unsigned)((unsigned long long)range*(sym+0ULL)>>nbits);
	x2=ec->x1+(unsigned)((unsigned long long)range*(sym+1ULL)>>nbits)-1;
	acval_enc(sym, sym, 1, ec->x1, ec->x2, x1, x2, 0, 0);
	ec->x1=x1;
	ec->x2=x2;
}
FORCEINLINE unsigned ac2_dec_bypass(AC2 *ec, int nbits)
{
	unsigned range, x1, x2, sym;
	ac2_dec_renorm(ec);
	sym=(unsigned)(((unsigned long long)(ec->code-ec->x1)<<nbits|((1LL<<nbits)-1))/(ec->x2-ec->x1));
	range=ec->x2-ec->x1;
	x1=ec->x1+(unsigned)((unsigned long long)range*(sym+0ULL)>>nbits);
	x2=ec->x1+(unsigned)((unsigned long long)range*(sym+1ULL)>>nbits)-1;
	acval_dec(sym, sym, 1, ec->x1, ec->x2, x1, x2, 0, 0, ec->code);
	ec->x1=x1;
	ec->x2=x2;
	return sym;
}
FORCEINLINE void ac2_enc_bypass_NPOT(AC2 *ec, unsigned sym, int nlevels)
{
	unsigned range, x1, x2;
	ac2_enc_renorm(ec);
	range=ec->x2-ec->x1;
	x1=ec->x1+(unsigned)((unsigned long long)range*(sym+0ULL)/nlevels);
	x2=ec->x1+(unsigned)((unsigned long long)range*(sym+1ULL)/nlevels)-1;
	acval_enc(sym, sym, 1, ec->x1, ec->x2, x1, x2, 0, 0);
	ec->x1=x1;
	ec->x2=x2;
}
FORCEINLINE unsigned ac2_dec_bypass_NPOT(AC2 *ec, int nlevels)
{
	unsigned range, x1, x2, sym;
	ac2_dec_renorm(ec);
	sym=(unsigned)(((unsigned long long)(ec->code-ec->x1)*nlevels+nlevels-1)/(ec->x2-ec->x1));
	range=ec->x2-ec->x1;
	x1=ec->x1+(unsigned)((unsigned long long)range*(sym+0ULL)/nlevels);
	x2=ec->x1+(unsigned)((unsigned long long)range*(sym+1ULL)/nlevels)-1;
	acval_dec(sym, sym, 1, ec->x1, ec->x2, x1, x2, 0, 0, ec->code);
	ec->x1=x1;
	ec->x2=x2;
	return sym;
}
FORCEINLINE void ac2_enc_bin(AC2 *ec, unsigned p1, int bit)
{
	unsigned mid;

	ac2_enc_renorm(ec);
#ifdef AC_SYMMETRIC
	if(p1<1ULL<<AC2_PROB_BITS>>1)
	{
		p1=(1<<AC2_PROB_BITS)-p1;
		bit^=1;
	}
#endif
#ifdef AC_VALIDATE
	unsigned x1=ec->x1, x2=ec->x2;
	if((unsigned long long)(p1-1LL)>=(1ULL<<AC2_PROB_BITS)-1)
		LOG_ERROR2("Invalid probability %08X", p1);
#endif
	mid=ec->x1+(unsigned)((unsigned long long)(ec->x2-ec->x1)*p1>>AC2_PROB_BITS);
	if(bit)
		ec->x2=mid;
	else
		ec->x1=mid+1;
	acval_enc(bit, 0, p1, x1, x2, ec->x1, ec->x2, 0, 0);
}
FORCEINLINE int ac2_dec_bin(AC2 *ec, unsigned p1)
{
	unsigned mid;
	int bit;

	ac2_dec_renorm(ec);
#ifdef AC_SYMMETRIC
	int flip=p1<1U<<AC2_PROB_BITS>>1;
	if(flip)
		p1=(1<<AC2_PROB_BITS)-p1;
#endif
#ifdef AC_VALIDATE
	unsigned x1=ec->x1, x2=ec->x2;
	if((unsigned long long)(p1-1LL)>=(1ULL<<AC2_PROB_BITS)-1)
		LOG_ERROR2("Invalid probability %08X", p1);
#endif
	mid=ec->x1+(unsigned)((unsigned long long)(ec->x2-ec->x1)*p1>>AC2_PROB_BITS);
	bit=ec->code<=mid;
	if(bit)
		ec->x2=mid;
	else
		ec->x1=mid+1;
	acval_dec(bit, 0, p1, x1, x2, ec->x1, ec->x2, 0, 0, ec->code);
#ifdef AC_SYMMETRIC
	bit^=flip;
#endif
	return bit;
}


//arithmetic coder (carry-less)
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
	BList *dst;
	const unsigned char *srcptr, *srcend;
	unsigned long long code;
} AC3;
FORCEINLINE void ac3_enc_init(AC3 *ec, BList *dst)
{
	memset(ec, 0, sizeof(*ec));
	ec->low=0;
	ec->range=~0LLU;
	ec->dst=dst;
}
FORCEINLINE void ac3_dec_init(AC3 *ec, const unsigned char *start, unsigned const char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->low=0;
	ec->range=~0LLU;
	ec->srcptr=start;
	ec->srcend=end;
	
	//if(ec->srcptr+8>ec->srcend)
	//{
	//	LOG_ERROR2("buffer overflow");
	//	return;
	//}
	int nbytes=8;
	if(ec->srcptr+8>ec->srcend)
		nbytes=(int)(ec->srcend-ec->srcptr);
	for(int k=0;k<nbytes;k+=AC3_RENORM/8)
		memcpy((unsigned char*)&ec->code+8-AC3_RENORM/8-k, ec->srcptr+k, AC3_RENORM/8);
	ec->srcptr+=nbytes;
}
FORCEINLINE void ac3_enc_renorm(AC3 *ec)//fast renorm by F. Rubin 1979
{
	unsigned long long rmax;

	blist_push_back(ec->dst, (unsigned char*)&ec->low+8-AC3_RENORM/8, AC3_RENORM/8);
	ec->range=ec->range<<AC3_RENORM|((1LL<<AC3_RENORM)-1);
	ec->low<<=AC3_RENORM;

	rmax=~ec->low;
	if(ec->range>rmax)//clamp hi to register size after renorm
		ec->range=rmax;
}
FORCEINLINE void ac3_dec_renorm(AC3 *ec)//fast renorm by F. Rubin 1979
{
	unsigned long long rmax;

	ec->range=ec->range<<AC3_RENORM|((1LL<<AC3_RENORM)-1);
	ec->code<<=AC3_RENORM;
	ec->low<<=AC3_RENORM;
	if(ec->srcptr+AC3_RENORM/8<=ec->srcend)
	{
		memcpy(&ec->code, ec->srcptr, AC3_RENORM/8);
		ec->srcptr+=AC3_RENORM/8;
	}

	rmax=~ec->low;
	if(ec->range>rmax)
		ec->range=rmax;
}
FORCEINLINE void ac3_enc_flush(AC3 *ec)
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
		blist_push_back(ec->dst, (unsigned char*)&code+8-(AC3_RENORM/8), AC3_RENORM/8);
		code<<=AC3_RENORM;
	}
}

FORCEINLINE void ac3_enc_update(AC3 *ec, unsigned cdf, unsigned freq)
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
FORCEINLINE unsigned ac3_dec_getcdf(AC3 *ec)
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
FORCEINLINE void ac3_dec_update(AC3 *ec, unsigned cdf, unsigned freq)
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

FORCEINLINE void ac3_enc_update_NPOT(AC3 *ec, unsigned cdf, unsigned freq, unsigned den)
{
	unsigned long long q;
#ifdef AC3_PREC
	unsigned long long r;
#endif
#ifdef AC_VALIDATE
	unsigned long long lo0, r0;
	if(!freq)
		LOG_ERROR2("ZPS");
	if(cdf+freq<cdf||cdf+freq>den)
		LOG_ERROR2("Invalid CDF");
#endif
	AC3_RENORM_STATEMENT(ec->range<(unsigned)den)//only when freq=1 -> range=0, this loop runs twice
		ac3_enc_renorm(ec);
#ifdef AC_VALIDATE
	lo0=ec->low, r0=ec->range;
#endif
	q=ec->range/den;
#ifdef AC3_PREC
	r=ec->range%den;
	ec->low+=q*cdf+r*cdf/den;
	ec->range=q*freq+r*freq/den-1;
#else
	ec->low+=q*cdf;
	ec->range=q*freq-1;//must decrement hi because decoder fails when code == hi2
#endif
	acval_enc(0, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
FORCEINLINE unsigned ac3_dec_getcdf_NPOT(AC3 *ec, unsigned den)
{
	AC3_RENORM_STATEMENT(ec->range<(unsigned)den)
		ac3_dec_renorm(ec);
#ifdef AC3_PREC
	{
#ifdef AC_VALIDATE
		ACVAL *val=(ACVAL*)array_at(&acval, acval_idx);
		if(ec->code<val->lo2||ec->code>val->hi2)
			printf("");
#endif
		unsigned long long lo, hi, lo0;
		lo0=lo=_umul128(ec->code-ec->low, den, &hi);
		lo+=den-1LL;
		hi+=lo<lo0;
		return (unsigned)_udiv128(hi, lo, ec->range, &hi);
	}
#else
	return (unsigned)((ec->code-ec->low)/(ec->range/den));
#endif
}
FORCEINLINE void ac3_dec_update_NPOT(AC3 *ec, unsigned cdf, unsigned freq, unsigned den)
{
	unsigned long long q=ec->range/den;
#ifdef AC3_PREC
	unsigned long long r=ec->range%den;
#endif
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
	ec->low+=q*cdf;
	ec->range=q*freq-1;//must decrement hi because decoder fails when code == hi2
#endif
	acval_dec(0, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
}

FORCEINLINE void ac3_enc_bypass(AC3 *ec, int bypass, int nbits)
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
FORCEINLINE int ac3_dec_bypass(AC3 *ec, int nbits)
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
FORCEINLINE void ac3_enc_bypass_NPOT(AC3 *ec, int bypass, int nlevels)
{
	unsigned long long q;
#ifdef AC3_PREC
	unsigned long long r;
#endif
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
#endif
	
	AC3_RENORM_STATEMENT(ec->range<(unsigned)nlevels)
		ac3_enc_renorm(ec);
	q=ec->range/nlevels;
#ifdef AC3_PREC
	r=ec->range%nlevels;
	ec->low+=q*bypass+r*bypass/nlevels;
	ec->range=q-1;
#else
	ec->low+=q*bypass;
	ec->range=q-1;
#endif
	acval_enc(bypass, bypass, 1, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
FORCEINLINE int ac3_dec_bypass_NPOT(AC3 *ec, int nlevels)
{
	unsigned long long q;
#ifdef AC3_PREC
	unsigned long long r;
#endif
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
	q=ec->range/nlevels;
	ec->low+=q*bypass;
	ec->range=q-1;
#endif
	acval_dec(bypass, bypass, 1, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return bypass;
}

FORCEINLINE void ac3_enc(AC3 *ec, int sym, const unsigned *CDF)
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
FORCEINLINE int ac3_dec(AC3 *ec, const unsigned *CDF, int nlevels)
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
FORCEINLINE void ac3_enc_bin(AC3 *ec, int bit, unsigned p0, int probbits)
{
	unsigned long long mid;

#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!p0||p0>=1ULL<<probbits)
		LOG_ERROR2("ZPS");
#endif
	mid=ec->range>>probbits;
	AC3_RENORM_STATEMENT(!mid)
	{
		ac3_enc_renorm(ec);
		mid=ec->range>>probbits;
	}
#ifdef AC_SYMMETRIC
	if(p0>1U<<probbits>>1)
	{
		p0=(1<<probbits)-p0;
		bit=!bit;
	}
#endif
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
FORCEINLINE int ac3_dec_bin(AC3 *ec, unsigned p0, int probbits)
{
	unsigned long long mid, t2;
	int bit;
	
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!p0||p0>=1ULL<<probbits)
		LOG_ERROR2("ZPS");
#endif
	mid=ec->range>>probbits;
	AC3_RENORM_STATEMENT(!mid)
	{
		ac3_dec_renorm(ec);
		mid=ec->range>>probbits;
	}
	
#ifdef AC_SYMMETRIC
	int flip=p0>1U<<probbits>>1;
	if(flip)
		p0=(1<<probbits)-p0;
#endif
	mid*=p0;
#ifdef AC3_PREC
	mid+=(ec->range&((1LL<<probbits)-1))*p0>>probbits;
#endif
	t2=ec->low+mid;
	ec->range-=mid;
	bit=ec->code>=t2;
	if(bit)
		ec->low=t2;
	else
		ec->range=mid-1;
	acval_dec(bit, bit?p0:0, bit?(1<<probbits)-p0:p0, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
#ifdef AC_SYMMETRIC
	bit^=flip;
#endif
	return bit;
}
#define AC3_ENC_BIN(EC, BIT, P0, PROB_BITS)\
	do\
	{\
		unsigned long long _mid=(EC)->range>>PROB_BITS;\
		if(!_mid)\
		{\
			blist_push_back((EC)->dst, (unsigned char*)&(EC)->low+8-AC3_RENORM/8, AC3_RENORM/8);\
			(EC)->range=(EC)->range<<AC3_RENORM|((1LL<<AC3_RENORM)-1);\
			(EC)->low<<=AC3_RENORM;\
			_mid=~(EC)->low;\
			if((EC)->range>_mid)\
				(EC)->range=_mid;\
			_mid=(EC)->range>>PROB_BITS;\
		}\
		_mid*=P0;\
		_mid+=((EC)->range&((1LL<<PROB_BITS)-1))*P0>>PROB_BITS;\
		if(BIT)\
		{\
			(EC)->low+=_mid;\
			(EC)->range-=_mid;\
		}\
		else\
			(EC)->range=_mid-1;\
	}while(0)
#define AC3_DEC_BIN(EC, DSTBIT, P0, PROB_BITS)\
	do\
	{\
		unsigned long long _mid=(EC)->range>>PROB_BITS, _t2;\
		if(!_mid)\
		{\
			(EC)->range=(EC)->range<<AC3_RENORM|((1LL<<AC3_RENORM)-1);\
			(EC)->code<<=AC3_RENORM;\
			(EC)->low<<=AC3_RENORM;\
			if((EC)->srcptr+AC3_RENORM/8<=(EC)->srcend)\
			{\
				memcpy(&(EC)->code, (EC)->srcptr, AC3_RENORM/8);\
				(EC)->srcptr+=AC3_RENORM/8;\
			}\
			_mid=~(EC)->low;\
			if((EC)->range>_mid)\
				(EC)->range=_mid;\
			_mid=(EC)->range>>PROB_BITS;\
		}\
		_mid*=P0;\
		_mid+=((EC)->range&((1LL<<PROB_BITS)-1))*P0>>PROB_BITS;\
		_t2=(EC)->low+_mid;\
		(EC)->range-=_mid;\
		DSTBIT=(EC)->code>=_t2;\
		(EC)->low=DSTBIT?_t2:(EC)->low;\
		(EC)->range=DSTBIT?(EC)->range:_mid-1;\
	}while(0)


//arithmetic coder (zpaq)	https://mattmahoney.net/dc/dce.html#Section_32
typedef struct _AC4
{
	unsigned lo, hi;
	BList *dst;
	const unsigned char *srcptr, *srcend;
	unsigned code;
} AC4;
FORCEINLINE void ac4_enc_init(AC4 *ec, BList *dst)
{
	memset(ec, 0, sizeof(*ec));
	ec->lo=0;
	ec->hi=0xFFFFFFFF;
	ec->dst=dst;
}
FORCEINLINE void ac4_dec_init(AC4 *ec, const unsigned char *srcstart, const unsigned char *srcend)
{
	memset(ec, 0, sizeof(*ec));
	ec->lo=0;
	ec->hi=0xFFFFFFFF;
	ec->srcptr=srcstart;
	ec->srcend=srcend;
	
	ec->code=0;
	for(int k=0;k<4;++k)
	{
		ec->code<<=8;
		if(ec->srcptr<ec->srcend)
			ec->code|=*ec->srcptr++;
	}
}
FORCEINLINE void ac4_enc_flush(AC4 *ec)
{
	unsigned mid=(unsigned)(((unsigned long long)ec->lo+ec->hi)>>1);
	blist_push_back1(ec->dst, (unsigned char*)&mid+3);//big-endian
	blist_push_back1(ec->dst, (unsigned char*)&mid+2);
	blist_push_back1(ec->dst, (unsigned char*)&mid+1);
	blist_push_back1(ec->dst, (unsigned char*)&mid+0);
}
FORCEINLINE void ac4_enc_bin(AC4 *ec, unsigned short p1, int bit)
{
	unsigned mid;
	
	while((ec->hi^ec->lo)<0x1000000)
	{
		blist_push_back1(ec->dst, (unsigned char*)&ec->lo+3);
		ec->hi=ec->hi<<8|255;
		ec->lo<<=8;
	}
#ifdef AC_VALIDATE
	unsigned lo=ec->lo, hi=ec->hi;
	if(!p1)
		LOG_ERROR("ZPS");
#endif
	mid=ec->lo+(unsigned)((unsigned long long)(ec->hi-ec->lo)*p1>>16);
	if(bit)
		ec->hi=mid;
	else
		ec->lo=mid+1;
	acval_enc(bit, bit, p1, lo, hi, ec->lo, ec->hi, 0, 0);
}
FORCEINLINE int ac4_dec_bin(AC4 *ec, unsigned short p1)
{
	unsigned mid;
	int bit;
	
	while((ec->hi^ec->lo)<0x1000000)
	{
		ec->hi=ec->hi<<8|255;
		ec->lo<<=8;
		ec->code<<=8;
		if(ec->srcptr<ec->srcend)
			ec->code|=*ec->srcptr++;
	}
#ifdef AC_VALIDATE
	unsigned lo=ec->lo, hi=ec->hi;
	if(!p1)
		LOG_ERROR("ZPS");
#endif
	mid=ec->lo+(unsigned)((unsigned long long)(ec->hi-ec->lo)*p1>>16);
	bit=ec->code<=mid;
	if(bit)
		ec->hi=mid;
	else
		ec->lo=mid+1;
	acval_dec(bit, bit, p1, lo, hi, ec->lo, ec->hi, 0, 0, ec->code);
	return bit;
}
FORCEINLINE void ac4_enc_update_NPOT(AC4 *ec, unsigned cdf_curr, unsigned cdf_next, unsigned den)
{
	unsigned lo, hi;

	while((ec->hi^ec->lo)<0x1000000)
	{
		blist_push_back1(ec->dst, (unsigned char*)&ec->lo+3);
		ec->hi=ec->hi<<8|255;
		ec->lo<<=8;
	}
	lo=ec->lo+(unsigned)((unsigned long long)(ec->hi-ec->lo)*cdf_curr/den)+1;
	hi=ec->lo+(unsigned)((unsigned long long)(ec->hi-ec->lo)*cdf_next/den);
	acval_enc(0, cdf_curr, cdf_next-cdf_curr, ec->lo, ec->hi, lo, hi, 0, 0);
	ec->lo=lo;
	ec->hi=hi;
}
FORCEINLINE unsigned ac4_dec_getcdf_NPOT(AC4 *ec, unsigned den)
{
	unsigned cdf;
	
	while((ec->hi^ec->lo)<0x1000000)
	{
		ec->hi=ec->hi<<8|255;
		ec->lo<<=8;
		ec->code<<=8;
		if(ec->srcptr<ec->srcend)
			ec->code|=*ec->srcptr++;
	}
	cdf=(unsigned)(((unsigned long long)(ec->code-ec->lo)*den+den-1)/(ec->hi-ec->lo));
	return cdf;
}
FORCEINLINE void ac4_dec_update_NPOT(AC4 *ec, unsigned cdf_curr, unsigned cdf_next, unsigned den)
{
	unsigned lo, hi;

	lo=ec->lo+(unsigned)((unsigned long long)(ec->hi-ec->lo)*cdf_curr/den)+1;
	hi=ec->lo+(unsigned)((unsigned long long)(ec->hi-ec->lo)*cdf_next/den);
	acval_dec(0, cdf_curr, cdf_next-cdf_curr, ec->lo, ec->hi, lo, hi, 0, 0, ec->code);
	ec->lo=lo;
	ec->hi=hi;
}


#define AC5_PROB_BITS 16	//10
typedef struct _AC5
{
	unsigned lo, hi, code, enc;
	FILE *f;
} AC5;
FORCEINLINE void ac5_enc_init(AC5 *ec, FILE *fdst)
{
	ec->lo=0;
	ec->hi=0xFFFFFFFF;
	ec->code=0;
	ec->enc=1;
	ec->f=fdst;
}
FORCEINLINE void ac5_dec_init(AC5 *ec, FILE *fsrc)
{
	ec->lo=0;
	ec->hi=0xFFFFFFFF;
	ec->enc=0;
	ec->f=fsrc;
	ec->code=fgetc(fsrc)&255;
	ec->code=ec->code<<8|(fgetc(fsrc)&255);
	ec->code=ec->code<<8|(fgetc(fsrc)&255);
	ec->code=ec->code<<8|(fgetc(fsrc)&255);
}
FORCEINLINE void ac5_enc_flush(AC5 *ec)
{
	unsigned code=ec->lo;
	fputc(code>>24, ec->f);	code<<=8;
	fputc(code>>24, ec->f);	code<<=8;
	fputc(code>>24, ec->f);	code<<=8;
	fputc(code>>24, ec->f);
}
FORCEINLINE void ac5_enc_bin(AC5 *ec, int p1, int bit)
{
	while((ec->lo^ec->hi)<0x1000000)
	{
		fputc(ec->lo>>24, ec->f);
		ec->lo<<=8;
		ec->hi=ec->hi<<8|255;
	}
	{
		unsigned mid=ec->lo+(unsigned)((unsigned long long)(ec->hi-ec->lo)*p1>>AC5_PROB_BITS);
	//	unsigned mid=ec->lo+((ec->hi-ec->lo)>>AC5_PROB_BITS)*p1;
		ec->lo=bit?ec->lo:mid+1;
		ec->hi=bit?mid:ec->hi;
	}
}
FORCEINLINE int ac5_dec_bin(AC5 *ec, int p1)
{
	while((ec->lo^ec->hi)<0x1000000)
	{
		ec->code=ec->code<<8|(fgetc(ec->f)&255);
		ec->lo<<=8;
		ec->hi=ec->hi<<8|255;
	}
	{
		unsigned mid=ec->lo+(unsigned)((unsigned long long)(ec->hi-ec->lo)*p1>>AC5_PROB_BITS);
	//	unsigned mid=ec->lo+((ec->hi-ec->lo)>>AC5_PROB_BITS)*p1;
		int bit=ec->code<=mid;
		ec->lo=bit?ec->lo:mid+1;
		ec->hi=bit?mid:ec->hi;
		return bit;
	}
}


//Golomb-Rice Coder

	typedef unsigned long long GREmit_t;//0.001% larger, 18% faster enc
//	typedef unsigned GREmit_t;

typedef struct GolombRiceCoderStruct
{
	GREmit_t cache;
	int nbits;//enc: number of free bits in cache, dec: number of unread bits in cache
	int is_enc;//for padding
	const unsigned char *srcptr, *srcend, *srcstart;
	unsigned char *dstptr, *dstend, *dststart;
	BList *dst;
} GolombRiceCoder;
FORCEINLINE size_t gr_enc_flush(GolombRiceCoder *ec)
{
#ifdef EC_USE_ARRAY
	if(ec->dstptr+sizeof(ec->cache)>ec->dstend)//compression failed
		return 0;
	memcpy(ec->dstptr, &ec->cache, sizeof(ec->cache));
	ec->dstptr+=sizeof(ec->cache);
	return 1;
#else
	blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));//size is qword-aligned
	return 1;
#endif
}
FORCEINLINE int gr_dec_impl_read(GolombRiceCoder *ec)
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
FORCEINLINE void gr_enc_init(GolombRiceCoder *ec,
#ifdef EC_USE_ARRAY
	unsigned char *start, unsigned char *end
#else
	BList *dst
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
FORCEINLINE void gr_dec_init(GolombRiceCoder *ec, const unsigned char *start, const unsigned char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->cache=0;
	ec->nbits=0;
	ec->is_enc=0;
	ec->srcptr=start;
	ec->srcend=end;
	ec->srcstart=start;
}

FORCEINLINE int gr_enc_NPOT(GolombRiceCoder *ec, unsigned sym, unsigned magnitude)
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
		//blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		if(nzeros>=(int)(sizeof(ec->cache)<<3))//just flush zeros
		{
			ec->cache=0;
			do
			{
				nzeros-=(sizeof(ec->cache)<<3);
				if(!gr_enc_flush(ec))
					return 0;
				//blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
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
		ec->cache|=(GREmit_t)bypass>>nbypass;
		bypass&=(1<<nbypass)-1;
		if(!gr_enc_flush(ec))
			return 0;
		//blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		ec->cache=0;
		ec->nbits=sizeof(ec->cache)<<3;
	}
	//now there is room for bypass:  0 <= nbypass < nbits <= 64
	ec->nbits-=nbypass;//emit remaining bypass to cache
	ec->cache|=(GREmit_t)bypass<<ec->nbits;
	return 1;
}
FORCEINLINE unsigned gr_dec_NPOT(GolombRiceCoder *ec, unsigned magnitude)
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

FORCEINLINE int gr_enc(GolombRiceCoder *ec, int sym, int nbypass)
{
	//buffer: {c,c,c,b,b,a,a,a, f,f,f,e,e,e,d,c}, cache: MSB gg[hhh]000 LSB	nbits 6->3, code h is about to be emitted
	//written 64-bit words are byte-reversed because the CPU is little-endian

	int nzeros=sym>>nbypass, bypass=sym&((1<<nbypass)-1);
	if(nzeros>=ec->nbits)//fill the rest of cache with zeros, and flush
	{
		nzeros-=ec->nbits;
		if(!gr_enc_flush(ec))
			return 0;
		//blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		if(nzeros>=(int)(sizeof(ec->cache)<<3))//just flush zeros
		{
			ec->cache=0;
			do
			{
				nzeros-=(sizeof(ec->cache)<<3);
				if(!gr_enc_flush(ec))
					return 0;
				//blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
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
		ec->cache|=(GREmit_t)bypass>>nbypass;
		bypass&=(1<<nbypass)-1;
		if(!gr_enc_flush(ec))
			return 0;
		//blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		ec->cache=0;
		ec->nbits=sizeof(ec->cache)<<3;
	}
	//now there is room for bypass:  0 <= nbypass < nbits <= 64
	ec->nbits-=nbypass;//emit remaining bypass to cache
	ec->cache|=(GREmit_t)bypass<<ec->nbits;
	return 1;
}
FORCEINLINE unsigned gr_dec(GolombRiceCoder *ec, int nbypass)
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
