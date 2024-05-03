#pragma once
#ifndef INC_AC_H
#define INC_AC_H
#include"util.h"
#include<stdio.h>
#include<string.h>
#ifdef _MSC_VER
#include<intrin.h>
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

	if(acval_idx>=(int)acval->count)
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
void debug_dec_update(unsigned state, unsigned cdf, int freq, int kx, int ky, int kq, int kc, unsigned char sym)
{
	if(freq<=0)
		LOG_ERROR2("ANS: invalid frequency %d", freq);
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
#define PROB_BITS 16//can't be changed
#ifdef EC_USE_ARRAY
typedef ArrayHandle EC_DSTStructure;
#else
typedef DList EC_DSTStructure;
#endif
typedef struct ArithmeticCoderStruct
{
	unsigned long long low, range;
	EC_DSTStructure *dst;
	const unsigned char *srcptr, *srcend;
	unsigned long long code;
} ArithmeticCoder;
INLINE void ac_enc_init(ArithmeticCoder *ec, EC_DSTStructure *dst)
{
	memset(ec, 0, sizeof(*ec));
	ec->low=0;
	ec->range=~0LLU>>PROB_BITS;
#ifdef EC_USE_ARRAY
	array_append(dst, 0, 1, 0, 0, 0, 0);
#endif
	ec->dst=dst;
}
INLINE void ac_dec_init(ArithmeticCoder *ec, const unsigned char *start, unsigned const char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->low=0;
	ec->range=~0LLU>>PROB_BITS;
	ec->srcptr=start;
	ec->srcend=end;
	
	if(ec->srcptr+6>ec->srcend)
	{
		LOG_ERROR2("buffer overflow");
		return;
	}
	memcpy((unsigned char*)&ec->code+4, ec->srcptr+0, 2);
	memcpy((unsigned char*)&ec->code+2, ec->srcptr+2, 2);
	memcpy((unsigned char*)&ec->code+0, ec->srcptr+4, 2);
	ec->srcptr+=6;
}
INLINE void ac_enc_renorm(ArithmeticCoder *ec)//fast renorm by F. Rubin 1979
{
#ifdef EC_USE_ARRAY
	array_append(ec->dst, (unsigned char*)&ec->low+4, 1, 2, 1, 0, 0);
#else
	dlist_push_back(ec->dst, (unsigned char*)&ec->low+4, 2);
#endif
	ec->range<<=PROB_BITS;
	ec->low<<=PROB_BITS;
	ec->range|=(1LL<<PROB_BITS)-1;
	ec->low&=~0LLU>>PROB_BITS;

	unsigned long long rmax=ec->low^(~0LLU>>PROB_BITS);
	if(ec->range>rmax)//clamp hi to register size after renorm
		ec->range=rmax;
}
INLINE void ac_dec_renorm(ArithmeticCoder *ec)//fast renorm by F. Rubin 1979
{
//	if(ec->srcptr+2>ec->srcend)
//	{
//#ifdef AC_VALIDATE
//		printf("buffer overflow\n");
//		acval_dump();
//#endif
//		LOG_ERROR2("buffer overflow");
//		return;
//	}
	ec->code<<=PROB_BITS;
	ec->range<<=PROB_BITS;
	ec->low<<=PROB_BITS;
	ec->range|=(1LL<<PROB_BITS)-1;
	if(ec->srcptr+2<=ec->srcend)
		memcpy(&ec->code, ec->srcptr, 2);
	ec->srcptr+=2;
	ec->low&=~0LLU>>PROB_BITS;
	ec->code&=~0LLU>>PROB_BITS;

	unsigned long long rmax=ec->low^(~0LLU>>PROB_BITS);
	if(ec->range>rmax)
		ec->range=rmax;
}
INLINE void ac_enc_flush(ArithmeticCoder *ec)
{
	unsigned long long code=ec->low+ec->range;
	int n=floor_log2(ec->low^code)-1;
	if(n<0)
		n=0;
	code&=~((1LL<<n)-1);		//minimize final code
#if 0
	int flushbits=get_lsb_index(code);//FIXME tail-chaining parallel decoders
#endif
	for(int k=4;k>=0&&code;k-=2)
	{
#ifdef EC_USE_ARRAY
		array_append(ec->dst, (unsigned char*)&code+4, 1, 2, 1, 0, 0);
#else
		dlist_push_back(ec->dst, (unsigned char*)&code+4, 2);
#endif
		code<<=16;
		code&=(~0LLU>>PROB_BITS);
	}
}

INLINE void ac_enc_av2(ArithmeticCoder *ec, int sym, const unsigned *CDF1, const unsigned *CDF2, int alpha)//CDF is 16 bit
{
	unsigned cdf=(CDF1[sym]+CDF2[sym])>>1;
	int freq=((CDF1[sym+1]+CDF2[sym+1])>>1)-cdf;
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
INLINE int ac_dec_packedsign_av2(ArithmeticCoder *ec, const unsigned *CDF1, const unsigned *CDF2, int alpha, int nlevels)//preferred for skewed distributions as with 'pack sign'
{
	unsigned cdf=(unsigned)(((ec->code-ec->low)<<16|0xFFFF)/ec->range);
	int freq;
	int sym;

	(void)nlevels;

	//_mm_prefetch((char*)CDF, _MM_HINT_T0);
	
	for(sym=0;cdf>=((CDF1[sym+2]+CDF2[sym+2])>>1);sym+=2);//must shift right to remove LSB which can corrupt decoded symbol
	sym+=cdf>=((CDF1[sym+1]+CDF2[sym+1])>>1);
	cdf=(CDF1[sym]+CDF2[sym])>>1;
	freq=((CDF1[sym+1]+CDF2[sym+1])>>1)-cdf;
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
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return sym;
}

#define AV4_CDFs(IDX) ((CDF1[IDX]+CDF2[IDX]+CDF3[IDX]+CDF4[IDX])>>2)
INLINE void ac_enc_av4(ArithmeticCoder *ec, int sym, const unsigned *CDF1, const unsigned *CDF2, const unsigned *CDF3, const unsigned *CDF4, int alpha, int beta)//CDF is 16 bit
{
	unsigned cdf=AV4_CDFs(sym);
	int freq=AV4_CDFs(sym+1)-cdf;
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
INLINE int ac_dec_packedsign_av4(ArithmeticCoder *ec, const unsigned *CDF1, const unsigned *CDF2, const unsigned *CDF3, const unsigned *CDF4, int alpha, int beta, int nlevels)//preferred for skewed distributions as with 'pack sign'
{
	unsigned cdf=(unsigned)(((ec->code-ec->low)<<16|0xFFFF)/ec->range);
	int freq;
	int sym;

	(void)nlevels;

	//_mm_prefetch((char*)CDF, _MM_HINT_T0);
	
	for(sym=0;cdf>=AV4_CDFs(sym+2);sym+=2);
	sym+=cdf>=AV4_CDFs(sym+1);
	cdf=AV4_CDFs(sym);
	freq=AV4_CDFs(sym+1)-cdf;
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
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return sym;
}

#define MIX2(X1, X2, ALPHA) (unsigned)(X1+((int)(X2-X1)*(ALPHA)>>15))
INLINE void ac_enc_mix2(ArithmeticCoder *ec, int sym, const unsigned *CDF1, const unsigned *CDF2, int alpha)//CDF is 16 bit
{
	unsigned cdf=MIX2(CDF1[sym], CDF2[sym], alpha);
	int freq=MIX2(CDF1[sym+1], CDF2[sym+1], alpha)-cdf;
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
INLINE int ac_dec_packedsign_mix2(ArithmeticCoder *ec, const unsigned *CDF1, const unsigned *CDF2, int alpha, int nlevels)//preferred for skewed distributions as with 'pack sign'
{
	unsigned cdf=(unsigned)(((ec->code-ec->low)<<16|0xFFFF)/ec->range);
	int freq;
	int sym;

	(void)nlevels;

	//_mm_prefetch((char*)CDF, _MM_HINT_T0);
	
	for(sym=0;cdf>=MIX2(CDF1[sym+2], CDF2[sym+2], alpha);sym+=2);
	sym+=cdf>=MIX2(CDF1[sym+1], CDF2[sym+1], alpha);
	cdf=MIX2(CDF1[sym], CDF2[sym], alpha);
	freq=MIX2(CDF1[sym+1], CDF2[sym+1], alpha)-cdf;
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
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return sym;
}

INLINE unsigned mix4(unsigned x0, unsigned x1, unsigned x2, unsigned x3, int alpha, int beta)
{
	x0=MIX2(x0, x1, alpha);
	x2=MIX2(x2, x3, alpha);
	x0=MIX2(x0, x2, beta);
	return x0;
}
#define MIX4_CDFs(IDX) mix4(CDF1[IDX], CDF2[IDX], CDF3[IDX], CDF4[IDX], alpha, beta)
INLINE void ac_enc_mix4(ArithmeticCoder *ec, int sym, const unsigned *CDF1, const unsigned *CDF2, const unsigned *CDF3, const unsigned *CDF4, int alpha, int beta)//CDF is 16 bit
{
	unsigned cdf=MIX4_CDFs(sym);
	int freq=MIX4_CDFs(sym+1)-cdf;
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
INLINE int ac_dec_packedsign_mix4(ArithmeticCoder *ec, const unsigned *CDF1, const unsigned *CDF2, const unsigned *CDF3, const unsigned *CDF4, int alpha, int beta, int nlevels)//preferred for skewed distributions as with 'pack sign'
{
	unsigned cdf=(unsigned)(((ec->code-ec->low)<<16|0xFFFF)/ec->range);
	int freq;
	int sym;

	(void)nlevels;

	//_mm_prefetch((char*)CDF, _MM_HINT_T0);
	
	for(sym=0;cdf>=MIX4_CDFs(sym+2);sym+=2);
	sym+=cdf>=MIX4_CDFs(sym+1);
	cdf=MIX4_CDFs(sym);
	freq=MIX4_CDFs(sym+1)-cdf;
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
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return sym;
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
INLINE int ac_dec_uniform(ArithmeticCoder *ec, const unsigned *CDF, int nlevels)
{
	unsigned cdf=(unsigned)(((ec->code-ec->low)<<16|0xFFFF)/ec->range);
	int freq;
	int sym;
	
	int range=nlevels;
	sym=0;
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
INLINE int ac_dec(ArithmeticCoder *ec, const unsigned *CDF, int nlevels)//preferred for skewed distributions as with 'pack sign'
{
	unsigned cdf=(unsigned)(((ec->code-ec->low)<<16|0xFFFF)/ec->range);
	int freq;
	int sym=0;
	
	int range=nlevels;
#ifdef UNIFORM_DIST
	sym=0;
#else
	for(range=1;;)//exponential search
	{
		if(range>nlevels)
		{
			sym=range>>1;
			range=nlevels;
			break;
		}
		if(cdf<CDF[range])
		{
			range>>=1;
			sym=range;
			break;
		}
		range<<=1;
	}
#endif
	while(range)//binary search		lg(nlevels) memory accesses per symbol
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
	if(freq<=0||cdf>0x10000||cdf+freq>0x10000)
		LOG_ERROR2("ZPS");
#endif
	ec->low+=ec->range*cdf>>16;
	ec->range*=freq;
	ec->range>>=16;
	--ec->range;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))
		ac_dec_renorm(ec);
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return sym;
}
INLINE int ac_dec_packedsign(ArithmeticCoder *ec, const unsigned *CDF, int nlevels)//preferred for skewed distributions as with 'pack sign'
{
	unsigned cdf=(unsigned)(((ec->code-ec->low)<<16|0xFFFF)/ec->range);
	int freq;
	int sym;

	(void)nlevels;

	//_mm_prefetch((char*)CDF, _MM_HINT_T0);

	for(sym=0;cdf>=CDF[sym+2];sym+=2);
	sym+=cdf>=CDF[sym+1];
	
	//for(sym=0;cdf>=CDF[sym+1];++sym);//slower

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
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return sym;
}
INLINE int ac_dec_POT(ArithmeticCoder *ec, const unsigned *CDF, int nbits)
{
	unsigned cdf;
	int freq;
	int sym=0;
	
	unsigned c=(unsigned)(((ec->code-ec->low)<<16|0xFFFF)/ec->range);
	switch(nbits)
	{
	default:
		LOG_ERROR2("Unsupported bit depth");
		break;
	case 9:sym|=(c>=CDF[sym|256])<<8;
	case 8:sym|=(c>=CDF[sym|128])<<7;
	case 7:sym|=(c>=CDF[sym| 64])<<6;
	case 6:sym|=(c>=CDF[sym| 32])<<5;
	case 5:sym|=(c>=CDF[sym| 16])<<4;
	case 4:sym|=(c>=CDF[sym|  8])<<3;
	case 3:sym|=(c>=CDF[sym|  4])<<2;
	case 2:sym|=(c>=CDF[sym|  2])<<1;
	case 1:sym|= c>=CDF[sym|  1];
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
INLINE int ac_dec_POT_permuted(ArithmeticCoder *ec, const unsigned *pCDF, const unsigned *CDF, int nbits)//~5% faster because more cache-friendly,  pCDF[0] & pCDF[1<<nbits] are not accessed
{
	unsigned cdf;
	int freq;
	unsigned sym=1;
	
	unsigned c=(unsigned)(((ec->code-ec->low)<<16|0xFFFF)/ec->range);
	switch(nbits)
	{
	default:
		LOG_ERROR2("Unsupported bit depth");
		break;
	case 16:sym=sym<<1|(c>=pCDF[sym]);
	case 15:sym=sym<<1|(c>=pCDF[sym]);
	case 14:sym=sym<<1|(c>=pCDF[sym]);
	case 13:sym=sym<<1|(c>=pCDF[sym]);
	case 12:sym=sym<<1|(c>=pCDF[sym]);
	case 11:sym=sym<<1|(c>=pCDF[sym]);
	case 10:sym=sym<<1|(c>=pCDF[sym]);
	case  9:sym=sym<<1|(c>=pCDF[sym]);
	case  8:sym=sym<<1|(c>=pCDF[sym]);
	case  7:sym=sym<<1|(c>=pCDF[sym]);
	case  6:sym=sym<<1|(c>=pCDF[sym]);
	case  5:sym=sym<<1|(c>=pCDF[sym]);
	case  4:sym=sym<<1|(c>=pCDF[sym]);
	case  3:sym=sym<<1|(c>=pCDF[sym]);
	case  2:sym=sym<<1|(c>=pCDF[sym]);
	case  1:sym=sym<<1|(c>=pCDF[sym]);
		break;
	}
	sym-=1<<nbits;
	//sym&=(1<<nbits)-1;

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
INLINE int ac_dec_CDF2sym(ArithmeticCoder *ec, const unsigned *CDF, const unsigned char *CDF2sym)
{
	unsigned cdf;
	int freq;
	unsigned sym;
	
	unsigned c=(unsigned)(((ec->code-ec->low)<<16|0xFFFF)/ec->range);
	sym=CDF2sym[c];

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
INLINE int ac_dec_4bit(ArithmeticCoder *ec, const unsigned *CDF)
{
	unsigned cdf;
	int freq;
	unsigned short sym=0;
	
	unsigned long long c=((ec->code-ec->low)<<16|0xFFFF)/ec->range;
	sym|=8&-(c>=CDF[sym|8]);
	sym|=4&-(c>=CDF[sym|4]);
	sym|=2&-(c>=CDF[sym|2]);
	sym|=1&-(c>=CDF[sym|1]);

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
	
	unsigned long long c=((ec->code-ec->low)<<16|0xFFFF)/ec->range;
	sym|=16&-(c>=CDF[sym|16]);
	sym|= 8&-(c>=CDF[sym| 8]);
	sym|= 4&-(c>=CDF[sym| 4]);
	sym|= 2&-(c>=CDF[sym| 2]);
	sym|= 1&-(c>=CDF[sym| 1]);

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

INLINE void ac_enc_packedCDF(ArithmeticCoder *ec, int sym, const unsigned short *CDF, int nlevels)//CDF is 16 bit
{
	unsigned cdf=CDF[sym];
	int freq=(sym>=nlevels-1?0x10000:CDF[sym+1])-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(freq<=0)
		LOG_ERROR2("ZPS");
	if(cdf>0x10000||cdf+freq>0x10000)
		LOG_ERROR2("");
	if(ec->range==~0LLU)
		LOG_ERROR2("");
	//if(cdf>>15||freq>>15)
	//	LOG_ERROR2("LOL_1");
#endif
	ec->low+=ec->range*cdf>>16;
	ec->range=ec->range*freq>>16;
	--ec->range;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))//only when freq=1 -> range=0, this loop runs twice
		ac_enc_renorm(ec);
	acval_enc(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
INLINE int ac_dec_packedCDF_POT(ArithmeticCoder *ec, const unsigned short *CDF, int nbits)
{
	unsigned cdf;
	int freq;
	int sym=0;
	
	unsigned long long c=((ec->code-ec->low)<<16|0xFFFF)/ec->range;
	switch(nbits)
	{
	default:
		LOG_ERROR2("Unsupported bit depth");
		break;
	case 9:sym|=(c>=CDF[sym|256])<<8;
	case 8:sym|=(c>=CDF[sym|128])<<7;
	case 7:sym|=(c>=CDF[sym| 64])<<6;
	case 6:sym|=(c>=CDF[sym| 32])<<5;
	case 5:sym|=(c>=CDF[sym| 16])<<4;
	case 4:sym|=(c>=CDF[sym|  8])<<3;
	case 3:sym|=(c>=CDF[sym|  4])<<2;
	case 2:sym|=(c>=CDF[sym|  2])<<1;
	case 1:sym|= c>=CDF[sym|  1];
		break;
	}

	cdf=CDF[sym];
	freq=(sym>=(1<<nbits)-1?0x10000:CDF[sym+1])-cdf;
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

INLINE void ac_enc_packedCDF_8x3(ArithmeticCoder *ec, const unsigned char *sym, const unsigned short *CDF0, const unsigned short *CDF1, const unsigned short *CDF2)//CDF is 16 bit
{
	unsigned cdf[]={CDF0[sym[0]], CDF1[sym[1]], CDF2[sym[2]]};
	int freq[]=
	{
		(int)((sym[0]==255?0x10000:CDF0[sym[0]+1])-cdf[0]),
		(int)((sym[1]==255?0x10000:CDF1[sym[1]+1])-cdf[1]),
		(int)((sym[2]==255?0x10000:CDF2[sym[2]+1])-cdf[2]),
	};
#ifdef AC_VALIDATE
	unsigned long long
		lo0[]={ec[0].low, ec[1].low, ec[2].low},
		r0[]={ec[0].range, ec[1].range, ec[2].range};
	if(freq[0]<=0||cdf[0]>0x10000||cdf[0]+freq[0]>0x10000)
		LOG_ERROR2("ZPS");
	if(freq[1]<=0||cdf[1]>0x10000||cdf[1]+freq[1]>0x10000)
		LOG_ERROR2("ZPS");
	if(freq[2]<=0||cdf[2]>0x10000||cdf[2]+freq[2]>0x10000)
		LOG_ERROR2("ZPS");
#endif
	ec[0].low+=ec[0].range*cdf[0]>>16;
	ec[1].low+=ec[1].range*cdf[1]>>16;
	ec[2].low+=ec[2].range*cdf[2]>>16;
	ec[0].range=ec[0].range*freq[0]>>16;
	ec[1].range=ec[1].range*freq[1]>>16;
	ec[2].range=ec[2].range*freq[2]>>16;
	--ec[0].range;//must decrement hi because decoder fails when code == hi2
	--ec[1].range;
	--ec[2].range;
	while(ec[0].range<(1LL<<PROB_BITS))//only when freq=1 -> range=0, this loop runs twice
		ac_enc_renorm(ec+0);
	while(ec[1].range<(1LL<<PROB_BITS))
		ac_enc_renorm(ec+1);
	while(ec[2].range<(1LL<<PROB_BITS))
		ac_enc_renorm(ec+2);
	acval_enc(sym[0], cdf[0], freq[0], lo0[0], lo0[0]+r0[0], ec[0].low, ec[0].low+ec[0].range, 0, 0);
	acval_enc(sym[1], cdf[1], freq[1], lo0[1], lo0[1]+r0[1], ec[1].low, ec[1].low+ec[1].range, 0, 0);
	acval_enc(sym[2], cdf[2], freq[2], lo0[2], lo0[2]+r0[2], ec[2].low, ec[2].low+ec[2].range, 0, 0);
}
INLINE void ac_dec_packedCDF_8x3(ArithmeticCoder *ec, const unsigned short *CDF0, const unsigned short *CDF1, const unsigned short *CDF2, unsigned char *ret_sym)
{
	unsigned cdf[3];
	int freq[3];
	
	unsigned long long c[]=
	{
		((ec[0].code-ec[0].low)<<16|0xFFFF)/ec[0].range,
		((ec[1].code-ec[1].low)<<16|0xFFFF)/ec[1].range,
		((ec[2].code-ec[2].low)<<16|0xFFFF)/ec[2].range,
	};
	ret_sym[0] =(c[0]>=CDF0[          128])<<7;
	ret_sym[1] =(c[1]>=CDF1[          128])<<7;
	ret_sym[2] =(c[2]>=CDF2[          128])<<7;
	ret_sym[0]|=(c[0]>=CDF0[ret_sym[0]|64])<<6;
	ret_sym[1]|=(c[1]>=CDF1[ret_sym[1]|64])<<6;
	ret_sym[2]|=(c[2]>=CDF2[ret_sym[2]|64])<<6;
	ret_sym[0]|=(c[0]>=CDF0[ret_sym[0]|32])<<5;
	ret_sym[1]|=(c[1]>=CDF1[ret_sym[1]|32])<<5;
	ret_sym[2]|=(c[2]>=CDF2[ret_sym[2]|32])<<5;
	ret_sym[0]|=(c[0]>=CDF0[ret_sym[0]|16])<<4;
	ret_sym[1]|=(c[1]>=CDF1[ret_sym[1]|16])<<4;
	ret_sym[2]|=(c[2]>=CDF2[ret_sym[2]|16])<<4;
	ret_sym[0]|=(c[0]>=CDF0[ret_sym[0]| 8])<<3;
	ret_sym[1]|=(c[1]>=CDF1[ret_sym[1]| 8])<<3;
	ret_sym[2]|=(c[2]>=CDF2[ret_sym[2]| 8])<<3;
	ret_sym[0]|=(c[0]>=CDF0[ret_sym[0]| 4])<<2;
	ret_sym[1]|=(c[1]>=CDF1[ret_sym[1]| 4])<<2;
	ret_sym[2]|=(c[2]>=CDF2[ret_sym[2]| 4])<<2;
	ret_sym[0]|=(c[0]>=CDF0[ret_sym[0]| 2])<<1;
	ret_sym[1]|=(c[1]>=CDF1[ret_sym[1]| 2])<<1;
	ret_sym[2]|=(c[2]>=CDF2[ret_sym[2]| 2])<<1;
	ret_sym[0]|= c[0]>=CDF0[ret_sym[0]| 1];
	ret_sym[1]|= c[1]>=CDF1[ret_sym[1]| 1];
	ret_sym[2]|= c[2]>=CDF2[ret_sym[2]| 1];

	cdf[0]=CDF0[ret_sym[0]];
	cdf[1]=CDF1[ret_sym[1]];
	cdf[2]=CDF2[ret_sym[2]];
	freq[0]=(ret_sym[0]==255?0x10000:CDF0[ret_sym[0]+1])-cdf[0];
	freq[1]=(ret_sym[1]==255?0x10000:CDF1[ret_sym[1]+1])-cdf[1];
	freq[2]=(ret_sym[2]==255?0x10000:CDF2[ret_sym[2]+1])-cdf[2];
#ifdef AC_VALIDATE
	unsigned long long
		lo0[]={ec[0].low, ec[1].low, ec[2].low},
		r0[]={ec[0].range, ec[1].range, ec[2].range};
	if(freq[0]<=0||cdf[0]>0x10000||cdf[0]+freq[0]>0x10000)
		LOG_ERROR2("ZPS");
	if(freq[1]<=0||cdf[1]>0x10000||cdf[1]+freq[1]>0x10000)
		LOG_ERROR2("ZPS");
	if(freq[2]<=0||cdf[2]>0x10000||cdf[2]+freq[2]>0x10000)
		LOG_ERROR2("ZPS");
#endif
	ec[0].low+=ec[0].range*cdf[0]>>16;
	ec[1].low+=ec[1].range*cdf[1]>>16;
	ec[2].low+=ec[2].range*cdf[2]>>16;
	ec[0].range=ec[0].range*freq[0]>>16;
	ec[1].range=ec[1].range*freq[1]>>16;
	ec[2].range=ec[2].range*freq[2]>>16;
	--ec[0].range;//must decrement hi because decoder fails when code == hi2
	--ec[1].range;
	--ec[2].range;
	while(ec[0].range<(1LL<<PROB_BITS))
		ac_dec_renorm(ec+0);
	while(ec[1].range<(1LL<<PROB_BITS))
		ac_dec_renorm(ec+1);
	while(ec[2].range<(1LL<<PROB_BITS))
		ac_dec_renorm(ec+2);
	acval_dec(ret_sym[0], cdf[0], freq[0], lo0[0], lo0[0]+r0[0], ec[0].low, ec[0].low+ec[0].range, 0, 0, ec[0].code);
	acval_dec(ret_sym[1], cdf[1], freq[1], lo0[1], lo0[1]+r0[1], ec[1].low, ec[1].low+ec[1].range, 0, 0, ec[1].code);
	acval_dec(ret_sym[2], cdf[2], freq[2], lo0[2], lo0[2]+r0[2], ec[2].low, ec[2].low+ec[2].range, 0, 0, ec[2].code);
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
	acval_enc(bit, bit?p0:0, bit?0x10000-p0:p0, ec->low, ec->low+ec->range, bit?ec->low+r2:ec->low, bit?ec->low+ec->range:ec->low+r2-1, 0, 0);
	//acval_enc(bit, bit?p0:0, bit?0x10000-p0:p0, ec->low, ec->low+ec->range, bit?ec->low+r2:ec->low, bit?ec->low+ec->range:ec->low+r2-1, ec->cache, ec->cidx);
#endif
	ec->low+=r2&-bit;
	ec->range=bit?ec->range-r2:r2-1;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))
		ac_enc_renorm(ec);
}
INLINE int ac_dec_bin(ArithmeticCoder *ec, unsigned short p0)//binary AC decoder doesn't do binary search
{
	unsigned long long r2=ec->range*p0>>16;
	int bit=ec->code>=ec->low+r2;
#ifdef AC_VALIDATE
	if(!p0)//reject degenerate distribution
		LOG_ERROR2("ZPS");
	acval_dec(bit, bit?p0:0, bit?0x10000-p0:p0, ec->low, ec->low+ec->range, bit?ec->low+r2:ec->low, bit?ec->low+ec->range:ec->low+r2-1, 0, 0, ec->code);
	//acval_dec(bit, bit?p0:0, bit?0x10000-p0:p0, ec->low, ec->low+ec->range, bit?ec->low+r2:ec->low, bit?ec->low+ec->range:ec->low+r2-1, ec->cache, ec->cidx, ec->code);
#endif
	ec->low+=r2&-bit;
	ec->range=bit?ec->range-r2:r2-1;
	//if(bit)
	//{
	//	ec->low+=r2;
	//	ec->range-=r2;
	//}
	//else
	//	ec->range=r2-1;//must decrement hi because decoder fails when code == hi2
	while(ec->range<(1LL<<PROB_BITS))
		ac_dec_renorm(ec);
	return bit;
}


//asymmetric numeral systems coder
typedef struct ANSCoderStruct
{
	unsigned state, unused;
	EC_DSTStructure *dst;
	const unsigned char *srcptr, *srcstart;
} ANSCoder;
INLINE void ans_enc_init(ANSCoder *ec, EC_DSTStructure *dst)
{
	memset(ec, 0, sizeof(*ec));
	ec->state=0x10000;
	ec->dst=dst;
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
	ARRAY_APPEND(*ec->dst, &ec->state, 4, 1, 0);
#else
	dlist_push_back(ec->dst, &ec->state, 4);
#endif
}

INLINE void ans_enc(ANSCoder *ec, int sym, const unsigned *CDF, int nlevels)
{
	int cdf, freq;
	if(CDF)
		cdf=CDF[sym], freq=CDF[sym+1]-cdf;
	else//bypass
		cdf=(sym<<16)/nlevels, freq=((sym+1)<<16)/nlevels-cdf;
#ifdef DEBUG_ANS
	if(!freq)
		LOG_ERROR2("ZPS");
#endif
	if((ec->state>>16)>=(unsigned)freq)//renorm
	{
#ifdef EC_USE_ARRAY
		ARRAY_APPEND(*ec->dst, &ec->state, 2, 1, 0);
#else
		dlist_push_back(ec->dst, &ec->state, 2);
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
		while(R)//binary search		lg(nlevels) memory accesses per symbol
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
#ifdef DEBUG_ANS
	if(!freq)
		LOG_ERROR2("ZPS");
#endif
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
	case 1:sym|= CDF[sym|  1]<=c;
		break;
	}

	cdf=CDF[sym], freq=CDF[sym+1]-cdf;
#ifdef DEBUG_ANS
	if(!freq)
		LOG_ERROR2("ZPS");
	debug_dec_update(ec->state, cdf, freq, 0, 0, 0, 0, sym);
#endif
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
		ARRAY_APPEND(*ec->dst, &ec->state, 2, 1, 0);
#else
		dlist_push_back(ec->dst, &ec->state, 2);
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
		ARRAY_APPEND(*ec->dst, &ec->state, 2, 1, 0);
#else
		dlist_push_back(ec->dst, &ec->state, 2);
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
	DList *dst;
} GolombRiceCoder;
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
	dlist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));//size is qword-aligned
	return 1;
#endif
	//while(ec->nbits<(sizeof(ec->cache)<<3))//loops up to 2 times
	//{
	//	dlist_push_back(ec->dst, (unsigned*)&ec->cache+1, 4);
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
		//dlist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		if(nzeros>=(sizeof(ec->cache)<<3))//just flush zeros
		{
			ec->cache=0;
			do
			{
				nzeros-=(sizeof(ec->cache)<<3);
				if(!gr_enc_flush(ec))
					return 0;
				//dlist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
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
		//dlist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
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
