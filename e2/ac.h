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
void    acval_enc(int sym, int cdf, int freq, unsigned lo1, unsigned hi1, unsigned lo2, unsigned hi2, unsigned cache, int nbits);
void    acval_dec(int sym, int cdf, int freq, unsigned lo1, unsigned hi1, unsigned lo2, unsigned hi2, unsigned cache, int nbits, unsigned code);
void    acval_dump();
#ifdef AC_IMPLEMENTATION
typedef struct ACVALStruct
{
	int sym, cdf, freq;
	unsigned
		lo1, hi1,
		lo2, hi2,
		cache;
	int nbits;
	unsigned code;
} ACVAL;
#define ACVAL_CBUFSIZE 128
ArrayHandle acval=0;
ACVAL acval_cbuf[ACVAL_CBUFSIZE]={0};
int acval_idx=0;
void acval_enc(int sym, int cdf, int freq, unsigned lo1, unsigned hi1, unsigned lo2, unsigned hi2, unsigned cache, int nbits)
{
	ACVAL val=
	{
		sym, cdf, freq,
		lo1, hi1,
		lo2, hi2,
		cache, nbits,
		0,
	};
	if(!acval)
		ARRAY_ALLOC(ACVAL, acval, 0, 0, 0, 0);
	ARRAY_APPEND(acval, &val, 1, 1, 0);
}
void acval_print(int idx, ACVAL const *p, int dec)
{
	printf("%4d sym 0x%02X cdf 0x%04X~0x%04X range 0x%08X~0x%08X->0x%08X~0x%08X cache 0x%08X:%2d %s", idx, p->sym, p->cdf, p->freq, p->lo1, p->hi1, p->lo2, p->hi2, p->cache, p->nbits, dec?"<-|":"|->");
	if(dec)
	{
		printf(" code 0x%08X ", p->code);
		print_bin32(p->code);
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
		if(k==end)
			printf("\n");
		acval_print(k2, acval_cbuf+k, 1);
		if(k==end)
			break;
	}
}
void acval_dec(int sym, int cdf, int freq, unsigned lo1, unsigned hi1, unsigned lo2, unsigned hi2, unsigned cache, int nbits, unsigned code)
{
	ACVAL *val, val2=
	{
		sym, cdf, freq,
		lo1, hi1,
		lo2, hi2,
		cache, nbits,
		code,
	};

	acval_cbuf[acval_idx%ACVAL_CBUFSIZE]=val2;
	
	//if(acval_idx>=acval_start&&acval_idx<acval_end)
	//{
	//	printf("%4d sym 0x%02X cdf 0x%04X freq 0x%04X  0x%08X~0x%08X->0x%08X~0x%08X cache 0x%08X:%2d right  code 0x%08X ", acval_idx, sym, cdf, freq, lo1, hi1, lo2, hi2, cache, nbits, code);
	//	print_bin32(code);
	//	printf("\n");
	//	//printf("%4d sym 0x%02X cdf 0x%04X freq 0x%04X  0x%08X~0x%08X->0x%08X~0x%08X cache 0x%08X:%2d code 0x%08X\n", acval_idx, sym, cdf, freq, lo1, hi1, lo2, hi2, cache, nbits, code);
	//}

	if(acval_idx>=acval->count)
		LOG_ERROR2("AC validation index error");

	val=(ACVAL*)array_at(&acval, acval_idx);
	if(sym!=val->sym||cdf!=val->cdf||freq!=val->freq)
	{
		acval_dump();
		//printf("%4d sym 0x%02X cdf 0x%04X freq 0x%04X  0x%08X~0x%08X->0x%08X~0x%08X cache 0x%08X:%2d right  code 0x%08X ", acval_idx, sym, cdf, freq, lo1, hi1, lo2, hi2, cache, nbits, code);
		//print_bin32(code);
		//printf("\n");
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


//adaptive binary arithmetic coder
typedef struct ABACEncContextStruct
{
	unsigned lo, hi, cache;
	int nbits;
	DList *list;
} ABACEncContext;
inline void abac_enc_init(ABACEncContext *ctx, DList *list)
{
	ctx->lo=0;
	ctx->hi=0xFFFFFFFF;
	ctx->cache=0;
	ctx->nbits=0;
	ctx->list=list;
}
inline void abac_enc(ABACEncContext *ctx, unsigned short p0, int bit)
{
	unsigned mid;

	for(;;)//renorm
	{
		mid=(ctx->hi-ctx->lo)>>16;
		if(ctx->lo+3>ctx->lo&&ctx->lo+3<=ctx->hi&&mid>0)
			break;
		if(ctx->nbits>=32)
		{
			dlist_push_back(ctx->list, &ctx->cache, 4);
			ctx->cache=0;
			ctx->nbits=0;
		}
		ctx->cache|=(ctx->lo&0x80000000)>>ctx->nbits;//cache is written MSB -> LSB
		++ctx->nbits;

		ctx->lo<<=1;//shift out MSB
		ctx->hi<<=1;
		ctx->hi|=1;
	}

	if(!p0)//reject degenerate distribution
		LOG_ERROR2("ZPS");

	mid=ctx->lo+(int)((unsigned long long)(ctx->hi-ctx->lo)*p0>>16);

	if(bit)
		ctx->lo=mid;
	else
		ctx->hi=mid-1;//OBLIGATORY range leak guard
}
inline void abac_enc_flush(ABACEncContext *ctx)
{
	int k2=0;
	do//flush
	{
		while(ctx->nbits<32)
		{
			ctx->cache|=(ctx->lo&0x80000000)>>ctx->nbits;//cache is written MSB -> LSB
			++ctx->nbits;
			++k2;

			ctx->lo<<=1;//shift out MSB
			ctx->hi<<=1;
			ctx->hi|=1;
		}
		dlist_push_back(ctx->list, &ctx->cache, 4);
		ctx->cache=0;
		ctx->nbits=0;
	}while(k2<32);
}
typedef struct ABACDecContextStruct
{
	unsigned lo, hi, cache;
	int nbits;
	unsigned code;
	const unsigned char *srcptr, *srcend;
} ABACDecContext;
inline void abac_dec_init(ABACDecContext *ctx, const unsigned char *start, unsigned const char *end)
{
	ctx->lo=0;
	ctx->hi=0xFFFFFFFF;
	ctx->nbits=32;
	ctx->srcptr=start;
	ctx->srcend=end;

	if(ctx->srcend-ctx->srcptr<4)
		LOG_ERROR2("buffer overflow");
	memcpy(&ctx->code, ctx->srcptr, 4);
	ctx->srcptr+=4;

	if(ctx->srcend-ctx->srcptr<4)
		LOG_ERROR2("buffer overflow");
	memcpy(&ctx->cache, ctx->srcptr, 4);
	ctx->srcptr+=4;
}
inline int abac_dec(ABACDecContext *ctx, unsigned short p0)
{
	unsigned mid;
	for(;;)//renorm
	{
		mid=(ctx->hi-ctx->lo)>>16;
		if(ctx->lo<ctx->hi&&mid>0)
			break;
		if(!ctx->nbits)
		{
			if(ctx->srcend-ctx->srcptr<4)
			{
#ifdef AC_VALIDATE
				printf("buffer overflow\n");
				acval_dump();
#endif
				LOG_ERROR2("buffer overflow");
			}
			memcpy(&ctx->cache, ctx->srcptr, 4);
			ctx->srcptr+=4;

			ctx->nbits=32;
		}
		--ctx->nbits;
		ctx->code<<=1;//shift out MSB		cache is read MSB -> LSB
		ctx->code|=(unsigned)(ctx->cache>>ctx->nbits&1);

		ctx->lo<<=1;
		ctx->hi<<=1;
		ctx->hi|=1;
	}

	if(!p0)//reject degenerate distribution
		LOG_ERROR2("ZPS");

	mid=ctx->lo+(int)((unsigned long long)(ctx->hi-ctx->lo)*p0>>16);
	int bit=ctx->code>mid;
	
	if(bit)
		ctx->lo=mid;
	else
		ctx->hi=mid-1;//OBLIGATORY range leak guard

	return bit;
}


//arithmetic coder
typedef struct ACEncContextStruct
{
	unsigned lo, hi, cache;
	int nbits;
	unsigned const *CDF2;
	DList *list;
} ACEncContext;
inline void ac_enc_init(ACEncContext *ctx, unsigned const *CDF2, DList *list)
{
	ctx->lo=0;
	ctx->hi=0xFFFFFFFF;
	ctx->cache=0;
	ctx->nbits=0;
	ctx->CDF2=CDF2;
	ctx->list=list;
}
inline void ac_enc(ACEncContext *ctx, int sym)
{
	unsigned mingap, lo2, hi2;
	int cdf_curr, cdf_next;

	for(;;)//renorm
	{
		mingap=(ctx->hi-ctx->lo)>>16;
		if(ctx->lo<ctx->hi&&mingap>0)
			break;
		if(ctx->nbits>=32)
		{
			dlist_push_back(ctx->list, &ctx->cache, 4);
			ctx->cache=0;
			ctx->nbits=0;
		}
		ctx->cache|=(ctx->lo&0x80000000)>>ctx->nbits;//cache is written MSB -> LSB
		++ctx->nbits;

		ctx->lo<<=1;//shift out MSB
		ctx->hi<<=1;
		ctx->hi|=1;
	}

	cdf_curr=ctx->CDF2[sym];
	cdf_next=ctx->CDF2[sym+1];

	if(cdf_curr>=cdf_next)
		LOG_ERROR2("ZPS");

	lo2=ctx->lo+((unsigned long long)(ctx->hi-ctx->lo)*cdf_curr>>16);
	hi2=ctx->lo+((unsigned long long)(ctx->hi-ctx->lo)*cdf_next>>16);
	acval_enc(sym, cdf_curr, cdf_next, ctx->lo, ctx->hi, lo2, hi2, ctx->cache, ctx->nbits);//
	ctx->lo=lo2;
	ctx->hi=hi2-1;//OBLIGATORY range leak guard
}
inline void ac_enc_flush(ACEncContext *ctx)
{
	int k2=0;
	do//flush
	{
		while(ctx->nbits<32)
		{
			ctx->cache|=(ctx->lo&0x80000000)>>ctx->nbits;//cache is written MSB -> LSB
			++ctx->nbits;
			++k2;

			ctx->lo<<=1;//shift out MSB
			ctx->hi<<=1;
			ctx->hi|=1;
		}
		dlist_push_back(ctx->list, &ctx->cache, 4);
		ctx->cache=0;
		ctx->nbits=0;
	}while(k2<32);
}
typedef struct ACDecContextStruct
{
	unsigned lo, hi, cache;
	int nbits;
	unsigned code;
	const unsigned *CDF2;
	const unsigned char *srcptr, *srcend;
} ACDecContext;
inline void ac_dec_init(ACDecContext *ctx, const unsigned *CDF2, const unsigned char *start, unsigned const char *end)
{
	ctx->lo=0;
	ctx->hi=0xFFFFFFFF;
	ctx->nbits=32;
	ctx->CDF2=CDF2;
	ctx->srcptr=start;
	ctx->srcend=end;

	if(ctx->srcend-ctx->srcptr<4)
		LOG_ERROR2("buffer overflow");
	memcpy(&ctx->code, ctx->srcptr, 4);
	ctx->srcptr+=4;

	if(ctx->srcend-ctx->srcptr<4)
		LOG_ERROR2("buffer overflow");
	memcpy(&ctx->cache, ctx->srcptr, 4);
	ctx->srcptr+=4;
}
inline int ac_dec(ACDecContext *ctx)
{
	unsigned mingap;
	for(;;)//renorm
	{
		mingap=(ctx->hi-ctx->lo)>>16;
		if(ctx->lo<ctx->hi&&mingap>0)
			break;
		if(!ctx->nbits)
		{
			if(ctx->srcend-ctx->srcptr<4)
			{
#ifdef AC_VALIDATE
				printf("buffer overflow\n");
				acval_dump();
#endif
				LOG_ERROR2("buffer overflow");
			}
			memcpy(&ctx->cache, ctx->srcptr, 4);
			ctx->srcptr+=4;

			ctx->nbits=32;
		}
		--ctx->nbits;
		ctx->code<<=1;//shift out MSB		cache is read MSB -> LSB
		ctx->code|=(unsigned)(ctx->cache>>ctx->nbits&1);

		ctx->lo<<=1;
		ctx->hi<<=1;
		ctx->hi|=1;
	}

	unsigned lo2, hi2;
	int sym=0;
	int L=0, R=256, found=0;
	unsigned code2;
	while(L<=R)
	{
		sym=(L+R)>>1;
		code2=ctx->lo+((unsigned long long)(ctx->hi-ctx->lo)*ctx->CDF2[sym]>>16);
		if(code2<ctx->code)
			L=sym+1;
		else if(code2>ctx->code)
			R=sym-1;
		else
		{
			found=1;
			break;
		}
	}
	if(found)
		for(;sym<256-1&&ctx->lo+((unsigned long long)(ctx->hi-ctx->lo)*ctx->CDF2[sym+1]>>16)==ctx->code;++sym);
	else
		sym=L+(L<256&&ctx->lo+((unsigned long long)(ctx->hi-ctx->lo)*ctx->CDF2[sym+1]>>16)<ctx->code)-(L!=0);

	//buf[(iw*ky+kx2)<<2|kc]=(unsigned char)sym;

	unsigned cdf_start=ctx->CDF2[sym], cdf_end=ctx->CDF2[sym+1];
					
	lo2=ctx->lo+((unsigned long long)(ctx->hi-ctx->lo)*cdf_start>>16);
	hi2=ctx->lo+((unsigned long long)(ctx->hi-ctx->lo)*cdf_end  >>16);
	acval_dec(sym, cdf_start, cdf_end, ctx->lo, ctx->hi, lo2, hi2, ctx->cache, ctx->nbits, ctx->code);//
	ctx->lo=lo2;
	ctx->hi=hi2-1;//OBLIGATORY range leak guard
	return sym;
}


//asymmetric numeral systems coder
typedef struct ANSEncContextStruct
{
	unsigned state;
	const unsigned *CDF2;
	DList *list;
} ANSEncContext;
inline void ans_enc_init(ANSEncContext *ctx, const unsigned *CDF2, DList *list)
{
	ctx->state=0x10000;
	ctx->CDF2=CDF2;
	ctx->list=list;
}
inline void ans_enc(ANSEncContext *ctx, int sym, int kc)
{
	int cdf=ctx->CDF2[sym], freq=ctx->CDF2[sym+1]-cdf;
	if(ctx->state>=(unsigned)(freq<<16))//renorm
	{
		dlist_push_back(ctx->list, &ctx->state, 2);
		ctx->state>>=16;
	}
	debug_enc_update(ctx->state, cdf, freq, 0, 0, 0, kc, sym);
	ctx->state=ctx->state/freq<<16|(cdf+ctx->state%freq);//update
}
inline void ans_enc_flush(ANSEncContext *ctx)
{
	dlist_push_back(ctx->list, &ctx->state, 4);
}
typedef struct ANSDecContextStruct
{
	unsigned state;
	const unsigned *CDF2;
	const unsigned char *srcptr, *srcstart;
} ANSDecContext;
inline void ans_dec_init(ANSDecContext *ctx, const unsigned *CDF2, const unsigned char *start, const unsigned char *end)
{
	ctx->CDF2=CDF2;
	ctx->srcptr=end;
	ctx->srcstart=start;
	
	ctx->srcptr-=4;
	if(ctx->srcptr<ctx->srcstart)
		LOG_ERROR2("ANS buffer overflow");
	memcpy(&ctx->state, ctx->srcptr, 4);
}
inline int ans_dec(ANSDecContext *ctx, int kc)
{
	unsigned c=(unsigned short)ctx->state;
	int sym=0;

	int L=0, R=256, found=0;
	while(L<=R)
	{
		sym=(L+R)>>1;
		if(ctx->CDF2[sym]<c)
			L=sym+1;
		else if(ctx->CDF2[sym]>c)
			R=sym-1;
		else
		{
			found=1;
			break;
		}
	}
	if(found)
		for(;sym<256-1&&ctx->CDF2[sym+1]==c;++sym);
	else
		sym=L+(L<256&&ctx->CDF2[L]<c)-(L!=0);

	//buf[(iw*ky+kx)<<2|kc]=(unsigned char)sym;

	unsigned cdf=ctx->CDF2[sym], freq=ctx->CDF2[sym+1]-cdf;
						
	debug_dec_update(ctx->state, cdf, freq, 0, 0, 0, kc, sym);
	ctx->state=freq*(ctx->state>>16)+c-cdf;//update
	if(ctx->state<0x10000)//renorm
	{
		ctx->state<<=16;
		if(ctx->srcptr-2>=ctx->srcstart)
		{
			ctx->srcptr-=2;
			memcpy(&ctx->state, ctx->srcptr, 2);
		}
	}
	return sym;
}


#ifdef __cplusplus
}
#endif
#endif
