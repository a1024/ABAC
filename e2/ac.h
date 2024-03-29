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
void    acval_enc(int sym, int cdf, int freq, unsigned lo1, unsigned hi1, unsigned lo2, unsigned hi2, unsigned long long cache, int nbits);
void    acval_dec(int sym, int cdf, int freq, unsigned lo1, unsigned hi1, unsigned lo2, unsigned hi2, unsigned long long cache, int nbits, unsigned code);
void    acval_dump();
#ifdef AC_IMPLEMENTATION
typedef struct ACVALStruct
{
	int sym, cdf, freq;
	unsigned
		lo1, hi1,
		lo2, hi2;
	unsigned long long cache;
	int nbits;
	unsigned code;
} ACVAL;
#define ACVAL_CBUFSIZE 128
ArrayHandle acval=0;
ACVAL acval_cbuf[ACVAL_CBUFSIZE]={0};
int acval_idx=0;
void acval_enc(int sym, int cdf, int freq, unsigned lo1, unsigned hi1, unsigned lo2, unsigned hi2, unsigned long long cache, int nbits)
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
	printf("%4d sym %02X cdf %04X freq %04X range %08X~%08X->%08X~%08X cache %08X %08X nbits %2d",
		idx, p->sym, p->cdf, p->freq, p->lo1, p->hi1, p->lo2, p->hi2, (int)(p->cache>>32), (int)p->cache, p->nbits
	);
	//printf("%4d sym 0x%02X cdf 0x%04X freq 0x%04X range 0x%08X~0x%08X->0x%08X~0x%08X cache 0x%016llX nbits %2d",
	//	idx, p->sym, p->cdf, p->freq, p->lo1, p->hi1, p->lo2, p->hi2, p->cache, p->nbits
	//);
	if(dec)
	{
		printf(" code %08X ", p->code);
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
		if(k2<0)
			continue;
		if(k==end)
			printf("\n");
		acval_print(k2, acval_cbuf+k, 1);
		if(k==end)
			break;
	}
}
void acval_dec(int sym, int cdf, int freq, unsigned lo1, unsigned hi1, unsigned lo2, unsigned hi2, unsigned long long cache, int nbits, unsigned code)
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
	if(sym!=val->sym||cdf!=val->cdf||freq!=val->freq||lo1!=val->lo1||hi1!=val->hi1||lo2!=val->lo2||hi1!=val->hi1)
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


//arithmetic coder
typedef struct ArithmeticCoderStruct
{
	unsigned lo, hi;
	unsigned long long cache;
	int nbits;//enc: number of free bits in cache [0 ~ 64], dec: number of unread bits in cache [0 ~ 32]
	int is_enc;
	union
	{
		struct
		{
			const unsigned char *srcptr, *srcend;
			unsigned code;
		};
#ifdef EC_USE_ARRAY
		ArrayHandle *arr;
#else
		DList *list;
#endif
	};
} ArithmeticCoder;
static void ac_enc_init(ArithmeticCoder *ec,
#ifdef EC_USE_ARRAY
	ArrayHandle *arr
#else
	DList *list
#endif
)
{
	ec->lo=0;
	ec->hi=0xFFFFFFFF;
	ec->cache=0;
	ec->nbits=64;
	ec->is_enc=1;
#ifdef EC_USE_ARRAY
	array_append(arr, 0, 1, 0, 0, 0, 0);
	ec->arr=arr;
#else
	ec->list=list;
#endif
}
static void ac_dec_init(ArithmeticCoder *ec, const unsigned char *start, unsigned const char *end)
{
	ec->lo=0;
	ec->hi=0xFFFFFFFF;
	ec->cache=0;
	ec->nbits=0;
	ec->is_enc=0;
	ec->srcptr=start;
	ec->srcend=end;

	if(ec->srcptr+4>ec->srcend)
	{
		LOG_ERROR2("buffer overflow");
		return;
	}
	memcpy(&ec->code, ec->srcptr, 4);
	ec->srcptr+=4;
}
static void ac_renorm(ArithmeticCoder *ec, unsigned lg_fmin)//one-time loopless renorm		this keeps hi & lo as far apart as possible from each other in the ALU
{
	//((hi-lo)<<n_emit)*fmin/0x10000 should be >= 1, where fmin is the current smallest nonzero freq
	//floor_log2(hi-lo)+n_emit+floor_log2(fmin)-16 >= 1
	//32-n_keep = n_emit >= 16-floor_log2(fmin)-floor_log2(hi-lo)
	//n_keep <= 16+floor_log2(fmin)+floor_log2(hi-lo)
	
	int n_keep=floor_log2_32(ec->hi^ec->lo)+1, n_emit;
	int range_guard=floor_log2_32(ec->lo<ec->hi?ec->hi-ec->lo:1)+lg_fmin+16;
	if(n_keep>range_guard)//a solution for the rare "Schrodinger's half" problem
	{
		unsigned lo2, hi2;
		n_emit=32-range_guard;
		lo2=ec->lo<<n_emit;
		hi2=ec->hi<<n_emit|((1<<n_emit)-1);
		while(hi2<lo2)
		{
			lo2<<=1;
			hi2<<=1;
			hi2|=1;
			++n_emit;
		}
		n_keep=32-n_emit;
	}
	else
		n_emit=32-n_keep;
	if(n_emit)
	{
		unsigned emit_mask=(unsigned)((1LL<<n_emit)-1);
		ec->nbits-=n_emit;//new nbits
		if(ec->is_enc)//encoding
		{
			//buffer: {b,a,a,a, c,c,c,b, e,e,d,c, f,f,f,e}, cache: MSB gg[hhh]000 LSB	nbits 6->3, h is about to be emitted
			//written 32-bit words are reversed because ACs are inherently big-endian (significance decreases during coding)

			if(ec->nbits<0)//number of free bits in cache will be negative, flush 32 MSBs to buffer
			{
#ifdef EC_USE_ARRAY
				ARRAY_APPEND(*ec->arr, (unsigned*)&ec->cache+1, 4, 1, 0);
#else
				dlist_push_back(ec->list, (unsigned*)&ec->cache+1, 4);
#endif
				ec->nbits+=32;
				ec->cache<<=32;
			}
			ec->cache|=(unsigned long long)(ec->lo>>n_keep)<<ec->nbits;
		}
		else//decoding
		{
			//same operations triger same renorms                            v srcptr
			//buffer: {b,a,a,a, c,c,c,b, e,e,d,c, f,f,f,e, h,h,g,g, j,j,i,h, | l,k,k,j, n,n,m,l, q,p,o,n, ...}
			//cache: MSB gg[hhh]ijj LSB		nbits 6->3, h is about to be read (g is now useless old code)

			if(ec->nbits<0)//number of bits in cache will be negative, insert 32 LSBs from buffer
			{
				if(ec->srcptr+4>ec->srcend)
				{
#ifdef AC_VALIDATE
					printf("buffer overflow\n");
					acval_dump();
#endif
					LOG_ERROR2("buffer overflow");
					return;
				}
				ec->nbits+=32;
				ec->cache<<=32;

				memcpy(&ec->cache, ec->srcptr, 4);
				ec->srcptr+=4;
			}
			if(n_emit==32)
				ec->code=(unsigned)(ec->cache>>ec->nbits);
			else
			{
				ec->code<<=n_emit;
				ec->code|=(unsigned)(ec->cache>>ec->nbits&emit_mask);
			}
		}
		if(n_emit==32)//shift left uint32 by 32 is UB, shift ammount is 5 bit (0 ~ 31)
		{
			ec->lo=0;
			ec->hi=0xFFFFFFFF;
		}
		else
		{
			ec->lo<<=n_emit;
			ec->hi<<=n_emit;
			ec->hi|=emit_mask;
		}
	}
}
static void ac_enc_flush(ArithmeticCoder *ec)
{
	ec->hi=ec->lo;//this will cause all remaining 32 lo bits to be written to the cache
	ac_renorm(ec, 0x10000);
	//now all remaining bits are in the cache (up to 64)
	while(64-ec->nbits>0)//loops up to 2 times
	{
#ifdef EC_USE_ARRAY
		ARRAY_APPEND(*ec->arr, (unsigned*)&ec->cache+1, 4, 1, 0);
#else
		dlist_push_back(ec->list, (unsigned*)&ec->cache+1, 4);
#endif
		ec->nbits+=32;
		ec->cache<<=32;
	}
}

static void ac_enc(ArithmeticCoder *ec, int sym, const unsigned *CDF, int nlevels, unsigned lg_fmin)//CDF is 16 bit
{
	unsigned lo2, hi2;
	int cdf_curr, cdf_next;
	
	ac_renorm(ec, lg_fmin);
	if(CDF)
	{
		cdf_curr=CDF[sym];
		cdf_next=CDF[sym+1];
	}
	else//bypass
	{
		cdf_curr=(sym<<16)/nlevels;
		cdf_next=((sym+1)<<16)/nlevels;
	}

	if(cdf_curr>=cdf_next)
		LOG_ERROR2("ZPS");

	lo2=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*cdf_curr>>16);
	hi2=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*cdf_next>>16);
	acval_enc(sym, cdf_curr, cdf_next-cdf_curr, ec->lo, ec->hi, lo2, hi2, ec->cache, ec->nbits);//
	ec->lo=lo2;
	ec->hi=hi2-1;//must decrement hi because decoder fails when code == hi2
}
static int ac_dec(ArithmeticCoder *ec, const unsigned *CDF, int nlevels, unsigned lg_fmin)
{
	unsigned lo2, hi2;
	int sym=0;
	unsigned cdf_start, cdf_end;
	
	ac_renorm(ec, lg_fmin);
	if(CDF)
	{
		int L=0, R=nlevels-1, found=0;
		unsigned code2;
		while(L<=R)//binary search
		{
			sym=(L+R)>>1;
			code2=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*CDF[sym]>>16);
			if(code2<ec->code)
				L=sym+1;
			else if(code2>ec->code)
				R=sym-1;
			else
			{
				found=1;
				break;
			}
		}
		if(found)
			for(;sym<nlevels-1&&ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*CDF[sym+1]>>16)==ec->code;++sym);
		else
			sym=L+(L<nlevels&&ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*CDF[sym+1]>>16)<ec->code)-(L!=0);
		cdf_start=CDF[sym];
		cdf_end=CDF[sym+1];
	}
	else//bypass
	{
		cdf_start=(int)(((unsigned long long)(ec->code-ec->lo)<<16)/(ec->hi-ec->lo));
		sym=(cdf_start*nlevels>>16)+1;//this is to handle the case when code == lo2
		cdf_start=(sym<<16)/nlevels;
		lo2=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*cdf_start>>16);
		if(lo2>ec->code)
		{
			--sym;
			cdf_start=(sym<<16)/nlevels;
		}
		cdf_end=((sym+1)<<16)/nlevels;
	}
	if(cdf_start>=cdf_end)
		LOG_ERROR2("ZPS");

	lo2=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*cdf_start>>16);
	hi2=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*cdf_end  >>16);
	acval_dec(sym, cdf_start, cdf_end-cdf_start, ec->lo, ec->hi, lo2, hi2, ec->cache, ec->nbits, ec->code);//
	ec->lo=lo2;
	ec->hi=hi2-1;//must decrement hi because decoder fails when code == hi2
	return sym;
}

static void ac_enc15(ArithmeticCoder *ec, int sym, const unsigned short *CDF, int nlevels)//CDF is 15 bit
{
	unsigned lo2, hi2;
	int cdf_curr, cdf_next;
	
	ac_renorm(ec, 0);
	if(CDF)
	{
		cdf_curr=CDF[sym]<<1;
		cdf_next=CDF[sym+1]<<1;
	}
	else//bypass
	{
		cdf_curr=(sym<<16)/nlevels;
		cdf_next=((sym+1)<<16)/nlevels;
	}

	if(cdf_curr>=cdf_next)
		LOG_ERROR2("ZPS");

	lo2=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*cdf_curr>>16);
	hi2=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*cdf_next>>16);
	acval_enc(sym, cdf_curr, cdf_next-cdf_curr, ec->lo, ec->hi, lo2, hi2, ec->cache, ec->nbits);//
	ec->lo=lo2;
	ec->hi=hi2-1;//must decrement hi because decoder fails when code == hi2
}
static int ac_dec15(ArithmeticCoder *ec, const unsigned short *CDF, int nlevels)
{
	unsigned lo2, hi2;
	int sym=0;
	unsigned cdf_start, cdf_end;
	
	ac_renorm(ec, 0);
	if(CDF)
	{
		int L=0, R=nlevels-1, found=0;
		unsigned code2;
		while(L<=R)//binary search
		{
			sym=(L+R)>>1;
			code2=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*(CDF[sym]<<1)>>16);
			if(code2<ec->code)
				L=sym+1;
			else if(code2>ec->code)
				R=sym-1;
			else
			{
				found=1;
				break;
			}
		}
		if(found)
			for(;sym<nlevels-1&&ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*(CDF[sym+1]<<1)>>16)==ec->code;++sym);
		else
			sym=L+(L<nlevels&&ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*(CDF[sym+1]<<1)>>16)<ec->code)-(L!=0);
		cdf_start=CDF[sym]<<1;
		cdf_end=CDF[sym+1]<<1;
	}
	else//bypass
	{
		cdf_start=(int)(((unsigned long long)(ec->code-ec->lo)<<16)/(ec->hi-ec->lo));
		sym=(cdf_start*nlevels>>16)+1;//this is to handle the case when code == lo2
		cdf_start=(sym<<16)/nlevels;
		lo2=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*cdf_start>>16);
		if(lo2>ec->code)
		{
			--sym;
			cdf_start=(sym<<16)/nlevels;
		}
		cdf_end=((sym+1)<<16)/nlevels;
	}
	if(cdf_start>=cdf_end)
		LOG_ERROR2("ZPS");

	lo2=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*cdf_start>>16);
	hi2=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*cdf_end  >>16);
	acval_dec(sym, cdf_start, cdf_end-cdf_start, ec->lo, ec->hi, lo2, hi2, ec->cache, ec->nbits, ec->code);//
	ec->lo=lo2;
	ec->hi=hi2-1;//must decrement hi because decoder fails when code == hi2
	return sym;
}

static void ac_enc_bin(ArithmeticCoder *ec, unsigned short p0, int bit)
{
	unsigned mid;
	
	ac_renorm(ec, 0);

	if(!p0)//reject degenerate distribution
		LOG_ERROR2("ZPS");

	mid=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*p0>>16);

	acval_enc(bit, bit?p0:0, bit?0x10000-p0:p0, ec->lo, ec->hi, bit?mid:ec->lo, bit?ec->hi:mid-1, ec->cache, ec->nbits);

	if(bit)
		ec->lo=mid;
	else
		ec->hi=mid-1;//must decrement hi because decoder fails when code == hi2
}
static int ac_dec_bin(ArithmeticCoder *ec, unsigned short p0)//binary AC decoder doesn't do binary search
{
	unsigned mid;

	ac_renorm(ec, 0);

	if(!p0)//reject degenerate distribution
		LOG_ERROR2("ZPS");

	mid=ec->lo+(int)((unsigned long long)(ec->hi-ec->lo)*p0>>16);
	int bit=ec->code>=mid;

	acval_dec(bit, bit?p0:0, bit?0x10000-p0:p0, ec->lo, ec->hi, bit?mid:ec->lo, bit?ec->hi:mid-1, ec->cache, ec->nbits, ec->code);
	
	if(bit)
		ec->lo=mid;
	else
		ec->hi=mid-1;//must decrement hi because decoder fails when code == hi2

	return bit;
}


//asymmetric numeral systems coder
typedef struct ANSCoderStruct
{
	unsigned state;
	union
	{
		struct
		{
			const unsigned char *srcptr, *srcstart;
		};
#ifdef EC_USE_ARRAY
		ArrayHandle *arr;
#else
		DList *list;
#endif
	};
} ANSCoder;
static void ans_enc_init(ANSCoder *ec,
#ifdef EC_USE_ARRAY
	ArrayHandle *arr
#else
	DList *list
#endif
)
{
	ec->state=0x10000;
#ifdef EC_USE_ARRAY
	ec->arr=arr;
#else
	ec->list=list;
#endif
}
static void ans_dec_init(ANSCoder *ec, const unsigned char *start, const unsigned char *end)
{
	ec->srcptr=end;
	ec->srcstart=start;

	ec->srcptr-=4;
	if(ec->srcptr<ec->srcstart)
		LOG_ERROR2("ANS buffer overflow");
	memcpy(&ec->state, ec->srcptr, 4);
}
static void ans_enc_flush(ANSCoder *ec)
{
#ifdef EC_USE_ARRAY
	ARRAY_APPEND(*ec->arr, &ec->state, 4, 1, 0);
#else
	dlist_push_back(ec->list, &ec->state, 4);
#endif
}

static void ans_enc(ANSCoder *ec, int sym, const unsigned *CDF, int nlevels)
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
static int ans_dec(ANSCoder *ec, const unsigned *CDF, int nlevels)
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

static void ans_enc15(ANSCoder *ec, int sym, const unsigned short *CDF, int nlevels)
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
static int ans_dec15(ANSCoder *ec, const unsigned short *CDF, int nlevels)
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

static void ans_enc_bin(ANSCoder *ec, unsigned short p0, int bit)
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
static int ans_dec_bin(ANSCoder *ec, unsigned short p0)
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


#ifdef __cplusplus
}
#endif
#endif
