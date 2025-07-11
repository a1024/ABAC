#pragma once
#ifndef INC_AC_H
#define INC_AC_H
#include"util.h"
#include<stdio.h>
#include<string.h>
#include<immintrin.h>
//#ifdef _MSC_VER
//#include<intrin.h>//_udiv128, _umul128 -> _mulx_u64
//#endif
#include"blist.h"
#ifdef __cplusplus
extern "C"
{
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
		printf(" code %016llX ", p->code);
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

#ifdef ANS_VAL
#define ANS_VAL_HISTSIZE 128
typedef struct _ANSVALHeader
{
	unsigned short esize, count;
	unsigned idx;
} ANSVALHeader;
static ArrayHandle debugstack=0;
static void ansval_push(const void *data, int esize, int count)
{
	static int idx=0;
	ANSVALHeader header={esize, count, idx};
	++idx;
	if(!debugstack)
		ARRAY_ALLOC(char, debugstack, 0, 0, 1024, 0);
	ARRAY_APPEND(debugstack, &header, sizeof(header), 1, 0);//lo header
	ARRAY_APPEND(debugstack, data, (ptrdiff_t)count*esize, 1, 0);
	ARRAY_APPEND(debugstack, &header, sizeof(header), 1, 0);//hi header
}
static void ansval_printr(const void *data, int esize, int count, const void *xdata)//print elements in reverse because little-endian
{
	const unsigned char *p=(const unsigned char*)data;
	const unsigned char *p2=(const unsigned char*)xdata;
	int size=count*esize;
	for(int k=0;k<size;k+=esize)
	{
		printf(" ");
		for(int k2=esize-1;k2>=0;--k2)
		{
			int val=p[k+k2];
			if(p2)
				val^=p2[k+k2];
			if(p2&&!val)
				printf("--");
			else
				printf("%02X", val);
		}
	}
	printf("\n");
	//for(int k=size-1;k>=0;--k)
	//{
	//	if((k&3)==3)
	//		printf(" ");
	//	printf("%02X", ((unsigned char*)data)[k]);
	//}
	//printf("\n");
}
static void* ansval_ptrguard(const void *start, const void *end, const void *ptr, ptrdiff_t nbytes)
{
	size_t istart=(size_t)start, iend=(size_t)end;
	ptrdiff_t size=iend-istart;
	size_t ip1=(size_t)ptr, ip2=ip1+nbytes;
	int problems[]=
	{
		size<0,
		(size_t)(ip1-istart)>=(size_t)size,
		(size_t)(ip2-istart)>=(size_t)size,
	};
	if(problems[0]||problems[1]||problems[2])
	{
		printf("\nOOB\n");
		printf("  inc     %+16td bytes\n", nbytes);
		printf("  start   %016zd  %16d\n", istart, 0);
		if(nbytes<0)
		{
			printf("  after   %016zd  %16td%s\n", ip2, ip2-istart, problems[2]?"  <-":"");
			printf("  before  %016zd  %16td%s\n", ip1, ip1-istart, problems[1]?"  <-":"");
		}
		else
		{
			printf("  before  %016zd  %16td%s\n", ip1, ip1-istart, problems[1]?"  <-":"");
			printf("  after   %016zd  %16td%s\n", ip2, ip2-istart, problems[2]?"  <-":"");
		}
		printf("  end     %016zd  %16td%s\n", iend, size, problems[0]?"  <-":"");
		LOG_ERROR("\n");
		return 0;
	}
	return (void*)(nbytes<0?ip2:ip1);
}
static void ansval_check(const void *data, int esize, int count)
{
	static const unsigned char *endptr=0;
	static int totalcount=0, popcount=0;
	int firstpop=!endptr;
	if(firstpop)
		endptr=debugstack->data+debugstack->count;
	const ANSVALHeader *hiheader=(const ANSVALHeader*)(debugstack->data+debugstack->count)-1;
	debugstack->count-=sizeof(ANSVALHeader);
	const unsigned char *data0=debugstack->data+debugstack->count-hiheader->count*hiheader->esize;
	debugstack->count-=hiheader->count*hiheader->esize;
	const ANSVALHeader *loheader=(const ANSVALHeader*)(debugstack->data+debugstack->count)-1;
	debugstack->count-=sizeof(ANSVALHeader);
	if(firstpop)
		totalcount=hiheader->idx+1;
	if(memcmp(hiheader, loheader, sizeof(*hiheader)))
	{
		printf("\n");
		printf("Validation Header Mismatch  idx,esize,count: loheader %d,%d,%d != hiheader %d,%d,%d\n",
			loheader->idx, loheader->esize, loheader->count,
			hiheader->idx, hiheader->esize, hiheader->count
		);
		LOG_ERROR("\n");
	}
	if(esize!=loheader->esize||count!=loheader->count||memcmp(data, data0, count*esize))
	{
		printf("\n");
		printf("Validation Error    pop #%d / total %d,  remaining %d,  using %8.2lf/%8.2lf MB\n",
			popcount,
			totalcount,
			totalcount-popcount-1,
			(double)debugstack->count/(1024*1024),
			(double)debugstack->cap/(1024*1024)
		);
		printf("\n");

		const unsigned char *verptr=debugstack->data+debugstack->count+loheader->count*loheader->esize+sizeof(ANSVALHeader[2]);
		const unsigned char *ptrstack[ANS_VAL_HISTSIZE]={0};
		const ANSVALHeader *loheader2=0, *hiheader2=0;
		const unsigned char *verdata=0, *unverdata=0;
		int nptrs=0;
		printf("Verified pops:\n");
		for(int k=0;k<ANS_VAL_HISTSIZE;++k)
		{
			if(verptr>=endptr)
				break;

			loheader2=(const ANSVALHeader*)ansval_ptrguard(debugstack->data, endptr, verptr, +sizeof(ANSVALHeader));
			verptr+=sizeof(ANSVALHeader);
			verdata=(const unsigned char*)ansval_ptrguard(debugstack->data, endptr, verptr, +loheader2->count*loheader2->esize);
			verptr+=loheader2->count*loheader2->esize;
			hiheader2=(const ANSVALHeader*)verptr;
			verptr+=sizeof(ANSVALHeader);

			ptrstack[nptrs++]=(const unsigned char*)loheader2;
			(void)verdata;
			(void)hiheader2;
		}
		for(int k=nptrs-1;k>=0;--k)
		{
			const unsigned char *ptr=ptrstack[k];
			loheader2=(const ANSVALHeader*)ptr;
			ptr+=sizeof(ANSVALHeader);
			verdata=ptr;
			ptr+=loheader2->count*loheader2->esize;
			hiheader2=(const ANSVALHeader*)ptr;
			ptr+=sizeof(ANSVALHeader);

			printf("  [%7d] %7d B    ", loheader2->idx, loheader2->count*loheader2->esize);
			ansval_printr(verdata, loheader2->esize, loheader2->count, 0);
			(void)hiheader2;
		}
		if(!nptrs)
			printf("  No data\n");
		printf("\n");

		printf("The error:\n");
		printf("  [%7d] Original %7d B    ", loheader->idx, loheader->count*loheader->esize);
		ansval_printr(data0, loheader->esize, loheader->count, 0);
		printf("  [%7d] Corrupt  %7d B    ", loheader->idx, count*esize);
		ansval_printr(data, esize, count, 0);
		printf("  [%7d] XOR      %7s      ", loheader->idx, "");
		ansval_printr(data0, esize, count, data);
		printf("\n");
		
		const unsigned char *unverptr=debugstack->data+debugstack->count;
		printf("Remaining pops:\n");
		for(int k=0;k<ANS_VAL_HISTSIZE;++k)
		{
			if(unverptr<=debugstack->data)
			{
				printf("  No data\n");
				break;
			}

			unverptr=(const unsigned char*)ansval_ptrguard(debugstack->data, endptr, unverptr, -(ptrdiff_t)sizeof(ANSVALHeader));
			hiheader2=(const ANSVALHeader*)unverptr;
			unverptr=(const unsigned char*)ansval_ptrguard(debugstack->data, endptr, unverptr, -(ptrdiff_t)hiheader2->count*hiheader2->esize);
			unverdata=unverptr;
			unverptr=(const unsigned char*)ansval_ptrguard(debugstack->data, endptr, unverptr, -(ptrdiff_t)sizeof(ANSVALHeader));
			loheader2=(const ANSVALHeader*)unverptr;
			
			if(memcmp(hiheader2, loheader2, sizeof(*hiheader2)))
			{
				printf("  Header Mismatch:  (idx,esize,count)\n");
				printf("    hiheader  0x%08X, 0x%08X, 0x%08X    %d, %d, %d\n",
					hiheader2->idx, hiheader2->esize, hiheader2->count,
					hiheader2->idx, hiheader2->esize, hiheader2->count
				);
				printf("    loheader  0x%08X, 0x%08X, 0x%08X    %d, %d, %d\n",
					loheader2->idx, loheader2->esize, loheader2->count,
					loheader2->idx, loheader2->esize, loheader2->count
				);
				break;
			}
			printf("  [%7d] %7d B    ", hiheader2->idx, hiheader2->count*hiheader2->esize);
			ansval_printr(unverdata, hiheader2->esize, hiheader2->count, 0);
		}
		printf("\n");
		LOG_ERROR("");
	}
	++popcount;
}
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
AWM_INLINE void ac2_bitwrite(AC2 *ec, int bit)
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
AWM_INLINE void ac2_bitwrite_withpending(AC2 *ec, int bit)
{
	ac2_bitwrite(ec, bit);
	bit^=1;
	while(ec->pending_bits)
	{
		ac2_bitwrite(ec, bit);
		--ec->pending_bits;
	}
}
AWM_INLINE int ac2_bitread(AC2 *ec)
{
	if(!ec->nbits)
	{
		ec->cache=ec->srcptr<ec->srcend?*ec->srcptr++:0;
		ec->nbits=8;
	}
	--ec->nbits;
	return ec->cache>>ec->nbits&1;
}
AWM_INLINE void ac2_enc_init(AC2 *ec, BList *dst)
{
	memset(ec, 0, sizeof(*ec));
	ec->x1=0;
	ec->x2=~0;
	ec->dst=dst;
}
AWM_INLINE void ac2_dec_init(AC2 *ec, const unsigned char *start, const unsigned char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->x1=0;
	ec->x2=~0;
	ec->srcptr=start;
	ec->srcend=end;
	for(int k=0;k<32;++k)
		ec->code=ec->code<<1|ac2_bitread(ec);
}
AWM_INLINE void ac2_enc_flush(AC2 *ec)
{
	do
	{
		ac2_bitwrite_withpending(ec, ec->x1>>31);
		ec->x1<<=1;
	}while(ec->nbits||ec->x1);//while(ec->nbits);
}
AWM_INLINE void ac2_enc_renorm(AC2 *ec)
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
AWM_INLINE void ac2_dec_renorm(AC2 *ec)
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

AWM_INLINE void ac2_enc_update(AC2 *ec, unsigned cdf_curr, unsigned cdf_next)
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
AWM_INLINE unsigned ac2_dec_getcdf(AC2 *ec)
{
	unsigned cdf;

	ac2_dec_renorm(ec);
	cdf=(unsigned)(((unsigned long long)(ec->code-ec->x1)<<AC2_PROB_BITS|((1LL<<AC2_PROB_BITS)-1))/(ec->x2-ec->x1));
	return cdf;
}
AWM_INLINE void ac2_dec_update(AC2 *ec, unsigned cdf_curr, unsigned cdf_next)
{
	unsigned range, x1, x2;

	range=ec->x2-ec->x1;
	x1=ec->x1+(unsigned)((unsigned long long)range*cdf_curr>>AC2_PROB_BITS);
	x2=ec->x1+(unsigned)((unsigned long long)range*cdf_next>>AC2_PROB_BITS)-1;//must decrement hi because decoder fails when code == hi2
	acval_dec(0, cdf_curr, cdf_next-cdf_curr, ec->x1, ec->x2, x1, x2, 0, 0, ec->code);
	ec->x1=x1;
	ec->x2=x2;
}
AWM_INLINE void ac2_enc_bypass(AC2 *ec, unsigned sym, int nbits)
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
AWM_INLINE unsigned ac2_dec_bypass(AC2 *ec, int nbits)
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
AWM_INLINE void ac2_enc_bypass_NPOT(AC2 *ec, unsigned sym, int nlevels)
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
AWM_INLINE unsigned ac2_dec_bypass_NPOT(AC2 *ec, int nlevels)
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
AWM_INLINE void ac2_enc_bin(AC2 *ec, unsigned p1, int bit)
{
	unsigned mid;

	ac2_enc_renorm(ec);
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
AWM_INLINE int ac2_dec_bin(AC2 *ec, unsigned p1)
{
	unsigned mid;
	int bit;

	ac2_dec_renorm(ec);
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
	return bit;
}


//arithmetic coder (carry-less)
#define AC3_PROB_BITS 16
#define AC3_RENORM 32	//multiple of 8!	//32-bit renorm is best
#if AC3_RENORM>AC3_PROB_BITS
#define AC3_RENORM_STATEMENT if
#else
#define AC3_RENORM_STATEMENT while
#endif
typedef struct _AC3
{
	unsigned long long low, range, code;
	BList *dst;
	unsigned char *dstptr, *dststart, *dstend;
	const unsigned char *srcptr, *srcstart, *srcend;
} AC3;
AWM_INLINE void ac3_encbuf_init(AC3 *ec, unsigned char *start, unsigned char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->low=0;
	ec->range=0xFFFFFFFFFFFF;
	ec->dstptr=start;
	ec->dststart=start;
	ec->dstend=end;
}
#if 0
AWM_INLINE void ac3_encbuf_renorm(AC3 *ec)
{
	unsigned long long rmax;

	*(unsigned*)ec->dstptr=(unsigned)(ec->low>>32);
	ec->dstptr+=4;
	ec->range=ec->range<<AC3_RENORM|((1LL<<AC3_RENORM)-1);
	ec->low<<=AC3_RENORM;

	rmax=~ec->low;
	if(ec->range>rmax)//clamp hi to register size after renorm
		ec->range=rmax;
}
#endif
AWM_INLINE unsigned char* ac3_encbuf_flush(AC3 *ec)
{
#ifndef AC_DISABLE_WRITE
	*(unsigned*)ec->dstptr=(unsigned)(ec->low>>32);
	ec->dstptr+=4;
	*(unsigned*)ec->dstptr=(unsigned)ec->low;
	ec->dstptr+=4;
#endif
	return ec->dstptr;
}
AWM_INLINE void ac3_encbuf_bypass_NPOT(AC3 *ec, int bypass, int nlevels)
{
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if((unsigned)bypass>=(unsigned)nlevels)
		LOG_ERROR("Bypass OOB");
#endif
	if(ec->range<(unsigned)nlevels)
	{
		unsigned long long rmax;

		*(unsigned*)ec->dstptr=(unsigned)(ec->low>>32);
		ec->dstptr+=4;
		ec->range=ec->range<<AC3_RENORM|((1LL<<AC3_RENORM)-1);
		ec->low<<=AC3_RENORM;

		rmax=~ec->low;
		if(ec->range>rmax)//clamp hi to register size after renorm
			ec->range=rmax;
	}
	ec->low+=ec->range*bypass/nlevels;
	ec->range=ec->range/nlevels-1;
	acval_enc(nlevels, bypass, 1, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
AWM_INLINE void ac3_encbuf_update_N(AC3 *ec, unsigned cdf, unsigned freq, int probbits)//probbits <= 16
{
	//if(ec->range<(1ULL<<16))
	//{
	//	*(unsigned*)ec->dstptr=(unsigned)(ec->low>>32);
	//	ec->dstptr+=4;
	//	ec->range=ec->range<<32|0xFFFFFFFF;
	//	ec->low<<=32;
	//
	//	if(ec->range>~ec->low)
	//		ec->range=~ec->low;
	//}
	//ec->low+=ec->range*cdf>>16;
	//ec->range=(ec->range*freq>>16)-1;

#ifdef _DEBUG
	if(cdf>=(unsigned)(1<<probbits)||cdf+freq>(unsigned)(1<<probbits))
		LOG_ERROR2("Invalid stats");
#endif
#ifdef AC3_ENC_BRANCHLESSRENORM
	unsigned emit=(unsigned)(ec->low>>32);
	unsigned char *nextptr=ec->dstptr+4;
	unsigned long long r2=ec->range<<AC3_RENORM|((1LL<<AC3_RENORM)-1);
	unsigned long long low2=ec->low<<AC3_RENORM;
	unsigned long long norenormflag=ec->range>>probbits;

	if(!norenormflag)			//FIXME this condition can be replaced with direct assignment
		*(unsigned*)ec->dstptr=emit;
	if(!norenormflag)			//FIXME ensure CMOVs
		ec->dstptr=nextptr;
	if(!norenormflag)
		ec->range=r2;
	if(!norenormflag)
		ec->low=low2;

	unsigned long long rmax=~ec->low;
	if(ec->range>rmax)
		ec->range=rmax;
#else
	if(ec->range<(1ULL<<probbits))
	{
#ifndef AC_DISABLE_WRITE
		*(unsigned*)ec->dstptr=(unsigned)(ec->low>>32);
		ec->dstptr+=4;
#endif
		ec->range=ec->range<<32|0xFFFFFFFF;
		ec->low<<=32;
		if(ec->range>~ec->low)//clamp hi to register size after renorm
			ec->range=~ec->low;
	}
#endif
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if((unsigned)(freq-1)>=(unsigned)((1<<probbits)-1)||cdf+freq<cdf)
		LOG_ERROR2("Invalid stats");
#endif
	ec->low+=ec->range*cdf>>probbits;
	ec->range=(ec->range*freq>>probbits)-1;
	acval_enc(probbits, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
AWM_INLINE void ac3_encbuf_bin(AC3 *ec, int bit, unsigned p0, int probbits)//probbits <= 16
{
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!p0||p0>=1ULL<<probbits)
		LOG_ERROR2("ZPS");
#endif
	if(ec->range<(1ULL<<probbits))
	{
#ifndef AC_DISABLE_WRITE
		*(unsigned*)ec->dstptr=(unsigned)(ec->low>>32);
		ec->dstptr+=4;
#endif
		ec->range=ec->range<<32|0xFFFFFFFF;
		ec->low<<=32;
		if(ec->range>~ec->low)//clamp hi to register size after renorm
			ec->range=~ec->low;
	}
	unsigned long long mid=ec->range*p0>>probbits;
	if(bit)
	{
		ec->low+=mid;
		ec->range-=mid;
	}
	else
		ec->range=mid-1;
	acval_enc(bit, bit?p0:0, bit?(1<<probbits)-p0:p0, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}

AWM_INLINE void ac3_enc_init(AC3 *ec, BList *dst)
{
	memset(ec, 0, sizeof(*ec));
	ec->low=0;
	ec->range=0xFFFFFFFFFFFF;
	ec->dst=dst;
}
AWM_INLINE void ac3_dec_init(AC3 *ec, const unsigned char *start, const unsigned char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->low=0;
	ec->range=0xFFFFFFFFFFFF;
	ec->srcptr=start;
	ec->srcstart=start;
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
AWM_INLINE void ac3_enc_renorm(AC3 *ec)
{
	unsigned long long rmax;

	blist_push_back(ec->dst, (unsigned char*)&ec->low+8-AC3_RENORM/8, AC3_RENORM/8);
	ec->range=ec->range<<AC3_RENORM|((1LL<<AC3_RENORM)-1);
	ec->low<<=AC3_RENORM;

	rmax=~ec->low;
	if(ec->range>rmax)//clamp hi to register size after renorm
		ec->range=rmax;
}
AWM_INLINE void ac3_dec_renorm(AC3 *ec)
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
AWM_INLINE void ac3_enc_flush(AC3 *ec)
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

AWM_INLINE void ac3_enc_update(AC3 *ec, unsigned cdf, unsigned freq)
{
	AC3_RENORM_STATEMENT(!(ec->range>>AC3_PROB_BITS))
		ac3_enc_renorm(ec);
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if((unsigned)(freq-1)>=(0x10000-1))
		LOG_ERROR2("ZPS");
	if(cdf>(1<<AC3_PROB_BITS)||cdf+freq<cdf)
		LOG_ERROR2("Invalid CDF");
#endif
	ec->low+=ec->range*cdf>>AC3_PROB_BITS;
	ec->range=(ec->range*freq>>AC3_PROB_BITS)-1;//must decrement hi because decoder fails when code == hi2
	acval_enc(0, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
AWM_INLINE unsigned ac3_dec_getcdf(AC3 *ec)
{
	AC3_RENORM_STATEMENT(!(ec->range>>AC3_PROB_BITS))
		ac3_dec_renorm(ec);
#ifdef _DEBUG
	if(ec->code<ec->low||ec->code>ec->low+ec->range)
		LOG_ERROR("CODE OOB");
#endif
	return (unsigned)(((ec->code-ec->low)<<AC3_PROB_BITS|((1ULL<<AC3_PROB_BITS)-1))/ec->range);
}
AWM_INLINE void ac3_dec_update(AC3 *ec, unsigned cdf, unsigned freq)
{
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if((unsigned)(freq-1)>=(0x10000-1))
		LOG_ERROR2("ZPS");
	if(cdf+freq<cdf)
		LOG_ERROR2("Invalid CDF");
#endif
	ec->low+=ec->range*cdf>>AC3_PROB_BITS;
	ec->range=(ec->range*freq>>AC3_PROB_BITS)-1;
	acval_dec(0, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
}

AWM_INLINE void ac3_enc_update_N(AC3 *ec, unsigned cdf, unsigned freq, int probbits)//probbits <= 16
{
#ifdef _DEBUG
	if(cdf>=(unsigned)(1<<probbits)||cdf+freq>(unsigned)(1<<probbits))
		LOG_ERROR2("Invalid stats");
#endif
	AC3_RENORM_STATEMENT(!(ec->range>>probbits))
		ac3_enc_renorm(ec);
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if((unsigned)(freq-1)>=(0x10000-1)||cdf+freq<cdf)
		LOG_ERROR2("Invalid stats");
#endif
	ec->low+=ec->range*cdf>>probbits;
	ec->range=(ec->range*freq>>probbits)-1;
	acval_enc(probbits, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
AWM_INLINE unsigned ac3_dec_getcdf_N(AC3 *ec, int probbits)
{
	//unsigned cdf=0, freq=0;
	//
	//if(!(ec->range>>16))
	//{
	//	ec->code=ec->code<<32|*(unsigned*)ec->srcptr;
	//	ec->srcptr+=4;
	//	ec->range=ec->range<<32|0xFFFFFFFF;
	//	ec->low<<=32;
	//
	//	if(ec->range>~ec->low)
	//		ec->range=~ec->low;
	//}
	//unsigned c=(unsigned)(((ec->code-ec->low)<<16|0xFFFF)/ec->range);
	//
	////linear search...
	//
	//ec->low+=ec->range*cdf>>16;
	//ec->range=(ec->range*freq>>16)-1;

#ifdef AC3_DEC_BRANCHLESSRENORM
	const unsigned char *nextptr=ec->srcptr+4;
	unsigned long long r2=ec->range<<AC3_RENORM|((1LL<<AC3_RENORM)-1);
	unsigned long long low2=ec->low<<AC3_RENORM;
	unsigned long long code2=ec->code<<AC3_RENORM|*(unsigned*)ec->srcptr;
	unsigned long long norenormflag=(unsigned long long)(ec->range>>probbits);
	if(!norenormflag)		//FIXME ensure CMOVs
		ec->range=r2;
	if(!norenormflag)
		ec->low=low2;
	if(!norenormflag)
		ec->code=code2;
	if(!norenormflag)
		ec->srcptr=nextptr;

	unsigned long long rmax=~ec->low;
	if(ec->range>rmax)
		ec->range=rmax;
#else
	AC3_RENORM_STATEMENT(!(ec->range>>probbits))
	{
		unsigned long long rmax;

		ec->code=ec->code<<AC3_RENORM|*(unsigned*)ec->srcptr;
		ec->srcptr+=AC3_RENORM/8;
		ec->range=ec->range<<AC3_RENORM|((1LL<<AC3_RENORM)-1);
		ec->low<<=AC3_RENORM;

		rmax=~ec->low;
		if(ec->range>rmax)
			ec->range=rmax;
	}
#endif
#if 0
	unsigned long long num=((ec->code-ec->low+1)<<probbits)-1, den=ec->range;
	int nzeros=(int)_lzcnt_u64(den);
	unsigned long long
		x1=-(long long)(den<<nzeros),
		x2=__umulh(x1, x1),
		x3=__umulh(x2, x2),
		x4=__umulh(x3, x3);
	//	x5=__umulh(x4, x4),
	//	x6=__umulh(x5, x5);
	x1>>=1;
	x2>>=1;
	x3>>=1;
	x4>>=1;
	//x5>>=1;
	//x6>>=1;
	x1|=1ULL<<63;
	x2|=1ULL<<63;
	x3|=1ULL<<63;
	x4|=1ULL<<63;
	//x5|=1ULL<<63;
	//x6|=1ULL<<63;
	//unsigned long long debug[]={x1, x2, x3, x4, x5, x6};
	x1=__umulh(x1, x2);
	x3=__umulh(x3, x4);
	x1=__umulh(x1, x3);
	//x1=__umulh(x1, x5);
	//x1=__umulh(x1, x6);
	x1=__umulh(num, x1);
	//x1+=1ULL<<(64-4-nzeros)>>1;
	x1>>=64-4-nzeros;
	//if(num/den!=x1)
	//	LOG_ERROR("");
	return x1&((1ULL<<probbits)-1);
#endif
	return (unsigned)((((ec->code-ec->low+1)<<probbits)-1)/ec->range);

	//return (unsigned)(((ec->code-ec->low)<<probbits|((1LL<<probbits)-1))/ec->range);
}
AWM_INLINE void ac3_dec_update_N(AC3 *ec, unsigned cdf, unsigned freq, int probbits)
{
#ifdef _DEBUG
	if(cdf>=(unsigned)(1<<probbits)||cdf+freq>(unsigned)(1<<probbits))
		LOG_ERROR2("Invalid stats");
#endif
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if((unsigned)(freq-1)>=(unsigned)((1<<probbits)-1)||cdf+freq<cdf)
		LOG_ERROR2("Invalid stats");
#endif
	ec->low+=ec->range*cdf>>probbits;
	ec->range=(ec->range*freq>>probbits)-1;//must decrement hi because decoder fails when code == hi2
	acval_dec(probbits, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
}

AWM_INLINE void ac3_enc_update_NPOT(AC3 *ec, unsigned cdf, unsigned freq, unsigned den)
{
#ifdef AC_VALIDATE
	unsigned long long lo0, r0;
	if((unsigned)(freq-1)>=(den-1))
		LOG_ERROR2("ZPS");
	if(cdf+freq<cdf||cdf+freq>den)
		LOG_ERROR2("Invalid CDF");
#endif
	AC3_RENORM_STATEMENT(ec->range<(unsigned)den)
		ac3_enc_renorm(ec);
#ifdef AC_VALIDATE
	lo0=ec->low, r0=ec->range;
#endif
	//unsigned	//X
	//	ncdf=(cdf<<16)/den,
	//	nfreq=(freq<<16)/den;
	//ec->low+=ec->range*ncdf>>16;
	//ec->range=(ec->range*nfreq>>16)-1;

	//__m128d t0=_mm_set1_pd((double)ec->range/den);	//X
	//__m128d t1=_mm_cvtepi32_pd(_mm_set_epi32(0, 0, freq, cdf));
	//t0=_mm_mul_pd(t0, t1);
	//unsigned long long r0=_mm_cvtsd_si64(t0);
	//unsigned long long r1=_mm_cvtsd_si64(_mm_shuffle_pd(t0, t0, 1));
	//ec->low+=r0;
	//ec->range=r1;

	ec->low+=ec->range*cdf/den;
	ec->range=ec->range*freq/den-1;
	acval_enc(den, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
AWM_INLINE unsigned ac3_dec_getcdf_NPOT(AC3 *ec, unsigned den)
{
	AC3_RENORM_STATEMENT(ec->range<(unsigned)den)
		ac3_dec_renorm(ec);
#ifdef _DEBUG
	if(ec->code<ec->low||ec->code>ec->low+ec->range)
		LOG_ERROR("CODE OOB");
#endif
#ifdef AC_VALIDATE
	ACVAL *val=(ACVAL*)array_at(&acval, acval_idx);
	if(ec->code<val->lo2||ec->code>val->hi2)
		LOG_ERROR("CODE OOB");
#endif
	return (unsigned)(((ec->code-ec->low)*den+den-1)/ec->range);
}
AWM_INLINE void ac3_dec_update_NPOT(AC3 *ec, unsigned cdf, unsigned freq, unsigned den)
{
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if((unsigned)(freq-1)>=(den-1))
		LOG_ERROR2("ZPS");
	if(cdf+freq<cdf)
		LOG_ERROR2("Invalid CDF");
#endif
	//unsigned	//X
	//	ncdf=(cdf<<16)/den,
	//	nfreq=(freq<<16)/den;
	//ec->low+=ec->range*ncdf>>16;
	//ec->range=(ec->range*nfreq>>16)-1;

	ec->low+=ec->range*cdf/den;
	ec->range=ec->range*freq/den-1;
	acval_dec(den, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
}

AWM_INLINE void ac3_enc(AC3 *ec, int sym, const unsigned *CDF)
{
	unsigned cdf, freq;

	AC3_RENORM_STATEMENT(!(ec->range>>AC3_PROB_BITS))
		ac3_enc_renorm(ec);
	cdf=CDF[sym];
	freq=CDF[sym+1]-cdf;
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!freq)
		LOG_ERROR2("ZPS");
	if(cdf+freq<cdf)
		LOG_ERROR2("Invalid CDF");
#endif
	ec->low+=ec->range*cdf>>AC3_PROB_BITS;
	ec->range=(ec->range*freq>>AC3_PROB_BITS)-1;
	acval_enc(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
AWM_INLINE int ac3_dec(AC3 *ec, const unsigned *CDF, int nbits)
{
	unsigned cdf, freq;
	int sym;

	AC3_RENORM_STATEMENT(!(ec->range>>AC3_PROB_BITS))
		ac3_dec_renorm(ec);
	cdf=(unsigned)(((ec->code-ec->low)<<AC3_PROB_BITS|((1ULL<<AC3_PROB_BITS)-1))/ec->range);
	sym=0;
	switch(nbits)
	{
	default:
		LOG_ERROR2("Unsupported bit depth");
		break;
	case 9:sym|=(cdf>=CDF[sym|256])<<8;
	case 8:sym|=(cdf>=CDF[sym|128])<<7;
	case 7:sym|=(cdf>=CDF[sym| 64])<<6;
	case 6:sym|=(cdf>=CDF[sym| 32])<<5;
	case 5:sym|=(cdf>=CDF[sym| 16])<<4;
	case 4:sym|=(cdf>=CDF[sym|  8])<<3;
	case 3:sym|=(cdf>=CDF[sym|  4])<<2;
	case 2:sym|=(cdf>=CDF[sym|  2])<<1;
	case 1:sym|= cdf>=CDF[sym|  1];
		break;
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
	ec->low+=ec->range*cdf>>AC3_PROB_BITS;
	ec->range=(ec->range*freq>>AC3_PROB_BITS)-1;
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return sym;
}
AWM_INLINE int ac3_dec_NPOT(AC3 *ec, const unsigned *CDF, int nlevels)
{
	unsigned cdf, freq;
	int range, sym;

	AC3_RENORM_STATEMENT(!(ec->range>>AC3_PROB_BITS))
		ac3_dec_renorm(ec);
	cdf=(unsigned)((((ec->code-ec->low+1)<<AC3_PROB_BITS)-1)/ec->range);
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
	ec->low+=ec->range*cdf>>AC3_PROB_BITS;
	ec->range=(ec->range*freq>>AC3_PROB_BITS)-1;
	acval_dec(sym, cdf, freq, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return sym;
}

AWM_INLINE void ac3_enc_bypass(AC3 *ec, int bypass, int nbits)//up to 16 bits
{
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if((unsigned)bypass>=(unsigned)(1<<nbits))
		LOG_ERROR("Bypass OOB");
#endif
	AC3_RENORM_STATEMENT(!(ec->range>>nbits))
		ac3_enc_renorm(ec);
	ec->low+=ec->range*bypass>>nbits;
	ec->range=(ec->range>>nbits)-1;
	acval_enc(nbits, bypass, 1, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
AWM_INLINE int ac3_dec_bypass(AC3 *ec, int nbits)
{
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
#endif
	AC3_RENORM_STATEMENT(!(ec->range>>nbits))
		ac3_dec_renorm(ec);
#ifdef _DEBUG
	if(ec->code<ec->low||ec->code>ec->low+ec->range)
		LOG_ERROR("CODE OOB");
#endif
	int bypass=(int)((((ec->code-ec->low+1)<<nbits)-1)/ec->range);
	ec->low+=ec->range*bypass>>nbits;
	ec->range=(ec->range>>nbits)-1;
	acval_dec(nbits, bypass, 1, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return bypass;
}

AWM_INLINE void ac3_enc_bypass_NPOT(AC3 *ec, int bypass, int nlevels)
{
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if((unsigned)bypass>=(unsigned)nlevels)
		LOG_ERROR("Bypass OOB");
#endif
	AC3_RENORM_STATEMENT(ec->range<(unsigned)nlevels)
		ac3_enc_renorm(ec);
	ec->low+=ec->range*bypass/nlevels;
	ec->range=ec->range/nlevels-1;
	acval_enc(nlevels, bypass, 1, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
AWM_INLINE int ac3_dec_bypass_NPOT(AC3 *ec, int nlevels)
{
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
#endif
	
	AC3_RENORM_STATEMENT(ec->range<(unsigned)nlevels)
		ac3_dec_renorm(ec);
	int bypass=(int)(((ec->code-ec->low)*nlevels+nlevels-1)/ec->range);
	ec->low+=ec->range*bypass/nlevels;
	ec->range=ec->range/nlevels-1;
	acval_dec(nlevels, bypass, 1, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
	return bypass;
}

AWM_INLINE void ac3_enc_bin(AC3 *ec, int bit, unsigned p0, int probbits)//probbits <= 16
{
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!p0||p0>=1ULL<<probbits)
		LOG_ERROR2("ZPS");
#endif
	AC3_RENORM_STATEMENT(!(ec->range>>probbits))
		ac3_enc_renorm(ec);
	unsigned long long mid=ec->range*p0>>probbits;
	if(bit)
	{
		ec->low+=mid;
		ec->range-=mid;
	}
	else
		ec->range=mid-1;
	acval_enc(bit, bit?p0:0, bit?(1<<probbits)-p0:p0, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0);
}
AWM_INLINE int ac3_dec_bin(AC3 *ec, unsigned p0, int probbits)
{
#ifdef AC_VALIDATE
	unsigned long long lo0=ec->low, r0=ec->range;
	if(!p0||p0>=1ULL<<probbits)
		LOG_ERROR2("ZPS");
#endif
	AC3_RENORM_STATEMENT(!(ec->range>>probbits))
		ac3_dec_renorm(ec);
	unsigned long long mid=ec->range*p0>>probbits;
	unsigned long long t2=ec->low+mid;
	ec->range-=mid;
	int bit=ec->code>=t2;
	if(bit)
		ec->low=t2;
	else
		ec->range=mid-1;
	acval_dec(bit, bit?p0:0, bit?(1<<probbits)-p0:p0, lo0, lo0+r0, ec->low, ec->low+ec->range, 0, 0, ec->code);
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
AWM_INLINE void ac4_enc_init(AC4 *ec, BList *dst)
{
	memset(ec, 0, sizeof(*ec));
	ec->lo=0;
	ec->hi=0xFFFFFFFF;
	ec->dst=dst;
}
AWM_INLINE void ac4_dec_init(AC4 *ec, const unsigned char *srcstart, const unsigned char *srcend)
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
AWM_INLINE void ac4_enc_flush(AC4 *ec)
{
	unsigned mid=(unsigned)(((unsigned long long)ec->lo+ec->hi)>>1);
	blist_push_back1(ec->dst, (unsigned char*)&mid+3);//big-endian
	blist_push_back1(ec->dst, (unsigned char*)&mid+2);
	blist_push_back1(ec->dst, (unsigned char*)&mid+1);
	blist_push_back1(ec->dst, (unsigned char*)&mid+0);
}
AWM_INLINE void ac4_enc_bin(AC4 *ec, unsigned short p1, int bit)
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
AWM_INLINE int ac4_dec_bin(AC4 *ec, unsigned short p1)
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
AWM_INLINE void ac4_enc_update_NPOT(AC4 *ec, unsigned cdf_curr, unsigned cdf_next, unsigned den)
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
AWM_INLINE unsigned ac4_dec_getcdf_NPOT(AC4 *ec, unsigned den)
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
AWM_INLINE void ac4_dec_update_NPOT(AC4 *ec, unsigned cdf_curr, unsigned cdf_next, unsigned den)
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
AWM_INLINE void ac5_enc_init(AC5 *ec, FILE *fdst)
{
	ec->lo=0;
	ec->hi=0xFFFFFFFF;
	ec->code=0;
	ec->enc=1;
	ec->f=fdst;
}
AWM_INLINE void ac5_dec_init(AC5 *ec, FILE *fsrc)
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
AWM_INLINE void ac5_enc_flush(AC5 *ec)
{
	unsigned code=ec->lo;
	fputc(code>>24, ec->f);	code<<=8;
	fputc(code>>24, ec->f);	code<<=8;
	fputc(code>>24, ec->f);	code<<=8;
	fputc(code>>24, ec->f);
}
AWM_INLINE void ac5_enc_bin(AC5 *ec, int p1, int bit)
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
AWM_INLINE int ac5_dec_bin(AC5 *ec, int p1)
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


//Bypass Coder (POT)		rbuf with branchless renorm requires read/write access 4 bytes past the end
typedef struct _BypassCoder
{
	unsigned long long state;
	int enc_nfree, dec_navailable;
	BList *dst;
	unsigned char *dstptr, *dststart, *dstend;
	const unsigned char *srcptr, *srcstart, *srcend;
} BypassCoder;
AWM_INLINE void bypass_encbuf_init(BypassCoder *ec, unsigned char *start, unsigned char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->enc_nfree=64;
	ec->dstptr=start;
	ec->dststart=start;
	ec->dstend=end;
}
AWM_INLINE void bypass_decbuf_init(BypassCoder *ec, const unsigned char *start, const unsigned char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->srcptr=start;
	ec->srcstart=start;
	ec->srcend=end;
}
AWM_INLINE unsigned char* bypass_encbuf_flush(BypassCoder *ec)
{
#ifndef AC_DISABLE_WRITE
	*(unsigned*)ec->dstptr=(unsigned)(ec->state>>32);
	ec->dstptr+=4;
	*(unsigned*)ec->dstptr=(unsigned)ec->state;
	ec->dstptr+=4;
#endif
	return ec->dstptr;
}
AWM_INLINE void bypass_encbuf(BypassCoder *ec, int bypass, int nbits)
{
	//int nfree=ec->enc_nfree-nbits;
#ifdef BYPASS_ENC_BRANCHLESSRENORM
	unsigned char *dstptr=ec->dstptr+4;
	unsigned emit=(unsigned)(ec->state>>32);
	int renorm=ec->enc_nfree<32;
	if(renorm)				//FIXME can be replaced with direct assignment
		*(unsigned*)ec->dstptr=emit;
	if(renorm)				//FIXME ensure CMOV
		ec->dstptr=dstptr;
	ec->state<<=32*renorm;
	ec->enc_nfree+=32*renorm;
#else
	if(ec->enc_nfree<32)
	{
#ifdef _DEBUG
		if(ec->dstptr>=ec->dstend)
		{
			LOG_ERROR("Encoder bypass buffer overflow");
			return;
		}
#endif
#ifndef AC_DISABLE_WRITE
		*(unsigned*)ec->dstptr=(unsigned)(ec->state>>32);
		ec->dstptr+=4;
#endif
		ec->state<<=32;
		//nfree+=32;
		ec->enc_nfree+=32;
	}
#endif
	ec->enc_nfree-=nbits;
	ec->state|=(unsigned long long)bypass<<ec->enc_nfree;
	//ec->enc_nfree=nfree;
}
AWM_INLINE int bypass_decbuf(BypassCoder *ec, int nbits)
{
#ifdef BYPASS_DEC_BRANCHLESSRENORM
	const unsigned char *srcptr=ec->srcptr+4;
	unsigned long long state2=ec->state<<32|*(const unsigned*)ec->srcptr;
	int renorm=ec->dec_navailable<32;
	if(renorm)				//FIXME ensure CMOV
		ec->state=state2;
	if(renorm)
		ec->srcptr=srcptr;
	ec->dec_navailable+=32*renorm-nbits;

	int bypass=_bextr_u64(ec->state, ec->dec_navailable, nbits);//BMI1 (2013)
//	int bypass=ec->state>>navailable&((1LL<<nbits)-1);
#else
	//int navailable=ec->dec_navailable-nbits;
	if(ec->dec_navailable<32)
	{
#ifdef _DEBUG
		if(ec->srcptr>=ec->srcend)
		{
			LOG_ERROR("Decoder bypass buffer overflow");
			return 0;
		}
#endif
		ec->state=ec->state<<32|*(unsigned*)ec->srcptr;
		ec->srcptr+=4;
		//navailable+=32;
		ec->dec_navailable+=32;
	}
	ec->dec_navailable-=nbits;

	int bypass=(int)_bextr_u64(ec->state, ec->dec_navailable, nbits);//BMI1 (2013)
//	int bypass=ec->state>>ec->dec_navailable&((1LL<<nbits)-1);
	//ec->dec_navailable=navailable;
#endif
	return bypass;
}

AWM_INLINE void bypass_enc_init(BypassCoder *ec, BList *dst)
{
	memset(ec, 0, sizeof(*ec));
	ec->enc_nfree=64;
	ec->dst=dst;
}
AWM_INLINE void bypass_dec_init(BypassCoder *ec, const unsigned char *start, const unsigned char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->srcptr=start;
	ec->srcend=end;
}
AWM_INLINE void bypass_enc_flush(BypassCoder *ec)
{
	blist_push_back(ec->dst, (unsigned char*)&ec->state+(sizeof(ec->state)>>1)*1, sizeof(ec->state)>>1);//big-endian
	blist_push_back(ec->dst, (unsigned char*)&ec->state+(sizeof(ec->state)>>1)*0, sizeof(ec->state)>>1);
}
AWM_INLINE void bypass_enc(BypassCoder *ec, int bypass, int nbits)
{
	if(ec->enc_nfree<(sizeof(ec->state)<<3>>1))	//32-bit renorm
	{
		blist_push_back(ec->dst, (unsigned char*)&ec->state+(sizeof(ec->state)>>1), sizeof(ec->state)>>1);
		ec->state<<=sizeof(ec->state)<<3>>1;
		ec->enc_nfree+=sizeof(ec->state)<<3>>1;
	}
	ec->enc_nfree-=nbits;
	ec->state|=(unsigned long long)bypass<<ec->enc_nfree;

	//ec->nbits+=nbits;		//64-bit renorm
	//if(ec->nbits>(sizeof(ec->state)<<3))
	//{
	//	int hicount=ec->nbits&((sizeof(ec->state)<<3)-1);
	//	ec->state<<=hicount;
	//	ec->state|=(unsigned long long)bypass>>(nbits-hicount);
	//	blist_push_back(ec->dst, &ec->state, sizeof(ec->state)<<3);
	//	ec->state=bypass&((1LL<<(nbits-hicount))-1);
	//	return;
	//}
	//ec->state=ec->state<<nbits|bypass;

	//ec->state=ec->state<<nbits|bypass;				//convention #1: shift state
	//ec->state|=(unsigned long long)bypass<<ec->nbits;		//convention #2: shift bypass
	//ec->nbits+=nbits;
}
AWM_INLINE int bypass_dec(BypassCoder *ec, int nbits)
{
	if(ec->dec_navailable<(sizeof(ec->state)<<3>>1))		//32-bit renorm
	{
		ec->state=ec->state<<(sizeof(ec->state)<<3>>1)|*(unsigned*)ec->srcptr;
		ec->srcptr+=sizeof(ec->state)>>1;
		ec->dec_navailable+=sizeof(ec->state)<<3>>1;
	}
	ec->dec_navailable-=nbits;
	int bypass=ec->state>>ec->dec_navailable&((1LL<<nbits)-1);
	return bypass;

	//int bypass=ec->state&((1LL<<nbits)-1);		//convention #1: shift state
	//ec->state>>=nbits;
	//return bypass;
	
	//ec->nbits-=nbits;
	//int bypass=ec->state>>ec->nbits&((1LL<<nbits)-1);	//convention #2: shift bypass
	//return bypass;
}


//LIFO Bypass Coder
typedef struct _BitPackerLIFO//bwd enc / fwd dec
{
	unsigned long long state;
	int enc_nwritten, dec_navailable;//bitcounts, only for tracking renorms
	unsigned char *dstbwdptr;
	const unsigned char *srcfwdptr, *streamend;
} BitPackerLIFO;
AWM_INLINE void bitpacker_enc_init(BitPackerLIFO *ec, const unsigned char *bufstart, unsigned char *bufptr0_OOB)
{
	memset(ec, 0, sizeof(*ec));
	ec->state=1ULL<<32;
	ec->enc_nwritten=33;
	ec->streamend=bufstart;
	ec->dstbwdptr=bufptr0_OOB;
}
AWM_INLINE void bitpacker_dec_init(BitPackerLIFO *ec, const unsigned char *bufptr0_start, const unsigned char *bufend)
{
	memset(ec, 0, sizeof(*ec));
	ec->srcfwdptr=bufptr0_start+8;
	ec->streamend=bufend;
	ec->state=*(const unsigned long long*)bufptr0_start;
	ec->dec_navailable=FLOOR_LOG2_P1(ec->state);
}
AWM_INLINE void bitpacker_enc_flush(BitPackerLIFO *ec)
{
	ec->dstbwdptr-=8;
#ifdef _DEBUG
	if(ec->dstbwdptr<ec->streamend)
		LOG_ERROR("IntPacker Encoder OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
	*(unsigned long long*)ec->dstbwdptr=ec->state;
}
AWM_INLINE void bitpacker_enc(BitPackerLIFO *ec, int inbits, int sym)
{
#ifdef _DEBUG
	if(!inbits)
		LOG_ERROR("BitPacker inbits=0");
#endif
	//renorm then push inbits
	ec->enc_nwritten+=inbits;
	if(ec->enc_nwritten>64)//renorm on overflow
	{
		ec->enc_nwritten-=32;
		ec->dstbwdptr-=4;
#ifdef _DEBUG
		if(ec->dstbwdptr<ec->streamend)
			LOG_ERROR("IntPacker OOB:  dstbwdptr = 0x%016zX < 0x%016zX", ec->dstbwdptr, ec->streamend);
#endif
		*(unsigned*)ec->dstbwdptr=(unsigned)ec->state;
		ec->state>>=32;
#ifdef ANS_VAL
		ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
	}
	ec->state=ec->state<<inbits|sym;
#ifdef ANS_VAL
	ansval_push(&ec->state, sizeof(ec->state), 1);
#endif
}
AWM_INLINE int bitpacker_dec(BitPackerLIFO *ec, int outbits)
{
#ifdef _DEBUG
	if(!outbits)
		LOG_ERROR("BitPacker outbits=0");
#endif
	int sym=ec->state&((1ULL<<outbits)-1);

	//pop outbits then renorm
#ifdef ANS_VAL
	ansval_check(&ec->state, sizeof(ec->state), 1);
#endif
	ec->dec_navailable-=outbits;
	ec->state>>=outbits;
	if(ec->dec_navailable<=32)
	{
#ifdef ANS_VAL
		ansval_check(&ec->state, sizeof(ec->state), 1);
#endif
		ec->dec_navailable+=32;
#ifdef _DEBUG
		if(ec->srcfwdptr>=ec->streamend)
			LOG_ERROR("IntPacker OOB:  srcfwdptr = 0x%016zX >= 0x%016zX", ec->srcfwdptr, ec->streamend);
#endif
		ec->state=ec->state<<32|*(const unsigned*)ec->srcfwdptr;
		ec->srcfwdptr+=4;
	}
	return sym;
}


//Golomb-Rice Coder

//	typedef unsigned GREmit_t;
	typedef unsigned long long GREmit_t;//0.001% larger, 18% faster enc

typedef struct GolombRiceCoderStruct
{
	GREmit_t cache;
	int nbits;//enc: number of free bits in cache, dec: number of unread bits in cache
	int is_enc;//for padding
	const unsigned char *srcptr, *srcend, *srcstart;
	unsigned char *dstptr, *dstend, *dststart;
	BList *dst;
} GolombRiceCoder;
AWM_INLINE void gr_enc_flush(GolombRiceCoder *ec)
{
#ifdef GR_USE_ARRAY
#ifdef _DEBUG
	if(ec->dstptr+sizeof(ec->cache)>ec->dstend)//compression failed
		LOG_ERROR("GR buffer overflow");
#endif
	*(GREmit_t*)ec->dstptr=ec->cache;
	//memcpy(ec->dstptr, &ec->cache, sizeof(ec->cache));
	ec->dstptr+=sizeof(ec->cache);
#else
	blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));//size is qword-aligned
#endif
}
AWM_INLINE int gr_dec_impl_read(GolombRiceCoder *ec)
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
AWM_INLINE void gr_enc_init(GolombRiceCoder *ec,
#ifdef GR_USE_ARRAY
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
#ifdef GR_USE_ARRAY
	ec->dstptr=start;
	ec->dstend=end;
	ec->dststart=start;
#else
	ec->dst=dst;
#endif
}
AWM_INLINE void gr_dec_init(GolombRiceCoder *ec, const unsigned char *start, const unsigned char *end)
{
	memset(ec, 0, sizeof(*ec));
	ec->cache=0;
	ec->nbits=0;
	ec->is_enc=0;
	ec->srcptr=start;
	ec->srcend=end;
	ec->srcstart=start;
}

AWM_INLINE int gr_enc_NPOT(GolombRiceCoder *ec, unsigned sym, unsigned magnitude)
{
	//buffer: {c,c,c,b,b,a,a,a, f,f,f,e,e,e,d,c}, cache: MSB gg[hhh]000 LSB	nbits 6->3, code h is about to be emitted
	//written 64-bit words are byte-reversed because the CPU is little-endian

	//magnitude+=!magnitude;
	int nbypass=32-_lzcnt_u32(magnitude);//FLOOR_LOG2(magnitude)+1
	int nzeros=sym/magnitude, bypass=sym%magnitude;
	if(nzeros>=ec->nbits)//fill the rest of cache with zeros, and flush
	{
		nzeros-=ec->nbits;
		gr_enc_flush(ec);
		//if(!gr_enc_flush(ec))
		//	return 0;
		//blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		if(nzeros>=(int)(sizeof(ec->cache)<<3))//just flush zeros
		{
			ec->cache=0;
			do
			{
				nzeros-=(sizeof(ec->cache)<<3);
				gr_enc_flush(ec);
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
		gr_enc_flush(ec);
		//blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		ec->cache=0;
		ec->nbits=sizeof(ec->cache)<<3;
	}
	//now there is room for bypass:  0 <= nbypass < nbits <= 64
	ec->nbits-=nbypass;//emit remaining bypass to cache
	ec->cache|=(GREmit_t)bypass<<ec->nbits;
	return 1;
}
AWM_INLINE unsigned gr_dec_NPOT(GolombRiceCoder *ec, unsigned magnitude)
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
	unsigned bypass=0, nunused=(1<<(nbypass+1))-magnitude;//this reads k=FLOOR_LOG2(magnitude) bits, then possibly reads one more bit
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
	if(bypass>=nunused)//read one more bit
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

AWM_INLINE int gr_enc(GolombRiceCoder *ec, int sym, int nbypass)
{
	//buffer: {c,c,c,b,b,a,a,a, f,f,f,e,e,e,d,c}, cache: MSB gg[hhh]000 LSB	nbits 6->3, code h is about to be emitted
	//written 64-bit words are byte-reversed because the CPU is little-endian

	int nzeros=sym>>nbypass, bypass=sym&((1<<nbypass)-1);
	if(nzeros>=ec->nbits)//fill the rest of cache with zeros, and flush
	{
		nzeros-=ec->nbits;
		gr_enc_flush(ec);
		//blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		ec->cache=0;
		if(nzeros>=(int)(sizeof(ec->cache)<<3))//just flush zeros
		{
			do
			{
				nzeros-=(sizeof(ec->cache)<<3);
				gr_enc_flush(ec);
				//blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
			}
			while(nzeros>(int)(sizeof(ec->cache)<<3));
		}
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
		gr_enc_flush(ec);
		//blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		ec->cache=0;
		ec->nbits=sizeof(ec->cache)<<3;
	}
	//now there is room for bypass:  0 <= nbypass < nbits <= 64
	ec->nbits-=nbypass;//emit remaining bypass to cache
	ec->cache|=(GREmit_t)bypass<<ec->nbits;
	return 1;
}
AWM_INLINE unsigned gr_dec(GolombRiceCoder *ec, int nbypass)
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

AWM_INLINE int gr_enc_bypass(GolombRiceCoder *ec, int sym, int nbits)
{
	//buffer: {c,c,c,b,b,a,a,a, f,f,f,e,e,e,d,c}, cache: MSB gg[hhh]000 LSB	nbits 6->3, code h is about to be emitted
	//written 64-bit words are byte-reversed because the CPU is little-endian

	if(nbits>=ec->nbits)//not enough free bits in cache:  fill cache, write to list, and repeat
	{
		nbits-=ec->nbits;
		ec->cache|=(unsigned long long)sym>>nbits;
		sym&=(1<<nbits)-1;
		gr_enc_flush(ec);
		//blist_push_back(ec->dst, &ec->cache, sizeof(ec->cache));
		ec->cache=0;
		ec->nbits=sizeof(ec->cache)<<3;
	}
	//now there is room for bypass:  0 <= nbypass < nbits <= 64
	ec->nbits-=nbits;//emit remaining bypass to cache
	ec->cache|=(GREmit_t)sym<<ec->nbits;
	return 1;
}
AWM_INLINE unsigned gr_dec_bypass(GolombRiceCoder *ec, int nbits)
{
	//cache: MSB 00[hhh]ijj LSB		nbits 6->3, h is about to be read (past codes must be cleared from cache)

	unsigned sym=0;
	if(ec->nbits<nbits)
	{
		nbits-=ec->nbits;
		sym|=(int)(ec->cache<<nbits);
		if(gr_dec_impl_read(ec))
			return 0;
	}
	if(nbits)		//FIXME minimize branches & quick return for common case
	{
		ec->nbits-=nbits;
		sym|=(int)(ec->cache>>ec->nbits);
		ec->cache&=(1ULL<<ec->nbits)-1;
	}
	return sym;
}


#ifdef __cplusplus
}
#endif
#endif
