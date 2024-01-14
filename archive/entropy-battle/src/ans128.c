#include"battle.h"
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#include<string.h>
#include<stdio.h>//for debugging
static const char file[]=__FILE__;

//128-bit ANS demonstration
#if 1
#define ADD128(LARGE, SMALL)	*(LARGE)+=(SMALL), (LARGE)[1]+=*(LARGE)<(SMALL)
void ans128_enc_update(unsigned long long *state, unsigned long long CDF_val, unsigned long long freq, DList *list)
{
	if(!freq)
		LOG_ERROR("ZPS");

	if(state[1]>freq||state[1]==freq&&state[0])
	{
		dlist_push_back(list, state, 8);
		state[0]=state[1];
		state[1]=0;
	}

	state[1]=_udiv128(state[1], state[0], freq, state);
	state[0]+=CDF_val;

	//for(int k=127;k>=0;--k)
	//{
	//	int bit=(&state->lo)[k>=64]>>(k&63)&1;
	//	r<<=1;
	//	r|=bit;
	//	r-=freq&-(r>freq);
	//	(&q.lo)[k>=64]|=(unsigned long long)bit<<(k&63);
	//}
}
void ans128_dec_update(unsigned long long *state, const unsigned char *srcstart, const unsigned char **srcptr, unsigned long long CDF_val, unsigned long long freq)
{
	unsigned long long delta=state[0]-CDF_val;
	state[0]=_mul128(state[0], freq, state+1);
	ADD128(state, delta);

	if(!state[1])
	{
		state[1]=state[0];
		if(*srcptr-8<srcstart)
		{
			*srcptr-=8;
			memcpy(state, *srcptr, 8);
		}
	}
}
#endif

typedef struct SymbolStruct
{
	unsigned freq, sym;
} Symbol;
int cmp_uint32(const void *p1, const void *p2)
{
	unsigned u1=*(const unsigned*)p1, u2=*(const unsigned*)p2;
	return (u1>u2)-(u1<u2);
}
int cmp_uint64(const void *p1, const void *p2)
{
	unsigned long long u1=*(const unsigned long long*)p1, u2=*(const unsigned long long*)p2;
	return (u1>u2)-(u1<u2);
}
size_t ans128_encode(const void *src, ptrdiff_t symcount, int *depths, int nch, int bytestride, ArrayHandle *data)//depth up to 32 bit per channel, channels can be packed, each symbol must start at new byte
{
	if(depths)
	{
		for(int kc=0;kc<nch;++kc)
		{
			if(depths[kc]>32)
			{
				LOG_ERROR("Channel depth > 32 is not supported");
				return 0;
			}
		}
	}
	const unsigned char *buf=(const unsigned char*)src;
	unsigned
		*ch=(unsigned*)malloc(symcount*sizeof(unsigned)),
		*ans_offsets=(unsigned*)malloc(nch*sizeof(unsigned));
	Symbol *pal=(Symbol*)malloc(symcount*sizeof(Symbol));
	if(!ch||!pal||!ans_offsets)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	size_t overheadsize=0;
	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, 4LL*nch);
	memset(pal, 0, symcount*sizeof(Symbol));
	for(int kc=0, bitoffset=0;kc<nch;++kc)
	{
		int depth=depths?depths[kc]:8;
		int bytestart=bitoffset>>3, chbytes=(depth+7)>>3;
		memset(ch, 0, symcount*sizeof(unsigned));
		for(ptrdiff_t ks=0;ks<symcount;++ks)//copy channel
		{
			memcpy(ch+ks, buf+bytestride*ks+bytestart, chbytes);
			ch[ks]>>=bitoffset&7;
		}

		//memset(hist, 0, symcount*sizeof(unsigned));
		ptrdiff_t nlevels=0;
		for(ptrdiff_t ks=0;ks<symcount;++ks)//build palette, substitute symbols & calculate histogram
		{
			unsigned val=ch[ks];
			size_t idx=0;
			int found=binary_search(pal, nlevels, 4, cmp_uint32, &val, &idx);
			if(!found)
			{
				memmove(pal+idx+1, pal+idx, (nlevels-idx)*sizeof(Symbol));
				++nlevels;
				pal[idx].sym=val;
				pal[idx].freq=0;
			}
			++pal[idx].freq;
			ch[ks]=(unsigned)idx;//replace value with palette idx
		}
		
		dlist_push_back1(&list, "P");
		dlist_push_back(&list, &nlevels, 4);
		for(int idx=0;idx<nlevels;++idx)
			dlist_push_back(&list, pal+idx, 4LL+chbytes);
		overheadsize+=nlevels*(4LL+chbytes)+1;
#if 0
		unsigned char tag='P'+chbytes;
		dlist_push_back1(&list, &tag);
		//dlist_push_back1(&list, "P");
		//dlist_push_back1(&list, &chbytes);
		dlist_push_back(&list, &nlevels, 4);
		for(int idx=0;idx<nlevels;++idx)//bypass palette
			dlist_push_back(&list, pal+idx, chbytes);
		overheadsize+=nlevels*chbytes+5;

		dlist_push_back1(&list, "H");
		dlist_push_back(&list, hist, nlevels*sizeof(unsigned));//bypass histogram
		overheadsize+=nlevels*sizeof(unsigned)+1;
#endif

		unsigned long long *CDF=(unsigned long long*)malloc(nlevels*sizeof(unsigned long long));
		unsigned long long sum=0;
		for(int idx=0;idx<nlevels;++idx)//normalize histogram
		{
			unsigned long long qfreq, r=0;
			qfreq=_div128(pal[idx].freq, 0, symcount, &r);
			if(!qfreq)
			{
				LOG_ERROR("ZPS");
				return 0;
			}
			CDF[idx]=sum;
			sum+=qfreq;
		}

		//dlist_push_back1(&list, "A");
		unsigned long long state[2]={1, 0};
		for(ptrdiff_t ks=symcount-1;ks>=0;--ks)//encode
		{
			unsigned idx=ch[ks];
			unsigned long long cdf=CDF[idx], freq=(idx<nlevels-1?CDF[idx+1]:0)-cdf;
			ans128_enc_update(state, cdf, freq, &list);
		}
		free(CDF);
		dlist_push_back(&list, state, 16);
		ans_offsets[kc]=(unsigned)list.nobj;
		if(kc>0)
			ans_offsets[kc]-=ans_offsets[kc-1];

		bitoffset+=depth;
	}
	free(ch);
	free(pal);
	size_t puresize=list.nobj-overheadsize;
	
	size_t dststart=0;
	if(*data)
	{
		if(data[0]->esize!=1)
			LOG_ERROR("Invalid destination array");
		dststart=data[0]->count;
	}
	dlist_appendtoarray(&list, data);
	memcpy(data[0]->data+dststart, ans_offsets, (size_t)nch*4);
	//memcpy(data[0]->data+dststart, &list.nobj, 4);

	dlist_clear(&list);
	free(ans_offsets);
	return puresize;
}
int    ans128_decode(const unsigned char *data, ptrdiff_t srclen, ptrdiff_t symcount, int *depths, int nch, int bytestride, void *dst)
{
	if(srclen<4LL*nch)
	{
		LOG_ERROR("Incomplete file");
		return 0;
	}
	unsigned *ans_offsets=(unsigned*)malloc(4LL*nch);
	if(!ans_offsets)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	const unsigned char *srcptr=data, *srcstart;
	memcpy(ans_offsets, srcptr, 4LL*nch);
	srcptr+=4*nch;
	if(srclen<ans_offsets[nch-1])
	{
		LOG_ERROR("Incomplete file");
		return 0;
	}
	unsigned char *buf=(unsigned char*)dst;
	for(int kc=0, bitoffset=0;kc<nch;++kc)
	{
		if(kc>0)
			srcptr=data+ans_offsets[kc-1];
		int depth=depths?depths[kc]:8;
		int bytestart=bitoffset>>3, o2=bitoffset&7, chbytes=(depth+7)>>3;

		if(srcptr+5>data+srclen)
		{
			LOG_ERROR("Incomplete file");
			return 0;
		}
		if(*srcptr!='P')
		{
			LOG_ERROR("Wrong tag");
			return 0;
		}
		int nlevels;
		memcpy(&nlevels, srcptr+1, 4);
		if(nlevels<0||nlevels>0x1000000)
		{
			LOG_ERROR("Wrong nlevels");
			return 0;
		}
		srcptr+=5;

		Symbol *pal=(Symbol*)malloc(nlevels*sizeof(Symbol));
		unsigned long long *CDF=(unsigned long long*)malloc(nlevels*sizeof(unsigned long long));
		if(!CDF||!pal)
		{
			LOG_ERROR("Allocation error");
			return 0;
		}

		for(int idx=0;idx<nlevels;++idx)
		{
			pal[idx].sym=0;
			memcpy(pal+idx, srcptr, (size_t)(4+chbytes));
			srcptr+=(size_t)(4+chbytes);
			//pal[idx]=0;
			//memcpy(pal+idx, srcptr, chbytes);
			//srcptr+=chbytes;
		}

		unsigned long long sum=0;
		for(int idx=0;idx<nlevels;++idx)//normalize histogram
		{
			unsigned long long freq, r=0;
			//memcpy(&freq, srcptr, 4);
			//srcptr+=4;
			freq=_div128(pal[idx].freq, 0, symcount, &r);
			CDF[idx]=sum;
			sum+=freq;
		}
		srcstart=srcptr;
		srcptr=data+ans_offsets[kc];
		unsigned long long state[2];
		srcptr-=16;
		memcpy(state, srcptr, 16);
		for(ptrdiff_t ks=0;ks<symcount;++ks)
		{
			size_t idx;
			binary_search(CDF, nlevels, 8, cmp_uint64, state, &idx);//here CDF has unique values
			idx-=state[0]<CDF[idx];

			unsigned sym=pal[idx].sym;
			unsigned char *sstart=buf+bytestride*ks+bytestart;
			for(int kb=0;kb<depth;++kb)
				sstart[(o2+kb)>>3]|=(sym>>kb&1)<<(o2+kb);

			unsigned long long cdf=CDF[idx], freq=(idx<nlevels-1?CDF[idx+1]:0)-cdf;
			ans128_dec_update(state, srcstart, &srcptr, cdf, freq);
		}
		free(CDF);
		free(pal);
		bitoffset+=depth;
	}
	return 1;
}

size_t ans16_encode(const unsigned char *src, ptrdiff_t res, ArrayHandle *data)
{
	unsigned *hist=(unsigned*)malloc(256*sizeof(unsigned));
	if(!hist)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	DList list;
	dlist_init(&list, 1, 1024, 0);
	dlist_push_back(&list, 0, 4);
	for(int kc=0;kc<3;++kc)
	{
		memset(hist, 0, 256*sizeof(unsigned));
		for(ptrdiff_t ks=0;ks<res;++ks)
		{
			unsigned char sym=src[ks<<2|kc]-64;
			++hist[sym];
		}

		unsigned sum=0;
		for(int sym=0;sym<256;++sym)
		{
			unsigned freq=(unsigned)(((long long)hist[sym]<<16)/res);
			hist[sym]=sum;
			sum+=freq;
		}

		//for(int sym=0;sym<256;++sym)
		//	dlist_push_back(&list, hist+sym, 2);

		unsigned short state=0x0100;
		for(ptrdiff_t ks=res-1;ks>=0;--ks)
		{
			unsigned char sym=src[ks<<2|kc]-64;
			unsigned char cdf[2], freq[2];

			unsigned char start=sym&0xF0, end=start+0x10;
			unsigned c0=hist[start], c1=end<256?hist[end]:0x10000;
			cdf[0]=c0>>8&0xF0|start>>4;
			freq[0]=(c1-c0)>>8&0xF0|1;
			unsigned c2=hist[sym], c3=sym<255?hist[sym+1]:0x10000;
			cdf[1]=(c2<<8)/freq[0]&0xF0|sym&15;
			freq[1]=((c3-c2)<<8)/freq[0]&0xF0|1;

			for(int kp=1;kp>=0;--kp)
			{
				if(state>freq[kp]<<8)//renorm
				{
					dlist_push_back1(&list, &state);
					state>>=8;
				}
				state=state/freq[kp]<<8|(cdf[kp]+state%freq[kp]);//update
			}
		}
		dlist_push_back(&list, &state, 2);
	}
	free(hist);

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
	return 1;
}