//util.c - Utilities
//Copyright (C) 2023  Ayman Wagih Mohsen, unless source link provided
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include"util.h"
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#include<sys/stat.h>
#include<errno.h>
#include<time.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#include<unistd.h>
#endif
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>//QueryPerformance...
#include<processthreadsapi.h>
#include<conio.h>
#include<psapi.h>
#else
#include<dirent.h>
#include<sys/types.h>
#include<sys/sysinfo.h>
#define sprintf_s	snprintf
#define vsprintf_s	vsnprintf
#ifndef _HUGE
#define _HUGE	HUGE_VAL
#endif
#endif
#if defined _MSC_VER && defined _WIN32
#include<process.h>
#define THREAD_CALL __stdcall
typedef unsigned THREAD_RET;
#else
#include<pthread.h>
#define THREAD_CALL
typedef void *THREAD_RET;
#endif
#include<immintrin.h>
static const char file[]=__FILE__;

char g_buf[G_BUF_SIZE]={0};

void memfill(void *dst, const void *src, size_t dstbytes, size_t srcbytes)
{
	size_t copied;
	char *d=(char*)dst;
	const char *s=(const char*)src;
	if(dstbytes<srcbytes)
	{
		memcpy(dst, src, dstbytes);
		return;
	}
	copied=srcbytes;
	memcpy(d, s, copied);
	while(copied<<1<=dstbytes)
	{
		memcpy(d+copied, d, copied);
		copied<<=1;
	}
	if(copied<dstbytes)
		memcpy(d+copied, d, dstbytes-copied);
}
void memswap_slow(void *p1, void *p2, size_t size)
{
	unsigned char *s1=(unsigned char*)p1, *s2=(unsigned char*)p2, *end=s1+size;
	for(;s1<end;++s1, ++s2)
	{
		const unsigned char t=*s1;
		*s1=*s2;
		*s2=t;
	}
}
void memswap(void *p1, void *p2, size_t size, void *temp)
{
	memcpy(temp, p1, size);
	memcpy(p1, p2, size);
	memcpy(p2, temp, size);
}
void memreverse(void *p, size_t count, size_t esize)
{
	size_t totalsize=count*esize;
	unsigned char *s1=(unsigned char*)p, *s2=s1+totalsize-esize;
	void *temp=malloc(esize);
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	while(s1<s2)
	{
		memswap(s1, s2, esize, temp);
		s1+=esize, s2-=esize;
	}
	free(temp);
}
void reverse16(void *start, void *end)
{
	__m256i reverse=_mm256_set_epi8(
		 1,  0,  3,  2,  5,  4,  7,  6,  9,  8, 11, 10, 13, 12, 15, 14,
		 1,  0,  3,  2,  5,  4,  7,  6,  9,  8, 11, 10, 13, 12, 15, 14
		// 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
		// 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15
	);
	unsigned short
		*p1=(unsigned short*)start,
		*p2=(unsigned short*)end;
	while(p1<p2-(sizeof(__m256i)/sizeof(short)-1))
	{
		p2-=sizeof(__m256i)/sizeof(short);
		__m256i v1=_mm256_loadu_si256((__m256i*)p1);
		__m256i v2=_mm256_loadu_si256((__m256i*)p2);
		v1=_mm256_shuffle_epi8(v1, reverse);
		v2=_mm256_shuffle_epi8(v2, reverse);
		v1=_mm256_permute2x128_si256(v1, v1, 1);
		v2=_mm256_permute2x128_si256(v2, v2, 1);
		_mm256_store_si256((__m256i*)p1, v2);
		_mm256_store_si256((__m256i*)p2, v1);
		p1+=sizeof(__m256i)/sizeof(short);
	}
	while(p1<p2-1)
	{
		--p2;
		unsigned short temp=*p1;
		*p1=*p2;
		*p2=temp;
		++p1;
	}
}
void memrotate(void *p, size_t byteoffset, size_t bytesize, void *temp)//temp buffer is min(byteoffset, bytesize-byteoffset)
{
	unsigned char *buf=(unsigned char*)p;

	if(byteoffset<bytesize-byteoffset)
	{
		memcpy(temp, buf, byteoffset);
		memmove(buf, buf+byteoffset, bytesize-byteoffset);
		memcpy(buf+bytesize-byteoffset, temp, byteoffset);
	}
	else
	{
		memcpy(temp, buf+byteoffset, bytesize-byteoffset);
		memmove(buf+bytesize-byteoffset, buf, byteoffset);
		memcpy(buf, temp, bytesize-byteoffset);
	}
}
int binary_search(const void *base, size_t count, size_t esize, int (*threeway)(const void*, const void*), const void *val, size_t *idx)
{
	const unsigned char *buf=(const unsigned char*)base;
#if 1
	ptrdiff_t low=0, range=count, mid;
	while(range)//binary search		log2(nlevels) memory accesses per symbol
	{
		ptrdiff_t floorhalf=range>>1;
		mid=low+floorhalf;
		low+=(range-floorhalf)&-(ptrdiff_t)(threeway(buf+mid*esize, val)<0);
		range=floorhalf;
	}
	if(idx)
		*idx=low;
	return !threeway(buf+low*esize, val);
#endif
#if 0
	ptrdiff_t left=0, right=(ptrdiff_t)count-1, mid;
	int ret;
	while(left<=right)
	{
		mid=(left+right)>>1;
		ret=threeway(buf+mid*esize, val);
		if(ret<0)
			left=mid+1;
		else if(ret>0)
			right=mid-1;
		else
		{
			if(idx)
				*idx=mid;
			return 1;
		}
	}
	if(idx)
		*idx=left+(left<(ptrdiff_t)count&&threeway(buf+left*esize, val)<0);
	return 0;
#endif
}
void isort(void *base, size_t count, size_t esize, int (*threeway)(const void*, const void*))
{
	unsigned char *buf=(unsigned char*)base;
	size_t k;
	void *temp;

	if(count<2)
		return;

	temp=malloc((count>>1)*esize);
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(k=1;k<count;++k)
	{
		size_t idx=0;
		binary_search(buf, k, esize, threeway, buf+k*esize, &idx);
		if(idx<k)
			memrotate(buf+idx*esize, (k-idx)*esize, (k+1-idx)*esize, temp);
	}
	free(temp);
}
int strcmp_ci(const char *s1, const char *s2)
{
#ifdef _MSC_VER
	return _stricmp(s1, s2);
#else
	return strcasecmp(s1, s2);
#endif
}
int acme_getopt(int argc, char **argv, int *start, const char **keywords, int kw_count)
{
	int k;
	size_t len;
	const char *arg, *cand;

	if(*start>=argc)
		return OPT_ENDOFARGS;
	
	arg=argv[*start];
	len=strlen(arg);
	if(len<=0)
		return OPT_INVALIDARG;
	//len>=1
	if(arg[0]!='-')
		return OPT_NOMATCH;
	++arg, --len;
	if(len<=0)
		return OPT_INVALIDARG;
	//len>=1
	if(arg[0]!='-')//short form (single dash followed by one character)
	{
		if(len!=1)
			return OPT_INVALIDARG;
		//len==1
		for(k=0;k<kw_count;++k)
		{
			cand=keywords[k];
			if(cand[0]==arg[0])
				return k;
		}
	}
	else//long form (double dash followed by a word)
	{
		++arg, --len;
		if(len<=0)
			return OPT_INVALIDARG;
		//len>=1
		for(k=0;k<kw_count;++k)
		{
			cand=keywords[k];
			if(!strcmp(arg, cand+1))
				return k;
		}
	}
	return OPT_NOMATCH;
}

int hammingweight16(unsigned short x)
{
#ifdef _MSC_VER
	return __popcnt16(x);
#else
	return __builtin_popcount(x);
#endif
}
int hammingweight32(unsigned x)
{
#ifdef _MSC_VER
	return __popcnt(x);
#else
	return __builtin_popcount(x);
#endif
}
int hammingweight64(unsigned long long x)
{
#ifdef _MSC_VER
	return (int)_mm_popcnt_u64(x);
#else
	return __builtin_popcountll(x);
#endif
}
#if 0
int floor_log2_p1(unsigned long long n)
{
#ifdef _MSC_VER
	return (sizeof(n)<<3)-(int)__lzcnt64(n);

	//unsigned long logn=0;
	//int success=_BitScanReverse64(&logn, n);
	//logn=success?logn+1:0;
	//return logn;
#elif defined __GNUC__
	int logn=64-__builtin_clzll(n);
	return logn;
#else
	int	logn=n!=0;
	int	sh=(n>=1ULL<<32)<<5;	logn+=sh, n>>=sh;
		sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
#endif
}
int floor_log2(unsigned long long n)
{
#if defined _MSC_VER || defined __GNUC__
	return (sizeof(n)<<3)-1-(int)_lzcnt_u64(n);//since Haswell
#else
	int	logn=-!n;
	int	sh=(n>=1ULL<<32)<<5;	logn+=sh, n>>=sh;
		sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
#endif
}
int floor_log2_32(unsigned n)
{
#if defined _MSC_VER || defined __GNUC__
	return (sizeof(n)<<3)-1-(int)_lzcnt_u32(n);//SSE4
#else
	//binary search
#if 0
	int	logn=-!n;
	int	sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
#endif

	//https://github.com/Cyan4973/FiniteStateEntropy/blob/d41d8be8e7955787ce486b90708d7b8de53137bd/lib/bitstream.h#L187
#if 1
	static const int table[]=
	{
		 0,  9,  1, 10, 13, 21,  2, 29,
		11, 14, 16, 18, 22, 25,  3, 30,
		 8, 12, 20, 28, 15, 17, 24,  7,
		19, 27, 23,  6, 26,  5,  4, 31,
	};
	n|=n>>1;
	n|=n>>2;
	n|=n>>4;
	n|=n>>8;
	n|=n>>16;
	n*=0x07C4ACDD;//0b0000_0111_1100_0100_1010_1100_1101_1101
	n>>=27;
	n=table[n];
	return n;
#endif

	//https://en.wikipedia.org/wiki/De_Bruijn_sequence?useskin=monobook
	//memory-bound, doesn't return -1 for zero
#if 0
	static const int table[]=
	{
		 0,  1, 16,  2, 29, 17,  3, 22, 30, 20, 18, 11, 13,  4,  7, 23,
		31, 15, 28, 21, 19, 10, 12,  6, 14, 27,  9,  5, 26,  8, 25, 24,
	};
	n|=n>>1;
	n|=n>>2;
	n|=n>>4;
	n|=n>>8;
	n|=n>>16;
	n-=n>>1;
	n*=0x06EB14F9;//0b0000_0110_1110_1011_0001_0100_1111_1001
	n>>=27;
	n=table[n];
	return n;
#endif
#endif
}
#endif
int ceil_log2(unsigned long long n)
{
	int lgn=FLOOR_LOG2(n);
	lgn+=(1ULL<<lgn)<n;
	return lgn;
}
int ceil_log2_32(unsigned n)
{
	int lgn=FLOOR_LOG2(n);
	lgn+=(1U<<lgn)<n;
	return lgn;
}
#if 0
int get_lsb_index(unsigned long long n)
{
#ifdef _MSC_VER
	return (int)_tzcnt_u64(n);//BMI1 (2013)

	//unsigned long idx;
	//_BitScanForward64(&idx, n);
	//return idx;
#elif defined __GNUC__
	return __builtin_ctzll(n);
	//return __builtin_ffsll(n)-1;
#else
	int cond, lsb;
	if(!n)
		return sizeof(n)<<3;
	lsb=0;
	cond=(n>>32<<32==n)<<5, lsb+=cond, n>>=cond;
	cond=(n>>16<<16==n)<<4, lsb+=cond, n>>=cond;
	cond=(n>> 8<< 8==n)<<3, lsb+=cond, n>>=cond;
	cond=(n>> 4<< 4==n)<<2, lsb+=cond, n>>=cond;
	cond=(n>> 2<< 2==n)<<1, lsb+=cond, n>>=cond;
	cond= n>> 1<< 1==n    , lsb+=cond;
	return lsb;
#endif
}
int get_lsb_index32(unsigned n)
{
#ifdef _MSC_VER
	return (int)_tzcnt_u32(n);

	//unsigned long idx;
	//_BitScanForward(&idx, n);
	//return idx;
#elif defined __GNUC__
	return __builtin_ctz(n);
	//return __builtin_ffs(n)-1;
#else
	int cond, lsb;
	if(!n)
		return sizeof(n)<<3;
	lsb=0;
	cond=(n>>16<<16==n)<<4, lsb+=cond, n>>=cond;
	cond=(n>> 8<< 8==n)<<3, lsb+=cond, n>>=cond;
	cond=(n>> 4<< 4==n)<<2, lsb+=cond, n>>=cond;
	cond=(n>> 2<< 2==n)<<1, lsb+=cond, n>>=cond;
	cond= n>> 1<< 1==n    , lsb+=cond;
	return lsb;
#endif
}
int get_lsb_index16(unsigned short n)
{
#ifdef _MSC_VER
	return (int)_tzcnt_u16(n);
#elif defined __GNUC__
	return __builtin_ctzs(n);
#else
	int cond, lsb;
	if(!n)
		return sizeof(n)<<3;
	lsb=0;
	cond=(n>> 8<< 8==n)<<3, lsb+=cond, n>>=cond;
	cond=(n>> 4<< 4==n)<<2, lsb+=cond, n>>=cond;
	cond=(n>> 2<< 2==n)<<1, lsb+=cond, n>>=cond;
	cond= n>> 1<< 1==n    , lsb+=cond;
	return lsb;
#endif
}
#endif
int floor_log10(double x)
{
	static const double pmask[]=//positive powers
	{
		1, 10,		//10^2^0
		1, 100,		//10^2^1
		1, 1e4,		//10^2^2
		1, 1e8,		//10^2^3
		1, 1e16,	//10^2^4
		1, 1e32,	//10^2^5
		1, 1e64,	//10^2^6
		1, 1e128,	//10^2^7
		1, 1e256	//10^2^8
	};
	static const double nmask[]=//negative powers
	{
		1, 0.1,		//1/10^2^0
		1, 0.01,	//1/10^2^1
		1, 1e-4,	//1/10^2^2
		1, 1e-8,	//1/10^2^3
		1, 1e-16,	//1/10^2^4
		1, 1e-32,	//1/10^2^5
		1, 1e-64,	//1/10^2^6
		1, 1e-128,	//1/10^2^7
		1, 1e-256	//1/10^2^8
	};
	int logn, sh;
	if(x<=0)
		return 0x80000000;
	if(x>=1)
	{
		logn=0;
		sh=(x>=pmask[17])<<8; logn+=sh, x*=nmask[16+(sh!=0)];
		sh=(x>=pmask[15])<<7; logn+=sh, x*=nmask[14+(sh!=0)];
		sh=(x>=pmask[13])<<6; logn+=sh, x*=nmask[12+(sh!=0)];
		sh=(x>=pmask[11])<<5; logn+=sh, x*=nmask[10+(sh!=0)];
		sh=(x>=pmask[9])<<4;  logn+=sh, x*=nmask[8+(sh!=0)];
		sh=(x>=pmask[7])<<3;  logn+=sh, x*=nmask[6+(sh!=0)];
		sh=(x>=pmask[5])<<2;  logn+=sh, x*=nmask[4+(sh!=0)];
		sh=(x>=pmask[3])<<1;  logn+=sh, x*=nmask[2+(sh!=0)];
		sh= x>=pmask[1];      logn+=sh;
		return logn;
	}
	logn=-1;
	sh=(x<nmask[17])<<8; logn-=sh;	x*=pmask[16+(sh!=0)];
	sh=(x<nmask[15])<<7; logn-=sh;	x*=pmask[14+(sh!=0)];
	sh=(x<nmask[13])<<6; logn-=sh;	x*=pmask[12+(sh!=0)];
	sh=(x<nmask[11])<<5; logn-=sh;	x*=pmask[10+(sh!=0)];
	sh=(x<nmask[9])<<4;  logn-=sh;	x*=pmask[8+(sh!=0)];
	sh=(x<nmask[7])<<3;  logn-=sh;	x*=pmask[6+(sh!=0)];
	sh=(x<nmask[5])<<2;  logn-=sh;	x*=pmask[4+(sh!=0)];
	sh=(x<nmask[3])<<1;  logn-=sh;	x*=pmask[2+(sh!=0)];
	sh= x<nmask[1];      logn-=sh;
	return logn;
}
unsigned floor_sqrt(unsigned long long x)
{
	unsigned long long low=0, range=x;
	while(range)
	{
		unsigned long long floorhalf=range>>1;
		unsigned long long level=low+floorhalf+1;
		if(level*level<=x)
			low+=range-floorhalf;
		range=floorhalf;
	}
	return (unsigned)low;
#if 0
	if(x<2)
		return (unsigned)x;
	int lg_sqrtx_p1=(floor_log2(x)>>1)+1;
	unsigned long long
		U=1ULL<<lg_sqrtx_p1,
		L=U>>1,
		temp=x>>(lg_sqrtx_p1-1);
	--U;
	U=MINVAR(U, temp);
	temp>>=1;
	L=MAXVAR(L, temp);
	//unsigned long long
	//	L=(1ULL<<lg_sqrtx_p1)>>1,
	//	U=(1ULL<<lg_sqrtx_p1)-1,
	//	L2=x>>lg_sqrtx_p1,	//where is the proof that these are also bounds?
	//	U2=x>>(lg_sqrtx_p1-1);
	//L=MAXVAR(L, L2);
	//U=MINVAR(U, U2);
	//int nmuls=0;//
	while(L<=U)//binary search
	{
		unsigned long long level=(L+U)>>1, sq=level*level;
		//++nmuls;//
		if(x>sq)
			L=level+1;
		else if(x<sq)
			U=level-1;
		else
		{
			U=level;
			break;
		}
	}
	//printf("%d muls\n", nmuls);//
	return (unsigned)U;//U <= L, we want floor, so return U
#endif

	//https://stackoverflow.com/questions/1100090/looking-for-an-efficient-integer-square-root-algorithm-for-arm-thumb2
#if 0
	unsigned res=0, add=0x80000000;
	int i;
	for(i=0;i<32;++i)
	{
		unsigned temp=res|add;
		if(x>=(unsigned long long)temp*temp)
			res=temp;
		add>>=1;
	}
	return res;
#endif
}
#define FRAC_BITS 24
unsigned exp2_fix24_neg(unsigned x)
{
	//return (unsigned)(exp2(-(x/16777216.))*0x1000000);//53% slower
	/*
	transcendental fractional powers of two
	x					inv(x)
	2^-0x0.000001 = 0x0.FFFFFF4F...		0x1.000000B1... = 2^0x0.000001
	2^-0x0.000002 = 0x0.FFFFFE9D...		0x1.00000163... = 2^0x0.000002
	2^-0x0.000004 = 0x0.FFFFFD3A...		0x1.000002C6... = 2^0x0.000004
	2^-0x0.000008 = 0x0.FFFFFA74...		0x1.0000058C... = 2^0x0.000008
	2^-0x0.000010 = 0x0.FFFFF4E9...		0x1.00000B17... = 2^0x0.000010
	2^-0x0.000020 = 0x0.FFFFE9D2...		0x1.0000162E... = 2^0x0.000020
	2^-0x0.000040 = 0x0.FFFFD3A3...		0x1.00002C5D... = 2^0x0.000040
	2^-0x0.000080 = 0x0.FFFFA747...		0x1.000058B9... = 2^0x0.000080
	2^-0x0.000100 = 0x0.FFFF4E8E...		0x1.0000B172... = 2^0x0.000100
	2^-0x0.000200 = 0x0.FFFE9D1D...		0x1.000162E5... = 2^0x0.000200
	2^-0x0.000400 = 0x0.FFFD3A3B...		0x1.0002C5CD... = 2^0x0.000400
	2^-0x0.000800 = 0x0.FFFA747F...		0x1.00058BA0... = 2^0x0.000800
	2^-0x0.001000 = 0x0.FFF4E91C...		0x1.000B175F... = 2^0x0.001000
	2^-0x0.002000 = 0x0.FFE9D2B3...		0x1.00162F39... = 2^0x0.002000
	2^-0x0.004000 = 0x0.FFD3A752...		0x1.002C605E... = 2^0x0.004000
	2^-0x0.008000 = 0x0.FFA75652...		0x1.0058C86E... = 2^0x0.008000
	2^-0x0.010000 = 0x0.FF4ECB59...		0x1.00B1AFA6... = 2^0x0.010000
	2^-0x0.020000 = 0x0.FE9E115C...		0x1.0163DAA0... = 2^0x0.020000
	2^-0x0.040000 = 0x0.FD3E0C0D...		0x1.02C9A3E7... = 2^0x0.040000
	2^-0x0.080000 = 0x0.FA83B2DB...		0x1.059B0D32... = 2^0x0.080000
	2^-0x0.100000 = 0x0.F5257D15...		0x1.0B5586D0... = 2^0x0.100000
	2^-0x0.200000 = 0x0.EAC0C6E8...		0x1.172B83C8... = 2^0x0.200000
	2^-0x0.400000 = 0x0.D744FCCB...		0x1.306FE0A3... = 2^0x0.400000
	2^-0x0.800000 = 0x0.B504F334...		0x1.6A09E667... = 2^0x0.800000
	*/
	static const unsigned long long frac_pots[]=
	{
		0x100000000,
		0x0FFFFFF4F,//extra 8 bits of precision
		0x0FFFFFE9D,
		0x0FFFFFD3A,
		0x0FFFFFA74,
		0x0FFFFF4E9,
		0x0FFFFE9D2,
		0x0FFFFD3A3,
		0x0FFFFA747,
		0x0FFFF4E8E,
		0x0FFFE9D1D,
		0x0FFFD3A3B,
		0x0FFFA747F,
		0x0FFF4E91C,
		0x0FFE9D2B3,
		0x0FFD3A752,
		0x0FFA75652,
		0x0FF4ECB59,
		0x0FE9E115C,
		0x0FD3E0C0D,
		0x0FA83B2DB,
		0x0F5257D15,
		0x0EAC0C6E8,
		0x0D744FCCB,
		0x0B504F334,
	};
#if 0
	unsigned long long x2=(unsigned long long)x<<1;
	unsigned long long r0=0x1000000, r1=0x1000000, r2=0x1000000, r3=0x1000000;
	for(int k=1;k<=FRAC_BITS;k+=8)
	{
		r0=r0*frac_pots[(k+0)&-(int)(x2>>(k+0)&1)]>>32;
		r1=r1*frac_pots[(k+1)&-(int)(x2>>(k+1)&1)]>>32;
		r2=r2*frac_pots[(k+2)&-(int)(x2>>(k+2)&1)]>>32;
		r3=r3*frac_pots[(k+3)&-(int)(x2>>(k+3)&1)]>>32;
		r0=r0*frac_pots[(k+4)&-(int)(x2>>(k+4)&1)]>>32;
		r1=r1*frac_pots[(k+5)&-(int)(x2>>(k+5)&1)]>>32;
		r2=r2*frac_pots[(k+6)&-(int)(x2>>(k+6)&1)]>>32;
		r3=r3*frac_pots[(k+7)&-(int)(x2>>(k+7)&1)]>>32;
	}
	r2=r2*r3>>24;
	r0=r0*r1>>24;
	r0=r0*r2>>24;
	r0>>=x>>FRAC_BITS;
	return (unsigned)r0;
#endif
#if 0
	unsigned long long x2=(unsigned long long)x<<1;
	unsigned long long r0=0x1000000, r1=0x1000000;
	for(int k=1;k<=FRAC_BITS;k+=2)
	{
		//unsigned long long t0=r0*frac_pots[k+0], t1=r1*frac_pots[k+1];//continuous access
		//t0>>=32;
		//t1>>=32;
		//if(x2>>(k+0)&1)
		r0=r0*frac_pots[(k+0)&-(int)(x2>>(k+0)&1)]>>32;
		r1=r1*frac_pots[(k+1)&-(int)(x2>>(k+1)&1)]>>32;
		//r0=r0*frac_pots[(k+2)&-(int)(x2>>(k+2)&1)]>>32;
		//r1=r1*frac_pots[(k+3)&-(int)(x2>>(k+3)&1)]>>32;
		//r0=r0*frac_pots[(k+4)&-(int)(x2>>(k+4)&1)]>>32;
		//r1=r1*frac_pots[(k+5)&-(int)(x2>>(k+5)&1)]>>32;
		//r0=r0*frac_pots[(k+6)&-(int)(x2>>(k+6)&1)]>>32;
		//r1=r1*frac_pots[(k+7)&-(int)(x2>>(k+7)&1)]>>32;
	}
	r0*=r1;
	r0>>=(x>>FRAC_BITS)+24;
	return (unsigned)r0;
#endif
#if 1
	unsigned long long result=0x1000000;
	for(int k=0;k<FRAC_BITS;)//up to 24 muls
	{
		int bit=x>>k&1;
		++k;
		result=result*frac_pots[k&-bit]>>32;
	}
	result>>=x>>FRAC_BITS;
	return (unsigned)result;
#endif
}
unsigned exp2_neg_fix24_avx2(unsigned x)
{
	/*
	transcendental fractional powers of two
	x					inv(x)
	2^-0x0.000001 = 0x0.FFFFFF4F...		0x1.000000B1... = 2^0x0.000001
	2^-0x0.000002 = 0x0.FFFFFE9D...		0x1.00000163... = 2^0x0.000002
	2^-0x0.000004 = 0x0.FFFFFD3A...		0x1.000002C6... = 2^0x0.000004
	2^-0x0.000008 = 0x0.FFFFFA74...		0x1.0000058C... = 2^0x0.000008
	2^-0x0.000010 = 0x0.FFFFF4E9...		0x1.00000B17... = 2^0x0.000010
	2^-0x0.000020 = 0x0.FFFFE9D2...		0x1.0000162E... = 2^0x0.000020
	2^-0x0.000040 = 0x0.FFFFD3A3...		0x1.00002C5D... = 2^0x0.000040
	2^-0x0.000080 = 0x0.FFFFA747...		0x1.000058B9... = 2^0x0.000080
	2^-0x0.000100 = 0x0.FFFF4E8E...		0x1.0000B172... = 2^0x0.000100
	2^-0x0.000200 = 0x0.FFFE9D1D...		0x1.000162E5... = 2^0x0.000200
	2^-0x0.000400 = 0x0.FFFD3A3B...		0x1.0002C5CD... = 2^0x0.000400
	2^-0x0.000800 = 0x0.FFFA747F...		0x1.00058BA0... = 2^0x0.000800
	2^-0x0.001000 = 0x0.FFF4E91C...		0x1.000B175F... = 2^0x0.001000
	2^-0x0.002000 = 0x0.FFE9D2B3...		0x1.00162F39... = 2^0x0.002000
	2^-0x0.004000 = 0x0.FFD3A752...		0x1.002C605E... = 2^0x0.004000
	2^-0x0.008000 = 0x0.FFA75652...		0x1.0058C86E... = 2^0x0.008000
	2^-0x0.010000 = 0x0.FF4ECB59...		0x1.00B1AFA6... = 2^0x0.010000
	2^-0x0.020000 = 0x0.FE9E115C...		0x1.0163DAA0... = 2^0x0.020000
	2^-0x0.040000 = 0x0.FD3E0C0D...		0x1.02C9A3E7... = 2^0x0.040000
	2^-0x0.080000 = 0x0.FA83B2DB...		0x1.059B0D32... = 2^0x0.080000
	2^-0x0.100000 = 0x0.F5257D15...		0x1.0B5586D0... = 2^0x0.100000
	2^-0x0.200000 = 0x0.EAC0C6E8...		0x1.172B83C8... = 2^0x0.200000
	2^-0x0.400000 = 0x0.D744FCCB...		0x1.306FE0A3... = 2^0x0.400000
	2^-0x0.800000 = 0x0.B504F334...		0x1.6A09E667... = 2^0x0.800000
	*/
	ALIGN(32) static const unsigned long long frac_pots[]=
	{
		(0xFFFFFF4F+128)>>8,//bit  0
		(0xFFFFD3A3+128)>>8,//bit  6
		(0xFFF4E91C+128)>>8,//bit 12
		(0xFD3E0C0D+128)>>8,//bit 18
		
		0xFFFFFE9D,//bit  1
		0xFFFFA747,//bit  7
		0xFFE9D2B3,//bit 13
		0xFA83B2DB,//bit 19
		
		0xFFFFFD3A,//bit  2
		0xFFFF4E8E,//bit  8
		0xFFD3A752,//bit 14
		0xF5257D15,//bit 20
		
		0xFFFFFA74,//bit  3
		0xFFFE9D1D,//bit  9
		0xFFA75652,//bit 15
		0xEAC0C6E8,//bit 21
		
		0xFFFFF4E9,//bit  4
		0xFFFD3A3B,//bit 10
		0xFF4ECB59,//bit 16
		0xD744FCCB,//bit 22
		
		0xFFFFE9D2,//bit  5
		0xFFFA747F,//bit 11
		0xFE9E115C,//bit 17
		0xB504F334,//bit 23

		0x1000000,//initialization
		0x1000000,
		0x1000000,
		0x1000000,
	};
	__m256i sr32=_mm256_set_epi8(
		-1, -1, -1, -1, 15, 14, 13, 12,
		-1, -1, -1, -1,  7,  6,  5,  4,
		-1, -1, -1, -1, 15, 14, 13, 12,
		-1, -1, -1, -1,  7,  6,  5,  4
		//15, 14, 13, 12, 11, 10,  9,  8,
		// 7,  6,  5,  4,  3,  2,  1,  0
	);
	__m128i sr24=_mm_set_epi8(
		-1, -1, -1, 15, 14, 13, 12, 11,
		-1, -1, -1,  7,  6,  5,  4,  3
	);
	__m256i ones=_mm256_set_epi32(0, 1, 0, 1, 0, 1, 0, 1);

	__m256i x2=_mm256_set_epi32(0, x>>18&63, 0, x>>12&63, 0, x>>6&63, 0, x&63);
	__m256i result=_mm256_load_si256((const __m256i*)frac_pots+6);

	__m256i factor=_mm256_load_si256((const __m256i*)frac_pots+0);
	__m256i mask=_mm256_and_si256(x2, ones);
	x2=_mm256_srli_epi32(x2, 1);
	mask=_mm256_cmpeq_epi64(mask, ones);
	result=_mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(result), _mm256_castsi256_pd(factor), _mm256_castsi256_pd(mask)));
	for(int k=1;k<FRAC_BITS/4;++k)
	{
		factor=_mm256_load_si256((const __m256i*)frac_pots+k);
		factor=_mm256_mul_epu32(factor, result);
		mask=_mm256_and_si256(x2, ones);
		x2=_mm256_srli_epi32(x2, 1);
		factor=_mm256_shuffle_epi8(factor, sr32);
		mask=_mm256_cmpeq_epi64(mask, ones);
		result=_mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(result), _mm256_castsi256_pd(factor), _mm256_castsi256_pd(mask)));
	}
	__m128i r0=_mm256_extracti128_si256(result, 0);
	__m128i r1=_mm256_extracti128_si256(result, 1);
	r0=_mm_mul_epu32(r0, r1);
	r0=_mm_shuffle_epi8(r0, sr24);
	unsigned sh=(x>>FRAC_BITS)+24;
	unsigned res=(unsigned)((unsigned long long)_mm_extract_epi64(r0, 0)*(unsigned long long)_mm_extract_epi64(r0, 1)>>sh);
	return res;
}
unsigned long long exp2_fix24(int x)
{
	/*
	transcendental fractional powers of two
	x					inv(x)
	2^-0x0.000001 = 0x0.FFFFFF4F...		0x1.000000B1... = 2^0x0.000001
	2^-0x0.000002 = 0x0.FFFFFE9D...		0x1.00000163... = 2^0x0.000002
	2^-0x0.000004 = 0x0.FFFFFD3A...		0x1.000002C6... = 2^0x0.000004
	2^-0x0.000008 = 0x0.FFFFFA74...		0x1.0000058C... = 2^0x0.000008
	2^-0x0.000010 = 0x0.FFFFF4E9...		0x1.00000B17... = 2^0x0.000010
	2^-0x0.000020 = 0x0.FFFFE9D2...		0x1.0000162E... = 2^0x0.000020
	2^-0x0.000040 = 0x0.FFFFD3A3...		0x1.00002C5D... = 2^0x0.000040
	2^-0x0.000080 = 0x0.FFFFA747...		0x1.000058B9... = 2^0x0.000080
	2^-0x0.000100 = 0x0.FFFF4E8E...		0x1.0000B172... = 2^0x0.000100
	2^-0x0.000200 = 0x0.FFFE9D1D...		0x1.000162E5... = 2^0x0.000200
	2^-0x0.000400 = 0x0.FFFD3A3B...		0x1.0002C5CD... = 2^0x0.000400
	2^-0x0.000800 = 0x0.FFFA747F...		0x1.00058BA0... = 2^0x0.000800
	2^-0x0.001000 = 0x0.FFF4E91C...		0x1.000B175F... = 2^0x0.001000
	2^-0x0.002000 = 0x0.FFE9D2B3...		0x1.00162F39... = 2^0x0.002000
	2^-0x0.004000 = 0x0.FFD3A752...		0x1.002C605E... = 2^0x0.004000
	2^-0x0.008000 = 0x0.FFA75652...		0x1.0058C86E... = 2^0x0.008000
	2^-0x0.010000 = 0x0.FF4ECB59...		0x1.00B1AFA6... = 2^0x0.010000
	2^-0x0.020000 = 0x0.FE9E115C...		0x1.0163DAA0... = 2^0x0.020000
	2^-0x0.040000 = 0x0.FD3E0C0D...		0x1.02C9A3E7... = 2^0x0.040000
	2^-0x0.080000 = 0x0.FA83B2DB...		0x1.059B0D32... = 2^0x0.080000
	2^-0x0.100000 = 0x0.F5257D15...		0x1.0B5586D0... = 2^0x0.100000
	2^-0x0.200000 = 0x0.EAC0C6E8...		0x1.172B83C8... = 2^0x0.200000
	2^-0x0.400000 = 0x0.D744FCCB...		0x1.306FE0A3... = 2^0x0.400000
	2^-0x0.800000 = 0x0.B504F334...		0x1.6A09E667... = 2^0x0.800000
	*/
	static const unsigned long long frac_pots[]=
	{
		0x1000000B1,//round(2^(0x000001*2^-24)*2^32)		//extra 8 bits of precision
		0x100000163,
		0x1000002C6,
		0x10000058C,
		0x100000B17,
		0x10000162E,
		0x100002C5D,
		0x1000058B9,
		0x10000B172,
		0x1000162E5,
		0x10002C5CC,
		0x100058BA0,
		0x1000B175F,
		0x100162F39,
		0x1002C605E,
		0x10058C86E,
		0x100B1AFA6,
		0x10163DAA0,
		0x102C9A3E7,
		0x1059B0D31,
		0x10B5586D0,
		0x1172B83C8,
		0x1306FE0A3,
		0x16A09E668,//round(2^(0x800000*2^-24)*2^32)
	};
	int neg=x<0;
	x=abs(x);
	if(x>=0x27000000)
		return neg?0:0xFFFFFFFFFFFFFFFF;
	unsigned long long result=0x1000000;
	for(int k=0;k<FRAC_BITS;++k)//up to 24 muls
	{
		if(x&1)
		{
			result*=frac_pots[k];
			result>>=32;
		}
		x>>=1;
	}
	result=SHIFT_LEFT_SIGNED(result, x);
	if(neg)
		result=0x1000000000000/result;
	return result;
}
int log2_fix24(unsigned long long x)
{
	//https://stackoverflow.com/questions/4657468/fast-fixed-point-pow-log-exp-and-sqrt
	//and YouChad
	if(!x)
		return -((FRAC_BITS+1)<<FRAC_BITS);
	int lgx=FLOOR_LOG2(x)-24;
	x=SHIFT_RIGHT_SIGNED(x, lgx);
	lgx<<=24;
	if(!(x&(x-1)))//no need to do 24 muls if x is a power of two
		return lgx;
	for(int k=FRAC_BITS-1;k>=0;--k)//constant 24 muls, hopefully will be unrolled
	{
		x=x*x>>FRAC_BITS;
		int cond=x>=(2LL<<FRAC_BITS);
		x>>=cond;
		lgx|=cond<<k;
	}
	return lgx;
}
#undef  FRAC_BITS
double power(double x, int y)
{
	double mask[]={1, 0}, product=1;
	if(y<0)
		mask[1]=1/x, y=-y;
	else
		mask[1]=x;
	for(;;)
	{
		product*=mask[y&1], y>>=1;	//67.7
		if(!y)
			return product;
		mask[1]*=mask[1];
	}
	return product;
}
double _10pow(int n)
{
	static double mask[616]={0};
//	const double _ln10=log(10.);
	if(!mask[308])
	{
		int k;
		for(k=-308;k<308;++k)		//23.0
			mask[k+308]=power(10., k);
		//	mask[k+308]=exp(k*_ln10);//inaccurate
	}
	if(n<-308)
		return 0;
	if(n>307)
		return _HUGE;
	return mask[n+308];
}
int acme_isdigit(char c, char base)
{
	switch(base)
	{
	case 2:		return BETWEEN_INC('0', c, '1');
	case 8:		return BETWEEN_INC('0', c, '7');
	case 10:	return BETWEEN_INC('0', c, '9');
	case 16:	return BETWEEN_INC('0', c, '9')||BETWEEN_INC('A', c&0xDF, 'F');
	}
	return 0;
}
#ifdef __GNUC__
unsigned long long _udiv128(unsigned long long hi, unsigned long long lo, unsigned long long den, unsigned long long *rem)
{
	//FIXME get rid of __int128
	//https://stackoverflow.com/questions/1870158/unsigned-128-bit-division-on-64-bit-machine
	unsigned __int128 num=(unsigned __int128)hi<<64|lo;
	*rem=(unsigned long long)(num%den);
	return (unsigned long long)(num/den);
}
unsigned long long _umul128(unsigned long long a, unsigned long long b, unsigned long long *hi)
{
	unsigned __int128 num=(unsigned __int128)a*b;
	*hi=(unsigned long long)(num>>64);
	return (unsigned long long)num;
}
#endif

double time_ms(void)
{
#ifdef _MSC_VER
	static long long t0=0;
	LARGE_INTEGER li;
	double t;
	QueryPerformanceCounter(&li);
	if(!t0)
		t0=li.QuadPart;
	t=(double)(li.QuadPart-t0);
	QueryPerformanceFrequency(&li);
	t/=(double)li.QuadPart;
	t*=1000;
	return t;
#else
	struct timespec t;
	clock_gettime(CLOCK_REALTIME, &t);//<time.h>
	return t.tv_sec*1000+t.tv_nsec*1e-6;
#endif
}
double time_sec(void)
{
#ifdef _MSC_VER
	static long long t0=0;
	LARGE_INTEGER li;
	double t;
	QueryPerformanceCounter(&li);
	if(!t0)
		t0=li.QuadPart;
	t=(double)(li.QuadPart-t0);
	QueryPerformanceFrequency(&li);
	t/=(double)li.QuadPart;
	return t;
#else
	struct timespec t;
	clock_gettime(CLOCK_REALTIME, &t);//<time.h>
	return t.tv_sec+t.tv_nsec*1e-9;
#endif
}
void parsetimedelta(double secs, TimeInfo *ti)
{
	ti->days=(int)floor(secs/(60*60*24));
	secs-=ti->days*(60*60*24);

	ti->hours=(int)floor(secs/(60*60));
	secs-=ti->hours*(60*60);

	ti->mins=(int)floor(secs/60);
	secs-=ti->mins*60;

	ti->secs=(float)secs;
}
int timedelta2str(char *buf, size_t len, double secs)
{
	int printed;
	TimeInfo ti;

	parsetimedelta(secs, &ti);

	printed=0;
	if(buf)
	{
		if(ti.days)
			printed+=snprintf(buf, len, "%02dD-", ti.days);
		printed+=snprintf(buf, len, "%02d-%02d-%09.6lf", ti.hours, ti.mins, ti.secs);
	}
	else
	{
		if(ti.days)
			printed+=printf("%02dD-", ti.days);
		printed+=printf("%02d-%02d-%09.6lf", ti.hours, ti.mins, ti.secs);
	}
	return printed;
}
int acme_strftime(char *buf, size_t len, const char *format)
{
#ifdef _MSC_VER
	time_t tstamp;
	struct tm tformat;

	tstamp=time(0);
	localtime_s(&tformat, &tstamp);
	return (int)strftime(buf, len, format, &tformat);
#else
	time_t tstamp=time(0);
	struct tm *tformat=localtime(&tstamp);
	return (int)strftime(buf, len, format, tformat);
#endif
}
int print_bin8(int x)
{
	//printf("0b");
	for(int k=7;k>=0;--k)
	{
		int bit=x>>k&1;
		printf("%c", '0'+bit);
	}
	return 34;
}
int print_bin32(unsigned x)
{
	//printf("0b");
	for(int k=31;k>=0;--k)
	{
		int bit=x>>k&1;
		printf("%c", '0'+bit);
	}
	return 34;
}
int print_binn(unsigned long long x, int nbits)
{
	for(int k=nbits-1;k>=0;--k)
	{
		int bit=x>>k&1;
		printf("%c", '0'+bit);
	}
	return nbits;
}
void print_nan(double x, int total, int decimal)
{
	if(x!=x)
		printf("%*s", total, "");
	else
		printf("%*.*lf", total, decimal, x);
}
double convert_size(double bytesize, int *log1024)
{
	int k=0;
	double p=1;
	for(;k<4&&bytesize>p;++k, p*=1024);
	*log1024=k-1;
	return bytesize*1024/p;
}
int print_size(double bytesize, int ndigits, int pdigits, char *str, int len)
{
	static const char *labels[]=
	{
		"bytes",
		"KB",
		"MB",
		"GB",
		"TB",
		"PB",
		"EB",
	};
	int idx=0;
	bytesize=convert_size(bytesize, &idx);
	if(str)
		return snprintf(str, len-1LL, "%*.*lf %s", ndigits, pdigits, bytesize, labels[idx]);
	return printf("%*.*lf %s", ndigits, pdigits, bytesize, labels[idx]);
}

//error handling
char first_error_msg[G_BUF_SIZE]={0}, latest_error_msg[G_BUF_SIZE]={0};
int log_error(const char *fn, int line, int quit, const char *format, ...)
{
	int firsttime=first_error_msg[0]=='\0';

	ptrdiff_t size=strlen(fn), start=size-1;
	for(;start>=0&&fn[start]!='/'&&fn[start]!='\\';--start);
	start+=start==-1||fn[start]=='/'||fn[start]=='\\';

	int printed=sprintf_s(latest_error_msg, G_BUF_SIZE, "\n%s(%d): ", fn+start, line);
	va_list args;
	va_start(args, format);
	printed+=vsprintf_s(latest_error_msg+printed, G_BUF_SIZE-printed-1, format, args);
	va_end(args);

	if(firsttime)
		memcpy(first_error_msg, latest_error_msg, printed+1LL);
	fprintf(stderr, "%s\n", latest_error_msg);
	//messagebox(MBOX_OK, "Error", latest_error_msg);
	if(quit)
	{
		pause();
		exit(0);
	}
	return firsttime;
}
int pause(void)
{
	int k;

	printf("Enter 0 to continue: ");
	while(!scanf("%d", &k));
	return k;
}
int pause_abort(const char *fn, int lineno, const char *extraInfo)
{
	printf("INTERNAL ERROR %s(%d)\nABORTING\n", fn, lineno);
	if(extraInfo)
		printf("%s\n\n", extraInfo);
	pause();
	abort();
	return 0;
}


//C array
#if 1
static void array_realloc(ArrayHandle *arr, size_t count, size_t pad)//CANNOT be nullptr, array must be initialized with array_alloc()
{
	size_t size, newcap;

	size=(count+pad)*arr[0]->esize, newcap=arr[0]->esize;
	for(;newcap<size;newcap<<=1);
	if(newcap>arr[0]->cap)
	{
		void *p2=realloc(*arr, sizeof(ArrayHeader)+newcap);
		if(!p2)
		{
			LOG_ERROR("Alloc error %016llX", sizeof(ArrayHeader)+newcap);
			return;
		}
		*arr=(ArrayHandle)p2;
		if(arr[0]->cap<newcap)
			memset(arr[0]->data+arr[0]->cap, 0, newcap-arr[0]->cap);
		arr[0]->cap=newcap;
	}
	arr[0]->count=count;
}

//Array API
ArrayHandle array_construct(const void *src, size_t esize, size_t count, size_t rep, size_t pad, void (*destructor)(void*))
{
	ArrayHandle arr;
	size_t srcsize, dstsize, cap;
	
	srcsize=count*esize;
	dstsize=rep*srcsize;
	cap=dstsize+pad*esize;
	arr=(ArrayHandle)malloc(sizeof(ArrayHeader)+cap);
	if(!arr)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	arr->count=count*rep;
	arr->esize=esize;
	arr->cap=cap;
	arr->destructor=destructor;
	if(src)
		memfill(arr->data, src, dstsize, srcsize);
	else
		memset(arr->data, 0, dstsize);
		
	if(cap-dstsize>0)//zero pad
		memset(arr->data+dstsize, 0, cap-dstsize);
	return arr;
}
ArrayHandle array_copy(ArrayHandle *arr)
{
	ArrayHandle a2;
	size_t bytesize;

	if(!*arr)
		return 0;
	bytesize=sizeof(ArrayHeader)+arr[0]->cap;
	a2=(ArrayHandle)malloc(bytesize);
	if(!a2)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memcpy(a2, *arr, bytesize);
	return a2;
}
void array_clear(ArrayHandle *arr)//can be nullptr
{
	if(*arr)
	{
		if(arr[0]->destructor)
		{
			for(size_t k=0;k<arr[0]->count;++k)
				arr[0]->destructor(array_at(arr, k));
		}
		arr[0]->count=0;
	}
}
void array_free(ArrayHandle *arr)//can be nullptr
{
	if(*arr&&arr[0]->destructor)
	{
		for(size_t k=0;k<arr[0]->count;++k)
			arr[0]->destructor(array_at(arr, k));
	}
	free(*arr);
	*arr=0;
}
void array_fit(ArrayHandle *arr, size_t pad)//can be nullptr
{
	void *p2;

	if(!*arr)
		return;
	arr[0]->cap=(arr[0]->count+pad)*arr[0]->esize;
	p2=realloc(*arr, sizeof(ArrayHeader)+arr[0]->cap);
	if(!p2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	*arr=(ArrayHandle)p2;
}

void* array_insert(ArrayHandle *arr, size_t idx, const void *data, size_t count, size_t rep, size_t pad)//cannot be nullptr
{
	size_t start, srcsize, dstsize, movesize;
	
	start=idx*arr[0]->esize;
	srcsize=count*arr[0]->esize;
	dstsize=rep*srcsize;
	movesize=arr[0]->count*arr[0]->esize-start;
	array_realloc(arr, arr[0]->count+rep*count, pad);
	memmove(arr[0]->data+start+dstsize, arr[0]->data+start, movesize);
	if(data)
		memfill(arr[0]->data+start, data, dstsize, srcsize);
	else
		memset(arr[0]->data+start, 0, dstsize);
	return arr[0]->data+start;
}
void* array_erase(ArrayHandle *arr, size_t idx, size_t count)//does not reallocate
{
	size_t k;

	if(arr[0]->count<idx+count)
	{
		LOG_ERROR("array_erase() out of bounds: idx=%lld count=%lld size=%lld", (long long)idx, (long long)count, (long long)arr[0]->count);
		if(arr[0]->count<idx)
			return 0;
		count=arr[0]->count-idx;//erase till end of array if just idx+count is OOB
	}
	if(arr[0]->destructor)
	{
		for(k=0;k<count;++k)
			arr[0]->destructor(array_at(arr, idx+k));
	}
	memmove(arr[0]->data+idx*arr[0]->esize, arr[0]->data+(idx+count)*arr[0]->esize, (arr[0]->count-(idx+count))*arr[0]->esize);
	arr[0]->count-=count;
	return arr[0]->data+idx*arr[0]->esize;
}
void* array_replace(ArrayHandle *arr, size_t idx, size_t rem_count, const void *data, size_t ins_count, size_t rep, size_t pad)
{
	size_t k, c0, c1, start, srcsize, dstsize;

	if(arr[0]->count<idx+rem_count)
	{
		LOG_ERROR("array_replace() out of bounds: idx=%lld rem_count=%lld size=%lld ins_count=%lld", (long long)idx, (long long)rem_count, (long long)arr[0]->count, (long long)ins_count);
		if(arr[0]->count<idx)
			return 0;
		rem_count=arr[0]->count-idx;//erase till end of array if just idx+count is OOB
	}
	if(arr[0]->destructor)
	{
		for(k=0;k<rem_count;++k)//destroy removed objects
			arr[0]->destructor(array_at(arr, idx+k));
	}
	start=idx*arr[0]->esize;
	srcsize=ins_count*arr[0]->esize;
	dstsize=rep*srcsize;
	c0=arr[0]->count;//copy original count
	c1=arr[0]->count+rep*ins_count-rem_count;//calculate new count

	if(ins_count!=rem_count||(c1+pad)*arr[0]->esize>arr[0]->cap)//resize array
		array_realloc(arr, c1, pad);

	if(ins_count!=rem_count)//shift objects
		memmove(arr[0]->data+(idx+ins_count)*arr[0]->esize, arr[0]->data+(idx+rem_count)*arr[0]->esize, (c0-(idx+rem_count))*arr[0]->esize);

	if(dstsize)//initialize inserted range
	{
		if(data)
			memfill(arr[0]->data+start, data, dstsize, srcsize);
		else
			memset(arr[0]->data+start, 0, dstsize);
	}
	return arr[0]->data+start;//return start of inserted range
}

void* array_at(ArrayHandle *arr, size_t idx)
{
	if(!arr[0])
	{
		LOG_ERROR("nullptr exception");
		return 0;
	}
	if(idx>=arr[0]->count)
	{
		LOG_ERROR("OOB");
		return 0;
	}
	return arr[0]->data+idx*arr[0]->esize;
}
void* array_back(ArrayHandle *arr)
{
	if(!*arr||!arr[0]->count)
		return 0;
	return arr[0]->data+(arr[0]->count-1)*arr[0]->esize;
}

int str_append(ArrayHandle *str, const char *format, ...)
{
	size_t reqlen;
	va_list args;
	va_start(args, format);
	reqlen=vsnprintf(0, 0, format, args);//requires C99
	if(str[0]->count+reqlen+1>str[0]->cap)
	{
		size_t c0=str[0]->count;
		array_realloc(str, str[0]->count+reqlen, 1);
		str[0]->count=c0;
	}
	reqlen=vsnprintf((char*)str[0]->data+str[0]->count, str[0]->cap-str[0]->count, format, args);
	str[0]->count+=reqlen;
	va_end(args);
	return (int)reqlen;
}

size_t array_append(ArrayHandle *dst, const void *src, size_t esize, size_t count, size_t rep, size_t pad, void (*destructor)(void*))//arr can be 0, returns original array size
{
	size_t dststart=0;
	if(!*dst)
		*dst=array_construct(src, esize, count, rep, pad, destructor);
	else
	{
		dststart=dst[0]->count*dst[0]->esize;
		if(dst[0]->esize!=esize)
			LOG_ERROR("Array element size mismatch");
		else
			ARRAY_APPEND(*dst, src, count, rep, pad);
	}
	return dststart;
}
#endif


//double-linked array list
#if 1
void dlist_init(DListHandle list, size_t objsize, size_t objpernode, void (*destructor)(void*))
{
	list->i=list->f=0;
	list->objsize=objsize;
	list->objpernode=objpernode;
	list->nnodes=list->nobj=0;//empty
	list->destructor=destructor;
}
#define DLIST_COPY_NODE(DST, PREV, NEXT, SRC, PAYLOADSIZE)\
	do\
	{\
		DST=(DNodeHandle)malloc(sizeof(DNodeHeader)+(PAYLOADSIZE));\
		if(!(DST))\
		{\
			LOG_ERROR("Alloc error");\
			return;\
		}\
		DST->prev=PREV;\
		DST->next=NEXT;\
		memcpy(DST->data, SRC->data, PAYLOADSIZE);\
	}\
	while(0)
void dlist_copy(DListHandle dst, DListHandle src)
{
	DNodeHandle it;
	size_t payloadsize;

	dlist_init(dst, src->objsize, src->objpernode, src->destructor);
	it=dst->i;
	if(it)
	{
		payloadsize=src->objpernode*src->objsize;

		DLIST_COPY_NODE(dst->f, 0, 0, it, payloadsize);
		//dst->f=(DNodeHandle)malloc(sizeof(DNodeHeader)+payloadsize);
		//memcpy(dst->f->data, it->data, payloadsize);

		dst->i=dst->f;

		it=it->next;

		for(;it;it=it->next)
		{
			DLIST_COPY_NODE(dst->f->next, dst->f, 0, it, payloadsize);
			dst->f=dst->f->next;
		}
	}
	dst->nnodes=src->nnodes;
	dst->nobj=src->nobj;
}
void dlist_clear(DListHandle list)
{
	DNodeHandle it;

	it=list->i;
	if(it)
	{
		while(it->next)
		{
			if(list->destructor)
			{
				for(size_t k=0;k<list->objpernode;++k)
					list->destructor(it->data+k*list->objsize);
				list->nobj-=list->objpernode;
			}
			it=it->next;
			free(it->prev);
		}
		if(list->destructor)
		{
			for(size_t k=0;k<list->nobj;++k)
				list->destructor(it->data+k*list->objsize);
		}
		free(it);
		list->i=list->f=0;
		list->nobj=list->nnodes=0;
	}
}
size_t dlist_appendtoarray(DListHandle list, ArrayHandle *dst)
{
	DNodeHandle it;
	size_t start, payloadsize;

	if(!*dst)
	{
		start=0;
		*dst=array_construct(0, list->objsize, 0, 0, list->nnodes*list->objpernode, list->destructor);
	}
	else
	{
		if(dst[0]->esize!=list->objsize)
		{
			LOG_ERROR("dlist_appendtoarray(): dst->esize=%d, list->objsize=%d", dst[0]->esize, list->objsize);
			return 0;
		}
		start=dst[0]->count;
		ARRAY_APPEND(*dst, 0, 0, 0, list->nnodes*list->objpernode);
	}
	it=list->i;
	payloadsize=list->objpernode*list->objsize;
	for(size_t offset=dst[0]->count;it;)
	{
		memcpy(dst[0]->data+offset*list->objsize, it->data, payloadsize);
		offset+=list->objpernode;
		it=it->next;
	}
	dst[0]->count+=list->nobj;
	return start;
}

static void dlist_append_node(DListHandle list)
{
	DNodeHandle temp;

	temp=(DNodeHandle)malloc(sizeof(DNodeHeader)+list->objpernode*list->objsize);
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	temp->next=0;
	if(list->nnodes)
	{
		temp->prev=list->f;
		list->f->next=temp;
	}
	else
	{
		temp->prev=0;
		list->i=temp;
	}
	list->f=temp;
	++list->nnodes;
}
#if 0
#define dlist_fill_node(LIST, COPYSIZE, SRC, DST)\
	LIST->nobj+=COPYSIZE;\
	COPYSIZE*=LIST->objsize;\
	if(SRC)\
		memcpy(DST, SRC, COPYSIZE), SRC+=COPYSIZE;\
	else\
		memset(DST, 0, COPYSIZE);
#else
static void dlist_fill_node(DListHandle list, size_t copysize, const char **src, void *dst)
{
	list->nobj+=copysize;
	copysize*=list->objsize;

	if(*src)
		memcpy(dst, *src, copysize), *src+=copysize;
	else
		memset(dst, 0, copysize);
}
#endif
static void dlist_rdestroy(DListHandle list, ptrdiff_t rstart, ptrdiff_t rend)
{
	if(list->destructor)
	{
		--rstart;
		rstart*=list->objsize;
		rend  *=list->objsize;
		for(;rstart>=rend;rstart-=list->objsize)
			list->destructor(list->f->data+rstart);
	}
}
void* dlist_push_back1(DListHandle list, const void *obj)
{
	size_t obj_idx=list->nobj%list->objpernode;//index of next object
	if(!obj_idx)//need a new node
		dlist_append_node(list);
	void *p=list->f->data+obj_idx*list->objsize;
	if(obj)
		memcpy(p, obj, list->objsize);
	else
		memset(p, 0, list->objsize);
	++list->nobj;
	return p;
}
void* dlist_push_back(DListHandle list, const void *data, size_t count)
{
	size_t obj_idx, copysize;
	const char *buffer;
	void *ret;
	
	buffer=(const char*)data;
	ret=0;
	obj_idx=list->nobj%list->objpernode;
	if(obj_idx)
	{
		copysize=list->objpernode<obj_idx+count?list->objpernode-obj_idx:count;
		count-=copysize;
		ret=list->f->data+obj_idx*list->objsize;

		dlist_fill_node(list, copysize, &buffer, ret);
	}
	while(count)
	{
		dlist_append_node(list);
		
		copysize=list->objpernode<count?list->objpernode:count;
		count-=copysize;

		if(!ret)
			ret=list->f->data;
		dlist_fill_node(list, copysize, &buffer, list->f->data);
	}
	return ret;
}
void* dlist_back(DListHandle list)
{
	size_t obj_idx;

	if(!list->nobj)
	{
		LOG_ERROR("dlist_back() called on empty list");
		return 0;
	}
	obj_idx=(list->nobj-1)%list->objpernode;
	return list->f->data+obj_idx*list->objsize;
}
void  dlist_pop_back1(DListHandle list)
{
	size_t obj_idx;

	if(!list->nobj)
		LOG_ERROR("dlist_pop_back1() called on empty list");
	if(list->destructor)
		list->destructor(dlist_back(list));
	obj_idx=(list->nobj-1)%list->objpernode;
	if(!obj_idx)//last object is first in the last block
	{
		DNodeHandle last=list->f;
		list->f=last->prev;
		free(last);
		--list->nnodes;
		if(!list->nnodes)//last object was popped out
			list->i=0;
	}
	--list->nobj;
}
void  dlist_pop_back(DListHandle list, size_t count)
{
	DNodeHandle last;
	size_t
		q1, r1,
		q2, r2;

	if(list->nobj<count)
		LOG_ERROR("dlist_pop_back()  pop count %lld > list->nobj %lld", count, list->nobj);

	q2=list->nobj/list->objpernode;
	r2=list->nobj%list->objpernode;
	list->nobj-=count;
	q1=list->nobj/list->objpernode;
	r1=list->nobj%list->objpernode;
	list->nobj+=count;
	
	if(q1==q2)
	{
		dlist_rdestroy(list, r2, r1);
		list->nobj-=count;
	}
	else
	{
		if(r2)
		{
			dlist_rdestroy(list, r2, 0);

			list->nobj-=r2;
			count-=r2;

			last=list->f;
			list->f=last->prev;
			free(last);
			--list->nnodes;
		}
		while(count>=list->objpernode)
		{
			dlist_rdestroy(list, list->objpernode, 0);

			list->nobj-=list->objpernode;
			count-=list->objpernode;

			last=list->f;
			list->f=last->prev;
			free(last);
			--list->nnodes;
		}
		if(count)
		{
			dlist_rdestroy(list, list->nobj%list->objpernode, r1);
			list->nobj-=count;
		}
	}
	if(!list->nnodes)//last object was popped out
		list->i=0;
}

void  dlist_first(DListHandle list, DListItHandle it)
{
	it->list=list;
	it->node=list->i;
	it->obj_idx=0;
}
void  dlist_last(DListHandle list, DListItHandle it)
{
	it->list=list;
	it->node=list->f;
	it->obj_idx=(list->nobj-1)%list->objpernode;
}
void* dlist_it_deref(DListItHandle it)
{
	if(it->obj_idx>=it->list->nobj)
	{
		LOG_ERROR("dlist_it_deref() out of bounds");
		return 0;
	}
	if(!it->node)
	{
		LOG_ERROR("dlist_it_deref() node == nullptr");
		return 0;
	}
	return it->node->data+it->obj_idx%it->list->objpernode*it->list->objsize;
}
int   dlist_it_inc(DListItHandle it)
{
	++it->obj_idx;
	if(it->obj_idx>=it->list->objpernode)
	{
		it->obj_idx=0;
		if(!it->node||!it->node->next)
			return 0;
			//LOG_ERROR("dlist_it_inc() attempting to read node == nullptr");
		it->node=it->node->next;
	}
	return 1;
}
int   dlist_it_dec(DListItHandle it)
{
	if(it->obj_idx)
		--it->obj_idx;
	else
	{
		if(!it->node||!it->node->prev)
			return 0;
			//LOG_ERROR("dlist_it_dec() attempting to read node == nullptr");
		it->node=it->node->prev;
		it->obj_idx=it->list->objpernode-1;
	}
	return 1;
}
#endif

//map/set (red-black tree)
#if 1
void map_init(MapHandle map, size_t esize, MapCmpFn comparator, void (*destructor)(void*))
{
	map->esize=esize;
	map->nnodes=0;
	map->root=0;
	map->comparator=comparator;
	map->destructor=destructor;
}
void map_clear_r(MapHandle map, RBNodeHandle node)
{
	if(node)
	{
		map_clear_r(map, node->left);
		map_clear_r(map, node->right);
		if(map->destructor)
			map->destructor(node->data);
		free(node);
	}
}
static RBNodeHandle* get_node_addr(MapHandle map, RBNodeHandle node)
{
	if(!node->parent)
		return &map->root;
	if(node==node->parent->left)
		return &node->parent->left;
	if(node==node->parent->right)
		return &node->parent->right;
	return 0;//inconsistent, should be unreachable
}
RBNodeHandle*	map_find(MapHandle map, const void *key)
{
	RBNodeHandle node;
	CmpRes result;

	for(node=map->root;node;)
	{
		result=map->comparator(key, node->data);
		switch(result)
		{
		case RESULT_LESS:
			node=node->left;
			break;
		case RESULT_GREATER:
			node=node->right;
			break;
		case RESULT_EQUAL://found: return node address
			return get_node_addr(map, node);
		}
	}
	return 0;//not found
}
static void rb_rotateleft(MapHandle map, RBNodeHandle x)
{
	RBNodeHandle y;

	y=x->right;
	x->right=y->left;
	if(y->left)
		y->left->parent=x;
	y->parent=x->parent;
	if(map->root==x)
		map->root=y;
	else if(x==x->parent->left)
		x->parent->left=y;
	else
		x->parent->right=y;
	y->left=x;
	x->parent=y;
}
static void rb_rotateright(MapHandle map, RBNodeHandle x)
{
	RBNodeHandle y;

	y=x->left;
	x->left=y->right;
	if(y->right)
		y->right->parent=x;
	y->parent=x->parent;
	if(map->root==x)
		map->root=y;
	else if(x==x->parent->right)
		x->parent->right=y;
	else
		x->parent->left=y;
	y->right=x;
	x->parent=y;
}
RBNodeHandle* map_insert(MapHandle map, const void *key, int *found)
{
	RBNodeHandle x, y, z;
	CmpRes result;

	//rb-insert - Cormen page 315
	x=map->root;
	y=0;
	while(x)
	{
		y=x;
		result=map->comparator(key, x->data);
		switch(result)
		{
		case RESULT_LESS:
			x=x->left;
			break;
		case RESULT_GREATER:
			x=x->right;
			break;
		case RESULT_EQUAL:
			if(found)
				*found=1;
			return get_node_addr(map, x);
		default:
			return 0;
		}
	}
	if(found)
		*found=0;

	z=(RBNodeHandle)malloc(sizeof(RBNodeHeader)+map->esize);
	if(!z)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	z->parent=y;
	z->left=z->right=0;
	z->is_red=1;
	memset(z->data, 0, map->esize);
	//memcpy(z->data, key, map->esize);//X  what if sizeof(key) != sizeof(data)
	x=z;

	if(!y)
		map->root=z;
	else
	{
		result=map->comparator(key, y->data);
		switch(result)
		{
		case RESULT_LESS:
			y->left=z;
			break;
		case RESULT_GREATER:
			y->right=z;
			break;
		case RESULT_EQUAL:
		default:
			free(z);
			return 0;
		}
	}

	//rb-fixup - Cormen page 316
	//https://github.com/gcc-mirror/gcc/blob/master/libstdc%2B%2B-v3/src/c%2B%2B98/tree.cc		line 195
	while(x!=map->root&&x->parent->is_red)//all nullptr's are black
	{
		RBNodeHandle xpp;

		xpp=x->parent->parent;
		if(xpp)
		{
			if(x->parent==xpp->left)
			{
				y=xpp->right;//uncle
				if(y&&y->is_red)//case 1
				{
					x->parent->is_red=0;
					y->is_red=0;
					xpp->is_red=1;
					x=xpp;
				}
				else//cases 2 & 3
				{
					if(x==x->parent->right)//case 2
					{
						x=x->parent;
						rb_rotateleft(map, x);
					}
					x->parent->is_red=0;
					xpp->is_red=1;
					rb_rotateright(map, xpp);
				}
			}
			else
			{
				y=xpp->left;//uncle
				if(y&&y->is_red)//case 1
				{
					x->parent->is_red=0;
					y->is_red=0;
					xpp->is_red=1;
					x=xpp;
				}
				else//cases 2 & 3
				{
					if(x==x->parent->left)//case 2
					{
						x=x->parent;
						rb_rotateright(map, x);
					}
					x->parent->is_red=0;
					xpp->is_red=1;
					rb_rotateleft(map, xpp);
				}
			}
		}
	}
	map->root->is_red=0;

	++map->nnodes;
	return get_node_addr(map, z);
}
//static void rb_transplant(MapHandle map, RBNodeHandle u, RBNodeHandle v)
//{
//	if(!u->parent)
//		map->root=v;
//	else if(u==u->parent->left)
//		u->parent->left=v;
//	else
//		u->parent->right=v;
//	v->parent=u->parent;
//}
static RBNodeHandle tree_minimum(RBNodeHandle root)
{
	if(!root)
		return 0;
	while(root->left)
		root=root->left;
	return root;
}
static RBNodeHandle tree_maximum(RBNodeHandle root)
{
	if(!root)
		return 0;
	while(root->right)
		root=root->right;
	return root;
}
int  map_erase(MapHandle map, const void *data, RBNodeHandle node)
{
	//https://github.com/gcc-mirror/gcc/blob/master/libstdc%2B%2B-v3/src/c%2B%2B98/tree.cc		line 286
	RBNodeHandle *leftmost, *rightmost, x, xp, y, z, *r2;
	size_t y_is_red;

	if(node)
		z=node;
	else if(data)
	{
		r2=map_find(map, data);
		if(!r2)
			return 0;
		z=*r2;
	}
	else
	{
		LOG_ERROR("map_erase() usage error: nullptr args");
		return 0;
	}

	leftmost=&map->root->left;
	rightmost=&map->root->right;
	y=z;
	x=0;
	xp=0;
	if(!y->left)		//z has at most one non-null child. y == z.
		x=y->right;		//x might be null
	else if(!y->right)	//z has exactly one non-null child. y == z.
		x=y->left;		//x is not null
	else
	{
		y=y->right;//z has two non-null children. Set y to z's successor.  x might be null
		while(y->left)
			y=y->left;
		x=y->right;
	}
	if(y!=z)
	{
		//relink y in place of z.  y is z's successor
		z->left->parent=y;
		y->left=z->left;
		if(y!=z->right)
		{
			xp=y->parent;
			if(x)
				x->parent=y->parent;
			y->parent->left=x;//y must be a child of left
			y->right=z->right;
			z->right->parent=y;
		}
		else
			xp=y;
		if(map->root==z)
			map->root=y;
		else if(z->parent->left==z)
			z->parent->left=y;
		else
			z->parent->right=y;
		y->parent=z->parent;
		y_is_red=y->is_red, y->is_red=z->is_red, z->is_red=y_is_red;
		y=z;
		//y now points to node to be actually deleted
	}
	else//y==z
	{
		xp=y->parent;

		if(x)
			x->parent=y->parent;

		if(map->root==z)
			map->root=x;
		else if(z->parent->left==z)
			z->parent->left=x;
		else
			z->parent->right=x;

		if(*leftmost==z)
		{
			if(!z->right)//z->left must be null also
				leftmost=&z->parent;//makes __leftmost == _M_header if __z == __root
			else
				leftmost=get_node_addr(map, tree_minimum(x));
		}
		if(*rightmost==z)
		{
			if(!z->left)//z->right must be null also
				rightmost=&z->parent;//makes __rightmost == _M_header if __z == __root
			else//x == z->left
				rightmost=get_node_addr(map, tree_maximum(x));
		}
	}
	if(!y->is_red)
	{
		RBNodeHandle wnode;

		while(x!=map->root&&(!x||!x->is_red))
		{
			if(x==xp->left)
			{
				wnode=xp->right;
				if(wnode->is_red)
				{
					wnode->is_red=0;
					xp->is_red=1;
					rb_rotateleft(map, xp);
					wnode=xp->right;
				}
				if((!wnode->left||!wnode->left->is_red)&&(!wnode->right||!wnode->right->is_red))
				{
					wnode->is_red=1;
					x=xp;
					xp=xp->parent;
				}
				else
				{
					if(!wnode->right||!wnode->right->is_red)
					{
						wnode->left->is_red=0;
						wnode->is_red=1;
						rb_rotateright(map, wnode);
						wnode=xp->right;
					}
					wnode->is_red=xp->is_red;
					xp->is_red=0;
					if(wnode->right)
						wnode->right->is_red=0;
					rb_rotateleft(map, xp);
					break;
				}
			}
			else//same as above with right <-> left
			{
				wnode=xp->left;
				if(wnode->is_red)
				{
					wnode->is_red=0;
					xp->is_red=1;
					rb_rotateright(map, xp);
					wnode=xp->left;
				}
				if((!wnode->right||!wnode->right->is_red)&&(!wnode->left||!wnode->left->is_red))
				{
					wnode->is_red=1;
					x=xp;
					xp=xp->parent;
				}
				else
				{
					if(!wnode->left||!wnode->left->is_red)
					{
						wnode->right->is_red=0;
						wnode->is_red=1;
						rb_rotateleft(map, wnode);
						wnode=xp->left;
					}
					wnode->is_red=xp->is_red;
					xp->is_red=0;
					if(wnode->left)
						wnode->left->is_red=0;
					rb_rotateright(map, xp);
					break;
				}
			}
		}
		if(x)
			x->is_red=0;
	}
	//return y;
	return 1;
}
void map_debugprint_r(RBNodeHandle *node, int depth, void (*printer)(RBNodeHandle *node, int depth))
{
	if(*node)
	{
		map_debugprint_r(&node[0]->left, depth+1, printer);
		printer(node, depth);
		map_debugprint_r(&node[0]->right, depth+1, printer);
	}
}
#endif


//single-linked list implementation
#if 1
void slist_init(SListHandle list, size_t esize, void (*destructor)(void*))
{
	list->esize=esize;
	list->count=0;
	list->destructor=destructor;
	list->front=list->back=0;
}
void slist_clear(SListHandle list)
{
	SNodeHandle temp;

	while(list->front)
	{
		temp=list->front;//copy pointer
		list->front=temp->prev;//advance

		if(list->destructor)//destroy & free
			list->destructor(temp->data);
		free(temp);
	}
	list->back=0;
	list->count=0;
}
static SNodeHandle slist_alloc_node(SListHandle list, SNodeHandle prev, const void *data)
{
	SNodeHandle temp;

	//allocate new node
	temp=(SNodeHandle)malloc(sizeof(SNode)+list->esize);
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	temp->prev=prev;
	if(data)
		memcpy(temp->data, data, list->esize);
	else
		memset(temp->data, 0, list->esize);
	return temp;
}
void* slist_push_front(SListHandle list, const void *data)
{
	SNodeHandle temp;

	//allocate new node
	temp=slist_alloc_node(list, list->front, data);

	if(list->count)//if front is not nullptr
		list->front=temp;//assign the new front
	else//list was empty
		list->front=list->back=temp;//initialize front and back with same node
	
	++list->count;
	return temp->data;
}
void* slist_push_back(SListHandle list, const void *data)
{
	SNodeHandle temp;

	//allocate new node
	temp=slist_alloc_node(list, 0, data);//back->prev is always nullptr

	if(list->count)//if back is not nullptr
	{
		list->back->prev=temp;//make old back point to new note (new back)
		list->back=temp;//assign the new back
	}
	else//list was empty
		list->front=list->back=temp;//initialize front and back with same node
	
	++list->count;
	return temp->data;
}
void *slist_front(SListHandle list)
{
	if(!list->front)
		return 0;
	return list->front->data;
}
void *slist_back(SListHandle list)
{
	if(!list->back)
		return 0;
	return list->back->data;
}
void slist_pop_front(SListHandle list)
{
	SNodeHandle temp;

	if(!list->count)
		return;
	temp=list->front;//copy front pointer
	list->front=temp->prev;//advance

	if(list->destructor)//destroy & free
		list->destructor(temp->data);
	free(temp);
	
	--list->count;
}
void slist_print(SListHandle list, void (*printer)(const void*))
{
	SNodeHandle node=list->front;
	while(node)
	{
		printer(node->data);
		node=node->prev;
	}
}
#endif


//implementation of bit-string
#if 1
BitstringHandle bitstring_construct(const void *src, size_t bitCount, size_t bitOffset, size_t bytePad)
{
	BitstringHandle str;
	size_t srcBytes, cap;

	srcBytes=(bitCount+7)>>3;
	cap=srcBytes+bytePad;
	str=(BitstringHandle)malloc(sizeof(BitstringHeader)+cap);
	if(!str)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	str->bitCount=bitCount;
	str->byteCap=cap;

	memset(str->data, 0, cap);
	if(src)
	{
		const unsigned char *srcbytes=(const unsigned char*)src;
		for(size_t b=0;b<bitCount;++b)
		{
			int bit=srcbytes[(bitOffset+b)>>3]>>((bitOffset+b)&7)&1;
			str->data[b>>3]|=bit<<(b&7);
		}
	}
	return str;
}
void bitstring_free(BitstringHandle *str)
{
	free(*str);
	*str=0;
}
void bitstring_append(BitstringHandle *str, const void *src, size_t bitCount, size_t bitOffset)
{
	size_t reqcap, newcap;
	size_t byteIdx;

	newcap=str[0]->byteCap;
	newcap+=!newcap;
	reqcap=(str[0]->bitCount+bitCount+7)/8;
	for(;newcap<reqcap;newcap<<=1);
	if(str[0]->byteCap<newcap)
	{
		void *p=realloc(*str, sizeof(BitstringHeader)+newcap);
		if(!p)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		*str=(BitstringHandle)p;
		str[0]->byteCap=newcap;
	}
	byteIdx=(str[0]->bitCount+7)/8;
	memset(str[0]->data+byteIdx, 0, newcap-byteIdx);
	if(src)
	{
		const unsigned char *srcbytes=(const unsigned char*)src;
		for(size_t b=0;b<bitCount;++b)
		{
			int bit=srcbytes[(bitOffset+b)>>3]>>((bitOffset+b)&7)&1;
			str[0]->data[(str[0]->bitCount+b)>>3]|=bit<<((str[0]->bitCount+b)&7);
		}
	}
	str[0]->bitCount+=bitCount;
}
int bitstring_get(BitstringHandle *str, size_t bitIdx)
{
	if(!*str)
	{
		LOG_ERROR("bitstring_get str=%p", *str);
		return 0;
	}
	if(bitIdx>=str[0]->bitCount)
	{
		LOG_ERROR("bitstring_get OOB: bitCount=%lld, bitIdx=%lld", (long long)str[0]->bitCount, (long long)bitIdx);
		return 0;
	}
	return str[0]->data[bitIdx>>3]>>(bitIdx&7)&1;
}
void bitstring_set(BitstringHandle *str, size_t bitIdx, int bit)
{
	if(!*str)
	{
		LOG_ERROR("bitstring_get str=%p", *str);
		return;
	}
	if(bitIdx>=str[0]->bitCount)
	{
		LOG_ERROR("bitstring_get OOB: bitCount=%lld, bitIdx=%lld", (long long)str[0]->bitCount, (long long)bitIdx);
		return;
	}
	if(bit)
		str[0]->data[bitIdx>>3]|=1<<(bitIdx&7);
	else
		str[0]->data[bitIdx>>3]&=0<<(bitIdx&7);
}
void bitstring_print(BitstringHandle str)
{
	for(int i=0;i<(int)str->bitCount;++i)
		printf("%d", bitstring_get(&str, i));
	//printf("\n");
}
#endif


//implementation of max-heap-based priority queue
#if 1
static void pqueue_realloc(PQueueHandle *pq, size_t count, size_t pad)//CANNOT be nullptr, array must be initialized with array_alloc()
{
	size_t size, newcap;

	size=(count+pad)*pq[0]->esize, newcap=pq[0]->esize;
	for(;newcap<size;newcap<<=1);
	if(newcap>pq[0]->byteCap)
	{
		//printf("realloc(%p, %lld)\n", *pq, sizeof(PQueueHeader)+newcap);//

		void *p2=realloc(*pq, sizeof(PQueueHeader)+newcap);
		if(!p2)
		{
			LOG_ERROR("Alloc error");
			return;
		}
		*pq=(PQueueHandle)p2;
		if(pq[0]->byteCap<newcap)
			memset(pq[0]->data+pq[0]->byteCap, 0, newcap-pq[0]->byteCap);
		pq[0]->byteCap=newcap;
	}
	pq[0]->count=count;
}
static void pqueue_heapifyup(PQueueHandle *pq, size_t idx, void *temp)
{
	for(;idx!=0;)
	{
		size_t parent=(idx-1)/2;
		if(pq[0]->less(pq[0]->data+parent*pq[0]->esize, pq[0]->data+idx*pq[0]->esize))
			memswap(pq[0]->data+parent*pq[0]->esize, pq[0]->data+idx*pq[0]->esize, pq[0]->esize, temp);
		else
			break;
		idx=parent;
	}
}
void        pqueue_heapifydown(PQueueHandle *pq, size_t idx, void *temp)
{
	size_t left, right, largest;

	for(;idx<pq[0]->count;)
	{
		left=(idx<<1)+1, right=left+1;

		if(left<pq[0]->count&&pq[0]->less(pq[0]->data+idx*pq[0]->esize, pq[0]->data+left*pq[0]->esize))//if [idx] < [left]
			largest=left;
		else
			largest=idx;

		if(right<pq[0]->count&&pq[0]->less(pq[0]->data+largest*pq[0]->esize, pq[0]->data+right*pq[0]->esize))//if [largest] < [right]
			largest=right;

		if(largest==idx)
			break;
		memswap(pq[0]->data+idx*pq[0]->esize, pq[0]->data+largest*pq[0]->esize, pq[0]->esize, temp);
		idx=largest;
	}
}
void        pqueue_buildheap(PQueueHandle *pq)
{
	void *temp;

	if(pq[0]->count<2)
		return;
	temp=malloc(pq[0]->esize);
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(ptrdiff_t i=pq[0]->count/2-1;i>=0;--i)
		pqueue_heapifydown(pq, i, temp);
	free(temp);
}

//Priority Queue
PQueueHandle pqueue_construct(
	size_t esize,
	size_t pad,
	int (*less)(const void*, const void*),
	void (*destructor)(void*)
)
{
	PQueueHandle pq;
	size_t cap;
	
	cap=esize+pad*esize;
	pq=(PQueueHandle)malloc(sizeof(PQueueHeader)+cap);
	if(!pq)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	pq->count=0;
	pq->esize=esize;
	pq->byteCap=cap;
	pq->destructor=destructor;
	pq->less=less;

	memset(pq->data, 0, cap);
	return pq;
}
void pqueue_free(PQueueHandle *pq)//can be nullptr
{
	if(*pq&&pq[0]->destructor)
	{
		for(size_t k=0;k<pq[0]->count;++k)
			pq[0]->destructor(pq[0]->data+k*pq[0]->esize);
	}
	free(*pq);
	*pq=0;
}

void  pqueue_enqueue(PQueueHandle *pq, const void *src)//src cannot be nullptr
{
	size_t start;
	void *temp;
	
	start=pq[0]->count*pq[0]->esize;
	pqueue_realloc(pq, pq[0]->count+1, 0);

	memcpy(pq[0]->data+start, src, pq[0]->esize);

	temp=malloc(pq[0]->esize);
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	pqueue_heapifyup(pq, pq[0]->count-1, temp);
	free(temp);
}
void* pqueue_front(PQueueHandle *pq)
{
	return pq[0]->data;
}
void  pqueue_dequeue(PQueueHandle *pq)
{
	void *temp;

	if(!pq[0]->count)
	{
		LOG_ERROR("pqueue_erase() out of bounds: size=%lld", (long long)pq[0]->count);
		return;
	}

	if(pq[0]->destructor)
		pq[0]->destructor(pq[0]->data);

	temp=malloc(pq[0]->esize);
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memswap(pq[0]->data, pq[0]->data+(pq[0]->count-1)*pq[0]->esize, pq[0]->esize, temp);
	--pq[0]->count;
	pqueue_heapifydown(pq, 0, temp);
	free(temp);
}
void  pqueue_print(PQueueHandle *pq, void (*printer)(const void*))
{
	for(int k=0;k<(int)pq[0]->count;++k)
		printer(pq[0]->data+k*pq[0]->esize);
}
void  pqueue_print_heap(PQueueHandle *pq, void (*printer)(const void*))
{
	for(int x=1, start=0;start<(int)pq[0]->count;x<<=1)
	{
		for(int i=0;i<x&&start+i<(int)pq[0]->count;++i)
			printer(pq[0]->data+(start+i)*pq[0]->esize);
		printf("\n");
		start+=x;
	}
	printf("\n");
}
#endif

ptrdiff_t get_filesize(const char *filename)//-1 not found,  0: not a file,  ...: regular file size
{
	struct stat info={0};
	int e2=stat(filename, &info);
	if(e2)
		return -1;
#if defined _MSC_VER || defined _WIN32
	if((info.st_mode&S_IFMT)==S_IFREG)
#else
	if(S_ISREG(info.st_mode))
#endif
		return info.st_size;
	return 0;
}

int acme_stricmp(const char *a, const char *b)//case insensitive strcmp
{
	if(!a||!b)
		return !a&&!b;
	while(*a&&tolower(*a)==tolower(*b))
		++a, ++b;
	return (*a>*b)-(*a<*b);
}
ptrdiff_t acme_strrchr(const char *str, ptrdiff_t len, char c)//find last occurrence, with known length for backward search
{
	ptrdiff_t k;

	for(k=len-1;k>=0;--k)
		if(str[k]==c)
			return k;
	return -1;
}
ArrayHandle filter_path(const char *path, int len)//replaces back slashes with slashes, adds trailing slash if missing, as ArrayHandle
{
	ArrayHandle path2;

	if(len<0)
		len=(int)strlen(path);
	STR_COPY(path2, path, len);
	for(ptrdiff_t k=0;k<(ptrdiff_t)path2->count;++k)//replace back slashes
	{
		if(path2->data[k]=='\\')
			path2->data[k]='/';
	}
	if(path2->data[path2->count-1]!='/')//ensure trailing slash
		STR_APPEND(path2, "/", 1, 1);
	return path2;
}
void get_filetitle(const char *fn, int len, int *idx_start, int *idx_end)//pass -1 for len if unknown
{
	int kpoint, kslash;
	if(len<0)
		len=(int)strlen(fn);
	for(kpoint=(int)len-1;kpoint>=0&&fn[kpoint]!='.';--kpoint);
	if(kpoint<0)
		kpoint=len;
	for(kslash=kpoint-1;kslash>=0&&fn[kslash]!='/'&&fn[kslash]!='\\';--kslash);
	++kslash;
	if(idx_start)
		*idx_start=kslash;
	if(idx_end)
		*idx_end=kpoint;
}
static const char* get_extension(const char *filename, ptrdiff_t len)//excludes the dot
{
	ptrdiff_t idx;

	idx=acme_strrchr(filename, len, '.');
	if(idx==-1)
		return 0;
	return filename+idx+1;
#if 0
	const char *dot=strrchr(filename, '.');//https://stackoverflow.com/questions/5309471/getting-file-extension-in-c
	if(!dot||dot==filename)
		return "";
	return dot+1;
#endif
}
static void free_str(void *p)
{
	ArrayHandle *str;
	
	str=(ArrayHandle*)p;
	array_free(str);
}
#ifdef __linux__
static int cmp_str(const void *p1, const void *p2)
{
	ArrayHandle const
		*s1=(ArrayHandle const*)p1,
		*s2=(ArrayHandle const*)p2;
	return strcmp_ci((char*)s1[0]->data, (char*)s2[0]->data);
}
#endif
ArrayHandle get_filenames(const char *path, const char **extensions, int extCount, int fullyqualified)
{
#if defined _MSC_VER || defined _WIN32
	ArrayHandle searchpath, filename, filenames;
	char c;
	WIN32_FIND_DATAA data={0};
	void *hSearch;
	int success;
	const char *extension;
	ptrdiff_t len;
	int found;
	
	//prepare searchpath
	searchpath=filter_path(path, -1);
	c='*';
	STR_APPEND(searchpath, &c, 1, 1);

	hSearch=FindFirstFileA((char*)searchpath->data, &data);//skip .
	if(hSearch==INVALID_HANDLE_VALUE)
		return 0;
	success=FindNextFileA(hSearch, &data);//skip ..

	STR_POPBACK(searchpath, 1);//pop the '*'
	ARRAY_ALLOC(ArrayHandle, filenames, 0, 0, 0, free_str);

	for(;;)
	{
		success=FindNextFileA(hSearch, &data);
		if(!success)
			break;
		len=strlen(data.cFileName);
		extension=get_extension(data.cFileName, len);
		if(!(data.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY))
		{
			found=0;
			for(int k=0;k<extCount;++k)
			{
				if(!acme_stricmp(extension, extensions[k]))
				{
					found=1;
					break;
				}
			}
			if(found)
			{
				STR_ALLOC(filename, 0);
				if(fullyqualified)
					STR_APPEND(filename, searchpath->data, searchpath->count, 1);
				STR_APPEND(filename, data.cFileName, len, 1);
				ARRAY_APPEND(filenames, &filename, 1, 1, 0);
			}
		}
	}
	success=FindClose(hSearch);
	array_free(&searchpath);
	return filenames;
#elif defined __linux__
	ArrayHandle searchpath, filename, filenames;
	searchpath=filter_path(path, -1);
	struct dirent *dir;
	DIR *d=opendir((char*)searchpath->data);
	if(!d)
		return 0;
	ARRAY_ALLOC(ArrayHandle, filenames, 0, 0, 0, free_str);
	while((dir=readdir(d)))
	{
		if(dir->d_type==DT_REG)//regular file
		{
			const char *name=dir->d_name;
			ptrdiff_t len=strlen(name);
			const char *extension=get_extension(name, len);
			int found=0;
			for(int k=0;k<extCount;++k)
			{
				if(!acme_stricmp(extension, extensions[k]))
				{
					found=1;
					break;
				}
			}
			if(found)
			{
				STR_ALLOC(filename, 0);
				if(fullyqualified)
					STR_APPEND(filename, searchpath->data, searchpath->count, 1);
				STR_APPEND(filename, name, len, 1);
				ARRAY_APPEND(filenames, &filename, 1, 1, 0);
			}
		}
	}
	closedir(d);
	array_free(&searchpath);
	qsort(filenames->data, filenames->count, filenames->esize, cmp_str);
	return filenames;
#endif
}

ArrayHandle load_file(const char *filename, int bin, int pad, int erroronfail)
{
	struct stat info={0};
	FILE *f;
	ArrayHandle str;
	char mode[3]={'r', bin?'b':'\0', 0};

	int e2=stat(filename, &info);
	if(e2)
	{
		if(erroronfail)
		{
#ifdef _MSC_VER
			strerror_s(g_buf, G_BUF_SIZE, errno);
			LOG_ERROR("Cannot open %s\n%s", filename, g_buf);
#else
			LOG_ERROR("Cannot open %s\n%s", filename, strerror(errno));
#endif
		}
		return 0;
	}
#ifdef _MSC_VER
	fopen_s(&f, filename, mode);
#else
	f=fopen(filename, mode);
#endif
	if(!f)
	{
		if(erroronfail)
		{
#ifdef _MSC_VER
			strerror_s(g_buf, G_BUF_SIZE, errno);
			LOG_ERROR("Cannot open %s\n%s", filename, g_buf);
#else
			LOG_ERROR("Cannot open %s\n%s", filename, strerror(errno));
#endif
		}
		return 0;
	}

	str=array_construct(0, 1, info.st_size, 1, pad+1, 0);
	str->count=fread(str->data, 1, info.st_size, f);
	fclose(f);
	memset(str->data+str->count, 0, str->cap-str->count);
	return str;
}
int save_file(const char *filename, const unsigned char *src, size_t srcsize, int is_bin)
{
	FILE *f;
	size_t bytesRead;
	char mode[]={'w', is_bin?'b':'\0', 0};

#ifdef _MSC_VER
	fopen_s(&f, filename, mode);
#else
	f=fopen(filename, mode);
#endif
	if(!f)
	{
		//printf("Failed to save %s\n", filename);
		return 0;
	}
	bytesRead=fwrite(src, 1, srcsize, f);
	fclose(f);
	if(is_bin&&bytesRead!=srcsize)
	{
		//printf("Failed to save %s\n", filename);
		return 0;
	}
	return 1;
}

ArrayHandle searchfor_file(const char *searchpath, const char *filetitle)
{
	ArrayHandle filename;
	ptrdiff_t size;

	STR_COPY(filename, filetitle, strlen(filetitle));
	size=get_filesize((char*)filename->data);
	if(size==-1)
	{
		array_insert(&filename, 0, searchpath, strlen(searchpath), 1, 1);
		size=get_filesize((char*)filename->data);
		if(size==-1)
			array_free(&filename);
	}
	return filename;
}

int get_cpu_features(void)//returns  0: old CPU,  1: AVX2,  3: AVX-512
{
#ifdef _MSC_VER
	int avx2=IsProcessorFeaturePresent(PF_AVX2_INSTRUCTIONS_AVAILABLE)!=0;
	int avx512=IsProcessorFeaturePresent(PF_AVX512F_INSTRUCTIONS_AVAILABLE)!=0;
#else
	__builtin_cpu_init();
	int avx2=__builtin_cpu_supports("avx2")!=0;
	int avx512=__builtin_cpu_supports("avx512f")!=0;
#endif
	return avx512<<1|avx2;
}
int query_cpu_cores(void)
{
#ifdef _WIN32
	SYSTEM_INFO info;
	GetNativeSystemInfo(&info);
	return info.dwNumberOfProcessors;
#else
	return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}
size_t query_mem_usage()
{
#ifdef _WIN32
	PROCESS_MEMORY_COUNTERS_EX pmc={0};
	K32GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
	return pmc.PrivateUsage;
#else
	FILE *file=fopen("/proc/self/status", "r");
	size_t result=0;
	char line[128];

	while(fgets(line, 128, file))
	{
		if(!strncmp(line, "VmSize:", 7))
		{
			//const char *p=line+7;
			//while((unsigned)(*p-'0')>=10)++p;
			//result=atoi(p);
			result=atoi(line+7);
			break;
		}
	}
	fclose(file);
	return result<<10;//KB -> bytes
#endif
}

typedef struct _ThreadParam
{
	void (*func)(void*);
	void *args;
	void *handle;
} ThreadParam;
static THREAD_RET THREAD_CALL thread_caller(void *param)
{
	ThreadParam *args=(ThreadParam*)param;

	args->func(args->args);
	
	return 0;
}
void* mt_exec(void (*func)(void*), void *args, int argbytes, int nthreads)
{
	ArrayHandle handles;

	ARRAY_ALLOC(ThreadParam, handles, 0, nthreads, 0, 0);
	if(!handles)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	for(int k=0;k<nthreads;++k)
	{
		ThreadParam *h;
		int error;
		
		h=(ThreadParam*)array_at(&handles, k);
		h->func=func;
		h->args=(char*)args+argbytes*k;
#ifdef _MSC_VER
		h->handle=(void*)_beginthreadex(0, 0, thread_caller, h, 0, 0);
		error=!h->handle;
#else
		error=pthread_create((pthread_t*)&h->handle, 0, thread_caller, h);
#endif
		if(error)
		{
			LOG_ERROR("Alloc error");
			return 0;
		}
	}
	return handles;
}
void mt_finish(void *ctx)
{
	ArrayHandle handles=(ArrayHandle)ctx;
#ifdef _MSC_VER
	ArrayHandle h2;
	ARRAY_ALLOC(HANDLE, h2, 0, handles->count, 0, 0);
	if(!h2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(int k=0;k<(int)handles->count;++k)
	{
		ThreadParam *src;
		HANDLE *dst;
		
		src=(ThreadParam*)array_at(&handles, k);
		dst=(HANDLE*)array_at(&h2, k);
		*dst=src->handle;
	}
	WaitForMultipleObjects((int)h2->count, (HANDLE*)h2->data, TRUE, INFINITE);
	for(int k=0;k<(int)handles->count;++k)
	{
		HANDLE *h=(HANDLE*)array_at(&h2, k);
		CloseHandle(*h);
	}
	array_free(&h2);
#else
	for(int k=0;k<(int)handles->count;++k)
	{
		ThreadParam *h=(ThreadParam*)array_at(&handles, k);
		pthread_join((pthread_t)h->handle, 0);
	}
#endif
	array_free(&handles);
}