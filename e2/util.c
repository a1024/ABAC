//util.c - Utilities implementation
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
#include<unistd.h>
#endif
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include<Windows.h>//QueryPerformance...
#include<conio.h>
#else
#define sprintf_s	snprintf
#define vsprintf_s	vsnprintf
#ifndef _HUGE
#define _HUGE	HUGE_VAL
#endif
#endif
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
	while(s1<s2)
	{
		memswap(s1, s2, esize, temp);
		s1+=esize, s2-=esize;
	}
	free(temp);
}
void memrotate(void *p, size_t byteoffset, size_t bytesize, void *temp)
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
	ptrdiff_t L=0, R=(ptrdiff_t)count-1, mid;
	int ret;

	while(L<=R)
	{
		mid=(L+R)>>1;
		ret=threeway(buf+mid*esize, val);
		if(ret<0)
			L=mid+1;
		else if(ret>0)
			R=mid-1;
		else
		{
			if(idx)
				*idx=mid;
			return 1;
		}
	}
	if(idx)
		*idx=L+(L<(ptrdiff_t)count&&threeway(buf+L*esize, val)<0);
	return 0;
}
void isort(void *base, size_t count, size_t esize, int (*threeway)(const void*, const void*))
{
	unsigned char *buf=(unsigned char*)base;
	size_t k;
	void *temp;

	if(count<2)
		return;

	temp=malloc((count>>1)*esize);
	for(k=1;k<count;++k)
	{
		size_t idx=0;
		binary_search(buf, k, esize, threeway, buf+k*esize, &idx);
		if(idx<k)
			memrotate(buf+idx*esize, (k-idx)*esize, (k+1-idx)*esize, temp);
	}
	free(temp);
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

int floor_log2(unsigned long long n)
{
#ifdef _MSC_VER
	unsigned long logn=0;
	int success=_BitScanReverse64(&logn, n);
	logn=success?logn:-1;
	return logn;
#elif defined __GNUC__
	int logn=63-__builtin_clzll(n);
	return logn;
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
#ifdef _MSC_VER
	unsigned long logn=0;
	int success=_BitScanReverse(&logn, n);
	logn=success?logn:-1;
	return logn;
#elif defined __GNUC__
	int logn=31-__builtin_clz(n);
	return logn;
#else
	int	logn=-!n;
	int	sh=(n>=1<<16)<<4;	logn+=sh, n>>=sh;
		sh=(n>=1<< 8)<<3;	logn+=sh, n>>=sh;
		sh=(n>=1<< 4)<<2;	logn+=sh, n>>=sh;
		sh=(n>=1<< 2)<<1;	logn+=sh, n>>=sh;
		sh= n>=1<< 1;		logn+=sh;
	return logn;
#endif
}
int ceil_log2(unsigned long long n)
{
	int lgn=floor_log2(n);
	lgn+=(1ULL<<lgn)<n;
	return lgn;
}
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
int minimum(int a, int b){return a<b?a:b;}
int maximum(int a, int b){return a>b?a:b;}
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

double time_sec()
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

	ti->secs=(float)(secs);
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
int		acme_strftime(char *buf, size_t len, const char *format)
{
	time_t tstamp;
	struct tm tformat;

	tstamp=time(0);
	localtime_s(&tformat, &tstamp);
	return (int)strftime(buf, len, format, &tformat);
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

//error handling
char first_error_msg[G_BUF_SIZE]={0}, latest_error_msg[G_BUF_SIZE]={0};
int log_error(const char *file, int line, int quit, const char *format, ...)
{
	int firsttime=first_error_msg[0]=='\0';

	ptrdiff_t size=strlen(file), start=size-1;
	for(;start>=0&&file[start]!='/'&&file[start]!='\\';--start);
	start+=start==-1||file[start]=='/'||file[start]=='\\';

	int printed=sprintf_s(latest_error_msg, G_BUF_SIZE, "\n%s(%d): ", file+start, line);
	va_list args;
	va_start(args, format);
	printed+=vsprintf_s(latest_error_msg+printed, G_BUF_SIZE-printed, format, args);
	va_end(args);

	if(firsttime)
		memcpy(first_error_msg, latest_error_msg, printed+1);
	fprintf(stderr, "%s\n", latest_error_msg);
	//messagebox(MBOX_OK, "Error", latest_error_msg);
	if(quit)
	{
		pause();
		exit(0);
	}
	return firsttime;
}
int valid(const void *p)//only makes sense with MSVC debugger
{
	size_t val=(size_t)p;

	if(sizeof(size_t)==4)
	{
		switch(val)
		{
		case 0:
		case 0xCCCCCCCC:
		case 0xFEEEFEEE:
		case 0xEEFEEEFE:
		case 0xCDCDCDCD:
		case 0xFDFDFDFD:
		case 0xBAADF00D:
		case 0xBAAD0000:
			return 0;
		}
	}
	else
	{
		if(val==0xCCCCCCCCCCCCCCCC)
			return 0;
		if(val==0xFEEEFEEEFEEEFEEE)
			return 0;
		if(val==0xEEFEEEFEEEFEEEFE)
			return 0;
		if(val==0xCDCDCDCDCDCDCDCD)
			return 0;
		if(val==0xBAADF00DBAADF00D)
			return 0;
		if(val==0xADF00DBAADF00DBA)
			return 0;
	}
	return 1;
}
int pause()
{
	int k;

	printf("Enter 0 to continue: ");
	while(!scanf("%d", &k));
	return k;
}
#ifdef _MSC_VER
int pause1()
{
	return _getch();
}
#endif
int pause_abort(const char *file, int lineno, const char *extraInfo)
{
	printf("INTERNAL ERROR %s(%d)\nABORTING\n", file, lineno);
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
	ArrayHandle p2;
	size_t size, newcap;

	ASSERT_P(*arr);
	size=(count+pad)*arr[0]->esize, newcap=arr[0]->esize;
	for(;newcap<size;newcap<<=1);
	if(newcap>arr[0]->cap)
	{
		p2=(ArrayHandle)realloc(*arr, sizeof(ArrayHeader)+newcap);
		ASSERT_P(p2);
		*arr=p2;
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
	{
		ASSERT_P(src);
		memfill(arr->data, src, dstsize, srcsize);
	}
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
	ArrayHandle p2;
	if(!*arr)
		return;
	arr[0]->cap=(arr[0]->count+pad)*arr[0]->esize;
	p2=(ArrayHandle)realloc(*arr, sizeof(ArrayHeader)+arr[0]->cap);
	ASSERT_P(p2);
	*arr=p2;
}

void* array_insert(ArrayHandle *arr, size_t idx, const void *data, size_t count, size_t rep, size_t pad)//cannot be nullptr
{
	size_t start, srcsize, dstsize, movesize;
	
	ASSERT_P(*arr);
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

	ASSERT_P(*arr);
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

	ASSERT_P(*arr);
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
	DST=(DNodeHandle)malloc(sizeof(DNodeHeader)+(PAYLOADSIZE)),\
	DST->prev=PREV,\
	DST->next=NEXT,\
	memcpy(DST->data, SRC->data, PAYLOADSIZE)
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
		PANIC();
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
static void dlist_fill_node(DListHandle list, size_t copysize, char **src, void *dst)
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
	char *buffer;
	void *ret;
	
	buffer=(char*)data;
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
		LOG_ERROR("dlist_it_deref() out of bounds");
	if(!it->node)
		LOG_ERROR("dlist_it_deref() node == nullptr");
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
#if 0
static void rb_transplant(MapHandle map, RBNodeHandle u, RBNodeHandle v)
{
	if(!u->parent)
		map->root=v;
	else if(u==u->parent->left)
		u->parent->left=v;
	else
		u->parent->right=v;
	v->parent=u->parent;
}
#endif
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
		RBNodeHandle w;

		while(x!=map->root&&(!x||!x->is_red))
		{
			if(x==xp->left)
			{
				w=xp->right;
				if(w->is_red)
				{
					w->is_red=0;
					xp->is_red=1;
					rb_rotateleft(map, xp);
					w=xp->right;
				}
				if((!w->left||!w->left->is_red)&&(!w->right||!w->right->is_red))
				{
					w->is_red=1;
					x=xp;
					xp=xp->parent;
				}
				else
				{
					if(!w->right||!w->right->is_red)
					{
						w->left->is_red=0;
						w->is_red=1;
						rb_rotateright(map, w);
						w=xp->right;
					}
					w->is_red=xp->is_red;
					xp->is_red=0;
					if(w->right)
						w->right->is_red=0;
					rb_rotateleft(map, xp);
					break;
				}
			}
			else//same as above with right <-> left
			{
				w=xp->left;
				if(w->is_red)
				{
					w->is_red=0;
					xp->is_red=1;
					rb_rotateright(map, xp);
					w=xp->left;
				}
				if((!w->right||!w->right->is_red)&&(!w->left||!w->left->is_red))
				{
					w->is_red=1;
					x=xp;
					xp=xp->parent;
				}
				else
				{
					if(!w->left||!w->left->is_red)
					{
						w->right->is_red=0;
						w->is_red=1;
						rb_rotateleft(map, w);
						w=xp->left;
					}
					w->is_red=xp->is_red;
					xp->is_red=0;
					if(w->left)
						w->left->is_red=0;
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
	unsigned char *srcbytes=(unsigned char*)src;

	srcBytes=(bitCount+7)>>3;
	cap=srcBytes+bytePad;
	str=(BitstringHandle)malloc(sizeof(BitstringHeader)+cap);
	str->bitCount=bitCount;
	str->byteCap=cap;

	memset(str->data, 0, cap);
	if(src)
	{
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
	void *p;
	size_t byteIdx;
	unsigned char *srcbytes=(unsigned char*)src;

	newcap=str[0]->byteCap;
	newcap+=!newcap;
	reqcap=(str[0]->bitCount+bitCount+7)/8;
	for(;newcap<reqcap;newcap<<=1);
	if(str[0]->byteCap<newcap)
	{
		p=realloc(*str, sizeof(BitstringHeader)+newcap);
		if(!p)
			LOG_ERROR("realloc returned nullptr");
		*str=p;
		str[0]->byteCap=newcap;
	}
	byteIdx=(str[0]->bitCount+7)/8;
	memset(str[0]->data+byteIdx, 0, newcap-byteIdx);
	if(src)
	{
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
	void *p2;
	size_t size, newcap;

	ASSERT_P(*pq);
	size=(count+pad)*pq[0]->esize, newcap=pq[0]->esize;
	for(;newcap<size;newcap<<=1);
	if(newcap>pq[0]->byteCap)
	{
		//printf("realloc(%p, %lld)\n", *pq, sizeof(PQueueHeader)+newcap);//

		p2=realloc(*pq, sizeof(PQueueHeader)+newcap);
		ASSERT_P(p2);
		*pq=(PQueueHandle)p2;
		if(pq[0]->byteCap<newcap)
			memset(pq[0]->data+pq[0]->byteCap, 0, newcap-pq[0]->byteCap);
		pq[0]->byteCap=newcap;
	}
	pq[0]->count=count;
}
void        pqueue_heapifyup(PQueueHandle *pq, size_t idx, void *temp)
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
	size_t L, R, largest;

	for(;idx<pq[0]->count;)
	{
		L=(idx<<1)+1, R=L+1;

		if(L<pq[0]->count&&pq[0]->less(pq[0]->data+idx*pq[0]->esize, pq[0]->data+L*pq[0]->esize))//if [idx] < [L]
			largest=L;
		else
			largest=idx;

		if(R<pq[0]->count&&pq[0]->less(pq[0]->data+largest*pq[0]->esize, pq[0]->data+R*pq[0]->esize))//if [largest] < [R]
			largest=R;

		if(largest==idx)
			break;
		memswap(pq[0]->data+idx*pq[0]->esize, pq[0]->data+largest*pq[0]->esize, pq[0]->esize, temp);
		idx=largest;
	}
}
void        pqueue_buildheap(PQueueHandle *pq)
{
	void *temp;

	temp=malloc(pq[0]->esize);
	for(size_t i=pq[0]->count/2-1;i>=0;--i)
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
	
	ASSERT_P(*pq);
	start=pq[0]->count*pq[0]->esize;
	pqueue_realloc(pq, pq[0]->count+1, 0);

	memcpy(pq[0]->data+start, src, pq[0]->esize);

	temp=malloc(pq[0]->esize);
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

	ASSERT_P(*pq);
	if(!pq[0]->count)
	{
		LOG_ERROR("pqueue_erase() out of bounds: size=%lld", (long long)pq[0]->count);
		return;
	}

	if(pq[0]->destructor)
		pq[0]->destructor(pq[0]->data);

	temp=malloc(pq[0]->esize);
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
	int error=stat(filename, &info);
	if(error)
		return -1;
	if((info.st_mode&S_IFMT)==S_IFREG)
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
ArrayHandle get_filetitle(const char *fn, int len)
{
	ArrayHandle title;
	int kpoint, kslash;
	if(len<0)
		len=(int)strlen(fn);
	for(kpoint=(int)len-1;kpoint>=0&&fn[kpoint]!='.';--kpoint);
	if(kpoint<0)
		kpoint=len;
	for(kslash=kpoint-1;kslash>=0&&fn[kslash]!='/'&&fn[kslash]!='\\';--kslash);
	++kslash;
	STR_COPY(title, fn+kslash, kpoint-kslash);
	return title;
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
void free_str(void *p)
{
	ArrayHandle *str;
	
	str=(ArrayHandle*)p;
	array_free(str);
}
ArrayHandle get_filenames(const char *path, const char **extensions, int extCount, int fullyqualified)
{
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

	for(;(success=FindNextFileA(hSearch, &data));)
	{
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
				STR_APPEND(filename, searchpath->data, searchpath->count, 1);
				STR_APPEND(filename, data.cFileName, strlen(data.cFileName), 1);
				ARRAY_APPEND(filenames, &filename, 1, 1, 0);
			}
		}
	}
	success=FindClose(hSearch);
	array_free(&searchpath);
	return filenames;
}

ArrayHandle load_file(const char *filename, int bin, int pad, int erroronfail)
{
	struct stat info={0};
	FILE *f;
	ArrayHandle str;
	char mode[3]={'r', bin?'b':0, 0};

	int error=stat(filename, &info);
	if(error)
	{
		if(erroronfail)
		{
			strerror_s(g_buf, G_BUF_SIZE, errno);
			LOG_ERROR("Cannot open %s\n%s", filename, g_buf);
		}
		return 0;
	}
	fopen_s(&f, filename, mode);
	//f=fopen(filename, "r");
	//f=fopen(filename, "r, ccs=UTF-8");//gets converted to UTF-16 on Windows
	if(!f)
	{
		if(erroronfail)
		{
			strerror_s(g_buf, G_BUF_SIZE, errno);
			LOG_ERROR("Cannot open %s\n%s", filename, g_buf);
		}
		return 0;
	}

	str=array_construct(0, 1, info.st_size, 1, pad+1, 0);
	str->count=fread(str->data, 1, info.st_size, f);
	fclose(f);
	memset(str->data+str->count, 0, str->cap-str->count);
	return str;
}
int save_file(const char *filename, const unsigned char *src, size_t srcSize, int is_bin)
{
	FILE *f;
	size_t bytesRead;
	char mode[]={'w', is_bin?'b':0, 0};

	fopen_s(&f, filename, mode);
	if(!f)
	{
		printf("Failed to save %s\n", filename);
		return 0;
	}
	bytesRead=fwrite(src, 1, srcSize, f);
	fclose(f);
	if(is_bin&&bytesRead!=srcSize)
	{
		printf("Failed to save %s\n", filename);
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

int query_cpu_cores()
{
#ifdef _WIN32
	SYSTEM_INFO info;
	GetNativeSystemInfo(&info);
	return info.dwNumberOfProcessors;
#else
	return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}
