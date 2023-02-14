#include"battle.h"
#include<stdio.h>
#include<string.h>
#include<math.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
static const char file[]=__FILE__;

//	#define UABS_PRINT_RENORM
	#define UABS_PRINT_PROB
	#define RIGHTMOST_CODER

const int tag_uabs='U'|'B'<<8|'0'<<16|'1'<<24;
int uabs_encode_ch(const void *src, size_t nbytes, int bitoffset, int bitdepth, int bytestride, ArrayHandle *out)
{
	if(!src||!nbytes||!bitdepth||!bytestride||!out)
		return -2;//invalid func params

	size_t dststart;
	if(*out)
	{
		if(out[0]->esize!=1)
			return -1;//invalid out array
		dststart=out[0]->count;
	}
	else
		dststart=0;

	DList list;
	const unsigned char *buf=(const unsigned char*)src;
	dlist_init(&list, 1, 1024, 0);

	dlist_push_back(&list, &tag_uabs, 4);
	dlist_push_back(&list, 0, 8);

	unsigned short *prob=(unsigned short*)malloc(bitdepth*sizeof(short));

	size_t nsym=nbytes/bytestride;
	for(int kp=0;kp<bitdepth;++kp)
	{
		size_t freq=0;
		for(int ks=0;ks<nbytes;ks+=bytestride)
		{
			size_t bit_idx=(ks<<3)+bitoffset+kp, byte_idx=bit_idx>>3, sh=bit_idx&7;
			unsigned char bit=buf[byte_idx]>>sh&1;
			freq+=bit;
		}
		int qfreq=(int)((freq<<16)/nbytes);
		if(qfreq>0xFFFE)
			qfreq=0xFFFE;
		qfreq+=!qfreq;
		prob[kp]=qfreq;//adaptive prediction must be opposite to encoding direction
#ifdef UABS_PRINT_PROB
		printf("bit %d: %7d/%7d = %04X\n", kp, (int)freq, (int)nsym, (int)qfreq);
#endif
	}
	dlist_push_back(&list, prob, bitdepth*sizeof(short));

	unsigned state=0x10000;
	for(ptrdiff_t ks=nbytes-bytestride;ks>=0;ks-=bytestride)
	{
		for(int kp=bitdepth-1;kp>=0;--kp)
		{
			size_t bit_idx=(ks<<3)+bitoffset+kp, byte_idx=bit_idx>>3, sh=bit_idx&7;
			unsigned char bit=buf[byte_idx]>>sh&1;
			unsigned short p1=prob[kp];

			//update
			unsigned long long temp=(unsigned long long)state<<16;
#ifdef RIGHTMOST_CODER
			temp|=p1&-bit;
			p1^=-!bit;
			p1+=!bit;
			temp/=p1;//TODO Barett division table
			temp-=bit;
#else
			bit=!bit;
			temp|=p1&-bit;
			p1^=-bit;
			p1+=bit;
			temp/=p1;
			temp-=bit;
#endif

			//renorm
			//if(temp>=p1<<16)
			if(((unsigned*)&temp)[1])
			{
#ifdef UABS_PRINT_RENORM
				printf("ENC.RENORM ks %07lld kp%2d 0x%016llX ->16\n", ks, kp, temp);
#endif
				dlist_push_back(&list, &temp, 2);
				temp>>=16;
			}
			state=(unsigned)temp;
#if 0
#ifdef RIGHTMOST_CODER
			if(bit)
				state=(state<<16|p1)/p1-1;//TODO renorm
			else
				state=(state<<16)/(unsigned short)-p1;
#else
			if(bit)
				state=(state<<16)/p1;
			else
				state=(state<<16|p1)/(unsigned short)-p1-1;
#endif
#endif
		}
	}
	dlist_push_back(&list, &state, 4);
	free(prob);
	dlist_appendtoarray(&list, out);
	memcpy(out[0]->data+dststart+4, &list.nobj, 8);
	dlist_clear(&list);
	return 0;
}
const unsigned char* uabs_decode_ch(const unsigned char *data, size_t srclen, void *dst, size_t nbytes, int bitoffset, int bitdepth, int bytestride
#ifdef ENABLE_GUIDE
	, const unsigned char *guide
#endif
)
{
	if(!data||srclen<12+bitdepth*sizeof(short)||!dst||!nbytes)
	{
		LOG_ERROR("uABS: Invalid parameters");
		return 0;
	}

	int tag;
	memcpy(&tag, data, 4);
	if(tag!=tag_uabs)
	{
		LOG_ERROR("uABS: Invalid tag 0x%08X != 0x%08X", tag, tag_uabs);
		return 0;
	}

	size_t csize;
	memcpy(&csize, data+4, 8);

	if(srclen<csize)
	{
		LOG_ERROR("uABS: OOB %lld < %lld", srclen, 12+bitdepth*sizeof(short)+csize-4);
		return 0;
	}

	unsigned short *prob=(unsigned short*)malloc(bitdepth*sizeof(short));
	memcpy(prob, data+12, bitdepth*sizeof(short));

	unsigned char *buf=(unsigned char*)dst;
	unsigned state;
	const unsigned char
		*startptr=data+12+bitdepth*sizeof(short),
		*endptr=data+csize,
		*ptr=endptr;
	ptr-=4;
	if(ptr<startptr)
	{
		LOG_ERROR("OOB: start=%p, ptr=%p", startptr, ptr);
		return 0;
	}
	memcpy(&state, ptr, 4);
#ifdef UABS_PRINT_RENORM
	printf("DEC.START 0x%08X <-32\n", state);
#endif
	for(ptrdiff_t ks=0;ks<(ptrdiff_t)nbytes;ks+=bytestride)
	{
		for(int kp=0;kp<bitdepth;++kp)
		{
			unsigned short p1=prob[kp];
			unsigned long long temp1=(unsigned long long)state*p1;
#ifndef RIGHTMOST_CODER
			temp1+=0xFFFF;
#endif
			unsigned long long temp2=temp1+p1;
			temp1>>=16;
			temp2>>=16;
			unsigned char bit=(unsigned)temp2-(unsigned)temp1;
			
			size_t bit_idx=(ks<<3)+bitoffset+kp, byte_idx=bit_idx>>3, sh=bit_idx&7;
#ifdef ENABLE_GUIDE
			unsigned char b0=guide[byte_idx]>>sh&1;
			if(bit!=b0)
				LOG_ERROR("uABS: at %lld:%d dec %d instead of %d", byte_idx, sh, bit, b0);
#endif
			buf[byte_idx]=bit<<sh;

			bit=!bit;
			temp1^=-bit;
			temp1+=bit;
			state&=-bit;
			state+=(unsigned)temp1;
			//if(bit)
			//	state=(unsigned)temp1;
			//else
			//	state-=(unsigned)temp1;

			//renorm
			if(state<0x10000)
			{
				ptr-=2;
				if(ptr<startptr)
					LOG_ERROR("OOB: start=%p, ptr=%p", startptr, ptr);
				state<<=16;
				memcpy(&state, ptr, 2);
#ifdef UABS_PRINT_RENORM
				printf("DEC.RENORM ks%7lld kp%2d 0x%08X <-16\n", ks, kp, state);
#endif
			}
		}
	}

	free(prob);
	return endptr;
}
