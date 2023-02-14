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


//	#define AC_PRINT_HITCOUNT
//	#define AC_PRINT_PROB
//	#define AC_DISABLE_BYPASS


static const int tag_ac00='A'|'C'<<8|'0'<<16|'0'<<24;

#ifdef DEBUG_ABAC2
static int examined_plane=9, examined_start=0, examined_end=100;
#endif
long long ac0_encode(const void *src, ptrdiff_t nbytes, int bitoffset, int bitdepth, int bytestride, ArrayHandle *out, int loud)
{
	long long cycles=__rdtsc();

	if(!src||!nbytes||!bitdepth||!bytestride||!out)
		LOG_ERROR("binac0_encode(src=%p, nbytes=%d, bitdepth=%d, stride=%d, out=%p)", src, nbytes, bitdepth, bytestride, out);

	const unsigned char *buf=(const unsigned char*)src;
	unsigned char *prob=(unsigned char*)malloc(bitdepth);
	size_t nsym=nbytes/bytestride;
	for(int kp=0;kp<bitdepth;++kp)//analyze bitplanes
	{
		int startbitidx=bitoffset+kp;
		const unsigned char
			*srcptr=buf+(startbitidx>>3),
			*srcend=buf+nbytes;
		size_t freq=nsym;
		//ptrdiff_t step=bytestride<<3;//needed when bytestride < 1 byte
		for(ptrdiff_t sh=startbitidx&7;srcptr<srcend;srcptr+=bytestride)
		{
			unsigned char bit=*srcptr>>sh&1;
			freq-=bit;//freq=P(0)

			//sh=(sh+step)&7;//needed when bytestride < 1 byte
		}
		int qfreq=(int)(((freq<<8)+(nsym>>1))/nsym);//round(p0/size)
		qfreq+=(qfreq<0x80)-(qfreq>0x80);//0x80 means bypass

		//qfreq-=(freq<nsym)&(qfreq>255);
		//--qfreq;//add 1 when assigning p0
		//if(qfreq<1)//0x80 qfreq means bypass, which is triggered when ratio < 1
		//	qfreq=1;

		prob[kp]=qfreq;//adaptive prediction must be opposite to encoding direction
#ifdef AC_PRINT_PROB
		if(loud)
			printf("bit %d: P(0) = %7d/%7d = 0x%02X\n", kp, (int)freq, (int)nsym, qfreq);
#endif
	}

	size_t dststartidx;
	if(*out)
	{
		if(out[0]->esize!=1)
			return 0;
		ARRAY_APPEND(*out, 0, 0, 1, nbytes);
		dststartidx=out[0]->count;
	}
	else
	{
		ARRAY_ALLOC(char, *out, 0, 0, nbytes, 0);
		dststartidx=0;
	}
	ARRAY_APPEND(*out, 0, 12+bitdepth, 1, 0);//{tag, total written bytes, prob array}
	memcpy(out[0]->data+dststartidx, &tag_ac00, 4);

#ifdef AC_PRINT_HITCOUNT
	size_t hitcount=0;
#endif
	DList list={0};
	for(int kp=bitdepth-1;kp>=0;--kp)//bit-plane loop		encode MSB first
	{
#ifdef AC_PRINT_HITCOUNT
		size_t hitcount_p=0;
#endif
		unsigned char p0=prob[kp];
		int startidx=(bitoffset+kp)>>3, sh=(bitoffset+kp)&7;
#ifndef AC_DISABLE_BYPASS
		if(p0!=0x80)
#endif
		{
			unsigned r_start=0, r_end=0xFFFFFFFF;
			dlist_init(&list, 1, 1024, 0);
			for(ptrdiff_t ks=0;ks<nbytes;ks+=bytestride)
			{
				unsigned middle=r_start+(unsigned)((unsigned long long)(r_end-r_start)*p0>>8);

				int bit=buf[startidx+ks]>>sh&1;
				if(bit)
					r_start=middle+1;
				else
					r_end=middle-1;
#ifdef AC_PRINT_HITCOUNT
				int correct=bit^(p0>=0x80);
				hitcount_p+=correct;
#endif

				while((r_start^r_end)<0x1000000)
				{
					dlist_push_back1(&list, (char*)&r_start+3);
					r_start<<=8;
					r_end=r_end<<8|0xFF;
				}

				if(r_start+3<r_start||r_start+3>r_end)
				{
					dlist_push_back1(&list, (char*)&r_start+3);//big endian
					dlist_push_back1(&list, (char*)&r_start+2);
					dlist_push_back1(&list, (char*)&r_start+1);
					dlist_push_back1(&list, (char*)&r_start);

					r_start=0, r_end=0xFFFFFFFF;//because 1=0.9999...
				}
			}
			dlist_push_back1(&list, (char*)&r_start+3);//big endian
			dlist_push_back1(&list, (char*)&r_start+2);
			dlist_push_back1(&list, (char*)&r_start+1);
			dlist_push_back1(&list, (char*)&r_start);
		}
		size_t csize;
#ifndef AC_DISABLE_BYPASS
		if(p0==0x80||list.nobj>nsym)//ratio > 1, bypass
		{
			prob[kp]=0x80;
			size_t bypass_idx=out[0]->count, bypass_len=(nsym+7)>>3;
			ARRAY_APPEND(*out, 0, bypass_len, 1, 0);
			unsigned char frag;
			ptrdiff_t ks=0, end=nbytes-bytestride*7, kb=0;
			for(;ks<end;++kb)
			{
				frag=0;
				frag|= buf[startidx+ks]>>sh&1,     ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<1, ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<2, ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<3, ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<4, ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<5, ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<6, ks+=bytestride;
				frag|=(buf[startidx+ks]>>sh&1)<<7, ks+=bytestride;
				out[0]->data[bypass_idx+kb]=frag;
			}
			if(ks<nbytes)
			{
				frag=0;
				int sh2=0;
				do
				{
					frag|=(buf[startidx+ks]>>sh&1)<<sh2;
					ks+=bytestride;
					++sh2;
				}
				while(ks<nbytes);
				out[0]->data[bypass_idx+kb]=frag;
			}
			csize=nsym;
		}
		else
#endif
		{
			csize=list.nobj;
			dlist_appendtoarray(&list, out);
		}
		if(p0!=0x80)
			dlist_clear(&list);

#ifdef AC_PRINT_HITCOUNT
		if(loud)
			printf("bit %d: r =%6d /%6d = %lf, hit=%6d=%lf%%\n", kp, (int)nsym, (int)(csize<<3), (double)nsym/(csize<<3), (int)hitcount_p, 100.*hitcount_p/nsym);
		hitcount+=hitcount_p;
#else
		if(loud)
			printf("bit %d: r =%6d /%6d = %lf\n", kp, nsym, csize, (double)nsym/csize);
#endif
	}
	size_t byteswritten=out[0]->count-dststartidx;
	memcpy(out[0]->data+dststartidx+4, &byteswritten, 8);
	memcpy(out[0]->data+dststartidx+12, prob, bitdepth);
	cycles=__rdtsc()-cycles;

	if(loud)
	{
		size_t original_bitsize=nsym*bitdepth, compressed_bitsize=byteswritten<<3;
		printf("AC encode:  %lld cycles, %lf c/byte\n", cycles, (double)(cycles<<3)/original_bitsize);
		printf("Size: %d -> %d, ratio: %lf, %lf bpp\n", (int)original_bitsize>>3, (int)compressed_bitsize>>3, (double)original_bitsize/compressed_bitsize, (double)compressed_bitsize/nsym);
#ifdef AC_PRINT_HITCOUNT
		printf("Predicted: %6lld / %6lld = %lf%%\n", hitcount, original_bitsize, 100.*hitcount/original_bitsize);
#endif
		
		printf("Preview:\n");
		int kprint=200;
		if(byteswritten<kprint)
			kprint=(int)byteswritten;
		for(int k=0;k<kprint;++k)
			printf("%02X-", out[0]->data[dststartidx+k]&0xFF);
		printf("\n");
	}
	return byteswritten;
}
const void* ac0_decode(const void *srcstart, const void *srcend, void *dst, size_t nbytes, int bitoffset, int bitdepth, int bytestride, int loud
#ifdef ENABLE_GUIDE
	, unsigned char *guide
#endif
)//set the dst buffer to zero
{
	if(!srcstart||!srcend||!dst||!nbytes||!bitdepth||!bytestride)
		LOG_ERROR("abac4_decode(in_start=%p, in_end=%p, dst=%p, imsize=%d, bitdepth=%d, stride=%d)", srcstart, srcend, dst, nbytes, bitdepth, bytestride);

	long long cycles=__rdtsc();

	const unsigned char
		*data=(const unsigned char*)srcstart,
		*srcptr=data+12+bitdepth,
		*dataend=(const unsigned char*)srcend;
	size_t datalen=dataend-data;
	unsigned char *buf=(unsigned char*)dst;

	if(datalen<12)
		LOG_ERROR("AC: datalen = %p", datalen);
	int tag;
	memcpy(&tag, data, 4);
	if(tag!=tag_ac00)
		LOG_ERROR("AC: invalid tag 0x%08X, expected 0x%08X", tag, tag_ac00);

	size_t csize;
	memcpy(&csize, data+4, 8);
	if(datalen<csize)
		LOG_ERROR("AC: Unexpected EOF, datalen = %p", datalen);
	dataend=data+csize;

	unsigned char *prob=(unsigned char*)malloc(bitdepth);
	memcpy(prob, data+12, bitdepth);

	size_t nsym=nbytes/bytestride;
	for(int kp=bitdepth-1;kp>=0;--kp)
	{
		int startidx=(bitoffset+kp)>>3, sh=(bitoffset+kp)&7;
		unsigned char p0=prob[kp];
#ifndef AC_DISABLE_BYPASS
		if(p0==0x80)//bypass
		{
			ptrdiff_t ks=0, end=nbytes-bytestride*7;
			for(;ks<end;++srcptr)
			{
				unsigned char frag=*srcptr;
				buf[startidx+ks]|=(frag>>0&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>1&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>2&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>3&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>4&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>5&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>6&1)<<sh, ks+=bytestride;
				buf[startidx+ks]|=(frag>>7&1)<<sh, ks+=bytestride;
			}
			if(ks<nbytes)
			{
				unsigned char frag=*srcptr;
				int sh2=0;
				do
				{
					buf[startidx+ks]|=(frag>>sh2&1)<<sh;
					ks+=bytestride;
					++sh2;
				}while(ks<nbytes);
				++srcptr;
			}
#ifdef ENABLE_GUIDE
			for(ptrdiff_t ks=0;ks<nbytes;ks+=bytestride)
			{
				unsigned char
					b1=guide[startidx+ks]>>sh&1,
					b2=buf[startidx+ks]>>sh&1;
				if(b2!=b1)
					LOG_ERROR("AC bypass error kp%2d ks%7d dec %d != %d", kp, ks, b2, b1);
			}
#endif
		}
		else
#endif
		{
			unsigned r_start=0, r_end=0xFFFFFFFF, code;

#ifdef ENABLE_GUIDE
			unsigned r_start2=0, r_end2=0xFFFFFFFF;
			size_t enc_renorms=0, dec_renorms=0;
#endif
			srcptr+=4;
			if(srcptr>dataend)
				LOG_ERROR("AC OOB [start] srcptr %p dataend %p", srcptr, dataend);
			code=srcptr[-4]<<24|srcptr[-3]<<16|srcptr[-2]<<8|srcptr[-1];//big endian

			//for(ptrdiff_t ks=0;;)
			for(ptrdiff_t ks=0;ks<(ptrdiff_t)nbytes;ks+=bytestride)
			{
				unsigned middle=r_start+(unsigned)((unsigned long long)(r_end-r_start)*p0>>8);

				unsigned char bit=code>middle;
#ifdef ENABLE_GUIDE
				unsigned mid2=r_start2+(unsigned)((unsigned long long)(r_end2-r_start2)*p0>>8);
				unsigned char b0=guide[startidx+ks]>>sh&1;
				if(bit!=b0||enc_renorms!=dec_renorms||r_start!=r_start2||r_end!=r_end2)
					LOG_ERROR("kp%2d ks%7d dec %d != %d", kp, ks, bit, b0);
				if(b0)
					r_start2=mid2+1;
				else
					r_end2=mid2-1;
				while((r_start2^r_end2)<0x1000000)
				{
					r_start2<<=8;
					r_end2=r_end2<<8|0xFF;
					++enc_renorms;
				}
				if(r_start2+3<r_start2||r_start2+3>r_end2)
				{
					r_start2=0, r_end2=0xFFFFFFFF;//because 1=0.9999...
					enc_renorms+=4;
				}
#endif

				buf[startidx+ks]|=bit<<sh;

				//ks+=bytestride;
				//if(ks>=nbytes)
				//	break;

				if(bit)
					r_start=middle+1;
				else
					r_end=middle-1;

				//if(srcptr<dataend)
				//{
					while((r_start^r_end)<0x1000000)
					{
						++srcptr;
						if(srcptr>dataend)
							LOG_ERROR("AC OOB kp%2d ks%7d srcptr %p dataend %p", kp, ks, srcptr, dataend);

						//++srcptr;
						code=code<<8|srcptr[-1];
						r_start<<=8;
						r_end=r_end<<8|0xFF;
#ifdef ENABLE_GUIDE
						++dec_renorms;
#endif
					}
				//}
				//else
				//	srcptr=dataend;

				if(r_start+3<r_start||r_start+3>r_end)
				{
					srcptr+=4;
					if(srcptr>dataend)
						LOG_ERROR("AC OOB kp%2d ks%7d srcptr %p dataend %p", kp, ks, srcptr, dataend);

					code=srcptr[-4]<<24|srcptr[-3]<<16|srcptr[-2]<<8|srcptr[-1];//big endian

					r_start=0, r_end=0xFFFFFFFF;//because 1=0.9999...
#ifdef ENABLE_GUIDE
					dec_renorms+=4;
#endif
				}
			}
		}
	}
	free(prob);
	cycles=__rdtsc()-cycles;

	if(loud)
		printf("AC decode:  %lld cycles, %lf c/byte\n", cycles, (double)(cycles<<3)/(nsym*bitdepth));
	return srcptr;
}