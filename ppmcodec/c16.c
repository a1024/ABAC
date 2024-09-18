#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#endif
#include<sys/stat.h>
static const char file[]=__FILE__;


//	#define ENABLE_FILEGUARD	//makes using scripts harder

//	#define ESTIMATE_SIZE		//FOR DEBUGGING

//	#define AC_VALIDATE

#ifdef AC_VALIDATE
#define AC_IMPLEMENTATION
#include"entropy.h"
#endif
typedef unsigned short Emit_t;
static void update_CDF(const unsigned *hist, int nlevels, unsigned *CDF)
{
	const int hsum=hist[nlevels];
	for(int ks=0, c=0;ks<nlevels;++ks)
	{
		int freq=hist[ks];
		CDF[ks]=(int)(c*((1LL<<16)-nlevels)/hsum)+ks;
		c+=freq;
	}
	CDF[nlevels]=1<<16;
}
static void rescale_hist(unsigned *hist, int nlevels)
{
	int hsum=0;
	for(int ks=0;ks<nlevels;++ks)
		hsum+=hist[ks]=(hist[ks]+1)>>1;
	hist[nlevels]=hsum;
}
int c16_codec(const char *srcfn, const char *dstfn)
{
#ifdef ESTIMATE_SIZE
	double t=time_sec();
	double csizes[3]={0};
#endif
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
#ifdef ENABLE_FILEGUARD
	ptrdiff_t dstsize;
	dstsize=get_filesize(dstfn);
	if(dstsize>=0)
	{
		LOG_ERROR("Destination file already exists");
		return 1;
	}
#endif
	unsigned char *srcbuf=0;
	size_t srcsize=0;
	{
		struct stat info={0};
		int e2=stat(srcfn, &info);
		if(e2)
		{
			LOG_ERROR("Not found: \"%s\"", srcfn);
			return 1;
		}
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Not found: \"%s\"", srcfn);
			return 1;
		}
		srcsize=info.st_size;
		srcbuf=(unsigned char*)malloc(srcsize+1LL);
		size_t nread=fread(srcbuf, 1, srcsize, fsrc);
		if(nread!=srcsize)
			printf("Warning: file size %zd, read %zd\n", srcsize, nread);
		fclose(fsrc);
	}
	unsigned char *srcptr=srcbuf;
	int tag=0;
	memcpy(&tag, srcptr, 2); srcptr+=2;
	int fwd=tag==('P'|'6'<<8);
	if(!fwd&&tag!=('C'|'H'<<8))
	{
		LOG_ERROR("Unsupported source file");
		return 1;
	}
	int iw=0, ih=0;
	if(fwd)
	{
		int temp=0;
		if(*srcptr!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
		++srcptr;

		while((unsigned)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++-'0';
		if(iw<=0||*srcptr!=' ')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
		++srcptr;

		while((unsigned)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++-'0';
		if(ih<=0||*srcptr!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
		++srcptr;

		while((unsigned)(*srcptr-'0')<10)
			temp=10*temp+*srcptr++-'0';
		if(temp!=255||*srcptr!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
		++srcptr;
	}
	else
	{
		memcpy(&iw, srcptr, 4); srcptr+=4;
		memcpy(&ih, srcptr, 4); srcptr+=4;
	}
	int psize=(iw+16LL)*sizeof(short[4*4*2]);//4 padded rows * 4 channels max * {pixels, errors}
	short *pixels=(short*)_mm_malloc(psize, sizeof(__m128i));
	int hsize=(int)sizeof(int[3][257]);
	unsigned *hist=(unsigned*)malloc(hsize);
	unsigned *CDF=(unsigned*)malloc(hsize);
	if(!pixels||!hist||!CDF)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	FILLMEM(hist, 1, hsize, sizeof(int));
	hist[257*0+256]=256;
	hist[257*1+256]=256;
	hist[257*2+256]=256;
	update_CDF(hist+257*0, 256, CDF+257*0);
	update_CDF(hist+257*1, 256, CDF+257*1);
	update_CDF(hist+257*2, 256, CDF+257*2);
#ifdef AC_VALIDATE
	unsigned long long lo0, r0;
#endif
	if(fwd)
	{
		unsigned long long low=0, range=0xFFFFFFFFFFFF;
		unsigned short *dstbuf=(unsigned short*)malloc(sizeof(short[2])*iw*ih);//4/3 of image size
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		unsigned short *dstptr=dstbuf;
		memset(pixels, 0, psize);
		for(int ky=0, idx=0;ky<ih;++ky)//enc
		{
			ALIGN(32) short *rows[]=
			{
				pixels+((iw+16LL)*((ky-0LL)&1)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-1LL)&1)+8LL)*4*2,
			//	pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			//	pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
			};
			ALIGN(16) short preds[8]={0};
			for(int kx=0;kx<iw;++kx, idx+=3)
			{
				short
					*NW	=rows[1]-1*4*2,
					*N	=rows[1]+0*4*2,
					*W	=rows[0]-1*4*2,
					*curr	=rows[0]+0*4*2;
				__m128i mNW	=_mm_load_si128((__m128i*)NW);
				__m128i mN	=_mm_load_si128((__m128i*)N);
				__m128i mW	=_mm_load_si128((__m128i*)W);
				__m128i vmin=_mm_min_epi16(mN, mW);
				__m128i vmax=_mm_max_epi16(mN, mW);
				__m128i mp=_mm_sub_epi16(_mm_add_epi16(mN, mW), mNW);
				mp=_mm_max_epi16(mp, vmin);
				mp=_mm_min_epi16(mp, vmax);
				_mm_store_si128((__m128i*)preds, mp);

				//if(ky==3&&kx==114)//
				//	printf("");
				
				char comp[]=
				{
					srcptr[idx+1]-128,//y	g
					srcptr[idx+2]-128,//u	b
					srcptr[idx+0]-128,//v	r
				};

				curr[0]=comp[0];
				comp[0]-=preds[0];

				preds[1]+=curr[0];
				CLAMP2(preds[1], -128, 127);
				curr[1]=comp[1]-curr[0];
				comp[1]-=preds[1];
				
				preds[2]+=curr[0];
				CLAMP2(preds[2], -128, 127);
				curr[2]=comp[2]-curr[0];
				comp[2]-=preds[2];

				srcptr[idx+0]=comp[0]<<1^comp[0]>>7;//pack sign for fast decoding
				srcptr[idx+1]=comp[1]<<1^comp[1]>>7;
				srcptr[idx+2]=comp[2]<<1^comp[2]>>7;
				
				unsigned *curr_CDF;
				int sym[3], cdfs[3], freqs[3];

				sym[0]=srcptr[idx+0];
				sym[1]=srcptr[idx+1];
				sym[2]=srcptr[idx+2];

				curr_CDF=CDF+257*0;
				cdfs[0]=curr_CDF[sym[0]];
				freqs[0]=curr_CDF[sym[0]+1]-cdfs[0];

				curr_CDF=CDF+257*1;
				cdfs[1]=curr_CDF[sym[1]];
				freqs[1]=curr_CDF[sym[1]+1]-cdfs[1];
			
				curr_CDF=CDF+257*2;
				cdfs[2]=curr_CDF[sym[2]];
				freqs[2]=curr_CDF[sym[2]+1]-cdfs[2];
#ifdef ESTIMATE_SIZE
				csizes[0]-=log2((double)freqs[0]/0x10000);
				csizes[1]-=log2((double)freqs[1]/0x10000);
				csizes[2]-=log2((double)freqs[2]/0x10000);
#endif
				while(range<0x10000)
				{
					*dstptr++=(unsigned short)(low>>32);
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				low+=range*cdfs[0]>>16;
				range=(range*freqs[0]>>16)-1;
#ifdef AC_VALIDATE
				acval_enc(0, cdfs[0], freqs[0], lo0, lo0+r0, low, low+range, 0, 0);//
#endif

				while(range<0x10000)
				{
					*dstptr++=(unsigned short)(low>>32);
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				low+=range*cdfs[1]>>16;
				range=(range*freqs[1]>>16)-1;
#ifdef AC_VALIDATE
				acval_enc(0, cdfs[1], freqs[1], lo0, lo0+r0, low, low+range, 0, 0);//
#endif

				while(range<0x10000)
				{
					*dstptr++=(unsigned short)(low>>32);
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				low+=range*cdfs[2]>>16;
				range=(range*freqs[2]>>16)-1;
#ifdef AC_VALIDATE
				acval_enc(0, cdfs[2], freqs[2], lo0, lo0+r0, low, low+range, 0, 0);//
#endif
				++hist[257*0+sym[0]];
				++hist[257*0+256];
				++hist[257*1+sym[1]];
				++hist[257*1+256];
				++hist[257*2+sym[2]];
				++hist[257*2+256];

				rows[0]+=4*2;
				rows[1]+=4*2;
			//	rows[2]+=4*2;
			//	rows[3]+=4*2;
			}
			update_CDF(hist+257*0, 256, CDF+257*0);
			update_CDF(hist+257*1, 256, CDF+257*1);
			update_CDF(hist+257*2, 256, CDF+257*2);
			rescale_hist(hist+257*0, 256);
			rescale_hist(hist+257*1, 256);
			rescale_hist(hist+257*2, 256);
		}
		*dstptr++=(unsigned short)(low>>32);
		*dstptr++=(unsigned short)(low>>16);
		*dstptr++=(unsigned short)low;
		int streamsize=(int)(dstptr-dstbuf);//number of emits
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			LOG_ERROR("Cannot open \"%s\"", fdst);
			return 1;
		}
		fwrite("CH", 1, 2, fdst);
		fwrite(&iw, 1, 4, fdst);
		fwrite(&ih, 1, 4, fdst);
		fwrite(dstbuf, sizeof(short), streamsize, fdst);
		fclose(fdst);
		free(dstbuf);
	}
	else
	{
		int usize=(int)sizeof(char[3])*iw*ih;
		unsigned char *dstbuf=(unsigned char*)malloc(usize);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		memset(dstbuf, 0, usize);
		unsigned long long low=0, range=0xFFFFFFFFFFFF, code=0;
		unsigned short *sptr=(unsigned short*)srcptr;
		code=*sptr++;
		code=*sptr++|code<<16;
		code=*sptr++|code<<16;
		memset(pixels, 0, psize);
		for(int ky=0, idx=0;ky<ih;++ky)//dec
		{
			ALIGN(32) short *rows[]=
			{
				pixels+((iw+16LL)*((ky-0LL)&1)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-1LL)&1)+8LL)*4*2,
			//	pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			//	pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
			};
			ALIGN(16) short preds[8]={0};
			for(int kx=0;kx<iw;++kx, idx+=3)
			{
				unsigned *curr_CDF;
				unsigned long long code2;
				int cdf, freq, sym[3];
				
				curr_CDF=CDF+257*0;
				while(!(range>>16))
				{
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
				code2=(code-low)<<16|0xFFFF;
				sym[0]=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym[0]+2];
					if(range*freq>code2)
					{
						sym[0]+=range*curr_CDF[sym[0]+1]<=code2;
						break;
					}
					sym[0]+=2;
				}
				cdf=curr_CDF[sym[0]];
				freq=curr_CDF[sym[0]+1]-cdf;
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
#ifdef AC_VALIDATE
				acval_dec(0, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				dstbuf[idx+0]=sym[0];
				
				curr_CDF=CDF+257*1;
				while(!(range>>16))
				{
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
				code2=(code-low)<<16|0xFFFF;
				sym[1]=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym[1]+2];
					if(range*freq>code2)
					{
						sym[1]+=range*curr_CDF[sym[1]+1]<=code2;
						break;
					}
					sym[1]+=2;
				}
				cdf=curr_CDF[sym[1]];
				freq=curr_CDF[sym[1]+1]-cdf;
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
#ifdef AC_VALIDATE
				acval_dec(0, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				dstbuf[idx+1]=sym[1];
				
				curr_CDF=CDF+257*2;
				while(!(range>>16))
				{
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
				code2=(code-low)<<16|0xFFFF;
				sym[2]=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym[2]+2];
					if(range*freq>code2)
					{
						sym[2]+=range*curr_CDF[sym[2]+1]<=code2;
						break;
					}
					sym[2]+=2;
				}
				cdf=curr_CDF[sym[2]];
				freq=curr_CDF[sym[2]+1]-cdf;
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
#ifdef AC_VALIDATE
				acval_dec(0, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				dstbuf[idx+2]=sym[2];

				++hist[257*0+sym[0]];
				++hist[257*0+256];
				++hist[257*1+sym[1]];
				++hist[257*1+256];
				++hist[257*2+sym[2]];
				++hist[257*2+256];
				
				short
					*NW	=rows[1]-1*4*2,
					*N	=rows[1]+0*4*2,
					*W	=rows[0]-1*4*2,
					*curr	=rows[0]+0*4*2;
				__m128i mNW	=_mm_load_si128((__m128i*)NW);
				__m128i mN	=_mm_load_si128((__m128i*)N);
				__m128i mW	=_mm_load_si128((__m128i*)W);
				__m128i vmin=_mm_min_epi16(mN, mW);
				__m128i vmax=_mm_max_epi16(mN, mW);
				__m128i mp=_mm_sub_epi16(_mm_add_epi16(mN, mW), mNW);
				mp=_mm_max_epi16(mp, vmin);
				mp=_mm_min_epi16(mp, vmax);
				_mm_store_si128((__m128i*)preds, mp);

				char comp[]=
				{
					dstbuf[idx+0],
					dstbuf[idx+1],
					dstbuf[idx+2],
				};
				comp[0]=(unsigned char)comp[0]>>1^-(comp[0]&1);//unpack sign for fast decoding
				comp[1]=(unsigned char)comp[1]>>1^-(comp[1]&1);
				comp[2]=(unsigned char)comp[2]>>1^-(comp[2]&1);

				comp[0]+=preds[0];
				curr[0]=comp[0];

				preds[1]+=curr[0];
				CLAMP2(preds[1], -128, 127);
				comp[1]+=preds[1];
				curr[1]=comp[1]-curr[0];

				preds[2]+=curr[0];
				CLAMP2(preds[2], -128, 127);
				comp[2]+=preds[2];
				curr[2]=comp[2]-curr[0];

				dstbuf[idx+1]=comp[0]+128;//y	g
				dstbuf[idx+2]=comp[1]+128;//u	b
				dstbuf[idx+0]=comp[2]+128;//v	r

				rows[0]+=4*2;
				rows[1]+=4*2;
			//	rows[2]+=4*2;
			//	rows[3]+=4*2;
			}
			update_CDF(hist+257*0, 256, CDF+257*0);
			update_CDF(hist+257*1, 256, CDF+257*1);
			update_CDF(hist+257*2, 256, CDF+257*2);
			rescale_hist(hist+257*0, 256);
			rescale_hist(hist+257*1, 256);
			rescale_hist(hist+257*2, 256);
		}
		FILE *fdst=fopen(dstfn, "wb");
		if(!fdst)
		{
			LOG_ERROR("Cannot open \"%s\"", fdst);
			return 1;
		}
		fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		fwrite(dstbuf, sizeof(char), usize, fdst);
		fclose(fdst);
		free(dstbuf);
	}
	_mm_free(pixels);
	free(hist);
	free(CDF);
	return 0;
}