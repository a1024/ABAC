#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#endif
//#include"blist.h"
static const char file[]=__FILE__;


	#define BYPASS_ON_INFLATION
//	#define ENABLE_FILEGUARD	//makes using scripts harder

//	#define ESTIMATE_SIZE
	#define USE_RCTJ2K
	#define USE_AC48

//	#include"entropy.h"
#ifdef USE_AC48
typedef unsigned short Emit_t;
#else
typedef unsigned Emit_t;
#endif
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
int c10_codec(const char *srcfn, const char *dstfn)
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
	ptrdiff_t dstsize;
#ifdef ENABLE_FILEGUARD
	dstsize=get_filesize(dstfn);
	if(dstsize>=0)
	{
		LOG_ERROR("Destination file already exists");
		return 1;
	}
#endif
	FILE *fsrc=fopen(srcfn, "rb");
	int tag=0;
	size_t nread=fread(&tag, 1, 2, fsrc), nwritten=0;
	if(nread!=2)
	{
		LOG_ERROR("File is empty");
		return 1;
	}
	int fwd=tag==('P'|'6'<<8);
	int bypass=tag==('B'|'P'<<8);
	if(!fwd&&tag!=('C'|'H'<<8)&&!bypass)
	{
		LOG_ERROR("Unsupported source file");
		return 1;
	}
	int iw=0, ih=0;
	int usize=0;
	if(fwd)
	{
		int temp=fgetc(fsrc);
		if(temp!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
		nread=fscanf(fsrc, "%d %d", &iw, &ih);
		if(nread!=2)
		{
			LOG_ERROR("Unsupported PPM file");
			return 1;
		}
		int vmax=0;
		nread=fscanf(fsrc, "%d", &vmax);
		if(nread!=1||vmax!=255)
		{
			LOG_ERROR("Unsupported PPM file");
			return 1;
		}
		temp=fgetc(fsrc);
		if(temp!='\n')
		{
			LOG_ERROR("Invalid PPM file");
			return 1;
		}
	}
	else
	{
		nread=fread(&iw, 1, 4, fsrc);
		nread+=fread(&ih, 1, 4, fsrc);
		if(nread!=4LL+4)
		{
			LOG_ERROR("Unsupported archive");
			return 1;
		}
	}
	usize=iw*ih*3;
	unsigned char *buffer=(unsigned char*)malloc(usize);
	if(!buffer)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	FILE *fdst=fopen(dstfn, "wb");
	nwritten=0;
	if(fwd)
	{
		nread=fread(buffer, 1, usize, fsrc);
		if(nread!=usize)
		{
			LOG_ERROR("Corrupt PPM");
			return 1;
		}

		nwritten+=fwrite("CH", 1, 2, fdst);
		nwritten+=fwrite(&iw, 1, 4, fdst);
		nwritten+=fwrite(&ih, 1, 4, fdst);
	}
	else
	{
		nwritten+=fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
		if(bypass)
		{
			printf("Bypass decode\n");
			nread=fread(buffer, 1, usize, fsrc);
			if(nread!=usize)
				printf("Error  read %d/%d bytes\n", (int)nread, usize);
			fwrite(buffer, 1, usize, fdst);
			fclose(fsrc);
			fclose(fdst);
			free(buffer);
			return 0;
		}
		memset(buffer, 0, usize);
	}
	int psize=(iw+16LL)*sizeof(short[4*4*2]);//4 padded rows * 4 channels max * {pixels, errors}
	short *pixels=(short*)_mm_malloc(psize, sizeof(__m128i));
	int hsize=(int)sizeof(int[257+513+513]);
	unsigned *hist=(unsigned*)malloc(hsize);
	unsigned *CDF=(unsigned*)malloc(hsize);
	int sbufsize=iw*16;//rowsize = iw*3, allocate iw*16 just in case
	Emit_t *sbuf=(Emit_t*)malloc(sbufsize);//stream buffer
	if(!pixels||!hist||!CDF||!sbuf)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(sbuf, 0, sbufsize);
	memset(hist, 0, hsize);
	for(int k=0;k<257;++k)
		CDF[k]=k<<8;
	for(int ks=0;ks<513;++ks)
	{
		CDF[257+ks]=ks<<7;
		CDF[257+513+ks]=ks<<7;
	}
	memset(pixels, 0, psize);
	int rowsize=iw*3;
	for(int ky=0, idx=0;ky<ih;++ky)
	{
		ALIGN(32) short *rows[]=
		{
			pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		ALIGN(16) short preds[8]={0};
		unsigned *curr_CDF=0;
		unsigned cdfs[3]={0}, freqs[3]={0}, c2;
		int yuv[3]={0}, errors[3]={0};
		unsigned long long low=0, range=0xFFFFFFFFFFFF, code=0;
#ifndef USE_AC48
		unsigned long long r2;
#endif
#ifdef AC_VALIDATE
		unsigned long long lo0, r0;
#endif
		Emit_t *sptr=sbuf;
		int streamsize=0;
		if(!fwd)
		{
			nread=fread(&streamsize, 1, 4, fsrc);
			if(streamsize)
			{
				nread=fread(sbuf, sizeof(Emit_t), streamsize, fsrc);
#ifdef USE_AC48
				code=*sptr++;
				code=*sptr++|code<<16;
				code=*sptr++|code<<16;
#else
				code=*sptr++;
				code=*sptr++|code<<32;
#endif
			}
			else
				nread=fread(buffer+rowsize*ky, 1, rowsize, fsrc);
		}
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
			
			//if(ky==1&&kx==39)//
			//if(ky==1&&kx==5)//
			//if(ky==1&&kx==743)//
			//	printf("");
			if(fwd||!streamsize)
			{
				yuv[0]=buffer[idx+1]-128;//g
				yuv[1]=buffer[idx+2]-128;//b
				yuv[2]=buffer[idx+0]-128;//r
#ifdef USE_RCTJ2K
				yuv[1]-=yuv[0];
				yuv[2]-=yuv[0];
				yuv[0]+=(yuv[1]+yuv[2])>>2;
#else
				//Pei09 RCT		b-=(87*r+169*g+128)>>8; r-=g; g+=(86*r+29*b+128)>>8;
				yuv[1]-=(87*yuv[2]+169*yuv[0]+128)>>8;
				yuv[2]-=yuv[0];
				yuv[0]+=(86*yuv[2]+29*yuv[1]+128)>>8;
#endif
				curr[0]=yuv[0];
				curr[1]=yuv[1];
				curr[2]=yuv[2];

				errors[0]=yuv[0]-preds[0];
				errors[1]=yuv[1]-preds[1];
				errors[2]=yuv[2]-preds[2];
				curr[4]=errors[0];
				curr[5]=errors[1];
				curr[6]=errors[2];
				errors[0]=errors[0]<<(32-8)>>(32-8);
				errors[1]=errors[1]<<(32-9)>>(32-9);
				errors[2]=errors[2]<<(32-9)>>(32-9);
				errors[0]=errors[0]<<1^errors[0]>>31;//pack sign
				errors[1]=errors[1]<<1^errors[1]>>31;
				errors[2]=errors[2]<<1^errors[2]>>31;

				if(fwd)
				{
					curr_CDF=CDF;
					cdfs[0]=curr_CDF[errors[0]];
					freqs[0]=curr_CDF[errors[0]+1]-cdfs[0];

					curr_CDF=CDF+257;
					cdfs[1]=curr_CDF[errors[1]];
					freqs[1]=curr_CDF[errors[1]+1]-cdfs[1];
			
					curr_CDF=CDF+257+513;
					cdfs[2]=curr_CDF[errors[2]];
					freqs[2]=curr_CDF[errors[2]+1]-cdfs[2];
#ifdef ESTIMATE_SIZE
					csizes[0]-=log2((double)freqs[0]/0x10000);
					csizes[1]-=log2((double)freqs[1]/0x10000);
					csizes[2]-=log2((double)freqs[2]/0x10000);
#endif
#ifdef USE_AC48
					while(range<0x10000)
					{
						*sptr++=(Emit_t)(low>>32);
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
						*sptr++=(Emit_t)(low>>32);
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
						*sptr++=(Emit_t)(low>>32);
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
#else
					r2=range>>16;
					if(!r2)
					{
						*sptr++=(unsigned)(low>>32);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
						r2=range>>16;
					}
					low+=r2*cdfs[0]+((range&0xFFFF)*cdfs[0]>>16);
					range=r2*freqs[0]+((range&0xFFFF)*freqs[0]>>16)-1;
				
					r2=range>>16;
					if(!r2)
					{
						*sptr++=(unsigned)(low>>32);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
						r2=range>>16;
					}
					low+=r2*cdfs[1]+((range&0xFFFF)*cdfs[1]>>16);
					range=r2*freqs[1]+((range&0xFFFF)*freqs[1]>>16)-1;

					r2=range>>16;
					if(!r2)
					{
						*sptr++=(unsigned)(low>>32);
						low<<=32;
						range=range<<32|0xFFFFFFFF;
						if(range>~low)
							range=~low;
						r2=range>>16;
					}
					low+=r2*cdfs[2]+((range&0xFFFF)*cdfs[2]>>16);
					range=r2*freqs[2]+((range&0xFFFF)*freqs[2]>>16)-1;
#endif
				}
			}
			else
			{
				unsigned cdf, freq, sym;
#ifdef USE_AC48
				//if(!ky&&kx==513)//
				//	printf("");
				while(!(range>>16))
				{
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
				c2=(unsigned)(((code-low)<<16|0xFFFF)/range);
				curr_CDF=CDF;
				sym=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym+1];
					if((unsigned)freq>c2)
						break;
					++sym;
				}
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				freq-=cdf;
				errors[0]=sym;
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
#ifdef AC_VALIDATE
				acval_dec(0, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				
				while(!(range>>16))
				{
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
				c2=(unsigned)(((code-low)<<16|0xFFFF)/range);
				curr_CDF=CDF+257;
				sym=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym+1];
					if((unsigned)freq>c2)
						break;
					++sym;
				}
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				freq-=cdf;
				errors[1]=sym;
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
#ifdef AC_VALIDATE
				acval_dec(0, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				
				while(!(range>>16))
				{
					low=low<<16&0xFFFFFFFFFFFF;
					range=range<<16|0xFFFF;
					code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
					if(range>(low^0xFFFFFFFFFFFF))
						range=low^0xFFFFFFFFFFFF;
				}
				c2=(unsigned)(((code-low)<<16|0xFFFF)/range);
				curr_CDF=CDF+257+513;
				sym=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym+1];
					if((unsigned)freq>c2)
						break;
					++sym;
				}
#ifdef AC_VALIDATE
				lo0=low; r0=range;
#endif
				freq-=cdf;
				errors[2]=sym;
				low+=range*cdf>>16;
				range=(range*freq>>16)-1;
#ifdef AC_VALIDATE
				acval_dec(0, cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
#else
				if(!(range>>16))
				{
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					code=code<<32|*sptr++;
					if(range>~low)
						range=~low;
				}
				r2=code-low;
				c2=(unsigned)_udiv128(r2>>48, r2<<16|0xFFFF, range, &r2);
				curr_CDF=CDF;
				sym=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym+1];
					if((unsigned)freq>c2)
						break;
					++sym;
				}
				freq-=cdf;
				errors[0]=sym;
				r2=range>>16;
				low+=r2*cdf+((range&0xFFFF)*cdf>>16);
				range=r2*freq+((range&0xFFFF)*freq>>16)-1;

				if(!(range>>16))
				{
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					code=code<<32|*sptr++;
					if(range>~low)
						range=~low;
				}
				r2=code-low;
				c2=(unsigned)_udiv128(r2>>48, r2<<16|0xFFFF, range, &r2);
				curr_CDF=CDF+257;
				sym=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym+1];
					if((unsigned)freq>c2)
						break;
					++sym;
				}
				freq-=cdf;
				errors[1]=sym;
				r2=range>>16;
				low+=r2*cdf+((range&0xFFFF)*cdf>>16);
				range=r2*freq+((range&0xFFFF)*freq>>16)-1;

				if(!(range>>16))
				{
					low<<=32;
					range=range<<32|0xFFFFFFFF;
					code=code<<32|*sptr++;
					if(range>~low)
						range=~low;
				}
				r2=code-low;
				c2=(unsigned)_udiv128(r2>>48, r2<<16|0xFFFF, range, &r2);
				curr_CDF=CDF+257+513;
				sym=0;
				for(freq=0;;)
				{
					cdf=freq;
					freq=curr_CDF[sym+1];
					if((unsigned)freq>c2)
						break;
					++sym;
				}
				freq-=cdf;
				errors[2]=sym;
				r2=range>>16;
				low+=r2*cdf+((range&0xFFFF)*cdf>>16);
				range=r2*freq+((range&0xFFFF)*freq>>16)-1;
#endif
				yuv[0]=errors[0]>>1^-(errors[0]&1);//unpack sign
				yuv[1]=errors[1]>>1^-(errors[1]&1);
				yuv[2]=errors[2]>>1^-(errors[2]&1);
				//if(!ky&&kx&&yuv[0])
				//	printf("");
				yuv[0]+=preds[0];
				yuv[1]+=preds[1];
				yuv[2]+=preds[2];
				yuv[0]=yuv[0]<<(32-8)>>(32-8);
				yuv[1]=yuv[1]<<(32-9)>>(32-9);
				yuv[2]=yuv[2]<<(32-9)>>(32-9);
				curr[0]=yuv[0];
				curr[1]=yuv[1];
				curr[2]=yuv[2];
				curr[4]=yuv[0]-preds[0];
				curr[5]=yuv[1]-preds[1];
				curr[6]=yuv[2]-preds[2];
#ifdef USE_RCTJ2K
				yuv[0]-=(yuv[1]+yuv[2])>>2;
				yuv[2]+=yuv[0];
				yuv[1]+=yuv[0];
#else
				//inverse Pei09 RCT
				yuv[0]-=(86*yuv[2]+29*yuv[1]+128)>>8;
				yuv[2]+=yuv[0];
				yuv[1]+=(87*yuv[2]+169*yuv[0]+128)>>8;
#endif
				buffer[idx+1]=(unsigned char)(yuv[0]+128);
				buffer[idx+2]=(unsigned char)(yuv[1]+128);
				buffer[idx+0]=(unsigned char)(yuv[2]+128);
			}

			++hist[errors[0]];
			++hist[256];
			++hist[257+errors[1]];
			++hist[257+512];
			++hist[257+513+errors[2]];
			++hist[257+513+512];
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
		if(fwd)
		{
#ifdef USE_AC48
			*sptr++=(Emit_t)(low>>32);
			*sptr++=(Emit_t)(low>>16);
			*sptr++=(Emit_t)low;
#else
			*sptr++=(unsigned)(low>>32);//AC is big-endian
			*sptr++=(unsigned)low;
#endif

			streamsize=(int)(sptr-sbuf);//number of emits
			if(streamsize<<1>rowsize+4&&0)//bypass row
			{
				streamsize=0;
				fwrite(&streamsize, 1, 4, fdst);
				fwrite(buffer+rowsize*ky, 1, rowsize, fdst);
			}
			else
			{
				fwrite(&streamsize, 1, 4, fdst);
				fwrite(sbuf, sizeof(Emit_t), streamsize, fdst);
			}
		}

		update_CDF(hist, 256, CDF);
		update_CDF(hist+257, 512, CDF+257);
		update_CDF(hist+257+513, 512, CDF+257+513);
	}
	if(!fwd)
		fwrite(buffer, 1, usize, fdst);
	fclose(fsrc);
	fclose(fdst);
	_mm_free(pixels);
	free(hist);
	free(CDF);
	free(sbuf);
	dstsize=get_filesize(dstfn);
#ifdef ESTIMATE_SIZE
	t=time_sec()-t;
	if(fwd)
	{
		csizes[0]/=8;
		csizes[1]/=8;
		csizes[2]/=8;
		printf("%.2lf + %.2lf + %.2lf = %.2lf  ->  %td\n",
			csizes[0],
			csizes[1],
			csizes[2],
			csizes[0]+csizes[1]+csizes[2],
			dstsize
		);
	}
	printf("%c  %lf sec  %lf MB/s\n", 'D'+fwd, t, usize/(t*1024*1024));
#endif
#ifdef BYPASS_ON_INFLATION
	if(fwd&&dstsize>usize)
	{
		printf("Bypass encode\n");
		fdst=fopen(dstfn, "wb");
		fwrite("BP", 1, 2, fdst);
		fwrite(&iw, 1, 4, fdst);
		fwrite(&ih, 1, 4, fdst);
		fwrite(buffer, 1, usize, fdst);
		fclose(fdst);
	}
#endif
	free(buffer);
	(void)nwritten;
	return 0;
}