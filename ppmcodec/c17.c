#include"codec.h"
#include"util.h"
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


int c17_codec(const char *srcfn, const char *dstfn, int nthreads0)
{
#if 0
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
	FILE *fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		LOG_ERROR("Cannot open \"%s\"", fdst);
		return 1;
	}
	fwrite(srcbuf, sizeof(char), srcsize, fdst);
	fclose(fdst);
	free(srcbuf);
#endif
#if 1
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
	if(fwd)
	{
		unsigned long long *dstbuf=(unsigned long long*)malloc(sizeof(long long)*iw*ih);//4/3 of image size
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		unsigned long long *dstptr=dstbuf;
		unsigned long long cache;
		int nbits;//enc: number of free bits in cache, dec: number of unread bits in cache
		
		cache=0;
		nbits=sizeof(cache)<<3;
		for(int ky=0, idx=0;ky<ih;++ky)//enc
		{
			int pred=0, error=0;
			for(int kx=0;kx<iw;++kx, idx+=3)
			{
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
				for(int kc=0;kc<3;++kc)
				{
					int nbypass=FLOOR_LOG2(error+1);
					int sym=srcptr[idx+kc]-pred;
					sym=sym<<1^sym>>31;

					int nzeros=sym>>nbypass, bypass=sym&((1<<nbypass)-1);
					if(nzeros>=nbits)//fill the rest of cache with zeros, and flush
					{
						nzeros-=nbits;
						*dstptr++=cache;//flush
						if(nzeros>=(int)(sizeof(cache)<<3))//just flush zeros
						{
							cache=0;
							do
							{
								nzeros-=(sizeof(cache)<<3);
								*dstptr++=cache;//flush
							}
							while(nzeros>(int)(sizeof(cache)<<3));
						}
						cache=0;
						nbits=(sizeof(cache)<<3);
					}
					//now there is room for zeros:  0 <= nzeros < nbits <= 64
					nbits-=nzeros;//emit remaining zeros to cache

					bypass|=1<<nbypass;//append 1 stop bit
					++nbypass;
					if(nbypass>=nbits)//not enough free bits in cache:  fill cache, write to list, and repeat
					{
						nbypass-=nbits;
						cache|=(unsigned long long)bypass>>nbypass;
						bypass&=(1<<nbypass)-1;
						*dstptr++=cache;//flush
						cache=0;
						nbits=sizeof(cache)<<3;
					}
					//now there is room for bypass:  0 <= nbypass < nbits <= 64
					nbits-=nbypass;//emit remaining bypass to cache
					cache|=(unsigned long long)bypass<<nbits;
					
					error=(3*error+sym)>>2;
					pred=srcptr[idx+kc];
				}
			}
		}
		*dstptr++=cache;

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
		fwrite(dstbuf, sizeof(long long), streamsize, fdst);
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
		unsigned long long *sptr=(unsigned long long*)srcptr;
		unsigned long long cache;
		int nbits;//enc: number of free bits in cache, dec: number of unread bits in cache
		
		cache=0;
		nbits=0;
		for(int ky=0, idx=0;ky<ih;++ky)//dec
		{
			int pred=0, error=0;
			for(int kx=0;kx<iw;++kx, idx+=3)
			{
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
				for(int kc=0;kc<3;++kc)
				{
					int nbypass=FLOOR_LOG2(error+1);

					//cache: MSB 00[hhh]ijj LSB		nbits 6->3, h is about to be read (past codes must be cleared from cache)
	
					int sym=0, nleadingzeros=0;
					if(!nbits)//cache is empty
						goto read;
					for(;;)//cache reading loop
					{
						nleadingzeros=nbits-FLOOR_LOG2_P1(cache);//count leading zeros
						nbits-=nleadingzeros;//remove accessed zeros
						sym+=nleadingzeros;

						if(nbits)
							break;
					read://cache is empty
						nbits=sizeof(cache)<<3;
						cache=*sptr++;
						//if(gr_dec_impl_read(ec))
						//	return 0;
					}
					//now  0 < nbits <= 64
					--nbits;
					//now  0 <= nbits < 64
					cache-=1ULL<<nbits;//remove stop bit

					unsigned bypass=0;
					sym<<=nbypass;
					if(nbits<nbypass)
					{
						nbypass-=nbits;
						bypass|=(int)(cache<<nbypass);
						nbits=sizeof(cache)<<3;
						cache=*sptr++;
						//if(gr_dec_impl_read(ec))
						//	return 0;
					}
					if(nbypass)
					{
						nbits-=nbypass;
						bypass|=(int)(cache>>nbits);
						cache&=(1ULL<<nbits)-1;
					}
					sym|=bypass;
					
					error=(3*error+sym)>>2;
					sym=sym>>1^-(sym&1);
					dstbuf[idx+kc]=sym+pred;
					pred=dstbuf[idx+kc];
				}
			}
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
	free(srcbuf);
#endif
	return 0;
}