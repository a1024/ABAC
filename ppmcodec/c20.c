#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
//#include"blist.h"
static const char file[]=__FILE__;


//	#define ESTIMATE_SIZE
//	#define ENABLE_GUIDE
//	#define ENABLE_VALIDATION
//	#define COPY_TEST

//	#define USE_CTR


#define ABAC_SH 6
//#define ABAC_RISK 16

#ifdef ENABLE_GUIDE
static int g_iw=0, g_ih=0;
static unsigned char *g_image=0;
static void guide_save(unsigned char *image, int iw, int ih)
{
	int size=3*iw*ih;
	g_iw=iw;
	g_ih=ih;
	g_image=(unsigned char*)malloc(size);
	if(!g_image)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(g_image, image, size);
}
static void guide_check(unsigned char *image, int kx, int ky)
{
	int idx=3*(g_iw*ky+kx);
	if(memcmp(image+idx, g_image+idx, 3))
		LOG_ERROR("");
}
#else
#define guide_save(...)
#define guide_check(...)
#endif


#ifdef ENABLE_VALIDATION
static ArrayHandle g_probs=0;
static ptrdiff_t g_probidx=0;
static void val_append(unsigned short prob)
{
	if(!g_probs)
		ARRAY_ALLOC(short, g_probs, 0, 0, 0, 0);
	ARRAY_APPEND(g_probs, &prob, 1, 1, 0);
}
static void val_check(unsigned short prob)
{
	unsigned short *correctprob=(unsigned short*)array_at(&g_probs, g_probidx);
	if(prob!=*correctprob)
	{
		printf("%8td %04X %04X\n", g_probidx, prob, *correctprob);
		LOG_ERROR("Val Error");
	}
	++g_probidx;
}
#else
#define val_append(...)
#define val_check(...)
#endif


int c20_codec(const char *srcfn, const char *dstfn)
{
#ifdef _MSC_VER
	double elapsed=time_sec();
	unsigned long long cycles=__rdtsc();
#endif
#ifdef ESTIMATE_SIZE
	double csizes[3]={0};
#endif
#if defined ESTIMATE_SIZE || defined _MSC_VER
	ptrdiff_t csize_actual=0;
#endif
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	ptrdiff_t srcsize;
	unsigned char *srcbuf=0, *srcptr=0, *srcend=0;
	{
		srcsize=get_filesize(srcfn);
		if(srcsize<1)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		srcbuf=(unsigned char*)malloc(srcsize+16);
		if(!srcbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		fread(srcbuf, 1, srcsize, fsrc);
		fclose(fsrc);
		srcbuf[srcsize]=0;
	}
	srcptr=srcbuf;
	srcend=srcbuf+srcsize;
	if(srcsize<=2)
	{
		LOG_ERROR("File is empty");
		return 1;
	}
	int tag=*(unsigned short*)srcptr;
	srcptr+=2;
	int fwd=tag==('P'|'6'<<8);
	if(!fwd&&tag!=('C'|'H'<<8))
	{
		LOG_ERROR("Unsupported source file");
		return 1;
	}
	int iw=0, ih=0;
	//unsigned char *image=0;
	//ptrdiff_t dstbufsize=0;
	//unsigned char *dstbuf;
	if(fwd)//encode
	{
		if(*srcptr++ != '\n')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		while((unsigned)(*srcptr-'0')<10)
			iw=10*iw+*srcptr++ - '0';
		if(iw<1||*srcptr++ != ' ')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		while((unsigned)(*srcptr-'0')<10)
			ih=10*ih+*srcptr++ - '0';
		if(ih<1||*srcptr++ != '\n')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		int temp=0;
		while((unsigned)(*srcptr-'0')<10)
			temp=10*temp+*srcptr++ - '0';
		if(temp!=255||*srcptr++ != '\n')
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
	}
	else//decode
	{
		if(srcptr+4*3>=srcend)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		memcpy(&iw, srcptr, 4); srcptr+=4;
		memcpy(&ih, srcptr, 4); srcptr+=4;
		if(iw<1||ih<1)
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
	}
	
	int psize=(iw+32LL)*sizeof(short[4*3]);//4 padded rows * 3 channels
	short *pixels=(short*)_mm_malloc(psize, sizeof(__m256i));
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(pixels, 0, psize);
	
	unsigned short stats[3][256]={0};
#ifndef USE_CTR
	FILLMEM((unsigned short*)stats, 0x8000, sizeof(stats), sizeof(short));
#endif
	
	unsigned long long low=0, range=0xFFFFFFFFFFFF;
	if(fwd)//encode
	{
#ifndef COPY_TEST
		ptrdiff_t dstbufsize=4LL*iw*ih;
		unsigned short *dstbuf=0, *dstptr=0;
		dstbuf=(unsigned short*)malloc(dstbufsize);
		if(!dstbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		dstptr=dstbuf;
		unsigned char *image=srcptr;
		guide_save(image, iw, ih);
		for(int ky=0, idx=0;ky<ih;++ky)
		{
			ALIGN(16) short *rows[]=
			{
				pixels+((iw+32LL)*((ky-0LL)&1)+16LL)*4,
				pixels+((iw+32LL)*((ky-1LL)&1)+16LL)*4,
			};
			for(int kx=0;kx<iw;++kx, idx+=3)
			{
				int yuv[]=
				{
					image[idx+1]-128,
					image[idx+2]-128,
					image[idx+0]-128,
				};
				int offset=0;
				for(int kc=0;kc<3;++kc)
				{
					//if(ky==0&&kx==6&&kc==0)//
					//	printf("");

					short
						NW	=rows[1][-1*4+kc],
						N	=rows[1][+0*4+kc],
						W	=rows[0][-1*4+kc],
						*pcurr	=rows[0] +0*4 ;
					int vmax=N, vmin=W; if(N<W)vmin=N, vmax=W;
					int pred=N+W-NW; CLAMP2(pred, vmin, vmax);
					pred+=offset; CLAMP2(pred, -128, 127);
					int error=yuv[kc]-pred;
#ifdef __GNUC__
#pragma GCC unroll 8
#elif defined __clang__
#pragma clang loop unroll_count(8)
#endif
					for(int kb=7, tidx=1;kb>=0;--kb)
					{
						int bit=error>>kb&1;
#ifdef USE_CTR
						unsigned char *ctr=(unsigned char*)(stats[kc]+tidx);
						int p0=(ctr[0]+ctr[1])*3+2;
						p0=(((ctr[0]*3+1)<<16)+(p0>>1))/p0;
						CLAMP2(p0, 1, 0xFFFF);
#else
						unsigned short *pp0=stats[kc]+tidx, p0=*pp0;
						//CLAMP2(p0, ABAC_RISK, 0x10000-ABAC_RISK);
#endif
						val_append(p0);
						while(range<0x10000)
						{
							*dstptr++=(unsigned short)(low>>32);
							low=low<<16&0xFFFFFFFF0000;
							range=range<<16|0xFFFF;
							unsigned long long rmax=low^0xFFFFFFFFFFFF;
							if(range>rmax)
								range=rmax;
						}
						unsigned long long mid=range*p0>>16;
						if(bit)
						{
							low+=mid;
							range-=mid;
						}
						else
							range=mid-1;
#ifdef ESTIMATE_SIZE
						csizes[kc]-=log2((double)(bit?0x10000-p0:p0)/0x10000);
#endif
#ifdef USE_CTR
						if(ctr[bit]>=255)
						{
							ctr[0]>>=1;
							ctr[1]>>=1;
						}
						++ctr[bit];
#else
						p0+=((!bit<<16)-p0+(1<<ABAC_SH>>1))>>ABAC_SH;
						p0+=((!bit<<16)-p0+(1<<ABAC_SH>>1))>>ABAC_SH;
						*pp0=p0;
#endif
						tidx=2*tidx+bit;
					}
					pcurr[kc]=yuv[kc]-offset;
					offset=yuv[0];
				}
				rows[0]+=4;
				rows[1]+=4;
			}
		}
		*dstptr++=(unsigned short)(low>>32);
		*dstptr++=(unsigned short)(low>>16);
		*dstptr++=(unsigned short)low;
#endif
		{
			ptrdiff_t csize=0;
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
				return 1;
			}
			csize+=fwrite("CH", 1, 2, fdst);
			csize+=fwrite(&iw, 1, 4, fdst);
			csize+=fwrite(&ih, 1, 4, fdst);
#ifndef COPY_TEST
			ptrdiff_t size=dstptr-dstbuf;
			csize+=fwrite(&size, 1, 4, fdst);
			csize+=fwrite(dstbuf, 1, 2LL*size, fdst);
#else
			csize+=fwrite(srcptr, 1, 3LL*iw*ih, fdst);
#endif
#if defined ESTIMATE_SIZE || defined _MSC_VER
			csize_actual=csize;
#endif
			fclose(fdst);
		}
#ifndef COPY_TEST
		free(dstbuf);
#endif
	}
	else//decode
	{
#ifndef COPY_TEST
		int imsize=3*iw*ih;
		unsigned char *image=(unsigned char*)malloc(imsize);
		if(!image)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		ptrdiff_t nemitts=0;
		memcpy(&nemitts, srcptr, 4); srcptr+=4;
		unsigned short *sptr=(unsigned short*)srcptr;
		unsigned long long code=(unsigned long long)sptr[0]<<32|(unsigned long long)sptr[1]<<16|sptr[2];
		sptr+=3;
		//code=code<<16|*sptr++;
		//code=code<<16|*sptr++;
		for(int ky=0, idx=0;ky<ih;++ky)
		{
			ALIGN(16) short *rows[]=
			{
				pixels+((iw+32LL)*((ky-0LL)&1)+16LL)*4,
				pixels+((iw+32LL)*((ky-1LL)&1)+16LL)*4,
			};
			int yuv[3]={0};
			for(int kx=0;kx<iw;++kx, idx+=3)
			{
				int offset=0;
				for(int kc=0;kc<3;++kc)
				{
					//if(ky==0&&kx==6&&kc==0)//
					//	printf("");

					short
						NW	=rows[1][-1*4+kc],
						N	=rows[1][+0*4+kc],
						W	=rows[0][-1*4+kc],
						*pcurr	=rows[0] +0*4 ;
					int vmax=N, vmin=W; if(N<W)vmin=N, vmax=W;
					int pred=N+W-NW; CLAMP2(pred, vmin, vmax);
					pred+=offset; CLAMP2(pred, -128, 127);
					int error=0;
#ifdef __GNUC__
#pragma GCC unroll 8
#elif defined __clang__
#pragma clang loop unroll_count(8)
#endif
					for(int kb=7, tidx=1;kb>=0;--kb)
					{
						int bit;
#ifdef USE_CTR
						unsigned char *ctr=(unsigned char*)(stats[kc]+tidx);
						int p0=(ctr[0]+ctr[1])*3+2;
						p0=(((ctr[0]*3+1)<<16)+(p0>>1))/p0;
						CLAMP2(p0, 1, 0xFFFF);
#else
						unsigned short *pp0=stats[kc]+tidx, p0=*pp0;
						//CLAMP2(p0, ABAC_RISK, 0x10000-ABAC_RISK);
#endif
						val_check(p0);
						while(range<0x10000)
						{
							code=(code<<16&0xFFFFFFFF0000)|*sptr++;
							low=low<<16&0xFFFFFFFF0000;
							range=range<<16|0xFFFF;
							unsigned long long rmax=low^0xFFFFFFFFFFFF;
							if(range>rmax)
								range=rmax;
						}
						unsigned long long mid=range*p0>>16;
						unsigned long long t2=low+mid;
						range-=mid;
						bit=code>=t2;
						if(bit)
							low=t2;
						else
							range=mid-1;
						error|=bit<<kb;
#ifdef USE_CTR
						if(ctr[bit]>=255)
						{
							ctr[0]>>=1;
							ctr[1]>>=1;
						}
						++ctr[bit];
#else
						p0+=((!bit<<16)-p0+(1<<ABAC_SH>>1))>>ABAC_SH;
						p0+=((!bit<<16)-p0+(1<<ABAC_SH>>1))>>ABAC_SH;
						*pp0=p0;
#endif
						tidx=2*tidx+bit;
					}
					error+=pred;
					error<<=24;
					error>>=24;
					yuv[kc]=error;
					pcurr[kc]=yuv[kc]-offset;
					offset=yuv[0];
				}
				image[idx+1]=yuv[0]+128;
				image[idx+2]=yuv[1]+128;
				image[idx+0]=yuv[2]+128;
				guide_check(image, kx, ky);
				rows[0]+=4;
				rows[1]+=4;
			}
		}
#endif
		{
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
				return 1;
			}
			fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
#ifndef COPY_TEST
			fwrite(image, 1, imsize, fdst);
#else
			fwrite(srcptr, 1, 3LL*iw*ih, fdst);
#endif
			fclose(fdst);
		}
#ifndef COPY_TEST
		free(image);
#endif
	}
	_mm_free(pixels);
	free(srcbuf);
#ifdef _MSC_VER
	cycles=__rdtsc()-cycles;
	elapsed=time_sec()-elapsed;
	ptrdiff_t usize=3LL*iw*ih;
#endif
#ifdef ESTIMATE_SIZE
	if(fwd)
	{
		printf("T%12.2lf\n", (csizes[0]+csizes[1]+csizes[2])/8);
		for(int kc=0;kc<3;++kc)
			printf("%c%12.2lf\n", "YUV"[kc], csizes[kc]/8);
		printf(" %9td   /%9td %12.6lf:1\n", csize_actual, usize, (double)usize/csize_actual);
	}
#elif defined _MSC_VER
	if(fwd)
		printf(" %9td   /%9td %12.6lf:1\n", csize_actual, usize, (double)usize/csize_actual);
#endif
#ifdef _MSC_VER
	printf("%c%12lf sec %12lf MB/s  %12lld cycles %12lf C/B\n",
		'D'+fwd,
		elapsed,
		usize/(elapsed*1024*1024),
		cycles,
		(double)cycles/usize
	);
#endif
	return 0;
}