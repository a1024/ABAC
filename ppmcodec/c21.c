#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
#include<immintrin.h>
static const char file[]=__FILE__;

//	#define LOUD
//	#define ENABLE_GUIDE

	#define USE_RANGECODER

#define NSTREAMS 24

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
	{
		LOG_ERROR("");
		printf("");
	}
}
#else
#define guide_save(...)
#define guide_check(...)
#endif
int c21_codec(const char *srcfn, const char *dstfn)
{
#ifdef LOUD
	double elapsed=time_sec();
	unsigned long long cycles=__rdtsc();
#endif
	if(!srcfn||!dstfn)
	{
		LOG_ERROR("Codec requires both source and destination filenames");
		return 1;
	}
	int iw=0, ih=0, fwd=0;
	ptrdiff_t srcsize, res;
	unsigned char *srcbuf=0, *srcptr=0;
//	unsigned char *srcend=0;
	unsigned char *decbuf=0;
	unsigned char *image=0;
	int headersize=0;
	{
		srcsize=get_filesize(srcfn);
		if(srcsize<1)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		if(srcsize<=2)
		{
			LOG_ERROR("File is empty");
			return 1;
		}
		FILE *fsrc=fopen(srcfn, "rb");
		if(!fsrc)
		{
			LOG_ERROR("Cannot open \"%s\"", srcfn);
			return 1;
		}
		unsigned char header[128];
		fread(header, 1, sizeof(header)-1, fsrc);
		srcptr=header;
		fwd=!memcmp(srcptr, "P6\n", 3);
		if(!fwd&&memcmp(srcptr, "CH", 2))
		{
			LOG_ERROR("Unsupported source file");
			return 1;
		}
		if(fwd)
		{
			srcptr+=3;
			while((unsigned)(*srcptr-'0')<10)
				iw=10*iw+*srcptr++ - '0';
			if(iw<1||*srcptr++ != ' ')
			{
				LOG_ERROR("Unsupported source file");
				return 1;
			}
			while((unsigned)(*srcptr-'0')<10)
				ih=10*ih+*srcptr++ - '0';
			if(memcmp(srcptr, "\n255\n", 5))
			{
				LOG_ERROR("Unsupported source file");
				return 1;
			}
			srcptr+=5;
			if(iw<1||ih<1)
			{
				LOG_ERROR("Unsupported source file");
				return 1;
			}
			res=(ptrdiff_t)iw*ih;

			headersize=(int)(srcptr-header);
			fseek(fsrc, 0, SEEK_SET);
			srcbuf=(unsigned char*)malloc(headersize+sizeof(char[3])*iw*(ih+3LL));
			if(!srcbuf)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			image=srcbuf+headersize+sizeof(char[3])*iw*3;
			memset(srcbuf, 0, image-srcbuf);
			fread(image-headersize, 1, headersize+3LL*res, fsrc);
			srcptr=image;
		//	srcend=image+3LL*res;
		}
		else
		{
			srcptr+=2;
			memcpy(&iw, srcptr, 4); srcptr+=4;
			memcpy(&ih, srcptr, 4); srcptr+=4;
			if(iw<1||ih<1)
			{
				LOG_ERROR("Unsupported source file");
				return 1;
			}
			res=(ptrdiff_t)iw*ih;

			fseek(fsrc, 0, SEEK_SET);
			srcbuf=(unsigned char*)malloc(srcsize+16);
			if(!srcbuf)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			fread(srcbuf, 1, srcsize, fsrc);
			srcbuf[srcsize]=0;
			srcptr=srcbuf+10;//CH header size = 10
		//	srcend=srcbuf+srcsize;
			
			headersize=snprintf((char*)header, sizeof(header)-1, "P6\n%d %d\n255\n", iw, ih);
			decbuf=(unsigned char*)malloc(headersize+sizeof(char[3])*iw*(ih+3LL));
			if(!decbuf)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			memset(decbuf, 0, headersize+iw*9LL);
			image=decbuf+headersize+iw*9LL;
			memcpy(image-headersize, header, headersize);
		}
		fclose(fsrc);
	}
	unsigned char *NNNptr	=image-3*3*iw;
	unsigned char *NNptr	=image-2*3*iw;
	unsigned char *Nptr	=image-1*3*iw;
	unsigned char *currptr	=image+0*3*iw;

	unsigned short hist[3][257]={0}, CDF[3][257]={0};
	for(int ks=0;ks<257;++ks)
	{
		int c=ks<<(12-8);
		CDF[0][ks]=c;
		CDF[1][ks]=c;
		CDF[2][ks]=c;
	}

	ALIGN(32) unsigned low[NSTREAMS]={0}, range[NSTREAMS]={0}, code[NSTREAMS]={0}, code2[NSTREAMS]={0};
	memset(range, -1, sizeof(range));
	//for(int k=0;k<NSTREAMS;++k)
	//	range[k]=0xFFFFFFFF;
	//__m256i mlow0=_mm256_load_si256((__m256i*)low+0);
	//__m256i mlow1=_mm256_load_si256((__m256i*)low+1);
	//__m256i mlow2=_mm256_load_si256((__m256i*)low+2);
	//__m256i mrange0=_mm256_load_si256((__m256i*)range+0);
	//__m256i mrange1=_mm256_load_si256((__m256i*)range+1);
	//__m256i mrange2=_mm256_load_si256((__m256i*)range+2);

	ALIGN(32) unsigned weights[3*8]={0}, errors[3*8]={0}, circlebuf[3*8]={0};
	FILLMEM(weights, 0x10000/8, sizeof(weights), sizeof(int));

	ptrdiff_t maxemitts=res/12;
	unsigned short *dstbufs[NSTREAMS]={0}, *streamptrs[NSTREAMS]={0};
	int streamsizes[NSTREAMS]={0};
	if(fwd)
	{
		guide_save(image, iw, ih);
		for(int k=0;k<NSTREAMS;++k)
		{
			dstbufs[k]=(unsigned short*)malloc(sizeof(short)*maxemitts);
			if(!dstbufs[k])
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			streamptrs[k]=dstbufs[k];
		}
	}
	else
	{
		memcpy(streamsizes, srcptr, sizeof(streamsizes)); srcptr+=sizeof(streamsizes);
		for(int k=0;k<NSTREAMS;++k)
		{
			streamptrs[k]=(unsigned short*)srcptr;
			srcptr+=streamsizes[k]*sizeof(short);
		}
		for(int k=0;k<NSTREAMS;++k)
		{
			code[k]=*streamptrs[k]++;
			code[k]=code[k]<<16|*streamptrs[k]++;
		}
	}
	for(int kp=0;kp<res;++kp)
	{
		unsigned *curr_weights=weights, *curr_errors=errors;
		int update=fwd?((kp&7)==7||kp==res-1):!(kp&7);
		if(update)//decode circlebuf
		{
			if(fwd)
			{
#ifdef __GNUC__
#pragma GCC unroll 48
#endif
				for(int k=0;k<NSTREAMS;++k)
				{
					if(range[k]<0x1000)//enc renorm
					{
						*streamptrs[k]++=low[k]>>16;
						low[k]<<=16;
						range[k]=range[k]<<16|0xFFFF;
						unsigned rmax=~low[k];
						if(range[k]>rmax)
							range[k]=rmax;
					}
				}
			}
			else
			{
				//if(kp==5760)//
				//if(kp==5752)//
				//	printf("");//
#ifdef __GNUC__
#pragma GCC unroll 48
#endif
				for(int k=0;k<NSTREAMS;++k)
				{
					if(range[k]<0x1000)//dec renorm
					{
						code[k]=code[k]<<16|*streamptrs[k]++;
						low[k]<<=16;
						range[k]=range[k]<<16|0xFFFF;
						unsigned rmax=~low[k];
						if(range[k]>rmax)
							range[k]=rmax;
					}
				}
#ifdef __GNUC__
#pragma GCC unroll 48
#endif
				for(int k=0;k<NSTREAMS;++k)
#ifdef USE_RANGECODER
				//	code2[k]=code[k]-low[k];
					code2[k]=(code[k]-low[k])/(range[k]>>12);
#else
					code2[k]=(unsigned)(((unsigned long long)(code[k]-low[k])<<12|0xFFF)/range[k]);
#endif
#ifdef __GNUC__
#pragma GCC unroll 48
#endif
				for(int k=0;k<NSTREAMS;++k)
				{
					int sym;

					sym=0;
					for(;;)//dec search
					{
					//	if((range[k]>>12)*CDF[k>>3][sym+2]>code2[k])
						if(CDF[k>>3][sym+2]>code2[k])
						{
						//	sym+=(range[k]>>12)*CDF[k>>3][sym+1]<=code2[k];
							sym+=CDF[k>>3][sym+1]<=code2[k];
							break;
						}
						sym+=2;
#ifdef _DEBUG
						if(sym>255)
							LOG_ERROR("Decode error at %d", kp);
#endif
					}
					circlebuf[k]=sym;
				}
			}
		}
		int offset=0;
		for(int kc=0;kc<3;++kc)
		{
			int
				NNN	=NNNptr		[kc+0*3],
				NN	=NNptr		[kc+0*3],
				NNE	=NNptr		[kc+1*3],
				NW	=Nptr		[kc-1*3],
				N	=Nptr		[kc+0*3],
				NE	=Nptr		[kc+1*3],
				NEE	=Nptr		[kc+2*3],
				NEEE	=Nptr		[kc+3*3],
				NEEEE	=Nptr		[kc+4*3],
				WWWW	=currptr	[kc-4*3],
				WWW	=currptr	[kc-3*3],
				WW	=currptr	[kc-2*3],
				W	=currptr	[kc-1*3];
			if(kc)
			{
				NNN	-=NNNptr	[kc-1+0*3];
				NN	-=NNptr		[kc-1+0*3];
				NNE	-=NNptr		[kc-1+1*3];
				NW	-=Nptr		[kc-1-1*3];
				N	-=Nptr		[kc-1+0*3];
				NE	-=Nptr		[kc-1+1*3];
				NEE	-=Nptr		[kc-1+2*3];
				NEEE	-=Nptr		[kc-1+3*3];
				NEEEE	-=Nptr		[kc-1+4*3];
				WWWW	-=currptr	[kc-1-4*3];
				WWW	-=currptr	[kc-1-3*3];
				WW	-=currptr	[kc-1-2*3];
				W	-=currptr	[kc-1-1*3];
			}
			int preds[]=
			{
				N,
				W,
				3*(N-NN)+NNN,
				3*(W-WW)+WWW,
				W+NE-N,
				(WWWW+WWW+NNN+NEE+NEEE+NEEEE-2*NW+2)>>2,
				N+W-NW,
				N+NE-NNE,
			};
			int vmax=N, vmin=W;
			if(N<W)vmin=N, vmax=W;
			if(vmin<NE)vmin=NE;
			if(vmax<NE)vmax=NE;
			int pred=(
				+curr_weights[0]*preds[0]
				+curr_weights[1]*preds[1]
				+curr_weights[2]*preds[2]
				+curr_weights[3]*preds[3]
				+curr_weights[4]*preds[4]
				+curr_weights[5]*preds[5]
				+curr_weights[6]*preds[6]
				+curr_weights[7]*preds[7]
				+0x8000
			)>>16;
			CLAMP2(pred, vmin, vmax);
			pred+=offset;
			CLAMP2(pred, 0, 255);

			int curr, error;
			//if(kp==5762)//
			//	printf("");
			if(fwd)
			{
				curr=currptr[kc];
				error=(char)(curr-pred);
				error=error<<1^error>>31;//pack sign
				circlebuf[kc<<3|(kp&7)]=error;
			}
			else
			{
				error=circlebuf[kc<<3|(kp&7)];
				curr=error>>1^-(error&1);//unpack sign
				curr=(unsigned char)(curr+pred);
				currptr[kc]=curr;
			}
			++hist[kc][error];
			++hist[kc][256];

			curr_errors[0]+=abs(preds[0]-curr);
			curr_errors[1]+=abs(preds[1]-curr);
			curr_errors[2]+=abs(preds[2]-curr);
			curr_errors[3]+=abs(preds[3]-curr);
			curr_errors[4]+=abs(preds[4]-curr);
			curr_errors[5]+=abs(preds[5]-curr);
			curr_errors[6]+=abs(preds[6]-curr);
			curr_errors[7]+=abs(preds[7]-curr);
			if(!(kp&127))
			{
				unsigned w0, w1, w2, w3, w4, w5, w6, w7, wsum;
				w0=0x8000/(curr_errors[0]+1);//FIXME try curr_weights[0]/(curr_errors[0]+1)
				w1=0x8000/(curr_errors[1]+1);
				w2=0x8000/(curr_errors[2]+1);
				w3=0x8000/(curr_errors[3]+1);
				w4=0x8000/(curr_errors[4]+1);
				w5=0x8000/(curr_errors[5]+1);
				w6=0x8000/(curr_errors[6]+1);
				w7=0x8000/(curr_errors[7]+1);
				wsum=w0+w1+w2+w3+w4+w5+w6+w7+1;
				curr_weights[0]=(w0<<16)/wsum;
				curr_weights[1]=(w1<<16)/wsum;
				curr_weights[2]=(w2<<16)/wsum;
				curr_weights[3]=(w3<<16)/wsum;
				curr_weights[4]=(w4<<16)/wsum;
				curr_weights[5]=(w5<<16)/wsum;
				curr_weights[6]=(w6<<16)/wsum;
				curr_weights[7]=0x10000-(
					+curr_weights[0]
					+curr_weights[1]
					+curr_weights[2]
					+curr_weights[3]
					+curr_weights[4]
					+curr_weights[5]
					+curr_weights[6]
				);
#ifdef _DEBUG
				if(
					(int)curr_weights[0]<0||
					(int)curr_weights[1]<0||
					(int)curr_weights[2]<0||
					(int)curr_weights[3]<0||
					(int)curr_weights[4]<0||
					(int)curr_weights[5]<0||
					(int)curr_weights[6]<0||
					(int)curr_weights[7]<0
				)
					LOG_ERROR("");
#endif
				memset(curr_errors, 0, sizeof(int[8]));//FIXME try halving errors
			}

			offset=curr;
			curr_weights+=8;
			curr_errors+=8;
		}
		if(!fwd)
			guide_check(image, kp%iw, kp/iw);//

		if(update)
		{
			//if(kp==5759)//
			//	printf("");

			int cdf[3*8]=
			{
				CDF[0][circlebuf[0+0*8]],
				CDF[0][circlebuf[1+0*8]],
				CDF[0][circlebuf[2+0*8]],
				CDF[0][circlebuf[3+0*8]],
				CDF[0][circlebuf[4+0*8]],
				CDF[0][circlebuf[5+0*8]],
				CDF[0][circlebuf[6+0*8]],
				CDF[0][circlebuf[7+0*8]],
				CDF[1][circlebuf[0+1*8]],
				CDF[1][circlebuf[1+1*8]],
				CDF[1][circlebuf[2+1*8]],
				CDF[1][circlebuf[3+1*8]],
				CDF[1][circlebuf[4+1*8]],
				CDF[1][circlebuf[5+1*8]],
				CDF[1][circlebuf[6+1*8]],
				CDF[1][circlebuf[7+1*8]],
				CDF[2][circlebuf[0+2*8]],
				CDF[2][circlebuf[1+2*8]],
				CDF[2][circlebuf[2+2*8]],
				CDF[2][circlebuf[3+2*8]],
				CDF[2][circlebuf[4+2*8]],
				CDF[2][circlebuf[5+2*8]],
				CDF[2][circlebuf[6+2*8]],
				CDF[2][circlebuf[7+2*8]],
			};
			int freq[3*8]=
			{
				CDF[0][circlebuf[0+0*8]+1]-cdf[0+0*8],
				CDF[0][circlebuf[1+0*8]+1]-cdf[1+0*8],
				CDF[0][circlebuf[2+0*8]+1]-cdf[2+0*8],
				CDF[0][circlebuf[3+0*8]+1]-cdf[3+0*8],
				CDF[0][circlebuf[4+0*8]+1]-cdf[4+0*8],
				CDF[0][circlebuf[5+0*8]+1]-cdf[5+0*8],
				CDF[0][circlebuf[6+0*8]+1]-cdf[6+0*8],
				CDF[0][circlebuf[7+0*8]+1]-cdf[7+0*8],
				CDF[1][circlebuf[0+1*8]+1]-cdf[0+1*8],
				CDF[1][circlebuf[1+1*8]+1]-cdf[1+1*8],
				CDF[1][circlebuf[2+1*8]+1]-cdf[2+1*8],
				CDF[1][circlebuf[3+1*8]+1]-cdf[3+1*8],
				CDF[1][circlebuf[4+1*8]+1]-cdf[4+1*8],
				CDF[1][circlebuf[5+1*8]+1]-cdf[5+1*8],
				CDF[1][circlebuf[6+1*8]+1]-cdf[6+1*8],
				CDF[1][circlebuf[7+1*8]+1]-cdf[7+1*8],
				CDF[2][circlebuf[0+2*8]+1]-cdf[0+2*8],
				CDF[2][circlebuf[1+2*8]+1]-cdf[1+2*8],
				CDF[2][circlebuf[2+2*8]+1]-cdf[2+2*8],
				CDF[2][circlebuf[3+2*8]+1]-cdf[3+2*8],
				CDF[2][circlebuf[4+2*8]+1]-cdf[4+2*8],
				CDF[2][circlebuf[5+2*8]+1]-cdf[5+2*8],
				CDF[2][circlebuf[6+2*8]+1]-cdf[6+2*8],
				CDF[2][circlebuf[7+2*8]+1]-cdf[7+2*8],
			};
#ifdef USE_RANGECODER
#ifdef __GNUC__
#pragma GCC unroll 48
#endif
			for(int k=0;k<NSTREAMS;++k)
				low[k]+=(range[k]>>12)*cdf[k];
#ifdef __GNUC__
#pragma GCC unroll 48
#endif
			for(int k=0;k<NSTREAMS;++k)
				range[k]=(range[k]>>12)*freq[k]-1;
#else
#ifdef __GNUC__
#pragma GCC unroll 48
#endif
			for(int k=0;k<NSTREAMS;++k)
				low[k]+=(unsigned)((unsigned long long)range[k]*cdf[k]>>12);
#ifdef __GNUC__
#pragma GCC unroll 48
#endif
			for(int k=0;k<NSTREAMS;++k)
				range[k]=(unsigned)((unsigned long long)range[k]*freq[k]>>12)-1;
#endif
		}
		for(int kc=0;kc<3;++kc)
		{
			unsigned short *curr_hist=hist[kc];
			if(curr_hist[256]==4096-256)//snapshot-CDF
			{
				//if(kp==5765&&kc==1)//
				//	printf("");
				unsigned short *curr_CDF=CDF[kc];
				//int hsum=0;
				for(int ks=0, sum=0;ks<256;++ks)
				{
					curr_CDF[ks]=sum+ks;
					sum+=curr_hist[ks];
					curr_hist[ks]=0;
					//hsum+=curr_hist[ks]>>=1;
				}
				curr_hist[256]=0;
				//curr_hist[256]=hsum;
			}
		}
		NNptr	+=3;
		Nptr	+=3;
		currptr	+=3;
	}
	ptrdiff_t usize=res*3, csize_actual=0;
	FILE *fdst=fopen(dstfn, "wb");
	if(!fdst)
	{
		LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
		return 1;
	}
	if(fwd)
	{
		for(int k=0;k<NSTREAMS;++k)//flush
		{
			*streamptrs[k]++=low[k]>>16;
			*streamptrs[k]++=low[k];
		}
		csize_actual+=fwrite("CH", 1, 2, fdst);
		csize_actual+=fwrite(&iw, 1, 4, fdst);
		csize_actual+=fwrite(&ih, 1, 4, fdst);
		for(int k=0;k<NSTREAMS;++k)
			streamsizes[k]=(int)(streamptrs[k]-dstbufs[k]);
		csize_actual+=fwrite(streamsizes, 1, sizeof(int[NSTREAMS]), fdst);
		for(int k=0;k<NSTREAMS;++k)
			csize_actual+=fwrite(dstbufs[k], 1, streamsizes[k]*sizeof(short), fdst);

		for(int k=0;k<NSTREAMS;++k)
			free(dstbufs[k]);
	}
	else
	{
		fwrite(image-headersize, 1, headersize+3*res, fdst);
		free(decbuf);
	}
	fclose(fdst);
	free(srcbuf);

#ifdef LOUD
	cycles=__rdtsc()-cycles;
	elapsed=time_sec()-elapsed;
	if(fwd)
		printf(" %9td   /%9td %12.6lf:1\n", csize_actual, usize, (double)usize/csize_actual);
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