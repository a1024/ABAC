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

//	#define ENABLE_MT
//	#define ENABLE_GUIDE

//	#define USE_FASTPRED	//useless-fast
	#define USE_ROWPRED	//good
	#define USE_O1		//good


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


//https://github.com/samtools/htscodecs
unsigned char *rans_compress_O0_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_O0_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int out_sz);
unsigned char *rans_compress_O1_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_O1_32x16_avx2(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int out_sz);
#ifdef USE_O1
#define ENTROPY_ENC rans_compress_O1_32x16_avx2
#define ENTROPY_DEC rans_uncompress_O1_32x16_avx2
#else
#define ENTROPY_ENC rans_compress_O0_32x16_avx2
#define ENTROPY_DEC rans_uncompress_O0_32x16_avx2
#endif

#ifdef ENABLE_MT
typedef struct _ThreadArgs
{
	unsigned char *in, *out, *ret;
	unsigned insize, outsize;
} ThreadArgs;
static void c20_enc(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	args->ret=ENTROPY_ENC(args->in, args->insize, args->out, &args->outsize);
}
static void c20_dec(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	args->ret=ENTROPY_DEC(args->in, args->insize, args->out, args->outsize);
}
#endif
int c20_codec(const char *srcfn, const char *dstfn)
{
#if defined _MSC_VER && defined LOUD
	double ptime=0, etime=0;
	double elapsed=time_sec();
	unsigned long long cycles=__rdtsc();
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
	
	if(fwd)//encode
	{
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
		unsigned char *dptr=(unsigned char*)dstptr;
		int res=iw*ih;
#if defined _MSC_VER && defined LOUD
		ptime=time_sec();
#endif
		for(int ky=0, idx=0, idx2=0;ky<ih;++ky)
		{
			ALIGN(16) short *rows[]=
			{
				pixels+((iw+32LL)*((ky-0LL)&1)+16LL)*4,
				pixels+((iw+32LL)*((ky-1LL)&1)+16LL)*4,
			};
			int kx=0;
#if 1
			__m128i getyuv=_mm_set_epi8(
			//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
				-1,  9, 11, 10, -1,  6,  8,  7, -1,  3,  5,  4, -1,  0,  2,  1
			);
			__m128i half=_mm_set_epi8(
				0, -128, -128, -128,
				0, -128, -128, -128,
				0, -128, -128, -128,
				0, -128, -128, -128
			);
			__m256i getoffset=_mm256_set_epi8(
			//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
				-1, -1,  9,  8,  9,  8, -1, -1, -1, -1,  1,  0,  1,  0, -1, -1,
				-1, -1,  9,  8,  9,  8, -1, -1, -1, -1,  1,  0,  1,  0, -1, -1
			);
			__m256i deinter=_mm256_set_epi8(
			//	15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
				-1, -1, -1, -1, 12,  4, -1, -1, 10,  2, -1, -1,  8,  0, -1, -1,
				-1, -1, -1, -1, -1, -1, 12,  4, -1, -1, 10,  2, -1, -1,  8,  0
			);
			__m256i mmin=_mm256_set1_epi16(-128);
			__m256i mmax=_mm256_set1_epi16(127);
			ALIGN(16) int errors[4]={0};
			for(;kx<iw-3;kx+=4, idx+=3*4, idx2+=4)
			{
				//if(kx==4&&ky==1)//
				//	printf("");

				__m128i myuv8=_mm_loadu_si128((__m128i*)(image+idx));
				myuv8=_mm_shuffle_epi8(myuv8, getyuv);
				myuv8=_mm_xor_si128(myuv8, half);
				__m256i myuv=_mm256_cvtepi8_epi16(myuv8);
				__m256i moffset=_mm256_shuffle_epi8(myuv, getoffset);
				_mm256_store_si256((__m256i*)rows[0], _mm256_sub_epi16(myuv, moffset));
				
#ifdef USE_FASTPRED
				__m256i mp=_mm256_loadu_si256((__m256i*)(rows[1]+0*4));//N
#elif defined USE_ROWPRED
				//(2*N+NW+NE+2)>>2
				__m256i mNW	=_mm256_loadu_si256((__m256i*)(rows[1]-1*4));
				__m256i mN	=_mm256_loadu_si256((__m256i*)(rows[1]+0*4));
				__m256i mNE	=_mm256_loadu_si256((__m256i*)(rows[1]+1*4));
				__m256i mp=_mm256_add_epi16(mNW, mNE);
				mp=_mm256_add_epi16(mp, _mm256_slli_epi16(mN, 1));
				mp=_mm256_add_epi16(mp, _mm256_set1_epi16(2));
				mp=_mm256_srai_epi16(mp, 2);
#else
				__m256i mNW	=_mm256_loadu_si256((__m256i*)(rows[1]-1*4));
				__m256i mN	=_mm256_loadu_si256((__m256i*)(rows[1]+0*4));
				__m256i mNE	=_mm256_loadu_si256((__m256i*)(rows[1]+1*4));
				__m256i mW	=_mm256_loadu_si256((__m256i*)(rows[0]-1*4));
				__m256i vmin=_mm256_min_epi16(mN, mW);
				__m256i vmax=_mm256_max_epi16(mN, mW);
				vmin=_mm256_min_epi16(vmin, mNE);
				vmax=_mm256_max_epi16(vmax, mNE);
				__m256i mp=_mm256_slli_epi16(_mm256_add_epi16(mN, mW), 2);
				mp=_mm256_add_epi16(mp, mNE);
				mp=_mm256_sub_epi16(mp, mNW);
				mp=_mm256_add_epi16(mp, _mm256_set1_epi16(4));
				mp=_mm256_srai_epi16(mp, 3);
				mp=_mm256_max_epi16(mp, vmin);
				mp=_mm256_min_epi16(mp, vmax);
#endif

				mp=_mm256_add_epi16(mp, moffset);
				mp=_mm256_max_epi16(mp, mmin);
				mp=_mm256_min_epi16(mp, mmax);

				myuv=_mm256_sub_epi16(myuv, mp);

				myuv=_mm256_shuffle_epi8(myuv, deinter);
				myuv8=_mm_or_si128(_mm256_extracti128_si256(myuv, 0), _mm256_extracti128_si256(myuv, 1));
				_mm_store_si128((__m128i*)errors, myuv8);
				*(unsigned*)(dptr+res*0+idx2)=errors[0];
				*(unsigned*)(dptr+res*1+idx2)=errors[1];
				*(unsigned*)(dptr+res*2+idx2)=errors[2];

				rows[0]+=4*4;
				rows[1]+=4*4;
			}
#endif
			for(;kx<iw;++kx, idx+=3, ++idx2)
			{
				//if(ky==0&&kx==6&&kc==0)//
				//if(ky==399&&kx==714)//
				//if(kx==4&&ky==1)//
				//	printf("");

				int yuv[]=
				{
					image[idx+1]-128,
					image[idx+2]-128,
					image[idx+0]-128,
				};
#ifdef USE_FASTPRED
				int ypred=rows[1][+0*4+0];
				int upred=rows[1][+0*4+1];
				int vpred=rows[1][+0*4+2];
				short *pcurr=rows[0]+0*4;
#elif defined USE_ROWPRED
				int ypred=(2*rows[1][+0*4+0]+rows[1][-1*4+0]+rows[1][+1*4+0]+2)>>2;
				int upred=(2*rows[1][+0*4+1]+rows[1][-1*4+1]+rows[1][+1*4+1]+2)>>2;
				int vpred=(2*rows[1][+0*4+2]+rows[1][-1*4+2]+rows[1][+1*4+2]+2)>>2;
				short *pcurr=rows[0]+0*4;
#else
				short
					yNW	=rows[1][-1*4+0], uNW	=rows[1][-1*4+1], vNW	=rows[1][-1*4+2],
					yN	=rows[1][+0*4+0], uN	=rows[1][+0*4+1], vN	=rows[1][+0*4+2],
					yNE	=rows[1][+1*4+0], uNE	=rows[1][+1*4+1], vNE	=rows[1][+1*4+2],
					yW	=rows[0][-1*4+0], uW	=rows[0][-1*4+1], vW	=rows[0][-1*4+2],
					*pcurr	=rows[0] +0*4;
				int ymax=yN, ymin=yW;
				int umax=uN, umin=uW;
				int vmax=vN, vmin=vW;
				if(yN<yW)ymin=yN, ymax=yW;
				if(ymin>yNE)ymin=yNE;
				if(ymax<yNE)ymax=yNE;
				if(uN<uW)umin=uN, umax=uW;
				if(umin>uNE)umin=uNE;
				if(umax<uNE)umax=uNE;
				if(vN<vW)vmin=vN, vmax=vW;
				if(vmin>vNE)vmin=vNE;
				if(vmax<vNE)vmax=vNE;

				int ypred=(4*(yN+yW)+yNE-yNW+4)>>3;
				CLAMP2(ypred, ymin, ymax);
				int upred=(4*(uN+uW)+uNE-uNW+4)>>3;
				CLAMP2(upred, umin, umax);
				int vpred=(4*(vN+vW)+vNE-vNW+4)>>3;
				CLAMP2(vpred, vmin, vmax);
#endif
				pcurr[0]=yuv[0];
				pcurr[1]=yuv[1]-yuv[0];
				pcurr[2]=yuv[2]-yuv[0];
				upred+=yuv[0];
				CLAMP2(upred, -128, 127);
				vpred+=yuv[0];
				CLAMP2(vpred, -128, 127);

				dptr[res*0+idx2]=yuv[0]-ypred;
				dptr[res*1+idx2]=yuv[1]-upred;
				dptr[res*2+idx2]=yuv[2]-vpred;

				rows[0]+=4;
				rows[1]+=4;
			}
		}
		unsigned char *cdata[3]={0};
		unsigned cdatasize[3]={0};
#if defined _MSC_VER && defined LOUD
		ptime=time_sec()-ptime;
		etime=time_sec();
#endif
#ifdef ENABLE_MT
		ThreadArgs args[]=
		{
			{dptr+res*0, 0, 0, res, 0},
			{dptr+res*1, 0, 0, res, 0},
			{dptr+res*2, 0, 0, res, 0},
		};
		void *mt=mt_exec(c20_enc, args, sizeof(*args), _countof(args));
		mt_finish(mt);
		cdata[0]=args[0].ret;	cdatasize[0]=args[0].outsize;
		cdata[1]=args[1].ret;	cdatasize[1]=args[1].outsize;
		cdata[2]=args[2].ret;	cdatasize[2]=args[2].outsize;
#else
		cdata[0]=ENTROPY_ENC(dptr+res*0, res, cdata[0], cdatasize+0);
		cdata[1]=ENTROPY_ENC(dptr+res*1, res, cdata[1], cdatasize+1);
		cdata[2]=ENTROPY_ENC(dptr+res*2, res, cdata[2], cdatasize+2);
#endif
		if(!cdata[0]||!cdata[1]||!cdata[2])
		{
			LOG_ERROR("HTS Encode Error");
			return 1;
		}
#if defined _MSC_VER && defined LOUD
		etime=time_sec()-etime;
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

			csize+=fwrite(cdatasize+0, 1, 4, fdst);
			csize+=fwrite(cdatasize+1, 1, 4, fdst);
			csize+=fwrite(cdata[0], 1, cdatasize[0], fdst);
			csize+=fwrite(cdata[1], 1, cdatasize[1], fdst);
			csize+=fwrite(cdata[2], 1, cdatasize[2], fdst);
			
#if defined _MSC_VER && defined LOUD
			csize_actual=csize;
#else
			(void)csize;
#endif
			fclose(fdst);
		}
		free(cdata[0]);
		free(cdata[1]);
		free(cdata[2]);
		//free(cdata);

		free(dstbuf);
	}
	else//decode
	{
		int imsize=3*iw*ih;
		unsigned char *image=(unsigned char*)malloc(imsize);
		if(!image)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		unsigned char *im2=(unsigned char*)malloc(imsize);
		if(!im2)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		int res=iw*ih;
		unsigned char *planes[]=
		{
			im2+res*0,
			im2+res*1,
			im2+res*2,
		};
		unsigned char *ret[3]={0};
		unsigned cdatasize[2]={0};
		memcpy(cdatasize+0, srcptr, 4); srcptr+=4;
		memcpy(cdatasize+1, srcptr, 4); srcptr+=4;
		unsigned char *cdata[3]={0};
		cdata[0]=srcptr;
		cdata[1]=cdata[0]+cdatasize[0];
		cdata[2]=cdata[1]+cdatasize[1];
#if defined _MSC_VER && defined LOUD
		etime=time_sec();
#endif
#ifdef ENABLE_MT
		ThreadArgs args[]=
		{
			{cdata[0], planes[0], 0, cdatasize[0],			res},
			{cdata[1], planes[1], 0, cdatasize[1],			res},
			{cdata[2], planes[2], 0, (unsigned)(srcend-cdata[2]),	res},
		};
		void *mt=mt_exec(c20_dec, args, sizeof(*args), _countof(args));
		mt_finish(mt);
		ret[0]=args[0].ret;
		ret[1]=args[1].ret;
		ret[2]=args[2].ret;
#else
		ret[0]=ENTROPY_DEC(cdata[0], cdatasize[0], planes[0], res);
		ret[1]=ENTROPY_DEC(cdata[1], cdatasize[1], planes[1], res);
		ret[2]=ENTROPY_DEC(cdata[2], (unsigned)(srcend-cdata[2]), planes[2], res);
#endif
		if(!ret[0]||!ret[1]||!ret[2])
		{
			LOG_ERROR("HTS Decode Error");
			return 1;
		}
#if defined _MSC_VER && defined LOUD
		etime=time_sec()-etime;
		ptime=time_sec();
#endif
		for(int ky=0, idx=0, idx2=0;ky<ih;++ky)
		{
			ALIGN(16) short *rows[]=
			{
				pixels+((iw+32LL)*((ky-0LL)&1)+16LL)*4,
				pixels+((iw+32LL)*((ky-1LL)&1)+16LL)*4,
			};
			int yuv[3]={0};
			for(int kx=0;kx<iw;++kx, idx+=3, ++idx2)
			{
				//if(ky==0&&kx==6&&kc==0)//
				//if(ky==399&&kx==714)//
				//	printf("");
#ifdef USE_FASTPRED
				int ypred=rows[1][+0*4+0];
				int upred=rows[1][+0*4+1];
				int vpred=rows[1][+0*4+2];
				short *pcurr=rows[0]+0*4;
#elif defined USE_ROWPRED
				int ypred=(2*rows[1][+0*4+0]+rows[1][-1*4+0]+rows[1][+1*4+0]+2)>>2;
				int upred=(2*rows[1][+0*4+1]+rows[1][-1*4+1]+rows[1][+1*4+1]+2)>>2;
				int vpred=(2*rows[1][+0*4+2]+rows[1][-1*4+2]+rows[1][+1*4+2]+2)>>2;
				short *pcurr=rows[0]+0*4;
#else
				short
					yNW	=rows[1][-1*4+0], uNW	=rows[1][-1*4+1], vNW	=rows[1][-1*4+2],
					yN	=rows[1][+0*4+0], uN	=rows[1][+0*4+1], vN	=rows[1][+0*4+2],
					yNE	=rows[1][+1*4+0], uNE	=rows[1][+1*4+1], vNE	=rows[1][+1*4+2],
					yW	=rows[0][-1*4+0], uW	=rows[0][-1*4+1], vW	=rows[0][-1*4+2],
					*pcurr	=rows[0] +0*4;

				int ymax=yN, ymin=yW;
				int umax=uN, umin=uW;
				int vmax=vN, vmin=vW;
				if(yN<yW)ymin=yN, ymax=yW;
				if(ymin>yNE)ymin=yNE;
				if(ymax<yNE)ymax=yNE;
				if(uN<uW)umin=uN, umax=uW;
				if(umin>uNE)umin=uNE;
				if(umax<uNE)umax=uNE;
				if(vN<vW)vmin=vN, vmax=vW;
				if(vmin>vNE)vmin=vNE;
				if(vmax<vNE)vmax=vNE;

				int ypred=(4*(yN+yW)+yNE-yNW+4)>>3;
				int upred=(4*(uN+uW)+uNE-uNW+4)>>3;
				int vpred=(4*(vN+vW)+vNE-vNW+4)>>3;
				CLAMP2(ypred, ymin, ymax);
				CLAMP2(upred, umin, umax);
				CLAMP2(vpred, vmin, vmax);
#endif

				pcurr[0]=yuv[0]=(char)(planes[0][idx2]+ypred);
				upred+=yuv[0];
				CLAMP2(upred, -128, 127);
				vpred+=yuv[0];
				CLAMP2(vpred, -128, 127);

				yuv[1]=(char)(planes[1][idx2]+upred);
				yuv[2]=(char)(planes[2][idx2]+vpred);
				pcurr[1]=yuv[1]-yuv[0];
				pcurr[2]=yuv[2]-yuv[0];

				image[idx+1]=yuv[0]+128;
				image[idx+2]=yuv[1]+128;
				image[idx+0]=yuv[2]+128;
				guide_check(image, kx, ky);
				rows[0]+=4;
				rows[1]+=4;
			}
		}
#if defined _MSC_VER && defined LOUD
		ptime=time_sec()-ptime;
#endif
		free(im2);
		{
			FILE *fdst=fopen(dstfn, "wb");
			if(!fdst)
			{
				LOG_ERROR("Cannot open \"%s\" for writing", dstfn);
				return 1;
			}
			fprintf(fdst, "P6\n%d %d\n255\n", iw, ih);
			fwrite(image, 1, imsize, fdst);
			fclose(fdst);
		}
		free(image);
	}
	_mm_free(pixels);
	free(srcbuf);
#if defined _MSC_VER && defined LOUD
	cycles=__rdtsc()-cycles;
	elapsed=time_sec()-elapsed;
	ptrdiff_t usize=3LL*iw*ih;
	if(fwd)
		printf(" %9td   /%9td %12.6lf:1\n", csize_actual, usize, (double)usize/csize_actual);
	printf("%c%12lf sec %12lf MB/s  %12lld cycles %12lf C/B  pred %lf sec  EC %lf sec\n",
		'D'+fwd,
		elapsed,
		usize/(elapsed*1024*1024),
		cycles,
		(double)cycles/usize,
		ptime,
		etime
	);
#endif
	return 0;
}