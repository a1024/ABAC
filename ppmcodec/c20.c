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


	#define SILENT

//	#define ENABLE_MT
//	#define ENABLE_GUIDE

//	#define ERROR_CORRECTION
//	#define ERROR_SHIFT 8

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

#ifdef ENABLE_MT
typedef struct _ThreadArgs
{
	unsigned char *in, *out, *ret;
	unsigned insize, outsize;
} ThreadArgs;
static void c20_enc(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	args->ret=rans_compress_O1_32x16_avx2(args->in, args->insize, args->out, &args->outsize);
}
static void c20_dec(void *param)
{
	ThreadArgs *args=(ThreadArgs*)param;
	args->ret=rans_uncompress_O1_32x16_avx2(args->in, args->insize, args->out, args->outsize);
}
#endif
int c20_codec(const char *srcfn, const char *dstfn)
{
#ifdef _MSC_VER
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
		for(int ky=0, idx=0, idx2=0;ky<ih;++ky)
		{
			ALIGN(16) short *rows[]=
			{
				pixels+((iw+32LL)*((ky-0LL)&1)+16LL)*4,
				pixels+((iw+32LL)*((ky-1LL)&1)+16LL)*4,
			};
#ifdef ERROR_CORRECTION
			int eWWW[3]={0}, eWW[3]={0}, eW[3]={0};
#endif
			for(int kx=0;kx<iw;++kx, idx+=3, ++idx2)
			{
				int yuv[]=
				{
					image[idx+1]-128,
					image[idx+2]-128,
					image[idx+0]-128,
				};
				int offset=0;
				//if(ky==0&&kx==6&&kc==0)//
				//if(ky==399&&kx==714)//
				//	printf("");
				{
					short
						NW	=rows[1][-1*4+0],
						N	=rows[1][+0*4+0],
						NE	=rows[1][+1*4+0],
						WW	=rows[0][-2*4+0],
						W	=rows[0][-1*4+0],
						*pcurr	=rows[0] +0*4 ;
					int vmax=N, vmin=W; if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE; if(vmax<NE)vmax=NE; int pred=(4*(N+W)+NE-NW+4)>>3;
					//int pred=N+W-NW;
					//if(vmin>NE)vmin=NE; if(vmax<NE)vmax=NE; int pred=W+((5*(N-NW)+NE-WW+4)>>3);
#ifdef ERROR_CORRECTION
					pred-=(3*(eW[0]+eWW[0])+eWWW[0]+(1<<5>>1))>>5; //CLAMP2(pred, -128, 127);
#endif
					CLAMP2(pred, vmin, vmax);
				//	pred+=offset; CLAMP2(pred, -128, 127);
					int error=yuv[0]-pred;

					dptr[res*0+idx2]=error;
					//*dptr++=error;
					
#ifdef ERROR_CORRECTION
					eWWW[0]=eWW[0];
					eWW[0]=eW[0];
					eW[0]=error;
#endif
					pcurr[0]=yuv[0]-offset;
				}
				offset=yuv[0];
				{
					short
						NW	=rows[1][-1*4+1],
						N	=rows[1][+0*4+1],
						NE	=rows[1][+1*4+1],
						WW	=rows[0][-2*4+1],
						W	=rows[0][-1*4+1],
						*pcurr	=rows[0] +0*4 ;
					int vmax=N, vmin=W; if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE; if(vmax<NE)vmax=NE; int pred=(4*(N+W)+NE-NW+4)>>3;
					//int pred=N+W-NW;
					//if(vmin>NE)vmin=NE; if(vmax<NE)vmax=NE; int pred=W+((5*(N-NW)+NE-WW+4)>>3);
#ifdef ERROR_CORRECTION
					pred-=(3*(eW[1]+eWW[1])+eWWW[1]+(1<<5>>1))>>5; //CLAMP2(pred, -128, 127);
#endif
					CLAMP2(pred, vmin, vmax);
					pred+=offset; CLAMP2(pred, -128, 127);
					int error=yuv[1]-pred;

					dptr[res*1+idx2]=error;
					//*dptr++=error;
					
#ifdef ERROR_CORRECTION
					eWWW[1]=eWW[1];
					eWW[1]=eW[1];
					eW[1]=error;
#endif
					pcurr[1]=yuv[1]-offset;
				}
				{
					short
						NW	=rows[1][-1*4+2],
						N	=rows[1][+0*4+2],
						NE	=rows[1][+1*4+2],
						WW	=rows[0][-2*4+2],
						W	=rows[0][-1*4+2],
						*pcurr	=rows[0] +0*4 ;
					int vmax=N, vmin=W; if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE; if(vmax<NE)vmax=NE; int pred=(4*(N+W)+NE-NW+4)>>3;
					//int pred=N+W-NW;
					//if(vmin>NE)vmin=NE; if(vmax<NE)vmax=NE; int pred=W+((5*(N-NW)+NE-WW+4)>>3);
#ifdef ERROR_CORRECTION
					pred-=(3*(eW[2]+eWW[2])+eWWW[2]+(1<<5>>1))>>5; //CLAMP2(pred, -128, 127);
#endif
					CLAMP2(pred, vmin, vmax);
					pred+=offset; CLAMP2(pred, -128, 127);
					int error=yuv[2]-pred;

					dptr[res*2+idx2]=error;
					//*dptr++=error;
					
#ifdef ERROR_CORRECTION
					eWWW[2]=eWW[2];
					eWW[2]=eW[2];
					eW[2]=error;
#endif
					pcurr[2]=yuv[2]-offset;
				}
				rows[0]+=4;
				rows[1]+=4;
			}
		}
		unsigned char *cdata[3]={0};
		unsigned cdatasize[3]={0};
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
		cdata[0]=rans_compress_O1_32x16_avx2(dptr+res*0, res, cdata[0], cdatasize+0);
		cdata[1]=rans_compress_O1_32x16_avx2(dptr+res*1, res, cdata[1], cdatasize+1);
		cdata[2]=rans_compress_O1_32x16_avx2(dptr+res*2, res, cdata[2], cdatasize+2);
#endif
		if(!cdata[0]||!cdata[1]||!cdata[2])
		{
			LOG_ERROR("HTS Encode Error");
			return 1;
		}
		//unsigned char *cdata=0;
		//unsigned cdatasize=0;
		//cdata=rans_compress_O1_32x16_avx2((unsigned char*)dstbuf, 3*iw*ih, cdata, &cdatasize);//FIXME planar
		//if(!cdata)
		//{
		//	LOG_ERROR("HTS Encode Error");
		//	return 1;
		//}
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
			//csize+=fwrite(cdata, 1, cdatasize, fdst);

#ifdef _MSC_VER
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
		ret[0]=rans_uncompress_O1_32x16_avx2(cdata[0], cdatasize[0], planes[0], res);
		ret[1]=rans_uncompress_O1_32x16_avx2(cdata[1], cdatasize[1], planes[1], res);
		ret[2]=rans_uncompress_O1_32x16_avx2(cdata[2], (unsigned)(srcend-cdata[2]), planes[2], res);
#endif
		if(!ret[0]||!ret[1]||!ret[2])
		{
			LOG_ERROR("HTS Decode Error");
			return 1;
		}
		//unsigned char *ret=rans_uncompress_O1_32x16_avx2(srcptr, (unsigned)(srcend-srcptr), image, imsize);
		//if(!ret)
		//{
		//	LOG_ERROR("HTS Decode Error");
		//	return 1;
		//}
		for(int ky=0, idx=0, idx2=0;ky<ih;++ky)
		{
			ALIGN(16) short *rows[]=
			{
				pixels+((iw+32LL)*((ky-0LL)&1)+16LL)*4,
				pixels+((iw+32LL)*((ky-1LL)&1)+16LL)*4,
			};
#ifdef ERROR_CORRECTION
			int eWWW[3]={0}, eWW[3]={0}, eW[3]={0};
#endif
			int yuv[3]={0};
			for(int kx=0;kx<iw;++kx, idx+=3, ++idx2)
			{
				int offset=0;
				//if(ky==0&&kx==6&&kc==0)//
				//if(ky==399&&kx==714)//
				//	printf("");
				{
					short
						NW	=rows[1][-1*4+0],
						N	=rows[1][+0*4+0],
						NE	=rows[1][+1*4+0],
						WW	=rows[0][-2*4+0],
						W	=rows[0][-1*4+0],
						*pcurr	=rows[0] +0*4 ;
					int vmax=N, vmin=W; if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE; if(vmax<NE)vmax=NE; int pred=(4*(N+W)+NE-NW+4)>>3;
					//int pred=N+W-NW;
					//if(vmin>NE)vmin=NE; if(vmax<NE)vmax=NE; int pred=W+((5*(N-NW)+NE-WW+4)>>3);
#ifdef ERROR_CORRECTION
					pred-=(3*(eW[0]+eWW[0])+eWWW[0]+(1<<5>>1))>>5; //CLAMP2(pred, -128, 127);
#endif
					CLAMP2(pred, vmin, vmax);
				//	pred+=eW[0]>>2; CLAMP2(pred, -128, 127);
				//	pred+=offset; CLAMP2(pred, -128, 127);
					int error=0;

					error=planes[0][idx2];
					//error=image[idx+0];

					error+=pred;
					error<<=24;
					error>>=24;
					yuv[0]=error;
					pcurr[0]=yuv[0]-offset;
					
#ifdef ERROR_CORRECTION
					eWWW[0]=eWW[0];
					eWW[0]=eW[0];
					eW[0]=error-pred;
				//	eW[0]+=(error-pred-eW[0]+(1<<ERROR_SHIFT>>1))>>ERROR_SHIFT;
#endif
				}
				offset=yuv[0];
				{
					//if(ky==0&&kx==6&&kc==0)//
					//	printf("");

					short
						NW	=rows[1][-1*4+1],
						N	=rows[1][+0*4+1],
						NE	=rows[1][+1*4+1],
						WW	=rows[0][-2*4+1],
						W	=rows[0][-1*4+1],
						*pcurr	=rows[0] +0*4 ;
					int vmax=N, vmin=W; if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE; if(vmax<NE)vmax=NE; int pred=(4*(N+W)+NE-NW+4)>>3;
					//int pred=N+W-NW;
					//if(vmin>NE)vmin=NE; if(vmax<NE)vmax=NE; int pred=W+((5*(N-NW)+NE-WW+4)>>3);
#ifdef ERROR_CORRECTION
					pred-=(3*(eW[1]+eWW[1])+eWWW[1]+(1<<5>>1))>>5; //CLAMP2(pred, -128, 127);
#endif
					CLAMP2(pred, vmin, vmax);
				//	pred+=eW[1]>>2; CLAMP2(pred, -128, 127);
					pred+=offset; CLAMP2(pred, -128, 127);
					int error=0;

					error=planes[1][idx2];
					//error=image[idx+1];

					error+=pred;
					error<<=24;
					error>>=24;
					yuv[1]=error;
					pcurr[1]=yuv[1]-offset;
					
#ifdef ERROR_CORRECTION
					eWWW[1]=eWW[1];
					eWW[1]=eW[1];
					eW[1]=error-pred;
				//	eW[1]+=(error-pred-eW[1]+(1<<ERROR_SHIFT>>1))>>ERROR_SHIFT;
#endif
				}
				{
					//if(ky==0&&kx==6&&kc==0)//
					//	printf("");

					short
						NW	=rows[1][-1*4+2],
						N	=rows[1][+0*4+2],
						NE	=rows[1][+1*4+2],
						WW	=rows[0][-2*4+2],
						W	=rows[0][-1*4+2],
						*pcurr	=rows[0] +0*4 ;
					int vmax=N, vmin=W; if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE; if(vmax<NE)vmax=NE; int pred=(4*(N+W)+NE-NW+4)>>3;
					//int pred=N+W-NW;
					//if(vmin>NE)vmin=NE; if(vmax<NE)vmax=NE; int pred=W+((5*(N-NW)+NE-WW+4)>>3);
#ifdef ERROR_CORRECTION
					pred-=(3*(eW[2]+eWW[2])+eWWW[2]+(1<<5>>1))>>5; //CLAMP2(pred, -128, 127);
#endif
					CLAMP2(pred, vmin, vmax);
				//	pred+=eW[2]>>2; CLAMP2(pred, -128, 127);
					pred+=offset; CLAMP2(pred, -128, 127);
					int error=0;

					error=planes[2][idx2];
					//error=image[idx+2];

					error+=pred;
					error<<=24;
					error>>=24;
					yuv[2]=error;
					pcurr[2]=yuv[2]-offset;
					
#ifdef ERROR_CORRECTION
					eWWW[2]=eWW[2];
					eWW[2]=eW[2];
					eW[2]=error-pred;
				//	eW[2]+=(error-pred-eW[2]+(1<<ERROR_SHIFT>>1))>>ERROR_SHIFT;
#endif
				}
				image[idx+1]=yuv[0]+128;
				image[idx+2]=yuv[1]+128;
				image[idx+0]=yuv[2]+128;
				guide_check(image, kx, ky);
				rows[0]+=4;
				rows[1]+=4;
			}
		}
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
#if defined _MSC_VER && !defined SILENT
	cycles=__rdtsc()-cycles;
	elapsed=time_sec()-elapsed;
	ptrdiff_t usize=3LL*iw*ih;
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