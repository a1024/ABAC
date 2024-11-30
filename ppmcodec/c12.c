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


//	#define ENABLE_GUIDE

	#define BYPASS_ON_INFLATION
//	#define ENABLE_FILEGUARD	//makes using scripts harder

//	#define ESTIMATE_SIZE

//	#define ENABLE_O2

//	#include"entropy.h"
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
typedef unsigned short Emit_t;
#ifdef ENABLE_O2
ALIGN(32) static unsigned short stats1[3][128][128][256];
#else
ALIGN(32) static unsigned short stats0[3][256];
#endif
int c12_codec(const char *srcfn, const char *dstfn)
{
#ifdef ESTIMATE_SIZE
	double t=time_sec();
	double csizes[6]={0};
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

		guide_save(buffer, iw, ih);//
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
	//int hsize=(int)sizeof(int[257+513+513]);
	//unsigned *hist=(unsigned*)malloc(hsize);
	//unsigned *CDF=(unsigned*)malloc(hsize);
	int sbufsize=iw*16;//rowsize = iw*3, allocate iw*16 just in case
	Emit_t *sbuf=(Emit_t*)malloc(sbufsize);//stream buffer
	if(!pixels||!sbuf)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	memset(sbuf, 0, sbufsize);
	
#ifdef ENABLE_O2
	FILLMEM((unsigned short*)stats1, 0x8000, sizeof(stats1), sizeof(short));
#else
	FILLMEM((unsigned short*)stats0, 0x8000, sizeof(stats0), sizeof(short));
#endif
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
		char yuv[3]={0};
		int errors[3]={0};
		unsigned long long low=0, range=0xFFFFFFFFFFFF, code=0;
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
				code=*sptr++;
				code=*sptr++|code<<16;
				code=*sptr++|code<<16;
			}
			else
				nread=fread(buffer+rowsize*ky, 1, rowsize, fsrc);
		}
		for(int kx=0;kx<iw;++kx, idx+=3)
		{
			short
				*NNE	=rows[2]+1*4*2,
				*NW	=rows[1]-1*4*2,
				*N	=rows[1]+0*4*2,
				*NE	=rows[1]+1*4*2,
				*W	=rows[0]-1*4*2,
				*curr	=rows[0]+0*4*2;
#ifdef ENABLE_O2
			memset(preds, 0, sizeof(preds));
#else
			__m128i mNW	=_mm_load_si128((__m128i*)NW);
			__m128i mN	=_mm_load_si128((__m128i*)N);
			__m128i mW	=_mm_load_si128((__m128i*)W);
			__m128i vmin=_mm_min_epi16(mN, mW);
			__m128i vmax=_mm_max_epi16(mN, mW);
			__m128i mp=_mm_sub_epi16(_mm_add_epi16(mN, mW), mNW);
			mp=_mm_max_epi16(mp, vmin);
			mp=_mm_min_epi16(mp, vmax);
			_mm_store_si128((__m128i*)preds, mp);
#endif
			//if(kx==1128&&ky==154)//
			//	printf("");
#ifdef ENABLE_O2
			unsigned short *stats0[]=
			{
				stats1[0][(N[0]+NE[0]-NNE[0])>>1&127][(N[0]+W[0]-NW[0])>>1&127],
				stats1[1][(N[1]+NE[1]-NNE[1])>>1&127][(N[1]+W[1]-NW[1])>>1&127],
				stats1[2][(N[2]+NE[2]-NNE[2])>>1&127][(N[2]+W[2]-NW[2])>>1&127],
			//	stats1[0][N[4]>>1&127][W[4]>>1&127],
			//	stats1[1][N[5]>>1&127][W[5]>>1&127],
			//	stats1[2][N[6]>>1&127][W[6]>>1&127],
			//	stats1[0][(N[0]-W[0])>>2&127][(N[4]-W[4])>>2&127],
			//	stats1[1][(N[1]-W[1])>>2&127][(N[5]-W[5])>>2&127],
			//	stats1[2][(N[2]-W[2])>>2&127][(N[6]-W[6])>>2&127],
			//	stats1[0][(N[4]+W[4]+1)>>1&255],
			//	stats1[1][(N[5]+W[5]+1)>>1&255],
			//	stats1[2][(N[6]+W[6]+1)>>1&255],
			//	stats1[0][(abs(W[4])+abs(N[4]))&255],
			//	stats1[1][(abs(W[5])+abs(N[5]))&255],
			//	stats1[2][(abs(W[6])+abs(N[6]))&255],
			};
#endif
			
			//if(ky==1&&kx==39)//
			//if(ky==1&&kx==5)//
			//if(ky==1&&kx==743)//
			//if(ky==10&&kx==10)//
			//	printf("");
			if(fwd||!streamsize)
			{
				yuv[0]=buffer[idx+1]-128;//g
				yuv[1]=buffer[idx+2]-128;//b
				yuv[2]=buffer[idx+0]-128;//r
				
				preds[1]+=yuv[0]; CLAMP2(preds[1], -128, 127);
				preds[2]+=yuv[0]; CLAMP2(preds[2], -128, 127);
				errors[0]=yuv[0]-preds[0];
				errors[1]=yuv[1]-preds[1];
				errors[2]=yuv[2]-preds[2];
				curr[0]=yuv[0];
				curr[1]=yuv[1]-yuv[0];
				curr[2]=yuv[2]-yuv[0];
				curr[4]=(char)errors[0];
				curr[5]=(char)errors[1];
				curr[6]=(char)errors[2];
			}
			else
			{
				yuv[0]=0;
				yuv[1]=0;
				yuv[2]=0;
			}
			if(fwd||streamsize)
			{
#ifdef __GNUC__
#pragma GCC unroll 8
#endif
				for(int kb=7, tidx[]={1, 1, 1};kb>=0;--kb)
				{
					unsigned long long r2;
					int bit[3];
					unsigned short *cell0[]=
					{
						stats0[0]+tidx[0],
						stats0[1]+tidx[1],
						stats0[2]+tidx[2],
					};
					unsigned short p0[]=
					{
						*cell0[0],
						*cell0[1],
						*cell0[2],
					};

					if(fwd)
					{
						bit[0]=errors[0]>>kb&1;
						bit[1]=errors[1]>>kb&1;
						bit[2]=errors[2]>>kb&1;
						
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
						for(int kc=0;kc<3;++kc)
						{
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
							r2=range*p0[kc]>>16;
							if(bit[kc])
							{
								low+=r2;
								range-=r2;
							}
							else
								range=r2-1;
#ifdef AC_VALIDATE
							acval_enc(bit[kc], 0, p0[kc], lo0, lo0+r0, low, low+range, 0, 0);
#endif
						}
					}
					else
					{
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
						for(int kc=0;kc<3;++kc)
						{
							unsigned long long mid;

							while(!(range>>16))
							{
								low=low<<16&0xFFFFFFFFFFFF;
								range=range<<16|0xFFFF;
								code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
								if(range>(low^0xFFFFFFFFFFFF))
									range=low^0xFFFFFFFFFFFF;
							}
#ifdef AC_VALIDATE
							lo0=low; r0=range;
#endif
							r2=range*p0[kc]>>16;
							//bit[kc]=code>=low+r2;
							//if(bit[kc])
							//{
							//	low+=r2;
							//	range-=r2;
							//}
							//else
							//	range=r2-1;
							mid=low+r2;
							range-=r2;
							bit[kc]=code>=mid;
							if(bit[kc])
								low=mid;
							else
								range=r2-1;
#ifdef AC_VALIDATE
							acval_dec(bit[kc], 0, p0[kc], lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
						}
						yuv[0]|=bit[0]<<kb;
						yuv[1]|=bit[1]<<kb;
						yuv[2]|=bit[2]<<kb;
					}
					*cell0[0]=p0[0]+(((!bit[0]<<16)-p0[0]+(1<<7>>1))>>7);
					*cell0[1]=p0[1]+(((!bit[1]<<16)-p0[1]+(1<<7>>1))>>7);
					*cell0[2]=p0[2]+(((!bit[2]<<16)-p0[2]+(1<<7>>1))>>7);
					tidx[0]=tidx[0]*2+bit[0];
					tidx[1]=tidx[1]*2+bit[1];
					tidx[2]=tidx[2]*2+bit[2];
				}
			}
			if(!fwd)
			{
				curr[4]=(char)yuv[0];
				curr[5]=(char)yuv[1];
				curr[6]=(char)yuv[2];
				yuv[0]+=preds[0]; preds[1]+=yuv[0]; CLAMP2(preds[1], -128, 127);
				yuv[1]+=preds[1]; preds[2]+=yuv[0]; CLAMP2(preds[2], -128, 127);
				yuv[2]+=preds[2];
				curr[0]=yuv[0];
				curr[1]=yuv[1]-yuv[0];
				curr[2]=yuv[2]-yuv[0];
				buffer[idx+1]=(unsigned char)(yuv[0]+128);
				buffer[idx+2]=(unsigned char)(yuv[1]+128);
				buffer[idx+0]=(unsigned char)(yuv[2]+128);

				guide_check(buffer, kx, ky);//
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
		if(fwd)
		{
			*sptr++=(Emit_t)(low>>32);
			*sptr++=(Emit_t)(low>>16);
			*sptr++=(Emit_t)low;

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
	}
	if(!fwd)
		fwrite(buffer, 1, usize, fdst);
	fclose(fsrc);
	fclose(fdst);
	_mm_free(pixels);
	free(sbuf);
	dstsize=get_filesize(dstfn);
#ifdef ESTIMATE_SIZE
	t=time_sec()-t;
	if(fwd)
	{
		csizes[0]/=8;
		csizes[1]/=8;
		csizes[2]/=8;
		csizes[3]/=8;
		csizes[4]/=8;
		csizes[5]/=8;
		printf("%12.2lf +\n%12.2lf +\n%12.2lf +\n%12.2lf +\n%12.2lf +\n%12.2lf =\n%12.2lf ->\n%9td\n",
			csizes[0],
			csizes[1],
			csizes[2],
			csizes[3],
			csizes[4],
			csizes[5],
			csizes[0]+csizes[1]+csizes[2]+csizes[3]+csizes[4]+csizes[5],
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