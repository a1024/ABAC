#include"util.h"
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

//	#include"entropy.h"
typedef unsigned short Emit_t;
ALIGN(32) static unsigned short stats0[3][32], stats1[3][16][32], mixinCDFs[16][16];
int c11_codec(int argc, char **argv)
{
	if(argc!=3)
	{
		printf(
			"Usage: \"%s\"  input  output    Encode/decode.\n"
			, argv[0]
		);
		return 1;
	}
	const char *srcfn=argv[1], *dstfn=argv[2];
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

	for(int k=0;k<16;++k)
		stats0[0][k]=k<<(14-4);
	stats0[0][16]=0x4000;
	memfill(stats0[0]+17, stats0[0]+16, sizeof(short[15]), sizeof(short));
	memfill(stats0[1], stats0[0], sizeof(short[2][32]), sizeof(short[32]));
	memfill(stats1[0], stats0[0], sizeof(stats1), sizeof(stats0));
	for(int ks=0;ks<16;++ks)
	{
		for(int kv=0;kv<16;++kv)
			mixinCDFs[ks][kv]=((0x4000-16)&-(ks<kv))+kv;
	}
	//memset(hist, 0, hsize);
	//for(int k=0;k<257;++k)
	//	CDF[k]=k<<8;
	//for(int ks=0;ks<513;++ks)
	//{
	//	CDF[257+ks]=ks<<7;
	//	CDF[257+513+ks]=ks<<7;
	//}
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
		int errors[3]={0}, syms[6]={0}, cdfs[6]={0}, freqs[6]={0};
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
				curr[4]=errors[0];
				curr[5]=errors[1];
				curr[6]=errors[2];
				if(fwd)
				{
					syms[0]=errors[0]>>4&15;
					syms[1]=errors[1]>>4&15;
					syms[2]=errors[2]>>4&15;
					syms[3]=errors[0]&15;
					syms[4]=errors[1]&15;
					syms[5]=errors[2]&15;
					cdfs[0]=stats0[0][syms[0]];
					cdfs[1]=stats0[1][syms[1]];
					cdfs[2]=stats0[2][syms[2]];
					cdfs[3]=stats1[0][syms[0]][syms[3]];
					cdfs[4]=stats1[1][syms[1]][syms[4]];
					cdfs[5]=stats1[2][syms[2]][syms[5]];
					freqs[0]=stats0[0][syms[0]+1]-cdfs[0];
					freqs[1]=stats0[1][syms[1]+1]-cdfs[1];
					freqs[2]=stats0[2][syms[2]+1]-cdfs[2];
					freqs[3]=stats1[0][syms[0]][syms[3]+1]-cdfs[3];
					freqs[4]=stats1[1][syms[1]][syms[4]+1]-cdfs[4];
					freqs[5]=stats1[2][syms[2]][syms[5]+1]-cdfs[5];
#ifdef ESTIMATE_SIZE
					csizes[0]-=log2((double)freqs[0]/0x4000);
					csizes[1]-=log2((double)freqs[1]/0x4000);
					csizes[2]-=log2((double)freqs[2]/0x4000);
					//if(isinf(csizes[0])||isinf(csizes[1])||isinf(csizes[2]))
					//	LOG_ERROR("");
					csizes[3]-=log2((double)freqs[3]/0x4000);
					csizes[4]-=log2((double)freqs[4]/0x4000);
					csizes[5]-=log2((double)freqs[5]/0x4000);
					//if(isinf(csizes[0])||isinf(csizes[1])||isinf(csizes[2]))
					//	LOG_ERROR("");
#endif
#ifdef __GNUC__
#pragma GCC unroll 6
#endif
					for(int k=0;k<6;++k)
					{
						while(range<0x4000)
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
						low+=range*cdfs[k]>>14;
						range=(range*freqs[k]>>14)-1;
#ifdef AC_VALIDATE
						acval_enc(syms[k], cdfs[k], freqs[k], lo0, lo0+r0, low, low+range, 0, 0);//
#endif
					}
				}
			}
			else
			{
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
				for(int kc=0;kc<3;++kc)
				{
					unsigned cdf, freq;

					while(!(range>>14))
					{
						low=low<<16&0xFFFFFFFFFFFF;
						range=range<<16|0xFFFF;
						code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
						if(range>(low^0xFFFFFFFFFFFF))
							range=low^0xFFFFFFFFFFFF;
					}
					__m256i ms0=_mm256_load_si256((__m256i*)stats0[kc]);
					__m256i mc2=_mm256_set1_epi16((unsigned short)(((code-low)<<14|0x3FFF)/range));
					int mask=_mm256_movemask_epi8(_mm256_cmpgt_epi16(ms0, mc2));
					syms[kc]=FLOOR_LOG2(~mask)>>1;
					cdf=stats0[kc][syms[kc]];
					freq=stats0[kc][syms[kc]+1]-cdf;
#ifdef AC_VALIDATE
					lo0=low; r0=range;
#endif
					low+=range*cdf>>14;
					range=(range*freq>>14)-1;
#ifdef AC_VALIDATE
					acval_dec(syms[kc], cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				}
#ifdef __GNUC__
#pragma GCC unroll 3
#endif
				for(int kc=0;kc<3;++kc)
				{
					unsigned cdf, freq;

					while(!(range>>14))
					{
						low=low<<16&0xFFFFFFFFFFFF;
						range=range<<16|0xFFFF;
						code=(code<<16&0xFFFFFFFFFFFF)|*sptr++;
						if(range>(low^0xFFFFFFFFFFFF))
							range=low^0xFFFFFFFFFFFF;
					}
					__m256i ms0=_mm256_load_si256((__m256i*)stats1[kc][syms[kc]]);
					__m256i mc2=_mm256_set1_epi16((unsigned short)(((code-low)<<14|0x3FFF)/range));
					int mask=_mm256_movemask_epi8(_mm256_cmpgt_epi16(ms0, mc2));
					syms[kc+3]=FLOOR_LOG2(~mask)>>1;
					cdf=stats1[kc][syms[kc]][syms[kc+3]];
					freq=stats1[kc][syms[kc]][syms[kc+3]+1]-cdf;
#ifdef AC_VALIDATE
					lo0=low; r0=range;
#endif
					low+=range*cdf>>14;
					range=(range*freq>>14)-1;
#ifdef AC_VALIDATE
					acval_dec(syms[kc+3], cdf, freq, lo0, lo0+r0, low, low+range, 0, 0, code);
#endif
				}
				yuv[0]=syms[0]<<4|syms[3];
				yuv[1]=syms[1]<<4|syms[4];
				yuv[2]=syms[2]<<4|syms[5];
				yuv[0]+=preds[0]; preds[1]+=yuv[0]; CLAMP2(preds[1], -128, 127);
				yuv[1]+=preds[1]; preds[2]+=yuv[0]; CLAMP2(preds[2], -128, 127);
				yuv[2]+=preds[2];
				curr[0]=yuv[0];
				curr[1]=yuv[1]-yuv[0];
				curr[2]=yuv[2]-yuv[0];
				_mm_store_si128((__m128i*)preds, mp);
				curr[4]=curr[0]-preds[0];
				curr[5]=curr[1]-preds[1];
				curr[6]=curr[2]-preds[2];
				buffer[idx+1]=(unsigned char)(yuv[0]+128);
				buffer[idx+2]=(unsigned char)(yuv[1]+128);
				buffer[idx+0]=(unsigned char)(yuv[2]+128);
			}
			__m256i mixin0=_mm256_load_si256((__m256i*)mixinCDFs[syms[0]]);
			__m256i mixin1=_mm256_load_si256((__m256i*)mixinCDFs[syms[1]]);
			__m256i mixin2=_mm256_load_si256((__m256i*)mixinCDFs[syms[2]]);
			__m256i mixin3=_mm256_load_si256((__m256i*)mixinCDFs[syms[3]]);
			__m256i mixin4=_mm256_load_si256((__m256i*)mixinCDFs[syms[4]]);
			__m256i mixin5=_mm256_load_si256((__m256i*)mixinCDFs[syms[5]]);
			__m256i ms0=_mm256_load_si256((__m256i*)stats0[0]);
			__m256i ms1=_mm256_load_si256((__m256i*)stats0[1]);
			__m256i ms2=_mm256_load_si256((__m256i*)stats0[2]);
			__m256i ms3=_mm256_load_si256((__m256i*)stats1[0][syms[0]]);
			__m256i ms4=_mm256_load_si256((__m256i*)stats1[1][syms[1]]);
			__m256i ms5=_mm256_load_si256((__m256i*)stats1[2][syms[2]]);
			mixin0=_mm256_sub_epi16(mixin0, ms0);
			mixin1=_mm256_sub_epi16(mixin1, ms1);
			mixin2=_mm256_sub_epi16(mixin2, ms2);
			mixin3=_mm256_sub_epi16(mixin3, ms3);
			mixin4=_mm256_sub_epi16(mixin4, ms4);
			mixin5=_mm256_sub_epi16(mixin5, ms5);
			mixin0=_mm256_srai_epi16(mixin0, 7);
			mixin1=_mm256_srai_epi16(mixin1, 7);
			mixin2=_mm256_srai_epi16(mixin2, 7);
			mixin3=_mm256_srai_epi16(mixin3, 8);
			mixin4=_mm256_srai_epi16(mixin4, 8);
			mixin5=_mm256_srai_epi16(mixin5, 8);
			ms0=_mm256_add_epi16(ms0, mixin0);
			ms1=_mm256_add_epi16(ms1, mixin1);
			ms2=_mm256_add_epi16(ms2, mixin2);
			ms3=_mm256_add_epi16(ms3, mixin3);
			ms4=_mm256_add_epi16(ms4, mixin4);
			ms5=_mm256_add_epi16(ms5, mixin5);
			_mm256_store_si256((__m256i*)stats0[0], ms0);
			_mm256_store_si256((__m256i*)stats0[1], ms1);
			_mm256_store_si256((__m256i*)stats0[2], ms2);
			_mm256_store_si256((__m256i*)stats1[0][syms[0]], ms3);
			_mm256_store_si256((__m256i*)stats1[1][syms[1]], ms4);
			_mm256_store_si256((__m256i*)stats1[2][syms[2]], ms5);
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