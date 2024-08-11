#include"codec.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
//#include"blist.h"
static const char file[]=__FILE__;


	#define BYPASS_ON_INFLATION
//	#define ENABLE_FILEGUARD	//makes using scripts harder

//	#define ESTIMATE_SIZE
	#define USE_RCTJ2K
	#define ANS_INTERLEAVE3


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
int c09_codec(const char *srcfn, const char *dstfn)
{
#ifdef ESTIMATE_SIZE
	double t=time_sec();
	double csizes[3]={0};
	size_t nqueries=0;
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
	unsigned short *sbuf=(unsigned short*)malloc(sbufsize);//stream buffer
	if(!pixels||!hist||!CDF||!sbuf)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
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
	if(fwd)//encode
	{
#ifdef ANS_INTERLEAVE3
		int statbufsize=iw*sizeof(int[8]);
		unsigned *statbuf=(unsigned*)_mm_malloc(statbufsize, sizeof(__m128i));
#else
		unsigned *statbuf=(unsigned*)malloc(rowsize*sizeof(int[2]));
#endif
		if(!statbuf)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
#ifdef ANS_INTERLEAVE3
		memset(statbuf, 0, statbufsize);
#endif
		//BList list;
		//blist_init(&list);
		for(int ky=0, idx=0;ky<ih;++ky)
		{
			ALIGN(32) short *rows[]=
			{
				pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
			};
			unsigned short *sptr=sbuf, *sbufend=sbuf+sbufsize;
			(void)sbufend;
			for(int kx=0;kx<iw;++kx, idx+=3)
			{
				short
					*NW	=rows[1]-1*4*2,
					*N	=rows[1]+0*4*2,
					*W	=rows[0]-1*4*2,
					*curr	=rows[0]+0*4*2;
				ALIGN(16) short preds[8];
				__m128i mNW	=_mm_load_si128((__m128i*)NW);
				__m128i mN	=_mm_load_si128((__m128i*)N);
				__m128i mW	=_mm_load_si128((__m128i*)W);
				__m128i vmin=_mm_min_epi16(mN, mW);
				__m128i vmax=_mm_max_epi16(mN, mW);
				__m128i mp=_mm_sub_epi16(_mm_add_epi16(mN, mW), mNW);
				mp=_mm_max_epi16(mp, vmin);
				mp=_mm_min_epi16(mp, vmax);
				_mm_store_si128((__m128i*)preds, mp);

				int yuv[]=
				{
					buffer[idx+1]-128,//g
					buffer[idx+2]-128,//b
					buffer[idx+0]-128,//r
				};
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

				int errors[]=
				{
					yuv[0]-preds[0],
					yuv[1]-preds[1],
					yuv[2]-preds[2],
				};
				curr[4]=errors[0];
				curr[5]=errors[1];
				curr[6]=errors[2];
				errors[0]=errors[0]<<(32-8)>>(32-8);
				errors[1]=errors[1]<<(32-9)>>(32-9);
				errors[2]=errors[2]<<(32-9)>>(32-9);
				errors[0]=errors[0]<<1^errors[0]>>31;//pack sign
				errors[1]=errors[1]<<1^errors[1]>>31;
				errors[2]=errors[2]<<1^errors[2]>>31;
				
				unsigned *curr_CDF;
				
#ifdef ANS_INTERLEAVE3
				curr_CDF=CDF;
				statbuf[kx*8+0]=curr_CDF[errors[0]];
				statbuf[kx*8+4]=curr_CDF[errors[0]+1]-curr_CDF[errors[0]];
				curr_CDF=CDF+257;
				statbuf[kx*8+1]=curr_CDF[errors[1]];
				statbuf[kx*8+5]=curr_CDF[errors[1]+1]-curr_CDF[errors[1]];
				curr_CDF=CDF+257+513;
				statbuf[kx*8+2]=curr_CDF[errors[2]];
				statbuf[kx*8+6]=curr_CDF[errors[2]+1]-curr_CDF[errors[2]];
#else
				curr_CDF=CDF;
				statbuf[kx*6+0]=curr_CDF[errors[0]];
				statbuf[kx*6+1]=curr_CDF[errors[0]+1]-curr_CDF[errors[0]];
				curr_CDF=CDF+257;
				statbuf[kx*6+2]=curr_CDF[errors[1]];
				statbuf[kx*6+3]=curr_CDF[errors[1]+1]-curr_CDF[errors[1]];
				curr_CDF=CDF+257+513;
				statbuf[kx*6+4]=curr_CDF[errors[2]];
				statbuf[kx*6+5]=curr_CDF[errors[2]+1]-curr_CDF[errors[2]];
#endif

				++hist[errors[0]];
				++hist[256];
				++hist[257+errors[1]];
				++hist[257+512];
				++hist[257+513+errors[2]];
				++hist[257+513+512];
#ifdef ESTIMATE_SIZE
				csizes[0]-=log2((double)statbuf[kx*6+1]/0x10000);
				csizes[1]-=log2((double)statbuf[kx*6+3]/0x10000);
				csizes[2]-=log2((double)statbuf[kx*6+5]/0x10000);
#endif
				rows[0]+=4*2;
				rows[1]+=4*2;
				rows[2]+=4*2;
				rows[3]+=4*2;
			}
#ifdef ANS_INTERLEAVE3
			ALIGN(16) unsigned states[]={0x10000, 0x10000, 0x10000, 0x10000};
			for(int kx=iw-1;kx>=0;--kx)
			{
				if(!ky&&kx==1)//
					printf("");
				__m128i mstate=_mm_load_si128((__m128i*)states);
				__m128i mfreq=_mm_load_si128((__m128i*)(statbuf+kx*8+4));
				int mask=_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(mfreq, _mm_srli_epi32(mstate, 16))));
				if(!(mask&1))
				{
					*sptr++=(unsigned short)states[0];
					states[0]>>=16;
				}
				if(!(mask&2))
				{
					*sptr++=(unsigned short)states[1];
					states[1]>>=16;
				}
				if(!(mask&4))
				{
					*sptr++=(unsigned short)states[2];
					states[2]>>=16;
				}
				mstate=_mm_load_si128((__m128i*)states);
				__m256d dstate=_mm256_cvtepi32_pd(mstate);
				__m256d dfreq=_mm256_cvtepi32_pd(mfreq);
				dstate=_mm256_div_pd(dstate, dfreq);
				__m128i q=_mm256_cvttpd_epi32(dstate);
				__m128i r=_mm_mullo_epi32(q, mfreq);
				__m128i mcdf=_mm_load_si128((__m128i*)(statbuf+kx*8+0));
				r=_mm_sub_epi32(mstate, r);
				r=_mm_add_epi32(r, mcdf);
				q=_mm_slli_epi32(q, 16);
				q=_mm_or_si128(q, r);
				_mm_store_si128((__m128i*)states, q);
			}
			*sptr++=(unsigned short)states[0];
			*sptr++=(unsigned short)(states[0]>>16);
			*sptr++=(unsigned short)states[1];
			*sptr++=(unsigned short)(states[1]>>16);
			*sptr++=(unsigned short)states[2];
			*sptr++=(unsigned short)(states[2]>>16);
#else
			unsigned state=0x10000;
			for(int kx=rowsize-1;kx>=0;--kx)
			{
				unsigned cdf=statbuf[kx<<1|0], freq=statbuf[kx<<1|1];
				if((state>>16)>=(unsigned)freq)//renorm
				{
					*sptr++=(unsigned short)state;
					state>>=16;
				}
				state=state/freq<<16|(cdf+state%freq);//update
			}
			*sptr++=(unsigned short)state;//ANS is little-endian
			*sptr++=(unsigned short)(state>>16);
#endif

			int streamsize=(int)(sptr-sbuf);//number of emits
			if(streamsize<<1>rowsize+4)//bypass row
			{
				streamsize=0;
				fwrite(&streamsize, 1, 4, fdst);
				fwrite(buffer+rowsize*ky, 1, rowsize, fdst);
			}
			else
			{
				fwrite(&streamsize, 1, 4, fdst);
				fwrite(sbuf, sizeof(short), streamsize, fdst);
			}

			update_CDF(hist, 256, CDF);
			update_CDF(hist+257, 512, CDF+257);
			update_CDF(hist+257+513, 512, CDF+257+513);
		}
#ifdef ANS_INTERLEAVE3
		_mm_free(statbuf);
#else
		free(statbuf);
#endif
	}
	else//decode
	{
		//unsigned short *CDF2sym=(unsigned short*)malloc(sizeof(short[3<<16]));
		for(int ky=0, idx=0;ky<ih;++ky)
		{
			ALIGN(32) short *rows[]=
			{
				pixels+((iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
				pixels+((iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
			};
			int streamsize=0;
			nread=fread(&streamsize, 1, 4, fsrc);
			if(!streamsize)//bypass decode row
			{
				fread(buffer+rowsize*ky, 1, rowsize, fsrc);
				for(int kx=0;kx<iw;++kx, idx+=3)//dummy process to set neighbors & stats
				{
					short
						*NW	=rows[1]-1*4*2,
						*N	=rows[1]+0*4*2,
						*W	=rows[0]-1*4*2,
						*curr	=rows[0]+0*4*2;
					ALIGN(16) short preds[8];
					__m128i mNW	=_mm_load_si128((__m128i*)NW);
					__m128i mN	=_mm_load_si128((__m128i*)N);
					__m128i mW	=_mm_load_si128((__m128i*)W);
					__m128i vmin=_mm_min_epi16(mN, mW);
					__m128i vmax=_mm_max_epi16(mN, mW);
					__m128i mp=_mm_sub_epi16(_mm_add_epi16(mN, mW), mNW);
					mp=_mm_max_epi16(mp, vmin);
					mp=_mm_min_epi16(mp, vmax);
					_mm_store_si128((__m128i*)preds, mp);

					int yuv[]=
					{
						buffer[idx+1]-128,//g
						buffer[idx+2]-128,//b
						buffer[idx+0]-128,//r
					};
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

					int errors[]=
					{
						yuv[0]-preds[0],
						yuv[1]-preds[1],
						yuv[2]-preds[2],
					};
					curr[4]=errors[0];
					curr[5]=errors[1];
					curr[6]=errors[2];
					errors[0]=errors[0]<<(32-8)>>(32-8);
					errors[1]=errors[1]<<(32-9)>>(32-9);
					errors[2]=errors[2]<<(32-9)>>(32-9);
					errors[0]=errors[0]<<1^errors[0]>>31;//pack sign
					errors[1]=errors[1]<<1^errors[1]>>31;
					errors[2]=errors[2]<<1^errors[2]>>31;

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
			}
			else//normal decode row
			{
				fread(sbuf, sizeof(short), streamsize, fsrc);
				unsigned short *sptr=sbuf+streamsize-1;
#ifdef ANS_INTERLEAVE3
				ALIGN(16) unsigned states[4]={0}, freqs[4]={0}, cdfs[4]={0};
				states[2]=*sptr--;
				states[2]=*sptr--|states[2]<<16;
				states[1]=*sptr--;
				states[1]=*sptr--|states[1]<<16;
				states[0]=*sptr--;
				states[0]=*sptr--|states[0]<<16;
#else
				unsigned state=0;
				state=*sptr--;
				state=*sptr--|state<<16;
#endif
				for(int kx=0;kx<iw;++kx, idx+=3)
				{
					short
						*NW	=rows[1]-1*4*2,
						*N	=rows[1]+0*4*2,
						*W	=rows[0]-1*4*2,
						*curr	=rows[0]+0*4*2;

					ALIGN(16) short preds[8];
					__m128i mNW	=_mm_load_si128((__m128i*)NW);
					__m128i mN	=_mm_load_si128((__m128i*)N);
					__m128i mW	=_mm_load_si128((__m128i*)W);
					__m128i vmin=_mm_min_epi16(mN, mW);
					__m128i vmax=_mm_max_epi16(mN, mW);
					__m128i mp=_mm_sub_epi16(_mm_add_epi16(mN, mW), mNW);
					mp=_mm_max_epi16(mp, vmin);
					mp=_mm_min_epi16(mp, vmax);
					_mm_store_si128((__m128i*)preds, mp);

					ALIGN(16) int errors[3];
#ifdef ANS_INTERLEAVE3
					int cdf, freq, sym;
					unsigned code, *curr_CDF;

					curr_CDF=CDF;
					code=(unsigned short)states[0];
					sym=0;
					for(freq=0;;)
					{
						cdf=freq;
						freq=curr_CDF[sym+1];
						if((unsigned)freq>code)
							break;
						++sym;
					}
					freq-=cdf;
					errors[0]=sym;
					freqs[0]=freq;
					cdfs[0]=cdf;
					
					curr_CDF=CDF+257;
					code=(unsigned short)states[1];
					sym=0;
					for(freq=0;;)
					{
						cdf=freq;
						freq=curr_CDF[sym+1];
						if((unsigned)freq>code)
							break;
						++sym;
					}
					freq-=cdf;
					errors[1]=sym;
					freqs[1]=freq;
					cdfs[1]=cdf;
					
					curr_CDF=CDF+257+513;
					code=(unsigned short)states[2];
					sym=0;
					for(freq=0;;)
					{
						cdf=freq;
						freq=curr_CDF[sym+1];
						if((unsigned)freq>code)
							break;
						++sym;
					}
					freq-=cdf;
					errors[2]=sym;
					freqs[2]=freq;
					cdfs[2]=cdf;

					__m128i mstate=_mm_load_si128((__m128i*)states);
					__m128i mfreq=_mm_load_si128((__m128i*)freqs);
					__m128i mcdf=_mm_load_si128((__m128i*)cdfs);
					__m128i mlo=_mm_and_si128(mstate, _mm_set1_epi32(0xFFFF));
					__m128i mhi=_mm_srli_epi32(mstate, 16);
					mhi=_mm_mullo_epi32(mhi, mfreq);
					mhi=_mm_add_epi32(mhi, mlo);
					mhi=_mm_sub_epi32(mhi, mcdf);
					int mask=_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(_mm_set1_epi32(0x10000), mhi)));
					_mm_store_si128((__m128i*)states, mhi);
					if(mask&4)
						states[2]=*sptr--|states[2]<<16;
					if(mask&2)
						states[1]=*sptr--|states[1]<<16;
					if(mask&1)
						states[0]=*sptr--|states[0]<<16;
#else
					unsigned code;
					int cdf, freq, sym;
					unsigned *curr_CDF;

					curr_CDF=CDF;
					code=(unsigned short)state;
					sym=0;

					for(freq=0;;)
					{
						cdf=freq;
						freq=curr_CDF[sym+1];
						if((unsigned)freq>code)
							break;
						++sym;
					}
					freq-=cdf;
					//for(;curr_CDF[sym+1]<=code;++sym);
					//cdf=curr_CDF[sym], freq=curr_CDF[sym+1]-cdf;
					state=freq*(state>>16)+code-cdf;//update
					if(state<0x10000)//renorm
						state=*sptr--|state<<16;
					errors[0]=sym;

					curr_CDF=CDF+257;
					code=(unsigned short)state;
					sym=0;

					for(freq=0;;)
					{
						cdf=freq;
						freq=curr_CDF[sym+1];
						if((unsigned)freq>code)
							break;
						++sym;
					}
					freq-=cdf;
					//for(;curr_CDF[sym+1]<=code;++sym);
					//cdf=curr_CDF[sym], freq=curr_CDF[sym+1]-cdf;
					state=freq*(state>>16)+code-cdf;//update
					if(state<0x10000)//renorm
						state=*sptr--|state<<16;
					errors[1]=sym;

					curr_CDF=CDF+257+513;
					code=(unsigned short)state;
					sym=0;

					for(freq=0;;)
					{
						cdf=freq;
						freq=curr_CDF[sym+1];
						if((unsigned)freq>code)
							break;
						++sym;
					}
					freq-=cdf;
					//for(;curr_CDF[sym+1]<=code;++sym);
					//cdf=curr_CDF[sym], freq=curr_CDF[sym+1]-cdf;
					state=freq*(state>>16)+code-cdf;//update
					if(state<0x10000)//renorm
						state=*sptr--|state<<16;
					errors[2]=sym;
#endif

					++hist[errors[0]];
					++hist[256];
					++hist[257+errors[1]];
					++hist[257+512];
					++hist[257+513+errors[2]];
					++hist[257+513+512];

					errors[0]=errors[0]>>1^-(errors[0]&1);//unpack sign
					errors[1]=errors[1]>>1^-(errors[1]&1);
					errors[2]=errors[2]>>1^-(errors[2]&1);
					int yuv[]=
					{
						errors[0]+preds[0],
						errors[1]+preds[1],
						errors[2]+preds[2],
					};
					yuv[0]=yuv[0]<<(32-8)>>(32-8);
					yuv[1]=yuv[1]<<(32-9)>>(32-9);
					yuv[2]=yuv[2]<<(32-9)>>(32-9);
					curr[0]=yuv[0];
					curr[1]=yuv[1];
					curr[2]=yuv[2];
					curr[4]=errors[0]-preds[0];
					curr[5]=errors[1]-preds[1];
					curr[6]=errors[2]-preds[2];
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

					rows[0]+=4*2;
					rows[1]+=4*2;
					rows[2]+=4*2;
					rows[3]+=4*2;
				}
#ifdef ANS_INTERLEAVE3
				if(states[0]!=0x10000||states[1]!=0x10000||states[2]!=0x10000)
				{
					LOG_ERROR("Decode error  Y %d", ky);
					return 1;
				}
#else
				if(state!=0x10000)
				{
					LOG_ERROR("Decode error  Y %d", ky);
					return 1;
				}
#endif
			}
			update_CDF(hist, 256, CDF);
			update_CDF(hist+257, 512, CDF+257);
			update_CDF(hist+257+513, 512, CDF+257+513);
		}
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
		printf("%.2lf  ->  %.2lf + %.2lf + %.2lf = %.2lf  ->  %td\n",
			nqueries/8.,
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