#include"fast.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;


//	#define ENABLE_GUIDE
//	#define SIMD_FLOOR_LOG2		//X
//	#define PROFILE_CSIZE
//	#define SAVE_K


#include"ac.h"
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
static void calc_ctx(const __m256i *N, const __m256i *W, const __m256i *NW, const __m256i *eW, const __m256i *mminmag, __m256i *pred, __m256i *mag)
{
#ifdef SIMD_FLOOR_LOG2
	__m256i hmasks[]=
	{
		_mm256_set1_epi32(0x55555555),
		_mm256_set1_epi32(0x33333333),
		_mm256_set1_epi32(0x0F0F0F0F),
		_mm256_set1_epi32(0x00FF00FF),
	};
#endif
	//pred = median3(N, W, N+W-NW)
	//mag = eW
	//mag = FLOOR_LOG2_P1(PACK_SIGN(mag)+1)

	*mag=_mm256_max_epi32(*eW, *mminmag);
	*mag=_mm256_add_epi32(*mag, _mm256_set1_epi32(1));
#ifdef SIMD_FLOOR_LOG2
	*mag=_mm256_or_si256(*mag, _mm256_srli_epi32(*mag, 1));
	*mag=_mm256_or_si256(*mag, _mm256_srli_epi32(*mag, 2));
	*mag=_mm256_or_si256(*mag, _mm256_srli_epi32(*mag, 4));
	*mag=_mm256_or_si256(*mag, _mm256_srli_epi32(*mag, 8));
#endif
	__m256i mmin=_mm256_min_epi32(*N, *W);
	__m256i mmax=_mm256_max_epi32(*N, *W);
	*pred=_mm256_add_epi32(*N, *W);
	*pred=_mm256_sub_epi32(*pred, *NW);
	*pred=_mm256_max_epi32(*pred, mmin);
	*pred=_mm256_min_epi32(*pred, mmax);
#ifdef SIMD_FLOOR_LOG2
	*mag=_mm256_sub_epi32(*mag, _mm256_and_si256(_mm256_srli_epi32(*mag, 1), hmasks[0]));//hamming weight
	*mag=_mm256_add_epi32(_mm256_and_si256(*mag, hmasks[1]), _mm256_and_si256(_mm256_srli_epi32(*mag, 2), hmasks[1]));
	*mag=_mm256_add_epi32(_mm256_and_si256(*mag, hmasks[2]), _mm256_and_si256(_mm256_srli_epi32(*mag, 4), hmasks[2]));
	*mag=_mm256_add_epi32(_mm256_and_si256(*mag, hmasks[3]), _mm256_and_si256(_mm256_srli_epi32(*mag, 8), hmasks[3]));
#endif
}
static void pack_sign(__m256i *data)
{
	//x = curr-pred
	//x = x<<1^-(x<0)
	__m256i negmask=_mm256_cmpgt_epi32(_mm256_setzero_si256(), *data);
	*data=_mm256_slli_epi32(*data, 1);
	*data=_mm256_xor_si256(*data, negmask);
}
static void unpack_sign(__m256i *curr)
{
	//curr = (curr>>1^-(curr&1)) + pred
	__m256i negmask=_mm256_and_si256(*curr, _mm256_set1_epi32(1));
	negmask=_mm256_cmpeq_epi32(negmask, _mm256_set1_epi32(1));
	*curr=_mm256_srai_epi32(*curr, 1);
	*curr=_mm256_xor_si256(*curr, negmask);
}
static void calc_update(const __m256i *eW, const __m256i *eNEEE, __m256i *ecurr)//writes to hi64 bits
{
	//x = (2*eW+eNEEE+ecurr)>>2
	*ecurr=_mm256_add_epi32(*ecurr, _mm256_slli_epi32(*eW, 1));
	*ecurr=_mm256_add_epi32(*ecurr, *eNEEE);
	*ecurr=_mm256_srai_epi32(*ecurr, 2);//doesn't matter logic or arithmetic shift
}
int f22_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;

	size_t ebufsize=sizeof(int[4*4*2*2])*(image->iw+16LL);//4 padded rows * 4 channels max * {pixels, errors} * 2 threads (interleaved)
	int *pixels=(int*)_mm_malloc(ebufsize, sizeof(__m256i));
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	int minmag=1<<image->depth>>7;
	pixels[ 0]=0;
	pixels[ 1]=0;
	pixels[ 2]=0;
	pixels[ 3]=0;
	pixels[ 4]=0;
	pixels[ 5]=0;
	pixels[ 6]=0;
	pixels[ 7]=0;
	pixels[ 8]=minmag;
	pixels[ 9]=minmag;
	pixels[10]=minmag;
	pixels[11]=minmag;
	pixels[12]=minmag;
	pixels[13]=minmag;
	pixels[14]=minmag;
	pixels[15]=minmag;
	memfill(pixels+16, pixels, ebufsize-sizeof(int[16]), sizeof(int[16]));
	__m256i mminmag=_mm256_set1_epi32(minmag);
	if(image->depth==8)
	{
		mminmag=_mm256_srli_epi32(mminmag, 2);
		minmag>>=2;
	}
	ptrdiff_t usize=(ptrdiff_t)image->iw*image->ih*image->nch*image->depth>>3;
	if(fwd)
	{
#ifdef ENABLE_GUIDE
		guide=image;
#endif
#ifdef PROFILE_CSIZE
		ptrdiff_t csizes[4*3]={0};//{unary, bypass}
#endif
#ifdef SAVE_K
		int half=1<<image->depth>>1;
		Image kimage={0};
		if(loud)
		{
			image_copy(&kimage, image);
			if(!kimage.data)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
		}
#endif
		GolombRiceCoder ec0, ec1;
		DList list0, list1;
		dlist_init(&list0, 1, 1024, 0);
		dlist_init(&list1, 1, 1024, 0);
		gr_enc_init(&ec0, &list0);
		gr_enc_init(&ec1, &list1);
		ALIGN(16) short c2[8]={0};
		int floorhalfheight=image->ih>>1, ceilhalfheight=image->ih-floorhalfheight;
		for(int ky=0, idx0=0, idx1=floorhalfheight*image->iw*image->nch;ky<ceilhalfheight;++ky)
		{
			ALIGN(32) int *rows[]=
			{
				pixels+(((image->iw+16LL)*((ky-0LL)&3)+8)<<4),
				pixels+(((image->iw+16LL)*((ky-1LL)&3)+8)<<4),
				pixels+(((image->iw+16LL)*((ky-2LL)&3)+8)<<4),
				pixels+(((image->iw+16LL)*((ky-3LL)&3)+8)<<4),
			};
			for(int kx=0;kx<image->iw;++kx, idx0+=image->nch, idx1+=image->nch)
			{
				__m256i
					NW	=_mm256_load_si256((__m256i*)rows[1+0]-1*2+0),
					N	=_mm256_load_si256((__m256i*)rows[1+0]+0*2+0),
					W	=_mm256_load_si256((__m256i*)rows[0+0]-1*2+0),
					eNEEE	=_mm256_load_si256((__m256i*)rows[1+0]+3*2+1),
					eW	=_mm256_load_si256((__m256i*)rows[0+0]-1*2+1);
				int	*curr	=rows[0+0]+0*16;
				__m256i mcurr, mpred, mmag;
				calc_ctx(&N, &W, &NW, &eW, &mminmag, &mpred, &mmag);
				//if(ky==1&&kx==5)//
				//if(idx==1170)//
				//if(idx==6936)//
				//if(idx==4722)//
				//if(idx==4731)//
				//	printf("");
				ALIGN(32) int mag[8];
				_mm256_store_si256((__m256i*)mag, mmag);
				switch(image->nch)
				{
				case 4:
					c2[3+0]=image->data[idx0+3];
					c2[3+4]=image->data[idx1+3];
				case 3:
					c2[2+0]=image->data[idx0+2];
					c2[2+4]=image->data[idx1+2];
				case 2:
					c2[1+0]=image->data[idx0+1];
					c2[1+4]=image->data[idx1+1];
				case 1:
					c2[0+0]=image->data[idx0+0];
					c2[0+4]=image->data[idx1+0];
					break;
				}
				if(image->nch>=3)
				{
					//JPEG2000 RCT
					c2[0+0]-=c2[1+0];
					c2[0+4]-=c2[1+4];
					c2[2+0]-=c2[1+0];
					c2[2+4]-=c2[1+4];
					c2[1+0]+=(c2[0+0]+c2[2+0])>>2;
					c2[1+4]+=(c2[0+4]+c2[2+4])>>2;
				}
				mcurr=_mm256_cvtepi16_epi32(_mm_load_si128((__m128i*)c2));
				__m256i delta=_mm256_sub_epi32(mcurr, mpred);
				pack_sign(&delta);
				_mm256_store_si256((__m256i*)curr+0, mcurr);
				_mm256_store_si256((__m256i*)curr+1, delta);
#ifdef SIMD_FLOOR_LOG2
				switch(image->nch)
				{
				case 4:
					gr_enc_POT(&ec0, curr[8+0+3], mag[3+0]);
					gr_enc_POT(&ec1, curr[8+4+3], mag[3+4]);
				case 3:
					gr_enc_POT(&ec0, curr[8+0+2], mag[2+0]);
					gr_enc_POT(&ec1, curr[8+4+2], mag[2+4]);
				case 2:
					gr_enc_POT(&ec0, curr[8+0+1], mag[1+0]);
					gr_enc_POT(&ec1, curr[8+4+1], mag[1+4]);
				case 1:
					gr_enc_POT(&ec0, curr[8+0+0], mag[0+0]);
					gr_enc_POT(&ec1, curr[8+4+0], mag[0+4]);
					break;
				}
#else
				switch(image->nch)
				{
				case 4:
					gr_enc_POT(&ec0, curr[8+0+3], FLOOR_LOG2(mag[3+0]));
					gr_enc_POT(&ec1, curr[8+4+3], FLOOR_LOG2(mag[3+4]));
				case 3:
					gr_enc_POT(&ec0, curr[8+0+2], FLOOR_LOG2(mag[2+0]));
					gr_enc_POT(&ec1, curr[8+4+2], FLOOR_LOG2(mag[2+4]));
				case 2:
					gr_enc_POT(&ec0, curr[8+0+1], FLOOR_LOG2(mag[1+0]));
					gr_enc_POT(&ec1, curr[8+4+1], FLOOR_LOG2(mag[1+4]));
				case 1:
					gr_enc_POT(&ec0, curr[8+0+0], FLOOR_LOG2(mag[0+0]));
					gr_enc_POT(&ec1, curr[8+4+0], FLOOR_LOG2(mag[0+4]));
					break;
				}
#endif
#ifdef SAVE_K
				if(loud)
				{
					int temp;
					temp=1<<FLOOR_LOG2(mag[0]+1), kimage.data[idx+0]=CLAMP(-half, temp, half-1)-half;
					temp=1<<FLOOR_LOG2(mag[1]+1), kimage.data[idx+1]=CLAMP(-half, temp, half-1)-half;
					temp=1<<FLOOR_LOG2(mag[2]+1), kimage.data[idx+2]=CLAMP(-half, temp, half-1)-half;
				}
#endif
#ifdef PROFILE_CSIZE
				for(int kc=image->nch-1;kc>=0;--kc)
				{
					int sym=curr[kc+4];
					int magnitude=mag[kc]+1;
					csizes[0<<2|kc]+=sym/magnitude;
					++csizes[1<<2|kc];
					int bypass=sym%magnitude;
					int nbypass=floor_log2_32(magnitude)+1;
					csizes[2<<2|kc]+=nbypass-(bypass<(1<<nbypass)-magnitude);
				}
#endif
				calc_update(&eW, &eNEEE, &delta);
				_mm256_store_si256((__m256i*)curr+1, delta);
				
				__m256i step=_mm256_set1_epi64x(sizeof(int[8*2]));
				__m256i mrows=_mm256_load_si256((__m256i*)rows);
				mrows=_mm256_add_epi64(mrows, step);
				_mm256_store_si256((__m256i*)rows, mrows);
			}
		}
		gr_enc_flush(&ec0);
		gr_enc_flush(&ec1);
		array_append(data, &list0.nobj, 1, 4, 1, 0, 0);
		dlist_appendtoarray(&list0, data);
		dlist_appendtoarray(&list1, data);
		if(loud)
		{
#ifdef SAVE_K
			image_snapshot8(&kimage);
			image_clear(&kimage);
#endif
			ptrdiff_t csize=list0.nobj+list1.nobj;
			t0=time_sec()-t0;
			printf("%14td/%14td = %10.6lf%%  CR %lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
#ifdef PROFILE_CSIZE
			printf("C         unary         stop_bit          bypass\n");
			double total_unary=0, total_stop=0, total_bypass=0;
			for(int kc=0;kc<image->nch;++kc)
			{
				double size_unary=csizes[0<<2|kc]/8., size_stop=csizes[1<<2|kc]/8., size_bypass=csizes[2<<2|kc]/8.;
				printf("%d %16.2lf %16.2lf %16.2lf\n", kc, size_unary, size_stop, size_bypass);
				total_unary+=size_unary;
				total_stop+=size_stop;
				total_bypass+=size_bypass;
			}
			printf("T %16.2lf %16.2lf %16.2lf\n", total_unary, total_stop, total_bypass);
#endif
			printf("E  %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
		dlist_clear(&list0);
		dlist_clear(&list1);
	}
	else
	{
		GolombRiceCoder ec0, ec1;
		unsigned halfsize;
		memcpy(&halfsize, cbuf, 4);
		cbuf+=4;
		clen-=4;
		gr_dec_init(&ec0, cbuf, cbuf+halfsize);
		gr_dec_init(&ec1, cbuf+halfsize, cbuf+clen);
		__m256i packpixels=_mm256_set_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1, 13, 12,  9,  8,  5,  4,  1,  0,
			-1, -1, -1, -1, -1, -1, -1, -1, 13, 12,  9,  8,  5,  4,  1,  0
			//15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
		);
		int floorhalfheight=image->ih>>1, ceilhalfheight=image->ih-floorhalfheight;
		for(int ky=0, idx0=0, idx1=floorhalfheight*image->iw*image->nch;ky<ceilhalfheight;++ky)
		{
			ALIGN(32) int *rows[]=
			{
				pixels+(((image->iw+16LL)*((ky-0LL)&3)+8)<<4),
				pixels+(((image->iw+16LL)*((ky-1LL)&3)+8)<<4),
				pixels+(((image->iw+16LL)*((ky-2LL)&3)+8)<<4),
				pixels+(((image->iw+16LL)*((ky-3LL)&3)+8)<<4),
			};
			for(int kx=0;kx<image->iw;++kx, idx0+=image->nch, idx1+=image->nch)
			{
				__m256i
					NW	=_mm256_load_si256((__m256i*)rows[1+0]-1*2+0),
					N	=_mm256_load_si256((__m256i*)rows[1+0]+0*2+0),
					W	=_mm256_load_si256((__m256i*)rows[0+0]-1*2+0),
					eNEEE	=_mm256_load_si256((__m256i*)rows[1+0]+3*2+1),
					eW	=_mm256_load_si256((__m256i*)rows[0+0]-1*2+1);
				int	*curr	=rows[0+0]+0*16;
				__m256i mcurr, mpred, mmag;
				calc_ctx(&N, &W, &NW, &eW, &mminmag, &mpred, &mmag);
				//if(ky==1&&kx==5)//
				//if(idx==1161)//
				//if(idx==6936)//
				//if(idx==4722)//
				//	printf("");
				ALIGN(32) int mag[8];
				_mm256_store_si256((__m256i*)mag, mmag);
#ifdef SIMD_FLOOR_LOG2
				switch(image->nch)
				{
				case 4:
					curr[8+0+3]=gr_dec_POT(&ec0, mag[3+0]);
					curr[8+4+3]=gr_dec_POT(&ec1, mag[3+4]);
				case 3:
					curr[8+0+2]=gr_dec_POT(&ec0, mag[2+0]);
					curr[8+4+2]=gr_dec_POT(&ec1, mag[2+4]);
				case 2:
					curr[8+0+1]=gr_dec_POT(&ec0, mag[1+0]);
					curr[8+4+1]=gr_dec_POT(&ec1, mag[1+4]);
				case 1:
					curr[8+0+0]=gr_dec_POT(&ec0, mag[0+0]);
					curr[8+4+0]=gr_dec_POT(&ec1, mag[0+4]);
					break;
				}
#else
				switch(image->nch)
				{
				case 4:
					curr[8+0+3]=gr_dec_POT(&ec0, FLOOR_LOG2(mag[3+0]));
					curr[8+4+3]=gr_dec_POT(&ec1, FLOOR_LOG2(mag[3+4]));
				case 3:
					curr[8+0+2]=gr_dec_POT(&ec0, FLOOR_LOG2(mag[2+0]));
					curr[8+4+2]=gr_dec_POT(&ec1, FLOOR_LOG2(mag[2+4]));
				case 2:
					curr[8+0+1]=gr_dec_POT(&ec0, FLOOR_LOG2(mag[1+0]));
					curr[8+4+1]=gr_dec_POT(&ec1, FLOOR_LOG2(mag[1+4]));
				case 1:
					curr[8+0+0]=gr_dec_POT(&ec0, FLOOR_LOG2(mag[0+0]));
					curr[8+4+0]=gr_dec_POT(&ec1, FLOOR_LOG2(mag[0+4]));
					break;
				}
#endif
				__m256i delta=_mm256_load_si256((__m256i*)curr+1);
				mcurr=delta;
				unpack_sign(&mcurr);
				calc_update(&eW, &eNEEE, &delta);
				mcurr=_mm256_add_epi32(mcurr, mpred);
				_mm256_store_si256((__m256i*)curr+0, mcurr);
				_mm256_store_si256((__m256i*)curr+1, delta);

				mcurr=_mm256_shuffle_epi8(mcurr, packpixels);
				short *c2_0=image->data+idx0;
				short *c2_1=image->data+idx1;
				ALIGN(16) short c3[16];
				_mm256_store_si256((__m256i*)c3, mcurr);
				switch(image->nch)
				{
				case 4:
					c2_0[3]=c3[3+0];
					c2_1[3]=c3[3+8];
				case 3:
					c2_0[2]=c3[2+0];
					c2_1[2]=c3[2+8];
				case 2:
					c2_0[1]=c3[1+0];
					c2_1[1]=c3[1+8];
				case 1:
					c2_0[0]=c3[0+0];
					c2_1[0]=c3[0+8];
					break;
				}
				if(image->nch>=3)
				{
					//JPEG2000 RCT
					c2_0[1]-=(c2_0[0]+c2_0[2])>>2;
					c2_1[1]-=(c2_1[0]+c2_1[2])>>2;
					c2_0[2]+=c2_0[1];
					c2_1[2]+=c2_1[1];
					c2_0[0]+=c2_0[1];
					c2_1[0]+=c2_1[1];
				}
#ifdef ENABLE_GUIDE
				if(guide&&memcmp(c2_0, guide->data+idx0, sizeof(short)*image->nch))
				{
					short c0[4]={0};
					memcpy(c0, guide->data+idx0, sizeof(short)*image->nch);
					c0[0]-=c0[1];
					c0[2]-=c0[1];
					c0[1]+=(c0[0]+c0[2])>>2;
					c2_0[0]-=c2_0[1];
					c2_0[2]-=c2_0[1];
					c2_0[1]+=(c2_0[0]+c2_0[2])>>2;
					LOG_ERROR("Guide error IDX %d/%d", idx0, image->nch*image->iw*image->ih);
					printf("");//
				}
#endif
				__m256i step=_mm256_set1_epi64x(sizeof(int[8*2]));
				__m256i mrows=_mm256_load_si256((__m256i*)rows);
				mrows=_mm256_add_epi64(mrows, step);
				_mm256_store_si256((__m256i*)rows, mrows);
			}
		}
		if(loud)
		{
			t0=time_sec()-t0;
			printf("D  %16.6lf sec  %16.6lf MB/s\n", t0, usize/(t0*1024*1024));
		}
	}
	_mm_free(pixels);
	return 0;
}