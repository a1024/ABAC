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


//	#define PROFILE_CSIZE
//	#define ENABLE_GUIDE
//	#define ALLOW_IDIV
//	#define SAVE_K
//	#define MEASURE_GR_PENALTY	//it's ~1.1%
//	#define INSTRUCTION_SPEED


#include"ac.h"
#ifdef ENABLE_GUIDE
static const Image *guide=0;
#endif
static void calc_ctx(const __m128i *N, const __m128i *W, const __m128i *NW, const __m128i *eN, const __m128i *eW, const __m128i *eNE, const __m128i *eNW, const __m128i *mminmag, __m128i *pred, __m128i *mag)
{
	//pred = median3(N, W, N+W-NW)
	//mag = eW
	__m128i mmin=_mm_min_epi32(*N, *W);
	__m128i mmax=_mm_max_epi32(*N, *W);
	*pred=_mm_add_epi32(*N, *W);
	*mag=*eW;
	//*mag=_mm_add_epi32(*eN, *eW);
	//*mag=_mm_add_epi32(*mag, *eNE);
	//*mag=_mm_add_epi32(*mag, *eNW);
	*pred=_mm_sub_epi32(*pred, *NW);
	//*mag=_mm_srai_epi32(*mag, 2);//mag = (eN+eW)>>1		//mag = (4*(eN+eW)+eNE-eNW)>>2
	*pred=_mm_max_epi32(*pred, mmin);
	*mag=_mm_max_epi32(*mag, *mminmag);
	*pred=_mm_min_epi32(*pred, mmax);

	(void)eN;
	(void)eNE;
	(void)eNW;
}
static void pack_sign(__m128i *data)
{
	//x = curr-pred
	//x = x<<1^-(x<0)
	__m128i negmask=_mm_cmpgt_epi32(_mm_setzero_si128(), *data);
	*data=_mm_slli_epi32(*data, 1);
	*data=_mm_xor_si128(*data, negmask);
}
static void unpack_sign(__m128i *curr)
{
	//curr = (curr>>1^-(curr&1)) + pred
	__m128i negmask=_mm_and_si128(*curr, _mm_set1_epi32(1));
	negmask=_mm_cmpeq_epi32(negmask, _mm_set1_epi32(1));
	*curr=_mm_srai_epi32(*curr, 1);
	*curr=_mm_xor_si128(*curr, negmask);
}
static void calc_update(const __m128i *eWW, const __m128i *eW, const __m128i *eNEEE, __m128i *ecurr)//writes to hi64 bits
{
	//x = (2*eW+eNEEE+ecurr)>>2
	*ecurr=_mm_add_epi32(*ecurr, _mm_slli_epi32(*eW, 1));
	*ecurr=_mm_add_epi32(*ecurr, *eNEEE);
	*ecurr=_mm_srli_epi32(*ecurr, 2);

	//x = (eWW+eW+eNEEE+ecurr)>>2
	//*ecurr=_mm_add_epi32(*ecurr, *eW);
	//*ecurr=_mm_add_epi32(*ecurr, *eWW);
	//*ecurr=_mm_add_epi32(*ecurr, *eNEEE);
	//*ecurr=_mm_srli_epi32(*ecurr, 2);
	
	//x = (eW+eNEEE+2*ecurr)>>2
	//*ecurr=_mm_slli_epi32(*ecurr, 1);
	//*ecurr=_mm_add_epi32(*ecurr, *eW);
	//*ecurr=_mm_add_epi32(*ecurr, *eNEEE);
	//*ecurr=_mm_srli_epi32(*ecurr, 2);

	//x = (eN+eW+eNEEE+ecurr)>>2
	//*ecurr=_mm_add_epi32(*ecurr, *eN);
	//*ecurr=_mm_add_epi32(*ecurr, *eW);
	//*ecurr=_mm_add_epi32(*ecurr, *eNEEE);
	//*ecurr=_mm_srli_epi32(*ecurr, 2);
	
	//x = (2*eW-eN+eNEEE+ecurr)/3		X
	//*ecurr=_mm_sub_epi32(*ecurr, *eN);
	//*ecurr=_mm_add_epi32(*ecurr, _mm_slli_epi32(*eW, 1));
	//*ecurr=_mm_add_epi32(*ecurr, *eNEEE);
	//*ecurr=_mm_add_epi32(*ecurr, _mm_slli_epi32(*ecurr, 2));//x/3 ~= x*0x55>>8 = x*5*17>>8
	//*ecurr=_mm_add_epi32(*ecurr, _mm_slli_epi32(*ecurr, 4));
	//*ecurr=_mm_srai_epi32(*ecurr, 8);//need arithmetic shift because of subtraction

	//x = (eN+2*eW+eNEEE+ecurr)/5		X
	//*ecurr=_mm_add_epi32(*ecurr, *eN);
	//*ecurr=_mm_add_epi32(*ecurr, _mm_slli_epi32(*eW, 1));
	//*ecurr=_mm_add_epi32(*ecurr, *eNEEE);
	//*ecurr=_mm_add_epi32(*ecurr, _mm_slli_epi32(*ecurr, 1));//x/5 ~= x*0x33>>8 = x*3*17>>8
	//*ecurr=_mm_add_epi32(*ecurr, _mm_slli_epi32(*ecurr, 4));
	//*ecurr=_mm_srli_epi32(*ecurr, 8);

	//(void)eN;
}
int f20_codec(Image const *src, ArrayHandle *data, const unsigned char *cbuf, size_t clen, Image *dst, int loud)
{
	double t0=time_sec();
	int fwd=src!=0;
	Image const *image=fwd?src:dst;

	size_t ebufsize=sizeof(int[4*8])*(image->iw+16LL);//4 padded rows * 4 channels max * {pixels, errors}
	int *pixels=(int*)_mm_malloc(ebufsize, sizeof(__m128i));
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return 1;
	}
	int minmag=1<<image->depth>>7;
	pixels[0]=0;
	pixels[1]=0;
	pixels[2]=0;
	pixels[3]=0;
	pixels[4]=minmag;
	pixels[5]=minmag;
	pixels[6]=minmag;
	pixels[7]=minmag;
	memfill(pixels+8, pixels, ebufsize-sizeof(int[8]), sizeof(int[8]));
	__m128i mminmag=_mm_set1_epi32(minmag);
	if(image->depth==8)
	{
		mminmag=_mm_srli_epi32(mminmag, 2);
		minmag>>=2;
	}
	ptrdiff_t usize=(ptrdiff_t)image->iw*image->ih*image->nch*image->depth>>3;
	if(fwd)
	{
#ifdef INSTRUCTION_SPEED
		if(loud)
		{
			Image im2={0};
			image_copy(&im2, image);
			if(!im2.data)
			{
				LOG_ERROR("Alloc error");
				return 1;
			}
			ptrdiff_t nvals=(ptrdiff_t)image->iw*image->ih*image->nch;
			size_t result=0;
			volatile double t1=time_sec();
			for(ptrdiff_t k=0;k<nvals-1;++k)
				result+=im2.data[k]=im2.data[k]*im2.data[k+1]>>im2.depth;
			t1=time_sec()-t1;
			printf("MUL>>d       %16lf sec %16lf MB/s  result %16zd\n", t1, usize/(t1*1024*1024), result);
			
			t1=time_sec();
			memcpy(im2.data, image->data, usize);
			t1=time_sec()-t1;
			printf("memcpy       %16lf sec %16lf MB/s\n", t1, usize/(t1*1024*1024));
			
			t1=time_sec();
			result=0;
			for(ptrdiff_t k=0;k<nvals-1;++k)
				result+=(int)(im2.data[k]<<im2.depth)/(abs(im2.data[k+1])+1);
			t1=time_sec()-t1;
			printf("DIV<<d       %16lf sec %16lf MB/s  result %16zd\n", t1, usize/(t1*1024*1024), result);
			
			//t1=time_sec();
			//memcpy(im2.data, image->data, usize);
			//t1=time_sec()-t1;
			//printf("memcpy  %16lf sec %16lf MB/s", t1, usize/(t1*1024*1024));
			
			t1=time_sec();
			result=0;
			for(ptrdiff_t k=0;k<nvals-1;++k)
				result+=FLOOR_LOG2_P1(im2.data[k]);
			t1=time_sec()-t1;
			printf("LZCNT        %16lf sec %16lf MB/s  result %16zd\n", t1, usize/(t1*1024*1024), result);
			
			t1=time_sec();
			result=0;
			for(ptrdiff_t k=0;k<nvals-1;++k)
				result^=im2.data[k]-=im2.data[k+1];
			t1=time_sec()-t1;
			printf("SUB          %16lf sec %16lf MB/s  result %16zd\n", t1, usize/(t1*1024*1024), result);
			
			{
				t1=time_sec();
				__m128i mres=_mm_setzero_si128();
				for(ptrdiff_t k=0;k<nvals-(ptrdiff_t)(2*sizeof(__m128i)/sizeof(short)-1);k+=sizeof(__m128i)/sizeof(short))
				{
					__m128i v0=_mm_loadu_si128((__m128i*)(im2.data+k));
					__m128i v1=_mm_loadu_si128((__m128i*)(im2.data+k+sizeof(__m128i)/sizeof(short)));
					v0=_mm_mullo_epi16(v0, v1);
					mres=_mm_add_epi16(mres, v0);
				}
				ALIGN(16) short r2[sizeof(__m128i)/sizeof(short)]={0};
				_mm_store_si128((__m128i*)r2, mres);
				t1=time_sec()-t1;
				printf("MULLO16_128  %16lf sec %16lf MB/s  result %16zd %16zd\n",
					t1, usize/(t1*1024*1024),
					((long long*)&r2)[0],
					((long long*)&r2)[1]
				);
			}
			
			{
				t1=time_sec();
				__m256i mres=_mm256_setzero_si256();
				for(ptrdiff_t k=0;k<nvals-(ptrdiff_t)(2*sizeof(__m256i)/sizeof(short)-1);k+=sizeof(__m256i)/sizeof(short))
				{
					__m256i v0=_mm256_loadu_si256((__m256i*)(im2.data+k));
					__m256i v1=_mm256_loadu_si256((__m256i*)(im2.data+k+sizeof(__m256i)/sizeof(short)));
					v0=_mm256_mullo_epi16(v0, v1);
					mres=_mm256_add_epi16(mres, v0);
				}
				ALIGN(16) short r2[sizeof(__m256i)/sizeof(short)]={0};
				_mm256_store_si256((__m256i*)r2, mres);
				t1=time_sec()-t1;
				printf("MULLO16_256  %16lf sec %16lf MB/s  result %16zd %16zd %16zd %16zd\n",
					t1, usize/(t1*1024*1024),
					((long long*)&r2)[0],
					((long long*)&r2)[1],
					((long long*)&r2)[2],
					((long long*)&r2)[3]
				);
			}
			
			int cpu_features=get_cpu_features();
			if(cpu_features>>1)
			{
				t1=time_sec();
				__m512i mres=_mm512_setzero_si512();
				for(ptrdiff_t k=0;k<nvals-(ptrdiff_t)(2*sizeof(__m512i)/sizeof(short)-1);k+=sizeof(__m512i)/sizeof(short))
				{
					__m512i v0=_mm512_loadu_si512((__m512i*)(im2.data+k));
					__m512i v1=_mm512_loadu_si512((__m512i*)(im2.data+k+sizeof(__m512i)/sizeof(short)));
					v0=_mm512_mullo_epi16(v0, v1);
					mres=_mm512_add_epi16(mres, v0);
				}
				ALIGN(16) short r2[sizeof(__m512i)/sizeof(short)]={0};
				_mm512_store_si512((__m512i*)r2, mres);
				t1=time_sec()-t1;
				printf("MULLO16_512  %16lf sec %16lf MB/s  result %16zd %16zd %16zd %16zd %16zd %16zd %16zd %16zd\n",
					t1, usize/(t1*1024*1024),
					((long long*)&r2)[0],
					((long long*)&r2)[1],
					((long long*)&r2)[2],
					((long long*)&r2)[3],
					((long long*)&r2)[4],
					((long long*)&r2)[5],
					((long long*)&r2)[6],
					((long long*)&r2)[7]
				);
			}
			image_clear(&im2);
		}
#endif
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
#ifdef MEASURE_GR_PENALTY
		size_t statssize=sizeof(short)*image->nch<<image->depth;
		unsigned short *stats=(unsigned short*)malloc(statssize);
		if(!stats)
		{
			LOG_ERROR("Alloc error");
			return 1;
		}
		*stats=0x8000;
		memfill(stats+1, stats, statssize-sizeof(short), sizeof(short));
		unsigned long long low=0, range=0xFFFFFFFFFFFF;
		size_t abacsize=0;
#endif
		GolombRiceCoder ec;
		DList list;
		dlist_init(&list, 1, 1024, 0);
		gr_enc_init(&ec, &list);
		ALIGN(16) short c2[8]={0};
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			int *rows[]=
			{
				pixels+(((image->iw+16LL)*((ky-0LL)&3)+8)<<3),
				pixels+(((image->iw+16LL)*((ky-1LL)&3)+8)<<3),
				pixels+(((image->iw+16LL)*((ky-2LL)&3)+8)<<3),
				pixels+(((image->iw+16LL)*((ky-3LL)&3)+8)<<3),
			};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				__m128i
					NW	=_mm_load_si128((__m128i*)rows[1]-1*2+0),
					N	=_mm_load_si128((__m128i*)rows[1]+0*2+0),
					W	=_mm_load_si128((__m128i*)rows[0]-1*2+0),
					eNW	=_mm_load_si128((__m128i*)rows[1]-1*2+1),
					eN	=_mm_load_si128((__m128i*)rows[1]+0*2+1),
					eNE	=_mm_load_si128((__m128i*)rows[1]+1*2+1),
					eNEEE	=_mm_load_si128((__m128i*)rows[1]+3*2+1),
					eWW	=_mm_load_si128((__m128i*)rows[0]-2*2+1),
					eW	=_mm_load_si128((__m128i*)rows[0]-1*2+1);
				int	*curr	=rows[0]+0*8;
				__m128i mcurr, mpred, mmag;
				calc_ctx(&N, &W, &NW, &eN, &eW, &eNE, &eNW, &mminmag, &mpred, &mmag);
				//if(ky==1&&kx==5)//
				//if(idx==1170)//
				//if(idx==6936)//
				//if(idx==4722)//
				//if(idx==4731)//
				//	printf("");
				ALIGN(16) int mag[4];
				_mm_store_si128((__m128i*)mag, mmag);
				switch(image->nch)
				{
				case 4:c2[3]=image->data[idx+3];
				case 3:c2[2]=image->data[idx+2];
				case 2:c2[1]=image->data[idx+1];
				case 1:c2[0]=image->data[idx+0];
					break;
				}
				if(image->nch>=3)
				{
					//JPEG2000 RCT
					c2[0]-=c2[1];
					c2[2]-=c2[1];
					c2[1]+=(c2[0]+c2[2])>>2;
				}
				mcurr=_mm_load_si128((__m128i*)c2);
				mcurr=_mm_cvtepi16_epi32(mcurr);
				__m128i delta=_mm_sub_epi32(mcurr, mpred);
				pack_sign(&delta);
				_mm_store_si128((__m128i*)curr+0, mcurr);
				_mm_store_si128((__m128i*)curr+1, delta);
				switch(image->nch)
				{
#ifdef ALLOW_IDIV
				case 4:gr_enc(&ec, curr[4+3], mag[3]+1);
				case 3:gr_enc(&ec, curr[4+2], mag[2]+1);
				case 2:gr_enc(&ec, curr[4+1], mag[1]+1);
				case 1:gr_enc(&ec, curr[4+0], mag[0]+1);
#else
				case 4:gr_enc_POT(&ec, curr[4+3], FLOOR_LOG2(mag[3]+1));
				case 3:gr_enc_POT(&ec, curr[4+2], FLOOR_LOG2(mag[2]+1));
				case 2:gr_enc_POT(&ec, curr[4+1], FLOOR_LOG2(mag[1]+1));
				case 1:gr_enc_POT(&ec, curr[4+0], FLOOR_LOG2(mag[0]+1));
#endif
					break;
				}
#ifdef SAVE_K
				if(loud)
				{
					int temp;
					temp=1<<FLOOR_LOG2(mag[0]+1), kimage.data[idx+0]=CLAMP(-half, temp, half-1)-half;
					temp=1<<FLOOR_LOG2(mag[1]+1), kimage.data[idx+1]=CLAMP(-half, temp, half-1)-half;
					temp=1<<FLOOR_LOG2(mag[2]+1), kimage.data[idx+2]=CLAMP(-half, temp, half-1)-half;
				}
#endif
#ifdef MEASURE_GR_PENALTY
				for(int kc=0;kc<image->nch;++kc)
				{
					int tidx=1;
					unsigned short *curr_stats=stats+((size_t)kc<<image->depth);
					for(int kb=image->depth-1;kb>=0;--kb)
					{
						int p0=curr_stats[tidx];
						int bit=curr[kc+4]>>kb&1;

						unsigned long long r2=range*p0>>16;
						low+=r2&-bit;
						range=bit?range-r2:r2-1;
						while(range<0x10000)
						{
							abacsize+=2;

							range<<=16;
							low<<=16;
							range|=0xFFFF;
							low&=0xFFFFFFFFFFFF;
							unsigned long long rmax=low^0xFFFFFFFFFFFF;
							if(range>rmax)//clamp hi to register size after renorm
								range=rmax;
						}

						p0+=((!bit<<16)-p0)>>7;
						curr_stats[tidx]=CLAMP(1, p0, 0xFFFF);
						tidx=tidx<<1|bit;
					}
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
				calc_update(&eWW, &eW, &eNEEE, &delta);
				_mm_store_si128((__m128i*)curr+1, delta);
				rows[0]+=8;
				rows[1]+=8;
				rows[2]+=8;
				rows[3]+=8;
			}
		}
		gr_enc_flush(&ec);
		dlist_appendtoarray(&list, data);
#ifdef MEASURE_GR_PENALTY
		unsigned long long code=low+range;
		int n=FLOOR_LOG2_P1(low^code);
		code&=~((1LL<<n)-1);		//minimize final code
		int flushbits=LSB_IDX_64(code);
		abacsize+=(flushbits+7LL)>>3;
		free(stats);
#endif
		if(loud)
		{
#ifdef SAVE_K
			image_snapshot8(&kimage);
			image_clear(&kimage);
#endif
			ptrdiff_t csize=list.nobj;
			t0=time_sec()-t0;
			printf("%14td/%14td = %10.6lf%%  CR %lf\n", csize, usize, 100.*csize/usize, (double)usize/csize);
#ifdef MEASURE_GR_PENALTY
			printf("%14td/%14td = %10.6lf%%  CR %lf  ABAC\n", abacsize, usize, 100.*abacsize/usize, (double)usize/abacsize);
#endif
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
		dlist_clear(&list);
	}
	else
	{
		GolombRiceCoder ec;
		gr_dec_init(&ec, cbuf, cbuf+clen);
		__m128i packpixels=_mm_set_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1, 13, 12,  9,  8,  5,  4,  1,  0
			//15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
		);
		for(int ky=0, idx=0;ky<image->ih;++ky)
		{
			int *rows[]=
			{
				pixels+(((image->iw+16LL)*((ky-0LL)&3)+8)<<3),
				pixels+(((image->iw+16LL)*((ky-1LL)&3)+8)<<3),
				pixels+(((image->iw+16LL)*((ky-2LL)&3)+8)<<3),
				pixels+(((image->iw+16LL)*((ky-3LL)&3)+8)<<3),
			};
			for(int kx=0;kx<image->iw;++kx, idx+=image->nch)
			{
				__m128i
					NW	=_mm_load_si128((__m128i*)rows[1]-1*2+0),
					N	=_mm_load_si128((__m128i*)rows[1]+0*2+0),
					W	=_mm_load_si128((__m128i*)rows[0]-1*2+0),
					eNW	=_mm_load_si128((__m128i*)rows[1]-1*2+1),
					eN	=_mm_load_si128((__m128i*)rows[1]+0*2+1),
					eNE	=_mm_load_si128((__m128i*)rows[1]+1*2+1),
					eNEEE	=_mm_load_si128((__m128i*)rows[1]+3*2+1),
					eWW	=_mm_load_si128((__m128i*)rows[0]-2*2+1),
					eW	=_mm_load_si128((__m128i*)rows[0]-1*2+1);
				int	*curr	=rows[0]+0*8;
				__m128i mpred, mmag;
				calc_ctx(&N, &W, &NW, &eN, &eW, &eNE, &eNW, &mminmag, &mpred, &mmag);
				//if(ky==1&&kx==5)//
				//if(idx==1161)//
				//if(idx==6936)//
				//if(idx==4722)//
				//	printf("");
				ALIGN(16) int mag[4];
				_mm_store_si128((__m128i*)mag, mmag);
				switch(image->nch)
				{
#ifdef ALLOW_IDIV
				case 4:curr[3+4]=gr_dec(&ec, mag[3]+1);
				case 3:curr[2+4]=gr_dec(&ec, mag[2]+1);
				case 2:curr[1+4]=gr_dec(&ec, mag[1]+1);
				case 1:curr[0+4]=gr_dec(&ec, mag[0]+1);
#else
				case 4:curr[3+4]=gr_dec_POT(&ec, FLOOR_LOG2(mag[3]+1));
				case 3:curr[2+4]=gr_dec_POT(&ec, FLOOR_LOG2(mag[2]+1));
				case 2:curr[1+4]=gr_dec_POT(&ec, FLOOR_LOG2(mag[1]+1));
				case 1:curr[0+4]=gr_dec_POT(&ec, FLOOR_LOG2(mag[0]+1));
#endif
					break;
				}
				__m128i delta=_mm_load_si128((__m128i*)curr+1);
				__m128i mcurr=delta;
				unpack_sign(&mcurr);
				calc_update(&eWW, &eW, &eNEEE, &delta);
				mcurr=_mm_add_epi32(mcurr, mpred);
				_mm_store_si128((__m128i*)curr+0, mcurr);
				_mm_store_si128((__m128i*)curr+1, delta);

				mcurr=_mm_shuffle_epi8(mcurr, packpixels);
				short *c2=image->data+idx;
				ALIGN(16) short c3[8];
				_mm_store_si128((__m128i*)c3, mcurr);
				switch(image->nch)
				{
				case 4:c2[3]=c3[3];
				case 3:c2[2]=c3[2];
				case 2:c2[1]=c3[1];
				case 1:c2[0]=c3[0];
					break;
				}
				if(image->nch>=3)
				{
					//JPEG2000 RCT
					c2[1]-=(c2[0]+c2[2])>>2;
					c2[2]+=c2[1];
					c2[0]+=c2[1];
				}
#ifdef ENABLE_GUIDE
				if(guide&&memcmp(c2, guide->data+idx, sizeof(short)*image->nch))
				{
					short c0[4]={0};
					memcpy(c0, guide->data+idx, sizeof(short)*image->nch);
					c0[0]-=c0[1];
					c0[2]-=c0[1];
					c0[1]+=(c0[0]+c0[2])>>2;
					c2[0]-=c2[1];
					c2[2]-=c2[1];
					c2[1]+=(c2[0]+c2[2])>>2;
					LOG_ERROR("Guide error IDX %d/%d", idx, image->nch*image->iw*image->ih);
					printf("");//
				}
#endif
				rows[0]+=8;
				rows[1]+=8;
				rows[2]+=8;
				rows[3]+=8;
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