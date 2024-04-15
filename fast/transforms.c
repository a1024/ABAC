#include"fast.h"
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<intrin.h>
#include<immintrin.h>
static const char file[]=__FILE__;

void rct_JPEG2000_32(Image *image, int fwd)
{
	if(image->nch<3)
		return;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*image->nch;k<len;k+=image->nch)
		{
			short
				r=image->data[k+0],
				g=image->data[k+1],
				b=image->data[k+2];
			
			r-=g;       //r-g				[1     -1     0  ].RGB
			b-=g;       //b-g				[0     -1     1  ].RGB
			g+=(r+b)>>2;//g+(r-g+b-g)/4 = r/4+g/2+b/4	[1/4    1/2   1/4].RGB

			image->data[k+0]=g;//Y
			image->data[k+1]=b;//Cb
			image->data[k+2]=r;//Cr
		}
	}
	else
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*image->nch;k<len;k+=image->nch)
		{
			short
				Y =image->data[k+0],
				Cb=image->data[k+1],
				Cr=image->data[k+2];
			
			Y-=(Cr+Cb)>>2;
			Cb+=Y;
			Cr+=Y;

			image->data[k+0]=Cr;
			image->data[k+1]=Y;
			image->data[k+2]=Cb;
		}
	}
}
#if 0
void rct_JPEG2000_32_simd(Image *image, int fwd)
{
	if(image->nch<3)
		return;
	ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*image->nch;
	ALIGN(32) short buf[64];
	short *pixels=image->data;
	__m256i mask=_mm256_set_epi16(
		-1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1
	);
	__m256i crosslane=_mm256_set_epi32(6, 5, 4, 3, 2, 1, 0, 2);
	__m256i inlane=_mm256_set_epi8(
		13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0,  1,  0,
		13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0, 15, 14
		//15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0,
		//15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
	);
	if(fwd)
	{
		for(;k<len-80;k+=72)
		{
			__m256i r0=_mm256_loadu_si256((__m256i*)(pixels+18*0+0));
			__m256i g0=_mm256_loadu_si256((__m256i*)(pixels+18*0+1));
			__m256i b0=_mm256_loadu_si256((__m256i*)(pixels+18*0+2));
			__m256i r1=_mm256_loadu_si256((__m256i*)(pixels+18*1+0));
			__m256i g1=_mm256_loadu_si256((__m256i*)(pixels+18*1+1));
			__m256i b1=_mm256_loadu_si256((__m256i*)(pixels+18*1+2));

			r0=_mm256_sub_epi16(r0, g0);
			r1=_mm256_sub_epi16(r1, g1);
			b0=_mm256_sub_epi16(b0, g0);
			b1=_mm256_sub_epi16(b1, g1);
			__m256i t0=_mm256_add_epi16(r0, b0);
			__m256i t1=_mm256_add_epi16(r1, b1);
			t0=_mm256_srai_epi16(t0, 2);
			t1=_mm256_srai_epi16(t1, 2);
			g0=_mm256_add_epi16(g0, t0);
			g1=_mm256_add_epi16(g1, t1);

			r0=_mm256_and_si256(r0, mask);
			g0=_mm256_and_si256(g0, mask);
			b0=_mm256_and_si256(b0, mask);
			r1=_mm256_and_si256(r1, mask);
			g1=_mm256_and_si256(g1, mask);
			b1=_mm256_and_si256(b1, mask);
			_mm256_store_si256((__m256i*)buf+0, b0);
			_mm256_store_si256((__m256i*)buf+1, r0);
			_mm256_store_si256((__m256i*)buf+2, b1);
			_mm256_store_si256((__m256i*)buf+3, r1);
			r0=_mm256_permutevar8x32_epi32(r0, crosslane);
			r1=_mm256_permutevar8x32_epi32(r1, crosslane);
			b0=_mm256_shuffle_epi8(b0, inlane);
			b1=_mm256_shuffle_epi8(b1, inlane);
			g0=_mm256_or_si256(g0, b0);
			g1=_mm256_or_si256(g1, b1);
			g0=_mm256_or_si256(g0, r0);
			g1=_mm256_or_si256(g1, r1);
			_mm256_storeu_si256((__m256i*)(pixels+18*0+0), g0);
			_mm256_storeu_si256((__m256i*)(pixels+18*1+0), g1);
			pixels[];
			//_mm256_storeu_si256((__m256i*)(pixels+18*0+0), g0);//Y
			//_mm256_store_si256((__m256i*)buf+0, b0);//Cb
			//_mm256_store_si256((__m256i*)buf+1, r0);//Cr
			//for(int k2=0;k2<6;++k2)
			//{
			//	pixels[3*k2+1]=buf[16*0+3*k2];
			//	pixels[3*k2+2]=buf[16*1+3*k2];
			//}
			//pixels[3*0+1]=buf[16*0+3*0];
			//pixels[3*0+2]=buf[16*1+3*0];
			//pixels[3*1+1]=buf[16*0+3*1];
			//pixels[3*1+2]=buf[16*1+3*1];
			//pixels[3*2+1]=buf[16*0+3*2];
			//pixels[3*2+2]=buf[16*1+3*2];
			//pixels[3*3+1]=buf[16*0+3*3];
			//pixels[3*3+2]=buf[16*1+3*3];
		}
	}
	else
	{
	}
	for(;k<len;k+=image->nch)
	{
		short
			r=image->data[k+0],
			g=image->data[k+1],
			b=image->data[k+2];
			
		r-=g;       //r-g				[1     -1     0  ].RGB
		b-=g;       //b-g				[0     -1     1  ].RGB
		g+=(r+b)>>2;//g+(r-g+b-g)/4 = r/4+g/2+b/4	[1/4    1/2   1/4].RGB

		image->data[k+0]=g;//Y
		image->data[k+1]=b;//Cb
		image->data[k+2]=r;//Cr
	}
}
#endif

//clamped gradient predictor, aka LOCO-I / Median Edge Detector (MED) predictor from JPEG-LS
void pred_clampgrad(Image *src, int fwd, const char *depths)
{
	Image dst={0};
	image_copy(&dst, src);
	if(!dst.data)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	const short *pixels=fwd?dst.data:src->data;
	int dy=src->nch*src->iw;
	int dx=src->nch;
	for(int kc=0;kc<src->nch;++kc)
	{
		int nlevels=1<<(depths?depths[kc]:src->depth);
		for(int ky=0, idx=kc;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, idx+=src->nch)
			{
				int
					NW=kx&&ky	?pixels[idx-dy	-dx	]:0,
					N =ky		?pixels[idx	-dy	]:0,
					W =kx		?pixels[idx-dx		]:0;
				int pred=N+W-NW;
				pred=MEDIAN3(N, W, pred);

				pred^=-fwd;
				pred+=fwd;

				pred+=src->data[idx];

				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;

				src->data[idx]=pred;
			}
		}
	}
	image_clear(&dst);
}
void pred_clampgrad_fast(Image *src, int fwd, const char *depths)
{
	short *pixels=(short*)malloc((src->iw+2LL)*sizeof(short[2*4]));//2 padded rows * 4 channels max
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, (src->iw+2LL)*sizeof(short[2*4]));
	int nlevels[]=
	{
		1<<(depths?depths[0]:src->depth),
		1<<(depths?depths[1]:src->depth),
		1<<(depths?depths[2]:src->depth),
		1<<(depths?depths[3]:src->depth),
	};
	int fwdmask=-fwd;
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((src->iw+2LL)*((ky-0)&1)+1LL)<<2),
			pixels+(((src->iw+2LL)*((ky-1)&1)+1LL)<<2),
		};
		for(int kx=0;kx<src->iw;++kx, idx+=src->nch)
		{
			for(int kc=0;kc<src->nch;++kc)
			{
				int
					NW	=rows[1][kc-4],
					N	=rows[1][kc+0],
					W	=rows[0][kc-4];
				int pred=N+W-NW;
				int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
				pred=CLAMP(vmin, pred, vmax);

				int curr=src->data[idx+kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=curr;

				pred+=nlevels[kc]>>1;
				pred&=nlevels[kc]-1;
				pred-=nlevels[kc]>>1;

				src->data[idx+kc]=pred;
				rows[0][kc]=fwd?curr:pred;
			}

			rows[0]+=4;
			rows[1]+=4;
		}
	}
	free(pixels);
}
void pred_simd(Image *src, int fwd, const char *depths)
{
	short *pixels=(short*)malloc((src->iw+4LL)*sizeof(short[2*4]));//2 padded rows * 4 channels max
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, (src->iw+4LL)*sizeof(short[2*4]));
	int nlevels[]=
	{
		1<<(depths?depths[0]:src->depth),
		1<<(depths?depths[1]:src->depth),
		1<<(depths?depths[2]:src->depth),
		1<<(depths?depths[3]:src->depth),
	};
	int fwdmask=-fwd;
	__m128i mfwd=_mm_set1_epi16(fwdmask);
	__m128i mhalf=_mm_set_epi16(
		0, 0, 0, 0,
		0,
		nlevels[2]>>1,
		nlevels[1]>>1,
		nlevels[0]>>1
	);
	__m128i symmask=_mm_set_epi16(
		0, 0, 0, 0,
		0,
		nlevels[2]-1,
		nlevels[1]-1,
		nlevels[0]-1
	);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((src->iw+2LL)*((ky-0)&1)+1LL)<<2),
			pixels+(((src->iw+2LL)*((ky-1)&1)+1LL)<<2),
		};
		int kx=0;
		for(;kx<src->iw;++kx, idx+=src->nch)
		{
			__m128i NW	=_mm_loadu_si128((__m128i*)(rows[1]-4));
			__m128i N	=_mm_loadu_si128((__m128i*)(rows[1]+0));
			__m128i W	=_mm_loadu_si128((__m128i*)(rows[0]-4));
			__m128i pred=_mm_add_epi16(N, W);
			pred=_mm_sub_epi16(pred, NW);
			__m128i vmin=_mm_min_epi16(N, W);
			__m128i vmax=_mm_max_epi16(N, W);
			pred=_mm_min_epi16(pred, vmax);
			pred=_mm_max_epi16(pred, vmin);

			__m128i mc=_mm_set_epi16(
				0, 0, 0, 0,
				0,
				src->data[idx+2],
				src->data[idx+1],
				src->data[idx+0]
			);
			ALIGN(16) short curr[8]={0};
			pred=_mm_xor_si128(pred, mfwd);
			pred=_mm_sub_epi16(pred, mfwd);
			pred=_mm_add_epi16(pred, mc);
			pred=_mm_add_epi16(pred, mhalf);
			pred=_mm_and_si128(pred, symmask);
			pred=_mm_sub_epi16(pred, mhalf);
			_mm_store_si128((__m128i*)curr, pred);
			src->data[idx+0]=curr[0];
			src->data[idx+1]=curr[1];
			src->data[idx+2]=curr[2];
			_mm_store_si128((__m128i*)(rows[0]+0), fwd?mc:pred);

			rows[0]+=4;
			rows[1]+=4;
		}
	}
	free(pixels);
}


	#define WP_RCT

#define NWP 4LL		//multiple of 4
#define LGBLOCKSIZE 0
#define BLOCKSIZE (1<<LGBLOCKSIZE)
#define WP_PBITS 5
#define WP_WBITS 10
void pred_wp_deferred(Image *src, int fwd)
{
	//if(src->nch!=3)
	//{
	//	LOG_ERROR("Expected 3 channels, got %d", src->nch);
	//	return;
	//}
	int nblocks=(src->iw+BLOCKSIZE-1)>>LGBLOCKSIZE;
	int *weights=(int*)_mm_malloc((nblocks+NWP)*sizeof(int[4*NWP]), sizeof(__m256i));//interleaved channels
	int *pixels=(int*)malloc((src->iw+16LL)*sizeof(int[4*4*2]));//4 padded rows * 4 channels max * {pixel, error}		int, to accelerate _mullo_epi32
	if(!weights||!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	weights[0]=(int)((1LL<<WP_WBITS)/NWP);
	memfill(weights+1, weights, (nblocks+NWP)*sizeof(int[4*NWP])-sizeof(*weights), sizeof(*weights));
	//memset(weights, 0, (nblocks+NWP)*sizeof(int[4*NWP]));
	memset(pixels, 0, (src->iw+16LL)*sizeof(int[4*4*2]));
	int nlevels[]=
	{
#ifdef WP_RCT
		2<<src->depth,
		1<<src->depth,
		2<<src->depth,
#else
		1<<src->depth,
		1<<src->depth,
		1<<src->depth,
#endif
	};
	__m128i pixelmask=_mm_set_epi32(0, nlevels[2]-1, nlevels[1]-1, nlevels[0]-1);
	__m128i pixelhalf=_mm_set_epi32(0, nlevels[2]>>1, nlevels[1]>>1, nlevels[0]>>1);
	__m128i roundoffset=_mm_set1_epi32(1<<WP_PBITS>>1);
	__m128i fwdmask=_mm_set1_epi32(-fwd);
	__m128i invmask=_mm_set1_epi32(-!fwd);
	ALIGN(16) int pred_errors[4*NWP]={0};//interleaved channels
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *wp=weights;
		int *rows[]=
		{
			pixels+(((src->iw+16LL)*((ky-0LL)&3)+4LL)<<3),
			pixels+(((src->iw+16LL)*((ky-1LL)&3)+4LL)<<3),
			pixels+(((src->iw+16LL)*((ky-2LL)&3)+4LL)<<3),
			pixels+(((src->iw+16LL)*((ky-3LL)&3)+4LL)<<3),
		};
		for(int kx=0;kx<src->iw;++kx, idx+=3)
		{
			__m128i NN	=_mm_loadu_si128((__m128i*)(rows[2]+0*8+0));
			__m128i NW	=_mm_loadu_si128((__m128i*)(rows[1]-1*8+0));
			__m128i N	=_mm_loadu_si128((__m128i*)(rows[1]+0*8+0));
			__m128i NE	=_mm_loadu_si128((__m128i*)(rows[1]+1*8+0));
			__m128i W	=_mm_loadu_si128((__m128i*)(rows[0]-1*8+0));
			__m128i eNN	=_mm_loadu_si128((__m128i*)(rows[2]+0*8+4));
			__m128i eNW	=_mm_loadu_si128((__m128i*)(rows[1]-1*8+4));
			__m128i eN	=_mm_loadu_si128((__m128i*)(rows[1]+0*8+4));
			__m128i eNE	=_mm_loadu_si128((__m128i*)(rows[1]+1*8+4));
			__m128i eW	=_mm_loadu_si128((__m128i*)(rows[0]-1*8+4));
			
			__m128i vmin=_mm_min_epi32(N, W);
			__m128i vmax=_mm_max_epi32(N, W);
			vmin=_mm_min_epi32(vmin, NE);
			vmax=_mm_max_epi32(vmax, NE);
			//wp1 = W+NE-N
			//wp2 = W+((eN+eW+eNW)*5>>4)
			//wp3 = N+((eN+eW+eNE)*5>>4)
			//wp4 = N+NW-W+((-eNN*13-eN*5-eNE*11+(N-NN)*5)>>4)
#if 0
			ALIGN(16) int spred[16]={0};
			for(int kc=0;kc<4;++kc)
			{
				spred[kc|0*4]=W.m128i_i32[kc]+NE.m128i_i32[kc]-N.m128i_i32[kc];
				spred[kc|1*4]=W.m128i_i32[kc]+((eN.m128i_i32[kc]+eW.m128i_i32[kc]+eNW.m128i_i32[kc])*5>>4);
				spred[kc|2*4]=N.m128i_i32[kc]+((eN.m128i_i32[kc]+eW.m128i_i32[kc]+eNE.m128i_i32[kc])*5>>4);
				spred[kc|3*4]=N.m128i_i32[kc]+NW.m128i_i32[kc]-W.m128i_i32[kc]
					+((-eNN.m128i_i32[kc]*13-eN.m128i_i32[kc]*5-eNE.m128i_i32[kc]*11+(N.m128i_i32[kc]-NN.m128i_i32[kc])*5)>>4);
			}
#endif
			__m128i e=_mm_add_epi32(eN, eW);
			__m128i wp1=_mm_sub_epi32(_mm_add_epi32(W, NE), N);
			__m128i wp2=_mm_add_epi32(e, eNW);
			__m128i wp3=_mm_add_epi32(e, eNE);
			__m128i wp4=_mm_sub_epi32(_mm_add_epi32(N, NW), W);
			__m128i e2=_mm_add_epi32(eNE, _mm_add_epi32(_mm_slli_epi32(eNE, 1), _mm_slli_epi32(eNE, 3)));
			e=_mm_add_epi32(eNN, _mm_add_epi32(_mm_slli_epi32(eNN, 2), _mm_slli_epi32(eNN, 3)));
			e=_mm_add_epi32(e, _mm_add_epi32(eN, _mm_slli_epi32(eN, 2)));
			e=_mm_add_epi32(e, e2);
			e2=_mm_sub_epi32(N, NN);
			e2=_mm_add_epi32(e2, _mm_slli_epi32(e2, 2));
			e2=_mm_sub_epi32(e2, e);
			e2=_mm_srai_epi32(e2, 4);
			wp2=_mm_add_epi32(W, _mm_srai_epi32(_mm_add_epi32(wp2, _mm_slli_epi32(wp2, 2)), 4));
			wp3=_mm_add_epi32(N, _mm_srai_epi32(_mm_add_epi32(wp3, _mm_slli_epi32(wp3, 2)), 4));
			wp4=_mm_add_epi32(wp4, e2);
#if 0
			//for(int kc=0;kc<4;++kc)
			//{
			//	if(wp1.m128i_i32[kc]!=spred[kc|0*4])
			//		LOG_ERROR("");
			//	if(wp2.m128i_i32[kc]!=spred[kc|1*4])
			//		LOG_ERROR("");
			//	if(wp3.m128i_i32[kc]!=spred[kc|2*4])
			//		LOG_ERROR("");
			//	if(wp4.m128i_i32[kc]!=spred[kc|3*4])
			//		LOG_ERROR("");
			//}
			wp1=_mm_load_si128((__m128i*)spred+0);//
			wp2=_mm_load_si128((__m128i*)spred+1);//
			wp3=_mm_load_si128((__m128i*)spred+2);//
			wp4=_mm_load_si128((__m128i*)spred+3);//
#endif

			__m128i w1=_mm_load_si128((__m128i*)wp+0);
			__m128i w2=_mm_load_si128((__m128i*)wp+1);
			__m128i w3=_mm_load_si128((__m128i*)wp+2);
			__m128i w4=_mm_load_si128((__m128i*)wp+3);
			w1=_mm_mullo_epi32(w1, wp1);//FIXME use AVX2
			w2=_mm_mullo_epi32(w2, wp2);
			w3=_mm_mullo_epi32(w3, wp3);
			w4=_mm_mullo_epi32(w4, wp4);
			w1=_mm_add_epi32(w1, w2);
			w1=_mm_add_epi32(w1, w3);
			w1=_mm_add_epi32(w1, w4);
			w1=_mm_srai_epi32(w1, WP_WBITS);
			w1=_mm_min_epi32(w1, vmax);
			w1=_mm_max_epi32(w1, vmin);

			ALIGN(16) int curr[]=
			{
				src->data[idx+0],
				src->data[idx+1],
				src->data[idx+2],
				//src->data[idx+1],//g		X  can't permute VYU -> YUV
				//src->data[idx+2],//b
				//src->data[idx+0],//r
				0,
			};
#ifdef WP_RCT
			if(fwd)
			{
				curr[0]-=curr[1];			//curr[0]=((curr[0]+(nlevels[0]>>1))&(nlevels[0]-1))-(nlevels[0]>>1);
				curr[2]-=curr[1];			//curr[1]=((curr[1]+(nlevels[1]>>1))&(nlevels[1]-1))-(nlevels[1]>>1);
				curr[1]+=(curr[0]+curr[2])>>2;		//curr[2]=((curr[2]+(nlevels[2]>>1))&(nlevels[2]-1))-(nlevels[2]>>1);

				//curr[2]-=curr[0];//Cr
				//curr[1]-=curr[0];//Cb
				//curr[0]+=(curr[2]+curr[1])>>2;//Y
			}
#endif
			__m128i mc=_mm_load_si128((__m128i*)curr);
			__m128i pixel=_mm_and_si128(mc, fwdmask);
			__m128i crisppred=_mm_add_epi32(w1, roundoffset);
			crisppred=_mm_srai_epi32(crisppred, WP_PBITS);
			crisppred=_mm_xor_si128(crisppred, fwdmask);
			crisppred=_mm_sub_epi32(crisppred, fwdmask);
			crisppred=_mm_add_epi32(crisppred, mc);
			crisppred=_mm_add_epi32(crisppred, pixelhalf);
			crisppred=_mm_and_si128(crisppred, pixelmask);
			crisppred=_mm_sub_epi32(crisppred, pixelhalf);
			_mm_store_si128((__m128i*)curr, crisppred);
#ifdef WP_RCT
			if(!fwd)
			{
				curr[1]-=(curr[0]+curr[2])>>2;		//curr[2]=((curr[2]+(nlevels[2]>>1))&(nlevels[2]-1))-(nlevels[2]>>1);
				curr[2]+=curr[1];			//curr[1]=((curr[1]+(nlevels[1]>>1))&(nlevels[1]-1))-(nlevels[1]>>1);
				curr[0]+=curr[1];			//curr[0]=((curr[0]+(nlevels[0]>>1))&(nlevels[0]-1))-(nlevels[0]>>1);

				//curr[0]-=(curr[2]+curr[1])>>2;//g
				//curr[1]+=curr[0];//b
				//curr[2]+=curr[0];//r
			}
#endif
			src->data[idx+0]=curr[0];
			src->data[idx+1]=curr[1];
			src->data[idx+2]=curr[2];
			//src->data[idx+0]=curr[2];//r
			//src->data[idx+1]=curr[0];//g
			//src->data[idx+2]=curr[1];//b
			pixel=_mm_or_si128(pixel, _mm_and_si128(crisppred, invmask));
			pixel=_mm_slli_epi32(pixel, WP_PBITS);
			w1=_mm_sub_epi32(pixel, w1);
			_mm_storeu_si128((__m128i*)(rows[0]+8*0+4), w1);
			_mm_storeu_si128((__m128i*)(rows[0]+8*0+0), pixel);

			__m128i pe0=_mm_load_si128((__m128i*)pred_errors+0);
			__m128i pe1=_mm_load_si128((__m128i*)pred_errors+1);
			__m128i pe2=_mm_load_si128((__m128i*)pred_errors+2);
			__m128i pe3=_mm_load_si128((__m128i*)pred_errors+3);
			wp1=_mm_sub_epi32(pixel, wp1);
			wp2=_mm_sub_epi32(pixel, wp2);
			wp3=_mm_sub_epi32(pixel, wp3);
			wp4=_mm_sub_epi32(pixel, wp4);
			wp1=_mm_abs_epi32(wp1);
			wp2=_mm_abs_epi32(wp2);
			wp3=_mm_abs_epi32(wp3);
			wp4=_mm_abs_epi32(wp4);
			pe0=_mm_add_epi32(pe0, wp1);
			pe1=_mm_add_epi32(pe1, wp2);
			pe2=_mm_add_epi32(pe2, wp3);
			pe3=_mm_add_epi32(pe3, wp4);
			_mm_store_si128((__m128i*)pred_errors+0, pe0);
			_mm_store_si128((__m128i*)pred_errors+1, pe1);
			_mm_store_si128((__m128i*)pred_errors+2, pe2);
			_mm_store_si128((__m128i*)pred_errors+3, pe3);

			if(kx&&!(kx&(BLOCKSIZE-1)))//update WP
			{
				for(int kc=0;kc<3;++kc)
				{
					long long
						w0=0x100000000/((long long)pred_errors[kc+0*4]+1),
						w1=0x100000000/((long long)pred_errors[kc+1*4]+1),
						w2=0x100000000/((long long)pred_errors[kc+2*4]+1),
						w3=0x100000000/((long long)pred_errors[kc+3*4]+1),
						sum=w0+w1+w2+w3+1;
					wp[kc+0*4]=(int)((w0<<WP_WBITS)/sum);
					wp[kc+1*4]=(int)((w1<<WP_WBITS)/sum);
					wp[kc+2*4]=(int)((w2<<WP_WBITS)/sum);
					wp[kc+3*4]=(int)((w3<<WP_WBITS)/sum);
				}
				memset(pred_errors, 0, sizeof(pred_errors));
				wp+=4*NWP;
			}
			rows[0]+=8;
			rows[1]+=8;
			rows[2]+=8;
			rows[3]+=8;
		}
	}
	_mm_free(weights);
	free(pixels);
}

void calc_csize(Image const *src, const char *depths, double *ret_csizes, double *ret_invCR)//ret: TRGBA
{
	char bitdepths[]=
	{
		depths?depths[0]:src->depth,
		depths?depths[1]:src->depth,
		depths?depths[2]:src->depth,
		depths?depths[3]:src->depth,
	};
	int nlevels[]=
	{
		1<<bitdepths[0],
		1<<bitdepths[1],
		1<<bitdepths[2],
		1<<bitdepths[3],
	};
	int maxdepth=bitdepths[0];
	UPDATE_MAX(maxdepth, bitdepths[1]);
	UPDATE_MAX(maxdepth, bitdepths[2]);
	UPDATE_MAX(maxdepth, bitdepths[3]);
	unsigned *hist=(unsigned*)malloc(sizeof(int)*src->nch<<maxdepth);
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(hist, 0, sizeof(int)*src->nch<<maxdepth);
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih, nvals=res*src->nch;
	for(ptrdiff_t k=0;k<nvals;k+=src->nch)
	{
		for(int kc=0;kc<src->nch;++kc)
		{
			int val=src->data[k+kc];
			val+=nlevels[kc]>>1;
			val&=nlevels[kc]-1;
			++hist[kc<<maxdepth|val];
		}
	}
	if(ret_invCR)
		*ret_invCR=0;
	if(ret_csizes)
		*ret_csizes=0;
	for(int kc=0;kc<src->nch;++kc)
	{
		double csize=0;
		for(int ks=0;ks<nlevels[kc];++ks)
		{
			unsigned freq=hist[kc<<maxdepth|ks];
			if(freq)
				csize-=freq*log2((double)freq/res);
		}
		if(ret_invCR)
		{
			ret_invCR[kc+1]=csize/((double)res*src->depth);
			*ret_invCR+=ret_invCR[kc+1];
		}
		csize/=8;
		if(ret_csizes)
		{
			ret_csizes[kc+1]=csize;
			*ret_csizes+=ret_csizes[kc+1];
		}
	}
	if(ret_invCR)
		*ret_invCR/=src->nch;
	free(hist);
}
void calc_csize_vlc(Image const *src, const char *depths, double *ret_csizes, double *ret_csizes_vlc)
{
	char bitdepths[]=
	{
		depths?depths[0]:src->depth,
		depths?depths[1]:src->depth,
		depths?depths[2]:src->depth,
		depths?depths[3]:src->depth,
	};
	int nlevels[]=
	{
		1<<bitdepths[0],
		1<<bitdepths[1],
		1<<bitdepths[2],
		1<<bitdepths[3],
	};
	int maxdepth=bitdepths[0];
	UPDATE_MAX(maxdepth, bitdepths[1]);
	UPDATE_MAX(maxdepth, bitdepths[2]);
	UPDATE_MAX(maxdepth, bitdepths[3]);
	unsigned *hist=(unsigned*)malloc(sizeof(int)*src->nch<<maxdepth);
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	*hist=1;
	memfill(hist+1, hist, (sizeof(int)*src->nch<<maxdepth)-sizeof(int), sizeof(int));
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih, nvals=res*src->nch;
	double csizes[4]={0}, csizes_vlc[4]={0};
	for(ptrdiff_t k=0, idx=1;k<nvals;k+=src->nch, ++idx)
	{
		for(int kc=0;kc<src->nch;++kc)
		{
			int val=src->data[k+kc];
			val+=nlevels[kc]>>1;
			val&=nlevels[kc]-1;
			int idx2=kc<<maxdepth|val;
			unsigned freq=hist[idx2];
			csizes[kc]-=log2((double)freq/idx);
			csizes_vlc[kc]+=(double)ceil_log2_32(((unsigned)idx))-ceil_log2_32((unsigned)freq)+1;//prefix bit
			++freq;
			hist[idx2]=freq;
		}
	}
	for(int kc=0;kc<4;++kc)
	{
		csizes[kc]/=8;
		csizes_vlc[kc]/=8;
	}
	if(ret_csizes)
	{
		ret_csizes[0]=csizes[0]+csizes[1]+csizes[2]+csizes[3];
		memcpy(ret_csizes+1, csizes, sizeof(csizes));
	}
	if(ret_csizes_vlc)
	{
		ret_csizes_vlc[0]=csizes_vlc[0]+csizes_vlc[1]+csizes_vlc[2]+csizes_vlc[3];
		memcpy(ret_csizes_vlc+1, csizes_vlc, sizeof(csizes_vlc));
	}
	free(hist);
}
void calc_csize_bin(Image const *src, const char *depths, double *ret_csizes)
{
	char bitdepths[]=
	{
		depths?depths[0]:src->depth,
		depths?depths[1]:src->depth,
		depths?depths[2]:src->depth,
		depths?depths[3]:src->depth,
	};
	int maxdepth=bitdepths[0];
	UPDATE_MAX(maxdepth, bitdepths[1]);
	UPDATE_MAX(maxdepth, bitdepths[2]);
	UPDATE_MAX(maxdepth, bitdepths[3]);
	unsigned short *stats=(unsigned short*)malloc(sizeof(short)*src->nch<<maxdepth);
	if(!stats)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	*stats=0x8000;
	memfill(stats+1, stats, (sizeof(short)*src->nch<<maxdepth)-sizeof(short), sizeof(short));
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih, nvals=res*src->nch;
	double csizes[4]={0};
	for(ptrdiff_t k=0, idx=1;k<nvals;k+=src->nch, ++idx)
	{
		for(int kc=0;kc<src->nch;++kc)
		{
			unsigned short *curr_stats=stats+((size_t)kc<<maxdepth);
			int idx2=1;
			for(int kb=bitdepths[kc]-1;kb>=0;--kb)
			{
				int p0=curr_stats[idx2];
				int bit=src->data[k+kc]>>kb&1;
				int p=bit?0x10000-p0:p0;
				csizes[kc]-=log2((double)p/0x10000);
				int update=((!bit<<16)-p0);
				p0+=(update>>7)+(update<0);
				//if(!(unsigned short)p0)
				//	LOG_ERROR("");
				curr_stats[idx2]=p0;
				idx2=idx2<<1|bit;
			}
		}
	}
	*ret_csizes=0;
	for(int kc=0;kc<src->nch;++kc)
		*ret_csizes+=csizes[kc]/=8;
	memcpy(ret_csizes+1, csizes, sizeof(csizes));
	free(stats);
}
size_t calc_csize_ABAC(Image const *src, const char *depths)
{
	char bitdepths[]=
	{
		depths?depths[0]:src->depth,
		depths?depths[1]:src->depth,
		depths?depths[2]:src->depth,
		depths?depths[3]:src->depth,
	};
	int maxdepth=bitdepths[0];
	UPDATE_MAX(maxdepth, bitdepths[1]);
	UPDATE_MAX(maxdepth, bitdepths[2]);
	UPDATE_MAX(maxdepth, bitdepths[3]);
	unsigned short *stats=(unsigned short*)malloc(sizeof(short)*src->nch<<maxdepth);
	if(!stats)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	*stats=0x8000;
	memfill(stats+1, stats, (sizeof(short)*src->nch<<maxdepth)-sizeof(short), sizeof(short));
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih, nvals=res*src->nch;
	size_t csize=0;
	//unsigned long long low=0, range=0xFFFFFFFF;
	unsigned range=0xFFFFFFFF;
	for(ptrdiff_t k=0, idx=1;k<nvals;k+=src->nch, ++idx)
	{
		for(int kc=0;kc<src->nch;++kc)
		{
			int val=src->data[k+kc];
			unsigned short *curr_stats=stats+((size_t)kc<<maxdepth);
			unsigned long long idx2=1;
			int kb=bitdepths[kc]-1;

			unsigned short p0;
			int update;
			unsigned long long bit, r2;
#if 1
			for(int kb2=0;kb2<8;++kb2)//compiler unrolls this
			{
				p0=curr_stats[idx2];
				bit=(unsigned long long)(val>>kb&1);
				r2=(unsigned long long)range*p0>>16;
				//low+=r2&-bit;
				//range=bit?range-r2:r2-1;
				if(bit)
				{
					//low+=r2;
					range-=r2;
				}
				else
					range=r2-1;
				bit^=1;
				if(range<0x10000)
				{
					range<<=16;
					//low<<=16;
					range|=0xFFFF;
					csize+=2;
				}
				update=(bit<<16)-p0;
				curr_stats[idx2]=(unsigned short)(p0+(update>>8)+(update>>31));
				idx2=idx2<<1|bit;
				--kb;
			}
#endif
			do
			{
				p0=curr_stats[idx2];
				bit=(unsigned long long)(val>>kb&1);
				r2=(unsigned long long)range*p0>>16;
				//low+=r2&-bit;
				//range=bit?range-r2:r2-1;
				if(bit)
				{
					//low+=r2;
					range-=r2;
				}
				else
					range=r2-1;
				bit^=1;
				if(range<0x10000)
				{
					range<<=16;
					//low<<=16;
					range|=0xFFFF;
					csize+=2;
				}
				update=(bit<<16)-p0;
				curr_stats[idx2]=(unsigned short)(p0+(update>>8)+(update>>31));//rounding towards zero
				idx2=idx2<<1|bit;
				--kb;
			}while(kb>=0);
		}
	}
	csize+=8;
	free(stats);
	return csize;
}