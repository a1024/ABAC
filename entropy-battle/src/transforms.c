#include"battle.h"
#include<stdio.h>//for debugging
#include<string.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#include<tmmintrin.h>
static const char file[]=__FILE__;

static void print_buf(short *buf, int bw, int bh, int kc, int bytestride)
{
	static int call=0;
	//for(int kc=0;kc<nch;++kc)
	{
		printf("call %d ch %d\n", call, kc);
		for(int ky=0;ky<bh;++ky)
		{
			for(int kx=0;kx<bw;++kx)
				printf("%d\t", buf[bytestride*(bw*ky+kx)+kc]);
			printf("\n");
		}
		printf("\n");
	}
	++call;
}

//Haar DWT: 8 bit -> 9 bit
void haar_1d_fwd(short *buffer, int count, int stride, short *b2)
{
	int floorhalf=count>>1, ceilhalf=count-floorhalf;//ceil(count/2)			//ordinarily (but not here), nOdd <= nEven			nOdd = nEven - (count&1)
	short *even=b2, *odd=b2+ceilhalf;
	
	for(int k=0, ks=0;k<floorhalf;++k, ks+=stride<<1)//lazy wavelet: split into even (low frequency) & odd (high frequency), contrary to DWT notation			TODO _mm_shuffle_epi8?
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(count&1)
		even[floorhalf]=buffer[stride*(count-1)];

	for(int k=0;k<floorhalf;++k)
	{
		short
			av=(even[k]+odd[k]+(even[k]>odd[k]))>>1,//average is rounded towards the even (first) value		from JPEG XL
			diff=even[k]-odd[k];

		//int neg=diff<0;
		//diff^=-neg;
		//diff+=neg;
		//diff<<=1;
		//diff|=neg;

		even[k]=av;
		odd[k]=diff;
	}
	//dwt2_1d_9_7(even, odd, halfsize);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void haar_2d_fwd(const unsigned char *buf, int bw, int bh, int nch, int bytestride, int nstages, short **ret)//lossless DWT, don't forget to free ret if was zero
{
	int tsize=bw>bh?bw:bh;
	short *temp=(short*)malloc(tsize*sizeof(short));

	size_t dstlen=(size_t)bw*bh*bytestride;
	if(!*ret)
	{
		*ret=(short*)malloc(dstlen*sizeof(short));
		if(!*ret)
		{
			free(temp);
			return;
		}
	}
	for(int k=0;k<dstlen;++k)
		ret[0][k]=buf[k];
	//if(*ret!=buf)
	//	memcpy(*ret, buf, dstlen);

	int rowlen=bytestride*bw;
	for(int kc=0;kc<nch;++kc)
	{
		//print_buf(*ret, bw, bh, kc, bytestride);//

		for(int w2=bw, h2=bh, it=0;w2>=3&&h2>=3&&(!nstages||it<nstages);++it)
		{
			for(int ky=0;ky<h2;++ky)//horizontal DWT
				haar_1d_fwd(*ret+rowlen*ky+kc, w2, bytestride, temp);
			
			//print_buf(*ret, bw, bh, kc, bytestride);//

			for(int kx=0;kx<w2;++kx)//vertical DWT
				haar_1d_fwd(*ret+bytestride*kx+kc, h2, rowlen, temp);
			
			//print_buf(*ret, bw, bh, kc, bytestride);//

			w2-=w2>>1;//w=ceil(w/2)
			h2-=h2>>1;//h=ceil(h/2)
		}
	}
	free(temp);
}

void haar_1d_inv(short *buffer, int count, int stride, short *b2)
{
	int floorhalf=count>>1, ceilhalf=count-floorhalf;//ceil(count/2)
	short *even=b2, *odd=b2+ceilhalf;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];

	for(int k=0;k<floorhalf;++k)
	{
		short av=even[k], diff=odd[k];
		
		//int neg=diff&1;
		//diff>>=1;
		//diff^=-neg;
		//diff+=neg;
		//diff|=neg<<7&-!diff;

		even[k]=((av<<1)+diff+(diff>0?-(diff&1):(diff&1)))>>1;
		odd[k]=even[k]-diff;
	}
	//dwt2_1d_9_7_inv(even, odd, ceilhalf);

	for(int k=0, ks=0;k<floorhalf;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(count&1)
		buffer[stride*(count-1)]=even[floorhalf];
}
void haar_2d_inv(short *buf, int bw, int bh, int nch, int bytestride, int nstages, unsigned char **ret)//buf is destroyed
{
	int lw=floor_log2(bw), lh=floor_log2(bh);
	int *sizes=(int*)malloc(((size_t)(lw<lh?lw:lh)<<1)*sizeof(int));
	if(!sizes)
		return;
	int nsizes=0;
	for(int w2=bw, h2=bh;w2>=3&&h2>=3&&(!nstages||nsizes<nstages);++nsizes)//calculate dimensions of each stage
	{
		sizes[nsizes<<1]=w2;
		sizes[nsizes<<1|1]=h2;
		w2-=w2>>1;//w=ceil(w/2)
		h2-=h2>>1;//h=ceil(h/2)
	}

	int tsize=maximum(bw, bh);
	short *temp=(short*)malloc(tsize*sizeof(short));
	if(!temp)
	{
		free(sizes);
		return;
	}

	int rowlen=bytestride*bw;
	for(int kc=0;kc<nch;++kc)
	{
		//print_buf(buf, bw, bh, kc, bytestride);//

		for(int it=nsizes-1;it>=0;--it)
		{
			int w2=sizes[it<<1], h2=sizes[it<<1|1];

			for(int kx=0;kx<w2;++kx)//vertical IDWT
				haar_1d_inv(buf+bytestride*kx+kc, h2, rowlen, temp);

			//print_buf(buf, bw, bh, kc, bytestride);//

			for(int ky=0;ky<h2;++ky)//horizontal IDWT
				haar_1d_inv(buf+rowlen*ky+kc, w2, bytestride, temp);

			//print_buf(buf, bw, bh, kc, bytestride);//
		}
	}
	free(sizes);
	free(temp);

	size_t dstlen=(size_t)bw*bh*bytestride;
	if(!*ret)
	{
		*ret=(unsigned char*)malloc(dstlen);
		if(!*ret)
			return;
	}
	for(int k=0;k<dstlen;++k)
		ret[0][k]=(unsigned char)buf[k];
}


//integer DCT-II
void intDCTblock8_fwd(int *buf)//buf is aligned packed 32 bit as fixed precision
{
	//https://fgiesen.wordpress.com/2013/11/04/bink-2-2-integer-dct-design-part-1/
	__m128i half=_mm_set1_epi32(1<<5);

	__m128i
		a0lo=_mm_load_si128((__m128i*)buf+ 0), a0hi=_mm_load_si128((__m128i*)buf+ 1),
		a1lo=_mm_load_si128((__m128i*)buf+ 2), a1hi=_mm_load_si128((__m128i*)buf+ 3),
		a2lo=_mm_load_si128((__m128i*)buf+ 4), a2hi=_mm_load_si128((__m128i*)buf+ 5),
		a3lo=_mm_load_si128((__m128i*)buf+ 6), a3hi=_mm_load_si128((__m128i*)buf+ 7),
		a4lo=_mm_load_si128((__m128i*)buf+ 8), a4hi=_mm_load_si128((__m128i*)buf+ 9),
		a5lo=_mm_load_si128((__m128i*)buf+10), a5hi=_mm_load_si128((__m128i*)buf+11),
		a6lo=_mm_load_si128((__m128i*)buf+12), a6hi=_mm_load_si128((__m128i*)buf+13),
		a7lo=_mm_load_si128((__m128i*)buf+14), a7hi=_mm_load_si128((__m128i*)buf+15);

	__m128i
		b0lo=_mm_add_epi32(a0lo, a7lo), b0hi=_mm_add_epi16(a0hi, a7hi),
		b1lo=_mm_add_epi32(a1lo, a6lo), b1hi=_mm_add_epi32(a1hi, a6hi),
		b2lo=_mm_add_epi32(a2lo, a5lo), b2hi=_mm_add_epi32(a2hi, a5hi),
		b3lo=_mm_add_epi32(a3lo, a4lo), b3hi=_mm_add_epi32(a3hi, a4hi),
		b4lo=_mm_sub_epi32(a3lo, a4lo), b4hi=_mm_sub_epi32(a3hi, a4hi),
		b5lo=_mm_sub_epi32(a2lo, a5lo), b5hi=_mm_sub_epi32(a2hi, a5hi),
		b6lo=_mm_sub_epi32(a1lo, a6lo), b6hi=_mm_sub_epi32(a1hi, a6hi),
		b7lo=_mm_sub_epi32(a0lo, a7lo), b7hi=_mm_sub_epi32(a0hi, a7hi);

#define MUL_C7(REG)		_mm_add_epi32(REG, _mm_add_epi32(_mm_slli_epi32(REG, 2), _mm_slli_epi32(REG, 3)))
#define MUL_S7(REG)		_mm_add_epi32(REG, _mm_slli_epi32(REG, 6))
#define MUL_C5(REG)		_mm_add_epi32(REG, _mm_add_epi32(_mm_slli_epi32(REG, 2), _mm_slli_epi32(REG, 5)))
#define MUL_S5(REG)		_mm_sub_epi32(_mm_slli_epi32(REG, 6), _mm_add_epi32(_mm_slli_epi32(REG, 3), REG))
#define MUL_C6(REG)		_mm_sub_epi32(_mm_slli_epi32(REG, 5), _mm_add_epi32(_mm_slli_epi32(REG, 1), REG))
#define MUL_S6(REG)		_mm_add_epi32(_mm_slli_epi32(REG, 6), _mm_add_epi32(_mm_slli_epi32(REG, 2), _mm_slli_epi32(REG, 1)))
#define ROUND(REG)		_mm_srai_epi32(_mm_add_epi32(REG, half), 6)
	a0lo=_mm_add_epi32(b0lo, b3lo); a0hi=_mm_add_epi32(b0hi, b3hi);
	a1lo=_mm_add_epi32(b1lo, b2lo); a1hi=_mm_add_epi32(b1hi, b2hi);
	a2lo=_mm_sub_epi32(b1lo, b2lo); a2hi=_mm_sub_epi32(b1hi, b2hi);
	a3lo=_mm_sub_epi32(b0lo, b3lo); a3hi=_mm_sub_epi32(b0hi, b3hi);
	a4lo=_mm_add_epi32(MUL_C7(b4lo), MUL_S7(b7lo)); a4hi=_mm_add_epi32(MUL_C7(b4hi), MUL_S7(b7hi));//c7*b4+s7*b7
	a5lo=_mm_add_epi32(MUL_C5(b5lo), MUL_S5(b6lo)); a5hi=_mm_add_epi32(MUL_C5(b5hi), MUL_S5(b6hi));//c5*b5+s5*b6
	a6lo=_mm_sub_epi32(MUL_C5(b6lo), MUL_S7(b5lo)); a6hi=_mm_sub_epi32(MUL_C5(b6hi), MUL_S5(b5hi));//c5*b6-s5*b5
	a7lo=_mm_sub_epi32(MUL_C7(b7lo), MUL_S5(b4lo)); a7hi=_mm_sub_epi32(MUL_C7(b7hi), MUL_S7(b4hi));//c7*b7-s7*b4

	b0lo=_mm_add_epi32(a0lo, a1lo); b0hi=_mm_add_epi32(a0hi, a1hi);
	b1lo=_mm_sub_epi32(a0lo, a1lo); b1hi=_mm_sub_epi32(a0hi, a1hi);
	b2lo=_mm_add_epi32(MUL_C6(a2lo), MUL_S6(a3lo)); b2hi=_mm_add_epi32(MUL_C6(a2hi), MUL_S6(a3hi));//c6*a2+s6*a3
	b3lo=_mm_sub_epi32(MUL_C6(a3lo), MUL_S6(a2lo)); b3hi=_mm_sub_epi32(MUL_C6(a3hi), MUL_S6(a2hi));//c6*a3-s6*a2
	b4lo=_mm_add_epi32(a4lo, a5lo); b4hi=_mm_add_epi32(a4hi, a5hi);
	b5lo=_mm_sub_epi32(a4lo, a5lo); b5hi=_mm_sub_epi32(a4hi, a5hi);
	b6lo=_mm_add_epi32(a6lo, a7lo); b6hi=_mm_add_epi32(a6hi, a7hi);
	b7lo=_mm_sub_epi32(a6lo, a7lo); b7hi=_mm_sub_epi32(a6hi, a7hi);

	b2lo=ROUND(b2lo); b2hi=ROUND(b2hi);//divide by 64
	b3lo=ROUND(b3lo); b3hi=ROUND(b3hi);
	b4lo=ROUND(b4lo); b4hi=ROUND(b4hi);
	b5lo=ROUND(b5lo); b5hi=ROUND(b5hi);
	b6lo=ROUND(b6lo); b6hi=ROUND(b6hi);
	b7lo=ROUND(b7lo); b7hi=ROUND(b7hi);

	a5lo=_mm_add_epi32(b5lo, b6lo); a5hi=_mm_add_epi32(b5hi, b6hi);
	a6lo=_mm_sub_epi32(b5lo, b6lo); a6hi=_mm_sub_epi32(b5hi, b6hi);

	_mm_store_si128((__m128i*)buf+ 0, b0lo); _mm_store_si128((__m128i*)buf+ 1, b0hi);
	_mm_store_si128((__m128i*)buf+ 2, b4lo); _mm_store_si128((__m128i*)buf+ 3, b4hi);
	_mm_store_si128((__m128i*)buf+ 4, b2lo); _mm_store_si128((__m128i*)buf+ 5, b2hi);
	_mm_store_si128((__m128i*)buf+ 6, a6lo); _mm_store_si128((__m128i*)buf+ 7, a6hi);
	_mm_store_si128((__m128i*)buf+ 8, b1lo); _mm_store_si128((__m128i*)buf+ 9, b1hi);
	_mm_store_si128((__m128i*)buf+10, b3lo); _mm_store_si128((__m128i*)buf+11, b3hi);
	_mm_store_si128((__m128i*)buf+12, a5lo); _mm_store_si128((__m128i*)buf+13, a5hi);
	_mm_store_si128((__m128i*)buf+14, b7lo); _mm_store_si128((__m128i*)buf+15, b7hi);
#undef	MUL_C7
#undef	MUL_S7
#undef	MUL_C5
#undef	MUL_S5
#undef	MUL_C6
#undef	MUL_S6
}
void transpose_block8(int *buf)
{
	__m128
		a0lo=_mm_load_ps((float*)buf+ 0*4), a0hi=_mm_load_ps((float*)buf+ 1*4),
		a1lo=_mm_load_ps((float*)buf+ 2*4), a1hi=_mm_load_ps((float*)buf+ 3*4),
		a2lo=_mm_load_ps((float*)buf+ 4*4), a2hi=_mm_load_ps((float*)buf+ 5*4),
		a3lo=_mm_load_ps((float*)buf+ 6*4), a3hi=_mm_load_ps((float*)buf+ 7*4),
		a4lo=_mm_load_ps((float*)buf+ 8*4), a4hi=_mm_load_ps((float*)buf+ 9*4),
		a5lo=_mm_load_ps((float*)buf+10*4), a5hi=_mm_load_ps((float*)buf+11*4),
		a6lo=_mm_load_ps((float*)buf+12*4), a6hi=_mm_load_ps((float*)buf+13*4),
		a7lo=_mm_load_ps((float*)buf+14*4), a7hi=_mm_load_ps((float*)buf+15*4);

	_MM_TRANSPOSE4_PS(a0lo, a1lo, a2lo, a3lo);
	_MM_TRANSPOSE4_PS(a0hi, a1hi, a2hi, a3hi);
	_MM_TRANSPOSE4_PS(a4lo, a5lo, a6lo, a7lo);
	_MM_TRANSPOSE4_PS(a4hi, a5hi, a6hi, a7hi);
	
	_mm_store_ps((float*)buf+ 0*4, a0lo); _mm_store_ps((float*)buf+ 1*4, a4lo);
	_mm_store_ps((float*)buf+ 2*4, a1lo); _mm_store_ps((float*)buf+ 3*4, a5lo);
	_mm_store_ps((float*)buf+ 4*4, a2lo); _mm_store_ps((float*)buf+ 5*4, a6lo);
	_mm_store_ps((float*)buf+ 6*4, a3lo); _mm_store_ps((float*)buf+ 7*4, a7lo);
	_mm_store_ps((float*)buf+ 8*4, a0hi); _mm_store_ps((float*)buf+ 9*4, a4hi);
	_mm_store_ps((float*)buf+10*4, a1hi); _mm_store_ps((float*)buf+11*4, a5hi);
	_mm_store_ps((float*)buf+12*4, a2hi); _mm_store_ps((float*)buf+13*4, a6hi);
	_mm_store_ps((float*)buf+14*4, a3hi); _mm_store_ps((float*)buf+15*4, a7hi);
}
void intDCT8x8_fwd(const unsigned char *buf, int bw, int bh, int nch, int bytestride, unsigned short **ret)
{
	if(bw&7||bh&7)
		return;
	size_t len=(size_t)bw*bh*bytestride;
	if(!*ret)
	{
		*ret=(unsigned short*)malloc(len*sizeof(short));
		if(!*ret)
			return;
	}
	ALIGN(16) int temp[64];
	for(int ky=0;ky<bh;ky+=8)
	{
		for(int kx=0;kx<bw;kx+=8)
		{
			for(int kc=0;kc<nch;++kc)
			{
				for(int ky2=0;ky2<8;++ky2)
				{
					for(int kx2=0;kx2<8;++kx2)
						temp[ky2<<3|kx2]=buf[bytestride*(bw*(ky|ky2)+(kx|kx2))+kc];
					intDCTblock8_fwd(temp);
					transpose_block8(temp);
					intDCTblock8_fwd(temp);
					transpose_block8(temp);
					for(int kx2=0;kx2<8;++kx2)
						ret[0][bytestride*(bw*(ky|ky2)+(kx|kx2))+kc]=temp[ky2<<3|kx2];
				}
			}
		}
	}
}
#if 1
void DCTtest()
{
	ALIGN(16) int block[64];
	int b2[64];
	long long result[64];
	for(int ky=0;ky<8;++ky)
	{
		for(int kx=0;kx<8;++kx)
			block[ky<<3|kx]=ky==kx;
	}
	intDCTblock8_fwd(block);
	memcpy(b2, block, 64*sizeof(int));
	transpose_block8(block);

	for(int ky=0;ky<8;++ky)//matmul
	{
		for(int kx=0;kx<8;++kx)
		{
			long long sum=0;
			for(int ki=0;ki<8;++ki)
				sum+=(long long)b2[ky<<3|ki]*block[ki<<3|kx];
			result[ky<<3|kx]=sum;
		}
	}
	for(int ky=0;ky<8;++ky)
	{
		for(int kx=0;kx<8;++kx)
			printf("%lld\t", result[ky<<3|kx]);
		printf("\n");
	}
}
#endif