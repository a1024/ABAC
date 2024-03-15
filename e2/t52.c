#include"e2.h"
#include"ac.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctype.h>
#ifdef _MSC_VER
#include<intrin.h>
#else
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;

#define FIXBITS 12
#define ROUNDCONST ((1<<FIXBITS)>>1)
#define FLOAT2FIX(X) (int)((X)*(1<<FIXBITS)+0.5)
void ct_YCbCr_lossy(Image *src, int fwd)
{
	if(fwd)
	{
		for(ptrdiff_t k=0, res=(ptrdiff_t)src->iw*src->ih;k<res;++k)
		{
			int r=src->data[k<<2|0], g=src->data[k<<2|1], b=src->data[k<<2|2];

			int
				Y=(FLOAT2FIX(0.299)*r+FLOAT2FIX(0.587)*g+FLOAT2FIX(0.114)*b+ROUNDCONST)>>FIXBITS,
				Cb=(-FLOAT2FIX(0.168736)*r-FLOAT2FIX(0.331264)*g+FLOAT2FIX(0.5)*b+ROUNDCONST)>>FIXBITS,
				Cr=(FLOAT2FIX(0.5)*r-FLOAT2FIX(0.418688)*g-FLOAT2FIX(0.081312)*b+ROUNDCONST)>>FIXBITS;

			src->data[k<<2|0]=Y;
			src->data[k<<2|1]=Cb;
			src->data[k<<2|2]=Cr;
		}
	}
	else
	{
		for(ptrdiff_t k=0, res=(ptrdiff_t)src->iw*src->ih;k<res;++k)
		{
			int Y=src->data[k<<2|0]<<FIXBITS|ROUNDCONST, Cb=src->data[k<<2|1], Cr=src->data[k<<2|2];
			
			int
				r=(Y+FLOAT2FIX(1.402)*Cr+ROUNDCONST)>>FIXBITS,
				g=(Y-FLOAT2FIX(0.344136)*Cb-FLOAT2FIX(0.714136)*Cr+ROUNDCONST)>>FIXBITS,
				b=(Y+FLOAT2FIX(1.772)*Cb+ROUNDCONST)>>FIXBITS;
			
			src->data[k<<2|0]=r;
			src->data[k<<2|1]=g;
			src->data[k<<2|2]=b;
		}
	}
}

void resample_YUV420(Image *src, int down)
{
	if(down)
	{
		for(int kc=1;kc<3;++kc)
		{
			for(int ky=0;ky<src->ih;ky+=2)
			{
				for(int kx=0;kx<src->iw;kx+=2)
				{
					int idx=src->iw*ky+kx;
					int v0=src->data[idx<<2|kc];
					int v1=src->data[(idx+1)<<2|kc];
					int v2=src->data[(idx+src->iw)<<2|kc];
					int v3=src->data[(idx+src->iw+1)<<2|kc];
					idx=src->iw*(ky>>1)+(kx>>1);
					src->data[idx<<2|kc]=(v0+v1+v2+v3)>>2;
				}
			}
			//for(int ky=0;ky<src->ih;ky+=2)
			//{
			//	for(int kx=0;kx<src->iw;kx+=2)
			//	{
			//		int idx=src->iw*ky+kx;
			//		int v0=src->data[idx<<2|kc];
			//		int v1=src->data[(idx+1)<<2|kc];
			//		int v2=src->data[(idx+src->iw)<<2|kc];
			//		int v3=src->data[(idx+src->iw+1)<<2|kc];
			//		src->data[idx<<2|kc]=(v0+v1+v2+v3)>>2;
			//	}
			//}
		}
	}
	else
	{
		for(int kc=1;kc<3;++kc)
		{
			for(int ky=src->ih-2;ky>=0;ky-=2)
			{
				for(int kx=src->iw-2;kx>=0;kx-=2)
				{
					int idx=src->iw*(ky>>1)+(kx>>1);
					int v0=src->data[idx<<2|kc];
					idx=src->iw*ky+kx;
					src->data[idx<<2|kc]=v0;
					src->data[(idx+1)<<2|kc]=v0;
					src->data[(idx+src->iw)<<2|kc]=v0;
					src->data[(idx+src->iw+1)<<2|kc]=v0;
				}
			}
			//for(int ky=0;ky<src->ih;ky+=2)
			//{
			//	for(int kx=0;kx<src->iw;kx+=2)
			//	{
			//		int idx=src->iw*ky+kx;
			//		int v0=src->data[idx<<2|kc];
			//		src->data[(idx+1)<<2|kc]=v0;
			//		src->data[(idx+src->iw)<<2|kc]=v0;
			//		src->data[(idx+src->iw+1)<<2|kc]=v0;
			//	}
			//}
		}
	}
}

static void block8_dct_lossy_fwd(int *data)//https://github.com/rockcarry/ffjpeg
{
	__m256i c45=_mm256_set1_epi32(FLOAT2FIX(0.707106781));//cosd45 = 1/sqrt2
	__m256i c675=_mm256_set1_epi32(FLOAT2FIX(0.382683433));//cosd67.5
	__m256i c2=_mm256_set1_epi32(FLOAT2FIX(0.541196100));//sqrt(1-1/sqrt2)
	__m256i c3=_mm256_set1_epi32(FLOAT2FIX(1.306562965));//sqrt(1+1/sqrt2)

	__m256i r0=_mm256_load_si256((__m256i*)data+0);
	__m256i r1=_mm256_load_si256((__m256i*)data+1);
	__m256i r2=_mm256_load_si256((__m256i*)data+2);
	__m256i r3=_mm256_load_si256((__m256i*)data+3);
	__m256i r4=_mm256_load_si256((__m256i*)data+4);
	__m256i r5=_mm256_load_si256((__m256i*)data+5);
	__m256i r6=_mm256_load_si256((__m256i*)data+6);
	__m256i r7=_mm256_load_si256((__m256i*)data+7);

	__m256i t0=_mm256_add_epi32(r0, r7);
	__m256i t7=_mm256_sub_epi32(r0, r7);
	__m256i t1=_mm256_add_epi32(r1, r6);
	__m256i t6=_mm256_sub_epi32(r1, r6);
	__m256i t2=_mm256_add_epi32(r2, r5);
	__m256i t5=_mm256_sub_epi32(r2, r5);
	__m256i t3=_mm256_add_epi32(r3, r4);
	__m256i t4=_mm256_sub_epi32(r3, r4);

	r0=_mm256_add_epi32(t0, t3);
	r1=_mm256_sub_epi32(t0, t3);
	r2=_mm256_add_epi32(t1, t2);
	r3=_mm256_sub_epi32(t1, t2);
	r4=_mm256_add_epi32(t4, t5);
	r5=_mm256_add_epi32(t5, t6);
	r6=_mm256_add_epi32(t6, t7);

	__m256i d0=_mm256_add_epi32(r1, r3);
	__m256i d1=_mm256_add_epi32(r4, r6);

	__m256i p0=_mm256_mullo_epi32(d0, c45);	//z1
	__m256i p1=_mm256_mullo_epi32(d1, c675);//z5
	__m256i p2=_mm256_mullo_epi32(r4, c2);	//z2
	__m256i p3=_mm256_mullo_epi32(r6, c3);	//z4
	__m256i p4=_mm256_mullo_epi32(r5, c45);	//z3

	p0=_mm256_srai_epi32(p0, FIXBITS);
	p1=_mm256_srai_epi32(p1, FIXBITS);
	p2=_mm256_srai_epi32(p2, FIXBITS);
	p3=_mm256_srai_epi32(p3, FIXBITS);
	p4=_mm256_srai_epi32(p4, FIXBITS);

	__m256i s0=_mm256_add_epi32(p2, p1);//z2
	__m256i s1=_mm256_add_epi32(p3, p1);//z4
	__m256i s2=_mm256_add_epi32(t7, p4);//z11
	__m256i s3=_mm256_sub_epi32(t7, p4);//z13

	__m256i f0=_mm256_add_epi32(r0, r2);
	__m256i f4=_mm256_sub_epi32(r0, r2);
	__m256i f2=_mm256_add_epi32(r1, p0);
	__m256i f6=_mm256_sub_epi32(r1, p0);
	__m256i f5=_mm256_add_epi32(s0, s3);
	__m256i f3=_mm256_sub_epi32(s0, s3);
	__m256i f1=_mm256_add_epi32(s2, s1);
	__m256i f7=_mm256_sub_epi32(s2, s1);

	_mm256_store_si256((__m256i*)data+0, f0);
	_mm256_store_si256((__m256i*)data+1, f1);
	_mm256_store_si256((__m256i*)data+2, f2);
	_mm256_store_si256((__m256i*)data+3, f3);
	_mm256_store_si256((__m256i*)data+4, f4);
	_mm256_store_si256((__m256i*)data+5, f5);
	_mm256_store_si256((__m256i*)data+6, f6);
	_mm256_store_si256((__m256i*)data+7, f7);
}
static void block8_dct_lossy_inv(int *data)
{
	__m256i sqrt2=_mm256_set1_epi32(FLOAT2FIX(1.414213562));
	__m256i c2=_mm256_set1_epi32(FLOAT2FIX(1.847759065));
	__m256i c3=_mm256_set1_epi32(FLOAT2FIX(1.082392200));
	__m256i c4=_mm256_set1_epi32(FLOAT2FIX(-2.613125930));

	__m256i r0=_mm256_load_si256((__m256i*)data+0);
	__m256i r4=_mm256_load_si256((__m256i*)data+1);
	__m256i r1=_mm256_load_si256((__m256i*)data+2);
	__m256i r5=_mm256_load_si256((__m256i*)data+3);
	__m256i r2=_mm256_load_si256((__m256i*)data+4);
	__m256i r6=_mm256_load_si256((__m256i*)data+5);
	__m256i r3=_mm256_load_si256((__m256i*)data+6);
	__m256i r7=_mm256_load_si256((__m256i*)data+7);

	__m256i tmp10=_mm256_add_epi32(r0, r2);
	__m256i tmp11=_mm256_sub_epi32(r0, r2);
	__m256i tmp13=_mm256_add_epi32(r1, r3);
	__m256i tmp12=_mm256_sub_epi32(r1, r3);
	tmp12=_mm256_mullo_epi32(tmp12, sqrt2);
	tmp12=_mm256_srai_epi32(tmp12, FIXBITS);
	tmp12=_mm256_sub_epi32(tmp12, tmp13);

	__m256i tmp0=_mm256_add_epi32(tmp10, tmp13);
	__m256i tmp3=_mm256_sub_epi32(tmp10, tmp13);
	__m256i tmp1=_mm256_add_epi32(tmp11, tmp12);
	__m256i tmp2=_mm256_sub_epi32(tmp11, tmp12);

	__m256i z13=_mm256_add_epi32(r6, r5);
	__m256i z10=_mm256_sub_epi32(r6, r5);
	__m256i z11=_mm256_add_epi32(r4, r7);
	__m256i z12=_mm256_sub_epi32(r4, r7);

	__m256i tmp7=_mm256_add_epi32(z11, z13);
	tmp11=_mm256_sub_epi32(z11, z13);
	tmp11=_mm256_mullo_epi32(tmp11, sqrt2);
	tmp11=_mm256_srai_epi32(tmp11, FIXBITS);

	__m256i z5=_mm256_add_epi32(z10, z12);
	z5=_mm256_mullo_epi32(z5, c2);
	z5=_mm256_srai_epi32(z5, FIXBITS);

	tmp10=_mm256_mullo_epi32(z12, c3);
	tmp12=_mm256_mullo_epi32(z10, c4);
	tmp10=_mm256_srai_epi32(tmp10, FIXBITS);
	tmp12=_mm256_srai_epi32(tmp12, FIXBITS);

	tmp10=_mm256_sub_epi32(tmp10, z5);
	tmp12=_mm256_add_epi32(tmp12, z5);

	__m256i tmp6=_mm256_sub_epi32(tmp12, tmp7);
	__m256i tmp5=_mm256_sub_epi32(tmp11, tmp6);
	__m256i tmp4=_mm256_add_epi32(tmp10, tmp5);

	r0=_mm256_add_epi32(tmp0, tmp7);
	r7=_mm256_sub_epi32(tmp0, tmp7);
	r1=_mm256_add_epi32(tmp1, tmp6);
	r6=_mm256_sub_epi32(tmp1, tmp6);
	r2=_mm256_add_epi32(tmp2, tmp5);
	r5=_mm256_sub_epi32(tmp2, tmp5);
	r4=_mm256_add_epi32(tmp3, tmp4);
	r3=_mm256_sub_epi32(tmp3, tmp4);

	_mm256_store_si256((__m256i*)data+0, r0);
	_mm256_store_si256((__m256i*)data+1, r1);
	_mm256_store_si256((__m256i*)data+2, r2);
	_mm256_store_si256((__m256i*)data+3, r3);
	_mm256_store_si256((__m256i*)data+4, r4);
	_mm256_store_si256((__m256i*)data+5, r5);
	_mm256_store_si256((__m256i*)data+6, r6);
	_mm256_store_si256((__m256i*)data+7, r7);
}

//#define BUTTERFLY_FWD(A, B, TEMP) TEMP=A, A=_mm256_add_epi32(A, B), B=_mm256_sub_epi32(TEMP, B)
//#define BUTTERFLY_INV(A, B, TEMP) TEMP=A, A=_mm256_add_epi32(A, B), B=_mm256_sub_epi32(TEMP, B)
#define BUTTERFLY_FWD(A, B, TEMP) B=_mm256_sub_epi32(B, A), A=_mm256_add_epi32(A, _mm256_srai_epi32(B, 1))
#define BUTTERFLY_INV(A, B, TEMP) A=_mm256_sub_epi32(A, _mm256_srai_epi32(B, 1)), B=_mm256_add_epi32(B, A)
#define HALF(X) _mm256_srai_epi32(X, 1)
#define MUL15(X, SR) _mm256_srai_epi32(_mm256_sub_epi32(_mm256_slli_epi32(X, 4), X), SR)
#define MUL13(X, SR) _mm256_srai_epi32(_mm256_add_epi32(_mm256_add_epi32(X, _mm256_slli_epi32(X, 2)), _mm256_slli_epi32(X, 3)), SR)
#define MUL11(X, SR) _mm256_srai_epi32(_mm256_add_epi32(_mm256_add_epi32(X, _mm256_slli_epi32(X, 1)), _mm256_slli_epi32(X, 3)), SR)
#define MUL03(X, SR) _mm256_srai_epi32(_mm256_add_epi32(_mm256_slli_epi32(X, 1), X), SR)
static void block8_dct_ls_inv(int *data);
static void block8_dct_ls_fwd(int *data)
{
	//approximate lossless DCT
	//https://fgiesen.wordpress.com/2013/11/04/bink-2-2-integer-dct-design-part-1/
	//http://thanglong.ece.jhu.edu/Tran/Pub/binDCT.pdf
	//binDCT-C1
	__m256i r0=_mm256_load_si256((__m256i*)data+0);
	__m256i r1=_mm256_load_si256((__m256i*)data+1);
	__m256i r2=_mm256_load_si256((__m256i*)data+2);
	__m256i r3=_mm256_load_si256((__m256i*)data+3);
	__m256i r4=_mm256_load_si256((__m256i*)data+4);
	__m256i r5=_mm256_load_si256((__m256i*)data+5);
	__m256i r6=_mm256_load_si256((__m256i*)data+6);
	__m256i r7=_mm256_load_si256((__m256i*)data+7);
	__m256i temp;

	BUTTERFLY_FWD(r0, r7, temp);
	BUTTERFLY_FWD(r1, r6, temp);
	BUTTERFLY_FWD(r2, r5, temp);
	BUTTERFLY_FWD(r3, r4, temp);

	r5=_mm256_sub_epi32(r5, MUL13(r6, 5));//p4 = p1 ~= 13/32
	r6=_mm256_add_epi32(r6, MUL11(r5, 4));//u4 = 0.7071067811 ~= 11/16
	r5=_mm256_sub_epi32(MUL13(r6, 5), r5);//p5 = p1 ~= 13/32

	BUTTERFLY_FWD(r0, r3, temp);
	BUTTERFLY_FWD(r1, r2, temp);
	BUTTERFLY_FWD(r4, r5, temp);
	BUTTERFLY_FWD(r7, r6, temp);
	
	r0=_mm256_add_epi32(r0, r1);
	r1=_mm256_sub_epi32(HALF(r0), r1);
	r2=_mm256_sub_epi32(MUL13(r3, 5), r2);//p1 = 0.4142135623 ~= 13/32
	r3=_mm256_sub_epi32(r3, MUL11(r2, 5));//u1 = 0.3535533905 ~= 11/32
	r5=_mm256_add_epi32(r5, MUL11(r6, 4));//p2 = 0.6681786379 ~= 11/16
	r6=_mm256_sub_epi32(r6, MUL15(r5, 5));//u2 = 0.4619397662 ~= 15/32
	r4=_mm256_sub_epi32(MUL03(r7, 4), r4);//p3 = 0.1989123673 ~= 3/16
	r7=_mm256_sub_epi32(r7, MUL03(r4, 4));//u3 = 0.1913417161 ~= 3/16
	
#if 0
	__m256i b0=r0;
	__m256i b1=r1;
	__m256i b2=r2;
	__m256i b3=r3;
	__m256i b4=r4;
	__m256i b5=r5;
	__m256i b6=r6;
	__m256i b7=r7;
	r7=_mm256_add_epi32(r7, MUL03(r4, 4));//u3 = 0.1913417161 ~= 3/16
	r4=_mm256_sub_epi32(MUL03(r7, 4), r4);//p3 = 0.1989123673 ~= 3/16
	r6=_mm256_add_epi32(r6, MUL15(r5, 5));//u2 = 0.4619397662 ~= 15/32
	r5=_mm256_sub_epi32(r5, MUL11(r6, 4));//p2 = 0.6681786379 ~= 11/16
	r3=_mm256_add_epi32(r3, MUL11(r2, 5));//u1 = 0.3535533905 ~= 11/32
	r2=_mm256_sub_epi32(MUL13(r3, 5), r2);//p1 = 0.4142135623 ~= 13/32
	r1=_mm256_sub_epi32(HALF(r0), r1);
	r0=_mm256_sub_epi32(r0, r1);
	
	BUTTERFLY_INV(r0, r3, temp);
	BUTTERFLY_INV(r1, r2, temp);
	BUTTERFLY_INV(r4, r5, temp);
	BUTTERFLY_INV(r7, r6, temp);
	
	r5=_mm256_sub_epi32(MUL13(r6, 5), r5);//p5 = p1 ~= 13/32
	r6=_mm256_sub_epi32(r6, MUL11(r5, 4));//u4 = 0.7071067811 ~= 11/16
	r5=_mm256_add_epi32(r5, MUL13(r6, 5));//p4 = p1 ~= 13/32

	BUTTERFLY_INV(r0, r7, temp);
	BUTTERFLY_INV(r1, r6, temp);
	BUTTERFLY_INV(r2, r5, temp);
	BUTTERFLY_INV(r3, r4, temp);

	//r0=_mm256_srai_epi32(r0, 2);
	//r1=_mm256_srai_epi32(r1, 2);
	//r2=_mm256_srai_epi32(r2, 2);
	//r3=_mm256_srai_epi32(r3, 2);
	//r4=_mm256_srai_epi32(r4, 2);
	//r5=_mm256_srai_epi32(r5, 2);
	//r6=_mm256_srai_epi32(r6, 2);
	//r7=_mm256_srai_epi32(r7, 2);
	
	int d2[64];
	_mm256_store_si256((__m256i*)d2+0, r0);
	_mm256_store_si256((__m256i*)d2+1, r1);
	_mm256_store_si256((__m256i*)d2+2, r2);
	_mm256_store_si256((__m256i*)d2+3, r3);
	_mm256_store_si256((__m256i*)d2+4, r4);
	_mm256_store_si256((__m256i*)d2+5, r5);
	_mm256_store_si256((__m256i*)d2+6, r6);
	_mm256_store_si256((__m256i*)d2+7, r7);
	//block8_dct_ls_inv(d2);
	if(memcmp(d2, data, sizeof(d2)))
		LOG_ERROR("");
	_mm256_store_si256((__m256i*)data+0, b0);
	_mm256_store_si256((__m256i*)data+1, b7);
	_mm256_store_si256((__m256i*)data+2, b3);
	_mm256_store_si256((__m256i*)data+3, b6);
	_mm256_store_si256((__m256i*)data+4, b1);
	_mm256_store_si256((__m256i*)data+5, b5);
	_mm256_store_si256((__m256i*)data+6, b2);
	_mm256_store_si256((__m256i*)data+7, b4);
#else
	_mm256_store_si256((__m256i*)data+0, r0);
	_mm256_store_si256((__m256i*)data+1, r7);
	_mm256_store_si256((__m256i*)data+2, r3);
	_mm256_store_si256((__m256i*)data+3, r6);
	_mm256_store_si256((__m256i*)data+4, r1);
	_mm256_store_si256((__m256i*)data+5, r5);
	_mm256_store_si256((__m256i*)data+6, r2);
	_mm256_store_si256((__m256i*)data+7, r4);
#endif
#if 0
	r7=_mm256_sub_epi32(r0, r7);
	r6=_mm256_sub_epi32(r1, r6);
	r5=_mm256_sub_epi32(r2, r5);
	r4=_mm256_sub_epi32(r3, r4);
	r0=_mm256_sub_epi32(r0, _mm256_srai_epi32(r7, 1));
	r1=_mm256_sub_epi32(r1, _mm256_srai_epi32(r6, 1));
	r2=_mm256_sub_epi32(r2, _mm256_srai_epi32(r5, 1));
	r3=_mm256_sub_epi32(r3, _mm256_srai_epi32(r4, 1));

	r3=_mm256_sub_epi32(r0, r3);
	r2=_mm256_sub_epi32(r1, r2);
	r0=_mm256_sub_epi32(r0, _mm256_srai_epi32(r3, 1));
	r1=_mm256_sub_epi32(r1, _mm256_srai_epi32(r2, 1));

	r1=_mm256_sub_epi32(r0, r1);
	r0=_mm256_sub_epi32(r0, _mm256_srai_epi32(r1, 1));
	r2=_mm256_sub_epi32(MUL13(r3), r2);
	r3=_mm256_sub_epi32(r3, MUL11(r2));

	r5=_mm256_sub_epi32(r5, MUL13(r6));
	r6=_mm256_add_epi32(r6, MUL11(r5));
	r5=_mm256_sub_epi32(MUL15(r6), r5);

	r5=_mm256_sub_epi32(r4, r5);
	r6=_mm256_sub_epi32(r7, r6);
	r4=_mm256_sub_epi32(r4, _mm256_srai_epi32(r5, 1));
	r7=_mm256_sub_epi32(r7, _mm256_srai_epi32(r6, 1));

	r4=_mm256_sub_epi32(MUL03(r7), r4);
	r7=_mm256_sub_epi32(r7, MUL03(r4));

	r5=_mm256_add_epi32(r5, MUL11(r6));
	r6=_mm256_sub_epi32(r6, MUL15(r5));

	_mm256_store_si256((__m256i*)data+0, r0);
	_mm256_store_si256((__m256i*)data+0, r0);
#endif
}
static void block8_dct_ls_inv(int *data)
{
	__m256i r0=_mm256_load_si256((__m256i*)data+0);
	__m256i r7=_mm256_load_si256((__m256i*)data+1);
	__m256i r3=_mm256_load_si256((__m256i*)data+2);
	__m256i r6=_mm256_load_si256((__m256i*)data+3);
	__m256i r1=_mm256_load_si256((__m256i*)data+4);
	__m256i r5=_mm256_load_si256((__m256i*)data+5);
	__m256i r2=_mm256_load_si256((__m256i*)data+6);
	__m256i r4=_mm256_load_si256((__m256i*)data+7);
	__m256i temp;

	r7=_mm256_add_epi32(r7, MUL03(r4, 4));//u3 = 0.1913417161 ~= 3/16
	r4=_mm256_sub_epi32(MUL03(r7, 4), r4);//p3 = 0.1989123673 ~= 3/16
	r6=_mm256_add_epi32(r6, MUL15(r5, 5));//u2 = 0.4619397662 ~= 15/32
	r5=_mm256_sub_epi32(r5, MUL11(r6, 4));//p2 = 0.6681786379 ~= 11/16
	r3=_mm256_add_epi32(r3, MUL11(r2, 5));//u1 = 0.3535533905 ~= 11/32
	r2=_mm256_sub_epi32(MUL13(r3, 5), r2);//p1 = 0.4142135623 ~= 13/32
	r1=_mm256_sub_epi32(HALF(r0), r1);
	r0=_mm256_sub_epi32(r0, r1);
	
	BUTTERFLY_INV(r0, r3, temp);
	BUTTERFLY_INV(r1, r2, temp);
	BUTTERFLY_INV(r4, r5, temp);
	BUTTERFLY_INV(r7, r6, temp);
	
	r5=_mm256_sub_epi32(MUL13(r6, 5), r5);//p5 = p1 ~= 13/32
	r6=_mm256_sub_epi32(r6, MUL11(r5, 4));//u4 = 0.7071067811 ~= 11/16
	r5=_mm256_add_epi32(r5, MUL13(r6, 5));//p4 = p1 ~= 13/32

	BUTTERFLY_INV(r0, r7, temp);
	BUTTERFLY_INV(r1, r6, temp);
	BUTTERFLY_INV(r2, r5, temp);
	BUTTERFLY_INV(r3, r4, temp);

	//r0=_mm256_srai_epi32(r0, 2);
	//r1=_mm256_srai_epi32(r1, 2);
	//r2=_mm256_srai_epi32(r2, 2);
	//r3=_mm256_srai_epi32(r3, 2);
	//r4=_mm256_srai_epi32(r4, 2);
	//r5=_mm256_srai_epi32(r5, 2);
	//r6=_mm256_srai_epi32(r6, 2);
	//r7=_mm256_srai_epi32(r7, 2);
	
	_mm256_store_si256((__m256i*)data+0, r0);
	_mm256_store_si256((__m256i*)data+1, r1);
	_mm256_store_si256((__m256i*)data+2, r2);
	_mm256_store_si256((__m256i*)data+3, r3);
	_mm256_store_si256((__m256i*)data+4, r4);
	_mm256_store_si256((__m256i*)data+5, r5);
	_mm256_store_si256((__m256i*)data+6, r6);
	_mm256_store_si256((__m256i*)data+7, r7);
}
#undef  MUL15
#undef  MUL13
#undef  MUL11
#undef  MUL03

static void block8_transpose(int *data)//https://stackoverflow.com/questions/73977998/simd-transposition-of-8x8-matrix-of-32-bit-values-in-java
{
	__m256i a[8], b[8], c[8];
	a[0]=_mm256_load_si256((__m256i*)data+0);
	a[1]=_mm256_load_si256((__m256i*)data+1);
	a[2]=_mm256_load_si256((__m256i*)data+2);
	a[3]=_mm256_load_si256((__m256i*)data+3);
	a[4]=_mm256_load_si256((__m256i*)data+4);
	a[5]=_mm256_load_si256((__m256i*)data+5);
	a[6]=_mm256_load_si256((__m256i*)data+6);
	a[7]=_mm256_load_si256((__m256i*)data+7);

	b[0]=_mm256_unpacklo_epi32(a[0], a[1]);
	b[1]=_mm256_unpackhi_epi32(a[0], a[1]);
	b[2]=_mm256_unpacklo_epi32(a[2], a[3]);
	b[3]=_mm256_unpackhi_epi32(a[2], a[3]);
	b[4]=_mm256_unpacklo_epi32(a[4], a[5]);
	b[5]=_mm256_unpackhi_epi32(a[4], a[5]);
	b[6]=_mm256_unpacklo_epi32(a[6], a[7]);
	b[7]=_mm256_unpackhi_epi32(a[6], a[7]);

	c[0]=_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(b[0]), _mm256_castsi256_ps(b[2]), _MM_SHUFFLE(1, 0, 1, 0)));
	c[1]=_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(b[0]), _mm256_castsi256_ps(b[2]), _MM_SHUFFLE(3, 2, 3, 2)));
	c[2]=_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(b[1]), _mm256_castsi256_ps(b[3]), _MM_SHUFFLE(1, 0, 1, 0)));
	c[3]=_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(b[1]), _mm256_castsi256_ps(b[3]), _MM_SHUFFLE(3, 2, 3, 2)));
	c[4]=_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(b[4]), _mm256_castsi256_ps(b[6]), _MM_SHUFFLE(1, 0, 1, 0)));
	c[5]=_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(b[4]), _mm256_castsi256_ps(b[6]), _MM_SHUFFLE(3, 2, 3, 2)));
	c[6]=_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(b[5]), _mm256_castsi256_ps(b[7]), _MM_SHUFFLE(1, 0, 1, 0)));
	c[7]=_mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(b[5]), _mm256_castsi256_ps(b[7]), _MM_SHUFFLE(3, 2, 3, 2)));

	a[0]=_mm256_permute2f128_si256(c[0], c[4], _MM_SHUFFLE(0, 2, 0, 0));
	a[1]=_mm256_permute2f128_si256(c[1], c[5], _MM_SHUFFLE(0, 2, 0, 0));
	a[2]=_mm256_permute2f128_si256(c[2], c[6], _MM_SHUFFLE(0, 2, 0, 0));
	a[3]=_mm256_permute2f128_si256(c[3], c[7], _MM_SHUFFLE(0, 2, 0, 0));
	a[4]=_mm256_permute2f128_si256(c[0], c[4], _MM_SHUFFLE(0, 3, 0, 1));
	a[5]=_mm256_permute2f128_si256(c[1], c[5], _MM_SHUFFLE(0, 3, 0, 1));
	a[6]=_mm256_permute2f128_si256(c[2], c[6], _MM_SHUFFLE(0, 3, 0, 1));
	a[7]=_mm256_permute2f128_si256(c[3], c[7], _MM_SHUFFLE(0, 3, 0, 1));

	_mm256_store_si256((__m256i*)data+0, a[0]);
	_mm256_store_si256((__m256i*)data+1, a[1]);
	_mm256_store_si256((__m256i*)data+2, a[2]);
	_mm256_store_si256((__m256i*)data+3, a[3]);
	_mm256_store_si256((__m256i*)data+4, a[4]);
	_mm256_store_si256((__m256i*)data+5, a[5]);
	_mm256_store_si256((__m256i*)data+6, a[6]);
	_mm256_store_si256((__m256i*)data+7, a[7]);
}
static void block8_sra3(int *data)
{
	__m256i a[8];
	a[0]=_mm256_load_si256((__m256i*)data+0);
	a[1]=_mm256_load_si256((__m256i*)data+1);
	a[2]=_mm256_load_si256((__m256i*)data+2);
	a[3]=_mm256_load_si256((__m256i*)data+3);
	a[4]=_mm256_load_si256((__m256i*)data+4);
	a[5]=_mm256_load_si256((__m256i*)data+5);
	a[6]=_mm256_load_si256((__m256i*)data+6);
	a[7]=_mm256_load_si256((__m256i*)data+7);

	a[0]=_mm256_srai_epi32(a[0], 3);
	a[1]=_mm256_srai_epi32(a[1], 3);
	a[2]=_mm256_srai_epi32(a[2], 3);
	a[3]=_mm256_srai_epi32(a[3], 3);
	a[4]=_mm256_srai_epi32(a[4], 3);
	a[5]=_mm256_srai_epi32(a[5], 3);
	a[6]=_mm256_srai_epi32(a[6], 3);
	a[7]=_mm256_srai_epi32(a[7], 3);
	
	_mm256_store_si256((__m256i*)data+0, a[0]);
	_mm256_store_si256((__m256i*)data+1, a[1]);
	_mm256_store_si256((__m256i*)data+2, a[2]);
	_mm256_store_si256((__m256i*)data+3, a[3]);
	_mm256_store_si256((__m256i*)data+4, a[4]);
	_mm256_store_si256((__m256i*)data+5, a[5]);
	_mm256_store_si256((__m256i*)data+6, a[6]);
	_mm256_store_si256((__m256i*)data+7, a[7]);
}
void dct_8x8(Image *dst, Image const *src, int subsampled, int *ymatrix, int *uvmatrix, int fwd)
{
	ALIGN(16) int data[64];
	int iw[]=
	{
		src->iw,
		src->iw>>subsampled,
		src->iw>>subsampled,
		src->iw,
	};
	int ih[]=
	{
		src->ih,
		src->ih>>subsampled,
		src->ih>>subsampled,
		src->ih,
	};
	for(int kc=0;kc<src->nch;++kc)
	{
		int *matrix=kc?uvmatrix:ymatrix;
		for(int ky=0;ky<ih[kc]-7;ky+=8)
		{
			for(int kx=0;kx<iw[kc]-7;kx+=8)
			{
				if(fwd)
				{
					for(int ky2=0;ky2<8;++ky2)
					{
						for(int kx2=0;kx2<8;++kx2)
							data[ky2<<3|kx2]=src->data[(src->iw*(ky+ky2)+kx+kx2)<<2|kc];
					}
					//int d2[64];//
					//memcpy(d2, data, sizeof(d2));//

					block8_dct_ls_fwd(data);
					block8_transpose(data);
					block8_dct_ls_fwd(data);
					if(matrix)
					{
						for(int k=0;k<64;++k)
							data[k]/=matrix[k];
					}

					//block8_dct_ls_inv(data);//
					//block8_transpose(data);//
					//block8_dct_ls_inv(data);//
					//if(memcmp(d2, data, sizeof(d2)))//
					//	LOG_ERROR("");

					//block8_sra3(data);
					for(int ky2=0;ky2<8;++ky2)
					{
						for(int kx2=0;kx2<8;++kx2)
							dst->data[(src->iw*((ih[kc]>>3)*ky2+(ky>>3))+(iw[kc]>>3)*kx2+(kx>>3))<<2|kc]=data[ky2<<3|kx2];
					}
				}
				else
				{
					for(int ky2=0;ky2<8;++ky2)
					{
						for(int kx2=0;kx2<8;++kx2)
							data[ky2<<3|kx2]=src->data[(src->iw*((ih[kc]>>3)*ky2+(ky>>3))+(iw[kc]>>3)*kx2+(kx>>3))<<2|kc];
					}
					if(matrix)
					{
						for(int k=0;k<64;++k)
							data[k]*=matrix[k];
					}
					block8_dct_ls_inv(data);
					block8_transpose(data);
					block8_dct_ls_inv(data);
					//block8_sra3(data);
					for(int ky2=0;ky2<8;++ky2)
					{
						for(int kx2=0;kx2<8;++kx2)
							dst->data[(src->iw*(ky+ky2)+kx+kx2)<<2|kc]=data[ky2<<3|kx2];
					}
				}

#if 0
				__m128i tmp;
				__m128i row0=_mm_load_si128((const __m128i*)block+0);
				__m128i row1=_mm_load_si128((const __m128i*)block+1);
				__m128i row2=_mm_load_si128((const __m128i*)block+2);
				__m128i row3=_mm_load_si128((const __m128i*)block+3);
				__m128i row4=_mm_load_si128((const __m128i*)block+4);
				__m128i row5=_mm_load_si128((const __m128i*)block+5);
				__m128i row6=_mm_load_si128((const __m128i*)block+6);
				__m128i row7=_mm_load_si128((const __m128i*)block+7);
				
				__m128i rot0_0 = dct_const(stbi__f2f(0.5411961f), stbi__f2f(0.5411961f) + stbi__f2f(-1.847759065f));
				__m128i rot0_1 = dct_const(stbi__f2f(0.5411961f) + stbi__f2f( 0.765366865f), stbi__f2f(0.5411961f));
				__m128i rot1_0 = dct_const(stbi__f2f(1.175875602f) + stbi__f2f(-0.899976223f), stbi__f2f(1.175875602f));
				__m128i rot1_1 = dct_const(stbi__f2f(1.175875602f), stbi__f2f(1.175875602f) + stbi__f2f(-2.562915447f));
				__m128i rot2_0 = dct_const(stbi__f2f(-1.961570560f) + stbi__f2f( 0.298631336f), stbi__f2f(-1.961570560f));
				__m128i rot2_1 = dct_const(stbi__f2f(-1.961570560f), stbi__f2f(-1.961570560f) + stbi__f2f( 3.072711026f));
				__m128i rot3_0 = dct_const(stbi__f2f(-0.390180644f) + stbi__f2f( 2.053119869f), stbi__f2f(-0.390180644f));
				__m128i rot3_1 = dct_const(stbi__f2f(-0.390180644f), stbi__f2f(-0.390180644f) + stbi__f2f( 1.501321110f));
				{
				__m128i rot0_0lo = _mm_unpacklo_epi16((row2), (row6));
				__m128i rot0_0hi = _mm_unpackhi_epi16((row2), (row6));
				__m128i t2e_l = _mm_madd_epi16(rot0_0lo, rot0_0);
				__m128i t2e_h = _mm_madd_epi16(rot0_0hi, rot0_0);
				__m128i t3e_l = _mm_madd_epi16(rot0_0lo, rot0_1);
				__m128i t3e_h = _mm_madd_epi16(rot0_0hi, rot0_1);
				__m128i sum04 = _mm_add_epi16(row0, row4);
				__m128i dif04 = _mm_sub_epi16(row0, row4);
				__m128i t0e_l = _mm_srai_epi32(_mm_unpacklo_epi16(_mm_setzero_si128(), (sum04)), 4);
				__m128i t0e_h = _mm_srai_epi32(_mm_unpackhi_epi16(_mm_setzero_si128(), (sum04)), 4);
				__m128i t1e_l = _mm_srai_epi32(_mm_unpacklo_epi16(_mm_setzero_si128(), (dif04)), 4);
				__m128i t1e_h = _mm_srai_epi32(_mm_unpackhi_epi16(_mm_setzero_si128(), (dif04)), 4);
				__m128i x0_l = _mm_add_epi32(t0e_l, t3e_l); __m128i x0_h = _mm_add_epi32(t0e_h, t3e_h);
				__m128i x3_l = _mm_sub_epi32(t0e_l, t3e_l); __m128i x3_h = _mm_sub_epi32(t0e_h, t3e_h);
				__m128i x1_l = _mm_add_epi32(t1e_l, t2e_l); __m128i x1_h = _mm_add_epi32(t1e_h, t2e_h);
				__m128i x2_l = _mm_sub_epi32(t1e_l, t2e_l); __m128i x2_h = _mm_sub_epi32(t1e_h, t2e_h);
				__m128i rot2_0lo = _mm_unpacklo_epi16((row7), (row3));
				__m128i rot2_0hi = _mm_unpackhi_epi16((row7), (row3));
				__m128i y0o_l = _mm_madd_epi16(rot2_0lo, rot2_0);
				__m128i y0o_h = _mm_madd_epi16(rot2_0hi, rot2_0);
				__m128i y2o_l = _mm_madd_epi16(rot2_0lo, rot2_1);
				__m128i y2o_h = _mm_madd_epi16(rot2_0hi, rot2_1);
				__m128i rot3_0lo = _mm_unpacklo_epi16((row5), (row1));
				__m128i rot3_0hi = _mm_unpackhi_epi16((row5), (row1));
				__m128i y1o_l = _mm_madd_epi16(rot3_0lo, rot3_0);
				__m128i y1o_h = _mm_madd_epi16(rot3_0hi, rot3_0);
				__m128i y3o_l = _mm_madd_epi16(rot3_0lo, rot3_1);
				__m128i y3o_h = _mm_madd_epi16(rot3_0hi, rot3_1);
				__m128i sum17 = _mm_add_epi16(row1, row7);
				__m128i sum35 = _mm_add_epi16(row3, row5);
				__m128i rot1_0lo = _mm_unpacklo_epi16((sum17), (sum35));
				__m128i rot1_0hi = _mm_unpackhi_epi16((sum17), (sum35));
				__m128i y4o_l = _mm_madd_epi16(rot1_0lo, rot1_0);
				__m128i y4o_h = _mm_madd_epi16(rot1_0hi, rot1_0);
				__m128i y5o_l = _mm_madd_epi16(rot1_0lo, rot1_1);
				__m128i y5o_h = _mm_madd_epi16(rot1_0hi, rot1_1);
				__m128i x4_l = _mm_add_epi32(y0o_l, y4o_l);
				__m128i x4_h = _mm_add_epi32(y0o_h, y4o_h);
				__m128i x5_l = _mm_add_epi32(y1o_l, y5o_l);
				__m128i x5_h = _mm_add_epi32(y1o_h, y5o_h);
				__m128i x6_l = _mm_add_epi32(y2o_l, y5o_l);
				__m128i x6_h = _mm_add_epi32(y2o_h, y5o_h);
				__m128i x7_l = _mm_add_epi32(y3o_l, y4o_l);
				__m128i x7_h = _mm_add_epi32(y3o_h, y4o_h);
				{
				__m128i abiased_l = _mm_add_epi32(x0_l, bias_1);
				__m128i abiased_h = _mm_add_epi32(x0_h, bias_1);
				__m128i sum_l = _mm_add_epi32(abiased_l, x7_l);
				__m128i sum_h = _mm_add_epi32(abiased_h, x7_h);
				__m128i dif_l = _mm_sub_epi32(abiased_l, x7_l);
				__m128i dif_h = _mm_sub_epi32(abiased_h, x7_h);
				row0 = _mm_packs_epi32(_mm_srai_epi32(sum_l, 17), _mm_srai_epi32(sum_h, 17));
				row7 = _mm_packs_epi32(_mm_srai_epi32(dif_l, 17), _mm_srai_epi32(dif_h, 17));
				}
				{
				__m128i abiased_l = _mm_add_epi32(x1_l, bias_1);
				__m128i abiased_h = _mm_add_epi32(x1_h, bias_1);
				__m128i sum_l = _mm_add_epi32(abiased_l, x6_l);
				__m128i sum_h = _mm_add_epi32(abiased_h, x6_h);
				__m128i dif_l = _mm_sub_epi32(abiased_l, x6_l);
				__m128i dif_h = _mm_sub_epi32(abiased_h, x6_h);
				row1 = _mm_packs_epi32(_mm_srai_epi32(sum_l, 17), _mm_srai_epi32(sum_h, 17));
				row6 = _mm_packs_epi32(_mm_srai_epi32(dif_l, 17), _mm_srai_epi32(dif_h, 17));
				}
				{
				__m128i abiased_l = _mm_add_epi32(x2_l, bias_1);
				__m128i abiased_h = _mm_add_epi32(x2_h, bias_1);
				__m128i sum_l = _mm_add_epi32(abiased_l, x5_l);
				__m128i sum_h = _mm_add_epi32(abiased_h, x5_h);
				__m128i dif_l = _mm_sub_epi32(abiased_l, x5_l);
				__m128i dif_h = _mm_sub_epi32(abiased_h, x5_h);
				row2 = _mm_packs_epi32(_mm_srai_epi32(sum_l, 17), _mm_srai_epi32(sum_h, 17));
				row5 = _mm_packs_epi32(_mm_srai_epi32(dif_l, 17), _mm_srai_epi32(dif_h, 17));
				}
				{
				__m128i abiased_l = _mm_add_epi32(x3_l, bias_1);
				__m128i abiased_h = _mm_add_epi32(x3_h, bias_1);
				__m128i sum_l = _mm_add_epi32(abiased_l, x4_l);
				__m128i sum_h = _mm_add_epi32(abiased_h, x4_h);
				__m128i dif_l = _mm_sub_epi32(abiased_l, x4_l);
				__m128i dif_h = _mm_sub_epi32(abiased_h, x4_h);
				row3 = _mm_packs_epi32(_mm_srai_epi32(sum_l, 17), _mm_srai_epi32(sum_h, 17));
				row4 = _mm_packs_epi32(_mm_srai_epi32(dif_l, 17), _mm_srai_epi32(dif_h, 17));
				}
				}
#endif
			}
		}
	}
}
#if 0
void jpeg_quantize(Image *src, int quality, int fwd)//quality [1, 100]
{
	int YQT[]=
	{
		16,11,10,16,24,40,51,61,
		12,12,14,19,26,58,60,55,
		14,13,16,24,40,57,69,56,
		14,17,22,29,51,87,80,62,
		18,22,37,56,68,109,103,77,
		24,35,55,64,81,104,113,92,
		49,64,78,87,103,121,120,101,
		72,92,95,98,112,100,103,99,
	};
	int UVQT[]=
	{
		17,18,24,47,99,99,99,99,
		18,21,26,66,99,99,99,99,
		24,26,56,99,99,99,99,99,
		47,66,99,99,99,99,99,99,
		99,99,99,99,99,99,99,99,
		99,99,99,99,99,99,99,99,
		99,99,99,99,99,99,99,99,
		99,99,99,99,99,99,99,99,
	};
	int aasf[]=
	{
		(int)(0x10000*2.828427125),
		(int)(0x10000*2.828427125),
		(int)(0x10000*2.828427125),
		(int)(0x10000*2.828427125),
		(int)(0x10000*2.828427125),
		(int)(0x10000*2.828427125),
		(int)(0x10000*2.828427125),
		(int)(0x10000*2.828427125),

		//(int)(0x10000*2.828427125),
		//(int)(0x10000*2.828427125*1.387039845),
		//(int)(0x10000*2.828427125*1.306562965),
		//(int)(0x10000*2.828427125*1.175875602),
		//(int)(0x10000*2.828427125),
		//(int)(0x10000*2.828427125*0.785694958),
		//(int)(0x10000*2.828427125*0.541196100),
		//(int)(0x10000*2.828427125*0.275899379),
	};
	int subsampled=quality<=90;
	int iw[]=
	{
		src->iw,
		src->iw>>subsampled,
		src->iw>>subsampled,
		src->iw,
	};
	int ih[]=
	{
		src->ih,
		src->ih>>subsampled,
		src->ih>>subsampled,
		src->ih,
	};
	quality=CLAMP(1, quality, 100);
	int s=quality<50?5000/quality:200-2*quality;
	for(int ky=0;ky<8;++ky)
	{
		for(int kx=0;kx<8;++kx)
		{
			int idx=ky<<3|kx, val;
			
			val=(YQT[idx]*s+50)/100;
			val=CLAMP(1, val, 255);
			if(fwd)
				val=(int)((1LL<<56)/((long long)val*aasf[kx]*aasf[ky]));
			else
				val=(int)((long long)val*aasf[kx]*aasf[ky]>>16);
			YQT[idx]=val;

			val=(UVQT[idx]*s+50)/100;
			val=CLAMP(1, val, 255);
			if(fwd)
				val=(int)((1LL<<56)/((long long)val*aasf[kx]*aasf[ky]));
			else
				val=(int)((long long)val*aasf[kx]*aasf[ky]>>16);
			UVQT[idx]=val;
		}
	}
	for(int kc=0;kc<src->nch;++kc)
	{
		int *qt=kc?UVQT:YQT;
		for(int ky=0;ky<ih[kc];ky+=8)
		{
			for(int kx=0;kx<iw[kc];kx+=8)
			{
				if(fwd)
				{
					for(int ky2=0;ky2<8;++ky2)
					{
						for(int kx2=0;kx2<8;++kx2)
						{
							int idx=(src->iw*(ky+ky2)+kx+kx2)<<2|kc;
							int val=src->data[idx];
							val=(long long)val*qt[ky2<<3|kx2]>>24;
							src->data[idx]=val;
						}
					}
				}
				else
				{
					for(int ky2=0;ky2<8;++ky2)
					{
						for(int kx2=0;kx2<8;++kx2)
						{
							int idx=(src->iw*(ky+ky2)+kx+kx2)<<2|kc;
							int val=src->data[idx];
							val=(long long)val*qt[ky2<<3|kx2]>>16;
							src->data[idx]=val;
						}
					}
				}
			}
		}
	}
}
#endif

#define T52_SYM_EXP 5
#define T52_SYM_MSB 2
#define T52_SYM_LSB 0

#define T52_CTX_EXP 2
#define T52_CTX_MSB 1
#define T52_CTX_LSB 0
#define T52_CTX_DIM 64
typedef struct HybridUintStruct
{
	unsigned short token, nbits;
	unsigned bypass;
} HybridUint;
static void hybriduint_encode(int val, int exp, int msb, int lsb, HybridUint *hu)
{
	int token, bypass, nbits;
	val=val<<1^-(val<0);//pack sign
	if(val<(1<<exp))
	{
		token=val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=floor_log2_32((unsigned)val);
		int mantissa=val-(1<<lgv);
		token = (1<<exp) + (
				(lgv-exp)<<(msb+lsb)|
				(mantissa>>(lgv-msb))<<lsb|
				(mantissa&((1<<lsb)-1))
			);
		nbits=lgv-(msb+lsb);
		bypass=val>>lsb&((1LL<<nbits)-1);
	}
	hu->token=token;
	hu->bypass=bypass;
	hu->nbits=nbits;
}

static void update_hist(unsigned short *curr_hist, unsigned *curr_CDF, int cdfsize, int token)
{
	++curr_hist[cdfsize];
	++curr_hist[token];
	if(curr_hist[cdfsize]>cdfsize&&!(curr_hist[cdfsize]&63))//update CDF
	{
		int sum=curr_hist[cdfsize];
		long long c=0;
		for(int ks=0;ks<cdfsize;++ks)
		{
			long long freq=curr_hist[ks];
			curr_CDF[ks]=(int)(c*((1LL<<16)-cdfsize)/sum)+ks;
			c+=freq;
		}
		curr_CDF[cdfsize]=1<<16;
	}
	if(curr_hist[cdfsize]>=0x1800)//rescale hist
	{
		int sum=0;
		for(int k=0;k<cdfsize;++k)
		{
			curr_hist[k]=(curr_hist[k]+1)>>1;
			sum+=curr_hist[k];
		}
		curr_hist[cdfsize]=sum;
	}
}
int t52_encode(Image const *src, ArrayHandle *data, int loud)
{
	double t_start=time_sec();
	int nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);//TODO
	UPDATE_MIN(nch, src->nch);
	if(loud)
	{
		int maxdepth=calc_maxdepth(src, 0);
		acme_strftime(g_buf, G_BUF_SIZE, "%Y-%m-%d_%H-%M-%S");
		printf("T52 JPEG-recompress  Enc %s  CWHD %d*%d*%d*%d/8\n", g_buf, nch, src->iw, src->ih, maxdepth);
	}
	if(nch==3&&(src->depth[0]!=8||src->depth[1]!=8||src->depth[2]!=8)||nch==1&&src->depth[0]!=8||nch!=3&&nch!=1)
	{
		if(loud)
			printf("Not a JPEG\n");
		return 0;
	}
	if(src->iw&15||src->ih&15)
	{
		if(loud)
			printf("Dimensions must be divisible by 16\n");
		return 0;
	}
	Image *im2=0, *im3=0;
	image_copy(&im2, src);

	image_copy(&im3, src);//
	//image_copy_nodata(&im3, src);
	
	int cdfsize=127;
	unsigned short *hist=(unsigned short*)malloc((cdfsize+1LL)*sizeof(short[3*T52_CTX_DIM*T52_CTX_DIM+3]));
	unsigned *CDF=(unsigned*)malloc((cdfsize+1LL)*sizeof(int[3*T52_CTX_DIM*T52_CTX_DIM+3]));
	if(!im2||!im3||!hist||!CDF)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(hist, 0, (cdfsize+1LL)*sizeof(short[3*T52_CTX_DIM*T52_CTX_DIM+3]));
	int sum=0;
	for(int ks=0;ks<cdfsize;++ks)
	{
		int freq=0x10000/(ks+1);
		CDF[ks]=freq;
		sum+=freq;
		//freq=(int)((long long)freq*freq>>16);
	}
	for(int ks=0, c=0;ks<cdfsize;++ks)
	{
		int freq=CDF[ks];
		CDF[ks]=(int)((long long)c*((1LL<<16)-cdfsize)/sum)+ks;
		c+=freq;
	}
	CDF[cdfsize]=1<<16;
	memfill(CDF+cdfsize+1, CDF, (cdfsize+1LL)*sizeof(int[3*T52_CTX_DIM*T52_CTX_DIM+3-1]), (cdfsize+1LL)*sizeof(int));

	int matrix[]=
	{
		16,11,10,16,24,40,51,61,
		12,12,14,19,26,58,60,55,
		14,13,16,24,40,57,69,56,
		14,17,22,29,51,87,80,62,
		18,22,37,56,68,109,103,77,
		24,35,55,64,81,104,113,92,
		49,64,78,87,103,121,120,101,
		72,92,95,98,112,100,103,99,

		17,18,24,47,99,99,99,99,
		18,21,26,66,99,99,99,99,
		24,26,56,99,99,99,99,99,
		47,66,99,99,99,99,99,99,
		99,99,99,99,99,99,99,99,
		99,99,99,99,99,99,99,99,
		99,99,99,99,99,99,99,99,
		99,99,99,99,99,99,99,99,
	};
	int quality=90;
	quality=CLAMP(1, quality, 100);
	int s=quality<50?5000/quality:200-2*quality;
	for(int k=0;k<128;++k)
	{
		int val=matrix[k];
		val=(val*s+50)/100;
		val=CLAMP(1, val, 255);
		matrix[k]=val;
	}
	int subsampled=quality<=90;
	if(nch==3)
	{
		ct_YCbCr_lossy(im2, 1);
		if(subsampled)
			resample_YUV420(im2, 1);
	}
	dct_8x8(im3, im2, subsampled, matrix, matrix+64, 1);

	dct_8x8(im2, im3, subsampled, matrix, matrix+64, 0);
	if(nch==3)
	{
		if(subsampled)
			resample_YUV420(im2, 0);
		ct_YCbCr_lossy(im2, 0);
	}
	int iw[]=
	{
		src->iw,
		src->iw>>subsampled,
		src->iw>>subsampled,
		src->iw,
	};
	int ih[]=
	{
		src->ih,
		src->ih>>subsampled,
		src->ih>>subsampled,
		src->ih,
	};
	double t_transforms=time_sec()-t_start;
	//compare_bufs_32(src->data, im2->data, src->iw, src->ih, src->nch, 4, "DCT", 0, 1);//
	//calc_depthfromdata(im2->data, im2->iw, im2->ih, im2->depth, im2->src_depth);
	//image_snapshot(im2);//
#if 0
	int vmax[4]={0}, vmin[4]={0};
	for(int kc=0;kc<nch;++kc)
	{
		for(int ky=0, idx=0;ky<ih[kc];++ky)
		{
			for(int kx=0;kx<iw[kc];++kx, ++idx)
			{
				int val=im3->data[idx<<2|kc];
				if((unsigned)val==0xBEBEBEBE)
					LOG_ERROR("");
				UPDATE_MIN(vmin[kc], val);
				UPDATE_MAX(vmax[kc], val);
			}
		}
	}
	for(int kc=0;kc<nch;++kc)
		printf("%c  %d ~ %d\n", "YUVA"[kc], vmin[kc], vmax[kc]);
#endif
	HybridUint hu;
	ArithmeticCoder ec;
	DList list;
	dlist_init(&list, 1, 0x10000, 0);
	ac_enc_init(&ec, &list);
	for(int kc=0;kc<nch;++kc)
	{
		for(int ky=0;ky<ih[kc];++ky)
		{
			for(int kx=0;kx<iw[kc];++kx)
			{
				int
					idx=src->iw*ky+kx,
					N=ky?im3->data[(idx-src->iw)<<2|kc]:0,
					W=kx?im3->data[(idx-1)<<2|kc]:0,
					NW=ky&&kx?im3->data[(idx-1)<<2|kc]:0;
				int val=im3->data[idx<<2|kc];
				if(kx<(iw[kc]>>3)&&ky<(ih[kc]>>3))
				{
					int pred=N+W-NW;
					pred=MEDIAN3(N, W, pred);
					val-=pred;
				}
				hybriduint_encode(N, T52_CTX_EXP, T52_CTX_MSB, T52_CTX_LSB, &hu);
				int ctx=CLAMP(0, hu.token, T52_CTX_DIM-1);
				hybriduint_encode(W, T52_CTX_EXP, T52_CTX_MSB, T52_CTX_LSB, &hu);
				ctx=T52_CTX_DIM*ctx+CLAMP(0, hu.token, T52_CTX_DIM-1);

				hybriduint_encode(val, T52_SYM_EXP, T52_SYM_MSB, T52_SYM_LSB, &hu);
				int ctx_idx=(cdfsize+1)*(T52_CTX_DIM*T52_CTX_DIM*kc+ctx);
				unsigned *curr_CDF=CDF+ctx_idx;
				unsigned short *curr_hist=hist+ctx_idx;
				ac_enc(&ec, hu.token, curr_CDF, cdfsize, 0);
				if(hu.nbits)
				{
					int bypass=hu.bypass, nbits=hu.nbits;
					while(nbits>8)
					{
						ac_enc(&ec, bypass>>(nbits-8)&0xFF, 0, 1<<8, 16-8);
						nbits-=8;
					}
					ac_enc(&ec, bypass&((1<<nbits)-1), 0, 1<<nbits, 16-nbits);
				}
				update_hist(curr_hist, curr_CDF, cdfsize, hu.token);
			}
		}
	}
	int csize_lossy=(int)list.nobj;
	double t_enc1=time_sec()-t_start;

	int vmin[4]={0}, vmax[4]={0};

	unsigned short *curr_hist=hist+(cdfsize+1LL)*3*T52_CTX_DIM*T52_CTX_DIM;
	unsigned *curr_CDF=CDF+(cdfsize+1LL)*3*T52_CTX_DIM*T52_CTX_DIM;
	for(int kc=0;kc<nch;++kc)
	{
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				int delta=src->data[idx<<2|kc]-im2->data[idx<<2|kc];
				//im2->data[idx<<2|kc]=src->data[idx<<2|kc]-im2->data[idx<<2|kc];//
				UPDATE_MIN(vmin[kc], delta);//
				UPDATE_MAX(vmax[kc], delta);//
				hybriduint_encode(delta, T52_SYM_EXP, T52_SYM_MSB, T52_SYM_LSB, &hu);
				ac_enc(&ec, hu.token, curr_CDF, cdfsize, 0);
				if(hu.nbits)
				{
					int bypass=hu.bypass, nbits=hu.nbits;
					while(nbits>8)
					{
						ac_enc(&ec, bypass>>(nbits-8)&0xFF, 0, 1<<8, 16-8);
						nbits-=8;
					}
					ac_enc(&ec, bypass&((1<<nbits)-1), 0, 1<<nbits, 16-nbits);
				}
				update_hist(curr_hist, curr_CDF, cdfsize, hu.token);
			}
		}
	}
	//image_snapshot(im2);
	ac_enc_flush(&ec);
	double t_enc2=time_sec()-t_start;
	if(loud)
	{
		t_enc2-=t_enc1;
		t_enc1-=t_transforms;
		double usize=image_getBMPsize(src);
		printf("\n");
		printf("Transforms ");
		timedelta2str(0, 0, t_transforms);
		printf("\n");
		printf("ENC Lossy  ");
		timedelta2str(0, 0, t_enc1);
		printf("\n");
		printf("ENC LS     ");
		timedelta2str(0, 0, t_enc2);
		printf("\n");

		printf("csize %8d  invCR %10.6lf%%  usize %lf\n", (int)list.nobj, 100.*list.nobj/usize, usize);
		printf("lossy %8d  invCR %10.6lf%%\n", csize_lossy, 100.*csize_lossy/usize);
		printf("ls    %8d  invCR %10.6lf%%\n", (int)list.nobj-csize_lossy, 100.*(list.nobj-csize_lossy)/usize);

		for(int kc=0;kc<nch;++kc)
			printf("%c  %d ~ %d\n", "YUVA"[kc], vmin[kc], vmax[kc]);
	}
	free(im2);
	free(im3);
	return 1;
}
int t52_decode(const unsigned char *data, size_t srclen, Image *dst, int loud)
{
	return 0;
}