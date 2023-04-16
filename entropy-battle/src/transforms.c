#include"battle.h"
#include<stdio.h>//for debugging
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#ifdef __AVX__
#include<immintrin.h>
#else
#include<tmmintrin.h>
#endif
static const char file[]=__FILE__;

static void print_sbuf(short *buf, int bw, int bh, int kc, int bytestride)
{
	static int call=0;
	printf("call %d ch %d\n", call, kc);
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
			printf("%d\t", buf[bytestride*(bw*ky+kx)+kc]);
		printf("\n");
	}
	printf("\n");
	++call;
}
static void print_fbuf(float *buf, int bw, int bh, int nplaces)
{
	static int call=0;
	printf("call %d\n", call);
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
			printf("%.*f\t", nplaces, buf[bw*ky+kx]);
		printf("\n");
	}
	printf("\n");
	++call;
}
void compare_bufs_uint8(unsigned char *b1, unsigned char *b0, int iw, int ih, int symbytes, int bytestride, const char *name, int backward)
{
	ptrdiff_t len=(ptrdiff_t)bytestride*iw*ih;
	int inc=bytestride*(1-(backward<<1));
	for(ptrdiff_t k=backward?len-bytestride:0;k>=0&&k<len;k+=inc)
	{
		//{//
		//	ptrdiff_t idx=k>>2, kx=idx%1920, ky=idx/1920;
		//	if(kx==5&&ky==1)
		//		kx=5;
		//}//
		if(memcmp(b1+k, b0+k, symbytes))
		{
			ptrdiff_t idx=k/bytestride, kx=idx%iw, ky=idx/iw;
			printf("%s error XY (%5lld, %5lld) / %5d x %5d  b1 != b0\n", name, kx, ky, iw, ih);
			for(int kc=0;kc<symbytes;++kc)
				printf("C%d  0x%02X != 0x%02X    %d != %d\n", kc, (unsigned)b1[k+kc], (unsigned)b0[k+kc], (unsigned)b1[k+kc], (unsigned)b0[k+kc]);
			//if(backward)
			//	printf("%s error at %d - %d: 0x%02X != 0x%02X\n", name, (int)len-1, (int)(len-1-k), b1[k], b0[k]);
			//else
			//	printf("%s error at %d: 0x%02X != 0x%02X\n", name, (int)k, b1[k], b0[k]);
			return;
		}
	}
	printf("%s:\tSUCCESS\n", name);
}
void compare_bufs_ps(float *b1, float *b0, int iw, int ih, const char *name, int backward)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	int inc=1-(backward<<1);
	for(ptrdiff_t k=backward?res-1:0;k>=0&&k<res;k+=inc)
	{
		if(fabsf(b1[k]-b0[k])>1e-6)
		{
			ptrdiff_t kx=k%iw, ky=k/iw;
			printf("%s error XY (%5lld, %5lld) / %5d x %5d  b1 != b0  %f != %f  error=%f\n", name, kx, ky, iw, ih, b1[k], b0[k], b1[k]-b0[k]);
			return;
		}
	}
	printf("%s:\tSUCCESS\n", name);
}


void addbuf(unsigned char *buf, int iw, int ih, int nch, int bytestride, int ammount)
{
	for(int kp=0, len=iw*ih*bytestride;kp<len;kp+=bytestride)
	{
		for(int kc=0;kc<nch;++kc)
			buf[kp+kc]+=ammount;
	}
}


//YCoCg-R (lifting-based YCoCg)		8 bit <-> ps,  3 channels, pixel stride 4 bytes,  Y in [0, 1], Co/Cg in [-1, 1],  used by AVC/HEVC/VVC
#ifndef CLAMP
#define CLAMP(LO, X, HI)	(X<LO?LO:(X>HI?HI:X))
#endif
void YCoCg_8bit_ps_fwd(const unsigned char *src, ptrdiff_t res, float *bufY, float *bufCo, float *bufCg)
{
	const float gain=1.f/255;
	for(ptrdiff_t k=0;k<res;++k)
	{
		ptrdiff_t idx=k<<2;
		unsigned R=src[idx], G=src[idx|1], B=src[idx|2];
		int Co=R-B,
			temp=B+(Co>>1),
			Cg=G-temp,
			Y=temp+(Cg>>1);
		bufY[k]=(float)Y*gain;
		bufCo[k]=(float)Co*gain;
		bufCg[k]=(float)Cg*gain;
	}
}
void YCoCg_8bit_ps_inv(const float *bufY, const float *bufCo, const float *bufCg, ptrdiff_t res, unsigned char *dst)
{
	const float gain=255.f;
	for(ptrdiff_t k=0;k<res;++k)
	{
		ptrdiff_t idx=k<<2;
		int Y=(int)(gain*bufY[k]), Co=(int)(gain*bufCo[k]), Cg=(int)(gain*bufCg[k]);
		int temp=Y-(Cg>>1), G=Cg+temp, B=temp-(Co>>1), R=B+Co;
		R=CLAMP(0, R, 255);
		G=CLAMP(0, G, 255);
		B=CLAMP(0, B, 255);
		dst[idx]=R, dst[idx|1]=G, dst[idx|3]=B;
	}
}

void colortransform_xgz_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;
		b-=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_xgz_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		b+=g;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_xyz_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;//what is this transform?
		b-=g;
		g-=(r+b)>>1;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;

		//buf[k  ]=r-g+128;//XGZ
		//buf[k|1]=g;
		//buf[k|2]=b-g+128;

		//buf[k  ]=r;//RYZ
		//buf[k|1]=g-r+128;
		//buf[k|2]=b-r+128;

		//buf[k  ]=r-b+128;//XYBdash
		//buf[k|1]=g-b+128;
		//buf[k|2]=b;
	}
}
void colortransform_xyz_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		g+=(r+b)>>1;//what is this transform???
		b+=g;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocg_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=b;
		b+=r>>1;
		g-=b;
		r+=g>>1;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocg_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];
		
		r-=g>>1;
		g+=b;
		b-=r>>1;
		r+=b;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}


//lifting-based 8bit Haar DWT
void dwt1d_haar_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(int k=0;k<nodd;++k)//predict
		even[k]-=odd[k];
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(int k=0;k<nodd;++k)//update
		odd[k]+=even[k]>>1;


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_haar_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	
	for(int k=0;k<nodd;++k)//unupdate
		odd[k]-=even[k]>>1;
	
	for(int k=0;k<nodd;++k)//unpredict
		even[k]+=odd[k];
	if(extraeven)
		even[nodd]+=odd[nodd-1];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_haar_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_haar_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_haar_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_haar_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_haar_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_haar_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit CDF 5/3 DWT
void dwt1d_cdf53_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	even[0]-=odd[0];
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]>>1;

#if 0
	for(int k=0;k<nodd-!extraeven;++k)//predict
		odd[k]+=128-((even[k]+even[k+1])>>1);
	if(!extraeven)
		odd[nodd-1]+=128-even[nodd-1];
	
	even[0]+=(odd[0]-128)>>1;
	for(int k=1;k<nodd-1;++k)//update
		even[k]+=(odd[k]+odd[k+1])>>2;
	if(extraeven)
		even[nodd]+=(odd[nodd-1]-128)>>1;
#endif


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_cdf53_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];

	
	for(int k=0;k<nodd-!extraeven;++k)//un-update
		odd[k]-=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]>>1;
	
	even[0]+=odd[0];
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]+=odd[nodd-1];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf53_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf53_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf53_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_cdf53_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_cdf53_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_cdf53_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit CDF 9/7 DWT
static const int cdf97_coeff[]=
{
	//'factring wavelet transforms into lifting steps' - page 19
	//'a wavelet tour of signal processing - the sparse way' - page 376
	-0x1960C,	//-1.58613434342059f,	//alpha
	-0x00D90,	//-0.0529801185729f,	//beta
	 0x0E206,	// 0.8829110755309f,	//gamma
	 0x07189,	// 0.4435068520439f,	//delta
				// 1.1496043988602f,	//zeta		output gain is 1.89
};
static void dwt1d_u8_predict(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	even[0]-=odd[0]*coeff>>15;
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(odd[k-1]+odd[k])*coeff>>16;
	if(extraeven)
		even[nodd]-=odd[nodd-1]*coeff>>15;
}
static void dwt1d_u8_update(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(even[k]+even[k+1])*coeff>>16;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]*coeff>>15;
}
static void dwt1d_u8_unpredict(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	even[0]+=odd[0]*coeff>>15;
	for(int k=1;k<nodd;++k)//unpredict
		even[k]+=(odd[k-1]+odd[k])*coeff>>16;
	if(extraeven)
		even[nodd]+=odd[nodd-1]*coeff>>15;
}
static void dwt1d_u8_unupdate(char *odd, char *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//unupdate
		odd[k]-=(even[k]+even[k+1])*coeff>>16;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]*coeff>>15;
}
void dwt1d_cdf97_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

	dwt1d_u8_predict(odd, even, nodd, extraeven, cdf97_coeff[3]);
	dwt1d_u8_update (odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_u8_predict(odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_u8_update (odd, even, nodd, extraeven, cdf97_coeff[0]);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_cdf97_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	dwt1d_u8_unupdate (odd, even, nodd, extraeven, cdf97_coeff[0]);
	dwt1d_u8_unpredict(odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_u8_unupdate (odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_u8_unpredict(odd, even, nodd, extraeven, cdf97_coeff[3]);

	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf97_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf97_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf97_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_cdf97_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_cdf97_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_cdf97_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}


//lifting-based 9/7 DWT
ArrayHandle dwt2d_gensizes(int iw, int ih, int wstop, int hstop, int nstages_override)//calculate dimensions of each DWT stage in descending order
{
	ArrayHandle sizes;
	int lw=floor_log2(iw), lh=floor_log2(ih), lmin=MINVAR(lw, lh);
	ARRAY_ALLOC(DWTSize, sizes, 0, 0, lmin, 0);
	if(wstop<3)
		wstop=3;
	if(hstop<3)
		hstop=3;
	int nstages=0;
	DWTSize *p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);
	p->w=iw;
	p->h=ih;
	for(int w2=iw, h2=ih;w2>=wstop&&h2>=hstop&&(!nstages_override||nstages<nstages_override);++nstages)
	{
		p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);

		w2>>=1;//w=floor(w/2)
		h2>>=1;//h=floor(h/2)
		//w2-=w2>>1;//w=ceil(w/2)
		//h2-=h2>>1;//h=ceil(h/2)

		p->w=w2;
		p->h=h2;
	}
	return sizes;
}
static const float dwt97_coeff[]=
{
	//'factring wavelet transforms into lifting steps' - page 19
	//'a wavelet tour of signal processing - the sparse way' - page 376
	-1.58613434342059f,	//alpha
	-0.0529801185729f,	//beta
	 0.8829110755309f,	//gamma
	 0.4435068520439f,	//delta
	 1.1496043988602f,	//zeta		output gain is 1.89
};
void        dwt1d_ps_scale(float *odd, float *even, int nodd, int neven, float ce, float co)
{
	int k, c2;
	__m128 factor;

	factor=_mm_set1_ps(ce);
	for(k=0, c2=nodd-3;k<c2;k+=4)//x[2k] *= ce
	{
		__m128 val=_mm_load_ps(odd+k);
		val=_mm_mul_ps(val, factor);
		_mm_store_ps(odd+k, val);
	}
	for(;k<nodd;++k)
		odd[k]*=ce;

	factor=_mm_set1_ps(co);
	for(k=0, c2=neven-3;k<c2;k+=4)//x[2k+1] *= co
	{
		__m128 val=_mm_load_ps(even+k);
		val=_mm_mul_ps(val, factor);
		_mm_store_ps(even+k, val);
	}
	for(;k<neven;++k)
		even[k]*=co;
}
void        dwt1d_ps_predict_next(float *odd, float *even, int nodd, int neven, float z10)
{
	int k, c2;
	__m128 f0;
	__m128i sh1, sh2;
	int extraeven=nodd<neven;

	//  _mm_set_epi8(15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0);
	sh1=_mm_set_epi8(-1, -1, -1, -1, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4);
	sh2=_mm_set_epi8( 3,  2,  1,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
	f0=_mm_set1_ps(z10);

	//c2=neven-4;//same
	//c2=nodd-3-!extraeven;

	for(k=0, c2=neven-4;k<c2;k+=4)//x[2k+1] += z10*(x[2k]+x[2k+2])		odd[k]+=z10*(even[k]+even[k+1])
	{
		__m128 vo=_mm_load_ps(odd+k);
		__m128 ve0=_mm_load_ps(even+k);
		
		__m128 ve1=_mm_load_ss(even+k+4);
		__m128 temp=_mm_castsi128_ps(_mm_shuffle_epi8(_mm_castps_si128(ve0), sh1));
		ve1=_mm_castsi128_ps(_mm_shuffle_epi8(_mm_castps_si128(ve1), sh2));
		ve1=_mm_or_ps(ve1, temp);
		
		//__m128 ve1=_mm_loadu_ps(even+k+1);//unaligned

		ve0=_mm_add_ps(ve0, ve1);
		ve0=_mm_mul_ps(ve0, f0);
		ve0=_mm_add_ps(ve0, vo);
		_mm_store_ps(odd+k, ve0);
	}
	for(;k<nodd;++k)
		odd[k]+=z10*(even[k]+even[k+(k+1<neven)]);

	//for(c2=nodd-1-!extraeven;k<c2;++k)
	//	odd[k]+=z10*(even[k]+even[k+1]);
	//--nodd;
	//odd[nodd]+=z10*(even[nodd]+even[nodd+extraeven]);//symmetric pad if no extra even

	//if(nodd==neven)//count was even: use symmetric padding
	//	odd[nodd-1]+=z10*(even[nodd-1]+even[nodd-1]);
	//else//count was odd
	//	odd[nodd-1]+=z10*(even[nodd-1]+even[nodd]);
}
void        dwt1d_ps_update_prev(float *odd, float *even, int nodd, int neven, float z0m1)
{
	int k, c2;
	__m128 f0;
	__m128i shift=_mm_set_epi8(11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0, -1, -1, -1, -1);

	//even[0]+=z0m1*(odd[0]+odd[0]);//always pad at start
	//for(k=1;k<3;++k)
	//	even[k]+=z0m1*(odd[k-1]+odd[k]);
	
	f0=_mm_set1_ps(z0m1);
	for(k=0, c2=nodd-3;k<c2;k+=4)//x[2k] += z0m1*(x[2k-1]+x[2k+1])		even[k]+=z0m1*(odd[k-1]+odd[k])
	{
		__m128 ve=_mm_load_ps(even+k);
		__m128 vo0=_mm_load_ps(odd+k);

		__m128 vo1=_mm_load_ss(odd+k-(k>0));//always pad symmetric at start
		__m128 temp=_mm_castsi128_ps(_mm_shuffle_epi8(_mm_castps_si128(vo0), shift));
		vo1=_mm_or_ps(vo1, temp);

		//__m128 vo1=_mm_loadu_ps(odd+k-1);//unaligned

		vo0=_mm_add_ps(vo0, vo1);
		vo0=_mm_mul_ps(vo0, f0);
		vo0=_mm_add_ps(vo0, ve);
		_mm_storeu_ps(even+k, vo0);
	}
	for(;k<neven;++k)
		even[k]+=z0m1*(odd[k-1]+odd[k-(k>=nodd)]);

#if 0
	__m128 f0=_mm_set1_ps(z0m1);
	int k=halfsize-5;
	for(;k>0;k-=4)//x[2k] += z0m1*(x[2k-1]+x[2k+1])
	{
		__m128 ve=_mm_loadu_ps(even+k);
		__m128 vo0=_mm_loadu_ps(odd+k);
		__m128 vo1=_mm_loadu_ps(odd+k-1);
		vo0=_mm_add_ps(vo0, vo1);
		vo0=_mm_mul_ps(vo0, f0);
		vo0=_mm_add_ps(vo0, ve);
		_mm_storeu_ps(even+k, vo0);
	}
	k+=4;
	for(;k>0;--k)
		even[k]+=z0m1*(odd[k]+odd[k-1]);
	even[0]+=z0m1*(odd[0]+odd[1]);
#endif
}
void        dwt1d_ps_fwd(float *buffer, int count, int stride, float *b2)
{
	int nodd=count>>1, neven=count-nodd;//counting from 0:  nOdd <= nEven			nOdd = nEven - (count&1)
	float *odd=b2, *even=odd+nodd;
	even=(float*)(((size_t)even+15)&~15);//16 byte alignment
	
	int k, ks;
	for(k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(count&1)//extra even if count is odd
		even[k]=buffer[ks];
	
	dwt1d_ps_scale       (odd, even, nodd, neven, dwt97_coeff[4], 1/dwt97_coeff[4]);
	dwt1d_ps_predict_next(odd, even, nodd, neven, dwt97_coeff[3]);
	dwt1d_ps_update_prev (odd, even, nodd, neven, dwt97_coeff[2]);
	dwt1d_ps_predict_next(odd, even, nodd, neven, dwt97_coeff[1]);
	dwt1d_ps_update_prev (odd, even, nodd, neven, dwt97_coeff[0]);

	for(k=0, ks=0;k<nodd;++k, ks+=stride)
		buffer[ks]=odd[k];
	for(k=0;k<neven;++k, ks+=stride)//skip align pad
		buffer[ks]=even[k];
}
void        dwt1d_ps_inv(float *buffer, int count, int stride, float *b2)
{
	int nodd=count>>1, neven=count-nodd;
	float *odd=b2, *even=b2+nodd;
	even=(float*)(((size_t)even+15)&~15);//16 byte alignment

	int k, ks;
	for(k=0, ks=0;k<nodd;++k, ks+=stride)
		odd[k]=buffer[ks];
	for(k=0;k<neven;++k, ks+=stride)//skip align pad
		even[k]=buffer[ks];
	
	dwt1d_ps_update_prev (odd, even, nodd, neven, -dwt97_coeff[0]);
	dwt1d_ps_predict_next(odd, even, nodd, neven, -dwt97_coeff[1]);
	dwt1d_ps_update_prev (odd, even, nodd, neven, -dwt97_coeff[2]);
	dwt1d_ps_predict_next(odd, even, nodd, neven, -dwt97_coeff[3]);
	dwt1d_ps_scale       (odd, even, nodd, neven, 1/dwt97_coeff[4], dwt97_coeff[4]);

	for(k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(count&1)//extra even if count is odd
		buffer[ks]=even[k];
}
void        dwt2d_ps_fwd(float *buffer, DWTSize *sizes, int nsizes)//sizes in descending order
{
	if(!nsizes)
		return;
	int iw=sizes->w, ih=sizes->h,
		tsize=MAXVAR(iw, ih)+16;//16 is for align pad
	float *temp=(float*)_mm_malloc(tsize*sizeof(float), 16);
	for(int it=0;it<nsizes-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_ps_fwd(buffer+iw*ky, w2, 1, temp);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_ps_fwd(buffer+kx, h2, iw, temp);
	}
	_mm_free(temp);
}
void        dwt2d_ps_inv(float *buffer, DWTSize *sizes, int nsizes)//sizes in descending order
{
	if(!nsizes)
		return;
	int iw=sizes->w, ih=sizes->h,
		tsize=MAXVAR(iw, ih)+16;
	float *temp=(float*)_mm_malloc(tsize*sizeof(float), 16);
	for(int it=nsizes-2;it>=0;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_ps_inv(buffer+kx, h2, iw, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_ps_inv(buffer+iw*ky, w2, 1, temp);
	}
	_mm_free(temp);
}


//Haar (squeeze) DWT: n bit -> {n, n+1, n+1, n+2} bit
void squeeze_1d_fwd(short *buffer, int count, int stride, int vmax, short *b2)//nOdd <= nEven			nOdd = nEven - (count&1)
{
	int floorhalf=count>>1;
	short *odd=b2, *even=b2+floorhalf;
	
	for(int k=0, ks=0;k<floorhalf;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
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
			diff=even[k]-odd[k]+vmax;

		//int neg=diff<0;
		//diff^=-neg;
		//diff+=neg;
		//diff<<=1;
		//diff|=neg;

		//if(diff<0)
		//	LOG_ERROR("Wrong vmax %d", vmax);

		odd[k]=av;
		even[k]=diff;
	}
	//dwt2_1d_9_7_fwd(even, odd, halfsize);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void squeeze_2d_fwd(short *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int vmax, short *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			squeeze_1d_fwd(buffer+rowlen*ky, w2, stride, vmax, temp);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			squeeze_1d_fwd(buffer+stride*kx, h2, rowlen, vmax<<(kx>=(w2>>1)), temp);
	}
}
void squeeze_1d_inv(short *buffer, int count, int stride, int vmax, short *b2)
{
	int floorhalf=count>>1;//ceil(count/2)
	short *odd=b2, *even=b2+floorhalf;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];

	for(int k=0;k<floorhalf;++k)
	{
		short av=odd[k], diff=even[k]-vmax;
		
		//int neg=diff&1;
		//diff>>=1;
		//diff^=-neg;
		//diff+=neg;
		//diff|=neg<<7&-!diff;

		even[k]=((av<<1)+diff+(diff>0?-(diff&1):(diff&1)))>>1;//JPEG XL
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
void squeeze_2d_inv(short *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int vmax, short *temp)//temp size is maxdim*sizeof(short)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, tsize=MAXVAR(iw, ih), rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			squeeze_1d_inv(buffer+stride*kx, h2, rowlen, vmax, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			squeeze_1d_inv(buffer+rowlen*ky, w2, stride, vmax, temp);
	}
}


//Haar DWT: 8 bit -> 9 bit
void haar_1d_fwd(short *buffer, int count, int stride, short *b2)//nOdd <= nEven			nOdd = nEven - (count&1)
{
	int floorhalf=count>>1,
		ceilhalf=count-floorhalf;//ceil(count/2)
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
	//dwt2_1d_9_7_fwd(even, odd, halfsize);

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
		//print_sbuf(*ret, bw, bh, kc, bytestride);//

		for(int w2=bw, h2=bh, it=0;w2>=3&&h2>=3&&(!nstages||it<nstages);++it)
		{
			for(int ky=0;ky<h2;++ky)//horizontal DWT
				haar_1d_fwd(*ret+rowlen*ky+kc, w2, bytestride, temp);
			
			//print_sbuf(*ret, bw, bh, kc, bytestride);//

			for(int kx=0;kx<w2;++kx)//vertical DWT
				haar_1d_fwd(*ret+bytestride*kx+kc, h2, rowlen, temp);
			
			//print_sbuf(*ret, bw, bh, kc, bytestride);//

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
		//print_sbuf(buf, bw, bh, kc, bytestride);//

		for(int it=nsizes-1;it>=0;--it)
		{
			int w2=sizes[it<<1], h2=sizes[it<<1|1];

			for(int kx=0;kx<w2;++kx)//vertical IDWT
				haar_1d_inv(buf+bytestride*kx+kc, h2, rowlen, temp);

			//print_sbuf(buf, bw, bh, kc, bytestride);//

			for(int ky=0;ky<h2;++ky)//horizontal IDWT
				haar_1d_inv(buf+rowlen*ky+kc, w2, bytestride, temp);

			//print_sbuf(buf, bw, bh, kc, bytestride);//
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


//integer DCT-II		TODO benchmark vs FCT, where is the inverse?
void transpose_block8x8_ps(float *ptr)//32 byte aligned
{
	//AVX implementation is broken
#if 0
//#ifdef __AVX__
	//https://stackoverflow.com/questions/25622745/transpose-an-8x8-float-using-avx-avx2
	__m256
		r0=_mm256_load_ps(ptr   ),
		r4=_mm256_load_ps(ptr+ 8),
		r6=_mm256_load_ps(ptr+16),
		r2=_mm256_load_ps(ptr+24),
		r7=_mm256_load_ps(ptr+32),
		r3=_mm256_load_ps(ptr+40),
		r5=_mm256_load_ps(ptr+48),
		r1=_mm256_load_ps(ptr+56);

	__m256 t0, t1, t2, t3, t4, t5, t6, t7;
	__m256 tt0, tt1, tt2, tt3, tt4, tt5, tt6, tt7;
	t0=_mm256_unpacklo_ps(r0, r1);
	t1=_mm256_unpackhi_ps(r0, r1);
	t2=_mm256_unpacklo_ps(r2, r3);
	t3=_mm256_unpackhi_ps(r2, r3);
	t4=_mm256_unpacklo_ps(r4, r5);
	t5=_mm256_unpackhi_ps(r4, r5);
	t6=_mm256_unpacklo_ps(r6, r7);
	t7=_mm256_unpackhi_ps(r6, r7);
	tt0=_mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(1,0,1,0));
	tt1=_mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(3,2,3,2));
	tt2=_mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1,0,1,0));
	tt3=_mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(3,2,3,2));
	tt4=_mm256_shuffle_ps(t4, t6, _MM_SHUFFLE(1,0,1,0));
	tt5=_mm256_shuffle_ps(t4, t6, _MM_SHUFFLE(3,2,3,2));
	tt6=_mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(1,0,1,0));
	tt7=_mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(3,2,3,2));
	r0=_mm256_permute2f128_ps(tt0, tt4, 0x20);
	r1=_mm256_permute2f128_ps(tt1, tt5, 0x20);
	r2=_mm256_permute2f128_ps(tt2, tt6, 0x20);
	r3=_mm256_permute2f128_ps(tt3, tt7, 0x20);
	r4=_mm256_permute2f128_ps(tt0, tt4, 0x31);
	r5=_mm256_permute2f128_ps(tt1, tt5, 0x31);
	r6=_mm256_permute2f128_ps(tt2, tt6, 0x31);
	r7=_mm256_permute2f128_ps(tt3, tt7, 0x31);

	_mm256_store_ps((float*)ptr   , r0);
	_mm256_store_ps((float*)ptr+ 8, r1);
	_mm256_store_ps((float*)ptr+16, r2);
	_mm256_store_ps((float*)ptr+24, r3);
	_mm256_store_ps((float*)ptr+32, r4);
	_mm256_store_ps((float*)ptr+40, r5);
	_mm256_store_ps((float*)ptr+48, r6);
	_mm256_store_ps((float*)ptr+56, r7);
#else
	__m128
		a0lo=_mm_load_ps(ptr   ), a0hi=_mm_load_ps(ptr+ 4),
		a1lo=_mm_load_ps(ptr+ 8), a1hi=_mm_load_ps(ptr+12),
		a2lo=_mm_load_ps(ptr+16), a2hi=_mm_load_ps(ptr+20),
		a3lo=_mm_load_ps(ptr+24), a3hi=_mm_load_ps(ptr+28),
		a4lo=_mm_load_ps(ptr+32), a4hi=_mm_load_ps(ptr+36),
		a5lo=_mm_load_ps(ptr+40), a5hi=_mm_load_ps(ptr+44),
		a6lo=_mm_load_ps(ptr+48), a6hi=_mm_load_ps(ptr+52),
		a7lo=_mm_load_ps(ptr+56), a7hi=_mm_load_ps(ptr+60);

	_MM_TRANSPOSE4_PS(a0lo, a1lo, a2lo, a3lo);
	_MM_TRANSPOSE4_PS(a0hi, a1hi, a2hi, a3hi);
	_MM_TRANSPOSE4_PS(a4lo, a5lo, a6lo, a7lo);
	_MM_TRANSPOSE4_PS(a4hi, a5hi, a6hi, a7hi);
	
	_mm_store_ps(ptr   , a0lo); _mm_store_ps(ptr+ 4, a4lo);
	_mm_store_ps(ptr+ 8, a1lo); _mm_store_ps(ptr+12, a5lo);
	_mm_store_ps(ptr+16, a2lo); _mm_store_ps(ptr+20, a6lo);
	_mm_store_ps(ptr+24, a3lo); _mm_store_ps(ptr+28, a7lo);
	_mm_store_ps(ptr+32, a0hi); _mm_store_ps(ptr+36, a4hi);
	_mm_store_ps(ptr+40, a1hi); _mm_store_ps(ptr+44, a5hi);
	_mm_store_ps(ptr+48, a2hi); _mm_store_ps(ptr+52, a6hi);
	_mm_store_ps(ptr+56, a3hi); _mm_store_ps(ptr+60, a7hi);
#endif
}
void intDCT2_8x8_p32_rows(int *buf)//buf is aligned packed 32 bit as fixed precision
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
void intDCT2_8x8_p32(const unsigned char *buf, int bw, int bh, int nch, int bytestride, unsigned short **ret)
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
					intDCT2_8x8_p32_rows(temp);
					transpose_block8x8_ps((float*)temp);
					intDCT2_8x8_p32_rows(temp);
					transpose_block8x8_ps((float*)temp);
					for(int kx2=0;kx2<8;++kx2)
						ret[0][bytestride*(bw*(ky|ky2)+(kx|kx2))+kc]=temp[ky2<<3|kx2];
				}
			}
		}
	}
}
#if 0
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
	intDCT2_8x8_p32_rows(block);
	memcpy(b2, block, 64*sizeof(int));
	transpose_block8x8_ps((float*)block);

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


//	#define DCT_8x8_PS_SCALED

//fast 8 point DCT-II/III, call transpose then call again, ptr to 64 packed floats is aligned by 32 bytes
void DCT2_8x8_ps(float *ptr)
{
	//Loeffler's factorization using lifting		http://www.sfu.ca/~jiel/papers/c003-bindct-ieee.pdf
#ifdef __AVX__
	__m256
		half=_mm256_set1_ps(0.5f),
		p1=_mm256_set1_ps(13.f/32),//binDCT-L5
		u1=_mm256_set1_ps(11.f/32),
		p2=_mm256_set1_ps(19.f/64),
		u2=_mm256_set1_ps( 9.f/16),
		p3=_mm256_set1_ps(19.f/64),
		p4=_mm256_set1_ps( 3.f/32),
		u3=_mm256_set1_ps( 3.f/16),
		p5=_mm256_set1_ps( 3.f/32);
	
	__m256
		x0=_mm256_load_ps(ptr   ),
		x1=_mm256_load_ps(ptr+ 8),
		x2=_mm256_load_ps(ptr+16),
		x3=_mm256_load_ps(ptr+24),
		x4=_mm256_load_ps(ptr+32),
		x5=_mm256_load_ps(ptr+40),
		x6=_mm256_load_ps(ptr+48),
		x7=_mm256_load_ps(ptr+56);
	
	//new stage after each extra newline
	x7=_mm256_sub_ps(x0, x7);
	x6=_mm256_sub_ps(x1, x6);
	x5=_mm256_sub_ps(x2, x5);
	x4=_mm256_sub_ps(x3, x4);

	x0=_mm256_fnmadd_ps(x7, half, x0);
	x1=_mm256_fnmadd_ps(x6, half, x1);
	x2=_mm256_fnmadd_ps(x5, half, x2);
	x3=_mm256_fnmadd_ps(x4, half, x3);

	x3=_mm256_sub_ps(x0, x3);
	x2=_mm256_sub_ps(x1, x2);
	x7=_mm256_fnmadd_ps(x4, p2, x7);
	x6=_mm256_fnmadd_ps(x5, p4, x6);
	
	x0=_mm256_fnmadd_ps(x3, half, x0);
	x1=_mm256_fnmadd_ps(x2, half, x1);
	x4=_mm256_fmadd_ps(x7, u2, x4);
	x5=_mm256_fmadd_ps(x6, u3, x5);

	x0=_mm256_add_ps(x0, x1);
	x2=_mm256_fmsub_ps(x3, p1, x2);
	x7=_mm256_fnmadd_ps(x4, p3, x7);
	x6=_mm256_fnmadd_ps(x5, p5, x6);

	x1=_mm256_fmsub_ps(x0, half, x1);
	x3=_mm256_fnmadd_ps(x2, u1, x3);
	x6=_mm256_sub_ps(x4, x6);
	x5=_mm256_sub_ps(x7, x5);

	x4=_mm256_fnmadd_ps(x6, half, x4);
	x7=_mm256_fnmadd_ps(x5, half, x7);

	x7=_mm256_add_ps(x7, x4);
	
	x4=_mm256_fmsub_ps(x7, half, x4);
	
	_mm256_store_ps(ptr   , x0);
	_mm256_store_ps(ptr+ 8, x4);
	_mm256_store_ps(ptr+16, x6);
	_mm256_store_ps(ptr+24, x2);
	_mm256_store_ps(ptr+32, x7);
	_mm256_store_ps(ptr+40, x3);
	_mm256_store_ps(ptr+48, x5);
	_mm256_store_ps(ptr+56, x1);
#else
#error TODO
#endif

#if 0
	//Chen's factorization using lifting		http://www.sfu.ca/~jiel/papers/c003-bindct-ieee.pdf
#ifdef __AVX__
	__m256
		half=_mm256_set1_ps(0.5f),
		p1=_mm256_set1_ps(13.f/32),//binDCT-C7
		u1=_mm256_set1_ps(11.f/32),
		p2=_mm256_set1_ps(11.f/16),
		u2=_mm256_set1_ps(15.f/32),
		p3=_mm256_set1_ps( 3.f/16),
		u3=_mm256_set1_ps( 3.f/16),
		p4=_mm256_set1_ps(13.f/32),
		u4=_mm256_set1_ps(11.f/16),
		p5=_mm256_set1_ps(13.f/32);
	
	__m256
		x0=_mm256_load_ps(ptr   ),
		x1=_mm256_load_ps(ptr+ 8),
		x2=_mm256_load_ps(ptr+16),
		x3=_mm256_load_ps(ptr+24),
		x4=_mm256_load_ps(ptr+32),
		x5=_mm256_load_ps(ptr+40),
		x6=_mm256_load_ps(ptr+48),
		x7=_mm256_load_ps(ptr+56);

	//new stage after each extra newline
	x7=_mm256_sub_ps(x0, x7);
	x6=_mm256_sub_ps(x1, x6);
	x5=_mm256_sub_ps(x2, x5);
	x4=_mm256_sub_ps(x3, x4);

	x0=_mm256_sub_ps(x0, _mm256_mul_ps(x7, half));
	x1=_mm256_sub_ps(x1, _mm256_mul_ps(x6, half));
	x2=_mm256_sub_ps(x2, _mm256_mul_ps(x5, half));
	x3=_mm256_sub_ps(x3, _mm256_mul_ps(x4, half));

	x3=_mm256_sub_ps(x0, x3);
	x2=_mm256_sub_ps(x1, x2);
	x5=_mm256_sub_ps(x5, _mm256_mul_ps(x6, p4));

	x0=_mm256_sub_ps(x0, _mm256_mul_ps(x3, half));
	x1=_mm256_sub_ps(x1, _mm256_mul_ps(x2, half));
	x6=_mm256_add_ps(x6, _mm256_mul_ps(x5, u4));

	x5=_mm256_sub_ps(_mm256_mul_ps(x6, p5), x5);

	x5=_mm256_sub_ps(x4, x5);
	x6=_mm256_sub_ps(x7, x6);
	x1=_mm256_sub_ps(x0, x1);
	x2=_mm256_sub_ps(_mm256_mul_ps(x3, p1), x2);

	x4=_mm256_sub_ps(x4, _mm256_mul_ps(x5, half));
	x7=_mm256_sub_ps(x7, _mm256_mul_ps(x6, half));

	x4=_mm256_sub_ps(_mm256_mul_ps(x7, p3), x4);
	x5=_mm256_add_ps(_mm256_mul_ps(x6, p2), x5);
	
	x0=_mm256_sub_ps(x0, _mm256_mul_ps(x1, half));
	x2=_mm256_sub_ps(x2, _mm256_mul_ps(x1, half));

#else
	__m128
		half=_mm_set1_ps(0.5f),
		p1=_mm_set1_ps(13.f/32),//binDCT-C7
		u1=_mm_set1_ps(11.f/32),
		p2=_mm_set1_ps(11.f/16),
		u2=_mm_set1_ps(15.f/32),
		p3=_mm_set1_ps( 3.f/16),
		u3=_mm_set1_ps( 3.f/16),
		p4=_mm_set1_ps(13.f/32),
		u4=_mm_set1_ps(11.f/16),
		p5=_mm_set1_ps(13.f/32);
	
	__m128
		x0lo=_mm_load_ps(ptr   ), x0hi=_mm_load_ps(ptr+ 4),
		x1lo=_mm_load_ps(ptr+ 8), x1hi=_mm_load_ps(ptr+12),
		x2lo=_mm_load_ps(ptr+16), x2hi=_mm_load_ps(ptr+20),
		x3lo=_mm_load_ps(ptr+24), x3hi=_mm_load_ps(ptr+28),
		x4lo=_mm_load_ps(ptr+32), x4hi=_mm_load_ps(ptr+36),
		x5lo=_mm_load_ps(ptr+40), x5hi=_mm_load_ps(ptr+44),
		x6lo=_mm_load_ps(ptr+48), x6hi=_mm_load_ps(ptr+52),
		x7lo=_mm_load_ps(ptr+56), x7hi=_mm_load_ps(ptr+60);

	x7lo=_mm_sub_ps(x0lo, x7lo), x7hi=_mm_sub_ps(x0hi, x7hi);
#endif
#endif

#if 0
	//Loeffler's factorization using lifting		http://www.sfu.ca/~jiel/papers/c003-bindct-ieee.pdf
	__m128
		half=_mm_set1_ps(0.5f),
		p1=_mm_set1_ps(13.f/32),//binDCT-L5
		u1=_mm_set1_ps(11.f/32),
		p2=_mm_set1_ps(19.f/64),
		u2=_mm_set1_ps( 9.f/16),
		p3=_mm_set1_ps(19.f/64),
		p4=_mm_set1_ps( 3.f/32),
		u3=_mm_set1_ps( 3.f/16),
		p5=_mm_set1_ps( 3.f/32);
	
#ifdef DCT_8x8_PS_SCALED
	__m128
		g0=_mm_set1_ps(0.3535533905932737622004221810524f),//sin(pi/4)/2
		g1=_mm_set1_ps(0.7071067811865475244008443621048f),//sin(pi/4)
		g2=_mm_set1_ps(0.4619397662556433780640915946984f),//sin(3pi/8)/2
		g3=_mm_set1_ps(0.5411961001461969843997232053664f),//1/(2sin(3pi/8))
		g4=_mm_set1_ps(0.7071067811865475244008443621048f),//1/sqrt2
		g5=half,//1/2
		g6=half,//1/2
		g7=_mm_set1_ps(0.3535533905932737622004221810524f);//1/sqrt8
#endif

	__m128
		a0lo=_mm_load_ps(ptr   ), a0hi=_mm_load_ps(ptr+ 4),
		a1lo=_mm_load_ps(ptr+ 8), a1hi=_mm_load_ps(ptr+12),
		a2lo=_mm_load_ps(ptr+16), a2hi=_mm_load_ps(ptr+20),
		a3lo=_mm_load_ps(ptr+24), a3hi=_mm_load_ps(ptr+28),
		a4lo=_mm_load_ps(ptr+32), a4hi=_mm_load_ps(ptr+36),
		a5lo=_mm_load_ps(ptr+40), a5hi=_mm_load_ps(ptr+44),
		a6lo=_mm_load_ps(ptr+48), a6hi=_mm_load_ps(ptr+52),
		a7lo=_mm_load_ps(ptr+56), a7hi=_mm_load_ps(ptr+60);

	__m128
		b0lo=_mm_add_ps(a0lo, a7lo), b0hi=_mm_add_ps(a0hi, a7hi),
		b1lo=_mm_add_ps(a1lo, a6lo), b1hi=_mm_add_ps(a1hi, a6hi),
		b2lo=_mm_add_ps(a2lo, a5lo), b2hi=_mm_add_ps(a2hi, a5hi),
		b3lo=_mm_add_ps(a3lo, a4lo), b3hi=_mm_add_ps(a3hi, a4hi),
		b4lo=_mm_sub_ps(a0lo, a7lo), b4hi=_mm_sub_ps(a0hi, a7hi),
		b5lo=_mm_sub_ps(a1lo, a6lo), b5hi=_mm_sub_ps(a1hi, a6hi),
		b6lo=_mm_sub_ps(a2lo, a5lo), b6hi=_mm_sub_ps(a2hi, a5hi),
		b7lo=_mm_sub_ps(a3lo, a4lo), b7hi=_mm_sub_ps(a3hi, a4hi);

	__m128 t0lo, t0hi;

	t0lo=_mm_mul_ps(b4lo, p2  ), t0hi=_mm_mul_ps(b4hi, p2  );//pred2:   x[7]-=p2*x[4]
	b7lo=_mm_sub_ps(b7lo, t0lo), b7hi=_mm_sub_ps(b7hi, t0hi);

	b7lo=_mm_mul_ps(b7lo, u2  ), b7hi=_mm_mul_ps(b7hi, u2  );//update2: x[4]+=u2*x[7]
	b4lo=_mm_add_ps(b4lo, b7lo), b4hi=_mm_add_ps(b4hi, b7hi);

	t0lo=_mm_mul_ps(b4lo, p3  ), t0hi=_mm_mul_ps(b4hi, p3  );//pred3:   x[7]-=p3*x[4]
	b7lo=_mm_sub_ps(b7lo, t0lo), b7hi=_mm_sub_ps(b7hi, t0hi);


	t0lo=_mm_mul_ps(b5lo, p4  ), t0hi=_mm_mul_ps(b5hi, p4  );//pred4:   x[6]-=p4*x[5]
	b6lo=_mm_sub_ps(b6lo, t0lo), b6hi=_mm_sub_ps(b6hi, t0hi);

	t0lo=_mm_mul_ps(b6lo, u3  ), t0hi=_mm_mul_ps(b6hi, u3  );//update3: x[5]+=u3*x[6]
	b5lo=_mm_add_ps(b5lo, t0lo), b5hi=_mm_add_ps(b5hi, t0hi);

	t0lo=_mm_mul_ps(b5lo, p5  ), t0hi=_mm_mul_ps(b5hi, p5  );//pred5:   x[6]-=p5*x[5]
	b6lo=_mm_sub_ps(b6lo, t0lo), b6hi=_mm_sub_ps(b6hi, t0hi);

	a0lo=_mm_add_ps(b0lo, b3lo), a0hi=_mm_add_ps(b0hi, b3hi);
	a1lo=_mm_add_ps(b1lo, b2lo), a1hi=_mm_add_ps(b1hi, b2hi);
	a2lo=_mm_sub_ps(b0lo, b3lo), a2hi=_mm_sub_ps(b0hi, b3hi);
	a3lo=_mm_sub_ps(b1lo, b2lo), a3hi=_mm_sub_ps(b1hi, b2hi);
	a4lo=_mm_add_ps(b4lo, b6lo), a4hi=_mm_add_ps(b4hi, b6hi);
	a5lo=_mm_add_ps(b5lo, b7lo), a5hi=_mm_add_ps(b5hi, b7hi);
	a6lo=_mm_sub_ps(b4lo, b6lo), a6hi=_mm_sub_ps(b4hi, b6hi);
	a7lo=_mm_sub_ps(b5lo, b7lo), a7hi=_mm_sub_ps(b5hi, b7hi);

	a0lo=_mm_add_ps(a0lo, a1lo), a0hi=_mm_add_ps(a0hi, a1hi);//x[0]+=x[1]
	a7lo=_mm_add_ps(a7lo, a4lo), a7hi=_mm_add_ps(a7hi, a4hi);//x[7]+=x[4]

	t0lo=_mm_mul_ps(a0lo, half), t0hi=_mm_mul_ps(a0hi, half);//x[1] = x[0]/2 - x[1]
	a1lo=_mm_sub_ps(t0lo, a1lo), a1hi=_mm_sub_ps(t0hi, a1hi);

	t0lo=_mm_mul_ps(a3lo, p1  ), t0hi=_mm_mul_ps(a3hi, p1  );//pred1: x[2] = p1*x[3] - x[2]
	a2lo=_mm_sub_ps(t0lo, a2lo), a2hi=_mm_sub_ps(t0hi, a2hi);

	t0lo=_mm_mul_ps(a2lo, u1  ), t0hi=_mm_mul_ps(a2hi, u1  );//update1: x[3]-=u1*x[2]
	a3lo=_mm_sub_ps(a3lo, t0lo), a3hi=_mm_sub_ps(a3hi, t0hi);

	t0lo=_mm_mul_ps(a7lo, half), t0hi=_mm_mul_ps(a7hi, half);//x[4] = x[7]/2 - x[4]
	a4lo=_mm_sub_ps(t0lo, a4lo), a4hi=_mm_sub_ps(t0hi, a4hi);

#ifdef DCT_8x8_PS_SCALED
	a0lo=_mm_mul_ps(a0lo, g0), a0hi=_mm_mul_ps(a0hi, g0);
	a1lo=_mm_mul_ps(a1lo, g1), a1hi=_mm_mul_ps(a1hi, g1);
	a2lo=_mm_mul_ps(a2lo, g2), a2hi=_mm_mul_ps(a2hi, g2);
	a3lo=_mm_mul_ps(a3lo, g3), a3hi=_mm_mul_ps(a3hi, g3);
	a4lo=_mm_mul_ps(a4lo, g4), a4hi=_mm_mul_ps(a4hi, g4);
	a5lo=_mm_mul_ps(a5lo, g5), a5hi=_mm_mul_ps(a5hi, g5);
	a6lo=_mm_mul_ps(a6lo, g6), a6hi=_mm_mul_ps(a6hi, g6);
	a7lo=_mm_mul_ps(a7lo, g7), a7hi=_mm_mul_ps(a7hi, g7);
#endif
	
	_mm_store_ps(ptr   , a0lo), _mm_load_ps(ptr+ 4, a0hi);
	_mm_store_ps(ptr+ 8, a4lo), _mm_load_ps(ptr+12, a4hi);
	_mm_store_ps(ptr+16, a6lo), _mm_load_ps(ptr+20, a6hi);
	_mm_store_ps(ptr+24, a2lo), _mm_load_ps(ptr+28, a2hi);
	_mm_store_ps(ptr+32, a7lo), _mm_load_ps(ptr+36, a7hi);
	_mm_store_ps(ptr+40, a3lo), _mm_load_ps(ptr+44, a3hi);
	_mm_store_ps(ptr+48, a5lo), _mm_load_ps(ptr+52, a5hi);
	_mm_store_ps(ptr+56, a1lo), _mm_load_ps(ptr+60, a1hi);
#endif
}
void DCT3_8x8_ps(float *ptr)
{
	//INVERSE OF Loeffler's factorization using lifting		http://www.sfu.ca/~jiel/papers/c003-bindct-ieee.pdf
#ifdef __AVX__
	__m256
		half=_mm256_set1_ps(0.5f),
		p1=_mm256_set1_ps(13.f/32),//binDCT-L5
		u1=_mm256_set1_ps(11.f/32),
		p2=_mm256_set1_ps(19.f/64),
		u2=_mm256_set1_ps( 9.f/16),
		p3=_mm256_set1_ps(19.f/64),
		p4=_mm256_set1_ps( 3.f/32),
		u3=_mm256_set1_ps( 3.f/16),
		p5=_mm256_set1_ps( 3.f/32);
	
	__m256
		x0=_mm256_load_ps(ptr   ),
		x4=_mm256_load_ps(ptr+ 8),
		x6=_mm256_load_ps(ptr+16),
		x2=_mm256_load_ps(ptr+24),
		x7=_mm256_load_ps(ptr+32),
		x3=_mm256_load_ps(ptr+40),
		x5=_mm256_load_ps(ptr+48),
		x1=_mm256_load_ps(ptr+56);
	
	//new stage after each extra newline, to undo lifting:
	//- permuted load instead of store
	//- reverse operation order
	//- replace each operation with its inverse
	x4=_mm256_fmsub_ps(x7, half, x4);

	x7=_mm256_sub_ps(x7, x4);//add

	x4=_mm256_fmadd_ps(x6, half, x4);//_mm256_fnmadd_ps
	x7=_mm256_fmadd_ps(x5, half, x7);//_mm256_fnmadd_ps

	x1=_mm256_fmsub_ps(x0, half, x1);
	x3=_mm256_fmadd_ps(x2, u1, x3);//_mm256_fnmadd_ps
	x6=_mm256_sub_ps(x4, x6);
	x5=_mm256_sub_ps(x7, x5);

	x0=_mm256_sub_ps(x0, x1);//add
	x2=_mm256_fmsub_ps(x3, p1, x2);
	x7=_mm256_fmadd_ps(x4, p3, x7);//_mm256_fnmadd_ps
	x6=_mm256_fmadd_ps(x5, p5, x6);//_mm256_fnmadd_ps
	
	x0=_mm256_fmadd_ps(x3, half, x0);//_mm256_fnmadd_ps
	x1=_mm256_fmadd_ps(x2, half, x1);//_mm256_fnmadd_ps
	x4=_mm256_fnmadd_ps(x7, u2, x4);//_mm256_fmadd_ps
	x5=_mm256_fnmadd_ps(x6, u3, x5);//_mm256_fmadd_ps

	x3=_mm256_sub_ps(x0, x3);
	x2=_mm256_sub_ps(x1, x2);
	x7=_mm256_fmadd_ps(x4, p2, x7);//_mm256_fnmadd_ps
	x6=_mm256_fmadd_ps(x5, p4, x6);//_mm256_fnmadd_ps

	x0=_mm256_fmadd_ps(x7, half, x0);//_mm256_fnmadd_ps
	x1=_mm256_fmadd_ps(x6, half, x1);//_mm256_fnmadd_ps
	x2=_mm256_fmadd_ps(x5, half, x2);//_mm256_fnmadd_ps
	x3=_mm256_fmadd_ps(x4, half, x3);//_mm256_fnmadd_ps

	x7=_mm256_sub_ps(x0, x7);
	x6=_mm256_sub_ps(x1, x6);
	x5=_mm256_sub_ps(x2, x5);
	x4=_mm256_sub_ps(x3, x4);
	
	_mm256_store_ps(ptr   , x0);
	_mm256_store_ps(ptr+ 8, x1);
	_mm256_store_ps(ptr+16, x2);
	_mm256_store_ps(ptr+24, x3);
	_mm256_store_ps(ptr+32, x4);
	_mm256_store_ps(ptr+40, x5);
	_mm256_store_ps(ptr+48, x6);
	_mm256_store_ps(ptr+56, x7);
#else
#error TODO
#endif
}
#if 1
void DCTtest2()
{
	ALIGN(32) float buf[64], b2[64];
	
	float
		srcmin=0, srcmax=0,//src [0.000000 1.000000]		Y/Co/Cg	src [-1.000000 1.000061]
		dstmin=0, dstmax=0;//DCT [-4.664061 3.825694]				DCT [-7.948018 10.167978]
	for(int it=0;;++it)
	{
		int k;
		srand((unsigned)__rdtsc());

		for(k=0;k<64;++k)
			//buf[k]=roundf((float)rand()*(10000.f/RAND_MAX))*0.0001f;//[0, 1] 4 digits
			//buf[k]=(float)rand()/RAND_MAX;//[0, 1]
			buf[k]=(float)(rand()-(RAND_MAX>>1))/(RAND_MAX>>1);//[-1, 1]
		memcpy(b2, buf, 64*sizeof(float));
		for(k=0;k<64;++k)
		{
			if(!it&&!k||srcmin>b2[k])
				srcmin=b2[k];
			if(!it&&!k||srcmax<b2[k])
				srcmax=b2[k];
		}

		DCT2_8x8_ps(b2);
		transpose_block8x8_ps(b2);
		DCT2_8x8_ps(b2);
		
		for(k=0;k<64;++k)
		{
			if(!it&&!k||dstmin>b2[k])
				dstmin=b2[k];
			if(!it&&!k||dstmax<b2[k])
				dstmax=b2[k];
		}

		DCT3_8x8_ps(b2);
		transpose_block8x8_ps(b2);
		DCT3_8x8_ps(b2);

		printf("src [%f %f] DCT [%f %f]\t\t", srcmin, srcmax, dstmin, dstmax);
		//int nplaces=4;
		//print_fbuf(buf, 8, 8, nplaces);
		//print_fbuf(b2, 8, 8, nplaces);
		compare_bufs_ps(b2, buf, 8, 8, "DCT", 0);
		pause1();
	}

	printf("Done.\n");
	pause();
	exit(0);
}
#endif


//FCT-II/III, O(N*lg(N)), POT size (uses FFT size N)
void fft_dit_ps(float *re, float *im, int lgsize, int forward)//DIT vs DIF: simply reverse the order of stage matrices along with the permutation, and transpose the stages
{
	int size=1<<lgsize, sign=!forward-(forward!=0);
	
	for(int k=0;k<size;k+=2)//separate mul-free sum & difference (m==2)
	{
		float
			re0=re[k],
			im0=im[k],
			re1=re[k+1],
			im1=im[k+1];

		re[k]=re0+re1;
		im[k]=im0+im1;
		re[k+1]=re0-re1;
		im[k+1]=im0-im1;
	}
	for(int m=4;m<=size;m<<=1)
	{
		int m2=m>>1;
		float angle=(float)((2*M_PI)*sign/m);
		float rewm=(float)cosf(angle), imwm=(float)sinf(angle);//exp(-+ 2*pi*i*(N/m)/N)
		for(int k=0;k<size;k+=m)
		{

			float
				re0=re[k],
				im0=im[k],
				re1=re[k+m2],
				im1=im[k+m2];

			re[k]=re0+re1;
			im[k]=im0+im1;
			re[k+m2]=re0-re1;
			im[k+m2]=im0-im1;

			float rew=rewm, imw=imwm;
			for(int j=1;;)
			{
				int idx=k+j;
				re0=re[idx];
				im0=im[idx];
				re1=re[idx+m2];
				im1=im[idx+m2];
				
				float
					re2=rew*re1-imw*im1,
					im2=rew*im1+imw*re1;

				re[idx]=re0+re2;
				im[idx]=im0+im2;
				re[idx+m2]=re0-re2;
				im[idx+m2]=im0-im2;

				//nmuls+=4;//

				++j;
				if(j>=m2)
					break;

				re0=rewm*rew-imwm*imw;//w*=wm
				im0=rewm*imw+imwm*rew;
				rew=re0;
				imw=im0;

				//nmuls+=4;//
			}
		}
	}
}
void fft_dif_ps(float *re, float *im, int lgsize, int forward)
{
	int size=1<<lgsize, sign=!forward-(forward!=0);
	
	for(int m=size;m>2;m>>=1)
	{
		int m2=m>>1;
		float angle=(float)((2*M_PI)*sign/m);
		float rewm=(float)cosf(angle), imwm=(float)sinf(angle);//exp(-+ 2*pi*i*(N/m)/N)
		for(int k=0;k<size;k+=m)
		{
			float
				re0=re[k],
				im0=im[k],
				re1=re[k+m2],
				im1=im[k+m2];

			re[k]=re0+re1;
			im[k]=im0+im1;
			re[k+m2]=re0-re1;
			im[k+m2]=im0-im1;
			
			float rew=rewm, imw=imwm;
			for(int j=1;;)
			{
				int idx=k+j;
				float
					re0=re[idx],
					im0=im[idx],
					re1=re[idx+m2],
					im1=im[idx+m2];

				re[idx]=re0+re1;
				im[idx]=im0+im1;

				re0-=re1;
				im0-=im1;
				
				re[idx+m2]=rew*re0-imw*im0;
				im[idx+m2]=rew*im0+imw*re0;

				//nmuls+=4;//
				
				++j;
				if(j>=m2)
					break;
				
				re0=rewm*rew-imwm*imw;//w*=wm
				im0=rewm*imw+imwm*rew;
				rew=re0;
				imw=im0;

				//nmuls+=4;//
			}
		}
	}
	for(int k=0;k<size;k+=2)//separate mul-free sum & difference (m==2)
	{
		float
			re0=re[k],
			im0=im[k],
			re1=re[k+1],
			im1=im[k+1];

		re[k]=re0+re1;
		im[k]=im0+im1;
		re[k+1]=re0-re1;
		im[k+1]=im0-im1;
	}
}

//A fast cosine transform (FCT) in one and two dimensions - page 103
#if 0
typedef struct FCT1D_PS_Params
{
	int lgsize;
	int *fctp;
	float *re, *im, *rew4N, *imw4N;
} FCT1D_PS_Params;
#endif
#define BIT_REVERSE(N, NBITS)\
	N=N<<1&0xAAAA|N>>1&0x5555;\
	N=N<<2&0xCCCC|N>>2&0x3333;\
	N=N<<4&0xF0F0|N>>4&0x0F0F;\
	N=N<<8&0xFF00|N>>8&0x00FF;\
	N>>=16-NBITS;
void FCT1D_ps_gen(int lgsize, FCT1D_PS_Params *p)//up to size 65536
{
	int size=1<<lgsize, sh=16-lgsize, mask=(1<<lgsize)-1;
	p->lgsize=lgsize;
	p->fctp=(int*)malloc(size*sizeof(int));
	p->re=(float*)malloc(size*sizeof(float));
	p->im=(float*)malloc(size*sizeof(float));
	p->rew4N=(float*)malloc(size*sizeof(float));
	p->imw4N=(float*)malloc(size*sizeof(float));
	for(int k=0;k<size;++k)//FCT permutation: f(n) = (bit_reverse(n>>1)^-(n&1)) & ((1<<lgsize)-1)
	{
		unsigned val=k>>1;
		BIT_REVERSE(val, lgsize)
		val^=-(k&1);
		val&=mask;
		p->fctp[k]=val;
	}
	float invsqrtN=(float)(1/sqrt(size));
	for(int k=0;k<size;++k)
	{
		float temp=(float)(-2*M_PI*k/((size_t)size<<2));
		p->rew4N[k]=invsqrtN*cosf(temp);
		p->imw4N[k]=invsqrtN*sinf(temp);
	}
}
void FCT1D_ps_free(FCT1D_PS_Params *p)
{
	free(p->fctp);
	free(p->re);
	free(p->im);
	free(p->rew4N);
	free(p->imw4N);
	memset(p, 0, sizeof(*p));
}
void FCT1D_ps_fwd(FCT1D_PS_Params *p, float *data, int stride)//data, re, im, rew, imw are all of size (1<<lgsize), result is in data buffer
{
	int size=1<<p->lgsize;

	for(int ks=0, kd=0;kd<size;ks+=stride, ++kd)
		p->re[p->fctp[kd]]=data[ks];
	memset(p->im, 0, size*sizeof(float));

	fft_dit_ps(p->re, p->im, p->lgsize, 1);//V[k] = sum n=0 to N-1: v[n] WN^nk

	for(int ks=0, kd=0;ks<size;++ks, kd+=stride)//C[k] = 2 Re{W4N^k V[k]}
	{
		float r=p->re[ks], i=p->im[ks];
		//data[kd]=p->rew4N[ks]*r-p->imw4N[ks]*i;
		data[kd]=2*(p->rew4N[ks]*r-p->imw4N[ks]*i);
		//imaginary part would be 2(imw4N[ks]*r+rew4N[ks]*i) = data[(N*stride-kd)%N*stride] or data[(N-1)*stride-kd] ?
	}
}
void FCT1D_ps_inv(FCT1D_PS_Params *p, float *data, int stride)//data, re, im, rew, imw are all of size (1<<lgsize), result is in data buffer
{
	int size=1<<p->lgsize, offset=size*stride;

	for(int ks=0, kd=0;kd<size;ks+=stride, ++kd)
	{
		float r=data[ks], i=ks?-data[offset-ks]:0;//V[k] = (1/2) W4N^-k (C[k] - i*C[N-k])
		//p->re[kd]= p->rew4N[kd]*r-p->imw4N[kd]*i;
		//p->im[kd]=-p->rew4N[kd]*i-p->imw4N[kd]*r;
		p->re[kd]=0.5f*(p->rew4N[kd]*r+p->imw4N[kd]*i);	//(1/2) (cos(2pik/4N) *  C[k]   - -sin(2pik/4N) *(-C[N-k]))
		p->im[kd]=0.5f*(p->rew4N[kd]*i-p->imw4N[kd]*r);	//(1/2) (cos(2pik/4N) *(-C[N-k])+ -sin(2pik/4N) *  C[k])
	}

	fft_dif_ps(p->re, p->im, p->lgsize, 0);

	for(int ks=0, kd=0;kd<size;ks+=stride, ++kd)
		data[p->fctp[kd]]=p->re[ks];
}


void printmat_pd(double *matrix, int bw, int bh, int prec)
{
	for(int ky=0;ky<bh;++ky)
	{
		for(int kx=0;kx<bw;++kx)
			printf("%*.*lf", prec+5, prec, matrix[bw*ky+kx]);
		printf("\n");
	}
	printf("\n");
}
void matmul_pd(double *res, const double *A, const double *B, int h1, int w1h2, int w2)
{
	for(int ky=0;ky<h1;++ky)
	{
		for(int kx=0;kx<w2;++kx)
		{
			double sum=0;
			for(int ki=0;ki<w1h2;++ki)
				sum+=A[w1h2*ky+ki]*B[w2*ki+kx];
			res[w2*ky+kx]=sum;
		}
	}
}
void calc_DCT_matrix(int lgn, double *DCT2, double *DCT3)//arrays size (1<<(lgn<<1))
{
	int n=1<<lgn;
	double gain=sqrt(2./n), pi_n=M_PI/n;
	for(int ko=0;ko<n;++ko)
	{
		for(int ki=0;ki<n;++ki)
		{
			DCT2[ko<<2|ki]=gain*cos(pi_n*(ki+0.5)*ko);
			if(ki)
				DCT3[ko<<2|ki]=gain*cos(pi_n*(ko+0.5)*ki);
			else
				DCT3[ko<<2|ki]=gain*0.5;
		}
	}
}
#define lgN 1
#define N (1<<lgN)
double DCT2[N*N], DCT3[N*N], prod[N*N];
void test4()
{
	calc_DCT_matrix(lgN, DCT2, DCT3);

	printf("DCT2 (fwd):\n");
	printmat_pd(DCT2, N, N, 10);

	printf("DCT3 (inv):\n");
	printmat_pd(DCT3, N, N, 10);

	matmul_pd(DCT3, DCT2, prod, N, N, N);

	printf("Product:\n");
	printmat_pd(prod, N, N, 10);

	printf("Done.\n");
	pause();
	exit(0);
}