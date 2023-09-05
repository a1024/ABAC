#include"lossy.h"
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
#include<immintrin.h>
static const char file[]=__FILE__;


#ifndef __AVX__
#define __AVX__ 1
#endif

void measure_distortion(const unsigned char *b1, const unsigned char *b2, int iw, int ih, double *ret_rmse, double *ret_psnr)
{
	int res=iw*ih;
	long long sum[3]={0};
	for(int k=0;k<res;++k)
	{
		int dr=b1[k<<2  ]-b2[k<<2  ],
			dg=b1[k<<2|1]-b2[k<<2|1],
			db=b1[k<<2|2]-b2[k<<2|2];
		sum[0]+=dr*dr;
		sum[1]+=dg*dg;
		sum[2]+=db*db;
	}
	double rmse[]=
	{
		sqrt((double)sum[0]/res),
		sqrt((double)sum[1]/res),
		sqrt((double)sum[2]/res),
		sqrt((double)(sum[0]+sum[1]+sum[2])/(res*3)),
	};
	double psnr[]=
	{
		20*log10(255/rmse[0]),
		20*log10(255/rmse[1]),
		20*log10(255/rmse[2]),
		20*log10(255/rmse[3]),
	};
	if(ret_rmse)
		memcpy(ret_rmse, rmse, 4*sizeof(double));
	if(ret_psnr)
		memcpy(ret_psnr, psnr, 4*sizeof(double));
}

void cvt_u8_ps(float *dst, const unsigned char *src, int res, int srcstride, int dststride)
{
	for(int ks=0, kd=0, slen=res*srcstride;ks<slen;ks+=srcstride, kd+=dststride)
		dst[kd]=(float)(src[ks]-128);
}
void cvt_ps_u8(unsigned char *dst, const float *src, int res, int srcstride, int dststride)
{
	for(int ks=0, kd=0, slen=res*srcstride;ks<slen;ks+=srcstride, kd+=dststride)
	{
		float val=src[ks];
		val+=128;
		CLAMP(0, val, 255);
		dst[kd]=(unsigned char)val;
	}
}

void dnsample(const float *src, int iw, int ih, float *dst)//src & dst can be the same buffer
{
	for(int ky=0;ky<ih;ky+=2)
	{
		for(int kx=0;kx<iw;kx+=2)
		{
			int idx=iw*ky+kx;
			float val=src[idx]+src[idx+1]+src[idx+iw]+src[idx+iw+1];
			val*=0.25f;
			dst[(iw>>1)*(ky>>1)+(kx>>1)]=val;
		}
	}
}
void upsample(const float *src, int iw, int ih, float *dst)//src & dst CANNOT be the same buffer
{
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
			float v0=src[idx];
			float v1=kx+1<iw?src[idx+1]:v0;
			float v2=ky+1<ih?src[idx+iw]:v0;
			float v3=kx+1<iw&&ky+1<ih?src[idx+iw+1]:v0;
			idx=(iw<<1)*(ky<<1)+(kx<<1);
			dst[idx]=v0;
			dst[idx+1]=(v0+v1)*0.5f;
			dst[idx+(iw<<1)]=(v0+v2)*0.5f;
			dst[idx+(iw<<1)+1]=(v0+v1+v2+v3)*0.25f;
		}
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
		b+=g>>1;

		buf[k  ]=r;//Co
		buf[k|1]=g;//Cg
		buf[k|2]=b;//Y
	}
}
void colortransform_ycocg_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];
		
		b-=g>>1;
		g+=b;
		b-=r>>1;
		r+=b;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_ycocb_fwd(char *buf, int iw, int ih)//YCoCb is YCoCg but with swapped green and blue
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];
		
		//Y    1/4   1/4  1/2   r
		//Co   1    -1    0     g
		//Cb  -1/2   1   -1/2   b

		r-=g;       //diff(r, g) = Co
		g+=r>>1;    //g+floor((r-g)/2) = av(r, g)
		b-=g;       //b-av(r, g) = Cb
		g+=b>>1;    //av(r, g)+floor((b-av(r, g))/2) = av(av(r, g), b) = Y

		buf[k  ]=r;//Co
		buf[k|1]=g;//Y
		buf[k|2]=b;//Cb
	}
}
void colortransform_ycocb_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		//char r=buf[k|1], g=buf[k], b=buf[k|2];//Y Co Cb
		//char r=buf[k], g=buf[k|2], b=buf[k|1];//Co Cb Y
		char r=buf[k], g=buf[k|1], b=buf[k|2];//Co Y Cb	original
		
		g-=b>>1;
		b+=g;
		g-=r>>1;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}

//#define YCOCB_FIX 0//4
void ycocb_fwd_subsample_separate(const unsigned char *buf, int iw, int ih, short *luma, short *co, short *cb)
{
	int ky, kx, idx;
	short p[16]={0};
	int w2=(iw+1)>>1;
	for(ky=0;ky<ih;ky+=2)
	{
		for(kx=0;kx<iw;kx+=2)
		{
			idx=iw*ky+kx;
			
			for(int k=0;k<3;++k)
				p[k]=(char)(buf[(idx<<2)+k]-128);
			if(kx+1<iw)
			{
				for(int k=4;k<7;++k)
					p[k]=(char)(buf[(idx<<2)+k]-128);
			}
			else
				memcpy(p+4, p, 4);
			if(ky+1<ih)
			{
				for(int k=0;k<3;++k)
					p[8+k]=(char)(buf[((idx+iw)<<2)+k]-128);
				if(kx+1<iw)
				{
					for(int k=4;k<7;++k)
						p[8+k]=(char)(buf[((idx+iw)<<2)+k]-128);
				}
				else
					memcpy(p+12, p+8, 4);
			}
			else
				memcpy(p+8, p, 8);
			//for(int k=0;k<8;++k)
			//	p[k]=(char)(buf[(idx<<2)+k]-128)<<4;
			//for(int k=0;k<8;++k)
			//	p[k+8]=(char)(buf[((idx+iw)<<2)+k]-128)<<4;
			//char p[16];
			//memcpy(p, buf+(idx<<2), 8);
			//memcpy(p, buf+((idx+iw)<<2), 8);

			for(int k=0;k<16;k+=4)
			{
				p[k]-=p[k|1];       //r-=g      Co	9 bit	fwd
				p[k|1]+=p[k]>>1;    //g+=r>>1
				p[k|2]-=p[k|1];     //b-=g      Cb	9 bit
				p[k|1]+=p[k|2]>>1;  //g+=b>>1   Y	8 bit
			}

			if((unsigned)(p[1]+128)>=256||(unsigned)(p[4+1]+128)>=256||(unsigned)(p[8+1]+128)>=256||(unsigned)(p[12+1]+128)>=256)//
				LOG_ERROR("Range error");//

			luma[idx]=p[1];
			if(kx+1<iw)
				luma[idx+1]=p[4+1];
			if(ky+1<ih)
			{
				luma[idx+iw]=p[8+1];
				if(kx+1<iw)
					luma[idx+iw+1]=p[12+1];
			}
			idx=w2*(ky>>1)+(kx>>1);
			co[idx]=(p[0  ]+p[4  ]+p[8  ]+p[12  ])>>2;
			cb[idx]=(p[0+2]+p[4+2]+p[8+2]+p[12+2])>>2;
		}
	}
}
void ycocb_inv_upsample_separate(const short *luma, const short *co, const short *cb, int iw, int ih, unsigned char *buf)
{
	int ky, kx, idx;
	short p[16]={0};
	int w2=(iw+1)>>1, h2=(ih+1)>>1;
	for(ky=0;ky<ih;ky+=2)
	{
		for(kx=0;kx<iw;kx+=2)
		{
			idx=w2*(ky>>1)+(kx>>1);
			p[ 0]=co[idx];
			p[ 4]=(kx>>1)+1<w2              ?(co[idx   +1]+p[0])>>1:p[0];
			p[ 8]=              (ky>>1)+1<h2?(co[idx+iw  ]+p[0])>>1:p[0];
			p[12]=(kx>>1)+1<w2&&(ky>>1)+1<h2?(co[idx+iw+1]+p[0])>>1:p[0];
			p[ 0+2]=cb[idx];
			p[ 4+2]=(kx>>1)+1<w2              ?(cb[idx   +1]+p[2])>>1:p[2];
			p[ 8+2]=              (ky>>1)+1<h2?(cb[idx+iw  ]+p[2])>>1:p[2];
			p[12+2]=(kx>>1)+1<w2&&(ky>>1)+1<h2?(cb[idx+iw+1]+p[2])>>1:p[2];
			idx=iw*ky+kx;
			p[0+1]=luma[idx];
			if(kx+1<iw)
				p[4+1]=luma[idx+1];
			else
				p[4+1]=p[1];
			if(ky+1<ih)
			{
				p[8+1]=luma[idx+iw];
				if(kx+1<iw)
					p[12+1]=luma[idx+iw+1];
				else
					p[12+1]=p[8+1];
			}
			else
				p[8+1]=p[1], p[12+1]=p[4+1];
			for(int k=0;k<16;k+=4)
			{
				p[k|1]-=p[k|2]>>1;  //g-=b>>1   Y		inv
				p[k|2]+=p[k|1];     //b+=g      Cb
				p[k|1]-=p[k]>>1;    //g-=r>>1
				p[k]+=p[k|1];       //r+=g      Co
			}
			for(int k=0;k<3;++k)
				buf[idx<<2|k]=CLAMP(-128, p[k], 127)+128;
			buf[idx<<2|3]=0xFF;

			if(kx+1<iw)
			{
				for(int k=0;k<3;++k)
					buf[(idx+1)<<2|k]=CLAMP(-128, p[4|k], 127)+128;
				buf[(idx+1)<<2|3]=0xFF;
			}
			if(ky+1<ih)
			{
				for(int k=0;k<3;++k)
					buf[(idx+iw)<<2|k]=CLAMP(-128, p[8|k], 127)+128;
				buf[(idx+iw)<<2|3]=0xFF;
				if(kx+1<iw)
				{
					for(int k=0;k<3;++k)
						buf[(idx+iw+1)<<2|k]=CLAMP(-128, p[12|k], 127)+128;
					buf[(idx+iw+1)<<2|3]=0xFF;
				}
			}
		}
	}
}
void ycocb_inv_upsample_i8(const short *luma, const short *co, const short *cb, int iw, int ih, unsigned char *buf)
{
	int ky, kx, idx;
	short p[16]={0};
	int w2=(iw+1)>>1, h2=(ih+1)>>1;
	for(ky=0;ky<ih;ky+=2)
	{
		for(kx=0;kx<iw;kx+=2)
		{
			idx=w2*(ky>>1)+(kx>>1);
			p[ 0]=co[idx];
			p[ 4]=(kx>>1)+1<w2              ?(co[idx   +1]+p[0])>>1:p[0];
			p[ 8]=              (ky>>1)+1<h2?(co[idx+iw  ]+p[0])>>1:p[0];
			p[12]=(kx>>1)+1<w2&&(ky>>1)+1<h2?(co[idx+iw+1]+p[0])>>1:p[0];
			p[ 0+2]=cb[idx];
			p[ 4+2]=(kx>>1)+1<w2              ?(cb[idx   +1]+p[2])>>1:p[2];
			p[ 8+2]=              (ky>>1)+1<h2?(cb[idx+iw  ]+p[2])>>1:p[2];
			p[12+2]=(kx>>1)+1<w2&&(ky>>1)+1<h2?(cb[idx+iw+1]+p[2])>>1:p[2];
			idx=iw*ky+kx;
			p[0+1]=luma[idx];
			if(kx+1<iw)
				p[4+1]=luma[idx+1];
			else
				p[4+1]=p[1];
			if(ky+1<ih)
			{
				p[8+1]=luma[idx+iw];
				if(kx+1<iw)
					p[12+1]=luma[idx+iw+1];
				else
					p[12+1]=p[8+1];
			}
			else
				p[8+1]=p[1], p[12+1]=p[4+1];
			for(int k=0;k<16;k+=4)
			{
				p[k|1]-=p[k|2]>>1;  //g-=b>>1   Y		inv
				p[k|2]+=p[k|1];     //b+=g      Cb
				p[k|1]-=p[k]>>1;    //g-=r>>1
				p[k]+=p[k|1];       //r+=g      Co
			}
			for(int k=0;k<3;++k)
				buf[idx<<2|k]=CLAMP(-128, p[k], 127)+128;
			buf[idx<<2|3]=0xFF;

			if(kx+1<iw)
			{
				for(int k=0;k<3;++k)
					buf[(idx+1)<<2|k]=CLAMP(-128, p[4|k], 127)+128;
				buf[(idx+1)<<2|3]=0xFF;
			}
			if(ky+1<ih)
			{
				for(int k=0;k<3;++k)
					buf[(idx+iw)<<2|k]=CLAMP(-128, p[8|k], 127)+128;
				buf[(idx+iw)<<2|3]=0xFF;
				if(kx+1<iw)
				{
					for(int k=0;k<3;++k)
						buf[(idx+iw+1)<<2|k]=CLAMP(-128, p[12|k], 127)+128;
					buf[(idx+iw+1)<<2|3]=0xFF;
				}
			}
		}
	}
}

void ycocb_fwd_subsample_ps(const unsigned char *buf, int iw, int ih, float *luma, float *co, float *cb)
{
	int ky, kx, idx;
	float p[16]={0};
	int w2=(iw+1)>>1;
	for(ky=0;ky<ih;ky+=2)
	{
		for(kx=0;kx<iw;kx+=2)
		{
			//if(kx==(iw>>1)&&ky==(ih>>1))
			//	printf("");

			idx=iw*ky+kx;
			
			for(int k=0;k<3;++k)
				p[k]=(char)(buf[(idx<<2)+k]-128);
			if(kx+1<iw)
			{
				for(int k=4;k<7;++k)
					p[k]=(char)(buf[(idx<<2)+k]-128);
			}
			else
				memcpy(p+4, p, 4*sizeof(float));
			if(ky+1<ih)
			{
				for(int k=0;k<3;++k)
					p[8|k]=(char)(buf[((idx+iw)<<2)+k]-128);
				if(kx+1<iw)
				{
					for(int k=4;k<7;++k)
						p[8|k]=(char)(buf[((idx+iw)<<2)+k]-128);
				}
				else
					memcpy(p+12, p+8, 4*sizeof(float));
			}
			else
				memcpy(p+8, p, 8*sizeof(float));

			for(int k=0;k<16;k+=4)
			{
				//float p0[16];//
				//memcpy(p0, p, 16*sizeof(float));//
				//memcpy(p, p0, 16*sizeof(float));//

				//if(p[k]<-1||p[k]>1||p[k|1]<-1||p[k|1]>1||p[k|2]<-1||p[k|2]>1)
				//	LOG_ERROR("Overflow");
				p[k]-=p[k|1];       //r-=g      Co	9 bit	fwd
				p[k|1]+=p[k]*0.5f;  //g+=r>>1
				p[k|2]-=p[k|1];     //b-=g      Cb	9 bit
				p[k|1]+=p[k|2]*0.5f;//g+=b>>1   Y	8 bit

				//if(p[k]<-1||p[k]>1||p[k|1]<-1||p[k|1]>1||p[k|2]<-1||p[k|2]>1)
				//	LOG_ERROR("Overflow");
			}

			//if((unsigned)(p[1]+128)>=256||(unsigned)(p[4+1]+128)>=256||(unsigned)(p[8+1]+128)>=256||(unsigned)(p[12+1]+128)>=256)//
			//	LOG_ERROR("Range error");//

			luma[idx]=p[1];
			if(kx+1<iw)
				luma[idx+1]=p[4|1];
			if(ky+1<ih)
			{
				luma[idx+iw]=p[8|1];
				if(kx+1<iw)
					luma[idx+iw+1]=p[12|1];
			}
			idx=w2*(ky>>1)+(kx>>1);
			co[idx]=(p[0  ]+p[4  ]+p[8  ]+p[12  ])*0.25f;
			cb[idx]=(p[0|2]+p[4|2]+p[8|2]+p[12|2])*0.25f;
		}
	}
}
void ycocb_inv_upsample_ps(const float *luma, const float *co, const float *cb, int iw, int ih, unsigned char *buf)
{
	int ky, kx, idx;
	float p[16]={0};
	int w2=(iw+1)>>1, h2=(ih+1)>>1;
	for(ky=0;ky<ih;ky+=2)
	{
		for(kx=0;kx<iw;kx+=2)
		{
			//if(kx==(iw>>1)&&ky==(ih>>1))
			//	printf("");

			idx=w2*(ky>>1)+(kx>>1);
			p[ 0]=co[idx];
			p[ 4]=(kx>>1)+1<w2              ?(co[idx   +1]+p[0])*0.5f:p[0];
			p[ 8]=              (ky>>1)+1<h2?(co[idx+w2  ]+p[0])*0.5f:p[0];
			p[12]=(kx>>1)+1<w2&&(ky>>1)+1<h2?(co[idx+w2+1]+p[0]+p[4]+p[8])*0.25f:p[0];
			p[ 0|2]=cb[idx];
			p[ 4|2]=(kx>>1)+1<w2              ?(cb[idx   +1]+p[2])*0.5f:p[2];
			p[ 8|2]=              (ky>>1)+1<h2?(cb[idx+w2  ]+p[2])*0.5f:p[2];
			p[12|2]=(kx>>1)+1<w2&&(ky>>1)+1<h2?(cb[idx+w2+1]+p[2]+p[4|2]+p[8|2])*0.25f:p[2];
			idx=iw*ky+kx;
			p[ 0|1]=luma[idx];
			p[ 4|1]=kx+1<iw?luma[idx+1]:p[1];
			p[ 8|1]=ky+1<ih?luma[idx+iw]:p[1];
			p[12|1]=kx+1<iw&&ky+1<ih?luma[idx+iw+1]:p[1];

			for(int k=0;k<16;k+=4)
			{
				//float p0[16];//
				//memcpy(p0, p, 16*sizeof(float));//
				//memcpy(p, p0, 16*sizeof(float));//

				p[k|1]-=p[k|2]*0.5f;//g-=b>>1   Y		inv
				p[k|2]+=p[k|1];     //b+=g      Cb
				p[k|1]-=p[k]*0.5f;  //g-=r>>1
				p[k]+=p[k|1];       //r+=g      Co

				//if(p[k]<-1||p[k]>1||p[k|1]<-1||p[k|1]>1||p[k|2]<-1||p[k|2]>1)
				//	LOG_ERROR("Overflow");
			}
			for(int k=0;k<3;++k)
				buf[idx<<2|k]=CLAMP(-128, p[k], 127)+128;
			//{
			//	//if(p[k]<-1||p[k]>1)
			//	//	LOG_ERROR("Overflow");
			//	buf[idx<<2|k]=(unsigned char)((CLAMP(-1, p[k], 1)+1)*127.5f);
			//}
			buf[idx<<2|3]=0xFF;

			if(kx+1<iw)
			{
				for(int k=0;k<3;++k)
					buf[(idx+1)<<2|k]=CLAMP(-128, p[4|k], 127)+128;
				buf[(idx+1)<<2|3]=0xFF;
			}
			if(ky+1<ih)
			{
				for(int k=0;k<3;++k)
					buf[(idx+iw)<<2|k]=CLAMP(-128, p[8|k], 127)+128;
				buf[(idx+iw)<<2|3]=0xFF;
				if(kx+1<iw)
				{
					for(int k=0;k<3;++k)
						buf[(idx+iw+1)<<2|k]=CLAMP(-128, p[12|k], 127)+128;
					buf[(idx+iw+1)<<2|3]=0xFF;
				}
			}
		}
	}
}

#if 0
void colortransform_ycocb_ps_fwd(float *c0, float *c1, float *c2, int iw, int ih, int stride)
{
	for(ptrdiff_t k=0, res=(ptrdiff_t)iw*ih, idx=0;k<res;++k, idx+=stride)
	{
		float r=c0[idx], g=c1[idx], b=c2[idx];

		r-=g;
		g+=r*0.5f;
		b-=g;
		g+=b*0.5f;

		c0[idx]=r;//Co
		c1[idx]=g;//Y
		c2[idx]=b;//Cb
	}
}
void colortransform_ycocb_ps_inv(float *c0, float *c1, float *c2, int iw, int ih, int stride)
{
	for(ptrdiff_t k=0, res=(ptrdiff_t)iw*ih, idx=0;k<res;++k, idx+=stride)
	{
		float r=c0[idx], g=c1[idx], b=c2[idx];//Co Y Cb
		
		g-=b*0.5f;
		b+=g;
		g-=r*0.5f;
		r+=g;

		c0[idx]=r;
		c1[idx]=g;
		c2[idx]=b;
	}
}
#endif

//	#define DCT_8x8_PS_SCALED

//fast 8 point DCT-II/III, call transpose then call again, ptr to 64 packed floats is aligned by 32 bytes
void DCT2_8x8_ps(float *ptr)
{
	//Loeffler's factorization using lifting	https://thanglong.ece.jhu.edu/Tran/Pub/bindct-IEEESP.pdf	http://www.sfu.ca/~jiel/papers/c003-bindct-ieee.pdf
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
#if 0
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
#if 0
void DCT2_8x8_ps_buf(float *buf, int iw, int ih)
{
	ALIGN(32) float temp[64];
	for(int ky=0;ky<ih-7;ky+=8)
	{
		for(int kx=0;kx<iw-7;kx+=8)
		{
			int idx=iw*ky+kx;
			for(int k=0;k<8;++k)
				memcpy(temp+(k<3), buf+idx+iw*k, 8*sizeof(float));

			DCT2_8x8_ps(temp);
			transpose_block8x8_ps(temp);
			DCT2_8x8_ps(temp);

			for(int k=0;k<8;++k)
				memcpy(buf+idx+iw*k, temp+(k<3), 8*sizeof(float));
		}
	}
}
void DCT3_8x8_ps_buf(float *buf, int iw, int ih)
{
	ALIGN(32) float temp[64];
	for(int ky=0;ky<ih-7;ky+=8)
	{
		for(int kx=0;kx<iw-7;kx+=8)
		{
			int idx=iw*ky+kx;
			for(int k=0;k<8;++k)
				memcpy(temp+(k<3), buf+idx+iw*k, 8*sizeof(float));

			DCT3_8x8_ps(temp);
			transpose_block8x8_ps(temp);
			DCT3_8x8_ps(temp);

			for(int k=0;k<8;++k)
				memcpy(buf+idx+iw*k, temp+(k<3), 8*sizeof(float));
		}
	}
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

	for(int k=0;k<size;++k)
		data[p->fctp[k]*stride]=p->re[k];
	//for(int ks=0, kd=0;kd<size;ks+=stride, ++kd)//X
	//	data[p->fctp[kd]]=p->re[ks];
}


void DCT2_8x8_ps_buf(float *buf, int iw, int ih)
{
	FCT1D_PS_Params params;
	FCT1D_ps_gen(3, &params);
	ALIGN(32) float temp[64];
	for(int ky=0;ky<ih-7;ky+=8)
	{
		for(int kx=0;kx<iw-7;kx+=8)
		{
			int idx=iw*ky+kx;
			for(int k=0;k<8;++k)
				memcpy(temp+(k<<3), buf+idx+iw*k, 8*sizeof(float));
			
			for(int k=0;k<8;++k)
				FCT1D_ps_fwd(&params, temp+(k<<3), 1);
			transpose_block8x8_ps(temp);
			for(int k=0;k<8;++k)
				FCT1D_ps_fwd(&params, temp+(k<<3), 1);

			for(int k=0;k<8;++k)
				memcpy(buf+idx+iw*k, temp+(k<<3), 8*sizeof(float));
		}
	}
	FCT1D_ps_free(&params);
}
void DCT3_8x8_ps_buf(float *buf, int iw, int ih)
{
	FCT1D_PS_Params params;
	FCT1D_ps_gen(3, &params);
	ALIGN(32) float temp[64];
	for(int ky=0;ky<ih-7;ky+=8)
	{
		for(int kx=0;kx<iw-7;kx+=8)
		{
			int idx=iw*ky+kx;
			for(int k=0;k<8;++k)
				memcpy(temp+(k<<3), buf+idx+iw*k, 8*sizeof(float));
			
			for(int k=0;k<8;++k)
				FCT1D_ps_inv(&params, temp+(k<<3), 1);
			transpose_block8x8_ps(temp);
			for(int k=0;k<8;++k)
				FCT1D_ps_inv(&params, temp+(k<<3), 1);

			for(int k=0;k<8;++k)
				memcpy(buf+idx+iw*k, temp+(k<<3), 8*sizeof(float));
		}
	}
	FCT1D_ps_free(&params);
}


static void DCT_col_SSE2(short *data)//google/webmproject/sjpeg/fdct.cc
{
	__m128i m0, m1, m2, m3, m4, m5, m6, m7;

#define LOAD(REG, ADDR) REG=_mm_load_si128((__m128i*)(ADDR))
#define STORE(ADDR, REG) _mm_store_si128((__m128i*)(ADDR), REG)
#define LOAD1(REG, VAL) REG=_mm_set1_epi16(VAL)
#define ADD(A, B) A=_mm_add_epi16(A, B)
#define SUB(A, B) A=_mm_sub_epi16(A, B)
#define BUTTERFLY(A, B)\
	SUB(A, B),\
	ADD(B, B),\
	ADD(B, A)
#define SLLI(REG, SL) REG=_mm_slli_epi16(REG, SL)
#define MUL(A, B) A=_mm_mulhi_epi16(A, B)

	LOAD(m0, data);
	LOAD(m2, data+2*8);
	LOAD(m7, data+7*8);
	LOAD(m5, data+5*8);

	BUTTERFLY(m0, m7);
	BUTTERFLY(m2, m5);

	LOAD(m3, data+3*8);
	LOAD(m4, data+4*8);
	BUTTERFLY(m3, m4);

	LOAD(m6, data+6*8);
	LOAD(m1, data+1*8);
	BUTTERFLY(m1, m6);

	BUTTERFLY(m7, m4);
	BUTTERFLY(m6, m5);

	SLLI(m4, 3);//multiply by 8 to use mulhi here and in row DCT
	SLLI(m5, 3);

	BUTTERFLY(m4, m5);
	STORE(data+0*8, m5);
	STORE(data+4*8, m4);
	
	SLLI(m7, 3);
	SLLI(m6, 3);
	SLLI(m3, 3);
	SLLI(m0, 3);

	LOAD1(m4, 27146);//tan(2*pi/16) = sqrt(2)-1
	m5=m4;
	MUL(m4, m7);
	MUL(m5, m6);
	SUB(m4, m6);
	ADD(m5, m7);
	STORE(data+2*8, m5);
	STORE(data+6*8, m4);

	LOAD1(m6, 23170);//1/(2sqrt(2))		half of 1/sqrt2, to fit in 15bit
	SLLI(m2, 3+1);
	SLLI(m1, 3+1);
	BUTTERFLY(m1, m2);
	MUL(m2, m6);
	MUL(m1, m6);
	BUTTERFLY(m3, m1);
	BUTTERFLY(m0, m2);

	LOAD1(m4, -21746);//tan(3*pi/16)-1
	LOAD1(m5, 13036);//tan(pi/16)
	m7=m3;
	m6=m1;
	MUL(m3, m4);
	MUL(m1, m5);

	ADD(m3, m7);
	ADD(m1, m2);

	__m128i bit=_mm_set1_epi16(1);//correct LSB
	ADD(m1, bit);
	ADD(m3, bit);

	MUL(m4, m0);
	MUL(m5, m2);
	ADD(m4, m0);
	SUB(m0, m3);
	ADD(m7, m4);
	SUB(m5, m6);

	STORE(data+1*8, m1);
	STORE(data+3*8, m0);
	STORE(data+5*8, m7);
	STORE(data+7*8, m5);

#undef LOAD
#undef STORE
#undef LOAD1
#undef ADD
#undef SUB
#undef BUTTERFLY
#undef SLLI
#undef MUL
}
static void DCT_row_SSE2(short *data, const short *table1, const short *table2)//google/webmproject/sjpeg/fdct.cc
{
#if defined __AVX__ && 0
	__m256i m0=_mm256_load_si256((__m256i*)data);//TODO process 4 rows at a time
#else
	//process two rows in parallel
	__m128i m0=_mm_load_si128((__m128i*)(data+0*8));//x0 ~ x15
	__m128i m2=_mm_load_si128((__m128i*)(data+1*8));
	m0=_mm_shufflehi_epi16(m0, _MM_SHUFFLE(0, 1, 2, 3));//[ 0  1  2  3| 7  6  5  4]
	m2=_mm_shufflehi_epi16(m2, _MM_SHUFFLE(0, 1, 2, 3));//[ 8  9 10 11|15 14 13 12]
	__m128i m4=m0;
	m0=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(m0), _mm_castsi128_ps(m2), _MM_SHUFFLE(1, 0, 1, 0)));//[ 0  1  2  3| 8  9 10 11]
	m4=_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(m4), _mm_castsi128_ps(m2), _MM_SHUFFLE(3, 2, 3, 2)));//[ 7  6  5  4|15 14 13 12]
	
	//initial butterfly
	m2=m0;
	m0=_mm_add_epi16(m0, m4);// a0=x0+x7 | a1=x1+x6 | a2=x2+x5 | a3=x3+x4
	m2=_mm_sub_epi16(m2, m4);// b0=x0-x7 | b1=x1-x6 | b2=x2-x5 | b3=x3-x4

	__m128i m6;
	m4=m0;
	m0=_mm_unpacklo_epi32(m0, m2);
	m4=_mm_unpackhi_epi32(m4, m2);
	m2=_mm_shuffle_epi32(m0, _MM_SHUFFLE(1, 0, 3, 2));
	m6=_mm_shuffle_epi32(m4, _MM_SHUFFLE(1, 0, 3, 2));
	
  __m128i m1, m3, m5, m7;
  m1=_mm_madd_epi16(m2, _mm_load_si128((const __m128i*)table1+1));
  m3=_mm_madd_epi16(m0, _mm_load_si128((const __m128i*)table1+2));
  m5=_mm_madd_epi16(m6, _mm_load_si128((const __m128i*)table2+1));
  m7=_mm_madd_epi16(m4, _mm_load_si128((const __m128i*)table2+2));

  m2=_mm_madd_epi16(m2, _mm_load_si128((const __m128i*)table1+3));
  m0=_mm_madd_epi16(m0, _mm_load_si128((const __m128i*)table1+0));
  m6=_mm_madd_epi16(m6, _mm_load_si128((const __m128i*)table2+3));
  m4=_mm_madd_epi16(m4, _mm_load_si128((const __m128i*)table2+0));
  
  m0=_mm_add_epi32(m0, m1);
  m4=_mm_add_epi32(m4, m5);
  m2=_mm_add_epi32(m2, m3);
  m6=_mm_add_epi32(m6, m7);
  
  m0=_mm_srai_epi32(m0, 16);
  m4=_mm_srai_epi32(m4, 16);
  m2=_mm_srai_epi32(m2, 16);
  m6=_mm_srai_epi32(m6, 16);

  m0=_mm_packs_epi32(m0, m2);
  m4=_mm_packs_epi32(m4, m6);
  
  _mm_store_si128((__m128i*)(data+0*8), m0);
  _mm_store_si128((__m128i*)(data+1*8), m4);
#endif
}
static ALIGN(16) const short DCT_SSE2_table[]=//google/webmproject/sjpeg/fdct.cc
{
	// Tables for fdct, roughly the transposed of the above, shuffled
	0x4000, 0x4000, 0x58c5, 0x4b42, 0xdd5d, 0xac61, 0xa73b, 0xcdb7,
	0x4000, 0x4000, 0x3249, 0x11a8, 0x539f, 0x22a3, 0x4b42, 0xee58,
	0x4000, 0xc000, 0x3249, 0xa73b, 0x539f, 0xdd5d, 0x4b42, 0xa73b,
	0xc000, 0x4000, 0x11a8, 0x4b42, 0x22a3, 0xac61, 0x11a8, 0xcdb7,
	0x58c5, 0x58c5, 0x7b21, 0x6862, 0xcff5, 0x8c04, 0x84df, 0xba41,
	0x58c5, 0x58c5, 0x45bf, 0x187e, 0x73fc, 0x300b, 0x6862, 0xe782,
	0x58c5, 0xa73b, 0x45bf, 0x84df, 0x73fc, 0xcff5, 0x6862, 0x84df,
	0xa73b, 0x58c5, 0x187e, 0x6862, 0x300b, 0x8c04, 0x187e, 0xba41,
	0x539f, 0x539f, 0x73fc, 0x6254, 0xd2bf, 0x92bf, 0x8c04, 0xbe4d,
	0x539f, 0x539f, 0x41b3, 0x1712, 0x6d41, 0x2d41, 0x6254, 0xe8ee,
	0x539f, 0xac61, 0x41b3, 0x8c04, 0x6d41, 0xd2bf, 0x6254, 0x8c04,
	0xac61, 0x539f, 0x1712, 0x6254, 0x2d41, 0x92bf, 0x1712, 0xbe4d,
	0x4b42, 0x4b42, 0x6862, 0x587e, 0xd746, 0x9dac, 0x979e, 0xc4df,
	0x4b42, 0x4b42, 0x3b21, 0x14c3, 0x6254, 0x28ba, 0x587e, 0xeb3d,
	0x4b42, 0xb4be, 0x3b21, 0x979e, 0x6254, 0xd746, 0x587e, 0x979e,
	0xb4be, 0x4b42, 0x14c3, 0x587e, 0x28ba, 0x9dac, 0x14c3, 0xc4df,
};
void DCT2_8x8_i16_buf(const short *src, int iw, int ih, short *dst, float quantize, const unsigned char *qmatrix)
{
	//exact DCT
#if 1
	FCT1D_PS_Params p;
	float temp[64];
	FCT1D_ps_gen(3, &p);
	for(int ky=0, kb=0;ky<ih-7;ky+=8)
	{
		for(int kx=0;kx<iw-7;kx+=8, ++kb)
		{
			for(int ky2=0;ky2<8;++ky2)
			{
				for(int kx2=0;kx2<8;++kx2)
					temp[ky2<<3|kx2]=(float)src[iw*(ky+ky2)+kx+kx2];
			}

			for(int k=0;k<8;++k)//row DCT2
				FCT1D_ps_fwd(&p, temp+(k<<3), 1);

			for(int k=0;k<8;++k)//column DCT2
				FCT1D_ps_fwd(&p, temp+k, 8);

			if(quantize)
			{
				if(qmatrix)//quantize
				{
					for(int k=0;k<64;++k)
						temp[k]/=qmatrix[k]*quantize;
				}
				else
				{
					for(int k=0;k<64;++k)
						temp[k]/=quantize;
				}
			}
			
			for(int k=0;k<64;++k)
			{
				float val=roundf(temp[k]);
				if(val<-0x7FFF||val>0x7FFF)
					LOG_ERROR("DCT2 Overflow");
				dst[kb<<6|k]=(short)val;
			}
		}
	}
	FCT1D_ps_free(&p);
#endif

	//fast DCT2 from sjpeg
#if 0
	ALIGN(16) short temp[64];
	for(int ky=0, kb=0;ky<ih-7;ky+=8)
	{
		for(int kx=0;kx<iw-7;kx+=8, ++kb)
		{
			//if(ky==0&&kx==8)//
			//	printf("");

			int idx=iw*ky+kx;
			for(int ky2=0;ky2<8;++ky2)
				memcpy(temp+(ky2<<3), src+idx+iw*ky2, 8*sizeof(short));

			DCT_col_SSE2(temp);
			DCT_row_SSE2(temp+0*8, DCT_SSE2_table+0*32, DCT_SSE2_table+1*32);
			DCT_row_SSE2(temp+2*8, DCT_SSE2_table+2*32, DCT_SSE2_table+3*32);
			DCT_row_SSE2(temp+4*8, DCT_SSE2_table+0*32, DCT_SSE2_table+3*32);
			DCT_row_SSE2(temp+6*8, DCT_SSE2_table+2*32, DCT_SSE2_table+1*32);

			//for(int k=0;k<64;++k)//
			//{
			//	if((unsigned)(temp[k]+(1<<(quantize-1)))>=(unsigned)(1<<quantize))//
			//		LOG_ERROR("Range error");//
			//}
			if(quantize)
			{
				if(qmatrix)//quantize
				{
					for(int k=0;k<64;++k)
						temp[k]/=qmatrix[k];
				}
				for(int k=0;k<64;++k)
					temp[k]>>=7;//15 -> 8 bit
			}

			memcpy(dst+(kb<<6), temp, 64*sizeof(short));
		}
	}
#endif
#if 0
	int vmin=0, vmax=0;
	for(int k=0, res=iw*ih;k<res;++k)
	{
		if(vmin>dst[k])
			vmin=dst[k];
		if(vmax<dst[k])
			vmax=dst[k];
	}
	printf("DCT [%d, %d]\n", vmin, vmax);
#endif
}


#define stbi__f2f(x)  ((int) (((x) * 4096 + 0.5)))
static void stbi__idct_simd(char *out, int out_stride, const short *data)//stb_image.h(2415)
{
   // This is constructed to match our regular (generic) integer IDCT exactly.
   __m128i row0, row1, row2, row3, row4, row5, row6, row7;
   __m128i tmp;

   // dot product constant: even elems=x, odd elems=y
   #define dct_const(x,y)  _mm_setr_epi16((x),(y),(x),(y),(x),(y),(x),(y))

   // out(0) = c0[even]*x + c0[odd]*y   (c0, x, y 16-bit, out 32-bit)
   // out(1) = c1[even]*x + c1[odd]*y
   #define dct_rot(out0,out1, x,y,c0,c1) \
	  __m128i c0##lo = _mm_unpacklo_epi16((x),(y)); \
	  __m128i c0##hi = _mm_unpackhi_epi16((x),(y)); \
	  __m128i out0##_l = _mm_madd_epi16(c0##lo, c0); \
	  __m128i out0##_h = _mm_madd_epi16(c0##hi, c0); \
	  __m128i out1##_l = _mm_madd_epi16(c0##lo, c1); \
	  __m128i out1##_h = _mm_madd_epi16(c0##hi, c1)

   // out = in << 12  (in 16-bit, out 32-bit)
   #define dct_widen(out, in) \
	  __m128i out##_l = _mm_srai_epi32(_mm_unpacklo_epi16(_mm_setzero_si128(), (in)), 4); \
	  __m128i out##_h = _mm_srai_epi32(_mm_unpackhi_epi16(_mm_setzero_si128(), (in)), 4)

   // wide add
   #define dct_wadd(out, a, b) \
	  __m128i out##_l = _mm_add_epi32(a##_l, b##_l); \
	  __m128i out##_h = _mm_add_epi32(a##_h, b##_h)

   // wide sub
   #define dct_wsub(out, a, b) \
	  __m128i out##_l = _mm_sub_epi32(a##_l, b##_l); \
	  __m128i out##_h = _mm_sub_epi32(a##_h, b##_h)

   // butterfly a/b, add bias, then shift by "s" and pack
   #define dct_bfly32o(out0, out1, a,b,bias,s) \
	  { \
		 __m128i abiased_l = _mm_add_epi32(a##_l, bias); \
		 __m128i abiased_h = _mm_add_epi32(a##_h, bias); \
		 dct_wadd(sum, abiased, b); \
		 dct_wsub(dif, abiased, b); \
		 out0 = _mm_packs_epi32(_mm_srai_epi32(sum_l, s), _mm_srai_epi32(sum_h, s)); \
		 out1 = _mm_packs_epi32(_mm_srai_epi32(dif_l, s), _mm_srai_epi32(dif_h, s)); \
	  }

   // 8-bit interleave step (for transposes)
   #define dct_interleave8(a, b) \
	  tmp = a; \
	  a = _mm_unpacklo_epi8(a, b); \
	  b = _mm_unpackhi_epi8(tmp, b)

   // 16-bit interleave step (for transposes)
   #define dct_interleave16(a, b) \
	  tmp = a; \
	  a = _mm_unpacklo_epi16(a, b); \
	  b = _mm_unpackhi_epi16(tmp, b)

   #define dct_pass(bias,shift) \
	  { \
		 /* even part */ \
		 dct_rot(t2e,t3e, row2,row6, rot0_0,rot0_1); \
		 __m128i sum04 = _mm_add_epi16(row0, row4); \
		 __m128i dif04 = _mm_sub_epi16(row0, row4); \
		 dct_widen(t0e, sum04); \
		 dct_widen(t1e, dif04); \
		 dct_wadd(x0, t0e, t3e); \
		 dct_wsub(x3, t0e, t3e); \
		 dct_wadd(x1, t1e, t2e); \
		 dct_wsub(x2, t1e, t2e); \
		 /* odd part */ \
		 dct_rot(y0o,y2o, row7,row3, rot2_0,rot2_1); \
		 dct_rot(y1o,y3o, row5,row1, rot3_0,rot3_1); \
		 __m128i sum17 = _mm_add_epi16(row1, row7); \
		 __m128i sum35 = _mm_add_epi16(row3, row5); \
		 dct_rot(y4o,y5o, sum17,sum35, rot1_0,rot1_1); \
		 dct_wadd(x4, y0o, y4o); \
		 dct_wadd(x5, y1o, y5o); \
		 dct_wadd(x6, y2o, y5o); \
		 dct_wadd(x7, y3o, y4o); \
		 dct_bfly32o(row0,row7, x0,x7,bias,shift); \
		 dct_bfly32o(row1,row6, x1,x6,bias,shift); \
		 dct_bfly32o(row2,row5, x2,x5,bias,shift); \
		 dct_bfly32o(row3,row4, x3,x4,bias,shift); \
	  }

   __m128i rot0_0 = dct_const(stbi__f2f(0.5411961f), stbi__f2f(0.5411961f) + stbi__f2f(-1.847759065f));
   __m128i rot0_1 = dct_const(stbi__f2f(0.5411961f) + stbi__f2f( 0.765366865f), stbi__f2f(0.5411961f));
   __m128i rot1_0 = dct_const(stbi__f2f(1.175875602f) + stbi__f2f(-0.899976223f), stbi__f2f(1.175875602f));
   __m128i rot1_1 = dct_const(stbi__f2f(1.175875602f), stbi__f2f(1.175875602f) + stbi__f2f(-2.562915447f));
   __m128i rot2_0 = dct_const(stbi__f2f(-1.961570560f) + stbi__f2f( 0.298631336f), stbi__f2f(-1.961570560f));
   __m128i rot2_1 = dct_const(stbi__f2f(-1.961570560f), stbi__f2f(-1.961570560f) + stbi__f2f( 3.072711026f));
   __m128i rot3_0 = dct_const(stbi__f2f(-0.390180644f) + stbi__f2f( 2.053119869f), stbi__f2f(-0.390180644f));
   __m128i rot3_1 = dct_const(stbi__f2f(-0.390180644f), stbi__f2f(-0.390180644f) + stbi__f2f( 1.501321110f));

   // rounding biases in column/row passes, see stbi__idct_block for explanation.
   __m128i bias_0 = _mm_set1_epi32(512);
   __m128i bias_1 = _mm_set1_epi32(65536 + (128<<17));

   // load
   row0 = _mm_load_si128((const __m128i *) (data + 0*8));
   row1 = _mm_load_si128((const __m128i *) (data + 1*8));
   row2 = _mm_load_si128((const __m128i *) (data + 2*8));
   row3 = _mm_load_si128((const __m128i *) (data + 3*8));
   row4 = _mm_load_si128((const __m128i *) (data + 4*8));
   row5 = _mm_load_si128((const __m128i *) (data + 5*8));
   row6 = _mm_load_si128((const __m128i *) (data + 6*8));
   row7 = _mm_load_si128((const __m128i *) (data + 7*8));

   // column pass
   dct_pass(bias_0, 10);

   {
	  // 16bit 8x8 transpose pass 1
	  dct_interleave16(row0, row4);
	  dct_interleave16(row1, row5);
	  dct_interleave16(row2, row6);
	  dct_interleave16(row3, row7);

	  // transpose pass 2
	  dct_interleave16(row0, row2);
	  dct_interleave16(row1, row3);
	  dct_interleave16(row4, row6);
	  dct_interleave16(row5, row7);

	  // transpose pass 3
	  dct_interleave16(row0, row1);
	  dct_interleave16(row2, row3);
	  dct_interleave16(row4, row5);
	  dct_interleave16(row6, row7);
   }

   // row pass
   dct_pass(bias_1, 17);

   {
	  // pack
	  __m128i p0 = _mm_packus_epi16(row0, row1); // a0a1a2a3...a7b0b1b2b3...b7
	  __m128i p1 = _mm_packus_epi16(row2, row3);
	  __m128i p2 = _mm_packus_epi16(row4, row5);
	  __m128i p3 = _mm_packus_epi16(row6, row7);

	  // 8bit 8x8 transpose pass 1
	  dct_interleave8(p0, p2); // a0e0a1e1...
	  dct_interleave8(p1, p3); // c0g0c1g1...

	  // transpose pass 2
	  dct_interleave8(p0, p1); // a0c0e0g0...
	  dct_interleave8(p2, p3); // b0d0f0h0...

	  // transpose pass 3
	  dct_interleave8(p0, p2); // a0b0c0d0...
	  dct_interleave8(p1, p3); // a4b4c4d4...

	  // store
	  _mm_storel_epi64((__m128i *) out, p0); out += out_stride;
	  _mm_storel_epi64((__m128i *) out, _mm_shuffle_epi32(p0, 0x4e)); out += out_stride;
	  _mm_storel_epi64((__m128i *) out, p2); out += out_stride;
	  _mm_storel_epi64((__m128i *) out, _mm_shuffle_epi32(p2, 0x4e)); out += out_stride;
	  _mm_storel_epi64((__m128i *) out, p1); out += out_stride;
	  _mm_storel_epi64((__m128i *) out, _mm_shuffle_epi32(p1, 0x4e)); out += out_stride;
	  _mm_storel_epi64((__m128i *) out, p3); out += out_stride;
	  _mm_storel_epi64((__m128i *) out, _mm_shuffle_epi32(p3, 0x4e));
   }

#undef dct_const
#undef dct_rot
#undef dct_widen
#undef dct_wadd
#undef dct_wsub
#undef dct_bfly32o
#undef dct_interleave8
#undef dct_interleave16
#undef dct_pass
}
void DCT3_8x8_i16_buf(const short *src, int iw, int ih, short *dst, float dequantize, const unsigned char *qmatrix)
{
	//exact DCT
#if 1
	FCT1D_PS_Params p;
	float temp[64];
	FCT1D_ps_gen(3, &p);
	for(int ky=0, kb=0;ky<ih-7;ky+=8)
	{
		for(int kx=0;kx<iw-7;kx+=8, ++kb)
		{
			for(int k=0;k<64;++k)
				temp[k]=(float)src[kb<<6|k];
			
			if(dequantize)
			{
				if(qmatrix)//quantize
				{
					for(int k=0;k<64;++k)
						temp[k]*=qmatrix[k]*dequantize;
				}
				else
				{
					for(int k=0;k<64;++k)
						temp[k]*=dequantize;
				}
			}

			for(int k=0;k<8;++k)//row DCT3
				FCT1D_ps_inv(&p, temp+(k<<3), 1);

			for(int k=0;k<8;++k)//column DCT3
				FCT1D_ps_inv(&p, temp+k, 8);

			for(int ky2=0;ky2<8;++ky2)
			{
				for(int kx2=0;kx2<8;++kx2)
				{
					float val=roundf(temp[ky2<<3|kx2]);
					if(val<-0x7FFF||val>0x7FFF)
						LOG_ERROR("DCT3 Overflow");
					dst[iw*(ky+ky2)+kx+kx2]=(short)val;
				}
			}
		}
	}
	FCT1D_ps_free(&p);
#endif

	//fast DCT3 from stb_image.h
#if 0
	ALIGN(16) short temp[64];
	for(int ky=0, kb=0;ky<ih-7;ky+=8)
	{
		for(int kx=0;kx<iw-7;kx+=8, ++kb)
		{
			memcpy(temp, dst+(kb<<6), 64*sizeof(short));

			if(dequantize)
			{
				if(qmatrix)//quantize
				{
					for(int k=0;k<64;++k)
						temp[k]*=qmatrix[k];
				}
				//for(int k=0;k<64;++k)
				//	temp[k]>>=7;//15 -> 8 bit
			}

			stbi__idct_simd(dst+iw*ky+kx, iw, temp);

			//int idx=iw*ky+kx;
			//for(int ky2=0;ky2<8;++ky2)
			//	memcpy(temp+(ky2<<3), src+idx+iw*ky2, 8*sizeof(short));
		}
	}
#endif
#if 0
	int vmin=0, vmax=0;
	for(int k=0, res=iw*ih;k<res;++k)
	{
		if(vmin>dst[k])
			vmin=dst[k];
		if(vmax<dst[k])
			vmax=dst[k];
	}
	printf("DCT [%d, %d]\n", vmin, vmax);
#endif
}


void quantize(float *buf, int iw, int ih, const float *qmatrix)
{
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=iw*ky+kx;
			float val=buf[idx];
			val/=qmatrix[(ky&7)<<3|kx&7];
			val=roundf(val);
			buf[idx]=CLAMP(-128, val, 128);
		}
	}
}




//DWTs

//	#define CHECK_OVERFLOW

#ifdef CHECK_OVERFLOW
int FIX_MUL(int a, int b, int sr)
{
	long long x=(long long)a*b>>sr;
	if(x<-0x8000||x>0x7FFF)
		LOG_ERROR("Overflow");
	return (int)x;
}
int FIX_DIV(int num, int den, int sl)
{
	long long x=((long long)num<<16)/den;
	if(x<-0x8000||x>0x7FFF)
		LOG_ERROR("Overflow");
	return (int)x;
}
#else
#define FIX_MUL(A, B, SR) (int)((long long)(A)*(B)>>SR)
#define FIX_DIV(N, D, SL) (int)(((long long)(N)<<SL)/(D))
#endif
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

//lifting-based 16bit CDF 9/7 DWT
static const int cdf97_coeff[]=
{
	//'factring wavelet transforms into lifting steps' - page 19
	//'a wavelet tour of signal processing - the sparse way' - page 376
	//Lifting-based Discrete Wavelet Transform for Real-Time Signal Detection 2015-10
	 0x1960C,	//-1.58613434342059f,	//alpha
	 0x00D90,	//-0.0529801185729f,	//beta
	-0x0E206,	// 0.8829110755309f,	//gamma
	-0x07189,	// 0.4435068520439f,	//delta
	 0x1264C,	// 1.1496043988602f,	//zeta		output gain is 1.89
};
static void dwt1d_i16_predict(short *odd, short *even, int nodd, int extraeven, int coeff)
{
	even[0]-=FIX_MUL(odd[0], coeff, 15);
	for(int k=1;k<nodd;++k)//predict
		even[k]-=FIX_MUL(odd[k-1]+odd[k], coeff, 16);
	if(extraeven)
		even[nodd]-=FIX_MUL(odd[nodd-1], coeff, 15);
}
static void dwt1d_i16_update(short *odd, short *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=FIX_MUL(even[k]+even[k+1], coeff, 16);
	if(!extraeven)
		odd[nodd-1]+=FIX_MUL(even[nodd-1], coeff, 15);
}
static void dwt1d_i16_unpredict(short *odd, short *even, int nodd, int extraeven, int coeff)
{
	even[0]+=FIX_MUL(odd[0], coeff, 15);
	for(int k=1;k<nodd;++k)//unpredict
		even[k]+=FIX_MUL(odd[k-1]+odd[k], coeff, 16);
	if(extraeven)
		even[nodd]+=FIX_MUL(odd[nodd-1], coeff, 15);
}
static void dwt1d_i16_unupdate(short *odd, short *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//unupdate
		odd[k]-=FIX_MUL(even[k]+even[k+1], coeff, 16);
	if(!extraeven)
		odd[nodd-1]-=FIX_MUL(even[nodd-1], coeff, 15);
}
static void dwt1d_i16_mul(short *buf, int count, int coeff)//scale
{
	for(int k=0;k<count;++k)
		buf[k]=FIX_MUL(buf[k], coeff, 16);
}
static void dwt1d_i16_div(short *buf, int count, int coeff)//unscale
{
	for(int k=0;k<count;++k)
		buf[k]=FIX_DIV(buf[k], coeff, 16);
}
void dwt1d_cdf97_fwd(short *buffer, int count, int stride, short *b2, int qfactor)
{
	int nodd=count>>1, extraeven=count&1;
	short *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

	dwt1d_i16_predict(odd, even, nodd, extraeven, cdf97_coeff[0]);
	dwt1d_i16_update (odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_i16_predict(odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_i16_update (odd, even, nodd, extraeven, cdf97_coeff[3]);
	//dwt1d_i16_div(even, nodd+extraeven, qfactor);
	//dwt1d_i16_mul(odd, nodd, cdf97_coeff[4]);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_cdf97_inv(short *buffer, int count, int stride, short *b2, int qfactor)
{
	int nodd=count>>1, extraeven=count&1;
	short *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	//dwt1d_i16_mul(even, nodd+extraeven, qfactor);
	//dwt1d_i16_div(odd, nodd, cdf97_coeff[4]);
	dwt1d_i16_unupdate (odd, even, nodd, extraeven, cdf97_coeff[3]);
	dwt1d_i16_unpredict(odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_i16_unupdate (odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_i16_unpredict(odd, even, nodd, extraeven, cdf97_coeff[0]);

	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf97_fwd(short *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, short *temp, int qfactor)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf97_fwd(buffer+rowlen*ky, w2, stride, temp, qfactor);

		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dA.PNG", it);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf97_fwd(buffer+stride*kx, h2, rowlen, temp, qfactor);
		
		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_cdf97_inv(short *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, short *temp, int qfactor)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_cdf97_inv(buffer+stride*kx, h2, rowlen, temp, qfactor);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_cdf97_inv(buffer+rowlen*ky, w2, stride, temp, qfactor);
	}
}


static const float cdf97ps_coeff[]=
{
	//'factring wavelet transforms into lifting steps' - page 19
	//'a wavelet tour of signal processing - the sparse way' - page 376
	//Lifting-based Discrete Wavelet Transform for Real-Time Signal Detection 2015-10
	+1.58613434342059f,	//alpha			sign is flipped to make predict subtract
	-0.0529801185729f,	//beta
	-0.8829110755309f,	//gamma			sign is flipped to make predict subtract
	 0.4435068520439f,	//delta
	 //1.1496043988602f,	//zeta		unused		//output gain is 1.89 ?

	 1.230174105f,//K	openjpeg/src/lib/openjp2/dwt.c(113)
};
static void dwtps1d_predict(float *odd, float *even, int nodd, int extraeven, float coeff)
{
	even[0]-=odd[0]*2*coeff;
	for(int k=1;k<nodd;++k)//predict				//nonstandard: in JPEG2000's CDF 9/7, predict modifies odd samples by adding
		even[k]-=(odd[k-1]+odd[k])*coeff;
	if(extraeven)
		even[nodd]-=odd[nodd-1]*2*coeff;
}
static void dwtps1d_update(float *odd, float *even, int nodd, int extraeven, float coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//update		//nonstandard: in JPEG2000's CDF 9/7, update modifies even samples
		odd[k]+=(even[k]+even[k+1])*coeff;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]*2*coeff;
}
static void dwtps1d_unpredict(float *odd, float *even, int nodd, int extraeven, float coeff)
{
	even[0]+=odd[0]*2*coeff;
	for(int k=1;k<nodd;++k)//unpredict
		even[k]+=(odd[k-1]+odd[k])*coeff;
	if(extraeven)
		even[nodd]+=odd[nodd-1]*2*coeff;
}
static void dwtps1d_unupdate(float *odd, float *even, int nodd, int extraeven, float coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//unupdate
		odd[k]-=(even[k]+even[k+1])*coeff;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]*2*coeff;
}
static void dwtps1d_mul(float *buf, int count, float coeff)//scale
{
	for(int k=0;k<count;++k)
		buf[k]*=coeff;
}
static void dwtps1d_div(float *buf, int count, float coeff)//unscale
{
	float gain=1/coeff;
	for(int k=0;k<count;++k)
		buf[k]*=gain;
}
void cdf97ps1d_fwd(float *buffer, int count, int stride, float *b2, float qfactor)
{
	int nodd=count>>1, extraeven=count&1;
	float *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

	dwtps1d_predict(odd, even, nodd, extraeven, cdf97ps_coeff[0]);
	dwtps1d_update (odd, even, nodd, extraeven, cdf97ps_coeff[1]);
	dwtps1d_predict(odd, even, nodd, extraeven, cdf97ps_coeff[2]);
	dwtps1d_update (odd, even, nodd, extraeven, cdf97ps_coeff[3]);
	dwtps1d_div(even, nodd+extraeven, cdf97ps_coeff[4]);
	dwtps1d_mul(odd, nodd, cdf97ps_coeff[4]);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void cdf97ps1d_inv(float *buffer, int count, int stride, float *b2, float qfactor)
{
	int nodd=count>>1, extraeven=count&1;
	float *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	dwtps1d_div(odd, nodd, cdf97ps_coeff[4]);
	dwtps1d_mul(even, nodd+extraeven, cdf97ps_coeff[4]);
	dwtps1d_unupdate (odd, even, nodd, extraeven, cdf97ps_coeff[3]);
	dwtps1d_unpredict(odd, even, nodd, extraeven, cdf97ps_coeff[2]);
	dwtps1d_unupdate (odd, even, nodd, extraeven, cdf97ps_coeff[1]);
	dwtps1d_unpredict(odd, even, nodd, extraeven, cdf97ps_coeff[0]);

	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void cdf97ps2d_fwd(float *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, float *temp, float qfactor)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			cdf97ps1d_fwd(buffer+rowlen*ky, w2, stride, temp, qfactor);

		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dA.PNG", it);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			cdf97ps1d_fwd(buffer+stride*kx, h2, rowlen, temp, qfactor);
		
		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dB.PNG", it);
	}
}
void cdf97ps2d_inv(float *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, float *temp, float qfactor)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			cdf97ps1d_inv(buffer+stride*kx, h2, rowlen, temp, qfactor);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			cdf97ps1d_inv(buffer+rowlen*ky, w2, stride, temp, qfactor);
	}
}

void lg53ps1d_fwd(float *buffer, int count, int stride, float *b2, float qfactor)
{
	int nodd=count>>1, extraeven=count&1;
	float *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

	dwtps1d_predict(odd, even, nodd, extraeven, 0.5f);
	dwtps1d_update (odd, even, nodd, extraeven, 0.25f);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void lg53ps1d_inv(float *buffer, int count, int stride, float *b2, float qfactor)
{
	int nodd=count>>1, extraeven=count&1;
	float *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	dwtps1d_unupdate (odd, even, nodd, extraeven, 0.25f);
	dwtps1d_unpredict(odd, even, nodd, extraeven, 0.5f);

	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void lg53ps2d_fwd(float *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, float *temp, float qfactor)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			lg53ps1d_fwd(buffer+rowlen*ky, w2, stride, temp, qfactor);

		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dA.PNG", it);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			lg53ps1d_fwd(buffer+stride*kx, h2, rowlen, temp, qfactor);
		
		//save_channel(buffer, iw, ih, 4, 128, "cdf53-stage%02dB.PNG", it);
	}
}
void lg53ps2d_inv(float *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, float *temp, float qfactor)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, ih=sizes->h, rowlen=stride*iw;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			lg53ps1d_inv(buffer+stride*kx, h2, rowlen, temp, qfactor);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			lg53ps1d_inv(buffer+rowlen*ky, w2, stride, temp, qfactor);
	}
}


void quantize_dwt_i16_fwd(short *buf, int llw, int llh, int iw, int ih, int Q)
{
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			if(ky<llh&&kx<llh)
				continue;
			int idx=iw*ky+kx;
			buf[idx]=FIX_DIV(buf[idx], Q, 16);
		}
	}
}
void quantize_dwt_i16_inv(short *buf, int llw, int llh, int iw, int ih, int Q)
{
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			if(ky<llh&&kx<llh)
				continue;
			int idx=iw*ky+kx;
			buf[idx]=FIX_MUL(buf[idx], Q, 16);
		}
	}
}

void quantize_dwt_ps_i16(const float *src, short *dst, int llw, int llh, int iw, int ih, float Q, float Z)//dead zone Z is at least 1
{
	if(Z<1)
		Z=1;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			//if(ky<llh&&kx<llw)
			//	continue;

			//if(kx==(iw*2/3)&&ky==(ih*2/3))
			//	printf("");

			int idx=iw*ky+kx;
			float val=src[idx];
			int neg=val<0;
			val=fabsf(val)-(Z-1);
			if(val<0)
				val=0;
			val/=Q;
			if(neg)
				val=-val;
			if(val<-0x8000||val>0x7FFF)
				LOG_ERROR("Overflow");
			val=CLAMP(-0x8000, val, 0x7FFF);
			dst[idx]=(short)val;
		}
	}
}
void dequantize_dwt_i16_ps(const short *src, float *dst, int llw, int llh, int iw, int ih, float Q, float Z)
{
	if(Z<1)
		Z=1;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			//if(ky<llh&&kx<llw)
			//	continue;

			//if(kx==(iw*2/3)&&ky==(ih*2/3))
			//	printf("");

			int idx=iw*ky+kx;
			float val=src[idx];
			int neg=val<0;
			val*=Q;
			val=fabsf(val);
			val+=Z-1;
			if(neg)
				val=-val;
			dst[idx]=val;
		}
	}
}