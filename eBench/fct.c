#include"ebench.h"
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#if defined _MSC_VER
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif

//FCT-II/III, O(N*lg(N)), POT size (uses FFT size N)
static void fft_dit_ps(float *re, float *im, int lgsize, int forward)//DIT vs DIF: simply reverse the order of stage matrices along with the permutation, and transpose the stages
{
	int size=1<<lgsize, sign=!forward-(forward!=0);
	int k, m, j;
	
	for(k=0;k<size;k+=2)//separate mul-free sum & difference (m==2)
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
	for(m=4;m<=size;m<<=1)
	{
		int m2=m>>1;
		float angle=(float)((2*M_PI)*sign/m);
		float rewm=(float)cosf(angle), imwm=(float)sinf(angle);//exp(-+ 2*pi*i*(N/m)/N)
		for(k=0;k<size;k+=m)
		{

			float
				re0=re[k],
				im0=im[k],
				re1=re[k+m2],
				im1=im[k+m2];
			float rew=rewm, imw=imwm;

			re[k]=re0+re1;
			im[k]=im0+im1;
			re[k+m2]=re0-re1;
			im[k+m2]=im0-im1;
			for(j=1;;)
			{
				float re2, im2;
				int idx=k+j;
				re0=re[idx];
				im0=im[idx];
				re1=re[idx+m2];
				im1=im[idx+m2];
				
				re2=rew*re1-imw*im1;
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
static void fft_dif_ps(float *re, float *im, int lgsize, int forward)
{
	int size=1<<lgsize, sign=!forward-(forward!=0);
	int k, m, j;
	
	for(m=size;m>2;m>>=1)
	{
		int m2=m>>1;
		float angle=(float)((2*M_PI)*sign/m);
		float rewm=(float)cosf(angle), imwm=(float)sinf(angle);//exp(-+ 2*pi*i*(N/m)/N)
		for(k=0;k<size;k+=m)
		{
			float
				re0=re[k],
				im0=im[k],
				re1=re[k+m2],
				im1=im[k+m2];
			float rew=rewm, imw=imwm;

			re[k]=re0+re1;
			im[k]=im0+im1;
			re[k+m2]=re0-re1;
			im[k+m2]=im0-im1;
			
			for(j=1;;)
			{
				int idx=k+j;

				re0=re[idx];
				im0=im[idx];
				re1=re[idx+m2];
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
	for(k=0;k<size;k+=2)//separate mul-free sum & difference (m==2)
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
#if 1
typedef struct _FCT1D_PS_Params
{
	int lgsize;
	int *fctp;
	float *re, *im, *rew4N, *imw4N;
} FCT1D_PS_Params;
#endif
#define BIT_REVERSE(N, NBITS)\
	N=(N<<1&0xAAAA)|(N>>1&0x5555),\
	N=(N<<2&0xCCCC)|(N>>2&0x3333),\
	N=(N<<4&0xF0F0)|(N>>4&0x0F0F),\
	N=(N<<8&0xFF00)|(N>>8&0x00FF),\
	N>>=16-NBITS
static void FCT1D_ps_gen(int lgsize, FCT1D_PS_Params *p)//up to size 65536
{
	int size=1<<lgsize, mask=(1<<lgsize)-1;
	p->lgsize=lgsize;
	p->fctp=(int*)malloc(size*sizeof(int));
	p->re=(float*)malloc(size*sizeof(float));
	p->im=(float*)malloc(size*sizeof(float));
	p->rew4N=(float*)malloc(size*sizeof(float));
	p->imw4N=(float*)malloc(size*sizeof(float));
	if(!p->fctp||!p->re||!p->im||!p->rew4N||!p->imw4N)
	{
		CRASH("Alloc error");
		return;
	}
	for(int k=0;k<size;++k)//FCT permutation: f(n) = (bit_reverse(n>>1)^-(n&1)) & ((1<<lgsize)-1)
	{
		uint32_t val=k>>1;
		BIT_REVERSE(val, lgsize);
		val^=-(k&1);
		val&=mask;
		p->fctp[k]=val;
	}
	{
		float invsqrtN=(float)(1/sqrt(size));
		for(int k=0;k<size;++k)
		{
			float temp=(float)(-2*M_PI*k/((size_t)size<<2));
			p->rew4N[k]=invsqrtN*cosf(temp);
			p->imw4N[k]=invsqrtN*sinf(temp);
		}
	}
}
static void FCT1D_ps_free(FCT1D_PS_Params *p)
{
	free(p->fctp);
	free(p->re);
	free(p->im);
	free(p->rew4N);
	free(p->imw4N);
	memset(p, 0, sizeof(*p));
}
static void FCT1D_ps_fwd(FCT1D_PS_Params *p, float *data, int stride)//data, re, im, rew, imw are all of size (1<<lgsize), result is in data buffer
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
static void FCT1D_ps_inv(FCT1D_PS_Params *p, float *data, int stride)//data, re, im, rew, imw are all of size (1<<lgsize), result is in data buffer
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

	for(int ks=0, kd=0;kd<size;++ks, ++kd)
		data[p->fctp[kd]*stride]=p->re[ks];
}

void pred_fct(Image *src, int fwd)
{
	enum
	{
		INTGAIN=16
	};
	size_t bsize;
	float *buf;
	FCT1D_PS_Params px={0}, py={0};
	int lx, ly, kc, kx, ky;

	lx=31^LZCNT32(src->iw);
	ly=31^LZCNT32(src->ih);
	bsize=sizeof(float)<<(lx+ly);
	buf=(float*)malloc(bsize);
	if(!buf)
	{
		CRASH("Alloc error");
		return;
	}
	FCT1D_ps_gen(lx, &px);
	FCT1D_ps_gen(ly, &py);
//	g=1.f;
//	g=1.f/sqrtf((float)(src->iw*src->ih));
	if(fwd)
	{
	//	g*=INTGAIN;
		for(kc=0;kc<3;++kc)
		{
			float *p;

			for(p=buf, ky=0;ky<1<<ly;++ky)//load
			{
				int32_t *s=((int32_t*)(src+1))+4*src->iw*ky+kc;
				for(kx=0;kx<1<<lx;++kx)
				{
					*p++=(float)*s;
					s+=4;
				}
			}
			for(ky=0;ky<1<<ly;++ky)
				FCT1D_ps_fwd(&px, buf+((intptr_t)ky<<lx), 1);
			for(kx=0;kx<1<<lx;++kx)
				FCT1D_ps_fwd(&py, buf+kx, 1<<lx);
		//	for(kx=0;kx<1<<(lx+ly);++kx)
		//		buf[kx]*=g;
			for(p=buf, ky=0;ky<1<<ly;++ky)//store
			{
				int32_t *s=((int32_t*)(src+1))+4*src->iw*ky+kc;
				for(kx=0;kx<1<<lx;++kx)
				{
					float x=*p++;
				//	int neg=x<0;
				//	x=logf(fabsf(x));
				//	if(neg)
				//		x=-x;
				//	x*=INTGAIN;
					*s=(int32_t)x;
					s+=4;
				}
			}
		//	src->depth[kc]+=4*(src->depth[kc]<=31-4);
		}
	}
	else
	{
	//	g/=INTGAIN;
		for(kc=0;kc<3;++kc)
		{
			float *p;

			for(p=buf, ky=0;ky<1<<ly;++ky)//load
			{
				int32_t *s=((int32_t*)(src+1))+4*src->iw*ky+kc;
				for(kx=0;kx<1<<lx;++kx)
				{
					float x=(float)*s;
				//	int neg=x<0;
				//	x*=1.f/INTGAIN;
				//	x=expf(fabsf(x));
				//	if(neg)
				//		x=-x;
					*p++=x;
					s+=4;
				}
			}
		//	for(kx=0;kx<1<<(lx+ly);++kx)
		//		buf[kx]*=g;
			for(kx=0;kx<1<<lx;++kx)
				FCT1D_ps_inv(&py, buf+kx, 1<<lx);
			for(ky=0;ky<1<<ly;++ky)
				FCT1D_ps_inv(&px, buf+((intptr_t)ky<<lx), 1);
			for(p=buf, ky=0;ky<1<<ly;++ky)//store
			{
				int32_t *s=((int32_t*)(src+1))+4*src->iw*ky+kc;
				for(kx=0;kx<1<<lx;++kx)
				{
					*s=(int32_t)*p++;
					s+=4;
				}
			}
		//	src->depth[kc]-=4*(src->depth[kc]>=src->src_depth[kc]+4);
		}
	}
	FCT1D_ps_free(&px);
	FCT1D_ps_free(&py);
	free(buf);
}
