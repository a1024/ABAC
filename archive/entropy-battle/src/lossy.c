#include"battle.h"
#include<stdio.h>//for debugging
#include<stdlib.h>
#include<string.h>
#ifdef __GNUC__
#include<x86intrin.h>
#elif defined _MSC_VER
#include<intrin.h>
#endif
#include<math.h>
#include"lodepng.h"//for testing
static const char file[]=__FILE__;

void calc_distortion(unsigned char *b0, unsigned char *b1, int iw, int ih, int symbytes, int bytestride, RateDistortion *ret)//returned info is for each channel separately
{
	ptrdiff_t res=(ptrdiff_t)iw*ih, len=bytestride*res;
	if(!ret)
		return;
	RateDistortion *p;
	for(int kc=0;kc<symbytes;++kc)
	{
		long long mean=0;
		for(ptrdiff_t k=0;k<len;k+=bytestride)
		{
			unsigned char diff=b0[k+kc]-b1[k+kc];
			mean+=diff*diff;
		}
		double RMSE=sqrt((double)mean/res), PSNR=20*log10(255/RMSE);
		if(ret)
		{
			p=ret+kc;
			p->RMSE=RMSE;
			p->PSNR=PSNR;
		}
	}
}

int lossy1_encode(const void *src, int iw, int ih, int symbytes, int bytestride, ArrayHandle *data, int bitrate)
{
	if(bytestride!=4||symbytes!=3)
	{
		LOG_ERROR("Expected 8 bit 3 channels with stride 4 bytes");
		return 0;
	}
	const unsigned char *srcbuf=(const unsigned char*)src;
	ptrdiff_t res=(ptrdiff_t)iw*ih;

	float *tbuf=(float*)malloc(res*3*sizeof(float));
	if(!tbuf)
	{
		LOG_ERROR("Allocation error");
		return 0;
	}
	float *channels[3]={tbuf, tbuf+res, tbuf+(res<<1)};
	YCoCg_8bit_ps_fwd(srcbuf, res, channels[0], channels[1], channels[2]);

	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 8, 8, 0);
	if(!sizes)
	{
		LOG_ERROR("Allocation error");
		free(tbuf);
		return 0;
	}
	float vmin[3]={0}, vmax[3]={0};
	//float Ymean=0.5f;
	for(int kc=0;kc<3;++kc)
	{
		float *buf=channels[kc];

		if(!kc)//scale Y to [-1, 1]
		{
			for(ptrdiff_t kp=0;kp<res;++kp)
			{
				float val=buf[kp];
				val-=0.5;
				val+=val;
				buf[kp]=val;
			}
#if 0
			double sum=0;
			for(int kp=0;kp<res;++kp)//X  apply per block and differentiate
				sum+=buf[kp];
			Ymean=sum/res;
			for(int kp=0;kp<res;++kp)
				buf[kp]-=Ymean;
#endif
		}

		dwt2d_ps_fwd(buf, (DWTSize*)sizes->data, (int)sizes->count);

		//normalize buffer
		float *pmin=vmin+kc, *pmax=vmax+kc;
		for(ptrdiff_t kp=0;kp<res;++kp)
		{
			float val=buf[kp];
			if(*pmin>val)
				*pmin=val;
			if(*pmax<val)
				*pmax=val;
		}
		float gain=*pmax-*pmin;
		if(gain>0)//X  normalization should only depend on DWT type
		{
			gain=1/gain;
			for(ptrdiff_t kp=0;kp<res;++kp)
			{
				float val=buf[kp];
				val-=*pmin;
				val*=gain;
				buf[kp]=val;
			}
		}

		ptrdiff_t ks=sizes->count-1;
		for(;ks>0;--ks)
		{

		}
	}

	array_free(&sizes);
	free(tbuf);
	return 1;
}
int    lossy1_decode(const void *src, ptrdiff_t srclen, int iw, int ih, int symbytes, int bytestride, void *dst, int bitrate)
{
	if(bytestride!=4||symbytes!=3)
	{
		LOG_ERROR("Expected 8 bit 3 channels with stride 4 bytes");
		return 0;
	}
	return 1;
}