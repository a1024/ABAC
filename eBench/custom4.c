#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
//#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>
static const char file[]=__FILE__;

int lossyconv_page=0;//[0~15]: 4 channels max * 4 stages
char lossyconv_params[5*5*4*4]={0};//(r2 = 5*5) * 4 channels max * 4 stages		precision: sign_bit.2.4.clamp_bit
unsigned char lossyconv_stride[2*4]={0}, lossyconv_offset[2*4]={0};//2 dimensions * 4 stages
#define LOAD(BUF, X, Y) ((unsigned)(ky+(Y))<(unsigned)src->ih&&(unsigned)(kx+(X))<(unsigned)src->iw?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
void pred_lossyconv(Image *src)
{
	Image *dst=0;
	image_copy(&dst, src);
	if(!dst->data)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	int *p1=src->data, *p2=dst->data;
	for(int kb=0;kb<4;++kb)//no need for memcpy due to even number of stages:  src->dst->src
	{
		int
			xstart=lossyconv_offset[kb<<1|0], xstride=lossyconv_stride[kb<<1|0]+1,
			ystart=lossyconv_offset[kb<<1|1], ystride=lossyconv_stride[kb<<1|1]+1;
		for(int kc=0;kc<4;++kc)
		{
			int depth=src->depth[kc], half=1<<depth>>1;
			if(!depth)
				continue;
			char *curr_filt=lossyconv_params+5*5*(4*kb+kc);
			int clampsize=0;
			for(int k=0;k<5*5;++k)
				clampsize+=curr_filt[k]&1;
			for(int ky=ystart, idx=0;ky<src->ih;ky+=ystride)
			{
				for(int kx=xstart;kx<src->iw;kx+=xstride)
				{
					int idx=src->iw*ky+kx;
					int nb[]=
					{
						LOAD(p1, -2, -2),
						LOAD(p1, -1, -2),
						LOAD(p1,  0, -2),
						LOAD(p1,  1, -2),
						LOAD(p1,  2, -2),
						LOAD(p1, -2, -1),
						LOAD(p1, -1, -1),
						LOAD(p1,  0, -1),
						LOAD(p1,  1, -1),
						LOAD(p1,  2, -1),
						LOAD(p1, -2,  0),
						LOAD(p1, -1,  0),
						LOAD(p1,  0,  0),
						LOAD(p1,  1,  0),
						LOAD(p1,  2,  0),
						LOAD(p1, -2,  1),
						LOAD(p1, -1,  1),
						LOAD(p1,  0,  1),
						LOAD(p1,  1,  1),
						LOAD(p1,  2,  1),
						LOAD(p1, -2,  2),
						LOAD(p1, -1,  2),
						LOAD(p1,  0,  2),
						LOAD(p1,  1,  2),
						LOAD(p1,  2,  2),
					};
					int pred=0;
					int vmin=-half, vmax=half-1, uninit=0;
					for(int k=0;k<5*5;++k)
					{
						pred+=curr_filt[k]*(nb[k]>>1);
						if(curr_filt[k]&1)
						{
							if(uninit)
								vmin=nb[k], vmax=nb[k];
							else
							{
								UPDATE_MIN(vmin, nb[k]);
								UPDATE_MAX(vmax, nb[k]);
							}
						}
					}
					pred+=8;
					pred>>=4;
					pred=CLAMP(vmin, pred, vmax);

					int curr=p2[idx<<2|kc];
					curr-=pred;
					curr<<=32-depth;//signed-MA
					curr>>=32-depth;
					p2[idx<<2|kc]=curr;
				}
			}
		}
		int *temp;
		SWAPVAR(p1, p2, temp);
	}
	free(dst);
}