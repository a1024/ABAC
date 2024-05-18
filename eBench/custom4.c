#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
//#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include<immintrin.h>
//static const char file[]=__FILE__;

short lossyconv_clipboard=0;
int lossyconv_page=0;//[0~15]: 4 channels max * 4 layers	layer<<2|ch
short lossyconv_params[5*5*4*4]={0};//(r2 = 5*5) * 4 channels max * 4 layers			9-bit precision: sign_bit.2.5.clamp_bit
unsigned char lossyconv_stride[2*4]={0}, lossyconv_offset[2*4]={0};//2 dimensions * 4 layers
unsigned char lossyconv_causalRCT[4]={0};
#define LOAD(BUF, X, Y) ((unsigned)(ky+(Y))<(unsigned)src->ih&&(unsigned)(kx+(X))<(unsigned)src->iw?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
void pred_lossyconv(Image *src)
{
	Image *dst=0;
	image_copy(&dst, src);
	//if(!dst->data)
	//{
	//	LOG_ERROR("Alloc error");
	//	return;
	//}
	for(int kb=0;kb<4;++kb)//no need for memcpy due to even number of stages:  src->dst->src
	{
		int *p1=src->data, *p2=dst->data;
		int
			xstart=lossyconv_offset[kb<<1|0], xstride=lossyconv_stride[kb<<1|0]+1,
			ystart=lossyconv_offset[kb<<1|1], ystride=lossyconv_stride[kb<<1|1]+1;
		int skip_layer=1;
		for(int k=0;k<5*5*4;++k)
		{
			if(lossyconv_params[5*5*4*kb+k])
				skip_layer=0;
		}
		if(skip_layer)
			memcpy(p2, p1, sizeof(int[4])*src->iw*src->ih);
		else
		{
			for(int kc=0;kc<4;++kc)
			{
				int depth, half;
				short *curr_filt;

				depth=src->depth[kc];
				if(!depth)
					continue;
				half=1<<depth>>1;
				curr_filt=lossyconv_params+5*5*(4*kb+kc);
				for(int ky=ystart;ky<src->ih-(ystride-1);ky+=ystride)
				{
					for(int kx=xstart;kx<src->iw-(xstride-1);kx+=xstride)
					{
						int idx=src->iw*ky+kx;
						int offset=0, pred=0, vmin=-half, vmax=half-1, uninit=1;
						for(int ky2=-2, idx2=0;ky2<=2;++ky2)
						{
							for(int kx2=-2;kx2<=2;++kx2, ++idx2)
							{
								if((unsigned)(ky+ky2)<(unsigned)src->ih&&(unsigned)(kx+kx2)<(unsigned)src->iw)
								{
									short param=curr_filt[idx2];
									int idx3=src->iw*(ky+ky2)+kx+kx2;
									int nb=p1[idx3<<2|kc];
									if(kc!=1&&lossyconv_causalRCT[kb])
										nb-=p1[idx3<<2|1];
									
									pred+=(param>>1)*nb;
									if(param&1)
									{
										if(uninit)
										{
											uninit=0;
											vmin=nb, vmax=nb;
										}
										else
										{
											UPDATE_MIN(vmin, nb);
											UPDATE_MAX(vmax, nb);
										}
									}
								}
							}
						}
						if(kc!=1&&lossyconv_causalRCT[kb])
							offset+=p1[idx<<2|1];
						pred+=16;
						pred>>=5;
						pred=CLAMP(vmin, pred, vmax);
						pred+=offset;
						pred=CLAMP(-half, pred, half);
#if 0
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
						//if(ky==10&&kx==10)//
						//	printf("");
						int pred=0;
						int vmin=-half, vmax=half-1, uninit=1;
						for(int k=0;k<5*5;++k)
						{
							pred+=(curr_filt[k]>>1)*nb[k];
							if(curr_filt[k]&1)
							{
								if(uninit)
								{
									uninit=0;
									vmin=nb[k], vmax=nb[k];
								}
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
#endif
						{
							int curr=p1[idx<<2|kc];
							curr-=pred;
							curr<<=32-depth;//signed-MA
							curr>>=32-depth;
							p2[idx<<2|kc]=curr;
						}
					}
				}
			}
		}
		{
			int *temp;
			SWAPVAR(p1, p2, temp);
		}
	}
	free(dst);
}