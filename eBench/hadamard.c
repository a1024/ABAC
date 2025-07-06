#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#ifdef _MSC_VER
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif
#include<immintrin.h>
static const char file[]=__FILE__;

static void wht4_fwd(int *x)
{
	x[1]-=x[0];
	x[3]-=x[2];
	x[0]+=x[1]>>1;
	x[2]+=x[3]>>1;

	x[2]-=x[0];
	x[3]-=x[1];
	x[0]+=x[2]>>1;
	x[1]+=x[3]>>1;
}
static void wht4_inv(int *x)
{
	x[1]-=x[3]>>1;
	x[0]-=x[2]>>1;
	x[3]+=x[1];
	x[2]+=x[0];
	
	x[2]-=x[3]>>1;
	x[0]-=x[1]>>1;
	x[3]+=x[2];
	x[1]+=x[0];
}
void image_wht4_fwd(Image *image)
{
	int *temp=(int*)malloc(MAXVAR(image->iw, image->ih)*sizeof(int));
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(temp, 0, MAXVAR(image->iw, image->ih)*sizeof(int));
	for(int kc=0;kc<4;++kc)
	{
		if(!image->depth[kc])
			continue;
#if 1
		for(int ky=0;ky<image->ih;++ky)
		{
			for(int kx=0;kx<image->iw-3;kx+=4)
			{
				int idx=image->iw*ky+kx;
				int x[]=
				{
					image->data[(idx+0)<<2|kc],
					image->data[(idx+1)<<2|kc],
					image->data[(idx+2)<<2|kc],
					image->data[(idx+3)<<2|kc],
				};

				wht4_fwd(x);

				temp[(kx>>2)+(image->iw>>2)*0]=x[0];
				temp[(kx>>2)+(image->iw>>2)*1]=x[1];
				temp[(kx>>2)+(image->iw>>2)*2]=x[2];
				temp[(kx>>2)+(image->iw>>2)*3]=x[3];
			}
			for(int kx=0;kx<image->iw;++kx)
				image->data[(image->iw*ky+kx)<<2|kc]=temp[kx];
		}
#endif
#if 1
		for(int kx=0;kx<image->iw;++kx)
		{
			for(int ky=0;ky<image->ih-3;ky+=4)
			{
				int idx=image->iw*ky+kx;
				int x[]=
				{
					image->data[(idx+image->iw*0)<<2|kc],
					image->data[(idx+image->iw*1)<<2|kc],
					image->data[(idx+image->iw*2)<<2|kc],
					image->data[(idx+image->iw*3)<<2|kc],
				};

				wht4_fwd(x);

				temp[(ky>>2)+(image->ih>>2)*0]=x[0];
				temp[(ky>>2)+(image->ih>>2)*1]=x[1];
				temp[(ky>>2)+(image->ih>>2)*2]=x[2];
				temp[(ky>>2)+(image->ih>>2)*3]=x[3];
			}
			for(int ky=0;ky<image->ih;++ky)
				image->data[(image->iw*ky+kx)<<2|kc]=temp[ky];
		}
#endif
	}
	free(temp);
	for(int kc=0;kc<3;++kc)
	{
		int d=image->depth[kc];
		if(d<24-2)
			d+=2;
		image->depth[kc]=d;
	}
}
void image_wht4_inv(Image *image)
{
	int *temp=(int*)malloc(MAXVAR(image->iw, image->ih)*sizeof(int));
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(temp, 0, MAXVAR(image->iw, image->ih)*sizeof(int));
	for(int kc=0;kc<4;++kc)
	{
		if(!image->depth[kc])
			continue;
#if 1
		for(int kx=0;kx<image->iw;++kx)
		{
			for(int ky=0;ky<image->ih;++ky)
				temp[ky]=image->data[(image->iw*ky+kx)<<2|kc];
			for(int ky=0;ky<image->ih-3;ky+=4)
			{
				int idx=image->iw*ky+kx;
				int x[]=
				{
					temp[(ky>>2)+(image->ih>>2)*0],
					temp[(ky>>2)+(image->ih>>2)*1],
					temp[(ky>>2)+(image->ih>>2)*2],
					temp[(ky>>2)+(image->ih>>2)*3],
				};

				wht4_inv(x);
				
				image->data[(idx+image->iw*0)<<2|kc]=x[0];
				image->data[(idx+image->iw*1)<<2|kc]=x[1];
				image->data[(idx+image->iw*2)<<2|kc]=x[2];
				image->data[(idx+image->iw*3)<<2|kc]=x[3];
			}
		}
#endif
#if 1
		for(int ky=0;ky<image->ih;++ky)
		{
			for(int kx=0;kx<image->iw;++kx)
				temp[kx]=image->data[(image->iw*ky+kx)<<2|kc];
			for(int kx=0;kx<image->iw-3;kx+=4)
			{
				int idx=image->iw*ky+kx;
				int x[]=
				{
					temp[(kx>>2)+(image->iw>>2)*0],
					temp[(kx>>2)+(image->iw>>2)*1],
					temp[(kx>>2)+(image->iw>>2)*2],
					temp[(kx>>2)+(image->iw>>2)*3],
				};

				wht4_inv(x);
				
				image->data[(idx+0)<<2|kc]=x[0];
				image->data[(idx+1)<<2|kc]=x[1];
				image->data[(idx+2)<<2|kc]=x[2];
				image->data[(idx+3)<<2|kc]=x[3];
			}
		}
#endif
	}
	free(temp);
	for(int kc=0;kc<3;++kc)
	{
		int d=image->depth[kc];
		if(d>=image->src_depth[kc]+2)
			d-=2;
		image->depth[kc]=d;
	}
}


static void wht8_fwd(int *x)
{
	x[1]-=x[0];
	x[3]-=x[2];
	x[5]-=x[4];
	x[7]-=x[6];
	x[0]+=x[1]>>1;
	x[2]+=x[3]>>1;
	x[4]+=x[5]>>1;
	x[6]+=x[7]>>1;

	x[2]-=x[0];
	x[3]-=x[1];
	x[6]-=x[4];
	x[7]-=x[5];
	x[0]+=x[2]>>1;
	x[1]+=x[3]>>1;
	x[4]+=x[6]>>1;
	x[5]+=x[7]>>1;

	x[4]-=x[0];
	x[5]-=x[1];
	x[6]-=x[2];
	x[7]-=x[3];
	x[0]+=x[4]>>1;
	x[1]+=x[5]>>1;
	x[2]+=x[6]>>1;
	x[3]+=x[7]>>1;
}
static void wht8_inv(int *x)
{
	x[3]-=x[7]>>1;
	x[2]-=x[6]>>1;
	x[1]-=x[5]>>1;
	x[0]-=x[4]>>1;
	x[7]+=x[3];
	x[6]+=x[2];
	x[5]+=x[1];
	x[4]+=x[0];

	x[5]-=x[7]>>1;
	x[4]-=x[6]>>1;
	x[1]-=x[3]>>1;
	x[0]-=x[2]>>1;
	x[7]+=x[5];
	x[6]+=x[4];
	x[3]+=x[1];
	x[2]+=x[0];

	x[6]-=x[7]>>1;
	x[4]-=x[5]>>1;
	x[2]-=x[3]>>1;
	x[0]-=x[1]>>1;
	x[7]+=x[6];
	x[5]+=x[4];
	x[3]+=x[2];
	x[1]+=x[0];
}
void image_wht8_fwd(Image *image)
{
	int *temp=(int*)malloc(MAXVAR(image->iw, image->ih)*sizeof(int));
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(temp, 0, MAXVAR(image->iw, image->ih)*sizeof(int));
	for(int kc=0;kc<4;++kc)
	{
		if(!image->depth[kc])
			continue;
#if 1
		for(int ky=0;ky<image->ih;++ky)
		{
			for(int kx=0;kx<image->iw-7;kx+=8)
			{
				int idx=image->iw*ky+kx;
				int x[]=
				{
					image->data[(idx+0)<<2|kc],
					image->data[(idx+1)<<2|kc],
					image->data[(idx+2)<<2|kc],
					image->data[(idx+3)<<2|kc],
					image->data[(idx+4)<<2|kc],
					image->data[(idx+5)<<2|kc],
					image->data[(idx+6)<<2|kc],
					image->data[(idx+7)<<2|kc],
				};

				//char y[8];
				//memcpy(y, x, 8);
				//dct8_fwd_i8(y);
				//dct8_inv_i8(y);
				//if(memcmp(x, y, 8))
				//	x[0]=y[0];

				wht8_fwd(x);

				temp[(kx>>3)+(image->iw>>3)*0]=x[0];
				temp[(kx>>3)+(image->iw>>3)*1]=x[1];
				temp[(kx>>3)+(image->iw>>3)*2]=x[2];
				temp[(kx>>3)+(image->iw>>3)*3]=x[3];
				temp[(kx>>3)+(image->iw>>3)*4]=x[4];
				temp[(kx>>3)+(image->iw>>3)*5]=x[5];
				temp[(kx>>3)+(image->iw>>3)*6]=x[6];
				temp[(kx>>3)+(image->iw>>3)*7]=x[7];
			}
			for(int kx=0;kx<image->iw;++kx)
				image->data[(image->iw*ky+kx)<<2|kc]=temp[kx];
		}
#endif
#if 1
		for(int kx=0;kx<image->iw;++kx)
		{
			for(int ky=0;ky<image->ih-7;ky+=8)
			{
				int idx=image->iw*ky+kx;
				int x[]=
				{
					image->data[(idx+image->iw*0)<<2|kc],
					image->data[(idx+image->iw*1)<<2|kc],
					image->data[(idx+image->iw*2)<<2|kc],
					image->data[(idx+image->iw*3)<<2|kc],
					image->data[(idx+image->iw*4)<<2|kc],
					image->data[(idx+image->iw*5)<<2|kc],
					image->data[(idx+image->iw*6)<<2|kc],
					image->data[(idx+image->iw*7)<<2|kc],
				};

				wht8_fwd(x);

				temp[(ky>>3)+(image->ih>>3)*0]=x[0];
				temp[(ky>>3)+(image->ih>>3)*1]=x[1];
				temp[(ky>>3)+(image->ih>>3)*2]=x[2];
				temp[(ky>>3)+(image->ih>>3)*3]=x[3];
				temp[(ky>>3)+(image->ih>>3)*4]=x[4];
				temp[(ky>>3)+(image->ih>>3)*5]=x[5];
				temp[(ky>>3)+(image->ih>>3)*6]=x[6];
				temp[(ky>>3)+(image->ih>>3)*7]=x[7];
			}
			for(int ky=0;ky<image->ih;++ky)
				image->data[(image->iw*ky+kx)<<2|kc]=temp[ky];
		}
#endif
	}
	free(temp);
	for(int kc=0;kc<3;++kc)
	{
		int d=image->depth[kc];
		if(d<24-3)
			d+=3;
		image->depth[kc]=d;
	}
}
void image_wht8_inv(Image *image)
{
	int *temp=(int*)malloc(MAXVAR(image->iw, image->ih)*sizeof(int));
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(temp, 0, MAXVAR(image->iw, image->ih)*sizeof(int));
	for(int kc=0;kc<4;++kc)
	{
		if(!image->depth[kc])
			continue;
#if 1
		for(int kx=0;kx<image->iw;++kx)
		{
			for(int ky=0;ky<image->ih;++ky)
				temp[ky]=image->data[(image->iw*ky+kx)<<2|kc];
			for(int ky=0;ky<image->ih-7;ky+=8)
			{
				int idx=image->iw*ky+kx;
				int x[]=
				{
					temp[(ky>>3)+(image->ih>>3)*0],
					temp[(ky>>3)+(image->ih>>3)*1],
					temp[(ky>>3)+(image->ih>>3)*2],
					temp[(ky>>3)+(image->ih>>3)*3],
					temp[(ky>>3)+(image->ih>>3)*4],
					temp[(ky>>3)+(image->ih>>3)*5],
					temp[(ky>>3)+(image->ih>>3)*6],
					temp[(ky>>3)+(image->ih>>3)*7],
				};

				wht8_inv(x);
				
				image->data[(idx+image->iw*0)<<2|kc]=x[0];
				image->data[(idx+image->iw*1)<<2|kc]=x[1];
				image->data[(idx+image->iw*2)<<2|kc]=x[2];
				image->data[(idx+image->iw*3)<<2|kc]=x[3];
				image->data[(idx+image->iw*4)<<2|kc]=x[4];
				image->data[(idx+image->iw*5)<<2|kc]=x[5];
				image->data[(idx+image->iw*6)<<2|kc]=x[6];
				image->data[(idx+image->iw*7)<<2|kc]=x[7];
			}
		}
#endif
#if 1
		for(int ky=0;ky<image->ih;++ky)
		{
			for(int kx=0;kx<image->iw;++kx)
				temp[kx]=image->data[(image->iw*ky+kx)<<2|kc];
			for(int kx=0;kx<image->iw-7;kx+=8)
			{
				int idx=image->iw*ky+kx;
				int x[]=
				{
					temp[(kx>>3)+(image->iw>>3)*0],
					temp[(kx>>3)+(image->iw>>3)*1],
					temp[(kx>>3)+(image->iw>>3)*2],
					temp[(kx>>3)+(image->iw>>3)*3],
					temp[(kx>>3)+(image->iw>>3)*4],
					temp[(kx>>3)+(image->iw>>3)*5],
					temp[(kx>>3)+(image->iw>>3)*6],
					temp[(kx>>3)+(image->iw>>3)*7],
				};

				wht8_inv(x);
				
				image->data[(idx+0)<<2|kc]=x[0];
				image->data[(idx+1)<<2|kc]=x[1];
				image->data[(idx+2)<<2|kc]=x[2];
				image->data[(idx+3)<<2|kc]=x[3];
				image->data[(idx+4)<<2|kc]=x[4];
				image->data[(idx+5)<<2|kc]=x[5];
				image->data[(idx+6)<<2|kc]=x[6];
				image->data[(idx+7)<<2|kc]=x[7];
			}
		}
#endif
	}
	free(temp);
	for(int kc=0;kc<3;++kc)
	{
		int d=image->depth[kc];
		if(d>=image->src_depth[kc]+3)
			d-=3;
		else
			d=image->src_depth[kc];
		image->depth[kc]=d;
	}
}