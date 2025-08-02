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
			for(int kx=0;kx<(image->iw&~3);++kx)
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
			for(int ky=0;ky<(image->ih&~3);++ky)
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
	for(int kc=0;kc<3;++kc)
	{
		int d=image->depth[kc];
		if(d>=image->src_depth[kc]+2)
			d-=2;
		image->depth[kc]=d;
	}
	for(int kc=0;kc<4;++kc)
	{
		int vmin=-(1<<image->depth[kc]>>1), vmax=(1<<image->depth[kc]>>1)-1;
		if(!image->depth[kc])
			continue;
#if 1
		for(int kx=0;kx<image->iw;++kx)
		{
			for(int ky=0;ky<(image->ih&~3);++ky)
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
			for(int kx=0;kx<(image->iw&~3);++kx)
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
				CLAMP2(x[0], vmin, vmax);
				CLAMP2(x[1], vmin, vmax);
				CLAMP2(x[2], vmin, vmax);
				CLAMP2(x[3], vmin, vmax);
				
				image->data[(idx+0)<<2|kc]=x[0];
				image->data[(idx+1)<<2|kc]=x[1];
				image->data[(idx+2)<<2|kc]=x[2];
				image->data[(idx+3)<<2|kc]=x[3];
			}
		}
#endif
	}
	free(temp);
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
			for(int kx=0;kx<(image->iw&~7);++kx)
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
			for(int ky=0;ky<(image->ih&~7);++ky)
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
	for(int kc=0;kc<3;++kc)
	{
		int d=image->depth[kc];
		if(d>=image->src_depth[kc]+3)
			d-=3;
		else
			d=image->src_depth[kc];
		image->depth[kc]=d;
	}
	for(int kc=0;kc<4;++kc)
	{
		int vmin=-(1<<image->depth[kc]>>1), vmax=(1<<image->depth[kc]>>1)-1;
		if(!image->depth[kc])
			continue;
#if 1
		for(int kx=0;kx<image->iw;++kx)
		{
			for(int ky=0;ky<(image->ih&~7);++ky)
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
			for(int kx=0;kx<(image->iw&~7);++kx)
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
				CLAMP2(x[0], vmin, vmax);
				CLAMP2(x[1], vmin, vmax);
				CLAMP2(x[2], vmin, vmax);
				CLAMP2(x[3], vmin, vmax);
				CLAMP2(x[4], vmin, vmax);
				CLAMP2(x[5], vmin, vmax);
				CLAMP2(x[6], vmin, vmax);
				CLAMP2(x[7], vmin, vmax);
				
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
}


static void haar8_fwd(int *x)
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
	x[6]-=x[4];
	x[0]+=x[2]>>1;
	x[4]+=x[6]>>1;

	x[4]-=x[0];
	x[0]+=x[4]>>1;
}
static void haar8_inv(int *x)
{
	x[0]-=x[4]>>1;
	x[4]+=x[0];

	x[4]-=x[6]>>1;
	x[0]-=x[2]>>1;
	x[6]+=x[4];
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
void image_haar8_fwd(Image *image)
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

				haar8_fwd(x);

				temp[(kx>>3)+(image->iw>>3)*0]=x[0];
				temp[(kx>>3)+(image->iw>>3)*1]=x[1];
				temp[(kx>>3)+(image->iw>>3)*2]=x[2];
				temp[(kx>>3)+(image->iw>>3)*3]=x[3];
				temp[(kx>>3)+(image->iw>>3)*4]=x[4];
				temp[(kx>>3)+(image->iw>>3)*5]=x[5];
				temp[(kx>>3)+(image->iw>>3)*6]=x[6];
				temp[(kx>>3)+(image->iw>>3)*7]=x[7];
			}
			for(int kx=0;kx<(image->iw&~7);++kx)
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

				haar8_fwd(x);

				temp[(ky>>3)+(image->ih>>3)*0]=x[0];
				temp[(ky>>3)+(image->ih>>3)*1]=x[1];
				temp[(ky>>3)+(image->ih>>3)*2]=x[2];
				temp[(ky>>3)+(image->ih>>3)*3]=x[3];
				temp[(ky>>3)+(image->ih>>3)*4]=x[4];
				temp[(ky>>3)+(image->ih>>3)*5]=x[5];
				temp[(ky>>3)+(image->ih>>3)*6]=x[6];
				temp[(ky>>3)+(image->ih>>3)*7]=x[7];
			}
			for(int ky=0;ky<(image->ih&~7);++ky)
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
void image_haar8_inv(Image *image)
{
	int *temp=(int*)malloc(MAXVAR(image->iw, image->ih)*sizeof(int));
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(temp, 0, MAXVAR(image->iw, image->ih)*sizeof(int));
	for(int kc=0;kc<3;++kc)
	{
		int d=image->depth[kc];
		if(d>=image->src_depth[kc]+2)
			d-=2;
		else
			d=image->src_depth[kc];
		image->depth[kc]=d;
	}
	for(int kc=0;kc<4;++kc)
	{
		int vmin=-(1<<image->depth[kc]>>1), vmax=(1<<image->depth[kc]>>1)-1;
		if(!image->depth[kc])
			continue;
#if 1
		for(int kx=0;kx<image->iw;++kx)
		{
			for(int ky=0;ky<(image->ih&~7);++ky)
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

				haar8_inv(x);
				
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
			for(int kx=0;kx<(image->iw&~7);++kx)
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

				haar8_inv(x);
				CLAMP2(x[0], vmin, vmax);
				CLAMP2(x[1], vmin, vmax);
				CLAMP2(x[2], vmin, vmax);
				CLAMP2(x[3], vmin, vmax);
				CLAMP2(x[4], vmin, vmax);
				CLAMP2(x[5], vmin, vmax);
				CLAMP2(x[6], vmin, vmax);
				CLAMP2(x[7], vmin, vmax);
				
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
}


//#define DCT_SH 7	// <= 7
static void fdct8_fwd(int *x)
{
	//The BinDCT: Fast Multiplierless Approximation of the DCT - Trac D Tran - IEEE SPL 2000-06
#if 1
	float
		t0=(float)x[0],
		t1=(float)x[1],
		t2=(float)x[2],
		t3=(float)x[3],
		t4=(float)x[4],
		t5=(float)x[5],
		t6=(float)x[6],
		t7=(float)x[7];

	t7=t0-t7;
	t6=t1-t6;
	t5=t2-t5;
	t4=t3-t4;
	t0=2*t0-t7;//2*a-(a-b) = a+b
	t1=2*t1-t6;
	t2=2*t2-t5;
	t3=2*t3-t4;
	{
		float a=t5, b=t6;
		t5=(-a+b)*0.7071067811865475244008443621048f;//cos(pi/4)
		t6=(+a+b)*0.7071067811865475244008443621048f;
	}
	t3=t0-t3;
	t2=t1-t2;
	t5=t4-t5;
	t6=t7-t6;
	t0=2*t0-t3;
	t1=2*t1-t2;
	t4=2*t4-t5;
	t7=2*t7-t6;
	{
		float a0=t0, b0=t1;
		float a1=t2, b1=t3;
		float a2=t4, b2=t7;
		float a3=t5, b3=t6;

		t0=(a0+b0)*0.7071067811865475244008443621048f;
		t4=(a0-b0)*0.7071067811865475244008443621048f;

		t2=a1*+0.9238795325112867561281831893968f + b1*+0.3826834323650897717284599840304f;//+cos(3pi/8), +sin(3pi/8)
		t6=a1*-0.3826834323650897717284599840304f + b1*+0.9238795325112867561281831893968f;//-sin(3pi/8), +cos(3pi/8)

		t1=a2*+0.1950903220161282678482848684770f + b2*+0.9807852804032304491261822361342f;//+cos(7pi/16), +sin(7pi/16)
		t7=a2*-0.9807852804032304491261822361342f + b2*+0.1950903220161282678482848684770f;//-sin(7pi/16), +cos(7pi/16)

		t5=a3*+0.8314696123025452370787883776179f + b3*+0.5555702330196022247428308139485f;//+cos(3pi/16), +sin(3pi/16)
		t3=a3*-0.5555702330196022247428308139485f + b3*+0.8314696123025452370787883776179f;//-sin(3pi/16), +cos(3pi/16)
	}

	x[0]=CVTFP32_I32(t0);
	x[1]=CVTFP32_I32(t1);
	x[2]=CVTFP32_I32(t2);
	x[3]=CVTFP32_I32(t3);
	x[4]=CVTFP32_I32(t4);
	x[5]=CVTFP32_I32(t5);
	x[6]=CVTFP32_I32(t6);
	x[7]=CVTFP32_I32(t7);
#endif
#if 0
	int
		t0=x[0]<<DCT_SH,
		t1=x[1]<<DCT_SH,
		t2=x[2]<<DCT_SH,
		t3=x[3]<<DCT_SH,
		t4=x[4]<<DCT_SH,
		t5=x[5]<<DCT_SH,
		t6=x[6]<<DCT_SH,
		t7=x[7]<<DCT_SH;
#if 1
	t7=t0-t7;
	t6=t1-t6;
	t5=t2-t5;
	t4=t3-t4;
	t0=2*t0-t7;//2*a-(a-b) = a+b
	t1=2*t1-t6;
	t2=2*t2-t5;
	t3=2*t3-t4;
	{
		int a=t5, b=t6;
		t5=(-a+b)*(int)(0.7071067811865475244008443621048*(1<<12))>>12;//cos(pi/4)
		t6=(+a+b)*(int)(0.7071067811865475244008443621048*(1<<12))>>12;
	}
	t3=t0-t3;
	t2=t1-t2;
	t5=t4-t5;
	t6=t7-t6;
	t0=2*t0-t3;
	t1=2*t1-t2;
	t4=2*t4-t5;
	t7=2*t7-t6;
#endif
#if 1
	{
		int a0=t0, b0=t1;
		int a1=t2, b1=t3;
		int a2=t4, b2=t7;
		int a3=t5, b3=t6;

		t0=(a0+b0)*(int)(0.7071067811865475244008443621048*(1<<12))>>12;
		t4=(a0-b0)*(int)(0.7071067811865475244008443621048*(1<<12))>>12;

		t2=(a1*(int)(+0.9238795325112867561281831893968*(1<<12)+0.5) + b1*(int)(+0.3826834323650897717284599840304*(1<<12)+0.5))>>12;//+cos(3pi/8), +sin(3pi/8)
		t6=(a1*(int)(-0.3826834323650897717284599840304*(1<<12)+0.5) + b1*(int)(+0.9238795325112867561281831893968*(1<<12)+0.5))>>12;//-sin(3pi/8), +cos(3pi/8)

		t1=(a2*(int)(+0.1950903220161282678482848684770*(1<<12)+0.5) + b2*(int)(+0.9807852804032304491261822361342*(1<<12)+0.5))>>12;//+cos(7pi/16), +sin(7pi/16)
		t7=(a2*(int)(-0.9807852804032304491261822361342*(1<<12)+0.5) + b2*(int)(+0.1950903220161282678482848684770*(1<<12)+0.5))>>12;//-sin(7pi/16), +cos(7pi/16)

		t5=(a3*(int)(+0.8314696123025452370787883776179*(1<<12)+0.5) + b3*(int)(+0.5555702330196022247428308139485*(1<<12)+0.5))>>12;//+cos(3pi/16), +sin(3pi/16)
		t3=(a3*(int)(-0.5555702330196022247428308139485*(1<<12)+0.5) + b3*(int)(+0.8314696123025452370787883776179*(1<<12)+0.5))>>12;//-sin(3pi/16), +cos(3pi/16)
	}
#endif

	x[0]=t0>>DCT_SH;
	x[1]=t1>>DCT_SH;
	x[2]=t2>>DCT_SH;
	x[3]=t3>>DCT_SH;
	x[4]=t4>>DCT_SH;
	x[5]=t5>>DCT_SH;
	x[6]=t6>>DCT_SH;
	x[7]=t7>>DCT_SH;
#endif
#if 0
	x[7]=x[0]-x[7];
	x[6]=x[1]-x[6];
	x[5]=x[2]-x[5];
	x[4]=x[3]-x[4];
	x[0]=2*x[0]-x[7];//2*a-(a-b) = a+b
	x[1]=2*x[1]-x[6];
	x[2]=2*x[2]-x[5];
	x[3]=2*x[3]-x[4];

	{
		int a=x[5], b=x[6];
		x[5]=(b-a)*(int)(0.7071067811865475244008443621048*(1<<12))>>12;//cos(pi/4)
		x[6]=(b+a)*(int)(0.7071067811865475244008443621048*(1<<12))>>12;
	}

	x[3]=x[0]-x[3];
	x[2]=x[1]-x[2];
	x[5]=x[4]-x[5];
	x[6]=x[7]-x[6];
	x[0]=2*x[0]-x[3];
	x[1]=2*x[1]-x[2];
	x[4]=2*x[4]-x[5];
	x[7]=2*x[7]-x[6];

	{
		int a0=x[0], b0=x[1];
		int a1=x[2], b1=x[3];
		int a2=x[4], b2=x[7];
		int a3=x[5], b3=x[6];

		x[0]=(a0+b0)*(int)(0.7071067811865475244008443621048*(1<<12))>>12;
		x[4]=(a0-b0)*(int)(0.7071067811865475244008443621048*(1<<12))>>12;
		
		x[2]=(a1*(int)(+0.9238795325112867561281831893968*(1<<12)+0.5) + b1*(int)(+0.3826834323650897717284599840304*(1<<12)+0.5))>>12;//+cos(3pi/8), +sin(3pi/8)
		x[6]=(a1*(int)(-0.3826834323650897717284599840304*(1<<12)+0.5) + b1*(int)(+0.9238795325112867561281831893968*(1<<12)+0.5))>>12;//-sin(3pi/8), +cos(3pi/8)

		x[1]=(a2*(int)(+0.1950903220161282678482848684770*(1<<12)+0.5) + b2*(int)(+0.9807852804032304491261822361342*(1<<12)+0.5))>>12;//+cos(7pi/16), +sin(7pi/16)
		x[7]=(a2*(int)(-0.9807852804032304491261822361342*(1<<12)+0.5) + b2*(int)(+0.1950903220161282678482848684770*(1<<12)+0.5))>>12;//-sin(7pi/16), +cos(7pi/16)

		x[5]=(a3*(int)(+0.8314696123025452370787883776179*(1<<12)+0.5) + b3*(int)(+0.5555702330196022247428308139485*(1<<12)+0.5))>>12;//+cos(3pi/16), +sin(3pi/16)
		x[3]=(a3*(int)(-0.5555702330196022247428308139485*(1<<12)+0.5) + b3*(int)(+0.8314696123025452370787883776179*(1<<12)+0.5))>>12;//-sin(3pi/16), +cos(3pi/16)
	}
#endif
}
static void fdct8_inv(int *x)
{
	float
		a0=(float)x[0], b0=(float)x[4],
		a1=(float)x[2], b1=(float)x[6],
		a2=(float)x[1], b2=(float)x[7],
		a3=(float)x[5], b3=(float)x[3];

	float t0, t1, t2, t3, t4, t5, t6, t7;

	t0=(a0+b0)*0.7071067811865475244008443621048f;
	t1=(a0-b0)*0.7071067811865475244008443621048f;

	t2=a1*+0.9238795325112867561281831893968f + b1*-0.3826834323650897717284599840304f;//+cos(3pi/8), +sin(3pi/8)
	t3=a1*+0.3826834323650897717284599840304f + b1*+0.9238795325112867561281831893968f;//-sin(3pi/8), +cos(3pi/8)

	t4=a2*+0.1950903220161282678482848684770f + b2*-0.9807852804032304491261822361342f;//+cos(7pi/16), +sin(7pi/16)
	t7=a2*+0.9807852804032304491261822361342f + b2*+0.1950903220161282678482848684770f;//-sin(7pi/16), +cos(7pi/16)

	t5=a3*+0.8314696123025452370787883776179f + b3*-0.5555702330196022247428308139485f;//+cos(3pi/16), +sin(3pi/16)
	t6=a3*+0.5555702330196022247428308139485f + b3*+0.8314696123025452370787883776179f;//-sin(3pi/16), +cos(3pi/16)

	t0+=t3;//a = (s+d)/2
	t1+=t2;
	t4+=t5;
	t7+=t6;
	t0*=0.5;
	t1*=0.5;
	t4*=0.5;
	t7*=0.5;
	t3=t0-t3;//b = a - (a-b)
	t2=t1-t2;
	t5=t4-t5;
	t6=t7-t6;
	{
		float a=t5, b=t6;
		t5=(-a+b)*0.7071067811865475244008443621048f;//cos(pi/4)
		t6=(+a+b)*0.7071067811865475244008443621048f;
	}
	t0+=t7;//a = (s+d)/2
	t1+=t6;
	t2+=t5;
	t3+=t4;
	t0*=0.5;
	t1*=0.5;
	t2*=0.5;
	t3*=0.5;
	t7=t0-t7;//b = a - (a-b)
	t6=t1-t6;
	t5=t2-t5;
	t4=t3-t4;

	x[0]=CVTFP32_I32(t0);
	x[1]=CVTFP32_I32(t1);
	x[2]=CVTFP32_I32(t2);
	x[3]=CVTFP32_I32(t3);
	x[4]=CVTFP32_I32(t4);
	x[5]=CVTFP32_I32(t5);
	x[6]=CVTFP32_I32(t6);
	x[7]=CVTFP32_I32(t7);
#if 0
	int
		a0=x[0]<<DCT_SH, b0=x[4]<<DCT_SH,
		a1=x[2]<<DCT_SH, b1=x[6]<<DCT_SH,
		a2=x[1]<<DCT_SH, b2=x[7]<<DCT_SH,
		a3=x[5]<<DCT_SH, b3=x[3]<<DCT_SH;

	int t0, t1, t2, t3, t4, t5, t6, t7;
#if 1
	t0=(a0+b0)*(int)(0.7071067811865475244008443621048*(1<<12))>>12;
	t1=(a0-b0)*(int)(0.7071067811865475244008443621048*(1<<12))>>12;

	t2=(a1*(int)(+0.9238795325112867561281831893968*(1<<12)+0.5) + b1*(int)(-0.3826834323650897717284599840304*(1<<12)+0.5))>>12;//+cos(3pi/8), +sin(3pi/8)
	t3=(a1*(int)(+0.3826834323650897717284599840304*(1<<12)+0.5) + b1*(int)(+0.9238795325112867561281831893968*(1<<12)+0.5))>>12;//-sin(3pi/8), +cos(3pi/8)

	t4=(a2*(int)(+0.1950903220161282678482848684770*(1<<12)+0.5) + b2*(int)(-0.9807852804032304491261822361342*(1<<12)+0.5))>>12;//+cos(7pi/16), +sin(7pi/16)
	t7=(a2*(int)(+0.9807852804032304491261822361342*(1<<12)+0.5) + b2*(int)(+0.1950903220161282678482848684770*(1<<12)+0.5))>>12;//-sin(7pi/16), +cos(7pi/16)

	t5=(a3*(int)(+0.8314696123025452370787883776179*(1<<12)+0.5) + b3*(int)(-0.5555702330196022247428308139485*(1<<12)+0.5))>>12;//+cos(3pi/16), +sin(3pi/16)
	t6=(a3*(int)(+0.5555702330196022247428308139485*(1<<12)+0.5) + b3*(int)(+0.8314696123025452370787883776179*(1<<12)+0.5))>>12;//-sin(3pi/16), +cos(3pi/16)
#else
	t0=x[0]<<DCT_SH;
	t1=x[1]<<DCT_SH;
	t2=x[2]<<DCT_SH;
	t3=x[3]<<DCT_SH;
	t4=x[4]<<DCT_SH;
	t5=x[5]<<DCT_SH;
	t6=x[6]<<DCT_SH;
	t7=x[7]<<DCT_SH;
#endif
#if 1
	t0+=t3;//a = (s+d)/2
	t1+=t2;
	t4+=t5;
	t7+=t6;
	t0>>=1;
	t1>>=1;
	t4>>=1;
	t7>>=1;
	t3=t0-t3;//b = a - (a-b)
	t2=t1-t2;
	t5=t4-t5;
	t6=t7-t6;
	{
		int a=t5, b=t6;
		t5=(-a+b)*(int)(0.7071067811865475244008443621048*(1<<12))>>12;//cos(pi/4)
		t6=(+a+b)*(int)(0.7071067811865475244008443621048*(1<<12))>>12;
	}
	t0+=t7;//a = (s+d)/2
	t1+=t6;
	t2+=t5;
	t3+=t4;
	t0>>=1;
	t1>>=1;
	t2>>=1;
	t3>>=1;
	t7=t0-t7;//b = a - (a-b)
	t6=t1-t6;
	t5=t2-t5;
	t4=t3-t4;
#endif

	x[0]=t0>>DCT_SH;
	x[1]=t1>>DCT_SH;
	x[2]=t2>>DCT_SH;
	x[3]=t3>>DCT_SH;
	x[4]=t4>>DCT_SH;
	x[5]=t5>>DCT_SH;
	x[6]=t6>>DCT_SH;
	x[7]=t7>>DCT_SH;
#endif
}
void image_fdct8_fwd(Image *image)
{
	//{
	//	int x1[8], x2[8];
	//	for(int k=0;k<8;++k)
	//		x1[k]=(rand()&255)-128;
	//	memcpy(x2, x1, sizeof(x2));
	//	fdct8_fwd(x2);
	//	fdct8_inv(x2);
	//	if(memcmp(x1, x2, sizeof(x1)))
	//		LOG_ERROR("");
	//}
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

				fdct8_fwd(x);

				temp[(kx>>3)+(image->iw>>3)*0]=x[0];
				temp[(kx>>3)+(image->iw>>3)*1]=x[1];
				temp[(kx>>3)+(image->iw>>3)*2]=x[2];
				temp[(kx>>3)+(image->iw>>3)*3]=x[3];
				temp[(kx>>3)+(image->iw>>3)*4]=x[4];
				temp[(kx>>3)+(image->iw>>3)*5]=x[5];
				temp[(kx>>3)+(image->iw>>3)*6]=x[6];
				temp[(kx>>3)+(image->iw>>3)*7]=x[7];
			}
			for(int kx=0;kx<(image->iw&~7);++kx)
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

				fdct8_fwd(x);
				x[0]>>=3;
				x[1]>>=3;
				x[2]>>=3;
				x[3]>>=3;
				x[4]>>=3;
				x[5]>>=3;
				x[6]>>=3;
				x[7]>>=3;

				temp[(ky>>3)+(image->ih>>3)*0]=x[0];
				temp[(ky>>3)+(image->ih>>3)*1]=x[1];
				temp[(ky>>3)+(image->ih>>3)*2]=x[2];
				temp[(ky>>3)+(image->ih>>3)*3]=x[3];
				temp[(ky>>3)+(image->ih>>3)*4]=x[4];
				temp[(ky>>3)+(image->ih>>3)*5]=x[5];
				temp[(ky>>3)+(image->ih>>3)*6]=x[6];
				temp[(ky>>3)+(image->ih>>3)*7]=x[7];
			}
			for(int ky=0;ky<(image->ih&~7);++ky)
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
void image_fdct8_inv(Image *image)
{
	int *temp=(int*)malloc(MAXVAR(image->iw, image->ih)*sizeof(int));
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(temp, 0, MAXVAR(image->iw, image->ih)*sizeof(int));
	for(int kc=0;kc<3;++kc)
	{
		int d=image->depth[kc];
		if(d>=image->src_depth[kc]+3)
			d-=3;
		else
			d=image->src_depth[kc];
		image->depth[kc]=d;
	}
	for(int kc=0;kc<4;++kc)
	{
		int vmin=-(1<<image->depth[kc]>>1), vmax=(1<<image->depth[kc]>>1)-1;
		if(!image->depth[kc])
			continue;
#if 1
		for(int kx=0;kx<image->iw;++kx)
		{
			for(int ky=0;ky<(image->ih&~7);++ky)
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
				x[0]<<=3;
				x[1]<<=3;
				x[2]<<=3;
				x[3]<<=3;
				x[4]<<=3;
				x[5]<<=3;
				x[6]<<=3;
				x[7]<<=3;

				fdct8_inv(x);
				
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
			for(int kx=0;kx<(image->iw&~7);++kx)
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

				fdct8_inv(x);
				CLAMP2(x[0], vmin, vmax);
				CLAMP2(x[1], vmin, vmax);
				CLAMP2(x[2], vmin, vmax);
				CLAMP2(x[3], vmin, vmax);
				CLAMP2(x[4], vmin, vmax);
				CLAMP2(x[5], vmin, vmax);
				CLAMP2(x[6], vmin, vmax);
				CLAMP2(x[7], vmin, vmax);
				
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
}