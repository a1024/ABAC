#include"best.h"
#include<stdio.h>
#include<string.h>
int get_nch(const char *buf, int res)//returns nch = {0 degenerate, 1 gray, 2 gray_alpha, 3, rgb, 4, rgb_alpha}
{
	int has_image=0, gray=1, has_alpha=0;
	int initial=buf[0];
	for(int k=0;k<res;++k)
	{
		char r=buf[k<<2|0], g=buf[k<<2|1], b=buf[k<<2|2], a=buf[k<<2|3];
		if(r!=initial)
			has_image=1;
		if(r!=g||g!=b)
			gray=0;
		if(a!=-1)
			has_alpha=1;
	}
	if(!has_image)
		return 0;
	int nch=gray+!gray*3+has_alpha;
	return nch;
}
int get_nch32(const int *buf, int res)//returns nch = {0 degenerate, 1 gray, 2 gray_alpha, 3, rgb, 4, rgb_alpha}
{
	int has_image=0, gray=1, has_alpha=0;
	int initial=buf[0];
	for(int k=0;k<res;++k)
	{
		int r=buf[k<<2|0], g=buf[k<<2|1], b=buf[k<<2|2], a=buf[k<<2|3];
		if(r!=initial)
			has_image=1;
		if(r!=g||g!=b)
			gray=0;
		if(a)
			has_alpha=1;
	}
	if(!has_image)
		return 0;
	int nch=1+!gray*2+has_alpha;
	return nch;
}
void addbuf(unsigned char *buf, int iw, int ih, int nch, int bytestride, int ammount)
{
	for(int kp=0, len=iw*ih*bytestride;kp<len;kp+=bytestride)
	{
		for(int kc=0;kc<nch;++kc)
			buf[kp+kc]+=ammount;
	}
}
int compare_bufs_32(const int *b1, const int *b0, int iw, int ih, int nch, int chstride, const char *name, int backward, int loud)
{
	ptrdiff_t len=(ptrdiff_t)chstride*iw*ih;
	int inc=chstride*(1-(backward<<1));
	for(ptrdiff_t k=backward?len-chstride:0;k>=0&&k<len;k+=inc)
	{
		if(memcmp(b1+k, b0+k, nch*sizeof(int)))
		{
			if(loud)
			{
				ptrdiff_t idx=k/chstride, kx=idx%iw, ky=idx/iw;
				printf("%s error XY (%5lld, %5lld) / %5d x %5d  b1 != b0\n", name, kx, ky, iw, ih);
				for(int kc=0;kc<nch;++kc)
					printf("C%d  0x%08X != 0x%08X    %d != %d\n",
						kc, (unsigned)b1[k+kc], (unsigned)b0[k+kc], (unsigned)b1[k+kc], (unsigned)b0[k+kc]
					);
			}
			return 1;
		}
	}
	if(loud)
		printf("%s:\tSUCCESS\n", name);
	return 0;
}
int compare_bufs_uint8(const unsigned char *b1, const unsigned char *b0, int iw, int ih, int symbytes, int bytestride, const char *name, int backward, int loud)
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
			if(loud)
			{
				ptrdiff_t idx=k/bytestride, kx=idx%iw, ky=idx/iw;
				printf("%s error XY (%5lld, %5lld) / %5d x %5d  b1 != b0\n", name, kx, ky, iw, ih);
				for(int kc=0;kc<symbytes;++kc)
					printf("C%d  0x%02X != 0x%02X    %d != %d\n", kc, (unsigned)b1[k+kc], (unsigned)b0[k+kc], (unsigned)b1[k+kc], (unsigned)b0[k+kc]);
			}
			//if(backward)
			//	printf("%s error at %d - %d: 0x%02X != 0x%02X\n", name, (int)len-1, (int)(len-1-k), b1[k], b0[k]);
			//else
			//	printf("%s error at %d: 0x%02X != 0x%02X\n", name, (int)k, b1[k], b0[k]);
			return 1;
		}
	}
	if(loud)
		printf("%s:\tSUCCESS\n", name);
	return 0;
}

void rct_JPEG2000_32(Image *image, int fwd)
{
	//char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			r-=g;       //r-g				[1     -1     0  ].RGB
			b-=g;       //b-g				[0     -1     1  ].RGB
			g+=(r+b)>>2;//g+(r-g+b-g)/4 = r/4+g/2+b/4	[1/4    1/2   1/4].RGB

			image->data[k  ]=g;//Y
			image->data[k|1]=b;//Cb
			image->data[k|2]=r;//Cr
		}
		//ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		//image->depth[1]+=image->depth[1]<24;
		//image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		//image->depth[1]-=image->depth[1]>image->src_depth[1];
		//image->depth[2]-=image->depth[2]>image->src_depth[2];
		//ROTATE3(image->depth[2], image->depth[1], image->depth[0], temp);
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y=image->data[k], Cb=image->data[k|1], Cr=image->data[k|2];
			
			Y-=(Cr+Cb)>>2;
			Cb+=Y;
			Cr+=Y;

			image->data[k  ]=Cr;
			image->data[k|1]=Y;
			image->data[k|2]=Cb;
		}
	}
}
void colortransform_JPEG2000_fwd(char *buf, int iw, int ih)
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;
		b-=g;
		g+=(r+b)>>2;

		buf[k  ]=g;//luma
		buf[k|1]=b;
		buf[k|2]=r;
	}
}
void colortransform_JPEG2000_inv(char *buf, int iw, int ih)
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char g=buf[k], b=buf[k|1], r=buf[k|2];
		
		g-=(r+b)>>2;
		b+=g;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_YCbCr_R_v0_fwd(char *buf, int iw, int ih)
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;   //Co = diff(r, g)
		g+=r>>1;//     av(r, g)				g+floor((r-g)/2) = floor((r+g)/2)
		b-=g;   //Cb = diff(b, av(r, g))
		g+=b>>1;//Y  = av(b, av(r, g))

		buf[k  ]=r;//Co
		buf[k|1]=g;//Y
		buf[k|2]=b;//Cb
	}
}
void colortransform_YCbCr_R_v0_inv(char *buf, int iw, int ih)
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];//Co Y Cb
		
		g-=b>>1;
		b+=g;
		g-=r>>1;
		r+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}
void colortransform_subgreen_fwd(char *buf, int iw, int ih)
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=buf[k], g=buf[k|1], b=buf[k|2];

		r-=g;
		b-=g;

		buf[k  ]=g;//luma
		buf[k|1]=b;
		buf[k|2]=r;
	}
}
void colortransform_subgreen_inv(char *buf, int iw, int ih)
{
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char g=buf[k], b=buf[k|1], r=buf[k|2];
		
		r+=g;
		b+=g;

		buf[k  ]=r;
		buf[k|1]=g;
		buf[k|2]=b;
	}
}

void pack3_fwd(char *buf, int res)
{
	for(int k=1, idx=3;k<res;++k, idx+=3)
	{
		char r=buf[k<<2|0], g=buf[k<<2|1], b=buf[k<<2|2];
		buf[idx+0]=r;
		buf[idx+1]=g;
		buf[idx+2]=b;
	}
}
void pack3_inv(char *buf, int res)
{
	for(int k=res-1, idx=3*(res-1);k>0;--k, idx-=3)
	{
		char r=buf[idx+0], g=buf[idx+1], b=buf[idx+2];
		buf[k<<2|0]=r;
		buf[k<<2|1]=g;
		buf[k<<2|2]=b;
		buf[k<<2|3]=(char)-1;
	}
}