#include<stddef.h>
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
void addbuf(unsigned char *buf, int iw, int ih, int nch, int bytestride, int ammount)
{
	for(int kp=0, len=iw*ih*bytestride;kp<len;kp+=bytestride)
	{
		for(int kc=0;kc<nch;++kc)
			buf[kp+kc]+=ammount;
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