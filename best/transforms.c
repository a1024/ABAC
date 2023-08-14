#include<stddef.h>
void addbuf(unsigned char *buf, int iw, int ih, int nch, int bytestride, int ammount)
{
	for(int kp=0, len=iw*ih*bytestride;kp<len;kp+=bytestride)
	{
		for(int kc=0;kc<nch;++kc)
			buf[kp+kc]+=ammount;
	}
}
void colortransform_ycocb_fwd(char *buf, int iw, int ih)//3 channels, stride 4 bytes
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
void colortransform_ycocb_inv(char *buf, int iw, int ih)//3 channels, stride 4 bytes
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