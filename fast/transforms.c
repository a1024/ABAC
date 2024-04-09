#include"fast.h"
#include<stdlib.h>
#include<string.h>
static const char file[]=__FILE__;

void rct_JPEG2000_32(Image *image, int fwd)
{
	if(image->nch<3)
		return;
	//char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*image->nch;k<len;k+=image->nch)
		{
			short
				r=image->data[k+0],
				g=image->data[k+1],
				b=image->data[k+2];
			
			r-=g;       //r-g				[1     -1     0  ].RGB
			b-=g;       //b-g				[0     -1     1  ].RGB
			g+=(r+b)>>2;//g+(r-g+b-g)/4 = r/4+g/2+b/4	[1/4    1/2   1/4].RGB

			image->data[k+0]=g;//Y
			image->data[k+1]=b;//Cb
			image->data[k+2]=r;//Cr
		}
	}
	else
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*image->nch;k<len;k+=image->nch)
		{
			short
				Y =image->data[k+0],
				Cb=image->data[k+1],
				Cr=image->data[k+2];
			
			Y-=(Cr+Cb)>>2;
			Cb+=Y;
			Cr+=Y;

			image->data[k+0]=Cr;
			image->data[k+1]=Y;
			image->data[k+2]=Cb;
		}
	}
}
//clamped gradient predictor, aka LOCO-I / Median Edge Detector (MED) predictor from JPEG-LS
void pred_clampgrad(Image *src, int fwd, char *depths)
{
	Image dst={0};
	image_copy(&dst, src);
	if(!dst.data)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	const short *pixels=fwd?dst.data:src->data;
	int dy=src->nch*src->iw;
	int dx=src->nch;
	for(int kc=0;kc<src->nch;++kc)
	{
		int nlevels=1<<depths[kc];
		for(int ky=0, idx=kc;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, idx+=src->nch)
			{
				int
					NW=kx&&ky	?pixels[idx-dy	-dx	]:0,
					N =ky		?pixels[idx	-dy	]:0,
					W =kx		?pixels[idx-dx		]:0;
				int pred=N+W-NW;
				pred=MEDIAN3(N, W, pred);

				pred^=-fwd;
				pred+=fwd;

				pred+=src->data[idx];

				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;

				src->data[idx]=pred;
			}
		}
	}
	image_clear(&dst);
}