#include"ebench.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<process.h>
#include"intercept_malloc.h"
#include"lodepng.h"//
#include<immintrin.h>
#ifdef _MSC_VER
#include<intrin.h>
#elif defined __GNUC__
#include<x86intrin.h>
#endif
static const char file[]=__FILE__;

void calc_histogram(const int *buf, int iw, int ih, int kc, int x1, int x2, int y1, int y2, int depth, int *hist, int *hist8)
{
	int nlevels=1<<depth;
	int shift=MAXVAR(0, depth-8);

	(void)ih;

	memset(hist, 0, nlevels*sizeof(int));
	if(hist8)
		memset(hist8, 0, 256*sizeof(int));
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			int sym=buf[(iw*ky+kx)<<2|kc]+(nlevels>>1);
			sym&=nlevels-1;
			//sym=CLAMP(0, sym, nlevels-1);
			++hist[sym];
			if(hist8)
				++hist8[sym>>shift];
		}
	}
}
double calc_entropy(const int *hist, int nlevels, int sum)//nlevels = 1<<current_depth
{
	double entropy=0;
	if(sum<0)
	{
		sum=0;
		for(int k=0;k<nlevels;++k)
			sum+=hist[k];
	}
	if(!sum)
		return 0;//degenerate distribution
	for(int k=0;k<nlevels;++k)
	{
		int freq=hist[k];
		if(freq)
		{
			double p=(double)freq/sum;
			entropy-=p*log2(p);		//Shannon's law
		}
	}
	return entropy;//invCR = entropy/original_depth
}
int calc_maxdepth(Image const *image, int *inflation)
{
	int maxdepth=image->depth[0]+(inflation?inflation[0]:0);
	int next=image->depth[1]+(inflation?inflation[1]:0);
	if(maxdepth<next)
		maxdepth=next;
	next=image->depth[2]+(inflation?inflation[2]:0);
	if(maxdepth<next)
		maxdepth=next;
	return maxdepth;
}
void calc_depthfromdata(int *image, int iw, int ih, char *depths, const char *src_depths)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	for(int kc=0;kc<3;++kc)
	{
		int vmin=image[0<<2|kc], vmax=image[0<<2|kc];
		for(ptrdiff_t k=1;k<res;++k)
		{
			int sym=image[k<<2|kc];
			if(vmin>sym)
				vmin=sym;
			if(vmax<sym)
				vmax=sym;
		}
		++vmax;
		vmin=abs(vmin);
		if(vmax<vmin)
			vmax=vmin;
		{
			int nlevels=vmax<<1;
			depths[kc]=(char)ceil_log2_32(nlevels);
		}
		depths[kc]=MAXVAR(depths[kc], src_depths[kc]);
	}
}


//color transforms
void colortransform_YCoCg_R(Image *image, int fwd)
{
	char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];

			r-=b;		//co = r-b			diff(r, b)
			b+=r>>1;	//(r+b)/2
			g-=b;		//cg = g-(r+b)/2		diff(g, av(r, b))
			b+=g>>1;	//Y  = (r+b)/2 + (g-(r+b)/2)/2 = r/4+g/2+b/4	av(g, av(r, b))

			image->data[k  ]=b;//Y
			image->data[k|1]=r;//Co
			image->data[k|2]=g;//Cg
		}
		ROTATE3(image->depth[2], image->depth[1], image->depth[0], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y=image->data[k], Co=image->data[k|1], Cg=image->data[k|2];
			
			Y-=Cg>>1;
			Cg+=Y;
			Y-=Co>>1;
			Co+=Y;

			image->data[k  ]=Co;//r
			image->data[k|1]=Cg;//g
			image->data[k|2]=Y;//b
		}
	}
}
void colortransform_YCbCr_R_v1(Image *image, int fwd)
{
	char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			//g+=r;	//r+g		X
			//r-=g>>1;//r-(r+g)/2 = (r-g)/2
			//b-=g>>1;
			//g+=b;

			r-=g;		//diff(r, g)            [ 1      -1      0  ].RGB	Cr
			g+=r>>1;
			b-=g;		//diff(b, av(r, g))     [-1/2    -1/2    1  ].RGB	Cb
			g+=b>>1;	//av(b, av(r, g))       [ 1/4     1/4    1/2].RGB	Y

			image->data[k  ]=g;//Y
			image->data[k|1]=b;//Cb
			image->data[k|2]=r;//Cr
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		ROTATE3(image->depth[2], image->depth[1], image->depth[0], temp);
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y=image->data[k], Cb=image->data[k|1], Cr=image->data[k|2];
			
			Y-=Cb>>1;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			image->data[k  ]=Cr;
			image->data[k|1]=Y;
			image->data[k|2]=Cb;
		}
	}
}
void colortransform_YCbCr_R_v2(Image *image, int fwd)
{
	char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			r-=g;		//Cr =	[1	-1	0].RGB
			g+=r>>1;	//	[1/2	1/2	0]
			b-=g;		//Cb =	[-1/2	-1/2	1]
			g+=(2*b-r)>>3;	//Y  =	[1/4	1/2	1/4]	v2

			image->data[k  ]=g;//Y
			image->data[k|1]=b;//Cb
			image->data[k|2]=r;//Cr
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		ROTATE3(image->depth[2], image->depth[1], image->depth[0], temp);
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y=image->data[k], Cb=image->data[k|1], Cr=image->data[k|2];
			
			Y-=(2*Cb-Cr)>>3;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			image->data[k  ]=Cr;
			image->data[k|1]=Y;
			image->data[k|2]=Cb;
		}
	}
}
void colortransform_YCbCr_R_v3(Image *image, int fwd)
{
	char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			r-=g;		//Cr =	[1	-1	0].RGB
			g+=r>>1;	//	[1/2	1/2	0]
			b-=g;		//Cb =	[-1/2	-1/2	1]
			g+=(2*b+r)>>3;	//Y  =	[1/2	1/4	1/4]	v3

			image->data[k  ]=g;//Y
			image->data[k|1]=b;//Cb
			image->data[k|2]=r;//Cr
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		ROTATE3(image->depth[2], image->depth[1], image->depth[1], temp);
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y=image->data[k], Cb=image->data[k|1], Cr=image->data[k|2];
			
			Y-=(2*Cb+Cr)>>3;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			image->data[k  ]=Cr;
			image->data[k|1]=Y;
			image->data[k|2]=Cb;
		}
	}
}
void colortransform_YCbCr_R_v4(Image *image, int fwd)
{
	char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			r-=g;		//Cr =	[1	-1	0].RGB
			g+=r>>1;	//	[1/2	1/2	0]
			b-=g;		//Cb =	[-1/2	-1/2	1]
			g+=b/3;		//Y  =	[1/3	1/3	1/3]	v4

			image->data[k  ]=g;//Y
			image->data[k|1]=b;//Cb
			image->data[k|2]=r;//Cr
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		ROTATE3(image->depth[2], image->depth[1], image->depth[1], temp);
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y=image->data[k], Cb=image->data[k|1], Cr=image->data[k|2];
			
			Y-=Cb/3;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			image->data[k  ]=Cr;
			image->data[k|1]=Y;
			image->data[k|2]=Cb;
		}
	}
}
void colortransform_YCbCr_R_v5(Image *image, int fwd)
{
	char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			r-=g;		//Cr =	[1	-1	0].RGB
			g+=r>>1;	//	[1/2	1/2	0]
			b-=g;		//Cb =	[-1/2	-1/2	1]
			g+=6*b>>4;	//Y  =	[5/16	5/16	6/16]	v5

			image->data[k  ]=g;//Y
			image->data[k|1]=b;//Cb
			image->data[k|2]=r;//Cr
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y=image->data[k], Cb=image->data[k|1], Cr=image->data[k|2];
			
			Y-=6*Cb>>4;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			image->data[k  ]=Cr;
			image->data[k|1]=Y;
			image->data[k|2]=Cb;
		}
		ROTATE3(image->depth[2], image->depth[1], image->depth[1], temp);
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
	}
}
void colortransform_YCbCr_R_v6(Image *image, int fwd)
{
	char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			r-=g;		//Cr =	[1	-1	0].RGB
			g+=r>>1;	//	[1/2	1/2	0]
			b-=g;		//Cb =	[-1/2	-1/2	1]
			g+=14*b>>5;	//Y  =	[9/32	9/32	14/32]	v6

			image->data[k  ]=g;//Y
			image->data[k|1]=b;//Cb
			image->data[k|2]=r;//Cr
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y=image->data[k], Cb=image->data[k|1], Cr=image->data[k|2];
			
			Y-=14*Cb>>5;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			image->data[k  ]=Cr;
			image->data[k|1]=Y;
			image->data[k|2]=Cb;
		}
		ROTATE3(image->depth[2], image->depth[1], image->depth[1], temp);
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
	}
}
void colortransform_YCbCr_R_v7(Image *image, int fwd)
{
	char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			r-=g;		//Cr =	[1	-1	0].RGB
			g+=r>>1;	//	[1/2	1/2	0]
			b-=g;		//Cb =	[-1/2	-1/2	1]
			g+=(10*b-r)>>5;	//Y  =	[5/16	 6/16	5/16]	v7

			image->data[k  ]=g;//Y
			image->data[k|1]=b;//Cb
			image->data[k|2]=r;//Cr
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y=image->data[k], Cb=image->data[k|1], Cr=image->data[k|2];
			
			Y-=(10*Cb-Cr)>>5;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			image->data[k  ]=Cr;
			image->data[k|1]=Y;
			image->data[k|2]=Cb;
		}
		ROTATE3(image->depth[2], image->depth[1], image->depth[1], temp);
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
	}
}
void colortransform_CrCgCb_R(Image *image, int fwd)
{
	if(fwd)
	{
		for(ptrdiff_t k=(ptrdiff_t)image->iw*image->ih*4-4;k>=0;k-=4)
		{
			//if(k==(ptrdiff_t)image->iw*image->ih*4-4)
			//	printf("");

			int prev_b=k?image->data[k-2]:0, r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			//b-=g;
			//g+=b>>1;
			//g-=r;
			//r+=g>>1;
			//r-=prev_b;
			//prev_b+=r>>1;

			//b-=g;	//diff	LPF
			//g+=b>>1;//av	HPF

			b-=g;
			g-=r;
			r-=prev_b;

			if(k)
				image->data[k-2]=prev_b;
			image->data[k|0]=r;//Cr
			image->data[k|1]=g;//Cg
			image->data[k|2]=b;//Cb
		}
		image->depth[0]+=image->depth[0]<24;
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			//if(k==(ptrdiff_t)image->iw*image->ih*4-4)
			//	printf("");

			int prev_b=k?image->data[k-2]:0, r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			//g-=b>>1;
			//b+=g;

			r+=prev_b;
			g+=r;
			b+=g;

			image->data[k|0]=r;
			image->data[k|1]=g;
			image->data[k|2]=b;
		}
		image->depth[0]-=image->depth[0]>image->src_depth[1];
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
	}
}
void colortransform_Pei09(Image *image, int fwd)//Pei09 RCT
{
	char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			b-=(87*r+169*g+128)>>8;	//Cb = [-87/256  -169/256  1]
			r-=g;			//Cr = [1  -1  0].RGB
			g+=(86*r+29*b+128)>>8;	//Y  = [19493/65536  38619/65536  29/256]	g+86/256*(r-g)+29/256*(b-87/256*r-169/256*g) = 19493/65536*r + 38619/65536*g + 29/256*b

			image->data[k  ]=g;//Y
			image->data[k|1]=b;//Cb
			image->data[k|2]=r;//Cr
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		ROTATE3(image->depth[2], image->depth[1], image->depth[0], temp);
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y=image->data[k], Cb=image->data[k|1], Cr=image->data[k|2];
			
			Y-=(86*Cr+29*Cb+128)>>8;
			Cr+=Y;
			Cb+=(87*Cr+169*Y+128)>>8;

			image->data[k  ]=Cr;
			image->data[k|1]=Y;
			image->data[k|2]=Cb;
		}
	}
}
void colortransform_JPEG2000(Image *image, int fwd)
{
	char temp;
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
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		ROTATE3(image->depth[2], image->depth[1], image->depth[0], temp);
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
void colortransform_subtractgreen(Image *image, int fwd)
{
	char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			r-=g;
			b-=g;

			image->data[k  ]=g;//Y
			image->data[k|1]=b;//Cb
			image->data[k|2]=r;//Cr
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		ROTATE3(image->depth[2], image->depth[1], image->depth[0], temp);
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y=image->data[k], Cb=image->data[k|1], Cr=image->data[k|2];
			
			Cb+=Y;
			Cr+=Y;

			image->data[k  ]=Cr;
			image->data[k|1]=Y;
			image->data[k|2]=Cb;
		}
	}
}
void colortransform_J2K2(Image *image, int fwd)//apply twice
{
	int sh[]=
	{
		32-image->depth[0],
		32-image->depth[1],
		32-image->depth[2],
		32-image->depth[3],
	};
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int
				r=image->data[k|0],
				g=image->data[k|1],
				b=image->data[k|2];
			
			r-=g>>1;	r=r<<sh[0]>>sh[0];
			b-=g>>1;	b=b<<sh[2]>>sh[2];
			g+=(r+b)>>4;	g=g<<sh[1]>>sh[1];

			image->data[k|0]=r;//Cr
			image->data[k|1]=g;//Y
			image->data[k|2]=b;//Cb
		}
	}
	else
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int
				r=image->data[k|0],
				g=image->data[k|1],
				b=image->data[k|2];
			
			g-=(r+b)>>4;	g=g<<sh[1]>>sh[1];
			b+=g>>1;	b=b<<sh[2]>>sh[2];
			r+=g>>1;	r=r<<sh[0]>>sh[0];

			image->data[k|0]=r;
			image->data[k|1]=g;
			image->data[k|2]=b;
		}
	}
}
void rct_yrgb_v1(Image *image, int fwd)
{
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];

			//int base=r>>1;
			//r-=base;

			int base=r;//C={r, g, b}, base=C[argmin(abs(C))]
			if(abs(base)>abs(g))base=g;
			if(abs(base)>abs(b))base=b;
			r-=base;
			g-=base;
			b-=base;

			//int base=image->data[k|3];
			//r>>=1;
			//g>>=1;
			//b>>=1;

			image->data[k  ]=r;//Y
			image->data[k|1]=g;//Cb
			image->data[k|2]=b;//Cr
			image->data[k|3]=base;//Cg
		}
		{
			int d0=image->depth[0], d1=image->depth[1], d2=image->depth[2];
			int dbase=MINVAR(d0, d1);
			dbase=MINVAR(dbase, d2);
			d0+=d0<24;
			d1+=d1<24;
			d2+=d2<24;
			image->depth[0]=(char)dbase;
			image->depth[1]=(char)d2;
			image->depth[2]=(char)d0;
			image->depth[3]=(char)d1;
		}
	}
	else
	{
		image->depth[0]-=image->depth[0]>image->src_depth[0];
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		image->depth[3]=0;
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int base=image->data[k], b=image->data[k|1], r=image->data[k|2], g=image->data[k|3];
			
			r+=base;
			b+=base;
			g+=base;

			image->data[k  ]=r;
			image->data[k|1]=g;
			image->data[k|2]=b;
			image->data[k|3]=0xFF;
		}
	}
}
void rct_yrgb_v2(Image *image, int fwd)
{
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];

			int base=(r+g+b+1)/3;
			//int base=(r+g+g+b+2)>>2;
			r-=base;
			g-=base;
			b-=base;

			image->data[k  ]=base;//Y
			image->data[k|1]=b;//Cb
			image->data[k|2]=r;//Cr
			image->data[k|3]=g;//Cg
		}
		{
			int d0=image->depth[0], d1=image->depth[1], d2=image->depth[2];
			int dbase=MINVAR(d0, d1);
			dbase=MINVAR(dbase, d2);
			d0+=d0<24;
			d1+=d1<24;
			d2+=d2<24;
			image->depth[0]=(char)dbase;
			image->depth[1]=(char)d2;
			image->depth[2]=(char)d0;
			image->depth[3]=(char)d1;
		}
	}
	else
	{
		image->depth[0]-=image->depth[0]>image->src_depth[0];
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		image->depth[3]=0;
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int base=image->data[k], b=image->data[k|1], r=image->data[k|2], g=image->data[k|3];
			
			r+=base;
			b+=base;
			g+=base;

			image->data[k  ]=r;
			image->data[k|1]=g;
			image->data[k|2]=b;
			image->data[k|3]=0xFF;
		}
	}
}
void ct_cmyk_fwd(Image *image)//irreversible
{
	int nlevels[]=
	{
		1<<image->depth[0],
		1<<image->depth[1],
		1<<image->depth[2],
	};
	for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
	{
		double
			r=(double)(image->data[k|0]+(nlevels[0]>>1))/(nlevels[0]-1),//[0, 1]
			g=(double)(image->data[k|1]+(nlevels[1]>>1))/(nlevels[1]-1),
			b=(double)(image->data[k|2]+(nlevels[2]>>1))/(nlevels[2]-1);

		double luma=MAXVAR(r, g);
		luma=MAXVAR(luma, b);
		if(fabs(luma)<1e-9)
		{
			r=(luma-r)/luma;
			g=(luma-g)/luma;
			b=(luma-b)/luma;
			luma=1-luma;
		}
		else
		{
			luma=1;
			r=g=b=0;
		}

		if((unsigned)(int)(r*255)>255||(unsigned)(int)(g*255)>255||(unsigned)(int)(b*255)>255||(unsigned)(int)(luma*255)>255)
			LOG_ERROR("");
		image->data[k|0]=(unsigned char)(255*r)-128;//C
		image->data[k|1]=(unsigned char)(255*g)-128;//M
		image->data[k|2]=(unsigned char)(255*b)-128;//Y
		image->data[k|3]=(unsigned char)(255*luma)-128;//K
	}
	image->depth[0]=8;
	image->depth[1]=8;
	image->depth[2]=8;
	image->depth[3]=8;
}

void colortransform_lossy_YCbCr(Image *image, int fwd)//for demonstration purposes only
{
	if(fwd)
	{
		for(ptrdiff_t k=0, res=image->iw*image->ih;k<res;++k)
		{
			int
				r=image->data[k<<2|0],
				g=image->data[k<<2|1],
				b=image->data[k<<2|2];

			double
				Y=0.299*r+0.587*g+0.114*b,
				Cb=-0.168736*r-0.331264*g+0.5*b,
				Cr=+0.5*r-0.418688*g-0.081312*b;

			image->data[k<<2|0]=(int)Y;
			image->data[k<<2|1]=(int)Cb;
			image->data[k<<2|2]=(int)Cr;
		}
		//image->depth[1]+=image->depth[1]<24;
		//image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		//image->depth[1]-=image->depth[1]>image->src_depth[1];
		//image->depth[2]-=image->depth[2]>image->src_depth[2];
		for(ptrdiff_t k=0, res=image->iw*image->ih;k<res;++k)
		{
			int
				Y =image->data[k<<2|0],
				Cb=image->data[k<<2|1],
				Cr=image->data[k<<2|2];

			double
				r=Y+1.402*Cr,
				g=Y-0.344136*Cb-0.714136*Cr,
				b=Y+1.772*Cb;

			image->data[k<<2|0]=(int)r;
			image->data[k<<2|1]=(int)g;
			image->data[k<<2|2]=(int)b;
		}
	}
}
void colortransform_lossy_XYB(Image *image, int fwd)
{
	const double bias=0.00379307325527544933;
	double cbrt_bias=cbrt(bias);
	int nlevels[3];
	char temp;
	if(fwd)
	{
		nlevels[0]=1<<image->depth[0];
		nlevels[1]=1<<image->depth[1];
		nlevels[2]=1<<image->depth[2];
		for(ptrdiff_t k=0, res=(ptrdiff_t)image->iw*image->ih;k<res;++k)
		{
			double
				r=(double)(image->data[k<<2|0]+(nlevels[0]>>1))/(nlevels[0]-1),//input: RGB in [0, 1]
				g=(double)(image->data[k<<2|1]+(nlevels[1]>>1))/(nlevels[1]-1),
				b=(double)(image->data[k<<2|2]+(nlevels[2]>>1))/(nlevels[2]-1);

			double
				L=0.3*r+0.622*g+0.078*b,
				M=0.23*r+0.692*g+0.078*b,
				S=0.24342268924547819*r+0.20476744424496821*g+0.55180986650955360*b;
			L=cbrt(L+bias)-cbrt_bias;
			M=cbrt(M+bias)-cbrt_bias;
			S=cbrt(S+bias)-cbrt_bias;
			r=(L-M)*0.5;	//X
			g=(L+M)*0.5;	//Y
			b=S-g;		//B-Y

			r*=(int)(4096.*((long long)nlevels[0]<<(nlevels[0]<0x1000000))/512);//customparam_ct[0]*1000
			g*=(int)(144.*nlevels[1]/256);
			b*=(int)(1024.*((long long)nlevels[2]<<(nlevels[2]<0x1000000))/512);
			r=CLAMP(-nlevels[0], r, nlevels[0]-1);
			g=CLAMP(-(nlevels[1]>>1), g, (nlevels[1]>>1)-1);
			b=CLAMP(-nlevels[2], b, nlevels[2]-1);

			image->data[k<<2|0]=(int)g;//Y
			image->data[k<<2|1]=(int)b;//Cb
			image->data[k<<2|2]=(int)r;//Cr
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		ROTATE3(image->depth[2], image->depth[1], image->depth[0], temp);
		nlevels[0]=1<<image->depth[0];
		nlevels[1]=1<<image->depth[1];
		nlevels[2]=1<<image->depth[2];
		for(ptrdiff_t k=0, res=(ptrdiff_t)image->iw*image->ih;k<res;++k)
		{
			double
				Y =(double)image->data[k<<2|0]/(144.*nlevels[1]/256),
				Cb=(double)image->data[k<<2|1]/(1024.*((long long)nlevels[2]<<(nlevels[2]<0x1000000))/512),
				Cr=(double)image->data[k<<2|2]/(4096.*((long long)nlevels[0]<<(nlevels[0]<0x1000000))/512);

			double
				L=Y+Cr,
				M=Y-Cr,//g is the sum
				S=Cb+Y;
			L+=cbrt_bias, L*=L*L, L-=bias;
			M+=cbrt_bias, M*=M*M, M-=bias;
			S+=cbrt_bias, S*=S*S, S-=bias;
			{
				double
					r=11.03160*L-9.86694*M-0.164623*S,
					g=-3.25415*L+4.41877*M-0.164623*S,
					b=-3.65885*L+2.71292*M+1.945930*S;

				image->data[k<<2|0]=(int)((r-0.5)*nlevels[0]);//output: RGB in [0, 1]
				image->data[k<<2|1]=(int)((g-0.5)*nlevels[1]);
				image->data[k<<2|2]=(int)((b-0.5)*nlevels[2]);
			}
		}
	}
}
void colortransform_lossy_matrix(Image *image, int fwd)//for demonstration purposes only
{
	if(fwd)
	{
		for(ptrdiff_t k=0, res=(ptrdiff_t)image->iw*image->ih;k<res;++k)
		{
			int
				r=image->data[k<<2|0],
				g=image->data[k<<2|1],
				b=image->data[k<<2|2];

			double
				Y=0.299*r+0.587*g+0.114*b,
				Cb=(-0.168736*r-0.331264*g+0.5*b)/2,
				Cr=(0.5*r-0.418688*g-0.081312*b)/2;

			image->data[k<<2|0]=(int)Y;
			image->data[k<<2|1]=(int)Cb;
			image->data[k<<2|2]=(int)Cr;
		}
		//image->depth[1]+=image->depth[1]<24;
		//image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		//image->depth[1]-=image->depth[1]>image->src_depth[1];
		//image->depth[2]-=image->depth[2]>image->src_depth[2];
		for(ptrdiff_t k=0, res=image->iw*image->ih;k<res;++k)
		{
			int
				Y =image->data[k<<2|0],
				Cb=image->data[k<<2|1],
				Cr=image->data[k<<2|2];

			double
				r=Y+(1.402*Cr)*2,
				g=Y+(-0.344136*Cb-0.714136*Cr)*2,
				b=Y+(1.772*Cb)*2;

			image->data[k<<2|0]=(int)r;
			image->data[k<<2|1]=(int)g;
			image->data[k<<2|2]=(int)b;
		}
	}
}

short rct_custom_params[RCT_CUSTOM_NPARAMS]={0};
void rct_custom_unpackpermutation(short p, unsigned char *permutation)
{
	int temp=0;
	switch(p)
	{
	case 0:temp=0x020100;break;//rgb
	case 1:temp=0x020001;break;//grb
	case 2:temp=0x000201;break;//gbr
	case 3:temp=0x010200;break;//rbg
	case 4:temp=0x000102;break;//bgr
	case 5:temp=0x010002;break;//brg
	default:
		LOG_ERROR("Invalid RGB permutation");
		return;
	}
	memcpy(permutation, &temp, sizeof(char[3]));
}
void rct_custom(Image *image, int fwd, const short *params)//4 params	fixed 15.16
{
	int temp;
	unsigned char permutation[4]={0};
	rct_custom_unpackpermutation(params[8], permutation);
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int comp[]=
			{
				image->data[k|permutation[0]],
				image->data[k|permutation[1]],
				image->data[k|permutation[2]],
			};
			temp=params[0]*comp[1]+params[1]*comp[2], comp[0]+=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			temp=params[2]*comp[0]+params[3]*comp[2], comp[1]+=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			temp=params[4]*comp[0]+params[5]*comp[1], comp[2]+=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			temp=params[6]*comp[0]+params[7]*comp[2], comp[1]+=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			//temp=params[0]*comp[1]+params[1]*comp[2], comp[0]+=(temp>>RCT_CUSTOM_PARAMBITS)+(temp<0);
			//temp=params[2]*comp[0]+params[3]*comp[2], comp[1]+=(temp>>RCT_CUSTOM_PARAMBITS)+(temp<0);
			//temp=params[4]*comp[0]+params[5]*comp[1], comp[2]+=(temp>>RCT_CUSTOM_PARAMBITS)+(temp<0);
			//temp=params[6]*comp[0]+params[7]*comp[2], comp[1]+=(temp>>RCT_CUSTOM_PARAMBITS)+(temp<0);
			//image->data[k|0]=comp[1];//Y
			//image->data[k|1]=comp[2];//Cb
			//image->data[k|2]=comp[0];//Cr
			memcpy(image->data+k, comp, sizeof(comp));
		}
		image->depth[0]+=image->depth[0]<24;
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
		//ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		//image->depth[1]+=image->depth[1]<24;
		//image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		image->depth[0]-=image->depth[0]>image->src_depth[0];
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		//image->depth[1]-=image->depth[1]>image->src_depth[1];
		//image->depth[2]-=image->depth[2]>image->src_depth[2];
		//ROTATE3(image->depth[2], image->depth[1], image->depth[0], temp);
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int comp[3];
			//int comp[]=
			//{
			//	image->data[k|2],
			//	image->data[k|0],
			//	image->data[k|1],
			//};
			memcpy(comp, image->data+k, sizeof(comp));
			temp=params[6]*comp[0]+params[7]*comp[2], comp[1]-=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			temp=params[4]*comp[0]+params[5]*comp[1], comp[2]-=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			temp=params[2]*comp[0]+params[3]*comp[2], comp[1]-=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			temp=params[0]*comp[1]+params[1]*comp[2], comp[0]-=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			image->data[k|permutation[0]]=comp[0];
			image->data[k|permutation[1]]=comp[1];
			image->data[k|permutation[2]]=comp[2];
		}
	}
#if 0
	char temp;
	if(fwd)
	{
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=image->data[k], g=image->data[k|1], b=image->data[k|2];
			
			r-=g;
			g+=(r*params[0]+b*params[1])>>16;
			b-=g;
			g+=(r*params[2]+b*params[3])>>16;

			image->data[k  ]=g;//Y
			image->data[k|1]=b;//Cb
			image->data[k|2]=r;//Cr
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		ROTATE3(image->depth[2], image->depth[1], image->depth[0], temp);
		for(ptrdiff_t k=0, len=(ptrdiff_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y=image->data[k], Cb=image->data[k|1], Cr=image->data[k|2];
			
			Y-=(Cr*params[2]+Cb*params[3])>>16;
			Cb+=Y;
			Y-=(Cr*params[0]+Cb*params[1])>>16;
			Cr+=Y;

			image->data[k  ]=Cr;
			image->data[k|1]=Y;
			image->data[k|2]=Cb;
		}
	}
#endif
}
static void rct_custom_calcloss(Image const *src, Image *dst, int *hist, double *loss)
{
	int depths[]=
	{
		src->depth[1],
		src->depth[2]+1,
		src->depth[0]+1,
	};
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	memcpy(dst->data, src->data, res*sizeof(int[4]));

	apply_selected_transforms(dst, 0);
	//rct_custom(dst, 1, params);
	//pred_grad2(dst, 1);
	//pred_clampedgrad(dst, 1, 1);
	
	for(int kc=0;kc<3;++kc)
	{
		calc_histogram(dst->data, dst->iw, dst->ih, kc, 0, dst->iw, 0, dst->ih, depths[kc], hist, 0);
		{
			double entropy=calc_entropy(hist, 1<<depths[kc], (int)res);
			loss[kc]=entropy/src->src_depth[(kc+1)%3];
		}
	}
	loss[3]=(loss[0]+loss[1]+loss[2])/3;
}
#define RCT_CUSTOM_NITER 200
#define RCT_CUSTOM_DELTAGROUP 2
void rct_custom_optimize(Image const *image, short *params)
{
	static int call_idx=0;

	int inflation[]={1, 0, 1};//original image is expected, RCT permutes depths {r+1, g, b+1} -> {Y:g:+0, Cb:b:+1, Cr:r:+1}
	int maxdepth=calc_maxdepth(image, inflation);
	int maxlevels=1<<maxdepth;
	int *hist;
	double loss_bestsofar[4], loss_prev[4], loss_curr[4];
	short params2[RCT_CUSTOM_NPARAMS];
	int shakethreshold=RCT_CUSTOM_NPARAMS;
	Image *im2;

	++call_idx;
	if(call_idx==1)
		DisableProcessWindowsGhosting();

	im2=0;
	image_copy(&im2, image);
	hist=(int*)malloc(maxlevels*sizeof(int));
	if(!im2||!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}

	memcpy(params2, params, sizeof(params2));

#define CALC_LOSS(L) rct_custom_calcloss(image, im2, hist, L)
#ifndef _DEBUG
	srand((unsigned)__rdtsc());
#endif
	CALC_LOSS(loss_bestsofar);
	memcpy(loss_prev, loss_bestsofar, sizeof(loss_prev));

	for(int it=0, watchdog=0;it<RCT_CUSTOM_NITER;++it)
	{
		int idx[RCT_CUSTOM_DELTAGROUP]={0}, stuck=0;
		short params_original_selected[RCT_CUSTOM_DELTAGROUP]={0};

		if(watchdog>=shakethreshold)//bump if stuck
		{
			memcpy(params, params2, sizeof(params2));
			for(int k=0;k<RCT_CUSTOM_NPARAMS-1;++k)
			{
				//if(k==8)
				//	params[8]=rand()%6;
				//else
					params[k]+=(rand()&63)-32;
			}
			watchdog=0;
			stuck=1;
		}
		else
		{
			for(int k=0;k<RCT_CUSTOM_DELTAGROUP;++k)//increment params
			{
				//int inc=0;
				idx[k]=rand()%RCT_CUSTOM_NPARAMS;
				//while(!(inc=rand()-(RAND_MAX>>1)));//reject zero delta
		
				params_original_selected[k]=params[idx[k]];
				//params[idx[k]]+=(inc<<2)/RAND_MAX;
				if(idx[k]==8)
					params[8]=rand()%6;
				//{
				//	params[8]+=((rand()&1)<<1)-1;
				//	MODVAR(params[8], params[8], 6);
				//}
				else
					params[idx[k]]+=(rand()&63)-32;
			}
		}

		CALC_LOSS(loss_curr);
		
		if(loss_prev[3]<loss_curr[3])//revert if worse
		{
			if(stuck)//a bad branch may surpass the local minimum
				memcpy(loss_prev, loss_curr, sizeof(loss_prev));
			else
			{
				memcpy(loss_curr, loss_prev, sizeof(loss_prev));
				for(int k=0;k<RCT_CUSTOM_DELTAGROUP;++k)
					params[idx[k]]=params_original_selected[k];
			}
			++watchdog;
		}
		else//save if better
		{
			if(loss_curr[3]<loss_bestsofar[3])//publish if record best
			{
				memcpy(params2, params, sizeof(params2));
				memcpy(loss_bestsofar, loss_curr, sizeof(loss_bestsofar));
				--it;//again
			}
			memcpy(loss_prev, loss_curr, sizeof(loss_prev));
			watchdog=0;
		}

		if(loud_transforms)
			set_window_title(
				"%d %4d/%4d,%d/%d: %lf%% RGB %lf %lf %lf%s",
				call_idx,
				it+1,
				RCT_CUSTOM_NITER,
				watchdog,
				shakethreshold,
				100.*loss_bestsofar[3],
				100.*loss_bestsofar[0],
				100.*loss_bestsofar[1],
				100.*loss_bestsofar[2],
				it+1<RCT_CUSTOM_NITER?"...":" Done."
			);

		//preview
#if 1
		{
			ch_entropy[0]=(float)(loss_bestsofar[0]*image->src_depth[0]);
			ch_entropy[1]=(float)(loss_bestsofar[1]*image->src_depth[1]);
			ch_entropy[2]=(float)(loss_bestsofar[2]*image->src_depth[2]);
			ch_entropy[3]=0;
			//ch_entropy[3]=(float)(loss_bestsofar[3]*image->src_depth[1]);//X
			//unsigned char *ptr;
			//addhalf(temp, iw, ih, 3, 4);
			//SWAPVAR(image, temp, ptr);
			io_render();
			//SWAPVAR(image, temp, ptr);
		}
#endif
	}
	memcpy(params, params2, sizeof(params2));
#undef  CALC_LOSS
	free(hist);
	free(im2);
}
void rct_custom_getmatrix(double *matrix, int fwd)
{
	unsigned char per[4]={0};
	rct_custom_unpackpermutation(rct_custom_params[8], per);
	memset(matrix, 0, sizeof(double[9]));
	matrix[0]=1;
	matrix[4]=1;
	matrix[8]=1;
	for(int k=0;k<3;++k)
	{
		double vrtx[]=
		{
			matrix[k+3*(fwd?per[0]:0)],
			matrix[k+3*(fwd?per[1]:1)],
			matrix[k+3*(fwd?per[2]:2)],
			//matrix[k+3*0],
			//matrix[k+3*1],
			//matrix[k+3*2],
		};
		if(fwd)
		{
			vrtx[0]+=(rct_custom_params[0]*vrtx[1]+rct_custom_params[1]*vrtx[2])*(1/4096.);
			vrtx[1]+=(rct_custom_params[2]*vrtx[0]+rct_custom_params[3]*vrtx[2])*(1/4096.);
			vrtx[2]+=(rct_custom_params[4]*vrtx[0]+rct_custom_params[5]*vrtx[1])*(1/4096.);
			vrtx[1]+=(rct_custom_params[6]*vrtx[0]+rct_custom_params[7]*vrtx[2])*(1/4096.);
		}
		else
		{
			vrtx[1]-=(rct_custom_params[6]*vrtx[0]+rct_custom_params[7]*vrtx[2])*(1/4096.);
			vrtx[2]-=(rct_custom_params[4]*vrtx[0]+rct_custom_params[5]*vrtx[1])*(1/4096.);
			vrtx[1]-=(rct_custom_params[2]*vrtx[0]+rct_custom_params[3]*vrtx[2])*(1/4096.);
			vrtx[0]-=(rct_custom_params[0]*vrtx[1]+rct_custom_params[1]*vrtx[2])*(1/4096.);
		}
		matrix[k+3*(fwd?0:per[0])]=vrtx[0];
		matrix[k+3*(fwd?1:per[1])]=vrtx[1];
		matrix[k+3*(fwd?2:per[2])]=vrtx[2];
		//matrix[k+3*0]=vrtx[0];
		//matrix[k+3*1]=vrtx[1];
		//matrix[k+3*2]=vrtx[2];
	}
}



//void addhalf(unsigned char *buf, int iw, int ih, int nch, int bytestride)
//{
//	for(int kp=0, len=iw*ih*bytestride;kp<len;kp+=bytestride)
//	{
//		for(int kc=0;kc<nch;++kc)
//			buf[kp+kc]+=128;
//	}
//}
static void c_mul_add(double *dst, const double *a, const double *b)
{
	double
		r=a[0]*b[0]-a[1]*b[1],
		i=a[0]*b[1]+a[1]*b[0];
	dst[0]+=r, dst[1]+=i;
}
static void c_mul_sub(double *dst, const double *a, const double *b)
{
	double
		r=a[0]*b[0]-a[1]*b[1],
		i=a[0]*b[1]+a[1]*b[0];
	dst[0]-=r, dst[1]-=i;
}
static double c_abs2(const double *z)
{
	return z[0]*z[0]+z[1]*z[1];
}
static void c_div(double *dst, double const *a, double const *b)
{
	double invabsb2=1/c_abs2(b);
	double
		r=(a[0]*b[0]+a[1]*b[1])*invabsb2,
		i=(a[1]*b[0]-a[0]*b[1])*invabsb2;
	dst[0]=r, dst[1]=i;
}
static void c_exp(double *dst, double const *x)
{
	double m=exp(x[0]);
	dst[0]=m*cos(x[1]);
	dst[1]=m*sin(x[1]);
}
static void c_ln(double *dst, double const *x)
{
	double
		r=log(sqrt(c_abs2(x))),
		i=atan2(x[1], x[0]);
	dst[0]=r, dst[1]=i;
}
static void c_cbrt(double *dst, const double *x)//sqrt(x)=exp(0.5lnx)
{
	if(fabs(x[0])<1e-9||fabs(x[1])<1e-9)
	{
		double temp[2];

		c_ln(temp, x);
		temp[0]*=1./3, temp[1]*=1./3;
		c_exp(dst, temp);
	}
	else
		dst[0]=x[0], dst[1]=x[1];
}
static void impl_solve_cubic(const double *coeffs, double *roots)//finds r[0], r[1], & r[2], the solutions of x^3 + c[2]x^2 + c[1]x + c[0] = 0
{
	//https://math.stackexchange.com/questions/15865/why-not-write-the-solutions-of-a-cubic-this-way/18873#18873
	double p=coeffs[2], q=coeffs[1], r=coeffs[0],
		p2, p3, q2, pq, A[2], B[2];
	double sqrt27=3*sqrt(3), inv3cbrt2=1./(3*cbrt(2)), ninth=1./9;
	double cm[]={-0.5, -sqrt(3)*0.5}, cp[]={-0.5, sqrt(3)*0.5};

	if(fabs(p)<1e-9&&fabs(q)<1e-9)
	{
		roots[4]=roots[2]=roots[0]=-cbrt(r);
		roots[5]=roots[3]=roots[1]=0;
		return;
	}

	p2=p*p;
	p3=p2*p;
	q2=q*q;
	pq=p*q;

	A[0]=(4*q-p2)*q2+(4*p3-18*pq+27*r)*r;
	if(A[0]<0)//complex sqrt
		A[1]=sqrt(fabs(A[0])), A[0]=0;
	else
		A[0]=sqrt(A[0]), A[1]=0;
	A[0]*=sqrt27;
	A[1]*=sqrt27;
	A[0]-=27*r;
	A[0]+=9*pq-2*p3;
	c_cbrt(A, A);
	A[0]*=inv3cbrt2;
	A[1]*=inv3cbrt2;
	B[0]=3*q-p2;
	B[1]=0;
	c_div(B, B, A);
	B[0]*=ninth;
	B[1]*=ninth;
	roots[4]=roots[2]=roots[0]=p*(-1./3);
	roots[5]=roots[3]=roots[1]=0;
	roots[0]+=A[0]-B[0];
	roots[1]+=A[1]-B[1];
	c_mul_add(roots+2, A, cm);
	c_mul_sub(roots+2, B, cp);
	c_mul_add(roots+4, A, cp);
	c_mul_sub(roots+4, B, cm);
}
static void impl_rref(double *m, short dx, short dy)
{
#ifdef _DEBUG
	double pivot;
#endif
	double coeff;
	int mindim=dx<dy?dx:dy, it, ky, kx, npivots, kpivot;

	for(it=0, npivots=0;it<mindim;++it)//iteration
	{
		kpivot=-1;
		for(ky=npivots;ky<dy;++ky)//find pivot
		{
			if(fabs(m[dx*ky+it])>1e-10)
			{
#ifdef _DEBUG
				pivot=m[dx*ky+it];
#endif
				kpivot=ky;
				++npivots;
				break;
			}
		}
		if(kpivot==-1)
			continue;
		if(kpivot>npivots-1)
		{
			for(kx=0;kx<dx;++kx)//swap rows
				coeff=m[dx*kpivot+kx], m[dx*kpivot+kx]=m[dx*(npivots-1)+kx], m[dx*(npivots-1)+kx]=coeff;
			kpivot=npivots-1;
		}
		for(ky=0;ky<dy;++ky)
		{
			if(ky==kpivot)//normalize pivot row
			{
				coeff=1/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*kpivot+kx]*=coeff;
			}
			else//subtract pivot row from all other rows
			{
				coeff=m[dx*ky+it]/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*ky+kx]-=coeff*m[dx*kpivot+kx];
			}
			//print_matrix_debug(m, dx, dy);//
		}
	}
}
static int impl_nullspace(double *M, int dx, int dy, double *solution, int sstart, char *dep_flags, short *row_idx)
{
	//M is rref'ed
	//solution allocated size dy*dy, pre-memset solution to zero
	//p_flags & row_idx both size dx
	//returns number of vectors in solution
	int kx, kxdst, ky, keq, kfree, idx, idx2, nvec;
#ifdef DEBUG_NULLSPACE
	printf("Before RREF:\n");
	print_matrix_debug(M, dx, dy);
#endif
	impl_rref(M, (short)dx, (short)dy);//first, get rref(M)
#ifdef DEBUG_NULLSPACE
	printf("After RREF:\n");
	print_matrix_debug(M, dx, dy);
#endif
	memset(dep_flags, 0, dx);
	memset(row_idx, 0, dx*sizeof(short));
	for(ky=0;ky<dy;++ky)//find pivots (dependent variables)
	{
		for(kx=ky;kx<dx;++kx)//find first nonzero element
		{
			idx=dx*ky+kx;
			if(fabs(M[idx])>1e-10)
				break;
		}
		if(kx<dx)//if found
			dep_flags[kx]=1, row_idx[ky]=(short)kx;
		else
			break;
	}
	nvec=dx-ky;
	for(ky=0, keq=0, kfree=0;ky<dx;++ky)
	{
		if(dep_flags[ky])//pivot, dependent variable
		{
			for(kx=0, kxdst=0;kx<dx;++kx)
			{
				if(dep_flags[kx])
					continue;
				idx=dx*ky+kxdst, idx2=dx*keq+kx;
				if(sstart+idx>=9)
					LOG_ERROR("");
				solution[sstart+idx]=-M[idx2];
				++kxdst;
			}
			++keq;
		}
		else//free variable
		{
			idx=dx*ky;
			if(sstart+idx+nvec>9)
				LOG_ERROR("");
			memset(solution+sstart+idx, 0, nvec*sizeof(double));
			if(sstart+idx+kfree>=9)
				LOG_ERROR("");
			solution[sstart+idx+kfree]=1;
			++kfree;
		}
#ifdef DEBUG_NULLSPACE
		//printf("Nullspace row %d:\n", ky);
		print_matrix_debug(solution+dx*ky, nvec, 1);//
#endif
	}
	return nvec;
}
static int impl_egvec(double const *M, int n, const double *lambdas, double *S)
{
	int kv, kx, nvec, size=n*n, again;
	double *temp=(double*)malloc(size*sizeof(double));
	char *dep_flags=(char*)malloc(n);
	short *row_idx=(short*)malloc(n*sizeof(short));
	if(!temp||!dep_flags||!row_idx)
	{
		LOG_ERROR("Alloc error");
		return 0;
	}
	memset(S, 0, size*sizeof(double));
	for(kv=0, nvec=0;kv<n&&nvec<n;++kv)//get eigenvectors
	{
		again=0;
		for(kx=0;kx<kv;++kx)//check for repeated eigenvalues
		{
			if(fabs(lambdas[kx]-lambdas[kv])<1e-10)
			{
				again=1;
				break;
			}
		}
		if(again)
			continue;
		memcpy(temp, M, size*sizeof(double));
		for(kx=0;kx<n;++kx)
			temp[(n+1)*kx]-=lambdas[kv];

		nvec+=impl_nullspace(temp, n, n, S, nvec, dep_flags, row_idx);
		//if(nvec>6)
		//	LOG_ERROR("OOB");
		//print_matrix_debug(S, n, n);//
		//for(ky=0;ky<n;++ky)
		//{
		//	for(kx=ky;kx<n;++kx)
		//		if(fabs(temp[n*ky+kx].r)>1e-10||fabs(temp[n*ky+kx].i)>1e-10)
		//			break;
		//}
	}
	free(dep_flags), free(row_idx);
	free(temp);
	return nvec;
}
static void impl_egval3(double const *M, double *lambdas)//finds the complex eigenvalues of a 3x3 matrix
{
	double C[3];

	C[2]=-(M[0]+M[4]+M[8]);//-tr(M)
	C[1]=M[4]*M[8]-M[5]*M[7] + M[0]*M[8]-M[2]*M[6] + M[0]*M[4]-M[1]*M[3];//cof0+cof4+cof8
	C[0]=-(M[0]*(M[4]*M[8]-M[5]*M[7]) - M[1]*(M[3]*M[8]-M[5]*M[6]) + M[2]*(M[3]*M[7]-M[4]*M[6]));//-det(M)

	impl_solve_cubic(C, lambdas);
}

#define ARCT_REACH 10
//#define ARCT_REACH 1
#define ARCT_NNB (2*(ARCT_REACH+1)*ARCT_REACH)
#define CROSS(DST, A, B)\
	DST[0]=A[1]*B[2]-A[2]*B[1],\
	DST[1]=A[2]*B[0]-A[0]*B[2],\
	DST[2]=A[0]*B[1]-A[1]*B[0]
#if 0
static void default_RCT(const int *src, int *dst, int fwd)
{
	int r=src[0], g=src[1], b=src[2];
	if(fwd)
	{
		r-=g;
		b-=g;
		g+=(r+b)>>2;
	}
	else
	{
		g-=(r+b)>>2;
		b+=g;
		r+=g;
	}
	dst[0]=r;
	dst[1]=g;
	dst[2]=b;
}
#endif
static void rct_adaptive_block(int *src, int iw, int ih, int x, int y, int blocksize, double *axes)
{
	double
		mean[3]={0}, cov[9]={0},
		lambdas[6],
		axis1[3], axis2[3], axis3[3],
		ev[9],
		d;
	int res, nv;
	char d1=0, d2=1, d3=2, temp;

	res=iw*ih;
	for(int ky=0;ky<blocksize&&y+ky<ih;++ky)
	{
		for(int kx=0;kx<blocksize&&x+kx<iw;++kx)
		{
			int idx=(iw*(y+ky)+x+kx)<<2;
			int r=src[idx|0], g=src[idx|1], b=src[idx|2];
			mean[0]+=r;
			mean[1]+=g;
			mean[2]+=b;
		}
	}
	mean[0]/=res;
	mean[1]/=res;
	mean[2]/=res;
	for(int ky=0;ky<blocksize&&y+ky<ih;++ky)
	{
		for(int kx=0;kx<blocksize&&x+kx<iw;++kx)
		{
			int idx=(iw*(y+ky)+x+kx)<<2;
			double
				r=(src[idx|0]-mean[0]),
				g=(src[idx|1]-mean[1]),
				b=(src[idx|2]-mean[2]);
			double rr=r*r, gg=g*g, bb=b*b, rg=r*g, gb=g*b, br=b*r;
			cov[0]+=rr;
			cov[1]+=rg;
			cov[2]+=br;
			//cov[3]+=rg;//cov[1]
			cov[4]+=gg;
			cov[5]+=gb;
			//cov[6]+=br;//cov[2]
			//cov[7]+=gb;//cov[5]
			cov[8]+=bb;
		}
	}
	cov[3]=cov[1];
	cov[6]=cov[2];
	cov[7]=cov[5];
	for(int k=0;k<9;++k)
		cov[k]/=res-1;
	impl_egval3(cov, lambdas);//get eigenvalues
	lambdas[0]=sqrt(lambdas[0]*lambdas[0]+lambdas[1]*lambdas[1]);
	lambdas[1]=sqrt(lambdas[2]*lambdas[2]+lambdas[3]*lambdas[3]);
	lambdas[2]=sqrt(lambdas[4]*lambdas[4]+lambdas[5]*lambdas[5]);
	nv=impl_egvec(cov, 3, lambdas, ev);//get eigenvectors
	if(nv<2)
		return;

#define NORMALIZE(IDX) d=ev[IDX+0]*ev[IDX+0]+ev[IDX+1]*ev[IDX+1]+ev[IDX+2]*ev[IDX+2], d=fabs(d)>1e-9?1/sqrt(d):1, ev[IDX+0]*=d, ev[IDX+1]*=d, ev[IDX+2]*=d
	NORMALIZE(0);//normalize eigenvectors
	NORMALIZE(3);
	NORMALIZE(6);
#undef  NORMALIZE

#define COMPARE(A, B) if(fabs(lambdas[(int)A])<fabs(lambdas[(int)B]))SWAPVAR(A, B, temp)
	COMPARE(d1, d2);//sort eigenvalues
	COMPARE(d2, d3);
	COMPARE(d1, d2);
#undef  COMPARE
	d1*=3;
	d2*=3;
	d3*=3;

	//use only 2 most significant eigenvectors
	memcpy(axis1, ev+d1, sizeof(double[3]));
	memcpy(axis2, ev+d2, sizeof(double[3]));
	CROSS(axis3, axis2, axis1);
	CROSS(axis2, axis1, axis3);
	//#define DOT(A, B) A[0]*B[0]+A[1]*B[1]+A[2]*B[2]
	//double
	//	m1=DOT(axis1, axis2),
	//	m2=DOT(axis2, axis3),
	//	m3=DOT(axis3, axis1);
	//if(m1*m1+m2*m2+m3*m3>1e-3)
	//	LOG_ERROR("");
#define NORMALIZE(X) d=X[0]*X[0]+X[1]*X[1]+X[2]*X[2], d=fabs(d)<1e-9?1/sqrt(d):1, X[0]*=d, X[1]*=d, X[2]*=d
	NORMALIZE(axis1);
	NORMALIZE(axis2);
	NORMALIZE(axis3);
#undef  NORMALIZE
	memcpy(axes+0, axis1, sizeof(double[3]));
	memcpy(axes+3, axis2, sizeof(double[3]));
	memcpy(axes+6, axis3, sizeof(double[3]));
	memcpy(axes+9, mean, sizeof(double[3]));
#if 0
	for(int ky=0;ky<blocksize&&y+ky<ih;++ky)
	{
		for(int kx=0;kx<blocksize&&x+kx<iw;++kx)
		{
			int idx=(iw*(y+ky)+x+kx)<<2;
			double r=src->data[idx|0], g=src->data[idx|1], b=src->data[idx|2];
			r-=mean[0];
			g-=mean[1];
			b-=mean[2];
			double
				r2=axis3[0]*r+axis3[1]*g+axis3[2]*b,//Cr
				g2=axis1[0]*r+axis1[1]*g+axis1[2]*b,//Y
				b2=axis2[0]*r+axis2[1]*g+axis2[2]*b;//Cb
			//double
			//	r2=ev[d3+0]*r+ev[d3+1]*g+ev[d3+2]*b,//Cr
			//	g2=ev[d1+0]*r+ev[d1+1]*g+ev[d1+2]*b,//Y
			//	b2=ev[d2+0]*r+ev[d2+1]*g+ev[d2+2]*b;//Cb
			src->data[idx|0]=(int)r2;
			src->data[idx|1]=(int)g2;
			src->data[idx|2]=(int)b2;
		}
	}
#endif
	//set_window_title("%lf %lf %lf", lambdas[0], lambdas[1], lambdas[2]);
}
void rct_adaptive(Image *src, int fwd)
{
	int blocksize=4, margin=16;
	int bx=(src->iw+blocksize-1)/blocksize;
	int by=(src->ih+blocksize-1)/blocksize;
	int nblocks=bx*by;
	double *axes=(double*)malloc(nblocks*sizeof(double[12]));

	if(!axes)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(int ky=0;ky<by;++ky)
	{
		int ky2=blocksize*ky-margin;
		if(ky2<0)
			ky2=0;
		for(int kx=0;kx<bx;++kx)
		{
			int kx2=blocksize*kx-margin;
			if(kx2<0)
				kx2=0;
			rct_adaptive_block(src->data, src->iw, src->ih, kx2, ky2, blocksize+margin*2, axes+12*(bx*ky+kx));
		}
	}
	for(int ky=0;ky<src->ih;++ky)
	{
		int
			kby1, kby0,
			alphay;

		kby1=(ky+(blocksize>>1))/blocksize;
		kby0=kby1-(kby1>0);
		kby1-=kby1>=by;
		alphay=(ky+(blocksize>>1))%blocksize;
		for(int kx=0;kx<src->iw;++kx)
		{
			int kbx1, kbx0, alphax, idx;
			double *axesTL, *axesTR, *axesBL, *axesBR, localaxes[12], d, r, g, b, r2, g2, b2;

			kbx1=(kx+(blocksize>>1))/blocksize;
			kbx0=kbx1-(kbx1>0);
			kbx1-=kbx1>=bx;
			alphax=(kx+(blocksize>>1))%blocksize;
			axesTL=axes+12*(bx*kby0+kbx0);
			axesTR=axes+12*(bx*kby0+kbx1);
			axesBL=axes+12*(bx*kby1+kbx0);
			axesBR=axes+12*(bx*kby1+kbx1);
			//if(kby0!=kby1&&kbx0!=kbx1)
			//	printf("");
			for(int k=0;k<12;++k)//alpha blend
			{
				double vTL=axesTL[k], vTR=axesTR[k], vBL=axesBL[k], vBR=axesBR[k];
				double vT=vTL+(vTR-vTL)*alphax/blocksize;
				double vB=vBL+(vBR-vBL)*alphax/blocksize;
				double v=vT+(vB-vT)*alphay/blocksize;
				localaxes[k]=v;
			}
#define NORMALIZE(IDX) d=localaxes[IDX+0]*localaxes[IDX+0]+localaxes[IDX+1]*localaxes[IDX+1]+localaxes[IDX+2]*localaxes[IDX+2], d=fabs(d)<1e-9?1/sqrt(d):1, localaxes[IDX+0]*=d, localaxes[IDX+1]*=d, localaxes[IDX+2]*=d
			NORMALIZE(0);//normalize axes
			NORMALIZE(3);
			NORMALIZE(6);
#undef  NORMALIZE
			idx=(src->iw*ky+kx)<<2;
			r=src->data[idx|0], g=src->data[idx|1], b=src->data[idx|2];
			r-=localaxes[ 9];//subtract mean
			g-=localaxes[10];
			b-=localaxes[11];
			r2=localaxes[6]*r+localaxes[7]*g+localaxes[8]*b;//Cr	project on axes
			g2=localaxes[0]*r+localaxes[1]*g+localaxes[2]*b;//Y
			b2=localaxes[3]*r+localaxes[4]*g+localaxes[5]*b;//Cb
			src->data[idx|0]=(int)g2;
			src->data[idx|1]=(int)b2;
			src->data[idx|2]=(int)r2;
		}
	}
	{
		char temp;
		if(fwd)
		{
			ROTATE3(src->depth[0], src->depth[1], src->depth[2], temp);
			src->depth[1]+=src->depth[1]<24;
			src->depth[2]+=src->depth[2]<24;
		}
		else
		{
			src->depth[1]-=src->depth[1]>src->src_depth[1];
			src->depth[2]-=src->depth[2]>src->src_depth[2];
			ROTATE3(src->depth[2], src->depth[1], src->depth[0], temp);
		}
	}
	free(axes);
#if 0
	int res=iw*ih;
	char *dst=(char*)malloc((size_t)res<<2);
	if(!dst)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(dst, src, (size_t)res<<2);
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=(iw*ky+kx)<<2;
			//if(kx==20&&ky==20)//
			//	printf("");
			int nb[ARCT_NNB*3];
			for(int ky2=-ARCT_REACH, idx2=0;ky2<=0;++ky2)//fetch causal neighbors
			{
				for(int kx2=-ARCT_REACH;kx2<=ARCT_REACH;++kx2, idx2+=3)
				{
					if(idx2>=ARCT_NNB*3)
						break;
					if((unsigned)(ky+ky2)<(unsigned)ih&&(unsigned)(kx+kx2)<(unsigned)iw)
					{
						int idx3=(iw*(ky+ky2)+kx+kx2)<<2;
						nb[idx2+0]=pixels[idx3|0];
						nb[idx2+1]=pixels[idx3|1];
						nb[idx2+2]=pixels[idx3|2];
					}
					else
					{
						nb[idx2+0]=0;
						nb[idx2+1]=0;
						nb[idx2+2]=0;
					}
				}
			}
			int mean[3]={0};
			for(int k=0;k<ARCT_NNB*3-2;k+=3)//calculate mean
			{
				mean[0]+=nb[k+0];
				mean[1]+=nb[k+1];
				mean[2]+=nb[k+2];
			}
			mean[0]=(mean[0]+(ARCT_NNB>>1))/ARCT_NNB;
			mean[1]=(mean[1]+(ARCT_NNB>>1))/ARCT_NNB;
			mean[2]=(mean[2]+(ARCT_NNB>>1))/ARCT_NNB;
			default_RCT(src+idx, dst+idx, fwd);
			continue;

			for(int k=0;k<ARCT_NNB*3-2;k+=3)//subtract mean
			{
				nb[k+0]-=mean[0];
				nb[k+1]-=mean[1];
				nb[k+2]-=mean[2];
			}
#if 1
			double cov[9]={0};
			for(int k=0;k<ARCT_NNB*3-2;k+=3)//calculate covariance matrix
			{
				double r=nb[k+0], g=nb[k+1], b=nb[k+2];
				double rr=r*r, gg=g*g, bb=b*b, rg=r*g, gb=g*b, br=b*r;
				cov[0]+=rr;
				cov[1]+=rg;
				cov[2]+=br;
				//cov[3]+=rg;//cov[1]
				cov[4]+=gg;
				cov[5]+=gb;
				//cov[6]+=br;//cov[2]
				//cov[7]+=gb;//cov[5]
				cov[8]+=bb;
			}
			cov[3]=cov[1];
			cov[6]=cov[2];
			cov[7]=cov[5];
			for(int k=0;k<9;++k)
				cov[k]/=ARCT_NNB-1;
			
			//double det=cov[0]*(cov[4]*cov[8]-cov[5]*cov[7]) - cov[1]*(cov[3]*cov[8]-cov[5]*cov[6]) + cov[2]*(cov[3]*cov[7]-cov[4]*cov[6]);
			//if(det<1e-6)
			//{
			//	default_RCT(src+idx, dst+idx, fwd);
			//	continue;
			//}

			//calculate eigenvectors & eigenvalues
			double lambdas[6];
			impl_egval3(cov, lambdas);//get eigenvalues
			if(fabs(lambdas[1])>1e-6||fabs(lambdas[3])>1e-6||fabs(lambdas[5])>1e-6)//reject complex eigenvalues
			{
				default_RCT(src+idx, dst+idx, fwd);
				continue;
			}
			lambdas[0]=sqrt(lambdas[0]*lambdas[0]+lambdas[1]*lambdas[1]);
			lambdas[1]=sqrt(lambdas[2]*lambdas[2]+lambdas[3]*lambdas[3]);
			lambdas[2]=sqrt(lambdas[4]*lambdas[4]+lambdas[5]*lambdas[5]);
			if(fabs(lambdas[0])<1e-6||fabs(lambdas[1])<1e-6||fabs(lambdas[2])<1e-6)//reject zero eigenvalues
			{
				default_RCT(src+idx, dst+idx, fwd);
				continue;
			}
			double ev[9];
			int nv=impl_egvec(cov, 3, lambdas, ev);//get eigenvectors
			if(nv<2)
			{
				default_RCT(src+idx, dst+idx, fwd);
				continue;
			}
			double d;
#define NORMALIZE(IDX) d=ev[IDX+0]*ev[IDX+0]+ev[IDX+1]*ev[IDX+1]+ev[IDX+2]*ev[IDX+2], d=fabs(d)<1e-9?1/sqrt(d):1, ev[IDX+0]*=d, ev[IDX+1]*=d, ev[IDX+2]*=d
			NORMALIZE(0);//normalize eigenvectors
			NORMALIZE(3);
			NORMALIZE(6);
#undef  NORMALIZE
			char d1=0, d2=1, d3=2, temp;
#define COMPARE(A, B) if(fabs(lambdas[A])<fabs(lambdas[B]))SWAPVAR(A, B, temp)
			COMPARE(d1, d2);//sort eigenvalues
			COMPARE(d2, d3);
			COMPARE(d1, d2);
#undef  COMPARE
			d1*=3;
			d2*=3;
			d3*=3;

			//use only 2 most significant eigenvectors
			double axis1[3], axis2[3], axis3[3];
			memcpy(axis1, ev+d1, sizeof(double[3]));
			memcpy(axis2, ev+d2, sizeof(double[3]));
			CROSS(axis3, axis2, axis1);
			CROSS(axis2, axis1, axis3);
			//#define DOT(A, B) A[0]*B[0]+A[1]*B[1]+A[2]*B[2]
			//double
			//	m1=DOT(axis1, axis2),
			//	m2=DOT(axis2, axis3),
			//	m3=DOT(axis3, axis1);
			//if(m1*m1+m2*m2+m3*m3>1e-3)
			//	LOG_ERROR("");
#define NORMALIZE(X) d=X[0]*X[0]+X[1]*X[1]+X[2]*X[2], d=fabs(d)<1e-9?1/sqrt(d):1, X[0]*=d, X[1]*=d, X[2]*=d
			NORMALIZE(axis1);
			NORMALIZE(axis2);
			NORMALIZE(axis3);
#undef  NORMALIZE
			char r=pixels[idx|0], g=pixels[idx|1], b=pixels[idx|2];
			r-=mean[0];
			g-=mean[1];
			b-=mean[2];
			double
				r2=ev[d3+0]*r+ev[d3+1]*g+ev[d3+2]*b,//Cr
				g2=ev[d1+0]*r+ev[d1+1]*g+ev[d1+2]*b,//Y
				b2=ev[d2+0]*r+ev[d2+1]*g+ev[d2+2]*b;//Cb
			dst[idx|0]=(char)r2;
			dst[idx|1]=(char)g2;
			dst[idx|2]=(char)b2;
#endif
#if 0
			//check for degenerate tetrahedron
			{
				double ab_ac[]=
				{
					(nb[1*3+0]-nb[0*3+0])*(nb[2*3+0]-nb[0*3+0]),
					(nb[1*3+1]-nb[0*3+1])*(nb[2*3+1]-nb[0*3+1]),
					(nb[1*3+2]-nb[0*3+2])*(nb[2*3+2]-nb[0*3+2]),
				};
				double ad[]=
				{
					nb[3*3+0]-nb[0*3+0],
					nb[3*3+1]-nb[0*3+1],
					nb[3*3+2]-nb[0*3+2],
				};
				double temp[3];
				CROSS(temp, ab_ac, ad);
				double volume=temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2];
				if(volume<1e-10)//shortcut for degenerate neighbors	JPEG2000
				{
					if(fwd)
					{
						r-=g;
						b-=g;
						g+=(r+b)>>2;
					}
					else
					{
						g-=(r+b)>>2;
						b+=g;
						r+=g;
					}
					dst[idx|0]=r;
					dst[idx|1]=g;
					dst[idx|2]=b;
					continue;
				}
			}
			int dist[ARCT_NNB]={0};
			for(int k=0, k2=0;k<ARCT_NNB*3-2;k+=3, ++k2)//calculate vertex distances from mean
			{
				int r=nb[k+0], g=nb[k+1], b=nb[k+2];
				dist[k2]=r*r+b*b+g*g;
			}
			int d1=0, d2=0;
			for(int k=1;k<ARCT_NNB;++k)//select farthest vertex
			{
				if(dist[d1]<dist[k])
					d1=k;
			}
			for(int k=1;k<ARCT_NNB;++k)//select 2nd farthest vertex
			{
				if(k!=d1&&dist[d2]<dist[k])
					d2=k;
			}

			//define new axes
			double
				inv_d1=dist[d1]?1./sqrt(dist[d1]):1,
				inv_d2=dist[d2]?1./sqrt(dist[d2]):1;
			d1*=3;
			d2*=3;
			double
				axis1[]={nb[d1+0]*inv_d1, nb[d1+1]*inv_d1, nb[d1+2]*inv_d1},
				axis2[]={nb[d2+0]*inv_d2, nb[d2+1]*inv_d2, nb[d2+2]*inv_d2},
				axis3[3];
			//int axis1[3], axis2[3], axis3[3];
			//memcpy(axis1, nb+d1, sizeof(int[3]));
			//memcpy(axis2, nb+d2, sizeof(int[3]));
			CROSS(axis3, axis2, axis1);//axis3 = d2 cross d1
			CROSS(axis2, axis1, axis3);//axis2 = d1 cross axis3

			inv_d1=axis3[0]*axis3[0]+axis3[1]*axis3[1]+axis3[2]*axis3[2];
			inv_d1=inv_d1?1/sqrt(inv_d1):1;
			axis3[0]*=inv_d1;
			axis3[1]*=inv_d1;
			axis3[2]*=inv_d1;

			inv_d1=axis2[0]*axis2[0]+axis2[1]*axis2[1]+axis2[2]*axis2[2];
			inv_d1=inv_d1?1/sqrt(inv_d1):1;
			axis2[0]*=inv_d1;
			axis2[1]*=inv_d1;
			axis2[2]*=inv_d1;
			
			double
				r2=axis3[0]*r+axis3[1]*g+axis3[2]*b,
				g2=axis1[0]*r+axis1[1]*g+axis1[2]*b,
				b2=axis2[0]*r+axis2[1]*g+axis2[2]*b;
			dst[idx|0]=(char)r2;
			dst[idx|1]=(char)g2;
			dst[idx|2]=(char)b2;
#endif
		}
	}
	memcpy(src, dst, (size_t)res<<2);
	free(dst);
#endif
	//const int blocksize=4;
	//for(int ky=0;ky<ih-(blocksize-1);ky+=blocksize)
	//{
	//	for(int kx=0;kx<iw-(blocksize-1);kx+=blocksize)
	//	{
	//	}
	//}
#if 0
	int A=0x8000, C=0x8000;
	int LR=(int)(customparam_st[0]*1000);
	for(ptrdiff_t k=0, len=(ptrdiff_t)iw*ih*4;k<len;k+=4)
	{
		char r=src[k], g=src[k|1], b=src[k|2];
		char r0, g0, b0, Cr, Y, Cb;

		if(fwd)
		{
			r0=r, g0=g, b0=b;
			
			r-=g;
			b-=g;
			g+=(A*r+C*b+0x8000)>>16;//minimize abs(Y) = abs(g+A*(r-g)+C*(b-g))

			Cr=r, Y=g, Cb=b;
		}
		else
		{
			Cr=r, Y=g, Cb=b;

			g-=(A*r+C*b+0x8000)>>16;
			b+=g;
			r+=g;

			r0=r, g0=g, b0=b;
		}
		
		src[k  ]=r;//Cr
		src[k|1]=g;//Y
		src[k|2]=b;//Cb

		//backprop	d/dA abs(Y) = sgn(Y)*(r-g)	d/dB abs(Y) = sgn(Y)*(b-g)
		int sgnY=(Y>0)-(Y<0);
		A-=sgnY*Cr*LR/1000;
		C-=sgnY*Cb*LR/1000;
	}
#endif
	//int res=iw*ih;
	//char *buf=(char*)malloc((size_t)res<<2);
	//if(!buf)
	//{
	//	LOG_ERROR("Alloc error");
	//	return;
	//}
	//memcpy(buf, src, (size_t)res<<2);
#if 0
	for(int ky=0, idx=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx, ++idx)
		{
			char r=src[idx<<2|0], g=src[idx<<2|1], b=src[idx<<2|2];
			char lumas[]=
			{
				(r+g+b)/3,
				(r+r+g+b)>>2,
				(r+g+g+b)>>2,
				(r+g+b+b)>>2,
				(5*r+5*g+6*b)>>4,
				(9*r*9*g+14*b)>>4,
			};
			int kmin=0;
			for(int k=0;k<(int)_countof(lumas);++k)
			{
				if(abs(lumas[kmin])>abs(lumas[k]))
					kmin=k;
			}
			int color=0;
			switch(kmin)
			{
			case 0:color=0x0000FF;break;
			case 1:color=0x00C0C0;break;
			case 2:color=0x00FF00;break;
			case 3:color=0xC0C000;break;
			case 4:color=0xFF0000;break;
			case 5:color=0xC000C0;break;
			}
			char *p=(char*)&color;
			src[idx<<2|0]=p[0]-128;
			src[idx<<2|1]=p[1]-128;
			src[idx<<2|2]=p[2]-128;
		}
	}
#endif
}



//chroma from luma (CfL)
#if 0
short cfl_params[4]={0};//{alpha, beta} x2 chroma channels
void pred_cfl_calcparams(char *buf, int iw, int ih, int x1, int x2, int y1, int y2, int kl, int kc, short *params)
{
	double alpha, beta, term[4]={0};
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			int idx=(iw*ky+kx)<<2;
			char L=buf[idx|kl], C=buf[idx|kc];
			term[0]+=L*C;
			term[1]+=L;
			term[2]+=C;
			term[3]+=L*L;
		}
	}
	int bw=x2-x1, bh=y2-y1, count=bw*bh;
	if(!count||!term[1]&&!term[3])
	{
		params[0]=0, params[1]=0;
		return;
	}
	alpha=(count*term[0]-term[1]*term[2])/(count*term[3]-term[1]*term[1]);
	beta=(term[2]-alpha*term[1])/count;
	params[0]=(short)(alpha*256);
	params[1]=(short)(beta*256);
}
void pred_cfl(char *buf, int iw, int ih, int fwd)
{
#if 0
	if(fwd)
	{
		pred_cfl_calcparams(buf, iw, ih, 0, iw, 0, ih, 1, 0, cfl_params);
		pred_cfl_calcparams(buf, iw, ih, 0, iw, 0, ih, 1, 2, cfl_params+2);
	}
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
			int idx=(iw*ky+kx)<<2;
			char Co=buf[idx], Y=buf[idx|1], Cb=buf[idx|2];
			int pred_o=(cfl_params[0]*Y+cfl_params[1]+127)>>8;
			int pred_b=(cfl_params[2]*Y+cfl_params[3]+127)>>8;
			if(fwd)
			{
				Co-=pred_o;
				Cb-=pred_b;
			}
			else
			{
				Co+=pred_o;
				Cb+=pred_b;
			}
			buf[idx  ]=Co;
			buf[idx|2]=Cb;
		}
	}
#endif
}
#endif


//spatial transforms
void filt_deint422(Image *src)
{
	for(int ky=0;ky<src->ih;++ky)
	{
		for(int kx=1;kx<src->iw-1;kx+=2)
		{
			int
				duW=src->data[(src->iw*(ky+0)+kx+0)<<2|1],
				duE=src->data[(src->iw*(ky+0)+kx+1)<<2|1],
				dvW=src->data[(src->iw*(ky+0)+kx+0)<<2|2],
				dvE=src->data[(src->iw*(ky+0)+kx+1)<<2|2];
			int
				suW=(+3*duW-1*duE+1)>>1,
				suE=(-1*duW+3*duE+1)>>1;
			int
				svW=(+3*dvW-1*dvE+1)>>1,
				svE=(-1*dvW+3*dvE+1)>>1;
			src->data[(src->iw*(ky+0)+kx+0)<<2|1]=suW;
			src->data[(src->iw*(ky+0)+kx+1)<<2|1]=suE;
			src->data[(src->iw*(ky+0)+kx+0)<<2|2]=svW;
			src->data[(src->iw*(ky+0)+kx+1)<<2|2]=svE;
		}
	}
}
void filt_deint420(Image *src)
{
	for(int ky=1;ky<src->ih-1;ky+=2)
	{
		for(int kx=1;kx<src->iw-1;kx+=2)
		{
			int
				duNW=src->data[(src->iw*(ky+0)+kx+0)<<2|1],
				duNE=src->data[(src->iw*(ky+0)+kx+1)<<2|1],
				duSW=src->data[(src->iw*(ky+1)+kx+0)<<2|1],
				duSE=src->data[(src->iw*(ky+1)+kx+1)<<2|1],
				dvNW=src->data[(src->iw*(ky+0)+kx+0)<<2|2],
				dvNE=src->data[(src->iw*(ky+0)+kx+1)<<2|2],
				dvSW=src->data[(src->iw*(ky+1)+kx+0)<<2|2],
				dvSE=src->data[(src->iw*(ky+1)+kx+1)<<2|2];
			int
				suNW=(+9*duNW-3*duNE-3*duSW+1*duSE+2)>>2,
				suNE=(-3*duNW+9*duNE+1*duSW-3*duSE+2)>>2,
				suSW=(-3*duNW+1*duNE+9*duSW-3*duSE+2)>>2,
				suSE=(+1*duNW-3*duNE-3*duSW+9*duSE+2)>>2;
			int
				svNW=(+9*dvNW-3*dvNE-3*dvSW+1*dvSE+2)>>2,
				svNE=(-3*dvNW+9*dvNE+1*dvSW-3*dvSE+2)>>2,
				svSW=(-3*dvNW+1*dvNE+9*dvSW-3*dvSE+2)>>2,
				svSE=(+1*dvNW-3*dvNE-3*dvSW+9*dvSE+2)>>2;
			src->data[(src->iw*(ky+0)+kx+0)<<2|1]=suNW;
			src->data[(src->iw*(ky+0)+kx+1)<<2|1]=suNE;
			src->data[(src->iw*(ky+1)+kx+0)<<2|1]=suSW;
			src->data[(src->iw*(ky+1)+kx+1)<<2|1]=suSE;
			src->data[(src->iw*(ky+0)+kx+0)<<2|2]=svNW;
			src->data[(src->iw*(ky+0)+kx+1)<<2|2]=svNE;
			src->data[(src->iw*(ky+1)+kx+0)<<2|2]=svSW;
			src->data[(src->iw*(ky+1)+kx+1)<<2|2]=svSE;
		}
	}
}

void packsign(Image *src, int fwd)
{
	ptrdiff_t nvals=(ptrdiff_t)src->iw*src->ih<<2;
	int half[]=
	{
		1<<src->depth[0]>>1,
		1<<src->depth[1]>>1,
		1<<src->depth[2]>>1,
		1<<src->depth[3]>>1,
	};
	if(fwd)
	{
		for(ptrdiff_t k=0;k<nvals;k+=4)
		{
			for(int kc=0;kc<src->nch;++kc)
			{
				int val=src->data[k|kc];
				val=val<<1^-(val<0);
				val-=half[kc];
				src->data[k|kc]=val;
			}
		}
	}
	else
	{
		for(ptrdiff_t k=0;k<nvals;k+=4)
		{
			for(int kc=0;kc<src->nch;++kc)
			{
				int val=src->data[k|kc];
				val+=half[kc];
				val=val>>1^-(val&1);
				src->data[k|kc]=val;
			}
		}
	}
}
void pred_MTF(Image *src, int fwd)
{
	int maxdepth, *dict;

	maxdepth=MAXVAR(src->depth[0], src->depth[1]);
	UPDATE_MAX(maxdepth, src->depth[2]);
	UPDATE_MAX(maxdepth, src->depth[3]);
	dict=(int*)malloc(sizeof(int[4])<<maxdepth);
	if(!dict)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(int kc=0;kc<4;++kc)
	{
		int depth, *curr_dict, nlevels;

		depth=src->depth[kc];
		if(!depth)
			continue;
		curr_dict=dict+((size_t)kc<<maxdepth);
		nlevels=1<<src->depth[kc];
		for(int ks=0;ks<nlevels;++ks)
			curr_dict[ks]=ks;
	}
	if(fwd)
	{
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, idx+=4)
			{
				for(int kc=0;kc<4;++kc)
				{
					int depth, nlevels, half, *curr_dict, curr, ks;

					depth=src->depth[kc];
					if(!depth)
						continue;
					nlevels=1<<depth;
					half=nlevels>>1;
					curr_dict=dict+((size_t)kc<<maxdepth);
					curr=src->data[idx|kc]+half;
					ks=0;
					for(;ks<nlevels&&curr_dict[ks]!=curr;++ks);
					//if(ks>=nlevels)
					//	LOG_ERROR("");
					src->data[idx|kc]=ks-half;
					if(ks)
					{
						memmove(curr_dict+1, curr_dict, ks*sizeof(int));
						curr_dict[0]=curr;
					}
				}
			}
		}
	}
	else
	{
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, idx+=4)
			{
				for(int kc=0;kc<4;++kc)
				{
					int depth, nlevels, half, *curr_dict, curr, ks;

					depth=src->depth[kc];
					if(!depth)
						continue;
					nlevels=1<<depth;
					half=nlevels>>1;
					curr_dict=dict+((size_t)kc<<maxdepth);
					ks=src->data[idx|kc]+half;
					curr=curr_dict[ks];
					src->data[idx|kc]=curr-half;
					if(ks)
					{
						memmove(curr_dict+1, curr_dict, ks*sizeof(int));
						curr_dict[0]=curr;
					}
				}
			}
		}
	}
	free(dict);
}

//clamped gradient / LOCO-I / Median Edge Detector (MED) predictor from JPEG-LS
void pred_clampgrad(Image *src, int fwd, int enable_ma)
{
#if 1
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	__m128i mhalf=_mm_set_epi32(
		nlevels[3]>>1,
		nlevels[2]>>1,
		nlevels[1]>>1,
		nlevels[0]>>1
	);
	__m128i symmask=_mm_set_epi32(
		nlevels[3]-1,
		nlevels[2]-1,
		nlevels[1]-1,
		nlevels[0]-1
	);
	ALIGN(16) int curr[4]={0};
	int *pixels=(int*)_mm_malloc((src->iw+4LL)*sizeof(int[2*4]), sizeof(__m128i));//2 padded rows * 4 channels max

	(void)enable_ma;

	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, (src->iw+4LL)*sizeof(int[2*4]));
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+(((src->iw+4LL)*((ky-0LL)&1)+1)<<2),
			pixels+(((src->iw+4LL)*((ky-1LL)&1)+1)<<2),
		};
		if(fwd)
		{
			for(int kx=0;kx<src->iw;++kx, idx+=4)
			{
				__m128i NW	=_mm_load_si128((__m128i*)rows[1]-1);
				__m128i N	=_mm_load_si128((__m128i*)rows[1]+0);
				__m128i W	=_mm_load_si128((__m128i*)rows[0]-1);
				__m128i pred=_mm_add_epi32(N, W);
				pred=_mm_sub_epi32(pred, NW);
				{
					__m128i vmin=_mm_min_epi32(N, W);
					__m128i vmax=_mm_max_epi32(N, W);
					pred=_mm_min_epi32(pred, vmax);
					pred=_mm_max_epi32(pred, vmin);
				}
				{
					__m128i mc=_mm_loadu_si128((__m128i*)(src->data+idx));
					pred=_mm_sub_epi32(mc, pred);
					pred=_mm_add_epi32(pred, mhalf);
					pred=_mm_and_si128(pred, symmask);
					pred=_mm_sub_epi32(pred, mhalf);
					_mm_store_si128((__m128i*)curr, pred);
					src->data[idx|0]=curr[0];
					src->data[idx|1]=curr[1];
					src->data[idx|2]=curr[2];
					_mm_store_si128((__m128i*)rows[0], mc);
				}
				rows[0]+=4;
				rows[1]+=4;
			}
		}
		else
		{
			for(int kx=0;kx<src->iw;++kx, idx+=4)
			{
				__m128i NW	=_mm_load_si128((__m128i*)rows[1]-1);
				__m128i N	=_mm_load_si128((__m128i*)rows[1]+0);
				__m128i W	=_mm_load_si128((__m128i*)rows[0]-1);
				__m128i pred=_mm_add_epi32(N, W);
				pred=_mm_sub_epi32(pred, NW);
				{
					__m128i vmin=_mm_min_epi32(N, W);
					__m128i vmax=_mm_max_epi32(N, W);
					pred=_mm_min_epi32(pred, vmax);
					pred=_mm_max_epi32(pred, vmin);
				}
				{
					__m128i mc=_mm_loadu_si128((__m128i*)(src->data+idx));
					pred=_mm_add_epi32(mc, pred);
					pred=_mm_add_epi32(pred, mhalf);
					pred=_mm_and_si128(pred, symmask);
					pred=_mm_sub_epi32(pred, mhalf);
					_mm_store_si128((__m128i*)curr, pred);
					src->data[idx|0]=curr[0];
					src->data[idx|1]=curr[1];
					src->data[idx|2]=curr[2];
					_mm_store_si128((__m128i*)rows[0], pred);
				}
				rows[0]+=4;
				rows[1]+=4;
			}
		}
#if 0
		for(;kx<src->iw;++kx, idx+=4)
		{
			//if(ky==src->ih-1&&kx==src->iw-1)//
			//	printf("");

			__m128i NW	=_mm_load_si128((__m128i*)rows[1]-1);
			__m128i N	=_mm_load_si128((__m128i*)rows[1]+0);
			__m128i W	=_mm_load_si128((__m128i*)rows[0]-1);
			__m128i pred=_mm_add_epi32(N, W);
			pred=_mm_sub_epi32(pred, NW);
			__m128i vmin=_mm_min_epi32(N, W);
			__m128i vmax=_mm_max_epi32(N, W);
			pred=_mm_min_epi32(pred, vmax);
			pred=_mm_max_epi32(pred, vmin);

			__m128i mc=_mm_loadu_si128((__m128i*)(src->data+idx));
			pred=_mm_xor_si128(pred, mfwd);
			pred=_mm_sub_epi32(pred, mfwd);
			pred=_mm_add_epi32(pred, mc);
			pred=_mm_add_epi32(pred, mhalf);
			pred=_mm_and_si128(pred, symmask);
			pred=_mm_sub_epi32(pred, mhalf);
			_mm_store_si128((__m128i*)curr, pred);
			src->data[idx|0]=curr[0];
			src->data[idx|1]=curr[1];
			src->data[idx|2]=curr[2];
			_mm_store_si128((__m128i*)rows[0], fwd?mc:pred);

			rows[0]+=4;
			rows[1]+=4;
		}
#endif
	}
	_mm_free(pixels);
#endif
#if 0
	short *pixels=(short*)malloc((src->iw+2LL)*sizeof(short[2*4]));//2 padded rows * 4 channels max
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, (src->iw+2LL)*sizeof(short[2*4]));
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int fwdmask=-fwd;
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		short *rows[]=
		{
			pixels+(((src->iw+2LL)*((ky-0LL)&1)+1)<<2),
			pixels+(((src->iw+2LL)*((ky-1LL)&1)+1)<<2),
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(int kc=0;kc<src->nch;++kc)
			{
				//if(ky==100&&kx==100)//
				//	printf("");

				int
					NW	=rows[1][kc-4],
					N	=rows[1][kc+0],
					W	=rows[0][kc-4];
				int pred=N+W-NW;
				int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
				pred=CLAMP(vmin, pred, vmax);

				int curr=src->data[idx+kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=curr;

				pred+=nlevels[kc]>>1;
				pred&=nlevels[kc]-1;
				pred-=nlevels[kc]>>1;

				src->data[idx+kc]=pred;
				rows[0][kc]=fwd?curr:pred;
			}

			rows[0]+=4;
			rows[1]+=4;
		}
	}
	free(pixels);
#endif
#if 0
	Image *dst=0;
	image_copy(&dst, src);
	if(!dst)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	const int *pixels=fwd?src->data:dst->data;
//	const int *errors=fwd?dst->data:src->data;
	for(int kc=0;kc<4;++kc)
	{
		if(!src->depth[kc])
			continue;
		int nlevels=1<<src->depth[kc];
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				//if(ky==100&&kx==100)//
				//	printf("");

				int
					NW=kx&&ky?pixels[(idx-src->iw-1)<<2|kc]:0,
					N =ky    ?pixels[(idx-src->iw  )<<2|kc]:0,
					W =kx    ?pixels[(idx        -1)<<2|kc]:0;
				int pred=N+W-NW;
				pred=MEDIAN3(N, W, pred);

				//int p=N+W-NW;
				//int pred;
				//if(p>=-(nlevels>>1)&&p<(nlevels>>1))
				//	pred=p;
				//else
				//{
				//	int pN=abs(p-N), pW=abs(p-W), pNW=abs(p-NW);
				//	if(pN<=pW&&pN<=pNW)
				//		pred=N;
				//	else if(pW<=pNW)
				//		pred=W;
				//	else
				//		pred=NW;
				//}
				
				pred^=-fwd;
				pred+=fwd;
				pred+=src->data[idx<<2|kc];
				if(enable_ma)
				{
					pred+=nlevels>>1;
					pred&=nlevels-1;
					pred-=nlevels>>1;
				}
				dst->data[idx<<2|kc]=pred;
			}
		}
	}
	memcpy(src->data, dst->data, (size_t)src->iw*src->ih*sizeof(int[4]));
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}
	free(dst);
#endif
}
void pred_CG420(Image *src, int fwd)
{
	int nch;
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int fwdmask=-fwd;

	int bufsize=(src->iw+16LL)*sizeof(int[4*4]);//4 padded rows * 4 channels max
	int *pixels=(int*)malloc(bufsize);

	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8)*4,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8)*4,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8)*4,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8)*4,
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(int kc=0;kc<3;++kc)
			{
				int pred, curr;
				int
					NW	=rows[1][kc-1*4],
					N	=rows[1][kc+0*4],
					W	=rows[0][kc-1*4];

				if(kc)
				{
					if(ky&1)
					{
						if(kx&1)
						{
							MEDIAN3_32(pred, N, W, N+W-NW);
						}
						else
							pred=N;
					}
					else
					{
						if(kx&1)
							pred=W;
						else
						{
							MEDIAN3_32(pred, N, W, N+W-NW);
						}
					}
				}
				else
					MEDIAN3_32(pred, N, W, N+W-NW);

				curr=src->data[idx+kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=curr;

				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];

				src->data[idx+kc]=pred;
				rows[0][kc]=fwd?curr:pred;
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	free(pixels);
}
void pred_CG422(Image *src, int fwd)
{
	int nch;
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int fwdmask=-fwd;

	int bufsize=(src->iw+16LL)*sizeof(int[4*4]);//4 padded rows * 4 channels max
	int *pixels=(int*)malloc(bufsize);

	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8)*4,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8)*4,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8)*4,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8)*4,
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(int kc=0;kc<3;++kc)
			{
				int pred, curr;
				int
					NW	=rows[1][kc-1*4],
					N	=rows[1][kc+0*4],
					W	=rows[0][kc-1*4];

				if(kc)
				{
					if(kx&1)
						pred=W;
					else
					{
						MEDIAN3_32(pred, N, W, N+W-NW);
					}
				}
				else
					MEDIAN3_32(pred, N, W, N+W-NW);

				curr=src->data[idx+kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=curr;

				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];

				src->data[idx+kc]=pred;
				rows[0][kc]=fwd?curr:pred;
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	free(pixels);
}


//	#define CG3D_ENABLE_MA

void pred_CG3D(Image *src, int fwd, int enable_ma)
{
	int nch;
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int chpermutations[]=
	{
		0, 1, 2, 3,//subR
		1, 2, 0, 3,//subG
	};
	int
		nlumas=1,
		chpoints[3]={0},
		*perm,
		fwdmask=-fwd;

	int *pixels=(int*)malloc((src->iw+2LL)*sizeof(int[2*4]));//2 padded rows * 4 channels max

	(void)enable_ma;

	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, (src->iw+2LL)*sizeof(int[2*4]));
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
#ifdef CG3D_ENABLE_MA
	for(int kc=0;kc<4;++kc)
	{
		if(src->depth[kc])
			src->depth[kc]+=fwd;
	}
#else
	//src->depth[1]-=src->depth[1]>src->src_depth[1];
#endif
	for(int kc=1;kc<nch;++kc)
		nlumas+=src->depth[kc]==src->depth[0];
	chpoints[1]=nlumas;
	chpoints[2]=nch;
	perm=nlumas==nch?chpermutations+4:chpermutations;//permute channels RGB->GBR if no RCT was applied earlier
	//int pred[3]={0}, vmin[3]={0}, vmax[3]={0};
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+(((src->iw+2LL)*((ky-0LL)&1)+1)<<2),
			pixels+(((src->iw+2LL)*((ky-1LL)&1)+1)<<2),
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
#if 1
			for(int kk=0;kk<2;++kk)
			{
				for(int kc0=chpoints[kk];kc0<chpoints[kk+1];++kc0)
				{
					int kc=perm[kc0], pred, curr;
					int
						NW	=rows[1][kc-4],
						N	=rows[1][kc+0],
						W	=rows[0][kc-4],
						offset	=0;
					if(kc0>chpoints[kk])
						offset+=rows[0][perm[chpoints[kk]]];

					MEDIAN3_32(pred, N, W, N+W-NW);
					//pred=(N+W)>>1;
					//pred=N+W-NW;//X

					pred+=offset;
					CLAMP2_32(pred, pred, -halfs[kc], halfs[kc]-1);

					curr=src->data[idx+kc];
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=curr;

					//if(kc!=1)
					//{
						pred<<=32-src->depth[kc];
						pred>>=32-src->depth[kc];
					//}
					//pred+=nlevels[kc]>>1;
					//pred&=nlevels[kc]-1;
					//pred-=nlevels[kc]>>1;

					src->data[idx+kc]=pred;
					rows[0][kc]=(fwd?curr:pred)-offset;
				}
			}
#endif
#if 0
			if(ky==100&&kx==100)
				printf("");
			int kc=0;
			if(nch>=3)
			{
				int
					vNW	=rows[1][0-4],
					vN	=rows[1][0+0],
					vW	=rows[0][0-4];
				int
					yNW	=rows[1][1-4],
					yN	=rows[1][1+0],
					yW	=rows[0][1-4];
				int
					uNW	=rows[1][2-4],
					uN	=rows[1][2+0],
					uW	=rows[0][2-4];
				pred[0]=vN+vW-vNW;
				vmin[0]=MINVAR(vN, vW), vmax[0]=MAXVAR(vN, vW);
				pred[0]=CLAMP(vmin[0], pred[0], vmax[0]);

				pred[1]=yN+yW-yNW;
				vmin[1]=MINVAR(yN, yW), vmax[1]=MAXVAR(yN, yW);
				pred[1]=CLAMP(vmin[1], pred[1], vmax[1]);

				pred[2]=uN+uW-uNW;
				vmin[2]=MINVAR(uN, uW), vmax[2]=MAXVAR(uN, uW);
				pred[2]=CLAMP(vmin[2], pred[2], vmax[2]);
				
				//pred[1]-=(pred[0]+pred[2])>>2;
				pred[2]+=pred[1];
				pred[0]+=pred[1];//X  adding pred, not prev
				pred[0]=CLAMP(-halfs[kc], pred[0], halfs[kc]-1);
				pred[1]=CLAMP(-halfs[kc], pred[1], halfs[kc]-1);
				pred[2]=CLAMP(-halfs[kc], pred[2], halfs[kc]-1);
				int curr[]=
				{
					src->data[idx+0],
					src->data[idx+1],
					src->data[idx+2],
				};
				pred[0]^=fwdmask;
				pred[1]^=fwdmask;
				pred[2]^=fwdmask;
				pred[0]-=fwdmask;
				pred[1]-=fwdmask;
				pred[2]-=fwdmask;
				pred[0]+=curr[0];
				pred[1]+=curr[1];
				pred[2]+=curr[2];
				pred[0]+=nlevels[0]>>1;
				pred[1]+=nlevels[1]>>1;
				pred[2]+=nlevels[2]>>1;
				pred[0]&=nlevels[0]-1;
				pred[1]&=nlevels[1]-1;
				pred[2]&=nlevels[2]-1;
				pred[0]-=nlevels[0]>>1;
				pred[1]-=nlevels[1]>>1;
				pred[2]-=nlevels[2]>>1;
				src->data[idx+0]=pred[0];
				src->data[idx+1]=pred[1];
				src->data[idx+2]=pred[2];
				if(fwd)
				{
					curr[0]-=curr[1];
					curr[2]-=curr[1];
					//curr[1]+=(curr[0]+curr[2])>>2;
					rows[0][0]=curr[0];
					rows[0][1]=curr[1];
					rows[0][2]=curr[2];
				}
				else
				{
					pred[0]-=pred[1];
					pred[2]-=pred[1];
					//pred[1]+=(pred[0]+pred[2])>>2;
					rows[0][0]=pred[0];
					rows[0][1]=pred[1];
					rows[0][2]=pred[2];
				}
				kc=3;
			}
			while(kc<nch)
			{
				int
					NW	=rows[1][kc-4],
					N	=rows[1][kc+0],
					W	=rows[0][kc-4];
				int pred=N+W-NW;
				int vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
				pred=CLAMP(vmin, pred, vmax);

				int curr=src->data[idx+kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=curr;

				pred+=nlevels[kc]>>1;
				pred&=nlevels[kc]-1;
				pred-=nlevels[kc]>>1;

				src->data[idx+kc]=pred;
				rows[0][kc]=fwd?curr:pred;
				++kc;
			}
#endif
			rows[0]+=4;
			rows[1]+=4;
		}
	}
#ifdef CG3D_ENABLE_MA
	for(int kc=0;kc<4;++kc)
	{
		if(src->depth[kc])
			src->depth[kc]-=!fwd;
	}
#else
	//src->depth[1]+=src->depth[1]<24;
#endif
	free(pixels);
}


	#define PU_UPDATE_LUMA

void pred_PU(Image *src, int fwd)
{
	int nch;
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int perm[]={1, 2, 0, 3};
	int fwdmask=-fwd;

	size_t bufsize=(src->iw+4LL)*sizeof(int[4*4*2]);//4 padded rows * 4 channels max * {pixels, errors}
	int *pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+(((src->iw+4LL)*((ky-0LL)&3)+2)<<3),
			pixels+(((src->iw+4LL)*((ky-1LL)&3)+2)<<3),
			pixels+(((src->iw+4LL)*((ky-2LL)&3)+2)<<3),
			pixels+(((src->iw+4LL)*((ky-3LL)&3)+2)<<3),
		};
		//__m128i vmin=_mm_set_epi32(-halfs[3], -halfs[2], -halfs[1], -halfs[0]);
		//__m128i vmax=_mm_set_epi32(halfs[3]-1, halfs[2]-1, halfs[1]-1, halfs[0]-1);
		//ALIGN(16) int preds[4]={0};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
#if 1
#ifdef PU_UPDATE_LUMA
			if(!fwd)
			{
				{
					int *ptr=src->data+idx, temp;
					ROTATE3(ptr[2], ptr[1], ptr[0], temp);
				}
				int
					//NN	=rows[2][1+0*8],
					//NW	=rows[1][1-1*8],
					//N	=rows[1][1+0*8],
					//NE	=rows[1][1+1*8],
					//WW	=rows[0][1-2*8],
					//W	=rows[0][1-1*8],
					cb	=src->data[idx+2],
					cr	=src->data[idx+0];
				//int update=(2*(NN+WW)-(N+W+NW+NE)+4*(cb+cr))>>4;
				//update=CLAMP(-halfs[1], update, halfs[1]-1);
				int update=(cb+cr)>>2;
				int luma=src->data[idx+1]-update;//subtract update
				luma<<=32-src->depth[1];
				luma>>=32-src->depth[1];
				src->data[idx+1]=luma;
			}
#endif
			for(int kc0=0;kc0<src->nch;++kc0)
			{
				int kc=perm[kc0], pred;
				int
					NNWW	=rows[2][kc-2*8+0],
					NNW	=rows[2][kc-1*8+0],
					NN	=rows[2][kc+0*8+0],
					NNE	=rows[2][kc+1*8+0],
					NNEE	=rows[2][kc+2*8+0],
					NWW	=rows[1][kc-2*8+0],
					NW	=rows[1][kc-1*8+0],
					N	=rows[1][kc+0*8+0],
					NE	=rows[1][kc+1*8+0],
					NEE	=rows[1][kc+2*8+0],
					WW	=rows[0][kc-2*8+0],
					W	=rows[0][kc-1*8+0],
					eNNWW	=rows[2][kc-2*8+4],
					eNNW	=rows[2][kc-1*8+4],
					eNN	=rows[2][kc+0*8+4],
					eNNE	=rows[2][kc+1*8+4],
					eNNEE	=rows[2][kc+2*8+4],
					eNWW	=rows[1][kc-2*8+4],
					eNW	=rows[1][kc-1*8+4],
					eN	=rows[1][kc+0*8+4],
					eNE	=rows[1][kc+1*8+4],
					eNEE	=rows[1][kc+2*8+4],
					eWW	=rows[0][kc-2*8+4],
					eW	=rows[0][kc-1*8+4],
					offset	=0;
				//(void)NNWW;	(void)eNNWW;
				//(void)NNW;	(void)eNNW;
				//(void)NN;	(void)eNN;
				//(void)NNE;	(void)eNNE;
				//(void)NNEE;	(void)eNNEE;
				//(void)NWW;	(void)eNWW;
				//(void)NW;	(void)eNW;
				//(void)N;	(void)eN;
				//(void)NE;	(void)eNE;
				//(void)NEE;	(void)eNEE;
				//(void)WW;	(void)eWW;
				//(void)W;	(void)eW;
				//if(kx==10&&ky==10)//
				//	printf("");
				//if(kc==1)//luma
				//	pred=(4*(N+W)+NE-NW)>>3;
				//else
				{
					static const int update_sh[]={2, 2, 2, 2};
					int
						vx=abs(eW-eWW)+abs(eN-eNW)+abs(eNE-eN)+1,
						vy=abs(eW-eNW)+abs(eN-eNN)+abs(eNE-eNNE)+1;
					int update=vx<vy?eN:eW;
					update/=4;
					//update=(update+(((1<<update_sh[kc])-1)&update>>31))>>update_sh[kc];
					//update=(update>>update_sh[kc])+(update>>31&((update&((1<<update_sh[kc])-1))!=0));
					//update/=abs(eN)+abs(eW)+1;

					//int update2=(update+(((1<<2)-1)&update>>31))>>2;
					//int update2=(update>>update_sh[kc])+(update>>31&((update&3)!=0));
					//update/=4;
					//if(update!=update2)
					//	LOG_ERROR("");//ceil(floor(5.13)*-1.3) = -6

					MEDIAN3_32(pred, N, W, N+W-NW);
					//pred=(N*vx+W*vy)/(vx+vy);//X
					//pred=(N+W)>>1;//X
					//pred=W;//X
					pred+=update;
					//CLAMP2_32(pred, pred, -halfs[kc], halfs[kc]-1);//X

					//if(kc==1)
					//	pred=((512-4)*pred+N+W+NW+NE)>>9;
				}

				//pred=(252*pred+N+W+NW+NE)>>8;
				//pred=((512-48)*pred+16*(N+W)+3*(NW+NE+NN+WW)+2*(NNW+NNE+NWW+NEE+NNWW+NNEE))>>9;
				
				if(kc0>0)
				{
					offset=rows[0][1];
					//if(kc0>1)
					//	offset=(2*offset+rows[0][2])>>1;
					pred+=offset;
					CLAMP2_32(pred, pred, -halfs[kc], halfs[kc]-1);
				}
				{
					int curr=src->data[idx+kc];
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=curr;
				
					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];

					src->data[idx+kc]=pred;
					rows[0][kc  ]=(fwd?curr:pred)-offset;
					rows[0][kc+4]=(fwd?pred:curr);
				}
			}
#ifdef PU_UPDATE_LUMA
			if(fwd)
			{
				int
					//NN	=rows[2][1+0*8],
					//NW	=rows[1][1-1*8],
					//N	=rows[1][1+0*8],
					//NE	=rows[1][1+1*8],
					//WW	=rows[0][1-2*8],
					//W	=rows[0][1-1*8],
					cb	=src->data[idx+2],
					cr	=src->data[idx+0];
				//int update=(2*(NN+WW)-(N+W+NW+NE)+4*(cb+cr))>>4;
				//update=CLAMP(-halfs[1], update, halfs[1]-1);
				int update=(cb+cr)>>2;
				{
					int luma=src->data[idx+1]+update;//add update
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					src->data[idx+1]=luma;
				}
				{
					int *ptr=src->data+idx, temp;
					ROTATE3(ptr[0], ptr[1], ptr[2], temp);
				}
			}
#endif
#endif
#if 0
			int
				*NNWW	=rows[2]-2*8,
				*NNW	=rows[2]-1*8,
				*NN	=rows[2]+0*8,
				*NNE	=rows[2]+1*8,
				*NNEE	=rows[2]+2*8,
				*NWW	=rows[1]-2*8,
				*NW	=rows[1]-1*8,
				*N	=rows[1]+0*8,
				*NE	=rows[1]+1*8,
				*NEE	=rows[1]+2*8,
				*WW	=rows[0]-2*8,
				*W	=rows[0]-1*8;
			for(int kc=0;kc<src->nch;++kc)
			{
				int
					vx=abs(W[kc+4]-WW[kc+4])+abs(N[kc+4]-NW[kc+4])+abs(NE[kc+4]-N[kc+4])+1,
					vy=abs(W[kc+4]-NW[kc+4])+abs(N[kc+4]-NN[kc+4])+abs(NE[kc+4]-NNE[kc+4])+1;
				MEDIAN3_32(preds[kc], N[kc], W[kc], N[kc]+W[kc]-NW[kc]);
				preds[kc]+=(vx<vy?N[kc+4]:W[kc+4])/4;
			}
			if(src->nch>=3)
			{
				preds[1]-=preds[2]>>1;
				preds[2]+=preds[1];
				preds[1]-=preds[0]>>1;
				preds[0]+=preds[1];
				{
					__m128i mpred=_mm_load_si128((__m128i*)preds);
					mpred=_mm_max_epi32(mpred, vmin);
					mpred=_mm_min_epi32(mpred, vmax);
					_mm_store_si128((__m128i*)preds, mpred);
				}
			}
			for(int kc=0;kc<src->nch;++kc)
			{
				int curr=src->data[idx+kc], pred=preds[kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=curr;
				
				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];

				src->data[idx+kc]=pred;
				rows[0][kc  ]=(fwd?curr:pred);
				rows[0][kc+4]=(fwd?pred:curr);
			}
			if(src->nch>=3)
			{
				int *curr=rows[0];
				curr[0]-=curr[1];
				curr[1]+=curr[0]>>1;
				curr[2]-=curr[1];
				curr[1]+=curr[2]>>1;
			}
#endif
			rows[0]+=8;
			rows[1]+=8;
			rows[2]+=8;
			rows[3]+=8;
		}
	}
	free(pixels);
}

#define DFWP_NPREDS 8
void pred_divfreeWP(Image *src, int fwd)//not DIV-free
{
	int nch;
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int perm[]={1, 2, 0, 3};
	int fwdmask=-fwd;
	//int weights[4][DFWP_NPREDS*2]={0};

	size_t bufsize=(src->iw+16LL)*sizeof(int[4*4*DFWP_NPREDS*2]);//4 padded rows * 4 channels max * {pixel, error, 2 perrors, 4 perrors, 8 perrors}
	int *pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	//FILLMEM((int*)weights, 0x8000, sizeof(weights), sizeof(int));
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8)*DFWP_NPREDS*2*4,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8)*DFWP_NPREDS*2*4,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8)*DFWP_NPREDS*2*4,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8)*DFWP_NPREDS*2*4,
		};
		int preds[DFWP_NPREDS*2]={0};
		int weights[DFWP_NPREDS*2]={0};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(int kc=0;kc<src->nch;++kc)
			{
				int kc2=kc*DFWP_NPREDS*2;
				int
					*NNN	=rows[3]+kc2+0*DFWP_NPREDS*2*4,
					*NNWW	=rows[2]+kc2-2*DFWP_NPREDS*2*4,
					*NNW	=rows[2]+kc2-1*DFWP_NPREDS*2*4,
					*NN	=rows[2]+kc2+0*DFWP_NPREDS*2*4,
					*NNE	=rows[2]+kc2+1*DFWP_NPREDS*2*4,
					*NNEE	=rows[2]+kc2+2*DFWP_NPREDS*2*4,
					*NWW	=rows[1]+kc2-2*DFWP_NPREDS*2*4,
					*NW	=rows[1]+kc2-1*DFWP_NPREDS*2*4,
					*N	=rows[1]+kc2+0*DFWP_NPREDS*2*4,
					*NE	=rows[1]+kc2+1*DFWP_NPREDS*2*4,
					*NEE	=rows[1]+kc2+2*DFWP_NPREDS*2*4,
					*NEEE	=rows[1]+kc2+3*DFWP_NPREDS*2*4,
					*WWW	=rows[0]+kc2-3*DFWP_NPREDS*2*4,
					*WW	=rows[0]+kc2-2*DFWP_NPREDS*2*4,
					*W	=rows[0]+kc2-1*DFWP_NPREDS*2*4,
					*curr	=rows[0]+kc2+0*DFWP_NPREDS*2*4;
				int pred=0;
				//int *curr_weights=weights[kc];
				int j=DFWP_NPREDS;

				//if(ky==100&&kx==100)//
				//	printf("");

				preds[j++]=N[0];
				preds[j++]=W[0];
				preds[j++]=3*(N[0]-NN[0])+NNN[0];
				preds[j++]=3*(W[0]-WW[0])+WWW[0];
				preds[j++]=N[0]+W[0]-NW[0];
				preds[j++]=N[0]+NE[0]-NNE[0];
				preds[j++]=W[0]+NE[0]-N[0];
				preds[j++]=(W[0]+NEEE[0])>>1;

				for(int nsrc=DFWP_NPREDS, srcidx=DFWP_NPREDS;nsrc>1;)
				{
					int ndst=nsrc>>1, dstidx=srcidx-ndst;
					for(int kp=0;kp<ndst;++kp)
					{
						int w1=NNWW[srcidx+kp*2+0]+NNW[srcidx+kp*2+0]+NN[srcidx+kp*2+0]+NNE[srcidx+kp*2+0]+NNEE[srcidx+kp*2+0]+NWW[srcidx+kp*2+0]+NW[srcidx+kp*2+0]+N[srcidx+kp*2+0]*3+NE[srcidx+kp*2+0]+NEE[srcidx+kp*2+0]+WW[srcidx+kp*2+0]+W[srcidx+kp*2+0]*3+1;
						int w2=NNWW[srcidx+kp*2+1]+NNW[srcidx+kp*2+1]+NN[srcidx+kp*2+1]+NNE[srcidx+kp*2+1]+NNEE[srcidx+kp*2+1]+NWW[srcidx+kp*2+1]+NW[srcidx+kp*2+1]+N[srcidx+kp*2+1]*3+NE[srcidx+kp*2+1]+NEE[srcidx+kp*2+1]+WW[srcidx+kp*2+1]+W[srcidx+kp*2+1]*3+1;

						//weights[dstidx+kp]=(w2<<16)/(w1+w2);
						//preds[dstidx+kp]=preds[srcidx+kp*2+1]+((preds[srcidx+kp*2+0]-preds[srcidx+kp*2+1])*weights[dstidx+kp]>>16);

						preds[dstidx+kp]=(preds[srcidx+kp*2+0]*w2+preds[srcidx+kp*2+1]*w1)/(w1+w2);//NPREDs DIVs
					}
					nsrc=ndst;
					srcidx=dstidx;
				}
				pred=preds[1];
				
				{
					int c2=src->data[idx+kc];
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=c2;
				
					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];

					src->data[idx+kc]=pred;
					curr[0]=fwd?c2:pred;
					curr[1]=fwd?pred:c2;
				}
				for(int k=2;k<DFWP_NPREDS*2;++k)
					curr[k]=abs(curr[0]-preds[k]);
			}
			rows[0]+=DFWP_NPREDS*2*4;
			rows[1]+=DFWP_NPREDS*2*4;
			rows[2]+=DFWP_NPREDS*2*4;
			rows[3]+=DFWP_NPREDS*2*4;
		}
	}
	free(pixels);
}


void pred_nblic(Image *src, int fwd)//https://github.com/WangXuan95/NBLIC-Image-Compression
{
	int nch;
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int fwdmask=-fwd;
	static const int thresholds[]=
	{
		1, 3, 9, 20, 50, 110, 300, 800,
	};
	const static int q_mid[]={0, 2, 4, 7, 10, 14, 20, 26, 34, 42, 52, 64, 78, 95, 135, 200};
	int bufsize=(src->iw+16LL)*sizeof(int[4*4*2]);//4 padded rows * 4 channels max
	int *pixels=(int*)malloc(bufsize);
	int ssesize=(int)sizeof(int[4<<4>>1<<8]);
	int *sse=(int*)malloc(ssesize);
	if(!pixels||!sse)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	memset(sse, 0, ssesize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(int ky=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<src->nch;++kc)
			{
				//NNWW NNW NN NNE NNEE
				//NWW  NW  N  NE  NEE
				//WW   W   ?
				int
					idx	=(src->iw*ky+kx)<<2|kc,
					NNWW	=rows[2][kc-2*4*2+0],
					NNW	=rows[2][kc-1*4*2+0],
					NN	=rows[2][kc+0*4*2+0],
					NNE	=rows[2][kc+1*4*2+0],
					NNEE	=rows[2][kc+2*4*2+0],
					NWW	=rows[1][kc-2*4*2+0],
					NW	=rows[1][kc-1*4*2+0],
					N	=rows[1][kc+0*4*2+0],
					NE	=rows[1][kc+1*4*2+0],
					NEE	=rows[1][kc+2*4*2+0],
					WW	=rows[0][kc-2*4*2+0],
					W	=rows[0][kc-1*4*2+0],
					*curr	=rows[0]+kc+0*4*2+0,
					eN	=rows[1][kc+0*4*2+4],
					eW	=rows[0][kc-1*4*2+4];

				//simple predict:

				//      -1
				//   -2  9  2
				//-1  9  ?
				int pred_linear=9*(N+W)+2*(NE-NW)-NN-WW;
				int pred_ang=0, cost, csum=0, cmin=0xFFFFFF, wt;
				CLAMP2_32(pred_linear, pred_linear, -halfs[kc]<<4, halfs[kc]<<4);

				cost=2*(abs(W-WW)+abs(NW-NWW)+abs(N-NW)+abs(NE-N));//180 degrees
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=2*W;
				}

				cost=2*(abs(W-NW)+abs(NW-NNW)+abs(N-NN)+abs(NE-NNE));//90 degrees
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=2*N;
				}

				cost=2*(abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW)+abs(NE-NN));//135 degrees
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=2*NW;
				}

				cost=2*(abs(W-N)+abs(NW-NN)+abs(N-NNE)+abs(NE-NNEE));//45 degrees
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=2*NE;
				}

				cost=abs(2*W-WW-NWW)+abs(2*NW-NNW-NNWW)+abs(2*N-NW-NNW)+abs(2*NE-N-NN);//180-atand(1/2) ~= 153 degrees
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=W+NW;
				}

				cost=abs(2*W-NWW-NW)+abs(2*NW-NNWW-NNW)+abs(2*N-NNW-NN)+abs(2*NE-NN-NNE);//90+atand(1/2) ~= 116 degrees
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=NW+N;
				}

				cost=abs(2*W-NW-N)+abs(2*NW-NNW-NN)+abs(2*N-NN-NNE)+abs(2*NE-NNE-NNEE);//90-atand(1/2) ~= 63 degrees
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=N+NE;
				}
				
				csum-=7*cmin;
				for(wt=0;wt<_countof(thresholds)&&(thresholds[wt]<<src->depth[kc])-thresholds[wt]<=(csum<<3);++wt);
				int pred=(8*wt*pred_ang+(8-wt)*pred_linear+64)>>7, pred0=pred;

				//get quantized delta:
				int delta=abs(W-WW)+abs(N-NW)+abs(N-NE)+abs(W-NW)+abs(N-NN)+abs(NE-NNE)+2*abs(eW), qd, qu, qv, qw;
				for(qd=0;qd<_countof(q_mid)-1&&delta>q_mid[qd];++qd);
				qu=qv=qd;
				qw=0;
				if(delta<q_mid[qd])
				{
					qw=((delta-q_mid[qd])<<5)/(q_mid[qd]-q_mid[qd-1]);
					if(qw<(1<<5>>1))
						qu=qd-1;
					else
					{
						qv=qd-1;
						qw=(1<<5)-qw;
					}
				}

				//get context address:
				int ctx=kc;
				ctx=ctx<<(4-1)|qu>>1;
				ctx=ctx<<1|(pred>2*N-NN);
				ctx=ctx<<1|(pred>2*W-WW);
				ctx=ctx<<1|(pred>NN);
				ctx=ctx<<1|(pred>WW);
				ctx=ctx<<1|(pred>NE);
				ctx=ctx<<1|(pred>NW);
				ctx=ctx<<1|(pred>N);
				ctx=ctx<<1|(pred>W);

				//correct prediction by context:
				int *cell=sse+ctx;
				int corr=*cell;
				int ssesign=corr>>7&1;
				pred+=corr>>8;

				CLAMP2_32(pred, pred, -halfs[kc], halfs[kc]);
				
				//apply prediction:
				int p2=pred;
				int c2=src->data[idx];
				p2^=fwdmask;
				p2-=fwdmask;
				p2+=c2;
				p2<<=32-src->depth[kc];
				p2>>=32-src->depth[kc];
				src->data[idx]=p2;
				rows[0][kc+0]=fwd?c2:p2;
				rows[0][kc+4]=rows[0][kc+0]-pred0;//FIXME error uses prediction before SSE correction?

				//update context:
				int e2=rows[0][kc+4];
				*cell=(corr*((1<<7)-1)+(e2<<7)+(1<<7>>1))>>7;//curr = curr*127/128+error	FIXME compare with traditional SSE
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	free(pixels);
	free(sse);
}


static void sort_int32(int *data, int count)
{
	int temp=0;
	for(int k=0;k<count-1;++k)
	{
		for(int k2=k+1;k2<count;++k2)
		{
			if(data[k2]<data[k])
				SWAPVAR(data[k], data[k2], temp);
		}
	}
}

	#define WGRAD_UPDATE_LUMA
//	#define WGRAD_SSE
#ifdef WGRAD_SSE
static int quantize_nb(int x)
{
	int ax=abs(x), negmask=x<0;
	x=FLOOR_LOG2(ax);
	x^=negmask;
	x-=negmask;
	return x;
}
#define QUANTIZE_NB(X) quantize_nb(X)+maxdepth
#endif
void pred_wgrad(Image *src, int fwd, int hasRCT)
{
	//static int sh=8;
	//if(fwd)
	//{
	//	--sh;
	//	if(sh<0)
	//		sh=8;
	//}
	int nch;
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	const int perm0[]={1, 2, 0, 3, 0, 1, 2, 3}, *perm=perm0+hasRCT*4;
	int fwdmask=-fwd;

	int *pixels=(int*)malloc((src->iw+16LL)*sizeof(int[4*4]));//4 padded rows * 4 channels max
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, (src->iw+16LL)*sizeof(int[4*4]));
#ifdef WGRAD_SSE
	int maxdepth=MAXVAR(src->depth[0], src->depth[1]);
	UPDATE_MAX(maxdepth, src->depth[2]);
	UPDATE_MAX(maxdepth, src->depth[3]);
	int clevels=maxdepth<<1|1;
	size_t ssesize=sizeof(int[4*2])*clevels*clevels;//4 channels max * {sum, count} * clevels^2
	int *sse=(int*)malloc(ssesize);
	if(!sse)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(sse, 0, ssesize);
#endif
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int eW[16]={0};
		int *rows[]=
		{
			pixels+(((src->iw+16LL)*((ky-0LL)&3)+8LL)<<2),
			pixels+(((src->iw+16LL)*((ky-1LL)&3)+8LL)<<2),
			pixels+(((src->iw+16LL)*((ky-2LL)&3)+8LL)<<2),
			pixels+(((src->iw+16LL)*((ky-3LL)&3)+8LL)<<2),
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
#ifdef WGRAD_UPDATE_LUMA
			if(!hasRCT&&!fwd)
			{
				int
					cb	=src->data[idx+2],
					cr	=src->data[idx+0];
				int update=(cb+cr)>>2;
				update=CLAMP(-halfs[1], update, halfs[1]-1);
				{
					int luma=src->data[idx+1]-update;//subtract update
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					src->data[idx+1]=luma;
				}
			}
#endif
			for(int kc0=0;kc0<src->nch;++kc0)
			{
				int kc=perm[kc0], pred;
				int
					NNNN	=rows[0][kc+0*4],
					NNN	=rows[3][kc+0*4],
					NNWW	=rows[2][kc-2*4],
					NNW	=rows[2][kc-1*4],
					NN	=rows[2][kc+0*4],
					NNE	=rows[2][kc+1*4],
					NNEE	=rows[2][kc+2*4],
					NWW	=rows[1][kc-2*4],
					NW	=rows[1][kc-1*4],
					N	=rows[1][kc+0*4],
					NE	=rows[1][kc+1*4],
					NEE	=rows[1][kc+2*4],
					NEEE	=rows[1][kc+3*4],
					NEEEE	=rows[1][kc+4*4],
					NEEEEE	=rows[1][kc+5*4],
					WWWWW	=rows[0][kc-5*4],
					WWWW	=rows[0][kc-4*4],
					WWW	=rows[0][kc-3*4],
					WW	=rows[0][kc-2*4],
					W	=rows[0][kc-1*4],
					offset	=0;
				(void)NNNN;		//<- what insanity looks like
				(void)NNN;
				(void)NNWW;
				(void)NNW;
				(void)NN;
				(void)NNE;
				(void)NNEE;
				(void)NWW;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)NEEEE;
				(void)NEEEEE;
				(void)WWW;
				(void)WW;
				(void)W;
				if(!hasRCT&&kc0>0)
				{
					offset+=rows[0][1];
					if(kc0>1)
						offset=(2*offset+rows[0][2])>>1;
				}
				//int nb[]={N, W, NW, NE, NEE};
				//sort_int32(nb, (int)_countof(nb));
				//int pred=nb[_countof(nb)>>1];
				//pred=(pred+N+W-NW)>>1;
#if 1
				long long lpred, wsum;
				int
					gy=abs(N-NN)+abs(W-NW)+1,
					gx=abs(W-WW)+abs(N-NW)+1;
				int preds[]=
				{
					//(N*gx+W*gy)/(gy+gx),
					N+W-NW,
					W+NE-N,
					N+NE-NNE,
					N,
					W,
					3*(N-NN)+NNN,
					3*(W-WW)+WWW,
					(W+NEEE)>>1,

					//(N+NE)>>1,
					//W+((WWWWW-WWWW+NEEEE-NEEEEE)>>2),
					//(W+N+NE+NEE+NEEE)/5,
					//(N+2*NE-NNE)>>1,
				};
				//MEDIAN3_32(preds[0], N, W, preds[0]);

				lpred=0;
				wsum=0;
				for(int k=0;k<_countof(preds);++k)
				{
					int weight=0x1000000/(eW[k]+1);
					lpred+=(long long)weight*preds[k];
					wsum+=weight;
				}
				lpred/=wsum+1;
				pred=(int)lpred;
				//offset=0;
#endif
#if 0
				int cgrad;
				MEDIAN3_32(cgrad, N, W, N+W-NW);
				pred=(eW[1]*W+eW[0]*cgrad)/(eW[0]+eW[1]);
#endif
#if 0
				{
					int
						gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1,
						gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1,
						g135=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW)+abs(NE-NN)+1,
						g45=abs(W-N)+abs(N-NNE)+abs(NE-NNEE)+abs(WW-NW)+1;
					//int gx=abs(W-WW)+1, gy=abs(N-NN)+1;//X
					(void)gx;
					(void)gy;
					(void)g135;
					(void)g45;

					//NW, N, NE, W
					//if(NW>N)
					//{
					//	if(N>NE)
					//	{
					//		if(W>N)
					//	}
					//	else
					//	{
					//	}
					//}
					//else
					//{
					//	if(NW>NE)
					//	{
					//	}
					//	else
					//	{
					//	}
					//}
					
					pred=(int)(((long long)gx*N+(long long)gy*W)/(gx+gy));		//<- NOT better than CG: LPCB, DSLR2
					
					//pred=(gx*N+gy*W)/((gx+gy)+((gx+gy)>>2));

					//if(kx==10&&ky==10)
					//	printf("");
					//pred=(NW+((2*N-NW)*gy+(2*W-NW)*gx)/(gx+gy))>>1;
					//pred=NW+((N-NW)*gy+(W-NW)*gx)/(gx+gy);
					//MEDIAN3_32(pred, N, W, pred);

					//MEDIAN3_32(pred, N, W, N+W-NW);
				}
#endif
				//const int sh=8;
				//int grad=N+W-NW, pred;
				//MEDIAN3_32(pred, N, W, grad);
				//pred=(pred<<sh)-pred+grad;
				//pred=(pred>>sh)+(pred<0);
				
#if 0
				int sh=src->depth[kc]-8;
				int diff=(gy-gx)>>sh, diff2=(g45-g135)>>sh, diff3=NE-NW;
				int paper_GAP;
				if(gy+gx>32)
					paper_GAP=(int)(((long long)gx*N+(long long)gy*W)/((long long)gy+gx));
				else if(diff>12)
					paper_GAP=(N+2*W)/3;
				else if(diff<-12)
					paper_GAP=(2*N+W)/3;
				else
					paper_GAP=(N+W)>>1;

				if(diff2>32)
					paper_GAP+=diff3>>2;
				else if(diff2>16)
					paper_GAP+=diff3*3>>4;
				else if(diff2>=-16)
					paper_GAP+=diff3>>3;
				else if(diff2>=-32)
					paper_GAP+=diff3>>4;
				int pred=paper_GAP;
#endif
				//int pred=(3*(N+W)-2*(NN+WW))>>1;
				//pred=(pred+NE+NEE+NNE-NW+NWW+WW)/5;
				//int pred=(gx*N+gy*W+g135*((3*NW+NNWW)>>2)+g45*((3*NE+NNEE)>>2))/(gx+gy+g45+g135);
				//int
				//	w45=0x10000/g45,
				//	w90=0x10000/gy,
				//	w135=0x10000/g135,
				//	w180=0x10000/gx,
				//	wsum=w45+w90+w135+w180;
				//wsum+=!wsum;
				//int pred=(w45*NE+w90*N+w135*NW+w180*W)/wsum;
				//int pred=(gx*(2*N-NN)+gy*(2*W-WW))/(gx+gy);
				//int pred=(W+WW+NEE-NW+3*(N-NN)+NNN+NE)>>2;
				//int pred=MEDIAN3(N, W, NW);
#ifdef WGRAD_SSE
				int sse_idx=kc;
				{
					int x=N;
					int ax=abs(x), negmask=x<0;
					x=FLOOR_LOG2(ax);
					x^=negmask;
					x-=negmask;
					x+=maxdepth;
					sse_idx*=clevels;
					sse_idx+=x;
				}
				{
					int x=W;
					int ax=abs(x), negmask=x<0;
					x=FLOOR_LOG2(ax);
					x^=negmask;
					x-=negmask;
					x+=maxdepth;
					sse_idx*=clevels;
					sse_idx+=x;
				}
				int *curr_sse=sse+((size_t)sse_idx<<1);
				pred+=curr_sse[0]/(curr_sse[1]+1);
#endif

				pred+=offset;
				pred=CLAMP(-halfs[kc], pred, halfs[kc]-1);
				{
					int curr=src->data[idx+kc];
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=curr;

					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];

					src->data[idx+kc]=pred;
					rows[0][kc]=(fwd?curr:pred)-offset;
#ifdef WGRAD_SSE
					int error=fwd?pred:curr;
					curr_sse[0]+=error;
					++curr_sse[1];
					if(curr_sse[1]>=640)
					{
						curr_sse[0]>>=1;
						curr_sse[1]>>=1;
					}
#endif
				}
#if 1
				{
					int curr=rows[0][kc];
					for(int k=0;k<_countof(preds);++k)
					{
						int err=abs(curr-preds[k]);
						eW[k]+=err;
						eW[k]-=(eW[k]+3)>>2;
						//eW[k]*=3;
						//eW[k]>>=2;
					}
				}
#endif
#if 0
				{
					int curr=rows[0][kc];
					eW[0]=(eW[0]+abs(curr-W)+1)>>1;
					eW[1]=(eW[1]+abs(curr-cgrad)+1)>>1;
					//eW[0]=abs(curr-W);
					//eW[1]=abs(curr-cgrad);
				}
#endif
			}
#ifdef WGRAD_UPDATE_LUMA
			if(!hasRCT&&fwd)
			{
				int
					cb	=src->data[idx+2],
					cr	=src->data[idx+0];
				int update=(cb+cr)>>2;
				update=CLAMP(-halfs[1], update, halfs[1]-1);
				{
					int luma=src->data[idx+1]+update;//add update
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					src->data[idx+1]=luma;
				}
			}
#endif
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	free(pixels);
	//set_window_title("sh=%d", sh);
}

	#define WG2_USE_FP64
//	#define WG2_DIVFREE
//	#define WG2_NPREDS 8
void pred_wgrad2(Image *src, int fwd)
{
	int nch;
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int perm[]={1, 2, 0, 3};
	int fwdmask=-fwd;

	int *pixels;
	//int maxdepth, histsize, *hist;

	pixels=(int*)malloc((src->iw+16LL)*sizeof(int[4*4]));//4 padded rows * 4 channels max
	//maxdepth=MAXVAR(src->depth[0], src->depth[1]);
	//UPDATE_MAX(maxdepth, src->depth[2]);
	//UPDATE_MAX(maxdepth, src->depth[3]);
	//histsize=sizeof(int[4])<<maxdepth;
	//hist=(int*)malloc(histsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, (src->iw+16LL)*sizeof(int[4*4]));
	//memset(hist, 0, histsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int eprev[16*4]={0};
		int *rows[]=
		{
			pixels+(((src->iw+16LL)*((ky-0LL)&3)+8LL)<<2),
			pixels+(((src->iw+16LL)*((ky-1LL)&3)+8LL)<<2),
			pixels+(((src->iw+16LL)*((ky-2LL)&3)+8LL)<<2),
			pixels+(((src->iw+16LL)*((ky-3LL)&3)+8LL)<<2),
		};
#ifdef WG2_DIVFREE
		int mixinCDF[WG2_NPREDS]={0}, tempCDF[WG2_NPREDS+1]={0};
		for(int k=0;k<WG2_NPREDS+1;++k)
		{
			int val=0x10000*k/WG2_NPREDS;
			eprev[k<<2|0]=val;
			eprev[k<<2|1]=val;
			eprev[k<<2|2]=val;
			eprev[k<<2|3]=val;
		}
		for(int k=0;k<WG2_NPREDS;++k)
			mixinCDF[k]=(0x10000-16)*k/(WG2_NPREDS-1);
#endif
//#ifdef WG2_DIVFREE
//		eprev[0]=2<<src->depth[0];
//		eprev[1]=2<<src->depth[1];
//		eprev[2]=2<<src->depth[2];
//		eprev[3]=2<<src->depth[3];
//		memfill(eprev+4, eprev, sizeof(eprev)-sizeof(int[4]), sizeof(int[4]));
//#endif
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
//#ifdef WGRAD_UPDATE_LUMA
//			if(!fwd)
//			{
//				int
//					cb	=src->data[idx+2],
//					cr	=src->data[idx+0];
//				int update=(cb+cr)>>2;
//				update=CLAMP(-halfs[1], update, halfs[1]-1);
//				{
//					int luma=src->data[idx+1]-update;//subtract update
//					luma<<=32-src->depth[1];
//					luma>>=32-src->depth[1];
//					src->data[idx+1]=luma;
//				}
//			}
//#endif
			for(int kc0=0;kc0<src->nch;++kc0)
			{
				int kc=perm[kc0], pred;
				int
					NNNN	=rows[0][kc+0*4],
					NNN	=rows[3][kc+0*4],
					NNWW	=rows[2][kc-2*4],
					NNW	=rows[2][kc-1*4],
					NN	=rows[2][kc+0*4],
					NNE	=rows[2][kc+1*4],
					NNEE	=rows[2][kc+2*4],
					NWW	=rows[1][kc-2*4],
					NW	=rows[1][kc-1*4],
					N	=rows[1][kc+0*4],
					NE	=rows[1][kc+1*4],
					NEE	=rows[1][kc+2*4],
					NEEE	=rows[1][kc+3*4],
					NEEEE	=rows[1][kc+4*4],
					NEEEEE	=rows[1][kc+5*4],
					WWWWW	=rows[0][kc-5*4],
					WWWW	=rows[0][kc-4*4],
					WWW	=rows[0][kc-3*4],
					WW	=rows[0][kc-2*4],
					W	=rows[0][kc-1*4],
					offset	=0;
				//int *curr_hist=hist+((size_t)kc<<maxdepth);
				//(void)NNNN;
				//(void)NNN;
				//(void)NNWW;
				//(void)NNW;
				//(void)NN;
				//(void)NNE;
				//(void)NNEE;
				//(void)NWW;
				//(void)NW;
				//(void)N;
				//(void)NE;
				//(void)NEE;
				//(void)NEEE;
				//(void)NEEEE;
				//(void)NEEEEE;
				//(void)WWW;
				//(void)WW;
				//(void)W;
				//if(kc0>0)
				//{
				//	offset+=rows[0][1];
				//	if(kc0>1)
				//		offset=(2*offset+rows[0][2])>>1;
				//}
				long long lpred, wsum;
				int preds[]=
				{
					N+W-NW,
					W+NE-N,
					N+NE-NNE,
					N,
					W,
					3*(N-NN)+NNN,
					3*(W-WW)+WWW,
					(W+NEEE)>>1,
				};

				lpred=0;
				wsum=0;
#ifdef WG2_DIVFREE
				for(int k=0;k<_countof(preds);++k)
				{
					int weight=eprev[(k+1)<<2|kc]-eprev[k<<2|kc];
					lpred+=(long long)weight*preds[k];
				}
				lpred>>=16;
#elif defined WG2_USE_FP64
				{
					double fpred=0, fsum=0;
					for(int k=0;k<_countof(preds);++k)
					{
						double weight=1./(eprev[k<<2|kc]+1);
						fpred+=weight*preds[k];
						fsum+=weight;
					}
					fpred/=fsum;
					pred=_cvt_dtoi_fast(fpred);
				}
#else
				for(int k=0;k<_countof(preds);++k)
				{
//#ifdef WG2_DIVFREE
//					int weight=eprev[k<<2|kc];
//#else
					int weight=0x1000000/(eprev[k<<2|kc]+1);
//#endif
					lpred+=(long long)weight*preds[k];
					wsum+=weight;
				}
				lpred/=wsum+1;
#endif
#ifndef WG2_USE_FP64
				pred=(int)lpred;
#endif
#if 0
				{
					int val=pred+halfs[kc];
					int start=val, end=val;
					//if(ky==10&&kx==10)
					//	printf("");

					for(int k=val;k>=0;--k)
					{
						if(curr_hist[k])
						{
							start=k;
							break;
						}
					}
					for(int k=val;k<nlevels[kc];++k)
					{
						if(curr_hist[k])
						{
							end=k;
							break;
						}
					}
					if(val-start<end-val)
						pred=start-halfs[kc];
					else if(val-start>end-val)
						pred=end-halfs[kc];
					
					//const int windowsize=10;
					//int val=pred+halfs[kc], start, end;
					//start=MAXVAR(val-(windowsize>>1), 0);
					//end=MINVAR(val+(windowsize>>1), nlevels[kc]-1);
					//for(int k=start;k<end;++k)
					//{
					//	if(curr_hist[val]<curr_hist[k])
					//		val=k;
					//}
					//pred=val-halfs[kc];
				}
#endif
				//pred+=offset;
				//pred=CLAMP(-halfs[kc], pred, halfs[kc]-1);
				{
					int curr=src->data[idx+kc];
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=curr;

					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];

					src->data[idx+kc]=pred;
					rows[0][kc]=(fwd?curr:pred)-offset;
				}
				{
					int curr=rows[0][kc];
#ifdef WG2_DIVFREE
					memset(tempCDF, 0, sizeof(tempCDF));
					for(int k=0;k<_countof(preds);++k)
					{
						int err=abs(curr-preds[k]);
						int *ptr=mixinCDF;
						for(int k2=0;k2<_countof(preds);++k2)
						{
							tempCDF[k2]+=(*ptr-tempCDF[k2])*err>>src->depth[kc];
							ptr+=k2!=k;
						}
					}
					for(int k=0;k<_countof(preds);++k)
					{
						int update=tempCDF[k]>>3;
						eprev[k<<2|kc]=update;
						//eprev[k<<2|kc]+=(update-eprev[k<<2|kc])>>1;
					}
#else
//#ifdef WG2_DIVFREE
//					int minerr=0;
//#endif
					for(int k=0;k<_countof(preds);++k)
					{
						int err=abs(curr-preds[k]);
//#ifdef WG2_DIVFREE
//						eprev[k<<2|kc]-=err;
//						if(eprev[k<<2|kc]<1)
//							eprev[k<<2|kc]=1;
//						if(!k||minerr>eprev[k<<2|kc])
//							minerr=eprev[k<<2|kc];
//#else
						eprev[k<<2|kc]+=err;
						eprev[k<<2|kc]-=(eprev[k<<2|kc]+3)>>2;
						//eprev[k]*=3;
						//eprev[k]>>=2;
//#endif
					}
//#ifdef WG2_DIVFREE
//					if(minerr<nlevels[kc])
//					{
//						for(int k=0;k<_countof(preds);++k)
//						{
//							if(eprev[k<<2|kc]>(0x7FFFFFFF>>8))
//								eprev[k<<2|kc]=0x7FFFFFFF>>8;
//							eprev[k<<2|kc]<<=8;
//						}
//					}
//#endif
#endif
				}
				//{
				//	int curr=rows[0][kc];
				//	curr+=halfs[kc];
				//	curr&=nlevels[kc]-1;
				//	++curr_hist[curr];
				//	if(curr_hist[curr]>6144)
				//	{
				//		for(int ks=0;ks<nlevels[kc];++ks)
				//			curr_hist[ks]>>=1;
				//	}
				//}
			}
//#ifdef WGRAD_UPDATE_LUMA
//			if(fwd)
//			{
//				int
//					cb	=src->data[idx+2],
//					cr	=src->data[idx+0];
//				int update=(cb+cr)>>2;
//				update=CLAMP(-halfs[1], update, halfs[1]-1);
//				{
//					int luma=src->data[idx+1]+update;//add update
//					luma<<=32-src->depth[1];
//					luma>>=32-src->depth[1];
//					src->data[idx+1]=luma;
//				}
//			}
//#endif
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	free(pixels);
	//free(hist);
}
void pred_wgrad3(Image *src, int fwd, int hasRCT)
{
	//static int callidx=0;

	int nch;
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	const int perm0[]={1, 2, 0, 3, 0, 1, 2, 3}, *perm=perm0+hasRCT*4;
	int fwdmask=-fwd;
	int *pixels;

	//++callidx;

	pixels=(int*)malloc((src->iw+16LL)*sizeof(int[4*4]));//4 padded rows * 4 channels max
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, (src->iw+16LL)*sizeof(int[4*4]));
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		double fprev[16*4]={0};
		int eprev[16*4]={0};
		int *rows[]=
		{
			pixels+(((src->iw+16LL)*((ky-0LL)&3)+8LL)<<2),
			pixels+(((src->iw+16LL)*((ky-1LL)&3)+8LL)<<2),
			pixels+(((src->iw+16LL)*((ky-2LL)&3)+8LL)<<2),
			pixels+(((src->iw+16LL)*((ky-3LL)&3)+8LL)<<2),
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			if(!hasRCT&&!fwd)
			{
				int
					cb	=src->data[idx+2],
					cr	=src->data[idx+0];
				int update=(cb+cr)>>2;
				update=CLAMP(-halfs[1], update, halfs[1]-1);
				{
					int luma=src->data[idx+1]-update;//subtract update
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					src->data[idx+1]=luma;
				}
			}
			for(int kc0=0;kc0<src->nch;++kc0)
			{
				int kc=perm[kc0], pred;
				int
					NNNN	=rows[0][kc+0*4],
					NNN	=rows[3][kc+0*4],
					NNWW	=rows[2][kc-2*4],
					NNW	=rows[2][kc-1*4],
					NN	=rows[2][kc+0*4],
					NNE	=rows[2][kc+1*4],
					NNEE	=rows[2][kc+2*4],
					NWW	=rows[1][kc-2*4],
					NW	=rows[1][kc-1*4],
					N	=rows[1][kc+0*4],
					NE	=rows[1][kc+1*4],
					NEE	=rows[1][kc+2*4],
					NEEE	=rows[1][kc+3*4],
					NEEEE	=rows[1][kc+4*4],
					NEEEEE	=rows[1][kc+5*4],
					WWWWW	=rows[0][kc-5*4],
					WWWW	=rows[0][kc-4*4],
					WWW	=rows[0][kc-3*4],
					WW	=rows[0][kc-2*4],
					W	=rows[0][kc-1*4],
					offset	=0;
				int preds[]=
				{
					N+W-NW,
					W+NE-N,
					N+NE-NNE,
					N,
					W,
					3*(N-NN)+NNN,
					3*(W-WW)+WWW,
					(W+NEEE)>>1,
				};
				double fpred=0, fsum=0;
				for(int k=0;k<_countof(preds);++k)
				{
					double weight=1./(eprev[k<<2|kc]+1);
					//double weight=1./(callidx&1?eprev[k<<2|kc]+1:fprev[k<<2|kc]+1);
					fpred+=weight*preds[k];
					fsum+=weight;
				}
				fpred/=fsum;
				pred=_cvt_dtoi_fast(fpred);
				if(!hasRCT&&kc0>0)
				{
					offset+=rows[0][1];
					if(kc0>1)
						offset=(2*offset+rows[0][2])>>1;
					pred+=offset;
					pred=CLAMP(-halfs[kc], pred, halfs[kc]-1);
				}
				{
					int curr=src->data[idx+kc];
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=curr;

					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];

					src->data[idx+kc]=pred;
					rows[0][kc]=(fwd?curr:pred)-offset;
				}
				{
					int curr=rows[0][kc];
					for(int k=0;k<_countof(preds);++k)
					{
						int err=abs(curr-preds[k]);
						//if(callidx&1)
						//{
							eprev[k<<2|kc]+=err;
							eprev[k<<2|kc]=eprev[k<<2|kc]*0xE6>>8;
							//eprev[k<<2|kc]-=(eprev[k<<2|kc]+3)>>2;
						//}
						//else
						//{
						//	fprev[k<<2|kc]+=err;
						//	fprev[k<<2|kc]*=0.9;
						//}
					}
				}
			}
			if(!hasRCT&&fwd)
			{
				int
					cb	=src->data[idx+2],
					cr	=src->data[idx+0];
				int update=(cb+cr)>>2;
				update=CLAMP(-halfs[1], update, halfs[1]-1);
				{
					int luma=src->data[idx+1]+update;//add update
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					src->data[idx+1]=luma;
				}
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	free(pixels);

	//messagebox(MBOX_OK, "Info", "%s", callidx&1?"int":"double");
}


#define LWAV_NPARAMS 5
//#define LWAV_LPARAMS 3
#define LWAV_PBITS 8
#define LWAV_NBITS 16
void pred_lwav(Image *src, int fwd)
{
	int nch;
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int fwdmask=-fwd;

	int bufsize=(src->iw+16LL)*sizeof(int[4*4*(2+LWAV_NPARAMS)]);
	int *pixels=(int*)malloc(bufsize);//4 padded rows * 4 channels max * {pixels, errors}
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	int params[4][LWAV_NPARAMS]={0};
	FILLMEM(params[0], 1<<LWAV_NBITS, sizeof(int[LWAV_NPARAMS]), sizeof(int));
	FILLMEM(params[1], 1<<LWAV_NBITS, sizeof(int[LWAV_NPARAMS]), sizeof(int));
	FILLMEM(params[2], 1<<LWAV_NBITS, sizeof(int[LWAV_NPARAMS]), sizeof(int));
	FILLMEM(params[3], 1<<LWAV_NBITS, sizeof(int[LWAV_NPARAMS]), sizeof(int));
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int eW[16]={0};
		int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8LL)*4*(2+LWAV_NPARAMS),
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8LL)*4*(2+LWAV_NPARAMS),
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8LL)*4*(2+LWAV_NPARAMS),
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8LL)*4*(2+LWAV_NPARAMS),
		};
		int pred=0;
		//for(int k=0;k<LWAV_NPARAMS*4;++k)
		//	((int*)params)[k]>>=1;
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(int kc=0;kc<src->nch;++kc)
			{
				int
					*NNNN	=rows[0]+(kc+0*4)*(2+LWAV_NPARAMS),
					*NNN	=rows[3]+(kc+0*4)*(2+LWAV_NPARAMS),
					*NNWW	=rows[2]+(kc-2*4)*(2+LWAV_NPARAMS),
					*NNW	=rows[2]+(kc-1*4)*(2+LWAV_NPARAMS),
					*NN	=rows[2]+(kc+0*4)*(2+LWAV_NPARAMS),
					*NNE	=rows[2]+(kc+1*4)*(2+LWAV_NPARAMS),
					*NNEE	=rows[2]+(kc+2*4)*(2+LWAV_NPARAMS),
					*NWW	=rows[1]+(kc-2*4)*(2+LWAV_NPARAMS),
					*NW	=rows[1]+(kc-1*4)*(2+LWAV_NPARAMS),
					*N	=rows[1]+(kc+0*4)*(2+LWAV_NPARAMS),
					*NE	=rows[1]+(kc+1*4)*(2+LWAV_NPARAMS),
					*NEE	=rows[1]+(kc+2*4)*(2+LWAV_NPARAMS),
					*NEEE	=rows[1]+(kc+3*4)*(2+LWAV_NPARAMS),
					*NEEEE	=rows[1]+(kc+4*4)*(2+LWAV_NPARAMS),
					*NEEEEE	=rows[1]+(kc+5*4)*(2+LWAV_NPARAMS),
					*WWWWW	=rows[0]+(kc-5*4)*(2+LWAV_NPARAMS),
					*WWWW	=rows[0]+(kc-4*4)*(2+LWAV_NPARAMS),
					*WWW	=rows[0]+(kc-3*4)*(2+LWAV_NPARAMS),
					*WW	=rows[0]+(kc-2*4)*(2+LWAV_NPARAMS),
					*W	=rows[0]+(kc-1*4)*(2+LWAV_NPARAMS),
					*curr	=rows[0]+(kc+0*4)*(2+LWAV_NPARAMS);
				(void)NNNN;
				(void)NNN;
				(void)NNWW;
				(void)NNW;
				(void)NN;
				(void)NNE;
				(void)NNEE;
				(void)NWW;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)NEEEE;
				(void)NEEEEE;
				(void)WWW;
				(void)WW;
				(void)W;

				//if(ky==10&&kx==10)//
				//	printf("");
				//int errors[]=
				//{
				//	abs(eN)+abs(eW)+abs(eNW),
				//	abs(eW)+abs(eNE)+abs(eN),
				//	abs(eN),
				//	abs(eW),
				//	abs(eN)*2+abs(eNN),
				//	abs(eW)*2+abs(eWW),
				//	abs(eN)+abs(eW),
				//	abs(eW)+abs(eNEE),
				//	//abs(eN)+abs(eNN)+abs(eNNN),
				//	//abs(eW)+abs(eWW)+abs(eWWW),
				//};
				int ctx[]=
				{
					0x2100,	N[0]+W[0]-NW[0],
				//	0x0DB8, (NNW[0]-NN[0]-2*N[0]+2*NE[0]-2*WW[0]+5*W[0])/3,
					0x0DB8, W[0]+NE[0]-N[0],
					0x0E22, N[0]-((N[1]+W[1]+NE[1])*-0x005C>>8),
					0x181F, W[0]-((N[1]+W[1]+NW[1])*-0x005B>>8),
					0x0BF3, N[0]-((0x00DF*NW[1]+0x0051*N[1]+0x00BD*NE[1]+0x005C*(NN[0]-N[0])-0x0102*(NW[0]-W[0]))>>8),

				//	13,	N[0]+W[0]-NW[0],
				//	1,	N[0]+W[0]-((NN[0]+WW[0])>>1),
				//	1,	2*N[0]-NN[0],
				//	1,	2*W[0]-WW[0],

				//	(W+NE)>>1,
				//	W,
				//	N,
				//	W+NE-N,
				//	W+((5*(N-NW)+NE-WW)>>3),//AV5
				//	W+((N+9*(N-NW)+2*(NE-NN-WW)+NNW-NNE-NWW)>>4),//AV9

					//3*(N-NN)+NNN,//X
					//3*(W-WW)+WWW,
					//(N+W+NE+NW)>>2,//X

					//rows[1][kc-1*4],//NW
					//rows[1][kc+0*4],//N
					//rows[1][kc+1*4],//NE
					//rows[0][kc-1*4],//W
				};
				MEDIAN3_32(ctx[0<<1|1], N[0], W[0], ctx[0<<1|1]);
			//	ctx[6]=ctx[5]=ctx[4]=ctx[3]=ctx[2]=ctx[1]=ctx[0];
			//	MEDIAN3_32(ctx[1], W, NE, ctx[1]);//X
				//if(NW&&N&&NE&&W)//
				//	printf("");
				long long lpred=0;
				int wsum=0;
				for(int k=0;k<LWAV_NPARAMS;++k)
				{
					int weight=params[kc][k]*ctx[k<<1|0]/(NW[k+2]+N[k+2]+NE[k+2]+WW[k+2]+W[k+2]+1);
					lpred+=(long long)weight*ctx[k<<1|1];
					wsum+=weight;
				}
				if(wsum)
					pred=(int)((lpred+(wsum>>1))/wsum);
				else
					pred=ctx[0<<1|1];
				//pred=(int)(lpred>>(LWAV_LPARAMS+LWAV_NBITS));
				//pred=(int)((
				//	+(long long)params[kc][0]*ctx[0]
				//	+(long long)params[kc][1]*ctx[1]
				//	+(long long)params[kc][2]*ctx[2]
				//	+(long long)params[kc][3]*ctx[3]
				//)>>(LWAV_LPARAMS+LWAV_NBITS));

				CLAMP3_32(pred, pred, N[0], W[0], NE[0]);

				{
					int p2=(pred+(1<<LWAV_PBITS>>1))>>LWAV_PBITS;
					int c=src->data[idx+kc];
					p2^=fwdmask;
					p2-=fwdmask;
					p2+=c;

					p2<<=32-src->depth[kc];
					p2>>=32-src->depth[kc];

					src->data[idx+kc]=p2;
					curr[0]=(fwd?c:p2)<<LWAV_PBITS;
					curr[1]=curr[0]-pred;
				}

				int error=(2*W[1]+curr[1]+N[1]+NE[1]+NEEE[1]);
				//int error=
				//	+rows[0][kc+4-1*4]*2
				//	+rows[0][kc+4+0*4]
				//	+rows[1][kc+4+0*4]
				//	+rows[1][kc+4+1*4]
				//	+rows[1][kc+4+3*4];
				//int error=
				//	+rows[0][kc+4-1*4]
				//	+rows[0][kc+4+0*4]*2
				//	+rows[1][kc+4+0*4]
				//	+rows[1][kc+4+1*4];
				//int error=rows[0][kc+4];//X
				int vmax=0;
				for(int k=0;k<LWAV_NPARAMS;++k)
				{
					curr[k+2]=abs(curr[0]-ctx[k<<1|1]);
					params[kc][k]+=(int)((long long)error*ctx[k<<1|1]>>src->depth[kc]);
					//params[kc][k]+=(int)((long long)error*errors[k]*ctx[k]>>src->depth[kc]>>src->depth[kc]);
					//if(abs(params[kc][k])>1<<LWAV_NBITS)//
					//	printf("");
					CLAMP2_32(params[kc][k], params[kc][k], -(1<<LWAV_NBITS), 1<<LWAV_NBITS);
					//int val=abs(params[kc][k]);
					//if(vmax<val)
					//	vmax=val;
				}
				//if(vmax>=0x100000)
				//{
				//	for(int k=0;k<LWAV_NPARAMS;++k)
				//		params[kc][k]>>=1;
				//}
			//	params[kc][0]+=error*NW	<<LWAV_NBITS>>src->depth[kc]*2;	CLAMP2_32(params[kc][0], params[kc][0], -(1<<LWAV_NBITS), 1<<LWAV_NBITS);
			//	params[kc][1]+=error*N	<<LWAV_NBITS>>src->depth[kc]*2;	CLAMP2_32(params[kc][1], params[kc][1], -(1<<LWAV_NBITS), 1<<LWAV_NBITS);
			//	params[kc][2]+=error*NE	<<LWAV_NBITS>>src->depth[kc]*2;	CLAMP2_32(params[kc][2], params[kc][2], -(1<<LWAV_NBITS), 1<<LWAV_NBITS);
			//	params[kc][3]+=error*W	<<LWAV_NBITS>>src->depth[kc]*2;	CLAMP2_32(params[kc][3], params[kc][3], -(1<<LWAV_NBITS), 1<<LWAV_NBITS);

			//	rows[0][kc+4+0*4]=error;
			}
			rows[0]+=4*(2+LWAV_NPARAMS);
			rows[1]+=4*(2+LWAV_NPARAMS);
			rows[2]+=4*(2+LWAV_NPARAMS);
			rows[3]+=4*(2+LWAV_NPARAMS);
		}
	}
	free(pixels);
	//set_window_title("sh=%d", sh);
}


//	#define WPU_TRACE
	#define WPU_UPDATE_LUMA

#define WPU_NPREDS 8
#define WPU_PREDLIST\
	WPU_PRED((N+W+NE-NW)>>1)\
	WPU_PRED(W+(eW>>1))\
	WPU_PRED(N+(eN>>1))\
	WPU_PRED(W+NW-NWW)\
	WPU_PRED(W+NE-N)\
	WPU_PRED((3*N+W-(NNW+NNE))>>1)\
	WPU_PRED(N+W-NW+((eN+eW-eNW)>>1))\
	WPU_PRED(N+NE-NNE)
//	WPU_PRED(grad)
//	WPU_PRED((W+NEE)>>1)
//	WPU_PRED((3*W+NEEE)>>2)
//	WPU_PRED(((abs(eW)+1)*N+(abs(eN)+1)*W)/(abs(eN)+abs(eW)+2))
//	WPU_PRED(((abs(eNW)+1)*NE+(abs(eNE)+1)*NW)/(abs(eNE)+abs(eNW)+2))
//	WPU_PRED(median1)
//	WPU_PRED(median2)
//	WPU_PRED(cgrad)
#define WPU_NBITS 8
#define WPU_PIXEL_STRIDE (4*(WPU_NPREDS+2LL))
void pred_WPU(Image *src, int fwd)
{
	int nch;
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	int halfs[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
		nlevels[3]>>1,
	};
	int boosthalfs[]=
	{
		nlevels[0]<<WPU_NBITS>>1,
		nlevels[1]<<WPU_NBITS>>1,
		nlevels[2]<<WPU_NBITS>>1,
		nlevels[3]<<WPU_NBITS>>1,
	};
	int perm[]={1, 2, 0, 3};
#ifdef WPU_TRACE
	long long trace[WPU_NPREDS]={0};
#endif
	int eprev[WPU_NPREDS]={0};
	int weights[WPU_NPREDS<<2]={0};
	size_t bufsize;
	int *pixels;

	for(int k=0;k<(WPU_NPREDS<<2);++k)
		weights[k]=0x8000;
	bufsize=(src->iw+8LL)*sizeof(int[4*WPU_PIXEL_STRIDE]);
	pixels=(int*)malloc(bufsize);//4 padded rows * 4 channels max * {pixels, errors, pred_errors...}
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);//{...padding    Y U V A  eY eU eV eA  a0Y a0U a0V a0A  a1Y a1U a1V a1A ...} * 4 rows
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+((src->iw+8LL)*((ky-0LL)&3)+4)*WPU_PIXEL_STRIDE,
			pixels+((src->iw+8LL)*((ky-1LL)&3)+4)*WPU_PIXEL_STRIDE,
			pixels+((src->iw+8LL)*((ky-2LL)&3)+4)*WPU_PIXEL_STRIDE,
			pixels+((src->iw+8LL)*((ky-3LL)&3)+4)*WPU_PIXEL_STRIDE,
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
#ifdef WPU_UPDATE_LUMA
			if(!fwd)
			{
				int
					cb	=src->data[idx+2],
					cr	=src->data[idx+0];
				int update=(cb+cr)>>3;
				update=CLAMP(-halfs[1], update, halfs[1]-1);
				{
					int luma=src->data[idx+1]-update;//subtract update
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					src->data[idx+1]=luma;
				}
			}
#endif
			for(int kc0=0;kc0<src->nch;++kc0)
			{
				int kc=perm[kc0], pred, vmin, vmax, curr, val;
				int
					NNN	=rows[3][kc+0*WPU_PIXEL_STRIDE+0],
					NNWW	=rows[2][kc-2*WPU_PIXEL_STRIDE+0],
					NNW	=rows[2][kc-1*WPU_PIXEL_STRIDE+0],
					NN	=rows[2][kc+0*WPU_PIXEL_STRIDE+0],
					NNE	=rows[2][kc+1*WPU_PIXEL_STRIDE+0],
					NNEE	=rows[2][kc+2*WPU_PIXEL_STRIDE+0],
					NNEEE	=rows[2][kc+3*WPU_PIXEL_STRIDE+0],
					NWW	=rows[1][kc-2*WPU_PIXEL_STRIDE+0],
					NW	=rows[1][kc-1*WPU_PIXEL_STRIDE+0],
					N	=rows[1][kc+0*WPU_PIXEL_STRIDE+0],
					NE	=rows[1][kc+1*WPU_PIXEL_STRIDE+0],
					NEE	=rows[1][kc+2*WPU_PIXEL_STRIDE+0],
					NEEE	=rows[1][kc+3*WPU_PIXEL_STRIDE+0],
					WWW	=rows[0][kc-3*WPU_PIXEL_STRIDE+0],
					WW	=rows[0][kc-2*WPU_PIXEL_STRIDE+0],
					W	=rows[0][kc-1*WPU_PIXEL_STRIDE+0],

					eNN	=rows[2][kc+0*WPU_PIXEL_STRIDE+4],
					eNW	=rows[1][kc-1*WPU_PIXEL_STRIDE+4],
					eN	=rows[1][kc+0*WPU_PIXEL_STRIDE+4],
					eNE	=rows[1][kc+1*WPU_PIXEL_STRIDE+4],
					eW	=rows[0][kc-1*WPU_PIXEL_STRIDE+4],
					offset	=0;
				(void)NNN;
				(void)NNWW;
				(void)NNW;
				(void)NN;
				(void)NNE;
				(void)NNEE;
				(void)NNEEE;
				(void)NWW;
				(void)NW;
				(void)N;
				(void)NE;
				(void)NEE;
				(void)NEEE;
				(void)WWW;
				(void)WW;
				(void)W;
				(void)eNN;
				(void)eNW;
				(void)eN;
				(void)eNE;
				(void)eW;
				if(kc0>0)
				{
					offset+=rows[0][1];
					if(kc0>1)
						offset=(2*offset+rows[0][2])>>1;//(2*g+[b-g])/2 = (g+b)/2
				}

				//int
				//	gx=abs(W-WW)+abs(N-NW)+abs(NE-N)+1,
				//	gy=abs(W-NW)+abs(N-NN)+abs(NE-NNE)+1,
				//	g135=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW)+abs(NE-NN)+1,
				//	g45=abs(W-N)+abs(N-NNE)+abs(NE-NNEE)+abs(WW-NW)+1;
				//int grad=(gx*N+gy*W)/(gx+gy);
				
				vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
				UPDATE_MIN(vmin, NE);
				UPDATE_MAX(vmax, NE);
#if 0
				int vmin=MINVAR(N, W), vmax=MAXVAR(N, W), pred, curr, val;
				int cgrad=N+W-NW;
				cgrad=CLAMP(vmin, cgrad, vmax);
				UPDATE_MIN(vmin, NE);
				UPDATE_MAX(vmax, NE);
				int nb[]=
				{
					N,
					W,
					NE,
					NW,
				};
				int temp;
				if(nb[0]>nb[1])SWAPVAR(nb[0], nb[1], temp);
				if(nb[2]>nb[3])SWAPVAR(nb[2], nb[3], temp);
				if(nb[0]>nb[2])SWAPVAR(nb[0], nb[2], temp);
				if(nb[1]>nb[3])SWAPVAR(nb[1], nb[3], temp);
				if(nb[1]>nb[2])SWAPVAR(nb[1], nb[2], temp);
				int median1=(nb[1]+nb[2])>>1;
				int median2=(nb[0]+nb[3])>>1;
#endif
				{
					int preds[]=
					{
#define WPU_PRED(X) X,
						WPU_PREDLIST
#undef  WPU_PRED
						//(N+W+NE-NW)>>1,
						//(3*N+W-(NNW+NNE))>>1,
						//N+W-NW,
						//W+(eW>>1),
						//N+(eN>>1),
						//W+NW-NWW,
						//W+NE-N,
						//N+NE-NNE,

						//(3*(N+W-NW)+((N+W)>>1))>>2,
						//(4*(N+W)+NE-NW)>>3,
						//NE+NEE-NEEE,

						//W+NE-N,
						//W+((eN+eW+eNW)>>2),
						//N+((eN+eW+eNE)>>2),
						//N+((4*eNW+0*eN+3*eNE-23*(NN-N)+2*(NW-W))>>5),

						//lib/jxl/modular/encoding/context_predict.h:	default lossless8
						//prediction[0] = W + NE - N;
						//prediction[1] = N - (((sumWN + teNE) * 8) >> 5);
						//prediction[2] = W - (((sumWN + teNW) * 8) >> 5);
						//prediction[3] = N - ((teNW * 4 + teN * 0 + teNE * 3 + (NN - N) * 23 + (NW - W) * 2) >> 5);
						//W+NE-N,
						//W-((eN+eW+eNW)>>2),
						//N-((eN+eW+eNE)>>2),
						//N-((4*eNW+0*eN+3*eNE+23*(N-NN)+2*(NW-W))>>5),

						//W+NE-N,
						//W-(int)(-13*((long long)eN+eW+eNW)>>8),
						//N+(int)(-7*((long long)eN+eW+eNE)>>8),
						//N-(int)((-0x85LL*eNW-0x63LL*eN-0xA2LL*eNE+0x17*((long long)N-NN)+0xF2*((long long)NW-W))>>8),
					};
					//if(ky==1&&kx==230)//
					//if(ky==0&&kx==141)//
					//if(ky==100&&kx==100)//
					//	printf("");

					long long llpred=0, wsum=1;
					for(int k=0;k<WPU_NPREDS;++k)
					{
						long long weight=0;
						weight+=(long long)rows[1][kc+0*WPU_PIXEL_STRIDE+8+((size_t)k<<2)];
#ifdef WPU_TRACE
						if(kc==1)
							trace[k]+=weight*weight;
#endif
						weight+=(long long)rows[1][kc-1*WPU_PIXEL_STRIDE+8+((size_t)k<<2)];
						weight+=(long long)rows[1][kc+1*WPU_PIXEL_STRIDE+8+((size_t)k<<2)];
						weight+=!weight;
						weight=((long long)weights[k<<2|kc]<<16)/weight;
						llpred+=weight*preds[k];
						wsum+=weight;
					}
					llpred+=(wsum>>1)-1;
					wsum+=!wsum;
					llpred/=wsum;
					pred=(int)llpred;
					pred=CLAMP(vmin, pred, vmax);
					pred=((64-4)*pred+N+W+NW+NE)>>6;

					pred+=offset;
					pred=CLAMP(-boosthalfs[kc], pred, boosthalfs[kc]-1);

					val=src->data[idx+kc];
					if(fwd)
					{
						curr=val;
						val-=(pred+(1<<WPU_NBITS>>1))>>WPU_NBITS;
						val+=halfs[kc];
						val&=nlevels[kc]-1;
						val-=halfs[kc];
					}
					else
					{
						val+=(pred+(1<<WPU_NBITS>>1))>>WPU_NBITS;
						val+=halfs[kc];
						val&=nlevels[kc]-1;
						val-=halfs[kc];
						curr=val;
					}
					src->data[idx+kc]=val;
					curr<<=WPU_NBITS;
					rows[0][kc+0]=curr-offset;
					rows[0][kc+4]=curr-pred;
					{
						int kbest=0;
						for(int k=0;k<WPU_NPREDS;++k)
						{
							int e2=abs(curr-preds[k]);
							eprev[k]=e2;
							rows[0][kc+0*WPU_PIXEL_STRIDE+8+((size_t)k<<2)]=e2;
							rows[1][kc+1*WPU_PIXEL_STRIDE+8+((size_t)k<<2)]+=e2;
							if(eprev[kbest]>eprev[k])
								kbest=k;
						}
						(void)kbest;
						//++weights[kbest<<2|kc];
						//if(weights[kbest<<2|kc]>0x8000)
						//{
						//	for(int k=0;k<WPU_NPREDS;++k)
						//		weights[k<<2|kc]>>=1;
						//}
						{
							int rescale=0;
							for(int k=0;k<WPU_NPREDS;++k)
							{
								weights[k<<2|kc]-=eprev[k];
								rescale|=weights[k<<2|kc]<0x800;
							}
							if(rescale)
							{
								for(int k=0;k<WPU_NPREDS;++k)
									weights[k<<2|kc]<<=1;
							}
						}
					}
				}
			}
#ifdef WPU_UPDATE_LUMA
			if(fwd)
			{
				int
					cb	=src->data[idx+2],
					cr	=src->data[idx+0];
				int update=(cb+cr)>>3;
				update=CLAMP(-halfs[1], update, halfs[1]-1);
				{
					int luma=src->data[idx+1]+update;//add update
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					src->data[idx+1]=luma;
				}
			}
#endif
			rows[0]+=WPU_PIXEL_STRIDE;
			rows[1]+=WPU_PIXEL_STRIDE;
			rows[2]+=WPU_PIXEL_STRIDE;
			rows[3]+=WPU_PIXEL_STRIDE;
		}
	}
#ifdef WPU_TRACE
	if(loud_transforms)
	{
		const char *prednames[]=
		{
#define WPU_PRED(X) #X,
			WPU_PREDLIST
#undef  WPU_PRED
		};
		char buf[2048]={0};
		int printed=0;
		for(int k=0;k<WPU_NPREDS;++k)
			printed+=snprintf(buf+printed, 1024LL-printed-1, "%20lld  %s\n", trace[k], prednames[k]);
			//printed+=snprintf(buf+printed, 1024LL-printed-1, "0x%016llX  %s\n", trace[k], prednames[k]);
		copy_to_clipboard(buf, printed);
		messagebox(MBOX_OK, "Copied to clipboard", buf);
	}
#endif
	free(pixels);
}

//(N+W)>>1
void pred_av2(Image *src, int fwd)
{
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
		1<<src->depth[3],
	};
	__m128i mhalf=_mm_set_epi32(
		nlevels[3]>>1,
		nlevels[2]>>1,
		nlevels[1]>>1,
		nlevels[0]>>1
	);
	__m128i symmask=_mm_set_epi32(
		nlevels[3]-1,
		nlevels[2]-1,
		nlevels[1]-1,
		nlevels[0]-1
	);
	ALIGN(16) int curr[4]={0};

	int *pixels=(int*)_mm_malloc((src->iw+4LL)*sizeof(int[2*4]), sizeof(__m128i));//2 padded rows * 4 channels max
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, (src->iw+4LL)*sizeof(int[2*4]));
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *rows[]=
		{
			pixels+(((src->iw+4LL)*((ky-0LL)&1)+1)<<2),
			pixels+(((src->iw+4LL)*((ky-1LL)&1)+1)<<2),
		};
		int kx=0;
		if(fwd)
		{
			for(;kx<src->iw;++kx, idx+=4)
			{
				__m128i N	=_mm_load_si128((__m128i*)rows[1]+0);
				__m128i W	=_mm_load_si128((__m128i*)rows[0]-1);
				__m128i pred=_mm_add_epi32(N, W);
				pred=_mm_srai_epi32(pred, 1);
				{
					__m128i mc=_mm_loadu_si128((__m128i*)(src->data+idx));
					pred=_mm_sub_epi32(mc, pred);
					pred=_mm_add_epi32(pred, mhalf);
					pred=_mm_and_si128(pred, symmask);
					pred=_mm_sub_epi32(pred, mhalf);
					_mm_store_si128((__m128i*)curr, pred);//error
					_mm_store_si128((__m128i*)rows[0], mc);//pixel
				}
				src->data[idx|0]=curr[0];
				src->data[idx|1]=curr[1];
				src->data[idx|2]=curr[2];

				rows[0]+=4;
				rows[1]+=4;
			}
		}
		else
		{
			for(;kx<src->iw;++kx, idx+=4)
			{
				__m128i N	=_mm_load_si128((__m128i*)rows[1]+0);
				__m128i W	=_mm_load_si128((__m128i*)rows[0]-1);
				__m128i pred=_mm_add_epi32(N, W);
				pred=_mm_srai_epi32(pred, 1);
				{
					__m128i mc=_mm_loadu_si128((__m128i*)(src->data+idx));
					pred=_mm_add_epi32(mc, pred);
					pred=_mm_add_epi32(pred, mhalf);
					pred=_mm_and_si128(pred, symmask);
					pred=_mm_sub_epi32(pred, mhalf);
					_mm_store_si128((__m128i*)curr, pred);//pixel
					_mm_store_si128((__m128i*)rows[0], pred);//pixel
				}
				src->data[idx|0]=curr[0];
				src->data[idx|1]=curr[1];
				src->data[idx|2]=curr[2];

				rows[0]+=4;
				rows[1]+=4;
			}
		}
	}
	_mm_free(pixels);
}


	#define WP_RCT
	#define ETOTAL_MBOX

#define NWP 8LL		//multiple of 4
#define LGBLOCKSIZE 0
#define BLOCKSIZE (1<<LGBLOCKSIZE)
#define WP_PBITS 5
#define WP_WBITS 10
void pred_wp_deferred(Image *src, int fwd)
{
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
	};
	__m128i pixelmask=_mm_set_epi32(0, nlevels[2]-1, nlevels[1]-1, nlevels[0]-1);
	__m128i pixelhalf=_mm_set_epi32(0, nlevels[2]>>1, nlevels[1]>>1, nlevels[0]>>1);
	__m128i roundoffset=_mm_set1_epi32(1<<WP_PBITS>>1);
	__m128i fwdmask=_mm_set1_epi32(-fwd);
	__m128i invmask=_mm_set1_epi32(-!fwd);
	ALIGN(16) int pred_errors[4*NWP]={0};//interleaved channels
#ifdef ETOTAL_MBOX
	long long etotal[4*NWP]={0};
#endif
	//if(src->nch!=3)
	//{
	//	LOG_ERROR("Expected 3 channels, got %d", src->nch);
	//	return;
	//}
	int nblocks=(src->iw+BLOCKSIZE-1)>>LGBLOCKSIZE;
	int *weights=(int*)_mm_malloc((nblocks+NWP)*sizeof(int[4*NWP]), sizeof(__m256i));//interleaved channels
	int *pixels=(int*)malloc((src->iw+16LL)*sizeof(int[4*4*2]));//4 padded rows * 4 channels max * {pixel, error}		int, to accelerate _mullo_epi32
	if(!weights||!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	weights[0]=(int)((1LL<<WP_WBITS)/NWP);
	memfill(weights+1, weights, (nblocks+NWP)*sizeof(int[4*NWP])-sizeof(*weights), sizeof(*weights));
	//memset(weights, 0, (nblocks+NWP)*sizeof(int[4*NWP]));
	memset(pixels, 0, (src->iw+16LL)*sizeof(int[4*4*2]));
#ifdef WP_RCT
	src->depth[0]+=(char)fwd;
	src->depth[2]+=(char)fwd;
#endif
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int *wp=weights;
		int *rows[]=
		{
			pixels+(((src->iw+16LL)*((ky-0LL)&3)+4LL)<<3),
			pixels+(((src->iw+16LL)*((ky-1LL)&3)+4LL)<<3),
			pixels+(((src->iw+16LL)*((ky-2LL)&3)+4LL)<<3),
			pixels+(((src->iw+16LL)*((ky-3LL)&3)+4LL)<<3),
		};
		for(int kx=0;kx<src->iw;++kx, idx+=4)
		{
			__m128i NN	=_mm_loadu_si128((__m128i*)(rows[2]+0*8+0));
			__m128i NNE	=_mm_loadu_si128((__m128i*)(rows[2]+1*8+0));
			__m128i NW	=_mm_loadu_si128((__m128i*)(rows[1]-1*8+0));
			__m128i N	=_mm_loadu_si128((__m128i*)(rows[1]+0*8+0));
			__m128i NE	=_mm_loadu_si128((__m128i*)(rows[1]+1*8+0));
			__m128i NEE	=_mm_loadu_si128((__m128i*)(rows[1]+2*8+0));
			__m128i NEEE	=_mm_loadu_si128((__m128i*)(rows[1]+3*8+0));
			__m128i W	=_mm_loadu_si128((__m128i*)(rows[0]-1*8+0));
			__m128i eNN	=_mm_loadu_si128((__m128i*)(rows[2]+0*8+4));
			__m128i eNW	=_mm_loadu_si128((__m128i*)(rows[1]-1*8+4));
			__m128i eN	=_mm_loadu_si128((__m128i*)(rows[1]+0*8+4));
			__m128i eNE	=_mm_loadu_si128((__m128i*)(rows[1]+1*8+4));
			__m128i eW	=_mm_loadu_si128((__m128i*)(rows[0]-1*8+4));
			
			__m128i vmin=_mm_min_epi32(N, W);
			__m128i vmax=_mm_max_epi32(N, W);
			vmin=_mm_min_epi32(vmin, NE);
			vmax=_mm_max_epi32(vmax, NE);
			//wp0 = W+NE-N
			//wp1 = W		//W+((eN+eW+eNW)*5>>4)
			//wp2 = N		//N+((eN+eW+eNE)*5>>4)
			//wp3 = (N+W)>>1	//N+NW-W+((-eNN*13-eN*5-eNE*11+(N-NN)*5)>>4)
			//wp4 = N+NE-NNE	//(N+W)>>1
			//wp5 = N+W-NW
			//wp6 = (W+NEE)>>1
			//wp7 = (3*W+NEEE)>>2
#if 0
			ALIGN(16) int spred[16]={0};
			for(int kc=0;kc<4;++kc)
			{
				spred[kc|0*4]=W.m128i_i32[kc]+NE.m128i_i32[kc]-N.m128i_i32[kc];
				spred[kc|1*4]=W.m128i_i32[kc]+((eN.m128i_i32[kc]+eW.m128i_i32[kc]+eNW.m128i_i32[kc])*5>>4);
				spred[kc|2*4]=N.m128i_i32[kc]+((eN.m128i_i32[kc]+eW.m128i_i32[kc]+eNE.m128i_i32[kc])*5>>4);
				spred[kc|3*4]=N.m128i_i32[kc]+NW.m128i_i32[kc]-W.m128i_i32[kc]
					+((-eNN.m128i_i32[kc]*13-eN.m128i_i32[kc]*5-eNE.m128i_i32[kc]*11+(N.m128i_i32[kc]-NN.m128i_i32[kc])*5)>>4);
			}
#endif
			{
				__m128i e=_mm_add_epi32(eN, eW);
				__m128i wp0=_mm_sub_epi32(_mm_add_epi32(W, NE), N);
				__m128i wp1=_mm_add_epi32(e, eNW);
				__m128i wp2=_mm_add_epi32(e, eNE);
				__m128i wp3=_mm_sub_epi32(_mm_add_epi32(N, NW), W);
				__m128i wp4=_mm_sub_epi32(_mm_add_epi32(N, NE), NNE);
				__m128i wp5=_mm_sub_epi32(_mm_add_epi32(N, W), NW);
				__m128i wp6=_mm_add_epi32(W, NEE);
				__m128i wp7=_mm_add_epi32(W, _mm_slli_epi32(W, 1));
				wp6=_mm_srai_epi32(wp6, 1);
				wp7=_mm_srai_epi32(_mm_add_epi32(wp7, NEEE), 2);
				{
					__m128i e2=_mm_add_epi32(eNE, _mm_add_epi32(_mm_slli_epi32(eNE, 1), _mm_slli_epi32(eNE, 3)));
					e=_mm_add_epi32(eNN, _mm_add_epi32(_mm_slli_epi32(eNN, 2), _mm_slli_epi32(eNN, 3)));
					e=_mm_add_epi32(e, _mm_add_epi32(eN, _mm_slli_epi32(eN, 2)));
					e=_mm_add_epi32(e, e2);
					e2=_mm_sub_epi32(N, NN);
					e2=_mm_add_epi32(e2, _mm_slli_epi32(e2, 2));
					e2=_mm_sub_epi32(e2, e);
					e2=_mm_srai_epi32(e2, 4);
				}

				//wp1=_mm_add_epi32(W, _mm_srai_epi32(_mm_add_epi32(wp1, _mm_slli_epi32(wp1, 2)), 4));
			//wp2=_mm_add_epi32(N, _mm_srai_epi32(_mm_add_epi32(wp2, _mm_slli_epi32(wp2, 2)), 4));
			//wp3=_mm_add_epi32(wp3, e2);
				wp1=W;
				//wp1=_mm_srai_epi32(_mm_add_epi32(_mm_slli_epi32(W, 2), _mm_add_epi32(_mm_add_epi32(N, NE), _mm_add_epi32(NEE, NEEE))), 3);
				wp2=N;
				wp3=_mm_srai_epi32(_mm_add_epi32(N, W), 1);
				//wp3=_mm_slli_epi32(_mm_add_epi32(N, W), 2);
			//wp3=_mm_add_epi32(wp3, _mm_sub_epi32(NE, NW));
			//wp3=_mm_srai_epi32(wp3, 3);
#if 0
				//for(int kc=0;kc<4;++kc)
				//{
				//	if(wp1.m128i_i32[kc]!=spred[kc|0*4])
				//		LOG_ERROR("");
				//	if(wp2.m128i_i32[kc]!=spred[kc|1*4])
				//		LOG_ERROR("");
				//	if(wp3.m128i_i32[kc]!=spred[kc|2*4])
				//		LOG_ERROR("");
				//	if(wp4.m128i_i32[kc]!=spred[kc|3*4])
				//		LOG_ERROR("");
				//}
				wp1=_mm_load_si128((__m128i*)spred+0);//
				wp2=_mm_load_si128((__m128i*)spred+1);//
				wp3=_mm_load_si128((__m128i*)spred+2);//
				wp4=_mm_load_si128((__m128i*)spred+3);//
#endif
				{
					ALIGN(16) int curr[]=
					{
						src->data[idx+0],
						src->data[idx+1],
						src->data[idx+2],
						//src->data[idx+1],//g		X  can't permute VYU -> YUV
						//src->data[idx+2],//b
						//src->data[idx+0],//r
						0,
					};

					__m128i w0=_mm_load_si128((__m128i*)wp+0);
					__m128i w1=_mm_load_si128((__m128i*)wp+1);
					__m128i w2=_mm_load_si128((__m128i*)wp+2);
					__m128i w3=_mm_load_si128((__m128i*)wp+3);
					__m128i w4=_mm_load_si128((__m128i*)wp+4);
					__m128i w5=_mm_load_si128((__m128i*)wp+5);
					__m128i w6=_mm_load_si128((__m128i*)wp+6);
					__m128i w7=_mm_load_si128((__m128i*)wp+7);
					w0=_mm_mullo_epi32(w0, wp0);//FIXME use AVX2
					w1=_mm_mullo_epi32(w1, wp1);
					w2=_mm_mullo_epi32(w2, wp2);
					w3=_mm_mullo_epi32(w3, wp3);
					w4=_mm_mullo_epi32(w4, wp4);
					w5=_mm_mullo_epi32(w5, wp5);
					w6=_mm_mullo_epi32(w6, wp6);
					w7=_mm_mullo_epi32(w7, wp7);
					w0=_mm_add_epi32(w0, w1);
					w0=_mm_add_epi32(w0, w2);
					w0=_mm_add_epi32(w0, w3);
					w0=_mm_add_epi32(w0, w4);
					w0=_mm_add_epi32(w0, w5);
					w0=_mm_add_epi32(w0, w6);
					w0=_mm_add_epi32(w0, w7);
					w0=_mm_srai_epi32(w0, WP_WBITS);
					w0=_mm_min_epi32(w0, vmax);
					w0=_mm_max_epi32(w0, vmin);
#ifdef WP_RCT
					if(fwd)
					{
						curr[0]-=curr[1];			//curr[0]=((curr[0]+(nlevels[0]>>1))&(nlevels[0]-1))-(nlevels[0]>>1);
						curr[2]-=curr[1];			//curr[1]=((curr[1]+(nlevels[1]>>1))&(nlevels[1]-1))-(nlevels[1]>>1);
						curr[1]+=(curr[0]+curr[2])>>2;		//curr[2]=((curr[2]+(nlevels[2]>>1))&(nlevels[2]-1))-(nlevels[2]>>1);

						//curr[2]-=curr[0];//Cr
						//curr[1]-=curr[0];//Cb
						//curr[0]+=(curr[2]+curr[1])>>2;//Y
					}
#endif
					{
						__m128i mc=_mm_load_si128((__m128i*)curr);
						__m128i pixel=_mm_and_si128(mc, fwdmask);
						__m128i crisppred=_mm_add_epi32(w0, roundoffset);
						crisppred=_mm_srai_epi32(crisppred, WP_PBITS);

						crisppred=_mm_xor_si128(crisppred, fwdmask);
						crisppred=_mm_sub_epi32(crisppred, fwdmask);
						crisppred=_mm_add_epi32(crisppred, mc);//error = curr - pred

						crisppred=_mm_add_epi32(crisppred, pixelhalf);
						crisppred=_mm_and_si128(crisppred, pixelmask);
						crisppred=_mm_sub_epi32(crisppred, pixelhalf);
						_mm_store_si128((__m128i*)curr, crisppred);
#ifdef WP_RCT
						if(!fwd)
						{
							curr[1]-=(curr[0]+curr[2])>>2;		//curr[2]=((curr[2]+(nlevels[2]>>1))&(nlevels[2]-1))-(nlevels[2]>>1);
							curr[2]+=curr[1];			//curr[1]=((curr[1]+(nlevels[1]>>1))&(nlevels[1]-1))-(nlevels[1]>>1);
							curr[0]+=curr[1];			//curr[0]=((curr[0]+(nlevels[0]>>1))&(nlevels[0]-1))-(nlevels[0]>>1);

							//curr[0]-=(curr[2]+curr[1])>>2;//g
							//curr[1]+=curr[0];//b
							//curr[2]+=curr[0];//r
						}
#endif
						src->data[idx+0]=curr[0];
						src->data[idx+1]=curr[1];
						src->data[idx+2]=curr[2];
						//src->data[idx+0]=curr[2];//r
						//src->data[idx+1]=curr[0];//g
						//src->data[idx+2]=curr[1];//b
						pixel=_mm_or_si128(pixel, _mm_and_si128(crisppred, invmask));
						pixel=_mm_slli_epi32(pixel, WP_PBITS);
						w1=_mm_sub_epi32(pixel, w1);
						_mm_storeu_si128((__m128i*)(rows[0]+8*0+4), w1);
						_mm_storeu_si128((__m128i*)(rows[0]+8*0+0), pixel);

						{
							__m128i pe0=_mm_load_si128((__m128i*)pred_errors+0);
							__m128i pe1=_mm_load_si128((__m128i*)pred_errors+1);
							__m128i pe2=_mm_load_si128((__m128i*)pred_errors+2);
							__m128i pe3=_mm_load_si128((__m128i*)pred_errors+3);
							__m128i pe4=_mm_load_si128((__m128i*)pred_errors+4);
							__m128i pe5=_mm_load_si128((__m128i*)pred_errors+5);
							__m128i pe6=_mm_load_si128((__m128i*)pred_errors+6);
							__m128i pe7=_mm_load_si128((__m128i*)pred_errors+7);
							wp0=_mm_sub_epi32(pixel, wp0);
							wp1=_mm_sub_epi32(pixel, wp1);
							wp2=_mm_sub_epi32(pixel, wp2);
							wp3=_mm_sub_epi32(pixel, wp3);
							wp4=_mm_sub_epi32(pixel, wp4);
							wp5=_mm_sub_epi32(pixel, wp5);
							wp6=_mm_sub_epi32(pixel, wp6);
							wp7=_mm_sub_epi32(pixel, wp7);
							wp0=_mm_abs_epi32(wp0);
							wp1=_mm_abs_epi32(wp1);
							wp2=_mm_abs_epi32(wp2);
							wp3=_mm_abs_epi32(wp3);
							wp4=_mm_abs_epi32(wp4);
							wp5=_mm_abs_epi32(wp5);
							wp6=_mm_abs_epi32(wp6);
							wp7=_mm_abs_epi32(wp7);
							pe0=_mm_add_epi32(pe0, wp0);
							pe1=_mm_add_epi32(pe1, wp1);
							pe2=_mm_add_epi32(pe2, wp2);
							pe3=_mm_add_epi32(pe3, wp3);
							pe4=_mm_add_epi32(pe4, wp4);
							pe5=_mm_add_epi32(pe5, wp5);
							pe6=_mm_add_epi32(pe6, wp6);
							pe7=_mm_add_epi32(pe7, wp7);
							_mm_store_si128((__m128i*)pred_errors+0, pe0);
							_mm_store_si128((__m128i*)pred_errors+1, pe1);
							_mm_store_si128((__m128i*)pred_errors+2, pe2);
							_mm_store_si128((__m128i*)pred_errors+3, pe3);
							_mm_store_si128((__m128i*)pred_errors+4, pe4);
							_mm_store_si128((__m128i*)pred_errors+5, pe5);
							_mm_store_si128((__m128i*)pred_errors+6, pe6);
							_mm_store_si128((__m128i*)pred_errors+7, pe7);
						}
					}
				}
			}
			
			if(kx&&!(kx&(BLOCKSIZE-1)))//update WP
			{
				for(int kc=0;kc<3;++kc)
				{
					long long
						weight0=0x100000000/((long long)pred_errors[kc+0*4]+1),
						weight1=0x100000000/((long long)pred_errors[kc+1*4]+1),
						weight2=0x100000000/((long long)pred_errors[kc+2*4]+1),
						weight3=0x100000000/((long long)pred_errors[kc+3*4]+1),
						weight4=0x100000000/((long long)pred_errors[kc+4*4]+1),
						weight5=0x100000000/((long long)pred_errors[kc+5*4]+1),
						weight6=0x100000000/((long long)pred_errors[kc+6*4]+1),
						weight7=0x100000000/((long long)pred_errors[kc+7*4]+1),
						sum=weight0+weight1+weight2+weight3+weight4+weight5+weight6+weight7+1;
					wp[kc+0*4]=(int)((weight0<<WP_WBITS)/sum);
					wp[kc+1*4]=(int)((weight1<<WP_WBITS)/sum);
					wp[kc+2*4]=(int)((weight2<<WP_WBITS)/sum);
					wp[kc+3*4]=(int)((weight3<<WP_WBITS)/sum);
					wp[kc+4*4]=(int)((weight4<<WP_WBITS)/sum);
					wp[kc+5*4]=(int)((weight5<<WP_WBITS)/sum);
					wp[kc+6*4]=(int)((weight6<<WP_WBITS)/sum);
					wp[kc+7*4]=(int)((weight7<<WP_WBITS)/sum);
#ifdef ETOTAL_MBOX
					etotal[kc+0*4]+=(long long)pred_errors[kc+0*4];
					etotal[kc+1*4]+=(long long)pred_errors[kc+1*4];
					etotal[kc+2*4]+=(long long)pred_errors[kc+2*4];
					etotal[kc+3*4]+=(long long)pred_errors[kc+3*4];
					etotal[kc+4*4]+=(long long)pred_errors[kc+4*4];
					etotal[kc+5*4]+=(long long)pred_errors[kc+5*4];
					etotal[kc+6*4]+=(long long)pred_errors[kc+6*4];
					etotal[kc+7*4]+=(long long)pred_errors[kc+7*4];
#endif
				}
				memset(pred_errors, 0, sizeof(pred_errors));
				wp+=4*NWP;
			}
			rows[0]+=8;
			rows[1]+=8;
			rows[2]+=8;
			rows[3]+=8;
		}
	}
#ifdef WP_RCT
	src->depth[0]-=!fwd;
	src->depth[2]-=!fwd;
#endif
#ifdef ETOTAL_MBOX
	{
#define BUFSIZE 1024
		char buf[1024]={0};
		int printed=0;
		long long esum[3]={0};
		for(int kp=0;kp<NWP;++kp)
		{
			esum[0]+=etotal[kp<<2|0];
			esum[1]+=etotal[kp<<2|1];
			esum[2]+=etotal[kp<<2|2];
		}
		for(int kp=0;kp<NWP;++kp)
		{
			for(int kc=0;kc<3;++kc)
				printed+=snprintf(buf+printed, BUFSIZE-printed-1, " %6.2lf", 100.*etotal[kp<<2|kc]/esum[kc]);
			printed+=snprintf(buf+printed, BUFSIZE-printed-1, "\n");
		}
		copy_to_clipboard(buf, printed);
		messagebox(MBOX_OK, "Copied to clipboard", "%s", buf);
#undef BUFSIZE
	}
#endif
	_mm_free(weights);
	free(pixels);
}

void pred_average(Image *src, int fwd, int enable_ma)
{
	Image *dst;
	const int *pixels;
//	const int *errors;

	dst=0;
	image_copy(&dst, src);
	if(!dst)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	pixels=fwd?src->data:dst->data;
//	errors=fwd?dst->data:src->data;
	for(int kc=0;kc<4;++kc)
	{
		int nlevels;

		if(!src->depth[kc])
			continue;
		nlevels=1<<src->depth[kc];
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
#define LOAD(BUF, X, Y) ((unsigned)(ky+(Y))<(unsigned)src->ih&&(unsigned)(kx+(X))<(unsigned)src->iw?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
				int
				//	NNE	=LOAD(pixels,  1, -2),
					NW	=LOAD(pixels, -1, -1),
					N	=LOAD(pixels,  0, -1),
					NE	=LOAD(pixels,  1, -1),
					W	=LOAD(pixels, -1,  0);
#undef  LOAD
				int pred=(4*(N+W)+NE-NW)>>3;
			//	int pred=(W+2*NE-NNE+1)>>1;
				
				pred^=-fwd;
				pred+=fwd;
				pred+=src->data[idx<<2|kc];
				if(enable_ma)
				{
					pred+=nlevels>>1;
					pred&=nlevels-1;
					pred-=nlevels>>1;
				}
				dst->data[idx<<2|kc]=pred;
				//dst->data[idx<<2|kc]=((src->data[idx<<2|kc]+pred+(nlevels>>1))&(nlevels-1))-(nlevels>>1);
			}
		}
	}
	memcpy(src->data, dst->data, (size_t)src->iw*src->ih*sizeof(int[4]));
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}
	free(dst);
}

void pred_ecoeff(Image *src, int fwd, int enable_ma)
{
	const int *pixels, *errors;
	Image *dst;

	dst=0;
	image_copy(&dst, src);
	if(!dst)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	pixels=fwd?dst->data:src->data;
	errors=fwd?src->data:dst->data;
	for(int kc=0;kc<4;++kc)
	{
		int nlevels, half;

		if(!src->depth[kc])
			continue;
		nlevels=1<<src->depth[kc];
		half=nlevels>>1;
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			int ecoeffs[12]={0};
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
#define LOAD(BUF, X, Y) ((unsigned)(ky+(Y))<(unsigned)src->ih&&(unsigned)(kx+(X))<(unsigned)src->iw?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
				int
				//	NNE	=LOAD(pixels,  1, -2),
					NW	=LOAD(pixels, -1, -1),
					N	=LOAD(pixels,  0, -1),
				//	NE	=LOAD(pixels,  1, -1),
					W	=LOAD(pixels, -1,  0),
					enb[]=
				{
					LOAD(errors, -2, -2),
					LOAD(errors, -1, -2),
					LOAD(errors,  0, -2),
					LOAD(errors,  1, -2),
					LOAD(errors,  2, -2),
					LOAD(errors, -2, -1),
					LOAD(errors, -1, -1),
					LOAD(errors,  0, -1),
					LOAD(errors,  1, -1),
					LOAD(errors,  2, -1),
					LOAD(errors, -2,  0),
					LOAD(errors, -1,  0),
				};
				int pred, correction, e;
#undef  LOAD
				correction=0;
				for(int k=0;k<(int)_countof(enb);++k)
					correction+=enb[k]*ecoeffs[k];
				//correction+=((1<<8)>>1)-1;
				//correction>>=8;
				correction+=12*256/2-1;
				correction/=12*256;

				//if(correction)
				//	printf("");

				pred=N+W-NW;
				pred=MEDIAN3(N, W, pred);
				//int pred=(4*(N+W)+NE-NW)>>3;
				//int pred=(W+2*NE-NNE+1)>>1;

				pred+=correction;
				pred=CLAMP(-half, pred, half-1);
				
				pred^=-fwd;
				pred+=fwd;
				pred+=dst->data[idx<<2|kc];
				if(enable_ma)
				{
					pred+=nlevels>>1;
					pred&=nlevels-1;
					pred-=nlevels>>1;
				}
				src->data[idx<<2|kc]=pred;

				//L = sq(curr - pred)
				//dL/dcoeff = (curr - pred) * -dpred/dcoeff
				//	= (pred - curr)*srcval
				//coeff -= dL/dcoeff * lr
				//coeff += (curr-pred)*srcval*lr
				e=errors[idx<<2|kc];
				for(int k=0;k<(int)_countof(enb);++k)
				{
					int c=ecoeffs[k];
					int update=e*enb[k];
					c+=(update>>((src->depth[kc]<<1)-12))+(update<0);
					ecoeffs[k]=CLAMP(-256, c, 256);
				}
			}
		}
	}
	//memcpy(src->data, dst->data, (size_t)src->iw*src->ih*sizeof(int[4]));
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}
	free(dst);
}

#if 0
static void zipper_predict(Image *src, int fwd)
{
	int negmask=-fwd;
	for(int ky=fwd?src->ih-1:0, step=fwd?-1:1;(unsigned)ky<(unsigned)src->ih;ky+=step)
	{
		for(int kx=fwd?src->iw-1:0;(unsigned)kx<(unsigned)src->iw;kx+=step)
		{
			int idx=src->iw*ky+kx;
			for(int kc=0;kc<src->nch;++kc)
			{
				int
					NW	=ky&&kx?src->data[(idx-src->iw-1)<<2|kc]:0,
					N	=ky?src->data[(idx-src->iw)<<2|kc]:0,
					W	=kx?src->data[(idx-1)<<2|kc]:0;
				int pred;
				MEDIAN3_32(pred, N, W, N+W-NW);
				
				pred^=negmask;
				pred-=negmask;
				pred+=src->data[idx<<2|kc];
				//pred<<=32-src->depth[kc];
				//pred>>=32-src->depth[kc];
				src->data[idx<<2|kc]=pred;
			}
		}
	}
}
static void zipper_update(Image *src, int fwd)
{
	int negmask=-fwd;
	for(int ky=fwd?src->ih-1:0, step=fwd?-1:1;(unsigned)ky<(unsigned)src->ih;ky+=step)
	{
		for(int kx=fwd?src->iw-1:0;(unsigned)kx<(unsigned)src->iw;kx+=step)
		{
			int idx=src->iw*ky+kx;
			for(int kc=0;kc<src->nch;++kc)
			{
				int pred;
				if(kc>1)
				{
					int
						NW	=ky<src->ih-1&&kx<src->iw-1?src->data[(idx+src->iw+1)<<2|kc]:0,
						N	=ky<src->ih-1?src->data[(idx+src->iw)<<2|kc]:0,
						W	=kx<src->iw-1?src->data[(idx+1)<<2|kc]:0;
					pred=(N+W-NW)>>4;
				}
				else
					pred=0;
				//if(kc)
				//	pred-=NW;
				//pred>>=4;
				//int pred=kc?(N+W-NW)>>4:(N+W)>>4;
				
				pred^=negmask;
				pred-=negmask;
				pred+=src->data[idx<<2|kc];
				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];
				src->data[idx<<2|kc]=pred;
			}
		}
	}
}
#endif
void pred_zipper(Image *src, int fwd, int enable_ma)
{
	int negmask=-fwd;
	int half[]=
	{
		1<<src->depth[0]>>1,
		1<<src->depth[1]>>1,
		1<<src->depth[2]>>1,
		1<<src->depth[3]>>1,
	};
	for(int ky=fwd?src->ih-1:0, step=fwd?-1:1;(unsigned)ky<(unsigned)src->ih;ky+=step)
	{
		for(int kx=fwd?src->iw-1:0;(unsigned)kx<(unsigned)src->iw;kx+=step)
		{
			int idx=src->iw*ky+kx;
			for(int kc=0;kc<src->nch;++kc)
			{
#define LOAD(X, Y) ((unsigned)(ky+(Y))<(unsigned)src->ih&&(unsigned)(kx+(X))<(unsigned)src->iw?src->data[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
				int
					NW	=LOAD(-1, -1),
					N	=LOAD( 0, -1),
					W	=LOAD(-1,  0),
					eE	=LOAD( 1,  0),
					eEE	=LOAD( 2,  0),
					eSWW	=LOAD(-2,  1),
					eSW	=LOAD(-1,  1),
					eS	=LOAD( 0,  1),
					eSE	=LOAD( 1,  1),
					eSEE	=LOAD( 2,  1),
					eSSWW	=LOAD(-2,  2),
					eSSW	=LOAD(-1,  2),
					eSS	=LOAD( 0,  2),
					eSSE	=LOAD( 1,  2),
					eSSEE	=LOAD( 2,  2);
				//int
				//	NW	=ky&&kx?src->data[(idx-src->iw-1)<<2|kc]:0,
				//	N	=ky?src->data[(idx-src->iw)<<2|kc]:0,
				//	W	=kx?src->data[(idx-1)<<2|kc]:0,
				//	eE	=kx<src->iw-1?src->data[(idx+1)<<2|kc]:0,
				//	eEE	=kx<src->iw-2?src->data[(idx+2)<<2|kc]:0,
				//	eSWW	=ky<src->ih-1&&kx>1?src->data[(idx+src->iw-2)<<2|kc]:0,
				//	eSW	=ky<src->ih-1&&kx>0?src->data[(idx+src->iw-1)<<2|kc]:0,
				//	eS	=ky<src->ih-1?src->data[(idx+src->iw)<<2|kc]:0,
				//	eSE	=ky<src->ih-1&&kx<src->iw-1?src->data[(idx+src->iw+1)<<2|kc]:0,
				//	eSEE	=ky<src->ih-1&&kx<src->iw-2?src->data[(idx+src->iw+2)<<2|kc]:0,
				//	eSS	=ky<src->ih-2?src->data[(idx+src->iw*2)<<2|kc]:0;
#undef  LOAD
				int pred;
				int
					vx=abs(eE-eEE)+abs(eS-eSE)+abs(eS-eSW),
					vy=abs(eSSW-eSW)+abs(eSS-eS)+abs(eSE-eE);
				int update=(vy>vx?eS:eE)/4;
				//int update=(abs(eS-eSS)>abs(eE-eEE)?eS:eE)/4;
				//int update=(abs(eS)>abs(eE)?eS:eE)/4;
				MEDIAN3_32(pred, N, W, N+W-NW);
				pred+=update;
				CLAMP2_32(pred, pred, -half[kc], half[kc]-1);
				
				pred^=negmask;
				pred-=negmask;
				pred+=src->data[idx<<2|kc];
				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];
				src->data[idx<<2|kc]=pred;
			}
		}
	}
	//if(fwd)
	//{
	//	zipper_predict(src, 1);
	//	for(int kc=2;kc<src->nch;++kc)
	//		src->depth[kc]+=src->depth[kc]!=0&&src->depth[kc]<24;
	//	zipper_update(src, 0);
	//}
	//else
	//{
	//	zipper_update(src, 1);
	//	for(int kc=2;kc<src->nch;++kc)
	//		src->depth[kc]-=src->depth[kc]>0;
	//	zipper_predict(src, 0);
	//}

#if 0
	Image *src=*psrc;
	int *row=(int*)malloc(src->iw*sizeof(int));
	//Image *dst=(Image*)malloc(sizeof(Image)+(src->iw+1LL)*(src->ih+1LL)*sizeof(int[4]));

	(void)fwd;//WIP
	(void)enable_ma;

	if(!row)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(int kc=0;kc<4;++kc)
	{
		int depth=src->depth[kc], nlevels=1<<depth;
		if(!depth)
			continue;
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			//int WW=0;
			int W=0;
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				int E=kx<src->iw-1?src->data[(idx+1)<<2|kc]:0;
				//int EE=kx<src->iw-2?src->data[(idx+2)<<2|kc]:0;
				int pred=(W+E)>>1;
				//int pred=W+E-((WW+EE)>>1);

				int curr=src->data[idx<<2|kc];
				curr-=pred;
				curr+=nlevels>>1;
				curr&=nlevels-1;
				curr-=nlevels>>1;
				src->data[idx<<2|kc]=curr;

				W+=curr;
				//WW+=kx?src->data[(idx-1)<<2|kc]:0;
			}
		}
	}
	free(row);
	//free(src);
	//*psrc=dst;
#endif
#if 0
	//reference cheating (W+E)/2  (fwd-only)
	Image *src=*psrc, *dst=0;
	image_copy(&dst, src);
	if(!dst)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(int kc=0;kc<4;++kc)
	{
		int depth=src->depth[kc], nlevels=1<<depth;
		if(!depth)
			continue;
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				int
					W=kx?dst->data[(idx-1)<<2|kc]:0,
					E=kx<src->iw-1?dst->data[(idx+1)<<2|kc]:0;
				int pred=(W+E)>>1;

				int curr=src->data[idx<<2|kc];
				curr-=pred;
				curr+=nlevels>>1;
				curr&=nlevels-1;
				curr-=nlevels>>1;
				src->data[idx<<2|kc]=curr;
			}
		}
	}
	free(dst);
#endif
}

#define MULTISTAGE_CALCCSIZE
#define LOAD(BUF, X, Y) ((unsigned)(ky+(Y))<(unsigned)src->ih&&(unsigned)(kx+(X))<(unsigned)src->iw?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
static void pred_multistage_stage1(Image *src, Image *dst, int fwd, int enable_ma, int kc, int *hist)
{
	int nlevels=1<<src->depth[kc];
	const int *pixels=fwd?src->data:dst->data;
//	const int *errors=fwd?dst->data:src->data;
#ifdef MULTISTAGE_CALCCSIZE
	memset(hist, 0, (nlevels+1LL)*sizeof(int));//
#endif
	for(int ky=0;ky<src->ih;ky+=2)
	{
		for(int kx=0;kx<src->iw;kx+=2)
		{
			int
				idx=src->iw*ky+kx,
				NNWW	=LOAD(pixels, -2, -2),
				NN	=LOAD(pixels,  0, -2),
				WW	=LOAD(pixels, -2,  0);
			int pred=NN+WW-NNWW;
			pred=MEDIAN3(NN, WW, pred);

			pred^=-fwd;
			pred+=fwd;
			pred+=src->data[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
#ifdef MULTISTAGE_CALCCSIZE
				++hist[pred+(nlevels>>1)+1];//
				++hist[0];
#endif
			}
			dst->data[idx<<2|kc]=pred;
		}
	}
}
#if 0
static void pred_multistage_stage2(Image *src, Image *dst, int fwd, int enable_ma, int kc, int *hist)
{
	int nlevels=1<<src->depth[kc];
	const int *pixels=fwd?src->data:dst->data;
//	const int *errors=fwd?dst->data:src->data;
#ifdef MULTISTAGE_CALCCSIZE
	memset(hist, 0, (nlevels+1LL)*sizeof(int));//
#endif
	for(int ky=0;ky<src->ih;ky+=2)
	{
		for(int kx=1;kx<src->iw;kx+=2)
		{
			int
				idx=src->iw*ky+kx,
			//	NNWW	=LOAD(pixels, -2, -2),
			//	NN	=LOAD(pixels,  0, -2),
			//	WW	=LOAD(pixels, -2,  0),
				W	=LOAD(pixels, -1,  0),
				E	=LOAD(pixels,  1,  0);
			//int pred2=NN+WW-NNWW;
			//pred2=MEDIAN3(NN, WW, pred2);
			//int pred=(7*(W+E)+2*pred2+8)>>4;
			int pred=(W+E+1)>>1;

			pred^=-fwd;
			pred+=fwd;
			pred+=src->data[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
#ifdef MULTISTAGE_CALCCSIZE
				++hist[pred+(nlevels>>1)+1];//
				++hist[0];
#endif
			}
			dst->data[idx<<2|kc]=pred;
		}
	}
}
static void pred_multistage_stage3(Image *src, Image *dst, int fwd, int enable_ma, int kc, int *hist)
{
	int nlevels=1<<src->depth[kc];
	const int *pixels=fwd?src->data:dst->data;
//	const int *errors=fwd?dst->data:src->data;
#ifdef MULTISTAGE_CALCCSIZE
	memset(hist, 0, (nlevels+1LL)*sizeof(int));//
#endif
	for(int ky=1;ky<src->ih;ky+=2)
	{
		for(int kx=0;kx<src->iw;kx+=2)
		{
			int
				idx=src->iw*ky+kx,
				N	=LOAD(pixels,  0, -1),
				S	=LOAD(pixels,  0,  1);
			//	NW	=LOAD(pixels, -1, -1),
			//	SW	=LOAD(pixels, -1,  1),
			//	NE	=LOAD(pixels,  1, -1),
			//	SE	=LOAD(pixels,  1,  1);
			int pred=(N+S+1)>>1;
			//int pred=(4*(N+S)-(NW+NE+SW+SE)+2)>>2;
			//int pred=(30*(N+S)+NW+SW+NE+SE+32)>>6;
			//int pred=(14*(N+S)+NW+SW+NE+SE+16)>>5;
			//int pred=(6*(N+S)+NW+SW+NE+SE+8)>>4;
			//int pred=(34*(N+S)-(NW+SW+NE+SE)+32)>>6;

			pred^=-fwd;
			pred+=fwd;
			pred+=src->data[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
#ifdef MULTISTAGE_CALCCSIZE
				++hist[pred+(nlevels>>1)+1];//
				++hist[0];
#endif
			}
			dst->data[idx<<2|kc]=pred;
		}
	}
}
static void pred_multistage_stage4(Image *src, Image *dst, int fwd, int enable_ma, int kc, int *hist)
{
	int nlevels=1<<src->depth[kc];
	const int *pixels=fwd?src->data:dst->data;
//	const int *errors=fwd?dst->data:src->data;
#ifdef MULTISTAGE_CALCCSIZE
	memset(hist, 0, (nlevels+1LL)*sizeof(int));//
#endif
	for(int ky=1;ky<src->ih;ky+=2)
	{
		for(int kx=1;kx<src->iw;kx+=2)
		{
			int
				idx=src->iw*ky+kx,
				N	=LOAD(pixels,  0, -1),
				S	=LOAD(pixels,  0,  1),
				W	=LOAD(pixels, -1,  0),
				E	=LOAD(pixels,  1,  0),
				NW	=LOAD(pixels, -1, -1),
				SW	=LOAD(pixels, -1,  1),
				NE	=LOAD(pixels,  1, -1),
				SE	=LOAD(pixels,  1,  1);
			int pred=(7*(N+S+W+E)+NW+SW+NE+SE+16)>>5;
			//int pred=(3*(N+S+W+E)+NW+SW+NE+SE+8)>>4;
			//int pred=(N+S+W+E+2)>>2;

			pred^=-fwd;
			pred+=fwd;
			pred+=src->data[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
#ifdef MULTISTAGE_CALCCSIZE
				++hist[pred+(nlevels>>1)+1];//
				++hist[0];
#endif
			}
			dst->data[idx<<2|kc]=pred;
		}
	}
}
#endif

static void pred_multistage_stage2c(Image *src, Image *dst, int fwd, int enable_ma, int kc, int *hist)
{
	int nlevels=1<<src->depth[kc], half=nlevels>>1;
	const int *pixels=fwd?src->data:dst->data;
//	const int *errors=fwd?dst->data:src->data;
#ifdef MULTISTAGE_CALCCSIZE
	memset(hist, 0, (nlevels+1LL)*sizeof(int));//
#endif
	for(int ky=1;ky<src->ih;ky+=2)
	{
		for(int kx=1;kx<src->iw;kx+=2)
		{
			int
				idx=src->iw*ky+kx,
				NW	=LOAD(pixels, -1, -1),
				SW	=LOAD(pixels, -1,  1),
				NE	=LOAD(pixels,  1, -1),
				SE	=LOAD(pixels,  1,  1),
			//	NWWW	=LOAD(pixels, -3, -1),
			//	NEEE	=LOAD(pixels,  3, -1),
			//	SWWW	=LOAD(pixels, -3,  1),
			//	SEEE	=LOAD(pixels,  3,  1),
			//	NNNW	=LOAD(pixels, -1, -3),
			//	NNNE	=LOAD(pixels,  1, -3),
			//	SSSW	=LOAD(pixels, -1,  3),
			//	SSSE	=LOAD(pixels,  1,  3),
				NNNWWW	=LOAD(pixels, -3, -3),
				SSSWWW	=LOAD(pixels, -3,  3),
				NNNEEE	=LOAD(pixels,  3, -3),
				SSSEEE	=LOAD(pixels,  3,  3);
			//	NNWW	=LOAD(pixels, -2, -2),
			//	NNEE	=LOAD(pixels,  2, -2),
			//	NN	=LOAD(pixels,  0, -2),
			//	WW	=LOAD(pixels, -2,  0);

			int pred=(9*(NW+NE+SW+SE)-(NNNWWW+NNNEEE+SSSWWW+SSSEEE)+16)>>5;
#if 0
			int
				dNW=abs(NW-NNNWWW),
				dNE=abs(NE-NNNEEE),
				dSW=abs(SW-SSSWWW),
				dSE=abs(SE-SSSEEE),
				sum=dNW+dNE+dSW+dSE;
			int pred;
			//if(kx==5&&ky==5)
			//	printf("");
			if(!sum)
				pred=NW;
			else
				pred=(
					dSE*(9*NW-(NWWW+NNNW)/2+NNNWWW)+
					dSW*(9*NE-(NEEE+NNNE)/2+NNNEEE)+
					dNE*(9*SW-(SWWW+SSSW)/2+SSSWWW)+
					dNW*(9*SE-(SEEE+SSSE)/2+SSSEEE)
				)/(9*sum);
				//pred=(dSE*(9*NW-NNNWWW)+dSW*(9*NE-NNNEEE)+dNE*(9*SW-SSSWWW)+dNW*(9*SE-SSSEEE))/(8*sum);
				//pred=(dSE*NW+dSW*NE+dNE*SW+dNW*SE)/sum;
#endif

			//int pred=(15*(NW+NE+SW+SE)-3*(NNNW+NNNE+SSSW+SSSE+NWWW+NEEE+SWWW+SEEE)+NNNWWW+NNNEEE+SSSWWW+SSSEEE+20)/40;

			//int pred=(6*(NW+NE+SW+SE)-(NNNW+NNNE+SSSW+SSSE+NWWW+NEEE+SWWW+SEEE+NNNWWW+SSSWWW+NNNEEE+SSSEEE))/12;
			//int pred=(NW+NE+SW+SE-(NN+WW))>>1;
			//int pred=(3*(NW+NE+SW+SE)+2*(NN+WW)+8)>>4;
			//int pred=(9*(NW+NE+SW+SE)-(NNNWWW+NNNEEE+SSSWWW+SSSEEE)+16)>>5;

			//int pred2=(2*(NW+NE+SW+SE)-(NNNWWW+SSSWWW+NNNEEE+SSSEEE)+2)>>2;
			//int pred=(4*(NW+NE+SW+SE)-(NNNW+NNNE+SSSW+SSSE+NWWW+NEEE+SWWW+SEEE)+4)>>3;
			//pred=(pred+pred2)>>1;

			//int pred3=NW+NE-NN;
			//pred3=MEDIAN3(NW, NE, pred3);
			//int pred2=NW+SW-WW;
			//pred2=MEDIAN3(NW, SW, pred2);
			//int pred=NN+WW-NNWW;
			//pred=MEDIAN3(NN, WW, pred);
			//pred=(2*pred+3*pred2+3*pred3+4)>>3;

			//int pN=(NW+NE)>>1, pW=(NW+SW)>>1;
			//int pred=pN+pW-NW;
			//pred=MEDIAN3(pN, pW, pred);
			//int pred;
			//if(abs(NW-SE)<abs(NE-SW))
			//	pred=(NW+SE)>>1;
			//else
			//	pred=(NE+SW)>>1;

			//int grad=NN+WW-NNWW;
			//int grad=WW+NNEE-NN;
			//grad=MEDIAN3(NN, WW, grad);
			//int vmin=WW, vmax=WW;
			//UPDATE_MIN(vmin, NN);
			//UPDATE_MAX(vmax, NN);
			//UPDATE_MIN(vmin, NNEE);
			//UPDATE_MAX(vmax, NNEE);

			//int pred=(27*(NW+SW+SE+NE)-(NWWW+NEEE+SWWW+SEEE+NNNW+SSSW+NNNE+SSSE+NNNWWW+NNNEEE+SSSWWW+SSSEEE))/96;
			//int pred=(NW+SW+NE+SE+2)>>2;
			//int grad=NN+WW-NNWW;
			//grad=MEDIAN3(NN, WW, grad);
			//int pred=(3*(NW+SW+NE+SE)+4*grad+8)>>4;
			//pred=(31*pred+grad+16)>>5;
			pred=CLAMP(-half, pred, half-1);
				
			pred^=-fwd;
			pred+=fwd;
			pred+=src->data[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
#ifdef MULTISTAGE_CALCCSIZE
				++hist[pred+(nlevels>>1)+1];//
				++hist[0];
#endif
			}
			dst->data[idx<<2|kc]=pred;
		}
	}
}
static void pred_multistage_stage3c(Image *src, Image *dst, int fwd, int enable_ma, int kc, int *hist)
{
	int nlevels=1<<src->depth[kc], half=nlevels>>1;
	const int *pixels=fwd?src->data:dst->data;
//	const int *errors=fwd?dst->data:src->data;
#ifdef MULTISTAGE_CALCCSIZE
	memset(hist, 0, (nlevels+1LL)*sizeof(int));//
#endif
	for(int ky=1;ky<src->ih;ky+=2)
	{
		for(int kx=0;kx<src->iw;kx+=2)
		{
			int
				idx=src->iw*ky+kx,
				N	=LOAD(pixels,  0, -1),
				S	=LOAD(pixels,  0,  1),
				W	=LOAD(pixels, -1,  0),
				E	=LOAD(pixels,  1,  0),
				NNW	=LOAD(pixels, -1, -2),
				NNE	=LOAD(pixels,  1, -2),
				SSW	=LOAD(pixels, -1,  2),
				SSE	=LOAD(pixels,  1,  2),
				NWW	=LOAD(pixels, -2, -1),
				NEE	=LOAD(pixels,  2, -1),
				SWW	=LOAD(pixels, -2,  1),
				SEE	=LOAD(pixels,  2,  1),
				WWW	=LOAD(pixels, -3,  0),
				EEE	=LOAD(pixels,  3,  0),
				NNN	=LOAD(pixels,  0, -3),
				SSS	=LOAD(pixels,  0,  3);
			int pred=(27*(N+W+S+E)-(WWW+EEE+NNN+SSS+NNW+SEE+NWW+SSE+NNE+SWW+SSW+NEE))/96;
			//int pred=(N+W+S+E+2)>>2;
			pred=CLAMP(-half, pred, half-1);
				
			pred^=-fwd;
			pred+=fwd;
			pred+=src->data[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
#ifdef MULTISTAGE_CALCCSIZE
				++hist[pred+(nlevels>>1)+1];//
				++hist[0];
#endif
			}
			dst->data[idx<<2|kc]=pred;
		}
	}
}
static void pred_multistage_stage4c(Image *src, Image *dst, int fwd, int enable_ma, int kc, int *hist)
{
	int nlevels=1<<src->depth[kc], half=nlevels>>1;
	const int *pixels=fwd?src->data:dst->data;
//	const int *errors=fwd?dst->data:src->data;
#ifdef MULTISTAGE_CALCCSIZE
	memset(hist, 0, (nlevels+1LL)*sizeof(int));//
#endif
	for(int ky=0;ky<src->ih;ky+=2)
	{
		for(int kx=1;kx<src->iw;kx+=2)
		{
			int
				idx=src->iw*ky+kx,
				N	=LOAD(pixels,  0, -1),
				S	=LOAD(pixels,  0,  1),
				W	=LOAD(pixels, -1,  0),
				E	=LOAD(pixels,  1,  0),
			//	NW	=LOAD(pixels, -1, -1),
			//	SW	=LOAD(pixels, -1,  1),
			//	NE	=LOAD(pixels,  1, -1),
			//	SE	=LOAD(pixels,  1,  1),
				NNW	=LOAD(pixels, -1, -2),
				NNE	=LOAD(pixels,  1, -2),
				SSW	=LOAD(pixels, -1,  2),
				SSE	=LOAD(pixels,  1,  2),
				NWW	=LOAD(pixels, -2, -1),
				NEE	=LOAD(pixels,  2, -1),
				SWW	=LOAD(pixels, -2,  1),
				SEE	=LOAD(pixels,  2,  1),
				WWW	=LOAD(pixels, -3,  0),
				EEE	=LOAD(pixels,  3,  0),
				NNN	=LOAD(pixels,  0, -3),
				SSS	=LOAD(pixels,  0,  3);
			//	NNNWWW	=LOAD(pixels, -3, -3),
			//	NNNEEE	=LOAD(pixels,  3, -3),
			//	SSSWWW	=LOAD(pixels, -3,  3),
			//	SSSEEE	=LOAD(pixels,  3,  3);
			int pred=(27*(N+W+S+E)-(WWW+EEE+NNN+SSS+NNW+SEE+NWW+SSE+NNE+SWW+SSW+NEE))/96;
			//int pred=(9*(N+W+S+E+NW+NE+SW+SE)-(NNN+WWW+SSS+EEE+NNNWWW+NNNEEE+SSSWWW+SSSEEE))>>6;
			//int pred=(N+W+S+E+2)>>2;
			//int pred=(7*(N+S+W+E)+NW+SW+NE+SE+16)>>5;
			pred=CLAMP(-half, pred, half-1);
				
			pred^=-fwd;
			pred+=fwd;
			pred+=src->data[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
#ifdef MULTISTAGE_CALCCSIZE
				++hist[pred+(nlevels>>1)+1];//
				++hist[0];
#endif
			}
			dst->data[idx<<2|kc]=pred;
		}
	}
}
void pred_multistage(Image *src, int fwd, int enable_ma)
{
#ifdef MULTISTAGE_CALCCSIZE
	double csizes[16]={0};
	int maxdepth=calc_maxdepth(src, 0), maxlevels=1<<maxdepth;
	int *hist=(int*)malloc((maxlevels+1LL)*sizeof(int));
#define CALC(X) csizes[kc<<2|X]=calc_entropy(hist+1, 1<<src->depth[kc], hist[0])
#else
#define CALC(X)
#endif
	Image *dst=0;
	image_copy_nodata(&dst, src);
	if(!dst
#ifdef MULTISTAGE_CALCCSIZE
		||!hist
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return;
	}
#ifdef MULTISTAGE_CALCCSIZE
	memset(dst->data, 0, (size_t)src->iw*src->ih*sizeof(int[4]));
#endif
	for(int kc=0;kc<4;++kc)
	{
		if(!src->depth[kc])
			continue;
		pred_multistage_stage1(src, dst, fwd, enable_ma, kc, hist);	CALC(0);
		pred_multistage_stage2c(src, dst, fwd, enable_ma, kc, hist);	CALC(1);
		pred_multistage_stage3c(src, dst, fwd, enable_ma, kc, hist);	CALC(2);
		pred_multistage_stage4c(src, dst, fwd, enable_ma, kc, hist);	CALC(3);
	}
#ifdef MULTISTAGE_CALCCSIZE
	//double csize=0;
	//for(int k=0;k<16;++k)//
	//{
	//	csizes[k]*=(double)(src->iw*src->ih>>2)/8;
	//	csize+=csizes[k];
	//}
	for(int k=0;k<16;++k)
		csizes[k]*=100./8;
	if(fwd&&GET_KEY_STATE(KEY_SHIFT))
	{
		int printed=snprintf(g_buf, G_BUF_SIZE,
			"Y0123 %7.4lf  %7.4lf %7.4lf %7.4lf %7.4lf\n"
			"U0123 %7.4lf  %7.4lf %7.4lf %7.4lf %7.4lf\n"
			"V0123 %7.4lf  %7.4lf %7.4lf %7.4lf %7.4lf\n",
			(csizes[ 0]+csizes[ 1]+csizes[ 2]+csizes[ 3])/4, csizes[ 0], csizes[ 1], csizes[ 2], csizes[ 3],
			(csizes[ 4]+csizes[ 5]+csizes[ 6]+csizes[ 7])/4, csizes[ 4], csizes[ 5], csizes[ 6], csizes[ 7],
			(csizes[ 8]+csizes[ 9]+csizes[10]+csizes[11])/4, csizes[ 8], csizes[ 9], csizes[10], csizes[11]
		);
		//int printed=snprintf(g_buf, G_BUF_SIZE,
		//	"%d | %d\n"
		//	"Y0123 %12.2lf | %12.2lf %12.2lf %12.2lf %12.2lf\n"
		//	"U0123 %12.2lf | %12.2lf %12.2lf %12.2lf %12.2lf\n"
		//	"V0123 %12.2lf | %12.2lf %12.2lf %12.2lf %12.2lf\n",
		//	src->iw*src->ih, src->iw*src->ih>>2,
		//	csizes[0]+csizes[1]+csizes[2]+csizes[3],
		//	csizes[0], csizes[1], csizes[2], csizes[3],
		//	csizes[4]+csizes[5]+csizes[6]+csizes[7],
		//	csizes[4], csizes[5], csizes[6], csizes[7],
		//	csizes[8]+csizes[9]+csizes[10]+csizes[11],
		//	csizes[8], csizes[9], csizes[10], csizes[11]
		//);
		copy_to_clipboard(g_buf, printed);
		messagebox(MBOX_OK, "Copied to clipboard", g_buf);
	}
	//set_window_title("Y0123 %lf %lf %lf %lf %lf",
	//	csizes[0]+csizes[1]+csizes[2]+csizes[3],
	//	csizes[0], csizes[1], csizes[2], csizes[3]
	//);
	free(hist);
#endif
#undef  CALC
	memcpy(src->data, dst->data, (size_t)src->iw*src->ih*sizeof(int[4]));
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}
	free(dst);
}
#undef  LOAD

//	#define P3_CALCCSIZE
#define P3_MAXPREDS 128
#define LOAD(BUF, X, Y) ((unsigned)(ky+(Y))<(unsigned)src->ih&&(unsigned)(kx+(X))<(unsigned)src->iw?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
//static Image *guide=0;
static void pred_separate_x(Image const *src, Image *dst, int fwd, int enable_ma, int kc
#ifdef P3_CALCCSIZE
	, int *hist, double *csizes
#endif
)
{
	int nlevels=1<<src->depth[kc];
//	int half=nlevels>>1;
	const int *pixels=fwd?src->data:dst->data;
//	const int *errors=fwd?dst->data:src->data;
#ifdef P3_CALCCSIZE
	memset(hist, 0, nlevels*sizeof(int[P3_MAXPREDS]));
#endif
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		long long errors[P3_MAXPREDS]={nlevels<<8};
		memfill(errors+1, errors, sizeof(errors)-sizeof(long long), sizeof(long long));
		for(int kx=0;kx<src->iw;++kx, ++idx)
		{
			long long lpred, lsum;
			int pred, curr;

			int
				NNNNN	=LOAD(pixels,  0, -5),
				NNNNWW	=LOAD(pixels, -2, -4),
				NNNNW	=LOAD(pixels, -1, -4),
				NNNN	=LOAD(pixels,  0, -4),
				NNNNE	=LOAD(pixels,  1, -4),
				NNNNEE	=LOAD(pixels,  2, -4),
			//	NNNNEEE	=LOAD(pixels,  3, -4),
				NNNNEEEE=LOAD(pixels,  4, -4),
			//	NNNWWWW	=LOAD(pixels, -4, -3),
				NNNWWW	=LOAD(pixels, -3, -3),
			//	NNNWW	=LOAD(pixels, -2, -3),
			//	NNNW	=LOAD(pixels, -1, -3),
				NNN	=LOAD(pixels,  0, -3),
			//	NNNE	=LOAD(pixels,  1, -3),
				NNNEE	=LOAD(pixels,  2, -3),
				NNNEEE	=LOAD(pixels,  3, -3),
				NNNEEEE	=LOAD(pixels,  4, -3),
				NNWWWW	=LOAD(pixels, -4, -2),
				NNWWW	=LOAD(pixels, -3, -2),
				NNWW	=LOAD(pixels, -2, -2),
				NNW	=LOAD(pixels, -1, -2),
				NN	=LOAD(pixels,  0, -2),
				NNE	=LOAD(pixels,  1, -2),
				NNEE	=LOAD(pixels,  2, -2),
				NNEEE	=LOAD(pixels,  3, -2),
				NNEEEE	=LOAD(pixels,  4, -2),
			//	NWWWWW	=LOAD(pixels, -5, -1),
				NWWWW	=LOAD(pixels, -4, -1),
				NWWW	=LOAD(pixels, -3, -1),
				NWW	=LOAD(pixels, -2, -1),
				NW	=LOAD(pixels, -1, -1),
				N	=LOAD(pixels,  0, -1),
				NE	=LOAD(pixels,  1, -1),
				NEE	=LOAD(pixels,  2, -1),
				NEEE	=LOAD(pixels,  3, -1),
				NEEEE	=LOAD(pixels,  4, -1),
				NEEEEE	=LOAD(pixels,  5, -1),
				NEEEEEE	=LOAD(pixels,  6, -1),
				NEEEEEEE	=LOAD(pixels,  7, -1),
				NEEEEEEEE	=LOAD(pixels,  8, -1),
				WWWWWWWWW	=LOAD(pixels, -9,  0),
				WWWWWWWW	=LOAD(pixels, -8,  0),
				WWWWWWW	=LOAD(pixels, -7,  0),
				WWWWWW	=LOAD(pixels, -6,  0),
				WWWWW	=LOAD(pixels, -5,  0),
				WWWW	=LOAD(pixels, -4,  0),
				WWW	=LOAD(pixels, -3,  0),
				WW	=LOAD(pixels, -2,  0),
				W	=LOAD(pixels, -1,  0);
			int preds[]=//P3
			{
				 0x2E00000, N+W-NW,
				//0x0100000, (2*N+W+NE-NW-NNE)>>1,//X
				 0x0400000, N+W-((NW+NN+WW+NE)>>2),
				-0x0300000, ((2*(N+W)-(NW+NN+WW+NE))*9+(WWW+NWW+NNW+NNN+NNE+NEE)*2)/12,
				 0x0C00000, 2*(N+W-NW)-(NN+WW-NNWW),
				 0x0200000, 3*(N+W-NW-(NN+WW-NNWW))+NNN+WWW-NNNWWW,
				-0x0300000, N+W-NW+NEE-(NE+NNEE-NNE),
				//-0x0100000, N+W-NW+WW-(WWW+NWW-NWWW),//X
				 0x0100000, N+W-NW+WW+NEE-(WWW+NWW-NWWW+NE+NNEE-NNE),
				//0x0100000, N+W-NW+((NE-(N+NNE-NN))<<1),//X
				-0x0800000, (N+W)>>1,
				-0x0100000, N+W-((NN+WW)>>1),
				 0x0200000, (4*(N+W)+NE-NW)>>2,//shift right 3?
				 0x0800000, (4*(N+W)-NE+NW)>>2,//shift right 3?
				 0x0800000, (4*(N+W)-(NE+NW))>>2,
				-0x0200000, (4*(N+W)-(NE+NW))/6,
				-0x0100000, (4*N-(NW+NE))>>2,
				-0x0E00000, (W+NE)>>1,
				-0x0E00000, (2*W+NEE)/3,
				 0x0F00000, (3*W+NEEE)>>2,
				 0x1000000, 4*W+N+NE+NEE+NEEE,//shift right 3?
				-0x0500000, (0x289*W+0x080*N+0x06A*NE+0x04F*NEE+0x03E*NEEE)>>10,//proportional to distance
				-0x0500000, (192*W+33*N+22*NE+11*NEE+6*NEEE)/264,//proportional to square distance
				 0x0500000, W+NE-N,
				 0x0200000, 2*(W+NE-N)-(WW+NNEE-NN),
				 0x0200000, 3*(W+NE-N-(WW+NNEE-NN))+WWW+NNNEEE-NNN,
				 0x0500000, W+NE-N+NW-(NWW+NN-NNW),
				-0x0100000, W+NE-((3*N+NW)>>2),
				 0x0800000, (2*W+NEE-N)>>1,
				 0x0800000, (3*W+NEEE-N)>>1,
				 0x1400000, (3*(3*W+NE+NEE)-10*N)/5,
				 0x0300000, (WWWWWWWWW+WWWWWWWW)>>1,
				 0x0400000, (WWWWWWW+WWWWWW)>>1,
				-0x0100000, (WWWWW+WWWW+NWWW+NEEE)>>2,
				 0x0200000, (NWWWW+NWWW+NEEEE+NEEEEE)>>2,
				 0x0200000, W+NW-NWW,
				 0x0700000, NW+NWW-NNWWW,
				 0x0600000, N+NW-NNW,
				 0x1800000, N+NE-NNE,
				-0x0500000, (9*(NW-NNWW)+NNNNWW+NNNWWW+NNWWWW)/3,
				 0x0300000, (9*(NE-NNEE)+NNNNEE+NNNEEE+NNEEEE)/3,
				 0x0500000, NN+NNW-(NNNN+NNNNW+NNNNWW)/3,
				 0x0500000, 2*NN-(NNNNW+NNNN+NNNNE)/3,
				 0x0500000, NN+NNE-(NNNN+NNNNE+NNNNEE)/3,
				 0x0B00000, NNEE+NEE-(NNNNEEEE+NNNEEEE+NNEEEE)/3,
				 0x0500000, 3*NE-(NNE+NEE),
				//-0x0700000, NE+NEE-(NNEE+NNEEE+NNEEEE)/3,
				-0x0600000, (14*NE-(NNEE+NNNEE+NNEEE))/11,
				-0x0500000, ((NW+N+W)*10-(NNWW+NNW+NN+NWW+WW)*3)/15,
				-0x0300000, (75*W+53*NW)>>7,
				-0x0300000, (75*N+53*NW)>>7,
				-0x0300000, (75*N+53*NE)>>7,
				 0x0F00000, (2*(N+NE+NEE)-(NNE+NNEE+NNEEE))/3,
				 0x0F00000, NE+NEE-(NNEE+NNEEE+NNEEEE)/3,
				 0x1400000, (NEEE+NEEEE)>>1,
				 0x0E00000, (NEEEEE+NEEEEEE)>>1,
				 0x0C00000, (NEEEEEEE+NEEEEEEEE)>>1,
				 0x1800000, (4*N-2*NN+NW+NE)>>2,
				 0x0800000, (N+W+2*(NW+NE)-(NNWW+NNEE))>>2,
				
				 0x0A00000, 2*N-W,
				 0x0A00000, 2*W-N,
				-0x0500000, (2*(NW+2*N+NE)-(NNW+2*NN+NNE))>>2,
				 0x0500000, 4*(N-W)-NW,
				 0x0100000, N+W+NW+NE,
				-0x0400000, 2*NE-NN,
				-0x0400000, 2*NE-NNE,
				-0x0600000, (2*(W+NE)-(WW+NNEE))>>1,
				//0x0100000, 2*NNW-NNNWW,

				0x0F00000, N,
				0x1400000, 2*N-NN,
				0x1000000, 3*(N-NN)+NNN,
				0x1800000, 4*(N+NNN)-6*NN-(NNNNW+NNNN+NNNNE)/3,
				0x1000000, (N+NN)>>1,
				0x0800000, (N+NN+NNN)>>1,		//why 3/2?
				0x1200000, (NN+NNN+NNNN+NNNNN)>>2,
				
				0x0F00000, W,
				0x1400000, 2*W-WW,
				0x1000000, 3*(W-WW)+WWW,
				0x1800000, 4*(W+WWW)-6*WW-(WWWW+NWWWW)/2,
				0x1000000, (W+WW)>>1,
				0x0800000, (W+WW+WWW)>>1,		//why 3/2?
				0x1200000, (WW+WWW+WWWW+WWWWW)>>2,
				
				-0x0800000, NW,
				 0x0400000, 2*NW-NNWW,
				 0x0200000, 3*(NW-NNWW)+NNNWWW,
				 0x0400000, (NW+NNWW)>>1,
				 0x0400000, (NW+NNWW+NNNWWW)>>1,	//why 3/2?
				-0x0300000, (NW+NNW+NWW+NNWW)>>1,
				
				-0x0800000, NE,
				 0x0400000, 2*NE-NNEE,
				 0x0200000, 3*(NE-NNEE)+NNNEEE,
				 0x0400000, (NE+NNEE)>>1,
				 0x0400000, (NE+NNEE+NNNEEE)>>1,	//why 3/2?
				-0x0300000, (NE+NNE+NEE+NNEE)>>1,
#if 0
				N+W-NW,
				(4*(N+W)-(NE+NW))>>2,
				(2*W+NEE-N)>>1,
				(3*W+NEEE)>>2,
				(3*(3*W+NE+NEE)-10*N+2)/5,

				//directional predictors:
				//-0x0100000, (3*(N+W-(NN+WW))+NNN+WWW)>>1,//X
				//0x0100000, (14*(NE+NNEE)-3*(NNNNEEE+NNNNEEEE+NNNEEE+NNNEEEE))>>4,//X
				//0x0100000, (2*(9*W+NW)-(NWWWWW+NWWWW+WWWWW+WWWW))>>4,//X
				//-0x0100000, W+WW-((NWWWW+NWWW+WWWW+WWW)>>2),//X
				//2*NE-((NNEE+NNEEE+NNNEE+NNNEEE)>>2),
				//2*NW-((NNWW+NNWWW+NNNWW+NNNWWW)>>2),
				//(10*NE-(NNEE+NNEEE+NNNEE+NNNEEE))/6,
				//WWWWWW,
				//(W+NW+N)/3,
				//(N+NE+NEE)/3,
				//(WWWW+WWWWW+WWWWWW)/3,
				//(WWWWWWW+WWWWWW)>>1,
				//W,			//180.00 degrees
				//((WW+W+NW+NWW)>>1)-NWWW,//165.96 degrees
				//W+NW-NWW,		//153.43 degrees
				//NW+NWW-NNWWW,		//123.69 degrees
				//NW,			//135.00 degrees
				//N+NW-NNW,		//116.57 degrees
				//N,			// 90.00 degrees
				N+NE-NNE,		// 63.43 degrees
				(2*(N+NE+NEE)-(NNE+NNEE+NNEEE))/3,
				//N+(N+NE+NEE)/3-(NNNE+NNNEE+NNNEEE+NNE+NNEE+NNEEE+NE+NEE+NEEE)/9,
				//(10*(N+NE+NEE)-3*(NNE+NNEE+NNEEE+NEE+NNNEE))/15,
				//NE+NNE-NNNEE,		// 56.31 degrees
				//NE,			// 45.00 degrees
				NE+NEE-(NNEE+NNEEE+NNEEEE)/3,
				//NE+NEE-(NEEE+NNEE+NNEEE+NNEEEE+NNNEEE)/5,//X
				//NE+NEE-NNEEE,		// 33.69 degrees
				//NEE,
				//NEEEE,		// 14.04 degrees
				(NEEE+NEEEE)>>1,
				(NEEEEE+NEEEEEE)>>1,
				//(NEEEEEE+NEEEEEEE)>>1,
				//(NEEEEEEE+NEEEEEEEE)>>1,
				//NEEEEEEEE,

				(4*W+N+NE+NEE+NEEE)>>2,
				(4*N-2*NN+NW+NE)>>2,
				(N+W+2*(NW+NE)-(NNWW+NNEE))>>2,

				//W,
				//(75*W+53*NW)>>7,
				//(75*N+53*NW)>>7,
				//N,
				//(75*N+53*NE)>>7,

				//N,
				2*N-NN,
				//5*N-(NN+NNN),
				3*(N-NN)+NNN,
				//3*(N-NN)+(NNN+NNNE)/2,
				//3*(N-NN)+(NNNW+NNN+NNNE)/3,
				//4*(N+NNN)-6*NN-NNNN,
				4*(N+NNN)-6*NN-(NNNNW+NNNN+NNNNE)/3,
				(N+NN)>>1,
				(N+NN+NNN)>>1,		//why 3/2?
				(NN+NNN+NNNN+NNNNN)>>2,

				//W,
				2*W-WW,
				3*(W-WW)+WWW,
				//3*(W-WW)+(WWW+WWWW)/2,
				4*(W+WWW)-6*WW-(WWWW+NWWWW)/2,
				(W+WW)>>1,
				(W+WW+WWW)>>1,		//why 3/2?
				(WW+WWW+WWWW+WWWWW)>>2,
				
				//NW,
				//2*NW-NNWW,
				//3*(NW-NNWW)+NNNWWW,
				//(NW+NNWW)>>1,
				//(NW+NNWW+NNNWWW)>>1,		//why 3/2?
				
				//NE,
				//2*NE-NNEE,
				//3*(NE-NNEE)+NNNEEE,
				//(NE+NNEE)>>1,
				//(NE+NNEE+NNNEEE)>>1,		//why 3/2?

				//negative effect:
				//0,
				//W+NE-N,		//superseded
				//(6*(N+W)-(NW+NWW+NNW+NNWW))>>3,
				//(4*(N+W)+NE-NW)>>3,	//superseded
				//(N+W)>>1,		//redundant bias
				//NE+WW-NW,
				//W+NEE-NE,
				//(WW+W+NEE-NW)>>1,
				//(9*W+4*NEE)/13,
				//(7*W+NEEE+NEEEE)/9,
				//NW+NWW-NNWWW,
				//(4*(N+W+NW+NE)-(NN+WW+NNWW+NNEE)+6)/12,
				//floor_sqrt(((long long)N+half)*((long long)W+half)),
				//5*(N-NNNN+2*(NNN-NN))+NNNNN,
				//5*(W-WWWW+2*(WWW-WW))+WWWWW,
				//(N+NN+NNN+NNNN+NNNNN)/5,
				//(W+WW+WWW+WWWW+WWWWW)/5,
#endif
			};
			//preds[1]=MEDIAN3(N, W, preds[1]);
#if 0
			int preds[]=
			{
				-W,
				0,
				W,
				2*W-WW,
				3*(W-WW)+WWW,
				//4*(W+WWW)-6*WW-WWWW,
				(W+WW)>>1,
				(W+WW+WWW)>>1,
				//(W+WW+WWW+WWWW)>>2,
			};
#endif
			//ALIGN(16) int cmin[4], cmax[4];

			//if(kc==1&&ky==5&&kx==6)
			//	printf("");
			lpred=0, lsum=0;
			for(int k=0;k<(int)_countof(preds)/2;++k)
			{
				long long weight=((long long)preds[k<<1|0]<<8)/errors[k];
				lpred+=weight*preds[k<<1|1];
				lsum+=weight;
			}
			pred=(int)(lpred/lsum);
			CLAMP3_32(pred, pred, N, W, NE);
#if 0
			{
				__m128i mpred	=_mm_set_epi32(0, 0, 0, pred);
				__m128i mW	=_mm_set_epi32(0, 0, 0, W);
				__m128i mN	=_mm_set_epi32(0, 0, 0, N);
				__m128i mNE	=_mm_set_epi32(0, 0, 0, NE);
				__m128i mmin=_mm_min_epi32(mN, mW);
				__m128i mmax=_mm_max_epi32(mN, mW);
				mmin=_mm_min_epi32(mmin, mNE);
				mmax=_mm_max_epi32(mmax, mNE);
				mpred=_mm_max_epi32(mpred, mmin);
				mpred=_mm_min_epi32(mpred, mmax);
				{
					ALIGN(16) int apred[4];
					_mm_store_si128((__m128i*)apred, mpred);
					pred=apred[0];
				}

				//cmin=W, cmax=W;
				//UPDATE_MIN(cmin, N);
				//UPDATE_MAX(cmax, N);
				//UPDATE_MIN(cmin, NE);
				//UPDATE_MAX(cmax, NE);
				//pred=CLAMP(cmin[0], pred, cmax[0]);
			}
#endif

			pred^=-fwd;
			pred+=fwd;
			pred+=src->data[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
			}
			dst->data[idx<<2|kc]=pred;

			//if(!fwd&&guide)
			//{
			//	if(dst->data[idx<<2|kc]!=guide->data[idx<<2|kc])
			//		LOG_ERROR("");
			//}
			//const int count=(int)_countof(preds)/2;//

			curr=pixels[idx<<2|kc];
			for(int k=0;k<(int)_countof(preds)/2;++k)
			{
				long long e=errors[k];
				e+=(((long long)abs(curr-preds[k<<1|1])<<8)-e)/3;
				UPDATE_MAX(e, 1);
				errors[k]=e;
#ifdef P3_CALCCSIZE
				e=curr-preds[k<<1|1];
				e+=nlevels>>1;
				e&=nlevels-1;
				++hist[k<<src->depth[kc]|(int)e];
#endif
			}
		}
	}
#ifdef P3_CALCCSIZE
	for(int k=0;k<P3_MAXPREDS;++k)
		csizes[k<<2|kc]=calc_entropy(hist+((size_t)k<<src->depth[kc]), nlevels, src->iw*src->ih)*100./8;
#endif
}
#if 0
static void pred_separate_y(Image const *src, Image *dst, int fwd, int enable_ma, int kc)
{
	int nlevels=1<<src->depth[kc], half=nlevels>>1;
	const int *pixels=fwd?src->data:dst->data;
//	const int *errors=fwd?dst->data:src->data;
#if 1
#define LG_PREDLEVELS 7
	int preds[1<<LG_PREDLEVELS];
	for(int k=0;k<(int)_countof(preds);++k)
		preds[k]=(k*nlevels>>LG_PREDLEVELS)-half;
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		long long errors[1<<LG_PREDLEVELS]={1};
		memfill(errors+1, errors, sizeof(errors)-sizeof(long long), sizeof(long long));
		for(int kx=0;kx<src->iw;++kx, ++idx)
		{
			int
				//NW	=LOAD(pixels, -1, -1),
				N	=LOAD(pixels,  0, -1),
				NE	=LOAD(pixels,  1, -1),
				W	=LOAD(pixels, -1,  0);
			int vmin=W, vmax=W;
			UPDATE_MIN(vmin, N);
			UPDATE_MAX(vmax, N);
			UPDATE_MIN(vmin, NE);
			UPDATE_MAX(vmax, NE);
			long long lpred=0, lsum=0;
			for(int k=0;k<(int)_countof(errors);++k)
			{
				//int p=(k*nlevels>>7)-half;
				long long w=(((long long)kx+1)<<24)/errors[k];
				lpred+=w*preds[k];
				lsum+=w;
			}
			int pred=(int)(lpred/lsum);
			pred=CLAMP(vmin, pred, vmax);
			
			pred^=-fwd;
			pred+=fwd;
			pred+=src->data[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
			}
			dst->data[idx<<2|kc]=pred;

			int curr=pixels[idx<<2|kc];
			for(int k=0;k<(int)_countof(errors);++k)
			{
				long long e=errors[k];
				e+=((abs(curr-preds[k])<<8)-e)/3;
				UPDATE_MAX(e, 1);
				errors[k]=e;
			}
		}
	}
#endif
#if 0
	for(int kx=0, idx=0;kx<src->iw;++kx)
	{
		long long errors[128]={1};
		memfill(errors+1, errors, sizeof(errors)-sizeof(long long), sizeof(long long));
		for(int ky=0;ky<src->ih;++ky, ++idx)
		{
			int
				WWWW	=LOAD(pixels, 0, -4),
				WWW	=LOAD(pixels, 0, -3),
				WW	=LOAD(pixels, 0, -2),
				W	=LOAD(pixels, 0, -1);
			int preds[]=
			{
				-W,
				0,
				W,
				2*W-WW,
				3*(W-WW)+WWW,
				//4*(W+WWW)-6*WW-WWWW,
				(W+WW)>>1,
				(W+WW+WWW)>>1,
				//(W+WW+WWW+WWWW)>>2,
			};
			long long lpred=0, lsum=0;
			for(int k=0;k<(int)_countof(preds);++k)
			{
				long long w=(((long long)ky+1)<<16)/errors[k];
				lpred+=w*preds[k];
				lsum+=w;
			}
			int pred=(int)(lpred/lsum);

			pred^=-fwd;
			pred+=fwd;
			pred+=src->data[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
			}
			dst->data[idx<<2|kc]=pred;

			int curr=pixels[idx<<2|kc];
			for(int k=0;k<(int)_countof(preds);++k)
				errors[k]+=abs(curr-preds[k]);
		}
	}
#endif
}
#endif
#undef  LOAD
void pred_separate(Image *src, int fwd, int enable_ma)
{
	double t=time_sec();
	Image *dst;
#ifdef P3_CALCCSIZE
	int maxdepth, *hist;
	double *csizes;
	char *str;
#endif

	//if(fwd)
	//	image_copy(&guide, src);
	dst=0;
	image_copy_nodata(&dst, src);
#ifdef P3_CALCCSIZE
	maxdepth=calc_maxdepth(src, 0);
	hist=(int*)malloc(sizeof(int[P3_MAXPREDS])<<maxdepth);
	csizes=(double*)malloc(sizeof(double[P3_MAXPREDS*4]));
	str=(char*)malloc(0x100000);
#endif
	if(!dst
#ifdef P3_CALCCSIZE
		||!hist||!csizes||!str
#endif
	)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(dst->data, 0, (size_t)src->iw*src->ih*sizeof(int[4]));
	for(int kc=0;kc<4;++kc)
	{
		if(!src->depth[kc])
			continue;
		//if(fwd)
		//{
		//	pred_separate_x(src, dst, fwd, enable_ma, kc);
		//	pred_separate_y(dst, src, fwd, enable_ma, kc);
		//}
		//else
		//{
		//	pred_separate_y(src, dst, fwd, enable_ma, kc);
		//	pred_separate_x(dst, src, fwd, enable_ma, kc);
		//}
		pred_separate_x(src, dst, fwd, enable_ma, kc
#ifdef P3_CALCCSIZE
			, hist, csizes
#endif
		);
	}
	memcpy(src->data, dst->data, (size_t)src->iw*src->ih*sizeof(int[4]));
#ifdef P3_CALCCSIZE
	if(GET_KEY_STATE(KEY_SHIFT)&&fwd)
	{
		int printed=0;
		for(int k=0;k<P3_MAXPREDS;++k)
		{
			printed+=snprintf(str+printed, 0x100000LL-printed-1, "%d", k);
			for(int kc=0;kc<4;++kc)
			{
				if(!src->depth[kc])
					continue;
				printed+=snprintf(str+printed, 0x100000LL-printed-1, "\t%8.4lf", csizes[k<<2|kc]);
			}
			printed+=snprintf(str+printed, 0x100000LL-printed-1, "\n");
		}
		copy_to_clipboard(str, printed);
		messagebox(MBOX_OK, "Copied to clipboard", "Too much to print here...");
	}
	free(hist);
	free(csizes);
	free(str);
#endif
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}
	free(dst);
	t=time_sec()-t;
	if(loud_transforms)
	{
		if(fn)
			set_window_title("%s - %lf sec", (char*)fn->data, t);
		else
			set_window_title("%lf sec", t);
	}
}

void pred_dir(Image *src, int fwd, int enable_ma)
{
	Image *dst;
	int maxdepth, ssesize, *ssearr;
	const int *pixels;
//	const int *errors;

	dst=0;
	image_copy(&dst, src);
	maxdepth=calc_maxdepth(src, 0);
	ssesize=(1<<maxdepth)<<1;
	ssearr=(int*)malloc(ssesize*sizeof(int));
	if(!ssearr||!dst)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	pixels=fwd?src->data:dst->data;
//	errors=fwd?dst->data:src->data;
	for(int kc=0;kc<4;++kc)
	{
		int nlevels;

		if(!src->depth[kc])
			continue;
		memset(ssearr, 0, ssesize*sizeof(int));
		nlevels=1<<src->depth[kc];
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				int pred, curr, arg;
#define LOAD(X, Y) ((unsigned)(ky+(Y))<(unsigned)src->ih&&(unsigned)(kx+(X))<(unsigned)src->iw?pixels[(idx+src->iw*(Y)+(X))<<2|kc]:0)
				int nb[]=
				{
					//LOAD(-2, -2),
					//LOAD(-1, -2),
					//LOAD( 0, -2),
					//LOAD( 1, -2),
					//LOAD( 2, -2),
					//LOAD(-2, -1),
					//LOAD(-1, -1),
					//LOAD( 0, -1),
					//LOAD( 1, -1),
					//LOAD( 2, -1),
					//LOAD(-2,  0),
					//LOAD(-1,  0),

					LOAD(-1,  0),//W	180
					LOAD(-1, -1),//NW	135
					LOAD( 0, -1),//N	90
					LOAD( 1, -1),//NE	45
				};
#undef  LOAD
				curr=pixels[idx<<2|kc];
				arg=0;
				for(int k=1;k<(int)_countof(nb);++k)
				{
					if(abs(curr-nb[arg])>abs(curr-nb[k]))
						arg=k;
				}
				pred=nb[arg];
				
				pred^=-fwd;
				pred+=fwd;
				pred+=src->data[idx<<2|kc];
				if(enable_ma)
				{
					pred+=nlevels>>1;
					pred&=nlevels-1;
					pred-=nlevels>>1;
				}
				dst->data[idx<<2|kc]=pred;
				//dst->data[idx<<2|kc]=127-arg*64;
#if 0
				int
					N =ky    ?pixels[(idx-src->iw  )<<2|kc]:0,
					W =kx    ?pixels[(idx        -1)<<2|kc]:0,
					NW=kx&&ky?pixels[(idx-src->iw-1)<<2|kc]:0,
					NE=kx+1<src->iw&&ky?pixels[(idx-src->iw+1)<<2|kc]:0;

				int
					eN =ky    ?errors[(idx-src->iw  )<<2|kc]:0,
					eW =kx    ?errors[(idx        -1)<<2|kc]:0,
					eNW=kx&&ky?errors[(idx-src->iw-1)<<2|kc]:0,
					eNE=kx+1<src->iw&&ky?errors[(idx-src->iw-1)<<2|kc]:0;

				//int eprev[]=
				//{
				//	kc>0?errors[idx<<2|0]:0,
				//	kc>1?errors[idx<<2|1]:0,
				//	kc>2?errors[idx<<2|2]:0,
				//};
				long long weights[]=
				{
					0x10000LL/(abs(eN)+1),
					0x10000LL/(abs(eW)+1),
					0x10000LL/(abs(eNW)+1),
					0x10000LL/(abs(eNE)+1),
				};
				if(kx==100&&ky==100)
					printf("");
				long long sum=weights[0]+weights[1]+weights[2]+weights[3];
				int pred=(int)((weights[0]*N+weights[1]*W+weights[2]*NW+weights[3]*NE+(sum>>1))/sum);

				int *curr_cell=ssearr+(eN+eW+nlevels);
				//int *curr_cell=ssearr+(((eN+0+(nlevels>>1))<<src->depth[kc]|(eW+0+(nlevels>>1)))&((1<<(src->depth[kc]<<1))-1));
				int sse_sum=*curr_cell>>12, sse_count=*curr_cell&0xFFF;
				//int *cell[]=
				//{
				//	//ssearr+(((eN-1+(nlevels>>1))<<src->depth[kc]|(eW-1+(nlevels>>1)))&((1<<(src->depth[kc]<<1))-1)),
				//	ssearr+(((eN-1+(nlevels>>1))<<src->depth[kc]|(eW+0+(nlevels>>1)))&((1<<(src->depth[kc]<<1))-1)),
				//	//ssearr+(((eN-1+(nlevels>>1))<<src->depth[kc]|(eW+1+(nlevels>>1)))&((1<<(src->depth[kc]<<1))-1)),
				//	ssearr+(((eN+0+(nlevels>>1))<<src->depth[kc]|(eW-1+(nlevels>>1)))&((1<<(src->depth[kc]<<1))-1)),
				//	curr_cell,
				//	ssearr+(((eN+0+(nlevels>>1))<<src->depth[kc]|(eW+1+(nlevels>>1)))&((1<<(src->depth[kc]<<1))-1)),
				//	//ssearr+(((eN+1+(nlevels>>1))<<src->depth[kc]|(eW-1+(nlevels>>1)))&((1<<(src->depth[kc]<<1))-1)),
				//	ssearr+(((eN+1+(nlevels>>1))<<src->depth[kc]|(eW+0+(nlevels>>1)))&((1<<(src->depth[kc]<<1))-1)),
				//	//ssearr+(((eN+1+(nlevels>>1))<<src->depth[kc]|(eW+1+(nlevels>>1)))&((1<<(src->depth[kc]<<1))-1)),
				//};
				//int sse_sum=0, sse_count=0;
				//for(int k=0;k<(int)_countof(cell);++k)
				//{
				//	sse_sum+=*cell[k]>>12;
				//	sse_count+=*cell[k]&0xFFF;
				//}
				if(sse_count)
					pred+=sse_sum/sse_count;

				pred=MEDIAN3(N, W, pred);
				
				pred^=-fwd;
				pred+=fwd;
				pred+=src->data[idx<<2|kc];
				if(enable_ma)
				{
					pred+=nlevels>>1;
					pred&=nlevels-1;
					pred-=nlevels>>1;
				}
				dst->data[idx<<2|kc]=pred;

				int error=errors[idx<<2|kc];
				++sse_count;
				sse_sum+=error;
				if(sse_count>640)
				{
					sse_count>>=1;
					sse_sum>>=1;
				}
				*curr_cell=sse_sum<<12|sse_count;
#endif
			}
		}
	}
	memcpy(src->data, dst->data, (size_t)src->iw*src->ih*sizeof(int[4]));
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}
	free(dst);
}

static int invert_matrix(double *matrix, int size, double *temprow)
{
	int success;

	success=1;
	for(int it=0;it<size;++it)
	{
		int kp;
		double pivot;

		kp=it;
		for(;kp<size;++kp)
		{
			if(fabs(matrix[((size_t)size<<1)*kp+it])>1e-6)
				break;
		}
		if(kp==size)
		{
			success=0;
			break;
		}
		if(kp!=it)
		{
			memcpy(temprow, matrix+((size_t)size<<1)*it, size*sizeof(double[2]));
			memcpy(matrix+((size_t)size<<1)*it, matrix+((size_t)size<<1)*kp, size*sizeof(double[2]));
			memcpy(matrix+((size_t)size<<1)*kp, temprow, size*sizeof(double[2]));
		}
		pivot=matrix[((size_t)size<<1)*it+it];
		for(int kx=it;kx<(size<<1);++kx)
			matrix[((size_t)size<<1)*it+kx]/=pivot;
		for(int ky=0;ky<size;++ky)
		{
			double factor;

			if(ky==it)
				continue;
			factor=matrix[((size_t)size<<1)*ky+it];
			if(fabs(factor)>1e-6)
			{
				for(int kx=it;kx<(size<<1);++kx)
					matrix[((size_t)size<<1)*ky+kx]-=matrix[((size_t)size<<1)*it+kx]*factor;
			}
		}
	}
	//_freea(temp);
	return success;
}
#define OLS_NPARAMS 4	//8
#define OLS_NSAMPLES 17
void pred_ols(Image *src, int fwd, int enable_ma)
{
	int successcount;
	Image *dst;
	const int *pixels;
//	const int *errors;
	int c1_init[4]={0};
	double c1_params[OLS_NPARAMS*4]={0};

	successcount=0;
	dst=0;
	image_copy(&dst, src);
	if(!dst)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	pixels=fwd?src->data:dst->data;
//	errors=fwd?dst->data:src->data;
	for(int kc=0;kc<4;++kc)
	{
		if(!src->depth[kc])
			continue;
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				int pred, success;
#define LOAD(BUF, X, Y) ((unsigned)(ky-(Y))<(unsigned)src->ih&&(unsigned)(kx-(X))<(unsigned)src->iw?BUF[(src->iw*(ky-(Y))+kx-(X))<<2|kc]:0)
				int
					NNNWWWW	=LOAD(pixels,  4, 3),
					NNNWWW	=LOAD(pixels,  3, 3),
					NNNWW	=LOAD(pixels,  2, 3),
					NNNW	=LOAD(pixels,  1, 3),
					NNN	=LOAD(pixels,  0, 3),
					NNNE	=LOAD(pixels, -1, 3),
					NNNEE	=LOAD(pixels, -2, 3),
					NNNEEE	=LOAD(pixels, -3, 3),
					NNNEEEE	=LOAD(pixels, -4, 3),
					NNWWWW	=LOAD(pixels,  4, 2),
					NNWWW	=LOAD(pixels,  3, 2),
					NNWW	=LOAD(pixels,  2, 2),
					NNW	=LOAD(pixels,  1, 2),
					NN	=LOAD(pixels,  0, 2),
					NNE	=LOAD(pixels, -1, 2),
					NNEE	=LOAD(pixels, -2, 2),
					NNEEE	=LOAD(pixels, -3, 2),
					NNEEEE	=LOAD(pixels, -4, 2),
					NWWWW	=LOAD(pixels,  4, 1),
					NWWW	=LOAD(pixels,  3, 1),
					NWW	=LOAD(pixels,  2, 1),
					NW	=LOAD(pixels,  1, 1),
					N	=LOAD(pixels,  0, 1),
					NE	=LOAD(pixels, -1, 1),
					NEE	=LOAD(pixels, -2, 1),
					NEEE	=LOAD(pixels, -3, 1),
				//	NEEEE	=LOAD(pixels, -4, 1),
					WWWW	=LOAD(pixels,  4, 0),
					WWW	=LOAD(pixels,  3, 0),
					WW	=LOAD(pixels,  2, 0),
					W	=LOAD(pixels,  1, 0);
#if 0
				int
					eNNNWWWW	=LOAD(errors,  4, 3),
					eNNNWWW		=LOAD(errors,  3, 3),
					eNNNWW		=LOAD(errors,  2, 3),
					eNNNW		=LOAD(errors,  1, 3),
					eNNN		=LOAD(errors,  0, 3),
					eNNNE		=LOAD(errors, -1, 3),
					eNNNEE		=LOAD(errors, -2, 3),
					eNNNEEE		=LOAD(errors, -3, 3),
					eNNNEEEE	=LOAD(errors, -4, 3),
					eNNWWWW		=LOAD(errors,  4, 2),
					eNNWWW		=LOAD(errors,  3, 2),
					eNNWW		=LOAD(errors,  2, 2),
					eNNW		=LOAD(errors,  1, 2),
					eNN		=LOAD(errors,  0, 2),
					eNNE		=LOAD(errors, -1, 2),
					eNNEE		=LOAD(errors, -2, 2),
					eNNEEE		=LOAD(errors, -3, 2),
					eNNEEEE		=LOAD(errors, -4, 2),
					eNWWWW		=LOAD(errors,  4, 1),
					eNWWW		=LOAD(errors,  3, 1),
					eNWW		=LOAD(errors,  2, 1),
					eNW		=LOAD(errors,  1, 1),
					eN		=LOAD(errors,  0, 1),
					eNE		=LOAD(errors, -1, 1),
					eNEE		=LOAD(errors, -2, 1),
					eNEEE		=LOAD(errors, -3, 1),
					eNEEEE		=LOAD(errors, -4, 1),
					eWWWW		=LOAD(errors,  4, 0),
					eWWW		=LOAD(errors,  3, 0),
					eWW		=LOAD(errors,  2, 0),
					eW		=LOAD(errors,  1, 0);
#endif
				//int
				//	eNNWW	=LOAD(errors,  2, 2),//error = (curr<<8) - pred
				//	eNNW	=LOAD(errors,  1, 2),
				//	eNN	=LOAD(errors,  0, 2),
				//	eNNE	=LOAD(errors, -1, 2),
				//	eNNEE	=LOAD(errors, -2, 2),
				//	eNWW	=LOAD(errors,  2, 1),
				//	eNW	=LOAD(errors,  1, 1),
				//	eN	=LOAD(errors,  0, 1),
				//	eNE	=LOAD(errors, -1, 1),
				//	eNEE	=LOAD(errors, -2, 1),
				//	eWW	=LOAD(errors,  2, 0),
				//	eW	=LOAD(errors,  1, 0);
#undef  LOAD
				//NNNWWWW NNNWWW NNNWW NNNW NNN  NNNE NNNEE NNNEEE NNNEEEE
				//NNWWWW  NNWWW  NNWW  NNW  NN   NNE  NNEE  NNEEE  NNEEEE
				//NWWWW   NWWW   NWW   NW   N    NE   NEE   NEEE   NEEEE
				//WWWW    WWW    WW    W    curr
				int nb[]=
				{
					//slightly better than clampgrad
#if 1
					//NW	N	NE	W	curr
					NNNWWWW,NNNWWW,	NNNWW,	NNWWWW,	//NNWWW
					NNNWWW,	NNNWW,	NNNW,	NNWWW,	//NNWW
					NNNWW,	NNNW,	NNN,	NNWW,	//NNW
					NNNW,	NNN,	NNNE,	NNW,	//NN
					NNN,	NNNE,	NNNEE,	NN,	//NNE
					NNNE,	NNNEE,	NNNEEE,	NNE,	//NNEE
					NNNEE,	NNNEEE,	NNNEEEE,NNEE,	//NNEEE
					NNWWWW,	NNWWW,	NNWW,	NWWWW,	//NWWW
					NNWWW,	NNWW,	NNW,	NWWW,	//NWW
					NNWW,	NNW,	NN,	NWW,	//NW
					NNW,	NN,	NNE,	NW,	//N
					NN,	NNE,	NNEE,	N,	//NE
					NNE,	NNEE,	NNEEE,	NE,	//NEE
					NNEE,	NNEEE,	NNEEEE,	NEE,	//NEEE
					NWWWW,	NWWW,	NWW,	WWWW,	//WWW
					NWWW,	NWW,	NW,	WWW,	//WW
					NWW,	NW,	N,	WW,	//W
#endif

					//worse than clampgrad
#if 0
					//N+eN,			W+eW,			(N+W)/2,		(4*(N+W)+NW-NW)/8,					curr
					NNNWWW	+eNNNWWW,	NNWWWW	+eNNWWWW,	(NNNWWW	+NNWWWW	)>>1,	(4*(NNNWWW	+NNWWWW	)+NNNWW		-NNNWWWW)>>3,//NNWWW
					NNNWW	+eNNNWW,	NNWWW	+eNNWWW,	(NNNWW	+NNWWW	)>>1,	(4*(NNNWW	+NNWWW	)+NNNW		-NNNWWW	)>>3,//NNWW
					NNNW	+eNNNW,		NNWW	+eNNWW,		(NNNW	+NNWW	)>>1,	(4*(NNNW	+NNWW	)+NNN		-NNNWW	)>>3,//NNW
					NNN	+eNNN,		NNW	+eNNW,		(NNN	+NNW	)>>1,	(4*(NNN		+NNW	)+NNNE		-NNNW	)>>3,//NN
					NNNE	+eNNNE,		NN	+eNN,		(NNNE	+NN	)>>1,	(4*(NNNE	+NN	)+NNNEE		-NNN	)>>3,//NNE
					NNNEE	+eNNNEE,	NNE	+eNNE,		(NNNEE	+NNE	)>>1,	(4*(NNNEE	+NNE	)+NNNEEE	-NNNE	)>>3,//NNEE
					NNNEEE	+eNNNEEE,	NNEE	+eNNEE,		(NNNEEE	+NNEE	)>>1,	(4*(NNNEEE	+NNEE	)+NNNEEEE	-NNNEE	)>>3,//NNEEE
					NNWWW	+eNNWWW,	NWWWW	+eNWWWW,	(NNWWW	+NWWWW	)>>1,	(4*(NNWWW	+NWWWW	)+NNWW		-NNWWWW	)>>3,//NWWW
					NNWW	+eNNWW,		NWWW	+eNWWW,		(NNWW	+NWWW	)>>1,	(4*(NNWW	+NWWW	)+NNW		-NNWWW	)>>3,//NWW
					NNW	+eNNW,		NWW	+eNWW,		(NNW	+NWW	)>>1,	(4*(NNW		+NWW	)+NN		-NNWW	)>>3,//NW
					NN	+eNN,		NW	+eNW,		(NN	+NW	)>>1,	(4*(NN		+NW	)+NNE		-NNW	)>>3,//N
					NNE	+eNNE,		N	+eN,		(NNE	+N	)>>1,	(4*(NNE		+N	)+NNEE		-NN	)>>3,//NE
					NNEE	+eNNEE,		NE	+eNE,		(NNEE	+NE	)>>1,	(4*(NNEE	+NE	)+NNEEE		-NNE	)>>3,//NEE
					NNEEE	+eNNEEE,	NEE	+eNEE,		(NNEEE	+NEE	)>>1,	(4*(NNEEE	+NEE	)+NNEEEE	-NNEE	)>>3,//NEEE
					NWWW	+eNWWW,		WWWW	+eWWWW,		(NWWW	+WWWW	)>>1,	(4*(NWWW	+WWWW	)+NWW		-NWWWW	)>>3,//WWW
					NWW	+eNWW,		WWW	+eWWW,		(NWW	+WWW	)>>1,	(4*(NWW		+WWW	)+NW		-NWWW	)>>3,//WW
					NW	+eNW,		WW	+eWW,		(NW	+WW	)>>1,	(4*(NW		+WW	)+N		-NWW	)>>3,//W
#endif

					//X
#if 0
					//NW	N	NE	W		eNW		eN		eNE		eW		curr
					NNNWWWW,NNNWWW,	NNNWW,	NNWWWW,		eNNNWWWW,	eNNNWWW,	eNNNWW,		eNNWWWW,	//NNWWW
					NNNWWW,	NNNWW,	NNNW,	NNWWW,		eNNNWWW,	eNNNWW,		eNNNW,		eNNWWW,		//NNWW
					NNNWW,	NNNW,	NNN,	NNWW,		eNNNWW,		eNNNW,		eNNN,		eNNWW,		//NNW
					NNNW,	NNN,	NNNE,	NNW,		eNNNW,		eNNN,		eNNNE,		eNNW,		//NN
					NNN,	NNNE,	NNNEE,	NN,		eNNN,		eNNNE,		eNNNEE,		eNN,		//NNE
					NNNE,	NNNEE,	NNNEEE,	NNE,		eNNNE,		eNNNEE,		eNNNEEE,	eNNE,		//NNEE
					NNNEE,	NNNEEE,	NNNEEEE,NNEE,		eNNNEE,		eNNNEEE,	eNNNEEEE,	eNNEE,		//NNEEE
					NNWWWW,	NNWWW,	NNWW,	NWWWW,		eNNWWWW,	eNNWWW,		eNNWW,		eNWWWW,		//NWWW
					NNWWW,	NNWW,	NNW,	NWWW,		eNNWWW,		eNNWW,		eNNW,		eNWWW,		//NWW
					NNWW,	NNW,	NN,	NWW,		eNNWW,		eNNW,		eNN,		eNWW,		//NW
					NNW,	NN,	NNE,	NW,		eNNW,		eNN,		eNNE,		eNW,		//N
					NN,	NNE,	NNEE,	N,		eNN,		eNNE,		eNNEE,		eN,		//NE
					NNE,	NNEE,	NNEEE,	NE,		eNNE,		eNNEE,		eNNEEE,		eNE,		//NEE
					NNEE,	NNEEE,	NNEEEE,	NEE,		eNNEE,		eNNEEE,		eNNEEEE,	eNEE,		//NEEE
					NWWWW,	NWWW,	NWW,	WWWW,		eNWWWW,		eNWWW,		eNWW,		eWWWW,		//WWW
					NWWW,	NWW,	NW,	WWW,		eNWWW,		eNWW,		eNW,		eWWW,		//WW
					NWW,	NW,	N,	WW,		eNWW,		eNW,		eN,		eWW,		//W
#endif
				};
				const int priority[]=
				{
					1, 1, 1, 1, 1, 1, 1,
					1, 1, 2, 4, 2, 1, 1,
					1, 1, 4,
		
					//1, 1, 1, 1, 1, 1, 1,
					//1, 1, 1, 3, 1, 1, 1,
					//1, 1, 3,

					//1, 2, 4,  8, 4, 2, 1,
					//2, 4, 8, 16, 8, 4, 2,
					//4, 8, 16,
				};
				//for(int k=0;k<OLS_NSAMPLES;++k)//TODO store preds instead of re-calculating
				//{
				//	int xNW=nb[k<<2|0], xN=nb[k<<2|1], xNE=nb[k<<2|2], xW=nb[k<<2|3];
				//	int grad=xN+xW-xNW;
				//	grad=MEDIAN3(xN, xW, grad);
				//	nb[k<<2|0]=xN;
				//	nb[k<<2|1]=xW;
				//	nb[k<<2|2]=grad;
				//	nb[k<<2|3]=(4*(xN+xW)+xNE-xNW)>>3;
				//}
				double nb2[_countof(nb)];
				double sqm[OLS_NPARAMS*(OLS_NPARAMS<<1)]={0};
				double temp[OLS_NPARAMS<<1];

				for(int ky2=0;ky2<4;++ky2)
				{
					for(int kx2=0;kx2<OLS_NSAMPLES;++kx2)
						nb2[OLS_NSAMPLES*ky2+kx2]=(double)nb[kx2<<2|ky2]/256;///(1<<24);
				}
				for(int ky2=0;ky2<OLS_NPARAMS;++ky2)
				{
					for(int kx2=0;kx2<OLS_NPARAMS;++kx2)
					{
						double sum=0;
						for(int j=0;j<OLS_NSAMPLES;++j)
							sum+=priority[j]*nb2[OLS_NSAMPLES*ky2+j]*nb2[OLS_NSAMPLES*kx2+j];
						sqm[(OLS_NPARAMS<<1)*ky2+kx2]=sum;
					}
					sqm[((OLS_NPARAMS<<1)+1)*ky2+OLS_NPARAMS]=1;
				}
				success=invert_matrix(sqm, OLS_NPARAMS, temp);
				if(success||c1_init[kc])
				{
					//pred = nb * params,		params += ((inv(NB * NBT) * NB * Y) - params)*lr
					if(success)
					{
						int targets[]=
						{
							NNWWW,
							NNWW,
							NNW,
							NN,
							NNE,
							NNEE,
							NNEEE,
							NWWW,
							NWW,
							NW,
							N,
							NE,
							NEE,
							NEEE,
							WWW,
							WW,
							W,
						};
						for(int ky2=0;ky2<OLS_NPARAMS;++ky2)
						{
							double sum=0;
							for(int kx2=0;kx2<OLS_NSAMPLES;++kx2)
								sum+=priority[kx2]*nb2[OLS_NSAMPLES*ky2+kx2]*targets[kx2]/256;///(1<<24);
							temp[ky2]=sum;
						}
						for(int ky2=0;ky2<OLS_NPARAMS;++ky2)
						{
							double sum=0;
							for(int j=0;j<OLS_NPARAMS;++j)
								sum+=sqm[(OLS_NPARAMS<<1)*ky2+j+OLS_NPARAMS]*temp[j];
							if(c1_init[kc])
								c1_params[kc<<2|ky2]+=(sum-c1_params[kc<<2|ky2])*0.2;
							else
								c1_params[kc<<2|ky2]=sum;
						}
						c1_init[kc]=1;
					}
					{
						int nb3[]=
						{
							NW,	N,	NE,	W,//slightly better than clampgrad

							//N+eN,	W+eW,	(N+W)>>1,	(4*(N+W)+NW-NW)>>3,//worse than clampgrad

							//NW,	N,	NE,	W,	eNW,	eN,	eNE,	eW,//X
						};
						//{
						//	int xNW=nb3[0], xN=nb3[1], xNE=nb3[2], xW=nb3[3];
						//	int grad=xN+xW-xNW;
						//	grad=MEDIAN3(xN, xW, grad);
						//	nb3[0]=xN;
						//	nb3[1]=xW;
						//	nb3[2]=grad;
						//	nb3[3]=(4*(xN+xW)+xNE-xNW)>>3;
						//}
						double fpred=0;
						for(int k=0;k<OLS_NPARAMS;++k)
							fpred+=nb3[k]*c1_params[kc<<2|k];
#ifdef _MSC_VER
						pred=_cvt_dtoi_fast(fpred);
#else
						pred=(int)round(fpred);
#endif
					}
					++successcount;
				}
				else//singular matrix
				{
					pred=N+W-NW;
					pred=MEDIAN3(N, W, pred);
				}
				{
					int curr=pred;
					curr^=-fwd;
					curr+=fwd;
					curr+=src->data[idx<<2|kc];
					if(enable_ma)
					{
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
					}
					dst->data[idx<<2|kc]=curr;
				}
			}
		}
	}
	memcpy(src->data, dst->data, (size_t)src->iw*src->ih*sizeof(int[4]));
	if(loud_transforms)
		set_window_title("OLS %lf%%", 100.*successcount/(src->nch*src->iw*src->ih));
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}
	free(dst);
}
#undef  OLS_NPARAMS
#undef  OLS_NSAMPLES

#define OLS2_REACH 2
#define OLS2_STEP 2	//power-of-two
#define OLS2_SAMPLEREACH 7
#define OLS2_NPARAMS (2*(OLS2_REACH+1)*OLS2_REACH+1)
#define OLS2_NSAMPLES ((2*OLS2_SAMPLEREACH+OLS2_STEP+1)*OLS2_SAMPLEREACH)
static int ols2_fallbackpred(const int *pixels, int iw, int kc, int kx, int ky, int idx)
{
	int
		N =ky?pixels[(idx-iw)<<2|kc]:0,
		W =kx?pixels[(idx-1)<<2|kc]:0,
		NW=ky&&kx?pixels[(idx-iw-1)<<2|kc]:0;
	int pred=N+W-NW;
	pred=MEDIAN3(N, W, pred);
	return pred;
}
static void ols2_add_sample(const int *pixels, int depth, int iw, int kc, int kx, int ky, double *sample, double *matrix, int remove)
{
	if(remove)
	{
		for(int ky3=0;ky3<OLS2_NPARAMS;++ky3)
		{
			double vy=sample[ky3];
			for(int kx3=0;kx3<OLS2_NPARAMS;++kx3)
			{
				double vx=sample[kx3];
				matrix[(OLS2_NPARAMS<<1)*ky3+kx3]-=vy*vx;
			}
		}
	}
	{
		int half=(1<<depth)>>1;
		double gain=1./half;
		sample[0]=1;
		for(int ky3=-OLS2_REACH, idx2=1;ky3<=0;++ky3)
		{
			for(int kx3=-OLS2_REACH;kx3<=OLS2_REACH;++kx3, ++idx2)
			{
				//if(idx2>=OLS2_NPARAMS+1)
				//	LOG_ERROR("");
				sample[idx2]=(double)pixels[(iw*(ky+ky3)+kx+kx3)<<2|kc]*gain;
				if(!ky3&&!kx3)//last element is target
					break;
			}
		}
	}
	if(!remove)
	{
		for(int ky3=0;ky3<OLS2_NPARAMS;++ky3)
		{
			double vy=sample[ky3];
			for(int kx3=0;kx3<OLS2_NPARAMS;++kx3)
			{
				double vx=sample[kx3];
				matrix[(OLS2_NPARAMS<<1)*ky3+kx3]+=vy*vx;
			}
		}
	}
}
void pred_ols2(Image *src, int fwd, int enable_ma)
{
	int successcount;
	Image *dst;
	double *samples, *matrix, *matrix2, *temp;
	const int *pixels;
//	const int *errors;

	successcount=0;
	dst=0;
	image_copy(&dst, src);
	samples=(double*)malloc(sizeof(double[OLS2_NSAMPLES*(OLS2_NPARAMS+1)]));
	matrix=(double*)malloc(sizeof(double[OLS2_NPARAMS*OLS2_NPARAMS<<1]));
	matrix2=(double*)malloc(sizeof(double[OLS2_NPARAMS*OLS2_NPARAMS<<1]));
	temp=(double*)malloc(sizeof(double[OLS2_NPARAMS<<1]));
	if(!dst||!samples||!matrix||!matrix2||!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(samples, 0, sizeof(double[OLS2_NSAMPLES*(OLS2_NPARAMS+1)]));
	pixels=fwd?src->data:dst->data;
//	errors=fwd?dst->data:src->data;
	for(int kc=0;kc<4;++kc)
	{
		int initialized,
			nlevels, half;
		double params[OLS2_NPARAMS]={0};

		if(!src->depth[kc])
			continue;
		initialized=0;
		nlevels=1<<src->depth[kc], half=nlevels>>1;
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				ALIGN(16) int pred[4];

				if(kx<OLS2_SAMPLEREACH+OLS2_REACH||kx>src->iw-OLS2_SAMPLEREACH-OLS2_REACH-OLS2_STEP||ky<OLS2_SAMPLEREACH+OLS2_REACH)
					pred[0]=ols2_fallbackpred(pixels, src->iw, kc, kx, ky, idx);
				else
				{
					int success=0;
					if(kx==OLS2_SAMPLEREACH+OLS2_REACH)//initialize row
					{
						memset(matrix, 0, sizeof(double[OLS2_NPARAMS*OLS2_NPARAMS<<1]));
						for(int k=0;k<OLS2_NPARAMS;++k)
						{
							matrix[(OLS2_NPARAMS<<1|1)*k]=0.00005;//add small lambda for regularization
							matrix[(OLS2_NPARAMS<<1|1)*k+OLS2_NPARAMS]=1;//identity
						}
						for(int ky2=-OLS2_SAMPLEREACH, idx2=0;ky2<=0;++ky2)
						{
							for(int kx2=-OLS2_SAMPLEREACH;kx2<=OLS2_SAMPLEREACH+OLS2_STEP-1&&(ky2||kx2);++kx2, ++idx2)
							{
								//if(idx2>=OLS2_NSAMPLES)
								//	LOG_ERROR("");
								ols2_add_sample(pixels, src->depth[kc], src->iw, kc, kx+kx2, ky+ky2, samples+(OLS2_NPARAMS+1)*idx2, matrix, 0);
							}
						}

						memcpy(matrix2, matrix, sizeof(double[OLS2_NPARAMS*OLS2_NPARAMS<<1]));
						success=invert_matrix(matrix2, OLS2_NPARAMS, temp);
					}
					else if(!initialized||!(kx%OLS2_STEP))
					{
						int idx2;
						for(int ky2=0;ky2<OLS2_SAMPLEREACH-1;++ky2)
						{
							idx2=(2*OLS2_SAMPLEREACH+OLS2_STEP)*ky2+(kx-OLS2_SAMPLEREACH-1)%(2*OLS2_SAMPLEREACH+OLS2_STEP);
							ols2_add_sample(pixels, src->depth[kc], src->iw, kc, kx-OLS2_SAMPLEREACH-1, ky-OLS2_SAMPLEREACH+ky2, samples+(OLS2_NPARAMS+1)*idx2, matrix, 1);
							ols2_add_sample(pixels, src->depth[kc], src->iw, kc, kx+OLS2_SAMPLEREACH-1, ky-OLS2_SAMPLEREACH+ky2, samples+(OLS2_NPARAMS+1)*idx2, matrix, 0);
						}
						idx2=(2*OLS2_SAMPLEREACH+OLS2_STEP)*(OLS2_SAMPLEREACH-1)+(kx-OLS2_SAMPLEREACH-1)%OLS2_SAMPLEREACH;
						ols2_add_sample(pixels, src->depth[kc], src->iw, kc, kx-OLS2_SAMPLEREACH, ky, samples+(OLS2_NPARAMS+1)*idx2, matrix, 1);
						ols2_add_sample(pixels, src->depth[kc], src->iw, kc, kx-1, ky, samples+(OLS2_NPARAMS+1)*idx2, matrix, 0);

						memcpy(matrix2, matrix, sizeof(double[OLS2_NPARAMS*OLS2_NPARAMS<<1]));
						success=invert_matrix(matrix2, OLS2_NPARAMS, temp);
					}
					//NNWW NNW NN NNE NNEE
					//NWW  NW  N  NE  NEE
					//WW   W   ?
					if(!success&&!initialized)
						pred[0]=ols2_fallbackpred(pixels, src->iw, kc, kx, ky, idx);
					else
					{
						double fpred;
						int
							N =pixels[(idx-src->iw)<<2|kc],
							W =pixels[(idx-1)<<2|kc],
							NE=pixels[(idx-src->iw+1)<<2|kc];

						if(success)
						{
							//params = ((inv(NB * NBT) * NB * Targets) - params)*lr
							for(int ky2=0;ky2<OLS2_NPARAMS;++ky2)//temp = NB * Targets
							{
								double sum=0;
								for(int kx2=0;kx2<OLS2_NSAMPLES;++kx2)
									sum+=samples[(OLS2_NPARAMS+1)*kx2+ky2]*samples[(OLS2_NPARAMS+1)*kx2+OLS2_NPARAMS];
								temp[ky2]=sum;
							}
							for(int ky2=0;ky2<OLS2_NPARAMS;++ky2)//params = matrix * temp
							{
								double sum=0;
								for(int j=0;j<OLS2_NPARAMS;++j)
									sum+=matrix2[(OLS2_NPARAMS<<1)*ky2+j+OLS2_NPARAMS]*temp[j];
								if(initialized)
									params[ky2]+=(sum-params[ky2])*0.2;
								else
									params[ky2]=sum;
							}
							initialized=1;
						}
						//pred = nb * params
						fpred=params[0]*half;
						for(int ky3=-OLS2_REACH, idx2=1;ky3<=0;++ky3)
						{
							for(int kx3=-OLS2_REACH;kx3<=OLS2_REACH;++kx3, ++idx2)
							{
								if(!ky3&&!kx3)//exclude current pixel
									break;
								fpred+=pixels[(src->iw*(ky+ky3)+kx+kx3)<<2|kc]*params[idx2];
							}
						}
						{
							__m128i mpred=_mm_cvtpd_epi32(_mm_set_sd(fpred));
							__m128i mN	=_mm_set_epi32(0, 0, 0, N);
							__m128i mW	=_mm_set_epi32(0, 0, 0, W);
							__m128i mNE	=_mm_set_epi32(0, 0, 0, NE);
							__m128i vmin=_mm_min_epi32(mN, mW);
							__m128i vmax=_mm_max_epi32(mN, mW);
							vmin=_mm_min_epi32(vmin, mNE);
							vmax=_mm_max_epi32(vmin, mNE);
							mpred=_mm_max_epi32(mpred, vmin);
							mpred=_mm_min_epi32(mpred, vmax);
							_mm_store_si128((__m128i*)pred, mpred);
						}
						//int cmin=W, cmax=W;
						//UPDATE_MIN(cmin, N);
						//UPDATE_MAX(cmax, N);
						//UPDATE_MIN(cmin, NE);
						//UPDATE_MAX(cmax, NE);
						//fpred=CLAMP(cmin, fpred, cmax);
						//pred=(int)round(fpred);
						++successcount;
					}
				}
				{
					int curr=pred[0];
					curr^=-fwd;
					curr+=fwd;
					curr+=src->data[idx<<2|kc];
					if(enable_ma)
					{
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
					}
					dst->data[idx<<2|kc]=curr;
				}
			}
		}
	}
	memcpy(src->data, dst->data, (size_t)src->iw*src->ih*sizeof(int[4]));
	if(loud_transforms)
		set_window_title("OLS2 %lf%%", 100.*successcount/(src->nch*src->iw*src->ih));
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}
	free(dst);
	free(samples);
	free(matrix);
	free(matrix2);
	free(temp);
}
#undef  OLS2_REACH
#undef  OLS2_STEP
#undef  OLS2_SAMPLEREACH
#undef  OLS2_NPARAMS
#undef  OLS2_NSAMPLES

#define OLS3_REACH 1
#define OLS3_STEP 1
#define OLS3_SAMPLEREACH 8
#define OLS3_NPARAMS0 (2*(OLS3_REACH+1)*OLS3_REACH*3+0+1)
#define OLS3_NPARAMS1 (2*(OLS3_REACH+1)*OLS3_REACH*3+1+1)
#define OLS3_NPARAMS2 (2*(OLS3_REACH+1)*OLS3_REACH*3+2+1)
//#define OLS3_NPARAMS_TOTAL (2*(OLS3_REACH+1)*OLS3_REACH*9+3+3)
#define OLS3_NSAMPLES ((2*OLS3_SAMPLEREACH+OLS3_STEP+1)*OLS3_SAMPLEREACH)
static int ols3_fallbackpred(const int *pixels, int iw, int kc, int kx, int ky, int idx)
{
	int
		N =ky?pixels[(idx-iw)<<2|kc]:0,
		W =kx?pixels[(idx-1)<<2|kc]:0,
		NW=ky&&kx?pixels[(idx-iw-1)<<2|kc]:0;
	int pred=N+W-NW;
	pred=MEDIAN3(N, W, pred);
	return pred;
}
static void ols3_add_sample(const int *pixels, double gain, int iw, int kc, int kx, int ky, double *sample, double *matrix, int remove)
{
	int nparams=OLS3_NPARAMS0+kc;
	if(remove)
	{
		for(int ky3=0;ky3<nparams;++ky3)
		{
			double vy=sample[ky3];
			for(int kx3=0;kx3<nparams;++kx3)
			{
				double vx=sample[kx3];
				matrix[(nparams<<1)*ky3+kx3]-=vy*vx;
			}
		}
	}
	sample[0]=1;
	for(int ky3=-OLS3_REACH, idx2=1;ky3<=0;++ky3)
	{
		for(int kx3=-OLS3_REACH;kx3<=OLS3_REACH;++kx3)
		{
			for(int kc2=0;kc2<3;++kc2, ++idx2)
			{
				//if(idx2>=nparams+1)
				//	LOG_ERROR("");
				sample[idx2]=(double)pixels[(iw*(ky+ky3)+kx+kx3)<<2|kc2]*gain;
				if(!ky3&&!kx3&&kc2==kc)//last element is target
					goto finish;
			}
		}
	}
finish:
	if(!remove)
	{
		for(int ky3=0;ky3<nparams;++ky3)
		{
			double vy=sample[ky3];
			for(int kx3=0;kx3<nparams;++kx3)
			{
				double vx=sample[kx3];
				matrix[(nparams<<1)*ky3+kx3]+=vy*vx;
			}
		}
	}
}
void pred_ols3(Image *src, int fwd, int enable_ma)
{
	double t_start;
	int successcount;
	Image *dst;
	double *samples, *matrix, *matrix2, *temp;
	const int *pixels;
//	const int *errors;
	double avparams[OLS3_NPARAMS2*3]={0};
	double params[OLS3_NPARAMS2*3]={0};
	int nlevels[]=
	{
		1<<src->depth[0],
		1<<src->depth[1],
		1<<src->depth[2],
	};
	int half[]=
	{
		nlevels[0]>>1,
		nlevels[1]>>1,
		nlevels[2]>>1,
	};
	double gains[]=
	{
		1./half[0],
		1./half[1],
		1./half[2],
		//1,
		//1,
		//1,
	};

	t_start=time_sec();
	if(loud_transforms)
		DisableProcessWindowsGhosting();
	successcount=0;
	dst=0;
	image_copy(&dst, src);
	samples=(double*)malloc(sizeof(double[OLS3_NSAMPLES*(OLS3_NPARAMS2+1)*3]));
	matrix=(double*)malloc(sizeof(double[OLS3_NPARAMS2*OLS3_NPARAMS2*3<<1]));
	matrix2=(double*)malloc(sizeof(double[OLS3_NPARAMS2*OLS3_NPARAMS2*3<<1]));
	temp=(double*)malloc(sizeof(double[OLS3_NPARAMS2<<1]));
	if(!dst||!samples||!matrix||!matrix2||!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(samples, 0, sizeof(double[OLS3_NSAMPLES*(OLS3_NPARAMS2+1)*3]));
	pixels=fwd?dst->data:src->data;
//	errors=fwd?src->data:dst->data;
	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		int initialized[3]={0};
		for(int kx=0;kx<src->iw;++kx, ++idx)
		{
			for(int kc=0;kc<3;++kc)
			{
				double *cparams=params+OLS3_NPARAMS2*kc;
				double *csamples=samples+(OLS3_NPARAMS2+1)*OLS3_NSAMPLES*kc;
				int nparams=OLS3_NPARAMS0+kc;
				ALIGN(16) int pred[4];

				if(kx<OLS3_SAMPLEREACH+OLS3_REACH||kx>src->iw-OLS3_SAMPLEREACH-OLS3_REACH-OLS3_STEP||ky<OLS3_SAMPLEREACH+OLS3_REACH)
					pred[0]=ols3_fallbackpred(pixels, src->iw, kc, kx, ky, idx);
				else
				{
					int success=0;
					double
						*mat1=matrix+(OLS3_NPARAMS2*OLS3_NPARAMS2<<1)*kc,
						*mat2=matrix2+(OLS3_NPARAMS2*OLS3_NPARAMS2<<1)*kc;
					//if(kx==OLS3_SAMPLEREACH+OLS3_REACH)//initialize row
					if(!(kx%OLS3_STEP))
					{
						memset(mat1, 0, sizeof(double[2])*nparams*nparams);
						for(int k=0;k<nparams;++k)
						{
							mat1[(nparams<<1|1)*k]=0.00005;//add small lambda for regularization
							mat1[(nparams<<1|1)*k+nparams]=1;//identity
						}
						for(int ky2=-OLS3_SAMPLEREACH, idx2=0;ky2<=0;++ky2)
						{
							for(int kx2=-OLS3_SAMPLEREACH;kx2<=OLS3_SAMPLEREACH+OLS3_STEP-1&&(ky2||kx2);++kx2, ++idx2)
								ols3_add_sample(pixels, gains[kc], src->iw, kc, kx+kx2, ky+ky2, csamples+(OLS3_NPARAMS2+1)*idx2, mat1, 0);
						}

						memcpy(mat2, mat1, sizeof(double[2])*nparams*nparams);
						success=invert_matrix(mat2, nparams, temp);
					}
#if 0
					else if(!initialized[kc]||!(kx%OLS3_STEP))//INCOMPATIBLE WITH STEP>1
					{
						int idx2;
						for(int ky2=0;ky2<OLS3_SAMPLEREACH-1;++ky2)
						{
							idx2=(2*OLS3_SAMPLEREACH+OLS3_STEP)*ky2+(kx-OLS3_SAMPLEREACH-1)%(2*OLS3_SAMPLEREACH+OLS3_STEP);
							ols3_add_sample(pixels, gains[kc], src->iw, kc, kx-OLS3_SAMPLEREACH-1, ky-OLS3_SAMPLEREACH+ky2, csamples+(OLS3_NPARAMS2+1)*idx2, mat1, 1);
							ols3_add_sample(pixels, gains[kc], src->iw, kc, kx+OLS3_SAMPLEREACH-1, ky-OLS3_SAMPLEREACH+ky2, csamples+(OLS3_NPARAMS2+1)*idx2, mat1, 0);
						}
						idx2=(2*OLS3_SAMPLEREACH+OLS3_STEP)*(OLS3_SAMPLEREACH-1)+(kx-OLS3_SAMPLEREACH-1)%OLS3_SAMPLEREACH;
						ols3_add_sample(pixels, gains[kc], src->iw, kc, kx-OLS3_SAMPLEREACH, ky, csamples+(OLS3_NPARAMS2+1)*idx2, mat1, 1);
						ols3_add_sample(pixels, gains[kc], src->iw, kc, kx-1, ky, csamples+(OLS3_NPARAMS2+1)*idx2, mat1, 0);

						memcpy(mat2, mat1, sizeof(double[2])*nparams*nparams);
						success=invert_matrix(mat2, nparams, temp);
					}
#endif
					//NNWW NNW NN NNE NNEE
					//NWW  NW  N  NE  NEE
					//WW   W   ?
					if(!success&&!initialized[kc])
						pred[0]=ols3_fallbackpred(pixels, src->iw, kc, kx, ky, idx);
					else
					{
						if(success)
						{
							//params = ((inv(NB * NBT) * NB * Targets) - params)*lr
							for(int ky2=0;ky2<nparams;++ky2)//temp = NB * Targets
							{
								double sum=0;
								for(int kx2=0;kx2<OLS3_NSAMPLES;++kx2)
									sum+=csamples[(OLS3_NPARAMS2+1)*kx2+ky2]*csamples[(OLS3_NPARAMS2+1)*kx2+nparams];
								temp[ky2]=sum;
							}
							for(int ky2=0;ky2<nparams;++ky2)//params = matrix * temp
							{
								double sum=0;
								for(int j=0;j<nparams;++j)
									sum+=mat2[(nparams<<1)*ky2+j+nparams]*temp[j];
								if(initialized[kc])
									cparams[ky2]+=(sum-cparams[ky2])*0.2;
								else
									cparams[ky2]=sum;
								avparams[OLS3_NPARAMS2*kc+ky2]+=cparams[ky2];
							}
							initialized[kc]=1;
						}
						{
							//pred = nb * params
							int
								N =pixels[(idx-src->iw)<<2|kc],
								W =pixels[(idx-1)<<2|kc],
								NE=pixels[(idx-src->iw+1)<<2|kc];
							double fpred=cparams[0]/gains[kc];
							for(int ky3=-OLS3_REACH, idx2=1;ky3<=0;++ky3)
							{
								for(int kx3=-OLS3_REACH;kx3<=OLS3_REACH;++kx3)
								{
									for(int kc2=0;kc2<3;++kc2, ++idx2)
									{
										if(!ky3&&!kx3&&kc2==kc)//exclude current pixel
											goto finish;
										fpred+=pixels[(src->iw*(ky+ky3)+kx+kx3)<<2|kc2]*cparams[idx2];
									}
								}
							}
						finish:
							{
								__m128i mpred=_mm_cvtpd_epi32(_mm_set_sd(fpred));
								__m128i mN	=_mm_set_epi32(0, 0, 0, N);
								__m128i mW	=_mm_set_epi32(0, 0, 0, W);
								__m128i mNE	=_mm_set_epi32(0, 0, 0, NE);
								__m128i vmin=_mm_min_epi32(mN, mW);
								__m128i vmax=_mm_max_epi32(mN, mW);
								vmin=_mm_min_epi32(vmin, mNE);
								vmax=_mm_max_epi32(vmax, mNE);
								mpred=_mm_max_epi32(mpred, vmin);
								mpred=_mm_min_epi32(mpred, vmax);
								_mm_store_si128((__m128i*)pred, mpred);
							}
							//int cmin=W, cmax=W;
							//UPDATE_MIN(cmin, N);
							//UPDATE_MAX(cmax, N);
							//UPDATE_MIN(cmin, NE);
							//UPDATE_MAX(cmax, NE);
							//fpred=CLAMP(cmin, fpred, cmax);
							//pred=(int)round(fpred);
							++successcount;
						}
					}
				}
				{
					int curr=pred[0];
					curr^=-fwd;
					curr+=fwd;
					curr+=src->data[idx<<2|kc];
					if(enable_ma)
					{
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
					}
					src->data[idx<<2|kc]=curr;
				}
			}
		}
		if(loud_transforms)
			set_window_title("%d/%d = %7.3lf%%  OLS3 rate %lf%%  %lf sec", ky+1, src->ih, 100.*(ky+1)/src->ih, 100.*successcount/(src->nch*src->iw*(ky+1)), time_sec()-t_start);
	}
	//memcpy(src->data, dst->data, (size_t)src->iw*src->ih*sizeof(int[4]));
	if(loud_transforms)
	{
		const char *nbnames[]=
		{
			"bias",
			"NW",
			"N",
			"NE",
			"W",
			"curr",
		};
		set_window_title("OLS3 %lf%%  %lf sec", 100.*successcount/(src->nch*src->iw*src->ih), time_sec()-t_start);
		for(int k=0;k<(int)_countof(avparams);++k)
			avparams[k]*=(double)src->nch/successcount;
		int printed=0;
		printed+=snprintf(g_buf+printed, G_BUF_SIZE-1-printed, "Y U V\n");
		for(int kg=0, idx=0;kg<3;++kg)
		{
			int idx2=0;
			for(int k=0;k<OLS3_NPARAMS2;++k)
			{
				printed+=snprintf(g_buf+printed, G_BUF_SIZE-1-printed, " %+12.8lf", avparams[idx++]);
				if(!(k%3))
					printed+=snprintf(g_buf+printed, G_BUF_SIZE-1-printed, "  %s\n", nbnames[idx2++]);
			}
			printed+=snprintf(g_buf+printed, G_BUF_SIZE-1-printed, "  %s", nbnames[idx2]);
			//for(int kx=0;kx<5;++kx)
			//{
			//	for(int kc=0;kc<3;++kc)
			//		printed+=snprintf(g_buf+printed, G_BUF_SIZE-1-printed, " %+12.8lf", avparams[(kg*3+kx)*3+kc]);
			//	printed+=snprintf(g_buf+printed, G_BUF_SIZE-1-printed, "\n");
			//}
			if(kg<2)
				printed+=snprintf(g_buf+printed, G_BUF_SIZE-1-printed, "\n\n");
		}
		copy_to_clipboard(g_buf, printed);
		messagebox(MBOX_OK, "Copied to clipboard", g_buf);
	}
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}
	free(dst);
	free(samples);
	free(matrix);
	free(matrix2);
	free(temp);
}
#undef  OLS3_REACH
#undef  OLS3_STEP
#undef  OLS3_SAMPLEREACH
#undef  OLS3_NPARAMS
#undef  OLS3_NSAMPLES

void pred_select(Image const *src, Image *dst, int fwd, int enable_ma)
{
	const int *pixels=fwd?src->data:dst->data;
	for(int kc=0;kc<4;++kc)
	{
		int nlevels;

		if(!src->depth[kc])
			continue;
		nlevels=1<<src->depth[kc];
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				int
					N=ky?pixels[(idx-src->iw)<<2|kc]:0,
					W=kx?pixels[(idx-1)<<2|kc]:0,
					NW=kx&&ky?pixels[(idx-src->iw-1)<<2|kc]:0;
				int pred=abs(W-NW)<abs(N-NW)?N:W;
				
				pred^=-fwd;
				pred+=fwd;
				pred+=src->data[idx<<2|kc];
				if(enable_ma)
				{
					pred+=nlevels>>1;
					pred&=nlevels-1;
					pred-=nlevels>>1;
				}
				dst->data[idx<<2|kc]=pred;
			}
		}
	}
	dst->depth[0]=src->depth[0]+!enable_ma;
	dst->depth[1]=src->depth[1]+!enable_ma;
	dst->depth[2]=src->depth[2]+!enable_ma;
	dst->depth[3]=src->depth[3]+(!enable_ma&(src->depth[3]!=0));
}

#define LIN_NNB (2*(LIN_REACH+1)*LIN_REACH)
void pred_linear(Image const *src, Image *dst, const int *coeffs, int lgden, int fwd, int enable_ma)
{
	const int *pixels=fwd?src->data:dst->data;
	for(int kc=0;kc<4;++kc)
	{
		if(!src->depth[kc])
			continue;
		for(int ky=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx)
			{
				int nb[LIN_NNB], idx;
				int pred;

				idx=0;
				for(int ky2=-LIN_REACH;ky2<=0;++ky2)
				{
					for(int kx2=-LIN_REACH;kx2<=LIN_REACH;++kx2, ++idx)
					{
						if(!ky2&&!kx2)
							break;
						if((unsigned)(ky+ky2)<(unsigned)src->ih&&(unsigned)(kx+kx2)<(unsigned)src->iw)
							nb[idx]=pixels[(src->iw*(ky+ky2)+kx+kx2)<<2|kc];
						else
							nb[idx]=0;
					}
				}
				pred=(1<<lgden)>>1;
				for(int k=0;k<LIN_NNB;++k)
					pred+=coeffs[k]*nb[k];
				pred>>=lgden;
				idx=(src->iw*ky+kx)<<2|kc;
				
				//if(kx==5&&ky==5)//
				//	printf("");

				pred^=-fwd;
				pred+=fwd;
				pred+=src->data[idx];
				if(enable_ma)
				{
					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];
				}
				dst->data[idx]=pred;
			}
		}
	}
	dst->depth[0]=src->depth[0]+!enable_ma;
	dst->depth[1]=src->depth[1]+!enable_ma;
	dst->depth[2]=src->depth[2]+!enable_ma;
	dst->depth[3]=src->depth[3]+(!enable_ma&(src->depth[3]!=0));
}


#define SLIC5_NPREDS 32
#define SLIC5_PREDLIST\
	SLIC5_PRED(W+NE-N-((2*(eN+eW)+eNE-eNW+4)>>3))\
	SLIC5_PRED(N-(int)(((long long)eN+eW+eNE)*-0x05C>>8))\
	SLIC5_PRED(W-(int)(((long long)eN+eW+eNW)*-0x05B>>8))\
	SLIC5_PRED(N+(int)((-eNN*0x0DFLL-eN*0x051LL-eNE*0x0BDLL+((long long)N-NN)*0x05C+((long long)NW-W)*0x102)>>8))\
	SLIC5_PRED(3*(N-NN)+NNN)\
	SLIC5_PRED((N+W)>>1)\
	SLIC5_PRED(N+W-NW)\
	SLIC5_PRED((W+NEE)>>1)\
	SLIC5_PRED((3*W+NEEE)>>2)\
	SLIC5_PRED((3*(3*W+NE+NEE)-10*N+2)/5)\
	SLIC5_PRED((3*(3*W+NE+NEE)-10*N)/5)\
	SLIC5_PRED((4*N-2*NN+NW+NE)>>2)\
	SLIC5_PRED(N+NE-NNE-eNNE)\
	SLIC5_PRED((4*(N+W+NW+NE)-(NN+WW+NNWW+NNEE)+6)/12)\
	SLIC5_PRED(W+((eW-eWW)>>1))\
	SLIC5_PRED(paper_GAP)\
	SLIC5_PRED(calic_GAP)\
	SLIC5_PRED(N+W-((NW+NN+WW+NE)>>2))\
	SLIC5_PRED(((2*(N+W)-(NW+NN+WW+NE))*9+(WWW+NWW+NNW+NNN+NNE+NEE)*2)/12)\
	SLIC5_PRED(3*(N+W-NW-(NN+WW-NNWW))+NNN+WWW-NNNWWW)\
	SLIC5_PRED(2*(W+NE-N)-(WW+NNEE-NN))\
	SLIC5_PRED((2*W+NEE-N)>>1)\
	SLIC5_PRED(NW+NWW-NNWWW)\
	SLIC5_PRED((14*NE-(NNEE+NNNEE+NNEEE))/11)\
	SLIC5_PRED((NEEE+NEEEE)>>1)\
	SLIC5_PRED((NNNEEEE+NNEEE)>>1)\
	SLIC5_PRED(NNEEEE)\
	SLIC5_PRED((NNWWWW+NNNWWWW)>>1)\
	SLIC5_PRED((WWW+WWWW)>>1)\
	SLIC5_PRED((N+NN)>>1)\
	SLIC5_PRED((NE+NNEE)>>1)\
	SLIC5_PRED((NE+NNE+NEE+NNEE)>>2)
//	SLIC5_PRED(geomean)
void pred_t47(Image *src, int fwd, int enable_ma)
{
	Image *dst;
	int *pred_errors;
	const int *pixels, *errors;

	dst=0;
	image_copy(&dst, src);
	pred_errors=(int*)malloc((src->iw+4LL)*sizeof(int[SLIC5_NPREDS*4]));//4 padded rows
	if(!dst||!pred_errors)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	pixels=fwd?dst->data:src->data;
	errors=fwd?src->data:dst->data;
	for(int kc=0;kc<4;++kc)
	{
		int sh;
		int params[SLIC5_NPREDS], param0=176;
		long long bias_sum=0;
		int bias_count=0;

		if(!src->depth[kc])
			continue;
		sh=src->depth[kc];
		memfill(params, &param0, sizeof(params), sizeof(int));
		memset(pred_errors, 0, (src->iw+4LL)*sizeof(int[SLIC5_NPREDS*4]));
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			int kym[]=
			{
				(src->iw+4)*(ky&3),
				(src->iw+4)*((ky-1)&3),
				(src->iw+4)*((ky-2)&3),
				(src->iw+4)*((ky-3)&3),
			};
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
#define LOAD(BUF, X, Y) ((unsigned)(ky+(Y))<(unsigned)src->ih&&(unsigned)(kx+(X))<(unsigned)src->iw?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]<<8:0)
				int
				//	NNNNWW	=LOAD(pixels, -2, -4),
				//	NNNNW	=LOAD(pixels, -1, -4),
				//	NNNN	=LOAD(pixels,  0, -4),
				//	NNNNE	=LOAD(pixels,  1, -4),
				//	NNNNEE	=LOAD(pixels,  2, -4),
				//	NNNNEEEE=LOAD(pixels,  4, -4),

					NNNWWWW	=LOAD(pixels, -4, -3),
					NNNWWW	=LOAD(pixels, -3, -3),
				//	NNNWW	=LOAD(pixels, -2, -3),
				//	NNNW	=LOAD(pixels, -1, -3),
					NNN	=LOAD(pixels,  0, -3),
				//	NNNE	=LOAD(pixels,  1, -3),
					NNNEE	=LOAD(pixels,  2, -3),
				//	NNNEEE	=LOAD(pixels,  3, -3),
					NNNEEEE	=LOAD(pixels,  4, -3),
		
					NNWWWW	=LOAD(pixels, -4, -2),
					NNWWW	=LOAD(pixels, -3, -2),
					NNWW	=LOAD(pixels, -2, -2),
					NNW	=LOAD(pixels, -1, -2),
					NN	=LOAD(pixels,  0, -2),
					NNE	=LOAD(pixels,  1, -2),
					NNEE	=LOAD(pixels,  2, -2),
					NNEEE	=LOAD(pixels,  3, -2),
					NNEEEE	=LOAD(pixels,  4, -2),
				//	NWWWW	=LOAD(pixels, -4, -1),
				//	NWWW	=LOAD(pixels, -3, -1),
					NWW	=LOAD(pixels, -2, -1),
					NW	=LOAD(pixels, -1, -1),
					N	=LOAD(pixels,  0, -1),
					NE	=LOAD(pixels,  1, -1),
					NEE	=LOAD(pixels,  2, -1),
					NEEE	=LOAD(pixels,  3, -1),
					NEEEE	=LOAD(pixels,  4, -1),
					WWWW	=LOAD(pixels, -4,  0),
					WWW	=LOAD(pixels, -3,  0),
					WW	=LOAD(pixels, -2,  0),
					W	=LOAD(pixels, -1,  0),

				//	eNNWW	=LOAD(errors, -2, -2),//error = (curr<<8) - pred
				//	eNNW	=LOAD(errors, -1, -2),
					eNN	=LOAD(errors,  0, -2),
					eNNE	=LOAD(errors,  1, -2),
				//	eNNEE	=LOAD(errors,  2, -2),
				//	eNWW	=LOAD(errors, -2, -1),
					eNW	=LOAD(errors, -1, -1),
					eN	=LOAD(errors,  0, -1),
					eNE	=LOAD(errors,  1, -1),
				//	eNEE	=LOAD(errors,  2, -1),
					eWW	=LOAD(errors, -2,  0),
					eW	=LOAD(errors, -1,  0);
#undef  LOAD
				int paper_GAP, calic_GAP;
				int geomean;
				{
					int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
					int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
					int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
					int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
					int diff=(dy-dx)>>sh, diff2=(d45-d135)>>sh, diff3=NE-NW;
					if(dy+dx>32)
						paper_GAP=(int)(((long long)dx*N+(long long)dy*W)/((long long)dy+dx));
					else if(diff>12)
						paper_GAP=(N+2*W)/3;
					else if(diff<-12)
						paper_GAP=(2*N+W)/3;
					else
						paper_GAP=(N+W)>>1;

					if(diff2>32)
						paper_GAP+=diff3>>2;
					else if(diff2>16)
						paper_GAP+=diff3*3>>4;
					else if(diff2>=-16)
						paper_GAP+=diff3>>3;
					else if(diff2>=-32)
						paper_GAP+=diff3>>4;

					if(diff>80)
						calic_GAP=W;
					else if(diff<-80)
						calic_GAP=N;
					else if(diff>32)
						calic_GAP=(2*N+6*W+NE-NW)>>3;		//c1	[1/4  3/4  1/8  -1/8].[N W NE NW]
					else if(diff>8)
						calic_GAP=(6*N+10*W+3*(NE-NW))>>4;	//c2	[3/8  5/8  3/16  -3/16]
					else if(diff<-32)
						calic_GAP=(6*N+2*W+NE-NW)>>3;		//c3	[3/4  1/4  1/8  -1/8]
					else if(diff<-8)
						calic_GAP=(10*N+6*W+3*(NE-NW))>>4;	//c4	[5/8  3/8  3/16  -3/16]
					else
						calic_GAP=(((N+W)<<1)+NE-NW)>>2;	//c5	[1/2  1/2  1/4  -1/4]
				}
				{
					long long aN=N+0x800000LL, aW=W+0x800000LL;
					aN=CLAMP(0, aN, 0xFFFFFF);
					aW=CLAMP(0, aW, 0xFFFFFF);
					geomean=(int)floor_sqrt(aN*aW)-0x800000;
				}
				{
					int val, curr, kbest, besterr;
					int preds[]=
					{
#define SLIC5_PRED(X) X,
						SLIC5_PREDLIST
#undef  SLIC5_PRED
					};
					long long pred=0, wsum=0;
#define LOAD(mX, mY) pred_errors[SLIC5_NPREDS*(kym[mY]+kx+2-(mX))+k]
					for(int k=0;k<(int)_countof(preds);++k)
					{
						long long weight=1+
							(long long)LOAD(-2,  2)+
							(long long)LOAD(-1,  2)+
							(long long)LOAD( 0,  2)+
							(long long)LOAD( 1,  2)+
							(long long)LOAD(-2,  1)+
							(long long)LOAD(-1,  1)+
							(long long)LOAD( 0,  1)*3+
							(long long)LOAD( 1,  1);
						weight=((long long)params[k]<<29)/wndw;
						pred+=weight*preds[k];
						wsum+=weight;
					}
					if(wsum)
					{
						pred+=(wsum>>1)-1;
						pred/=wsum;
					}
					else
						pred=preds[4];
					pred+=bias_sum/(bias_count+1LL);
					CLAMP3_32(pred, (int)pred, N, W, NE);
					//{
					//	int clamp_lo=N, clamp_hi=N;
					//	UPDATE_MIN(clamp_lo, W);
					//	UPDATE_MAX(clamp_hi, W);
					//	UPDATE_MIN(clamp_lo, NE);
					//	UPDATE_MAX(clamp_hi, NE);
					//	pred=CLAMP(clamp_lo, pred, clamp_hi);
					//}
					val=(int)((pred+127)>>8);

					//if(!kc&&ky==0&&kx==1)
					//	printf("");
					val^=-fwd;
					val+=fwd;
					val+=src->data[idx<<2|kc];
					if(enable_ma)
					{
						val<<=32-sh;
						val>>=32-sh;
					}
					src->data[idx<<2|kc]=val;

					curr=pixels[idx<<2|kc]<<8;
					kbest=-1, besterr=0;
					for(int k=0;k<(int)_countof(preds);++k)
					{
						int e=abs(curr-preds[k]);
						LOAD(0, 0)=e;
						if(ky&&kx<src->iw-1)
							LOAD(-1, 1)+=e;
						if(kbest==-1||besterr>e)
							kbest=k, besterr=e;
					}
					++params[kbest];
					if(params[kbest]>352)
					{
						for(int k=0;k<(int)_countof(preds);++k)
							params[k]=(params[k]+1)>>1;
					}
					++bias_count;
					bias_sum+=curr-pred;
				}
#undef  LOAD
			}
		}
	}
	free(dst);
	free(pred_errors);
}

//CUSTOM reach-2 predictor
int custom_params[CUSTOM_NPARAMS]={0};
int custom_clamp[4]={0};
void pred_custom(Image *src, int fwd, int enable_ma, const int *params)
{
	int bufsize=(src->iw+16)*(int)sizeof(int[4*4*2]);//4 padded rows * 4 channels max * {pixel, pred}
	int *pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	//Image *dst;
	//int *pixels, *errors;

	//dst=0;
	//image_copy(&dst, src);
	//pixels=fwd?dst->data:src->data;
	//errors=fwd?src->data:dst->data;
	//int nlevels;

	for(int ky=0, idx=0;ky<src->ih;++ky)
	{
		ALIGN(16) int *rows[]=
		{
			pixels+((src->iw+16LL)*((ky-0LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-1LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-2LL)&3)+8LL)*4*2,
			pixels+((src->iw+16LL)*((ky-3LL)&3)+8LL)*4*2,
		};
		for(int kx=0;kx<src->iw;++kx, ++idx)
		{
			for(int kc=0;kc<3;++kc)//there are params only for 3 channels in CUSTOM predictor
			{
				int nlevels=1<<src->depth[kc];
				if(!src->depth[kc])
					continue;
				int nb[CUSTOM_NNB*2]=
				{
					rows[2][kc-2*8+0], rows[2][kc-2*8+4],
					rows[2][kc-1*8+0], rows[2][kc-1*8+4],
					rows[2][kc+0*8+0], rows[2][kc+0*8+4],
					rows[2][kc+1*8+0], rows[2][kc+1*8+4],
					rows[2][kc+2*8+0], rows[2][kc+2*8+4],
					rows[1][kc-2*8+0], rows[1][kc-2*8+4],
					rows[1][kc-1*8+0], rows[1][kc-1*8+4],
					rows[1][kc+0*8+0], rows[1][kc+0*8+4],
					rows[1][kc+1*8+0], rows[1][kc+1*8+4],
					rows[1][kc+2*8+0], rows[1][kc+2*8+4],
					rows[0][kc-2*8+0], rows[0][kc-2*8+4],
					rows[0][kc-1*8+0], rows[0][kc-1*8+4],
				};
				int clamp_idx[]={22, 12, 14, 16};//W NW N NE
				int vmin, vmax, uninitialized;
				long long pred;
				int curr=src->data[idx<<2|kc];

				//for(int ky2=-CUSTOM_REACH, idx2=0;ky2<=0;++ky2)
				//{
				//	for(int kx2=-CUSTOM_REACH;kx2<=CUSTOM_REACH&&idx2<CUSTOM_NNB*2;++kx2, idx2+=2)
				//	{
				//		if((unsigned)(ky+ky2)<(unsigned)src->ih&&(unsigned)(kx+kx2)<(unsigned)src->iw)
				//		{
				//			nb[idx2  ]=pixels[(src->iw*(ky+ky2)+kx+kx2)<<2|kc];
				//			nb[idx2+1]=errors[(src->iw*(ky+ky2)+kx+kx2)<<2|kc];
				//		}
				//		else
				//		{
				//			nb[idx2  ]=0;
				//			nb[idx2+1]=0;
				//		}
				//	}
				//}
				vmin=-(nlevels>>1);
				vmax=(nlevels>>1);
				uninitialized=1;
				for(int k=0;k<4;++k)
				{
					if(custom_clamp[k])
					{
						int val=nb[clamp_idx[k]];
						if(uninitialized)
						{
							vmin=val;
							vmax=val;
							uninitialized=0;
						}
						else
						{
							UPDATE_MIN(vmin, val);
							UPDATE_MAX(vmax, val);
						}
					}
				}

				pred=0;
				for(int k=0;k<CUSTOM_NNB*2;++k)
					pred+=(long long)nb[k]*params[kc*CUSTOM_NNB*2+k];
				pred+=1LL<<15;
				pred>>=16;
				pred=CLAMP(vmin, pred, vmax);

				rows[0][kc+4*!fwd]=curr;
				pred^=-fwd;
				pred+=fwd;
				pred+=(long long)curr;
				//pred+=(long long)dst->data[(size_t)idx<<2|kc];
				if(enable_ma)
				{
					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];
				}
				src->data[idx<<2|kc]=(int)pred;
				rows[0][kc+4* fwd]=(int)pred;
			}
			rows[0]+=8;
			rows[1]+=8;
			rows[2]+=8;
			rows[3]+=8;
		}
		//params+=24;
	}
	//memcpy(src->data, dst->data, (size_t)src->iw*src->ih*sizeof(int[4]));
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}
	free(pixels);
	//free(dst);
}
static void pred_custom_calcloss(Image const *src, Image *dst, int *hist, const int *params, double *loss)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	memcpy(dst->data, src->data, res*sizeof(int[4]));

	pred_custom(dst, 1, 1, params);

	for(int kc=0;kc<3;++kc)
	{
		int nlevels=1<<src->depth[kc];
		calc_histogram(dst->data, dst->iw, dst->ih, kc, 0, dst->iw, 0, dst->ih, src->depth[kc], hist, 0);
		{
			double entropy=calc_entropy(hist, nlevels, (int)res);
			loss[kc]=entropy/src->src_depth[(kc+1)%3];
		}
	}
	loss[3]=(loss[0]+loss[1]+loss[2])/3;
}
#define CUSTOM_NITER 128
#define CUSTOM_DELTAGROUP 2//must be 2		//must be even
void pred_custom_optimize(Image const *image, int *params)
{
	static int call_idx=0;

	Image *im2;
	int maxdepth, maxlevels, *hist;
	double loss_bestsofar[4], loss_prev[4], loss_curr[4];
	int params2[CUSTOM_NPARAMS];
	int shakethreshold;

	++call_idx;
	if(call_idx==1)
		DisableProcessWindowsGhosting();

	im2=0;
	image_copy(&im2, image);
	maxdepth=calc_maxdepth(image, 0);
	maxlevels=1<<maxdepth;
	hist=(int*)malloc(maxlevels*sizeof(int));
	if(!im2||!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}

	memcpy(params2, params, sizeof(params2));

#define CALC_LOSS(L) pred_custom_calcloss(image, im2, hist, params, L)
#ifndef _DEBUG
	srand((unsigned)__rdtsc());
#endif
	CALC_LOSS(loss_bestsofar);
	//for(int kc=0;kc<3;++kc)
	{
		int sum=0, idx;
		for(int k=0;k<CUSTOM_NNB*2;k+=2)
			sum+=params[CUSTOM_NNB*2*custom_pred_ch_idx+k];
		idx=rand()%CUSTOM_NNB;
		params[CUSTOM_NNB*2*custom_pred_ch_idx+idx*2]+=0x10000-sum;
	}
	CALC_LOSS(loss_prev);

	shakethreshold=CUSTOM_NPARAMS;
	for(int it=0, watchdog=0;it<CUSTOM_NITER;++it)
	{
		int idx[CUSTOM_DELTAGROUP]={0}, stuck=0;
		int params_original_selected[CUSTOM_DELTAGROUP]={0};

		if(watchdog>=shakethreshold)//bump if stuck?
		{
			//memcpy(params2, params, sizeof(params2));
			//for(int k=0;k<CUSTOM_NPARAMS;++k)
			//	params2[k]+=((rand()&1)<<1)-1;
			watchdog=0;
			stuck=1;
		}
		else
		{
			for(int k=0;k<CUSTOM_DELTAGROUP;++k)//increment params
			{
				int inc=0;
#if 1
				//inc=(rand()&0xFFFF)-0x8000;
				//inc=(rand()&0x0FF0)-0x0800;
				inc=rand()&1?0x100:-0x100;
				
			//	idx[k]=(CUSTOM_NNB*2*custom_pred_ch_idx+rand()%(CUSTOM_NNB*2));
				idx[k]=(CUSTOM_NNB*custom_pred_ch_idx+rand()%CUSTOM_NNB)*2;
				++k;
			//	idx[k]=(CUSTOM_NNB*2*custom_pred_ch_idx+rand()%(CUSTOM_NNB*2));
				idx[k]=(CUSTOM_NNB*custom_pred_ch_idx+rand()%CUSTOM_NNB)*2;

				params_original_selected[k-1]=params[idx[k-1]];
				params[idx[k-1]]+=inc;
				params_original_selected[k+0]=params[idx[k+0]];
				params[idx[k+0]]-=inc;
				//{
				//	int sum=0;
				//	for(int k=0;k<CUSTOM_NNB*2;k+=2)
				//		sum+=params[CUSTOM_NNB*2*custom_pred_ch_idx+k];
				//	if(sum!=0x10000)
				//		LOG_ERROR("");
				//}
#endif
#if 0
				while(!(inc=rand()-(RAND_MAX>>1)));//reject zero delta
				idx[k]=rand()%CUSTOM_NPARAMS;
				if(!(idx[k]&1))//maintain pixel coefficient sum at 1
				{
					if(k<CUSTOM_DELTAGROUP-1)
					{
						params_original_selected[k]=params2[idx[k]];
						params2[idx[k]]+=(int)(((long long)inc*CUSTOM_NITER<<8)/((it+1)*RAND_MAX));
						++k;
						idx[k]=rand()%CUSTOM_NPARAMS;
						params_original_selected[k]=params2[idx[k]];
						params2[idx[k]]-=(int)(((long long)inc*CUSTOM_NITER<<8)/((it+1)*RAND_MAX));
					}
					else
					{
						++idx[k];
						params_original_selected[k]=params2[idx[k]];
						params2[idx[k]]+=(int)(((long long)inc*CUSTOM_NITER<<8)/((it+1)*RAND_MAX));
					}
				}
				else
				{
					params_original_selected[k]=params2[idx[k]];
					params2[idx[k]]+=(int)(((long long)inc*CUSTOM_NITER<<8)/((it+1)*RAND_MAX));
					//params2[idx[k]]+=inc*(16<<1)/RAND_MAX;
				}
#endif
			}
		}

		CALC_LOSS(loss_curr);
		
		if(loss_prev[3]<loss_curr[3])//revert if worse
		{
			if(stuck)//a bad branch may surpass the local minimum
				memcpy(loss_prev, loss_curr, sizeof(loss_prev));
			else
			{
				memcpy(loss_curr, loss_prev, sizeof(loss_prev));
				for(int k=0;k<CUSTOM_DELTAGROUP;++k)
					params[idx[k]]=params_original_selected[k];
			}
			++watchdog;
		}
		else//save if better
		{
			if(loss_curr[3]<loss_bestsofar[3])//publish if record best
			{
				memcpy(params2, params, sizeof(params2));
				memcpy(loss_bestsofar, loss_curr, sizeof(loss_bestsofar));
				--it;//again
			}
			memcpy(loss_prev, loss_curr, sizeof(loss_prev));
			watchdog=0;
		}
		if(loud_transforms)
			set_window_title(
				"%d %4d/%4d,%d/%d: %lf%% RGB %lf %lf %lf%s",
				call_idx,
				it+1,
				CUSTOM_NITER,
				watchdog,
				shakethreshold,
				100.*loss_bestsofar[3],
				100.*loss_bestsofar[0],
				100.*loss_bestsofar[1],
				100.*loss_bestsofar[2],
				it+1<CUSTOM_NITER?"...":" Done."
			);

		//preview
#if 1
		{
			ch_entropy[0]=(float)(loss_bestsofar[0]*image->src_depth[0]);
			ch_entropy[1]=(float)(loss_bestsofar[1]*image->src_depth[1]);
			ch_entropy[2]=(float)(loss_bestsofar[2]*image->src_depth[2]);
			ch_entropy[3]=0;
			//ch_entropy[3]=(float)(loss_bestsofar[3]*image->src_depth[3]);//X
			//unsigned char *ptr;
			//addhalf(temp, iw, ih, 3, 4);
			//SWAPVAR(image, temp, ptr);
			io_render();
			//SWAPVAR(image, temp, ptr);
		}
#endif
	}
	memcpy(params, params2, sizeof(params2));
#undef  CALC_LOSS
	free(hist);
	free(im2);
}

void conv_custom(Image *src)
{
	Image *dst=0;
	image_copy(&dst, src);
	if(!dst)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(int ky=0;ky<src->ih;++ky)
	{
		for(int kx=0;kx<src->iw;++kx)
		{
			for(int kc=0;kc<3;++kc)
			{
				int nb[]=
				{
#define LOAD(BUF, X, Y) ((unsigned)(kx+(X))<(unsigned)src->iw&&(unsigned)(ky+(Y))<(unsigned)src->ih?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
					LOAD(src->data, -2, -2),
					LOAD(src->data, -1, -2),
					LOAD(src->data,  0, -2),
					LOAD(src->data,  1, -2),
					LOAD(src->data,  2, -2),
					LOAD(src->data, -2, -1),
					LOAD(src->data, -1, -1),
					LOAD(src->data,  0, -1),
					LOAD(src->data,  1, -1),
					LOAD(src->data,  2, -1),
					LOAD(src->data, -2,  0),
					LOAD(src->data, -1,  0),
#undef  LOAD
				};
				int *curr_params=custom_params+kc*CUSTOM_NNB*2;
				if(!curr_params[0])
					continue;
				int idx=(src->iw*ky+kx)<<2|kc;
				int val=dst->data[idx]*curr_params[0];
				for(int k=1;k<_countof(nb);++k)
					val+=nb[k]*curr_params[k<<1];
				val+=1<<15;
				val>>=16;
				src->data[idx]=val;
			}
		}
	}
	free(dst);
}




double pw2_errors[PW2_NPRED]={0};//
short pw2_params[PW2_NPARAM*3]=
{
	//0
	
	 0x00EA, 0x01C8, 0x00A2, 0x005E, 0x01F4, 0x0045, 0x0091, 0x0066, 0x003B, 0x0027,-0x0011, 0x001B, 0x00FF, 0x007E, 0x00D1, 0x00F3, 0x008F, 0x0130, 0x018E,-0x00AC, 0x010C, 0x0008,-0x007E, 0x00A2, 0x000E,-0x0069,-0x0073,-0x0125,-0x0092, 0x0000, 0x0078,//Y
	 0x0006, 0x003D, 0x0031, 0x002F, 0x003F, 0x0015, 0x0011, 0x0036, 0x002E,-0x0022, 0x0011, 0x0034,-0x0007, 0x0012,-0x0018, 0x0012, 0x002F, 0x0000, 0x0000, 0x001C, 0x00A2, 0x02E1, 0x00C9,-0x00E0,-0x0068,-0x004E,-0x013E,-0x0012, 0x0001, 0x0000,-0x0046,//Cb
	 0x003D, 0x0036, 0x0006, 0x007E, 0x0012, 0x0007, 0x0007, 0x0005, 0x001E, 0x0000, 0x0028, 0x0055,-0x0020, 0x0020, 0x0005, 0x0011, 0x0034, 0x0000, 0x0004, 0x003E,-0x0100, 0x0001,-0x0086,-0x0041, 0x0051,-0x0080, 0x0004, 0x0002,-0x0003,-0x0003, 0x00D9,//Cr

	// 0x007D, 0x0040, 0x0039, 0x004A, 0x0007, 0x0003, 0x001D, 0x0007, 0x0000, 0x0007, 0x0002, 0x000F, 0x001B, 0x0018, 0x000A, 0x0008, 0x001C,-0x0008, 0x0004, 0x0005, 0x0006,-0x0022,-0x003B,-0x0041,-0x00C8,-0x0040,-0x0085,-0x0050, 0x0060,
	// 0x007A, 0x00F9, 0x0165, 0x00C2, 0x0036, 0x0100, 0x0054, 0x0000,-0x0081,-0x0078, 0x0020, 0x004D,-0x0010, 0x0028, 0x00BD, 0x009D, 0x0020,-0x0082,-0x003F, 0x0060, 0x002A, 0x0161,-0x004E,-0x001D, 0x0123,-0x0008,-0x0080, 0x0020, 0x003C,
	// 0x0010, 0x0039, 0x002E, 0x0037, 0x000E,-0x0010, 0x0014, 0x0008,-0x0007,-0x001C, 0x0074, 0x0019, 0x0010, 0x001B, 0x000D, 0x0047, 0x000A, 0x001C, 0x0008, 0x0004, 0x0023,-0x0012,-0x0156,-0x0074,-0x00A0,-0x0002,-0x0088,-0x0060, 0x0102,
};
//static int pred_w2_paeth2(int T, int L, int TL, int TR)
//{
//	int p=T+L-TL, closest=T;
//	if(abs(closest-p)>abs(L-p))
//		closest=L;
//	if(abs(closest-p)>abs(TL-p))
//		closest=TL;
//	if(abs(closest-p)>abs(TR-p))
//		closest=TR;
//	return closest;
//}
//static int pred_w2_select(int T, int L, int TL)
//{
//	int p=T+L-TL, pT=abs(p-T), pL=abs(p-L);
//	return pT<pL?L:T;
//}
//static pred_w2_cgrad(int T, int L, int TL)
//{
//	int vmin, vmax, grad;
//
//	if(T<L)
//		vmin=T, vmax=L;
//	else
//		vmin=L, vmax=T;
//	grad=T+L-TL;
//	grad=CLAMP(vmin, grad, vmax);
//	return grad;
//}
static int clamp4(int p, int a, int b, int c, int d)
{
	int vmin=a, vmax=a;
	if(vmin>b)vmin=b;
	if(vmin>c)vmin=c;
	if(vmin>d)vmin=d;
	if(vmax<b)vmax=b;
	if(vmax<c)vmax=c;
	if(vmax<d)vmax=d;
	p=CLAMP(vmin, p, vmax);
	return p;
}
//static int clip(int x, int nlevels)
//{
//	x=CLAMP(-nlevels, x, nlevels-1);
//	return x;
//}
static void pred_w2_prealloc(const int *src, int iw, int ih, int depth, int kc, short *params, int fwd, int enable_ma, int *dst, int *temp)//temp is (PW2_NPRED+1)*2w
{
	int
		errorbuflen, rowlen,
		*hireserror, *pred_errors[PW2_NPRED],
		idx, nlevels;
	const int *src2;

	errorbuflen=iw<<1;
	rowlen=iw<<2;
	hireserror=temp;
	for(int k=0;k<PW2_NPRED;++k)
		pred_errors[k]=temp+errorbuflen*(k+1);
	idx=kc;

	nlevels=1<<depth;
	src2=fwd?src:dst;
	memset(pw2_errors, 0, sizeof(pw2_errors));
	for(int ky=0;ky<ih;++ky)
	{
		//int pred_left=0, pred_left2=0;

		int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
		for(int kx=0;kx<iw;++kx, idx+=4)
		{
			int eT, eL, eTL, eTR, eT_L;
			//           T3
			//   T2L2    T2    T2R2
			//        TL T  TR TR2
			//L3 L2   L  X
			int
				cT6  =         ky-6>=0?src2[idx-rowlen*6   ]<<8:0,
			//	cT5  =         ky-5>=0?src2[idx-rowlen*5   ]<<8:0,

			//	cT4L3=kx-3>=0&&ky-4>=0?src2[idx-rowlen*4-12]<<8:0,
			//	cT4  =         ky-4>=0?src2[idx-rowlen*4   ]<<8:0,
			//	cT4R3=kx+3<iw&&ky-4>=0?src2[idx-rowlen*4+12]<<8:0,
				
			//	cT3L5=kx-5>=0&&ky-3>=0?src2[idx-rowlen*3-20]<<8:0,
			//	cT3L4=kx-4>=0&&ky-3>=0?src2[idx-rowlen*3-16]<<8:0,
			//	cT3L2=kx-2>=0&&ky-3>=0?src2[idx-rowlen*3- 8]<<8:0,
			//	cT3L =kx-1>=0&&ky-3>=0?src2[idx-rowlen*3- 4]<<8:0,
				cT3  =         ky-3>=0?src2[idx-rowlen*3   ]<<8:0,
			//	cT3R =kx+1<iw&&ky-3>=0?src2[idx-rowlen*3+ 4]<<8:0,
			//	cT3R2=kx+2<iw&&ky-3>=0?src2[idx-rowlen*3+ 8]<<8:0,
			//	cT3R3=kx+3<iw&&ky-3>=0?src2[idx-rowlen*3+12]<<8:0,
			//	cT3R4=kx+4<iw&&ky-3>=0?src2[idx-rowlen*3+16]<<8:0,
				
			//	cT2L3=kx-3>=0&&ky-2>=0?src2[idx-rowlen*2-12]<<8:0,
			//	cT2L2=kx-2>=0&&ky-2>=0?src2[idx-rowlen*2- 8]<<8:0,
				cT2L =kx-1>=0&&ky-2>=0?src2[idx-rowlen*2- 4]<<8:0,
				cT2  =         ky-2>=0?src2[idx-rowlen*2   ]<<8:0,
				cT2R =kx+1<iw&&ky-2>=0?src2[idx-rowlen*2+ 4]<<8:0,
			//	cT2R2=kx+2<iw&&ky-2>=0?src2[idx-rowlen*2+ 8]<<8:0,
			//	cT2R3=kx+3<iw&&ky-2>=0?src2[idx-rowlen*2+12]<<8:0,
			//	cT2R4=kx+4<iw&&ky-2>=0?src2[idx-rowlen*2+16]<<8:0,
				
			//	cTL3 =kx-3>=0&&ky-1>=0?src2[idx-rowlen  -12]<<8:0,
			//	cTL2 =kx-2>=0&&ky-1>=0?src2[idx-rowlen  - 8]<<8:0,
				cTL  =kx-1>=0&&ky-1>=0?src2[idx-rowlen  - 4]<<8:0,
				cT   =         ky-1>=0?src2[idx-rowlen     ]<<8:0,
				cTR  =kx+1<iw&&ky-1>=0?src2[idx-rowlen  + 4]<<8:0,
				cTR2 =kx+2<iw&&ky-1>=0?src2[idx-rowlen  + 8]<<8:0,
			//	cTR3 =kx+3<iw&&ky-1>=0?src2[idx-rowlen  +12]<<8:0,
				cTR4 =kx+4<iw&&ky-1>=0?src2[idx-rowlen  +16]<<8:0,
				cTR5 =kx+5<iw&&ky-1>=0?src2[idx-rowlen  +20]<<8:0,
				cTR6 =kx+6<iw&&ky-1>=0?src2[idx-rowlen  +24]<<8:0,
				cTR7 =kx+7<iw&&ky-1>=0?src2[idx-rowlen  +28]<<8:0,

				cL6  =kx-6>=0         ?src2[idx         -24]<<8:0,
			//	cL5  =kx-5>=0         ?src2[idx         -20]<<8:0,
				cL4  =kx-4>=0         ?src2[idx         -16]<<8:0,
			//	cL3  =kx-2>=0         ?src2[idx         -12]<<8:0,
			//	cL2  =kx-2>=0         ?src2[idx         - 8]<<8:0,
				cL   =kx-1>=0         ?src2[idx         - 4]<<8:0;

			//w0   w1   w2   w3
			//p3Ca p3Cb p3Cc p3Cd p3Ce
			//p1C  p2c

			//if(kx==(iw>>1)&&ky==(ih>>1))//
			//	kx=iw>>1;
			
			int weights[PW2_NPRED];//fixed 23.8 bit
			for(int k=0;k<PW2_NPRED;++k)
			{
				//eNW + eN + eNE
				int e2=
					 (ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0)
					+(ky-1>=0?pred_errors[k][prevrow+kx]:0)
					+(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0);
				weights[k]=(params[k]<<8)/(e2+1);
			}

			//TL T TR
			//L  X
			eT=ky-1>=0?hireserror[prevrow+kx]:0;
			eL=kx-1>=0?hireserror[currrow+kx-1]:0;
			eTL=ky-1>=0&&kx-1>=0?hireserror[prevrow+kx-1]:0;
			eTR=ky-1>=0&&kx+1<iw?hireserror[prevrow+kx+1]:0;
			eT_L=eT+eL;

			//pred_left=pred_left+((cL-pred_left)*params[PW2_NPRED+11]>>8);
			//pred_left2=pred_left2+((cL-pred_left2)*abs(eL)>>16);

			//int pred_w=(cT*(0x10000-eT)+cL*(0x10000-eL)+cTL*(0x10000-eTL)+cTR*(0x10000-eTR))>>16;
			{
				int vmin, vmax, curr;

				int predictions[PW2_NPRED]=//fixed 23.8 bit
				{
					//from jxl
					cT-((eTL*params[PW2_NPRED]+eT*params[PW2_NPRED+1]+eTR*params[PW2_NPRED+2]+(cT2-cT)*params[PW2_NPRED+3]+(cTL-cL)*params[PW2_NPRED+4])>>8),//k13: 1.998458 how many optimizations?
					//k13: 1.737220, 1.837736, 1.842377, 1.847278, 1.865093, 1.861586, 1.866601, 1.872176, 1.878888, 1.883146, 1.883862, 1.883863, <

					cL-((eT_L+eTL)*params[PW2_NPRED+5]>>8),
					//k13: 1.737220, 1.932766, 1.960469, 1.966698, 1.970106, 1.971381, 1.971891, 1.972187, 1.972368, 1.972492, 1.972565, 1.972580, 1.972594, 1.972602, <

					cT-((eT_L+eTR)*params[PW2_NPRED+6]>>8),
					//k13: 1.737220, 1.934911, 1.951553, 1.962075, 1.973068, 1.979630, 1.983403, 1.985205, 1.986109, 1.986488, 1.986609, 1.986677, 1.986688, 1.986689, <

					cL+cTR-cT,
					//k13: 1.909206, 1.918954, 1.946385, 1.963439, 1.970334, 1.981971, 1.983898, 1.984527, 1.985102, 1.985773, 1.986008, 1.986331, 1.986678, 1.986722, 1.986756			...1.998458


					cL-(eL*params[PW2_NPRED+7]>>8),
					cT-(eT*params[PW2_NPRED+8]>>8),
					//k13: 1.909206, 1.922112, 1.931268, 1.954690, 1.964123, 1.978742, 1.981578, 1.984138, 1.985535, 1.986711, 1.987659, 1.988190, 1.988474, 1.988532, 1.988550
					cTL-(eTL*params[PW2_NPRED+9]>>8),
					cTR-(eTR*params[PW2_NPRED+10]>>8),
					//k13: 1.909206, 1.921490, 1.932631, 1.949766, 1.950930, 1.951645, 1.951977, 1.960758, 1.967595, 1.969669, 1.972408, 1.973050, 1.973506, 1.974268, 1.975184			...1.977183


					//paq8px by Matt Mahoney

					//k13: 1.737220, 1.958513, 1.973267, 1.979685, 1.983374, 1.985860, 1.987622, 1.989731, 1.991147, 1.992018, 1.992707, 1.993444, 1.994374, 1.995238, 1.996056,   1.996876, 1.997423, 1.997708, 1.997946, 1.998162, 1.998320, 1.998364, 1.998611, 1.998815, 1.998948, 1.999125, 1.999207, 1.999222, 1.999229, 1.999235, 1.999241, 1.999242, 1.999247, 1.999248, 1.999250, 1.999251, <

					clamp4(cL+cT-cTL, cL, cTL, cT, cTR),//0
					//clip(cL+cT-cTL),//1
					clamp4(cL+cTR-cT, cL, cTL, cT, cTR),//2
					//clip(cL+cTR-cT),//3
					clamp4(cT+cTL-cT2L, cL, cTL, cT, cTR),//4
					//clip(cT+cTL-cT2L),//5
					clamp4(cT+cTR-cT2R, cL, cT, cTR, cTR2),//6
					//clip(cT+cTR-cT2R),//7
					(cL+cTR2)>>1,//8
					//clip(cT*3-cT2*3+cT3),//9
					//clip(cL*3-cL2*3+cL3),//10
					//(cL+clip(cTR*3-cT2R*3+cT3R))>>1,//11
					//(cL+clip(cTR2*3-cT2R3*3+cT3R4))>>1,//12
					//clip(cT2+cT4-cT6),//13
					//clip(cL2+cL4-cL6),//14
					//clip((cT5-6*cT4+15*cT3-20*cT2+15*cT+clamp4(cL*2-cTL2, cL, cTL, cT, cT2))/6),//15
					//clip((-3*cL2+8*cL+clamp4(3*cTR2-3*cT2R2+cT3R2, cTR, cTR2, cTR3, cTR4))/6),//16
					//clip(cT2+cTL-cT3L),//17
					//clip(cT2+cTR-cT3R),//18
					//clip((cL*2+cTL) - (cL2*2+cTL2) + cL3),//19
					//clip(3*(cTL+cTL2)/2-cT2L3*3+(cT3L4+cT3L5)/2),//20
					//clip(cTR2+cTR-cT2R3),//21
					//clip(cTL2+cL2-cL4),//22
					//clip(((cL+cTL)*3-cTL2*6+cTL3+cT2L3)/2),//23
					//clip((cTR*2+cTR2) - (cT2R2+cT3R2*2) + cT4R3),//24
					cT6,//25
					(cTR4+cTR6)>>1,//26
					(cL4+cL6)>>1,//27
					(cL+cT+cTR5+cTR7)>>2,//28
					//clip(cTR3+cL-cTR2),//29
					//clip(4*cT3-3*cT4),//30
					//clip(cT+cT2-cT3),//31
					//clip(cL+cL2-cL3),//32
					//clip(cL+cTR2-cTR),//33
					//clip(cL2+cTR2-cT),//34
					//(clip(cL*2-cTL)+clip(cL*2-cTL2)+cT+cTR)>>2,//35
					clamp4(cT*2-cT2, cL, cT, cTR, cTR2),//36
					(cT+cT3)>>1,//37
					//clip(cT2+cL-cT2L),//38
					//clip(cTR2+cT-cT2R2),//39
					//clip((4*cL3-15*cL2+20*cL+clip(cTR2*2-cT2R2))/10),//40
					//clip((cT3R3-4*cT2R2+6*cTR+clip(cL*3-cTL*3+cT2L))/4),//41
					//clip((cT*2+cTR) - (cT2+2*cT2R) + cT3R),//42
					//clip((cTL*2+cT2L) - (cT2L2+cT3L2*2) + cT4L3),//43
					//clip(cT2L2+cL-cT2L3),//44
					//clip((-cT4+5*cT3-10*cT2+10*cT+clip(cL*4-cTL2*6+cT2L3*4-cT3L4))/5),//45
					//clip(cTR2+clip(cTR3*2-cT2R4-cTR4)),//46
					//clip(cTL+cL-cTL2),//47
					//clip((cT*2+cTL) - (cT2+2*cT2L) + cT3L),//48
					//clip(cT2+clip(cTR2*2-cT2R3) - cT2R),//49
					//clip((-cL4+5*cL3-10*cL2+10*cL+clip(cTR*2-cT2R))/5),//50
					//clip((-cL5+4*cL4-5*cL3+5*cL+clip(cTR*2-cT2R))>>2),//51
					//clip((cL3-4*cL2+6*cL+clip(cTR*3-cT2R*3+cT3R))>>2),//52
					//clip((-cT2R2+3*cTR+clip(4*cL-6*cTL+4*cT2L-cT3L))/3),//53
					((cL+cT)*3-cTL*2)>>2,//54
				};

				long long pred, sum=0;
				for(int k=0;k<PW2_NPRED;++k)
					sum+=weights[k];
				if(sum)
				{
					pred=(sum>>1)-1;
					for(int k=0;k<PW2_NPRED;++k)
						pred+=predictions[k]*weights[k];
					pred/=sum;
				}
				else
					pred=predictions[0];

				vmin=cL;
				vmax=cL;
				if(vmin>cT)vmin=cT;
				if(vmax<cT)vmax=cT;
				if(vmin>cTR)vmin=cTR;
				if(vmax<cTR)vmax=cTR;
				pred=CLAMP(vmin, pred, vmax);


				if(fwd)
				{
					curr=src[idx]<<8;
					dst[idx]=src[idx]-(int)((pred+127)>>8);
					if(enable_ma)
						dst[idx]=((dst[idx]+(nlevels>>1))&(nlevels-1))-(nlevels>>1);
				}
				else
				{
					dst[idx]=src[idx]+(int)((pred+127)>>8);
					if(enable_ma)
						dst[idx]=((dst[idx]+(nlevels>>1))&(nlevels-1))-(nlevels>>1);
					curr=dst[idx]<<8;
				}

				hireserror[currrow+kx]=curr-(int)pred;
				for(int k=0;k<PW2_NPRED;++k)
				{
					int e=abs(curr-predictions[k]);
					pw2_errors[k]+=e;//
					pred_errors[k][currrow+kx]=e;
					if(kx+1<iw)//add current error to eNE, such that eN = eN0 + eW
						pred_errors[k][prevrow+kx+1]+=e;
				}
			}
		}
	}
	for(int k=0;k<PW2_NPRED;++k)
		pw2_errors[k]/=iw*ih*256;
}
static double pred_w2_calcloss(const int *src, int iw, int ih, int depth, int src_depth, int kc, short *params, int *temp, int *dst, int *hist)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	pred_w2_prealloc(src, iw, ih, depth, kc, params, 1, 1, dst, temp);
	calc_histogram(dst, iw, ih, kc, 0, iw, 0, ih, depth, hist, 0);
	{
		double entropy=calc_entropy(hist, 1<<depth, (int)res);

		double invCR=entropy/src_depth, csize=res*invCR;
		return csize;
	}
}
void pred_w2_opt_v2(Image *src, short *params, int loud)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int *buf3=(int*)malloc((size_t)res*sizeof(int[4]));
	int *temp=(int*)malloc((size_t)src->iw*(PW2_NPRED+1)*2*sizeof(int));
	int maxdepth=calc_maxdepth(src, 0);
	int nlevels=1<<maxdepth;
	int *hist=(int*)malloc(nlevels*sizeof(int));
	char title0[256];
	int steps[]={128, 64, 32, 16, 8, 4, 2, 1};

	(void)loud;

	if(!buf3||!temp||!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	get_window_title(title0, 256);
	for(int kc=0;kc<3;++kc)
	{
		short *param=params+kc*PW2_NPARAM;
		double csize0=pred_w2_calcloss(src->data, src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
		for(int ks=0;ks<(int)_countof(steps);++ks)
		{
			int step=steps[ks];
			double bestcsize=csize0;
			int bestidx=0, beststep=0;
			for(int idx=0;idx<PW2_NPARAM;++idx)
			{
				double csize;
				short prev;

				prev=param[idx];
				param[idx]+=(short)step;
				csize=pred_w2_calcloss(src->data, src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=step;

				prev=param[idx];
				param[idx]-=(short)step;
				if(idx<4&&param[idx]<1)
					param[idx]=1;
				csize=pred_w2_calcloss(src->data, src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=-step;

				if(loud_transforms)
					set_window_title("Ch%d csize %lf [%d/%d %d/%d]...", kc, csize0, (int)(kc*_countof(steps)+ks+1), (int)(_countof(steps)*3), idx+1, PW2_NPARAM);//
			}
			if(csize0>bestcsize)
			{
				csize0=bestcsize;

				param[bestidx]+=(short)beststep;
				if(bestidx<4&&param[bestidx]<1)
					param[bestidx]=1;
			}
			//set_window_title("Ch%d csize %lf [%d/%d]...", kc, csize0, (int)(kc*_countof(steps)+ks+1), (int)(_countof(steps)*3));//
		}
	}
	free(hist);
	free(temp);
	free(buf3);
	if(loud_transforms)
		set_window_title("%s", title0);//
}
void pred_w2_apply(Image *src, int fwd, int enable_ma, short *params)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int *temp=(int*)malloc((size_t)src->iw*(PW2_NPRED+1)*2*sizeof(int));
	int *buf2=(int*)malloc(res*sizeof(int[4]));
	if(!temp||!buf2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	
	pred_w2_prealloc(src->data, src->iw, src->ih, src->depth[0], 0, params             , fwd, enable_ma, buf2, temp);
	pred_w2_prealloc(src->data, src->iw, src->ih, src->depth[1], 1, params+PW2_NPARAM  , fwd, enable_ma, buf2, temp);
	pred_w2_prealloc(src->data, src->iw, src->ih, src->depth[2], 2, params+PW2_NPARAM*2, fwd, enable_ma, buf2, temp);
	if(src->depth[3])
		pred_w2_prealloc(src->data, src->iw, src->ih, src->depth[3], 3, params+PW2_NPARAM, fwd, enable_ma, buf2, temp);

	for(ptrdiff_t k=0;k<res;++k)
	{
		src->data[k<<2|0]=buf2[k<<2|0];
		src->data[k<<2|1]=buf2[k<<2|1];
		src->data[k<<2|2]=buf2[k<<2|2];
		if(src->depth[3])
			src->data[k<<2|3]=buf2[k<<2|3];
	}
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}

	free(temp);
	free(buf2);
}


short jxlparams_i16[33]=//signed fixed 7.8 bit
{
	0x0DB8,  0x0E22,  0x181F,  0x0BF3, -0x005C, -0x005B,  0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,//Y
	0x064C,  0x0F31,  0x1040,  0x0BF8, -0x0007, -0x000D, -0x0085, -0x0063, -0x00A2, -0x0017,  0x00F2,//Cb
	0x0B37,  0x110B,  0x121B,  0x0BFC, -0x0001,  0x000E, -0x0188, -0x00E7, -0x00BB, -0x004A,  0x00BA,//Cr

	//0x0DD3, 0x1002, 0x11F9, 0x0BEC, 0x0002,-0x0029,-0x010E,-0x00E8,-0x00A2,-0x0089, 0x0064,
	//0x0FAF, 0x0C9E, 0x139F, 0x0E13,-0x009D,-0x005C, 0x0041, 0x0013, 0x00DD,-0x0004,-0x00C3,
	//0x09A8, 0x0F00, 0x1191, 0x0B9E,-0x0026,-0x0027,-0x0086,-0x0084,-0x005F,-0x0038, 0x00B2,

	//0x1000, 0x1000, 0x1000, 0x1000, 0, 0, 0, 0, 0, 0, 0,
	//0x1000, 0x1000, 0x1000, 0x1000, 0, 0, 0, 0, 0, 0, 0,
	//0x1000, 0x1000, 0x1000, 0x1000, 0, 0, 0, 0, 0, 0, 0,

	//0x0C7B, 0x0B5C, 0x0A14, 0x0B33,  -0x0148,-0x0029, 0, 0, 0, 0, 0,
	//0x0A14, 0x0829, 0x1548, 0x0CA4,   0x047B, 0x0052, 0, 0, 0, 0, 0,
	//0x0B33, 0x0C29, 0x0DC3, 0x119A,  -0x0148,-0x00F6, 0, 0, 0, 0, 0,

	//0x0C7B, 0x0B5C, 0x0A14, 0x0B33,-0x0148,-0x0029, 0x0971, 0x01EC,-0x01C3, 0x047B, 0x0AB8,
	//0x0A14, 0x0829, 0x1548, 0x0CA4, 0x047B, 0x0052,-0x011F, 0x0000, 0x0029, 0x063D, 0x0266,
	//0x0B33, 0x0C29, 0x0DC3, 0x119A,-0x0148,-0x00F6, 0x0614, 0x00A4,-0x007B, 0x019A, 0x0E8F,

	//0x63, 0x5A, 0x50, 0x59,    -0x0A, -0x01,  0x4B, 0x0F, -0x0E, -0x24, 0x55,
	//0x50, 0x41, 0xA9, 0x64,     0x24,  0x03, -0x09,    0,  0x01,  0x32, 0x13,
	//0x59, 0x61, 0x6D, 0x8C,    -0x0A, -0x08,  0x30, 0x05, -0x04,  0x0D, 0x74,
};
static void pred_jxl_prealloc(const int *src, int iw, int ih, int depth, int kc, short *params, int fwd, int enable_ma, int *dst, int *temp_w10)
{
	int errorbuflen=iw<<1, rowlen=iw<<2;
	int *hireserror=temp_w10, *pred_errors[]=
	{
		temp_w10+errorbuflen*1,
		temp_w10+errorbuflen*2,
		temp_w10+errorbuflen*3,
		temp_w10+errorbuflen*4,
	};
	int idx=kc;
	const int *pixels=fwd?src:dst;
//	const int *errors=fwd?dst:src;

	memset(temp_w10, 0, iw*sizeof(int[10]));
	for(int ky=0;ky<ih;++ky)
	{
		int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
		for(int kx=0;kx<iw;++kx, idx+=4)
		{
			int pred, curr;
			
			int
				ctt      =         ky-2>=0?pixels[idx-rowlen*2]:0,
				ctopleft =kx-1>=0&&ky-1>=0?pixels[idx-rowlen-4]:0,
				ctop     =kx  <iw&&ky-1>=0?pixels[idx-rowlen  ]:0,
				ctopright=kx+1<iw&&ky-1>=0?pixels[idx-rowlen+4]:0,
				cleft    =kx-1>=0         ?pixels[idx       -4]:0;

			//if(kx==(iw>>1)&&ky==(ih>>1))
			//	kx=iw>>1;

			//w0   w1   w2   w3
			//p3Ca p3Cb p3Cc p3Cd p3Ce
			//p1C  p2c
			
			int weights[4];//fixed 23.8 bit
			for(int k=0;k<4;++k)
			{
				//	eNW	eN	eNE
				//eWW	eW	?
				int e2=
					+(ky-1>=0?pred_errors[k][prevrow+kx]:0)			//eN + eW
					+(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0)	//eNE
					+(ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0);	//eNW + eWW
				weights[k]=(params[k]<<8)/(e2+1);
			}
			{
				int
					etop=ky-1>=0?hireserror[prevrow+kx]:0,
					eleft=kx-1>=0?hireserror[currrow+kx-1]:0,
					etopleft=ky-1>=0&&kx-1>=0?hireserror[prevrow+kx-1]:0,
					etopright=ky-1>=0&&kx+1<iw?hireserror[prevrow+kx+1]:0,
					etopplusleft=etop+eleft;
				long long predictions[]=//fixed 23.8 bit
				{
					(cleft+ctopright-ctop)<<8,
					(ctop<<8)-((etopplusleft+etopright)*params[4]>>8),
					(cleft<<8)-((etopplusleft+etopleft)*params[5]>>8),
					(ctop<<8)-(((etopleft*params[6]+etop*params[7]+etopright*params[8])>>8)+(ctt-ctop)*params[9]+(ctopleft-cleft)*params[10]),
				};

				int sum=weights[0]+weights[1]+weights[2]+weights[3];
				if(sum)
					pred=(int)((predictions[0]*weights[0]+predictions[1]*weights[1]+predictions[2]*weights[2]+predictions[3]*weights[3]+(sum>>1)-1)/sum);
				else
					pred=(int)predictions[0];

				//int vmin=cleft, vmax=cleft;
				//if(vmin>ctopright)vmin=ctopright;
				//if(vmax<ctopright)vmax=ctopright;
				//if(vmin>ctop)vmin=ctop;
				//if(vmax<ctop)vmax=ctop;
				//vmin<<=8;
				//vmax<<=8;

				//if(kc==0&&kx==0&&ky==1)//
				//if(kc==0&&kx==256&&ky==256)//
				//	printf("");

				//pred=CLAMP(vmin, pred, vmax);
				CLAMP3_32(pred, pred, ctop<<8, cleft<<8, ctopright<<8);
				{
					int pred_final=(pred+127)>>8;
					pred_final^=-fwd;
					pred_final+=fwd;
					pred_final+=src[idx];
					if(enable_ma)
					{
						pred_final<<=32-depth;
						pred_final>>=32-depth;
					}
					dst[idx]=pred_final;
					curr=pixels[idx]<<8;
					//if(fwd)
					//{
					//	dst[idx]=src[idx]-((pred+127)>>8);
					//	if(enable_ma)
					//		dst[idx]=((dst[idx]+(nlevels>>1))&(nlevels-1))-(nlevels>>1);
					//	curr=src[idx]<<8;
					//}
					//else
					//{
					//	dst[idx]=src[idx]+((pred+127)>>8);
					//	if(enable_ma)
					//		dst[idx]=((dst[idx]+(nlevels>>1))&(nlevels-1))-(nlevels>>1);
					//	curr=dst[idx]<<8;
					//}

					hireserror[currrow+kx]=curr-pred;
					for(int k=0;k<4;++k)
					{
						int e=abs(curr-(int)predictions[k]);
						pred_errors[k][currrow+kx]=e;
						if(kx+1<iw)
							pred_errors[k][prevrow+kx+1]+=e;//eNE+=ecurr
					}
				}
			}
		}
	}
}
static double pred_jxl_calcloss(const int *src, int iw, int ih, int depth, int src_depth, int kc, short *params, int *temp, int *dst, int *hist)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	pred_jxl_prealloc(src, iw, ih, depth, kc, params, 1, 1, dst, temp);
	//addhalf((unsigned char*)dst+kc, iw, ih, 1, 4);
	calc_histogram(dst, iw, ih, kc, 0, iw, 0, ih, depth, hist, 0);
	{
		double entropy=calc_entropy(hist, 1<<depth, (int)res);
		double invCR=entropy/src_depth, csize=res*invCR;
		return csize;
	}
}
void pred_jxl_opt_v2(Image *src, short *params, int loud)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int *buf3=(int*)malloc((size_t)res*sizeof(int[4]));
	int *temp=(int*)malloc((size_t)src->iw*10*sizeof(int));
	int maxdepth=calc_maxdepth(src, 0);
	int nlevels=1<<maxdepth;
	int *hist=(int*)malloc(nlevels*sizeof(int));
	int steps[]={256, 128, 64, 32, 16, 8, 4, 2, 1};

	(void)loud;

	if(!buf3||!temp||!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(int kc=0;kc<3;++kc)
	{
		short *param=params+kc*11;
		double csize0=pred_jxl_calcloss(src->data, src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
		for(int ks=0;ks<(int)_countof(steps);++ks)
		{
			int step=steps[ks];
			double bestcsize=csize0;
			int bestidx=0, beststep=0;
			for(int idx=0;idx<11;++idx)
			{
				double csize;
				short prev;

				prev=param[idx];
				param[idx]+=(short)step;
				csize=pred_jxl_calcloss(src->data, src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=step;

				prev=param[idx];
				param[idx]-=(short)step;
				if(idx<4&&param[idx]<1)
					param[idx]=1;
				csize=pred_jxl_calcloss(src->data, src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=-step;
			}
			if(csize0>bestcsize)
			{
				csize0=bestcsize;

				param[bestidx]+=(short)beststep;
				if(bestidx<4&&param[bestidx]<1)
					param[bestidx]=1;
			}
			//set_window_title("%lf", csize0);//
		}
	}
	free(hist);
	free(temp);
	free(buf3);
}
void pred_jxl_apply(Image *src, int fwd, int enable_ma, short *params)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int *temp=(int*)malloc((size_t)src->iw*10*sizeof(int));
	int *buf2=(int*)malloc(res*sizeof(int[4]));
	if(!temp||!buf2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	
	pred_jxl_prealloc(src->data, src->iw, src->ih, src->depth[0], 0, params   , fwd, enable_ma, buf2, temp);
	pred_jxl_prealloc(src->data, src->iw, src->ih, src->depth[1], 1, params+11, fwd, enable_ma, buf2, temp);
	pred_jxl_prealloc(src->data, src->iw, src->ih, src->depth[2], 2, params+22, fwd, enable_ma, buf2, temp);
	if(src->depth[3])
		pred_jxl_prealloc(src->data, src->iw, src->ih, src->depth[3], 3, params+11, fwd, enable_ma, buf2, temp);

	for(int k=0;k<res;++k)
	{
		src->data[k<<2  ]=buf2[k<<2  ];
		src->data[k<<2|1]=buf2[k<<2|1];
		src->data[k<<2|2]=buf2[k<<2|2];
		if(src->depth[3])
			src->data[k<<2|3]=buf2[k<<2|3];
	}

	free(temp);
	free(buf2);
}

void pred_jmj_apply(Image *src, int fwd, int enable_ma)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int *temp=(int*)malloc((size_t)src->iw*(PW2_NPRED+1)*2*sizeof(int));
	int *buf2=(int*)malloc(res*sizeof(int[4]));
	if(!temp||!buf2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	
	pred_jxl_prealloc(src->data, src->iw, src->ih, src->depth[0], 0, jxlparams_i16        , fwd, enable_ma, buf2, temp);
	pred_w2_prealloc (src->data, src->iw, src->ih, src->depth[1], 1, pw2_params+PW2_NPARAM, fwd, enable_ma, buf2, temp);
	pred_jxl_prealloc(src->data, src->iw, src->ih, src->depth[2], 2, jxlparams_i16+22     , fwd, enable_ma, buf2, temp);
	if(src->depth[3])
		pred_jxl_prealloc(src->data, src->iw, src->ih, src->depth[3], 3, jxlparams_i16+22, fwd, enable_ma, buf2, temp);

	for(int k=0;k<res;++k)
	{
		src->data[k<<2  ]=buf2[k<<2  ];
		src->data[k<<2|1]=buf2[k<<2|1];
		src->data[k<<2|2]=buf2[k<<2|2];
		if(src->depth[3])
			src->data[k<<2|3]=buf2[k<<2|3];
	}

	free(temp);
	free(buf2);
}


//CUSTOM3

	#define C3_3D

#define C3_OPT_NCOMP (C3_NPARAMS>>5)
Custom3Params c3_params={0};
static int fast_dot(const short *a, const short *b, int count)
{
	int k;
	__m256i sum=_mm256_setzero_si256();
	for(k=0;k<count-31;k+=16)//https://stackoverflow.com/questions/62041400/inner-product-of-two-16bit-integer-vectors-with-avx2-in-c
	{
		__m256i va=_mm256_loadu_si256((const __m256i*)(a+k));
		__m256i vb=_mm256_loadu_si256((const __m256i*)(b+k));
		va=_mm256_madd_epi16(va, vb);
		sum=_mm256_add_epi32(sum, va);
	}
	{
		__m128i s2=_mm_add_epi32(_mm256_extracti128_si256(sum, 1), _mm256_castsi256_si128(sum));
		__m128i hi=_mm_shuffle_epi32(s2, _MM_SHUFFLE(2, 1, 3, 2));
		s2=_mm_add_epi32(s2, hi);
		s2=_mm_hadd_epi32(s2, s2);
		{
			int s3=_mm_extract_epi32(s2, 0);
			for(;k<count;++k)
				s3+=a[k]*b[k];
			return s3;
		}
	}
}
static int custom3_loadnb(const int *pixels, const int *errors, int iw, int ih, int kc, int kx, int ky, short *nb)
{
	int idx=-1;
	for(int ky2=-C3_REACH;ky2<0;++ky2)
	{
		for(int kx2=-C3_REACH;kx2<=C3_REACH;++kx2)
		{
			if((unsigned)(kx+kx2)<(unsigned)iw&&(unsigned)(ky+ky2)<(unsigned)ih)
			{
				int idx2=(iw*(ky+ky2)+kx+kx2)<<2|kc;
				nb[++idx]=(short)pixels[idx2];
				nb[++idx]=(short)errors[idx2];
			}
			else
			{
				nb[++idx]=0;
				nb[++idx]=0;
			}
		}
	}
	for(int kx2=-C3_REACH;kx2<0;++kx2)
	{
		if((unsigned)(kx+kx2)<(unsigned)iw)
		{
			int idx2=(iw*ky+kx+kx2)<<2|kc;
			nb[++idx]=(short)pixels[idx2];
			nb[++idx]=(short)errors[idx2];
		}
		else
		{
			nb[++idx]=0;
			nb[++idx]=0;
		}
	}
	return ++idx;
}
static void custom3_prealloc(const int *src, int iw, int ih, const char *depths, int fwd, int enable_ma, Custom3Params const *params, int *dst)
{
	const int *pixels=fwd?src:dst, *errors=fwd?dst:src;
	int pred, idx;
	const short *coeffs[]=
	{
		params->c00, params->c01, params->c02,
		params->c10, params->c11, params->c12,
		params->c20, params->c21, params->c22,
	};
	short nb[3][C3_NNB+2]={0};
	int nlevels[]=
	{
		1<<depths[0],
		1<<depths[1],
		1<<depths[2],
	};
	for(int ky=0;ky<ih;++ky)
	{
		for(int kx=0;kx<iw;++kx)
		{
#if 1
			int clampers[6]=
			{
				(nlevels[0]>>1)-1, -(nlevels[0]>>1),
				(nlevels[1]>>1)-1, -(nlevels[1]>>1),
				(nlevels[2]>>1)-1, -(nlevels[2]>>1),
			};
			int val;
			int count[3], idx2;
#define UPDATE_CLAMPER(X, A, B) val=X; UPDATE_MIN(clampers[A], val); UPDATE_MAX(clampers[B], val);
			if(kx)
			{
				UPDATE_CLAMPER(pixels[(iw*ky+kx-1)<<2|0], 0, 1)
				UPDATE_CLAMPER(pixels[(iw*ky+kx-1)<<2|1]-pixels[(iw*ky+kx-1)<<2|0], 2, 3)
				UPDATE_CLAMPER(pixels[(iw*ky+kx-1)<<2|2]-pixels[(iw*ky+kx-1)<<2|0], 4, 5)
			}
			if(ky)
			{
				UPDATE_CLAMPER(pixels[(iw*(ky-1)+kx)<<2|0], 0, 1)
				UPDATE_CLAMPER(pixels[(iw*(ky-1)+kx)<<2|1]-pixels[(iw*(ky-1)+kx)<<2|0], 2, 3)
				UPDATE_CLAMPER(pixels[(iw*(ky-1)+kx)<<2|2]-pixels[(iw*(ky-1)+kx)<<2|0], 4, 5)
				if(kx)
				{
					UPDATE_CLAMPER(pixels[(iw*(ky-1)+kx-1)<<2|0], 0, 1)
					UPDATE_CLAMPER(pixels[(iw*(ky-1)+kx-1)<<2|1]-pixels[(iw*(ky-1)+kx-1)<<2|0], 2, 3)
					UPDATE_CLAMPER(pixels[(iw*(ky-1)+kx-1)<<2|2]-pixels[(iw*(ky-1)+kx-1)<<2|0], 4, 5)
				}
			}
#undef  UPDATE_CLAMPER
			if(clampers[0]>clampers[1])
				clampers[0]=-(nlevels[0]>>1), clampers[1]=(nlevels[0]>>1)-1;
			if(clampers[2]>clampers[3])
				clampers[2]=-(nlevels[1]>>1), clampers[3]=(nlevels[1]>>1)-1;
			if(clampers[4]>clampers[5])
				clampers[4]=-(nlevels[2]>>1), clampers[5]=(nlevels[2]>>1)-1;
#endif
			for(int kc=0;kc<3;++kc)
				count[kc]=custom3_loadnb(pixels, errors, iw, ih, kc, kx, ky, nb[kc]);
#ifdef C3_3D
			for(int kc=1;kc<3;++kc)
			{
				for(int kn=0;kn<count[0];++kn)
					nb[kc][kn]-=nb[0][kn];
			}
#endif
			idx=(iw*ky+kx)<<2;
			idx2=0;
			
			for(int kdst=0;kdst<3;++kdst)
			{
				pred=0;
				for(int kc=0;kc<3;++kc)
					pred+=fast_dot(coeffs[idx2+kc], nb[kc], count[kc]);
				pred+=1<<13;
				pred>>=14;
				pred=CLAMP(clampers[kdst<<1|0], pred, clampers[kdst<<1|1]);
#ifdef C3_3D
				if(kdst)
					pred+=nb[0][C3_NNB];
				pred=CLAMP(-(nlevels[kdst]>>1), pred, (nlevels[kdst]>>1)-1);
#endif

				pred^=-fwd;
				pred+=fwd;
				pred+=src[idx];
				if(enable_ma)
				{
					pred+=nlevels[kdst]>>1;
					pred&=nlevels[kdst]-1;
					pred-=nlevels[kdst]>>1;
				}
				dst[idx]=pred;

				nb[kdst][C3_NNB  ]=(short)pixels[idx];
				nb[kdst][C3_NNB+1]=(short)errors[idx];
				count[kdst]+=2;
				++idx;
				idx2+=3;
			}
		}
	}
}

void custom3_apply(Image *src, int fwd, int enable_ma, Custom3Params const *params)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int *temp=(int*)malloc(res*sizeof(int[4]));
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(temp, src->data, res*sizeof(int[4]));//copy alpha

	custom3_prealloc(src->data, src->iw, src->ih, src->depth, fwd, enable_ma, params, temp);

	memcpy(src->data, temp, res*sizeof(int[4]));
	free(temp);
}

typedef struct Custom3OptInfoStruct
{
	double invCR[4];//loss == invCR[3]=(invCR[0]+invCR[1]+invCR[2])/3
	short params[C3_NPARAMS];
} Custom3OptInfo;
static void custom3_calcloss(const int *src, int iw, int ih, const char *depths, const char *src_depths, Custom3OptInfo *info, int *temp, int *hist)
{
	ptrdiff_t res=(ptrdiff_t)iw*ih;
	custom3_prealloc(src, iw, ih, depths, 1, 1, (Custom3Params*)info->params, temp);
	
	info->invCR[3]=0;
	for(int kc=0;kc<3;++kc)
	{
		calc_histogram(temp, iw, ih, kc, 0, iw, 0, ih, depths[kc], hist, 0);
		{
			double entropy=calc_entropy(hist, 1<<depths[kc], (int)res);
			info->invCR[kc]=entropy/src_depths[kc];
			info->invCR[3]+=info->invCR[kc];
		}
	}
	info->invCR[3]/=3;
}
void custom3_opt(Image const *src, Custom3Params *srcparams, int niter, int maskbits, int loud, double *loss)
{
	static int call_idx=0;

	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int maxdepth=calc_maxdepth(src, 0), maxlevels=1<<maxdepth;
	int *temp=(int*)malloc(res*sizeof(int[4]));
	int *hist=(int*)malloc(maxlevels*sizeof(int));
	Custom3OptInfo info;
	double invCR[4], loss0;
	int shakethreshold=C3_NPARAMS;

	++call_idx;
	if(call_idx==1)
		DisableProcessWindowsGhosting();
	if(!temp||!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(temp, 0xFF, res*sizeof(int[4]));//set alpha for preview
	memcpy(info.params, srcparams, sizeof(info.params));
#define CALC_LOSS() custom3_calcloss(src->data, src->iw, src->ih, src->depth, src->src_depth, &info, temp, hist)
	if(loud)
		srand((unsigned)__rdtsc());//

	CALC_LOSS();
	memcpy(invCR, info.invCR, sizeof(invCR));
	loss0=info.invCR[3];
	if(!niter)
		niter=C3_NPARAMS*10;
	for(int it=0, watchdog=0;it<niter;++it)
	{
		int idx[C3_OPT_NCOMP]={0}, inc[C3_OPT_NCOMP]={0}, stuck=0;
		if(loud_transforms&&loud)
			set_window_title("%d %4d/%4d,%d/%d: %.4lf%% RGB %.4lf %.4lf %.4lf%%",
				call_idx, it+1, niter, watchdog, shakethreshold,
				100.*info.invCR[3],
				100.*info.invCR[0],
				100.*info.invCR[1],
				100.*info.invCR[2]
			);
		if(watchdog>=shakethreshold)//bump if stuck
		{
			memcpy(info.params, srcparams, sizeof(info.params));
			for(int k=0;k<C3_NPARAMS;++k)
				info.params[k]+=((rand()&1)<<1)-1;
			watchdog=0;
			stuck=1;
		}
		else
		{
			for(int k=0;k<C3_OPT_NCOMP;++k)
			{
				idx[k]=rand()%C3_NPARAMS;
				while(!(inc[k]=(rand()&((1<<maskbits)-1))-(1<<(maskbits-1))));//reject zero delta
				info.params[idx[k]]+=(short)inc[k];
			}
		}

		CALC_LOSS();
		
		if(info.invCR[3]>invCR[3])//revert if worse
		{
			if(stuck)//a bad branch may surpass the local minimum
				memcpy(invCR, info.invCR, sizeof(info.invCR));
			else
			{
				memcpy(info.invCR, invCR, sizeof(info.invCR));
				for(int k=0;k<C3_OPT_NCOMP;++k)
					info.params[idx[k]]-=(short)inc[k];
			}
			++watchdog;
		}
		else//save if better
		{
			if(loss0>info.invCR[3])
			{
				memcpy(srcparams, info.params, sizeof(info.params));
				loss0=info.invCR[3];
				--it;//again
			}
			memcpy(invCR, info.invCR, sizeof(invCR));
			watchdog=0;
		}

		//preview
#if 1
		if(loud)
		{
			ch_entropy[0]=(float)(info.invCR[0]*src->src_depth[0]);
			ch_entropy[1]=(float)(info.invCR[1]*src->src_depth[1]);
			ch_entropy[2]=(float)(info.invCR[2]*src->src_depth[2]);
			//ch_entropy[3]=(float)(info.invCR[3]*(src->src_depth[0]+src->src_depth[1]+src->src_depth[2])/3);//X
			//int *ptr;
			//SWAPVAR(image, temp, ptr);
			io_render();
			//SWAPVAR(image, temp, ptr);
		}
#endif
	}
#undef  CALC_LOSS
	if(loss)
		memcpy(loss, invCR, sizeof(invCR));
	
	free(temp);
	free(hist);
}
#if 0
void custom3_opt_batch(Custom3Params *srcparams, int niter, int maskbits, int loud, double *loss)
{
	ArrayHandle folder=dialog_open_folder();
	if(!folder)
		return;
	const char *extensions[]=
	{
		"PNG",
		"JPG",
		"JPEG",
	};
	ArrayHandle filenames=get_filenames((char*)folder->data, extensions, _countof(extensions), 1);
	array_free(&folder);
	if(!filenames)
		return;

	int iw=0, ih=0;
	ArrayHandle bmp;
	ARRAY_ALLOC(int, bmp, 0, 0, 0, 0);
	for(int ks=0;ks<(int)filenames->count;++ks)
	{
		ArrayHandle *fn=(ArrayHandle*)array_at(&filenames, ks);
		int iw2=0, ih2=0;
		char *buf=(char*)stbi_load(fn[0]->data, &iw2, &ih2, 0, 4);
		if(!buf)
			continue;
		if(!iw||iw==iw2)
		{
			iw=iw2;
			ih+=ih2;
			ARRAY_APPEND(bmp, buf, iw2*ih2, 1, 0);
		}
		free(buf);
	}
	array_free(&filenames);
	if(bmp->count)
	{
		addhalf(bmp->data, iw, ih, 3, 4);
		colortransform_YCbCr_R_fwd((char*)bmp->data, iw, ih);//
		custom3_opt((char*)bmp->data, iw, ih, srcparams, niter, maskbits, loud, 0);
	}
	array_free(&bmp);
}
#endif


//CALIC - optimized for 8-bit	https://github.com/play-co/gcif
void pred_calic(Image *src, int fwd, int enable_ma)
{
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	//const int thresholds[]={5, 15, 25, 42, 60, 85, 140};
	int idx, pred;
	int *b2=(int*)malloc(res*sizeof(int[4]));
	int *arrN=(int*)malloc(2048*sizeof(int));//count
	int *arrS=(int*)malloc(2048*sizeof(int));//sum
	const int *pixels=fwd?src->data:b2, *errors=fwd?b2:src->data;
	if(!b2||!arrN||!arrS)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(b2, src->data, (size_t)res*sizeof(int[4]));//copy alpha
	for(int kc=0;kc<3;++kc)//process each channel separately
	{
		int depth=src->depth[kc], nlevels=1<<depth, half=nlevels>>1, sh=depth-8;
		memset(arrN, 0, 2048*sizeof(int));
		memset(arrS, 0, 2048*sizeof(int));
		for(int ky=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx)
			{
				int energy, pattern, context;
#define LOAD(BUF, X, Y) ((unsigned)(kx+(X))<(unsigned)src->iw&&(unsigned)(ky+(Y))<(unsigned)src->ih?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
				int
				//	NNWW=LOAD(pixels, -2, -2),
				//	NNW =LOAD(pixels, -1, -2),
					NN  =LOAD(pixels,  0, -2),
					NNE =LOAD(pixels,  1, -2),
				//	NNEE=LOAD(pixels,  2, -2),
				//	NWW =LOAD(pixels, -2, -1),
					NW  =LOAD(pixels, -1, -1),
					N   =LOAD(pixels,  0, -1),
					NE  =LOAD(pixels,  1, -1),
					WW  =LOAD(pixels, -2,  0),
					W   =LOAD(pixels, -1,  0),

					eN  =LOAD(errors,  0, -1),
					eW  =LOAD(errors, -1,  0);
#undef  LOAD
				//NNWW NNW NN NNE NNEE
				//NWW  NW  N  NE
				//WW   W   ?
				int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
				int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
				//int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
				//int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);

				//'CALIC - A context-based adaptive lossless image codec'	https://github.com/play-co/gcif/tree/master/refs/calic
#if 1
				int temp=dy-dx;
				temp=SHIFT_RIGHT_SIGNED(temp, sh);
				if(temp>80)
					pred=W;
				else if(temp<-80)
					pred=N;
				else
				{
					pred=(W+N)/2+(NE-NW)/4;
					//pred=((W+N)*3+(NE+NW))/8;
					//pred=((W+N)*5+(NE+NW)*3)/16;
					if(temp>32)
						pred=(pred+W)/2;
					else if(temp>8)
						pred=(pred*3+W)/4;
					else if(temp<-32)
						pred=(pred+N)/2;
					else if(temp<-8)
						pred=(pred*3+N)/4;
				}
#endif

				//clamped grad
#if 0
				int vmin, vmax;
				if(N<W)
					vmin=N, vmax=W;
				else
					vmin=W, vmax=N;
				pred=N+W-NW;
				pred=CLAMP(vmin, pred, vmax);
#endif

				//'A context-based adaptive lossless/nearly-lossless coding scheme for continuous-tone images'
#if 0
				if(dy+dx>32)//sharp edge
					pred=(dy*W+dx*N)/(dy+dx)+(NE-NW)/8;
				else if(dy-dx>12)//horizontal edge
					pred=(2*W+N)/3+(NE-NW)/8;
				else if(dy-dx<-12)//vertical edge
					pred=(W+2*N)/3+(NE-NW)/8;
				else//shallow area
					pred=(W+N)/2+(NE-NW)/8;

				if(d45-d135>32)//sharp 135-deg diagonal edge
					pred+=(NE-NW)/8;
				else if(d45-d135>16)//135-deg diagonal edge
					pred+=(NE-NW)/16;
				else if(d45-d135<-32)//sharp 45-deg diagonal edge
					pred-=(NE-NW)/8;
				else if(d45-d135<-16)//45-deg diagonal edge
					pred-=(NE-NW)/16;
#endif

				energy=dx+dy+abs(eN)+abs(eW);
				energy=SHIFT_RIGHT_SIGNED(energy, sh);
				if(energy<42)//quantization using binary search
				{
					if(energy<15)
						energy=energy<5?0:1;
					else
						energy=energy<25?2:3;
				}
				else
				{
					if(energy<85)
						energy=energy<60?4:5;
					else
						energy=energy<140?6:7;
				}

				pattern=(2*W-WW<pred)<<7|(2*N-NN<pred)<<6|(WW<pred)<<5|(NN<pred)<<4|(NE<pred)<<3|(NW<pred)<<2|(W<pred)<<1|(N<pred);

				context=pattern<<3|energy;

				if(arrN[context]>0)
					pred+=(arrS[context]+arrN[context]/2)/arrN[context];//sum of encountered errors / number of occurrences
				pred=CLAMP(-half, pred, half-1);

				idx=(src->iw*ky+kx)<<2|kc;
				pred^=-fwd;
				pred+=fwd;
				pred+=src->data[idx];
				if(enable_ma)
				{
					pred+=nlevels>>1;
					pred&=nlevels-1;
					pred-=nlevels>>1;
				}
				b2[idx]=pred;

				//update
				++arrN[context];
				arrS[context]+=errors[idx];
				if(arrN[context]>255)
				{
					arrN[context]>>=1;
					arrS[context]>>=1;
				}
			}
		}
	}
	memcpy(src->data, b2, res*sizeof(int[4]));
	free(b2);
	free(arrN);
	free(arrS);
}


#define G2_NPRED 21
#define G2_SSE_SIZE 1024
//#define G2_NPRED (21+12)
short g2_weights[]=
{
	//g	Y
	 0x00EA, 0x01C8, 0x00A2, 0x005E, 0x01F4, 0x0045, 0x0091, 0x0066, 0x003B, 0x0027,-0x0011, 0x001B, 0x00FF, 0x007E, 0x00D1, 0x00F3, 0x008F, 0x0130, 0x018E,-0x00AC,
	 0x0004,
	 //0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010,
	 0x010C, 0x0008,-0x007E, 0x00A2, 0x000E,-0x0069,-0x0073,-0x0125,-0x0092, 0x0000, 0x0078,

	//b	Cb
	 0x0006, 0x003D, 0x0031, 0x002F, 0x003F, 0x0015, 0x0011, 0x0036, 0x002E,-0x0022, 0x0011, 0x0034,-0x0007, 0x0012,-0x0018, 0x0012, 0x002F, 0x0000, 0x0000, 0x001C,
	 0x0004,
	 //0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010,
	 0x00A2, 0x02E1, 0x00C9,-0x00E0,-0x0068,-0x004E,-0x013E,-0x0012, 0x0001, 0x0000,-0x0046,
	 
	//r	Cr
	 0x003D, 0x0036, 0x0006, 0x007E, 0x0012, 0x0007, 0x0007, 0x0005, 0x001E, 0x0000, 0x0028, 0x0055,-0x0020, 0x0020, 0x0005, 0x0011, 0x0034, 0x0000, 0x0004, 0x003E,
	 0x0004,
	 //0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010,
	 -0x0100, 0x0001,-0x0086,-0x0041, 0x0051,-0x0080, 0x0004, 0x0002,-0x0003,-0x0003, 0x00D9,
};
typedef long long SseCell;
//typedef struct SseCellStruct
//{
//	int count, sum;//, sum7, sum8;
//} SseCell;
//SseCell g2_SSE[12][256];
//SseCell g2_SSE0[128], g2_SSE1[128], g2_SSE2[128], g2_SSE3[128];
//char g2_SSE_debug[768*512*4];
//static unsigned char g2_hash(unsigned x)//24 -> 7 bit
//{
//	unsigned short *p=(unsigned short*)&x;
//	unsigned char *p2=(unsigned char*)&x;
//	x=(p[0]^p[1])*0x45D9F3B;//https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
//	x=(p[0]^p[1])*0x45D9F3B;
//	x=p[0]^p[1];
//	x=p2[0]^p2[1];
//	return x;
//}
void pred_grad2(Image *src, int fwd, int enable_ma)
{
	//const int *guide=im0->data;

	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int *b2=(int*)malloc((size_t)res*sizeof(int[4]));
	int *perrors=(int*)malloc((size_t)src->iw*(G2_NPRED+1)*2*sizeof(int));
	SseCell *sse=(SseCell*)malloc(sizeof(SseCell[G2_SSE_SIZE]));
	//int *SSE_count=(int*)malloc(256*sizeof(int));
	//int *SSE_sum=(int*)malloc(256*sizeof(int));
	const int *pixels=fwd?src->data:b2, *errors=fwd?b2:src->data;
	if(!b2||!perrors||!sse)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(b2, src->data, res*sizeof(int[4]));//copy alpha
	for(int kc=0;kc<4;++kc)
	{
		int depth;
		short *params;
		int *hireserror;
		//int maxerror=0;

		depth=src->depth[kc];
		if(!depth)
			continue;
		params=g2_weights+(_countof(g2_weights)/3)*kc;
		hireserror=perrors+src->iw*2*G2_NPRED;
		memset(perrors, 0, 2LL*src->iw*(G2_NPRED+1)*sizeof(int));

		memset(sse, 0, sizeof(SseCell[G2_SSE_SIZE]));
		//memset(g2_SSE0, 0, sizeof(g2_SSE0));
		//memset(g2_SSE1, 0, sizeof(g2_SSE1));
		//memset(g2_SSE2, 0, sizeof(g2_SSE2));
		//memset(g2_SSE3, 0, sizeof(g2_SSE3));

		//memset(SSE_count, 0, 256*sizeof(int));
		//memset(SSE_sum, 0, 256*sizeof(int));

		//int weights[3]=
		//{
		//	0x8000,
		//	0x8000,
		//	0x8000,
		//};
		//int certainty=512;
		for(int ky=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx)
			{
#define LOAD(BUF, X, Y) (unsigned)(kx+(X))<(unsigned)src->iw&&(unsigned)(ky+(Y))<(unsigned)src->ih?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]<<8:0
				int
					NNNNNN  =LOAD(pixels,  0, -6),
				//	NNNNWWWW=LOAD(pixels, -4, -4),
				//	NNNN    =LOAD(pixels,  0, -4),
				//	NNNNEEEE=LOAD(pixels,  4, -4),
				//	NNNWWW  =LOAD(pixels, -3, -3),
					NNN     =LOAD(pixels,  0, -3),
				//	NNNEEE  =LOAD(pixels,  3, -3),
					NNWW    =LOAD(pixels, -2, -2),
					NNW     =LOAD(pixels, -1, -2),
					NN      =LOAD(pixels,  0, -2),
					NNE     =LOAD(pixels,  1, -2),
					NNEE    =LOAD(pixels,  2, -2),
					NWW     =LOAD(pixels, -2, -1),
					NW      =LOAD(pixels, -1, -1),
					N       =LOAD(pixels,  0, -1),
					NE      =LOAD(pixels,  1, -1),
					NEE     =LOAD(pixels,  2, -1),
					NEEE    =LOAD(pixels,  3, -1),
					NEEEE   =LOAD(pixels,  4, -1),
					NEEEEE  =LOAD(pixels,  5, -1),
					NEEEEEE =LOAD(pixels,  6, -1),
					NEEEEEEE=LOAD(pixels,  7, -1),
					WWWWWW  =LOAD(pixels, -6,  0),
					WWWW    =LOAD(pixels, -4,  0),
				//	WWW     =LOAD(pixels, -3,  0),
					WW      =LOAD(pixels, -2,  0),
					W       =LOAD(pixels, -1,  0);
#undef  LOAD
#define LOAD(BUF, X, Y) (unsigned)(kx+(X))<(unsigned)src->iw&&(unsigned)(ky+(Y))<(unsigned)src->ih?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0
				int
				//	dNNNWWW=LOAD(errors, -3, -3),
				//	dNNNWW =LOAD(errors, -2, -3),
				//	dNNNW  =LOAD(errors, -1, -3),
				//	dNNN   =LOAD(errors,  0, -3),
				//	dNNNE  =LOAD(errors,  1, -3),
				//	dNNNEE =LOAD(errors,  2, -3),
				//	dNNNEEE=LOAD(errors,  3, -3),
				//	dNNWWW =LOAD(errors, -3, -2),
					dNNWW  =LOAD(errors, -2, -2),
					dNNW   =LOAD(errors, -1, -2),
					dNN    =LOAD(errors,  0, -2),
					dNNE   =LOAD(errors,  1, -2),
					dNNEE  =LOAD(errors,  2, -2),
				//	dNNEEE =LOAD(errors,  3, -2),
				//	dNWWW  =LOAD(errors, -3, -1),
					dNWW   =LOAD(errors, -2, -1),
					dNW    =LOAD(errors, -1, -1),
					dN     =LOAD(errors,  0, -1),
					dNE    =LOAD(errors,  1, -1),
					dNEE   =LOAD(errors,  2, -1),
				//	dNEEE  =LOAD(errors,  3, -1),
				//	dWWW   =LOAD(errors, -3,  0),
					dWW    =LOAD(errors, -2,  0),
					dW     =LOAD(errors, -1,  0);
#undef  LOAD
#define LOAD(X, Y) (unsigned)(kx+(X))<(unsigned)src->iw&&(unsigned)(ky+(Y))<(unsigned)src->ih?hireserror[src->iw*((ky+(Y))&1)+kx+(X)]:0
				int
					eNW=LOAD(-1, -1),
					eN =LOAD( 0, -1),
					eNE=LOAD( 1, -1),
					eW =LOAD(-1,  0);
#undef  LOAD
				int preds[G2_NPRED], pred;
				int prederrs[G2_NPRED];
				int paper_GAP, calic_GAP;
				int j;
				long long num;
				int weights[G2_NPRED], den;
				int idx;

				for(int k=0;k<G2_NPRED;++k)
				{
					//eNW + eN + eNE
					prederrs[k]=
						((unsigned)(kx-1)<(unsigned)src->iw?perrors[src->iw*2*k+src->iw*((ky-1)&1)+kx-1]:0)+
						perrors[src->iw*2*k+src->iw*((ky-1)&1)+kx]+
						((unsigned)(kx+1)<(unsigned)src->iw?perrors[src->iw*2*k+src->iw*((ky-1)&1)+kx+1]:0);
				}
				//int clamp_min=N, clamp_max=N;
				//UPDATE_MIN(clamp_min, W);
				//UPDATE_MAX(clamp_max, W);
				//UPDATE_MIN(clamp_min, NE);
				//UPDATE_MAX(clamp_max, NE);
				{
					int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
					int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
					int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
					int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
					int diff=(dy-dx)>>8, diff2=(d45-d135)>>8, diff3=NE-NW;
					if(dy+dx>32)
						paper_GAP=(int)(((long long)dx*N+(long long)dy*W)/((long long)dy+dx));
					else if(diff>12)
						paper_GAP=(N+2*W)/3;
					else if(diff<-12)
						paper_GAP=(2*N+W)/3;
					else
						paper_GAP=(N+W)>>1;

					if(diff2>32)
						paper_GAP+=diff3>>2;
					else if(diff2>16)
						paper_GAP+=diff3*3>>4;
					else if(diff2>=-16)
						paper_GAP+=diff3>>3;
					else if(diff2>=-32)
						paper_GAP+=diff3>>4;

					if(diff>80)
						calic_GAP=W;
					else if(diff<-80)
						calic_GAP=N;
					else if(diff>32)
						calic_GAP=(2*N+6*W+NE-NW)>>3;		//c1	[1/4  3/4  1/8  -1/8].[N W NE NW]
					else if(diff>8)
						calic_GAP=(6*N+10*W+3*(NE-NW))>>4;	//c2	[3/8  5/8  3/16  -3/16]
					else if(diff<-32)
						calic_GAP=(6*N+2*W+NE-NW)>>3;		//c3	[3/4  1/4  1/8  -1/8]
					else if(diff<-8)
						calic_GAP=(10*N+6*W+3*(NE-NW))>>4;	//c4	[5/8  3/8  3/16  -3/16]
					else
						calic_GAP=(((N+W)<<1)+NE-NW)>>2;	//c5	[1/2  1/2  1/4  -1/4]
				}
				j=-1;
				//the 4 predictors from JPEG XL:
				++j, preds[j]=W+NE-N;
				++j, preds[j]=N-((eN+eW+eNE)*params[G2_NPRED+6]>>8);
				++j, preds[j]=W-((eN+eW+eNW)*params[G2_NPRED+5]>>8);
				++j, preds[j]=N-((eNW*params[G2_NPRED]+eN*params[G2_NPRED+1]+eNE*params[G2_NPRED+2]+(NN-N)*params[G2_NPRED+3]+(NW-W)*params[G2_NPRED+4])>>8);
				//++j, preds[j]=N+((eN+eW+eNE)*10>>5);
				//++j, preds[j]=W+((eN+eW+eNW)*10>>5);
				//++j, preds[j]=N+((eNW*5+eN*5+eNE*5+(NN-N)*12+(NW-W)*4)>>5);
				
				++j, preds[j]=W -(eW *params[G2_NPRED+ 7]>>8);
				++j, preds[j]=N -(eN *params[G2_NPRED+ 8]>>8);
				++j, preds[j]=NW-(eNW*params[G2_NPRED+ 9]>>8);
				++j, preds[j]=NE-(eNE*params[G2_NPRED+10]>>8);
				//++j, preds[j]=W+(eW*10>>5);
				//++j, preds[j]=N+(eN*10>>5);
				//++j, preds[j]=NW+(eNW*8>>5);
				//++j, preds[j]=NE+(eNE*8>>5);
				//++j, preds[j]=W+(prederrs[j]*10>>5);
				//++j, preds[j]=N+(prederrs[j]*10>>5);
				//++j, preds[j]=NW+(prederrs[j]*8>>5);
				//++j, preds[j]=NE+(prederrs[j]*8>>5);
				//++j, preds[j]=W+correction;
				//++j, preds[j]=N+correction;
				//++j, preds[j]=NW+correction;
				//++j, preds[j]=NE+correction;
				
				++j, preds[j]=calic_GAP;
				//++j, preds[j]=(N+W)>>1;
				++j, preds[j]=(4*(N+W+NW+NE)-(NN+WW+NNWW+NNEE)+6)/12;
				++j, preds[j]=(3*W+NEEE)>>2;
				++j, preds[j]=paper_GAP;
				//++j, preds[j]=N+NE-NNE;
				//++j, preds[j]=clamp4(N+W -NW , N, W, NW, NE);
				//++j, preds[j]=clamp4(W+NE-N  , N, W, NW, NE);
				//++j, preds[j]=clamp4(N+NW-NNW, N, W, NW, NE);
				//++j, preds[j]=clamp4(N+NE-NNE, N, W, NE, NEE);
				
				++j, preds[j]=(W+NEE)>>1;		//this predicts E0.5
				++j, preds[j]=NNNNNN;			//this predicts N6
				++j, preds[j]=(NEEEE+NEEEEEE)>>1;	//this predicts NE5
				++j, preds[j]=(WWWW+WWWWWW)>>1;		//this predicts W5
				++j, preds[j]=(N+W+NEEEEE+NEEEEEEE)>>2;	//this predicts NE2.75

				//++j, preds[j]=W*3-WW*3+WWW;//parabolic (3rd derivative)
				//++j, preds[j]=N*3-NN*3+NNN;
				//++j, preds[j]=NW*3-NNWW*3+NNNWWW;
				//++j, preds[j]=NE*3-NNEE*3+NNNEEE;
				//++j, preds[j]=W*2-WW;//linear
				//++j, preds[j]=N*2-NN;
				//++j, preds[j]=NW*2-NNWW;
				//++j, preds[j]=NE*2-NNEE;
				//++j, preds[j]=4*W-6*WW+4*WWW-WWWW, preds[j]=CLAMP(vmin, preds[j], vmax);//cubic
				//++j, preds[j]=4*N-6*NN+4*NNN-NNNN, preds[j]=CLAMP(vmin, preds[j], vmax);
				//++j, preds[j]=4*NW-6*NNWW+4*NNNWWW-NNNNWWWW, preds[j]=CLAMP(vmin, preds[j], vmax);
				//++j, preds[j]=4*NE-6*NNEE+4*NNNEEE-NNNNEEEE, preds[j]=CLAMP(vmin, preds[j], vmax);

				++j, preds[j]=3*(N-NN)+NNN;
				//++j, preds[j]=clamp4(N*2-NN, N, W, NE, NEE);
				++j, preds[j]=(N+NNN)>>1;		//this predicts NN
				++j, preds[j]=((N+W)*3-NW*2)>>2;

				++j, preds[j]=N+W-NW;//, preds[j]=MEDIAN3(N, W, preds[j]);
				//++j; GRAD(preds[j], N, W, NW, vmin, vmax);
				//++j, preds[j]=W+NE-N, preds[j]=CLAMP(vmin, preds[j], vmax);

				//++j;
				//int gx1=N-NW, gx2=W-WW, gy1=W-NW, gy2=N-NN;
				//if(abs(gx1-gx2)<abs(gy1-gy2))
				//	preds[j]=W+(gx1*11+gx2*5)/16;
				//else
				//	preds[j]=N+(gy1*11+gy2*5)/16;

				//++j, preds[j]=(W+((N-NW)*11+(W-WW)*5)/16 + N+((W-NW)*11+(N-NN)*5)/16)/2;
				//++j, preds[j]=NW+(N-NW)+(W-NW);//==N+W-NW
				//++j, preds[j]=((N+W)*5+(NE+NW)*3)/16;
				//++j, preds[j]=(N+W)/2+(NE-NW)/4;
				//++j, preds[j]=N;
				//++j, preds[j]=W;
				//++j; GRAD(preds[j], NE, NW, NN, vmin, vmax);
				//++j; GRAD(preds[j], NN, WW, NNWW, vmin, vmax);

#if 0
				++j, preds[j]=2*W -WW  , preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);//linear
				++j, preds[j]=2*N -NN  , preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
				++j, preds[j]=2*NW-NNWW, preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
				++j, preds[j]=2*NE-NNEE, preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);

				++j, preds[j]=3*W -3*WW  +WWW   , preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);//parabolic
				++j, preds[j]=3*N -3*NN  +NNN   , preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
				++j, preds[j]=3*NW-3*NNWW+NNNWWW, preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
				++j, preds[j]=3*NE-3*NNEE+NNNEEE, preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);

				++j, preds[j]=4*W -6*WW  -4*WWW   +WWWW    , preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);//cubic
				++j, preds[j]=4*N -6*NN  -4*NNN   +NNNN    , preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
				++j, preds[j]=4*NW-6*NNWW-4*NNNWWW+NNNNWWWW, preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
				++j, preds[j]=4*NE-6*NNEE-4*NNNEEE+NNNNEEEE, preds[j]=CLAMP(-(128<<8), preds[j], 127<<8);
#endif
				
				num=0;
				den=0;
				for(int k=0;k<G2_NPRED;++k)
				{
					weights[k]=(params[k]<<8)/(prederrs[k]+1);
					//if(maxerror<prederrs[k])
					//	maxerror=prederrs[k];
					//weights[k]=maxerror-prederrs[k]+1;

					num+=(long long)preds[k]*weights[k];
					den+=weights[k];
				}
				pred=den?(int)(num/den):preds[0];
				
				//pred=CLAMP(vmin, pred, vmax);
				
				
				//int pred0=pred;
				//int energy=abs(eN)+abs(eW)+abs(eNW)+abs(eNE);
				//pred-=pred*energy/certainty;//pred=pred0*(1-energy/certainty)

				//pred-=pred*(abs(eN)+abs(eW)+abs(eNW))/1108;//improves only with kodim13

				//pred*=(1108-(abs(eN)+abs(eW)+abs(eNW)));
				//pred/=1108;

				//pred/=(abs(eN)+abs(eW)+abs(eNW))/128+1;

				//if(abs(eN)+abs(eW)+abs(eNW)>512)pred/=2;

				//pred+=(eN>0&&eW>0&&eNW>0)-(eN<0&&eW<0&&eNW<0);

				//if(eN>0&&eW>0&&eNW>0)pred=MAXVAR(vmax, NW);
				//else if(eN<0&&eW<0&&eNW<0)pred=MINVAR(vmin, NW);

				//pred=CLAMP(-128, pred, 127);

				idx=(src->iw*ky+kx)<<2|kc;
				{
					int val, curr, delta;
				//Secondary Symbol Estimation (SSE)
#if 1
				int dnb[]=//error neighbors
				{
					dN   ,//order matters slightly
					dW   ,
					dNW  ,
					dNE  ,
					dNN  ,
					dWW  ,
					dNNWW,
					dNNW ,
					dNNE ,
					dNNEE,
					dNWW ,
					dNEE ,
					
					//dNNWW,
					//dNNW ,
					//dNN  ,
					//dNNE ,
					//dNNEE,
					//dNWW ,
					//dNW  ,
					//dN   ,
					//dNE  ,
					//dNEE ,
					//dWW  ,
					//dW   ,
					
					//dN     ,
					//dW     ,
					//dNW    ,
					//dNE    ,
					//dNN    ,
					//dWW    ,
					//dNNW   ,
					//dNNE   ,
					//dNWW   ,
					//dNEE   ,
					//dNNN   ,
					//dWWW   ,
					//dNNNW  ,
					//dNNNE  ,
					//dNNWW  ,
					//dNNEE  ,
					//dNEEE  ,
					//dNWWW  ,
					//dNNNWW ,
					//dNNNEE ,
					//dNNWWW ,
					//dNNEEE ,
					//dNNNWWW,
					//dNNNEEE,
					
					//dNNNWWW,
					//dNNNWW ,
					//dNNNW  ,
					//dNNN   ,
					//dNNNE  ,
					//dNNNEE ,
					//dNNNEEE,
					//dNNWWW ,
					//dNNWW  ,
					//dNNW   ,
					//dNN    ,
					//dNNE   ,
					//dNNEE  ,
					//dNNEEE ,
					//dNWWW  ,
					//dNWW   ,
					//dNW    ,
					//dN     ,
					//dNE    ,
					//dNEE   ,
					//dNEEE  ,
					//dWWW   ,
					//dWW    ,
					//dW     ,
					
					//eNW,
					//eN ,
					//eNE,
					//eW ,
				};
				SseCell *pc[12];
				int hashval;

				for(int k2=3;k2>=0;--k2)
				{
					//int *dnbk=dnb+((size_t)k2<<2);//X
					//hashval=pred>>(depth+8-4)<<8;
					//int temp=dnbk[0]+dnbk[1]+dnbk[2]+dnbk[3];
					//hashval+=floor_log2_32(abs(temp))<<1|(temp<0);
					//hashval&=0xFF;
					//pc[k2]=sse+((size_t)k2<<8)+hashval;

					hashval=(pred+128)>>(8+6)&3;
					for(int k=0;k<3;++k)
					{
						hashval<<=3;
						hashval+=(int)(dnb[3*k2+k]*0xFC28>>(16+6)&3);
					}
					hashval&=0xFF;
					pc[k2]=sse+((size_t)k2<<8)+hashval;
					//pc[k2]=g2_SSE[k2]+hashval;

					//int c=2, d=pc[k2]->sum8+pc[k2]->sum7;
					//if(pc[k2]->count)
					//{
					//	++c;
					//	d+=(int)((((long long)pc[k2]->sum<<8)+(pc[k2]->count>>1))/pc[k2]->count);
					//}
					//pred+=d/c;
					{
						long long sum=pc[k2][0]>>12;
						int count=pc[k2][0]&0xFFF;
						pred+=(int)((sum+(count>>1))/(count+1LL));
						//pred=CLAMP(-(nlevels<<7), pred, (nlevels<<7)-1);
					}
				}
				CLAMP3_32(pred, pred, N, W, NE);
				//pred=CLAMP(clamp_min, pred, clamp_max);

				val=(pred+128)>>8;
				val^=-fwd;
				val+=fwd;
				val+=src->data[idx];
				if(enable_ma)
				{
					val<<=32-depth;
					val>>=32-depth;
				}
				b2[idx]=val;

				curr=pixels[idx]<<8;
				delta=curr-pred;
				for(int k2=0;k2<4;++k2)
				{
					long long sum=pc[k2][0]>>12;
					int count=pc[k2][0]&0xFFF;
					if(count+1>640)
					{
						count>>=1;
						sum>>=1;
					}
					++count;
					sum+=delta;
					pc[k2][0]=sum<<12|count;
					//pc[k]->sum8=(pc[k]->sum8*255+(delta<<8))>>8;
					//pc[k]->sum7=(pc[k]->sum7*127+(delta<<8))>>7;
				}
#endif
#if 0
				//if(kc==0&&kx==255&&ky==1)//
				//	printf("");

				SseCell *pc[3];
				pc[0]=g2_SSE[0]+((dNW>>7&1)<<6|(dNE>>7&1)<<5|(dN>>7&1)<<4|(dW>>7&1)<<3|pred>>(5+8)&7);
				pred+=pc[0]->count?(int)((((long long)pc[0]->sum<<8)+(pc[0]->count>>1))/pc[0]->count):0;
				pred=CLAMP(-(128<<8), pred, 127<<8);
				
				pc[1]=g2_SSE[1]+((dNNWW>>7&1)<<4|(dNNW>>7&1)<<3|(dNWW>>7&1)<<2|(dWW>>7&1)<<1|pred>>(7+8)&1);
				pred+=pc[1]->count?(int)((((long long)pc[1]->sum<<8)+(pc[1]->count>>1))/pc[1]->count):0;
				pred=CLAMP(-(128<<8), pred, 127<<8);
				
				pc[2]=g2_SSE[2]+((dNN>>7&1)<<4|(dNNE>>7&1)<<3|(dNNEE>>7&1)<<2|(dNEE>>7&1)<<1|pred>>(7+8)&1);
				pred+=pc[2]->count?(pc[2]->sum+(pc[2]->count>>1))/pc[2]->count:0;
				pred=CLAMP(-(128<<8), pred, 127<<8);

				int delta;
				if(fwd)
				{
					delta=buf[idx]-((pred+128)>>8);
					b2[idx]=delta;
				}
				else
				{
					delta=buf[idx];
					b2[idx]=delta+((pred+128)>>8);
				}
				for(int k=0;k<3;++k)
				{
					if(pc[k]->count+1>640)
					{
						pc[k]->count>>=1;
						pc[k]->sum>>=1;
					}
					++pc[k]->count;
					pc[k]->sum+=delta;
				}
#endif
#if 0
				SseCell *pc[3];
				pc[0]=g2_SSE0+((dNW>>7&1)<<6|(dNE>>7&1)<<5|(dN>>7&1)<<4|(dW>>7&1)<<3|pred>>(5+8)&7);
				pred+=pc[0]->count?(int)((((long long)pc[0]->sum<<8)+(pc[0]->count>>1))/pc[0]->count):0;
				pred=CLAMP(-(128<<8), pred, 127<<8);
				
				pc[1]=g2_SSE1+((dNNWW>>7&1)<<4|(dNNW>>7&1)<<3|(dNWW>>7&1)<<2|(dWW>>7&1)<<1|pred>>(7+8)&1);
				pred+=pc[1]->count?(int)((((long long)pc[1]->sum<<8)+(pc[1]->count>>1))/pc[1]->count):0;
				pred=CLAMP(-(128<<8), pred, 127<<8);
				
				//pc[2]=g2_SSE2+((dNN>>7&1)<<7|(dNNE>>7&1)<<6|(dNNEE>>7&1)<<5|(dNEE>>7&1)<<4|ipred>>4&15);
				//ipred+=pc[2]->count?(pc[2]->sum+(pc[2]->count>>1))/pc[2]->count:0;
				//ipred=CLAMP(-128, ipred, 127);

				int delta;
				if(fwd)
				{
					delta=buf[idx]-((pred+128)>>8);
					b2[idx]=delta;
				}
				else
				{
					delta=buf[idx];
					b2[idx]=delta+((pred+128)>>8);
				}
				//if(fwd)//
				//	g2_SSE_debug[idx]=pixels[idx];
				//else if((pixels[idx]&0xFF)!=(g2_SSE_debug[idx]&0xFF))//
				//	LOG_ERROR("SSE error");
				for(int k=0;k<2;++k)
				{
					if(pc[k]->count+1>640)
					{
						pc[k]->count>>=1;
						pc[k]->sum>>=1;
					}
					++pc[k]->count;
					pc[k]->sum+=delta;
				}
#endif
#if 0
				int ctx=(dNW>>6&3)<<6|(dN>>6&3)<<4|(dNE>>6&3)<<2|(dW>>6&3);
				//int ctx=(dNNWW>>7&1)<<7|(dNNW>>7&1)<<6|(dNN>>7&1)<<5|(dNNE>>7&1)<<4|(dNNEE>>7&1)<<3|(dNWW>>7&1)<<2|(dNEE>>7&1)<<1|(dWW>>7&1);
				//int ctx=
				//	(dNNWW>>7&1)<<11|
				//	(dNNW >>7&1)<<10|
				//	(dNN  >>7&1)<<9|
				//	(dNNE >>7&1)<<8|
				//	(dNNEE>>7&1)<<7|
				//	(dNWW >>7&1)<<6|
				//	(dNW  >>7&1)<<5|
				//	(dN   >>7&1)<<4|
				//	(dNE  >>7&1)<<3|
				//	(dNEE >>7&1)<<2|
				//	(dWW  >>7&1)<<1|
				//	(dW   >>7&1);
				SseCell *pc=g2_SSE0+ctx;
				int SSE_bias=pc->count?(pc->sum+(pc->count>>1))/pc->count:0;
				ipred+=SSE_bias;
				ipred=CLAMP(-128, ipred, 127);
				int delta;
				if(fwd)
				{
					delta=buf[idx]-ipred;
					b2[idx]=delta;
				}
				else
				{
					delta=buf[idx];
					b2[idx]=delta+ipred;
				}
				if(pc->count+1>256)
				{
					pc->count>>=1;
					pc->sum>>=1;
				}
				++pc->count;
				pc->sum+=delta;
#endif
#if 0
				SseCell *pc[4];
				int ctx, SSE_b[4];

				ipred=0;

				pc[0]=g2_SSE0+(((dW&0xF8)<<2|ipred>>3)&0xFF);
				SSE_b[0]=pc[0]->count?(pc[0]->sum+(pc[0]->count>>1))/pc[0]->count:0;

				pc[1]=g2_SSE1+(((dN&0xF8)<<2|SSE_b[0]>>3)&0xFF);
				SSE_b[1]=pc[1]->count?(pc[1]->sum+(pc[1]->count>>1))/pc[1]->count:0;

				pc[2]=g2_SSE2+(((dNW&0xF8)<<2|SSE_b[1]>>3)&0xFF);
				SSE_b[2]=pc[2]->count?(pc[2]->sum+(pc[2]->count>>1))/pc[2]->count:0;

				pc[3]=g2_SSE3+(((dNE&0xF8)<<2|SSE_b[2]>>3)&0xFF);
				SSE_b[3]=pc[3]->count?(pc[3]->sum+(pc[3]->count>>1))/pc[3]->count:0;

				ipred+=SSE_b[3];
				ipred=CLAMP(-128, ipred, 127);
				int delta;
				if(fwd)
				{
					delta=buf[idx]-ipred;
					b2[idx]=delta;
				}
				else
				{
					delta=buf[idx];
					b2[idx]=delta+ipred;
				}
				for(int k=0;k<4;++k)
				{
					if(pc[k]->count+1>256)
					{
						pc[k]->count>>=1;
						pc[k]->sum>>=1;
					}
					++pc[k]->count;
					pc[k]->sum+=delta;
				}
#endif
#if 0
				//if(kc==0&&kx==553&&ky==4)//
				//if(kc==1&&kx==684&&ky==7)//
				//	printf("");

				//int ctx=((dNW>>6&3)<<6|(dN>>6&3)<<4|(dNE>>6&3)<<2|(dW>>6&3))+(kx<<4)/iw+(ky<<4)/ih;
				int ctx=((dNW>>6&3)<<6|(dN>>6&3)<<4|(dNE>>6&3)<<2|(dW>>6&3));
				int SSE_bias=SSE_count[ctx]?(SSE_sum[ctx]+(SSE_count[ctx]>>1))/SSE_count[ctx]:0;

				//if(SSE_bias)//
				//	printf("");

				ipred+=SSE_bias;
				ipred=CLAMP(-128, ipred, 127);
				
				int delta;//modular arithmetic
				if(fwd)
				{
					delta=buf[idx]-ipred;
					b2[idx]=delta;
				}
				else
				{
					delta=buf[idx];
					b2[idx]=delta+ipred;
				}

				//if(fwd)//
				//	g2_SSE_debug[idx]=SSE_bias;
				//else if((SSE_bias&0xFF)!=(g2_SSE_debug[idx]&0xFF))//
				//	LOG_ERROR("SSE error");

				if(SSE_count[ctx]+1>256)
				{
					SSE_count[ctx]>>=1;
					SSE_sum[ctx]>>=1;
				}
				++SSE_count[ctx];
				SSE_sum[ctx]+=delta;
#endif

				//no SSE
#if 0
				int spred=pred;
				spred^=-fwd;
				spred+=fwd;
				b2[idx]=src->data[idx]+((spred+128)>>8);
#endif

				//int curr=pixels[idx]<<8;
				hireserror[src->iw*(ky&1)+kx]=curr-pred;
				for(int k=0;k<G2_NPRED;++k)
				{
					int e=abs(curr-preds[k]);
					perrors[src->iw*2*k+src->iw*(ky&1)+kx]=e;
					if(kx<src->iw&&ky>0)
						perrors[src->iw*2*k+src->iw*((ky-1)&1)+kx+1]+=e;
				}
#if 0
				int curr=pixels[idx];
				//curr=pred0*(1-energy/temperature)	->	temperature = energy/(1 - abs(curr/pred0))

				if(abs(curr)<abs(pred))
				{
					certainty-=10;
					if(certainty<512)
						certainty=512;
				}
				else
					++certainty;

				//certainty+=abs(curr)<abs(pred)?-1:1;

				//int den=pred0?0x10000-abs((curr<<16)/pred0):1;
				//if(den>0)
				//{
				//	int t2=(energy<<16)/den;
				//	if(t2>512)
				//		temperature=t2;
				//}

				//int num=pred0*energy, den=pred0-curr;
				//if(num&&den&&abs(num)>abs(den))
				//	temperature=abs(num/den);
#endif
				}
			}
		}
	}
	memcpy(src->data, b2, res*sizeof(int[4]));
	free(perrors);
	free(sse);
	free(b2);
	//free(SSE_count);
	//free(SSE_sum);
}




//DWTs

ArrayHandle dwt2d_gensizes(int iw, int ih, int wstop, int hstop, int nstages_override)//calculate dimensions of each DWT stage in descending order
{
	int nstages;
	DWTSize *p;
	ArrayHandle sizes;

	int lw=FLOOR_LOG2(iw), lh=FLOOR_LOG2(ih), lmin=MINVAR(lw, lh);
	ARRAY_ALLOC(DWTSize, sizes, 0, 0, lmin, 0);
	if(wstop<3)
		wstop=3;
	if(hstop<3)
		hstop=3;
	nstages=0;
	p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);
	p->w=(unsigned short)iw;
	p->h=(unsigned short)ih;
	for(int w2=iw, h2=ih;w2>=wstop&&h2>=hstop&&(!nstages_override||nstages<nstages_override);++nstages)
	{
		p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);

		w2>>=1;//w=floor(w/2)
		h2>>=1;//h=floor(h/2)
		//w2-=w2>>1;//w=ceil(w/2)
		//h2-=h2>>1;//h=ceil(h/2)

		p->w=(unsigned short)w2;
		p->h=(unsigned short)h2;
	}
	return sizes;
}

#if 0
void dwt2d_grad_fwd(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	for(int it=sizes_start;it<sizes_end-1;++it)
	//for(int it=sizes_start;it<1;++it)
	{
		for(int ky=0;ky<sizes[it].h-1;ky+=2)
		{
			for(int kx=0;kx<sizes[it].w-1;kx+=2)
			{
				int idx=sizes->w*ky+kx;
				int v[]=
				{
					buffer[ idx            *stride],
					buffer[(idx         +1)*stride],
					buffer[(idx+sizes->w  )*stride],
					buffer[(idx+sizes->w+1)*stride],
				};

				//if((unsigned char)v[0]==0xFF||(unsigned char)v[1]==0xFF||(unsigned char)v[2]==0xFF||(unsigned char)v[3]==0xFF)
				//	v[0]=0xFF;

				char vmin, vmax;
				if(v[1]<v[2])
					vmin=v[1], vmax=v[2];
				else
					vmin=v[2], vmax=v[1];
				if(v[0]<vmin)
					v[3]-=vmax;
				else if(v[0]>vmax)
					v[3]-=vmin;
				else
					v[3]-=v[1]+v[2]-v[0];

				v[2]-=v[0];
				v[1]-=v[0];
				v[0]+=(v[1]+v[2])>>2;

				buffer[ idx            *stride]=v[3];//grad
				buffer[(idx         +1)*stride]=v[2];//diffy
				buffer[(idx+sizes->w  )*stride]=v[1];//diffx
				buffer[(idx+sizes->w+1)*stride]=v[0];//av

				//buffer[ idx            <<2|kc]=v[0];//av
				//buffer[(idx         +1)<<2|kc]=v[1];//diffx
				//buffer[(idx+sizes->w  )<<2|kc]=v[2];//diffy
				//buffer[(idx+sizes->w+1)<<2|kc]=v[3];//grad
			}
		}
		dwt2d_lazy_fwd(buffer, sizes, it, it+2, stride, temp);
	}
}
void dwt2d_grad_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		dwt2d_lazy_inv(buffer, sizes, it, it+2, stride, temp);
		for(int ky=0;ky<sizes[it].h-1;ky+=2)
		{
			for(int kx=0;kx<sizes[it].w-1;kx+=2)
			{
				int idx=sizes->w*ky+kx;
				char v[]=
				{
					buffer[ idx            *stride],
					buffer[(idx         +1)*stride],
					buffer[(idx+sizes->w  )*stride],
					buffer[(idx+sizes->w+1)*stride],
				};

				v[0]-=(v[1]+v[2])>>2;
					
				v[1]+=v[0];
				v[2]+=v[0];

				char vmin, vmax;
				if(v[1]<v[2])
					vmin=v[1], vmax=v[2];
				else
					vmin=v[2], vmax=v[1];
				if(v[0]<vmin)
					v[3]+=vmax;
				else if(v[0]>vmax)
					v[3]+=vmin;
				else
					v[3]+=v[1]+v[2]-v[0];

				buffer[ idx            *stride]=v[3];
				buffer[(idx         +1)*stride]=v[2];
				buffer[(idx+sizes->w  )*stride]=v[1];
				buffer[(idx+sizes->w+1)*stride]=v[0];
			}
		}
	}
}

double customdwtparams[12]={0};
void dwt1d_custom_fwd(char *buffer, int count, int stride, char *b2, const double *params)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	even[0]-=(char)floor(odd[0]*(params[0]+params[5]));
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(char)floor(odd[k-1]*params[0]+odd[k]*params[5]);
	if(extraeven)
		even[nodd]-=(char)floor(odd[nodd-1]*(params[0]+params[5]));
	
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(char)floor(even[k]*params[1]+even[k+1]*params[6]);
	if(!extraeven)
		odd[nodd-1]+=(char)floor(even[nodd-1]*(params[1]+params[6]));


	even[0]-=(char)(odd[0]*(params[2]+params[7]));
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(char)floor(odd[k-1]*params[2]+odd[k]*params[7]);
	if(extraeven)
		even[nodd]-=(char)floor(odd[nodd-1]*(params[2]+params[7]));
	
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(char)floor(even[k]*params[3]+even[k+1]*params[8]);
	if(!extraeven)
		odd[nodd-1]+=(char)floor(even[nodd-1]*(params[3]+params[8]));


	even[0]-=(char)floor(odd[0]*(params[4]+params[9]));
	for(int k=1;k<nodd;++k)//predict
		even[k]-=(char)floor(odd[k-1]*params[4]+odd[k]*params[9]);
	if(extraeven)
		even[nodd]-=(char)floor(odd[nodd-1]*(params[4]+params[9]));


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_custom_inv(char *buffer, int count, int stride, char *b2, const double *params)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	even[0]+=(char)floor(odd[0]*(params[4]+params[9]));
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(char)floor(odd[k-1]*params[4]+odd[k]*params[9]);
	if(extraeven)
		even[nodd]+=(char)floor(odd[nodd-1]*(params[4]+params[9]));

	
	for(int k=0;k<nodd-!extraeven;++k)//un-update
		odd[k]-=(char)floor(even[k]*params[3]+even[k+1]*params[8]);
	if(!extraeven)
		odd[nodd-1]-=(char)floor(even[nodd-1]*(params[3]+params[8]));
	
	even[0]+=(char)floor(odd[0]*(params[2]+params[7]));
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(char)floor(odd[k-1]*params[2]+odd[k]*params[7]);
	if(extraeven)
		even[nodd]+=(char)floor(odd[nodd-1]*(params[2]+params[7]));

	
	for(int k=0;k<nodd-!extraeven;++k)//un-update
		odd[k]-=(char)floor(even[k]*params[1]+even[k+1]*params[6]);
	if(!extraeven)
		odd[nodd-1]-=(char)floor(even[nodd-1]*(params[1]+params[6]));
	
	even[0]+=(char)floor(odd[0]*(params[0]+params[5]));
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(char)floor(odd[k-1]*params[0]+odd[k]*params[5]);
	if(extraeven)
		even[nodd]+=(char)floor(odd[nodd-1]*(params[0]+params[5]));


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_custom_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp, const double *params)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_custom_fwd(buffer+rowlen*ky, w2, stride, temp, params);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_custom_fwd(buffer+stride*kx, h2, rowlen, temp, params);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_custom_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp, const double *params)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_custom_inv(buffer+stride*kx, h2, rowlen, temp, params);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_custom_inv(buffer+rowlen*ky, w2, stride, temp, params);
	}
}

void dwt1d_exp_fwd(char *buffer, int count, int stride, char *b2, const double *paramsx12)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(int k=0;k<nodd+extraeven;++k)//predict
	{
		char
			prev5=k-4>=0?odd[k-4]:0,
			prev4=k-3>=0?odd[k-3]:0,
			prev3=k-2>=0?odd[k-2]:0,
			prev2=k-1>=0?odd[k-1]:0,
			prev=k<nodd?odd[k]:0,
			next=k+1<nodd?odd[k+1]:0,
			next2=k+2<nodd?odd[k+2]:0,
			next3=k+3<nodd?odd[k+3]:0,
			next4=k+4<nodd?odd[k+4]:0,
			next5=k+5<nodd?odd[k+5]:0;
		char pred=(char)(0.5*floor(
			paramsx12[0]*(prev+next)+
			paramsx12[1]*(prev2+next2)+
			paramsx12[2]*(prev3+next3)+
			paramsx12[3]*(prev4+next4)+
			paramsx12[4]*(prev5+next5)
		));
		even[k]-=pred;
		//even[k]-=(9*(prev+next)+prev2+next2)>>4;
	}
	for(int k=0;k<nodd;++k)//update
	{
		char
			prev5=k-5>=0?even[k-5]:0,
			prev4=k-4>=0?even[k-4]:0,
			prev3=k-3>=0?even[k-3]:0,
			prev2=k-2>=0?even[k-2]:0,
			prev =k-1>=0?even[k-1]:0,
			next =even[k],
			next2=k+1<nodd+extraeven?even[k+1]:0,
			next3=k+2<nodd+extraeven?even[k+2]:0,
			next4=k+3<nodd+extraeven?even[k+3]:0,
			next5=k+4<nodd+extraeven?even[k+4]:0;
		char update=(char)(0.5*floor(
			paramsx12[5]*(prev+next)+
			paramsx12[6]*(prev2+next2)+
			paramsx12[7]*(prev3+next3)+
			paramsx12[8]*(prev4+next4)+
			paramsx12[9]*(prev5+next5)
		));
		odd[k]+=update;
	}


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_exp_inv(char *buffer, int count, int stride, char *b2, const double *paramsx12)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	
	for(int k=0;k<nodd;++k)//unupdate
	{
		char
			prev5=k-5>=0?even[k-5]:0,
			prev4=k-4>=0?even[k-4]:0,
			prev3=k-3>=0?even[k-3]:0,
			prev2=k-2>=0?even[k-2]:0,
			prev =k-1>=0?even[k-1]:0,
			next =even[k],
			next2=k+1<nodd+extraeven?even[k+1]:0,
			next3=k+2<nodd+extraeven?even[k+2]:0,
			next4=k+3<nodd+extraeven?even[k+3]:0,
			next5=k+4<nodd+extraeven?even[k+4]:0;
		char update=(char)(0.5*floor(
			paramsx12[5]*(prev+next)+
			paramsx12[6]*(prev2+next2)+
			paramsx12[7]*(prev3+next3)+
			paramsx12[8]*(prev4+next4)+
			paramsx12[9]*(prev5+next5)
		));
		odd[k]-=update;
	}
	for(int k=0;k<nodd+extraeven;++k)//unpredict
	{
		char
			prev5=k-4>=0?odd[k-4]:0,
			prev4=k-3>=0?odd[k-3]:0,
			prev3=k-2>=0?odd[k-2]:0,
			prev2=k-1>=0?odd[k-1]:0,
			prev=k<nodd?odd[k]:0,
			next=k+1<nodd?odd[k+1]:0,
			next2=k+2<nodd?odd[k+2]:0,
			next3=k+3<nodd?odd[k+3]:0,
			next4=k+4<nodd?odd[k+4]:0,
			next5=k+5<nodd?odd[k+5]:0;
		char pred=(char)(0.5*floor(
			paramsx12[0]*(prev+next)+
			paramsx12[1]*(prev2+next2)+
			paramsx12[2]*(prev3+next3)+
			paramsx12[3]*(prev4+next4)+
			paramsx12[4]*(prev5+next5)
		));
		even[k]+=pred;
		//even[k]+=(9*(prev+next)+prev2+next2)>>4;
	}


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_exp_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp, const double *params)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_exp_fwd(buffer+rowlen*ky, w2, stride, temp, params);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_exp_fwd(buffer+stride*kx, h2, rowlen, temp, params);
	}
}
void dwt2d_exp_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp, const double *params)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_exp_inv(buffer+stride*kx, h2, rowlen, temp, params);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_exp_inv(buffer+rowlen*ky, w2, stride, temp, params);
	}
}

//lifting-based 8bit lazy DWT
void dwt1d_lazy_fwd(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
void dwt1d_lazy_inv(char *buffer, int count, int stride, char *b2)
{
	int nodd=count>>1, extraeven=count&1;
	char *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_lazy_fwd(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_lazy_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_lazy_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_lazy_inv(char *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, char *temp)
{
	if(sizes_start>=sizes_end-1)
		return;
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_lazy_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_lazy_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}
#endif

//lazy DWT (1 step)
static void dwt1d_lazy_fwd(int *buffer, int count, int stride, int *b2)
{
	int nodd=count>>1, extraeven=count&1;
	int *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
static void dwt1d_lazy_inv(int *buffer, int count, int stride, int *b2)
{
	int nodd=count>>1, extraeven=count&1;
	int *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_lazy_fwd(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	if(sizes_start>=sizes_end-1)
		return;
	//for(int it=sizes_start;it<sizes_end-1;++it)
	{
		//int w2=sizes[it].w, h2=sizes[it].h;
		int w2=sizes[sizes_start].w, h2=sizes[sizes_start].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_lazy_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_lazy_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
		//break;
	}
}
void dwt2d_lazy_inv(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	if(sizes_start>=sizes_end-1)
		return;
	//for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		//int w2=sizes[it].w, h2=sizes[it].h;
		int w2=sizes[sizes_start].w, h2=sizes[sizes_start].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_lazy_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_lazy_inv(buffer+rowlen*ky, w2, stride, temp);
		//break;
	}
}

//lifting-based 8bit Haar DWT
static void dwt1d_haar_fwd(int *buffer, int count, int stride, int *b2)
{
	int nodd=count>>1, extraeven=count&1;
	int *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(int k=0;k<nodd;++k)//const predictor (difference)
		even[k]-=odd[k];
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(int k=0;k<nodd;++k)//update (O+(E-O)/2 = average)
		odd[k]+=even[k]>>1;


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
static void dwt1d_haar_inv(int *buffer, int count, int stride, int *b2)
{
	int nodd=count>>1, extraeven=count&1;
	int *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	
	for(int k=0;k<nodd;++k)//unupdate
		odd[k]-=even[k]>>1;
	
	for(int k=0;k<nodd;++k)//unpredict
		even[k]+=odd[k];
	if(extraeven)
		even[nodd]+=odd[nodd-1];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_haar_fwd(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	if(sizes_start>=sizes_end-1)
		return;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_haar_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_haar_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_haar_inv(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	if(sizes_start>=sizes_end-1)
		return;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_haar_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_haar_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit squeeze DWT  from JPEG XL
static int smoothtendency(int B, int a, int n)
{
	int diff=0;
	if(B>=a&&a>=n)
	{
		diff=(4*B-3*n-a+6)/12;
		//      2C = a<<1 + diff - diff&1 <= 2B  so diff - diff&1 <= 2B - 2a
		//      2D = a<<1 - diff - diff&1 >= 2n  so diff + diff&1 <= 2a - 2n
		if(diff-(diff&1)>2*(B-a))
			diff=2*(B-a)+1;
		if(diff+(diff&1)>2*(a-n))
			diff=2*(a-n);
	}
	else if(B<=a&&a<=n)
	{
		diff=(4*B-3*n-a-6)/12;
		//      2C = a<<1 + diff + diff&1 >= 2B  so diff + diff&1 >= 2B - 2a
		//      2D = a<<1 - diff + diff&1 <= 2n  so diff - diff&1 >= 2a - 2n
		if(diff+(diff&1)<2*(B-a))
			diff=2*(B-a)-1;
		if(diff-(diff&1)<2*(a-n))
			diff=2*(a-n);
	}
	return diff;
}
static void dwt1d_squeeze_fwd(int *buffer, int count, int stride, int *b2)//nOdd <= nEven			nOdd = nEven - (count&1)
{
	int nodd=count>>1, extraeven=count&1;
	int *odd=b2, *even=b2+nodd;

	for(int ks=0, kd=0;kd<nodd+extraeven;ks+=stride<<1, ++kd)
	{
		int
			o1=kd>0?buffer[ks-stride]:0,
			e =buffer[ks],
			o =(kd<<1)+1<count?buffer[ks+stride]:0,
			e2=(kd<<1)+2<count?buffer[ks+(stride<<1)]:0,//n-1-(idx-(n-1)) = ((n-1)<<1)-idx
			o2=(kd<<1)+3<count?buffer[ks+ stride*3  ]:0;

		e-=o;		//diff
		o+=e>>1;	//av
		e2-=o2;
		o2+=e2>>1;
		e-=smoothtendency(o1, o, o2);
		
		if(kd<nodd)
			odd[kd]=o;
		even[kd]=e;
	}


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
static void dwt1d_squeeze_inv(int *buffer, int count, int stride, int *b2)
{
	int nodd=count>>1, extraeven=count&1;
	int *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	for(int ks=0, kd=0;ks<nodd+extraeven;++ks, kd+=stride<<1)
	{
		int
			o1=ks>0?buffer[kd-stride]:0,
			e =even[ks],
			o=ks<nodd?odd[ks]:0,
			o2=ks+1<nodd?odd[ks+1]:0;

		e+=smoothtendency(o1, o, o2);
		o-=e>>1;
		e+=o;

		buffer[kd]=e;
		if(ks<nodd)
			buffer[kd+stride]=o;
	}
}
void dwt2d_squeeze_fwd(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	if(sizes_start>=sizes_end-1)
		return;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_squeeze_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, 128, "squeeze-stage%02dA.PNG", it);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_squeeze_fwd(buffer+stride*kx, h2, rowlen, temp);

		//save_channel(buffer, iw, ih, 4, 128, "squeeze-stage%02dB.PNG", it);
	}
}
void dwt2d_squeeze_inv(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	if(sizes_start>=sizes_end-1)
		return;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_squeeze_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_squeeze_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit CDF 5/3 DWT
static void dwt1d_cdf53_fwd(int *buffer, int count, int stride, int *b2)
{
	int nodd=count>>1, extraeven=count&1;
	int *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

#if 0
	if(fabs(customparam_st[0])<1.5)//linear
	{
		even[0]-=odd[0];
		for(int k=1;k<nodd;++k)//linear predictor (deviation from linearity)
			even[k]-=(odd[k-1]+odd[k]+1)>>1;
		if(extraeven)
			even[nodd]-=odd[nodd-1];
	
		for(int k=0;k<nodd-!extraeven;++k)//update (smoothing?)
			odd[k]+=(even[k]+even[k+1])>>2;
		if(!extraeven)
			odd[nodd-1]+=even[nodd-1]>>1;
	}
	else if(fabs(customparam_st[0])<3.5)//cubic
	{
		for(int k=0;k<nodd+extraeven;++k)
		{
			int prev3=k-2>=0?odd[k-2]:0, prev1=k-1>=0?odd[k-1]:0, next1=k<nodd?odd[k]:0, next3=k+1<nodd?odd[k+1]:0;
			even[k]-=(-(prev3+next3)+(prev1+next1)*9+8)>>4;
		}
	
		for(int k=0;k<nodd-!extraeven;++k)//update (smoothing)
			odd[k]+=(even[k]+even[k+1])>>2;
		if(!extraeven)
			odd[nodd-1]+=even[nodd-1]>>1;

		//for(int k=0;k<nodd;++k)
		//{
		//	int prev3=k-1>=0?even[k-1]:0, prev1=even[k], next1=k+1<nodd+extraeven?even[k+1]:0, next3=k+2<nodd+extraeven?even[k+2]:0;
		//	odd[k]+=((prev3+next3)+(prev1+next1)*3+4)>>3;
		//}
	}
	else if(fabs(customparam_st[0])<5.5)//power 5
	{
		for(int k=0;k<nodd+extraeven;++k)
		{
			int prev5=k-3>=0?odd[k-3]:0,
				prev3=k-2>=0?odd[k-2]:0,
				prev1=k-1>=0?odd[k-1]:0,
				next1=k<nodd?odd[k]:0,
				next3=k+1<nodd?odd[k+1]:0,
				next5=k+2<nodd?odd[k+2]:0;
			even[k]-=(3*(prev5+next5)-25*(prev3+next3)+150*(prev1+next1)+128)>>8;
		}
	
		for(int k=0;k<nodd-!extraeven;++k)//update (smoothing?)
			odd[k]+=(even[k]+even[k+1])>>2;
		if(!extraeven)
			odd[nodd-1]+=even[nodd-1]>>1;
	}
	else if(fabs(customparam_st[0])<7.5)//power 7
	{
		for(int k=0;k<nodd+extraeven;++k)
		{
			int prev7=k-4>=0?odd[k-4]:0,
				prev5=k-3>=0?odd[k-3]:0,
				prev3=k-2>=0?odd[k-2]:0,
				prev1=k-1>=0?odd[k-1]:0,
				next1=k<nodd?odd[k]:0,
				next3=k+1<nodd?odd[k+1]:0,
				next5=k+2<nodd?odd[k+2]:0,
				next7=k+3<nodd?odd[k+3]:0;
			even[k]-=(-5*(prev7+next7)+49*(prev5+next5)-245*(prev3+next3)+1225*(prev1+next1)+1024)>>11;
		}
	
		for(int k=0;k<nodd-!extraeven;++k)//update (smoothing?)
			odd[k]+=(even[k]+even[k+1])>>2;
		if(!extraeven)
			odd[nodd-1]+=even[nodd-1]>>1;
	}
	else if(fabs(customparam_st[0])<9.5)//power 9
	{
		for(int k=0;k<nodd+extraeven;++k)
		{
			int prev9=k-5>=0?odd[k-5]:0,
				prev7=k-4>=0?odd[k-4]:0,
				prev5=k-3>=0?odd[k-3]:0,
				prev3=k-2>=0?odd[k-2]:0,
				prev1=k-1>=0?odd[k-1]:0,
				next1=k<nodd?odd[k]:0,
				next3=k+1<nodd?odd[k+1]:0,
				next5=k+2<nodd?odd[k+2]:0,
				next7=k+3<nodd?odd[k+3]:0,
				next9=k+4<nodd?odd[k+4]:0;
			even[k]-=(35*(prev9+next9)-405*(prev7+next7)+2268*(prev5+next5)-8820*(prev3+next3)+39690*(prev1+next1)+0x8000)>>16;
		}
	
		for(int k=0;k<nodd-!extraeven;++k)//update (smoothing?)
			odd[k]+=(even[k]+even[k+1])>>2;
		if(!extraeven)
			odd[nodd-1]+=even[nodd-1]>>1;
	}
#endif

	//orginal CDF 5/3 DWT
#if 1
	even[0]-=odd[0];
	for(int k=1;k<nodd;++k)//linear predictor (deviation from linearity)
		even[k]-=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(int k=0;k<nodd-!extraeven;++k)//update (smoothing)
		odd[k]+=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]>>1;
#endif


	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
static void dwt1d_cdf53_inv(int *buffer, int count, int stride, int *b2)
{
	int nodd=count>>1, extraeven=count&1;
	int *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];

	
	for(int k=0;k<nodd-!extraeven;++k)//un-update
		odd[k]-=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]>>1;
	
	even[0]+=odd[0];
	for(int k=1;k<nodd;++k)//un-predict
		even[k]+=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]+=odd[nodd-1];


	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf53_fwd(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	if(sizes_start>=sizes_end-1)
		return;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf53_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);
		//snprintf(g_buf, G_BUF_SIZE, "cdf53-stage%02dA.PNG", it);
		//lodepng_encode_file(g_buf, buffer, iw, ih, LCT_RGBA, 8);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf53_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_cdf53_inv(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	if(sizes_start>=sizes_end-1)
		return;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_cdf53_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_cdf53_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

//lifting-based 8bit CDF 9/7 DWT
static const int cdf97_coeff[]=
{
	//'factring wavelet transforms into lifting steps' - page 19
	//'a wavelet tour of signal processing - the sparse way' - page 376
	//Lifting-based Discrete Wavelet Transform for Real-Time Signal Detection 2015-10
	-0x1960C,	//-1.58613434342059f,	//alpha
	-0x00D90,	//-0.0529801185729f,	//beta
	 0x0E206,	// 0.8829110755309f,	//gamma
	 0x07189,	// 0.4435068520439f,	//delta
	 0x1264C,	// 1.1496043988602f,	//zeta		output gain is 1.89
};
static void dwt1d_u8_predict(int *odd, int *even, int nodd, int extraeven, int coeff)
{
	even[0]+=odd[0]*coeff>>15;
	for(int k=1;k<nodd;++k)//predict
		even[k]+=(odd[k-1]+odd[k])*coeff>>16;
	if(extraeven)
		even[nodd]+=odd[nodd-1]*coeff>>15;
}
static void dwt1d_u8_unpredict(int *odd, int *even, int nodd, int extraeven, int coeff)
{
	even[0]-=odd[0]*coeff>>15;
	for(int k=1;k<nodd;++k)//unpredict
		even[k]-=(odd[k-1]+odd[k])*coeff>>16;
	if(extraeven)
		even[nodd]-=odd[nodd-1]*coeff>>15;
}
static void dwt1d_u8_update(int *odd, int *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//update
		odd[k]+=(even[k]+even[k+1])*coeff>>16;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]*coeff>>15;
}
static void dwt1d_u8_unupdate(int *odd, int *even, int nodd, int extraeven, int coeff)
{
	for(int k=0;k<nodd-!extraeven;++k)//unupdate
		odd[k]-=(even[k]+even[k+1])*coeff>>16;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]*coeff>>15;
}
//static void dwt1d_u8_scale(int *buf, int count, int coeff)
//{
//	for(int k=0;k<count;++k)
//		buf[k]=buf[k]*coeff>>16;
//}
//static void dwt1d_u8_unscale(int *buf, int count, int coeff)
//{
//	for(int k=0;k<count;++k)
//		buf[k]=(buf[k]<<16)/coeff;
//}
static void dwt1d_cdf97_fwd(int *buffer, int count, int stride, int *b2)
{
	int nodd=count>>1, extraeven=count&1;
	int *odd=b2, *even=b2+nodd;
	
	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//lazy wavelet: split into odd (low frequency) & even (high frequency)
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

	dwt1d_u8_predict(odd, even, nodd, extraeven, cdf97_coeff[0]);
	dwt1d_u8_update (odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_u8_predict(odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_u8_update (odd, even, nodd, extraeven, cdf97_coeff[3]);
	//dwt1d_u8_scale(b2, count, cdf97_coeff[4]);

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
static void dwt1d_cdf97_inv(int *buffer, int count, int stride, int *b2)
{
	int nodd=count>>1, extraeven=count&1;
	int *odd=b2, *even=b2+nodd;

	for(int k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	//dwt1d_u8_unscale(b2, count, cdf97_coeff[4]);
	dwt1d_u8_unupdate (odd, even, nodd, extraeven, cdf97_coeff[3]);
	dwt1d_u8_unpredict(odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_u8_unupdate (odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_u8_unpredict(odd, even, nodd, extraeven, cdf97_coeff[0]);

	for(int k=0, ks=0;k<nodd;++k, ks+=stride<<1)//inv lazy wavelet: join even & odd
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf97_fwd(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	if(sizes_start>=sizes_end-1)
		return;
	for(int it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf97_fwd(buffer+rowlen*ky, w2, stride, temp);

		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dA.PNG", it);

		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf97_fwd(buffer+stride*kx, h2, rowlen, temp);
		
		//save_channel(buffer, iw, ih, 4, "cdf53-stage%02dB.PNG", it);
	}
}
void dwt2d_cdf97_inv(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw=sizes->w, rowlen=stride*iw;
//	int ih=sizes->h;
	if(sizes_start>=sizes_end-1)
		return;
	for(int it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(int kx=0;kx<w2;++kx)//vertical IDWT
			dwt1d_cdf97_inv(buffer+stride*kx, h2, rowlen, temp);

		for(int ky=0;ky<h2;++ky)//horizontal IDWT
			dwt1d_cdf97_inv(buffer+rowlen*ky, w2, stride, temp);
	}
}

#if 0
void dwt2d_dec_fwd(char *buffer, int iw, int ih)
{
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	int w2=iw, h2=ih;
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf97_fwd(buffer+4*iw*ky+kc, w2, 4, temp);
		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf97_fwd(buffer+4*kx+kc, h2, 4*iw, temp);
	}
	w2>>=1, h2>>=1;
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_cdf53_fwd(buffer+4*iw*ky+kc, w2, 4, temp);
		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_cdf53_fwd(buffer+4*kx+kc, h2, 4*iw, temp);
	}
	w2>>=1, h2>>=1;
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<h2;++ky)//horizontal DWT
			dwt1d_haar_fwd(buffer+4*iw*ky+kc, w2, 4, temp);
		for(int kx=0;kx<w2;++kx)//vertical DWT
			dwt1d_haar_fwd(buffer+4*kx+kc, h2, 4*iw, temp);
	}
	free(temp);
}
void dwt2d_dec_inv(char *buffer, int iw, int ih)
{
	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 3, 3, 3);
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	DWTSize *psizes=(DWTSize*)sizes->data;
	int nsizes=(int)sizes->count;
	for(int kc=0;kc<3;++kc)
	{
		for(int kx=0;kx<psizes[2].w;++kx)//vertical invDWT
			dwt1d_haar_inv(buffer+4*kx+kc, psizes[2].h, 4*iw, temp);
		for(int ky=0;ky<psizes[2].h;++ky)//horizontal invDWT
			dwt1d_haar_inv(buffer+4*iw*ky+kc, psizes[2].w, 4, temp);
	}
	for(int kc=0;kc<3;++kc)
	{
		for(int kx=0;kx<psizes[1].w;++kx)//vertical invDWT
			dwt1d_cdf53_inv(buffer+4*kx+kc, psizes[1].h, 4*iw, temp);
		for(int ky=0;ky<psizes[1].h;++ky)//horizontal invDWT
			dwt1d_cdf53_inv(buffer+4*iw*ky+kc, psizes[1].w, 4, temp);
	}
	for(int kc=0;kc<3;++kc)
	{
		for(int kx=0;kx<psizes[0].w;++kx)//vertical invDWT
			dwt1d_cdf97_inv(buffer+4*kx+kc, psizes[0].h, 4*iw, temp);
		for(int ky=0;ky<psizes[0].h;++ky)//horizontal invDWT
			dwt1d_cdf97_inv(buffer+4*iw*ky+kc, psizes[0].w, 4, temp);
	}
	free(temp);
	array_free(&sizes);
}

#define LOAD(BUF, X, Y) ((unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?BUF[(iw*(ky+(Y))+kx+(X))<<2|kc]:0)
static void pred_wu97_pass3_row(const char *src, char *dst, int iw, int ih, int fwd, int kc, int kx0, int ky)
{
	const char *pixels=fwd?src:dst, *errors=fwd?dst:src;
	for(int kx=kx0;kx<iw;kx+=2)
	{
		char
			NW=LOAD(pixels, -1, -1),
			N =LOAD(pixels,  0, -1),
			NE=LOAD(pixels,  1, -1),
			W =LOAD(pixels, -1,  0),
			E =LOAD(pixels,  1,  0),
			S =LOAD(pixels,  0,  1);
		int pred=((N+W+S+E)*3-(NW+NE)*2+4)>>3;
		pred=CLAMP(-128, pred, 127);

		int idx=(iw*ky+kx)<<2|kc;
		if(fwd)
			dst[idx]=src[idx]-pred;
		else
			dst[idx]=src[idx]+pred;
	}
}
static void pred_wu97_pass3(const char *src, char *dst, int iw, int ih, int fwd)
{
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<ih;ky+=2)
		{
			pred_wu97_pass3_row(src, dst, iw, ih, fwd, kc, 1, ky);
			pred_wu97_pass3_row(src, dst, iw, ih, fwd, kc, 0, ky+1);
		}
	}
}
static void pred_wu97_pass2_fwd(const char *pixels, char *errors, int iw, int ih)//diamonds to sticks
{
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int ky=0;ky<ih;ky+=2)//pass2
		{
			for(int kx=0;kx<iw;kx+=2)
			{
				char
					NN  =LOAD(pixels,  0, -2),
					NW  =LOAD(pixels, -1, -1),
					NE  =LOAD(pixels,  1, -1),
					WW  =LOAD(pixels, -2,  0),
					curr=LOAD(pixels,  0,  0),
					EE  =LOAD(pixels,  2,  0),
					SW  =LOAD(pixels, -1,  1),
					SE  =LOAD(pixels,  1,  1),
					SEEE=LOAD(pixels,  3,  1),
					SS  =LOAD(pixels,  0,  2),
					SSSE=LOAD(pixels,  1,  3);
				
				char
					dcurr=curr-SE,
					acurr=SE+(dcurr>>1),
					aE=SEEE+((char)(EE-SEEE)>>1),
					aS=SSSE+((char)(SS-SSSE)>>1);

				//if(kx==(iw>>1)&&ky==(ih>>1))//
				//if(kc==0&&kx==0&&ky==0)//
				//	printf("");

				int pred=(acurr*0xE666+(NE+NW+SW)*0x2AAB-(NN+WW)*0x0CCD-(aE+aS)*0x2666+0x8000)>>16;
				pred=CLAMP(-128, pred, 127);
				pred-=acurr;
				pred<<=1;
				//pred=CLAMP(-128, pred, 127);

				int idx=(iw*ky+kx)<<2|kc;
				errors[idx]=acurr;
				errors[idx+((iw+1)<<2)]=dcurr-pred;
				//errors[idx+((iw+1)<<2)]=dcurr;
			}
		}
#endif
#if 1
		for(int ky=(ih-2)&(-2);ky>=0;ky-=2)//pass1
		{
			for(int kx=(iw-2)&(-2);kx>=0;kx-=2)
			{
				char
					NNWW=LOAD(errors, -2, -2),
					NN  =LOAD(errors,  0, -2),
					NNEE=LOAD(errors,  2, -2),
					WW  =LOAD(errors, -2,  0);
				//int pred=(NN+WW)/2+(NNEE-NNWW)/4;

				int vmin, vmax;
				if(NN<WW)
					vmin=NN, vmax=WW;
				else
					vmin=WW, vmax=NN;
				int pred=NN+WW-NNWW;
				pred=CLAMP(vmin, pred, vmax);

				//if(kc==0&&kx==2&&ky==0)//
				//	printf("");
				//if(kc==0&&kx==0&&ky==2)//
				//	printf("");

				int idx=(iw*ky+kx)<<2|kc;
				errors[idx]-=pred;
			}
		}
#endif
	}
}
static void pred_wu97_pass2_inv(const char *errors, char *pixels, int iw, int ih)//sticks to diamonds
{
	for(int kc=0;kc<3;++kc)
	{
#if 1
		for(int ky=0;ky<ih;ky+=2)//pass1
		{
			for(int kx=0;kx<iw;kx+=2)
			{
				char
					NNWW=LOAD(pixels, -2, -2),
					NN  =LOAD(pixels,  0, -2),
					NNEE=LOAD(pixels,  2, -2),
					WW  =LOAD(pixels, -2,  0);
				//int pred=(NN+WW)/2+(NNEE-NNWW)/4;
				
				int vmin, vmax;
				if(NN<WW)
					vmin=NN, vmax=WW;
				else
					vmin=WW, vmax=NN;
				int pred=NN+WW-NNWW;
				pred=CLAMP(vmin, pred, vmax);

				//if(kc==0&&kx==2&&ky==0)//
				//	printf("");
				//if(kc==0&&kx==0&&ky==2)//
				//	printf("");

				int idx=(iw*ky+kx)<<2|kc;
				pixels[idx]=errors[idx]+pred;

				//idx+=(iw+1)<<2;
				//pixels[idx]=errors[idx];//copy diff
			}
		}
#endif
#if 1
		for(int ky=0;ky<ih;ky+=2)//pass2
		{
			for(int kx=0;kx<iw;kx+=2)
			{
				char
					NN  =LOAD(pixels,  0, -2),
					NW  =LOAD(pixels, -1, -1),
					NE  =LOAD(pixels,  1, -1),
					WW  =LOAD(pixels, -2,  0),
					SW  =LOAD(pixels, -1,  1);
#if 1
				char//original
					acurr=LOAD(pixels, 0, 0),
					dcurr=LOAD(errors, 1, 1),
					aE   =LOAD(pixels, 2, 0),
					aS   =LOAD(pixels, 0, 2);
#endif
#if 0
				char//pass2 only
					acurr=LOAD(errors, 0, 0),
					dcurr=LOAD(errors, 1, 1),
					aE   =LOAD(errors, 2, 0),
					aS   =LOAD(errors, 0, 2);
#endif

				//if(kc==0&&kx==0&&ky==0)//
				//	printf("");

				int pred=(acurr*0xE666+(NE+NW+SW)*0x2AAB-(NN+WW)*0x0CCD-(aE+aS)*0x2666+0x8000)>>16;
				pred=CLAMP(-128, pred, 127);
				pred-=acurr;
				pred<<=1;
				//pred=CLAMP(-128, pred, 127);

				dcurr+=pred;

				char curr=dcurr, SE=acurr;
				SE-=curr>>1;	//SE = av - floor(diff/2)
				curr+=SE;		//curr = diff + SE

				int idx=(iw*ky+kx)<<2|kc;
				pixels[idx]=curr;
				pixels[idx+((iw+1)<<2)]=SE;
			}
		}
#endif
	}
}
#undef  LOAD
void pred_wu97(char *buf, int iw, int ih, int fwd)//'Lossless Compression of Continuous-Tone Images via Context Selection, Quantization, and Modeling'
{
	int res=iw*ih;
	char *b2=(char*)malloc((size_t)res<<2);
	char *temp=(char*)malloc(MAXVAR(iw, ih));
	if(!b2||!temp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	//int black=0xFF000000;
	//memfill(b2, &black, (size_t)res<<2, sizeof(int));
	memcpy(b2, buf, (size_t)res<<2);//copy alpha
	ArrayHandle sizes=dwt2d_gensizes(iw, ih, 0, 0, 1);
	if(fwd)
	{
		pred_wu97_pass3(buf, b2, iw, ih, fwd);
		pred_wu97_pass2_fwd(buf, b2, iw, ih);
		dwt2d_lazy_fwd(b2  , (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
		dwt2d_lazy_fwd(b2+1, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
		dwt2d_lazy_fwd(b2+2, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);

	}
	else
	{
		dwt2d_lazy_inv(buf  , (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
		dwt2d_lazy_inv(buf+1, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
		dwt2d_lazy_inv(buf+2, (DWTSize*)sizes->data, 0, (int)sizes->count, 4, temp);
		pred_wu97_pass2_inv(buf, b2, iw, ih);
		pred_wu97_pass3(buf, b2, iw, ih, fwd);
	}
	array_free(&sizes);
	memcpy(buf, b2, (size_t)res<<2);
	free(temp);
	free(b2);
}
#endif


//DCTs
static void dct4_fwd_i8(int *x)
{
	x[3]=x[0]-x[3];
	x[0]-=x[3]>>1;
	x[2]=x[1]-x[2];
	x[1]-=x[2]>>1;
	x[1]=x[0]-x[1];
	x[0]-=x[1]>>1;
	x[2]=(x[3]*13>>5)-x[2];
	x[3]-=x[2]*11>>5;
}
static void dct4_inv_i8(int *x)
{
	x[3]+=x[2]*11>>5;
	x[2]=(x[3]*13>>5)-x[2];
	x[0]+=x[1]>>1;
	x[1]=x[0]-x[1];
	x[1]+=x[2]>>1;
	x[2]=x[1]-x[2];
	x[0]+=x[3]>>1;
	x[3]=x[0]-x[3];
}
void image_dct4_fwd(Image *image)
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
					image->data[ idx   <<2|kc],
					image->data[(idx+1)<<2|kc],
					image->data[(idx+2)<<2|kc],
					image->data[(idx+3)<<2|kc],
				};

				dct4_fwd_i8(x);

				temp[ kx>>2                  ]=x[0];
				temp[(kx>>2)+(image->iw>>2)  ]=x[1];
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
					image->data[ idx             <<2|kc],
					image->data[(idx+image->iw  )<<2|kc],
					image->data[(idx+image->iw*2)<<2|kc],
					image->data[(idx+image->iw*3)<<2|kc],
				};

				dct4_fwd_i8(x);

				temp[(ky>>2)                 ]=x[0];
				temp[(ky>>2)+(image->ih>>2)  ]=x[1];
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
		image->depth[kc]+=(image->depth[kc]+2<=24)<<1;
}
void image_dct4_inv(Image *image)
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
					temp[(ky>>2)                 ],
					temp[(ky>>2)+(image->ih>>2)  ],
					temp[(ky>>2)+(image->ih>>2)*2],
					temp[(ky>>2)+(image->ih>>2)*3],
				};

				dct4_inv_i8(x);
				
				image->data[ idx             <<2|kc]=x[0];
				image->data[(idx+image->iw  )<<2|kc]=x[1];
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
					temp[(kx>>2)                 ],
					temp[(kx>>2)+(image->iw>>2)  ],
					temp[(kx>>2)+(image->iw>>2)*2],
					temp[(kx>>2)+(image->iw>>2)*3],
				};

				dct4_inv_i8(x);
				
				image->data[ idx   <<2|kc]=x[0];
				image->data[(idx+1)<<2|kc]=x[1];
				image->data[(idx+2)<<2|kc]=x[2];
				image->data[(idx+3)<<2|kc]=x[3];
			}
		}
#endif
	}
	free(temp);
	for(int kc=0;kc<3;++kc)
		image->depth[kc]-=(image->depth[kc]-2>=image->src_depth[kc])<<1;
}

static void dct8_fwd_i8(int *x)
{
	//binDCT-C7
	x[7]=x[0]-x[7];
	x[6]=x[1]-x[6];
	x[5]=x[2]-x[5];
	x[4]=x[3]-x[4];
	x[0]-=x[7]>>1;
	x[1]-=x[6]>>1;
	x[2]-=x[5]>>1;
	x[3]-=x[4]>>1;

	x[3]=x[0]-x[3];
	x[2]=x[1]-x[2];
	x[0]-=x[3]>>1;
	x[1]-=x[2]>>1;

	x[1]=x[0]-x[1];
	x[0]-=x[1]>>1;

	x[2]=(x[3]*13>>5)-x[2];
	x[3]-=x[2]*11>>5;

	x[5]-=x[6]*13>>5;
	x[6]+=x[5]*11>>4;
	x[5]=(x[6]*15>>5)-x[5];

	x[5]=x[4]-x[5];
	x[6]=x[7]-x[6];
	x[4]-=x[5]>>1;
	x[7]-=x[6]>>1;

	x[4]=(x[7]*3>>4)-x[4];
	x[7]-=x[4]*3>>4;

	x[5]+=x[6]*11>>4;
	x[6]-=x[5]*15>>5;
}
static void dct8_inv_i8(int *x)
{
	//invBinDCT-C7
	x[6]+=x[5]*15>>5;
	x[5]-=x[6]*11>>4;

	x[7]+=x[4]*3>>4;
	x[4]=(x[7]*3>>4)-x[4];

	x[7]+=x[6]>>1;
	x[4]+=x[5]>>1;
	x[6]=x[7]-x[6];
	x[5]=x[4]-x[5];

	x[5]=(x[6]*15>>5)-x[5];
	x[6]-=x[5]*11>>4;
	x[5]+=x[6]*13>>5;

	x[3]+=x[2]*11>>5;
	x[2]=(x[3]*13>>5)-x[2];

	x[0]+=x[1]>>1;
	x[1]=x[0]-x[1];

	x[1]+=x[2]>>1;
	x[0]+=x[3]>>1;
	x[2]=x[1]-x[2];
	x[3]=x[0]-x[3];

	x[3]+=x[4]>>1;
	x[2]+=x[5]>>1;
	x[1]+=x[6]>>1;
	x[0]+=x[7]>>1;
	x[4]=x[3]-x[4];
	x[5]=x[2]-x[5];
	x[6]=x[1]-x[6];
	x[7]=x[0]-x[7];
}
void image_dct8_fwd(Image *image)
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
					image->data[ idx   <<2|kc],
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

				dct8_fwd_i8(x);

				temp[ kx>>3                  ]=x[0];
				temp[(kx>>3)+(image->iw>>3)  ]=x[1];
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
					image->data[ idx             <<2|kc],
					image->data[(idx+image->iw  )<<2|kc],
					image->data[(idx+image->iw*2)<<2|kc],
					image->data[(idx+image->iw*3)<<2|kc],
					image->data[(idx+image->iw*4)<<2|kc],
					image->data[(idx+image->iw*5)<<2|kc],
					image->data[(idx+image->iw*6)<<2|kc],
					image->data[(idx+image->iw*7)<<2|kc],
				};

				dct8_fwd_i8(x);

				temp[(ky>>3)                 ]=x[0];
				temp[(ky>>3)+(image->ih>>3)  ]=x[1];
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
		image->depth[kc]+=(image->depth[kc]+3<=24)*3;
}
void image_dct8_inv(Image *image)
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
					temp[(ky>>3)                 ],
					temp[(ky>>3)+(image->ih>>3)  ],
					temp[(ky>>3)+(image->ih>>3)*2],
					temp[(ky>>3)+(image->ih>>3)*3],
					temp[(ky>>3)+(image->ih>>3)*4],
					temp[(ky>>3)+(image->ih>>3)*5],
					temp[(ky>>3)+(image->ih>>3)*6],
					temp[(ky>>3)+(image->ih>>3)*7],
				};

				dct8_inv_i8(x);
				
				image->data[ idx             <<2|kc]=x[0];
				image->data[(idx+image->iw  )<<2|kc]=x[1];
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
					temp[(kx>>3)                 ],
					temp[(kx>>3)+(image->iw>>3)  ],
					temp[(kx>>3)+(image->iw>>3)*2],
					temp[(kx>>3)+(image->iw>>3)*3],
					temp[(kx>>3)+(image->iw>>3)*4],
					temp[(kx>>3)+(image->iw>>3)*5],
					temp[(kx>>3)+(image->iw>>3)*6],
					temp[(kx>>3)+(image->iw>>3)*7],
				};

				dct8_inv_i8(x);
				
				image->data[ idx   <<2|kc]=x[0];
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
		image->depth[kc]-=(image->depth[kc]-3>=image->src_depth[kc])*3;
}


#if 0
void predict_dct3_prep(float *dct3, float *dct4)
{
	dct3[0]=2,          dct3[1]= 2, dct3[2]=2;
	dct3[3]=sqrtf(3.f), dct3[4]= 0, dct3[5]=-sqrtf(3.f);
	dct3[6]=1,          dct3[7]=-2, dct3[8]=1;

	for(int k=0;k<9;++k)
		dct3[k]/=3.f;

	float a=(float)cos(M_PI/8), b=(float)cos(M_PI/4), c=(float)cos(M_PI*3/8);
	dct4[ 0]=0.5f, dct4[ 1]= a, dct4[ 2]= b, dct4[ 3]= c;
	dct4[ 4]=0.5f, dct4[ 5]= a, dct4[ 6]=-b, dct4[ 7]=-c;
	dct4[ 8]=0.5f, dct4[ 9]=-a, dct4[10]=-b, dct4[11]= c;
	dct4[12]=0.5f, dct4[13]=-a, dct4[14]= b, dct4[15]=-c;
}
int  predict_dct3(const int *buf, int iw, int kx, int ky, int idx, int bytestride, int rowlen, const float *dct3, const float *dct4)
{
	float
		left[]=
	{
		kx-3>=0?(float)buf[idx-bytestride*3]:0,
		kx-2>=0?(float)buf[idx-bytestride*2]:0,
		kx-1>=0?(float)buf[idx-bytestride  ]:0,
	},
		top []=
	{
		ky-3>=0?(float)buf[idx-rowlen*3]:0,
		ky-2>=0?(float)buf[idx-rowlen*2]:0,
		ky-1>=0?(float)buf[idx-rowlen  ]:0,
	},
		topleft[]=
	{
		kx-3>=0&&ky-3>=0?(float)buf[idx-(rowlen+bytestride)*3]:0,
		kx-2>=0&&ky-2>=0?(float)buf[idx-(rowlen+bytestride)*2]:0,
		kx-1>=0&&ky-1>=0?(float)buf[idx-(rowlen+bytestride)  ]:0,
	},
		topright[]=
	{
		kx+3<iw&&ky-3>=0?(float)buf[idx-(rowlen-bytestride)*3]:0,
		kx+2<iw&&ky-2>=0?(float)buf[idx-(rowlen-bytestride)*2]:0,
		kx+1<iw&&ky-1>=0?(float)buf[idx-(rowlen-bytestride)  ]:0,
	};
	float x[]=
	{
		(dct3[0]*left[0]+dct3[1]*left[1]+dct3[2]*left[2] + dct3[0]*top[0]+dct3[1]*top[1]+dct3[2]*top[2] + dct3[0]*topleft[0]+dct3[1]*topleft[1]+dct3[2]*topleft[2] + dct3[0]*topright[0]+dct3[1]*topright[1]+dct3[2]*topright[2])*0.25f,
		(dct3[3]*left[0]+dct3[4]*left[1]+dct3[5]*left[2] + dct3[3]*top[0]+dct3[4]*top[1]+dct3[5]*top[2] + dct3[3]*topleft[0]+dct3[4]*topleft[1]+dct3[5]*topleft[2] + dct3[3]*topright[0]+dct3[4]*topright[1]+dct3[5]*topright[2])*0.25f,
		0,
		(dct3[6]*left[0]+dct3[7]*left[1]+dct3[8]*left[2] + dct3[6]*top[0]+dct3[7]*top[1]+dct3[8]*top[2] + dct3[6]*topleft[0]+dct3[7]*topleft[1]+dct3[8]*topleft[2] + dct3[6]*topright[0]+dct3[7]*topright[1]+dct3[8]*topright[2])*0.25f,
	};
	x[2]=x[1];
	float pred=dct4[12]*x[0]+dct4[13]*x[1]+dct4[14]*x[2]+dct4[15]*x[3];
	pred=roundf(pred);
	
	pred=CLAMP(customparam_clamp[0], pred, customparam_clamp[1]);
	return (int)pred;
}
void pred_dct3_fwd(int *buf, int iw, int ih, int nch, int bytestride)
{
	float coeff[25];
	predict_dct3_prep(coeff, coeff+9);

	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=(iw*ih-1)*bytestride+kc;
		for(int ky=ih-1;ky>=0;--ky)
		{
			for(int kx=iw-1;kx>=0;--kx, idx-=bytestride)
			{
				int pred=predict_dct3(buf, iw, kx, ky, idx, bytestride, rowlen, coeff, coeff+9);

				buf[idx]-=pred;
			}
		}
	}
}
void pred_dct3_inv(int *buf, int iw, int ih, int nch, int bytestride)
{
	float coeff[25];
	predict_dct3_prep(coeff, coeff+9);

	int rowlen=iw*bytestride;
	for(int kc=0;kc<nch;++kc)
	{
		int idx=kc;
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx, idx+=bytestride)
			{
				int pred=predict_dct3(buf, iw, kx, ky, idx, bytestride, rowlen, coeff, coeff+9);

				buf[idx]+=pred;
			}
		}
	}
}


//other
void image_split_fwd(char *image, int iw, int ih)
{
	char *b2=(char*)malloc((size_t)iw*ih<<2);
	if(!b2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(b2, image, (size_t)iw*ih<<2);
	//int maxdim=MAXVAR(iw, ih);
	//char *temp=(char*)malloc(maxdim);
	//if(!temp)
	//{
	//	LOG_ERROR("Alloc error");
	//	return;
	//}
	//memset(temp, 0, maxdim);
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<ih-1;ky+=2)
		{
			for(int kx=0;kx<iw-1;kx+=2)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					b2[ idx      <<2|kc],
					b2[(idx+1   )<<2|kc],
					b2[(idx  +iw)<<2|kc],
					b2[(idx+1+iw)<<2|kc],
				};
#if 0
				for(int k=1;k<4;++k)//mini-CDF (insertion sort)
				{
					int L=0, R=k-1, mid, found=0;
					while(L<=R)
					{
						mid=(L+R)>>1;
						if(x[mid]<x[k])
							L=mid+1;
						else if(x[mid]>x[k])
							R=mid-1;
						else
						{
							found=1;
							break;
						}
					}
					if(!found)
						mid=L+(L<k&&x[L]<x[k]);
				}
#endif
#if 1
				char temp;
				if(x[0]>x[1])temp=x[0], x[0]=x[1], x[1]=temp;//mini-CDF (dedicated sort)
				if(x[2]>x[3])temp=x[2], x[2]=x[3], x[3]=temp;//https://stackoverflow.com/questions/6145364/sort-4-number-with-few-comparisons
				if(x[0]>x[2])temp=x[0], x[0]=x[2], x[2]=temp;
				if(x[1]>x[3])temp=x[1], x[1]=x[3], x[3]=temp;
				if(x[1]>x[2])temp=x[1], x[1]=x[2], x[2]=temp;
#endif
				int idx2=iw*(ky>>1)+(kx>>1);
				image[ idx2                    <<2|kc]=x[0];
				image[(idx2           +(iw>>1))<<2|kc]=x[1];
				image[(idx2+iw*(ih>>1)        )<<2|kc]=x[2];
				image[(idx2+iw*(ih>>1)+(iw>>1))<<2|kc]=x[3];
			}
		}
#if 0
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw-1;kx+=2)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx   <<2|kc],
					image[(idx+1)<<2|kc],
				};

				temp[ kx>>1           ]=x[0];
				temp[(kx>>1)+(iw>>1)  ]=x[1];
			}
			for(int kx=0;kx<iw;++kx)
				image[(iw*ky+kx)<<2|kc]=temp[kx];
		}
#endif
#if 0
		for(int kx=0;kx<iw;++kx)
		{
			for(int ky=0;ky<ih-3;ky+=4)
			{
				int idx=iw*ky+kx;
				char x[]=
				{
					image[ idx      <<2|kc],
					image[(idx+iw  )<<2|kc],
				};

				temp[(ky>>1)          ]=x[0];
				temp[(ky>>1)+(ih>>1)  ]=x[1];
			}
			for(int ky=0;ky<ih;++ky)
				image[(iw*ky+kx)<<2|kc]=temp[ky];
		}
#endif
	}
	free(b2);
	//free(temp);
}
void image_split_inv(char *image, int iw, int ih)
{
	char *b2=(char*)malloc((size_t)iw*ih<<2);
	if(!b2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memcpy(b2, image, (size_t)iw*ih<<2);
	for(int kc=0;kc<3;++kc)
	{
		for(int ky=0;ky<ih-1;ky+=2)
		{
			for(int kx=0;kx<iw-1;kx+=2)
			{
				int idx=iw*ky+kx, idx2=iw*(ky>>1)+(kx>>1);
				char x[]=
				{
					b2[ idx2                    <<2|kc],
					b2[(idx2           +(iw>>1))<<2|kc],
					b2[(idx2+iw*(ih>>1)        )<<2|kc],
					b2[(idx2+iw*(ih>>1)+(iw>>1))<<2|kc],
				};
				
#if 0
				char temp;
				if(x[0]>x[1])temp=x[0], x[0]=x[1], x[1]=temp;//mini-CDF (dedicated sort)
				if(x[2]>x[3])temp=x[2], x[2]=x[3], x[3]=temp;//https://stackoverflow.com/questions/6145364/sort-4-number-with-few-comparisons
				if(x[0]>x[2])temp=x[0], x[0]=x[2], x[2]=temp;
				if(x[1]>x[3])temp=x[1], x[1]=x[3], x[3]=temp;
				if(x[1]>x[2])temp=x[1], x[1]=x[2], x[2]=temp;
#endif
				
				image[ idx      <<2|kc]=x[0];
				image[(idx+1   )<<2|kc]=x[1];
				image[(idx  +iw)<<2|kc]=x[2];
				image[(idx+1+iw)<<2|kc]=x[3];
			}
		}
	}
	free(b2);
}

static unsigned qhist[256]={0};
void channel_entropy(unsigned char *buf, int resolution, int nch, int bytestride, float *cr, int *usage)
{
#if 0
	if(debug_buf)
	{
		lodepng_encode_file("buf_debug.PNG", debug_buf, 768, resolution/768, LCT_RGBA, 8);
		lodepng_encode_file("buf.PNG", buf, 768, resolution/768, LCT_RGBA, 8);
		for(int k=0;k<resolution;++k)
		{
			debug_buf[k<<2  ]-=buf[k<<2  ]-128;
			debug_buf[k<<2|1]-=buf[k<<2|1]-128;
			debug_buf[k<<2|2]-=buf[k<<2|2]-128;
		}
		lodepng_encode_file("buf_diff.PNG", debug_buf, 768, resolution/768, LCT_RGBA, 8);
		for(int k=0;k<resolution;++k)
		{
			debug_buf[k<<2  ]+=buf[k<<2  ]-128;
			debug_buf[k<<2|1]+=buf[k<<2|1]-128;
			debug_buf[k<<2|2]+=buf[k<<2|2]-128;
		}
		//for(int k=0;k<resolution;++k)
		//{
		//	if(buf[k]!=(unsigned char)debug_buf[k])
		//	{
		//		LOG_ERROR("Error");
		//	}
		//}
	}
#endif

	double entropy[4]={0};
	memset(usage, 0, 4*sizeof(int));
	for(int kc=0;kc<nch;++kc)
	{
		memset(qhist, 0, 256*sizeof(unsigned));
		for(int k=0, end=resolution*bytestride;k<end;k+=bytestride)
		{
			unsigned char val=buf[k+kc];
			++qhist[val];
		}
		for(int ks=0;ks<256;++ks)
		{
			unsigned freq=qhist[ks];
			if(freq)
			{
				double p=(double)freq/resolution;
				//p*=0xFF00;
				//++p;
				//p/=0x10000;
				entropy[kc]-=p*log2(p);
				++usage[kc];
			}
		}
		cr[kc]=(float)(8/entropy[kc]);
	}
	
	//calculate csize with joint histogram
	unsigned *h2=(unsigned*)malloc(0x1000000*sizeof(unsigned));
	if(!h2)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(h2, 0, 0x1000000*sizeof(unsigned));
	for(int k=0;k<resolution;++k)
	{
		unsigned color=((int*)buf)[k]&0xFFFFFF;
		++h2[color];
	}
	for(int k=0;k<0x1000000;++k)
		usage[3]+=h2[k]!=0;
	double csize=0;
	for(int k=0;k<resolution;++k)
	{
		unsigned color=((int*)buf)[k]&0xFFFFFF;
		double p=(double)h2[color]/resolution;
		//p*=0xFFFFFFFF;
		//++p;
		//p/=0x100000000;
		double bitsize=-log2(p);
		csize+=bitsize;
	}
	free(h2);
	csize/=8;
	cr[3]=(float)(resolution*3/csize);
}
void jointhistogram(unsigned char *buf, int iw, int ih, int nbits, ArrayHandle *hist, int space_not_color)
{
	int nlevels=1<<nbits, hsize=nlevels*nlevels*nlevels;
	if(*hist)
		array_free(hist);
	unsigned *htemp=(unsigned*)malloc(hsize*sizeof(unsigned));
	if(!htemp)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(htemp, 0, hsize*sizeof(unsigned));

	int res=iw*ih;
	switch(space_not_color)
	{
	case 0://show correlation in color
		for(int k=0;k<res;++k)
		{
			unsigned char r=buf[k<<2]>>(8-nbits), g=buf[k<<2|1]>>(8-nbits), b=buf[k<<2|2]>>(8-nbits);
			int color=b<<(nbits<<1)|g<<nbits|r;

			++htemp[color];
		}
		break;
	case 1://show correlation in space x (CURR, W, WW)
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
				unsigned char
					v2=kx-2>=0?buf[(iw*ky+kx-2)<<2|1]>>(8-nbits):0,//WW
					v1=kx-1>=0?buf[(iw*ky+kx-1)<<2|1]>>(8-nbits):0,//W
					v0=kx-0>=0?buf[(iw*ky+kx-0)<<2|1]>>(8-nbits):0;//curr
				int color=v2<<(nbits<<1)|v1<<nbits|v0;

				++htemp[color];
			}
		}
		break;
	case 2://show correlation in space x (CURR, N, NN)
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
				unsigned char
					v2=ky-2>=0?buf[(iw*(ky-2)+kx)<<2|1]>>(8-nbits):0,//NN
					v1=ky-1>=0?buf[(iw*(ky-1)+kx)<<2|1]>>(8-nbits):0,//N
					v0=ky-0>=0?buf[(iw*(ky-0)+kx)<<2|1]>>(8-nbits):0;//curr
				int color=v2<<(nbits<<1)|v1<<nbits|v0;

				++htemp[color];
			}
		}
		break;
	case 3://show correlation in space x (CURR, N, W)
		for(int ky=0;ky<ih;++ky)
		{
			for(int kx=0;kx<iw;++kx)
			{
				unsigned char
					v2=kx-1>=0?buf[(iw* ky   +kx-1)<<2|1]>>(8-nbits):0,//W
					v1=ky-1>=0?buf[(iw*(ky-1)+kx  )<<2|1]>>(8-nbits):0,//N
					v0=        buf[(iw* ky   +kx  )<<2|1]>>(8-nbits)  ;//curr
				int color=v2<<(nbits<<1)|v1<<nbits|v0;

				++htemp[color];
			}
		}
		break;
	}

	//don't calculate csize from downsampled histogram

	unsigned histmin=0;
	unsigned histmax=0;
	for(int k=0;k<hsize;++k)//get min & max
	{
		if(histmin>htemp[k])
			histmin=htemp[k];
		if(histmax<htemp[k])
			histmax=htemp[k];
	}

	ARRAY_ALLOC(char, *hist, 0, hsize, 0, 0);
	unsigned char *h2=hist[0]->data;
	for(int k=0;k<hsize;++k)//normalize
		h2[k]=htemp[k]*255/histmax;
	free(htemp);
}


ArrayHandle bayes_mem[8]={0};
void bayes_estimate(unsigned char *src, int iw, int ih, int x1, int x2, int y1, int y2, int kc)
{
	if(!*bayes_mem)
	{
		for(int kb=7;kb>=0;--kb)
			ARRAY_ALLOC(BayesCounter, bayes_mem[kb], 0, 256LL<<(7-kb), 0, 0);
	}
	int val=1;
	for(int kb=7;kb>=0;--kb)//reset memory
		memfill(bayes_mem[kb]->data, &val, bayes_mem[kb]->count*bayes_mem[kb]->esize, sizeof(int));
	for(int ky=y1;ky<y2;++ky)
	{
		for(int kx=x1;kx<x2;++kx)
		{
			char nb[]=
			{
#define LOAD(X, Y) (unsigned)(kx+(X))<(unsigned)iw&&(unsigned)(ky+(Y))<(unsigned)ih?src[(iw*(ky+(Y))+kx+(X))<<2|kc]-128:0
				LOAD(-2, -2), LOAD(-1, -2), LOAD( 0, -2), LOAD( 1, -2), LOAD( 2, -2),
				LOAD(-2, -1), LOAD(-1, -1), LOAD( 0, -1), LOAD( 1, -1), LOAD( 2, -1),
				LOAD(-2,  0), LOAD(-1,  0),
#undef  LOAD
			};
			int pred=0;
			{
				double fpred=0;
				for(int k=0;k<12;++k)
					fpred+=(double)nb[k]/12;
					//fpred+=customparam_st[12*(kc<<1)+k]*nb[k];
				pred=(int)CLAMP(-128, fpred, 127)+128;
			}
			int context=pred;
			int idx=(iw*ky+kx)<<2|kc;
			unsigned char sym=src[idx];
			for(int kb=7;kb>=0;--kb)
			{
				int bit=sym>>kb&1;
				ArrayHandle mem=bayes_mem[kb];
				BayesCounter *ctr=(BayesCounter*)array_at(&mem, context);
				++ctr->n[bit];
				context|=bit<<(8+7-kb);
			}
		}
	}
}
#endif


typedef struct HybridUintStruct
{
	unsigned short token, nbits;
	unsigned bypass;
} HybridUint;
static void quantize_signed(int val, int expbits, int msb, int lsb, HybridUint *hu)
{
	int token, bypass, nbits;
	val=val<<1^-(val<0);
	if(val<(1<<expbits))
	{
		token=val;//token
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=FLOOR_LOG2((unsigned)val);
		int mantissa=val-(1<<lgv);
		token = (1<<expbits) + (
				(lgv-expbits)<<(msb+lsb)|
				(mantissa>>(lgv-msb))<<lsb|
				(mantissa&((1<<lsb)-1))
			);
		nbits=lgv-(msb+lsb);
		bypass=val>>lsb&((1LL<<nbits)-1);
	}
	hu->token=(unsigned short)token;
	hu->nbits=(unsigned short)nbits;
	hu->bypass=bypass;
}
const char* ec_method_label(EContext ec_method)
{
	const char *label=0;
	switch(ec_method)
	{
#define CASE(L) case L:label=#L;break;
	CASE(ECTX_HIST)
	CASE(ECTX_ABAC)
	CASE(ECTX_ZERO)
	CASE(ECTX_YUV422)
	CASE(ECTX_YUV420)
	CASE(ECTX_QNW)
	CASE(ECTX_MIN_QN_QW)
	CASE(ECTX_MAX_QN_QW)
	CASE(ECTX_MIN_N_W_NW_NE)
	CASE(ECTX_ARGMIN_N_W_NW_NE)
#undef  CASE
	default:
		break;
	}
	return label;
}
typedef enum NBIndexEnum
{
	NB_NNWW, NB_NNW, NB_NN, NB_NNE, NB_NNEE,
	NB_NWW,  NB_NW,  NB_N,  NB_NE,  NB_NEE,
	NB_WW,   NB_W,
} NBIndex;
static void getnb(Image const *src, int kc, int kx, int ky, int *nb)
{
	for(int ky2=-2, idx2=0;ky2<=0;++ky2)
	{
		for(int kx2=-2;kx2<=2;++kx2, ++idx2)
		{
			if(!ky2&&!kx2)
				break;
			if((unsigned)(ky+ky2)<(unsigned)src->ih&&(unsigned)(kx+kx2)<(unsigned)src->iw)
				nb[idx2]=src->data[(src->iw*(ky+ky2)+kx+kx2)<<2|kc];
			else
				nb[idx2]=0;
		}
	}
}
static int getctx_zero(const int *nb, int expbits, int msb, int lsb)
{
	(void)nb;
	(void)expbits;
	(void)msb;
	(void)lsb;
	return 0;
}
static int getctx_QNW(const int *nb, int expbits, int msb, int lsb)
{
	int ctx=(nb[NB_N]+nb[NB_W])>>1;
	HybridUint hu;
	quantize_signed(ctx, expbits, msb, lsb, &hu);
	return hu.token;
}
static int getctx_min_QN_QW(const int *nb, int expbits, int msb, int lsb)
{
	int dy=(nb[NB_N]-nb[NB_NN])>>1, dx=(nb[NB_W]-nb[NB_WW])>>1;
	HybridUint hu1, hu2;
	quantize_signed(dx, expbits, msb, lsb, &hu1);
	quantize_signed(dy, expbits, msb, lsb, &hu2);
	return MAXVAR(hu1.token, hu2.token);
}
static int getctx_max_QN_QW(const int *nb, int expbits, int msb, int lsb)
{
	int dy=(nb[NB_N]-nb[NB_NN])>>1, dx=(nb[NB_W]-nb[NB_WW])>>1;
	HybridUint hu1, hu2;
	quantize_signed(dx, expbits, msb, lsb, &hu1);
	quantize_signed(dy, expbits, msb, lsb, &hu2);
	return MINVAR(hu1.token, hu2.token);
}
static int getctx_min_N_W_NW_NE(const int *nb, int expbits, int msb, int lsb)
{
	int dy=(nb[NB_N]-nb[NB_NN])>>1, dx=(nb[NB_W]-nb[NB_WW])>>1;
	int d135=(nb[NB_NW]-nb[NB_NNWW])>>1, d45=(nb[NB_NE]-nb[NB_NNEE])>>1;
	HybridUint hu[4];
	quantize_signed(dx, expbits, msb, lsb, &hu[0]);
	quantize_signed(dy, expbits, msb, lsb, &hu[1]);
	quantize_signed(d45, expbits, msb, lsb, &hu[2]);
	quantize_signed(d135, expbits, msb, lsb, &hu[3]);
	{
		int ctx=hu->token;
		UPDATE_MIN(ctx, hu[1].token);
		UPDATE_MIN(ctx, hu[2].token);
		UPDATE_MIN(ctx, hu[3].token);
		return ctx;
	}
}
static int getctx_argmin_N_W_NW_NE(const int *nb, int expbits, int msb, int lsb)
{
	int dy=abs(nb[NB_N]-nb[NB_NN]), dx=abs(nb[NB_W]-nb[NB_WW]);
	int d135=abs(nb[NB_NW]-nb[NB_NNWW]), d45=abs(nb[NB_NE]-nb[NB_NNEE]);
	int idx=0, bestval=dx;
	if(bestval>d135)
		bestval=d135, idx=1;
	if(bestval>dy)
		bestval=dy, idx=2;
	if(bestval>d45)
		bestval=d45, idx=3;
	{
		int nb2[]={nb[NB_W], nb[NB_NW], nb[NB_N], nb[NB_NE]};
		int val=nb2[idx];
		HybridUint hu;
		quantize_signed(val, expbits, msb, lsb, &hu);
		return hu.token;
	}
}
void calc_csize_ec(Image const *src, EContext method, int adaptive, int expbits, int msb, int lsb, double *entropy)
{
	int (*const getctx[])(const int *nb, int expbits, int msb, int lsb)=
	{
		getctx_zero,
		getctx_zero,//unused
		getctx_zero,//unused
		getctx_zero,//unused
		getctx_zero,//unused
		getctx_QNW,
		getctx_min_QN_QW,
		getctx_max_QN_QW,
		getctx_min_N_W_NW_NE,
		getctx_argmin_N_W_NW_NE,
	};
	int maxdepth;
	HybridUint hu;
	int cdfsize;
	int *hist, *hsum;
	double bitsizes[4]={0}, bypasssizes[4]={0};
	int res=src->iw*src->ih;

	maxdepth=calc_maxdepth(src, 0);
	quantize_signed(1<<maxdepth, expbits, msb, lsb, &hu);
	cdfsize=hu.token+1;
	hist=(int*)malloc(sizeof(int)*cdfsize*cdfsize);
	hsum=(int*)malloc(sizeof(int)*cdfsize);
	if(!hist||!hsum)
	{
		LOG_ERROR("Alloc error");
		return;
	}
//#define LOAD(X, Y) ((unsigned)(ky+(Y))<src->ih&&(unsigned)(kx+(X))<src->iw?src->data[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
	for(int kc=0;kc<4;++kc)
	{
		int depth=src->depth[kc];
		int constant=1;
		if(!depth)
			continue;
		for(int k=1;k<res;++k)//check for constant value
		{
			if(src->data[k<<2|kc]!=src->data[kc])
			{
				constant=0;
				break;
			}
		}
		if(constant)
		{
			bypasssizes[kc]+=depth;
			continue;
		}
		if(adaptive)
		{
			int fillval=1;
			memfill(hist, &fillval, sizeof(int)*cdfsize*cdfsize, sizeof(int));
			for(int ky=0, idx=0;ky<src->ih;++ky)
			{
				for(int kx=0;kx<src->iw;++kx, ++idx)
				{
					int nb[12], ctx, val, *curr_hist, den;

					getnb(src, kc, kx, ky, nb);
					ctx=getctx[method](nb, expbits, msb, lsb);
					val=src->data[idx<<2|kc];
					quantize_signed(val, expbits, msb, lsb, &hu);

					curr_hist=hist+cdfsize*ctx;
					den=0;
					for(int k=0;k<cdfsize;++k)
						den+=curr_hist[k];
					{
						double p=(double)curr_hist[hu.token]/den;
						double bitsize=-log2(p);
						bitsizes[kc]+=bitsize;
						bypasssizes[kc]+=hu.nbits;
					}

					++curr_hist[hu.token];
					if(curr_hist[hu.token]>=adaptive)
					{
						for(int k=0;k<cdfsize;++k)
							curr_hist[k]=(curr_hist[k]+1)>>1;
					}
				}
			}
			continue;
		}
		memset(hist, 0, sizeof(int)*cdfsize*cdfsize);
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				int nb[12], ctx, val;

				getnb(src, kc, kx, ky, nb);
				ctx=getctx[method](nb, expbits, msb, lsb);

				val=src->data[idx<<2|kc];
				//if(kc==3)
				//	val-=nb[NB_W];
				quantize_signed(val, expbits, msb, lsb, &hu);
				val=hu.token;
				++hist[cdfsize*ctx+val];
			}
		}
		for(int k=0;k<cdfsize;++k)
		{
			int sum=0;
			for(int ks=0;ks<cdfsize;++ks)
				sum+=hist[cdfsize*k+ks];
			hsum[k]=sum;
		}
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				int nb[12], ctx, val;

				getnb(src, kc, kx, ky, nb);
				ctx=getctx[method](nb, expbits, msb, lsb);

				val=src->data[idx<<2|kc];
				//if(kc==3)
				//	val-=nb[NB_W];
				quantize_signed(val, expbits, msb, lsb, &hu);
				val=hu.token;
				{
					int freq=hist[cdfsize*ctx+val], sum=hsum[ctx];
					if(!freq||!sum)
						LOG_ERROR("ZPS");
					{
						double p=(double)freq/sum;
						double bitsize=-log2(p);
						bitsizes[kc]+=bitsize;
						bypasssizes[kc]+=hu.nbits;
					}
				}
			}
		}
	}
//#undef  LOAD
	{
		int nch=(src->src_depth[0]!=0)+(src->src_depth[1]!=0)+(src->src_depth[2]!=0)+(src->src_depth[3]!=0);
		double chubitsize=image_getBMPsize(src)*8/nch;
		for(int kc=0;kc<4;++kc)
		{
			entropy[kc]=(bitsizes[kc]+bypasssizes[kc])*src->src_depth[kc]/chubitsize;
			bitsizes[kc]/=8;
			bypasssizes[kc]/=8;
		}
	}
	free(hist);
	free(hsum);
}


	#define ABAC_DUMMY_EC 48
//	#define ABAC_DUMMY_EC 32

static void calc_qlevels(Image const *src, int kc, int *qdivs, int nqlevels)
{
	int depth=src->depth[kc], nlevels=1<<depth, half=nlevels>>1;
	int *hist=(int*)malloc((nlevels+1LL)*sizeof(int));
	ptrdiff_t res=(ptrdiff_t)src->iw*src->ih;
	int sum=0;
	if(!hist)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(hist, 0, (nlevels+1LL)*sizeof(int));
	for(ptrdiff_t k=0;k<res;++k)
	{
		int val=src->data[k<<2|kc];
		val+=half;
		val&=nlevels-1;
		++hist[val];
	}
	for(int ks=0;ks<nlevels;++ks)
	{
		int freq=hist[ks];
		hist[ks]=sum;
		sum+=freq;
	}
	hist[nlevels]=sum;
	for(int kq=1, ks=0;kq<nqlevels;++kq)
	{
		int div=kq*sum/nqlevels;
		for(;hist[ks]<div;++ks);
		div=(ks-half)<<1;
		qdivs[kq-1]=CLAMP(-half, div, half-1);
	}
#if 0
	for(int kq=1;kq<nqlevels-1;++kq)
	{
		if(qdivs[kq]==qdivs[kq-1])
		{
			if(qdivs[kq-1]<=0)
			{
				for(int kq2=0;kq2<kq;++kq2)
					--qdivs[kq2];
			}
			else
			{
				for(int kq2=kq;kq2<nqlevels-1;++kq2)
					++qdivs[kq2];
			}
		}
	}
#endif
	free(hist);
}
#if 0
static int abac_quantize(int x, int *qdivs, int ndivs)
{
	int sym=0;
#if 0
	int L=0, R=qlevels-1;
	while(R)//binary search		lg(nlevels) memory accesses per symbol
	{
		int floorhalf=R>>1;
		sym=L+floorhalf;
		L+=(R-floorhalf)&-(qdivs[sym]<x);
		R=floorhalf;
	}
#endif
#if 1
	for(;sym<ndivs-1&&x<=qdivs[sym];++sym);//linear search
#endif
	return sym;
}
#endif
#define ABAC_NCTX 4
#define ABAC_NCTRS 1
#define ABAC_QLEVELS 8
#define ABAC_PROBSHIFT 21
#define ABAC_STATSSHIFT 25
#define ABAC_MIXERSHIFT 21
#define MINABS(A, B) (abs(A)<abs(B)?(A):(B))
void calc_csize_abac(Image const *src, double *entropy)
{
	int maxdepth=calc_maxdepth(src, 0);
	//HybridUint hu;
	//quantize_signed(1<<maxdepth, expbits, msb, lsb, &hu);
	//int qlevels=floor_log2_32(hu.token)+2;
	//int qlevels=hu.token+1;
	int qdivs[ABAC_QLEVELS-1]={0};
	int treesize=maxdepth*(maxdepth+1)>>1;
	size_t statssize=sizeof(short[ABAC_QLEVELS*ABAC_NCTX*ABAC_NCTRS])*treesize;
	short *stats=(short*)malloc(statssize);
	int *mixer=(int*)malloc(sizeof(int[ABAC_NCTX*ABAC_NCTRS])*maxdepth);
	double bitsizes[4]={0};
	if(!stats||!mixer)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	for(int kc=0;kc<4;++kc)
	{
#ifdef ABAC_DUMMY_EC
#if ABAC_DUMMY_EC==48
		unsigned long long low=0, range=0xFFFFFFFFFFFF;
#else
		unsigned low=0, range=0xFFFFFFFF;
#endif
#endif
		int depth=src->depth[kc], nlevels=1<<depth;
		//int half=nlevels>>1;
		if(!depth)
			continue;
		calc_qlevels(src, kc, qdivs, ABAC_QLEVELS);
		memset(stats, 0, statssize);
		for(int ks=0;ks<nlevels;++ks)
		{
			int sh=depth+depth+7;
			int weight=nlevels-1-ks;
			weight*=weight;
			for(int kb=0, MSBidx=0;kb<depth;++kb)
			{
				int bit=ks>>(depth-1-kb)&1;
				short *ctr=stats+(kb*(kb+1LL)>>1)+kb-MSBidx;

				int val=*ctr+0x8000;
				val+=(int)((0x10000LL-val)*weight>>sh);
				*ctr=(short)(CLAMP(1, val, 0xFFFF)-0x8000);

				MSBidx+=(!bit)&-(MSBidx==kb);
			}
		}
		memfill(stats+treesize, stats, statssize-treesize*sizeof(short), treesize*sizeof(short));
		//memset(stats, 0, sizeof(short[ABAC_QLEVELS*ABAC_NCTX])*treesize);
		//*stats=0x8000;
		//memfill(stats+1, stats, sizeof(short[ABAC_QLEVELS*ABAC_NCTX])*treesize-sizeof(short), sizeof(short));
		*mixer=0x80000;
		memfill(mixer+1, mixer, sizeof(int[ABAC_NCTX*ABAC_NCTRS])*maxdepth-sizeof(int), sizeof(int));
		for(int ky=0, idx=0;ky<src->ih;++ky)
		{
			for(int kx=0;kx<src->iw;++kx, ++idx)
			{
				int nb[12];
				getnb(src, kc, kx, ky, nb);
				{
					int val, tidx;
					short *curr_stats[ABAC_NCTX]={0}, *p0a[ABAC_NCTX]={0};
					int *curr_mixer;
#if 1
					int preds[]=
					{
						nb[NB_N]>>1,
						nb[NB_W]>>1,
						nb[NB_NW]>>1,
						nb[NB_NE]>>1,
						//nb[NB_N]-nb[NB_W],
						//MINABS(nb[NB_N], nb[NB_W]),
						//MINABS(nb[NB_NW], nb[NB_NE]),
						//MINVAR(abs(nb[NB_WW]), abs(nb[NB_NWW]))*3,
						//MINVAR(abs(nb[NB_NNWW]), abs(nb[NB_NNW]))*3,
						//MINVAR(abs(nb[NB_NN]), abs(nb[NB_NNE]))*3,
						//MINVAR(abs(nb[NB_NNEE]), abs(nb[NB_NEE]))*3,

						//nb[NB_N],
						//nb[NB_W],
						//nb[NB_NW]*3,
						//nb[NB_NE]*3,
						//nb[NB_NEE]*3,
						//nb[NB_NNE]*3,
						//nb[NB_NWW]*3,

						//2*nb[NB_N]-nb[NB_NN],
						//2*nb[NB_W]-nb[NB_WW],
						//nb[NB_N]+nb[NB_W]-nb[NB_NW],
						//nb[NB_N]+nb[NB_W]-(nb[NB_NW]+nb[NB_NE]),
						//nb[NB_N]+nb[NB_W]-nb[NB_NW],
						//(nb[NB_W]+nb[NB_N]+nb[NB_NE]+nb[NB_NEE])>>2,
						//(nb[NB_W]+nb[NB_NE])>>1,

						//abac_quantize(nb[NB_N	], qdivs, ABAC_QLEVELS-1),
						//abac_quantize(nb[NB_W	], qdivs, ABAC_QLEVELS-1),
						//abac_quantize(nb[NB_NW	], qdivs, ABAC_QLEVELS-1),
						//abac_quantize(nb[NB_NE	], qdivs, ABAC_QLEVELS-1),
						//abac_quantize(nb[NB_NN	], qdivs, ABAC_QLEVELS-1),
						//abac_quantize(nb[NB_WW	], qdivs, ABAC_QLEVELS-1),
						//abac_quantize(nb[NB_NEE	], qdivs, ABAC_QLEVELS-1),
					};
#endif
					//ctx=ABAC_QLEVELS*ctx+abac_quantize(nb[NB_W], qdivs, ABAC_QLEVELS-1);
					//int ctx=getctx[method](nb, expbits, msb, lsb);
					//ctx=floor_log2(ctx)+1;
					//if(ctx>=qlevels)
					//	LOG_ERROR("");
					//if(kx==5&&ky==5)
					//	printf("");

					val=src->data[idx<<2|kc];
					val=val<<1^-(val<0);
				
					tidx=0;
					for(int kt=0;kt<ABAC_NCTX;++kt)
					{
						int v2=preds[kt];
						v2=THREEWAY(v2, 0)+1;
						//v2=(v2>1)-(v2<-1)+1;//less efficient
						curr_stats[kt]=stats+treesize*ABAC_NCTRS*v2;
					}
						//curr_stats[kt]=stats+treesize*abac_quantize(CLAMP(-half, preds[kt], half-1), qdivs, ABAC_QLEVELS-1);
					curr_mixer=mixer;
					for(int kb=0;kb<depth;++kb)
					{
						int bit;

						int idx2=ABAC_NCTRS*((kb*(kb+1LL)>>1)+kb-tidx);
						long long p0=0;
						for(int kt=0;kt<ABAC_NCTX;++kt)
						{
							for(int k2=0;k2<ABAC_NCTRS;++k2)
							{
								p0a[kt]=curr_stats[kt]+idx2;
								p0+=(long long)curr_mixer[ABAC_NCTRS*kt+k2]*p0a[kt][k2];
							}
						}
						p0+=p0>>3;
						p0+=(1LL<<ABAC_PROBSHIFT>>1)-1;
						p0>>=ABAC_PROBSHIFT;
						//p0/=ABAC_NCTX;
						p0+=0x8000;
						p0=CLAMP(1, p0, 0xFFFF);
						//int p0=curr_stats[idx2];

						bit=val>>(depth-1-kb)&1;
					
#ifdef ABAC_DUMMY_EC
#if ABAC_DUMMY_EC==48
						{
							unsigned long long r2=range*p0>>16;
							low+=r2&-bit;
							range=bit?range-r2:r2-1;
						}
						while(range<0x10000)
						{
							bitsizes[kc]+=16;

							range<<=16;
							low<<=16;
							range|=0xFFFF;
							low&=0xFFFFFFFFFFFF;
							{
								unsigned long long rmax=0xFFFFFFFFFFFF-low;
								if(range>rmax)//clamp range
									range=rmax;
							}
						}
#else
						{
#ifdef _MSC_VER
							unsigned r2=(unsigned)(__emulu(range, p0)>>16);
#else
							unsigned r2=(int)((unsigned long long)range*p0>>16);
#endif
							low+=r2&-bit;
							range=bit?range-r2:r2-1;
							while(range<0x10000)
							{
								bitsizes[kc]+=16;

								range<<=16;
								low<<=16;
								range|=0xFFFF;
								unsigned rmax=~low;
								if(range>rmax)//clamp range
									range=rmax;
							}
						}
#endif
#else
						{
							int prob=(int)p0;
							prob^=-bit;
							prob+=0x10001&-bit;
							//if(!prob)
							//	LOG_ERROR("");
							bitsizes[kc]-=log2((double)prob/0x10000);
						}
#endif
						{
							//L = (1/2)sq(p0_collapse-p0)		p0 = sum: m[k]*s[k]
							//dL/ds[k] = (p0-p0_collapse)*m[k]
							//dL/dm[k] = (p0-p0_collapse)*s[k]
							int collapse=!bit<<16;
							int diff=collapse-(int)p0;
							for(int kt=0;kt<ABAC_NCTX;++kt)
							{
								for(int k2=0;k2<ABAC_NCTRS;++k2)
								{
									int temp=curr_mixer[ABAC_NCTRS*kt+k2];
									curr_mixer[ABAC_NCTRS*kt+k2]+=(int)((((long long)diff*p0a[kt][k2])+(1LL<<ABAC_MIXERSHIFT>>1)-1)>>ABAC_MIXERSHIFT);
									temp=p0a[kt][k2]+(int)((((long long)diff*temp)+(1LL<<(ABAC_STATSSHIFT-k2)>>1)-1)>>(ABAC_STATSSHIFT-k2));
									p0a[kt][k2]=(short)CLAMP(-0x7FFF, temp, 0x7FFF);
									//*p0a[kt]+=(collapse-*p0a[kt])>>8;
								}
							}
							//p0+=((!bit<<16)-p0)>>8;
							//curr_stats[idx2]=CLAMP(1, p0, 0xFFFF);

							curr_mixer+=ABAC_NCTX*ABAC_NCTRS;
							tidx+=(!bit)&-(tidx==kb);
						}
					}
				}
			}
		}
#ifdef ABAC_DUMMY_EC
#if ABAC_DUMMY_EC==48
		{
			unsigned long long code=low+range;
			int nzeros=FLOOR_LOG2(low^code)-1;
			if(nzeros<0)
				nzeros=0;
			code&=~((1LL<<nzeros)-1);
			bitsizes[kc]+=LSB_IDX_64(code);
		}
#else
		{
			unsigned code=low+range;//coming from renorm, cannot overflow
			int nzeros=FLOOR_LOG2(low^code)-1;
			if(nzeros<0)
				nzeros=0;
			code&=~((1LL<<nzeros)-1);
			bitsizes[kc]+=LSB_IDX_32(code);
		}
#endif
#endif
	}
#if 0
	const size_t strsize=0x10000;
	char *str=(char*)malloc(strsize);
	if(!str)
	{
		LOG_ERROR("Alloc error");
		return;
	}
	memset(str, 0, strsize);
	int printed=0;
	int *curr_mixer=mixer;
	for(int kb=0;kb<maxdepth;++kb)
	{
		for(int kt=0;kt<ABAC_NCTX;++kt)
		{
			int val=curr_mixer[kt];
			//printed+=snprintf(str+printed, strsize-printed-1, " %c0x%08X,", val<0?'-':' ', abs(val));
			printed+=snprintf(str+printed, strsize-printed-1, " %8d,", val);
		}
		printed+=snprintf(str+printed, strsize-printed-1, "\n");
		curr_mixer+=ABAC_NCTX;
	}
	copy_to_clipboard(str, printed);
	messagebox(MBOX_OK, "Copied to clipboard", "Paste into a text editor");
	free(str);
#endif
	{
		int nch=(src->src_depth[0]!=0)+(src->src_depth[1]!=0)+(src->src_depth[2]!=0)+(src->src_depth[3]!=0);
		double chubitsize=image_getBMPsize(src)*8/nch;
		for(int kc=0;kc<4;++kc)
			entropy[kc]=bitsizes[kc]*src->src_depth[kc]/chubitsize;
	}
	free(stats);
	free(mixer);
}