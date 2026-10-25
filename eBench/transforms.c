#include"ebench.h"
#include<stdint.h>
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
void calc_histogram(const int *buf, int iw, int ih, int kc, int x1, int x2, int y1, int y2, int depth, int *hist, int *hist8)
{
	int nlevels, shift, ky, kx;

	(void)ih;
	nlevels=1<<depth;
	shift=MAXVAR(0, depth-8);
	memset(hist, 0, nlevels*sizeof(int));
	if(hist8)
		memset(hist8, 0, 256*sizeof(int));
	for(ky=y1;ky<y2;++ky)
	{
		for(kx=x1;kx<x2;++kx)
		{
			int sym=buf[(iw*ky+kx)<<2|kc]+(nlevels>>1);
			sym&=nlevels-1;
			++hist[sym];
			if(hist8)
				++hist8[sym>>shift];
		}
	}
}
double calc_entropy(const int *hist, int nlevels, int sum)/*nlevels = 1<<current_depth*/
{
	double entropy, norm;
	int k;

	entropy=0;
	if(sum<0)
	{
		sum=0;
		for(k=0;k<nlevels;++k)
			sum+=hist[k];
	}
	if(!sum)/*degenerate distribution*/
		return 0;
	norm=1./sum;
	for(k=0;k<nlevels;++k)
	{
		int f=hist[k];
		if(f)
		{
			double p=(double)f*norm;
			entropy-=p*log(p);
		}
	}
	return entropy*(1./M_LN2);
}
int calc_maxdepth(Image const *image, const int *inflation)
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
	intptr_t res, k;
	int kc;

	res=(intptr_t)iw*ih;
	for(kc=0;kc<3;++kc)
	{
		int vmin=image[0<<2|kc], vmax=image[0<<2|kc];
		for(k=1;k<res;++k)
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
			depths[kc]=31^LZCNT32(nlevels);
			depths[kc]+=1<<depths[kc]<nlevels;
		}
		depths[kc]=MAXVAR(depths[kc], src_depths[kc]);
	}
}


/*color transforms*/
void colortransform_YCoCg_R(Image *image, int fwd)
{
	intptr_t k, len;
	char temp;

	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];

			r-=b;		/*co = r-b			diff(r, b)			*/
			b+=r>>1;	/*(r+b)/2							*/
			g-=b;		/*cg = g-(r+b)/2		diff(g, av(r, b))		*/
			b+=g>>1;	/*Y  = (r+b)/2 + (g-(r+b)/2)/2 = r/4+g/2+b/4	av(g, av(r, b))	*/

			((int32_t*)HDR_DATA(image))[k  ]=b;/*Y	*/
			((int32_t*)HDR_DATA(image))[k|1]=r;/*Co	*/
			((int32_t*)HDR_DATA(image))[k|2]=g;/*Cg	*/
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
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Co=((int32_t*)HDR_DATA(image))[k+1];
			int Cg=((int32_t*)HDR_DATA(image))[k+2];
			
			Y-=Cg>>1;
			Cg+=Y;
			Y-=Co>>1;
			Co+=Y;

			((int32_t*)HDR_DATA(image))[k  ]=Co;	/*r*/
			((int32_t*)HDR_DATA(image))[k|1]=Cg;	/*g*/
			((int32_t*)HDR_DATA(image))[k|2]=Y;	/*b*/
		}
	}
}
void colortransform_YCbCr_R_v1(Image *image, int fwd)
{
	intptr_t k, len;
	char temp;

	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];

			r-=g;		/*diff(r, g)            [ 1      -1      0  ].RGB	Cr	*/
			g+=r>>1;
			b-=g;		/*diff(b, av(r, g))     [-1/2    -1/2    1  ].RGB	Cb	*/
			g+=b>>1;	/*av(b, av(r, g))       [ 1/4     1/4    1/2].RGB	Y	*/

			((int32_t*)HDR_DATA(image))[k+0]=g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k+1]=b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k+2]=r;/*Cr	*/
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
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Cb=((int32_t*)HDR_DATA(image))[k+1];
			int Cr=((int32_t*)HDR_DATA(image))[k+2];
			
			Y-=Cb>>1;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			((int32_t*)HDR_DATA(image))[k+0]=Cr;
			((int32_t*)HDR_DATA(image))[k+1]=Y;
			((int32_t*)HDR_DATA(image))[k+2]=Cb;
		}
	}
}
void colortransform_YCbCr_R_v2(Image *image, int fwd)
{
	intptr_t k, len;
	char temp;

	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];
			
			r-=g;		/*Cr =	[1	-1	0].RGB		*/
			g+=r>>1;	/*	[1/2	1/2	0]		*/
			b-=g;		/*Cb =	[-1/2	-1/2	1]		*/
			g+=(2*b-r)>>3;	/*Y  =	[1/4	1/2	1/4]	v2	*/

			((int32_t*)HDR_DATA(image))[k  ]=g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k|1]=b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k|2]=r;/*Cr	*/
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
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Cb=((int32_t*)HDR_DATA(image))[k+1];
			int Cr=((int32_t*)HDR_DATA(image))[k+2];
			
			Y-=(2*Cb-Cr)>>3;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			((int32_t*)HDR_DATA(image))[k+0]=Cr;
			((int32_t*)HDR_DATA(image))[k+1]=Y;
			((int32_t*)HDR_DATA(image))[k+2]=Cb;
		}
	}
}
void colortransform_YCbCr_R_v3(Image *image, int fwd)
{
	intptr_t k, len;
	char temp;

	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];
			
			r-=g;		/*Cr =	[1	-1	0].RGB		*/
			g+=r>>1;	/*	[1/2	1/2	0]		*/
			b-=g;		/*Cb =	[-1/2	-1/2	1]		*/
			g+=(2*b+r)>>3;	/*Y  =	[1/2	1/4	1/4]	v3	*/

			((int32_t*)HDR_DATA(image))[k+0]=g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k+1]=b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k+2]=r;/*Cr	*/
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
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Cb=((int32_t*)HDR_DATA(image))[k+1];
			int Cr=((int32_t*)HDR_DATA(image))[k+2];
			
			Y-=(2*Cb+Cr)>>3;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			((int32_t*)HDR_DATA(image))[k+0]=Cr;
			((int32_t*)HDR_DATA(image))[k+1]=Y;
			((int32_t*)HDR_DATA(image))[k+2]=Cb;
		}
	}
}
void colortransform_YCbCr_R_v4(Image *image, int fwd)
{
	intptr_t k, len;
	char temp;

	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];
			
			r-=g;		/*Cr =	[1	-1	0].RGB		*/
			g+=r>>1;	/*	[1/2	1/2	0]		*/
			b-=g;		/*Cb =	[-1/2	-1/2	1]		*/
			g+=b/3;		/*Y  =	[1/3	1/3	1/3]	v4	*/

			((int32_t*)HDR_DATA(image))[k+0]=g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k+1]=b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k+2]=r;/*Cr	*/
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
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Cb=((int32_t*)HDR_DATA(image))[k+1];
			int Cr=((int32_t*)HDR_DATA(image))[k+2];
			
			Y-=Cb/3;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			((int32_t*)HDR_DATA(image))[k  ]=Cr;
			((int32_t*)HDR_DATA(image))[k|1]=Y;
			((int32_t*)HDR_DATA(image))[k|2]=Cb;
		}
	}
}
void colortransform_YCbCr_R_v5(Image *image, int fwd)
{
	intptr_t k, len;
	char temp;

	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];
			
			r-=g;		/*Cr =	[1	-1	0].RGB		*/
			g+=r>>1;	/*	[1/2	1/2	0]		*/
			b-=g;		/*Cb =	[-1/2	-1/2	1]		*/
			g+=6*b>>4;	/*Y  =	[5/16	5/16	6/16]	v5	*/

			((int32_t*)HDR_DATA(image))[k+0]=g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k+1]=b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k+2]=r;/*Cr	*/
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Cb=((int32_t*)HDR_DATA(image))[k+1];
			int Cr=((int32_t*)HDR_DATA(image))[k+2];
			
			Y-=6*Cb>>4;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			((int32_t*)HDR_DATA(image))[k+0]=Cr;
			((int32_t*)HDR_DATA(image))[k+1]=Y;
			((int32_t*)HDR_DATA(image))[k+2]=Cb;
		}
		ROTATE3(image->depth[2], image->depth[1], image->depth[1], temp);
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
	}
}
void colortransform_YCbCr_R_v6(Image *image, int fwd)
{
	intptr_t k, len;
	char temp;

	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];
			
			r-=g;		/*Cr =	[1	-1	0].RGB		*/
			g+=r>>1;	/*	[1/2	1/2	0]		*/
			b-=g;		/*Cb =	[-1/2	-1/2	1]		*/
			g+=14*b>>5;	/*Y  =	[9/32	9/32	14/32]	v6	*/

			((int32_t*)HDR_DATA(image))[k+0]=g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k+1]=b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k+2]=r;/*Cr	*/
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Cb=((int32_t*)HDR_DATA(image))[k+1];
			int Cr=((int32_t*)HDR_DATA(image))[k+2];
			
			Y-=14*Cb>>5;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			((int32_t*)HDR_DATA(image))[k+0]=Cr;
			((int32_t*)HDR_DATA(image))[k+1]=Y;
			((int32_t*)HDR_DATA(image))[k+2]=Cb;
		}
		ROTATE3(image->depth[2], image->depth[1], image->depth[1], temp);
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
	}
}
void colortransform_YCbCr_R_v7(Image *image, int fwd)
{
	intptr_t k, len;
	char temp;

	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];
			
			r-=g;		/*Cr =	[1	-1	0].RGB		*/
			g+=r>>1;	/*	[1/2	1/2	0]		*/
			b-=g;		/*Cb =	[-1/2	-1/2	1]		*/
			g+=(10*b-r)>>5;	/*Y  =	[5/16	 6/16	5/16]	v7	*/

			((int32_t*)HDR_DATA(image))[k+0]=g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k+1]=b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k+2]=r;/*Cr	*/
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Cb=((int32_t*)HDR_DATA(image))[k+1];
			int Cr=((int32_t*)HDR_DATA(image))[k+2];
			
			Y-=(10*Cb-Cr)>>5;
			Cb+=Y;
			Y-=Cr>>1;
			Cr+=Y;

			((int32_t*)HDR_DATA(image))[k+0]=Cr;
			((int32_t*)HDR_DATA(image))[k+1]=Y;
			((int32_t*)HDR_DATA(image))[k+2]=Cb;
		}
		ROTATE3(image->depth[2], image->depth[1], image->depth[1], temp);
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
	}
}
void colortransform_Pei09(Image *image, int fwd)/*Pei09 RCT*/
{
	intptr_t k, len;
	char temp;
	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];

			b-=(87*r+169*g+128)>>8;	/*Cb = [-87/256  -169/256  1]	*/
			r-=g;			/*Cr = [1  -1  0].RGB		*/
			g+=(86*r+29*b+128)>>8;	/*Y  = [19493/65536  38619/65536  29/256]	g+86/256*(r-g)+29/256*(b-87/256*r-169/256*g) = 19493/65536*r + 38619/65536*g + 29/256*b	*/

			((int32_t*)HDR_DATA(image))[k+0]=g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k+1]=b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k+2]=r;/*Cr	*/
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
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Cb=((int32_t*)HDR_DATA(image))[k+1];
			int Cr=((int32_t*)HDR_DATA(image))[k+2];
			
			Y-=(86*Cr+29*Cb+128)>>8;
			Cr+=Y;
			Cb+=(87*Cr+169*Y+128)>>8;

			((int32_t*)HDR_DATA(image))[k+0]=Cr;
			((int32_t*)HDR_DATA(image))[k+1]=Y;
			((int32_t*)HDR_DATA(image))[k+2]=Cb;
		}
	}
}
void colortransform_optsub(Image *image, int fwd)
{
	intptr_t k, len;
	const uint8_t *combination;
	int yidx, uidx, vidx, uc0, vc0, vc1, yuv[3];

	if(fwd)
		image->rct=crct_analysis(image);
	combination=rct_combinations[image->rct];
	yidx=combination[II_PERM_Y];
	uidx=combination[II_PERM_U];
	vidx=combination[II_PERM_V];
	uc0=combination[II_COEFF_U_SUB_Y];
	vc0=combination[II_COEFF_V_SUB_Y];
	vc1=combination[II_COEFF_V_SUB_U];
	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			yuv[0]=((int32_t*)HDR_DATA(image))[k+yidx];
			yuv[1]=((int32_t*)HDR_DATA(image))[k+uidx];
			yuv[2]=((int32_t*)HDR_DATA(image))[k+vidx];

			yuv[2]-=(vc0*yuv[0]+vc1*yuv[1])>>2;
			yuv[1]-=uc0*yuv[0]>>2;

			((int32_t*)HDR_DATA(image))[k+0]=yuv[0];
			((int32_t*)HDR_DATA(image))[k+1]=yuv[1];
			((int32_t*)HDR_DATA(image))[k+2]=yuv[2];
		}
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		int halfs[3];

		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		halfs[0]=1<<image->depth[0]>>1;
		halfs[1]=1<<image->depth[1]>>1;
		halfs[2]=1<<image->depth[2]>>1;
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			yuv[0]=((int32_t*)HDR_DATA(image))[k+0];
			yuv[1]=((int32_t*)HDR_DATA(image))[k+1];
			yuv[2]=((int32_t*)HDR_DATA(image))[k+2];

			yuv[1]+=uc0*yuv[0]>>2;
			yuv[2]+=(vc0*yuv[0]+vc1*yuv[1])>>2;
			CLAMP2(yuv[0], -halfs[0], halfs[0]-1);
			CLAMP2(yuv[1], -halfs[1], halfs[1]-1);
			CLAMP2(yuv[2], -halfs[2], halfs[2]-1);
			((int32_t*)HDR_DATA(image))[k+yidx]=yuv[0];
			((int32_t*)HDR_DATA(image))[k+uidx]=yuv[1];
			((int32_t*)HDR_DATA(image))[k+vidx]=yuv[2];
		}
	}
}
void colortransform_JPEG2000(Image *image, int fwd)
{
	intptr_t k, len;
	char temp;

	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];
			
			r-=g;       /*r-g				[1     -1     0  ].RGB*/
			b-=g;       /*b-g				[0     -1     1  ].RGB*/
			g+=(r+b)>>2;/*g+(r-g+b-g)/4 = r/4+g/2+b/4	[1/4    1/2   1/4].RGB*/

			((int32_t*)HDR_DATA(image))[k+0]=g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k+1]=b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k+2]=r;/*Cr	*/
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
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Cb=((int32_t*)HDR_DATA(image))[k+1];
			int Cr=((int32_t*)HDR_DATA(image))[k+2];
			
			Y-=(Cr+Cb)>>2;
			Cb+=Y;
			Cr+=Y;

			((int32_t*)HDR_DATA(image))[k+0]=Cr;
			((int32_t*)HDR_DATA(image))[k+1]=Y;
			((int32_t*)HDR_DATA(image))[k+2]=Cb;
		}
	}
}
void colortransform_JPEG2000_MA(Image *image, int fwd)
{
	intptr_t k, len;
	char temp;

	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];
			
			r-=g;       /*r-g				[1     -1     0  ].RGB*/
			r<<=32-image->depth[0];
			r>>=32-image->depth[0];
			b-=g;       /*b-g				[0     -1     1  ].RGB*/
			b<<=32-image->depth[2];
			b>>=32-image->depth[2];
			g+=(r+b)>>2;/*g+(r-g+b-g)/4 = r/4+g/2+b/4	[1/4    1/2   1/4].RGB*/
			g<<=32-image->depth[1];
			g>>=32-image->depth[1];

			((int32_t*)HDR_DATA(image))[k+0]=g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k+1]=b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k+2]=r;/*Cr	*/
		}
		ROTATE3(image->depth[0], image->depth[1], image->depth[2], temp);
	}
	else
	{
		ROTATE3(image->depth[2], image->depth[1], image->depth[0], temp);
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Cb=((int32_t*)HDR_DATA(image))[k+1];
			int Cr=((int32_t*)HDR_DATA(image))[k+2];
			
			Y-=(Cr+Cb)>>2;
			Y<<=32-image->depth[1];
			Y>>=32-image->depth[1];
			Cb+=Y;
			Cb<<=32-image->depth[2];
			Cb>>=32-image->depth[2];
			Cr+=Y;
			Cr<<=32-image->depth[0];
			Cr>>=32-image->depth[0];

			((int32_t*)HDR_DATA(image))[k+0]=Cr;
			((int32_t*)HDR_DATA(image))[k+1]=Y;
			((int32_t*)HDR_DATA(image))[k+2]=Cb;
		}
	}
}
void colortransform_NBLI(Image *image, int fwd)
{
	intptr_t k, len;
	char temp;

	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];
			
			r-=g;       /*r-g				[1     -1     0  ].RGB*/
			b-=g;       /*b-g				[0     -1     1  ].RGB*/
			g+=(r+b)>>2;/*g+(r-g+b-g)/4 = r/4+g/2+b/4	[1/4    1/2   1/4].RGB*/
			b-=r>>2;    /*b-g-(r-g)/4 = b-(r+3*g)/4		[1     -3/4  -1/4].RGB*/

			((int32_t*)HDR_DATA(image))[k+0]=g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k+1]=b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k+2]=r;/*Cr	*/
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
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Cb=((int32_t*)HDR_DATA(image))[k+1];
			int Cr=((int32_t*)HDR_DATA(image))[k+2];
			
			Cb+=Cr>>2;
			Y-=(Cr+Cb)>>2;
			Cb+=Y;
			Cr+=Y;

			((int32_t*)HDR_DATA(image))[k+0]=Cr;
			((int32_t*)HDR_DATA(image))[k+1]=Y;
			((int32_t*)HDR_DATA(image))[k+2]=Cb;
		}
	}
}
void colortransform_subtractgreen(Image *image, int fwd)
{
	intptr_t k, len;
	char temp;

	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int r=((int32_t*)HDR_DATA(image))[k+0];
			int g=((int32_t*)HDR_DATA(image))[k+1];
			int b=((int32_t*)HDR_DATA(image))[k+2];
			
			r-=g;
			b-=g;

			((int32_t*)HDR_DATA(image))[k+0]=g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k+1]=b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k+2]=r;/*Cr	*/
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
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int Y =((int32_t*)HDR_DATA(image))[k+0];
			int Cb=((int32_t*)HDR_DATA(image))[k+1];
			int Cr=((int32_t*)HDR_DATA(image))[k+2];
			
			Cb+=Y;
			Cr+=Y;

			((int32_t*)HDR_DATA(image))[k+0]=Cr;
			((int32_t*)HDR_DATA(image))[k+1]=Y;
			((int32_t*)HDR_DATA(image))[k+2]=Cb;
		}
	}
}

void colortransform_lossy_YCbCr(Image *image, int fwd)
{
	intptr_t k, res;

	if(fwd)
	{
		for(k=0, res=image->iw*image->ih;k<res;++k)
		{
			int
				r=((int32_t*)HDR_DATA(image))[k*4+0],
				g=((int32_t*)HDR_DATA(image))[k*4+1],
				b=((int32_t*)HDR_DATA(image))[k*4+2];

			double
				y=0.299*r+0.587*g+0.114*b,
				u=-0.168736*r-0.331264*g+0.5*b,
				v=+0.5*r-0.418688*g-0.081312*b;

			CLAMP2(y, -128, 127);
			CLAMP2(u, -128, 127);
			CLAMP2(v, -128, 127);
			((int32_t*)HDR_DATA(image))[k*4+0]=(int)y;
			((int32_t*)HDR_DATA(image))[k*4+1]=(int)u;
			((int32_t*)HDR_DATA(image))[k*4+2]=(int)v;
		}
	}
	else
	{
		for(k=0, res=image->iw*image->ih;k<res;++k)
		{
			int
				Y =((int32_t*)HDR_DATA(image))[k*4+0],
				Cb=((int32_t*)HDR_DATA(image))[k*4+1],
				Cr=((int32_t*)HDR_DATA(image))[k*4+2];

			double
				r=Y+1.402*Cr,
				g=Y-0.344136*Cb-0.714136*Cr,
				b=Y+1.772*Cb;

			CLAMP2(r, -128, 127);
			CLAMP2(g, -128, 127);
			CLAMP2(b, -128, 127);
			((int32_t*)HDR_DATA(image))[k*4+0]=(int)r;
			((int32_t*)HDR_DATA(image))[k*4+1]=(int)g;
			((int32_t*)HDR_DATA(image))[k*4+2]=(int)b;
		}
	}
}
void colortransform_lossy_XYB(Image *image, int fwd)
{
	const double bias=0.00379307325527544933;
	double cbrt_bias=pow(bias, 1./3);
	int nlevels[3];
	intptr_t k, res;
	char temp;
	if(fwd)
	{
		nlevels[0]=1<<image->depth[0];
		nlevels[1]=1<<image->depth[1];
		nlevels[2]=1<<image->depth[2];
		for(k=0, res=(intptr_t)image->iw*image->ih;k<res;++k)
		{
			double
				r=(double)(((int32_t*)HDR_DATA(image))[k*4+0]+(nlevels[0]>>1))/(nlevels[0]-1),/*input: RGB in [0, 1]*/
				g=(double)(((int32_t*)HDR_DATA(image))[k*4+1]+(nlevels[1]>>1))/(nlevels[1]-1),
				b=(double)(((int32_t*)HDR_DATA(image))[k*4+2]+(nlevels[2]>>1))/(nlevels[2]-1);

			double
				L=0.3*r+0.622*g+0.078*b,
				M=0.23*r+0.692*g+0.078*b,
				S=0.24342268924547819*r+0.20476744424496821*g+0.55180986650955360*b;
			L=pow(L+bias, 1./3)-cbrt_bias;
			M=pow(M+bias, 1./3)-cbrt_bias;
			S=pow(S+bias, 1./3)-cbrt_bias;
			r=(L-M)*0.5;	/*X	*/
			g=(L+M)*0.5;	/*Y	*/
			b=S-g;		/*B-Y	*/

			r*=(int)(4096.*((int64_t)nlevels[0]<<(nlevels[0]<0x1000000))/512);/*customparam_ct[0]*1000*/
			g*=(int)(144.*nlevels[1]/256);
			b*=(int)(1024.*((int64_t)nlevels[2]<<(nlevels[2]<0x1000000))/512);
			r=CLAMP(-nlevels[0], r, nlevels[0]-1);
			g=CLAMP(-(nlevels[1]>>1), g, (nlevels[1]>>1)-1);
			b=CLAMP(-nlevels[2], b, nlevels[2]-1);

			((int32_t*)HDR_DATA(image))[k*4+0]=(int)g;/*Y	*/
			((int32_t*)HDR_DATA(image))[k*4+1]=(int)b;/*Cb	*/
			((int32_t*)HDR_DATA(image))[k*4+2]=(int)r;/*Cr	*/
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
		for(k=0, res=(intptr_t)image->iw*image->ih;k<res;++k)
		{
			double
				Y =(double)((int32_t*)HDR_DATA(image))[k*4+0]/(144.*nlevels[1]/256),
				Cb=(double)((int32_t*)HDR_DATA(image))[k*4+1]/(1024.*((int64_t)nlevels[2]<<(nlevels[2]<0x1000000))/512),
				Cr=(double)((int32_t*)HDR_DATA(image))[k*4+2]/(4096.*((int64_t)nlevels[0]<<(nlevels[0]<0x1000000))/512);

			double
				L=Y+Cr,
				M=Y-Cr,/*g is the sum*/
				S=Cb+Y;
			L+=cbrt_bias, L*=L*L, L-=bias;
			M+=cbrt_bias, M*=M*M, M-=bias;
			S+=cbrt_bias, S*=S*S, S-=bias;
			{
				double
					r=11.03160*L-9.86694*M-0.164623*S,
					g=-3.25415*L+4.41877*M-0.164623*S,
					b=-3.65885*L+2.71292*M+1.945930*S;

				((int32_t*)HDR_DATA(image))[k*4+0]=(int)((r-0.5)*nlevels[0]);/*output: RGB in [0, 1]*/
				((int32_t*)HDR_DATA(image))[k*4+1]=(int)((g-0.5)*nlevels[1]);
				((int32_t*)HDR_DATA(image))[k*4+2]=(int)((b-0.5)*nlevels[2]);
			}
		}
	}
}

int16_t rct_custom_params[RCT_CUSTOM_NPARAMS]={0};
void rct_custom_unpackpermutation(int16_t p, uint8_t *permutation)
{
	int temp=0;

	switch(p)
	{
	case 0:temp=0x020100;break;/*rgb*/
	case 1:temp=0x020001;break;/*grb*/
	case 2:temp=0x000201;break;/*gbr*/
	case 3:temp=0x010200;break;/*rbg*/
	case 4:temp=0x000102;break;/*bgr*/
	case 5:temp=0x010002;break;/*brg*/
	default:
		CRASH("Invalid RGB permutation");
		return;
	}
	memcpy(permutation, &temp, sizeof(char[3]));
}
void rct_custom(Image *image, int fwd, const int16_t *params)/*4 params	fixed 15.16*/
{
	int temp;
	uint8_t permutation[4]={0};
	intptr_t k, len;

	rct_custom_unpackpermutation(params[8], permutation);
	if(fwd)
	{
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int comp[3];

			comp[0]=((int32_t*)HDR_DATA(image))[k+permutation[0]];
			comp[1]=((int32_t*)HDR_DATA(image))[k+permutation[1]];
			comp[2]=((int32_t*)HDR_DATA(image))[k+permutation[2]];
			temp=params[0]*comp[1]+params[1]*comp[2], comp[0]+=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			temp=params[2]*comp[0]+params[3]*comp[2], comp[1]+=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			temp=params[4]*comp[0]+params[5]*comp[1], comp[2]+=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			temp=params[6]*comp[0]+params[7]*comp[2], comp[1]+=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			memcpy(((int32_t*)HDR_DATA(image))+k, comp, sizeof(comp));
		}
		image->depth[0]+=image->depth[0]<24;
		image->depth[1]+=image->depth[1]<24;
		image->depth[2]+=image->depth[2]<24;
	}
	else
	{
		image->depth[0]-=image->depth[0]>image->src_depth[0];
		image->depth[1]-=image->depth[1]>image->src_depth[1];
		image->depth[2]-=image->depth[2]>image->src_depth[2];
		for(k=0, len=(intptr_t)image->iw*image->ih*4;k<len;k+=4)
		{
			int comp[3];

			memcpy(comp, ((int32_t*)HDR_DATA(image))+k, sizeof(comp));
			temp=params[6]*comp[0]+params[7]*comp[2], comp[1]-=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			temp=params[4]*comp[0]+params[5]*comp[1], comp[2]-=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			temp=params[2]*comp[0]+params[3]*comp[2], comp[1]-=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			temp=params[0]*comp[1]+params[1]*comp[2], comp[0]-=(temp+(((1<<RCT_CUSTOM_PARAMBITS)-1)&-(temp<0)))>>RCT_CUSTOM_PARAMBITS;
			((int32_t*)HDR_DATA(image))[k+permutation[0]]=comp[0];
			((int32_t*)HDR_DATA(image))[k+permutation[1]]=comp[1];
			((int32_t*)HDR_DATA(image))[k+permutation[2]]=comp[2];
		}
	}
}
static void rct_custom_calcloss(Image const *src, Image **pdst, int *hist, double *loss)
{
	Image *dst;
	double entropy;
	intptr_t res;
	int depths[3], kc;

	dst=*pdst;
	res=(intptr_t)src->iw*src->ih;
	depths[0]=src->depth[1];
	depths[1]=src->depth[2]+1;
	depths[2]=src->depth[0]+1;
	memcpy(HDR_DATA(dst), HDR_DATC(src), res*sizeof(int[4]));
	apply_selected_transforms(pdst, 0, 1, 1);
	dst=*pdst;
	for(kc=0;kc<3;++kc)
	{
		calc_histogram((int*)HDR_DATA(dst), dst->iw, dst->ih, kc, 0, dst->iw, 0, dst->ih, depths[kc], hist, 0);
		entropy=calc_entropy(hist, 1<<depths[kc], (int)res);
		loss[kc]=entropy/src->src_depth[(kc+1)%3];
	}
	loss[3]=(loss[0]+loss[1]+loss[2])/3;
}
enum
{
	RCT_CUSTOM_NITER=200,
	RCT_CUSTOM_DELTAGROUP=2
};
void rct_custom_optimize(Image const *image, int16_t *params)
{
	static int call_idx=0;
	const int inflation[]={1, 0, 1};/*original image is expected, RCT permutes depths {r+1, g, b+1} -> {Y:g:+0, Cb:b:+1, Cr:r:+1}*/
	double loss_bestsofar[4], loss_prev[4], loss_curr[4];
	int maxdepth, maxlevels, *hist, it, watchdog, k, loud;
	int16_t params2[RCT_CUSTOM_NPARAMS];
	int shakethreshold=RCT_CUSTOM_NPARAMS;
	Image *im2;

	++call_idx;
	if(call_idx==1)
		DisableProcessWindowsGhosting();
#define CALC_LOSS(L) rct_custom_calcloss(image, &im2, hist, L)
#ifndef _DEBUG
	srand((uint32_t)__rdtsc());
#endif
	maxdepth=calc_maxdepth(image, inflation);
	maxlevels=1<<maxdepth;
	loud=loud_transforms;
	loud_transforms=0;
	im2=0;
	image_copy(&im2, image);
	hist=(int*)malloc(maxlevels*sizeof(int));
	if(!im2||!hist)
	{
		CRASH("Alloc error");
		return;
	}
	memcpy(params2, params, sizeof(params2));
	CALC_LOSS(loss_bestsofar);
	memcpy(loss_prev, loss_bestsofar, sizeof(loss_prev));

	for(it=0, watchdog=0;it<RCT_CUSTOM_NITER;++it)
	{
		int idx[RCT_CUSTOM_DELTAGROUP]={0}, stuck=0;
		int16_t params_original_selected[RCT_CUSTOM_DELTAGROUP]={0};

		if(watchdog>=shakethreshold)/*bump if stuck*/
		{
			memcpy(params, params2, sizeof(params2));
			for(k=0;k<RCT_CUSTOM_NPARAMS-1;++k)
				params[k]+=(rand()&63)-32;
			watchdog=0;
			stuck=1;
		}
		else
		{
			for(k=0;k<RCT_CUSTOM_DELTAGROUP;++k)/*increment params*/
			{
				idx[k]=rand()%RCT_CUSTOM_NPARAMS;
		
				params_original_selected[k]=params[idx[k]];
				if(idx[k]==8)
					params[8]=rand()%6;
				else
					params[idx[k]]+=(rand()&63)-32;
			}
		}

		CALC_LOSS(loss_curr);
		
		if(loss_prev[3]<loss_curr[3])/*revert if worse*/
		{
			if(stuck)/*a bad branch may surpass the local minimum*/
				memcpy(loss_prev, loss_curr, sizeof(loss_prev));
			else
			{
				memcpy(loss_curr, loss_prev, sizeof(loss_prev));
				for(k=0;k<RCT_CUSTOM_DELTAGROUP;++k)
					params[idx[k]]=params_original_selected[k];
			}
			++watchdog;
		}
		else/*save if better*/
		{
			if(loss_curr[3]<loss_bestsofar[3])/*publish if record best*/
			{
				memcpy(params2, params, sizeof(params2));
				memcpy(loss_bestsofar, loss_curr, sizeof(loss_bestsofar));
				--it;/*again*/
			}
			memcpy(loss_prev, loss_curr, sizeof(loss_prev));
			watchdog=0;
		}
		
		/*if(loud_transforms)*/
			set_window_title(
				"%d %4d/%4d,%d/%d: %f%% RGB %f %f %f%s",
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

		/*preview*/
#if 1
		{
			ch_entropy[0]=(float)loss_bestsofar[0];
			ch_entropy[1]=(float)loss_bestsofar[1];
			ch_entropy[2]=(float)loss_bestsofar[2];
			ch_entropy[3]=0;
			io_render();
		}
#endif
	}
	loud_transforms=loud;
	memcpy(params, params2, sizeof(params2));
#undef  CALC_LOSS
	free(hist);
	free(im2);
}
void rct_custom_getmatrix(double *matrix, int fwd)
{
	uint8_t per[4]={0};
	int k;

	rct_custom_unpackpermutation(rct_custom_params[8], per);
	memset(matrix, 0, sizeof(double[9]));
	matrix[0]=1;
	matrix[4]=1;
	matrix[8]=1;
	for(k=0;k<3;++k)
	{
		double vrtx[3];

		vrtx[0]=matrix[k+3*(fwd?per[0]:0)];
		vrtx[1]=matrix[k+3*(fwd?per[1]:1)];
		vrtx[2]=matrix[k+3*(fwd?per[2]:2)];
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
	}
}

#if 0
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
static void c_cbrt(double *dst, const double *x)/*sqrt(x)=exp(0.5*lnx)*/
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
static void impl_solve_cubic(const double *coeffs, double *roots)/*finds r[0], r[1], & r[2], the solutions of x^3 + c[2]x^2 + c[1]x + c[0] = 0*/
{
	/*
	https://math.stackexchange.com/questions/15865/why-not-write-the-solutions-of-a-cubic-this-way/18873#18873
	*/
	double sqrt27, inv3cbrt2, ninth, cm[2], cp[2];
	double p=coeffs[2], q=coeffs[1], r=coeffs[0], p2, p3, q2, pq, A[2], B[2];

	sqrt27=3*sqrt(3);
	inv3cbrt2=1./(3*pow(2, 1./3));
	ninth=1./9;
	cm[0]=-0.5; cm[1]=-0.5*sqrt(3);
	cp[0]=-0.5; cp[1]=+0.5*sqrt(3);
	if(fabs(p)<1e-9&&fabs(q)<1e-9)
	{
		roots[4]=roots[2]=roots[0]=-pow(r, 1./3);
		roots[5]=roots[3]=roots[1]=0;
		return;
	}

	p2=p*p;
	p3=p2*p;
	q2=q*q;
	pq=p*q;

	A[0]=(4*q-p2)*q2+(4*p3-18*pq+27*r)*r;
	if(A[0]<0)/*complex sqrt*/
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
static void impl_rref(double *m, int16_t dx, int16_t dy)
{
#ifdef _DEBUG
	double pivot;
#endif
	double coeff;
	int mindim=dx<dy?dx:dy, it, ky, kx, npivots, kpivot;

	for(it=0, npivots=0;it<mindim;++it)/*iteration*/
	{
		kpivot=-1;
		for(ky=npivots;ky<dy;++ky)/*find pivot*/
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
			for(kx=0;kx<dx;++kx)/*swap rows*/
				coeff=m[dx*kpivot+kx], m[dx*kpivot+kx]=m[dx*(npivots-1)+kx], m[dx*(npivots-1)+kx]=coeff;
			kpivot=npivots-1;
		}
		for(ky=0;ky<dy;++ky)
		{
			if(ky==kpivot)/*normalize pivot row*/
			{
				coeff=1/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*kpivot+kx]*=coeff;
			}
			else/*subtract pivot row from all other rows*/
			{
				coeff=m[dx*ky+it]/m[dx*kpivot+it];
				for(kx=it;kx<dx;++kx)
					m[dx*ky+kx]-=coeff*m[dx*kpivot+kx];
			}
			/*print_matrix_debug(m, dx, dy);*/
		}
	}
}
static int impl_nullspace(double *M, int dx, int dy, double *solution, int sstart, char *dep_flags, int16_t *row_idx)
{
	/*
	M is rref'ed
	solution allocated size dy*dy, pre-memset solution to zero
	p_flags & row_idx both size dx
	returns number of vectors in solution
	*/
	int kx, kxdst, ky, keq, kfree, idx, idx2, nvec;
#ifdef DEBUG_NULLSPACE
	printf("Before RREF:\n");
	print_matrix_debug(M, dx, dy);
#endif
	impl_rref(M, (int16_t)dx, (int16_t)dy);/*first, get rref(M)*/
#ifdef DEBUG_NULLSPACE
	printf("After RREF:\n");
	print_matrix_debug(M, dx, dy);
#endif
	memset(dep_flags, 0, dx);
	memset(row_idx, 0, dx*sizeof(int16_t));
	for(ky=0;ky<dy;++ky)/*find pivots (dependent variables)*/
	{
		for(kx=ky;kx<dx;++kx)/*find first nonzero element*/
		{
			idx=dx*ky+kx;
			if(fabs(M[idx])>1e-10)
				break;
		}
		if(kx<dx)/*if found*/
			dep_flags[kx]=1, row_idx[ky]=(int16_t)kx;
		else
			break;
	}
	nvec=dx-ky;
	for(ky=0, keq=0, kfree=0;ky<dx;++ky)
	{
		if(dep_flags[ky])/*pivot, dependent variable*/
		{
			for(kx=0, kxdst=0;kx<dx;++kx)
			{
				if(dep_flags[kx])
					continue;
				idx=dx*ky+kxdst, idx2=dx*keq+kx;
				if(sstart+idx>=9)
					CRASH("");
				solution[sstart+idx]=-M[idx2];
				++kxdst;
			}
			++keq;
		}
		else/*free variable*/
		{
			idx=dx*ky;
			if(sstart+idx+nvec>9)
				CRASH("");
			memset(solution+sstart+idx, 0, nvec*sizeof(double));
			if(sstart+idx+kfree>=9)
				CRASH("");
			solution[sstart+idx+kfree]=1;
			++kfree;
		}
#ifdef DEBUG_NULLSPACE
		/*printf("Nullspace row %d:\n", ky);*/
		print_matrix_debug(solution+dx*ky, nvec, 1);/**/
#endif
	}
	return nvec;
}
static int impl_egvec(double const *M, int n, const double *lambdas, double *S)
{
	int kv, kx, nvec, size=n*n, again;
	double *temp=(double*)malloc(size*sizeof(double));
	char *dep_flags=(char*)malloc(n);
	int16_t *row_idx=(int16_t*)malloc(n*sizeof(int16_t));
	if(!temp||!dep_flags||!row_idx)
	{
		CRASH("Alloc error");
		return 0;
	}
	memset(S, 0, size*sizeof(double));
	for(kv=0, nvec=0;kv<n&&nvec<n;++kv)/*get eigenvectors*/
	{
		again=0;
		for(kx=0;kx<kv;++kx)/*check for repeated eigenvalues*/
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
	}
	free(dep_flags), free(row_idx);
	free(temp);
	return nvec;
}
static void impl_egval3(double const *M, double *lambdas)/*finds the complex eigenvalues of a 3x3 matrix*/
{
	double C[3];

	C[2]=-(M[0]+M[4]+M[8]);/*-tr(M)*/
	C[1]=M[4]*M[8]-M[5]*M[7] + M[0]*M[8]-M[2]*M[6] + M[0]*M[4]-M[1]*M[3];/*cof0+cof4+cof8*/
	C[0]=-(M[0]*(M[4]*M[8]-M[5]*M[7]) - M[1]*(M[3]*M[8]-M[5]*M[6]) + M[2]*(M[3]*M[7]-M[4]*M[6]));/*-det(M)*/

	impl_solve_cubic(C, lambdas);
}
#endif




/*spatial transforms*/
void filt_deint422(Image *src)
{
	int ky, kx;

	for(ky=0;ky<src->ih;++ky)
	{
		for(kx=1;kx<src->iw-1;kx+=2)
		{
			int
				duW=((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+0)<<2|1],
				duE=((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|1],
				dvW=((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+0)<<2|2],
				dvE=((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|2];
			int
				suW=(+3*duW-1*duE+1)>>1,
				suE=(-1*duW+3*duE+1)>>1;
			int
				svW=(+3*dvW-1*dvE+1)>>1,
				svE=(-1*dvW+3*dvE+1)>>1;
			((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+0)<<2|1]=suW;
			((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|1]=suE;
			((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+0)<<2|2]=svW;
			((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|2]=svE;
		}
	}
}
void filt_deint420(Image *src)
{
	int ky, kx;

	for(ky=1;ky<src->ih-1;ky+=2)
	{
		for(kx=1;kx<src->iw-1;kx+=2)
		{
			int
				duNW=((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+0)<<2|1],
				duNE=((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|1],
				duSW=((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+0)<<2|1],
				duSE=((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+1)<<2|1],
				dvNW=((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+0)<<2|2],
				dvNE=((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|2],
				dvSW=((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+0)<<2|2],
				dvSE=((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+1)<<2|2];
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
			((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+0)<<2|1]=suNW;
			((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|1]=suNE;
			((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+0)<<2|1]=suSW;
			((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+1)<<2|1]=suSE;
			((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+0)<<2|2]=svNW;
			((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|2]=svNE;
			((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+0)<<2|2]=svSW;
			((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+1)<<2|2]=svSE;
		}
	}
}

void packsign(Image *src, int fwd)
{
	intptr_t nvals, k;
	int half[4], kc;

	nvals=(intptr_t)src->iw*src->ih<<2;
	half[0]=1<<src->depth[0]>>1;
	half[1]=1<<src->depth[1]>>1;
	half[2]=1<<src->depth[2]>>1;
	half[3]=1<<src->depth[3]>>1;
	if(fwd)
	{
		for(k=0;k<nvals;k+=4)
		{
			for(kc=0;kc<src->nch;++kc)
			{
				int val=((int32_t*)HDR_DATA(src))[k+kc];
				val=val<<1^-(val<0);
				val-=half[kc];
				((int32_t*)HDR_DATA(src))[k+kc]=val;
			}
		}
	}
	else
	{
		for(k=0;k<nvals;k+=4)
		{
			for(kc=0;kc<src->nch;++kc)
			{
				int val=((int32_t*)HDR_DATA(src))[k+kc];
				val+=half[kc];
				val=val>>1^-(val&1);
				((int32_t*)HDR_DATA(src))[k+kc]=val;
			}
		}
	}
}
void pred_MTF(Image *src, int fwd)
{
	int maxdepth, *dict, idx, ky, kx, kc, ks;

	maxdepth=MAXVAR(src->depth[0], src->depth[1]);
	UPDATE_MAX(maxdepth, src->depth[2]);
	UPDATE_MAX(maxdepth, src->depth[3]);
	dict=(int*)malloc(sizeof(int[4])<<maxdepth);
	if(!dict)
	{
		CRASH("Alloc error");
		return;
	}
	for(kc=0;kc<4;++kc)
	{
		int depth, *curr_dict, nlevels;

		depth=src->depth[kc];
		if(!depth)
			continue;
		curr_dict=dict+((size_t)kc<<maxdepth);
		nlevels=1<<src->depth[kc];
		for(ks=0;ks<nlevels;++ks)
			curr_dict[ks]=ks;
	}
	if(fwd)
	{
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx, idx+=4)
			{
				for(kc=0;kc<4;++kc)
				{
					int depth, nlevels, half, *curr_dict, curr;

					depth=src->depth[kc];
					if(!depth)
						continue;
					nlevels=1<<depth;
					half=nlevels>>1;
					curr_dict=dict+((size_t)kc<<maxdepth);
					curr=((int32_t*)HDR_DATA(src))[idx|kc]+half;
					for(ks=0;ks<nlevels&&curr_dict[ks]!=curr;++ks);
					((int32_t*)HDR_DATA(src))[idx|kc]=ks-half;
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
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx, idx+=4)
			{
				for(kc=0;kc<4;++kc)
				{
					int depth, nlevels, half, *curr_dict, curr;

					depth=src->depth[kc];
					if(!depth)
						continue;
					nlevels=1<<depth;
					half=nlevels>>1;
					curr_dict=dict+((size_t)kc<<maxdepth);
					ks=((int32_t*)HDR_DATA(src))[idx|kc]+half;
					curr=curr_dict[ks];
					((int32_t*)HDR_DATA(src))[idx|kc]=curr-half;
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
void pred_signbit(Image *image, int fwd)
{
	int k, kc, res;

	if(fwd)
	{
		for(kc=0;kc<4;++kc)
		{
			if(!image->depth[kc])
				continue;
			for(k=0, res=image->iw*image->ih;k<res;++k)
			{
				int val=((int32_t*)HDR_DATA(image))[k<<2|kc];
				((int32_t*)HDR_DATA(image))[k<<2|kc]=val<0?-128:val>0?127:0;
			}
		}
	}
	else
	{
		for(kc=0;kc<4;++kc)
		{
			if(!image->depth[kc])
				continue;
			for(k=0, res=image->iw*image->ih;k<res;++k)
			{
				int val=((int32_t*)HDR_DATA(image))[k<<2|kc];
				((int32_t*)HDR_DATA(image))[k<<2|kc]=abs(val)-128;
			}
		}
	}
}

/*clamped gradient / LOCO-I / Median Edge Detector (MED) predictor from JPEG-LS*/
void pred_clampgrad(Image *src, int fwd, int enable_ma)
{
	int rstride, sh0, sh1, sh2, ky, kx;

	/*0.123 0.121 0.117 sec		0.157 sec*/
	rstride=4*src->iw;
	sh0=32-src->depth[0];
	sh1=32-src->depth[1];
	sh2=32-src->depth[2];
	if(g_dist>1)
		goto cg_near;
	if(fwd)
	{
		int *ptr=((int32_t*)HDR_DATA(src))+4*src->iw*src->ih-4;
		for(ky=src->ih-1;ky>=1;--ky)
		{
			for(kx=src->iw-1;kx>=1;--kx)
			{
				int
					NW0	=ptr[-1*rstride-1*4+0],
					NW1	=ptr[-1*rstride-1*4+1],
					NW2	=ptr[-1*rstride-1*4+2],
					N0	=ptr[-1*rstride+0*4+0],
					N1	=ptr[-1*rstride+0*4+1],
					N2	=ptr[-1*rstride+0*4+2],
					W0	=ptr[+0*rstride-1*4+0],
					W1	=ptr[+0*rstride-1*4+1],
					W2	=ptr[+0*rstride-1*4+2],
					c0	=ptr[0],
					c1	=ptr[1],
					c2	=ptr[2];
				int p0=N0+W0-NW0;
				int p1=N1+W1-NW1;
				int p2=N2+W2-NW2;
				int ymax=N0, ymin=W0;
				int umax=N1, umin=W1;
				int vmax=N2, vmin=W2;
				if(N0<W0)ymin=N0, ymax=W0;
				if(N1<W1)umin=N1, umax=W1;
				if(N2<W2)vmin=N2, vmax=W2;
				CLAMP2(p0, ymin, ymax);
				CLAMP2(p1, umin, umax);
				CLAMP2(p2, vmin, vmax);
				c0-=p0;
				c1-=p1;
				c2-=p2;
				c0<<=sh0;
				c1<<=sh1;
				c2<<=sh2;
				c0>>=sh0;
				c1>>=sh1;
				c2>>=sh2;
				ptr[0]=c0;
				ptr[1]=c1;
				ptr[2]=c2;
				ptr-=4;
			}
			{
				int
					N0	=ptr[-1*rstride+0*4+0],
					N1	=ptr[-1*rstride+0*4+1],
					N2	=ptr[-1*rstride+0*4+2],
					c0	=ptr[0],
					c1	=ptr[1],
					c2	=ptr[2];
				c0-=N0;
				c1-=N1;
				c2-=N2;
				c0<<=sh0;
				c1<<=sh1;
				c2<<=sh2;
				c0>>=sh0;
				c1>>=sh1;
				c2>>=sh2;
				ptr[0]=c0;
				ptr[1]=c1;
				ptr[2]=c2;
				ptr-=4;
			}
		}
		for(kx=src->iw-1;kx>=1;--kx)
		{
			int
				W0	=ptr[+0*rstride-1*4+0],
				W1	=ptr[+0*rstride-1*4+1],
				W2	=ptr[+0*rstride-1*4+2],
				c0	=ptr[0],
				c1	=ptr[1],
				c2	=ptr[2];
			c0-=W0;
			c1-=W1;
			c2-=W2;
			c0<<=sh0;
			c1<<=sh1;
			c2<<=sh2;
			c0>>=sh0;
			c1>>=sh1;
			c2>>=sh2;
			ptr[0]=c0;
			ptr[1]=c1;
			ptr[2]=c2;
			ptr-=4;
		}
	}
	else
	{
		int *ptr=((int32_t*)HDR_DATA(src))+4;
		for(kx=1;kx<src->iw;++kx)
		{
			int
				W0	=ptr[+0*rstride-1*4+0],
				W1	=ptr[+0*rstride-1*4+1],
				W2	=ptr[+0*rstride-1*4+2],
				c0	=ptr[0],
				c1	=ptr[1],
				c2	=ptr[2];
			c0+=W0;
			c1+=W1;
			c2+=W2;
			c0<<=sh0;
			c1<<=sh1;
			c2<<=sh2;
			c0>>=sh0;
			c1>>=sh1;
			c2>>=sh2;
			ptr[0]=c0;
			ptr[1]=c1;
			ptr[2]=c2;
			ptr+=4;
		}
		for(ky=1;ky<src->ih;++ky)
		{
			{
				int
					N0	=ptr[-1*rstride+0*4+0],
					N1	=ptr[-1*rstride+0*4+1],
					N2	=ptr[-1*rstride+0*4+2],
					c0	=ptr[0],
					c1	=ptr[1],
					c2	=ptr[2];
				c0+=N0;
				c1+=N1;
				c2+=N2;
				c0<<=sh0;
				c1<<=sh1;
				c2<<=sh2;
				c0>>=sh0;
				c1>>=sh1;
				c2>>=sh2;
				ptr[0]=c0;
				ptr[1]=c1;
				ptr[2]=c2;
				ptr+=4;
			}
			for(kx=1;kx<src->iw;++kx)
			{
				int
					NW0	=ptr[-1*rstride-1*4+0],
					NW1	=ptr[-1*rstride-1*4+1],
					NW2	=ptr[-1*rstride-1*4+2],
					N0	=ptr[-1*rstride+0*4+0],
					N1	=ptr[-1*rstride+0*4+1],
					N2	=ptr[-1*rstride+0*4+2],
					W0	=ptr[+0*rstride-1*4+0],
					W1	=ptr[+0*rstride-1*4+1],
					W2	=ptr[+0*rstride-1*4+2],
					c0	=ptr[0],
					c1	=ptr[1],
					c2	=ptr[2];
				int p0=N0+W0-NW0;
				int p1=N1+W1-NW1;
				int p2=N2+W2-NW2;
				int ymax=N0, ymin=W0;
				int umax=N1, umin=W1;
				int vmax=N2, vmin=W2;
				if(N0<W0)ymin=N0, ymax=W0;
				if(N1<W1)umin=N1, umax=W1;
				if(N2<W2)vmin=N2, vmax=W2;
				CLAMP2(p0, ymin, ymax);
				CLAMP2(p1, umin, umax);
				CLAMP2(p2, vmin, vmax);
				c0+=p0;
				c1+=p1;
				c2+=p2;
				c0<<=sh0;
				c1<<=sh1;
				c2<<=sh2;
				c0>>=sh0;
				c1>>=sh1;
				c2>>=sh2;
				ptr[0]=c0;
				ptr[1]=c1;
				ptr[2]=c2;
				ptr+=4;
			}
		}
	}
	return;

	/*0.218 0.220 0.214 sec*/
cg_near:
	{
		int nlevels[4], amin[4], amax[4], invdist, psize, kc, idx;
		int16_t *pixels;

		nlevels[0]=1<<src->depth[0];
		nlevels[1]=1<<src->depth[1];
		nlevels[2]=1<<src->depth[2];
		nlevels[3]=1<<src->depth[3];
		amin[0]=-(nlevels[0]>>1);
		amin[1]=-(nlevels[1]>>1);
		amin[2]=-(nlevels[2]>>1);
		amin[3]=-(nlevels[3]>>1);
		amax[0]=(nlevels[0]>>1)-1;
		amax[1]=(nlevels[1]>>1)-1;
		amax[2]=(nlevels[2]>>1)-1;
		amax[3]=(nlevels[3]>>1)-1;
		invdist=((1<<16)+g_dist-1)/g_dist;
		psize=(src->iw+2)*(int)sizeof(int16_t[2*4]);
		pixels=(int16_t*)malloc(psize);/*2 padded rows * 4 channels max*/
		if(!pixels)
		{
			CRASH("Alloc error");
			return;
		}
		memset(pixels, 0, psize);
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			int16_t *rows[2];

			rows[0]=pixels+(((src->iw+(intptr_t)2)*((ky-(intptr_t)0)&1)+1)<<2);
			rows[1]=pixels+(((src->iw+(intptr_t)2)*((ky-(intptr_t)1)&1)+1)<<2);
			for(kx=0;kx<src->iw;++kx, idx+=4)
			{
				for(kc=0;kc<src->nch;++kc)
				{
					int
						NW	=rows[1][kc-4],
						N	=rows[1][kc+0],
						W	=rows[0][kc-4];
					int pred, vmin, vmax, curr, e;

					pred=N+W-NW;
					vmax=N; vmin=W;
					if(N<W)vmin=N, vmax=W;
					CLAMP2(pred, vmin, vmax);

					curr=((int32_t*)HDR_DATA(src))[idx+kc];
					if(g_dist>1)
					{
						if(fwd)
						{
							curr-=(int)pred;
							curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));/*curr/=g_dist*/
							((int32_t*)HDR_DATA(src))[idx+kc]=curr;

							curr=g_dist*curr+(int)pred;
							CLAMP2(curr, amin[kc], amax[kc]);
						}
						else
						{
							curr=g_dist*curr+(int)pred;
							CLAMP2(curr, amin[kc], amax[kc]);

							((int32_t*)HDR_DATA(src))[idx+kc]=curr;
						}
					}
					else
					{
						if(fwd)
						{
							e=curr-pred;
							e<<=32-src->depth[kc];
							e>>=32-src->depth[kc];
							((int32_t*)HDR_DATA(src))[idx+kc]=e;
						}
						else
						{
							curr+=pred;
							curr<<=32-src->depth[kc];
							curr>>=32-src->depth[kc];
							((int32_t*)HDR_DATA(src))[idx+kc]=curr;
						}
					}
					rows[0][kc]=curr;
				}

				rows[0]+=4;
				rows[1]+=4;
			}
		}
		free(pixels);
	}
	(void)enable_ma;
}
void pred_cg_crct(Image *src, int fwd, int enable_ma)
{
	int nlevels[4], amin[4], amax[4], invdist, psize, yidx, uidx, vidx, uc0, vc0, vc1, ky, kx, kc, idx;
	int32_t *pixels;
	const uint8_t *combination;

	if(fwd)
		src->rct=crct_analysis(src);
	nlevels[0]=1<<src->depth[0];
	nlevels[1]=1<<src->depth[1];
	nlevels[2]=1<<src->depth[2];
	nlevels[3]=1<<src->depth[3];
	amin[0]=-(nlevels[0]>>1);
	amin[1]=-(nlevels[1]>>1);
	amin[2]=-(nlevels[2]>>1);
	amin[3]=-(nlevels[3]>>1);
	amax[0]=(nlevels[0]>>1)-1;
	amax[1]=(nlevels[1]>>1)-1;
	amax[2]=(nlevels[2]>>1)-1;
	amax[3]=(nlevels[3]>>1)-1;
	invdist=((1<<16)+g_dist-1)/g_dist;
	psize=(src->iw+2)*(int)sizeof(int32_t[2*4]);
	pixels=(int32_t*)malloc(psize);/*2 padded rows * 4 channels max*/
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	combination=rct_combinations[src->rct];
	yidx=combination[II_PERM_Y];
	uidx=combination[II_PERM_U];
	vidx=combination[II_PERM_V];
	uc0=combination[II_COEFF_U_SUB_Y];
	vc0=combination[II_COEFF_V_SUB_Y];
	vc1=combination[II_COEFF_V_SUB_U];
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int32_t *rows[2];

		rows[0]=pixels+((src->iw+(intptr_t)2)*((ky-(intptr_t)0)&1)+1)*4-1;
		rows[1]=pixels+((src->iw+(intptr_t)2)*((ky-(intptr_t)1)&1)+1)*4-1;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			int offset=0, yuv[3];

			yuv[0]=((int32_t*)HDR_DATA(src))[idx+yidx];
			yuv[1]=((int32_t*)HDR_DATA(src))[idx+uidx];
			yuv[2]=((int32_t*)HDR_DATA(src))[idx+vidx];
			for(kc=0;kc<4;++kc)
			{
				int NW, N, W, pred, vmin, vmax, curr, e;

				++rows[0];
				++rows[1];
				if(!src->depth[kc])
					continue;
				NW	=rows[1][-1*4];
				N	=rows[1][+0*4];
				W	=rows[0][-1*4];
				pred=N+W-NW;
				vmax=N; vmin=W;
				if(N<W)vmin=N, vmax=W;
				CLAMP2(pred, vmin, vmax);
				if(kc)
				{
					pred+=offset;
					CLAMP2(pred, amin[kc], amax[kc]);
				}
				curr=yuv[kc];
				if(g_dist>1)
				{
					if(fwd)
					{
						curr-=pred;
						curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));/*curr/=g_dist*/
						((int32_t*)HDR_DATA(src))[idx+kc]=curr;

						curr=g_dist*curr+pred;
					}
					else
						curr=g_dist*((int32_t*)HDR_DATA(src))[idx+kc]+pred;
					CLAMP2(curr, amin[kc], amax[kc]);
					yuv[kc]=curr;
				}
				else
				{
					if(fwd)
					{
						e=curr-pred;
						e<<=32-src->depth[kc];
						e>>=32-src->depth[kc];
						((int32_t*)HDR_DATA(src))[idx+kc]=e;
					}
					else
					{
						curr=((int32_t*)HDR_DATA(src))[idx+kc]+pred;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						yuv[kc]=curr;
					}
				}
				curr-=offset;
				rows[0][0]=curr;
				offset=kc?(vc0*yuv[0]+vc1*yuv[1])>>2:uc0*yuv[0]>>2;
			}
			if(!fwd)
			{
				((int32_t*)HDR_DATA(src))[idx+yidx]=yuv[0];
				((int32_t*)HDR_DATA(src))[idx+uidx]=yuv[1];
				((int32_t*)HDR_DATA(src))[idx+vidx]=yuv[2];
			}
		}
	}
	free(pixels);

	(void)enable_ma;
}
void pred_cleartype(Image *src, int fwd)
{
	int ky, kx, idx, psize;
	int16_t *pixels;

	psize=(src->iw+2)*(int)sizeof(int16_t[2*4]);
	pixels=(int16_t*)malloc(psize);/*2 padded rows * 4 channels max*/
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int16_t *rows[2];

		rows[0]=pixels+((src->iw+(intptr_t)2)*((ky-(intptr_t)0)&1)+1)*4;
		rows[1]=pixels+((src->iw+(intptr_t)2)*((ky-(intptr_t)1)&1)+1)*4;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			int
				rNW	=rows[1][0-1*4],
				gNW	=rows[1][1-1*4],
				bNW	=rows[1][2-1*4],
				rN	=rows[1][0+0*4],
				gN	=rows[1][1+0*4],
				bN	=rows[1][2+0*4],
				rW	=rows[0][0-1*4],
				gW	=rows[0][1-1*4],
				bW	=rows[0][2-1*4],
				r	=((int32_t*)HDR_DATA(src))[idx+0],
				g	=((int32_t*)HDR_DATA(src))[idx+1],
				b	=((int32_t*)HDR_DATA(src))[idx+2];
			int pred;

			(void)rNW;
			(void)gNW;
			(void)bNW;
			(void)rN;
			(void)gN;
			(void)bN;
			(void)rW;
			(void)gW;
			(void)bW;
			(void)r;
			(void)g;
			(void)b;
			/*
			rNW gNW bNW rN gN bN
			rW  gW  bW  r  g  b
			pred(r) = SELECT(rN, bW, bNW)
			pred(g) = SELECT(gN, r, rN)
			pred(b) = SELECT(bN, g, gN)
			*/
#define PRED_SELECT(aN, aW, aNW) (abs((aN)-(aNW))<abs((aW)-(aNW))?aW:aN)

			pred=PRED_SELECT(rN, bW, bNW);
			if(fwd)
				pred=-pred;
			pred=(pred+r)<<(32-src->depth[0])>>(32-src->depth[0]);
			((int32_t*)HDR_DATA(src))[idx+0]=pred;
			if(!fwd)
				r=pred;

			pred=PRED_SELECT(gN, r, rN);
			if(fwd)
				pred=-pred;
			pred=(pred+g)<<(32-src->depth[1])>>(32-src->depth[1]);
			((int32_t*)HDR_DATA(src))[idx+1]=pred;
			if(!fwd)
				g=pred;
			pred=PRED_SELECT(bN, g, gN);
			if(fwd)
				pred=-pred;
			pred=(pred+b)<<(32-src->depth[2])>>(32-src->depth[2]);
			((int32_t*)HDR_DATA(src))[idx+2]=pred;
			if(!fwd)
				b=pred;
			rows[0][0]=r;
			rows[0][1]=g;
			rows[0][2]=b;
#undef  PRED_SELECT

			rows[0]+=4;
			rows[1]+=4;
		}
	}
	free(pixels);
}
void pred_sse(Image *src, int fwd)
{
	enum
	{
		cellsize=9
	};
	int nlevels[4], mask[4], ssesize, *sse, bufsize, ky, kx, kc, idx;
	int16_t *pixels;

	nlevels[0]=1<<src->depth[0];
	nlevels[1]=1<<src->depth[1];
	nlevels[2]=1<<src->depth[2];
	nlevels[3]=1<<src->depth[3];
	mask[0]=nlevels[0]-1;
	mask[1]=nlevels[1]-1;
	mask[2]=nlevels[2]-1;
	mask[3]=nlevels[3]-1;
	ssesize=(int)sizeof(int)*cellsize*(nlevels[0]+nlevels[1]+nlevels[2]+nlevels[3]);
	sse=(int*)malloc(ssesize);
	bufsize=(src->iw+16)*(int)sizeof(int16_t[4*4*3]);/*4 padded rows * 4 channels max * {pixels, errors}*/
	pixels=(int16_t*)malloc(bufsize);
	if(!sse||!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	memset(sse, 0, ssesize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int16_t *rows[4];

		rows[0]=pixels+((src->iw+(intptr_t)8)*((ky-(intptr_t)0)&3)+1)*4*3;
		rows[1]=pixels+((src->iw+(intptr_t)8)*((ky-(intptr_t)1)&3)+1)*4*3;
		rows[2]=pixels+((src->iw+(intptr_t)8)*((ky-(intptr_t)2)&3)+1)*4*3;
		rows[3]=pixels+((src->iw+(intptr_t)8)*((ky-(intptr_t)3)&3)+1)*4*3;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			int *currsse=sse;
			for(kc=0;kc<src->nch;++kc)
			{
				int
					bNW	=rows[1][kc-1*4*3+4],
					bN	=rows[1][kc+0*4*3+4],
					bNE	=rows[1][kc+1*4*3+4],
					bNEE	=rows[1][kc+2*4*3+4],
					bW	=rows[0][kc-1*4*3+4],
					yNW	=rows[1][kc-1*4*3+8],
					yN	=rows[1][kc+0*4*3+8],
					yNE	=rows[1][kc+1*4*3+8],
					yNEE	=rows[1][kc+2*4*3+8],
					yW	=rows[0][kc-1*4*3+8],
					NNNW	=rows[3][kc-1*4*3],
					NNN	=rows[3][kc+0*4*3],
					NNNE	=rows[3][kc+1*4*3],
					NNWW	=rows[2][kc-2*4*3],
					NNW	=rows[2][kc-1*4*3],
					NN	=rows[2][kc+0*4*3],
					NNE	=rows[2][kc+1*4*3],
					NNEE	=rows[2][kc+2*4*3],
					NWW	=rows[1][kc-2*4*3],
					NW	=rows[1][kc-1*4*3],
					N	=rows[1][kc+0*4*3],
					NE	=rows[1][kc+1*4*3],
					NEE	=rows[1][kc+2*4*3],
					NEEE	=rows[1][kc+3*4*3],
					NEEEE	=rows[1][kc+4*4*3],
					NEEEEE	=rows[1][kc+5*4*3],
					WWWWW	=rows[0][kc-5*4*3],
					WWWW	=rows[0][kc-4*4*3],
					WWW	=rows[0][kc-3*4*3],
					WW	=rows[0][kc-2*4*3],
					W	=rows[0][kc-1*4*3];
				int64_t psum;
				int preds[4], bias, vmax, vmin, pred, target, e, sseidx, *cells, wsum, weight;

#if 1
				(void)bNW;
				(void)bN;
				(void)bNE;
				(void)bNEE;
				(void)bW;
				(void)yNW;
				(void)yN;
				(void)yNE;
				(void)yNEE;
				(void)yW;
				(void)NNNW;
				(void)NNN;
				(void)NNNE;
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
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
				(void)bias;
#endif
				bias=(yNEE+bW)>>(6+1);
				pred=0;
				target=((int32_t*)HDR_DATA(src))[idx+kc];
				e=target;
				preds[0]=N;
				preds[1]=W;
				preds[2]=N+W-NW;
				preds[3]=(NN+WW+NEE+NEEE)>>2;
				sseidx=(int)(((N-W)>>1&mask[kc])*cellsize);
				cells=currsse+sseidx;
				psum=0;
				wsum=0;
				preds[0]+=cells[4+0]>>8;
				preds[1]+=cells[4+1]>>8;
				preds[2]+=cells[4+2]>>8;
				preds[3]+=cells[4+3]>>8;
				weight=0x1000000/(cells[0]+1); psum+=(int64_t)weight*preds[0]; wsum+=weight;
				weight=0x1000000/(cells[1]+1); psum+=(int64_t)weight*preds[1]; wsum+=weight;
				weight=0x1000000/(cells[2]+1); psum+=(int64_t)weight*preds[2]; wsum+=weight;
				weight=0x1000000/(cells[3]+1); psum+=(int64_t)weight*preds[3]; wsum+=weight;
				pred=(int)(psum/wsum);
				pred+=cells[8]>>8;
				vmax=N; vmin=W;
				if(N<W)vmin=N, vmax=W;
				if(vmin>NE)vmin=NE;
				if(vmax<NE)vmax=NE;
				CLAMP2(pred, vmin, vmax);
				currsse+=nlevels[kc]*cellsize;
				if(fwd)
				{
					e-=pred;
					e<<=32-src->depth[kc];
					e>>=32-src->depth[kc];
					((int32_t*)HDR_DATA(src))[idx+kc]=e;
				}
				else
				{
					target+=pred;
					target<<=32-src->depth[kc];
					target>>=32-src->depth[kc];
					((int32_t*)HDR_DATA(src))[idx+kc]=target;
				}
				cells[0]+=((abs(target-preds[0])<<8)-cells[0]+(1<<3>>1))>>3;
				cells[1]+=((abs(target-preds[1])<<8)-cells[1]+(1<<3>>1))>>3;
				cells[2]+=((abs(target-preds[2])<<8)-cells[2]+(1<<3>>1))>>3;
				cells[3]+=((abs(target-preds[3])<<8)-cells[3]+(1<<3>>1))>>3;
				cells[4+0]+=(((target-preds[0])<<8)-cells[4+0]+(1<<7>>1))>>7;
				cells[4+1]+=(((target-preds[1])<<8)-cells[4+1]+(1<<7>>1))>>7;
				cells[4+2]+=(((target-preds[2])<<8)-cells[4+2]+(1<<7>>1))>>7;
				cells[4+3]+=(((target-preds[3])<<8)-cells[4+3]+(1<<7>>1))>>7;
				cells[8]+=(((target-pred)<<8)-cells[8]+(1<<9>>1))>>9;
				rows[0][kc]=target;
			}
			rows[0]+=4*3;
			rows[1]+=4*3;
			rows[2]+=4*3;
			rows[3]+=4*3;
		}
	}
	free(pixels);
	free(sse);
}
void pred_quantize(Image *src, int fwd)
{
	enum
	{
		XPAD=8,
		NROWS=4,
		NCH=4,
		NVAL=1
	};
	int amin[4], amax[4], ky, kx, kc, idx, psize;
	int16_t *pixels, *rows[4];

	(void)amax;
	(void)amin;
	psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i));/*4 padded rows * 4 channels max * {pixels, drift}*/
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;/*sub 1 channel for pre-increment*/
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				{
					int
						NN	=rows[2][0-1*NCH*NROWS*NVAL],
						NW	=rows[1][0-1*NCH*NROWS*NVAL],
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						NE	=rows[1][0+1*NCH*NROWS*NVAL],
						WW	=rows[0][0-2*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL];
					int val;

					(void)NN;
					(void)NW;
					(void)N;
					(void)NE;
					(void)WW;
					(void)W;
					val=((int32_t*)HDR_DATA(src))[idx];
					if(fwd)
					{
						int smask=val>>31;

						val^=smask;
						val-=smask;
						val>>=1;
						val^=smask;
						val-=smask;
					}
					else
					{
						val<<=1;
#if 0
						val=val<<1|(val-amin[kc])>>(src->depth[kc]-1);
						CLAMP2(val, amin[kc], amax[kc]);
#endif
					}
					rows[0][0]=val;
					((int32_t*)HDR_DATA(src))[idx]=val;
				}
			}
		}
	}
	_mm_free(pixels);
}
static void pred_LeGall_1d_fwd(int *data, int *temp, int stride, int count)
{
	int odd, end, *ptr, ks, kd, kx, *loptr, *hiptr;

	for(ks=0, kd=0;kd<count;ks+=stride, ++kd)
		temp[kd]=data[ks];

	odd=count&1, end=count+odd-3;
	ptr=temp;
	ptr[0]-=ptr[1];
	++ptr;
	for(kx=1;kx<end;kx+=2)
	{
		ptr[ 1]-=(ptr[ 0]+ptr[ 2])>>1;
		ptr[ 0]+=(ptr[-1]+ptr[ 1])>>2;
		ptr+=2;
	}
	if(odd)
	{
		ptr[ 1]-=ptr[0];
		ptr[ 0]+=(ptr[-1]+ptr[ 1])>>2;
	}
	else
	{
		ptr[ 1]-=(ptr[ 0]+ptr[ 2])>>1;
		ptr[ 0]+=(ptr[-1]+ptr[ 1])>>2;
		ptr[ 2]+=ptr[1]>>1;
	}
	loptr=data;
	hiptr=loptr+(size_t)(count>>1)*stride;
	ptr=temp;
	for(kx=0;kx<end;kx+=2)
	{
		*loptr=ptr[1];
		*hiptr=ptr[0];
		ptr+=2;
		loptr+=stride;
		hiptr+=stride;
	}
	*hiptr=ptr[0];
	if(!odd)
		*loptr=ptr[1];
}
static void pred_LeGall_1d_inv(int *data, int *temp, int stride, int count)
{
	int odd, end, *ptr, *loptr, *hiptr, kx, ks, kd;

	odd=count&1, end=count+odd-3;
	ptr=temp;
	loptr=data;
	hiptr=loptr+(count>>1)*stride;
	for(kx=0;kx<end;kx+=2)
	{
		ptr[1]=*loptr;
		ptr[0]=*hiptr;
		ptr+=2;
		loptr+=stride;
		hiptr+=stride;
	}
	ptr[0]=*hiptr;
	if(!odd)
		ptr[1]=*loptr;

	ptr=temp;
	ptr[ 1]-=(ptr[ 0]+ptr[ 2])>>2;
	ptr[ 0]+=ptr[1];
	ptr+=2;
	for(kx=2;kx<end;kx+=2)
	{
		ptr[ 1]-=(ptr[ 0]+ptr[ 2])>>2;
		ptr[ 0]+=(ptr[-1]+ptr[ 1])>>1;
		ptr+=2;
	}
	if(odd)
		ptr[ 0]+=ptr[-1];
	else
	{
		ptr[ 1]-=ptr[ 0]>>1;
		ptr[ 0]+=(ptr[-1]+ptr[ 1])>>1;
	}

	for(ks=0, kd=0;kd<count;ks+=stride, ++kd)
		data[ks]=temp[kd];
}
void pred_LeGallCG(Image *src, int fwd)
{
	int nstages, psize, *pixels, k, kx, ky, kc;

	nstages=31^LZCNT32(src->iw);
	psize=MAXVAR(src->iw+16, src->ih)*(int)sizeof(int[2*4]);/*2 padded rows * 4 channels max*/
	pixels=(int*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	{
		int temp=31^LZCNT32(src->ih);
		if(nstages>temp)
			nstages=temp;
		nstages-=2;
	}
	if(nstages<0)
		return;
	if(fwd)
	{
		for(kc=0;kc<src->nch;++kc)
			src->depth[kc]+=src->depth[kc]<24;
		for(k=0;k<nstages;++k)
		{
			int w2=src->iw>>k, h2=src->ih>>k;
			for(ky=0;ky<h2;++ky)
				for(kc=0;kc<src->nch;++kc)
					pred_LeGall_1d_fwd(((int32_t*)HDR_DATA(src))+((size_t)src->iw*ky<<2|kc), pixels, 4, w2);
			for(kx=0;kx<w2;++kx)
				for(kc=0;kc<src->nch;++kc)
					pred_LeGall_1d_fwd(((int32_t*)HDR_DATA(src))+((size_t)kx<<2|kc), pixels, src->iw<<2, h2);
		}
	}
	else
	{
		for(k=nstages-1;k>=0;--k)
		{
			int w2=src->iw>>k, h2=src->ih>>k;
			for(kx=0;kx<w2;++kx)
				for(kc=0;kc<src->nch;++kc)
					pred_LeGall_1d_inv(((int32_t*)HDR_DATA(src))+((size_t)kx<<2|kc), pixels, src->iw<<2, h2);
			for(ky=0;ky<h2;++ky)
				for(kc=0;kc<src->nch;++kc)
					pred_LeGall_1d_inv(((int32_t*)HDR_DATA(src))+((size_t)src->iw*ky<<2|kc), pixels, 4, w2);
		}
		for(kc=0;kc<src->nch;++kc)
			src->depth[kc]-=src->depth[kc]>src->src_depth[kc];
	}
	free(pixels);
}
void pred_cgplus(Image *src, int fwd)
{
	int kx, ky, kc;

	for(kc=0;kc<src->nch;++kc)
	{
		if(fwd)
		{
			for(ky=0;ky<src->ih;ky+=2)
			{
				for(kx=1;kx<src->iw;kx+=2)	/*AV8*/
				{
					int
						N	=ky-1>=0			?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx+0)<<2|kc]:0,
						W	=kx-1>=0			?((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx-1)<<2|kc]:0,
						S	=ky+1<src->ih			?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+0)<<2|kc]:0,
						E	=kx+1<src->iw			?((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|kc]:0,
						NW	=ky-1>=0&&kx-1>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx-1)<<2|kc]:0,
						NE	=ky-1>=0&&kx+1<src->iw		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx+1)<<2|kc]:0,
						SW	=ky+1<src->ih&&kx-1>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx-1)<<2|kc]:0,
						SE	=ky+1<src->ih&&kx+1<src->iw	?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+1)<<2|kc]:0;
					int pred, *curr;

					pred=(3*(N+W+S+E)-(NW+NE+SW+SE)+4)>>3;

					curr=((int32_t*)HDR_DATA(src))+((src->iw*ky+kx)<<2|kc);
					*curr=(*curr-pred)<<(32-src->depth[kc])>>(32-src->depth[kc]);
				}
			}
			for(ky=1;ky<src->ih;ky+=2)
			{
				for(kx=0;kx<src->iw;kx+=2)	/*AV12A*/
				{
					int
						N	=ky-1>=0			?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx+0)<<2|kc]:0,
						W	=kx-1>=0			?((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx-1)<<2|kc]:0,
						S	=ky+1<src->ih			?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+0)<<2|kc]:0,
						E	=kx+1<src->iw			?((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|kc]:0,
						NNW	=ky-2>=0&&kx-1>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-2)+kx-1)<<2|kc]:0,
						NNE	=ky-2>=0&&kx+1<src->iw		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-2)+kx+1)<<2|kc]:0,
						NWW	=ky-1>=0&&kx-2>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx-2)<<2|kc]:0,
						NEE	=ky-1>=0&&kx+2<src->iw		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx+2)<<2|kc]:0,
						SWW	=ky+1<src->ih&&kx-2>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx-2)<<2|kc]:0,
						SEE	=ky+1<src->ih&&kx+2<src->iw	?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+2)<<2|kc]:0,
						SSW	=ky+2<src->ih&&kx-1>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky+2)+kx-1)<<2|kc]:0,
						SSE	=ky+2<src->ih&&kx+1<src->iw	?((int32_t*)HDR_DATA(src))[(src->iw*(ky+2)+kx+1)<<2|kc]:0;
					int vmax, vmin, pred, *curr;

					vmax=N; vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>S)vmin=S;
					if(vmax<S)vmax=S;
					if(vmin>E)vmin=E;
					if(vmax<E)vmax=E;
					pred=(6*(N+W+S+E)-(NNW+NNE+NWW+NEE+SWW+SEE+SSW+SSE)+8)>>4;
					CLAMP2(pred, vmin, vmax);

					curr=((int32_t*)HDR_DATA(src))+((src->iw*ky+kx)<<2|kc);
					*curr=(*curr-pred)<<(32-src->depth[kc])>>(32-src->depth[kc]);
				}
			}
			for(ky=1;ky<src->ih;ky+=2)
			{
				for(kx=1;kx<src->iw;kx+=2)	/*AV4B*/
				{
					int
						NW	=ky-1>=0&&kx-1>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx-1)<<2|kc]:0,
						NE	=ky-1>=0&&kx+1<src->iw		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx+1)<<2|kc]:0,
						SW	=ky+1<src->ih&&kx-1>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx-1)<<2|kc]:0,
						SE	=ky+1<src->ih&&kx+1<src->iw	?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+1)<<2|kc]:0;
					int pred, *curr;

					pred=(NW+NE+SW+SE+2)>>2;

					curr=((int32_t*)HDR_DATA(src))+((src->iw*ky+kx)<<2|kc);
					*curr=(*curr-pred)<<(32-src->depth[kc])>>(32-src->depth[kc]);
				}
			}
			for(ky=(src->ih-1)&~1;ky>=0;ky-=2)
			{
				for(kx=(src->iw-1)&~1;kx>=0;kx-=2)	/*CG2*/
				{
					int
						NNWW	=ky-2>=0&&kx-2>=0	?((int32_t*)HDR_DATA(src))[(src->iw*(ky-2)+kx-2)<<2|kc]:0,
						NN	=ky-2>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-2)+kx+0)<<2|kc]:0,
						WW	=kx-2>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx-2)<<2|kc]:0;
					int pred, vmin, vmax, *curr;

					vmax=NN; vmin=WW;
					if(NN<WW)vmin=NN, vmax=WW;
					pred=NN+WW-NNWW;
					CLAMP2(pred, vmin, vmax);

					curr=((int32_t*)HDR_DATA(src))+((src->iw*ky+kx)<<2|kc);
					*curr=(*curr-pred)<<(32-src->depth[kc])>>(32-src->depth[kc]);
				}
			}
		}
		else/*inv*/
		{
			for(ky=0;ky<src->ih;ky+=2)
			{
				for(kx=0;kx<src->iw;kx+=2)	/*CG2*/
				{
					int
						NNWW	=ky-2>=0&&kx-2>=0	?((int32_t*)HDR_DATA(src))[(src->iw*(ky-2)+kx-2)<<2|kc]:0,
						NN	=ky-2>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-2)+kx+0)<<2|kc]:0,
						WW	=kx-2>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx-2)<<2|kc]:0;
					int pred, vmin, vmax, *curr;

					vmax=NN; vmin=WW;
					if(NN<WW)vmin=NN, vmax=WW;
					pred=NN+WW-NNWW;
					CLAMP2(pred, vmin, vmax);

					curr=((int32_t*)HDR_DATA(src))+((src->iw*ky+kx)<<2|kc);
					*curr=(*curr+pred)<<(32-src->depth[kc])>>(32-src->depth[kc]);
				}
			}
			for(ky=1;ky<src->ih;ky+=2)
			{
				for(kx=1;kx<src->iw;kx+=2)	/*AV4B*/
				{
					int
						NW	=ky-1>=0&&kx-1>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx-1)<<2|kc]:0,
						NE	=ky-1>=0&&kx+1<src->iw		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx+1)<<2|kc]:0,
						SW	=ky+1<src->ih&&kx-1>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx-1)<<2|kc]:0,
						SE	=ky+1<src->ih&&kx+1<src->iw	?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+1)<<2|kc]:0;
					int pred, *curr;

					pred=(NW+NE+SW+SE+2)>>2;

					curr=((int32_t*)HDR_DATA(src))+((src->iw*ky+kx)<<2|kc);
					*curr=(*curr+pred)<<(32-src->depth[kc])>>(32-src->depth[kc]);
				}
			}
			for(ky=1;ky<src->ih;ky+=2)
			{
				for(kx=0;kx<src->iw;kx+=2)	/*AV12A*/
				{
					int
						N	=ky-1>=0			?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx+0)<<2|kc]:0,
						W	=kx-1>=0			?((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx-1)<<2|kc]:0,
						S	=ky+1<src->ih			?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+0)<<2|kc]:0,
						E	=kx+1<src->iw			?((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|kc]:0,
						NNW	=ky-2>=0&&kx-1>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-2)+kx-1)<<2|kc]:0,
						NNE	=ky-2>=0&&kx+1<src->iw		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-2)+kx+1)<<2|kc]:0,
						NWW	=ky-1>=0&&kx-2>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx-2)<<2|kc]:0,
						NEE	=ky-1>=0&&kx+2<src->iw		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx+2)<<2|kc]:0,
						SWW	=ky+1<src->ih&&kx-2>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx-2)<<2|kc]:0,
						SEE	=ky+1<src->ih&&kx+2<src->iw	?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+2)<<2|kc]:0,
						SSW	=ky+2<src->ih&&kx-1>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky+2)+kx-1)<<2|kc]:0,
						SSE	=ky+2<src->ih&&kx+1<src->iw	?((int32_t*)HDR_DATA(src))[(src->iw*(ky+2)+kx+1)<<2|kc]:0;
					int vmax, vmin, pred, *curr;

					vmax=N; vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>S)vmin=S;
					if(vmax<S)vmax=S;
					if(vmin>E)vmin=E;
					if(vmax<E)vmax=E;
					pred=(6*(N+W+S+E)-(NNW+NNE+NWW+NEE+SWW+SEE+SSW+SSE)+8)>>4;
					CLAMP2(pred, vmin, vmax);

					curr=((int32_t*)HDR_DATA(src))+((src->iw*ky+kx)<<2|kc);
					*curr=(*curr+pred)<<(32-src->depth[kc])>>(32-src->depth[kc]);
				}
			}
			for(ky=0;ky<src->ih;ky+=2)
			{
				for(kx=1;kx<src->iw;kx+=2)	/*AV8*/
				{
					int
						N	=ky-1>=0			?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx+0)<<2|kc]:0,
						W	=kx-1>=0			?((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx-1)<<2|kc]:0,
						S	=ky+1<src->ih			?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+0)<<2|kc]:0,
						E	=kx+1<src->iw			?((int32_t*)HDR_DATA(src))[(src->iw*(ky+0)+kx+1)<<2|kc]:0,
						NW	=ky-1>=0&&kx-1>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx-1)<<2|kc]:0,
						NE	=ky-1>=0&&kx+1<src->iw		?((int32_t*)HDR_DATA(src))[(src->iw*(ky-1)+kx+1)<<2|kc]:0,
						SW	=ky+1<src->ih&&kx-1>=0		?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx-1)<<2|kc]:0,
						SE	=ky+1<src->ih&&kx+1<src->iw	?((int32_t*)HDR_DATA(src))[(src->iw*(ky+1)+kx+1)<<2|kc]:0;
					int pred, *curr;

					pred=(3*(N+W+S+E)-(NW+NE+SW+SE)+4)>>3;

					curr=((int32_t*)HDR_DATA(src))+((src->iw*ky+kx)<<2|kc);
					*curr=(*curr+pred)<<(32-src->depth[kc])>>(32-src->depth[kc]);
				}
			}
		}
	}
}
void pred_CG420(Image *src, int fwd)
{
	int nch, fwdmask, bufsize, *pixels, *rows[4], ky, kx, kc, idx;

	fwdmask=-fwd;
	bufsize=(src->iw+16)*(int)sizeof(int[4*4]);/*4 padded rows * 4 channels max*/
	pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<3;++kc)
			{
				int NW, N, W, pred, curr;

				NW	=rows[1][kc-1*4];
				N	=rows[1][kc+0*4];
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

				curr=((int32_t*)HDR_DATA(src))[idx+kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=curr;

				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];

				((int32_t*)HDR_DATA(src))[idx+kc]=pred;
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
	int nch, fwdmask, bufsize, *pixels, *rows[4], ky, kx, kc, idx;

	fwdmask=-fwd;
	bufsize=(src->iw+16)*(int)sizeof(int[4*4]);/*4 padded rows * 4 channels max*/
	pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<3;++kc)
			{
				int NW, N, W, pred, curr;

				NW	=rows[1][kc-1*4];
				N	=rows[1][kc+0*4];
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

				curr=((int32_t*)HDR_DATA(src))[idx+kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=curr;

				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];

				((int32_t*)HDR_DATA(src))[idx+kc]=pred;
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
void pred_sel4(Image *src, int fwd)
{
	int nch, fwdmask, bufsize, *pixels, *rows[4], ky, kx, kc, idx;

	fwdmask=-fwd;
	bufsize=(src->iw+16)*(int)sizeof(int[4*4]);/*4 padded rows * 4 channels max*/
	pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<3;++kc)
			{
				int
					NNWW	=rows[2][kc-2*4],
					NN	=rows[2][kc+0*4],
					NNEE	=rows[2][kc+2*4],
					NW	=rows[1][kc-1*4],
					N	=rows[1][kc+0*4],
					NE	=rows[1][kc+1*4],
					WW	=rows[0][kc-2*4],
					W	=rows[0][kc-1*4];
				int curr, diff, mindiff, pred;

				mindiff=abs(W-WW);
				pred=W;

				diff=abs(N-NN);
				if(mindiff>diff)
					mindiff=diff, pred=N;

				diff=abs(NW-NNWW);
				if(mindiff>diff)
					mindiff=diff, pred=NW;

				diff=abs(NE-NNEE);
				if(mindiff>diff)
					mindiff=diff, pred=NE;

				curr=((int32_t*)HDR_DATA(src))[idx+kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=curr;

				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];

				((int32_t*)HDR_DATA(src))[idx+kc]=pred;
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
void pred_select(Image *src, int fwd)
{
	int amin[4], amax[4], nch, bufsize, *pixels, *rows[4], ky, kx, kc, idx;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	bufsize=(src->iw+16)*(int)sizeof(int[4*4]);/*4 padded rows * 4 channels max*/
	pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<3;++kc)
			{
				int NW, N, W, pred, curr, val;

				NW	=rows[1][kc-1*4];
				N	=rows[1][kc+0*4];
				W	=rows[0][kc-1*4];
				pred=abs(N-NW)>abs(W-NW)?N:W;/*WebP select (Paeth without NW)*/
				curr=((int32_t*)HDR_DATA(src))[idx+kc];
				if(g_dist>1)
				{
					if(fwd)
					{
						val=(curr-(int)pred+g_dist/2)/g_dist;
						curr=g_dist*val+(int)pred;
					}
					else
					{
						val=g_dist*curr+(int)pred;
						curr=val;
						CLAMP2(val, amin[kc], amax[kc]);
					}
				}
				else if(fwd)
				{
					val=curr-pred;
					val<<=32-src->depth[kc];
					val>>=32-src->depth[kc];
				}
				else
				{
					val=curr+pred;
					val<<=32-src->depth[kc];
					val>>=32-src->depth[kc];
					curr=val;
				}
				((int32_t*)HDR_DATA(src))[idx+kc]=val;

				rows[0][kc]=curr;
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	free(pixels);
}
void pred_paeth(Image *src, int fwd)
{
	int amin[4], amax[4], nch, bufsize, *pixels, *rows[4], ky, kx, kc, idx;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	bufsize=(src->iw+16)*(int)sizeof(int[4*4]);/*4 padded rows * 4 channels max*/
	pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<3;++kc)
			{
				int NW, N, W, pred, curr, p, score, s2, val;

				NW	=rows[1][kc-1*4];
				N	=rows[1][kc+0*4];
				W	=rows[0][kc-1*4];

				p=N+W-NW;
				pred=W; score=abs(pred-p);
				s2=abs(N-p);
				if(score>s2)pred=N, score=s2;
				s2=abs(NW-p);
				if(score>s2)pred=NW;

				curr=((int32_t*)HDR_DATA(src))[idx+kc];
				if(g_dist>1)
				{
					if(fwd)
					{
						val=(curr-(int)pred+g_dist/2)/g_dist;
						curr=g_dist*val+(int)pred;
					}
					else
					{
						val=g_dist*curr+(int)pred;
						curr=val;
						CLAMP2(val, amin[kc], amax[kc]);
					}
				}
				else if(fwd)
				{
					val=curr-pred;
					val<<=32-src->depth[kc];
					val>>=32-src->depth[kc];
				}
				else
				{
					val=curr+pred;
					val<<=32-src->depth[kc];
					val>>=32-src->depth[kc];
					curr=val;
				}
				((int32_t*)HDR_DATA(src))[idx+kc]=val;

				rows[0][kc]=curr;
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	free(pixels);
}
void pred_select4(Image *src, int fwd)
{
	enum
	{
		XPAD=8,
		NCH=4,
		NROWS=4,
		NVAL=3
	};
	int amin[4], amax[4], psize, invdist, uc0, vc0, vc1, ky, kx, kc, idx, yuv[3];
	int16_t *pixels, *rows[4];
	const uint8_t *perm;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	invdist=((1<<16)+g_dist-1)/g_dist;
	uc0=vc0=vc1=0;
	if(fwd)
		src->rct=crct2_analysis(src);
	perm=crct2_unpack(src->rct, &uc0, &vc0, &vc1);
	psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NCH*NROWS*NVAL]);/*4 padded rows * 4 channels max*/
	pixels=(int16_t*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		rows[2]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)2+NROWS)%NROWS)*NVAL;
		rows[3]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)3+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			yuv[0]=((int32_t*)HDR_DATA(src))[idx+perm[0]];
			yuv[1]=((int32_t*)HDR_DATA(src))[idx+perm[1]];
			yuv[2]=((int32_t*)HDR_DATA(src))[idx+perm[2]];
			for(kc=0;kc<4;++kc)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				rows[2]+=NROWS*NVAL;
				rows[3]+=NROWS*NVAL;
				if(kc==3)
					continue;
				{
					int
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL],

						xNW	=rows[1][1-1*NCH*NROWS*NVAL],
						xN	=rows[1][1+0*NCH*NROWS*NVAL],
						xW	=rows[0][1-1*NCH*NROWS*NVAL],

						yNW	=rows[1][2-1*NCH*NROWS*NVAL],
						yN	=rows[1][2+0*NCH*NROWS*NVAL],
						yW	=rows[0][2-1*NCH*NROWS*NVAL];
					int pred, curr, offset, e;
					
					offset=0;
					if(kc==1)offset=uc0*yuv[0]>>RCTBITS;
					if(kc==2)offset=(vc0*yuv[0]+vc1*yuv[1])>>RCTBITS;

					pred = xNW+xN+xW>yNW+yN+yW ? N : W;

					pred+=offset;
					CLAMP2(pred, amin[kc], amax[kc]);

					curr=yuv[kc];
					if(g_dist>1)
					{
						if(fwd)
						{
							curr-=pred;
							curr=(curr*invdist>>16)-(curr>>31);
							((int32_t*)HDR_DATA(src))[idx+kc]=curr;

							curr=g_dist*curr+pred;
						}
						else
							curr=g_dist*((int32_t*)HDR_DATA(src))[idx+kc]+pred;
						CLAMP2(curr, amin[kc], amax[kc]);
						yuv[kc]=curr;
					}
					else if(fwd)
					{
						e=curr-pred;
						e<<=32-src->depth[kc];
						e>>=32-src->depth[kc];
						((int32_t*)HDR_DATA(src))[idx+kc]=e;
					}
					else
					{
						curr=((int32_t*)HDR_DATA(src))[idx+kc]+pred;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						yuv[kc]=curr;
					}
					curr-=offset;
					rows[0][0]=curr;
					rows[0][1]=abs(curr-W);
					rows[0][2]=abs(curr-N);
				}
			}
			if(!fwd)
			{
				((int32_t*)HDR_DATA(src))[idx+perm[0]]=yuv[0];
				((int32_t*)HDR_DATA(src))[idx+perm[1]]=yuv[1];
				((int32_t*)HDR_DATA(src))[idx+perm[2]]=yuv[2];
			}
		}
	}
	free(pixels);
}


void pred_CG3D(Image *src, int fwd, int enable_ma)
{
	const int chpermutations[]=
	{
		0, 1, 2, 3,/*subR*/
		1, 2, 0, 3,/*subG*/
	};
	const int *perm;
	int nlevels[4], halfs[4], nch, nlumas, chpoints[3]={0}, fwdmask, psize, *pixels, *rows[2], ky, kx, kc, kc0, kk, idx;

	(void)enable_ma;
	nlevels[0]=1<<src->depth[0];
	nlevels[1]=1<<src->depth[1];
	nlevels[2]=1<<src->depth[2];
	nlevels[3]=1<<src->depth[3];
	halfs[0]=nlevels[0]>>1;
	halfs[1]=nlevels[1]>>1;
	halfs[2]=nlevels[2]>>1;
	halfs[3]=nlevels[3]>>1;
	nlumas=1;
	fwdmask=-fwd;
	psize=(src->iw+2)*(int)sizeof(int[2*4]);
	pixels=(int*)malloc(psize);/*2 padded rows * 4 channels max*/
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
#ifdef CG3D_ENABLE_MA
	for(kc=0;kc<4;++kc)
	{
		if(src->depth[kc])
			src->depth[kc]+=fwd;
	}
#else
	/*src->depth[1]-=src->depth[1]>src->src_depth[1];*/
#endif
	for(kc=1;kc<nch;++kc)
		nlumas+=src->depth[kc]==src->depth[0];
	chpoints[1]=nlumas;
	chpoints[2]=nch;
	perm=nlumas==nch?chpermutations+4:chpermutations;/*permute channels RGB->GBR if no RCT was applied earlier*/
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+(((src->iw+(intptr_t)2)*((ky-(intptr_t)0)&1)+1)<<2);
		rows[1]=pixels+(((src->iw+(intptr_t)2)*((ky-(intptr_t)1)&1)+1)<<2);
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kk=0;kk<2;++kk)
			{
				for(kc0=chpoints[kk];kc0<chpoints[kk+1];++kc0)
				{
					int NW, N, W, offset, pred, curr;

					kc=perm[kc0];
					NW	=rows[1][kc-4];
					N	=rows[1][kc+0];
					W	=rows[0][kc-4];
					offset	=0;
					if(kc0>chpoints[kk])
						offset+=rows[0][perm[chpoints[kk]]];

					MEDIAN3_32(pred, N, W, N+W-NW);

					pred+=offset;
					CLAMP2(pred, -halfs[kc], halfs[kc]-1);

					curr=((int32_t*)HDR_DATA(src))[idx+kc];
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=curr;

					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];

					((int32_t*)HDR_DATA(src))[idx+kc]=pred;
					rows[0][kc]=(fwd?curr:pred)-offset;
				}
			}
			rows[0]+=4;
			rows[1]+=4;
		}
	}
#ifdef CG3D_ENABLE_MA
	for(kc=0;kc<4;++kc)
	{
		if(src->depth[kc])
			src->depth[kc]-=!fwd;
	}
#else
	/*src->depth[1]+=src->depth[1]<24;*/
#endif
	free(pixels);
}


void pred_table(Image *src, int fwd)
{
	int nlevels[4], halfs[4], nch, maxdepth, tablesize, *table, psize, *pixels, *rows[2], fwdmask, ky, kx, kc, idx;

	nlevels[0]=1<<src->depth[0];
	nlevels[1]=1<<src->depth[1];
	nlevels[2]=1<<src->depth[2];
	nlevels[3]=1<<src->depth[3];
	halfs[0]=nlevels[0]>>1;
	halfs[1]=nlevels[1]>>1;
	halfs[2]=nlevels[2]>>1;
	halfs[3]=nlevels[3]>>1;
	fwdmask=-fwd;
	psize=(src->iw+2)*(int)sizeof(int[2*4]);
	pixels=(int*)malloc(psize);/*2 padded rows * 4 channels max*/
	maxdepth=src->depth[0];
	UPDATE_MAX(maxdepth, src->depth[1]);
	UPDATE_MAX(maxdepth, src->depth[2]);
	UPDATE_MAX(maxdepth, src->depth[3]);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	tablesize=(int)sizeof(int[4])*nch<<maxdepth;/*{N, W, NE, NW} tables * nch*/
	table=(int*)malloc(tablesize);
	if(!pixels||!table)
	{
		CRASH("Alloc error");
		return;
	}
	memset(table, 0, tablesize);
	memset(pixels, 0, psize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+(((src->iw+(intptr_t)2)*((ky-(intptr_t)0)&1)+1)<<2);
		rows[1]=pixels+(((src->iw+(intptr_t)2)*((ky-(intptr_t)1)&1)+1)<<2);
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<nch;++kc)
			{
				int
					NW	=rows[1][kc-4],
					N	=rows[1][kc+0],
					NE	=rows[1][kc+1],
					W	=rows[0][kc-4];
				int
					*curr_table=table+((intptr_t)4*kc<<maxdepth),
					*tNW	=curr_table+(0<<maxdepth|((NW	+halfs[kc])&(nlevels[kc]-1))),
					*tN	=curr_table+(1<<maxdepth|((N	+halfs[kc])&(nlevels[kc]-1))),
					*tNE	=curr_table+(2<<maxdepth|((NE	+halfs[kc])&(nlevels[kc]-1))),
					*tW	=curr_table+(3<<maxdepth|((W	+halfs[kc])&(nlevels[kc]-1)));
				int
					NW2	=*tNW	,/*==(1<<maxdepth)?NW	:*tNW,	*/
					N2	=*tN	,/*==(1<<maxdepth)?N	:*tN,	*/
					NE2	=*tNE	,/*==(1<<maxdepth)?NE	:*tNE,	*/
					W2	=*tW	;/*==(1<<maxdepth)?W	:*tW;	*/
				int pred, curr;

				pred=(N2+W2+NE2+NW2+2)>>2;

				curr=((int32_t*)HDR_DATA(src))[idx+kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=curr;

				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];

				((int32_t*)HDR_DATA(src))[idx+kc]=pred;
				rows[0][kc]=(fwd?curr:pred);

				curr=rows[0][kc];
				*tNW	=curr;
				*tN	=curr;
				*tNE	=curr;
				*tW	=curr;
			}
			rows[0]+=4;
			rows[1]+=4;
		}
	}
	free(pixels);
	free(table);
}


	#define PU_UPDATE_LUMA

void pred_PU(Image *src, int fwd)
{
	const int perm[]={1, 2, 0, 3};
	int nlevels[4], halfs[4], nch, fwdmask, psize, *pixels, *rows[4], ky, kx, idx, kc0;

	nlevels[0]=1<<src->depth[0];
	nlevels[1]=1<<src->depth[1];
	nlevels[2]=1<<src->depth[2];
	nlevels[3]=1<<src->depth[3];
	halfs[0]=nlevels[0]>>1;
	halfs[1]=nlevels[1]>>1;
	halfs[2]=nlevels[2]>>1;
	halfs[3]=nlevels[3]>>1;
	fwdmask=-fwd;
	psize=(src->iw+4)*(int)sizeof(int[4*4*2]);/*4 padded rows * 4 channels max * {pixels, errors}*/
	pixels=(int*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+(((src->iw+(intptr_t)4)*((ky-(intptr_t)0)&3)+2)<<3);
		rows[1]=pixels+(((src->iw+(intptr_t)4)*((ky-(intptr_t)1)&3)+2)<<3);
		rows[2]=pixels+(((src->iw+(intptr_t)4)*((ky-(intptr_t)2)&3)+2)<<3);
		rows[3]=pixels+(((src->iw+(intptr_t)4)*((ky-(intptr_t)3)&3)+2)<<3);
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
#if 1
#ifdef PU_UPDATE_LUMA
			if(!fwd)
			{
				int *ptr=((int32_t*)HDR_DATA(src))+idx, temp;
				ROTATE3(ptr[2], ptr[1], ptr[0], temp);
				{
					int
						cb	=((int32_t*)HDR_DATA(src))[idx+2],
						cr	=((int32_t*)HDR_DATA(src))[idx+0];
					int update=(cb+cr)>>2;
					int luma=((int32_t*)HDR_DATA(src))[idx+1]-update;/*subtract update*/
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					((int32_t*)HDR_DATA(src))[idx+1]=luma;
				}
			}
#endif
			for(kc0=0;kc0<src->nch;++kc0)
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
				(void)WW;
				(void)W;
				(void)eNNWW;
				(void)eNNW;
				(void)eNN;
				(void)eNNE;
				(void)eNNEE;
				(void)eNWW;
				(void)eNW;
				(void)eN;
				(void)eNE;
				(void)eNEE;
				(void)eWW;
				(void)eW;
				{
					int
						vx=abs(eW-eWW)+abs(eN-eNW)+abs(eNE-eN)+1,
						vy=abs(eW-eNW)+abs(eN-eNN)+abs(eNE-eNNE)+1;
					int update=vx<vy?eN:eW;
					update/=4;
					MEDIAN3_32(pred, N, W, N+W-NW);
					pred+=update;
				}
				if(kc0>0)
				{
					offset=rows[0][1];
					pred+=offset;
					CLAMP2(pred, -halfs[kc], halfs[kc]-1);
				}
				{
					int curr=((int32_t*)HDR_DATA(src))[idx+kc];
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=curr;
				
					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];

					((int32_t*)HDR_DATA(src))[idx+kc]=pred;
					rows[0][kc+0]=(fwd?curr:pred)-offset;
					rows[0][kc+4]=(fwd?pred:curr);
				}
			}
#ifdef PU_UPDATE_LUMA
			if(fwd)
			{
				int
					cb	=((int32_t*)HDR_DATA(src))[idx+2],
					cr	=((int32_t*)HDR_DATA(src))[idx+0];
				int update=(cb+cr)>>2;
				{
					int luma=((int32_t*)HDR_DATA(src))[idx+1]+update;/*add update*/
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					((int32_t*)HDR_DATA(src))[idx+1]=luma;
				}
				{
					int *ptr=((int32_t*)HDR_DATA(src))+idx, temp;
					ROTATE3(ptr[0], ptr[1], ptr[2], temp);
				}
			}
#endif
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
void pred_divfreeWP(Image *src, int fwd)/*not DIV-free*/
{
	int nch, fwdmask, bufsize, *pixels, *rows[4], ky, kx, kc, idx, preds[DFWP_NPREDS*2]={0};

	fwdmask=-fwd;
	bufsize=(src->iw+16)*(int)sizeof(int[4*4*DFWP_NPREDS*2]);/*4 padded rows * 4 channels max * {pixel, error, 2 perrors, 4 perrors, 8 perrors}*/
	pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*DFWP_NPREDS*2*4;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*DFWP_NPREDS*2*4;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*DFWP_NPREDS*2*4;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*DFWP_NPREDS*2*4;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<src->nch;++kc)
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
				int pred, j, nsrc, srcidx, kp, c2, k;

				pred=0;
				j=DFWP_NPREDS;

				preds[j++]=N[0];
				preds[j++]=W[0];
				preds[j++]=3*(N[0]-NN[0])+NNN[0];
				preds[j++]=3*(W[0]-WW[0])+WWW[0];
				preds[j++]=N[0]+W[0]-NW[0];
				preds[j++]=N[0]+NE[0]-NNE[0];
				preds[j++]=W[0]+NE[0]-N[0];
				preds[j++]=(W[0]+NEEE[0])>>1;

				for(nsrc=DFWP_NPREDS, srcidx=DFWP_NPREDS;nsrc>1;)
				{
					int ndst=nsrc>>1, dstidx=srcidx-ndst;
					for(kp=0;kp<ndst;++kp)
					{
						int w1=NNWW[srcidx+kp*2+0]+NNW[srcidx+kp*2+0]+NN[srcidx+kp*2+0]+NNE[srcidx+kp*2+0]+NNEE[srcidx+kp*2+0]+NWW[srcidx+kp*2+0]+NW[srcidx+kp*2+0]+N[srcidx+kp*2+0]*3+NE[srcidx+kp*2+0]+NEE[srcidx+kp*2+0]+WW[srcidx+kp*2+0]+W[srcidx+kp*2+0]*3+1;
						int w2=NNWW[srcidx+kp*2+1]+NNW[srcidx+kp*2+1]+NN[srcidx+kp*2+1]+NNE[srcidx+kp*2+1]+NNEE[srcidx+kp*2+1]+NWW[srcidx+kp*2+1]+NW[srcidx+kp*2+1]+N[srcidx+kp*2+1]*3+NE[srcidx+kp*2+1]+NEE[srcidx+kp*2+1]+WW[srcidx+kp*2+1]+W[srcidx+kp*2+1]*3+1;

						preds[dstidx+kp]=(preds[srcidx+kp*2+0]*w2+preds[srcidx+kp*2+1]*w1)/(w1+w2);/*NPREDs DIVs*/
					}
					nsrc=ndst;
					srcidx=dstidx;
				}
				pred=preds[1];
				
				c2=((int32_t*)HDR_DATA(src))[idx+kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=c2;
			
				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];

				((int32_t*)HDR_DATA(src))[idx+kc]=pred;
				curr[0]=fwd?c2:pred;
				curr[1]=fwd?pred:c2;
				for(k=2;k<DFWP_NPREDS*2;++k)
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


/*
https://github.com/WangXuan95/NBLIC-Image-Compression
*/
void pred_nblic(Image *src, int fwd)
{
	const int thresholds[]={1, 3, 9, 20, 50, 110, 300, 800};
	const int q_mid[]={0, 2, 4, 7, 10, 14, 20, 26, 34, 42, 52, 64, 78, 95, 135, 200};
	int amin[4], amax[4], nlevels[4], halfs[4], nch, invdist, bufsize, *pixels, *rows[4], ssesize, *sse, ky, kx, kc, idx;

	nlevels[0]=1<<src->depth[0];
	nlevels[1]=1<<src->depth[1];
	nlevels[2]=1<<src->depth[2];
	nlevels[3]=1<<src->depth[3];
	halfs[0]=nlevels[0]>>1;
	halfs[1]=nlevels[1]>>1;
	halfs[2]=nlevels[2]>>1;
	halfs[3]=nlevels[3]>>1;
	amin[0]=-halfs[0];
	amin[1]=-halfs[1];
	amin[2]=-halfs[2];
	amin[3]=-halfs[3];
	amax[0]=halfs[0]-1;
	amax[1]=halfs[1]-1;
	amax[2]=halfs[2]-1;
	amax[3]=halfs[3]-1;
	invdist=((1<<16)+g_dist-1)/g_dist;
	bufsize=(src->iw+16)*(int)sizeof(int[4*4*2]);/*4 padded rows * 4 channels max*/
	pixels=(int*)malloc(bufsize);
	ssesize=(int)sizeof(int[4<<4>>1<<8]);
	sse=(int*)malloc(ssesize);
	if(!pixels||!sse)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	memset(sse, 0, ssesize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4*2;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4*2;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4*2;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4*2;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<src->nch;++kc)
			{
				int
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
					eN	=rows[1][kc+0*4*2+4],
					eW	=rows[0][kc-1*4*2+4];
				int pred_linear, pred_ang, cost, csum, cmin, wt, pred, delta, qd, qu, qv, qw, ctx, *cell, corr, curr, e, e2;

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
				(void)WW;
				(void)W;
				(void)eN;
				(void)eW;
				/*
				NNWW NNW NN NNE NNEE
				NWW  NW  N  NE  NEE
				WW   W   ?
				*/
				/*
				|simple predict:
				|
				|      -1
				|   -2  9  2
				|-1  9  ?
				*/
				pred_linear=9*(N+W)+2*(NE-NW)-NN-WW;
				pred_ang=0;
				csum=0;
				cmin=0xFFFFFF;
				CLAMP2(pred_linear, -halfs[kc]<<4, halfs[kc]<<4);

				cost=2*(abs(W-WW)+abs(NW-NWW)+abs(N-NW)+abs(NE-N));/*180 degrees*/
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=2*W;
				}

				cost=2*(abs(W-NW)+abs(NW-NNW)+abs(N-NN)+abs(NE-NNE));/*90 degrees*/
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=2*N;
				}

				cost=2*(abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW)+abs(NE-NN));/*135 degrees*/
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=2*NW;
				}

				cost=2*(abs(W-N)+abs(NW-NN)+abs(N-NNE)+abs(NE-NNEE));/*45 degrees*/
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=2*NE;
				}

				cost=abs(2*W-WW-NWW)+abs(2*NW-NNW-NNWW)+abs(2*N-NW-NNW)+abs(2*NE-N-NN);/*180-atand(1/2) ~= 153 degrees*/
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=W+NW;
				}

				cost=abs(2*W-NWW-NW)+abs(2*NW-NNWW-NNW)+abs(2*N-NNW-NN)+abs(2*NE-NN-NNE);/*90+atand(1/2) ~= 116 degrees*/
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=NW+N;
				}

				cost=abs(2*W-NW-N)+abs(2*NW-NNW-NN)+abs(2*N-NN-NNE)+abs(2*NE-NNE-NNEE);/*90-atand(1/2) ~= 63 degrees*/
				csum+=cost;
				if(cmin>cost)
				{
					cmin=cost;
					pred_ang=N+NE;
				}
				
				csum-=7*cmin;
				for(wt=0;wt<(int)_countof(thresholds)&&(thresholds[wt]<<src->depth[kc])-thresholds[wt]<=(csum<<3);++wt);
				pred=(8*wt*pred_ang+(8-wt)*pred_linear+64)>>7;

				/*get quantized delta:*/
				delta=abs(W-WW)+abs(N-NW)+abs(N-NE)+abs(W-NW)+abs(N-NN)+abs(NE-NNE)+2*abs(eW);
				for(qd=0;qd<(int)_countof(q_mid)-1&&delta>q_mid[qd];++qd);
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

				/*get context address:*/
				ctx=kc;
				ctx=ctx<<(4-1)|qu>>1;
				ctx=ctx<<1|(pred>2*N-NN);
				ctx=ctx<<1|(pred>2*W-WW);
				ctx=ctx<<1|(pred>NN);
				ctx=ctx<<1|(pred>WW);
				ctx=ctx<<1|(pred>NE);
				ctx=ctx<<1|(pred>NW);
				ctx=ctx<<1|(pred>N);
				ctx=ctx<<1|(pred>W);

				/*correct prediction by context:*/
				cell=sse+ctx;
				corr=*cell;
				pred+=corr>>8;
				CLAMP2(pred, -halfs[kc], halfs[kc]);

				/*apply prediction:*/
				curr=((int32_t*)HDR_DATA(src))[idx];
				if(g_dist>1)
				{
					if(fwd)
					{
						curr-=pred;
						curr=(curr*invdist>>16)-(curr>>31&-(g_dist>1));/*curr/=g_dist*/
						((int32_t*)HDR_DATA(src))[idx]=curr;

						curr=g_dist*curr+pred;
						CLAMP2(curr, amin[kc], amax[kc]);
					}
					else
					{
						curr=g_dist*curr+pred;
						CLAMP2(curr, amin[kc], amax[kc]);

						((int32_t*)HDR_DATA(src))[idx]=curr;
					}
				}
				else
				{
					if(fwd)
					{
						e=curr-pred;
						e<<=32-src->depth[kc];
						e>>=32-src->depth[kc];
						((int32_t*)HDR_DATA(src))[idx]=e;
					}
					else
					{
						curr+=pred;
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
						((int32_t*)HDR_DATA(src))[idx]=curr;
					}
				}
				rows[0][kc+0]=curr;
				rows[0][kc+4]=curr-pred;
				/*update context:*/
				e2=rows[0][kc+4];
				*cell=corr+(((e2<<7)-corr+(1<<7>>1))>>7);/*FIXME compare with traditional SSE*/
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


#ifdef WGRAD_SSE
static int quantize_nb(int x)
{
	int ax=abs(x), negmask=x<0;
	x=31^LZCNT32(ax);
	x^=negmask;
	x-=negmask;
	return x;
}
#define QUANTIZE_NB(X) quantize_nb(X)+maxdepth
#endif
void pred_wgrad(Image *src, int fwd, int hasRCT)
{
	int nch, fwdmask, psize, *pixels, *rows[4], ky, kx, kc, idx;

	fwdmask=-fwd;
	psize=(src->iw+16)*(int)sizeof(int[4*4*3]);/*4 padded rows * 4 channels max * {pixels, gx, gy}*/
	pixels=(int*)malloc(psize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4*3;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4*3;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4*3;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4*3;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<nch;++kc)
			{
				int
					NNNN	=rows[0][kc+0*4*3],
					NNN	=rows[3][kc+0*4*3],
					NNWW	=rows[2][kc-2*4*3],
					NNW	=rows[2][kc-1*4*3],
					NN	=rows[2][kc+0*4*3],
					NNE	=rows[2][kc+1*4*3],
					NNEE	=rows[2][kc+2*4*3],
					NWW	=rows[1][kc-2*4*3],
					NW	=rows[1][kc-1*4*3],
					N	=rows[1][kc+0*4*3],
					NE	=rows[1][kc+1*4*3],
					NEE	=rows[1][kc+2*4*3],
					NEEE	=rows[1][kc+3*4*3],
					NEEEE	=rows[1][kc+4*4*3],
					NEEEEE	=rows[1][kc+5*4*3],
					WWWWW	=rows[0][kc-5*4*3],
					WWWW	=rows[0][kc-4*4*3],
					WWW	=rows[0][kc-3*4*3],
					WW	=rows[0][kc-2*4*3],
					W	=rows[0][kc-1*4*3];
				int
					gy=abs(N-NN)+abs(W-NW)+abs(NE-NNE)+abs(NW-NNW)+1,
					gx=abs(N-NW)+abs(W-WW)+abs(NE-N)+abs(NW-NWW)+1;
				int pred;

#if 1
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
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
#endif
				pred=(gx*N+gy*W)/(gx+gy);
				{
					int curr=((int32_t*)HDR_DATA(src))[idx+kc];

					{
						pred^=fwdmask;
						pred-=fwdmask;
						pred+=curr;
						pred<<=32-src->depth[kc];
						pred>>=32-src->depth[kc];
					}
					((int32_t*)HDR_DATA(src))[idx+kc]=pred;
					if(!fwd)
						curr=pred;
					rows[0][kc+0*4]=curr;
				}
			}
			rows[0]+=4*3;
			rows[1]+=4*3;
			rows[2]+=4*3;
			rows[3]+=4*3;
		}
	}
	free(pixels);

	(void)hasRCT;
}

	#define WG2_USE_FP64
#if 0
#	define WG2_DIVFREE
#	define WG2_NPREDS 8
#endif
void pred_wgrad2(Image *src, int fwd)
{
	const int perm[]={1, 2, 0, 3};
	int nch, fwdmask, psize, *pixels, *rows[4], k, ky, kx, kc0, idx;

	fwdmask=-fwd;
	psize=(src->iw+16)*sizeof(int[4*4]);
	pixels=(int*)malloc(psize);/*4 padded rows * 4 channels max*/
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int eprev[16*4]={0};

		rows[0]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)<<2);
		rows[1]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)<<2);
		rows[2]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)<<2);
		rows[3]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)<<2);
#ifdef WG2_DIVFREE
		int mixinCDF[WG2_NPREDS]={0}, tempCDF[WG2_NPREDS+1]={0};
		for(k=0;k<WG2_NPREDS+1;++k)
		{
			int val=0x10000*k/WG2_NPREDS;
			eprev[k<<2|0]=val;
			eprev[k<<2|1]=val;
			eprev[k<<2|2]=val;
			eprev[k<<2|3]=val;
		}
		for(k=0;k<WG2_NPREDS;++k)
			mixinCDF[k]=(0x10000-16)*k/(WG2_NPREDS-1);
#endif
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc0=0;kc0<src->nch;++kc0)
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
				int preds[8], j;

#if 1
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
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
#endif
				j=0;
				preds[j++]=N+W-NW;
				preds[j++]=W+NE-N;
				preds[j++]=N+NE-NNE;
				preds[j++]=N;
				preds[j++]=W;
				preds[j++]=3*(N-NN)+NNN;
				preds[j++]=3*(W-WW)+WWW;
				preds[j++]=(W+NEEE)>>1;

#ifdef WG2_DIVFREE
				for(k=0;k<(int)_countof(preds);++k)
				{
					int weight=eprev[(k+1)<<2|kc]-eprev[k<<2|kc];
					lpred+=(int64_t)weight*preds[k];
				}
				lpred>>=16;
#elif defined WG2_USE_FP64
				{
					double fpred=0, fsum=0;
					for(k=0;k<(int)_countof(preds);++k)
					{
						double weight=1./(eprev[k<<2|kc]+1);
						fpred+=weight*preds[k];
						fsum+=weight;
					}
					fpred/=fsum;
					pred=(int)TRUNC64(fpred);
				}
#else
				for(k=0;k<(int)_countof(preds);++k)
				{
					int weight=0x1000000/(eprev[k<<2|kc]+1);
					lpred+=(int64_t)weight*preds[k];
					wsum+=weight;
				}
				lpred/=wsum+1;
#endif
#ifndef WG2_USE_FP64
				pred=(int)lpred;
#endif
				{
					int curr=((int32_t*)HDR_DATA(src))[idx+kc];
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=curr;

					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];

					((int32_t*)HDR_DATA(src))[idx+kc]=pred;
					rows[0][kc]=(fwd?curr:pred)-offset;
				}
				{
					int curr=rows[0][kc];
#ifdef WG2_DIVFREE
					memset(tempCDF, 0, sizeof(tempCDF));
					for(k=0;k<(int)_countof(preds);++k)
					{
						int err, *ptr, k2;

						err=abs(curr-preds[k]);
						ptr=mixinCDF;
						for(k2=0;k2<(int)_countof(preds);++k2)
						{
							tempCDF[k2]+=(*ptr-tempCDF[k2])*err>>src->depth[kc];
							ptr+=k2!=k;
						}
					}
					for(k=0;k<(int)_countof(preds);++k)
					{
						int update=tempCDF[k]>>3;
						eprev[k<<2|kc]=update;
					}
#else
					for(k=0;k<(int)_countof(preds);++k)
					{
						int err=abs(curr-preds[k]);
						eprev[k<<2|kc]+=err;
						eprev[k<<2|kc]-=(eprev[k<<2|kc]+3)>>2;
					}
#endif
				}
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	free(pixels);
}
void pred_wgrad3(Image *src, int fwd, int hasRCT)
{
	const int perm0[]={1, 2, 0, 3, 0, 1, 2, 3}, *perm;
	int nlevels[4], halfs[4], nch, fwdmask, psize, *pixels, *rows[4], ky, kx, kc0, idx;

	nlevels[0]=1<<src->depth[0];
	nlevels[1]=1<<src->depth[1];
	nlevels[2]=1<<src->depth[2];
	nlevels[3]=1<<src->depth[3];
	halfs[0]=nlevels[0]>>1;
	halfs[1]=nlevels[1]>>1;
	halfs[2]=nlevels[2]>>1;
	halfs[3]=nlevels[3]>>1;
	perm=perm0+hasRCT*4;
	fwdmask=-fwd;
	psize=(src->iw+16)*(int)sizeof(int[4*4]);
	pixels=(int*)malloc(psize);/*4 padded rows * 4 channels max*/
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int eprev[16*4]={0};

		rows[0]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)<<2);
		rows[1]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)<<2);
		rows[2]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)<<2);
		rows[3]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)<<2);
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			if(!hasRCT&&!fwd)
			{
				int
					cb	=((int32_t*)HDR_DATA(src))[idx+2],
					cr	=((int32_t*)HDR_DATA(src))[idx+0];
				int update=(cb+cr)>>2;
				update=CLAMP(-halfs[1], update, halfs[1]-1);
				{
					int luma=((int32_t*)HDR_DATA(src))[idx+1]-update;/*subtract update*/
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					((int32_t*)HDR_DATA(src))[idx+1]=luma;
				}
			}
			for(kc0=0;kc0<src->nch;++kc0)
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
				double fpred=0, fsum=0;
				int preds[8], j, k;

#if 1
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
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
#endif
				j=0;
				preds[j++]=N+W-NW;
				preds[j++]=W+NE-N;
				preds[j++]=N+NE-NNE;
				preds[j++]=N;
				preds[j++]=W;
				preds[j++]=3*(N-NN)+NNN;
				preds[j++]=3*(W-WW)+WWW;
				preds[j++]=(W+NEEE)>>1;
				for(k=0;k<(int)_countof(preds);++k)
				{
					double weight=1./(eprev[k<<2|kc]+1);
					fpred+=weight*preds[k];
					fsum+=weight;
				}
				fpred/=fsum;
				pred=(int)TRUNC64(fpred);
				if(!hasRCT&&kc0>0)
				{
					offset+=rows[0][1];
					if(kc0>1)
						offset=(2*offset+rows[0][2])>>1;
					pred+=offset;
					pred=CLAMP(-halfs[kc], pred, halfs[kc]-1);
				}
				{
					int curr=((int32_t*)HDR_DATA(src))[idx+kc];
					pred^=fwdmask;
					pred-=fwdmask;
					pred+=curr;

					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];

					((int32_t*)HDR_DATA(src))[idx+kc]=pred;
					rows[0][kc]=(fwd?curr:pred)-offset;
				}
				{
					int curr=rows[0][kc];
					for(k=0;k<(int)_countof(preds);++k)
					{
						int err=abs(curr-preds[k]);
						eprev[k<<2|kc]+=err;
						eprev[k<<2|kc]=eprev[k<<2|kc]*0xE6>>8;
					}
				}
			}
			if(!hasRCT&&fwd)
			{
				int
					cb	=((int32_t*)HDR_DATA(src))[idx+2],
					cr	=((int32_t*)HDR_DATA(src))[idx+0];
				int update=(cb+cr)>>2;
				update=CLAMP(-halfs[1], update, halfs[1]-1);
				{
					int luma=((int32_t*)HDR_DATA(src))[idx+1]+update;/*add update*/
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					((int32_t*)HDR_DATA(src))[idx+1]=luma;
				}
			}
			rows[0]+=4;
			rows[1]+=4;
			rows[2]+=4;
			rows[3]+=4;
		}
	}
	free(pixels);
}


#define LWAV_NPARAMS 5	/*3*/
#define LWAV_PBITS 8
#define LWAV_NBITS 16
void pred_lwav(Image *src, int fwd)
{
	int nch, fwdmask, bufsize, *pixels, *rows[4], params[4][LWAV_NPARAMS]={0}, ky, kx, kc, idx;

	fwdmask=-fwd;
	bufsize=(src->iw+16)*(int)sizeof(int[4*4*(2+LWAV_NPARAMS)]);
	pixels=(int*)malloc(bufsize);/*4 padded rows * 4 channels max * {pixels, errors}*/
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	FILLMEM(params[0], 1<<LWAV_NBITS, sizeof(int[LWAV_NPARAMS]), sizeof(int));
	FILLMEM(params[1], 1<<LWAV_NBITS, sizeof(int[LWAV_NPARAMS]), sizeof(int));
	FILLMEM(params[2], 1<<LWAV_NBITS, sizeof(int[LWAV_NPARAMS]), sizeof(int));
	FILLMEM(params[3], 1<<LWAV_NBITS, sizeof(int[LWAV_NPARAMS]), sizeof(int));
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int pred=0;

		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4*(2+LWAV_NPARAMS);
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4*(2+LWAV_NPARAMS);
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4*(2+LWAV_NPARAMS);
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4*(2+LWAV_NPARAMS);
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<src->nch;++kc)
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
				int64_t lpred;
				int ctx[10], e, j, k, wsum;

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
				(void)WWWWW;
				(void)WWWW;
				(void)WWW;
				(void)WW;
				(void)W;
				(void)curr;

				j=0;
				ctx[j++]=0x2100; ctx[j++]=N[0]+W[0]-NW[0];
				ctx[j++]=0x0DB8; ctx[j++]=W[0]+NE[0]-N[0];
				ctx[j++]=0x0E22; ctx[j++]=N[0]-((N[1]+W[1]+NE[1])*-0x005C>>8);
				ctx[j++]=0x181F; ctx[j++]=W[0]-((N[1]+W[1]+NW[1])*-0x005B>>8);
				ctx[j++]=0x0BF3; ctx[j++]=N[0]-((0x00DF*NW[1]+0x0051*N[1]+0x00BD*NE[1]+0x005C*(NN[0]-N[0])-0x0102*(NW[0]-W[0]))>>8);
				MEDIAN3_32(ctx[0<<1|1], N[0], W[0], ctx[0<<1|1]);
				lpred=0;
				wsum=0;
				for(k=0;k<LWAV_NPARAMS;++k)
				{
					int weight=params[kc][k]*ctx[k<<1|0]/(NW[k+2]+N[k+2]+NE[k+2]+WW[k+2]+W[k+2]+1);
					lpred+=(int64_t)weight*ctx[k<<1|1];
					wsum+=weight;
				}
				if(wsum)
					pred=(int)((lpred+(wsum>>1))/wsum);
				else
					pred=ctx[0<<1|1];

				CLAMP3_32(pred, pred, N[0], W[0], NE[0]);

				{
					int p2=(pred+(1<<LWAV_PBITS>>1))>>LWAV_PBITS;
					int c=((int32_t*)HDR_DATA(src))[idx+kc];
					p2^=fwdmask;
					p2-=fwdmask;
					p2+=c;

					p2<<=32-src->depth[kc];
					p2>>=32-src->depth[kc];

					((int32_t*)HDR_DATA(src))[idx+kc]=p2;
					curr[0]=(fwd?c:p2)<<LWAV_PBITS;
					curr[1]=curr[0]-pred;
				}

				e=(2*W[1]+curr[1]+N[1]+NE[1]+NEEE[1]);
				for(k=0;k<LWAV_NPARAMS;++k)
				{
					curr[k+2]=abs(curr[0]-ctx[k<<1|1]);
					params[kc][k]+=(int)((int64_t)e*ctx[k<<1|1]>>src->depth[kc]);
					CLAMP2(params[kc][k], -(1<<LWAV_NBITS), 1<<LWAV_NBITS);
				}
			}
			rows[0]+=4*(2+LWAV_NPARAMS);
			rows[1]+=4*(2+LWAV_NPARAMS);
			rows[2]+=4*(2+LWAV_NPARAMS);
			rows[3]+=4*(2+LWAV_NPARAMS);
		}
	}
	free(pixels);
}


#if 0
#	define WPU_TRACE
#endif
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
	WPU_PRED(N+NE-NNE)\

#if 0
	WPU_PRED(grad)\
	WPU_PRED((W+NEE)>>1)\
	WPU_PRED((3*W+NEEE)>>2)\
	WPU_PRED(((abs(eW)+1)*N+(abs(eN)+1)*W)/(abs(eN)+abs(eW)+2))\
	WPU_PRED(((abs(eNW)+1)*NE+(abs(eNE)+1)*NW)/(abs(eNE)+abs(eNW)+2))\
	WPU_PRED(median1)\
	WPU_PRED(median2)\
	WPU_PRED(cgrad)\

#endif
#define WPU_NBITS 8
#define WPU_PIXEL_STRIDE (4*(WPU_NPREDS+(intptr_t)2))
void pred_WPU(Image *src, int fwd)
{
	enum
	{
#define WPU_PRED(E) +1
		NPREDS2=WPU_PREDLIST
#undef  WPU_PRED
	};
	const int perm[]={1, 2, 0, 3};
	int nlevels[4], halfs[4], boosthalfs[4], nch, bufsize, *pixels, *rows[4], k, ky, kx, kc0, idx;
	int eprev[WPU_NPREDS]={0}, weights[WPU_NPREDS<<2]={0};
#ifdef WPU_TRACE
	int64_t trace[WPU_NPREDS]={0};
#endif

	nlevels[0]=1<<src->depth[0],
	nlevels[1]=1<<src->depth[1],
	nlevels[2]=1<<src->depth[2],
	nlevels[3]=1<<src->depth[3],
	halfs[0]=nlevels[0]>>1;
	halfs[1]=nlevels[1]>>1;
	halfs[2]=nlevels[2]>>1;
	halfs[3]=nlevels[3]>>1;
	boosthalfs[0]=nlevels[0]<<WPU_NBITS>>1;
	boosthalfs[1]=nlevels[1]<<WPU_NBITS>>1;
	boosthalfs[2]=nlevels[2]<<WPU_NBITS>>1;
	boosthalfs[3]=nlevels[3]<<WPU_NBITS>>1;
	for(k=0;k<(WPU_NPREDS<<2);++k)
		weights[k]=0x8000;
	bufsize=(src->iw+8)*(int)sizeof(int[4*WPU_PIXEL_STRIDE]);
	pixels=(int*)malloc(bufsize);/*4 padded rows * 4 channels max * {pixels, errors, pred_errors...}*/
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);/*{...padding    Y U V A  eY eU eV eA  a0Y a0U a0V a0A  a1Y a1U a1V a1A ...} * 4 rows*/
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)8)*((ky-(intptr_t)0)&3)+4)*WPU_PIXEL_STRIDE;
		rows[1]=pixels+((src->iw+(intptr_t)8)*((ky-(intptr_t)1)&3)+4)*WPU_PIXEL_STRIDE;
		rows[2]=pixels+((src->iw+(intptr_t)8)*((ky-(intptr_t)2)&3)+4)*WPU_PIXEL_STRIDE;
		rows[3]=pixels+((src->iw+(intptr_t)8)*((ky-(intptr_t)3)&3)+4)*WPU_PIXEL_STRIDE;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
#ifdef WPU_UPDATE_LUMA
			if(!fwd)
			{
				int
					cb	=((int32_t*)HDR_DATA(src))[idx+2],
					cr	=((int32_t*)HDR_DATA(src))[idx+0];
				int update=(cb+cr)>>3;
				update=CLAMP(-halfs[1], update, halfs[1]-1);
				{
					int luma=((int32_t*)HDR_DATA(src))[idx+1]-update;/*subtract update*/
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					((int32_t*)HDR_DATA(src))[idx+1]=luma;
				}
			}
#endif
			for(kc0=0;kc0<src->nch;++kc0)
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

#if 1
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
#endif
				if(kc0>0)
				{
					offset+=rows[0][1];
					if(kc0>1)
						offset=(2*offset+rows[0][2])>>1;/*(2*g+[b-g])/2 = (g+b)/2*/
				}

				vmin=MINVAR(N, W), vmax=MAXVAR(N, W);
				UPDATE_MIN(vmin, NE);
				UPDATE_MAX(vmax, NE);
				{
					int64_t llpred, wsum;
					int preds[NPREDS2], j;

#define WPU_PRED(E) preds[j++]=E;
					j=0;
					WPU_PREDLIST;
#undef  WPU_PRED
					llpred=0;
					wsum=1;
					for(k=0;k<WPU_NPREDS;++k)
					{
						int64_t weight=0;
						weight+=(int64_t)rows[1][kc+0*WPU_PIXEL_STRIDE+8+((size_t)k<<2)];
#ifdef WPU_TRACE
						if(kc==1)
							trace[k]+=weight*weight;
#endif
						weight+=(int64_t)rows[1][kc-1*WPU_PIXEL_STRIDE+8+((size_t)k<<2)];
						weight+=(int64_t)rows[1][kc+1*WPU_PIXEL_STRIDE+8+((size_t)k<<2)];
						weight+=!weight;
						weight=((int64_t)weights[k<<2|kc]<<16)/weight;
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

					val=((int32_t*)HDR_DATA(src))[idx+kc];
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
					((int32_t*)HDR_DATA(src))[idx+kc]=val;
					curr<<=WPU_NBITS;
					rows[0][kc+0]=curr-offset;
					rows[0][kc+4]=curr-pred;
					{
						int kbest=0;
						for(k=0;k<WPU_NPREDS;++k)
						{
							int e2=abs(curr-preds[k]);
							eprev[k]=e2;
							rows[0][kc+0*WPU_PIXEL_STRIDE+8+((size_t)k<<2)]=e2;
							rows[1][kc+1*WPU_PIXEL_STRIDE+8+((size_t)k<<2)]+=e2;
							if(eprev[kbest]>eprev[k])
								kbest=k;
						}
						(void)kbest;
						{
							int rescale=0;
							for(k=0;k<WPU_NPREDS;++k)
							{
								weights[k<<2|kc]-=eprev[k];
								rescale|=weights[k<<2|kc]<0x800;
							}
							if(rescale)
							{
								for(k=0;k<WPU_NPREDS;++k)
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
					cb	=((int32_t*)HDR_DATA(src))[idx+2],
					cr	=((int32_t*)HDR_DATA(src))[idx+0];
				int update=(cb+cr)>>3;
				update=CLAMP(-halfs[1], update, halfs[1]-1);
				{
					int luma=((int32_t*)HDR_DATA(src))[idx+1]+update;/*add update*/
					luma<<=32-src->depth[1];
					luma>>=32-src->depth[1];
					((int32_t*)HDR_DATA(src))[idx+1]=luma;
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
		for(k=0;k<WPU_NPREDS;++k)
			printed+=snprintf(buf+printed, (intptr_t)1024-printed-1, "%20lld  %s\n", trace[k], prednames[k]);
		copy_to_clipboard(buf, printed);
		messagebox(MBOX_OK, "Copied to clipboard", buf);
	}
#endif
	free(pixels);
}

void pred_sub(Image *src, int fwd)/*W*/
{
	intptr_t k, size;
	int p0=0, p1=0, p2=0, p3=0;

	size=(intptr_t)4*src->iw*src->ih;
	if(g_dist>1)
	{
		int amin[4], amax[4], invdist;

		amin[0]=-(1<<src->depth[0]>>1);
		amin[1]=-(1<<src->depth[1]>>1);
		amin[2]=-(1<<src->depth[2]>>1);
		amin[3]=-(1<<src->depth[3]>>1);
		amax[0]=(1<<src->depth[0]>>1)-1;
		amax[1]=(1<<src->depth[1]>>1)-1;
		amax[2]=(1<<src->depth[2]>>1)-1;
		amax[3]=(1<<src->depth[3]>>1)-1;
		invdist=((1<<16)+g_dist-1)/g_dist;
		for(k=0;k<size;k+=4)
		{
			int c0=((int32_t*)HDR_DATA(src))[k|0];
			int c1=((int32_t*)HDR_DATA(src))[k|1];
			int c2=((int32_t*)HDR_DATA(src))[k|2];
			int c3=((int32_t*)HDR_DATA(src))[k|3];
			int e0, e1, e2, e3;
			if(fwd)
			{
				e0=c0-p0;
				e1=c1-p1;
				e2=c2-p2;
				e3=c3-p3;
				e0=(e0*invdist>>16)-(e0>>31);/*curr/=g_dist*/
				e1=(e1*invdist>>16)-(e1>>31);
				e2=(e2*invdist>>16)-(e2>>31);
				e3=(e3*invdist>>16)-(e3>>31);
				if(src->nch>=0)((int32_t*)HDR_DATA(src))[k|0]=e0;
				if(src->nch>=1)((int32_t*)HDR_DATA(src))[k|1]=e1;
				if(src->nch>=2)((int32_t*)HDR_DATA(src))[k|2]=e2;
				if(src->nch>=3)((int32_t*)HDR_DATA(src))[k|3]=e3;
			}
			else
			{
				e0=c0;
				e1=c1;
				e2=c2;
				e3=c3;
			}
			c0=g_dist*e0+p0;
			c1=g_dist*e1+p1;
			c2=g_dist*e2+p2;
			c3=g_dist*e3+p3;
			CLAMP2(c0, amin[0], amax[0]);
			CLAMP2(c1, amin[1], amax[1]);
			CLAMP2(c2, amin[2], amax[2]);
			CLAMP2(c3, amin[3], amax[3]);
			if(!fwd)
			{
				if(src->nch>=0)((int32_t*)HDR_DATA(src))[k|0]=c0;
				if(src->nch>=1)((int32_t*)HDR_DATA(src))[k|1]=c1;
				if(src->nch>=2)((int32_t*)HDR_DATA(src))[k|2]=c2;
				if(src->nch>=3)((int32_t*)HDR_DATA(src))[k|3]=c3;
			}
			p0=c0;
			p1=c1;
			p2=c2;
			p3=c3;
		}
	}
	else if(fwd)
	{
		for(k=0;k<size;k+=4)
		{
			int c0=((int32_t*)HDR_DATA(src))[k|0];
			int c1=((int32_t*)HDR_DATA(src))[k|1];
			int c2=((int32_t*)HDR_DATA(src))[k|2];
			int c3=((int32_t*)HDR_DATA(src))[k|3];
			int d0=c0-p0;
			int d1=c1-p1;
			int d2=c2-p2;
			int d3=c3-p3;
			d0<<=32-src->depth[0];
			d1<<=32-src->depth[1];
			d2<<=32-src->depth[2];
			d3<<=32-src->depth[3];
			d0>>=32-src->depth[0];
			d1>>=32-src->depth[1];
			d2>>=32-src->depth[2];
			d3>>=32-src->depth[3];
			p0=c0;
			p1=c1;
			p2=c2;
			p3=c3;
			if(src->nch>=0)((int32_t*)HDR_DATA(src))[k|0]=d0;
			if(src->nch>=1)((int32_t*)HDR_DATA(src))[k|1]=d1;
			if(src->nch>=2)((int32_t*)HDR_DATA(src))[k|2]=d2;
			if(src->nch>=3)((int32_t*)HDR_DATA(src))[k|3]=d3;
		}
	}
	else
	{
		for(k=0;k<size;k+=4)
		{
			int c0=((int32_t*)HDR_DATA(src))[k|0];
			int c1=((int32_t*)HDR_DATA(src))[k|1];
			int c2=((int32_t*)HDR_DATA(src))[k|2];
			int c3=((int32_t*)HDR_DATA(src))[k|3];
			int d0=c0+p0;
			int d1=c1+p1;
			int d2=c2+p2;
			int d3=c3+p3;
			d0<<=32-src->depth[0];
			d1<<=32-src->depth[1];
			d2<<=32-src->depth[2];
			d3<<=32-src->depth[3];
			d0>>=32-src->depth[0];
			d1>>=32-src->depth[1];
			d2>>=32-src->depth[2];
			d3>>=32-src->depth[3];
			p0=d0;
			p1=d1;
			p2=d2;
			p3=d3;
			if(src->nch>=0)((int32_t*)HDR_DATA(src))[k|0]=d0;
			if(src->nch>=1)((int32_t*)HDR_DATA(src))[k|1]=d1;
			if(src->nch>=2)((int32_t*)HDR_DATA(src))[k|2]=d2;
			if(src->nch>=3)((int32_t*)HDR_DATA(src))[k|3]=d3;
		}
	}
}
void pred_av2(Image *src, int fwd)/*(N+W)>>1*/
{
	int amin[4], amax[4], psize, ky, kx, kc, idx;
	int16_t *pixels, *rows[2];

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	psize=(src->iw+2)*(int)sizeof(int16_t[2*4]);
	pixels=(int16_t*)malloc(psize);/*2 padded rows * 4 channels max*/
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+(((src->iw+(intptr_t)2)*((ky-(intptr_t)0)&1)+1)<<2);
		rows[1]=pixels+(((src->iw+(intptr_t)2)*((ky-(intptr_t)1)&1)+1)<<2);
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<src->nch;++kc)
			{
				int N, W, pred, curr, val;

				N	=rows[1][kc+0];
				W	=rows[0][kc-4];
				pred=(N+W)>>1;
				curr=((int32_t*)HDR_DATA(src))[idx+kc];
				if(g_dist>1)
				{
					if(fwd)
					{
						val=(curr-(int)pred+g_dist/2)/g_dist;
						curr=g_dist*val+(int)pred;
					}
					else
					{
						val=g_dist*curr+(int)pred;
						curr=val;
						CLAMP2(val, amin[kc], amax[kc]);
					}
				}
				else if(fwd)
				{
					val=curr-pred;
					val<<=32-src->depth[kc];
					val>>=32-src->depth[kc];
				}
				else
				{
					val=curr+pred;
					val<<=32-src->depth[kc];
					val>>=32-src->depth[kc];
					curr=val;
				}
				((int32_t*)HDR_DATA(src))[idx+kc]=val;
				rows[0][kc]=curr;
			}
			rows[0]+=4;
			rows[1]+=4;
		}
	}
	free(pixels);
}
void pred_mix2(Image *src, int fwd)
{
	int nch, fwdmask, bufsize, *pixels, *rows[4], ky, kx, kc, idx;

	fwdmask=-fwd;
	bufsize=(src->iw+16)*(int)sizeof(int[4*4*2]);/*4 padded rows * 4 channels max*/
	pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	nch=(src->depth[0]!=0)+(src->depth[1]!=0)+(src->depth[2]!=0)+(src->depth[3]!=0);
	UPDATE_MAX(nch, src->nch);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int asum=0;

		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4*2;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4*2;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4*2;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4*2;
		for(kx=0;kx<src->iw;++kx, idx+=4)
		{
			for(kc=0;kc<3;++kc)
			{
				int
					NW	=rows[1][kc-1*4*2+0],
					N	=rows[1][kc+0*4*2+0],
					NE	=rows[1][kc+1*4*2+0],
					W	=rows[0][kc-1*4*2+0],
					aNNWW	=rows[2][kc-2*4*2+4],
					aNNW	=rows[2][kc-1*4*2+4],
					aNN	=rows[2][kc+0*4*2+4],
					aNNE	=rows[2][kc+1*4*2+4],
					aNNEE	=rows[2][kc+2*4*2+4],
					aNWW	=rows[1][kc-2*4*2+4],
					aNW	=rows[1][kc-1*4*2+4],
					aN	=rows[1][kc+0*4*2+4],
					aNE	=rows[1][kc+1*4*2+4],
					aNEE	=rows[1][kc+2*4*2+4],
					aNEEE	=rows[1][kc+3*4*2+4],
					aNEEEE	=rows[1][kc+4*4*2+4],
					aWWWW	=rows[0][kc-4*4*2+4],
					aWWW	=rows[0][kc-3*4*2+4],
					aWW	=rows[0][kc-2*4*2+4],
					aW	=rows[0][kc-1*4*2+4];
				int pred, curr;

#if 1
				(void)NW;
				(void)N;
				(void)NE;
				(void)W;
				(void)aNNWW;
				(void)aNNW;
				(void)aNN;
				(void)aNNE;
				(void)aNNEE;
				(void)aNWW;
				(void)aNW;
				(void)aN;
				(void)aNE;
				(void)aNEE;
				(void)aNEEE;
				(void)aNEEEE;
				(void)aWWWW;
				(void)aWWW;
				(void)aWW;
				(void)aW;
#endif
				/*predict*/
				if(kc)/*chroma*/
					pred=W+((N-W)*(asum+(1<<8>>1))>>8)+(NE-NW)/4;
				else/*luma*/
				{
					MEDIAN3_32(pred, N, W, N+W-NW);
				}

				curr=((int32_t*)HDR_DATA(src))[idx+kc];
				pred^=fwdmask;
				pred-=fwdmask;
				pred+=curr;

				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];

				((int32_t*)HDR_DATA(src))[idx+kc]=pred;
				rows[0][kc]=fwd?curr:pred;

				/*update*/
				curr=rows[0][kc];
				if(N!=W)
				{
					int alpha=((curr-W)<<8|1<<8>>1)/(N-W);
					CLAMP2(alpha, 0, 1<<8);
					alpha-=1<<8>>1;

					asum+=(alpha-asum)>>6;
				}
			}
			rows[0]+=4*2;
			rows[1]+=4*2;
			rows[2]+=4*2;
			rows[3]+=4*2;
		}
	}
	free(pixels);
}


	#define WP_RCT
	#define ETOTAL_MBOX

#define NWP (intptr_t)8		/*multiple of 4*/
#define LGBLOCKSIZE 0
#define BLOCKSIZE (1<<LGBLOCKSIZE)
#define WP_PBITS 5
#define WP_WBITS 10
void pred_wp_deferred(Image *src, int fwd)
{
	__m128i pixelmask, pixelhalf, roundoffset, fwdmask, invmask;
	ALIGN(16) int pred_errors[4*NWP]={0};/*interleaved channels*/
	int nlevels[4], nblocks, wsize, *weights, psize, *pixels, *rows[4], ky, kx, kc, idx;
#ifdef ETOTAL_MBOX
	int64_t etotal[4*NWP]={0};
#endif

	nlevels[0]=1<<src->depth[0];
	nlevels[1]=1<<src->depth[1];
	nlevels[2]=1<<src->depth[2];
	nlevels[3]=1<<src->depth[3];
	pixelmask=_mm_set_epi32(0, nlevels[2]-1, nlevels[1]-1, nlevels[0]-1);
	pixelhalf=_mm_set_epi32(0, nlevels[2]>>1, nlevels[1]>>1, nlevels[0]>>1);
	roundoffset=_mm_set1_epi32(1<<WP_PBITS>>1);
	fwdmask=_mm_set1_epi32(-fwd);
	invmask=_mm_set1_epi32(-!fwd);
	nblocks=(src->iw+BLOCKSIZE-1)>>LGBLOCKSIZE;
	wsize=(nblocks+NWP)*(int)sizeof(int[4*NWP]);
	weights=(int*)_mm_malloc(wsize, sizeof(__m256i));/*interleaved channels*/
	psize=(src->iw+16)*(int)sizeof(int[4*4*2]);
	pixels=(int*)malloc(psize);/*4 padded rows * 4 channels max * {pixel, error}		int, to accelerate _mullo_epi32*/
	if(!weights||!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	weights[0]=(int)(((intptr_t)1<<WP_WBITS)/NWP);
	memfill(weights+1, weights, wsize-sizeof(*weights), sizeof(*weights));
	memset(pixels, 0, psize);
#ifdef WP_RCT
	src->depth[0]+=(char)fwd;
	src->depth[2]+=(char)fwd;
#endif
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int *wp=weights;

		rows[0]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+4)<<3);
		rows[1]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+4)<<3);
		rows[2]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+4)<<3);
		rows[3]=pixels+(((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+4)<<3);
		for(kx=0;kx<src->iw;++kx, idx+=4)
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

				{
					ALIGN(16) int curr[4];
					__m128i w0=_mm_load_si128((__m128i*)wp+0);
					__m128i w1=_mm_load_si128((__m128i*)wp+1);
					__m128i w2=_mm_load_si128((__m128i*)wp+2);
					__m128i w3=_mm_load_si128((__m128i*)wp+3);
					__m128i w4=_mm_load_si128((__m128i*)wp+4);
					__m128i w5=_mm_load_si128((__m128i*)wp+5);
					__m128i w6=_mm_load_si128((__m128i*)wp+6);
					__m128i w7=_mm_load_si128((__m128i*)wp+7);

					w0=_mm_mullo_epi32(w0, wp0);/*FIXME use AVX2*/
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
					curr[0]=((int32_t*)HDR_DATA(src))[idx+0];
					curr[1]=((int32_t*)HDR_DATA(src))[idx+1];
					curr[2]=((int32_t*)HDR_DATA(src))[idx+2];
					curr[3]=0;
#ifdef WP_RCT
					if(fwd)
					{
						curr[0]-=curr[1];
						curr[2]-=curr[1];
						curr[1]+=(curr[0]+curr[2])>>2;
					}
#endif
					{
						__m128i mc=_mm_load_si128((__m128i*)curr);
						__m128i pixel=_mm_and_si128(mc, fwdmask);
						__m128i crisppred=_mm_add_epi32(w0, roundoffset);
						crisppred=_mm_srai_epi32(crisppred, WP_PBITS);

						crisppred=_mm_xor_si128(crisppred, fwdmask);
						crisppred=_mm_sub_epi32(crisppred, fwdmask);
						crisppred=_mm_add_epi32(crisppred, mc);/*error = curr - pred*/

						crisppred=_mm_add_epi32(crisppred, pixelhalf);
						crisppred=_mm_and_si128(crisppred, pixelmask);
						crisppred=_mm_sub_epi32(crisppred, pixelhalf);
						_mm_store_si128((__m128i*)curr, crisppred);
#ifdef WP_RCT
						if(!fwd)
						{
							curr[1]-=(curr[0]+curr[2])>>2;
							curr[2]+=curr[1];
							curr[0]+=curr[1];
						}
#endif
						((int32_t*)HDR_DATA(src))[idx+0]=curr[0];
						((int32_t*)HDR_DATA(src))[idx+1]=curr[1];
						((int32_t*)HDR_DATA(src))[idx+2]=curr[2];
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
			
			if(kx&&!(kx&(BLOCKSIZE-1)))/*update WP*/
			{
				for(kc=0;kc<3;++kc)
				{
					int64_t
						weight0=((int64_t)1<<32)/((int64_t)pred_errors[kc+0*4]+1),
						weight1=((int64_t)1<<32)/((int64_t)pred_errors[kc+1*4]+1),
						weight2=((int64_t)1<<32)/((int64_t)pred_errors[kc+2*4]+1),
						weight3=((int64_t)1<<32)/((int64_t)pred_errors[kc+3*4]+1),
						weight4=((int64_t)1<<32)/((int64_t)pred_errors[kc+4*4]+1),
						weight5=((int64_t)1<<32)/((int64_t)pred_errors[kc+5*4]+1),
						weight6=((int64_t)1<<32)/((int64_t)pred_errors[kc+6*4]+1),
						weight7=((int64_t)1<<32)/((int64_t)pred_errors[kc+7*4]+1),
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
					etotal[kc+0*4]+=(int64_t)pred_errors[kc+0*4];
					etotal[kc+1*4]+=(int64_t)pred_errors[kc+1*4];
					etotal[kc+2*4]+=(int64_t)pred_errors[kc+2*4];
					etotal[kc+3*4]+=(int64_t)pred_errors[kc+3*4];
					etotal[kc+4*4]+=(int64_t)pred_errors[kc+4*4];
					etotal[kc+5*4]+=(int64_t)pred_errors[kc+5*4];
					etotal[kc+6*4]+=(int64_t)pred_errors[kc+6*4];
					etotal[kc+7*4]+=(int64_t)pred_errors[kc+7*4];
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
		int64_t esum[3]={0};
		int printed, kp;

		printed=0;
		for(kp=0;kp<NWP;++kp)
		{
			esum[0]+=etotal[kp<<2|0];
			esum[1]+=etotal[kp<<2|1];
			esum[2]+=etotal[kp<<2|2];
		}
		for(kp=0;kp<NWP;++kp)
		{
			for(kc=0;kc<3;++kc)
				printed+=snprintf(buf+printed, BUFSIZE-printed-1, " %6.2f", 100.*etotal[kp<<2|kc]/esum[kc]);
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

void pred_av4(Image *src, int fwd, int enable_ma)
{
	enum
	{
		XPAD=8,
		NROWS=2,
		NCH=4,
		NVAL=1
	};
	int amin[4], amax[4], psize, ky, kx, kc, idx;
	int16_t *pixels, *rows[2];

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	psize=(src->iw+2*XPAD)*(int)sizeof(int16_t[NROWS*NCH*NVAL]);
	pixels=(int16_t*)_mm_malloc(psize, sizeof(__m128i));
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, psize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)0+NROWS)%NROWS)*NVAL;/*sub 1 channel for pre-increment*/
		rows[1]=pixels+(XPAD*NCH*NROWS-NROWS+(ky-(intptr_t)1+NROWS)%NROWS)*NVAL;
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<4;++kc, ++idx)
			{
				rows[0]+=NROWS*NVAL;
				rows[1]+=NROWS*NVAL;
				if(!src->depth[kc])
					continue;
				{
					int16_t
						NW	=rows[1][0-1*NCH*NROWS*NVAL],
						N	=rows[1][0+0*NCH*NROWS*NVAL],
						NE	=rows[1][0+1*NCH*NROWS*NVAL],
						NEEE	=rows[1][0+3*NCH*NROWS*NVAL],
						W	=rows[0][0-1*NCH*NROWS*NVAL];
					int pred, vmin, vmax, curr, val;

					pred=(4*(N+W)+NE-NW)>>3;
					vmax=N; vmin=W;
					if(N<W)vmin=N, vmax=W;
					if(vmin>NE)vmin=NE;
					if(vmax<NE)vmax=NE;
					if(vmin>NEEE)vmin=NEEE;
					if(vmax<NEEE)vmax=NEEE;
					CLAMP2(pred, vmin, vmax);

					curr=((int32_t*)HDR_DATA(src))[idx];
					if(g_dist>1)
					{
						if(fwd)
						{
							val=(curr-(int)pred+g_dist/2)/g_dist;
							curr=g_dist*val+(int)pred;
						}
						else
						{
							val=g_dist*curr+(int)pred;
							curr=val;
							CLAMP2(val, amin[kc], amax[kc]);
						}
					}
					else if(fwd)
					{
						val=curr-pred;
						val<<=32-src->depth[kc];
						val>>=32-src->depth[kc];
					}
					else
					{
						val=curr+pred;
						val<<=32-src->depth[kc];
						val>>=32-src->depth[kc];
						curr=val;
					}
					((int32_t*)HDR_DATA(src))[idx]=val;
					rows[0][0]=curr;
				}
			}
		}
	}
	_mm_free(pixels);

	(void)enable_ma;
}

void pred_ecoeff(Image *src, int fwd, int enable_ma)
{
	Image *dst;
	const int *pixels, *errors;
	int k, ky, kx, kc, idx;

	dst=0;
	image_copy(&dst, src);
	if(!dst)
	{
		CRASH("Alloc error");
		return;
	}
	pixels=fwd?((int32_t*)HDR_DATA(dst)):((int32_t*)HDR_DATA(src));
	errors=fwd?((int32_t*)HDR_DATA(src)):((int32_t*)HDR_DATA(dst));
	for(kc=0;kc<4;++kc)
	{
		int nlevels, half;

		if(!src->depth[kc])
			continue;
		nlevels=1<<src->depth[kc];
		half=nlevels>>1;
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			int ecoeffs[12]={0};
			for(kx=0;kx<src->iw;++kx, ++idx)
			{
				int NW, N, W, enb[12], j, pred, correction, e;
#define LOAD(BUF, X, Y) ((uint32_t)(ky+(Y))<(uint32_t)src->ih&&(uint32_t)(kx+(X))<(uint32_t)src->iw?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
				NW	=LOAD(pixels, -1, -1);
				N	=LOAD(pixels,  0, -1);
				W	=LOAD(pixels, -1,  0);
				j=0;
				enb[j++]=LOAD(errors, -2, -2);
				enb[j++]=LOAD(errors, -1, -2);
				enb[j++]=LOAD(errors,  0, -2);
				enb[j++]=LOAD(errors,  1, -2);
				enb[j++]=LOAD(errors,  2, -2);
				enb[j++]=LOAD(errors, -2, -1);
				enb[j++]=LOAD(errors, -1, -1);
				enb[j++]=LOAD(errors,  0, -1);
				enb[j++]=LOAD(errors,  1, -1);
				enb[j++]=LOAD(errors,  2, -1);
				enb[j++]=LOAD(errors, -2,  0);
				enb[j++]=LOAD(errors, -1,  0);
#undef  LOAD
				correction=0;
				for(k=0;k<(int)_countof(enb);++k)
					correction+=enb[k]*ecoeffs[k];
				correction+=12*256/2-1;
				correction/=12*256;

				pred=N+W-NW;
				pred=MEDIAN3(N, W, pred);

				pred+=correction;
				pred=CLAMP(-half, pred, half-1);
				
				pred^=-fwd;
				pred+=fwd;
				pred+=((int32_t*)HDR_DATA(dst))[idx<<2|kc];
				if(enable_ma)
				{
					pred+=nlevels>>1;
					pred&=nlevels-1;
					pred-=nlevels>>1;
				}
				((int32_t*)HDR_DATA(src))[idx<<2|kc]=pred;

				/*
				L = sq(curr - pred)
				dL/dcoeff = (curr - pred) * -dpred/dcoeff
					= (pred - curr)*srcval
				coeff -= dL/dcoeff * lr
				coeff += (curr-pred)*srcval*lr
				*/
				e=errors[idx<<2|kc];
				for(k=0;k<(int)_countof(enb);++k)
				{
					int c=ecoeffs[k];
					int update=e*enb[k];
					c+=(update>>((src->depth[kc]<<1)-12))+(update<0);
					ecoeffs[k]=CLAMP(-256, c, 256);
				}
			}
		}
	}
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}
	free(dst);
}
void pred_zipper(Image *src, int fwd, int enable_ma)
{
	int half[4], negmask, step, ky, kx, kc;

	negmask=-fwd;
	half[0]=1<<src->depth[0]>>1;
	half[1]=1<<src->depth[1]>>1;
	half[2]=1<<src->depth[2]>>1;
	half[3]=1<<src->depth[3]>>1;
	for(ky=fwd?src->ih-1:0, step=fwd?-1:1;(uint32_t)ky<(uint32_t)src->ih;ky+=step)
	{
		for(kx=fwd?src->iw-1:0;(uint32_t)kx<(uint32_t)src->iw;kx+=step)
		{
			int idx=src->iw*ky+kx;
			for(kc=0;kc<src->nch;++kc)
			{
#define LOAD(X, Y) ((uint32_t)(ky+(Y))<(uint32_t)src->ih&&(uint32_t)(kx+(X))<(uint32_t)src->iw?((int32_t*)HDR_DATA(src))[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
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
#undef  LOAD
				int pred;
				int
					vx=abs(eE-eEE)+abs(eS-eSE)+abs(eS-eSW),
					vy=abs(eSSW-eSW)+abs(eSS-eS)+abs(eSE-eE);
				int update=(vy>vx?eS:eE)/4;
				MEDIAN3_32(pred, N, W, N+W-NW);
				pred+=update;
				CLAMP2(pred, -half[kc], half[kc]-1);
				
				pred^=negmask;
				pred-=negmask;
				pred+=((int32_t*)HDR_DATA(src))[idx<<2|kc];
				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];
				((int32_t*)HDR_DATA(src))[idx<<2|kc]=pred;

				(void)NW;
				(void)N;
				(void)W;
				(void)eE;
				(void)eEE;
				(void)eSWW;
				(void)eSW;
				(void)eS;
				(void)eSE;
				(void)eSEE;
				(void)eSSWW;
				(void)eSSW;
				(void)eSS;
				(void)eSSE;
				(void)eSSEE;
			}
		}
	}
	(void)enable_ma;
}

#define MULTISTAGE_CALCCSIZE
#define LOAD(BUF, X, Y) ((uint32_t)(ky+(Y))<(uint32_t)src->ih&&(uint32_t)(kx+(X))<(uint32_t)src->iw?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
static void pred_multistage_stage1(Image *src, Image *dst, int fwd, int enable_ma, int kc, int *hist)
{
	const int *pixels;
	int nlevels, ky, kx;

	nlevels=1<<src->depth[kc];
	pixels=fwd?((int32_t*)HDR_DATA(src)):((int32_t*)HDR_DATA(dst));
#ifdef MULTISTAGE_CALCCSIZE
	memset(hist, 0, (nlevels+(intptr_t)1)*sizeof(int));
#endif
	for(ky=0;ky<src->ih;ky+=2)
	{
		for(kx=0;kx<src->iw;kx+=2)
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
			pred+=((int32_t*)HDR_DATA(src))[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
#ifdef MULTISTAGE_CALCCSIZE
				++hist[pred+(nlevels>>1)+1];
				++hist[0];
#endif
			}
			((int32_t*)HDR_DATA(dst))[idx<<2|kc]=pred;
		}
	}
}
static void pred_multistage_stage2c(Image *src, Image *dst, int fwd, int enable_ma, int kc, int *hist)
{
	const int *pixels;
	int nlevels, half, ky, kx;

	nlevels=1<<src->depth[kc];
	half=nlevels>>1;
	pixels=fwd?((int32_t*)HDR_DATA(src)):((int32_t*)HDR_DATA(dst));
#ifdef MULTISTAGE_CALCCSIZE
	memset(hist, 0, (nlevels+(intptr_t)1)*sizeof(int));
#endif
	for(ky=1;ky<src->ih;ky+=2)
	{
		for(kx=1;kx<src->iw;kx+=2)
		{
			int
				idx=src->iw*ky+kx,
				NW	=LOAD(pixels, -1, -1),
				SW	=LOAD(pixels, -1,  1),
				NE	=LOAD(pixels,  1, -1),
				SE	=LOAD(pixels,  1,  1),
				NNNWWW	=LOAD(pixels, -3, -3),
				SSSWWW	=LOAD(pixels, -3,  3),
				NNNEEE	=LOAD(pixels,  3, -3),
				SSSEEE	=LOAD(pixels,  3,  3);

			int pred=(9*(NW+NE+SW+SE)-(NNNWWW+NNNEEE+SSSWWW+SSSEEE)+16)>>5;

			pred=CLAMP(-half, pred, half-1);
				
			pred^=-fwd;
			pred+=fwd;
			pred+=((int32_t*)HDR_DATA(src))[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
#ifdef MULTISTAGE_CALCCSIZE
				++hist[pred+(nlevels>>1)+1];
				++hist[0];
#endif
			}
			((int32_t*)HDR_DATA(dst))[idx<<2|kc]=pred;
		}
	}
}
static void pred_multistage_stage3c(Image *src, Image *dst, int fwd, int enable_ma, int kc, int *hist)
{
	const int *pixels;
	int nlevels, half, ky, kx;

	nlevels=1<<src->depth[kc];
	half=nlevels>>1;
	pixels=fwd?((int32_t*)HDR_DATA(src)):((int32_t*)HDR_DATA(dst));
#ifdef MULTISTAGE_CALCCSIZE
	memset(hist, 0, (nlevels+(intptr_t)1)*sizeof(int));
#endif
	for(ky=1;ky<src->ih;ky+=2)
	{
		for(kx=0;kx<src->iw;kx+=2)
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
			pred=CLAMP(-half, pred, half-1);
				
			pred^=-fwd;
			pred+=fwd;
			pred+=((int32_t*)HDR_DATA(src))[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
#ifdef MULTISTAGE_CALCCSIZE
				++hist[pred+(nlevels>>1)+1];
				++hist[0];
#endif
			}
			((int32_t*)HDR_DATA(dst))[idx<<2|kc]=pred;
		}
	}
}
static void pred_multistage_stage4c(Image *src, Image *dst, int fwd, int enable_ma, int kc, int *hist)
{
	const int *pixels;
	int nlevels, half, ky, kx;

	nlevels=1<<src->depth[kc];
	half=nlevels>>1;
	pixels=fwd?((int32_t*)HDR_DATA(src)):((int32_t*)HDR_DATA(dst));
#ifdef MULTISTAGE_CALCCSIZE
	memset(hist, 0, (nlevels+(intptr_t)1)*sizeof(int));
#endif
	for(ky=0;ky<src->ih;ky+=2)
	{
		for(kx=1;kx<src->iw;kx+=2)
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
			pred=CLAMP(-half, pred, half-1);
				
			pred^=-fwd;
			pred+=fwd;
			pred+=((int32_t*)HDR_DATA(src))[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
#ifdef MULTISTAGE_CALCCSIZE
				++hist[pred+(nlevels>>1)+1];
				++hist[0];
#endif
			}
			((int32_t*)HDR_DATA(dst))[idx<<2|kc]=pred;
		}
	}
}
void pred_multistage(Image *src, int fwd, int enable_ma)
{
	Image *dst;
	int kc;
#ifdef MULTISTAGE_CALCCSIZE
	double csizes[16]={0};
	int maxdepth, maxlevels, hsize, *hist, k;
#endif

#ifdef MULTISTAGE_CALCCSIZE
	maxdepth=calc_maxdepth(src, 0);
	maxlevels=1<<maxdepth;
	hsize=(maxlevels+1)*(int)sizeof(int);
	hist=(int*)malloc(hsize);
#define CALC(X) csizes[kc<<2|X]=calc_entropy(hist+1, 1<<src->depth[kc], hist[0])
#else
#define CALC(X)
#endif
	dst=0;
	image_copy_nodata(&dst, src);
	if(!dst
#ifdef MULTISTAGE_CALCCSIZE
		||!hist
#endif
	)
	{
		CRASH("Alloc error");
		return;
	}
#ifdef MULTISTAGE_CALCCSIZE
	memset(((int32_t*)HDR_DATA(dst)), 0, (size_t)src->iw*src->ih*sizeof(int[4]));
#endif
	for(kc=0;kc<4;++kc)
	{
		if(!src->depth[kc])
			continue;
		pred_multistage_stage1(src, dst, fwd, enable_ma, kc, hist);	CALC(0);
		pred_multistage_stage2c(src, dst, fwd, enable_ma, kc, hist);	CALC(1);
		pred_multistage_stage3c(src, dst, fwd, enable_ma, kc, hist);	CALC(2);
		pred_multistage_stage4c(src, dst, fwd, enable_ma, kc, hist);	CALC(3);
	}
#ifdef MULTISTAGE_CALCCSIZE
	for(k=0;k<16;++k)
		csizes[k]*=100./8;
	if(fwd&&GET_KEY_STATE(KEY_SHIFT))
	{
		int printed=snprintf(g_buf, G_BUF_SIZE,
			"Y0123 %7.4f  %7.4f %7.4f %7.4f %7.4f\n"
			"U0123 %7.4f  %7.4f %7.4f %7.4f %7.4f\n"
			"V0123 %7.4f  %7.4f %7.4f %7.4f %7.4f\n",
			(csizes[ 0]+csizes[ 1]+csizes[ 2]+csizes[ 3])/4, csizes[ 0], csizes[ 1], csizes[ 2], csizes[ 3],
			(csizes[ 4]+csizes[ 5]+csizes[ 6]+csizes[ 7])/4, csizes[ 4], csizes[ 5], csizes[ 6], csizes[ 7],
			(csizes[ 8]+csizes[ 9]+csizes[10]+csizes[11])/4, csizes[ 8], csizes[ 9], csizes[10], csizes[11]
		);
		copy_to_clipboard(g_buf, printed);
		messagebox(MBOX_OK, "Copied to clipboard", g_buf);
	}
	free(hist);
#endif
#undef  CALC
	memcpy(((int32_t*)HDR_DATA(src)), ((int32_t*)HDR_DATA(dst)), (size_t)src->iw*src->ih*sizeof(int[4]));
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

enum
{
	P3_MAXPREDS=256
};
#define LOAD(BUF, X, Y) ((uint32_t)(ky+(Y))<(uint32_t)src->ih&&(uint32_t)(kx+(X))<(uint32_t)src->iw?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
static void pred_separate_x(Image const *src, Image *dst, int fwd, int enable_ma, int kc
#ifdef P3_CALCCSIZE
	, int *hist, double *csizes
#endif
)
{
	const int *pixels;
	int nlevels, ky, kx, idx;

	nlevels=1<<src->depth[kc];
	pixels=fwd?((const int32_t*)HDR_DATC(src)):((const int32_t*)HDR_DATC(dst));
#ifdef P3_CALCCSIZE
	memset(hist, 0, nlevels*sizeof(int[P3_MAXPREDS]));
#endif
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int64_t errors[P3_MAXPREDS]={0};

		*errors=nlevels<<8;
		memfill(errors+1, errors, sizeof(errors)-sizeof(int64_t), sizeof(int64_t));
		for(kx=0;kx<src->iw;++kx, ++idx)
		{
			int
				NNNNN		=LOAD(pixels,  0, -5),
				NNNNWW		=LOAD(pixels, -2, -4),
				NNNNW		=LOAD(pixels, -1, -4),
				NNNN		=LOAD(pixels,  0, -4),
				NNNNE		=LOAD(pixels,  1, -4),
				NNNNEE		=LOAD(pixels,  2, -4),
				NNNNEEEE	=LOAD(pixels,  4, -4),
				NNNWWW		=LOAD(pixels, -3, -3),
				NNN		=LOAD(pixels,  0, -3),
				NNNEE		=LOAD(pixels,  2, -3),
				NNNEEE		=LOAD(pixels,  3, -3),
				NNNEEEE		=LOAD(pixels,  4, -3),
				NNWWWW		=LOAD(pixels, -4, -2),
				NNWWW		=LOAD(pixels, -3, -2),
				NNWW		=LOAD(pixels, -2, -2),
				NNW		=LOAD(pixels, -1, -2),
				NN		=LOAD(pixels,  0, -2),
				NNE		=LOAD(pixels,  1, -2),
				NNEE		=LOAD(pixels,  2, -2),
				NNEEE		=LOAD(pixels,  3, -2),
				NNEEEE		=LOAD(pixels,  4, -2),
				NWWWW		=LOAD(pixels, -4, -1),
				NWWW		=LOAD(pixels, -3, -1),
				NWW		=LOAD(pixels, -2, -1),
				NW		=LOAD(pixels, -1, -1),
				N		=LOAD(pixels,  0, -1),
				NE		=LOAD(pixels,  1, -1),
				NEE		=LOAD(pixels,  2, -1),
				NEEE		=LOAD(pixels,  3, -1),
				NEEEE		=LOAD(pixels,  4, -1),
				NEEEEE		=LOAD(pixels,  5, -1),
				NEEEEEE		=LOAD(pixels,  6, -1),
				NEEEEEEE	=LOAD(pixels,  7, -1),
				NEEEEEEEE	=LOAD(pixels,  8, -1),
				WWWWWWWWW	=LOAD(pixels, -9,  0),
				WWWWWWWW	=LOAD(pixels, -8,  0),
				WWWWWWW		=LOAD(pixels, -7,  0),
				WWWWWW		=LOAD(pixels, -6,  0),
				WWWWW		=LOAD(pixels, -5,  0),
				WWWW		=LOAD(pixels, -4,  0),
				WWW		=LOAD(pixels, -3,  0),
				WW		=LOAD(pixels, -2,  0),
				W		=LOAD(pixels, -1,  0);
			int64_t lpred, lsum;
			int preds[P3_MAXPREDS], j, k, pred, curr;

			j=0;/*P3*/
			preds[j++]= 0x2E00000; preds[j++]=N+W-NW;
			preds[j++]= 0x0400000; preds[j++]=N+W-((NW+NN+WW+NE)>>2);
			preds[j++]=-0x0300000; preds[j++]=((2*(N+W)-(NW+NN+WW+NE))*9+(WWW+NWW+NNW+NNN+NNE+NEE)*2)/12;
			preds[j++]= 0x0C00000; preds[j++]=2*(N+W-NW)-(NN+WW-NNWW);
			preds[j++]= 0x0200000; preds[j++]=3*(N+W-NW-(NN+WW-NNWW))+NNN+WWW-NNNWWW;
			preds[j++]=-0x0300000; preds[j++]=N+W-NW+NEE-(NE+NNEE-NNE);
			preds[j++]= 0x0100000; preds[j++]=N+W-NW+WW+NEE-(WWW+NWW-NWWW+NE+NNEE-NNE);
			preds[j++]=-0x0800000; preds[j++]=(N+W)>>1;
			preds[j++]=-0x0100000; preds[j++]=N+W-((NN+WW)>>1);
			preds[j++]= 0x0200000; preds[j++]=(4*(N+W)+NE-NW)>>2;/*shift right 3?*/
			preds[j++]= 0x0800000; preds[j++]=(4*(N+W)-NE+NW)>>2;/*shift right 3?*/
			preds[j++]= 0x0800000; preds[j++]=(4*(N+W)-(NE+NW))>>2;
			preds[j++]=-0x0200000; preds[j++]=(4*(N+W)-(NE+NW))/6;
			preds[j++]=-0x0100000; preds[j++]=(4*N-(NW+NE))>>2;
			preds[j++]=-0x0E00000; preds[j++]=(W+NE)>>1;
			preds[j++]=-0x0E00000; preds[j++]=(2*W+NEE)/3;
			preds[j++]= 0x0F00000; preds[j++]=(3*W+NEEE)>>2;
			preds[j++]= 0x1000000; preds[j++]=4*W+N+NE+NEE+NEEE;/*shift right 3?*/
			preds[j++]=-0x0500000; preds[j++]=(0x289*W+0x080*N+0x06A*NE+0x04F*NEE+0x03E*NEEE)>>10;/*proportional to distance*/
			preds[j++]=-0x0500000; preds[j++]=(192*W+33*N+22*NE+11*NEE+6*NEEE)/264;/*proportional to square distance*/
			preds[j++]= 0x0500000; preds[j++]=W+NE-N;
			preds[j++]= 0x0200000; preds[j++]=2*(W+NE-N)-(WW+NNEE-NN);
			preds[j++]= 0x0200000; preds[j++]=3*(W+NE-N-(WW+NNEE-NN))+WWW+NNNEEE-NNN;
			preds[j++]= 0x0500000; preds[j++]=W+NE-N+NW-(NWW+NN-NNW);
			preds[j++]=-0x0100000; preds[j++]=W+NE-((3*N+NW)>>2);
			preds[j++]= 0x0800000; preds[j++]=(2*W+NEE-N)>>1;
			preds[j++]= 0x0800000; preds[j++]=(3*W+NEEE-N)>>1;
			preds[j++]= 0x1400000; preds[j++]=(3*(3*W+NE+NEE)-10*N)/5;
			preds[j++]= 0x0300000; preds[j++]=(WWWWWWWWW+WWWWWWWW)>>1;
			preds[j++]= 0x0400000; preds[j++]=(WWWWWWW+WWWWWW)>>1;
			preds[j++]=-0x0100000; preds[j++]=(WWWWW+WWWW+NWWW+NEEE)>>2;
			preds[j++]= 0x0200000; preds[j++]=(NWWWW+NWWW+NEEEE+NEEEEE)>>2;
			preds[j++]= 0x0200000; preds[j++]=W+NW-NWW;
			preds[j++]= 0x0700000; preds[j++]=NW+NWW-NNWWW;
			preds[j++]= 0x0600000; preds[j++]=N+NW-NNW;
			preds[j++]= 0x1800000; preds[j++]=N+NE-NNE;
			preds[j++]=-0x0500000; preds[j++]=(9*(NW-NNWW)+NNNNWW+NNNWWW+NNWWWW)/3;
			preds[j++]= 0x0300000; preds[j++]=(9*(NE-NNEE)+NNNNEE+NNNEEE+NNEEEE)/3;
			preds[j++]= 0x0500000; preds[j++]=NN+NNW-(NNNN+NNNNW+NNNNWW)/3;
			preds[j++]= 0x0500000; preds[j++]=2*NN-(NNNNW+NNNN+NNNNE)/3;
			preds[j++]= 0x0500000; preds[j++]=NN+NNE-(NNNN+NNNNE+NNNNEE)/3;
			preds[j++]= 0x0B00000; preds[j++]=NNEE+NEE-(NNNNEEEE+NNNEEEE+NNEEEE)/3;
			preds[j++]= 0x0500000; preds[j++]=3*NE-(NNE+NEE);
			preds[j++]=-0x0600000; preds[j++]=(14*NE-(NNEE+NNNEE+NNEEE))/11;
			preds[j++]=-0x0500000; preds[j++]=((NW+N+W)*10-(NNWW+NNW+NN+NWW+WW)*3)/15;
			preds[j++]=-0x0300000; preds[j++]=(75*W+53*NW)>>7;
			preds[j++]=-0x0300000; preds[j++]=(75*N+53*NW)>>7;
			preds[j++]=-0x0300000; preds[j++]=(75*N+53*NE)>>7;
			preds[j++]= 0x0F00000; preds[j++]=(2*(N+NE+NEE)-(NNE+NNEE+NNEEE))/3;
			preds[j++]= 0x0F00000; preds[j++]=NE+NEE-(NNEE+NNEEE+NNEEEE)/3;
			preds[j++]= 0x1400000; preds[j++]=(NEEE+NEEEE)>>1;
			preds[j++]= 0x0E00000; preds[j++]=(NEEEEE+NEEEEEE)>>1;
			preds[j++]= 0x0C00000; preds[j++]=(NEEEEEEE+NEEEEEEEE)>>1;
			preds[j++]= 0x1800000; preds[j++]=(4*N-2*NN+NW+NE)>>2;
			preds[j++]= 0x0800000; preds[j++]=(N+W+2*(NW+NE)-(NNWW+NNEE))>>2;
			preds[j++]= 0x0A00000; preds[j++]=2*N-W;
			preds[j++]= 0x0A00000; preds[j++]=2*W-N;
			preds[j++]=-0x0500000; preds[j++]=(2*(NW+2*N+NE)-(NNW+2*NN+NNE))>>2;
			preds[j++]= 0x0500000; preds[j++]=4*(N-W)-NW;
			preds[j++]= 0x0100000; preds[j++]=N+W+NW+NE;
			preds[j++]=-0x0400000; preds[j++]=2*NE-NN;
			preds[j++]=-0x0400000; preds[j++]=2*NE-NNE;
			preds[j++]=-0x0600000; preds[j++]=(2*(W+NE)-(WW+NNEE))>>1;
			preds[j++]= 0x0F00000; preds[j++]=N;
			preds[j++]= 0x1400000; preds[j++]=2*N-NN;
			preds[j++]= 0x1000000; preds[j++]=3*(N-NN)+NNN;
			preds[j++]= 0x1800000; preds[j++]=4*(N+NNN)-6*NN-(NNNNW+NNNN+NNNNE)/3;
			preds[j++]= 0x1000000; preds[j++]=(N+NN)>>1;
			preds[j++]= 0x0800000; preds[j++]=(N+NN+NNN)>>1;		/*why 3/2?*/
			preds[j++]= 0x1200000; preds[j++]=(NN+NNN+NNNN+NNNNN)>>2;
			preds[j++]= 0x0F00000; preds[j++]=W;
			preds[j++]= 0x1400000; preds[j++]=2*W-WW;
			preds[j++]= 0x1000000; preds[j++]=3*(W-WW)+WWW;
			preds[j++]= 0x1800000; preds[j++]=4*(W+WWW)-6*WW-(WWWW+NWWWW)/2;
			preds[j++]= 0x1000000; preds[j++]=(W+WW)>>1;
			preds[j++]= 0x0800000; preds[j++]=(W+WW+WWW)>>1;		/*why 3/2?*/
			preds[j++]= 0x1200000; preds[j++]=(WW+WWW+WWWW+WWWWW)>>2;
			preds[j++]=-0x0800000; preds[j++]=NW;
			preds[j++]= 0x0400000; preds[j++]=2*NW-NNWW;
			preds[j++]= 0x0200000; preds[j++]=3*(NW-NNWW)+NNNWWW;
			preds[j++]= 0x0400000; preds[j++]=(NW+NNWW)>>1;
			preds[j++]= 0x0400000; preds[j++]=(NW+NNWW+NNNWWW)>>1;	/*why 3/2?*/
			preds[j++]=-0x0300000; preds[j++]=(NW+NNW+NWW+NNWW)>>1;
			preds[j++]=-0x0800000; preds[j++]=NE;
			preds[j++]= 0x0400000; preds[j++]=2*NE-NNEE;
			preds[j++]= 0x0200000; preds[j++]=3*(NE-NNEE)+NNNEEE;
			preds[j++]= 0x0400000; preds[j++]=(NE+NNEE)>>1;
			preds[j++]= 0x0400000; preds[j++]=(NE+NNEE+NNNEEE)>>1;	/*why 3/2?*/
			preds[j++]=-0x0300000; preds[j++]=(NE+NNE+NEE+NNEE)>>1;
			lpred=0, lsum=0;
			for(k=0;k<j/2;++k)
			{
				int64_t weight=((int64_t)preds[k<<1|0]<<8)/errors[k];
				lpred+=weight*preds[k<<1|1];
				lsum+=weight;
			}
			pred=(int)(lpred/lsum);
			CLAMP3_32(pred, pred, N, W, NE);
			pred^=-fwd;
			pred+=fwd;
			pred+=((const int32_t*)HDR_DATC(src))[idx<<2|kc];
			if(enable_ma)
			{
				pred+=nlevels>>1;
				pred&=nlevels-1;
				pred-=nlevels>>1;
			}
			((int32_t*)HDR_DATA(dst))[idx<<2|kc]=pred;
			curr=pixels[idx<<2|kc];
			for(k=0;k<j/2;++k)
			{
				int64_t e=errors[k];
				e+=(((int64_t)abs(curr-preds[k<<1|1])<<8)-e)/3;
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
	for(k=0;k<P3_MAXPREDS;++k)
		csizes[k<<2|kc]=calc_entropy(hist+((size_t)k<<src->depth[kc]), nlevels, src->iw*src->ih)*100./8;
#endif
}
#undef  LOAD
void pred_separate(Image *src, int fwd, int enable_ma)
{
	Image *dst;
	int kc;
	double t;
#ifdef P3_CALCCSIZE
	int maxdepth, *hist, k;
	double *csizes;
	char *str;
#endif

	t=time_sec();
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
		CRASH("Alloc error");
		return;
	}
	memset(((int32_t*)HDR_DATA(dst)), 0, (size_t)src->iw*src->ih*sizeof(int[4]));
	for(kc=0;kc<4;++kc)
	{
		if(!src->depth[kc])
			continue;
		pred_separate_x(src, dst, fwd, enable_ma, kc
#ifdef P3_CALCCSIZE
			, hist, csizes
#endif
		);
	}
	memcpy(((int32_t*)HDR_DATA(src)), ((int32_t*)HDR_DATA(dst)), (size_t)src->iw*src->ih*sizeof(int[4]));
#ifdef P3_CALCCSIZE
	if(GET_KEY_STATE(KEY_SHIFT)&&fwd)
	{
		int printed=0;
		for(k=0;k<P3_MAXPREDS;++k)
		{
			printed+=snprintf(str+printed, (intptr_t)0x100000-printed-1, "%d", k);
			for(kc=0;kc<4;++kc)
			{
				if(!src->depth[kc])
					continue;
				printed+=snprintf(str+printed, (intptr_t)0x100000-printed-1, "\t%8.4f", csizes[k<<2|kc]);
			}
			printed+=snprintf(str+printed, (intptr_t)0x100000-printed-1, "\n");
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
		if(currfn)
			set_window_title("%s - %f sec", (char*)HDR_DATA(currfn), t);
		else
			set_window_title("%f sec", t);
	}
}

void pred_dir(Image *src, int fwd, int enable_ma)
{
	Image *dst;
	const int *pixels;
	int maxdepth, ssesize, *ssearr, k, ky, kx, kc, idx;

	dst=0;
	image_copy(&dst, src);
	maxdepth=calc_maxdepth(src, 0);
	ssesize=(1<<maxdepth)<<1;
	ssearr=(int*)malloc(ssesize*sizeof(int));
	if(!ssearr||!dst)
	{
		CRASH("Alloc error");
		return;
	}
	pixels=fwd?((int32_t*)HDR_DATA(src)):((int32_t*)HDR_DATA(dst));
	for(kc=0;kc<4;++kc)
	{
		int nlevels;

		if(!src->depth[kc])
			continue;
		memset(ssearr, 0, ssesize*sizeof(int));
		nlevels=1<<src->depth[kc];
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx, ++idx)
			{
				int pred, curr, arg, nb[4];
#define LOAD(X, Y) ((uint32_t)(ky+(Y))<(uint32_t)src->ih&&(uint32_t)(kx+(X))<(uint32_t)src->iw?pixels[(idx+src->iw*(Y)+(X))<<2|kc]:0)
				nb[0]=LOAD(-1,  0);/*W	180	*/
				nb[1]=LOAD(-1, -1);/*NW	135	*/
				nb[2]=LOAD( 0, -1);/*N	90	*/
				nb[3]=LOAD( 1, -1);/*NE	45	*/
#undef  LOAD
				curr=pixels[idx<<2|kc];
				arg=0;
				for(k=1;k<(int)_countof(nb);++k)
				{
					if(abs(curr-nb[arg])>abs(curr-nb[k]))
						arg=k;
				}
				pred=nb[arg];
				
				pred^=-fwd;
				pred+=fwd;
				pred+=((int32_t*)HDR_DATA(src))[idx<<2|kc];
				if(enable_ma)
				{
					pred+=nlevels>>1;
					pred&=nlevels-1;
					pred-=nlevels>>1;
				}
				((int32_t*)HDR_DATA(dst))[idx<<2|kc]=pred;
			}
		}
	}
	memcpy(((int32_t*)HDR_DATA(src)), ((int32_t*)HDR_DATA(dst)), (size_t)src->iw*src->ih*sizeof(int[4]));
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
	int success, it, kx, ky;

	success=1;
	for(it=0;it<size;++it)
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
		for(kx=it;kx<(size<<1);++kx)
			matrix[((size_t)size<<1)*it+kx]/=pivot;
		for(ky=0;ky<size;++ky)
		{
			double factor;

			if(ky==it)
				continue;
			factor=matrix[((size_t)size<<1)*ky+it];
			if(fabs(factor)>1e-6)
			{
				for(kx=it;kx<(size<<1);++kx)
					matrix[((size_t)size<<1)*ky+kx]-=matrix[((size_t)size<<1)*it+kx]*factor;
			}
		}
	}
	return success;
}
#define OLS_NPARAMS 4	/*8*/
#define OLS_NSAMPLES 17
void pred_ols(Image *src, int fwd, int enable_ma)
{
	double c1_params[OLS_NPARAMS*4]={0};
	Image *dst;
	const int *pixels;
	int successcount, c1_init[4]={0}, kc, ky, kx, idx, ky2, kx2, j, k;

	successcount=0;
	dst=0;
	image_copy(&dst, src);
	if(!dst)
	{
		CRASH("Alloc error");
		return;
	}
	pixels=fwd?((int32_t*)HDR_DATA(src)):((int32_t*)HDR_DATA(dst));
	for(kc=0;kc<4;++kc)
	{
		if(!src->depth[kc])
			continue;
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx, ++idx)
			{
				const int priority[]=
				{
					1, 1, 1, 1, 1, 1, 1,
					1, 1, 2, 4, 2, 1, 1,
					1, 1, 4,
#if 0
					1, 1, 1, 1, 1, 1, 1,
					1, 1, 1, 3, 1, 1, 1,
					1, 1, 3,
#endif
#if 0
					1, 2, 4,  8, 4, 2, 1,
					2, 4, 8, 16, 8, 4, 2,
					4, 8, 16,
#endif
				};
				int pred, success;
#define LOAD(BUF, X, Y) ((uint32_t)(ky-(Y))<(uint32_t)src->ih&&(uint32_t)(kx-(X))<(uint32_t)src->iw?BUF[(src->iw*(ky-(Y))+kx-(X))<<2|kc]:0)
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
					WWWW	=LOAD(pixels,  4, 0),
					WWW	=LOAD(pixels,  3, 0),
					WW	=LOAD(pixels,  2, 0),
					W	=LOAD(pixels,  1, 0);
				int nb[68];
				double nb2[_countof(nb)], sqm[OLS_NPARAMS*(OLS_NPARAMS<<1)]={0}, temp[OLS_NPARAMS<<1];
#undef  LOAD
				/*
				NNNWWWW NNNWWW NNNWW NNNW NNN  NNNE NNNEE NNNEEE NNNEEEE
				NNWWWW  NNWWW  NNWW  NNW  NN   NNE  NNEE  NNEEE  NNEEEE
				NWWWW   NWWW   NWW   NW   N    NE   NEE   NEEE   NEEEE
				WWWW    WWW    WW    W    curr
				*/
				/*slightly better than clampgrad*/
#if 1
				j=0;
				/*NW		N		NE		W		curr*/
				nb[j++]=NNNWWWW;nb[j++]=NNNWWW;	nb[j++]=NNNWW;	nb[j++]=NNWWWW;	/*NNWWW	*/
				nb[j++]=NNNWWW;	nb[j++]=NNNWW;	nb[j++]=NNNW;	nb[j++]=NNWWW;	/*NNWW	*/
				nb[j++]=NNNWW;	nb[j++]=NNNW;	nb[j++]=NNN;	nb[j++]=NNWW;	/*NNW	*/
				nb[j++]=NNNW;	nb[j++]=NNN;	nb[j++]=NNNE;	nb[j++]=NNW;	/*NN	*/
				nb[j++]=NNN;	nb[j++]=NNNE;	nb[j++]=NNNEE;	nb[j++]=NN;	/*NNE	*/
				nb[j++]=NNNE;	nb[j++]=NNNEE;	nb[j++]=NNNEEE;	nb[j++]=NNE;	/*NNEE	*/
				nb[j++]=NNNEE;	nb[j++]=NNNEEE;	nb[j++]=NNNEEEE;nb[j++]=NNEE;	/*NNEEE	*/
				nb[j++]=NNWWWW;	nb[j++]=NNWWW;	nb[j++]=NNWW;	nb[j++]=NWWWW;	/*NWWW	*/
				nb[j++]=NNWWW;	nb[j++]=NNWW;	nb[j++]=NNW;	nb[j++]=NWWW;	/*NWW	*/
				nb[j++]=NNWW;	nb[j++]=NNW;	nb[j++]=NN;	nb[j++]=NWW;	/*NW	*/
				nb[j++]=NNW;	nb[j++]=NN;	nb[j++]=NNE;	nb[j++]=NW;	/*N	*/
				nb[j++]=NN;	nb[j++]=NNE;	nb[j++]=NNEE;	nb[j++]=N;	/*NE	*/
				nb[j++]=NNE;	nb[j++]=NNEE;	nb[j++]=NNEEE;	nb[j++]=NE;	/*NEE	*/
				nb[j++]=NNEE;	nb[j++]=NNEEE;	nb[j++]=NNEEEE;	nb[j++]=NEE;	/*NEEE	*/
				nb[j++]=NWWWW;	nb[j++]=NWWW;	nb[j++]=NWW;	nb[j++]=WWWW;	/*WWW	*/
				nb[j++]=NWWW;	nb[j++]=NWW;	nb[j++]=NW;	nb[j++]=WWW;	/*WW	*/
				nb[j++]=NWW;	nb[j++]=NW;	nb[j++]=N;	nb[j++]=WW;	/*W	*/
#endif

				/*worse than clampgrad*/
#if 0
				/*N+eN,			W+eW,			(N+W)/2,		(4*(N+W)+NW-NW)/8,					curr*/
				NNNWWW	+eNNNWWW,	NNWWWW	+eNNWWWW,	(NNNWWW	+NNWWWW	)>>1,	(4*(NNNWWW	+NNWWWW	)+NNNWW		-NNNWWWW)>>3,/*NNWWW	*/
				NNNWW	+eNNNWW,	NNWWW	+eNNWWW,	(NNNWW	+NNWWW	)>>1,	(4*(NNNWW	+NNWWW	)+NNNW		-NNNWWW	)>>3,/*NNWW	*/
				NNNW	+eNNNW,		NNWW	+eNNWW,		(NNNW	+NNWW	)>>1,	(4*(NNNW	+NNWW	)+NNN		-NNNWW	)>>3,/*NNW	*/
				NNN	+eNNN,		NNW	+eNNW,		(NNN	+NNW	)>>1,	(4*(NNN		+NNW	)+NNNE		-NNNW	)>>3,/*NN	*/
				NNNE	+eNNNE,		NN	+eNN,		(NNNE	+NN	)>>1,	(4*(NNNE	+NN	)+NNNEE		-NNN	)>>3,/*NNE	*/
				NNNEE	+eNNNEE,	NNE	+eNNE,		(NNNEE	+NNE	)>>1,	(4*(NNNEE	+NNE	)+NNNEEE	-NNNE	)>>3,/*NNEE	*/
				NNNEEE	+eNNNEEE,	NNEE	+eNNEE,		(NNNEEE	+NNEE	)>>1,	(4*(NNNEEE	+NNEE	)+NNNEEEE	-NNNEE	)>>3,/*NNEEE	*/
				NNWWW	+eNNWWW,	NWWWW	+eNWWWW,	(NNWWW	+NWWWW	)>>1,	(4*(NNWWW	+NWWWW	)+NNWW		-NNWWWW	)>>3,/*NWWW	*/
				NNWW	+eNNWW,		NWWW	+eNWWW,		(NNWW	+NWWW	)>>1,	(4*(NNWW	+NWWW	)+NNW		-NNWWW	)>>3,/*NWW	*/
				NNW	+eNNW,		NWW	+eNWW,		(NNW	+NWW	)>>1,	(4*(NNW		+NWW	)+NN		-NNWW	)>>3,/*NW	*/
				NN	+eNN,		NW	+eNW,		(NN	+NW	)>>1,	(4*(NN		+NW	)+NNE		-NNW	)>>3,/*N	*/
				NNE	+eNNE,		N	+eN,		(NNE	+N	)>>1,	(4*(NNE		+N	)+NNEE		-NN	)>>3,/*NE	*/
				NNEE	+eNNEE,		NE	+eNE,		(NNEE	+NE	)>>1,	(4*(NNEE	+NE	)+NNEEE		-NNE	)>>3,/*NEE	*/
				NNEEE	+eNNEEE,	NEE	+eNEE,		(NNEEE	+NEE	)>>1,	(4*(NNEEE	+NEE	)+NNEEEE	-NNEE	)>>3,/*NEEE	*/
				NWWW	+eNWWW,		WWWW	+eWWWW,		(NWWW	+WWWW	)>>1,	(4*(NWWW	+WWWW	)+NWW		-NWWWW	)>>3,/*WWW	*/
				NWW	+eNWW,		WWW	+eWWW,		(NWW	+WWW	)>>1,	(4*(NWW		+WWW	)+NW		-NWWW	)>>3,/*WW	*/
				NW	+eNW,		WW	+eWW,		(NW	+WW	)>>1,	(4*(NW		+WW	)+N		-NWW	)>>3,/*W	*/
#endif

				/*X*/
#if 0
				*/NW	N	NE	W		eNW		eN		eNE		eW		curr*/
				NNNWWWW,NNNWWW,	NNNWW,	NNWWWW,		eNNNWWWW,	eNNNWWW,	eNNNWW,		eNNWWWW,	/*NNWWW	*/
				NNNWWW,	NNNWW,	NNNW,	NNWWW,		eNNNWWW,	eNNNWW,		eNNNW,		eNNWWW,		/*NNWW	*/
				NNNWW,	NNNW,	NNN,	NNWW,		eNNNWW,		eNNNW,		eNNN,		eNNWW,		/*NNW	*/
				NNNW,	NNN,	NNNE,	NNW,		eNNNW,		eNNN,		eNNNE,		eNNW,		/*NN	*/
				NNN,	NNNE,	NNNEE,	NN,		eNNN,		eNNNE,		eNNNEE,		eNN,		/*NNE	*/
				NNNE,	NNNEE,	NNNEEE,	NNE,		eNNNE,		eNNNEE,		eNNNEEE,	eNNE,		/*NNEE	*/
				NNNEE,	NNNEEE,	NNNEEEE,NNEE,		eNNNEE,		eNNNEEE,	eNNNEEEE,	eNNEE,		/*NNEEE	*/
				NNWWWW,	NNWWW,	NNWW,	NWWWW,		eNNWWWW,	eNNWWW,		eNNWW,		eNWWWW,		/*NWWW	*/
				NNWWW,	NNWW,	NNW,	NWWW,		eNNWWW,		eNNWW,		eNNW,		eNWWW,		/*NWW	*/
				NNWW,	NNW,	NN,	NWW,		eNNWW,		eNNW,		eNN,		eNWW,		/*NW	*/
				NNW,	NN,	NNE,	NW,		eNNW,		eNN,		eNNE,		eNW,		/*N	*/
				NN,	NNE,	NNEE,	N,		eNN,		eNNE,		eNNEE,		eN,		/*NE	*/
				NNE,	NNEE,	NNEEE,	NE,		eNNE,		eNNEE,		eNNEEE,		eNE,		/*NEE	*/
				NNEE,	NNEEE,	NNEEEE,	NEE,		eNNEE,		eNNEEE,		eNNEEEE,	eNEE,		/*NEEE	*/
				NWWWW,	NWWW,	NWW,	WWWW,		eNWWWW,		eNWWW,		eNWW,		eWWWW,		/*WWW	*/
				NWWW,	NWW,	NW,	WWW,		eNWWW,		eNWW,		eNW,		eWWW,		/*WW	*/
				NWW,	NW,	N,	WW,		eNWW,		eNW,		eN,		eWW,		/*W	*/
#endif
				for(ky2=0;ky2<4;++ky2)
				{
					for(kx2=0;kx2<OLS_NSAMPLES;++kx2)
						nb2[OLS_NSAMPLES*ky2+kx2]=(double)nb[kx2<<2|ky2]/256;/*/(1<<24);*/
				}
				for(ky2=0;ky2<OLS_NPARAMS;++ky2)
				{
					for(kx2=0;kx2<OLS_NPARAMS;++kx2)
					{
						double sum=0;

						for(j=0;j<OLS_NSAMPLES;++j)
							sum+=priority[j]*nb2[OLS_NSAMPLES*ky2+j]*nb2[OLS_NSAMPLES*kx2+j];
						sqm[(OLS_NPARAMS<<1)*ky2+kx2]=sum;
					}
					sqm[((OLS_NPARAMS<<1)+1)*ky2+OLS_NPARAMS]=1;
				}
				success=invert_matrix(sqm, OLS_NPARAMS, temp);
				if(success||c1_init[kc])
				{
					/*pred = nb * params,		params += ((inv(NB * NBT) * NB * Y) - params)*lr*/
					if(success)
					{
						int targets[17];
						
						j=0;
						targets[j++]=NNWWW;
						targets[j++]=NNWW;
						targets[j++]=NNW;
						targets[j++]=NN;
						targets[j++]=NNE;
						targets[j++]=NNEE;
						targets[j++]=NNEEE;
						targets[j++]=NWWW;
						targets[j++]=NWW;
						targets[j++]=NW;
						targets[j++]=N;
						targets[j++]=NE;
						targets[j++]=NEE;
						targets[j++]=NEEE;
						targets[j++]=WWW;
						targets[j++]=WW;
						targets[j++]=W;
						for(ky2=0;ky2<OLS_NPARAMS;++ky2)
						{
							double sum=0;

							for(kx2=0;kx2<OLS_NSAMPLES;++kx2)
								sum+=priority[kx2]*nb2[OLS_NSAMPLES*ky2+kx2]*targets[kx2]/256;/*/(1<<24);*/
							temp[ky2]=sum;
						}
						for(ky2=0;ky2<OLS_NPARAMS;++ky2)
						{
							double sum=0;

							for(j=0;j<OLS_NPARAMS;++j)
								sum+=sqm[(OLS_NPARAMS<<1)*ky2+j+OLS_NPARAMS]*temp[j];
							if(c1_init[kc])
								c1_params[kc<<2|ky2]+=(sum-c1_params[kc<<2|ky2])*0.2;
							else
								c1_params[kc<<2|ky2]=sum;
						}
						c1_init[kc]=1;
					}
					{
						double fpred=0;
						int nb3[4];
						
						nb3[0]=NW;	nb3[1]=N;	nb3[2]=NE;	nb3[3]=W;	/*slightly better than clampgrad*/
#if 0
						N+eN,	W+eW,	(N+W)>>1,	(4*(N+W)+NW-NW)>>3,		/*worse than clampgrad*/

						NW,	N,	NE,	W,	eNW,	eN,	eNE,	eW,	/*X*/
#endif
						for(k=0;k<OLS_NPARAMS;++k)
							fpred+=nb3[k]*c1_params[kc<<2|k];
						pred=(int)TRUNC64(fpred);
					}
					++successcount;
				}
				else/*singular matrix*/
				{
					pred=N+W-NW;
					pred=MEDIAN3(N, W, pred);
				}
				{
					int curr=pred;
					curr^=-fwd;
					curr+=fwd;
					curr+=((int32_t*)HDR_DATA(src))[idx<<2|kc];
					if(enable_ma)
					{
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
					}
					((int32_t*)HDR_DATA(dst))[idx<<2|kc]=curr;
				}
			}
		}
	}
	memcpy(((int32_t*)HDR_DATA(src)), ((int32_t*)HDR_DATA(dst)), (size_t)src->iw*src->ih*sizeof(int[4]));
	if(loud_transforms)
		set_window_title("OLS %f%%", 100.*successcount/(src->nch*src->iw*src->ih));
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
#define OLS2_STEP 2	/*power-of-two*/
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
	int ky3, kx3;

	if(remove)
	{
		for(ky3=0;ky3<OLS2_NPARAMS;++ky3)
		{
			double vy=sample[ky3];
			for(kx3=0;kx3<OLS2_NPARAMS;++kx3)
			{
				double vx=sample[kx3];
				matrix[(OLS2_NPARAMS<<1)*ky3+kx3]-=vy*vx;
			}
		}
	}
	{
		double gain;
		int half, idx2;

		half=(1<<depth)>>1;
		gain=1./half;
		sample[0]=1;
		for(ky3=-OLS2_REACH, idx2=1;ky3<=0;++ky3)
		{
			for(kx3=-OLS2_REACH;kx3<=OLS2_REACH;++kx3, ++idx2)
			{
				sample[idx2]=(double)pixels[(iw*(ky+ky3)+kx+kx3)<<2|kc]*gain;
				if(!ky3&&!kx3)/*last element is target*/
					break;
			}
		}
	}
	if(!remove)
	{
		for(ky3=0;ky3<OLS2_NPARAMS;++ky3)
		{
			double vy=sample[ky3];
			for(kx3=0;kx3<OLS2_NPARAMS;++kx3)
			{
				double vx=sample[kx3];
				matrix[(OLS2_NPARAMS<<1)*ky3+kx3]+=vy*vx;
			}
		}
	}
}
void pred_ols2(Image *src, int fwd, int enable_ma)
{
	Image *dst;
	double *samples, *matrix, *matrix2, *temp;
	const int *pixels;
	int successcount, j, k, kc, ky, kx, idx, ky2, kx2, idx2, ky3, kx3;

	successcount=0;
	dst=0;
	image_copy(&dst, src);
	samples=(double*)malloc(sizeof(double[OLS2_NSAMPLES*(OLS2_NPARAMS+1)]));
	matrix=(double*)malloc(sizeof(double[OLS2_NPARAMS*OLS2_NPARAMS<<1]));
	matrix2=(double*)malloc(sizeof(double[OLS2_NPARAMS*OLS2_NPARAMS<<1]));
	temp=(double*)malloc(sizeof(double[OLS2_NPARAMS<<1]));
	if(!dst||!samples||!matrix||!matrix2||!temp)
	{
		CRASH("Alloc error");
		return;
	}
	memset(samples, 0, sizeof(double[OLS2_NSAMPLES*(OLS2_NPARAMS+1)]));
	pixels=fwd?((int32_t*)HDR_DATA(src)):((int32_t*)HDR_DATA(dst));
	for(kc=0;kc<4;++kc)
	{
		int initialized,
			nlevels, half;
		double params[OLS2_NPARAMS]={0};

		if(!src->depth[kc])
			continue;
		initialized=0;
		nlevels=1<<src->depth[kc], half=nlevels>>1;
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx, ++idx)
			{
				ALIGN(16) int pred[4];

				if(kx<OLS2_SAMPLEREACH+OLS2_REACH||kx>src->iw-OLS2_SAMPLEREACH-OLS2_REACH-OLS2_STEP||ky<OLS2_SAMPLEREACH+OLS2_REACH)
					pred[0]=ols2_fallbackpred(pixels, src->iw, kc, kx, ky, idx);
				else
				{
					int success=0;
					if(kx==OLS2_SAMPLEREACH+OLS2_REACH)/*initialize row*/
					{
						memset(matrix, 0, sizeof(double[OLS2_NPARAMS*OLS2_NPARAMS<<1]));
						for(k=0;k<OLS2_NPARAMS;++k)
						{
							matrix[(OLS2_NPARAMS<<1|1)*k]=0.00005;/*add small lambda for regularization*/
							matrix[(OLS2_NPARAMS<<1|1)*k+OLS2_NPARAMS]=1;/*identity*/
						}
						for(ky2=-OLS2_SAMPLEREACH, idx2=0;ky2<=0;++ky2)
						{
							for(kx2=-OLS2_SAMPLEREACH;kx2<=OLS2_SAMPLEREACH+OLS2_STEP-1&&(ky2||kx2);++kx2, ++idx2)
								ols2_add_sample(pixels, src->depth[kc], src->iw, kc, kx+kx2, ky+ky2, samples+(OLS2_NPARAMS+1)*idx2, matrix, 0);
						}

						memcpy(matrix2, matrix, sizeof(double[OLS2_NPARAMS*OLS2_NPARAMS<<1]));
						success=invert_matrix(matrix2, OLS2_NPARAMS, temp);
					}
					else if(!initialized||!(kx%OLS2_STEP))
					{
						for(ky2=0;ky2<OLS2_SAMPLEREACH-1;++ky2)
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
					/*
					NNWW NNW NN NNE NNEE
					NWW  NW  N  NE  NEE
					WW   W   ?
					*/
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
							/*params = ((inv(NB * NBT) * NB * Targets) - params)*lr*/
							for(ky2=0;ky2<OLS2_NPARAMS;++ky2)/*temp = NB * Targets*/
							{
								double sum=0;
								for(kx2=0;kx2<OLS2_NSAMPLES;++kx2)
									sum+=samples[(OLS2_NPARAMS+1)*kx2+ky2]*samples[(OLS2_NPARAMS+1)*kx2+OLS2_NPARAMS];
								temp[ky2]=sum;
							}
							for(ky2=0;ky2<OLS2_NPARAMS;++ky2)/*params = matrix * temp*/
							{
								double sum=0;
								for(j=0;j<OLS2_NPARAMS;++j)
									sum+=matrix2[(OLS2_NPARAMS<<1)*ky2+j+OLS2_NPARAMS]*temp[j];
								if(initialized)
									params[ky2]+=(sum-params[ky2])*0.2;
								else
									params[ky2]=sum;
							}
							initialized=1;
						}
						/*pred = nb * params*/
						fpred=params[0]*half;
						for(ky3=-OLS2_REACH, idx2=1;ky3<=0;++ky3)
						{
							for(kx3=-OLS2_REACH;kx3<=OLS2_REACH;++kx3, ++idx2)
							{
								if(!ky3&&!kx3)/*exclude current pixel*/
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
						++successcount;
					}
				}
				{
					int curr=pred[0];
					curr^=-fwd;
					curr+=fwd;
					curr+=((int32_t*)HDR_DATA(src))[idx<<2|kc];
					if(enable_ma)
					{
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
					}
					((int32_t*)HDR_DATA(dst))[idx<<2|kc]=curr;
				}
			}
		}
	}
	memcpy(((int32_t*)HDR_DATA(src)), ((int32_t*)HDR_DATA(dst)), (size_t)src->iw*src->ih*sizeof(int[4]));
	if(loud_transforms)
		set_window_title("OLS2 %f%%", 100.*successcount/(src->nch*src->iw*src->ih));
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
	int nparams, ky3, kx3, kc2, idx2;

	nparams=OLS3_NPARAMS0+kc;
	if(remove)
	{
		for(ky3=0;ky3<nparams;++ky3)
		{
			double vy=sample[ky3];
			for(kx3=0;kx3<nparams;++kx3)
			{
				double vx=sample[kx3];
				matrix[(nparams<<1)*ky3+kx3]-=vy*vx;
			}
		}
	}
	sample[0]=1;
	for(ky3=-OLS3_REACH, idx2=1;ky3<=0;++ky3)
	{
		for(kx3=-OLS3_REACH;kx3<=OLS3_REACH;++kx3)
		{
			for(kc2=0;kc2<3;++kc2, ++idx2)
			{
				sample[idx2]=(double)pixels[(iw*(ky+ky3)+kx+kx3)<<2|kc2]*gain;
				if(!ky3&&!kx3&&kc2==kc)/*last element is target*/
					goto finish;
			}
		}
	}
finish:
	if(!remove)
	{
		for(ky3=0;ky3<nparams;++ky3)
		{
			double vy=sample[ky3];
			for(kx3=0;kx3<nparams;++kx3)
			{
				double vx=sample[kx3];
				matrix[(nparams<<1)*ky3+kx3]+=vy*vx;
			}
		}
	}
}
void pred_ols3(Image *src, int fwd, int enable_ma)
{
	double t_start, *samples, *matrix, *matrix2, *temp, gains[4];
	double avparams[OLS3_NPARAMS2*3]={0};
	double params[OLS3_NPARAMS2*3]={0};
	Image *dst;
	const int *pixels;
	int nlevels[4], half[4], successcount, j, k, ky, kx, kc, idx, ky2, kx2, kc2, idx2, ky3, kx3;

	t_start=time_sec();
	if(loud_transforms)
		DisableProcessWindowsGhosting();
	nlevels[0]=1<<src->depth[0];
	nlevels[1]=1<<src->depth[1];
	nlevels[2]=1<<src->depth[2];
	nlevels[3]=1<<src->depth[3];
	half[0]=nlevels[0]>>1;
	half[1]=nlevels[1]>>1;
	half[2]=nlevels[2]>>1;
	half[3]=nlevels[3]>>1;
	gains[0]=1./half[0];
	gains[1]=1./half[1];
	gains[2]=1./half[2];
	gains[3]=1./half[3];
	successcount=0;
	dst=0;
	image_copy(&dst, src);
	samples=(double*)malloc(sizeof(double[OLS3_NSAMPLES*(OLS3_NPARAMS2+1)*3]));/*FIXME FIXME FIXME FIXME arena*/
	matrix=(double*)malloc(sizeof(double[OLS3_NPARAMS2*OLS3_NPARAMS2*3<<1]));
	matrix2=(double*)malloc(sizeof(double[OLS3_NPARAMS2*OLS3_NPARAMS2*3<<1]));
	temp=(double*)malloc(sizeof(double[OLS3_NPARAMS2<<1]));
	if(!dst||!samples||!matrix||!matrix2||!temp)
	{
		CRASH("Alloc error");
		return;
	}
	memset(samples, 0, sizeof(double[OLS3_NSAMPLES*(OLS3_NPARAMS2+1)*3]));
	pixels=fwd?((int32_t*)HDR_DATA(dst)):((int32_t*)HDR_DATA(src));
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		int initialized[3]={0};
		for(kx=0;kx<src->iw;++kx, ++idx)
		{
			for(kc=0;kc<3;++kc)
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
					if(!(kx%OLS3_STEP))
					{
						memset(mat1, 0, sizeof(double[2])*nparams*nparams);
						for(k=0;k<nparams;++k)
						{
							mat1[(nparams<<1|1)*k]=0.00005;/*add small lambda for regularization*/
							mat1[(nparams<<1|1)*k+nparams]=1;/*identity*/
						}
						for(ky2=-OLS3_SAMPLEREACH, idx2=0;ky2<=0;++ky2)
						{
							for(kx2=-OLS3_SAMPLEREACH;kx2<=OLS3_SAMPLEREACH+OLS3_STEP-1&&(ky2||kx2);++kx2, ++idx2)
								ols3_add_sample(pixels, gains[kc], src->iw, kc, kx+kx2, ky+ky2, csamples+(OLS3_NPARAMS2+1)*idx2, mat1, 0);
						}

						memcpy(mat2, mat1, sizeof(double[2])*nparams*nparams);
						success=invert_matrix(mat2, nparams, temp);
					}
					/*
					NNWW NNW NN NNE NNEE
					NWW  NW  N  NE  NEE
					WW   W   ?
					*/
					if(!success&&!initialized[kc])
						pred[0]=ols3_fallbackpred(pixels, src->iw, kc, kx, ky, idx);
					else
					{
						if(success)
						{
							/*params = ((inv(NB * NBT) * NB * Targets) - params)*lr*/
							for(ky2=0;ky2<nparams;++ky2)/*temp = NB * Targets*/
							{
								double sum=0;
								for(kx2=0;kx2<OLS3_NSAMPLES;++kx2)
									sum+=csamples[(OLS3_NPARAMS2+1)*kx2+ky2]*csamples[(OLS3_NPARAMS2+1)*kx2+nparams];
								temp[ky2]=sum;
							}
							for(ky2=0;ky2<nparams;++ky2)/*params = matrix * temp*/
							{
								double sum=0;
								for(j=0;j<nparams;++j)
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
							/*pred = nb * params*/
							int
								N =pixels[(idx-src->iw)<<2|kc],
								W =pixels[(idx-1)<<2|kc],
								NE=pixels[(idx-src->iw+1)<<2|kc];
							double fpred=cparams[0]/gains[kc];
							for(ky3=-OLS3_REACH, idx2=1;ky3<=0;++ky3)
							{
								for(kx3=-OLS3_REACH;kx3<=OLS3_REACH;++kx3)
								{
									for(kc2=0;kc2<3;++kc2, ++idx2)
									{
										if(!ky3&&!kx3&&kc2==kc)/*exclude current pixel*/
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
							++successcount;
						}
					}
				}
				{
					int curr=pred[0];
					curr^=-fwd;
					curr+=fwd;
					curr+=((int32_t*)HDR_DATA(src))[idx<<2|kc];
					if(enable_ma)
					{
						curr<<=32-src->depth[kc];
						curr>>=32-src->depth[kc];
					}
					((int32_t*)HDR_DATA(src))[idx<<2|kc]=curr;
				}
			}
		}
		if(loud_transforms)
			set_window_title("%d/%d = %7.3f%%  OLS3 rate %f%%  %f sec", ky+1, src->ih, 100.*(ky+1)/src->ih, 100.*successcount/(src->nch*src->iw*(ky+1)), time_sec()-t_start);
	}
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
		int kg, printed;

		printed=0;
		set_window_title("OLS3 %f%%  %f sec", 100.*successcount/(src->nch*src->iw*src->ih), time_sec()-t_start);
		for(k=0;k<(int)_countof(avparams);++k)
			avparams[k]*=(double)src->nch/successcount;
		printed+=snprintf(g_buf+printed, G_BUF_SIZE-1-printed, "Y U V\n");
		for(kg=0, idx=0;kg<3;++kg)
		{
			idx2=0;
			for(k=0;k<OLS3_NPARAMS2;++k)
			{
				printed+=snprintf(g_buf+printed, G_BUF_SIZE-1-printed, " %+12.8f", avparams[idx++]);
				if(!(k%3))
					printed+=snprintf(g_buf+printed, G_BUF_SIZE-1-printed, "  %s\n", nbnames[idx2++]);
			}
			printed+=snprintf(g_buf+printed, G_BUF_SIZE-1-printed, "  %s", nbnames[idx2]);
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

#define LIN_NNB (2*(LIN_REACH+1)*LIN_REACH)
void pred_linear(Image const *src, Image *dst, const int *coeffs, int lgden, int fwd, int enable_ma)
{
	const int *pixels;
	int k, kc, ky, kx, ky2, kx2;

	pixels=fwd?((const int32_t*)HDR_DATC(src)):((const int32_t*)HDR_DATC(dst));
	for(kc=0;kc<4;++kc)
	{
		if(!src->depth[kc])
			continue;
		for(ky=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx)
			{
				int nb[LIN_NNB], idx;
				int pred;

				idx=0;
				for(ky2=-LIN_REACH;ky2<=0;++ky2)
				{
					for(kx2=-LIN_REACH;kx2<=LIN_REACH;++kx2, ++idx)
					{
						if(!ky2&&!kx2)
							break;
						if((uint32_t)(ky+ky2)<(uint32_t)src->ih&&(uint32_t)(kx+kx2)<(uint32_t)src->iw)
							nb[idx]=pixels[(src->iw*(ky+ky2)+kx+kx2)<<2|kc];
						else
							nb[idx]=0;
					}
				}
				pred=(1<<lgden)>>1;
				for(k=0;k<LIN_NNB;++k)
					pred+=coeffs[k]*nb[k];
				pred>>=lgden;
				idx=(src->iw*ky+kx)<<2|kc;

				pred^=-fwd;
				pred+=fwd;
				pred+=((const int32_t*)HDR_DATC(src))[idx];
				if(enable_ma)
				{
					pred<<=32-src->depth[kc];
					pred>>=32-src->depth[kc];
				}
				((int32_t*)HDR_DATA(dst))[idx]=pred;
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
	SLIC5_PRED(N-(int)(((int64_t)eN+eW+eNE)*-0x05C>>8))\
	SLIC5_PRED(W-(int)(((int64_t)eN+eW+eNW)*-0x05B>>8))\
	SLIC5_PRED(N+(int)((-eNN*(intptr_t)0x0DF-eN*(intptr_t)0x051-eNE*(intptr_t)0x0BD+((int64_t)N-NN)*0x05C+((int64_t)NW-W)*0x102)>>8))\
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
	SLIC5_PRED((NE+NNE+NEE+NNEE)>>2)\

#if 0
	SLIC5_PRED(geomean)\

#endif
void pred_t47(Image *src, int fwd, int enable_ma)
{
	Image *dst;
	const int *pixels, *errors;
	int pesize, *pred_errors, k, kc, ky, kx, idx;

	dst=0;
	image_copy(&dst, src);
	pesize=(src->iw+4)*(int)sizeof(int[SLIC5_NPREDS*4]);
	pred_errors=(int*)malloc(pesize);/*4 padded rows*/
	if(!dst||!pred_errors)
	{
		CRASH("Alloc error");
		return;
	}
	pixels=fwd?((int32_t*)HDR_DATA(dst)):((int32_t*)HDR_DATA(src));
	errors=fwd?((int32_t*)HDR_DATA(src)):((int32_t*)HDR_DATA(dst));
	for(kc=0;kc<4;++kc)
	{
		int sh;
		int params[SLIC5_NPREDS], param0=176;
		int64_t bias_sum=0;
		int bias_count=0;

		if(!src->depth[kc])
			continue;
		sh=src->depth[kc];
		memfill(params, &param0, sizeof(params), sizeof(int));
		memset(pred_errors, 0, pesize);
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			int kym[4];

			kym[0]=(src->iw+4)*((ky+0)&3);
			kym[1]=(src->iw+4)*((ky-1)&3);
			kym[2]=(src->iw+4)*((ky-2)&3);
			kym[3]=(src->iw+4)*((ky-3)&3);
			for(kx=0;kx<src->iw;++kx, ++idx)
			{
#define LOAD(BUF, X, Y) ((uint32_t)(ky+(Y))<(uint32_t)src->ih&&(uint32_t)(kx+(X))<(uint32_t)src->iw?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]<<8:0)
				int
					NNNWWWW	=LOAD(pixels, -4, -3),
					NNNWWW	=LOAD(pixels, -3, -3),
					NNN	=LOAD(pixels,  0, -3),
					NNNEE	=LOAD(pixels,  2, -3),
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

					eNN	=LOAD(errors,  0, -2),
					eNNE	=LOAD(errors,  1, -2),
					eNW	=LOAD(errors, -1, -1),
					eN	=LOAD(errors,  0, -1),
					eNE	=LOAD(errors,  1, -1),
					eWW	=LOAD(errors, -2,  0),
					eW	=LOAD(errors, -1,  0);
#undef  LOAD
				int paper_GAP, calic_GAP;
				/*int geomean;*/
				{
					int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
					int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
					int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
					int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
					int diff=(dy-dx)>>sh, diff2=(d45-d135)>>sh, diff3=NE-NW;
					if(dy+dx>32)
						paper_GAP=(int)(((int64_t)dx*N+(int64_t)dy*W)/((int64_t)dy+dx));
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
						calic_GAP=(2*N+6*W+NE-NW)>>3;		/*c1	[1/4  3/4  1/8  -1/8].[N W NE NW]*/
					else if(diff>8)
						calic_GAP=(6*N+10*W+3*(NE-NW))>>4;	/*c2	[3/8  5/8  3/16  -3/16]*/
					else if(diff<-32)
						calic_GAP=(6*N+2*W+NE-NW)>>3;		/*c3	[3/4  1/4  1/8  -1/8]*/
					else if(diff<-8)
						calic_GAP=(10*N+6*W+3*(NE-NW))>>4;	/*c4	[5/8  3/8  3/16  -3/16]*/
					else
						calic_GAP=(((N+W)<<1)+NE-NW)>>2;	/*c5	[1/2  1/2  1/4  -1/4]*/
				}
				{
					int64_t pred=0, wsum=0;
					int val, curr, kbest, besterr, preds[32], j;

#define SLIC5_PRED(X) preds[j++]=X;
					j=0;
					SLIC5_PREDLIST;
#undef  SLIC5_PRED
#define LOAD(mX, mY) pred_errors[SLIC5_NPREDS*(kym[mY]+kx+2-(mX))+k]
					for(k=0;k<(int)_countof(preds);++k)
					{
						int64_t weight=1+
							(int64_t)LOAD(-2,  2)+
							(int64_t)LOAD(-1,  2)+
							(int64_t)LOAD( 0,  2)+
							(int64_t)LOAD( 1,  2)+
							(int64_t)LOAD(-2,  1)+
							(int64_t)LOAD(-1,  1)+
							(int64_t)LOAD( 0,  1)*3+
							(int64_t)LOAD( 1,  1);
						weight=((int64_t)params[k]<<29)/wndw;
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
					pred+=bias_sum/(bias_count+(intptr_t)1);
					CLAMP3_32(pred, (int)pred, N, W, NE);
					val=(int)((pred+127)>>8);

					val^=-fwd;
					val+=fwd;
					val+=((int32_t*)HDR_DATA(src))[idx<<2|kc];
					if(enable_ma)
					{
						val<<=32-sh;
						val>>=32-sh;
					}
					((int32_t*)HDR_DATA(src))[idx<<2|kc]=val;

					curr=pixels[idx<<2|kc]<<8;
					kbest=-1, besterr=0;
					for(k=0;k<(int)_countof(preds);++k)
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
						for(k=0;k<(int)_countof(preds);++k)
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

/*CUSTOM reach-2 predictor*/
int custom_params[CUSTOM_NPARAMS]={0};
int custom_clamp[4]={0};
void pred_custom(Image *src, int fwd, int enable_ma, const int *params)
{
	ALIGN(16) int *rows[4];
	int amin[4], amax[4], bufsize, *pixels, k, ky, kx, kc, idx;

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	bufsize=(src->iw+16)*(int)sizeof(int[4*4*2]);/*4 padded rows * 4 channels max * {pixel, pred}*/
	pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4*2;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4*2;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4*2;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4*2;
		for(kx=0;kx<src->iw;++kx, ++idx)
		{
			for(kc=0;kc<3;++kc)/*there are params only for 3 channels in CUSTOM predictor*/
			{
				const int clamp_idx[]={22, 12, 14, 16};/*W NW N NE*/
				int64_t pred;
				int nb[CUSTOM_NNB*2], j, nlevels, vmin, vmax, uninitialized, curr, val;

				if(!src->depth[kc])
					continue;
				nlevels=1<<src->depth[kc];
				j=0;
				nb[j++]=rows[2][kc-2*8+0]; nb[j++]=rows[2][kc-2*8+4];
				nb[j++]=rows[2][kc-1*8+0]; nb[j++]=rows[2][kc-1*8+4];
				nb[j++]=rows[2][kc+0*8+0]; nb[j++]=rows[2][kc+0*8+4];
				nb[j++]=rows[2][kc+1*8+0]; nb[j++]=rows[2][kc+1*8+4];
				nb[j++]=rows[2][kc+2*8+0]; nb[j++]=rows[2][kc+2*8+4];
				nb[j++]=rows[1][kc-2*8+0]; nb[j++]=rows[1][kc-2*8+4];
				nb[j++]=rows[1][kc-1*8+0]; nb[j++]=rows[1][kc-1*8+4];
				nb[j++]=rows[1][kc+0*8+0]; nb[j++]=rows[1][kc+0*8+4];
				nb[j++]=rows[1][kc+1*8+0]; nb[j++]=rows[1][kc+1*8+4];
				nb[j++]=rows[1][kc+2*8+0]; nb[j++]=rows[1][kc+2*8+4];
				nb[j++]=rows[0][kc-2*8+0]; nb[j++]=rows[0][kc-2*8+4];
				nb[j++]=rows[0][kc-1*8+0]; nb[j++]=rows[0][kc-1*8+4];
				curr=((int32_t*)HDR_DATA(src))[idx<<2|kc];
				vmin=-(nlevels>>1);
				vmax=(nlevels>>1);
				uninitialized=1;
				for(k=0;k<4;++k)
				{
					if(custom_clamp[k])
					{
						val=nb[clamp_idx[k]];
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
				for(k=0;k<CUSTOM_NNB*2;++k)
					pred+=(int64_t)nb[k]*params[kc*CUSTOM_NNB*2+k];
				pred+=(intptr_t)1<<15;
				pred>>=16;
				pred=CLAMP(vmin, pred, vmax);
				if(fwd)
				{
					val=(curr-(int)pred+g_dist/2)/g_dist;
					curr=g_dist*val+(int)pred;
				}
				else
				{
					val=g_dist*curr+(int)pred;
					curr=val;
					if(g_dist==1)
					{
						val<<=32-src->depth[kc];
						val>>=32-src->depth[kc];
					}
					else
						CLAMP2(val, amin[kc], amax[kc]);
				}
				((int32_t*)HDR_DATA(src))[idx<<2|kc]=val;
				rows[0][kc+0]=curr;
				rows[0][kc+4]=curr-(int)pred;
			}
			rows[0]+=8;
			rows[1]+=8;
			rows[2]+=8;
			rows[3]+=8;
		}
	}
	if(!enable_ma)
	{
		++src->depth[0];
		++src->depth[1];
		++src->depth[2];
		src->depth[3]+=src->depth[3]!=0;
	}
	free(pixels);
}
static void pred_custom_calcloss(Image const *src, Image *dst, int *hist, const int *params, double *loss)
{
	double entropy;
	intptr_t res;
	int kc;

	res=(intptr_t)src->iw*src->ih;
	memcpy(((int32_t*)HDR_DATA(dst)), ((const int32_t*)HDR_DATC(src)), res*sizeof(int[4]));

	pred_custom(dst, 1, 1, params);

	for(kc=0;kc<3;++kc)
	{
		int nlevels=1<<src->depth[kc];
		calc_histogram(((int32_t*)HDR_DATA(dst)), dst->iw, dst->ih, kc, 0, dst->iw, 0, dst->ih, src->depth[kc], hist, 0);
		entropy=calc_entropy(hist, nlevels, (int)res);
		loss[kc]=entropy/src->src_depth[(kc+1)%3];
	}
	loss[3]=(loss[0]+loss[1]+loss[2])/3;
}
enum
{
	CUSTOM_NITER=128,
	CUSTOM_DELTAGROUP=2	/*must be 2*/
};
void pred_custom_optimize(Image const *image, int *params, int reach)
{
	static int call_idx=0;
	const int r1_idx[]=
	{
		6, 7, 8,	/*NW N NE*/
		11,		/*W  ?*/
	};
	Image *im2;
	double loss_bestsofar[4], loss_prev[4], loss_curr[4];
	int maxdepth, maxlevels, *hist, params2[CUSTOM_NPARAMS], shakethreshold, k, it, watchdog;

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
		CRASH("Alloc error");
		return;
	}
	memcpy(params2, params, sizeof(params2));
#define CALC_LOSS(L) pred_custom_calcloss(image, im2, hist, params, L)
#ifndef _DEBUG
	srand((uint32_t)__rdtsc());
#endif
	CALC_LOSS(loss_bestsofar);
	{
		int sum=0, idx;
		for(k=0;k<CUSTOM_NNB*2;k+=2)
			sum+=params[CUSTOM_NNB*2*custom_pred_ch_idx+k];
		if(reach==1)
			idx=r1_idx[rand()&3];
		else
			idx=rand()%CUSTOM_NNB;
		params[CUSTOM_NNB*2*custom_pred_ch_idx+idx*2]+=0x10000-sum;
	}
	CALC_LOSS(loss_prev);

	shakethreshold=CUSTOM_NPARAMS;
	for(it=0, watchdog=0;it<CUSTOM_NITER;++it)
	{
		int idx[CUSTOM_DELTAGROUP]={0}, stuck=0;
		int params_original_selected[CUSTOM_DELTAGROUP]={0};

		if(watchdog>=shakethreshold)/*bump if stuck?*/
		{
			watchdog=0;
			stuck=1;
		}
		else if(reach==1)
		{
			int inc;

			inc=(rand()&1?1:-1)<<(rand()&7);
			idx[0]=(CUSTOM_NNB*custom_pred_ch_idx+r1_idx[rand()&3])*2;
			idx[1]=(CUSTOM_NNB*custom_pred_ch_idx+r1_idx[rand()&3])*2;
			params_original_selected[0]=params[idx[0]];
			params[idx[0]]+=inc;
			params_original_selected[1]=params[idx[1]];
			params[idx[1]]-=inc;
		}
		else
		{
			for(k=0;k<CUSTOM_DELTAGROUP;++k)/*increment params*/
			{
				int inc=0;

				inc=rand()&1?0x100:-0x100;

				idx[k]=(CUSTOM_NNB*custom_pred_ch_idx+rand()%CUSTOM_NNB)*2;
				++k;
				idx[k]=(CUSTOM_NNB*custom_pred_ch_idx+rand()%CUSTOM_NNB)*2;

				params_original_selected[k-1]=params[idx[k-1]];
				params[idx[k-1]]+=inc;
				params_original_selected[k+0]=params[idx[k+0]];
				params[idx[k+0]]-=inc;
			}
		}

		CALC_LOSS(loss_curr);
		
		if(loss_prev[3]<loss_curr[3])/*revert if worse*/
		{
			if(stuck)/*a bad branch may surpass the local minimum*/
				memcpy(loss_prev, loss_curr, sizeof(loss_prev));
			else
			{
				memcpy(loss_curr, loss_prev, sizeof(loss_prev));
				for(k=0;k<CUSTOM_DELTAGROUP;++k)
					params[idx[k]]=params_original_selected[k];
			}
			++watchdog;
		}
		else/*save if better*/
		{
			if(loss_curr[3]<loss_bestsofar[3])/*publish if record best*/
			{
				memcpy(params2, params, sizeof(params2));
				memcpy(loss_bestsofar, loss_curr, sizeof(loss_bestsofar));
				--it;/*again*/
			}
			memcpy(loss_prev, loss_curr, sizeof(loss_prev));
			watchdog=0;
		}
		if(loud_transforms)
			set_window_title(
				"%d %4d/%4d,%d/%d: %f%% RGB %f %f %f%s",
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

		/*preview*/
#if 1
		{
			ch_entropy[0]=(float)loss_bestsofar[0];
			ch_entropy[1]=(float)loss_bestsofar[1];
			ch_entropy[2]=(float)loss_bestsofar[2];
			ch_entropy[3]=0;
			io_render();
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
	Image *dst;
	int k, ky, kx, kc;

	dst=0;
	image_copy(&dst, src);
	if(!dst)
	{
		CRASH("Alloc error");
		return;
	}
	for(ky=0;ky<src->ih;++ky)
	{
		for(kx=0;kx<src->iw;++kx)
		{
			for(kc=0;kc<3;++kc)
			{
				int nb[12], j, *curr_params, idx, val;
#define LOAD(BUF, X, Y) ((uint32_t)(kx+(X))<(uint32_t)src->iw&&(uint32_t)(ky+(Y))<(uint32_t)src->ih?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0)
				j=0;
				nb[j++]=LOAD(((int32_t*)HDR_DATA(src)), -2, -2);
				nb[j++]=LOAD(((int32_t*)HDR_DATA(src)), -1, -2);
				nb[j++]=LOAD(((int32_t*)HDR_DATA(src)),  0, -2);
				nb[j++]=LOAD(((int32_t*)HDR_DATA(src)),  1, -2);
				nb[j++]=LOAD(((int32_t*)HDR_DATA(src)),  2, -2);
				nb[j++]=LOAD(((int32_t*)HDR_DATA(src)), -2, -1);
				nb[j++]=LOAD(((int32_t*)HDR_DATA(src)), -1, -1);
				nb[j++]=LOAD(((int32_t*)HDR_DATA(src)),  0, -1);
				nb[j++]=LOAD(((int32_t*)HDR_DATA(src)),  1, -1);
				nb[j++]=LOAD(((int32_t*)HDR_DATA(src)),  2, -1);
				nb[j++]=LOAD(((int32_t*)HDR_DATA(src)), -2,  0);
				nb[j++]=LOAD(((int32_t*)HDR_DATA(src)), -1,  0);
#undef  LOAD
				curr_params=custom_params+kc*CUSTOM_NNB*2;
				if(!curr_params[0])
					continue;
				idx=(src->iw*ky+kx)<<2|kc;
				val=((int32_t*)HDR_DATA(dst))[idx]*curr_params[0];
				for(k=1;k<(int)_countof(nb);++k)
					val+=nb[k]*curr_params[k<<1];
				val+=1<<15;
				val>>=16;
				((int32_t*)HDR_DATA(src))[idx]=val;
			}
		}
	}
	free(dst);
}
void pred_convtest(Image *src)
{
	ALIGN(16) int *rows[4];
	int fwd, bufsize, *pixels, k, ky, kx, kc, idx, nlevels;

	fwd=0;
	bufsize=(src->iw+16)*(int)sizeof(int[4*4*2]);/*4 padded rows * 4 channels max * {pixel, pred}*/
	pixels=(int*)malloc(bufsize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	memset(pixels, 0, bufsize);
	for(ky=0, idx=0;ky<src->ih;++ky)
	{
		rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0)&3)+8)*4*2;
		rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1)&3)+8)*4*2;
		rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2)&3)+8)*4*2;
		rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3)&3)+8)*4*2;
		for(kx=0;kx<src->iw;++kx, ++idx)
		{
			for(kc=0;kc<3;++kc)/*there are params only for 3 channels in CUSTOM predictor*/
			{
				const int clamp_idx[]={22, 12, 14, 16};/*W NW N NE*/
				int nb[24], j, vmin, vmax, uninitialized, curr;
				int64_t pred;

				if(!src->depth[kc])
					continue;
				nlevels=1<<src->depth[kc];
				j=0;
				nb[j++]=rows[2][kc-2*8+0]; nb[j++]=rows[2][kc-2*8+4];
				nb[j++]=rows[2][kc-1*8+0]; nb[j++]=rows[2][kc-1*8+4];
				nb[j++]=rows[2][kc+0*8+0]; nb[j++]=rows[2][kc+0*8+4];
				nb[j++]=rows[2][kc+1*8+0]; nb[j++]=rows[2][kc+1*8+4];
				nb[j++]=rows[2][kc+2*8+0]; nb[j++]=rows[2][kc+2*8+4];
				nb[j++]=rows[1][kc-2*8+0]; nb[j++]=rows[1][kc-2*8+4];
				nb[j++]=rows[1][kc-1*8+0]; nb[j++]=rows[1][kc-1*8+4];
				nb[j++]=rows[1][kc+0*8+0]; nb[j++]=rows[1][kc+0*8+4];
				nb[j++]=rows[1][kc+1*8+0]; nb[j++]=rows[1][kc+1*8+4];
				nb[j++]=rows[1][kc+2*8+0]; nb[j++]=rows[1][kc+2*8+4];
				nb[j++]=rows[0][kc-2*8+0]; nb[j++]=rows[0][kc-2*8+4];
				nb[j++]=rows[0][kc-1*8+0]; nb[j++]=rows[0][kc-1*8+4];
				curr=!(kx&63)&&!(ky&63)?64:0;
				vmin=-(nlevels>>1);
				vmax=(nlevels>>1);
				uninitialized=1;
				for(k=0;k<4;++k)
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
				for(k=0;k<CUSTOM_NNB*2;++k)
					pred+=(int64_t)nb[k]*custom_params[k];
				pred+=(intptr_t)1<<15;
				pred>>=16;
				pred=CLAMP(vmin, pred, vmax);

				rows[0][kc+4*!fwd]=curr;
				pred^=-fwd;
				pred+=fwd;
				pred+=(int64_t)curr;
				pred<<=32-src->depth[kc];
				pred>>=32-src->depth[kc];
				((int32_t*)HDR_DATA(src))[idx<<2|kc]=(int)pred;
				rows[0][kc+4* fwd]=(int)pred;
			}
			rows[0]+=8;
			rows[1]+=8;
			rows[2]+=8;
			rows[3]+=8;
		}
	}
	free(pixels);
}




double pw2_errors[PW2_NPRED]={0};/**/
int16_t pw2_params[PW2_NPARAM*3]=
{
	/*0*/
	
	 0x00EA, 0x01C8, 0x00A2, 0x005E, 0x01F4, 0x0045, 0x0091, 0x0066, 0x003B, 0x0027,-0x0011, 0x001B, 0x00FF, 0x007E, 0x00D1, 0x00F3, 0x008F, 0x0130, 0x018E,-0x00AC, 0x010C, 0x0008,-0x007E, 0x00A2, 0x000E,-0x0069,-0x0073,-0x0125,-0x0092, 0x0000, 0x0078,/*Y	*/
	 0x0006, 0x003D, 0x0031, 0x002F, 0x003F, 0x0015, 0x0011, 0x0036, 0x002E,-0x0022, 0x0011, 0x0034,-0x0007, 0x0012,-0x0018, 0x0012, 0x002F, 0x0000, 0x0000, 0x001C, 0x00A2, 0x02E1, 0x00C9,-0x00E0,-0x0068,-0x004E,-0x013E,-0x0012, 0x0001, 0x0000,-0x0046,/*Cb	*/
	 0x003D, 0x0036, 0x0006, 0x007E, 0x0012, 0x0007, 0x0007, 0x0005, 0x001E, 0x0000, 0x0028, 0x0055,-0x0020, 0x0020, 0x0005, 0x0011, 0x0034, 0x0000, 0x0004, 0x003E,-0x0100, 0x0001,-0x0086,-0x0041, 0x0051,-0x0080, 0x0004, 0x0002,-0x0003,-0x0003, 0x00D9,/*Cr	*/

#if 0
	 0x007D, 0x0040, 0x0039, 0x004A, 0x0007, 0x0003, 0x001D, 0x0007, 0x0000, 0x0007, 0x0002, 0x000F, 0x001B, 0x0018, 0x000A, 0x0008, 0x001C,-0x0008, 0x0004, 0x0005, 0x0006,-0x0022,-0x003B,-0x0041,-0x00C8,-0x0040,-0x0085,-0x0050, 0x0060,
	 0x007A, 0x00F9, 0x0165, 0x00C2, 0x0036, 0x0100, 0x0054, 0x0000,-0x0081,-0x0078, 0x0020, 0x004D,-0x0010, 0x0028, 0x00BD, 0x009D, 0x0020,-0x0082,-0x003F, 0x0060, 0x002A, 0x0161,-0x004E,-0x001D, 0x0123,-0x0008,-0x0080, 0x0020, 0x003C,
	 0x0010, 0x0039, 0x002E, 0x0037, 0x000E,-0x0010, 0x0014, 0x0008,-0x0007,-0x001C, 0x0074, 0x0019, 0x0010, 0x001B, 0x000D, 0x0047, 0x000A, 0x001C, 0x0008, 0x0004, 0x0023,-0x0012,-0x0156,-0x0074,-0x00A0,-0x0002,-0x0088,-0x0060, 0x0102,
#endif
};
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
static void pred_w2_prealloc(const int *src, int iw, int ih, int depth, int kc, int16_t *params, int fwd, int enable_ma, int *dst, int *temp)/*temp is (PW2_NPRED+1)*2w*/
{
	const int *src2;
	int errorbuflen, rowlen, *hireserror, *pred_errors[PW2_NPRED], idx, nlevels, k, ky, kx;

	errorbuflen=iw<<1;
	rowlen=iw<<2;
	hireserror=temp;
	for(k=0;k<PW2_NPRED;++k)
		pred_errors[k]=temp+errorbuflen*(k+1);
	idx=kc;

	nlevels=1<<depth;
	src2=fwd?src:dst;
	memset(pw2_errors, 0, sizeof(pw2_errors));
	for(ky=0;ky<ih;++ky)
	{
		int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
		for(kx=0;kx<iw;++kx, idx+=4)
		{
			int
				cT6  =         ky-6>=0?src2[idx-rowlen*6   ]<<8:0,

				cT3  =         ky-3>=0?src2[idx-rowlen*3   ]<<8:0,

				cT2L =kx-1>=0&&ky-2>=0?src2[idx-rowlen*2- 4]<<8:0,
				cT2  =         ky-2>=0?src2[idx-rowlen*2   ]<<8:0,
				cT2R =kx+1<iw&&ky-2>=0?src2[idx-rowlen*2+ 4]<<8:0,

				cTL  =kx-1>=0&&ky-1>=0?src2[idx-rowlen  - 4]<<8:0,
				cT   =         ky-1>=0?src2[idx-rowlen     ]<<8:0,
				cTR  =kx+1<iw&&ky-1>=0?src2[idx-rowlen  + 4]<<8:0,
				cTR2 =kx+2<iw&&ky-1>=0?src2[idx-rowlen  + 8]<<8:0,
				cTR4 =kx+4<iw&&ky-1>=0?src2[idx-rowlen  +16]<<8:0,
				cTR5 =kx+5<iw&&ky-1>=0?src2[idx-rowlen  +20]<<8:0,
				cTR6 =kx+6<iw&&ky-1>=0?src2[idx-rowlen  +24]<<8:0,
				cTR7 =kx+7<iw&&ky-1>=0?src2[idx-rowlen  +28]<<8:0,

				cL6  =kx-6>=0         ?src2[idx         -24]<<8:0,
				cL4  =kx-4>=0         ?src2[idx         -16]<<8:0,
				cL   =kx-1>=0         ?src2[idx         - 4]<<8:0;
			int eT, eL, eTL, eTR, eT_L;
			int weights[PW2_NPRED];/*fixed 23.8 bit*/
			/*
			|				NNN
			|		NNWW		NN		NNEE
			|			NW	N	NE	NEE
			|	WWW	WW	W	?
			|
			|           T3
			|   T2L2    T2    T2R2
			|        TL T  TR TR2
			|L3 L2   L  X
			*/

			/*
			w0   w1   w2   w3
			p3Ca p3Cb p3Cc p3Cd p3Ce
			p1C  p2c
			*/
			
			for(k=0;k<PW2_NPRED;++k)
			{
				/*eNW + eN + eNE*/
				int e2=
					 (ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0)
					+(ky-1>=0?pred_errors[k][prevrow+kx]:0)
					+(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0);
				weights[k]=(params[k]<<8)/(e2+1);
			}

			/*
			TL T TR
			L  X
			*/
			eT=ky-1>=0?hireserror[prevrow+kx]:0;
			eL=kx-1>=0?hireserror[currrow+kx-1]:0;
			eTL=ky-1>=0&&kx-1>=0?hireserror[prevrow+kx-1]:0;
			eTR=ky-1>=0&&kx+1<iw?hireserror[prevrow+kx+1]:0;
			eT_L=eT+eL;
			{
				int64_t pred, sum;
				int vmin, vmax, curr;
				int predictions[PW2_NPRED];/*fixed 23.8 bit*/
				int j;

				j=0;
				/*from jxl*/
				predictions[j++]=cT-((eTL*params[PW2_NPRED]+eT*params[PW2_NPRED+1]+eTR*params[PW2_NPRED+2]+(cT2-cT)*params[PW2_NPRED+3]+(cTL-cL)*params[PW2_NPRED+4])>>8);/*k13: 1.998458 how many optimizations?*/
				/*k13: 1.737220, 1.837736, 1.842377, 1.847278, 1.865093, 1.861586, 1.866601, 1.872176, 1.878888, 1.883146, 1.883862, 1.883863, <*/

				predictions[j++]=cL-((eT_L+eTL)*params[PW2_NPRED+5]>>8);
				/*k13: 1.737220, 1.932766, 1.960469, 1.966698, 1.970106, 1.971381, 1.971891, 1.972187, 1.972368, 1.972492, 1.972565, 1.972580, 1.972594, 1.972602, <*/

				predictions[j++]=cT-((eT_L+eTR)*params[PW2_NPRED+6]>>8);
				/*k13: 1.737220, 1.934911, 1.951553, 1.962075, 1.973068, 1.979630, 1.983403, 1.985205, 1.986109, 1.986488, 1.986609, 1.986677, 1.986688, 1.986689, <*/

				predictions[j++]=cL+cTR-cT;
				/*k13: 1.909206, 1.918954, 1.946385, 1.963439, 1.970334, 1.981971, 1.983898, 1.984527, 1.985102, 1.985773, 1.986008, 1.986331, 1.986678, 1.986722, 1.986756			...1.998458*/


				predictions[j++]=cL-(eL*params[PW2_NPRED+7]>>8);
				predictions[j++]=cT-(eT*params[PW2_NPRED+8]>>8);
				/*k13: 1.909206, 1.922112, 1.931268, 1.954690, 1.964123, 1.978742, 1.981578, 1.984138, 1.985535, 1.986711, 1.987659, 1.988190, 1.988474, 1.988532, 1.988550*/
				predictions[j++]=cTL-(eTL*params[PW2_NPRED+9]>>8);
				predictions[j++]=cTR-(eTR*params[PW2_NPRED+10]>>8);
				/*k13: 1.909206, 1.921490, 1.932631, 1.949766, 1.950930, 1.951645, 1.951977, 1.960758, 1.967595, 1.969669, 1.972408, 1.973050, 1.973506, 1.974268, 1.975184			...1.977183*/


				/*paq8px by Matt Mahoney*/

				/*k13: 1.737220, 1.958513, 1.973267, 1.979685, 1.983374, 1.985860, 1.987622, 1.989731, 1.991147, 1.992018, 1.992707, 1.993444, 1.994374, 1.995238, 1.996056,   1.996876, 1.997423, 1.997708, 1.997946, 1.998162, 1.998320, 1.998364, 1.998611, 1.998815, 1.998948, 1.999125, 1.999207, 1.999222, 1.999229, 1.999235, 1.999241, 1.999242, 1.999247, 1.999248, 1.999250, 1.999251, <*/

				predictions[j++]=clamp4(cL+cT-cTL, cL, cTL, cT, cTR);/*0*/
				/*clip(cL+cT-cTL),*/ /*1*/
				predictions[j++]=clamp4(cL+cTR-cT, cL, cTL, cT, cTR);/*2*/
				/*clip(cL+cTR-cT),*/ /*3*/
				predictions[j++]=clamp4(cT+cTL-cT2L, cL, cTL, cT, cTR);/*4*/
				/*clip(cT+cTL-cT2L),*/ /*5*/
				predictions[j++]=clamp4(cT+cTR-cT2R, cL, cT, cTR, cTR2);/*6*/
				/*clip(cT+cTR-cT2R),*/ /*7*/
				predictions[j++]=(cL+cTR2)>>1;/*8*/
				/*clip(cT*3-cT2*3+cT3),*/ /*9*/
				/*clip(cL*3-cL2*3+cL3),*/ /*10*/
				/*(cL+clip(cTR*3-cT2R*3+cT3R))>>1,*/ /*11*/
				/*(cL+clip(cTR2*3-cT2R3*3+cT3R4))>>1,*/ /*12*/
				/*clip(cT2+cT4-cT6),*/ /*13*/
				/*clip(cL2+cL4-cL6),*/ /*14*/
				/*clip((cT5-6*cT4+15*cT3-20*cT2+15*cT+clamp4(cL*2-cTL2, cL, cTL, cT, cT2))/6),*/ /*15*/
				/*clip((-3*cL2+8*cL+clamp4(3*cTR2-3*cT2R2+cT3R2, cTR, cTR2, cTR3, cTR4))/6),*/ /*16*/
				/*clip(cT2+cTL-cT3L),*/ /*17*/
				/*clip(cT2+cTR-cT3R),*/ /*18*/
				/*clip((cL*2+cTL) - (cL2*2+cTL2) + cL3),*/ /*19*/
				/*clip(3*(cTL+cTL2)/2-cT2L3*3+(cT3L4+cT3L5)/2),*/ /*20*/
				/*clip(cTR2+cTR-cT2R3),*/ /*21*/
				/*clip(cTL2+cL2-cL4),*/ /*22*/
				/*clip(((cL+cTL)*3-cTL2*6+cTL3+cT2L3)/2),*/ /*23*/
				/*clip((cTR*2+cTR2) - (cT2R2+cT3R2*2) + cT4R3);*/ /*24*/
				predictions[j++]=cT6;/*25*/
				predictions[j++]=(cTR4+cTR6)>>1;/*26*/
				predictions[j++]=(cL4+cL6)>>1;/*27*/
				predictions[j++]=(cL+cT+cTR5+cTR7)>>2;/*28*/
				/*clip(cTR3+cL-cTR2),*/ /*29*/
				/*clip(4*cT3-3*cT4),*/ /*30*/
				/*clip(cT+cT2-cT3),*/ /*31*/
				/*clip(cL+cL2-cL3),*/ /*32*/
				/*clip(cL+cTR2-cTR),*/ /*33*/
				/*clip(cL2+cTR2-cT),*/ /*34*/
				/*(clip(cL*2-cTL)+clip(cL*2-cTL2)+cT+cTR)>>2,*/ /*35*/
				predictions[j++]=clamp4(cT*2-cT2, cL, cT, cTR, cTR2);/*36*/
				predictions[j++]=(cT+cT3)>>1;/*37*/
				/*clip(cT2+cL-cT2L),*/ /*38*/
				/*clip(cTR2+cT-cT2R2),*/ /*39*/
				/*clip((4*cL3-15*cL2+20*cL+clip(cTR2*2-cT2R2))/10),*/ /*40*/
				/*clip((cT3R3-4*cT2R2+6*cTR+clip(cL*3-cTL*3+cT2L))/4),*/ /*41*/
				/*clip((cT*2+cTR) - (cT2+2*cT2R) + cT3R),*/ /*42*/
				/*clip((cTL*2+cT2L) - (cT2L2+cT3L2*2) + cT4L3),*/ /*43*/
				/*clip(cT2L2+cL-cT2L3),*/ /*44*/
				/*clip((-cT4+5*cT3-10*cT2+10*cT+clip(cL*4-cTL2*6+cT2L3*4-cT3L4))/5),*/ /*45*/
				/*clip(cTR2+clip(cTR3*2-cT2R4-cTR4)),*/ /*46*/
				/*clip(cTL+cL-cTL2),*/ /*47*/
				/*clip((cT*2+cTL) - (cT2+2*cT2L) + cT3L),*/ /*48*/
				/*clip(cT2+clip(cTR2*2-cT2R3) - cT2R),*/ /*49*/
				/*clip((-cL4+5*cL3-10*cL2+10*cL+clip(cTR*2-cT2R))/5),*/ /*50*/
				/*clip((-cL5+4*cL4-5*cL3+5*cL+clip(cTR*2-cT2R))>>2),*/ /*51*/
				/*clip((cL3-4*cL2+6*cL+clip(cTR*3-cT2R*3+cT3R))>>2),*/ /*52*/
				/*clip((-cT2R2+3*cTR+clip(4*cL-6*cTL+4*cT2L-cT3L))/3),*/ /*53*/
				predictions[j++]=((cL+cT)*3-cTL*2)>>2;/*54*/

				sum=0;
				for(k=0;k<PW2_NPRED;++k)
					sum+=weights[k];
				if(sum)
				{
					pred=(sum>>1)-1;
					for(k=0;k<PW2_NPRED;++k)
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
				for(k=0;k<PW2_NPRED;++k)
				{
					int e=abs(curr-predictions[k]);
					pw2_errors[k]+=e;
					pred_errors[k][currrow+kx]=e;
					if(kx+1<iw)/*add current error to eNE, such that eN = eN0 + eW*/
						pred_errors[k][prevrow+kx+1]+=e;
				}
			}
		}
	}
	for(k=0;k<PW2_NPRED;++k)
		pw2_errors[k]/=iw*ih*256;
}
static double pred_w2_calcloss(const int *src, int iw, int ih, int depth, int src_depth, int kc, int16_t *params, int *temp, int *dst, int *hist)
{
	double entropy, invCR, csize;
	intptr_t res;

	res=(intptr_t)iw*ih;
	pred_w2_prealloc(src, iw, ih, depth, kc, params, 1, 1, dst, temp);
	calc_histogram(dst, iw, ih, kc, 0, iw, 0, ih, depth, hist, 0);
	entropy=calc_entropy(hist, 1<<depth, (int)res);
	invCR=entropy/src_depth;
	csize=res*invCR;
	return csize;
}
void pred_w2_opt_v2(Image *src, int16_t *params, int loud)
{
	const int steps[]={128, 64, 32, 16, 8, 4, 2, 1};
	intptr_t res;
	int *buf3, *temp, maxdepth, nlevels, *hist, kc, ks, idx;
	char title0[256];

	(void)loud;
	res=(intptr_t)src->iw*src->ih;
	buf3=(int*)malloc((size_t)res*sizeof(int[4]));
	temp=(int*)malloc((size_t)src->iw*(PW2_NPRED+1)*2*sizeof(int));
	maxdepth=calc_maxdepth(src, 0);
	nlevels=1<<maxdepth;
	hist=(int*)malloc(nlevels*sizeof(int));
	if(!buf3||!temp||!hist)
	{
		CRASH("Alloc error");
		return;
	}
	get_window_title(title0, 256);
	for(kc=0;kc<3;++kc)
	{
		int16_t *param=params+kc*PW2_NPARAM;
		double csize0=pred_w2_calcloss(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);

		for(ks=0;ks<(int)_countof(steps);++ks)
		{
			int step=steps[ks];
			double bestcsize=csize0;
			int bestidx=0, beststep=0;

			for(idx=0;idx<PW2_NPARAM;++idx)
			{
				double csize;
				int16_t prev;

				prev=param[idx];
				param[idx]+=(int16_t)step;
				csize=pred_w2_calcloss(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=step;

				prev=param[idx];
				param[idx]-=(int16_t)step;
				if(idx<4&&param[idx]<1)
					param[idx]=1;
				csize=pred_w2_calcloss(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=-step;

				if(loud_transforms)
					set_window_title("Ch%d csize %f [%d/%d %d/%d]...", kc, csize0, (int)(kc*_countof(steps)+ks+1), (int)(_countof(steps)*3), idx+1, PW2_NPARAM);
			}
			if(csize0>bestcsize)
			{
				csize0=bestcsize;

				param[bestidx]+=(int16_t)beststep;
				if(bestidx<4&&param[bestidx]<1)
					param[bestidx]=1;
			}
		}
	}
	free(hist);
	free(temp);
	free(buf3);
	if(loud_transforms)
		set_window_title("%s", title0);
}
void pred_w2_apply(Image *src, int fwd, int enable_ma, int16_t *params)
{
	intptr_t k, res;
	int *temp, *buf2;

	res=(intptr_t)src->iw*src->ih;
	temp=(int*)malloc((size_t)src->iw*(PW2_NPRED+1)*2*sizeof(int));
	buf2=(int*)malloc(res*sizeof(int[4]));
	if(!temp||!buf2)
	{
		CRASH("Alloc error");
		return;
	}
	pred_w2_prealloc(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[0], 0, params             , fwd, enable_ma, buf2, temp);
	pred_w2_prealloc(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[1], 1, params+PW2_NPARAM  , fwd, enable_ma, buf2, temp);
	pred_w2_prealloc(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[2], 2, params+PW2_NPARAM*2, fwd, enable_ma, buf2, temp);
	if(src->depth[3])
		pred_w2_prealloc(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[3], 3, params+PW2_NPARAM, fwd, enable_ma, buf2, temp);
	for(k=0;k<res;++k)
	{
		((int32_t*)HDR_DATA(src))[k<<2|0]=buf2[k<<2|0];
		((int32_t*)HDR_DATA(src))[k<<2|1]=buf2[k<<2|1];
		((int32_t*)HDR_DATA(src))[k<<2|2]=buf2[k<<2|2];
		if(src->depth[3])
			((int32_t*)HDR_DATA(src))[k<<2|3]=buf2[k<<2|3];
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


int16_t jxlparams_i16[33]=/*signed fixed 7.8 bit*/
{
	0x0DB8,  0x0E22,  0x181F,  0x0BF3, -0x005C, -0x005B,  0x00DF,  0x0051,  0x00BD,  0x005C, -0x0102,/*Y	*/
	0x064C,  0x0F31,  0x1040,  0x0BF8, -0x0007, -0x000D, -0x0085, -0x0063, -0x00A2, -0x0017,  0x00F2,/*Cb	*/
	0x0B37,  0x110B,  0x121B,  0x0BFC, -0x0001,  0x000E, -0x0188, -0x00E7, -0x00BB, -0x004A,  0x00BA,/*Cr	*/

#if 0
	0x0DD3, 0x1002, 0x11F9, 0x0BEC, 0x0002,-0x0029,-0x010E,-0x00E8,-0x00A2,-0x0089, 0x0064,
	0x0FAF, 0x0C9E, 0x139F, 0x0E13,-0x009D,-0x005C, 0x0041, 0x0013, 0x00DD,-0x0004,-0x00C3,
	0x09A8, 0x0F00, 0x1191, 0x0B9E,-0x0026,-0x0027,-0x0086,-0x0084,-0x005F,-0x0038, 0x00B2,

	0x1000, 0x1000, 0x1000, 0x1000, 0, 0, 0, 0, 0, 0, 0,
	0x1000, 0x1000, 0x1000, 0x1000, 0, 0, 0, 0, 0, 0, 0,
	0x1000, 0x1000, 0x1000, 0x1000, 0, 0, 0, 0, 0, 0, 0,

	0x0C7B, 0x0B5C, 0x0A14, 0x0B33,  -0x0148,-0x0029, 0, 0, 0, 0, 0,
	0x0A14, 0x0829, 0x1548, 0x0CA4,   0x047B, 0x0052, 0, 0, 0, 0, 0,
	0x0B33, 0x0C29, 0x0DC3, 0x119A,  -0x0148,-0x00F6, 0, 0, 0, 0, 0,

	0x0C7B, 0x0B5C, 0x0A14, 0x0B33,-0x0148,-0x0029, 0x0971, 0x01EC,-0x01C3, 0x047B, 0x0AB8,
	0x0A14, 0x0829, 0x1548, 0x0CA4, 0x047B, 0x0052,-0x011F, 0x0000, 0x0029, 0x063D, 0x0266,
	0x0B33, 0x0C29, 0x0DC3, 0x119A,-0x0148,-0x00F6, 0x0614, 0x00A4,-0x007B, 0x019A, 0x0E8F,

	0x63, 0x5A, 0x50, 0x59,    -0x0A, -0x01,  0x4B, 0x0F, -0x0E, -0x24, 0x55,
	0x50, 0x41, 0xA9, 0x64,     0x24,  0x03, -0x09,    0,  0x01,  0x32, 0x13,
	0x59, 0x61, 0x6D, 0x8C,    -0x0A, -0x08,  0x30, 0x05, -0x04,  0x0D, 0x74,
#endif
};
static void pred_jxl_prealloc(const int *src, int iw, int ih, int depth, int kc, int16_t *params, int fwd, int enable_ma, int *dst, int *temp_w10)
{
	const int *pixels;
	int errorbuflen, rowlen, *hireserror, *pred_errors[4], idx, k, ky, kx;

	errorbuflen=iw<<1;
	rowlen=iw<<2;
	hireserror=temp_w10;
	pred_errors[0]=temp_w10+errorbuflen*1;
	pred_errors[1]=temp_w10+errorbuflen*2;
	pred_errors[2]=temp_w10+errorbuflen*3;
	pred_errors[3]=temp_w10+errorbuflen*4;
	idx=kc;
	pixels=fwd?src:dst;
	memset(temp_w10, 0, iw*sizeof(int[10]));
	for(ky=0;ky<ih;++ky)
	{
		int currrow=ky&1?0:iw, prevrow=ky&1?iw:0;
		for(kx=0;kx<iw;++kx, idx+=4)
		{
			int
				ctt      =         ky-2>=0?pixels[idx-rowlen*2]:0,
				ctopleft =kx-1>=0&&ky-1>=0?pixels[idx-rowlen-4]:0,
				ctop     =kx  <iw&&ky-1>=0?pixels[idx-rowlen  ]:0,
				ctopright=kx+1<iw&&ky-1>=0?pixels[idx-rowlen+4]:0,
				cleft    =kx-1>=0         ?pixels[idx       -4]:0;
			int weights[4];/*fixed 23.8 bit*/
			int pred, curr, pred_final, e;

			/*
			w0   w1   w2   w3
			p3Ca p3Cb p3Cc p3Cd p3Ce
			p1C  p2c
			*/
			for(k=0;k<4;++k)
			{
				/*
				|	eNW	eN	eNE
				|eWW	eW	?
				*/
				int e2=
					+(ky-1>=0?pred_errors[k][prevrow+kx]:0)			/*eN + eW	*/
					+(ky-1>=0&&kx+1<iw?pred_errors[k][prevrow+kx+1]:0)	/*eNE		*/
					+(ky-1>=0&&kx-1>=0?pred_errors[k][prevrow+kx-1]:0);	/*eNW + eWW	*/
				weights[k]=(params[k]<<8)/(e2+1);
			}
			{
				int
					etop=ky-1>=0?hireserror[prevrow+kx]:0,
					eleft=kx-1>=0?hireserror[currrow+kx-1]:0,
					etopleft=ky-1>=0&&kx-1>=0?hireserror[prevrow+kx-1]:0,
					etopright=ky-1>=0&&kx+1<iw?hireserror[prevrow+kx+1]:0,
					etopplusleft=etop+eleft;
				int sum;
				int64_t predictions[4];/*fixed 23.8 bit*/

				predictions[0]=(cleft+ctopright-ctop)<<8;
				predictions[1]=(ctop<<8)-((etopplusleft+etopright)*params[4]>>8);
				predictions[2]=(cleft<<8)-((etopplusleft+etopleft)*params[5]>>8);
				predictions[3]=(ctop<<8)-(((etopleft*params[6]+etop*params[7]+etopright*params[8])>>8)+(ctt-ctop)*params[9]+(ctopleft-cleft)*params[10]);
				sum=weights[0]+weights[1]+weights[2]+weights[3];

				if(sum)
					pred=(int)((predictions[0]*weights[0]+predictions[1]*weights[1]+predictions[2]*weights[2]+predictions[3]*weights[3]+(sum>>1)-1)/sum);
				else
					pred=(int)predictions[0];

				CLAMP3_32(pred, pred, ctop<<8, cleft<<8, ctopright<<8);
				pred_final=(pred+127)>>8;
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

				hireserror[currrow+kx]=curr-pred;
				for(k=0;k<4;++k)
				{
					e=abs(curr-(int)predictions[k]);
					pred_errors[k][currrow+kx]=e;
					if(kx+1<iw)
						pred_errors[k][prevrow+kx+1]+=e;/*eNE+=ecurr*/
				}
			}
		}
	}
}
static double pred_jxl_calcloss(const int *src, int iw, int ih, int depth, int src_depth, int kc, int16_t *params, int *temp, int *dst, int *hist)
{
	double entropy, invCR, csize;
	intptr_t res;

	res=(intptr_t)iw*ih;
	pred_jxl_prealloc(src, iw, ih, depth, kc, params, 1, 1, dst, temp);
	calc_histogram(dst, iw, ih, kc, 0, iw, 0, ih, depth, hist, 0);
	entropy=calc_entropy(hist, 1<<depth, (int)res);
	invCR=entropy/src_depth;
	csize=res*invCR;
	return csize;
}
void pred_jxl_opt_v2(Image *src, int16_t *params, int loud)
{
	const int steps[]={256, 128, 64, 32, 16, 8, 4, 2, 1};
	intptr_t res;
	int *buf3, *temp, maxdepth, nlevels, *hist, kc, ks, idx;

	(void)loud;
	res=(intptr_t)src->iw*src->ih;
	buf3=(int*)malloc((size_t)res*sizeof(int[4]));
	temp=(int*)malloc((size_t)src->iw*10*sizeof(int));
	maxdepth=calc_maxdepth(src, 0);
	nlevels=1<<maxdepth;
	hist=(int*)malloc(nlevels*sizeof(int));
	if(!buf3||!temp||!hist)
	{
		CRASH("Alloc error");
		return;
	}
	for(kc=0;kc<3;++kc)
	{
		int16_t *param=params+kc*11;
		double csize0=pred_jxl_calcloss(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);

		for(ks=0;ks<(int)_countof(steps);++ks)
		{
			int step=steps[ks];
			double bestcsize=csize0;
			int bestidx=0, beststep=0;

			for(idx=0;idx<11;++idx)
			{
				double csize;
				int16_t prev;

				prev=param[idx];
				param[idx]+=(int16_t)step;
				csize=pred_jxl_calcloss(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=step;

				prev=param[idx];
				param[idx]-=(int16_t)step;
				if(idx<4&&param[idx]<1)
					param[idx]=1;
				csize=pred_jxl_calcloss(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[kc], src->src_depth[kc], kc, param, temp, buf3, hist);
				param[idx]=prev;
				if(bestcsize>csize)
					bestcsize=csize, bestidx=idx, beststep=-step;
			}
			if(csize0>bestcsize)
			{
				csize0=bestcsize;

				param[bestidx]+=(int16_t)beststep;
				if(bestidx<4&&param[bestidx]<1)
					param[bestidx]=1;
			}
		}
	}
	free(hist);
	free(temp);
	free(buf3);
}
void pred_jxl_apply(Image *src, int fwd, int enable_ma, int16_t *params)
{
	int k;
	intptr_t res=(intptr_t)src->iw*src->ih;
	int *temp=(int*)malloc((size_t)src->iw*10*sizeof(int));
	int *buf2=(int*)malloc(res*sizeof(int[4]));

	if(!temp||!buf2)
	{
		CRASH("Alloc error");
		return;
	}
	pred_jxl_prealloc(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[0], 0, params   , fwd, enable_ma, buf2, temp);
	pred_jxl_prealloc(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[1], 1, params+11, fwd, enable_ma, buf2, temp);
	pred_jxl_prealloc(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[2], 2, params+22, fwd, enable_ma, buf2, temp);
	if(src->depth[3])
		pred_jxl_prealloc(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[3], 3, params+11, fwd, enable_ma, buf2, temp);

	for(k=0;k<res;++k)
	{
		((int32_t*)HDR_DATA(src))[k<<2  ]=buf2[k<<2  ];
		((int32_t*)HDR_DATA(src))[k<<2|1]=buf2[k<<2|1];
		((int32_t*)HDR_DATA(src))[k<<2|2]=buf2[k<<2|2];
		if(src->depth[3])
			((int32_t*)HDR_DATA(src))[k<<2|3]=buf2[k<<2|3];
	}

	free(temp);
	free(buf2);
}

void pred_jmj_apply(Image *src, int fwd, int enable_ma)
{
	int k;
	intptr_t res=(intptr_t)src->iw*src->ih;
	int *temp=(int*)malloc((size_t)src->iw*(PW2_NPRED+1)*2*sizeof(int));
	int *buf2=(int*)malloc(res*sizeof(int[4]));

	if(!temp||!buf2)
	{
		CRASH("Alloc error");
		return;
	}
	pred_jxl_prealloc(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[0], 0, jxlparams_i16        , fwd, enable_ma, buf2, temp);
	pred_w2_prealloc (((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[1], 1, pw2_params+PW2_NPARAM, fwd, enable_ma, buf2, temp);
	pred_jxl_prealloc(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[2], 2, jxlparams_i16+22     , fwd, enable_ma, buf2, temp);
	if(src->depth[3])
		pred_jxl_prealloc(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth[3], 3, jxlparams_i16+22, fwd, enable_ma, buf2, temp);

	for(k=0;k<res;++k)
	{
		((int32_t*)HDR_DATA(src))[k<<2  ]=buf2[k<<2  ];
		((int32_t*)HDR_DATA(src))[k<<2|1]=buf2[k<<2|1];
		((int32_t*)HDR_DATA(src))[k<<2|2]=buf2[k<<2|2];
		if(src->depth[3])
			((int32_t*)HDR_DATA(src))[k<<2|3]=buf2[k<<2|3];
	}

	free(temp);
	free(buf2);
}


/*CUSTOM3*/

	#define C3_3D

#define C3_OPT_NCOMP (C3_NPARAMS>>5)
Custom3Params c3_params={0};
static int fast_dot(const int16_t *a, const int16_t *b, int count)
{
	/*
	https://stackoverflow.com/questions/62041400/inner-product-of-two-16bit-integer-vectors-with-avx2-in-c
	*/
	int k;

	__m256i sum=_mm256_setzero_si256();
	for(k=0;k<count-31;k+=16)
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
static int custom3_loadnb(const int *pixels, const int *errors, int iw, int ih, int kc, int kx, int ky, int16_t *nb)
{
	int idx, ky2, kx2;

	idx=-1;
	for(ky2=-C3_REACH;ky2<0;++ky2)
	{
		for(kx2=-C3_REACH;kx2<=C3_REACH;++kx2)
		{
			if((uint32_t)(kx+kx2)<(uint32_t)iw&&(uint32_t)(ky+ky2)<(uint32_t)ih)
			{
				int idx2=(iw*(ky+ky2)+kx+kx2)<<2|kc;
				nb[++idx]=(int16_t)pixels[idx2];
				nb[++idx]=(int16_t)errors[idx2];
			}
			else
			{
				nb[++idx]=0;
				nb[++idx]=0;
			}
		}
	}
	for(kx2=-C3_REACH;kx2<0;++kx2)
	{
		if((uint32_t)(kx+kx2)<(uint32_t)iw)
		{
			int idx2=(iw*ky+kx+kx2)<<2|kc;
			nb[++idx]=(int16_t)pixels[idx2];
			nb[++idx]=(int16_t)errors[idx2];
		}
		else
		{
			nb[++idx]=0;
			nb[++idx]=0;
		}
	}
	return ++idx;
}
static void custom3_prealloc(const int *src, int iw, int ih, const int8_t *depths, int fwd, int enable_ma, Custom3Params const *params, int *dst)
{
	const int *pixels=fwd?src:dst, *errors=fwd?dst:src;
	int pred, idx, j, nlevels[4], ky, kx, kc, kn, kdst;
	const int16_t *coeffs[9];
	int16_t nb[3][C3_NNB+2]={0};

	j=0;
	coeffs[j++]=params->c00; coeffs[j++]=params->c01; coeffs[j++]=params->c02;
	coeffs[j++]=params->c10; coeffs[j++]=params->c11; coeffs[j++]=params->c12;
	coeffs[j++]=params->c20; coeffs[j++]=params->c21; coeffs[j++]=params->c22;
	nlevels[0]=1<<depths[0];
	nlevels[1]=1<<depths[1];
	nlevels[2]=1<<depths[2];
	nlevels[3]=1<<depths[3];
	for(ky=0;ky<ih;++ky)
	{
		for(kx=0;kx<iw;++kx)
		{
#if 1
			int clampers[6];
			int val;
			int count[3], idx2;

			j=0;
			clampers[j++]=(nlevels[0]>>1)-1; clampers[j++]=-(nlevels[0]>>1);
			clampers[j++]=(nlevels[1]>>1)-1; clampers[j++]=-(nlevels[1]>>1);
			clampers[j++]=(nlevels[2]>>1)-1; clampers[j++]=-(nlevels[2]>>1);
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
			for(kc=0;kc<3;++kc)
				count[kc]=custom3_loadnb(pixels, errors, iw, ih, kc, kx, ky, nb[kc]);
#ifdef C3_3D
			for(kc=1;kc<3;++kc)
			{
				for(kn=0;kn<count[0];++kn)
					nb[kc][kn]-=nb[0][kn];
			}
#endif
			idx=(iw*ky+kx)<<2;
			idx2=0;
			
			for(kdst=0;kdst<3;++kdst)
			{
				pred=0;
				for(kc=0;kc<3;++kc)
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

				nb[kdst][C3_NNB  ]=(int16_t)pixels[idx];
				nb[kdst][C3_NNB+1]=(int16_t)errors[idx];
				count[kdst]+=2;
				++idx;
				idx2+=3;
			}
		}
	}
}

void custom3_apply(Image *src, int fwd, int enable_ma, Custom3Params const *params)
{
	intptr_t res=(intptr_t)src->iw*src->ih;
	int *temp=(int*)malloc(res*sizeof(int[4]));
	if(!temp)
	{
		CRASH("Alloc error");
		return;
	}
	memcpy(temp, ((int32_t*)HDR_DATA(src)), res*sizeof(int[4]));/*copy alpha*/

	custom3_prealloc(((int32_t*)HDR_DATA(src)), src->iw, src->ih, src->depth, fwd, enable_ma, params, temp);

	memcpy(((int32_t*)HDR_DATA(src)), temp, res*sizeof(int[4]));
	free(temp);
}

typedef struct Custom3OptInfoStruct
{
	double invCR[4];/*loss == invCR[3]=(invCR[0]+invCR[1]+invCR[2])/3*/
	int16_t params[C3_NPARAMS];
} Custom3OptInfo;
static void custom3_calcloss(const int *src, int iw, int ih, const int8_t *depths, const int8_t *src_depths, Custom3OptInfo *info, int *temp, int *hist)
{
	double entropy;
	intptr_t res;
	int kc;

	res=(intptr_t)iw*ih;
	custom3_prealloc(src, iw, ih, depths, 1, 1, (Custom3Params*)info->params, temp);
	
	info->invCR[3]=0;
	for(kc=0;kc<3;++kc)
	{
		calc_histogram(temp, iw, ih, kc, 0, iw, 0, ih, depths[kc], hist, 0);
		{
			entropy=calc_entropy(hist, 1<<depths[kc], (int)res);
			info->invCR[kc]=entropy/src_depths[kc];
			info->invCR[3]+=info->invCR[kc];
		}
	}
	info->invCR[3]/=3;
}
void custom3_opt(Image const *src, Custom3Params *srcparams, int niter, int maskbits, int loud, double *loss)
{
	static int call_idx=0;
	int it, watchdog, k;
	intptr_t res=(intptr_t)src->iw*src->ih;
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
		CRASH("Alloc error");
		return;
	}
	memset(temp, 0xFF, res*sizeof(int[4]));/*set alpha for preview*/
	memcpy(info.params, srcparams, sizeof(info.params));
#define CALC_LOSS() custom3_calcloss(((const int32_t*)HDR_DATC(src)), src->iw, src->ih, src->depth, src->src_depth, &info, temp, hist)
	if(loud)
		srand((uint32_t)__rdtsc());

	CALC_LOSS();
	memcpy(invCR, info.invCR, sizeof(invCR));
	loss0=info.invCR[3];
	if(!niter)
		niter=C3_NPARAMS*10;
	for(it=0, watchdog=0;it<niter;++it)
	{
		int idx[C3_OPT_NCOMP]={0}, inc[C3_OPT_NCOMP]={0}, stuck=0;
		if(loud_transforms&&loud)
			set_window_title("%d %4d/%4d,%d/%d: %.4f%% RGB %.4f %.4f %.4f%%",
				call_idx, it+1, niter, watchdog, shakethreshold,
				100.*info.invCR[3],
				100.*info.invCR[0],
				100.*info.invCR[1],
				100.*info.invCR[2]
			);
		if(watchdog>=shakethreshold)/*bump if stuck*/
		{
			memcpy(info.params, srcparams, sizeof(info.params));
			for(k=0;k<C3_NPARAMS;++k)
				info.params[k]+=((rand()&1)<<1)-1;
			watchdog=0;
			stuck=1;
		}
		else
		{
			for(k=0;k<C3_OPT_NCOMP;++k)
			{
				idx[k]=rand()%C3_NPARAMS;
				while(!(inc[k]=(rand()&((1<<maskbits)-1))-(1<<(maskbits-1))));/*reject zero delta*/
				info.params[idx[k]]+=(int16_t)inc[k];
			}
		}

		CALC_LOSS();
		
		if(info.invCR[3]>invCR[3])/*revert if worse*/
		{
			if(stuck)/*a bad branch may surpass the local minimum*/
				memcpy(invCR, info.invCR, sizeof(info.invCR));
			else
			{
				memcpy(info.invCR, invCR, sizeof(info.invCR));
				for(k=0;k<C3_OPT_NCOMP;++k)
					info.params[idx[k]]-=(int16_t)inc[k];
			}
			++watchdog;
		}
		else/*save if better*/
		{
			if(loss0>info.invCR[3])
			{
				memcpy(srcparams, info.params, sizeof(info.params));
				loss0=info.invCR[3];
				--it;/*again*/
			}
			memcpy(invCR, info.invCR, sizeof(invCR));
			watchdog=0;
		}

		/*preview*/
#if 1
		if(loud)
		{
			ch_entropy[0]=(float)(info.invCR[0]*src->src_depth[0]);
			ch_entropy[1]=(float)(info.invCR[1]*src->src_depth[1]);
			ch_entropy[2]=(float)(info.invCR[2]*src->src_depth[2]);
			io_render();
		}
#endif
	}
#undef  CALC_LOSS
	if(loss)
		memcpy(loss, invCR, sizeof(invCR));
	
	free(temp);
	free(hist);
}


/*
CALIC - optimized for 8-bit
https://github.com/play-co/gcif
*/
void pred_calic(Image *src, int fwd, int enable_ma)
{
	int amin[4], amax[4], ky, kx, kc, idx, bufsize, pred, *arrN, *arrS;
	int16_t *pixels, *rows[4];

	amin[0]=-(1<<src->depth[0]>>1);
	amin[1]=-(1<<src->depth[1]>>1);
	amin[2]=-(1<<src->depth[2]>>1);
	amin[3]=-(1<<src->depth[3]>>1);
	amax[0]=(1<<src->depth[0]>>1)-1;
	amax[1]=(1<<src->depth[1]>>1)-1;
	amax[2]=(1<<src->depth[2]>>1)-1;
	amax[3]=(1<<src->depth[3]>>1)-1;
	bufsize=(src->iw+8*2)*(int)sizeof(int16_t[4*1*2]);/*4 padded rows * 1 channel * {pixels, errors}*/
	pixels=(int16_t*)malloc(bufsize);
	if(!pixels)
	{
		CRASH("Alloc error");
		return;
	}
	arrN=(int*)malloc(2048*sizeof(int));/*count*/
	arrS=(int*)malloc(2048*sizeof(int));/*sum*/
	if(!arrN||!arrS)
	{
		CRASH("Alloc error");
		return;
	}
	for(kc=0;kc<3;++kc)/*process each channel separately*/
	{
		int depth=src->depth[kc], nlevels=1<<depth, half=nlevels>>1, sh=depth-8;

		memset(arrN, 0, 2048*sizeof(int));
		memset(arrS, 0, 2048*sizeof(int));
		memset(pixels, 0, bufsize);
		for(ky=0;ky<src->ih;++ky)
		{
			rows[0]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)0+4)%4)+8)*2;
			rows[1]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)1+4)%4)+8)*2;
			rows[2]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)2+4)%4)+8)*2;
			rows[3]=pixels+((src->iw+(intptr_t)16)*((ky-(intptr_t)3+4)%4)+8)*2;
			for(kx=0;kx<src->iw;++kx)
			{
				rows[0]+=2;
				rows[1]+=2;
				rows[2]+=2;
				rows[3]+=2;
				{
					int
						NN	=rows[2][0+0*2],
						NNE	=rows[2][0+1*2],
						NW	=rows[1][0-1*2],
						N	=rows[1][0+0*2],
						NE	=rows[1][0+1*2],
						WW	=rows[0][0-2*2],
						W	=rows[0][0-1*2],
						eN	=rows[1][1+0*2],
						eW	=rows[0][1-1*2];
					int energy, pattern, context, dx, dy, temp, curr, val, e;
				/*	int vmin, vmax;*/
				/*	int d45, d135;*/

					/*
					NNWW NNW NN NNE NNEE
					NWW  NW  N  NE
					WW   W   ?
					*/
					dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
					dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
#if 0
					d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
					d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
#endif

					/*
					'CALIC - A context-based adaptive lossless image codec'
					https://github.com/play-co/gcif/tree/master/refs/calic
					*/
#if 1
					temp=dy-dx;
					temp=SHIFT_RIGHT_SIGNED(temp, sh);
					if(temp>80)
						pred=W;
					else if(temp<-80)
						pred=N;
					else
					{
						pred=(W+N)/2+(NE-NW)/4;
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

					/*clamped grad*/
#if 0
					if(N<W)
						vmin=N, vmax=W;
					else
						vmin=W, vmax=N;
					pred=N+W-NW;
					pred=CLAMP(vmin, pred, vmax);
#endif

					/*'A context-based adaptive lossless/nearly-lossless coding scheme for continuous-tone images'*/
#if 0
					if(dy+dx>32)/*sharp edge*/
						pred=(dy*W+dx*N)/(dy+dx)+(NE-NW)/8;
					else if(dy-dx>12)/*horizontal edge*/
						pred=(2*W+N)/3+(NE-NW)/8;
					else if(dy-dx<-12)/*vertical edge*/
						pred=(W+2*N)/3+(NE-NW)/8;
					else/*shallow area*/
						pred=(W+N)/2+(NE-NW)/8;

					if(d45-d135>32)/*sharp 135-deg diagonal edge*/
						pred+=(NE-NW)/8;
					else if(d45-d135>16)/*135-deg diagonal edge*/
						pred+=(NE-NW)/16;
					else if(d45-d135<-32)/*sharp 45-deg diagonal edge*/
						pred-=(NE-NW)/8;
					else if(d45-d135<-16)/*45-deg diagonal edge*/
						pred-=(NE-NW)/16;
#endif

					energy=dx+dy+abs(eN)+abs(eW);
					energy=SHIFT_RIGHT_SIGNED(energy, sh);
					if(energy<42)/*quantization using binary search*/
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
						pred+=(arrS[context]+arrN[context]/2)/arrN[context];/*sum of encountered errors / number of occurrences*/
					pred=CLAMP(-half, pred, half-1);

					idx=(src->iw*ky+kx)<<2|kc;
					curr=((int32_t*)HDR_DATA(src))[idx];
					if(fwd)
					{
						val=(curr-(int)pred+g_dist/2)/g_dist;
						curr=g_dist*val+(int)pred;
					}
					else
					{
						val=g_dist*curr+(int)pred;
						curr=val;
						CLAMP2(val, amin[kc], amax[kc]);
					}
					((int32_t*)HDR_DATA(src))[idx]=val;
					e=curr-pred;
					rows[0][0]=curr;
					rows[0][1]=e;

					/*update*/
					++arrN[context];
					arrS[context]+=e;
					if(arrN[context]>255)
					{
						arrN[context]>>=1;
						arrS[context]>>=1;
					}
				}
			}
		}
	}
	free(pixels);
	free(arrN);
	free(arrS);

	(void)enable_ma;
}


#define G2_NPRED 21
#define G2_SSE_SIZE 1024
int16_t g2_weights[]=
{
	/*g	Y*/
	 0x00EA, 0x01C8, 0x00A2, 0x005E, 0x01F4, 0x0045, 0x0091, 0x0066, 0x003B, 0x0027,-0x0011, 0x001B, 0x00FF, 0x007E, 0x00D1, 0x00F3, 0x008F, 0x0130, 0x018E,-0x00AC,
	 0x0004,
	/*0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010,*/
	 0x010C, 0x0008,-0x007E, 0x00A2, 0x000E,-0x0069,-0x0073,-0x0125,-0x0092, 0x0000, 0x0078,

	/*b	Cb*/
	 0x0006, 0x003D, 0x0031, 0x002F, 0x003F, 0x0015, 0x0011, 0x0036, 0x002E,-0x0022, 0x0011, 0x0034,-0x0007, 0x0012,-0x0018, 0x0012, 0x002F, 0x0000, 0x0000, 0x001C,
	 0x0004,
	/*0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010,*/
	 0x00A2, 0x02E1, 0x00C9,-0x00E0,-0x0068,-0x004E,-0x013E,-0x0012, 0x0001, 0x0000,-0x0046,
	 
	/*r	Cr*/
	 0x003D, 0x0036, 0x0006, 0x007E, 0x0012, 0x0007, 0x0007, 0x0005, 0x001E, 0x0000, 0x0028, 0x0055,-0x0020, 0x0020, 0x0005, 0x0011, 0x0034, 0x0000, 0x0004, 0x003E,
	 0x0004,
	/*0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010, 0x0010,*/
	 -0x0100, 0x0001,-0x0086,-0x0041, 0x0051,-0x0080, 0x0004, 0x0002,-0x0003,-0x0003, 0x00D9,
};
typedef int64_t SseCell;
void pred_grad2(Image *src, int fwd, int enable_ma)
{
	int j, k, k2, kc, ky, kx;
	intptr_t res=(intptr_t)src->iw*src->ih;
	int *b2=(int*)malloc((size_t)res*sizeof(int[4]));
	int *perrors=(int*)malloc((size_t)src->iw*(G2_NPRED+1)*2*sizeof(int));
	SseCell *sse=(SseCell*)malloc(sizeof(SseCell[G2_SSE_SIZE]));
	const int *pixels=fwd?((int32_t*)HDR_DATA(src)):b2, *errors=fwd?b2:((int32_t*)HDR_DATA(src));
	if(!b2||!perrors||!sse)
	{
		CRASH("Alloc error");
		return;
	}
	memcpy(b2, ((int32_t*)HDR_DATA(src)), res*sizeof(int[4]));/*copy alpha*/
	for(kc=0;kc<4;++kc)
	{
		int depth;
		int16_t *params;
		int *hireserror;

		depth=src->depth[kc];
		if(!depth)
			continue;
		params=g2_weights+(_countof(g2_weights)/3)*kc;
		hireserror=perrors+src->iw*2*G2_NPRED;
		memset(perrors, 0, (intptr_t)2*src->iw*(G2_NPRED+1)*sizeof(int));

		memset(sse, 0, sizeof(SseCell[G2_SSE_SIZE]));
		for(ky=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx)
			{
#define LOAD(BUF, X, Y) (uint32_t)(kx+(X))<(uint32_t)src->iw&&(uint32_t)(ky+(Y))<(uint32_t)src->ih?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]<<8:0
				int
					NNNNNN  =LOAD(pixels,  0, -6),
					NNN     =LOAD(pixels,  0, -3),
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
					WW      =LOAD(pixels, -2,  0),
					W       =LOAD(pixels, -1,  0);
#undef  LOAD
#define LOAD(BUF, X, Y) (uint32_t)(kx+(X))<(uint32_t)src->iw&&(uint32_t)(ky+(Y))<(uint32_t)src->ih?BUF[(src->iw*(ky+(Y))+kx+(X))<<2|kc]:0
				int
					dNNWW  =LOAD(errors, -2, -2),
					dNNW   =LOAD(errors, -1, -2),
					dNN    =LOAD(errors,  0, -2),
					dNNE   =LOAD(errors,  1, -2),
					dNNEE  =LOAD(errors,  2, -2),
					dNWW   =LOAD(errors, -2, -1),
					dNW    =LOAD(errors, -1, -1),
					dN     =LOAD(errors,  0, -1),
					dNE    =LOAD(errors,  1, -1),
					dNEE   =LOAD(errors,  2, -1),
					dWW    =LOAD(errors, -2,  0),
					dW     =LOAD(errors, -1,  0);
#undef  LOAD
#define LOAD(X, Y) (uint32_t)(kx+(X))<(uint32_t)src->iw&&(uint32_t)(ky+(Y))<(uint32_t)src->ih?hireserror[src->iw*((ky+(Y))&1)+kx+(X)]:0
				int
					eNW=LOAD(-1, -1),
					eN =LOAD( 0, -1),
					eNE=LOAD( 1, -1),
					eW =LOAD(-1,  0);
#undef  LOAD
				int preds[G2_NPRED], pred;
				int prederrs[G2_NPRED];
				int paper_GAP, calic_GAP;
				int64_t num;
				int weights[G2_NPRED], den;
				int idx;

				for(k=0;k<G2_NPRED;++k)
				{
					/*eNW + eN + eNE*/
					prederrs[k]=
						((uint32_t)(kx-1)<(uint32_t)src->iw?perrors[src->iw*2*k+src->iw*((ky-1)&1)+kx-1]:0)+
						perrors[src->iw*2*k+src->iw*((ky-1)&1)+kx]+
						((uint32_t)(kx+1)<(uint32_t)src->iw?perrors[src->iw*2*k+src->iw*((ky-1)&1)+kx+1]:0);
				}
				{
					int dx=abs(W-WW)+abs(N-NW)+abs(NE-N);
					int dy=abs(W-NW)+abs(N-NN)+abs(NE-NNE);
					int d45=abs(W-NWW)+abs(NW-NNWW)+abs(N-NNW);
					int d135=abs(NE-NNEE)+abs(N-NNE)+abs(W-N);
					int diff=(dy-dx)>>8, diff2=(d45-d135)>>8, diff3=NE-NW;
					if(dy+dx>32)
						paper_GAP=(int)(((int64_t)dx*N+(int64_t)dy*W)/((int64_t)dy+dx));
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
						calic_GAP=(2*N+6*W+NE-NW)>>3;		/*c1	[1/4  3/4  1/8  -1/8].[N W NE NW]*/
					else if(diff>8)
						calic_GAP=(6*N+10*W+3*(NE-NW))>>4;	/*c2	[3/8  5/8  3/16  -3/16]*/
					else if(diff<-32)
						calic_GAP=(6*N+2*W+NE-NW)>>3;		/*c3	[3/4  1/4  1/8  -1/8]*/
					else if(diff<-8)
						calic_GAP=(10*N+6*W+3*(NE-NW))>>4;	/*c4	[5/8  3/8  3/16  -3/16]*/
					else
						calic_GAP=(((N+W)<<1)+NE-NW)>>2;	/*c5	[1/2  1/2  1/4  -1/4]*/
				}
				j=0;
				/*the 4 predictors from JPEG XL:*/
				preds[j++]=W+NE-N;
				preds[j++]=N-((eN+eW+eNE)*params[G2_NPRED+6]>>8);
				preds[j++]=W-((eN+eW+eNW)*params[G2_NPRED+5]>>8);
				preds[j++]=N-((eNW*params[G2_NPRED]+eN*params[G2_NPRED+1]+eNE*params[G2_NPRED+2]+(NN-N)*params[G2_NPRED+3]+(NW-W)*params[G2_NPRED+4])>>8);

				/*spam some more*/
				preds[j++]=W -(eW *params[G2_NPRED+ 7]>>8);
				preds[j++]=N -(eN *params[G2_NPRED+ 8]>>8);
				preds[j++]=NW-(eNW*params[G2_NPRED+ 9]>>8);
				preds[j++]=NE-(eNE*params[G2_NPRED+10]>>8);
				
				preds[j++]=calic_GAP;
				preds[j++]=(4*(N+W+NW+NE)-(NN+WW+NNWW+NNEE)+6)/12;
				preds[j++]=(3*W+NEEE)>>2;
				preds[j++]=paper_GAP;
				
				preds[j++]=(W+NEE)>>1;			/*E0.5		*/
				preds[j++]=NNNNNN;			/*N6		*/
				preds[j++]=(NEEEE+NEEEEEE)>>1;		/*NE5		*/
				preds[j++]=(WWWW+WWWWWW)>>1;		/*W5		*/
				preds[j++]=(N+W+NEEEEE+NEEEEEEE)>>2;	/*NE2.75	*/

				preds[j++]=3*(N-NN)+NNN;
				preds[j++]=(N+NNN)>>1;		/*this predicts NN*/
				preds[j++]=((N+W)*3-NW*2)>>2;

				preds[j++]=N+W-NW;
				
				num=0;
				den=0;
				for(k=0;k<G2_NPRED;++k)
				{
					weights[k]=(params[k]<<8)/(prederrs[k]+1);
					num+=(int64_t)preds[k]*weights[k];
					den+=weights[k];
				}
				pred=den?(int)(num/den):preds[0];

				idx=(src->iw*ky+kx)<<2|kc;
				{
					SseCell *pc[12];
					int hashval, val, curr, delta, dnb[12];

					/*Secondary Symbol Estimation (SSE)*/
#if 1
					j=0;
					dnb[j++]=dN   ;/*error neighbors		order matters slightly*/
					dnb[j++]=dW   ;
					dnb[j++]=dNW  ;
					dnb[j++]=dNE  ;
					dnb[j++]=dNN  ;
					dnb[j++]=dWW  ;
					dnb[j++]=dNNWW;
					dnb[j++]=dNNW ;
					dnb[j++]=dNNE ;
					dnb[j++]=dNNEE;
					dnb[j++]=dNWW ;
					dnb[j++]=dNEE ;
					for(k2=3;k2>=0;--k2)
					{
						hashval=(pred+128)>>(8+6)&3;
						for(k=0;k<3;++k)
						{
							hashval<<=3;
							hashval+=(int)(dnb[3*k2+k]*0xFC28>>(16+6)&3);
						}
						hashval&=0xFF;
						pc[k2]=sse+((size_t)k2<<8)+hashval;
						{
							int64_t sum=pc[k2][0]>>12;
							int count=pc[k2][0]&0xFFF;
							pred+=(int)((sum+(count>>1))/(count+(intptr_t)1));
						}
					}
					CLAMP3_32(pred, pred, N, W, NE);
					val=(pred+128)>>8;
					val^=-fwd;
					val+=fwd;
					val+=((int32_t*)HDR_DATA(src))[idx];
					if(enable_ma)
					{
						val<<=32-depth;
						val>>=32-depth;
					}
					b2[idx]=val;

					curr=pixels[idx]<<8;
					delta=curr-pred;
					for(k2=0;k2<4;++k2)
					{
						int64_t sum=pc[k2][0]>>12;
						int count=pc[k2][0]&0xFFF;
						if(count+1>640)
						{
							count>>=1;
							sum>>=1;
						}
						++count;
						sum+=delta;
						pc[k2][0]=sum<<12|count;
					}
#endif
					hireserror[src->iw*(ky&1)+kx]=curr-pred;
					for(k=0;k<G2_NPRED;++k)
					{
						int e=abs(curr-preds[k]);
						perrors[src->iw*2*k+src->iw*(ky&1)+kx]=e;
						if(kx<src->iw&&ky>0)
							perrors[src->iw*2*k+src->iw*((ky-1)&1)+kx+1]+=e;
					}
				}
			}
		}
	}
	memcpy(((int32_t*)HDR_DATA(src)), b2, res*sizeof(int[4]));
	free(perrors);
	free(sse);
	free(b2);
}




/*DWTs*/

ArrayHandle dwt2d_gensizes(int iw, int ih, int wstop, int hstop, int nstages_override)/*calculate dimensions of each DWT stage in descending order	FIXME FIXME FIXME FIXME remove this garbage*/
{
	int w2, h2, nstages;
	DWTSize *p;
	ArrayHandle sizes;

	int lw=31^LZCNT32(iw), lh=31^LZCNT32(ih), lmin=MINVAR(lw, lh);
	ARRAY_ALLOC(DWTSize, sizes, 0, 0, lmin, 0);
	if(wstop<3)
		wstop=3;
	if(hstop<3)
		hstop=3;
	nstages=0;
	p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);
	p->w=(uint16_t)iw;
	p->h=(uint16_t)ih;
	for(w2=iw, h2=ih;w2>=wstop&&h2>=hstop&&(!nstages_override||nstages<nstages_override);++nstages)
	{
		p=(DWTSize*)ARRAY_APPEND(sizes, 0, 1, 1, 0);

		w2>>=1;
		h2>>=1;

		p->w=(uint16_t)w2;
		p->h=(uint16_t)h2;
	}
	return sizes;
}

/*lazy DWT*/
static void dwt1d_lazy_fwd(int *buffer, int count, int stride, int *b2)
{
	int k, ks, nodd, extraeven, *odd, *even;

	nodd=count>>1;
	extraeven=count&1;
	odd=b2;
	even=b2+nodd;
	for(k=0, ks=0;k<nodd;++k, ks+=stride<<1)/*lazy wavelet: split into odd (low frequency) & even (high frequency)*/
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
static void dwt1d_lazy_inv(int *buffer, int count, int stride, int *b2)
{
	int k, ks, nodd, extraeven, *odd, *even;

	nodd=count>>1;
	extraeven=count&1;
	odd=b2;
	even=b2+nodd;
	for(k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];


	for(k=0, ks=0;k<nodd;++k, ks+=stride<<1)/*inv lazy wavelet: join even & odd*/
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_lazy_fwd(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw, rowlen, it, kx, ky;

	iw=sizes->w;
	rowlen=stride*iw;
	if(sizes_start>=sizes_end-1)
		return;
	for(it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;
		/*int w2=sizes[sizes_start].w, h2=sizes[sizes_start].h;*/

		for(ky=0;ky<h2;++ky)/*horizontal DWT*/
			dwt1d_lazy_fwd(buffer+rowlen*ky, w2, stride, temp);

		for(kx=0;kx<w2;++kx)/*vertical DWT*/
			dwt1d_lazy_fwd(buffer+stride*kx, h2, rowlen, temp);
	}
}
void dwt2d_lazy_inv(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp, int vmin, int vmax)
{
	int iw, rowlen, it, kx, ky;

	iw=sizes->w;
	rowlen=stride*iw;
	if(sizes_start>=sizes_end-1)
		return;
	for(it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;
		/*int w2=sizes[sizes_start].w, h2=sizes[sizes_start].h;*/

		for(kx=0;kx<w2;++kx)/*vertical IDWT*/
			dwt1d_lazy_inv(buffer+stride*kx, h2, rowlen, temp);

		for(ky=0;ky<h2;++ky)/*horizontal IDWT*/
			dwt1d_lazy_inv(buffer+rowlen*ky, w2, stride, temp);
	}

	(void)vmin;
	(void)vmax;
}

/*lifting-based 8bit Haar DWT*/
static void dwt1d_haar_fwd(int *buffer, int count, int stride, int *b2)
{
	int k, ks, nodd, extraeven, *odd, *even;

	nodd=count>>1;
	extraeven=count&1;
	odd=b2;
	even=b2+nodd;
	for(k=0, ks=0;k<nodd;++k, ks+=stride<<1)/*lazy wavelet: split into odd (low frequency) & even (high frequency)*/
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];


	for(k=0;k<nodd;++k)/*const predictor (difference)*/
		even[k]-=odd[k];
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(k=0;k<nodd;++k)/*update (O+(E-O)/2 = average)*/
		odd[k]+=even[k]>>1;


	for(k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
static void dwt1d_haar_inv(int *buffer, int count, int stride, int *b2)
{
	int k, ks, nodd, extraeven, *odd, *even;

	nodd=count>>1;
	extraeven=count&1;
	odd=b2;
	even=b2+nodd;
	for(k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	
	for(k=0;k<nodd;++k)/*unupdate*/
		odd[k]-=even[k]>>1;
	
	for(k=0;k<nodd;++k)/*unpredict*/
		even[k]+=odd[k];
	if(extraeven)
		even[nodd]+=odd[nodd-1];


	for(k=0, ks=0;k<nodd;++k, ks+=stride<<1)/*inv lazy wavelet: join even & odd*/
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_haar_fwd(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw, rowlen, it, kx, ky;

	iw=sizes->w;
	rowlen=stride*iw;
	if(sizes_start>=sizes_end-1)
		return;
	for(it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(ky=0;ky<h2;++ky)/*horizontal DWT*/
			dwt1d_haar_fwd(buffer+rowlen*ky, w2, stride, temp);

		for(kx=0;kx<w2;++kx)/*vertical DWT*/
			dwt1d_haar_fwd(buffer+stride*kx, h2, rowlen, temp);
	}
}
void dwt2d_haar_inv(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp, int vmin, int vmax)
{
	intptr_t k, size;
	int iw, rowlen, it, kx, ky;

	iw=sizes->w;
	rowlen=stride*iw;
	if(sizes_start>=sizes_end-1)
		return;
	for(it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(kx=0;kx<w2;++kx)/*vertical IDWT*/
			dwt1d_haar_inv(buffer+stride*kx, h2, rowlen, temp);

		for(ky=0;ky<h2;++ky)/*horizontal IDWT*/
			dwt1d_haar_inv(buffer+rowlen*ky, w2, stride, temp);
	}
	for(k=0, size=(intptr_t)stride*sizes->w*sizes->h;k<size;k+=stride)
	{
		CLAMP2(buffer[k], vmin, vmax);
	}
}

/*lifting-based 8bit squeeze DWT  from JPEG XL*/
static int smoothtendency(int B, int a, int n)
{
	int diff=0;
	if(B>=a&&a>=n)
	{
		diff=(4*B-3*n-a+6)/12;
		/*
		2C = a<<1 + diff - diff&1 <= 2B  so diff - diff&1 <= 2B - 2a
		2D = a<<1 - diff - diff&1 >= 2n  so diff + diff&1 <= 2a - 2n
		*/
		if(diff-(diff&1)>2*(B-a))
			diff=2*(B-a)+1;
		if(diff+(diff&1)>2*(a-n))
			diff=2*(a-n);
	}
	else if(B<=a&&a<=n)
	{
		diff=(4*B-3*n-a-6)/12;
		/*
		2C = a<<1 + diff + diff&1 >= 2B  so diff + diff&1 >= 2B - 2a
		2D = a<<1 - diff + diff&1 <= 2n  so diff - diff&1 >= 2a - 2n
		*/
		if(diff+(diff&1)<2*(B-a))
			diff=2*(B-a)-1;
		if(diff-(diff&1)<2*(a-n))
			diff=2*(a-n);
	}
	return diff;
}
static void dwt1d_squeeze_fwd(int *buffer, int count, int stride, int *b2)/*nOdd <= nEven			nOdd = nEven - (count&1)*/
{
	int k, ks, kd, nodd, extraeven, *odd, *even;

	nodd=count>>1;
	extraeven=count&1;
	odd=b2;
	even=b2+nodd;
	for(ks=0, kd=0;kd<nodd+extraeven;ks+=stride<<1, ++kd)
	{
		int
			o1=kd>0?buffer[ks-stride]:0,
			e =buffer[ks],
			o =(kd<<1)+1<count?buffer[ks+stride]:0,
			e2=(kd<<1)+2<count?buffer[ks+(stride<<1)]:0,/*n-1-(idx-(n-1)) = ((n-1)<<1)-idx*/
			o2=(kd<<1)+3<count?buffer[ks+ stride*3  ]:0;

		e-=o;		/*diff*/
		o+=e>>1;	/*av*/
		e2-=o2;
		o2+=e2>>1;
		e-=smoothtendency(o1, o, o2);
		
		if(kd<nodd)
			odd[kd]=o;
		even[kd]=e;
	}


	for(k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
static void dwt1d_squeeze_inv(int *buffer, int count, int stride, int *b2)
{
	int k, ks, kd, nodd, extraeven, *odd, *even;

	nodd=count>>1;
	extraeven=count&1;
	odd=b2;
	even=b2+nodd;
	for(k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];
	
	for(ks=0, kd=0;ks<nodd+extraeven;++ks, kd+=stride<<1)
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
	int iw, rowlen, it, kx, ky;

	iw=sizes->w;
	rowlen=stride*iw;
	if(sizes_start>=sizes_end-1)
		return;
	for(it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(ky=0;ky<h2;++ky)/*horizontal DWT*/
			dwt1d_squeeze_fwd(buffer+rowlen*ky, w2, stride, temp);

		for(kx=0;kx<w2;++kx)/*vertical DWT*/
			dwt1d_squeeze_fwd(buffer+stride*kx, h2, rowlen, temp);
	}
}
void dwt2d_squeeze_inv(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp, int vmin, int vmax)
{
	intptr_t k, size;
	int iw, rowlen, it, kx, ky;

	iw=sizes->w;
	rowlen=stride*iw;
	if(sizes_start>=sizes_end-1)
		return;
	for(it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(kx=0;kx<w2;++kx)/*vertical IDWT*/
			dwt1d_squeeze_inv(buffer+stride*kx, h2, rowlen, temp);

		for(ky=0;ky<h2;++ky)/*horizontal IDWT*/
			dwt1d_squeeze_inv(buffer+rowlen*ky, w2, stride, temp);
	}
	for(k=0, size=(intptr_t)stride*sizes->w*sizes->h;k<size;k+=stride)
	{
		CLAMP2(buffer[k], vmin, vmax);
	}
}

/*lifting-based 8bit LeGall 5/3 DWT*/
static void dwt1d_legall53_fwd(int *buffer, int count, int stride, int *b2)
{
	int k, ks, nodd, extraeven, *odd, *even;

	nodd=count>>1;
	extraeven=count&1;
	odd=b2;
	even=b2+nodd;
	for(k=0, ks=0;k<nodd;++k, ks+=stride<<1)/*lazy wavelet: split into odd (low frequency) & even (high frequency)*/
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

	even[0]-=odd[0];
	for(k=1;k<nodd;++k)/*predict linear (deviation from linearity)*/
		even[k]-=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]-=odd[nodd-1];
	
	for(k=0;k<nodd-!extraeven;++k)/*update (smoothing)*/
		odd[k]+=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]+=even[nodd-1]>>1;

	for(k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
static void dwt1d_legall53_inv(int *buffer, int count, int stride, int *b2)
{
	int k, ks, nodd, extraeven, *odd, *even;

	nodd=count>>1;
	extraeven=count&1;
	odd=b2;
	even=b2+nodd;
	for(k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];

	
	for(k=0;k<nodd-!extraeven;++k)/*un-update*/
		odd[k]-=(even[k]+even[k+1])>>2;
	if(!extraeven)
		odd[nodd-1]-=even[nodd-1]>>1;
	
	even[0]+=odd[0];
	for(k=1;k<nodd;++k)/*un-predict*/
		even[k]+=(odd[k-1]+odd[k])>>1;
	if(extraeven)
		even[nodd]+=odd[nodd-1];


	for(k=0, ks=0;k<nodd;++k, ks+=stride<<1)/*inv lazy wavelet: join even & odd*/
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_legall53_fwd(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw, rowlen, it, kx, ky;

	iw=sizes->w;
	rowlen=stride*iw;
	if(sizes_start>=sizes_end-1)
		return;
	for(it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(ky=0;ky<h2;++ky)/*horizontal DWT*/
			dwt1d_legall53_fwd(buffer+rowlen*ky, w2, stride, temp);

		for(kx=0;kx<w2;++kx)/*vertical DWT*/
			dwt1d_legall53_fwd(buffer+stride*kx, h2, rowlen, temp);
	}
}
void dwt2d_legall53_inv(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp, int vmin, int vmax)
{
	intptr_t k, size;
	int iw, rowlen, it, kx, ky;

	iw=sizes->w;
	rowlen=stride*iw;
	if(sizes_start>=sizes_end-1)
		return;
	for(it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(kx=0;kx<w2;++kx)/*vertical IDWT*/
			dwt1d_legall53_inv(buffer+stride*kx, h2, rowlen, temp);

		for(ky=0;ky<h2;++ky)/*horizontal IDWT*/
			dwt1d_legall53_inv(buffer+rowlen*ky, w2, stride, temp);
	}
	for(k=0, size=(intptr_t)stride*sizes->w*sizes->h;k<size;k+=stride)
	{
		CLAMP2(buffer[k], vmin, vmax);
	}
}

/*lifting-based CDF 9/7 DWT*/
static const int cdf97_coeff[]=
{
	/*
	'Factoring wavelet transforms into lifting steps' - page 19
	'A wavelet tour of signal processing - the sparse way' - page 376
	Lifting-based Discrete Wavelet Transform for Real-Time Signal Detection 2015-10
	*/
	 0x1960C,	/* 1.58613434342059f,	alpha		predict is negative	*/
	-0x00D90,	/*-0.0529801185729f,	beta		update is positive	*/
	-0x0E206,	/*-0.8829110755309f,	gamma		predict			*/
	 0x07189,	/* 0.4435068520439f,	delta		update			*/
	 0x1264C,	/* 1.1496043988602f,	zeta		output gain is 1.89	*/
};
static void dwt1d_cdf97_predict(int *odd, int *even, int nodd, int extraeven, int coeff)
{
	int k, p;

	p=odd[0]*coeff;
	even[0]-=p<0?-(-p>>15):p>>15;
	for(k=1;k<nodd;++k)/*predict*/
	{
		p=(odd[k-1]+odd[k])*coeff;
		even[k]-=p<0?-(-p>>16):p>>16;
	}
	if(extraeven)
	{
		p=odd[nodd-1]*coeff;
		even[nodd]-=p<0?-(-p>>15):p>>15;
	}
}
static void dwt1d_cdf97_update(int *odd, int *even, int nodd, int extraeven, int coeff)
{
	int k, p;

	for(k=0;k<nodd-!extraeven;++k)/*update*/
	{
		p=(even[k]+even[k+1])*coeff;
		odd[k]+=p<0?-(-p>>16):p>>16;
	}
	if(!extraeven)
	{
		p=even[nodd-1]*coeff;
		odd[nodd-1]+=p<0?-(-p>>15):p>>15;
	}
}
static void dwt1d_cdf97_fwd(int *buffer, int count, int stride, int *b2)
{
	int k, ks, nodd, extraeven, *odd, *even;

	nodd=count>>1;
	extraeven=count&1;
	odd=b2;
	even=b2+nodd;
	for(k=0, ks=0;k<nodd;++k, ks+=stride<<1)/*lazy wavelet: split into odd (low frequency) & even (high frequency)*/
	{
		even[k]=buffer[ks];
		odd[k]=buffer[ks+stride];
	}
	if(extraeven)
		even[nodd]=buffer[stride*(count-1)];

	dwt1d_cdf97_predict(odd, even, nodd, extraeven, cdf97_coeff[0]);
	dwt1d_cdf97_update (odd, even, nodd, extraeven, cdf97_coeff[1]);
	dwt1d_cdf97_predict(odd, even, nodd, extraeven, cdf97_coeff[2]);
	dwt1d_cdf97_update (odd, even, nodd, extraeven, cdf97_coeff[3]);

	for(k=0, ks=0;k<count;++k, ks+=stride)
		buffer[ks]=b2[k];
}
static void dwt1d_cdf97_inv(int *buffer, int count, int stride, int *b2)
{
	int k, ks, nodd, extraeven, *odd, *even;

	nodd=count>>1;
	extraeven=count&1;
	odd=b2;
	even=b2+nodd;
	for(k=0, ks=0;k<count;++k, ks+=stride)
		b2[k]=buffer[ks];

	dwt1d_cdf97_update (odd, even, nodd, extraeven, -cdf97_coeff[3]);
	dwt1d_cdf97_predict(odd, even, nodd, extraeven, -cdf97_coeff[2]);
	dwt1d_cdf97_update (odd, even, nodd, extraeven, -cdf97_coeff[1]);
	dwt1d_cdf97_predict(odd, even, nodd, extraeven, -cdf97_coeff[0]);

	for(k=0, ks=0;k<nodd;++k, ks+=stride<<1)/*inv lazy wavelet: join even & odd*/
	{
		buffer[ks]=even[k];
		buffer[ks+stride]=odd[k];
	}
	if(extraeven)
		buffer[stride*(count-1)]=even[nodd];
}
void dwt2d_cdf97_fwd(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp)
{
	int iw, rowlen, it, kx, ky;

	iw=sizes->w;
	rowlen=stride*iw;
	if(sizes_start>=sizes_end-1)
		return;
	for(it=sizes_start;it<sizes_end-1;++it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(ky=0;ky<h2;++ky)/*horizontal DWT*/
			dwt1d_cdf97_fwd(buffer+rowlen*ky, w2, stride, temp);

		for(kx=0;kx<w2;++kx)/*vertical DWT*/
			dwt1d_cdf97_fwd(buffer+stride*kx, h2, rowlen, temp);
	}
	if(g_dist>1)
	{
		int invdist, *ptr, w0, h0;

		invdist=(0x10000+g_dist-1)/g_dist;
		ptr=buffer;
		w0=sizes[sizes_end-1].w;
		h0=sizes[sizes_end-1].h;
		for(ky=0;ky<sizes->h;++ky)
		{
			for(kx=0;kx<sizes->w;++kx, ptr+=stride)
			{
				if(ky>=h0&&kx>=w0)
				{
					int val=*ptr;
					val=val*invdist>>16;
					val-=val>>31;
					*ptr=val;
				}
			}
		}
	}
}
void dwt2d_cdf97_inv(int *buffer, DWTSize *sizes, int sizes_start, int sizes_end, int stride, int *temp, int vmin, int vmax)
{
	intptr_t k, size;
	int iw, rowlen, it, kx, ky;

	iw=sizes->w;
	rowlen=stride*iw;
	if(sizes_start>=sizes_end-1)
		return;
	if(g_dist>1)
	{
		int *ptr=buffer;
		int w0=sizes[sizes_end-1].w;
		int h0=sizes[sizes_end-1].h;
		for(ky=0;ky<sizes->h;++ky)
		{
			for(kx=0;kx<sizes->w;++kx, ptr+=stride)
			{
				if(ky>=h0&&kx>=w0)
				{
					int val=*ptr;
					val=val*g_dist;
					val-=val>>31;
					*ptr=val;
				}
			}
		}
	}
	for(it=sizes_end-2;it>=sizes_start;--it)
	{
		int w2=sizes[it].w, h2=sizes[it].h;

		for(kx=0;kx<w2;++kx)/*vertical IDWT*/
			dwt1d_cdf97_inv(buffer+stride*kx, h2, rowlen, temp);

		for(ky=0;ky<h2;++ky)/*horizontal IDWT*/
			dwt1d_cdf97_inv(buffer+rowlen*ky, w2, stride, temp);
	}
	for(k=0, size=(intptr_t)stride*sizes->w*sizes->h;k<size;k+=stride)
	{
		CLAMP2(buffer[k], vmin, vmax);
	}
}


/*DCTs*/
static void dct4_fwd_i8(int *x)
{
#if 0
	int a[4], b[4];

	a[0]=x[0]+x[3];
	a[1]=x[1]+x[2];
	a[2]=x[1]-x[2];
	a[3]=x[0]-x[3];

	b[0]=a[0]+a[1];
	b[1]=(a[0]-a[1])*46341>>16;
	b[2]=(a[3]*60547+a[2]*25080)>>16;
	b[3]=(a[3]*25080-a[2]*60547)>>16;

	x[0]=b[0]>>2;
	x[1]=b[2];
	x[2]=b[1];
	x[3]=b[3];
#endif
#if 1
	x[3]=x[0]-x[3];
	x[0]-=x[3]>>1;
	x[2]=x[1]-x[2];
	x[1]-=x[2]>>1;
	x[1]=x[0]-x[1];
	x[0]-=x[1]>>1;
	x[2]=(x[3]*13>>5)-x[2];
	x[3]-=x[2]*11>>5;
#endif
}
static void dct4_inv_i8(int *x)
{
#if 0
	int a[4], b[4];

	a[0]=x[0]<<1;/*half*/
	a[2]=x[1];
	a[1]=x[2];
	a[3]=x[3];

	a[1]=a[1]*46341>>16;

	b[0]=a[0]+a[1];
	b[1]=a[0]-a[1];
	b[2]=(a[2]*25080-a[3]*60547)>>16;
	b[3]=(a[2]*60547+a[3]*25080)>>16;

	x[0]=(b[0]+b[3])>>1;
	x[1]=(b[1]+b[2])>>1;
	x[2]=(b[1]-b[2])>>1;
	x[3]=(b[0]-b[3])>>1;
#endif
#if 1
	x[3]+=x[2]*11>>5;
	x[2]=(x[3]*13>>5)-x[2];
	x[0]+=x[1]>>1;
	x[1]=x[0]-x[1];
	x[1]+=x[2]>>1;
	x[2]=x[1]-x[2];
	x[0]+=x[3]>>1;
	x[3]=x[0]-x[3];
#endif
}
void image_dct4_fwd(Image *image)
{
	int tsize, *temp, j, kc, ky, kx;

	tsize=MAXVAR(image->iw, image->ih)*(int)sizeof(int);
	temp=(int*)malloc(tsize);
	if(!temp)
	{
		CRASH("Alloc error");
		return;
	}
	memset(temp, 0, tsize);
	for(kc=0;kc<4;++kc)
	{
		if(!image->depth[kc])
			continue;
#if 1
		for(ky=0;ky<image->ih;++ky)
		{
			for(kx=0;kx<image->iw-3;kx+=4)
			{
				int idx, x[4];

				idx=image->iw*ky+kx;
				j=0;
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+0)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+1)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+2)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+3)<<2|kc];

				dct4_fwd_i8(x);

				temp[(kx>>2)+(image->iw>>2)*0]=x[0];
				temp[(kx>>2)+(image->iw>>2)*1]=x[1];
				temp[(kx>>2)+(image->iw>>2)*2]=x[2];
				temp[(kx>>2)+(image->iw>>2)*3]=x[3];
			}
			for(kx=0;kx<(image->iw&~3);++kx)
				((int32_t*)HDR_DATA(image))[(image->iw*ky+kx)<<2|kc]=temp[kx];
		}
#endif
#if 1
		for(kx=0;kx<image->iw;++kx)
		{
			for(ky=0;ky<image->ih-3;ky+=4)
			{
				int idx, x[4];

				idx=image->iw*ky+kx;
				j=0;
				x[j++]=((int32_t*)HDR_DATA(image))[ idx             <<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+image->iw  )<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+image->iw*2)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+image->iw*3)<<2|kc];

				dct4_fwd_i8(x);

				temp[(ky>>2)                 ]=x[0];
				temp[(ky>>2)+(image->ih>>2)  ]=x[1];
				temp[(ky>>2)+(image->ih>>2)*2]=x[2];
				temp[(ky>>2)+(image->ih>>2)*3]=x[3];
			}
			for(ky=0;ky<(image->ih&~3);++ky)
				((int32_t*)HDR_DATA(image))[(image->iw*ky+kx)<<2|kc]=temp[ky];
		}
#endif
	}
	free(temp);
	for(kc=0;kc<3;++kc)
		image->depth[kc]+=(image->depth[kc]+2<=24)<<1;
}
void image_dct4_inv(Image *image)
{
	int amin[4], amax[4], tsize, *temp, j, kc, ky, kx;

	tsize=MAXVAR(image->iw, image->ih)*(int)sizeof(int);
	temp=(int*)malloc(tsize);
	if(!temp)
	{
		CRASH("Alloc error");
		return;
	}
	memset(temp, 0, tsize);
	for(kc=0;kc<3;++kc)
		image->depth[kc]-=(image->depth[kc]-2>=image->src_depth[kc])<<1;
	amin[0]=-(1<<image->depth[0]>>1);
	amin[1]=-(1<<image->depth[1]>>1);
	amin[2]=-(1<<image->depth[2]>>1);
	amin[3]=-(1<<image->depth[3]>>1);
	amax[0]=(1<<image->depth[0]>>1)-1;
	amax[1]=(1<<image->depth[1]>>1)-1;
	amax[2]=(1<<image->depth[2]>>1)-1;
	amax[3]=(1<<image->depth[3]>>1)-1;
	for(kc=0;kc<4;++kc)
	{
		if(!image->depth[kc])
			continue;
#if 1
		for(kx=0;kx<image->iw;++kx)
		{
			for(ky=0;ky<(image->ih&~3);++ky)
				temp[ky]=((int32_t*)HDR_DATA(image))[(image->iw*ky+kx)<<2|kc];
			for(ky=0;ky<image->ih-3;ky+=4)
			{
				int idx, x[4];

				idx=image->iw*ky+kx;
				j=0;
				x[j++]=temp[(ky>>2)                 ];
				x[j++]=temp[(ky>>2)+(image->ih>>2)  ];
				x[j++]=temp[(ky>>2)+(image->ih>>2)*2];
				x[j++]=temp[(ky>>2)+(image->ih>>2)*3];

				dct4_inv_i8(x);
				
				((int32_t*)HDR_DATA(image))[ idx             <<2|kc]=x[0];
				((int32_t*)HDR_DATA(image))[(idx+image->iw  )<<2|kc]=x[1];
				((int32_t*)HDR_DATA(image))[(idx+image->iw*2)<<2|kc]=x[2];
				((int32_t*)HDR_DATA(image))[(idx+image->iw*3)<<2|kc]=x[3];
			}
		}
#endif
#if 1
		for(ky=0;ky<image->ih;++ky)
		{
			for(kx=0;kx<(image->iw&~3);++kx)
				temp[kx]=((int32_t*)HDR_DATA(image))[(image->iw*ky+kx)<<2|kc];
			for(kx=0;kx<image->iw-3;kx+=4)
			{
				int idx, x[4];

				idx=image->iw*ky+kx;
				j=0;
				x[j++]=temp[(kx>>2)                 ];
				x[j++]=temp[(kx>>2)+(image->iw>>2)  ];
				x[j++]=temp[(kx>>2)+(image->iw>>2)*2];
				x[j++]=temp[(kx>>2)+(image->iw>>2)*3];

				dct4_inv_i8(x);
				CLAMP2(x[0], amin[kc], amax[kc]);
				CLAMP2(x[1], amin[kc], amax[kc]);
				CLAMP2(x[2], amin[kc], amax[kc]);
				CLAMP2(x[3], amin[kc], amax[kc]);
				
				((int32_t*)HDR_DATA(image))[ idx   <<2|kc]=x[0];
				((int32_t*)HDR_DATA(image))[(idx+1)<<2|kc]=x[1];
				((int32_t*)HDR_DATA(image))[(idx+2)<<2|kc]=x[2];
				((int32_t*)HDR_DATA(image))[(idx+3)<<2|kc]=x[3];
			}
		}
#endif
	}
	free(temp);
}

static void dct8_fwd_i8(int *x)
{
	/*binDCT-C7*/
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
	{
		int a[8];

		a[0]=x[0];
		a[4]=x[1];
		a[6]=x[2];
		a[2]=x[3];
		a[7]=x[4];
		a[5]=x[5];
		a[3]=x[6];
		a[1]=x[7];

		x[0]=a[0];
		x[1]=a[1];
		x[2]=a[2];
		x[3]=a[3];
		x[4]=a[4];
		x[5]=a[5];
		x[6]=a[6];
		x[7]=a[7];
	}
}
static void dct8_inv_i8(int *x)
{
	{
		int a[8];

		a[0]=x[0];
		a[1]=x[1];
		a[2]=x[2];
		a[3]=x[3];
		a[4]=x[4];
		a[5]=x[5];
		a[6]=x[6];
		a[7]=x[7];

		x[0]=a[0];
		x[1]=a[4];
		x[2]=a[6];
		x[3]=a[2];
		x[4]=a[7];
		x[5]=a[5];
		x[6]=a[3];
		x[7]=a[1];
	}

	/*invBinDCT-C7*/
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
	int tsize, *temp, j, kc, kx, ky;

	tsize=MAXVAR(image->iw, image->ih)*(int)sizeof(int);
	temp=(int*)malloc(tsize);
	if(!temp)
	{
		CRASH("Alloc error");
		return;
	}
	memset(temp, 0, tsize);
	for(kc=0;kc<4;++kc)
	{
		if(!image->depth[kc])
			continue;
#if 1
		for(ky=0;ky<image->ih;++ky)
		{
			for(kx=0;kx<image->iw-7;kx+=8)
			{
				int idx, x[8];

				idx=image->iw*ky+kx;
				j=0;
				x[j++]=((int32_t*)HDR_DATA(image))[ idx   <<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+1)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+2)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+3)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+4)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+5)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+6)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+7)<<2|kc];


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
			for(kx=0;kx<(image->iw&~7);++kx)
				((int32_t*)HDR_DATA(image))[(image->iw*ky+kx)<<2|kc]=temp[kx];
		}
#endif
#if 1
		for(kx=0;kx<image->iw;++kx)
		{
			for(ky=0;ky<image->ih-7;ky+=8)
			{
				int idx, x[8];

				idx=image->iw*ky+kx;
				j=0;
				x[j++]=((int32_t*)HDR_DATA(image))[ idx             <<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+image->iw  )<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+image->iw*2)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+image->iw*3)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+image->iw*4)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+image->iw*5)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+image->iw*6)<<2|kc];
				x[j++]=((int32_t*)HDR_DATA(image))[(idx+image->iw*7)<<2|kc];

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
			for(ky=0;ky<(image->ih&~7);++ky)
				((int32_t*)HDR_DATA(image))[(image->iw*ky+kx)<<2|kc]=temp[ky];
		}
#endif
	}
	free(temp);
	for(kc=0;kc<3;++kc)
		image->depth[kc]+=(image->depth[kc]+3<=24)*3;
}
void image_dct8_inv(Image *image)
{
	int tsize, *temp, j, kc, kx, ky;

	tsize=MAXVAR(image->iw, image->ih)*sizeof(int);
	temp=(int*)malloc(tsize);
	if(!temp)
	{
		CRASH("Alloc error");
		return;
	}
	memset(temp, 0, tsize);
	for(kc=0;kc<3;++kc)
		image->depth[kc]-=(image->depth[kc]-3>=image->src_depth[kc])*3;
	for(kc=0;kc<4;++kc)
	{
		int vmin=-(1<<image->depth[kc]>>1), vmax=(1<<image->depth[kc]>>1)-1;
		if(!image->depth[kc])
			continue;
#if 1
		for(kx=0;kx<image->iw;++kx)
		{
			for(ky=0;ky<(image->ih&~7);++ky)
				temp[ky]=((int32_t*)HDR_DATA(image))[(image->iw*ky+kx)<<2|kc];
			for(ky=0;ky<image->ih-7;ky+=8)
			{
				int idx, x[8];

				idx=image->iw*ky+kx;
				j=0;
				x[j++]=temp[(ky>>3)                 ];
				x[j++]=temp[(ky>>3)+(image->ih>>3)  ];
				x[j++]=temp[(ky>>3)+(image->ih>>3)*2];
				x[j++]=temp[(ky>>3)+(image->ih>>3)*3];
				x[j++]=temp[(ky>>3)+(image->ih>>3)*4];
				x[j++]=temp[(ky>>3)+(image->ih>>3)*5];
				x[j++]=temp[(ky>>3)+(image->ih>>3)*6];
				x[j++]=temp[(ky>>3)+(image->ih>>3)*7];

				dct8_inv_i8(x);
				
				((int32_t*)HDR_DATA(image))[ idx             <<2|kc]=x[0];
				((int32_t*)HDR_DATA(image))[(idx+image->iw  )<<2|kc]=x[1];
				((int32_t*)HDR_DATA(image))[(idx+image->iw*2)<<2|kc]=x[2];
				((int32_t*)HDR_DATA(image))[(idx+image->iw*3)<<2|kc]=x[3];
				((int32_t*)HDR_DATA(image))[(idx+image->iw*4)<<2|kc]=x[4];
				((int32_t*)HDR_DATA(image))[(idx+image->iw*5)<<2|kc]=x[5];
				((int32_t*)HDR_DATA(image))[(idx+image->iw*6)<<2|kc]=x[6];
				((int32_t*)HDR_DATA(image))[(idx+image->iw*7)<<2|kc]=x[7];
			}
		}
#endif
#if 1
		for(ky=0;ky<image->ih;++ky)
		{
			for(kx=0;kx<(image->iw&~7);++kx)
				temp[kx]=((int32_t*)HDR_DATA(image))[(image->iw*ky+kx)<<2|kc];
			for(kx=0;kx<image->iw-7;kx+=8)
			{
				int idx, x[8];

				idx=image->iw*ky+kx;
				j=0;
				x[j++]=temp[(kx>>3)                 ];
				x[j++]=temp[(kx>>3)+(image->iw>>3)  ];
				x[j++]=temp[(kx>>3)+(image->iw>>3)*2];
				x[j++]=temp[(kx>>3)+(image->iw>>3)*3];
				x[j++]=temp[(kx>>3)+(image->iw>>3)*4];
				x[j++]=temp[(kx>>3)+(image->iw>>3)*5];
				x[j++]=temp[(kx>>3)+(image->iw>>3)*6];
				x[j++]=temp[(kx>>3)+(image->iw>>3)*7];

				dct8_inv_i8(x);
				CLAMP2(x[0], vmin, vmax);
				CLAMP2(x[1], vmin, vmax);
				CLAMP2(x[2], vmin, vmax);
				CLAMP2(x[3], vmin, vmax);
				CLAMP2(x[4], vmin, vmax);
				CLAMP2(x[5], vmin, vmax);
				CLAMP2(x[6], vmin, vmax);
				CLAMP2(x[7], vmin, vmax);
				
				((int32_t*)HDR_DATA(image))[ idx   <<2|kc]=x[0];
				((int32_t*)HDR_DATA(image))[(idx+1)<<2|kc]=x[1];
				((int32_t*)HDR_DATA(image))[(idx+2)<<2|kc]=x[2];
				((int32_t*)HDR_DATA(image))[(idx+3)<<2|kc]=x[3];
				((int32_t*)HDR_DATA(image))[(idx+4)<<2|kc]=x[4];
				((int32_t*)HDR_DATA(image))[(idx+5)<<2|kc]=x[5];
				((int32_t*)HDR_DATA(image))[(idx+6)<<2|kc]=x[6];
				((int32_t*)HDR_DATA(image))[(idx+7)<<2|kc]=x[7];
			}
		}
#endif
	}
	free(temp);
}


typedef struct HybridUintStruct
{
	uint16_t token, nbits;
	uint32_t bypass;
} HybridUint;
static void quantize_signed(int val, int expbits, int msb, int lsb, HybridUint *hu)
{
	int token, bypass, nbits;

	val=val<<1^-(val<0);
	if(val<(1<<expbits))
	{
		token=val;/*token*/
		nbits=0;
		bypass=0;
	}
	else
	{
		int lgv=31^LZCNT32((uint32_t)val);
		int mantissa=val-(1<<lgv);
		token = (1<<expbits) + (
				(lgv-expbits)<<(msb+lsb)|
				(mantissa>>(lgv-msb))<<lsb|
				(mantissa&((1<<lsb)-1))
			);
		nbits=lgv-(msb+lsb);
		bypass=val>>lsb&(((intptr_t)1<<nbits)-1);
	}
	hu->token=(uint16_t)token;
	hu->nbits=(uint16_t)nbits;
	hu->bypass=bypass;
}
const char* ec_method_label(EContext ec_method)
{
	const char *label="ERROR";

	switch(ec_method)
	{
#define CASE(L) case L:label=#L;break;
	CASE(ECTX_QNW)
	CASE(ECTX_MIN_QN_QW)
	CASE(ECTX_MAX_QN_QW)
	CASE(ECTX_MIN_N_W_NW_NE)
	CASE(ECTX_ARGMIN_N_W_NW_NE)
	CASE(ECTX_ABAC0)
	CASE(ECTX_ABAC1)
	CASE(ECTX_C2)
	CASE(ECTX_GRCTX)
	CASE(ECTX_GR)
	CASE(ECTX_STATIC_O0)
	CASE(ECTX_STATIC_O1)
	CASE(ECTX_INTERLEAVED)
	CASE(ECTX_DWT)
	CASE(ECTX_YUV422)
	CASE(ECTX_YUV420)
	CASE(ECTX_ZERO)
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
	NB_WW,   NB_W
} NBIndex;
static void getnb(Image const *src, int kc, int kx, int ky, int *nb)
{
	int ky2, kx2, idx2;

	for(ky2=-2, idx2=0;ky2<=0;++ky2)
	{
		for(kx2=-2;kx2<=2;++kx2, ++idx2)
		{
			if(!ky2&&!kx2)
				break;
			if((uint32_t)(ky+ky2)<(uint32_t)src->ih&&(uint32_t)(kx+kx2)<(uint32_t)src->iw)
				nb[idx2]=((const int32_t*)HDR_DATC(src))[(src->iw*(ky+ky2)+kx+kx2)<<2|kc];
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
		HybridUint hu;
		int nb2[4], val;

		nb2[0]=nb[NB_W];
		nb2[1]=nb[NB_NW];
		nb2[2]=nb[NB_N];
		nb2[3]=nb[NB_NE];
		val=nb2[idx];
		quantize_signed(val, expbits, msb, lsb, &hu);
		return hu.token;
	}
}
void calc_csize_ec(Image const *src, EContext method, int adaptive, int expbits, int msb, int lsb, double *entropy)
{
	int (*const getctx[])(const int *nb, int expbits, int msb, int lsb)=
	{
		getctx_QNW,
		getctx_min_QN_QW,
		getctx_max_QN_QW,
		getctx_min_N_W_NW_NE,
		getctx_argmin_N_W_NW_NE,
		getctx_zero,/*unused*/
		getctx_zero,/*unused*/
		getctx_zero,/*unused*/
		getctx_zero,/*unused*/
		getctx_zero,/*unused*/
		getctx_zero,/*unused*/
		getctx_zero,/*unused*/
		getctx_zero,/*unused*/
		getctx_zero,/*unused*/
		getctx_zero,/*unused*/
		getctx_zero,/*unused*/
		getctx_zero,
	};
	int k, ks, ky, kx, kc, idx, maxdepth;
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
		CRASH("Alloc error");
		return;
	}
	for(kc=0;kc<4;++kc)
	{
		int depth=src->depth[kc];
		int constant=1;

		if(!depth)
			continue;
		for(k=1;k<res;++k)/*check for constant value*/
		{
			if(((const int32_t*)HDR_DATC(src))[k<<2|kc]!=((const int32_t*)HDR_DATC(src))[kc])
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
			for(ky=0, idx=0;ky<src->ih;++ky)
			{
				for(kx=0;kx<src->iw;++kx, ++idx)
				{
					int nb[12], ctx, val, *curr_hist, den;

					getnb(src, kc, kx, ky, nb);
					ctx=getctx[method](nb, expbits, msb, lsb);
					val=((const int32_t*)HDR_DATC(src))[idx<<2|kc];
					quantize_signed(val, expbits, msb, lsb, &hu);

					curr_hist=hist+cdfsize*ctx;
					den=0;
					for(k=0;k<cdfsize;++k)
						den+=curr_hist[k];
					{
						double p=(double)curr_hist[hu.token]/den;
						double bitsize=-log(p)*(1./M_LN2);
						bitsizes[kc]+=bitsize;
						bypasssizes[kc]+=hu.nbits;
					}

					++curr_hist[hu.token];
					if(curr_hist[hu.token]>=adaptive)
					{
						for(k=0;k<cdfsize;++k)
							curr_hist[k]=(curr_hist[k]+1)>>1;
					}
				}
			}
			continue;
		}
		memset(hist, 0, sizeof(int)*cdfsize*cdfsize);
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx, ++idx)
			{
				int nb[12], ctx, val;

				getnb(src, kc, kx, ky, nb);
				ctx=getctx[method](nb, expbits, msb, lsb);

				val=((const int32_t*)HDR_DATC(src))[idx<<2|kc];
				quantize_signed(val, expbits, msb, lsb, &hu);
				val=hu.token;
				++hist[cdfsize*ctx+val];
			}
		}
		for(k=0;k<cdfsize;++k)
		{
			int sum=0;
			for(ks=0;ks<cdfsize;++ks)
				sum+=hist[cdfsize*k+ks];
			hsum[k]=sum;
		}
		for(ky=0, idx=0;ky<src->ih;++ky)
		{
			for(kx=0;kx<src->iw;++kx, ++idx)
			{
				int nb[12], ctx, val;

				getnb(src, kc, kx, ky, nb);
				ctx=getctx[method](nb, expbits, msb, lsb);

				val=((const int32_t*)HDR_DATC(src))[idx<<2|kc];
				quantize_signed(val, expbits, msb, lsb, &hu);
				val=hu.token;
				{
					int freq=hist[cdfsize*ctx+val], sum=hsum[ctx];
					if(!freq||!sum)
						CRASH("ZPS");
					{
						double p=(double)freq/sum;
						double bitsize=-log(p)*(1./M_LN2);
						bitsizes[kc]+=bitsize;
						bypasssizes[kc]+=hu.nbits;
					}
				}
			}
		}
	}
	{
		int nch=(src->src_depth[0]!=0)+(src->src_depth[1]!=0)+(src->src_depth[2]!=0)+(src->src_depth[3]!=0);
		double chubitsize=image_getBMPsize(src)*8/nch;
		for(kc=0;kc<4;++kc)
		{
			entropy[kc]=(bitsizes[kc]+bypasssizes[kc])/chubitsize;
			bitsizes[kc]/=8;
			bypasssizes[kc]/=8;
		}
	}
	free(hist);
	free(hsum);
}


void calc_csize_abac(Image const *src, int order, double *entropy)
{
	int64_t cbitsizes[4]={0};
	int ky, kx, kc, kb, idx, tidx, treesizes[4], maxtreesize, statssize;
	int32_t *stats;

	treesizes[0]=1<<src->depth[0];
	treesizes[1]=1<<src->depth[1];
	treesizes[2]=1<<src->depth[2];
	treesizes[3]=1<<src->depth[3];
	maxtreesize=treesizes[0];
	if(maxtreesize<treesizes[1])maxtreesize=treesizes[1];
	if(maxtreesize<treesizes[2])maxtreesize=treesizes[2];
	if(maxtreesize<treesizes[3])maxtreesize=treesizes[3];
	statssize=(int)sizeof(uint32_t)*maxtreesize;
	if(order)
		statssize*=maxtreesize;
	stats=(int32_t*)malloc(statssize);
	if(!stats)
	{
		CRASH("Alloc error");
		return;
	}
	memset(stats, 0, statssize);
	for(kc=0;kc<4;++kc)
	{
		uint64_t lo=0, hi=0xFFFFFFFFFFFF;
		int depth=src->depth[kc];
		if(!depth)
			continue;
	/*	FILLMEM(stats, 0x8000, statssize, sizeof(int16_t));*/
	/*	memset(stats, 0, statssize);*/
		for(ky=0, idx=kc;ky<src->ih;++ky)
		{
			int NW=0, N=0, WWW=0, WW=0, W=0;
			for(kx=0;kx<src->iw;++kx, idx+=4)
			{
				int32_t *curr_stats;
				int NE, val;

				curr_stats=stats;
				NE=0;
				if(order)
				{
					int ctx=W;

					NE=ky&&kx<src->iw-1?((const int32_t*)HDR_DATC(src))[idx-(src->iw<<2)+4]:0;
					if(ky)
						ctx+=N-NW;
					ctx=(ctx&((1<<depth)-1))<<depth;
					if(ctx>=(int)(statssize/sizeof(uint32_t)))
						CRASH("ABAC says hi");
					curr_stats+=ctx;
				}
				val=((const int32_t*)HDR_DATC(src))[idx];
				for(kb=depth-1, tidx=1;kb>=0;--kb)
				{
					enum
					{
						USEB=12,
						STOREB=20,
						AGEBITS=9,
						AGEMASK=(1<<AGEBITS)-1
					};
					uint64_t x;
					int32_t *pp, cell, prob, ctr, sh, p0, bit;

					pp=curr_stats+tidx;
					cell=*pp;
					p0=cell>>(AGEBITS+STOREB-USEB);
					prob=cell>>AGEBITS;
					ctr=cell&AGEMASK;
					sh=31^LZCNT32(ctr+3);
					ctr+=ctr<AGEMASK;
					p0+=p0<0;
					p0+=1<<USEB>>1;
					bit=val>>kb&1;
					x=hi-lo;
					if(x<=0xFFFF)
					{
						cbitsizes[kc]+=32;

						lo<<=32;
						hi=hi<<32|(uint32_t)~0;
						if(hi<lo)hi=~(uint64_t)0;
						x=hi-lo;
					}
					x=lo+(x*p0>>USEB);
					*(bit?&lo:&hi)=x-(bit^1);
					prob+=(((bit^1)<<STOREB)-(1<<STOREB>>1)-prob+(1<<sh>>1))>>sh;
				/*	p0+=((!bit<<16)-p0+(1<<7>>1))>>7;*/
					*pp=prob<<AGEBITS|ctr;
					tidx=tidx*2+bit;
				}
				NW=N;
				N=NE;
				WWW=WW;
				WW=W;
				W=val;
				(void)WWW;
			}
		}
		cbitsizes[kc]+=32;/*flush*/
	}
	free(stats);
	for(kc=0;kc<4;++kc)
	{
		double usize=(double)src->iw*src->ih*src->src_depth[kc];
		entropy[kc]=src->src_depth[kc]?cbitsizes[kc]/usize:0;/*entropy = cbitsize/ubitsize*/
	}
}